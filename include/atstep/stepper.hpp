// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/bdf.hpp"
#include "detail/controller.hpp"
#include "detail/output.hpp"
#include "detail/state.hpp"
#include "settings.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace atstep {

  /// What a step's solve is told about the step.
  ///
  /// A first-order equation uses the affine derivative `cj*y + known`, with
  /// `cj = derivative_weight()` and `known = derivative_history()`. Thus
  /// `M y' + K y = f` at `request.time` is
  ///
  ///     (cj M + K) y_{n+1} = f(t_{n+1}) - M known.
  ///
  /// `derivative(y, out)` evaluates it in correction form. Direct
  /// higher-derivative equations retain the classical level-based formula and
  /// its weights.
  ///
  /// A request borrows the stepper's states. It is valid for the solve
  /// callback's own call and never after: a solve that stores a reference to
  /// it, or to a state it returns, reads the next step's values.
  template <class State>
  class StepRequest {
  public:
    struct ClassicalDerivative {
      const detail::StepFormula& formula;
    };

    struct BdfDerivative {
      const State& predictor;
      const State& derivative_predictor;
      double leading;
    };

    using Derivative = std::variant<ClassicalDerivative, BdfDerivative>;

    StepRequest(double new_time,
      double step_size,
      int step_order,
      Derivative derivative_formula,
      const detail::Levels<State>& levels,
      const State& derivative_history,
      State& derivative_scratch,
      const State& absolute_scale,
      std::size_t accepted_steps,
      bool initialization = false)
      : time(new_time)
      , step(step_size)
      , order(step_order)
      , step_index(accepted_steps)
      , initializing(initialization)
      , derivative_(std::move(derivative_formula))
      , levels_(levels)
      , derivative_history_(derivative_history)
      , derivative_scratch_(derivative_scratch)
      , absolute_scale_(absolute_scale)
    {
    }

    const double time;
    const double step;
    const int order;
    const std::size_t step_index;
    const bool initializing;

    const State& absolute_tolerance_reference() const noexcept
    {
      return absolute_scale_;
    }

    /// The weight `y^(derivative)` carries on the state being solved for.
    double derivative_weight(int derivative = 1) const
    {
      if (const auto* bdf = std::get_if<BdfDerivative>(&derivative_)) {
        require_first_derivative(derivative);
        return bdf->leading;
      }
      return std::get<ClassicalDerivative>(derivative_)
        .formula.weights(derivative)[0];
    }

    /// The known term of the first derivative: `y' = cj*y + known`.
    const State& derivative_history() const noexcept
    {
      return derivative_history_;
    }

    /// The known term of the derivative level `derivative`, written to `out`.
    /// A level above the first is the classical stencil's, which a first-order
    /// run does not have.
    void derivative_history(int derivative, State& out) const
    {
      if (std::holds_alternative<BdfDerivative>(derivative_)) {
        require_first_derivative(derivative);
        out = derivative_history_;
        return;
      }
      if (derivative == 1) {
        out = derivative_history_;
        return;
      }
      const auto weights = std::get<ClassicalDerivative>(derivative_)
                             .formula.weights(derivative);
      out = previous();
      detail::axpy(out, weights[1] - 1.0, previous());
      for (std::size_t j = 2; j < weights.size(); ++j)
        detail::axpy(out, weights[j], previous(j - 1));
    }

    /// `out = y^(derivative)` of the state `y`. `out` is a state of the
    /// caller's and may not be one of this request's.
    ///
    /// A first-order equation evaluates its derivative in correction form,
    /// `predictor' + cj (y - predictor)`, rather than as `cj y + known`: the
    /// correction is small where the predictor is good, which is where the
    /// cancellation of the direct form would bite.
    void derivative(const State& y, State& out, int derivative = 1) const
    {
      if (const auto* bdf = std::get_if<BdfDerivative>(&derivative_)) {
        require_first_derivative(derivative);
        derivative_scratch_ = y;
        detail::axpy(derivative_scratch_, -1.0, bdf->predictor);
        out = bdf->derivative_predictor;
        detail::axpy(out, bdf->leading, derivative_scratch_);
        return;
      }
      derivative_history(derivative, out);
      detail::axpy(out, derivative_weight(derivative), y);
    }

    /// Views are valid for the solve callback; the stepper owns their
    /// states.
    ///
    /// Level 0 is the state the run has accepted. A first-order equation
    /// keeps no level behind it — the fixed-leading table holds them — so an
    /// older level is there only for a higher-derivative equation, which
    /// starts this request's classical stencil.
    const State& previous(std::size_t level = 0) const
    {
      return levels_.state(level);
    }

  private:
    static void require_first_derivative(int derivative)
    {
      if (derivative != 1)
        throw std::invalid_argument("atstep: this BDF request supports "
                                    "only the first derivative");
    }

    Derivative derivative_;
    const detail::Levels<State>& levels_;
    const State& derivative_history_;
    /// The difference `y - predictor` of the correction form, borrowed from
    /// the stepper.
    State& derivative_scratch_;
    const State& absolute_scale_;
  };

  /// What the caller's solve reports back.
  struct SolveOutcome {
    /// Whether the nonlinear solve converged. A step that did not is
    /// retried at a quarter of its size.
    bool converged = false;
    /// The nonlinear iterations the solve took, for the log.
    int iterations = 0;
  };

  /// An accepted step, as the log sees it.
  struct StepRecord {
    double time = 0.0;
    double step = 0.0;
    int order = 0;
    int iterations = 0;
    /// The error estimate of the step: held at most to one.
    double error = 0.0;
  };

  /// What a run reports as it goes.
  template <class State>
  struct RunHooks {
    /// Every accepted step, in order, with the state it left: the level the
    /// history now holds. The state is the stepper's own buffer and is only
    /// valid during the call.
    std::function<void(const StepRecord&, const State&)> on_step;
    /// Every stored value: the time, the state, and whether it was interpolated
    /// from the steps around it rather than being a level. The state is the
    /// stepper's own buffer and is only valid during the call — copy it, or
    /// evaluate the derivative first, before calling `interpolate` again.
    std::function<void(double, const State&, bool)> on_store;
  };

  /// What a run did.
  struct RunReport {
    int steps = 0;
    int rejected = 0;
    int failed_solves = 0;
    int stored = 0;
    int max_order = 1;
    double first_step = 0.0;
    double last_step = 0.0;
    double final_time = 0.0;
  };

  /// Advances a state in time with COMSOL's adaptive BDF method.
  ///
  /// The caller owns the equations: it is given the weights of the step and
  /// solves the system they define, with its own operator, its own linear
  /// solver and its own nonlinear iteration. The stepper owns everything
  /// else: the history, the predictor the solve starts from, the local error
  /// estimate, the step size and the order, and the placement of the steps
  /// relative to the output times. A step it takes is one the error test
  /// accepts; a step the solve or the test refuses is retried smaller, inside
  /// the call, without the caller seeing it.
  ///
  ///     atstep::Settings settings;              // COMSOL's defaults
  ///     settings.max_order = 5;
  ///     atstep::TimeStepper<Vector> stepper(settings,
  ///         [&](const atstep::StepRequest<Vector>& request, Vector& y) {
  ///             assemble_and_solve(request, y); // M y' + K y = f
  ///             return atstep::SolveOutcome {true, iterations};
  ///         });
  ///     stepper.start(0.0, y0, 10.0);
  ///     while (stepper.time() < stepper.end_time()) {
  ///         const auto step = stepper.advance();
  ///         const Vector& y = stepper.state(); // the solution at step.time
  ///     }
  ///
  /// `run` is the same loop for a whole schedule of output times.
  ///
  /// `y` comes in holding the predictor of the step — the extrapolation of
  /// the history, which is a good initial guess — and must leave holding the
  /// solution of the step.
  ///
  /// An integrator is neither copyable nor movable: a run is one object's
  /// state from `start` to its last step, and `start` begins it again. The
  /// states a run hands out are its own buffers, valid until the next step.
  template <class State, class Norm = detail::WeightedRmsNorm>
    requires detail::StateVector<State>
  class TimeStepper {
  public:
    using Solve = std::function<SolveOutcome(const StepRequest<State>&, State&)>;

    TimeStepper(const TimeStepper&) = delete;
    TimeStepper& operator=(const TimeStepper&) = delete;
    TimeStepper(TimeStepper&&) = delete;
    TimeStepper& operator=(TimeStepper&&) = delete;

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one step (see the class comment).
    /// @param[in] norm the norm the error of a step is measured in. The
    /// default is COMSOL's weighted RMS, whose tolerances are taken from
    /// `settings`.
    explicit TimeStepper(const Settings& settings,
      Solve solve,
      Norm norm = Norm {})
      : settings_(settings)
      , solve_(std::move(solve))
      , norm_(std::move(norm))
    {
      validate(settings_);
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>) {
        norm_.relative = settings_.relative_tolerance;
        norm_.absolute = settings_.relative_tolerance * settings_.absolute_factor;
        norm_.scaled = settings_.scale_absolute_tolerance;
        norm_.scale = settings_.absolute_scale;
      }
      // The formulation does not change, so the levels a run keeps are known
      // before it starts: the newest alone, or the whole stencil's.
      levels_.reserve(fixed_leading() ? 1U : static_cast<std::size_t>(max_levels));
    }

    /// Prepare the run: the state at `time`, integrated to `end_time`.
    ///
    /// The state is advanced by two artificial backward-Euler steps of
    /// `init_step_fraction * first_step_fraction` of the span (of the
    /// manual step in manual mode) — COMSOL's measured consistent
    /// initialization, which reconciles the initial values with the
    /// constraints of the equations and gives the time derivative the first
    /// step's size is read from. That step is a solve of the caller's, and
    /// its result is the level the run starts at.
    void start(double time, const State& initial, double end_time)
    {
      if (not std::isfinite(time) or not std::isfinite(end_time) or not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
      start_time_ = time;
      end_time_ = end_time;
      const double span = end_time - time;
      controller_ = detail::StepController(settings_, span);
      started_ = false;
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>)
        norm_.step_index = 0;
      has_previous_correction_ = false;
      previous_order_ = 0;
      last_order_ = 1;
      work_.derivative = initial;
      detail::axpy(work_.derivative, -1.0, initial);
      probe_ = 1.0;

      levels_.reset(time, initial);
      work_.scale_reference = initial;

      // An equation of the second order in time has no derivative inside
      // the state, so there is nothing for the initialization step to
      // measure: the run starts from the span rule and a constant
      // predictor, and the equation should be given the derivative of its
      // state, or written as a first-order system in `(y, y')`, which is
      // how COMSOL states it.
      const double probe = settings_.derivative_orders >= 2
        ? 0.0
        : (settings_.step_mode == StepMode::manual
              ? settings_.manual_step
              : span)
          * settings_.first_step_fraction * settings_.init_step_fraction;
      if (probe <= 0.0) {
        controller_.begin(controller_.first_step(0.0));
        started_ = true;
        return;
      }

      const detail::StepFormula probe_formula(1,
        std::span<const double>(&probe, 1),
        settings_.derivative_orders);
      // COMSOL's backward-Euler initialization advances twice. Its stored
      // initial state is the second result, labeled with the start time.
      work_.solution = initial;
      for (int j = 1; j <= 2; ++j) {
        classical_derivative_history(probe_formula, work_.derivative_history);
        StepRequest<State> request {
          time + static_cast<double>(j) * probe,
          probe,
          1,
          typename StepRequest<State>::ClassicalDerivative {
            probe_formula},
          levels_,
          work_.derivative_history,
          work_.scratch,
          levels_.state(0),
          0,
          true};
        const SolveOutcome outcome = solve_(request, work_.solution);
        if (not outcome.converged)
          throw std::runtime_error(
            "atstep: the initialization step of the run at t = " + std::to_string(time) + " did not converge");
        work_.scratch = work_.solution;
        detail::axpy(work_.scratch, -1.0, levels_.state(0));
        levels_.reset(request.time, work_.solution);
      }

      // The time derivative the first step is read from and predicted by:
      // the increment of the initialization step over its size.
      const double derivative_norm = norm_(work_.scratch, work_.solution) / probe;

      // The initialization step is not a level of the scheme: its
      // interval is a thousandth of the first step, and a polynomial
      // through it would be meaningless. It leaves the state it produced
      // and the derivative it measured, which is what the first step is
      // predicted — and sized — from. IDA starts the same way, from a
      // state and a time derivative and no history at all.
      work_.derivative = work_.scratch;
      probe_ = probe;
      levels_.reset(time + 2.0 * probe, work_.solution);
      work_.scale_reference = work_.solution;
      controller_.begin(controller_.first_step(derivative_norm));
      if (fixed_leading()) {
        work_.last_derivative = work_.derivative;
        detail::axpy(work_.last_derivative, 1.0 / probe - 1.0, work_.derivative);
        // Zero only the columns the configured order can reach: the scheme
        // reads `phi[0..max_order+1]`, and a state allocated here is a
        // state held for the whole run.
        for (int j = 0; j <= settings_.max_order + 1; ++j) {
          State& column = table_.phi[static_cast<std::size_t>(j)];
          column = work_.solution;
          detail::axpy(column, -1.0, work_.solution);
        }
        table_.phi[0] = work_.solution;
        detail::axpy(
          table_.phi[1], controller_.step_size() / probe, work_.derivative);
        table_.psi.fill(0.0);
        table_.psi[0] = controller_.step_size();
      }
      started_ = true;
    }

    const Settings& settings() const noexcept
    {
      return settings_;
    }
    double start_time() const noexcept
    {
      return start_time_;
    }
    double end_time() const noexcept
    {
      return end_time_;
    }

    /// The time the run has reached: the newest accepted level's.
    double time() const
    {
      require_started();
      return levels_.time(0);
    }

    /// The current solution: the newest accepted level's state. It is the
    /// stepper's own buffer, overwritten by the next step.
    const State& state() const
    {
      require_started();
      return levels_.state(0);
    }

    /// The number of accepted levels the run keeps: the newest alone for a
    /// first-order equation, the whole stencil's for a higher-derivative one.
    std::size_t history_size() const noexcept
    {
      return levels_.size();
    }

    /// The step size and order of the next step.
    double step_size() const
    {
      return controller_.step_size();
    }
    int order() const
    {
      return controller_.order();
    }

    /// Take one step the error test accepts, from the newest level, and
    /// report it.
    ///
    /// It proposes the step the controller asks for, retries it smaller when
    /// the solve does not converge or the error test fails, and accepts it.
    /// It knows nothing of the output times: `run` is what cuts a step short
    /// to land on one.
    ///
    /// @throws std::runtime_error when the step has fallen to the smallest
    /// one the controller proposes and still does not converge or still
    /// misses the tolerance.
    StepRecord advance()
    {
      require_started();
      return take_step(nullptr, nullptr);
    }

    /// Integrate from the start time to the end time, and report the values
    /// the store mode asks for.
    ///
    /// `output_times` is increasing, starts at the start time and ends at
    /// the end time the run was started with.
    RunReport run(std::span<const double> output_times,
      RunHooks<State> hooks = {})
    {
      require_started();
      detail::StepSchedule schedule(
        settings_, output_times, start_time_, end_time_);

      RunReport report;
      double time = levels_.time(0);

      // COMSOL labels the state left by consistent initialization with
      // the requested start time; the initialization solves are not
      // accepted steps.
      store(start_time_, levels_.state(0), false, hooks, report);

      while (time < end_time_) {
        const StepRecord record = take_step(&schedule, &report);
        time = record.time;
        ++report.steps;
        report.max_order = std::max(report.max_order, record.order);
        if (report.steps == 1)
          report.first_step = record.step;
        report.last_step = record.step;
        report.final_time = time;
        if (hooks.on_step)
          hooks.on_step(record, levels_.state(0));
        for (const auto& stored : schedule.stores(time))
          store(stored.time,
            stored.interpolated ? interpolate(stored.time)
                                : levels_.state(0),
            stored.interpolated,
            hooks,
            report);
      }
      return report;
    }

    /// The dense-output value (derivative=0) or first derivative
    /// (derivative=1). At the newest level the stored derivative is the
    /// actual solve derivative; inside a step it is the derivative of the
    /// BDF dense-output polynomial. The result is kept in the stepper:
    /// another call can overwrite it.
    const State& interpolate(double time, int derivative = 0)
    {
      require_started();
      if (derivative < 0 or derivative > 1)
        throw std::invalid_argument(
          "atstep: interpolation supports value or first derivative");
      if (fixed_leading()) {
        if (derivative == 1 and time == levels_.time(0))
          return work_.last_derivative;
        const double delta = time - levels_.time(0);
        double coefficient = 1.0;
        double slope = 0.0;
        double gamma = delta / table_.psi[0];
        work_.interpolated = table_.phi[0];
        if (derivative == 1)
          detail::axpy(work_.interpolated, -1.0, table_.phi[0]);
        for (int j = 1; j <= last_order_; ++j) {
          slope = slope * gamma + coefficient / table_.psi[static_cast<std::size_t>(j - 1)];
          coefficient *= gamma;
          detail::axpy(work_.interpolated,
            derivative == 0 ? coefficient : slope,
            table_.phi[static_cast<std::size_t>(j)]);
          gamma = (delta + table_.psi[static_cast<std::size_t>(j - 1)]) / table_.psi[static_cast<std::size_t>(j)];
        }
        return work_.interpolated;
      }
      const std::size_t count = std::min<std::size_t>(
        static_cast<std::size_t>(last_order_) + 1, levels_.size());
      std::array<double, max_levels> nodes {};
      for (std::size_t j = 0; j < count; ++j)
        nodes[j] = levels_.time(j) - time;
      const auto weights = detail::lagrange_derivative(
        std::span<const double>(nodes.data(), count),
        0.0,
        derivative);
      combine(work_.interpolated,
        std::span<const double>(weights.data(), count));
      return work_.interpolated;
    }

  private:
    /// One attempt: the step, its order, the coefficients the trial reads the
    /// accepted table through, and the estimates it produced. It borrows the
    /// workspace's states, so it is created and consumed inside one step and
    /// never outlives it.
    struct Trial {
      double step = 0.0;
      double time = 0.0;
      int order = 1;
      int iterations = 0;
      bool solved = false;
      bool passed = false;
      detail::StepErrors errors;
      /// The weight `y'` carries on the state being solved for, the factor the
      /// error test is held to, the columns' factors at this step size, the
      /// distances back to the levels, and IDA's `sigma`.
      double cj = 0.0;
      double ck = 0.5;
      std::array<double, max_levels> beta {};
      std::array<double, max_levels> psi {};
      std::array<double, max_levels> sigma {};
      /// The classical stencil, which a first-order run does not build.
      detail::StepFormula formula;
    };

    /// The columns of the first-order table: `phi[j]` is the j-th scaled
    /// difference of the accepted steps, and `phi[0]` the level the predictor
    /// extrapolates from — which the solve's own state need not equal in its
    /// last bits. Only `accept` writes a column of it beyond the first step; an
    /// attempt reads it. `psi[j]` is the distance back to the level `j`.
    struct Table {
      std::array<State, max_levels> phi {};
      std::array<double, max_levels> psi {};
    };

    /// The states of a step. `scratch` is shared by the derivative's
    /// correction form and the error estimates, which never run at the same
    /// time.
    struct Workspace {
      State derivative {};
      State predictor {};
      State solution {};
      State correction {};
      State previous_correction {};
      State scratch {};
      State interpolated {};
      State derivative_history {};
      State derivative_predictor {};
      State last_derivative {};
      State scale_reference {};
    };

    /// Take one step the error test accepts, retrying it smaller while the
    /// solve or the test refuses it.
    ///
    /// @param[in] schedule the output times the step is placed against, or
    /// nullptr for a step of the controller's own.
    /// @param[in,out] report the counters of the run, or nullptr.
    StepRecord take_step(detail::StepSchedule* schedule, RunReport* report)
    {
      const double from = levels_.time(0);
      for (;;) {
        double step = controller_.step_size();
        std::optional<double> boundary;
        if (schedule) {
          const detail::StepSchedule::Placement placement
            = schedule->place(step, from);
          step = placement.step;
          boundary = placement.boundary;
        }
        if (not(step > 0.0))
          throw std::runtime_error(
            "atstep: no step is left at t = " + std::to_string(from));

        Trial trial = attempt(step);
        if (not trial.solved) {
          if (report)
            ++report->failed_solves;
          const double tried = trial.step;
          controller_.reject_convergence(tried);
          if (tried <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from) + " did not converge at the smallest step "
              + std::to_string(tried));
          continue;
        }
        if (not trial.passed) {
          if (report)
            ++report->rejected;
          const double error = trial.errors.test;
          controller_.reject(trial.errors, trial.step);
          if (controller_.step_size() <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller_.step_size()) + " and the local error is still " + std::to_string(error) + "; the tolerance cannot be met on this interval");
          continue;
        }
        if (boundary)
          trial.time = *boundary;
        const StepRecord record {trial.time,
          trial.step,
          trial.order,
          trial.iterations,
          trial.errors.test};
        accept(trial);
        return record;
      }
    }

    /// Solve the step of size `step` and estimate its error. The accepted
    /// levels and table are read, not written, apart from the first step's
    /// column, which the attempt states: only `accept` commits a step.
    Trial attempt(double step)
    {
      const double next_time = levels_.time(0) + step;
      if (not std::isfinite(step) or not(step > 0.0)
        or not std::isfinite(next_time) or not(next_time > levels_.time(0)))
        throw std::invalid_argument("atstep: step must be finite, positive and advance time");

      // Both formulations use the same history bound on startup order: an
      // order-k step reads k levels, and the run has one level per accepted
      // step on top of the one it started from.
      const int available = std::min(
        controller_.steps() + 1, max_levels);
      const int order = std::min(controller_.order(), std::max(1, available - 1));
      controller_.begin_attempt(step, order);

      Trial trial;
      trial.step = step;
      trial.time = next_time;
      trial.order = order;

      auto derivative = prepare(trial);
      work_.solution = work_.predictor;
      StepRequest<State> request {
        next_time,
        step,
        order,
        std::move(derivative),
        levels_,
        work_.derivative_history,
        work_.scratch,
        work_.scale_reference,
        static_cast<std::size_t>(controller_.steps())};
      const SolveOutcome outcome = solve_(request, work_.solution);
      trial.iterations = outcome.iterations;
      if (not outcome.converged)
        return trial;
      trial.solved = true;

      if (fixed_leading() or controller_.error_controlled()) {
        work_.correction = work_.solution;
        detail::axpy(work_.correction, -1.0, work_.predictor);
      }
      if (not controller_.error_controlled()) {
        trial.passed = true;
        return trial;
      }
      trial.errors = estimate(trial);
      trial.passed = trial.errors.passed();
      return trial;
    }

    /// Commit the step: the level, the table, and the controller's next size
    /// and order.
    void accept(const Trial& trial)
    {
      if (fixed_leading()) {
        work_.last_derivative = work_.derivative_predictor;
        detail::axpy(work_.last_derivative, trial.cj, work_.correction);
      }
      levels_.push(trial.time, work_.solution);
      last_order_ = trial.order;
      if (fixed_leading()) {
        const auto k = static_cast<std::size_t>(trial.order);
        // The factors the trial read are the factors the table keeps: the
        // rescaled column of the order it is at is a scalar multiple of the
        // accepted one, so a dropped attempt has nothing to put back.
        for (std::size_t j = 1; j <= k; ++j)
          detail::axpy(
            table_.phi[j], trial.beta[j] - 1.0, table_.phi[j]);
        table_.phi[k + 1] = work_.correction;
        detail::axpy(table_.phi[k], 1.0, work_.correction);
        for (std::size_t j = k; j > 0; --j)
          detail::axpy(table_.phi[j - 1], 1.0, table_.phi[j]);
        table_.psi = trial.psi;
      }
      controller_.accept(trial.errors, trial.step);
      work_.scale_reference = work_.solution;
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>)
        norm_.step_index = static_cast<std::size_t>(controller_.steps());
      if (not fixed_leading() and controller_.error_controlled()) {
        work_.previous_correction = work_.correction;
        previous_order_ = trial.order;
        has_previous_correction_ = true;
      }
    }

    bool fixed_leading() const
    {
      return settings_.derivative_orders == 1;
    }

    void require_started() const
    {
      if (not started_)
        throw std::logic_error(
          "atstep: the run has not been started (see start)");
    }

    /// The known term of the first derivative of the classical stencil: the
    /// combination of the accepted levels, without the level being solved for.
    void classical_derivative_history(const detail::StepFormula& formula,
      State& out) const
    {
      const auto weights = formula.weights(1);
      out = levels_.state(0);
      detail::axpy(out, weights[1] - 1.0, levels_.state(0));
      for (std::size_t j = 2; j < weights.size(); ++j)
        detail::axpy(out, weights[j], levels_.state(j - 1));
    }

    /// Fill the trial's coefficients, the predictor the solve starts from, the
    /// derivative's predictor and the known term of the derivative, and return
    /// the derivative formulation the solve reads.
    typename StepRequest<State>::Derivative prepare(Trial& trial)
    {
      if (fixed_leading()) {
        prepare_fixed_leading(trial);
        return typename StepRequest<State>::BdfDerivative {
          work_.predictor, work_.derivative_predictor, trial.cj};
      }
      const std::size_t count = levels_.size();
      std::array<double, max_levels> steps {};
      steps[0] = trial.step;
      for (std::size_t j = 1; j < count; ++j)
        steps[j] = levels_.time(j - 1) - levels_.time(j);
      trial.formula = detail::StepFormula(trial.order,
        std::span<const double>(steps.data(), count),
        settings_.derivative_orders);
      predict(trial);
      classical_derivative_history(trial.formula, work_.derivative_history);
      return typename StepRequest<State>::ClassicalDerivative {
        trial.formula};
    }

    /// The coefficients of a first-order trial, from the accepted table.
    void prepare_fixed_leading(Trial& trial)
    {
      trial.psi = table_.psi;
      if (controller_.steps() == 0) {
        // The first step of a run has no accepted column behind it: the one the
        // scheme reads is the derivative the initialization measured, at the
        // size this attempt takes. It is a function of the attempt, so the
        // attempt states it, and a retry states it again.
        trial.psi[0] = trial.step;
        table_.phi[1] = work_.derivative;
        detail::axpy(
          table_.phi[1], trial.step / probe_ - 1.0, work_.derivative);
      }
      // Beyond that column the trial reads the accepted table and writes none
      // of it: the rescaled column of the order it is at is `beta_j phi[j]`,
      // which is a scalar multiple of the accepted column, so the factors
      // `beta_j` are all it keeps. A refused attempt therefore has nothing to
      // put back, and an accepted one writes the factors in.
      work_.predictor = table_.phi[0];
      work_.derivative_predictor = table_.phi[0];
      detail::axpy(work_.derivative_predictor, -1.0, table_.phi[0]);
      double beta = 1.0;
      double distance = trial.step;
      double alpha = 1.0;
      double gamma = 0.0;
      double alpha0 = -1.0;
      double harmonic = 0.0;
      trial.sigma[0] = 1.0;
      for (int i = 1; i <= trial.order; ++i) {
        const auto j = static_cast<std::size_t>(i);
        const double old_distance = trial.psi[j - 1];
        trial.psi[j - 1] = distance;
        beta *= distance / old_distance;
        distance = old_distance + trial.step;
        gamma += alpha / trial.step;
        alpha = trial.step / distance;
        if (i < trial.order)
          alpha0 -= alpha;
        harmonic += 1.0 / static_cast<double>(i);
        trial.sigma[j] = static_cast<double>(i) * trial.sigma[j - 1] * alpha;
        trial.beta[j] = beta;
        detail::axpy(work_.predictor, beta, table_.phi[j]);
        detail::axpy(work_.derivative_predictor, gamma * beta, table_.phi[j]);
      }
      trial.psi[static_cast<std::size_t>(trial.order)] = distance;
      trial.ck = std::max(alpha, std::abs(alpha - harmonic - alpha0));
      trial.cj = harmonic / trial.step;
      work_.derivative_history = work_.derivative_predictor;
      detail::axpy(work_.derivative_history, -trial.cj, work_.predictor);
    }

    /// The extrapolation of the history to the new time, in `predictor`. A
    /// history of a single level is the initial state, and the step from it
    /// is predicted by the derivative the initialization step measured.
    void predict(const Trial& trial)
    {
      if (levels_.size() == 1) {
        // `derivative` is the increment of the initialization step, so the
        // extrapolation of the step being taken is `step / probe` of it.
        work_.predictor = levels_.state(0);
        detail::axpy(work_.predictor, trial.step / probe_, work_.derivative);
        return;
      }
      combine(work_.predictor, trial.formula.predictor());
    }

    /// The error estimates of the trial, whose solved state is `solution` and
    /// whose predictor is `predictor`.
    detail::StepErrors estimate(const Trial& trial)
    {
      detail::StepErrors errors;
      // IDA fixes the weights from the last accepted state before
      // attempting a step; every estimate and retry uses those same
      // weights.
      const double enorm = error_norm(work_.correction);
      if (fixed_leading()) {
        const auto k = static_cast<std::size_t>(trial.order);
        errors.test = trial.ck * enorm;
        errors.at = trial.sigma[k] * enorm;
        work_.scratch = work_.correction;
        if (trial.order >= 2) {
          detail::axpy(work_.scratch, trial.beta[k], table_.phi[k]);
          errors.below = trial.sigma[k - 1] * error_norm(work_.scratch);
        }
        if (trial.order >= 3) {
          detail::axpy(
            work_.scratch, trial.beta[k - 1], table_.phi[k - 1]);
          errors.below2 = trial.sigma[k - 2] * error_norm(work_.scratch);
        }
        work_.scratch = work_.correction;
        detail::axpy(work_.scratch, -1.0, table_.phi[k + 1]);
        errors.above = error_norm(work_.scratch) / static_cast<double>(trial.order + 2);
        return errors;
      }
      errors.test = trial.formula.ck() * enorm;
      errors.at = trial.formula.sigma(trial.order) * enorm;

      // The orders below read the polynomial through one level fewer
      // each, and the estimate above compares this step's correction with
      // the one the previous step left.
      if (trial.order >= 2)
        errors.below = trial.formula.sigma(trial.order - 1) * extrapolation_error(trial.formula.extrapolation(trial.order));
      if (trial.order >= 3)
        errors.below2 = trial.formula.sigma(trial.order - 2) * extrapolation_error(trial.formula.extrapolation(trial.order - 1));
      if (has_previous_correction_ and previous_order_ == trial.order) {
        work_.scratch = work_.correction;
        detail::axpy(work_.scratch, -1.0, work_.previous_correction);
        errors.above = error_norm(work_.scratch) / static_cast<double>(trial.order + 2);
      }
      return errors;
    }

    /// The weighted norm of `solution - p`, `p` being the polynomial through
    /// the levels with `weights`.
    double extrapolation_error(std::span<const double> weights)
    {
      combine(work_.scratch, weights);
      detail::axpy(work_.scratch, -1.0, work_.solution);
      return error_norm(work_.scratch);
    }

    double error_norm(const State& error)
    {
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>)
        return norm_(error, levels_.state(0), work_.scale_reference);
      else
        return norm_(error, levels_.state(0));
    }

    /// `out = weights[0] y_0 + ... + weights[count - 1] y_{count - 1}`,
    /// from the levels of the history.
    void combine(State& out, std::span<const double> weights) const
    {
      const std::size_t count = std::min(weights.size(), levels_.size());
      out = levels_.state(0);
      detail::axpy(out, weights[0] - 1.0, levels_.state(0));
      for (std::size_t j = 1; j < count; ++j)
        detail::axpy(out, weights[j], levels_.state(j));
    }

    void store(double time,
      const State& state,
      bool interpolated,
      const RunHooks<State>& hooks,
      RunReport& report)
    {
      ++report.stored;
      work_.scale_reference = state;
      if (hooks.on_store)
        hooks.on_store(time, state, interpolated);
    }

    Settings settings_;
    Solve solve_;
    Norm norm_;

    detail::StepController controller_ {Settings {}, 1.0};
    detail::Levels<State> levels_;
    Table table_;
    Workspace work_;

    double start_time_ = 0.0;
    double end_time_ = 0.0;
    /// The size of the initialization step (see `start`).
    double probe_ = 1.0;
    bool started_ = false;

    /// The order of the last accepted step, which the dense output reads.
    int last_order_ = 1;
    /// The correction of the last accepted step and its order: the estimate at
    /// the order above compares against them.
    bool has_previous_correction_ = false;
    int previous_order_ = 0;
  };

} // namespace atstep
