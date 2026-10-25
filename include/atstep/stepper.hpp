// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

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

  template <class State, class Norm>
    requires detail::StateVector<State>
  class TimeStepper;

  /// Borrowed coefficients for one first-order solve.
  ///
  /// The equation is first order in time, and the step holds its derivative to
  /// the affine form
  ///
  ///     y' = cj*y + known,   cj = derivative_weight(),
  ///                          known = derivative_history(),
  ///
  /// so that `M y' + K y = f` at `request.time` is
  ///
  ///     (cj M + K) y_{n+1} = f(t_{n+1}) - M known.
  ///
  /// `derivative(y, out)` evaluates `predictor' + cj (y - predictor)` to
  /// reduce cancellation. The solve receives the predictor in its mutable
  /// state argument and replaces it with the solution.
  ///
  /// The request and its state references are valid only during the solve
  /// callback. Do not retain pointers or references to them.
  template <class State>
  class StepRequest {
  public:
    StepRequest(const StepRequest&) = delete;
    StepRequest& operator=(const StepRequest&) = delete;
    StepRequest(StepRequest&&) = delete;
    StepRequest& operator=(StepRequest&&) = delete;

    const double time;
    const double step;
    const int order;
    const std::size_t step_index;
    const bool initializing;

    const State& absolute_tolerance_reference() const noexcept
    {
      return absolute_scale_;
    }

    /// The weight `y'` carries on the state being solved for.
    double derivative_weight() const noexcept
    {
      return leading_;
    }

    /// The known term of the derivative: `y' = cj*y + known`.
    const State& derivative_history() const noexcept
    {
      return derivative_history_;
    }

    /// `out = y'`. Input and output may alias; output must not be a state
    /// borrowed from this request.
    void derivative(const State& y, State& out) const
    {
      derivative_scratch_ = y;
      detail::axpy(derivative_scratch_, -1.0, predictor_);
      out = derivative_predictor_;
      detail::axpy(out, leading_, derivative_scratch_);
    }

    /// The accepted state used for relative error weights.
    const State& previous() const noexcept
    {
      return previous_;
    }

  private:
    template <class S, class N>
      requires detail::StateVector<S>
    friend class TimeStepper;
    template <class>
    friend class SecondOrderRequest;

    StepRequest(double new_time,
      double step_size,
      int step_order,
      const State& predictor,
      const State& derivative_predictor,
      double leading,
      const State& previous,
      const State& derivative_history,
      State& derivative_scratch,
      const State& absolute_scale,
      std::size_t accepted_steps,
      bool initialization)
      : time(new_time)
      , step(step_size)
      , order(step_order)
      , step_index(accepted_steps)
      , initializing(initialization)
      , predictor_(predictor)
      , derivative_predictor_(derivative_predictor)
      , leading_(leading)
      , previous_(previous)
      , derivative_history_(derivative_history)
      , derivative_scratch_(derivative_scratch)
      , absolute_scale_(absolute_scale)
    {
    }

    const State& predictor_;
    const State& derivative_predictor_;
    double leading_ = 0.0;
    const State& previous_;
    const State& derivative_history_;
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
    /// stepper's own buffer and is only valid during the call: copy it to
    /// retain it.
    std::function<void(double, const State&, bool)> on_store;
  };

  /// What a run did.
  struct RunReport {
    int steps = 0;
    int rejected = 0;
    int failed_solves = 0;
    int stored = 0;
    int max_order = 1;
    double final_time = 0.0;
  };

  /// Adaptive BDF integration of a first-order state. The caller supplies the
  /// algebraic solve; the stepper owns history, work buffers, error control
  /// and output. `advance` retries internally until one step is accepted.
  /// `run` performs the same loop over a complete output schedule.
  ///
  /// The integrator is neither copyable nor movable. `start` resets it. The
  /// states a hook is handed are borrowed for the call; `state()` hands out
  /// the run's own buffer, while `interpolate` writes into a state of the
  /// caller's. See doc/design.md for the numerical formulation.
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
      configure_norm();
    }

    /// Prepare the run: the state at `time`, integrated to `end_time`.
    ///
    /// The state is advanced by two artificial backward-Euler steps, each
    /// `initialization.fraction` of the run's initial step: the span times
    /// `AdaptiveSteps::initial_step_fraction`, or the manual step times the
    /// library's default fraction. That is COMSOL's measured consistent
    /// initialization, which reconciles the initial values with the constraints
    /// of the equations and gives the time derivative the first step's size is
    /// read from. Each step is a solve of the caller's, and the result of the
    /// second is the level the run starts at.
    void start(double time, const State& initial, double end_time)
    {
      require_bounds(time, end_time);
      started_ = false;
      controller_.emplace(settings_, end_time - time);
      configure_norm();
      const auto* backward = std::get_if<BackwardEuler>(&settings_.initialization);
      if (not backward)
        throw std::invalid_argument(
          "atstep: this run takes its initial values as they are, so it needs the "
          "initial derivative: call start(time, initial, derivative, end_time)");
      const double span = end_time - time;

      // The probe is the initial step of the alternative the run steps in, times
      // the fraction of it the initialization takes. A manual run has no initial
      // step of its own, and is scaled by the library default, as it always was
      // (doc/evidence/stepping-sweep.md, the manual startup curve). A named
      // initial step is the one the caller gave, and the probe is a fraction of
      // it.
      const auto [probe_step, probe_fraction] = std::visit(
        [span](const auto& steps) -> std::pair<double, double> {
          using Alternative = std::decay_t<decltype(steps)>;
          if constexpr (std::is_same_v<Alternative, ManualSteps>)
            return {steps.step, default_initial_step_fraction};
          else if (steps.initial_step)
            return {*steps.initial_step, 1.0};
          else
            return {span, default_initial_step_fraction};
        },
        settings_.steps);
      const double probe = probe_step * probe_fraction * backward->fraction;
      const double cj = 1.0 / probe;
      const double initialized_time = time + 2.0 * probe;
      if (not std::isfinite(probe) or not(probe > 0.0)
        or not std::isfinite(cj) or not std::isfinite(initialized_time)
        or not(time + probe > time) or not(initialized_time > time + probe)
        or not(initialized_time < end_time))
        throw std::invalid_argument(
          "atstep: initialization steps must be finite and advance within the run");

      // COMSOL's backward-Euler initialization advances twice, and neither is
      // a step of the scheme: the level is a local until the second one, so a
      // run that cannot initialize accepts nothing.
      State level = initial;
      State increment {};
      for (int j = 1; j <= 2; ++j) {
        // A backward-Euler step of the level: the level itself is the
        // predictor, and its derivative is zero.
        work_.predictor = level;
        work_.derivative_predictor = work_.predictor;
        detail::axpy(work_.derivative_predictor, -1.0, work_.predictor);
        const SolveOutcome outcome = solve(
          time + static_cast<double>(j) * probe, probe, 1, cj, level, level, true);
        if (not outcome.converged)
          throw std::runtime_error(
            "atstep: the initialization step of the run at t = " + std::to_string(time) + " did not converge");
        increment = work_.solution;
        detail::axpy(increment, -1.0, level);
        level = work_.solution;
      }

      // The time derivative the first step is read from and predicted by:
      // the increment of the initialization step over its size.
      State derivative = increment;
      detail::axpy(derivative, cj - 1.0, increment);
      begin_run(time, level, derivative, norm_(increment, level) / probe,
        initialized_time, end_time);
    }

    /// Start from a state and a derivative the caller gives, without the
    /// consistent initialization: COMSOL's "Consistent initialization: Off".
    void start(double time, const State& initial, const State& derivative, double end_time)
    {
      require_bounds(time, end_time);
      if (std::holds_alternative<BackwardEuler>(settings_.initialization))
        throw std::invalid_argument(
          "atstep: this run initializes consistently, so it takes no derivative: "
          "call start(time, initial, end_time)");
      started_ = false;
      controller_.emplace(settings_, end_time - time);
      configure_norm();
      begin_run(time, initial, derivative, norm_(derivative, initial), time, end_time);
    }

    /// The time the run has reached: the newest accepted level's.
    double time() const
    {
      require_started();
      return accepted_.time;
    }

    /// The bounds of a run, and the fields of its error test.
    void require_bounds(double time, double end_time) const
    {
      if (not std::isfinite(time) or not std::isfinite(end_time) or not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
    }

    /// The tail both starts share: the level and the derivative the run begins
    /// at, the first step, and the difference table a level and a derivative
    /// imply.
    void begin_run(double time,
      const State& level,
      const State& derivative,
      double derivative_norm,
      double initialized_time,
      double end_time)
    {
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>) {
        if (not norm_.fields.empty()) {
          std::size_t dofs = 0;
          for (const std::size_t field : norm_.fields)
            dofs += field;
          if (dofs != level.size())
            throw std::invalid_argument(
              "atstep: the fields of the error test do not add up to the size of the state");
        }
      }
      accepted_.state = level;
      accepted_.derivative = derivative;
      accepted_.time = initialized_time;
      accepted_.order = 1;
      output_.reference = level;
      const double step = controller_->first_step(derivative_norm, end_time - initialized_time);
      controller_->begin(step, initialized_time);
      // Zero only the columns the configured order can reach: the scheme
      // reads `phi[0..max_order+1]`, and a state allocated here is a
      // state held for the whole run.
      for (int j = 0; j <= controller_->order_limit() + 1; ++j) {
        State& column = accepted_.phi[static_cast<std::size_t>(j)];
        column = level;
        detail::axpy(column, -1.0, level);
      }
      // The first difference the scheme reads is the initial derivative at the
      // size the first step proposes — IDA's `phi[2] = h yp_0`. A step
      // rescaled to another size reads it through its own factor.
      accepted_.phi[0] = level;
      detail::axpy(accepted_.phi[1], controller_->step_size(), derivative);
      accepted_.psi.fill(0.0);
      accepted_.psi[0] = controller_->step_size();
      start_time_ = time;
      end_time_ = end_time;
      started_ = true;
    }

    /// The current solution: the newest accepted level's state. It is the
    /// stepper's own buffer, overwritten by the next step.
    const State& state() const
    {
      require_started();
      return accepted_.state;
    }

    /// The step size and order of the next step.
    double step_size() const
    {
      require_started();
      return controller_->step_size();
    }
    int order() const
    {
      require_started();
      return controller_->order();
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
    /// Call after `start`, before any `advance` or successful `run`.
    /// `output_times` increases from the start time to the end time.
    RunReport run(std::span<const double> output_times,
      RunHooks<State> hooks = {})
    {
      require_started();
      if (controller_->steps() != 0)
        throw std::logic_error("atstep: run requires a fresh start");
      detail::StepSchedule schedule(
        settings_, output_times, start_time_, end_time_);

      RunReport report;
      double time = accepted_.time;
      report.final_time = time;

      // COMSOL labels the state left by consistent initialization with
      // the requested start time; the initialization solves are not
      // accepted steps.
      store(start_time_, accepted_.state, false, hooks, report);

      while (time < end_time_) {
        const StepRecord record = take_step(&schedule, &report);
        time = record.time;
        ++report.steps;
        report.max_order = std::max(report.max_order, record.order);
        report.final_time = time;
        if (hooks.on_step)
          hooks.on_step(record, accepted_.state);
        schedule.for_each_store(time, static_cast<std::size_t>(report.steps),
          [&](double output_time, bool interpolated) {
            if (not interpolated) {
              store(output_time, accepted_.state, false, hooks, report);
              return;
            }
            interpolate(output_time, work_.scratch);
            store(output_time, work_.scratch, true, hooks, report);
          });
      }
      return report;
    }

    /// Write the dense-output value (derivative=0) or first derivative
    /// (derivative=1) at `time` into `out`. At the newest level the stored
    /// derivative is the actual solve derivative; inside a step it is the
    /// derivative of the BDF dense-output polynomial. `out` is the caller's,
    /// and may be the same state on every call.
    void interpolate(double time, State& out, int derivative = 0)
    {
      require_started();
      if (derivative < 0 or derivative > 1)
        throw std::invalid_argument(
          "atstep: interpolation supports value or first derivative");
      if (derivative == 1 and time == accepted_.time) {
        out = accepted_.derivative;
        return;
      }
      const double delta = time - accepted_.time;
      double coefficient = 1.0;
      double slope = 0.0;
      double gamma = delta / accepted_.psi[0];
      out = accepted_.phi[0];
      if (derivative == 1)
        detail::axpy(out, -1.0, accepted_.phi[0]);
      for (int j = 1; j <= accepted_.order; ++j) {
        slope = slope * gamma + coefficient / accepted_.psi[static_cast<std::size_t>(j - 1)];
        coefficient *= gamma;
        detail::axpy(out,
          derivative == 0 ? coefficient : slope,
          accepted_.phi[static_cast<std::size_t>(j)]);
        if (j < accepted_.order)
          gamma = (delta + accepted_.psi[static_cast<std::size_t>(j - 1)]) / accepted_.psi[static_cast<std::size_t>(j)];
      }
    }

  private:
    /// One attempt: the step, its order, the coefficients it reads the accepted
    /// table through, and the estimates it produced. It is created and consumed
    /// inside one step of `take_step` and never outlives it: a refused attempt
    /// is a `continue` of that loop, not a state the stepper keeps.
    struct Trial {
      double step = 0.0;
      double time = 0.0;
      int order = 1;
      int iterations = 0;
      detail::StepErrors errors;
      /// The weight `y'` carries on the state being solved for, the factor the
      /// error test is held to, the columns' factors at this step size, the
      /// distances back to the levels, and IDA's `sigma`.
      double cj = 0.0;
      double ck = 0.5;
      std::array<double, max_levels> beta {};
      std::array<double, max_levels> psi {};
      std::array<double, max_levels> sigma {};
    };

    /// Accepted state and scaled differences. A trial only reads them.
    /// `phi[0]` is updated by accumulation and can differ from the solved state
    /// in its last bits; keeping both preserves the BDF arithmetic and weights.
    struct History {
      State state {};
      State derivative {};
      double time = 0.0;
      int order = 1;
      std::array<State, max_levels> phi {};
      std::array<double, max_levels> psi {};
    };

    /// Reused across attempts. During solve, scratch holds the known term and
    /// correction is the request's derivative workspace. After solve returns,
    /// correction holds the solved-state increment; scratch serves estimates,
    /// the accepted derivative, and the dense output a store reports. No
    /// request survives that phase boundary.
    struct Work {
      State predictor {};
      State derivative_predictor {};
      State solution {};
      State correction {};
      State scratch {};
    };

    /// The state of the output: the value the last store reported, which the
    /// scaled absolute tolerance is measured against. The dense output itself
    /// is the caller's, so the stepper keeps no buffer of it.
    struct Output {
      State reference {};
    };

    /// Construct the request only for the duration of the callback. This is
    /// also the buffer-lifetime boundary between solving and error estimation.
    SolveOutcome solve(double time, double step, int order, double cj,
      const State& previous, const State& reference, bool initializing = false)
    {
      work_.scratch = work_.derivative_predictor;
      detail::axpy(work_.scratch, -cj, work_.predictor);
      work_.solution = work_.predictor;
      const StepRequest<State> request {time, step, order,
        work_.predictor, work_.derivative_predictor, cj, previous,
        work_.scratch, work_.correction, reference,
        initializing ? 0 : static_cast<std::size_t>(controller_->steps()),
        initializing};
      return solve_(request, work_.solution);
    }

    /// Take one step the error test accepts, retrying it smaller while the
    /// solve or the test refuses it.
    ///
    /// The accepted history is written at the end of a pass and nowhere else, so
    /// a refusal is a `continue` that leaves the run as the accepted steps left
    /// it.
    ///
    /// @param[in] schedule the output times the step is placed against, or
    /// nullptr for a step of the controller's own.
    /// @param[in,out] report the counters of the run, or nullptr.
    StepRecord take_step(detail::StepSchedule* schedule, RunReport* report)
    {
      const double from = accepted_.time;
      auto& controller = *controller_;
      for (;;) {
        double step = controller.step_size();
        std::optional<double> boundary;
        if (schedule) {
          const detail::StepSchedule::Landing placement
            = schedule->place(step, from);
          step = placement.step;
          boundary = placement.boundary;
        }
        if (not(step > 0.0))
          throw std::runtime_error(
            "atstep: no step is left at t = " + std::to_string(from));
        const double next_time = from + step;
        if (not std::isfinite(step) or not std::isfinite(next_time) or not(next_time > from))
          throw std::invalid_argument(
            "atstep: step must be finite, positive and advance time");

        // An order-k step reads k levels, and the run has one level per accepted
        // step on top of the one it started from.
        Trial trial;
        trial.step = step;
        trial.time = next_time;
        trial.order = std::min(controller.order(), std::max(1, controller.steps()));
        controller.begin_attempt(step, trial.order);
        prepare(trial);

        const SolveOutcome outcome = solve(next_time, step, trial.order,
          trial.cj, accepted_.state, output_.reference);
        trial.iterations = outcome.iterations;
        if (not outcome.converged) {
          if (report)
            ++report->failed_solves;
          controller.reject_convergence(step);
          if (step <= controller.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from) + " did not converge at the smallest step "
              + std::to_string(step));
          continue;
        }

        work_.correction = work_.solution;
        detail::axpy(work_.correction, -1.0, work_.predictor);
        if (controller.error_controlled()) {
          trial.errors = estimate(trial);
          if (not trial.errors.passed()) {
            if (report)
              ++report->rejected;
            controller.reject(trial.errors, step);
            if (controller.step_size() <= controller.smallest_step())
              throw std::runtime_error(
                "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller.step_size()) + " and the local error is still " + std::to_string(trial.errors.test) + "; the tolerance cannot be met on this interval");
            continue;
          }
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

    /// Update history in place, reusing its storage. If a State operation
    /// throws, the run stays unreadable until start rebuilds it. Ordinary
    /// rejected attempts never enter here and leave accepted history intact.
    void accept(const Trial& trial)
    {
      // The derivative the level pairs with: the predictor's derivative plus the
      // correction's, which is the formula's derivative at the new level.
      work_.scratch = work_.derivative_predictor;
      detail::axpy(work_.scratch, trial.cj, work_.correction);
      started_ = false;
      accepted_.state = work_.solution;
      accepted_.derivative = work_.scratch;
      accepted_.time = trial.time;
      accepted_.order = trial.order;
      const auto k = static_cast<std::size_t>(trial.order);
      // The factors the trial read are the factors the table keeps: the
      // rescaled column of the order it is at is a scalar multiple of the
      // accepted one, so a dropped attempt has nothing to put back.
      for (std::size_t j = 1; j <= k; ++j)
        detail::axpy(
          accepted_.phi[j], trial.beta[j] - 1.0, accepted_.phi[j]);
      accepted_.phi[k + 1] = work_.correction;
      detail::axpy(accepted_.phi[k], 1.0, work_.correction);
      for (std::size_t j = k; j > 0; --j)
        detail::axpy(accepted_.phi[j - 1], 1.0, accepted_.phi[j]);
      accepted_.psi = trial.psi;
      controller_->accept(trial.errors, trial.step);
      output_.reference = work_.solution;
      configure_norm();
      started_ = true;
    }

    void require_started() const
    {
      if (not started_)
        throw std::logic_error(
          "atstep: no valid run; call start before accessing or advancing it");
    }

    /// Fill the trial coefficients and predictors without changing history.
    void prepare(Trial& trial)
    {
      trial.psi = accepted_.psi;
      work_.predictor = accepted_.phi[0];
      work_.derivative_predictor = accepted_.phi[0];
      detail::axpy(work_.derivative_predictor, -1.0, accepted_.phi[0]);
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
        detail::axpy(work_.predictor, beta, accepted_.phi[j]);
        detail::axpy(work_.derivative_predictor, gamma * beta, accepted_.phi[j]);
      }
      trial.psi[static_cast<std::size_t>(trial.order)] = distance;
      trial.ck = std::max(alpha, std::abs(alpha - harmonic - alpha0));
      trial.cj = harmonic / trial.step;
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
      const auto k = static_cast<std::size_t>(trial.order);
      errors.test = trial.ck * enorm;
      errors.at = trial.sigma[k] * enorm;
      work_.scratch = work_.correction;
      if (trial.order >= 2) {
        detail::axpy(work_.scratch, trial.beta[k], accepted_.phi[k]);
        errors.below = trial.sigma[k - 1] * error_norm(work_.scratch);
      }
      if (trial.order >= 3) {
        detail::axpy(
          work_.scratch, trial.beta[k - 1], accepted_.phi[k - 1]);
        errors.below2 = trial.sigma[k - 2] * error_norm(work_.scratch);
      }
      work_.scratch = work_.correction;
      detail::axpy(work_.scratch, -1.0, accepted_.phi[k + 1]);
      errors.above = error_norm(work_.scratch) / static_cast<double>(trial.order + 2);
      return errors;
    }

    /// Point the norm at the run: the tolerances of its settings, copied once,
    /// and the steps it has accepted, which its scaled absolute tolerance
    /// follows. A norm of the caller's holds neither, and this is nothing for
    /// it.
    void configure_norm()
    {
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>) {
        if (not controller_)
          norm_ = detail::WeightedRmsNorm::from(settings_.error);
        norm_.step_index = controller_ ? static_cast<std::size_t>(controller_->steps()) : 0;
      }
    }

    /// The norm of the run's error, in the shape the run's norm has: COMSOL's
    /// scaled absolute tolerance follows the solution returned to the caller,
    /// and a norm of the caller's needs only the accepted state.
    double error_norm(const State& error)
    {
      if constexpr (requires { norm_(error, accepted_.state, output_.reference); })
        return norm_(error, accepted_.state, output_.reference);
      else
        return norm_(error, accepted_.state);
    }

    void store(double time,
      const State& state,
      bool interpolated,
      const RunHooks<State>& hooks,
      RunReport& report)
    {
      ++report.stored;
      output_.reference = state;
      if (hooks.on_store)
        hooks.on_store(time, state, interpolated);
    }

    Settings settings_;
    Solve solve_;
    Norm norm_;

    std::optional<detail::StepController> controller_;
    History accepted_;
    Work work_;
    Output output_;

    double start_time_ = 0.0;
    double end_time_ = 0.0;
    bool started_ = false;
  };

} // namespace atstep
