// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "controller.hpp"
#include "formula.hpp"
#include "settings.hpp"
#include "state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace atstep {

  /// The levels a run has passed through, newest first.
  ///
  /// The classical scheme reads the last `order + 1` of them, the predictor one
  /// more, and the error estimate at the order two below three fewer, so a run
  /// of a higher-derivative equation keeps `max_levels` of them. A first-order
  /// equation reads the difference table instead and keeps the newest level
  /// alone: the levels behind it are what the fixed-leading table holds.
  template <class State>
  class History {
  public:
    void reset(double time, const State& state)
    {
      std::deque<Level> replacement;
      replacement.push_front(Level {time, state});
      levels_.swap(replacement);
    }

    void push(double time, const State& state)
    {
      levels_.push_front(Level {time, state});
    }

    void pop()
    {
      if (not levels_.empty())
        levels_.pop_front();
    }

    std::size_t size() const noexcept
    {
      return levels_.size();
    }
    bool empty() const noexcept
    {
      return levels_.empty();
    }

    /// The time of the level `level`, 0 being the newest.
    double time(std::size_t level) const
    {
      return levels_.at(level).time;
    }

    /// The state of the level `level`, 0 being the newest.
    const State& state(std::size_t level) const
    {
      return levels_.at(level).state;
    }

    /// The times, newest first.
    std::vector<double> times() const
    {
      std::vector<double> result;
      result.reserve(size());
      for (const auto& level : levels_)
        result.push_back(level.time);
      return result;
    }

    /// Forget the levels past the newest `levels`.
    void keep(std::size_t levels)
    {
      while (size() > levels)
        levels_.pop_back();
    }

  private:
    struct Level {
      double time;
      State state;
    };
    std::deque<Level> levels_;
  };

  /// What a step's solve is told about the step.
  ///
  /// A first-order equation uses the affine derivative `cj*y + known`, with
  /// `cj = derivative_weight()` and `known = derivative_history()`. Thus
  /// `M y' + K y = f` at `request.time` is
  ///
  ///     (cj M + K) y_{n+1} = f(t_{n+1}) - M known.
  ///
  /// `derivative(y)` evaluates it in correction form. Direct
  /// higher-derivative equations retain the classical level-based formula and
  /// its weights.
  template <class State>
  class StepRequest {
  public:
    struct ClassicalDerivative {
      const StepFormula& formula;
    };

    struct BdfDerivative {
      const State& predictor;
      const State& derivative_predictor;
      const State& history;
      double leading;
    };

    using Derivative = std::variant<ClassicalDerivative, BdfDerivative>;

    StepRequest(double new_time,
      double step_size,
      int step_order,
      Derivative derivative_formula,
      const History<State>& history,
      const State& absolute_scale,
      std::size_t accepted_steps,
      bool initialization = false)
      : time(new_time)
      , step(step_size)
      , order(step_order)
      , step_index(accepted_steps)
      , initializing(initialization)
      , derivative_(std::move(derivative_formula))
      , history_(history)
      , absolute_scale_(absolute_scale)
    {
    }

    const double time;
    const double step;
    const int order;
    const std::size_t step_index;
    const bool initializing;

    const State& absolute_tolerance_reference() const
    {
      return absolute_scale_;
    }

    double derivative_weight(int derivative = 1) const
    {
      if (derivative == 0)
        return 1.0;
      if (const auto* bdf = std::get_if<BdfDerivative>(&derivative_)) {
        require_first_derivative(derivative);
        return bdf->leading;
      }
      return std::get<ClassicalDerivative>(derivative_)
        .formula.weights(derivative)[0];
    }

    State derivative_history(int derivative = 1) const
    {
      if (derivative == 0) {
        State result = previous();
        atstep::axpy(result, -1.0, previous());
        return result;
      }
      if (const auto* bdf = std::get_if<BdfDerivative>(&derivative_)) {
        require_first_derivative(derivative);
        return bdf->history;
      }
      const auto weights = std::get<ClassicalDerivative>(derivative_)
                             .formula.weights(derivative);
      State result = previous();
      atstep::axpy(result, weights[1] - 1.0, previous());
      for (std::size_t j = 2; j < weights.size(); ++j)
        atstep::axpy(result, weights[j], previous(j - 1));
      return result;
    }

    State derivative(const State& y, int derivative = 1) const
    {
      if (derivative == 0)
        return y;
      if (const auto* bdf = std::get_if<BdfDerivative>(&derivative_)) {
        require_first_derivative(derivative);
        State correction = y;
        atstep::axpy(correction, -1.0, bdf->predictor);
        State result = bdf->derivative_predictor;
        atstep::axpy(result, bdf->leading, correction);
        return result;
      }
      State result = derivative_history(derivative);
      atstep::axpy(result, derivative_weight(derivative), y);
      return result;
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
      return history_.state(level);
    }

  private:
    static void require_first_derivative(int derivative)
    {
      if (derivative != 1)
        throw std::invalid_argument("atstep: this BDF request supports "
                                    "only the first derivative");
    }

    Derivative derivative_;
    const History<State>& history_;
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

  template <class State, class Norm>
    requires StateVector<State>
  class TimeStepper;

  /// One attempted step: the time and order it was taken at, whether its
  /// solve converged, and the error estimate of it.
  ///
  /// It is the receipt of an attempt: the stepper creates it, it cannot be
  /// copied, and one call of `accept`, `reject` or `reject_convergence`
  /// settles it. A handle of an older attempt, or of an attempt another
  /// stepper took, is refused rather than committed. The candidate state it
  /// exposes is the stepper's own buffer, valid until the next attempt.
  template <class State>
  class TrialStep {
  public:
    TrialStep(const TrialStep&) = delete;
    TrialStep& operator=(const TrialStep&) = delete;
    TrialStep(TrialStep&&) noexcept = default;
    TrialStep& operator=(TrialStep&&) noexcept = default;

    /// Whether the solve converged.
    bool solved() const noexcept
    {
      return solved_;
    }

    /// Whether the step met the error test. A step that did not is rejected.
    bool passed() const noexcept
    {
      return passed_;
    }

    /// The time the step ends at.
    double time() const noexcept
    {
      return time_;
    }

    /// The size of the step.
    double step() const noexcept
    {
      return step_;
    }

    /// The order the step was taken at.
    int order() const noexcept
    {
      return order_;
    }

    /// The nonlinear iterations the solve reported.
    int iterations() const noexcept
    {
      return iterations_;
    }

    /// The error estimates of the step.
    const StepErrors& errors() const noexcept
    {
      return errors_;
    }

    /// The state the solve left: the candidate solution of this step.
    const State& state() const noexcept
    {
      return *state_;
    }

    /// Label the step as ending at `time` rather than at the time its size
    /// proposes. A landing mode snaps a step onto an output time this way,
    /// within a rounding error; any other use breaks the time the history
    /// records.
    void snap_time(double time) noexcept
    {
      time_ = time;
    }

  private:
    template <class S, class N>
      requires StateVector<S>
    friend class TimeStepper;

    TrialStep(const State& state, double time, double step, int order)
      : state_(&state)
      , time_(time)
      , step_(step)
      , order_(order)
    {
    }

    const State* state_;
    double time_;
    double step_;
    int order_;
    int iterations_ = 0;
    bool solved_ = false;
    bool passed_ = false;
    StepErrors errors_ {};
    /// The attempt of the stepper this handle describes; 0 is not an attempt.
    unsigned long long id_ = 0;
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
  /// relative to the output times.
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
  /// `run` is the same loop for a whole schedule of output times; `attempt`
  /// with `accept`, `reject` and `reject_convergence` is one step at a time,
  /// with the decision taken by the caller.
  ///
  /// `y` comes in holding the predictor of the step — the extrapolation of
  /// the history, which is a good initial guess — and must leave holding the
  /// solution of the step.
  template <class State, class Norm = WeightedRmsNorm>
    requires StateVector<State>
  class TimeStepper {
  public:
    using Solve = std::function<SolveOutcome(const StepRequest<State>&, State&)>;

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
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>) {
        norm_.relative = settings_.relative_tolerance;
        norm_.absolute = settings_.relative_tolerance * settings_.absolute_factor;
        norm_.scaled = settings_.scale_absolute_tolerance;
        norm_.scale = settings_.absolute_scale;
      }
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
      controller_ = StepController(settings_, span);
      started_ = false;
      pending_ = false;
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        norm_.step_index = 0;
      has_previous_ee_ = false;
      previous_order_ = 0;
      last_order_ = 1;
      derivative_ = initial;
      atstep::axpy(derivative_, -1.0, initial);
      probe_ = 1.0;

      history_.reset(time, initial);
      scale_reference_ = initial;

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

      const StepFormula probe_formula(1,
        std::span<const double>(&probe, 1),
        settings_.derivative_orders);
      // COMSOL's backward-Euler initialization advances twice. Its stored
      // initial state is the second result, labeled with the start time.
      work_ = initial;
      for (int j = 1; j <= 2; ++j) {
        StepRequest<State> request {
          time + static_cast<double>(j) * probe,
          probe,
          1,
          typename StepRequest<State>::ClassicalDerivative {
            probe_formula},
          history_,
          history_.state(0),
          0,
          true};
        const SolveOutcome outcome = solve_(request, work_);
        if (not outcome.converged)
          throw std::runtime_error(
            "atstep: the initialization step of the run at t = " + std::to_string(time) + " did not converge");
        scratch_ = work_;
        atstep::axpy(scratch_, -1.0, history_.state(0));
        history_.reset(request.time, work_);
      }

      // The time derivative the first step is read from and predicted by:
      // the increment of the initialization step over its size.
      const double derivative_norm = norm_(scratch_, work_) / probe;

      // The initialization step is not a level of the scheme: its
      // interval is a thousandth of the first step, and a polynomial
      // through it would be meaningless. It leaves the state it produced
      // and the derivative it measured, which is what the first step is
      // predicted — and sized — from. IDA starts the same way, from a
      // state and a time derivative and no history at all.
      derivative_ = scratch_;
      probe_ = probe;
      history_.reset(time + 2.0 * probe, work_);
      scale_reference_ = work_;
      controller_.begin(controller_.first_step(derivative_norm));
      if (fixed_leading()) {
        last_derivative_ = derivative_;
        atstep::axpy(last_derivative_, 1.0 / probe - 1.0, derivative_);
        // Zero only the columns the configured order can reach: the scheme
        // reads `phi_[0..max_order+1]`, and a state allocated here is a
        // state held for the whole run.
        for (int j = 0; j <= settings_.max_order + 1; ++j) {
          State& column = phi_[static_cast<std::size_t>(j)];
          column = work_;
          atstep::axpy(column, -1.0, work_);
        }
        phi_[0] = work_;
        atstep::axpy(
          phi_[1], controller_.step_size() / probe, derivative_);
        psi_.fill(0.0);
        psi_[0] = controller_.step_size();
      }
      started_ = true;
    }

    const Settings& settings() const noexcept
    {
      return settings_;
    }
    const StepController& controller() const noexcept
    {
      return controller_;
    }
    const History<State>& history() const noexcept
    {
      return history_;
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
      return history_.time(0);
    }

    /// The current solution: the newest accepted level's state. It is the
    /// stepper's own buffer, overwritten by the next attempt.
    const State& state() const
    {
      require_started();
      return history_.state(0);
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

    /// Take one attempt at `step`: solve the system the scheme defines, and
    /// estimate the error of the step. Nothing is committed: `accept`
    /// pushes the state onto the history, `reject` and `reject_convergence`
    /// drop it.
    ///
    /// Only the attempt this returns is pending, and only until it is
    /// settled: another attempt, or a second settlement, invalidates it.
    TrialStep<State> attempt(double step)
    {
      require_started();
      pending_ = false; // an attempt that throws leaves none pending
      const double next_time = history_.time(0) + step;
      if (not std::isfinite(step) or not(step > 0.0)
        or not std::isfinite(next_time) or not(next_time > history_.time(0)))
        throw std::invalid_argument("atstep: step must be finite, positive and advance time");

      // Both formulations use the same history bound on startup order: an
      // order-k step reads k levels, and the run has one level per accepted
      // step on top of the one it started from.
      const int levels = std::min(controller_.steps() + 1, max_levels);
      const int order = std::min(controller_.order(), std::max(1, levels - 1));
      controller_.begin_attempt(step, order);

      auto derivative = prepare_derivative(step, order);
      work_ = predictor_;
      TrialStep<State> result {work_, next_time, step, order};
      result.id_ = ++attempt_id_;
      pending_ = true;

      StepRequest<State> request {
        next_time,
        step,
        order,
        std::move(derivative),
        history_,
        scale_reference_,
        static_cast<std::size_t>(controller_.steps())};
      const SolveOutcome outcome = solve_(request, work_);
      result.iterations_ = outcome.iterations;
      if (not outcome.converged)
        return result;
      result.solved_ = true;

      if (fixed_leading() or controller_.error_controlled()) {
        ee_ = work_;
        atstep::axpy(ee_, -1.0, predictor_);
      }
      if (not controller_.error_controlled()) {
        result.passed_ = true;
        return result;
      }
      result.errors_ = estimate(order);
      result.passed_ = result.errors_.passed();
      return result;
    }

    /// Commit the step `attempt` produced.
    ///
    /// @param[in] attempt the pending attempt, whose solve converged and
    /// whose error test passed. The call settles it.
    void accept(TrialStep<State> attempt)
    {
      require_pending(attempt);
      if (not attempt.solved_)
        throw std::logic_error("atstep: a step whose solve did not converge "
                               "cannot be accepted (see reject_convergence)");
      if (not attempt.passed_)
        throw std::logic_error("atstep: the step did not meet the error test ("
          + std::to_string(attempt.errors_.test) + "); reject it");
      if (not std::isfinite(attempt.time_) or not(attempt.time_ > history_.time(0)))
        throw std::invalid_argument(
          "atstep: an accepted step must advance the time of the run");
      pending_ = false;

      if (fixed_leading()) {
        last_derivative_ = derivative_predictor_;
        atstep::axpy(last_derivative_, cj_, ee_);
      }
      history_.push(attempt.time_, work_);
      history_.keep(history_depth());
      last_order_ = attempt.order_;
      if (fixed_leading()) {
        const auto k = static_cast<std::size_t>(attempt.order_);
        commit_table(k);
        phi_[k + 1] = ee_;
        atstep::axpy(phi_[k], 1.0, ee_);
        for (std::size_t j = k; j > 0; --j)
          atstep::axpy(phi_[j - 1], 1.0, phi_[j]);
        psi_ = trial_psi_;
      }
      controller_.accept(attempt.errors_, attempt.step_);
      scale_reference_ = work_;
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        norm_.step_index = static_cast<std::size_t>(controller_.steps());
      if (not fixed_leading() and controller_.error_controlled()) {
        previous_ee_ = ee_;
        previous_order_ = attempt.order_;
        has_previous_ee_ = true;
      }
    }

    /// Drop the step `attempt` produced: it missed the tolerance.
    void reject(TrialStep<State> attempt)
    {
      require_pending(attempt);
      if (not attempt.solved_)
        throw std::logic_error("atstep: the solve of the step did not converge "
                               "(see reject_convergence)");
      pending_ = false;
      controller_.reject(attempt.errors_, attempt.step_);
    }

    /// Drop the step `attempt` produced: its solve did not converge.
    void reject_convergence(TrialStep<State> attempt)
    {
      require_pending(attempt);
      if (attempt.solved_)
        throw std::logic_error(
          "atstep: the solve of the step converged; reject it for its error "
          "test instead");
      pending_ = false;
      controller_.reject_convergence(attempt.step_);
    }

    /// Take one step the error test accepts, and report it.
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
      for (;;) {
        const double from = history_.time(0);
        const double step = controller_.step_size();
        if (not(step > 0.0))
          throw std::runtime_error(
            "atstep: no step is left at t = " + std::to_string(from));
        TrialStep<State> piece = attempt(step);
        if (not piece.solved_) {
          const double tried = piece.step_;
          if (tried <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from) + " did not converge at the smallest step "
              + std::to_string(tried));
          reject_convergence(std::move(piece));
          continue;
        }
        if (not piece.passed_) {
          const double error = piece.errors_.test;
          reject(std::move(piece));
          if (controller_.step_size() <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller_.step_size()) + " and the local error is still " + std::to_string(error) + "; the tolerance cannot be met on this interval");
          continue;
        }
        const StepRecord record = record_of(piece);
        accept(std::move(piece));
        return record;
      }
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
      if (output_times.size() < 2)
        throw std::invalid_argument(
          "atstep: a run needs at least two output times");
      if (output_times.front() != start_time_ or output_times.back() != end_time_)
        throw std::invalid_argument("atstep: the output times must "
                                    "start at the start time and end "
                                    "at the end time of the run");
      for (std::size_t i = 1; i < output_times.size(); ++i)
        if (not(output_times[i] > output_times[i - 1]))
          throw std::invalid_argument(
            "atstep: the output times must increase");

      RunReport report;
      double time = history_.time(0);

      // COMSOL labels the state left by consistent initialization with
      // the requested start time; the initialization solves are not
      // accepted steps.
      std::size_t next = 0;
      if (output_times[0] == start_time_) {
        store(start_time_, history_.state(0), false, hooks, report);
        next = 1;
      }

      bool landed = false; // a step has landed inside the current subinterval
      while (time < end_time_) {
        const double target = next < output_times.size() ? output_times[next] : end_time_;
        double step = controller_.step_size();
        const std::size_t intermediate_index = next + (landed ? 1U : 0U);
        const double intermediate_target = intermediate_index < output_times.size()
          ? output_times[intermediate_index]
          : std::numeric_limits<double>::infinity();
        const bool intermediate_landing = settings_.step_mode == StepMode::intermediate and time + step >= intermediate_target;
        switch (settings_.step_mode) {
        case StepMode::strict:
          // COMSOL avoids a tiny remainder: split the last two steps
          // equally if one full proposed step would leave less than
          // one.
          step = landing_step(step, time, target);
          break;
        case StepMode::intermediate:
          // Every subinterval holds at least one step: this one is
          // cut short to land inside the current one.
          if (intermediate_landing)
            step = (intermediate_target - time) * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
          break;
        case StepMode::free:
        case StepMode::manual:
          break;
        }
        if (not settings_.interpolate_end_time) {
          step = landing_step(step, time, end_time_);
        }
        if (not(step > 0.0))
          throw std::runtime_error("atstep: no step is left at t = " + std::to_string(time));

        const double from = time;
        TrialStep<State> piece = attempt(step);
        const double boundary = intermediate_landing
          ? intermediate_target
          : (settings_.step_mode == StepMode::strict ? target
                                                     : end_time_);
        if ((settings_.step_mode == StepMode::strict or intermediate_landing or not settings_.interpolate_end_time) and std::abs(piece.time() - boundary) <= 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(time) + std::abs(boundary)))
          piece.snap_time(boundary);
        if (not piece.solved()) {
          ++report.failed_solves;
          const double tried = piece.step();
          reject_convergence(std::move(piece));
          if (tried <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from) + " did not converge at the "
                                                                                 "smallest step "
              + std::to_string(tried));
          continue;
        }
        if (not piece.passed()) {
          ++report.rejected;
          const double error = piece.errors().test;
          reject(std::move(piece));
          if (controller_.step_size() <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller_.step_size()) + " and the local error is still " + std::to_string(error) + "; the tolerance cannot be met on this interval");
          continue;
        }

        const StepRecord record = record_of(piece);
        accept(std::move(piece));
        time = record.time;
        ++report.steps;
        report.max_order = std::max(report.max_order, record.order);
        if (report.steps == 1)
          report.first_step = record.step;
        report.last_step = record.step;
        report.final_time = time;
        if (hooks.on_step)
          hooks.on_step(record, history_.state(0));

        if (settings_.store_mode == StoreMode::steps) {
          if (time > end_time_ and settings_.interpolate_end_time)
            store(end_time_,
              interpolate(end_time_),
              true,
              hooks,
              report);
          else
            store(time, history_.state(0), false, hooks, report);
        }
        next = store_outputs(next, time, output_times, hooks, report);
        // A step crossing an output time can already lie inside the
        // next interval. Count that landing instead of forcing another
        // step.
        landed = next > 0 and time > output_times[next - 1];
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
        if (derivative == 1 and time == history_.time(0))
          return last_derivative_;
        const double delta = time - history_.time(0);
        double coefficient = 1.0;
        double slope = 0.0;
        double gamma = delta / psi_[0];
        interpolated_ = phi_[0];
        if (derivative == 1)
          atstep::axpy(interpolated_, -1.0, phi_[0]);
        for (int j = 1; j <= last_order_; ++j) {
          slope = slope * gamma + coefficient / psi_[static_cast<std::size_t>(j - 1)];
          coefficient *= gamma;
          atstep::axpy(interpolated_,
            derivative == 0 ? coefficient : slope,
            phi_[static_cast<std::size_t>(j)]);
          gamma = (delta + psi_[static_cast<std::size_t>(j - 1)]) / psi_[static_cast<std::size_t>(j)];
        }
        return interpolated_;
      }
      const std::size_t count = std::min<std::size_t>(
        static_cast<std::size_t>(last_order_) + 1, history_.size());
      std::array<double, max_levels> nodes {};
      for (std::size_t j = 0; j < count; ++j)
        nodes[j] = history_.time(j) - time;
      const auto weights = lagrange_derivative(
        std::span<const double>(nodes.data(), count),
        0.0,
        derivative);
      combine(interpolated_,
        std::span<const double>(weights.data(), count));
      return interpolated_;
    }

  private:
    /// Refuse a handle that is not the pending attempt of this stepper.
    void require_pending(const TrialStep<State>& attempt) const
    {
      if (not pending_ or attempt.id_ != attempt_id_ or attempt.state_ != &work_)
        throw std::logic_error(
          "atstep: the step is not the pending attempt of this stepper");
    }

    static StepRecord record_of(const TrialStep<State>& attempt)
    {
      return StepRecord {attempt.time_,
        attempt.step_,
        attempt.order_,
        attempt.iterations_,
        attempt.errors_.test};
    }

    typename StepRequest<State>::Derivative prepare_derivative(double step,
      int order)
    {
      if (fixed_leading()) {
        prepare_fixed_leading(step, order);
        return typename StepRequest<State>::BdfDerivative {
          predictor_,
          derivative_predictor_,
          derivative_history_,
          cj_};
      }
      const std::size_t count = history_.size();
      std::array<double, max_levels> steps {};
      steps[0] = step;
      for (std::size_t j = 1; j < count; ++j)
        steps[j] = history_.time(j - 1) - history_.time(j);
      formula_ = StepFormula(order,
        std::span<const double>(steps.data(), count),
        settings_.derivative_orders);
      predictor(step);
      return typename StepRequest<State>::ClassicalDerivative {formula_};
    }

    double error_norm(const State& error)
    {
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        return norm_(error, history_.state(0), scale_reference_);
      else
        return norm_(error, history_.state(0));
    }

    static double landing_step(double step, double time, double target)
    {
      const double remaining = target - time;
      const double fuzz = 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(time) + std::abs(target));
      if (time + step > target) {
        return remaining * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
      }
      if (std::abs(time + step - target) <= fuzz)
        return step;
      return 2.0 * step > remaining ? 0.5 * remaining : step;
    }

    bool fixed_leading() const
    {
      return settings_.derivative_orders == 1;
    }

    /// The levels the history keeps: the fixed-leading table stands in for
    /// every level but the newest, so a first-order equation keeps one.
    std::size_t history_depth() const
    {
      return fixed_leading() ? 1U : static_cast<std::size_t>(max_levels);
    }

    void prepare_fixed_leading(double step, int order)
    {
      trial_psi_ = psi_;
      if (controller_.steps() == 0) {
        trial_psi_[0] = step;
        phi_[1] = derivative_;
        atstep::axpy(phi_[1], step / probe_ - 1.0, derivative_);
      }
      // The trial reads the accepted table and writes no column of it: the
      // rescaled column of the order it is at is `beta_j phi_[j]`, which is
      // a scalar multiple of the accepted column, so the factors `beta_j`
      // are all it keeps. A dropped attempt therefore has nothing to put
      // back, and an accepted one writes the factors in.
      predictor_ = phi_[0];
      derivative_predictor_ = phi_[0];
      atstep::axpy(derivative_predictor_, -1.0, phi_[0]);
      double beta = 1.0;
      double distance = step;
      double alpha = 1.0;
      double gamma = 0.0;
      double alpha0 = -1.0;
      double harmonic = 0.0;
      sigma_[0] = 1.0;
      for (int i = 1; i <= order; ++i) {
        const auto j = static_cast<std::size_t>(i);
        const double old_distance = trial_psi_[j - 1];
        trial_psi_[j - 1] = distance;
        beta *= distance / old_distance;
        distance = old_distance + step;
        gamma += alpha / step;
        alpha = step / distance;
        if (i < order)
          alpha0 -= alpha;
        harmonic += 1.0 / static_cast<double>(i);
        sigma_[j] = static_cast<double>(i) * sigma_[j - 1] * alpha;
        beta_[j] = beta;
        atstep::axpy(predictor_, beta, phi_[j]);
        atstep::axpy(derivative_predictor_, gamma * beta, phi_[j]);
      }
      trial_psi_[static_cast<std::size_t>(order)] = distance;
      ck_ = std::max(alpha, std::abs(alpha - harmonic - alpha0));
      cj_ = harmonic / step;
      derivative_history_ = derivative_predictor_;
      atstep::axpy(derivative_history_, -cj_, predictor_);
    }

    /// Write the factors the trial read into the table. This is the whole of
    /// what an attempt owes it: the columns the predictor and the estimates
    /// read, in the form they read them in. The accepted column is rescaled
    /// where it stands, so the table holds one column per order throughout.
    void commit_table(std::size_t order)
    {
      for (std::size_t j = 1; j <= order; ++j)
        atstep::axpy(phi_[j], beta_[j] - 1.0, phi_[j]);
    }

    void require_started() const
    {
      if (not started_)
        throw std::logic_error(
          "atstep: the run has not been started (see start)");
    }

    /// `out = weights[0] y_0 + ... + weights[count - 1] y_{count - 1}`,
    /// from the levels of the history.
    void combine(State& out, std::span<const double> weights) const
    {
      const std::size_t count = std::min(weights.size(), history_.size());
      out = history_.state(0);
      atstep::axpy(out, weights[0] - 1.0, history_.state(0));
      for (std::size_t j = 1; j < count; ++j)
        atstep::axpy(out, weights[j], history_.state(j));
    }

    /// The extrapolation of the history to the new time, in `predictor_`. A
    /// history of a single level is the initial state, and the step from it
    /// is predicted by the derivative the initialization step measured.
    void predictor(double step)
    {
      if (history_.size() == 1) {
        // `derivative_` is the increment of the initialization step, so
        // the extrapolation of the step being taken is `step / probe`
        // of it.
        predictor_ = history_.state(0);
        atstep::axpy(predictor_, step / probe_, derivative_);
        return;
      }
      combine(predictor_, formula_.predictor());
    }

    /// The error estimates of the step at `order`, whose state is `work_`
    /// and whose predictor is `predictor_`.
    StepErrors estimate(int order)
    {
      StepErrors errors;
      // IDA fixes the weights from the last accepted state before
      // attempting a step; every estimate and retry uses those same
      // weights.
      const double enorm = error_norm(ee_);
      if (fixed_leading()) {
        const auto k = static_cast<std::size_t>(order);
        errors.test = ck_ * enorm;
        errors.at = sigma_[k] * enorm;
        scratch_ = ee_;
        if (order >= 2) {
          atstep::axpy(scratch_, beta_[k], phi_[k]);
          errors.below = sigma_[k - 1] * error_norm(scratch_);
        }
        if (order >= 3) {
          atstep::axpy(scratch_, beta_[k - 1], phi_[k - 1]);
          errors.below2 = sigma_[k - 2] * error_norm(scratch_);
        }
        scratch_ = ee_;
        atstep::axpy(scratch_, -1.0, phi_[k + 1]);
        errors.above = error_norm(scratch_) / static_cast<double>(order + 2);
        return errors;
      }
      errors.test = formula_.ck() * enorm;
      errors.at = formula_.sigma(order) * enorm;

      // The orders below read the polynomial through one level fewer
      // each, and the estimate above compares this step's correction with
      // the one the previous step left.
      if (order >= 2)
        errors.below = formula_.sigma(order - 1) * correction(formula_.extrapolation(order));
      if (order >= 3)
        errors.below2 = formula_.sigma(order - 2) * correction(formula_.extrapolation(order - 1));
      if (has_previous_ee_ and previous_order_ == order) {
        scratch_ = ee_;
        atstep::axpy(scratch_, -1.0, previous_ee_);
        errors.above = error_norm(scratch_) / static_cast<double>(order + 2);
      }
      return errors;
    }

    /// The weighted norm of `y - p`, `p` being the polynomial through the
    /// levels with `weights`.
    double correction(std::span<const double> weights)
    {
      combine(scratch_, weights);
      atstep::axpy(scratch_, -1.0, work_);
      return error_norm(scratch_);
    }

    void store(double time,
      const State& state,
      bool interpolated,
      const RunHooks<State>& hooks,
      RunReport& report)
    {
      ++report.stored;
      scale_reference_ = state;
      if (hooks.on_store)
        hooks.on_store(time, state, interpolated);
    }

    std::size_t store_outputs(std::size_t next,
      double to,
      std::span<const double> output_times,
      const RunHooks<State>& hooks,
      RunReport& report)
    {
      while (next < output_times.size() and output_times[next] <= to) {
        const double target = output_times[next];
        switch (settings_.store_mode) {
        case StoreMode::interpolate:
          store(target, interpolate(target), true, hooks, report);
          break;
        case StoreMode::closest: {
          // COMSOL's tstepsclosest stores the first step that reaches
          // an output time, even if the preceding step is nearer to
          // it.
          if (target == end_time_ and to > end_time_ and settings_.interpolate_end_time)
            store(target, interpolate(target), true, hooks, report);
          else
            store(to, history_.state(0), false, hooks, report);
          break;
        }
        case StoreMode::steps:
          break; // the output times are only counted, not stored
        }
        ++next;
      }
      return next;
    }

    Settings settings_;
    Solve solve_;
    Norm norm_;

    StepController controller_ {Settings {}, 1.0};
    History<State> history_;

    double start_time_ = 0.0;
    double end_time_ = 0.0;
    /// The size of the initialization step (see `start`).
    double probe_ = 1.0;
    bool started_ = false;

    /// The attempt a `TrialStep` describes, and whether one is still to be
    /// settled. A handle whose id is not the current one is not that attempt.
    unsigned long long attempt_id_ = 0;
    bool pending_ = false;

    /// The buffers of a step: the initial derivative, the predictor the solve
    /// starts from, the state it produces, the correction of the step, the
    /// previous step's correction, a scratch state, and the interpolated
    /// state.
    State derivative_ {};
    State predictor_ {};
    State work_ {};
    State ee_ {};
    State previous_ee_ {};
    State scratch_ {};
    State interpolated_ {};
    State derivative_history_ {};
    State derivative_predictor_ {};
    State last_derivative_ {};
    State scale_reference_ {};
    /// The columns of the first-order table: `phi_[j]` is the j-th scaled
    /// difference of the accepted steps, and `phi_[0]` the level the predictor
    /// extrapolates from — which the solve's own state need not equal in its
    /// last bits. Only `accept` writes a column of it; an attempt reads it.
    std::array<State, max_levels> phi_ {};
    /// The factor of `phi_[j]` at the step size of the pending attempt, which
    /// `accept` commits and the estimates read back. The scaled column the
    /// scheme reads is this factor times the accepted one, so a trial costs
    /// no state at all.
    std::array<double, max_levels> beta_ {};
    std::array<double, max_levels> psi_ {};
    std::array<double, max_levels> trial_psi_ {};
    std::array<double, max_levels> sigma_ {};
    double cj_ = 0.0;
    double ck_ = 0.5;
    bool has_previous_ee_ = false;
    int previous_order_ = 0;
    int last_order_ = 1;
    StepFormula formula_ {};
  };

} // namespace atstep
