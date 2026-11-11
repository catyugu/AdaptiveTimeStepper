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
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace atstep {

  /// Error norm context: estimate, accepted state, last returned state and step count.
  template <class State>
  using ErrorNorm = std::function<double(const State& error,
    const State& accepted,
    const State& returned,
    std::size_t accepted_steps)>;

  /// Coefficients for `y' = cj*y + known` at `time`.
  /// For `M y' + K y = f`, solve `(cj*M + K)*y = f - M*known`.
  /// The solve receives the predictor and replaces it with the solution.
  /// Requests and their state references are valid only during the callback.
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
    /// The accepted steps behind the run.
    const std::size_t accepted_steps;
    const bool initializing;

    /// The weight `y'` carries on the state being solved for.
    double derivative_weight() const noexcept
    {
      return weight_;
    }

    /// The known term of the derivative: `y' = cj*y + known`.
    const State& derivative_history() const noexcept
    {
      return known_;
    }

    /// Write `predicted_derivative + cj*(y - predictor)` to caller-owned `out`.
    /// This correction form reduces cancellation at large `cj`; `out` may alias `y`.
    void derivative(const State& y, State& out) const
    {
      out = y;
      detail::axpy(out, -1.0, predictor_);
      detail::axpy(out, weight_ - 1.0, out);
      detail::axpy(out, 1.0, derivative_predictor_);
    }

    /// The last accepted level, which the error weights' relative part is
    /// measured against.
    const State& accepted_state() const noexcept
    {
      return accepted_;
    }

    /// The value last returned to the caller, which may be interpolated; the
    /// scaled absolute tolerance follows it. A model that repeats the error
    /// weight for its own nonlinear convergence test reads both states here.
    const State& returned_state() const noexcept
    {
      return returned_;
    }

    StepRequest(double new_time,
      double step_size,
      int step_order,
      std::size_t steps,
      bool initialization,
      double weight,
      const State& known,
      const State& predictor,
      const State& derivative_predictor,
      const State& accepted,
      const State& returned)
      : time(new_time)
      , step(step_size)
      , order(step_order)
      , accepted_steps(steps)
      , initializing(initialization)
      , weight_(weight)
      , known_(known)
      , predictor_(predictor)
      , derivative_predictor_(derivative_predictor)
      , accepted_(accepted)
      , returned_(returned)
    {
    }

  private:
    double weight_ = 0.0;
    const State& known_;
    const State& predictor_;
    const State& derivative_predictor_;
    const State& accepted_;
    const State& returned_;
  };

  /// An accepted step, as the log sees it.
  struct StepRecord {
    double time = 0.0;
    double step = 0.0;
    int order = 0;
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
    /// Stored time, state and interpolation flag. State references last only for the call.
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

  /// BDF integration owning history, work buffers, error control and output.
  /// Noncopyable and nonmovable; `start` resets the run. See doc/design.md.
  template <class State>
  class TimeStepper {
    static_assert(detail::StateVector<State>,
      "atstep: State must be default-constructible, copyable and support axpy");

  public:
    /// The solve of one step: whether its nonlinear iteration converged. A step
    /// that did not is retried at a quarter of its size.
    using Solve = std::function<bool(const StepRequest<State>&, State&)>;
    using Norm = ErrorNorm<State>;

    TimeStepper(const TimeStepper&) = delete;
    TimeStepper& operator=(const TimeStepper&) = delete;
    TimeStepper(TimeStepper&&) = delete;
    TimeStepper& operator=(TimeStepper&&) = delete;

    /// Default weighted RMS with the configured variable partition, scales and tolerances.
    static Norm default_norm(const ErrorTest& error)
    {
      return [weights = detail::WeightedRmsNorm::from(error)](
               const State& estimate,
               const State& accepted,
               const State& returned,
               std::size_t accepted_steps) {
        return weights(estimate, accepted, returned, accepted_steps);
      };
    }

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one step (see the class comment).
    /// @param[in] norm the norm the error of a step is measured in. Empty means
    /// `default_norm(settings.error)`, which a state must expose its dofs to
    /// accept; a state of its own arithmetic brings its own norm.
    explicit TimeStepper(const Settings& settings,
      Solve solve,
      Norm norm = {})
      : settings_ {}
      , solve_(std::move(solve))
      , norm_(resolve_norm(std::move(norm), settings.error))
    {
      // Copy active alternatives directly. Copying the nested variants as one
      // object triggers GCC 13/14's maybe-uninitialized diagnostic in optimized
      // consumers, on the inactive std::function storage.
      settings_.error = settings.error;
      settings_.initialization = settings.initialization;
      settings_.initialization_fraction = settings.initialization_fraction;
      settings_.storage = settings.storage;
      settings_.store_every = settings.store_every;
      if (const auto* manual = std::get_if<ManualSteps>(&settings.steps))
        settings_.steps.emplace<ManualSteps>(*manual);
      else {
        const auto& source = std::get<AdaptiveSteps>(settings.steps);
        auto& target = std::get<AdaptiveSteps>(settings_.steps);
        target.placement = source.placement;
        target.interpolate_end_time = source.interpolate_end_time;
        target.min_order = source.min_order;
        target.max_order = source.max_order;
        target.initial_step = source.initial_step;
        target.min_step_fraction = source.min_step_fraction;
        std::visit([&target](const auto& maximum) {
          target.max_step.template emplace<std::decay_t<decltype(maximum)>>(maximum);
        },
          source.max_step);
      }
      validate(settings_);
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
    }

    /// Start a run over [time, end_time]. Default initialization performs two
    /// backward-Euler probes; initialization off uses the supplied state and zero derivative.
    void start(double time, const State& initial, double end_time)
    {
      prepare_run(time, end_time);
      if (settings_.initialization == Initialization::off) {
        State derivative = initial;
        detail::axpy(derivative, -1.0, initial);
        begin_run(
          time, initial, derivative, norm_(derivative, initial, initial, 0), time, end_time);
        return;
      }
      const double span = end_time - time;

      // Probe basis: explicit initial step, or 0.1% of span/manual step.
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
      const double probe = probe_step * probe_fraction * settings_.initialization_fraction;
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
        if (not solve(time + static_cast<double>(j) * probe, probe, 1, cj, level,
              level, true))
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
      begin_run(time, level, derivative, norm_(increment, level, level, 0) / probe,
        initialized_time, end_time);
    }

    /// Start from a state and a derivative the caller gives, as a restart from a
    /// recorded level does. The initialization setting is not read: the derivative
    /// is the caller's, so a probe has none to measure.
    void start(double time, const State& initial, const State& derivative, double end_time)
    {
      prepare_run(time, end_time);
      begin_run(
        time, initial, derivative, norm_(derivative, initial, initial, 0), time, end_time);
    }

    /// The time the run has reached: the newest accepted level's.
    double time() const
    {
      require_started();
      return accepted_.time;
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

    /// Advance one accepted step, retrying failed attempts and respecting end placement.
    /// @throws std::runtime_error if the minimum step cannot converge or meet tolerance.
    StepRecord advance()
    {
      require_started();
      const StepRecord record = take_step(std::nullopt);
      note_step(record);
      return record;
    }

    /// Whether the run has reached its end time.
    bool finished() const
    {
      require_started();
      return accepted_.time >= end_time_;
    }

    /// The counts of the run so far.
    const RunReport& report() const noexcept
    {
      return report_;
    }

    /// The settings the run was built from: the copy the constructor took, which an
    /// adapter layered on this kernel reads instead of keeping its own.
    const Settings& settings() const noexcept
    {
      return settings_;
    }

    /// Integrate from the start time to the end time, reporting the values the
    /// storage asks for. Call after `start`, before any `advance`; `output_times`
    /// increases from the start time to the end time, and the schedule copies it.
    RunReport run(const std::vector<double>& output_times,
      RunHooks<State> hooks = {})
    {
      require_started();
      if (controller_->steps() != 0)
        throw std::logic_error("atstep: run requires a fresh start");
      detail::StepSchedule schedule(
        settings_, output_times, start_time_, end_time_);

      // COMSOL labels the state left by consistent initialization with
      // the requested start time; the initialization solves are not
      // accepted steps.
      store(start_time_, accepted_.state, false, hooks);

      // The schedule only says where a step may land and what it stores. One
      // driver, two callers.
      while (not finished()) {
        const StepRecord record = take_step(schedule);
        note_step(record);
        if (hooks.on_step)
          hooks.on_step(record, accepted_.state);
        schedule.for_each_store(record.time, static_cast<std::size_t>(report_.steps),
          [&](double output_time, bool interpolated) {
            if (not interpolated) {
              store(output_time, accepted_.state, false, hooks);
              return;
            }
            interpolate(output_time, work_.scratch);
            store(output_time, work_.scratch, true, hooks);
          });
      }
      report_.final_time = accepted_.time;
      return report_;
    }

    /// Write dense value (derivative=0) or first derivative (derivative=1) to `out`.
    /// At the newest node use the solve derivative; inside a step use the polynomial
    /// derivative. The caller owns `out` and may reuse it.
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
    /// The tail both starts share: the level, its derivative, the first step and
    /// the difference table they imply.
    void begin_run(double time,
      const State& level,
      const State& derivative,
      double derivative_norm,
      double initialized_time,
      double end_time)
    {
      accepted_.state = level;
      accepted_.derivative = derivative;
      accepted_.time = initialized_time;
      accepted_.order = 1;
      returned_ = level;
      controller_->begin(derivative_norm, end_time - initialized_time, initialized_time);
      // Zero only the columns the configured order can reach: the scheme
      // reads `phi[0..max_order]` and the column above it, and a state
      // allocated here is a state held for the whole run.
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

    /// The preamble both starts share: the bounds, the counts, the controller and
    /// the placement rule.
    void prepare_run(double time, double end_time)
    {
      require_bounds(time, end_time);
      started_ = false;
      report_ = {};
      controller_.emplace(settings_, end_time - time);
      const detail::PlacementRule rule = detail::placement_rule(settings_);
      interpolate_end_time_ = rule.interpolate_end_time;
    }

    /// The bounds of a run.
    void require_bounds(double time, double end_time) const
    {
      if (not std::isfinite(time) or not std::isfinite(end_time) or not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
    }

    /// One attempt, created and consumed inside one pass of `take_step`: a refused
    /// attempt is a `continue` of that loop, never a state the stepper keeps.
    struct Trial {
      double step = 0.0;
      double time = 0.0;
      int order = 1;
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

    /// The accepted level and the scaled differences. A trial only reads them;
    /// `phi[0]` is accumulated and can differ from the solved state in its last
    /// bits, which the weights must keep.
    struct History {
      State state {};
      State derivative {};
      double time = 0.0;
      int order = 1;
      std::array<State, max_levels> phi {};
      std::array<double, max_levels> psi {};
    };

    /// Reusable predictor, predicted derivative, solution, correction and scratch.
    /// Scratch holds `known` during the solve, then serves estimates and output.
    struct Work {
      State predictor {};
      State derivative_predictor {};
      State solution {};
      State correction {};
      State scratch {};
    };

    /// Construct the request only for the duration of the callback. This is
    /// also the buffer-lifetime boundary between solving and error estimation.
    bool solve(double time, double step, int order, double cj,
      const State& previous, const State& reference, bool initializing = false)
    {
      work_.scratch = work_.derivative_predictor;
      detail::axpy(work_.scratch, -cj, work_.predictor);
      work_.solution = work_.predictor;
      const StepRequest<State> request {time, step, order,
        initializing ? 0 : static_cast<std::size_t>(controller_->steps()),
        initializing, cj, work_.scratch, work_.predictor,
        work_.derivative_predictor, previous, reference};
      return solve_(request, work_.solution);
    }

    /// Take one step the error test accepts. Accepted history is written at the
    /// end of a pass and nowhere else, so a refusal leaves the run as the accepted
    /// steps left it; a schedule only says where the step may land.
    StepRecord take_step(std::optional<std::reference_wrapper<detail::StepSchedule>> schedule)
    {
      const double from = accepted_.time;
      auto& controller = *controller_;
      for (;;) {
        double step = controller.step_size();
        std::optional<double> boundary;
        if (schedule) {
          const detail::StepSchedule::Landing placement
            = schedule->get().place(step, from);
          step = placement.step;
          boundary = placement.boundary;
        }
        else if (not interpolate_end_time_) {
          // A run of the caller's own lands on the end time exactly as a
          // scheduled one does when the end is not interpolated.
          step = detail::landing_step(step, from, end_time_);
          if (detail::lands_on(from, step, end_time_))
            boundary = end_time_;
        }
        // A step that would stop within the clock's own roundings of the end time
        // ends on it, so the end is a time the run reaches and not one it misses by
        // a few ulps.
        if (not boundary)
          step = detail::end_landing_step(step, from, end_time_);
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

        if (not solve(next_time, step, trial.order, trial.cj, accepted_.state,
              returned_)) {
          ++report_.failed_solves;
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
            ++report_.rejected;
            controller.reject(trial.errors, step);
            if (controller.step_size() <= controller.smallest_step())
              throw std::runtime_error(
                "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller.step_size()) + " and the local error is still " + std::to_string(trial.errors.test) + "; the tolerance cannot be met on this interval");
            continue;
          }
        }

        if (boundary)
          trial.time = *boundary;
        const StepRecord record {trial.time, trial.step, trial.order, trial.errors.test};
        accept(trial);
        return record;
      }
    }

    /// Update history in place, reusing its storage. Ordinary rejected attempts
    /// never enter here and leave accepted history intact.
    void accept(const Trial& trial)
    {
      // The derivative the level pairs with: the predictor's derivative plus the
      // correction's, which is the formula's derivative at the new level.
      work_.scratch = work_.derivative_predictor;
      detail::axpy(work_.scratch, trial.cj, work_.correction);
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
      returned_ = work_.solution;
    }

    /// Record an accepted step in the run's counts.
    void note_step(const StepRecord& record)
    {
      ++report_.steps;
      report_.max_order = std::max(report_.max_order, record.order);
      report_.final_time = record.time;
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

    /// The error estimates of the trial: the accepted step's, and the ones the
    /// order selection compares. IDA fixes the weights from the last accepted
    /// state before attempting a step, so every estimate and retry uses them.
    detail::StepErrors estimate(const Trial& trial)
    {
      detail::StepErrors errors;
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

    /// The norm a run uses: the caller's, or COMSOL's when the state can be read
    /// dof by dof.
    static Norm resolve_norm(Norm norm, const ErrorTest& error)
    {
      if constexpr (detail::Indexable<State>) {
        return norm ? std::move(norm) : default_norm(error);
      }
      else {
        if (not norm)
          throw std::invalid_argument(
            "atstep: this state has no size() and operator[], so the run needs a norm");
        return std::move(norm);
      }
    }

    /// The norm of the run's error, with the context COMSOL's weight reads.
    double error_norm(const State& error)
    {
      return norm_(error, accepted_.state, returned_,
        controller_ ? static_cast<std::size_t>(controller_->steps()) : 0);
    }

    void store(double time, const State& state, bool interpolated, const RunHooks<State>& hooks)
    {
      ++report_.stored;
      returned_ = state;
      if (hooks.on_store)
        hooks.on_store(time, state, interpolated);
    }

    Settings settings_;
    Solve solve_;
    Norm norm_;

    std::optional<detail::StepController> controller_;
    History accepted_;
    Work work_;
    /// The value last returned to the caller: the scaled absolute tolerance
    /// follows it, and it may be interpolated. The dense output itself is the
    /// caller's, so the stepper keeps no buffer of it.
    State returned_ {};

    double start_time_ = 0.0;
    double end_time_ = 0.0;
    /// Whether a step may pass the end time, whose value is then interpolated.
    /// The schedule resolves the same rule for its own use.
    bool interpolate_end_time_ = true;
    bool started_ = false;
    RunReport report_ {};
  };

} // namespace atstep
