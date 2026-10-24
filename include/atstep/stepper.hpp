// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/controller.hpp"
#include "detail/output.hpp"
#include "detail/state.hpp"
#include "norm.hpp"
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

namespace atstep {

  template <class State, class Norm>
    requires detail::StateVector<State>
  class TimeStepper;

  /// Borrowed coefficients for one first-order solve. The step holds the
  /// derivative to
  ///
  ///     y' = cj*y + known,   cj = derivative_weight(),
  ///                          known = derivative_history(),
  ///
  /// so `M y' + K y = f` at `time` is `(cj M + K) y = f - M known`.
  /// `derivative(y, out)` evaluates `y'` in correction form. The request and
  /// its references live only for the solve callback.
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
      return absolute_reference_;
    }

    double derivative_weight() const noexcept
    {
      return leading_;
    }

    const State& derivative_history() const noexcept
    {
      return derivative_history_;
    }

    /// `out = y'`. Input and output may alias; `out` must not be a state
    /// borrowed from this request.
    void derivative(const State& y, State& out) const
    {
      derivative_scratch_ = y;
      detail::axpy(derivative_scratch_, -1.0, predictor_);
      out = derivative_predictor_;
      detail::axpy(out, leading_, derivative_scratch_);
    }

    /// The accepted state, the relative error-weight reference.
    const State& previous() const noexcept
    {
      return previous_;
    }

  private:
    template <class S, class N>
      requires detail::StateVector<S>
    friend class TimeStepper;

    StepRequest(double new_time,
      double step_size,
      int step_order,
      const State& predictor,
      const State& derivative_predictor,
      double leading,
      const State& previous,
      const State& derivative_history,
      State& derivative_scratch,
      const State& absolute_reference,
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
      , absolute_reference_(absolute_reference)
    {
    }

    const State& predictor_;
    const State& derivative_predictor_;
    double leading_ = 0.0;
    const State& previous_;
    const State& derivative_history_;
    State& derivative_scratch_;
    const State& absolute_reference_;
  };

  struct SolveOutcome {
    bool converged = false;
    int iterations = 0;
  };

  struct StepRecord {
    double time = 0.0;
    double step = 0.0;
    int order = 0;
    int iterations = 0;
    /// The error estimate of the step, at most one when it passed.
    double error = 0.0;
  };

  template <class State>
  struct RunHooks {
    std::function<void(const StepRecord&, const State&)> on_step;
    std::function<void(double, const State&, bool)> on_store;
  };

  struct RunReport {
    int steps = 0;
    int rejected = 0;
    int failed_solves = 0;
    int stored = 0;
    int max_order = 1;
    double final_time = 0.0;
  };

  /// Adaptive BDF integration of a first-order state. The caller supplies the
  /// algebraic solve; the stepper owns history, work buffers, error control and
  /// output. `advance` retries internally until one step is accepted; `run`
  /// performs the same loop over a complete output schedule.
  ///
  /// The integrator is neither copyable nor movable; `start` resets it. State
  /// operations must not throw: history is published in place.
  template <class State, class Norm = WeightedRmsNorm>
    requires detail::StateVector<State>
  class TimeStepper {
  public:
    using Solve = std::function<SolveOutcome(const StepRequest<State>&, State&)>;

    TimeStepper(const TimeStepper&) = delete;
    TimeStepper& operator=(const TimeStepper&) = delete;
    TimeStepper(TimeStepper&&) = delete;
    TimeStepper& operator=(TimeStepper&&) = delete;

    /// The default norm, configured from `settings`.
    TimeStepper(const Settings& settings, Solve solve)
      : settings_(settings)
      , solve_(std::move(solve))
      , norm_(Norm {})
    {
      validate(settings_);
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
      configure_default_norm();
    }

    /// A norm of the caller's, used as given.
    TimeStepper(const Settings& settings, Solve solve, Norm norm)
      : settings_(settings)
      , solve_(std::move(solve))
      , norm_(std::move(norm))
    {
      validate(settings_);
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
    }

    /// Consistent initialization by two artificial backward-Euler steps. The
    /// resulting state is labelled with `time`, and its derivative is what the
    /// first step is read from.
    void start(double time, const State& initial, double end_time)
    {
      if (not std::isfinite(time) or not std::isfinite(end_time)
        or not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
      started_ = false;
      const double span = end_time - time;
      controller_.emplace(settings_, span, end_time);
      update_norm_step_index();

      const double probe = controller_->initial_step_base()
        * settings_.initial_step_fraction * settings_.backward_euler_fraction;
      const double cj = 1.0 / probe;
      const double first_time = time + probe;
      const double initialized_time = time + 2.0 * probe;
      if (not std::isfinite(probe) or not(probe > 0.0)
        or not std::isfinite(cj) or not std::isfinite(initialized_time)
        or not(first_time > time) or not(initialized_time > first_time)
        or not(initialized_time < end_time))
        throw std::invalid_argument(
          "atstep: initialization steps must be finite and advance within the run");

      State level = initial;
      State increment {};
      for (int j = 1; j <= 2; ++j) {
        work_.predictor = level;
        work_.derivative_predictor = work_.predictor;
        detail::axpy(work_.derivative_predictor, -1.0, work_.predictor);
        const SolveOutcome outcome = solve(time + static_cast<double>(j) * probe,
          probe, 1, cj, level, level, true);
        if (not outcome.converged)
          throw std::runtime_error("atstep: the initialization step at t = "
            + std::to_string(time) + " did not converge");
        increment = work_.solution;
        detail::axpy(increment, -1.0, level);
        level = work_.solution;
      }

      const double derivative_norm = norm_(increment, level) / probe;
      State derivative = increment;
      detail::axpy(derivative, cj - 1.0, increment);

      accepted_.state = level;
      accepted_.derivative = derivative;
      accepted_.time = initialized_time;
      accepted_.order = 1;
      output_.reference = &accepted_.state;
      controller_->begin(controller_->first_step(derivative_norm));
      reset_table(level, controller_->step_size(), derivative);
      start_time_ = time;
      end_time_ = end_time;
      started_ = true;
    }

    /// Consistent initial values supplied by the caller: no artificial step,
    /// no perturbation of the state. `derivative` is the time derivative of the
    /// state at `time`.
    void start(double time,
      const State& initial,
      const State& derivative,
      double end_time)
    {
      if (not std::isfinite(time) or not std::isfinite(end_time)
        or not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
      started_ = false;
      controller_.emplace(settings_, end_time - time, end_time);
      update_norm_step_index();

      const double derivative_norm = norm_(derivative, initial);
      accepted_.state = initial;
      accepted_.derivative = derivative;
      accepted_.time = time;
      accepted_.order = 1;
      output_.reference = &accepted_.state;
      controller_->begin(controller_->first_step(derivative_norm));
      reset_table(initial, controller_->step_size(), derivative);
      start_time_ = time;
      end_time_ = end_time;
      started_ = true;
    }

    double time() const
    {
      require_started();
      return accepted_.time;
    }

    const State& state() const
    {
      require_started();
      return accepted_.state;
    }

    double step_size() const
    {
      require_started();
      return controller_->step_size();
    }

    /// The order the next attempt is taken at.
    int order() const
    {
      require_started();
      return std::min(controller_->order(),
        std::max(1, controller_->steps()));
    }

    /// Take one step the error test accepts, retrying it smaller when the solve
    /// or the test refuses it; it knows nothing of the output times.
    StepRecord advance()
    {
      require_started();
      return take_step(nullptr, nullptr);
    }

    /// Integrate over the output times, which are the start time or at/after
    /// the initialized time and increase, and report the store-mode values.
    RunReport run(std::span<const double> output_times, RunHooks<State> hooks = {})
    {
      require_started();
      if (controller_->steps() != 0)
        throw std::logic_error("atstep: run requires a fresh start");
      detail::StepSchedule schedule(
        settings_, output_times, start_time_, accepted_.time, end_time_);

      RunReport report;
      double time = accepted_.time;
      report.final_time = time;

      if (schedule.stores_initial())
        store(start_time_, accepted_.state, false, hooks, report);

      while (time < end_time_) {
        const StepRecord record = take_step(&schedule, &report);
        time = record.time;
        ++report.steps;
        report.max_order = std::max(report.max_order, record.order);
        report.final_time = time;
        if (hooks.on_step)
          hooks.on_step(record, accepted_.state);
        schedule.for_each_store(time, [&](double output_time, bool interpolated) {
          if (not interpolated) {
            store(output_time, accepted_.state, false, hooks, report);
            return;
          }
          interpolate(output_time, output_.interpolated);
          store(output_time, output_.interpolated, true, hooks, report);
        });
      }
      return report;
    }

    /// Write the dense-output value (derivative=0) or first derivative
    /// (derivative=1) at `time` into `out`. The value is the BDF polynomial
    /// through the last `order + 1` levels: inside the last step it is
    /// interpolated, older times are extrapolated. At the newest level the
    /// stored derivative is the solve derivative.
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
        slope = slope * gamma
          + coefficient / accepted_.psi[static_cast<std::size_t>(j - 1)];
        coefficient *= gamma;
        detail::axpy(out,
          derivative == 0 ? coefficient : slope,
          accepted_.phi[static_cast<std::size_t>(j)]);
        if (j < accepted_.order)
          gamma = (delta + accepted_.psi[static_cast<std::size_t>(j - 1)])
            / accepted_.psi[static_cast<std::size_t>(j)];
      }
    }

  private:
    struct Trial {
      double step = 0.0;
      double time = 0.0;
      int order = 1;
      int iterations = 0;
      detail::StepErrors errors;
      double cj = 0.0;
      double ck = 0.5;
      std::array<double, detail::max_levels> beta {};
      std::array<double, detail::max_levels> psi {};
      std::array<double, detail::max_levels> sigma {};
    };

    struct History {
      State state {};
      State derivative {};
      double time = 0.0;
      int order = 1;
      std::array<State, detail::max_levels> phi {};
      std::array<double, detail::max_levels> psi {};
    };

    struct Work {
      State predictor {};
      State derivative_predictor {};
      State solution {};
      State correction {};
      State scratch {};
    };

    /// The scaled absolute tolerance follows the last state returned to the
    /// caller. A level is the accepted state itself; an interpolated value is
    /// one buffer, overwritten by the next store.
    struct Output {
      const State* reference = nullptr;
      State interpolated {};
    };

    SolveOutcome solve(double time,
      double step,
      int order,
      double cj,
      const State& previous,
      const State& reference,
      bool initializing = false)
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

    StepRecord take_step(detail::StepSchedule* schedule, RunReport* report)
    {
      const double from = accepted_.time;
      auto& controller = *controller_;
      for (;;) {
        double step = controller.step_size();
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
        const double next_time = from + step;
        if (not std::isfinite(step) or not std::isfinite(next_time)
          or not(next_time > from))
          throw std::invalid_argument(
            "atstep: step must be finite, positive and advance time");

        Trial trial;
        trial.step = step;
        trial.time = next_time;
        trial.order = std::min(controller.order(), std::max(1, controller.steps()));
        controller.begin_attempt(step, trial.order);
        prepare(trial);

        const SolveOutcome outcome = solve(
          next_time, step, trial.order, trial.cj, accepted_.state,
          *output_.reference);
        trial.iterations = outcome.iterations;
        if (not outcome.converged) {
          if (report)
            ++report->failed_solves;
          controller.reject_convergence(step);
          if (step <= controller.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from)
              + " did not converge at the smallest step " + std::to_string(step));
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
                "atstep: the time step at t = " + std::to_string(from)
                + " fell to " + std::to_string(controller.step_size())
                + " and the local error is still "
                + std::to_string(trial.errors.test));
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

    /// Update history in place. State operations must not throw; the accepted
    /// tables are written before the controller and output are updated.
    void accept(const Trial& trial)
    {
      work_.scratch = work_.derivative_predictor;
      detail::axpy(work_.scratch, trial.cj, work_.correction);
      accepted_.state = work_.solution;
      accepted_.derivative = work_.scratch;
      accepted_.time = trial.time;
      accepted_.order = trial.order;
      const auto k = static_cast<std::size_t>(trial.order);
      for (std::size_t j = 1; j <= k; ++j)
        detail::axpy(accepted_.phi[j], trial.beta[j] - 1.0, accepted_.phi[j]);
      accepted_.phi[k + 1] = work_.correction;
      detail::axpy(accepted_.phi[k], 1.0, work_.correction);
      for (std::size_t j = k; j > 0; --j)
        detail::axpy(accepted_.phi[j - 1], 1.0, accepted_.phi[j]);
      accepted_.psi = trial.psi;
      controller_->accept(trial.errors, trial.step);
      output_.reference = &accepted_.state;
      update_norm_step_index();
    }

    void require_started() const
    {
      if (not started_)
        throw std::logic_error(
          "atstep: no valid run; call start before accessing or advancing it");
    }

    /// Zero the table the configured order can read and seed the first
    /// difference of the initialized derivative: `phi[1] = h * y'`.
    void reset_table(
      const State& level, double first_step, const State& derivative)
    {
      for (int j = 0; j <= settings_.max_order + 1; ++j) {
        State& column = accepted_.phi[static_cast<std::size_t>(j)];
        column = level;
        detail::axpy(column, -1.0, level);
      }
      accepted_.phi[0] = level;
      detail::axpy(accepted_.phi[1], first_step, derivative);
      accepted_.psi.fill(0.0);
      accepted_.psi[0] = first_step;
    }

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
        detail::axpy(work_.scratch,
          trial.beta[k - 1],
          accepted_.phi[k - 1]);
        errors.below2 = trial.sigma[k - 2] * error_norm(work_.scratch);
      }
      work_.scratch = work_.correction;
      detail::axpy(work_.scratch, -1.0, accepted_.phi[k + 1]);
      errors.above = error_norm(work_.scratch)
        / static_cast<double>(trial.order + 2);
      return errors;
    }

    /// Configure a default norm from the run's settings. A norm of the
    /// caller's, passed to the three-argument constructor, is used as it is.
    void configure_default_norm()
    {
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>) {
        norm_.relative = settings_.relative_tolerance;
        norm_.absolute = settings_.relative_tolerance
          * settings_.absolute_tolerance_factor;
        norm_.scaled = settings_.scale_absolute_tolerance;
        norm_.field_sizes = settings_.field_sizes;
      }
    }

    void update_norm_step_index()
    {
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        norm_.step_index = controller_
          ? static_cast<std::size_t>(controller_->steps())
          : 0;
    }

    double error_norm(const State& error)
    {
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        return norm_(error, accepted_.state, *output_.reference);
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
      if (&state == &accepted_.state)
        output_.reference = &accepted_.state;
      else {
        if (&state != &output_.interpolated)
          output_.interpolated = state;
        output_.reference = &output_.interpolated;
      }
      if (hooks.on_store)
        hooks.on_store(time, *output_.reference, interpolated);
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
