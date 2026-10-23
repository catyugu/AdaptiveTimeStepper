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
  /// The scheme reads the last `order + 1` of them, the predictor one more,
  /// and the error estimate at the order two below three fewer, so a run
  /// keeps `max_levels` of them.
  template <class State>
  class History {
  public:
    void reset(double time, const State& state)
    {
      states_.clear();
      times_.clear();
      push(time, state);
    }

    void push(double time, const State& state)
    {
      states_.push_front(state);
      times_.push_front(time);
    }

    void pop()
    {
      if (not states_.empty()) {
        states_.pop_front();
        times_.pop_front();
      }
    }

    std::size_t size() const noexcept
    {
      return states_.size();
    }
    bool empty() const noexcept
    {
      return states_.empty();
    }

    /// The time of the level `level`, 0 being the newest.
    double time(std::size_t level) const
    {
      return times_[level];
    }

    /// The state of the level `level`, 0 being the newest.
    const State& state(std::size_t level) const
    {
      return states_[level];
    }

    /// The times, newest first.
    std::vector<double> times() const
    {
      return std::vector<double>(times_.begin(), times_.end());
    }

    /// Forget the levels past the newest `levels`.
    void keep(std::size_t levels)
    {
      while (states_.size() > levels) {
        states_.pop_back();
        times_.pop_back();
      }
    }

  private:
    std::deque<State> states_;
    std::deque<double> times_;
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
      if (std::holds_alternative<BdfDerivative>(derivative_)) {
        require_first_derivative(derivative);
        return std::get<BdfDerivative>(derivative_).leading;
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
      if (std::holds_alternative<BdfDerivative>(derivative_)) {
        require_first_derivative(derivative);
        return std::get<BdfDerivative>(derivative_).history;
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
      if (std::holds_alternative<BdfDerivative>(derivative_)) {
        require_first_derivative(derivative);
        const auto& bdf = std::get<BdfDerivative>(derivative_);
        State correction = y;
        atstep::axpy(correction, -1.0, bdf.predictor);
        State result = bdf.derivative_predictor;
        atstep::axpy(result, bdf.leading, correction);
        return result;
      }
      State result = derivative_history(derivative);
      atstep::axpy(result, derivative_weight(derivative), y);
      return result;
    }

    /// Views are valid for the solve callback; the stepper owns their
    /// states.
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
    bool converged = true;
    /// The nonlinear iterations the solve took, for the log.
    int iterations = 0;
  };

  /// One attempted step: the state it produced, and the error estimate of it.
  struct Attempt {
    bool solved = false;
    bool passed = false;
    double time = 0.0;
    double step = 0.0;
    int order = 0;
    int iterations = 0;
    StepErrors errors;
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
    /// Every accepted step, in order.
    std::function<void(const StepRecord&)> on_step;
    /// Every stored value: the time, the state, and whether it was
    /// interpolated from the steps around it rather than being a level.
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
  ///     const auto report = stepper.run(times, hooks);
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
      if (not(end_time > time))
        throw std::invalid_argument(
          "atstep: the end time must be past the start time");
      start_time_ = time;
      end_time_ = end_time;
      const double span = end_time - time;
      controller_ = StepController(settings_, span);
      started_ = false;
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
        for (State& column : phi_) {
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
    /// pushes the state onto the history, `reject` drops it.
    Attempt attempt(double step)
    {
      require_started();
      Attempt result;
      result.step = step;
      result.time = history_.time(0) + step;

      // Both formulations use the same history bound on startup order.
      const int available = static_cast<int>(
        std::min<std::size_t>(history_.size(), max_levels));
      const int order = std::min(controller_.order(), std::max(1, available - 1));
      controller_.begin_attempt(step, order);
      result.order = order;

      auto derivative = prepare_derivative(step, order);
      work_ = predictor_;
      StepRequest<State> request {
        result.time,
        step,
        order,
        std::move(derivative),
        history_,
        scale_reference_,
        static_cast<std::size_t>(controller_.steps())};
      const SolveOutcome outcome = solve_(request, work_);
      result.iterations = outcome.iterations;
      if (not outcome.converged)
        return result;
      result.solved = true;

      if (fixed_leading() or controller_.error_controlled()) {
        ee_ = work_;
        atstep::axpy(ee_, -1.0, predictor_);
      }
      if (not controller_.error_controlled()) {
        result.passed = true;
        return result;
      }
      result.errors = estimate(order);
      result.passed = result.errors.passed();
      return result;
    }

    /// Commit the step `attempt` produced.
    void accept(const Attempt& attempt)
    {
      if (fixed_leading()) {
        last_derivative_ = derivative_predictor_;
        atstep::axpy(last_derivative_, cj_, ee_);
      }
      history_.push(attempt.time, work_);
      history_.keep(static_cast<std::size_t>(max_levels));
      last_order_ = attempt.order;
      if (fixed_leading()) {
        const auto k = static_cast<std::size_t>(attempt.order);
        trial_phi_[k + 1] = ee_;
        atstep::axpy(trial_phi_[k], 1.0, ee_);
        for (std::size_t j = k; j > 0; --j)
          atstep::axpy(trial_phi_[j - 1], 1.0, trial_phi_[j]);
        phi_ = trial_phi_;
        psi_ = trial_psi_;
      }
      controller_.accept(attempt.errors, attempt.step);
      scale_reference_ = work_;
      if constexpr (std::is_same_v<Norm, WeightedRmsNorm>)
        norm_.step_index = static_cast<std::size_t>(controller_.steps());
      if (not fixed_leading() and controller_.error_controlled()) {
        previous_ee_ = ee_;
        previous_order_ = attempt.order;
        has_previous_ee_ = true;
      }
    }

    /// Drop the step `attempt` produced: it missed the tolerance.
    void reject(const Attempt& attempt)
    {
      controller_.reject(attempt.errors, attempt.step);
    }

    /// Drop the step `attempt` produced: its solve did not converge.
    void reject_convergence(const Attempt& attempt)
    {
      controller_.reject_convergence(attempt.step);
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
        Attempt piece = attempt(step);
        const double boundary = intermediate_landing
          ? intermediate_target
          : (settings_.step_mode == StepMode::strict ? target
                                                     : end_time_);
        if ((settings_.step_mode == StepMode::strict or intermediate_landing or not settings_.interpolate_end_time) and std::abs(piece.time - boundary) <= 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(time) + std::abs(boundary)))
          piece.time = boundary;
        if (not piece.solved) {
          ++report.failed_solves;
          if (piece.step <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the solve of the step from t = " + std::to_string(from) + " did not converge at the "
                                                                                 "smallest step "
              + std::to_string(piece.step));
          reject_convergence(piece);
          continue;
        }
        if (not piece.passed) {
          ++report.rejected;
          reject(piece);
          if (controller_.step_size() <= controller_.smallest_step())
            throw std::runtime_error(
              "atstep: the time step at t = " + std::to_string(from) + " fell to " + std::to_string(controller_.step_size()) + " and the local error is still " + std::to_string(piece.errors.test) + "; the tolerance cannot be met on this interval");
          continue;
        }

        accept(piece);
        time = piece.time;
        ++report.steps;
        report.max_order = std::max(report.max_order, piece.order);
        if (report.steps == 1)
          report.first_step = piece.step;
        report.last_step = piece.step;
        report.final_time = time;
        if (hooks.on_step)
          hooks.on_step(StepRecord {time,
            piece.step,
            piece.order,
            piece.iterations,
            piece.errors.test});

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
        derivative,
        static_cast<int>(count));
      combine(interpolated_,
        std::span<const double>(weights.data(), count));
      return interpolated_;
    }

  private:
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

    void prepare_fixed_leading(double step, int order)
    {
      trial_phi_ = phi_;
      trial_psi_ = psi_;
      if (controller_.steps() == 0) {
        trial_psi_[0] = step;
        trial_phi_[1] = derivative_;
        atstep::axpy(trial_phi_[1], step / probe_ - 1.0, derivative_);
      }
      predictor_ = phi_[0];
      derivative_history_ = phi_[0];
      atstep::axpy(derivative_history_, -1.0, phi_[0]);
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
        State column = trial_phi_[j];
        atstep::axpy(trial_phi_[j], beta - 1.0, column);
        atstep::axpy(predictor_, 1.0, trial_phi_[j]);
        atstep::axpy(derivative_history_, gamma, trial_phi_[j]);
      }
      trial_psi_[static_cast<std::size_t>(order)] = distance;
      ck_ = std::max(alpha, std::abs(alpha - harmonic - alpha0));
      cj_ = harmonic / step;
      derivative_predictor_ = derivative_history_;
      atstep::axpy(derivative_history_, -cj_, predictor_);
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
          atstep::axpy(scratch_, 1.0, trial_phi_[k]);
          errors.below = sigma_[k - 1] * error_norm(scratch_);
        }
        if (order >= 3) {
          atstep::axpy(scratch_, 1.0, trial_phi_[k - 1]);
          errors.below2 = sigma_[k - 2] * error_norm(scratch_);
        }
        scratch_ = ee_;
        atstep::axpy(scratch_, -1.0, trial_phi_[k + 1]);
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
    std::array<State, max_levels> phi_ {};
    std::array<State, max_levels> trial_phi_ {};
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
