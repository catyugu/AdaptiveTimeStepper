// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/state.hpp"
#include "settings.hpp"
#include "stepper.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <utility>

namespace atstep {

  namespace detail {

    /// The state and its first derivative, advanced by the same BDF formula: the
    /// pair the kernel integrates for an equation that is second order in time.
    template <class State>
    struct Augmented {
      State value {};
      State derivative {};

      void axpy(double a, const Augmented& x)
      {
        detail::axpy(value, a, x.value);
        detail::axpy(derivative, a, x.derivative);
      }
    };

  } // namespace detail

  /// Coefficients for a solve in the state alone, of an equation that is second
  /// order in time: `M u'' + C u' + K u = f(t)`, or any residual in `u`, `u'`
  /// and `u''`.
  ///
  /// The step holds the two derivatives to the affine forms
  ///
  ///     u'  = derivative_weight(1) u + derivative_history(1),
  ///     u'' = derivative_weight(2) u + derivative_history(2),
  ///
  /// with `derivative_weight(2) == derivative_weight(1)^2`, which leaves the
  /// condensed system `(K + cj C + cj^2 M) u = f - M known(2) - C known(1)` and
  /// the nonlinear Jacobian `dF/du + cj dF/du' + cj^2 dF/du''`.
  ///
  /// The stepper builds the request for one solve callback and its state
  /// references are borrowed for that call; it cannot be copied or retained.
  template <class State>
  class SecondOrderRequest {
  public:
    SecondOrderRequest(const SecondOrderRequest&) = delete;
    SecondOrderRequest& operator=(const SecondOrderRequest&) = delete;
    SecondOrderRequest(SecondOrderRequest&&) = delete;
    SecondOrderRequest& operator=(SecondOrderRequest&&) = delete;

    const double time;
    const double step;
    const int order;
    /// The accepted steps behind the run.
    const std::size_t accepted_steps;
    const bool initializing;

    /// The weight `u^(derivative)` carries on the state being solved for.
    double derivative_weight(int derivative = 1) const
    {
      require_supported(derivative);
      return derivative == 1 ? weight_ : weight_ * weight_;
    }

    /// The known term of the derivative's affine form.
    const State& derivative_history(int derivative = 1) const
    {
      require_supported(derivative);
      return derivative == 1 ? known1_ : known2_;
    }

    SecondOrderRequest(double new_time,
      double step_size,
      int step_order,
      std::size_t steps,
      bool initialization,
      double weight,
      const State& known1,
      const State& known2)
      : time(new_time)
      , step(step_size)
      , order(step_order)
      , accepted_steps(steps)
      , initializing(initialization)
      , weight_(weight)
      , known1_(known1)
      , known2_(known2)
    {
    }

  private:
    static void require_supported(int derivative)
    {
      if (derivative < 1 or derivative > 2)
        throw std::invalid_argument(
          "atstep: a second-order request supports the first and second derivatives");
    }

    double weight_ = 0.0;
    const State& known1_;
    const State& known2_;
  };

  /// What a second-order run reports as it goes.
  template <class State>
  struct SecondOrderHooks {
    /// Every accepted step, with the state and the first derivative it left.
    /// The states are the stepper's own buffers and are only valid during the
    /// call.
    std::function<void(const StepRecord&, const State&, const State&)> on_step;
    /// Every stored value: the time, the state, its first derivative, and
    /// whether the value was interpolated from the steps around it. The states
    /// are the stepper's own buffers and are only valid during the call.
    std::function<void(double, const State&, const State&, bool)> on_store;
  };

  /// `M u'' + C u' + K u = f` integrated by the first-order kernel over the pair
  /// `(u, v = u')`, with every solve condensed back to the state: COMSOL's IDA
  /// path. The error test measures the state alone — COMSOL excludes the
  /// auxiliary velocity except that its documentation includes it while the
  /// consistent initialization is on, which is still an open measurement
  /// (doc/refactor-plan.md, E05).
  ///
  /// The integrator is neither copyable nor movable.
  template <class State>
    requires detail::StateVector<State>
  class SecondOrderStepper {
    using Pair = detail::Augmented<State>;
    using Inner = TimeStepper<Pair>;

  public:
    using Solve
      = std::function<SolveOutcome(const SecondOrderRequest<State>&, State&)>;
    using Norm = ErrorNorm<State>;

    SecondOrderStepper(const SecondOrderStepper&) = delete;
    SecondOrderStepper& operator=(const SecondOrderStepper&) = delete;
    SecondOrderStepper(SecondOrderStepper&&) = delete;
    SecondOrderStepper& operator=(SecondOrderStepper&&) = delete;

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one condensed step (see the class comment).
    /// @param[in] norm the norm the error of a step is measured in, over the
    /// state. Empty means COMSOL's, over the state alone.
    explicit SecondOrderStepper(const Settings& settings,
      Solve solve,
      Norm norm = {})
      : solve_(std::move(solve))
      , inner_(settings,
          [this](const StepRequest<Pair>& request, Pair& y) {
            return solve_step(request, y);
          },
          pair_norm(settings.error, std::move(norm)))
    {
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
    }

    /// Prepare the run: the state and its first derivative at `time`, integrated
    /// to `end_time`. The two artificial backward-Euler steps apply to both
    /// initial values.
    void start(double time, const State& initial, const State& derivative, double end_time)
    {
      inner_.start(time, Pair {initial, derivative}, end_time);
    }

    /// The time the run has reached: the newest accepted level's.
    double time() const
    {
      return inner_.time();
    }

    /// The current state, and the first derivative the level pairs it with.
    /// They are the stepper's own buffers, overwritten by the next step.
    const State& state() const
    {
      return inner_.state().value;
    }
    const State& derivative() const
    {
      return inner_.state().derivative;
    }

    /// The step size and order of the next step.
    double step_size() const
    {
      return inner_.step_size();
    }
    int order() const
    {
      return inner_.order();
    }

    /// Take the next step the error test accepts, cut to the end time when it is
    /// not interpolated.
    StepRecord advance()
    {
      return inner_.advance();
    }

    /// Whether the run has reached its end time.
    bool finished() const
    {
      return inner_.finished();
    }

    /// Integrate the rest of the run, reporting the values the store mode asks
    /// for, with the derivative each of them pairs with.
    RunReport run(std::span<const double> output_times,
      SecondOrderHooks<State> hooks = {})
    {
      RunHooks<Pair> inner_hooks;
      if (hooks.on_step)
        inner_hooks.on_step = [&hooks](const StepRecord& record, const Pair& y) {
          hooks.on_step(record, y.value, y.derivative);
        };
      if (hooks.on_store)
        inner_hooks.on_store = [&hooks](double time, const Pair& y, bool interpolated) {
          hooks.on_store(time, y.value, y.derivative, interpolated);
        };
      return inner_.run(output_times, std::move(inner_hooks));
    }

    /// The counts of the run so far.
    const RunReport& report() const noexcept
    {
      return inner_.report();
    }

    /// Write the dense-output state (derivative=0) or its first derivative
    /// (derivative=1) at `time` into `out`. The pair's dense output of order
    /// `derivative` is `(u, u')^(derivative)`, so its first half is the state's
    /// derivative of that order. At the newest level it is the derivative the
    /// solve left; between the levels it is the dense-output polynomial's, which
    /// need not be the pair's own `u'` there. `out` is the caller's, and may be
    /// the same state on every call.
    void interpolate(double time, State& out, int derivative = 0)
    {
      inner_.interpolate(time, dense_, derivative);
      out = dense_.value;
    }

  private:
    /// The norm the kernel applies to the pair: the caller's over the state, or
    /// COMSOL's, which measures the state and not the velocity the library added.
    static ErrorNorm<Pair> pair_norm(const ErrorTest& error, Norm norm)
    {
      if (not norm) {
        if constexpr (detail::Indexable<State>)
          norm = TimeStepper<State>::default_norm(error);
        else
          throw std::invalid_argument(
            "atstep: this state has no size() and operator[], so the run needs a norm");
      }
      return [norm](const Pair& estimate,
               const Pair& accepted,
               const Pair& returned,
               std::size_t accepted_steps) {
        return norm(estimate.value, accepted.value, returned.value, accepted_steps);
      };
    }

    /// Turn the kernel's step of the pair into the condensed step of the state:
    /// name the two affine forms, hand the solve the state alone, and derive the
    /// derivative the level pairs with from the kinematic constraint.
    SolveOutcome solve_step(const StepRequest<Pair>& request, Pair& y)
    {
      const double cj = request.derivative_weight();
      const State& known1 = request.derivative_history().value;
      // `known(2) = cj known(1) + known(v)`: the velocity's own known term,
      // shifted to the second derivative's level.
      known2_ = request.derivative_history().derivative;
      detail::axpy(known2_, cj, known1);
      const SecondOrderRequest<State> condensed {request.time, request.step,
        request.order, request.accepted_steps, request.initializing, cj, known1,
        known2_};
      const SolveOutcome outcome = solve_(condensed, y.value);
      if (outcome.converged) {
        // Restore the kinematic constraint only for a solved candidate.
        y.derivative = known1;
        detail::axpy(y.derivative, cj, y.value);
      }
      return outcome;
    }

    Solve solve_;
    State known2_ {};
    Inner inner_;
    /// The pair the dense output is written into before one half of it is handed
    /// to the caller: the kernel interpolates both halves together.
    Pair dense_ {};
  };

} // namespace atstep
