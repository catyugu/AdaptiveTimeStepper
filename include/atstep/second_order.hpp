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
#include <type_traits>
#include <utility>

namespace atstep {

  namespace detail {

    /// Displacement and velocity, advanced by the same BDF formula. The
    /// default error norm measures only displacement, as COMSOL's IDA path
    /// does; the solve enforces the velocity's kinematic constraint.
    template <class State> struct Augmented {
      State value {};
      State derivative {};

      void axpy(double a, const Augmented& x)
      {
        detail::axpy(value, a, x.value);
        detail::axpy(derivative, a, x.derivative);
      }

      std::size_t size() const
      {
        return value.size();
      }

      double operator[](std::size_t i) const
      {
        return value[i];
      }
    };

    /// Apply a caller-supplied norm to displacement only.
    template <class State, class Norm> struct ValueNorm {
      Norm norm;

      double operator()(const Augmented<State>& error,
        const Augmented<State>& reference) const
      {
        return norm(error.value, reference.value);
      }
    };

    /// Owned by the adapter and overwritten on each attempt. The displacement
    /// predictor, first history and derivative scratch belong to the kernel.
    template <class State> struct SecondOrderWork {
      State first_predictor {};
      State second_predictor {};
      State second_history {};
    };

  } // namespace detail

  template <class State, class Norm>
    requires detail::StateVector<State>
  class SecondOrderStepper;

  /// Coefficients for a solve in displacement alone.
  ///
  /// The equation is `M u'' + C u' + K u = f(t)`, or any residual in `u`, `u'`
  /// and `u''`, and the step holds its derivatives to the affine forms
  ///
  ///     u'  = cj    u + known(1),   cj     = derivative_weight(),
  ///     u'' = cj^2  u + known(2),   known(d) = derivative_history(d).
  ///
  /// Substitution gives the condensed linear system
  ///
  ///     (K + cj C + cj^2 M) u_{n+1} = f(t_{n+1}) - M known(2) - C known(1),
  ///
  /// or the nonlinear Jacobian `dF/du + cj dF/du' + cj^2 dF/du''`.
  ///
  /// `derivative(u, out, d)` evaluates `u^(d)` in correction form,
  /// `u^(d)_predictor + weight(d) (u - predictor)`, which is small where the
  /// predictor is good and free of the direct form's cancellation.
  ///
  /// The request and its state references are valid only during the solve
  /// callback. Do not retain pointers or references to them.
  template <class State> class SecondOrderRequest {
  public:
    SecondOrderRequest(const SecondOrderRequest&) = delete;
    SecondOrderRequest& operator=(const SecondOrderRequest&) = delete;
    SecondOrderRequest(SecondOrderRequest&&) = delete;
    SecondOrderRequest& operator=(SecondOrderRequest&&) = delete;

    const double time;
    const double step;
    const int order;
    const std::size_t step_index;
    const bool initializing;

    const State& absolute_tolerance_reference() const noexcept
    {
      return step_.absolute_tolerance_reference().value;
    }

    /// The weight `u^(derivative)` carries on the state being solved for:
    /// `cj` for the first derivative, `cj^2` for the second.
    double derivative_weight(int derivative = 1) const
    {
      require_supported(derivative);
      const double cj = step_.derivative_weight();
      return derivative == 1 ? cj : cj * cj;
    }

    /// The known term of the derivative level `derivative` of the affine
    /// forms above.
    const State& derivative_history(int derivative = 1) const
    {
      require_supported(derivative);
      return derivative == 1 ? step_.derivative_history().value
                             : work_.second_history;
    }

    /// `out = u^(derivative)` of the displacement `u`. `out` is a state of the
    /// caller's and may not be one of this request's.
    void derivative(const State& u, State& out, int derivative = 1) const
    {
      require_supported(derivative);
      State& scratch = step_.derivative_scratch_.value;
      scratch = u;
      detail::axpy(scratch, -1.0, step_.predictor_.value);
      out = derivative == 1 ? work_.first_predictor : work_.second_predictor;
      const double cj = step_.derivative_weight();
      detail::axpy(out, derivative == 1 ? cj : cj * cj, scratch);
    }

  private:
    template <class S, class N>
      requires detail::StateVector<S>
    friend class SecondOrderStepper;

    SecondOrderRequest(
      const StepRequest<detail::Augmented<State>>& step_request,
      const detail::SecondOrderWork<State>& work)
      : time(step_request.time)
      , step(step_request.step)
      , order(step_request.order)
      , step_index(step_request.step_index)
      , initializing(step_request.initializing)
      , step_(step_request)
      , work_(work)
    {
    }

    static void require_supported(int derivative)
    {
      if (derivative < 1 or derivative > 2)
        throw std::invalid_argument(
          "atstep: a second-order request supports the first and second "
          "derivatives");
    }

    const StepRequest<detail::Augmented<State>>& step_;
    const detail::SecondOrderWork<State>& work_;
  };

  /// What a second-order run reports as it goes.
  template <class State> struct SecondOrderHooks {
    /// Every accepted step, with the displacement and the velocity it left.
    /// The states are the stepper's own buffers and are only valid during the
    /// call.
    std::function<void(const StepRecord&, const State&, const State&)>
      on_step;
    /// Every stored value: the time, the displacement, the velocity, and
    /// whether the value was interpolated from the steps around it. The states
    /// are the stepper's own buffers and are only valid during the call.
    std::function<void(double, const State&, const State&, bool)> on_store;
  };

  /// Integrates `(displacement, velocity)` through TimeStepper, with each
  /// solve condensed to displacement (see SecondOrderRequest). The error test
  /// measures displacement only. `start(t0, u0, v0, end)` applies the kernel's
  /// two artificial backward-Euler steps to both initial conditions.
  ///
  /// An integrator is neither copyable nor movable.
  template <class State, class Norm = detail::WeightedRmsNorm>
    requires detail::StateVector<State>
  class SecondOrderStepper {
    using InnerNorm = std::conditional_t<std::is_same_v<Norm, detail::WeightedRmsNorm>,
      Norm,
      detail::ValueNorm<State, Norm>>;
    using Inner = TimeStepper<detail::Augmented<State>, InnerNorm>;

  public:
    using Solve
      = std::function<SolveOutcome(const SecondOrderRequest<State>&, State&)>;

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one condensed step (see the class
    /// comment).
    /// @param[in] norm the norm the error of a step is measured in, over the
    /// displacement. The default is COMSOL's weighted RMS, whose tolerances
    /// are taken from `settings`.
    explicit SecondOrderStepper(const Settings& settings,
      Solve solve,
      Norm norm = Norm {})
      : solve_(std::move(solve))
      , inner_(settings,
          [this](const StepRequest<detail::Augmented<State>>& request,
            detail::Augmented<State>& y) { return solve_step(request, y); },
          inner_norm(std::move(norm)))
    {
      if (not solve_)
        throw std::invalid_argument("atstep: a run needs a solve");
    }

    const Settings& settings() const noexcept
    {
      return inner_.settings();
    }
    double start_time() const noexcept
    {
      return inner_.start_time();
    }
    double end_time() const noexcept
    {
      return inner_.end_time();
    }

    /// Prepare the run: the displacement and the velocity at `time`,
    /// integrated to `end_time`.
    void start(double time,
      const State& initial,
      const State& velocity,
      double end_time)
    {
      inner_.start(
        time, detail::Augmented<State> {initial, velocity}, end_time);
    }

    /// The time the run has reached: the newest accepted level's.
    double time() const
    {
      return inner_.time();
    }

    /// The current displacement, and the velocity the level pairs it with.
    /// They are the stepper's own buffers, overwritten by the next step.
    const State& state() const
    {
      return inner_.state().value;
    }
    const State& velocity() const
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

    /// Take one step the error test accepts, from the newest level, and
    /// report it.
    StepRecord advance()
    {
      return inner_.advance();
    }

    /// Integrate from the start time to the end time, and report the values
    /// the store mode asks for, with the velocity each of them pairs with.
    RunReport run(std::span<const double> output_times,
      SecondOrderHooks<State> hooks = {})
    {
      RunHooks<detail::Augmented<State>> inner_hooks;
      if (hooks.on_step)
        inner_hooks.on_step = [&hooks](const StepRecord& record,
                                 const detail::Augmented<State>& y) {
          hooks.on_step(record, y.value, y.derivative);
        };
      if (hooks.on_store)
        inner_hooks.on_store = [&hooks](double time,
                                  const detail::Augmented<State>& y,
                                  bool interpolated) {
          hooks.on_store(time, y.value, y.derivative, interpolated);
        };
      return inner_.run(output_times, std::move(inner_hooks));
    }

    /// The dense-output displacement (derivative=0) or velocity
    /// (derivative=1). At the newest level the velocity is the derivative the
    /// solve left, which agrees with `velocity()` to round-off; between the
    /// levels it is the derivative of the dense-output polynomial, which is not
    /// the pair's own `u'` there. The result is kept in the stepper: another
    /// call can overwrite it.
    const State& interpolate(double time, int derivative = 0)
    {
      return inner_.interpolate(time, derivative).value;
    }

  private:
    static InnerNorm inner_norm(Norm norm)
    {
      if constexpr (std::is_same_v<Norm, detail::WeightedRmsNorm>)
        return norm;
      else
        return InnerNorm {std::move(norm)};
    }

    /// Turn the kernel's step of the pair into the condensed step of the
    /// displacement: name the three affine forms, hand the solve the
    /// displacement alone, and derive the velocity the level pairs with it.
    SolveOutcome solve_step(
      const StepRequest<detail::Augmented<State>>& request,
      detail::Augmented<State>& y)
    {
      const double cj = request.derivative_weight();
      const State& known1 = request.derivative_history().value;
      work_.second_history = request.derivative_history().derivative;
      detail::axpy(work_.second_history, cj, known1);
      work_.first_predictor = known1;
      detail::axpy(work_.first_predictor, cj, y.value);
      work_.second_predictor = work_.second_history;
      detail::axpy(work_.second_predictor, cj * cj, y.value);
      const SecondOrderRequest<State> condensed(request, work_);
      const SolveOutcome outcome = solve_(condensed, y.value);
      if (outcome.converged) {
        // Restore the kinematic constraint only for a solved candidate.
        y.derivative = known1;
        detail::axpy(y.derivative, cj, y.value);
      }
      return outcome;
    }

    Solve solve_;
    detail::SecondOrderWork<State> work_;
    Inner inner_;
  };

} // namespace atstep
