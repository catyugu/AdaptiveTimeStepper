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

    /// The state of an equation of the second order in time, written the way
    /// IDA integrates it: the displacement, and the velocity as a variable of
    /// its own.
    ///
    /// The two blocks are advanced by the same BDF formula, at the same order
    /// and step, and the step's solve keeps them consistent: the velocity is
    /// the derivative the formula gives the displacement. That is the state
    /// augmentation COMSOL documents for a second-order system under IDA —
    /// "extra variables are internally introduced so that it is possible to
    /// form a first-order time-derivative system".
    ///
    /// `size()` and `operator[]` are the variables the error test measures,
    /// and they are the displacement alone: the auxiliary velocity carries no
    /// independent error, being the derivative of the displacement's, and a
    /// norm over it would hold the step to the derivative's accuracy rather
    /// than the solution's. COMSOL's error test excludes the augmented
    /// velocity for the same reason.
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

    /// A norm of the caller's over the variables of the model, applied to the
    /// pair the augmentation carries.
    template <class State, class Norm> struct ValueNorm {
      Norm norm;

      double operator()(const Augmented<State>& error,
        const Augmented<State>& reference) const
      {
        return norm(error.value, reference.value);
      }
    };

    /// The coefficients of one condensed step, which the request reads (see
    /// `SecondOrderRequest`). The states are the stepper's; the known terms
    /// are formed once per attempt.
    template <class State> struct SecondOrderWork {
      /// The weight of the first derivative, which the second carries squared.
      double cj = 0.0;
      /// The known term of the first derivative, `u' = cj u + known1`.
      const State& known1;
      /// The known term of the second, `u'' = cj^2 u + known2`.
      const State& known2;
      /// The predictor of the displacement, and of its two derivatives.
      const State& predictor;
      const State& first_predictor;
      const State& second_predictor;
      /// The difference `u - predictor` of the correction form, borrowed from
      /// the stepper.
      State& scratch;
    };

  } // namespace detail

  /// What a step of a second-order-in-time equation is told about the step.
  ///
  /// The equation is `M u'' + C u' + K u = f(t)`, or any residual in `u`, `u'`
  /// and `u''`, and the step holds its derivatives to the affine forms
  ///
  ///     u'  = cj    u + known(1),   cj     = derivative_weight(),
  ///     u'' = cj^2  u + known(2),   known(d) = derivative_history(d).
  ///
  /// `cj` is the weight the BDF step gives the first derivative; the second
  /// carries it squared, because the velocity is itself the derivative of the
  /// displacement. Substituting the two into the equation leaves the condensed
  /// system of the displacement alone,
  ///
  ///     (K + cj C + cj^2 M) u_{n+1} = f(t_{n+1}) - M known(2) - C known(1),
  ///
  /// which is what COMSOL's IDA path solves: the velocity never enters the
  /// linear system, so a step costs a matrix of the displacement's size, not
  /// twice it. A nonlinear residual is condensed the same way, and its
  /// Jacobian is `dF/du + cj dF/du' + cj^2 dF/du''`.
  ///
  /// `derivative(u, out, d)` evaluates `u^(d)` in correction form,
  /// `u^(d)_predictor + weight(d) (u - predictor)`, which is small where the
  /// predictor is good and free of the direct form's cancellation.
  ///
  /// A request borrows the stepper's states. It is valid for the solve
  /// callback's own call and never after.
  template <class State> class SecondOrderRequest {
  public:
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
      return derivative == 1 ? work_.cj : work_.cj * work_.cj;
    }

    /// The known term of the derivative level `derivative` of the affine
    /// forms above.
    const State& derivative_history(int derivative = 1) const
    {
      require_supported(derivative);
      return derivative == 1 ? work_.known1 : work_.known2;
    }

    /// `out = u^(derivative)` of the displacement `u`. `out` is a state of the
    /// caller's and may not be one of this request's.
    void derivative(const State& u, State& out, int derivative = 1) const
    {
      require_supported(derivative);
      work_.scratch = u;
      detail::axpy(work_.scratch, -1.0, work_.predictor);
      out = derivative == 1 ? work_.first_predictor : work_.second_predictor;
      detail::axpy(out, derivative_weight(derivative), work_.scratch);
    }

  private:
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

  /// Advances an equation of the second order in time with COMSOL's adaptive
  /// BDF method.
  ///
  /// The equation is handed to the library as it is — `M u'' + C u' + K u =
  /// f(t)`, or any residual in `u`, `u'` and `u''` — and the caller's solve is
  /// a solve of the displacement alone. The library introduces the velocity
  /// the first-order kernel needs, integrates it with the same formula and
  /// keeps it consistent with the displacement, and gives the solve the
  /// condensed coefficients of the substitution (see `SecondOrderRequest`).
  /// Nothing of the velocity reaches the linear system, the settings or the
  /// tolerances; it is a companion state, not another unknown of the model.
  ///
  ///     atstep::SecondOrderStepper<Vector> stepper(settings,
  ///         [&](const atstep::SecondOrderRequest<Vector>& request, Vector& u) {
  ///             // u' = cj u + k1, u'' = cj^2 u + k2
  ///             solve((K + cj*C + cj*cj*M), u, f(t) - M*k2 - C*k1);
  ///             return atstep::SolveOutcome {true, iterations};
  ///         });
  ///     stepper.start(0.0, u0, v0, 10.0);  // displacement and velocity
  ///     const auto report = stepper.run(times, hooks);
  ///
  /// The initial velocity is the second initial condition, and `start` runs the
  /// same two artificial backward-Euler steps of consistent initialization as a
  /// first-order run does: the pair is advanced, the displacement's increment
  /// is what the first step is sized from, and the first returned state is at
  /// `t0 + 2 step_probe`, as it is there.
  ///
  /// The error test measures the displacement alone. A third or higher
  /// derivative in time is the caller's to write as a first-order system in
  /// `TimeStepper`; there is no formula for it here.
  ///
  /// An integrator is neither copyable nor movable.
  template <class State, class Norm = detail::WeightedRmsNorm>
    requires detail::StateVector<State>
  class SecondOrderStepper {
  public:
    using Solve
      = std::function<SolveOutcome(const SecondOrderRequest<State>&, State&)>;
    using InnerNorm = std::conditional_t<std::is_same_v<Norm, detail::WeightedRmsNorm>,
      Norm,
      detail::ValueNorm<State, Norm>>;
    using Inner = TimeStepper<detail::Augmented<State>, InnerNorm>;
    using Work = detail::SecondOrderWork<State>;

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
      work_.value_predictor = y.value;
      work_.second_history = request.derivative_history().derivative;
      detail::axpy(work_.second_history, cj, known1);
      work_.first_predictor = known1;
      detail::axpy(work_.first_predictor, cj, work_.value_predictor);
      work_.second_predictor = work_.second_history;
      detail::axpy(work_.second_predictor, cj * cj, work_.value_predictor);
      const Work view {cj,
        known1,
        work_.second_history,
        work_.value_predictor,
        work_.first_predictor,
        work_.second_predictor,
        work_.scratch};
      const SecondOrderRequest<State> condensed(request, view);
      const SolveOutcome outcome = solve_(condensed, y.value);
      // The kinematic constraint of the augmentation: the level's velocity is
      // the derivative the formula gives the level's displacement.
      y.derivative = known1;
      detail::axpy(y.derivative, cj, y.value);
      return outcome;
    }

    /// The states of one condensed step. `scratch` is the correction form's.
    struct Workspace {
      State value_predictor {};
      State first_predictor {};
      State second_predictor {};
      State second_history {};
      State scratch {};
    };

    Solve solve_;
    Workspace work_;
    Inner inner_;
  };

} // namespace atstep
