// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/state.hpp"
#include "settings.hpp"
#include "stepper.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace atstep {

  /// Displacement and velocity: the solved unknowns of a second-order run.
  template <class State>
  struct SecondOrderState {
    State value {};
    State derivative {};

    void axpy(double a, const SecondOrderState& x)
    {
      detail::axpy(value, a, x.value);
      detail::axpy(derivative, a, x.derivative);
    }
  };

  /// Custom error norm over the full pair. The default measures displacement / sqrt(2).
  template <class State>
  using SecondOrderNorm = std::function<double(const SecondOrderState<State>& error,
    const SecondOrderState<State>& accepted,
    const SecondOrderState<State>& returned,
    std::size_t accepted_steps)>;

  /// What a second-order run reports as it goes.
  template <class State>
  struct SecondOrderHooks {
    /// Every accepted step, with the state and the first derivative it left.
    /// The states are the stepper's own buffers and are only valid during the
    /// call.
    std::function<void(const StepRecord&, const State&, const State&)> on_step;
    /// Stored time, displacement, first derivative and interpolation flag.
    /// Interpolated derivatives come from dense output. References last only for the call.
    std::function<void(double, const State&, const State&, bool)> on_store;
  };

  /// Second-order integration over `(u, v = u')` using the first-order BDF kernel.
  /// The callback solves both components; the default norm measures displacement
  /// error / sqrt(2). See doc/design.md for the formulation and roundoff rationale.
  template <class State>
  class SecondOrderStepper {
    static_assert(detail::StateVector<State>,
      "atstep: State must be default-constructible, copyable and support axpy");

    using Pair = SecondOrderState<State>;
    using Inner = TimeStepper<Pair>;

  public:
    /// Solve both unknowns and return convergence. For `M u'' + C u' + K u = f`:
    ///
    ///     cj*u - v = -known.value
    ///     K*u + (cj*M + C)*v = f - M*known.derivative
    using Solve = std::function<bool(const StepRequest<Pair>&, Pair&)>;
    using Norm = SecondOrderNorm<State>;

    SecondOrderStepper(const SecondOrderStepper&) = delete;
    SecondOrderStepper& operator=(const SecondOrderStepper&) = delete;
    SecondOrderStepper(SecondOrderStepper&&) = delete;
    SecondOrderStepper& operator=(SecondOrderStepper&&) = delete;

    /// Default displacement weighted RMS / sqrt(2); requires size() and operator[].
    static Norm default_norm(const ErrorTest& error)
    {
      if constexpr (detail::Indexable<State>) {
        return [weights = detail::WeightedRmsNorm::from(error)](
                 const Pair& estimate,
                 const Pair& accepted,
                 const Pair& returned,
                 std::size_t accepted_steps) {
          return std::sqrt(0.5)
            * weights(estimate.value, accepted.value, returned.value, accepted_steps);
        };
      }
      else {
        throw std::invalid_argument(
          "atstep: this state has no size() and operator[], so the run needs a norm");
      }
    }

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one step of the pair (see `Solve`).
    /// @param[in] norm the norm the error of a step is measured in, over the
    /// pair. Empty means `default_norm(settings.error)`.
    explicit SecondOrderStepper(const Settings& settings,
      Solve solve,
      Norm norm = {})
      : inner_(settings, std::move(solve), resolve_norm(std::move(norm), settings.error))
    {
    }

    /// Start from displacement and velocity, with two initialization probes.
    /// With initialization off, the pair derivative starts at `(v, 0)`.
    void start(double time, const Pair& initial, double end_time)
    {
      if (inner_.settings().initialization == Initialization::off) {
        Pair derivative {};
        derivative.value = initial.derivative;
        derivative.derivative = initial.derivative;
        detail::axpy(derivative.derivative, -1.0, initial.derivative);
        inner_.start(time, initial, derivative, end_time);
        return;
      }
      inner_.start(time, initial, end_time);
    }

    /// Prepare the run from a pair and its time derivative, `(u', u'')`, as a
    /// restart from a recorded pair does. The initialization setting is not read:
    /// the derivative is the caller's, so a probe has none to measure.
    void start(double time, const Pair& initial, const Pair& derivative, double end_time)
    {
      inner_.start(time, initial, derivative, end_time);
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
    /// for, with the derivative each of them pairs with. `output_times` is any
    /// range of `double` that increases from the start time to the end time.
    template <class Times, std::enable_if_t<detail::TimeRange<Times>, int> = 0>
    RunReport run(const Times& output_times,
      SecondOrderHooks<State> hooks = {})
    {
      RunHooks<Pair> inner_hooks;
      if (hooks.on_step)
        inner_hooks.on_step = [&hooks](const StepRecord& record, const Pair& y) {
          hooks.on_step(record, y.value, y.derivative);
        };
      if (hooks.on_store)
        inner_hooks.on_store = [this, &hooks](double time, const Pair& y, bool interpolated) {
          if (not interpolated) {
            hooks.on_store(time, y.value, y.derivative, false);
            return;
          }
          // An interpolated value's derivative is the dense output's, not the
          // pair's own velocity there: `u'` of the polynomial that produced the
          // value is what COMSOL reports (doc/comsol-alignment.md#second-order-error-and-output).
          inner_.interpolate(time, dense_, 1);
          hooks.on_store(time, y.value, dense_.value, true);
        };
      return inner_.run(output_times, std::move(inner_hooks));
    }

    /// The counts of the run so far.
    const RunReport& report() const noexcept
    {
      return inner_.report();
    }

    /// Write displacement (derivative=0) or its first derivative (derivative=1).
    /// Between nodes the derivative comes from the displacement polynomial,
    /// which can differ from the interpolated velocity component. `out` is caller-owned.
    void interpolate(double time, State& out, int derivative = 0)
    {
      inner_.interpolate(time, dense_, derivative);
      out = dense_.value;
    }

  private:
    /// The norm a run uses: the caller's, or COMSOL's over the pair.
    static Norm resolve_norm(Norm norm, const ErrorTest& error)
    {
      if (norm)
        return norm;
      return default_norm(error);
    }

    Inner inner_;
    /// The pair the dense output is written into before one half of it is handed
    /// to the caller: the kernel interpolates both halves together.
    Pair dense_ {};
  };

} // namespace atstep
