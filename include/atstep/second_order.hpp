// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/state.hpp"
#include "settings.hpp"
#include "stepper.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <utility>

namespace atstep {

  /// The pair a second-order equation is integrated as: the state and its first
  /// derivative. It is public because it is the solve's unknown and the error
  /// test's argument, not an implementation detail.
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

  /// The error norm of a second-order run: COMSOL's weighted RMS over both
  /// halves of the pair, averaged over the two.
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
    /// Every stored value: the time, the state, its first derivative, and
    /// whether the value was interpolated from the steps around it. An
    /// interpolated value's derivative is the dense output's, which is what
    /// COMSOL reports there. The states are the stepper's own buffers and are
    /// only valid during the call.
    std::function<void(double, const State&, const State&, bool)> on_store;
  };

  /// `M u'' + C u' + K u = f(t)` as the first-order system COMSOL's IDA path
  /// integrates,
  ///
  ///     u' = v,     M v' + C v + K u = f,
  ///
  /// over the pair `(u, v)`: the kernel's state is `SecondOrderState`, and the
  /// caller's solve receives the pair's two affine forms
  ///
  ///     u' = cj u + known.value,
  ///     v' = cj v + known.derivative,
  ///
  /// to assemble the step for both halves. The velocity is a solved unknown
  /// here, not something this class reconstructs from the displacement: a
  /// reconstruction is a divided difference, and its round-off is the state's
  /// divided by the step, so it doubles as the step shrinks until a tight
  /// absolute tolerance on the velocity cannot be met at all. With the pair
  /// solved, each unknown comes from the equation that determines it best, and
  /// what accuracy is left over belongs to the caller's solver.
  ///
  /// The error test measures the state alone, as COMSOL's does: the velocity is
  /// carried by the system but left out of the local error, while the average
  /// over the variables the augmented system solves for still counts it, so the
  /// state's weighted norm enters divided by `sqrt(2)`
  /// (doc/evidence/second-order.md).
  ///
  /// The integrator is neither copyable nor movable.
  template <class State>
    requires detail::StateVector<State>
  class SecondOrderStepper {
    using Pair = SecondOrderState<State>;
    using Inner = TimeStepper<Pair>;

  public:
    /// The solve of one step of the pair: the caller fills both halves, and
    /// reports whether its nonlinear iteration converged. For a linear
    /// `M u'' + C u' + K u = f` the step is
    ///
    ///     cj u - v = -known.value
    ///     u + (cj M + C) v = f - M known.derivative
    ///
    /// and `v` is read from the second row, whose right-hand side is `O(cj |v|)`
    /// rather than `O(|u|)`; reading it back from `cj u + known.value` instead
    /// throws away the digits the velocity has.
    using Solve = std::function<bool(const StepRequest<Pair>&, Pair&)>;
    using Norm = SecondOrderNorm<State>;

    SecondOrderStepper(const SecondOrderStepper&) = delete;
    SecondOrderStepper& operator=(const SecondOrderStepper&) = delete;
    SecondOrderStepper(SecondOrderStepper&&) = delete;
    SecondOrderStepper& operator=(SecondOrderStepper&&) = delete;

    /// The norm COMSOL's settings describe: the weighted RMS of the state alone,
    /// divided by the two variables the augmented system solves for. COMSOL adds
    /// the velocity to the system its integrator carries but leaves it out of the
    /// error test, while the average over `M` still counts it, so with `M = 2` the
    /// state's norm enters as `rms / sqrt(2)`. A state must expose its dofs.
    static Norm default_norm(const ErrorTest& error)
    {
      if constexpr (detail::Indexable<State>) {
        detail::WeightedRmsNorm weights = detail::WeightedRmsNorm::from(error);
        return [weights](const Pair& estimate,
                 const Pair& accepted,
                 const Pair& returned,
                 std::size_t accepted_steps) mutable {
          weights.accepted_steps = accepted_steps;
          return std::sqrt(0.5) * weights(estimate.value, accepted.value, returned.value);
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

    /// Prepare the run from the pair at `time`, integrated to `end_time`, with
    /// COMSOL's consistent initialization: two artificial backward-Euler steps
    /// that reconcile both initial values with the equations.
    void start(double time, const Pair& initial, double end_time)
    {
      inner_.start(time, initial, end_time);
    }

    /// Prepare the run from a pair and its time derivative, `(u', u'')`, without
    /// the consistent initialization: COMSOL's "Consistent initialization: Off".
    /// The derivative completes the affine form, so the probes have nothing to
    /// derive.
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
        inner_hooks.on_store = [this, &hooks](double time, const Pair& y, bool interpolated) {
          if (not interpolated) {
            hooks.on_store(time, y.value, y.derivative, false);
            return;
          }
          // An interpolated value's derivative is the dense output's, not the
          // pair's own velocity there: `u'` of the polynomial that produced the
          // value is what COMSOL reports (doc/evidence/second-order.md).
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
