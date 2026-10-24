// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "detail/state.hpp"
#include "stepper.hpp"

#include <cstddef>

namespace atstep {

  /// The pair a second-order-in-time equation is integrated as:
  /// `(u, v = u')`. `M u'' + C u' + K u = f` becomes
  ///
  ///     u' - v            = 0
  ///     M v' + C v + K u  = f,
  ///
  /// the first-order system COMSOL's IDA path forms. The pair is a state of
  /// `TimeStepper`; `size()` and `operator[]` expose the displacement, so the
  /// default error norm measures the displacement alone. A custom norm that
  /// reads `value` and `derivative` can include the velocity.
  template <class State>
  struct SecondOrderState {
    State value {};
    State derivative {};

    void axpy(double a, const SecondOrderState& x)
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

  /// `known2 = cj known1 + known(v)`, the known term of
  /// `u'' = cj^2 u + known2`. `known2` is the caller's buffer, reused across
  /// steps.
  template <class State>
  void second_derivative_history(
    const StepRequest<SecondOrderState<State>>& request, State& known2)
  {
    known2 = request.derivative_history().derivative;
    detail::axpy(
      known2, request.derivative_weight(), request.derivative_history().value);
  }

  /// Enforce the kinematic constraint at the solved state:
  /// `v = cj u + known1`. Call it after the displacement solve succeeded; a
  /// failed solve must leave the pair alone so the step can be retried.
  template <class State>
  void set_velocity(SecondOrderState<State>& state,
    double cj,
    const State& known1)
  {
    state.derivative = known1;
    detail::axpy(state.derivative, cj, state.value);
  }

} // namespace atstep
