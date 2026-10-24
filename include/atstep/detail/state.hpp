// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace atstep::detail {

  /// The default arithmetic of a `std::vector<double>`, declared before the
  /// concept so that a plain vector models a state without an adapter.
  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x);

  /// `y += a * x`, the only arithmetic the stepper does on a state. `x` may be
  /// `y`, which rescales a state. The default forwards to a member `axpy`.
  template <class State>
    requires requires(State& y, const State& x, double a) { y.axpy(a, x); }
  void axpy(State& y, double a, const State& x)
  {
    y.axpy(a, x);
  }

  /// A default-constructible, copyable value with the arithmetic above. State
  /// operations must not throw: the stepper publishes history in place.
  template <class State>
  concept StateVector
    = std::is_copy_assignable_v<State> && std::is_copy_constructible_v<State>
    && std::is_default_constructible_v<State>
    && requires(State& y, const State& x, double a) {
         atstep::detail::axpy(y, a, x);
       };

  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x)
  {
    if (y.size() != x.size())
      throw std::invalid_argument("atstep: axpy state sizes differ");
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] += a * x[i];
  }

} // namespace atstep::detail
