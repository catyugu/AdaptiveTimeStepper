// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace atstep {

  /// The default arithmetic of a `std::vector<double>`, declared before the
  /// concept so that a plain vector models a state without an adapter.
  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x);

  /// The customization point that makes a state a vector of the stepper:
  ///
  ///     y += a * x
  ///
  /// It is the only arithmetic the stepper does on a state. Everything else
  /// it needs — the operator, the matrix, the assembly — is the caller's, and
  /// a state of any shape works: a vector of one value, a block vector of a
  /// coupled system, a distributed vector with a layout of its own.
  ///
  /// `x` may be `y` itself, which rescales a state: the stepper scales stored
  /// columns this way rather than copying them. An implementation must
  /// therefore read an element of `x` before it writes the same element of
  /// `y`, as `y += a * x` elementwise does.
  ///
  /// The default forwards to a member `axpy(a, x)`. A type that has no such
  /// member can be wrapped in a type providing that member:
  ///
  ///     struct State {
  ///         std::vector<float> values;
  ///         void axpy(double a, const State& x);
  ///     };
  template <class State>
    requires requires(State& y, const State& x, double a) { y.axpy(a, x); }
  void axpy(State& y, double a, const State& x)
  {
    y.axpy(a, x);
  }

  /// The states the stepper advances: a copyable value the arithmetic above
  /// is defined on. A state is held by value, so a run copies it once per
  /// level.
  template <class State>
  concept StateVector = std::is_copy_assignable_v<State> && std::is_copy_constructible_v<State> && std::is_default_constructible_v<State> && requires(State& y, const State& x, double a) { atstep::axpy(y, a, x); };

  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x)
  {
    if (y.size() != x.size())
      throw std::invalid_argument("atstep: axpy state sizes differ");
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] += a * x[i];
  }

  /// Weighted RMS: sqrt(mean((error_i / (A_i + relative*abs(y_i)))^2)).
  /// With scaled=false, A_i=absolute. Otherwise, for accepted-step index j:
  ///
  ///   A_i = absolute * ((2-exp(-j/5))*abs(returned_i) + exp(-j/5)*sigma_i).
  ///
  /// y is the accepted internal state; returned is the most recently stored
  /// state, possibly interpolated. Missing scale entries use sigma_i=1.
  /// States expose size() and operator[]. Other layouts can supply a custom
  /// norm as TimeStepper<State, Norm>'s second template parameter.
  struct WeightedRmsNorm {
    double relative = 1e-2;
    double absolute = 1e-3;
    bool scaled = false;
    std::size_t step_index = 0;
    /// `sigma_i`; an index past the end is one.
    std::vector<double> scale;

    template <class State>
      requires requires(const State& s) {
        { s.size() } -> std::convertible_to<std::size_t>;
        { s[std::size_t {}] } -> std::convertible_to<double>;
      }
    double operator()(const State& error, const State& reference) const
    {
      return (*this)(error, reference, reference);
    }

    /// COMSOL updates the scaled absolute tolerance from the solution
    /// returned to the caller, which can be interpolated. The relative part
    /// continues to use the last accepted internal state.
    template <class State>
    double operator()(const State& error,
      const State& reference,
      const State& scale_reference) const
    {
      const std::size_t n = reference.size();
      if (error.size() != n or scale_reference.size() != n)
        throw std::invalid_argument("atstep: norm state sizes differ");
      if (n == 0)
        return 0.0;
      const double decay = std::exp(-0.2 * static_cast<double>(step_index));
      double sum = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        const double magnitude = std::abs(scale_reference[i]);
        const double sigma = i < scale.size() ? scale[i] : 1.0;
        const double limit
          = scaled ? absolute * ((2.0 - decay) * magnitude + decay * sigma)
                   : absolute;
        const double weight = limit + relative * std::abs(reference[i]);
        if (not std::isfinite(error[i]) or not std::isfinite(weight) or weight < 0.0)
          return std::numeric_limits<double>::infinity();
        if (weight == 0.0) {
          if (error[i] != 0.0)
            return std::numeric_limits<double>::infinity();
        }
        else {
          const double ratio = error[i] / weight;
          sum += ratio * ratio;
        }
      }
      return std::sqrt(sum / static_cast<double>(n));
    }
  };

} // namespace atstep
