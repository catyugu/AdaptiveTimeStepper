// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
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
  /// The default forwards to a member `axpy(a, x)`. A type that has no such
  /// member gets one by overloading this function in this namespace:
  ///
  ///     namespace atstep {
  ///         inline void axpy(std::vector<float>& y, double a,
  ///                          const std::vector<float>& x)
  ///         {
  ///             for (std::size_t i = 0; i < y.size(); ++i) y[i] += a * x[i];
  ///         }
  ///     }
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
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] += a * x[i];
  }

  /// The norm a step's error is measured in: COMSOL's weighted root mean square,
  ///
  ///     ||e|| = sqrt( (1/N) sum_i ( e_i / (A_i + R |y_i|) )^2 ),
  ///
  /// with `R` the relative tolerance, `A_i` the absolute one, `N` the number of
  /// degrees of freedom, and `y` the state the error is taken against. It is the
  /// norm of COMSOL's error test, whose every variable of a coupled system
  /// shares one mean square: `N` counts all of them.
  ///
  /// `A_i` has COMSOL's two methods of its absolute tolerance. With `scaled` off
  /// it is `R * absolute_factor`: the plain reading of its "Factor" method, and
  /// the form a run of a physical field is held to. With it on it is its
  /// *scaled* absolute tolerance,
  ///
  ///     A_i = R * absolute_factor * (|y_i| + s_i),
  ///     s_i = (1 - e^(-j/5)) |y_i| + e^(-j/5) sigma_i,
  ///
  /// where `j` is the number of accepted time steps (zero at initialization) and
  /// `sigma_i` the scale of the degree of freedom — `scale` below, COMSOL's
  /// "scale for dependent variable", one by default. So the absolute tolerance
  /// follows the magnitude of each variable rather than a fixed floor, and it
  /// tightens as the run proceeds. At the first step `A_i` is
  /// `R * absolute_factor * (|y_i| + sigma_i)`; a scalar state of magnitude one
  /// is therefore held twice as tightly as the unscaled method holds it, and
  /// that factor is the difference between the first step of a run and COMSOL's
  /// own (see doc/comsol-alignment.md).
  ///
  /// A state is expected to expose its degrees of freedom as `size()` and
  /// `operator[]`. A distributed or blocked state passes a norm of its own
  /// instead — the stepper takes any callable object with this shape, so a
  /// reduction over the ranks is a lambda:
  ///
  ///     atstep::TimeStepper<State> stepper(settings, solve,
  ///         [](const State& e, const State& y) { return global_rms(e, y);
  ///         });
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
        if (weight > 0.0) {
          const double ratio = error[i] / weight;
          sum += ratio * ratio;
        }
      }
      return std::sqrt(sum / static_cast<double>(n));
    }
  };

} // namespace atstep
