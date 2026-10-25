// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "../settings.hpp"

namespace atstep::detail {

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

  /// A default-constructible, copyable value with the arithmetic above.
  /// Assignment should reuse storage when dimensions are unchanged. State
  /// operations may throw; if a history update fails, start must rebuild the
  /// run before it can be read or advanced again.
  template <class State>
  concept StateVector = std::is_copy_assignable_v<State> && std::is_copy_constructible_v<State> && std::is_default_constructible_v<State> && requires(State& y, const State& x, double a) { atstep::detail::axpy(y, a, x); };

  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x)
  {
    if (y.size() != x.size())
      throw std::invalid_argument("atstep: axpy state sizes differ");
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] += a * x[i];
  }

  /// Weighted RMS: sqrt(mean((error_i / (A_i + relative*abs(y_i)))^2)).
  /// `A_i` is `absolute` everywhere when `scaled` is false, and when it is true,
  /// for accepted-step index j:
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
    /// The scaled form of COMSOL's absolute tolerance, which blends the magnitude
    /// of the returned solution with `scale`. `TimeStepper` sets it from
    /// `ErrorTest`; the plain form is what a norm on its own starts with.
    bool scaled = false;
    /// `sigma_i`; an empty vector, and an index past its end, are one.
    std::vector<double> scale;
    /// The partition of the state into fields, as their sizes: COMSOL averages
    /// the squared ratios over the fields first and over the degrees of freedom
    /// of a field second. An empty vector is one field over every degree of
    /// freedom.
    std::vector<std::size_t> fields;
    std::size_t step_index = 0;

    /// The norm of COMSOL's settings: the tolerances, and the scales of the
    /// fields laid over the degrees of freedom in order. An empty layout is one
    /// field of scale one over every degree of freedom.
    static WeightedRmsNorm from(const ErrorTest& error)
    {
      WeightedRmsNorm norm;
      norm.relative = error.relative;
      norm.absolute = error.absolute.kind == AbsoluteTolerance::Kind::factor
        ? error.relative * error.absolute.value
        : error.absolute.value;
      norm.scaled = error.scaled;
      for (const Field& field : error.layout.fields) {
        if (field.dofs == 0)
          throw std::invalid_argument(
            "atstep: a field of the error test needs at least one degree of freedom");
        norm.fields.push_back(field.dofs);
        norm.scale.insert(norm.scale.end(), field.dofs, field.scale);
      }
      return norm;
    }

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
      std::vector<double> squares(n, 0.0);
      for (std::size_t i = 0; i < n; ++i) {
        double limit = absolute;
        if (scaled) {
          const double magnitude = std::abs(scale_reference[i]);
          const double sigma = i < scale.size() ? scale[i] : 1.0;
          limit = absolute * ((2.0 - decay) * magnitude + decay * sigma);
        }
        const double weight = limit + relative * std::abs(reference[i]);
        if (not std::isfinite(error[i]) or not std::isfinite(weight) or weight < 0.0)
          return std::numeric_limits<double>::infinity();
        if (weight == 0.0) {
          if (error[i] != 0.0)
            return std::numeric_limits<double>::infinity();
          continue;
        }
        const double ratio = error[i] / weight;
        squares[i] = ratio * ratio;
      }
      if (fields.empty())
        return std::sqrt(std::accumulate(squares.begin(), squares.end(), 0.0)
          / static_cast<double>(n));
      // COMSOL averages over the fields first, so that a field of many degrees
      // of freedom does not outvote a field of few.
      double total = 0.0;
      std::size_t first = 0;
      for (const std::size_t dofs : fields) {
        const std::size_t last = std::min(first + dofs, n);
        total += std::accumulate(squares.begin() + static_cast<std::ptrdiff_t>(first),
                   squares.begin() + static_cast<std::ptrdiff_t>(last),
                   0.0)
          / static_cast<double>(dofs);
        first = last;
      }
      return std::sqrt(total / static_cast<double>(fields.size()));
    }
  };

} // namespace atstep::detail
