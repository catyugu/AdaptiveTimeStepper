// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
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
  /// Assignment should reuse storage when dimensions are unchanged.
  template <class State>
  concept StateVector = std::is_copy_assignable_v<State> && std::is_copy_constructible_v<State> && std::is_default_constructible_v<State> && requires(State& y, const State& x, double a) { atstep::detail::axpy(y, a, x); };

  /// A state whose dofs can be read one by one: what the default norm needs and
  /// a state of its own arithmetic may not have.
  template <class State>
  concept Indexable = requires(const State& s) {
    { s.size() } -> std::convertible_to<std::size_t>;
    { s[std::size_t {}] } -> std::convertible_to<double>;
  };

  inline void
  axpy(std::vector<double>& y, double a, const std::vector<double>& x)
  {
    if (y.size() != x.size())
      throw std::invalid_argument("atstep: axpy state sizes differ");
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] += a * x[i];
  }

  /// COMSOL's local error test: the weighted RMS of the error estimate over the
  /// state, averaged over the dependent variables first and over the degrees of
  /// freedom of a variable second.
  ///
  /// For a variable `j` with `N_j` degrees of freedom, its scale `sigma_j` and
  /// its error estimate `E`:
  ///
  ///     Y_i     = returned_i / sigma_j
  ///     A_s,i   = A_k (beta_j + rms(Y_j))     automatic scale, update on
  ///     A_s,i   = A_k                          otherwise
  ///     beta_j  = (1 - exp(-j/5)) max(abs(Y_j)) + exp(-j/5)
  ///     W_i     = A_s,i + R abs(Y_i)           scaled method
  ///     W_i     = A_k   + R abs(U_i)           unscaled method
  ///     E^2     = (1/M) sum_j (1/N_j) sum_i (E_i / W_i)^2
  ///
  /// where `j` is the accepted-step index and `A_k` is the entered absolute
  /// tolerance. Measured: doc/evidence/tolerance-axes.md. A single variable of
  /// one dof at scale one reduces to COMSOL's scalar form, which is the only
  /// case the recorded trajectories cover.
  ///
  /// States expose `size()` and `operator[]`; a layout that does not is given
  /// its norm by its caller (the `Norm` parameter of `TimeStepper`).
  struct WeightedRmsNorm {
    /// The relative tolerance `R`.
    double relative = 1e-2;

    /// The absolute tolerance `A_k`: a factor of `R`, or a value.
    double absolute = 1e-3;

    /// COMSOL's "Tolerance method": the criterion on `Y = U / scale`, or on `U`.
    bool scaled = false;

    /// The variables of the state, as their sizes: COMSOL averages the squared
    /// ratios over the variables first and over the degrees of freedom of one
    /// second. An empty vector is one variable over every degree of freedom.
    std::vector<std::size_t> variables;

    /// `sigma_i`, the scale of the variable each degree of freedom belongs to.
    std::vector<double> scale;

    /// Whether each variable's absolute tolerance follows the solution: the
    /// `A_s` branch above. An empty vector means every variable follows it,
    /// which is COMSOL's default.
    std::vector<bool> scaled_tolerance;

    /// The accepted steps the run has behind it: COMSOL's `j`.
    std::size_t accepted_steps = 0;

    /// The norm of COMSOL's settings. An empty `variables` is one variable of
    /// scale one over the whole state, whose size the caller supplies.
    static WeightedRmsNorm from(const ErrorTest& error)
    {
      WeightedRmsNorm norm;
      norm.relative = error.relative;
      norm.absolute = error.absolute.kind == AbsoluteTolerance::Kind::factor
        ? error.relative * error.absolute.value
        : error.absolute.value;
      norm.scaled = error.scaled;
      for (const Variable& variable : error.variables) {
        if (variable.dofs == 0)
          throw std::invalid_argument(
            "atstep: a variable of the error test needs at least one degree of freedom");
        norm.variables.push_back(variable.dofs);
        norm.scale.insert(norm.scale.end(), variable.dofs, variable.scale);
        norm.scaled_tolerance.push_back(variable.scaled_tolerance);
      }
      return norm;
    }

    /// COMSOL updates the scaled absolute tolerance from the solution returned
    /// to the caller, which can be interpolated. The relative part and the
    /// variable's own magnitude continue to use the last accepted internal
    /// state.
    template <Indexable State>
    double operator()(const State& error,
      const State& reference,
      const State& returned) const
    {
      const std::size_t n = reference.size();
      if (error.size() != n or returned.size() != n)
        throw std::invalid_argument("atstep: norm state sizes differ");
      if (n == 0)
        return 0.0;
      // An empty list is one variable over the whole state.
      const std::size_t count = variables.empty() ? 1 : variables.size();
      if (not variables.empty()) {
        std::size_t dofs = 0;
        for (const std::size_t variable : variables)
          dofs += variable;
        if (dofs != n)
          throw std::invalid_argument(
            "atstep: the variables of the error test do not add up to the size of the state");
      }
      // The scales and the partition, walked together.
      const auto scale_at = [this](std::size_t i) {
        return i < scale.size() ? scale[i] : 1.0;
      };
      const double decay = std::exp(-0.2 * static_cast<double>(accepted_steps));
      double total = 0.0;
      std::size_t first = 0;
      for (std::size_t variable = 0; variable < count; ++variable) {
        const std::size_t dofs = variables.empty() ? n : variables[variable];
        const bool update = scaled_tolerance.empty()
          or (variable < scaled_tolerance.size() and scaled_tolerance[variable]);
        // The variable's own norm, from the solution the caller was given.
        double square = 0.0;
        double largest = 0.0;
        for (std::size_t i = first; i < first + dofs; ++i) {
          const double y = returned[std::size_t {i}] / scale_at(std::size_t {i});
          square += y * y;
          largest = std::max(largest, std::abs(y));
        }
        const double beta = (largest > 0.0 and largest < 1.0)
          ? (1.0 - decay) * largest + decay
          : 1.0;
        const double limit = scaled
          ? (update ? absolute * (beta + std::sqrt(square / static_cast<double>(dofs)))
                    : absolute)
          : absolute;
        // The variable's ratios, averaged over its own degrees of freedom.
        double sum = 0.0;
        for (std::size_t i = first; i < first + dofs; ++i) {
          const double s = scale_at(std::size_t {i});
          // The relative part follows the last accepted internal state; only
          // the absolute one follows the solution the caller was given.
          const double reference_value = scaled ? reference[i] / s : reference[i];
          const double estimate = scaled ? error[i] / s : error[i];
          const double weight = limit + relative * std::abs(reference_value);
          if (not std::isfinite(estimate) or not std::isfinite(weight) or weight < 0.0)
            return std::numeric_limits<double>::infinity();
          if (weight == 0.0) {
            if (estimate != 0.0)
              return std::numeric_limits<double>::infinity();
            continue;
          }
          const double ratio = estimate / weight;
          sum += ratio * ratio;
        }
        // COMSOL averages over the variables, so that a variable of many
        // degrees of freedom does not outvote one of few.
        total += sum / static_cast<double>(dofs);
        first += dofs;
      }
      return std::sqrt(total / static_cast<double>(count));
    }
  };

} // namespace atstep::detail
