// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "settings.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>

namespace atstep {

  /// The weights `c` with `sum_j c[j] * v[j]` equal to the value at `x` of
  /// the polynomial through `(nodes[j], v[j])`.
  std::array<double, max_levels>
  lagrange_value(std::span<const double> nodes, double x);

  /// The weights `c` with `sum_j c[j] * v[j]` equal to the `order`-th
  /// derivative at `x` of the polynomial through `(nodes[j], v[j])`.
  std::array<double, max_levels> lagrange_derivative(
    std::span<const double> nodes, double x, int order);

  namespace detail {

    inline void validate_nodes(std::span<const double> nodes, double x)
    {
      if (nodes.empty() or nodes.size() > static_cast<std::size_t>(max_levels))
        throw std::invalid_argument("atstep: invalid interpolation node count");
      if (not std::isfinite(x))
        throw std::invalid_argument("atstep: interpolation time must be finite");
      for (std::size_t index = 0; index < nodes.size(); ++index) {
        if (not std::isfinite(nodes[index]))
          throw std::invalid_argument("atstep: interpolation nodes must be finite");
        for (std::size_t i = 0; i < index; ++i)
          if (nodes[index] == nodes[i])
            throw std::invalid_argument("atstep: interpolation nodes must be distinct");
      }
    }

    /// The coefficients of `prod_{i != skip} (u - nodes[i])` in the
    /// monomial basis of `u`, up to the degree `count - 1`, and that
    /// product's value at `nodes[skip]` — the denominator of the Lagrange
    /// basis polynomial.
    struct Basis {
      std::array<double, max_levels> coefficients {};
      double denominator = 1.0;
    };

    inline Basis
    basis_polynomial(std::span<const double> nodes, int count, int skip)
    {
      Basis basis;
      basis.coefficients[0] = 1.0;
      int degree = 0;
      for (int i = 0; i < count; ++i) {
        if (i == skip)
          continue;
        const double node = nodes[static_cast<std::size_t>(i)];
        for (int m = degree + 1; m >= 1; --m)
          basis.coefficients[static_cast<std::size_t>(m)] = basis.coefficients[static_cast<std::size_t>(m - 1)] - node * basis.coefficients[static_cast<std::size_t>(m)];
        basis.coefficients[0] = -node * basis.coefficients[0];
        ++degree;
        basis.denominator *= nodes[static_cast<std::size_t>(skip)] - node;
      }
      return basis;
    }

  } // namespace detail

  inline std::array<double, max_levels>
  lagrange_value(std::span<const double> nodes, double x)
  {
    return lagrange_derivative(nodes, x, 0);
  }

  inline std::array<double, max_levels> lagrange_derivative(
    std::span<const double> nodes, double x, int order)
  {
    detail::validate_nodes(nodes, x);
    const int count = static_cast<int>(nodes.size());
    std::array<double, max_levels> weights {};
    if (order < 0 or order >= count)
      return weights;
    for (int j = 0; j < count; ++j) {
      const detail::Basis basis = detail::basis_polynomial(nodes, count, j);
      // The order-th derivative of the monomial expansion at x is
      //   sum_{m >= order} c_m * m!/(m - order)! * x^(m - order).
      double value = 0.0;
      double falling = 1.0;
      for (int i = 1; i <= order; ++i)
        falling *= static_cast<double>(i);
      double power = 1.0;
      for (int m = order; m < count; ++m) {
        if (m > order) {
          falling *= static_cast<double>(m) / static_cast<double>(m - order);
          power *= x;
        }
        value += basis.coefficients[static_cast<std::size_t>(m)] * falling * power;
      }
      weights[static_cast<std::size_t>(j)] = value / basis.denominator;
    }
    return weights;
  }

  /// The differentiation weights of one BDF step, and the coefficients its
  /// error estimate is formed with.
  ///
  /// The state is known at the levels `t_{n+1}` (the level being solved for) and
  /// `t_{n+1-j}`, `j = 1..k`, and `h_j = t_{n+2-j} - t_{n+1-j}` is the step that
  /// reached the level `j - 1` from the level `j`. The p-th time derivative at
  /// `t_{n+1}` of the polynomial through those `k + 1` levels is
  ///
  ///     y^(p)(t_{n+1}) = sum_{j=0..k} weights(p)[j] * y_{n+1-j},   p = 0..d,
  ///
  /// with `weights(0) = {1, 0, ..., 0}` (the polynomial reproduces the level
  /// itself) and `weights(1)` the classical BDF formula. The weights of a
  /// variable step are what keeps a run of order `k` at order `k` when the step
  /// size changes: they differentiate that polynomial exactly, at any steps.
  ///
  /// Because a BDF step is implicit, a residual `F(t, y, y', y'')` assembled
  /// with these weights is linear in the unknown `y_{n+1}` through
  /// `weights(p)[0]`.
  ///
  /// The error-test coefficients are IDA's, the solver COMSOL's BDF is: `sigma`
  /// scales the correction of a step into the error estimate at a given order,
  /// and `ck` is the scale of the test the step is actually held to.
  class StepFormula {
  public:
    /// The first-order formula of the step `h`: backward Euler, with the
    /// derivative weights `{1/h, -1/h}`.
    explicit StepFormula(double h = 1.0)
    {
      const double steps[1] = {h};
      fill(1, std::span<const double>(steps, 1), 1);
    }

    /// The formula of the step into `t_{n+1}`, of BDF order `order`, with
    /// the derivative levels `0..derivatives`.
    ///
    /// `steps[j - 1] = h_j` for `j = 1..k + 1`, the newest first:
    /// `steps[0]` is the step into `t_{n+1}`, and `steps[j]` the step into
    /// `t_{n+2-j}`. The formula needs `derivatives` of them and the
    /// predictor of the step `order + 1`, so an order the steps do not
    /// reach is lowered to the one they do — a run whose history is short
    /// starts at a lower order rather than extrapolating on a polynomial it
    /// has no points for.
    StepFormula(int order,
      std::span<const double> steps,
      int derivatives = 1)
    {
      fill(order, steps, derivatives);
    }

    int order() const noexcept
    {
      return order_;
    }
    int derivatives() const noexcept
    {
      return derivatives_;
    }

    /// The weights of the derivative level `p`: `order() + 1` of them, the
    /// coefficient of `y_{n+1-j}` at index `j`.
    std::span<const double> weights(int derivative) const
    {
      if (derivative < 0 or derivative > derivatives_)
        throw std::invalid_argument("atstep: derivative level " + std::to_string(derivative) + " is outside the formula's 0.." + std::to_string(derivatives_));
      return std::span<const double>(
        weights_[static_cast<std::size_t>(derivative)].data(),
        static_cast<std::size_t>(order_) + 1);
    }

    /// The weights of the polynomial through the `order() + 1` levels
    /// *before* `t_{n+1}`, evaluated at `t_{n+1}`: the extrapolation the
    /// step starts from. Their sum is one.
    std::span<const double> predictor() const
    {
      return std::span<const double>(
        predictor_.data(), static_cast<std::size_t>(predictor_count_));
    }

    /// The weights of the polynomial through the `levels` most recent
    /// levels before `t_{n+1}`, evaluated at `t_{n+1}`: the estimate at a
    /// lower order is the correction against one level fewer. `levels` is
    /// clamped to what the formula's steps reach.
    std::array<double, max_levels> extrapolation(int levels) const
    {
      std::array<double, max_levels> weights {};
      if (levels > available_)
        levels = available_;
      if (levels <= 1) {
        weights[0] = 1.0;
        return weights;
      }
      std::array<double, max_levels> nodes {};
      for (int j = 0; j < levels; ++j)
        nodes[static_cast<std::size_t>(j)] = -psi_[static_cast<std::size_t>(j + 1)];
      return lagrange_value(
        std::span<const double>(nodes.data(),
          static_cast<std::size_t>(levels)),
        0.0);
    }

    /// IDA's `sigma[j]`: the factor that turns the correction of a step
    /// into the error estimate at the order `j`.
    double sigma(int order) const noexcept
    {
      return order >= 0 and order <= max_order
        ? sigma_[static_cast<std::size_t>(order)]
        : 0.0;
    }

    /// IDA's `ck`: the factor the correction of this step is held to.
    double ck() const noexcept
    {
      return ck_;
    }

    /// The weight `y^(1)` carries on `y_{n+1}`: the `cj` of IDA.
    double leading() const noexcept
    {
      return weights_[1][0];
    }

  private:
    void fill(int order, std::span<const double> steps, int derivatives)
    {
      const int available = static_cast<int>(std::min<std::size_t>(steps.size(), max_order + 1));
      if (available < 1)
        throw std::invalid_argument(
          "atstep: a step formula needs at least one step size");
      for (int j = 0; j < available; ++j)
        if (not std::isfinite(steps[static_cast<std::size_t>(j)])
          or not(steps[static_cast<std::size_t>(j)] > 0.0))
          throw std::invalid_argument("atstep: step size " + std::to_string(j) + " is not positive");

      // The formula reads `order` steps, the predictor `order + 1`.
      order_ = std::clamp(order, 1, std::max(1, available - 1));
      derivatives_ = std::clamp(derivatives, 0, max_derivative);

      const int count = order_ + 1;
      // The predictor is the polynomial through the levels *before*
      // t_{n+1}, of which there are one more than the scheme reads; a
      // history of a single level has no polynomial, and the predictor is
      // that level: the constant extrapolation a first step starts from.
      available_ = available;
      predictor_count_ = std::min(count, available);

      // psi[i] = t_{n+1} - t_{n+1-i}: the distance back to the level i. A
      // step the history does not reach — the levels before the first one
      // of a run — is taken to be the oldest step it does, which is what
      // IDA's own initialization does with the step before the first.
      for (int i = 1; i <= std::max(count, available); ++i) {
        const int known = i - 1 < available ? i - 1 : available - 1;
        psi_[static_cast<std::size_t>(i)] = psi_[static_cast<std::size_t>(i - 1)] + steps[static_cast<std::size_t>(known)];
      }

      // The derivative weights: the levels, evaluated at t_{n+1}.
      std::array<double, max_levels> nodes {};
      for (int j = 0; j < count; ++j)
        nodes[static_cast<std::size_t>(j)] = -psi_[static_cast<std::size_t>(j)];
      for (int p = 0; p <= derivatives_; ++p)
        weights_[static_cast<std::size_t>(p)] = lagrange_derivative(
          std::span<const double>(nodes.data(),
            static_cast<std::size_t>(count)),
          0.0,
          p);

      // The predictor: the polynomial through the levels before t_{n+1} —
      // one more level than the scheme reads — evaluated at t_{n+1}.
      if (predictor_count_ == 1) {
        predictor_[0] = 1.0;
      }
      else {
        std::array<double, max_levels> older {};
        for (int j = 0; j < predictor_count_; ++j)
          older[static_cast<std::size_t>(j)] = -psi_[static_cast<std::size_t>(j + 1)];
        predictor_ = lagrange_value(
          std::span<const double>(
            older.data(),
            static_cast<std::size_t>(predictor_count_)),
          0.0);
      }

      // IDA's coefficients of the error test (see IDASetCoeffs). Its
      // `psi[i]` is the distance back to the level `i` — psi[0] being the
      // step into t_n — while psi_[0] here is the level being solved for,
      // so IDA's index runs one ahead of this one.
      sigma_[0] = 1.0;
      for (int i = 1; i <= order_; ++i) {
        const double alpha_i = steps[0] / psi_[static_cast<std::size_t>(i + 1)];
        sigma_[static_cast<std::size_t>(i)] = static_cast<double>(i) * sigma_[static_cast<std::size_t>(i - 1)] * alpha_i;
      }
      double alphas = 0.0;
      double alpha0 = 0.0;
      for (int i = 0; i < order_; ++i) {
        alphas -= 1.0 / static_cast<double>(i + 1);
        alpha0 -= (i == 0 ? 1.0
                          : steps[0] / psi_[static_cast<std::size_t>(i + 1)]);
      }
      const double alpha_k = steps[0] / psi_[static_cast<std::size_t>(order_ + 1)];
      ck_ = std::abs(alpha_k + alphas - alpha0);
      if (alpha_k > ck_)
        ck_ = alpha_k;
    }

    int order_ = 1;
    int derivatives_ = 1;
    int predictor_count_ = 1;
    int available_ = 1;
    std::array<double, max_levels> psi_ {};
    std::array<std::array<double, max_levels>, max_derivative + 1>
      weights_ {};
    std::array<double, max_levels> predictor_ {};
    std::array<double, max_levels + 1> sigma_ {};
    double ck_ = 0.5;
  };

} // namespace atstep
