// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace atstep {

  /// COMSOL's weighted RMS over fields and degrees of freedom:
  ///
  ///     sqrt( (1/M) * sum_j (1/N_j) * sum_i (E_i / (A_i + R abs(Y_i)))^2 )
  ///
  /// with `M` fields, `N_j` dofs in field `j`, `R` the relative tolerance and
  /// `A_i` the absolute tolerance. A field whose size is not given by
  /// `field_sizes` is the whole state.
  ///
  /// With scaled absolute tolerance, for a field of scale reference `S`:
  ///
  ///     norm2 = sqrt(sum_i S_i^2)
  ///     norm_inf = max_i abs(S_i)
  ///     d = exp(-j/5)          j = accepted steps
  ///     beta = (1-d) norm_inf + d   if 0 < norm_inf < 1
  ///     beta = 1                    otherwise
  ///     A_i = absolute * (beta + norm2)
  ///
  /// and `A_i = absolute` when unscaled. The documented COMSOL form of the
  /// scalar case, `absolute * ((2-d) abs(S) + d)`, is evaluated directly so it
  /// keeps its rounding.
  struct WeightedRmsNorm {
    double relative = 1e-2;
    double absolute = 1e-3;
    bool scaled = false;
    /// Accepted steps, `j` in the scaled absolute tolerance.
    std::size_t step_index = 0;
    /// Dof count of each field. Empty is one field of the state's size.
    std::vector<std::size_t> field_sizes;

    template <class State>
      requires requires(const State& s) {
        { s.size() } -> std::convertible_to<std::size_t>;
        { s[std::size_t {}] } -> std::convertible_to<double>;
      }
    double operator()(const State& error, const State& reference) const
    {
      return (*this)(error, reference, reference);
    }

    /// `reference` is the accepted internal state, `scale_reference` the last
    /// state returned to the caller; COMSOL's scaled tolerance follows the
    /// latter, its relative part the former.
    template <class State>
    double
    operator()(const State& error, const State& reference, const State& scale_reference) const
    {
      const std::size_t n = reference.size();
      if (error.size() != n or scale_reference.size() != n)
        throw std::invalid_argument("atstep: norm state sizes differ");
      if (n == 0)
        return 0.0;

      std::size_t fields = 1;
      std::size_t dof_sum = n;
      if (not field_sizes.empty()) {
        fields = field_sizes.size();
        dof_sum = 0;
        for (std::size_t size : field_sizes)
          dof_sum += size;
        if (dof_sum != n)
          throw std::invalid_argument(
            "atstep: field_sizes does not cover the state");
      }

      const double decay = std::exp(-0.2 * static_cast<double>(step_index));
      double total = 0.0;
      std::size_t offset = 0;
      for (std::size_t f = 0; f < fields; ++f) {
        const std::size_t size = field_sizes.empty() ? n : field_sizes[f];
        if (size == 0)
          throw std::invalid_argument("atstep: field_sizes contains an empty field");
        double limit = absolute;
        if (scaled) {
          double norm2 = 0.0;
          double norm_inf = 0.0;
          for (std::size_t i = offset; i < offset + size; ++i) {
            const double magnitude = std::abs(scale_reference[i]);
            if (not std::isfinite(magnitude))
              return std::numeric_limits<double>::infinity();
            norm2 += magnitude * magnitude;
            if (magnitude > norm_inf)
              norm_inf = magnitude;
          }
          norm2 = std::sqrt(norm2);
          if (size == 1 and norm_inf > 0.0 and norm_inf < 1.0)
            limit = absolute * ((2.0 - decay) * norm_inf + decay);
          else {
            const double beta = (norm_inf > 0.0 and norm_inf < 1.0)
              ? (1.0 - decay) * norm_inf + decay
              : 1.0;
            limit = absolute * (beta + norm2);
          }
        }

        double field = 0.0;
        for (std::size_t i = offset; i < offset + size; ++i) {
          const double weight = limit + relative * std::abs(reference[i]);
          if (not std::isfinite(error[i]) or not std::isfinite(weight)
            or weight < 0.0)
            return std::numeric_limits<double>::infinity();
          if (weight == 0.0) {
            if (error[i] != 0.0)
              return std::numeric_limits<double>::infinity();
          }
          else {
            const double ratio = error[i] / weight;
            field += ratio * ratio;
          }
        }
        total += field / static_cast<double>(size);
        offset += size;
      }
      return std::sqrt(total / static_cast<double>(fields));
    }
  };

} // namespace atstep
