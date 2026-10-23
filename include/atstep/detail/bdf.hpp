// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace atstep::detail {

  /// The local error estimates of one step, and the value the step is held to.
  ///
  /// They are IDA's (`IDATestError`), which is the solver COMSOL's BDF is.
  /// With `k` the order of the step and `ee` its correction — the difference
  /// between the corrected state and the value the predictor extrapolated to
  /// it — the estimate at the order `j` is `sigma_j ||ee_j||`, the weighted
  /// norm of the correction a polynomial of order `j` would have left:
  ///
  ///     at     = sigma_k   ||ee_k||      the order the step was taken at
  ///     below  = sigma_k-1 ||ee_k-1||    one order lower, where it exists
  ///     below2 = sigma_k-2 ||ee_k-2||    two orders lower, where it exists
  ///     above  = ||ee_k - ee_k-1|| / (k + 2), the estimate at the order above
  ///
  /// `test` is `ck ||ee_k||`, and the step is accepted while it is at most
  /// one: `ck` differs from `sigma_k` only when the step size changed, and
  /// IDA tests the former. `passed()` is that test.
  struct StepErrors {
    double test = 0.0;
    double at = 0.0;
    std::optional<double> below;
    std::optional<double> below2;
    std::optional<double> above;

    /// Whether the step met the tolerance.
    bool passed() const noexcept
    {
      return std::isfinite(test) and test >= 0.0 and test <= 1.0;
    }
  };

  /// The level a run has accepted: the state, and the time it was accepted at.
  ///
  /// A level pairs its state with its time in one object, so a copy that throws
  /// in the member-wise assignment leaves the time alone. The new level is
  /// constructed before the member is written, which is what leaves the
  /// accepted level alone when the copy throws.
  ///
  /// Nothing behind it is kept: the scaled differences of the accepted steps
  /// are the history the scheme reads, and they are differences of states, not
  /// older levels.
  template <class State> class AcceptedLevel {
  public:
    /// Take `state` as the accepted level, at `time`.
    void commit(double time, const State& state)
    {
      Level level {state, time};
      level_ = level;
    }

    double time() const noexcept
    {
      return level_.time;
    }

    const State& state() const noexcept
    {
      return level_.state;
    }

  private:
    /// The state comes first so that a copy that throws in the member-wise
    /// assignment leaves the time alone.
    struct Level {
      State state {};
      double time = 0.0;
    };

    Level level_ {};
  };

} // namespace atstep::detail
