// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

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

  /// The step a run has accepted: the level, the time it was accepted at, and
  /// the time derivative the solve left there.
  ///
  /// The three are one object because `accept` publishes them together: a run
  /// either has the whole step or none of it, and no caller can see a level
  /// that is not paired with its own derivative. The new step is constructed
  /// before the accepted one is written, so a copy that throws while the step
  /// is being published leaves the accepted step alone.
  ///
  /// Nothing behind it is kept: the scaled differences of the accepted steps
  /// are the history the scheme reads, and they are differences of states, not
  /// older levels.
  template <class State> class AcceptedStep {
  public:
    /// Publish `state` and its derivative `derivative` as the accepted step at
    /// `time`.
    void commit(double time, const State& state, const State& derivative)
    {
      Step step {state, derivative, time};
      step_ = std::move(step);
    }

    double time() const noexcept
    {
      return step_.time;
    }

    const State& state() const noexcept
    {
      return step_.state;
    }

    /// The derivative the solve of the step evaluated at the level, which the
    /// dense output returns at the newest time.
    const State& derivative() const noexcept
    {
      return step_.derivative;
    }

  private:
    /// The state and its derivative come first: a copy that throws in the
    /// member-wise assignment then leaves the time alone.
    struct Step {
      State state {};
      State derivative {};
      double time = 0.0;
    };

    Step step_ {};
  };

} // namespace atstep::detail
