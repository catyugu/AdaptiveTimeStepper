// The algebraic solve of COMSOL's recorded decay models, emulated for the
// tests.
//
// COMSOL's default Fully Coupled node reuses its Jacobian instead of solving
// the system exactly, and that choice is visible in the recorded
// trajectories: the library's error estimate is the Newton correction
// (TimeStepper::estimate), so a different solve gives a different step
// sequence.  The tests that compare step sequences therefore have to reproduce
// the solve, not just the equation.
//
// This adapter is specific to the recorded models and is not part of the
// library.  It reaches into two accessors the rewrite removes
// (`request.previous()` and `request.absolute_tolerance_reference()`), so it
// changes together with the request boundary (doc/refactor-plan.md, F08).
#pragma once

#include "atstep/stepper.hpp"

#include <cmath>
#include <cstddef>

namespace comsol_decay {

  /// The state of the recorded model: one degree of freedom, `u`.
  struct State {
    double u = 0.0;
    void axpy(double a, const State& x)
    {
      u += a * x.u;
    }
    std::size_t size() const
    {
      return 1;
    }
    double operator[](std::size_t) const
    {
      return u;
    }
  };

  /// The solve of `u' + u = 0` with COMSOL's reused Jacobian.
  struct Solve {
    atstep::Settings settings;
    double jacobian_cj = 0.0;
    double last_cj = 0.0;
    double convergence_factor = 20.0;
    /// The derivative of the residual's evaluation: the only state the solve
    /// writes besides `y`.
    State derivative;

    atstep::SolveOutcome
    operator()(const atstep::StepRequest<State>& request, State& y)
    {
      const double cj = request.derivative_weight();
      if (request.initializing) {
        y.u = -request.derivative_history().u / (cj + 1.0);
        jacobian_cj = 0.0;
        return {true, 1};
      }
      if (cj != last_cj)
        convergence_factor = 100.0;
      if (jacobian_cj == 0.0 or cj / jacobian_cj < 0.6 or cj / jacobian_cj > 1.0 / 0.6) {
        jacobian_cj = cj;
        convergence_factor = 20.0;
      }
      atstep::detail::WeightedRmsNorm norm;
      norm.relative = settings.error.relative;
      norm.absolute = settings.error.relative * settings.error.absolute.value;
      norm.scaled = settings.error.scaled;
      norm.step_index = request.step_index;
      double first_norm = 0.0;
      last_cj = cj;
      for (int iteration = 0; iteration < 4; ++iteration) {
        request.derivative(y, derivative);
        const double residual = y.u + derivative.u;
        const double correction = -residual / (jacobian_cj + 1.0) * 2.0 / (1.0 + cj / jacobian_cj);
        y.u += correction;
        const double size = norm(State {correction},
          request.previous(),
          request.absolute_tolerance_reference());
        if (iteration == 0)
          first_norm = size;
        else {
          const double rate = std::pow(size / first_norm, 1.0 / static_cast<double>(iteration));
          if (rate > 0.9)
            return {false, iteration + 1};
          convergence_factor = rate / (1.0 - rate);
        }
        if (size <= 3.3e-9 or convergence_factor * size <= 0.33)
          return {true, iteration + 1};
      }
      return {false, 4};
    }
  };

} // namespace comsol_decay
