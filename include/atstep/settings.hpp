// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace atstep {

  /// The highest BDF order, the highest derivative level, and the number of
  /// time levels a formula reads. The orders are COMSOL's: its BDF runs at 1
  /// to 5.
  inline constexpr int max_order = 5;
  inline constexpr int max_derivative = 3;
  inline constexpr int max_levels = max_order + 2;

  /// How the solver places its steps relative to the output times of the
  /// study. This is COMSOL's "Steps taken by solver" setting, with the same
  /// four values.
  enum class StepMode {
    /// The step is the solver's own. An output time a step stepped over is
    /// reported by interpolating the scheme's polynomial there.
    free,
    /// The step is the solver's own, but every subinterval between two output
    /// times holds at least one step.
    intermediate,
    /// Land on every output time, taking smaller internal steps as needed.
    strict,
    /// A step of the model's own, with the time-step error test switched off.
    manual,
  };

  /// Which times the result is stored at: COMSOL's "Times to store".
  enum class StoreMode {
    /// The output times, interpolated from the steps that stepped over them.
    interpolate,
    /// The first solver step reaching each output time.
    closest,
    /// The solver's own steps.
    steps,
  };

  /// The time stepping of a run.
  ///
  /// The defaults are COMSOL's, read back from a COMSOL 6.2 Time-Dependent
  /// Solver node (see doc/comsol-alignment.md): the relative tolerance 0.01,
  /// the absolute tolerance `0.1` times it, the first step `1e-3` of the
  /// span, the largest step `0.1` of the span, the orders 1 to 5, and the
  /// interpolation of the end time.
  struct Settings {
    /// Lowest order the order selection may take. COMSOL's solver accepts 1
    /// and 2 here (its "Minimum BDF order" keeps a solver off the first
    /// order).
    int min_order = 1;

    /// Highest order the order selection may take: 1 to 5.
    int max_order = 5;

    /// The order a step of `StepMode::manual` is taken at: COMSOL's "BDF
    /// order" of a manual step, 2 by default.
    int manual_order = 2;

    StepMode step_mode = StepMode::free;
    StoreMode store_mode = StoreMode::interpolate;

    /// Whether the stepping may pass the last output time, whose value is
    /// then interpolated from the steps around it.
    bool interpolate_end_time = true;

    /// The relative tolerance `R`: a dof is held to the local error
    /// `A + R |y_i|`, with `A` the absolute tolerance below. It bounds the
    /// error of a step, so it is the accuracy of the time discretization.
    double relative_tolerance = 1e-2;

    /// The factor of the absolute tolerance `A = R * absolute_factor`: what
    /// keeps a dof whose value is near zero from being held to a relative
    /// accuracy on a value that carries nothing. COMSOL's "Factor" method
    /// of the absolute tolerance, 0.1 by default.
    double absolute_factor = 0.1;

    /// Whether the absolute tolerance is scaled by the magnitude of each
    /// degree of freedom, COMSOL's "Scaled" method:
    ///
    ///     A_i = R * absolute_factor * (|y_i| + s_i),
    ///     s_i = (1 - e^(-j/5)) |y_i| + e^(-j/5) sigma_i,
    ///
    /// with `j` the number of accepted steps and `sigma_i` the scale below. It
    /// is COMSOL's own default; off, `A` is the fixed `R * absolute_factor`,
    /// which is what a run of a physical field whose magnitude the mesh changes
    /// is easier to reason about with.
    bool scale_absolute_tolerance = true;

    /// The scale `sigma_i` of each degree of freedom in the scaled absolute
    /// tolerance above. It is COMSOL's "scale for dependent variable", which
    /// its solver log prints: one for a variable whose equation moves at a
    /// moderate rate, and below one for a fast one — a variable obeying
    /// `v' = -100 v` carries 0.32. COMSOL derives it from the model, and the
    /// library cannot, so a run of a model COMSOL would scale differently
    /// passes it here. An empty vector is `sigma_i = 1`, which is the value of
    /// every configuration of the recorded reference.
    std::vector<double> absolute_scale;

    /// The first step, as a fraction of the span integrated. COMSOL's
    /// "Initial step: Automatic", 0.1% of the end time. The initial
    /// derivative of the state may take the first step below it (see
    /// `StepController::first_step()`).
    double first_step_fraction = 1e-3;

    /// The largest step, as a fraction of the span integrated: COMSOL's
    /// "Maximum step constraint: Automatic", 10% of the span. A finite
    /// nonpositive value disables the maximum-step constraint.
    double max_step_fraction = 0.1;

    /// The smallest step, as a fraction of the span integrated: below it
    /// the estimated error cannot be met by shrinking the step any further.
    double min_step_fraction = 1e-10;

    /// The step of `StepMode::manual`, in the time unit of the run.
    double manual_step = 0.0;

    /// The artificial backward-Euler steps of consistent initialization, as a
    /// fraction of the first step: COMSOL's "Fraction of initial step for
    /// Backward Euler", 0.001. The run takes two of them; the increment of the
    /// second is the time derivative the initial-step rule reads.
    double init_step_fraction = 1e-3;

    /// How many derivative levels the step formula reports (see
    /// `StepFormula::weights`): 1 for a first-order-in-time equation, 2 for one
    /// with a second time derivative (`M u'' + C u' + K u = f`), up to 3.
    ///
    /// It also selects the formulation, not only how many weights the solve is
    /// given: one uses the fixed-leading BDF history, above one the classical
    /// stencil of `StepFormula`, whose convergence on a directly discretized
    /// second derivative is first order (see the scope note in the README).
    int derivative_orders = 1;
  };

  /// Reject a setting a run cannot be built from. The message names the
  /// field.
  void validate(const Settings& settings);

  namespace detail {
    [[noreturn]] inline void
    reject_setting(const char* field, const char* why, double value)
    {
      throw std::invalid_argument(std::string("atstep: ") + field + " is " + why + " (" + std::to_string(value) + ")");
    }
  } // namespace detail

  inline void validate(const Settings& settings)
  {
    // Check finiteness once, before range checks and step-bound arithmetic.
    const auto finite = [](const char* field, double value) {
      if (not std::isfinite(value))
        detail::reject_setting(field, "not finite", value);
    };
    finite("relative_tolerance", settings.relative_tolerance);
    finite("absolute_factor", settings.absolute_factor);
    finite("first_step_fraction", settings.first_step_fraction);
    finite("max_step_fraction", settings.max_step_fraction);
    finite("min_step_fraction", settings.min_step_fraction);
    finite("init_step_fraction", settings.init_step_fraction);
    finite("manual_step", settings.manual_step);
    finite("absolute_tolerance", settings.relative_tolerance * settings.absolute_factor);
    for (double scale : settings.absolute_scale)
      finite("absolute_scale", scale);
    if (settings.max_step_fraction > 0.0
      and settings.min_step_fraction > settings.max_step_fraction)
      detail::reject_setting("min_step_fraction", "greater than max_step_fraction", settings.min_step_fraction);
    if (settings.max_order < 1 or settings.max_order > max_order)
      detail::reject_setting("max_order",
        "outside 1..5",
        static_cast<double>(settings.max_order));
    if (settings.min_order < 1 or settings.min_order > settings.max_order)
      detail::reject_setting("min_order",
        "outside 1..max_order",
        static_cast<double>(settings.min_order));
    if (settings.manual_order < 1 or settings.manual_order > max_order)
      detail::reject_setting("manual_order",
        "outside 1..5",
        static_cast<double>(settings.manual_order));
    if (not(settings.relative_tolerance > 0.0))
      detail::reject_setting("relative_tolerance",
        "not positive",
        settings.relative_tolerance);
    if (not(settings.absolute_factor >= 0.0))
      detail::reject_setting(
        "absolute_factor", "negative", settings.absolute_factor);
    for (double scale : settings.absolute_scale)
      if (not(scale >= 0.0))
        detail::reject_setting(
          "absolute_scale", "negative", scale);
    if (not(settings.first_step_fraction > 0.0))
      detail::reject_setting("first_step_fraction",
        "not positive",
        settings.first_step_fraction);
    if (not(settings.min_step_fraction > 0.0))
      detail::reject_setting("min_step_fraction",
        "not positive",
        settings.min_step_fraction);
    if (not(settings.init_step_fraction > 0.0))
      detail::reject_setting("init_step_fraction",
        "not positive",
        settings.init_step_fraction);
    if (settings.derivative_orders < 1 or settings.derivative_orders > max_derivative)
      detail::reject_setting(
        "derivative_orders",
        "outside 1..3",
        static_cast<double>(settings.derivative_orders));
    if (settings.step_mode == StepMode::manual and not(settings.manual_step > 0.0))
      detail::reject_setting("manual_step",
        "not a positive step of a manual run",
        settings.manual_step);
  }

} // namespace atstep
