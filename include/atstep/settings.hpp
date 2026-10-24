// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace atstep {

  /// The highest BDF order, and the number of time levels the scaled-difference
  /// table holds: the scheme reads the columns `0..max_order + 1`. The orders
  /// are COMSOL's: its BDF runs at 1 to 5.
  inline constexpr int max_order = 5;
  inline constexpr int max_levels = max_order + 2;

  /// Where a step the solver chose lands relative to the output times of the
  /// study: COMSOL's "Steps taken by solver", less the manual step, which is a
  /// step of the model's own (see `ManualStep`).
  enum class StepPlacement {
    /// The step is the solver's own. An output time a step stepped over is
    /// reported by interpolating the scheme's polynomial there.
    free,
    /// The step is the solver's own, but every subinterval between two output
    /// times holds at least one step.
    intermediate,
    /// Land on every output time, taking smaller internal steps as needed.
    strict,
  };

  /// A step the solver chooses: its size from the local error, its order from
  /// the order selection, and `placement` where it may land.
  struct AutomaticStep {
    StepPlacement placement = StepPlacement::free;
  };

  /// A step of the model's own. The time-step error test is off, so neither
  /// the size nor the order follows from the tolerance.
  struct ManualStep {
    /// The step, in the time unit of the run: positive, which `validate`
    /// enforces.
    double step = 0.0;

    /// The order the step is taken at: COMSOL's "BDF order" of a manual step,
    /// 1 to 5.
    int order = 2;
  };

  /// How a run takes its steps. The two cases share no parameter: a run is
  /// either the solver's or the model's own, and neither holds a field the
  /// other reads.
  using StepControl = std::variant<AutomaticStep, ManualStep>;

  /// Which times the result is stored at: COMSOL's "Times to store".
  enum class StoreMode {
    /// The output times, interpolated from the steps that stepped over them.
    interpolate,
    /// The first solver step reaching each output time.
    closest,
    /// The solver's own steps.
    steps,
  };

  /// The same absolute tolerance on every degree of freedom, `A = R *
  /// absolute_factor`. It is COMSOL's "Unscaled" method, and the one a run of a
  /// physical field whose magnitude the mesh changes is easier to reason about
  /// with.
  struct FixedAbsoluteTolerance { };

  /// The absolute tolerance `A_i` of each degree of freedom, scaled by the
  /// magnitude of the solution, COMSOL's "Scaled" method:
  ///
  ///     A_i = R * absolute_factor * (|y_i| + s_i),
  ///     s_i = (1 - e^(-j/5)) |y_i| + e^(-j/5) sigma_i,
  ///
  /// with `j` the number of accepted steps and `sigma_i` the scale below. It
  /// keeps a degree of freedom whose value is near zero from being held to a
  /// relative accuracy on a value that carries nothing, and is COMSOL's own
  /// default.
  struct ScaledAbsoluteTolerance {
    /// The scale `sigma_i` of each degree of freedom. It is COMSOL's "scale for
    /// dependent variable", which its solver log prints: one for a variable
    /// whose equation moves at a moderate rate, and below one for a fast one — a
    /// variable obeying `v' = -100 v` carries 0.32. COMSOL derives it from the
    /// model, and the library cannot, so a run of a model COMSOL would scale
    /// differently passes it here. An empty vector is `sigma_i = 1`, and so is
    /// an index past the end of it; one is the value of every configuration of
    /// the recorded reference.
    std::vector<double> scale;
  };

  /// How the absolute tolerance of the error test is shaped: one value for
  /// every degree of freedom, or one per degree of freedom that the solution's
  /// own magnitude scales.
  using AbsoluteTolerance = std::variant<FixedAbsoluteTolerance, ScaledAbsoluteTolerance>;

  /// The time stepping of a run.
  ///
  /// The defaults are COMSOL's, read back from a COMSOL 6.2 Time-Dependent
  /// Solver node (see doc/comsol-alignment.md): the relative tolerance 0.01,
  /// the absolute tolerance `0.1` times it, scaled by the magnitude of each
  /// degree of freedom, the first step `1e-3` of the span, the largest step
  /// `0.1` of the span, the orders 1 to 5, and the interpolation of the end
  /// time.
  struct Settings {
    /// Lowest order the order selection may take. COMSOL's solver accepts 1
    /// and 2 here (its "Minimum BDF order" keeps a solver off the first
    /// order).
    int min_order = 1;

    /// Highest order the order selection may take: 1 to 5.
    int max_order = 5;

    /// How the run takes its steps: the solver's, from the local error, or the
    /// model's own at a fixed step and order.
    StepControl stepping = AutomaticStep {};

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

    /// The shape of the absolute tolerance `A`: the fixed `R * absolute_factor`
    /// on every degree of freedom, or the scaled form of
    /// `ScaledAbsoluteTolerance`.
    AbsoluteTolerance absolute_tolerance = ScaledAbsoluteTolerance {};

    /// The first step, as a fraction of the span integrated. COMSOL's
    /// "Initial step: Automatic" is 0.1% of the end time for BDF, subject to
    /// DT * yp_norm <= 1/2; the span and the end time differ when the run
    /// starts after zero, and which one COMSOL uses there is unverified. The
    /// initial derivative may take the first step below it (see
    /// `StepController::first_step()`).
    double first_step_fraction = 1e-3;

    /// The largest step, as a fraction of the span integrated. COMSOL's
    /// "Maximum step constraint: Automatic" has no documented value; this
    /// 10% of the span is the library's own choice, and no recorded run makes
    /// it bind. A finite nonpositive value disables the constraint.
    double max_step_fraction = 0.1;

    /// The smallest step, as a fraction of the span integrated: below it
    /// the estimated error cannot be met by shrinking the step any further.
    double min_step_fraction = 1e-10;

    /// The artificial backward-Euler steps of consistent initialization, as a
    /// fraction of the first step: COMSOL's "Fraction of initial step for
    /// Backward Euler", 0.001. The run takes two of them; the increment of the
    /// second is the time derivative the initial-step rule reads.
    double init_step_fraction = 1e-3;

    /// The step of a manual run, or nothing when the solver chooses the steps.
    std::optional<ManualStep> manual() const
    {
      if (const auto* step = std::get_if<ManualStep>(&stepping))
        return *step;
      return std::nullopt;
    }

    /// Where the steps the solver chooses may land. A manual step is placed as
    /// a free one: the step of the model's own is not cut to an output time.
    StepPlacement placement() const
    {
      if (const auto* step = std::get_if<AutomaticStep>(&stepping))
        return step->placement;
      return StepPlacement::free;
    }
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
    const auto negative = [](const char* field, double value) {
      if (not(value >= 0.0))
        detail::reject_setting(field, "negative", value);
    };
    // The settings a run reads are the ones of the alternative it is in.
    const auto* scaled = std::get_if<ScaledAbsoluteTolerance>(&settings.absolute_tolerance);
    const auto* manual = std::get_if<ManualStep>(&settings.stepping);
    finite("relative_tolerance", settings.relative_tolerance);
    finite("absolute_factor", settings.absolute_factor);
    finite("first_step_fraction", settings.first_step_fraction);
    finite("max_step_fraction", settings.max_step_fraction);
    finite("min_step_fraction", settings.min_step_fraction);
    finite("init_step_fraction", settings.init_step_fraction);
    finite("relative_tolerance * absolute_factor",
      settings.relative_tolerance * settings.absolute_factor);
    if (scaled)
      for (double scale : scaled->scale)
        finite("absolute_scale", scale);
    if (manual)
      finite("manual_step", manual->step);
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
    if (manual and (manual->order < 1 or manual->order > max_order))
      detail::reject_setting("manual_order",
        "outside 1..5",
        static_cast<double>(manual->order));
    if (not(settings.relative_tolerance > 0.0))
      detail::reject_setting("relative_tolerance",
        "not positive",
        settings.relative_tolerance);
    negative("absolute_factor", settings.absolute_factor);
    if (scaled)
      for (double scale : scaled->scale)
        negative("absolute_scale", scale);
    if (manual and not(manual->step > 0.0))
      detail::reject_setting("manual_step",
        "not a positive step of a manual run",
        manual->step);
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
  }

} // namespace atstep
