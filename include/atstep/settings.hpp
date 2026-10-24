// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace atstep {

  /// Highest BDF order the scheme implements.
  inline constexpr int max_bdf_order = 5;

  namespace detail {
    /// The scaled-difference table holds the columns `0..max_order + 1`.
    inline constexpr int max_levels = max_bdf_order + 2;
  } // namespace detail

  /// COMSOL's "Steps taken by solver".
  enum class StepMode {
    /// The solver chooses the step. An output time a step stepped over is
    /// interpolated from the dense output.
    free,
    /// Every subinterval between two output times holds at least one step.
    intermediate,
    /// Every output time is a step end point.
    strict,
    /// A step of the model's own, with the time-step error test switched off.
    manual,
  };

  /// COMSOL's "Times to store".
  enum class StoreMode {
    /// The output times, interpolated from the steps that stepped over them.
    interpolate,
    /// The first solver step reaching each output time.
    closest,
    /// The solver's own steps.
    steps,
  };

  /// The time stepping of a run.
  struct Settings {
    /// Lowest order the order selection may take. The first step of a run is
    /// always taken at order one: it has one level behind it.
    int min_order = 1;
    /// Highest order the order selection may take: 1 to `max_bdf_order`.
    int max_order = max_bdf_order;
    /// The order of a `StepMode::manual` step.
    int manual_order = 2;

    StepMode step_mode = StepMode::free;
    StoreMode store_mode = StoreMode::interpolate;

    /// Whether a step may pass the final output time, whose value is then
    /// interpolated from the steps around it.
    bool interpolate_end_time = true;

    /// The relative tolerance `R`: a dof is held to the local error
    /// `A + R abs(y_i)`, with `A` the absolute tolerance below.
    double relative_tolerance = 1e-2;

    /// The absolute tolerance as a factor of the relative one:
    /// `A = R * absolute_tolerance_factor`.
    double absolute_tolerance_factor = 0.1;

    /// Whether the absolute tolerance is scaled by the magnitude of each
    /// field: COMSOL's "Scaled" method.
    bool scale_absolute_tolerance = true;

    /// The number of degrees of freedom in each field of the state, in state
    /// order. Empty means one field holding the whole state. The error norm
    /// weights fields equally and dofs within a field by their count.
    std::vector<std::size_t> field_sizes;

    /// The automatic first step, as a fraction of the final time (its
    /// magnitude, and at least the span of the run). In a manual run the base
    /// is `manual_step`, not the final time.
    double initial_step_fraction = 1e-3;

    /// The largest step, in the time unit of the run. Zero means no limit.
    double max_step = 0.0;

    /// The smallest step the controller proposes, as a fraction of the span.
    double min_step_fraction = 1e-10;

    /// The step of `StepMode::manual`, in the time unit of the run.
    double manual_step = 0.0;

    /// The artificial backward-Euler steps of consistent initialization, as a
    /// fraction of the automatic initial step. COMSOL's default is 1e-3.
    double backward_euler_fraction = 1e-3;
  };

  /// Reject a setting a run cannot be built from. The message names the field.
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
    const auto finite = [](const char* field, double value) {
      if (not std::isfinite(value))
        detail::reject_setting(field, "not finite", value);
    };
    finite("relative_tolerance", settings.relative_tolerance);
    finite("absolute_tolerance_factor", settings.absolute_tolerance_factor);
    finite("initial_step_fraction", settings.initial_step_fraction);
    finite("max_step", settings.max_step);
    finite("min_step_fraction", settings.min_step_fraction);
    finite("backward_euler_fraction", settings.backward_euler_fraction);
    finite("manual_step", settings.manual_step);
    finite(
      "absolute_tolerance",
      settings.relative_tolerance * settings.absolute_tolerance_factor);
    for (std::size_t size : settings.field_sizes)
      if (size == 0)
        throw std::invalid_argument(
          "atstep: field_sizes contains an empty field");
    if (settings.max_order < 1 or settings.max_order > max_bdf_order)
      detail::reject_setting("max_order",
        "outside 1..max_bdf_order",
        static_cast<double>(settings.max_order));
    if (settings.min_order < 1 or settings.min_order > settings.max_order)
      detail::reject_setting("min_order",
        "outside 1..max_order",
        static_cast<double>(settings.min_order));
    if (settings.manual_order < 1 or settings.manual_order > max_bdf_order)
      detail::reject_setting("manual_order",
        "outside 1..max_bdf_order",
        static_cast<double>(settings.manual_order));
    if (not(settings.relative_tolerance > 0.0))
      detail::reject_setting("relative_tolerance",
        "not positive",
        settings.relative_tolerance);
    if (not(settings.absolute_tolerance_factor >= 0.0))
      detail::reject_setting("absolute_tolerance_factor",
        "negative",
        settings.absolute_tolerance_factor);
    if (not(settings.initial_step_fraction > 0.0))
      detail::reject_setting("initial_step_fraction",
        "not positive",
        settings.initial_step_fraction);
    if (not(settings.max_step >= 0.0))
      detail::reject_setting("max_step", "negative", settings.max_step);
    if (not(settings.min_step_fraction > 0.0))
      detail::reject_setting("min_step_fraction",
        "not positive",
        settings.min_step_fraction);
    if (not(settings.backward_euler_fraction > 0.0))
      detail::reject_setting("backward_euler_fraction",
        "not positive",
        settings.backward_euler_fraction);
    if (settings.step_mode == StepMode::manual and not(settings.manual_step > 0.0))
      detail::reject_setting("manual_step",
        "not a positive step of a manual run",
        settings.manual_step);
  }

} // namespace atstep
