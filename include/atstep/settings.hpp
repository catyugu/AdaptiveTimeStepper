// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
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

  /// The automatic initial step of a run, as a fraction of the span integrated:
  /// COMSOL's "Initial step fraction", measured as 0.1% of the span and not of
  /// the end time (doc/evidence/stepping-sweep.md, Q05). A manual run has no such
  /// field, and its two backward-Euler steps advance the clock by
  /// `manual_step * this * initialization.fraction`, 2e-7 s at a manual step of
  /// 0.1 in every measured manual case (Q03).
  inline constexpr double default_initial_step_fraction = 1e-3;

  /// The largest step of a run whose maximum step constraint is automatic, as a
  /// fraction of the span: COMSOL's Automatic constraint, measured as exactly 0.1
  /// of the span (doc/evidence/stepping-sweep.md, Q04) — 10.0 s over six
  /// consecutive steps at a span of 100, and 1.0 s rather than 1.2 s for
  /// t0 = 2, end = 12.
  inline constexpr double automatic_max_step_fraction = 0.1;

  /// The smallest step a run may be reduced to before it gives up, as a fraction
  /// of the span: the library's own guard, with no COMSOL counterpart.
  inline constexpr double default_min_step_fraction = 1e-10;

  /// The solver's own steps. An output time a step stepped over is reported by
  /// interpolating the scheme's polynomial there.
  struct Free {
    /// Whether a step may pass the last output time, whose value is then
    /// interpolated from the steps around it.
    bool interpolate_end_time = true;
  };

  /// The solver's own steps, with at least one in every subinterval between two
  /// output times.
  struct Intermediate {
    bool interpolate_end_time = true;
  };

  /// Land on every output time, taking smaller internal steps as needed. COMSOL
  /// offers no end interpolation here, and none is needed: a step of this
  /// placement ends at the end time.
  struct Strict { };

  /// Where a step the solver chose lands relative to the output times of the
  /// study: COMSOL's "Steps taken by solver", less the manual step, which is a
  /// step of the model's own (see `ManualSteps`).
  using Placement = std::variant<Free, Intermediate, Strict>;

  /// The largest step is a fraction of the span the library chooses.
  struct AutomaticMaxStep { };

  /// The largest step is `fraction` of the span.
  struct ConstantMaxStep {
    double fraction = automatic_max_step_fraction;
  };

  /// No largest step beyond the arithmetic of the run.
  struct UnlimitedMaxStep { };

  /// The largest step of a run: COMSOL's "Maximum step constraint", less the
  /// expression, which the model layer evaluates.
  using MaximumStep = std::variant<AutomaticMaxStep, ConstantMaxStep, UnlimitedMaxStep>;

  /// A step the run chooses: its size from the local error, its order from the
  /// order selection.
  struct AdaptiveSteps {
    /// Where a step may land.
    Placement placement {};

    /// The largest step it may take.
    MaximumStep max_step {};

    /// Lowest order the order selection may take: COMSOL's "Minimum BDF order",
    /// which keeps a solver off the first order, so 1 or 2.
    int min_order = 1;

    /// Highest order the order selection may take: 1 to 5.
    int max_order = 5;

    /// The first step, as a fraction of the span integrated. COMSOL's
    /// "Initial step: Automatic" is 0.1% of the span, not of the end time: at
    /// t0 = 2, end = 12 the initialization advances the clock by 2.0e-5 s rather
    /// than 2.4e-5, and with the error test out of the way the first step is
    /// 0.00999998 s at a span of 10 and 0.0999998 s at a span of 100
    /// (doc/evidence/stepping-sweep.md, Q05). The initial derivative may take the
    /// first step below it, COMSOL's `DT * yp_norm <= 1/2` (see
    /// `StepController::first_step()`).
    double initial_step_fraction = default_initial_step_fraction;

    /// The smallest step, as a fraction of the span integrated: below it the
    /// estimated error cannot be met by shrinking the step any further.
    double min_step_fraction = default_min_step_fraction;
  };

  /// A step of the model's own: the time-step error test is off, so neither the
  /// size nor the order follows from the tolerance.
  struct ManualSteps {
    /// The step, in the time unit of the run: positive, which `validate`
    /// enforces.
    double step = 0.0;

    /// The order the step is taken at: COMSOL's "BDF order" of a manual step,
    /// 1 to 5.
    int order = 2;

    /// Whether the run may pass the last output time, whose value is then
    /// interpolated from the steps around it.
    bool interpolate_end_time = true;
  };

  /// The local error test, and the tolerance it holds a degree of freedom to.
  struct ErrorTest {
    /// The relative tolerance `R`: a dof is held to the local error
    /// `A + R |y_i|`, with `A` the absolute tolerance below. It bounds the error
    /// of a step, so it is the accuracy of the time discretization.
    double relative = 1e-2;

    /// Whether the absolute tolerance is scaled by the magnitude of each degree
    /// of freedom, COMSOL's "Scaled" method, or the same for every dof, its
    /// "Unscaled" one.
    bool scaled = true;

    /// The factor of the absolute tolerance `A = R * absolute_factor`: what keeps
    /// a dof whose value is near zero from being held to a relative accuracy on
    /// a value that carries nothing. COMSOL's "Factor" method of the absolute
    /// tolerance, 0.1 by default.
    double absolute_factor = 0.1;

    /// The scale `sigma_i` of each degree of freedom, which the scaled form
    /// blends the magnitude of the solution with:
    ///
    ///     A_i = R * absolute_factor * ((2 - e^(-j/5)) |returned_i|
    ///                                  + e^(-j/5) sigma_i),
    ///
    /// with `j` the number of accepted steps. It is COMSOL's "scale for
    /// dependent variable", which its solver log prints: one for a variable
    /// whose equation moves at a moderate rate, and below one for a fast one — a
    /// variable obeying `v' = -100 v` carries 0.32. COMSOL derives it from the
    /// model, and the library cannot, so a run of a model COMSOL would scale
    /// differently passes it here. An empty vector is `sigma_i = 1`, and so is
    /// an index past the end of it; one is the value of every configuration of
    /// the recorded reference.
    std::vector<double> scale;
  };

  /// The artificial backward-Euler steps of consistent initialization, which
  /// reconcile the initial values with the constraints of the equations and give
  /// the time derivative the first step's size is read from. The run takes two
  /// of them; the increment of the second is that derivative.
  struct BackwardEuler {
    /// The size of one, as a fraction of the first step: COMSOL's "Fraction of
    /// initial step for Backward Euler".
    double fraction = 1e-3;
  };

  /// The output times, interpolated from the steps that stepped over them: the
  /// values at the requested times.
  struct Interpolate { };

  /// The first solver step reaching each output time.
  struct Closest { };

  /// The solver's own steps: the initial state and every accepted one.
  struct Steps { };

  /// Which times the result is stored at: COMSOL's "Times to store".
  using Storage = std::variant<Interpolate, Closest, Steps>;

  /// The time stepping of a run.
  ///
  /// The defaults are COMSOL's, read back from a COMSOL 6.2 Time-Dependent
  /// Solver node (see doc/comsol-alignment.md): the relative tolerance 0.01, the
  /// absolute tolerance `0.1` times it, scaled by the magnitude of each degree of
  /// freedom, the first step `1e-3` of the span, the largest step `0.1` of the
  /// span, the orders 1 to 5, and the interpolation of the end time.
  struct Settings {
    /// How the run takes its steps: the solver's, from the local error, or the
    /// model's own at a fixed step and order.
    std::variant<AdaptiveSteps, ManualSteps> steps {};

    /// The tolerance the local error is held to.
    ErrorTest error {};

    /// The consistent initialization every run starts with.
    BackwardEuler initialization {};

    /// Which times the run reports.
    Storage storage = Interpolate {};
  };

  /// Reject a setting a run cannot be built from. Each message names the field.
  void validate(const Settings& settings);

  inline void validate(const Settings& settings)
  {
    const auto reject = [](const char* field, const char* why, double value) {
      throw std::invalid_argument(
        std::string("atstep: ") + field + " is " + why + " (" + std::to_string(value) + ")");
    };
    const auto finite = [&reject](const char* field, double value) {
      if (not std::isfinite(value))
        reject(field, "not finite", value);
    };
    const auto positive = [&reject](const char* field, double value) {
      if (not(value > 0.0))
        reject(field, "not positive", value);
    };
    const auto non_negative = [&reject](const char* field, double value) {
      if (not(value >= 0.0))
        reject(field, "negative", value);
    };

    finite("error.relative", settings.error.relative);
    positive("error.relative", settings.error.relative);
    finite("error.absolute_factor", settings.error.absolute_factor);
    non_negative("error.absolute_factor", settings.error.absolute_factor);
    finite("error.relative * error.absolute_factor",
      settings.error.relative * settings.error.absolute_factor);
    for (double scale : settings.error.scale) {
      finite("error.scale", scale);
      non_negative("error.scale", scale);
    }
    finite("initialization.fraction", settings.initialization.fraction);
    positive("initialization.fraction", settings.initialization.fraction);

    if (const auto* adaptive = std::get_if<AdaptiveSteps>(&settings.steps)) {
      if (adaptive->min_order < 1 or adaptive->min_order > 2)
        reject("steps.adaptive.min_order",
          "outside 1..2",
          static_cast<double>(adaptive->min_order));
      if (adaptive->max_order < adaptive->min_order or adaptive->max_order > max_order)
        reject("steps.adaptive.max_order",
          "outside steps.adaptive.min_order..5",
          static_cast<double>(adaptive->max_order));
      finite("steps.adaptive.initial_step_fraction",
        adaptive->initial_step_fraction);
      positive("steps.adaptive.initial_step_fraction",
        adaptive->initial_step_fraction);
      finite("steps.adaptive.min_step_fraction",
        adaptive->min_step_fraction);
      positive("steps.adaptive.min_step_fraction",
        adaptive->min_step_fraction);
      if (const auto* constant = std::get_if<ConstantMaxStep>(&adaptive->max_step)) {
        finite("steps.adaptive.max_step.fraction", constant->fraction);
        positive("steps.adaptive.max_step.fraction", constant->fraction);
        if (adaptive->min_step_fraction > constant->fraction)
          reject("steps.adaptive.min_step_fraction",
            "greater than steps.adaptive.max_step.fraction",
            adaptive->min_step_fraction);
      }
    }
    else {
      const auto& manual = std::get<ManualSteps>(settings.steps);
      finite("steps.manual.step", manual.step);
      positive("steps.manual.step", manual.step);
      if (manual.order < 1 or manual.order > max_order)
        reject("steps.manual.order",
          "outside 1..5",
          static_cast<double>(manual.order));
    }
  }

} // namespace atstep
