// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
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

  /// The automatic initial step, as a fraction of the interval that remains
  /// after the consistent initialization. Measured: 0.1% of what is left, not
  /// of the span and not of the end time (doc/evidence/stepping-sweep.md, Q05
  /// and Q02).
  inline constexpr double default_initial_step_fraction = 1e-3;

  /// The automatic maximum step, as a fraction of the span: COMSOL's Automatic
  /// constraint, measured as exactly 0.1 of the span
  /// (doc/evidence/stepping-sweep.md, Q04).
  inline constexpr double automatic_max_step_fraction = 0.1;

  /// The smallest step a run may be reduced to before it gives up, as a fraction
  /// of the span: the library's own guard, with no COMSOL counterpart.
  inline constexpr double default_min_step_fraction = 1e-10;

  /// Where the solver's own steps land: COMSOL's "Steps taken by solver".
  /// `free` lets the stepping choose, `intermediate` puts at least one step in
  /// every interval between two output times, `strict` ends every step at an
  /// output time.
  enum class Placement { free, intermediate, strict };

  /// What the run reports: COMSOL's "Times to store". `interpolate` stores the
  /// values at the output times, `closest` the first solver step reaching each
  /// of them, `steps` the solver's own steps.
  enum class Storage { interpolate, closest, steps };

  /// The largest step: COMSOL's "Maximum step constraint".
  struct AutomaticMaxStep { };

  struct ConstantMaxStep {
    /// The maximum step, in the run's time unit. COMSOL's "maxstepbdf" is an
    /// absolute time value, not a fraction of the span.
    double seconds = 0.0;
  };

  struct ExpressionMaxStep {
    /// The maximum step as a function of the time, evaluated while solving;
    /// its absolute value is used.
    std::function<double(double)> value {};
  };

  using MaximumStep =
    std::variant<AutomaticMaxStep, ConstantMaxStep, ExpressionMaxStep>;

  /// A step the run chooses: its size from the local error, its order from the
  /// order selection.
  struct AdaptiveSteps {
    /// Where a step may land.
    Placement placement = Placement::free;

    /// Whether a step may pass the last output time, whose value is then
    /// interpolated from the steps around it. COMSOL offers this for every
    /// placement but Strict.
    bool interpolate_end_time = true;

    /// The largest step it may take.
    MaximumStep max_step = AutomaticMaxStep {};

    /// Lowest order the order selection may take: COMSOL's "Minimum BDF order",
    /// which keeps a solver off the first order, so 1 or 2.
    int min_order = 1;

    /// Highest order the order selection may take: 1 to 5.
    int max_order = 5;

    /// COMSOL's "Initial step": unset means automatic, a fraction of what
    /// remains of the interval after the consistent initialization.
    std::optional<double> initial_step {};

    /// The smallest step, as a fraction of the span integrated.
    double min_step_fraction = default_min_step_fraction;
  };

  /// A step of the model's own: the time-step error test is off, so neither the
  /// size nor the order follows from the tolerance.
  struct ManualSteps {
    /// The step, in the time unit of the run.
    double step = 0.0;

    /// The order the step is taken at: COMSOL's "BDF order" of a manual step,
    /// 1 to 5.
    int order = 2;

    /// COMSOL's per-order "Initial step fraction" and "Initial step growth
    /// rate": unset means the default for `order`.
    std::optional<double> startup_fraction {};
    std::optional<double> startup_growth {};

    /// Whether a step may pass the last output time.
    bool interpolate_end_time = true;
  };

  /// One dependent variable: the degrees of freedom it holds, COMSOL's scale
  /// of it, and whether it is algebraic.
  struct Field {
    /// The number of degrees of freedom of the field.
    std::size_t dofs = 0;

    /// COMSOL's "scale for dependent variable", which its solver log prints.
    double scale = 1.0;

    /// Whether the field is algebraic, for COMSOL's "Exclude algebraic".
    bool algebraic = false;
  };

  /// The partition of the state into dependent variables. An empty list means
  /// one field over every degree of freedom.
  struct Layout {
    std::vector<Field> fields;
  };

  /// The absolute tolerance: COMSOL's "Tolerance method". `factor` makes it a
  /// factor of the relative tolerance, `manual` a value of its own.
  struct AbsoluteTolerance {
    enum class Kind { factor, manual };
    Kind kind = Kind::factor;
    /// The factor of the relative tolerance, or the manual value.
    double value = 0.1;
  };

  /// The local error test, and the tolerance it holds a degree of freedom to.
  struct ErrorTest {
    /// The relative tolerance `R`.
    double relative = 1e-2;

    /// The absolute tolerance.
    AbsoluteTolerance absolute {};

    /// Whether the tolerance is applied to scaled variables, COMSOL's "Scaled"
    /// method, or to the variables as they are, its "Unscaled" one.
    bool scaled = true;

    /// COMSOL's "Update scaled absolute tolerance": whether the absolute
    /// tolerance follows the magnitude of the solution.
    bool update_scaled = true;

    /// Whether the algebraic degrees of freedom are held to the tolerance:
    /// COMSOL's "Error estimation".
    bool include_algebraic = true;

    /// The partition of the state into dependent variables.
    Layout layout {};
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

  /// No consistent initialization: the initial values are taken as they are,
  /// and the caller supplies the initial derivative.
  struct AssumeConsistent { };

  /// The time stepping of a run.
  ///
  /// The defaults are COMSOL's, read back from a COMSOL 6.2 Time-Dependent
  /// Solver node (see doc/comsol-alignment.md): the relative tolerance 0.01, the
  /// absolute tolerance `0.1` times it, scaled by the magnitude of each degree
  /// of freedom, the first step `1e-3` of what remains of the interval, the
  /// largest step `0.1` of the span, the orders 1 to 5, and the interpolation
  /// of the end time.
  struct Settings {
    /// How the run takes its steps: the solver's, from the local error, or the
    /// model's own at a fixed step and order.
    std::variant<AdaptiveSteps, ManualSteps> steps {};

    /// The tolerance the local error is held to.
    ErrorTest error {};

    /// The consistent initialization every run starts with.
    std::variant<BackwardEuler, AssumeConsistent> initialization = BackwardEuler {};

    /// Which times the run reports.
    Storage storage = Storage::interpolate;

    /// COMSOL's "Store every Nth step"; the initial state is step 0 and the
    /// modulus is anchored on it (doc/evidence/stepping-sweep.md, Q07).
    int store_every = 1;
  };

  /// Reject a setting a run cannot be built from. Each message names the field.
  void validate(const Settings& settings);

  inline void validate(const Settings& settings)
  {
    const auto require = [](const char* name, bool ok) {
      if (not ok)
        throw std::invalid_argument(std::string("atstep: ") + name + " is out of range");
    };
    const auto finite_positive = [](double value) {
      return std::isfinite(value) and value > 0.0;
    };

    require("error.relative", finite_positive(settings.error.relative));
    require("error.absolute.value", std::isfinite(settings.error.absolute.value));
    require("error.absolute.value",
      settings.error.absolute.kind == AbsoluteTolerance::Kind::factor
        ? settings.error.absolute.value >= 0.0
        : finite_positive(settings.error.absolute.value));
    require("error.relative * error.absolute.value",
      std::isfinite(settings.error.relative * settings.error.absolute.value));
    for (const Field& field : settings.error.layout.fields)
      require("error.layout.fields[].scale", finite_positive(field.scale));

    if (const auto* manual = std::get_if<ManualSteps>(&settings.steps)) {
      require("steps.manual.step", finite_positive(manual->step));
      require("steps.manual.order", manual->order >= 1 and manual->order <= max_order);
    }
    else {
      const auto& adaptive = std::get<AdaptiveSteps>(settings.steps);
      require("steps.adaptive.min_order",
        adaptive.min_order >= 1 and adaptive.min_order <= 2);
      require("steps.adaptive.max_order",
        adaptive.max_order >= adaptive.min_order and adaptive.max_order <= max_order);
      require("steps.adaptive.min_step_fraction", finite_positive(adaptive.min_step_fraction));
      if (adaptive.initial_step)
        require("steps.adaptive.initial_step", finite_positive(*adaptive.initial_step));
      if (const auto* constant = std::get_if<ConstantMaxStep>(&adaptive.max_step))
        require("steps.adaptive.max_step.seconds", finite_positive(constant->seconds));
    }
    if (const auto* backward = std::get_if<BackwardEuler>(&settings.initialization))
      require("initialization.fraction", finite_positive(backward->fraction));
    require("store_every", settings.store_every >= 1);
  }

} // namespace atstep
