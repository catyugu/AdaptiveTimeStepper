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
  /// after the consistent initialization: 0.1% of what is left, not of the span
  /// and not of the end time (doc/evidence/stepping-sweep.md).
  inline constexpr double default_initial_step_fraction = 1e-3;

  /// The automatic maximum step, as a fraction of the span: COMSOL's Automatic
  /// constraint, measured as exactly 0.1 of the span
  /// (doc/evidence/stepping-sweep.md, the automatic maximum step).
  inline constexpr double automatic_max_step_fraction = 0.1;

  /// The smallest step a run may be reduced to before it gives up, as a fraction
  /// of the span: the library's own guard, with no COMSOL counterpart.
  inline constexpr double default_min_step_fraction = 1e-10;

  /// Where the solver's own steps land: COMSOL's "Steps taken by solver".
  /// `free` lets the stepping choose, `intermediate` puts at least one step in
  /// every interval between two output times, `strict` ends every step at an
  /// output time.
  enum class Placement { free,
    intermediate,
    strict };

  /// What the run reports: COMSOL's "Times to store". `interpolate` stores the
  /// values at the output times, `closest` the first solver step at or past each
  /// of them — one value per step, so a step that reaches several stores them as
  /// one (doc/evidence/dofs.md) — and `steps` the solver's own steps.
  enum class Storage { interpolate,
    closest,
    steps };

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

  using MaximumStep
    = std::variant<AutomaticMaxStep, ConstantMaxStep, ExpressionMaxStep>;

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

  /// One dependent variable of the model: its share of the state, the scale
  /// COMSOL's solver uses for it, and whether its absolute tolerance follows the
  /// solution.
  ///
  /// The scale divides the solution in the scaled criterion, `Y = U / scale`
  /// (doc/evidence/tolerance-axes.md). The library cannot derive it — COMSOL
  /// analyses the model — so the caller passes the value the solver's log prints.
  struct Variable {
    /// The number of degrees of freedom of the variable.
    std::size_t dofs = 0;

    /// COMSOL's "scale for dependent variable", a divisor of the solution.
    double scale = 1.0;

    /// Whether the scaled absolute tolerance follows the solution. COMSOL sets
    /// this for a variable whose scaling method is Automatic while "Update
    /// scaled absolute tolerance" is on; every other combination leaves
    /// `A_s = A_k` (doc/evidence/tolerance-axes.md).
    bool scaled_tolerance = true;
  };

  /// The absolute tolerance: COMSOL's "Tolerance method". `factor` makes it a
  /// factor of the relative tolerance, `manual` a value of its own.
  struct AbsoluteTolerance {
    enum class Kind { factor,
      manual };
    Kind kind = Kind::factor;
    /// The factor of the relative tolerance, or the manual value.
    double value = 0.1;
  };

  /// The local error test: COMSOL's absolute-tolerance settings.
  struct ErrorTest {
    /// The relative tolerance `R`.
    double relative = 1e-2;

    /// The absolute tolerance `A_k`.
    AbsoluteTolerance absolute {};

    /// Whether the criterion is applied to the scaled variables `Y = U / scale`,
    /// COMSOL's "Tolerance method: Scaled", or to the variables as they are.
    bool scaled = true;

    /// The dependent variables of the state, in the order their degrees of
    /// freedom appear. An empty list means one variable over every degree of
    /// freedom.
    std::vector<Variable> variables {};
  };

  /// The time stepping of a run.
  ///
  /// The defaults are COMSOL's, read back from a COMSOL 6.2 Time-Dependent
  /// Solver node (see doc/comsol-alignment.md): the relative tolerance 0.01, the
  /// absolute tolerance `0.1` times it, scaled by the magnitude of each variable,
  /// the first step `1e-3` of what remains of the interval, the largest step
  /// `0.1` of the span, the orders 1 to 5, and the interpolation of the end time.
  struct Settings {
    /// How the run takes its steps: the solver's, from the local error, or the
    /// model's own at a fixed step and order.
    std::variant<AdaptiveSteps, ManualSteps> steps {};

    /// The tolerance the local error is held to.
    ErrorTest error {};

    /// COMSOL's "Fraction of initial step for Backward Euler": the size of each
    /// of the two artificial steps a consistent initialization takes, as a
    /// fraction of the run's first step. A start that is given the initial
    /// derivative takes none.
    double initialization_fraction = 1e-3;

    /// Which times the run reports.
    Storage storage = Storage::interpolate;

    /// COMSOL's "Store every Nth step"; the initial state is step 0 and the
    /// modulus is anchored on it (doc/evidence/stepping-sweep.md).
    int store_every = 1;
  };

  /// Reject a setting a run cannot be built from. Each message names the field.
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
    for (const Variable& variable : settings.error.variables)
      require("error.variables[].scale", finite_positive(variable.scale));

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
    require("initialization_fraction", finite_positive(settings.initialization_fraction));
    require("store_every", settings.store_every >= 1);
  }

} // namespace atstep
