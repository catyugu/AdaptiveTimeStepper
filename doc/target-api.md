# The target interface of AdaptiveTimeStepper

The refactor plan says what to change; this says what it is changed *to*. It
is the contract the rewrite is written against and the acceptance tests are
written against, so it is written down before the code moves.

Nothing here is implemented yet. The tests in `test/reference/test_comsol_*.cpp`
are red until the corresponding part lands, and each failure names the setting
that is missing.

## The configuration

```cpp
namespace atstep {

  inline constexpr int max_order = 5;
  inline constexpr int max_levels = max_order + 2;
  inline constexpr double default_initial_step_fraction = 1e-3;
  inline constexpr double automatic_max_step_fraction = 0.1;
  inline constexpr double default_min_step_fraction = 1e-10;

  /// Where the solver's own steps land: COMSOL's "Steps taken by solver".
  enum class Placement { free, intermediate, strict };

  /// What the run reports: COMSOL's "Times to store".
  enum class Storage { interpolate, closest, steps };

  /// The largest step: COMSOL's "Maximum step constraint".
  struct AutomaticMaxStep { };

  struct ConstantMaxStep {
    /// The maximum step, in the run's time unit.  COMSOL's "maxstepbdf" is an
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

  struct AdaptiveSteps {
    Placement placement = Placement::free;
    bool interpolate_end_time = true;
    MaximumStep max_step = AutomaticMaxStep {};
    int min_order = 1;    // COMSOL's "Minimum BDF order": 1 or 2
    int max_order = 5;    // 1 to 5
    /// COMSOL's "Initial step": unset means automatic.
    std::optional<double> initial_step {};
    double min_step_fraction = default_min_step_fraction;
  };

  struct ManualSteps {
    double step = 0.0;
    int order = 2;
    /// COMSOL's per-order defaults; unset means the default for `order`.
    std::optional<double> startup_fraction {};
    std::optional<double> startup_growth {};
    bool interpolate_end_time = true;
  };

  /// One dependent variable.
  struct Field {
    std::size_t dofs = 0;
    double scale = 1.0;
    bool algebraic = false;
  };

  struct Layout {
    std::vector<Field> fields;
  };

  struct AbsoluteTolerance {
    enum class Kind { factor, manual };
    Kind kind = Kind::factor;
    double value = 0.1;
  };

  struct ErrorTest {
    double relative = 1e-2;
    AbsoluteTolerance absolute {};
    bool scaled = true;
    /// COMSOL's "Update scaled absolute tolerance".
    bool update_scaled = true;
    /// COMSOL's "Error estimation": whether the algebraic dofs are held to the
    /// tolerance.
    bool include_algebraic = true;
    Layout layout {};
  };

  struct BackwardEuler {
    double fraction = 1e-3;   // COMSOL's "bwinitstepfrac"
  };

  /// No consistent initialization: the initial values are taken as they are,
  /// and the caller supplies the initial derivative.
  struct AssumeConsistent { };

  struct Settings {
    std::variant<AdaptiveSteps, ManualSteps> steps {};
    ErrorTest error {};
    std::variant<BackwardEuler, AssumeConsistent> initialization = BackwardEuler {};
    Storage storage = Storage::interpolate;
    /// COMSOL's "Store every Nth step"; the initial state is step 0 and the
    /// modulus is anchored on it.
    int store_every = 1;
  };

  void validate(const Settings& settings);

} // namespace atstep
```

## What each setting means

| Setting | COMSOL's name | Rule |
| --- | --- | --- |
| `Placement` | "Steps taken by solver" | `free` lets the stepping choose, `intermediate` puts a step in every interval between two output times, `strict` ends every step at an output time. |
| `Storage` | "Times to store" | `interpolate` stores the values at the output times, `closest` the first solver step reaching each of them, `steps` the solver's own steps. |
| `ConstantMaxStep::seconds` | "maxstepbdf" | An absolute time value. The library's present field is a *fraction of the span*, and the two differ by the span (ten at the recorded span of 10), which is why the recorded reference cannot see the difference. |
| `ExpressionMaxStep` | "maxstepbdf" as an expression | Evaluated while solving; the absolute value is the cap. |
| `AdaptiveSteps::initial_step` | "Initial step" | Unset means automatic: a fraction of the interval that remains *after* the consistent initialization, further limited by `DT * yp_norm <= 1/2`. The library's present rule takes the fraction of the whole span and so folds the initialization in twice; the measured difference is 2e-8 s at a span of 10 and 2e-7 s at a span of 100. |
| `AdaptiveSteps::min_order` | "Minimum BDF order" | 1 or 2. |
| `ManualSteps::startup_fraction`, `startup_growth` | "Initial step fraction", "Initial step growth rate" | The step grows as `min(fraction * step * growth^(k-1), step)`. Unset means the default for the order: 1 and 1.5 at the orders 1 and 2, 0.1 and 1.2 at the order 3, 0.001 and 1.1 at the order 5. |
| `ErrorTest::absolute` | "Tolerance method" | `factor` makes the absolute tolerance a factor of the relative one, `manual` a value of its own. |
| `ErrorTest::scaled` | "Scaled" / "Unscaled" | Whether the tolerance is applied to scaled variables. |
| `ErrorTest::update_scaled` | "Update scaled absolute tolerance" | Whether the absolute tolerance follows the magnitude of the solution. |
| `ErrorTest::include_algebraic` | "Error estimation" | Whether the algebraic degrees of freedom are held to the tolerance. |
| `ErrorTest::layout` | the dependent variables of the model | The partition of the state into fields. The error test averages over the fields first and over the degrees of freedom of a field second; a flat average over all of them disagrees by five to seven orders of magnitude in the step sequence. |
| `Settings::store_every` | "Store every Nth step" | The initial state is step 0 and the modulus is anchored on it. |
| `Settings::initialization` | "Consistent initialization" | `BackwardEuler` takes the two artificial steps; `AssumeConsistent` takes the initial values as they are and needs the initial derivative from the caller. |

## The current interface, and what becomes of it

| Now | Becomes | Note |
| --- | --- | --- |
| `atstep::Free`, `atstep::Intermediate`, `atstep::Strict` | `Placement` | The variants carry no data of their own; the enum says the same thing in one line. |
| `atstep::Interpolate`, `atstep::Closest`, `atstep::Steps` | `Storage` | The same. |
| `atstep::ConstantMaxStep{fraction}` | `ConstantMaxStep{seconds}` | The unit changes from a fraction of the span to an absolute time value. The name and the aggregate form stay, so a caller that writes the target value compiles against both. |
| — | `ExpressionMaxStep` | New. |
| `AdaptiveSteps::max_step_fraction` | `AdaptiveSteps::initial_step` | The automatic rule is a fraction of what remains of the interval; the optional holds an absolute value when the caller gives one. |
| — | `AdaptiveSteps::min_order` | New; the present run always starts at the order 1. |
| — | `ManualSteps::startup_fraction`, `startup_growth` | New. |
| `ErrorTest::absolute_factor` | `ErrorTest::absolute` | The kind joins the factor: a factor of the relative tolerance, or a value of its own. |
| `ErrorTest::scale` (a scale per degree of freedom) | `ErrorTest::layout` (a scale per field) | COMSOL scales a variable, not a degree of freedom; the field's degrees of freedom share the scale. |
| `ErrorTest::scaled` | `ErrorTest::scaled` | Unchanged, with `update_scaled` beside it. |
| — | `ErrorTest::include_algebraic` | New. |
| `Settings::initialization.fraction` | `Settings::initialization` (a variant) | The two artificial steps, or none. |
| — | `Settings::store_every` | New. |

## What the acceptance tests ask for

| Test case | Needs | Today |
| --- | --- | --- |
| `comsol stepping: the automatic maximum step and the initial step` | the automatic initial step over the interval that remains | the counts agree; the first stored time is 2e-8 s late and the offset compounds to 1.2e-5 s |
| `comsol stepping: the constant maximum step is an absolute value` | `ConstantMaxStep::seconds` | 37 steps instead of 204: the value is read as a fraction of the span |
| `comsol stepping: the expression maximum step` | `ExpressionMaxStep` | not expressible |
| `comsol stepping: a specified initial step` | `AdaptiveSteps::initial_step` | not expressible |
| `comsol stepping: store every Nth step` | `Settings::store_every` | not expressible |
| `comsol stepping: the manual startup curve` | `ManualSteps::startup_fraction`, `startup_growth` | not expressible |
| `comsol two fields: the error test averages over the fields` | `ErrorTest::layout` | not expressible |
| `comsol two fields: a single PDE field reproduces COMSOL` | `AssumeConsistent` | COMSOL shrinks its backward-Euler probes by the safety factor of the nonlinear solver, so the case has to be replayed from the recorded state |

The cases that pass today are the ones the rewrite must not break: the five
automatic cases of the sweep (37 steps, 1 rejection, the stored times to
2.4e-11 s), the four manual cases whose startup curve is the trivial one
(100 steps, the stored times exact), and the Global Equation case of the
two-field sweep (60 steps, the stored times to 4.0e-11 s).

## What is not part of the interface

The algebraic solve of the model, the Jacobian of the nonlinear solver, and
the consistent initialization's *step size* belong to the model layer. The
library takes the initial derivative from the caller when the initialization
is switched off, and the size of the backward-Euler probes is a fraction of
the first step, which COMSOL may shrink by the safety factor of the
nonlinear solver (the recorded PDE model does, the recorded Global Equation
model does not).
