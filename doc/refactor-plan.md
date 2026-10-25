# AdaptiveTimeStepper: the rewrite

The interface the library is being brought to, what has already moved, and what
is left. The measurements are in [evidence/](evidence/); the alignment with
COMSOL is in [comsol-alignment.md](comsol-alignment.md).

One rule outranks the others: the 21+23 recorded single-variable configurations,
the 11 two-field configurations and the 161 thermoelastic samples are the
regression net. They stay green to the last digit at every stage, and a change
that moves them without a measurement behind it is wrong.

## The interface as it stands

### Settings: COMSOL's Time-Dependent Solver node

```cpp
namespace atstep {

  inline constexpr int max_order = 5;

  enum class Placement { free, intermediate, strict };
  enum class Storage { interpolate, closest, steps };

  struct AutomaticMaxStep { };
  struct ConstantMaxStep { double seconds; };
  struct ExpressionMaxStep { std::function<double(double)> seconds; };
  using MaximumStep = std::variant<AutomaticMaxStep, ConstantMaxStep, ExpressionMaxStep>;

  struct AdaptiveSteps {
    Placement placement = Placement::free;
    MaximumStep max_step = AutomaticMaxStep {};
    int min_order = 1;                       // COMSOL's "Minimum BDF order"
    int max_order = 5;
    std::optional<double> initial_step {};   // unset: the automatic rule
    bool interpolate_end_time = true;
    double min_step_fraction = 1e-10;        // the library's own guard
  };

  struct ManualSteps {
    double step = 0.0;
    int order = 2;
    std::optional<double> startup_fraction {};
    std::optional<double> startup_growth {};
    bool interpolate_end_time = true;
  };

  struct Variable {                  // one dependent variable
    std::size_t dofs = 0;
    double scale = 1.0;              // COMSOL's printed scale; divides the solution
    bool scaled_tolerance = true;    // automatic scale *and* "Update scaled absolute tolerance"
  };

  struct ErrorTest {
    double relative = 1e-2;
    struct Absolute { enum class Kind { factor, manual }; Kind kind; double value; };
    Absolute absolute {};            // "Factor" of the relative tolerance, or a value
    bool scaled = true;              // "Tolerance method"
    std::vector<Variable> variables; // empty: one automatic variable over the state
  };

  struct Settings {
    std::variant<AdaptiveSteps, ManualSteps> steps {};
    ErrorTest error {};
    double initialization_fraction = 1e-3;
    Storage storage = Storage::interpolate;
    int store_every = 1;
  };
}
```

`start(time, initial, end)` takes the consistent initialization;
`start(time, initial, derivative, end)` takes the caller's derivative and takes
none. The overload is the decision, so the settings hold only the fraction.

### The solve request

```cpp
template <class State> class StepRequest {
public:
  const double time;  const double step;  const int order;
  const std::size_t accepted_steps;  const bool initializing;

  double derivative_weight() const;            // cj
  const State& derivative_history() const;     // known
  const State& predictor() const;
  const State& derivative_predictor() const;
  void derivative(const State& y, State& out) const;   // out = yypredict + cj*(y - ypred)
  const State& accepted_state() const;  const State& returned_state() const;
};
```

The request borrows five states and holds no buffer of its own: the two the
correction form needs, the known term, and the two the error weights are
measured against — which a model that repeats the weight for its own nonlinear
convergence test also needs. `derivative` writes the correction form *into the
caller's* argument, so it needs no scratch. That form is not a convenience: IDA
solves for the correction and forms `y' = yypredict + cj*ee`, which is why the
recorded `irregular_intermediate` trajectory only reproduces with it.

### The run

```cpp
template <class State> class TimeStepper {
public:
  using Solve = std::function<SolveOutcome(const StepRequest<State>&, State&)>;
  using Norm  = ErrorNorm<State>;   // (error, accepted, returned, accepted_steps) -> double
  static Norm default_norm(const ErrorTest&);

  explicit TimeStepper(const Settings& settings, Solve solve, Norm norm = {});

  void start(double time, const State& initial, double end_time);
  void start(double time, const State& initial, const State& derivative, double end_time);

  StepRecord advance();          // one accepted step, cut to the placement and the end
  bool finished() const;
  RunReport run(std::span<const double> times, RunHooks<State> hooks = {});

  double time() const;  const State& state() const;
  double step_size() const;  int order() const;
  void interpolate(double time, State& out, int derivative = 0);
  const RunReport& report() const;
};
```

One driver: `run` is `while (!finished()) { advance(); store(); }`. `advance`
cuts the step to the next placement boundary and, when the end is not
interpolated, to the end time, so a caller-driven loop lands where the
scheduled one does.

The norm is a value the caller may replace. The default is COMSOL's weight,
built from the settings; a state that exposes its dofs gets it for free, and a
state of its own arithmetic passes its own. There is one protocol, no type is
compared, and the kernel never mutates the norm.

### Second order

`SecondOrderRequest` is a value the adapter builds: `time`, `step`, `order`,
`accepted_steps`, `initializing`, and the two affine forms. It borrows the two
known terms and nothing else; no friendship anywhere. The adapter keeps the pair
`(u, v = u')` as the kernel's state — that is what keeps one discretization —
and the caller's norm is applied to the state's half.

## What has landed

| Stage | What | The record |
| --- | --- | --- |
| 0 sessions | `advance` lands on the end, `run` is its loop, `finished()` and `report()` exist, `start` picks the initialization by arity | `test/acceptance/test_session.cpp` |
| 0 settings | `update_scaled`, `include_algebraic` and `Field::algebraic` are gone; `initialization` is a number | all suites |
| 1 the weight | `WeightedRmsNorm` is the measured law: per variable RMS and maximum, `Y = U/scale`, the `A_s` branch | 17 acceptance cases: every recorded first step, six series to 9e-16, the three controls |
| 2 the request | no friend, no scratch, no `step_index` field, one accessor per borrowed state, the correction form in place | 21+23+11+161 recorded samples, unchanged |
| 2 the norm | one `ErrorNorm` protocol, defaulted by capability, `ValueNorm` and `configure_norm` gone | all suites |
| 2 second order | value request, no `SecondOrderWork`, no `Augmented::size` | all suites |
| 3 the scheduler | the times are owned; the crossing flag is named and its reason recorded | the 23 stepping configurations, including the two non-uniform schedules |
| 4 the controller | the first step is computed by `begin`, the dead manual branch is gone | `test_controller.cpp`, the 10 manual configurations |

## What is left

### The second-order error test (E05)

Whether the auxiliary velocity joins the error test while consistent
initialization is on: COMSOL 6.4 says it does, the adapter excludes it. The
decision waits for the run that settles it: one dof, `u'' + 0.1 u' + u = 0`,
physics forced to `timemethod = bdf` with the readback asserted, replayed with
both norms.

### The duplicate sweep

`test/reference/comsol_trajectories.hpp` (1,401 lines) replays the same model as
`comsol_stepping.hpp` with the same Newton emulation: I replayed the 23 stepping
cases with `comsol_decay_solve.hpp` and they agree to 1.75e-12 in value and
3.7e-10 in derivative, so the value comparison can move there. Its only unique
content is its `logged_orders`, which are wall-clock-sampled log rows the
replacement already documents as irreproducible. Delete the header, its test,
`GeDecayTransient.java` and `run_decay.py` once the comparison has moved.

### The documents

This file is the last one that describes work in progress. When E05 and the
sweep land, it goes: the code and `evidence/` are then the record.

## Open questions

Each waits for a measurement; none is implemented on a guess.

| ID | Question | What settles it |
| --- | --- | --- |
| Q01 | `beta`'s argument for a multi-dof variable at `j > 0`: the variable's maximum norm (documented) or its RMS (indistinguishable at `j = 0` and on one-dof variables). Also whether `A_s` is updated on every accepted step. | An automatic-scale model whose scale is exact, or a readback of the computed scale after solving, replayed for its whole series. |
| Q02 | Does the auxiliary velocity join the second-order error test while consistent initialization is on? | E05. |
| Q03 | Does COMSOL count physics interfaces or dependent variables as `M`? | A model with two variables in one interface, or one variable in two (`field-norm.md`, section 6). |
| Q04 | Does a rejected attempt count toward "Store every Nth step"? No recorded rejection happens before a stored step. | A case whose rejection falls between two stored steps. |
| Q05 | `Closest` on a non-uniform schedule: the first step that reaches the time, or the nearest step end? | One run with `tout = tlist` and uneven times. |
| Q06 | Can the multi-variable tolerance trajectories be replayed exactly? The record carries COMSOL's inexact Newton from the first size change. | Emulating the Jacobian reuse for the field model, as `comsol_decay_solve.hpp` does for the scalar one. |

## What the audit found, and what stands

The audit of commit `65b9f1f` is in the history. Its findings that still shape
the tree:

- **The weight was one special case.** Fixed by stage 1; the 69 regression cases
  did not move, which is the check that they were that special case.
- **Three settings did nothing.** Gone in stage 0.
- **The caller-driven run had no end.** Fixed in stage 0.
- **The second-order error test contradicts the manual.** Open, E05.
- **Documents contradicted the code** about the per-field average and the
  tolerance law. Fixed.
- The over-design it named — the friend-based request, the norm configured by
  type in three places, the exception protocol, the dead `first_step` branch,
  the two `start` overloads policing a variant — is gone.

Two of its judgments the measurements overturned, which is why the code differs
from that audit:

- It proposed deleting `StepRequest::derivative`. IDA forms `y'` from the
  correction, and the recorded `irregular_intermediate` trajectory reproduces
  only that way. The method stays; what went is the mutable scratch it used.
- It called the output schedule's `landed_` flag avoidable. It is not: the level
  a consistent initialization leaves sits `2 * probe` past the start time and is
  not a step, so "this interval already has one" is state, not a function of the
  clock. It survives as `StepSchedule::interval_has_step_`, with the reason
  written down.

## What is not a defect

Worth keeping in view, because most of the numerics is right and the defects
were at the edges:

- The variable-step BDF machinery is faithful to IDA: the `phi/psi/beta/sigma`
  table, the predictor, `ck = max(alpha, abs(alpha + alphas - alpha0))`, the
  `k-1`/`k-2`/`k+1` estimates, the order selection, the startup phase, the
  first-failure reduction. Each was checked against `ida.c` and `ida_impl.h`, and
  the `ida::` constants are IDA's own, not fitted.
- The COMSOL stepping settings are measured and implemented: the automatic
  maximum step, the automatic initial step with its derivative bound, the manual
  startup curve with its per-order defaults, the constant and expression maximum
  step, the three placements, the three storage modes, store-every-Nth anchored
  on the initial state, and the end-time interpolation.
- The solve boundary is right: the caller owns assembly and the algebraic solve,
  the library owns prediction, history, error control and output, and the state
  needs one operation (`axpy`).
- The end-of-run interpolation, the initial-state label at the start time, and
  the modulus anchored on step 0 are measured and right.
- The recorded trajectories are reproduced because the library keeps the phi
  accumulation rather than the solved state for its weights, and the same
  arithmetic order. Do not "simplify" those away.

## What the library does not do

No integrator framework: no base class, no virtual dispatch, no policy matrix.
No equation layer: no Newton, no sparse assembly, no expression parser, no unit
system — the solve callback is the boundary. No copy of COMSOL's nonlinear
solver strategy: that belongs to the model, and the reference models emulate
their own. No explicit methods, no generalized alpha, no events, no local time
stepping, no remeshing. No field reserved for a path that does not exist: every
branch has an implementation, a test and a COMSOL case, or it is not merged.
