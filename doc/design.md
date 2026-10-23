# Design

Two headers are public, and everything they need is under `detail`:

| Header                  | Responsibility                                        |
| ----------------------- | ----------------------------------------------------- |
| `settings.hpp`          | Configuration and validation                          |
| `stepper.hpp`           | The solve request, the run, and its output            |
| `detail/state.hpp`      | State arithmetic and the default weighted RMS norm    |
| `detail/bdf.hpp`        | Accepted levels, polynomial stencils, error estimates |
| `detail/controller.hpp` | Scalar step-size and order decisions                  |
| `detail/output.hpp`     | Step placement against the output times               |

The `detail` headers are installed like the rest, because a header-only library
cannot do without them, and they carry no compatibility promise: they are the
implementation of the two public headers, in the namespace `atstep::detail`.

## Ownership and solve boundary

`TimeStepper` owns its settings, its solve callback, its norm, and the state of
a run in four parts:

| Member        | What it is                                             |
| ------------- | ------------------------------------------------------ |
| `levels_`     | the accepted steps, newest first                       |
| `table_`      | the scaled differences a first-order run reads instead |
| `controller_` | the step size and order proposed for the next step     |
| `work_`       | the eleven buffers of a step                           |

A step's own coefficients — its size, order, the factors the trial reads the
table through, the distances back to the levels, IDA's `sigma`, and its error
estimates — are a `Trial`, created by the attempt that computes them and
consumed by the same call. No attempt outlives its step, so the lifetime of
every temporary is a C++ scope rather than a runtime flag: the stepper keeps no
"pending attempt", no attempt identity, and no handle for a caller to settle.
A step whose solve or error test is refused is retried inside `advance` and
`run`, and the caller sees only accepted steps.

An integrator is therefore neither copyable nor movable: it is one run's state
from `start` to its last step, and `start` begins it again.

A `StepRequest` borrows its references for the duration of the solve, and a
state it exposes is the stepper's buffer. Its derivative formulation is a
variant: fixed-leading BDF or a classical polynomial stencil. Callers use the
same affine derivative interface for both:

```text
y' = cj*y + known
M y' + K y = f  =>  (cj*M + K)*y = f - M*known
```

`derivative_history()` returns `known` by reference, and `derivative(y, out)`
evaluates `cj*y + known` in correction form into a buffer of the caller's,
reducing cancellation where the predictor is good. A derivative level above the
first is written by `derivative_history(level, out)`. The solve starts from the
predictor and returns convergence and iteration information; matrix assembly and
algebraic solvers belong to the caller. Custom states implement `axpy`; custom
norms implement `norm(error, reference)`.

## Accepted history

First-order integration stores scaled divided-difference columns and step
distances. A step size rescales every column the scheme reads by one factor: the
rescaled column of order `j` is `beta_j * phi_j`, the accepted column scaled. An
attempt therefore keeps the factors and no column at all — it reads the accepted
table and writes nothing, so a step the solve or the error test refuses leaves it
exactly as the accepted steps left it. Accepting rescales the columns in place
and writes the factors in. The lead coefficient is `H_k/h`, where
`H_k = 1 + 1/2 + ... + 1/k`. Error estimates use the solved-state correction and
the same factors.

The table is therefore never copied and never written speculatively: an attempt
costs a few doubles, not one state. Dense output reads the accepted columns.

The level behind it is the table alone: a first-order run keeps one accepted
level, which `previous()` and the error weights read, and no older one. The
levels are a ring of the depth the formulation needs — one, or `max_levels` for
a higher-derivative run, whose classical stencil reads them — so an accepted
step costs one state copy, no allocation and no shifting of older levels.
`previous(level)` is what starts that stencil.

A level pairs its state with its time in one object, so a copy that throws
cannot leave the two out of step; a level is constructed before the ring is
written, which is what leaves the accepted history alone when it throws.

`start` resets the levels, the table, the controller, the tolerance indexing,
and the correction the estimate at the order above compares against.

The stepper holds the eleven buffers of a step, the seven columns and the
accepted levels as states: nineteen objects for a first-order run whatever the
order, and twenty-five for a higher-derivative one. A state's storage is
allocated when a step first writes it, so a first-order run at order `k`
allocates the columns `0..k+1` and the buffers the estimates read.

`StepController` consumes scalar estimates independently of state and equation
types. It implements startup growth, settled-step order selection and retries.
Maximum-step arithmetic preserves rounding at the transition to a constant
step: a one-bit overshoot can reset the settled-step counter.

`StepFormula` serves classical stencils, including higher derivatives. It is
not reconstructed for first-order fixed-leading steps. Lagrange values and
derivatives share one evaluator. Counts, node uniqueness and finiteness are
checked before accessing fixed-size buffers; extra step history is bounded
to the supported order range.

## Output and error weights

Step placement and stored output are separate settings. Strict mode lands on
each requested time; intermediate mode ensures a step in each output interval.
Landing logic preserves a step already reaching its target, otherwise shortens
it with a rounding margin or splits the final two steps. `StepSchedule` holds
the output times and the two modes, and answers the two questions a step asks of
them: the size it must be cut to, and the values to report once it ended.

Dense output uses committed BDF columns. At the newest node the returned
derivative is the solve derivative; it can differ from the dense polynomial's
derivative after a step change. Returned state references use internal buffers
and must be copied if retained across callbacks or interpolation calls.

The default norm uses the last accepted internal state for relative weights
and the most recently stored state for scaled absolute weights. Interpolated
output therefore affects later error weights, as covered by the COMSOL tests.
Custom norms receive only the error and accepted internal state.

Invalid dimensions and nonfinite configuration are rejected. A nonzero error
with zero weight has infinite norm. A step must be finite, positive and large
enough to advance floating-point time.

## Validation boundaries

Unit tests cover polynomial reproduction, controller rules, convergence,
rejection, restart, invalid input and level-copy failure. CTest also runs
four examples. The [scalar COMSOL regressions](comsol-alignment.md) compare 21 recorded
scalar configurations, totaling 1,336 output samples, including times, states,
derivatives, failure counts and logged orders. Reference data is never injected
into the numerical solve.

The [thermoelastic regression](multiphysics.md) additionally compares 161
COMSOL states of a reciprocal three-variable coupled system, including all
internal times, states, derivatives and solver counts. These regressions
validate the recorded configurations, rather than arbitrary nonlinear or
continuum multiphysics models. Regenerating references requires COMSOL;
running the C++ tests does not.
