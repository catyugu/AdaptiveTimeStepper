# Design

The library has five components:

| Header           | Responsibility                                             |
| ---------------- | :--------------------------------------------------------- |
| `settings.hpp`   | Configuration and validation                               |
| `state.hpp`      | State arithmetic and the default weighted RMS norm         |
| `formula.hpp`    | Classical polynomial interpolation and differentiation     |
| `controller.hpp` | Scalar step-size and order decisions                       |
| `stepper.hpp`    | Solve requests, accepted/trial history and output handling |

## Ownership and solve boundary

`TimeStepper` owns its settings, solve callback, norm and state buffers by value.
Each history level pairs a time and state in one container; insertion cannot
leave time and state sequences out of sync if copying a state throws.
Reset constructs replacement history before discarding existing levels.

A `StepRequest` borrows required references for the duration of the solve.
Its derivative formulation is a variant: fixed-leading BDF or a classical
polynomial stencil. Callers use the same affine derivative interface for both:

```text
y' = cj*y + known
M y' + K y = f  =>  (cj*M + K)*y = f - M*known
```

`derivative(y)` evaluates the correction form for fixed-leading BDF, reducing
cancellation. The solve starts from the predictor and returns convergence and
iteration information. Matrix assembly and algebraic solvers belong to the caller.
Custom states implement `axpy`; custom norms implement `norm(error, reference)`.

## Accepted and trial history

First-order integration stores scaled divided-difference columns and step
distances. A step size rescales every column the scheme reads by one factor:
the rescaled column of order `j` is `beta_j * phi_j`, the accepted column
scaled. A trial therefore holds the factors and no column at all — it reads the
accepted table and writes nothing, so a rejection, and an attempt a caller walks
away from, leave it exactly as the accepted steps left it. `accept` rescales the
columns in place and commits the factors. The lead coefficient is `H_k/h`, where
`H_k = 1 + 1/2 + ... + 1/k`. Error estimates use the solved-state correction and
the same factors.

The table is therefore never copied and never written speculatively: a trial
costs a few doubles, not one state. Dense output reads the accepted columns
whether an attempt is pending or not.

The level behind it is the table alone: a first-order run keeps one accepted
level, which is what `previous()` and the error weights read, and no older one.
A higher-derivative run keeps `max_levels` levels, because its classical stencil
reads them; `previous(level)` is then what starts that stencil.

Acceptance commits trial history; rejection preserves accepted history.
One attempt is pending at a time:
`attempt` returns a move-only `TrialStep` carrying the attempt's identity, and
`accept`, `reject` or `reject_convergence` settle that identity once. A handle
of an older attempt, of another stepper, or of an attempt already settled is
refused, and so is accepting a step whose solve failed or whose error test did
not pass. `start` resets numerical history, tolerance indexing and the pending
attempt.

The stepper's own buffers are the eleven of a step and the seven columns, plus
the levels the history keeps. A state is allocated when a step first writes it,
and the columns follow `max_order`: a first-order run at order 5 holds 18 states
with a dense output allocated, at order 2 15 and at order 1 14. A
higher-derivative run keeps all `max_levels` levels and holds 15 at order 5.

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
it with a rounding margin or splits the final two steps.

Dense output uses committed BDF columns. At the newest node the returned
derivative is the solve derivative; it can differ from the dense polynomial's
derivative after a step change. Returned state references use internal buffers
and must be copied if retained across callbacks or interpolation calls.

The default norm uses the last accepted internal state for relative weights
and the most recently stored state for scaled absolute weights. Interpolated
output therefore affects later error weights, as covered by the COMSOL tests.
Custom norms receive only the error and accepted internal state.

Invalid dimensions and nonfinite configuration are rejected. A nonzero error
with zero weight has infinite norm. An attempted step must be finite, positive
and large enough to advance floating-point time.

## Validation boundaries

Unit tests cover polynomial reproduction, controller rules, convergence,
rejection, restart, invalid input and history-copy failure. CTest also runs
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
