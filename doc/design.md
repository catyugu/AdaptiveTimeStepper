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
distances. Each attempt rescales a copy, forms the predictor and derivative
predictor, then solves. The leading derivative coefficient is `H_k/h`, where
`H_k = 1 + 1/2 + ... + 1/k`. Error estimates use the solved-state correction
and trial columns.

Acceptance commits trial history; rejection preserves accepted history.
Only one set of trial buffers exists. Handle the latest `Attempt` before
creating another: attempts are diagnostics, not saved state snapshots.
`start` resets numerical history and tolerance indexing.

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
three examples. The [COMSOL regressions](comsol-alignment.md) compare 21 recorded
scalar configurations, totaling 1,336 output samples, including times, states,
derivatives, failure counts and logged orders. Reference data is never injected
into the numerical solve.

These regressions validate the recorded configurations, rather than arbitrary
nonlinear or coupled models. Regenerating references requires COMSOL; running
the C++ tests does not.
