# Design

`settings.hpp`, `stepper.hpp` and `second_order.hpp` are the public headers, with
`atstep.hpp` as their umbrella, and everything they need is under `detail`:

| Header                  | Responsibility                                        |
| ----------------------- | ----------------------------------------------------- |
| `settings.hpp`          | Configuration and validation                          |
| `stepper.hpp`           | The solve request, the run, and its output            |
| `second_order.hpp`      | The second-order-in-time adapter                      |
| `detail/state.hpp`      | State arithmetic and the default weighted RMS norm    |
| `detail/controller.hpp` | Scalar error estimates, step-size and order decisions |
| `detail/output.hpp`     | Step placement against the output times               |

The `detail` headers are installed like the rest, because a header-only library
cannot do without them, and they carry no compatibility promise: they are the
implementation of the public headers, in the namespace `atstep::detail`.

## Ownership and solve boundary

`TimeStepper` owns its settings, its solve callback, its norm, and the state of
a run in four parts, grouped by lifetime:

| Member        | What it is                                                         |
| ------------- | ------------------------------------------------------------------ |
| `accepted_`   | the accepted state, derivative, time, order and difference columns |
| `controller_` | the step size and order proposed for the next step                 |
| `work_`       | five buffers reused by successive attempts                         |
| `returned_`   | the value last returned, which the absolute weight is measured on  |

All four are owned by the stepper. Accepted history, controller decisions and
the value last returned persist across steps; work-buffer contents do not.
`returned_` is separate because COMSOL measures the scaled absolute tolerance
against the value it last returned to the caller, which may be interpolated.
`controller_` is built by `start`, once the run's span is known. State, time,
step size, order and interpolation require a successfully initialized run.

A step's own coefficients — its size, order, the factors the trial reads the
table through, the distances back to the levels, IDA's `sigma`, and its error
estimates — are a `Trial`, created and consumed inside one pass of `take_step`.
A refused step is a `continue` of that loop, so the stepper keeps no pending
attempt and the caller sees only accepted steps.

An integrator is therefore neither copyable nor movable: it is one run's state
from `start` to its last step, and `start` begins it again.

A `StepRequest` is constructed only by the stepper and borrows its buffers for
one solve callback. It is neither copyable nor movable; storing a pointer or
reference still does not extend its lifetime. There is one formulation, the
fixed-leading BDF history of a first-order equation:

```text
y' = cj*y + known
M y' + K y = f  =>  (cj*M + K)*y = f - M*known
```

`derivative_history()` returns `known` by reference, and `derivative(y, out)`
evaluates the same level in the form IDA forms it, `yypredict + cj*(y - ypred)`,
into a buffer of the caller's: the correction is small where the predictor is
good, which is where the direct form cancels. The request therefore borrows the
predictor and its derivative as well, and keeps no buffer of its own. The solve
starts from the predictor and returns convergence and iteration information.

The private `solve` function bounds the request's lifetime. During that call,
`work_.scratch` holds the known derivative term. After it returns, `work_.correction`
holds the final solution minus predictor, and scratch is reused for the error
estimates, the accepted derivative and the dense output a store reports. No
buffer needs a separate owning view or a second allocation for each use.

A second time derivative is not a second formulation: it is an equation the
adapter writes as a first-order system over the same kernel (see below). A third
is the caller's to write that way itself.

## Accepted history

Integration stores scaled divided-difference columns and step distances. A step
size rescales every column the scheme reads by one factor: the
rescaled column of order `j` is `beta_j * phi_j`, the accepted column scaled. An
attempt therefore keeps the factors and no column at all — it reads the accepted
table and writes nothing, so a step the solve or the error test refuses leaves it
exactly as the accepted steps left it. Accepting rescales the columns in place
and writes the factors in. The lead coefficient is `H_k/h`, where
`H_k = 1 + 1/2 + ... + 1/k`. Error estimates use the solved-state correction and
the same factors.

The first step's column is the derivative consistent initialization measured, at
the size the first step proposes — IDA's `phi[2] = h yp_0`. It belongs to the
accepted history, not to an attempt, so a refused first step changes nothing.

The table is therefore never copied and never written speculatively: an attempt
costs a few doubles, not one state. Dense output reads the accepted columns.

The private `History` owns the accepted state and derivative alongside the
table, time and interpolation order. Older solution levels are not stored.
The solved state and `phi[0]` remain separate: the latter is accumulated from
differences and can differ in the last bits, which matters for reproducible
error weights and trajectories.

Accepting assigns into existing buffers instead of constructing temporary
copies of the state and derivative. With fixed-size `std::vector` states, these
assignments reuse capacity. There is no rollback table and no exception handler:
if a state operation throws while the level, the derivative and the columns are
written, the history is partly written, and only `start` rebuilds it. That
differs from a failed solve or an error-test rejection, which the attempt loop
covers with a `continue` and which leaves accepted history untouched.

`start` resets the level, the derivative, the table, the controller and the
counts. Its two backward-Euler probes accept nothing until the
second. Their times must be finite, advance time and remain inside the run;
their derivative weight must be finite. A failed restart cannot expose old
history paired with newly assigned run bounds.

The stepper holds five work buffers, seven difference columns, and the accepted
state, derivative and the value last returned: fifteen state objects, whatever the
order and whatever the equation. A state's storage is allocated when a step
first writes it. Initialization allocates columns `0..max_order+1`; unused higher
columns remain empty. Dense output allocates only the buffer its caller passes.

`StepController` consumes scalar estimates independently of state and equation
types. It implements startup growth, settled-step order selection and retries.
Maximum-step arithmetic preserves rounding at the transition to a constant
step: a one-bit overshoot can reset the settled-step counter.

## Second-order equations

`SecondOrderStepper` advances `M u'' + C u' + K u = f(t)` with the first-order
kernel, the way COMSOL's IDA path does: "extra variables are internally
introduced so that it is possible to form a first-order time-derivative system",
and "the linearized system is reduced to the original size with the matrix
E + σD + σ²K" ([COMSOL 6.4, About the Time-Dependent
Solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html)).

The pair the adapter integrates in the kernel's state is `(u, v = u')`:

```text
u' - v            = 0
M v' + C v + K u  = f
```

The caller's state is the state, and the first derivative is a companion the
adapter keeps in a `detail::Augmented` state: prediction, levels, table, dense
output and the error test of the kernel all run on the pair without knowing
which half is which.

The adapter owns only three additional state-sized buffers: the first
and second derivative predictors and the second derivative's known term. Its
request borrows the state predictor, first known term and derivative
scratch from the kernel request. These borrowed references live for the same
solve callback; no separate reference bundle is stored. A failed solve returns
immediately, without trying to reconstruct the derivative from a failed candidate.

The step's solve is then the condensed system of the state alone. With
`cj` the step's first-derivative weight, the two derivative levels are the
affine forms

```text
u'  = cj    u + known(1)
u'' = cj^2  u + known(2),   known(2) = cj known(1) + known(v),
```

where `known(v)` is the known term of the first derivative's own formula.

and substituting them into the equation leaves

```text
(K + cj C + cj^2 M) u_{n+1} = f(t_{n+1}) - M known(2) - C known(1).
```

A nonlinear residual is condensed the same way, and its Jacobian is
`dF/du + cj dF/du' + cj² dF/du''`. The linear system is the size it was before
the augmentation, and the adapter is the only place that knows about it. The
solve returns the state; the adapter derives the level's first derivative from
the kinematic constraint, `v = cj u + known(1)`, which is what makes the pair at
every accepted level a consistent one.

`run` reports the first derivative with every stored state, `derivative()` is the
newest level's, and `interpolate(t, 1)` is the dense first derivative: the kernel's dense
derivative of the pair, whose derivative block is the state's derivative.

Two consequences of the layering are deliberate:

- The error test measures the state alone: the adapter wraps the caller's norm
  so that the pair reaches it as its state's half. The first derivative's error
  is the derivative of the state's, so a norm over it would hold the step to the
  accuracy of a derivative — a factor `1/h` tighter — rather than to the
  accuracy of the solution. State augmentation and error-norm selection are
  therefore separate decisions: the caller's norm describes its own variables,
  and the pairs never reach it. COMSOL excludes its auxiliary velocity for the
  same reason, except that its documentation includes it while the consistent
  initialization is on; that sentence is still unmeasured
  (doc/refactor-plan.md, E05).
- No equation ever carries two derivatives into the kernel, so there is one
  integration path, one order selection and one error estimate. Differencing the
  state levels twice — the classical stencil of version 1.0 — is a
  different and worse discretization: the BDF polynomial through the levels has
  a parasitic root of multiplicity two, and the run is first order however high
  the scheme's order is.

The initialization is the kernel's: two artificial backward-Euler steps of the
pair. The initial derivative is a second initial value, and the increment of the
state over the probe is what the first step is sized from, as the
measured derivative is for a first-order run.

## Output and error weights

Step placement and stored output are separate settings. Strict mode lands on
each requested time; intermediate mode ensures a step in each output interval.
Landing logic preserves a step already reaching its target, otherwise shortens
it with a rounding margin or splits the final two steps. `StepSchedule` holds
the output times and the two modes, and answers the two questions a step asks of
them: the size it must be cut to, and the values to report once it ended.
It owns a copy of the times and emits each output synchronously, without
allocating a list of pending output records. Its one piece of state is whether
the interval that starts at the last stored output already holds a step: an
accepted step that crossed an output time counts as the landing of the next
interval, and the level a consistent initialization leaves sits `2 * probe` past
the start time without being a step, so that flag cannot be read off the clock.

`run` consumes a complete schedule immediately after `start`. Calling it after
`advance` or after a completed `run` is rejected, so a later state cannot be
mislabelled as the initial output. Caller-driven stepping uses `advance`;
starting a new scheduled run requires `start` again.

Dense output uses committed BDF columns. At the newest node the returned
derivative is the solve derivative; it can differ from the dense polynomial's
derivative after a step change. `interpolate` writes into a state the caller
passes, so it keeps nothing and invalidates nothing: the same buffer can serve
every call. The states a hook is handed do use the run's own buffers and are
valid only during the call.

The default norm is COMSOL's error test: the weighted RMS over the state,
averaged per dependent variable first and per degree of freedom second, with
`Y = U / scale` and `A_s = A_k (beta + rms(Y))` for a variable whose absolute
tolerance follows the solution and `A_s = A_k` otherwise
(doc/evidence/tolerance-axes.md). It uses the last accepted internal state for
the relative weights and the most recently stored state for the scaled absolute
weight, so interpolated output affects later error weights, as the COMSOL tests
cover. A custom norm receives the error, the accepted state, the returned state
and the accepted-step count; the kernel never mutates it.

Invalid dimensions and nonfinite configuration are rejected. A nonzero error
with zero weight has infinite norm. A step must be finite, positive and large
enough to advance floating-point time.

## Validation boundaries

Unit tests cover controller rules, convergence, the second-order adapter,
rejection, restart, invalid input and aliased derivative output. CTest also runs
the three examples and `test/acceptance`, the COMSOL conformance target of the
rewrite. The [scalar COMSOL regressions](comsol-alignment.md) compare 1,336
recorded output samples of 21 configurations, including times, states,
derivatives and failure counts, 1,904 samples of 23 stepping configurations, 682
samples of 11 two-field configurations, and the 17 tolerance configurations of
[the tolerance evidence](evidence/tolerance-axes.md) at the level that evidence
establishes. The orders in COMSOL's solver log are printed as a diagnostic.
Reference data is never injected into the numerical solve.

The [thermoelastic regression](multiphysics.md) additionally compares 161
COMSOL states of a reciprocal three-variable coupled system, including all
internal times, states, derivatives and solver counts. These regressions
validate the recorded configurations, rather than arbitrary nonlinear or
continuum multiphysics models. Regenerating references requires COMSOL;
running the C++ tests does not.
