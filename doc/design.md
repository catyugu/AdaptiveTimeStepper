# Design

`settings.hpp`, `stepper.hpp`, `norm.hpp` and `second_order.hpp` are the public
headers; `atstep.hpp` is their umbrella. Everything else is under `detail` and
carries no compatibility promise.

| Header | Responsibility |
| :--- | :--- |
| `settings.hpp` | Configuration and validation |
| `norm.hpp` | COMSOL's field-and-DOF weighted RMS |
| `stepper.hpp` | The solve request, the run, and its output |
| `second_order.hpp` | The pair `(u, u')` and the condensed coefficients |
| `detail/state.hpp` | The state concept and its `axpy` arithmetic |
| `detail/controller.hpp` | Scalar error estimates, step size and order |
| `detail/output.hpp` | Step placement against the output times |

## Ownership and the solve boundary

`TimeStepper` owns its settings, solve, norm and the state of a run:

| Member | What it is |
| :--- | :--- |
| `accepted_` | the accepted state, derivative, time, order and difference columns |
| `controller_` | the step size and order proposed for the next step |
| `work_` | five buffers reused by successive attempts |
| `output_` | a pointer to the last state returned, and one interpolated buffer |

The accepted history, controller and output reference persist across steps; the
work buffers do not. `output_` points at `accepted_.state` for a level and at
its own buffer for an interpolated value, so a stored level costs no copy.

A step's coefficients — size, order, the table factors, the distances back to
the levels, IDA's `sigma` and the error estimates — are a `Trial`, created and
consumed inside one pass of `take_step`. A refused step is a `continue`; the
caller sees only accepted steps.

A `StepRequest` is constructed only by the stepper and borrows its buffers for
one solve callback. The formulation is the fixed-leading BDF history of a
first-order equation:

```text
y' = cj*y + known
M y' + K y = f  =>  (cj*M + K)*y = f - M*known
```

The solve starts from the predictor and returns convergence and iterations;
matrix assembly and algebraic solvers belong to the caller. Custom states
implement `axpy`; custom norms implement `norm(error, reference)`.

## Start and consistent initialization

`start(t0, initial, end)` performs two artificial backward-Euler solves of size
`backward_euler_fraction` of the initial step, the COMSOL default. The resulting
level and derivative are the value `run` stores at `t0`; the probes advance the
internal time to `t0 + 2*probe`, which `time()` reports. The initialization
solves are not accepted time steps. The initial step is
`initial_step_fraction * max(abs(end), span)` (the manual step in manual mode),
limited by `0.5 / norm(y')`.

`start(t0, initial, derivative, end)` takes consistent initial values from the
caller and performs no artificial step. It is the restart path: the initial
level is not perturbed and no solve is spent before the first step. The same
history setup follows: `phi[0] = initial`, `phi[1] = h*derivative`.

## Accepted history

Integration stores scaled divided-difference columns and step distances. An
attempt reads the accepted table and writes nothing; a refusal leaves it exactly
as the accepted steps left it. Accepting rescales the columns in place, writes
the correction and folds the columns down. The lead coefficient is `H_k/h`,
where `H_k = 1 + 1/2 + ... + 1/k`.

The first step's column is the derivative the initialization measured, at the
size the first step proposes. Older solution levels are not stored. `phi[0]` is
updated by accumulation and can differ from the solved state in its last bits;
both are kept.

The kernel holds fifteen state objects (five work buffers, seven difference
columns, the accepted state, derivative and output reference) whatever the
order. A second-order run is a `SecondOrderState`, whose two halves make that
thirty; there is no adapter layer.

## Second-order equations

`SecondOrderState<State>` is `(u, v = u')`; the kernel integrates it as the
first-order system

```text
u' - v            = 0
M v' + C v + K u  = f
```

With `cj` the step's first-derivative weight:

```text
u'  = cj u + known1
u'' = cj^2 u + known2,   known2 = cj known1 + known(v)
(K + cj C + cj^2 M) u = f - M known2 - C known1
```

`second_derivative_history` builds `known2` into a caller buffer;
`set_velocity` enforces `v = cj u + known1` after a successful solve. The pair's
`size()` and `operator[]` expose the displacement, so the default norm measures
the displacement alone. A custom norm can read both halves. The COMSOL rule
that includes the auxiliary velocity when consistent initialization is on, and
the derivative tolerance that weights it, are not implemented; no COMSOL
second-order model has been recorded.

Differencing displacement levels twice is a different and worse discretization:
the BDF polynomial through the levels has a parasitic root of multiplicity two
and the run is first order. The pair does not do that.

## Error weights

The default norm is the documented COMSOL field-and-DOF weighted RMS:

```text
sqrt( (1/M) * sum_j (1/N_j) * sum_i (E_i / (A_i + R*abs(Y_i)))^2 )
```

`field_sizes` gives `N_j`; empty is one field. With the scaled method, a field's
scale reference `S` gives

```text
norm2 = sqrt(sum S_i^2),  norm_inf = max abs(S_i)
d = exp(-j/5),  j = accepted steps
beta = (1-d) norm_inf + d   if 0 < norm_inf < 1, else 1
A_i = absolute * (beta + norm2)
```

and `A_i = absolute` unscaled. The relative part uses the accepted internal
state; the scaled part uses the last state returned to the caller, including an
interpolated one.

A norm of a caller's own type receives the error and the accepted internal
state; a `WeightedRmsNorm` passed explicitly also receives the scale reference.
A caller that needs a per-variable scaling or another reference implements it
in the norm.

## Output and step placement

`StepSchedule` holds the output times and the two modes. An output time is the
start time or after the initialized time, and the list increases; the initial
state is stored if the list starts at the start time, and always in `steps`
mode. `strict` lands on each output time; `intermediate` ensures a step in each
interval; `free` and `manual` let the step run.

A stored level is the accepted state; a value between levels is the BDF dense
polynomial, whose derivative at the newest node is the solve derivative. The
`interpolated` flag is false when the step ended on the output time. The dense
output writes into a state of the caller's; the same buffer can serve every
call, and times older than the last step are extrapolated.

## Validation boundaries

Invalid dimensions, non-finite configuration and output times, and a schedule
outside the run are rejected. A nonzero error with zero weight has infinite
norm. A step must be finite, positive and advance floating-point time. State
operations must not throw; a throwing one leaves the run in an unspecified
state, and `start` must rebuild it.

The [COMSOL regressions](comsol-alignment.md) cover 21 scalar configurations and
1,336 samples, plus the 161-sample thermoelastic run. Reference data never
enters a solve.
