# COMSOL BDF alignment

The first-order BDF implementation is checked against COMSOL 6.2.0.290 on
`u' + u = 0`, `u(0) = 1`, over `[0,10]`. All 21 recorded configurations agree
on accepted-step counts, rejection counts, stored times, values and
derivatives. The tests compare actual trajectories, rather than accepting a
band around the step count. The orders in COMSOL's solver log are printed as
a diagnostic: the log is sampled on the wall clock, so its rows are not
reproducible (doc/evidence/stepping-sweep.md, "The log is sampled unless told
otherwise").

## Evidence and acceptance

`tools/comsol/GeSteppingSweep.java` builds and solves the COMSOL models.
It reads time, state and derivative directly through `getPVals`, `getU` and
`getUDot`. Data-export time labels were unsuitable: they round the time, and
the final label `10` may be an interpolated value from a step ending past 10.
`tools/comsol/run_stepping.py` validates the exports against `exp(-t)` and the
solver log, then writes `test/reference/comsol_stepping.hpp`; `run_tolerance.py`,
`run_twofield.py`, `run_thermoelastic.py`, `run_second_order.py` and
`run_dofs.py` do the same for the other five models.
The 1,231 samples in that file come from COMSOL, not from this library.

The regression checks each stored time, value and derivative, the accepted,
error-test-failure and failed-solve counts, and the order of every logged step.
One tolerance covers every recorded quantity: `1e-8`. The measured worst case
over all 21 configurations is `5e-10` in time and `7e-10` in the derivative — the
latter is a consistent-initialization derivative, reconstructed by dividing by a
backward-Euler step as small as `1e-7 s` — so the single number keeps an order of
magnitude of margin without a special case per quantity.

| Configuration                                      | COMSOL accepted steps |
| -------------------------------------------------- | :-------------------- |
| Global-equation baseline, maximum order 5          | 37                    |
| Maximum orders 1, 2, 3, 4                          | 108, 58, 39, 39       |
| Relative tolerances 1e-3, 1e-5                     | 60, 98                |
| Strict, intermediate                               | 49, 37                |
| Maximum step 0.05 s                                | 204                   |
| Land exactly at the end                            | 37                    |
| Exact algebraic solve, default / tight / strict    | 44, 101, 45           |
| Nonuniform output times, strict / intermediate     | 47, 42                |
| Manual step 0.1 s, default / exact algebraic solve | 100, 100              |
| Manual step 0.2 s                                  | 50                    |
| Closest-step / interpolated output                 | 37, 35                |

The nonuniform schedule is `0 .001 .01 .05 .2 .3 .7 1 1.1 1.4 2 3 5 8 10`.
It forces the mode constraints to act; checking only equally spaced outputs
with a smaller maximum step would not exercise intermediate stepping.

## Default order depends on the physics interface

There is no single maximum-order default for all COMSOL transient models.
The generic Time solver's API default is `maxorder = 5`, but physics-controlled
solver generation can select a lower maximum.

A local COMSOL 6.2.0.290 readback on 2026-10-23, without setting `maxorder`,
`bdforder` or `timemethod`, gave:

| Configuration                                    | Method | Maximum BDF order (`maxorder`) | Manual BDF order (`bdforder`) |
| ------------------------------------------------ | ------ | ------------------------------ | ----------------------------- |
| Directly created Time solver node                | BDF    | 5                              | 2                             |
| Global Equations, automatically generated solver | BDF    | 5                              | 2                             |
| Heat Transfer, automatically generated solver    | BDF    | 2                              | 2                             |

The [COMSOL 6.2 Time API](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/comsol_api_solver.49.48.html)
lists `maxorder = 5` and the separate manual-mode `bdforder = 2`.
The [Heat Transfer documentation](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_heattransfer.30.65.html)
specifies second-order BDF as its physics default. The
[Navier–Stokes solver documentation](https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_single.06.103.html)
also specifies a maximum order of two; that interface was not included in
the local readback above.

`Settings::max_order = 5` follows the generic solver default. Set it to 2 to
match the automatic Heat Transfer or Navier–Stokes configuration. This limits
adaptive orders to 1–2; it does not force every step to use order 2.

The recorded case named `default` uses Global Equations and explicitly sets
`maxorder = 5` in `GeSteppingSweep.java`. It is a configured regression
baseline, not evidence that all COMSOL physics interfaces default to 5.

## What caused the mismatch

### Fixed-leading BDF history

COMSOL's state/derivative pairs distinguish its method from a classical BDF
stencil reconstructed from recent time levels. The confirmed first-derivative
formula is

```text
cj = (1 + 1/2 + ... + 1/k) / h
y' = predicted_derivative + cj * (y - predicted_state).
```

The predictor and error estimates use persistent scaled divided differences.
Each attempt scales every column it reads by the factor its step size gives it,
and writes none of them; a rejection leaves the accepted table as it is.
Lower-order estimates use the corresponding difference columns plus the
correction. Reconstructing a lower-degree Lagrange polynomial is not equivalent
after a step or order change. This is a change to the numerical formula, not a
constant fitted to the reference step counts.

Use `request.derivative_weight()`, `request.derivative_history()` and
`request.derivative(y, dy)` when assembling first-order systems. A second time
derivative is not a weight of this formula: `SecondOrderStepper` hands the same
first-derivative weight to the pair `(u, u')` and lets the solve fill both halves,
which is what COMSOL's IDA path does. The classical stencil
version 1.0 offered for a second derivative is gone: it is a different
discretization, first order whatever the scheme's order.

### Tolerances and output

The recorded references are one scalar equation, and for a single dof the
weight of the error test is

```text
A_i = rtol * factor * ( |y_i| + s_i ),
s_i = (1 - exp(-j/5)) |y_i| + exp(-j/5) * sigma_i,
```

with `j` the accepted-step count and `sigma_i` the scale of the variable — the
"Scales for dependent variables" line of the solver log. Three measurements fix
that form, each of them a model the recorded reference cannot see:

| Measurement                                            | Result                                                                                                                                                      |
| ------------------------------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `u' = -1`, `u(0) = c`, four magnitudes                 | The first step is `0.5 (A + rtol \|y\|)` with `A = 2 rtol factor c` to five digits. The first term is `\|y_i\|`; a single dof's scale is its own magnitude. |
| `u' + u = 0`, the recorded 44-step exact run           | A constant `sigma = 1` gives 31 steps; `s_i = \|y_i\|` gives 45. The `exp(-j/5)` blend gives 44.                                                            |
| `u' = -1e-9 u, v' = -1` and `u' = -u + v, v' = -100 v` | The first step needs `sigma = (1, 1)` and `(1, 0.32)`, respectively: scale is per variable and comes from the model.                                        |

That form is COMSOL's general rule for a variable of one dof whose scale is
chosen automatically and updated. The general rule is per variable:

```text
Y_i   = U_i / sigma_variable
A_s,i = A_k (beta + rms(Y over the variable))   automatic scale, update on
A_s,i = A_k                                     otherwise
beta  = (1 - exp(-j/5)) max(abs(Y over the variable)) + exp(-j/5)
```

so the single-dof form is what its name says — a special case, not the rule.
The measured three-branch table, the scale's role as a divisor of the solution,
and the two recorded switches ("Update scaled absolute tolerance" off, and a
scale the model gives) are in [evidence/tolerance-axes.md](evidence/tolerance-axes.md).
The kernel implements the special case only: for a state of more than one
variable, or for a scale the model gives, its first step misses the record by up
to 9% before the repair; every recorded first step now reproduces (doc/evidence/tolerance-axes.md).

Two readings of that formula are decided on a model whose variables are not one
dof each ([evidence/dofs.md](evidence/dofs.md)). `M` counts dependent variables,
not physics interfaces: two global equations in one interface and the same two in
an interface each are one trajectory to the last bit, while a partition of
unequal size changes the step count. And `beta`'s argument is the maximum of
`abs(Y)` over the variable's dofs, not their RMS: the RMS reading agrees with the
44-sample record of a six-dof field for eleven samples and then leaves it by
1.9e-3. A variable of six dofs replays over its whole run once the algebraic
solve is exact, which is what the tolerance sweep's own records could not do.

COMSOL also updates the *magnitude inside* `A_i` from the solution returned to
the caller, including an interpolated solution. Its relative reference remains
the internal state. This was tested against three alternatives: unchanged
references, updating only the absolute one, and updating both. Only updating the
absolute one reproduced the 35-step interpolated-output run and its stored
values. Assuming that output cannot affect the internal step sequence would miss
this COMSOL behavior.

### Initialization and mode constraints

The recorded consistent initialization takes two short backward-Euler solves
and labels the resulting state with the requested start time. Automatic runs
use 0.1% of the span times `initialization_fraction` for each solve; manual
runs use the manual step instead of the span. The independent 0.2 s manual
case confirms that dependence. The initialization solves are not counted as
accepted time steps.

Strict stepping and end-time landing avoid a small remainder by splitting the
last two steps. When a proposed step overshoots a required landing time, the
clipped step includes the floating-point safety margin. A step already landing
on that time is retained. These details matter because changes of step size
reset the equal-step counter and can reset nonlinear convergence estimates.

Intermediate mode counts an accepted step that crossed an output time as a
landing in the next interval. It also prevents a jump across an entire empty
interval. The first shortened step becomes the base for startup growth.

COMSOL's `tstepsclosest` stores the first accepted step reaching an output
time; it can choose a later point even when the previous point is nearer. A step
stores one value however many output times it reaches, so several outputs inside
one step collapse into that step's time (evidence/dofs.md). At an internal node
the stored derivative is the derivative from the solve. Between nodes it is the
dense-output polynomial's derivative.

## Separate the algebraic solve from time stepping

COMSOL's default Fully Coupled node uses `jtech=minimal`, reusing a Jacobian.
Even this scalar linear equation can therefore have a nonzero residual at an
accepted step. An exact linear solve is not an equivalent comparison.

The test-only `ComsolDecaySolve` reproduces the scalar Jacobian-reuse and
iteration-termination behavior. It contains no reference trajectory values.
It remains outside the library, where the caller owns algebraic solution.
Independent COMSOL runs set `jtech=onevery` and `ntolfact=1e-8`; those runs are
compared with a separate exact-solve callback. Changing the COMSOL algebraic
strategy changes the trajectory, and both strategies are reproduced.

IDA source supplied hypotheses and a control-law reference; COMSOL experiments
decided whether those hypotheses held. The validation does not establish
bitwise equality for arbitrary models or undocumented COMSOL configurations.
Coupled-field weighting has its own reference run (the two-field sweep) and the
kernel implements it; the scaled tolerance's branch structure is measured in
doc/evidence/tolerance-axes.md and implemented.
The recorded references still cover one scalar equation, so a state of several
variables is validated by the acceptance cases of the tolerance sweep rather
than by a long trajectory. Nonlinear physics, events and other COMSOL versions
need their own reference runs. A second-order equation keeps the validated BDF path, through
the state augmentation COMSOL documents, rather than a separate discretization:
the pair `(u, u')` is the solve's unknown, and the error test measures the state
alone divided by the two variables the augmentation adds — the velocity is in the
system, not in the local error (doc/evidence/second-order.md).

## Reproduce

```text
python tools/comsol/run_stepping.py --directory .cache/comsol-reference
python tools/comsol/run_second_order.py --directory .cache/comsol-reference
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

COMSOL is required only to regenerate the reference. Normal tests use the
recorded data and have no COMSOL dependency. Put `comsolcompile` and
`comsolbatch` on PATH, or set `COMSOL_ROOT`. Compilation and batch logs are
checked because these launchers can report errors while returning exit code 0.

Primary references: [COMSOL time-dependent solver API](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_solver.51.50.html),
[COMSOL time-step control](https://www.comsol.com/support/knowledgebase/1254),
and [SUNDIALS IDA source](https://github.com/LLNL/sundials/tree/main/src/ida).
The measurements above are from the locally executed COMSOL 6.2 model.
