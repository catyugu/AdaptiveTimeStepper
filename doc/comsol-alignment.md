# COMSOL BDF alignment

The first-order BDF implementation is checked against COMSOL 6.2.0.290 on
`u' + u = 0`, `u(0) = 1`, over `[0,10]`. All 21 recorded configurations agree
on accepted-step counts, rejection counts, stored times, values, derivatives
and every order printed in COMSOL's solver log. The tests compare trajectories,
not a band around the step count.

## Evidence and acceptance

`tools/comsol/GeDecayTransient.java` builds and solves the COMSOL models. It
reads time, state and derivative through `getPVals`, `getU` and `getUDot`.
`tools/comsol/run_decay.py` validates the exports against `exp(-t)` and the
solver log, then writes `test/reference/comsol_trajectories.hpp`. The 1,336
samples come from COMSOL, not from this library.

The regression checks each stored time, value and derivative, the accepted,
error-test-failure and failed-solve counts, and the order of every logged step.
One tolerance covers every recorded quantity: `1e-8`. The measured worst case
is about `5e-10` in time and `4e-10` in the derivative, so the number keeps an
order of magnitude of margin.

| Configuration | COMSOL accepted steps |
| :--- | ---: |
| Global-equation baseline, maximum order 5 | 37 |
| Maximum orders 1, 2, 3, 4 | 108, 58, 39, 39 |
| Relative tolerances 1e-3, 1e-5 | 60, 98 |
| Strict, intermediate | 49, 37 |
| Maximum step 0.05 s | 204 |
| Land exactly at the end | 37 |
| Exact algebraic solve, default / tight / strict | 44, 101, 45 |
| Nonuniform output times, strict / intermediate | 47, 42 |
| Manual step 0.1 s, default / exact algebraic solve | 100, 100 |
| Manual step 0.2 s | 50 |
| Closest-step / interpolated output | 37, 35 |

The nonuniform schedule is `0 .001 .01 .05 .2 .3 .7 1 1.1 1.4 2 3 5 8 10`.

## Default order depends on the physics interface

The generic Time solver's API default is `maxorder = 5`; physics-controlled
solver generation can select a lower maximum. A local COMSOL 6.2.0.290 readback
on 2026-10-23, without setting `maxorder`, `bdforder` or `timemethod`, gave:

| Configuration | Method | Maximum BDF order | Manual BDF order |
| :--- | :--- | ---: | ---: |
| Directly created Time solver node | BDF | 5 | 2 |
| Global Equations, automatically generated solver | BDF | 5 | 2 |
| Heat Transfer, automatically generated solver | BDF | 2 | 2 |

`Settings::max_order = 5` follows the generic solver. Set it to 2 for the
automatic Heat Transfer or Navier–Stokes configuration. This limits adaptive
orders to 1–2; it does not force every step to order 2.

## Fixed-leading BDF history

COMSOL's state/derivative pairs distinguish its method from a classical BDF
stencil reconstructed from recent time levels. The first-derivative formula is

```text
cj = (1 + 1/2 + ... + 1/k) / h
y' = predicted_derivative + cj * (y - predicted_state).
```

The predictor and error estimates use persistent scaled divided differences.
Each attempt scales every column it reads and writes none; a rejection leaves
the accepted table as it is. Reconstructing a lower-degree Lagrange polynomial
is not equivalent after a step or order change.

## Initial step and maximum step

COMSOL 6.4 states the automatic BDF first step as `min(0.1% of the end time,
DT)` with `DT * norm(y') <= 1/2`. The library uses

```text
initial_step_fraction * max(abs(end time), span)
```

so a run that starts at zero matches COMSOL's `0.1% of the end time`; a run that
starts later uses the final time rather than the shorter span. The recorded runs
all start at zero, so the later-start rule is not measured against COMSOL.

COMSOL's automatic maximum step is **not documented**. The 10% figure in its
manuals belongs to the Particle Tracing ODE solver, not the FEM Time-Dependent
Solver. The library therefore has no automatic maximum: `max_step = 0` is no
limit, and `max_step` is an absolute time. The recorded `max_step_5e3` case sets
COMSOL's `maxstepbdf = 0.05` and the library's `max_step = 0.05`; the other
recorded runs leave the maximum unset and their largest step is below it, so
they do not measure a default.

## Tolerances and the error norm

COMSOL's error test is a field-and-DOF weighted RMS:

```text
sqrt( (1/M) * sum_j (1/N_j) * sum_i (E_i / (A_i + R*abs(Y_i)))^2 )
```

with `M` fields and `N_j` dofs in field `j`. `WeightedRmsNorm` implements it;
`Settings::field_sizes` supplies `N_j`. With the scaled method, a field of scale
reference `S` gives

```text
norm2 = sqrt(sum S_i^2),  norm_inf = max abs(S_i)
d = exp(-j/5)
beta = (1-d) norm_inf + d   if 0 < norm_inf < 1, else 1
A_i = R*factor * (beta + norm2)
```

The 21 recorded runs are a single scalar equation: one field, one dof, and
`abs(u) <= 1`. In that case the formula reduces to the scalar
`R*factor*((2-d) abs(u) + d)` that the recorded runs measured. The field-aware
formula for multiple fields or multi-dof fields is the documented one, but no
COMSOL run with several fields was recorded here.

COMSOL's per-variable scaling and its derivative tolerance are not implemented.
A caller that needs them supplies a custom norm; the default norm uses the
documented formula only.

COMSOL also updates the magnitude inside `A_i` from the solution returned to the
caller, including an interpolated solution, while its relative part stays the
internal state. The library does the same: `output_` points at the last stored
value.

## Initialization and mode constraints

The recorded consistent initialization takes two short backward-Euler solves
and labels the resulting state with the requested start time. Automatic runs
use `initial_step_fraction * backward_euler_fraction` of the initial-step base;
manual runs use `manual_step * initial_step_fraction * backward_euler_fraction`.
The initialization solves are not counted as accepted time steps.

`start(t0, initial, derivative, end)` is the COMSOL `Consistent initialization:
Off` analog: it takes the initial values as consistent, uses the supplied
derivative, and spends no artificial solve. COMSOL's `On` (index-1) routine,
`Safety factor for Backward Euler`, and `Rescale after initialization` are not
implemented.

Strict stepping and end-time landing avoid a small remainder by splitting the
last two steps. Intermediate mode counts an accepted step that crossed an
output time as a landing in the next interval. The first shortened step becomes
the base for startup growth.

COMSOL's `tstepsclosest` stores the first accepted step reaching an output time.
At an internal node the stored derivative is the derivative from the solve;
between nodes it is the dense-output polynomial's derivative. `steps` stores
the solver steps; when the last step passes the end and end interpolation is on,
the end is stored by interpolation instead of that step.

## Second-order equations

COMSOL's IDA path introduces `Uv = u'` and condenses the linearized system to
`E + sigma D + sigma^2 K`. The library exposes the pair `(u, v)` and the
condensed coefficients; the caller assembles and solves the displacement
system.

COMSOL's error test **excludes** `Uv` only while consistent initialization is
off; with it on the differential `Uv` variables are included, with a derivative
tolerance. The library's pair exposes the displacement to the default norm, so
it measures the displacement alone; a custom norm can include the velocity, but
the derivative tolerance is not implemented. No COMSOL second-order model has
been recorded, so the pair is not claimed as COMSOL-validated.

Structural dynamics in COMSOL uses generalized alpha by default, not BDF; the
library does not provide generalized alpha.

## What is not implemented

COMSOL's algebraic-variable error-test options (`Include`/`Exclude algebraic`,
singular mass detection), per-variable tolerances, the derivative tolerance,
`Store every Nth step`, the nonlinear controller (STAB), the time-step-increase
delay, events, generalized alpha and Runge–Kutta are out of scope. The library
has no automatic maximum step.

## Reproduce

```text
python tools/comsol/run_decay.py --directory .cache/comsol-reference
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

COMSOL is needed only to regenerate the reference. Put `comsolcompile` and
`comsolbatch` on `PATH`, or set `COMSOL_ROOT`. Compilation and batch logs are
checked because these launchers can report errors while returning exit code 0.

Primary references: [COMSOL time-dependent solver API](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_solver.51.50.html),
[COMSOL time-step control](https://www.comsol.com/support/knowledgebase/1254),
and [SUNDIALS IDA source](https://github.com/LLNL/sundials/tree/main/src/ida).
The measurements above are from the locally executed COMSOL 6.2 model.
