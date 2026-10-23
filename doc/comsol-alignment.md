# COMSOL BDF alignment

The first-order BDF implementation is checked against COMSOL 6.2.0.290 on
`u' + u = 0`, `u(0) = 1`, over `[0,10]`. All 21 recorded configurations agree
on accepted-step counts, rejection counts, stored times, values, derivatives,
and every order printed in COMSOL's solver log. The tests compare actual
trajectories, rather than accepting a band around the step count.

## Evidence and acceptance

`tools/comsol/GeDecayTransient.java` builds and solves the COMSOL models.
It reads time, state and derivative directly through `getPVals`, `getU` and
`getUDot`. Data-export time labels were unsuitable: they round the time, and
the final label `10` may be an interpolated value from a step ending past 10.
`tools/comsol/run.py` validates the exports against `exp(-t)`
and the solver log, then writes `test/reference/comsol_trajectories.hpp`.
The 1,336 samples in that file come from COMSOL, not from this library.

The regression checks each stored time, value and derivative, the accepted,
error-test-failure and failed-solve counts, and the order of every logged step.
One tolerance covers every recorded quantity: `1e-8`. The measured worst case
over all 21 configurations is `5e-10` in time and `7e-10` in the derivative — the
latter is a consistent-initialization derivative, reconstructed by dividing by a
backward-Euler step as small as `1e-7 s` — so the single number keeps an order of
magnitude of margin without a special case per quantity.

| Configuration                                      | COMSOL accepted steps |
| -------------------------------------------------- | :-------------------- |
| Default, maximum order 5                           | 37                    |
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
Each attempt rescales a copy of the accepted history; a rejection discards the
copy. Lower-order estimates use the corresponding difference columns plus the
correction. Reconstructing a lower-degree Lagrange polynomial is not equivalent
after a step or order change. This is a change to the numerical formula, not a
constant fitted to the reference step counts.

Use `request.derivative_weight()`, `request.derivative_history()` and
`request.derivative(y)` when assembling first-order systems. The older
`request.formula->weights(1)` interface describes the classical stencil and
cannot express the confirmed fixed-leading history.

### Tolerances and output

The weight of the error test is per degree of freedom, and the scaled absolute
tolerance is

```text
A_i = rtol * factor * ( |y_i| + s_i ),
s_i = (1 - exp(-j/5)) |y_i| + exp(-j/5) * sigma_i,
```

with `j` the accepted-step count and `sigma_i` a scale of the variable — the
"Scales for dependent variables" line of the solver log. Three measurements fix
it, and each of them is a model the recorded reference cannot see:

| Measurement                                            | Result                                                                                                                                                      |
| ------------------------------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `u' = -1`, `u(0) = c`, four magnitudes                 | The first step is `0.5 (A + rtol \|y\|)` with `A = 2 rtol factor c` to five digits. The first term is `\|y_i\|`; a single dof's scale is its own magnitude. |
| `u' + u = 0`, the recorded 44-step exact run           | A constant `sigma = 1` gives 31 steps; `s_i = \|y_i\|` gives 45. The `exp(-j/5)` blend gives 44.                                                            |
| `u' = -1e-9 u, v' = -1` and `u' = -u + v, v' = -100 v` | The first step needs `sigma = (1, 1)` and `(1, 0.32)`, respectively: scale is per variable and comes from the model.                                        |

An earlier form substituted the state's 2-norm for `|y_i|` and the constant one
for `sigma_i`. It agrees with the measured form only when `sigma_i = 1`, which is
every configuration of the reference — the 21 recorded runs are one scalar
equation, so the difference was invisible to them. Replaying the two models whose
scale is not one, with the scale passed in, reproduces COMSOL's accepted steps,
every logged order and its error-test-failure count exactly.

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
use `span * first_step_fraction * init_step_fraction` for each solve; manual
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
time; it can choose a later point even when the previous point is nearer.
At an internal node the stored derivative is the derivative from the solve.
Between nodes it is the dense-output polynomial's derivative.

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
Coupled-field weighting, nonlinear physics, events and other COMSOL versions
need their own reference runs. Direct second derivatives use the separately
documented classical formula; first-order systems use the validated BDF path.

## Reproduce

```text
python tools/comsol/run.py --directory .cache/comsol-reference
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
