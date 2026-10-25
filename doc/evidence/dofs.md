# COMSOL's error test over variables of several degrees of freedom

What the local COMSOL 6.2.0.290 solver does when a model's dependent variables
carry more than one degree of freedom, recorded by
`tools/comsol/GeDofSweep.java` and `tools/comsol/run_dofs.py`, and stored as raw
doubles in `test/reference/comsol_dofs.hpp`.

Every earlier recording is a model whose dependent variables each hold one
degree of freedom — the scalar decay, the four tolerance models, the two-field
sweep's global equation — or a field whose dofs all carry the same value. None of
them can see how the norm treats a variable of many dofs, whether `M` counts
variables or physics interfaces, or whether a multi-variable run is replayable.
This sweep is the same model under six settings, with the algebraic solve made
exact so that the trajectory is a property of the time integration alone.

- Regenerate with `python tools/comsol/run_dofs.py --directory .cache/comsol-dofs`,
  or rewrite the header from a recorded run with `--header-only`. The whole sweep
  is under a minute including the Java compilation and the COMSOL start; raw logs
  and exports stay in `.cache/comsol-dofs/dofs/`.
- Six runs, 169 stored samples recorded and 110 kept: the three controls are
  compared against the kept runs by the generator, so storing them again would
  be a second copy of the same numbers.

## The model as built

One component, one 1D geometry, a Coefficient Form PDE and up to two Global
Equations:

```text
c   Coefficient Form PDE on the default unit interval, five linear elements,
    da = 1, a = 1, c = 0, f = 0, u(x,0) = 1+x
    -> elements + 1 = 6 solved dofs, and two internal reaction-force dofs the
       solver reports apart from them
ge  Global Equations, one degree of freedom each:
    vt + 0.5*v = 0, v(0) = 0.01     and    wt + 0.1*w = 0, w(0) = 0.1
```

`da = a = 1` with no spatial operator gives every field dof `u_i' = -u_i`, so the
exported solution has an analytic reference, `u_i(t) = (1+x_i) exp(-(t-t0))`, and
the initial value is deliberately not uniform: the six dofs hold 1, 1.2 … 2, so a
per-dof weight, a per-variable average and the maximum of the scaled solution are
three different numbers.

Every run integrates `[0,10]` with `tstepsbdf = free`, `tout = tsteps`,
`tstepsstore = 1`, `keeplog = on`, `maxorder = 5`, `rtol = 0.01`, the solver's
default absolute-tolerance factor (0.1, so `A_k = 1e-3`), every accepted step in
the log (`logsampling = 0`), and `atolglobalmethod = unscaled` except where the
scale's role is the question. The variable scales are pinned to one
(`scalemethod = manual`, `scaleval = 1`) except in the one run that asks for the
solver's own.

Each run also sets the Fully Coupled node's `jtech = onevery` and `ntolfact =
1e-8`, so that every Newton iteration builds its Jacobian. That is what makes the
trajectory replayable: with COMSOL's default Jacobian reuse the linear solve
terminates on a small increment instead of a small residual, and the recorded
state stops being the exact solve's after the first step-size change
(evidence/tolerance-axes.md measures that on the seven-dof model).

| Case                   | variables (dofs)  | criterion | scale | stored steps | rejected | stored samples |
| ---------------------- | ----------------- | --------- | ----- | ------------ | -------- | -------------- |
| `field_only`           | `u` (6)           | unscaled  | 1     | 30           | 2        | 31             |
| `field_only_scaled`    | `u` (6)           | scaled    | 1     | 30           | 2        | 31             |
| `field_only_auto`      | `u` (6)           | scaled    | 2     | 43           | 5        | 44             |
| `field_global2_one`    | `u` (6), `v`, `w` | unscaled  | 1     | 27           | 1        | 28             |
| `field_global2_two`    | `u` (6), `v`, `w` | unscaled  | 1     | 27           | 1        | 28             |
| `field_closest_uneven` | `u` (6)           | unscaled  | 1     | 30           | 2        | 7              |

## 1. `M` counts dependent variables, not physics interfaces

`field_global2_one` puts both global equations in one Global Equations interface;
`field_global2_two` gives each its own. The dependent variables, their equations,
their size, the tolerance, the scale, the mesh and the solver settings are
identical, so a norm that averages over physics interfaces predicts different
step sizes and a norm that averages over dependent variables predicts the same
run:

| Pair                          | samples  | worst difference over all times and dofs |
| ----------------------------- | -------- | ---------------------------------------- |
| `field_global2_one` vs `_two` | 28 vs 28 | 0 (byte-identical)                       |
| `field_only` vs `_scaled`     | 31 vs 31 | 0 (byte-identical)                       |

Zero, not small: the two models share one trajectory to the last bit. So `M = 3`
here, the field and the two global equations, and an interface is not a unit of
the average. That is also what the solver's own log says: it prints one scale per
dependent variable, `comp1.u`, `comp1.ODE1`, `comp1.ODE2`, and the count of
degrees of freedom it reports (`8 (plus 2 internal DOFs)`) is the sum over them.

The second row is the control for the scaled criterion: with a given scale its
absolute tolerance stays `A_k`, so `scaled` and `unscaled` are the same
criterion, and the two records are the same run.

Which partitions a model can expose at all is worth stating. With every variable
of one dof, averaging over variables and averaging over the dofs of one interface
are the same arithmetic, so a model of scalar variables cannot tell the readings
apart whatever it does. `field_global2_*` has a field of six dofs beside two
one-dof variables, so the two readings do differ — `(1/3)(r_u^2 + r_v^2 + r_w^2)`
against `(1/2)(r_u^2 + (r_v^2 + r_w^2)/2)` — and the two records are identical,
which is what rules the interface reading out.

## 2. A variable of six degrees of freedom replays over its whole run

Replayed from the model's own initial state, with the two backward-Euler probes
of the consistent initialization and an exact solve of `y' = -rate*y`, the two
unscaled records reproduce every decision and every number:

| Case                | steps   | rejections | samples | worst `abs(dt)` | worst `abs(du)` |
| ------------------- | ------- | ---------- | ------- | --------------- | --------------- |
| `field_only`        | 30 / 30 | 2 / 2      | 31 / 31 | 2.2e-12         | 5.3e-13         |
| `field_global2_one` | 27 / 27 | 1 / 1      | 28 / 28 | 1.2e-12         | 3.3e-13         |

The field's six dofs carry six different error-to-weight ratios, so the
per-variable average is exercised by every step, and the three-variable case puts
that average beside two variables of one dof. Twenty-seven and thirty steps of
agreement is not a first-step calibration: the multi-dof average, the per-variable
partition and the retry after a rejection all have to be right, and the residual
is the arithmetic of the two initialization probes.

## 3. The automatic scale's `beta` takes the maximum, not the RMS

`field_only_auto` leaves the scale to the solver, which prints it as `2` — the
maximum of the initial field, exactly representable and printed in full. With
`scaled` and the solution-following absolute tolerance, the weight of a dof is
`A_k(beta + rms(Y)) + R abs(Y_i)` with `Y = U/2`, and the field decays below 1 from
the first step, so `beta` is in its non-trivial branch for the whole run.

The two readings of the argument of `beta` — the maximum of `abs(Y)` over the
variable's dofs, and their RMS — differ by 29% at `t = 0` and less later, and the
44-sample record separates them:

| `beta`'s argument           | steps   | samples | worst difference from the record                            |
| --------------------------- | ------- | ------- | ----------------------------------------------------------- |
| `max(abs(Y))` over the dofs | 43 / 43 | 44 / 44 | 1.1e-12 in time, 4.6e-13 in value                           |
| `rms(Y)` over the dofs      | 43 / 43 | 44 / 44 | agrees for 11 samples, then leaves it by 1.9e-3 at the 12th |

The step count is the same under both readings and the trajectory is not, which
is why the test compares every stored value and not the decisions alone: the
growth-limited first steps hide the norm until the error test binds, and from
there the two readings part.

## 4. The closest step stores one value per solver step

`field_closest_uneven` asks for the closest-step storage on the output schedule
`0, 0.013, 0.5, 0.53, 0.56, 2, 2.05, 7, 10`. Nine requested times, seven stored
samples, and the run takes the same 30 steps as `field_only`: the schedule does
not move the steps, and every stored value is one of the run's own levels, bit
for bit — no interpolation except at the end time, which a step passed.

Reading the record against the step ends, with `n_k = min {n : t_n >= tau_k}`:

| requested time | first step at or after it | stored |
| -------------- | ------------------------- | ------ |
| 0              | the initial state         | 0      |
| 0.013          | 0.02142136                | yes    |
| 0.5            | 0.51365274                | yes    |
| 0.53           | 0.68486365                | yes    |
| 0.56           | 0.68486365                | **no** |
| 2              | 2.14918126                | yes    |
| 2.05           | 2.14918126                | **no** |
| 7              | 7.26294398                | yes    |
| 10             | the end, interpolated     | yes    |

So the rule has two halves, and the second one is invisible on a uniform
schedule, which is the only schedule the earlier sweeps recorded:

1. The stored time is the first solver step that reaches the requested time, not
   the nearest step end. At `0.5` the previous step ends at `0.34244183`, which is
   0.158 s away, and the next at `0.51365274`, which is 0.0137 s away; the record
   stores the later one.
2. A solver step stores one value, however many requested times it reaches: the
   three outputs in `[0.5, 0.56]` collapse into the step at `0.68486365`, and the
   two in `[2, 2.05]` into the step at `2.14918126`.

The single-variable uniform `store_closest` case of the stepping sweep is
consistent with the first half and cannot see the second.

## What this settles, and what it does not

Settled for the local COMSOL 6.2.0.290, by the records above:

- The norm's `M` is a count of dependent variables, and the average is over the
  dofs of each variable (§1).
- A run of one field of six dofs and of three variables with seven and eight dofs
  replays over its whole length, decisions and values (§2) — which answers the
  question the tolerance sweep had to leave open, because its own records kept
  COMSOL's inexact Newton.
- `beta` takes the maximum of the scaled solution over the variable's dofs (§3).
- The closest-step storage stores the first step reaching each requested time,
  once per step (§4).

Not settled:

- Q07, the second-order divergence at its nineteenth accepted step: nothing here
  involves the pair.
- What a run with COMSOL's default Jacobian reuse does. This sweep removes the
  inexact solve deliberately; reproducing it is the caller's business, and the
  test-only `ComsolDecaySolve` covers the scalar case only.
- One dimension, linear elements, decoupled linear equations. The partition
  question is answered where the partition is observable, which needs variables
  of unequal dof counts in one interface; every COMSOL interface that holds
  several variables holds ones of equal dofs (the global equations, the fields of
  a vector unknown), so where it matters the two readings coincide anyway.
- An automatic scale that the log prints to fewer digits than it holds. The one
  measured here is `2`, and the replay depends on it being exact.
- `A_k`: the sweep left COMSOL's default factor at `0.1` and never varied it, so
  the record fixes the weight but not the axis.
