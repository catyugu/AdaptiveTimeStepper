# Variables of several degrees of freedom

What COMSOL 6.2 does when a dependent variable carries more than one degree of
freedom, recorded by `tools/comsol/GeDofSweep.java` into
`test/reference/comsol_dofs.hpp` (4 runs kept of 6). The model is one
Coefficient Form PDE and up to two Global Equations:

```text
c   u(x,0) = 1 + x, da = a = 1, c = f = 0   ->  six solved dofs 1, 1.2, ... 2,
                                                every one with u_i' = -u_i
ge  v' + 0.5 v = 0, v(0) = 0.01  and  w' + 0.1 w = 0, w(0) = 0.1   one dof each
```

The non-uniform initial value is what makes the questions decidable: a per-dof
weight, a per-variable average and the maximum of the scaled solution are three
different numbers. Every run pins the scales to one and takes the unscaled
criterion except where the scale is the question, and every run rebuilds the
Jacobian at every iteration (`jtech = onevery`, `ntolfact = 1e-8`) so that its
trajectory is the time integration's alone — with COMSOL's default reuse a
multi-variable record is not the exact solve's after its first step-size change.

## `M` counts dependent variables, not physics interfaces

Both global equations in one interface and the same two in an interface each are
one trajectory to the last bit, with the variables, their sizes, the mesh and
every setting identical. A partition of unequal size then matters: replaying the
three-variable record with `(6,1,1)` gives 27 accepted steps and with `(7,1)`,
which is the interface reading, 28.

## A variable of six dofs replays over its whole run

| Case                | steps   | rejections | worst difference                  |
| ------------------- | ------- | ---------- | --------------------------------- |
| `field_only`        | 30 / 30 | 2 / 2      | 2.2e-12 in time, 5.3e-13 in value |
| `field_global2_one` | 27 / 27 | 1 / 1      | 1.2e-12 in time, 3.3e-13 in value |

Six different error-to-weight ratios in one variable, three variables in the
other case: the per-variable average, the partition and the retry after a
rejection all have to be right, and what is left is the arithmetic of the two
initialization probes.

## `beta`'s argument is the maximum, not the RMS

With the solver's own scale (printed `2`, the maximum of the initial field) the
field decays below 1 from the first step, so `beta` is in its non-trivial branch
for the whole run and the two readings differ by 29% at `t = 0`. The 44-sample
record reproduces with the maximum; with the RMS it agrees for eleven samples and
then leaves the record by 1.9e-3.

## The closest-step storage stores one value per step

On the schedule `0, 0.013, 0.5, 0.53, 0.56, 2, 2.05, 7, 10` the run stores seven
samples: the first step at or past each requested time, with the three outputs in
`[0.5, 0.56]` collapsing into the step at 0.68486365 and two in `[2, 2.05]` into
2.14918126. Every stored value is one of the run's own levels, bit for bit, and
the end time is interpolated. At `0.5` the step before ends at 0.3424, 0.158 away,
and the one after at 0.5137, 0.0137 away: the later one is stored, so the rule is
crossing and not nearest.

## The weight at a near-zero magnitude

A run started at the divergence's state (`a = 92.9388`, `b = 0.2801`, manual
scales) has the automatic first step `h1 = 2.8769e-5`. That step is
`0.5 / norm(y')` with `y' = (b, -a - 0.1b)`, so it reads the near-zero variable's
weight out as `3.78e-3`, against `A_k + R abs(b) = 3.80e-3`; a weight of `A_k`
alone would give `7.6e-6`.

## What this does not settle

- One dimension, linear elements, decoupled linear equations. An interface that
  holds several variables of unequal size would be needed to separate the
  partition from the interface reading; every COMSOL interface that holds several
  variables holds ones of equal dofs, so where it matters the readings coincide.
- An automatic scale the log prints to fewer digits than it holds: the replay
  needs it exact.
- `A_k`: every run leaves COMSOL's default factor at 0.1, so the record fixes the
  weight and not that axis.
- COMSOL's default Jacobian reuse, which the replay removes deliberately: the
  algebraic solve is the caller's, and a record taken with it is reproducible only
  by reproducing COMSOL's iteration.
