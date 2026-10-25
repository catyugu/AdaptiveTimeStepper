# The average over fields

Whether COMSOL's error test averages per field or over every degree of freedom,
recorded by `tools/comsol/GeFieldNormSweep.java` into
`test/reference/comsol_twofield.hpp` (3 cases kept of 11 recorded). The model is
one Global Equation and one Coefficient Form PDE field, decoupled:

```text
A: uA' + rateA uA = 0, uA(0) = 1        one dof
B: da u' + rateB u = 0, u(x,0) = 1      elements + 1 dofs, every one identical
```

Both scales are pinned to one, so the weights are the entered absolute tolerance
plus `rtol abs(y)`, and the two candidate means are

```text
per field:  E^2 = (1/M) sum_j (1/N_j) sum_{i in j} (e_i / W_i)^2
flat:       E^2 = (1/N) sum over every dof (e_i / W_i)^2
```

## The per-field average is the rule

The first step is a pure ratio of norms, calibrated against the same model with
A alone, so the weight cancels:

| Pair of rates          | measured `h1 / h1(A alone)` | per-field prediction | flat prediction |
| ---------------------- | --------------------------- | -------------------- | --------------- |
| `rateB/rateA` = 0.2236 | 1.380133024                 | 1.380133117          | 4.106           |
| `rateB/rateA` = 0.1    | 1.407195067                 | 1.407195089          | 7.124           |

The same rule survives a fourfold refinement of field B's mesh: 26 and 101 field
dofs give one trajectory to 4.6e-16, where a flat average would move every step
size by 34%. When the two rates are equal both formulas coincide, and the
combined run reproduces the single-field one to 2.1e-16 — so the insensitivity to
`N_B` is a property of the mean, not of the model.

## What this does not settle, and why it is not needed

- Whether `M` counts physics interfaces or dependent variables: this model has one
  variable per interface, so the two readings coincide. A model where they differ
  is in [dofs.md](dofs.md).
- The values of a two-field trajectory. These runs keep COMSOL's default Jacobian
  reuse, so their states are not the exact solve's after the first step-size
  change; the accepted counts and the step sequence are what they pin, and
  [dofs.md](dofs.md) replays a multi-variable model with the solve made exact.
