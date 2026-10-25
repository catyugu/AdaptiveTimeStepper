# COMSOL field-norm evidence: does the error test average per field or per dof?

What the local COMSOL 6.2 solver actually does when one model holds two fields of very different
size and very different error-to-weight ratio, recorded by `tools/comsol/GeFieldNormSweep.java` and
`tools/comsol/run_twofield.py`, and stored as raw doubles in `test/reference/comsol_twofield.hpp`.

The two candidates, for M fields with N_j dofs each and per-dof error e_i against weight W_i:

```text
per field:  E^2 = (1/M) * sum_j (1/N_j) * sum_{i in j} (e_i / W_i)^2
flat:       E^2 = (1/N) * sum_i (e_i / W_i)^2,  N = sum_j N_j
```

`WeightedRmsNorm` in `include/atstep/detail/state.hpp` computes the flat one. The same norm is used twice in a
run: the spring step `DT * yp_norm <= 1/2` (`include/atstep/detail/controller.hpp`, `first_step`) and
the acceptance criterion. Sections 2 and 3 measure the two uses separately, so the claim does not
rest on the assumption that they share one function.

- COMSOL version and build, read from `batch.log` of the recorded run: `COMSOL 6.2.0.290 progress`,
  i.e. `COMSOL Multiphysics 6.2 (Build: 290)`. No documentation is installed with this build
  (`doc/pdf` holds only the license component), so the two formulas above are the candidates as
  stated for this experiment; the manual page was not re-read here.
- 11 cases, 682 stored samples. Driver wall time about 3 min including the Java compilation; the
  eleven solves report `Total time` of 0-2 s each, so the whole evidence set is cheap to regenerate:
  `python tools/comsol/run_twofield.py --directory .cache/comsol-twofield`, or rewrite the header
  from an existing run with `--header-only`. Raw logs and exports stay in
  `.cache/comsol-twofield/twofield/`, next to the driver's own `batch.log`.
- Every case exports `time,a,b,spread` from `getPVals()/getU()/getUDot()`, its own solver log, and
  one `manifest.csv` row. `spread` is field B's max minus min over its dofs, i.e. the record that the
  field really is spatially uniform; it is checked when the header is written and not stored, because
  a uniform field is one number, so `Sample` in the header holds `time,a,b` only. The driver rejects a
  run whose exported values are not the analytic exponentials, whose field is not spatially uniform,
  whose element count, dof count, variable scale or log completeness is not the one the case asked
  for, or whose control pair does not reproduce.
- Controls: `two_r100_repeat` reproduces `two_r100` byte for byte (all 57 samples, every double),
  `two_equal100` reproduces `single_a` to 2.1e-16 relative, and `single_b100`/`single_b25` are
  byte-identical.

## The model as built

One component, one 1D geometry, two decoupled physics interfaces, both integrated from 1:

```text
A: Global Equations "ge"/"ge1",   uA_t + rateA*uA = 0,  uA(0) = 1
   -> exactly 1 dof.
B: Coefficient Form PDE "c" on the default unit interval, linear elements,
   da = 1, a = rateB, c = 0, f = 0, u(0) = 1, zero-flux boundaries
   -> elements + 1 dofs: 26 for 25 elements, 101 for 100 elements.
```

Field B is `M*u' + rateB*M*u = 0` with the same mass matrix on both sides and no spatial operator,
so every dof of B obeys `u' = -rateB*u` exactly: the spatial discretization is not an approximation
for this state, and the exported `spread` (max minus min over the field's dofs) is 0 for 100
elements and 4.4e-16 for 25. That is what makes B a field whose N_j dofs carry one identical
error-to-weight ratio, and it is why a change of N_j is a change of nothing else.

The mesh is the default `Interval` (`[0,1]`; the 1D geometry feature exposes no `min`/`max`
properties in 6.2 - `properties()` is empty and `set("min", ...)` fails with `Unknown property`), so
the element count is set with `mesh1/size` `custom = on` and `hmax = 1/elements`, and then verified:
`getNumElem()` must equal the case's count, and the log's dof count must equal `elements + 1` plus
one for `ge`. `elements = 25` meshed 25 elements (26 field dofs), `elements = 100` meshed 100 (101
field dofs).

Dof counts as printed by the solver, `Number of degrees of freedom solved for: ...`:

| Case | `field_a` | `field_b` | elements | field dofs | dofs in the log | samples | accepted | rejected | failed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `single_a` | global | none | 100 | 0 | `1.` | 61 | 60 | 0 | 0 |
| `single_b25` | none | pde | 25 | 26 | `26 (plus 2 internal DOFs)` | 37 | 36 | 0 | 0 |
| `single_b100` | none | pde | 100 | 101 | `101 (plus 2 internal DOFs)` | 37 | 36 | 0 | 0 |
| `two_equal100` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 61 | 60 | 0 | 0 |
| `two_r25` | global | pde | 25 | 26 | `27 (plus 2 internal DOFs)` | 57 | 56 | 0 | 0 |
| `two_r100` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 57 | 56 | 0 | 0 |
| `two_r100_repeat` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 57 | 56 | 0 | 0 |
| `two_slow25` | global | pde | 25 | 26 | `27 (plus 2 internal DOFs)` | 52 | 51 | 0 | 0 |
| `two_slow100` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 52 | 51 | 0 | 0 |
| `two_r100_order5` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 29 | 28 | 1 | 0 |
| `two_r100_rtol1e3` | global | pde | 100 | 101 | `102 (plus 2 internal DOFs)` | 182 | 181 | 0 | 0 |

The two internal dofs are the reaction-force slots at the ends of the interval; the solver reports
them apart from the solved dofs and they carry `getUDot == 0` in every sample, which is how the
exporter finds the field's dofs.

Every case runs on `[0,10]` with `tstepsbdf = free`, `tout = tsteps`, `tstepsstore = 1`,
`keeplog = on`, `minorder = 1` and a full solver log. `rateA = 1` in every two-field case;
`rateB = 0.2236` (`two_r*`), `0.1` (`two_slow*`) or `1` (`two_equal100`). `maxorder = 1` except
`two_r100_order5` (5), because the ratio between the two fields' per-dof errors is
`(rateB/rateA)^(p+1)` - `0.05` at order 1 but `0.0025` at order 5 - so pinning the order to 1 makes
that ratio one number for every step instead of a function of the order ramp.

## 0. The setting names, read back from the solver

No name below was guessed: each was found by reading `SolverFeature.properties()`, then set and read
back, and the driver rejects a case whose readback differs from what it asked for
(`max_order_readback`, `rtol_readback`, `tout_readback`, `elements_readback`, `element_order`,
`field_scale`, `global_scale`).

| Property | Meaning | Value used |
| --- | --- | --- |
| `sol1/t1` `tout` | "Times to store" | `tsteps` |
| `sol1/t1` `tstepsstore` | "Store every Nth step" | `1` |
| `sol1/t1` `keeplog` | keep the solver log | `on` |
| `sol1/t1` `maxorder` | maximum BDF order | `1`, and `5` in `two_r100_order5` |
| `sol1/t1` `rtol` | relative tolerance | `0.01`, `1e-3` in `two_r100_rtol1e3` |
| `sol1/t1` `tstepsbdf` | "Steps taken by solver" | `free` |
| `sol1/t1/aDef` `logsampling` | log-row sampling interval | `0`, print every accepted step |
| `sol1/v1/comp1_u` `scalemethod`, `scaleval` | scale of field B | `manual`, `1` |
| `sol1/v1/comp1_ODE1` `scalemethod`, `scaleval` | scale of field A | `manual`, `1` |
| `mesh1/size` `custom`, `hmax` | element size | `on`, `1/elements` |
| `c/cfeq1` `da`, `a`, `c`, `f` | the PDE coefficients | `1`, `rateB`, `0`, `0` |
| `c/init1` `u` | initial value | `1` |
| `c` `prop("ShapeProperty")` `order` | element order | `1` (linear) |
| `ge/ge1` `equation`, `name`, `initialValueU` | the global equation | `uAt+rateA*uA`, `uA`, `1` |

The scale override is not cosmetic. Left to `auto`, the solver picks a scale by itself - it read
back `13` for 26 field dofs and `25` for 101 - and the scale enters the absolute-tolerance term, so
two meshes of the same case would no longer be weighed alike. With `scalemethod = manual` and
`scaleval = 1` the log prints `1` for both variables, which the driver checks in the
`Scales for dependent variables:` block of every case:

```text
Dependent variable u (comp1.u): 1
Global Equations 1 (comp1.ODE1): 1
```

Also found while building this: the Coefficient Form PDE's feature tags are `cfeq1` (the equation),
`zflx1` (`ZeroFluxBoundary`), `init1`, `dcont1`; `t1` accepts `estrat` with only `include` and
`exclude`; and `getU` returns, per sample, the `ge` dof at index 0 when that interface exists, then
the field's dofs interleaved with the two reaction-force slots - which is what `GeFieldNormSweep.java`
walks and what the driver's analytic check re-verifies on every sample.

## 1. The log has no error-estimate column

The per-step error estimate is not in the log. A time-step log row is, exactly:

```text
Step        Time    Stepsize      Res  Jac  Sol Order Tfail NLfail   LinErr   LinRes
   0           0           - out    2    3    2                  0  8.5e-18  4.9e-18
   1   0.0076107   0.0075907 out    3    4    3     1     0      0  1.5e-17  8.8e-18
   2    0.022792    0.015181 out    4    5    4     1     0      0  1.3e-17  8.2e-18
   -          10           - out
```

`LinErr` and `LinRes` are the linear-solver error and residual, not the local error estimate: they
sit at 1e-17 to 1e-16 for these runs (both fields are linear and the Newton solve is exact to
roundoff) while the steps they accompany are the very steps the norm accepted. `Tfail` and `NLfail`
are cumulative counters, `out` marks a stored state, and the `- 10 - out` row is the interpolation
to the end time, not a step. So the norm has to be measured indirectly, which sections 2 and 3 do.

`sol1/t1` (131 properties) and `sol1/t1/aDef` (35 properties) were both enumerated looking for a
switch that prints the error estimate or the norm; there is none, the closest being `convinfo`, which
prints Newton convergence information. `estrat` is not it either: it takes only `include` and
`exclude`.

## 2. E01, first discriminator: the first step is an exact readout of the norm

The first accepted step of a run is `min(0.1% of the span, 0.5 / yp_norm)`, where `yp_norm` is the
weighted norm of `u'` at `t = 0`. Both fields start at `u = 1` with the same `rtol` and the same
manual scale, so both carry the same weight `W`, and their per-dof ratios at `t = 0` are exactly

```text
u'_A = -rateA = -1        u'_B = -rateB        rho = |u'_B| / |u'_A| = rateB / rateA
```

so the two candidate criteria predict different norms for the same state, hence different first
steps. Taking `single_a`'s first step as the calibration removes the weight from the prediction: it
becomes a pure ratio of norms.

```text
per field:  0.5 / (W * sqrt((1 + rho^2)/2))
flat:       0.5 / (W * sqrt((1 + N_B*rho^2)/(1 + N_B)))
```

Measured, with the first step recovered from the export as `t1 - 2e-5` (the two backward-Euler
consistency steps advance the clock by `2 * 1e-3 * 0.1% of the span`, see `stepping-sweep.md`):

| Case | field dofs | first step `h1` | `h1/h1(single_a)` | per-field prediction | flat prediction |
| --- | --- | --- | --- | --- | --- |
| `single_a` | 0 | 0.00550001000005 | 1 | 1 | 1 |
| `two_equal100` | 101 | 0.00550001000005 | 1.000000000 | 1 | 1 |
| `two_r100` | 101 | 0.00759074543243 | 1.380133024 | 1.380133117 | 4.106137 |
| `two_r25` | 26 | 0.00759074543243 | 1.380133024 | 1.380133117 | 3.426300 |
| `two_slow100` | 101 | 0.00773958693866 | 1.407195067 | 1.407195089 | 7.123642 |
| `two_slow25` | 26 | 0.00773958693866 | 1.407195067 | 1.407195089 | 4.629100 |
| `two_r100_rtol1e3` | 101 | 0.00075907454324 | 0.138013302 | 0.138013312 | 0.410614 |

For `rho = 0.2236` the measured 1.380133024 agrees with the per-field prediction 1.380133117 to
6.7e-8 relative, the roundoff of the difference that recovers `h1`, and it is 2.98x away from the
flat prediction. For `rho = 0.1` the measured 1.407195067 agrees with 1.407195089 to 1.6e-8. The
same factor is unchanged at `rtol = 1e-3` (0.138013302 measured against 0.138013312 predicted,
exactly one tenth of the `rtol = 0.01` number), so it is a property of the norm and not of the
tolerance, and it is identical at 26 and 101 field dofs.

The flat reading cannot be rescued by the initial-step rule either. Flat would allow a first step of
`0.01884` (26 dofs) or `0.02258` (101 dofs) for `rho = 0.2236`, and `0.02546` / `0.03918` for
`rho = 0.1` - all above the `0.1% of span` cap of `0.00999998`, so under a flat norm the clamp
would stop binding and the first step would be the `0.00999998` that the `single_b*` cases do show.
It is `0.00759074543` and `0.00773958694` instead. Even reading flat in its most favourable light -
the cap binding, i.e. a ratio of `1.8182` to `single_a` - it is 29% away from the measurement.

The solver log agrees at its own print precision: the `Stepsize` column reads `0.0055`, `0.0075907`
and `0.0077396`, ratios of 1.38013 and 1.40720.

## 3. E01, second discriminator: the same case at 26 and 101 field dofs

The members of each pair are the same model, the same two rates, the same tolerance, order, solver
settings and export schedule; only `hmax` differs. Everything about field B except its size is
bitwise identical, so a per-field average has the same number of terms per field and cannot move,
while a flat average weights field B's 101 (or 26) identical ratios by `N_B`.

There are two pairs, at two error ratios (`rateB/rateA = 0.2236` and `0.1`). They are not two
discriminators but the same one twice: the slow pair gives section 2's first-step measurement a
second, independent instance at `rho = 0.1` (where the flat clamp value is 3.3x to 5.1x away instead
of 3x), and it runs the mesh comparison in the regime where field A's single dof dominates the flat
sum and the dof count dilutes it - the regime the `sqrt(N)` heuristic describes - instead of in the
regime where field B's own ratio dominates.

| Pair | samples | worst relative difference over all times and both fields |
| --- | --- | --- |
| `two_slow100` vs `two_slow25` | 52 vs 52 | 1.5e-15 |
| `two_r100` vs `two_r25` | 57 vs 57 | 4.6e-16 |
| `two_r100` vs `two_r100_repeat` | 57 vs 57 | 0 (byte-identical) |
| `single_b100` vs `single_b25` | 37 vs 37 | 0 (byte-identical) |
| `two_equal100` vs `single_a` | 61 vs 61 | 2.1e-16 |

The step sequence, the order ramp, the accepted count (51 and 56) and the rejected count (0) are the
same for both meshes. That equality is what rules out a flat norm, and the argument is exact rather
than a magnitude comparison: with `alpha` field A's per-dof ratio and `beta` field B's at the same
state, the flat norm of the two meshes is equal only if

```text
(alpha^2 + 26*beta^2)/27 = (alpha^2 + 101*beta^2)/102   <=>   alpha = beta
```

and `alpha = beta` cannot hold here at any step. Both fields are `u' = -rate*u` with the same weight
structure, so at the order-1 error test their per-dof ratios are in the ratio
`(rateA/rateB)^2` - 20 for the `two_r*` pair and 100 for the `two_slow*` pair - and to make them
equal the rates would have to be equal. What a flat norm would therefore predict for the two
meshes, step for step, with `E` proportional to `h^2` at order 1:

| Pair | flat norm ratio, 26 vs 101 dofs | flat prediction for every step size | measured |
| --- | --- | --- | --- |
| `two_r*` | 1.792 | 1.34x smaller at 26 dofs | equal to 4.6e-16 |
| `two_slow*` | 1.539 | 1.24x smaller at 26 dofs | equal to 1.5e-15 |

The control in section 4 closes the loop: when the two rates *are* equal, both formulas coincide,
and the combined run does reproduce `single_a` to 2.1e-16 - so the insensitivity to `N_B` measured
here is a property of the norm, not of the model.

## 4. Controls: where the two formulas must coincide, they do

- `single_a` (1 dof) and `single_b100` (101 dofs) are single-field models, so M = 1 and the two
  formulas are the same formula. `single_b100` equals `single_b25` byte for byte: changing the
  number of dofs inside one field changes nothing, exactly as both formulas say.
- `two_equal100` (`rateA = rateB = 1`, 102 dofs) must reproduce `single_a` under either formula,
  because `rho = 1` and the two averages coincide whenever the fields are identical. It does, to
  2.1e-16 over all 61 samples: the extra 101 dofs neither add error nor dilute it.
- `two_r100_repeat` reproduces `two_r100` byte for byte, so the pairs of section 3 are not two lucky
  draws.

## 5. The step sequences, and why the plateau is not a counter-argument

Measured first six step sizes, from the log's `Stepsize` column:

| Case | first six step sizes | plateau | accepted |
| --- | --- | --- | --- |
| `single_a` | 0.0055, 0.011, 0.022, 0.044, 0.088, 0.088 | 1.0 (the 0.1*span cap) | 60 |
| `single_b100` | 0.01, 0.02, 0.04, 0.08, 0.16, 0.32 | 0.32 | 36 |
| `two_r100` | 0.0075907, 0.015181, 0.030363, 0.060726, 0.060726, 0.12145 | 0.3935 | 56 |
| `two_r25` | 0.0075907, 0.015181, 0.030363, 0.060726, 0.060726, 0.12145 | 0.3935 | 56 |
| `two_slow100` | 0.0077396, 0.015479, 0.030958, 0.061917, 0.061917, 0.12383 | 0.80244 | 51 |
| `two_slow25` | 0.0077396, 0.015479, 0.030958, 0.061917, 0.061917, 0.12383 | 0.80244 | 51 |
| `two_r100_order5` | 0.0075907, 0.0075907, 0.015181, 0.030363, 0.060726, 0.12145 | 1.0 (cap) | 28 |
| `two_r100_rtol1e3` | 0.00075907, 0.0015181, 0.0030363, 0.0060726, 0.012145, 0.02429 | 0.13859 | 181 |

The first four steps of every run are the initial step doubled (the growth limit), so they carry the
ratio of section 2 unchanged; from the fifth step the error test takes over and the doubling pauses.
The measured early ratio between `two_r100` and `single_a` is 1.380127, 1.380091, 1.380136,
1.380136 - the clamp's number again.

Late in the run the sequences stop being a pure scaling, and this is the one place where a naive
reading of the flat idea could look attractive: the fast field A decays to 4.5e-5 while the slow
field B stays near 0.1-0.4, so the absolute-tolerance floor in the weight is no longer proportional
to each field's own magnitude and the fields' relative weights drift. For the `two_r*` family the
three plateaus are 0.32 s (B alone), 0.3935 s (both) and 0.6336 s (A alone, before the cap): the
combined run sits between the two single-field runs, which is what a per-field average of two
comparable norms gives. With `E` proportional to `h^2` those two single-field plateaus predict
0.3746 s for the per-field average and 0.3208 s for a flat average dominated by field B's 101
identical dofs; the measurement, 0.3935 s, is 5% from the first and 23% from the second. Because
that plateau depends on how the weight floor cancels between the fields, this is supporting evidence
only; sections 2 and 3 carry the claim.

## 6. What this settles for the library

**Defended.** In COMSOL 6.2.0.290, for a model with two fields, the weighted norm behind both the
spring step and the acceptance criterion is a **per-field average**: with `M = 2`, field A
contributing 1 dof and field B contributing 26 or 101, the measured norm behaves as
`sqrt((N_A^2 + N_B_field^2)/M)` and not as `sqrt((N_A^2 + N_field*N_B^2)/(1 + N_field))`, to 7
significant digits in the first step of the run and exactly (4.6e-16 and 1.5e-15) under a fourfold
refinement of field B.

`WeightedRmsNorm` in `include/atstep/detail/state.hpp` is the flat average, so it is the wrong rule for a model
with more than one field, and the discrepancy is not a rounding effect: for these models the two
rules differ by a factor of 2.5 to 5 in the first step, by 24-34% in every step of a mesh comparison,
and the accepted step count moves from 60 to 56 instead of to roughly 20. The repair is not a
different exponent or a constant: the norm needs to know the *partition* of the state into fields -
which dofs belong to which dependent variable, and how many fields there are - which it currently
never receives, since it sees one flat `error`/`accepted` pair and a single `ErrorTest`. The
partition has to enter through `ErrorTest`/`ErrorContext` or through a norm constructed with it.

**Not defended.**

- Whether COMSOL counts *physics interfaces* or *dependent variables* as `M`. This model has one
  variable per interface, so the two readings coincide; nothing here separates them.
- The tolerance formula itself. The measured first step implies an effective weight of
  `0.011 = rtol*(1 + atolfactor)` at step 0 for both fields, where the `scaled` branch of
  `WeightedRms` would give `0.012` and the `unscaled` branch `0.011`. That is a 9% difference in an
  absolute step size, and it was not measured independently of the backward-Euler advance used to
  recover `h1`, so no claim is made about `scaled` versus `unscaled`, about `ewtrescale`, or about
  how the weight's absolute term decays with the accepted-step index.
- Anything at order above 1. The quantitative ratio argument is an order-1 argument
  (`beta/alpha = (rateB/rateA)^2`); `two_r100_order5` reproduces the same first step and a
  per-field-looking sequence, but no number here depends on it.
- What a rejected step measures: there is exactly one rejection in the whole sweep
  (`two_r100_order5`, `Tfail = 1`), so no case exercises rejection of a two-field state.

## 7. Not run: E02, the absolute-tolerance axes

E02 was left out. The property names its six cases need are already read back from `sol1/t1`, so the
work is a case table rather than a discovery: `atolglobalmethod` (`scaled`), `atolglobalvaluemethod`
(`factor`), `atolglobal` (`0.001`), `atolglobalfactor` (`0.1`), `atolfactor` (`0.1`), `atolmethod`,
`atolvaluemethod`, `ewtrescale` (`on`), `tderglobalmethod` (`auto`), `tderglobalfactor` (`1`),
`tdermethod`, `tderfactor`. Whether `ewtrescale` off, a manual `atolglobal` and an unscaled
`atolglobalmethod` reproduce the `0.011`-versus-`0.012` question above is exactly what those axes
would settle.

## 8. Limits of this evidence

- One physics pairing (Global Equations + Coefficient Form PDE), one spatial dimension, one span
  (`[0,10]`) and one export schedule. The claim is about the averaging rule, which these cases pin
  down, not about any other setting.
- Field B's dofs are identical by construction, so the experiment measures the averaging over fields
  and the scaling with N_j, and says nothing about how a field whose own dofs disagree is averaged.
  For that, a field with a spatially non-uniform error is needed.
- The log prints `Stepsize` with 5 significant digits; the 7-digit agreement in section 2 comes from
  the exported times, not from the log.
