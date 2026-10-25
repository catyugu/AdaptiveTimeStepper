# The stepping and output settings

What COMSOL 6.2 does with the settings of its Time-Dependent Solver, recorded by
`tools/comsol/GeSteppingSweep.java` into `test/reference/comsol_stepping.hpp`
(19 cases, 1231 samples). The model is the Global Equations interface
`u' + u = 0`, `u(0) = 1`, over `[0,10]` with `maxorder = 5`, `rtol = 0.01`,
`tout = tsteps` and every accepted step in the log; nothing else is named per
claim.

## The automatic maximum step is 0.1 of the span

Exactly 1.0 s at a span of 10 and 10.0 s at a span of 100, where six consecutive
steps hold it, and 1.0 s for `start 2, end 12` — so it is the span and not the
end time. A constant maximum step is an absolute time (0.05 s is held), and an
expression one is evaluated while solving (`0.05 (1 + 0.5 sin t)` reaches
0.075 s at `t = pi/2`).

## The automatic initial step is 0.1% of what remains

`h1 = min(0.1% of the interval left after the initialization, 0.5 / norm(y'))`.
The first stored interval is not that step: two backward-Euler probes of
`1e-3 × "Fraction of initial step" (1e-3)` of the run's first step precede it, so
`h1 = h_first + 2e-5` at a span of 10, and that offset is what makes the rule
readable — a named initial step of 0.02 s shows `h1 = 0.02004` s.

For `start 2, end 12` the probe reads 0.01 s, i.e. 0.1% of 10 and not of 12; and
the fraction is of the interval that remains: at `rtol = 1` the derivative
condition stops binding and the first step is `0.00999998 = 1e-3 (10 - 2e-5)`,
not `0.01`.

## A manual step follows a per-order startup curve

The order ramps `1, 1, 2, 3, …` to `bdforder` and then holds, one step per
accepted step. The size is `h_k = min(fraction × step × growth^(k-1), step)`,
with `(fraction, growth)` defaults `(1, 1.5)`, `(1, 1.5)`, `(0.1, 1.2)`,
`(0.01, 1.15)`, `(0.001, 1.1)` for orders 1 to 5, and the caller's own values
override them. At orders 1 and 2 the first step already is the manual step, so
the growth rate never acts.

## "Store every Nth step" anchors on the initial state

For N = 3 the stored steps are `t0` and then every third accepted step, plus the
end time, which is stored whatever the modulus is. An anchor on step 1 or 2 would
store a different set and neither appears. A rejected attempt does not move the
phase.

## The log is sampled unless told otherwise

`sol1/t1/aDef` `logsampling` defaults to 0.005 and prints a subset of the accepted
steps that varies between runs of the same input; `logsampling = 0` prints all of
them, and only the log's timestamps then differ between runs. The recorded header
therefore stores the samples and the last row's cumulative counters, and no
per-step log list.

## What this does not settle

- A rejection inside a store-N phase: every rejection here falls in the last step,
  after all stored samples.
- The anchor for output-time store modes (`tlist`, `tstepsclosest`), which no case
  here exercises.
- The order 4 startup curve, which is only read back (0.01 and 1.15), and a manual
  step other than 0.1 s.
- How a physics interface whose own default differs from the Global Equations one
  behaves: nothing here covers that.
