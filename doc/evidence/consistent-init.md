# Consistent initialization off

COMSOL's Time-Dependent Solver can skip its consistent initialization, and the
library's `Settings::initialization` is that switch. What "off" leaves behind was
measured on two models, each one Global Equation, both integrated from 0 to 10
with `maxorder = 5`, `rtol = 0.01`, the scaled tolerance at its defaults, every
accepted step in the log, and the Fully Coupled node's Jacobian rebuilt at every
iteration, so each record is the time integration's alone:

```text
u' + u = 0,                 u(0) = 1
qtt + 0.1 qt + q = 0,       q(0) = 1,   q'(0) = 100
```

Recorded by `tools/comsol/GeConsistentInit.java` into
`test/reference/comsol_consistent_init.hpp`, six runs on the two models and the
switch on and off.

## The setting

The solver property is `consistent` on `sol1/t1`. Read back from COMSOL 6.2.0.290
and run under the sweep:

| Value     | Meaning                                                    |
| --------- | ---------------------------------------------------------- |
| `bweuler` | the default: the two artificial backward-Euler steps       |
| `on`      | the same run, bit for bit: 45 samples, the same first step |
| `off`     | no initialization                                          |
| `none`    | refused: `Invalid_property_value`                          |

`on` is a synonym the solver echoes back rather than normalizes, so the sweep sets
only `bweuler` and `off` and requires the property to read them back.

## What the record starts from

| Case               | `consistent` | first sample: time, state, derivative        | first step |
| ------------------ | ------------ | -------------------------------------------- | ---------- |
| `on`               | `bweuler`    | `0, 0.9999800002999961, -0.9999800002999961` | 0.00602    |
| `off`              | `off`        | `0, 1, 0`                                    | 0.01       |
| `off_initial_step` | `off`        | `0, 1, 0`                                    | 0.02       |
| `off_manual`       | `off`        | `0, 1, 0`                                    | 0.1        |
| `second_on`        | `bweuler`    | `0, 1.0019999966996045, 99.99977997033008`   | 1.05e-4    |
| `second_off`       | `off`        | `0, 1, 100`                                  | 8.485e-5   |

Three things follow, and the library now does all three:

- The level is the initial state itself, labelled with the start time. Nothing
  reconciles it with the equations: `u(0)` stays 1, where the two probes leave
  `0.99998`.
- The derivative the run starts from is zero for a first-order equation. A
  second-order form keeps the second initial value it was given: `the state's u'`
  is 100, the acceleration is 0, and the pair therefore starts at `(1, 100)` with
  the derivative `(100, 0)`. The first step is the one that derivative justifies,
  `0.5 / E(y')` with the pair's `sqrt(1/2)`, which is 8.485e-5 here.
- The first step is the rule's own, uncapped: 0.1% of the span (0.01), the step a
  case names (0.02), or the manual step (0.1). With the initialization on, the
  measured derivative caps the automatic step at 0.00602 instead.

The `off_manual` step is 0.1 s, whose clock accumulates to `9.99999999999998`: the
record stops there after 100 steps, and the library's step that would stop within
the clock's own roundings lands on the end instead, so its last stored time is
`10`, `1.95e-14` away.

## The derivative a stored sample carries

COMSOL reports the derivative of the polynomial it stores, not the one the solve
left. At a level the two differ in the last digits the step gives them: over
`off`'s 35 samples the dense derivative is within `6.2e-15` of the record and the
solve's own within `1.7e-5`, and `on`'s worst solve-derivative gap is `4.0e-6` at
its last sample, where the derivative is `-4.5e-5`. The replay in
`test/reference/test_comsol_consistent_init.cpp` therefore compares the record with
`interpolate(time, out, 1)`, which is also what a second-order run reports at an
interpolated output time ([second-order.md](second-order.md)).

## The replay

Every case is replayed from the model's own initial values and the settings the
record carries:

| Case               | accepted  | rejected | samples   | worst `dt` | worst `du` |
| ------------------ | --------- | -------- | --------- | ---------- | ---------- |
| `on`               | 44 / 44   | 5 / 5    | 45 / 45   | 1.0e-11    | 4.7e-12    |
| `off`              | 34 / 34   | 5 / 5    | 35 / 35   | 1.1e-13    | 6.2e-15    |
| `off_initial_step` | 34 / 34   | 4 / 4    | 35 / 35   | 2.2e-13    | 2.4e-15    |
| `off_manual`       | 100 / 100 | 0 / 0    | 101 / 101 | 2.0e-14    | 7.8e-16    |
| `second_on`        | 49 / 49   | 5 / 5    | 50 / 50   | 1.5e-10    | 7.1e-09    |
| `second_off`       | 48 / 48   | 3 / 3    | 49 / 49   | 7.1e-15    | 4.0e-13    |

Every stored time, state and derivative of the six runs is compared at `1e-8`, so
`second_off` is the switch's own arithmetic to round-off, and `second_on` keeps
the residual the second-order record already has.

## Not settled

A first step the tolerance refuses many times over was measured and is **not**
reproduced. `off` at `rtol = 1e-5` proposes the automatic 0.01 s step, refuses it,
and walks down: the record takes 5 rejections to 9.765625e-06 s (0.01 times
`0.25^5`) before accepting, while the library accepts after 3, at 1.5625e-4 s.
The difference is the first step's error estimate when the step is far below the
distance back to the level it starts from — the library's `ck` shrinks with that
ratio and the record's does not — and it is not the initialization: the same
measurement with the initialization on reproduces (`rtol_1e5`, 98 accepted steps,
3 rejections), its first step being the small one the derivative justifies. That
case is not in the record; a rule for it needs its own sweep.

## Reproduce

```text
python tools/comsol/run_consistent_init.py --directory .cache/comsol-consistent
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

COMSOL is required only to regenerate the reference; the test replays the
recorded data and has no COMSOL dependency.
