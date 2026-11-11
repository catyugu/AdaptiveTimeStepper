# Settings

Construct `atstep::Settings`, adjust the relevant fields, then pass it to the
stepper. The stepper owns a copy; later changes to the original do not affect it.
All time values use the same unit as the model's time coordinate.

```cpp
atstep::Settings settings;
auto& adaptive = std::get<atstep::AdaptiveSteps>(settings.steps);
adaptive.max_order = 2;
adaptive.max_step = atstep::ConstantMaxStep {0.1};
settings.error.relative = 1e-5;
settings.storage = atstep::Storage::steps;
```

## Adaptive steps

`settings.steps` is a variant of `AdaptiveSteps` (the default) and `ManualSteps`.
Adaptive settings apply only while it holds `AdaptiveSteps`.

| `AdaptiveSteps` field | Default | Meaning |
| --- | --- | --- |
| `placement` | `Placement::free` | Step placement, described below |
| `min_order` | 1 | Minimum selected order: 1 or 2; startup uses available history |
| `max_order` | 5 | Maximum order, from `min_order` through 5 |
| `initial_step` | unset | Automatic first step; set a positive value to override |
| `max_step` | `AutomaticMaxStep {}` | Maximum-step alternative, described below |
| `min_step_fraction` | 1e-10 | Positive minimum-step fraction of the run span |
| `interpolate_end_time` | true | Allow a step past the end and interpolate its value; ignored for strict placement |

Free placement lets error control choose steps. `Placement::intermediate` adds
at least one step per output interval; `Placement::strict` lands on each requested
time, using smaller internal steps as needed.

| `max_step` alternative | Meaning |
| --- | --- |
| `AutomaticMaxStep {}` | 0.1 of the run span |
| `ConstantMaxStep {h}` | Positive absolute limit `h` |
| `ExpressionMaxStep {function}` | Time-dependent limit; the function's absolute value is used |

The automatic first step is bounded by 0.1% of the interval remaining after
initialization and `0.5 / norm(y')`. An explicit `initial_step` replaces that
rule, while maximum-step and placement constraints still apply.

Maximum order 5 follows COMSOL's generic Time solver. Heat Transfer's automatic
configuration uses a maximum of 2. This permits adaptive orders 1–2; it does
not fix every step at order 2. See the [default-order evidence](comsol-alignment.md#default-order).

## Manual steps

Switch the variant to select manual stepping:

```cpp
atstep::ManualSteps manual;
manual.step = 0.01;
manual.order = 2;
settings.steps = manual;
```

| `ManualSteps` field | Default | Meaning |
| --- | --- | --- |
| `step` | 0.0 | Required positive step in the model's time unit |
| `order` | 2 | Target BDF order, 1–5 |
| `startup_fraction` | unset | Initial fraction of `step`; unset uses the per-order default |
| `startup_growth` | unset | Startup growth multiplier; unset uses the per-order default |
| `interpolate_end_time` | true | Allow an internal step past the end |

Manual stepping disables the local error acceptance test. It ramps toward the
configured size and order using the [startup curve](comsol-alignment.md#stepping-and-output).
Failed algebraic solves still reduce the step; after such a reduction the run
can remain below the configured size. Tolerances also remain relevant to the
error norm used during initialization.

## Error tolerances and variables

| `settings.error` field | Default | Meaning |
| --- | --- | --- |
| `relative` | 1e-2 | Positive relative tolerance `R` |
| `absolute.kind` | `AbsoluteTolerance::Kind::factor` | Factor of `R`, or `manual` absolute value |
| `absolute.value` | 0.1 | Nonnegative factor, or positive manual absolute tolerance |
| `scaled` | true | Apply the criterion to variables divided by their scales |
| `variables` | empty | Contiguous variable blocks in state order; empty means one block at scale one |

For an absolute tolerance independent of `relative`, set both its kind and value:

```cpp
settings.error.absolute.kind = atstep::AbsoluteTolerance::Kind::manual;
settings.error.absolute.value = 1e-8;
```

Each `Variable` has `dofs` (default 0), `scale` (default 1.0) and
`scaled_tolerance` (default true). Block sizes must cover the state exactly;
scales must be finite and positive. For example, a six-dof temperature field
followed by a scalar state can use:

```cpp
settings.error.variables = {
    atstep::Variable {6, 300.0, false},
    atstep::Variable {1, 1.0, true}
};
```

`scale` divides the solution in the scaled criterion. `scaled_tolerance` controls
whether the scaled absolute weight follows the solution; disable it to match
manual scaling or disabled tolerance updates in COMSOL. The library does not
derive scales from the model. Unscaled error testing ignores these scales.

The default norm averages weighted errors per variable, then over its dofs.
Relative weights follow the accepted internal state; absolute weights follow
the most recently stored state, including interpolated output. The complete
[error-norm formula](comsol-alignment.md#error-norm) is maintained with its evidence.
A custom norm supplied to the stepper replaces this default criterion.

## Initialization and storage

| `Settings` field | Default | Meaning |
| --- | --- | --- |
| `initialization` | `Initialization::backward_euler` | Two initialization probes; `Initialization::off` skips them |
| `initialization_fraction` | 1e-3 | Positive probe-size multiplier; read only when probes run |
| `storage` | `Storage::interpolate` | Output selection, described below |
| `store_every` | 1 | Positive N for internal-step storage |

The probe basis is an explicit initial step, 0.1% of the span for automatic
stepping, or 0.1% of the manual step. Each of the two probes uses that basis
multiplied by `initialization_fraction`. Initial reported values can therefore
differ slightly from the supplied state. Passing an explicit initial derivative
to `start` skips probes regardless of the initialization setting. See
[initialization evidence and its known limitation](comsol-alignment.md#initialization).

| Storage | Reported values |
| --- | --- |
| `Storage::interpolate` | Values at requested output times |
| `Storage::closest` | First internal step reaching each request; requests reached by one step merge |
| `Storage::steps` | Accepted internal steps, filtered by `store_every` |

Internal-step storage counts the initial state as step 0 and always stores the
endpoint. Step placement and storage are independent: strict placement controls
where integration lands, while storage controls which values callbacks receive.

Settings are validated by the stepper. `atstep::validate(settings)` also checks
configuration before constructing a run; state partitions and time-dependent
limits receive further checks when used.
