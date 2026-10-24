# AdaptiveTimeStepper

[![CI](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml/badge.svg)](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Header-only C++20 adaptive BDF integration, variable order 1–5, no external
dependencies. The caller solves the equations; the library owns prediction,
history, the error test, step size, order and output. The first-order kernel is
checked against recorded COMSOL 6.2 scalar runs. A second-order-in-time equation
is written as the first-order pair `(u, v = u')`; the library supplies the pair
and the condensed coefficients, and the caller owns the solve.

## Example

Integrate `y' + y = 0` from `y(0) = 1`:

```cpp
#include <atstep/atstep.hpp>
#include <iostream>
#include <vector>

int main()
{
    using State = std::vector<double>;
    atstep::TimeStepper<State> stepper(atstep::Settings {},
        [](const atstep::StepRequest<State>& request, State& y) {
            const State& known = request.derivative_history();
            y[0] = -known[0] / (request.derivative_weight() + 1.0);
            return atstep::SolveOutcome {true, 1};
        });

    atstep::RunHooks<State> hooks;
    hooks.on_store = [](double time, const State& y, bool) {
        std::cout << time << ' ' << y[0] << '\n';
    };
    const std::vector<double> times {0.0, 0.5, 1.0};
    stepper.start(0.0, State {1.0}, 1.0);
    const auto report = stepper.run(times, hooks);
    return report.steps > 0 ? 0 : 1;
}
```

The solve receives the predictor in `y` and replaces it with the solved state.
For `M y' + K y = f`, assemble `(cj*M + K)*y = f - M*known`; for a nonlinear
residual, `request.derivative(y, dy)` evaluates `y'` in correction form. Return
`{false, iterations}` when the solve fails; the step is retried smaller.

## State and norm

A state is a default-constructible, copyable value with
`void axpy(double a, const State& x)` (`*this += a*x`). `std::vector<double>`
is supported directly. The default error norm is COMSOL's field-and-DOF
weighted RMS:

```text
sqrt( (1/M) * sum_j (1/N_j) * sum_i (E_i / (A_i + R*abs(Y_i)))^2 )
```

`R` is `relative_tolerance`; `field_sizes` gives `N_j` for each of the `M`
fields (empty means one field holding the state). With
`scale_absolute_tolerance`, a field of scale reference `S` has

```text
A_i = R*factor * (beta + norm2(S))
beta = (1 - d)*norm_inf(S) + d   if 0 < norm_inf(S) < 1
beta = 1                         otherwise
d = exp(-j/5),  j = accepted steps
```

and `A_i = R*factor` unscaled. A custom norm receives `(error, reference)` and
may be passed to the three-argument constructor; it can, for example, scale a
velocity field or a variable the documented COMSOL formula does not cover.

## Settings

| Setting | Default | Meaning |
| :--- | :--- | :--- |
| `min_order`, `max_order` | 1, 5 | Adaptive order range; the first attempts use the history they have |
| `relative_tolerance` | 1e-2 | Relative error weight `R` |
| `absolute_tolerance_factor` | 0.1 | Absolute tolerance as a factor of `R` |
| `scale_absolute_tolerance` | true | COMSOL's scaled absolute tolerance |
| `field_sizes` | empty | DOF count of each field; empty is one field |
| `initial_step_fraction` | 1e-3 | Automatic first step, of the final time's magnitude (at least the span) |
| `max_step` | 0 | Absolute maximum step in run time units; 0 is no limit |
| `min_step_fraction` | 1e-10 | Smallest step the controller proposes, as a fraction of the span |
| `step_mode` | `free` | Step placement |
| `store_mode` | `interpolate` | Output selection |
| `interpolate_end_time` | true | Allow a step past the final output time |
| `manual_step`, `manual_order` | 0, 2 | Fixed step and order in manual mode |
| `backward_euler_fraction` | 1e-3 | Artificial initialization steps, as a fraction of the initial step |

Step modes: `free`, `intermediate` (a step in every output subinterval),
`strict` (end every step at an output time), `manual` (the time-step error test
off; a failed algebraic solve still reduces the step).

Store modes: `interpolate`, `closest` (the first internal step reaching each
output time), `steps` (the internal steps; a final step past the end is stored
at the end when `interpolate_end_time` is on).

## Running

`start(t0, initial, end)` performs two artificial backward-Euler solves, the
consistent initialization; the returned state is labelled `t0` and its
derivative is what the first step is read from. `start(t0, initial, derivative,
end)` uses consistent initial values supplied by the caller and takes no
artificial step. The second overload is the way to restart from a solution
without perturbing it.

`run(times, hooks)` takes any strictly increasing set of times that are the
start time or after the initialized time; the initial state is stored when the list
starts at `t0`, and always in `steps` mode. `advance()` takes
one accepted step for a caller-driven loop. `RunHooks::on_step` sees every
accepted step; `on_store` sees every stored value and whether it was
interpolated. `time()` and `state()` are the newest accepted level.
`interpolate(t, out, derivative)` writes the dense value or first derivative
into a buffer of the caller's; it is the BDF polynomial through the last
`order + 1` levels, so times older than the last step are extrapolated.

An integrator is neither copyable nor movable; `start` begins a run again. State
operations must not throw: history is published in place.

## Second-order equations

`M u'' + C u' + K u = f` is integrated as the pair `(u, v = u')`:

```text
u' - v            = 0
M v' + C v + K u  = f
```

`SecondOrderState<State>` is that pair. The solve is the condensed system of
the displacement:

```text
u'  = cj u + known1
u'' = cj^2 u + known2,   known2 = cj known1 + known(v)
(K + cj C + cj^2 M) u = f - M known2 - C known1
```

`second_derivative_history(request, known2)` builds `known2`;
`set_velocity(state, cj, known1)` enforces `v = cj u + known1` after a
successful solve. The default error norm measures the displacement, because the
pair exposes the displacement to the norm; a custom norm that reads `value` and
`derivative` can include the velocity. The COMSOL rule that includes the
auxiliary velocity when consistent initialization is on, and the derivative
tolerance that weights it, are not implemented.

A third or higher derivative in time is the caller's to write as a first-order
system.

## Build and test

```text
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix .cache/install
```

Tests cover the controller, the norm, convergence, the pair, rejection, restart,
output modes, invalid input and the recorded COMSOL trajectories. Examples:
[heat](examples/heat_equation.cpp), [oscillator](examples/second_order.cpp),
[Brusselator](examples/vector_state.cpp), [motor](examples/coupled_motor.cpp).

## COMSOL alignment and limits

The 21 recorded scalar configurations agree on accepted steps, rejections,
stored times, values, derivatives and logged orders to `1e-8`; the thermoelastic
regression agrees on 161 samples. See [COMSOL alignment](doc/comsol-alignment.md)
for the measured settings and the [research note](doc/research/comsol-time-dependent-solver.md)
for the documentation quotes.

The library does not provide algebraic solvers, events, generalized alpha,
Runge–Kutta, COMSOL's algebraic-variable error-test options, per-variable
tolerances, a derivative tolerance, COMSOL's `Safety factor for Backward Euler`
or `Rescale after initialization`, store-every-Nth-step, a nonlinear controller,
or an automatic maximum step of its own. The default maximum step is no limit;
COMSOL's automatic value is not documented. The scaled absolute tolerance is
validated on a one-field, one-dof scalar run with `abs(y) <= 1`; the field-aware
formula is the documented one, but multi-field COMSOL runs were not recorded.
The second-order pair is not validated against a COMSOL second-order model.

## Documentation

- [Design](doc/design.md): ownership, history, error weights and output.
- [COMSOL alignment](doc/comsol-alignment.md): measurements, assumptions, gaps.
- [Coupled multiphysics](doc/multiphysics.md): thermoelastic and motor cases.
- [Research notes](doc/research/comsol-time-dependent-solver.md): source quotes.

Licensed under [MIT](LICENSE).
