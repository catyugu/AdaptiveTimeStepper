# Coupled multiphysics cases

Both cases advance all fields in one state, using one step size, order and
error test. Their solve callbacks assemble the coupled equations at the new
time level. They are lumped models with no spatial mesh discretization.

## COMSOL reference: reciprocal thermoelasticity

[test_comsol_thermoelastic.cpp](../test/reference/test_comsol_thermoelastic.cpp)
compares a thermoelastic oscillator against an independently solved COMSOL
6.2.0.290 model. Mechanical deformation changes temperature, and thermal
expansion feeds back into motion. The nondimensional state is `(x,v,theta)`:

```text
x'     = v
v'     = -4*x - 4.5*v + 0.5*theta
theta' = -0.5*v - 1.5*theta
```

Initial values are `(1.75,-2.75,-0.25)`, integrated over `[0,5]`. The solution
contains three decay modes:

```text
x     = exp(-t) + 0.5*exp(-2*t) + 0.25*exp(-3*t)
v     = -exp(-t) - exp(-2*t) - 0.75*exp(-3*t)
theta = exp(-t) - exp(-2*t) - 0.25*exp(-3*t)
```

Its quadratic energy `E = 2*x^2 + v^2/2 + theta^2/2` satisfies
`E' = -4.5*v^2 - 1.5*theta^2`. The reciprocal coupling terms cancel.

### Reference generation and acceptance

[GeThermoelasticTransient.java](../tools/comsol/GeThermoelasticTransient.java)
implements the three equations in COMSOL Global Equations, with `q0=x`,
`q1=v`, `q2=theta`. It uses COMSOL's BDF solver, not library-produced states.
The numerical settings match the C++ run:

| Setting                                | Value                                                  |
| -------------------------------------- | :----------------------------------------------------- |
| Maximum BDF order                      | 5                                                      |
| Relative tolerance                     | 1e-7                                                   |
| Absolute tolerance                     | unscaled, factor 0.1                                   |
| Backward-Euler initialization fraction | 1e-5                                                   |
| Output                                 | every internal solver step; land at the final time     |
| Algebraic solve                        | Jacobian updated each iteration; tolerance factor 1e-8 |

The smaller initialization fraction limits the initial-state time shift when
checking the analytic solution. It is applied to both solvers.

`getPVals`, `getU` and `getUDot` export full-precision times, states and
derivatives. [run_thermoelastic.py](../tools/comsol/run_thermoelastic.py)
checks finiteness, increasing times, analytic accuracy, energy decay and the
solver log before generating
[comsol_thermoelastic.hpp](../test/reference/comsol_thermoelastic.hpp).
The raw CSV, model and logs stay in the selected run directory.

The recorded run has 160 accepted steps, one error-test rejection, no failed
solves and 161 stored states. The C++ regression compares every internal time,
all three states and all three derivatives, plus those solver counts. The
absolute comparison limit is `1e-7` for each quantity. Separate analytic and
energy checks validate the model and the export's degree-of-freedom ordering.
Reference values are only read after the C++ solve; they never enter its callback.

The local MSVC Release comparison measured maximum differences of `3.79e-9`
in time, `7.13e-10` in state and `1.28e-8` in derivative. COMSOL's own
maximum analytic errors were `5.25e-7` in state and `1.18e-6` in derivative.

This verifies a reduced coupled system implemented with Global Equations.
It does not validate COMSOL's spatial Solid Mechanics/Heat Transfer interfaces
or arbitrary continuum thermoelastic models.

Regenerate with COMSOL installed:

```text
python tools/comsol/run_thermoelastic.py --directory .cache/comsol-thermoelastic
```

To rebuild the header from the same raw evidence without running COMSOL:

```text
python tools/comsol/run_thermoelastic.py --header-only --directory .cache/comsol-thermoelastic
```

The scalar decay reference has its own entrypoint, `tools/comsol/run_decay.py`.
Both scripts use `tools/comsol/common.py` only for tool discovery, compilation
and batch execution; their model configuration and export checks are independent.
Normal CTest runs use the committed reference and require no COMSOL installation.

## Example: electrical, thermal and mechanical motor dynamics

[coupled_motor.cpp](../examples/coupled_motor.cpp) models a DC motor with
temperature-dependent winding resistance. The state is current `i` in amperes,
rotor speed `omega` in radians per second and temperature rise `theta` in kelvin:

```text
L*i'     + R(theta)*i + k*omega = V
J*omega' + b*omega             = k*i
C*theta' + G*theta             = R(theta)*i^2
R(theta) = R0*(1 + alpha*theta)
```

Current generates torque and heat. Rotation creates back EMF, and temperature
changes electrical resistance. Every Newton update solves the full coupled
Jacobian by block elimination. A failed nonlinear solve reports nonconvergence
to the stepper for retry.

The example uses `V=2 V`, `L=0.5 H`, `J=0.2 kg m^2`, `C=2 J/K`, `R0=1 ohm`,
`alpha=0.1/K`, `k=0.5` in consistent SI torque/back-EMF units,
`b=0.1 N m s` and `G=0.5 W/K`. It starts from rest at ambient temperature
and integrates to 60 seconds, using maximum order 2 and relative tolerance 1e-6.

It prints the trajectory and checks finite, nonnegative states, an independent
stationary solution and the power balance

```text
d/dt (L*i^2/2 + J*omega^2/2 + C*theta)
    = V*i - b*omega^2 - G*theta.
```

Joule heating transfers energy into the winding; mechanical friction is lost
to the surroundings. At equilibrium, the reference solves the scalar thermal
balance by bisection, then derives current and speed. The example requires
maximum final-state error below `2e-5` and power defect below `1e-7 W`.
This example has no COMSOL comparison.

The local Release run took 880 steps without rejection or failed solves.
At 60 seconds, `(i,omega,theta)` was approximately
`(0.56068099,2.80340494,0.67090812)`. Its maximum stationary-reference error
was `1.79e-8`; the worst power defect was `5.23e-12 W`.

Build and run both cases through CTest:

```text
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```
