# Design

The library is header-only C++20. `settings.hpp` holds configuration,
`formula.hpp` classical polynomial stencils, `controller.hpp` the scalar
step/order decisions, `state.hpp` state arithmetic and error norms, and
`stepper.hpp` the run, BDF history and output handling. The controller does
not depend on the state, equation or norm.

## The solve boundary

The caller supplies the state type, algebraic solve and optionally the norm.
The state needs copying and `axpy(y,a,x)`, meaning `y += a*x`; the library
does not inspect its matrix, operator or storage layout.

For a first-order residual `M y' + K y = f`, the request supplies
`y' = cj*y + known`. The caller assembles `(cj*M + K)*y = f - M*known` using
`derivative_weight()` and `derivative_history()`. `derivative(y)` evaluates
the equivalent correction form, reducing cancellation in a nonlinear loop.
The solve receives the predictor as its initial guess and returns convergence
and iteration information.

The standalone `StepFormula` still computes classical level-based stencils
for directly discretized higher derivatives. It is not the first-derivative
formula of the fixed-leading BDF method after a step change. Inertia is best
written as a first-order system `(u,v=u')`, as in `examples/second_order.cpp`.

## Accepted and trial history

First-order integration retains scaled divided-difference columns and their
step distances. A trial copies and rescales this history for the proposed
step, forms the predictor and derivative predictor, and solves. Its correction
is the difference between the solved state and predictor. The error test and
order comparisons use this correction and the trial difference columns.

Acceptance updates and commits the trial columns. Rejection leaves the
accepted columns and time levels intact. This also makes a failed first trial
safe: its retry rescales the initialization derivative for the smaller step.
`start` resets both the numerical history and accepted-step tolerance index.

`StepController` consumes scalar error estimates and returns the next step and
order. Startup, settled-step order selection, rejection and convergence failure
are independently tested. Maximum-step arithmetic preserves the measured
transition onto a constant step; an overshoot of one floating-point bit can
otherwise keep resetting the settled-step counter.

## Output is part of the measured contract

The run places steps according to free, intermediate, strict or manual mode.
Landing constraints preserve a full step when it already reaches the target,
otherwise shorten it with a rounding margin and, when needed, split the last
two steps. Intermediate mode tracks the interval in which the accepted state
actually lies, including a state that crossed an output time.

Dense output uses the committed BDF difference columns. The derivative stored
at a new node is the derivative from the solve, which need not be the derivative
of the dense-output polynomial at that node after a step change.

The default norm takes separate relative and absolute-scale references.
The relative reference is the last accepted internal state; the absolute scale
is the last returned state, which can be interpolated. That distinction is
measured in COMSOL and changes its trajectory in interpolated-output mode.
Custom norms remain owned by the caller.

## Validation

Unit tests check polynomial differentiation, convergence, accuracy, restart and
rejection behavior. The examples compile and run finite-element heat conduction,
a first-order oscillator system and a nonlinear reaction system.

The [COMSOL comparison](comsol-alignment.md) uses 21 externally recorded
configurations and 1,336 full-precision samples. It checks each stored time,
state and derivative, failure counts and logged orders. A model-specific
Jacobian-reuse callback and independent exact-solve configurations keep the
algebraic-solver assumptions visible. No reference samples are injected into
the numerical solve.
