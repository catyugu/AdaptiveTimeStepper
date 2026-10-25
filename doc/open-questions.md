# Open questions

What the library does not know yet about COMSOL's time-dependent solver, and the
measurement that would settle each one. The measurements already taken are in
[evidence/](evidence/); what they established about the solver is in
[comsol-alignment.md](comsol-alignment.md).

Nothing here is implemented on a guess: an open question is a place the code
makes a choice the recordings do not separate.

## Open

None. Every rule the kernel implements was decided by a recorded run; the limits
the evidence documents name — the automatic scale's value, a state that mixes the
two orders in time, nonlinear models and other COMSOL versions — are outside the
behaviour the kernel claims, not questions it answers by choice.

## Settled

Each of these was answered by a measurement, not by reading the manual; the
evidence document named in the answer is the record.

| ID  | Question                                                                                                                        | Answer                                                                                                                                                                                                                                                                                                                                                                             |
| --- | ------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Q01 | `beta`'s argument for a variable of several dofs at `j > 0`, and whether the term is observable at all                          | The maximum of `abs(Y)` over the variable, not its RMS; freezing it at 1 fails the unit and the acceptance suites. A six-dof field at an automatic scale replays 44 samples that way, and with the RMS for 11 before parting by 1.9e-3 (doc/evidence/dofs.md).                                                                                                                     |
| Q02 | Does the auxiliary velocity join the second-order error test?                                                                   | No. The test measures the state alone, divided by the two variables the augmented system solves for: the exact record reproduces all 50 samples that way and 20 with the velocity in the norm. The first step cannot tell the two apart, which is why the earlier answer held (doc/evidence/second-order.md).                                                                      |
| Q03 | Does `M` count dependent variables or physics interfaces?                                                                       | Dependent variables, which COMSOL also calls fields: a vector unknown is one of them with several dofs, the library's `Variable::dofs`. Two global equations in one interface and the same two in an interface each are one trajectory to the last bit, and a partition of unequal size (`(6,1,1)` against `(7,1)`) moves the accepted count from 27 to 28 (doc/evidence/dofs.md). |
| Q04 | Does a rejected attempt count toward "Store every Nth step"?                                                                    | No. The recorded phase does not move when a step is refused (doc/evidence/stepping-sweep.md).                                                                                                                                                                                                                                                                                      |
| Q05 | Which step does the closest-step storage store, and how many values does a step report?                                         | The first step at or past each requested time, once per step: nine requests on an uneven schedule store seven samples, at the step ends bit for bit (doc/evidence/dofs.md).                                                                                                                                                                                                        |
| Q06 | Can a multi-variable tolerance trajectory be replayed exactly?                                                                  | Yes, when the algebraic solve is exact: with `jtech = onevery` the multi-variable records replay over their whole length (doc/evidence/dofs.md). COMSOL's default Jacobian reuse is the caller's business, not the kernel's.                                                                                                                                                       |
| Q07 | The second-order record diverged at its twentieth accepted step, a 0.9 reduction where COMSOL kept its size.                    | The error test measured the pair instead of the state. With the state alone over `M = 2` the exact record reproduces all 50 samples; the earlier hypothesis, COMSOL's inexact Newton, was refuted first. Also measured there: an interpolated value's derivative is the dense output's, not the pair's own velocity (doc/evidence/second-order.md).                                |
| Q08 | Can the condensed second-order step resolve the velocity it derives?                                                            | No, and it no longer has to: solving the pair makes the velocity an unknown of the caller's step, and the setting a derived velocity cannot meet becomes an ordinary tight tolerance -- `rtol = 1e-8` takes 261 steps and holds `2.4e-7` (doc/evidence/second-order.md).                                                                                                           |
| Q09 | Does the second-order error test divide by the two variables the augmentation adds, or by two however many variables there are? | By `2K`: `K` independent copies of the recorded oscillator, K = 1, 2 and 3, give the same stored times and states bit for bit, so the divisor grows with the number of variables (doc/evidence/second-order.md).                                                                                                                                                                   |
