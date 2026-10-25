# Open questions

What the library does not know yet about COMSOL's time-dependent solver, and the
measurement that would settle each one. The measurements already taken are in
[evidence/](evidence/); what they established about the solver is in
[comsol-alignment.md](comsol-alignment.md).

Nothing here is implemented on a guess: an open question is a place the code
makes a choice the recordings do not separate.

## Open

| ID  | Question                                                                                                                                                                                                                                                               | What settles it                                                                                                     |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Q01 | `beta`'s argument for a variable of several degrees of freedom at `j > 0`: the maximum of `abs(Y)` over the variable, or its RMS. The recorded one-variable models cannot tell the two apart, and the argument is `Q01a` below only for `j = 0`.                       | A variable of several dofs whose scale is exact, replayed over its whole series.                                    |
| Q03 | Does `M` count dependent variables or physics interfaces? One variable per interface is what both recorded multi-variable models have, so the two readings coincide there.                                                                                             | A model with two variables in one interface, or one variable in two.                                                |
| Q05 | `Closest` on a non-uniform output schedule: the first step that reaches the time, or the nearest step end. Only a uniform schedule is recorded (`store_closest`).                                                                                                      | One run with `tout = tstepsclosest` and uneven output times.                                                        |
| Q06 | Can a multi-variable tolerance trajectory be replayed exactly? The record carries COMSOL's inexact Newton from the first step size change, so only the first stored step is compared.                                                                                  | Emulating the Jacobian reuse for the field model, as `test/support/comsol_decay_solve.hpp` does for the scalar one. |
| Q07 | The second-order record agrees for eighteen accepted steps and then diverges: the library reduces the step by 0.9 at order 3 where COMSOL keeps its size. The norm is settled; this is one borderline decision, in a pair whose velocity is a hundred times the state. | A record whose velocity varies slowly against the state, or the velocity's `A_s` read out of the solver log.        |

## Settled

Each of these was answered by a measurement, not by reading the manual; the
evidence document named in the answer is the record.

| ID   | Question                                                             | Answer                                                                                                                                                                                                                                                          |
| ---- | -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Q01a | Is the scaled absolute tolerance's `beta` observable at all?         | Yes. Freezing it at 1 fails the unit and the acceptance suites (doc/evidence/tolerance-axes.md).                                                                                                                                                                |
| Q02  | Does the auxiliary velocity join the second-order error test?        | Yes, always, each half with its own weight. The recorded first step is the norm over both halves, and eighteen accepted steps reproduce from the model's own initial values, nineteen from the record's own state (doc/evidence/second-order.md).               |
| Q04  | Does a rejected attempt count toward "Store every Nth step"?         | No. `store_every3` rejects a step and its recorded phase does not move (doc/evidence/stepping-sweep.md).                                                                                                                                                        |
| Q08  | Can the condensed second-order step resolve the velocity it derives? | No, and it no longer has to: solving the pair makes the velocity an unknown of the caller's step, and the setting a derived velocity cannot meet becomes an ordinary tight tolerance -- 296 steps and `1.6e-7` at `rtol = 1e-8` (doc/evidence/second-order.md). |
