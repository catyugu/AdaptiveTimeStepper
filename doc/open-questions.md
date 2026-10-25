# Open questions

What the library does not know yet about COMSOL's time-dependent solver, and the
measurement that would settle each one. The measurements already taken are in
[evidence/](evidence/); what they established about the solver is in
[comsol-alignment.md](comsol-alignment.md).

Nothing here is implemented on a guess: an open question is a place the code
makes a choice the recordings do not separate.

## Open

| ID  | Question                                                                                                                                                                                                                                                               | What settles it                                                                                                                    |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Q07 | The second-order record agrees for eighteen accepted steps and then diverges: the library reduces the step by 0.9 at order 3 where COMSOL keeps its size. The norm is settled; this is one borderline decision, in a pair whose velocity is a hundred times the state. | A record whose velocity varies slowly against the state, or the velocity's `A_s` read out of the solver log, with the solve exact. |

## Settled

Each of these was answered by a measurement, not by reading the manual; the
evidence document named in the answer is the record.

| ID   | Question                                                                                | Answer                                                                                                                                                                                                                                                          |
| ---- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Q01a | Is the scaled absolute tolerance's `beta` observable at all?                            | Yes. Freezing it at 1 fails the unit and the acceptance suites (doc/evidence/tolerance-axes.md).                                                                                                                                                                |
| Q01  | `beta`'s argument for a variable of several dofs at `j > 0`                             | The maximum of `abs(Y)` over the variable, not its RMS. A six-dof field at an automatic scale replays 44 samples to 1.1e-12; the RMS reading agrees for 11 and parts from the record by 1.9e-3 at the 12th (doc/evidence/dofs.md).                              |
| Q02  | Does the auxiliary velocity join the second-order error test?                           | Yes, always, each half with its own weight. The recorded first step is the norm over both halves, and eighteen accepted steps reproduce from the model's own initial values, nineteen from the record's own state (doc/evidence/second-order.md).               |
| Q03  | Does `M` count dependent variables or physics interfaces?                               | Dependent variables. Two global equations in one interface and the same two in an interface each are one trajectory to the last bit, and a partition of unequal size (`(6,1,1)` against `(7,1)`) moves the accepted count from 27 to 28 (doc/evidence/dofs.md). |
| Q04  | Does a rejected attempt count toward "Store every Nth step"?                            | No. `store_every3` rejects a step and its recorded phase does not move (doc/evidence/stepping-sweep.md).                                                                                                                                                        |
| Q05  | Which step does the closest-step storage store, and how many values does a step report? | The first step at or past each requested time, once per step: nine requests on an uneven schedule store seven samples, at the step ends bit for bit (doc/evidence/dofs.md).                                                                                     |
| Q06  | Can a multi-variable tolerance trajectory be replayed exactly?                          | Yes, when the algebraic solve is exact: with `jtech = onevery` the seven- and eight-dof records replay over their whole length, values to 3.3e-13 (doc/evidence/dofs.md). COMSOL's default Jacobian reuse is the caller's business, not the kernel's.           |
| Q08  | Can the condensed second-order step resolve the velocity it derives?                    | No, and it no longer has to: solving the pair makes the velocity an unknown of the caller's step, and the setting a derived velocity cannot meet becomes an ordinary tight tolerance -- 296 steps and `1.6e-7` at `rtol = 1e-8` (doc/evidence/second-order.md). |
