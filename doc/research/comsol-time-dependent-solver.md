# COMSOL 6.x adaptive BDF time stepping — documented facts

Documentation research only. Every claim below is a quote or a direct reading of the released
COMSOL 6.x Reference Manual; page numbers are PDF page numbers of the cited PDF and are written
into the URL as `#page=N` so each one opens on the page being quoted.

## 0. Source situation (read this first)

- The HTML chapter pages do not exist: `https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/comsol_ref_solver.24.html`
  and the neighbouring `.23`/`.25` names return **HTTP 404** for 5.6, 6.0, 6.1 and 6.2 (also
  `comsol_ref_solver.html`, `comsol_ref.html`, `index.html`). `https://doc.comsol.com/` now redirects to
  `https://doc.comsol.com/6.4/docserver/`, a Vaadin single-page application (`docserver.js`, no static
  per-chapter HTML), so there is no quotable chapter URL for 6.x.
- Exception: **6.4 does serve per-chapter HTML** under `.../6.4/doc/com.comsol.help.comsol/comsol_ref_solver.<n>.<m>.html`,
  and `comsol_ref_solver.36.122.html` ("About the Time-Dependent Solver") returns HTTP 200 with the
  second-order passage of §9 verbatim.
- What *is* published and fetchable: the full Reference Manual PDFs, `.../doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf`,
  for 6.0, 6.1, 6.2, 6.3, 6.4 (HTTP 200, ~48 MB). All quotes below come from the **6.2** PDF unless marked,
  and the load-bearing passages were re-checked in the **6.0** and **6.4** PDFs.
- Drift between 6.0 / 6.2 / 6.4: the settings, their names and their defaults are stable across the three
  releases; the *only* substantive change found is the BDF initial-step rule, which 6.4 states more
  precisely than 6.0/6.2 (see §4).
- COMSOL writes the SUNDIALS BDF solver as "**IDA**" in the theory text but as "**IDAS**" in two
  sentences of the same section ("the built-in nonlinear solver IDAS", "the BDF (IDAS) solver"). This is
  a documentation inconsistency, not a second solver.

## 1. Where the settings live (study step vs. solver node)

- Study step **Time Dependent** (`Study > Time Dependent`) — `RM62` p. 1490–1491. Sections:
  `Study Settings` (Time unit; `Times`; `Tolerance` = `Physics controlled` (default) / `User controlled`
  + `Relative tolerance`), `Results While Solving`, `Adaptation`, `Study Extensions`,
  plus the common sections. There is **no** BDF/step-size setting on the study step.
- Solver node **Time-Dependent Solver** (`Solver Configurations > Time-Dependent Solver`) — `RM62`
  p. 1616–1624. Sections: `General` (Defined by study step, Time unit, Output times, Times to store,
  Relative tolerance), `Absolute Tolerance`, `Time Stepping (General Settings)`,
  `Time Stepping (BDF and Generalized Alpha)`, `Time Stepping (Runge–Kutta Methods)`,
  `Algebraic Variable Settings`, `Output`, `Advanced`, `Least-Squares Data`.
- Synchronisation: "Use the **Defined by study step** list to specify if the settings are synchronized
  with the corresponding study step (the default). If you select **User defined** ... you can specify the
  **Time unit**, **Output times**, and **Relative tolerance** settings." (p. 1616)
- The study step field is called `Times`; the solver-node field is called `Output times`; the
  `Steps taken by solver` bullets refer to "the **Times** field in the General section" (p. 1620).

## 2. Solver type and method list (`Time Stepping (General Settings)`, p. 1619)

- "From the **Solver type** list, choose **Implicit** (the default) to use an implicit time-stepping method
  or choose **Explicit** to use an explicit time-stepping method."
- Implicit `Method`: **BDF** ("to use a backward differentiation formula"), **Generalized alpha**
  ("to use the generalized-α method"), **Initialization only** ("to compute consistent initial values only
  and then stop").
- Explicit `Method`: **Runge–Kutta** (sub-list `RK34`, `Cash–Karp 5`, `Dormand–Prince 5`),
  **Classical Runge–Kutta** ("order 1–4 ... Runge–Kutta 4 is the default choice for the discontinuous
  Galerkin method"), **Adams–Bashforth 3**, **Adams–Bashforth 3 (local)**.
- "The **Generalized alpha** time-stepping method requires a Fully Coupled or Segregated attribute node.
  The BDF time-stepping method can be used without a Fully Coupled or Segregated attribute node. In such
  a situation, the BDF method uses an internal automatic nonlinear solver." (p. 1620)

## 3. `Steps taken by solver` — the exact four values (p. 1620)

Quoted verbatim, in list order:

- "**Free** to let the time-stepping method choose time steps freely. The times specified in the **Times**
  field in the General section are not considered when a time step is chosen."
- "**Intermediate** to force the time-stepping method to take at least one step in each subinterval of the
  times specified in the **Times** field in the General section."
- "**Strict** to force the time-stepping method to take steps that end at the times specified in the
  **Times** field in the General section. The solver takes additional steps in between these times if
  necessary."
- "**Manual** to override the automatic choice of time step with a manual choice. Manual time stepping can
  be useful in cases where the automatic time-step selection does not work; for example, in contact
  problems, rotating machinery, or fluid-structure interaction."
- Note attached to `Manual`: "**Manual** is available for BDF and Generalized alpha and overrides the
  local error estimation made in each time step." (p. 1620)

So "Intermediate" = at least one step per output subinterval, not necessarily landing on the times;
"Strict" = every output time is a step end point.

## 4. `Initial step` — default is automatic, and it is *not* purely geometric

- 6.2 (p. 1620): "**Initial step.** By default, the solver chooses an initial step automatically. Select the
  **Initial step** check box for manual specification of an initial step. **By default, the first step is
  0.1% of the end time**, which can affect your solution so that you do not get the same time-step history
  up to a certain earlier time when the end time is changed."
- Same paragraph continues: "Also, the results of the consistent initialization can be strongly dependent on
  the initial time step taken, if the initial conditions and the boundary conditions do not match".
- 6.4 (pp. 1779–1780) splits the rule per method and adds the derivative condition — this is the most
  important release-to-release change:
  - "By default, for generalized-α, the first step is 0.1% of the end time, ..."
  - "**By default, for BDF, the first step is less than 0.1% of the end time but also fulfills
    DT·ypnorm <= 1/2, where DT is the initial time step and ypnorm is the weighted norm of the time
    derivative vector.**"
- Consequence for a reimplementation: the automatic initial step is `min(1e-3 * t_end, DT)` with `DT`
  fixed by `DT * ||u'(t0)||_w <= 0.5` (weighted norm, i.e. the same scaling as the error test in §10).
  6.0 and 6.2 only state the `0.1% of the end time` half.

## 5. `Maximum step constraint` — automatic default, no documented formula (pp. 1620–1621)

- "**Maximum step constraint.** By default, the solver chooses a maximum time step automatically. Select
  **Constant** as the maximum step constraint for manual specification of a fixed maximum time step. A
  constant maximum step constraint is a positive scalar value, which can be an expression that evaluates to
  a numerical value before entering the solver. ... Select **Expression** as the maximum step constraint for
  more general expressions of the allowed maximum time step. These expressions are evaluated while solving
  and can, for instance, depend on the time parameter itself."
- "The solver uses the absolute value of the expression for the maximum step constraint."
- The value computed by the *automatic* setting is **not documented anywhere in the fetched docs** (6.0,
  6.2, 6.4). The only "10 %" rule in the manuals belongs to a different solver: "The automatic maximum step
  size is 10% of the total simulation time for time-dependent flows ..." — that sentence sits in the
  *Particle Tracing* ODE-solver settings (RM62 p. 1945), not in the FEM Time-Dependent Solver.

## 6. BDF order settings (`Additionally, if BDF is selected as the time stepping Method`, pp. 1620–1621)

- "**Maximum BDF order** (available if **Free**, **Intermediate**, or **Strict** is selected for the **Steps
  taken by solver**). This setting controls the maximum allowed degree of the interpolating polynomial of
  the BDF method."
- "**Minimum BDF order** ... This setting can be used to prevent the solver from decreasing the order of the
  BDF method below 2."
- "**BDF order** (available if **Manual** is selected ...). **The order of the BDF can be 1–5 (default order:
  2).**"
- "**Initial step fraction** (available if **Manual** ...). During the startup of the BDF method, a shorter
  time step will be used to compensate for the lower order that is used for the first handful of steps. The
  initial step is a fraction of the first step, and the solver then exponentially increases the step length
  until the requested step length is reached. This settings and the initial step growth rate below control
  that startup phase. The default values depend on the selected BDF order."
- "**Initial step growth rate** (available if **Manual** ...)." — no numeric default given.
- Defaults of `Maximum BDF order` / `Minimum BDF order` are **not stated numerically**; they are
  physics-controlled (see §8). `BDF order = 2` is the default only in Manual mode.

## 7. Manual (fixed-step) time stepping (pp. 1620–1621)

- "**Time step** (available if **Manual** is selected for the **Steps taken by solver**). Enter a manual time
  step specification as a scalar, a vector of times, or an expression containing global variables or
  parameters in the **Time step** field using the base unit for time. **The relative and absolute tolerances
  are still used to terminate the algebraic equations at each time step. Also, the requested time step will
  be reduced if the algebraic solver does not converge.**"
- Manual overrides the local error estimation (§3), but is not unconditionally rigid: it is cut when the
  nonlinear solve fails.
- `Time step increase delay` (Free/Intermediate/Strict): "enter a positive integer ... the number of time
  steps taken by the solver before the increase of the time step is actually performed, from the first step
  where the error estimator signals that the current step is too small. ... Entering 0 results in the same
  behavior as clearing the check box." (p. 1621)
- `Nonlinear controller` (BDF): "Select this check box to use a nonlinear controller for more efficient
  time-step control in the BDF method, especially for highly nonlinear problems such as multiphase flow and
  turbulence in fluid dynamics. When nonlinear failures occur, the nonlinear controller becomes active and
  uses a more careful time step control. The nonlinear controller acknowledges that the step size for Newton
  stability might be smaller than the step size for BDF accuracy." (p. 1621)

## 8. Which physics defaults to which method/order

- Heat transfer (Reference Manual, `CHAPTER 14`, p. 1238, "Time-Dependent Study Step — Default Settings for
  Heat Transfer Interfaces"): "**The default time-stepping method for the Heat Transfer interfaces is BDF at
  second order. It excludes algebraic variable from the error estimate.**"
- Laminar/fluid flow (Reference Manual, `CHAPTER 13`, p. 1165, "Time-Dependent Solvers"): "**The default
  time-dependent solver for Navier–Stokes is the BDF method with maximum order set to two.** Higher BDF
  orders are not stable for transport problems in general nor for Navier–Stokes in particular."
- Structural mechanics (Structural Mechanics Module User's Guide, p. 69): "**The Generalized alpha method,
  which is recommended for structural dynamics problems. This is the default solver if Structural Transient
  Behavior is set to Include inertial terms in the physics interface settings.**" / "**The BDF method, which
  is recommended for first-order problems. This is the default solver if Structural Transient Behavior is
  set to Quasistatic in the physics interface settings.**" Also p. 238: "the BDF solvers should be avoided
  for dynamic contact analyses. The default solver suggests the generalized-α solver when inertial terms are
  included"; and "The default solver for the time stepping is changed from Generalized alpha to BDF when
  Quasistatic is selected."
- Electromagnetic waves (Reference Manual, p. 1583): "some physics interfaces in COMSOL Multiphysics — for
  structural mechanics and electromagnetic waves, for example — use generalized-α as the default transient
  solver."
- Electric currents: **not found in the fetched docs.** The AC/DC Module User's Guide (6.2, 501 pages) was
  downloaded and searched for `BDF`, `time-stepping`, `generalized alpha` — no statement on the default
  transient method or BDF order. Its only related sentences are about formulation (quasistatic/transient
  continuity equations), not about the solver.

## 9. The BDF method as documented (`The Implicit Time-Dependent Solver Algorithms`, pp. 1580–1583)

- "In COMSOL Multiphysics, the **IDA** and generalized-α solvers are available to solve the above ODE or DAE
  system:" (p. 1581)
- "**IDA** was created at the Lawrence Livermore National Laboratory (Ref. 4) and is a modernized
  implementation of the DAE solver **DASPK** (Ref. 5), which uses **variable-order variable-step-size
  backward differentiation formulas (BDF)**. Optionally, a nonlinear controller can provide more careful
  time-step control, which can reduce the thrashing of time steps that increases the number of time steps
  and degrades performance for highly nonlinear problems in CFD, for example. See Ref. 7 for more
  information about the nonlinear controller (STAB controller)." (p. 1581)
- References (p. 1629 list, "References for the Solution Operation Nodes and Solvers"): Ref. 4 = Hindmarsh
  et al., "SUNDIALS: Suite of Nonlinear and Differential/Algebraic Equation Solvers", ACM TOMS 31 (2005);
  Ref. 5 = Brown, Hindmarsh, Petzold, SIAM J. Sci. Comput. 15 (1994); Ref. 6 = Brown/Hindmarsh/Petzold,
  "Consistent Initial Condition Calculation for Differential-Algebraic Systems", SIAM J. Sci. Comput. 19
  (1998); Ref. 7 = Moore & Petzold, "A Stepsize Control Strategy for Stiff Systems of ODEs", Appl. Numer.
  Math. 15 (1994).
- "The BDF solver uses backward differentiation formulas with **an order of accuracy varying from one (that
  is, backward Euler) to five**." (p. 1583) — consistent with `BDF order ... 1–5` in the settings.
- Damping character: "BDF methods have been used for a long time and are known for their stability. However,
  they can have severe damping effects, especially the lower-order methods. Backward Euler severely damps
  any high frequencies. Even if you are expecting a solution with sharp gradients, you might get a very
  smooth solution due to the damping in the backward Euler method." (p. 1583)
- "The generalized-α (alpha) solver has properties similar to the second-order BDF solver but the underlying
  technology is different. ... Compared to BDF (with maximum order two), generalized-α causes much less
  damping and is thereby more accurate. For the same reason it is also less stable." / "the BDF method is
  more robust and better at handling changes to the time step." (p. 1583)
- Second-order-in-time systems under IDA: "extra variables are internally introduced so that it is possible
  to form a first-order time-derivative system ... This procedure expands the original ODE or DAE system to
  double its original size, but the linearized system is reduced to the original size with the matrix
  E + σD + σ²K, where σ is a scalar inversely proportional to the time step. ... **The error test excludes
  the variable Uv unless consistent initialization is on**, in which case the differential Uv-variables are
  included in the error test and the error estimation strategy applies to the algebraic Uv-variables."
  (pp. 1581–1582)
- The 6.4 HTML page ([`comsol_ref_solver.36.122.html`](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html))
  states the same passage with the two parentheticals a library author needs: the augmentation happens
  "**(this does not happen when using generalized-α because it can integrate second-order equations)**",
  and "By the added equation, the original variable U is therefore always a differential variable
  (index-0)." On the mixed system: "The implementation of the generalized-α method in COMSOL Multiphysics
  detects which variables are first order in time and which variables are second order in time and applies
  the correct formulas to the variables."
- The documented statement about what selects the time step is only the error test of §10 plus the
  test in §5/§7. **No safety factor for the BDF step-size update is documented**, and no explicit
  `h_new = h * f(err)` formula is given for BDF. (`Step size safety factor`, default 0.9, and the PI
  controller with `h_{n+1} = S h_n err_n^α err_{n-1}^β` are documented only for the explicit
  Dormand–Prince 5 controller, p. 1623.)

## 10. Tolerances and the error test (`Absolute Tolerance Settings for the Time-Dependent Solver`, p. 1582)

Relative tolerance:

- "Use the **Relative tolerance** field to enter a positive scalar number (**default: 0.01**). The solver
  uses this number to control the relative error in each time step." (solver node, p. 1617)
- Study step: "From the **Tolerance** list, choose **Physics controlled** (the default) to use the tolerance
  suggested by the physics. Choose **User controlled** to override the suggested relative tolerance with a
  value that you enter in the **Relative tolerance** field." (p. 1491)

Absolute tolerance (`Absolute Tolerance` section, pp. 1617–1618):

- "**Global method**: **Scaled** to let the absolute tolerance be applied to scaled variables. **Unscaled** to
  let the absolute tolerance be applied to unscaled variables."
- "**Tolerance method**: **Factor** to specify a factor (**default: 0.1**) in the **Tolerance factor** field
  that makes the absolute tolerance proportional to the relative tolerance. With this method, the absolute
  error does not become dominating when you reduce the relative tolerance so that it is much smaller than
  the absolute tolerance. **Manual** to specify a value for the absolute tolerance."
- "**Derivative tolerance method**: **Automatic** to use an automatic derivative tolerance. The automatic
  method sets the time derivative factor (see below) to **one over the time span**. **Factor** to specify a
  factor (**default: 1 Hz**) in the **Time derivative factor** field. The time derivative factor acts
  between the tolerance on the variables and their time derivatives."
- Per-variable overrides: `Method` = `Scaled` / `Unscaled` / `Use global` (default), each with its own
  `Tolerance method`, `Tolerance factor` (0.1), `Tolerance`, `Tolerance for time derivatives`.

The accepted-step error test, **Unscaled** method (exact formula, p. 1582):

```text
( (1/M) * sum_j (1/N_j) * sum_i ( |E_i| / (A_us,i + R*|U_i|) )^2 )^(1/2)  <  1
```

- "where `A_us,i` is the unscaled absolute tolerance for DOF i, `R` is the relative tolerance, `M` is the
  number of fields, and `N_j` is the number of degrees of freedom in field j. ... for vector elements there
  is a field-to-DOF conversion factor involved."
- **Scaled** method (and `Update scaled absolute tolerance` selected):

```text
( (1/M) * sum_j (1/N_j) * sum_i ( |E_Yi| / (A_s,i + R*|Y_i|) )^2 )^(1/2)  <  1
```

- with `A_s,i = A_ki * (beta + ||Y||_{2,i})`, `beta = (1 - exp(-alpha*j)) * ||Y||_{inf,ki} + exp(-alpha*j)`
  when `0 < ||Y||_{inf,ki} < 1`, else `beta = 1`, `alpha = 1/5`, `j` = time-step iteration number
  `j = 0, 1, ...`; `||Y||_{2,ki}` and `||Y||_{inf,ki}` are the 2-norm and max norm of dependent variable
  `k_i`. For another scaling method or when the check box is cleared, `A_s,i = A_ki`.
- So the error weight per DOF is `A + R*|U_i|` (scaled or unscaled variant), and the norm is a
  **field-and-DOF-averaged weighted RMS** (not an infinity norm) required to be `< 1`.
- Warnings attached to this section: "**If the solution is smaller than the absolute tolerance, there is no
  accuracy at all.**" and "Note that with the **Tolerance method** set to **Factor**, the absolute tolerance
  is specified as a factor of the relative tolerance. If you use a specific absolute tolerance instead, that
  value can reduce the effect of further reductions of the relative tolerance (so that when `R|U_i| ≈ A_i`,
  the effect of further reductions of the relative tolerance diminishes)."
- DAE handling (`Algebraic Variable Settings`, p. 1622): "**Error estimation**: **Include algebraic** (the
  default) to include the algebraic degrees of freedom in the error norm. **Exclude algebraic** to exclude
  the algebraic degrees of freedom from the error norm. ... In general, excluding algebraic degrees of freedom
  is not recommended when solving DAE systems of index 1, whereas it can be generally encouraged for DAE
  systems of index 2."
- Also `Singular mass matrix`: `Maybe` ("look for zero-filled rows or columns in the mass matrix as a means of
  detecting a differential-algebraic equation") or `Yes, if ...` (p. 1622).

## 11. Consistent initialization / initialization of the time stepping (p. 1622)

- "Use the **Consistent initialization** list to control how the solver performs consistent initialization of
  differential-algebraic systems. The consistent initialization procedure tries to reconcile inconsistent
  initial values. In doing that, **it takes a small time step (by default, 0.1% of the initial step size)**."
- "**Backward Euler (the default)** to perform consistent initialization using a small **artificial step**
  with the backward Euler method. When this is selected, enter a value in the **Fraction of initial step for
  Backward Euler** field. This value is a dimensionless quantity that determines the size of the time step
  for the backward Euler method (in terms of the initial step). ... **The default value is 0.001 (that is,
  the small backward Euler step size is 0.1% of the initial step size).**"
- "Also, the **Safety factor for Backward Euler** is used in the algebraic equation termination during
  consistent initialization using the Backward Euler method. A larger value gives a stricter condition.
  **The default value of 20 corresponds to the normal behavior for any new time step.**"
- `Off`: "indicate that the initial values already are consistent, which means that the solver does not
  modify them." `On`: "use a consistent initialization routine that is preferable to Backward Euler for
  index-1 differential-algebraic equations. This method detect the algebraic and differential variables from
  the sparsity structure of the damping and mass matrices D and E (see Ref. 6)."
- `Rescale after initialization` (only with Backward Euler): "select ... if you want to perform a rescaling
  of the variables. In fluid flow simulations, for example, the scaling for velocity can be problematic
  because the velocity is often zero at the start but then changes dramatically after a consistent
  initialization. Rescaling of the velocity and pressure can provide a smoother start and **avoid using very
  small time steps**." (pp. 1622–1623)
- The BDF startup phase (lower order for the first steps, exponential growth to the requested step) is a
  separate mechanism, described in §6; the artificial backward-Euler step of §11 is used for consistency,
  not for accuracy.

## 12. Events (implicit / explicit) and their interaction with the step

Interface (Reference Manual, `The Events Interface`, pp. 1351–1354):

- "The **Events** (ev) interface ... is used to create solver events. An event can be explicit or implicit,
  and the difference is that for explicit events you must specify the exact time when the event occurs.
  **When an event occurs, the solver stops and provides a possibility to reinitialize the values of states
  and dependent variables.**" (p. 1351)
- Implicit event conditions: "When the condition changes its state from false to true at the zero crossing of
  one of the indicator states, this triggers the implicit event. ... The condition is only evaluated at zero
  crossings of the indicator states, so, for example, using z>1 as a condition never triggers the event."
  (p. 1353; indicator states satisfy `nojac(g(v,vt,vtt,t)) - u = 0`)
- `Event tolerance` (solver, p. 1621): "**Event tolerance.** This setting can be used to set the event
  tolerance (**default value: 0.01**), which is used for **root finding of event conditions** when using
  implicit events; see Explicit Event." — with the note "This does not apply to the startup phase of the
  simulation."
- Reinitialization process (p. 1353), verbatim:
  1. "Dependent variables are reinitialized sequentially. If there are N triggered events, then N
     reinitializations are carried out."
  2. "The reinitializations carried out in step 1 may have led to nonconsistent dependent variables ...
     By default, the solver will try to correct the reinitialized state to enforce consistency, **using a
     small artificial time step with the backward Euler method**. This can be bypassed by clearing the
     **Use consistent initialization** check box".
  3. "If some implicit events have been triggered by steps 1–2, then steps 1–2 are carried out again, to
     process these events. This is indicated in the solver log."
- "Even if several events (explicit or implicit) are triggered simultaneously, consistency is enforced at
  most once in step 2."
- The consistency method after reinitialization is taken from the solver node's `Consistent initialization`
  list: "If set to **Backward Euler**, then dependent variables are corrected with the backward Euler method
  using a small artificial time step. If set to **Off**, then consistency is not enforced after
  reinitialization. If set to **On** and a Fully Coupled or Segregated solver subnode is active, then
  dependent variables are corrected with a backward Euler step. If set to **On** and the built-in nonlinear
  IDAS solver is used, then consistency is enforced using either IDAS or a backward Euler step, depending on
  the equation solved." (p. 1354)
- Output around events: `Store solution before and after events` — "store two additional solutions every
  time an implicit or explicit event is triggered ... This stores the solutions before and after the
  reinitialization." (p. 1624, `Output` section). For the orbit study steps: "Additional output times will be
  stored immediately before and after every explicit and implicit events" (p. 1541).
- **Not documented**: how event location limits or restarts the step (step rejection/rollback to the event
  time, restart order, first-step-after-event size). The docs only give the root-finding tolerance, the
  stop-and-reinitialize procedure and the two extra stored solutions.

## 13. Output times, stored steps and interpolation

`General` section of the solver node (p. 1617):

- `Output times`: "enter a vector of times that define the time span for the simulation using the Range
  button if needed (**default: range(0,0.1,1)**)."
- `Times to store`, exactly three options:
  - "**Output times by interpolation** (the default) to store the output times using interpolation from the
    solver's time steps."
  - "**Steps taken by solver closest to output times** to use the time steps from the solver that are
    closest to a given output time among all the steps that have passed this output time."
  - "**Steps taken by solver** to store the actual time steps from the solver. Specify to store every Nth
    step in the **Store every Nth step** field (**default value: 1**; that is, all time steps are stored)."
- (The `Output` section, p. 1624, has a second, narrower list also called `Times to store`: `Specified
  values` — "store solutions at the values entered in the **Times** field in the General section" — or
  `Steps taken by solver`.)
- `Interpolate solution at end time` (p. 1619): "For all settings for **Steps taken by solver** except
  **Strict** for the implicit solvers and for all explicit solvers, the **Interpolate solution at end time**
  check box is selected by default. **The solver can then step past the last time in the Output times list,
  leading to the end time being interpolated.** If you clear this check box, the solver will not step past
  the last time in the Output times list, so that the time stepping includes the last time and takes no steps
  past that time."
- Interpolation used for postprocessing between stored times (study step, p. 1491): "When plotting the
  results from a time-dependent simulation, you can choose to plot the solution at any of the times specified
  in the **Output times** field. You can also plot an interpolated solution at any intermediate time. **The
  interpolation used between times is a cubic Hermite spline; that is, the interpolation uses both the
  solution values and their time derivatives at two points: the closest output times before and after the
  time for which the interpolated solution is computed.**"
- "The tolerance settings control the internal time steps taken by the solver, so **selecting large time
  steps for the output times does not affect the accuracy in the time stepping**." (p. 1491)
- Side effect of interpolated storage: "For interpolated output (the **Times to store** list set to
  **Specified values**), no data is stored for the residual (NaN)." (p. 320 / 6.4 p. 354)
- **Not documented**: the order of the polynomial the solver itself uses when it produces `Output times by
  interpolation` (only the cubic-Hermite rule for plotting/intermediate evaluation is stated).

## 14. Documented statements about number of steps and performance

- STAB/nonlinear controller (p. 1581): the controller "can reduce the **thrashing of time steps that
  increases the number of time steps and degrades performance** for highly nonlinear problems in CFD".
- Heat transfer tuning section (p. 1238): "The quality of the time-stepping influences the nonlinear solver
  convergence. **Tiny time steps usually lead to mildly nonlinear problems at each time step whereas large
  time steps can result in (fewer) highly nonlinear problems.**" and "The default solver settings for
  transient heat transfer defines the maximal number of nonlinear iterations to 5."
- Heat transfer, on controlling the step (p. 1238): lower the relative tolerance (and the absolute tolerance
  in the same proportion), or set a maximum time step, or "include times of interest in the **Times** field
  of the time-dependent study and to use the **Intermediate** option", or trigger an event.
- Fluid flow (p. 1165): "Higher BDF orders are not stable for transport problems in general nor for
  Navier–Stokes in particular"; "Both BDF and generalized-α are per default set to automatically adjust the
  time step."
- Structural mechanics / dynamic contact (SME guide p. 238): "If you rely on the (default) adaptive
  time-stepping algorithm, the solver typically also tries to resolve the wave propagation in the domains
  adjacent to the contact pairs. **This can cause unnecessarily small time steps and increase the computation
  time.** To avoid this, you can modify the solver to use a manual time-stepping algorithm".
- No documented step-count formula, no "18 steps" statement, and **"deadbeat" does not occur anywhere** in
  the 6.0, 6.2 or 6.4 Reference Manuals, in the AC/DC Module User's Guide, or in the Structural Mechanics
  Module User's Guide (full-text search of the extracted PDF text).

## 15. Facts requested but not found in the fetched docs

- `Maximum step constraint = Automatic`: the value/formula it computes — **not found**.
- `Maximum BDF order` / `Minimum BDF order` numeric defaults (per physics or global) — **not found**; only
  the heat-transfer ("BDF at second order") and Navier–Stokes ("maximum order set to two") statements.
- Electric currents default time-stepping method / BDF order — **not found** (AC/DC Module User's Guide 6.2
  has no statement).
- BDF step-size update law, safety factor for BDF, "local truncation error below tolerance" phrasing —
  **not found**; what exists is the accept/reject inequality of §10 and the STAB-controller reference.
- "18 steps", "deadbeat" — **not found** (see §14).
- Internal interpolation order for `Output times by interpolation` — **not found** (only cubic Hermite for
  plotting, §13).
- Physics default solver for `Electric Currents`, and per-physics `Relative tolerance` ("Physics controlled")
  numeric values — **not found**.

## Sources

- COMSOL Multiphysics Reference Manual 6.2 (chapter 20 "Studies and Solvers", chapter 16 "Equation-Based
  Modeling"), PDF, fetched 2026-10:
  [full PDF](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf)
  - [p. 1490 — Time Dependent study step](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1490)
  - [p. 1491 — cubic Hermite interpolation, Tolerance list](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1491)
  - [p. 1580 — The Implicit Time-Dependent Solver Algorithms](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1580)
  - [p. 1581 — IDA / DASPK / SUNDIALS / STAB controller / IDAS](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1581)
  - [p. 1582 — Absolute Tolerance Settings, error-test formulas](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1582)
  - [p. 1583 — BDF vs. generalized-α and explicit methods](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1583)
  - [p. 1616 — Time-Dependent Solver, General](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1616)
  - [p. 1617 — Output times, Times to store, Relative tolerance](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1617)
  - [p. 1618 — Absolute Tolerance section](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1618)
  - [p. 1619 — Time Stepping (General Settings)](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1619)
  - [p. 1620 — Steps taken by solver, Initial step, Maximum step constraint](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1620)
  - [p. 1621 — BDF order, Manual, Event tolerance, Nonlinear controller](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1621)
  - [p. 1622 — Consistent initialization, Error estimation](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1622)
  - [p. 1623 — Interpolate solution at end time (RK), step growth ratios](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1623)
  - [p. 1624 — Output section (Times to store, store before/after events)](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1624)
  - [p. 1629 — References for the Solution Operation Nodes and Solvers](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1629)
  - [p. 1351 — The Events Interface](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1351)
  - [p. 1353 — Implicit Event, About the Reinitialization Process](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1353)
  - [p. 1354 — Enforcement of consistency after reinitialization](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1354)
  - [p. 1238 — default heat transfer method: BDF at second order](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1238)
  - [p. 1165 — Navier–Stokes: BDF with maximum order two](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1165)
  - [p. 1945 — particle-tracing ODE solver: 10 % automatic maximum step](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1945)
- Reference Manual 6.0 (for older wording):
  [p. 1485](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1485),
  [p. 1451](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1451),
  [p. 1450](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1450)
- Reference Manual 6.4 (for the refined BDF initial-step rule):
  [p. 1779](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1779),
  [p. 1780](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1780),
  [p. 1727](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1727),
  [p. 1728](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf#page=1728)
- Structural Mechanics Module User's Guide 6.2:
  [p. 69 — generalized-α vs BDF default per Structural Transient Behavior](https://doc.comsol.com/6.2/doc/com.comsol.help.sme/StructuralMechanicsModuleUsersGuide.pdf#page=69),
  [p. 238 — dynamic contact, manual stepping, numerical damping](https://doc.comsol.com/6.2/doc/com.comsol.help.sme/StructuralMechanicsModuleUsersGuide.pdf#page=238)
- AC/DC Module User's Guide 6.2 (searched; no statement on the transient solver):
  [full PDF](https://doc.comsol.com/6.2/doc/com.comsol.help.acdc/ACDCModuleUsersGuide.pdf)
- Starting points that 404 (documented here so the dead ends are not repeated):
  `https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/comsol_ref_solver.24.html` (and `.23`, `.25`,
  `comsol_ref_solver.html`, `index.html`, `helpindex.html`) → HTTP 404;
  `https://doc.comsol.com/` → redirect to `https://doc.comsol.com/6.4/docserver/` (Vaadin SPA).
