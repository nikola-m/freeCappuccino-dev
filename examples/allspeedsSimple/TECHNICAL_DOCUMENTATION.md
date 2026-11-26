# All-speeds SIMPLE - technical documentation

## 1. Governing equations and unknowns
Ideal gas, `rho = p/(R T)`, `R = Rair`. Unknowns per cell: `u,v,w` (momentum), `p'` (pressure correction), `E` (energy variable of the energy module,
total enthalpy by default), `T` and `rho` from the equation of state. Momentum, energy and thermodynamic properties are those of the existing
code (`calcuvw`, `calc_energy`); this module adds the pressure-velocity coupling and the wall functions.

## 2. One outer iteration (`calc_simple_allspeeds`)
1. `as_time_levels` (energy module): at a new time step shift `rho_o, rho_oo, p_o, E_o` and set BDF coefficients `tc0,tc1,tc2` (steady: 0).
2. `update_wall_viscosity`: wall function -> `visw`.
3. `calcuvw`: momentum predictor, under-relaxation `urfU`, `apu = 1/diag` before relaxation.
4. `calcp_allspeeds`: assemble and solve p', correct fluxes, velocity, pressure, density.
5. in `main`: turbulence, then `calc_energy` (since `compressible = T`), which updates `T`, properties and `rho = p/(R T)`.

## 3. Pressure-correction equation
Face mass flux `m_f = rho_f V_f`, with the Rhie-Chow volume flux (existing `facefluxmass2`, central density `rho_c`)
`V_f = Uf.S - (Vol/ap)_f Df [ (p_N - p_P) - (grad p)_f.d ]`. Face density
`rho_f = beta rho_c + (1-beta) rho_up`, `beta = asimple_beta max(0, 1 - M_f^2)`, `M_f = |V_f|/(|S| c_f)`: central for low Mach numbers
(second order, no decoupling), upwind when `M_f >= 1`.

Perturbations (`T` and `V*` frozen): `m' = rho_f V' + rho'_f V*`, `V' = -(Vol/ap)_f Df (p'_N - p'_P)`, `rho' = p'/(R T)`:

    m'_f = cap (p'_N - p'_P) + cvP p'_P + cvN p'_N ,                cap = -rho_f (Vol/ap)_f Df
    cvP = V* [ beta (1-lam) + (1-beta) H(V*)   ] / (R T_P)
    cvN = V* [ beta lam     + (1-beta) (1-H(V*)) ] / (R T_N)

Cell equation with the BDF density term:

    tc0 Vol/(R T dt) p'_P + sum_f m'_f = - (tc0 rho* - tc1 rho_o + tc2 rho_oo) Vol/dt - sum_f m*_f

Matrix (owner row P, neighbour row N): `a_PN = cap + cvN`, `a_NP = cap - cvP`, `diag_P += -cap + cvP`, `diag_N += -cap - cvN`.
For `M < 1` `cap` dominates (elliptic, M-matrix up to the central part of `cv`); for `M > 1` `cap` is small and the `cv` terms are upwind (hyperbolic).
For `compressible = F` all `cv` and the transient term vanish: classical SIMPLE (verified: identical iteration count with `iccg`).

Steady compressible flow: optional pseudo-time term `(Vol/(R T dtau)) p'` on the diagonal, `dtau = CFL Vol / (0.5 sum_f (|V_f| + c |S_f|))`, CFL ramped geometrically
from 1 over `asimple_cflRamp` iterations. It only affects the p' equation (no momentum pseudo-time; `urfU` plays that role).

Boundaries: see QUICK_REFERENCE. Compressible `outlet`: `m_b = rho_P U_P.S`, diagonal `+ (U_P.S)+/(R T_P)`; compressible `inlet`: `m_b = rho_b U_b.S` with `rho_b = p_b/(R T_b)` (explicit,
the dependency on `p'` is not included: with fixed U it would reduce the diagonal). `pressure`: `facefluxmassPressBnd` / `facefluxmassCorrPressBnd` of the code (use `den(ijp)`).

Correction: `m = m* + m'(p')`, `u = u* + sum_f(-p' S)/ap`, `p = p + asimple_urfP p'` (compressible: no reference pressure; incompressible: `pRefCell`
or zero at pressure patches), `rho = p/(R T)`. Non-orthogonality: `npcor` passes with the existing `fluxmc2` correction.

Consistency check built in (`ltest = T`): the net boundary mass inflow equals `d(mass)/dt` of the BDF term to 4 digits at every iteration (checked on the nozzle).

## 4. AMG for p' (`asimple_amg`)
Setup per solve (matrix changes every iteration, setup cost small): strength `|a_ij| >= theta max_k |a_ik|`; greedy aggregation (seed with its free strong neighbours;
leftovers join the strongest neighbour aggregate); piecewise-constant `P`; Galerkin `A_c = P^T A P` by summing entries into `(agg(i), agg(j))`; V(1,1) cycle with forward/backward
Gauss-Seidel; dense LU below `pamg_coarse_size` unknowns; one cycle is the preconditioner of right-preconditioned flexible GMRES(20). A plain aggregation V-cycle is not used as a
stand-alone solver (piecewise-constant `P` converges slowly). Closed domains give a singular coarse matrix (constant null space): a vanishing pivot is replaced by a huge one, which pins the constant.
The first version without this produced NaN after ~100 iterations of the cavity - found and fixed. Measured: cavity at convergence 5 FGMRES iterations per outer iteration for a 1e-1.5 drop;
NACA dual-time 1.6 iterations on average (tolRel 0.05). No comparison of CPU time with `iccg`/`bicgstab` was made, so no speed claim is made.
`block_amg` is not used or modified; AMG parameter names carry the prefix `pamg_` because `amg_*` belong to `block_amg`.

## 5. Wall functions (`asimple_wallfun`)
For every wall face with tangential velocity `U_t` (relative to the wall) at distance `d_n`, `u_tau` is found by geometric bisection of the law of the wall, `visw = tau_w d_n / U_t`
(the momentum equation uses `max(viscos, visw)`: a wall shear below the laminar one or of the opposite sign - separation, reversed shear - cannot be represented). `ypl`, `tau` are filled.

* `standard`: `tau_w = rho C_mu^1/4 sqrt(k) kappa U_t / ln(E y*)`, `y* = rho C_mu^1/4 sqrt(k) d_n/mu > 11.63`, otherwise `mu U_t/d_n`.
* `reichardt`: `u+ = ln(1+kappa y+)/kappa + 7.8 [1 - exp(-y+/11) - (y+/11) exp(-y+/3)]`, `kappa = 0.41`.
* `spl` - **Shih, Povinelli, Liu, J. Turbulence 4 (2003) N15** (the paper is in the project files), eq. (9)-(11), (16), (20):
  `U_t/u_c = (tau_w/(rho u_tau^2)) (u_tau/u_c) f1(Y_tau+) + (nu/rho) G f2(Y_c+)/u_c^3`, `u_c = u_tau + u_p`, `u_p = ((nu/rho)|G|)^(1/3)`, `Y_tau+ = u_tau d_n/nu`, `Y_c+ = u_c d_n/nu`,
  `G` = tangential part of the cell-centre pressure gradient. `f1`, `f2` are the piecewise polynomial fits of Tables 1 and 2 (viscous, buffer, inertial; `kappa = 0.41`, `C = 5`, `alpha = 5`, `beta = 8`).
  Eq. (16) gives the wall shear vector `tau_w = rho u_tau^2 R/|R|`, `R = U_t - (nu/rho) G f2(Y_c+)/u_c^2`, so that `u_tau f1(Y_tau+) = |R|` closes the system (bisection for `u_tau`).
  This **replaces the earlier `spl`** (Reichardt + clipped adverse-gradient term with constants guessed from a garbled text), which was not the paper's formulation.
  Near zero shear the equation can have several roots (`u_tau = 0` is one when `U_t = u_p f2(Y_p+)`, checked); the largest is taken (attached branch, continuous with the zero-gradient limit).
* `unified` - **Ong & Chan, J. Turbulence 19 (2018) 414** (in the project files): `u_t = u1 + u2`, `u2 = (1/2) u2+ u_p G/|G|` with `(y_p+)^2 = u2+ + exp(-2 beta/alpha) [exp(u2+/alpha) - 1 - u2+/alpha]` (eq. 6, Newton),
  `u_tau` from `y_tau+ = u1+ + exp(-kappa C) [ E - sum_{n<=5} (kappa u1+)^n/n! ]`, `u1+ = |u1|/u_tau`, `E = exp(kappa u1+)` (incompressible) or
  `E = exp((kappa/sqrt(g)) [asin((2 g u1+ - b_q)/Q) - phi])` (isothermal wall, eq. 17), `E = exp((kappa/sqrt(g)) asin(sqrt(g) u1+))` (adiabatic wall, eq. 23);
  `g = r u_tau^2/(2 c_p T_w)`, `r = Pr^1/3`, `b_q = u_tau [(T/T_w - 1)/|u1| + r |u1|/(2 c_p T_w)]` (eq. 35), `Q = sqrt(b_q^2 + 4 g)`, `phi = asin(-b_q/Q)`; the heat-transfer-only limit (`g -> 0`) is `E = exp((2 kappa/b_q)(sqrt(1 + b_q u1+) - 1))`
  (derived from eq. 18, the printed eq. 19 is garbled in the available text). Isothermal wall (`bcDFraction(7,ib) = 1`, Dirichlet T in `0/T`): `T_w` = wall face value; adiabatic: `T_w = T + r u1^2/(2 c_p)`
  (two passes, because the wall properties depend on `u1`). `rho_w = p_w/(R T_w)` with `p_w` extrapolated linearly from the cell, `mu_w` Sutherland (if `sutherland = T`). With `compressible = F`: `g = b_q = 0`, `rho` of the cell.
  `tau_w = rho_w u_tau^2 u1/|u1|` (projected on the velocity direction), `q_w = b_q T_w c_p rho_w u_tau/Pr` into the wall (stored as `qwall`, positive into the fluid).
  Deviations from the paper: (1) bisection on `u_tau` instead of Steffensen's method (the residual is monotone); (2) the derivative of eq. (29) is the analytical one, `1 + exp(-2 beta/alpha)(exp(u2+/alpha) - 1)/alpha`
  - the printed eq. (29) contains an extra term `2 (y_p+)^2/u2+` that is not the derivative of eq. (28); (3) the start of the Newton iteration is `min((y_p+)^2, 2 alpha ln y_p+ + 2 beta + 1)` (right of the root, monotone convergence).
* `compressible`: van Driest transformed Reichardt law with the Crocco-Busemann relation (classical construction, **not** Ong & Chan)
  `T/T_w = rho_w/rho = 1 - B u+ - A^2 u+^2`, `A^2 = r u_tau^2/(2 c_p T_w)`, `r = Pr^1/3`, `B = q_w/(rho_w c_p u_tau T_w)`. With `Bc = -B` the van Driest velocity is
  `u_vD+ = int_0^u+ sqrt(rho/rho_w) du+ = (1/A) [asin((2A^2 u+ - Bc)/S) + asin(Bc/S)]`, `S = sqrt(Bc^2 + 4A^2)`; it is inverted analytically for `u+` (series for `A -> 0`) and `u_vD+` is Reichardt's law in the wall variables
  (`rho_w`, `mu_w`). `q_w = rho_w c_p u_tau (T_w - T_P)/T+` with Kader's law `T+(y+, Pr)`.
* The wall heat flux of `compressible` and `unified` is **not** passed to the energy equation (isothermal walls keep the laminar conduction flux `k Df (T_w - T_P)`).

Verification (this is what was and was not done):
* Stand-alone test of the routines (zero pressure gradient, incompressible, U = 30 m/s, y+ 4.5 ... 900): `u_tau` of Shih et al. / Ong & Chan within 0 ... +6 % / -4 ... +1.6 % of Reichardt (the log constants differ: C = 5.0 vs about 5.5).
  Stratford limit (adverse gradient, `U_t = u_p f2(Y_p+)`): Ong & Chan `u_tau = 0.028 u_p` (zero shear within the accuracy of the fits). Adverse-gradient sweep: `u_tau` decreases monotonically with `G` for both (until `u2` exceeds `U_t` in Ong & Chan,
  where the shear reverses and `visw` is floored). Compressible, `U = 500 m/s`, `T = 170 K`, `p = 40 kPa`, `d_n = 1e-4 m`: adiabatic wall `T_w = 282 K`, `tau_w` 20 % below the incompressible value; isothermal wall hotter than the recovery temperature -> heat flux into the fluid, colder -> into the wall (signs as expected).
* CFD: incompressible channel Re 20000 (README), NACA0012 with y+ < 1 (no change). **No compressible CFD case with a turbulent no-slip wall was run successfully**: the bump with a no-slip wall (Re about 3e7 on an unresolved mesh, impulsive start) diverged for all wall functions and for `none`,
  and the block-coupled solver diverged on the same set-up too, so the cause is the set-up and not the wall function; this was not pursued.

## 6. Time stepping
BDF1 (`bdf = T`) and BDF2 (`bdf2 = T`): momentum from `calcuvw_assemble` (the density of the unsteady momentum term is the current iterate), continuity and energy from the energy module's coefficients.
Dual-time stepping with `maxit` outer iterations per step. BDF2 tested on the NACA case (same dt 1e-4, 1500 steps): Cd 0.05637 against 0.05639 for BDF1. The first step of a new run is BDF1-like for the energy module;
for momentum the module copies `uo` into `uoo` at its first call, because `time_shift` copies the never-initialised `uo` into `uoo` (zero velocities at level n-1): without this copy the BDF2 run gave Cd 0.011 and T down to 249 K.

## 5b. Compressible boundary conditions of the namelist
The kinds of `implicit_bc_compressible.f90` (`stagnation`, `totalpressure`, `massflow`, `compinlet`, `compoutlet`, `open`, `riemann`) are reused through `bc_setup` (once) and `bc_prepare_face`
(every outer iteration, before the momentum predictor: boundary `p_b`, `T_b`, `rho_b`, face mode, Picard under-relaxation `cplbc_urf`). The block-coupled assembly routines of that module are not used; the segregated solver treats the modes as follows:
`MODE_PRESSURE`/`MODE_RIEMANN` - known `p_b` with the Rhie-Chow flux of the boundary density `rho_b = p_b/(R T_b)` (`pressure_bnd_flux_kind`, `p'_b = 0`; the linearised acoustic relations of the `riemann` kind are **not** implicit in the segregated solver, `p_b`, `T_b` are Picard-lagged);
`MODE_EXTRAP` - supersonic exit of a pressure patch, as the `outlet` patch (`rho_P U_P.S` with the upwind `p'` coupling); `MODE_FIXFLUX` (`massflow`) - flux `-mdot` explicit, `U_b` from `mdot/(rho_b A)` after each correction; `compinlet` - the usual inlet with the lagged `rho_b`.
After every correction (`post_bc_kinds`) inflow faces of the pressure-type kinds get the boundary velocity along the inflow direction with the magnitude from the mass flux. Only `public` statements were added to `implicit_bc_compressible.f90` (`mdot_face`, `edir`, `U_set`, `p_set`, `T_set`, `pbc0`, `Tbc0`).
Experience: with a fixed-pressure outlet these kinds need `asimple_urfP` about 0.2 (0.7 diverges at M_in = 0.3); at M_in = 0.5 the default central density (`asimple_beta = 1`) gives a limit cycle, `asimple_beta = 0.5` converges.

## 7. Verification
See README. Beyond it: with `compressible = F` all-speeds SIMPLE and stock SIMPLE take 391 iterations on the cavity with the same solver (`iccg`). The nozzle gives the isentropic solution to 5e-4 in p.
A transient run of the same nozzle with a fixed-velocity inlet and a fixed-pressure outlet (reflective acoustic boundaries, inviscid) does not settle (mass flow oscillates by +-2 kg/s around 8.2)
although continuity is satisfied exactly each iteration (checked via the mass-flow diagnostic); the steady pseudo-time formulation is the one to use there.

## 8. What was changed in the existing files
`main.f90`: `use simple_allspeeds`, branch `elseif ( asimple )`, `calc_energy` also if `asimple .and. compressible`. `read_input.f90`: `use` and namelist entries. `files`, `directories`: new sources.
Nothing else: `velocity.f90`, `pressure.f90`, `block_amg.f90`, energy module untouched. With `ASIMPLE = F` the run is as before.

## 9. Limitations and open points
* Periodic boundaries, Crank-Nicolson, BDF3 with `compressible`: stop with a message.
* Steady NACA M 0.5 on the C-grid (`input-asimple-st.nml`): Cd 0.057-0.060 for several settings but the momentum residual stalls at ~1e-2 (the existing coupled solver has the same
  plateau on this mesh, see `examples/coupledPU`); use the dual-time input. Pseudo CFL 30 with urfP 0.5, urfU 0.5 gave a wrong solution (Mach 0.87).
* Fixed-velocity compressible inlet only robust at low Mach number (M_in = 0.2 works, 0.3 diverges, as for the coupled solver); use `stagnation`, `totalpressure` or `massflow` inlets. Transonic and supersonic flow: bump channel M = 0.5, 0.85, 1.4 (README);
  the transonic case needs `asimple_beta = 0` and a pseudo CFL of 3, and does not converge to machine precision (shock position oscillates by about one cell). Not tested: shocks interacting with walls/boundary layers, 3D, unstructured meshes.
* The heat flux of the wall functions does not feed the energy equation; the compressible wall functions have no compressible CFD verification (section 5). Turbulence is coupled only through `visw` and the existing viscosity modification; no compressible turbulence tests.
* Not parallel (the code is serial).

## 10. Issues noticed in the existing code (not changed)
* `npcor` has no default value: an input file without it gives an undefined number of pressure corrections in the stock SIMPLE (the new module uses `max(npcor,1)`).
* Input, monitor and restart file names are limited to 20 characters; longer names give `End of file` in `read_input` without a hint.
* `continuityErrors` ignores `pressure` patches (the 'global' error of a case with a pressure outlet is the inflow).
* BDF2 first step: `time_shift` copies the never-initialised `uo` into `uoo` before the first step, so BDF2 of the stock solvers starts from zero velocity at level n-1 (workaround in the new module only).
* Concurrent runs in the same case directory overwrite each other's files (`statistics`, monitor, restart).

## 11. References
Przulj V., Basara B., AIAA 2001-2593 (bounded convection schemes); Demirdzic, Lilek, Peric, IJNMF 16 (1993); Ferziger & Peric; Kader 1981 (temperature law);
Reichardt 1951; van Driest 1951; Shih, Povinelli, Liu, J. Turbulence 4 (2003) N15; Ong & Chan, J. Turbulence 19 (2018) 414-430; Darwish, Manganai, Moukalled, Computers & Fluids 168 (2018) 54-66 (BC kinds); Delery / GAMM bump channel (test case description).
