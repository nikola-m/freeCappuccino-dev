# Technical documentation - turbulence models, patches 0031-0040

## 1. Architecture

### 1.1 Model registry (`turbulence_registry.f90`)

Before: the model name was compared in three places (flag setter in `TurbModelData`, `modify_viscosity_turbulence`, `modify_viscosity_inlet`).
Now `select_turbulence_model` (called from `read_input` after the namelist) points two procedure pointers, `turb_update` and
`turb_update_inlet`, to the model routines and sets `solveTKE/solveEpsilon/solveOmega/solveNuTilda` and `lasm`. The main loop calls
`turb_update`. An unknown name stops the run with a message. A new model is one `case` block, one entry in `src/cappuccino/files`
(module before its users), nothing else. Existing models are registered unchanged; a restart of `k_omega_sst` is bitwise identical
before/after the refactor.

`TurbModel%Scalar` was enlarged from 4 to 10 slots (urf, gds, cScheme, lSolver, ... per transported scalar); `TurbModel%name` is 40 characters
(`easm_rumsey_gatski_komega` was truncated at 20).

### 1.2 Scalar transport engine (`turb_transport.f90`)

Every turbulence transport equation has the form
`d(rho phi)/dt + div(rho U phi) - div(Gamma grad phi) = Su - Sp phi`. A model only evaluates the volume-integrated sources `suv`
(explicit), `spv` (implicit, >= 0) and the diffusivity field `dif` (cells and boundary faces) and calls

```
transport_scalar( fi, dfidxi, fio, fioo, ifi, fname, suv, spv, dif, lconv, ltrans, wallDirichlet, fixed, fixval, positive )
```

The engine does gradient, convection with any scheme of `facefluxsc`, diffusion with non-orthogonal correction, BDF/BDF2/CN, inlet/outlet/
pressure/wall/periodic boundaries, under-relaxation, linear solve, boundary update, clipping. Options: `lconv=.false., ltrans=.false.`
(purely elliptic: `f` of zeta-f, `alpha` of EB-RSM); `fixed/fixval` (rows replaced by identity: wall-function epsilon/omega, wall-cell stresses);
`wallDirichlet` (walls hold a Dirichlet value with diffusive flux, otherwise zero flux); `positive=.false.` (stresses). Helpers:
`shift_time_levels` (time levels of model-owned fields, detects a new step by `itime`), `wall_normal_tangent`, `wall_viscosity_loglaw`,
`copy_owner_to_boundary`, `report_range`. Caution: `grad` writes boundary values, so gradient arrays of engine fields must be `(3,numTotal)`.

Convention for solver parameters: `TurbModel%Scalar(i)` with `i` as listed per model below.

### 1.3 Momentum coupling for non-Boussinesq stresses

For models with a stress tensor that is not `2 mut S` the momentum equations get an explicit source. Convention, with `aij(6,numCells)` ordered
xx, xy, xz, yy, yz, zz:

```
T_ij = -rho <u_i u_j> = 2 mu_t S_ij - 2/3 rho k delta_ij - rho k a_ij
```

`asm_stress_terms` integrates `-rho k a_ij` over the faces (inner: linear interpolation with `facint`; inlet/outlet/pressure faces: owner value)
and is called from `calcuvw_assemble` when `lasm` is true. The matching eddy-viscosity part stays implicit in `vis` (stabilising).
Master had the module but with a wrong sign in `calcstress`, `aij` was not initialised and `lasm` not set.

- EASM (algebraic): `aij` is the closure's anisotropy in excess of Boussinesq, `mu_t = rho C*mu k tau`.
- RSM (transported): `mu_a = rho C_a k T` (C_a = 0.1, bounded by 1e3 mu) in `vis`, `a^ex_ij = R_ij/k - 2/3 delta_ij + 2 mu_a S_ij/(rho k)`.
  The solution does not depend on `mu_a` at convergence; it only sets how much of the stress is implicit. A bound on `mu_a/mu` was necessary:
  epsilon overshooting to ~0 gives `T -> infinity`, `mu_a` of 1e5 and divergence in the next momentum solve.

## 2. Models

### 2.1 `k_epsilon_zeta_f` (slots 1..4 = k, eps, zeta, f)

Hanjalic, Popovac, Hadziabdic (2004): `nut = Cmu zeta k T`, `Cmu = 0.22`; `Dzeta/Dt = f - zeta P/k`; `L^2 lapl f - f = (C1 - 1 + C2' P/eps)(zeta - 2/3)/T`;
`T = max(min(k/eps, 0.6/(sqrt6 Cmu S zeta)), 6 sqrt(nu/eps))`, `L = 0.36 max(min(k^1.5/eps, k^.5/(sqrt6 Cmu S zeta)), 85 (nu^3/eps)^.25)`,
`Ceps1 = 1.4(1 + 0.012/zeta)`, `Ceps2 = 1.9`. Walls: `k_w = 0`, `zeta_w = 0` (Dirichlet), wall-cell `eps = 2 nu k/y^2`, `f = -2 nu zeta/y^2`
(taken implicitly as a sink in zeta - explicit made the iteration oscillate); log-law branch for `y* > 11.63`.
Added later: production limiter `P <= 10 rho eps`, `zeta <= 2`. The airfoil startup still diverges (see README).
Defects of the previous module (patch 0033 replaced it, 1392 -> ~450 lines): `zeta` allocated but never initialised, wall cell of zeta overwritten with
`2 nu k/y^2`, factor `eps` twice in the epsilon sink, division by zero in the time scale; the result was NaN.

### 2.2 `transition_sst` (slots 1..4 = k, omega, gamma, Re_theta_t~)

Langtry-Menter (2009) on SST-2003: `Pk~ = gamma_eff Pk`, `Dk~ = min(max(gamma_eff, 0.1), 1) Dk`, `F1 = max(F1, F3)`, `F3 = exp(-(Ry/120)^8)`, intermittency
equation with `Fonset`, `Flength`, `Fturb`, separation-induced `gamma_sep`; `Re_theta_t` transport with the free-stream correlation (local fixed-point iteration for
the pressure-gradient parameter). Mean-flow side is SST (`modify_mu_eff` of `k_omega_SST.f90`, now public). Needs y+ ~ 1. Correlations are typed from the paper.

### 2.3 EASM (`easm_closure.f90`, `easm_rumsey_gatski_komega.f90`, `earsm_gatski_speziale.f90`)

Weak-equilibrium solution of the algebraic stress equation with SSG pressure-strain, three-term tensor basis (Gatski-Speziale 1993, Rumsey-Gatski form):

```
T_ij = 2 mut { S + a2 a4 (S W - W S) - 2 a3 a4 (S S - 1/3 {S S} I) } - 2/3 rho k I,   mut = rho C*mu k tau
a1 = (4/3 - C2)/2, a2 = (2 - C4)/2, a3 = (2 - C3)/2, a4 = tau/(g1 + 2 g0 C*mu eta^2 tau^2)
g0 = C1(1)/2, g1 = C1(0)/2 + (Ceps2 - Ceps1)/(Ceps1 - 1)
```
Constants: `C1(0) = 3.4, C1(1) = 1.8, C2 = 0.36, C3 = 1.25, C4 = 0.4, Ceps1 = 1.44, Ceps2 = 1.83`. `C*mu` is the root with lowest real part of a cubic (Cardano, trigonometric
branch; triple-root case handled - it gave NaN before; verified in Python: `C*mu(S tau = 3.3) = 0.0912`), floored at 5e-4. Excess anisotropy is scaled so that
`<u_i u_i>/k >= 0.04`.

`easm_rumsey_gatski_komega` (slots 1, 2 = k, omega): `omega = eps/k`; `Dk/Dt = P - f_b* rho omega k + div[(mu + mut) grad k]`;
`Domega/Dt = gamma (omega/k) P - beta rho omega^2 + div[(mu + sigma_w mut) grad omega]`; `Cmu = 0.0895, sigma_k = 1, gamma = 0.53, beta = 0.83,
sigma_w = sqrt(Cmu)(beta - gamma)/kappa^2`. `P = T_ij dU_i/dx_j` with the full algebraic stress, limited to `20 f_b* rho omega k`. Wall: `k_w = 0` Dirichlet, wall-cell
`omega = 10 * 6 nu/(beta y^2)`, y+ ~ 1.
**Bug found while verifying** (Cf 40 % low): the cross-diffusion switch of Wilcox, `chi_k = (1/omega_w^3) grad k . grad omega_w`, carries a factor `Cmu^2` when `omega` is `eps/k`
(`omega = Cmu omega_w`). Without it `f_b*` was ~1.7 over the whole outer half of the channel and k was destroyed too fast. Fixed: Cf 0.0038 -> 0.0064.
Variants tested before: wall omega factor 1 or 10, k wall Dirichlet or zero-flux, closure on/off (EVM): no effect - all gave the same low Cf, which pointed to the transport equations.

`earsm_gatski_speziale` (slots 1, 2 = k, eps): standard k-eps with `P` from the algebraic stress, time scale `max(k/eps, 6 sqrt(nu/eps))`, `sigma_k = 1, sigma_e = 1.3`; wall functions with
`Cmu = 0.09` (not `C*mu`). High-Re only.

### 2.4 `rsm_full`, `ebrsm` (`rsm_common.f90`; slots 1..6 stresses xx, xy, xz, yy, yz, zz; 7 = eps; 8 = alpha)

```
DR_ij/Dt = P_ij + Phi_ij - eps_ij + diffusion,     P_ij = -R_ik dU_j/dx_k - R_jk dU_i/dx_k
Phi_ij   = (1-alpha^2) Phi^w_ij + alpha^2 Phi^h_ij,   eps_ij = (1-alpha^2)(R_ij/k) eps + alpha^2 (2/3) eps delta_ij        (alpha = 1 in rsm_full)
Phi^h    = -(C1 + C1* P/eps) eps a + C2 eps (a a - 1/3 a:a I) + (C3 - C3* |a|) k S + C4 k (a S + S a - 2/3 a:S I) + C5 k (a W + W a)
Phi^w    = -5 (eps/k) ( R_ik n_j n_k + R_jk n_i n_k - 1/2 R_kl n_k n_l (n_i n_j + delta_ij) ),   n = grad alpha/|grad alpha|
Deps/Dt  = (Ceps1 P - Ceps2 eps)/T + diffusion
```
`a = R/k - 2/3 I`, constants (SSG in `a`, Manceau 2015): `C1 = 1.7, C1* = 0.9, C2 = 1.05, C3 = 0.8, C3* = 0.65, C4 = 0.625, C5 = 0.2`; `Ceps1 = 1.44, Ceps2 = 1.83`.

All source terms are evaluated once per iteration from the previous stresses (`sx`), and a positive multiple `c_x = (eps/k)(6(1-alpha^2) + alpha^2(C1 + C1* P/eps))` of the old stress is added to both
sides, which makes the matrix diagonally dominant but leaves the converged solution unchanged. After the solve: diagonal stresses >= small, Schwarz inequality on the shear stresses, `k = R_kk/2`.
Diffusion: Daly-Harlow `Cs T R_kl` (`Cs = 0.22`, epsilon `0.18`) contracted with the direction of the gradient of the diffused quantity (`ggdh`), scalar form for the engine.
The isotropic form `nu + Cs (2/3) k T` gave Cf 25 % low in the channel (over-diffusion normal to the wall since `R_nn < 2/3 k`).
epsilon is floored by `nu_t/nu <= 1e4` after each solve.

`rsm_full`: `alpha = 1`, `T = k/eps`, wall cells hold the log-layer stresses
`R = k_w(1.214 tt + 0.238 nn + 0.548 ss + sqrt(Cmu)(tn + nt))`, `k_w = u_tau^2/sqrt(Cmu)` (u_tau from the previous wall shear stress), wall-cell `eps = Cmu^.75 k^1.5/(kappa y)`, wall shear from the standard wall function.

`ebrsm`: `alpha - L^2 lapl alpha = 1`, `alpha_w = 0`, `L = 0.161 max(min(k^1.5/eps, k^.5/(sqrt6 * 0.22 * 2/3 * S)), 80 (nu^3/eps)^.25)` bounded by the largest wall distance of the domain
(without the bounds alpha blew up to 1e23 where epsilon had overshot), `T = max(k/eps, 6 sqrt(nu/eps))`, `Ceps1 = 1.44(1 + 0.03 (1-alpha^3) P/eps)` (Manceau 2015 as remembered - not checked against the paper text),
`R_w = 0`, wall-cell `eps = 2 nu k/y^2`. The alpha equation uses the engine without convection and time term.

### 2.5 `sgs_one_equation_dynamic` (slot 1 = k_sgs)

`Dk/Dt = P - Ceps k^1.5/Delta + div[(nu + nut) grad k]`, `nut = Ck Delta sqrt(k)`, `Delta = V^(1/3)`. Local dynamic `Ck`, `Ceps` with a face-based test filter
`phi^_P = (phi_P + sum_nb phi_nb)/(1 + n_nb)`; test-to-grid width ratio per cell from the second moment of the kernel, `Delta^2 = Delta^2 + 12 sigma^2` (2.1 on hexahedra). `Ck = <L^d_ij M_ij>/<M_ij M_ij>`,
`L_ij = (u_i u_j)^ - u^_i u^_j`, `M_ij = -2[Delta^ sqrt(k_test) S^_ij - Delta (sqrt(k) S_ij)^]`, `k_test = L_kk/2`; `Ceps` from the viscous dissipation of the band between the filters over the model
dissipation difference. The averaging `<.>` is the same face filter, then clipping (`Ck` in [0, 0.2], `Ceps` in [0.3, 2.5]) and relaxation 0.5.
**Simplification:** the term `|grad u^|^2` of the `Ceps` numerator is approximated by the cell-wise gradient energy, not by the gradient of the filtered field. `Ceps` sits on the lower clip most
of the time (mean 0.4, the literature values are ~1). Not validated against DNS or LES data. No backscatter, no wall model.

## 3. Verification numbers

| test | result |
|---|---|
| registry refactor | `k_omega_sst` restart bitwise identical |
| channel 200x50 (Dean 0.00730) | SST 0.00679, zeta-f 0.00677, EASM-RG 0.00640, EB-RSM 0.00745; zeta-f `u+ = y+` in the sublayer |
| channel 100x10, y+ ~ 28 | k-eps 0.00662, EARSM 0.00649, RSM 0.00545 |
| T3A-like plate | see README |
| lid-driven cavity 3D, LES | stable 1500 steps |
| airfoil NACA 0012 Re 3e6 | SST 0.00907, EASM-RG 0.00912, EB-RSM 0.01022, transition 0.00507; y+~18 mesh: k-eps 0.00660, EARSM 0.00634, RSM 0.00565; zeta-f diverges |

Every model was run to convergence on the channel. On the airfoil `k_epsilon_zeta_f` diverges (README, limitations). `k_omega_sst` and the other models need first-order convection and urf 0.5 to start the airfoil.

## 4. Bugs and weak points found in master (not fixed unless said)

1. `k_epsilon_zeta_f` module: NaN, see 2.1 (replaced).
2. `calcstress.f90`: sign of `aij`; `asm_stress_terms` unused/incorrect; `lasm` uninitialised; `aij` not initialised (fixed here).
3. `grad()` writes boundary values: a `(3,numCells)` gradient array overflows (we use `(3,numTotal)`; the interface of `calcsc`-style routines still declares `numCells`).
4. `k_omega_sst` plate and airfoil: NaN with default numerics (not fixed; workaround in the generators).
5. Time level arrays of model-owned fields are not shifted by `time_shift.f90` (handled in `shift_time_levels`); restart file lacks them.
6. `ncorr`, `npcor` have no default (declared in `parameters.f90`, uninitialised): an input without them silently runs SIMPLE with no pressure solve (the airfoil input of the coupled generator is such a file). Not fixed.

## 5. Adding a model

1. Write `my_model.f90`: module with `modify_viscosity_my_model` and `modify_viscosity_inlet_my_model`; use `turb_transport` for each equation.
2. Add `case ('my_model')` in `turbulence_registry.f90` (pointers and flags; `lasm = .true.` if you fill `aij`).
3. Add `my_model.f90` to `src/cappuccino/files` after `turb_transport.f90` and before `turbulence_registry.f90`.
