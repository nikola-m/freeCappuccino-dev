# CLSVOF - technical documentation

Coupled level set and volume-of-fluid method for incompressible two-phase flow on the unstructured, face-based meshes of
freeCappuccino. Source files:

| file | module | content |
|------|--------|---------|
| `src/finiteVolume/fvEqnDiscretization/TwoPhase/clsvof.f90` | `clsvof` | interface: PLIC geometry, volume fraction and level set advection, VOF-LS coupling, properties, surface tension, gravity, balanced-force correction, monitoring, restart |
| `src/finiteVolume/fvEqnDiscretization/TwoPhase/twophase_pressure.f90` | `twophase_pressure` | pressure-velocity coupling with the volume flux (`calcp_twophase`) |
| hooks | `main.f90`, `velocity.f90` (`calcuvw_assemble`), `read_input.f90`, `init.f90`, `write_restart_files.f90`, `writefiles.f90` | see section 8 |

The method follows Sussman & Puckett, J. Comput. Phys. 162 (2000) 301-337 (CLSVOF), Sussman, Smereka & Osher, J. Comput. Phys.
114 (1994) 146-159 (level set), with the balanced-force surface tension of Francois et al., J. Comput. Phys. 213 (2006) 141-173.
**What is and is not as in the paper** is listed in section 9 - please read it before comparing with the paper.

## 1. Governing equations

Phase 1 (`vof = 1`, `lsf > 0`) and phase 2 (`vof = 0`, `lsf < 0`), interface `Gamma = {lsf = 0}`, normal `n = grad(lsf)/|grad(lsf)|`
pointing into phase 1:

```
rho (du/dt + u.grad u) = - grad p + div( mu (grad u + grad u^T) ) + sigma kappa grad H(lsf) + (rho - rho_ref) g
div u = 0
d vof/dt + u.grad vof = 0 ,    d lsf/dt + u.grad lsf = 0
rho = rho_2 + (rho_1 - rho_2) vof ,   mu = mu_2 + (mu_1 - mu_2) vof ,   kappa = - div n
```

`p` is the pressure without the hydrostatic part of the reference density `rho_ref = densit` (input parameter of the code), i.e. the
force `(rho - densit) g` is added to the momentum equations. `H` is the smoothed Heaviside function of half width
`cl_eps * h` (`h`: distance to the nearest cell center), `H = 1/2 (1 + phi/e + sin(pi phi/e)/pi)`.

The key difference to the single-phase pressure equation of freeCappuccino: the stock SIMPLE/PISO routines (`calcp_simple`,
`calcp_piso`) enforce `div(rho u) = 0`. With a density jump this is wrong, the constraint is `div u = 0`, it is the **volume flux**
`Q_f = u_f . S_f` which has to be conservative. `calcp_twophase` therefore solves for the volume flux, the mass flux of the momentum
equations `flmass` is a separate quantity (section 5).

## 2. One time step

`clsvof_advance` (called once per time step from `main`, after `time_shift`) :

1. **Volume flux** of the step: `qeff = (1 + r/2) Q^n - r/2 Q^(n-1)` (Adams-Bashforth, `r = dt/dt_old`; first step: `Q^n`; `cl_ab2 = F`: `Q^n`).
   For the prescribed flows of the advection tests (`cl_flow /= 'solve'`) the flux at the middle of the step is evaluated analytically.
2. **Face velocities** for the geometric flux: velocity at the face (interpolated cell velocities, or analytic) with its normal component
   replaced by `qeff/|S|` (so that `S.u_f = qeff` exactly).
3. **Substeps**: `nsub = ceil(max_P (sum_out Q dt / V_P) / cl_cfl)`.
4. **Volume fraction** (section 3), **level set** (section 4), in each substep.
5. **Coupling** (section 6): new PLIC planes, level set reset to the distance to the planes in a band.
6. **Properties**: `den`, `vis` from `vof`; `kappa`; force `fst`; balanced-force correction `fbal` (section 7); consistent mass flux `flmass` (section 5).

Then the usual outer iteration loop of `main`: `calcuvw` (momentum, first order implicit `bdf`) and `calcp_twophase` (pressure
correction); `maxit` iterations per time step.

## 3. Volume fraction: unsplit geometric (PLIC) advection

Each cell has the piecewise linear interface (PLIC) `n_P . x = d_P`. The plane normal comes from the level set (section 6), `d_P` is
found by bisection such that the volume of the cell on the positive side equals `vof_P V_P` (`plane_constant`); the volume is the exact
clipping of the cell, decomposed in tetrahedra (cell center, face center, two face nodes) by `tet_pos_volume`, exact for planar faces
and one of 4 cases (number of positive vertices). No division by differences of values of equal sign.

For each face the **swept region** is the face polygon displaced against the flow by `u_f dt` (slanted prism, "Lagrangian flux region", volume
`|S.u_f| dt = |Q| dt` for any `u_f`), decomposed in tetrahedra and clipped by the plane of the **donor** cell: `q1_f` is the volume of phase 1 which
crosses the face. Then

```
vof_P^new = vof_P - (1/V_P) sum_f q1_f + vof_P (dt/V_P) sum_f Q_f      (last term: zero for div-free flux, keeps vof = const invariant)
```

Pure cells (`vof < cl_fmin`, `> 1 - cl_fmin`) donate 0 or 1 without clipping. The regions of neighbouring faces overlap/leave gaps at
corners for non-uniform flow, so the result may leave `[0,1]`. `redistribute_defects` clips and passes the volume defect
to face neighbours (those with room, or those with phase 1, proportional to it, 6 Jacobi passes). **What is left is added to `cl_clipVol`**, which
is printed in `clsvof_history.dat`. In all tests below the volume change is 1e-8 or less of the volume.

*Why slanted regions:* the first version extruded the face along its normal. Then the regions through x- and y-faces of a cell
overlap at the corner (double count) and the volume of a translating disk grew by 3 % in 0.3 time units, 12 % after one revolution.
With `u_f`-slanted regions, translation and rotation conserve the volume to 1e-10.

Accuracy (section 10): second order in space for smooth interfaces in the rotation test.
Limitations: planar faces assumed (non-planar faces are triangulated with the face center: the clipped volume is exact for the triangulated cell but the cell volume `vtet` can differ from `vol`).

## 4. Level set advection

`lsf_t + u.grad lsf = 0` by the 3rd order TVD Runge-Kutta (Shu & Osher), the face value is taken from the upwind cell with the van Leer
limited linear-upwind reconstruction (`r = 2 grad(phi)_C . d_CD / (phi_D - phi_C) - 1`), `rhs = -(1/V) sum Q_f (phi_f - phi_P)`.
The advected `lsf` is **only used for the normal of the planes** (it is replaced right afterwards, section 6), so its accuracy
matters near the interface only. The gradient is the one selected in `input.nml` (`gauss`, `lstsq`, ...).

## 5. Mass flux consistent with the density (high density ratios)

Momentum equations of freeCappuccino are in the form `rho^(n+1)(u^(n+1) - u^n)/dt + sum_f F_f (u_f - u_P)/V = ...` (this is what the
assembly gives: the diagonal is the minus the sum of the off-diagonal convection coefficients). It is equivalent to the conservative form, if the mass
flux `F` and the density satisfy the discrete mass balance `(rho^(n+1) - rho^n) V + dt sum_f F_f = 0`. If the density comes from the level set and the
flux from the pressure solver it does not, the error is of the size of `(rho_1 - rho_2) u` in interface cells and the solution blows up for density ratios of 1000
(this is what we observed). Here the mass flux is built from the advection:

```
F_f = rho_2 Q_f + (rho_1 - rho_2) q1_f / dt                     (advection_mass_flux)
```

with `q1_f` the volume of phase 1 which crossed the face in the step (section 3), and `rho = rho_2 + (rho_1 - rho_2) vof`.
Then the balance holds up to the clipped volume. `flmass` is set once per time step and **not** changed by the pressure correction (the
volume flux is). Viscosity is a (arithmetic or, `cl_harmonic`, harmonic) mixture of `vof` too. The level set enters only through the
normal and the curvature/surface tension.

## 6. VOF-LS coupling (Sussman & Puckett, section 2.4 - 2.5)

After the advection of `vof` and `lsf` (`reconstruct_planes`, `lsf_from_planes`):

1. `n_P = grad(lsf)/|grad(lsf)|` in interface cells (`vof` in `(cl_fmin, 1-cl_fmin)`); if `|grad lsf| < 1e-6` the Youngs normal of `vof` is used.
2. `d_P` from `vof_P` (bisection, above).
3. `lsf` is **reset to the distance to the PLIC planes**: in interface cells the exact signed distance of the cell center to the plane of the cell; in
   cells up to `cl_band` face-neighbour layers from an interface cell the smallest distance to the plane fragments (disks of radius `h/2` in the planes) of
   the interface cells that reach them, signed by the side of the nearest plane; outside the band `+-(cl_band + 1.5) h`.
   This replaces the PDE reinitialization of the level set (Sussman et al. 1994), which would move the interface; here the zero level set is exactly
   the PLIC interface and the level set is exactly a distance function near it, in a volume conserving way.
4. Optional (`cl_nReinit > 0`): `cl_nReinit` sweeps of the Godunov-type reinitialization of Sussman, Smereka & Osher in the cells outside the band with
   fixed values in the interface cells and their neighbours (Russo & Smereka).

## 7. Surface tension, gravity, balanced force

Surface tension is the continuum surface force `f = sigma kappa grad H`, with `kappa = - div( grad(lsf)/|grad(lsf)| )` from the
distance function (Gauss, normals interpolated linearly to the faces and renormalized), limited to `|kappa| <= 2/h`, and `grad H` by Gauss
from the face-interpolated `H`. The cell force `fst` (surface tension + `(rho - densit) g`) is a source of the three momentum equations (`velocity.f90`).

In `calcp_twophase` the Rhie-Chow flux is (`K_f` = linear interpolation of `V/A_P`, `Df = S.S/(S.d)` the geometric diffusion coefficient):

```
Q*_f = U_f . S_f  -  K_f Df ( p_N - p_P - (grad p)_f . d_PN - fbal_f )
fbal_f = [ sigma kappa_f (H_N - H_P) + (rho_f - rho_ref) g . d_PN ]  -  (interpolated cell force fst).d_PN
```

The first bracket is the force in the compact (face) form, the second the interpolation of the cell forces which is contained in the
interpolated velocity `U_f`: the surface tension and gravity enter the face flux in the same form as the pressure difference and a static
state (a drop in equilibrium, a hydrostatic column) is an exact solution of the discrete equations except for the errors of the curvature
(`fbal` computed in `update_properties`).

## 8. Pressure-velocity coupling `calcp_twophase`

SIMPLE type: the momentum predictor (`calcuvw`) is called before, with `apu = 1/A_P`.

```
sum_f [ Q*_f + a_f (p'_N - p'_P) ] = 0 ,    a_f = - K_f Df ,   Q_f = Q*_f + a_f (p'_N - p'_P)
u' = - (V/A_P) grad p'  (Gauss, gradp_and_sources),   p = p + urfP (p' - p'_ref)
```

Boundary faces have zero volume flux (closed domain: `wall`, `symmetry`, `empty`), `p'` has zero normal gradient. The linear solver and its
tolerances are those of the stock pressure equation (`lSolverP`, `maxiterP`, `tolAbsP`, `tolRelP`). `cl_contErr = max_P |sum_f Q_f| dt / V_P` is printed in each outer iteration and
written in `clsvof_history.dat`; it is 1e-7 or less in the cases below. `npcor` is ignored (no non-orthogonal correctors: orthogonal or mildly
non-orthogonal meshes).

### Hooks into the existing code

| file | change |
|------|--------|
| `src/cappuccino/files`, `directories` | `clsvof.f90` before `velocity.f90`, `twophase_pressure.f90` after `pressure.f90`, directory `TwoPhase` |
| `main.f90` | `clsvof_advance` after `CourantNo`; `calcp_twophase` instead of `calcp_simple/piso` if `twoPhase`; `clsvof_monitor` at the end of the time step |
| `velocity.f90` | in `calcuvw_assemble` the source `fst*vol` (the same loop where `lbuoy` is) |
| `read_input.f90` | namelist parameters `twoPhase`, `cl_*` |
| `init.f90` | `clsvof_init` at the end (after the restart is read and the gradient matrices are built) |
| `write_restart_files.f90`, `writefiles.f90` | `clsvof.restart` (stream file: `cl_time, dt_old, V0, clip`, `vof`, `lsf`, `flvol`, `flvolo`); `vof`, `lsf`, `kappa`, `rho` in the .vtu |

`clsvof_init` sets `viscos = min(mu1, mu2)`: the wall boundary viscosity of the momentum equations is `max(viscos, visw)`, and `visw` is
set by `update_properties` to the viscosity of the wall cell. Use `lbuoy = F` (gravity is in `fst`, `deno` of the buoyancy branch is not used).

## 9. Differences with the paper, assumptions

* Sussman & Puckett use directional splitting on structured grids with the VOF flux from the PLIC reconstruction in 1D sweeps and a level set advected by
  a 2nd order ENO scheme. On unstructured meshes there is no splitting: here unsplit swept-region fluxes with redistribution. Order of accuracy was measured (section 10), not proved.
* They solve the momentum equations in non-conservative form with the density from the smoothed level set (projection method, 2nd order in time). Here: freeCappuccino
  `bdf` (1st order in time), SIMPLE iterations, and the density from `vof` with the mass flux of section 5. The smoothed Heaviside function is used for surface tension only.
* The level set is reset to the exact distance to the planes in a band (as in the paper, but with the plane fragments of the neighbouring interface cells rather than only the own plane);
  the PDE reinitialization is optional and off in all tests.
* Not implemented: contact angles (the level set has zero normal gradient at walls: 90 degrees), open boundaries in the flow solver (the code stops: `wall`, `symmetry`, `empty` only), periodic boundaries,
  axisymmetry, adaptive time step based on the capillary constraint (the user chooses `timestep`: `dt < sqrt((rho_1 + rho_2) h^3 / (4 pi sigma))`), parallel runs.
* A very light gas with unresolved wake: a water-like drop (density ratio 1000, gas Reynolds number about 2e4) falling at about 1 - 3 length scales per second goes unstable when the gas wake develops
  (see README, section "Falling drop"). The same case with a laminar gas (Re of order 30) is stable. This is a limitation of the present implementation (no turbulence model, 1st order
  time integration, SIMPLE), not of the cases shown.

## 10. Verification results

See `README.md` for the tables and how to reproduce them (`run_validation.sh`).

## 11. Issues found in the existing code (not changed)

* `calcp_simple`/`calcp_piso` use `div(rho u) = 0`, unsuitable for variable density incompressible flow (hence `calcp_twophase`).
* `time_shift.f90`: `deno` is used with `lbuoy` but allocated only for compressible flows.
* The boundary viscosity of the wall in `calcuvw_assemble` is `max(viscos, visw(iWall))` (constant `viscos`, not the local one): worked around in `clsvof_init`.
* In `calcp_piso` the source of the momentum corrector (`rU`) seems to contain the old pressure gradient (not used here, not investigated).
