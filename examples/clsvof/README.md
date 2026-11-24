# CLSVOF - coupled level set and volume-of-fluid method for incompressible two-phase flow

Module `clsvof` (+ `twophase_pressure`) of freeCappuccino: 3D (and 2D one-cell-thick) incompressible two-phase flow with
a sharp, volume conserving interface on unstructured meshes of arbitrary polyhedra with planar faces.
The interface is tracked by the **volume fraction** `vof` (unsplit geometric PLIC advection: conserves the volume of each phase) and the **level set**
`lsf` (advected, then reset to the distance to the PLIC planes in a band: smooth normals and curvature), following
Sussman & Puckett, *J. Comput. Phys.* 162 (2000) 301-337 and Sussman, Smereka & Osher, *J. Comput. Phys.* 114 (1994) 146-159.
Surface tension is the CSF force with the balanced-force discretization of Francois et al. (2006), the pressure equation is for the **volume flux**
(`div u = 0`), not for the mass flux as in `calcp_simple`.

| document | content |
|----------|---------|
| this `README.md` | what it does, how to run, verification results, limitations |
| [`TECHNICAL_DOCUMENTATION.md`](TECHNICAL_DOCUMENTATION.md) | equations, algorithm, all discretizations, differences with the papers, code structure and hooks |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | cheat sheet: parameters, files, symptoms and remedies |

## Files

| file | content |
|------|---------|
| `src/finiteVolume/fvEqnDiscretization/TwoPhase/clsvof.f90` | module `clsvof` |
| `src/finiteVolume/fvEqnDiscretization/TwoPhase/twophase_pressure.f90` | module `twophase_pressure`, `calcp_twophase` |
| `examples/clsvof/mkcase.py` | generator of the cases below (native mesh, `0/U`, VTK templates, `input.nml`, `cells.npz`); reuses the writers of `examples/coupledPU/mkcase.py` |
| `examples/clsvof/clsvof_post.py` | evaluation: shape error, Laplace jump, Hysing benchmark quantities, free-fall acceleration |
| `examples/clsvof/run_validation.sh` | generates, runs and evaluates all cases, `quick` (about 10 min on one core) or `full` (adds the 128^2 advection and h = 1/80 bubble grids, about 30 min) |

## Quick start

```sh
make                                                           # repository root
python3 examples/clsvof/mkcase.py rising-bubble /tmp/bubble --n 40
cd /tmp/bubble && ~/freecappuccino-dev/bin/cappuccino input.nml monitor restart      # ~1.5 min
python3 ~/freecappuccino-dev/examples/clsvof/clsvof_post.py hysing .
```
Output: `clsvof_history.dat`, `vtk/<step>/` with `vof`, `lsf`, `kappa`, `rho`, `U`, `p` (open `vtk/solution_fields_<step>.vtm` in ParaView), `clsvof.restart`.
Parameters and everything else needed for a new case: [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md). Briefly: `twoPhase = T`, `ltransient = T`, `bdf = T`, `SIMPLE = T`, `lbuoy = F`,
`densit = rho_1`, `gravx/y/z`, `cl_rho1/2`, `cl_mu1/2`, `cl_sigma`, initial interface `cl_ic`, `cl_center`, `cl_radius`, `cl_inside`; boundaries `wall`, `symmetry`, `empty`.

## Verification (`run_validation.sh`, serial, gfortran 13.3 `-O3`, results below are copied from its output)

### 1. Advection of a disk (prescribed divergence free flow, no flow solver)

Volume fraction advected with the analytic volume flux (circulation of a vector potential along the face edges: exactly divergence free for every cell).
`E` is the L1 error `sum |f - f_exact| V` relative to the volume of the disk, after the flow returned the disk to its starting position.
Disk R = 0.15 at (0.5, 0.75), unit square, `n x n` cells.

*Solid body rotation, one revolution* (`omega = 1`, 10 n steps, Courant number of the step about 0.5):

| n | E / V_disk | order | volume change |
|---|-----------|-------|---------------|
| 32 | 2.39e-2 | | 4.2e-3 (initialization error of the disk, not advection) |
| 64 | 6.40e-3 | 1.90 | 4.7e-6 |
| 128 | 1.81e-3 | 1.82 | 3.0e-4 (same remark) |

The volume is conserved to 1e-9 (clipped volume `cl_clipVol` is 1e-12 in these runs).

*Reversing vortex* (LeVeque, T = 8, 8 n steps, `cl_cfl = 0.15`):

| n | E / V_disk | final volume / exact - 1 | clipped volume `cl_clipVol` / V_disk |
|---|-----------|--------------------------|-------------------------------------|
| 32 | 0.753 | +6.1e-3 | -1.9e-3 |
| 64 | 0.318 | +6.6e-3 | -6.6e-3 |
| 128 | 0.175 | +1.7e-2 | -1.7e-2 |

The error decreases slowly (factors 2.4 and 1.8) and the volume gained by clipping grows with refinement: the reversing vortex shows the weakness of the present unsplit scheme
for under-resolved filaments.

This is a severe test: the disk is stretched in a thin filament which becomes thinner than a cell. The volume error comes from the filament
(cells that get a negative or larger than 1 volume fraction from the overlapping swept regions,
clipped with redistribution to the neighbours, whatever cannot be redistributed is counted in `cl_clipVol`)
and depends on `cl_cfl`: for n = 64, `cl_cfl = 0.5`: volume + 3.2 %, `cl_cfl = 0.15`: +0.7 %.

### 2. Static drop: Laplace pressure jump and parasitic currents

Cylinder (2D) or sphere (3D) R = 0.25 in the unit box, no gravity, `rho = 1`, `mu = 0.05`, `sigma = 1`, equal densities, started at rest, velocity from curvature errors only
(200 steps of `dt = 0.5 sqrt(2 h^3/(4 pi sigma))`, `maxit = 4`). Exact jump `sigma/R = 4` (2D), `2 sigma/R = 8` (3D).

| | n | p jump | error | max parasitic velocity |
|---|---|--------|-------|-------|
| 2D | 16 | 4.1451 | +3.63 % | 2.4e-2 |
| 2D | 32 | 4.0166 | +0.42 % | 9.6e-3 |
| 2D | 64 | 3.9950 | -0.12 % | 3.5e-3 |
| 3D | 24 | 8.0176 | +0.22 % | 2.9e-2 |

Parasitic velocity scaled with `sigma/mu`: capillary number `mu U/sigma` = 4.8e-4 (n = 32) and 1.7e-4 (n = 64); the volume of the drop is constant to 1e-14 (n = 64).

### 3. Rising bubble: Hysing et al. (2009), test case 1

S. Hysing, S. Turek, D. Kuzmin, N. Parolini, E. Burman, S. Ganesan, L. Tobiska, *Quantitative benchmark computations of two-dimensional bubble dynamics*,
Int. J. Numer. Meth. Fluids 60 (2009) 1259-1288. Domain 1 x 2, bubble R = 0.25 at (0.5, 0.5), `rho_1 = 1000`, `rho_2 = 100`, `mu_1 = 10`, `mu_2 = 1`, `g = 0.98`, `sigma = 24.5`
(Re = 35, Eo = 10), no-slip top and bottom, free slip left and right (`symmetry`), up to t = 3, `dt = h/20` (1.25e-3 for h = 1/40), `maxit = 5`.
Reference values of the finest grids of the three groups of the paper:

| h | cells | y_c(3) | V_c,max | t(V_c,max) | min circularity | t(min circ.) | bubble area change | run time |
|---|-------|--------|---------|------------|------------------|--------------|--------------------|----------|
| 1/20 | 800 | 1.0693 | 0.2446 | 0.9283 | 0.303 (bubble breaks up) | 2.890 | 1e-6 | - |
| 1/40 | 3200 | 1.0837 | 0.2425 | 0.9191 | 0.8975 | 1.975 | 1e-8 | 80 s |
| 1/80 | 12800 | 1.0816 | 0.2419 | 0.9164 | 0.9003 | 1.950 | 4e-9 | 540 s |
| **reference** (Hysing et al.) | | **1.0813 / 1.0799 / 1.0817** | **0.2417 / 0.2421 / 0.2417** | 0.9213 / 0.9313 / 0.9239 | **0.9013 / 0.9011 / 0.9013** | 1.904 / 1.875 / 1.900 | | |

The centre of mass and the rise velocity are within 0.1 % of the reference at h = 1/80 (0.2 % and 0.3 % at h = 1/40), the minimum circularity within 0.1 % (0.4 % at h = 1/40).
A bubble resolved by 5 cells per radius (h = 1/20) is too coarse: it breaks up (this does not happen in the reference solutions on this grid); this is reported as it is.
Run times are wall clock on a 2-core machine with other jobs running at the same time.

### 4. Falling drop (density ratio 1000, closed box)

Drop (sphere in 3D, cylinder in 2D) `R = 0.2` in a box `1 x 2 (x 1)` falling from y = 1.5 under gravity `g = 1` in a gas: `rho_1 = 1000`, `rho_2 = 1`, `mu_1 = 0.1`, `mu_2 = 0.01`, `sigma = 1`
(gas Reynolds number about 30, Weber number of the gas about 1: laminar, nearly spherical drop). Everything is initially at rest.
Free fall of a heavy drop in a light fluid: the initial acceleration is `g (1 - r)/(1 + C_m r)`, `r = rho_gas/rho_drop`, `C_m = 1` (cylinder), `1/2` (sphere), that is the added mass of the displaced gas.
`clsvof_post.py fall` fits the first 10 % of the run:

| case | grid | measured | theory | difference | drop velocity at t = 1.2 | volume change |
|------|------|----------|--------|------------|-------------------------|---------------|
| 2D cylinder | 32 x 64 | 0.99626 | 0.99800 | -0.17 % | -1.1325 | 1e-8 |
| 3D sphere | 24 x 48 x 24 (27 648 cells) | 0.99696 | 0.99850 | -0.15 % | -1.1488 | 1e-8 |

(The 3D case: `python3 mkcase.py falling-drop dir` and `clsvof_post.py fall dir --added-mass 0.5`; the difference to the unbounded value is due to the walls and the viscous drag.)

**Limitation found**: the same drop with the properties of water and air (`mu_2 = 1.8e-5`, `sigma = 0.0728`, `g = 9.81`) is accurate in the first 0.3 s (initial acceleration 9.79 m/s2, theory 9.79) and then
blows up in every variant tried (3D and 2D; `sigma` 0.0728, 1 and 2; `dt` halved; `maxit` 10 and 12; `urfP` 0.3 and 0.8). The drop of 0.4 m diameter at 3 m/s has a gas Reynolds number of 2e4 with unresolved separation behind it, and the gas velocity
at the sides of the drop is 2-3 times the velocity of the drop. The cause was not found: it looks like an instability of the gas flow around the drop which the interface then amplifies (the interface length doubles before the blow-up), but this was not proven. The supported case is a drop of a size and velocity where the gas flow is laminar.
The harmonic interpolation of the Rhie-Chow coefficient was tried as a remedy and is worse (blows up at once).

## Limitations

* Closed domains only (`wall`, `symmetry`, `empty`); no contact angle model (90 degrees); no periodic, inlet, outlet or pressure boundaries in the flow solver (the advection tests allow open faces).
* First order implicit time integration (`bdf`), SIMPLE type iterations in each time step (`maxit` 4-5 with `urfU = 0.9`, `urfP = 0.5` were used), no non-orthogonal correctors: tested only on orthogonal hexahedral meshes (the geometry
  routines are for general polyhedra with planar faces, but only hexahedra were run).
* The time step is chosen by the user: no capillary or Courant based time step control. Interface advection is subcycled for `cl_cfl`.
* Volume conservation: exact up to `cl_clipVol` (printed), which is zero to rounding for rotation, translation, static drop, rising bubble, falling drop (1e-8 or less) but not for strongly stretched filaments (reversing vortex).
* Accuracy of the coupled method was measured, not proved: second order in space for the smooth disk rotation (observed order 1.8-1.9); the error of the static drop pressure jump drops by a factor 8.6 from n = 16 to 32 and 3.4 from 32 to 64.
* Parallel (MPI) not supported (the code of freeCappuccino is serial).

## Reproduce

```sh
sh examples/clsvof/run_validation.sh quick      # or: full
```
and apply with `git am 00*.patch` (see the commit messages). To remove: `git apply -R` the three patches in reverse order or `git revert`.
