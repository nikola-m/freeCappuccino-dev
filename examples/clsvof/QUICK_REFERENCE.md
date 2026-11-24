# CLSVOF quick reference

## Build and run
```sh
make                                              # in the repository root -> bin/cappuccino
python3 examples/clsvof/mkcase.py rising-bubble case --n 40      # mesh, 0/U, vtk templates, input.nml, cells.npz
cd case && ../../bin/cappuccino input.nml monitor restart
python3 ../clsvof_post.py hysing .                # numbers of the benchmark
```
Cases of `mkcase.py`: `advect-vortex`, `advect-rotate`, `advect-sphere` (prescribed flow, `calcU = calcP = F`), `static-drop`, `rising-bubble`, `falling-drop`.
`sh examples/clsvof/run_validation.sh quick|full` generates, runs and evaluates all of them.

## Required settings in `input.nml` for a flow case
```
 twoPhase = T   ltransient = T   bdf = T         ! transient, 1st order implicit time stepping
 SIMPLE = T     maxit = 5   urfU = 0.9,0.9,0.9   urfP = 0.5   ! outer iterations per step
 lbuoy = F                                        ! gravity is added by the module: gravx/gravy/gravz
 densit = <rho_1>                                 ! reference density: force is (rho - densit) g, p is the p - rho_ref g.x
 timestep = ..                                    ! capillary limit  dt < sqrt((rho1+rho2) h^3 / (4 pi sigma)),  CFL(flow) < 1
```
Boundaries: `wall`, `symmetry`, `empty` only (closed domain). Initial condition of `0/U` is the velocity (zero).

## Namelist parameters (defaults)
| parameter | default | meaning |
|-----------|---------|---------|
| `twoPhase` | F | activates the solver |
| `cl_rho1`, `cl_rho2` | 1000, 1 | density of phase 1 (`vof = 1`, `lsf > 0`), phase 2 |
| `cl_mu1`, `cl_mu2` | 1e-3, 1.8e-5 | dynamic viscosities |
| `cl_sigma` | 0 | surface tension coefficient |
| `cl_ic` | `'sphere'` | initial interface: `sphere`, `cylinder` (axis z), `plane` (phase 1 below `z = cl_center(3)`) |
| `cl_inside` | 1 | phase inside the sphere/cylinder: 1 drop, 2 bubble |
| `cl_center(3)`, `cl_radius` | 0, 0.1 | center and radius |
| `cl_eps` | 1.5 | half width of the smoothed Heaviside function [cells] (surface tension only) |
| `cl_band` | 4 | half width of the band where `lsf` is the distance to the PLIC planes [cells] |
| `cl_nReinit` | 0 | PDE reinitialization sweeps outside the band (Sussman et al. 1994) |
| `cl_cfl` | 0.25 | Courant number of an advection substep (the step is subcycled) |
| `cl_fmin` | 1e-8 | `vof < cl_fmin` -> 0, `> 1 - cl_fmin` -> 1 |
| `cl_ab2` | T | 2nd order extrapolation of the volume flux for the advection |
| `cl_harmonic` | F | harmonic mixing of viscosity |
| `cl_flow` | `'solve'` | `translate`, `rotate`, `vortex`: prescribed flow (advection tests), `calcU = calcP = F` |
| `cl_vel(3)`, `cl_omega`, `cl_vcenter(2)`, `cl_T` | | parameters of the prescribed flows |
| `cl_inflowF` | 0 | volume fraction of what enters through open faces (prescribed flows only) |

## Output
* `clsvof_history.dat`, one line per time step: `t, V1, (V1-V10)/V10, clipped volume, centroid of phase 1 (3), centroid of phase 2 (3), mean velocity of phase 1 (3), of phase 2 (3), max|U|, interface area, continuity error, no. interface cells`.
* `vtk/<step>/`: `vof`, `lsf`, `kappa`, `rho` (+ `U`, `p`, `mueff` only if `calcVis`).
* `clsvof.restart` (with the other restart files; `lread = T` continues, `NUMSTEP` counts total steps).

## Public routines (`use clsvof`)
`clsvof_init`, `clsvof_advance`, `clsvof_monitor`, `clsvof_write_restart`; fields `vof`, `lsf`, `kappa`, `hs`, `flvol` (volume flux through faces), `fbal`, `fst`.
`use twophase_pressure, only: calcp_twophase`.

## Typical problems
| symptom | remedy |
|---------|--------|
| `contErr` grows, pressure solver stalls | lower `tolRelP`, raise `maxiterP`, more outer iterations `maxit` |
| interface wrinkles near sharp features, `clip` large | refine, lower `timestep`, `cl_cfl = 0.2`; a bubble of 5 cells radius is under-resolved |
| blow-up with density ratio 1000 and unresolved gas wake (gas Re >~ 1e4) | not supported at present, see TECHNICAL_DOCUMENTATION, section 9 |
| stops with 'boundary type ... not supported' | only `wall`, `symmetry`, `empty` |
| `lread = T` fails | `clsvof.restart` must be in the working directory |
