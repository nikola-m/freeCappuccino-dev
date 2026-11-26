# Moving mesh quick reference

## Build and run
```sh
rm -f bin/cappuccino; make                         # repository root -> bin/cappuccino (make does not notice some source changes)
python3 examples/movingMesh/mkcase.py piston case --motion oscillate --amp 0.1 --freq 1.0 --steps 100   # mesh, 0/U, vtk templates, input.nml, case.npz
cd case && ../../bin/cappuccino input.nml monitor restart
python3 ../piston_check.py .                       # exact solution of the piston
```
Cases of `mkcase.py`: `piston` (`--ndim 2|3 --nx --ny --L --H --motion translate|oscillate --U --amp --freq --steps --dt --solver simple|piso --maxit --rbf`),
`body` (`--n --motion oscillate|rotate --amp --freq --omega --steps --dt`). `sh examples/movingMesh/run_validation.sh` runs all of them.

## Required settings in `input.nml`
```
 DynamicMesh = T                                  ! activates the module
 mm_patch = 'piston'                              ! up to 8 patch names, type wall in polyMesh/boundary
 mm_motion = 'translate'                          ! translate | oscillate | rotate | user
 ltransient = T   bdf = T                         ! transient, first order implicit time stepping (required)
 SIMPLE = T (maxit = 5..10)   or   PISO = T       ! pressure-velocity coupling
 gauss = T                                        ! gradients (least squares not supported)
 TurbModel = 'none'                               ! laminar
```
Other boundaries stay at rest. Moving patches: `Dirichlet` entry in `0/U` is replaced by the velocity of the moving face.

## Namelist parameters (defaults)
| parameter | default | meaning |
|-----------|---------|---------|
| `DynamicMesh` | F | activates the moving mesh |
| `mm_patch(8)` | `''` | names of the moving patches (type `wall`) |
| `mm_motion` | `'oscillate'` | `oscillate`: `d = mm_amp sin(2 pi mm_freq t)`; `translate`: `d = mm_vel t`; `rotate`: angle `mm_omega t` about the axis; `user`: routine `mm_user_displacement` |
| `mm_amp(3)`, `mm_freq` | 0, 1 | amplitude vector [m], frequency [1/s] |
| `mm_vel(3)` | 0 | translation velocity [m/s] |
| `mm_center(3)`, `mm_axis(3)`, `mm_omega` | 0, (0,0,1), 0 | point and direction of the axis of rotation, angular velocity [rad/s] |
| `mm_rbf` | `'wendland'` | `wendland` (C2, compact support, positive definite), `gaussian`, `multiquadric`, `inverse_multiquadric`, `thin_plate` |
| `mm_radius` | 0 | support radius of `wendland` [m]; `<= 0`: diagonal of the bounding box of the mesh |
| `mm_epsilon` | 1 | shape parameter of `gaussian`, `multiquadric`, `inverse_multiquadric` [1/m] |

## Output
* `monitor`, every time step: `Moving mesh: max displacement step, min vol, total vol, GCL error` (GCL error must be at the round-off level, ~1e-14).
* At start: number of control points, RBF type and radius, memory of the RBF matrix, RBF interpolation error at the control points (test displacement 1).
* `vtk/<step>/` written with the current node coordinates. Restart: `lread = T` moves the mesh to the time of the restart file (`NUMSTEP` = number of further steps).

## Public data and routines (`use moving_mesh`)
`DynamicMesh`, `xo,yo,zo` (old nodes), `volo` (old volumes), `flux_grid` (volume swept by a face, along `S`), `flmass_grid` (`den/dt*flux_grid`);
`mm_init`, `mm_advance`, `mm_pressure_source(su)`, `vol_unsteady(inp)` (= `volo` on a moving mesh, `vol` otherwise: use it in the first order unsteady term), `mm_write_points`, `mm_unsupported(what)`.

## Typical problems
| symptom | remedy |
|---------|--------|
| stops: `Moving mesh (DynamicMesh = T) cannot be used with: ...` | unsupported combination: use laminar, `bdf = T`, SIMPLE or PISO, `gauss = T`, no periodic |
| stops: `negative cell volume` | motion per step or amplitude too large for the RBF; reduce `timestep`/motion, raise `mm_radius`, refine near the moving wall |
| stops: `RBF matrix is singular` | coincident control points, or non-positive-definite RBF: use `wendland` |
| warning: `RBF matrix is ill-conditioned` | reduce `mm_radius` or use fewer boundary nodes |
| warning: `boundary faces that are not on moving patches sweep volume` | a fixed patch moves: list it in `mm_patch` or check the motion |
| `moving patch ... is not of type wall` / `does not exist` | patch type in `polyMesh/boundary`, spelling of `mm_patch` |
| `continuity errors global` equals the outflow | `pressure` patches are not included in `continuityErrors` (existing code), not a moving mesh error |
| velocity error at an accelerating wall | zero gradient wall pressure (`dp/dn = -rho a` not imposed): first order in the mesh size |
| large memory / slow start | dense RBF matrix `N_ctrl^2`; coarsen the boundary mesh or reduce control points |
