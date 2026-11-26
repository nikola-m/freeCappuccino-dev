# Moving (dynamic) mesh - RBF mesh deformation with the space conservation law

Module `moving_mesh` of freeCappuccino: laminar incompressible flow on meshes whose boundaries move with a prescribed motion.
The interior mesh follows the moving boundaries by **radial basis function (RBF) interpolation** of the boundary displacements
(Rendall & Allen, *Int. J. Numer. Meth. Engng.* 74 (2008) 1519-1559). The finite volume discretization is the arbitrary
Lagrangian-Eulerian one: the volume swept by every face is calculated **exactly** and the new cell volume is the old one plus
the sum of the swept volumes (space conservation law, Demirdzic & Peric 1988), the mass fluxes are made relative to the grid
and the pressure(-correction) equation gets the volume change of the cells as a source term, so a uniform flow stays uniform
on a mesh which breathes.

| document | content |
|----------|---------|
| this `README.md` | what it does, how to run, verification results, limitations |
| [`TECHNICAL_DOCUMENTATION.md`](TECHNICAL_DOCUMENTATION.md) | equations, algorithm, RBF, swept volumes, geometry update, code structure and hooks |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | cheat sheet: parameters, files, symptoms and remedies |

## Files

| file | content |
|------|---------|
| `src/mesh/moving_mesh.f90` | module `moving_mesh` (listed in `src/cappuccino/files` after `geometry.f90`) |
| `src/cappuccino/main.f90`, `init.f90` | `mm_advance` every time step, `mm_init` + checks of unsupported combinations |
| `src/finiteVolume/fvEqnDiscretization/Pressure/calcp_simple.f90`, `calcp_piso.f90`, `continuityErrors.f90` | relative mass flux, volume change source |
| `.../Velocity/velocity.f90`, `.../ScalarEqns/temperature.f90`, `concentration.f90` | unsteady term with `vol_unsteady` |
| `src/io/read_input.f90`, `writefiles.f90` | namelist parameters, VTK output with the current coordinates |
| `examples/movingMesh/mkcase.py` | generator of the cases (native mesh, `0/U`, VTK templates, `input.nml`, `case.npz`); reuses the writers of `examples/coupledPU/mkcase.py` |
| `examples/movingMesh/piston_check.py` | comparison of a piston run with the exact solution |
| `examples/movingMesh/run_validation.sh` | generates, runs and evaluates all cases below (about 5 s on one core) |

## Quick start

```sh
rm -f bin/cappuccino; make                 # repository root (see the note on make in QUICK_REFERENCE.md)
python3 examples/movingMesh/mkcase.py piston /tmp/piston --motion oscillate --amp 0.1 --freq 1.0 --steps 100
cd /tmp/piston && ~/freecappuccino-dev/bin/cappuccino input.nml monitor restart
python3 ~/freecappuccino-dev/examples/movingMesh/piston_check.py .
```
Output: the usual files, `vtk/<step>/` written with the **current** coordinates of the nodes (open `vtk/solution_fields_<step>.vtm` in ParaView), and in
`monitor` one line per time step: `Moving mesh: max displacement step = ..  min vol = ..  total vol = ..  GCL error = ..`.
All parameters: [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md). Briefly: `DynamicMesh = T`, `mm_patch = 'name'` (patches of type `wall`),
`mm_motion = 'translate' | 'oscillate' | 'rotate' | 'user'` and its parameters, `ltransient = T`, `bdf = T`, `SIMPLE = T` or `PISO = T`, `gauss = T`.

## Cases (`mkcase.py`)

* `piston` - channel `L x H` (2D one cell thick with `empty` faces, or a 3D box with `--ndim 3`), moving wall (piston) on the left, free slip walls (`symmetry`)
  on the other sides, `pressure` outlet (`p = 0`) on the right. Piston motion `translate` (`--U`) or `oscillate` (`--amp --freq`).
  The fluid is expelled through the outlet, the cells are squeezed by the mesh motion. Exact solution in `piston_check.py`.
* `body` - square (side 1) in a closed box (side 4, no-slip), the square is a moving wall: `oscillate` in y (amplitude 0.25, 0.5 Hz, the case of the Python RBF
  script this module is based on) or `rotate` (`--omega`). The box is closed and the body is rigid, so the fluid volume is constant: the checks are GCL and mass conservation.

## Verification (`run_validation.sh`, serial, gfortran `-O3`, numbers copied from its output)

Exact solution of the piston: with free slip walls the fluid moves as a plug. In every time step the velocity is the *mean velocity of the piston in the step*,
`u^n = (d(t_n) - d(t_n-1)) / dt` (continuity on the moving mesh), and the pressure gradient is the one of the acceleration of backward Euler,
`p^n(x) = rho (u^n - u^(n-1))/dt (L - x)`, `x` = laboratory coordinate of the cell centre (from the VTK points).

| case | max abs(u - u_exact) | max abs(v), abs(w) | p error | GCL error (max over run) | smallest cell volume |
|------|----------------------|--------------------|---------|--------------------------|----------------------|
| piston, constant velocity 0.5, 2D 40x10, SIMPLE, 40 steps | 4.4e-16 | 0 | 0 (p = 0) | 9.6e-15 | 4.0e-4 |
| the same, PISO | 1.2e-15 | 1.8e-17 | 0 | 9.6e-15 | 4.0e-4 |
| the same, 3D 20x6x6, 20 steps | 1.7e-15 | 3.6e-16 | 3.7e-14 (abs) | 5.6e-15 | 2.5e-3 |
| oscillating, x = 0.1 sin(2 pi t), 40x10, 100 steps (one period) | 1.3e-3 (bulk) | 3.6e-3 | 3.3e-3 (6.8e-3 rel.) | 9.3e-15 | 4.5e-4 |
| the same, 80x20 | 5.3e-4 (bulk) | 1.2e-3 | 2.2e-3 (4.4e-3 rel.) | 1.8e-14 | 1.1e-4 |

* Constant velocity: the plug flow is reproduced to the round-off, which tests the relative fluxes, the volume change source and the unsteady term together: all cells change their volume in every step.
* The positions of the piston in the VTK files agree with `d(t)` to 1e-16 (translate) and 1e-25 (oscillate, t = 1).
* **Accelerating piston**: the deviation is confined to the cells at the piston and halves with the mesh refinement (40 -> 80 cells: 0.030 -> 0.013 at the first cell). It comes from the
  zero gradient pressure condition of the code at the wall, the true one is `dp/dn = -rho a`. Therefore `piston_check.py` compares cells further than `0.25 L` from the piston
  (printed as bulk) and the deviation there is 2e-3 of the velocity (halves with refinement too).
* **Negative control**: with the `mm_pressure_source` call removed (binary built for this) the constant velocity case gives `max abs(u - u_exact) = 0.50`, so the test is sensitive to the volume change source.
* **Restart** (`lread = T`): 50 steps + 100 steps restarted reproduce the continuous run of 150 steps to 1.2e-12 in `u`, 4e-14 in `p` and 1e-14 in `flmass`; the mesh is moved to the time of the restart file.
* **Square in a closed box** (RBF, 8 cells per unit length, 100 steps, oscillation of 0.25): GCL error 4.6e-15, global continuity error of the absolute mass flux 6.5e-16 (max over the run),
  smallest volume 1.2e-3 (initial 1.6e-3). **Rotation** by 1 rad: GCL error 4.9e-15, continuity 1.3e-15, smallest volume 7.5e-4.

GCL error = `max_cells abs(vgeo - vol) / vol` with `vgeo` the cell volume calculated from the new nodal coordinates (`vol` is the sum of the swept volumes). It is at the round-off level, which is the point of the exact swept volumes.

## Limitations

* Laminar incompressible flow. `mm_init`/`init.f90` stop with a message for: steady runs, `bdf2`/`bdf3`/`cn`, turbulence models, compressible flow, periodic boundaries,
  the block-coupled p-U solver, two-phase flow, energy equation, MHD, least-squares and node based gradients (`gauss = T` only). SIMPLE and PISO only. Time scheme is first order backward Euler.
* The moving patches have to be of type `wall`; all other boundaries should be at rest (a warning is printed if a fixed boundary face sweeps volume).
* Prescribed motion only (no fluid-structure coupling); arbitrary motion by editing `mm_user_displacement`.
* **RBF deformation is dense**: matrix `N_ctrl x N_ctrl` (control points = nodes of all patches except `empty` and `symmetry`), LU factorization once, `O(N_ctrl^2)` memory (printed at start)
  and `O(N_ctrl * numNodes)` per time step. Fine up to roughly 10^4 control points; for larger boundaries a data reduction (greedy control point selection) is needed.
* RBF deformation fails for large rotations: the rotating square inverts cells at about 90 degrees (negative volume -> the run stops with a message). Reduce the motion, increase the support radius `mm_radius`, or use a larger gap.
* Nodal weights of the node interpolation used for output (`node_interpolation`) are not recomputed, so the point data in VTK files are slightly off on strongly deformed meshes; cell data are exact.
* Pressure at the accelerating moving wall: zero gradient condition (see above).

## Issues found in the existing code (not changed)

* `continuityErrors` does not include `pressure` patches, so the "global" error of a case with a `pressure` outlet is the outflow (it is 5.0e-2 = U*A for the constant velocity piston, not an error).
* `make` considers `bin/cappuccino` up to date after a source file is restored with `cp` (identical time stamps within the same second): `rm -f bin/cappuccino` before `make` when in doubt; the build is a single `gfortran` command.
* Lines longer than 132 characters are truncated by gfortran (errors "Line truncated").
* `calcp_simple_lis.f90` (LIS variant of the pressure solver) is not hooked up to the moving mesh.

## Reproduce

```sh
sh examples/movingMesh/run_validation.sh          # directory ./run is created
```
