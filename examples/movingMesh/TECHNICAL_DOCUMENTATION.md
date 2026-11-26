# Moving mesh - technical documentation

Implementation of `src/mesh/moving_mesh.f90` (module `moving_mesh`) and of its hooks in the flow solver. Notation: `P`, `N` owner and neighbour cell
of a face `f`, `S_f` the area vector (outward from the owner), `V` cell volume, `F_f` the mass flux through a face (`flmass`), `dt` the time step,
superscript `o` the previous time level. Reference papers:

* T.C.S. Rendall, C.B. Allen, *Unified fluid-structure interpolation and mesh motion using radial basis functions*, Int. J. Numer. Meth. Engng. 74 (2008) 1519-1559.
* I. Demirdzic, M. Peric, *Space conservation law in finite volume calculations of fluid flow*, Int. J. Numer. Meth. Fluids 8 (1988) 1037-1050.
* J.H. Ferziger, M. Peric, *Computational Methods for Fluid Dynamics*, 3rd ed., Springer 2002, Sec. 12.2.
* H. Wendland, *Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree*, Adv. Comput. Math. 4 (1995) 389-396.

## 1. Governing equations (ALE, finite volume)

For a cell with moving faces, the integral transport equation of a quantity `phi` with the flux relative to the moving surface reads

    d/dt ( int_V rho phi dV ) + sum_f F_f^rel phi_f = diffusion + sources,        F_f^rel = F_f - F_f^grid,    F_f^grid = rho_f flux_grid_f / dt

where `flux_grid_f` is the volume swept by face `f` in the time step (positive along `S_f`) and `F_f` is the absolute mass flux (through a face fixed in space).
The volume of a cell follows the **space conservation law** (GCL):

    V - V^o = sum_f flux_grid_f                                                   (1)

which must hold *identically* in the discrete form (not only approximately), otherwise a uniform flow is not preserved on a moving mesh: the
scheme would create or destroy mass in every cell proportionally to the error of (1). Continuity (absolute mass flux divergence free) in discrete form is

    sum_f F_f = 0   <=>   sum_f F_f^rel + sum_f F_f^grid = 0   <=>   sum_f F_f^rel = rho (V - V^o)/dt          (incompressible, uniform rho)

so the pressure correction equation, built from the relative fluxes, gets the source `rho (V^o - V)/dt` (section 6).

## 2. One time step

`mm_advance` is called from the time loop in `main.f90`, after the time level is shifted (`time_shift`), the time is advanced and before the equations are solved:

1. **Swap**: `xo = x`, `volo = vol` (also the old face centres of the moving patches). The old data are module data (`xo,yo,zo` `[numNodes]`, `volo` `[numCells]`).
2. **New nodal coordinates** by RBF deformation of the displacements of the control points (section 3), for the new time `t`.
3. **Swept volumes** `flux_grid` of all faces (section 4), new cell volumes `vol = volo + sum flux_grid` (1). A non-positive volume stops the run with a message.
4. **Mass fluxes of the grid** `flmass_grid = den_f/dt * flux_grid`, `den_f` linear interpolation with `facint` on inner faces (owner value on boundary faces).
5. **Geometry**: `xc,yc,zc,xf,yf,zf,arx,ary,arz,facint, Df`, the auxiliary points, wall and symmetry distances are recalculated from the new nodes (section 5).
6. **Relative fluxes**: the inner face mass fluxes are made relative, `flmass = flmass - flmass_grid`. At the beginning of the next step (`mm_step`, `shift_flux`) `flmass_grid` of the
   previous step is added back so that the absolute flux (of the old time level) enters the next relative flux. The flow solver recomputes `flmass` from the velocity in each outer iteration
   (`facefluxmass2`); there the grid flux is subtracted again (`calcp_simple`, `calcp_piso`).
7. **Pressure correction source** `su(inp) = su(inp) + den(inp)/dt (volo(inp) - vol(inp))` is added in `calcp_simple`, `calcp_piso` and `continuityErrors` by `mm_pressure_source`.

The monitor line reports `max |x - xo|`, the smallest and total volume and the GCL error `max abs(vgeo - vol)/vol`, where `vgeo` is the volume from the new nodal coordinates (divergence theorem with the face
triangulation of `geometry.f90`).

## 3. Mesh deformation with RBF

**Control points**: nodes of all patches except `empty` and `symmetry` (a node belonging to a moving and a fixed patch gets the moving prescription; their number is printed). The prescribed
displacement `d_j` of a control point is the motion of the patch (`oscillate`: `amp sin(2 pi f t)`, `translate`: `vel t`, `rotate`: Rodrigues rotation about `mm_axis` through `mm_center` by the angle `omega t`,
`user`: `mm_user_displacement`) on the moving patches and zero elsewhere. Interior nodes and the nodes of `symmetry`/`empty` patches are free.

**Interpolation** (Rendall & Allen): the displacement of any node `x` is

    d(x) = sum_j alpha_j phi( |x - x_j| ),         sum_j alpha_j phi( |x_i - x_j| ) = d_i   for all control points i

The matrix `M_ij = phi(|x_i - x_j|)` of the *initial* control point positions is assembled and LU factorized (dense, partial pivoting) once in `mm_init`; every time step only the right-hand side changes
(`lu_solve` for three components, then evaluation at all nodes). Displacements are measured from the **initial mesh** (`xini,yini,zini`), so the mesh is a function of the time only, and the restart is trivial: the mesh is moved to the time of the restart file.
The `mm_check_interpolation` routine reports how well a test displacement of 1 is reproduced at the control points.

**Radial functions** `phi(r)` (`mm_rbf`): `wendland` = `(1 - r/R)^4 (4 r/R + 1)` for `r < R`, else 0 (C2, compact support, positive definite; `R = mm_radius`, default is the diagonal of the
bounding box so the support is the whole mesh); `gaussian` `exp(-(eps r)^2)`; `multiquadric` `sqrt(1 + (eps r)^2)`; `inverse_multiquadric` `1/sqrt(1 + (eps r)^2)`; `thin_plate` `r^2 ln r`.
Only the positive definite ones (`wendland`, `gaussian`, `inverse_multiquadric`) guarantee a nonsingular matrix without a polynomial term, which is not included. The Python script that this follows (`_rbf_function`,
`_calculate_distance_matrix`, `update_mesh`) uses the same family of functions.

**Cost**: memory `8 N_ctrl^2` bytes (printed), factorization `O(N_ctrl^3)` once, per step `O(N_ctrl^2 + N_ctrl numNodes)`.

## 4. Swept volumes (exact GCL)

During a step the nodes move linearly, `x(t) = xo + t (x - xo)`, `t in [0,1]`. A face is split into the triangles `(1, i+1, i+2)` of `geometry.f90` (planar faces: same normal as the face, general faces: the same triangulation that defines `arx,ary,arz` and the volume).
The volume swept by a triangle is the time integral of the flux of the surface velocity through it. The velocity of the points of the triangle is the linear interpolation of the three nodal velocities; the triangle's
area vector `S(t)` is a **quadratic** function of time, so the flux

    dV = int_0^1 S(t) . dc dt,        dc = displacement of the centroid of the triangle

is integrated *exactly* by Simpson's rule: `dV = ( S(0) + 4 S(1/2) + S(1) ) . dc / 6` (`mm_swept_volumes`, with `S(1/2)` from the triangle of the mid-step node positions).
The volume of a polyhedron bounded by such triangles is a cubic function of time (and exactly the integral of the swept volumes), so (1) is satisfied to the round-off, whatever the deformation
(`GCL error` ~ 1e-14 in all cases of the README).

## 5. Recalculation of the geometry

`mm_update_geometry` replicates the calculation of the `geometry` module for the new nodes: cell centres and volumes (`vgeo` only for the GCL check; the actual `vol` is from (1)), face centres, area vectors, `facint`
(all three variants of the code), `Df`, the auxiliary points `xpp, xnp`, and the distances `dns, srds, dnw, srdw` to symmetry and wall faces. The routine is a copy of the formulas of `geometry.f90` (differences in the formulas would be
seen as the GCL error); when `geometry.f90` changes its formulas the routine has to follow.

The velocity of the fluid at the moving wall faces `(u,v,w)(ijb)` is the velocity of the face centre, corrected along `S_f` so that its normal component is `flux_grid/dt` exactly:

    U_b = v_f + S_f ( flux_grid/dt - v_f . S_f ) / |S_f|^2

so there is no flux through the moving wall relative to the grid and the no-slip condition holds (tangential part of the face motion).

## 6. Hooks in the flow solver

| place | change |
|-------|--------|
| `main.f90` | `if( DynamicMesh ) call mm_advance` inside `if( ltransient )` after `CourantNo` and the time update |
| `init.f90` | at the end: unsupported combinations stop (`mm_unsupported`), then `mm_init` (check of input, control points, RBF matrix, restart) |
| `calcp_simple.f90`, `calcp_piso.f90` | `flmass(i) = flmass(i) - flmass_grid(i)` after the Rhie-Chow flux of every inner face; `mm_pressure_source(su)` before the pressure equation is solved |
| `continuityErrors.f90` | `mm_pressure_source(su)`, so the reported error is that of the absolute flux |
| `velocity.f90`, `temperature.f90`, `concentration.f90` | `apotime = den*vol_unsteady/dt` in the first order unsteady term |
| `read_input.f90` | the namelist parameters `DynamicMesh`, `mm_*` |
| `writefiles.f90` | the `Points` of the VTK files are written from the current `x,y,z` |

**Unsteady term.** Convection in this code is assembled as `sum_f F_f (phi_f - phi_P)`, i.e. `sum F phi_f - (sum F) phi_P` with the continuity term subtracted. The conservative first order
scheme `[rho V phi - rho^o V^o phi^o]/dt + sum F^rel phi_f = 0` with (1) and discrete continuity `sum F^rel = rho (V - V^o)/dt` becomes
`rho V^o (phi - phi^o)/dt + sum_f F^rel (phi_f - phi_P) = 0`, hence `den * volo` (function `vol_unsteady`) multiplies the diagonal and the old value. It reduces to `den*vol` on a fixed mesh and preserves a uniform
field on a moving mesh *for any `dt` if the flux satisfies continuity* (tested by the plug flow of the piston).

**Pressure(-correction) equation.** Fluxes entering the pressure equation are the relative ones, `Rhie-Chow` flux of the absolute velocity minus the grid flux; the correction `a_f (p'_N - p'_P)` is the same.
The source `den (volo - vol)/dt` is added to `su` which, by the sign convention of the code, is `-sum_f F^rel` of the cell. At convergence `sum F^rel = rho (V - V^o)/dt`, i.e. the absolute flux is divergence free.
With uniform `den` the same equation is obtained for PISO (source added before the non-orthogonal correction loop).

## 7. Assumptions and limitations

* Incompressible laminar flow, bdf1 (the BDF2 and CN branches would need the volumes of several time levels: not implemented), SIMPLE and PISO, gauss gradients (the other gradient
  options were not adapted and not tested), no periodic patches, no energy equation (`calcEn`), `twoPhase`, MHD, coupled p-U, turbulence models (not tested; the wall distance `wallDistance` is not updated).
* Nothing else in the discretization needs the grid velocity except the wall velocity (section 5): all convective fluxes are the relative ones.
* Boundary conditions: the moving patches are `wall`. A moving `inlet`/`pressure`/`outlet` patch is not supported (the flux through it would need the relative velocity).
* Node interpolation data for the VTK point data (`node_interpolation`) is not updated, so the point data are approximate on strongly deformed meshes.
* Pressure at the moving wall is extrapolated with zero gradient (the code's default wall condition), which gives a first order error in the cells at an accelerating wall.
* RBF: the dense matrix limits `N_ctrl` to about 10^4; large rotations invert cells (RBF are good for translation and moderate rotation of a body inside a mesh).

## 8. Verification

See the README. The constant velocity piston (2D, 3D, SIMPLE, PISO) reproduces the plug flow to 1e-15 and the GCL error is 1e-14; the negative control without the volume source gives an error of 0.5.

## 9. Issues found in the existing code (not changed)

See the README: `continuityErrors` does not count `pressure` patches; `make` may not rebuild after a restored source file; 132 character lines; `calcp_simple_lis.f90` not hooked.
