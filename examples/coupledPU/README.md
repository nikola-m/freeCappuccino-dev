## Block-coupled pressure-velocity solver (`coupled = T`)

Instead of the segregated SIMPLE/PISO loop, momentum and continuity are assembled into **one** block-sparse system
with a 4x4 block (u, v, w, p) per cell, and solved simultaneously by **FGMRES preconditioned by an aggregation
algebraic multigrid** working on the 4x4 blocks. Source files:

| file | content |
|------|---------|
| `src/linearSolvers/block_amg.f90` | block CSR, 4x4 kernels, aggregation AMG (V/W cycle), block GS / block ILU(0) smoothers, dense coarsest-level LU, FGMRES(m) |
| `src/finiteVolume/fvEqnDiscretization/Pressure/coupled_pU.f90` | assembly, boundary conditions, solution driver `calc_coupled_pU` |
| `src/finiteVolume/fvEqnDiscretization/Velocity/velocity.f90` | `calcuvw` split into `calcuvw_assemble` (reused by the coupled solver) + the segregated solve |

### Discretization

* **Momentum**: convection, diffusion, unsteady, body-force and boundary terms come from the same code as the segregated
  solver (`calcuvw_assemble`). The pressure gradient is *implicit* (Gauss, linear interpolation `p_f = (1-λ)p_P + λp_N`),
  so every row has the 3x1 block `G` in its p column. Under-relaxation is implicit (diagonal divided by `urfU`).
* **Continuity**: `Σ_f ṁ_f = 0` with the Rhie-Chow flux
  `ṁ_f = ρ_f U_f·S_f − ρ_f (Vol/a)_f Df [ (p_N − p_P) − (∇p)_f·d_PN ]`.
  `U_f` (linear interpolation) and the compact pressure difference are implicit; the wide-stencil gradient term uses the
  previous iterate (deferred). It gives the 1x3 block `D` and the pressure Laplacian block `C`; for a linear pressure field the flux
  reduces to `ρ U_f·S`. `(Vol/a)` uses the *un-relaxed* diagonal as in SIMPLE, so the converged result does not depend on `urfU`.
  Mass fluxes after the solve are evaluated with the same formula, hence they satisfy continuity to the linear-solver tolerance
  and the next momentum assembly convects with a conservative flux.
* The saddle-point structure `[A G; D C]` with `D ≈ −ρ Gᵀ` is positive real, which is what the smoothers and AMG rely on.

### Boundary conditions (implicit where possible)

| type in `polyMesh/boundary` | momentum | continuity |
|-----------------------------|----------|-----------|
| `inlet`    | as in SIMPLE; `p_b = p_P` implicit | prescribed mass flux `flmass` |
| `wall`     | as in SIMPLE; `p_b = p_P + ∇p·d` (gradient lagged) | zero flux |
| `symmetry` | as in SIMPLE; `p_b = p_P` | zero flux |
| `outlet`   | zero-gradient velocity; `p_b = p_P` | extrapolated flux scaled to the inlet mass flow (as `adjustMassFlow`); **pressure level fixed by replacing the continuity row of cell `pRefCell` with `p = p_ref`** |
| `pressure` | fixed boundary pressure `p(ijb)` goes to the right-hand side | Rhie-Chow flux between cell centre and face, boundary velocity `U_b = U_P − S (Vol/a)/(S·d)[p_b − p_P − ∇p·d]` |
| `empty`    | handled as in the segregated solver | zero flux |

Not supported yet (the code stops with a message): periodic boundaries, Crank-Nicolson time stepping, compressible flow.
Memory is about 16 doubles per matrix non-zero (block values) on top of the scalar CSR.

### Input (namelist `input_parameters`)

```
 COUPLED = T,          ! use the block-coupled solver (SIMPLE and PISO are then ignored)
 urfU = 0.95, 0.95, 0.95   ! implicit velocity relaxation; much closer to 1 than for SIMPLE (0.9-0.99)
 maxit = 1             ! outer (Picard) iterations per time step; steady: 1, transient: 3-5
 cpl_tolRel = 0.05     ! relative reduction of the linear residual per outer iteration
 cpl_maxit = 40, cpl_restart = 20    ! FGMRES
 cpl_urfP = 1.0        ! pressure relaxation (1 = none)
 cpl_rebuild = 10      ! recompute AMG aggregates every N outer iterations (values are refreshed every iteration)
 amg_smoother = 'ilu0' ! or 'gs'
 amg_cycle = 2         ! 1 = V, 2 = W
 amg_theta = 0.25, amg_coarse_size = 16, amg_nu1 = 1, amg_nu2 = 1, amg_ccf = 1.0, amg_verbose = 0
```

`resor(1:3)` hold the L1 norms of the momentum residuals of the assembled system (this is what the convergence
check `tolerance` looks at), `resor(4)` the continuity residual. Each outer iteration writes a line
`Coupled p-U: FGMRES(AMG) iterations = ..` to the monitor file.

### Try it

```
cd examples/coupledPU
python3 mkcase.py cavity cav64 --nx 64 --Re 100        # also: channel (--outlet pressure|outlet), cavity3d
cd cav64
../../../bin/cappuccino input-simple.nml  monitor-simple  restart-simple
../../../bin/cappuccino input-coupled.nml monitor-coupled restart-coupled
python3 ../read_restart.py restart-coupled 64 64        # centreline u vs Ghia et al., Re=100
```

All case generators (`mkcase.py`, `mkcouette.py`, `mkairfoil.py`, `mknozzle.py`) write the ParaView mesh templates the VTK writer needs
(`vtk/mesh/mesh_0_0.vtu` and `mesh_<ib>_0.vtp` for every boundary region `ib` in the order of `polyMesh/boundary`), so the solver writes
`vtk/<step>/*.vtu`, `vtk/<step>/boundary/<patch>.vtp` and `vtk/boundary_<step>.vtm` for ParaView. `vtkmesh.py` builds them from `polyMesh`
(the cell connectivity that the native format lacks is reconstructed from `faces`, `owner`, `neighbour`; hexahedral cells only) and can also be
run on an existing case: `python3 vtkmesh.py <case>`.

### AMG stand-alone test

`test_block_amg.f90` builds a 2D saddle-point system with the same block structure and checks mesh-independent convergence:
```
gfortran -O2 ../../src/linearSolvers/block_amg.f90 test_block_amg.f90 -o t && ./t ilu0 && ./t gs
```

### Results obtained while developing (gfortran -O3, one core)

Outer iterations to a momentum L1 residual of 1e-9 (steady lid-driven cavity, Re 100; SIMPLE urfU/urfP = 0.8/0.3, coupled urfU = 0.98):

| mesh | SIMPLE its / s | coupled its / s |
|------|----------------|-----------------|
| 32²  | 390 / 0.5      | 142 / 0.4 |
| 64²  | 1272 / 6.8     | 200 / 1.5 |
| 128² | 4339 / 110     | 414 / 15  |
| 16³ (3D) | 133 / 1.0  | 47 / 0.8 |

FGMRES needs 1-3 iterations per outer iteration; for the synthetic test the AMG-FGMRES count for 1e-10 stays at 13-14 from 1k to 16k cells.
Fields from both solvers agree (max pointwise difference up to about 2e-2, located at the lid corners), also for a stretched
and sheared (non-orthogonal) cavity mesh, for channels with both outlet types (Poiseuille pressure gradient −0.2356 vs −0.24 analytic on 40x10 cells),
and for the transient 2D elbow of `examples/elbow2D` (rms velocity difference 1e-2 after 40 s).

## All-speeds (compressible) extension (`coupled = T`, `compressible = T`)

The same block-coupled solver handles compressible flow of an ideal gas, from low Mach number to supersonic
(Darwish & Moukalled type all-speeds pressure-based coupled algorithm). Switch on with `compressible = T` together with `coupled = T`.

### Discretization

* Unknowns of the coupled system are still (u,v,w,p); density is `rho = p/(R T)` with the temperature lagged in the p-U step
  (`Rair` from `parameters`, `c_p = 1006`). The continuity equation `d(rho)/dt + sum_f mdot_f = 0` has `rho = C p` (C = 1/(R T)).
* Face mass flux with density of the **upwind** cell: `mdot_f = rho_up Vn_f`, `Vn_f` the Rhie-Chow volume flux
  (`Uf.S - (Vol/ap)_f Df [(p_N-p_P) - (grad p)_f.d]`). The product is Newton-linearised around the current iterate:
  `mdot_f = rho_up* Vn_f + Vn_f* C_up (p_up - p_up*) + beta_f (rho_c - rho_up)* Vn_f*`.
  The middle term is the pressure-density coupling: a convection-like term of pressure in the continuity equation
  (it is what makes the system hyperbolic for M>1 and it also removes the incompressible-limit singularity of the pressure block).
  The last term is a deferred central-density correction, `beta_f = cpl_rhoBlend * max(0, 1 - M_f^2)`: second order for subsonic flow,
  pure upwind density in supersonic regions.
* The mass fluxes handed to the momentum and energy equations after each solve are the ones of the *linear* system that was solved,
  so continuity is satisfied exactly at every outer iteration.
* Unsteady term `(tc0 rho - tc1 rho_o + tc2 rho_oo) Vol/dt` (BDF1 for the first step, then BDF1/BDF2 as selected by `bdf`/`bdf2`), implicit in p.
* **Energy** (module `energy_equation`, `ScalarEqns/energy.f90`; `allspeeds.f90` is a thin layer on top of it): solved segregated after
  the p-U step. Formulation by namelist (priority in this order): `solveEnthalpy`, `solveInternalEnergy`, `solveTotalEnergy`,
  `solveTotalEnthalpy` (default). Implicit diffusion with `k/c_x` (c_x = c_p for the enthalpy forms, c_v for the energy forms) and a lagged
  explicit remainder `k grad T - (k/c_x) grad E`, so variable c_p and the kinetic-energy gradient are handled; convection is implicit upwind
  plus a deferred higher-order correction (`cSchemeEn`: uds, cds, sou, quick, muscl, anything else = linear upwind; `limiterEn`, `gdsEn`
  blending), non-conservative (default) or conservative (`conservativeEn = T`, which adds `m` to the diagonal and diverges in cold-start
  pseudo-time iterations; same converged solution). Further switches: `addViscDiss`, `cpModel` ('constant', 'polynomial', 'nasa', 'shomate'),
  `sutherland`/`mu_Tref`/`mu_S`, `urfEn`, `lSolverEn`, `maxiterEn`, `tolRelEn`. Temperature is recovered from the energy variable,
  density is updated. `calcEn` must stay `F`.
* **Bulk-viscosity term** `-(2/3) mu div(U)` is added explicitly to the momentum equations (`add_bulk_viscous_term` in `velocity.f90`) when
  `compressible = T`; it changes Cd of the NACA M = 0.5 case by less than 0.1 %.
* **Pseudo-time stabilisation** `cpl_pseudoCFL > 0` adds `(C Vol/dtau)(p - p*)` to the continuity rows with
  `dtau = CFL Vol / sum_f (|Vn|+c|S|)/2`. It vanishes at convergence and removes the start-up instability seen on finer meshes
  (starting impulsively from uniform flow). `cpl_cflRamp = n` grows the CFL geometrically from 1 to `cpl_pseudoCFL` over the first `n`
  iterations (needed for the impulsive start of transonic flow). Experience so far: subsonic M=0.5, 100-1000; transonic M=0.8, about 30 with
  a ramp of 150 iterations (CFL 200 diverged). The pseudo-time term lowers the convergence rate: do not use more than needed.
* **Sutherland's law** (`cpl_sutherland = T`, `cpl_Tref`, `cpl_Ssuth`): `viscos` is then the viscosity at `cpl_Tref`.

### Boundary conditions (compressible)

| patch type | velocity | pressure | temperature | mass flux |
|---|---|---|---|---|
| `inlet`    | prescribed (`0/U`) | zero gradient | prescribed (`0/T`, Dirichlet) | `p_P/(R T_b) U_b.S`, linear in p_P |
| `outlet`   | zero gradient | zero gradient, **not pinned** | zero gradient | `rho_P U_P.S`, linearised (supersonic or characteristic-free exit; no flux scaling) |
| `pressure` | zero gradient + Rhie-Chow correction | prescribed (`0/p`) | zero gradient | `rho_b [U_P.S + RC(p_b - p_P)]`, `rho_b = p_b/(R T_b)` |
| `wall`     | no slip | zero gradient (+ lagged linear extrapolation in the momentum force) | adiabatic (Neumann T) or isothermal (Dirichlet T) | 0 |
| `symmetry` | as in the incompressible solver | zero gradient | zero gradient | 0 |

Pressure level: fixed by a `pressure` patch; a closed steady domain is pinned at `pRefCell` with the value read from `0/p`.
Initial/boundary files `0/U`, `0/p`, `0/T` are all required.

#### Additional compressible boundary condition kinds (`implicit_bc_compressible.f90`)

A kind is attached to an existing patch by namelist arrays (the patch types in `polyMesh/boundary` do not change):

```
 cplbc_name(1) = 'inlet',   cplbc_kind(1) = 'stagnation', cplbc_p0(1) = 120183, cplbc_T0(1) = 300
 cplbc_name(2) = 'outlet',  cplbc_kind(2) = 'riemann'
 cplbc_urf = 0.3           ! Picard under-relaxation of the boundary state (default 0.3; 1.0 diverges for 'riemann')
```

| kind | patch type | prescribed | boundary state / flux |
|---|---|---|---|
| `stagnation` | `pressure` | `cplbc_p0`, `cplbc_T0`, direction (`cplbc_dir*`, default: normal) | inflow: `T_b = T0 - |U|^2/(2 c_p)`, `p_b = p0 (T_b/T0)^(g/(g-1))`; backflow: `p_b = p0` |
| `totalpressure` | `pressure` | `cplbc_p0`, static `T` of `0/T` | inflow: `p_b = p0 (1 + (g-1)/2 M^2)^(-g/(g-1))`, `M = |U|/c(T_b)` |
| `massflow` | `inlet` | `cplbc_mdot` (patch total, kg/s), direction | flux fixed, `U_b = mdot/(rho_b A)`, `rho_b = p_P/(R T_b)` |
| `compinlet` | `inlet` | `U`, `T` of the 0/ files | the default velocity inlet, density `p_b/(R T_b)` with the Picard-lagged `p_b` |
| `compoutlet` | `pressure` | static `p` (`cplbc_p0` or `0/p`) | subsonic: fixed p; supersonic (M>=1): extrapolated, linearised flux |
| `open` | `pressure` | `cplbc_p0`, `T` of `0/T` | outflow as `compoutlet`; inflow `p_b = p0 - rho u_n^2/2` (entrainment) |
| `riemann` | `pressure` | free stream `U`, `p`, `T` of the 0/ files | subsonic: linearised acoustic Riemann problem, `u_nb = (u_nP + u_n,inf + (p_P - p_inf)/(rho c))/2`, `p_b = (p_P + p_inf + rho c (u_nP - u_n,inf))/2` - continuity and the pressure force of the momentum equation are implicit in (U_P, p_P) (the cross-component coupling is lagged); entropy and tangential velocity from the upwind side; supersonic: free-stream state (inflow) or extrapolation (outflow) |

Mass fluxes of the pressure-type kinds use the fixed-p Rhie-Chow face of the `pressure` patch. The boundary state (`p_b`, `T_b`, `rho_b`) is updated
from the lagged cell state at the start of each outer iteration.

### Limitations

Constant (or Sutherland) viscosity with constant Prandtl number, no turbulence modelling consistent with variable density, no ideal-gas mixtures;
periodic boundaries and Crank-Nicolson are not supported. The kinds above are Picard-lagged (not Newton) except the Riemann kind; supersonic
inflow of `riemann` fixes only `p`, `T` (not `U`) of the patch. Known behaviour:

* A fixed-velocity/temperature inlet (`inlet`, `compinlet`) against a fixed-pressure outlet is only robust at low Mach number (the nozzle test
  converges at M_in = 0.1 and diverges at 0.3 and above, with and without the BC module, for any pseudo-time CFL tried); use `stagnation`, `totalpressure`
  or `massflow` inlets for internal flows at higher Mach number.
* `stagnation`/`totalpressure` are inflow conditions for ducts; as the far field of an external flow they give a non-physical solution (no suction peak),
  use `riemann` there.
* The pseudo-time CFL must not be too small for the nozzle tests: CFL 20 diverges where CFL 100-1000 converges (acoustic modes are only damped by the
  implicit pseudo-time step).
* Energy: `solveTotalEnergy` diverges from a cold (uniform) start; unlimited `sou` and `quick` diverge in the cold start (use `uds`, linear upwind, `muscl`).

### NACA 0012 test case

```
cd examples/coupledPU
python3 mkairfoil.py naca --nt 201 --nr 96 --y1 1e-3 --R 20 --M 0.5 --Re 5000 --alpha 0 --Pr 0.72
cd naca                                   # edit input-allspeeds.nml if wanted (cpl_pseudoCFL ...)
../../../bin/cappuccino input-allspeeds.nml monitor restart
python3 ../postairfoil.py . restart --cp cp.dat      # Cd, Cl, x_sep, Cp(x), Cf(x)
```

`mkairfoil.py` builds an O-grid (finite trailing-edge base, `--te`), writes the case (far field: upstream half `inlet`, downstream half
`pressure`; `--farfield riemann|riemann-in|riemann-out|open|compoutlet` selects the boundary condition kinds) and also `input-incomp.nml` (same flow, incompressible solver) for low-Mach comparisons (`--M 0.02`).

### Validation done during development

* Low Mach limit: NACA 0012, Re 5000, M = 0.02, 201x96 cells: Cd 0.05568 (all-speeds) vs 0.05568 (incompressible coupled solver),
  max |dCp| = 1.3e-4, max |dCf| = 3.7e-5.
* Energy equation: compressible Couette flow with viscous heating (`mkcouette.py`, `couette_check.py`): temperature profile error 1e-3 K on
  a 2 K heating amplitude, u error 7e-7 U; BDF1 and BDF2 transient runs relax to the same steady state; total temperature and entropy
  function are constant to 1e-4 in the inviscid region of the airfoil flow.
* NACA 0012, M = 0.5, Re = 5000, Pr = 0.72, alpha = 0 (Swanson & Langer, NASA/TM-2016-219003: Cd = 0.05557, pressure part 0.02279, friction
  part 0.03279, x_sep = 0.808): our Cd = 0.0636 / 0.0594 / 0.0583 on 121x56 / 201x96 / 301x128 cells (O-grid, R = 20 chords, first cell
  1e-3 .. 5e-4), i.e. approaching the reference from above, +5 % on the finest mesh; the separation bubble of the reference is not produced
  (the wall-adjacent tangential velocity goes to zero at x/c = 0.95 but does not reverse). The reference uses a sharp trailing edge,
  Sutherland's law, and a far-field vortex correction; our airfoil has a 6e-4 thick base.
* NACA 0012, M = 0.8, Re = 500, alpha = 10 deg (reference Cl = 0.436, Cd = 0.275, x_sep = 0.362): the solver converges (pseudo CFL 30 with ramp)
  to a steady transonic solution (max Mach 1.2), but the forces are not in agreement with the reference (Cl about 0.59, Cd about 0.32); neither the domain size (R = 20 or
  40) nor the viscosity law changes this. Mesh refinement does not close the gap either: 121x56 / 161x76 / 201x96 give Cd = 0.352 / 0.336 / 0.327 and
  Cl = 0.647 / 0.598 / 0.595, i.e. the solution converges to about Cd 0.32, Cl 0.59. Treat the transonic forces as unvalidated.
* Compressible boundary condition kinds, NACA 0012, M = 0.5, Re = 5000, 201x96 (default velocity inlet + pressure outlet: Cd = 0.05938):
  `riemann` on both patches Cd = 0.05929 (-0.15 %), `riemann-in` 0.05929, `riemann-out` / `open` / `compoutlet` 0.05938.
* Duct (`mknozzle.py`, `nozzle_check.py`, area ratio 0.8 at the throat, M_in = 0.5, 120x10 cells), against the quasi-1D isentropic solution:
  `stagnation` (inlet) max |dM| 0.003, `totalpressure` 0.003, `massflow` 0.006 (mass flow exact), `stagnation` + `compoutlet` / `open` 0.003,
  `stagnation` + `riemann` outlet 0.004; `compinlet` at M_in = 0.1: 0.001. (`python3 mknozzle.py nz --inlet stagnation --outlet riemann`,
  then run and `python3 nozzle_check.py nz restart`.)
