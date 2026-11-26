# All-speeds SIMPLE (alternative pressure-velocity coupling)

`simple_allspeeds` is a second, independent implementation of the SIMPLE algorithm that works from incompressible to
supersonic flow. It follows the idea of Przulj & Basara: pressure-based, collocated, ideal gas `rho = p/(R T)`, pressure-correction
equation with a transient and a convective density-pressure term, momentum -> p' -> corrections -> energy. It lives next to the
existing SIMPLE, PISO and block-coupled solvers and **does not change them**: with `ASIMPLE = F` (default) the code runs exactly as before.
It is meant for testing and comparison.

| part | file | what |
|------|------|------|
| driver, p' equation | `src/finiteVolume/fvEqnDiscretization/SimpleAllSpeeds/simple_allspeeds.f90` | module `simple_allspeeds` |
| AMG for p' | `.../SimpleAllSpeeds/asimple_amg.f90` | module `asimple_amg` (scalar, aggregation, FGMRES; `block_amg` is not touched) |
| wall functions | `.../SimpleAllSpeeds/asimple_wallfun.f90` | module `asimple_wallfun`: standard, Reichardt, Shih et al. (2003) generalized (`spl`), Ong & Chan (2018) unified (`unified`), van Driest/Kader (`compressible`) |
| boundary conditions | `Pressure/implicit_bc_compressible.f90` (existing) | the compressible BC kinds of the namelist (`cplbc_name/cplbc_kind/cplbc_p0/...`) are shared with the coupled solver: `bc_setup`, `bc_prepare_face`; the face assembly of the segregated solver is in `simple_allspeeds` |
| hooks | `main.f90`, `read_input.f90`, `files`, `directories` | one `elseif ( asimple )` branch, namelist entries |

What is reused from the code: momentum assembly and solution (`calcuvw`, including stress-tensor diffusion, BDF1/2, bulk viscosity
term), Rhie-Chow flux (`facefluxmass2`), pressure boundary fluxes, `gradp_and_sources`, the energy equation and the time levels
(`energy_equation`, `allspeeds`). Documents: `TECHNICAL_DOCUMENTATION.md` (method, derivation, limits), `QUICK_REFERENCE.md` (parameters, recipes).

## Build and run
```sh
make                                                   # repository root -> bin/cappuccino
python3 examples/allspeedsSimple/mkcase.py naca case   # also: nozzle, bump, cavity, channel
cd case && ../bin/cappuccino input-asimple.nml mon rst # use <= 20 characters for the file names!
python3 ../examples/coupledPU/postairfoil.py . rst     # Cd, Cl, ...
sh examples/allspeedsSimple/run_validation.sh          # all of the below
```
The case generators are those of `examples/coupledPU`; `mkcase.py` only derives the input files. Switch on in any `input.nml`:
```
 SIMPLE = T   ASIMPLE = T   npcor = 1        ! all-speeds SIMPLE instead of calcp_simple
 compressible = T                            ! ideal gas, energy equation is solved by main (calc_energy)
 asimple_urfP = 0.3 .. 0.7                   ! pressure under-relaxation
```
Compressible boundary conditions are attached to the patches exactly as in the coupled solver (see `examples/coupledPU/README.md`):
```
 cplbc_name(1) = 'inlet'   cplbc_kind(1) = 'stagnation'  cplbc_p0(1) = 120183  cplbc_T0(1) = 300
 cplbc_name(2) = 'outlet'  cplbc_kind(2) = 'compoutlet'  cplbc_p0(2) = 101325     ! kinds: stagnation, totalpressure, massflow, compinlet, compoutlet, open, riemann
```
With `compressible = F` the density is constant and the algorithm reduces to classical SIMPLE (useful check, see below).

## Results (verified in this session, serial, gfortran -O3)
| case | result |
|------|--------|
| lid-driven cavity 32x32, Re 100, incompressible | stock SIMPLE 391 iterations, all-speeds SIMPLE with `iccg` 391 (identical), with own AMG 390 |
| subsonic nozzle (60x6, M_in = 0.2, steady, pseudo CFL 10 ramp 100) | vs. quasi-1D isentropic: max abs(M-M_1D) = 0.0013 (M_throat 0.2539 vs 0.2536), p error 5e-4, T error 1e-4, 3000 iterations ~1 s |
| NACA0012, M 0.5, Re 5000, coarse C-grid 101x40x40 (7200 cells), BDF1 dual time, dt = 1e-4 s, 6 outer its | Cd = 0.05639 (pressure 0.02370, viscous 0.03269), Cl = 0 within 1e-5, momentum residual 1e-9, no NaN; 2 min |
| the same, existing block-coupled solver (steady, 3000 its) | Cd = 0.05765, 1 min. Reference (Swanson & Langer) 0.05557 is for a fine mesh: the mesh is coarse, this is NOT a validation of either number |
| the same with BDF2 (`bdf2 = T`) | Cd = 0.05637 |
| channel, Re 20000, `asimple_wallFunc='reichardt'` | y+ = 39..63, visw/viscos 3.28 (wall function active); converges (239 its vs 202 without) |
| the same channel, `'spl'` (Shih et al.) / `'unified'` (Ong & Chan, incompressible limit) | visw/viscos 3.60 / 3.33, 247 / 234 iterations. `spl` is about 10 % above Reichardt because its log law has C = 5.0 (Reichardt: about 5.5) |
| NACA0012 as above with `'spl'` or `'unified'` (y+ < 1) | Cd identical to `none` in 5 digits (0.05639): the laws reduce to the viscous sublayer |
| nozzle, M_in = 0.3, BC kinds `stagnation`, `totalpressure`, `massflow` inlet; `compoutlet`, `open`, `riemann` outlet; `asimple_urfP = 0.2` | max abs(M-M_1D) 0.0015 / 0.0015 / 0.0019; 0.0015 / 0.0015 / 0.0018 (`riemann` in and out 0.0022) |
| nozzle, M_in = 0.5, `stagnation` inlet | `asimple_beta = 1` (default) does not converge (limit cycle, M_throat 0.81 vs 0.733); `asimple_beta = 0.5`: M_throat 0.7264, max abs(M-M_1D) 0.010 |
| bump channel, M = 0.5, t = 10 % (180x60, `stagnation`/`compoutlet`) | max Mach 0.688 (coupled solver 0.689), wall p/p_inf min 0.857 (0.857), entropy s = p/rho^g/(p/rho^g)_inf - 1: mean 1.0e-4, max 8.7e-3 at the leading/trailing edge |
| bump, M = 0.85, t = 4.2 % (transonic, `asimple_beta = 0`, CFL 3) | supersonic pocket and shock at x = 1.89-1.93 (90 % chord): max Mach 1.332 (coupled 1.345), wall p/p_inf min 0.546 (0.535); the run stops at the iteration limit (6000); the shock position moves by a cell, the p' residual stalls at 0.03-0.1 (not converged to machine precision) |
| bump, M = 1.4, t = 4 % (supersonic, `compinlet`/`compoutlet`) | converged to tolerance in 2537 iterations; wall p/p_inf at the leading edge 1.58 (sharp wedge of 9.15 deg: oblique shock theory 1.68 - the arc expands the flow at once, not comparable to better than ~10 %), max Mach 1.65. The coupled solver was started on the same mesh but did not finish in the time limit: no comparison |

Known limits: see section 9 of the technical documentation (no periodic boundaries, no CN/BDF3 with `compressible`, steady NACA stalls
at a residual of 1e-2, fixed-U inlets only for low Mach numbers, the compressible wall functions were verified in a stand-alone test but NOT in a compressible CFD case, heat flux of the wall functions not passed to the energy equation).
