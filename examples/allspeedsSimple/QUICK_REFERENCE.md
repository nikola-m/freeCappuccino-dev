# All-speeds SIMPLE - quick reference

## Switch and recipes
```
 SIMPLE = T   ASIMPLE = T   npcor = 1       ! required (npcor has no default in the code!)
 compressible = T  calcEn = T               ! all speeds; compressible = F -> classical SIMPLE with AMG
```
| situation | settings |
|-----------|----------|
| incompressible (check against stock SIMPLE) | `compressible = F`, `urfU = 0.8`, `asimple_urfP = 0.3` |
| compressible, transient or dual time (robust) | `ltransient = T  bdf = T (or bdf2 = T)  timestep = dt  maxit = 4..6`, `asimple_urfP = 0.3..0.7`; NACA 7200 cells: dt = 1e-4 s with 6 outer its was stable (other dt not explored) |
| compressible, steady | `ltransient = F`, `asimple_pseudoCFL = 3..10`, `asimple_cflRamp = 100..200`, `asimple_urfP = 0.3..0.7`, `urfU = 0.5..0.8` |
| walls with y+ > 1 | `asimple_wallFunc = 'reichardt'` (or `standard` with a k-eps model, `spl`, `unified`, `compressible`) |
| inlet/outlet of internal compressible flow | `cplbc_name/kind/p0/T0/mdot` as in the coupled solver: `stagnation`, `totalpressure`, `massflow` inlets; `compoutlet`, `open`, `riemann` outlets; `asimple_urfP = 0.2`, fixed-U inlet only at M_in <= 0.2 |
| transonic / supersonic | `asimple_beta = 0` (upwind density), `asimple_pseudoCFL = 3 .. 10` with `asimple_cflRamp = 200`, `cSchemeU/En = 'vanleer'`; see `mkbump.py` |

Input files names, monitor and restart names: at most 20 characters (otherwise `End of file` in `read_input`).

## Parameters (namelist `input_parameters`)
| parameter | default | meaning |
|-----------|---------|---------|
| `asimple` | F | use the all-speeds SIMPLE instead of SIMPLE/PISO/coupled |
| `asimple_urfP` | 0.3 | under-relaxation of the pressure update `p = p + urf p'` (velocities use the full p') |
| `asimple_beta` | 1.0 | weight of the central face density, `beta = asimple_beta max(0, 1-M_f^2)`; 0 = pure upwind density |
| `asimple_lSolverP` | `'amg'` | `'amg'` = own AMG + FGMRES, or any `csrsolve` name (`iccg`, `bicgstab`, ...) |
| `asimple_maxitP`, `asimple_tolRelP`, `asimple_tolAbsP` | 40, 0.05, 1e-14 | iterations, relative drop, absolute level of the p' solver |
| `asimple_pseudoCFL`, `asimple_cflRamp` | 0, 0 | steady compressible: pseudo-time term in the p' equation (0 = off); ramp length in iterations |
| `asimple_wallFunc` | `'none'` | `none`, `standard`, `reichardt`, `spl` (Shih et al.), `unified` (Ong & Chan), `compressible` |
| `pamg_theta`, `pamg_coarse_size`, `pamg_nu1`, `pamg_nu2`, `pamg_max_levels`, `pamg_restart`, `pamg_verbose` | 0.25, 40, 1, 1, 25, 20, F | AMG strength threshold, direct solve below this size, smoothing sweeps, levels, GMRES restart, print hierarchy |

`ltest = T` prints the wall function statistics (y+, `visw/viscos`), the mass flow through inlet/outlet/pressure patches and `d(mass)/dt`.
Convergence monitor line of the p' solver: `AMG-FGMRES: Solving for p, Initial residual = .., Final residual = .., No Iterations ..`;
the initial residual is the L1 mass imbalance and is stored in `resor(4)`.

## Boundary conditions used
`inlet` - U and T prescribed, flux `rho_b U_b.S` (explicit) | `outlet` - zero gradient, `rho_P U_P.S` with `p'` coupling (compressible) or flux scaling (incompressible)
| `pressure` - fixed p, Rhie-Chow boundary flux | `wall`, `symmetry`, `empty` - no flux | `periodic` - not supported (stop).

## Wall functions
| value | law | needs |
|-------|-----|-------|
| `standard` | Launder-Spalding, `u* = Cmu^1/4 sqrt(k)` | k-based turbulence model |
| `reichardt` | Reichardt's law, all y+ | nothing |
| `spl` | Shih, Povinelli, Liu (2003): `U_t/u_c = ... f1(Y_tau+) + (nu/rho) G f2(Y_c+)/u_c^3`, `u_c = u_tau + u_p`, fits of Tables 1-2, all signs of the pressure gradient | cell pressure gradient |
| `unified` | Ong & Chan (2018): `u_t = u1 + u2`, Newton for `u2+`, Spalding-type law with compressibility and heat transfer; isothermal wall if `T` of the wall is Dirichlet in `0/T`, else adiabatic; `compressible = F` -> incompressible limit | cell pressure gradient, `T`, `p` |
| `compressible` | van Driest transformed Reichardt, Crocco-Busemann, Kader heat flux (diagnostic `qwall`) | `compressible = T`, wall T boundary value |

The result is `visw`, used by the existing momentum equation as `max(viscos, visw)`.

## Files
`mkcase.py naca|nozzle|bump|cavity|channel <dir>` - cases (`mkbump.py`, `bump_check.py`: circular-arc bump channel and its post-processing); `run_validation.sh [bin] [workdir]` - all tests; postprocessing scripts are in `../coupledPU`.
