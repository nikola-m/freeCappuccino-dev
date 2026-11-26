# Quick reference

## Select a model
```
TurbModel = 'ebrsm'                     ! in the namelist; see table. Unknown name stops the run.
TurbModel%Scalar(1)%urf    = 0.5        ! per scalar slot: urf, gds (0 = first order convection), cScheme, lSolver, maxiter, tolAbs, tolRel
TurbModel%Scalar(1)%gds    = 0.0
TurbModel%urfVis           = 1.0        ! eddy viscosity relaxation
```

## Models (partial list), fields, meshes

| `TurbModel` | slots | needs 0/ files | y+ | notes |
|---|---|---|---|---|
| `k_epsilon_zeta_f` | 1 k, 2 eps, 3 zeta, 4 f | k, epsilon | ~1 | slot 4 (f): elliptic, urf 1; diverges on the airfoil start |
| `transition_sst` | 1 k, 2 omega, 3 gamma, 4 Re_theta | k, omega | ~1 | gamma/Re_theta start from 1/inlet correlation; `Scalar(4)%urf = 1` |
| `easm_rumsey_gatski_komega` | 1 k, 2 omega (= eps/k!) | k, omega | ~1 | omega 0.09x SST omega: `omega = eps/k` |
| `earsm_gatski_speziale` | 1 k, 2 eps | k, epsilon | 30-300 | wall functions |
| `rsm_full` | 1-6 R_ij, 7 eps | k, epsilon | 30-300 | stresses initialised 2/3 k delta_ij |
| `ebrsm` | 1-6 R_ij, 7 eps, 8 alpha | k, epsilon | < 1 | `Scalar(8)%urf = 1` |
| `sgs_one_equation_dynamic` | 1 k_sgs | k | wall resolved | transient LES, `ltransient = T` |

Inlet: `k = 1.5 (Tu U)^2`, `eps = Cmu k^2/(nut/nu nu)`, `omega = k/(nut/nu nu)` (SST units), `omega = eps/k` for the EASM-RG model.

## Practical settings
* Start every run with `gds = 0` and `urf = 0.5` on the turbulence scalars (with defaults `k_omega_sst` NaNs on the airfoil in 8 iterations).
* Stress models: `urf = 0.3-0.5` on slots 1-7. Elliptic scalars (`f`, `alpha`): urf 1.
* SIMPLE needs `ncorr = 1`, `npcor = 1` in the namelist (the airfoil input of the coupled generator has none).
* Restart does not hold the stresses, zeta, f, gamma, Re_theta, alpha: they are re-initialised on restart.
* VTK: k, epsilon/omega, viscosity; not the individual stresses.

## Commands
```
python3 mkchannel_turb.py ch --model <name> [--wallfunction]    # channel; run: cappuccino input-turb.nml monitor restart
python3 postchannel.py ch                                       # Cf, Dean
python3 mkairfoil_turb.py af --model <name>                     # NACA0012 Re=3e6; post: ../coupledPU/postairfoil.py af af/restart --pinf 0
python3 mkplate.py plate --model transition_sst; python3 platecf.py plate
```

## Reference numbers (channel Re_H = 1e4, Dean Cf = 0.00730)
SST 0.00679 | zeta-f 0.00677 | EASM-RG 0.00640 | EB-RSM 0.00745 | (y+ 28) k-eps 0.00662, EARSM 0.00649, RSM 0.00545 (known low).

## Add a model
`my_model.f90` (+ `case` in `turbulence_registry.f90`, + line in `src/cappuccino/files`). Use `transport_scalar` for each equation.
