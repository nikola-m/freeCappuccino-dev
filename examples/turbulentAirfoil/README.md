# turbulentAirfoil - test bench for the new turbulence-model layer

Cases, generators and post-processors that were used to verify the turbulence models added in patches 0031-0040
(`src/finiteVolume/fvEqnDiscretization/TurbulenceModels`). Despite the folder name the first thing to run is the **channel**, which takes
a minute per model and gives a quantitative number (Cf against Dean's correlation); the **airfoil** is the demonstration case; the
**flat plate** is for the transition model.

| file | what |
|---|---|
| `mkchannel_turb.py`, `postchannel.py` | developing channel Re_H = 1e4, low-Re mesh 200x50 (y+<1) or wall-function mesh 100x10 (y+~25); Cf versus Dean |
| `mkairfoil_turb.py` | NACA 0012, C-grid, Re = 3e6, alpha = 0, incompressible SIMPLE; uses `../coupledPU/mkairfoil.py` for the mesh |
| `mkplate.py`, `platecf.py` | ERCOFTAC T3A-like flat plate (Tu = 3.3 %) for `transition_sst`; Cf versus Blasius and the turbulent correlation |
| `README.md`, `TECHNICAL_DOCUMENTATION.md`, `QUICK_REFERENCE.md` | this file; theory, implementation, verification; one-page cheat sheet |

Postprocessing of the airfoil: `../coupledPU/postairfoil.py <dir> restart --pinf 0` (Cd, Cl, Cp, Cf). The generators need the python
helpers of `../coupledPU` (`mkcase.py`, `vtkmesh.py`), numpy.

## Models (partial list - see `src/finiteVolume/fvEqnDiscretization/TurbulenceModels/turbulence_registry.f90` for the full list) 

`TurbModel` in the namelist selects the model by name (checked once after reading the input, unknown name stops the run):

| `TurbModel` | closure | wall treatment / mesh | momentum coupling |
|---|---|---|---|
| `k_epsilon_zeta_f` | Hanjalic-Popovac-Hadziabdic k-eps-zeta-f (rewritten, was NaN) | low-Re, y+ ~ 1 (log-law branch if y+ > 11.63) | eddy viscosity |
| `transition_sst` | Langtry-Menter gamma-Re_theta_t on SST | low-Re, y+ ~ 1 | eddy viscosity |
| `easm_rumsey_gatski_komega` | Gatski-Speziale EASM, Rumsey-Gatski k-omega (omega = eps/k) | low-Re, y+ ~ 1 | `aij` source (`lasm`) |
| `earsm_gatski_speziale` | Gatski-Speziale EASM, standard k-eps | wall functions, 30 < y+ < 300 | `aij` source |
| `rsm_full` | SSG Reynolds-stress transport + eps, six stress equations | wall functions, 30 < y+ < 300 | stresses via artificial viscosity + `aij` |
| `ebrsm` | Manceau-Hanjalic elliptic-blending RSM (SSG + eps + alpha) | low-Re, y+ < 1 | as `rsm_full` |
| `sgs_one_equation_dynamic` | Kim-Menon localized dynamic one-equation SGS (LES) | wall-resolved only | eddy viscosity |

Other existing models (`k_epsilon_std`, `k_omega_sst`, `Spalart_Allmaras`, ...) are unchanged and selected the same way.

## Quick start

```
make                                           # in the repository root -> bin/cappuccino
cd examples/turbulentAirfoil

# channel (a few minutes)
python3 mkchannel_turb.py ch --model ebrsm
(cd ch && ../../../bin/cappuccino input-turb.nml monitor restart)
python3 postchannel.py ch

# airfoil (a few minutes)
python3 mkairfoil_turb.py af --model k_omega_sst
(cd af && ../../../bin/cappuccino input-turb.nml monitor restart)
python3 ../coupledPU/postairfoil.py af af/restart --pinf 0

# transition on the flat plate
python3 mkplate.py plate --model transition_sst
(cd plate && ../../../bin/cappuccino input.nml monitor restart)
python3 platecf.py plate
```

Use `--wallfunction` with `mkchannel_turb.py` (and for the airfoil the generator picks `y1 = 3e-4` automatically) for
`earsm_gatski_speziale` and `rsm_full`; the low-Re models need the fine mesh. Without the right mesh class a model runs but is wrong.

## Results (what was actually run)

Channel, developing flow, x = 98 (Re_H = 1e4, Dean Cf = 0.00730):

| model | mesh | Cf | vs Dean |
|---|---|---|---|
| `k_omega_sst` (existing) | 200x50 | 0.00679 | -7 % |
| `k_epsilon_zeta_f` | 200x50 | 0.00677 | -7 % |
| `easm_rumsey_gatski_komega` | 200x50 | 0.00640 | -12 % |
| `ebrsm` | 200x50 | 0.00745 | +2 % |
| `k_epsilon_std` (existing) | 100x10, y+ 28 | 0.00662 | -9 % |
| `earsm_gatski_speziale` | 100x10, y+ 28 | 0.00649 | -11 % |
| `rsm_full` | 100x10, y+ 25 | 0.00545 | **-25 %** (see limitations) |

Flat plate T3A-like (`transition_sst`): Cf follows Blasius (about 1.15-1.2x at Re_x = 1.2e5, the usual finite-thickness effect), minimum
near Re_x = 1.4e5, fully turbulent correlation (0.00433 at x = 1.1 m against 0.0043) by Re_x of 3-4e5.

NACA 0012, Re = 3e6, alpha = 0, C-grid 201 x (60+60 wake) x 80, y1 = 5e-6 (high-Re models: y1 = 3e-4):

| model | mesh | Cd (pressure + viscous) | comment |
|---|---|---|---|
| `k_omega_sst` (existing) | y1 = 5e-6 | 0.00907 (0.00166 + 0.00741) | reference run of the series |
| `easm_rumsey_gatski_komega` | y1 = 5e-6 | 0.00912 | within 0.5 % of SST |
| `ebrsm` | y1 = 5e-6 | 0.01022 | +13 % against SST, not investigated |
| `transition_sst` | y1 = 5e-6 | 0.00507 | free-stream Tu = 0.5 %: laminar run, not compared with data |
| `k_epsilon_zeta_f` | y1 = 5e-6 | **diverges** | NaN after ~150 iterations from the impulsive start, see limitations; fine in the channel |
| `k_epsilon_std` (existing) | y1 = 3e-4, y+ ~ 18 | 0.00660 | baseline for the wall-function mesh |
| `earsm_gatski_speziale` | y1 = 3e-4 | 0.00634 | -4 % against k-eps on the same mesh |
| `rsm_full` | y1 = 3e-4 | 0.00565 | -14 % against k-eps on the same mesh |

|Cl| < 1e-4 for every run that converged (symmetric case). 3000 SIMPLE iterations each; the monitored residual of U fell to 1e-6 or below
for the low-Re runs, the high-Re runs were not examined as carefully. No comparison with experiment was made - the wall-function mesh is
coarse (y+ ~ 18, at the lower edge of the log layer), so its absolute Cd (0.0066 against 0.0091) says more about the mesh than the models.

LES: `sgs_one_equation_dynamic` on the 3D lid-driven cavity Re = 1e4, 24^3, 1500 steps of 0.02: stable, mean Ck = 0.046, mean Ceps = 0.4.
There is no LES validation (no DNS comparison) - the model runs and behaves plausibly, nothing more.

## Limitations you should know about

* `rsm_full` gives Cf 25 % too low in the channel on a wall-function mesh. EB-RSM, with the same stress equations, is +2 %, so the
  stress core is fine and the fault is in the high-Re wall treatment (wall cell stresses prescribed from the log-layer ratios).
  Not resolved. Treat `rsm_full` as a working framework, not a validated model.
* `k_epsilon_zeta_f` diverges on the airfoil (impulsive start: zeta goes negative in the first iteration because the large negative `f` of the wall cells
  enters the zeta equation as an explicit source). A production limiter `P <= 10 rho eps` and `zeta <= 2` were added; this was not enough. Treating the negative part of
  `f` implicitly keeps the airfoil stable but moves the channel Cf from -7 % to -17 %, or diverges, depending on the floor of zeta; both variants were rejected. Starting from
  a converged `k_omega_sst` solution is untested. Fine in the channel and (patch 0033) it no longer gives NaN there.
* Plain `k_omega_sst` also diverges on the T3A plate (`mkplate.py --model k_omega_sst`, NaN) in master and in the new build; the cause was not found. `transition_sst` and `k_epsilon_zeta_f` run on it.
* Diffusion of the stresses and of epsilon is Daly-Harlow contracted with the gradient direction (scalar form), exact for 1D gradients only.
* On the airfoil mesh the turbulence scalars need first-order convection (`TurbModel%Scalar(i)%gds = 0`) and urf 0.5 to start: with the defaults
  (second-order deferred correction, urf 0.7) even the *existing* `k_omega_sst` gives NaN at iteration 8 (startup shock at the wall from the
  uniform initial field; k reaches 70 in the first iteration). The two changes were made together, so which one matters was not isolated.
  The generators write these settings. A run with `gds = 1` after convergence was not tried.
* Restart files do not store the fields owned by the new models (stresses, zeta, f, gamma, ...): a restarted run starts those from the
  isotropic/initial values and relaxes in a few hundred iterations.
* `rsm_full`/`ebrsm`/EASM write `k` and `epsilon`/`omega` to VTK, not the individual stress components.
* Not done: MPI, wall-function LES, the original Rumsey-Gatski cross-diffusion variants, a dynamic model for Prandtl numbers.
