#!/usr/bin/env python3
"""
mkairfoil_turb.py - turbulent flow around a NACA 0012 airfoil (incompressible, SIMPLE) for the RANS models of freeCappuccino.

  python3 mkairfoil_turb.py <dir> [--model k_omega_sst] [--Re 3e6 --alpha 0 --Tu 0.5 --nutr 10
                            --nt 201 --nw 60 --nr 80 --R 20 --y1 5e-6 --numstep 3000]

Builds the C-grid with ../coupledPU/mkairfoil.py (mesh, 0/U 0/p, info.txt), then adds the turbulence fields and writes
input-turb.nml (SIMPLE, TurbModel = <model>, under-relaxation suited to the model class).
Freestream: k = 1.5 (Tu U)^2, nut/nu = nutr, omega = k/(nutr nu), epsilon = Cmu k omega.
For 'easm_rumsey_gatski_komega' omega is eps/k (omega = Cmu k/(nutr nu)). Walls: k = 0, omega/epsilon are overwritten in the
wall cells by the models. Meshes: low-Re models (k_omega_sst, transition_sst, k_epsilon_zeta_f, ebrsm, easm_rumsey_gatski_komega)
need y+ < 1 (default y1 = 5e-6 at Re = 3e6); high-Re models (earsm_gatski_speziale, rsm_full) 30 < y+ < 300 (y1 ~ 3e-4).
Run:  cappuccino input-turb.nml monitor restart ;  post: python3 ../coupledPU/postairfoil.py <dir> restart --pinf 0
"""
import argparse, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
KOMEGA = ('k_omega_sst', 'transition_sst', 'easm_rumsey_gatski_komega')
LOWRE = KOMEGA + ('k_epsilon_zeta_f', 'ebrsm')


def scalar_file(root, name, val, wall):
    txt = ['internalField', '  uniform', '    %g' % val, 'boundaryField']
    for n in ('wall', 'far-in', 'far-out', 'front-back'):
        txt.append(n)
        if n == 'far-in':
            txt += ['  Dirichlet', '    uniform', '      %g' % val]
        elif n == 'wall':
            txt += ['  Dirichlet', '    uniform', '      %g' % wall]
        elif n == 'front-back':
            txt += ['  empty']
        else:
            txt += ['  Neumann', '    zeroGradient']
    open(os.path.join(root, '0', name), 'w').write('\n'.join(txt) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir'); ap.add_argument('--model', default='k_omega_sst')
    ap.add_argument('--Re', type=float, default=3e6); ap.add_argument('--alpha', type=float, default=0.0)
    ap.add_argument('--Tu', type=float, default=0.5, help='freestream turbulence intensity, percent')
    ap.add_argument('--nutr', type=float, default=10.0, help='freestream eddy viscosity ratio')
    ap.add_argument('--nt', type=int, default=201); ap.add_argument('--nw', type=int, default=60)
    ap.add_argument('--nr', type=int, default=80); ap.add_argument('--R', type=float, default=20.0)
    ap.add_argument('--y1', type=float, default=None)
    ap.add_argument('--numstep', type=int, default=3000)
    a = ap.parse_args()
    y1 = a.y1 if a.y1 is not None else (5e-6 if a.model in LOWRE else 3e-4)

    subprocess.check_call([sys.executable, os.path.join(HERE, '..', 'coupledPU', 'mkairfoil.py'), a.dir, '--grid', 'c',
                           '--nt', str(a.nt), '--nw', str(a.nw), '--nr', str(a.nr), '--R', str(a.R), '--y1', str(y1),
                           '--Re', str(a.Re), '--M', '0.15', '--alpha', str(a.alpha)])
    info = {l.split()[0]: float(l.split()[1]) for l in open(os.path.join(a.dir, 'info.txt'))}
    U = info['Uinf']; mu = info['mu']; rho = info['rhoinf']; nu = mu / rho
    k = 1.5 * (U * a.Tu / 100.0)**2
    om = k / (a.nutr * nu)                 # k-omega (SST units)
    cmu = 0.09
    eps = cmu * k * om
    if a.model == 'easm_rumsey_gatski_komega':
        om = eps / k                       # omega = eps/k
    scalar_file(a.dir, 'k', k, 1e-12)
    scalar_file(a.dir, 'omega', om, om)
    scalar_file(a.dir, 'epsilon', eps, eps)

    txt = open(os.path.join(a.dir, 'input-incomp.nml')).read()
    txt = re.sub(r"TurbModel = '[^']*'", "TurbModel = '%s'" % a.model, txt)
    txt = re.sub(r"SIMPLE = F", "SIMPLE = T\n              ncorr = 1\n              npcor = 1", txt); txt = re.sub(r"COUPLED = T", "COUPLED = F", txt)
    txt = re.sub(r"urfU = [0-9., ]+", "urfU = 0.7, 0.7, 0.7", txt); txt = re.sub(r"urfP = [0-9.]+", "urfP = 0.2", txt)
    txt = re.sub(r"NUMSTEP = \d+", "NUMSTEP = %d" % a.numstep, txt); txt = re.sub(r"NZAPIS = \d+", "NZAPIS = %d" % a.numstep, txt)
    # Turbulence scalars: first order convection (gds = 0) and urf 0.5 - the uniform initial field gives a huge startup
    # production at the wall, and second-order deferred correction of k/omega diverges within a few iterations.
    slots = {'k_omega_sst': (1, 2), 'k_epsilon_std': (1, 2), 'transition_sst': (1, 2, 3), 'easm_rumsey_gatski_komega': (1, 2),
             'k_epsilon_zeta_f': (1, 2, 3), 'earsm_gatski_speziale': (1, 2), 'rsm_full': range(1, 8), 'ebrsm': range(1, 8)}
    extra = ''.join("          TurbModel%%Scalar(%d)%%urf = 0.5\n          TurbModel%%Scalar(%d)%%gds = 0.0\n" % (i, i)
                    for i in slots.get(a.model, ()))
    if a.model == 'transition_sst':
        extra += "          TurbModel%Scalar(4)%urf = 1.0\n"
    if a.model == 'ebrsm':
        extra += "          TurbModel%Scalar(8)%urf = 1.0\n"
    txt = txt.replace("TurbModel = '%s'\n" % a.model, "TurbModel = '%s'\n%s" % (a.model, extra))
    open(os.path.join(a.dir, 'input-turb.nml'), 'w').write(txt)
    print('%s: model %s, U=%.2f nu=%.3e k=%g omega=%g eps=%g y1=%g' % (a.dir, a.model, U, nu, k, om, eps, y1))


if __name__ == '__main__':
    main()
