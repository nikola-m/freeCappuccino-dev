#!/usr/bin/env python3
"""
mkchannel_turb.py - developing plane channel (H = 1, U_b = 1, nu = 1e-4, Re_H = 1e4, Lx = 100) for the RANS models.
The quickest regression/verification test of a turbulence model: Cf in the developed part against Dean, Cf = 0.073 Re_H^-0.25.

  python3 mkchannel_turb.py <dir> [--model k_omega_sst] [--wallfunction] [--nx 200 --ny 50 --stretch 2.5 --numstep 4000]

Default mesh 200x50, tanh stretching, y+ < 1 (low-Re models). --wallfunction: uniform 100x10 mesh, y+ ~ 25 (earsm_gatski_speziale,
rsm_full, or any k-epsilon model). Inlet and initial k = 0.01, epsilon = 9e-4 (nut/nu = 100); omega = eps/(0.09 k) = 1, or eps/k for
easm_rumsey_gatski_komega. Run:  cappuccino input-turb.nml monitor restart ;  post:  python3 postchannel.py <dir>
Needs mkcase.py of ../coupledPU.
"""
import argparse, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def field(root, name, val, wall):
    txt = ['internalField', '  uniform', '    %g' % val, 'boundaryField',
           'inlet', '  Dirichlet', '    uniform', '      %g' % val,
           'outlet', '  Neumann', '    zeroGradient']
    for w in ('wall-t', 'wall-b'):
        txt += [w, '  Dirichlet', '    uniform', '      %g' % wall]
    txt += ['front-back', '  empty']
    open(os.path.join(root, '0', name), 'w').write('\n'.join(txt) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir'); ap.add_argument('--model', default='k_omega_sst')
    ap.add_argument('--wallfunction', action='store_true')
    ap.add_argument('--nx', type=int, default=None); ap.add_argument('--ny', type=int, default=None)
    ap.add_argument('--stretch', type=float, default=None); ap.add_argument('--Lx', type=float, default=100.0)
    ap.add_argument('--numstep', type=int, default=4000)
    a = ap.parse_args()
    nx = a.nx or (100 if a.wallfunction else 200); ny = a.ny or (10 if a.wallfunction else 50)
    st = a.stretch if a.stretch is not None else (0.0 if a.wallfunction else 2.5)
    subprocess.check_call([sys.executable, os.path.join(HERE, '..', 'coupledPU', 'mkcase.py'), 'channel', a.dir,
                           '--nx', str(nx), '--ny', str(ny), '--Lx', str(a.Lx), '--Re', '1e4', '--stretch', str(st),
                           '--outlet', 'outlet'])
    os.makedirs(os.path.join(a.dir, 'vtk'), exist_ok=True)
    k0, e0 = 0.01, 9e-4
    om = e0 / k0 if a.model == 'easm_rumsey_gatski_komega' else e0 / (0.09 * k0)
    field(a.dir, 'k', k0, 0.01 if a.wallfunction else 1e-10)
    field(a.dir, 'omega', om, om)
    field(a.dir, 'epsilon', e0, e0)
    txt = open(os.path.join(a.dir, 'input-simple.nml')).read()
    txt = re.sub(r"TurbModel = '[^']*'", "TurbModel = '%s'" % a.model, txt)
    txt = re.sub(r"viscos = [0-9.eE+-]+", "viscos = 1e-4", txt)
    txt = re.sub(r"NUMSTEP = \d+", "NUMSTEP = %d" % a.numstep, txt); txt = re.sub(r"NZAPIS = \d+", "NZAPIS = %d" % a.numstep, txt)
    if a.model in ('ebrsm', 'rsm_full'):
        extra = ''.join("          TurbModel%%Scalar(%d)%%urf = 0.5\n" % i for i in range(1, 8))
        txt = txt.replace("TurbModel = '%s'\n" % a.model, "TurbModel = '%s'\n%s" % (a.model, extra))
    open(os.path.join(a.dir, 'input-turb.nml'), 'w').write(txt)
    print('%s: %s, %dx%d, stretch %g' % (a.dir, a.model, nx, ny, st))


if __name__ == '__main__':
    main()
