#!/usr/bin/env python3
"""
mkplate.py - zero-pressure-gradient flat plate (ERCOFTAC T3A-like) for the transition / RANS models.

  python3 mkplate.py <dir> [--model k_omega_sst|transition_sst|k_epsilon_zeta_f] [--nx 400 --ny 70 --Lx 1.6
                           --H 0.2 --U 5.4 --nu 1.5e-5 --Tu 3.3 --nutr 12 --y1 6e-5]

Domain x in [-0.05, Lx]; bottom: symmetry for x<0 (lead-in), wall for x>=0; top: symmetry; inlet; outlet = fixed pressure.
Writes native mesh, 0/U 0/k 0/omega 0/epsilon, input.nml. Run:  cappuccino input.nml monitor restart
Post-process with platecf.py (skin friction versus Blasius / turbulent flat plate correlations).
Needs mkcase.py of ../coupledPU (mesh builder and native writer).
"""
import argparse, os, sys, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'coupledPU'))
from mkcase import build_mesh, write_native, write_U, write_input


def geometric_lines(H, n, y1):
    """n cells, first cell height y1, constant growth ratio r such that the sum is H."""
    lo, hi = 1.0 + 1e-9, 2.0
    f = lambda r: y1 * (r**n - 1.0) / (r - 1.0) - H
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    r = 0.5 * (lo + hi)
    return np.concatenate(([0.0], np.cumsum(y1 * r**np.arange(n)))), r


def scalar_file(root, name, val, patches, wall=None):
    out = ['internalField', '  uniform', '    %g' % val, 'boundaryField']
    for n, t in patches:
        out.append(n)
        if t == 'inlet':
            out += ['  Dirichlet', '    uniform', '      %g' % val]
        elif t == 'wall' and wall is not None:
            out += ['  Dirichlet', '    uniform', '      %g' % wall]
        elif t == 'empty':
            out += ['  empty']
        else:
            out += ['  Neumann', '    zeroGradient']
    open(os.path.join(root, '0', name), 'w').write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--model', default='transition_sst')
    ap.add_argument('--nx', type=int, default=400); ap.add_argument('--nlead', type=int, default=10)
    ap.add_argument('--ny', type=int, default=70)
    ap.add_argument('--Lx', type=float, default=1.6); ap.add_argument('--H', type=float, default=0.2)
    ap.add_argument('--U', type=float, default=5.4); ap.add_argument('--nu', type=float, default=1.5e-5)
    ap.add_argument('--Tu', type=float, default=3.3, help='free-stream turbulence intensity in percent')
    ap.add_argument('--nutr', type=float, default=12.0, help='free-stream eddy viscosity ratio')
    ap.add_argument('--y1', type=float, default=6e-5)
    ap.add_argument('--numstep', type=int, default=4000)
    a = ap.parse_args()

    os.makedirs(os.path.join(a.dir, 'vtk', 'mesh'), exist_ok=True)
    xn = np.concatenate((np.linspace(-0.05, 0.0, a.nlead + 1)[:-1], np.linspace(0.0, a.Lx, a.nx + 1)))
    yn, r = geometric_lines(a.H, a.ny, a.y1)
    mesh = build_mesh(xn, yn)
    # split bottom boundary into lead-in (symmetry) and plate (wall)
    xc = 0.5 * (xn[:-1] + xn[1:])
    bot = mesh['bnd']['bottom']
    mesh['bnd']['lead'] = [b for b, x in zip(bot, xc) if x < 0.0]
    mesh['bnd']['plate'] = [b for b, x in zip(bot, xc) if x >= 0.0]
    patches = [('inlet', 'inlet', 'left'), ('outlet', 'pressure', 'right'), ('top', 'symmetry', 'top'),
               ('lead', 'symmetry', 'lead'), ('plate', 'wall', 'plate'), ('front-back', 'empty', 'frontback')]
    info = write_native(a.dir, mesh, patches)
    write_U(a.dir, (a.U, 0, 0), [('inlet', (a.U, 0, 0)), ('outlet', None), ('top', None), ('lead', None),
                                 ('plate', (0, 0, 0)), ('front-back', None)])
    k = 1.5 * (a.U * a.Tu / 100.0)**2
    om = k / (a.nutr * a.nu)
    eps = 0.09 * k * om
    pt = [(n, t) for n, t, _ in patches]
    scalar_file(a.dir, 'k', k, pt, wall=1e-12)
    scalar_file(a.dir, 'omega', om, pt)
    scalar_file(a.dir, 'epsilon', eps, pt)
    extra = ("          TurbModel = '%s'\n" % a.model +
             "          TurbModel%Scalar(4)%urf = 1.0\n")
    write_input(a.dir, 'input.nml', title='flat plate T3A-like', nu=a.nu, simple='T', numstep=a.numstep,
                urfU=0.6, urfP=0.2, tol='1e-9')
    txt = open(os.path.join(a.dir, 'input.nml')).read().replace("          TurbModel = 'none'\n", extra)
    open(os.path.join(a.dir, 'input.nml'), 'w').write(txt)
    print(a.dir, info, 'growth ratio %.4f' % r, 'k=%g omega=%g eps=%g' % (k, om, eps))


if __name__ == '__main__':
    main()
