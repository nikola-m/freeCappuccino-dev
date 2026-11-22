#!/usr/bin/env python3
"""
mkcouette.py - compressible planar Couette flow with viscous heating: a check of the energy equation of the all-speeds solver.

  python3 mkcouette.py <dir> [--U 150 --nx 20 --ny 40 --Lx 2 --h 1 --Pr 0.72 --mu 1.8e-3]

Both walls are isothermal (T_w); the upper wall moves with U. For constant properties the exact solution is
    u = U y/h,     T = T_w + Pr U^2/(2 c_p) (y/h)(1 - y/h),     p = const
which is prescribed at the inlet; the outlet is a fixed-pressure patch. A correct solver reproduces the profile
along the channel: compare the temperature at the outlet column with the exact one using couette_check.py.
"""
import argparse, os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkcase import build_mesh, write_native

cp, Rg = 1006.0, 287.058


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir'); ap.add_argument('--U', type=float, default=150.0)
    ap.add_argument('--nx', type=int, default=20); ap.add_argument('--ny', type=int, default=40)
    ap.add_argument('--Lx', type=float, default=2.0); ap.add_argument('--h', type=float, default=1.0)
    ap.add_argument('--Pr', type=float, default=0.72); ap.add_argument('--mu', type=float, default=1.8e-3)
    ap.add_argument('--numstep', type=int, default=400)
    a = ap.parse_args()
    d = a.dir
    os.makedirs(os.path.join(d, 'vtk', 'mesh'), exist_ok=True)
    xn = np.linspace(0, a.Lx, a.nx + 1); yn = np.linspace(0, a.h, a.ny + 1)
    mesh = build_mesh(xn, yn, 0.0)
    write_native(d, mesh, [('inlet', 'inlet', 'left'), ('outlet', 'pressure', 'right'), ('wall-t', 'wall', 'top'),
                           ('wall-b', 'wall', 'bottom'), ('front-back', 'empty', 'frontback')])
    Tw, p0 = 300.0, 101325.0
    # table for 'interpolated' inlet values; the code assumes an implicit first point (y=0, value 0), so start at y=1e-9
    eta = np.concatenate([[1e-9], np.linspace(0, 1, 41)[1:]])
    heat = a.Pr * a.U**2 / (2 * cp)
    os.makedirs(os.path.join(d, '0'), exist_ok=True)
    with open(os.path.join(d, '0', 'U'), 'w') as f:
        f.write('internalField\n  uniform\n    %g 0 0\nboundaryField\n' % (0.5 * a.U))
        f.write('inlet\n  Dirichlet\n    interpolated\n      %d\n' % len(eta))
        for e in eta:
            f.write('      %.10g %.10g 0 0\n' % (e * a.h, a.U * e))
        f.write('outlet\n  Neumann\n    zeroGradient\n')
        f.write('wall-t\n  Dirichlet\n    uniform\n      %g 0 0\n' % a.U)
        f.write('wall-b\n  Dirichlet\n    uniform\n      0 0 0\nfront-back\n  Neumann\n    zeroGradient\n')
    with open(os.path.join(d, '0', 'p'), 'w') as f:
        f.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % p0)
        f.write('inlet\n  Neumann\n    zeroGradient\noutlet\n  Dirichlet\n    uniform\n      %.10g\n' % p0)
        f.write('wall-t\n  Neumann\n    zeroGradient\nwall-b\n  Neumann\n    zeroGradient\nfront-back\n  Neumann\n    zeroGradient\n')
    with open(os.path.join(d, '0', 'T'), 'w') as f:
        f.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % Tw)
        f.write('inlet\n  Dirichlet\n    interpolated\n      %d\n' % len(eta))
        for e in eta:
            f.write('      %.10g %.10g\n' % (e * a.h, Tw + heat * e * (1 - e)))
        f.write('outlet\n  Neumann\n    zeroGradient\nwall-t\n  Dirichlet\n    uniform\n      %g\n' % Tw)
        f.write('wall-b\n  Dirichlet\n    uniform\n      %g\nfront-back\n  Neumann\n    zeroGradient\n' % Tw)
    # reuse the airfoil namelist template with different physics
    from mkairfoil import NML
    nml = NML.format(M=0, Re=0, alpha=0, mu='%.10g' % a.mu, Pr=a.Pr, urfU=0.8, urfEn=0.9, blend=1.0,
                     numstep=a.numstep, maxit=40, cfl=100, ramp=0, suth='F').replace("NACA0012 M=0 Re=0 alpha=0 (all-speeds block-coupled p-U)",
                                                          "compressible Couette flow")
    open(os.path.join(d, 'input-couette.nml'), 'w').write(nml)
    open(os.path.join(d, 'info.txt'), 'w').write('U %g\nh %g\nTw %g\nPr %g\nnx %d\nny %d\nLx %g\n' % (a.U, a.h, Tw, a.Pr, a.nx, a.ny, a.Lx))
    print(d, 'viscous heating amplitude  Pr U^2/(2 cp) = %.4f K' % heat)


if __name__ == '__main__':
    main()
