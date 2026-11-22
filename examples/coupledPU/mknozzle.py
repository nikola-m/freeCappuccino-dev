#!/usr/bin/env python3
"""
mknozzle.py - 2D converging-diverging duct (subsonic or transonic isentropic flow) for testing the compressible boundary
conditions of the all-speeds block-coupled p-U solver (implicit_bc_compressible.f90).

  python3 mknozzle.py <dir> [--inlet stagnation|totalpressure|massflow|compinlet] [--outlet pressure|compoutlet|open|riemann]
                      [--nx 120 --ny 10 --L 3 --H 1 --a 0.2 --Min 0.5 --T0 300 --pexit 101325 --numstep 3000]

Geometry: height h(x) = H (1 - a sin^2(pi x/L)), throat at x = L/2, the exit has the same area as the inlet. Both the axis
and the contour are slip ('symmetry') boundaries, no viscous effects (tiny viscosity), so the flow is the quasi-1D isentropic
flow of the area ratio h(x)/h_in.
Flow:  inlet Mach number Min is the design value; T_0, p_0 and the exit pressure follow from the isentropic relations
(p_exit = p_in because the exit area equals the inlet area), the mass flow of the 'massflow' inlet is the isentropic one.
Inlet kinds:  stagnation     p_0, T_0 prescribed  (patch type pressure, 'stagnation')
              totalpressure  p_0 prescribed, static T_in (patch type pressure, 'totalpressure')
              massflow       mdot prescribed, static T_in (patch type inlet, 'massflow')
              compinlet      U_in, static T_in (patch type inlet, 'compinlet')
Outlet kinds: pressure (default pressure patch), compoutlet, open, riemann (free stream = design exit state).
Check the result with nozzle_check.py.
"""
import argparse, os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkcase import build_mesh, write_native
from mkairfoil import NML

Rg, cp = 287.058, 1006.0
gam = cp / (cp - Rg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--inlet', default='stagnation'); ap.add_argument('--outlet', default='pressure')
    ap.add_argument('--nx', type=int, default=120); ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--L', type=float, default=3.0); ap.add_argument('--H', type=float, default=1.0)
    ap.add_argument('--a', type=float, default=0.2); ap.add_argument('--Min', type=float, default=0.5)
    ap.add_argument('--T0', type=float, default=300.0); ap.add_argument('--pexit', type=float, default=101325.0)
    ap.add_argument('--cfl', type=float, default=100.0); ap.add_argument('--ramp', type=int, default=100)
    ap.add_argument('--urfU', type=float, default=0.8); ap.add_argument('--urfEn', type=float, default=0.8)
    ap.add_argument('--numstep', type=int, default=3000)
    a = ap.parse_args()
    d = a.dir
    os.makedirs(os.path.join(d, 'vtk', 'mesh'), exist_ok=True)

    h = lambda x: a.H * (1 - a.a * np.sin(np.pi * x / a.L) ** 2)
    xn = np.linspace(0, a.L, a.nx + 1); yn = np.linspace(0, 1.0, a.ny + 1)
    mesh = build_mesh(xn, yn, 0.0)
    mesh['pts'] = [(x, y * h(x), z) for (x, y, z) in mesh['pts']]
    inl_type = 'inlet' if a.inlet in ('massflow', 'compinlet') else 'pressure'
    write_native(d, mesh, [('inlet', inl_type, 'left'), ('outlet', 'pressure', 'right'), ('wall-t', 'symmetry', 'top'),
                           ('wall-b', 'symmetry', 'bottom'), ('front-back', 'empty', 'frontback')])

    # isentropic design state
    f = 1 + 0.5 * (gam - 1) * a.Min ** 2
    Tin = a.T0 / f
    pexit = a.pexit
    pin = pexit                                  # same area
    p0 = pin * f ** (gam / (gam - 1))
    cin = math.sqrt(gam * Rg * Tin)
    Uin = a.Min * cin
    rhoin = pin / (Rg * Tin)
    thick = 0.1
    mdot = rhoin * Uin * h(0.0) * thick
    mu = 1e-5

    os.makedirs(os.path.join(d, '0'), exist_ok=True)
    pu = a.inlet in ('massflow', 'compinlet')
    with open(os.path.join(d, '0', 'U'), 'w') as fo:
        fo.write('internalField\n  uniform\n    %.10g 0 0\nboundaryField\n' % Uin)
        fo.write('inlet\n  Dirichlet\n    uniform\n      %.10g 0 0\n' % Uin)
        fo.write('outlet\n  Neumann\n    zeroGradient\n')
        for p in ('wall-t', 'wall-b', 'front-back'):
            fo.write('%s\n  Neumann\n    zeroGradient\n' % p)
    with open(os.path.join(d, '0', 'p'), 'w') as fo:
        fo.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % pexit)
        if pu:
            fo.write('inlet\n  Neumann\n    zeroGradient\n')
        else:
            fo.write('inlet\n  Dirichlet\n    uniform\n      %.10g\n' % pin)
        fo.write('outlet\n  Dirichlet\n    uniform\n      %.10g\n' % pexit)
        for p in ('wall-t', 'wall-b', 'front-back'):
            fo.write('%s\n  Neumann\n    zeroGradient\n' % p)
    with open(os.path.join(d, '0', 'T'), 'w') as fo:
        fo.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % Tin)
        fo.write('inlet\n  Dirichlet\n    uniform\n      %.10g\n' % Tin)
        fo.write('outlet\n  Neumann\n    zeroGradient\n')
        for p in ('wall-t', 'wall-b', 'front-back'):
            fo.write('%s\n  Neumann\n    zeroGradient\n' % p)

    ent = []
    kind_in = {'stagnation': 'stagnation', 'totalpressure': 'totalpressure', 'massflow': 'massflow', 'compinlet': 'compinlet'}[a.inlet]
    ent.append(("inlet", kind_in, p0 if a.inlet in ('stagnation', 'totalpressure') else 0.0,
                a.T0 if a.inlet == 'stagnation' else 0.0, mdot if a.inlet == 'massflow' else 0.0))
    if a.outlet != 'pressure':
        ent.append(("outlet", a.outlet, pexit if a.outlet != 'riemann' else 0.0, 0.0, 0.0))
    lines = ''
    for k, (nm, kd, p0v, T0v, md) in enumerate(ent, 1):
        lines += "         cplbc_name(%d) = '%s'\n         cplbc_kind(%d) = '%s'\n" % (k, nm, k, kd)
        if p0v: lines += '           cplbc_p0(%d) = %.10g\n' % (k, p0v)
        if T0v: lines += '           cplbc_T0(%d) = %.10g\n' % (k, T0v)
        if md: lines += '         cplbc_mdot(%d) = %.10g\n' % (k, md)
    nml = NML.format(M=a.Min, Re=0, alpha=0, mu='%.10g' % mu, Pr=0.72, urfU=a.urfU, urfEn=a.urfEn, blend=1.0,
                     numstep=a.numstep, maxit=40, cfl=a.cfl, ramp=a.ramp, suth='F')
    nml = nml.replace("title = 'NACA0012 M=%s Re=0 alpha=0 (all-speeds block-coupled p-U)'" % a.Min,
                      "title = 'nozzle %s/%s'" % (a.inlet, a.outlet))
    nml = nml.replace('           pRefCell = 1', lines + '           pRefCell = 1')
    with open(os.path.join(d, 'input-allspeeds.nml'), 'w') as fo:
        fo.write(nml)
    with open(os.path.join(d, 'info.txt'), 'w') as fo:
        fo.write('nx %d\nny %d\nL %.10g\nH %.10g\na %.10g\nMin %.10g\nT0 %.10g\np0 %.10g\npexit %.10g\nmdot %.10g\nthick %.10g\n'
                 % (a.nx, a.ny, a.L, a.H, a.a, a.Min, a.T0, p0, pexit, mdot, thick))
    print(d, 'p0 = %.1f Pa, T0 = %.1f K, T_in = %.2f K, U_in = %.3f m/s, mdot = %.5f kg/s' % (p0, a.T0, Tin, Uin, mdot))


if __name__ == '__main__':
    main()
