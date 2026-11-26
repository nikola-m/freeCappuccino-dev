#!/usr/bin/env python3
"""
mkbump.py - circular-arc bump channel (GAMM / ONERA-type test) for the all-speeds SIMPLE solver, inviscid limit.

  python3 mkbump.py <dir> [--M 0.5 --t 0.10 --nx 180 --ny 60 --scheme vanleer --urfP 0.2 --beta 0.5 --cfl 10 --numstep 6000]

2D channel 0 < x < 3, height 1, circular-arc bump of chord 1 (1 < x < 2) and relative thickness t on the lower wall; the walls
are slip ('symmetry' patches), the viscosity is tiny, so the flow is the Euler flow of the GAMM test case.  Classical set-ups:

   M = 0.50   t = 0.10    subsonic, isentropic (the entropy production is a measure of the numerical dissipation)
   M = 0.85   t = 0.042   transonic: supersonic pocket terminated by a shock on the bump
   M = 1.40   t = 0.04    supersonic: oblique shocks from the leading and the trailing edge of the bump

Boundary conditions (implicit_bc_compressible kinds): M < 1: 'stagnation' inlet (p_0, T_0) and 'compoutlet' (static p = p_inf);
M > 1: velocity inlet ('compinlet', U and T of the free stream, p extrapolated) and 'compoutlet' (extrapolated at supersonic exit).
Writes input-asimple.nml (all-speeds SIMPLE, steady, pseudo-time) and input-coupled.nml (block-coupled solver, for comparison).
Check the result with bump_check.py.
"""
import argparse, os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'coupledPU'))
from mkcase import build_mesh, write_native
from mkairfoil import NML

Rg, cp = 287.058, 1006.0
gam = cp / (cp - Rg)


def bump_y(x, t):
    """lower wall: circular arc of chord 1 between x = 1 and 2, thickness t"""
    R = (0.25 + t * t) / (2 * t)
    x = np.asarray(x, float)
    return np.where((x > 1) & (x < 2), np.sqrt(np.maximum(R * R - (x - 1.5) ** 2, 0)) - (R - t), 0.0)


def grid(nx, ny, L=3.0):
    """uniform streamwise nodes, uniform wall-normal parameter (the mesh is mapped onto the bump)"""
    return np.linspace(0, L, nx + 1), np.linspace(0, 1, ny + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--M', type=float, default=0.5); ap.add_argument('--t', type=float, default=None)
    ap.add_argument('--nx', type=int, default=180); ap.add_argument('--ny', type=int, default=60)
    ap.add_argument('--T0', type=float, default=300.0); ap.add_argument('--pinf', type=float, default=101325.0)
    ap.add_argument('--scheme', default=None, help="convection scheme of U and T ('linearUpwind' for M<0.7, 'vanleer' above)")
    ap.add_argument('--urfP', type=float, default=0.2); ap.add_argument('--beta', type=float, default=None)
    ap.add_argument('--cfl', type=float, default=None); ap.add_argument('--ramp', type=int, default=200)
    ap.add_argument('--urfU', type=float, default=0.7); ap.add_argument('--urfEn', type=float, default=0.8)
    ap.add_argument('--numstep', type=int, default=6000)
    ap.add_argument('--wall', action='store_true', help='lower wall = no-slip wall (with --wf a wall function), upper wall stays slip')
    ap.add_argument('--wf', default='unified', help="wall function of the no-slip wall: none, reichardt, spl, unified, compressible")
    ap.add_argument('--Tw', type=float, default=0.0, help='isothermal wall temperature [K] (0 = adiabatic)')
    ap.add_argument('--mu', type=float, default=1e-5, help='dynamic viscosity (1.8e-5 for air)')
    ap.add_argument('--stretch', type=float, default=0.0, help='geometric stretching of the wall-normal grid towards the lower wall')
    a = ap.parse_args()
    d = a.dir
    if a.t is None:
        a.t = 0.10 if a.M < 0.7 else (0.042 if a.M < 1.0 else 0.04)
    if a.scheme is None:
        a.scheme = 'linearUpwind' if a.M < 0.7 else 'vanleer'
    if a.beta is None:
        a.beta = 0.5 if a.M < 0.7 else 0.0       # density weight of the central part: subsonic 0.5, transonic/supersonic upwind
    if a.cfl is None:
        a.cfl = 3.0 if 0.7 <= a.M < 1.0 else 10.0
    os.makedirs(os.path.join(d, 'vtk', 'mesh'), exist_ok=True)

    xn, eta = grid(a.nx, a.ny)
    if a.stretch > 0.0:
        eta = (np.exp(a.stretch * eta) - 1.0) / (np.exp(a.stretch) - 1.0)
    mesh = build_mesh(xn, eta, 0.0)
    pts = []
    for (x, y, z) in mesh['pts']:
        yb = float(bump_y(x, a.t))
        pts.append((x, yb + y * (1.0 - yb), z))
    mesh['pts'] = pts
    sup = a.M > 1.0
    write_native(d, mesh, [('inlet', 'inlet' if sup else 'pressure', 'left'), ('outlet', 'pressure', 'right'),
                           ('wall-t', 'symmetry', 'top'), ('wall-b', 'wall' if a.wall else 'symmetry', 'bottom'), ('front-back', 'empty', 'frontback')])

    f = 1 + 0.5 * (gam - 1) * a.M ** 2
    Tinf = a.T0 / f
    p0 = a.pinf * f ** (gam / (gam - 1))
    Uinf = a.M * math.sqrt(gam * Rg * Tinf)

    os.makedirs(os.path.join(d, '0'), exist_ok=True)
    order = ['inlet', 'outlet', 'wall-t', 'wall-b', 'front-back']      # the order of polyMesh/boundary (the files are read in this order)
    zg = '  Neumann\n    zeroGradient\n'

    def dirichlet(val):
        return '  Dirichlet\n    uniform\n      %s\n' % val

    def write_field(name, internal, special):
        with open(os.path.join(d, '0', name), 'w') as fo:
            fo.write('internalField\n  uniform\n    %s\nboundaryField\n' % internal)
            for pn in order:
                fo.write(pn + '\n' + special.get(pn, zg))

    write_field('U', '%.10g 0 0' % Uinf, {'inlet': dirichlet('%.10g 0 0' % Uinf), 'wall-b': dirichlet('0 0 0') if a.wall else zg})
    write_field('p', '%.10g' % a.pinf, {'inlet': zg if sup else dirichlet('%.10g' % a.pinf), 'outlet': dirichlet('%.10g' % a.pinf)})
    write_field('T', '%.10g' % Tinf, {'inlet': dirichlet('%.10g' % Tinf),
                                      'wall-b': dirichlet('%.10g' % a.Tw) if (a.wall and a.Tw > 0.0) else zg})

    lines = "         cplbc_name(1) = 'inlet'\n"
    if sup:
        lines += "         cplbc_kind(1) = 'compinlet'\n"
    else:
        lines += "         cplbc_kind(1) = 'stagnation'\n           cplbc_p0(1) = %.10g\n           cplbc_T0(1) = %.10g\n" % (p0, a.T0)
    lines += "         cplbc_name(2) = 'outlet'\n         cplbc_kind(2) = 'compoutlet'\n           cplbc_p0(2) = %.10g\n" % a.pinf

    nml = NML.format(M=a.M, Re=0, alpha=0, mu='%g' % a.mu, Pr=0.72, urfU=a.urfU, urfEn=a.urfEn, blend=1.0,
                     numstep=a.numstep, maxit=40, cfl=a.cfl, ramp=a.ramp, suth='F')
    nml = nml.replace("title = 'NACA0012 M=%s Re=0 alpha=0 (all-speeds block-coupled p-U)'" % a.M,
                      "title = 'bump M=%s t=%s'" % (a.M, a.t))
    nml = nml.replace('           pRefCell = 1', lines + '           pRefCell = 1')
    nml = nml.replace("cSchemeU = 'linearUpwind'", "cSchemeU = '%s'" % a.scheme)
    nml = nml.replace("cSchemeEn = 'linearUpwind'", "cSchemeEn = '%s'" % a.scheme)
    with open(os.path.join(d, 'input-coupled.nml'), 'w') as fo:
        fo.write(nml)

    asn = nml.replace('SIMPLE = F', 'SIMPLE = T').replace('COUPLED = T', 'COUPLED = F').replace('calcEn = F', 'calcEn = T')
    asn = asn.replace('MAXIT = 1\n/', 'MAXIT = 1\n  ASIMPLE = T\n  npcor = 1\n  asimple_urfP = %g\n  asimple_beta = %g\n  asimple_pseudoCFL = %g\n  asimple_cflRamp = %d\n/'
                      % (a.urfP, a.beta, a.cfl, a.ramp))
    if a.wall:
        # the stop test of the code looks at the momentum residuals only, which are tiny in the (nearly) inviscid flow
        asn = asn.replace('tolerance = 1e-10', 'tolerance = 1e-14')
        asn = asn.replace('  ASIMPLE = T\n', "  ASIMPLE = T\n  asimple_wallFunc = '%s'\n" % a.wf)
    with open(os.path.join(d, 'input-asimple.nml'), 'w') as fo:
        fo.write(asn)
    with open(os.path.join(d, 'info.txt'), 'w') as fo:
        fo.write('nx %d\nny %d\nM %.10g\nt %.10g\nT0 %.10g\np0 %.10g\npinf %.10g\nTinf %.10g\nUinf %.10g\n'
                 % (a.nx, a.ny, a.M, a.t, a.T0, p0, a.pinf, Tinf, Uinf))
    print(d, 'M = %g, t = %g, p0 = %.1f Pa, T_inf = %.2f K, U_inf = %.2f m/s' % (a.M, a.t, p0, Tinf, Uinf))


if __name__ == '__main__':
    main()
