#!/usr/bin/env python3
"""
mkairfoil.py - O-grid around a NACA 00xx airfoil (2D, one cell thick, 'empty' front/back) with a case set-up for
the all-speeds block-coupled p-U solver (compressible = T, coupled = T), in the native freeCappuccino mesh format.

  python3 mkairfoil.py <dir> [--nt 201 --nr 96 --R 20 --y1 1e-3 --M 0.5 --Re 5000 --alpha 0 --Pr 0.72 --te 3e-4 --cfl 200 --ramp 0 --sutherland]

Mesh   : nt cells around (odd), nr cells in the radial direction, geometric growth chosen to reach the far field at
         radius R (chords) from the first cell height y1; wall points cluster at LE and TE (cosine spacing).
         Near the wall the grid lines follow the wall normal, further out they become radial rays.
         The trailing edge has a small base of half thickness --te (0 = would be sharp; keep > 0, O-grid needs a base cell).
Patches: wall (airfoil), inlet (upstream half of the far field: T and U prescribed), outlet (downstream half,
         type 'pressure': p prescribed, zero gradient for U and T), front-back (empty).
Flow   : p_inf = 101325 Pa, T_inf = 288.15 K, air (R = 287.058, c_p = 1006, Pr as given), U_inf from the Mach number,
         mu from the Reynolds number (chord = 1):  mu = rho_inf U_inf / Re.
Writes polyMesh/, 0/U 0/p 0/T, input-allspeeds.nml (steady, coupled), input-incomp.nml (the same flow incompressible, for
low-Mach comparisons) and info.txt (reference quantities for postprocessing).
"""
import argparse, os, sys, math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkcase import write_native

Rgas, cp = 287.058, 1006.0


def airfoil_nodes(N, t=0.12, te=3e-4):
    """N+1 nodes: lower TE -> LE -> upper TE, cosine spacing. Half thickness at x=1 equals te."""
    a4 = 0.1036 - te / (5 * t)            # a4 = 0.1036 gives a sharp TE; smaller -> finite base
    th = 2 * np.pi * np.arange(N + 1) / N
    x = 0.5 * (1 + np.cos(th))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - a4 * x**4)
    y = np.where(th <= np.pi, -yt, yt)
    return np.column_stack([x, y])


def radial_distribution(nr, y1n):
    """s_0=0 ... s_nr=1 with first spacing y1n and constant growth ratio g."""
    lo, hi = 1.0 + 1e-9, 2.0
    f = lambda g: y1n * (g**nr - 1) / (g - 1) - 1
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        (lo, hi) = (mid, hi) if f(mid) < 0 else (lo, mid)
    g = 0.5 * (lo + hi)
    s = np.concatenate([[0], np.cumsum(y1n * g ** np.arange(nr))])
    return s / s[-1], g


def build(nt, nr, R, y1, te, alpha):
    N = nt - 1
    assert N % 2 == 0, 'nt must be odd'
    P = airfoil_nodes(N, te=te)                       # (nt,2) wall nodes, ring order
    c0 = np.array([0.5, 0.0])
    # wall normals (outward from the body), from the central-difference tangent
    tg = np.roll(P, -1, axis=0) - np.roll(P, 1, axis=0)
    tg /= np.linalg.norm(tg, axis=1)[:, None]
    nrm = np.column_stack([tg[:, 1], -tg[:, 0]])
    if np.dot(nrm[N // 2], P[N // 2] - c0) < 0:       # at the LE the normal must point forward (-x)
        nrm = -nrm
    # far-field points: blend of the geometric angle and a uniform angle (clockwise ring order, wake at angle 0)
    # monotone clockwise: lower surface angles go 0- -> -pi, upper surface (wrapped) -pi -> -2pi
    phi_g = np.zeros(nt)
    for i in range(nt):
        a = math.atan2(P[i, 1] - c0[1], P[i, 0] - c0[0])
        phi_g[i] = a if i <= N // 2 else a - 2 * np.pi
    phi_u = -(2 * np.pi / nt) * (np.arange(nt) + 0.5)
    phi = 0.5 * phi_g + 0.5 * phi_u
    Pf = c0 + R * np.column_stack([np.cos(phi), np.sin(phi)])

    sj, g = radial_distribution(nr, y1 / R)
    X = np.zeros((nr + 1, nt, 2))
    for i in range(nt):
        ray = Pf[i] - P[i]
        L = np.linalg.norm(ray)
        e = ray / L
        for j in range(nr + 1):
            t = sj[j] * L
            w = 1.0 - math.exp(-(t / 1.0) ** 2)       # normal near the wall -> ray far away
            d = (1 - w) * nrm[i] + w * e
            d /= np.linalg.norm(d)
            X[j, i] = P[i] + t * d
        X[nr, i] = Pf[i]
    return X, g


def hex_mesh(X, thick=0.1, alpha=0.0):
    nr1, nt, _ = X.shape
    nr = nr1 - 1
    nid = lambda i, j, k: 1 + (i % nt) + nt * (j + nr1 * k)
    pts = [(X[j, i, 0], X[j, i, 1], thick * k) for k in range(2) for j in range(nr1) for i in range(nt)]
    cid = lambda i, j: 1 + (i % nt) + nt * j
    P = np.array(pts)
    cen = {}
    for j in range(nr):
        for i in range(nt):
            ids = [nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i, j + 1, 0),
                   nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]
            cen[cid(i, j)] = P[np.array(ids) - 1].mean(0)

    def orient(nodes, c_from):
        p = P[np.array(nodes) - 1]
        c = p.mean(0)
        n = np.zeros(3)
        for q in range(len(nodes)):
            n += np.cross(p[q] - c, p[(q + 1) % len(nodes)] - c) / 2
        return nodes if np.dot(n, c - c_from) > 0 else nodes[::-1]

    inner, bnd = [], {'wall': [], 'far-in': [], 'far-out': [], 'frontback': []}
    for j in range(nr):
        for i in range(nt):
            c = cid(i, j)
            # face between (i,j) and (i+1,j): nodes at ring index i+1
            nb = cid(i + 1, j)
            nodes = [nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i + 1, j + 1, 1), nid(i + 1, j, 1)]
            L, Rr = min(c, nb), max(c, nb)
            inner.append((L, Rr, orient(nodes, cen[L])))
            # radial faces
            nodes = [nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]
            if j < nr - 1:
                nb = cid(i, j + 1)
                inner.append((c, nb, orient(nodes, cen[c])))
            else:
                fc = P[np.array(nodes) - 1].mean(0)
                key = 'far-in' if (fc[0] - 0.5) * math.cos(alpha) + fc[1] * math.sin(alpha) < 0 else 'far-out'
                bnd[key].append((c, orient(nodes, cen[c])))
            if j == 0:
                nodes = [nid(i, 0, 0), nid(i + 1, 0, 0), nid(i + 1, 0, 1), nid(i, 0, 1)]
                bnd['wall'].append((c, orient(nodes, cen[c])))
            # front/back
            bnd['frontback'].append((c, orient([nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i, j + 1, 0)], cen[c])))
            bnd['frontback'].append((c, orient([nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)], cen[c])))
    inner.sort(key=lambda t: (t[0], t[1]))
    return dict(pts=[tuple(p) for p in pts], inner=inner, bnd=bnd)


NML = """&input_parameters
              title = 'NACA0012 M={M} Re={Re} alpha={alpha} (all-speeds block-coupled p-U)'
        mesh_format = 'nativeMesh'
              lread = F
              ltest = F
              calcU = T
               urfU = {urfU}, {urfU}, {urfU}
               gdsU = 1.0
           cSchemeU = 'linearUpwind'
           dSchemeU = 'skewness'
            nrelaxU = 1
           lSolverU = 'bicgstab'
           maxiterU = 5
            tolAbsU = 1e-13
            tolRelU = 0.01
            pscheme = 'linear'
              calcP = T
               urfP = 0.3
           lSolverP = 'iccg'
           maxiterP = 50
            tolAbsP = 1e-12
            tolRelP = 0.025
              calcT = F
             calcEn = F
              urfEn = {urfEn}
           cSchemeEn = 'linearUpwind'
           lSolverEn = 'bicgstab'
           maxiterEn = 10
           tolRelEn = 0.01
         addViscDiss = T
          TurbModel = 'none'
             densit = 1.225
             viscos = {mu}
              pranl = {Pr}
               beta = 1e-3
         ltransient = F
                bdf = F
               bdf2 = F
               bdf3 = F
                 CN = F
              lstsq = F
           lstsq_qr = F
           lstsq_dm = F
         node_gauss = F
              gauss = T
            limiter = 'no-limit'
             SIMPLE = F
               PISO = F
          compressible = T
            COUPLED = T
         cpl_maxit = {maxit}
        cpl_restart = 30
        cpl_tolRel = 0.1
        cpl_rhoBlend = {blend}
        cpl_pseudoCFL = {cfl}
        cpl_cflRamp = {ramp}
        cpl_sutherland = {suth}
           pRefCell = 1
          tolerance = 1e-10
            NUMSTEP = {numstep}
           TIMESTEP = 1.0
             NZAPIS = {numstep}
              MAXIT = 1
/
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--nt', type=int, default=201); ap.add_argument('--nr', type=int, default=96)
    ap.add_argument('--R', type=float, default=20.0); ap.add_argument('--y1', type=float, default=1e-3)
    ap.add_argument('--M', type=float, default=0.5); ap.add_argument('--Re', type=float, default=5000.0)
    ap.add_argument('--alpha', type=float, default=0.0); ap.add_argument('--Pr', type=float, default=0.72)
    ap.add_argument('--te', type=float, default=3e-4)
    ap.add_argument('--urfU', type=float, default=0.7); ap.add_argument('--urfEn', type=float, default=0.8)
    ap.add_argument('--blend', type=float, default=1.0)
    ap.add_argument('--cfl', type=float, default=None, help='pseudo-time CFL (default 200 for M<=0.6, 30 above)')
    ap.add_argument('--ramp', type=int, default=None, help='iterations of the CFL ramp (default 0, 150 for M>0.6)')
    ap.add_argument('--sutherland', action='store_true', help='Sutherland viscosity law')
    ap.add_argument('--numstep', type=int, default=3000); ap.add_argument('--maxit', type=int, default=40)
    ap.add_argument('--farfield', default='default', choices=['default', 'riemann', 'riemann-in', 'riemann-out', 'open', 'compoutlet'],
                    help='boundary conditions of the far field (implicit_bc_compressible.f90): default = velocity inlet '
                         '(far-in) + static pressure outlet (far-out); riemann = both patches (-in / -out: one of them); '
                         'open, compoutlet = far-out')
    a = ap.parse_args()
    os.makedirs(os.path.join(a.dir, 'vtk', 'mesh'), exist_ok=True)

    al = math.radians(a.alpha)
    cfl = a.cfl if a.cfl is not None else (200.0 if a.M <= 0.6 else 30.0)
    ramp = a.ramp if a.ramp is not None else (0 if a.M <= 0.6 else 150)
    suth = 'T' if a.sutherland else 'F'
    X, g = build(a.nt, a.nr, a.R, a.y1, a.te, al)
    mesh = hex_mesh(X, 0.1, al)
    in_type = 'pressure' if a.farfield in ('riemann', 'riemann-in') else 'inlet'
    info = write_native(a.dir, mesh, [('wall', 'wall', 'wall'), ('far-in', in_type, 'far-in'),
                                      ('far-out', 'pressure', 'far-out'), ('front-back', 'empty', 'frontback')])

    pinf, Tinf = 101325.0, 288.15
    gamma = cp / (cp - Rgas)
    rho = pinf / (Rgas * Tinf)
    c = math.sqrt(gamma * Rgas * Tinf)
    U = a.M * c
    mu = rho * U / a.Re
    ux, uy = U * math.cos(al), U * math.sin(al)
    os.makedirs(os.path.join(a.dir, '0'), exist_ok=True)
    with open(os.path.join(a.dir, '0', 'U'), 'w') as f:
        f.write('internalField\n  uniform\n    %.10g %.10g 0\nboundaryField\n' % (ux, uy))
        f.write('wall\n  Dirichlet\n    uniform\n      0 0 0\n')
        f.write('far-in\n  Dirichlet\n    uniform\n      %.10g %.10g 0\n' % (ux, uy))
        if a.farfield in ('riemann', 'riemann-out', 'open'):
            f.write('far-out\n  Dirichlet\n    uniform\n      %.10g %.10g 0\n' % (ux, uy))
        else:
            f.write('far-out\n  Neumann\n    zeroGradient\n')
        f.write('front-back\n  Neumann\n    zeroGradient\n')
    with open(os.path.join(a.dir, '0', 'p'), 'w') as f:
        f.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % pinf)
        f.write('wall\n  Neumann\n    zeroGradient\n')
        if in_type == 'pressure':
            f.write('far-in\n  Dirichlet\n    uniform\n      %.10g\n' % pinf)
        else:
            f.write('far-in\n  Neumann\n    zeroGradient\n')
        f.write('far-out\n  Dirichlet\n    uniform\n      %.10g\n' % pinf)
        f.write('front-back\n  Neumann\n    zeroGradient\n')
    with open(os.path.join(a.dir, '0', 'T'), 'w') as f:
        f.write('internalField\n  uniform\n    %.10g\nboundaryField\n' % Tinf)
        f.write('wall\n  Neumann\n    zeroGradient\n')
        f.write('far-in\n  Dirichlet\n    uniform\n      %.10g\n' % Tinf)
        if a.farfield in ('riemann', 'riemann-out', 'open'):
            f.write('far-out\n  Dirichlet\n    uniform\n      %.10g\n' % Tinf)
        else:
            f.write('far-out\n  Neumann\n    zeroGradient\n')
        f.write('front-back\n  Neumann\n    zeroGradient\n')
    nml = NML.format(M=a.M, Re=a.Re, alpha=a.alpha, mu='%.10g' % mu, Pr=a.Pr, urfU=a.urfU, urfEn=a.urfEn,
                     blend=a.blend, numstep=a.numstep, maxit=a.maxit, cfl=cfl, ramp=ramp, suth=suth)
    cpl = {'default': [], 'riemann': [('far-in', 'riemann', 0, 0), ('far-out', 'riemann', 0, 0)],
           'riemann-in': [('far-in', 'riemann', 0, 0)], 'riemann-out': [('far-out', 'riemann', 0, 0)],
           'open': [('far-out', 'open', pinf, 0)],
           'compoutlet': [('far-out', 'compoutlet', pinf, 0)]}[a.farfield]
    lines = ''
    for k, (nm, kd, p0v, T0v) in enumerate(cpl, 1):
        lines += "         cplbc_name(%d) = '%s'\n         cplbc_kind(%d) = '%s'\n" % (k, nm, k, kd)
        if p0v: lines += '           cplbc_p0(%d) = %.10g\n' % (k, p0v)
        if T0v: lines += '           cplbc_T0(%d) = %.10g\n' % (k, T0v)
    nml = nml.replace('           pRefCell = 1', lines + '           pRefCell = 1')
    with open(os.path.join(a.dir, 'input-allspeeds.nml'), 'w') as f:
        f.write(nml)
    # Same flow, incompressible (density rho_inf, p_out = 0): for low-Mach consistency checks (use --M 0.02)
    inc = NML.format(M=a.M, Re=a.Re, alpha=a.alpha, mu='%.10g' % mu, Pr=a.Pr, urfU=0.95, urfEn=a.urfEn,
                     blend=a.blend, numstep=a.numstep, maxit=a.maxit, cfl=0, ramp=0, suth='F')
    inc = inc.replace('densit = 1.225', 'densit = %.10g' % rho).replace('compressible = T', 'compressible = F')
    inc = inc.replace('all-speeds block-coupled', 'incompressible block-coupled')
    with open(os.path.join(a.dir, 'input-incomp.nml'), 'w') as f:
        f.write(inc)
    with open(os.path.join(a.dir, 'info.txt'), 'w') as f:
        f.write('pinf %.10g\nTinf %.10g\nrhoinf %.10g\nUinf %.10g\nalpha %.10g\nmu %.10g\nchord 1\nthickness_z 0.1\n'
                % (pinf, Tinf, rho, U, a.alpha, mu))
    print(a.dir, info, 'growth ratio %.4f' % g, 'U_inf %.3f m/s, mu %.5g Pa s, rho %.4f' % (U, mu, rho))


if __name__ == '__main__':
    main()
