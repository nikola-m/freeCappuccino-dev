#!/usr/bin/env python3
"""
postairfoil.py - forces, Cp and separation point from a restart file of the airfoil case (mkairfoil.py).

  python3 postairfoil.py <casedir> <restart-file> [--cp cp.dat] [--pinf 0]

Force on the body from the wall faces (patch 'wall'):  F = sum (p_b - p_inf) S  +  mu_w |S| U_t,P / d_n
(pressure force plus viscous drag from the cell-centre tangential velocity; U_wall = 0, wall at distance d_n).
Coefficients use q_inf = rho_inf U_inf^2 / 2, chord 1 and the span of the mesh (z thickness in info.txt).
Cp(x), the skin friction C_f(x) and the separation point (first sign change of tangential velocity from the TE
side on the upper surface, x_sep) are written to --cp.
"""
import sys, struct, argparse, math
import numpy as np


def records(fn):
    b = open(fn, 'rb').read(); o = 0; out = []
    while o < len(b):
        n = struct.unpack('<i', b[o:o + 4])[0]; out.append(b[o + 4:o + 4 + n]); o += 8 + n
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case'); ap.add_argument('restart'); ap.add_argument('--cp', default=None)
    ap.add_argument('--pinf', type=float, default=None, help='reference pressure (default: p_inf of info.txt; use 0 for incompressible runs)')
    a = ap.parse_args()
    d = a.case + '/polyMesh/'
    info = {l.split()[0]: float(l.split()[1]) for l in open(a.case + '/info.txt')}
    P = np.loadtxt(d + 'points')
    size = [int(x) for x in open(d + 'size').read().split()]
    npts, nc, nI, nB, nF = size
    faces = [[int(t) - 1 for t in l.split()[1:]] for l in open(d + 'faces')]
    own = np.loadtxt(d + 'owner', dtype=int) - 1
    bnd = [l.split() for l in open(d + 'boundary') if not l.startswith('#')]
    wall = [b for b in bnd if b[1] == 'wall'][0]
    nw, sw = int(wall[2]), int(wall[3])

    cc = np.zeros((nc, 3)); cnt = np.zeros(nc)
    for f, nodes in enumerate(faces):
        c = P[nodes].mean(0); cc[own[f]] += c; cnt[own[f]] += 1
    nb = np.loadtxt(d + 'neighbour', dtype=int) - 1
    for f in range(nI):
        cc[nb[f]] += P[faces[f]].mean(0); cnt[nb[f]] += 1
    cc /= cnt[:, None]

    r = records(a.restart)
    arr = [np.frombuffer(x, dtype='<f8') for x in r[1:]]
    ntot = nc + nB
    big = [x for x in arr if len(x) == ntot]
    u, v, w, p, te, ed, T, vis = big[:8]
    mu = info['mu']
    if a.pinf is not None:
        info['pinf'] = a.pinf
    q = 0.5 * info['rhoinf'] * info['Uinf'] ** 2
    span = info['thickness_z']

    F = np.zeros(3); Fv = np.zeros(3)
    rows = []
    for k in range(nw):
        f = sw + k
        nodes = faces[f]
        pts = P[nodes]; c = pts.mean(0); S = np.zeros(3)
        for i in range(len(nodes)):
            S += np.cross(pts[i] - c, pts[(i + 1) % len(nodes)] - c) / 2
        ib = nc + (f - nI)
        ip = own[f]
        area = np.linalg.norm(S); n = S / area
        dn = np.dot(c - cc[ip], n)
        Up = np.array([u[ip], v[ip], w[ip]])
        Ut = Up - np.dot(Up, n) * n
        Fp = (p[ib] - info['pinf']) * S
        Ff = mu * area * Ut / dn
        F += Fp + Ff; Fv += Ff
        tvec = Ut / (np.linalg.norm(Ut) + 1e-300)
        cf = mu * np.linalg.norm(Ut) / dn / q
        rows.append((c[0], c[1], (p[ib] - info['pinf']) / q, cf, np.dot(Ut, [1, 0, 0]), np.dot(Fp, [1, 0, 0]),
                     np.dot(Ff, [1, 0, 0])))
    al = math.radians(info['alpha'])
    Dir = np.array([math.cos(al), math.sin(al), 0]); Lir = np.array([-math.sin(al), math.cos(al), 0])
    cd = np.dot(F, Dir) / (q * span); cl = np.dot(F, Lir) / (q * span)
    cdv = np.dot(Fv, Dir) / (q * span)
    print('Cd = %.5f (pressure %.5f, viscous %.5f)   Cl = %.5f' % (cd, cd - cdv, cdv, cl))
    R = np.array(rows)
    up = R[R[:, 1] > 0]; up = up[np.argsort(-up[:, 0])]           # upper surface, from the TE to the LE
    # separation point on the upper surface: tangential velocity (x-component) changes sign going from the LE downstream
    ux = up[::-1, 4]; xs = up[::-1, 0]
    sep = None
    for i in range(1, len(ux)):
        if xs[i] > 0.3 and ux[i - 1] > 0 and ux[i] <= 0:
            sep = xs[i - 1] + (xs[i] - xs[i - 1]) * ux[i - 1] / (ux[i - 1] - ux[i]); break
    print('x_sep/c (upper surface) =', 'none' if sep is None else '%.4f' % sep)
    print('Tmin/Tmax = %.2f / %.2f K   pmin/pmax = %.0f / %.0f Pa' % (T[:nc].min(), T[:nc].max(), p[:nc].min(), p[:nc].max()))
    c = np.sqrt(1.39926 * 287.058 * T[:nc]); M = np.sqrt(u[:nc]**2 + v[:nc]**2 + w[:nc]**2) / c
    print('max Mach = %.4f' % M.max())
    if a.cp:
        np.savetxt(a.cp, R[np.argsort(R[:, 0])][:, :4], header='x y Cp Cf')


if __name__ == '__main__':
    main()
