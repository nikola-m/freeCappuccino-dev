#!/usr/bin/env python3
"""
mkairfoil.py - C-grid (default, --grid c) or O-grid (--grid o) around a NACA 00xx airfoil (2D, one cell thick, 'empty' front/back) with a case set-up for
the all-speeds block-coupled p-U solver (compressible = T, coupled = T), in the native freeCappuccino mesh format.

  python3 mkairfoil.py <dir> [--grid c|o --nt 201 --nw 120 --nr 96 --R 20 --y1 1e-3 --M 0.5 --Re 5000 --alpha 0 --Pr 0.72 --te 3e-4 --cfl 200 --ramp 0 --sutherland]

C-grid : sharp trailing edge (shared node), nt cells on the surface, nw cells in each half of the wake (wake cut at y = 0 is
         an inner interface), nr cells radially; grid lines leave the wall along its normal and turn vertical towards the TE,
         so there is no fan of rays behind the trailing edge (the O-grid weakness).  Far field: x_out = 0.5 + R behind
         the profile, circle of radius R in front.  Patches: far-in (front half circle, 'inlet'), far-out (top/bottom aft part
         and the outlet, 'pressure').  Mesh quality: python3 meshquality.py <dir>.
O-grid : nt cells around (odd), nr cells in the radial direction, geometric growth chosen to reach the far field at
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


def wake_distribution(d0, length, nw):
    """nw+1 positions from 0 to `length` with first spacing d0 and constant growth ratio g (g >= 1)."""
    if d0 * nw >= length:                              # uniform is already too coarse a start: fall back to uniform
        return np.linspace(0.0, length, nw + 1), 1.0
    f = lambda g: d0 * (g**nw - 1) / (g - 1) - length
    lo, hi = 1.0 + 1e-9, 4.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        (lo, hi) = (mid, hi) if f(mid) < 0 else (lo, mid)
    g = 0.5 * (lo + hi)
    return np.concatenate([[0.0], np.cumsum(d0 * g ** np.arange(nw))]), g


def airfoil_nodes_c(N, t=0.12, lam=0.1):
    """Sharp NACA 00xx surface, N+1 nodes lower TE -> LE -> upper TE (N even). Cosine spacing blended with the power law
    x = (1 - tau)^2 by the weight lam, which keeps the leading-edge spacing but makes the trailing-edge spacing coarser
    (about the wall-normal first cell height instead of 4x finer)."""
    h = N // 2
    tau = np.arange(h + 1) / h                      # 0 at the TE, 1 at the LE
    xl = (1 - lam) * 0.5 * (1 + np.cos(math.pi * tau)) + lam * (1 - tau) ** 2
    x = np.concatenate([xl, xl[::-1][1:]])
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    y = np.concatenate([-yt[:h + 1], yt[h + 1:]])
    y[0] = y[-1] = 0.0
    return np.column_stack([x, y])


def build_c(nt, nw, nr, R, y1, wake_h1=0.0, te_lam=0.1):
    """C-grid. Node line j = 0 (I = 0..M, M = N + 2 nw, N = nt - 1 cells on the airfoil): lower wake (outlet -> TE),
    lower surface (TE -> LE), upper surface (LE -> TE), upper wake (TE -> outlet). The two wake halves lie on the
    cut y = 0 and share their nodes. Sharp trailing edge (NACA 0012 with a4 = 0.1036).
    Far field: bottom line y = -R from the outlet to x = 0.5, half circle of radius R about (0.5, 0), top line, outlet x = 0.5 + R.
    Returns X (nr+1, M+1, 2), N, growth ratios (radial, wake)."""
    N = nt - 1
    assert N % 2 == 0, 'nt must be odd'
    P = airfoil_nodes_c(N, lam=te_lam)
    P[0] = P[N] = (1.0, 0.0)                              # exactly sharp
    xout = 0.5 + R
    d0 = np.linalg.norm(P[1] - P[0])
    sw, gw = wake_distribution(d0, xout - 1.0, nw)        # distance from the TE
    xw = 1.0 + sw                                         # xw[q], q = 0 at the TE ... nw at the outlet
    M = N + 2 * nw
    # wall normals (outward); at the trailing edge the grid lines are vertical so that they stay inside the wake columns
    tg = np.empty_like(P)
    tg[1:-1] = P[2:] - P[:-2]; tg[0] = P[1] - P[0]; tg[-1] = P[-1] - P[-2]
    tg /= np.linalg.norm(tg, axis=1)[:, None]
    nrm = np.column_stack([tg[:, 1], -tg[:, 0]])
    ile = int(np.argmin(P[:, 0]))
    if np.dot(nrm[ile], P[ile] - np.array([0.5, 0.0])) < 0:
        nrm = -nrm
    # far-field points. Aft of mid-chord the lines are vertical (they continue as the wake columns and end on the
    # bottom / top line at the same x); in front they are rays from (0.5, 0) with the angle blended 50 % with a uniform
    # distribution, as in the O-grid. The wall-normal direction turns to vertical towards the trailing edge, so that
    # the grid lines of the last surface cells cannot cross the wake columns.
    c0 = np.array([0.5, 0.0])
    P0 = np.zeros((M + 1, 2)); Pf = np.zeros((M + 1, 2)); nn = np.zeros((M + 1, 2))
    front = [m for m in range(N + 1) if P[m, 0] < 0.5]                    # nodes in front of mid-chord (contiguous, around the LE)
    mlo, mhi = front[0], front[-1]                                        # lower side starts at mlo, upper side ends at mhi
    # angular position (from the vertical through mid-chord) of the front far-field points: geometric growth that starts
    # with the x-spacing of the aft columns at the junction and ends at 90 degrees at the leading edge
    dxj = abs(P[mlo - 1, 0] - P[mlo, 0]) / R
    # far-field x positions of the aft surface nodes and of the wake columns: geometric growth starting with the junction
    # spacing, so that the dense columns at the trailing edge fan out slightly (a few degrees over R chords) instead of
    # staying 1e-4 apart up to the far field. Y[k], k = 0 at mid-chord ... naft - 1 at the TE (m = 0) ... naft - 1 + nw at the outlet
    naft = mlo                                            # lower-side nodes m = 0 .. mlo-1 have x >= 0.5
    Xaft = [P[m, 0] for m in range(naft - 1, -1, -1)]     # increasing x, last = TE
    Ys, _ = wake_distribution(dxj * R, xout - Xaft[0], naft - 1 + nw)
    Yfar = Xaft[0] + Ys                                   # same x for the lower and the upper side
    far_aft = lambda m: Yfar[naft - 1 - m]                # m = lower-side index (upper side: N - m)
    far_wake = lambda q: Yfar[naft - 1 + q]               # q = 0 at the TE ... nw at the outlet
    ang_lo, _ = wake_distribution(dxj, 0.5 * math.pi, ile - mlo + 1)
    ang_up, _ = wake_distribution(dxj, 0.5 * math.pi, mhi - ile + 1)
    for I in range(M + 1):
        if I < nw:                                        # lower wake, I = 0 at the outlet
            q = nw - I
            P0[I] = (xw[q], 0.0); Pf[I] = (far_wake(q), -R); nn[I] = (0.0, -1.0)
        elif I <= nw + N:                                 # surface (I = nw: lower TE, I = nw + N: upper TE)
            m = I - nw
            P0[I] = P[m]
            lower = m <= ile
            sg = -1.0 if lower else 1.0
            if P[m, 0] >= 0.5:
                Pf[I] = (far_aft(m if lower else N - m), sg * R)
            else:
                ang = (-0.5 * math.pi - ang_lo[m - mlo + 1]) if lower else (0.5 * math.pi + ang_up[mhi - m + 1])
                Pf[I] = c0 + R * np.array([math.cos(ang), math.sin(ang)])
            u = min(max((P[m, 0] - 0.5) / 0.5, 0.0), 1.0)
            beta = u * u * (3 - 2 * u)
            d = (1 - beta) * nrm[m] + beta * np.array([0.0, sg])
            nn[I] = d / np.linalg.norm(d)
        else:                                             # upper wake
            q = I - (nw + N)
            P0[I] = (xw[q], 0.0); Pf[I] = (far_wake(q), R); nn[I] = (0.0, 1.0)
    sj, g = radial_distribution(nr, y1 / R)
    X = np.zeros((nr + 1, M + 1, 2))
    for I in range(M + 1):
        ray = Pf[I] - P0[I]
        L = np.linalg.norm(ray); e = ray / L
        h1 = y1
        if I < nw or I > nw + N:                          # wake: first cell height grows with the distance behind the TE
            h1 = max(y1, wake_h1 * (P0[I, 0] - 1.0))
        sI = sj if h1 == y1 else radial_distribution(nr, h1 / R)[0]
        for j in range(nr + 1):
            t = sI[j] * L
            w = 1.0 - math.exp(-t ** 2)
            d = (1 - w) * nn[I] + w * e
            d /= np.linalg.norm(d)
            X[j, I] = P0[I] + t * d
        X[nr, I] = Pf[I]
    return X, N, g, gw


def hex_mesh_c(X, N, thick=0.1, alpha=0.0):
    nr1, M1, _ = X.shape
    nr, M = nr1 - 1, M1 - 1
    nw = (M - N) // 2
    canon = lambda I, j: (M - I, 0) if (j == 0 and I >= nw + N) else (I, j)
    ids2d = {}
    pts2d = []
    for j in range(nr1):
        for I in range(M1):
            key = canon(I, j)
            if key not in ids2d:
                ids2d[key] = len(pts2d)
                pts2d.append(X[key[1], key[0]])
    n2 = len(pts2d)
    pts = [(p[0], p[1], thick * k) for k in range(2) for p in pts2d]
    P = np.array(pts)
    nid = lambda I, j, k: 1 + ids2d[canon(I, j)] + n2 * k
    cid = lambda i, j: 1 + i + M * j
    cen = {}
    for j in range(nr):
        for i in range(M):
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
    ca, sa = math.cos(alpha), math.sin(alpha)
    for j in range(nr):
        for i in range(M):
            c = cid(i, j)
            # streamwise faces at node line i (between cells i-1 and i)
            nodes = [nid(i, j, 0), nid(i, j + 1, 0), nid(i, j + 1, 1), nid(i, j, 1)]
            if i == 0:
                bnd['far-out'].append((c, orient(nodes, cen[c])))                       # outlet, lower half
            else:
                nb = cid(i - 1, j)
                inner.append((nb, c, orient(nodes, cen[nb])))
            if i == M - 1:
                nodes = [nid(M, j, 0), nid(M, j + 1, 0), nid(M, j + 1, 1), nid(M, j, 1)]
                bnd['far-out'].append((c, orient(nodes, cen[c])))                       # outlet, upper half
            # faces on the node line j + 1 (between cells (i,j) and (i,j+1)) or the far field
            nodes = [nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]
            if j < nr - 1:
                inner.append((c, cid(i, j + 1), orient(nodes, cen[c])))
            else:
                fc = P[np.array(nodes) - 1].mean(0)
                key = 'far-in' if (fc[0] - 0.5) * ca + fc[1] * sa < 0 else 'far-out'
                bnd[key].append((c, orient(nodes, cen[c])))
            if j == 0:
                nodes = [nid(i, 0, 0), nid(i + 1, 0, 0), nid(i + 1, 0, 1), nid(i, 0, 1)]
                if nw <= i < nw + N:
                    bnd['wall'].append((c, orient(nodes, cen[c])))
                elif i < nw:                                                           # wake cut: lower cell i | upper cell M-1-i
                    nb = cid(M - 1 - i, 0)
                    inner.append((c, nb, orient(nodes, cen[c])))
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
        amg_smoother = 'gs'
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
    ap.add_argument('--te', type=float, default=3e-4, help='half thickness of the trailing-edge base (O-grid only; the C-grid is sharp)')
    ap.add_argument('--grid', default='c', choices=['c', 'o'], help='C-grid with a wake cut (default, sharp trailing edge) or the O-grid')
    ap.add_argument('--tespace', type=float, default=0.1, help='C-grid: weight of the power-law distribution of the surface nodes (0 = cosine, coarser trailing-edge spacing for larger values)')
    ap.add_argument('--nw', type=int, default=120, help='C-grid: number of cells in the wake (each half) from the trailing edge to the outlet')
    ap.add_argument('--urfU', type=float, default=0.7); ap.add_argument('--urfEn', type=float, default=0.8)
    ap.add_argument('--blend', type=float, default=1.0)
    ap.add_argument('--cfl', type=float, default=None, help='pseudo-time CFL (O-grid: 200 for M<=0.6, 30 above; C-grid: 10 / 5)')
    ap.add_argument('--ramp', type=int, default=None, help='iterations of the CFL ramp (O-grid: 0, 150 for M>0.6; C-grid: 200)')
    ap.add_argument('--sutherland', action='store_true', help='Sutherland viscosity law')
    ap.add_argument('--numstep', type=int, default=3000); ap.add_argument('--maxit', type=int, default=40)
    ap.add_argument('--farfield', default='default', choices=['default', 'riemann', 'riemann-in', 'riemann-out', 'open', 'compoutlet'],
                    help='boundary conditions of the far field (implicit_bc_compressible.f90): default = velocity inlet '
                         '(far-in) + static pressure outlet (far-out); riemann = both patches (-in / -out: one of them); '
                         'open, compoutlet = far-out')
    a = ap.parse_args()
    os.makedirs(os.path.join(a.dir, 'vtk', 'mesh'), exist_ok=True)

    al = math.radians(a.alpha)
    if a.grid == 'c':    # the C-grid resolves the near wake; the pseudo-time march needs a moderate CFL (200 diverges)
        cfl = a.cfl if a.cfl is not None else (10.0 if a.M <= 0.6 else 5.0)
        ramp = a.ramp if a.ramp is not None else 200
    else:
        cfl = a.cfl if a.cfl is not None else (200.0 if a.M <= 0.6 else 30.0)
        ramp = a.ramp if a.ramp is not None else (0 if a.M <= 0.6 else 150)
    suth = 'T' if a.sutherland else 'F'
    if a.grid == 'c':
        X, N, g, gw = build_c(a.nt, a.nw, a.nr, a.R, a.y1, te_lam=a.tespace)
        mesh = hex_mesh_c(X, N, 0.1, al)
        gridinfo = 'C-grid, wake growth ratio %.4f, ' % gw
    else:
        X, g = build(a.nt, a.nr, a.R, a.y1, a.te, al)
        mesh = hex_mesh(X, 0.1, al)
        gridinfo = 'O-grid, '
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
    print(a.dir, info, gridinfo + 'radial growth ratio %.4f' % g, 'U_inf %.3f m/s, mu %.5g Pa s, rho %.4f' % (U, mu, rho))


if __name__ == '__main__':
    main()
