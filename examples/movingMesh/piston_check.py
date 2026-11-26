#!/usr/bin/env python3
"""
piston_check.py - compare a finished 'piston' run (mkcase.py piston) with the exact solution.

  python3 piston_check.py <case dir> [restart file name]

Free slip walls, piston at x = 0 moving with the displacement d(t) (translate: U t, oscillate: A sin(2 pi f t)), fixed
pressure p = 0 at the outlet x = L. All the fluid moves as a plug, every time step the velocity is the average velocity
of the piston in the step (continuity on the moving mesh) and the pressure gradient is that of the acceleration
(backward Euler of the solver):

    u^n = ( d(t_n) - d(t_n-1) ) / dt,     p^n(x) = rho (u^n - u^(n-1))/dt * (L - x)         x = lab. coordinate of the cell centre.

The cell centres are taken from the points of the VTK file of the last step (written with the current coordinates), as the mean
of the 8 nodes of the (hexahedral) cell. Prints the errors, exit code 1 if they are larger than the tolerances.
"""
import os, sys, struct, re, math
import numpy as np


def read_restart(fn):
    b = open(fn, 'rb').read(); o = 0; recs = []
    while o < len(b):
        n = struct.unpack('<i', b[o:o + 4])[0]; recs.append(b[o + 4:o + 4 + n]); o += 8 + n
    it, t = struct.unpack('<id', recs[0][:12])
    d = [np.frombuffer(x, dtype='<f8') for x in recs[1:]]
    return it, t, d       # d: flmass, u, v, w, p, ...


def read_vtu_points_cells(fn):
    txt = open(fn).read()
    def block(name):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array(m.group(1).split(), dtype=float)
    pts = block('Points').reshape(-1, 3)
    conn = block('connectivity').astype(int).reshape(-1, 8)
    return pts, conn


def main():
    case = sys.argv[1]
    rfile = sys.argv[2] if len(sys.argv) > 2 else 'restart'
    cs = np.load(os.path.join(case, 'case.npz'))
    L = float(cs['L']); dt = float(cs['dt']); motion = str(cs['motion']); U = float(cs['U'])
    A = float(cs['amp']); f = float(cs['freq']); nd = int(cs['ndim'])
    it, t, d = read_restart(os.path.join(case, rfile))
    fl, u, v, w, p = d[0], d[1], d[2], d[3], d[4]
    nc = int(cs['nx']) * int(cs['ny']) * (1 if nd == 2 else int(cs['ny']))
    u, v, w, p = u[:nc], v[:nc], w[:nc], p[:nc]
    disp = (lambda tt: U * tt) if motion == 'translate' else (lambda tt: A * math.sin(2 * math.pi * f * tt))
    un = (disp(t) - disp(t - dt)) / dt
    uo = (disp(t - dt) - disp(t - 2 * dt)) / dt if it > 1 else un
    acc = (un - uo) / dt
    if abs(acc) < 1e-9:
        acc = 0.0           # constant velocity (roundoff of the difference)

    vtu = os.path.join(case, 'vtk', str(it), 'interior.vtu')
    pts, conn = read_vtu_points_cells(vtu)
    xc = pts[conn][:, :, 0].mean(axis=1)
    px = pts[:, 0]
    piston_ok = abs(px.min() - disp(t))

    pex = acc * (L - xc)
    # The zero-gradient pressure condition on the accelerating piston wall (dp/dn should be -rho*a) gives a first order
    # (O(dx)) error in the cells next to the piston, so with acceleration the comparison is made in the bulk, further
    # than 'bulk' * L from the piston. Without acceleration (translate) all cells are compared.
    bulk = 0.0 if acc == 0.0 else 0.25
    m = (xc - disp(t)) > bulk * L
    eu = np.abs(u - un)[m].max()
    ev = max(np.abs(v[m]).max(), np.abs(w[m]).max())
    ep = np.abs(p - pex)[m].max()
    pscale = max(np.abs(pex).max(), 1e-30)
    print('step %d  t = %.5f   piston position error %.2e (displacement %.5f)' % (it, t, piston_ok, disp(t)))
    print('u: exact %.8f   max|u - u_exact| = %.3e   max|v|,|w| = %.3e   (cells x > %.2f L)' % (un, eu, ev, bulk))
    print('p: max|p_exact| = %.3e   max|p - p_exact| = %.3e  (%.2e relative)' % (pscale, ep, ep / pscale if pscale > 1e-20 else 0.0))
    tol_u = 1e-7 * max(abs(un), 1.0) if bulk == 0.0 else 5e-3 * max(abs(un), 1.0)
    tol_v = 1e-8 if bulk == 0.0 else 5e-3
    ok = eu < tol_u and ev < tol_v and (ep < 2e-2 * pscale if pscale > 1e-20 else ep < 1e-8) and piston_ok < 1e-12
    print('RESULT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
