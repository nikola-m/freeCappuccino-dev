#!/usr/bin/env python3
"""
clsvof_post.py - post-processing of CLSVOF runs (reads clsvof.restart, clsvof_history.dat, cells.npz of a case directory).

  python3 clsvof_post.py shape <dir> --ic cylinder|sphere --center x y z --radius R [--inside 1|2]
        L1 error of the volume fraction against the exact shape  E = sum |vof - vof_exact| V   and the level set error
        (in the cells with |lsf| < 3h) against the exact signed distance. Used for the advection tests after the flow returned
        to the initial state (reversing vortex, one revolution) and for the static drop.
  python3 clsvof_post.py history <dir>
        prints the monitored quantities of clsvof_history.dat: the volume loss, bubble centroid, rise velocity, max velocity.
"""
import sys, os, re, argparse
import numpy as np


def read_vtu_field(fn, name, ncomp=1, cell=True):
    """Cell data array 'name' of an ascii interior.vtu written by the solver (the last array of that name)."""
    txt = open(fn).read()
    m = re.findall(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, flags=re.S)
    return np.array(m[-1].split(), dtype=float).reshape(-1, ncomp).squeeze()


def read_restart(d):
    c = np.load(os.path.join(d, 'cells.npz'))
    n = len(c['xc'])
    with open(os.path.join(d, 'clsvof.restart'), 'rb') as f:
        hdr = np.fromfile(f, dtype='<f8', count=4)
        vof = np.fromfile(f, dtype='<f8', count=n)
        lsf = np.fromfile(f, dtype='<f8', count=n)
    return c, hdr, vof, lsf


def exact_vof(c, ic, center, R, inside, nsub=8):
    """Exact volume fraction of phase 1 by sub-sampling each cell with nsub^ndim points."""
    nx, ny, nz = int(c['nx']), int(c['ny']), int(c['nz'])
    hx, hy = float(c['Lx']) / nx, float(c['Ly']) / ny
    hz = float(c['Lz']) / nz if nz > 0 else float(c['Lz'])
    s = (np.arange(nsub) + 0.5) / nsub - 0.5
    ox, oy = np.meshgrid(s * hx, s * hy, indexing='ij')
    ox, oy = ox.ravel(), oy.ravel()
    if ic == 'sphere':
        oz = None
        sz = s * hz
        X, Y, Z = np.meshgrid(s * hx, s * hy, sz, indexing='ij')
        ox, oy, oz = X.ravel(), Y.ravel(), Z.ravel()
    f = np.zeros(len(c['xc']))
    for i in range(len(f)):
        dx = c['xc'][i] + ox - center[0]
        dy = c['yc'][i] + oy - center[1]
        r2 = dx * dx + dy * dy
        if ic == 'sphere':
            r2 = r2 + (c['zc'][i] + oz - center[2]) ** 2
        ins = (r2 < R * R).mean()
        f[i] = ins if inside == 1 else 1.0 - ins
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['shape', 'history', 'laplace', 'hysing', 'fall'])
    ap.add_argument('dir')
    ap.add_argument('--ic', default='cylinder')
    ap.add_argument('--center', type=float, nargs=3, default=[0.5, 0.75, 0.0])
    ap.add_argument('--radius', type=float, default=0.15)
    ap.add_argument('--inside', type=int, default=1)
    ap.add_argument('--sigma', type=float, default=1.0)
    ap.add_argument('--gravity', type=float, default=1.0)
    ap.add_argument('--density-ratio', type=float, default=1000.0)
    ap.add_argument('--added-mass', type=float, default=1.0, help='added mass coefficient: 1 cylinder, 0.5 sphere')
    a = ap.parse_args()
    if a.mode == 'shape':
        c, hdr, vof, lsf = read_restart(a.dir)
        fe = exact_vof(c, a.ic, a.center, a.radius, a.inside)
        V = c['vol']
        E = np.sum(np.abs(vof - fe) * V)
        V0 = np.sum(fe * V)
        print('time = %g' % hdr[0])
        print('L1 error of vof   E = sum|f - f_exact| V = %.6e   (relative to initial volume: %.4e)' % (E, E / V0))
        print('volume of phase 1: %.9e  exact initial %.9e  relative change %.3e' % (np.sum(vof * V), V0, np.sum(vof * V) / V0 - 1))
        print('clipped volume (removed - added): %.3e' % hdr[3])
        # level set error in the band
        dx = c['xc'] - a.center[0]; dy = c['yc'] - a.center[1]
        r = dx * dx + dy * dy
        if a.ic == 'sphere':
            r = r + (c['zc'] - a.center[2]) ** 2
        r = np.sqrt(r)
        phi = (a.radius - r) * (1 if a.inside == 1 else -1)
        h = float(c['Lx']) / int(c['nx'])
        m = np.abs(phi) < 3 * h
        print('level set: max |lsf - phi| in |phi| < 3h : %.3e (h = %g), L1 : %.3e' %
              (np.max(np.abs(lsf - phi)[m]), h, np.mean(np.abs(lsf - phi)[m])))
    elif a.mode == 'laplace':
        # pressure jump of a static drop: p_in - p_out = sigma/R (2D cylinder), 2 sigma/R (3D sphere)
        c = np.load(os.path.join(a.dir, 'cells.npz'))
        step = sorted(int(x) for x in os.listdir(os.path.join(a.dir, 'vtk')) if x.isdigit())[-1]
        fn = os.path.join(a.dir, 'vtk', str(step), 'interior.vtu')
        p = read_vtu_field(fn, 'p'); f = read_vtu_field(fn, 'vof'); U = read_vtu_field(fn, 'U', 3)
        V = c['vol']
        rin = (f > 0.999); rout = (f < 0.001)
        # mean pressure in the drop away from the interface and outside (volume weighted)
        pin = np.sum(p[rin] * V[rin]) / np.sum(V[rin]); pout = np.sum(p[rout] * V[rout]) / np.sum(V[rout])
        R = a.radius; sigma = a.sigma
        ex = sigma / R if a.ic == 'cylinder' else 2 * sigma / R
        print('step %d: p_in = %.5f  p_out = %.5f  jump = %.5f   exact = %.5f   error = %.3f %%' %
              (step, pin, pout, pin - pout, ex, 100 * ((pin - pout) / ex - 1)))
        print('max |U| = %.3e   (sigma = %g)' % (np.max(np.linalg.norm(U, axis=1)), sigma))
    elif a.mode == 'fall':
        # free fall of a heavy drop in a light gas: initial acceleration  g (1 - r)/(1 + Cm r),  r = rho_gas/rho_drop
        h = np.loadtxt(os.path.join(a.dir, 'clsvof_history.dat'))
        t = h[:, 0]; vy = h[:, 11]
        m = t < 0.1 * t[-1] if t[-1] < 1.0 else t < 0.1
        acc = -np.polyfit(t[m], vy[m], 1)[0]
        r = 1.0 / a.density_ratio
        th = a.gravity * (1 - r) / (1 + a.added_mass * r)
        print('initial acceleration: measured %.5f  theory (added mass) %.5f  (%.2f %%)' % (acc, th, 100 * (acc / th - 1)))
        print('t_end = %g  drop velocity %.4f  drop centroid y = %.4f  volume change %.2e' % (t[-1], vy[-1], h[-1, 5], h[-1, 2]))
    elif a.mode == 'hysing':
        # Hysing et al. (2009) benchmark quantities of the 2D bubble: centre of mass, rise velocity, circularity
        c = np.load(os.path.join(a.dir, 'cells.npz'))
        dz = float(c['Lz']); Ltot = float(c['Lx']) * float(c['Ly']) * dz
        h = np.loadtxt(os.path.join(a.dir, 'clsvof_history.dat'))
        t = h[:, 0]; yc = h[:, 8]; vy = h[:, 14]
        A = (Ltot - h[:, 1]) / dz                         # bubble area
        P = h[:, 17] / dz                                 # interface length (PLIC)
        circ = np.pi * np.sqrt(4 * A / np.pi) / P         # perimeter of the area-equivalent circle / perimeter
        i1 = np.argmax(vy); i2 = np.argmin(circ)
        print('final centre of mass  y_c(t=%.3f) = %.4f     [ref 1.0813 / 1.0799 / 1.0817]' % (t[-1], yc[-1]))
        print('max rise velocity     V_c,max = %.4f at t = %.4f   [ref 0.2417 / 0.2421, t = 0.9213 / 0.9313]' % (vy[i1], t[i1]))
        print('min circularity       = %.4f at t = %.4f   [ref 0.9013 / 0.9011, t = 1.9041 / 1.8750]' % (circ[i2], t[i2]))
        print('final circularity %.4f ; bubble area change %.3e' % (circ[-1], A[-1] / A[0] - 1))
    else:
        h = np.loadtxt(os.path.join(a.dir, 'clsvof_history.dat'))
        t = h[:, 0]; yc2 = h[:, 9]; vol = h[:, 2]
        print('t_end = %g  volume change = %.3e  clip = %.3e' % (t[-1], vol[-1], h[-1, 3]))
        print('min / max vertical centroid of phase 2: %.6f %.6f' % (yc2.min(), yc2.max()))
        print('final centroid of phase 2 = (%.6f %.6f %.6f)' % tuple(h[-1, 7:10]))


if __name__ == '__main__':
    main()
