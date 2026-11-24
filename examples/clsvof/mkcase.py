#!/usr/bin/env python3
"""
mkcase.py - generate CLSVOF test cases in the native freeCappuccino 'nativeMesh' format
(polyMesh/, 0/U, vtk/mesh templates, input.nml and cells.npz with the cell centres and volumes for the post-processing).

  python3 mkcase.py advect-vortex  <dir> [--n 64 --T 8 --steps 400]       reversing vortex in a box (LeVeque), 2D
  python3 mkcase.py advect-rotate  <dir> [--n 64 --steps 628]              solid body rotation of a disk, one revolution, 2D
  python3 mkcase.py advect-sphere  <dir> [--n 32 --T 2 --steps 200]        reversing 3D vortex (sphere in a box)  (use --ndim 3)
  python3 mkcase.py static-drop    <dir> [--n 32 --ndim 2|3 --steps 200]   Laplace pressure jump, no gravity
  python3 mkcase.py rising-bubble  <dir> [--n 40 --case 1|2 --steps ...]   Hysing et al. (2009) 2D benchmark, h = 1/n
  python3 mkcase.py falling-drop   <dir> [--n 32 --ndim 3]                 water drop falling in air in a closed box

2D cases are one cell thick with 'empty' front/back patches (as in examples/coupledPU), 3D cases are boxes of nx x ny x nz cells.
"""
import argparse, os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'coupledPU'))
from mkcase import build_mesh, write_native, write_U          # noqa: E402  (reused 2D mesh writer)


def box3d(nx, ny, nz, Lx, Ly, Lz):
    """Structured hexahedral box, one 'walls' patch. Cell index: 1 + i + nx*(j + ny*k)."""
    X = np.linspace(0, Lx, nx + 1); Y = np.linspace(0, Ly, ny + 1); Z = np.linspace(0, Lz, nz + 1)
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    cid = lambda i, j, k: 1 + i + nx * (j + ny * k)
    pts = [(X[i], Y[j], Z[k]) for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    fx = lambda i, j, k: [nid(i, j, k), nid(i, j + 1, k), nid(i, j + 1, k + 1), nid(i, j, k + 1)]
    fy = lambda i, j, k: [nid(i, j, k), nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j, k)]
    fz = lambda i, j, k: [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]
    inner, walls = [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = cid(i, j, k)
                (inner.append((c, cid(i + 1, j, k), fx(i + 1, j, k))) if i < nx - 1 else walls.append((c, fx(i + 1, j, k))))
                (inner.append((c, cid(i, j + 1, k), fy(i, j + 1, k))) if j < ny - 1 else walls.append((c, fy(i, j + 1, k))))
                (inner.append((c, cid(i, j, k + 1), fz(i, j, k + 1))) if k < nz - 1 else walls.append((c, fz(i, j, k + 1))))
                if i == 0: walls.append((c, fx(i, j, k)[::-1]))
                if j == 0: walls.append((c, fy(i, j, k)[::-1]))
                if k == 0: walls.append((c, fz(i, j, k)[::-1]))
    inner.sort(key=lambda t: (t[0], t[1]))
    return dict(pts=pts, inner=inner, bnd={'walls': walls})


def make_mesh(root, nx, ny, nz, Lx, Ly, Lz, side='wall'):
    os.makedirs(os.path.join(root, 'vtk', 'mesh'), exist_ok=True)
    if nz == 0:   # 2D: one cell thick slab, dz = Lz
        xn = np.linspace(0, Lx, nx + 1); yn = np.linspace(0, Ly, ny + 1)
        mesh = build_mesh(xn, yn)
        mesh['pts'] = [(p[0], p[1], p[2] / 0.1 * Lz) for p in mesh['pts']]
        patches = [('walls-l', side, 'left'), ('walls-r', side, 'right'), ('walls-b', 'wall', 'bottom'),
                   ('walls-t', 'wall', 'top'), ('front-back', 'empty', 'frontback')]
        info = write_native(root, mesh, patches)
        write_U(root, (0, 0, 0), [(n, (0, 0, 0)) if t == 'wall' else (n, None) for n, t, _ in patches])  # symmetry/empty: zero gradient
        i, j = np.meshgrid(np.arange(nx), np.arange(ny))
        xc = (i.ravel() + 0.5) * Lx / nx; yc = (j.ravel() + 0.5) * Ly / ny; zc = np.full(xc.shape, 0.5 * Lz)
        vol = np.full(xc.shape, Lx / nx * Ly / ny * Lz)
    else:
        mesh = box3d(nx, ny, nz, Lx, Ly, Lz)
        info = write_native(root, mesh, [('walls', 'wall', 'walls')])
        write_U(root, (0, 0, 0), [('walls', (0, 0, 0))])
        i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
        order = np.lexsort((i.ravel(), j.ravel(), k.ravel()))
        xc = ((i.ravel() + 0.5) * Lx / nx)[order]; yc = ((j.ravel() + 0.5) * Ly / ny)[order]
        zc = ((k.ravel() + 0.5) * Lz / nz)[order]
        vol = np.full(xc.shape, Lx / nx * Ly / ny * Lz / nz)
    np.savez(os.path.join(root, 'cells.npz'), xc=xc, yc=yc, zc=zc, vol=vol, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)
    print(root, info)
    return info


INPUT = """&input_parameters
              title = '{title}'
        mesh_format = 'nativeMesh'
              lread = F
              lwrite = T
              ltest = F
              calcU = {calcU}
               urfU = 0.9, 0.9, 0.9
               gdsU = 1.0
           cSchemeU = 'linearUpwind'
           lSolverU = 'bicgstab'
           maxiterU = 10
            tolAbsU = 1e-14
            tolRelU = 0.001
            pscheme = 'linear'
              calcP = {calcU}
               urfP = {urfP}
           lSolverP = 'iccg'
           maxiterP = 400
            tolAbsP = 1e-14
            tolRelP = 0.001
          TurbModel = 'none'
             densit = {densit}
             viscos = {viscos}
         ltransient = T
                bdf = T
              gauss = T
            limiter = 'no-limit'
             SIMPLE = T
              ncorr = 1
              npcor = 1
           pRefCell = 1
          tolerance = 1e-12
            NUMSTEP = {steps}
           TIMESTEP = {dt}
             NZAPIS = {nzapis}
              MAXIT = {maxit}
               gravx = 0.0, gravy = {gy}, gravz = {gz}
{extra}/
"""


def write_input(root, **kw):
    kw.setdefault('calcU', 'T'); kw.setdefault('urfP', 0.5); kw.setdefault('densit', 1.0); kw.setdefault('viscos', 1e-3)
    kw.setdefault('maxit', 4); kw.setdefault('gy', 0.0); kw.setdefault('gz', 0.0); kw.setdefault('extra', '')
    kw.setdefault('nzapis', kw['steps'])
    with open(os.path.join(root, 'input.nml'), 'w') as f:
        f.write(INPUT.format(**kw))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case', choices=['advect-vortex', 'advect-rotate', 'advect-sphere', 'static-drop', 'rising-bubble',
                                     'falling-drop'])
    ap.add_argument('dir')
    ap.add_argument('--n', type=int, default=None, help='cells per unit length (or per domain side)')
    ap.add_argument('--ndim', type=int, default=None)
    ap.add_argument('--T', type=float, default=None)
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--bcase', type=int, default=1, help='Hysing case 1 (sigma=24.5) or 2 (sigma=1.96)')
    ap.add_argument('--maxit', type=int, default=None)
    a = ap.parse_args()
    d = a.dir
    if a.case == 'advect-vortex':
        n = a.n or 64; T = a.T or 8.0; steps = a.steps or 400
        make_mesh(d, n, n, 0, 1.0, 1.0, 1.0 / n)
        ex = ("           twoPhase = T\n  cl_flow = 'vortex', cl_T = %g\n  cl_ic = 'cylinder', cl_inside = 1\n"
              "  cl_center = 0.5, 0.75, 0.0, cl_radius = 0.15\n  cl_rho1 = 1.0, cl_rho2 = 1.0, cl_mu1 = 1e-3, cl_mu2 = 1e-3\n"
              "  cl_cfl = 0.15\n" % T)
        write_input(d, title='CLSVOF reversing vortex', steps=steps, dt=T / steps, calcU='F', maxit=1, extra=ex)
    elif a.case == 'advect-rotate':
        n = a.n or 64; steps = a.steps or 628
        make_mesh(d, n, n, 0, 1.0, 1.0, 1.0 / n)
        ex = ("           twoPhase = T\n  cl_flow = 'rotate', cl_omega = 1.0, cl_vcenter = 0.5, 0.5\n  cl_ic = 'cylinder', cl_inside = 1\n"
              "  cl_center = 0.5, 0.75, 0.0, cl_radius = 0.15\n  cl_rho1 = 1.0, cl_rho2 = 1.0, cl_mu1 = 1e-3, cl_mu2 = 1e-3\n")
        write_input(d, title='CLSVOF solid body rotation', steps=steps, dt=2 * math.pi / steps, calcU='F', maxit=1, extra=ex)
    elif a.case == 'advect-sphere':
        n = a.n or 32; T = a.T or 3.0; steps = a.steps or 300
        make_mesh(d, n, n, n, 1.0, 1.0, 1.0)
        ex = ("           twoPhase = T\n  cl_flow = 'vortex', cl_T = %g\n  cl_ic = 'sphere', cl_inside = 1\n"
              "  cl_center = 0.35, 0.35, 0.35, cl_radius = 0.15\n  cl_rho1 = 1.0, cl_rho2 = 1.0, cl_mu1 = 1e-3, cl_mu2 = 1e-3\n" % T)
        write_input(d, title='CLSVOF 3D vortex', steps=steps, dt=T / steps, calcU='F', maxit=1, extra=ex)
    elif a.case == 'static-drop':
        n = a.n or 32; nd = a.ndim or 2; steps = a.steps or 200
        R = 0.25; sigma = 1.0
        if nd == 2: make_mesh(d, n, n, 0, 1.0, 1.0, 1.0 / n)
        else:       make_mesh(d, n, n, n, 1.0, 1.0, 1.0)
        ic = 'cylinder' if nd == 2 else 'sphere'
        # Capillary time step limit  dt < sqrt( (rho1+rho2) h^3 / (4 pi sigma) )
        h = 1.0 / n; dt = 0.5 * math.sqrt(2.0 * h ** 3 / (4 * math.pi * sigma))
        ex = ("           twoPhase = T\n  cl_ic = '%s', cl_inside = 1\n  cl_center = 0.5, 0.5, 0.5, cl_radius = %g\n"
              "  cl_rho1 = 1.0, cl_rho2 = 1.0, cl_mu1 = 0.05, cl_mu2 = 0.05, cl_sigma = %g\n" % (ic, R, sigma))
        write_input(d, title='CLSVOF static drop', steps=steps, dt=dt, maxit=a.maxit or 4, extra=ex, densit=1.0)
    elif a.case == 'rising-bubble':
        n = a.n or 40; h = 1.0 / n
        sigma = 24.5 if a.bcase == 1 else 1.96
        rho1, rho2, mu1, mu2 = (1000.0, 100.0, 10.0, 1.0) if a.bcase == 1 else (1000.0, 1.0, 10.0, 0.1)
        steps = a.steps or int(3.0 / (0.5 * h * 0.1) )
        make_mesh(d, n, 2 * n, 0, 1.0, 2.0, h, side='symmetry')   # Hysing: free slip on the vertical walls
        dt = 3.0 / steps
        ex = ("           twoPhase = T\n  cl_ic = 'cylinder', cl_inside = 2\n  cl_center = 0.5, 0.5, 0.0, cl_radius = 0.25\n"
              "  cl_rho1 = %g, cl_rho2 = %g, cl_mu1 = %g, cl_mu2 = %g, cl_sigma = %g\n" % (rho1, rho2, mu1, mu2, sigma))
        write_input(d, title='Hysing rising bubble case %d' % a.bcase, steps=steps, dt=dt, gy=-0.98, maxit=a.maxit or 5,
                    extra=ex, densit=rho1, nzapis=max(1, steps // 6))
    else:
        # Falling drop in a closed box (density ratio 1000, laminar gas: Re_gas ~ 30, We_gas ~ 1): g = 1, sigma = 1, R = 0.2.
        n = a.n or 24; nd = a.ndim or 3; h = 1.0 / n
        make_mesh(d, n, 2 * n, n, 1.0, 2.0, 1.0) if nd == 3 else make_mesh(d, n, 2 * n, 0, 1.0, 2.0, h)
        steps = a.steps or 300; dt = 0.004
        ic = 'sphere' if nd == 3 else 'cylinder'
        ex = ("           twoPhase = T\n  cl_ic = '%s', cl_inside = 1\n  cl_center = 0.5, 1.5, 0.5, cl_radius = 0.2\n"
              "  cl_rho1 = 1000.0, cl_rho2 = 1.0, cl_mu1 = 0.1, cl_mu2 = 0.01, cl_sigma = 1.0\n" % ic)
        write_input(d, title='CLSVOF falling drop', steps=steps, dt=dt, gy=-1.0, maxit=a.maxit or 5, extra=ex,
                    densit=1.0, nzapis=max(1, steps // 4))


if __name__ == '__main__':
    main()
