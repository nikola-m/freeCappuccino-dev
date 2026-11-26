#!/usr/bin/env python3
"""
mkcase.py - generate moving mesh test cases in the native freeCappuccino 'nativeMesh' format
(polyMesh/, 0/U, vtk/mesh templates, input.nml).

  python3 mkcase.py piston <dir> [--ndim 2|3 --nx 40 --ny 10 --L 2 --H 1 --motion translate|oscillate
                                  --U 0.5 --amp 0.1 --freq 1 --steps 40 --dt 0.01 --solver simple|piso --rbf wendland]
  python3 mkcase.py body   <dir> [--n 8 --motion oscillate|rotate --amp 0.25 --freq 0.5 --omega 1.0 --steps 100 --dt 0.02]

piston  Channel (L x H, one cell thick or a box in 3D) closed by a piston wall on the left (moves with the prescribed motion),
        free-slip walls (symmetry) top/bottom/sides, fixed pressure (p = 0) outlet on the right. The fluid is expelled
        through the outlet, the cells are squeezed by the RBF mesh deformation, volumes change in every time step.
        Exact solution (free slip walls, no viscous effects): plug flow, u = piston velocity (at the time step average),
        p = rho*du/dt*(L - x) in the laboratory frame. See piston_check.py.
body    Square body (side 1) in a closed box (side 4, no-slip walls). The body is a moving wall with prescribed motion:
        oscillation in y (the example of the RBF mesh deformation Python script), or rotation. The box is closed, the
        body is rigid: the volume of the fluid is constant. The check is the space conservation law (GCL error) and the
        mass conservation, see the monitor file.

Cells of the 2D meshes are one cell thick (z = 0 .. 0.1) with 'empty' front/back patches (as in examples/coupledPU).
"""
import argparse, os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'coupledPU'))
from mkcase import build_mesh, write_native, write_U          # noqa: E402  (2D mesh writer reused from coupledPU)


def box3d(nx, ny, nz, Lx, Ly, Lz):
    """Structured hexahedral box, six boundary lists xmin, xmax, ymin, ymax, zmin, zmax. Cell index: 1 + i + nx*(j + ny*k)."""
    X = np.linspace(0, Lx, nx + 1); Y = np.linspace(0, Ly, ny + 1); Z = np.linspace(0, Lz, nz + 1)
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    cid = lambda i, j, k: 1 + i + nx * (j + ny * k)
    pts = [(X[i], Y[j], Z[k]) for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    fx = lambda i, j, k: [nid(i, j, k), nid(i, j + 1, k), nid(i, j + 1, k + 1), nid(i, j, k + 1)]   # normal +x
    fy = lambda i, j, k: [nid(i, j, k), nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j, k)]   # normal +y
    fz = lambda i, j, k: [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]   # normal +z
    inner = []
    b = {k: [] for k in ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = cid(i, j, k)
                (inner.append((c, cid(i + 1, j, k), fx(i + 1, j, k))) if i < nx - 1 else b['xmax'].append((c, fx(i + 1, j, k))))
                (inner.append((c, cid(i, j + 1, k), fy(i, j + 1, k))) if j < ny - 1 else b['ymax'].append((c, fy(i, j + 1, k))))
                (inner.append((c, cid(i, j, k + 1), fz(i, j, k + 1))) if k < nz - 1 else b['zmax'].append((c, fz(i, j, k + 1))))
                if i == 0: b['xmin'].append((c, fx(i, j, k)[::-1]))
                if j == 0: b['ymin'].append((c, fy(i, j, k)[::-1]))
                if k == 0: b['zmin'].append((c, fz(i, j, k)[::-1]))
    inner.sort(key=lambda t: (t[0], t[1]))
    return dict(pts=pts, inner=inner, bnd=b)


def square_in_box(n, half=2.0, hole=0.5, thick=0.1):
    """Cartesian grid on [-half,half]^2 (cells of size 1/n) without the square [-hole,hole]^2; one cell thick.
    Boundary lists: 'outer' (box walls), 'body' (faces of the removed square), 'frontback'."""
    N = int(round(2 * half * n))
    lo = int(round((half - hole) * n)); hi = int(round((half + hole) * n))
    h = 1.0 / n
    nid = lambda i, j, k: 1 + i + (N + 1) * (j + (N + 1) * k)
    pts = [(-half + i * h, -half + j * h, thick * k) for k in range(2) for j in range(N + 1) for i in range(N + 1)]
    exists = lambda i, j: 0 <= i < N and 0 <= j < N and not (lo <= i < hi and lo <= j < hi)
    ids = {}
    for j in range(N):
        for i in range(N):
            if exists(i, j):
                ids[(i, j)] = len(ids) + 1
    xface = lambda i, j: [nid(i, j, 0), nid(i, j + 1, 0), nid(i, j + 1, 1), nid(i, j, 1)]      # normal +x
    yface = lambda i, j: [nid(i, j, 0), nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j, 0)]      # normal +y
    inner, outer, body, fb = [], [], [], []
    for (i, j), c in ids.items():
        for (di, dj, mk, rev) in ((1, 0, lambda: xface(i + 1, j), False), (-1, 0, lambda: xface(i, j), True),
                                  (0, 1, lambda: yface(i, j + 1), False), (0, -1, lambda: yface(i, j), True)):
            ni, nj = i + di, j + dj
            nodes = mk()[::-1] if rev else mk()
            if (ni, nj) in ids:
                if not rev:                                  # each inner face once, from the cell with the lower (i,j)
                    inner.append((c, ids[(ni, nj)], nodes))
            elif 0 <= ni < N and 0 <= nj < N:
                body.append((c, nodes))
            else:
                outer.append((c, nodes))
        fb.append((c, [nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0)]))
    for (i, j), c in ids.items():
        fb.append((c, [nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]))
    inner = [t if t[0] < t[1] else (t[1], t[0], t[2][::-1]) for t in inner]
    inner.sort(key=lambda t: (t[0], t[1]))
    # frontback list: the 'front' faces first, then the 'back' faces (as in build_mesh): order is not important here
    return dict(pts=pts, inner=inner, bnd={'outer': outer, 'body': body, 'frontback': fb})


INPUT = """&input_parameters
              title = '{title}'
        mesh_format = 'nativeMesh'
              lread = F
              lwrite = T
              ltest = F
              calcU = T
               urfU = 0.9, 0.9, 0.9
               gdsU = 1.0
           cSchemeU = 'linearUpwind'
           lSolverU = 'bicgstab'
           maxiterU = 10
            tolAbsU = 1e-14
            tolRelU = 0.001
            pscheme = 'linear'
              calcP = T
               urfP = {urfP}
           lSolverP = 'iccg'
           maxiterP = 400
            tolAbsP = 1e-14
            tolRelP = 0.001
          TurbModel = 'none'
             densit = 1.0
             viscos = {viscos}
         ltransient = T
                bdf = T
              gauss = T
            limiter = 'no-limit'
             SIMPLE = {simple}
               PISO = {piso}
              ncorr = 2
              npcor = 1
           pRefCell = 1
          tolerance = 1e-12
            NUMSTEP = {steps}
           TIMESTEP = {dt}
             NZAPIS = {nzapis}
              MAXIT = {maxit}
        DynamicMesh = T
           mm_patch = '{patch}'
          mm_motion = '{motion}'
{motion_params}             mm_rbf = '{rbf}'
{extra}/
"""


def write_input(root, **kw):
    kw.setdefault('urfP', 0.5); kw.setdefault('viscos', 1e-3); kw.setdefault('extra', '')
    kw.setdefault('nzapis', kw['steps'])
    kw['simple'] = 'T' if kw['solver'] == 'simple' else 'F'
    kw['piso'] = 'T' if kw['solver'] == 'piso' else 'F'
    kw['maxit'] = kw['maxit'] if kw['solver'] == 'simple' else 1
    with open(os.path.join(root, 'input.nml'), 'w') as f:
        f.write(INPUT.format(**kw))


def vec(v):
    return ', '.join('%g' % x for x in v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case', choices=['piston', 'body'])
    ap.add_argument('dir')
    ap.add_argument('--ndim', type=int, default=2)
    ap.add_argument('--nx', type=int, default=40)
    ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--n', type=int, default=8, help='body: cells per unit length')
    ap.add_argument('--L', type=float, default=2.0)
    ap.add_argument('--H', type=float, default=1.0)
    ap.add_argument('--motion', default=None)
    ap.add_argument('--U', type=float, default=0.5, help='piston: velocity of the translation')
    ap.add_argument('--amp', type=float, default=None)
    ap.add_argument('--freq', type=float, default=None)
    ap.add_argument('--omega', type=float, default=1.0, help='body: angular velocity of the rotation')
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--dt', type=float, default=None)
    ap.add_argument('--maxit', type=int, default=8)
    ap.add_argument('--solver', choices=['simple', 'piso'], default='simple')
    ap.add_argument('--rbf', default='wendland')
    ap.add_argument('--viscos', type=float, default=None)
    a = ap.parse_args()
    root = a.dir
    os.makedirs(os.path.join(root, 'vtk', 'mesh'), exist_ok=True)

    if a.case == 'piston':
        motion = a.motion or 'translate'
        steps = a.steps or 40; dt = a.dt or 0.01
        if a.ndim == 2:
            xn = np.linspace(0, a.L, a.nx + 1); yn = np.linspace(0, a.H, a.ny + 1)
            mesh = build_mesh(xn, yn)
            mesh['pts'] = [(p[0], p[1], p[2]) for p in mesh['pts']]
            patches = [('piston', 'wall', 'left'), ('outlet', 'pressure', 'right'), ('bottom', 'symmetry', 'bottom'),
                       ('top', 'symmetry', 'top'), ('frontback', 'empty', 'frontback')]
        else:
            nz = a.ny
            mesh = box3d(a.nx, a.ny, nz, a.L, a.H, a.H)
            patches = [('piston', 'wall', 'xmin'), ('outlet', 'pressure', 'xmax'), ('bottom', 'symmetry', 'ymin'),
                       ('top', 'symmetry', 'ymax'), ('front', 'symmetry', 'zmin'), ('back', 'symmetry', 'zmax')]
        if motion == 'translate':
            u0 = a.U
            mp = '           mm_vel = %s\n' % vec((a.U, 0, 0))
        else:
            amp = a.amp if a.amp is not None else 0.1; freq = a.freq if a.freq is not None else 1.0
            u0 = amp * 2 * math.pi * freq          # u(0) = d/dt A sin(2 pi f t)
            mp = '           mm_amp = %s\n          mm_freq = %g\n' % (vec((amp, 0, 0)), freq)
        info = write_native(root, mesh, patches)
        write_U(root, (u0, 0, 0), [(n, (u0, 0, 0)) if t == 'wall' else (n, None) for n, t, _ in patches])
        write_input(root, title='Moving mesh: piston, %s' % motion, steps=steps, dt=dt, patch='piston', motion=motion,
                    motion_params=mp, rbf=a.rbf, solver=a.solver, maxit=a.maxit, viscos=a.viscos or 1e-3,
                    urfP=0.5 if a.solver == 'simple' else 0.3)
        np.savez(os.path.join(root, 'case.npz'), L=a.L, H=a.H, nx=a.nx, ny=a.ny, ndim=a.ndim, dt=dt, steps=steps,
                 motion=motion, U=a.U, amp=(a.amp if a.amp is not None else 0.1), freq=(a.freq if a.freq is not None else 1.0))
    else:
        motion = a.motion or 'oscillate'
        # rotate: 50 steps of 0.02 with omega = 1 is 1 rad. The RBF deformation inverts the cells near the body
        # at about 90 degrees (large rotations need a better mesh motion), see TECHNICAL_DOCUMENTATION.md.
        steps = a.steps or (50 if motion == 'rotate' else 100); dt = a.dt or 0.02
        mesh = square_in_box(a.n)
        patches = [('outer', 'wall', 'outer'), ('body', 'wall', 'body'), ('frontback', 'empty', 'frontback')]
        if motion == 'oscillate':
            amp = a.amp if a.amp is not None else 0.25; freq = a.freq if a.freq is not None else 0.5
            mp = '           mm_amp = %s\n          mm_freq = %g\n' % (vec((0, amp, 0)), freq)
        else:
            mp = '          mm_omega = %g\n         mm_center = 0, 0, 0\n           mm_axis = 0, 0, 1\n' % a.omega
        info = write_native(root, mesh, patches)
        write_U(root, (0, 0, 0), [('outer', (0, 0, 0)), ('body', (0, 0, 0)), ('frontback', None)])
        write_input(root, title='Moving mesh: square in a closed box, %s' % motion, steps=steps, dt=dt, patch='body',
                    motion=motion, motion_params=mp, rbf=a.rbf, solver=a.solver, maxit=a.maxit, viscos=a.viscos or 0.01)
    print(root, info)


if __name__ == '__main__':
    main()
