#!/usr/bin/env python3
"""
mkcase.py - generate small structured test cases for the block-coupled p-U solver, in the native
freeCappuccino 'nativeMesh' format (polyMesh/, 0/U, input-simple.nml, input-coupled.nml).

  python3 mkcase.py cavity   <dir> [--nx 32 --ny 32 --Re 100 --stretch 0 --shear 0]
  python3 mkcase.py channel  <dir> [--nx 40 --ny 10 --Lx 10 --Re 50 --outlet pressure|outlet]
  python3 mkcase.py cavity3d <dir> [--nx 16 --Re 100]

cavity    lid-driven square cavity (2D, one cell thick, 'empty' front/back patches)
channel   plane channel with velocity inlet; outlet is either a fixed-pressure ('pressure') or a
          zero-gradient, mass-flux-scaled ('outlet') boundary
cavity3d  lid-driven cube
--stretch clusters cells towards the walls (tanh), --shear skews the mesh (x += shear*y), which makes it
          non-orthogonal: useful to test robustness.
"""
import argparse, os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtkmesh import write_vtk_templates


def tanh_stretch(n, s):
    xi = np.linspace(0.0, 1.0, n + 1)
    if s == 0:
        return xi
    return 0.5 * (1.0 + np.tanh(s * (2 * xi - 1)) / math.tanh(s))


def build_mesh(xn, yn, skew_shear=0.0, patches=None):
    """xn, yn: 1-D arrays of grid lines. Returns dict with points, faces, owner, neighbour, boundary."""
    nx, ny = len(xn) - 1, len(yn) - 1
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    pts = []
    for k in range(2):
        for j in range(ny + 1):
            for i in range(nx + 1):
                pts.append((xn[i] + skew_shear * yn[j], yn[j], 0.1 * k))
    cid = lambda i, j: 1 + i + nx * j

    inner = []   # (owner, neighbour, nodes) ; owner < neighbour
    # x-faces between cells (i-1,j) and (i,j)
    for j in range(ny):
        for i in range(1, nx):
            L, R = cid(i - 1, j), cid(i, j)
            nodes = [nid(i, j, 0), nid(i, j + 1, 0), nid(i, j + 1, 1), nid(i, j, 1)]
            inner.append((L, R, nodes))
    # y-faces between (i,j-1) and (i,j)
    for j in range(1, ny):
        for i in range(nx):
            L, R = cid(i, j - 1), cid(i, j)
            nodes = [nid(i, j, 0), nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j, 0)]
            inner.append((L, R, nodes))
    inner.sort(key=lambda t: (t[0], t[1]))

    bnd = {'left': [], 'right': [], 'bottom': [], 'top': [], 'frontback': []}
    for j in range(ny):
        i = 0
        bnd['left'].append((cid(i, j), [nid(i, j, 0), nid(i, j, 1), nid(i, j + 1, 1), nid(i, j + 1, 0)]))
        i = nx
        bnd['right'].append((cid(i - 1, j), [nid(i, j, 0), nid(i, j + 1, 0), nid(i, j + 1, 1), nid(i, j, 1)]))
    for i in range(nx):
        j = 0
        bnd['bottom'].append((cid(i, j), [nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j, 1), nid(i, j, 1)]))
        j = ny
        bnd['top'].append((cid(i, j - 1), [nid(i, j, 0), nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j, 0)]))
    for j in range(ny):
        for i in range(nx):
            bnd['frontback'].append((cid(i, j), [nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0)]))
    for j in range(ny):
        for i in range(nx):
            bnd['frontback'].append((cid(i, j), [nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]))
    return dict(nx=nx, ny=ny, pts=pts, inner=inner, bnd=bnd)


def write_native(root, mesh, patch_list):
    """patch_list: [(name, type, key)] in order; 'key' indexes mesh['bnd']"""
    os.makedirs(os.path.join(root, 'polyMesh'), exist_ok=True)
    pm = os.path.join(root, 'polyMesh')
    pts = mesh['pts']
    inner = mesh['inner']
    faces, owner, neigh = [], [], []
    for (L, R, nodes) in inner:
        faces.append(nodes); owner.append(L); neigh.append(R)
    n_inner = len(inner)
    bnd_lines = []
    start = n_inner
    for (name, typ, key) in patch_list:
        lst = mesh['bnd'][key]
        bnd_lines.append((name, typ, len(lst), start))
        for (c, nodes) in lst:
            faces.append(nodes); owner.append(c)
        start += len(lst)
    nfaces = len(faces)
    with open(os.path.join(pm, 'points'), 'w') as f:
        for p in pts:
            f.write('%.16g %.16g %.16g\n' % p)
    with open(os.path.join(pm, 'faces'), 'w') as f:
        for nds in faces:
            f.write('%d %s\n' % (len(nds), ' '.join(map(str, nds))))
    with open(os.path.join(pm, 'owner'), 'w') as f:
        f.write('\n'.join(map(str, owner)) + '\n')
    with open(os.path.join(pm, 'neighbour'), 'w') as f:
        f.write('\n'.join(map(str, neigh)) + '\n')
    with open(os.path.join(pm, 'boundary'), 'w') as f:
        f.write('# bcName bcType nFaces startFace\n')
        for (name, typ, n, s) in bnd_lines:
            f.write('%s %s %d %d\n' % (name, typ, n, s))
    ncells = max(owner)
    with open(os.path.join(pm, 'size'), 'w') as f:
        f.write('%d\n%d\n%d\n%d\n%d\n' % (len(pts), ncells, n_inner, nfaces - n_inner, nfaces))
    # ParaView templates vtk/mesh/mesh_0_0.vtu, mesh_<ib>_0.vtp needed by the VTK writer of the solver
    write_vtk_templates(root)
    return dict(ncells=ncells, n_inner=n_inner, nfaces=nfaces)


def write_U(root, init, patch_bc):
    os.makedirs(os.path.join(root, '0'), exist_ok=True)
    with open(os.path.join(root, '0', 'U'), 'w') as f:
        f.write('internalField\n  uniform\n    %g %g %g\nboundaryField\n' % tuple(init))
        for name, bc in patch_bc:
            if bc is None:
                f.write('%s\n  Neumann\n    zeroGradient\n' % name)
            else:
                f.write('%s\n  Dirichlet\n    uniform\n      %g %g %g\n' % ((name,) + tuple(bc)))


INPUT_TEMPLATE = """&input_parameters
              title = '{title}'
        mesh_format = 'nativeMesh'
              lread = F
              ltest = F
              calcU = T
               urfU = {urfU}, {urfU}, {urfU}
               gdsU = 1.0
           cSchemeU = '{cscheme}'
           lSolverU = 'bicgstab'
           maxiterU = 5
            tolAbsU = 1e-13
            tolRelU = 0.01
            pscheme = 'linear'
              calcP = T
               urfP = {urfP}
           lSolverP = 'iccg'
           maxiterP = 50
            tolAbsP = 1e-12
            tolRelP = 0.025
          TurbModel = 'none'
             densit = 1.0
             viscos = {nu}
              pranl = 1.0
               beta = 1e-3
         ltransient = F
              gauss = T
            limiter = 'no-limit'
             SIMPLE = {simple}
              ncorr = 1
              npcor = 1
           pRefCell = 1
          tolerance = {tol}
            NUMSTEP = {numstep}
           TIMESTEP = 1.0
             NZAPIS = {numstep}
              MAXIT = 1
{extra}/
"""


def write_input(root, name, **kw):
    kw.setdefault('urfU', 0.8); kw.setdefault('urfP', 0.3); kw.setdefault('cscheme', 'linearUpwind')
    kw.setdefault('simple', 'T'); kw.setdefault('tol', '1e-9'); kw.setdefault('numstep', 2000)
    kw.setdefault('extra', '')
    with open(os.path.join(root, name), 'w') as f:
        f.write(INPUT_TEMPLATE.format(**kw))


def build_mesh3d(n):
    nx = ny = nz = n
    x = [i / n for i in range(n + 1)]
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    cid = lambda i, j, k: 1 + i + nx * (j + ny * k)
    pts = [(x[i], x[j], x[k]) for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    fx = lambda i, j, k: [nid(i, j, k), nid(i, j + 1, k), nid(i, j + 1, k + 1), nid(i, j, k + 1)]   # normal +x
    fy = lambda i, j, k: [nid(i, j, k), nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j, k)]   # normal +y
    fz = lambda i, j, k: [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]   # normal +z
    inner, lid, walls = [], [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = cid(i, j, k)
                (inner.append((c, cid(i + 1, j, k), fx(i + 1, j, k))) if i < nx - 1 else walls.append((c, fx(i + 1, j, k))))
                (inner.append((c, cid(i, j + 1, k), fy(i, j + 1, k))) if j < ny - 1 else lid.append((c, fy(i, j + 1, k))))
                (inner.append((c, cid(i, j, k + 1), fz(i, j, k + 1))) if k < nz - 1 else walls.append((c, fz(i, j, k + 1))))
                if i == 0: walls.append((c, fx(i, j, k)[::-1]))
                if j == 0: walls.append((c, fy(i, j, k)[::-1]))
                if k == 0: walls.append((c, fz(i, j, k)[::-1]))
    inner.sort(key=lambda t: (t[0], t[1]))
    return dict(pts=pts, inner=inner, bnd={'lid': lid, 'walls': walls})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case', choices=['cavity', 'channel', 'cavity3d'])
    ap.add_argument('dir')
    ap.add_argument('--nx', type=int, default=None); ap.add_argument('--ny', type=int, default=None)
    ap.add_argument('--stretch', type=float, default=0.0)
    ap.add_argument('--shear', type=float, default=0.0)
    ap.add_argument('--Lx', type=float, default=10.0)
    ap.add_argument('--outlet', default='pressure', choices=['pressure', 'outlet'])
    ap.add_argument('--Re', type=float, default=100.0)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.dir, 'vtk', 'mesh'), exist_ok=True)
    nu = 1.0 / a.Re
    if a.case == 'cavity':
        nx, ny = a.nx or 32, a.ny or a.nx or 32
        mesh = build_mesh(tanh_stretch(nx, a.stretch), tanh_stretch(ny, a.stretch), a.shear)
        patches = [('lid', 'wall', 'top'), ('wall-l', 'wall', 'left'), ('wall-r', 'wall', 'right'),
                   ('wall-b', 'wall', 'bottom'), ('front-back', 'empty', 'frontback')]
        info = write_native(a.dir, mesh, patches)
        write_U(a.dir, (0, 0, 0), [('lid', (1, 0, 0)), ('wall-l', (0, 0, 0)), ('wall-r', (0, 0, 0)),
                                   ('wall-b', (0, 0, 0)), ('front-back', None)])
        kw = dict(urfU=0.8)
    elif a.case == 'channel':
        nx, ny = a.nx or 40, a.ny or 10
        mesh = build_mesh(a.Lx * tanh_stretch(nx, 0.0), tanh_stretch(ny, a.stretch), a.shear)
        patches = [('inlet', 'inlet', 'left'), ('outlet', a.outlet, 'right'), ('wall-t', 'wall', 'top'),
                   ('wall-b', 'wall', 'bottom'), ('front-back', 'empty', 'frontback')]
        info = write_native(a.dir, mesh, patches)
        write_U(a.dir, (1, 0, 0), [('inlet', (1, 0, 0)), ('outlet', None), ('wall-t', (0, 0, 0)),
                                   ('wall-b', (0, 0, 0)), ('front-back', None)])
        kw = dict(urfU=0.5, urfP=0.2)
    else:
        mesh = build_mesh3d(a.nx or 16)
        info = write_native(a.dir, mesh, [('lid', 'wall', 'lid'), ('walls', 'wall', 'walls')])
        write_U(a.dir, (0, 0, 0), [('lid', (1, 0, 0)), ('walls', (0, 0, 0))])
        kw = dict(urfU=0.8)
    print(a.dir, info, 'viscosity = %g' % nu)
    write_input(a.dir, 'input-simple.nml', title=a.case, nu=nu, simple='T', numstep=20000, **kw)
    kw['urfU'] = 0.95
    write_input(a.dir, 'input-coupled.nml', title=a.case + ' (block-coupled p-U)', nu=nu, simple='F',
                numstep=20000, extra='             COUPLED = T\n', **kw)


if __name__ == '__main__':
    main()
