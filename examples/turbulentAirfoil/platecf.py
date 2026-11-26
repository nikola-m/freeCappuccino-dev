#!/usr/bin/env python3
"""
platecf.py - skin friction of the flat plate case (mkplate.py) from the restart file.

  python3 platecf.py <case-dir> [restart-file]

Cf(x) = 2 nu u_1 / (y_1 U^2) from the first cell (the cell height is taken from the mesh), compared with
Blasius, Cf = 0.664/sqrt(Re_x), and the turbulent correlation Cf = 0.0576 Re_x^-0.2.
Cell ordering of the generated mesh: i (streamwise) fastest.
"""
import sys, struct, os
import numpy as np

case = sys.argv[1]
rst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(case, 'restart')
nu = None; U = None
for l in open(os.path.join(case, 'input.nml')):
    if 'viscos' in l: nu = float(l.split('=')[1])
    if 'internalField' in l: pass
U = float(open(os.path.join(case, '0', 'U')).read().split()[2])

pts = np.loadtxt(os.path.join(case, 'polyMesh', 'points'))
size = [int(x) for x in open(os.path.join(case, 'polyMesh', 'size')).read().split()]
ncells = size[1]
xs = np.unique(np.round(pts[:, 0], 10)); ys = np.unique(np.round(pts[:, 1], 10))
nx, ny = len(xs) - 1, len(ys) - 1
assert nx * ny == ncells
xc = 0.5 * (xs[1:] + xs[:-1]); y1 = 0.5 * (ys[1] - ys[0])

b = open(rst, 'rb').read(); o = 0; rec = []
while o < len(b):
    n = struct.unpack('<i', b[o:o + 4])[0]; rec.append(np.frombuffer(b[o + 4:o + 4 + n], '<f8') if n % 8 == 0 else np.zeros(0)); o += 8 + n
arr = [r for r in rec[1:] if len(r) >= ncells]
u = arr[1][:ncells].reshape(ny, nx)
cf = 2 * nu * u[0, :] / (y1 * U**2)
print('   x[m]     Re_x      Cf(num)    Cf Blasius   Cf turb.')
for i in range(0, nx, max(1, nx // 25)):
    if xc[i] <= 0: continue
    rex = U * xc[i] / nu
    print('%8.4f %9.3e  %9.5f   %9.5f   %9.5f' % (xc[i], rex, cf[i], 0.664 / np.sqrt(rex), 0.0576 * rex**-0.2))
