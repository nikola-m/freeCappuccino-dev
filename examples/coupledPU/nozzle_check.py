#!/usr/bin/env python3
"""nozzle_check.py <casedir> <restart>: Mach number, p/p0 and T/T0 along the duct of mknozzle.py against the quasi-1D
isentropic solution (area ratio h(x)/h_in, subsonic branch), and the mass flow through the first and the last cell column."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from postairfoil import records
Rg, cp = 287.058, 1006.0
gam = cp / (cp - Rg)
case, rst = sys.argv[1], sys.argv[2]
info = {l.split()[0]: float(l.split()[1]) for l in open(case + '/info.txt')}
nx, ny = int(info['nx']), int(info['ny'])
size = [int(x) for x in open(case + '/polyMesh/size').read().split()]
nc, nB = size[1], size[3]
arr = [np.frombuffer(x, dtype='<f8') for x in records(rst)[1:]]
big = [x for x in arr if len(x) == nc + nB]
u, v, w, p, te, ed, T = big[:7]
x = (np.arange(nx) + 0.5) / nx * info['L']
h = lambda x: info['H'] * (1 - info['a'] * np.sin(np.pi * x / info['L']) ** 2)
Af = lambda M: (1 / M) * (2 / (gam + 1) * (1 + 0.5 * (gam - 1) * M * M)) ** ((gam + 1) / (2 * (gam - 1)))
Ast = h(0.0) / Af(info['Min'])


def Msub(A):
    lo, hi = 1e-4, 1.0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if Af(m) > A / Ast: lo = m
        else: hi = m
    return 0.5 * (lo + hi)


Mex = np.array([Msub(h(xx)) for xx in x])
col = lambda q: q[:nc].reshape(ny, nx).mean(0)
M = col(np.sqrt(u**2 + v**2) / np.sqrt(gam * Rg * T))
f = 1 + 0.5 * (gam - 1) * Mex ** 2
pex = info['p0'] * f ** (-gam / (gam - 1))
Tex = info['T0'] / f
print('max |M - M_1D| = %.4f  (M_throat: %.4f vs %.4f)   max |p/p_1D - 1| = %.4f   max |T/T_1D - 1| = %.4f'
      % (abs(M - Mex).max(), M.min() if False else M[nx // 2], Mex[nx // 2], abs(col(p) / pex - 1).max(), abs(col(T) / Tex - 1).max()))
rho = p[:nc] / (Rg * T[:nc])
for name, c in (('first', 0), ('last', nx - 1)):
    idx = np.arange(ny) * nx + c
    print('mass flow through %s column: %.5f kg/s (design %.5f)' % (name, (rho[idx] * u[idx]).sum() / ny * h(x[c]) * info['thick'], info['mdot']))
