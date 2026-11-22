#!/usr/bin/env python3
"""couette_check.py <casedir> <restart>: compare u and T in the last cell column with the exact Couette solution."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from postairfoil import records
cp = 1006.0
case, rst = sys.argv[1], sys.argv[2]
info = {l.split()[0]: float(l.split()[1]) for l in open(case + '/info.txt')}
nx, ny = int(info['nx']), int(info['ny'])
size = [int(x) for x in open(case + '/polyMesh/size').read().split()]
nc, nB = size[1], size[3]
arr = [np.frombuffer(x, dtype='<f8') for x in records(rst)[1:]]
big = [x for x in arr if len(x) == nc + nB]
u, v, w, p, te, ed, T = big[:7]
y = (np.arange(ny) + 0.5) / ny * info['h']
eta = y / info['h']
Tex = info['Tw'] + info['Pr'] * info['U']**2 / (2 * cp) * eta * (1 - eta)
for name, col in (('outlet column', nx - 1), ('mid column', nx // 2)):
    idx = np.arange(ny) * nx + col                       # cell id = i + nx*j
    print('%s:  max|u-u_ex|/U = %.2e   max|T-T_ex| = %.2e K  (heating amplitude %.3f K, T_ex(mid) = %.4f, T(mid) = %.4f)'
          % (name, abs(u[idx] - info['U'] * eta).max() / info['U'], abs(T[idx] - Tex).max(),
             info['Pr'] * info['U']**2 / (8 * cp), Tex[ny // 2], T[idx][ny // 2]))
print('p range', p[:nc].min(), p[:nc].max())
