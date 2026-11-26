#!/usr/bin/env python3
"""
postchannel.py - skin friction of the developing channel (mkchannel_turb.py) from the restart file, at the second last column.

  python3 postchannel.py <case-dir> [restart-file]

Cf = 2 tau_w/(rho U_b^2). tau_w from the first cell: viscous sublayer, tau = nu u_1/y_1 (y+ < 11), otherwise from the log law
u_1/u_tau = ln(9.8 y_1 u_tau/nu)/0.41 (wall functions). Compared with Dean's correlation Cf = 0.073 Re_H^-0.25. Prints u+ against y+
for the low-Re meshes.
"""
import sys, struct, os
import numpy as np

case = sys.argv[1]; rst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(case, 'rst')
if not os.path.exists(rst): rst = os.path.join(case, 'restart')
nu = float([l for l in open(os.path.join(case, 'input-turb.nml')) if 'viscos' in l][0].split('=')[1])
pts = np.loadtxt(os.path.join(case, 'polyMesh', 'points'))
xs = np.unique(np.round(pts[:, 0], 10)); ys = np.unique(np.round(pts[:, 1], 10))
nx, ny = len(xs) - 1, len(ys) - 1; nc = nx * ny
b = open(rst, 'rb').read(); o = 0; rec = []
while o < len(b):
    n = struct.unpack('<i', b[o:o + 4])[0]; rec.append(b[o + 4:o + 4 + n]); o += 8 + n
arr = [np.frombuffer(r, '<f8') for r in rec[1:] if len(r) % 8 == 0 and len(r) // 8 >= nc]
u = arr[1][:nc].reshape(ny, nx)[:, nx - 2]
H = ys[-1] - ys[0]; Ub = np.sum(u * np.diff(ys)) / H; y1 = 0.5 * (ys[1] - ys[0])
ut = np.sqrt(nu * u[0] / y1)
if ut * y1 / nu > 11.0:
    ut = 0.05
    for _ in range(100): ut = u[0] * 0.41 / np.log(9.8 * y1 * ut / nu)
cf = 2 * ut**2 / Ub**2
print('Ub %.4f  u_tau %.5f  y1+ %.2f  Cf %.5f   Dean Cf %.5f   (%+.1f %%)' %
      (Ub, ut, y1 * ut / nu, cf, 0.073 * (Ub * H / nu)**-0.25, 100 * (cf / (0.073 * (Ub * H / nu)**-0.25) - 1)))
if y1 * ut / nu < 3:
    yc = 0.5 * (ys[1:] + ys[:-1])
    for k in range(0, ny // 2, 6):
        print('  y+ = %8.2f   u+ = %7.3f   log law %7.3f' % (yc[k] * ut / nu, u[k] / ut, np.log(yc[k] * ut / nu) / 0.41 + 5.2))
