#!/usr/bin/env python3
"""bump_check.py <casedir> <restart>: summary of the bump channel result (mkbump.py).

 - entropy  s = (p/p_inf)/(rho/rho_inf)^gamma - 1   (zero for the exact isentropic flow; for M < 1 its maximum is the numerical entropy production)
 - Mach number: maximum, and its location on the lower wall
 - lower wall pressure: minimum p/p_inf, the strongest compression (shock) on the bump and its x position, p/p_inf at the exit
 - M > 1: oblique shock at the leading edge of the bump (wedge angle of the arc): p/p_inf after the shock from the theta-beta-M relation
"""
import sys, os, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'coupledPU'))
from postairfoil import records
Rg, cp = 287.058, 1006.0
gam = cp / (cp - Rg)
case, rst = sys.argv[1], sys.argv[2]
info = {l.split()[0]: float(l.split()[1]) for l in open(case + '/info.txt')}
nx, ny, M, t = int(info['nx']), int(info['ny']), info['M'], info['t']
size = [int(x) for x in open(case + '/polyMesh/size').read().split()]
nc, nB = size[1], size[3]
arr = [np.frombuffer(x, dtype='<f8') for x in records(rst)[1:]]
big = [x for x in arr if len(x) == nc + nB]
u, v, w, p, te, ed, T = [b[:nc] for b in big[:7]]
pinf, Tinf = info['pinf'], info['Tinf']
rho, rhoinf = p / (Rg * T), pinf / (Rg * Tinf)
Mach = np.sqrt(u * u + v * v) / np.sqrt(gam * Rg * T)
s = (p / pinf) / (rho / rhoinf) ** gam - 1.0
xn = np.linspace(0, 3.0, nx + 1); xc = 0.5 * (xn[1:] + xn[:-1])
ok = np.isfinite(p).all()
print('finite: %s' % ok)
if not ok: sys.exit(1)
pw = p[:nx] / pinf; Mw = Mach[:nx]                       # first cell row = lower wall
print('Mach: max %.4f (lower wall: max %.4f at x = %.3f, min %.4f)' % (Mach.max(), Mw.max(), xc[Mw.argmax()], Mw.min()))
print('entropy s: max %.2e, min %.2e  (inlet-to-exit mean %.2e)' % (s.max(), s.min(), s.mean()))
dp = np.diff(pw) / np.diff(xc); k = int(dp.argmin())
print('wall p/p_inf: min %.4f at x = %.3f, max %.4f at x = %.3f, exit %.4f' % (pw.min(), xc[pw.argmin()], pw.max(), xc[pw.argmax()], pw[-1]))
print('steepest wall pressure rise %.2f /m at x = %.3f' % (np.diff(pw)[np.diff(pw).argmax()] / (xc[1] - xc[0]), xc[int(np.diff(pw).argmax())]))
if M > 1.0:
    R = (0.25 + t * t) / (2 * t); th = math.asin(0.5 / R)
    # weak oblique shock: theta-beta-M relation by bisection on beta
    def theta(beta): return math.atan(2 / math.tan(beta) * (M * M * math.sin(beta) ** 2 - 1) / (M * M * (gam + math.cos(2 * beta)) + 2))
    bs = np.linspace(math.asin(1 / M) + 1e-6, math.radians(80), 4000)
    bmax = bs[int(np.argmax([theta(b) for b in bs]))]          # weak branch: beta below the one of the maximum deflection
    lo, hi = math.asin(1 / M) + 1e-9, bmax
    if th > theta(bmax): print('wedge angle exceeds the maximum deflection: detached shock'); sys.exit(0)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if theta(mid) < th: lo = mid
        else: hi = mid
    beta = 0.5 * (lo + hi); Mn = M * math.sin(beta)
    pr = 1 + 2 * gam / (gam + 1) * (Mn * Mn - 1)
    print('leading-edge wedge angle %.2f deg, oblique shock beta = %.2f deg, p2/p1 = %.4f (theory); computed wall p/p_inf at x = %.3f: %.4f'
          % (math.degrees(th), math.degrees(beta), pr, xc[np.searchsorted(xc, 1.0) + 2], pw[np.searchsorted(xc, 1.0) + 2]))
