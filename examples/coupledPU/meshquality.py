#!/usr/bin/env python3
"""
meshquality.py <case dir> [--plot mesh.png] [--zoom x0 x1 y0 y1]

Quality report of the native mesh in <case>/polyMesh (hexahedra, e.g. the 2D cases of the examples):
  * non-orthogonality  : angle between the face normal and the line joining the two cell centres (inner faces)
  * skewness           : distance of the intersection of that line with the face plane from the face centre, / |d|
  * aspect ratio       : longest / shortest edge of the cell
  * volume ratio       : larger / smaller volume of the two cells sharing a face
  * non-positive cells : cells with a non-positive Jacobian in the VTK node ordering
--plot writes a picture of the mesh (z = 0 layer; needs matplotlib), --zoom limits it to a window.
"""
import argparse, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtkmesh import read_polymesh, cell_connectivity, hex_is_positive


def quality(case):
    pts, faces, owner, neigh, bnd = read_polymesh(case)
    P = np.array(pts)
    cells = cell_connectivity(pts, faces, owner, neigh)
    C = np.array([P[c].mean(0) for c in cells])
    # cell volumes: divergence theorem with the outward faces
    vol = np.zeros(len(cells))
    fc, fa = [], []
    for i, nds in enumerate(faces):
        q = P[np.array(nds) - 1]
        c0 = q.mean(0)
        a = np.zeros(3)
        for k in range(len(nds)):
            a += 0.5 * np.cross(q[k] - c0, q[(k + 1) % len(nds)] - c0)
        fc.append(c0); fa.append(a)
        vol[owner[i] - 1] += np.dot(a, c0) / 3.0
        if i < len(neigh):
            vol[neigh[i] - 1] -= np.dot(a, c0) / 3.0
    fc, fa = np.array(fc), np.array(fa)
    ni = len(neigh)
    O = np.array(owner[:ni]) - 1
    N = np.array(neigh) - 1
    d = C[N] - C[O]
    dn = np.linalg.norm(d, axis=1)
    an = np.linalg.norm(fa[:ni], axis=1)
    cosang = np.sum(d * fa[:ni], axis=1) / (dn * an)
    nonorth = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
    # intersection of the line C_O -> C_N with the face plane
    t = np.sum((fc[:ni] - C[O]) * fa[:ni], axis=1) / np.sum(d * fa[:ni], axis=1)
    xi = C[O] + t[:, None] * d
    skew = np.linalg.norm(xi - fc[:ni], axis=1) / dn
    # aspect ratio from the 12 edges
    asp = np.zeros(len(cells))
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    for ic, c in enumerate(cells):
        e = [np.linalg.norm(P[c[a]] - P[c[b]]) for a, b in edges]
        asp[ic] = max(e) / min(e)
    vr = np.maximum(vol[O], vol[N]) / np.minimum(vol[O], vol[N])
    bad = [ic + 1 for ic, c in enumerate(cells) if not hex_is_positive([P[n] for n in c])]
    return dict(P=P, cells=cells, vol=vol, nonorth=nonorth, skew=skew, asp=asp, vr=vr, bad=bad, C=C, O=O, N=N)


def report(q, name=''):
    print('%scells %d  volume min/max %.3g / %.3g  non-positive cells: %d' % (name, len(q['vol']), q['vol'].min(), q['vol'].max(), len(q['bad'])))
    for key, lab in (('nonorth', 'non-orthogonality [deg]'), ('skew', 'skewness'), ('asp', 'aspect ratio'), ('vr', 'neighbour volume ratio')):
        v = q[key]
        print('  %-26s max %9.4g   mean %9.4g   99.9%% %9.4g' % (lab, v.max(), v.mean(), np.percentile(v, 99.9)))
    return q


def plot(q, fname, zoom=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    P = q['P']
    segs = []
    zmin = P[:, 2].min()
    for c in q['cells']:
        xy = np.array([P[n][:2] for n in c if abs(P[n][2] - zmin) < 1e-12])       # the 4 nodes of the z = min layer
        if len(xy) != 4:
            continue
        ctr = xy.mean(0)
        xy = xy[np.argsort(np.arctan2(xy[:, 1] - ctr[1], xy[:, 0] - ctr[0]))]
        segs.append([tuple(p) for p in np.vstack([xy, xy[:1]])])
    fig, ax = plt.subplots(figsize=(10, 6), dpi=130)
    ax.add_collection(LineCollection(segs, linewidths=0.25, colors='k'))
    ax.set_aspect('equal')
    if zoom:
        ax.set_xlim(zoom[0], zoom[1]); ax.set_ylim(zoom[2], zoom[3])
    else:
        ax.autoscale()
    fig.savefig(fname, bbox_inches='tight')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case'); ap.add_argument('--plot'); ap.add_argument('--zoom', type=float, nargs=4)
    a = ap.parse_args()
    q = report(quality(a.case), a.case + ': ')
    if a.plot:
        plot(q, a.plot, a.zoom)
