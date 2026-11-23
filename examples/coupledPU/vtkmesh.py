#!/usr/bin/env python3
"""
vtkmesh.py - writes the ParaView mesh templates needed by the VTK writer of freeCappuccino (src/io/writefiles.f90)
from the native mesh in <case>/polyMesh:

    vtk/mesh/mesh_0_0.vtu         unstructured grid: points + cell connectivity + offsets + types (the 'cells' file that
                                  the native format does not have is reconstructed from faces/owner/neighbour)
    vtk/mesh/mesh_<ib>_0.vtp      one polydata file per boundary region ib = 1..numBoundaries (order of polyMesh/boundary):
                                  the points of the region (local numbering) and its faces (as listed in 'faces')

The solver reads these files line by line: the first 4 (.vtu) / 6 (.vtp) lines are copied before the fields are written
(up to and including '<PointData>' / '<CellData>'), the remainder after them, so the layout below must not be changed.
All lines are shorter than the 200 characters the writer can read.

Usage:   python3 vtkmesh.py <case dir>            (an existing case: writes <case>/vtk/mesh/...)
         from vtkmesh import write_vtk_templates   (called by mkcase.write_native for every generated case)

Cell types: hexahedra only (VTK type 12, structured and O/C-type hexahedral meshes). Other cell shapes stop with a message.
Node, cell and face numbers in polyMesh are 1-based (as everywhere in the native format); VTK is 0-based.
"""
import os
import sys

VTK_HEXAHEDRON = 12


def _read_ints(path):
    with open(path) as f:
        return [int(t) for t in f.read().split()]


def read_polymesh(case):
    pm = os.path.join(case, 'polyMesh')
    pts = []
    with open(os.path.join(pm, 'points')) as f:
        for line in f:
            t = line.split()
            if len(t) >= 3:
                pts.append((float(t[0]), float(t[1]), float(t[2])))
    faces = []
    with open(os.path.join(pm, 'faces')) as f:
        for line in f:
            t = line.split()
            if t:
                n = int(t[0])
                faces.append([int(x) for x in t[1:1 + n]])
    owner = _read_ints(os.path.join(pm, 'owner'))
    neigh = _read_ints(os.path.join(pm, 'neighbour'))
    bnd = []
    with open(os.path.join(pm, 'boundary')) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith('#'):
                t = line.split()
                bnd.append((t[0], t[1], int(t[2]), int(t[3])))
    return pts, faces, owner, neigh, bnd


def hex_connectivity(cell_faces):
    """cell_faces: six quads (node lists, outward orientation) of one cell -> 8 nodes in VTK_HEXAHEDRON order."""
    if len(cell_faces) != 6 or any(len(f) != 4 for f in cell_faces):
        raise ValueError('not a hexahedron (%d faces, sizes %s)' % (len(cell_faces), [len(f) for f in cell_faces]))
    bottom = cell_faces[0]
    bset = set(bottom)
    opposite = [f for f in cell_faces[1:] if not (set(f) & bset)]
    if len(opposite) != 1:
        raise ValueError('not a hexahedron: no unique face opposite to the first one')
    oset = set(opposite[0])
    # edges of the four side faces: consecutive nodes of the cyclic node lists
    up = {}
    for f in cell_faces[1:]:
        if f is opposite[0]:
            continue
        for k in range(4):
            a, b = f[k], f[(k + 1) % 4]
            if a in bset and b in oset:
                up[a] = b
            elif b in bset and a in oset:
                up[b] = a
    bot = bottom[::-1]       # VTK: nodes 0-3 counter-clockwise seen from the top face, i.e. normal into the cell
    if any(n not in up for n in bot):
        raise ValueError('inconsistent hexahedron topology')
    return bot + [up[n] for n in bot]


def _sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def _cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def _dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def hex_is_positive(p):
    """Jacobian sign at the 8 corners of a VTK hexahedron (p: 8 points). True if all are positive."""
    for i in range(4):
        n, pr = (i + 1) % 4, (i + 3) % 4
        if _dot(_cross(_sub(p[n], p[i]), _sub(p[pr], p[i])), _sub(p[i + 4], p[i])) <= 0.0:
            return False
        if _dot(_cross(_sub(p[4 + pr], p[4 + i]), _sub(p[4 + n], p[4 + i])), _sub(p[i], p[4 + i])) <= 0.0:
            return False
    return True


def cell_connectivity(pts, faces, owner, neigh):
    """List of 0-based node lists (VTK hexahedra), one per cell, in cell order 1..numCells."""
    ncells = max(owner)
    cf = [[] for _ in range(ncells)]
    for i, nds in enumerate(faces):
        cf[owner[i] - 1].append(nds)                  # normal points out of the owner
        if i < len(neigh):
            cf[neigh[i] - 1].append(nds[::-1])        # reversed for the neighbour: outward again
    cells = []
    for c, fl in enumerate(cf):
        try:
            h = hex_connectivity(fl)
        except ValueError as e:
            raise SystemExit('vtkmesh: cell %d: %s. Only hexahedral cells are supported.' % (c + 1, e))
        cells.append([n - 1 for n in h])
        if not hex_is_positive([pts[n] for n in cells[-1]]):
            print('vtkmesh: warning: cell %d has a non-positive Jacobian in the VTK ordering' % (c + 1))
    return cells


def _w_ints(f, vals, per_line):
    for i in range(0, len(vals), per_line):
        f.write('          ' + ' '.join(str(v) for v in vals[i:i + per_line]) + '\n')


def _w_points(f, plist):
    f.write('        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">\n')
    for p in plist:
        f.write('          %.9g %.9g %.9g\n' % p)
    f.write('        </DataArray>\n')


def write_vtu(path, pts, cells):
    with open(path, 'w') as f:
        f.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write('    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n' % (len(pts), len(cells)))
        f.write('      <PointData>\n      </PointData>\n      <CellData>\n      </CellData>\n')   # lines 4-7: the solver inserts data here
        f.write('      <Points>\n')
        _w_points(f, pts)
        f.write('      </Points>\n      <Cells>\n')
        f.write('        <DataArray type="Int64" IdType="1" Name="connectivity" format="ascii">\n')
        for c in cells:
            f.write('          ' + ' '.join(map(str, c)) + '\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int64" IdType="1" Name="offsets" format="ascii">\n')
        _w_ints(f, [8 * (i + 1) for i in range(len(cells))], 10)
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        _w_ints(f, [VTK_HEXAHEDRON] * len(cells), 20)
        f.write('        </DataArray>\n      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n')


def write_vtp(path, pts, region_faces):
    """region_faces: list of 1-based node lists of the faces of one boundary region (outward orientation)."""
    used = sorted({n for fc in region_faces for n in fc})
    loc = {n: i for i, n in enumerate(used)}
    empty = ('      <%s>\n'
             '        <DataArray type="Int64" IdType="1" Name="connectivity" format="ascii">\n        </DataArray>\n'
             '        <DataArray type="Int64" IdType="1" Name="offsets" format="ascii">\n        </DataArray>\n'
             '      </%s>\n')
    with open(path, 'w') as f:
        f.write('<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
        f.write('  <PolyData>\n')
        f.write('    <Piece NumberOfPoints="%d" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="%d">\n'
                % (len(used), len(region_faces)))
        f.write('      <PointData>\n      </PointData>\n      <CellData>\n')                   # lines 4-6 are copied before the data
        f.write('      </CellData>\n      <Points>\n')
        _w_points(f, [pts[n - 1] for n in used])
        f.write('      </Points>\n')
        for tag in ('Verts', 'Lines', 'Strips'):
            f.write(empty % (tag, tag))
        f.write('      <Polys>\n')
        f.write('        <DataArray type="Int64" IdType="1" Name="connectivity" format="ascii">\n')
        for fc in region_faces:
            f.write('          ' + ' '.join(str(loc[n]) for n in fc) + '\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int64" IdType="1" Name="offsets" format="ascii">\n')
        off, offs = 0, []
        for fc in region_faces:
            off += len(fc); offs.append(off)
        _w_ints(f, offs, 10)
        f.write('        </DataArray>\n      </Polys>\n    </Piece>\n  </PolyData>\n</VTKFile>\n')


def write_vtk_templates(case):
    """Writes <case>/vtk/mesh/mesh_0_0.vtu and mesh_<ib>_0.vtp (ib = 1..numBoundaries) from <case>/polyMesh."""
    pts, faces, owner, neigh, bnd = read_polymesh(case)
    out = os.path.join(case, 'vtk', 'mesh')
    os.makedirs(out, exist_ok=True)
    cells = cell_connectivity(pts, faces, owner, neigh)
    write_vtu(os.path.join(out, 'mesh_0_0.vtu'), pts, cells)
    for ib, (name, typ, n, start) in enumerate(bnd, 1):
        write_vtp(os.path.join(out, 'mesh_%d_0.vtp' % ib), pts, faces[start:start + n])
    return dict(ncells=len(cells), npoints=len(pts), nregions=len(bnd))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    r = write_vtk_templates(sys.argv[1])
    print('%s/vtk/mesh: %d points, %d hexahedra, %d boundary regions' % (sys.argv[1], r['npoints'], r['ncells'], r['nregions']))
