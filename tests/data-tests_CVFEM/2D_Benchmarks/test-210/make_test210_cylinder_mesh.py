#!/usr/bin/env python3
# SPDX-FileCopyrightText: (c) 2026 The CSMP++ project
# SPDX-License-Identifier: LGPL-3.0-only
"""test-210-cylinder: the full axisymmetric section of a solid cylinder, x from -200 to 200 m about the axis x = 0,
400 m high (y from -400 to 0), gmsh Frontal-Delaunay triangles of target size 23 m. The axis is not a mesh line, so
elements straddle it (test-210c checks their split in AxisymmetricThickness.h). Families HOST, BOTTOM, RIGHT, LEFT,
TOP as the other test meshes. Needs the gmsh Python module. Usage: python3 make_test210_cylinder_mesh.py [outdir]
@author S. Scott (generated with Claude Code)
"""
import struct
import sys

import gmsh
import numpy as np

R, H, SIZE = 200.0, 400.0, 23.0
TITLE = "%s.asc: full section of a solid cylinder, x -200..200 m about the axis x = 0, 400 m high, gmsh triangles 23 m"


def build():
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("cylinder")
    g = gmsh.model.geo
    pts = [g.addPoint(-R, -H, 0, SIZE), g.addPoint(R, -H, 0, SIZE), g.addPoint(R, 0, 0, SIZE), g.addPoint(-R, 0, 0, SIZE)]
    lines = [g.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    g.addPlaneSurface([g.addCurveLoop(lines)]); g.synchronize()
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    pos = {t: i for i, t in enumerate(tags)}
    etypes, _, enodes = gmsh.model.mesh.getElements(2)
    tri = np.array([pos[t] for t in enodes[list(etypes).index(2)]]).reshape(-1, 3)
    gmsh.finalize()
    px, py = xyz[:, 0], xyz[:, 1]
    out = []
    for a, b, c in tri:   # clockwise, as the other test meshes
        area = (px[b] - px[a]) * (py[c] - py[a]) - (px[c] - px[a]) * (py[b] - py[a])
        out.append((a, b, c) if area < 0 else (a, c, b))
    return px, py, out


def families(px, py, tri):
    edges = {}
    for t in tri:
        for k in range(3):
            edges.setdefault(frozenset((t[k], t[(k + 1) % 3])), []).append(t)
    fam = {"BOTTOM": [], "RIGHT": [], "LEFT": [], "TOP": []}
    for e, owners in edges.items():
        if len(owners) != 1:
            continue
        a, b = sorted(e)
        if abs(py[a] + H) < 1e-6 and abs(py[b] + H) < 1e-6: fam["BOTTOM"].append((a, b))
        elif abs(px[a] - R) < 1e-6 and abs(px[b] - R) < 1e-6: fam["RIGHT"].append((a, b))
        elif abs(px[a] + R) < 1e-6 and abs(px[b] + R) < 1e-6: fam["LEFT"].append((a, b))
        else: fam["TOP"].append((a, b))
    return [("HOST", "TRI_3", tri)] + [(k, "BAR_2", v) for k, v in fam.items()]


def write(base, px, py, fams):
    """ANSYS .asc/.dat/-regions.txt in the layout of the other test meshes (as make_test211_obtuse_mesh.py)."""
    elements = [(8 if t == "TRI_3" else 2, e) for _, t, es in fams for e in es]
    owner = {}
    for gid, (code, nodes) in enumerate(elements):
        keys = [frozenset((nodes[(j + 1) % 3], nodes[(j + 2) % 3])) for j in range(3)] if code == 8 else [frozenset(nodes)]
        for k in keys:
            owner.setdefault(k, []).append(gid)
    pfverts = []
    for gid, (code, nodes) in enumerate(elements):
        if code == 8:
            for j in range(3):
                pfverts.append([g for g in owner[frozenset((nodes[(j + 1) % 3], nodes[(j + 2) % 3]))] if g != gid][0])
    name = base.split("/")[-1]
    lines = [".", TITLE % name, "%d # Number of families" % len(fams),
             "# Objectname        Elementtype        Material-ID   Number of elements"]
    for fname, ftype, es in fams:
        lines.append(f"{fname:<19} {ftype:<18} {0:>6} {len(es):>19}")
    lines.append("# now the elements which make up each object are listed in sequence")
    gid = 0
    for fname, ftype, es in fams:
        lines.append(f"{fname} {ftype} {len(es)}")
        ids = list(range(gid, gid + len(es)))
        gid += len(es)
        for k in range(0, len(ids), 10):
            lines.append(" ".join(str(v) for v in ids[k:k + 10]) + " ")
    open(base + ".asc", "w").write("\n".join(lines) + "\n")
    n = len(px)
    out = bytearray(struct.pack("<I", n))
    out += struct.pack(f"<{n}d", *px) + struct.pack(f"<{n}d", *py) + struct.pack(f"<{n}d", *([0.0] * n))
    out += struct.pack(f"<{n}i", *([-30] * n)) + struct.pack(f"<{n}d", *([0.0] * n))
    out += struct.pack("<I", len(elements)) + struct.pack(f"<{len(elements)}i", *[c for c, _ in elements])
    plist = [k for _, nodes in elements for k in nodes]
    out += struct.pack("<I", len(plist)) + struct.pack(f"<{len(plist)}I", *plist)
    out += struct.pack("<I", len(pfverts)) + struct.pack(f"<{len(pfverts)}i", *pfverts)
    out += struct.pack("<I", len(elements)) + struct.pack(f"<{len(elements)}i", *([0] * len(elements)))
    open(base + ".dat", "wb").write(out)
    with open(base + "-regions.txt", "w") as f:
        f.write("no properties\n\n" + "\n".join(fname for fname, _, _ in fams) + "\n")


if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    px, py, tri = build()
    cut = sum(1 for t in tri if min(px[list(t)]) < -1e-6 and max(px[list(t)]) > 1e-6)
    write(outdir + "/test-210-cylinder", px, py, families(px, py, tri))
    print("%d nodes, %d triangles, %d of them cut by the axis" % (len(px), len(tri), cut))
