#!/usr/bin/env python3
# SPDX-FileCopyrightText: (c) 2026 The CSMP++ project
# SPDX-License-Identifier: LGPL-3.0-only
"""test-211-obtuse: the test-211 box (600 m x 200 m, 30 x 20 rectangles of 20 m x 10 m, alternating diagonals) with the
nodes of the row just below the water table (y = -110 m; surface y = 0, base y = -200 m) shifted sideways by +-SHIFT m, alternately. The triangles that
span the water-table row then have angles above 90 degrees, and the edge opposite such an angle has a negative finite-
element coefficient -A gradN_a . gradN_b. The script prints the number of those edges that join a water-table node
(y = -100 m) to a node below it. Usage: python3 make_test211_obtuse_mesh.py [outdir]. Writes ANSYS .asc/.dat and
-regions.txt with the families HOST, BOTTOM, RIGHT, LEFT, TOP (as test-211-structured).
@author S. Scott (generated with Claude Code)
"""
import struct
import sys

import numpy as np

W, H, NX, NY = 600.0, 200.0, 30, 20
TITLE = ("%s.asc: box 600 m x 200 m, 30 x 20 rectangles (20 m x 10 m), alternating diagonals, row y = -110 m shifted "
         "sideways by +-8 m (obtuse triangles across the water table)")
WATER_TABLE, SHIFT = -100.0, 8.0   # water table 100 m above the base


def build():
    dx, dy = W / NX, H / NY
    px, py = [], []
    for i in range(NX + 1):
        for j in range(NY + 1):
            x, y = i * dx, j * dy - H
            if abs(y - (WATER_TABLE - dy)) < 1e-9 and 0 < i < NX:
                x += SHIFT if i % 2 else -SHIFT
            px.append(x)
            py.append(y)
    nid = lambda i, j: i * (NY + 1) + j
    tri = []
    for i in range(NX):
        for j in range(NY):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            tri += [(a, c, b), (a, d, c)] if (i + j) % 2 == 0 else [(a, d, b), (b, d, c)]
    px, py = np.array(px), np.array(py)
    out = []
    for a, b, c in tri:   # clockwise, as test-211-structured
        area = (px[b] - px[a]) * (py[c] - py[a]) - (px[c] - px[a]) * (py[b] - py[a])
        out.append((a, b, c) if area < 0 else (a, c, b))
    return px, py, out


def negative_edges(px, py, tri):
    found = []
    for t in tri:
        p = np.array([[px[k], py[k]] for k in t])
        for k in range(3):
            o, u, w = p[k], p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k]
            cot = np.dot(u, w) / abs(u[0] * w[1] - u[1] * w[0])
            if cot < -1e-12:
                found.append((t[(k + 1) % 3], t[(k + 2) % 3]))
    return found


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
        if abs(py[a] + H) < 1e-9 and abs(py[b] + H) < 1e-9: fam["BOTTOM"].append((a, b))
        elif abs(px[a] - W) < 1e-9 and abs(px[b] - W) < 1e-9: fam["RIGHT"].append((a, b))
        elif px[a] < 1e-9 and px[b] < 1e-9: fam["LEFT"].append((a, b))
        else: fam["TOP"].append((a, b))
    return [("HOST", "TRI_3", tri)] + [(k, "BAR_2", v) for k, v in fam.items()]


def write(base, px, py, fams):
    """ANSYS .asc/.dat/-regions.txt in the layout of the other test meshes."""
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
    neg = negative_edges(px, py, tri)
    spanning = [e for e in neg if (abs(py[e[0]] - WATER_TABLE) < 1e-9) != (abs(py[e[1]] - WATER_TABLE) < 1e-9)
                and max(py[e[0]], py[e[1]]) <= WATER_TABLE + 1e-9]
    write(outdir + "/test-211-obtuse", px, py, families(px, py, tri))
    print("%d nodes, %d triangles; %d negative edges, %d of them join a water-table node to a node below it"
          % (len(px), len(tri), len(neg), len(spanning)))
