#!/usr/bin/env python3
"""Prototype of the proposed solid_lineview interface: read the instanced .npz scene (g4poly/tgeopoly) directly.

    python3 lineview_npz.py scene.npz <any solid_lineview options>

Builds the same PrimScene that parse_prim() builds from a .prim (per-mesh weld/clean done once per unique mesh,
then instanced with the placement matrices) and runs solid_lineview's main() unchanged.
"""
import importlib.machinery
import importlib.util
import json
import sys
import time

import numpy as np

LV = '$HOME/Projects/solid/solid_dd4hep/bin/solid_lineview'
_spec = importlib.util.spec_from_loader('slv', importlib.machinery.SourceFileLoader('slv', LV))
slv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(slv)


def parse_npz(path):
    t0 = time.perf_counter()
    z = np.load(path)
    vo, fo = z['mesh_vertex_offset'], z['mesh_facet_offset']
    MV = z['mesh_vertices'].astype(np.float64)
    MF = z['mesh_facets']
    cache = {}

    def mesh(k):
        if k not in cache:
            V = MV[vo[k]:vo[k + 1]]
            F = np.abs(MF[fo[k]:fo[k + 1]]).astype(np.int64) - 1          # 0 (absent) -> -1
            Fvis = np.ones(F.shape, bool)                                   # DAWNFILE writes every edge visible
            V, F = slv._weld(V, F)
            F, Fvis = slv._clean_faces(F, Fvis)
            cache[k] = (V, F, Fvis)
        return cache[k]

    sc = slv.PrimScene()
    names, copy, depth = z['node_name'], z['node_copy'], z['node_depth']
    rgba, wire, M, nm = z['node_rgba'], z['node_wireframe'], z['node_matrix'], z['node_mesh']
    bV, bF, bE = [], [], []
    nv = nf = 0
    for i in np.nonzero(nm >= 0)[0]:
        V, F, Fvis = mesh(int(nm[i]))
        if len(F) == 0:
            continue
        R = M[i, :9].reshape(3, 3)
        W = V @ R.T + M[i, 9:]
        if np.linalg.det(R) < 0:                                            # reflection: keep normals outward
            F = slv._flip_faces(F)
        s = slv.Solid()
        s.index = len(sc.solids)
        s.name = '%s:%d' % (names[i], copy[i])
        s.depth = int(depth[i])
        s.color = tuple(float(x) for x in rgba[i][:3])
        s.wire = bool(wire[i])
        s.kind = 'Polyhedron'
        s.line = int(i)
        s.v0, s.v1, s.f0, s.f1 = nv, nv + len(W), nf, nf + len(F)
        s.bmin, s.bmax = W.min(0), W.max(0)
        bV.append(W)
        bF.append(np.where(F >= 0, F + nv, -1))
        bE.append(Fvis)
        nv += len(W)
        nf += len(F)
        sc.solids.append(s)
        sc.counts['Polyhedron'] += 1
    if bV:
        sc.V = np.vstack(bV)
        sc.F = np.vstack(bF)
        sc.Fvis = np.vstack(bE)
        sc.fsolid = np.repeat(np.arange(len(sc.solids), dtype=np.int32), [s.f1 - s.f0 for s in sc.solids])
        lo, hi = sc.V.min(0), sc.V.max(0)
        c, r = 0.5 * (lo + hi), 0.5 * np.linalg.norm(hi - lo)
        sc.bbox = np.concatenate([c - r, c + r])
    slv._build_hierarchy(sc)
    meta = json.loads(str(z['meta'])) if 'meta' in z.files else {}
    print('[lineview_npz] %s: %d solids from %d unique meshes in %.3f s (source: %s)' % (
        path, len(sc.solids), len(cache), time.perf_counter() - t0, meta.get('source', '?')), file=sys.stderr)
    return sc


_orig = slv.parse_prim
slv.parse_prim = lambda p: parse_npz(p) if p.endswith('.npz') else _orig(p)

if __name__ == '__main__':
    sys.exit(slv.main(sys.argv[1:]))
