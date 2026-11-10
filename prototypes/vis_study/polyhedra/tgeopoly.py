#!/usr/bin/env python3
"""Prototype (feasibility study): tessellate the DD4hep TGeo geometry with ROOT only (no Geant4 at all).

    python3 tgeopoly.py -c $DETECTOR_PATH/$DETECTOR_CONFIG.xml -o out/tgeo_cfg [--nseg 72]

TGeoShape::MakeBuffer3D -> RootCsg meshes (composites via RootCsg::BuildUnion/Difference/Intersection), DD4hep
VisAttr colours/visibility/showDaughters, global TGeoHMatrix.  Writes the same .prim / .npz as g4poly.py.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g4poly  # noqa: E402  (writers, IdEncoder, detelement_index; no Geant4 is loaded by this import)

log = g4poly.log


def cling():
    t0 = time.time()
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    import dd4hep  # noqa: F401  (loads DDCore + dictionaries)
    ROOT.gSystem.Load('libRCsg')
    if not ROOT.gInterpreter.Declare('#include "%s"' % os.path.join(HERE, 'solidvis_tgeo.h')):
        raise RuntimeError('cannot declare solidvis_tgeo.h')
    log('cling setup (DDCore + RootCsg + helper) %.2f s' % (time.time() - t0))
    return ROOT


# ------------------------------------------------------------------------------------------- polygon handling
def _newell(P):
    n = np.zeros(3)
    for a, b in zip(P, np.roll(P, -1, 0)):
        n += np.array([(a[1] - b[1]) * (a[2] + b[2]), (a[2] - b[2]) * (a[0] + b[0]), (a[0] - b[0]) * (a[1] + b[1])])
    return n


def _ear_clip(P2, idx):
    """Triangulate a simple 2-D polygon (counter-clockwise) given as P2 (n,2); returns index triples."""
    ids = list(range(len(P2)))
    tris = []

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    guard = 0
    while len(ids) > 3 and guard < 10 * len(P2) ** 2:
        guard += 1
        m = len(ids)
        found = False
        for k in range(m):
            i0, i1, i2 = ids[(k - 1) % m], ids[k], ids[(k + 1) % m]
            a, b, c = P2[i0], P2[i1], P2[i2]
            if cross(a, b, c) <= 1e-14:
                continue
            ok = True
            for j in ids:
                if j in (i0, i1, i2):
                    continue
                p = P2[j]
                if cross(a, b, p) >= -1e-12 and cross(b, c, p) >= -1e-12 and cross(c, a, p) >= -1e-12:
                    ok = False
                    break
            if ok:
                tris.append((idx[i0], idx[i1], idx[i2]))
                del ids[k]
                found = True
                break
        if not found:          # degenerate / self-touching: fall back to a fan of the rest
            for k in range(1, len(ids) - 1):
                tris.append((idx[ids[0]], idx[ids[k]], idx[ids[k + 1]]))
            return tris, False
    if len(ids) == 3:
        tris.append((idx[ids[0]], idx[ids[1]], idx[ids[2]]))
    return tris, True


def polys_to_faces(V, pol, stats):
    """RootCsg polygon list -> (Nf,4) 1-based facets (0 = no 4th vertex); n-gons triangulated properly."""
    F = []
    i = 0
    while i < len(pol):
        n = int(pol[i])
        idx = [int(x) for x in pol[i + 1:i + 1 + n]]
        i += 1 + n
        # drop repeated / coincident consecutive vertices
        clean = []
        for k in idx:
            if not clean or (k != clean[-1] and np.linalg.norm(V[k] - V[clean[-1]]) > 1e-9):
                clean.append(k)
        while len(clean) > 1 and (clean[0] == clean[-1] or np.linalg.norm(V[clean[0]] - V[clean[-1]]) <= 1e-9):
            clean.pop()
        n = len(clean)
        if n < 3:
            stats['degenerate_polygons'] += 1
            continue
        if n == 3:
            F.append([clean[0] + 1, clean[1] + 1, clean[2] + 1, 0])
            continue
        P = V[clean]
        nrm = _newell(P)
        ln = np.linalg.norm(nrm)
        if ln < 1e-14:
            stats['degenerate_polygons'] += 1
            continue
        nrm /= ln
        e1 = np.roll(P, -1, 0) - P
        e0 = P - np.roll(P, 1, 0)
        turn = np.einsum('ij,j->i', np.cross(e0, e1), nrm)
        convex = np.all(turn >= -1e-9 * max(1.0, np.abs(P).max()) ** 2)
        if n == 4 and convex:
            F.append([clean[0] + 1, clean[1] + 1, clean[2] + 1, clean[3] + 1])
            continue
        stats['ngons'] += 1
        if convex:
            for k in range(1, n - 1):
                F.append([clean[0] + 1, clean[k] + 1, clean[k + 1] + 1, 0])
            continue
        stats['concave_ngons'] += 1
        # project to the polygon plane, counter-clockwise w.r.t. nrm
        a = np.cross(nrm, [1.0, 0, 0] if abs(nrm[0]) < 0.9 else [0, 1.0, 0])
        a /= np.linalg.norm(a)
        b = np.cross(nrm, a)
        P2 = np.stack([P @ a, P @ b], 1)
        tris, ok = _ear_clip(P2, clean)
        if not ok:
            stats['earclip_fallback'] += 1
        for t in tris:
            F.append([t[0] + 1, t[1] + 1, t[2] + 1, 0])
    return np.array(F, np.int32).reshape(-1, 4)


def mesh_stats(V, F):
    """Signed volume and edge topology (boundary edges = cracks / T-junctions, non-manifold, orientation)."""
    q = F[:, 3] > 0
    tri = np.vstack([F[:, :3], F[q][:, [0, 2, 3]]]) - 1
    a, b, c = V[tri[:, 0]], V[tri[:, 1]], V[tri[:, 2]]
    vol = np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0
    nxt = np.where(q[:, None], F[:, [1, 2, 3, 0]], F[:, [1, 2, 0, 0]])
    valid = np.ones(F.shape, bool)
    valid[:, 3] = q
    ea, eb = F[valid] - 1, nxt[valid] - 1
    lo, hi = np.minimum(ea, eb), np.maximum(ea, eb)
    key = lo.astype(np.int64) * (len(V) + 1) + hi
    u, cnt = np.unique(key, return_counts=True)
    return dict(volume=float(vol), edges=int(len(u)), boundary_edges=int((cnt == 1).sum()),
                nonmanifold_edges=int((cnt > 2).sum()))


def _loops(F):
    return [[int(x) - 1 for x in f if x > 0] for f in F]


def _boundary_edges(loops):
    cnt = {}
    for fi, L in enumerate(loops):
        for k in range(len(L)):
            a, b = L[k], L[(k + 1) % len(L)]
            key = (a, b) if a < b else (b, a)
            cnt[key] = cnt.get(key, 0) + 1
    return [e for e, c in cnt.items() if c == 1]


def fix_tjunctions(V, F, tol_rel=1e-7, max_iter=6):
    """Make CSG output watertight: insert vertices that lie on boundary edges (T-junctions) into the faces
    owning those edges, re-triangulate touched faces as a centroid fan.  Returns V, F, remaining boundary edges."""
    V = np.asarray(V, float)
    loops = _loops(F)
    scale = float(np.ptp(V, 0).max()) if len(V) else 1.0
    tol = tol_rel * scale + 1e-9
    for _ in range(max_iter):
        bnd = _boundary_edges(loops)
        if not bnd:
            break
        E = np.array(bnd)
        bv = np.unique(E)
        P = V[bv]
        A, B = V[E[:, 0]], V[E[:, 1]]
        d = B - A
        L2 = np.maximum((d * d).sum(1), 1e-300)
        ins = {}
        step = max(1, 2_000_000 // max(1, len(bv)))
        for c0 in range(0, len(E), step):
            sl = slice(c0, c0 + step)
            t = np.einsum('ecj,ej->ec', P[None, :, :] - A[sl, None, :], d[sl]) / L2[sl, None]
            Q = A[sl, None, :] + t[..., None] * d[sl, None, :]
            dist = np.linalg.norm(P[None, :, :] - Q, axis=2)
            hit = (t > 1e-9) & (t < 1 - 1e-9) & (dist < tol)
            for ei, ci in zip(*np.nonzero(hit)):
                e = tuple(E[c0 + ei])
                ins.setdefault(e, []).append((float(t[ei, ci]), int(bv[ci])))
        if not ins:
            break
        Vn = [V]
        nv = len(V)
        new = []
        for L in loops:
            out = []
            touched = False
            for k in range(len(L)):
                a, b = L[k], L[(k + 1) % len(L)]
                out.append(a)
                key = (a, b) if a < b else (b, a)
                if key in ins:
                    # parameters were computed along min->max; walk them along a->b
                    pts = sorted(ins[key]) if a < b else sorted((1 - t, v) for t, v in ins[key])
                    out.extend(v for _, v in pts if v != a and v != b)
                    touched = True
            if not touched:
                new.append(L)
                continue
            c = V[out].mean(0)
            Vn.append(c[None])
            ci = nv
            nv += 1
            for k in range(len(out)):
                new.append([ci, out[k], out[(k + 1) % len(out)]])
        V = np.vstack(Vn)
        loops = new
    Fo = np.zeros((len(loops), 4), np.int32)
    for i, L in enumerate(loops):
        Fo[i, :len(L)] = np.array(L) + 1
    return V, Fo, len(_boundary_edges(loops))


def weld(V, F, tol=1e-7):
    """Merge coincident vertices (RootCsg output repeats vertices per fragment)."""
    if len(V) == 0:
        return V, F
    key = np.round(V / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    V2 = V[first]
    F2 = F.copy()
    m = F2 > 0
    F2[m] = inv[F2[m] - 1] + 1
    return V2, F2


# ------------------------------------------------------------------------------------------------ collect
def collect(ROOT, description, nseg=72, with_ids=True, do_weld=True, fix_tj=True):
    timings = {}
    t0 = time.time()
    w = ROOT.solidtgeo.Walker()
    w.walk(ROOT.gGeoManager.GetTopNode())
    nodes = w.nodes
    n = nodes.size()
    timings['walk'] = time.time() - t0

    t0 = time.time()
    names, copy = [], np.zeros(n, np.int32)
    depth, g4depth, parent = np.zeros(n, np.int32), np.zeros(n, np.int32), np.zeros(n, np.int32)
    assembly, vis, wire = np.zeros(n, bool), np.zeros(n, bool), np.zeros(n, bool)
    rgba, mat = np.zeros((n, 4), np.float32), np.zeros((n, 12))
    shape_ptr = np.zeros(n, np.int64)
    shapes = {}
    tnodes = []
    for i in range(n):
        nd = nodes[i]
        tn = nd.node
        tnodes.append(tn)
        names.append(str(tn.GetName()))
        copy[i] = tn.GetNumber()
        depth[i], g4depth[i], parent[i] = nd.depth, nd.g4depth, nd.parent
        assembly[i], vis[i], wire[i] = nd.assembly, nd.visible, nd.wire
        rgba[i] = list(nd.rgba)
        mat[i] = list(nd.m)
        a = ROOT.addressof(nd.shape)
        shape_ptr[i] = a
        shapes.setdefault(a, nd.shape)
    mat[:, 9:] *= 10.0                      # cm -> mm
    timings['nodes_to_numpy'] = time.time() - t0

    t0 = time.time()
    drawn = vis & ~assembly
    mesh_of_shape, mv, mf, mnames, mtypes, mfail, mstat = {}, [], [], [], [], [], []
    from collections import Counter
    pstats = Counter()
    cpp = 0.0
    flipped = 0
    for a in np.unique(shape_ptr[drawn]):
        s = shapes[int(a)]
        m = ROOT.solidtgeo.tessellate(s, nseg)
        cpp += m.seconds
        if not m.ok:
            mfail.append((str(s.GetName()), str(s.ClassName()), str(m.err)))
            mesh_of_shape[int(a)] = -1
            continue
        V = np.frombuffer(m.v.data(), dtype=np.float64, count=m.v.size()).reshape(-1, 3) * 10.0
        pol = np.frombuffer(m.pol.data(), dtype=np.int32, count=m.pol.size()).copy()
        F = polys_to_faces(V, pol, pstats)
        if do_weld:
            V, F = weld(V, F)
        nb_before = None
        if fix_tj and s.IsComposite():
            nb_before = mesh_stats(V, F)['boundary_edges']
            V, F, _ = fix_tjunctions(V, F)
        st = mesh_stats(V, F)
        if nb_before is not None:
            st['boundary_edges_before_fix'] = nb_before
        if st['volume'] < 0:                 # inside-out: flip the winding
            flipped += 1
            G = F.copy()
            q = F[:, 3] > 0
            G[q] = F[q][:, ::-1]
            G[~q, :3] = F[~q][:, 2::-1]
            F = G
            st['volume'] = -st['volume']
            st['flipped'] = True
        try:
            st['capacity'] = float(s.Capacity()) * 1000.0      # cm3 -> mm3 (exact for primitives, MC for booleans)
        except Exception:
            pass
        st['composite'] = bool(s.IsComposite())
        st['ms'] = m.seconds * 1e3
        mesh_of_shape[int(a)] = len(mv)
        mv.append(V)
        mf.append(F)
        mnames.append(str(s.GetName()))
        mtypes.append(str(s.ClassName()))
        mstat.append(st)
    node_mesh = np.array([mesh_of_shape.get(int(a), -1) if drawn[i] else -1 for i, a in enumerate(shape_ptr)],
                         np.int32)
    timings['tessellate_total'] = time.time() - t0
    timings['tessellate_cpp'] = cpp
    out = dict(names=names, depth=g4depth, tgeo_depth=depth, parent=parent, copy=copy, kind=np.zeros(n, np.int8),
               visible=vis, has_vis=np.ones(n, bool), wire=wire, rgba=rgba, matrix=mat, density=np.zeros(n, np.float32),
               node_mesh=node_mesh, assembly=assembly, mesh_v=mv, mesh_f=mf, mesh_names=mnames, mesh_types=mtypes,
               mesh_fail=mfail, mesh_stats=mstat, polygon_stats=dict(pstats), flipped=flipped, timings=timings)

    if with_ids:
        t0 = time.time()
        deidx = g4poly.detelement_index(ROOT, description)
        enc = g4poly.IdEncoder(ROOT, description)
        det_path, subdet, readout, tgeo_path = [''] * n, [''] * n, [''] * n, [''] * n
        sensitive = np.zeros(n, bool)
        volid, volmask = np.zeros(n, np.uint64), np.zeros(n, np.uint64)
        chains = {0: []}
        for i in range(1, n):
            chain = chains[int(parent[i])] + [tnodes[i]]
            chains[i] = chain
            tgeo_path[i] = '/'.join(names[j] for j in _path_idx(parent, i))
            sd_owner = ''
            for tn in reversed(chain):
                hit = deidx.get(ROOT.addressof(tn))
                if hit:
                    det_path[i], sd_owner = hit
                    break
            subdet[i] = sd_owner
            if not sd_owner:
                continue
            name, vid, mask, _ = enc.encode(sd_owner, chain)
            readout[i], volid[i], volmask[i] = name, vid, mask
            try:
                sensitive[i] = bool(ROOT.dd4hep.Volume(tnodes[i].GetVolume()).isSensitive())
            except Exception:
                pass
        out.update(det_path=det_path, subdet=subdet, readout=readout, tgeo_path=tgeo_path, sensitive=sensitive,
                   volid=volid, volmask=volmask,
                   readouts={k: dict(name=v[0], fields=v[2]) for k, v in enc.ro.items() if v[0]})
        timings['ids'] = time.time() - t0
    return out


def _path_idx(parent, i):
    p = []
    while i > 0:
        p.append(i)
        i = int(parent[i])
    return p[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--compact', required=True)
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--nseg', type=int, default=72)
    ap.add_argument('--formats', default='prim,npz')
    ap.add_argument('--bbox-from', default=None)
    ap.add_argument('--no-weld', action='store_true')
    ap.add_argument('--no-fix', action='store_true', help='keep RootCsg T-junctions')
    ap.add_argument('--honour-drawing-style', action='store_true',
                    help='VisAttr drawingStyle=wireframe -> /ForceWireframe 1 (DDG4 drops it)')
    ap.add_argument('--raw-tgeo', action='store_true',
                    help='plain gGeoManager->SetNsegments semantics (no phi scaling, seam and sphere fixes)')
    a = ap.parse_args()
    T0 = time.time()
    ROOT = cling()
    t0 = time.time()
    desc = ROOT.dd4hep.Detector.getInstance()
    desc.fromXML(a.compact)
    log('DD4hep load %.2f s' % (time.time() - t0))
    t0 = time.time()
    ROOT.solidtgeo.options().honour_drawing_style = a.honour_drawing_style
    if a.raw_tgeo:
        o = ROOT.solidtgeo.options()
        o.scale_phi = False
        o.full_seg_as_tube = False
        o.sphere_divisions = False
    d = collect(ROOT, desc, nseg=a.nseg, do_weld=not a.no_weld, fix_tj=not a.no_fix)
    log('collect %.2f s' % (time.time() - t0))
    if d['mesh_fail']:
        log('meshes failed for %d shapes:' % len(d['mesh_fail']), d['mesh_fail'][:10])
    bbox = None
    if a.bbox_from:
        with open(a.bbox_from) as fh:
            for line in fh:
                if line.startswith('/BoundingBox'):
                    bbox = [float(x) for x in line.split()[1:7]]
                    break
    fm = a.formats.split(',')
    if 'prim' in fm:
        g4poly.write_prim(d, a.out + '.prim', bbox=bbox)
    if 'npz' in fm:
        g4poly.write_npz(d, a.out + '.npz', dict(source='TGeo MakeBuffer3D + RootCsg', nseg=a.nseg,
                                                  readouts=d.get('readouts', {})))
    s = g4poly.summary(d)
    ms = d['mesh_stats']
    s.update(polygon_stats=d['polygon_stats'], flipped=d['flipped'],
             meshes_with_boundary_edges=sum(1 for x in ms if x['boundary_edges']),
             composite_meshes=sum(1 for x in ms if x['composite']),
             composite_with_boundary_edges=sum(1 for x in ms if x['composite'] and x['boundary_edges']),
             boundary_edges_total=sum(x['boundary_edges'] for x in ms),
             boundary_edges_before_fix=sum(x.get('boundary_edges_before_fix', 0) for x in ms),
             nonmanifold_edges_total=sum(x['nonmanifold_edges'] for x in ms))
    log(json.dumps(s))
    with open(a.out + '.summary.json', 'w') as fh:
        json.dump(dict(s, mesh_fail=d['mesh_fail'],
                       mesh_stats=[dict(st, name=nm, type=tp) for st, nm, tp in
                                   zip(ms, d['mesh_names'], d['mesh_types'])]), fh, indent=1)
    log('total wall %.2f s' % (time.time() - T0))


if __name__ == '__main__':
    main()
