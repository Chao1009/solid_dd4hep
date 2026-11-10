#!/usr/bin/env python3
"""Prototype (feasibility study): dump the Geant4 tessellation of a DD4hep geometry WITHOUT a vis driver.

Two ways to run it (inside eic_xl, after `source thissolid.sh <config>`):

  standalone (no run manager, no physics):   python3 g4poly.py -c $DETECTOR_PATH/$DETECTOR_CONFIG.xml -o out/geo
  inside ddsim (same session as the events): ddsim --steeringFile polydump_steer.py ...   (see that file)

Writes <out>.prim (DAWNFILE g4.prim 2.4 subset that bin/solid_lineview reads unchanged) and/or
<out>.npz (unique meshes + placements, with DD4hep DetElement path, readout and volume ID per placement).
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = None


def log(*a):
    print('[g4poly]', *a, flush=True)


def cling():
    """Load DDG4 + Geant4 headers and the helper into cling (once)."""
    global _ROOT
    if _ROOT is not None:
        return _ROOT
    t0 = time.time()
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gSystem.Load('libDDG4')
    try:
        flags = subprocess.check_output(['geant4-config', '--cflags'], text=True).split()
        incs = [f[2:] for f in flags if f.startswith('-I') and ('geant4' in f or 'clhep' in f)]
    except Exception:
        incs = ['/opt/local/include/Geant4']
    for inc in incs:
        ROOT.gInterpreter.AddIncludePath(inc)
    if not ROOT.gInterpreter.Declare('#include "%s"' % os.path.join(HERE, 'solidvis_g4.h')):
        raise RuntimeError('cannot declare solidvis_g4.h')
    log('cling setup (DDG4 + Geant4 headers + helper) %.2f s' % (time.time() - t0))
    _ROOT = ROOT
    return ROOT


# ----------------------------------------------------------------------------------------------------
def detelement_index(ROOT, description):
    """TGeoNode address -> (DetElement path, readout owner).  The readout owner is the nearest ancestor-or-self
    DetElement that has a SensitiveDetector of its name (e.g. GEMTrackerLA inside GEMTrackerLiteAssembly), else
    the top-level subdetector name."""
    idx = {}
    sdnames = set()
    try:
        for p in description.sensitiveDetectors():
            sdnames.add(str(p.first))
    except Exception:
        pass
    world = description.world()

    def rec(de, top, owner):
        nm = str(de.name())
        if nm in sdnames:
            owner = nm
        try:
            pl = de.placement()
            if pl.isValid():
                idx[ROOT.addressof(pl.ptr())] = (str(de.path()), owner or top)
        except Exception:
            pass
        for name, child in de.children():
            rec(child, top, owner)

    for name, child in world.children():
        rec(child, str(name), '')
    return idx


class IdEncoder:
    """Partial DD4hep volume IDs: fields set by the placements along a TGeo node chain, encoded with the
    readout of the subdetector, plus the mask of those fields (hit cellID & mask == volid -> this volume)."""

    def __init__(self, ROOT, description):
        self.ROOT = ROOT
        self.desc = description
        self.ro = {}          # subdet -> (readout name, {field: (offset, width)}, description string)
        self.pv_ids = {}      # TGeoNode address -> [(field, value)]

    def readout(self, subdet):
        if subdet not in self.ro:
            res = ('', {}, '')
            try:
                sd = self.desc.sensitiveDetector(subdet)
                if sd.isValid():
                    ro = sd.readout() if callable(sd.readout) else sd.readout
                    if ro.isValid():
                        spec = ro.idSpec()
                        fields = {}
                        for p in spec.fields():
                            fields[str(p.first)] = (int(p.second.offset()), int(p.second.width()))
                        res = (str(ro.name()), fields, str(spec.fieldDescription()))
            except Exception as exc:          # noqa: BLE001  (prototype)
                log('readout lookup failed for', subdet, exc)
            self.ro[subdet] = res
        return self.ro[subdet]

    def ids_of(self, tnode):
        a = self.ROOT.addressof(tnode)
        if a not in self.pv_ids:
            try:
                self.pv_ids[a] = [(str(p.first), int(p.second)) for p in self.ROOT.dd4hep.PlacedVolume(tnode).volIDs()]
            except Exception:
                self.pv_ids[a] = []
        return self.pv_ids[a]

    def encode(self, subdet, chain):
        name, fields, _ = self.readout(subdet)
        ids = {}
        for tn in chain:
            for f, v in self.ids_of(tn):
                ids[f] = v
        vid = mask = 0
        for f, v in ids.items():
            if f in fields:
                off, w = fields[f]
                m = (1 << w) - 1
                vid |= (v & m) << off
                mask |= m << off
        return name, vid, mask, ids


def collect(ROOT, world_pv, nseg=72, description=None, with_ids=True):
    """Walk the G4 tree; return a dict with unique meshes and placements (numpy arrays)."""
    timings = {}
    t0 = time.time()
    w = ROOT.solidvis.Walker()
    w.walk(world_pv)
    nodes = w.nodes
    n = nodes.size()
    timings['walk'] = time.time() - t0

    # --- per-node scalars (Python loop over the C++ vector; fine for 1e4 nodes)
    t0 = time.time()
    names, depth, parent, copy, kind = [], np.zeros(n, np.int32), np.zeros(n, np.int32), np.zeros(n, np.int32), \
        np.zeros(n, np.int8)
    vis = np.zeros(n, bool)
    hasvis = np.zeros(n, bool)
    dinv = np.zeros(n, bool)
    wire = np.zeros(n, bool)
    rgba = np.zeros((n, 4), np.float32)
    mat = np.zeros((n, 12))
    dens = np.zeros(n, np.float32)
    solid_ptr = np.zeros(n, np.int64)
    pv_ptr = np.zeros(n, np.int64)
    solids = {}
    for i in range(n):
        nd = nodes[i]
        names.append(str(nd.pv.GetName()))
        depth[i], parent[i], copy[i], kind[i] = nd.depth, nd.parent, nd.copy, nd.kind
        vis[i], hasvis[i], dinv[i], wire[i] = nd.visible, nd.has_vis, nd.daughters_invisible, nd.wire
        rgba[i] = list(nd.rgba)
        mat[i] = list(nd.m)
        dens[i] = nd.density
        a = ROOT.addressof(nd.solid)
        solid_ptr[i] = a
        pv_ptr[i] = ROOT.addressof(nd.pv)
        if a not in solids:
            solids[a] = nd.solid
    timings['nodes_to_numpy'] = time.time() - t0

    # --- unique meshes (only for solids of drawn = visible nodes, like the vis drivers)
    t0 = time.time()
    drawn = vis.copy()
    mesh_of_solid = {}
    mv, mf, mnames, mtypes, mfail = [], [], [], [], []
    cpp_time = 0.0
    for a in np.unique(solid_ptr[drawn]):
        s = solids[int(a)]
        m = ROOT.solidvis.tessellate(s, nseg)
        cpp_time += m.seconds
        if not m.ok:
            mfail.append((str(s.GetName()), str(s.GetEntityType())))
            mesh_of_solid[int(a)] = -1
            continue
        V = np.frombuffer(m.v.data(), dtype=np.float64, count=m.v.size()).reshape(-1, 3).copy()
        F = np.frombuffer(m.f.data(), dtype=np.int32, count=m.f.size()).reshape(-1, 4).copy()
        mesh_of_solid[int(a)] = len(mv)
        mv.append(V)
        mf.append(F)
        mnames.append(str(s.GetName()))
        mtypes.append(str(s.GetEntityType()))
    node_mesh = np.array([mesh_of_solid.get(int(a), -1) if drawn[i] else -1 for i, a in enumerate(solid_ptr)],
                         np.int32)
    timings['tessellate_total'] = time.time() - t0
    timings['tessellate_cpp'] = cpp_time

    out = dict(names=names, depth=depth, parent=parent, copy=copy, kind=kind, visible=vis, has_vis=hasvis,
               daughters_invisible=dinv, wire=wire, rgba=rgba, matrix=mat, density=dens, node_mesh=node_mesh,
               mesh_v=mv, mesh_f=mf, mesh_names=mnames, mesh_types=mtypes, mesh_fail=mfail, timings=timings)

    # --- DD4hep identification: TGeo node chain, DetElement path, readout, (partial) volume ID + mask
    if with_ids and description is not None:
        t0 = time.time()
        pvmap = ROOT.solidvis.PVMap()
        pvmap.build()
        place = {ROOT.addressof(p.first): p.second for p in pvmap.place}
        imprint = {ROOT.addressof(p.first): list(p.second) for p in pvmap.imprint}
        deidx = detelement_index(ROOT, description)
        enc = IdEncoder(ROOT, description)
        det_path = [''] * n
        subdet = [''] * n
        readout = [''] * n
        tgeo_path = [''] * n
        sensitive = np.zeros(n, bool)
        volid = np.zeros(n, np.uint64)
        volmask = np.zeros(n, np.uint64)
        g4volid = np.zeros(n, np.int64)
        vec_t = ROOT.std.vector['const G4VPhysicalVolume*']
        chains = {0: []}
        for i in range(n):
            if i == 0:
                continue
            a = int(pv_ptr[i])
            own = [place[a]] if a in place else imprint.get(a, [])
            chain = chains[int(parent[i])] + own
            chains[i] = chain
            tgeo_path[i] = '/'.join(str(tn.GetName()) for tn in chain)
            top = ''
            for tn in reversed(chain):
                hit = deidx.get(ROOT.addressof(tn))
                if hit:
                    det_path[i], top = hit
                    break
            subdet[i] = top
            if not chain or not top:
                continue
            name, vid, mask, ids = enc.encode(top, chain)
            readout[i], volid[i], volmask[i] = name, vid, mask
            try:
                sensitive[i] = bool(ROOT.dd4hep.PlacedVolume(chain[-1]).volume().isSensitive())
            except Exception:
                pass
            if sensitive[i]:
                path = vec_t()
                j = i
                while j > 0:               # deepest first, world excluded (DDG4 touchable convention)
                    path.push_back(nodes[j].pv)
                    j = int(parent[j])
                g4volid[i] = ROOT.solidvis.volume_id(path)
        # cross-check with what a sensitive detector gets during tracking (navigator touchable at the centroid)
        if os.environ.get('POLYDUMP_CHECK_IDS', '0') == '1':
            nchk = nok = nother = 0
            pvout = ROOT.std.vector['const G4VPhysicalVolume*'](1)
            for i in np.nonzero(sensitive)[0]:
                k = mesh_of_solid.get(int(solid_ptr[i]), -1)
                if k < 0:
                    continue
                c = mv[k].mean(0)
                M = mat[i]
                wc = M[:9].reshape(3, 3) @ c + M[9:]
                ptr = ROOT.std.vector['const G4VPhysicalVolume*'](1)
                tid = ROOT.solidvis.volume_id_at(world_pv, float(wc[0]), float(wc[1]), float(wc[2]), ptr.data())
                if ROOT.addressof(ptr[0]) != int(pv_ptr[i]):
                    nother += 1
                    continue
                nchk += 1
                nok += int((tid & ((1 << 64) - 1)) == int(volid[i]))
            timings['touchable_check'] = dict(checked=nchk, equal=nok, point_in_other_volume=nother)
            log('touchable volume-ID check: %d checked, %d equal to the DD4hep-side encoding, %d skipped' % (
                nchk, nok, nother))
        out.update(det_path=det_path, subdet=subdet, readout=readout, tgeo_path=tgeo_path, sensitive=sensitive,
                   volid=volid, volmask=volmask, g4volid=g4volid,
                   readouts={k: dict(name=v[0], fields=v[2]) for k, v in enc.ro.items() if v[0]})
        timings['ids'] = time.time() - t0
    return out


# ----------------------------------------------------------------------------------------------------
def write_prim(d, path, bbox=None, edge_flags='dawn'):
    """DAWNFILE-compatible .prim (only what solid_lineview needs): one /Polyhedron per drawn placement.

    edge_flags='dawn': all labels positive, as the Geant4 11.4 DAWNFILE driver writes them;
    'g4': keep the G4Polyhedron edge flags (negative label = invisible edge)."""
    t0 = time.time()
    lines = ['##G4.PRIM-FORMAT-2.4']
    if bbox is None:
        lo = np.full(3, np.inf)
        hi = -lo
        for i in np.nonzero(d['node_mesh'] >= 0)[0]:
            V = d['mesh_v'][d['node_mesh'][i]]
            M = d['matrix'][i]
            W = V @ M[:9].reshape(3, 3).T + M[9:]
            lo = np.minimum(lo, W.min(0))
            hi = np.maximum(hi, W.max(0))
        c, r = 0.5 * (lo + hi), 0.5 * np.linalg.norm(hi - lo)
        bbox = np.concatenate([c - r, c + r])
    lines.append('/BoundingBox ' + ' '.join('%.9g' % x for x in bbox))
    lines += ['!SetCamera', '!OpenDevice', '!BeginModeling']
    vtxt = {}
    ftxt = {}
    nw = 0
    for i in range(len(d['names'])):
        k = d['node_mesh'][i]
        if k < 0:
            continue
        M = d['matrix'][i]
        R = M[:9].reshape(3, 3)
        lines.append('#--------------------')
        lines.append('#/PVName' + ' ' * (5 + 5 * int(d['depth'][i])) + '%s:%d' % (d['names'][i], d['copy'][i]))
        c = d['rgba'][i]
        lines.append('/ColorRGB %.9g %.9g %.9g' % (c[0], c[1], c[2]))
        lines.append('/ForceWireframe  %d' % int(d['wire'][i]))
        if np.linalg.det(R) > 0:
            lines.append('/Origin %.9g %.9g %.9g' % tuple(M[9:]))
            # /BaseVector = images of the local x and y axes = first two columns of R
            lines.append('/BaseVector %.9g %.9g %.9g %.9g %.9g %.9g' % (R[0, 0], R[1, 0], R[2, 0],
                                                                        R[0, 1], R[1, 1], R[2, 1]))
            if k not in vtxt:
                vtxt[k] = '\n'.join('/Vertex %.9g %.9g %.9g' % tuple(v) for v in d['mesh_v'][k])
                if edge_flags == 'g4':
                    ftxt[k] = '\n'.join('/Facet ' + ' '.join(str(int(x)) for x in f if x != 0) for f in d['mesh_f'][k])
                else:
                    ftxt[k] = '\n'.join('/Facet ' + ' '.join(str(abs(int(x))) for x in f if x != 0)
                                         for f in d['mesh_f'][k])
            lines += ['/Polyhedron', vtxt[k], ftxt[k], '/EndPolyhedron']
        else:
            # reflection: DAWN /BaseVector cannot express it, write world coordinates and reverse the winding
            nw += 1
            W = d['mesh_v'][k] @ R.T + M[9:]
            lines.append('/Origin 0 0 0')
            lines.append('/BaseVector 1 0 0 0 1 0')
            lines.append('/Polyhedron')
            lines.append('\n'.join('/Vertex %.9g %.9g %.9g' % tuple(v) for v in W))
            for f in d['mesh_f'][k]:
                f = [int(x) for x in f if x != 0]
                # reverse vertex order; edge flag of edge (a->b) is carried by label a in G4 GetFacet
                idx = [abs(x) for x in f][::-1]
                flg = [x > 0 for x in f]
                m = len(f)
                # original edge e_j = (v_j -> v_{j+1}) with flag flg[j]; reversed polygon edge
                # (idx[t] -> idx[t+1]) = original e_{m-2-t}
                lab = [idx[t] if (flg[(m - 2 - t) % m] or edge_flags != 'g4') else -idx[t] for t in range(m)]
                lines.append('/Facet ' + ' '.join(map(str, lab)))
            lines.append('/EndPolyhedron')
    lines += ['!EndModeling', '!DrawAll', '!CloseDevice', '']
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines))
    log('wrote %s (%.1f MB, %d reflected placements baked) in %.2f s' % (
        path, os.path.getsize(path) / 1e6, nw, time.time() - t0))


def write_npz(d, path, meta):
    t0 = time.time()
    mv = d['mesh_v']
    mf = d['mesh_f']
    voff = np.cumsum([0] + [len(v) for v in mv]).astype(np.int64)
    foff = np.cumsum([0] + [len(f) for f in mf]).astype(np.int64)
    arrays = dict(
        mesh_vertices=(np.vstack(mv) if mv else np.zeros((0, 3))).astype(np.float32),
        mesh_vertex_offset=voff,
        mesh_facets=(np.vstack(mf) if mf else np.zeros((0, 4))).astype(np.int32),
        mesh_facet_offset=foff,
        mesh_name=np.array(d['mesh_names']), mesh_type=np.array(d['mesh_types']),
        node_name=np.array(d['names']), node_copy=d['copy'], node_depth=d['depth'], node_parent=d['parent'],
        node_mesh=d['node_mesh'], node_matrix=d['matrix'], node_rgba=d['rgba'],
        node_visible=d['visible'], node_wireframe=d['wire'], node_density=d['density'],
        meta=np.array(json.dumps(meta)))
    for k in ('det_path', 'subdet', 'readout', 'tgeo_path'):
        if k in d:
            arrays['node_' + k] = np.array(d[k])
    if 'volid' in d:
        arrays['node_volume_id'] = d['volid']
        arrays['node_volume_mask'] = d['volmask']
        arrays['node_sensitive'] = d['sensitive']
        if 'g4volid' in d:
            arrays['node_g4_volume_id'] = d['g4volid']
    np.savez_compressed(path, **arrays)
    log('wrote %s (%.2f MB) in %.2f s' % (path, os.path.getsize(path) / 1e6, time.time() - t0))


def summary(d):
    nm = len(d['mesh_v'])
    drawn = d['node_mesh'] >= 0
    nf_inst = sum(len(d['mesh_f'][k]) for k in d['node_mesh'][drawn])
    nv_inst = sum(len(d['mesh_v'][k]) for k in d['node_mesh'][drawn])
    s = dict(nodes=len(d['names']), drawn=int(drawn.sum()), invisible=int((~d['visible']).sum()),
             unique_meshes=nm, failed=len(d['mesh_fail']), vertices_instanced=int(nv_inst),
             facets_instanced=int(nf_inst),
             vertices_unique=int(sum(len(v) for v in d['mesh_v'])),
             facets_unique=int(sum(len(f) for f in d['mesh_f'])),
             kinds=dict(zip(*[x.tolist() for x in np.unique(d['kind'], return_counts=True)])),
             timings={k: (round(v, 3) if isinstance(v, float) else v) for k, v in d['timings'].items()})
    if 'readout' in d:
        dr = d['node_mesh'] >= 0
        s['drawn_with_readout'] = int(sum(1 for i, r in enumerate(d['readout']) if r and dr[i]))
        s['drawn_with_detelement'] = int(sum(1 for i, r in enumerate(d['det_path']) if r and dr[i]))
        s['drawn_sensitive'] = int((d['sensitive'] & dr).sum())
        if 'g4volid' in d:
            s['sensitive_with_g4_volume_id'] = int(((d['g4volid'] != 0) & d['sensitive']).sum())
            ok = d['sensitive'] & (d['g4volid'] != 0)
            s['g4_vs_dd4hep_volume_id_mismatch'] = int((d['g4volid'][ok].view(np.uint64) != d['volid'][ok]).sum())
    return s


def dump(world_pv, out, nseg=72, formats=('prim', 'npz'), description=None, meta=None, bbox=None,
         edge_flags='dawn'):
    ROOT = cling()
    t0 = time.time()
    d = collect(ROOT, world_pv, nseg=nseg, description=description)
    t_collect = time.time() - t0
    if d['mesh_fail']:
        log('polyhedron not available for %d solids:' % len(d['mesh_fail']), d['mesh_fail'][:20])
    meta = dict(meta or {}, readouts=d.get('readouts', {}), nseg=nseg, source='Geant4 G4VSolid::GetPolyhedron', geant4=str(ROOT.G4Version
                if hasattr(ROOT, 'G4Version') else ''))
    if world_pv:
        e = list(ROOT.solidvis.extent(world_pv))
        c, r = e[:3], e[3]
        meta['world_extent'] = e
        log('world extent centre %s radius %.6g, box %s' % (c, r, e[4:]))
        if bbox is None:
            bbox = [c[0] - r, c[1] - r, c[2] - r, c[0] + r, c[1] + r, c[2] + r]
    if 'prim' in formats:
        write_prim(d, out + '.prim', bbox=bbox, edge_flags=edge_flags)
        if edge_flags == 'both':
            write_prim(d, out + '.g4flags.prim', bbox=bbox, edge_flags='g4')
    if 'npz' in formats:
        write_npz(d, out + '.npz', meta)
    s = summary(d)
    s['collect_seconds'] = round(t_collect, 3)
    log(json.dumps(s))
    with open(out + '.summary.json', 'w') as fh:
        json.dump(dict(s, mesh_fail=d['mesh_fail']), fh, indent=1)
    return d


# ----------------------------------------------------------------------------------------------------
def install(SIM, out, nseg=72, formats=('prim', 'npz')):
    """ddsim steering hook: dump at /run/initialize (after ConstructSDandField), same session as the events."""

    def add_dump_action(kernel):
        import DDG4
        g4 = DDG4.Geant4(kernel)

        def cb(*args):
            ROOT = cling()
            t0 = time.time()
            info = ROOT.dd4hep.sim.Geant4Mapping.instance().data()
            dump(info.world(), out, nseg=nseg, formats=formats, description=kernel.detectorDescription(),
                 meta=dict(mode='ddsim-hook'))
            log('hook total %.2f s' % (time.time() - t0))
            return 1

        g4.addDetectorConstruction('Geant4PythonDetectorConstruction/SolidVisPolyDump', sensitives=cb,
                                   sensitives_args=())

    SIM.physics.setupUserPhysics(add_dump_action)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--compact', required=True)
    ap.add_argument('-o', '--out', required=True, help='output prefix')
    ap.add_argument('--nseg', type=int, default=72)
    ap.add_argument('--formats', default='prim,npz')
    ap.add_argument('--bbox-from', default=None, help='copy /BoundingBox from this .prim')
    ap.add_argument('--edge-flags', default='dawn', choices=('dawn', 'g4', 'both'))
    a = ap.parse_args()
    T0 = time.time()
    ROOT = cling()
    t0 = time.time()
    desc = ROOT.dd4hep.Detector.getInstance()
    desc.fromXML(a.compact)
    t_load = time.time() - t0
    t0 = time.time()
    world = ROOT.solidvis.convert(desc, 2)
    t_conv = time.time() - t0
    log('DD4hep load %.2f s, Geant4 conversion %.2f s' % (t_load, t_conv))
    bbox = None
    if a.bbox_from:
        with open(a.bbox_from) as fh:
            for line in fh:
                if line.startswith('/BoundingBox'):
                    bbox = [float(x) for x in line.split()[1:7]]
                    break
    dump(world, a.out, nseg=a.nseg, formats=tuple(a.formats.split(',')), description=desc,
         meta=dict(mode='standalone', compact=a.compact), bbox=bbox, edge_flags=a.edge_flags)
    log('total wall %.2f s' % (time.time() - T0))


if __name__ == '__main__':
    main()
