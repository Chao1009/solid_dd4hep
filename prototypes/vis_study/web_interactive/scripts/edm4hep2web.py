#!/usr/bin/env python3
"""Prototype: ddsim EDM4hep output -> (a) Firebird DEX 1.0 JSON, (b) JSROOT overlay ROOT file.

usage: edm4hep2web.py events.edm4hep.root display_geo.root out_prefix [--events 0-2]

(a) <out_prefix>.firebird.json  pieces per event:
      <TrackerColl>          BoxHit           SimTrackerHits (pos, dim, time, edep)
      <CaloColl>             BoxHit           SimCalorimeterHits = ECal layers (time = first contribution)
      <CaloColl>Modules      BoxHit           hits summed per ECal module, pos/dim of the module from TGeo
      McHitTrajectories      PointTrajectory  vertex + tracker hits + first calorimeter contribution, time-ordered
      NeutralParticles       PointTrajectory  straight lines vertex->endpoint for photons/neutrons, t at speed c
(b) <out_prefix>.jsroot.root  TObjArray ev<N>: TPolyLine3D tracks, TPolyMarker3D hits,
      TEveGeoShapeExtract with the hit ECal modules coloured by summed energy (log scale)  [cm units]
"""
import argparse, json, math, sys, time
from collections import defaultdict

import numpy as np
import uproot

C_MM_PER_NS = 299.792458


def parse_range(s, n):
    if not s:
        return list(range(n))
    out = []
    for part in s.split(','):
        if '-' in part:
            a, b = map(int, part.split('-')); out += list(range(a, b + 1))
        else:
            out.append(int(part))
    return [i for i in out if i < n]


def colormap(x):
    """x in [0,1] -> rgb, a 'hot'-like perceptual ramp (dark violet -> red -> orange -> yellow)."""
    stops = [(0.00, (0.23, 0.05, 0.42)), (0.35, (0.75, 0.10, 0.35)), (0.65, (0.98, 0.45, 0.10)),
             (1.00, (1.00, 0.95, 0.35))]
    x = min(max(x, 0.0), 1.0)
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if x <= x1:
            f = (x - x0) / (x1 - x0)
            return tuple(a + f * (b - a) for a, b in zip(c0, c1))
    return stops[-1][1]


def energy_fraction(e, emin=1e-3, emax=1.0):
    """log scale between emin and emax [GeV]"""
    if e <= emin:
        return 0.0
    return min(1.0, math.log(e / emin) / math.log(emax / emin))


class Geo:
    """Module lookup in the (pruned) display geometry: hit position -> module node, centre, matrix, bbox."""

    def __init__(self, path):
        import ROOT
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        self.ROOT = ROOT
        self.gm = ROOT.TGeoManager.Import(path)
        self.cache = {}

    def module(self, x_mm, y_mm, z_mm):
        gm = self.gm
        node = gm.FindNode(x_mm / 10, y_mm / 10, z_mm / 10)
        if not node:
            return None
        path = gm.GetPath()
        if path in self.cache:
            return self.cache[path]
        m = gm.GetCurrentMatrix()
        tr = [m.GetTranslation()[i] * 10 for i in range(3)]
        rot = [m.GetRotationMatrix()[i] for i in range(9)]
        shape = node.GetVolume().GetShape()
        bbox = [2 * shape.GetDX() * 10, 2 * shape.GetDY() * 10, 2 * shape.GetDZ() * 10]
        info = dict(path=path, name=node.GetVolume().GetName(), centre=tr, rot=rot, bbox=bbox, shape=shape)
        self.cache[path] = info
        return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input'); ap.add_argument('geometry'); ap.add_argument('out_prefix')
    ap.add_argument('--events', default=None)
    ap.add_argument('--hit-box', type=float, default=10.0, help='tracker hit box size [mm]')
    ap.add_argument('--neutral-emin', type=float, default=0.05, help='GeV, min energy of drawn neutrals')
    ap.add_argument('--no-jsroot', action='store_true')
    a = ap.parse_args()
    t0 = time.time()

    tree = uproot.open(a.input)['events']
    types = tree.typenames(recursive=False)
    trk_colls = sorted(k for k, v in types.items() if v == 'vector<edm4hep::SimTrackerHitData>')
    cal_colls = sorted(k for k, v in types.items() if v == 'vector<edm4hep::SimCalorimeterHitData>')
    keys = set(tree.keys())
    events = parse_range(a.events, tree.num_entries)
    geo = Geo(a.geometry)

    def arrays(names, i):
        return {n: np.asarray(tree[n].array(entry_start=i, entry_stop=i + 1, library='np')[0]) for n in names}

    def rel_branch(coll):
        for cand in (f'_{coll}_particle', f'_{coll}_MCParticle'):
            if cand in keys:
                return f'{cand}/{cand}.index'
        return None

    dex_events, jsroot_events, stats = [], [], []
    for iev in events:
        pieces = []
        mcf = ['PDG', 'charge', 'time', 'mass', 'vertex.x', 'vertex.y', 'vertex.z', 'endpoint.x', 'endpoint.y',
               'endpoint.z', 'momentum.x', 'momentum.y', 'momentum.z', 'generatorStatus']
        mc = arrays([f'MCParticles/MCParticles.{f}' for f in mcf], iev)
        mc = {f: mc[f'MCParticles/MCParticles.{f}'] for f in mcf}
        nmc = len(mc['PDG'])
        pts_by_particle = defaultdict(list)          # mc index -> [(t, x, y, z)]
        all_hits = []                                # for JSROOT markers (mm)

        # ---- tracker hits
        for coll in trk_colls:
            f = ['position.x', 'position.y', 'position.z', 'time', 'eDep']
            d = arrays([f'{coll}/{coll}.{x}' for x in f], iev)
            x, y, z, t, e = (d[f'{coll}/{coll}.{k}'] for k in f)
            n = len(x)
            if n == 0:
                continue
            pieces.append(dict(name=coll, type='BoxHit', version='1.0',
                               origin=dict(type='edm4hep::SimTrackerHitData', name=coll), count=n,
                               columns=dict(pos=np.column_stack([x, y, z]).ravel().round(2).tolist(),
                                            dim=[a.hit_box] * (3 * n), time=t.round(3).tolist(),
                                            edep=e.astype(float).tolist())))
            rb = rel_branch(coll)
            idx = arrays([rb], iev)[rb] if rb else np.full(n, -1)
            for k in range(n):
                all_hits.append((x[k], y[k], z[k]))
                if 0 <= idx[k] < nmc and coll != 'LiteCherenkovHits':
                    pts_by_particle[int(idx[k])].append((float(t[k]), float(x[k]), float(y[k]), float(z[k])))

        # ---- calorimeter hits: per layer cell and per module
        jsroot_modules = []
        for coll in cal_colls:
            f = ['position.x', 'position.y', 'position.z', 'energy', 'contributions_begin', 'contributions_end']
            d = arrays([f'{coll}/{coll}.{x}' for x in f], iev)
            x, y, z, e, cb, ce = (d[f'{coll}/{coll}.{k}'] for k in f)
            n = len(x)
            if n == 0:
                continue
            cc = coll + 'Contributions'
            cf = ['time', 'stepPosition.x', 'stepPosition.y', 'stepPosition.z']
            cd = arrays([f'{cc}/{cc}.{k}' for k in cf], iev)
            ct = cd[f'{cc}/{cc}.time']
            cx, cy, cz = (cd[f'{cc}/{cc}.stepPosition.{k}'] for k in 'xyz')
            crb = rel_branch(cc)
            cidx = arrays([crb], iev)[crb] if crb else np.full(len(ct), -1)
            tfirst = np.array([ct[cb[k]:ce[k]].min() if ce[k] > cb[k] else 0.0 for k in range(n)])
            # first calorimeter contribution of each particle -> extends its trajectory into the ECal
            # ddsim leaves stepPosition at 0 unless configured: then use the cell (layer) centre instead
            hit_of = np.zeros(len(ct), dtype=int)
            for k in range(n):
                hit_of[cb[k]:ce[k]] = k
            use_step = bool(np.any(cx != 0) or np.any(cy != 0) or np.any(cz != 0))
            first = {}
            for j in range(len(ct)):
                p = int(cidx[j])
                if 0 <= p < nmc and (p not in first or ct[j] < first[p][0]):
                    k = hit_of[j]
                    first[p] = (float(ct[j]), *((float(cx[j]), float(cy[j]), float(cz[j])) if use_step
                                                 else (float(x[k]), float(y[k]), float(z[k]))))
            for p, v in first.items():
                if pts_by_particle.get(p):
                    pts_by_particle[p].append(v)
            # modules
            mods = {}
            layer_dz = []
            for k in range(n):
                m = geo.module(x[k], y[k], z[k])
                if m is None:
                    continue
                mm = mods.setdefault(m['path'], dict(m=m, e=0.0, t=1e9, nl=0))
                mm['e'] += float(e[k]); mm['t'] = min(mm['t'], float(tfirst[k])); mm['nl'] += 1
            zs = np.unique(np.round(z, 2))
            layer_dz = float(np.min(np.diff(zs))) if len(zs) > 1 else 5.0
            some = next(iter(mods.values()))['m']['bbox'] if mods else [100, 100, 10]
            pieces.append(dict(name=coll, type='BoxHit', version='1.0',
                               origin=dict(type='edm4hep::SimCalorimeterHitData', name=coll), count=n,
                               columns=dict(pos=np.column_stack([x, y, z]).ravel().round(2).tolist(),
                                            dim=[some[0], some[1], layer_dz] * n, time=tfirst.round(3).tolist(),
                                            edep=e.astype(float).tolist())))
            ml = sorted(mods.values(), key=lambda v: -v['e'])
            pieces.append(dict(name=coll + 'Modules', type='BoxHit', version='1.0',
                               origin=dict(type='edm4hep::SimCalorimeterHitData', name=coll,
                                           note='energy summed per module'), count=len(ml),
                               columns=dict(pos=[round(c, 2) for v in ml for c in v['m']['centre']],
                                            dim=[round(c, 2) for v in ml for c in v['m']['bbox']],
                                            time=[round(v['t'], 3) for v in ml], edep=[v['e'] for v in ml],
                                            nlayers=[v['nl'] for v in ml])))
            jsroot_modules += [(coll, v) for v in ml]

        # ---- MC-truth trajectories from hits
        cols = {k: [] for k in ('pdg', 'charge', 'px', 'py', 'pz', 'p', 'mc_index')}
        tpoints = []
        for p in sorted(pts_by_particle):
            hits = sorted(pts_by_particle[p])
            if len(hits) < 1:
                continue
            pts = [[float(mc['vertex.x'][p]), float(mc['vertex.y'][p]), float(mc['vertex.z'][p]), float(mc['time'][p])]]
            pts += [[hx, hy, hz, ht] for ht, hx, hy, hz in hits]
            px, py, pz = (float(mc[f'momentum.{c}'][p]) for c in 'xyz')
            cols['pdg'].append(int(mc['PDG'][p])); cols['charge'].append(float(mc['charge'][p]))
            cols['px'].append(px * 1000); cols['py'].append(py * 1000); cols['pz'].append(pz * 1000)
            cols['p'].append(math.sqrt(px * px + py * py + pz * pz) * 1000); cols['mc_index'].append(p)
            tpoints.append([[round(v, 3) for v in q] for q in pts])
        if tpoints:
            pieces.append(dict(name='McHitTrajectories', type='PointTrajectory', version='1.0',
                               origin=dict(type=['edm4hep::SimTrackerHitData', 'edm4hep::CaloHitContributionData',
                                                 'edm4hep::MCParticleData']),
                               count=len(tpoints), columns=cols, pointColumns=['x', 'y', 'z', 't'], points=tpoints))

        # ---- neutral particles: straight lines, sampled every 0.2 ns for the time animation
        ncols = {k: [] for k in ('pdg', 'charge', 'px', 'py', 'pz', 'p', 'mc_index')}
        npoints = []
        for p in range(nmc):
            if mc['charge'][p] != 0:
                continue
            px, py, pz = (float(mc[f'momentum.{c}'][p]) for c in 'xyz')
            pp = math.sqrt(px * px + py * py + pz * pz)
            m = float(mc['mass'][p]); en = math.sqrt(pp * pp + m * m)
            if en - m < a.neutral_emin and en < a.neutral_emin:
                continue
            v0 = np.array([mc['vertex.x'][p], mc['vertex.y'][p], mc['vertex.z'][p]], float)
            v1 = np.array([mc['endpoint.x'][p], mc['endpoint.y'][p], mc['endpoint.z'][p]], float)
            L = float(np.linalg.norm(v1 - v0)); beta = pp / en if en > 0 else 1.0
            if L < 1.0:
                continue
            dt = L / (beta * C_MM_PER_NS); t0v = float(mc['time'][p])
            ns = max(1, min(200, int(dt / 0.2)))
            pts = [[*(v0 + (v1 - v0) * s / ns).round(2).tolist(), round(t0v + dt * s / ns, 3)] for s in range(ns + 1)]
            ncols['pdg'].append(int(mc['PDG'][p])); ncols['charge'].append(0.0)
            ncols['px'].append(px * 1000); ncols['py'].append(py * 1000); ncols['pz'].append(pz * 1000)
            ncols['p'].append(pp * 1000); ncols['mc_index'].append(p)
            npoints.append(pts)
        if npoints:
            pieces.append(dict(name='NeutralParticles', type='PointTrajectory', version='1.0',
                               origin=dict(type='edm4hep::MCParticleData', note='straight vertex->endpoint'),
                               count=len(npoints), columns=ncols, pointColumns=['x', 'y', 'z', 't'], points=npoints))

        dex_events.append(dict(id=iev, pieces=pieces))
        jsroot_events.append(dict(tracks=list(zip(cols['pdg'], tpoints)), neutrals=list(zip(ncols['pdg'], npoints)),
                                  hits=all_hits, modules=jsroot_modules))
        stats.append(dict(event=iev, pieces={pc['name']: pc['count'] for pc in pieces}))

    dex = {'type': 'firebird-dex-json', 'version': '1.0',
           'origin': {'file': a.input, 'producer': 'solid edm4hep2web.py prototype'}, 'events': dex_events}
    with open(a.out_prefix + '.firebird.json', 'w') as fo:
        json.dump(dex, fo, separators=(',', ':'))
    t1 = time.time()

    if not a.no_jsroot:
        write_jsroot(geo, jsroot_events, events, a.out_prefix + '.jsroot.root')
    t2 = time.time()
    print(json.dumps(stats))
    print(f'DEX written in {t1 - t0:.1f}s, JSROOT overlay in {t2 - t1:.1f}s')


def write_jsroot(geo, jsroot_events, events, out):
    ROOT = geo.ROOT
    ROOT.gSystem.Load('libEve')
    f = ROOT.TFile(out, 'RECREATE')
    pid_color = {11: ROOT.kRed, -11: ROOT.kBlue, 22: ROOT.kOrange + 1, 211: ROOT.kMagenta, -211: ROOT.kCyan + 2,
                 2212: ROOT.kViolet, 2112: ROOT.kGreen + 2, 13: ROOT.kRed + 2, -13: ROOT.kBlue + 2}
    keep = []
    for iev, ev in zip(events, jsroot_events):
        arr = ROOT.TObjArray(); arr.SetName(f'ev{iev}')
        ROOT.SetOwnership(arr, False)
        for coll, tracks in (('tracks', ev['tracks']), ('neutrals', ev['neutrals'])):
            for pdg, pts in tracks:
                pl = ROOT.TPolyLine3D(len(pts))
                for k, q in enumerate(pts):
                    pl.SetPoint(k, q[0] / 10, q[1] / 10, q[2] / 10)
                pl.SetLineColor(pid_color.get(pdg, ROOT.kGray + 2)); pl.SetLineWidth(2 if coll == 'tracks' else 1)
                arr.Add(pl); keep.append(pl)
        if ev['hits']:
            pm = ROOT.TPolyMarker3D(len(ev['hits']))
            for k, (x, y, z) in enumerate(ev['hits']):
                pm.SetPoint(k, x / 10, y / 10, z / 10)
            pm.SetMarkerColor(ROOT.kAzure + 1); pm.SetMarkerSize(0.4); pm.SetMarkerStyle(20)
            arr.Add(pm); keep.append(pm)
        if ev['modules']:
            top = ROOT.TEveGeoShapeExtract(f'calo_ev{iev}', 'ECal modules coloured by energy')
            ROOT.SetOwnership(top, False)  # extracts own (and would delete) the shared TGeo shapes
            top.SetRnrSelf(False); top.SetRnrElements(True)
            for coll, v in ev['modules']:
                m = v['m']
                ex = ROOT.TEveGeoShapeExtract(f"{coll}_{v['e'] * 1000:.1f}MeV", m['path'])
                ROOT.SetOwnership(ex, False)
                r = m['rot']; t = [c / 10 for c in m['centre']]
                trans = [r[0], r[3], r[6], 0, r[1], r[4], r[7], 0, r[2], r[5], r[8], 0, t[0], t[1], t[2], 1]
                ex.SetTrans(np.array(trans, dtype=np.float64))
                rgb = colormap(energy_fraction(v['e']))
                ex.SetRGBA(np.array([rgb[0], rgb[1], rgb[2], 0.85], dtype=np.float32))
                ex.SetShape(m['shape'])
                ex.SetRnrSelf(True); ex.SetRnrElements(True)
                top.AddElement(ex); keep.append(ex)
            arr.Add(top); keep.append(top)
        arr.Write(f'ev{iev}', ROOT.TObject.kSingleKey)
    f.Close()


if __name__ == '__main__':
    main()
    sys.stdout.flush()
    import os
    os._exit(0)  # skip ROOT/cppyy teardown (TEve extracts share TGeo shapes)
