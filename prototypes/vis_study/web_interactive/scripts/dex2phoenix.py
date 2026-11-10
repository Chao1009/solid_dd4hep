#!/usr/bin/env python3
"""Prototype DEX (edm4hep2web.py output) -> HSF Phoenix event JSON.
Tracks (MC-truth polylines), Hits (tracker hits), IrregularCaloCells (ECal modules, axis-aligned boxes
with RGB colour from the summed module energy, log scale 1 MeV..1 GeV).  Units: mm.
usage: dex2phoenix.py events.firebird.json out.phoenix.json"""
import json, math, sys

STOPS = [(0.00, (59, 13, 107)), (0.35, (191, 26, 89)), (0.65, (250, 115, 26)), (1.00, (255, 242, 89))]


def colour(e_gev, emin=1e-3, emax=1.0):
    x = 0.0 if e_gev <= emin else min(1.0, math.log(e_gev / emin) / math.log(emax / emin))
    for (x0, c0), (x1, c1) in zip(STOPS, STOPS[1:]):
        if x <= x1:
            f = (x - x0) / (x1 - x0)
            return [round(a + f * (b - a)) for a, b in zip(c0, c1)]
    return list(STOPS[-1][1])


dex = json.load(open(sys.argv[1]))
out = {}
for ev in dex['events']:
    e = {'event number': ev['id'], 'run number': 0, 'Tracks': {}, 'Hits': {}, 'IrregularCaloCells': {}}
    hits, cells = [], []
    for p in ev['pieces']:
        c = p['columns']
        if p['type'] == 'PointTrajectory':
            col = 'ff2020' if p['name'] == 'McHitTrajectories' else 'ffa500'
            e['Tracks'][p['name']] = [
                {'pos': [q[:3] for q in p['points'][k]], 'color': col, 'charge': c['charge'][k], 'mom': c['p'][k],
                 'linewidth': 3 if col == 'ff2020' else 1} for k in range(p['count'])]
        elif p['type'] == 'BoxHit' and p['name'].endswith('Modules'):
            for k in range(p['count']):
                x, y, z = c['pos'][3 * k:3 * k + 3]
                dx, dy, dz = (v / 2 for v in c['dim'][3 * k:3 * k + 3])
                vtx = [v for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1) for v in (x + sx * dx, y + sy * dy, z + sz * dz)]
                cells.append({'layer': 0 if 'PreShower' in p['name'] else 1, 'vtx': vtx, 'color': colour(c['edep'][k]),
                              'opacity': 0.9, 'energy': c['edep'][k] * 1000})
        elif p['type'] == 'BoxHit' and 'EC' not in p['name']:
            hits += [{'type': 'CircularPoint', 'pos': c['pos'][3 * k:3 * k + 3], 'color': '00a0ff', 'size': 12}
                     for k in range(p['count'])]
    e['Hits']['TrackerHits'] = hits
    e['IrregularCaloCells']['ECalModules'] = cells
    out[f'SoLID-Lite e- 4 GeV, event {ev["id"]}'] = e
json.dump(out, open(sys.argv[2], 'w'))
print({k: {t: {n: len(v) for n, v in d.items()} for t, d in ev_.items() if isinstance(d, dict)} for k, ev_ in out.items()})
