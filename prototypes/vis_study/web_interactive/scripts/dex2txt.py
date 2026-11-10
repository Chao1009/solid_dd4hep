#!/usr/bin/env python3
"""Flatten one event of the prototype DEX JSON (edm4hep2web.py output) to a text file for the REve macro.
lines:  T pdg charge vx vy vz px py pz  [cm, GeV]   |  H x y z [cm]  |  B x y z dx dy dz E_MeV [cm]
usage: dex2txt.py events.firebird.json EVENT out.txt"""
import json, sys
dex = json.load(open(sys.argv[1]))
ev = [e for e in dex['events'] if e['id'] == int(sys.argv[2])][0]
with open(sys.argv[3], 'w') as f:
    for p in ev['pieces']:
        c = p['columns']
        if p['type'] == 'PointTrajectory' and p['name'] == 'McHitTrajectories':
            for k in range(p['count']):
                v = p['points'][k][0]
                f.write(f"T {c['pdg'][k]} {c['charge'][k]:g} {v[0]/10:.3f} {v[1]/10:.3f} {v[2]/10:.3f} "
                        f"{c['px'][k]/1e3:.5f} {c['py'][k]/1e3:.5f} {c['pz'][k]/1e3:.5f}\n")
        elif p['type'] == 'BoxHit' and p['name'].endswith('Modules'):
            for k in range(p['count']):
                x, y, z = c['pos'][3*k:3*k+3]; dx, dy, dz = c['dim'][3*k:3*k+3]
                f.write(f"B {x/10:.3f} {y/10:.3f} {z/10:.3f} {dx/10:.3f} {dy/10:.3f} {dz/10:.3f} {c['edep'][k]*1e3:.3f}\n")
        elif p['type'] == 'BoxHit' and 'EC' not in p['name']:
            for k in range(p['count']):
                x, y, z = c['pos'][3*k:3*k+3]
                f.write(f"H {x/10:.3f} {y/10:.3f} {z/10:.3f}\n")
