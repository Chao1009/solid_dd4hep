#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check that every simulated hit maps to exactly one drawn volume of the geometry scene.

usage (inside eic-shell / solid-shell, after `source thissolid.sh <config>`):
    check_scene_hits.py [-c CONFIG] [-n EVENTS] [-o OUTDIR]            run ddsim, then check
    check_scene_hits.py --scene scene.npz --hits sim.edm4hep.root      check existing files

A drawn volume of the scene carries the readout it belongs to and the volume ID bits fixed by the placements
down to it, with their mask. A hit of that readout lies in the volume if (cellID & mask) == id. Several nested
volumes can match (a calorimeter module and the plate it stands on); the hit belongs to the most specific one,
the match with the most mask bits. The test fails if a hit matches no drawn volume, or if two volumes tie for
the most specific match. An event display relies on this rule to colour volumes by deposited energy; the
rule itself is implemented once, in python/solid_evdisplay/hitmap.py, which solid_evexport uses as well.

Without --scene/--hits the script runs ddsim with steering/solid_sim.py, SOLID_SCENE (so the scene comes from
the same Geant4 session as the hits), SOLID_OPTICAL=0 and the default particle gun, for EVENTS events
(default 3). Needs uproot.

Exit status: 0 = every hit has exactly one most specific volume, 1 = not, 2 = the run failed.
"""
import argparse
import os
import subprocess
import sys

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'python'))


def run_ddsim(config, events, outdir):
    path = os.environ.get('DETECTOR_PATH', '.')
    scene = os.path.join(outdir, 'scene.npz')
    hits = os.path.join(outdir, 'hits.edm4hep.root')
    cmd = ['ddsim', '--steeringFile', os.path.join(path, 'steering', 'solid_sim.py'), '--compactFile',
           os.path.join(path, config + '.xml'), '-N', str(events), '--outputFile', hits]
    env = dict(os.environ, SOLID_SCENE=scene, SOLID_OPTICAL='0', DETECTOR_CONFIG=config)
    with open(os.path.join(outdir, 'ddsim.log'), 'w') as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=outdir, check=False)
    if r.returncode != 0 or not os.path.isfile(scene) or not os.path.isfile(hits):
        print('check_scene_hits: ddsim failed or wrote no scene, see %s/ddsim.log' % outdir)
        return None, None
    return scene, hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--config', default=os.environ.get('DETECTOR_CONFIG', 'solid'))
    ap.add_argument('-n', '--events', type=int, default=3)
    ap.add_argument('-o', '--outdir', default='scene_hits_out')
    ap.add_argument('--scene', default=None)
    ap.add_argument('--hits', default=None, help='EDM4hep ROOT file of a ddsim run on the same geometry')
    a = ap.parse_args()

    import uproot
    scene, hits = a.scene, a.hits
    if not (scene and hits):
        os.makedirs(a.outdir, exist_ok=True)
        scene, hits = run_ddsim(a.config, a.events, os.path.abspath(a.outdir))
        if not scene:
            return 2

    from solid_evdisplay.hitmap import NO_MATCH, TIE, VolumeMap   # the rule solid_evexport uses
    vmap = VolumeMap(scene)
    tree = uproot.open(hits)['events']
    collections = sorted({k.split('/')[0] for k in tree.keys() if k.endswith('.cellID') and not k.startswith('_')})
    total = bad = 0
    print('%-22s %7s %9s %9s %9s  %s' % ('collection', 'hits', 'unique', 'no match', 'tie', 'matches per hit'))
    for c in collections:
        per_event = tree['%s/%s.cellID' % (c, c)].array(library='np')
        cid = np.concatenate(list(per_event)).astype(np.uint64) if len(per_event) else np.zeros(0, np.uint64)
        if len(cid) == 0:
            continue
        node, nmatch = vmap.match(c, cid)
        none = int((node == NO_MATCH).sum())
        tie = int((node == TIE).sum())
        unique = int((node >= 0).sum())
        dist = dict(zip(*[x.tolist() for x in np.unique(nmatch, return_counts=True)]))
        print('%-22s %7d %9d %9d %9d  %s' % (c, len(cid), unique, none, tie, dist))
        total += len(cid)
        bad += none + tie
    if total == 0:
        print('no hits in %s: nothing checked' % hits)
        return 1
    print('%d hits, %d without exactly one most specific drawn volume' % (total, bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
