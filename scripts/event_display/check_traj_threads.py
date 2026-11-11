#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check that the trajectory recorder gives the same events with one and with several threads.

usage (inside eic-shell / solid-shell, after `source thissolid.sh <config>`):
    check_traj_threads.py [-c CONFIG] [-n EVENTS] [-j THREADS] [-o OUTDIR] [--seed SEED]
    check_traj_threads.py --compare DIR1 DIR2          compare two directories of event.npz files

Runs ddsim twice with steering/solid_sim.py, SOLID_SCENE, SOLID_TRAJ, SOLID_OPTICAL=0, the default particle
gun, the same --random.seed and --random.enableEventSeed (so that an event has the same random numbers in
whichever thread it runs): once with -j 1 and once with -j THREADS (default 2).  Both runs are exported with
solid_evexport, and every array of every event must be identical: trajectory points, hits, energies per
volume, MCParticles.  This also shows that the single trajectory file written by several threads is complete.

Exit status: 0 = identical, 1 = different, 2 = a run or an export failed.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'python'))


def run(config, events, threads, seed, outdir):
    from solid_evdisplay import export
    path = os.environ.get('DETECTOR_PATH', '.')
    os.makedirs(outdir, exist_ok=True)
    cmd = ['ddsim', '--steeringFile', os.path.join(path, 'steering', 'solid_sim.py'), '--compactFile',
           os.path.join(path, config + '.xml'), '-N', str(events), '-j', str(threads), '--random.seed', str(seed),
           '--random.enableEventSeed', '--outputFile', 'sim.edm4hep.root']
    env = dict(os.environ, SOLID_SCENE='scene.npz', SOLID_TRAJ='traj.bin', SOLID_OPTICAL='0',
               DETECTOR_CONFIG=config)
    with open(os.path.join(outdir, 'ddsim.log'), 'w') as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=outdir, check=False)
    files = [os.path.join(outdir, f) for f in ('sim.edm4hep.root', 'scene.npz', 'traj.bin')]
    if r.returncode != 0 or not all(os.path.isfile(f) for f in files):
        print('check_traj_threads: ddsim -j %d failed, see %s/ddsim.log' % (threads, outdir))
        return None
    evdir = os.path.join(outdir, 'events')
    if export.main([files[0], '--scene', files[1], '--traj', files[2], '-o', evdir, '-q']) != 0:
        return None
    return evdir


def compare(dir1, dir2):
    """number of differing events (and prints them)"""
    f1 = sorted(f for f in os.listdir(dir1) if f.startswith('event_') and f.endswith('.npz'))
    f2 = sorted(f for f in os.listdir(dir2) if f.startswith('event_') and f.endswith('.npz'))
    bad = 0
    if f1 != f2:
        print('different event lists: %d and %d files' % (len(f1), len(f2)))
        bad += 1
    npts = nhits = 0
    for f in sorted(set(f1) & set(f2)):
        a, b = np.load(os.path.join(dir1, f)), np.load(os.path.join(dir2, f))
        diff = [k for k in a.files if k != 'meta' and not (a[k].shape == b[k].shape and np.array_equal(a[k], b[k]))]
        ma, mb = json.loads(str(a['meta'])), json.loads(str(b['meta']))
        if ma['collections'] != mb['collections']:
            diff.append('meta collections')
        npts += len(a['pt_time'])
        nhits += len(a['hit_node'])
        if diff:
            bad += 1
            print('%s differs in %s' % (f, ', '.join(diff)))
    print('%d events compared (%d trajectory points, %d hits): %d differ' % (len(set(f1) & set(f2)), npts, nhits,
                                                                            bad))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--config', default=os.environ.get('DETECTOR_CONFIG', 'solid'))
    ap.add_argument('-n', '--events', type=int, default=6)
    ap.add_argument('-j', '--threads', type=int, default=2)
    ap.add_argument('--seed', type=int, default=4242)
    ap.add_argument('-o', '--outdir', default='traj_threads_out')
    ap.add_argument('--compare', nargs=2, metavar='DIR', default=None)
    a = ap.parse_args()
    if a.compare:
        return 1 if compare(*a.compare) else 0
    out = os.path.abspath(a.outdir)
    d1 = run(a.config, a.events, 1, a.seed, os.path.join(out, 'j1'))
    dn = run(a.config, a.events, a.threads, a.seed, os.path.join(out, 'j%d' % a.threads))
    if not (d1 and dn):
        return 2
    return 1 if compare(d1, dn) else 0


if __name__ == '__main__':
    sys.exit(main())
