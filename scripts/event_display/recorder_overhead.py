#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure what the trajectory recorder (SOLID_TRAJ) costs in simulation time.

usage (inside eic-shell / solid-shell):
    recorder_overhead.py [-c CONFIG] [-n EVENTS] [--reps R] [--input EVENTS.hepmc] [--seed S] [-o DIR]

The same events (same generator file or gun, same seed) are simulated R times without and R times with the
recorder, alternating (off on, on off, off on, ...) so that a drift of the machine load hits both alike.
Per run it takes

    loop CPU    user + system time of the Geant4 event loop: the line "User=...s Real=...s Sys=...s" that
                Geant4 prints at the end of the run with /run/verbose 1 (start-up and geometry not included)
    loop real   the "Real" of the same line
    job CPU     user + system time of the whole ddsim process (os.wait4)

and prints mean and standard deviation per mode, and the overhead three ways: from the means, from the
minima (the least disturbed run of each mode), and as the mean of the R paired differences with its
standard error. CPU time is the number to trust on a shared machine; real time is given for comparison.
The trajectory file is deleted after each run; its size and the number of points are reported.
SOLID_OPTICAL=0 unless the environment says otherwise. Results also go to DIR/overhead.json.
"""
import argparse
import json
import os
import re
import statistics
import subprocess
import sys
import time

G4_TIMER = re.compile(r'User=([0-9.eE+-]+)s\s+Real=([0-9.eE+-]+)s\s+Sys=([0-9.eE+-]+)s')


def run(a, tag, with_traj):
    path = os.environ.get('DETECTOR_PATH', '.')
    out = os.path.join(a.outdir, 'run.edm4hep.root')
    traj = os.path.join(a.outdir, 'traj.bin')
    log = os.path.join(a.outdir, 'ddsim_%s.log' % tag)
    cmd = ['ddsim', '--steeringFile', os.path.join(path, 'steering', 'solid_sim.py'), '--compactFile',
           os.path.join(path, a.config + '.xml'), '-N', str(a.events), '--random.seed', str(a.seed),
           '--outputFile', out, '--ui.commandsPreRun', '/run/verbose 1']
    env = dict(os.environ, DETECTOR_CONFIG=a.config)
    env.setdefault('SOLID_OPTICAL', '0')
    env.pop('SOLID_TRAJ', None)
    env.pop('SOLID_SCENE', None)
    if a.input:
        cmd += ['--inputFiles'] + a.input
        env['SOLID_GUN'] = '0'
    if with_traj:
        env['SOLID_TRAJ'] = traj
    t0 = time.perf_counter()
    with open(log, 'w') as f:
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=a.outdir)
        _, status, ru = os.wait4(p.pid, 0)
    wall = time.perf_counter() - t0
    text = open(log, errors='replace').read()
    m = G4_TIMER.findall(text)
    saved = len(re.findall(r'Saving EDM4hep event', text))
    if status != 0 or not m or not os.path.isfile(out):
        sys.exit('recorder_overhead: ddsim failed or printed no run timer, see %s' % log)
    user, real, syst = (float(v) for v in m[-1])
    r = dict(tag=tag, traj=with_traj, loop_cpu=user + syst, loop_real=real, job_cpu=ru.ru_utime + ru.ru_stime,
             job_wall=wall, events=saved or a.events, edm4hep_bytes=os.path.getsize(out))
    if with_traj:
        r['traj_bytes'] = os.path.getsize(traj) if os.path.isfile(traj) else 0
        os.remove(traj)
    os.remove(out)
    print('%-8s loop CPU %8.2f s  loop real %8.2f s  job CPU %8.2f s  job wall %8.2f s%s' % (
        tag, r['loop_cpu'], r['loop_real'], r['job_cpu'], r['job_wall'],
        '  traj %.1f MB' % (r['traj_bytes'] / 1e6) if with_traj else ''), flush=True)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--config', default=os.environ.get('DETECTOR_CONFIG', 'solid_JPsi_lite'))
    ap.add_argument('-n', '--events', type=int, default=100)
    ap.add_argument('--reps', type=int, default=5, help='runs per mode (default 5)')
    ap.add_argument('--input', nargs='+', default=None, help='generator file(s) for --inputFiles (default: '
                    'the particle gun of the steering file)')
    ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('-o', '--outdir', default='recorder_overhead')
    a = ap.parse_args()
    a.outdir = os.path.abspath(a.outdir)
    if a.input:
        a.input = [os.path.abspath(f) for f in a.input]
    os.makedirs(a.outdir, exist_ok=True)
    runs = []
    for rep in range(a.reps):
        for with_traj in ((False, True) if rep % 2 == 0 else (True, False)):
            runs.append(dict(run(a, '%s%d' % ('on' if with_traj else 'off', rep), with_traj), rep=rep))
    off = sorted((r for r in runs if not r['traj']), key=lambda r: r['rep'])
    on = sorted((r for r in runs if r['traj']), key=lambda r: r['rep'])
    summary = dict(config=a.config, events=a.events, reps=a.reps, input=a.input, seed=a.seed, runs=runs)
    print('\n%s, %d events, %d runs per mode, %s' % (a.config, a.events, a.reps,
                                                   ', '.join(a.input) if a.input else 'particle gun'))
    print('%-10s %22s %22s %12s %12s %12s' % ('', 'without recorder [s]', 'with recorder [s]', 'means',
                                              'minima', 'pairs'))
    for key, name in (('loop_cpu', 'loop CPU'), ('loop_real', 'loop real'), ('job_cpu', 'job CPU')):
        x, y = [r[key] for r in off], [r[key] for r in on]
        sd = (lambda v: statistics.stdev(v) if len(v) > 1 else 0.0)
        pairs = [100.0 * (q - p) / p for p, q in zip(x, y)]
        res = dict(off_mean=statistics.mean(x), off_sd=sd(x), on_mean=statistics.mean(y), on_sd=sd(y),
                   overhead_means=100.0 * (statistics.mean(y) / statistics.mean(x) - 1.0),
                   overhead_minima=100.0 * (min(y) / min(x) - 1.0), overhead_pairs=statistics.mean(pairs),
                   overhead_pairs_sd=sd(pairs), overhead_pairs_sem=sd(pairs) / len(pairs) ** 0.5)
        summary[key] = res
        print('%-10s %12.2f +- %6.2f %12.2f +- %6.2f %+11.2f%% %+11.2f%% %+6.2f +- %4.2f%%' % (
            name, res['off_mean'], res['off_sd'], res['on_mean'], res['on_sd'], res['overhead_means'],
            res['overhead_minima'], res['overhead_pairs'], res['overhead_pairs_sem']))
    tb = statistics.mean(r['traj_bytes'] for r in on)
    summary.update(traj_mb_per_event=tb / 1e6 / a.events, points_per_event=(tb - 40) / 40.0 / a.events,
                   edm4hep_bytes_off=off[0]['edm4hep_bytes'], edm4hep_bytes_on=on[0]['edm4hep_bytes'])
    print('+- is the standard deviation of the runs; for the pairs, the standard error of their mean')
    print('trajectory file: %.2f MB and %.0f points per event; EDM4hep file %d bytes without, %d with' % (
        summary['traj_mb_per_event'], summary['points_per_event'], summary['edm4hep_bytes_off'],
        summary['edm4hep_bytes_on']))
    json.dump(summary, open(os.path.join(a.outdir, 'overhead.json'), 'w'), indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
