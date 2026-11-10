#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Smoke test of the event display chain: recorder -> solid_evexport -> solid_evview (1 event, 2 views).

usage (inside eic-shell / solid-shell; ctest runs it as test `event_display`):
    test_event_display.py [-c CONFIG] [-o OUTDIR] [--keep]

Steps, each checked:
  1. ddsim, 1 gun event with SOLID_SCENE and SOLID_TRAJ (SOLID_OPTICAL=0, fixed seed)
  2. the trajectory file: header, one event, a primary track with increasing times
  3. solid_evexport: exit status 0 (every hit has one drawn volume, energies add up), arrays of the event file
  4. solid_evview: side and beam views as SVG and PNG; the SVG is well-formed and holds the track colour
  5. solid_evview --animate with 8 frames: GIF with 8 frames whose last frame is the static figure, and a
     well-formed animated SVG; MP4 is requested too and must be skipped without an error when no ffmpeg exists

Needs the built plugins on LD_LIBRARY_PATH, DETECTOR_PATH, and the python packages of this repository
(found next to this script in a source checkout). About 30 s. Exit status 0 = passed.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'python'))

FAILED = []


def check(ok, what):
    print('%s  %s' % ('ok  ' if ok else 'FAIL', what), flush=True)
    if not ok:
        FAILED.append(what)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--config', default='solid_JPsi_lite')
    ap.add_argument('-o', '--outdir', default=None, help='work directory (default: a temporary one)')
    ap.add_argument('--seed', type=int, default=11)
    ap.add_argument('--keep', action='store_true')
    a = ap.parse_args()
    from solid_evdisplay import animate, export, load_event, trajfile, view

    out = os.path.abspath(a.outdir) if a.outdir else tempfile.mkdtemp(prefix='event_display_test_')
    os.makedirs(out, exist_ok=True)
    path = os.environ.get('DETECTOR_PATH', '.')
    scene, traj, sim = (os.path.join(out, f) for f in ('scene.npz', 'traj.bin', 'sim.edm4hep.root'))
    # 1 ---------------------------------------------------------------------------------------------
    cmd = ['ddsim', '--steeringFile', os.path.join(path, 'steering', 'solid_sim.py'), '--compactFile',
           os.path.join(path, a.config + '.xml'), '-N', '1', '--random.seed', str(a.seed), '--outputFile', sim]
    env = dict(os.environ, SOLID_SCENE=scene, SOLID_TRAJ=traj, SOLID_OPTICAL='0', DETECTOR_CONFIG=a.config)
    with open(os.path.join(out, 'ddsim.log'), 'w') as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=out, check=False)
    if not check(r.returncode == 0 and all(os.path.isfile(f) for f in (scene, traj, sim)),
                 'ddsim wrote the scene, the trajectory file and the EDM4hep file (log: %s/ddsim.log)' % out):
        return 1
    # 2 ---------------------------------------------------------------------------------------------
    hdr = trajfile.read_header(traj)
    pts = trajfile.read_trajectories(traj)
    check(hdr['version'] == 1 and abs(hdr['min_ekin'] - 5.0) < 1e-6, 'trajectory header: version 1, 5 MeV cut')
    check(len(pts) > 10 and set(pts['event']) == {0}, 'trajectory points of event 0 only (%d points)' % len(pts))
    prim = pts[pts['parent'] == 0]
    check(len(prim) >= 2 and bool(np.all(np.diff(prim[prim['track'] == prim['track'][0]]['t']) >= 0))
          and bool(prim['flags'][0] & 1), 'primary track: first-point flag, times increase')
    sec = pts[pts['parent'] != 0]
    first = sec[(sec['flags'] & 1) != 0]
    check(len(first) == 0 or float(first['ekin'].min()) >= 5.0, 'no secondary below the 5 MeV cut')
    # 3 ---------------------------------------------------------------------------------------------
    evdir = os.path.join(out, 'events')
    check(export.main([sim, '--scene', scene, '--traj', traj, '-o', evdir, '-q']) == 0,
          'solid_evexport: every hit has one most specific drawn volume, energies add up')
    evfile = os.path.join(evdir, 'event_000000.npz')
    if not check(os.path.isfile(evfile), 'event_000000.npz written'):
        return 1
    ev = load_event(evfile)
    check(len(ev.hit_node) > 0 and bool((ev.hit_node >= 0).all()), '%d hits, all mapped' % len(ev.hit_node))
    for ci, c in enumerate(ev.collections):
        ev_sum = float(ev.vol_energy[ev.vol_collection == ci].astype(float).sum())
        if not abs(ev_sum - c['energy']) <= 1e-4 * max(c['energy'], 1e-9):
            check(False, 'energy per volume of %s: %g against %g MeV' % (c['name'], ev_sum, c['energy']))
    check(len(ev.trk_id) > 0 and ev.trk_offset[-1] == len(ev.pt_time) == len(pts),
          'tracks: %d with %d points' % (len(ev.trk_id), len(ev.pt_time)))
    check(ev.calorimeter_energy() > 100.0, 'calorimeter deposit %.0f MeV' % ev.calorimeter_energy())
    # 4 ---------------------------------------------------------------------------------------------
    figs = os.path.join(out, 'figs')
    check(view.main([evdir, '--scene', scene, '-o', figs, '--views', 'side,beam', '--formats', 'svg,png',
                     '--name', 'test_{view}', '-q']) == 0, 'solid_evview side,beam')
    for v in ('side', 'beam'):
        svg, png = os.path.join(figs, 'test_%s.svg' % v), os.path.join(figs, 'test_%s.png' % v)
        ok = os.path.isfile(svg) and os.path.isfile(png) and os.path.getsize(png) > 20000
        if ok:
            ET.parse(svg)
            text = open(svg, encoding='utf-8').read()
            ok = '#2a78d6' in text and 'energy deposit per volume' in text      # e- track colour, colour bar
        check(ok, '%s view: SVG well-formed with track and colour bar, PNG written' % v)
    # 5 ---------------------------------------------------------------------------------------------
    env_ffmpeg = os.environ.pop('SOLID_FFMPEG', None)
    have_ffmpeg = animate.find_ffmpeg() is not None
    rc = view.main([evdir, '--scene', scene, '-o', figs, '--views', '', '--animate', 'side', '--frames', '8',
                    '--anim-formats', 'gif,svg,mp4', '--name', 'test_{view}', '-q'])
    if env_ffmpeg is not None:
        os.environ['SOLID_FFMPEG'] = env_ffmpeg
    check(rc == 0, 'solid_evview --animate (ffmpeg %s)' % ('found' if have_ffmpeg else 'not found: MP4 skipped'))
    gif, asvg = os.path.join(figs, 'test_side_anim.gif'), os.path.join(figs, 'test_side_anim.svg')
    ok = os.path.isfile(asvg)
    if ok:
        root = ET.parse(asvg).getroot()
        ok = sum(1 for e in root.iter() if e.tag.endswith('animate')) > 5
    check(ok, 'animated SVG well-formed, with animate elements')
    check(os.path.isfile(os.path.join(figs, 'test_side_anim.mp4')) == have_ffmpeg, 'MP4 written only with ffmpeg')
    try:
        from PIL import Image
    except ImportError:
        Image = None
    if Image is not None and check(os.path.isfile(gif), 'GIF written'):
        im = Image.open(gif)
        check(im.n_frames == 8, 'GIF has 8 frames')
        im.seek(im.n_frames - 1)
        last = np.asarray(im.convert('RGB')).astype(int)
        # the static figure at the frame size, through the same palette-free path
        a_ = view.build_parser().parse_args([evdir, '--scene', scene, '-q'])
        import solid_lineview as slv
        slv.core.VERBOSE = 0
        ev_view = view.EventView(slv.load_scene(scene), ev, 'side', a_, '%s, event 0' % a.config)
        ref = animate.FrameRenderer(ev_view, a_.frame_scale).static().astype(int)
        same = last.shape == ref.shape
        diff = float(np.abs(last - ref).mean()) if same else 255.0
        check(same and diff < 2.0, 'last GIF frame equals the static figure (mean difference %.2f of 255, palette '
              'rounding)' % diff)
    json.dump(dict(failed=FAILED), open(os.path.join(out, 'result.json'), 'w'))
    if not a.keep and not a.outdir and not FAILED:
        shutil.rmtree(out, ignore_errors=True)
    else:
        print('work directory: %s' % out)
    print('event display test: %s' % ('FAILED: ' + '; '.join(FAILED) if FAILED else 'passed'))
    return 1 if FAILED else 0


if __name__ == '__main__':
    sys.exit(main())
