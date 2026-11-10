#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare the geometry scene (solid_scene_dump) with Geant4's DAWNFILE output for one configuration.

usage (inside eic-shell / solid-shell, after `source thissolid.sh <config>`):
    compare_scene_dawnfile.py [-c CONFIG] [-o OUTDIR] [--scene FILE] [--prim FILE] [--views v1,v2,...]
                              [--max-diff PERCENT] [--no-render]

Both files describe the same Geant4 geometry: the scene is written by this repository's tree walk, the .prim by
the DAWNFILE vis driver in a ddsim run. The script checks, volume by volume in drawing order,
    count      the same number of drawn volumes
    names      the same name:copy (volumes DAWNFILE writes without a name are counted separately; Geant4 11.4
               does that for extruded solids, e.g. the 330 PVDIS baffle blades)
    depths     the same depth in the volume tree (for the named volumes)
    colours    equal within 1e-5
    boxes      bounding boxes equal within 1e-3 mm
    sphere     the bounding sphere of the scene equals DAWNFILE's /BoundingBox (this sets --mag framing)
and renders the views with make_solid_views from both files and counts the pixels that differ by more than
8/255 in any channel. The comparison fails if a check fails or a view differs in more than --max-diff percent
of its pixels (default 0.1). Diff images (red: only in the scene picture, blue: only in the DAWNFILE one) and
summary.json are written to OUTDIR. Without --scene / --prim the files are produced first (the DAWNFILE run takes
about 15 s to 1 min).

Exit status: 0 = equal within the limits, 1 = differences, 2 = a file could not be produced.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

MACRO = """/vis/open DAWNFILE
/vis/viewer/set/lineSegmentsPerCircle {segments}
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/flush
"""


def tool(name):
    for cand in (shutil.which(name) or '', os.path.join(REPO, 'bin', name), os.path.join(REPO, 'build', 'bin', name)):
        if cand and os.path.isfile(cand):
            return cand
    sys.exit('compare_scene_dawnfile: %s not found (source thissolid.sh)' % name)


def make_scene(compact, workdir, segments):
    scene = os.path.join(workdir, 'scene.npz')
    with open(os.path.join(workdir, 'solid_scene_dump.log'), 'w') as log:
        subprocess.run([tool('solid_scene_dump'), compact, '-o', scene, '--segments', str(segments)], stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    return scene if os.path.isfile(scene) else None


def make_prim(compact, workdir, segments):
    macro = os.path.join(workdir, 'dawnfile.mac')
    with open(macro, 'w') as f:
        f.write(MACRO.format(segments=segments))
    cmd = ['ddsim', '--compactFile', compact, '--runType', 'vis', '--macroFile', macro,
           '--outputFile', os.path.join(workdir, 'dawnfile_sim.edm4hep.root')]
    steering = os.path.join(os.environ.get('DETECTOR_PATH', '.'), 'steering', 'solid_sim.py')
    if os.path.exists(steering):
        cmd[1:1] = ['--steeringFile', steering]
    env = dict(os.environ, G4DAWNFILE_VIEWER='true', SOLID_OPTICAL='0')
    env.pop('SOLID_SCENE', None)
    with open(os.path.join(workdir, 'ddsim.log'), 'w') as log:
        subprocess.run(cmd, input=b'exit\n', stdout=log, stderr=subprocess.STDOUT, env=env, cwd=workdir,
                       timeout=3600, check=False)
    prims = sorted(p for p in os.listdir(workdir) if p.endswith('.prim'))
    return os.path.join(workdir, prims[-1]) if prims else None


def read_png(path):
    """(H, W, 3) float array in 0..1"""
    try:
        import matplotlib.image as mi
        img = np.asarray(mi.imread(path), float)
        if img.ndim == 2:                            # greyscale PNG (line drawings)
            img = np.repeat(img[..., None], 3, -1)
        return img[..., :3]
    except ImportError:
        import cairo
        s = cairo.ImageSurface.create_from_png(path)
        a = np.frombuffer(s.get_data(), np.uint8).reshape(s.get_height(), s.get_stride() // 4, 4)[:, :s.get_width()]
        return a[..., [2, 1, 0]] / 255.0


def write_diff(a, b, diff, path):
    """grey copy of the scene picture with the differing pixels in red (darker in the scene) or blue"""
    try:
        import matplotlib.image as mi
    except ImportError:
        return
    out = np.repeat((0.6 + 0.4 * a.mean(-1))[..., None], 3, -1)
    darker = a.sum(-1) <= b.sum(-1)
    out[diff & darker] = (0.9, 0.0, 0.0)
    out[diff & ~darker] = (0.0, 0.2, 0.9)
    mi.imsave(path, out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-c', '--config', default=os.environ.get('DETECTOR_CONFIG', 'solid'))
    ap.add_argument('-o', '--outdir', default='compare_scene_out')
    ap.add_argument('--scene', default=None, help='existing scene.npz (default: run solid_scene_dump)')
    ap.add_argument('--prim', default=None, help='existing DAWNFILE .prim (default: run ddsim with the driver)')
    ap.add_argument('--views', default='section_yz,section_xz,cutaway,lines_yz')
    ap.add_argument('--segments', type=int, default=72)
    ap.add_argument('--max-diff', type=float, default=0.1, help='largest accepted share of differing pixels, in '
                    'percent (default 0.1)')
    ap.add_argument('--no-render', action='store_true', help='compare the volume lists only')
    a = ap.parse_args()

    try:
        import solid_lineview as slv
    except ImportError:
        sys.path.insert(0, os.path.join(REPO, 'python'))
        import solid_lineview as slv
    slv.core.VERBOSE = 0

    out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True)
    compact = os.path.join(os.environ.get('DETECTOR_PATH', '.'), a.config + '.xml')
    scene = os.path.abspath(a.scene) if a.scene else make_scene(compact, out, a.segments)
    prim = os.path.abspath(a.prim) if a.prim else make_prim(compact, out, a.segments)
    if not scene or not prim:
        print('compare_scene_dawnfile: %s not produced, see the logs in %s' % ('scene' if not scene else '.prim', out))
        return 2

    A, B = slv.load_scene(scene), slv.load_scene(prim)
    res = dict(config=a.config, scene=scene, prim=prim, scene_bytes=os.path.getsize(scene),
               prim_bytes=os.path.getsize(prim), volumes_scene=len(A.solids), volumes_dawnfile=len(B.solids),
               vertices=(len(A.V), len(B.V)), facets=(len(A.F), len(B.F)))
    fail = []
    if len(A.solids) != len(B.solids):
        fail.append('count')
    n = min(len(A.solids), len(B.solids))
    anonymous = names = depths = colours = boxes = 0
    worst_box = worst_colour = 0.0
    for s, t in zip(A.solids[:n], B.solids[:n]):
        if t.depth < 0 and t.name.startswith('<'):       # written by DAWNFILE without a #/PVName line
            anonymous += 1
        else:
            names += s.name != t.name
            depths += s.depth != t.depth
        dc = float(np.abs(np.array(s.color) - np.array(t.color)).max())
        db = float(max(np.abs(s.bmin - t.bmin).max(), np.abs(s.bmax - t.bmax).max()))
        worst_colour, worst_box = max(worst_colour, dc), max(worst_box, db)
        colours += dc > 1e-5
        boxes += db > 1e-3
    res.update(anonymous_in_dawnfile=anonymous, name_mismatches=int(names), depth_mismatches=int(depths),
               colour_mismatches=int(colours), box_mismatches=int(boxes), max_colour_diff=worst_colour,
               max_box_diff_mm=worst_box)
    for key, bad in (('names', names), ('depths', depths), ('colours', colours), ('boxes', boxes)):
        if bad:
            fail.append(key)
    if A.bbox is not None and B.bbox is not None:
        res['sphere_diff_mm'] = float(np.abs(A.bbox - B.bbox).max())
        # DAWNFILE prints six significant digits
        if res['sphere_diff_mm'] > 2e-5 * float(np.abs(B.bbox).max()):
            fail.append('sphere')
    else:
        fail.append('sphere')

    print('%s: %d volumes in the scene, %d in the DAWNFILE file (%d without a name there)' % (
        a.config, len(A.solids), len(B.solids), anonymous))
    print('  mismatches: names %d, depths %d, colours %d (max %.1e), bounding boxes %d (max %.1e mm), '
          'bounding sphere %.1e mm' % (names, depths, colours, worst_colour, boxes, worst_box,
                                       res.get('sphere_diff_mm', float('nan'))))
    print('  vertices %d / %d, facets %d / %d, file size %.2f / %.2f MB (scene / DAWNFILE)' % (
        len(A.V), len(B.V), len(A.F), len(B.F), res['scene_bytes'] / 1e6, res['prim_bytes'] / 1e6))

    res['views'] = {}
    if not a.no_render:
        msv = tool('make_solid_views')
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        for label, src in (('scene', scene), ('dawnfile', prim)):
            cmd = [sys.executable, msv, '-c', a.config, '--scene', src, '--views', a.views, '--formats', 'png',
                   '--name', '{view}', '-o', os.path.join(out, label), '--', '-q']
            with open(os.path.join(out, 'render_%s.log' % label), 'w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, check=False)
        for view in [v.strip() for v in a.views.split(',') if v.strip()]:
            fa, fb = (os.path.join(out, lab, view + '.png') for lab in ('scene', 'dawnfile'))
            if not (os.path.isfile(fa) and os.path.isfile(fb)):
                res['views'][view] = None
                fail.append(view + ' (not rendered)')
                continue
            pa, pb = read_png(fa), read_png(fb)
            if pa.shape != pb.shape:
                res['views'][view] = dict(pixels=-1, percent=100.0, size=(pa.shape, pb.shape))
                fail.append(view + ' (canvas size)')
                print('  %-11s canvas sizes differ: %s / %s' % (view, pa.shape[:2], pb.shape[:2]))
                continue
            diff = np.abs(pa - pb).max(-1) > 8.0 / 255.0
            pct = 100.0 * float(diff.mean())
            res['views'][view] = dict(pixels=int(diff.sum()), percent=pct)
            print('  %-11s %7d of %d pixels differ (%.4f %%)' % (view, diff.sum(), diff.size, pct))
            if diff.any():
                write_diff(pa, pb, diff, os.path.join(out, 'diff_%s.png' % view))
            if pct > a.max_diff:
                fail.append(view)
    res['failed'] = fail
    with open(os.path.join(out, 'summary.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print('%s: %s' % (a.config, 'FAILED: ' + ', '.join(fail) if fail else 'scene and DAWNFILE agree'))
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
