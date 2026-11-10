#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check the geometry scene of each configuration against the stored reference (scene_hashes.json).

usage (inside eic-shell / solid-shell, after `source thissolid.sh`):
    check_scene_hash.py [CONFIG ...]            check the given configurations (default: all in the reference file)
    check_scene_hash.py --update [CONFIG ...]   store the current hashes (after an intended geometry change)
    check_scene_hash.py --scene FILE CONFIG     check an existing scene.npz instead of running solid_scene_dump

solid_scene_dump stores a SHA-256 of the scene content in the file: every mesh (vertices rounded to 0.1 um,
facets), every placement (name, copy number, depth, mother, mesh, transform, colour, flags, DetElement path,
readout, volume ID and mask) and the number of segments per circle. File names and software versions are not
hashed. A changed hash means that the drawn geometry, its colours or its identification changed: intended after
a geometry change (run --update and commit the file), a regression otherwise. The volume, drawn-volume and mesh
counts are stored next to the hash to make a difference easier to read.

The hash can also change with a new DD4hep or Geant4 version (volume names of assemblies, polyhedron code);
the versions that produced the reference are recorded in the file.

Exit status: 0 = all equal (or updated), 1 = a difference or a missing reference, 2 = a scene was not produced.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_FILE = os.path.join(HERE, 'scene_hashes.json')
KEYS = ('geometry_hash', 'nodes', 'drawn', 'meshes', 'segments')


def scene_meta(path):
    with np.load(path, allow_pickle=False) as z:
        return json.loads(str(z['meta']))


def dump(config, workdir):
    tool = shutil.which('solid_scene_dump') or os.path.join(REPO, 'build', 'bin', 'solid_scene_dump')
    compact = os.path.join(os.environ.get('DETECTOR_PATH', '.'), config + '.xml')
    scene = os.path.join(workdir, config + '.npz')
    r = subprocess.run([tool, compact, '-o', scene], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0 or not os.path.isfile(scene):
        print(r.stdout[-2000:])
        return None
    return scene


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('configs', nargs='*', help='configurations, e.g. solid_JPsi_lite (default: all stored ones)')
    ap.add_argument('--update', action='store_true', help='write the current values to the reference file')
    ap.add_argument('--scene', default=None, help='existing scene.npz of the one configuration given')
    ap.add_argument('--file', default=DEFAULT_FILE, help='reference file (default: scene_hashes.json next to this script)')
    a = ap.parse_args()

    ref = {'configurations': {}}
    if os.path.isfile(a.file):
        with open(a.file) as f:
            ref = json.load(f)
    configs = a.configs or sorted(ref['configurations'])
    if not configs:
        sys.exit('check_scene_hash: no configuration given and none stored in %s' % a.file)
    if a.scene and len(configs) != 1:
        sys.exit('check_scene_hash: --scene needs exactly one configuration')

    work = tempfile.mkdtemp(prefix='scene_hash_')
    status = 0
    versions = None
    try:
        for config in configs:
            scene = a.scene or dump(config, work)
            if not scene:
                print('%-18s no scene produced' % config)
                status = max(status, 2)
                continue
            meta = scene_meta(scene)
            now = {k: meta.get(k) for k in KEYS}
            versions = meta.get('versions', versions)
            old = ref['configurations'].get(config)
            if a.update:
                ref['configurations'][config] = now
                print('%-18s %s  %s' % (config, now['geometry_hash'][:16], 'unchanged' if old == now else 'updated'))
            elif old is None:
                print('%-18s %s  no reference (run with --update)' % (config, now['geometry_hash'][:16]))
                status = max(status, 1)
            elif old == now:
                print('%-18s %s  ok (%d volumes, %d drawn, %d meshes)' % (
                    config, now['geometry_hash'][:16], now['nodes'], now['drawn'], now['meshes']))
            else:
                diff = ', '.join('%s %s -> %s' % (k, str(old.get(k))[:16], str(now[k])[:16])
                                 for k in KEYS if old.get(k) != now[k])
                print('%-18s CHANGED: %s' % (config, diff))
                status = max(status, 1)
        if a.update and status == 0:
            if versions:
                ref['versions'] = versions
            ref['comment'] = ('Reference scene hashes, written by scripts/geometry/check_scene_hash.py --update. '
                              'Update after an intended change of the geometry, its colours or its volume IDs.')
            with open(a.file, 'w') as f:
                json.dump(ref, f, indent=1, sort_keys=True)
                f.write('\n')
            print('wrote', a.file)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return status


if __name__ == '__main__':
    sys.exit(main())
