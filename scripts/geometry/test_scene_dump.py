#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test of the scene producer on test/test_scene_dump.xml (ctest: scene_dump).

usage: test_scene_dump.py [--tool PATH/solid_scene_dump] [--compact FILE] [--workdir DIR]

The test geometry has one volume placed twice, a solid without mirror symmetry placed once rotated and once
reflected, and a parameterised volume. The script runs solid_scene_dump on it and checks the scene file:
counts, names, depths and parents, a transform, the mesh of the reflected placement, volume IDs and masks,
the bounding sphere, the hash, and that solid_lineview reads the .npz and the .prim to the same solids.
Needs DETECTOR_PATH and the plugin libraries on LD_LIBRARY_PATH (thissolid.sh, or the ctest environment).
Exit status 0 if every check passes.
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
FAILED = []


def check(ok, what):
    print('%s  %s' % ('ok  ' if ok else 'FAIL', what))
    if not ok:
        FAILED.append(what)


def signed_volume(V, F):
    """volume of a closed mesh with outward facets (1-based indices, 0 = no 4th vertex); negative if inside out"""
    vol = 0.0
    for f in F:
        idx = [abs(int(i)) - 1 for i in f if i != 0]
        for k in range(1, len(idx) - 1):
            vol += np.dot(V[idx[0]], np.cross(V[idx[k]], V[idx[k + 1]])) / 6.0
    return vol


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tool', default=shutil.which('solid_scene_dump') or 'solid_scene_dump')
    ap.add_argument('--compact', default=os.path.join(REPO, 'test', 'test_scene_dump.xml'))
    ap.add_argument('--workdir', default=None, help='keep the output files here (default: a temporary directory)')
    a = ap.parse_args()

    try:
        import solid_lineview as slv
    except ImportError:
        sys.path.insert(0, os.path.join(REPO, 'python'))
        import solid_lineview as slv

    work = a.workdir or tempfile.mkdtemp(prefix='scene_dump_test_')
    os.makedirs(work, exist_ok=True)
    npz, prim = os.path.join(work, 'test.npz'), os.path.join(work, 'test.prim')
    r = subprocess.run([a.tool, a.compact, '-o', npz, '--prim', prim], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0 or not os.path.isfile(npz):
        print(r.stdout)
        print('FAIL  solid_scene_dump exit status %d' % r.returncode)
        return 1

    z = np.load(npz)
    meta = json.loads(str(z['meta']))
    name, copy, depth, parent = z['node_name'], z['node_copy'], z['node_depth'], z['node_parent']
    mesh, kind, M = z['node_mesh'], z['node_kind'], z['node_matrix']
    vid, mask = z['node_volume_id'], z['node_volume_mask']

    def node(n, c):
        hit = np.nonzero((name == n) & (copy == c))[0]
        if len(hit) != 1:
            raise SystemExit('FAIL  expected one node %s:%d, found %d' % (n, c, len(hit)))
        return int(hit[0])

    def mesh_of(i):
        vo, fo = z['mesh_vertex_offset'], z['mesh_facet_offset']
        k = mesh[i]
        return z['mesh_vertices'][vo[k]:vo[k + 1]], z['mesh_facets'][fo[k]:fo[k + 1]]

    # ---- counts and structure
    check(meta.get('schema') == 1 and meta.get('config') == 'test_scene_dump', 'meta: schema 1, configuration name')
    check(len(name) == 11 and meta['nodes'] == 11, 'volumes: 11 (world, envelope, 2 pads, 2 wedges, frame, 4 strips)')
    check(int((mesh >= 0).sum()) == 8 and meta['drawn'] == 8, 'drawn: 8 (world, envelope and frame are invisible)')
    check(len(z['mesh_name']) == 4, 'meshes: 4 (pad, wedge, reflected wedge, strip)')
    check(z['mesh_vertices'].dtype == np.float64, 'vertices are float64')
    env = node('SceneTest_envelope_0', 0)
    check(depth[0] == 0 and parent[0] == -1 and depth[env] == 1 and parent[env] == 0, 'world and envelope: depth, parent')

    # ---- one volume placed twice
    p0, p1 = node('SceneTest_pad_0', 0), node('SceneTest_pad_1', 1)
    check(mesh[p0] == mesh[p1] and parent[p0] == env and parent[p1] == env and depth[p0] == 2,
          'pads: two placements share one mesh, mother = envelope, depth 2')
    check(np.allclose(M[p0], [1, 0, 0, 0, 1, 0, 0, 0, 1, -50, 0, 100], atol=1e-9) and
          np.allclose(M[p1, 9:], [50, 0, 100], atol=1e-9), 'pads: world transform (translation -50 / +50, 0, 100 mm)')
    V, F = mesh_of(p0)
    check(len(V) == 8 and len(F) == 6 and np.allclose(np.abs(V).max(0), [10, 20, 5]), 'pad mesh: box 20 x 40 x 10 mm')

    # ---- rotated and reflected placement of the same solid
    w0, w1 = node('SceneTest_wedge_0', 0), node('SceneTest_wedge_1', 1)
    (V0, F0), (V1, F1) = mesh_of(w0), mesh_of(w1)
    W0 = V0 @ M[w0, :9].reshape(3, 3).T + M[w0, 9:]
    W1 = V1 @ M[w1, :9].reshape(3, 3).T + M[w1, 9:]

    def key(P):
        return sorted(tuple(np.round(p, 6) + 0.0) for p in P)

    check(mesh[w0] != mesh[w1], 'wedges: the reflected placement has its own mesh')
    check(key(W1) == key(W0 * [-1, 1, 1]), 'reflected wedge: world vertices are the mirror image (x -> -x)')
    check(key(W1) != key(W0 * [1, 1, 1] - [40, 0, 0]), 'reflected wedge: not a mere translation')
    vol0, vol1 = signed_volume(V0, F0), signed_volume(V1, F1)
    det1 = np.linalg.det(M[w1, :9].reshape(3, 3))
    check(abs(vol0 - 1500.0) < 1e-6, 'wedge mesh: volume 1500 mm3, facets outward')
    check(abs(vol1 * det1 - 1500.0) < 1e-6, 'reflected wedge: volume 1500 mm3, facets outward after the transform')

    # ---- parameterised volume
    strips = np.nonzero(name == 'SceneTest_strip_0')[0]
    frame = node('SceneTest_strip_frame_4', 4)
    check(len(strips) == 4 and set(kind[strips]) == {2} and list(copy[strips]) == [0, 1, 2, 3],
          'strips: 4 copies of a parameterised volume, copy numbers 0..3')
    check(len(set(mesh[strips])) == 1 and set(parent[strips]) == {frame} and set(depth[strips]) == {3},
          'strips: one mesh, mother = frame, depth 3')
    check(np.allclose(M[strips, 10], [-82, -74, -66, -58], atol=1e-9) and np.allclose(M[strips, 11], 100),
          'strips: y = -82, -74, -66, -58 mm, z = 100 mm')

    # ---- identification: system:0:8, layer:8:4, module:12:8, strip:20:8
    check(all(z['node_det_path'][i] == '/world/SceneTest' for i in range(1, 11)) and z['node_det_path'][0] == '',
          'DetElement path /world/SceneTest on all volumes below the world')
    check(all(z['node_readout'][i] == 'SceneTestHits' for i in range(1, 11)), 'readout SceneTestHits')
    check(meta['readouts'] == {'SceneTestHits': 'system:0:8,layer:8:4,module:12:8,strip:20:8'}, 'readout ID fields')
    check(vid[env] == 0x7 and mask[env] == 0xff, 'envelope: volume ID 0x7, mask 0xff (system)')
    check([int(vid[i]) for i in (p0, p1, w0, w1)] == [0x1107, 0x2107, 0x3207, 0x4207] and
          all(mask[i] == 0xfffff for i in (p0, p1, w0, w1)), 'pads and wedges: IDs 0x1107 0x2107 0x3207 0x4207')
    check([int(v) for v in vid[strips]] == [0x307, 0x100307, 0x200307, 0x300307] and all(mask[strips] == 0xff00fff),
          'strips: ID = system 7, layer 3, strip = copy number; mask 0xff00fff')
    sens = z['node_sensitive']
    check(sens[[p0, p1, w0, w1]].all() and sens[strips].all() and not sens[[0, env, frame]].any(), 'sensitive flags')
    # every sensitive volume has a different (ID, mask): a hit maps to exactly one of them
    ids = {(int(vid[i]), int(mask[i])) for i in np.nonzero(sens)[0]}
    check(len(ids) == int(sens.sum()), 'sensitive volumes have distinct volume IDs')

    # ---- colour, bounding sphere, hash
    check(np.allclose(z['node_rgba'][p0], [0.2, 0.4, 0.8, 0.5], atol=1e-6), 'colour and alpha of the vis attributes')
    c, rad = np.array(meta['bounding_sphere'][:3]), meta['bounding_sphere'][3]
    inside = True
    for i in np.nonzero(mesh >= 0)[0]:
        Vi, _ = mesh_of(i)
        Wi = Vi @ M[i, :9].reshape(3, 3).T + M[i, 9:]
        inside &= bool((np.linalg.norm(Wi - c, axis=1) <= rad + 1e-9).all())
    check(inside and 80 < rad < 110, 'bounding sphere encloses the drawn volumes, not the 1 m world (r = %.1f mm)' % rad)
    check(len(meta['geometry_hash']) == 64, 'geometry hash present')
    r2 = subprocess.run([a.tool, a.compact, '-o', os.path.join(work, 'again.npz')], stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, text=True)
    again = np.load(os.path.join(work, 'again.npz')) if r2.returncode == 0 else None
    check(again is not None and json.loads(str(again['meta']))['geometry_hash'] == meta['geometry_hash'] and
          all(np.array_equal(z[k], again[k]) for k in z.files), 'a second run writes the same arrays and hash')

    # ---- solid_lineview reads both files to the same solids
    a_, b_ = slv.load_scene(npz), slv.load_scene(prim)
    same = len(a_.solids) == len(b_.solids) == 8
    for s, t in zip(a_.solids, b_.solids):
        same &= s.name == t.name and s.depth == t.depth and np.allclose(s.color, t.color, atol=1e-6)
        same &= bool(np.allclose(s.bmin, t.bmin, atol=1e-5) and np.allclose(s.bmax, t.bmax, atol=1e-5))
    check(same, 'solid_lineview: .npz and .prim give the same 8 solids (name, depth, colour, bounding box)')
    check(a_.solids[2].name == 'SceneTest_wedge_0:0' and all(s.closed for s in mesh_closed(slv, a_)),
          'solid_lineview: names are <volume>:<copy>, all meshes closed')
    svg = os.path.join(work, 'section.svg')
    rc = slv.main([npz, '--mode', 'section', '--section', 'z=100', '--view', 'beam', '--fill', '-q', '-o', svg])
    check(rc == 0 and os.path.getsize(svg) > 500, 'solid_lineview draws a section of the scene')

    if not a.workdir:
        shutil.rmtree(work, ignore_errors=True)
    print('%d check(s) failed' % len(FAILED) if FAILED else 'all checks passed')
    return 1 if FAILED else 0


def mesh_closed(slv, sc):
    slv.Mesh(sc, 35.0)          # sets Solid.closed
    return sc.solids


if __name__ == '__main__':
    sys.exit(main())
