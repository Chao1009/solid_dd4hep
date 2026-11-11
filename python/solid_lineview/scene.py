# SPDX-License-Identifier: LGPL-3.0-or-later
"""Reader of scene.npz files (written by solid_scene_dump or by ddsim with SOLID_SCENE=file).

A scene file holds the Geant4 polyhedra of a geometry once per distinct solid, and one record per
placed volume (schema version 1; all lengths in mm):

    meta                       JSON text: schema, config, mode, segments, bounding_sphere (x, y, z, r of the
                               drawn volumes), geometry_hash, versions, readouts {name: ID fields}
    mesh_vertices              (Nv, 3) float64, frame of the solid
    mesh_facets                (Nf, 4) int32: 1-based vertex index, negative when Geant4 flags the edge that
                               starts there as invisible, 0 = triangle (no 4th vertex)
    mesh_vertex_offset         (M + 1,) int64: mesh k owns vertices offset[k]:offset[k + 1]
    mesh_facet_offset          (M + 1,) int64, same for facets
    mesh_name, mesh_type       (M,) name and Geant4 type of the solid
    node_name, node_copy       (N,) physical-volume name and copy number, depth-first order
    node_depth, node_parent    (N,) depth below the world (0) and index of the mother node (-1 for the world)
    node_mesh                  (N,) mesh index, -1 = not drawn (invisible)
    node_kind                  (N,) 0 placement, 1 replica, 2 parameterised
    node_matrix                (N, 12) local -> world: rotation (3 x 3, row-major), then translation
    node_rgba                  (N, 4) colour and alpha of the vis attributes
    node_visible, node_daughters_invisible, node_wireframe   (N,) Geant4 vis attribute flags
    node_drawing_style         (N,) DD4hep VisAttr drawingStyle: 0 unknown, 1 solid, 2 wireframe
    node_det_path              (N,) path of the deepest DetElement at or above the volume
    node_readout               (N,) readout (hit collection) name, '' if none
    node_volume_id, node_volume_mask   (N,) uint64: volume ID bits fixed by the placements down to this
                               volume and their mask; a hit lies in the volume if (cellID & mask) == id
    node_sensitive             (N,) the volume is sensitive

load_scene_npz() turns this into the PrimScene the renderer works on, so a scene file and a DAWNFILE .prim
of the same geometry give the same drawing.
"""
import json
import zipfile

import numpy as np

from . import core

SCHEMA_VERSION = 1
_REQUIRED = ("mesh_vertices", "mesh_facets", "mesh_vertex_offset", "mesh_facet_offset", "node_name",
             "node_copy", "node_depth", "node_mesh", "node_matrix", "node_rgba")


def is_scene_file(path):
    """True for a .npz scene (a zip archive), whatever the file is called."""
    try:
        with open(path, "rb") as fh:
            return fh.read(4) == b"PK\x03\x04"
    except OSError:
        return False


class SceneData:
    """The arrays of a scene.npz, as numpy arrays under their file names, and the parsed `meta`."""

    def __init__(self, path):
        try:
            with np.load(path, allow_pickle=False) as z:
                self.arrays = {k: z[k] for k in z.files}
        except (OSError, ValueError, zipfile.BadZipFile) as exc:
            raise SystemExit("cannot read scene file %s: %s" % (path, exc))
        missing = [k for k in _REQUIRED if k not in self.arrays]
        if missing:
            raise SystemExit("%s is not a geometry scene: no array %s" % (path, ", ".join(missing)))
        try:
            self.meta = json.loads(str(self.arrays["meta"])) if "meta" in self.arrays else {}
        except ValueError:
            self.meta = {}
        if self.meta.get("schema", SCHEMA_VERSION) > SCHEMA_VERSION:
            core.log("warning: %s has scene schema %s, this reader knows %d" % (
                path, self.meta.get("schema"), SCHEMA_VERSION))

    def __getitem__(self, key):
        return self.arrays[key]

    def __contains__(self, key):
        return key in self.arrays

    def get(self, key, default=None):
        return self.arrays.get(key, default)

    @property
    def drawn(self):
        """indices of the drawn nodes, in file (depth-first) order"""
        return np.nonzero(self.arrays["node_mesh"] >= 0)[0]

    def mesh(self, k):
        """(V, F) of mesh k as stored: local vertices and signed 1-based facets"""
        vo, fo = self.arrays["mesh_vertex_offset"], self.arrays["mesh_facet_offset"]
        return self.arrays["mesh_vertices"][vo[k]:vo[k + 1]], self.arrays["mesh_facets"][fo[k]:fo[k + 1]]


def load_scene_npz(path, drawing_style=False):
    """Read a scene.npz into a PrimScene.

    Each mesh is welded and cleaned once, then placed with the matrix of every drawn node.  All edges are
    taken as visible, as Geant4's DAWNFILE driver writes them.  With drawing_style=True a volume whose
    DD4hep VisAttr says drawingStyle="wireframe" is marked as a wireframe volume (DDG4 does not pass that
    setting to Geant4, so by default the drawing equals the DAWNFILE one).

    The scene keeps the file content as sc.data (a SceneData) and sc.node (drawn solid -> node index), so
    that callers can look up DetElement paths and volume IDs of the solids.
    """
    data = SceneData(path)
    MV = np.asarray(data["mesh_vertices"], np.float64).reshape(-1, 3)
    cache = {}

    def mesh(k):
        hit = cache.get(k)
        if hit is None:
            vo, fo = data["mesh_vertex_offset"], data["mesh_facet_offset"]
            V = MV[vo[k]:vo[k + 1]]
            F = np.abs(data["mesh_facets"][fo[k]:fo[k + 1]]).astype(np.int64) - 1      # 0 (absent) -> -1
            V, F = core._weld(V, F)
            F, Fvis = core._clean_faces(F, np.ones(F.shape, bool))
            hit = cache[k] = (V, F, Fvis)
        return hit

    sc = core.PrimScene()
    names, copy, depth = data["node_name"], data["node_copy"], data["node_depth"]
    rgba, M, nmesh = data["node_rgba"], np.asarray(data["node_matrix"], np.float64), data["node_mesh"]
    wire = np.asarray(data.get("node_wireframe", np.zeros(len(names), bool)), bool)
    if drawing_style and "node_drawing_style" in data:
        wire = wire | (data["node_drawing_style"] == 2)
    bV, bF, bE, node = [], [], [], []
    nv = nf = 0
    for i in data.drawn:
        V, F, Fvis = mesh(int(nmesh[i]))
        if len(F) == 0:
            continue
        R = M[i, :9].reshape(3, 3)
        W = V @ R.T + M[i, 9:]
        if np.linalg.det(R) < 0:                 # a reflection in the matrix: keep the normals outward
            F = core._flip_faces(F)
        s = core.Solid()
        s.index = len(sc.solids)
        s.name = "%s:%d" % (names[i], copy[i])
        s.depth = int(depth[i])
        s.color = tuple(min(1.0, max(0.0, float(x))) for x in rgba[i][:3])
        s.wire = bool(wire[i])
        s.kind = "Polyhedron"
        s.line = int(i)
        s.v0, s.v1, s.f0, s.f1 = nv, nv + len(W), nf, nf + len(F)
        s.bmin, s.bmax = W.min(0), W.max(0)
        bV.append(W)
        bF.append(np.where(F >= 0, F + nv, -1))
        bE.append(Fvis)
        nv += len(W)
        nf += len(F)
        sc.solids.append(s)
        node.append(int(i))
        sc.counts["Polyhedron"] += 1
    if bV:
        sc.V = np.vstack(bV)
        sc.F = np.vstack(bF).astype(np.int64)
        sc.Fvis = np.vstack(bE)
        sc.fsolid = np.repeat(np.arange(len(sc.solids), dtype=np.int32), [s.f1 - s.f0 for s in sc.solids])
    sphere = data.meta.get("bounding_sphere")
    if sphere and len(sphere) == 4 and sphere[3] > 0:
        c, r = np.array(sphere[:3], float), float(sphere[3])
        sc.bbox = np.concatenate([c - r, c + r])   # what DAWNFILE writes as /BoundingBox
    core._build_hierarchy(sc)
    sc.data = data
    sc.node = np.array(node, np.int64)
    return sc


def write_prim(path_npz, path_prim):
    """Write the drawn volumes of a scene.npz as a DAWN g4.prim file (for dawn and dawncut).

    One /Polyhedron per drawn volume with its name, depth and colour; all edges visible.  The same file
    `solid_scene_dump --prim` writes.
    """
    data = SceneData(path_npz)
    names, copy, depth = data["node_name"], data["node_copy"], data["node_depth"]
    rgba, M, nmesh = data["node_rgba"], np.asarray(data["node_matrix"], np.float64), data["node_mesh"]
    wire = np.asarray(data.get("node_wireframe", np.zeros(len(names), bool)), bool)
    sphere = data.meta.get("bounding_sphere") or [0.0, 0.0, 0.0, 0.0]
    text = {}

    def fmt(v):
        return " ".join("%.9g" % x for x in v)

    with open(path_prim, "w") as out:
        out.write("##G4.PRIM-FORMAT-2.4\n")
        out.write("/BoundingBox %s\n" % fmt([sphere[k] - sphere[3] for k in range(3)] +
                                           [sphere[k] + sphere[3] for k in range(3)]))
        out.write("!SetCamera\n!OpenDevice\n!BeginModeling\n/Ndiv %d\n" % int(data.meta.get("segments", 72)))
        for i in data.drawn:
            k = int(nmesh[i])
            R, T = M[i, :9].reshape(3, 3), M[i, 9:]
            out.write("#--------------------\n#/PVName%s%s:%d\n" % (" " * (5 + 5 * int(depth[i])), names[i], copy[i]))
            out.write("/ColorRGB %s\n/ForceWireframe  %d\n" % (fmt(rgba[i][:3]), int(wire[i])))
            V, F = data.mesh(k)
            if np.linalg.det(R) > 0:
                if k not in text:
                    text[k] = "".join("/Vertex %s\n" % fmt(v) for v in V) + "".join(
                        "/Facet %s\n" % " ".join(str(abs(int(x))) for x in f if x != 0) for f in F)
                # /BaseVector: images of the local x and y axes; DAWN takes z = x cross y
                out.write("/Origin %s\n/BaseVector %s %s\n/Polyhedron\n%s/EndPolyhedron\n" % (
                    fmt(T), fmt(R[:, 0]), fmt(R[:, 1]), text[k]))
            else:
                # a reflection cannot be written as /BaseVector: world coordinates, reversed winding
                out.write("/Origin 0 0 0\n/BaseVector 1 0 0 0 1 0\n/Polyhedron\n")
                out.write("".join("/Vertex %s\n" % fmt(v) for v in V @ R.T + T))
                out.write("".join("/Facet %s\n" % " ".join(str(abs(int(x))) for x in f[::-1] if x != 0) for f in F))
                out.write("/EndPolyhedron\n")
        out.write("!EndModeling\n!DrawAll\n!CloseDevice\n")
