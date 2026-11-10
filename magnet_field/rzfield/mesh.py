# SPDX-License-Identifier: LGPL-3.0-or-later
"""Conforming quadratic (6-node) triangle mesh of an r-z deck with gmsh.

Every region polygon is a gmsh OCC surface; a boolean fragment splits overlapping polygons so that each piece
follows the region polygons exactly, and each piece takes the properties of the last region covering it
(Automesh rule).  Two outer boundaries:
  mode "poisson": the region-1 polygon of the deck (A = 0 on its edges, as in the Poisson decks);
  mode "open":    a half disk (polygonal, 128 sides) around the magnet with radius open_factor times the radius
                  of the magnet (iron + coils) seen from its centre, A = 0 on the arc.  Region 1 is then only used
                  as the refinement box for the map.
A = 0 also on the axis r = 0.  Sizes (cm, DEFAULT_SIZES, all multiplied by `scale`): h_boundary on iron / coil
edges growing to h_box over grow_dist, h_vertex at polygon vertices (singular corners), 0.7 x the local feature
size at vertices of thin layers / gaps, h_axis along the axis inside the box growing to h_box at axis_dist.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

import numpy as np

from .amdeck import Deck
from .geometry import seg_point_dist

DEFAULT_SIZES = {
    "h_boundary": 1.0,   # cm, at iron / coil boundaries
    "h_box": 4.0,        # cm, max size inside the refinement box (map range)
    "grow_dist": 40.0,   # cm, distance over which h grows from h_boundary to h_box
    "h_min": 0.05,       # cm, lower limit for thin features (gaps, thin layers)
    "far_growth": 0.12,  # open mode: size grows by this fraction of the distance outside the box
    "h_vertex": 0.3,     # cm, size at iron / coil polygon vertices (field singularities at corners)
    "h_axis": 1.0,       # cm, size on the axis inside the refinement box (smooth on-axis Bz = 2 dA/dr)
    "axis_dist": 20.0,   # cm, distance from the axis over which h grows from h_axis to h_box
}


@dataclass
class Mesh:
    nodes: np.ndarray            # (nn, 2) r, z in metres
    tri: np.ndarray              # (ne, 6) node indices: 3 vertices, then mid-edges 01, 12, 20
    zone: np.ndarray             # (ne,) zone index of each element
    zones: list                  # per zone: dict(region, mat, cur, name, area_cm2, J)
    dirichlet: np.ndarray        # node indices with A = 0
    info: dict = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def save(self, path: str):
        np.savez_compressed(path, nodes=self.nodes, tri=self.tri, zone=self.zone, dirichlet=self.dirichlet,
                            zones=json.dumps(self.zones), info=json.dumps(self.info))

    @classmethod
    def load(cls, path: str) -> "Mesh":
        d = np.load(path, allow_pickle=False)
        return cls(nodes=d["nodes"], tri=d["tri"], zone=d["zone"], dirichlet=d["dirichlet"],
                   zones=json.loads(str(d["zones"])), info=json.loads(str(d["info"])))


def _local_feature_size(pts, polys):
    """For every point: distance to the nearest polygon segment that does not pass through it."""
    A = np.vstack([P if len(P) > 2 else P[:1] for P in polys])
    Bv = np.vstack([np.roll(P, -1, axis=0) if len(P) > 2 else P[1:] for P in polys])
    d = seg_point_dist(pts, A, Bv)
    d[d < 1e-7] = np.inf                       # segments through the point itself
    return d.min(axis=1)


def build_mesh(deck: Deck, mode: str = "poisson", scale: float = 1.0, sizes: dict | None = None,
               open_factor: float = 10.0, refine_box=None, threads: int = 1, verbose: bool = False) -> Mesh:
    import gmsh  # only needed here

    if mode not in ("poisson", "open"):
        raise ValueError("mode must be 'poisson' or 'open'")
    S = dict(DEFAULT_SIZES)
    if sizes:
        S.update(sizes)
    hb = S["h_boundary"] * scale
    hbox = S["h_box"] * scale
    hmin = S["h_min"]
    t0 = time.time()

    regions = deck.regions
    box = deck.bbox("box") if refine_box is None else tuple(refine_box)
    mag = deck.bbox("magnet")
    zc = 0.5 * (mag[2] + mag[3])
    mag_pts = np.vstack([r.points for r in regions[1:] if r.mat >= 2 or r.cur != 0.0])
    r_mag = float(np.max(np.hypot(mag_pts[:, 0], mag_pts[:, 1] - zc)))
    R_out = open_factor * r_mag

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("General.Verbosity", 3 if verbose else 1)
        gmsh.option.setNumber("General.NumThreads", max(1, int(threads)))
        gmsh.model.add("rz")
        occ = gmsh.model.occ

        def add_poly(P):
            pt = [occ.addPoint(float(x), float(y), 0.0) for x, y in P]
            ln = [occ.addLine(pt[i], pt[(i + 1) % len(pt)]) for i in range(len(pt))]
            return occ.addPlaneSurface([occ.addCurveLoop(ln)])

        inputs = []          # (dimtag, region index); region index 0 = open domain
        if mode == "poisson":
            inputs.append(((2, add_poly(regions[0].points)), 1))
        else:
            nseg = 128
            th = np.linspace(-0.5 * np.pi, 0.5 * np.pi, nseg + 1)
            arc = np.column_stack([R_out * np.cos(th), zc + R_out * np.sin(th)])
            arc[0, 0] = arc[-1, 0] = 0.0
            inputs.append(((2, add_poly(arc[::-1])), 0))
        for r in regions[1:]:
            inputs.append(((2, add_poly(r.points)), r.index))
        out, out_map = occ.fragment([dt for dt, _ in inputs], [])
        occ.synchronize()

        owner = {}            # fragment surface tag -> list of region indices covering it
        for (dt, ridx), children in zip(inputs, out_map):
            for c in children:
                if c[0] == 2:
                    owner.setdefault(c[1], []).append(ridx)
        domain_idx = 1 if mode == "poisson" else 0
        keep, drop = {}, []
        for tag, owners in owner.items():
            if domain_idx not in owners:
                drop.append(tag)      # outside the problem domain
                continue
            regs = [o for o in owners if o >= 2]
            keep[tag] = max(regs) if regs else 1
        if drop:
            occ.remove([(2, t) for t in drop], recursive=True)
            occ.synchronize()

        # ---- mesh sizes
        all_curves = gmsh.model.getEntities(1)
        outer = gmsh.model.getBoundary([(2, t) for t in keep], combined=True, oriented=False, recursive=False)
        outer_curves = {abs(c[1]) for c in outer}
        inner_curves = [c[1] for c in all_curves if c[1] not in outer_curves]
        # thin features: size at the region vertices from the local feature size (layers, gaps)
        mag_polys = [r.points for r in regions[1:]]
        other = [regions[0].points] if mode == "poisson" else [np.array([[0.0, zc - R_out], [0.0, zc + R_out]])]
        vtx = np.vstack(mag_polys)
        lfs = _local_feature_size(vtx, mag_polys + other)
        h_far = max(hbox, R_out / 10.0) if mode == "open" else hbox
        for (_, ptag) in gmsh.model.getEntities(0):
            xyz = gmsh.model.getValue(0, ptag, [])
            d = np.hypot(vtx[:, 0] - xyz[0], vtx[:, 1] - xyz[1])
            j = int(np.argmin(d))
            h_pt = min(S["h_vertex"] * scale, max(hmin, 0.7 * lfs[j])) if d[j] < 1e-6 else h_far
            gmsh.model.mesh.setSize([(0, ptag)], float(h_pt))
        fields = []
        if inner_curves:
            fd = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(fd, "CurvesList", inner_curves)
            gmsh.model.mesh.field.setNumber(fd, "Sampling", 500)
            ft = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(ft, "InField", fd)
            gmsh.model.mesh.field.setNumber(ft, "SizeMin", hb)
            gmsh.model.mesh.field.setNumber(ft, "SizeMax", hbox)
            gmsh.model.mesh.field.setNumber(ft, "DistMin", 0.5 * hb)
            gmsh.model.mesh.field.setNumber(ft, "DistMax", S["grow_dist"])
            gmsh.model.mesh.field.setNumber(ft, "StopAtDistMax", 1)
            fields.append(ft)
        if mode == "poisson":
            fc = gmsh.model.mesh.field.add("Constant")
            gmsh.model.mesh.field.setNumber(fc, "VIn", hbox)
            fields.append(fc)
        else:
            # h_box inside the refinement box, growing linearly with the distance outside it
            x0, x1, y0, y1 = box
            expr = ("Min(%g, %g + %g*Sqrt((Max(0,Max(%g-x,x-%g)))^2 + (Max(0,Max(%g-y,y-%g)))^2))"
                    % (h_far, hbox, S["far_growth"], x0, x1, y0, y1))
            fm = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(fm, "F", expr)
            fields.append(fm)
        if S.get("h_axis"):
            # finer elements along the axis inside the box: Bz(0, z) is a derivative (2 dA/dr) of the first row
            x0, x1, y0, y1 = box
            ha = S["h_axis"] * scale
            # (no constraint beyond r = axis_dist or outside the box in z)
            expr = ("%g + %g*Min(1, x/%g) + 1000*(Max(0, x-%g) + Max(0, y-%g) + Max(0, %g-y))"
                    % (ha, max(hbox - ha, 0.0), S["axis_dist"], S["axis_dist"], y1, y0))
            fa = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(fa, "F", expr)
            fields.append(fa)
        fmin = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_far)
        gmsh.option.setNumber("Mesh.MeshSizeMin", hmin)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
        gmsh.model.mesh.generate(2)

        # ---- extract
        ntags, coords, _ = gmsh.model.mesh.getNodes()
        coords = coords.reshape(-1, 3)[:, :2]
        tag2idx = np.full(int(ntags.max()) + 1, -1, dtype=np.int64)
        tag2idx[ntags.astype(np.int64)] = np.arange(len(ntags))
        tri_list, zone_list = [], []
        zones = []
        zone_of_region = {}
        for stag, ridx in sorted(keep.items()):
            et, _, en = gmsh.model.mesh.getElements(2, stag)
            for t, conn in zip(et, en):
                if t != 9:
                    raise RuntimeError("unexpected element type %d (want 6-node triangles)" % t)
                conn = tag2idx[np.asarray(conn, dtype=np.int64)].reshape(-1, 6)
                if ridx not in zone_of_region:
                    reg = regions[ridx - 1] if ridx >= 1 else None
                    zone_of_region[ridx] = len(zones)
                    zones.append({"region": int(ridx), "mat": int(reg.mat) if reg else 1,
                                  "cur": float(reg.cur) if reg else 0.0,
                                  "name": (reg.name if reg else "open domain"),
                                  "polygon_area_cm2": float(reg.area_cm2) if reg else 0.0})
                tri_list.append(conn)
                zone_list.append(np.full(len(conn), zone_of_region[ridx], dtype=np.int32))
        dnodes = []
        for c in outer_curves:
            nt, _, _ = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)
            dnodes.append(tag2idx[nt.astype(np.int64)])
        n_curves = len(all_curves)
    finally:
        gmsh.finalize()

    tri = np.vstack(tri_list)
    zone = np.concatenate(zone_list)
    # compact node numbering (drop unused nodes, e.g. from removed fragments)
    used = np.unique(tri)
    remap = np.full(len(coords), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    nodes_cm = coords[used].copy()
    tri = remap[tri]
    dirichlet = np.unique(remap[np.concatenate(dnodes)])
    dirichlet = dirichlet[dirichlet >= 0]
    # straight-sided (affine) elements: put mid-edge nodes exactly at the edge midpoints
    for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
        nodes_cm[tri[:, 3 + k]] = 0.5 * (nodes_cm[tri[:, a]] + nodes_cm[tri[:, b]])
    # orient all triangles counter-clockwise
    v = nodes_cm[tri[:, :3]]
    det = (v[:, 1, 0] - v[:, 0, 0]) * (v[:, 2, 1] - v[:, 0, 1]) - (v[:, 2, 0] - v[:, 0, 0]) * (v[:, 1, 1] - v[:, 0, 1])
    neg = det < 0
    if np.any(neg):
        tri[neg] = tri[neg][:, [0, 2, 1, 5, 4, 3]]
    # axis nodes are always Dirichlet
    axis = np.where(np.abs(nodes_cm[:, 0]) < 1e-9)[0]
    dirichlet = np.union1d(dirichlet, axis)
    # assigned area and current density per zone
    v = nodes_cm[tri[:, :3]]
    area = 0.5 * np.abs((v[:, 1, 0] - v[:, 0, 0]) * (v[:, 2, 1] - v[:, 0, 1])
                        - (v[:, 2, 0] - v[:, 0, 0]) * (v[:, 1, 1] - v[:, 0, 1]))
    zarea = np.bincount(zone, weights=area, minlength=len(zones))
    for z, zd in enumerate(zones):
        zd["area_cm2"] = float(zarea[z])
        zd["J"] = zd["cur"] / (zarea[z] * 1e-4) if zd["cur"] != 0.0 else 0.0   # A/m^2, total current kept
    info = {"mode": mode, "scale": scale, "sizes_cm": S, "open_factor": open_factor if mode == "open" else None,
            "R_out_cm": R_out if mode == "open" else None, "z_center_cm": zc, "r_magnet_cm": r_mag,
            "refine_box_cm": list(box), "n_elements": int(len(tri)), "n_nodes": int(len(nodes_cm)),
            "n_curves": n_curves, "deck": deck.path, "title": deck.title, "mesh_seconds": time.time() - t0}
    return Mesh(nodes=nodes_cm * 1e-2, tri=tri.astype(np.int64), zone=zone, zones=zones,
                dirichlet=dirichlet.astype(np.int64), info=info)
