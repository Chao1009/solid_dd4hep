# SPDX-License-Identifier: LGPL-3.0-or-later
"""Solution container (mesh + nodal A) with B(r, z) evaluation; needs numpy and scipy only (no gmsh)."""
from __future__ import annotations

import json

import numpy as np
from scipy.spatial import cKDTree

from .fem import p2_basis


class RZSolution:
    """Quadratic-triangle solution A_phi on an r-z mesh (metres, Wb/m)."""

    def __init__(self, nodes, tri, A, zone=None, zones=None, info=None):
        self.nodes = np.asarray(nodes, dtype=float)
        self.tri = np.asarray(tri, dtype=np.int64)
        self.A = np.asarray(A, dtype=float)
        self.zone = None if zone is None else np.asarray(zone)
        self.zones = zones or []
        self.info = info or {}
        v = self.nodes[self.tri[:, :3]]
        self._x0 = v[:, 0, :]
        e1 = v[:, 1, :] - v[:, 0, :]
        e2 = v[:, 2, :] - v[:, 0, :]
        det = e1[:, 0] * e2[:, 1] - e2[:, 0] * e1[:, 1]
        # inverse affine map: (l1, l2) = M (p - x0); M = [dxi/dr, dxi/dz, deta/dr, deta/dz]
        self._inv = np.stack([e2[:, 1] / det, -e2[:, 0] / det, -e1[:, 1] / det, e1[:, 0] / det], -1)
        self._centroid = v.mean(axis=1)
        self._tree = cKDTree(self._centroid)
        self._hmax = float(np.max(np.linalg.norm(v - self._centroid[:, None, :], axis=2)))
        on_axis = np.abs(v[:, :, 0]) < 1e-12
        self._axis_edge = on_axis.sum(axis=1) >= 2
        self._rec = None
        self._air = None

    # ------------------------------------------------------------------------------------------------- I/O
    def save(self, path: str):
        np.savez_compressed(path, nodes=self.nodes, tri=self.tri, A=self.A,
                            zone=self.zone if self.zone is not None else np.zeros(0, int),
                            zones=json.dumps(self.zones), info=json.dumps(self.info, default=float))

    @classmethod
    def load(cls, path: str) -> "RZSolution":
        d = np.load(path, allow_pickle=False)
        zone = d["zone"] if d["zone"].size else None
        return cls(d["nodes"], d["tri"], d["A"], zone=zone, zones=json.loads(str(d["zones"])),
                   info=json.loads(str(d["info"])))

    # ------------------------------------------------------------------------------------------ point location
    def air_elements(self):
        """(ne,) bool: elements of the air class, i.e. mat = 1 (air and coils); all True without zone data."""
        if self._air is None:
            if self.zone is None or not self.zones:
                self._air = np.ones(len(self.tri), dtype=bool)
            else:
                self._air = np.array([int(zd.get("mat", 1)) == 1 for zd in self.zones], dtype=bool)[self.zone]
        return self._air

    def locate(self, r, z, k: int = 12, prefer_air: bool = True):
        """Element index and reference coordinates (l1, l2) of points (metres); -1 if outside the mesh.

        A point on an element edge belongs to several elements. On the axis an element with an edge on the axis is
        taken (A/r is regular there). On an iron/air interface, prefer_air (default) takes the air element, so that
        the field at such a point is the air-side value whatever the mesh (otherwise the side would depend on the
        mesh, i.e. on the gmsh version); the tolerance is 1e-7 in barycentric coordinates."""
        p = np.column_stack([np.ravel(r), np.ravel(z)]).astype(float)
        n = len(p)
        elem = np.full(n, -1, dtype=np.int64)
        l1 = np.zeros(n)
        l2 = np.zeros(n)
        todo = np.arange(n)
        on_axis = np.abs(p[:, 0]) < 1e-12
        air = self.air_elements() if prefer_air else None
        for kk in (k, 4 * k, 32 * k):
            if len(todo) == 0:
                break
            kk = min(kk, len(self.tri))
            _, cand = self._tree.query(p[todo], k=kk)
            cand = cand.reshape(len(todo), kk)
            d = p[todo][:, None, :] - self._x0[cand]                      # (m, kk, 2)
            m = self._inv[cand]
            a1 = m[..., 0] * d[..., 0] + m[..., 1] * d[..., 1]
            a2 = m[..., 2] * d[..., 0] + m[..., 3] * d[..., 1]
            a0 = 1.0 - a1 - a2
            score = np.minimum(np.minimum(a0, a1), a2)
            # on the axis use an element with an edge on the axis (A/r is regular there)
            ax = on_axis[todo]
            if np.any(ax):
                score[ax] = np.where(self._axis_edge[cand[ax]], score[ax], score[ax] - 10.0)
            if air is not None:
                score = score + np.where(air[cand], 1e-7, 0.0)
            best = np.argmax(score, axis=1)
            ok = score[np.arange(len(todo)), best] > -1e-9
            sel = todo[ok]
            elem[sel] = cand[ok, best[ok]]
            l1[sel] = a1[ok, best[ok]]
            l2[sel] = a2[ok, best[ok]]
            todo = todo[~ok]
        return elem, l1, l2

    # ---------------------------------------------------------------------------------------------- fields
    def B(self, r, z, smooth: bool = True, chunk: int = 200000, return_air: bool = False):
        """(Br, Bz) in tesla at points r, z in metres (arrays of any equal shape); NaN outside the mesh.
        return_air=True also returns a boolean array: the point lies in air or a coil (mat = 1), not in iron (points
        on an iron/air interface count as air and get the air-side field, see locate()).

        smooth=True (default): recovered field, i.e. B averaged at the quadratic nodes over the elements of the
        same material that share the node, then interpolated with the quadratic shape functions (continuous
        inside each material, discontinuous across iron/air interfaces as it should be).
        smooth=False: the raw finite-element field curl(A) of the element containing the point (linear in each
        element, small jumps between elements)."""
        r = np.asarray(r, dtype=float)
        z = np.asarray(z, dtype=float)
        shape = np.broadcast(r, z).shape
        rr = np.broadcast_to(r, shape).ravel()
        zz = np.broadcast_to(z, shape).ravel()
        Br = np.full(rr.size, np.nan)
        Bz = np.full(rr.size, np.nan)
        in_air = np.zeros(rr.size, dtype=bool)
        for s in range(0, rr.size, chunk):
            sl = slice(s, min(s + chunk, rr.size))
            e, l1, l2 = self.locate(rr[sl], zz[sl])
            ok = e >= 0
            if not np.any(ok):
                continue
            e = e[ok]
            N, dxi, deta = p2_basis(l1[ok], l2[ok])                       # (m, 6)
            idx = np.arange(sl.start, sl.stop)[ok]
            in_air[idx] = self.air_elements()[e]
            if smooth:
                rec = self.recovered_nodal_field()[e]                     # (m, 6, 2)
                Br[idx] = np.sum(N * rec[:, :, 0], axis=1)
                Bz[idx] = np.sum(N * rec[:, :, 1], axis=1)
                continue
            inv = self._inv[e]                                            # dxi/dr, dxi/dz, deta/dr, deta/dz
            dNdr = inv[:, 0:1] * dxi + inv[:, 2:3] * deta
            dNdz = inv[:, 1:2] * dxi + inv[:, 3:4] * deta
            Ae = self.A[self.tri[e]]
            rp = rr[sl][ok]
            dAdr = np.sum(Ae * dNdr, axis=1)
            Aval = np.sum(Ae * N, axis=1)
            small = rp < 1e-12
            A_over_r = np.where(small, dAdr, Aval / np.where(small, 1.0, rp))
            Br[idx] = np.where(small, 0.0, -np.sum(Ae * dNdz, axis=1))
            Bz[idx] = dAdr + A_over_r
        if return_air:
            return Br.reshape(shape), Bz.reshape(shape), in_air.reshape(shape)
        return Br.reshape(shape), Bz.reshape(shape)

    def B_cm_gauss(self, r_cm, z_cm, smooth: bool = True, return_air: bool = False):
        """(Br, Bz) in gauss at r, z in cm (plus the in-air flags with return_air=True, see B())."""
        out = self.B(np.asarray(r_cm, dtype=float) * 1e-2, np.asarray(z_cm, dtype=float) * 1e-2, smooth=smooth,
                     return_air=return_air)
        return (out[0] * 1e4, out[1] * 1e4) + tuple(out[2:])

    def material_class(self):
        """Element class for the field recovery: one class per distinct material (air and coils share one)."""
        if self.zone is None or not self.zones:
            return np.zeros(len(self.tri), dtype=np.int64)
        keys = []
        for zd in self.zones:
            k = zd.get("material") or ("mat%d" % int(zd.get("mat", 1)))
            if int(zd.get("mat", 1)) == 1:
                k = "air"
            keys.append(k)
        ids = {k: i for i, k in enumerate(dict.fromkeys(keys))}
        return np.array([ids[k] for k in keys], dtype=np.int64)[self.zone]

    def recovered_nodal_field(self):
        """(ne, 6, 2) nodal (Br, Bz) in tesla, averaged per node and material class (computed once)."""
        if self._rec is not None:
            return self._rec
        tri = self.tri
        cls = self.material_class()
        ncls = int(cls.max()) + 1
        ref = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5]])
        _, dxi, deta = p2_basis(ref[:, 0], ref[:, 1])                      # (6 nodes, 6 basis)
        Ae = self.A[tri]
        gxi = Ae @ dxi.T
        geta = Ae @ deta.T
        inv = self._inv
        dAdr = inv[:, 0:1] * gxi + inv[:, 2:3] * geta
        dAdz = inv[:, 1:2] * gxi + inv[:, 3:4] * geta
        rn = self.nodes[tri, 0]
        axis = rn < 1e-12
        Bz = np.where(axis, 2.0 * dAdr, dAdr + Ae / np.where(axis, 1.0, rn))
        Br = np.where(axis, 0.0, -dAdz)
        w = np.where(axis & ~self._axis_edge[:, None], 0.0, 1.0)   # A/r undefined at a lone axis vertex
        key = (tri * ncls + cls[:, None]).ravel()
        nk = len(self.nodes) * ncls
        ws = np.maximum(np.bincount(key, weights=w.ravel(), minlength=nk), 1e-300)
        br = np.bincount(key, weights=(w * Br).ravel(), minlength=nk) / ws
        bz = np.bincount(key, weights=(w * Bz).ravel(), minlength=nk) / ws
        self._rec = np.stack([br[key].reshape(tri.shape), bz[key].reshape(tri.shape)], -1)
        return self._rec

    def flux(self, r, z):
        """Flux through the circle of radius r at z: 2 pi r A (Wb)."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        z = np.atleast_1d(np.asarray(z, dtype=float))
        e, l1, l2 = self.locate(r, z)
        N, _, _ = p2_basis(l1, l2)
        out = 2 * np.pi * r * np.sum(self.A[self.tri[np.maximum(e, 0)]] * N, axis=1)
        out[e < 0] = np.nan
        return out
