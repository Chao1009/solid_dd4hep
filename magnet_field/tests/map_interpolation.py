#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""How well a GEMC map reproduces the field once SoLID_FieldMapB interpolates it: the bilinear interpolation of the
file against the rzfield solution, at random points binned by their distance to the nearest iron/air interface,
separately in air and in iron; air points within 2 cm of a coil edge (and inside the coil sheets) are reported
apart ("near coil"), since the 1 cm grid cannot resolve the 0.4 cm current sheets either.

  python3 tests/map_interpolation.py --solution work/maps/solutions/solid_v4c_open.npz \\
      --map work/maps/solenoid_SoLID_v4c_2D.dat [--map other_variant.dat ...] [--label v4c] [--json out.json]

Points (fixed seed): --n uniform points in the box around the iron (3 cm margin), plus --n/4 over the whole map
grid. The interface is taken from the mesh of the solution (edges between an air and an iron element), so it is
exact for the solved geometry. The bilinear interpolation follows src/FieldMapB.cpp. A POISSON map can be checked
against the rzfield solution of its own deck (e.g. solenoid_CLEOv10.dat against CLEOv10.am solved with
--mode poisson). Needs numpy and scipy only."""
import argparse
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from rzfield.field import RZSolution  # noqa: E402
from rzfield.gemcmap import read_map  # noqa: E402

BANDS = ((0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, np.inf))


def interface_segments(sol, cls):
    """mesh edges between elements of different class (cls: per-element labels): end points in cm, from the straight
    vertex edges of the elements"""
    tri = sol.tri[:, :3]
    e = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(tri)), 3)
    key = e[:, 0].astype(np.int64) * len(sol.nodes) + e[:, 1]
    order = np.argsort(key, kind="stable")
    same = key[order][1:] == key[order][:-1]
    i1, i2 = order[:-1][same], order[1:][same]
    iface = cls[owner[i1]] != cls[owner[i2]]
    seg = e[i1[iface]]
    return sol.nodes[seg[:, 0]] * 100.0, sol.nodes[seg[:, 1]] * 100.0


def coil_elements(sol):
    """elements of the current-carrying zones"""
    if sol.zone is None or not sol.zones:
        return np.zeros(len(sol.tri), dtype=bool)
    return np.array([float(zd.get("cur", 0.0)) != 0.0 for zd in sol.zones], dtype=bool)[sol.zone]


def distance_to_segments(p, a, b, k=16, chunk=200000):
    """distance (cm) from points p (n, 2) to the nearest of the segments a-b, among the k with the nearest
    midpoints (exact for the 0-3 cm bands used here: interface segments are at most a few cm long)"""
    tree = cKDTree(0.5 * (a + b))
    k = min(k, len(a))
    out = np.empty(len(p))
    for s in range(0, len(p), chunk):
        q = p[s:s + chunk]
        _, idx = tree.query(q, k=k)
        idx = idx.reshape(len(q), k)
        sa, sb = a[idx], b[idx]
        ab = sb - sa
        t = np.clip(np.einsum("nki,nki->nk", q[:, None, :] - sa, ab) / np.maximum(np.einsum("nki,nki->nk", ab, ab),
                                                                                  1e-300), 0.0, 1.0)
        d = np.linalg.norm(q[:, None, :] - (sa + t[..., None] * ab), axis=2)
        out[s:s + chunk] = d.min(axis=1)
    return out


def bilinear(m, r, z):
    """(Br, Bz) of map m at r, z (cm) as SoLID_FieldMapB: bilinear in the cell, 0 outside the grid"""
    ra, za = m["r"], m["z"]
    u = (r - ra[0]) / (ra[1] - ra[0])
    v = (z - za[0]) / (za[1] - za[0])
    inside = (u >= 0) & (u <= len(ra) - 1) & (v >= 0) & (v <= len(za) - 1)
    i = np.minimum(np.clip(u, 0, None).astype(int), len(ra) - 2)
    j = np.minimum(np.clip(v, 0, None).astype(int), len(za) - 2)
    a, b = u - i, v - j
    out = []
    for A in (m["Br"], m["Bz"]):
        val = A[i, j] * (1 - a) * (1 - b) + A[i, j + 1] * (1 - a) * b + A[i + 1, j] * a * (1 - b) \
            + A[i + 1, j + 1] * a * b
        out.append(np.where(inside, val, 0.0))
    return out


def band_stats(dB, B, d, r, z, Bmap):
    out = []
    for lo, hi in BANDS:
        k = (d >= lo) & (d < hi)
        if not k.any():
            out.append(dict(band_cm=[lo, hi], n=0))
            continue
        x = dB[k]
        rel = x / np.maximum(B[k], 1.0)
        w = int(np.argmax(x))
        out.append(dict(band_cm=[lo, None if np.isinf(hi) else hi], n=int(k.sum()), median_G=float(np.median(x)),
                        p90_G=float(np.percentile(x, 90)), p99_G=float(np.percentile(x, 99)), max_G=float(x.max()),
                        median_rel=float(np.median(rel)), p99_rel=float(np.percentile(rel, 99)),
                        worst=dict(r_cm=float(r[k][w]), z_cm=float(z[k][w]), map_G=float(Bmap[k][w]),
                                   solution_G=float(B[k][w]))))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solution", required=True, help="rzfield solution (.npz) the maps were made from")
    ap.add_argument("--map", action="append", required=True, help="GEMC map(s) to check (repeatable)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--n", type=int, default=3000000, help="points near the iron (default 3e6)")
    ap.add_argument("--seed", type=int, default=20261106)
    ap.add_argument("--json", default=None, help="write the numbers here")
    a = ap.parse_args(argv)

    sol = RZSolution.load(a.solution)
    sa, sb = interface_segments(sol, sol.air_elements())
    ca, cb = interface_segments(sol, coil_elements(sol))
    maps = [read_map(p) for p in a.map]
    m0 = maps[0]
    lo = np.minimum(sa.min(0), sb.min(0)) - 3.0
    hi = np.maximum(sa.max(0), sb.max(0)) + 3.0
    rng = np.random.default_rng(a.seed)
    rmin, rmax, zmin, zmax = m0["r"][0], m0["r"][-1], m0["z"][0], m0["z"][-1]
    n2 = a.n // 4
    r = np.concatenate([rng.uniform(max(lo[0], rmin), min(hi[0], rmax), a.n), rng.uniform(rmin, rmax, n2)])
    z = np.concatenate([rng.uniform(max(lo[1], zmin), min(hi[1], zmax), a.n), rng.uniform(zmin, zmax, n2)])
    Br, Bz, in_air = sol.B_cm_gauss(r, z, return_air=True)
    ok = np.isfinite(Br) & np.isfinite(Bz)
    r, z, Br, Bz, in_air = r[ok], z[ok], Br[ok], Bz[ok], in_air[ok]
    B = np.hypot(Br, Bz)
    P = np.column_stack([r, z])
    d = distance_to_segments(P, sa, sb)
    near_coil = distance_to_segments(P, ca, cb) < 2.0 if len(ca) else np.zeros(r.size, dtype=bool)
    air = in_air & ~near_coil
    res = dict(label=a.label or os.path.basename(a.solution), solution=os.path.basename(a.solution),
               n_points=int(r.size), n_interface_segments=int(len(sa)),
               bands_cm=[[x, None if np.isinf(y) else y] for x, y in BANDS], maps={})
    print("== %s: %d points, %d iron/air interface segments" % (res["label"], r.size, len(sa)))
    for path, m in zip(a.map, maps):
        mbr, mbz = bilinear(m, r, z)
        dB = np.hypot(mbr - Br, mbz - Bz)
        Bm = np.hypot(mbr, mbz)
        entry = {}
        for name, sel in (("air", air), ("iron", ~in_air)):
            entry[name] = band_stats(dB[sel], B[sel], d[sel], r[sel], z[sel], Bm[sel])
        far = air & (d >= 2.0)
        entry["air_ge_2cm"] = band_stats(dB[far], B[far], np.zeros(far.sum()), r[far], z[far], Bm[far])[0]
        nc = in_air & near_coil
        entry["near_coil"] = band_stats(dB[nc], B[nc], np.zeros(nc.sum()), r[nc], z[nc], Bm[nc])[0]
        for k in ("air_ge_2cm", "near_coil"):
            entry[k]["band_cm"] = None
        res["maps"][os.path.basename(path)] = entry
        print("  map %s" % os.path.basename(path))
        rows = [("air", s_) for s_ in entry["air"]] + [("iron", s_) for s_ in entry["iron"]] + \
            [("air >= 2 cm from iron and coils", entry["air_ge_2cm"]), ("air < 2 cm from a coil", entry["near_coil"])]
        for name, s_ in rows:
            if not s_["n"]:
                continue
            band = "" if s_["band_cm"] is None else "%.1f-%s cm" % (
                s_["band_cm"][0], "%.1f" % s_["band_cm"][1] if s_["band_cm"][1] else "inf")
            print("    %-4s %-10s n %8d  |dB| median %8.2f  p90 %8.1f  p99 %8.1f  max %8.1f G   rel median %.3f"
                  "  worst (%.2f, %.2f): map %.0f, solution %.0f G" % (
                      name, band, s_["n"], s_["median_G"], s_["p90_G"], s_["p99_G"], s_["max_G"], s_["median_rel"],
                      s_["worst"]["r_cm"], s_["worst"]["z_cm"], s_["worst"]["map_G"], s_["worst"]["solution_G"]))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
