#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 SoLID collaboration
"""2D axisymmetric (r, z) model of the SoLID CLEO-II magnet, extracted from the DD4hep geometry.

The DD4hep magnet (DetElements MagnetYoke and MagnetCoil, plugin SoLID_CLEOMagnet) is loaded from a compact file
(a magnet test top-level such as test/test_magnet_v4.xml, or a full configuration such as solid_SIDIS_He3.xml).
Every volume is turned into an r-z region from its TGeo shape, the iron pieces are merged, the coil current is
distributed from a parameter set, and the result is written as a POISSON Superfish deck (.am, CLEOv10.am
conventions: region 1 = problem box, mat=1 + cur= coil regions, mat=2 iron) plus a JSON sidecar.

Sub-commands
  extract COMPACT -o DECK.am   needs DD4hep python (eic_xl container after thissolid.sh); writes DECK.am, DECK.json
  compare MODEL REFERENCE      numpy + matplotlib only: union and region-by-region comparison with a reference deck
                               (.am) or the Stage 1/2 html (its JSON data block), text report and overlay plot
  plot-fraction MODEL NPZ      plot the sampled azimuthal iron fraction with the model outline
  check MODEL.json --repo DIR  no DD4hep needed: are the repository files a deck was extracted from unchanged (in
                               what the geometry reads)? exit 0 yes, 1 only warnings, 2 a magnet file changed

Conversion rules (details and rationale in models/README.md)
  * TGeoPcon / TGeoTube(Seg) / TGeoCone(Seg) with full phi on the axis: exact r-z polygon.
  * TGeoPgon on the axis: polygon at the apothem ("flat") radii (default, as the POISSON decks) or, with
    --convention equal-area, at the radii of the circle with the same area (x sqrt(n tan(pi/n) / pi)).
  * TGeoCompositeShape whose leaves are such coaxial shapes: the boolean is done in the r-z plane (each leaf in
    its own convention), then checked against geometric sampling.
  * Anything else (boxes placed off the axis such as the spacer bars, partial-phi solids): sampled with
    TGeoShape::Contains on an (r, phi, z) grid; radial and axial extents by ray casting (DistFromOutside /
    DistFromInside) through the centroid of each copy; replaced by an equivalent full ring of the same iron
    volume whose axial length is scaled by the azimuthal coverage (the POISSON decks' spacer convention).
  * Iron = materials with an Fe mass fraction >= --iron-min-fe (default 0.9; stainless steel is not iron).
    Coils = volumes matched by the parameter set (models/cleo2_coils.json); every other volume is non-magnetic.
  * Gaps and offsets below --snap (default 0.5 cm) between neighbouring iron pieces are closed by snapping
    coordinates to the earlier piece (reported); touching iron pieces are merged into simple polygons without
    holes (largest contact first).
"""

from __future__ import annotations

import argparse
import bisect
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PARAMS = os.path.join(HERE, "models", "cleo2_coils.json")
STEEL_DENSITY = 7.87  # g/cm3, POISSON decks' "1010 steel"

# =================================================================================================================
# r-z geometry engine: a region is a union of trapezoids (z0, z1, lo0, lo1, hi0, hi1) in the (z, r) half plane,
# z0 < z1, r = lo(z) .. hi(z) linear in z. Booleans are done slab by slab (exact for straight edges).
# =================================================================================================================
EPS = 1e-7   # cm, coordinate noise
VTOL = 1e-6  # cm, vertex identity when linking boundary edges


def _cluster(values, tol=EPS):
    out = []
    for v in sorted(values):
        if not out or v - out[-1] > tol:
            out.append(v)
    return out


def _nearest(sorted_vals, v):
    i = bisect.bisect_left(sorted_vals, v)
    cands = [sorted_vals[j] for j in (i - 1, i) if 0 <= j < len(sorted_vals)]
    return min(cands, key=lambda c: abs(c - v))


def _ev(line, t):
    if t == 0.0:
        return line[0]
    if t == 1.0:
        return line[1]
    return line[0] + (line[1] - line[0]) * t


def _at(za, zb, t):
    """z at the fraction t of [za, zb], exact at the ends (adjacent slabs must share the same float)"""
    return za if t == 0.0 else zb if t == 1.0 else za + (zb - za) * t


def _slab(traps, za, zb):
    out = []
    for z0, z1, l0, l1, h0, h1 in traps:
        if z0 <= za and z1 >= zb:
            ta, tb = (za - z0) / (z1 - z0), (zb - z0) / (z1 - z0)
            out.append(((l0 + (l1 - l0) * ta, l0 + (l1 - l0) * tb), (h0 + (h1 - h0) * ta, h0 + (h1 - h0) * tb)))
    return out


_PRED = {
    "union": lambda a, b: a or b,
    "sub": lambda a, b: a and not b,
    "inter": lambda a, b: a and b,
}


def _interval_bool(IA, IB, t, op):
    pred = _PRED[op]
    ev = []
    for k, ints in ((0, IA), (1, IB)):
        for lo, hi in ints:
            a, b = _ev(lo, t), _ev(hi, t)
            if b - a > EPS:
                ev.append((a, k, 1, lo))
                ev.append((b, k, -1, hi))
    ev.sort(key=lambda e: e[0])
    res, cnt, inside, start, i = [], [0, 0], False, None, 0
    while i < len(ev):
        x, j = ev[i][0], i
        while j < len(ev) and ev[j][0] - x <= EPS:
            j += 1
        for e in ev[i:j]:
            cnt[e[1]] += e[2]
        now = pred(cnt[0] > 0, cnt[1] > 0)
        if now and not inside:
            start = (x, ev[i][3])
        elif inside and not now and x - start[0] > EPS:
            res.append((start[1], ev[i][3]))
        inside = now
        i = j
    return res


def _coalesce(traps):
    traps = sorted(traps, key=lambda t: (t[0], t[2]))
    out, ends = [], {}
    for t in traps:
        key = (t[0], round(t[2] / VTOL), round(t[4] / VTOL))
        k = ends.pop(key, None)
        if k is not None:
            p = out[k]
            s = (t[1] - p[0]) / (p[1] - p[0])
            if abs(p[2] + (p[3] - p[2]) * s - t[3]) < VTOL and abs(p[4] + (p[5] - p[4]) * s - t[5]) < VTOL:
                out[k] = (p[0], t[1], p[2], t[3], p[4], t[5])
                ends[(t[1], round(t[3] / VTOL), round(t[5] / VTOL))] = k
                continue
            ends[key] = k
        out.append(t)
        ends[(t[1], round(t[3] / VTOL), round(t[5] / VTOL))] = len(out) - 1
    return out


def _combine(A, B, op):
    zs = _cluster([t[0] for t in A] + [t[1] for t in A] + [t[0] for t in B] + [t[1] for t in B])
    if len(zs) < 2:
        return []

    def snapz(T):
        out = []
        for t in T:
            z0, z1 = _nearest(zs, t[0]), _nearest(zs, t[1])
            if z1 > z0:
                out.append((z0, z1) + tuple(t[2:]))
        return out

    A, B = snapz(A), snapz(B)
    out = []
    for za, zb in zip(zs[:-1], zs[1:]):
        IA, IB = _slab(A, za, zb), _slab(B, za, zb)
        if not IA and (op != "union" or not IB):
            continue
        if op == "inter" and not IB:
            continue
        lines = [ln for iv in IA + IB for ln in iv]
        ts = [0.0, 1.0]
        for i in range(len(lines)):
            for j in range(i + 1, len(lines)):
                da, db = lines[i][0] - lines[j][0], lines[i][1] - lines[j][1]
                if (da > EPS and db < -EPS) or (da < -EPS and db > EPS):
                    ts.append(da / (da - db))
        ts = _cluster(ts, 1e-12)
        ts[-1] = 1.0
        for t0, t1 in zip(ts[:-1], ts[1:]):
            s0, s1 = _at(za, zb, t0), _at(za, zb, t1)
            if s1 - s0 <= EPS:
                continue
            for lo, hi in _interval_bool(IA, IB, 0.5 * (t0 + t1), op):
                out.append((s0, s1, _ev(lo, t0), _ev(lo, t1), _ev(hi, t0), _ev(hi, t1)))
    return _coalesce(out)


def _iunion(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1] + EPS:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _isub(X, Y):
    out, Y = [], _iunion(Y)
    for a, b in _iunion(X):
        cur = a
        for c, d in Y:
            if d <= cur + EPS or c >= b - EPS:
                continue
            if c > cur + EPS:
                out.append((cur, c))
            cur = max(cur, d)
        if b - cur > EPS:
            out.append((cur, b))
    return out


class _VertexPool:
    def __init__(self, tol):
        self.tol, self.grid, self.pts = tol, {}, []

    def id(self, p):
        kz, kr = math.floor(p[0] / self.tol), math.floor(p[1] / self.tol)
        for dz in (-1, 0, 1):
            for dr in (-1, 0, 1):
                for i in self.grid.get((kz + dz, kr + dr), ()):
                    q = self.pts[i]
                    if abs(q[0] - p[0]) <= self.tol and abs(q[1] - p[1]) <= self.tol:
                        return i
        self.pts.append(p)
        self.grid.setdefault((kz, kr), []).append(len(self.pts) - 1)
        return len(self.pts) - 1


def _simplify(loop, tol=VTOL):
    pts = list(loop)
    changed = True
    while changed and len(pts) > 3:
        changed = False
        for i in range(len(pts)):
            p, v, n = pts[i - 1], pts[i], pts[(i + 1) % len(pts)]
            dx, dy = n[0] - p[0], n[1] - p[1]
            L = math.hypot(dx, dy)
            if L < tol:
                continue
            dist = abs(dx * (v[1] - p[1]) - dy * (v[0] - p[0])) / L
            t = ((v[0] - p[0]) * dx + (v[1] - p[1]) * dy) / (L * L)
            if dist < tol and -1e-9 < t < 1 + 1e-9 or math.hypot(v[0] - p[0], v[1] - p[1]) < tol:
                del pts[i]
                changed = True
                break
    return pts


def _signed_area(loop):
    a = 0.0
    for i in range(len(loop)):
        (x0, y0), (x1, y1) = loop[i - 1], loop[i]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def _point_in_loop(p, loop):
    x, y, inside = p[0], p[1], False
    for i in range(len(loop)):
        (x0, y0), (x1, y1) = loop[i - 1], loop[i]
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            inside = not inside
    return inside


def _trace(traps):
    """boundary loops (z, r) of a canonical trapezoid list; outer loops counter-clockwise, holes clockwise"""
    pool, E = _VertexPool(VTOL), []

    def add(p, q):
        a, b = pool.id(p), pool.id(q)
        if a != b:
            E.append((a, b))

    zl = _cluster([t[0] for t in traps] + [t[1] for t in traps], VTOL)
    levels = {}
    for z0, z1, l0, l1, h0, h1 in traps:
        add((z0, l0), (z1, l1))
        add((z1, h1), (z0, h0))
        levels.setdefault(_nearest(zl, z1), ([], []))[0].append((l1, h1))
        levels.setdefault(_nearest(zl, z0), ([], []))[1].append((l0, h0))
    for z, (ends, starts) in levels.items():
        for a, b in _isub(ends, starts):
            add((z, a), (z, b))
        for a, b in _isub(starts, ends):
            add((z, b), (z, a))
    # cancel exactly opposite edges (shared boundaries that survived)
    cnt = defaultdict(int)
    for a, b in E:
        cnt[(a, b)] += 1
    E2 = []
    for a, b in E:
        if cnt[(b, a)] > 0 and cnt[(a, b)] > 0:
            cnt[(b, a)] -= 1
            cnt[(a, b)] -= 1
            continue
        E2.append((a, b))
    E = E2
    P = pool.pts
    out = defaultdict(list)
    for k, (a, b) in enumerate(E):
        out[a].append(k)
    used = [False] * len(E)
    loops = []
    for s in range(len(E)):
        if used[s]:
            continue
        loop, cur = [], s
        while True:
            used[cur] = True
            a, b = E[cur]
            loop.append(P[a])
            cands = [k for k in out[b] if not used[k]]
            if b == E[s][0]:
                cands.append(s)
            if not cands:
                raise RuntimeError("r-z engine: open boundary while tracing a region")
            if len(cands) > 1:
                dx, dy = P[b][0] - P[a][0], P[b][1] - P[a][1]

                def turn(k):
                    c = P[E[k][1]]
                    ex, ey = c[0] - P[b][0], c[1] - P[b][1]
                    return math.atan2(dx * ey - dy * ex, dx * ex + dy * ey)

                nxt = max(cands, key=turn)
            else:
                nxt = cands[0]
            if nxt == s:
                break
            cur = nxt
        loop = _simplify(loop)
        if len(loop) >= 3 and abs(_signed_area(loop)) > 1e-9:
            loops.append(loop)
    return loops


def _traps_from_loops(loops_zr):
    edges = []
    for L in loops_zr:
        for i in range(len(L)):
            p, q = L[i], L[(i + 1) % len(L)]
            if abs(p[0] - q[0]) > EPS:
                edges.append((p, q))
    zs = _cluster([p[0] for L in loops_zr for p in L])
    traps = []
    for za, zb in zip(zs[:-1], zs[1:]):
        lines = []
        for (z1, r1), (z2, r2) in edges:
            if min(z1, z2) <= za + EPS and max(z1, z2) >= zb - EPS:
                lines.append((r1 + (r2 - r1) * (za - z1) / (z2 - z1), r1 + (r2 - r1) * (zb - z1) / (z2 - z1)))
        ts = [0.0, 1.0]
        for i in range(len(lines)):
            for j in range(i + 1, len(lines)):
                da, db = lines[i][0] - lines[j][0], lines[i][1] - lines[j][1]
                if (da > EPS and db < -EPS) or (da < -EPS and db > EPS):
                    ts.append(da / (da - db))
        ts = _cluster(ts, 1e-12)
        ts[-1] = 1.0
        for t0, t1 in zip(ts[:-1], ts[1:]):
            tm = 0.5 * (t0 + t1)
            xs = sorted(lines, key=lambda ln: _ev(ln, tm))
            if len(xs) % 2:
                raise ValueError("r-z engine: polygon is not closed")
            s0, s1 = _at(za, zb, t0), _at(za, zb, t1)
            for k in range(0, len(xs), 2):
                traps.append((s0, s1, _ev(xs[k], t0), _ev(xs[k], t1), _ev(xs[k + 1], t0), _ev(xs[k + 1], t1)))
    return traps


class Region:
    """Union of trapezoids in the (z, r) half plane (canonical: disjoint per z slab)."""

    __slots__ = ("traps",)

    def __init__(self, traps=(), canonical=False):
        t = [tuple(float(x) for x in tr) for tr in traps]
        self.traps = t if canonical else _combine(t, [], "union")

    @classmethod
    def from_planes(cls, z, rmin, rmax):
        tr = [(z[i], z[i + 1], rmin[i], rmin[i + 1], rmax[i], rmax[i + 1])
              for i in range(len(z) - 1) if z[i + 1] - z[i] > EPS]
        return cls(tr)

    @classmethod
    def from_loops_rz(cls, loops_rz):
        return cls(_traps_from_loops([[(z, r) for r, z in L] for L in loops_rz]))

    @classmethod
    def rect(cls, r0, r1, z0, z1):
        return cls([(z0, z1, r0, r0, r1, r1)])

    def __or__(self, o):
        return Region(_combine(self.traps, o.traps, "union"), True)

    def __sub__(self, o):
        return Region(_combine(self.traps, o.traps, "sub"), True)

    def __and__(self, o):
        return Region(_combine(self.traps, o.traps, "inter"), True)

    def empty(self):
        return not self.traps

    def area(self):
        return sum(0.5 * (z1 - z0) * ((h0 - l0) + (h1 - l1)) for z0, z1, l0, l1, h0, h1 in self.traps)

    def volume(self):
        """volume of the body of revolution (cm3)"""
        v = 0.0
        for z0, z1, l0, l1, h0, h1 in self.traps:
            dh, dl = h1 - h0, l1 - l0
            v += (z1 - z0) * ((h0 * h0 + h0 * dh + dh * dh / 3.0) - (l0 * l0 + l0 * dl + dl * dl / 3.0))
        return math.pi * v

    def bbox(self):
        if not self.traps:
            return None
        rs = [x for t in self.traps for x in t[2:]]
        return (min(rs), max(rs), min(t[0] for t in self.traps), max(t[1] for t in self.traps))

    def loops(self):
        """boundary loops as lists of (r, z); outer loops counter-clockwise in (z, r), holes clockwise"""
        return [[(r, z) for z, r in L] for L in _trace(self.traps)]

    def polygons(self):
        """[(outer_rz, [holes_rz])] with each hole assigned to the smallest outer loop around it"""
        lz = _trace(self.traps)
        outers = [L for L in lz if _signed_area(L) > 0]
        holes = [L for L in lz if _signed_area(L) < 0]
        res = [[L, []] for L in sorted(outers, key=_signed_area)]
        for h in holes:
            for item in res:
                if _point_in_loop(h[0], item[0]) or _point_in_loop(
                        ((h[0][0] + h[1][0]) / 2, (h[0][1] + h[1][1]) / 2), item[0]):
                    item[1].append(h)
                    break
        return [([(r, z) for z, r in o], [[(r, z) for z, r in h] for h in hs]) for o, hs in res]

    def perimeter(self):
        p = 0.0
        for L in _trace(self.traps):
            for i in range(len(L)):
                p += math.hypot(L[i][0] - L[i - 1][0], L[i][1] - L[i - 1][1])
        return p


def union_all(regions):
    out = Region()
    for r in regions:
        out = out | r
    return out


# =================================================================================================================
# POISSON deck (.am) reading and writing
# =================================================================================================================
def parse_am(text):
    """minimal AUTOMESH/POISSON input reader: returns dict(title, control, regions=[dict(name, mat, cur, pts)])
    with pts as (r, z) = (x, y). Region 1 is the problem box. Arcs (nt=2) are not supported."""
    lines = text.splitlines()
    title = lines[0].strip() if lines else ""
    blocks, buf, last_comment, name_for = [], None, None, None
    for line in lines[1:]:
        code, _, comment = line.partition(";")
        if buf is None and "&" not in code:
            if comment.strip() and not code.strip():
                last_comment = comment.strip()
            continue
        s = code
        while s:
            if buf is None:
                i = s.find("&")
                if i < 0:
                    break
                buf, s, name_for = "", s[i + 1:], last_comment
            j = s.find("&")
            if j < 0:
                buf += s + " "
                s = ""
            else:
                buf += s[:j]
                blocks.append((buf, name_for))
                buf, s = None, s[j + 1:]
    num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"
    regions, control = [], {}
    for body, cname in blocks:
        m = re.match(r"\s*(\w+)(.*)", body, re.S)
        if not m:
            continue
        kind, rest = m.group(1).lower(), m.group(2)
        kv = {k.lower(): float(v.replace("d", "e").replace("D", "e"))
              for k, v in re.findall(r"(\w+)\s*=\s*(" + num + ")", rest)}
        if kind == "reg":
            if "kprob" in kv or not regions and not control:
                control = kv
                regions.append(dict(name="problem box", mat=0, cur=0.0, pts=[]))
            else:
                nm = re.sub(r"\s*Region\s+\d+.*$", "", cname or "").strip() or "region %d" % (len(regions) + 1)
                regions.append(dict(name=nm, mat=int(kv.get("mat", 1)), cur=kv.get("cur", 0.0), pts=[]))
        elif kind == "po":
            if "nt" in kv and int(kv["nt"]) != 1:
                raise ValueError("parse_am: arcs (nt=%d) are not supported" % int(kv["nt"]))
            if not regions:
                raise ValueError("parse_am: point before the first &reg")
            regions[-1]["pts"].append((kv["x"], kv["y"]))
    return dict(title=title, control=control, regions=regions)


def _fmt(v):
    s = "%.3f" % v
    return "0.000" if s == "-0.000" else s


BASE_ZONES = [(80.0, 2.0), (120.0, 1.0), (150.0, 0.5), (156.0, 0.1), (190.0, 0.5), (320.0, 1.0), (400.0, 2.0)]


def radial_zones(box_r, refine=()):
    """radial mesh zones [(x_end, spacing)]: those of CLEOv10.am (4 cm beyond 400 cm), plus refined zones
    (a, b, h) around thin features, aligned to the base grid of the zone they fall in"""
    zones = [z for z in BASE_ZONES if z[0] < box_r] + [(box_r, 4.0)]

    def base_h(x):
        for xe, h in zones:
            if x < xe:
                return h
        return zones[-1][1]

    for a, b, h in sorted(refine):
        ha = base_h(a)
        a = math.floor(a / ha) * ha
        hb = base_h(b)
        b = math.ceil(b / hb) * hb
        new, x0 = [], 0.0
        for xe, hz in zones:
            if xe <= a or x0 >= b:
                new.append((xe, hz))
            else:
                if x0 < a:
                    new.append((a, hz))
                if hz > h:
                    new.append((min(b, xe), h))
                else:
                    new.append((min(b, xe), hz))
                if xe > b:
                    new.append((xe, hz))
            x0 = xe
        zones = []
        for xe, hz in new:  # merge neighbours with the same spacing
            if zones and abs(zones[-1][1] - hz) < 1e-12:
                zones[-1] = (xe, hz)
            else:
                zones.append((xe, hz))
    return zones


def mesh_block(zones):
    """xreg/kreg lines for radial zones; uniform 1 cm in z (dy = dx)"""
    out, k, x0 = [], 1, 0.0
    for i, (x, h) in enumerate(zones):
        k += int(round((x - x0) / h))
        if i < len(zones) - 1:
            out.append("xreg%d=%7.1f,  kreg%d=%5d," % (i + 1, x, i + 1, k))
        x0 = x
    out.append("              kmax =%5d," % k)
    return out


def thin_radial_features(regions, zones, cells=2.0):
    """(a, b, h) for parts of regions thinner in r than `cells` radial mesh cells"""
    out = []
    for reg in regions:
        for z0, z1, l0, l1, h0, h1 in reg.traps:
            lo, hi = min(l0, l1), max(h0, h1)
            w = min(h0 - l0, h1 - l1)
            if w <= 0:
                continue
            hloc = [hz for (xe, hz), xs in zip(zones, [0.0] + [z[0] for z in zones[:-1]]) if xs < hi and xe > lo]
            if hloc and w < cells * min(hloc):
                h = max(0.05, round(w / 3.0 / 0.05) * 0.05)
                out.append((lo - 1.0, hi + 1.0, h))
    return out


def _comment(text, width=100):
    """comment line(s) '; text', wrapped at `width` characters (continuations indented)"""
    import textwrap
    if not text.strip():
        return [";"]
    lead = len(text) - len(text.lstrip(" "))
    return ["; " + s for s in textwrap.wrap(text, width - 2, subsequent_indent=" " * (lead + 2),
                                            break_long_words=False, break_on_hyphens=False)]


def write_am(path, model):
    L = []
    L.append(model["title"][:79])
    for c in model["header"]:
        L.extend(_comment(c))
    L.append(";")
    L.append("&reg kprob=0,           ; Declares a POISSON problem")
    L.append("")
    L.append("icylin=1,               ; Cylindrical symmetry: x = r, y = z (cm)")
    L.append("ienergy=1,              ; Calculate the full energy")
    L.append("xjfact=1.0,             ; The current multiplicative factor")
    L.append("dx=1.0,                 ; Mesh interval (dy = dx = 1 cm; radial zones below)")
    L.append("rhogam=0.02,            ; relaxation parameter (0.02 also converges in saturated steel)")
    L.append("mode=0,                 ; Variable permeability for material 2 (internal iron table, no MT)")
    L.append("")
    L.extend(mesh_block(model["mesh_zones"]))
    L.append("")
    L.append("nbsup=0,                ; Dirichlet boundary condition on upper edge")
    L.append("nbslf=0,                ; Dirichlet boundary condition on left edge")
    L.append("nbsrt=0,                ; Dirichlet boundary condition on right edge")
    L.append("nbslo=0,                ; Dirichlet boundary condition on bottom")
    L.append("ibound=0 &")
    L.append("")
    (r0, r1), (z0, z1) = model["box"]["r"], model["box"]["z"]
    L.append(";Problem box (hall)  Region 1")
    for r, z in ((r0, z0), (r0, z1), (r1, z1), (r1, z0), (r0, z0)):
        L.append("&po x=%9s, y=%9s &" % (_fmt(r), _fmt(z)))
    for k, reg in enumerate(model["regions"]):
        L.append("")
        for c in reg.get("comments", []):
            L.extend(_comment(c))
        L.append("; %s  Region %d" % (reg["name"], k + 2))
        if reg["mat"] == 1:
            L.append("&reg    mat=1,   cur=%.1f &" % reg["cur"])
        else:
            L.append("&reg    mat=%d  &" % reg["mat"])
        pts = list(reg["pts"]) + [reg["pts"][0]]
        for r, z in pts:
            L.append("&po x=%9s, y=%9s &" % (_fmt(r), _fmt(z)))
    for blk in model.get("commented_regions", []):
        L.append("")
        for c in blk:
            L.extend(_comment(c) if not c.startswith("&") else ["; " + c])
    L.append("")
    with open(path, "w") as f:
        f.write("\n".join(L))


def _oriented_rz(loop_rz):
    """deck vertex order: start at the vertex with the smallest (z, r), counter-clockwise in (z, r)"""
    pts = list(loop_rz)
    if _signed_area([(z, r) for r, z in pts]) < 0:
        pts.reverse()
    k = min(range(len(pts)), key=lambda i: (round(pts[i][1], 6), round(pts[i][0], 6)))
    return pts[k:] + pts[:k]


# =================================================================================================================
# provenance helpers
# =================================================================================================================
def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def compact_files(top):
    """top-level compact + every xml it includes (recursively), with ${VAR} expanded; the paths as written (made
    absolute), so that ${DD4hepINSTALL}/... stays recognisable when it is a view of symbolic links"""
    seen, order = set(), []

    def visit(p):
        p = os.path.abspath(p)
        rp = os.path.realpath(p)
        if rp in seen or not os.path.isfile(rp):
            return
        seen.add(rp)
        order.append(p)
        try:
            txt = open(rp, errors="replace").read()
        except OSError:
            return
        txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
        for ref in re.findall(r"<(?:include|file|gdmlFile)\b[^>]*\bref\s*=\s*\"([^\"]+)\"", txt):
            q = os.path.expandvars(ref)
            if not os.path.isabs(q):
                q = os.path.join(os.path.dirname(p), q)
            if q.endswith(".xml"):
                visit(q)

    visit(top)
    return order


def without_fields(compact):
    """Copy of a top-level compact without its <fields> blocks and compact/fields includes, in a new temporary
    directory (relative include paths are made absolute). Loading it needs no field map."""
    txt = re.sub(r"<!--.*?-->", "", open(compact).read(), flags=re.S)
    txt = _FIELDS_INCLUDE.sub("", _FIELDS_BLOCK.sub("", txt))

    def absref(m):
        ref = m.group(2)
        if not (ref.startswith("$") or os.path.isabs(ref) or "://" in ref):
            ref = os.path.join(os.path.dirname(compact), ref)
        return m.group(1) + ref + m.group(3)

    txt = re.sub(r"(\bref\s*=\s*\")([^\"]+)(\")", absref, txt)
    out = os.path.join(tempfile.mkdtemp(prefix="extract_rz_"), os.path.basename(compact))
    with open(out, "w") as f:
        f.write(txt)
    return out


_FIELDS_BLOCK = re.compile(r"<fields\b.*?</fields>", re.S)
_FIELDS_INCLUDE = re.compile(r"<include\b[^>]*\bref\s*=\s*\"[^\"]*/fields/[^\"]*\"[^>]*/>")
MAGNET_FILES = ("compact/magnet/", "test/test_magnet")        # a change here makes a committed deck stale
PLUGIN_SOURCES = ("src/CLEOMagnet_geo.cpp",)                  # recorded too (a change gives a warning)


def geometry_digest(path, no_field=False):
    """sha256 of what the geometry reads from an xml file: comments and <info> removed, whitespace collapsed, and,
    for a top-level loaded with --no-field, the <fields> blocks and compact/fields includes removed as well (the same
    rules as without_fields). Edits of comments, titles or the field map of a test top-level leave it unchanged."""
    txt = open(path, errors="replace").read()
    if not path.endswith(".xml"):
        return hashlib.sha256(txt.encode()).hexdigest()
    txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    txt = re.sub(r"<info\b[^>]*/>|<info\b.*?</info>", "", txt, flags=re.S)
    if no_field:
        txt = _FIELDS_INCLUDE.sub("", _FIELDS_BLOCK.sub("", txt))
    return hashlib.sha256(" ".join(txt.split()).encode()).hexdigest()


def run_check(args):
    """compare the repository files recorded in deck sidecars with a checkout (see make_maps.sh -n)"""
    stale = warn = 0
    for side in args.models:
        name = os.path.splitext(os.path.basename(side))[0]
        try:
            d = json.load(open(side))
            files = d["files"]
        except (OSError, ValueError, KeyError) as e:
            print("  %-18s cannot read %s (%s): not checked" % (name, side, e))
            warn = 1
            continue
        no_field = bool(d.get("options", {}).get("no_field"))
        changed, cosmetic, missing, plugin, n = [], [], [], [], 0
        for f in files:
            p = f["path"]
            if p.startswith(("${", ".../", "/")):
                continue                          # outside the repository (the container's DD4hep): not checked
            q = os.path.join(args.repo, p)
            if not os.path.isfile(q):
                missing.append(p)
                continue
            n += 1
            if _sha(q) == f["sha256"]:
                continue
            if f.get("geometry_sha256") and geometry_digest(q, no_field and p == d.get("compact")) == \
                    f["geometry_sha256"]:
                cosmetic.append(p)
            elif f.get("role") == "plugin source":
                plugin.append(p)
            else:
                changed.append(p)
        magnet = [p for p in changed if p.startswith(MAGNET_FILES)]
        other = [p for p in changed if p not in magnet]
        if magnet:
            stale = 1
            print("  %-18s STALE: changed since the deck was extracted: %s" % (name, " ".join(magnet)))
        if other or missing or plugin:
            warn = 1
            print("  %-18s warning: %s" % (name, "; ".join(x for x in (
                ("changed: " + " ".join(other) + " (the deck depends on them only if they hold magnet parameters)")
                if other else "", ("magnet plugin source changed: " + " ".join(plugin) + " (re-extract to be sure)")
                if plugin else "", ("not found: " + " ".join(missing)) if missing else "") if x)))
        if n == 0 and not missing:
            warn = 1
            print("  %-18s warning: no repository-relative file recorded in %s: not checked" % (name, side))
        elif not (magnet or other or missing or plugin):
            print("  %-18s ok: the %d repository files it was extracted from are unchanged%s" % (
                name, n, (" in what the geometry reads (only comments, titles or field includes changed: %s)"
                          % " ".join(cosmetic)) if cosmetic else ""))
    return 2 if stale else (1 if warn else 0)


def git_info(path):
    """root, commit, date and branch of the git checkout holding path (read only: no optional locks)"""
    d = os.path.dirname(os.path.realpath(path))
    git = ["git", "--no-optional-locks", "-C"]
    try:
        root = subprocess.run(git + [d, "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                              timeout=20).stdout.strip()
        if not root:
            return None
        head = subprocess.run(git + [root, "log", "-1", "--format=%H %ci"], capture_output=True, text=True,
                              timeout=20).stdout.strip()
        branch = subprocess.run(git + [root, "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
                                text=True, timeout=20).stdout.strip()
        return dict(root=root, commit=head.split()[0] if head else "", date=" ".join(head.split()[1:]),
                    branch=branch)
    except (OSError, subprocess.SubprocessError):
        return None


def git_dirty(root, files):
    root = os.path.realpath(root)
    rel = [os.path.relpath(os.path.realpath(f), root) for f in files
           if os.path.realpath(f).startswith(root + os.sep)]
    if not rel:
        return []
    try:
        out = subprocess.run(["git", "--no-optional-locks", "-C", root, "status", "--porcelain", "--"] + rel,
                             capture_output=True, text=True, timeout=20).stdout
        return [ln[3:] for ln in out.splitlines() if ln.strip()]
    except (OSError, subprocess.SubprocessError):
        return []


def portable_path(path, root=None):
    """path as recorded in the decks and sidecars, without host-specific directories: relative to the repository
    root (root, else the directory holding magnet_field/) when inside it; else with the prefix of an environment
    variable (${DD4hepINSTALL}/..., the container's DD4hep); else '.../<file name>'. make_maps.sh -n checks the
    sha256 of the repository-relative files against the checkout."""
    paths = (os.path.abspath(path), os.path.realpath(path))
    bases = [(root or os.path.dirname(HERE), None)]
    bases += [(os.environ[v], "${%s}" % v) for v in ("DD4hepINSTALL", "DETECTOR_PATH", "ROOTSYS") if os.environ.get(v)]
    for base, label in bases:
        for b in dict.fromkeys((os.path.abspath(base), os.path.realpath(base))):
            for p in paths:
                if p == b or p.startswith(b.rstrip(os.sep) + os.sep):
                    rel = os.path.relpath(p, b)
                    return rel if label is None else label + "/" + rel
    return ".../" + os.path.basename(path)


# =================================================================================================================
# DD4hep / TGeo extraction (imports ROOT and dd4hep lazily)
# =================================================================================================================
def _mat(m):
    rot, tr = m.GetRotationMatrix(), m.GetTranslation()
    return (np.array([rot[i] for i in range(9)], dtype=float).reshape(3, 3),
            np.array([tr[i] for i in range(3)], dtype=float))


def _coaxial(R, T):
    return (abs(abs(R[2, 2]) - 1) < 1e-9 and max(abs(R[0, 2]), abs(R[1, 2]), abs(R[2, 0]), abs(R[2, 1])) < 1e-9
            and abs(T[0]) < 1e-9 and abs(T[1]) < 1e-9)


def k_equal_area(n):
    return math.sqrt(n * math.tan(math.pi / n) / math.pi)


def _material_info(mat):
    """(name, density g/cm3, Fe mass fraction, Z list)"""
    name, dens = mat.GetName(), mat.GetDensity()
    fe, zs = 0.0, []
    if mat.IsMixture():
        mix = mat
        n = mix.GetNelements()
        w = mix.GetWmixt()
        for i in range(n):
            z = int(round(mix.GetElement(i).Z()))
            zs.append(z)
            if z == 26:
                fe += float(w[i])
    else:
        z = int(round(mat.GetZ()))
        zs.append(z)
        fe = 1.0 if z == 26 else 0.0
    return name, float(dens), fe, zs


def _shape_planes(shape):
    """(z, rmin, rmax, nsides or 0) in the shape frame for full-phi solids of revolution / prisms, else None"""
    cls = shape.ClassName()
    full = lambda dphi: dphi >= 360.0 - 1e-9
    if cls == "TGeoPcon" and full(shape.GetDphi()):
        n = shape.GetNz()
        return ([shape.GetZ(i) for i in range(n)], [shape.GetRmin(i) for i in range(n)],
                [shape.GetRmax(i) for i in range(n)], 0)
    if cls == "TGeoPgon" and full(shape.GetDphi()):
        n = shape.GetNz()
        return ([shape.GetZ(i) for i in range(n)], [shape.GetRmin(i) for i in range(n)],
                [shape.GetRmax(i) for i in range(n)], int(shape.GetNedges()))
    if cls == "TGeoTube" or (cls == "TGeoTubeSeg" and full((shape.GetPhi2() - shape.GetPhi1()) % 360.0 or 360.0)):
        dz = shape.GetDz()
        return [-dz, dz], [shape.GetRmin()] * 2, [shape.GetRmax()] * 2, 0
    if cls == "TGeoCone" or (cls == "TGeoConeSeg" and full((shape.GetPhi2() - shape.GetPhi1()) % 360.0 or 360.0)):
        dz = shape.GetDz()
        return [-dz, dz], [shape.GetRmin1(), shape.GetRmin2()], [shape.GetRmax1(), shape.GetRmax2()], 0
    return None


def shape_rz(shape, R, T, conv):
    """exact r-z regions of a coaxial shape: (model region, equal-area region, description) or None"""
    if not _coaxial(R, T):
        return None
    pl = _shape_planes(shape)
    if pl is not None:
        z, rmin, rmax, n = pl
        s = 1.0 if R[2, 2] > 0 else -1.0
        zg = [s * zz + T[2] for zz in z]
        if s < 0:
            zg, rmin, rmax = zg[::-1], rmin[::-1], rmax[::-1]
        k_model = k_equal_area(n) if (n and conv == "equal-area") else 1.0
        k_eq = k_equal_area(n) if n else 1.0
        desc = ("%s(%d sides, apothems)" % (shape.ClassName(), n)) if n else shape.ClassName()
        return (Region.from_planes(zg, [k_model * r for r in rmin], [k_model * r for r in rmax]),
                Region.from_planes(zg, [k_eq * r for r in rmin], [k_eq * r for r in rmax]), desc, bool(n))
    if shape.ClassName() == "TGeoCompositeShape":
        bn = shape.GetBoolNode()
        op = {0: "union", 1: "inter", 2: "sub"}.get(int(bn.GetBooleanOperator()))
        if op is None:
            return None
        RL, TL = _mat(bn.GetLeftMatrix())
        RR, TR = _mat(bn.GetRightMatrix())
        left = shape_rz(bn.GetLeftShape(), R @ RL, R @ TL + T, conv)
        right = shape_rz(bn.GetRightShape(), R @ RR, R @ TR + T, conv)
        if left is None or right is None:
            return None
        f = {"union": Region.__or__, "inter": Region.__and__, "sub": Region.__sub__}[op]
        sym = {"union": "+", "inter": "*", "sub": "-"}[op]
        return (f(left[0], right[0]), f(left[1], right[1]), "(%s %s %s)" % (left[2], sym, right[2]),
                left[3] or right[3])
    return None


def collect_leaves(desc, det_names):
    leaves = []
    for dname in det_names:
        de = desc.detector(dname)
        if not de.isValid():
            raise SystemExit("extract: DetElement %s not found in the compact" % dname)
        R0, T0 = _mat(de.nominal().worldTransformation())
        top = de.placement().ptr().GetVolume()

        def walk(vol, R, T):
            for i in range(vol.GetNdaughters()):
                node = vol.GetNode(i)
                Rl, Tl = _mat(node.GetMatrix())
                Rg, Tg = R @ Rl, R @ Tl + T
                dv = node.GetVolume()
                if dv.IsAssembly():
                    walk(dv, Rg, Tg)
                    continue
                if dv.GetNdaughters():
                    raise SystemExit("extract: volume %s has daughters (not supported: model them as separate "
                                     "volumes)" % dv.GetName())
                leaves.append(dict(det=dname, node=node.GetName(), volume=dv.GetName(), shape=dv.GetShape(),
                                   R=Rg, T=Tg, material=_material_info(dv.GetMaterial())))

        if top.IsAssembly() or top.GetNdaughters():
            walk(top, R0, T0)
        else:
            leaves.append(dict(det=dname, node=de.placement().ptr().GetName(), volume=top.GetName(),
                               shape=top.GetShape(), R=R0, T=T0, material=_material_info(top.GetMaterial())))
    return leaves


def _short(det, vol):
    return vol[len(det) + 1:] if vol.startswith(det + "_") else vol


def _bbox_corners(shape):
    dx, dy, dz = shape.GetDX(), shape.GetDY(), shape.GetDZ()
    o = shape.GetOrigin()
    o = np.array([o[0], o[1], o[2]], dtype=float)
    return np.array([o + np.array([sx * dx, sy * dy, sz * dz]) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])


def sample_fraction(leaves, r_edges, z_edges, nphi, max_points=2_000_000):
    """volume fraction f[i, j] of the annular cell [r_i, r_i+1] x [z_j, z_j+1] inside the union of the leaves,
    from TGeoShape::Contains_v at nphi points per cell: phi on a regular grid (golden-ratio phase per radius
    against aliasing with n-fold structures), r and z inside the cell from an R2/R3 low-discrepancy sequence
    (unbiased at boundaries that do not fall on the grid); points weighted by r"""
    dr, dz = np.diff(r_edges), np.diff(z_edges)
    nr, nz = len(dr), len(dz)
    k = np.arange(nphi)
    u = (0.5 + 0.7548776662466927 * k) % 1.0
    w = (0.5 + 0.5698402909980532 * k) % 1.0
    phase = (np.arange(nr) * 0.6180339887498949) % 1.0
    phi = (k[None, :] + phase[:, None]) * (2 * math.pi / nphi)
    rr = r_edges[:-1, None] + u[None, :] * dr[:, None]
    X = (rr * np.cos(phi)).ravel()
    Y = (rr * np.sin(phi)).ravel()
    rsum = rr.sum(axis=1)
    per_row = nr * nphi
    step = max(1, max_points // per_row)
    f = np.zeros((nr, nz))
    for j0 in range(0, nz, step):
        nj = min(step, nz - j0)
        n = per_row * nj
        Z = z_edges[j0:j0 + nj, None] + w[None, :] * dz[j0:j0 + nj, None]
        pts = np.empty((n, 3))
        pts[:, 0] = np.tile(X, nj)
        pts[:, 1] = np.tile(Y, nj)
        pts[:, 2] = np.broadcast_to(Z[:, None, :], (nj, nr, nphi)).ravel()
        inside = np.zeros(n, dtype=bool)
        for lf in leaves:
            R, T = lf["R"], lf["T"]
            if np.allclose(R, np.eye(3)) and not T.any():
                loc = pts
            else:
                loc = (pts - T) @ R
            loc = np.ascontiguousarray(loc.ravel())
            ins = np.zeros(n, dtype=np.bool_)
            lf["shape"].Contains_v(loc, ins, n)
            inside |= ins
        f[:, j0:j0 + nj] = ((inside.reshape(nj, nr, nphi) * rr[None]).sum(axis=2) / rsum[None]).T
    return f


def fraction_volume(f, r_edges, z_edges):
    return float(np.sum(f * (math.pi * np.diff(r_edges ** 2))[:, None] * np.diff(z_edges)[None, :]))


def _ray(shape, R, T, p, d, inside):
    pl = np.ascontiguousarray((p - T) @ R)
    dl = np.ascontiguousarray(d @ R)
    return shape.DistFromInside(pl, dl, 3) if inside else shape.DistFromOutside(pl, dl, 3)


def ring_from_copies(copies, conv, k_poly, anchor, bar_mode, nphi=2048, hr=0.1, hz=0.5):
    """equivalent ring of a set of copies of an off-axis volume (e.g. 8 spacer bars)"""
    ext = []
    for lf in copies:
        sh, R, T = lf["shape"], lf["R"], lf["T"]
        o = sh.GetOrigin()
        c = R @ np.array([o[0], o[1], o[2]], dtype=float) + T
        phi = math.atan2(c[1], c[0])
        d = np.array([math.cos(phi), math.sin(phi), 0.0])
        p0 = np.array([0.0, 0.0, c[2]])
        s_in = _ray(sh, R, T, p0, d, False)
        if not s_in < 1e29:
            raise SystemExit("extract: ray casting missed %s" % lf["node"])
        p1 = p0 + (s_in + 1e-7) * d
        r_out = s_in + 1e-7 + _ray(sh, R, T, p1, d, True)
        cl = np.ascontiguousarray((c - T) @ R)
        if not sh.Contains(cl):
            raise SystemExit("extract: centroid of %s is outside the shape (not a bar-like volume)" % lf["node"])
        z1 = c[2] + _ray(sh, R, T, c, np.array([0.0, 0.0, 1.0]), True)
        z0 = c[2] - _ray(sh, R, T, c, np.array([0.0, 0.0, -1.0]), True)
        corners = (_bbox_corners(sh) @ R.T) + T
        ext.append((s_in, r_out, z0, z1, float(np.max(np.hypot(corners[:, 0], corners[:, 1]))), phi))
    ext = np.array(ext)
    spread = np.max(ext[:, :4], axis=0) - np.min(ext[:, :4], axis=0)
    if np.max(spread) > 1e-4:
        raise SystemExit("extract: copies of %s have different radial/axial extents %s" % (copies[0]["volume"], spread))
    r_in, r_out, z0, z1 = ext[0, :4]
    r_far = float(np.max(ext[:, 4]))
    # azimuthal fraction on an (r, phi, z) grid aligned with the ray-cast extents
    r_lo = max(0.0, r_in - 1.0)
    r_edges = np.linspace(r_lo, r_far + 0.5, int(math.ceil((r_far + 0.5 - r_lo) / hr)) + 1)
    nzb = max(1, int(math.ceil((z1 - z0) / hz)))
    z_edges = np.linspace(z0, z1, nzb + 1)
    f = sample_fraction(copies, r_edges, z_edges, nphi)
    V = fraction_volume(f, r_edges, z_edges)
    # no material outside [z0, z1] and below r_in: one thin slab on each side and the first radial bins
    fz = sample_fraction(copies, r_edges, np.array([z0 - 0.02, z0 - 0.01, z1 + 0.01, z1 + 0.02]), 256)
    leak = float(fz[:, 0].max() + fz[:, 2].max() + f[r_edges[1:] <= r_in - 1e-6].sum())
    cap = float(sum(lf["shape"].Capacity() for lf in copies))
    k = k_poly if (conv == "equal-area" and k_poly) else 1.0
    rm_in, rm_out = k * r_in, k * r_out
    V_ring = math.pi * (rm_out ** 2 - rm_in ** 2) * (z1 - z0)
    cov = V / V_ring
    if bar_mode == "full":
        za, zb = z0, z1
    else:
        L = cov * (z1 - z0)
        if anchor == "center":
            zm = 0.5 * (z0 + z1)
            za, zb = zm - L / 2, zm + L / 2
        elif (anchor == "outer") == (abs(z0) > abs(z1)):
            za, zb = z0, z0 + L
        else:
            za, zb = z1 - L, z1
    info = dict(n_copies=len(copies), ray_r_in=r_in, ray_r_out=r_out, ray_z=[z0, z1], r_far=r_far,
                phi_deg=[math.degrees(p) for p in ext[:, 5]], sampled_volume_cm3=V, capacity_cm3=cap,
                sampled_over_capacity=V / cap if cap else None, leak_check=leak,
                model_radii=[rm_in, rm_out], coverage=cov, model_z=[za, zb], bar_mode=bar_mode, anchor=anchor,
                grid=dict(hr=hr, hz=(z1 - z0) / nzb, nphi=nphi))
    return Region.rect(rm_in, rm_out, za, zb), info


def _apply_snap(region, zvals, rvals, tol):
    """snap the outline vertices of a region to nearby canonical values; returns (region, list of shifts)"""
    loops = region.loops()
    zs = sorted({round(z, 9) for L in loops for _, z in L})
    rs = sorted({round(r, 9) for L in loops for r, _ in L})

    def mapping(vals, canon):
        m = {}
        for v in vals:
            if not canon:
                break
            c = _nearest(canon, v)
            if 1e-9 < abs(c - v) <= tol:
                m[v] = c
        # a piece must not lose a dimension: drop snaps that would merge two of its own values
        tgt = defaultdict(list)
        for v in vals:
            tgt[m.get(v, v)].append(v)
        for t, vs in tgt.items():
            if len(vs) > 1:
                for v in vs:
                    m.pop(v, None)
        return m

    mz, mr = mapping(zs, zvals), mapping(rs, rvals)
    if not mz and not mr:
        return region, []
    new = [[(mr.get(round(r, 9), r), mz.get(round(z, 9), z)) for r, z in L] for L in loops]
    shifts = [("z", v, c) for v, c in mz.items()] + [("r", v, c) for v, c in mr.items()]
    return Region.from_loops_rz(new), shifts


def _boxes_near(b1, b2, tol):
    return not (b1[1] < b2[0] - tol or b2[1] < b1[0] - tol or b1[3] < b2[2] - tol or b2[3] < b1[2] - tol)


def merge_iron(pieces, min_contact=0.1):
    """greedy merge of touching iron pieces as long as the union is one simple polygon: largest contact first,
    equal contacts -> smallest merged area first (keeps symmetric structures symmetric)"""
    groups = [dict(names=[p["name"]], region=p["model"]) for p in pieces]
    perim = [g["region"].perimeter() for g in groups]
    log = []
    while True:
        best = None
        for i in range(len(groups)):
            bi = groups[i]["region"].bbox()
            for j in range(i + 1, len(groups)):
                if not _boxes_near(bi, groups[j]["region"].bbox(), 1e-6):
                    continue
                u = groups[i]["region"] | groups[j]["region"]
                contact = 0.5 * (perim[i] + perim[j] - u.perimeter())
                if contact < min_contact:
                    continue
                polys = u.polygons()
                if len(polys) != 1 or polys[0][1]:
                    continue
                if (best is None or contact > best[0] + 1e-6
                        or (abs(contact - best[0]) <= 1e-6 and u.area() < best[3].area() - 1e-6)):
                    best = (contact, i, j, u)
        if best is None:
            break
        contact, i, j, u = best
        log.append(dict(merged=groups[i]["names"] + groups[j]["names"], contact_cm=contact))
        groups[i] = dict(names=groups[i]["names"] + groups[j]["names"], region=u)
        perim[i] = u.perimeter()
        del groups[j], perim[j]
    return groups, log


def build_main_coil(winding, par):
    """current sheets from the DD4hep winding volume and the parameter set"""
    polys = winding.polygons()
    if len(polys) != 1 or polys[0][1] or len(polys[0][0]) != 4:
        raise SystemExit("extract: the main coil winding is not a rectangle in r-z")
    r0, r1, z0, z1 = winding.bbox()
    nl = len(par["layer_names"])
    edges_ref = par["section_edges_ref_cm"]
    a, b = edges_ref[0], edges_ref[-1]
    zed = [z0 + (e - a) / (b - a) * (z1 - z0) for e in edges_ref]
    w = (r1 - r0) / nl
    t = par.get("sheet_thickness_cm")
    coils = []
    for i, lname in enumerate(par["layer_names"]):
        c = r0 + (i + 0.5) * w
        ra, rb = (c - t / 2, c + t / 2) if t else (r0 + i * w, r0 + (i + 1) * w)
        for j, sname in enumerate(par["section_names"]):
            turns = par["turns_per_layer"][j]
            coils.append(dict(name="%s coil %s" % (sname, lname), mat=1, cur=turns * par["current_A"],
                              region=Region.rect(ra, rb, zed[j], zed[j + 1]), turns=turns,
                              comments=["main coil: %d turns x %.1f A, layer %d of %d of the DD4hep winding "
                                        "(sheet %s cm at the layer mid-radius %.3f)"
                                        % (turns, par["current_A"], i + 1, nl, t if t else "full", c)]))
    info = dict(winding_r=[r0, r1], winding_z=[z0, z1], layer_width=w, sheet_thickness=t, section_edges=zed,
                total_turns=nl * sum(par["turns_per_layer"]),
                total_At=nl * sum(par["turns_per_layer"]) * par["current_A"])
    return coils, info


def run_extract(args):
    import ROOT  # noqa: F401  (needed before dd4hep for the TGeo classes)
    import dd4hep

    t_start = datetime.datetime.now()
    params = json.load(open(args.coil_params))
    compact = os.path.abspath(args.compact)
    desc = dd4hep.Detector.getInstance()
    if args.no_field:
        load = without_fields(compact)
        desc.fromXML(load)
        shutil.rmtree(os.path.dirname(load), ignore_errors=True)
    else:
        desc.fromXML(compact)
    leaves = collect_leaves(desc, args.detector)
    warnings, notes = [], []

    # ---------------- classification
    main_re = re.compile(params["main_coil"]["volume_regex"])
    aux = params.get("aux_coils", [])
    aux_ni = {}
    for s in args.aux_ni:
        k, _, v = s.partition("=")
        aux_ni[k] = float(v)
    iron, mains, auxs, other = [], [], [], []
    for lf in leaves:
        nm, dens, fe, zs = lf["material"]
        lf["name"] = _short(lf["det"], lf["volume"])
        if fe >= args.iron_min_fe:
            iron.append(lf)
        elif main_re.search(lf["volume"]):
            mains.append(lf)
        elif any(re.search(a["volume_regex"], lf["volume"]) for a in aux):
            auxs.append(lf)
        else:
            other.append(lf)
    for lf in other:
        if re.search(r"(?i)winding|coil$", lf["name"]):
            warnings.append("conductor-like volume %s (%s) has no current in the parameter set: treated as air"
                            % (lf["node"], lf["material"][0]))
    if len(mains) != 1:
        raise SystemExit("extract: expected exactly one main coil volume matching %r, found %d"
                         % (params["main_coil"]["volume_regex"], len(mains)))

    # polygon frame of each DetElement (for equal-area rings)
    k_poly = {}
    for lf in leaves:
        pl = _shape_planes(lf["shape"])
        if pl is not None and pl[3]:
            k_poly.setdefault(lf["det"], set()).add(pl[3])
    k_poly = {d: k_equal_area(sorted(s)[0]) for d, s in k_poly.items() if len(s) == 1}

    # ---------------- iron pieces
    pieces, groups = [], defaultdict(list)
    for lf in iron:
        rz = shape_rz(lf["shape"], lf["R"], lf["T"], args.convention)
        if rz is not None:
            model, eq, d, poly = rz
            pieces.append(dict(name=lf["name"], det=lf["det"], nodes=[lf["node"]], material=lf["material"][0],
                               density=lf["material"][1], method="analytic " + d, polygonal=poly, model=model,
                               eq=eq, leaves=[lf]))
        else:
            groups[(lf["det"], lf["volume"])].append(lf)
    for (det, vol), copies in groups.items():
        reg, info = ring_from_copies(copies, args.convention, k_poly.get(det), args.ring_anchor, args.bars)
        if info["leak_check"] > 1e-9:
            warnings.append("%s: material found outside the ray-cast extents (%.3g): ring is approximate"
                            % (vol, info["leak_check"]))
        pieces.append(dict(name=_short(det, vol), det=det, nodes=[c["node"] for c in copies],
                           material=copies[0]["material"][0], density=copies[0]["material"][1],
                           method="sampled ring (%d copies, %s)" % (len(copies), copies[0]["shape"].ClassName()),
                           polygonal=False, model=reg, eq=None, leaves=copies, ring=info))

    # ---------------- geometric sampling of every iron piece (volume check + azimuthal fraction map)
    hr, hz, nphi = args.sampling
    fmap = None
    if not args.no_sampling:
        rmax_all = 0.0
        zr = [1e9, -1e9]
        for p in pieces:
            b = p["model"].bbox()
            if p["eq"] is not None:
                be = p["eq"].bbox()
                b = (min(b[0], be[0]), max(b[1], be[1]), min(b[2], be[2]), max(b[3], be[3]))
            far = max(_rbound(lf["shape"], lf["R"], lf["T"]) for lf in p["leaves"])
            p["_sbox"] = (max(0.0, b[0] - 1.0), far + 1.0, b[2] - 1.0, b[3] + 1.0)
            if "ring" in p:
                p["_sbox"] = (max(0.0, p["ring"]["ray_r_in"] - 1.0), p["ring"]["r_far"] + 1.0,
                              p["ring"]["ray_z"][0] - 1.0, p["ring"]["ray_z"][1] + 1.0)
            rmax_all = max(rmax_all, p["_sbox"][1])
            zr = [min(zr[0], p["_sbox"][2]), max(zr[1], p["_sbox"][3])]
        R0, Z0 = 0.0, math.floor(zr[0] / hz) * hz
        nR = int(math.ceil(rmax_all / hr))
        nZ = int(math.ceil((zr[1] - Z0) / hz))
        fmap = np.zeros((nR, nZ), dtype=np.float32)
        for p in pieces:
            r_lo, r_hi, z_lo, z_hi = p["_sbox"]
            i0, i1 = int(math.floor(r_lo / hr)), int(math.ceil(r_hi / hr))
            j0, j1 = int(math.floor((z_lo - Z0) / hz)), int(math.ceil((z_hi - Z0) / hz))
            re_ = R0 + hr * np.arange(i0, i1 + 1)
            ze_ = Z0 + hz * np.arange(j0, j1 + 1)
            f = sample_fraction(p["leaves"], re_, ze_, nphi)
            fmap[i0:i1, j0:j1] += f.astype(np.float32)
            p["sampled_volume"] = fraction_volume(f, re_, ze_)
        if args.sample_out:
            np.savez_compressed(args.sample_out, f=fmap, r_edges=R0 + hr * np.arange(nR + 1),
                                z_edges=Z0 + hz * np.arange(nZ + 1), nphi=nphi,
                                note="azimuthal iron fraction from TGeoShape::Contains_v, sum over iron pieces")
    for p in pieces:
        p["real_volume"] = p["eq"].volume() if p["eq"] is not None else p["ring"]["sampled_volume_cm3"]
        p["model_volume_raw"] = p["model"].volume()

    # ---------------- snapping (close small gaps between neighbouring pieces) and overlap check
    snaps = []
    if args.snap > 0:
        done = []
        for p in pieces:
            b = p["model"].bbox()
            near = [q for q in done if _boxes_near(b, q["model"].bbox(), args.snap)]
            zc = sorted({round(z, 9) for q in near for L in q["model"].loops() for _, z in L})
            rc = sorted({round(r, 9) for q in near for L in q["model"].loops() for r, _ in L})
            reg, sh = _apply_snap(p["model"], zc, rc, args.snap)
            if sh:
                dv = reg.volume() - p["model"].volume()
                snaps.append(dict(piece=p["name"], shifts=[dict(coord=c, old=o, new=n) for c, o, n in sh],
                                  dvolume_cm3=dv))
                p["model"] = reg
            done.append(p)
    for i in range(len(pieces)):
        for j in range(i + 1, len(pieces)):
            a = (pieces[i]["model"] & pieces[j]["model"]).area()
            if a > 1e-6:
                warnings.append("iron pieces %s and %s overlap in r-z (%.4f cm2)" % (pieces[i]["name"],
                                                                                   pieces[j]["name"], a))

    # ---------------- coils
    winding = shape_rz(mains[0]["shape"], mains[0]["R"], mains[0]["T"], "flat")
    if winding is None:
        raise SystemExit("extract: the main coil volume is not a coaxial solid of revolution")
    coils, coil_info = build_main_coil(winding[0], params["main_coil"])
    for lf in auxs:
        a = next(a for a in aux if re.search(a["volume_regex"], lf["volume"]))
        rz = shape_rz(lf["shape"], lf["R"], lf["T"], "flat")
        if rz is None:
            raise SystemExit("extract: auxiliary coil %s is not coaxial" % lf["node"])
        ni = aux_ni.get(a["name"], a.get("NI_At", 0.0))
        b = rz[0].bbox()
        coils.append(dict(name="%s (DD4hep %s)" % (a["name"], lf["name"]), mat=1, cur=ni, region=rz[0],
                          comments=["auxiliary coil from DD4hep volume %s (%s), r %.3f-%.3f z %.3f..%.3f cm, "
                                    "NI = %.1f A-turns (J = %.2f A/cm2)" % (lf["node"], lf["material"][0], b[0],
                                                                           b[1], b[2], b[3], ni,
                                                                           ni / rz[0].area())]
                          + (["current from the parameter set: " + a.get("source", "")] if ni == 0 else [])))
        if ni == 0:
            notes.append("%s: NI = 0 (current unknown, region kept as air; set --aux-ni %s=<A-turns>)"
                         % (a["name"], a["name"]))
    commented = []
    for spec in args.addon_coil + args.commented_addon:
        nm, _, opt = spec.partition(":")
        if nm not in params.get("addon_coils", {}):
            raise SystemExit("extract: unknown add-on coil %s (parameter set has %s)"
                             % (nm, ", ".join(params.get("addon_coils", {}))))
        a = params["addon_coils"][nm]
        reg = Region.rect(a["r_cm"][0], a["r_cm"][1], a["z_cm"][0], a["z_cm"][1])
        ni = a["NI_At"]
        for kv in filter(None, opt.split(",")):
            k, _, v = kv.partition("=")
            if k == "NI":
                ni = float(v)
            elif k == "J":
                ni = float(v) * reg.area()
            else:
                raise SystemExit("extract: add-on option %s not understood (use NI=... or J=...)" % kv)
        c = dict(name="%s (add-on, not in DD4hep)" % nm, mat=1, cur=ni, region=reg,
                 comments=["add-on coil (NOT in the DD4hep geometry): r %.3f-%.3f z %.3f..%.3f cm, NI = %.1f "
                           "A-turns, J = %.2f A/cm2" % (a["r_cm"][0], a["r_cm"][1], a["z_cm"][0], a["z_cm"][1], ni,
                                                         ni / reg.area()),
                           "source: " + a.get("source", "")])
        if spec in args.addon_coil:
            coils.append(c)
        else:
            commented.append(c)

    iron_union = union_all([p["model"] for p in pieces])
    for c in coils + commented:
        ov = (c["region"] & iron_union).area()
        c["overlap_iron_cm2"] = ov
        if ov > 1e-6:
            hit = [p["name"] for p in pieces if (c["region"] & p["model"]).area() > 1e-6]
            msg = "coil %s overlaps iron (%s) by %.3f cm2" % (c["name"], ", ".join(hit), ov)
            if c in commented:
                notes.append("commented-out " + msg)
            elif args.allow_overlap:
                warnings.append(msg + " (allowed: the coil region is written after the iron and overrides it)")
            else:
                raise SystemExit("extract: " + msg + " (use --allow-overlap to force)")
    for i in range(len(coils)):
        for j in range(i + 1, len(coils)):
            if (coils[i]["region"] & coils[j]["region"]).area() > 1e-6:
                raise SystemExit("extract: coils %s and %s overlap" % (coils[i]["name"], coils[j]["name"]))

    # ---------------- merging
    if args.no_merge:
        groups_iron = [dict(names=[p["name"]], region=p["model"]) for p in pieces]
        merge_log = []
    else:
        groups_iron, merge_log = merge_iron(pieces)
    iron_regions = []
    for g in groups_iron:
        for outer, holes in g["region"].polygons():
            if holes:
                raise SystemExit("extract: internal error, merged iron region with a hole")
            iron_regions.append(dict(names=g["names"], pts=_oriented_rz(outer)))
    iron_regions.sort(key=lambda r: (min(z for _, z in r["pts"]), min(rr for rr, _ in r["pts"])))

    # ---------------- small features (report)
    small = []
    for reg in iron_regions:
        P = reg["pts"]
        for i in range(len(P)):
            L = math.hypot(P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1])
            if L < args.min_feature:
                small.append(dict(region=" + ".join(reg["names"]), edge=[P[i - 1], P[i]], length=L))

    # ---------------- problem box
    allb = [p["model"].bbox() for p in pieces] + [c["region"].bbox() for c in coils]
    zmin, zmax = min(b[2] for b in allb), max(b[3] for b in allb)
    rmax = max(b[1] for b in allb)
    if args.box_z:
        bz = list(args.box_z)
    else:
        bz = [min(-500.0, math.floor((zmin - 200.0) / 100.0) * 100.0),
              max(700.0, math.ceil((zmax + 200.0) / 100.0) * 100.0)]
    box = dict(r=[0.0, args.box_r], z=bz)
    if rmax > args.box_r - 1 or zmin < bz[0] + 1 or zmax > bz[1] - 1:
        raise SystemExit("extract: the magnet does not fit in the problem box %s" % box)

    # ---------------- provenance
    files = compact_files(compact)
    if args.no_field:
        files = [f for f in files if "/fields/" not in f]
    gi = git_info(compact)
    dirty = git_dirty(gi["root"], files) if gi else []
    try:
        ddv = str(ROOT.dd4hep.versionString())
    except Exception:  # noqa: BLE001
        ddv = "unknown"
    rootv = str(ROOT.gROOT.GetVersion())
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")

    total_real = sum(p["real_volume"] * p["density"] for p in pieces) / 1e6
    total_model = sum(p["model"].volume() for p in pieces) * STEEL_DENSITY / 1e6

    # ---------------- deck
    name = args.name or os.path.splitext(os.path.basename(args.output))[0]
    title = args.title or "SoLID CLEO-II r-z model %s (DD4hep %s)" % (name, os.path.basename(compact))
    H = []
    H.append("Generated by magnet_field/extract_rz_model.py on %s (DD4hep %s, ROOT %s)" % (now, ddv, rootv))
    root = gi["root"] if gi else None

    def rel(f):
        return portable_path(f, root)

    H.append("Compact: %s%s" % (rel(compact), " (loaded without its field fragments, --no-field)" if args.no_field
                                else ""))
    if gi:
        H.append("Repository: branch %s, commit %s (%s); %s." % (
            gi["branch"], gi["commit"][:12], gi["date"],
            "uncommitted changes in " + " ".join(dirty) if dirty else
            "the %d compact files read are unmodified" % len(files)))
    H.append("Magnet/field fragments read (sha256): " + "; ".join(
        "%s %s" % (rel(f), _sha(f)[:12]) for f in files
        if re.search(r"/(magnet|fields?)/", f) or re.search(r"magnet|field", os.path.basename(f))
        or os.path.realpath(f) == os.path.realpath(compact)) + ".")
    H.append("DetElements %s: %d iron pieces, %d coil regions. Non-magnetic (air here): %s." % (
        ", ".join(args.detector), len(pieces), len(coils),
        ", ".join("%s (%s)" % (lf["name"], lf["material"][0]) for lf in other) or "none"))
    H.append("")
    H.append("Conventions (magnet_field/models/README.md): x = r, y = z in cm, hall frame (origin at the coil "
             "centre). mat=1 with cur = coil region (cur = total ampere-turns of the region), mat=2 = iron with "
             "POISSON's internal permeability table (mode=0, no MT table). Polygonal iron at its %s. %s Gaps and "
             "offsets below %.2f cm between neighbouring iron pieces are closed by snapping; touching pieces are "
             "merged into simple polygons. Iron = Fe mass fraction >= %.2f; stainless steel, aluminium, the "
             "conductor and the cryostat are air." % (
                 "apothem (flat) radii, as CLEOv8/9/10.am" if args.convention == "flat" else
                 "equal-area radii (apothem x sqrt(n tan(pi/n)/pi))",
                 "Off-axis bars become a full ring of the same iron volume, its axial length scaled by the "
                 "azimuthal coverage (%s end fixed)." % args.ring_anchor if args.bars == "length" else
                 "Off-axis bars become a full ring (no coverage scaling).", args.snap, args.iron_min_fe))
    H.append("Main coil: %d turns, %.0f A-turns. The DD4hep winding (r %.2f-%.2f, z %.2f..%.2f cm) is split into "
             "%d equal layers; each carries its current in a %s cm sheet at the layer mid-radius; section edges "
             "at z = %s cm (CLEOv10.am grading)." % (
                 coil_info["total_turns"], coil_info["total_At"], coil_info["winding_r"][0],
                 coil_info["winding_r"][1], coil_info["winding_z"][0], coil_info["winding_z"][1],
                 len(params["main_coil"]["layer_names"]), coil_info["sheet_thickness"],
                 ", ".join("%.3f" % z for z in coil_info["section_edges"])))
    H.append("Iron mass: DD4hep %.1f t (analytic or sampled volumes x material densities); this r-z model "
             "%.1f t at %.2f g/cm3." % (total_real, total_model, STEEL_DENSITY))
    for p in pieces:
        if "ring" in p:
            ri = p["ring"]
            H.append("%s: %d copies -> ring r %.3f-%.3f, z %.3f..%.3f; coverage %.4f (sampled volume / TGeo "
                     "capacity %.5f)." % (p["name"], ri["n_copies"], ri["model_radii"][0], ri["model_radii"][1],
                                          ri["model_z"][0], ri["model_z"][1], ri["coverage"],
                                          ri["sampled_over_capacity"] or 0))
    for s_ in snaps:
        H.append("Snapped %s: %s." % (s_["piece"], ", ".join("%s %.3f -> %.3f" % (x["coord"], x["old"], x["new"])
                                                            for x in s_["shifts"])))
    for w in warnings:
        H.append("WARNING: " + w)
    for n_ in notes:
        H.append("NOTE: " + n_)
    regions = []
    for c in coils:
        polys = c["region"].polygons()
        regions.append(dict(name=c["name"], mat=1, cur=c["cur"], pts=_oriented_rz(polys[0][0]),
                            comments=c.get("comments", [])))
    ngroup = 0
    for reg in iron_regions:
        if len(reg["names"]) > 1:
            ngroup += 1
        short = reg["names"][0] if len(reg["names"]) == 1 else "Iron group %d" % ngroup
        cm = ["iron: " + ", ".join(reg["names"])] if len(reg["names"]) > 1 else []
        regions.append(dict(name=short, mat=2, cur=0.0, pts=reg["pts"], comments=cm, pieces=reg["names"]))
    commented_blocks = []
    for c in commented:
        polys = c["region"].polygons()
        pts = _oriented_rz(polys[0][0])
        blk = ["---- optional, NOT active: %s ----" % c["name"]] + c["comments"]
        blk += ["overlaps the DD4hep iron by %.2f cm2 in this geometry" % c["overlap_iron_cm2"]] \
            if c["overlap_iron_cm2"] > 1e-6 else []
        blk.append("&reg    mat=1,   cur=%.1f &" % c["cur"])
        blk += ["&po x=%9s, y=%9s &" % (_fmt(r), _fmt(z)) for r, z in pts + [pts[0]]]
        commented_blocks.append(blk)
    zones = radial_zones(box["r"][1])
    refine = thin_radial_features([p["model"] for p in pieces] + [c["region"] for c in coils], zones)
    zones = radial_zones(box["r"][1], refine)
    if refine:
        H.insert(H.index("") if "" in H else len(H), "Radial mesh refined around features thinner than 2 cells: %s."
                 % ", ".join("r %.2f-%.2f at %.2f cm" % r_ for r_ in refine))
    model = dict(title=title, header=H, box=box, regions=regions, commented_regions=commented_blocks,
                 mesh_zones=zones)
    write_am(args.output, model)

    # ---------------- sidecar
    side = dict(
        title=title, generated=now, generator="magnet_field/extract_rz_model.py", dd4hep=ddv, root=rootv,
        compact=rel(compact), git={k: v for k, v in (gi or {}).items() if k != "root"} or None, dirty_files=dirty,
        paths="relative to the repository root; ${VAR}/... = under that environment variable (container)",
        files=[dict(path=rel(f), sha256=_sha(f),
                    geometry_sha256=geometry_digest(f, args.no_field and os.path.realpath(f) == os.path.realpath(
                        compact))) for f in files]
        + [dict(path=src, sha256=_sha(os.path.join(root, src)), role="plugin source") for src in PLUGIN_SOURCES
           if root and os.path.isfile(os.path.join(root, src))],
        options=dict(no_field=args.no_field, convention=args.convention, bars=args.bars, ring_anchor=args.ring_anchor,
                     snap=args.snap,
                     merge=not args.no_merge, iron_min_fe=args.iron_min_fe, sampling=list(args.sampling),
                     coil_params=rel(args.coil_params), aux_ni=aux_ni, addon_coil=args.addon_coil,
                     commented_addon=args.commented_addon),
        box=box,
        main_coil=coil_info,
        regions=[dict(name=r["name"], mat=r["mat"], cur=r["cur"], pts=r["pts"], pieces=r.get("pieces"))
                 for r in regions],
        pieces=[dict(name=p["name"], det=p["det"], nodes=p["nodes"], material=p["material"], density=p["density"],
                     method=p["method"], polygonal=p["polygonal"],
                     loops=[_oriented_rz(L) for L in p["model"].loops()],
                     real_volume_cm3=p["real_volume"], real_mass_t=p["real_volume"] * p["density"] / 1e6,
                     model_volume_cm3=p["model"].volume(),
                     sampled_volume_cm3=p.get("sampled_volume"),
                     sampled_over_real=(p["sampled_volume"] / p["real_volume"]) if p.get("sampled_volume") else None,
                     ring=p.get("ring"))
                for p in pieces],
        coils=[dict(name=c["name"], cur=c["cur"], bbox=c["region"].bbox(), overlap_iron_cm2=c["overlap_iron_cm2"])
               for c in coils],
        commented_coils=[dict(name=c["name"], cur=c["cur"], bbox=c["region"].bbox(),
                              overlap_iron_cm2=c["overlap_iron_cm2"]) for c in commented],
        non_magnetic=[dict(node=lf["node"], material=lf["material"][0], fe_fraction=lf["material"][2])
                      for lf in other],
        snaps=snaps, merges=merge_log, small_features=small, warnings=warnings, notes=notes,
        iron_mass_real_t=total_real, iron_mass_model_t=total_model,
        runtime_s=(datetime.datetime.now() - t_start).total_seconds())
    jpath = args.json or os.path.splitext(args.output)[0] + ".json"
    with open(jpath, "w") as f:
        json.dump(side, f, indent=1, default=float)
    # ---------------- console summary
    print("extract_rz_model: %s -> %s (%d regions) and %s" % (os.path.basename(compact), args.output,
                                                               len(regions) + 1, jpath))
    print("  %-26s %-48s %11s %11s %11s %8s" % ("piece", "method", "V_real", "V_model", "V_sampled", "s/real"))
    for p in pieces:
        print("  %-26s %-48s %11.0f %11.0f %11s %8s" % (
            p["name"], p["method"][:48], p["real_volume"], p["model"].volume(),
            "%.0f" % p["sampled_volume"] if p.get("sampled_volume") else "-",
            "%.5f" % (p["sampled_volume"] / p["real_volume"]) if p.get("sampled_volume") else "-"))
    print("  iron mass: DD4hep %.2f t, model %.2f t" % (total_real, total_model))
    for s in snaps:
        print("  snap %s: %s" % (s["piece"], s["shifts"]))
    for m in merge_log:
        print("  merge %s (contact %.2f cm)" % (" + ".join(m["merged"]), m["contact_cm"]))
    for w in warnings:
        print("  WARNING: " + w)
    for n_ in notes:
        print("  NOTE: " + n_)
    for s in small:
        print("  small feature: %s edge %.3f cm at %s" % (s["region"], s["length"], s["edge"]))
    return 0


def _rbound(shape, R, T):
    """upper bound of the distance from the z axis of any point of a placed shape"""
    if _coaxial(R, T):
        pl = _shape_planes(shape)
        if pl is not None:
            return max(pl[2]) / (math.cos(math.pi / pl[3]) if pl[3] else 1.0)
        if shape.ClassName() == "TGeoCompositeShape":
            bn = shape.GetBoolNode()
            RL, TL = _mat(bn.GetLeftMatrix())
            RR, TR = _mat(bn.GetRightMatrix())
            left = _rbound(bn.GetLeftShape(), R @ RL, R @ TL + T)
            if int(bn.GetBooleanOperator()) == 2:
                return left
            right = _rbound(bn.GetRightShape(), R @ RR, R @ TR + T)
            return max(left, right) if int(bn.GetBooleanOperator()) == 0 else min(left, right)
    c = (_bbox_corners(shape) @ R.T) + T
    return float(np.max(np.hypot(c[:, 0], c[:, 1])))


# =================================================================================================================
# comparison with reference decks
# =================================================================================================================
def load_model(path, key=None):
    """-> dict(name, regions=[dict(name, mat, cur, region, pts)], pieces=[dict(name, region)] or None, box)"""
    if path.endswith(".am"):
        d = parse_am(open(path).read())
        regs = []
        for r in d["regions"][1:]:
            regs.append(dict(name=r["name"], mat=r["mat"], cur=r["cur"], region=Region.from_loops_rz([r["pts"]]),
                             pts=r["pts"]))
        box = d["regions"][0]["pts"]
        return dict(name=os.path.basename(path), title=d["title"], regions=regs, pieces=None, box=box)
    if path.endswith((".html", ".htm")):  # "CLEO-II Stage 1 and 2.html": data block <script id="data">
        m = re.search(r"<script[^>]*\bid=\"data\"[^>]*>(.*?)</script>", open(path, encoding="utf-8").read(), re.S)
        if not m:
            raise SystemExit("compare: no <script id=\"data\"> JSON block in %s" % path)
        j = json.loads(m.group(1))
    else:
        j = json.load(open(path))
    if "regions" in j and "pieces" in j:  # sidecar of this script
        regs = [dict(name=r["name"], mat=r["mat"], cur=r["cur"], region=Region.from_loops_rz([r["pts"]]),
                     pts=r["pts"]) for r in j["regions"]]
        pieces = [dict(name=p["name"], region=Region.from_loops_rz(p["loops"])) for p in j["pieces"]]
        return dict(name=os.path.basename(path), title=j["title"], regions=regs, pieces=pieces, box=None, side=j)
    # Stage 1/2 html data: list of regions with pts as (z, r)
    src = j
    for k in (key or "geometry").split("."):
        src = src[k]
    if isinstance(src, dict):
        src = src["geometry"]
    regs = []
    for r in src:
        pts = [(p[1], p[0]) for p in r["pts"]]
        if pts[0] == pts[-1]:
            pts = pts[:-1]
        regs.append(dict(name=r["name"], mat=int(r.get("mat", 1)), cur=float(r.get("cur", 0.0)),
                         region=Region.from_loops_rz([pts]), pts=pts))
    return dict(name="%s[%s]" % (os.path.basename(path), key or "geometry"), title=key, regions=regs, pieces=None,
                box=None)


def _components(region):
    out = []
    for outer, holes in region.polygons():
        reg = Region.from_loops_rz([outer] + holes)
        A, P = reg.area(), reg.perimeter()
        out.append(dict(region=reg, area=A, perimeter=P, width=2 * A / P if P else 0.0, bbox=reg.bbox(),
                        volume=reg.volume()))
    return sorted(out, key=lambda c: -c["area"])


def _hausdorff(la, lb, step=0.25):
    """symmetric Hausdorff distance between two sets of closed polylines (r, z)"""
    def pts(loops):
        out = []
        for L in loops:
            for i in range(len(L)):
                p, q = np.array(L[i - 1], float), np.array(L[i], float)
                n = max(1, int(math.ceil(np.linalg.norm(q - p) / step)))
                out.extend(p + (q - p) * t for t in np.linspace(0, 1, n, endpoint=False))
        return np.array(out)

    if not la or not lb:
        return float("nan")
    return float(max(_dist_to_segments(pts(la), _segments(lb)).max(),
                     _dist_to_segments(pts(lb), _segments(la)).max()))


def _fmt_bbox(b):
    return "r %7.2f-%7.2f  z %8.2f..%8.2f" % b if b else "-"


def _segments(loops):
    return np.array([(L[i - 1][0], L[i - 1][1], L[i][0], L[i][1]) for L in loops for i in range(len(L))])


def _dist_to_segments(P, S):
    a, b = S[:, :2], S[:, 2:]
    ab = b - a
    L2 = np.maximum((ab ** 2).sum(1), 1e-30)
    t = np.clip(((P[:, None, :] - a[None]) * ab[None]).sum(2) / L2[None], 0, 1)
    proj = a[None] + t[..., None] * ab[None]
    return np.sqrt(((P[:, None, :] - proj) ** 2).sum(2)).min(1)


def _inside(region_polys, p):
    for outer, holes in region_polys:
        if _point_in_loop(p, outer) and not any(_point_in_loop(p, h) for h in holes):
            return True
    return False


def edge_offsets(ref_pts, UR_polys, UM_polys, UM_segs, far=10.0):
    """signed offset of the model iron outline from each external edge of a reference region:
    + = model iron extends beyond the reference edge, - = model boundary lies inside the reference region"""
    pts = list(ref_pts)
    out = []
    for i in range(len(pts)):
        p, q = np.array(pts[i - 1], float), np.array(pts[i], float)
        d = q - p
        L = float(np.hypot(*d))
        if L < 1e-6:
            continue
        n = np.array([d[1], -d[0]]) / L
        m = 0.5 * (p + q)
        if _point_in_loop(tuple(m + 1e-3 * n), pts):
            n = -n
        if _inside(UR_polys, tuple(m + 1e-3 * n)):
            continue  # internal edge: the reference has iron on both sides
        ts = np.linspace(0.1, 0.9, 9) if L > 1.0 else np.array([0.5])
        S = p[None] + ts[:, None] * d[None]
        dist = _dist_to_segments(S, UM_segs)
        sign = np.array([1.0 if _inside(UM_polys, tuple(s)) else -1.0 for s in S])
        off = sign * dist
        if abs(d[0]) < 1e-6:
            kind = "r=%.2f [z %.1f..%.1f]" % (p[0], min(p[1], q[1]), max(p[1], q[1]))
        elif abs(d[1]) < 1e-6:
            kind = "z=%.2f [r %.1f-%.1f]" % (p[1], min(p[0], q[0]), max(p[0], q[0]))
        else:
            kind = "(%.1f,%.1f)-(%.1f,%.1f)" % (p[0], p[1], q[0], q[1])
        out.append(dict(edge=(tuple(p), tuple(q)), kind=kind, length=L, off_min=float(off.min()),
                        off_max=float(off.max()), far=bool(np.all(dist > far))))
    return out


def _fmt_off(e):
    if e["far"]:
        return "%s: no model iron within 10 cm" % e["kind"]
    a, b = e["off_min"], e["off_max"]
    return "%s: %+.2f" % (e["kind"], a) if abs(a - b) < 0.005 else "%s: %+.2f..%+.2f" % (e["kind"], a, b)


def run_compare(args):
    M = load_model(args.model, args.model_key)
    Rf = load_model(args.reference, args.ref_key)
    out = []
    p = out.append
    p("# r-z model comparison: %s vs reference %s" % (M["name"], Rf["name"]))
    p("")
    m_iron = [r for r in M["regions"] if r["mat"] >= 2]
    r_iron = [r for r in Rf["regions"] if r["mat"] >= 2]
    m_coil = [r for r in M["regions"] if r["mat"] == 1 and r["cur"] != 0]
    r_coil = [r for r in Rf["regions"] if r["mat"] == 1 and r["cur"] != 0]
    r_air = [r for r in Rf["regions"] if r["mat"] == 1 and r["cur"] == 0]
    m_air = [r for r in M["regions"] if r["mat"] == 1 and r["cur"] == 0]
    UM, UR = union_all([r["region"] for r in m_iron]), union_all([r["region"] for r in r_iron])
    UMp, URp = UM.polygons(), UR.polygons()
    UMs = _segments(UM.loops())
    pieces = M["pieces"] or [dict(name=r["name"], region=r["region"]) for r in m_iron]
    t = lambda v: v * STEEL_DENSITY / 1e6
    p("## Totals")
    p("")
    p("- iron (r-z body of revolution at %.2f g/cm3): model %.2f t in %d regions, reference %.2f t in %d regions, "
      "common %.2f t" % (STEEL_DENSITY, t(UM.volume()), len(m_iron), t(UR.volume()), len(r_iron),
                         t((UM & UR).volume())))
    p("- iron only in the model %.2f t, only in the reference %.2f t; Hausdorff distance of the outlines %.2f cm"
      % (t((UM - UR).volume()), t((UR - UM).volume()), _hausdorff(UM.loops(), UR.loops())))
    p("- coil ampere-turns: model %.1f, reference %.1f" % (sum(m["cur"] for m in m_coil),
                                                           sum(r["cur"] for r in r_coil)))
    if r_air or m_air:
        p("- non-magnetic regions (mat=1, cur=0): reference: %s; model: %s" % (
            "; ".join("%s [%s]" % (r["name"], _fmt_bbox(r["region"].bbox())) for r in r_air) or "-",
            "; ".join("%s [%s]" % (r["name"], _fmt_bbox(r["region"].bbox())) for r in m_air) or "-"))
    p("")
    p("## Coils (mat=1 with current)")
    p("")
    p("| reference region | reference r, z (cm) | A-turns | model region | model r, z (cm) | A-turns | max shift cm |")
    p("|---|---|---:|---|---|---:|---:|")
    used_m = set()
    for r in r_coil:
        best = max(m_coil, key=lambda m: (r["region"] & m["region"]).area(), default=None)
        if best is None or (r["region"] & best["region"]).area() < 1e-6:
            p("| %s | %s | %.1f | (none) | - | - | - |" % (r["name"], _fmt_bbox(r["region"].bbox()), r["cur"]))
            continue
        used_m.add(id(best))
        d = max(abs(a - b) for a, b in zip(r["region"].bbox(), best["region"].bbox()))
        p("| %s | %s | %.1f | %s | %s | %.1f | %.3f |" % (r["name"], _fmt_bbox(r["region"].bbox()), r["cur"],
                                                         best["name"], _fmt_bbox(best["region"].bbox()),
                                                         best["cur"], d))
    for m in m_coil:
        if id(m) not in used_m:
            p("| (none) | - | - | %s | %s | %.1f | - |" % (m["name"], _fmt_bbox(m["region"].bbox()), m["cur"]))
    p("")
    p("## Iron, region by region")
    p("")
    p("Model pieces overlapping at least half of the reference region or half of the piece. Offsets: signed distance "
      "of the model iron outline from each external edge of the reference region (+ = model iron extends beyond the "
      "edge, - = model outline inside the region); edges shared with other reference iron are skipped.")
    p("")
    p("| reference region | reference r, z (cm) | model piece(s) | covered | boundary offsets (cm) |")
    p("|---|---|---|---:|---|")
    for r in r_iron:
        A = r["region"].area()
        used = [q["name"] for q in pieces
                if (r["region"] & q["region"]).area() > 0.5 * min(A, q["region"].area())]
        offs = edge_offsets(r["pts"], URp, UMp, UMs)
        p("| %s | %s | %s | %.3f | %s |" % (r["name"], _fmt_bbox(r["region"].bbox()), ", ".join(used) or "(none)",
                                           (r["region"] & UM).area() / A, "; ".join(_fmt_off(e) for e in offs)))
    unmatched = [q["name"] for q in pieces if (q["region"] & UR).area() < 0.5 * q["region"].area()]
    if unmatched:
        p("")
        p("Model pieces with less than half of their area on reference iron: " + ", ".join(unmatched))
    p("")
    p("## Iron in one model only, by model piece (+) and by reference region (-)")
    p("")
    p("| | piece / region | area cm2 | mass t | extent of the difference (cm) |")
    p("|---|---|---:|---:|---|")
    for sign, items, other in (("+", pieces, UR), ("-", r_iron, UM)):
        for it in items:
            D = it["region"] - other
            if D.area() > 1e-3:
                p("| %s | %s | %.1f | %.3f | %s |" % (sign, it["name"], D.area(), t(D.volume()), _fmt_bbox(D.bbox())))
    rows = []
    for sign, D, own, other in (("+", UM - UR, pieces, r_iron), ("-", UR - UM, r_iron, pieces)):
        for c in _components(D):
            share = sorted(((c["region"] & it["region"]).area(), it["name"]) for it in own)
            share = [n for a, n in reversed(share) if a > 1e-6]
            near = [it["name"] for it in other if _boxes_near(c["bbox"], it["region"].bbox(), 0.6)]
            rows.append((sign, c, ", ".join(share) or "-", ", ".join(near) or "-"))
    rows.sort(key=lambda x: -x[1]["area"])
    big = [x for x in rows if x[1]["width"] >= args.tol]
    small = [x for x in rows if x[1]["width"] < args.tol]
    p("")
    p("## Iron differences (connected components of model - reference (+) and reference - model (-))")
    p("")
    p("width = 2 area / perimeter, i.e. the thickness of a sliver. %d components are %.2f cm thick or more, %d are "
      "thinner (boundary shifts, listed in the region table above)." % (len(big), args.tol, len(small)))
    p("")
    p("| | r, z extent (cm) | area cm2 | width cm | mass t | in (this side) | next to (other side) |")
    p("|---|---|---:|---:|---:|---|---|")
    for sign, c, own, oth in big + small:
        p("| %s | %s | %.1f | %.2f | %.3f | %s | %s |" % (sign, _fmt_bbox(c["bbox"]), c["area"], c["width"],
                                                         t(c["volume"]), own, oth))
    text = "\n".join(out)
    print(text)
    if args.report:
        with open(args.report, "w") as f:
            f.write(text + "\n")
    if args.plot:
        plot_compare(M, Rf, m_iron, r_iron, m_coil, r_coil, args.plot,
                     args.title or "%s vs %s" % (M["name"], Rf["name"]), args.zoom, m_air, r_air)
    return 0


# =================================================================================================================
# plots (matplotlib, light theme; palette: dataviz reference slots blue / orange / aqua + neutrals)
# =================================================================================================================
C_SURF, C_INK, C_INK2, C_MUTED, C_GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
C_COMMON, C_MODEL, C_REF, C_COIL = "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"


def _patch(polys, **kw):
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    verts, codes = [], []
    for outer, holes in polys:
        for L in [outer] + holes:
            pts = [(z, r) for r, z in L]
            verts += pts + [pts[0]]
            codes += [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]
    return PathPatch(Path(verts, codes), **kw)


def _style(ax):
    ax.set_facecolor(C_SURF)
    for s in ax.spines.values():
        s.set_color(C_GRID)
        s.set_linewidth(0.8)
    ax.tick_params(colors=C_INK2, labelsize=8, length=3, color=C_GRID)
    ax.grid(True, color=C_GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def _density_windows(regions, extent, wz, wr, n, cell=1.0):
    """n non-overlapping windows (wz x wr cm) with the largest difference area inside"""
    z0, z1, r0, r1 = extent
    nz, nr = int(math.ceil((z1 - z0) / cell)), int(math.ceil((r1 - r0) / cell))
    g = np.zeros((nr, nz))
    for reg in regions:
        for za, zb, l0, l1, h0, h1 in reg.traps:
            for j in range(max(0, int((za - z0) / cell)), min(nz, int(math.ceil((zb - z0) / cell)))):
                ca, cb = max(za, z0 + j * cell), min(zb, z0 + (j + 1) * cell)
                if cb <= ca:
                    continue
                tm = (0.5 * (ca + cb) - za) / (zb - za)
                lo, hi = l0 + (l1 - l0) * tm, h0 + (h1 - h0) * tm
                for i in range(max(0, int((lo - r0) / cell)), min(nr, int(math.ceil((hi - r0) / cell)))):
                    g[i, j] += (min(hi, r0 + (i + 1) * cell) - max(lo, r0 + i * cell)) * (cb - ca)
    kz, kr = max(1, int(wz / cell)), max(1, int(wr / cell))
    c = np.zeros((nr + 1, nz + 1))
    c[1:, 1:] = g.cumsum(0).cumsum(1)
    wins = []
    for _ in range(n):
        if kr > nr or kz > nz:
            break
        s = c[kr:, kz:] - c[:-kr, kz:] - c[kr:, :-kz] + c[:-kr, :-kz]
        if s.max() <= 1e-9:
            break
        i, j = np.unravel_index(np.argmax(s), s.shape)
        wins.append((z0 + j * cell, z0 + (j + kz) * cell, r0 + i * cell, r0 + (i + kr) * cell))
        ia, ib = max(0, i - kr + 1), min(nr, i + kr)
        ja, jb = max(0, j - kz + 1), min(nz, j + kz)
        g[ia:ib, ja:jb] = 0.0
        c[1:, 1:] = g.cumsum(0).cumsum(1)
    return wins


def plot_compare(M, Rf, m_iron, r_iron, m_coil, r_coil, path, title, zooms, m_air=(), r_air=()):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle

    UM, UR = union_all([r["region"] for r in m_iron]), union_all([r["region"] for r in r_iron])
    common, plus, minus = UM & UR, UM - UR, UR - UM
    b = [r["region"].bbox() for r in m_iron + r_iron + m_coil + r_coil]
    z0, z1 = min(x[2] for x in b) - 20, max(x[3] for x in b) + 20
    r1 = max(x[1] for x in b) + 15
    W = 12.0
    top_h = min(6.0, max(2.4, (W - 0.9) * r1 / (z1 - z0)))
    fig = plt.figure(figsize=(W, top_h + 4.3), dpi=80, facecolor=C_SURF)
    gs = fig.add_gridspec(2, 3, height_ratios=[top_h, 3.0], hspace=0.28, wspace=0.16, left=0.06, right=0.99,
                          top=1 - 0.95 / (top_h + 4.3), bottom=0.45 / (top_h + 4.3))
    axes = [fig.add_subplot(gs[0, :])] + [fig.add_subplot(gs[1, k]) for k in range(3)]

    def draw(ax, lw):
        _style(ax)
        if not common.empty():
            ax.add_patch(_patch(common.polygons(), facecolor=C_COMMON, edgecolor="none"))
        if not plus.empty():
            ax.add_patch(_patch(plus.polygons(), facecolor=C_MODEL, edgecolor=C_MODEL, linewidth=0.5 * lw))
        if not minus.empty():
            ax.add_patch(_patch(minus.polygons(), facecolor=C_REF, edgecolor=C_REF, linewidth=0.5 * lw))
        for r in r_iron:
            ax.add_patch(_patch(r["region"].polygons(), facecolor="none", edgecolor=C_REF, linewidth=0.6 * lw))
        for r in m_iron:
            ax.add_patch(_patch(r["region"].polygons(), facecolor="none", edgecolor=C_INK, linewidth=0.45 * lw))
        for c in r_coil:
            ax.add_patch(_patch(c["region"].polygons(), facecolor="none", edgecolor=C_REF, linewidth=0.5 * lw))
        for c in m_coil:
            ax.add_patch(_patch(c["region"].polygons(), facecolor=C_COIL, edgecolor=C_COIL, linewidth=0.4 * lw))
        for c in m_air:
            ax.add_patch(_patch(c["region"].polygons(), facecolor="none", edgecolor=C_COIL, linewidth=0.6 * lw))
        for c in r_air:
            ax.add_patch(_patch(c["region"].polygons(), facecolor="none", edgecolor=C_REF, linewidth=0.5 * lw))

    draw(axes[0], 1.0)
    axes[0].set_xlim(z0, z1)
    axes[0].set_ylim(0, r1)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("z (cm)", color=C_INK2, fontsize=9)
    axes[0].set_ylabel("r (cm)", color=C_INK2, fontsize=9)
    wz, wr = 64.0, 48.0
    wins = [tuple(zz) for zz in (zooms or [])][:3]
    if len(wins) < 3:
        for w_ in _density_windows([plus, minus], (z0, z1, 0.0, r1), wz, wr, 6):
            if len(wins) >= 3:
                break
            if all(w_[1] <= o[0] or w_[0] >= o[1] or w_[3] <= o[2] or w_[2] >= o[3] for o in wins):
                wins.append(w_)
    for k, ax in enumerate(axes[1:]):
        if k >= len(wins):
            ax.set_visible(False)
            continue
        draw(ax, 2.2)
        ax.set_xlim(wins[k][0], wins[k][1])
        ax.set_ylim(wins[k][2], wins[k][3])
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("zoom %d" % (k + 1), color=C_INK2, fontsize=8, loc="left")
        ax.set_xlabel("z (cm)", color=C_INK2, fontsize=8)
        axes[0].add_patch(Rectangle((wins[k][0], wins[k][2]), wins[k][1] - wins[k][0], wins[k][3] - wins[k][2],
                                    fill=False, edgecolor=C_INK2, linewidth=0.8))
        axes[0].text(wins[k][0] + 1, wins[k][3] - 1, "%d" % (k + 1), color=C_INK2, fontsize=8, va="top")
    handles = [Patch(facecolor=C_COMMON, label="iron in both"),
               Patch(facecolor=C_MODEL, label="iron only in the DD4hep r-z model"),
               Patch(facecolor=C_REF, label="iron only in the reference"),
               Line2D([], [], color=C_INK, linewidth=0.8, label="model regions"),
               Line2D([], [], color=C_REF, linewidth=0.8, label="reference regions and coils"),
               Patch(facecolor=C_COIL, label="model coils (no current: outline)")]
    fig.legend(handles=handles, loc="upper left", ncol=6, frameon=False, fontsize=8, labelcolor=C_INK2,
               bbox_to_anchor=(0.055, 1 - 0.42 / (top_h + 4.3)))
    fig.suptitle(title, color=C_INK, fontsize=11, x=0.06, ha="left", y=1 - 0.08 / (top_h + 4.3), va="top")
    fig.savefig(path, facecolor=C_SURF)
    plt.close(fig)


def run_plot_fraction(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    M = load_model(args.model)
    d = np.load(args.npz)
    f, re_, ze_ = d["f"], d["r_edges"], d["z_edges"]
    # downsample to <= 0.5 x 1 cm for a small image
    fr = max(1, int(round(0.5 / (re_[1] - re_[0]))))
    fz = max(1, int(round(1.0 / (ze_[1] - ze_[0]))))
    nr, nz = f.shape[0] // fr * fr, f.shape[1] // fz * fz
    f = f[:nr, :nz].reshape(nr // fr, fr, nz // fz, fz).mean(axis=(1, 3))
    re_, ze_ = re_[:nr + 1:fr], ze_[:nz + 1:fz]
    ramp = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    cmap = LinearSegmentedColormap.from_list("blue_seq", ramp)
    fig, ax = plt.subplots(figsize=(12, 4.8), dpi=80, facecolor=C_SURF)
    fig.subplots_adjust(left=0.06, right=0.93, top=0.84, bottom=0.12)
    _style(ax)
    ax.grid(False)
    im = ax.pcolormesh(ze_, re_, np.clip(f, 0, 1), cmap=cmap, vmin=0, vmax=1, shading="flat", rasterized=True)
    for r in M["regions"]:
        col = C_REF if r["mat"] >= 2 else C_COIL
        ax.add_patch(_patch(r["region"].polygons(), facecolor="none", edgecolor=col, linewidth=0.7))
    b = [r["region"].bbox() for r in M["regions"]]
    ax.set_xlim(min(x[2] for x in b) - 20, max(x[3] for x in b) + 20)
    ax.set_ylim(0, max(x[1] for x in b) + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("z (cm)", color=C_INK2, fontsize=9)
    ax.set_ylabel("r (cm)", color=C_INK2, fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("azimuthal iron fraction", color=C_INK2, fontsize=8)
    cb.ax.tick_params(labelsize=7, colors=C_INK2)
    cb.outline.set_edgecolor(C_GRID)
    fig.legend(handles=[Line2D([], [], color=C_REF, linewidth=1, label="r-z model iron regions (deck)"),
                        Line2D([], [], color=C_COIL, linewidth=1, label="coil regions (deck)")],
               loc="upper left", bbox_to_anchor=(0.055, 0.935), ncol=2, frameon=False, fontsize=8,
               labelcolor=C_INK2)
    fig.suptitle(args.title or M["title"], color=C_INK, fontsize=11, x=0.06, ha="left", y=0.99, va="top")
    fig.savefig(args.output, facecolor=C_SURF)
    plt.close(fig)
    return 0


# =================================================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract", help="DD4hep compact -> POISSON deck (+ JSON sidecar)")
    e.add_argument("compact")
    e.add_argument("-o", "--output", required=True, help="deck to write (.am)")
    e.add_argument("--json", help="sidecar JSON (default: deck name with .json)")
    e.add_argument("--name", help="model name for the title (default: deck basename)")
    e.add_argument("--title", help="deck title line")
    e.add_argument("--detector", action="append", help="DetElement(s) to read (default MagnetYoke, MagnetCoil)")
    e.add_argument("--convention", choices=("flat", "equal-area"), default="flat",
                   help="radius of polygonal (TGeoPgon) iron surfaces: apothem (default, as the POISSON decks) or "
                        "equal-area circle")
    e.add_argument("--bars", choices=("length", "full"), default="length",
                   help="off-axis bars: ring with the axial length scaled by the azimuthal coverage (default), or "
                        "a full ring")
    e.add_argument("--ring-anchor", choices=("outer", "inner", "center"), default="outer",
                   help="end of a shortened ring that stays in place (outer = larger |z|, as the decks)")
    e.add_argument("--snap", type=float, default=0.5, help="close gaps/offsets below this (cm) between "
                   "neighbouring iron pieces by snapping to the earlier piece (0 = off)")
    e.add_argument("--no-merge", action="store_true", help="one region per DD4hep iron piece")
    e.add_argument("--min-feature", type=float, default=0.5, help="report iron edges shorter than this (cm)")
    e.add_argument("--iron-min-fe", type=float, default=0.9, help="minimum Fe mass fraction of a magnetic material")
    e.add_argument("--coil-params", default=DEFAULT_PARAMS, help="current-distribution parameter set (json)")
    e.add_argument("--aux-ni", action="append", default=[], metavar="NAME=AT",
                   help="ampere-turns of an auxiliary coil of the parameter set (e.g. CorrectionCoil=-16000)")
    e.add_argument("--addon-coil", action="append", default=[], metavar="NAME[:NI=..|:J=..]",
                   help="add a coil that is not in DD4hep (parameter set 'addon_coils'); J in A/cm2")
    e.add_argument("--commented-addon", action="append", default=[], metavar="NAME[:NI=..|:J=..]",
                   help="write an add-on coil as a commented-out (inactive) region")
    e.add_argument("--allow-overlap", action="store_true", help="allow coils overlapping iron")
    e.add_argument("--box-r", type=float, default=500.0, help="problem box radius (cm)")
    e.add_argument("--box-z", type=float, nargs=2, metavar=("ZMIN", "ZMAX"),
                   help="problem box z range (default: magnet +- 200 cm rounded to 100 cm, at least -500..700)")
    e.add_argument("--sampling", type=float, nargs=3, default=(0.25, 0.5, 360), metavar=("HR", "HZ", "NPHI"),
                   help="grid of the volume check / fraction map (cm, cm, phi samples)")
    e.add_argument("--no-sampling", action="store_true", help="skip the sampling check of analytic pieces")
    e.add_argument("--sample-out", help="write the azimuthal iron fraction map (npz)")
    e.add_argument("--no-field", action="store_true",
                   help="load the compact without its <fields> and compact/fields includes (no field map needed)")
    c = sub.add_parser("compare", help="compare a model with a reference deck or the Stage 1/2 html data")
    c.add_argument("model", help="model sidecar .json (preferred) or .am")
    c.add_argument("reference", help="reference .am, or the html data .json (with --ref-key)")
    c.add_argument("--model-key")
    c.add_argument("--ref-key", help="key in the html data json: geometry (CLEOv11), v10, v12")
    c.add_argument("--tol", type=float, default=0.5, help="differences thinner than this (cm) are 'small'")
    c.add_argument("--report", help="write the text report here")
    c.add_argument("--plot", help="write an overlay plot (png/svg)")
    c.add_argument("--title")
    c.add_argument("--zoom", type=float, nargs=4, action="append", metavar=("Z0", "Z1", "R0", "R1"))
    k = sub.add_parser("check", help="are the repository files the decks were extracted from unchanged (no DD4hep)?")
    k.add_argument("models", nargs="+", help="deck sidecars (.json)")
    k.add_argument("--repo", required=True, help="solid_dd4hep checkout to compare with")
    f = sub.add_parser("plot-fraction", help="plot the sampled azimuthal iron fraction with the model outline")
    f.add_argument("model", help="model sidecar .json or .am")
    f.add_argument("npz", help="fraction map from extract --sample-out")
    f.add_argument("-o", "--output", required=True)
    f.add_argument("--title")
    args = ap.parse_args(argv)
    if args.cmd == "extract":
        args.detector = args.detector or ["MagnetYoke", "MagnetCoil"]
        args.sampling = (float(args.sampling[0]), float(args.sampling[1]), int(args.sampling[2]))
        return run_extract(args)
    if args.cmd == "compare":
        return run_compare(args)
    if args.cmd == "check":
        return run_check(args)
    return run_plot_fraction(args)


if __name__ == "__main__":
    sys.exit(main())
