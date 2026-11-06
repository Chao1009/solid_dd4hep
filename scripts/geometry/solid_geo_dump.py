#!/usr/bin/env python3
"""Dump hall-frame extents (z, r, phi) of a DD4hep geometry, per top-level detector and per component.

Run inside eic-shell / the eic_xl container with the plugin library on LD_LIBRARY_PATH:
    python3 solid_geo_dump.py <compact.xml> [--depth N] [--json out.json] [--md out.md] [--grep REGEX]

Conventions (chosen to match solid_gemc/script/geometry_outline.py so tables are comparable):
  * units: cm, deg; frame: the DD4hep world frame (= GEMC hall frame for SoLID: origin at the coil center)
  * extents come from the mesh vertices of every non-assembly volume (TGeoShape::SetPoints), transformed to the
    world frame. For polyhedra this means rmax is the OUTER CORNER (as geometry_outline.py), and rmin is the inner
    corner; the inner flat (apothem) is reported separately as 'rmin_flat' when a TGeoPgon is placed coaxially.
  * assemblies are traversed (they have no shape of their own); a non-assembly volume is a "leaf" for the extent
    unless --depth allows descending into its daughters as separate components.
"""
import argparse
import json
import math
import re
import sys

import numpy as np


def mat_to_np(m):
    r = np.array([m.GetRotationMatrix()[i] for i in range(9)], dtype=float).reshape(3, 3)
    t = np.array([m.GetTranslation()[i] for i in range(3)], dtype=float)
    return r, t


def shape_points(shape):
    n = shape.GetNmeshVertices()
    if n <= 0:
        return None
    buf = np.zeros(3 * n, dtype=np.float64)
    shape.SetPoints(buf)
    return buf.reshape(n, 3)


class Extent:
    def __init__(self):
        self.zmin = self.rmin = math.inf
        self.zmax = self.rmax = -math.inf
        self.rmin_flat = math.inf
        self.nvol = 0
        self.materials = {}

    def add_points(self, pts):
        r = np.hypot(pts[:, 0], pts[:, 1])
        self.zmin = min(self.zmin, float(pts[:, 2].min()))
        self.zmax = max(self.zmax, float(pts[:, 2].max()))
        self.rmin = min(self.rmin, float(r.min()))
        self.rmax = max(self.rmax, float(r.max()))

    def merge(self, o):
        self.zmin, self.zmax = min(self.zmin, o.zmin), max(self.zmax, o.zmax)
        self.rmin, self.rmax = min(self.rmin, o.rmin), max(self.rmax, o.rmax)
        self.rmin_flat = min(self.rmin_flat, o.rmin_flat)
        self.nvol += o.nvol
        for k, v in o.materials.items():
            self.materials[k] = self.materials.get(k, 0) + v

    def as_dict(self):
        f = lambda x: None if not math.isfinite(x) else round(x, 4)
        return dict(zmin=f(self.zmin), zmax=f(self.zmax), rmin=f(self.rmin), rmax=f(self.rmax),
                    rmin_flat=f(self.rmin_flat), nvol=self.nvol, materials=self.materials)


BIG = 1.0e5  # cm, larger than the world


def raycast_extent(shape, R, t, n_phi=72, n_z=24):
    """Exact-surface extents of any (e.g. boolean) shape by ray casting in the world frame.

    r_min: rays from the world z axis outward; r_max: rays from far outside inward; z range: rays along +-z.
    Returns (zmin, zmax, rmin, rmax) or None if no ray hits the shape."""
    import ROOT
    Rt = R.T
    # world-frame bounding box of the shape from its local bbox corners
    bb = shape  # TGeoBBox interface is available on every TGeoShape
    o = np.array([bb.GetOrigin()[i] for i in range(3)])
    d = np.array([bb.GetDX(), bb.GetDY(), bb.GetDZ()])
    corners = np.array([o + d * np.array([sx, sy, sz]) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    wc = corners @ Rt + t
    zlo, zhi = wc[:, 2].min(), wc[:, 2].max()
    rhi = np.hypot(wc[:, 0], wc[:, 1]).max() + 1.0
    buf_p, buf_d = ROOT.std.vector('double')(3), ROOT.std.vector('double')(3)

    def dist(pw, dw):
        pl = Rt @ (np.asarray(pw) - t)
        dl = Rt @ np.asarray(dw)
        for i in range(3):
            buf_p[i] = pl[i]
            buf_d[i] = dl[i]
        return shape.DistFromOutside(buf_p.data(), buf_d.data())

    rmin, rmax, zmin, zmax = math.inf, -math.inf, math.inf, -math.inf
    for iz in range(n_z):
        z = zlo + (zhi - zlo) * (iz + 0.5) / n_z
        for ip in range(n_phi):
            ph = 2 * math.pi * ip / n_phi
            u = np.array([math.cos(ph), math.sin(ph), 0.0])
            s = dist([0.0, 0.0, z], u)
            if s < BIG:
                rmin = min(rmin, s)
            s = dist(u * rhi + [0, 0, z], -u)
            if s < BIG:
                rmax = max(rmax, rhi - s)
    for ir in range(n_z):
        r = rhi * (ir + 0.5) / n_z
        for ip in range(n_phi // 2):
            ph = 2 * math.pi * ip / (n_phi // 2)
            x, y = r * math.cos(ph), r * math.sin(ph)
            s = dist([x, y, zlo - 1.0], [0, 0, 1])
            if s < BIG:
                zmin = min(zmin, zlo - 1.0 + s)
            s = dist([x, y, zhi + 1.0], [0, 0, -1])
            if s < BIG:
                zmax = max(zmax, zhi + 1.0 - s)
    if not math.isfinite(rmin) and not math.isfinite(zmin):
        return None
    # a ray starting on the axis inside the shape returns 0 -> the shape reaches the axis
    return zmin, zmax, (0.0 if not math.isfinite(rmin) else rmin), rmax


def leaf_extent(node, R, t):
    vol = node.GetVolume()
    shape = vol.GetShape()
    e = Extent()
    if shape.ClassName() == 'TGeoCompositeShape':
        # mesh points of boolean shapes do not describe subtracted bores/cuts: ray-cast instead
        rc = raycast_extent(shape, R, t)
        if rc is not None:
            e.zmin, e.zmax, e.rmin, e.rmax = rc
        pts = None
    else:
        pts = shape_points(shape)
    if pts is not None:
        e.add_points(pts @ R.T + t)
    # inner flat of a coaxial polyhedra (TGeoPgon) placed with only a z-rotation
    if shape.ClassName() == 'TGeoPgon' and abs(abs(R[2, 2]) - 1) < 1e-9 and np.hypot(t[0], t[1]) < 1e-6:
        e.rmin_flat = min(shape.GetRmin()[i] for i in range(shape.GetNz()))
    e.nvol = 1
    mat = vol.GetMaterial().GetName() if vol.GetMaterial() else '?'
    e.materials[mat] = 1
    return e


def walk(node, R, t, depth_left, out_components, path, grep):
    """Return the Extent of this node (recursing through assemblies)."""
    vol = node.GetVolume()
    is_assembly = vol.IsAssembly()
    total = Extent()
    if not is_assembly:
        total.merge(leaf_extent(node, R, t))
        if grep is None or re.search(grep, path):
            out_components.append((path, node.GetVolume().GetName(), total))
        if depth_left <= 0:
            return total
        depth_left -= 1
    for i in range(vol.GetNdaughters()):
        d = vol.GetNode(i)
        r, tt = mat_to_np(d.GetMatrix())
        sub = walk(d, R @ r, R @ tt + t, depth_left, out_components, f"{path}/{d.GetName()}", grep)
        if is_assembly:
            total.merge(sub)
    return total


def dump_detectors(desc, depth=0, grep=None):
    """Extents per top-level DetElement (named as in the compact file).

    World-volume daughters that do not belong to any top-level DetElement are reported under their node name."""
    world_vol = desc.manager().GetTopVolume()
    results, seen = [], set()
    for name, de in sorted(desc.world().children(), key=lambda kv: kv[1].id()):
        pv = de.placement()
        node = pv.ptr() if hasattr(pv, 'ptr') else pv
        if not node:
            continue
        seen.add(node.GetName())
        R, t = mat_to_np(node.GetMatrix())
        comps = []
        ext = walk(node, R, t, depth, comps, str(name), grep)
        results.append(dict(detector=str(name), id=int(de.id()), extent=ext.as_dict(),
                            components=[dict(path=p, volume=v, extent=e.as_dict()) for p, v, e in comps[:5000]]))
    for i in range(world_vol.GetNdaughters()):
        node = world_vol.GetNode(i)
        if node.GetName() in seen:
            continue
        R, t = mat_to_np(node.GetMatrix())
        comps = []
        ext = walk(node, R, t, depth, comps, node.GetName(), grep)
        results.append(dict(detector=node.GetName(), id=-1, extent=ext.as_dict(),
                            components=[dict(path=p, volume=v, extent=e.as_dict()) for p, v, e in comps[:5000]]))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('compact')
    ap.add_argument('--depth', type=int, default=0, help='levels to descend below the first non-assembly volume')
    ap.add_argument('--grep', default=None, help='only list components whose path matches REGEX')
    ap.add_argument('--json', default=None)
    ap.add_argument('--md', default=None)
    a = ap.parse_args()

    import ROOT  # noqa: F401
    import dd4hep
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(a.compact)
    mgr = desc.manager()
    world = mgr.GetTopVolume()

    results = dump_detectors(desc, a.depth, a.grep)

    # constants are handy for the comparison table
    consts = {}
    try:
        for name, c in desc.constants():
            try:
                consts[str(name)] = float(desc.constantAsDouble(str(name)))
            except Exception:
                pass
    except Exception:
        pass

    lines = ['| top-level volume | z from | z to | r from | r to | r flat | #vol | materials |',
             '|---|---:|---:|---:|---:|---:|---:|---|']
    for r in results:
        e = r['extent']
        mats = ', '.join(f"{k}x{v}" for k, v in sorted(e['materials'].items(), key=lambda kv: -kv[1])[:4])
        lines.append(f"| {r['detector']} | {e['zmin']} | {e['zmax']} | {e['rmin']} | {e['rmax']} | "
                     f"{e['rmin_flat'] if e['rmin_flat'] is not None else ''} | {e['nvol']} | {mats} |")
    md = '\n'.join(lines)
    print(md)
    if a.md:
        open(a.md, 'w').write(md + '\n')
    if a.json:
        json.dump(dict(compact=a.compact, detectors=results, constants=consts), open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
