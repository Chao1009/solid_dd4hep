# SPDX-License-Identifier: LGPL-3.0-or-later
"""Small polygon helpers on deck geometry (cm): point-in-polygon, region lookup, distance to edges."""
from __future__ import annotations

import numpy as np


def seg_point_dist(p, a, b):
    """Distances from points p (n, 2) to segments a-b (m, 2): array (n, m)."""
    ab = b - a
    L2 = np.maximum(np.einsum("mi,mi->m", ab, ab), 1e-300)
    ap = p[:, None, :] - a[None, :, :]
    t = np.clip(np.einsum("nmi,mi->nm", ap, ab) / L2, 0.0, 1.0)
    q = a[None, :, :] + t[..., None] * ab[None, :, :]
    return np.linalg.norm(p[:, None, :] - q, axis=2)


def points_in_polygon(x, y, poly):
    """Even-odd rule; x, y arrays of equal shape, poly (n, 2). Points on edges are ambiguous."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    inside = np.zeros(x.shape, dtype=bool)
    px, py = poly[:, 0], poly[:, 1]
    qx, qy = np.roll(px, -1), np.roll(py, -1)
    for x1, y1, x2, y2 in zip(px, py, qx, qy):
        if y1 == y2:
            continue
        cond = (y1 > y) != (y2 > y)
        xc = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
        inside ^= cond & (x < xc)
    return inside


def region_at(deck, r_cm, z_cm):
    """Index of the last deck region containing each point (Automesh rule), 0 outside region 1."""
    r_cm = np.asarray(r_cm, dtype=float)
    z_cm = np.asarray(z_cm, dtype=float)
    out = np.zeros(r_cm.shape, dtype=np.int32)
    for reg in deck.regions:
        out[points_in_polygon(r_cm, z_cm, reg.points)] = reg.index
    return out


def distance_to_edges(deck, r_cm, z_cm, which: str = "iron", chunk: int = 20000):
    """Distance (cm) from each point to the nearest edge of the iron (MAT >= 2) or of all non-box regions."""
    regs = [rg for rg in deck.regions[1:] if (rg.mat >= 2 or which != "iron")]
    shape = np.shape(r_cm)
    pts = np.column_stack([np.ravel(r_cm), np.ravel(z_cm)]).astype(float)
    if not regs:
        return np.full(shape, np.inf)
    a = np.vstack([rg.points for rg in regs])
    b = np.vstack([np.roll(rg.points, -1, axis=0) for rg in regs])
    out = np.empty(len(pts))
    for s in range(0, len(pts), chunk):
        out[s:s + chunk] = seg_point_dist(pts[s:s + chunk], a, b).min(axis=1)
    return out.reshape(shape)


def material_masks(deck, r_cm, z_cm, margin_cm: float = 1.0):
    """(iron, iron_far, air_far): points in MAT >= 2 regions (after overrides); iron points at least margin_cm
    from every iron edge; air points (MAT 1 without current) at least margin_cm from every region edge."""
    reg = region_at(deck, r_cm, z_cm)
    mats = np.array([0] + [rg.mat for rg in deck.regions])
    curs = np.array([0.0] + [rg.cur for rg in deck.regions])
    iron = mats[reg] >= 2
    air = (mats[reg] == 1) & (curs[reg] == 0.0)
    iron_far = iron & (distance_to_edges(deck, r_cm, z_cm, "iron") >= margin_cm)
    air_far = air & (distance_to_edges(deck, r_cm, z_cm, "all") >= margin_cm)
    return iron, iron_far, air_far
