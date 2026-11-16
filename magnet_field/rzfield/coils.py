# SPDX-License-Identifier: LGPL-3.0-or-later
"""Field of air-core coils: coaxial circular coils of rectangular cross-section with a uniform current density.

No iron and no mesh: every coil cross-section is replaced by circular current loops at the nodes of a
Gauss-Legendre rule, and the field of a loop is known in closed form (complete elliptic integrals). The order of
the rule grows as the field point comes closer to the coil, so that the far field costs little. Next to a winding
and inside it the field of the nearest loops diverges like 1/distance: there the cross-section is split into four
rectangles that meet at the point of the winding nearest to the field point, and the nodes are packed towards
that corner. tests/test_coils.py checks the result against the closed form on the axis, against Maxwell's
equations and against a rule of much higher order.

Used for the polarized NH3 target magnet (a split pair of superconducting coils on the hall x axis), whose field
does not depend on the iron of the SoLID solenoid in this approximation: see README.md, "Target magnet".

Units: lengths in cm, currents in ampere-turns, fields in gauss.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
from scipy.special import ellipe, ellipk

MU0_G_CM_PER_A = 4.0e-7 * np.pi * 1.0e4 * 1.0e2      # mu0 in gauss cm / A (T m / A = 1e4 G * 1e2 cm / A)

# Gauss-Legendre order by distance from the coil cross-section, in units of its larger side:
# (distance >= this, order). Closer than the last entry the corner rule below is used.
ORDERS = ((8.0, 3), (3.0, 5), (1.0, 8), (0.3, 16))
NEAR_ORDER = 24                  # per direction and per quarter, with nodes packed towards the common corner
NEAR_POWER = 3                   # node packing: x = corner + (edge - corner) * u**NEAR_POWER


@dataclass
class Coil:
    name: str
    r_in: float            # cm
    r_out: float
    z_min: float
    z_max: float
    ampere_turns: float    # signed: positive = field along +z at the coil centre

    @property
    def current_density(self) -> float:
        """A / cm2"""
        return self.ampere_turns / ((self.r_out - self.r_in) * (self.z_max - self.z_min))


@dataclass
class CoilSet:
    name: str
    coils: list
    axis: str = "z"                  # hall axis of the magnet axis (for the map header)
    source: str = ""

    @property
    def ampere_turns(self) -> float:
        return float(sum(abs(c.ampere_turns) for c in self.coils))


def read_coils(path: str) -> CoilSet:
    """Coil set from a JSON file (models/*.json with "coils": r_in_cm, r_out_cm, z_min_cm, z_max_cm, turns)."""
    with open(path) as f:
        d = json.load(f)
    current = float(d["current_A"])
    coils = []
    for c in d["coils"]:
        coils.append(Coil(str(c.get("name", "coil%d" % (len(coils) + 1))), float(c["r_in_cm"]), float(c["r_out_cm"]),
                          float(c["z_min_cm"]), float(c["z_max_cm"]), float(c["turns"]) * current))
        k = coils[-1]
        if not (0.0 <= k.r_in < k.r_out and k.z_min < k.z_max):
            raise ValueError("%s: coil %s needs 0 <= r_in < r_out and z_min < z_max" % (path, k.name))
    if not coils:
        raise ValueError("%s: no coils" % path)
    return CoilSet(str(d.get("name", "coils")), coils, str(d.get("axis", "z")), str(d.get("source", "")))


def loop_field(a, z0, rho, z):
    """Field of circular current loops, per ampere.

    a, z0: radius and axial position of the loops, arrays (M,), or (N, M) for loops of their own per point;
    rho, z: field points, arrays (N,); all in cm. Returns (B_rho, B_z), arrays (N, M), in gauss per ampere."""
    a, z0 = np.asarray(a, float), np.asarray(z0, float)
    if a.ndim == 1:
        a, z0 = a[None, :], z0[None, :]
    rho = np.abs(np.asarray(rho, float))[:, None]
    zeta = np.asarray(z, float)[:, None] - z0
    d2 = (a + rho) ** 2 + zeta ** 2
    e2 = np.maximum((a - rho) ** 2 + zeta ** 2, 1e-300)
    m = np.minimum(4.0 * a * rho / d2, 1.0 - 1e-15)            # k^2
    K, E = ellipk(m), ellipe(m)
    c = MU0_G_CM_PER_A / (2.0 * np.pi) / np.sqrt(d2)
    s2 = rho * rho + zeta * zeta
    bz = c * (K + (a * a - s2) / e2 * E)
    with np.errstate(divide="ignore", invalid="ignore"):
        br = np.where(rho > 1e-12, c * zeta / rho * ((a * a + s2) / e2 * E - K), 0.0)
    return br, bz


def _distance(coil: Coil, rho, z):
    """distance (cm) from the points to the coil cross-section, 0 inside"""
    dr = np.maximum(np.maximum(coil.r_in - rho, rho - coil.r_out), 0.0)
    dz = np.maximum(np.maximum(coil.z_min - z, z - coil.z_max), 0.0)
    return np.hypot(dr, dz)


def _near_field(coil: Coil, rho, z, order: int, power: int):
    """(B_rho, B_z) of one coil at points next to or inside its winding: four rectangles meeting at the point of
    the cross-section nearest to the field point, Gauss-Legendre in u with x = corner + (edge - corner) u**power"""
    x, w = np.polynomial.legendre.leggauss(order)
    u, wu = 0.5 * (x + 1.0), 0.5 * w                           # rule on (0, 1)
    t, wt = u ** power, wu * power * u ** (power - 1)           # fraction of the way from the corner to the edge
    rc = np.clip(rho, coil.r_in, coil.r_out)[:, None, None]
    zc = np.clip(z, coil.z_min, coil.z_max)[:, None, None]
    jd = coil.current_density
    br = np.zeros(len(rho))
    bz = np.zeros(len(rho))
    for redge in (coil.r_in, coil.r_out):
        for zedge in (coil.z_min, coil.z_max):
            dr, dz = redge - rc, zedge - zc                    # (N, 1, 1), signed
            a = np.broadcast_to(rc + dr * t[None, :, None], (len(rho), order, order))
            z0 = np.broadcast_to(zc + dz * t[None, None, :], (len(rho), order, order))
            cur = jd * np.abs(dr) * np.abs(dz) * wt[None, :, None] * wt[None, None, :]    # ampere per loop
            lr, lz = loop_field(a.reshape(len(rho), -1), z0.reshape(len(rho), -1), rho, z)
            cur = cur.reshape(len(rho), -1)
            br += (lr * cur).sum(1)
            bz += (lz * cur).sum(1)
    return br, bz


def coil_field(coils, rho, z, orders=ORDERS, near_order: int = NEAR_ORDER, near_power: int = NEAR_POWER,
               max_terms: int = 4_000_000):
    """Field (B_rho, B_z) in gauss of the coils at the points (rho, z) in cm; arrays of any common shape.

    `coils` is a CoilSet or a list of Coil. `orders` gives the Gauss-Legendre order by distance; closer than its
    last entry the corner rule of order `near_order` is used. `max_terms` bounds the size of the (points x loops)
    arrays worked on at once."""
    coils = coils.coils if isinstance(coils, CoilSet) else coils
    shape = np.broadcast(rho, z).shape
    rho = np.abs(np.broadcast_to(np.asarray(rho, float), shape)).ravel()
    z = np.broadcast_to(np.asarray(z, float), shape).ravel()
    br = np.zeros(len(rho))
    bz = np.zeros(len(rho))
    for coil in coils:
        size = max(coil.r_out - coil.r_in, coil.z_max - coil.z_min)
        d = _distance(coil, rho, z) / size
        todo = np.ones(len(rho), bool)
        for dmin, order in orders:
            sel = np.nonzero(todo & (d >= dmin))[0]
            todo[sel] = False
            if len(sel) == 0:
                continue
            x, w = np.polynomial.legendre.leggauss(order)
            ar = 0.5 * (coil.r_in + coil.r_out) + 0.5 * (coil.r_out - coil.r_in) * x
            az = 0.5 * (coil.z_min + coil.z_max) + 0.5 * (coil.z_max - coil.z_min) * x
            A, Z0 = np.meshgrid(ar, az, indexing="ij")
            cur = coil.ampere_turns * 0.25 * (w[:, None] * w[None, :]).ravel()    # ampere per loop
            A, Z0 = A.ravel(), Z0.ravel()
            step = max(1, max_terms // len(A))
            for i in range(0, len(sel), step):
                s = sel[i:i + step]
                lr, lz = loop_field(A, Z0, rho[s], z[s])
                br[s] += lr @ cur
                bz[s] += lz @ cur
        sel = np.nonzero(todo)[0]                              # next to the winding or inside it
        step = max(1, max_terms // (near_order * near_order))
        for i in range(0, len(sel), step):
            s = sel[i:i + step]
            lr, lz = _near_field(coil, rho[s], z[s], near_order, near_power)
            br[s] += lr
            bz[s] += lz
    return br.reshape(shape), bz.reshape(shape)


def axis_field(coils, z):
    """B_z (gauss) on the axis, from the closed form for a coil of rectangular cross-section. For tests."""
    coils = coils.coils if isinstance(coils, CoilSet) else coils
    z = np.asarray(z, float)
    out = np.zeros(z.shape)
    for c in coils:
        def part(u):
            return u * np.log((c.r_out + np.sqrt(c.r_out ** 2 + u * u)) / (c.r_in + np.sqrt(c.r_in ** 2 + u * u)))
        out += 0.5 * MU0_G_CM_PER_A * c.current_density * (part(z - c.z_min) - part(z - c.z_max))
    return out


def in_winding(coils, rho, z):
    """points inside a coil cross-section"""
    coils = coils.coils if isinstance(coils, CoilSet) else coils
    rho, z = np.abs(np.asarray(rho, float)), np.asarray(z, float)
    out = np.zeros(np.broadcast(rho, z).shape, bool)
    for c in coils:
        out |= (rho >= c.r_in) & (rho <= c.r_out) & (z >= c.z_min) & (z <= c.z_max)
    return out
