# SPDX-License-Identifier: LGPL-3.0-or-later
"""Magnetic material laws for the r-z solver.

SI units inside: B in tesla, H in A/m, reluctivity nu = H/B in m/H.  Every material returns, for an array of
|B|^2 values, the pair (nu, dnu_db_over_b) with dnu_db_over_b = (dnu/d|B|) / |B|, which is what the Newton
Jacobian needs: d(nu B)/dB = nu I + (dnu/d|B|)/|B| * B B^T.

Built-in table "poisson1010": the 1010-steel table that POISSON uses by default (data/bh_1010_poisson.txt).
"""
from __future__ import annotations

import os

import numpy as np

MU0 = 4.0e-7 * np.pi
GAUSS = 1.0e-4                    # tesla per gauss
OERSTED = 1000.0 / (4.0 * np.pi)  # A/m per oersted

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BUILTIN_TABLES = {"poisson1010": "bh_1010_poisson.txt"}

_B_UNITS = {"gauss": GAUSS, "g": GAUSS, "kg": 1000.0 * GAUSS, "kilogauss": 1000.0 * GAUSS, "tesla": 1.0, "t": 1.0}
_H_UNITS = {"oersted": OERSTED, "oe": OERSTED, "a/m": 1.0, "ka/m": 1000.0}

_BTINY = 1.0e-12  # |B| below this (tesla) is treated as zero in derivative terms


class LinearMaterial:
    """Constant relative permeability mu_r (mu_r = 1 is air / coil / aluminium)."""

    nonlinear = False

    def __init__(self, mu_r: float = 1.0, name: str | None = None):
        if not mu_r > 0:
            raise ValueError("mu_r must be positive")
        self.mu_r = float(mu_r)
        self.name = name or ("air" if mu_r == 1.0 else "linear mu_r=%g" % mu_r)

    def nu(self, b2):
        b2 = np.asarray(b2, dtype=float)
        return np.full(b2.shape, 1.0 / (MU0 * self.mu_r)), np.zeros(b2.shape)

    def H(self, b):
        return np.asarray(b, dtype=float) / (MU0 * self.mu_r)

    def energy_density(self, b):
        b = np.asarray(b, dtype=float)
        return 0.5 * b * b / (MU0 * self.mu_r)

    def describe(self) -> str:
        return self.name


class BHCurve:
    """Isotropic nonlinear B-H curve built from (B, H) points.

    interp = "poisson" (default): reluctivity gamma = mu0 H / B (= 1/mu_r) linear in B between the table points,
        as POISSON does (its fast tables are linear interpolations of the (B, gamma) data);
    interp = "pchip": monotone cubic (Fritsch-Carlson) H(B) through the same points (C1, for comparison).
    Below the first non-zero point gamma goes linearly to gamma0 at B = 0 (gamma0 from the table, else
    extrapolated from the first two points as Automesh does).  Beyond the last point the iron is saturated:
    B - mu0 H stays at its last tabulated value (POISSON's rule), so the differential permeability is mu0.
    """

    nonlinear = True

    def __init__(self, B, H, gamma0: float | None = None, name: str = "", interp: str = "poisson"):
        B = np.asarray(B, dtype=float)
        H = np.asarray(H, dtype=float)
        if B.ndim != 1 or B.shape != H.shape:
            raise ValueError("B and H must be 1-D arrays of equal length")
        if B[0] == 0.0:
            if H[0] != 0.0:
                raise ValueError("a B = 0 point must have H = 0")
            B, H = B[1:], H[1:]
        if len(B) < 2 or np.any(np.diff(B) <= 0) or np.any(np.diff(H) <= 0) or B[0] <= 0 or H[0] <= 0:
            raise ValueError("B-H table must have >= 2 points with B and H strictly increasing and positive")
        if np.any(np.diff(H) / np.diff(B) <= 0):
            raise ValueError("B-H table is not monotonic")
        gam = MU0 * H / B
        if gamma0 is None:  # Automesh: extrapolate gamma to B = 0 from the first two entries
            gamma0 = gam[0] - B[0] * (gam[1] - gam[0]) / (B[1] - B[0])
            gamma0 = max(gamma0, 0.2 * gam[0])
        self.name = name or "B-H table"
        self.interp = interp
        self.Bk = np.concatenate([[0.0], B])       # tesla
        self.Hk = np.concatenate([[0.0], H])       # A/m
        self.gk = np.concatenate([[gamma0], gam])  # gamma = mu0 H / B (dimensionless)
        self.Bmax = B[-1]
        self.Msat = B[-1] - MU0 * H[-1]            # tesla, magnetization at the last point
        if interp == "poisson":
            self.sk = np.diff(self.gk) / np.diff(self.Bk)        # dgamma/dB per segment
            # cumulative energy density at the nodes: w = int H dB = int b gamma(b) / mu0 db
            b0, b1 = self.Bk[:-1], self.Bk[1:]
            c0 = self.gk[:-1] - self.sk * b0
            seg = (c0 * (b1 ** 2 - b0 ** 2) / 2.0 + self.sk * (b1 ** 3 - b0 ** 3) / 3.0) / MU0
            self.wk = np.concatenate([[0.0], np.cumsum(seg)])
        elif interp == "pchip":
            from scipy.interpolate import PchipInterpolator

            # H(B) through (0, 0) and the table points (gamma0 is not used by this variant)
            self._pchip = PchipInterpolator(self.Bk, self.Hk, extrapolate=False)
            self._dpchip = self._pchip.derivative()
            self._ipchip = self._pchip.antiderivative()
            self.wk = np.array([0.0, float(self._ipchip(self.Bmax))])
        else:
            raise ValueError("interp must be 'poisson' or 'pchip'")

    # ---------------------------------------------------------------------------------------------------- loading
    @classmethod
    def from_file(cls, path: str, interp: str = "poisson") -> "BHCurve":
        """Read a table: '#' comments with optional 'key: value' keywords (units, gamma0, name, columns)."""
        keys = {"units": "gauss oersted", "columns": "B H"}
        rows = []
        with open(path) as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith("#"):
                    body = s[1:].strip()
                    if ":" in body:
                        k, v = body.split(":", 1)
                        if k.strip().lower() in ("units", "gamma0", "name", "columns"):
                            keys[k.strip().lower()] = v.strip()
                    continue
                s = s.split("#")[0].replace(",", " ")
                vals = [float(x) for x in s.split()]
                if len(vals) < 2:
                    raise ValueError("%s: bad line %r" % (path, line))
                rows.append(vals[:2])
        arr = np.array(rows, dtype=float)
        bu, hu = keys["units"].lower().split()[:2]
        cols = keys["columns"].lower().split()
        if bu not in _B_UNITS:
            raise ValueError("unknown B unit %r" % bu)
        B = arr[:, 0] * _B_UNITS[bu]
        if cols[1] == "h":
            if hu not in _H_UNITS:
                raise ValueError("unknown H unit %r" % hu)
            H = arr[:, 1] * _H_UNITS[hu]
        elif cols[1] in ("mu", "mu_r", "mur"):
            H = np.where(B > 0, B / (MU0 * np.where(arr[:, 1] > 0, arr[:, 1], 1.0)), 0.0)
        elif cols[1] in ("gamma", "nu_r"):
            H = B * arr[:, 1] / MU0
        else:
            raise ValueError("columns must be 'B H', 'B mu' or 'B gamma'")
        gamma0 = float(keys["gamma0"]) if "gamma0" in keys else None
        name = keys.get("name", os.path.basename(path))
        return cls(B, H, gamma0=gamma0, name=name, interp=interp)

    @classmethod
    def builtin(cls, key: str = "poisson1010", interp: str = "poisson") -> "BHCurve":
        return cls.from_file(os.path.join(DATA_DIR, BUILTIN_TABLES[key]), interp=interp)

    # ------------------------------------------------------------------------------------------------ evaluation
    def gamma(self, b):
        """Relative reluctivity gamma = 1/mu_r = mu0 H / B and d gamma / dB (1/T)."""
        b = np.asarray(b, dtype=float)
        g = np.empty_like(b)
        dg = np.empty_like(b)
        hi = b >= self.Bmax
        lo = ~hi
        if self.interp == "poisson":
            bl = b[lo]
            k = np.clip(np.searchsorted(self.Bk, bl, side="right") - 1, 0, len(self.Bk) - 2)
            g[lo] = self.gk[k] + self.sk[k] * (bl - self.Bk[k])
            dg[lo] = self.sk[k]
        else:
            bl = b[lo]
            small = bl < 1e-6 * self.Bk[1]
            bs = np.where(small, 1e-6 * self.Bk[1], bl)
            h = self._pchip(bs)
            dh = self._dpchip(bs)
            g[lo] = MU0 * h / bs
            dg[lo] = np.where(small, 0.0, MU0 * (dh * bs - h) / bs ** 2)
        bh = b[hi]
        g[hi] = 1.0 - self.Msat / bh
        dg[hi] = self.Msat / bh ** 2
        return g, dg

    def nu(self, b2):
        b2 = np.asarray(b2, dtype=float)
        b = np.sqrt(b2)
        g, dg = self.gamma(b)
        nu = g / MU0
        with np.errstate(divide="ignore", invalid="ignore"):
            dnu_over_b = np.where(b > _BTINY, dg / (MU0 * np.maximum(b, _BTINY)), 0.0)
        return nu, dnu_over_b

    def H(self, b):
        b = np.asarray(b, dtype=float)
        g, _ = self.gamma(b)
        return g * b / MU0

    def energy_density(self, b):
        """w(B) = int_0^B H dB (J/m^3)."""
        b = np.asarray(b, dtype=float)
        w = np.empty_like(b)
        hi = b >= self.Bmax
        lo = ~hi
        bl = b[lo]
        if self.interp == "poisson":
            k = np.clip(np.searchsorted(self.Bk, bl, side="right") - 1, 0, len(self.Bk) - 2)
            b0 = self.Bk[k]
            c0 = self.gk[k] - self.sk[k] * b0
            w[lo] = self.wk[k] + (c0 * (bl ** 2 - b0 ** 2) / 2.0 + self.sk[k] * (bl ** 3 - b0 ** 3) / 3.0) / MU0
        else:
            w[lo] = self._ipchip(bl)
        bh = b[hi]
        w[hi] = self.wk[-1] + ((bh ** 2 - self.Bmax ** 2) / 2.0 - self.Msat * (bh - self.Bmax)) / MU0
        return w

    def describe(self) -> str:
        return "%s [%d points to %.4g T, %s interpolation, saturated above]" % (
            self.name, len(self.Bk), self.Bmax, self.interp)


def material_from_spec(spec: str, interp: str = "poisson"):
    """'poisson1010' (built-in), 'mu=<value>' / 'linear:<mu_r>' (fixed), 'gamma=<value>' or a table file path."""
    s = spec.strip()
    low = s.lower()
    if low in BUILTIN_TABLES:
        return BHCurve.builtin(low, interp=interp)
    for prefix in ("mu=", "linear:", "mu_r="):
        if low.startswith(prefix):
            return LinearMaterial(float(s[len(prefix):]))
    if low.startswith("gamma="):
        return LinearMaterial(1.0 / float(s[len("gamma="):]), name="fixed gamma=%s" % s[len("gamma="):])
    if os.path.isfile(s):
        return BHCurve.from_file(s, interp=interp)
    raise ValueError("unknown material spec %r (use poisson1010, mu=<mu_r>, gamma=<1/mu_r> or a table file)" % spec)
