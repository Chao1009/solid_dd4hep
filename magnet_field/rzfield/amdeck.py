# SPDX-License-Identifier: LGPL-3.0-or-later
"""Reader for the subset of the Poisson Superfish Automesh input (.am) used by magnet decks.

Supported: title lines before the first namelist; comments after ';'; '&reg ... &' namelists (the first one, with
KPROB = 0, holds the problem constants and opens region 1 = the problem boundary); '&po x=.., y=.. &' points
(straight segments, NT = 1); '&mt ... &' material tables (MTID, BH / BMU / BGAM pairs, GAMMA / MU).
Cylindrical problems only (ICYLIN = 1): x is r and y is z.  Lengths are converted to cm with CONV.
Region semantics follow Automesh: a later region overrides earlier ones where they overlap.  CUR (or CURRENT) is
the total current of the region in ampere-turns, multiplied by XJFACT; DEN (or DENSITY) is a current density in
A per (length unit)^2.  Mesh parameters (DX, XREG, KREG, ...) and solver controls (RHOGAM, ...) are kept in
``Deck.params`` but not used.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from .bh import BHCurve, LinearMaterial, MU0, GAUSS, OERSTED

_NAMELIST = re.compile(r"&\s*([A-Za-z]+)(.*?)&", re.S)
_KEY = re.compile(r"([A-Za-z][A-Za-z0-9_]*(?:\(\s*\d+\s*\))?)\s*=")
_REGION_KEYS = {"mat", "cur", "current", "den", "density", "mtid", "mshape", "stacking", "ibound", "itri",
                "charge", "voltage", "npoint", "nreg"}
_PO_KEYS = {"x", "y", "nt"}


def _parse_value(tok: str):
    t = tok.strip()
    try:
        return float(t.lower().replace("d", "e"))
    except ValueError:
        return t.strip("'\"")


def parse_namelist_body(body: str) -> dict:
    """'key=value, key=v1,v2,...' -> {key: value or [values]} with lower-case keys."""
    out = {}
    matches = list(_KEY.finditer(body))
    if matches and body[:matches[0].start()].strip(" \t\r\n,"):
        raise ValueError("cannot parse namelist text %r" % body[:matches[0].start()].strip())
    for i, m in enumerate(matches):
        key = re.sub(r"\s+", "", m.group(1).lower())
        end = matches[i + 1].start() if i + 1 < len(matches) else len(body)
        toks = [t for t in re.split(r"[,\s]+", body[m.end():end].strip()) if t]
        vals = [_parse_value(t) for t in toks]
        out[key] = vals[0] if len(vals) == 1 else vals
    return out


def polygon_area(pts: np.ndarray) -> float:
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


@dataclass
class Region:
    index: int                 # Poisson region number (1 = problem boundary)
    mat: int = 1
    cur: float = 0.0           # total current, ampere-turns (XJFACT applied)
    mtid: int | None = None
    mshape: int | None = None
    name: str = ""
    params: dict = field(default_factory=dict)
    points: np.ndarray = None  # (n, 2) polygon in cm, columns (r, z), not closed

    @property
    def area_cm2(self) -> float:
        return abs(polygon_area(self.points))


@dataclass
class Deck:
    title: str
    params: dict
    regions: list
    mt: dict                   # MTID -> MT namelist dict
    path: str = ""

    # ------------------------------------------------------------------------------------------------ geometry
    @property
    def boundary(self) -> Region:
        return self.regions[0]

    def bbox(self, which: str = "all"):
        """(rmin, rmax, zmin, zmax) in cm of region 1 ('box'), of iron + coils ('magnet') or of everything."""
        if which == "box":
            regs = [self.regions[0]]
        elif which == "magnet":
            regs = [r for r in self.regions[1:] if r.mat >= 2 or r.cur != 0.0]
        else:
            regs = self.regions
        p = np.vstack([r.points for r in regs])
        return float(p[:, 0].min()), float(p[:, 0].max()), float(p[:, 1].min()), float(p[:, 1].max())

    # ----------------------------------------------------------------------------------------------- materials
    @property
    def mode(self) -> int:
        xj = float(self.params.get("xjfact", 1.0))
        return int(self.params.get("mode", -2 if xj != 0 else -1))

    def material_settings(self) -> dict:
        """mat -> (mshape, mtid): the last REG setting for each material number wins (Automesh rule)."""
        out = {}
        for r in self.regions[1:]:
            if r.mat < 2:
                continue
            ms, mt = out.get(r.mat, (None, None))
            if r.mshape is not None:
                ms = r.mshape
            if r.mtid is not None:
                mt = r.mtid
            out[r.mat] = (ms, mt)
        return out

    def material(self, mat: int, interp: str = "poisson"):
        """Material law for material number mat as Poisson would use it."""
        if mat <= 1:
            return LinearMaterial(1.0, name="air")
        mode = self.mode
        mshape, mtid = self.material_settings().get(mat, (None, None))
        if mshape is None:
            mshape = mode
        fixgam = float(self.params.get("fixgam", 0.004))
        if mode == -2:
            return LinearMaterial(1.0e5, name="infinite permeability (MODE=-2) as mu_r=1e5")
        mtab = self.mt.get(mtid) if (mtid is not None and mtid > 0) else None
        if mode == -1 or mshape == -1:
            if mtab is not None and "mu" in mtab:
                return LinearMaterial(float(mtab["mu"]), name="MT %d mu=%g" % (mtid, float(mtab["mu"])))
            gam = float(mtab["gamma"]) if (mtab is not None and "gamma" in mtab) else fixgam
            return LinearMaterial(1.0 / gam, name="fixed gamma=%g (mu_r=%g)" % (gam, 1.0 / gam))
        if mshape != 0:
            raise ValueError("material %d: MSHAPE=%s (anisotropic / permanent magnet) is not supported" % (mat, mshape))
        if mtid is None or mtid == 0:
            return BHCurve.builtin("poisson1010", interp=interp)
        if mtid < 0:
            raise ValueError("material %d: internal table MTID=%d is not available (only MTID=0, 1010 steel)"
                             % (mat, mtid))
        if mtab is None:
            raise ValueError("material %d: no &mt namelist with MTID=%d" % (mat, mtid))
        return _mt_table(mtab, mtid, interp)

    def summary(self) -> str:
        lines = ["%s: %s" % (self.path or "deck", self.title)]
        b = self.bbox("box")
        lines.append("  problem box r %.1f..%.1f z %.1f..%.1f cm, MODE=%d, XJFACT=%g, %d regions"
                     % (b[0], b[1], b[2], b[3], self.mode, float(self.params.get("xjfact", 1.0)), len(self.regions)))
        tot = 0.0
        for r in self.regions[1:]:
            tot += r.cur
            lines.append("  region %2d mat %d cur %12.1f A  area %9.3f cm2  %s"
                         % (r.index, r.mat, r.cur, r.area_cm2, r.name))
        lines.append("  total current %.6g A-turns" % tot)
        return "\n".join(lines)


def _mt_table(mtab: dict, mtid: int, interp: str):
    for key, kind in (("bh", "h"), ("bmu", "mu"), ("bgam", "gamma")):
        if key in mtab:
            v = np.asarray(mtab[key], dtype=float).reshape(-1, 2)
            B = v[:, 0] * GAUSS
            if kind == "h":
                H = v[:, 1] * OERSTED
            elif kind == "mu":
                H = np.where(B > 0, B / (MU0 * v[:, 1]), 0.0)
            else:
                H = B * v[:, 1] / MU0
            return BHCurve(B, H, name="MT %d (%s table)" % (mtid, key.upper()), interp=interp)
    if "mu" in mtab:
        return LinearMaterial(float(mtab["mu"]), name="MT %d mu=%g" % (mtid, float(mtab["mu"])))
    if "gamma" in mtab:
        return LinearMaterial(1.0 / float(mtab["gamma"]), name="MT %d gamma=%g" % (mtid, float(mtab["gamma"])))
    raise ValueError("&mt mtid=%d has no BH/BMU/BGAM table and no MU/GAMMA" % mtid)


def read_deck(path: str) -> Deck:
    with open(path, encoding="latin-1") as f:
        raw_lines = f.read().splitlines()
    # title = lines before the first namelist; the comment just above a &reg names the region
    title_lines, body_lines, comments = [], [], {}
    started = False
    last_comment = ""
    for ln in raw_lines:
        code, _, com = ln.partition(";")
        if not started:
            if "&" in code:
                started = True
            else:
                if code.strip():
                    title_lines.append(code.strip())
                continue
        if com.strip() and not code.strip():
            last_comment = com.strip()
        if re.search(r"&\s*reg\b", code, re.I):
            comments[len(body_lines)] = last_comment
            last_comment = ""
        body_lines.append(code)
    text = "\n".join(body_lines)
    # map character offsets back to the line index for region comments
    line_starts = np.cumsum([0] + [len(s) + 1 for s in body_lines])

    params, regions, mt = None, [], {}
    cur_region = None
    for m in _NAMELIST.finditer(text):
        group = m.group(1).lower()
        vals = parse_namelist_body(m.group(2))
        if group == "reg":
            line_idx = int(np.searchsorted(line_starts, m.start(), side="right") - 1)
            name = comments.get(line_idx, "")
            if params is None:
                if "kprob" not in vals:
                    raise ValueError("%s: the first &reg must declare KPROB" % path)
                params = vals
                if int(params.get("kprob", 0)) != 0:
                    raise ValueError("%s: only Poisson problems (KPROB = 0) are supported" % path)
                if int(params.get("icylin", 0)) != 1:
                    raise ValueError("%s: only cylindrical problems (ICYLIN = 1) are supported" % path)
                if float(params.get("xjfact", 1.0)) == 0.0:
                    raise ValueError("%s: XJFACT = 0 (electrostatics) is not supported" % path)
                cur_region = Region(index=1, mat=1, name=name or "problem boundary", params={})
            else:
                unknown = set(vals) - _REGION_KEYS
                if unknown:
                    raise ValueError("%s: unsupported REG keys %s" % (path, sorted(unknown)))
                cur_region = Region(index=len(regions) + 1, params=vals, name=name,
                                    mat=int(vals.get("mat", 1)),
                                    mtid=int(vals["mtid"]) if "mtid" in vals else None,
                                    mshape=int(vals["mshape"]) if "mshape" in vals else None)
            cur_region.points = []
            regions.append(cur_region)
        elif group == "po":
            if cur_region is None:
                raise ValueError("%s: &po before the first &reg" % path)
            unknown = set(vals) - _PO_KEYS
            if unknown or int(vals.get("nt", 1)) != 1:
                raise ValueError("%s: region %d: only straight-line points '&po x=.., y=.. &' are supported (got %s)"
                                 % (path, cur_region.index, sorted(vals)))
            if "x" not in vals or "y" not in vals:
                raise ValueError("%s: region %d: &po needs both x and y" % (path, cur_region.index))
            cur_region.points.append((float(vals["x"]), float(vals["y"])))
        elif group == "mt":
            mtid = int(vals.get("mtid", len(mt) + 1))
            mt[mtid] = vals
        else:
            raise ValueError("%s: unsupported namelist &%s" % (path, group))
    if params is None or not regions:
        raise ValueError("%s: no regions found" % path)

    conv = float(params.get("conv", 1.0))
    xjfact = float(params.get("xjfact", 1.0))
    for r in regions:
        pts = np.asarray(r.points, dtype=float) * conv
        if len(pts) > 1 and np.allclose(pts[0], pts[-1]):
            pts = pts[:-1]
        keep = np.ones(len(pts), bool)          # drop consecutive duplicates
        keep[1:] = np.any(np.abs(np.diff(pts, axis=0)) > 1e-9, axis=1)
        pts = pts[keep]
        if len(pts) < 3 or abs(polygon_area(pts)) < 1e-12:
            raise ValueError("%s: region %d (%s) is not a closed polygon with area" % (path, r.index, r.name))
        r.points = pts
        if r.index > 1:
            cur = r.params.get("cur", r.params.get("current"))
            den = r.params.get("den", r.params.get("density"))
            if cur is not None and den is not None:
                raise ValueError("%s: region %d has both CUR and DEN" % (path, r.index))
            if den is not None:
                cur = float(den) * abs(polygon_area(np.asarray(r.points) / conv))
            r.cur = float(cur or 0.0) * xjfact
            if r.cur != 0.0 and r.mat != 1:
                raise ValueError("%s: region %d carries current but has MAT=%d" % (path, r.index, r.mat))
    if np.min(regions[0].points[:, 0]) < -1e-9:
        raise ValueError("%s: the problem boundary extends to r < 0" % path)
    title = " / ".join(title_lines)
    return Deck(title=title, params=params, regions=regions, mt=mt, path=path)
