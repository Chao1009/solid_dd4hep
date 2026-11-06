#!/usr/bin/env python3
"""Parser for GEMC TEXT-factory geometry files (``<name>__geometry_<variation>.txt``).

Usage as a CLI:
    python3 -I gemc_txt.py <geometry.txt> [--roots] [--grep REGEX] [--all] [--world]

Library:
    import sys; sys.path.insert(0, '<this dir>'); import gemc_txt
    vols = gemc_txt.load('<file>')           # OrderedDict name -> Volume
    gemc_txt.world_transform(vols, name)     # (R_active 3x3, t[cm]) in hall frame

Column layout (gemc 2.x TEXT factory):
 0 name | 1 mother | 2 description | 3 pos | 4 rot | 5 color | 6 type | 7 dims | 8 material |
 9 magfield | 10 ncopy | 11 pMany | 12 exist | 13 visible | 14 style | 15 sensitivity |
 16 hit_type | 17 identifiers

Units: lengths are converted to cm, angles to deg.

ROTATION CONVENTION (verified in gemc 2.9 source, utilities/gemcUtils.cc + detector/detector.cc):
  rot string "a b c"  ->  CLHEP rot = I; rot.rotateX(a); rot.rotateY(b); rot.rotateZ(c)
                       ->  rot = Rz(c) * Ry(b) * Rx(a)   (rotateX left-multiplies)
  then G4PVPlacement(&rot, pos, ...) treats rot as the PASSIVE (frame) rotation, so the
  ACTIVE rotation applied to the daughter is  R_active = rot^-1 = Rx(-a) * Ry(-b) * Rz(-c).
  In DD4hep (active convention):  RotationX(-a)*RotationY(-b)*RotationZ(-c)
  "ordered: zxy a b c" applies rotateZ(a), rotateX(b), rotateY(c) in that order (same left-mult rule).
  Boolean "Operation:" second solid is transformed by T(pos2) * R(rot2^-1) (same active rule),
  "Operation:@" means pos/rot of solid2 is given in the common mother frame.
Pgon dims: phi0 dphi nsides nplanes rin[n] rout[n] z[n]  (radii = tangent distance / apothem,
  same meaning as G4Polyhedra and TGeoPgon/dd4hep::Polyhedra).
Polycone dims: phi0 dphi nplanes rin[n] rout[n] z[n]
Cons dims: rmin1 rmax1 rmin2 rmax2 hz phi0 dphi (hz = HALF length)
Tube dims: rmin rmax hz phi0 dphi (hz = HALF length)
Box dims: hx hy hz (HALF lengths)
Trd dims: hx1 hx2 hy1 hy2 hz (HALF lengths, G4Trd)
Sphere dims: rmin rmax phi0 dphi theta0 dtheta
"""
import collections
import math
import re
import sys

LEN = {'mm': 0.1, 'cm': 1.0, 'm': 100.0, 'um': 1e-4, 'nm': 1e-7, 'inch': 2.54, 'in': 2.54}
ANG = {'deg': 1.0, 'rad': 180.0 / math.pi, 'mrad': 0.18 / math.pi}


def num(tok):
    """Convert '12.5*cm' -> 12.5 (cm), '30*deg' -> 30 (deg), '3*counts' -> 3, '7' -> 7."""
    tok = tok.strip()
    m = re.match(r'^([-+0-9.eE]+)\s*\*\s*([A-Za-z]+)$', tok)
    if m:
        v, u = float(m.group(1)), m.group(2)
        if u in LEN:
            return v * LEN[u]
        if u in ANG:
            return v * ANG[u]
        return v  # counts, etc.
    try:
        return float(tok)
    except ValueError:
        return tok


class Volume:
    __slots__ = ('name', 'mother', 'desc', 'pos', 'rot_str', 'color', 'type', 'dims_str', 'material',
                 'magfield', 'ncopy', 'pmany', 'exist', 'visible', 'style', 'sensitivity', 'hit_type',
                 'identifiers', 'raw')

    def __init__(self, cols):
        cols = [c.strip() for c in cols] + [''] * 18
        (self.name, self.mother, self.desc, pos, self.rot_str, self.color, self.type, self.dims_str,
         self.material, self.magfield, self.ncopy, self.pmany, self.exist, self.visible, self.style,
         self.sensitivity, self.hit_type, self.identifiers) = cols[:18]
        self.pos = [num(x) for x in pos.split()] if pos else [0.0, 0.0, 0.0]
        self.raw = cols

    @property
    def dims(self):
        return [num(x) for x in self.dims_str.split()] if self.dims_str else []

    def __repr__(self):
        return f"<{self.name} in {self.mother}: {self.type} [{self.dims_str}] {self.material} pos={self.pos} rot='{self.rot_str}'>"


def load(path):
    vols = collections.OrderedDict()
    with open(path) as f:
        for row in f:
            if not row.strip():
                continue
            cols = row.rstrip('\n').split('|')
            v = Volume(cols)
            vols[v.name] = v
    return vols


# ---------------- rotation helpers (3x3 lists) -----------------
def _mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def _mv(a, v):
    return [sum(a[i][k] * v[k] for k in range(3)) for i in range(3)]


def _T(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def Rx(d):
    c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
    return [[1, 0, 0], [0, c, -s], [0, s, c]]


def Ry(d):
    c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]


def Rz(d):
    c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]


I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def gemc_rot_matrix(rot_str):
    """Return the CLHEP matrix gemc builds (PASSIVE convention). Active = transpose."""
    toks = rot_str.split()
    rot = I3
    axis = {'x': Rx, 'y': Ry, 'z': Rz}
    if not toks:
        return rot
    if toks[0] == 'ordered:':
        order, vals = toks[1], [num(t) for t in toks[2:5]]
        for ax, v in zip(order, vals):
            rot = _mm(axis[ax](v), rot)
    elif toks[0] == 'doubleRotation:':
        vals = [num(t) for t in toks[1:7]]
        for ax, v in zip('xyzxyz', vals):
            rot = _mm(axis[ax](v), rot)
    else:
        vals = [num(t) for t in toks[:3]]
        for ax, v in zip('xyz', vals):
            rot = _mm(axis[ax](v), rot)
    return rot


def active_rotation(rot_str):
    """Active rotation (as applied to daughter in mother frame), i.e. what DD4hep needs."""
    return _T(gemc_rot_matrix(rot_str))


def world_transform(vols, name):
    """Compose placements up to 'root'. Returns (R_active, t) mapping local->hall frame (cm)."""
    R, t = I3, [0.0, 0.0, 0.0]
    chain = []
    n = name
    while n in vols and n != 'root':
        chain.append(vols[n])
        n = vols[n].mother
    for v in reversed(chain):  # from outermost to innermost
        Ra = active_rotation(v.rot_str)
        t = [t[i] + x for i, x in enumerate(_mv(R, v.pos))]
        R = _mm(R, Ra)
    return R, t


def zxy_euler_from_active(R):
    """Decompose active rotation as RotationZYX(z,y,x) = Rz*Ry*Rx (deg) for DD4hep XML convenience."""
    sy = -R[2][0]
    sy = max(-1.0, min(1.0, sy))
    y = math.degrees(math.asin(sy))
    if abs(abs(sy) - 1) > 1e-9:
        x = math.degrees(math.atan2(R[2][1], R[2][2]))
        z = math.degrees(math.atan2(R[1][0], R[0][0]))
    else:
        x = 0.0
        z = math.degrees(math.atan2(-R[0][1], R[1][1]))
    return z, y, x


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('file')
    ap.add_argument('--roots', action='store_true', help='only volumes whose mother is root')
    ap.add_argument('--grep', default=None, help='regex filter on volume name')
    ap.add_argument('--all', action='store_true', help='include Component/exist=0 volumes')
    ap.add_argument('--world', action='store_true', help='print hall-frame position + ZYX euler of active rotation')
    ap.add_argument('--limit', type=int, default=200)
    a = ap.parse_args()
    vols = load(a.file)
    print(f"# {a.file}: {len(vols)} volumes")
    types = collections.Counter(v.type.split(':')[0] for v in vols.values())
    mats = collections.Counter(v.material for v in vols.values())
    print('# shape types:', dict(types))
    print('# materials  :', dict(mats))
    print('# sensitive  :', collections.Counter((v.sensitivity, v.hit_type) for v in vols.values()))
    n = 0
    for v in vols.values():
        if a.roots and v.mother != 'root':
            continue
        if a.grep and not re.search(a.grep, v.name):
            continue
        if not a.all and (v.material == 'Component' or v.exist == '0'):
            pass
        line = f"{v.name} | mother={v.mother} | {v.type} | dims=[{v.dims_str}] | mat={v.material} | pos={v.pos} | rot='{v.rot_str}' | sens={v.sensitivity} hit={v.hit_type} id={v.identifiers} exist={v.exist}"
        if a.world:
            R, t = world_transform(vols, v.name)
            z, y, x = zxy_euler_from_active(R)
            line += f" | HALL t=({t[0]:.4f},{t[1]:.4f},{t[2]:.4f}) cm ACTIVE RotZYX=({z:.4f},{y:.4f},{x:.4f}) deg"
        print(line)
        n += 1
        if n >= a.limit:
            print(f"# ... truncated at {a.limit} (use --limit)")
            break


if __name__ == '__main__':
    main()
