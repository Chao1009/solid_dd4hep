# SPDX-License-Identifier: LGPL-3.0-or-later
"""DD4hep integration check of a map: load a compact with a SoLID_FieldMapB fragment, query
desc.field().magneticField at reference and 20000 random points (random azimuth) and compare with (1) the bilinear
interpolation of the map file, as FieldMapB does, and (2) the rzfield solution B(r, z) (on the grid nodes they
agree to the 6-digit file format; off the nodes they differ by the interpolation error, reported in air at least
2 cm from any iron or coil edge). Needs DD4hep python (eic_xl, after thissolid.sh) with numpy and scipy:
  python3 tests/dd4hep_field_query.py test/test_magnet_v4_2D.xml work/maps/solutions/solid_v4c_open.npz \
          $DETECTOR_PATH/fieldmaps/solenoid_SoLID_v4c_2D.dat models/solid_v4c.am [result.json]"""
import json
import os
import sys

import numpy as np

compact, solution, mapfile, deckfile = sys.argv[1:5]
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dd4hep  # noqa: E402
from rzfield.field import RZSolution  # noqa: E402
from rzfield.amdeck import read_deck  # noqa: E402
from rzfield.gemcmap import read_map  # noqa: E402
from rzfield.geometry import distance_to_edges, region_at  # noqa: E402

desc = dd4hep.Detector.getInstance()
desc.fromXML(compact)
fld = desc.field()
cm, G = dd4hep.cm, dd4hep.gauss
sol = RZSolution.load(solution)
m = read_map(mapfile)
r_ax, z_ax, mBr, mBz = m["r"], m["z"], m["Br"], m["Bz"]


def dd(x, y, z):
    B = fld.magneticField(dd4hep.Position(x * cm, y * cm, z * cm))
    return np.array([B.X() / G, B.Y() / G, B.Z() / G])


def bilinear(r, z):
    out = []
    for A in (mBr, mBz):
        u = (r - r_ax[0]) / (r_ax[1] - r_ax[0])
        v = (z - z_ax[0]) / (z_ax[1] - z_ax[0])
        if u < 0 or u > len(r_ax) - 1 or v < 0 or v > len(z_ax) - 1:
            return np.zeros(2)
        i = min(int(u), len(r_ax) - 2)
        j = min(int(v), len(z_ax) - 2)
        du, dv = u - i, v - j
        out.append(A[i, j] * (1 - du) * (1 - dv) + A[i, j + 1] * (1 - du) * dv + A[i + 1, j] * du * (1 - dv)
                   + A[i + 1, j + 1] * du * dv)
    return np.array(out)


def cyl(Bxyz, phi):
    """(Br, Bphi, Bz) from (Bx, By, Bz) at azimuth phi"""
    c, s = np.cos(phi), np.sin(phi)
    return np.array([Bxyz[0] * c + Bxyz[1] * s, -Bxyz[0] * s + Bxyz[1] * c, Bxyz[2]])


res = {"compact": compact, "map": mapfile, "points": []}
print("== %s\n   map %s, solution %s" % (compact, mapfile, solution))
print("   %-26s %12s %12s %12s | %12s %12s | %9s" % ("point (r, phi, z)", "Br dd4hep", "Bz dd4hep", "Bphi dd4hep",
                                                     "Br solution", "Bz solution", "max diff"))
ref = [(0, 0, 0), (0, 0, -315), (0, 0, -350), (215, 0, 344), (215, 37, 344), (215, 200, 300), (100, 123, 100),
       (150.5, 45, -20.25), (60.3, 290, 450.7), (0, 0, 799), (0, 0, -599)]
for r, phid, z in ref:
    phi = np.radians(phid)
    b = cyl(dd(r * np.cos(phi), r * np.sin(phi), z), phi)
    sbr, sbz = sol.B_cm_gauss(np.array([float(r)]), np.array([float(z)]))
    d = max(abs(b[0] - sbr[0]), abs(b[2] - sbz[0]))
    print("   (%6.1f, %4d deg, %7.2f) %12.4f %12.4f %12.2e | %12.4f %12.4f | %9.2e G" % (
        r, phid, z, b[0], b[2], b[1], sbr[0], sbz[0], d))
    res["points"].append(dict(r=r, phi_deg=phid, z=z, dd4hep=[b[0], b[1], b[2]], solution=[sbr[0], sbz[0]]))
# outside the map: no field
out = [dd(0, 0, 801), dd(0, 0, -601), dd(501, 0, 0), dd(400, 400, 0)]
res["outside_max_G"] = float(np.max(np.abs(out)))
print("   outside the map (z 801, z -601, r 501, r 566): max |B| %.3g G" % res["outside_max_G"])
# random points: dd4hep vs the bilinear interpolation of the file, and vs the solution
rng = np.random.default_rng(12345)
n = 20000
r = rng.uniform(0, 500, n)
z = rng.uniform(-600, 800, n)
phi = rng.uniform(0, 2 * np.pi, n)
on = rng.random(n) < 0.25           # a quarter exactly on grid nodes
r[on] = np.round(r[on])
z[on] = np.round(z[on])
Bd = np.array([cyl(dd(r[k] * np.cos(phi[k]), r[k] * np.sin(phi[k]), z[k]), phi[k]) for k in range(n)])
Bi = np.array([bilinear(r[k], z[k]) for k in range(n)])
sBr, sBz = sol.B_cm_gauss(r, z)
Bmag = np.hypot(Bi[:, 0], Bi[:, 1])
d_file = np.hypot(Bd[:, 0] - Bi[:, 0], Bd[:, 2] - Bi[:, 1])
d_sol = np.hypot(Bd[:, 0] - sBr, Bd[:, 2] - sBz)
rel_file = d_file / np.maximum(Bmag, 1.0)
print("   %d random points: |Bphi| max %.2e G" % (n, np.max(np.abs(Bd[:, 1]))))
print("   dd4hep vs bilinear of the file: max |dB| %.2e G, max relative %.2e" % (d_file.max(), rel_file.max()))
print("   dd4hep vs solution on grid nodes (%d): max |dB| %.3e G, max relative %.2e" % (
    on.sum(), d_sol[on].max(), (d_sol[on] / np.maximum(Bmag[on], 1.0)).max()))
# off the nodes the map differs from the solution by the bilinear interpolation error; report it in air at least
# 2 cm from any iron or coil edge (inside iron and across iron surfaces a 1 cm map is not meant to be exact)
deck = read_deck(deckfile)
reg = region_at(deck, r, z)
mats = {rg.index: (rg.mat, rg.cur) for rg in deck.regions}
air = np.array([k_ <= 1 or (mats[k_][0] == 1 and mats[k_][1] == 0.0) for k_ in reg])
air &= distance_to_edges(deck, r, z, "all") >= 2.0
res.update(n_random=n, bphi_max_G=float(np.max(np.abs(Bd[:, 1]))), file_max_dB_G=float(d_file.max()),
           file_max_rel=float(rel_file.max()), nodes_max_dB_G=float(d_sol[on].max()),
           nodes_max_rel=float((d_sol[on] / np.maximum(Bmag[on], 1.0)).max()))
for name, sel in (("air", (~on) & air), ("tracking", (~on) & air & (r < 140) & (z > -350) & (z < 200)),
                  ("pmt", (~on) & air & (r >= 200) & (r <= 250) & (z >= 280) & (z <= 360)),
                  ("iron_or_edges", (~on) & ~air)):
    d = d_sol[sel]
    rel = d / np.maximum(Bmag[sel], 1.0)
    res["interp_" + name] = dict(n=int(sel.sum()), median_G=float(np.median(d)), p99_G=float(np.percentile(d, 99)),
                                 max_G=float(d.max()), max_rel=float(rel.max()))
    print("   dd4hep vs solution off the nodes, %-13s (%5d points): median %.4f G, p99 %.3f G, max %.3f G "
          "(max relative %.2e)" % (name, sel.sum(), np.median(d), np.percentile(d, 99), d.max(), rel.max()))
if len(sys.argv) > 5:
    json.dump(res, open(sys.argv[5], "w"), indent=1, default=float)
