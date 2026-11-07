# SPDX-License-Identifier: LGPL-3.0-or-later
"""Do tracks bend as the map says? Compares the GEM hits of the primary particle in a ddsim output (edm4hep) with
a trajectory integrated independently (RK4 in vacuum, bilinear field from the GEMC map), with a straight line and
with the trajectory in the reversed field. With charged geantinos (no material effects; run ddsim with
--gun.particle chargedgeantino --filter.tracker '') the hits must lie on the integrated trajectory; with electrons
multiple scattering and bremsstrahlung spread them by mm. Needs numpy and uproot (eic_xl):
  python3 tests/ddsim_bending.py out.edm4hep.root $DETECTOR_PATH/fieldmaps/solenoid_SoLID_v4c_2D.dat [result.json]"""
import json
import os
import sys

import numpy as np
import uproot

out_file, map_file = sys.argv[1:3]
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rzfield.gemcmap import read_map  # noqa: E402

m = read_map(map_file)
r_ax, z_ax, mBr, mBz = m["r"], m["z"], m["Br"] * 1e-4, m["Bz"] * 1e-4      # tesla
dr, dz = r_ax[1] - r_ax[0], z_ax[1] - z_ax[0]


def field(p_cm, sign=1.0):
    x, y, z = p_cm
    r = np.hypot(x, y)
    u, v = (r - r_ax[0]) / dr, (z - z_ax[0]) / dz
    if u < 0 or u > len(r_ax) - 1 or v < 0 or v > len(z_ax) - 1:
        return np.zeros(3)
    i, j = min(int(u), len(r_ax) - 2), min(int(v), len(z_ax) - 2)
    a, b = u - i, v - j
    w = np.array([(1 - a) * (1 - b), (1 - a) * b, a * (1 - b), a * b])
    br = w @ np.array([mBr[i, j], mBr[i, j + 1], mBr[i + 1, j], mBr[i + 1, j + 1]])
    bz = w @ np.array([mBz[i, j], mBz[i, j + 1], mBz[i + 1, j], mBz[i + 1, j + 1]])
    c, s = (x / r, y / r) if r > 0 else (0.0, 0.0)
    return sign * np.array([br * c, br * s, bz])


def track(x0_cm, p0_GeV, q, zs_cm, sign=1.0, ds=0.5):
    """positions (cm) where the trajectory crosses each z in zs (increasing z only), RK4 in path length (cm)"""
    k = 0.299792458e-2 * q        # dp/ds [GeV/cm] = k * (p_hat x B[T])
    pmag = np.linalg.norm(p0_GeV)

    def f(y):
        u = y[3:] / pmag
        return np.concatenate([u, k * np.cross(u, field(y[:3], sign))])

    y = np.concatenate([x0_cm, p0_GeV]).astype(float)
    res = {}
    targets = sorted(set(zs_cm))
    ti = 0
    for _ in range(200000):
        if ti >= len(targets):
            break
        k1 = f(y)
        k2 = f(y + 0.5 * ds * k1)
        k3 = f(y + 0.5 * ds * k2)
        k4 = f(y + ds * k3)
        yn = y + ds / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        while ti < len(targets) and y[2] <= targets[ti] < yn[2]:
            t = (targets[ti] - y[2]) / (yn[2] - y[2])
            res[targets[ti]] = y[:3] + t * (yn[:3] - y[:3])
            ti += 1
        y = yn
        if y[5] <= 0:
            break
    return res


t = uproot.open(out_file)["events"]
mc = t.arrays(["MCParticles.generatorStatus", "MCParticles.PDG", "MCParticles.charge", "MCParticles.vertex.x",
               "MCParticles.vertex.y", "MCParticles.vertex.z", "MCParticles.momentum.x", "MCParticles.momentum.y",
               "MCParticles.momentum.z"], library="np")
colls = sorted({k.split("/")[0] for k in t.keys() if k.startswith("GEMTracker") and k.split("/")[0].endswith("Hits")})
hits = {}
for c in colls:
    a = t.arrays(["%s.position.%s" % (c, x) for x in "xyz"] + ["%s.momentum.%s" % (c, x) for x in "xyz"]
                 + ["_%s_particle.index" % c], library="np")
    for x in "xyz":
        for w in ("position", "momentum"):
            hits.setdefault("GEMTrackerHits.%s.%s" % (w, x), []).append(a["%s.%s.%s" % (c, w, x)])
    hits.setdefault("_GEMTrackerHits_particle.index", []).append(a["_%s_particle.index" % c])
nev = len(mc["MCParticles.PDG"])
hits = {k: [np.concatenate([v[ic][ev] for ic in range(len(colls))]) for ev in range(nev)] for k, v in hits.items()}
rows = []
for ev in range(len(mc["MCParticles.PDG"])):
    gs = mc["MCParticles.generatorStatus"][ev]
    ip = int(np.where(gs == 1)[0][0])
    q = float(mc["MCParticles.charge"][ev][ip])
    x0 = np.array([mc["MCParticles.vertex." + c][ev][ip] for c in "xyz"]) / 10.0
    p0 = np.array([mc["MCParticles.momentum." + c][ev][ip] for c in "xyz"])
    sel = hits["_GEMTrackerHits_particle.index"][ev] == ip
    hx = np.column_stack([hits["GEMTrackerHits.position." + c][ev][sel] for c in "xyz"]) / 10.0
    hp = np.column_stack([hits["GEMTrackerHits.momentum." + c][ev][sel] for c in "xyz"])
    if len(hx) == 0:
        print("event %d: no GEM hit of the primary (theta %.1f deg)" % (ev, np.degrees(np.arctan2(np.hypot(*p0[:2]), p0[2]))))
        continue
    zs = list(hx[:, 2])
    pred = {s: track(x0, p0, q, zs, sign=s) for s in (1.0, 0.0, -1.0)}
    for h, hpp in zip(hx, hp):
        d = {s: np.hypot(*(pred[s][h[2]][:2] - h[:2])) if h[2] in pred[s] else np.nan for s in pred}
        phi_h = np.arctan2(h[1], h[0])
        phi_s = np.arctan2(pred[0.0][h[2]][1], pred[0.0][h[2]][0])
        dphi = (phi_h - phi_s + np.pi) % (2 * np.pi) - np.pi
        rows.append(dict(event=ev, q=q, theta_deg=float(np.degrees(np.arctan2(np.hypot(*p0[:2]), p0[2]))),
                         p_GeV=float(np.linalg.norm(p0)), p_hit_GeV=float(np.linalg.norm(hpp)), z_cm=float(h[2]),
                         r_cm=float(np.hypot(*h[:2])), d_field_cm=float(d[1.0]), d_nofield_cm=float(d[0.0]),
                         d_reversed_cm=float(d[-1.0]), dphi_vs_straight_deg=float(np.degrees(dphi))))
df = {k: np.array([r[k] for r in rows]) for k in rows[0]}
good = df["p_hit_GeV"] > 0.98 * df["p_GeV"]          # no hard bremsstrahlung before the hit
print("%s (%s): %d events with hits, %d GEM hits of the primary (%d with p_hit > 0.98 p_gun)" % (
    out_file, " ".join(colls), len(set(df["event"])), len(rows), good.sum()))
print("  hit - trajectory through the map:  median %.3f cm, 90%% %.3f cm, max %.3f cm (p_hit > 0.98 p: max %.3f cm)" % (
    np.median(df["d_field_cm"]), np.percentile(df["d_field_cm"], 90), np.max(df["d_field_cm"]),
    np.max(df["d_field_cm"][good])))
print("  hit - straight line (no field):    median %.2f cm, min %.2f cm, max %.2f cm" % (
    np.median(df["d_nofield_cm"]), np.min(df["d_nofield_cm"]), np.max(df["d_nofield_cm"])))
print("  hit - trajectory, reversed field:  median %.2f cm, min %.2f cm" % (
    np.median(df["d_reversed_cm"]), np.min(df["d_reversed_cm"])))
print("  azimuth of the hit minus the straight line: %.2f..%.2f deg (charge %s, Bz(0,0) > 0)" % (
    df["dphi_vs_straight_deg"].min(), df["dphi_vs_straight_deg"].max(), sorted(set(df["q"]))))
print("  hits at z %.0f..%.0f cm, r %.0f..%.0f cm; theta %.1f..%.1f deg" % (
    df["z_cm"].min(), df["z_cm"].max(), df["r_cm"].min(), df["r_cm"].max(), df["theta_deg"].min(),
    df["theta_deg"].max()))
if len(sys.argv) > 3:
    json.dump(dict(rows=rows), open(sys.argv[3], "w"), indent=1)
