#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare the new maps (make_maps.sh solutions) with the POISSON references they replace.

  table         v4c vs solenoid_CLEOv9.dat, stage 1 vs solenoid_CLEOv10.dat, stage 2 vs the CLEOv11 numbers of
                "CLEO-II Stage 1 and 2.html" (its <script id="data"> block). Optional extra rzfield solutions split
                each difference into solver, boundary, iron-model and geometry parts (see magnet_field/README.md).
  cleov11-deck  rebuild the CLEOv11 deck from the html data block (the block lacks 'Upstream coil Inner', which is
                added from CLEOv10.am; the problem box is taken as CLEOv10.am's, r 0..500, z -500..700 cm).

Examples:
  python3 compare_references.py cleov11-deck "$GEMC/CLEO-II Stage 1 and 2.html" -o work/CLEOv11_html.am
  python3 compare_references.py table --gemc $GEMC --solutions work/maps/solutions --extra work/ref \\
      --json validation/new_maps.json --plot validation/new_maps_axis.png
Extra solutions looked up by name in --extra: solid_lite_stage1_box.npz (stage-1 deck in the CLEOv10 box),
CLEOv10_poisson.npz / CLEOv10_open.npz (CLEOv10.am, box / open), CLEOv9_open.npz (CLEOv9.am, 1010 steel, open),
CLEOv11_html_box.npz / CLEOv11_html_open.npz (the rebuilt CLEOv11 deck, box / open).
"""
import argparse
import json
import os
import re
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from rzfield.amdeck import read_deck  # noqa: E402
from rzfield.field import RZSolution  # noqa: E402
from rzfield.gemcmap import read_map  # noqa: E402
from rzfield.geometry import material_masks  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
POINTS = [("Bz(0,0)", 0.0, 0.0, "Bz"), ("Bz(0,-315)", 0.0, -315.0, "Bz"), ("Bz(0,-350)", 0.0, -350.0, "Bz"),
          ("|B|(215,300)", 215.0, 300.0, "B"), ("|B|(215,344)", 215.0, 344.0, "B")]
REGIONS = {"tracking": (0, 140, -350, 200), "coil bore": (0, 140, -170, 170), "PMT": (200, 250, 280, 360),
           "target He3": (0, 2, -380, -320), "target LH2": (0, 2, -322.5, -307.5)}
UP_INNER = {"name": "Upstream coil Inner", "mat": 1, "cur": 542156.0,
            "pts": [[-173.8, 152.3], [-85.5, 152.3], [-85.5, 152.7], [-173.8, 152.7], [-173.8, 152.3]]}


def html_data(path):
    """the JSON data block of the Stage 1/2 html (data only, nothing is executed)"""
    m = re.search(r"<script[^>]*\bid=\"data\"[^>]*>(.*?)</script>", open(path, encoding="utf-8").read(), re.S)
    if not m:
        raise SystemExit("no <script id=\"data\"> block in %s" % path)
    return json.loads(m.group(1))


class Source:
    """B(r, z) in gauss from a GEMC map (bilinear, exact on its nodes) or an rzfield solution"""

    def __init__(self, path, label):
        self.path, self.label = path, label
        if path.endswith(".npz"):
            self.sol = RZSolution.load(path)
            self.map = None
        else:
            self.sol = None
            self.map = read_map(path)

    def B(self, r, z):
        r = np.asarray(r, float)
        z = np.asarray(z, float)
        if self.sol is not None:
            return self.sol.B_cm_gauss(r, z)
        m = self.map
        ra, za = m["r"], m["z"]
        u = np.clip((r - ra[0]) / (ra[1] - ra[0]), 0, len(ra) - 1 - 1e-9)
        v = np.clip((z - za[0]) / (za[1] - za[0]), 0, len(za) - 1 - 1e-9)
        i, j = u.astype(int), v.astype(int)
        a, b = u - i, v - j
        out = []
        for A in (m["Br"], m["Bz"]):
            out.append(A[i, j] * (1 - a) * (1 - b) + A[i, j + 1] * (1 - a) * b + A[i + 1, j] * a * (1 - b)
                       + A[i + 1, j + 1] * a * b)
        return out[0], out[1]

    def energy_MJ(self):
        return self.sol.info["energy_J"] / 1e6 if self.sol is not None else None


def point_values(src):
    r = np.array([p[1] for p in POINTS])
    z = np.array([p[2] for p in POINTS])
    br, bz = src.B(r, z)
    return {p[0]: float(bz[k] if p[3] == "Bz" else np.hypot(br[k], bz[k])) for k, p in enumerate(POINTS)}


def axis_integral(src, z0, z1):
    z = np.arange(z0, z1 + 0.5, 1.0)
    _, bz = src.B(np.zeros_like(z), z)
    return float(np.sum(0.5 * (bz[1:] + bz[:-1])) * 1e-6)          # 1 cm steps -> T m


def region_diffs(new, ref, decks, rr, zz):
    """|B_new - B_ref| on a 1 cm grid in each region, air points >= 1 cm from the edges of every deck; also the
    largest Bz difference (signed, with its relative size) and the largest one on the axis"""
    out = {}
    for name, (r0, r1, z0, z1) in REGIONS.items():
        r = np.arange(r0, r1 + 1e-9, 1.0)
        z = np.arange(np.ceil(z0), np.floor(z1) + 1e-9, 1.0)
        z = z[(z >= zz[0]) & (z <= zz[1])]
        R, Z = np.meshgrid(r, z, indexing="ij")
        air = np.ones(R.shape, bool)
        for d in decks:          # (points on the axis are nudged off the box edge for the point-in-polygon test)
            air &= material_masks(d, np.maximum(R, 1e-3), Z, margin_cm=1.0)[2]
        b1 = new.B(R, Z)
        b0 = ref.B(R, Z)
        B0 = np.hypot(*b0)
        B1 = np.hypot(*b1)
        dB = np.hypot(b1[0] - b0[0], b1[1] - b0[1])[air]
        at = (float(R[air][np.argmax(dB)]), float(Z[air][np.argmax(dB)])) if dB.size else None
        out[name] = dict(n=int(air.sum()), rms_dB_G=float(np.sqrt(np.mean(dB ** 2))) if dB.size else None,
                         max_dB_G=float(dB.max()) if dB.size else None, max_at_rz_cm=at,
                         ref_B_G=[float(B0[air].min()), float(B0[air].max())] if dB.size else None,
                         new_B_G=[float(B1[air].min()), float(B1[air].max())] if dB.size else None)
        if dB.size:
            dBz = (b1[1] - b0[1])[air]
            k = int(np.argmax(np.abs(dBz)))
            bz0 = b0[1][air][k]
            out[name].update(max_dBz_G=float(dBz[k]), max_dBz_rel=float(dBz[k] / bz0) if bz0 else None,
                             max_dBz_at_rz_cm=(float(R[air][k]), float(Z[air][k])))
            ax = air & (R == 0.0)
            if ax.any():
                dza = (b1[1] - b0[1])[ax]
                ka = int(np.argmax(np.abs(dza)))
                out[name].update(axis_max_dBz_G=float(dza[ka]), axis_max_dBz_at_z_cm=float(Z[ax][ka]))
    return out


def fmt(v, nd=1):
    return "-" if v is None else ("%.*f" % (nd, v))


def run_table(a):
    G = a.gemc
    sols = a.solutions
    extra = {}
    for d in a.extra:
        for f in os.listdir(d):
            if f.endswith(".npz"):
                extra.setdefault(f[:-4], os.path.join(d, f))
    new = {k: Source(os.path.join(sols, "%s_open.npz" % k), k) for k in ("solid_v4c", "solid_lite_stage1",
                                                                          "solid_lite_stage2")}
    decks = {k: read_deck(os.path.join(HERE, "models", k + ".am")) for k in new}
    res = {"points": [p[0] for p in POINTS], "regions": REGIONS, "comparisons": {}}

    def opt(name, label):
        return Source(extra[name], label) if name in extra else None

    def column(src):
        if src is None:
            return None
        c = dict(label=src.label, file=os.path.basename(src.path), points=point_values(src))
        if src.energy_MJ() is not None:
            c["energy_MJ"] = src.energy_MJ()
        return c

    # ---- stage 1 vs CLEOv10
    ref10 = Source(os.path.join(G, "solenoid_CLEOv10.dat"), "POISSON CLEOv10 map (box)")
    d10 = read_deck(os.path.join(G, "CLEOv10.am"))
    box1 = opt("solid_lite_stage1_box", "rzfield stage-1 deck, box")
    srcs = [s for s in (ref10, opt("CLEOv10_poisson", "rzfield CLEOv10.am, box"), box1,
                        opt("CLEOv10_open", "rzfield CLEOv10.am, open"), new["solid_lite_stage1"]) if s]
    cols = [column(s) for s in srcs]
    for c, s in zip(cols, srcs):
        c["int_Bz_dz_Tm"] = {"-500..700": axis_integral(s, -500, 700)}
        c["target"] = target_stats(s)
    comp = dict(columns=cols, regions_new_vs_ref=region_diffs(new["solid_lite_stage1"], ref10,
                                                              [decks["solid_lite_stage1"], d10], (0, 500), (-500, 700)))
    if box1:
        comp["regions_box_vs_ref"] = region_diffs(box1, ref10, [decks["solid_lite_stage1"], d10], (0, 500),
                                                  (-500, 700))
    res["comparisons"]["stage1_vs_CLEOv10"] = comp

    # ---- v4c vs CLEOv9
    ref9 = Source(os.path.join(G, "solenoid_CLEOv9.dat"), "POISSON CLEOv9 map (box, mu_r frozen ~250)")
    d9 = read_deck(os.path.join(G, "CLEOv9.am"))
    open9 = opt("CLEOv9_open", "rzfield CLEOv9.am, 1010 steel, open")
    srcs = [s for s in (ref9, open9, new["solid_v4c"]) if s]
    cols = [column(s) for s in srcs]
    for c, s in zip(cols, srcs):
        c["int_Bz_dz_Tm"] = {"-600..600": axis_integral(s, -600, 600)}
        c["target"] = target_stats(s)
    comp = dict(columns=cols, regions_new_vs_ref=region_diffs(new["solid_v4c"], ref9, [decks["solid_v4c"], d9],
                                                              (0, 500), (-600, 600)))
    if open9:
        comp["regions_new_vs_CLEOv9_open"] = region_diffs(new["solid_v4c"], open9, [decks["solid_v4c"], d9],
                                                          (0, 500), (-600, 600))
    res["comparisons"]["v4c_vs_CLEOv9"] = comp

    # ---- stage 2 vs CLEOv11 (html)
    h = html_data(a.html) if a.html else html_data(os.path.join(G, "CLEO-II Stage 1 and 2.html"))
    fin, stray, pmt = h["final"], h["stray"], h["pmt"]
    html_pts = {"Bz(0,0)": fin["bz0"] * 1e4, "Bz(0,-350)": dict(fin["axis"])[-350],
                "Bz(0,-315)": dict(fin["axis"])[-315], "|B|(215,300)": pmt["v11"][0][2],
                "|B|(215,344)": pmt["v11"][1][2]}
    srcs = [opt("CLEOv11_html_box", "rzfield CLEOv11 (html data), box"),
            opt("CLEOv11_html_open", "rzfield CLEOv11 (html data), open"), new["solid_lite_stage2"]]
    cols = [dict(label="POISSON CLEOv11 (html, 2026-10-31)", points=html_pts, energy_MJ=fin["energy"] / 1e6,
                 target={"He3": {"max_B_G": fin.get("maxB")}, "LH2": {"max_B_G": None}})]
    spts = [tuple(p) for p in stray["pts"]]
    srow = {r["case"]: r["B"] for r in stray["rows"]}
    cols[0]["stray"] = srow["final"]
    za = np.array([p[0] for p in fin["axis"]], float)
    bza = np.array([p[1] for p in fin["axis"]], float)
    for s in srcs:
        c = column(s)
        if c is None:
            continue
        br, bz = s.B(np.array([p[0] for p in spts], float), np.array([p[1] for p in spts], float))
        c["stray"] = [float(x) for x in np.hypot(br, bz)]
        c["target"] = target_stats(s)
        _, bz_ax = s.B(np.zeros_like(za), za)
        c["axis_-400..-300_vs_html"] = dict(rms_G=float(np.sqrt(np.mean((bz_ax - bza) ** 2))),
                                            max_G=float(np.max(np.abs(bz_ax - bza))))
        cols.append(c)
    res["comparisons"]["stage2_vs_CLEOv11"] = dict(columns=cols, stray_points=spts,
                                                   stray_v10_html=srow.get("v10"))
    print_tables(res)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    if a.plot:
        plot_axis(a.plot, new, ref9, ref10, za, bza)
    return 0


def target_stats(src):
    """on-axis Bz and max |B| (r <= 2 cm) over the He3 cell (z -350 +- 30) and the LH2 cell (-315 +- 7.5)"""
    out = {}
    for name, (r0, r1, z0, z1) in (("He3", REGIONS["target He3"]), ("LH2", REGIONS["target LH2"])):
        z = np.arange(np.ceil(z0), np.floor(z1) + 1e-9, 1.0)
        R, Z = np.meshgrid(np.arange(r0, r1 + 1e-9, 1.0), z, indexing="ij")
        br, bz = src.B(R, Z)
        out[name] = dict(Bz_axis_G=[float(bz[0].min()), float(bz[0].max())], max_B_G=float(np.hypot(br, bz).max()))
    return out


def print_tables(res):
    for key, comp in res["comparisons"].items():
        print("== %s" % key)
        cols = comp["columns"]
        print("  %-14s" % "" + "".join("%22s" % c["label"][:21] for c in cols))
        for p in res["points"]:
            print("  %-14s" % p + "".join("%22s" % fmt(c["points"].get(p), 2) for c in cols))
        if any("int_Bz_dz_Tm" in c for c in cols):
            rng = list(next(c for c in cols if "int_Bz_dz_Tm" in c)["int_Bz_dz_Tm"])[0]
            print("  %-14s" % ("intBz " + rng) + "".join("%22s" % fmt(c.get("int_Bz_dz_Tm", {}).get(rng), 4)
                                                         for c in cols))
        print("  %-14s" % "energy MJ" + "".join("%22s" % fmt(c.get("energy_MJ"), 3) for c in cols))
        for t in ("He3", "LH2"):
            if any("target" in c for c in cols):
                print("  %-14s" % ("max|B| " + t) + "".join(
                    "%22s" % fmt(c["target"][t]["max_B_G"] if "target" in c else None, 2) for c in cols))
        if "stray" in cols[0]:
            for k, p in enumerate(comp["stray_points"]):
                print("  %-14s" % ("|B|(%g,%g)" % p) + "".join("%22s" % fmt(c["stray"][k], 1) for c in cols))
        for rk in [k for k in comp if k.startswith("regions")]:
            print("  %s:" % rk)
            for name, s in comp[rk].items():
                if s["n"]:
                    print("    %-11s n %6d  rms |dB| %8.2f G  max %8.2f G at (%g, %g)  |B_ref| %s  |B_new| %s G" % (
                        name, s["n"], s["rms_dB_G"], s["max_dB_G"], s["max_at_rz_cm"][0], s["max_at_rz_cm"][1],
                        "%.1f..%.1f" % tuple(s["ref_B_G"]), "%.1f..%.1f" % tuple(s["new_B_G"])))
                    extra = "                max dBz %+.1f G (%s) at (%g, %g)" % (
                        s["max_dBz_G"], "%+.1f %%" % (100 * s["max_dBz_rel"]) if s.get("max_dBz_rel") is not None
                        else "-", s["max_dBz_at_rz_cm"][0], s["max_dBz_at_rz_cm"][1])
                    if "axis_max_dBz_G" in s:
                        extra += "; on the axis %+.1f G at z %g" % (s["axis_max_dBz_G"], s["axis_max_dBz_at_z_cm"])
                    print(extra)


def plot_axis(path, new, ref9, ref10, za, bza):
    """on-axis Bz of the new maps (solid) and of the POISSON references (dashed), same hue per magnet"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid, surf = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    col = {"solid_v4c": "#2a78d6", "solid_lite_stage1": "#eb6834", "solid_lite_stage2": "#1baf7a"}
    lab = {"solid_v4c": "full SoLID v4c (new)", "solid_lite_stage1": "SoLID-Lite stage 1 (new)",
           "solid_lite_stage2": "SoLID-Lite stage 2 (new)"}
    z = np.arange(-600.0, 800.5, 1.0)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 6.4), gridspec_kw={"height_ratios": [1.25, 1]},
                                   facecolor=surf)
    for ax in (ax1, ax2):
        ax.set_facecolor(surf)
        ax.grid(color=grid, lw=0.6)
        for s in ax.spines.values():
            s.set_color(grid)
        ax.tick_params(colors=ink2, labelsize=8)
    for k, src in new.items():
        _, bz = src.B(np.zeros_like(z), z)
        ax1.plot(z, bz / 1e4, color=col[k], lw=2, label=lab[k])
        ax2.semilogy(z, np.maximum(np.abs(bz), 1e-2), color=col[k], lw=2)
    for src, k, txt in ((ref9, "solid_v4c", "POISSON CLEOv9 map"), (ref10, "solid_lite_stage1", "POISSON CLEOv10 map")):
        zz = src.map["z"]
        bz = src.map["Bz"][0]
        ax1.plot(zz, bz / 1e4, color=col[k], lw=1.2, ls=(0, (4, 2)), label=txt)
        ax2.semilogy(zz, np.maximum(np.abs(bz), 1e-2), color=col[k], lw=1.2, ls=(0, (4, 2)))
    ax2.semilogy(za, np.maximum(np.abs(bza), 1e-2), color=col["solid_lite_stage2"], lw=1.2, ls=(0, (1, 1.5)),
                 label="POISSON CLEOv11 (html)")
    ax1.set_ylabel("Bz on axis [T]", color=ink2, fontsize=9)
    ax1.set_xlim(-600, 800)
    ax1.legend(frameon=False, fontsize=8, labelcolor=ink, loc="upper right")
    ax1.set_title("On-axis Bz: new 2D maps (solid) and the POISSON maps they replace (dashed)", color=ink,
                  fontsize=10, loc="left")
    ax2.set_xlim(-600, -150)
    ax2.set_ylim(1e-2, 2e4)
    ax2.axvspan(-380, -320, color="#cde2fb", alpha=0.6, lw=0)
    ax2.text(-350, 1.2e4, "3He cell", ha="center", va="top", fontsize=8, color=ink2)
    ax2.set_xlabel("z [cm]", color=ink2, fontsize=9)
    ax2.set_ylabel("|Bz| on axis [G]", color=ink2, fontsize=9)
    ax2.legend(frameon=False, fontsize=8, labelcolor=ink, loc="upper left")
    ax2.set_title("Upstream fringe field and target region (log scale)", color=ink, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=90, facecolor=surf)
    plt.close(fig)
    print("plot written to %s" % path)


def run_deck(a):
    h = html_data(a.html)
    regs = [UP_INNER] + list(h["geometry"])
    out = ["CLEOv11 (stage 2A) rebuilt from the data block of 'CLEO-II Stage 1 and 2.html'",
           "; 'Upstream coil Inner' added from CLEOv10.am; problem box assumed = CLEOv10.am (r 0..500, z -500..700)",
           "&reg kprob=0, icylin=1, ienergy=1, xjfact=1.0, dx=1.0, rhogam=0.02, mode=0,",
           "nbsup=0, nbslf=0, nbsrt=0, nbslo=0, ibound=0 &", "; Hall  Region 1"]
    out += ["&po x=%8.3f, y=%9.3f &" % (r, z) for z, r in [(-500, 0), (700, 0), (700, 500), (-500, 500), (-500, 0)]]
    for i, reg in enumerate(regs, start=2):
        out += ["", "; %s  Region %d" % (reg["name"], i),
                "&reg mat=1, cur=%.1f &" % reg["cur"] if reg["mat"] == 1 else "&reg mat=%d &" % reg["mat"]]
        pts = reg["pts"] if reg["pts"][0] == reg["pts"][-1] else reg["pts"] + [reg["pts"][0]]
        out += ["&po x=%8.3f, y=%9.3f &" % (r, z) for z, r in pts]
    with open(a.output, "w") as f:
        f.write("\n".join(out) + "\n")
    print("wrote %s: %d regions + box, %.0f A-turns" % (a.output, len(regs),
                                                       sum(r["cur"] for r in regs if r["mat"] == 1)))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("table", help="new maps vs the POISSON references")
    t.add_argument("--gemc", required=True, help="solid_gemc/field (CLEOv9/10 decks and maps, the Stage 1/2 html)")
    t.add_argument("--solutions", required=True, help="directory with solid_{v4c,lite_stage1,lite_stage2}_open.npz")
    t.add_argument("--extra", action="append", default=[], help="directories with the optional extra solutions")
    t.add_argument("--html", help="Stage 1/2 html (default: <gemc>/CLEO-II Stage 1 and 2.html)")
    t.add_argument("--json", help="write all numbers here")
    t.add_argument("--plot", help="write the on-axis plot here (png)")
    d = sub.add_parser("cleov11-deck", help="rebuild the CLEOv11 deck from the html data block")
    d.add_argument("html")
    d.add_argument("-o", "--output", required=True)
    a = ap.parse_args(argv)
    return run_table(a) if a.cmd == "table" else run_deck(a)


if __name__ == "__main__":
    sys.exit(main())
