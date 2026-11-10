#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare an r-z solution with a reference field (a Poisson GEMC map or another solution) on the map grid.

Reports the reference points (Bz on the axis at z = 0, -315, -350 cm, |B| at r 215 z 344 cm, integral of Bz dz on
the axis) and |dB| = |B - B_ref| statistics in
  (i)   the tracking region r < 140, -350 < z < 200 cm,
  (ii)  the Cherenkov PMT region r 200-250, z 280-360 cm (both: air points >= 1 cm from any iron or coil edge),
  (iii) the iron (grid points inside MAT >= 2 regions, >= 1 cm from their edges),
  (iv)  all air >= 1 cm from any region edge,
and draws the on-axis Bz overlay and the |dB| map.

Examples:
  python3 validate_poisson.py --deck CLEOv10.am --ref solenoid_CLEOv10.dat --solution work/CLEOv10_poisson.npz \\
      --label CLEOv10 --plots validation --json work/CLEOv10_validation.json
  # solution against solution (mesh refinement, open vs poisson boundary):
  python3 validate_poisson.py --deck CLEOv10.am --ref work/CLEOv10_fine.npz --solution work/CLEOv10_poisson.npz \\
      --label "CLEOv10 mesh" --ref-label "refined mesh"
"""
import argparse
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from rzfield.amdeck import read_deck  # noqa: E402
from rzfield.driver import join_negative_values, parse_range  # noqa: E402
from rzfield.field import RZSolution  # noqa: E402
from rzfield.gemcmap import grid_axis, read_map  # noqa: E402
from rzfield.geometry import material_masks  # noqa: E402

REF_POINTS = [("Bz(0,0)", 0.0, 0.0, "Bz"), ("Bz(0,-315)", 0.0, -315.0, "Bz"), ("Bz(0,-350)", 0.0, -350.0, "Bz"),
              ("|B|(215,344)", 215.0, 344.0, "B")]


def region_stats(dB, Bref, mask, rel_target=0.01, abs_target=5.0, low=1000.0):
    d = dB[mask]
    b = Bref[mask]
    if d.size == 0:
        return {"n": 0}
    hi = b >= low
    ok = np.where(hi, d <= rel_target * b, d <= abs_target)
    out = {"n": int(d.size), "max_dB_G": float(d.max()), "rms_dB_G": float(np.sqrt(np.mean(d ** 2))),
           "p99_dB_G": float(np.percentile(d, 99)), "median_dB_G": float(np.median(d)),
           "Bref_range_G": [float(b.min()), float(b.max())], "n_ge_1kG": int(hi.sum()),
           "frac_within_target": float(ok.mean())}
    if hi.any():
        rel = d[hi] / b[hi]
        out.update({"max_rel_ge_1kG": float(rel.max()), "rms_rel_ge_1kG": float(np.sqrt(np.mean(rel ** 2))),
                    "p99_rel_ge_1kG": float(np.percentile(rel, 99))})
    if (~hi).any():
        out.update({"max_dB_lt_1kG_G": float(d[~hi].max()), "p99_dB_lt_1kG_G": float(np.percentile(d[~hi], 99))})
    rel_all = d / np.maximum(b, 1e-6)
    out["max_rel_all"] = float(rel_all.max())
    out["rms_rel_all"] = float(np.sqrt(np.mean(rel_all ** 2)))
    return out


def load_reference(path, deck, grid_r=None, grid_z=None, smooth=True):
    if path.endswith(".npz"):
        sol = RZSolution.load(path)
        b = deck.bbox("box")
        r = grid_axis(*(parse_range(grid_r) if grid_r else (b[0], b[1], 1.0)))
        z = grid_axis(*(parse_range(grid_z) if grid_z else (b[2], b[3], 1.0)))
        R, Z = np.meshgrid(r, z, indexing="ij")
        Br, Bz = sol.B_cm_gauss(R, Z, smooth=smooth)
        return {"name": os.path.basename(path), "r": r, "z": z, "Br": Br, "Bz": Bz, "solution": sol}
    return read_map(path)


def axis_integral_Tm(z_cm, bz_G):
    return float(np.trapezoid(bz_G, z_cm) * 1e-6) if hasattr(np, "trapezoid") else float(np.trapz(bz_G, z_cm) * 1e-6)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deck", required=True, help="deck used for the iron outlines / iron mask")
    ap.add_argument("--ref", required=True, help="reference: GEMC map (.dat) or solution (.npz)")
    ap.add_argument("--solution", required=True, help="solution (.npz) to test")
    ap.add_argument("--label", default=None, help="label for plots and tables")
    ap.add_argument("--ref-label", default="POISSON map", help="label of the reference")
    ap.add_argument("--sol-label", default="rzfield", help="label of the tested solution")
    ap.add_argument("--grid-r", default=None, help="grid for .npz references, min:max:step cm")
    ap.add_argument("--grid-z", default=None, help="grid for .npz references, min:max:step cm")
    ap.add_argument("--tracking", default="0,140,-350,200", help="rmin,rmax,zmin,zmax of region (i), cm")
    ap.add_argument("--pmt", default="200,250,280,360", help="rmin,rmax,zmin,zmax of region (ii), cm")
    ap.add_argument("--plots", default=None, help="directory for the PNG plots")
    ap.add_argument("--json", default=None, help="write all numbers to this JSON file")
    ap.add_argument("--no-diffmap", action="store_true", help="only the axis plot")
    ap.add_argument("--raw", action="store_true", help="use the raw element field instead of the recovered one")
    a = ap.parse_args(join_negative_values(sys.argv[1:] if argv is None else argv, ("--grid-r", "--grid-z")))

    deck = read_deck(a.deck)
    ref = load_reference(a.ref, deck, a.grid_r, a.grid_z, smooth=not a.raw)
    sol = RZSolution.load(a.solution)
    label = a.label or os.path.splitext(os.path.basename(a.solution))[0]
    r, z = ref["r"], ref["z"]
    R, Z = np.meshgrid(r, z, indexing="ij")
    Br, Bz = sol.B_cm_gauss(R, Z, smooth=not a.raw)
    if not (np.all(np.isfinite(Br)) and np.all(np.isfinite(Bz))):
        raise SystemExit("the solution does not cover the reference grid")
    Bref = np.hypot(ref["Br"], ref["Bz"])
    dB = np.hypot(Br - ref["Br"], Bz - ref["Bz"])
    res = {"label": label, "reference": os.path.basename(a.ref), "solution": os.path.basename(a.solution),
           "deck": os.path.basename(a.deck),
           "grid": {"r": [float(r[0]), float(r[-1]), len(r)], "z": [float(z[0]), float(z[-1]), len(z)]},
           "field": "raw" if a.raw else "recovered",
           "solution_info": {k: sol.info.get(k) for k in ("mode", "scale", "n_elements", "n_nodes", "energy_J",
                                                            "converged", "interp", "overrides", "R_out_cm")}}
    if "newton" in sol.info:
        res["solution_info"]["newton_iterations"] = len(sol.info["newton"])
        res["solution_info"]["final_residual"] = sol.info["newton"][-1]["res"]
    # reference points (exactly on grid points for maps; evaluated directly for both if off-grid)
    pts = []
    for name, pr, pz, what in REF_POINTS:
        if not (r[0] <= pr <= r[-1] and z[0] <= pz <= z[-1]):
            continue
        i = int(np.argmin(np.abs(r - pr)))
        j = int(np.argmin(np.abs(z - pz)))
        on_grid = abs(r[i] - pr) < 1e-6 and abs(z[j] - pz) < 1e-6
        if on_grid:
            vref = ref["Bz"][i, j] if what == "Bz" else Bref[i, j]
        elif "solution" in ref:
            br0, bz0 = ref["solution"].B_cm_gauss(np.array([pr]), np.array([pz]), smooth=not a.raw)
            vref = bz0[0] if what == "Bz" else np.hypot(br0[0], bz0[0])
        else:
            continue
        br1, bz1 = sol.B_cm_gauss(np.array([pr]), np.array([pz]), smooth=not a.raw)
        v = bz1[0] if what == "Bz" else np.hypot(br1[0], bz1[0])
        pts.append({"name": name, "ref_G": float(vref), "sol_G": float(v), "diff_G": float(v - vref),
                    "rel": float((v - vref) / vref) if vref != 0 else None})
    i0 = int(np.argmin(np.abs(r)))
    res["points"] = pts
    res["axis_integral_Tm"] = {"ref": axis_integral_Tm(z, ref["Bz"][i0]), "sol": axis_integral_Tm(z, Bz[i0])}
    dz_axis = Bz[i0] - ref["Bz"][i0]
    res["axis"] = {"max_abs_dBz_G": float(np.max(np.abs(dz_axis))),
                   "z_of_max_G": float(z[int(np.argmax(np.abs(dz_axis)))]),
                   "max_rel_where_ge_1kG": float(np.max(np.abs(dz_axis)[np.abs(ref["Bz"][i0]) >= 1000]
                                                        / np.abs(ref["Bz"][i0])[np.abs(ref["Bz"][i0]) >= 1000]))}
    # regions
    t = [float(x) for x in a.tracking.split(",")]
    p = [float(x) for x in a.pmt.split(",")]
    m_trk = (R < t[1]) & (R >= t[0]) & (Z > t[2]) & (Z < t[3])
    m_pmt = (R >= p[0]) & (R <= p[1]) & (Z >= p[2]) & (Z <= p[3])
    iron, m_iron, m_air = material_masks(deck, R, Z, margin_cm=1.0)
    # detector regions: air points only (>= 1 cm from any iron or coil edge); excluded points are counted
    res["regions"] = {
        "tracking": dict(region_stats(dB, Bref, m_trk & m_air), box=t, n_excluded=int((m_trk & ~m_air).sum())),
        "pmt": dict(region_stats(dB, Bref, m_pmt & m_air), box=p, n_excluded=int((m_pmt & ~m_air).sum())),
        "iron_1cm": region_stats(dB, Bref, m_iron),
        "air_1cm": region_stats(dB, Bref, m_air),
        "near_edges_or_in_coils": region_stats(dB, Bref, ~(m_iron | m_air)),
    }
    # print summary
    print("== %s: %s vs %s" % (label, a.sol_label, a.ref_label))
    for q in pts:
        print("  %-14s ref %10.2f  sol %10.2f  diff %+9.2f G (%+.3f %%)" % (
            q["name"], q["ref_G"], q["sol_G"], q["diff_G"], 100 * q["rel"] if q["rel"] is not None else np.nan))
    ai = res["axis_integral_Tm"]
    print("  int Bz dz (axis, z %g..%g): ref %.4f  sol %.4f T m (%+.3f %%)" % (
        z[0], z[-1], ai["ref"], ai["sol"], 100 * (ai["sol"] / ai["ref"] - 1)))
    print("  axis: max |dBz| %.2f G at z %g; max rel (|Bz| >= 1 kG) %.3f %%" % (
        res["axis"]["max_abs_dBz_G"], res["axis"]["z_of_max_G"], 100 * res["axis"]["max_rel_where_ge_1kG"]))
    for k, s in res["regions"].items():
        if s.get("n", 0) == 0:
            continue
        line = "  %-9s n %6d  |dB| max %8.2f rms %7.2f p99 %7.2f G" % (k, s["n"], s["max_dB_G"], s["rms_dB_G"],
                                                                      s["p99_dB_G"])
        if s.get("n_excluded"):
            line += " (%d iron/edge points excluded)" % s["n_excluded"]
        if "max_rel_ge_1kG" in s:
            line += " | rel(>=1kG) max %.3f%% rms %.3f%%" % (100 * s["max_rel_ge_1kG"], 100 * s["rms_rel_ge_1kG"])
        if "max_dB_lt_1kG_G" in s:
            line += " | <1kG max %.2f G" % s["max_dB_lt_1kG_G"]
        line += " | within target %.1f%%" % (100 * s["frac_within_target"])
        print(line)
    if a.plots:
        os.makedirs(a.plots, exist_ok=True)
        stem = os.path.join(a.plots, label.replace(" ", "_").replace("/", "_"))
        files = make_plots(stem, label, a.ref_label, a.sol_label, deck, r, z, Br, Bz, ref, dB, Bref, t, p,
                           not a.no_diffmap)
        # recorded as <plot directory name>/<file>, e.g. validation/CLEOv10_axis.png (no host directories)
        res["plots"] = [os.path.join(os.path.basename(os.path.normpath(a.plots)), os.path.basename(f)) for f in files]
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    return 0


def make_plots(stem, label, ref_label, sol_label, deck, r, z, Br, Bz, ref, dB, Bref, t, p, diffmap=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.patches import Polygon, Rectangle

    files = []
    i0 = int(np.argmin(np.abs(r)))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 5.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax1.plot(z, ref["Bz"][i0] / 1e4, color="#1f4e79", lw=2.2, label=ref_label)
    ax1.plot(z, Bz[i0] / 1e4, color="#e07b00", lw=1.0, ls="--", label=sol_label)
    ax1.set_ylabel("Bz on axis [T]")
    ax1.legend(frameon=False, loc="upper right")
    ax1.set_title("%s: on-axis Bz" % label)
    ax1.grid(alpha=0.3)
    d = Bz[i0] - ref["Bz"][i0]
    ax2.plot(z, d, color="#444444", lw=1.0)
    ax2.axhline(0, color="#999999", lw=0.6)
    ax2.set_ylabel("difference [G]")
    ax2.set_xlabel("z [cm]")
    ax2.grid(alpha=0.3)
    axr = ax2.twinx()
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(np.abs(ref["Bz"][i0]) > 100, 100 * d / ref["Bz"][i0], np.nan)
    axr.plot(z, rel, color="#e07b00", lw=0.8, alpha=0.8)
    axr.set_ylabel("relative [%] (|Bz| > 100 G)", color="#e07b00")
    fig.tight_layout()
    fn = stem + "_axis.png"
    fig.savefig(fn, dpi=90)
    plt.close(fig)
    files.append(fn)
    if not diffmap:
        return files
    step = max(1, int(round(len(z) / 600)))
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    ext = [z[0], z[-1], r[0], r[-1]]
    with np.errstate(divide="ignore", invalid="ignore"):
        relmap = 100 * dB / np.maximum(Bref, 1.0)
    lev_abs = [0.0, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 1e9]
    lev_rel = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 1e9]
    panels = [(dB, lev_abs, "|B - B_ref| [G]"), (relmap, lev_rel, "|B - B_ref| / max(|B_ref|, 1 G) [%]")]
    for ax, (img, lev, title) in zip(axs, panels):
        cmap = plt.get_cmap("viridis", len(lev) - 1)
        norm = BoundaryNorm(lev, cmap.N)
        im = ax.imshow(img[::step, ::step], origin="lower", extent=ext, aspect="auto", norm=norm, cmap=cmap,
                       interpolation="nearest")
        for reg in deck.regions[1:]:
            col = "white" if reg.mat >= 2 else "#ff5050"
            ax.add_patch(Polygon(reg.points[:, ::-1], closed=True, fill=False, ec=col, lw=0.5))
        ax.add_patch(Rectangle((t[2], t[0]), t[3] - t[2], t[1] - t[0], fill=False, ec="#ffd000", lw=0.8, ls="--"))
        ax.add_patch(Rectangle((p[2], p[0]), p[3] - p[2], p[1] - p[0], fill=False, ec="#ff7f0e", lw=0.8))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, spacing="uniform")
        cb.set_ticks(lev[1:-1])
        cb.ax.set_yticklabels(["%g" % v for v in lev[1:-1]], fontsize=7)
        ax.set_xlabel("z [cm]")
        ax.set_title(title, fontsize=10)
    axs[0].set_ylabel("r [cm]")
    fig.suptitle("%s: %s vs %s (dashed: tracking region, solid: PMT region)" % (label, sol_label, ref_label),
                 fontsize=10)
    fig.tight_layout()
    fn = stem + "_diffmap.png"
    fig.savefig(fn, dpi=80)
    plt.close(fig)
    files.append(fn)
    return files


if __name__ == "__main__":
    sys.exit(main())
