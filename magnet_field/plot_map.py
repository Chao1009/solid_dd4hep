#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Draw a 2D (r-z) field map: |B| in colour bands, field lines, the magnet outline and the on-axis field.

  plot_map.py MAP.dat [--deck models/solid_v4c.am] [-o map.png]
              [--mark "3He target:-350:0" ...] [--zlim -600 800] [--rlim 0 500] [--lines 28] [--arrows 50]
              [--title TEXT]

MAP.dat is a GEMC map with <mfield> header (rows R Z Br Bz), as written by make_maps.sh or POISSON/SF7. The field
lines are contours of the flux function psi(r, z) = integral_0^r Bz r' dr' in equal steps, so the same flux runs
between neighbouring lines. --deck draws the iron and coil regions of the r-z model; give it for a map with iron:
psi is then carried through the air only (a 1 cm map cannot be integrated through iron), and the lines are drawn in
air. --mark LABEL:Z:R (cm) marks a
point and prints |B| there. Needs numpy and matplotlib only.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rzfield.gemcmap import read_map  # noqa: E402

C_SURF, C_INK, C_INK2, C_MUTED, C_GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
C_IRON, C_COIL, C_LINE = "#eb6834", "#1baf7a", "#2a78d6"
RAMP = ["#fcfcfb", "#e6f0fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
BANDS_G = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]      # |B| band edges (gauss)


def fmt_B(g):
    """gauss -> short text"""
    if g >= 1000:
        return "%.3g T" % (g / 1e4)
    return "%.3g G" % g


def bilinear(m, r, z):
    ri = np.interp(r, m["r"], np.arange(len(m["r"])))
    zi = np.interp(z, m["z"], np.arange(len(m["z"])))
    i, j = min(int(ri), len(m["r"]) - 2), min(int(zi), len(m["z"]) - 2)
    a, b = ri - i, zi - j
    out = []
    for F in (m["Br"], m["Bz"]):
        out.append(F[i, j] * (1 - a) * (1 - b) + F[i + 1, j] * a * (1 - b) + F[i, j + 1] * (1 - a) * b
                   + F[i + 1, j + 1] * a * b)
    return out


def flux_function(m, iron=None):
    """psi(r, z) = integral_0^r Bz r' dr' (G cm2).

    Without iron: trapezoid along r. With an iron mask (nR, nZ): the integral through iron on a 1 cm grid is poor
    (B jumps at the faces), so psi is carried from the axis through the air only, along r with r Bz and along z
    with -r Br, in alternating line sweeps; iron nodes and air pockets that cannot be reached stay NaN."""
    r = m["r"][:, None]
    fr = m["Bz"] * r                                                   # d psi / d r
    fz = -m["Br"] * r                                                  # d psi / d z
    dr = 0.5 * (fr[1:] + fr[:-1]) * np.diff(m["r"])[:, None]           # psi[i+1] - psi[i]
    dz = 0.5 * (fz[:, 1:] + fz[:, :-1]) * np.diff(m["z"])[None, :]     # psi[:, j+1] - psi[:, j]
    psi = np.zeros_like(fr)
    if iron is None:
        psi[1:] = np.cumsum(dr, axis=0)
        return psi
    air = ~iron
    psi[:] = np.nan
    i0 = int(np.argmin(np.abs(m["r"])))
    psi[i0] = np.where(air[i0], 0.0, np.nan)
    nr, nz = psi.shape
    for _ in range(12):
        before = int(np.isfinite(psi).sum())
        for i in range(1, nr):                                         # outwards
            new = air[i] & np.isnan(psi[i]) & np.isfinite(psi[i - 1])
            psi[i, new] = psi[i - 1, new] + dr[i - 1, new]
        for j in range(1, nz):                                         # downstream
            new = air[:, j] & np.isnan(psi[:, j]) & np.isfinite(psi[:, j - 1])
            psi[new, j] = psi[new, j - 1] + dz[new, j - 1]
        for j in range(nz - 2, -1, -1):                                # upstream
            new = air[:, j] & np.isnan(psi[:, j]) & np.isfinite(psi[:, j + 1])
            psi[new, j] = psi[new, j + 1] - dz[new, j]
        for i in range(nr - 2, -1, -1):                                # inwards
            new = air[i] & np.isnan(psi[i]) & np.isfinite(psi[i + 1])
            psi[i, new] = psi[i + 1, new] - dr[i, new]
        if int(np.isfinite(psi).sum()) == before:
            break
    return psi


def iron_mask(deck, m, grow=1):
    """nodes inside iron regions of the deck, grown by `grow` nodes (map values next to iron faces are poor)"""
    from rzfield.geometry import points_in_polygon
    R, Z = np.meshgrid(m["r"], m["z"], indexing="ij")
    mask = np.zeros(R.shape, dtype=bool)
    for reg in deck.regions[1:]:
        if reg.mat >= 2:
            mask |= points_in_polygon(R, Z, reg.points)
    for _ in range(grow):
        g = mask.copy()
        g[1:] |= mask[:-1]
        g[:-1] |= mask[1:]
        g[:, 1:] |= mask[:, :-1]
        g[:, :-1] |= mask[:, 1:]
        mask = g
    return mask


def deck_patches(deck):
    from matplotlib.patches import Polygon
    out = []
    for reg in deck.regions[1:]:
        iron = reg.mat >= 2
        if not iron and reg.cur == 0.0:
            continue
        pts = np.column_stack([reg.points[:, 1], reg.points[:, 0]])
        out.append(Polygon(pts, closed=True, facecolor="none", edgecolor=C_IRON if iron else C_COIL,
                           linewidth=0.9 if iron else 1.1, zorder=4))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map")
    ap.add_argument("--deck", help="r-z model (.am) whose iron and coil regions are outlined")
    ap.add_argument("-o", "--output", help="image file (default: <map name>.png)")
    ap.add_argument("--mark", action="append", default=[], metavar="LABEL:Z:R", help="mark a point (cm); repeatable")
    ap.add_argument("--zlim", nargs=2, type=float, metavar=("ZMIN", "ZMAX"))
    ap.add_argument("--rlim", nargs=2, type=float, metavar=("RMIN", "RMAX"))
    ap.add_argument("--lines", type=int, default=28, help="number of field lines (default 28; 0 = none)")
    ap.add_argument("--arrows", type=float, default=50.0, metavar="STEP",
                    help="spacing (cm) of the field-direction arrows, drawn in air where |B| >= 1 G (default 50; 0 = none)")
    ap.add_argument("--title")
    ap.add_argument("--dpi", type=int, default=130)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.lines import Line2D

    m = read_map(args.map)
    r, z = m["r"], m["z"]
    B = np.hypot(m["Br"], m["Bz"])
    zlim = args.zlim or (z[0], z[-1])
    rlim = args.rlim or (r[0], r[-1])
    marks = []
    for s in args.mark:
        lab, zz, rr = s.rsplit(":", 2)
        marks.append((lab, float(zz), float(rr)))

    wz, wr = zlim[1] - zlim[0], rlim[1] - rlim[0]
    W = 12.0
    h_map = (W * 0.86) * wr / wz
    fig = plt.figure(figsize=(W, h_map + 3.6), dpi=args.dpi, facecolor=C_SURF)
    H = h_map + 3.6
    ax = fig.add_axes([0.06, (2.75) / H, 0.86, h_map / H])
    cax = fig.add_axes([0.93, (2.75) / H, 0.012, h_map / H])
    ax2 = fig.add_axes([0.06, 0.55 / H, 0.86, 1.75 / H], sharex=ax)

    for a in (ax, ax2):
        a.set_facecolor(C_SURF)
        for s in a.spines.values():
            s.set_color(C_GRID)
            s.set_linewidth(0.8)
        a.tick_params(colors=C_INK2, labelsize=8, length=3, color=C_GRID)

    # |B| in bands
    cmap = ListedColormap(RAMP)
    cmap.set_under(RAMP[0])
    cmap.set_over(RAMP[-1])
    norm = BoundaryNorm(BANDS_G, cmap.N)
    im = ax.imshow(B, origin="lower", extent=(z[0], z[-1], r[0], r[-1]), aspect="equal", cmap=cmap, norm=norm,
                   interpolation="nearest", zorder=1)
    cb = fig.colorbar(im, cax=cax, extend="min", ticks=BANDS_G)
    cb.ax.set_yticklabels([fmt_B(v) for v in BANDS_G])
    cb.ax.tick_params(labelsize=7, colors=C_INK2, length=2, color=C_GRID)
    cb.outline.set_edgecolor(C_GRID)
    cb.set_label("|B|", color=C_INK2, fontsize=8)

    # field lines: equal flux between neighbours
    deck = None
    if args.deck:
        from rzfield.amdeck import read_deck
        deck = read_deck(args.deck)
    if args.lines > 0:
        psi = flux_function(m, iron_mask(deck, m) if deck is not None else None)
        sel = np.isfinite(psi) & (r[:, None] >= rlim[0]) & (r[:, None] <= rlim[1]) & (z[None, :] >= zlim[0]) & (z[None, :] <= zlim[1])
        lo, hi = psi[sel].min(), psi[sel].max()
        lev = lo + (hi - lo) * (np.arange(args.lines) + 0.5) / args.lines
        # white under ink: visible on the light and on the dark bands
        psi = np.ma.masked_invalid(psi)
        ax.contour(z, r, psi, levels=lev, colors="#ffffff", linewidths=1.5, alpha=0.55, zorder=2)
        ax.contour(z, r, psi, levels=lev, colors=C_INK, linewidths=0.55, alpha=0.8, zorder=3)

    # field direction on a coarse grid (unit arrows), where the field lines are too far apart to show it
    if args.arrows > 0:
        st = args.arrows
        za = np.arange(np.ceil(zlim[0] / st) * st + st / 2, zlim[1], st)
        ra = np.arange(np.ceil(rlim[0] / st) * st + st / 2, rlim[1], st)
        ii = np.clip(np.searchsorted(r, ra), 0, len(r) - 1)
        jj = np.clip(np.searchsorted(z, za), 0, len(z) - 1)
        br, bz = m["Br"][np.ix_(ii, jj)], m["Bz"][np.ix_(ii, jj)]
        bb = np.hypot(br, bz)
        ok = bb >= BANDS_G[0]
        if deck is not None:
            ok &= ~iron_mask(deck, m, grow=3)[np.ix_(ii, jj)]
        ZA, RA = np.meshgrid(z[jj], r[ii])
        bb = np.where(ok, bb, 1.0)
        ax.quiver(ZA[ok], RA[ok], (bz / bb)[ok], (br / bb)[ok], pivot="mid", angles="xy", scale_units="xy",
                  scale=1.0 / (0.42 * st), width=0.0011, headwidth=4.5, headlength=5, headaxislength=4.5,
                  color=C_INK2, zorder=3.5)

    handles = []
    if args.lines > 0:
        handles.append(Line2D([], [], color=C_INK, linewidth=0.8, label="field lines in air (equal flux between lines)"))
    if args.arrows > 0:
        handles.append(Line2D([], [], color=C_INK2, linewidth=0, marker=r"$\rightarrow$", markersize=9,
                              label="field direction"))
    if deck is not None:
        for p in deck_patches(deck):
            ax.add_patch(p)
        handles += [Line2D([], [], color=C_IRON, linewidth=1.2, label="iron (r-z model)"),
                    Line2D([], [], color=C_COIL, linewidth=1.2, label="coil")]

    for lab, zz, rr in marks:
        br, bz = bilinear(m, rr, zz)
        b = float(np.hypot(br, bz))
        ax.plot([zz], [rr], marker="o", ms=5, mfc=C_SURF, mec=C_INK, mew=1.1, zorder=6, clip_on=False)
        up = rr < rlim[0] + 0.75 * wr
        t = ax.annotate("%s\n%s" % (lab, fmt_B(b)), (zz, rr), xytext=(0, 9 if up else -9),
                        textcoords="offset points", ha="center", va="bottom" if up else "top", fontsize=7.5,
                        color=C_INK, zorder=7, linespacing=1.15)
        t.set_bbox(dict(boxstyle="round,pad=0.2", facecolor=C_SURF, edgecolor="none", alpha=0.85))
        print("%-24s z %8.2f r %8.2f cm   Br %11.3f  Bz %11.3f  |B| %11.3f G" % (lab, zz, rr, br, bz, b))

    ax.set_xlim(*zlim)
    ax.set_ylim(*rlim)
    ax.set_ylabel("r (cm)", color=C_INK2, fontsize=9)
    plt.setp(ax.get_xticklabels(), visible=False)
    if handles:
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.055, (2.75 + h_map + 0.36) / H),
                   ncol=len(handles), frameon=False, fontsize=8, labelcolor=C_INK2)

    # on-axis field
    j0 = int(np.argmin(np.abs(r)))
    bz = m["Bz"][j0]
    ax2.plot(z, np.abs(bz), color=C_LINE, linewidth=2, solid_capstyle="round")
    neg = bz < 0
    if neg.any():
        ax2.plot(z, np.where(neg, -bz, np.nan), color=C_IRON, linewidth=2)
        ax2.legend(handles=[Line2D([], [], color=C_LINE, linewidth=2, label="Bz > 0"),
                            Line2D([], [], color=C_IRON, linewidth=2, label="Bz < 0")],
                   loc="upper right", frameon=False, fontsize=8, labelcolor=C_INK2)
    ax2.set_yscale("log")
    inside = (z >= zlim[0]) & (z <= zlim[1])
    top = np.abs(bz[inside]).max()
    ax2.set_ylim(max(0.1, min(1.0, np.abs(bz[inside]).min() * 0.7)), top * 2.2)
    ax2.grid(True, which="major", color=C_GRID, linewidth=0.5)
    ax2.set_axisbelow(True)
    ax2.set_xlabel("z (cm)", color=C_INK2, fontsize=9)
    ax2.set_ylabel("|Bz| on the axis (G)", color=C_INK2, fontsize=9)
    k = int(np.argmin(np.abs(z)))
    integ = float(np.sum(0.5 * (bz[1:] + bz[:-1]) * np.diff(z))) * 1e-6          # G cm -> T m
    ax2.text(0.5, 0.10, "Bz(0, 0) = %s      integral of Bz dz on the axis = %.3f T m" % (fmt_B(bz[k]), integ),
             transform=ax2.transAxes, fontsize=8, color=C_INK2, va="bottom", ha="center")
    for lab, zz, rr in marks:
        if abs(rr) < 1e-9:
            b = abs(float(np.interp(zz, z, bz)))
            ax2.plot([zz], [b], marker="o", ms=5, mfc=C_SURF, mec=C_INK, mew=1.1, zorder=5)

    fig.text(0.06, (H - 0.22) / H, args.title or m["name"], color=C_INK, fontsize=11.5, ha="left", va="top")
    out = args.output or (m["name"] or os.path.splitext(os.path.basename(args.map))[0]) + ".png"
    fig.savefig(out, facecolor=C_SURF)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
