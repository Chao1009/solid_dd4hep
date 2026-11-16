#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests of the air-core coil field (rzfield/coils.py) and of the target magnet model (ctest: coil_field).

usage: test_coils.py [--reference MAP.dat]

  axis        on the axis the sum over loops equals the closed form for a coil of rectangular cross-section
  maxwell     div B = 0 everywhere and curl B = mu0 J: zero outside the windings, the current density inside
  near        next to a winding and inside it, the default rule agrees with a rule of four times the order
  symmetry    a split pair: Bz even and Br odd in z, Br = 0 in the mid plane and on the axis
  far         far away the field falls like a dipole
  model       the Oxford 2012 NH3 target magnet: 5 T at the centre, homogeneous over the target cell
  map         a map written with the cylindrical-x header reads back; writing it twice gives the same bytes
  reference   with --reference oxford_ptarget.dat (solid_gemc/field): agreement within that file's rounding more
              than 5 cm from the windings, and the nodes where its radial component has the wrong sign
Needs numpy and scipy. Exit status 0 if every check passes.
"""
import argparse
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from rzfield.coils import MU0_G_CM_PER_A, Coil, axis_field, coil_field, in_winding, read_coils  # noqa: E402
from rzfield.gemcmap import read_map, write_map  # noqa: E402

FAILED = []


def check(ok, what):
    print("%s  %s" % ("ok  " if ok else "FAIL", what))
    if not ok:
        FAILED.append(what)


def derivatives(coils, r, z, h=1e-3):
    """div B and (curl B)_phi by central differences, G/cm"""
    def f(rr, zz):
        a, b = coil_field(coils, np.array([rr]), np.array([zz]))
        return a[0], b[0]
    brp, bzp = f(r + h, z)
    brm, bzm = f(r - h, z)
    bru, bzu = f(r, z + h)
    brd, bzd = f(r, z - h)
    div = ((r + h) * brp - (r - h) * brm) / (2 * h * r) + (bzu - bzd) / (2 * h)
    curl = (bru - brd) / (2 * h) - (bzp - bzm) / (2 * h)
    return div, curl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", help="solid_gemc/field/oxford_ptarget.dat, for the comparison")
    args = ap.parse_args()

    one = [Coil("test", 10.0, 14.0, -3.0, 5.0, 2.0e5)]
    pair = [Coil("up", 10.0, 14.0, 6.0, 11.0, 2.0e5), Coil("down", 10.0, 14.0, -11.0, -6.0, 2.0e5)]

    # axis
    z = np.array([-40.0, -3.0, 0.0, 1.0, 5.0, 8.0, 30.0, 300.0])
    _, bz = coil_field(one, np.zeros_like(z), z)
    ref = axis_field(one, z)
    check(np.abs(bz / ref - 1.0).max() < 1e-9, "on the axis the loop sum equals the closed form (largest deviation %.1e)"
          % np.abs(bz / ref - 1.0).max())

    # Maxwell
    worst_div, worst_curl, scale = 0.0, 0.0, 0.0
    for r_, z_ in ((3.0, 1.0), (9.0, 6.0), (16.0, -2.0), (12.0, 9.0), (30.0, 20.0), (9.9, 1.0), (14.2, 5.1)):
        d, c = derivatives(one, r_, z_)
        a, b = coil_field(one, np.array([r_]), np.array([z_]))
        s = float(np.hypot(a[0], b[0]))
        worst_div, worst_curl = max(worst_div, abs(d) / s), max(worst_curl, abs(c) / s)
    check(worst_div < 1e-5 and worst_curl < 1e-5,
          "outside the winding div B = 0 and curl B = 0 (largest %.1e and %.1e of |B| per cm)" % (worst_div, worst_curl))
    d, c = derivatives(one, 12.0, 1.0)
    j = one[0].current_density
    check(abs(c / (MU0_G_CM_PER_A * j) - 1.0) < 1e-3 and abs(d) < 1e-3 * MU0_G_CM_PER_A * j,
          "inside the winding curl B = mu0 J (ratio %.5f) and div B = 0" % (c / (MU0_G_CM_PER_A * j)))

    # near field against a much higher order
    pts = np.array([(10.0, 1.0), (12.0, 1.0), (13.9, 4.9), (14.05, 2.0), (9.5, -3.2), (12.0, 5.3), (14.0, 5.0)])
    a1, b1 = coil_field(one, pts[:, 0], pts[:, 1])
    a2, b2 = coil_field(one, pts[:, 0], pts[:, 1], orders=(), near_order=96)
    dev = max(np.abs(a1 - a2).max(), np.abs(b1 - b2).max()) / np.hypot(a2, b2).max()
    check(dev < 2e-5, "at and inside the winding the default rule agrees with order 96 (largest deviation %.1e)" % dev)

    # symmetry of a split pair
    r = np.array([2.0, 8.0, 20.0, 50.0])
    brp, bzp = coil_field(pair, r, np.full(4, 4.0))
    brm, bzm = coil_field(pair, r, np.full(4, -4.0))
    br0, _ = coil_field(pair, r, np.zeros(4))
    bra, _ = coil_field(pair, np.zeros(3), np.array([1.0, 7.0, 30.0]))
    check(np.allclose(bzp, bzm, rtol=1e-12) and np.allclose(brp, -brm, rtol=1e-12) and np.abs(br0).max() < 1e-9
          and np.abs(bra).max() == 0.0, "split pair: Bz even and Br odd in z, Br = 0 in the mid plane and on the axis")

    # dipole far field: B_z(axis) = 2 m / d^3, B_z(mid plane) = -m / d^3
    d1, d2 = 2000.0, 4000.0
    ax = axis_field(pair, np.array([d1, d2]))
    _, mid = coil_field(pair, np.array([d1, d2]), np.zeros(2))
    check(abs(ax[0] / ax[1] / 8.0 - 1.0) < 1e-3 and abs(mid[0] / mid[1] / 8.0 - 1.0) < 1e-3 and abs(ax[1] / mid[1] + 2.0) < 1e-3,
          "far away the field is a dipole: 1/d^3, and twice as strong on the axis as in the mid plane with the other sign")

    # the target magnet model
    model = os.path.join(os.path.dirname(HERE), "models", "nh3_target_oxford2012.json")
    cs = read_coils(model)
    rr, zz = np.meshgrid(np.linspace(0.0, 1.5, 7), np.linspace(-1.5, 1.5, 13), indexing="ij")
    br, bz = coil_field(cs, rr, zz)
    b0 = float(axis_field(cs, np.array([0.0]))[0])
    hom = float(np.abs(np.hypot(br, bz) - b0).max() / b0)
    check(len(cs.coils) == 14 and abs(b0 - 49996.254) < 0.01 and hom < 1e-4 and cs.axis == "x",
          "Oxford 2012 NH3 target magnet: 14 coils, %.3f G at the centre, homogeneous to %.1e over the target cell" % (b0, hom))
    _, bz2 = coil_field(cs, np.array([0.0, 200.0]), np.array([200.0, 0.0]))
    check(abs(bz2[0] - 87.724) < 0.01 and abs(bz2[1] + 42.98) < 0.01,
          "2 m from the target: %.2f G on the magnet axis, %.2f G in the mid plane" % (bz2[0], bz2[1]))
    br, _ = coil_field(cs, np.array([10.0]), np.array([2.0]))
    check(br[0] < -3000.0, "in the gap between the coil halves the radial field points to the axis (r 10, z 2 cm: %.0f G)" % br[0])

    # map file
    with tempfile.TemporaryDirectory() as tmp:
        r, z = np.arange(0.0, 31.0, 1.0), np.arange(-20.0, 21.0, 1.0)
        R, Z = np.meshgrid(r, z, indexing="ij")
        br, bz = coil_field(cs, R, Z)
        f1, f2 = os.path.join(tmp, "a.dat"), os.path.join(tmp, "b.dat")
        write_map(f1, "test", r, z, br, bz, comment="test", symmetry="cylindrical-x")
        write_map(f2, "test", r, z, br, bz, comment="test", symmetry="cylindrical-x")
        m = read_map(f1)
        text = open(f1).read()
        check('type="cylindrical-x"' in text and np.allclose(m["Br"], br, rtol=1e-6, atol=1e-6) and np.allclose(m["Bz"], bz, rtol=1e-6)
              and open(f2).read() == text, "a cylindrical-x map reads back, and writing it twice gives the same file")

    if args.reference:
        ref = read_map(args.reference)
        R, Z = np.meshgrid(ref["r"], ref["z"], indexing="ij")
        br, bz = coil_field(cs, R, Z)
        ins = in_winding(cs, R, Z)
        from rzfield.coils import _distance
        dist = np.min([_distance(c, R, Z) for c in cs.coils], axis=0)
        far = ~ins & (dist > 5.0)
        d1 = float(np.abs(np.abs(br) - np.abs(ref["Br"]))[far].max())
        d2 = float(np.abs(bz - ref["Bz"])[far].max())
        check(d1 < 0.07 and d2 < 0.07, "more than 5 cm from the windings (%d nodes) the reference map agrees within its 0.1 G "
              "rounding: |Br| to %.3f G, Bz to %.3f G" % (far.sum(), d1, d2))
        sign = (np.sign(br) != np.sign(ref["Br"])) & (np.abs(br) > 0.06) & (np.abs(ref["Br"]) > 0.06)
        check((ref["Br"][Z > 0] >= 0).all() and (br[Z > 0] < -0.06).sum() > 300,
              "the reference map never has Br < 0 for z > 0; the field does, at %d of its nodes (the sign differs at %d "
              "nodes in all)" % ((br[Z > 0] < -0.06).sum(), sign.sum()))
    print("%d check%s failed" % (len(FAILED), "" if len(FAILED) == 1 else "s") if FAILED else "all checks passed")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
