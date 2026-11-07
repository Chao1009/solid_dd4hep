#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Evaluate B(r, z) of a solve_rz.py solution at points given in cm; prints Br, Bz, |B| in gauss.

Example:  python3 eval_field.py work/CLEOv10_open.npz 0,0 0,-350 215,344
          python3 eval_field.py work/CLEOv10_open.npz --axis -500:700:10
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rzfield.driver import join_negative_values, parse_range  # noqa: E402
from rzfield.field import RZSolution  # noqa: E402
from rzfield.gemcmap import grid_axis  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("solution")
    ap.add_argument("points", nargs="*", metavar="R,Z", help="points in cm")
    ap.add_argument("--axis", help="also print Bz on the axis for z = min:max:step (cm)")
    a = ap.parse_args(join_negative_values(sys.argv[1:] if argv is None else argv, ("--axis",)))
    sol = RZSolution.load(a.solution)
    pts = [tuple(float(v) for v in p.split(",")) for p in a.points]
    if a.axis:
        pts += [(0.0, z) for z in grid_axis(*parse_range(a.axis))]
    if not pts:
        ap.error("no points given")
    r = np.array([p[0] for p in pts])
    z = np.array([p[1] for p in pts])
    Br, Bz = sol.B_cm_gauss(r, z)
    print("%10s %10s %14s %14s %14s" % ("R[cm]", "Z[cm]", "Br[G]", "Bz[G]", "|B|[G]"))
    for i in range(len(r)):
        print("%10.3f %10.3f %14.6g %14.6g %14.6g" % (r[i], z[i], Br[i], Bz[i], np.hypot(Br[i], Bz[i])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
