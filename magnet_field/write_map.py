#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Write a GEMC 2D map (cm, gauss, rows 'R Z Br Bz', R-major with Z fastest) from a solve_rz.py solution.

Needs numpy and scipy only (no gmsh), so it also runs in the eic container.
Example:  python3 write_map.py work/CLEOv10_open.npz -o solenoid_CLEOv10_rz.dat --r 0:500:1 --z -500:700:1

Nodes exactly on an iron face get the air-side field. --air-side also replaces the values at the iron nodes of
the grid cells that contain air by an extrapolation of the air field: the interpolated map is then accurate in air
up to the iron faces, but wrong in the iron within about two cells of them (README.md, Limitations).
"""
import argparse
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rzfield.driver import join_negative_values, parse_range  # noqa: E402
from rzfield.field import RZSolution  # noqa: E402
from rzfield.gemcmap import grid_axis, map_from_solution, write_map  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("solution", help="solution .npz from solve_rz.py")
    ap.add_argument("-o", "--output", required=True, help="map file to write")
    ap.add_argument("--name", help="name in the header (default: output file stem)")
    ap.add_argument("--comment", help="comment attribute in the header (default: the name)")
    ap.add_argument("--r", default=None, help="R range min:max:step, cm (default: refinement box, 1 cm)")
    ap.add_argument("--z", default=None, help="Z range min:max:step, cm (default: refinement box, 1 cm)")
    ap.add_argument("--air-side", action="store_true",
                    help="extrapolate the air field to the iron nodes next to air (accurate in air up to the iron "
                         "faces, wrong in the iron within about 2 cells of them)")
    a = ap.parse_args(join_negative_values(sys.argv[1:] if argv is None else argv, ("--r", "--z")))
    sol = RZSolution.load(a.solution)
    box = sol.info.get("refine_box_cm")
    rr = parse_range(a.r) if a.r else (box[0], box[1], 1.0)
    zz = parse_range(a.z) if a.z else (box[2], box[3], 1.0)
    r = grid_axis(*rr)
    z = grid_axis(*zz)
    info = {}
    Br, Bz = map_from_solution(sol, r, z, air_side=a.air_side, info=info)
    name = a.name or os.path.splitext(os.path.basename(a.output))[0]
    comment = a.comment
    if a.air_side:
        comment = (comment or name) + "; air-side values at the iron nodes next to air"
    write_map(a.output, name, r, z, Br, Bz, comment=comment)
    print("wrote %s: R %g..%g step %g, Z %g..%g step %g cm (%d x %d points)%s"
          % (a.output, r[0], r[-1], rr[2], z[0], z[-1], zz[2], len(r), len(z),
             "; air-side values at %d iron nodes" % info["n_air_side_nodes"] if a.air_side else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
