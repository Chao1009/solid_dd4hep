#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Solve a Poisson-style r-z magnet deck (.am) with nonlinear iron and write the solution and a GEMC map.

Examples:
  # validation setup: Dirichlet box of the deck (region 1), POISSON's 1010 table
  python3 solve_rz.py CLEOv10.am --mode poisson -o work/CLEOv10_poisson.npz
  # production: open boundary (half disk 10x the magnet), 1 cm map over the deck's region-1 box
  python3 solve_rz.py CLEOv10.am --mode open -o work/CLEOv10_open.npz --map work/solenoid_CLEOv10_rz.dat
"""
import argparse
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rzfield.driver import join_negative_values, parse_range, solve_deck, standard_points  # noqa: E402
from rzfield.gemcmap import grid_axis, map_from_solution, write_map  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck", help="Poisson Automesh input (.am)")
    ap.add_argument("--mode", choices=("poisson", "open"), default="open",
                    help="'poisson': A = 0 on the deck's region-1 box (validation); 'open': large half disk "
                         "(default)")
    ap.add_argument("--open-factor", type=float, default=10.0,
                    help="open mode: outer radius / magnet radius (default 10)")
    ap.add_argument("--mesh-scale", type=float, default=1.0, help="multiply all mesh sizes (0.5 = refine)")
    ap.add_argument("--h-boundary", type=float, help="mesh size at iron/coil boundaries, cm (default 1.0)")
    ap.add_argument("--h-box", type=float, help="max mesh size in the map box, cm (default 4.0)")
    ap.add_argument("--iron", default=None,
                    help="material for all MAT>=2 regions: poisson1010 (POISSON default 1010 table), "
                         "mu=<mu_r>, gamma=<1/mu_r> or a B-H table file (default: as the deck says)")
    ap.add_argument("--material", action="append", default=[], metavar="MAT=SPEC",
                    help="material for one MAT number, e.g. 3=mu=1000 or 2=my_steel.txt (repeatable)")
    ap.add_argument("--bh-interp", choices=("poisson", "pchip"), default="poisson",
                    help="B-H interpolation: gamma linear in B as POISSON (default) or monotone cubic H(B)")
    ap.add_argument("--rtol", type=float, default=1e-9, help="Newton: ||R||/||F|| tolerance (default 1e-9)")
    ap.add_argument("--atol", type=float, default=1e-9, help="Newton: max|dA|/max|A| tolerance (default 1e-9)")
    ap.add_argument("--maxit", type=int, default=60)
    ap.add_argument("--linear-solver", choices=("lu", "amg"), default="lu",
                    help="Newton systems: sparse LU (SuperLU, default) or AMG-preconditioned CG (needs pyamg; "
                         "for meshes > ~1e6 nodes)")
    ap.add_argument("--threads", type=int, default=1,
                    help="gmsh threads (default 1: identical meshes run to run with the same gmsh version)")
    ap.add_argument("-o", "--output", help="solution file (.npz), reusable by write_map.py / eval_field.py")
    ap.add_argument("--map", help="write a GEMC map (cm, gauss, 'R Z Br Bz')")
    ap.add_argument("--map-name", help="name in the map header (default: map file stem)")
    ap.add_argument("--map-comment", help="comment attribute in the map header (default: the name)")
    ap.add_argument("--map-r", default=None, help="R range min:max:step in cm (default: region-1 box, 1 cm)")
    ap.add_argument("--map-z", default=None, help="Z range min:max:step in cm (default: region-1 box, 1 cm)")
    ap.add_argument("--map-air-side", action="store_true",
                    help="map: extrapolate the air field to the iron nodes next to air (write_map.py --air-side)")
    ap.add_argument("--point", action="append", default=[], metavar="R,Z",
                    help="print B at this point (cm); repeatable")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args(join_negative_values(sys.argv[1:] if argv is None else argv, ("--map-r", "--map-z", "--point")))

    log = (lambda *x, **k: None) if a.quiet else print
    overrides = {}
    if a.iron:
        overrides["all"] = a.iron
    for m in a.material:
        k, v = m.split("=", 1)
        overrides[int(k)] = v
    sizes = {}
    if a.h_boundary:
        sizes["h_boundary"] = a.h_boundary
    if a.h_box:
        sizes["h_box"] = a.h_box
    from rzfield.amdeck import read_deck
    deck = read_deck(a.deck)
    box = list(deck.bbox("box"))
    if a.map_r:
        r0, r1, _ = parse_range(a.map_r)
        box[0], box[1] = r0, r1
    if a.map_z:
        z0, z1, _ = parse_range(a.map_z)
        box[2], box[3] = z0, z1
    sol, deck = solve_deck(a.deck, mode=a.mode, mesh_scale=a.mesh_scale, open_factor=a.open_factor,
                           overrides=overrides, interp=a.bh_interp, refine_box=box, sizes=sizes or None,
                           rtol=a.rtol, atol=a.atol, maxit=a.maxit, threads=a.threads,
                           linear_solver=a.linear_solver, log=log)
    if a.output:
        sol.save(a.output)
        log("solution written to %s" % a.output)
    pts = [(0.0, 0.0)] + [tuple(float(v) for v in p.split(",")) for p in a.point]
    for r, z, br, bz in standard_points(sol, pts):
        print("B(r=%g cm, z=%g cm): Br %.4f G  Bz %.4f G  |B| %.4f G" % (r, z, br, bz, (br * br + bz * bz) ** 0.5))
    if a.map:
        rr = parse_range(a.map_r) if a.map_r else (box[0], box[1], 1.0)
        zz = parse_range(a.map_z) if a.map_z else (box[2], box[3], 1.0)
        r = grid_axis(*rr)
        z = grid_axis(*zz)
        info = {}
        Br, Bz = map_from_solution(sol, r, z, air_side=a.map_air_side, info=info)
        name = a.map_name or os.path.splitext(os.path.basename(a.map))[0]
        comment = a.map_comment
        if a.map_air_side:
            comment = (comment or name) + "; air-side values at the iron nodes next to air"
        write_map(a.map, name, r, z, Br, Bz, comment=comment)
        log("GEMC map %s: R %g..%g, Z %g..%g cm, %d x %d points%s" % (
            a.map, r[0], r[-1], z[0], z[-1], len(r), len(z),
            "; air-side values at %d iron nodes" % info["n_air_side_nodes"] if a.map_air_side else ""))
    return 0 if sol.info.get("converged") else 2


if __name__ == "__main__":
    sys.exit(main())
