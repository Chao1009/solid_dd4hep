#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Field map of an air-core coil set (models/*.json), e.g. the polarized NH3 target magnet.

  coil_map.py models/nh3_target_oxford2012.json -o fieldmaps/ptarget_oxford2012_6m.dat      # the map
  coil_map.py models/nh3_target_oxford2012.json --compare oxford_ptarget.dat                 # benchmark
  coil_map.py models/nh3_target_oxford2012.json --at 0:0 --at 10:2                           # field at points

The field is the Biot-Savart sum over the coils (rzfield/coils.py): no iron, no mesh, exact up to the quadrature
(about 1e-5 of the field or better, also inside the windings). The map is written in the GEMC format that
SoLID_FieldMapB reads: rows 'R Z Br Bz' in cm and gauss, R the distance from the magnet axis and Z the coordinate
along it, with the symmetry of the model's "axis" in the header (cylindrical-x for a magnet on the hall x axis).

  --r MIN:MAX:STEP, --z MIN:MAX:STEP   grid in cm (default 0:600:1 and -600:600:1)
  --compare REF.dat                    compare with another map on that map's grid: differences by distance from
                                       the windings, and the nodes where the radial component has the other sign
  --at R:Z                             print the field at a point (cm); repeatable
  --json FILE                          write the numbers of --compare and the reference values as JSON

Needs numpy and scipy. The whole default map (722 000 points) takes about 20 s.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rzfield.coils import _distance, axis_field, coil_field, in_winding, read_coils  # noqa: E402
from rzfield.driver import join_negative_values, parse_range  # noqa: E402
from rzfield.gemcmap import grid_axis, read_map, write_map  # noqa: E402

CLASSES = ((5.0, None), (2.0, 5.0), (1.0, 2.0), (0.5, 1.0), (0.0, 0.5))      # cm from the nearest winding


def reference_values(cs):
    """numbers that characterise the field, for the log and for the field fragment"""
    z = np.array([0.0, 100.0, 200.0, 300.0, 600.0])
    r = np.array([100.0, 200.0, 300.0, 600.0])
    bz_axis = axis_field(cs, z)
    _, bz_mid = coil_field(cs, r, np.zeros_like(r))
    # homogeneity over the target cell: +-1.5 cm along the axis and 1.5 cm radius
    rr, zz = np.meshgrid(np.linspace(0.0, 1.5, 7), np.linspace(-1.5, 1.5, 13), indexing="ij")
    br_c, bz_c = coil_field(cs, rr, zz)
    b0 = float(bz_axis[0])
    out = {"B_centre_G": b0, "ampere_turns": cs.ampere_turns,
           "homogeneity_r1.5_z1.5_cm": float(np.abs(np.hypot(br_c, bz_c) - b0).max() / abs(b0)),
           "Bz_axis_G": {"%g" % a: float(b) for a, b in zip(z[1:], bz_axis[1:])},
           "Bz_midplane_G": {"%g" % a: float(b) for a, b in zip(r, bz_mid)}}
    return out


def compare(cs, ref):
    """differences to a reference map on its grid"""
    R, Z = np.meshgrid(ref["r"], ref["z"], indexing="ij")
    br, bz = coil_field(cs, R, Z)
    ins = in_winding(cs, R, Z)
    dist = np.min([_distance(c, R, Z) for c in cs.coils], axis=0)
    dr_abs = np.abs(np.abs(br) - np.abs(ref["Br"]))
    dz = np.abs(bz - ref["Bz"])
    rows = []
    for lo, hi in CLASSES:
        sel = ~ins & (dist > lo) & ((dist <= hi) if hi is not None else True)
        if sel.any():
            rows.append({"distance_cm": [lo, hi], "nodes": int(sel.sum()), "max_abs_dBr_magnitude_G": float(dr_abs[sel].max()),
                         "max_abs_dBz_G": float(dz[sel].max()), "max_B_G": float(np.hypot(br, bz)[sel].max())})
    if ins.any():
        rows.append({"distance_cm": "inside a winding", "nodes": int(ins.sum()),
                     "max_abs_dBr_magnitude_G": float(dr_abs[ins].max()), "max_abs_dBz_G": float(dz[ins].max()),
                     "max_B_G": float(np.hypot(br, bz)[ins].max())})
    # sign of the radial component, where both are clearly non-zero
    eps = 0.06
    sign = (np.sign(br) != np.sign(ref["Br"])) & (np.abs(br) > eps) & (np.abs(ref["Br"]) > eps)
    out = {"grid": {"r_cm": [float(ref["r"][0]), float(ref["r"][-1])], "z_cm": [float(ref["z"][0]), float(ref["z"][-1])],
                    "nodes": int(R.size)},
           "by_distance_from_windings": rows,
           "Br_sign_differs": {"nodes": int(sign.sum()), "outside_windings": int((sign & ~ins).sum())}}
    if sign.any():
        out["Br_sign_differs"].update({"r_cm": [float(R[sign].min()), float(R[sign].max())],
                                       "abs_z_cm": [float(np.abs(Z[sign]).min()), float(np.abs(Z[sign]).max())],
                                       "largest_abs_Br_G": float(np.abs(br[sign]).max())})
        # the reference against itself: div B = 0 on the axis gives Br(r, z) = -(r/2) dBz/dz to first order in r
        j = np.nonzero((ref["z"] > 0) & (ref["z"] < ref["z"][-1]))[0]
        slope = (ref["Bz"][0, j + 1] - ref["Bz"][0, j - 1]) / (ref["z"][j + 1] - ref["z"][j - 1])
        pred = -0.5 * ref["r"][1] * slope
        clear = np.abs(pred) > 20.0 * eps
        out["Br_sign_differs"]["reference_against_its_own_axis_field"] = {
            "nodes_tested": int(clear.sum()),
            "reference_Br_has_other_sign_than_its_axis_field_requires": int((np.sign(ref["Br"][1, j][clear]) != np.sign(pred[clear])).sum()),
            "computed_Br_has_other_sign": int((np.sign(br[1, j][clear]) != np.sign(pred[clear])).sum())}
    return out


def main(argv=None):
    argv = join_negative_values(sys.argv[1:] if argv is None else argv, ("--r", "--z", "--at"))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="coil set (JSON)")
    ap.add_argument("-o", "--output", help="map file to write")
    ap.add_argument("--r", default="0:600:1", help="radial grid min:max:step in cm (default 0:600:1)")
    ap.add_argument("--z", default="-600:600:1", help="grid along the magnet axis min:max:step in cm (default -600:600:1)")
    ap.add_argument("--name", help="name in the header (default: output file stem)")
    ap.add_argument("--compare", metavar="REF.dat", help="GEMC map to compare with, on its grid")
    ap.add_argument("--at", action="append", default=[], metavar="R:Z", help="print the field at this point (cm)")
    ap.add_argument("--json", help="write reference values and comparison as JSON")
    args = ap.parse_args(argv)

    cs = read_coils(args.model)
    report = {"model": cs.name, "source": cs.source, "coils": len(cs.coils), "reference": reference_values(cs)}
    ref = report["reference"]
    print("%s: %d coils, %.4g ampere-turns; B at the centre %.3f G, homogeneity over the target cell %.1e"
          % (cs.name, len(cs.coils), cs.ampere_turns, ref["B_centre_G"], ref["homogeneity_r1.5_z1.5_cm"]))
    print("  on the axis:      " + ", ".join("%s cm: %.2f G" % kv for kv in ref["Bz_axis_G"].items()))
    print("  in the mid plane: " + ", ".join("%s cm: %.2f G" % kv for kv in ref["Bz_midplane_G"].items()))
    for item in args.at:
        r_, z_ = (float(v) for v in item.split(":"))
        br, bz = coil_field(cs, np.array([r_]), np.array([z_]))
        print("  r %g z %g cm: Br %.4f G, Bz %.4f G, |B| %.4f G" % (r_, z_, br[0], bz[0], np.hypot(br[0], bz[0])))

    if args.compare:
        cmp_ = compare(cs, read_map(args.compare))
        report["comparison"] = dict(cmp_, reference_map=os.path.basename(args.compare))
        print("comparison with %s (%d nodes): largest differences by distance from the windings"
              % (os.path.basename(args.compare), cmp_["grid"]["nodes"]))
        for row in cmp_["by_distance_from_windings"]:
            d = row["distance_cm"]
            lab = d if isinstance(d, str) else ("more than %g cm" % d[0] if d[1] is None else "%g to %g cm" % tuple(d))
            print("  %-18s %6d nodes: | |Br| - |ref Br| | <= %9.3f G, |Bz - ref Bz| <= %9.3f G (field up to %.0f G)"
                  % (lab, row["nodes"], row["max_abs_dBr_magnitude_G"], row["max_abs_dBz_G"], row["max_B_G"]))
        s = cmp_["Br_sign_differs"]
        if s["nodes"]:
            t = s["reference_against_its_own_axis_field"]
            print("  the radial component has the other sign at %d nodes (%d outside the windings): r %g to %g cm, |z| %g to "
                  "%g cm, up to %.0f G" % (s["nodes"], s["outside_windings"], s["r_cm"][0], s["r_cm"][1],
                                           s["abs_z_cm"][0], s["abs_z_cm"][1], s["largest_abs_Br_G"]))
            print("  div B = 0 next to the axis, from the reference's own axis field: the reference has the wrong sign "
                  "at %d of %d nodes, this calculation at %d" % (t["reference_Br_has_other_sign_than_its_axis_field_requires"],
                                                              t["nodes_tested"], t["computed_Br_has_other_sign"]))
        else:
            print("  the radial component has the same sign everywhere")

    if args.output:
        rmin, rmax, dr = parse_range(args.r)
        zmin, zmax, dz = parse_range(args.z)
        r, z = grid_axis(rmin, rmax, dr), grid_axis(zmin, zmax, dz)
        if rmin < 0:
            ap.error("--r starts at 0 or above")
        t0 = time.time()
        R, Z = np.meshgrid(r, z, indexing="ij")
        br, bz = coil_field(cs, R, Z)
        br[r == 0.0, :] = 0.0
        # exact zeros of the symmetry (Br in the mid plane) come out as rounding noise of about 1e-13 G, which differs
        # from one numpy or scipy version to the next: write them as zero, so that the file is reproducible
        br[np.abs(br) < 1e-9] = 0.0
        bz[np.abs(bz) < 1e-9] = 0.0
        name = args.name or os.path.splitext(os.path.basename(args.output))[0]
        comment = ("air-core field of %s (%s), Biot-Savart sum by solid_dd4hep magnet_field/coil_map.py; R = distance "
                   "from the magnet axis, Z = coordinate along it" % (cs.name, cs.source or "no source given"))
        write_map(args.output, name, r, z, br, bz, comment=comment, symmetry="cylindrical-%s" % cs.axis)
        edge = max(float(np.hypot(br[-1, :], bz[-1, :]).max()), float(np.hypot(br[:, 0], bz[:, 0]).max()),
                   float(np.hypot(br[:, -1], bz[:, -1]).max()))
        report["map"] = {"file": os.path.basename(args.output), "r_cm": [rmin, rmax, dr], "z_cm": [zmin, zmax, dz],
                         "points": int(R.size), "largest_field_on_the_edge_G": edge, "max_B_G": float(np.hypot(br, bz).max())}
        print("wrote %s: %d x %d points in %.0f s, largest field on the edge of the map %.2f G"
              % (args.output, len(r), len(z), time.time() - t0, edge))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(report, f, indent=1)
            f.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
