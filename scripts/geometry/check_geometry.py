#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Load a SoLID compact file and run the ROOT/TGeo overlap checker.

usage: check_geometry.py <compact.xml> [--tolerance CM] [--sampling N] [--max-print N] [--only REGEX]

Exit status: 0 if the geometry loads and no overlap/extrusion above the tolerance is found,
             1 if overlaps are found, 2 if loading fails.
Needs the plugin library on LD_LIBRARY_PATH and DETECTOR_PATH set (source thissolid.sh).

Notes:
  * TGeo's CheckOverlaps samples volume surfaces; with --sampling > 0 it additionally shoots
    random points (slower, finds more).  Geant4's checker (checkOverlaps / ddsim --runType=...)
    remains the reference before productions.
  * Tolerance is in cm (TGeo units).  Default 1 um.
"""
import argparse
import re
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('compact')
    ap.add_argument('--tolerance', type=float, default=1e-4, help='overlap tolerance in cm (default 1e-4 = 1 um)')
    ap.add_argument('--sampling', type=int, default=0, help='random sampling points per volume (0 = surface check only)')
    ap.add_argument('--max-print', type=int, default=50)
    ap.add_argument('--only', default=None, help='only report overlaps whose title matches REGEX')
    a = ap.parse_args()

    import ROOT  # noqa: F401
    import dd4hep
    t0 = time.time()
    desc = dd4hep.Detector.getInstance()
    try:
        desc.fromXML(a.compact)
    except Exception as e:  # noqa: BLE001
        print(f'LOAD FAILED: {e}')
        return 2
    t1 = time.time()
    world = desc.world()
    dets = sorted((name, de.id()) for name, de in world.children())
    print(f'loaded {a.compact} in {t1 - t0:.1f} s: {len(dets)} top-level detectors')
    for name, did in dets:
        print(f'  {did:4d}  {name}')

    mgr = desc.manager()
    mgr.SetVerboseLevel(0)
    if a.sampling > 0:
        mgr.CheckOverlaps(a.tolerance, 's%d' % a.sampling)
    else:
        mgr.CheckOverlaps(a.tolerance)
    ovl = mgr.GetListOfOverlaps()
    items = []
    for i in range(ovl.GetEntries()):
        o = ovl.At(i)
        title = o.GetTitle()
        if a.only and not re.search(a.only, title):
            continue
        items.append((o.GetOverlap(), title))
    items.sort(reverse=True)
    print(f'overlap check ({time.time() - t1:.1f} s, tolerance {a.tolerance} cm): {len(items)} overlaps/extrusions')
    for size, title in items[:a.max_print]:
        print(f'  OVERLAP {size * 10:.4f} mm  {title}')
    if len(items) > a.max_print:
        print(f'  ... {len(items) - a.max_print} more')
    return 1 if items else 0


if __name__ == '__main__':
    sys.exit(main())
