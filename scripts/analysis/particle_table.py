#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Make the particle table of a simulated file: one row per generated particle with what it left in every detector.

usage: particle_table.py SIM.edm4hep.root [...] -o TABLE.npz

Several input files are joined into one table (same configuration). The table is described in
python/solid_ana/table.py; scripts/analysis/acceptance.py reads it. Needs uproot, awkward, numpy.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))

from solid_ana.table import build_table, concat_tables, save_table  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()
    tables, infos = [], []
    for f in a.files:
        t, info = build_table(f)
        tables.append(t)
        infos.append(info)
        print("%s: %d events, %d generated particles" % (f, info["events"], info["rows"]))
    arrays = tables[0] if len(tables) == 1 else concat_tables(tables)
    info = dict(infos[0], file=[i["file"] for i in infos], events=sum(i["events"] for i in infos),
                rows=int(len(arrays["pdg"])))
    save_table(a.output, arrays, info)
    for name, c in info["collections"].items():
        if c["kind"] == "tracker" and c.get("layers"):
            n = int(arrays[name + "_hit"].any(axis=1).sum())
            print("  %-22s %d layers, element field %s: %d particles with a hit" % (name, c["layers"], c["element_field"], n))
        elif c["kind"] == "calo":
            print("  %-22s %d particles with energy" % (name, int((arrays[name + "_e"] > 0).sum())))
    print("wrote %s: %d rows" % (a.output, info["rows"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
