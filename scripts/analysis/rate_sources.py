#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Where the hit rates of a beam-on-target simulation come from: by place of origin and by kind of particle.

usage:
    rate_sources.py SIM.edm4hep.root [...] --electrons-per-event N --current UA [--emin NAME=GEV ...]
                    [--target-z -350] [-o SOURCES.json]

For every hit above threshold (as in rates.py) the particle that made it is followed back through its stored
ancestors to the beam electron. Three things are read off the chain:

    origin     where the first secondary of the chain was produced, i.e. where the beam electron did what led to
               the hit: the entrance window, the target cell's windows, wall and gas, the collimators, the exit
               window, the exit beam pipe, the beam stop, or material away from the beam line. The regions are
               those of the polarized 3He target (compact/target/target_He3.xml, beamline_SIDIS_He3.xml), placed
               by --target-z.
    first      what that first secondary was: a photon (bremsstrahlung), an electron (Moller scattering and
               knock-on electrons), a positron, a hadron (photo- and electro-nuclear reactions)
    at detector  what made the hit and how it got there: an electron or positron that came straight from the
               beam line, one that a photon made in or near the detector, a hadron

ddsim does not store every secondary, unless it is run with --part.keepAllParticles True (larger files). A hit
of a particle that was not stored is booked on its nearest stored ancestor; when that is the beam electron
itself the hit is counted as "not recorded". When a photon between the beam line and the detector was not
stored, the chain seems to start at the detector: origin "away from the beam line (chain not stored)". The
shares of these two bound what the table can say.

Rates are hits / beam electrons x electrons per second, summed over the whole detector (MHz). Statistical errors
are sqrt(hits).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))

from solid_ana.table import _collections, decode, parse_encoding, read_encodings  # noqa: E402

E_PER_UA = 6.241509e12


def origin(r, z, zt):
    """region of a production point (cm); zt = z of the target centre"""
    dz = z - zt
    if r < 3.2 and dz < -24.7:
        return "entrance pipe and window"
    if r < 1.15 and abs(dz + 20.0) < 0.3:
        return "cell window, upstream"
    if r < 1.15 and abs(dz - 20.0) < 0.3:
        return "cell window, downstream"
    if abs(dz) < 20.0 and r < 0.9:
        return "3He gas"
    if abs(dz) < 20.0 and r < 1.15:
        return "cell wall"
    if -10.5 < dz < 0.5 and 1.15 <= r < 11.0:
        return "collimator, upstream"
    if 29.5 < dz < 50.5 and 1.15 <= r < 5.6:
        return "collimator, downstream"
    if r < 1.2 and abs(dz - 25.0) < 0.4:
        return "exit window"
    if abs(dz) < 60.0 and r < 80.0:
        return "target region, other (N2, coils)"
    if z > 540.0 and r < 25.0:
        return "beam stop"
    if dz >= 25.0 and r < 20.0:
        return "exit beam pipe"
    return "away from the beam line (chain not stored)"


def kind(pdg):
    a = abs(int(pdg))
    if a == 22:
        return "photon"
    if pdg == 11:
        return "electron"
    if pdg == -11:
        return "positron"
    if a in (12, 14, 16):
        return "neutrino"
    if a == 13:
        return "muon"
    return "hadron"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--electrons-per-event", type=float, required=True)
    ap.add_argument("--current", type=float, required=True)
    ap.add_argument("--emin", action="append", default=[], metavar="NAME=GEV")
    ap.add_argument("--target-z", type=float, default=-350.0)
    ap.add_argument("-o", "--output")
    a = ap.parse_args()

    import awkward as ak
    import uproot

    emin = [("GEM", 26e-9), ("MRPC", 26e-9), ("FASPD", 0.6e-3), ("LASPD", 2.0e-3)]
    for item in a.emin:
        k, v = item.split("=")
        emin.insert(0, (k, float(v)))

    def threshold(name):
        return next((v for k, v in emin if k in name), 0.0)

    counts = {}        # detector -> {(origin, first, at detector): hits}
    nev = 0
    mc = "MCParticles/MCParticles."
    for path in a.files:
        f = uproot.open(path)
        tree = f["events"]
        enc = read_encodings(f)
        cols = [n for n, k in _collections(tree).items() if k == "tracker" and "Cherenkov" not in n]
        nev += tree.num_entries
        for start in range(0, tree.num_entries, 100):
            names = [mc + "PDG", mc + "generatorStatus", mc + "vertex.x", mc + "vertex.y", mc + "vertex.z", mc + "parents_begin",
                     mc + "parents_end", "_MCParticles_parents/_MCParticles_parents.index"]
            for c in cols:
                names += ["%s/%s.eDep" % (c, c), "%s/%s.cellID" % (c, c), "_%s_particle/_%s_particle.index" % (c, c)]
            arr = tree.arrays(names, entry_start=start, entry_stop=start + 100, library="ak")
            for ev in range(len(arr)):
                pdg = vx = None
                for c in cols:
                    ed = ak.to_numpy(arr["%s/%s.eDep" % (c, c)][ev])
                    sel = np.nonzero(ed >= threshold(c))[0]
                    if len(sel) == 0:
                        continue
                    if pdg is None:
                        pdg = ak.to_numpy(arr[mc + "PDG"][ev]); gs = ak.to_numpy(arr[mc + "generatorStatus"][ev])
                        vx = ak.to_numpy(arr[mc + "vertex.x"][ev]) / 10.0; vy = ak.to_numpy(arr[mc + "vertex.y"][ev]) / 10.0
                        vz = ak.to_numpy(arr[mc + "vertex.z"][ev]) / 10.0
                        pb = ak.to_numpy(arr[mc + "parents_begin"][ev]); pe = ak.to_numpy(arr[mc + "parents_end"][ev])
                        par = ak.to_numpy(arr["_MCParticles_parents/_MCParticles_parents.index"][ev])
                    pidx = ak.to_numpy(arr["_%s_particle/_%s_particle.index" % (c, c)][ev])
                    cid = ak.to_numpy(arr["%s/%s.cellID" % (c, c)][ev])
                    e = parse_encoding(enc[c]) if enc.get(c) else {}
                    layer = decode(cid, e, "layer") if "layer" in e else np.ones(len(cid), np.int64)
                    for h in sel:
                        det = "%s plane %d" % (c.replace("Hits", ""), layer[h]) if "GEM" in c else c.replace("Hits", "")
                        i = int(pidx[h])
                        if i < 0 or gs[i] != 0:
                            key = ("not recorded", "not recorded", "not recorded")
                        else:
                            hit_pdg, hit_r, hit_z = pdg[i], np.hypot(vx[i], vy[i]), vz[i]
                            parent_pdg = 0
                            j, first = i, i
                            for _ in range(200):
                                if pe[j] <= pb[j]:
                                    break
                                p = int(par[pb[j]])
                                if j == i:
                                    parent_pdg = pdg[p]
                                if gs[p] != 0:              # reached the beam electron
                                    first = j
                                    break
                                j = first = p
                            org = origin(np.hypot(vx[first], vy[first]), vz[first], a.target_z)
                            k = kind(hit_pdg)
                            if k in ("electron", "positron"):
                                born = origin(hit_r, hit_z, a.target_z)
                                if first == i or not born.startswith("away from"):
                                    how = k + " straight from the beam line"
                                elif abs(int(parent_pdg)) == 22:
                                    how = k + " made by a photon near the detector"
                                else:
                                    how = k + " from a shower near the detector"
                            else:
                                how = k
                            key = (org, kind(pdg[first]), how)
                        counts.setdefault(det, {})
                        counts[det][key] = counts[det].get(key, 0) + 1
    electrons = nev * a.electrons_per_event
    scale = a.current * E_PER_UA / electrons / 1e6
    print("%d events x %g electrons = %.3g beam electrons; rates at %g uA, whole detector, MHz" % (nev, a.electrons_per_event, electrons, a.current))

    def table(title, index, dets):
        rows = {}
        for d in dets:
            for key, n in counts[d].items():
                rows.setdefault(key[index], {})
                rows[key[index]][d] = rows[key[index]].get(d, 0) + n
        tot = {d: sum(counts[d].values()) for d in dets}
        order = sorted(rows, key=lambda r: -sum(rows[r].values()))
        print("\n%s (share of the hits in %%)" % title)
        print("  %-46s " % "" + " ".join("%10s" % d.replace("GEMTracker", "GEM").replace(" plane ", "-")[:10] for d in dets))
        print("  %-46s " % "hits" + " ".join("%10d" % tot[d] for d in dets))
        print("  %-46s " % "rate (MHz)" + " ".join("%10.0f" % (tot[d] * scale) for d in dets))
        for r in order:
            print("  %-46s " % r + " ".join("%10.0f" % (100.0 * rows[r].get(d, 0) / tot[d]) if tot[d] else "%10s" % "-" for d in dets))

    dets = sorted(counts)
    table("By origin: where the beam electron did what led to the hit", 0, dets)
    table("By first secondary", 1, dets)
    table("By what made the hit", 2, dets)
    if a.output:
        out = {d: {" | ".join(k): n for k, n in counts[d].items()} for d in dets}
        with open(a.output, "w") as f:
            json.dump({"electrons": electrons, "current_uA": a.current, "hits": out}, f, indent=1)
            f.write("\n")
        print("\nwrote %s" % a.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
