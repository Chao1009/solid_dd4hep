#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests of bin/solid_lite_config and of the SoLID-Lite coverage options (ctest: lite_config).

usage (inside eic-shell / solid-shell, with DETECTOR_PATH and the plugin library set): test_lite_config.py

  committed   solid_SIDIS_He3_lite.xml is what the tool writes with its defaults
  load        every option loads, builds detectors and has no overlap
  segments    with four 24-deg segments every forward detector has 8 of 30 sectors, or its share of modules and
              pads, and the large-angle detectors are unchanged; "none" removes the forward detectors
Exit status 0 if every check passes.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
TOOL = os.path.join(REPO, "bin", "solid_lite_config")
CHECK = os.path.join(HERE, "check_geometry.py")

OPTIONS = {   # name: arguments
    "he3": ["--experiment", "SIDIS_He3"],
    "he3_fullFA": ["--experiment", "SIDIS_He3", "--fa", "full"],
    "he3_stage1": ["--experiment", "SIDIS_He3", "--iron", "stage1"],
    "he3_mrpc": ["--experiment", "SIDIS_He3", "--mrpc", "on"],
    "jpsi_4fa": ["--experiment", "JPsi"],
    "jpsi_LAonly": ["--experiment", "JPsi", "--fa", "none"],
    "ddvcs_fullFA": ["--experiment", "DDVCS", "--fa", "full", "--mrpc", "on"],
}
FAILED = []


def check(ok, what):
    print("%s  %s" % ("ok  " if ok else "FAIL", what))
    if not ok:
        FAILED.append(what)


def children(compact):
    """number of daughter DetElements of every top-level detector, in a fresh process (one geometry per process)"""
    code = ("import dd4hep, sys\n"
            "d = dd4hep.Detector.getInstance(); d.fromXML(sys.argv[1])\n"
            "def n(e):\n"
            "    c = [x for _, x in e.children()]\n"
            "    return len(c) + sum(n(x) for x in c)\n"
            "for name, de in d.world().children():\n"
            "    print('COUNT', name, n(de))\n")
    out = subprocess.run([sys.executable, "-c", code, compact], capture_output=True, text=True).stdout
    return {l.split()[1]: int(l.split()[2]) for l in out.splitlines() if l.startswith("COUNT ")}


def main():
    with tempfile.TemporaryDirectory() as tmp:
        files = {}
        for name, args in OPTIONS.items():
            files[name] = os.path.join(tmp, name + ".xml")
            subprocess.run([sys.executable, TOOL] + args + ["-o", files[name]], check=True, capture_output=True)

        default = subprocess.run([sys.executable, TOOL, "--experiment", "SIDIS_He3"], check=True, capture_output=True,
                                 text=True).stdout
        committed = open(os.path.join(REPO, "solid_SIDIS_He3_lite.xml")).read()
        check(default == committed, "solid_SIDIS_He3_lite.xml is the tool's default output for SIDIS_He3")

        for name, f in files.items():
            r = subprocess.run([sys.executable, CHECK, f], capture_output=True, text=True)
            last = (r.stdout.strip().splitlines() or ["no output"])[-1]
            check(r.returncode == 0, "%s loads without overlap (%s)" % (name, last))

        seg, full = children(files["he3"]), children(files["he3_fullFA"])
        print("    DetElements below each detector, 4 x 24 deg: %s" % seg)
        print("    DetElements below each detector, full azimuth: %s" % full)
        for det in ("LAECShower", "LAECPreShower", "LASPD", "MagnetYoke", "Target"):
            check(seg.get(det) == full.get(det) and det in seg, "%s is the same in both (%s)" % (det, seg.get(det)))
        for det, lo, hi in (("FAECShower", 0.25, 0.29), ("FAECPreShower", 0.25, 0.29), ("FASPD", 0.26, 0.28)):
            ratio = seg[det] / float(full[det])
            check(lo < ratio < hi, "%s: %d of %d elements in the segments (%.3f; 4 x 24 / 360 = 0.267)"
                  % (det, seg[det], full[det], ratio))
        # GEM: the assembly holds the large-angle planes (unchanged) and the forward planes (8 of 30 sectors)
        gem_seg, gem_full = seg["GEMTrackerLiteAssembly"], full["GEMTrackerLiteAssembly"]
        check(gem_seg < gem_full, "GEM: fewer elements with segments (%d against %d)" % (gem_seg, gem_full))
        check(seg["LiteCherenkov"] < full["LiteCherenkov"], "Cherenkov: 8 sectors against 30 (%d against %d elements)"
              % (seg["LiteCherenkov"], full["LiteCherenkov"]))
        la = children(files["jpsi_LAonly"])
        check(not any(k.startswith("FA") or k == "LiteCherenkov" for k in la) and "LAECShower" in la,
              "fa = none: no forward detector is built (%s)" % ", ".join(sorted(la)))
        mr = children(files["he3_mrpc"])
        check("FAMRPC" in mr and "FAMRPC" not in seg, "mrpc = on adds the MRPC")
    print("%d check%s failed" % (len(FAILED), "" if len(FAILED) == 1 else "s") if FAILED else "all checks passed")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
