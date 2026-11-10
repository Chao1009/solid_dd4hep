#!/usr/bin/env python3
"""Build a web-display TGeo file from a DD4hep compact description.

 - loads the compact XML with DD4hep (all SoLID plugins),
 - prunes the 194+194 absorber/scintillator layers inside every ECal shower module
   (5.07 M of the 5.08 M physical nodes) so that browsers only see the module prisms,
 - maps DD4hep's run-time TColor indices (>= 1000, which are not stored in the ROOT file)
   to the closest standard ROOT colour that JSROOT knows,
 - exports TGeoManager to a ROOT file.

usage: make_display_geo.py compact.xml out.root [--keep-layers]
"""
import fnmatch, sys, time

import ROOT

t0 = time.time()
compact, out = sys.argv[1], sys.argv[2]
keep_layers = "--keep-layers" in sys.argv
ROOT.gErrorIgnoreLevel = ROOT.kWarning
import dd4hep  # noqa: E402

desc = dd4hep.Detector.getInstance()
desc.fromXML(compact)
gm = desc.manager()
t1 = time.time()


def unrolled(v, memo):
    k = v.GetName()
    if k not in memo:
        memo[k] = 1 + sum(unrolled(v.GetNode(i).GetVolume(), memo) for i in range(v.GetNdaughters()))
    return memo[k]


n_before = unrolled(gm.GetTopVolume(), {})

# --- prune calorimeter layers ------------------------------------------------------------
pruned = []
if not keep_layers:
    for v in gm.GetListOfVolumes():
        if fnmatch.fnmatch(v.GetName(), "*Shower_module") and v.GetNdaughters() > 0:
            pruned.append((v.GetName(), v.GetNdaughters()))
            v.ClearNodes()
n_after = unrolled(gm.GetTopVolume(), {})

# --- colour remap: run-time colours (index >= 1000 or > 100) -> nearest standard 1..50 ------
std = {}
for i in list(range(1, 50)):
    c = ROOT.gROOT.GetColor(i)
    if c:
        std[i] = (c.GetRed(), c.GetGreen(), c.GetBlue())
remap = {}
for v in gm.GetListOfVolumes():
    ci = v.GetLineColor()
    if ci <= 50 or ci in remap:
        if ci in remap:
            v.SetLineColor(remap[ci]); v.SetFillColor(remap[ci])
        continue
    c = ROOT.gROOT.GetColor(ci)
    if not c:
        continue
    rgb = (c.GetRed(), c.GetGreen(), c.GetBlue())
    best = min(std, key=lambda k: sum((a - b) ** 2 for a, b in zip(std[k], rgb)))
    remap[ci] = best
    v.SetLineColor(best); v.SetFillColor(best)

gm.Export(out)
t2 = time.time()
print(f"loaded compact in {t1 - t0:.1f}s; nodes {n_before} -> {n_after}; pruned {len(pruned)} volumes "
      f"{pruned}; colour remap {remap}; export {t2 - t1:.1f}s")
