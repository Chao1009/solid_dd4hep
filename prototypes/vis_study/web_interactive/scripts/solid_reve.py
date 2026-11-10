#!/usr/bin/env python3
"""Prototype ROOT EVE-7 (REve) scene for SoLID-Lite: TGeo shapes, helix tracks propagated in a
solenoid field (REveTrackPropagator), tracker hits, ECal modules as an REveBoxSet coloured by energy.

usage: solid_reve.py display_geo.root events.firebird.json EVENT PORT
Runs a web server (RWebDisplayArgs 'server' mode) and prints the URL; stays alive until killed.
"""
import json, sys, time

import ROOT

geo_file, dex_file, iev, port = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ROOT.gEnv.SetValue("WebGui.HttpPort", port)
ROOT.gEnv.SetValue("WebGui.HttpPortMin", port)
ROOT.gEnv.SetValue("WebGui.HttpPortMax", port + 5)
ROOT.gEnv.SetValue("WebGui.HttpLoopback", "yes")
ROOT.gEnv.SetValue("WebGui.OnetimeKey", "no")

ROOT.gInterpreter.Declare(r"""
#include <ROOT/REveTrackPropagator.hxx>
// crude CLEO-II-like solenoid: uniform Bz inside the coil (r < 145 cm, -180 < z < 180 cm), zero outside
class SolidField : public ROOT::Experimental::REveMagField {
public:
   double fB{-1.5};
   SolidField() { fFieldConstant = kFALSE; }
   ROOT::Experimental::REveVectorD GetField(Double_t x, Double_t y, Double_t z) const override {
      double r2 = x * x + y * y;
      if (r2 < 145. * 145. && z > -180. && z < 180.) return {0., 0., fB};
      return {0., 0., 0.};
   }
   Double_t GetMaxFieldMag() const override { return std::abs(fB); }
};
""")

REX = ROOT.Experimental
t0 = time.time()
eve = REX.REveManager.Create()

# ---------------- geometry: top-level volumes of the display TGeo file ---------------------------
gm = ROOT.TGeoManager.Import(geo_file)
geom = REX.REveElement("SoLID-Lite geometry")
eve.GetGlobalScene().AddElement(geom)
nshapes = 0


def add_node(node, mother_mat, depth, parent):
    global nshapes
    mat = ROOT.TGeoHMatrix(mother_mat)
    mat.Multiply(node.GetMatrix())
    vol = node.GetVolume()
    name = vol.GetName()
    # skip the hundreds of identical ECal/SPD/GEM modules; draw their envelopes instead
    skip_daughters = any(s in name for s in ("Shower_envelope", "SPD_assembly", "GEMTracker", "MRPC"))
    if vol.GetShape().ClassName() != "TGeoShapeAssembly" and (vol.IsVisible() or skip_daughters):
        sh = REX.REveGeoShape(name)
        sh.SetShape(vol.GetShape().Clone())
        col = vol.GetLineColor()
        sh.SetMainColor(col if col > 1 else ROOT.kGray + 1)
        sh.SetMainTransparency(70 if skip_daughters else 40)
        sh.RefMainTrans().SetFrom(mat)
        parent.AddElement(sh)
        nshapes += 1
    if depth < 4 and not skip_daughters:
        for i in range(vol.GetNdaughters()):
            add_node(vol.GetNode(i), mat, depth + 1, parent)


top = gm.GetTopVolume()
for i in range(top.GetNdaughters()):
    add_node(top.GetNode(i), ROOT.TGeoHMatrix(), 1, geom)

# ---------------- event data (from the DEX file written by edm4hep2web.py) ------------------------
dex = json.load(open(dex_file))
ev = [e for e in dex["events"] if e["id"] == iev][0]
pieces = {p["name"]: p for p in ev["pieces"]}

prop = REX.REveTrackPropagator("SolidProp")
field = ROOT.SolidField()
ROOT.SetOwnership(field, False)
prop.SetMagFieldObj(field, False)
prop.SetMaxR(400); prop.SetMaxZ(500); prop.SetMaxOrbs(2)
tracks = REX.REveTrackList("MC particles (helix in 1.5 T solenoid)", prop)
tracks.SetMainColor(ROOT.kRed)
tracks.SetLineWidth(2)
ntr = 0
mc = pieces.get("McHitTrajectories")
if mc:
    c = mc["columns"]
    for k in range(mc["count"]):
        if c["p"][k] < 50:
            continue
        v = mc["points"][k][0]
        rt = REX.REveRecTrackD()
        rt.fV.Set(v[0] / 10, v[1] / 10, v[2] / 10)
        rt.fP.Set(c["px"][k] / 1000, c["py"][k] / 1000, c["pz"][k] / 1000)
        rt.fSign = int(c["charge"][k])
        tr = REX.REveTrack(rt, prop)
        tr.SetMainColor(ROOT.kRed if c["pdg"][k] == 11 else ROOT.kBlue)
        tr.SetName(f"pdg {c['pdg'][k]} p={c['p'][k]:.0f} MeV")
        tracks.AddElement(tr)
        ntr += 1
tracks.MakeTracks()
eve.GetEventScene().AddElement(tracks)

hits = REX.REvePointSet("tracker hits")
hits.SetMarkerColor(ROOT.kAzure + 1); hits.SetMarkerSize(6)
nh = 0
for name, p in pieces.items():
    if p["type"] == "BoxHit" and "EC" not in name:
        pos = p["columns"]["pos"]
        for k in range(p["count"]):
            hits.SetNextPoint(pos[3 * k] / 10, pos[3 * k + 1] / 10, pos[3 * k + 2] / 10)
            nh += 1
eve.GetEventScene().AddElement(hits)

pal = REX.REveRGBAPalette(0, 300)  # MeV
boxes = REX.REveBoxSet("ECal modules (E sum, MeV)")
boxes.Reset(REX.REveBoxSet.kBT_InstancedScaled, False, 64)
boxes.SetPalette(pal)
nb = 0
for name, p in pieces.items():
    if name.endswith("ShowerHitsModules"):
        cc = p["columns"]
        for k in range(p["count"]):
            x, y, z = (cc["pos"][3 * k + j] / 10 for j in range(3))
            w, h, d = (cc["dim"][3 * k + j] / 10 for j in range(3))
            boxes.AddBox(x - w / 2, y - h / 2, z - d / 2, w, h, d)
            boxes.DigitValue(int(cc["edep"][k] * 1000))
            nb += 1
boxes.RefitPlex()
eve.GetEventScene().AddElement(boxes)

print(f"REve scene: {nshapes} geometry shapes, {ntr} tracks, {nh} hits, {nb} ECal boxes, built in {time.time() - t0:.1f}s",
      flush=True)
eve.Show(ROOT.RWebDisplayArgs("server"))
print("SHOWN", flush=True)
if len(sys.argv) > 6 and sys.argv[6] == "loop":
    t_end = time.time() + float(sys.argv[5])
    while time.time() < t_end:
        ROOT.gSystem.ProcessEvents()
        time.sleep(0.02)
else:
    ROOT.gApplication.Run(True)  # ROOT event loop (serves http/websocket requests); kill to stop
