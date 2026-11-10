import ROOT, sys, dd4hep
desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML(sys.argv[1])
top = ROOT.gGeoManager.GetTopNode()
n = 0
def rec(node, depth):
    global n
    pv = ROOT.dd4hep.PlacedVolume(node)
    vol = pv.volume()
    sd = vol.sensitiveDetector()
    if sd.isValid():
        n += 1
        if n <= 5:
            print(depth, node.GetName(), sd.name(), (sd.readout() if callable(sd.readout) else sd.readout).name(), [ (str(p.first), p.second) for p in pv.volIDs()])
    for i in range(node.GetNdaughters()):
        rec(node.GetDaughter(i), depth+1)
rec(top, 0)
print('sensitive nodes', n)
