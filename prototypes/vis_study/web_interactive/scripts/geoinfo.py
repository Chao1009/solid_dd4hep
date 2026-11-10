import sys, ROOT
f = ROOT.TFile.Open(sys.argv[1])
f.ls()
gm = None
for k in f.GetListOfKeys():
    o = k.ReadObj()
    if o.InheritsFrom("TGeoManager"):
        gm = o
print("TGeoManager:", gm.GetName(), "volumes", gm.GetListOfVolumes().GetEntries(), "shapes", gm.GetListOfShapes().GetEntries(),
      "nodes(phys)", gm.CountNodes(gm.GetTopVolume(), 100, 0) if False else "", "materials", gm.GetListOfMaterials().GetEntries())
top = gm.GetTopVolume()
print("top", top.GetName(), top.GetShape().ClassName(), "ndaughters", top.GetNdaughters())
def walk(v, depth, maxd, cnt):
    cnt[0] += 1
    for i in range(v.GetNdaughters()):
        n = v.GetNode(i)
        if depth < maxd:
            walk(n.GetVolume(), depth + 1, maxd, cnt)
        else:
            cnt[0] += 1
cnt=[0]; walk(top, 0, 30, cnt); print("total node count (unrolled)", cnt[0])
for i in range(top.GetNdaughters()):
    n = top.GetNode(i); v = n.GetVolume()
    print(f"  {n.GetName():40s} {v.GetShape().ClassName():20s} nd={v.GetNdaughters():4d} color={v.GetLineColor()} transp={v.GetTransparency()} vis={v.IsVisible()}")
from collections import Counter
c = Counter(s.ClassName() for s in gm.GetListOfShapes())
print(c.most_common())
