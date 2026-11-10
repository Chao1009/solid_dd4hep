import sys, ROOT
from functools import lru_cache
ROOT.gErrorIgnoreLevel = ROOT.kWarning
gm = ROOT.TGeoManager.Import(sys.argv[1])
memo = {}
def unrolled(v):
    k = v.GetName()
    if k in memo: return memo[k]
    n = 1
    for i in range(v.GetNdaughters()):
        n += unrolled(v.GetNode(i).GetVolume())
    memo[k] = n
    return n
top = gm.GetTopVolume()
print("unrolled", unrolled(top))
def show(v, depth, maxd):
    # group daughters by volume name
    from collections import Counter
    c = Counter(v.GetNode(i).GetVolume().GetName() for i in range(v.GetNdaughters()))
    for name, cnt in c.most_common():
        dv = gm.GetVolume(name)
        if memo[name] * cnt < 200 and depth > 0: continue
        print("  " * depth + f"{cnt:5d} x {name:45s} shape={dv.GetShape().ClassName():18s} unrolled={memo[name]:7d} tot={memo[name]*cnt:8d} col={dv.GetLineColor()} vis={dv.IsVisible()} transp={dv.GetTransparency()}")
        if depth < maxd: show(dv, depth + 1, maxd)
show(top, 0, int(sys.argv[2]) if len(sys.argv) > 2 else 4)
