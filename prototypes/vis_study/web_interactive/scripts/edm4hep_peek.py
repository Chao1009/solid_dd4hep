import sys, uproot, awkward as ak, numpy as np
t = uproot.open(sys.argv[1])["events"]
ev = int(sys.argv[2]) if len(sys.argv) > 2 else 0
def arr(b): return t[b].array(entry_start=ev, entry_stop=ev+1)[0]
for c in ["LAECShowerHits", "LAECPreShowerHits", "GEMTrackerLAHits", "LASPDHits"]:
    keys = [k for k in t.keys() if k.startswith(c + "/")]
    print(c, [k.split(".",1)[1] for k in keys])
    if not keys: continue
    if "Shower" in c:
        x, y, z, e = (np.asarray(arr(f"{c}/{c}.{f}")) for f in ["position.x", "position.y", "position.z", "energy"])
        cid = np.asarray(arr(f"{c}/{c}.cellID"))
        print(f"  n={len(x)} unique cellIDs={len(set(cid.tolist()))} unique xy={len(set(zip(np.round(x,1),np.round(y,1))))} unique z={len(set(np.round(z,1)))} Esum={e.sum():.4f} GeV  z range {z.min():.1f}..{z.max():.1f}")
        order = np.argsort(-e)[:5]
        for i in order: print(f"   cell {cid[i]:#x} pos=({x[i]:.1f},{y[i]:.1f},{z[i]:.1f}) E={e[i]*1000:.2f} MeV")
        # contributions
        cb = t[f"{c}/{c}.contributions_begin"].array(entry_start=ev, entry_stop=ev+1)[0]
        ce = t[f"{c}/{c}.contributions_end"].array(entry_start=ev, entry_stop=ev+1)[0]
        cc = c + "Contributions"
        ck = [k.split(".",1)[1] for k in t.keys() if k.startswith(cc + "/")]
        print("  contributions fields", ck)
        ct = np.asarray(arr(f"{cc}/{cc}.time"))
        print(f"  ncontrib={len(ct)} time range {ct.min():.2f}..{ct.max():.2f} ns; per-hit contrib counts {np.asarray(ce)[:10]-np.asarray(cb)[:10]}")
    else:
        x, y, z, tm = (np.asarray(arr(f"{c}/{c}.{f}")) for f in ["position.x", "position.y", "position.z", "time"])
        for i in range(len(x)): print(f"   hit pos=({x[i]:.1f},{y[i]:.1f},{z[i]:.1f}) t={tm[i]:.2f}")
mc = "MCParticles"
print([k.split(".",1)[1] for k in t.keys() if k.startswith(mc + "/")])
pdg = np.asarray(arr(f"{mc}/{mc}.PDG")); gs = np.asarray(arr(f"{mc}/{mc}.generatorStatus"))
px, py, pz = (np.asarray(arr(f"{mc}/{mc}.momentum.{a}")) for a in "xyz")
p = np.sqrt(px**2+py**2+pz**2)
print("n", len(pdg), "gen", (gs==1).sum(), "p>10MeV", (p>0.01).sum(), "p>100MeV", (p>0.1).sum())
from collections import Counter
print(Counter(pdg.tolist()).most_common(8))
