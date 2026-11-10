import sys, uproot, awkward as ak
f = uproot.open(sys.argv[1])
print(f.keys())
t = f["events"]
tn = t.typenames(recursive=False)
for name, typ in sorted(tn.items()):
    if name.startswith("_"): continue
    try:
        n = ak.num(t[f"{name}/{name}.{'PDG' if 'MCParticle' in typ else 'cellID'}"].array(), axis=1).tolist()
    except Exception as e:
        n = "?"
    print(f"{name:40s} {typ:45s} per-event={n}")
