import sys, importlib.machinery, importlib.util, collections
import numpy as np
LV = '$HOME/Projects/solid/solid_dd4hep/bin/solid_lineview'
spec = importlib.util.spec_from_loader('slv', importlib.machinery.SourceFileLoader('slv', LV))
slv = importlib.util.module_from_spec(spec); spec.loader.exec_module(slv)
a = slv.parse_prim(sys.argv[1]); b = slv.parse_prim(sys.argv[2])
z = np.load(sys.argv[3])
names = list(z['node_name']); mesh = z['node_mesh']; mtype = z['mesh_type']; mname = z['mesh_name']
drawn = [i for i in range(len(names)) if mesh[i] >= 0]
c = collections.Counter()
for k, (sa, sb) in enumerate(zip(a.solids, b.solids)):
    if sa.kind != 'Polyhedron': continue
    Ea = a.Fvis[sa.f0:sa.f1]; Eb = b.Fvis[sb.f0:sb.f1]
    if np.array_equal(a.F[sa.f0:sa.f1], b.F[sb.f0:sb.f1]) and not np.array_equal(Ea, Eb):
        i = drawn[k]
        c[(str(mtype[mesh[i]]), str(mname[mesh[i]]), int((Ea != Eb).sum()), int(Ea.sum()), int(Eb.sum()))] += 1
for k, v in c.items(): print(v, k)
