import sys, importlib.machinery, importlib.util, collections
import numpy as np
LV = '$HOME/Projects/solid/solid_dd4hep/bin/solid_lineview'
spec = importlib.util.spec_from_loader('slv', importlib.machinery.SourceFileLoader('slv', LV))
slv = importlib.util.module_from_spec(spec); spec.loader.exec_module(slv)
a = slv.parse_prim(sys.argv[1]); b = slv.parse_prim(sys.argv[2])
stats = collections.Counter(); ex = []
for sa, sb in zip(a.solids, b.solids):
    Va = a.V[sa.v0:sa.v1]; Vb = b.V[sb.v0:sb.v1]; Fa = a.F[sa.f0:sa.f1]; Fb = b.F[sb.f0:sb.f1]
    same_shape = len(Va) == len(Vb) and len(Fa) == len(Fb)
    if not same_shape:
        stats['count differs %s/%s' % (sa.kind, sb.kind)] += 1; ex.append((sa.name, sa.kind, len(Va), len(Vb), len(Fa), len(Fb))); continue
    dv = np.abs(Va - Vb).max() if len(Va) else 0
    sameF = np.array_equal(Fa, Fb); sameE = np.array_equal(a.Fvis[sa.f0:sa.f1], b.Fvis[sb.f0:sb.f1])
    if dv < 1e-6 and sameF and sameE: stats['identical'] += 1
    elif sameF and sameE: stats['same topology, max |dV| <= %.0e' % (10 ** np.ceil(np.log10(max(dv, 1e-12))))] += 1
    elif sameF: stats['same facets, edge flags differ (%s/%s)' % (sa.kind, sb.kind)] += 1
    else:
        # same vertex set but different facet order?
        stats['facets differ (%s/%s)' % (sa.kind, sb.kind)] += 1
        if len(ex) < 6: ex.append((sa.name, sa.kind, 'dv=%.2g' % dv))
for k, v in sorted(stats.items()): print('%6d  %s' % (v, k))
print(ex[:10])
