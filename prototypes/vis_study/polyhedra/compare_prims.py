"""Compare two .prim scenes volume by volume (as parsed by solid_lineview)."""
import sys, time, importlib.machinery, importlib.util, collections
import numpy as np
LV = '$HOME/Projects/solid/solid_dd4hep/bin/solid_lineview'
spec = importlib.util.spec_from_loader('slv', importlib.machinery.SourceFileLoader('slv', LV))
slv = importlib.util.module_from_spec(spec); spec.loader.exec_module(slv)

def load(p):
    t = time.time(); sc = slv.parse_prim(p); dt = time.time() - t
    return sc, dt

a, ta = load(sys.argv[1]); b, tb = load(sys.argv[2])
print('parse times: %.2f s (%s)  %.2f s (%s)' % (ta, sys.argv[1], tb, sys.argv[2]))
for name, sc in (('A', a), ('B', b)):
    print(name, 'solids', len(sc.solids), 'V', len(sc.V), 'F', len(sc.F), dict(sc.counts), 'warnings', dict(sc.warnings))
na = [s.name for s in a.solids]; nb = [s.name for s in b.solids]
print('same name sequence:', na == nb)
ca, cb = collections.Counter(na), collections.Counter(nb)
print('names only in A:', list((ca - cb).items())[:10]); print('names only in B:', list((cb - ca).items())[:10])
# per-solid comparison in order
nd = 0; maxd = 0; kinds = collections.Counter(); colour_diff = 0; depth_diff = 0; exact = 0; facet_diff = 0
worst = []
for sa, sb in zip(a.solids, b.solids):
    if sa.name != sb.name:
        nd += 1; continue
    if np.max(np.abs(np.array(sa.color) - np.array(sb.color))) > 1e-6: colour_diff += 1
    if sa.depth != sb.depth: depth_diff += 1
    Va = a.V[sa.v0:sa.v1]; Vb = b.V[sb.v0:sb.v1]
    d = max(np.max(np.abs(sa.bmin - sb.bmin)), np.max(np.abs(sa.bmax - sb.bmax)))
    maxd = max(maxd, d)
    if d > 1e-3: worst.append((d, sa.name, sa.kind, sb.kind))
    Fa = a.F[sa.f0:sa.f1]; Fb = b.F[sb.f0:sb.f1]
    if len(Va) == len(Vb) and len(Fa) == len(Fb) and np.allclose(Va, Vb, atol=1e-6) and np.array_equal(Fa, Fb) \
            and np.array_equal(a.Fvis[sa.f0:sa.f1], b.Fvis[sb.f0:sb.f1]):
        exact += 1
    else:
        kinds[(sa.kind, sb.kind)] += 1
        if sa.kind == sb.kind == 'Polyhedron': facet_diff += 1
print('name mismatches in order:', nd, ' colour diffs:', colour_diff, ' depth diffs:', depth_diff)
print('identical meshes (V, F, edge flags):', exact, ' differing by kind:', dict(kinds), ' polyhedron-vs-polyhedron diffs:', facet_diff)
print('max bbox difference [mm]: %.3g' % maxd); worst.sort(reverse=True); print('worst:', worst[:8])
