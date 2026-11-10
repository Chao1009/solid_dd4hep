"""Per-mesh topology: boundary edges (cracks/T-junctions), non-manifold edges, inconsistent orientation."""
import sys
import numpy as np
def stats(V, F):
    F = np.abs(F)
    q = F[:, 3] > 0
    nxt = np.where(q[:, None], F[:, [1, 2, 3, 0]], F[:, [1, 2, 0, 0]])
    valid = np.ones(F.shape, bool); valid[:, 3] = q
    ea, eb = F[valid] - 1, nxt[valid] - 1
    # weld exact duplicates first (lineview does the same)
    key = np.round(V / 1e-6).astype(np.int64)
    _, inv = np.unique(key, axis=0, return_inverse=True); inv = inv.reshape(-1)
    ea, eb = inv[ea], inv[eb]
    keep = ea != eb
    ea, eb = ea[keep], eb[keep]
    lo, hi = np.minimum(ea, eb), np.maximum(ea, eb)
    k = lo.astype(np.int64) * (len(V) + 1) + hi
    u, idx, cnt = np.unique(k, return_inverse=True, return_counts=True)
    fwd = (ea < eb).astype(int)
    s = np.bincount(idx, weights=fwd, minlength=len(u))
    incons = ((cnt == 2) & (s != 1)).sum()
    return len(u), int((cnt == 1).sum()), int((cnt > 2).sum()), int(incons)
for p in sys.argv[1:]:
    z = np.load(p)
    vo, fo = z['mesh_vertex_offset'], z['mesh_facet_offset']
    print(p)
    tot = np.zeros(4, int); bad = 0
    for k in range(len(vo) - 1):
        V = z['mesh_vertices'][vo[k]:vo[k + 1]].astype(float); F = z['mesh_facets'][fo[k]:fo[k + 1]]
        e, b, nm, inc = stats(V, F)
        tot += (e, b, nm, inc)
        if b or nm or inc:
            bad += 1
            print('   %-22s %-45s edges %6d boundary %5d nonmanifold %4d inconsistent %4d' % (z['mesh_type'][k], str(z['mesh_name'][k])[:45], e, b, nm, inc))
    print('   total edges %d boundary %d nonmanifold %d inconsistent %d; meshes with defects %d of %d' % (*tot, bad, len(vo) - 1))
