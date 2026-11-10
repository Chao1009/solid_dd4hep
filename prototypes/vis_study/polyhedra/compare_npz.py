"""Compare G4 (g4poly) and TGeo (tgeopoly) tessellations node by node (matched by TGeo node path)."""
import sys, json, collections
import numpy as np

def load(p):
    z = np.load(p)
    d = {k: z[k] for k in z.files}
    d['vo'] = d['mesh_vertex_offset']; d['fo'] = d['mesh_facet_offset']
    return d

def mesh(d, k):
    V = d['mesh_vertices'][d['vo'][k]:d['vo'][k + 1]].astype(float)
    F = d['mesh_facets'][d['fo'][k]:d['fo'][k + 1]]
    return V, np.abs(F)

def volume(V, F):
    q = F[:, 3] > 0
    tri = np.vstack([F[:, :3], F[q][:, [0, 2, 3]]]) - 1
    return np.einsum('ij,ij->i', V[tri[:, 0]], np.cross(V[tri[:, 1]], V[tri[:, 2]])).sum() / 6.0

def area(V, F):
    q = F[:, 3] > 0
    tri = np.vstack([F[:, :3], F[q][:, [0, 2, 3]]]) - 1
    return 0.5 * np.linalg.norm(np.cross(V[tri[:, 1]] - V[tri[:, 0]], V[tri[:, 2]] - V[tri[:, 0]]), axis=1).sum()

def wbox(d, i, k):
    V, _ = mesh(d, k); M = d['node_matrix'][i]
    W = V @ M[:9].reshape(3, 3).T + M[9:]
    return W.min(0), W.max(0)

g, t = load(sys.argv[1]), load(sys.argv[2])
gi = {str(p): i for i, p in enumerate(g['node_tgeo_path']) if g['node_mesh'][i] >= 0}
ti = {str(p): i for i, p in enumerate(t['node_tgeo_path']) if t['node_mesh'][i] >= 0}
both = sorted(set(gi) & set(ti))
print('drawn: G4 %d, TGeo %d, matched by TGeo path %d; only G4 %d, only TGeo %d' % (len(gi), len(ti), len(both),
      len(set(gi) - set(ti)), len(set(ti) - set(gi))))
print('  only G4:', sorted(set(gi) - set(ti))[:5]); print('  only TGeo:', sorted(set(ti) - set(gi))[:5])
pairs = collections.Counter()
maxbox = 0; badbox = []; colour = 0
for p in both:
    i, j = gi[p], ti[p]
    pairs[(int(g['node_mesh'][i]), int(t['node_mesh'][j]))] += 1
    a0, a1 = wbox(g, i, g['node_mesh'][i]); b0, b1 = wbox(t, j, t['node_mesh'][j])
    dd = max(np.abs(a0 - b0).max(), np.abs(a1 - b1).max()); maxbox = max(maxbox, dd)
    if dd > 0.5: badbox.append((round(float(dd), 2), p))
    if np.abs(g['node_rgba'][i][:3] - t['node_rgba'][j][:3]).max() > 1e-6: colour += 1
print('nodes with colour differences:', colour)
print('max world-bbox difference over matched nodes: %.3g mm; nodes > 0.5 mm: %d' % (maxbox, len(badbox)))
for x in sorted(badbox, reverse=True)[:10]: print('   ', x)
rows = []
for (kg, kt), n in sorted(pairs.items(), key=lambda x: -x[1]):
    Vg, Fg = mesh(g, kg); Vt, Ft = mesh(t, kt)
    vg, vt = volume(Vg, Fg), volume(Vt, Ft)
    rows.append((str(g['mesh_type'][kg]), str(t['mesh_type'][kt]), str(g['mesh_name'][kg])[:38], n, len(Vg), len(Vt), len(Fg), len(Ft),
                 vg / 1e3, vt / 1e3, 100.0 * (vt - vg) / vg if vg else float('nan'), area(Vg, Fg) / 1e2, area(Vt, Ft) / 1e2))
print('%-20s %-20s %-38s %5s %6s %6s %6s %6s %12s %12s %8s' % ('G4 solid', 'TGeo shape', 'name', 'inst', 'Vg4', 'Vtgeo', 'Fg4', 'Ftgeo', 'vol_g4[cm3]', 'vol_tgeo', 'dvol%'))
for r in rows:
    print('%-20s %-20s %-38s %5d %6d %6d %6d %6d %12.1f %12.1f %8.3f' % r[:11])
