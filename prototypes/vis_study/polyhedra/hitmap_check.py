"""Check that simulated hits map to drawn volumes via (cellID & node_volume_mask) == node_volume_id.
usage: python3 -I hitmap_check.py geometry.npz sim.edm4hep.root"""
import sys, json, collections
import numpy as np
import uproot
z = np.load(sys.argv[1])
ro = z['node_readout']; vid = z['node_volume_id'].astype(np.uint64); msk = z['node_volume_mask'].astype(np.uint64)
mesh = z['node_mesh']; M = z['node_matrix']
voff = z['mesh_vertex_offset']; V = z['mesh_vertices'].astype(float)
drawn = mesh >= 0
def world_bbox(i):
    k = mesh[i]; P = V[voff[k]:voff[k + 1]]
    W = P @ M[i, :9].reshape(3, 3).T + M[i, 9:]
    return W.min(0), W.max(0)
t = uproot.open(sys.argv[2])['events']
cols = sorted({k.split('/')[0] for k in t.keys() if k.endswith('.cellID') and not k.startswith('_')})
tot = collections.Counter()
for c in cols:
    cid = t[c + '/' + c + '.cellID'].array(library='np')
    if sum(len(x) for x in cid) == 0:
        continue
    cid = np.concatenate(cid).astype(np.uint64)
    try:
        pos = np.stack([np.concatenate(t[c + '/' + c + '.position.' + a].array(library='np')) for a in 'xyz'], 1)
    except Exception:
        pos = None
    cand = np.nonzero(drawn & (ro == c))[0]
    nmatch = np.zeros(len(cid), int); deepest = np.full(len(cid), -1)
    for i in cand:
        hit = (cid & msk[i]) == vid[i]
        nmatch += hit
        # prefer the node with the largest mask (deepest)
        better = hit & ((deepest < 0) | (msk[i] > msk[np.maximum(deepest, 0)]))
        deepest[better] = i
    inside = 0
    if pos is not None:
        for h in np.nonzero(deepest >= 0)[0]:
            lo, hi = world_bbox(deepest[h])
            if np.all(pos[h] >= lo - 1.0) and np.all(pos[h] <= hi + 1.0):
                inside += 1
    sample = z['node_det_path'][deepest[deepest >= 0][0]] if (deepest >= 0).any() else ''
    print('%-22s hits %5d  drawn candidates %5d  mapped %5d  unique-deepest %5d  hit position inside node bbox %5d   e.g. %s'
          % (c, len(cid), len(cand), (nmatch > 0).sum(), (deepest >= 0).sum(), inside, sample))
    tot['hits'] += len(cid); tot['mapped'] += int((nmatch > 0).sum()); tot['inside'] += inside
print(dict(tot))
