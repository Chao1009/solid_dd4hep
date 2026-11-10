"""How many drawn nodes match each hit by (cellID & mask) == volid? (reviewer check of 'exactly one')"""
import sys, collections
import numpy as np
import uproot
z = np.load(sys.argv[1])
ro = z['node_readout']; vid = z['node_volume_id'].astype(np.uint64); msk = z['node_volume_mask'].astype(np.uint64)
drawn = z['node_mesh'] >= 0
t = uproot.open(sys.argv[2])['events']
cols = sorted({k.split('/')[0] for k in t.keys() if k.endswith('.cellID') and not k.startswith('_')})
for c in cols:
    cid = t[c + '/' + c + '.cellID'].array(library='np')
    if sum(len(x) for x in cid) == 0:
        continue
    cid = np.concatenate(cid).astype(np.uint64)
    cand = np.nonzero(drawn & (ro == c))[0]
    n = np.zeros(len(cid), int)
    for i in cand:
        n += ((cid & msk[i]) == vid[i])
    print('%-20s hits %5d  nmatch distribution %s  zero-mask candidates %d' % (c, len(cid), dict(collections.Counter(n.tolist())), int((msk[cand] == 0).sum())))
