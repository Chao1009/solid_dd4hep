import sys
import numpy as np
import uproot
z = np.load(sys.argv[1])
ro = z['node_readout']; vid = z['node_volume_id'].astype(np.uint64); msk = z['node_volume_mask'].astype(np.uint64)
drawn = z['node_mesh'] >= 0
t = uproot.open(sys.argv[2])['events']
c = 'FAECPreShowerHits'
cid = np.concatenate(t[c + '/' + c + '.cellID'].array(library='np')).astype(np.uint64)
cand = np.nonzero(drawn & (ro == c))[0]
for h in cid[:2]:
    m = [i for i in cand if (h & msk[i]) == vid[i]]
    for i in m:
        print(hex(int(h)), z['node_name'][i], z['node_det_path'][i], hex(int(msk[i])), 'sensitive', bool(z['node_sensitive'][i]), 'depth', int(z['node_depth'][i]))
