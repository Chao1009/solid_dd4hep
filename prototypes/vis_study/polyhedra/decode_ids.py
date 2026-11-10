import sys, json, collections
import numpy as np
z = np.load(sys.argv[1])
ros = json.loads(str(z['meta']))['readouts']
def decode(desc, v):
    out = {}; off = 0
    for tok in desc.split(','):
        p = tok.split(':')
        if len(p) == 3: off = int(p[1]); w = int(p[2])
        else: w = int(p[1])
        signed = w < 0; w = abs(w)
        x = (int(v) >> off) & ((1 << w) - 1)
        if signed and x >= (1 << (w - 1)): x -= (1 << w)
        out[p[0]] = x; off += w
    return out
sens = np.nonzero(z['node_sensitive'] & (z['node_mesh'] >= 0))[0]
by = collections.Counter(); shown = collections.Counter()
for i in sens:
    sd = str(z['node_subdet'][i]); desc = ros[sd]['fields']
    mine = int(z['node_volume_id'][i]); g = int(z['node_g4_volume_id'][i]) & ((1 << 64) - 1)
    ok = mine == g
    by[(sd, ok)] += 1
    if not ok and shown[sd] < 2:
        shown[sd] += 1
        print(sd, z['node_tgeo_path'][i]); print('   mine', decode(desc, mine)); print('   ddg4', decode(desc, g))
print(dict(by))
