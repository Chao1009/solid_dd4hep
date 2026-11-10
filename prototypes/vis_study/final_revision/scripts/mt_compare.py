"""Compare EDM4hep outputs of two ddsim runs event by event (keyed by EventHeader.eventNumber).
usage: python3 -I mt_compare.py a.edm4hep.root b.edm4hep.root
Prints, per event: primary momentum, #MCParticles, #hits and summed energy per hit collection."""
import sys
import numpy as np
import uproot


def summary(path):
    t = uproot.open(path)['events']
    keys = set(t.keys())
    evn = t['EventHeader/EventHeader.eventNumber'].array(library='np')
    evn = [int(x[0]) for x in evn]
    px = t['MCParticles/MCParticles.momentum.x'].array(library='np')
    py = t['MCParticles/MCParticles.momentum.y'].array(library='np')
    pz = t['MCParticles/MCParticles.momentum.z'].array(library='np')
    cols = sorted({k.split('/')[0] for k in keys if k.endswith('.cellID') and not k.startswith('_')})
    out = {}
    for i, n in enumerate(evn):
        rec = {'p0': (round(float(px[i][0]), 9), round(float(py[i][0]), 9), round(float(pz[i][0]), 9)),
               'nMC': len(px[i])}
        for c in cols:
            ek = c + '/' + c + '.energy' if (c + '/' + c + '.energy') in keys else c + '/' + c + '.eDep'
            e = t[ek].array(library='np', entry_start=i, entry_stop=i + 1)[0]
            rec[c] = (len(e), round(float(np.sum(e)), 9))
        out[n] = rec
    return out


a, b = summary(sys.argv[1]), summary(sys.argv[2])
same = 0
for n in sorted(set(a) | set(b)):
    ra, rb = a.get(n), b.get(n)
    eq = ra == rb
    same += eq
    print('event %d: %s' % (n, 'IDENTICAL' if eq else 'DIFFERENT'))
    if not eq:
        for k in sorted(set(ra or {}) | set(rb or {})):
            va, vb = (ra or {}).get(k), (rb or {}).get(k)
            if va != vb:
                print('   %-22s %s  |  %s' % (k, va, vb))
    else:
        print('   p0 %s  nMC %d  %s' % (ra['p0'], ra['nMC'], {k: v for k, v in ra.items() if k not in ('p0', 'nMC')}))
print('identical events: %d of %d' % (same, len(set(a) | set(b))))
