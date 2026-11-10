import sys
from podio.root_io import Reader
r = Reader(sys.argv[1])
for i, ev in enumerate(r.get('events')):
    print('event', i)
    for name in ev.getAvailableCollections():
        c = ev.get(name)
        tn = c.getTypeName() if hasattr(c, 'getTypeName') else type(c).__name__
        if len(c) and 'Calorimeter' in str(tn):
            e = sum(h.getEnergy() for h in c)
            print(f'  {name:30s} {str(tn):40s} n={len(c):6d} Esum={e:.4f} GeV')
        elif len(c):
            print(f'  {name:30s} {str(tn):40s} n={len(c):6d}')
