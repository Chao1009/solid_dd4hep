#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Count simulated hits per readout collection in an EDM4hep file (smoke test for every sensitive detector).

usage: count_hits.py sim.edm4hep.root [--require NAME ...]   exit 1 if a required collection has no hits
"""
import argparse
import sys

import uproot


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('file')
    ap.add_argument('--require', nargs='*', default=[])
    a = ap.parse_args()
    tree = uproot.open(a.file)['events']
    n_ev = tree.num_entries
    counts = {}
    for key in tree.keys(recursive=False):
        if key.startswith('_') or not key.endswith('Hits'):
            continue
        branch = f'{key}/{key}.cellID'
        if branch not in tree:
            continue
        arr = tree[branch].array(library='np')
        counts[key] = (sum(len(x) for x in arr), sum(1 for x in arr if len(x) > 0))
    print(f'{a.file}: {n_ev} events')
    for k, (nh, ne) in sorted(counts.items()):
        print(f'  {k:28s} {nh:8d} hits in {ne:3d}/{n_ev} events')
    missing = [r for r in a.require if counts.get(r, (0, 0))[0] == 0]
    if missing:
        print('MISSING hits in:', ', '.join(missing))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
