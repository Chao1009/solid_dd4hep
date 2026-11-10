# SPDX-License-Identifier: LGPL-3.0-or-later
"""Hit cellID -> drawn volume of the geometry scene (scene.npz).

A drawn volume of the scene carries the readout (hit collection) it belongs to and the volume ID bits that
the placements down to it fix, with their mask.  A hit of that readout lies in the volume if
``(cellID & mask) == id``.  Several nested volumes can match (a calorimeter module and the plate it stands
on).  The rule used everywhere in this repository:

    the hit belongs to the matching drawn volume with the most mask bits (the most specific one);
    no match and two matches with the same number of bits (a tie) are errors.

Candidates are all drawn volumes of the hit's readout that fix at least one ID bit.  `sensitive_only` narrows
them to sensitive volumes; the SoLID geometries do not need that for uniqueness (checked by
scripts/geometry/check_scene_hits.py and by solid_evexport on every run).
"""
import numpy as np

NO_MATCH, TIE = -1, -2


def popcount(x):
    """number of set bits of each uint64"""
    x = np.asarray(x, np.uint64)
    n = np.zeros(x.shape, np.int64)
    for k in range(64):
        n += ((x >> np.uint64(k)) & np.uint64(1)).astype(np.int64)
    return n


class VolumeMap:
    """The volume-ID table of a scene, for mapping hits to scene nodes.

    scene: path of a scene.npz, or a mapping with its arrays (numpy.load result, solid_lineview SceneData).
    """

    def __init__(self, scene, sensitive_only=False):
        z = np.load(scene, allow_pickle=False) if isinstance(scene, (str, bytes)) or hasattr(scene, "__fspath__") \
            else scene
        self.readout = np.asarray(z["node_readout"]).astype(str)
        self.vid = np.asarray(z["node_volume_id"]).astype(np.uint64)
        self.mask = np.asarray(z["node_volume_mask"]).astype(np.uint64)
        cand = (np.asarray(z["node_mesh"]) >= 0) & (self.mask != 0)
        if sensitive_only:
            cand &= np.asarray(z["node_sensitive"]).astype(bool)
        self.candidate = cand
        self._tables = {}

    def readouts(self):
        """readout names that have candidate volumes"""
        return sorted(set(self.readout[self.candidate]) - {""})

    def _table(self, readout):
        """[(mask, bits, sorted ids, node of each id)] of a readout, most mask bits first"""
        t = self._tables.get(readout)
        if t is None:
            nodes = np.nonzero(self.candidate & (self.readout == readout))[0]
            t = []
            for m in np.unique(self.mask[nodes]):
                sel = nodes[self.mask[nodes] == m]
                order = np.argsort(self.vid[sel], kind="stable")
                t.append((m, int(popcount(m)), self.vid[sel][order], sel[order]))
            t.sort(key=lambda r: -r[1])
            self._tables[readout] = t
        return t

    def match(self, readout, cellid):
        """Map the cellIDs of hits of one readout.

        Returns (node, nmatch): node is the scene node index of the most specific drawn volume, NO_MATCH (-1)
        if no volume matches, TIE (-2) if several share the largest number of mask bits; nmatch is the
        number of matching volumes of any depth.
        """
        cid = np.asarray(cellid).astype(np.uint64)
        node = np.full(len(cid), NO_MATCH, np.int64)
        nmatch = np.zeros(len(cid), np.int64)
        best = np.full(len(cid), -1, np.int64)
        nbest = np.zeros(len(cid), np.int64)
        for m, bits, ids, nodes in self._table(readout):
            key = cid & m
            lo = np.searchsorted(ids, key, side="left")
            cnt = np.searchsorted(ids, key, side="right") - lo
            hit = cnt > 0
            nmatch += cnt
            better = hit & (bits > best)
            same = hit & (bits == best)
            node[better] = nodes[lo[better]]
            best[better] = bits
            nbest[better] = cnt[better]
            nbest[same] += cnt[same]
        node[nbest > 1] = TIE
        return node, nmatch
