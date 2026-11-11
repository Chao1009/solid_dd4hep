# SPDX-License-Identifier: LGPL-3.0-or-later
"""Reader of the trajectory file written by the SoLIDTrajectoryRecorder DDG4 action (SOLID_TRAJ in
steering/solid_sim.py; source src/sim/TrajectoryRecorder.cpp).

The file is little endian: a 40-byte header, then 40-byte records, one per trajectory point::

    header   char[8] "SOLIDTRJ", int32 version (1), int32 record size (40), float32 MinEkin [MeV],
             float32 MaxDistance [mm], float32 MaxAngle [rad], 12 bytes reserved
    record   int32 event, track, parent, pdg; float32 x, y, z [mm], t [ns], ekin [MeV]; int32 flags
             (bit 0: first point of the track, bit 1: last point)

so ``numpy.fromfile(path, TRAJ_DTYPE, offset=40)`` is the whole content.  The points of one event are
contiguous; events are in the order they finished (not sorted in a multi-threaded run).
"""
import numpy as np

HEADER_DTYPE = np.dtype([("magic", "S8"), ("version", "<i4"), ("record_size", "<i4"), ("min_ekin", "<f4"),
                         ("max_distance", "<f4"), ("max_angle", "<f4"), ("reserved", "S12")])
TRAJ_DTYPE = np.dtype([("event", "<i4"), ("track", "<i4"), ("parent", "<i4"), ("pdg", "<i4"),
                       ("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("t", "<f4"), ("ekin", "<f4"), ("flags", "<i4")])
FLAG_FIRST, FLAG_LAST = 1, 2


def read_header(path):
    """The header as a dict: version, record_size, min_ekin [MeV], max_distance [mm], max_angle [rad]."""
    h = np.fromfile(path, HEADER_DTYPE, count=1)
    if len(h) != 1 or h["magic"][0] != b"SOLIDTRJ":
        raise ValueError("%s is not a SoLID trajectory file (no SOLIDTRJ header)" % path)
    h = h[0]
    if int(h["version"]) != 1 or int(h["record_size"]) != TRAJ_DTYPE.itemsize:
        raise ValueError("%s: trajectory file version %d with %d-byte records is not supported" % (
            path, int(h["version"]), int(h["record_size"])))
    return dict(version=int(h["version"]), record_size=int(h["record_size"]), min_ekin=float(h["min_ekin"]),
                max_distance=float(h["max_distance"]), max_angle=float(h["max_angle"]))


def read_trajectories(path, event=None):
    """All trajectory points as a structured array (TRAJ_DTYPE), or those of one event.

    A file that ends in an incomplete record (a job that was killed) is read up to the last complete one.
    """
    read_header(path)
    pts = np.fromfile(path, TRAJ_DTYPE, offset=HEADER_DTYPE.itemsize)
    if event is not None:
        pts = pts[pts["event"] == event]
    return pts


def split_events(points):
    """{event number: points of that event}, each in file order."""
    order = np.argsort(points["event"], kind="stable")
    p = points[order]
    ev, first = np.unique(p["event"], return_index=True)
    bounds = list(first) + [len(p)]
    return {int(e): p[bounds[i]:bounds[i + 1]] for i, e in enumerate(ev)}


def tracks_of(points):
    """Arrange the points of one event by track.

    Returns (trk, offset, pts): trk is a structured array with one row per track (track, parent, pdg, ekin
    at the first point), pts the points sorted by track and, within a track, in file order (= time order);
    track k owns pts[offset[k]:offset[k + 1]].
    """
    order = np.argsort(points["track"], kind="stable")
    p = points[order]
    tid, first = np.unique(p["track"], return_index=True)
    offset = np.append(first, len(p)).astype(np.int64)
    trk = np.zeros(len(tid), np.dtype([("track", "<i4"), ("parent", "<i4"), ("pdg", "<i4"), ("ekin", "<f4")]))
    trk["track"] = tid
    for k in ("parent", "pdg", "ekin"):
        trk[k] = p[k][first]
    return trk, offset, p
