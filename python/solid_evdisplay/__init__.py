# SPDX-License-Identifier: LGPL-3.0-or-later
"""solid_evdisplay: event displays and animations on the SoLID line drawings.

The package turns a ddsim run into event figures in three steps::

    ddsim with SOLID_SCENE=scene.npz SOLID_TRAJ=traj.bin   -> sim.edm4hep.root, scene.npz, traj.bin
    solid_evexport sim.edm4hep.root --traj traj.bin --scene scene.npz -o events/   -> events/event_000000.npz ...
    solid_evview events/ --scene scene.npz --views side,beam,cutaway -o figures/  -> SVG, PNG, PDF, GIF ...

Modules:

    trajfile   reader of the trajectory file of the SoLIDTrajectoryRecorder DDG4 action
    hitmap     hit cellID -> most specific drawn volume of the geometry scene
    export     EDM4hep + trajectories + scene -> one event.npz per event (the ``solid_evexport`` command)
    event      event.npz reader, event selection, overlay of several events
    view       the views and their drawing on the solid_lineview renderer (the ``solid_evview`` command)
    animate    frames with cached geometry layers, GIF / WebP / animated SVG / ffmpeg output

From Python::

    from solid_evdisplay import load_event, VolumeMap, read_trajectories
    ev = load_event('events/event_000003.npz')
    ev.vol_energy                      # MeV per drawn volume, ev.vol_node are scene node indices
"""
__version__ = "0.1.0"

from .trajfile import TRAJ_DTYPE, read_header, read_trajectories  # noqa: F401,E402
from .hitmap import VolumeMap, popcount  # noqa: F401,E402
from .event import Event, load_event, merge_events, select_events  # noqa: F401,E402
