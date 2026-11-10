"""SoLID steering + Firebird trajectory writer (feasibility prototype).

ddsim --steeringFile solid_firebird_steering.py --compactFile $DETECTOR_PATH/solid_JPsi_lite.xml \
      -N 3 --outputFile out.edm4hep.root      # FB_TRAJ_OUT=out.traj.firebird.json (Geant4 rich trajectories)
"""
import os

exec(open(os.path.join(os.environ['DETECTOR_PATH'], 'steering', 'solid_sim.py')).read())

import DDG4  # noqa: E402

_out = os.environ.get('FB_TRAJ_OUT', 'trajectories.firebird.json')
SIM.ui.commandsConfigure = ['/tracking/storeTrajectory 3']   # G4RichTrajectory: points carry PreT/PostT times
_kernel = DDG4.Kernel()
_writer = DDG4.EventAction(_kernel, 'FirebirdTrajectoryWriterEventAction/TrajectoryWriter')
_writer.ComponentName = 'Geant4Trajectories'
_writer.OutputFile = _out
_writer.OnlyPrimary = False
_writer.MomentumMin = float(os.environ.get('FB_PMIN_MEV', '30'))   # MeV/c
_writer.StepCut = True             # drop points outside the SoLID-Lite hall region
_writer.StepZMin = -4000.0
_writer.StepZMax = 6000.0
_writer.StepRMax = 3500.0
_kernel.eventAction().add(_writer)
