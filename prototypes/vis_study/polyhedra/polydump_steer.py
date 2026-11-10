# ddsim steering add-on (prototype): the SoLID steering + a geometry tessellation dump at /run/initialize.
#   SOLID_OPTICAL=0 POLYDUMP_OUT=out/hook ddsim --steeringFile polydump_steer.py -N 2 --outputFile out/hook.edm4hep.root
import os
import sys
_steer = os.path.join(os.environ['DETECTOR_PATH'], 'steering', 'solid_sim.py')
exec(compile(open(_steer).read(), _steer, 'exec'))
sys.path.insert(0, '$SCRATCH/vis_study/polyhedra')
import g4poly  # noqa: E402
g4poly.install(SIM, os.environ.get('POLYDUMP_OUT', 'polydump'))  # noqa: F821
