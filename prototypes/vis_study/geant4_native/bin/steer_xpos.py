# Wrapper steering: SoLID default steering with a restricted gun (env overrides) and a fixed seed.
#   GUN_PHIMIN/GUN_PHIMAX (deg, default -50/50 = x>0 half), GUN_THETAMIN/GUN_THETAMAX (deg, default from solid_sim.py),
#   SEED (default 4242), NEVT (batch events)
import os
exec(open(os.path.join(os.environ['DETECTOR_PATH'], 'steering', 'solid_sim.py')).read())
from g4units import deg
SIM.gun.phiMin = float(os.environ.get('GUN_PHIMIN', '-50')) * deg
SIM.gun.phiMax = float(os.environ.get('GUN_PHIMAX', '50')) * deg
if 'GUN_THETAMIN' in os.environ:
    SIM.gun.thetaMin = float(os.environ['GUN_THETAMIN']) * deg
if 'GUN_THETAMAX' in os.environ:
    SIM.gun.thetaMax = float(os.environ['GUN_THETAMAX']) * deg
SIM.random.seed = int(os.environ.get('SEED', '4242'))
SIM.numberOfEvents = int(os.environ.get('NEVT', '3'))
