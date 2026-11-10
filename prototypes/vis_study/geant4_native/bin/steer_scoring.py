# Wrapper steering: steer_xpos.py + instantiate G4ScoringManager (Geant4 command-based scoring, /score/...)
# ddsim never calls G4ScoringManager::GetScoringManager(), so the /score/ commands do not exist by default.
import os
exec(open('$SCRATCH/vis_study/geant4_native/bin/steer_xpos.py').read())
import ROOT
_g4inc = os.popen('geant4-config --cflags').read().split()
for tok in _g4inc:
    if tok.startswith('-I'):
        ROOT.gInterpreter.AddIncludePath(tok[2:])
ROOT.gSystem.Load('libG4digits_hits')
ROOT.gInterpreter.Declare('#include "G4ScoringManager.hh"')
_scm = ROOT.G4ScoringManager.GetScoringManager()
print('steer_scoring: G4ScoringManager instance', _scm)
