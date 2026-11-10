// SolidStepRecorder: DDG4 stepping action that writes thinned trajectory points with global time
// to a flat binary side file (feasibility prototype for SoLID event displays / animations).
//
// Record (little endian, 48 bytes):
//   int32 event, track, parent, pdg, flags   (flags: 1 = first point of the track, 2 = last point)
//   float32 x, y, z [mm], t [ns], ekin [MeV], edep [MeV] (energy deposited since the previous point)
//   int32 volume code (0, reserved)
// Thinning: a post-step point is written when the track moved more than MaxDistance since the last
// written point, turned by more than MaxAngle, or ended. Tracks with vertex kinetic energy below MinEkin
// and optical photons are skipped.
//
// ddsim usage:
//   --action.step '{"name":"SolidStepRecorder/StepRec","parameter":{"OutputFile":"steps.bin"}}'
#include "DDG4/Factories.h"
#include "DDG4/Geant4SteppingAction.h"
#include "DD4hep/InstanceCount.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>

namespace dd4hep {
namespace sim {

class SolidStepRecorder : public Geant4SteppingAction {
  std::string m_file{"steps.bin"};
  double m_maxDist{20.0};   // mm
  double m_maxAngle{0.04};  // rad
  double m_minEkin{0.0};    // MeV
  FILE* m_fp{nullptr};
  int m_curEvent{-1}, m_curTrack{-1};
  bool m_skip{false};
  G4ThreeVector m_lastPos, m_lastDir;
  double m_edep{0};
  long m_nrec{0}, m_nsteps{0};

#pragma pack(push, 1)
  struct Rec {
    int32_t event, track, parent, pdg, flags;
    float x, y, z, t, ekin, edep;
    int32_t vol;
  };
#pragma pack(pop)
  static_assert(sizeof(Rec) == 48, "record size");

  void put(const G4Track* trk, int evt, const G4ThreeVector& p, double t, double ekin, int flags) {
    if (!m_fp) {
      m_fp = std::fopen(m_file.c_str(), "wb");
      if (!m_fp) except("cannot open %s", m_file.c_str());
    }
    Rec r{evt, trk->GetTrackID(), trk->GetParentID(), trk->GetDefinition()->GetPDGEncoding(), flags,
          float(p.x()), float(p.y()), float(p.z()), float(t), float(ekin), float(m_edep), 0};
    std::fwrite(&r, sizeof(r), 1, m_fp);
    m_edep = 0;
    ++m_nrec;
  }

public:
  SolidStepRecorder(Geant4Context* ctxt, const std::string& nam) : Geant4SteppingAction(ctxt, nam) {
    declareProperty("OutputFile", m_file);
    declareProperty("MaxDistance", m_maxDist);
    declareProperty("MaxAngle", m_maxAngle);
    declareProperty("MinEkin", m_minEkin);
    InstanceCount::increment(this);
  }
  virtual ~SolidStepRecorder() {
    if (m_fp) std::fclose(m_fp);
    printout(INFO, name(), "%ld steps seen, %ld points written to %s", m_nsteps, m_nrec, m_file.c_str());
    InstanceCount::decrement(this);
  }
  void operator()(const G4Step* step, G4SteppingManager*) override {
    const G4Track* trk = step->GetTrack();
    if (trk->GetDefinition() == G4OpticalPhoton::Definition()) return;
    ++m_nsteps;
    const G4Event* ev = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    const int evt = ev ? ev->GetEventID() : -1;
    const G4StepPoint* pre = step->GetPreStepPoint();
    const G4StepPoint* post = step->GetPostStepPoint();
    if (trk->GetCurrentStepNumber() == 1 || trk->GetTrackID() != m_curTrack || evt != m_curEvent) {
      m_curTrack = trk->GetTrackID();
      m_curEvent = evt;
      m_skip = trk->GetVertexKineticEnergy() < m_minEkin;
      m_edep = 0;
      if (m_skip) return;
      if (trk->GetCurrentStepNumber() == 1)
        put(trk, evt, pre->GetPosition(), pre->GetGlobalTime(), pre->GetKineticEnergy(), 1);
      m_lastPos = pre->GetPosition();
      m_lastDir = pre->GetMomentumDirection();
    }
    if (m_skip) return;
    m_edep += step->GetTotalEnergyDeposit();
    const G4ThreeVector& p = post->GetPosition();
    const G4ThreeVector& d = post->GetMomentumDirection();
    const bool last = trk->GetTrackStatus() != fAlive;
    const double dot = std::max(-1.0, std::min(1.0, d.dot(m_lastDir)));
    if (last || (p - m_lastPos).mag() > m_maxDist || std::acos(dot) > m_maxAngle) {
      put(trk, evt, p, post->GetGlobalTime(), post->GetKineticEnergy(), last ? 2 : 0);
      m_lastPos = p;
      m_lastDir = d;
    }
  }
};

}  // namespace sim
}  // namespace dd4hep

DECLARE_GEANT4ACTION(SolidStepRecorder)
