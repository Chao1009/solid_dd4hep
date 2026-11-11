// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// Plugins around the scene producer (SceneDump.h):
//
//   SoLIDSceneDump    DDG4 detector-construction action.  It writes the scene when Geant4 initialises the run,
//                     after DDG4 has built the geometry and the sensitive detectors.  steering/solid_sim.py adds
//                     it when SOLID_SCENE=<file.npz> is set.  Properties: Output (.npz), Prim (.prim, optional),
//                     Segments (line segments per circle, default 72), Config (name stored in the file).
//
//   SoLID_SceneDump   DD4hep apply plugin for use without ddsim: converts the loaded geometry to Geant4 (no run
//                     manager, no physics) and writes the scene.  bin/solid_scene_dump calls it.  Arguments:
//                     output=<file.npz> prim=<file.prim> segments=<n> config=<name> ids=0|1

#include "SceneDump.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Printout.h>
#include <DDG4/Factories.h>
#include <DDG4/Geant4DetectorConstruction.h>

#include <atomic>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace solid {

class SoLIDSceneDump : public dd4hep::sim::Geant4DetectorConstruction {
public:
  SoLIDSceneDump(dd4hep::sim::Geant4Context* context, const std::string& name)
    : dd4hep::sim::Geant4DetectorConstruction(context, name)
  {
    declareProperty("Output", m_output);
    declareProperty("Prim", m_prim);
    declareProperty("Segments", m_segments = 72);
    declareProperty("Config", m_config);
  }

  // Called from ConstructSDandField(), once per thread: the geometry is shared, write it once.
  void constructSensitives(dd4hep::sim::Geant4DetectorConstructionContext* ctxt) override
  {
    static std::atomic<bool> done{false};
    if (done.exchange(true) || (m_output.empty() && m_prim.empty())) {
      return;
    }
    scene::Options options;
    options.segments = m_segments;
    options.config   = m_config;
    options.mode     = "ddsim";
    try {
      scene::dump(ctxt->description, ctxt->world, m_output, m_prim, options);
    }
    catch (const std::exception& e) {
      // a drawing aid must not stop a simulation job
      error("scene not written: %s", e.what());
    }
  }

private:
  std::string m_output;
  std::string m_prim;
  std::string m_config;
  int         m_segments = 72;
};

} // namespace solid

DECLARE_GEANT4ACTION_NS(solid, SoLIDSceneDump)

static long scene_dump(dd4hep::Detector& description, int argc, char** argv)
{
  solid::scene::Options options;
  std::string           npz, prim;
  for (int i = 0; i < argc && argv[i] != nullptr; ++i) {
    const std::string arg = argv[i];
    const auto        eq  = arg.find('=');
    const std::string key = arg.substr(0, eq);
    const std::string val = eq == std::string::npos ? "" : arg.substr(eq + 1);
    if (key == "output") {
      npz = val;
    }
    else if (key == "prim") {
      prim = val;
    }
    else if (key == "segments") {
      options.segments = std::atoi(val.c_str());
    }
    else if (key == "config") {
      options.config = val;
    }
    else if (key == "ids") {
      options.ids = val != "0";
    }
    else {
      dd4hep::except("SoLID_SceneDump", "unknown argument '%s' (output= prim= segments= config= ids=)", arg.c_str());
    }
  }
  if (npz.empty() && prim.empty()) {
    dd4hep::except("SoLID_SceneDump", "nothing to write: give output=<file.npz> and/or prim=<file.prim>");
  }
  if (options.segments < 3) {
    dd4hep::except("SoLID_SceneDump", "segments must be at least 3");
  }
  solid::scene::dump(description, nullptr, npz, prim, options);
  return 1;
}
DECLARE_APPLY(SoLID_SceneDump, scene_dump)
