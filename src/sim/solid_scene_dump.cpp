// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// solid_scene_dump: write the drawing scene (scene.npz, optionally a DAWN .prim) of a DD4hep compact file.
// Loads the geometry, converts it to Geant4 and calls the SoLID_SceneDump plugin; no run manager, no physics.

#include <DD4hep/Detector.h>
#include <DD4hep/Printout.h>

#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

const char* const USAGE = R"(usage: solid_scene_dump [COMPACT.xml] [-o SCENE.npz] [--prim FILE.prim] [--segments N] [--no-ids] [-v]

Write the geometry scene that solid_lineview and make_solid_views draw: the Geant4 polyhedra of all visible
volumes with names, colours, transforms, DetElement paths and DD4hep volume IDs.

  COMPACT.xml     compact file (default: $DETECTOR_PATH/$DETECTOR_CONFIG.xml, set by thissolid.sh)
  -o, --output    scene file to write (default: <configuration>.scene.npz unless --prim is given)
  --prim FILE     also (or only) write a DAWN g4.prim file, for dawn and dawncut
  --segments N    line segments per full circle (default 72)
  --no-ids        skip DetElement paths, readouts and volume IDs
  -v, --verbose   print DD4hep messages of level INFO
  -h, --help      this text
)";

}

int main(int argc, char** argv)
{
  std::string compact, output, prim;
  int         segments = 72;
  bool        ids = true, verbose = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto value = [&](const char* what) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "solid_scene_dump: " << what << " needs a value\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "-h" || a == "--help") {
      std::cout << USAGE;
      return 0;
    }
    else if (a == "-o" || a == "--output") {
      output = value("--output");
    }
    else if (a == "--prim") {
      prim = value("--prim");
    }
    else if (a == "--segments") {
      segments = std::atoi(value("--segments").c_str());
    }
    else if (a == "--no-ids") {
      ids = false;
    }
    else if (a == "-v" || a == "--verbose") {
      verbose = true;
    }
    else if (!a.empty() && a[0] == '-') {
      std::cerr << "solid_scene_dump: unknown option " << a << "\n" << USAGE;
      return 2;
    }
    else if (compact.empty()) {
      compact = a;
    }
    else {
      std::cerr << "solid_scene_dump: more than one compact file given\n";
      return 2;
    }
  }
  if (compact.empty()) {
    const char* path   = std::getenv("DETECTOR_PATH");
    const char* config = std::getenv("DETECTOR_CONFIG");
    if (path == nullptr || config == nullptr) {
      std::cerr << "solid_scene_dump: no compact file given and DETECTOR_PATH/DETECTOR_CONFIG are not set\n" << USAGE;
      return 2;
    }
    compact = std::string(path) + "/" + config + ".xml";
  }
  if (segments < 3) {
    std::cerr << "solid_scene_dump: --segments must be at least 3\n";
    return 2;
  }
  const std::string config = std::filesystem::path(compact).stem().string();
  if (output.empty() && prim.empty()) {
    output = config + ".scene.npz";
  }

  try {
    dd4hep::setPrintLevel(verbose ? dd4hep::INFO : dd4hep::WARNING);
    dd4hep::Detector& description = dd4hep::Detector::getInstance();
    description.fromXML(compact);
    std::vector<std::string> args{"config=" + config, "segments=" + std::to_string(segments),
                                  std::string("ids=") + (ids ? "1" : "0")};
    if (!output.empty()) {
      args.push_back("output=" + output);
    }
    if (!prim.empty()) {
      args.push_back("prim=" + prim);
    }
    std::vector<char*> cargs;
    for (auto& s : args) {
      cargs.push_back(s.data());
    }
    dd4hep::setPrintLevel(dd4hep::INFO); // the one summary line of the plugin
    description.apply("SoLID_SceneDump", static_cast<int>(cargs.size()), cargs.data());
  }
  catch (const std::exception& e) {
    std::cerr << "solid_scene_dump: " << e.what() << "\n";
    return 1;
  }
  // Leave without running the destructors of the DD4hep and Geant4 geometry stores: nothing is left to save,
  // and tearing down two large geometries takes longer than writing the scene.
  std::cout.flush();
  std::cerr.flush();
  std::_Exit(0);
}
