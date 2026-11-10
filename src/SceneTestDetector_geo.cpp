// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_SceneTestDetector: a small detector for testing the drawing scene (test/test_scene_dump.xml, ctest
// "scene_dump"). It is not part of any SoLID configuration. It contains the kinds of placement that the SoLID
// detectors do not use, so that the scene producer is exercised on them:
//
//   pad      a sensitive box placed twice (one logical volume, two placements), layer 1, module 1 and 2
//   wedge    a sensitive extruded right triangle (no mirror symmetry in x), layer 2:
//            module 3 placed with a rotation, module 4 placed with a reflection x -> -x
//   strip    a sensitive box placed as a parameterised volume inside an invisible frame, `nstrip` copies along y,
//            layer 3 (set on the frame), ID field strip = copy number
// Replicas (Volume::replicate) are left out: DDG4 1.38 passes their width and offset to Geant4 in cm instead of mm.
//
// The envelope is an invisible box whose centre is at (0, 0, z).
//
// Compact:
//   <detector id="1" name="SceneTest" type="SoLID_SceneTestDetector" readout="SceneTestHits" z="100*mm"
//             nstrip="4" vis="..." material="..."/>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <string>
#include <vector>

using namespace dd4hep;

static Ref_t create_scene_test(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t         x_det    = e;
  const std::string name     = x_det.nameStr();
  const double      z        = x_det.attr<double>(_Unicode(z));
  const int         nstrip   = x_det.hasAttr(_Unicode(nstrip)) ? x_det.attr<int>(_Unicode(nstrip)) : 4;
  const Material    mat      = desc.material(x_det.attr<std::string>(_Unicode(material)));
  const Material    air      = desc.air();

  DetElement det(name, x_det.id());
  sens.setType("tracker");

  Volume envelope(name + "_envelope", Box(120 * mm, 120 * mm, 20 * mm), air);
  envelope.setVisAttributes(desc.invisible());

  // 1. one volume, two placements
  Volume pad(name + "_pad", Box(10 * mm, 20 * mm, 5 * mm), mat);
  pad.setVisAttributes(desc, x_det.visStr());
  pad.setSensitiveDetector(sens);
  for (int i = 0; i < 2; ++i) {
    PlacedVolume pv = envelope.placeVolume(pad, i, Position((i == 0 ? -50 : 50) * mm, 0, 0));
    pv.addPhysVolID("layer", 1).addPhysVolID("module", i + 1);
  }

  // 2. a solid without mirror symmetry, placed rotated and reflected
  const std::vector<double> px{0 * mm, 30 * mm, 0 * mm};
  const std::vector<double> py{0 * mm, 0 * mm, 10 * mm};
  Volume wedge(name + "_wedge", ExtrudedPolygon(px, py, {-5 * mm, 5 * mm}, {0, 0}, {0, 0}, {1, 1}), mat);
  wedge.setVisAttributes(desc, x_det.visStr());
  wedge.setSensitiveDetector(sens);
  {
    PlacedVolume pv = envelope.placeVolume(wedge, 0, Transform3D(RotationZ(0), Position(20 * mm, 60 * mm, 0)));
    pv.addPhysVolID("layer", 2).addPhysVolID("module", 3);
    const Rotation3D mirror_x(-1, 0, 0, 0, 1, 0, 0, 0, 1); // determinant -1: a reflection
    pv = envelope.placeVolume(wedge, 1, Transform3D(mirror_x, Position(-20 * mm, 60 * mm, 0)));
    pv.addPhysVolID("layer", 2).addPhysVolID("module", 4);
  }

  // 3. parameterised volume
  Volume strip(name + "_strip", Box(15 * mm, 2 * mm, 2 * mm), mat);
  strip.setVisAttributes(desc, x_det.visStr());
  strip.setSensitiveDetector(sens);
  {
    // Geant4 wants a parameterised volume to be the only daughter of its mother
    Volume frame(name + "_strip_frame", Box(20 * mm, 20 * mm, 3 * mm), air);
    frame.setVisAttributes(desc.invisible());
    PlacedVolume pv = frame.paramVolume1D(Transform3D(Position(0, -12 * mm, 0)), strip,
                                          static_cast<std::size_t>(nstrip), Transform3D(Position(0, 8 * mm, 0)));
    pv.addPhysVolID("strip", 0); // DD4hep replaces the 0 by the copy number of each strip
    pv = envelope.placeVolume(frame, Position(0, -70 * mm, 0));
    pv.addPhysVolID("layer", 3);
  }

  PlacedVolume pv = desc.pickMotherVolume(det).placeVolume(envelope, Position(0, 0, z));
  pv.addPhysVolID("system", x_det.id());
  det.setPlacement(pv);
  return det;
}

DECLARE_DETELEMENT(SoLID_SceneTestDetector, create_scene_test)
