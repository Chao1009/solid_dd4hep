// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// Geometry scene for drawings: the polygon meshes Geant4 would hand to a vis driver, taken straight from the
// Geant4 volume tree that DDG4 built, without a vis manager.
//
// The walk follows G4PhysicalVolumeModel: depth first, global transforms, replicas and parameterised volumes
// expanded, volumes drawn when their vis attributes are visible (and alpha > 0), descent stopped below volumes
// with "daughters invisible".  Meshes are G4VSolid::GetPolyhedron() with a chosen number of segments per circle,
// one mesh per distinct solid.  Each placement also gets its DD4hep DetElement path, readout and the part of the
// volume ID fixed by the placements above it (value and bit mask): a hit belongs to the placement when
// (cellID & mask) == id.
//
// Output: scene.npz (schema v1, see docs/visualization.md) and, optionally, a DAWN g4.prim file.

#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

class G4VPhysicalVolume;
namespace dd4hep {
class Detector;
}

namespace solid::scene {

inline constexpr int SCHEMA_VERSION = 1;

struct Mesh {
  std::string               name;   // name of the Geant4 solid
  std::string               type;   // G4VSolid::GetEntityType()
  std::vector<double>       v;      // 3 per vertex, frame of the solid [mm]
  std::vector<std::int32_t> f;      // 4 per facet: 1-based vertex index, negative = Geant4 marks the edge that
                                    // starts at this vertex as invisible, 0 = no 4th vertex
};

struct Node {
  std::string            name;                 // name of the Geant4 physical volume
  std::int32_t           copy   = 0;           // copy number (replica number for replicas and parameterised volumes)
  std::int32_t           depth  = 0;           // 0 = world
  std::int32_t           parent = -1;          // index of the mother node
  std::int32_t           mesh   = -1;          // index into Scene::meshes, -1 = not drawn
  std::int8_t            kind   = 0;           // 0 placement, 1 replica, 2 parameterised
  std::array<double, 12> matrix{1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}; // local -> world: R (3x3, row-major), T [mm]
  std::array<double, 4>  rgba{1, 1, 1, 1};
  bool                   visible             = true;  // vis attribute flag of the logical volume
  bool                   daughters_invisible = false;
  bool                   wireframe           = false; // Geant4 forced wireframe (DDG4 sets it from lineStyle)
  std::int8_t            drawing_style       = 0;     // DD4hep VisAttr drawingStyle: 0 unknown, 1 solid, 2 wireframe
  std::string            det_path;                    // deepest DetElement above or at this volume
  std::string            readout;                     // readout (hit collection) name, empty if none
  std::uint64_t          volume_id   = 0;             // volume ID bits fixed down to this placement
  std::uint64_t          volume_mask = 0;             // mask of those bits
  bool                   sensitive   = false;
};

struct Scene {
  std::vector<Mesh>                  meshes;
  std::vector<Node>                  nodes;                    // depth-first order
  std::string                        config;                   // configuration name (compact file stem)
  std::string                        mode;                     // "standalone" or "ddsim"
  int                                segments = 72;            // line segments per circle
  std::array<double, 4>              bounding_sphere{0, 0, 0, 0}; // centre x, y, z and radius of the drawn volumes
  std::map<std::string, std::string> readouts;                 // readout name -> ID field description
  std::string                        geometry_hash;            // SHA-256 over the rounded scene content
  int                                failed_meshes = 0;        // drawn volumes whose solid gave no polyhedron
};

struct Options {
  int         segments = 72;
  bool        ids      = true; // DD4hep DetElement paths, readouts and volume IDs
  std::string config;
  std::string mode = "standalone";
};

// Convert the DD4hep geometry to Geant4 the way DDG4 does, without a run manager (standalone use).
G4VPhysicalVolume* convertGeometry(dd4hep::Detector& description);

// Walk the Geant4 tree below `world` (DDG4's current geometry) and build the scene.
Scene collect(dd4hep::Detector& description, const G4VPhysicalVolume* world, const Options& options);

void writeNpz(const Scene& scene, const std::string& path);

// DAWN g4.prim (format 2.4): one /Polyhedron per drawn placement, all edges visible, as Geant4's DAWNFILE driver
// writes them.  For dawn and dawncut; solid_lineview reads the .npz directly.
void writePrim(const Scene& scene, const std::string& path);

// collect() and write the files whose names are not empty; prints one summary line.  A null `world` means
// standalone use: the geometry is converted with convertGeometry() first.
Scene dump(dd4hep::Detector& description, const G4VPhysicalVolume* world, const std::string& npz, const std::string& prim,
           const Options& options);

} // namespace solid::scene
