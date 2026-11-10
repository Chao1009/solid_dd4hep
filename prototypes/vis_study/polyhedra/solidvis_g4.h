// Prototype: tessellate the Geant4 geometry built by DDG4 without any vis driver.
// Declared into cling from Python (g4poly.py).  Mirrors what G4PhysicalVolumeModel +
// G4VSceneHandler::RequestPrimitives + the DAWNFILE scene handler do:
//   - depth-first walk of the G4 physical-volume tree, global transforms
//   - culling of invisible volumes (vis attributes), daughters-invisible stops the descent
//   - G4VSolid::GetPolyhedron() with G4Polyhedron::SetNumberOfRotationSteps(nseg)
//   - facets as 1-based vertex indices multiplied by the edge flag (DAWN /Facet convention)
#pragma once
#include <vector>
#include <string>
#include <map>
#include <chrono>

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVReplica.hh"
#include "G4ReplicaNavigation.hh"
#include "G4Material.hh"
#include "G4VisExtent.hh"
#include "G4Navigator.hh"
#include "G4TouchableHistory.hh"

#include "DD4hep/Detector.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4GeometryInfo.h"
#include "DDG4/Geant4VolumeManager.h"

namespace solidvis {

struct Mesh {
  std::vector<double> v;   // 3*nv, local frame of the solid [mm]
  std::vector<int>    f;   // 4*nf, 1-based vertex index * edge flag (+ visible, - invisible), 0 = no 4th vertex
  int  nv = 0, nf = 0;
  bool ok = false;
  double seconds = 0;
};

// Tessellate a solid exactly like the vis drivers do (cached polyhedron of the solid).
inline Mesh tessellate(const G4VSolid* solid, int nseg) {
  Mesh r;
  auto t0 = std::chrono::steady_clock::now();
  G4Polyhedron::SetNumberOfRotationSteps(nseg);
  G4Polyhedron* p = solid->GetPolyhedron();     // owned by the solid
  G4Polyhedron::ResetNumberOfRotationSteps();
  if (p && p->GetNoFacets() > 0) {
    r.nv = p->GetNoVertices();
    r.nf = p->GetNoFacets();
    r.v.reserve(3 * r.nv);
    r.f.reserve(4 * r.nf);
    for (int i = 1; i <= r.nv; ++i) {
      G4Point3D q = p->GetVertex(i);
      r.v.push_back(q.x()); r.v.push_back(q.y()); r.v.push_back(q.z());
    }
    for (int i = 1; i <= r.nf; ++i) {
      G4int n = 0, nodes[4] = {0, 0, 0, 0}, ef[4] = {1, 1, 1, 1};
      p->GetFacet(i, n, nodes, ef);
      for (int k = 0; k < 4; ++k) r.f.push_back(k < n ? nodes[k] * (ef[k] > 0 ? 1 : -1) : 0);
    }
    r.ok = true;
  }
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

struct Node {
  const G4VPhysicalVolume* pv = nullptr;
  const G4VSolid* solid = nullptr;
  int depth = 0, copy = 0, parent = -1;
  int kind = 0;               // 0 placement, 1 replica, 2 parameterised
  bool visible = true, has_vis = false, daughters_invisible = false, wire = false;
  double rgba[4] = {1, 1, 1, 1};
  double m[12] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};   // row-major R (3x3), then T
  double density = 0;         // g/cm3
};

class Walker {
public:
  std::vector<Node> nodes;
  int max_depth = -1;
  void walk(const G4VPhysicalVolume* world) {
    nodes.clear();
    descend(const_cast<G4VPhysicalVolume*>(world), 0, G4Transform3D(), -1);
  }

private:
  void add(G4VPhysicalVolume* pv, const G4VSolid* solid, int depth, int copy, int kind,
           const G4Transform3D& T, int parent, bool& descend_out) {
    Node n;
    n.pv = pv; n.solid = solid; n.depth = depth; n.copy = copy; n.kind = kind; n.parent = parent;
    const G4LogicalVolume* lv = pv->GetLogicalVolume();
    const G4VisAttributes* va = lv->GetVisAttributes();
    n.has_vis = va != nullptr;
    if (va) {
      n.visible = va->IsVisible();
      n.daughters_invisible = va->IsDaughtersInvisible();
      const G4Colour& c = va->GetColour();
      n.rgba[0] = c.GetRed(); n.rgba[1] = c.GetGreen(); n.rgba[2] = c.GetBlue(); n.rgba[3] = c.GetAlpha();
      n.wire = va->IsForceDrawingStyle() && va->GetForcedDrawingStyle() == G4VisAttributes::wireframe;
    }
    if (lv->GetMaterial()) n.density = lv->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);
    n.m[0] = T.xx(); n.m[1] = T.xy(); n.m[2] = T.xz();
    n.m[3] = T.yx(); n.m[4] = T.yy(); n.m[5] = T.yz();
    n.m[6] = T.zx(); n.m[7] = T.zy(); n.m[8] = T.zz();
    n.m[9] = T.dx(); n.m[10] = T.dy(); n.m[11] = T.dz();
    nodes.push_back(n);
    descend_out = !n.daughters_invisible;
  }

  void children(G4VPhysicalVolume* pv, int depth, const G4Transform3D& T, int me) {
    if (max_depth >= 0 && depth >= max_depth) return;
    G4LogicalVolume* lv = pv->GetLogicalVolume();
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
      descend(lv->GetDaughter(i), depth + 1, T, me);
  }

  void descend(G4VPhysicalVolume* pv, int depth, const G4Transform3D& parentT, int parent) {
    bool go = true;
    EVolume vt = pv->VolumeType();
    if (vt == kNormal) {
      G4Transform3D T = parentT * G4Transform3D(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
      int me = (int)nodes.size();
      add(pv, pv->GetLogicalVolume()->GetSolid(), depth, pv->GetCopyNo(), 0, T, parent, go);
      if (go) children(pv, depth, T, me);
    } else {
      // replica / parameterised: loop over copies as G4PhysicalVolumeModel does
      G4int nrep = 0; EAxis axis; G4double width, offset; G4bool consuming;
      pv->GetReplicationData(axis, nrep, width, offset, consuming);
      G4VPVParameterisation* param = pv->GetParameterisation();
      G4LogicalVolume* lv = pv->GetLogicalVolume();
      G4ReplicaNavigation repNav;
      for (G4int n = 0; n < nrep; ++n) {
        const G4VSolid* solid = lv->GetSolid();
        if (param) {
          param->ComputeTransformation(n, pv);
          G4VSolid* s = param->ComputeSolid(n, pv);
          if (s) { s->ComputeDimensions(param, n, pv); solid = s; }
        } else {
          repNav.ComputeTransformation(n, pv);
        }
        G4Transform3D T = parentT * G4Transform3D(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
        int me = (int)nodes.size();
        add(pv, solid, depth, n, param ? 2 : 1, T, parent, go);
        if (go) children(pv, depth, T, me);
      }
    }
  }
};

// Geometry without a run manager: the same conversion DDG4's Geant4DetectorGeometryConstruction does.
inline G4VPhysicalVolume* convert(dd4hep::Detector& description, int level = 3) {
  dd4hep::sim::Geant4Mapping& g4map = dd4hep::sim::Geant4Mapping::instance();
  dd4hep::sim::Geant4Converter conv(description, dd4hep::PrintLevel(level));
  dd4hep::sim::Geant4GeometryInfo* info = conv.create(description.world()).detach();
  g4map.attach(info);
  g4map.volumeManager();
  return info->world();
}

inline std::vector<double> extent(const G4VPhysicalVolume* pv) {
  G4VisExtent e = pv->GetLogicalVolume()->GetSolid()->GetExtent();
  G4Point3D c = e.GetExtentCentre();
  return {c.x(), c.y(), c.z(), e.GetExtentRadius(), e.GetXmin(), e.GetXmax(), e.GetYmin(), e.GetYmax(), e.GetZmin(), e.GetZmax()};
}

inline dd4hep::sim::Geant4GeometryInfo& info() { return dd4hep::sim::Geant4Mapping::instance().data(); }

// volume ID of a placement path (world excluded), DDG4 convention: deepest volume first.
inline long long volume_id(const std::vector<const G4VPhysicalVolume*>& path) {
  auto vm = dd4hep::sim::Geant4Mapping::instance().volumeManager();
  std::pair<dd4hep::VolumeID, std::vector<std::pair<const dd4hep::BitFieldElement*, dd4hep::VolumeID>>> d;
  vm.volumeDescriptor(path, d);
  return (long long)d.first;
}

// DDG4 volume ID the way a sensitive detector gets it during tracking: navigator touchable at a point.
// Returns the ID and (through pv_out) the physical volume found there.
inline long long volume_id_at(const G4VPhysicalVolume* world, double x, double y, double z,
                              const G4VPhysicalVolume** pv_out = nullptr) {
  static G4Navigator* nav = nullptr;
  if (!nav) { nav = new G4Navigator(); nav->SetWorldVolume(const_cast<G4VPhysicalVolume*>(world)); }
  nav->LocateGlobalPointAndSetup(G4ThreeVector(x, y, z), nullptr, false, true);
  G4TouchableHistory* t = nav->CreateTouchableHistory();
  if (pv_out) *pv_out = t->GetVolume(0);
  auto vm = dd4hep::sim::Geant4Mapping::instance().volumeManager();
  long long id = (long long)vm.volumeID(t);
  delete t;
  return id;
}

// G4 PV -> TGeo node for ordinary placements; imprints (assemblies) -> chain of TGeo nodes
struct PVMap {
  std::map<const G4VPhysicalVolume*, const TGeoNode*> place;
  std::map<const G4VPhysicalVolume*, std::vector<const TGeoNode*>> imprint;
  void build() {
    auto& d = info();
    for (auto& e : d.g4Placements) place[e.second] = e.first.ptr();
    for (auto& e : d.g4VolumeImprints)
      for (auto& im : e.second) imprint[im.second] = im.first;
  }
};

}  // namespace solidvis
