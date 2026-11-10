// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// Scene producer: see SceneDump.h.

#include "SceneDump.h"

#include "NpzWriter.h"
#include "SHA256.h"

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/IDDescriptor.h>
#include <DD4hep/Objects.h>
#include <DD4hep/Printout.h>
#include <DD4hep/Readout.h>
#include <DD4hep/Volumes.h>
#include <DDG4/Geant4Converter.h>
#include <DDG4/Geant4GeometryInfo.h>
#include <DDG4/Geant4Mapping.h>

#include <G4LogicalVolume.hh>
#include <G4Polyhedron.hh>
#include <G4ReflectionFactory.hh>
#include <G4ReplicaNavigation.hh>
#include <G4Transform3D.hh>
#include <G4VPVParameterisation.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4Version.hh>
#include <G4VisAttributes.hh>
#include <G4VisExtent.hh>

#include <RVersion.h>
#include <TGeoNode.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#ifndef SOLID_VERSION
#define SOLID_VERSION "unknown"
#endif

namespace solid::scene {

namespace {

using NodeChain = std::vector<const TGeoNode*>;

// ---------------------------------------------------------------------------------------------------------
// Geant4 polyhedron of a solid, as the vis drivers request it
bool tessellate(const G4VSolid* solid, int segments, Mesh& mesh)
{
  G4Polyhedron::SetNumberOfRotationSteps(segments);
  const G4Polyhedron* p = solid->GetPolyhedron(); // owned and cached by the solid
  G4Polyhedron::ResetNumberOfRotationSteps();
  if (p == nullptr || p->GetNoFacets() <= 0) {
    return false;
  }
  const int nv = p->GetNoVertices();
  const int nf = p->GetNoFacets();
  mesh.v.reserve(3 * static_cast<std::size_t>(nv));
  mesh.f.reserve(4 * static_cast<std::size_t>(nf));
  for (int i = 1; i <= nv; ++i) {
    const G4Point3D q = p->GetVertex(i);
    mesh.v.insert(mesh.v.end(), {q.x(), q.y(), q.z()});
  }
  for (int i = 1; i <= nf; ++i) {
    G4int n = 0, nodes[4] = {0, 0, 0, 0}, flags[4] = {1, 1, 1, 1};
    p->GetFacet(i, n, nodes, flags);
    for (int k = 0; k < 4; ++k) {
      mesh.f.push_back(k < n ? nodes[k] * (flags[k] > 0 ? 1 : -1) : 0);
    }
  }
  mesh.name = solid->GetName();
  // DD4hep puts the address of the TGeo shape into the name ("_0x55f7846b4bd0"): drop it, it changes from run to run
  for (auto p0 = mesh.name.find("_0x"); p0 != std::string::npos; p0 = mesh.name.find("_0x", p0)) {
    auto p1 = mesh.name.find_first_not_of("0123456789abcdefABCDEF", p0 + 3);
    if (p1 == p0 + 3) {
      p0 += 3; // not an address
      continue;
    }
    mesh.name.erase(p0, p1 == std::string::npos ? std::string::npos : p1 - p0);
  }
  mesh.type = solid->GetEntityType();
  return true;
}

// ---------------------------------------------------------------------------------------------------------
// DD4hep side of the identification: which TGeo placements a Geant4 physical volume stands for
struct Identifier {
  dd4hep::Detector&                                                        description;
  std::unordered_map<const G4VPhysicalVolume*, std::vector<const TGeoNode*>> placed; // ordinary placements
  std::unordered_map<const G4VPhysicalVolume*, NodeChain>                    imprinted; // volumes inside assemblies
  struct Element {
    std::string path;
    std::string owner; // name of the subdetector whose readout applies
  };
  std::unordered_map<const TGeoNode*, Element> elements;
  struct ReadoutInfo {
    std::string                                           name;
    std::vector<std::pair<std::string, const dd4hep::BitFieldElement*>> fields;
    std::string                                           description;
  };
  std::map<std::string, ReadoutInfo> readouts; // by subdetector name

  explicit Identifier(dd4hep::Detector& desc) : description(desc)
  {
    const auto& info = dd4hep::sim::Geant4Mapping::instance().data();
    for (const auto& [pv, g4pv] : info.g4Placements) {
      if (g4pv != nullptr) {
        placed[g4pv].push_back(pv.ptr());
      }
    }
    for (const auto& [volume, imprints] : info.g4VolumeImprints) {
      for (const auto& [chain, g4pv] : imprints) {
        imprinted[g4pv] = chain;
      }
    }
    for (const auto& [name, child] : desc.world().children()) {
      index(child, name, "");
    }
  }

  void index(dd4hep::DetElement de, const std::string& top, std::string owner)
  {
    const auto& sds = description.sensitiveDetectors();
    if (sds.find(de.name()) != sds.end()) {
      owner = de.name(); // e.g. a tracker inside an assembly DetElement
    }
    const dd4hep::PlacedVolume pv = de.placement();
    if (pv.isValid()) {
      elements[pv.ptr()] = {de.path(), owner.empty() ? top : owner};
    }
    for (const auto& [name, child] : de.children()) {
      index(child, top, owner);
    }
  }

  // TGeo nodes of Geant4 volume `pv`, daughter number `index` of `mother`, relative to the mother's node
  NodeChain chainOf(const G4VPhysicalVolume* pv, const G4LogicalVolume* mother, std::size_t index, int copy, int kind) const
  {
    if (auto it = placed.find(pv); it != placed.end() && !it->second.empty()) {
      const TGeoNode* node = it->second.front();
      if (kind != 0) {
        // replicas and parameterised volumes: one Geant4 volume for the n TGeo placements, in copy order
        const dd4hep::PlacedVolume first(node);
        const auto*                data = first.data();
        if (data != nullptr && data->params != nullptr && copy >= 0 &&
            static_cast<std::size_t>(copy) < data->params->placements.size()) {
          node = data->params->placements[static_cast<std::size_t>(copy)].ptr();
        }
      }
      return {node};
    }
    if (auto it = imprinted.find(pv); it != imprinted.end()) {
      return it->second;
    }
    // Daughters of a reflected logical volume are made by G4ReflectionFactory, not by DDG4: use the same
    // daughter of the volume it was reflected from.
    auto* factory = G4ReflectionFactory::Instance();
    auto* lv      = const_cast<G4LogicalVolume*>(mother);
    if (lv != nullptr && factory->IsReflected(lv)) {
      const G4LogicalVolume* original = factory->GetConstituentLV(lv);
      if (original != nullptr && original != mother) {
        const std::size_t n = original->GetNoDaughters();
        if (index < n && original->GetDaughter(index)->GetName() == pv->GetName()) {
          return chainOf(original->GetDaughter(index), original, index, copy, kind);
        }
        for (std::size_t i = 0; i < n; ++i) {
          const G4VPhysicalVolume* d = original->GetDaughter(i);
          if (d->GetName() == pv->GetName() && d->GetCopyNo() == pv->GetCopyNo()) {
            return chainOf(d, original, i, copy, kind);
          }
        }
      }
    }
    return {};
  }

  const ReadoutInfo& readoutOf(const std::string& subdetector)
  {
    auto it = readouts.find(subdetector);
    if (it != readouts.end()) {
      return it->second;
    }
    ReadoutInfo  r;
    const auto&  sds = description.sensitiveDetectors();
    if (auto is = sds.find(subdetector); is != sds.end()) {
      const dd4hep::SensitiveDetector sd(is->second);
      const dd4hep::Readout           ro = sd.isValid() ? sd.readout() : dd4hep::Readout();
      if (ro.isValid() && ro.idSpec().isValid()) {
        r.name        = ro.name();
        r.fields      = ro.idSpec().fields();
        r.description = ro.idSpec().fieldDescription();
      }
    }
    return readouts.emplace(subdetector, std::move(r)).first->second;
  }

  void identify(const NodeChain& chain, Node& node)
  {
    if (chain.empty()) {
      return;
    }
    const dd4hep::Volume volume(chain.back()->GetVolume());
    if (volume.isValid()) {
      const dd4hep::VisAttr vis = volume.visAttributes();
      if (vis.isValid()) {
        node.drawing_style = vis.drawingStyle() == dd4hep::VisAttr::WIREFRAME ? 2 : 1;
      }
    }
    const Element* element = nullptr;
    for (auto it = chain.rbegin(); it != chain.rend() && element == nullptr; ++it) {
      if (auto ie = elements.find(*it); ie != elements.end()) {
        element = &ie->second;
      }
    }
    if (element == nullptr) {
      return;
    }
    node.det_path      = element->path;
    node.sensitive     = volume.isValid() && volume.isSensitive();
    const auto& ro     = readoutOf(element->owner);
    if (ro.name.empty()) {
      return;
    }
    node.readout = ro.name;
    // volume IDs set by the placements from the top down; a deeper placement overrides a field set above it
    std::map<std::string, long long> ids;
    for (const TGeoNode* n : chain) {
      const dd4hep::PlacedVolume pv(n);
      if (pv.data() == nullptr) {
        continue;
      }
      for (const auto& id : pv.volIDs()) {
        ids[id.first] = id.second;
      }
    }
    for (const auto& [name, field] : ro.fields) {
      if (auto ii = ids.find(name); ii != ids.end() && field != nullptr) {
        const std::uint64_t mask = static_cast<std::uint64_t>(field->mask());
        node.volume_id |= (static_cast<std::uint64_t>(ii->second) << field->offset()) & mask;
        node.volume_mask |= mask;
      }
    }
  }
};

// ---------------------------------------------------------------------------------------------------------
class Walker {
public:
  Walker(Scene& scene, const Options& options, Identifier* identifier)
    : m_scene(scene), m_options(options), m_identifier(identifier)
  {
    for (double& x : m_lo) {
      x = std::numeric_limits<double>::max();
    }
    for (double& x : m_hi) {
      x = std::numeric_limits<double>::lowest();
    }
  }

  void walk(const G4VPhysicalVolume* world)
  {
    descend(const_cast<G4VPhysicalVolume*>(world), nullptr, 0, 0, G4Transform3D(), -1, {});
    if (m_extent) {
      double r2 = 0;
      for (int k = 0; k < 3; ++k) {
        m_scene.bounding_sphere[k] = 0.5 * (m_lo[k] + m_hi[k]);
        r2 += 0.25 * (m_hi[k] - m_lo[k]) * (m_hi[k] - m_lo[k]);
      }
      m_scene.bounding_sphere[3] = std::sqrt(r2);
    }
  }

private:
  Scene&         m_scene;
  const Options& m_options;
  Identifier*    m_identifier;
  bool           m_extent = false;
  double         m_lo[3], m_hi[3];
  std::unordered_map<const G4VSolid*, std::vector<std::int32_t>> m_meshes; // meshes made from each solid

  std::int32_t meshOf(const G4VSolid* solid, bool may_change)
  {
    auto& known = m_meshes[solid];
    if (!known.empty() && !may_change) {
      return known.front();
    }
    Mesh mesh;
    if (!tessellate(solid, m_options.segments, mesh)) {
      if (known.empty()) {
        dd4hep::printout(dd4hep::WARNING, "SoLIDScene", "no polyhedron for solid %s (%s)", solid->GetName().c_str(),
                         solid->GetEntityType().c_str());
      }
      known.push_back(-1);
      ++m_scene.failed_meshes;
      return -1;
    }
    // a parameterisation may resize the same solid object for each copy: reuse a mesh only if it is equal
    for (std::int32_t k : known) {
      if (k >= 0 && m_scene.meshes[static_cast<std::size_t>(k)].v == mesh.v &&
          m_scene.meshes[static_cast<std::size_t>(k)].f == mesh.f) {
        return k;
      }
    }
    m_scene.meshes.push_back(std::move(mesh));
    known.push_back(static_cast<std::int32_t>(m_scene.meshes.size() - 1));
    return known.back();
  }

  // bounding box of the drawn volumes, accumulated like G4BoundingExtentScene does (extent box of each solid,
  // transformed, then the box around all of them)
  void accrue(const G4VSolid* solid, const G4Transform3D& transform)
  {
    G4VisExtent e = solid->GetExtent();
    e.Transform(transform);
    const double lo[3] = {e.GetXmin(), e.GetYmin(), e.GetZmin()};
    const double hi[3] = {e.GetXmax(), e.GetYmax(), e.GetZmax()};
    for (int k = 0; k < 3; ++k) {
      m_lo[k] = std::min(m_lo[k], lo[k]);
      m_hi[k] = std::max(m_hi[k], hi[k]);
    }
    m_extent = true;
  }

  // returns true if the daughters are to be visited
  bool add(const G4VPhysicalVolume* pv, const G4VSolid* solid, int depth, int copy, int kind, const G4Transform3D& T,
           int parent, const NodeChain& chain)
  {
    Node n;
    n.name   = pv->GetName();
    n.copy   = copy;
    n.depth  = depth;
    n.parent = parent;
    n.kind   = static_cast<std::int8_t>(kind);
    n.matrix = {T.xx(), T.xy(), T.xz(), T.yx(), T.yy(), T.yz(), T.zx(), T.zy(), T.zz(), T.dx(), T.dy(), T.dz()};
    const G4VisAttributes* va = pv->GetLogicalVolume()->GetVisAttributes();
    if (va != nullptr) {
      const G4Colour& c     = va->GetColour();
      n.rgba                = {c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha()};
      n.visible             = va->IsVisible();
      n.daughters_invisible = va->IsDaughtersInvisible();
      n.wireframe           = va->IsForceDrawingStyle() && va->GetForcedDrawingStyle() == G4VisAttributes::wireframe;
    }
    if (m_identifier != nullptr) {
      m_identifier->identify(chain, n);
    }
    if (n.visible && n.rgba[3] > 0) { // G4PhysicalVolumeModel: culling of invisible volumes
      n.mesh = meshOf(solid, kind == 2);
      if (n.mesh >= 0) {
        accrue(solid, T);
      }
    }
    const bool go = !n.daughters_invisible;
    m_scene.nodes.push_back(std::move(n));
    return go;
  }

  void daughters(const G4VPhysicalVolume* pv, int depth, const G4Transform3D& T, int me, const NodeChain& chain)
  {
    const G4LogicalVolume* lv = pv->GetLogicalVolume();
    const std::size_t      n  = lv->GetNoDaughters();
    for (std::size_t i = 0; i < n; ++i) {
      descend(lv->GetDaughter(i), lv, i, depth + 1, T, me, chain);
    }
  }

  NodeChain extend(const NodeChain& above, const G4VPhysicalVolume* pv, const G4LogicalVolume* mother, std::size_t index,
                   int copy, int kind) const
  {
    NodeChain chain = above;
    if (m_identifier != nullptr && mother != nullptr) {
      const NodeChain own = m_identifier->chainOf(pv, mother, index, copy, kind);
      chain.insert(chain.end(), own.begin(), own.end());
    }
    return chain;
  }

  void descend(G4VPhysicalVolume* pv, const G4LogicalVolume* mother, std::size_t index, int depth,
               const G4Transform3D& parentT, int parent, const NodeChain& above)
  {
    if (pv->VolumeType() == kNormal) {
      const G4Transform3D T     = parentT * G4Transform3D(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
      const NodeChain     chain = extend(above, pv, mother, index, pv->GetCopyNo(), 0);
      const int           me    = static_cast<int>(m_scene.nodes.size());
      if (add(pv, pv->GetLogicalVolume()->GetSolid(), depth, pv->GetCopyNo(), 0, T, parent, chain)) {
        daughters(pv, depth, T, me, chain);
      }
      return;
    }
    // replica or parameterised volume: one node per copy, as G4PhysicalVolumeModel draws them
    EAxis    axis = kUndefined;
    G4int    count = 0;
    G4double width = 0, offset = 0;
    G4bool   consuming = false;
    pv->GetReplicationData(axis, count, width, offset, consuming);
    G4VPVParameterisation* param = pv->GetParameterisation();
    G4LogicalVolume*       lv    = pv->GetLogicalVolume();
    G4ReplicaNavigation    navigation;
    for (G4int copy = 0; copy < count; ++copy) {
      G4VSolid* solid = lv->GetSolid();
      if (param != nullptr) {
        solid = param->ComputeSolid(copy, pv);
        solid->ComputeDimensions(param, copy, pv);
        param->ComputeTransformation(copy, pv);
      }
      else {
        navigation.ComputeTransformation(copy, pv);
      }
      const int           kind  = param != nullptr ? 2 : 1;
      const G4Transform3D T     = parentT * G4Transform3D(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
      const NodeChain     chain = extend(above, pv, mother, index, copy, kind);
      const int           me    = static_cast<int>(m_scene.nodes.size());
      if (add(pv, solid, depth, copy, kind, T, parent, chain)) {
        daughters(pv, depth, T, me, chain);
      }
    }
  }
};

// ---------------------------------------------------------------------------------------------------------
// Hash of the scene content with lengths rounded to 0.1 um and matrix elements and colours to 1e-9, so that
// the last bits of floating-point results do not matter.  Software versions and file names are not hashed.
class Hasher {
public:
  void put(std::int64_t v)
  {
    unsigned char b[8];
    for (int i = 0; i < 8; ++i) {
      b[i] = static_cast<unsigned char>((static_cast<std::uint64_t>(v) >> (8 * i)) & 0xff);
    }
    m_sha.update(b, 8);
  }
  void put(const std::string& s)
  {
    put(static_cast<std::int64_t>(s.size()));
    m_sha.update(reinterpret_cast<const unsigned char*>(s.data()), s.size());
  }
  void length(double mm) { put(static_cast<std::int64_t>(std::llround(mm * 1e4))); }
  void number(double x) { put(static_cast<std::int64_t>(std::llround(x * 1e9))); }
  std::string hex() { return m_sha.hexdigest(); }

private:
  SHA256 m_sha;
};

std::string hashOf(const Scene& s)
{
  Hasher h;
  h.put(std::string("solid-scene-v") + std::to_string(SCHEMA_VERSION));
  h.put(static_cast<std::int64_t>(s.segments));
  h.put(static_cast<std::int64_t>(s.meshes.size()));
  for (const auto& m : s.meshes) {
    h.put(m.type);
    h.put(static_cast<std::int64_t>(m.v.size()));
    for (double x : m.v) {
      h.length(x);
    }
    h.put(static_cast<std::int64_t>(m.f.size()));
    for (std::int32_t i : m.f) {
      h.put(static_cast<std::int64_t>(i));
    }
  }
  h.put(static_cast<std::int64_t>(s.nodes.size()));
  for (const auto& n : s.nodes) {
    h.put(n.name);
    for (std::int64_t i : {static_cast<std::int64_t>(n.copy), static_cast<std::int64_t>(n.depth),
                           static_cast<std::int64_t>(n.parent), static_cast<std::int64_t>(n.mesh),
                           static_cast<std::int64_t>(n.kind), static_cast<std::int64_t>(n.visible),
                           static_cast<std::int64_t>(n.daughters_invisible), static_cast<std::int64_t>(n.wireframe),
                           static_cast<std::int64_t>(n.drawing_style), static_cast<std::int64_t>(n.sensitive),
                           static_cast<std::int64_t>(n.volume_id), static_cast<std::int64_t>(n.volume_mask)}) {
      h.put(i);
    }
    for (int k = 0; k < 9; ++k) {
      h.number(n.matrix[static_cast<std::size_t>(k)]);
    }
    for (int k = 9; k < 12; ++k) {
      h.length(n.matrix[static_cast<std::size_t>(k)]);
    }
    for (double c : n.rgba) {
      h.number(c);
    }
    h.put(n.det_path);
    h.put(n.readout);
  }
  return h.hex();
}

std::string jsonString(const std::string& s)
{
  std::string out = "\"";
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') {
      out += '\\';
      out += static_cast<char>(c);
    }
    else if (c < 0x20) {
      char buf[8];
      std::snprintf(buf, sizeof(buf), "\\u%04x", c);
      out += buf;
    }
    else {
      out += static_cast<char>(c);
    }
  }
  return out + "\"";
}

std::string number(double x)
{
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", x);
  return buf;
}

std::string geant4Version()
{
  // "$Name: geant4-11-04-patch-01 [MT]$" -> "geant4-11-04-patch-01 [MT]"
  std::string v = G4Version;
  if (const auto a = v.find(": "); a != std::string::npos) {
    v = v.substr(a + 2);
  }
  while (!v.empty() && (v.back() == '$' || v.back() == ' ')) {
    v.pop_back();
  }
  return v;
}

std::string metaJson(const Scene& s)
{
  std::size_t drawn = 0;
  for (const auto& n : s.nodes) {
    drawn += n.mesh >= 0 ? 1 : 0;
  }
  std::ostringstream o;
  o << "{\"schema\": " << SCHEMA_VERSION << ", \"config\": " << jsonString(s.config) << ", \"mode\": "
    << jsonString(s.mode) << ", \"source\": \"Geant4 G4VSolid::GetPolyhedron\", \"units\": \"mm\", \"segments\": "
    << s.segments << ", \"nodes\": " << s.nodes.size() << ", \"drawn\": " << drawn << ", \"meshes\": "
    << s.meshes.size() << ", \"failed_meshes\": " << s.failed_meshes << ", \"bounding_sphere\": ["
    << number(s.bounding_sphere[0]) << ", " << number(s.bounding_sphere[1]) << ", " << number(s.bounding_sphere[2])
    << ", " << number(s.bounding_sphere[3]) << "], \"geometry_hash\": " << jsonString(s.geometry_hash)
    << ", \"versions\": {\"solid_dd4hep\": " << jsonString(SOLID_VERSION) << ", \"dd4hep\": "
    << jsonString(dd4hep::versionString()) << ", \"geant4\": " << jsonString(geant4Version()) << ", \"root\": "
    << jsonString(ROOT_RELEASE) << "}, \"readouts\": {";
  bool first = true;
  for (const auto& [name, fields] : s.readouts) {
    o << (first ? "" : ", ") << jsonString(name) << ": " << jsonString(fields);
    first = false;
  }
  o << "}}";
  return o.str();
}

} // namespace

// ---------------------------------------------------------------------------------------------------------
G4VPhysicalVolume* convertGeometry(dd4hep::Detector& description)
{
  auto& mapping = dd4hep::sim::Geant4Mapping::instance();
  dd4hep::sim::Geant4Converter converter(description, dd4hep::VERBOSE); // level of its own messages
  converter.checkOverlaps = false;
  dd4hep::sim::Geant4GeometryInfo* info = converter.create(description.world()).detach();
  mapping.attach(info);
  return info->world();
}

Scene collect(dd4hep::Detector& description, const G4VPhysicalVolume* world, const Options& options)
{
  if (world == nullptr) {
    throw std::runtime_error("SoLIDScene: no Geant4 world volume");
  }
  Scene scene;
  scene.config   = options.config;
  scene.mode     = options.mode;
  scene.segments = options.segments;
  std::unique_ptr<Identifier> identifier;
  if (options.ids) {
    identifier = std::make_unique<Identifier>(description);
  }
  Walker walker(scene, options, identifier.get());
  walker.walk(world);
  if (identifier) {
    for (const auto& [subdetector, ro] : identifier->readouts) {
      if (!ro.name.empty()) {
        scene.readouts[ro.name] = ro.description;
      }
    }
  }
  scene.geometry_hash = hashOf(scene);
  return scene;
}

void writeNpz(const Scene& s, const std::string& path)
{
  NpzWriter npz;
  // meshes
  std::vector<double>       mv;
  std::vector<std::int32_t> mf;
  std::vector<std::int64_t> voff{0}, foff{0};
  std::vector<std::string>  mname, mtype;
  for (const auto& m : s.meshes) {
    mv.insert(mv.end(), m.v.begin(), m.v.end());
    mf.insert(mf.end(), m.f.begin(), m.f.end());
    voff.push_back(static_cast<std::int64_t>(mv.size() / 3));
    foff.push_back(static_cast<std::int64_t>(mf.size() / 4));
    mname.push_back(m.name);
    mtype.push_back(m.type);
  }
  npz.addScalar("meta", metaJson(s));
  npz.add("mesh_vertices", mv, {mv.size() / 3, 3});
  npz.add("mesh_vertex_offset", voff);
  npz.add("mesh_facets", mf, {mf.size() / 4, 4});
  npz.add("mesh_facet_offset", foff);
  npz.add("mesh_name", mname);
  npz.add("mesh_type", mtype);
  // placements
  const std::size_t          n = s.nodes.size();
  std::vector<std::string>   name, det, readout;
  std::vector<std::int32_t>  copy, depth, parent, mesh;
  std::vector<std::int8_t>   kind, style;
  std::vector<double>        matrix, rgba;
  std::vector<std::uint8_t>  visible, wire, dinv, sensitive;
  std::vector<std::uint64_t> volid, mask;
  for (const auto& nd : s.nodes) {
    name.push_back(nd.name);
    copy.push_back(nd.copy);
    depth.push_back(nd.depth);
    parent.push_back(nd.parent);
    mesh.push_back(nd.mesh);
    kind.push_back(nd.kind);
    matrix.insert(matrix.end(), nd.matrix.begin(), nd.matrix.end());
    rgba.insert(rgba.end(), nd.rgba.begin(), nd.rgba.end());
    visible.push_back(nd.visible ? 1 : 0);
    wire.push_back(nd.wireframe ? 1 : 0);
    dinv.push_back(nd.daughters_invisible ? 1 : 0);
    style.push_back(nd.drawing_style);
    det.push_back(nd.det_path);
    readout.push_back(nd.readout);
    volid.push_back(nd.volume_id);
    mask.push_back(nd.volume_mask);
    sensitive.push_back(nd.sensitive ? 1 : 0);
  }
  npz.add("node_name", name);
  npz.add("node_copy", copy);
  npz.add("node_depth", depth);
  npz.add("node_parent", parent);
  npz.add("node_mesh", mesh);
  npz.add("node_kind", kind);
  npz.add("node_matrix", matrix, {n, 12});
  npz.add("node_rgba", rgba, {n, 4});
  npz.addBool("node_visible", visible);
  npz.addBool("node_wireframe", wire);
  npz.addBool("node_daughters_invisible", dinv);
  npz.add("node_drawing_style", style);
  npz.add("node_det_path", det);
  npz.add("node_readout", readout);
  npz.add("node_volume_id", volid);
  npz.add("node_volume_mask", mask);
  npz.addBool("node_sensitive", sensitive);
  npz.write(path);
}

void writePrim(const Scene& s, const std::string& path)
{
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (f == nullptr) {
    throw std::runtime_error("SoLIDScene: cannot open " + path + " for writing");
  }
  const auto& b = s.bounding_sphere;
  std::fprintf(f, "##G4.PRIM-FORMAT-2.4\n");
  std::fprintf(f, "/BoundingBox %.9g %.9g %.9g %.9g %.9g %.9g\n", b[0] - b[3], b[1] - b[3], b[2] - b[3], b[0] + b[3],
               b[1] + b[3], b[2] + b[3]);
  std::fprintf(f, "!SetCamera\n!OpenDevice\n!BeginModeling\n");
  std::fprintf(f, "/Ndiv %d\n", s.segments);
  for (const auto& n : s.nodes) {
    if (n.mesh < 0) {
      continue;
    }
    const Mesh& m = s.meshes[static_cast<std::size_t>(n.mesh)];
    const auto& M = n.matrix;
    std::fprintf(f, "#--------------------\n#/PVName%*s%s:%d\n", 5 + 5 * n.depth, "", n.name.c_str(), n.copy);
    std::fprintf(f, "/ColorRGB %.9g %.9g %.9g\n", n.rgba[0], n.rgba[1], n.rgba[2]);
    std::fprintf(f, "/ForceWireframe  %d\n", n.wireframe ? 1 : 0);
    const double det = M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) +
                       M[2] * (M[3] * M[7] - M[4] * M[6]);
    const std::size_t nv = m.v.size() / 3, nf = m.f.size() / 4;
    if (det > 0) {
      // /BaseVector: images of the local x and y axes (first two columns of R); DAWN takes z = x cross y
      std::fprintf(f, "/Origin %.9g %.9g %.9g\n", M[9], M[10], M[11]);
      std::fprintf(f, "/BaseVector %.9g %.9g %.9g %.9g %.9g %.9g\n", M[0], M[3], M[6], M[1], M[4], M[7]);
      std::fprintf(f, "/Polyhedron\n");
      for (std::size_t i = 0; i < nv; ++i) {
        std::fprintf(f, "/Vertex %.9g %.9g %.9g\n", m.v[3 * i], m.v[3 * i + 1], m.v[3 * i + 2]);
      }
      for (std::size_t i = 0; i < nf; ++i) {
        std::fprintf(f, "/Facet");
        for (int k = 0; k < 4; ++k) {
          if (m.f[4 * i + k] != 0) {
            std::fprintf(f, " %d", std::abs(m.f[4 * i + k]));
          }
        }
        std::fprintf(f, "\n");
      }
    }
    else {
      // a transform with a reflection cannot be written as /BaseVector: world coordinates, reversed winding
      std::fprintf(f, "/Origin 0 0 0\n/BaseVector 1 0 0 0 1 0\n/Polyhedron\n");
      for (std::size_t i = 0; i < nv; ++i) {
        const double x = m.v[3 * i], y = m.v[3 * i + 1], z = m.v[3 * i + 2];
        std::fprintf(f, "/Vertex %.9g %.9g %.9g\n", M[0] * x + M[1] * y + M[2] * z + M[9],
                     M[3] * x + M[4] * y + M[5] * z + M[10], M[6] * x + M[7] * y + M[8] * z + M[11]);
      }
      for (std::size_t i = 0; i < nf; ++i) {
        std::fprintf(f, "/Facet");
        for (int k = 3; k >= 0; --k) {
          if (m.f[4 * i + k] != 0) {
            std::fprintf(f, " %d", std::abs(m.f[4 * i + k]));
          }
        }
        std::fprintf(f, "\n");
      }
    }
    std::fprintf(f, "/EndPolyhedron\n");
  }
  std::fprintf(f, "!EndModeling\n!DrawAll\n!CloseDevice\n");
  if (std::fclose(f) != 0) {
    throw std::runtime_error("SoLIDScene: error writing " + path);
  }
}

Scene dump(dd4hep::Detector& description, const G4VPhysicalVolume* world, const std::string& npz,
           const std::string& prim, const Options& options)
{
  if (world == nullptr) {
    world = convertGeometry(description);
  }
  const auto t0    = std::chrono::steady_clock::now();
  Scene      scene = collect(description, world, options);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::size_t drawn = 0;
  for (const auto& n : scene.nodes) {
    drawn += n.mesh >= 0 ? 1 : 0;
  }
  if (!npz.empty()) {
    writeNpz(scene, npz);
  }
  if (!prim.empty()) {
    writePrim(scene, prim);
  }
  dd4hep::printout(dd4hep::INFO, "SoLIDScene", "%s: %zu volumes, %zu drawn, %zu meshes in %.2f s, hash %.16s%s%s%s%s",
                   scene.config.c_str(), scene.nodes.size(), drawn, scene.meshes.size(), seconds,
                   scene.geometry_hash.c_str(), npz.empty() ? "" : " -> ", npz.c_str(), prim.empty() ? "" : " -> ", prim.c_str());
  return scene;
}

} // namespace solid::scene
