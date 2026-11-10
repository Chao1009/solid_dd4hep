// Prototype: tessellate the DD4hep TGeo geometry with ROOT only (no Geant4).
//   TGeoShape::MakeBuffer3D() (raw points/segments/polygons, gGeoManager->SetNsegments(n))
//   -> RootCsg::ConvertToMesh(); composite shapes -> RootCsg::BuildUnion/Difference/Intersection on the
//   transformed operand meshes (the same recursion ROOT's web geometry viewer RGeomDescription uses).
#pragma once
#include <vector>
#include <string>
#include <memory>
#include <chrono>

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"
#include "TGeoMatrix.h"
#include "TGeoScaledShape.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoPcon.h"
#include "TGeoSphere.h"
#include "TGeoTorus.h"
#include "TBuffer3D.h"
#include "CsgOps.h"

#include "DD4hep/Volumes.h"
#include "DD4hep/Objects.h"

namespace solidtgeo {

struct Mesh {
  std::vector<double> v;     // 3*nv local frame [cm, TGeo units]
  std::vector<int> pol;      // n, i0 .. i(n-1), n, ...  (0-based)
  int nv = 0, npol = 0;
  bool ok = false;
  double seconds = 0;
  std::string err;
};

// options of the prototype (set from Python)
struct Options {
  bool scale_phi = true;        // segments proportional to the phi extent (like G4 / lineSegmentsPerCircle)
  bool full_seg_as_tube = true; // TGeoTubeSeg/ConeSeg with dphi = 360: tessellate as TGeoTube/Cone (no seam faces)
  bool sphere_divisions = true; // TGeoSphere keeps its own divisions (20): set them from nseg
  bool honour_drawing_style = false; // VisAttr drawingStyle="wireframe" -> /ForceWireframe 1
};
inline Options& options() { static Options o; return o; }

// phi extent in degrees of the shapes whose TBuffer3D uses gGeoManager->GetNsegments() per *own* phi range
inline double phi_extent(const TGeoShape* s) {
  if (auto t = dynamic_cast<const TGeoTubeSeg*>(s)) return t->GetPhi2() - t->GetPhi1();
  if (auto c = dynamic_cast<const TGeoConeSeg*>(s)) return c->GetPhi2() - c->GetPhi1();
  if (auto p = dynamic_cast<const TGeoPcon*>(s)) if (!dynamic_cast<const TGeoPgon*>(s)) return p->GetDphi();
  if (auto t = dynamic_cast<const TGeoTorus*>(s)) return t->GetDphi();
  return 360.;
}

inline std::unique_ptr<RootCsg::TBaseMesh> make_mesh(const TGeoMatrix* matr, const TGeoShape* shape, std::string& err,
                                                     int nseg = 72) {
  std::unique_ptr<RootCsg::TBaseMesh> res;
  auto comp = dynamic_cast<const TGeoCompositeShape*>(shape);
  if (!comp) {
    const Options& o = options();
    std::unique_ptr<TGeoShape> tmp;
    double dphi = phi_extent(shape);
    if (o.full_seg_as_tube && std::abs(dphi - 360.) < 1e-9) {
      if (auto t = dynamic_cast<const TGeoTubeSeg*>(shape))
        { tmp.reset(new TGeoTube(t->GetRmin(), t->GetRmax(), t->GetDz())); shape = tmp.get(); }
      else if (auto c = dynamic_cast<const TGeoConeSeg*>(shape))
        { tmp.reset(new TGeoCone(c->GetDz(), c->GetRmin1(), c->GetRmax1(), c->GetRmin2(), c->GetRmax2())); shape = tmp.get(); }
    }
    if (o.sphere_divisions)
      if (auto sp = dynamic_cast<const TGeoSphere*>(shape)) {
        double dp = sp->GetPhi2() - sp->GetPhi1();
        int n = o.scale_phi ? std::max(3, (int)std::lround(nseg * dp / 360.)) : nseg;
        const_cast<TGeoSphere*>(sp)->SetNumberOfDivisions(n);
      }
    int old = gGeoManager->GetNsegments();
    if (o.scale_phi && dphi < 360. - 1e-9) gGeoManager->SetNsegments(std::max(3, (int)std::lround(nseg * dphi / 360.)));
    std::unique_ptr<TBuffer3D> b3d(shape->MakeBuffer3D());
    gGeoManager->SetNsegments(old);
    if (!b3d) { err += std::string("no MakeBuffer3D for ") + shape->ClassName() + "; "; return res; }
    if (matr) {
      Double_t* v = b3d->fPnts;
      Double_t buf[3];
      for (UInt_t i = 0; i < b3d->NbPnts(); ++i) {
        buf[0] = v[i * 3]; buf[1] = v[i * 3 + 1]; buf[2] = v[i * 3 + 2];
        matr->LocalToMaster(buf, &v[i * 3]);
      }
    }
    res.reset(RootCsg::ConvertToMesh(*b3d));
  } else {
    auto node = comp->GetBoolNode();
    TGeoHMatrix mleft, mright;
    if (matr) { mleft = *matr; mright = *matr; }
    mleft.Multiply(node->GetLeftMatrix());
    auto left = make_mesh(&mleft, node->GetLeftShape(), err, nseg);
    mright.Multiply(node->GetRightMatrix());
    auto right = make_mesh(&mright, node->GetRightShape(), err, nseg);
    if (!left || !right) return res;
    if (node->IsA() == TGeoUnion::Class())
      res.reset(RootCsg::BuildUnion(left.get(), right.get()));
    else if (node->IsA() == TGeoIntersection::Class())
      res.reset(RootCsg::BuildIntersection(left.get(), right.get()));
    else if (node->IsA() == TGeoSubtraction::Class())
      res.reset(RootCsg::BuildDifference(left.get(), right.get()));
  }
  return res;
}

inline Mesh tessellate(const TGeoShape* shape, int nseg) {
  Mesh r;
  auto t0 = std::chrono::steady_clock::now();
  int old = gGeoManager->GetNsegments();
  gGeoManager->SetNsegments(nseg);
  std::unique_ptr<RootCsg::TBaseMesh> m;
  try {
    m = make_mesh(nullptr, shape, r.err, nseg);
  } catch (std::exception& e) {
    r.err += e.what();
  }
  gGeoManager->SetNsegments(old);
  if (m && m->NumberOfPolys() > 0) {
    r.nv = m->NumberOfVertices();
    r.v.reserve(3 * r.nv);
    for (int i = 0; i < r.nv; ++i) {
      const Double_t* p = m->GetVertex(i);
      r.v.push_back(p[0]); r.v.push_back(p[1]); r.v.push_back(p[2]);
    }
    for (UInt_t k = 0; k < m->NumberOfPolys(); ++k) {
      UInt_t n = m->SizeOfPoly(k);
      if (n < 3) continue;
      r.pol.push_back(n);
      for (UInt_t j = 0; j < n; ++j) r.pol.push_back(m->GetVertexIndex(k, j));
      ++r.npol;
    }
    r.ok = true;
  }
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

struct Node {
  const TGeoNode* node = nullptr;
  const TGeoShape* shape = nullptr;
  int depth = 0, g4depth = 0, parent = -1;
  bool assembly = false, visible = true, has_vis = false, show_daughters = true, wire = false, draw_wire = false;
  double rgba[4] = {1, 1, 1, 1};
  double m[12] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};   // row-major R, then T [cm]
};

class Walker {
public:
  std::vector<Node> nodes;
  void walk(const TGeoNode* top) {
    nodes.clear();
    descend(top, 0, 0, TGeoHMatrix(), -1);
  }

private:
  void descend(const TGeoNode* tn, int depth, int g4depth, const TGeoHMatrix& parentM, int parent) {
    TGeoHMatrix M = parentM;
    if (depth > 0) M.Multiply(tn->GetMatrix());
    Node n;
    n.node = tn; n.depth = depth; n.g4depth = g4depth; n.parent = parent;
    TGeoVolume* tv = tn->GetVolume();
    n.shape = tv->GetShape();
    n.assembly = tv->IsAssembly();
    dd4hep::Volume vol(tv);
    dd4hep::VisAttr vis = vol.visAttributes();
    if (vis.isValid()) {
      n.has_vis = true;
      n.visible = vis.visible();
      n.show_daughters = vis.showDaughters();
      float r, g, b, a;
      vis.argb(a, r, g, b);
      n.rgba[0] = r; n.rgba[1] = g; n.rgba[2] = b; n.rgba[3] = vis.alpha();
      // DDG4's Geant4Converter maps the *line* style to G4 ForceWireframe (drawingStyle is not propagated);
      // mirror that by default so that the drawings match the Geant4/DAWNFILE ones
      n.draw_wire = vis.drawingStyle() == dd4hep::VisAttr::WIREFRAME;
      n.wire = options().honour_drawing_style ? n.draw_wire : vis.lineStyle() == dd4hep::VisAttr::WIREFRAME;
    }
    const Double_t* R = M.GetRotationMatrix();
    const Double_t* T = M.GetTranslation();
    for (int i = 0; i < 9; ++i) n.m[i] = R[i];
    for (int i = 0; i < 3; ++i) n.m[9 + i] = T[i];
    int me = (int)nodes.size();
    nodes.push_back(n);
    bool go = n.assembly ? true : n.show_daughters;   // assemblies are not volumes in Geant4: never stop there
    if (!go) return;
    int child_g4depth = n.assembly ? g4depth : g4depth + 1;
    for (int i = 0; i < tn->GetNdaughters(); ++i)
      descend(tn->GetDaughter(i), depth + 1, child_g4depth, M, me);
  }
};

}  // namespace solidtgeo
