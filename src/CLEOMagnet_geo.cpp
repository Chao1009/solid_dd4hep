// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_CLEOMagnet: passive CLEO-II solenoid magnet (iron yoke, or coil + support shell + cryostat).
//
// One generic, fully XML-driven builder. A <detector> of this type is a list of coaxial parts, each of
// them a solid of revolution with a round (circular) or octagonal (polyhedral) inner and outer surface,
// plus optional rings of discrete bars (the spacer bars between the two CLEO-II yoke layers).
// The yoke (MagnetYoke) and the coil/cryostat (MagnetCoil) are two <detector> elements of this same type.
//
//   <detector id="MagnetYoke_ID" name="MagnetYoke" type="SoLID_CLEOMagnet" vis="InvisibleWithDaughters">
//     <octagon nsides="8" phi0="-22.5*deg"/>     <!-- default polygon for "oct" surfaces; radii are APOTHEMS -->
//     <position x="0" y="0" z="0"/>               <!-- optional placement of the whole magnet (default origin) -->
//     <part name="InnerOctagon" inner="oct" outer="oct" material="Magnet_YokeSteel" vis="Magnet_YokeVis"
//           enable="Magnet_SomeSwitch">          <!-- optional: built only if the expression is != 0 -->
//       <zplane z="-266.70*cm" rmin="180.34*cm" rmax="212.62*cm"/>
//       ...                                       <!-- >= 2 planes, z non-decreasing (steps allowed) -->
//     </part>
//     <bars name="SpacerBarsUp" n="8" phi0="0*deg" rmin="212.62*cm" rmax="221.51*cm" width="149.86*cm"
//           zmin="-266.70*cm" zmax="-232.66*cm" material="Magnet_YokeSteel" vis="Magnet_YokeVis" enable="..."/>
//     <material_map from="Magnet_YokeSteel" to="Iron"/>  <!-- optional material substitution (also takes enable);
//                                                     from="*" maps every material (e.g. a GEMC "kill" absorber) -->
//   </detector>
//
// Solid per part (inner/outer surface kind):
//   oct/oct     -> Polyhedra(nsides, phi0, 2pi, z, rmin, rmax)
//   round/round -> Polycone(0, 2pi, rmin, rmax, z)
//   round/oct   -> Polyhedra(rmin = 0, rmax) - Polycone(0, rmin_ext)     (inner round, outer octagonal)
//   oct/round   -> Polycone(0, rmax) - Polyhedra(rmin = 0, rmin_ext)     (inner octagonal, outer round)
// For the booleans the subtracted solid has its first/last z planes pushed outward by `boolean_eps`
// (default 0.01 cm, attribute of <detector>) with the radius extrapolated linearly, so that no faces coincide.
// Polyhedra radii are apothems (distance to the flats), as in TGeoPgon and G4Polyhedra. With phi0 = -22.5 deg
// an octagon has its flats centred at phi = 0, 45, 90, ... deg (flats normal to x and y, as the real CLEO-II
// yoke; GEMC: Pgon phi0 = 0 with rot "0 0 22.5deg" = ACTIVE RotationZ(-22.5 deg)).
//
// Bars: n boxes of radial thickness (rmax - rmin), tangential width `width` and z length (zmax - zmin), centred
// at r = (rmin + rmax)/2 and phi_k = phi0 + k*360/n, ACTIVE rotation RotationZ(phi_k).
//
// Passive detector: no readout, no sensitive volume. All numbers come from XML.

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <cmath>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

  enum class Surface { Round, Oct };

  Surface parse_surface(const std::string& s, const std::string& part)
  {
    if (s == "round" || s == "circle" || s == "circular")
      return Surface::Round;
    if (s == "oct" || s == "octagon" || s == "polygon" || s == "poly")
      return Surface::Oct;
    throw std::runtime_error("SoLID_CLEOMagnet: part " + part + ": unknown surface kind '" + s +
                             "' (use round|oct)");
  }

  // Expand the identifiers of an expression that are compact constants by the constant's value as registered in
  // the Detector's constant table. Reason: if a constant is ever defined twice, DD4hep keeps the FIRST value in
  // Detector::constants() but the expression evaluator holds the LAST one ("Overwriting variable"); expanding from
  // the constant table makes the switch deterministic. (Variants are separate fragments, see docs/design.md sec. 4.)
  std::string expand_constants(const Detector& desc, std::string expr, int depth = 0)
  {
    static const std::regex ident(R"([A-Za-z_][A-Za-z0-9_]*)");
    const auto&             consts = desc.constants();
    std::string             out;
    auto                    begin = std::sregex_iterator(expr.begin(), expr.end(), ident);
    std::size_t             last  = 0;
    for (auto it = begin; it != std::sregex_iterator(); ++it) {
      const std::string name = it->str();
      out += expr.substr(last, it->position() - last);
      auto c = consts.find(name);
      if (c != consts.end() && depth < 8) {
        out += "(" + expand_constants(desc, c->second->GetTitle(), depth + 1) + ")";
      } else {
        out += name;
      }
      last = it->position() + it->length();
    }
    out += expr.substr(last);
    return out;
  }

  // evaluate an optional "enable" expression (default: enabled)
  bool is_enabled(const Detector& desc, xml_comp_t x)
  {
    if (!x.hasAttr(_Unicode(enable)))
      return true;
    const double v = _toDouble(expand_constants(desc, x.attr<std::string>(_Unicode(enable))));
    return std::fabs(v) > 0.5;
  }

  // copy of r(z) with the first and last plane pushed outward by eps, radius extrapolated linearly
  void extend_planes(const std::vector<double>& z, const std::vector<double>& r, double eps,
                     std::vector<double>& z_ext, std::vector<double>& r_ext)
  {
    z_ext = z;
    r_ext = r;
    const std::size_t n = z.size();
    auto slope = [&](std::size_t i0, std::size_t i1) {
      const double dz = z[i1] - z[i0];
      return (std::fabs(dz) > 1e-9) ? (r[i1] - r[i0]) / dz : 0.0;
    };
    // first non-degenerate segment from the front / back
    std::size_t f = 1;
    while (f < n - 1 && std::fabs(z[f] - z[0]) < 1e-9)
      ++f;
    std::size_t b = n - 2;
    while (b > 0 && std::fabs(z[n - 1] - z[b]) < 1e-9)
      --b;
    const double s0 = (std::fabs(z[f] - z[0]) < 1e-9) ? 0.0 : slope(0, f);
    const double s1 = (std::fabs(z[n - 1] - z[b]) < 1e-9) ? 0.0 : slope(b, n - 1);
    z_ext[0]     = z[0] - eps;
    r_ext[0]     = std::max(0.0, r[0] - s0 * eps);
    z_ext[n - 1] = z[n - 1] + eps;
    r_ext[n - 1] = std::max(0.0, r[n - 1] + s1 * eps);
  }

  Solid make_part_solid(const std::string& name, Surface inner, Surface outer, int nsides, double phi0,
                        const std::vector<double>& z, const std::vector<double>& rmin,
                        const std::vector<double>& rmax, double eps)
  {
    const double two_pi = 2.0 * M_PI;
    if (inner == Surface::Oct && outer == Surface::Oct) {
      return Polyhedra(name, nsides, phi0, two_pi, z, rmin, rmax);
    }
    if (inner == Surface::Round && outer == Surface::Round) {
      return Polycone(name, 0.0, two_pi, rmin, rmax, z);
    }
    std::vector<double> zero(z.size(), 0.0);
    std::vector<double> z_ext, rin_ext;
    extend_planes(z, rmin, eps, z_ext, rin_ext);
    std::vector<double> zero_ext(z_ext.size(), 0.0);
    if (inner == Surface::Round) { // outer octagonal
      Polyhedra outer_s(name + "_outer", nsides, phi0, two_pi, z, zero, rmax);
      Polycone  hole(name + "_bore", 0.0, two_pi, zero_ext, rin_ext, z_ext);
      return SubtractionSolid(name, outer_s, hole);
    }
    // inner octagonal, outer round
    Polycone  outer_s(name + "_outer", 0.0, two_pi, zero, rmax, z);
    Polyhedra hole(name + "_bore", nsides, phi0, two_pi, z_ext, zero_ext, rin_ext);
    return SubtractionSolid(name, outer_s, hole);
  }

} // namespace

static Ref_t create_cleo_magnet(Detector& desc, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t   x_det    = e;
  const int   det_id   = x_det.id();
  std::string det_name = x_det.nameStr();
  DetElement  det(det_name, det_id);

  // default polygon for octagonal surfaces
  int    def_nsides = 8;
  double def_phi0   = -22.5 * dd4hep::deg;
  if (x_det.hasChild(_Unicode(octagon))) {
    xml_comp_t x_oct = x_det.child(_Unicode(octagon));
    def_nsides       = getAttrOrDefault<int>(x_oct, _Unicode(nsides), def_nsides);
    def_phi0         = getAttrOrDefault<double>(x_oct, _Unicode(phi0), def_phi0);
  }
  const double eps = getAttrOrDefault<double>(x_det, _Unicode(boolean_eps), 0.01 * dd4hep::cm);

  // optional material substitutions (e.g. GEMC parity: G4_Fe instead of the yoke steel)
  std::map<std::string, std::string> mat_map;
  for (xml_coll_t im(x_det, _Unicode(material_map)); im; ++im) {
    xml_comp_t x_map = im;
    if (is_enabled(desc, x_map)) {
      mat_map[x_map.attr<std::string>(_Unicode(from))] = x_map.attr<std::string>(_Unicode(to));
    }
  }
  auto material = [&](const std::string& name) {
    auto it = mat_map.find(name);
    if (it == mat_map.end())
      it = mat_map.find("*");
    const std::string& use = (it == mat_map.end()) ? name : it->second;
    Material mat = desc.material(use);
    if (!mat.isValid()) {
      throw std::runtime_error("SoLID_CLEOMagnet: unknown material " + use);
    }
    return mat;
  };

  Assembly assembly(det_name + "_assembly");
  int      nparts = 0;

  // ---- coaxial parts
  for (xml_coll_t ip(x_det, _Unicode(part)); ip; ++ip) {
    xml_comp_t  x_part = ip;
    std::string pname  = x_part.nameStr();
    if (!is_enabled(desc, x_part)) {
      printout(DEBUG, "SoLID_CLEOMagnet", "%s: part %s disabled", det_name.c_str(), pname.c_str());
      continue;
    }
    const Surface inner =
        parse_surface(getAttrOrDefault<std::string>(x_part, _Unicode(inner), "round"), pname);
    const Surface outer =
        parse_surface(getAttrOrDefault<std::string>(x_part, _Unicode(outer), "round"), pname);
    const int    nsides = getAttrOrDefault<int>(x_part, _Unicode(nsides), def_nsides);
    const double phi0   = getAttrOrDefault<double>(x_part, _Unicode(phi0), def_phi0);

    std::vector<double> z, rmin, rmax;
    for (xml_coll_t iz(x_part, _Unicode(zplane)); iz; ++iz) {
      xml_comp_t x_z = iz;
      z.push_back(x_z.z());
      rmin.push_back(x_z.rmin());
      rmax.push_back(x_z.rmax());
    }
    if (z.size() < 2) {
      throw std::runtime_error("SoLID_CLEOMagnet: part " + pname + " needs at least 2 zplanes");
    }
    for (std::size_t i = 0; i < z.size(); ++i) {
      if (i > 0 && z[i] < z[i - 1]) {
        throw std::runtime_error("SoLID_CLEOMagnet: part " + pname + ": zplanes must be ordered in z");
      }
      if (rmin[i] < 0 || rmax[i] <= rmin[i]) {
        throw std::runtime_error("SoLID_CLEOMagnet: part " + pname + ": need 0 <= rmin < rmax at every plane");
      }
    }

    Solid  solid = make_part_solid(det_name + "_" + pname, inner, outer, nsides, phi0, z, rmin, rmax, eps);
    Volume vol(det_name + "_" + pname, solid, material(x_part.materialStr()));
    vol.setVisAttributes(desc.visAttributes(getAttrOrDefault<std::string>(x_part, _Unicode(vis), "")));
    assembly.placeVolume(vol);
    ++nparts;
    printout(DEBUG, "SoLID_CLEOMagnet", "%s: part %-24s %s/%s, %zu planes, z %.2f..%.2f cm, %s", det_name.c_str(),
             pname.c_str(), inner == Surface::Oct ? "oct" : "round", outer == Surface::Oct ? "oct" : "round",
             z.size(), z.front() / dd4hep::cm, z.back() / dd4hep::cm, x_part.materialStr().c_str());
  }

  // ---- rings of discrete bars (spacer bars between yoke layers)
  for (xml_coll_t ib(x_det, _Unicode(bars)); ib; ++ib) {
    xml_comp_t  x_bars = ib;
    std::string bname  = x_bars.nameStr();
    if (!is_enabled(desc, x_bars)) {
      printout(DEBUG, "SoLID_CLEOMagnet", "%s: bars %s disabled", det_name.c_str(), bname.c_str());
      continue;
    }
    const int    n      = x_bars.attr<int>(_Unicode(n));
    const double phi0   = getAttrOrDefault<double>(x_bars, _Unicode(phi0), 0.0);
    const double r0     = x_bars.rmin();
    const double r1     = x_bars.rmax();
    const double width  = x_bars.width();
    const double z0     = x_bars.zmin();
    const double z1     = x_bars.zmax();
    if (n <= 0 || r1 <= r0 || z1 <= z0 || width <= 0) {
      throw std::runtime_error("SoLID_CLEOMagnet: bars " + bname + ": invalid dimensions");
    }
    Box    box(det_name + "_" + bname, (r1 - r0) / 2.0, width / 2.0, (z1 - z0) / 2.0);
    Volume vol(det_name + "_" + bname, box, material(x_bars.materialStr()));
    vol.setVisAttributes(desc.visAttributes(getAttrOrDefault<std::string>(x_bars, _Unicode(vis), "")));
    const double rc = (r0 + r1) / 2.0;
    const double zc = (z0 + z1) / 2.0;
    for (int k = 0; k < n; ++k) {
      const double phi = phi0 + k * 2.0 * M_PI / n;
      assembly.placeVolume(vol, k, Transform3D(RotationZ(phi), Position(rc * std::cos(phi), rc * std::sin(phi), zc)));
    }
    ++nparts;
    printout(DEBUG, "SoLID_CLEOMagnet", "%s: bars %s x%d, r %.2f..%.2f, z %.2f..%.2f cm", det_name.c_str(),
             bname.c_str(), n, r0 / dd4hep::cm, r1 / dd4hep::cm, z0 / dd4hep::cm, z1 / dd4hep::cm);
  }

  if (nparts == 0) {
    printout(WARNING, "SoLID_CLEOMagnet", "%s: no part enabled, nothing is built", det_name.c_str());
  }

  // ---- placement in the mother (world)
  Position pos(0, 0, 0);
  if (x_det.hasChild(_Unicode(position))) {
    xml_comp_t x_pos = x_det.child(_Unicode(position));
    pos              = Position(getAttrOrDefault<double>(x_pos, _Unicode(x), 0.0),
                                getAttrOrDefault<double>(x_pos, _Unicode(y), 0.0),
                                getAttrOrDefault<double>(x_pos, _Unicode(z), 0.0));
  }
  assembly.setVisAttributes(desc.invisible());
  if (x_det.hasAttr(_Unicode(vis))) {
    assembly.setVisAttributes(desc.visAttributes(x_det.visStr()));
  }
  Volume        mother = desc.pickMotherVolume(det);
  PlacedVolume  pv     = mother.placeVolume(assembly, Transform3D(RotationZ(0.), pos));
  pv.addPhysVolID("system", det_id);
  det.setPlacement(pv);
  printout(INFO, "SoLID_CLEOMagnet", "%s (id %d): %d parts/bar sets built", det_name.c_str(), det_id, nparts);
  return det;
}

// SoLID_CLEOMagnet is the name used by the compact files; SoLID_Solenoid is the alias listed in docs/design.md.
DECLARE_DETELEMENT(SoLID_CLEOMagnet, create_cleo_magnet)
DECLARE_DETELEMENT(SoLID_Solenoid, create_cleo_magnet)
