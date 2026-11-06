// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_PolNH3Target: SIDIS polarized NH3 target (scattering chamber with windows, 30 K shield, LHe can, NH3 cup,
// split-coil target magnet with its aluminum boxes, lids and support).  Replaces the old NH3Target stub.
//
// The target is a tree of simple CSG parts, read entirely from XML, so that the two GEMC models are
// both data:  SM2023  (Scientific-Magnetics magnet, GEMC target/solid_SIDIS_target_NH3_SM, 2023-12, default) and
//             Oxford2020 (GEMC target/solid_SIDIS_target_NH3 Original, used by solid_SIDIS_NH3_moved_full).
//
//   <detector id="Target_ID" name="Target" type="SoLID_PolNH3Target" select="Target_NH3_Variant">
//     <position z="TargetCenter_z"/>                     origin of the part tree (hall frame)
//     <variant id="1" name="SM2023">
//       <part name=".." shape="Tube|Cons|Polycone|Box" material=".." vis=".." [dims]>
//         <position x y z/> <rotation x y z/>            placement in the mother part (ACTIVE rotation,
//                                                       R = RotationX(x) * RotationY(y) * RotationZ(z))
//         <zplane z rmin rmax/> ...                      Polycone planes
//         <subtract shape=".." [dims] frame="mother|self"> <position/> <rotation/> </subtract> ...
//         <part ...> ... </part>                         daughters
//       </part>
//     </variant>
//   </detector>
//
// Shapes / dims (G4 conventions, half lengths): Tube rmin rmax dz [phi0 dphi]; Cons rmin1 rmax1 rmin2 rmax2 dz
// [phi0 dphi]; Polycone [phi0 dphi] + <zplane>s; Box dx dy dz.
// A <subtract> with frame="mother" gives its transform in the MOTHER frame of the part (GEMC "Operation:@"):
// the solid is subtracted with the relative transform T_part^-1 * T_sub.  frame="self" (default) gives it in the
// part's own frame.  GEMC rotation "a b c" (passive) = <rotation x="-a" y="-b" z="-c"/> here.
// No sensitive volumes, no readout.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {
/// Value of a variant selector attribute: either a number/expression or the name of an OPTIONAL constant.
/// When the attribute names a constant that is not defined (the top level did not override it), use `def`.
/// (A fragment must not define its own default: with DD4hep 1.38 a second definition of a constant is ignored by
/// the constants table but overwrites the expression evaluator, i.e. the LAST definition wins in expressions.)
int selector_value(Detector& description, xml_comp_t x, const xml::XmlChar* attr, int def, const char* who)
{
  if (!x.hasAttr(attr)) {
    return def;
  }
  const std::string s = x.attr<std::string>(attr);
  bool identifier     = !s.empty() && (std::isalpha(static_cast<unsigned char>(s[0])) || s[0] == '_');
  for (char c : s) {
    identifier = identifier && (std::isalnum(static_cast<unsigned char>(c)) || c == '_');
  }
  if (identifier && description.constants().find(s) == description.constants().end()) {
    printout(INFO, who, "selector constant %s not defined: using the default %d", s.c_str(), def);
    return def;
  }
  return _toInt(s);
}
} // namespace

namespace {

Transform3D part_transform(xml_comp_t x)
{
  double px = 0., py = 0., pz = 0., rx = 0., ry = 0., rz = 0.;
  if (x.hasChild(_U(position))) {
    xml_comp_t p = x.child(_U(position));
    px           = p.x(0.);
    py           = p.y(0.);
    pz           = p.z(0.);
  }
  if (x.hasChild(_U(rotation))) {
    xml_comp_t r = x.child(_U(rotation));
    rx           = r.x(0.);
    ry           = r.y(0.);
    rz           = r.z(0.);
  }
  const Rotation3D rot = Rotation3D(RotationX(rx)) * Rotation3D(RotationY(ry)) * Rotation3D(RotationZ(rz));
  return Transform3D(rot, Position(px, py, pz));
}

Solid make_shape(xml_comp_t x)
{
  const std::string shape = x.attr<std::string>(_Unicode(shape));
  const double phi0       = getAttrOrDefault<double>(x, _Unicode(phi0), 0.);
  const double dphi       = getAttrOrDefault<double>(x, _Unicode(dphi), 2. * M_PI);
  if (shape == "Tube") {
    return Tube(getAttrOrDefault<double>(x, _U(rmin), 0.), x.rmax(), x.dz(), phi0, phi0 + dphi);
  }
  if (shape == "Cons") {
    return ConeSegment(x.dz(), x.attr<double>(_Unicode(rmin1)), x.attr<double>(_Unicode(rmax1)),
                       x.attr<double>(_Unicode(rmin2)), x.attr<double>(_Unicode(rmax2)), phi0, phi0 + dphi);
  }
  if (shape == "Polycone") {
    std::vector<double> z, rmin, rmax;
    for (xml_coll_t iz(x, _Unicode(zplane)); iz; ++iz) {
      xml_comp_t x_z = iz;
      z.push_back(x_z.z());
      rmin.push_back(x_z.rmin());
      rmax.push_back(x_z.rmax());
    }
    return Polycone(phi0, dphi, rmin, rmax, z);
  }
  if (shape == "Box") {
    return Box(x.dx(), x.dy(), x.dz());
  }
  throw std::runtime_error("SoLID_PolNH3Target: unknown shape " + shape);
}

int build_parts(Detector& description, xml_comp_t x_parent, Volume mother, const std::string& prefix)
{
  int n = 0;
  for (xml_coll_t ip(x_parent, _Unicode(part)); ip; ++ip) {
    xml_comp_t x_p          = ip;
    const std::string pname = prefix + "_" + x_p.nameStr();
    const Transform3D tr    = part_transform(x_p);
    Solid solid             = make_shape(x_p);
    for (xml_coll_t is(x_p, _Unicode(subtract)); is; ++is) {
      xml_comp_t x_s        = is;
      const std::string fr  = getAttrOrDefault<std::string>(x_s, _Unicode(frame), "self");
      const Transform3D tsub = (fr == "mother") ? tr.Inverse() * part_transform(x_s) : part_transform(x_s);
      solid                 = SubtractionSolid(solid, make_shape(x_s), tsub);
    }
    Volume vol(pname, solid, description.material(x_p.materialStr()));
    vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_p, _Unicode(vis), "InvisibleWithDaughters"));
    n += 1 + build_parts(description, x_p, vol, prefix);
    mother.placeVolume(vol, tr);
  }
  return n;
}

} // namespace

static Ref_t create_nh3_target(Detector& description, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  DetElement det(det_name, x_det.id());
  const double zc = x_det.hasChild(_U(position)) ? xml_comp_t(x_det.child(_U(position))).z(0.) : 0.;
  const int sel   = selector_value(description, x_det, _Unicode(select), 1, "SoLID_PolNH3Target");

  Assembly assembly(det_name + "_assembly");
  assembly.setVisAttributes(description.invisible());
  Assembly tree(det_name + "_parts");
  tree.setVisAttributes(description.invisible());

  bool found = false;
  for (xml_coll_t iv(x_det, _Unicode(variant)); iv; ++iv) {
    xml_comp_t x_v = iv;
    if (x_v.id() != sel) {
      continue;
    }
    found       = true;
    const int n = build_parts(description, x_v, tree, det_name);
    printout(INFO, "SoLID_PolNH3Target", "%s: variant %d (%s), %d volumes, origin at hall z %.3f cm", det_name.c_str(),
             sel, x_v.nameStr().c_str(), n, zc / cm);
  }
  if (!found) {
    throw std::runtime_error("SoLID_PolNH3Target: no <variant id=\"" + std::to_string(sel) + "\">");
  }
  assembly.placeVolume(tree, Position(0., 0., zc));

  PlacedVolume pv = description.pickMotherVolume(det).placeVolume(assembly);
  det.setPlacement(pv);
  return det;
}

DECLARE_DETELEMENT(SoLID_PolNH3Target, create_nh3_target)
