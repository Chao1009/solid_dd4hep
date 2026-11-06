// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_PolHe3Target: SIDIS polarized 3He target (GEMC target/solid_SIDIS_target_He3 + _collimator + _coil).
//
//   <detector id="Target_ID" name="Target" type="SoLID_PolHe3Target">
//     <position z="TargetCenter_z"/>
//     <cell length=".." rin=".." wall=".." window=".." material="Target_GE180" gas="Target_He3_10amg"
//           envelope_r=".." envelope_material="Vacuum"/>
//     <coolgas enable="1" material="Target_N2Gas" rmax="..">
//       <slab name="upstream" rmin="0" zmin=".." zmax=".."/> ...        (z relative to the cell center)
//     </coolgas>
//     <collimators select="expr">                                       (0 = none)
//       <set id="1" name="2022">
//         <collimator name="Upstream" material=".." [mother="<slab name>"]>
//           <zplane z=".." rmin=".." rmax=".."/> ...                    (polycone, z relative to the cell center)
//         </collimator> ...
//       </set> ...
//     </collimators>
//     <coil_pair name=".." axis="x|y|z" R=".." width=".." material="Copper" [enable="1"]/>
//   </detector>
//
// Cell (GEMC TRGB/TCEL/TLHE/TLW1/TLW2): vacuum envelope Tube(0, envelope_r, length/2 + window) holding the
// glass side wall Tube(rin, rin + wall, length/2), the gas Tube(0, rin, length/2) and two end windows
// Tube(0, rin, window/2) at +-(length/2 + window/2).
// Coil pair: two square-section rings Tube(R - width/2, R + width/2, width/2) at +-R/2 along "axis"
// (Helmholtz spacing) around the cell center.  No sensitive volumes, no readout.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <map>
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

static Ref_t create_he3_target(Detector& description, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  DetElement det(det_name, x_det.id());
  const double zc = x_det.hasChild(_U(position)) ? xml_comp_t(x_det.child(_U(position))).z(0.) : 0.;

  Assembly assembly(det_name + "_assembly");
  assembly.setVisAttributes(description.invisible());

  // ---------------------------------------------------------------- glass cell
  xml_comp_t x_cell   = x_det.child(_Unicode(cell));
  const double length = x_cell.attr<double>(_Unicode(length));
  const double rin    = x_cell.attr<double>(_Unicode(rin));
  const double wall   = x_cell.attr<double>(_Unicode(wall));
  const double window = x_cell.attr<double>(_Unicode(window));
  const double env_r  = x_cell.attr<double>(_Unicode(envelope_r));
  if (env_r < rin + wall) {
    throw std::runtime_error("SoLID_PolHe3Target: envelope_r smaller than the cell outer radius");
  }
  Material glass = description.material(x_cell.materialStr());
  Volume env_vol(det_name + "_CellEnvelope", Tube(0., env_r, 0.5 * length + window),
                 description.material(getAttrOrDefault<std::string>(x_cell, _Unicode(envelope_material), "Vacuum")));
  env_vol.setVisAttributes(description, "InvisibleWithDaughters");
  Volume wall_vol(det_name + "_CellWall", Tube(rin, rin + wall, 0.5 * length), glass);
  wall_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_cell, _Unicode(vis), "Target_GlassVis"));
  Volume gas_vol(det_name + "_Gas", Tube(0., rin, 0.5 * length), description.material(x_cell.attr<std::string>(_Unicode(gas))));
  gas_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_cell, _Unicode(gas_vis), "Target_GasVis"));
  Volume win_vol(det_name + "_CellWindow", Tube(0., rin, 0.5 * window), glass);
  win_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_cell, _Unicode(vis), "Target_GlassVis"));
  env_vol.placeVolume(wall_vol);
  env_vol.placeVolume(gas_vol);
  env_vol.placeVolume(win_vol, Position(0., 0., -0.5 * (length + window)));
  env_vol.placeVolume(win_vol, Position(0., 0., +0.5 * (length + window)));
  assembly.placeVolume(env_vol, Position(0., 0., zc));
  printout(INFO, "SoLID_PolHe3Target", "%s cell: %.2f cm, r_in %.3f cm, wall %.3f mm, windows %.4f mm, center z %.3f cm",
           det_name.c_str(), length / cm, rin / cm, wall / mm, window / mm, zc / cm);

  // ---------------------------------------------------------------- cooling gas slabs
  std::map<std::string, std::pair<Volume, double>> slabs; // name -> (volume, center z relative to the cell)
  if (x_det.hasChild(_Unicode(coolgas))) {
    xml_comp_t x_cg = x_det.child(_Unicode(coolgas));
    if (getAttrOrDefault<int>(x_cg, _Unicode(enable), 1) != 0) {
      Material cg_mat   = description.material(x_cg.materialStr());
      const double rmax = x_cg.attr<double>(_Unicode(rmax));
      for (xml_coll_t is(x_cg, _Unicode(slab)); is; ++is) {
        xml_comp_t x_s     = is;
        const double zmin  = x_s.attr<double>(_Unicode(zmin));
        const double zmax  = x_s.attr<double>(_Unicode(zmax));
        const double rmin  = getAttrOrDefault<double>(x_s, _Unicode(rmin), 0.);
        const std::string n = x_s.nameStr();
        Volume v(det_name + "_CoolGas_" + n, Tube(rmin, rmax, 0.5 * (zmax - zmin)), cg_mat);
        v.setVisAttributes(description, getAttrOrDefault<std::string>(x_cg, _Unicode(vis), "Target_CoolGasVis"));
        slabs[n] = {v, 0.5 * (zmin + zmax)};
        assembly.placeVolume(v, Position(0., 0., zc + 0.5 * (zmin + zmax)));
        printout(INFO, "SoLID_PolHe3Target", "%s cooling gas %s: %s r %.2f..%.2f cm, hall z %.3f..%.3f cm",
                 det_name.c_str(), n.c_str(), x_cg.materialStr().c_str(), rmin / cm, rmax / cm, (zc + zmin) / cm,
                 (zc + zmax) / cm);
      }
    }
  }

  // ---------------------------------------------------------------- collimators (selectable sets)
  if (x_det.hasChild(_Unicode(collimators))) {
    xml_comp_t x_cols = x_det.child(_Unicode(collimators));
    const int sel     = selector_value(description, x_cols, _Unicode(select), 1, "SoLID_PolHe3Target");
    bool found        = (sel == 0);
    for (xml_coll_t iset(x_cols, _Unicode(set)); iset; ++iset) {
      xml_comp_t x_set = iset;
      if (x_set.id() != sel) {
        continue;
      }
      found = true;
      for (xml_coll_t ic(x_set, _Unicode(collimator)); ic; ++ic) {
        xml_comp_t x_c = ic;
        std::vector<double> z, rmin, rmax;
        for (xml_coll_t iz(x_c, _Unicode(zplane)); iz; ++iz) {
          xml_comp_t x_z = iz;
          z.push_back(x_z.z());
          rmin.push_back(x_z.rmin());
          rmax.push_back(x_z.rmax());
        }
        const std::string cname = det_name + "_Collimator" + x_c.nameStr();
        Volume cv(cname, Polycone(0., 2. * M_PI, rmin, rmax, z), description.material(x_c.materialStr()));
        cv.setVisAttributes(description, getAttrOrDefault<std::string>(x_c, _Unicode(vis), "Target_CollimatorVis"));
        const std::string mother = getAttrOrDefault<std::string>(x_c, _Unicode(mother), "");
        if (!mother.empty() && slabs.count(mother)) {
          slabs[mother].first.placeVolume(cv, Position(0., 0., -slabs[mother].second));
        } else {
          assembly.placeVolume(cv, Position(0., 0., zc));
        }
        printout(INFO, "SoLID_PolHe3Target", "%s (set %s): hall z %.3f..%.3f cm, rmin %.3f..%.3f, rmax %.3f..%.3f cm%s",
                 cname.c_str(), x_set.nameStr().c_str(), (zc + z.front()) / cm, (zc + z.back()) / cm,
                 rmin.front() / cm, rmin.back() / cm, rmax.front() / cm, rmax.back() / cm,
                 mother.empty() ? "" : (" in cooling gas " + mother).c_str());
      }
    }
    if (!found) {
      throw std::runtime_error("SoLID_PolHe3Target: no collimator <set id=\"" + std::to_string(sel) + "\">");
    }
  }

  // ---------------------------------------------------------------- Helmholtz coil pairs
  for (xml_coll_t ip(x_det, _Unicode(coil_pair)); ip; ++ip) {
    xml_comp_t x_p = ip;
    if (getAttrOrDefault<int>(x_p, _Unicode(enable), 1) == 0) {
      continue;
    }
    const double R    = x_p.attr<double>(_Unicode(R));
    const double w    = x_p.attr<double>(_Unicode(width));
    const std::string axis = x_p.attr<std::string>(_Unicode(axis));
    const std::string pname = det_name + "_Coil" + x_p.nameStr();
    Volume cv(pname, Tube(R - 0.5 * w, R + 0.5 * w, 0.5 * w), description.material(x_p.materialStr()));
    cv.setVisAttributes(description, getAttrOrDefault<std::string>(x_p, _Unicode(vis), "Target_CoilVis"));
    for (int sgn : {-1, +1}) {
      const double d = sgn * 0.5 * R;
      if (axis == "z") {
        assembly.placeVolume(cv, Position(0., 0., zc + d));
      } else if (axis == "x") {
        // ring axis along hall x: local z -> hall x (GEMC rot "0 90 0" = active RotationY(-90 deg))
        assembly.placeVolume(cv, Transform3D(RotationY(-90. * degree), Position(d, 0., zc)));
      } else if (axis == "y") {
        // ring axis along hall y: local z -> hall y
        assembly.placeVolume(cv, Transform3D(RotationX(-90. * degree), Position(0., d, zc)));
      } else {
        throw std::runtime_error("SoLID_PolHe3Target: coil_pair axis must be x, y or z");
      }
    }
    printout(INFO, "SoLID_PolHe3Target", "%s: axis %s, R %.2f cm, %.1fx%.1f cm %s rings at +-%.2f cm", pname.c_str(),
             axis.c_str(), R / cm, w / cm, w / cm, x_p.materialStr().c_str(), 0.5 * R / cm);
  }

  PlacedVolume pv = description.pickMotherVolume(det).placeVolume(assembly);
  det.setPlacement(pv);
  return det;
}

DECLARE_DETELEMENT(SoLID_PolHe3Target, create_he3_target)
