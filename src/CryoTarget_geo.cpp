// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_CryoTarget: liquid H2/D2 cryotarget cell, optionally inside an envelope (J/psi: GEMC TACH, air)
// or a vacuum enclosure (PVDIS: GEMC TACH aluminum + TACV vacuum with beam openings).
//
//   <detector id="Target_ID" name="Target" type="SoLID_CryoTarget">
//     <position z="TargetCenter_z"/>                                      cell center (hall frame)
//     <cell radius=".." length=".." wall=".." window=".." material="Aluminum" fill="Target_LH2"
//           [fill_select="expr"]> [<fill id="1" material=".."/> ...] </cell>
//     <envelope rmax=".." halfz=".." material="Air"/>                      optional (no enclosure)
//     <enclosure rmax=".." halfz=".." wall=".." endwall=".." material="Aluminum" fill="Vacuum" [dz=".."]>
//       <opening side="upstream|downstream" r=".."/>                       vacuum hole through the end wall
//       <cavity rmin=".." rmax=".." halfz=".." z=".."/>                    slot cut from wall + vacuum
//     </enclosure>
//   </detector>
//
// Cell = solid wall-material cylinder r = radius + wall, half length = length/2 + window, with the liquid
// (r = radius, half length = length/2) as daughter: side wall "wall" thick, end windows "window" thick
// (GEMC TALU/TAH2 and TACL/TAW1/TAW2/TALH give the same material map).
// The <cavity> (cut from the enclosure wall AND its vacuum, z relative to the enclosure center) leaves room for
// a volume that another detector places in the world (the PVDIS baffle plate 1 slot).
// No sensitive volumes, no readout.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <cmath>
#include <cctype>
#include <stdexcept>
#include <string>

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

static Ref_t create_cryotarget(Detector& description, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  DetElement det(det_name, x_det.id());

  const double zc = x_det.hasChild(_U(position)) ? xml_comp_t(x_det.child(_U(position))).z(0.) : 0.;

  // ---------------------------------------------------------------- cell
  xml_comp_t x_cell = x_det.child(_Unicode(cell));
  const double radius = x_cell.attr<double>(_Unicode(radius));
  const double length = x_cell.attr<double>(_Unicode(length));
  const double wall   = x_cell.attr<double>(_Unicode(wall));
  const double window = x_cell.attr<double>(_Unicode(window));
  std::string fill    = getAttrOrDefault<std::string>(x_cell, _Unicode(fill), "Target_LH2");
  if (x_cell.hasAttr(_Unicode(fill_select))) {
    const int sel = selector_value(description, x_cell, _Unicode(fill_select), 1, "SoLID_CryoTarget");
    bool found    = false;
    for (xml_coll_t i(x_cell, _Unicode(fill)); i; ++i) {
      xml_comp_t x_f = i;
      if (x_f.id() == sel) {
        fill  = x_f.materialStr();
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error("SoLID_CryoTarget: no <fill id=\"" + std::to_string(sel) + "\"> in <cell>");
    }
  }
  const double cell_rmax  = radius + wall;
  const double cell_halfz = 0.5 * length + window;
  Volume cell_vol(det_name + "_Cell", Tube(0., cell_rmax, cell_halfz), description.material(x_cell.materialStr()));
  cell_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_cell, _Unicode(vis), "Target_CellVis"));
  Volume liq_vol(det_name + "_Liquid", Tube(0., radius, 0.5 * length), description.material(fill));
  liq_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_cell, _Unicode(fill_vis), "Target_LiquidVis"));
  cell_vol.placeVolume(liq_vol);
  printout(INFO, "SoLID_CryoTarget", "%s cell: %s %.3f cm long, r %.4f cm, wall %.4f mm, windows %.4f mm (%s), center z %.3f cm",
           det_name.c_str(), fill.c_str(), length / cm, radius / cm, wall / mm, window / mm,
           x_cell.materialStr().c_str(), zc / cm);

  Assembly assembly(det_name + "_assembly");
  assembly.setVisAttributes(description.invisible());

  if (x_det.hasChild(_Unicode(enclosure))) {
    // -------------------------------------------------------------- vacuum enclosure
    xml_comp_t x_enc     = x_det.child(_Unicode(enclosure));
    const double rmax    = x_enc.attr<double>(_Unicode(rmax));
    const double halfz   = x_enc.attr<double>(_Unicode(halfz));
    const double ewall   = x_enc.attr<double>(_Unicode(wall));
    const double endwall = x_enc.attr<double>(_Unicode(endwall));
    const double dz      = getAttrOrDefault<double>(x_enc, _Unicode(dz), 0.);   // enclosure center - cell center
    const std::string fmat = getAttrOrDefault<std::string>(x_enc, _Unicode(fill), "Vacuum");
    if (std::abs(dz) + cell_halfz > halfz - endwall || cell_rmax > rmax - ewall) {
      throw std::runtime_error("SoLID_CryoTarget: the cell does not fit inside the enclosure vacuum");
    }

    Solid enc_s = Tube(0., rmax, halfz);
    Solid vac_s = Tube(0., rmax - ewall, halfz - endwall);
    int ncav    = 0;
    for (xml_coll_t ic(x_enc, _Unicode(cavity)); ic; ++ic, ++ncav) {
      xml_comp_t x_c = ic;
      double crmin = 0., crmax = 0., chalfz = 0., cz = 0.;
      try {
        crmin  = x_c.attr<double>(_Unicode(rmin));
        crmax  = x_c.attr<double>(_Unicode(rmax));
        chalfz = x_c.attr<double>(_Unicode(halfz));
        cz     = x_c.attr<double>(_Unicode(z));
      } catch (const std::exception& ex) {
        printout(ERROR, "SoLID_CryoTarget",
                 "cannot evaluate the <cavity> of %s (%s). The PVDIS cavity uses the baffle constants "
                 "(BafflePlate1_*, BafflePlate_halfz): include compact/baffle/baffle_PVDIS.xml BEFORE the target.",
                 det_name.c_str(), ex.what());
        throw;
      }
      if (std::abs(cz) + chalfz > halfz - endwall || crmax > rmax - ewall) {
        throw std::runtime_error("SoLID_CryoTarget: cavity must lie inside the enclosure vacuum");
      }
      if (std::abs(cz + dz) < chalfz + cell_halfz && crmin < cell_rmax) {
        throw std::runtime_error("SoLID_CryoTarget: cavity intersects the target cell");
      }
      Tube slot(crmin, crmax, chalfz);
      enc_s = SubtractionSolid(enc_s, slot, Position(0., 0., cz));
      vac_s = SubtractionSolid(vac_s, slot, Position(0., 0., cz));
      printout(INFO, "SoLID_CryoTarget", "%s enclosure cavity: r %.4f..%.4f cm, hall z %.4f..%.4f cm", det_name.c_str(),
               crmin / cm, crmax / cm, (zc + dz + cz - chalfz) / cm, (zc + dz + cz + chalfz) / cm);
    }

    Volume enc_vol(det_name + "_Enclosure", enc_s, description.material(x_enc.materialStr()));
    enc_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_enc, _Unicode(vis), "Target_EnclosureVis"));
    Volume vac_vol(det_name + "_EnclosureVacuum", vac_s, description.material(fmat));
    vac_vol.setVisAttributes(description, "InvisibleWithDaughters");
    enc_vol.placeVolume(vac_vol);
    vac_vol.placeVolume(cell_vol, Position(0., 0., -dz));

    for (xml_coll_t io(x_enc, _Unicode(opening)); io; ++io) {
      xml_comp_t x_o         = io;
      const std::string side = x_o.attr<std::string>(_Unicode(side));
      const double r         = x_o.attr<double>(_Unicode(r));
      if (r > rmax - ewall) {
        throw std::runtime_error("SoLID_CryoTarget: opening radius larger than the enclosure vacuum");
      }
      const double sgn = (side == "upstream") ? -1. : 1.;
      Volume ov(det_name + "_Enclosure_" + side + "_opening", Tube(0., r, 0.5 * endwall), description.material(fmat));
      ov.setVisAttributes(description, "InvisibleNoDaughters");
      enc_vol.placeVolume(ov, Position(0., 0., sgn * (halfz - 0.5 * endwall)));
    }
    assembly.placeVolume(enc_vol, Position(0., 0., zc + dz));
    printout(INFO, "SoLID_CryoTarget", "%s enclosure: %s r %.3f cm, hall z %.3f..%.3f cm, %d cavity(ies)",
             det_name.c_str(), x_enc.materialStr().c_str(), rmax / cm, (zc + dz - halfz) / cm, (zc + dz + halfz) / cm,
             ncav);
  } else if (x_det.hasChild(_Unicode(envelope))) {
    // -------------------------------------------------------------- simple envelope (GEMC TACH, air)
    xml_comp_t x_env   = x_det.child(_Unicode(envelope));
    const double rmax  = x_env.attr<double>(_Unicode(rmax));
    const double halfz = getAttrOrDefault<double>(x_env, _Unicode(halfz), cell_halfz);
    if (rmax < cell_rmax || halfz < cell_halfz) {
      throw std::runtime_error("SoLID_CryoTarget: the envelope is smaller than the cell");
    }
    Volume env_vol(det_name + "_Envelope", Tube(0., rmax, halfz), description.material(x_env.materialStr()));
    env_vol.setVisAttributes(description, "InvisibleWithDaughters");
    env_vol.placeVolume(cell_vol);
    assembly.placeVolume(env_vol, Position(0., 0., zc));
  } else {
    assembly.placeVolume(cell_vol, Position(0., 0., zc));
  }

  PlacedVolume pv = description.pickMotherVolume(det).placeVolume(assembly);
  det.setPlacement(pv);
  return det;
}

DECLARE_DETELEMENT(SoLID_CryoTarget, create_cryotarget)
