// SPDX-License-Identifier: LGPL-3.0-or-later
//
// SoLID_MuonDetector: DDVCS muon system components (GEMC muon_moved, 2026-10).
//
// One plugin, used by several <detector> entries of compact/muon/muon_DDVCS.xml:
//   * MuonIron     - absorber layers made of boxes (CLEO-II layer-3 slabs), passive;
//   * MuonScint    - scintillator planes, segmented in N_phi azimuthal paddles, sensitive;
//   * MuonObserver - thin vacuum "observer" (virtual / flux) planes: disks, cylinders, polygon barrels, sensitive.
// Each entry may contain any of the elements below; all numbers come from the compact file.
//
//   <layers count="3" z0="638*cm" pitch="46*cm"/>       repeated layers, layer L (1..count) centred at
//                                                        z_L = z0 + (L-1)*pitch (default: one layer at z0 = 0)
//   <block name="..." x=".." y=".." z="0" dx=".." dy=".." dz=".." material=".." vis=".."/>
//                                                        Box (half lengths) at (x, y, z_L + z) in every layer
//   <disk z="0" rmin=".." rmax=".." thickness=".." nphi="60" phi0="90*deg" material=".." vis=".."
//         sensitive="true"/>                             annulus at z_L + z, cut in nphi paddles (Tube segments);
//                                                        if sensitive: volIDs layer, module (paddle 1..nphi),
//                                                        paddle m covers [phi0 + (m-1)*dphi, phi0 + m*dphi]
//   <plane id="0" enable="1" type="disk" z=".." rmin=".." rmax=".." thickness=".." material=".." vis=".."/>
//   <plane id="10" enable="1" type="cylinder" r=".." thickness=".." zmin=".." zmax=".." .../>
//   <plane id="11" enable="1" type="polygon" nsides="8" phi0="-22.5*deg" apothem=".." thickness=".."
//          zmin=".." zmax=".." .../>                     single planes (not repeated per layer), volID plane = id;
//                                                        polygon radii are apothems (flat distances), phi0 is the
//                                                        Polyhedra start angle (first corner)
// Sensitive parts use the detector readout (type tracker).

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <cctype>
#include <cmath>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

// ---- compact parameters with "first definition wins" semantics ----------------------------------------------
// DD4hep keeps the FIRST definition of a duplicated <constant> in Detector::constants(), but its expression evaluator
// keeps the LAST one. A switch defined in the top-level <define> (before this fragment is included) would therefore
// be replaced by the fragment default when attributes are evaluated the usual way. Attributes read with param()
// resolve every constant name through Detector::constants() (recursively), so the first definition always wins.
std::string expand_constants(Detector& desc, const std::string& expr, int depth) {
  if (depth > 32) {
    throw std::runtime_error("constant expansion too deep (recursive definition?): " + expr);
  }
  const auto& consts = desc.constants();
  std::string out;
  size_t i = 0;
  while (i < expr.size()) {
    const char c = expr[i];
    if (std::isdigit(static_cast<unsigned char>(c)) ||
        (c == '.' && i + 1 < expr.size() && std::isdigit(static_cast<unsigned char>(expr[i + 1])))) {
      size_t j = i;
      while (j < expr.size() && (std::isdigit(static_cast<unsigned char>(expr[j])) || expr[j] == '.')) {
        ++j;
      }
      if (j < expr.size() && (expr[j] == 'e' || expr[j] == 'E')) {
        size_t k = j + 1;
        if (k < expr.size() && (expr[k] == '+' || expr[k] == '-')) {
          ++k;
        }
        if (k < expr.size() && std::isdigit(static_cast<unsigned char>(expr[k]))) {
          j = k;
          while (j < expr.size() && std::isdigit(static_cast<unsigned char>(expr[j]))) {
            ++j;
          }
        }
      }
      out += expr.substr(i, j - i);
      i = j;
    } else if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
      size_t j = i;
      while (j < expr.size() && (std::isalnum(static_cast<unsigned char>(expr[j])) || expr[j] == '_')) {
        ++j;
      }
      const std::string tok = expr.substr(i, j - i);
      auto it = consts.find(tok);
      if (it != consts.end() && it->second.isValid()) {
        out += "(" + expand_constants(desc, it->second->GetTitle(), depth + 1) + ")";
      } else {
        out += tok;
      }
      i = j;
    } else {
      out += c;
      ++i;
    }
  }
  return out;
}

[[maybe_unused]] double param(Detector& desc, xml::Component x, const xml::XmlChar* attr) {
  return _toDouble(expand_constants(desc, x.attr<std::string>(attr), 0));
}

[[maybe_unused]] double param(Detector& desc, xml::Component x, const xml::XmlChar* attr, double def) {
  return x.hasAttr(attr) ? param(desc, x, attr) : def;
}

[[maybe_unused]] int iparam(Detector& desc, xml::Component x, const xml::XmlChar* attr) {
  return static_cast<int>(std::lround(param(desc, x, attr)));
}

[[maybe_unused]] int iparam(Detector& desc, xml::Component x, const xml::XmlChar* attr, int def) {
  return x.hasAttr(attr) ? iparam(desc, x, attr) : def;
}

} // namespace

static Ref_t create_muon_detector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml::DetElement x_det = e;
  const std::string det_name = x_det.nameStr();
  const int det_id = x_det.id();

  int n_layers = 1;
  double z0 = 0., pitch = 0.;
  if (x_det.hasChild(_Unicode(layers))) {
    xml::Component x_lay = x_det.child(_Unicode(layers));
    n_layers = iparam(desc, x_lay, _Unicode(count));
    z0 = param(desc, x_lay, _Unicode(z0));
    pitch = param(desc, x_lay, _Unicode(pitch), 0.);
  }

  const bool has_readout = sens.isValid() && sens.readout().isValid();
  if (has_readout) {
    sens.setType("tracker");
  }

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  // ---------------- shared volumes for the repeated layers
  struct BlockDef {
    Volume vol;
    Position pos;
  };
  std::vector<BlockDef> blocks;
  int ib = 0;
  for (xml_coll_t xb(x_det, _Unicode(block)); xb; ++xb, ++ib) {
    xml::Component x_b = xb;
    const std::string bname = x_b.hasAttr(_U(name)) ? x_b.nameStr() : "block" + std::to_string(ib + 1);
    Box shape(param(desc, x_b, _U(dx)), param(desc, x_b, _U(dy)), param(desc, x_b, _U(dz)));
    Volume vol(det_name + "_" + bname, shape, desc.material(x_b.attr<std::string>(_U(material))));
    vol.setVisAttributes(desc.visAttributes(x_b.visStr()));
    const double bz = param(desc, x_b, _U(z), 0.);
    blocks.push_back({vol, Position(param(desc, x_b, _U(x)), param(desc, x_b, _U(y)), bz)});
  }

  struct DiskDef {
    Volume pad;
    double z;
    int nphi;
    double phi0;
    bool sensitive;
  };
  std::vector<DiskDef> disks;
  int id = 0;
  for (xml_coll_t xd(x_det, _Unicode(disk)); xd; ++xd, ++id) {
    xml::Component x_d = xd;
    const int nphi = iparam(desc, x_d, _Unicode(nphi), 1);
    const double phi0 = param(desc, x_d, _Unicode(phi0), 0.);
    const double dphi = 2. * M_PI / nphi;
    const double thick = param(desc, x_d, _U(thickness));
    const bool sensitive = x_d.hasAttr(_U(sensitive)) && x_d.attr<bool>(_U(sensitive));
    Solid shape = (nphi == 1)
                      ? Solid(Tube(param(desc, x_d, _U(rmin)), param(desc, x_d, _U(rmax)), 0.5 * thick))
                      : Solid(Tube(param(desc, x_d, _U(rmin)), param(desc, x_d, _U(rmax)), 0.5 * thick, -0.5 * dphi, 0.5 * dphi));
    Volume pad(det_name + "_disk" + std::to_string(id + 1) + "_paddle", shape,
               desc.material(x_d.attr<std::string>(_U(material))));
    pad.setVisAttributes(desc.visAttributes(x_d.visStr()));
    if (sensitive) {
      if (!has_readout) {
        throw std::runtime_error(det_name + ": sensitive <disk> but no readout");
      }
      pad.setSensitiveDetector(sens);
    }
    disks.push_back({pad, param(desc, x_d, _U(z), 0.), nphi, phi0, sensitive});
  }

  // ---------------- layers
  int n_blocks = 0, n_paddles = 0;
  if (!blocks.empty() || !disks.empty()) {
    for (int L = 1; L <= n_layers; ++L) {
      const double zL = z0 + (L - 1) * pitch;
      Assembly lay(det_name + "_layer" + std::to_string(L));
      DetElement de_lay(sdet, "layer" + std::to_string(L), L);
      for (const auto& b : blocks) {
        lay.placeVolume(b.vol, b.pos);
        ++n_blocks;
      }
      for (const auto& d : disks) {
        const double dphi = 2. * M_PI / d.nphi;
        for (int ip = 1; ip <= d.nphi; ++ip) {
          const double phi_c = (d.nphi == 1) ? 0. : d.phi0 + (ip - 0.5) * dphi;
          PlacedVolume pv = lay.placeVolume(d.pad, Transform3D(RotationZ(phi_c), Position(0., 0., d.z)));
          if (d.sensitive) {
            pv.addPhysVolID("module", ip);
            DetElement de_pad(de_lay, "paddle" + std::to_string(ip), ip);
            de_pad.setPlacement(pv);
          }
          ++n_paddles;
        }
      }
      lay.setVisAttributes(desc.invisible());
      PlacedVolume pv_lay = assembly.placeVolume(lay, Position(0., 0., zL));
      if (has_readout) {
        pv_lay.addPhysVolID("layer", L);
      }
      de_lay.setPlacement(pv_lay);
    }
  }

  // ---------------- single (observer) planes
  int n_planes = 0;
  for (xml_coll_t xp(x_det, _Unicode(plane)); xp; ++xp) {
    xml::Component x_p = xp;
    if (iparam(desc, x_p, _Unicode(enable), 1) == 0) {
      continue;
    }
    const int pid = x_p.id();
    const std::string type = x_p.attr<std::string>(_U(type));
    const double thick = param(desc, x_p, _U(thickness));
    Solid shape;
    Position pos;
    if (type == "disk") {
      shape = Tube(param(desc, x_p, _U(rmin)), param(desc, x_p, _U(rmax)), 0.5 * thick);
      pos = Position(0., 0., param(desc, x_p, _U(z)));
    } else if (type == "cylinder") {
      const double zmin = param(desc, x_p, _U(zmin)), zmax = param(desc, x_p, _U(zmax));
      const double r = param(desc, x_p, _U(r));
      shape = Tube(r, r + thick, 0.5 * (zmax - zmin));
      pos = Position(0., 0., 0.5 * (zmin + zmax));
    } else if (type == "polygon") {
      const double zmin = param(desc, x_p, _U(zmin)), zmax = param(desc, x_p, _U(zmax));
      const double a = param(desc, x_p, _Unicode(apothem));
      const int nsides = iparam(desc, x_p, _Unicode(nsides));
      const double phi0 = param(desc, x_p, _Unicode(phi0));
      shape = Polyhedra(nsides, phi0, 2. * M_PI, {zmin, zmax}, {a, a}, {a + thick, a + thick});
      pos = Position(0., 0., 0.);
    } else {
      throw std::runtime_error(det_name + ": unknown plane type '" + type + "'");
    }
    Volume vol(det_name + "_plane" + std::to_string(pid), shape, desc.material(x_p.attr<std::string>(_U(material))));
    vol.setVisAttributes(desc.visAttributes(x_p.visStr()));
    PlacedVolume pv = assembly.placeVolume(vol, pos);
    if (has_readout) {
      vol.setSensitiveDetector(sens);
      pv.addPhysVolID("plane", pid);
    }
    DetElement de_pl(sdet, "plane" + std::to_string(pid), pid);
    de_pl.setPlacement(pv);
    ++n_planes;
  }

  printout(INFO, det_name, "%d layer(s) at z0 = %.2f cm, pitch %.2f cm: %d blocks, %d paddles; %d single planes",
           n_layers, z0 / cm, pitch / cm, n_blocks, n_paddles, n_planes);

  assembly.setVisAttributes(desc.visAttributes(x_det.visStr()));
  PlacedVolume pv = desc.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(SoLID_MuonDetector, create_muon_detector)
