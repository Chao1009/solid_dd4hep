// SPDX-License-Identifier: LGPL-3.0-or-later
//
// SoLID_ScintPadDisk: scintillator pad detector (FASPD / LASPD).
//
// A disk perpendicular to the beam, made of fan-shaped pads: N_phi azimuthal sectors x N_ring radial rings.
// Every pad is a Tube segment (rmin, rmax, thickness, dphi) and is sensitive; it gets its own DetElement and the
// volume IDs "sector" (1..N_phi) and, if the readout has that field, "ring" (1..N_ring).
//
// Sector convention (GEMC analysis/analysis_tree_solid_spd.C): sector s covers [phi0 + (s-1)*dphi, phi0 + s*dphi],
// with phi0 = 90 deg by default (sector 1 = 90..96 deg for 60 sectors).
//
// Optional azimuthal mask (SoLID-Lite forward-angle segments): only pads whose centre lies inside one of the
// windows are built (the sector numbering stays the global one).
//
// Compact:
//   <detector id="SPDFA_ID" name="FASPD" type="SoLID_ScintPadDisk" readout="FASPDHits" vis="...">
//     <dimensions z="..." thickness="..."/>                          z = centre of the pads (hall frame)
//     <pads nphi="60" phi0="90*deg" material="PlasticScint" vis="..."/>
//     <ring id="1" rmin="..." rmax="..."/> ...                         radial rings (at least one)
//     <phi_mask enable="0|1" count="N" width="24*deg">                 optional
//       <window center="42*deg"/> ...                                  the first `count` windows are used
//     </phi_mask>
//   </detector>

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

struct PhiWindow {
  double center;
  double width;
};

// read <phi_mask>; returns an empty vector when the mask is absent or disabled
std::vector<PhiWindow> read_phi_mask(Detector& desc, xml::Component x_det, const std::string& det_name) {
  std::vector<PhiWindow> wins;
  if (!x_det.hasChild(_Unicode(phi_mask))) {
    return wins;
  }
  xml::Component x_mask = x_det.child(_Unicode(phi_mask));
  const int enable = iparam(desc, x_mask, _Unicode(enable));
  if (enable == 0) {
    return wins;
  }
  const int count = iparam(desc, x_mask, _Unicode(count));
  const double width = param(desc, x_mask, _Unicode(width));
  for (xml_coll_t xw(x_mask, _Unicode(window)); xw && (int)wins.size() < count; ++xw) {
    xml::Component x_win = xw;
    wins.push_back({param(desc, x_win, _Unicode(center)), width});
  }
  if ((int)wins.size() != count) {
    printout(ERROR, det_name, "phi_mask: count = %d but only %zu <window> elements", count, wins.size());
    throw std::runtime_error(det_name + ": inconsistent phi_mask");
  }
  return wins;
}

bool inside_mask(const std::vector<PhiWindow>& wins, double phi) {
  if (wins.empty()) {
    return true;
  }
  for (const auto& w : wins) {
    const double d = std::remainder(phi - w.center, 2. * M_PI);
    if (std::fabs(d) < 0.5 * w.width) {
      return true;
    }
  }
  return false;
}

bool readout_has_field(SensitiveDetector sens, const std::string& field) {
  Readout ro = sens.readout();
  if (!ro.isValid()) {
    return false;
  }
  for (const auto& f : ro.idSpec().fields()) {
    if (f.first == field) {
      return true;
    }
  }
  return false;
}

} // namespace

static Ref_t create_scint_pad_disk(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml::DetElement x_det = e;
  const std::string det_name = x_det.nameStr();
  const int det_id = x_det.id();

  xml::Component x_dim = x_det.child(_U(dimensions));
  const double z_c = param(desc, x_dim, _Unicode(z));
  const double thick = param(desc, x_dim, _U(thickness));

  xml::Component x_pads = x_det.child(_Unicode(pads));
  const int nphi = iparam(desc, x_pads, _Unicode(nphi));
  const double phi0 = param(desc, x_pads, _Unicode(phi0));
  const double dphi = 2. * M_PI / nphi;
  Material mat = desc.material(x_pads.attr<std::string>(_U(material)));
  const std::string pad_vis = x_pads.hasAttr(_U(vis)) ? x_pads.visStr() : x_det.visStr();

  struct Ring {
    int id;
    double rmin, rmax;
  };
  std::vector<Ring> rings;
  for (xml_coll_t xr(x_det, _U(ring)); xr; ++xr) {
    xml::Component x_ring = xr;
    rings.push_back({x_ring.id(), param(desc, x_ring, _U(rmin)), param(desc, x_ring, _U(rmax))});
  }
  if (rings.empty()) {
    throw std::runtime_error(det_name + ": no <ring> defined");
  }

  const auto mask = read_phi_mask(desc, x_det, det_name);
  const bool use_ring_id = readout_has_field(sens, "ring");
  if (!use_ring_id && rings.size() > 1) {
    throw std::runtime_error(det_name + ": more than one ring but the readout has no 'ring' field");
  }

  sens.setType("tracker");

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  int n_pads = 0;
  for (const auto& ring : rings) {
    // one pad shape per ring, centred on phi = 0 and rotated into place
    Tube pad_shape(ring.rmin, ring.rmax, 0.5 * thick, -0.5 * dphi, 0.5 * dphi);
    Volume pad_vol(det_name + "_pad_ring" + std::to_string(ring.id), pad_shape, mat);
    pad_vol.setVisAttributes(desc.visAttributes(pad_vis));
    pad_vol.setSensitiveDetector(sens);

    for (int isec = 1; isec <= nphi; ++isec) {
      const double phi_c = phi0 + (isec - 0.5) * dphi;
      if (!inside_mask(mask, phi_c)) {
        continue;
      }
      PlacedVolume pv = assembly.placeVolume(pad_vol, Transform3D(RotationZ(phi_c), Position(0., 0., z_c)));
      pv.addPhysVolID("sector", isec);
      if (use_ring_id) {
        pv.addPhysVolID("ring", ring.id);
      }
      DetElement de_pad(sdet, "pad_s" + std::to_string(isec) + "_r" + std::to_string(ring.id),
                        1000 * ring.id + isec);
      de_pad.setPlacement(pv);
      ++n_pads;
    }
  }

  printout(INFO, det_name, "%d pads (%d sectors x %zu rings, mask %s), z = %.3f cm, thickness %.3f cm", n_pads,
           nphi, rings.size(), mask.empty() ? "off" : "on", z_c / cm, thick / cm);

  assembly.setVisAttributes(desc.visAttributes(x_det.visStr()));
  PlacedVolume pv = desc.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(SoLID_ScintPadDisk, create_scint_pad_disk)
