// SPDX-License-Identifier: LGPL-3.0-or-later
//
// SoLID_MRPCDisk: forward-angle multi-gap resistive plate chamber disk (FAMRPC).
//
// The disk is an annulus (rmin, rmax) centred at z, filled by a list of coaxial layers ("slices") given in the
// compact file from upstream to downstream. Several alternative stacks may be listed (<stack variant="k">); the one
// whose variant equals <dimensions stack="..."/> is built. Slices with sensitive="true" are gas gaps: they carry the
// volume IDs "stack" (attribute of the slice) and "gap" (counted 1..n inside each stack). The strip readout is
// a PolarGridRPhi segmentation of the readout, evaluated in the gas-gap frame (= hall frame orientation, the
// layers are never rotated).
//
// Optional azimuthal mask (SoLID-Lite forward-angle segments): the disk is then built as one phi segment per window
// (module = window index 1..N); without mask a single full annulus is built (module = 1). The segment shapes are
// built at their global phi (no rotation), so that the PolarGridRPhi phi index is the hall-frame phi.
//
// Compact:
//   <detector id="MRPCFA_ID" name="FAMRPC" type="SoLID_MRPCDisk" readout="FAMRPCHits" vis="...">
//     <dimensions z="..." rmin="..." rmax="..." stack="0"/>
//     <envelope material="Air" vis="..."/>
//     <stack variant="0" name="preCDR">
//       <slice material="..." thickness="..." vis="..." [sensitive="true" stack="1"]/> ...
//     </stack>
//     <stack variant="1" name="gemc"> ... </stack>
//     <phi_mask enable="0|1" count="N" width="..."> <window center="..."/> ... </phi_mask>
//   </detector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <cctype>
#include <cmath>
#include <map>
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

std::vector<PhiWindow> read_phi_mask(Detector& desc, xml::Component x_det, const std::string& det_name) {
  std::vector<PhiWindow> wins;
  if (!x_det.hasChild(_Unicode(phi_mask))) {
    return wins;
  }
  xml::Component x_mask = x_det.child(_Unicode(phi_mask));
  if (iparam(desc, x_mask, _Unicode(enable)) == 0) {
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

struct Slice {
  Material mat;
  double thick;
  std::string vis;
  bool sensitive;
  int stack;
  std::string name;
};

} // namespace

static Ref_t create_mrpc_disk(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml::DetElement x_det = e;
  const std::string det_name = x_det.nameStr();
  const int det_id = x_det.id();

  xml::Component x_dim = x_det.child(_U(dimensions));
  const double z_c = param(desc, x_dim, _Unicode(z));
  const double rmin = param(desc, x_dim, _U(rmin));
  const double rmax = param(desc, x_dim, _U(rmax));
  const int variant = iparam(desc, x_dim, _Unicode(stack));

  xml::Component x_env = x_det.child(_Unicode(envelope));
  Material env_mat = desc.material(x_env.attr<std::string>(_U(material)));
  const std::string env_vis = x_env.visStr();

  // ---- select the layer stack
  std::vector<Slice> slices;
  std::string stack_name;
  for (xml_coll_t xs(x_det, _Unicode(stack)); xs; ++xs) {
    xml::Component x_stack = xs;
    if (iparam(desc, x_stack, _Unicode(variant)) != variant) {
      continue;
    }
    stack_name = x_stack.hasAttr(_U(name)) ? x_stack.nameStr() : std::to_string(variant);
    int i = 0;
    for (xml_coll_t xl(x_stack, _U(slice)); xl; ++xl, ++i) {
      xml::Component x_sl = xl;
      Slice sl;
      sl.mat = desc.material(x_sl.attr<std::string>(_U(material)));
      sl.thick = param(desc, x_sl, _U(thickness));
      sl.vis = x_sl.hasAttr(_U(vis)) ? x_sl.visStr() : env_vis;
      sl.sensitive = x_sl.hasAttr(_U(sensitive)) && x_sl.attr<bool>(_U(sensitive));
      sl.stack = iparam(desc, x_sl, _Unicode(stack), 0);
      sl.name = "slice" + std::to_string(i + 1);
      slices.push_back(sl);
    }
  }
  if (slices.empty()) {
    throw std::runtime_error(det_name + ": no <stack> with variant " + std::to_string(variant));
  }
  double total = 0.;
  for (const auto& sl : slices) {
    total += sl.thick;
  }

  const auto mask = read_phi_mask(desc, x_det, det_name);
  // azimuthal pieces: (start, delta); the full annulus when no mask
  std::vector<std::pair<double, double>> pieces;
  if (mask.empty()) {
    pieces.emplace_back(0., 2. * M_PI);
  } else {
    for (const auto& w : mask) {
      pieces.emplace_back(w.center - 0.5 * w.width, w.width);
    }
  }

  sens.setType("tracker");

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  int n_gaps = 0;
  for (size_t im = 0; im < pieces.size(); ++im) {
    const int module = im + 1;
    const double phi_s = pieces[im].first;
    const double phi_e = pieces[im].first + pieces[im].second;
    const bool full = mask.empty();
    const std::string mname = det_name + "_module" + std::to_string(module);

    auto make_tube = [&](double dz) -> Solid {
      if (full) {
        return Tube(rmin, rmax, dz);
      }
      return Tube(rmin, rmax, dz, phi_s, phi_e);
    };

    Volume env_vol(mname, make_tube(0.5 * total), env_mat);
    env_vol.setVisAttributes(desc.visAttributes(env_vis));
    DetElement de_mod(sdet, "module" + std::to_string(module), module);

    std::map<int, int> gap_count;
    double zpos = -0.5 * total;
    for (const auto& sl : slices) {
      Volume sl_vol(mname + "_" + sl.name, make_tube(0.5 * sl.thick), sl.mat);
      sl_vol.setVisAttributes(desc.visAttributes(sl.vis));
      PlacedVolume pv = env_vol.placeVolume(sl_vol, Position(0., 0., zpos + 0.5 * sl.thick));
      if (sl.sensitive) {
        const int gap = ++gap_count[sl.stack];
        sl_vol.setSensitiveDetector(sens);
        pv.addPhysVolID("stack", sl.stack).addPhysVolID("gap", gap);
        DetElement de_gap(de_mod, "stack" + std::to_string(sl.stack) + "_gap" + std::to_string(gap),
                          100 * sl.stack + gap);
        de_gap.setPlacement(pv);
        ++n_gaps;
      }
      zpos += sl.thick;
    }

    PlacedVolume pv_mod = assembly.placeVolume(env_vol, Position(0., 0., z_c));
    pv_mod.addPhysVolID("module", module);
    de_mod.setPlacement(pv_mod);
  }

  printout(INFO, det_name,
           "stack '%s': %zu slices, %.3f mm total, %d sensitive gaps in %zu module(s) (mask %s), z = %.3f cm, "
           "r = %.1f-%.1f cm",
           stack_name.c_str(), slices.size(), total / mm, n_gaps, pieces.size(), mask.empty() ? "off" : "on",
           z_c / cm, rmin / cm, rmax / cm);

  assembly.setVisAttributes(desc.visAttributes(x_det.visStr()));
  PlacedVolume pv = desc.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(SoLID_MRPCDisk, create_mrpc_disk)
