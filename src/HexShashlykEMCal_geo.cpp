// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 SoLID collaboration
//
// SoLID_HexShashlykEMCal: SoLID electromagnetic calorimeter built from hexagonal modules on a fixed lattice.
//
// One <detector> builds ONE readout unit: the shower (194-layer shashlik modules) or the preshower (2 cm
// scintillator tile + 2 cm Al support tile per module, behind a full Pb radiator sheet). FAEC/LAEC and the
// SIDIS / PVDIS / lite variants differ only in compact constants (see compact/ecal/*.xml).
//
// Compact structure (all numbers from XML, units explicit):
//   <detector id= name= type="SoLID_HexShashlykEMCal" readout= vis=(envelope vis)>
//     <lattice type="anl_hex_rows" rows="k:a:b[,a:b];..." | string-constant name
//              pitch= row_pitch= row_center= sector_dphi=/>
//     <selection rmin= rmax=/>                              strict rmin < R(centre) < rmax (GEMC make_ec_module)
//     <phi_mask enable= overlap= nseg= dphi= phi0= centers="c1,c2,..." | string-constant name/>   (optional)
//     <envelope type="tube|assembly" margin= material=/>    (optional, default tube, 0.1 mm, Air)
//     <module type="hex_stack" nsides= apothem= phi0= z_front= material= vis=>
//       <layer repeat=>  <slice name= material= thickness= [sensitive="true"] vis=>
//                           [<rods enable= n= radius= r= phi0= material= vis=/>]  </slice> ... </layer>
//     </module>
//     <sheet name= material= rmin= rmax= z_front= thickness= vis=/>    (optional, any number; e.g. PS Pb radiator)
//   </detector>
//
// Volume / id structure (keeps the build fast for 1400 modules x 194 layers x 6 slices):
//   envelope (system)                     DetElement <name>
//     module volume (one logical volume shared by all modules)   volIDs sector, module (ANL id)
//                                         DetElement module<id> per module (no DetElements below)
//       slice placements directly in the module (no layer envelopes)  volIDs layer (1..N), slice (1..6)
//         rods are daughters of the slice logical volume (defined once, reused by every layer/module)
//   sheets are placed in the envelope (passive, no volIDs needed)
//
// Pluggable module type: <module type=...> selects a module builder from `module_builders()`. A future
// sPHENIX-block / STAR-FCS option adds a builder (returning its logical volume, length and footprint polygon)
// and, for a rectangular layout, a <lattice type="rect_grid"> branch in `make_sites()`; selection, phi mask,
// sector rule, envelope, DetElements and readout handling stay common.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/detail/ObjectsInterna.h"
#include "XML/Utilities.h"

#include "SoLIDHexLattice.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace dd4hep;
namespace ec = solid::ec;

namespace {

/// value of a string attribute, resolved through a string <constant> of that name if one exists
std::string string_attr(Detector& desc, xml_comp_t x, const char* name) {
  std::string s = x.attr<std::string>(xml::Strng_t(name));
  try {
    return desc.constantAsString(s);
  } catch (...) {
    return s;
  }
}

/// Evaluate an expression with "first definition wins" semantics for <constant>s.
/// DD4hep (1.38) keeps the FIRST definition of a duplicated constant in the Detector constant table
/// ("already defined. New value will be ignored") but OVERWRITES the expression-evaluator variable with the
/// LAST one, so a plain attr<double>() would silently use the fragment default instead of a value pre-defined
/// by the top-level file. Here every identifier that names a numeric constant is replaced (recursively) by the
/// value of its stored (first) definition, so documented constants can be overridden by defining them before
/// the fragment is included - including the constants derived from them inside the fragment.
double eval_first(Detector& desc, const std::string& expr, int depth = 0) {
  if (depth > 32) throw std::runtime_error("SoLID_HexShashlykEMCal: constant recursion too deep in '" + expr + "'");
  const auto& consts = desc.constants();
  std::string out;
  size_t i = 0;
  while (i < expr.size()) {
    const char c = expr[i];
    const bool id_start = std::isalpha(static_cast<unsigned char>(c)) || c == '_';
    const bool in_number = i > 0 && (std::isdigit(static_cast<unsigned char>(expr[i - 1])) || expr[i - 1] == '.');
    if (id_start && !in_number) {
      size_t j = i;
      while (j < expr.size() && (std::isalnum(static_cast<unsigned char>(expr[j])) || expr[j] == '_')) ++j;
      const std::string ident = expr.substr(i, j - i);
      auto it = consts.find(ident);
      Constant cst = (it != consts.end()) ? Constant(it->second) : Constant();
      if (cst.isValid() && cst->dataType != "string") {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "(%.17g)", eval_first(desc, cst->GetTitle(), depth + 1));
        out += buf;
      } else {
        out += ident;
      }
      i = j;
    } else {
      out += c;
      ++i;
    }
  }
  return _toDouble(out);
}

double num(Detector& desc, xml_comp_t x, const char* name) {
  return eval_first(desc, x.attr<std::string>(xml::Strng_t(name)));
}
double num_or(Detector& desc, xml_comp_t x, const char* name, double def) {
  return x.hasAttr(xml::Strng_t(name)) ? num(desc, x, name) : def;
}
std::string str_or(xml_comp_t x, const char* name, const std::string& def) {
  return x.hasAttr(xml::Strng_t(name)) ? x.attr<std::string>(xml::Strng_t(name)) : def;
}
int inum(double v) { return static_cast<int>(std::lround(v)); }

std::vector<double> parse_list(Detector& desc, const std::string& s) {
  std::vector<double> out;
  std::string tok;
  std::string t = s;
  for (auto& c : t)
    if (c == ',' || c == ';') c = ' ';
  std::stringstream ss(t);
  while (ss >> tok) out.push_back(eval_first(desc, tok));
  return out;
}

/// what a module builder returns
struct ModuleType {
  Volume volume;
  double length = 0;   ///< full length along z
  int nsides = 6;      ///< footprint polygon (for extents and the overlap-mode phi mask)
  double apothem = 0;
  double phi0 = 0;
};

/// hex_stack: regular-polygon prism (Air) filled by <layer repeat><slice/></layer> slices, slices placed
/// directly in the module with volIDs layer/slice; optional rods as daughters of slice logical volumes.
ModuleType build_hex_stack(Detector& desc, xml_comp_t x_mod, SensitiveDetector& sens, const std::string& prefix) {
  ModuleType mt;
  mt.nsides  = inum(num_or(desc, x_mod, "nsides", 6));
  mt.apothem = num(desc, x_mod, "apothem");
  mt.phi0    = num_or(desc, x_mod, "phi0", 0.);

  struct SliceDef {
    Volume vol;
    double thick;
    int num;
  };
  struct LayerDef {
    std::vector<SliceDef> slices;
    double thick = 0;
    int repeat = 1;
  };
  std::vector<LayerDef> layers;

  for (xml_coll_t li(x_mod, _U(layer)); li; ++li) {
    xml_comp_t x_layer = li;
    LayerDef ld;
    ld.repeat = inum(num(desc, x_layer, "repeat"));
    int s_num = 0;
    for (xml_coll_t si(x_layer, _U(slice)); si; ++si) {
      xml_comp_t x_slice = si;
      ++s_num;
      const double t = num(desc, x_slice, "thickness");
      if (t <= 0) continue;  // e.g. paint switched off: keep the slice numbering stable
      std::string sname = x_slice.hasAttr(_U(name)) ? x_slice.nameStr() : _toString(s_num, "slice%d");
      PolyhedraRegular shape(mt.nsides, mt.phi0, 0., mt.apothem, t);
      Volume vol(prefix + "_" + sname, shape, desc.material(x_slice.materialStr()));
      if (x_slice.isSensitive()) {
        vol.setSensitiveDetector(sens);
      }
      vol.setVisAttributes(desc, x_slice.visStr());
      // optional rods (e.g. 6 stainless-steel tie rods through Pb and scintillator)
      for (xml_coll_t ri(x_slice, _Unicode(rods)); ri; ++ri) {
        xml_comp_t x_rod = ri;
        if (!inum(num_or(desc, x_rod, "enable", 1))) continue;
        const int n       = inum(num_or(desc, x_rod, "n", 6));
        const double rr   = num(desc, x_rod, "radius");
        const double rpos = num(desc, x_rod, "r");
        const double p0   = num_or(desc, x_rod, "phi0", 0.);
        Volume rvol(prefix + "_" + sname + "_rod", Tube(0., rr, t / 2.), desc.material(x_rod.materialStr()));
        rvol.setVisAttributes(desc, x_rod.visStr());
        for (int i = 0; i < n; ++i) {
          double p = p0 + 2. * M_PI * i / n;
          vol.placeVolume(rvol, Position(rpos * std::cos(p), rpos * std::sin(p), 0.));
        }
      }
      ld.slices.push_back({vol, t, s_num});
      ld.thick += t;
    }
    layers.push_back(ld);
  }

  for (auto& ld : layers) mt.length += ld.thick * ld.repeat;
  if (mt.length <= 0) throw std::runtime_error(prefix + ": hex_stack module without layers");

  PolyhedraRegular mshape(mt.nsides, mt.phi0, 0., mt.apothem, mt.length);
  Volume mvol(prefix + "_module", mshape, desc.material(str_or(x_mod, "material", "Air")));
  mvol.setVisAttributes(desc, x_mod.visStr());

  // place slices; positions accumulated from the front face (upstream, -z)
  double z = -mt.length / 2.;
  int l_num = 0;
  for (auto& ld : layers) {
    for (int j = 0; j < ld.repeat; ++j) {
      ++l_num;
      for (auto& sd : ld.slices) {
        PlacedVolume pv = mvol.placeVolume(sd.vol, Position(0., 0., z + sd.thick / 2.));
        pv.addPhysVolID("layer", l_num).addPhysVolID("slice", sd.num);
        z += sd.thick;
      }
    }
  }
  mt.volume = mvol;
  return mt;
}

using ModuleBuilder = std::function<ModuleType(Detector&, xml_comp_t, SensitiveDetector&, const std::string&)>;

const std::map<std::string, ModuleBuilder>& module_builders() {
  static const std::map<std::string, ModuleBuilder> m{
      {"hex_stack", build_hex_stack},
      // {"spacal_block", build_spacal_block},   // future: sPHENIX W/SciFi blocks (needs a rect lattice)
  };
  return m;
}

std::vector<ec::LatticeSite> make_sites(Detector& desc, xml_comp_t x_lat) {
  const std::string type = str_or(x_lat, "type", "anl_hex_rows");
  if (type != "anl_hex_rows")
    throw std::runtime_error("SoLID_HexShashlykEMCal: unknown lattice type '" + type + "'");
  auto rows = ec::parse_rows(string_attr(desc, x_lat, "rows"));
  return ec::build_sites(rows, num(desc, x_lat, "pitch"), num(desc, x_lat, "row_pitch"), num(desc, x_lat, "row_center"),
                         num(desc, x_lat, "sector_dphi"));
}

ec::PhiMask make_phi_mask(Detector& desc, xml_det_t x_det) {
  ec::PhiMask mask;
  if (!x_det.hasChild(_Unicode(phi_mask))) return mask;
  xml_comp_t x_m = x_det.child(_Unicode(phi_mask));
  mask.enabled = inum(num_or(desc, x_m, "enable", 0)) != 0;
  if (!mask.enabled) return mask;
  mask.mode = inum(num_or(desc, x_m, "overlap", 0)) ? ec::PhiMask::Mode::Overlap : ec::PhiMask::Mode::Centre;
  mask.dphi = num(desc, x_m, "dphi");
  const int nseg = inum(num(desc, x_m, "nseg"));
  std::vector<double> c;
  if (x_m.hasAttr(_Unicode(centers))) c = parse_list(desc, string_attr(desc, x_m, "centers"));
  if (static_cast<int>(c.size()) >= nseg) {
    c.resize(nseg);
  } else {  // evenly spaced from phi0
    const double p0 = num_or(desc, x_m, "phi0", 0.);
    c.clear();
    for (int i = 0; i < nseg; ++i) c.push_back(p0 + 2. * M_PI * i / nseg);
  }
  mask.centres = c;
  return mask;
}

Ref_t create_detector(Detector& desc, xml_h handle, SensitiveDetector sens) {
  xml_det_t x_det = handle;
  const int det_id = x_det.id();
  const std::string det_name = x_det.nameStr();
  sens.setType("calorimeter");

  // ---- module type
  xml_comp_t x_mod = x_det.child(_U(module));
  const std::string mtype = str_or(x_mod, "type", "hex_stack");
  auto it = module_builders().find(mtype);
  if (it == module_builders().end())
    throw std::runtime_error(det_name + ": unknown module type '" + mtype + "'");
  ModuleType mt = it->second(desc, x_mod, sens, det_name);
  const double mod_zfront = num(desc, x_mod, "z_front");

  // ---- lattice, selection, phi mask
  auto sites = make_sites(desc, x_det.child(_Unicode(lattice)));
  xml_comp_t x_sel = x_det.child(_Unicode(selection));
  const double rsel_min = num(desc, x_sel, "rmin"), rsel_max = num(desc, x_sel, "rmax");
  const auto mask = make_phi_mask(desc, x_det);

  std::vector<ec::LatticeSite> chosen;
  double rmin_all = 1e300, rmax_all = 0;
  for (const auto& s : sites) {
    const double rc = std::hypot(s.x, s.y);
    if (!(rc > rsel_min && rc < rsel_max)) continue;
    auto fp = ec::polygon(s.x, s.y, mt.nsides, mt.apothem, mt.phi0);
    if (!mask.accept(s.x, s.y, fp)) continue;
    auto [r1, r2] = ec::polygon_rminmax(fp);
    rmin_all = std::min(rmin_all, r1);
    rmax_all = std::max(rmax_all, r2);
    chosen.push_back(s);
  }
  if (chosen.empty()) throw std::runtime_error(det_name + ": no module selected");

  // ---- sheets (full or phi-masked annuli, e.g. the preshower Pb radiator)
  struct SheetDef {
    std::string name;
    Material mat;
    double rmin, rmax, zf, t;
    std::string vis;
  };
  std::vector<SheetDef> sheets;
  double zmin = mod_zfront, zmax = mod_zfront + mt.length;
  for (xml_coll_t si(x_det, _Unicode(sheet)); si; ++si) {
    xml_comp_t x_s = si;
    SheetDef sd{x_s.nameStr(),          desc.material(x_s.materialStr()), num(desc, x_s, "rmin"),
                num(desc, x_s, "rmax"),  num(desc, x_s, "z_front"),        num(desc, x_s, "thickness"),
                x_s.visStr()};
    if (sd.t <= 0) continue;
    zmin = std::min(zmin, sd.zf);
    zmax = std::max(zmax, sd.zf + sd.t);
    rmin_all = std::min(rmin_all, sd.rmin);
    rmax_all = std::max(rmax_all, sd.rmax);
    sheets.push_back(sd);
  }

  // ---- envelope
  std::string env_type = "tube";
  double margin = 0.1 * dd4hep::mm;
  Material env_mat = desc.air();
  if (x_det.hasChild(_Unicode(envelope))) {
    xml_comp_t x_env = x_det.child(_Unicode(envelope));
    env_type = str_or(x_env, "type", env_type);
    margin = num_or(desc, x_env, "margin", margin);
    if (x_env.hasAttr(_U(material))) env_mat = desc.material(x_env.materialStr());
  }
  const double zc = 0.5 * (zmin + zmax), hz = 0.5 * (zmax - zmin);
  DetElement sdet(det_name, det_id);
  Volume env;
  if (env_type == "assembly") {
    env = Assembly(det_name + "_envelope");
  } else {
    env = Volume(det_name + "_envelope", Tube(std::max(0., rmin_all - margin), rmax_all + margin, hz), env_mat);
  }
  env.setAttributes(desc, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  // ---- sheets
  for (auto& sd : sheets) {
    const double zpos = sd.zf + sd.t / 2. - zc;
    if (!mask.enabled) {
      Volume v(det_name + "_" + sd.name, Tube(sd.rmin, sd.rmax, sd.t / 2.), sd.mat);
      v.setVisAttributes(desc, sd.vis);
      env.placeVolume(v, Position(0, 0, zpos));
    } else {  // one segment per phi window
      Volume v(det_name + "_" + sd.name, Tube(sd.rmin, sd.rmax, sd.t / 2., -mask.dphi / 2., mask.dphi / 2.), sd.mat);
      v.setVisAttributes(desc, sd.vis);
      for (double c : mask.centres) env.placeVolume(v, Transform3D(RotationZ(c), Position(0, 0, zpos)));
    }
  }

  // ---- modules: one shared logical volume, one DetElement per module (ANL id)
  const double zmod = mod_zfront + mt.length / 2. - zc;
  for (const auto& s : chosen) {
    PlacedVolume pv = env.placeVolume(mt.volume, Position(s.x, s.y, zmod));
    pv.addPhysVolID("sector", s.sector).addPhysVolID("module", s.id);
    DetElement de(sdet, _toString(s.id, "module%d"), s.id);
    de.setPlacement(pv);
  }

  Volume mother = desc.pickMotherVolume(sdet);
  PlacedVolume env_pv = mother.placeVolume(env, Position(0, 0, zc));
  env_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_pv);

  std::string mask_str = "off";
  if (mask.enabled) {
    std::stringstream ms;
    ms << mask.centres.size() << " x " << mask.dphi / dd4hep::deg << " deg ("
       << (mask.mode == ec::PhiMask::Mode::Centre ? "centre" : "overlap") << "), centres";
    for (double c : mask.centres) ms << " " << c / dd4hep::deg;
    mask_str = ms.str();
  }
  printout(INFO, "SoLID_HexShashlykEMCal",
           "%s: %zu of %zu lattice modules (centre R in (%.2f, %.2f) cm, phi mask %s); module %s length %.3f cm, "
           "z %.3f..%.3f cm; envelope r %.3f..%.3f cm",
           det_name.c_str(), chosen.size(), sites.size(), rsel_min / dd4hep::cm, rsel_max / dd4hep::cm,
           mask_str.c_str(), mtype.c_str(), mt.length / dd4hep::cm, zmin / dd4hep::cm, zmax / dd4hep::cm,
           rmin_all / dd4hep::cm, rmax_all / dd4hep::cm);
  return sdet;
}

} // namespace

DECLARE_DETELEMENT(SoLID_HexShashlykEMCal, create_detector)
