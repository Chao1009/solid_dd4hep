// SPDX-License-Identifier: LGPL-3.0-or-later
//
// SoLID_PVDISBaffle: the PVDIS "CLEO2" baffle of SoLID
//
//   * 11 slitted absorber plates (default Pb, 9 cm), 30 curved slits each. The lead between two slits ("blade") is
//     bounded by two edges phi0(r), phi1(r) that are piecewise linear in r (Archimedean-spiral pieces), with the
//     breakpoint r1 and the six edge angles phi00..phi12 given per plate by a quadratic function of the plate z
//     (the 21-coefficient "slit model" of solid_gemc/geometry/baffle_moved/solid_PVDIS_baffle_CLEO2_geometry.pl).
//       mode "polygon" (default): every blade is ONE G4ExtrudedSolid (exact edges, sampled every <= dr_sample in r and
//                                 <= dphi_arc on the arcs), placed 30 times inside an air (plate 1: vacuum) envelope
//                                 together with the inner/outer lead rings.
//       mode "gemc"              : the GEMC 20-step staircase (12 deg sector TubeSeg holding nblock lead TubeSegs),
//                                 for bit-level comparisons with GEMC.
//   * plate 1 sits inside the PVDIS target enclosure: with plate1_placement="cavity" (default) a vacuum slot
//     Tube(rinin1, routout1, halfz + gap) is placed in the world at the plate-1 z; the target enclosure must cut the
//     same cavity out of its wall and vacuum (constants BafflePlate1_* of compact/baffle/baffle_PVDIS.xml).
//   * the "12th plane" photon blocker in front of the forward-angle EC (30 lead strips),
//   * the borated-polyethylene neutron shielding (GEMC shldPOLY: inter-plate slabs = sector blocks + rings, inner and
//     outer absorber polycones),
//   * optional observer (flux) planes 20 um thick just upstream/downstream of each plate (readout BaffleObserverHits,
//     fields system, plate, side; GEMC id = (7000 + 10*plate + side)*10000).
//   * chirality = +1 reproduces GEMC (slits drift toward +phi with z: electrons in Bz > 0); -1 mirrors phi -> -phi
//     (reversed solenoid polarity).
// Every number comes from the compact XML. All rotations are active (DD4hep convention).
//
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

// ------------------------------------------------------------------------------------------------ xml helpers
double getD(xml_comp_t x, const char* name, double def) {
  xml::Strng_t tag(name);
  return x.hasAttr(tag) ? x.attr<double>(tag) : def;
}
int getI(xml_comp_t x, const char* name, int def) {
  xml::Strng_t tag(name);
  return x.hasAttr(tag) ? int(std::lround(x.attr<double>(tag))) : def;
}
std::string getS(xml_comp_t x, const char* name, const std::string& def) {
  xml::Strng_t tag(name);
  return x.hasAttr(tag) ? x.attr<std::string>(tag) : def;
}
// "true"/"false"/"on"/"off" or any expression (e.g. a constant) != 0
bool getB(xml_comp_t x, const char* name, bool def) {
  xml::Strng_t tag(name);
  if (!x.hasAttr(tag)) return def;
  std::string s = x.attr<std::string>(tag);
  if (s == "true" || s == "on" || s == "yes") return true;
  if (s == "false" || s == "off" || s == "no") return false;
  return std::fabs(_toDouble(s)) > 0.5;
}
// keyword or number: returns the index of the keyword in `keys`, or the evaluated number
int getEnum(xml_comp_t x, const char* name, const std::vector<std::string>& keys, int def) {
  xml::Strng_t tag(name);
  if (!x.hasAttr(tag)) return def;
  std::string s = x.attr<std::string>(tag);
  for (std::size_t i = 0; i < keys.size(); ++i)
    if (s == keys[i]) return int(i);
  return int(std::lround(_toDouble(s)));
}

// optional override: a constant <prefix><name> defined anywhere in the compact (normally the top-level <define>)
// replaces the attribute value. DD4hep constants cannot be "re-defined" by a later fragment, so the switches are
// not defined in the fragment; the top level defines only the ones it wants to change.
struct Overrides {
  Detector& desc;
  std::string prefix;
  bool has(const std::string& name) const {
    return !prefix.empty() && desc.constants().find(prefix + name) != desc.constants().end();
  }
  double value(const std::string& name) const { return _toDouble(prefix + name); }
  bool b(xml_comp_t x, const char* attr, const std::string& name, bool def) const {
    return has(name) ? std::fabs(value(name)) > 0.5 : getB(x, attr, def);
  }
  int e(xml_comp_t x, const char* attr, const std::string& name, const std::vector<std::string>& keys, int def) const {
    return has(name) ? int(std::lround(value(name))) : getEnum(x, attr, keys, def);
  }
  double d(xml_comp_t x, const char* attr, const std::string& name, double def) const {
    return has(name) ? value(name) : getD(x, attr, def);
  }
};

// ------------------------------------------------------------------------------------------------ slit model
// X(z) = c0 + (c1 + c2*(z - zc0))*(z - zc0) for X in {r1, phi00, phi01, phi02, phi10, phi11, phi12}
const char* const kModelNames[7] = {"r1", "phi00", "phi01", "phi02", "phi10", "phi11", "phi12"};

struct SlitModel {
  double zc0 = 0;
  double c[7][3] = {};
  double eval(int k, double z) const {
    double d = z - zc0;
    return c[k][0] + (c[k][1] + c[k][2] * d) * d;
  }
};

struct Plate {
  int n = 0;
  bool enabled = true;
  double z = 0, rinin = 0, rinout = 0, routin = 0, routout = 0;
  double rot = 0;      // hall angle of the sector-1 local frame (phi_rot0 + offset), before chirality
  double trim = 0;     // hot-edge trim (lead width reduction at the high edge phi1)
  double wadj = 0;     // 360/nsectors - pitch_ref (GEMC "adjust for non 30 sectors")
  int segments = 2;    // 1: one straight r-phi segment r0..r2 (GEMC plates 1, 2), 2: kink at r1
  std::string material, fill;
  // evaluated slit model (local frame of sector 1)
  double r1 = 0, p00 = 0, p01 = 0, p02 = 0, p10 = 0, p11 = 0, p12 = 0;

  // lead [phi0, phi1] at radius r, local frame (GEMC formula, valid for any r)
  std::pair<double, double> edges(double r) const {
    const double r0 = rinout, r2 = routin;
    double a, b;
    if (segments == 1 || r >= r1) {
      a = p01 + (r - r1) * (p02 - p01) / (r2 - r1);
      b = p11 + (r - r1) * (p12 - p11) / (r2 - r1);
    } else {
      a = p01 + (r - r1) * (p00 - p01) / (r0 - r1);
      b = p11 + (r - r1) * (p10 - p11) / (r0 - r1);
    }
    b += wadj - trim;
    return {a, b};
  }
  // hall-frame lead range of sector 1 (before chirality)
  std::pair<double, double> hall_edges(double r) const {
    auto e = edges(r);
    return {e.first + rot, e.second + rot};
  }
};

// mirror a phi range for chirality -1: [a, b] -> [-b, -a]
inline std::pair<double, double> chiral(double a, double b, int chi) {
  return chi > 0 ? std::make_pair(a, b) : std::make_pair(-b, -a);
}

// piecewise-linear interpolation in a (sorted) knot table
double interp(const std::vector<double>& zs, const std::vector<double>& vs, double z) {
  if (z <= zs.front()) return vs.front();
  if (z >= zs.back()) return vs.back();
  auto it = std::upper_bound(zs.begin(), zs.end(), z);
  std::size_t i = std::size_t(it - zs.begin());
  double t = (z - zs[i - 1]) / (zs[i] - zs[i - 1]);
  return vs[i - 1] + t * (vs[i] - vs[i - 1]);
}

} // namespace

// =================================================================================================================
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  const int det_id = x_det.id();
  DetElement det(det_name, det_id);

  xml_comp_t x_lay = x_det.child(_Unicode(layout));
  const Overrides ovr{desc, getS(x_lay, "override_prefix", "")};
  const bool enable = ovr.b(x_lay, "enable", "Enable", true);
  const int nplates = getI(x_lay, "nplates", 11);
  const int nsect = getI(x_lay, "nsectors", 30);
  const double pitch = 2. * M_PI / nsect;
  const double halfz = 0.5 * x_lay.attr<double>(_Unicode(thickness));
  const double zc0 = x_lay.attr<double>(_Unicode(zc0));
  const double dzc = x_lay.attr<double>(_Unicode(dzc));
  const double phi_rot0 = x_lay.attr<double>(_Unicode(phi_rot0));
  const double offset0 = getD(x_lay, "offset", 0.);
  const double pitch_ref = getD(x_lay, "pitch_ref", pitch);
  const double gap = getD(x_lay, "gap", 0.002 * dd4hep::cm);
  const std::string mat_lead = getS(x_lay, "material", "Lead");
  const std::string mat_fill = getS(x_lay, "fill", "Air");
  const std::string mat_fill1 = getS(x_lay, "plate1_fill", "Vacuum");
  const int plate1_place = ovr.e(x_lay, "plate1_placement", "Plate1_Placement", {"cavity", "world", "none"}, 0);
  const int mode = ovr.e(x_lay, "mode", "Mode", {"polygon", "gemc"}, 0);
  const int nblock = getI(x_lay, "nblock", 20);
  const double dr_sample = getD(x_lay, "dr_sample", 1.0 * dd4hep::cm);
  const double dphi_arc = getD(x_lay, "dphi_arc", 0.25 * dd4hep::deg);
  const double inset = getD(x_lay, "inset", 0.001 * dd4hep::cm);
  const double trim0 = getD(x_lay, "hot_edge_trim", 0.);
  int chi = int(std::lround(ovr.d(x_lay, "chirality", "Chirality", 1.)));
  chi = chi < 0 ? -1 : 1;
  const std::string vis_lead = getS(x_lay, "vis", "Baffle_LeadVis");

  // ---- slit model
  SlitModel model;
  model.zc0 = zc0;
  {
    xml_comp_t x_m = x_det.child(_Unicode(slit_model));
    for (int k = 0; k < 7; ++k)
      for (int j = 0; j < 3; ++j) {
        std::string key = "c" + std::to_string(j) + "_" + kModelNames[k];
        xml::Strng_t tag(key);
        if (!x_m.hasAttr(tag))
          throw std::runtime_error(det_name + ": slit_model attribute " + key + " missing");
        model.c[k][j] = x_m.attr<double>(tag);
      }
  }

  // ---- plates
  std::map<int, Plate> plates;
  for (xml_coll_t ip(x_det, _Unicode(plate)); ip; ++ip) {
    xml_comp_t xp = ip;
    Plate p;
    p.n = xp.attr<int>(_Unicode(n));
    p.enabled = getB(xp, "enable", true);
    p.z = getD(xp, "z", zc0 + dzc * (p.n - 0.5 * (nplates + 1)));
    p.rinin = xp.attr<double>(_Unicode(rinin));
    p.rinout = xp.attr<double>(_Unicode(rinout));
    p.routin = xp.attr<double>(_Unicode(routin));
    p.routout = xp.attr<double>(_Unicode(routout));
    p.rot = phi_rot0 + getD(xp, "offset", offset0);
    p.trim = getD(xp, "hot_edge_trim", trim0);
    p.wadj = pitch - pitch_ref;
    p.segments = getI(xp, "segments", 2);
    p.material = getS(xp, "material", mat_lead);
    p.fill = getS(xp, "fill", p.n == 1 ? mat_fill1 : mat_fill);
    p.r1 = model.eval(0, p.z);
    p.p00 = model.eval(1, p.z);
    p.p01 = model.eval(2, p.z);
    p.p02 = model.eval(3, p.z);
    p.p10 = model.eval(4, p.z);
    p.p11 = model.eval(5, p.z);
    p.p12 = model.eval(6, p.z);
    if (p.segments == 2 && !(p.r1 > p.rinout && p.r1 < p.routin)) {
      printout(WARNING, det_name, "plate %d: kink r1 = %.3f cm outside (%.3f, %.3f): using one segment", p.n,
               p.r1 / dd4hep::cm, p.rinout / dd4hep::cm, p.routin / dd4hep::cm);
      p.segments = 1;
    }
    if (plates.count(p.n)) throw std::runtime_error(det_name + ": plate defined twice");
    plates[p.n] = p;
  }

  // the readout (BaffleObserverHits) is declared in the compact even when the observers are off: always give the
  // sensitive detector a type, otherwise ddsim cannot map an action to it (the collection is then simply empty)
  if (sens.isValid()) sens.setType("tracker");

  Assembly assembly(det_name);
  Material vac = desc.material("Vacuum");
  Volume plate1_slot; // vacuum slot of plate 1 (cavity placement)

  // ---------------------------------------------------------------------------------------------- plates
  for (auto& [n, p] : plates) {
    const bool build = enable && p.enabled;
    const bool is_slot = (n == 1 && plate1_place == 0);
    if (n == 1 && plate1_place == 2) continue;
    if (!build && !is_slot) continue;

    const std::string pname = det_name + "_plate" + std::to_string(n);
    DetElement de_plate(det, "plate" + std::to_string(n), n);

    // the vacuum slot of plate 1 (fills the cavity of the target enclosure, holds the plate and its observers)
    Volume slot_vol;
    if (is_slot) {
      Tube slot_shape(p.rinin, p.routout, halfz + gap);
      slot_vol = Volume(pname + "_slot", slot_shape, vac);
      slot_vol.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
      PlacedVolume pv_slot = assembly.placeVolume(slot_vol, Position(0, 0, p.z));
      de_plate.setPlacement(pv_slot);
      plate1_slot = slot_vol;
    }
    if (!build) continue;

    Material m_lead = desc.material(p.material);
    Tube env_shape(p.rinin, p.routout, halfz);
    Volume env(pname, env_shape, desc.material(p.fill));
    env.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));

    // rings
    Volume oring(pname + "_outer_ring", Tube(p.routin, p.routout, halfz), m_lead);
    oring.setVisAttributes(desc.visAttributes(vis_lead));
    env.placeVolume(oring);
    if (p.rinout > p.rinin) {
      Volume iring(pname + "_inner_ring", Tube(p.rinin, p.rinout, halfz), m_lead);
      iring.setVisAttributes(desc.visAttributes(vis_lead));
      env.placeVolume(iring);
    }

    const double r0 = p.rinout, r2 = p.routin;
    if (mode == 0) {
      // ---- exact blade polygon (hall frame of sector 1, clockwise)
      const double ra = r0 + inset, rb = r2 - inset;
      if ((ra * (1. - std::cos(0.5 * dphi_arc))) >= inset)
        printout(WARNING, det_name, "plate %d: inner-arc chord sagitta %.4f mm >= inset %.4f mm (blade may touch the ring)",
                 n, ra * (1. - std::cos(0.5 * dphi_arc)) / dd4hep::mm, inset / dd4hep::mm);
      std::vector<double> knots{ra, rb};
      if (p.segments == 2) knots.push_back(p.r1);
      std::sort(knots.begin(), knots.end());
      std::vector<double> rs;
      for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
        int ns = std::max(1, int(std::ceil((knots[k + 1] - knots[k]) / dr_sample - 1e-9)));
        for (int j = 0; j < ns; ++j) rs.push_back(knots[k] + (knots[k + 1] - knots[k]) * j / ns);
      }
      rs.push_back(knots.back());
      std::vector<std::pair<double, double>> poly; // (r, phi)
      auto arc = [&](double r, double pa, double pb, bool skip_first) {
        int ns = std::max(1, int(std::ceil(std::fabs(pb - pa) / dphi_arc - 1e-9)));
        for (int j = skip_first ? 1 : 0; j <= ns; ++j) poly.emplace_back(r, pa + (pb - pa) * j / ns);
      };
      auto he_b = p.hall_edges(rb), he_a = p.hall_edges(ra);
      arc(rb, he_b.second, he_b.first, false);                         // outer arc, phi decreasing
      for (int k = int(rs.size()) - 2; k >= 0; --k)                     // low edge, rb -> ra
        poly.emplace_back(rs[k], p.hall_edges(rs[k]).first);
      arc(ra, he_a.first, he_a.second, true);                          // inner arc, phi increasing
      for (std::size_t k = 1; k + 1 < rs.size(); ++k)                   // high edge, ra -> rb (open)
        poly.emplace_back(rs[k], p.hall_edges(rs[k]).second);
      std::vector<double> px, py;
      for (auto& [r, ph] : poly) {
        px.push_back(r * std::cos(ph));
        py.push_back(r * std::sin(ph) * chi);
      }
      if (chi < 0) { // keep the clockwise order after the reflection
        std::reverse(px.begin(), px.end());
        std::reverse(py.begin(), py.end());
      }
      ExtrudedPolygon blade_shape(px, py, {-halfz, halfz}, {0., 0.}, {0., 0.}, {1., 1.});
      Volume blade(pname + "_blade", blade_shape, m_lead);
      blade.setVisAttributes(desc.visAttributes(vis_lead));
      for (int i = 0; i < nsect; ++i) env.placeVolume(blade, i + 1, Transform3D(RotationZ(chi * pitch * i)));
      printout(DEBUG, det_name, "plate %d: blade polygon with %zu vertices", n, px.size());
    } else {
      // ---- GEMC staircase: one sector (12 deg container) with nblock lead TubeSegs, placed nsect times
      const double dr = (r2 - r0) / nblock;
      double s0 = p.hall_edges(r0 + 0.5 * dr).first; // container start = phi0 of block 1
      auto cs = chiral(s0, s0 + pitch, chi);
      Volume sector(pname + "_sector", Tube(r0, r2, halfz, cs.first, cs.second), desc.material(p.fill));
      sector.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
      for (int j = 0; j < nblock; ++j) {
        double ra = r0 + j * dr, rb = ra + dr;
        auto he = p.hall_edges(ra + 0.5 * dr);
        auto c = chiral(he.first, he.second, chi);
        Volume blk(pname + "_block" + std::to_string(j + 1), Tube(ra, rb, halfz, c.first, c.second), m_lead);
        blk.setVisAttributes(desc.visAttributes(vis_lead));
        sector.placeVolume(blk, j + 1);
      }
      for (int i = 0; i < nsect; ++i) env.placeVolume(sector, i + 1, Transform3D(RotationZ(chi * pitch * i)));
    }

    if (is_slot) {
      slot_vol.placeVolume(env);
    } else {
      PlacedVolume pv = assembly.placeVolume(env, Position(0, 0, p.z));
      de_plate.setPlacement(pv);
    }
    printout(DEBUG, det_name, "plate %2d z=%.2f r %.2f/%.2f/%.2f/%.2f r1=%.3f cm", n, p.z / dd4hep::cm,
             p.rinin / dd4hep::cm, p.rinout / dd4hep::cm, p.routin / dd4hep::cm, p.routout / dd4hep::cm,
             p.r1 / dd4hep::cm);
  }

  // ---------------------------------------------------------------------------------------------- observers
  if (x_det.hasChild(_Unicode(observers))) {
    xml_comp_t x_obs = x_det.child(_Unicode(observers));
    if (ovr.b(x_obs, "enable", "Observers_Enable", false)) {
      const double oh = getD(x_obs, "halfthickness", 0.001 * dd4hep::cm);
      Material m_obs = desc.material(getS(x_obs, "material", "Vacuum"));
      const std::string vis_obs = getS(x_obs, "vis", "InvisibleNoDaughters");
      for (auto& [n, p] : plates) {
        if (n == 1 && plate1_place == 2) continue;
        if (!(enable && p.enabled) && !(n == 1 && plate1_place == 0)) continue;
        Volume ov(det_name + "_observer" + std::to_string(n), Tube(p.rinin, p.routout, oh), m_obs);
        ov.setSensitiveDetector(sens);
        ov.setVisAttributes(desc.visAttributes(vis_obs));
        const double dz = halfz + 0.5 * gap;
        for (int side = 1; side <= 2; ++side) {
          double zl = (side == 1 ? -dz : dz);
          PlacedVolume pv;
          if (n == 1 && plate1_place == 0) {
            pv = plate1_slot.placeVolume(ov, side, Position(0, 0, zl));
          } else {
            pv = assembly.placeVolume(ov, 10 * n + side, Position(0, 0, p.z + zl));
          }
          pv.addPhysVolID("plate", n).addPhysVolID("side", side);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------------------------- photon blocker
  if (x_det.hasChild(_Unicode(photon_blocker))) {
    xml_comp_t x_pb = x_det.child(_Unicode(photon_blocker));
    if (ovr.b(x_pb, "enable", "PhotonBlocker_Enable", true)) {
      const int ns = getI(x_pb, "nsectors", nsect);
      const double pb_pitch = 2. * M_PI / ns;
      const double z = ovr.d(x_pb, "z", "PhotonBlocker_z", 0.);
      const double hz = 0.5 * x_pb.attr<double>(_Unicode(thickness));
      const double rot = getD(x_pb, "phi_rot0", phi_rot0) + getD(x_pb, "offset", offset0);
      const double s0 = rot + x_pb.attr<double>(_Unicode(phi_start));
      auto c = chiral(s0, s0 + x_pb.attr<double>(_Unicode(dphi)), chi);
      Tube strip_shape(x_pb.attr<double>(_Unicode(rmin)), x_pb.attr<double>(_Unicode(rmax)), hz, c.first, c.second);
      Volume strip(det_name + "_photon_blocker", strip_shape, desc.material(getS(x_pb, "material", mat_lead)));
      strip.setVisAttributes(desc.visAttributes(getS(x_pb, "vis", vis_lead)));
      Assembly pb_asm(det_name + "_PhotonBlocker");
      for (int i = 0; i < ns; ++i) pb_asm.placeVolume(strip, i + 1, Transform3D(RotationZ(chi * pb_pitch * i)));
      PlacedVolume pv = assembly.placeVolume(pb_asm, Position(0, 0, z));
      DetElement de_pb(det, "photon_blocker", 100);
      de_pb.setPlacement(pv);
    }
  }

  // ---------------------------------------------------------------------------------------------- BorPoly shielding
  if (enable && x_det.hasChild(_Unicode(shielding))) {
    xml_comp_t x_sh = x_det.child(_Unicode(shielding));
    if (ovr.b(x_sh, "enable", "Shielding_Enable", true)) {
      Material m_sh = desc.material(x_sh.attr<std::string>(_Unicode(material)));
      Material m_shfill = desc.material(getS(x_sh, "fill", mat_fill));
      const std::string vis_sh = getS(x_sh, "vis", "Baffle_ShieldVis");
      const int first = getI(x_sh, "first_plate", 2);
      const int last = getI(x_sh, "last_plate", 8); // last plate n with a slab between n and n+1
      const int nb = getI(x_sh, "nblock", nblock);
      Assembly sh_asm(det_name + "_Shielding");

      auto ok = [&](int n) { return plates.count(n) && plates[n].enabled; };
      // ---- inter-plate slabs (GEMC staircase sector blocks + inner/outer rings)
      for (int n = first; n <= last; ++n) {
        if (!ok(n) || !ok(n + 1)) continue;
        const Plate &a = plates[n], &b = plates[n + 1];
        const double z0 = a.z + halfz + gap, z1 = b.z - halfz - gap, hz = 0.5 * (z1 - z0), zc = 0.5 * (z0 + z1);
        if (hz <= 0) continue;
        const std::string sn = det_name + "_slab" + std::to_string(n);
        // container (fill) between the rings
        Polycone cont_shape(0, 2 * M_PI, {a.rinout, b.rinout}, {a.routin, b.routin}, {-hz, hz});
        Volume cont(sn, cont_shape, m_shfill);
        cont.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
        // GEMC block radii (staircase) of both plates
        auto blocks = [&](const Plate& p) {
          std::vector<std::array<double, 4>> v; // r0b, r1b, phi0, phi1 (hall)
          const double dr = (p.routin - p.rinout) / nb;
          for (int j = 0; j < nb; ++j) {
            double ra = p.rinout + j * dr;
            auto he = p.hall_edges(ra + 0.5 * dr);
            v.push_back({ra, ra + dr, he.first, he.second});
          }
          return v;
        };
        auto ba = blocks(a), bb = blocks(b);
        const double s0 = 0.5 * (ba[0][2] + bb[0][2]);
        auto cs = chiral(s0, s0 + pitch, chi);
        Polycone sec_shape(cs.first, pitch, {a.rinout, b.rinout}, {a.routin, b.routin}, {-hz, hz});
        Volume sec(sn + "_sector", sec_shape, m_shfill);
        sec.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
        for (int j = 0; j < nb; ++j) {
          double sphi = 0.5 * (ba[j][2] + bb[j][2]);
          double dphi = 0.5 * ((ba[j][3] - ba[j][2]) + (bb[j][3] - bb[j][2]));
          if (sphi < s0 - 1e-9 || sphi + dphi > s0 + pitch + 1e-9)
            printout(WARNING, det_name, "slab %d block %d exceeds its sector container", n, j + 1);
          auto c = chiral(sphi, sphi + dphi, chi);
          Polycone bs(c.first, dphi, {ba[j][0], bb[j][0]}, {ba[j][1], bb[j][1]}, {-hz, hz});
          Volume blk(sn + "_block" + std::to_string(j + 1), bs, m_sh);
          blk.setVisAttributes(desc.visAttributes(vis_sh));
          sec.placeVolume(blk, j + 1);
        }
        for (int i = 0; i < nsect; ++i) cont.placeVolume(sec, i + 1, Transform3D(RotationZ(chi * pitch * i)));
        sh_asm.placeVolume(cont, n, Position(0, 0, zc));
        // rings
        if (a.rinout > a.rinin) {
          Volume ir(sn + "_inner_ring", Polycone(0, 2 * M_PI, {a.rinin, b.rinin}, {a.rinout, b.rinout}, {-hz, hz}), m_sh);
          ir.setVisAttributes(desc.visAttributes(vis_sh));
          sh_asm.placeVolume(ir, n, Position(0, 0, zc));
        }
        Volume orr(sn + "_outer_ring", Polycone(0, 2 * M_PI, {a.routin, b.routin}, {a.routout, b.routout}, {-hz, hz}),
                   m_sh);
        orr.setVisAttributes(desc.visAttributes(vis_sh));
        sh_asm.placeVolume(orr, n, Position(0, 0, zc));
      }

      // ---- absorber z planes: per plate n >= first: z0 = plate front - gap,
      //      z1 = back + gap (plates with a slab), plate z + ext (others), zend (last plate)
      xml_comp_t x_abs = x_sh.child(_Unicode(absorbers));
      if (getB(x_abs, "enable", true)) {
        const double ext = x_abs.attr<double>(_Unicode(downstream_extension));
        const double zend = x_abs.attr<double>(_Unicode(zend));
        std::vector<double> Z, rin_in, rin_out, rout_in;
        int nlast = plates.rbegin()->first;
        for (auto& [n, p] : plates) {
          if (n < first) continue;
          double za = p.z - halfz - gap;
          double zb = (n <= last) ? p.z + halfz + gap : (n < nlast ? p.z + ext : zend);
          Z.push_back(za);
          Z.push_back(zb);
          double rmx = p.rinin - getD(x_abs, "inner_gap", 0.1 * dd4hep::cm);
          rin_out.push_back(rmx);
          rin_out.push_back(rmx);
          rout_in.push_back(p.routout);
          rout_in.push_back(p.routout);
        }
        // inner absorber rmin: GEMC line (z_a, r_a) -> (z_b, r_b), constant beyond, evaluated at the planes;
        // raised to bore(z) + clearance where needed (bore = beam-pipe outer radius profile)
        xml_comp_t x_in = x_abs.child(_Unicode(inner));
        const double lz0 = x_in.attr<double>(_Unicode(rmin_z0)), lr0 = x_in.attr<double>(_Unicode(rmin_r0));
        const double lz1 = x_in.attr<double>(_Unicode(rmin_z1)), lr1 = x_in.attr<double>(_Unicode(rmin_r1));
        auto gline = [&](double z) { return z >= lz1 ? lr1 : lr0 + (lr1 - lr0) / (lz1 - lz0) * (z - lz0); };
        for (double z : Z) rin_in.push_back(gline(z));
        std::vector<double> bz, br;
        for (xml_coll_t ib(x_in, _Unicode(bore)); ib; ++ib) {
          xml_comp_t xb = ib;
          bz.push_back(xb.attr<double>(_Unicode(z)));
          br.push_back(xb.attr<double>(_Unicode(r)));
        }
        const double clr = getD(x_in, "bore_clearance", 0.1 * dd4hep::cm);
        std::vector<double> kz = Z;
        if (!bz.empty()) {
          for (double z : bz)
            if (z > Z.front() && z < Z.back()) kz.push_back(z);
          std::sort(kz.begin(), kz.end());
          kz.erase(std::unique(kz.begin(), kz.end(), [](double u, double v) { return std::fabs(u - v) < 1e-9; }),
                   kz.end());
          // insert crossings of (line - bore - clearance)
          std::vector<double> kz2;
          for (std::size_t k = 0; k + 1 < kz.size(); ++k) {
            kz2.push_back(kz[k]);
            double f0 = interp(Z, rin_in, kz[k]) - interp(bz, br, kz[k]) - clr;
            double f1 = interp(Z, rin_in, kz[k + 1]) - interp(bz, br, kz[k + 1]) - clr;
            if (f0 * f1 < 0) kz2.push_back(kz[k] + (kz[k + 1] - kz[k]) * f0 / (f0 - f1));
          }
          kz2.push_back(kz.back());
          kz = kz2;
        }
        std::vector<double> ri, ro;
        for (double z : kz) {
          double r = interp(Z, rin_in, z);
          if (!bz.empty()) r = std::max(r, interp(bz, br, z) + clr);
          ri.push_back(r);
          ro.push_back(interp(Z, rin_out, z));
        }
        if (getB(x_in, "enable", true)) {
          Volume ia(det_name + "_inner_absorber", Polycone(0, 2 * M_PI, ri, ro, kz), m_sh);
          ia.setVisAttributes(desc.visAttributes(vis_sh));
          sh_asm.placeVolume(ia);
          for (std::size_t k = 0; k < kz.size(); ++k)
            printout(DEBUG, det_name, "inner absorber plane z=%.4f rmin=%.5f rmax=%.5f", kz[k] / dd4hep::cm,
                     ri[k] / dd4hep::cm, ro[k] / dd4hep::cm);
        }
        // outer absorber: rmin = routout of the plate, rmax constant; stop where the thickness vanishes
        xml_comp_t x_out = x_abs.child(_Unicode(outer));
        if (getB(x_out, "enable", true)) {
          const double rmax = x_out.attr<double>(_Unicode(rmax));
          std::vector<double> oz, oi, oo;
          for (std::size_t k = 0; k < Z.size(); ++k) {
            oz.push_back(Z[k]);
            oi.push_back(std::min(rout_in[k], rmax));
            oo.push_back(rmax);
            if (rout_in[k] >= rmax) break;
          }
          Volume oa(det_name + "_outer_absorber", Polycone(0, 2 * M_PI, oi, oo, oz), m_sh);
          oa.setVisAttributes(desc.visAttributes(vis_sh));
          sh_asm.placeVolume(oa);
        }
      }
      PlacedVolume pv = assembly.placeVolume(sh_asm);
      DetElement de_sh(det, "shielding", 200);
      de_sh.setPlacement(pv);
    }
  }

  // ---------------------------------------------------------------------------------------------- place
  PlacedVolume pv = desc.pickMotherVolume(det).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);
  det.setPlacement(pv);
  printout(INFO, det_name, "%zu plates (%s, chirality %+d)%s", plates.size(), mode == 0 ? "polygon" : "gemc", chi,
           enable ? "" : " DISABLED");
  return det;
}

DECLARE_DETELEMENT(SoLID_PVDISBaffle, create_detector)
