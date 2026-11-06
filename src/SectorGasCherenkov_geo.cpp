// SPDX-License-Identifier: LGPL-3.0-or-later
//
// SoLID_SectorGasCherenkov: a sectored, axially symmetric gas Cherenkov counter with one spherical mirror,
// one light-collection funnel, one magnetic back shield and one photosensor array per sector.
// Used for the SoLID heavy gas Cherenkov (HGC, GEMC hgc_moved) and the SoLID-Lite CO2 Cherenkov
// (GEMC hgc_moved_ext / hgc_moved_extshift). All numbers come from the compact file.
//
// Volume tree (hall frame = detector frame shifted by <placement z>):
//   tank (Polycone, hall-frame z planes)                    DetElement <name>, volID system
//   +- radiator (Polycone)
//   |  +- block_i       (phi-wedge Polycone, half-tank separators)
//   |  +- window_i      (Tube, if mother="radiator")
//   |  +- mirror_s      (ConeSegment(phi_s +- dphi/2) ∩ Sphere shell)   skin surface
//   |  +- cone_s        (Trd2 - Trd2 funnel)                             skin surface
//   |  +- shield_s      (Box - Box tube behind the photosensors)
//   |  +- pmt_s         (Box container, radiator gas)                    DetElement sector_s, volID sector
//   |     +- absorber   (back half of the container)
//   |     +- tile_k     (nx x ny sensitive tiles, radiator gas)          DetElement pmt_k, volID pmt
//   +- window_i         (Tube, if mother="tank")
//
// Sector convention: sector s = 1..count is centred at phi_s = phi0 + (s-1)*dphi. A component with a
// sector-local pose T_loc (sector axis along +x) is placed with RotationZ(phi_s) * T_loc. The optical axis
// (local +z of cone, shield and photosensor plane) is RotationY(-tilt) applied to +z, i.e. it points to
// (-sin(tilt), 0, cos(tilt)) in (r, ., z): from the photosensor towards the mirror. All rotations are ACTIVE.
//
// Frames: the detector frame is the hall frame shifted by <placement z> (z shift, e.g. -40 cm for the GEMC
// hgc_moved_extshift lite variant). z planes and windows are given in the detector frame unless frame="hall".
// Options: an attribute that names a constant directly uses the FIRST definition of that constant (top-level override);
// string options (material, sector list) may be string constants. Sector list: "all" or numbers ("4 5 11 12").
//
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

struct ZPlanes {
  std::vector<double> z, rmin, rmax;
  bool empty() const { return z.empty(); }
};

// Strings may be given literally or as the name of a string constant (<constant ... type="string"/>).
std::string resolveString(const Detector& desc, const std::string& s) {
  const auto& consts = desc.constants();
  if (consts.find(s) != consts.end()) {
    try {
      return desc.constantAsString(s);
    } catch (...) {
    }
  }
  return s;
}

std::string strAttr(const Detector& desc, xml_comp_t x, const char* name, const std::string& def = "") {
  if (!x.hasAttr(Unicode(name))) return def;
  return resolveString(desc, x.attr<std::string>(Unicode(name)));
}

// Numeric attribute. When the attribute text is exactly the name of a constant, the value of the FIRST definition
// of that constant is used (DD4hep keeps the first <constant> object but its expression evaluator keeps the last
// one), so that a top-level file can override an option by defining the constant before including the fragment.
// Any other text is evaluated as an expression as usual.
double num(const Detector& desc, xml_comp_t x, const char* name) {
  std::string raw = x.attr<std::string>(Unicode(name));
  raw.erase(0, raw.find_first_not_of(" \t"));
  raw.erase(raw.find_last_not_of(" \t") + 1);
  const auto& consts = desc.constants();
  if (consts.find(raw) != consts.end()) return _toDouble(desc.constantAsString(raw));
  return x.attr<double>(Unicode(name));
}

double num(const Detector& desc, xml_comp_t x, const char* name, double def) {
  return x.hasAttr(Unicode(name)) ? num(desc, x, name) : def;
}

// z planes; frame="hall" planes are given in the hall frame and moved by -zshift into the detector frame
ZPlanes readZPlanes(const Detector& desc, xml_comp_t x, double zshift) {
  ZPlanes p;
  for (xml_coll_t iz(x, _U(zplane)); iz; ++iz) {
    xml_comp_t xz = iz;
    const bool hall = strAttr(desc, xz, "frame", "detector") == "hall";
    p.z.push_back(num(desc, xz, "z") - (hall ? zshift : 0.));
    p.rmin.push_back(num(desc, xz, "rmin"));
    p.rmax.push_back(num(desc, xz, "rmax"));
  }
  return p;
}

// "all", or a list of sector numbers separated by blanks/commas, possibly via a string constant.
std::vector<int> parseSectors(const Detector& desc, const std::string& raw, int count) {
  std::string s = resolveString(desc, raw);
  std::replace(s.begin(), s.end(), ',', ' ');
  std::vector<int> out;
  std::istringstream is(s);
  std::string tok;
  std::set<int> seen;
  while (is >> tok) {
    if (tok == "all" || tok == "ALL") {
      for (int i = 1; i <= count; ++i)
        if (seen.insert(i).second) out.push_back(i);
      continue;
    }
    int id = std::stoi(tok);
    if (id < 1 || id > count)
      throw std::runtime_error("SoLID_SectorGasCherenkov: sector " + tok + " outside 1.." + std::to_string(count));
    if (seen.insert(id).second) out.push_back(id);
  }
  std::sort(out.begin(), out.end());
  return out;
}

// normalize an angle to [0, 2pi)
double norm2pi(double a) {
  a = std::fmod(a, 2. * M_PI);
  return a < 0 ? a + 2. * M_PI : a;
}

} // namespace

static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {
  xml::DetElement x_det = handle;
  const std::string name = x_det.nameStr();
  const int det_id = x_det.id();
  DetElement det(name, det_id);
  sens.setType("tracker");

  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // ------------------------------------------------------------------ global placement (z shift)
  double zshift = 0.;
  if (x_det.hasChild(_U(placement))) zshift = num(desc, x_det.child(_U(placement)), "z", 0.);

  // ------------------------------------------------------------------ tank and radiator
  xml_comp_t x_tank = x_det.child(_Unicode(tank));
  ZPlanes tank_planes = readZPlanes(desc, x_tank, zshift);
  Volume v_tank(name + "_tank", Polycone(0., 2. * M_PI, tank_planes.rmin, tank_planes.rmax, tank_planes.z),
                desc.material(strAttr(desc, x_tank, "material")));
  v_tank.setVisAttributes(desc, strAttr(desc, x_tank, "vis", x_det.visStr()));

  xml_comp_t x_gas = x_det.child(_Unicode(radiator));
  ZPlanes gas_planes = readZPlanes(desc, x_gas, zshift);
  const std::string gas_mat_name = strAttr(desc, x_gas, "material");
  Material gas_mat = desc.material(gas_mat_name);
  Volume v_gas(name + "_radiator", Polycone(0., 2. * M_PI, gas_planes.rmin, gas_planes.rmax, gas_planes.z), gas_mat);
  v_gas.setVisAttributes(desc, strAttr(desc, x_gas, "vis", "InvisibleWithDaughters"));
  v_tank.placeVolume(v_gas, Position(0., 0., 0.));

  // ------------------------------------------------------------------ windows (any number)
  int iwin = 0;
  for (xml_coll_t iw(x_det, _Unicode(window)); iw; ++iw, ++iwin) {
    xml_comp_t xw = iw;
    const double thk = num(desc, xw, "thickness");
    Volume v_win(name + "_window_" + std::to_string(iwin), Tube(num(desc, xw, "rmin"), num(desc, xw, "rmax"), thk / 2.),
                 desc.material(strAttr(desc, xw, "material")));
    v_win.setVisAttributes(desc, strAttr(desc, xw, "vis", x_det.visStr()));
    const std::string mother = strAttr(desc, xw, "mother", "tank");
    // z = upstream face of the window
    const bool hall = strAttr(desc, xw, "frame", "detector") == "hall";
    Position pos(0., 0., num(desc, xw, "z") - (hall ? zshift : 0.) + thk / 2.);
    if (mother == "radiator")
      v_gas.placeVolume(v_win, pos);
    else
      v_tank.placeVolume(v_win, pos);
  }

  // ------------------------------------------------------------------ phi-wedge support blocks
  std::vector<std::pair<double, double>> blocks; // [start, start+dphi) normalized
  int iblk = 0;
  for (xml_coll_t ib(x_det, _Unicode(block)); ib; ++ib, ++iblk) {
    xml_comp_t xb = ib;
    ZPlanes bp = readZPlanes(desc, xb, zshift);
    if (bp.empty()) bp = gas_planes;
    const double phi = num(desc, xb, "phi"), dphi = num(desc, xb, "dphi");
    Volume v_blk(name + "_block_" + std::to_string(iblk), Polycone(phi, dphi, bp.rmin, bp.rmax, bp.z),
                 desc.material(strAttr(desc, xb, "material")));
    v_blk.setVisAttributes(desc, strAttr(desc, xb, "vis", x_det.visStr()));
    v_gas.placeVolume(v_blk, Position(0., 0., 0.));
    blocks.emplace_back(norm2pi(phi), dphi);
  }

  // ------------------------------------------------------------------ sectors
  xml_comp_t x_sec = x_det.child(_Unicode(sectors));
  const int nsec = int(std::lround(num(desc, x_sec, "count")));
  const double phi0 = num(desc, x_sec, "phi0");
  const double sec_dphi = num(desc, x_sec, "dphi", 2. * M_PI / nsec);
  const std::vector<int> sectors = parseSectors(desc, strAttr(desc, x_sec, "list", "all"), nsec);

  // optical image point (sector-local r, z) and tilt of the optical axis
  xml_comp_t x_img = x_det.child(_Unicode(image));
  const double img_r = num(desc, x_img, "r"), img_z = num(desc, x_img, "z");
  const double tilt = num(desc, x_img, "tilt");
  const RotationY rot_tilt(-tilt);
  // position along the optical axis (offset measured from the image point toward the mirror)
  auto axisPos = [&](double offset) {
    return Position(img_r - offset * std::sin(tilt), 0., img_z + offset * std::cos(tilt));
  };

  // ---- mirror: ConeSegment(phi range) ∩ spherical shell, sector-local, placed at (0,0,mirror z)
  xml_comp_t x_mir = x_det.child(_Unicode(mirror));
  const double mir_R = num(desc, x_mir, "R");
  const double mir_thk = num(desc, x_mir, "thickness");
  const double mir_cr = num(desc, x_mir, "center_r");
  const double mir_cz = num(desc, x_mir, "center_z");
  const double mir_z = num(desc, x_mir, "z");
  const double mir_dz = num(desc, x_mir, "dz");
  const double mir_rmin1 = num(desc, x_mir, "rmin1"), mir_rmax1 = num(desc, x_mir, "rmax1");
  const double mir_rmin2 = num(desc, x_mir, "rmin2"), mir_rmax2 = num(desc, x_mir, "rmax2");
  const double mir_dphi = num(desc, x_mir, "dphi", sec_dphi);
  const double blk_clear = num(desc, x_mir, "block_clearance", 0.);
  Material mir_mat = desc.material(strAttr(desc, x_mir, "material"));
  const std::string mir_vis = strAttr(desc, x_mir, "vis", x_det.visStr());
  const std::string mir_surf_name = strAttr(desc, x_mir, "surface");
  // shell: inner radius = reflecting surface (concave toward the target), polar angles 0-90 deg as GEMC
  Sphere mir_sphere(mir_R, mir_R + mir_thk, 0., M_PI / 2., 0., 2. * M_PI);
  const Position mir_sphere_off(mir_cr, 0., mir_cz - mir_z);
  std::map<std::pair<long, long>, Volume> mirror_vols; // keyed by (lo, hi) in micro-radians
  int n_trimmed = 0;
  auto mirrorVolume = [&](double lo, double hi) -> Volume {
    auto key = std::make_pair(std::lround(lo * 1e6), std::lround(hi * 1e6));
    auto it = mirror_vols.find(key);
    if (it != mirror_vols.end()) return it->second;
    ConeSegment cs(mir_dz, mir_rmin1, mir_rmax1, mir_rmin2, mir_rmax2, lo, hi);
    IntersectionSolid sol(cs, mir_sphere, mir_sphere_off);
    Volume v(name + "_mirror_" + std::to_string(mirror_vols.size()), sol, mir_mat);
    v.setVisAttributes(desc, mir_vis);
    if (!mir_surf_name.empty())
      SkinSurface(desc, det, name + "_mirror_skin_" + std::to_string(mirror_vols.size()),
                  surfMgr.opticalSurface(mir_surf_name), v);
    mirror_vols[key] = v;
    return v;
  };

  // ---- light-collection funnel: Trd2(outer) - Trd2(inner); inner half sizes given, wall added to all four
  xml_comp_t x_cone = x_det.child(_Unicode(cone));
  const double c_x1 = num(desc, x_cone, "x1"), c_x2 = num(desc, x_cone, "x2"), c_y1 = num(desc, x_cone, "y1"), c_y2 = num(desc, x_cone, "y2");
  const double c_dz = num(desc, x_cone, "dz"), c_wall = num(desc, x_cone, "thickness");
  const double c_off = num(desc, x_cone, "offset", c_dz); // cone centre along the axis from the image point
  const double eps = num(desc, x_cone, "boolean_epsilon", 0.01 * dd4hep::cm);
  const double sx = (c_x2 - c_x1) / (2. * c_dz), sy = (c_y2 - c_y1) / (2. * c_dz);
  Trd2 cone_out(c_x1 + c_wall, c_x2 + c_wall, c_y1 + c_wall, c_y2 + c_wall, c_dz);
  // subtracted solid made longer by eps keeping the same wall planes (no coincident faces)
  Trd2 cone_in(c_x1 - sx * eps, c_x2 + sx * eps, c_y1 - sy * eps, c_y2 + sy * eps, c_dz + eps);
  Volume v_cone(name + "_cone", SubtractionSolid(cone_out, cone_in), desc.material(strAttr(desc, x_cone, "material")));
  v_cone.setVisAttributes(desc, strAttr(desc, x_cone, "vis", x_det.visStr()));
  if (x_cone.hasAttr(_Unicode(surface)))
    SkinSurface(desc, det, name + "_cone_skin", surfMgr.opticalSurface(strAttr(desc, x_cone, "surface")), v_cone);

  // ---- back shield: Box(outer) - Box(inner), open square tube along the optical axis
  Volume v_shield;
  double s_off = 0.;
  const bool has_shield = x_det.hasChild(_Unicode(shield));
  if (has_shield) {
    xml_comp_t x_sh = x_det.child(_Unicode(shield));
    const double s_dx = num(desc, x_sh, "dx"), s_dy = num(desc, x_sh, "dy"), s_dz = num(desc, x_sh, "dz"), s_wall = num(desc, x_sh, "thickness");
    s_off = num(desc, x_sh, "offset");
    Box sh_out(s_dx + s_wall, s_dy + s_wall, s_dz);
    Box sh_in(s_dx, s_dy, s_dz + eps);
    v_shield = Volume(name + "_shield", SubtractionSolid(sh_out, sh_in), desc.material(strAttr(desc, x_sh, "material")));
    v_shield.setVisAttributes(desc, strAttr(desc, x_sh, "vis", x_det.visStr()));
  }

  // ---- photosensor array: container box (front face at the image point), absorber back half,
  //      nx x ny sensitive tiles in the front half (Level 0: one tile = whole plane)
  xml_comp_t x_pmt = x_det.child(_Unicode(pmt_array));
  const double p_hx = num(desc, x_pmt, "dx") / 2., p_hy = num(desc, x_pmt, "dy") / 2., p_hz = num(desc, x_pmt, "thickness") / 2.;
  const int p_nx = num(desc, x_pmt, "nx", 1), p_ny = num(desc, x_pmt, "ny", 1);
  const double p_pitchx = num(desc, x_pmt, "pitch_x", num(desc, x_pmt, "dx") / p_nx);
  const double p_pitchy = num(desc, x_pmt, "pitch_y", num(desc, x_pmt, "dy") / p_ny);
  const double t_hx = num(desc, x_pmt, "tile_x", p_pitchx) / 2., t_hy = num(desc, x_pmt, "tile_y", p_pitchy) / 2.;
  const std::string sens_mat_name = strAttr(desc, x_pmt, "sensitive_material", gas_mat_name);
  Volume v_pmt(name + "_pmt", Box(p_hx, p_hy, p_hz), desc.material(strAttr(desc, x_pmt, "material", gas_mat_name)));
  v_pmt.setVisAttributes(desc, strAttr(desc, x_pmt, "vis", x_det.visStr()));
  Volume v_abs(name + "_pmt_absorber", Box(p_hx, p_hy, p_hz / 2.),
               desc.material(strAttr(desc, x_pmt, "absorber_material", "Aluminum")));
  v_abs.setVisAttributes(desc, strAttr(desc, x_pmt, "absorber_vis", "InvisibleNoDaughters"));
  v_pmt.placeVolume(v_abs, Position(0., 0., -p_hz / 2.));
  Volume v_tile(name + "_pmt_tile", Box(t_hx, t_hy, p_hz / 2.), desc.material(sens_mat_name));
  v_tile.setVisAttributes(desc, strAttr(desc, x_pmt, "sensitive_vis", x_det.visStr()));
  v_tile.setSensitiveDetector(sens);
  std::vector<std::pair<int, PlacedVolume>> tile_pvs;
  for (int iy = 0; iy < p_ny; ++iy) {
    for (int ix = 0; ix < p_nx; ++ix) {
      const double x = (ix - 0.5 * (p_nx - 1)) * p_pitchx, y = (iy - 0.5 * (p_ny - 1)) * p_pitchy;
      PlacedVolume pv = v_pmt.placeVolume(v_tile, Position(x, y, p_hz / 2.));
      pv.addPhysVolID("pmt", ix + p_nx * iy);
      tile_pvs.emplace_back(ix + p_nx * iy, pv);
    }
  }

  // sector-local transforms of the per-sector components
  const Transform3D T_cone(rot_tilt, axisPos(c_off));
  const Transform3D T_shield(rot_tilt, axisPos(s_off));
  const Transform3D T_pmt(rot_tilt, axisPos(-p_hz)); // front face at the image point

  for (int isec : sectors) {
    const double phi_s = phi0 + (isec - 1) * sec_dphi;
    const Transform3D T_sec(RotationZ(phi_s), Position(0., 0., 0.));

    // mirror phi range, trimmed against the support blocks (GEMC overlap of mirrors 1/16 with block_1/16)
    double lo = -mir_dphi / 2., hi = mir_dphi / 2.;
    for (auto& [b0, bd] : blocks) {
      // block start/end relative to the sector centre, in (-pi, pi]
      double rel0 = norm2pi(b0 - phi_s + M_PI) - M_PI, rel1 = rel0 + bd;
      if (rel1 <= lo + 1e-9 || rel0 >= hi - 1e-9) continue; // touching only
      if (rel0 <= lo + 1e-9) lo = std::max(lo, rel1 + blk_clear);
      else if (rel1 >= hi - 1e-9) hi = std::min(hi, rel0 - blk_clear);
      else
        printout(WARNING, "SoLID_SectorGasCherenkov", "%s: block inside mirror of sector %d not handled", name.c_str(), isec);
      ++n_trimmed;
    }
    Volume v_mir = mirrorVolume(lo, hi);
    v_gas.placeVolume(v_mir, Transform3D(RotationZ(phi_s), Position(0., 0., mir_z)));

    v_gas.placeVolume(v_cone, T_sec * T_cone);
    if (has_shield) v_gas.placeVolume(v_shield, T_sec * T_shield);

    PlacedVolume pv_pmt = v_gas.placeVolume(v_pmt, T_sec * T_pmt);
    pv_pmt.addPhysVolID("sector", isec);
    DetElement de_sec(det, "sector" + std::to_string(isec), isec);
    de_sec.setPlacement(pv_pmt);
    for (auto& [k, tpv] : tile_pvs) {
      DetElement de_tile(de_sec, "pmt" + std::to_string(k), k);
      de_tile.setPlacement(tpv);
    }
  }

  // ------------------------------------------------------------------ place the tank
  Volume mother = desc.pickMotherVolume(det);
  PlacedVolume pv_tank = mother.placeVolume(v_tank, Position(0., 0., zshift));
  pv_tank.addPhysVolID("system", det_id);
  det.setPlacement(pv_tank);

  std::ostringstream os;
  for (int isec : sectors) os << " " << isec;
  printout(INFO, "SoLID_SectorGasCherenkov",
           "%s: radiator %s, %zu/%d sectors [%s ], phi0 %.2f deg, %d mirror volume(s) (%d trimmed at blocks), "
           "%dx%d sensitive tiles/sector, z shift %.2f cm",
           name.c_str(), gas_mat_name.c_str(), sectors.size(), nsec, os.str().c_str(), phi0 / dd4hep::deg,
           int(mirror_vols.size()), n_trimmed, p_nx, p_ny, zshift / dd4hep::cm);
  return det;
}

DECLARE_DETELEMENT(SoLID_SectorGasCherenkov, createDetector)
