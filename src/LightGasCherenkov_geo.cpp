// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_LightGasCherenkov: the full-SoLID light gas Cherenkov (LGC), SIDIS and PVDIS configurations.
//
// Reference geometry: solid_gemc/geometry/lgc_moved (extZ_smallAng optics, standard-nose tank), see
// compact/cherenkov/lgc_SIDIS.xml for the parameter documentation. All numbers come from the compact file.
//
// Volume / DetElement tree
//   LightGasCherenkov            gas Polycone placed in the world (system)
//   |- window_front/back         PVF disks at the tank faces (radii follow the Polycone, no DetElement)
//   `- sector<s>, s = 1..N       Assembly, RotationZ(phi1 + dphi (s-1) - 90 deg)          [sector = s]
//      |- mirror1, mirror2       Sphere shell ∩ ConeSegment (+ optional rigid tilt about a pivot), skin surface
//      `- pmt_assembly           Assembly, RotationZ(dphi_pmt) * [RotationX(rotx), image point]
//         |- winston_cone        straight cone, skin surface
//         |- shield(s)           tubes
//         `- pmt<p>, p = 1..9    gas Box envelope (DetElement)                             [pmt = nx*... see below]
//            |- window           UV glass
//            |- pixel x 64       photocathode, sensitive                                   [pixel]
//            `- lightstop        absorbing skin surface
// Numbering (identical to GEMC solid_lgc): pmt = j*ny + k + 1, pixel = jp*npy + kp + 1, where j (jp) runs along
// the local x axis and k (kp) along the local y axis of the PMT frame, starting at the most negative coordinate.
// The design frame is the sector centred on +y (phi = 90 deg); the sector DetElement is built once with all its
// children and cloned for every sector, so the VolumeManager resolves every pixel of every sector.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"

#include <fmt/core.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

struct ZPlane {
  double z, rmin, rmax;
};

// Radii of a polycone at z, taken on the first segment of non-zero length that contains z.
std::pair<double, double> polyconeRadii(const std::vector<ZPlane>& p, double z) {
  for (size_t i = 0; i + 1 < p.size(); ++i) {
    const auto& a = p[i];
    const auto& b = p[i + 1];
    if (b.z > a.z && z >= a.z && z <= b.z) {
      double f = (z - a.z) / (b.z - a.z);
      return {a.rmin + f * (b.rmin - a.rmin), a.rmax + f * (b.rmax - a.rmax)};
    }
  }
  throw std::runtime_error(fmt::format("SoLID_LightGasCherenkov: z = {} cm outside the tank", z / cm));
}

std::string strAttr(xml_comp_t x, const char* name, const std::string& def = "") {
  return x.hasAttr(Unicode(name)) ? x.attr<std::string>(Unicode(name)) : def;
}

double dblAttr(xml_comp_t x, const char* name, double def) {
  return x.hasAttr(Unicode(name)) ? x.attr<double>(Unicode(name)) : def;
}

} // namespace

static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {
  xml::DetElement x_det = handle;
  const std::string det_name = x_det.nameStr();
  const int det_id = x_det.id();
  DetElement det(det_name, det_id);
  sens.setType("tracker");
  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // ------------------------------------------------------------------------------------------------
  // tank: radiator gas Polycone (GEMC builds no walls) with PVF windows at its faces
  xml_comp_t x_tank = x_det.child(_Unicode(tank));
  std::vector<ZPlane> planes;
  for (xml_coll_t ip(x_tank, _Unicode(zplane)); ip; ++ip) {
    xml_comp_t xp = ip;
    planes.push_back({xp.z(), xp.rmin(), xp.rmax()});
  }
  if (planes.size() < 2) {
    throw std::runtime_error("SoLID_LightGasCherenkov: the tank needs at least two zplanes");
  }
  std::vector<double> pz, prmin, prmax;
  for (const auto& p : planes) {
    pz.push_back(p.z);
    prmin.push_back(p.rmin);
    prmax.push_back(p.rmax);
  }
  Material gas = desc.material(x_tank.materialStr());
  Polycone tank_solid(0., 2. * M_PI, prmin, prmax, pz);
  Volume v_gas(det_name + "_gas", tank_solid, gas);
  v_gas.setVisAttributes(desc.visAttributes(strAttr(x_tank, "vis")));

  PlacedVolume pv_det = desc.pickMotherVolume(det).placeVolume(v_gas, Position(0., 0., 0.));
  pv_det.addPhysVolID("system", det_id);
  det.setPlacement(pv_det);

  for (xml_coll_t iw(x_tank, _Unicode(window)); iw; ++iw) {
    xml_comp_t xw = iw;
    const std::string side = strAttr(xw, "side");
    const double t = xw.thickness();
    double z0 = 0., z1 = 0.;
    if (side == "front") {
      z0 = planes.front().z;
      z1 = z0 + t;
    } else if (side == "back") {
      z1 = planes.back().z;
      z0 = z1 - t;
    } else {
      throw std::runtime_error("SoLID_LightGasCherenkov: window side must be front or back");
    }
    // a cone whose radii follow the tank profile on both faces: covers the full face, no extrusion
    auto [rmin0, rmax0] = polyconeRadii(planes, z0);
    auto [rmin1, rmax1] = polyconeRadii(planes, z1);
    Cone w_solid(t / 2., rmin0, rmax0, rmin1, rmax1);
    Volume v_w(det_name + "_" + strAttr(xw, "name", "window_" + side), w_solid, desc.material(xw.materialStr()));
    v_w.setVisAttributes(desc.visAttributes(strAttr(xw, "vis")));
    v_gas.placeVolume(v_w, Position(0., 0., (z0 + z1) / 2.));
  }

  // ------------------------------------------------------------------------------------------------
  // sector prototype (design frame: sector centred on +y)
  xml_comp_t x_sec = x_det.child(_Unicode(sectors));
  const int nsec = x_sec.number();
  const double phi1 = x_sec.attr<double>(_Unicode(phi1));
  const double dphi = x_sec.attr<double>(_Unicode(dphi));
  Assembly v_sector(det_name + "_sector");
  v_sector.setVisAttributes(desc.visAttributes(strAttr(x_sec, "vis")));
  DetElement de_proto(det_name + "_sector_proto", 0);

  // ---- mirrors
  for (xml_coll_t im(x_det, _Unicode(mirror)); im; ++im) {
    xml_comp_t xm = im;
    const int mid = xm.id();
    const std::string mname = strAttr(xm, "name", fmt::format("mirror{}", mid));

    xml_comp_t xs = xm.child(_Unicode(sphere));
    const double R = xs.rmin();
    const double thk = xs.thickness();
    const Position C(xs.x(0.), xs.y(), xs.z());
    Sphere sph(R, R + thk, dblAttr(xs, "theta0", 0.), dblAttr(xs, "theta1", M_PI / 2.));

    // cut: 'dphi' wedge centred on +y; inner/outer cone surfaces are the lines from apex_z at theta_min/theta_max
    xml_comp_t xc = xm.child(_Unicode(cut));
    const double zc = xc.z();
    const double hz = xc.dz();
    const double apex = xc.attr<double>(_Unicode(apex_z));
    const double tmin = std::tan(xc.attr<double>(_Unicode(theta_min)));
    const double tmax = std::tan(xc.attr<double>(_Unicode(theta_max)));
    const double l1 = zc - hz - apex, l2 = zc + hz - apex;
    const double cphi0 = M_PI / 2. - dphi / 2., cphi1 = M_PI / 2. + dphi / 2.;
    ConeSegment cut(hz, l1 * tmin, l1 * tmax, l2 * tmin, l2 * tmax, cphi0, cphi1);
    // GEMC "Operation:@": both solids are given in the tank frame; here the solid lives in the sphere frame
    Solid msolid = IntersectionSolid(sph, cut, Position(0., 0., zc) - C);

    // optional rigid tilt about the axis parallel to x through the pivot (0, py, pz) (GEMC mirrorRot)
    double tilt = 0., py = 0., pzp = 0.;
    bool clip = false;
    if (xm.hasChild(_Unicode(tilt))) {
      xml_comp_t xt = xm.child(_Unicode(tilt));
      tilt = xt.attr<double>(_Unicode(angle));
      clip = xt.hasAttr(_Unicode(clip)) ? xt.attr<bool>(_Unicode(clip)) : false;
      py = xt.y();
      pzp = xt.z();
    }
    RotationX rtilt(tilt);
    const Position P(0., py, pzp);
    const Position Cp = (std::abs(tilt) > 0.) ? P + rtilt * (C - P) : C;
    const Transform3D tr_mirror(rtilt, Cp);
    if (std::abs(tilt) > 0. && clip) {
      // a tilted wedge no longer ends on the sector boundary planes; optionally clip it with the untilted sector
      // wedge so that neighbouring sectors cannot overlap (for the PVDIS -8 deg tilt the wedge shrinks by
      // 0.04-0.28 deg on each side, so no clipping is needed: default off)
      const double rbig = 2. * (std::abs(C.y()) + R + thk);
      Tube wedge(0., rbig, rbig, cphi0, cphi1);
      msolid = IntersectionSolid(msolid, wedge, tr_mirror.Inverse());
    }
    Volume v_m(det_name + "_" + mname, msolid, desc.material(xm.materialStr()));
    v_m.setVisAttributes(desc.visAttributes(strAttr(xm, "vis")));
    PlacedVolume pv_m = v_sector.placeVolume(v_m, tr_mirror);
    DetElement de_m(de_proto, mname, mid);
    de_m.setPlacement(pv_m);
    SkinSurface skin(desc, de_m, det_name + "_" + mname + "_skin", surfMgr.opticalSurface(strAttr(xm, "surface")), v_m);
    skin.isValid();
    printout(INFO, "SoLID_LightGasCherenkov",
             "%s: R = %.3f cm, C = (0, %.3f, %.3f) cm, tilt %.2f deg -> placed at (0, %.3f, %.3f) cm; cut %.3f/%.3f -> %.3f/%.3f cm",
             mname.c_str(), R / cm, C.y() / cm, C.z() / cm, tilt / deg, Cp.y() / cm, Cp.z() / cm, l1 * tmin / cm,
             l1 * tmax / cm, l2 * tmin / cm, l2 * tmax / cm);
  }

  // ---- PMT assembly: frame origin at the image point, local z towards the mirrors
  xml_comp_t x_pa = x_det.child(_Unicode(pmt_assembly));
  Assembly v_pa(det_name + "_pmt_assembly");
  const Transform3D tr_pa = Transform3D(RotationZ(dblAttr(x_pa, "dphi", 0.))) *
                            Transform3D(RotationX(dblAttr(x_pa, "rotx", 0.)), Position(x_pa.x(0.), x_pa.y(), x_pa.z()));
  v_sector.placeVolume(v_pa, tr_pa);

  // Winston cone (narrow end at local z = 0)
  {
    xml_comp_t xc = x_pa.child(_Unicode(cone));
    const double len = xc.length();
    Cone c_solid(len / 2., xc.rmin1(), xc.rmax1(), xc.rmin2(), xc.rmax2());
    Volume v_c(det_name + "_winston_cone", c_solid, desc.material(xc.materialStr()));
    v_c.setVisAttributes(desc.visAttributes(strAttr(xc, "vis")));
    v_pa.placeVolume(v_c, Position(0., 0., len / 2.));
    SkinSurface skin(desc, det, det_name + "_winston_cone_skin", surfMgr.opticalSurface(strAttr(xc, "surface")), v_c);
    skin.isValid();
  }
  // magnetic shields
  for (xml_coll_t is(x_pa, _Unicode(shield)); is; ++is) {
    xml_comp_t xs = is;
    const double z0 = xs.attr<double>(_Unicode(zmin)), z1 = xs.attr<double>(_Unicode(zmax));
    Tube s_solid(xs.rmin(), xs.rmax(), (z1 - z0) / 2.);
    Volume v_s(det_name + "_" + strAttr(xs, "name", "shield"), s_solid, desc.material(xs.materialStr()));
    v_s.setVisAttributes(desc.visAttributes(strAttr(xs, "vis")));
    v_pa.placeVolume(v_s, Position(0., 0., (z0 + z1) / 2.));
  }

  // MaPMT: gas envelope [light stop back face, 0] holding window, pixels, light stop
  xml_comp_t x_pmt = x_pa.child(_Unicode(pmt));
  const int nx = x_pmt.attr<int>(_Unicode(nx)), ny = x_pmt.attr<int>(_Unicode(ny));
  const double pitch = x_pmt.attr<double>(_Unicode(pitch));
  const double size = x_pmt.attr<double>(_Unicode(size));
  xml_comp_t x_win = x_pmt.child(_Unicode(window));
  xml_comp_t x_pix = x_pmt.child(_Unicode(pixel));
  xml_comp_t x_stop = x_pmt.child(_Unicode(lightstop));
  const double t_win = x_win.thickness();
  const double t_stop = x_stop.thickness();
  const double z_stop = x_stop.z();
  const double env_z0 = z_stop - t_stop / 2.;
  const double env_hz = -env_z0 / 2., env_zc = env_z0 / 2.;
  if (env_z0 >= -t_win) {
    throw std::runtime_error("SoLID_LightGasCherenkov: the light stop must sit behind the PMT window");
  }

  Box env_solid(pitch / 2., pitch / 2., env_hz);
  Volume v_pmt(det_name + "_pmt", env_solid, gas);
  v_pmt.setVisAttributes(desc.visAttributes(strAttr(x_pmt, "vis")));

  Volume v_win(det_name + "_pmt_window", Box(size / 2., size / 2., t_win / 2.), desc.material(x_win.materialStr()));
  v_pmt.placeVolume(v_win, Position(0., 0., -t_win / 2. - env_zc));

  Volume v_stop(det_name + "_pmt_lightstop", Box(size / 2., size / 2., t_stop / 2.), desc.material(x_stop.materialStr()));
  v_pmt.placeVolume(v_stop, Position(0., 0., z_stop - env_zc));
  {
    SkinSurface skin(desc, det, det_name + "_lightstop_skin", surfMgr.opticalSurface(strAttr(x_stop, "surface")), v_stop);
    skin.isValid();
  }

  const int npx = x_pix.attr<int>(_Unicode(nx)), npy = x_pix.attr<int>(_Unicode(ny));
  const double ppitch = x_pix.attr<double>(_Unicode(pitch));
  const double psize = x_pix.attr<double>(_Unicode(size));
  const double pthk = x_pix.thickness();
  const double pzc = x_pix.z();
  if (pzc + pthk / 2. >= -t_win || pzc - pthk / 2. <= z_stop + t_stop / 2.) {
    throw std::runtime_error("SoLID_LightGasCherenkov: the pixels must sit between the PMT window and the light stop");
  }
  Volume v_pix(det_name + "_pixel", Box(psize / 2., psize / 2., pthk / 2.), desc.material(x_pix.materialStr()));
  v_pix.setVisAttributes(desc.visAttributes(strAttr(x_pix, "vis")));
  if (x_pix.isSensitive()) {
    v_pix.setSensitiveDetector(sens);
  }
  for (int jp = 0; jp < npx; ++jp) {
    for (int kp = 0; kp < npy; ++kp) {
      const double px = (jp - (npx - 1) / 2.) * ppitch, py = (kp - (npy - 1) / 2.) * ppitch;
      PlacedVolume pv = v_pmt.placeVolume(v_pix, Position(px, py, pzc - env_zc));
      pv.addPhysVolID("pixel", jp * npy + kp + 1);
    }
  }
  for (int j = 0; j < nx; ++j) {
    for (int k = 0; k < ny; ++k) {
      const int pid = j * ny + k + 1;
      const double x = (j - (nx - 1) / 2.) * pitch, y = (k - (ny - 1) / 2.) * pitch;
      PlacedVolume pv = v_pa.placeVolume(v_pmt, Position(x, y, env_zc));
      pv.addPhysVolID("pmt", pid);
      DetElement de_pmt(de_proto, fmt::format("pmt{}", pid), pid);
      de_pmt.setPlacement(pv);
    }
  }

  // ------------------------------------------------------------------------------------------------
  // sectors: GEMC numbering, sector s centred at phi1 + dphi (s - 1)
  for (int isec = 1; isec <= nsec; ++isec) {
    const double phi = phi1 + dphi * (isec - 1);
    PlacedVolume pv = v_gas.placeVolume(v_sector, Transform3D(RotationZ(phi - M_PI / 2.)));
    pv.addPhysVolID("sector", isec);
    DetElement de_s = de_proto.clone(fmt::format("sector{}", isec), isec);
    de_s.setPlacement(pv);
    det.add(de_s);
  }
  printout(INFO, "SoLID_LightGasCherenkov", "%s: %d sectors (sector 1 at %.2f deg), %d x %d PMTs of %d x %d pixels",
           det_name.c_str(), nsec, phi1 / deg, nx, ny, npx, npy);
  return det;
}

DECLARE_DETELEMENT(SoLID_LightGasCherenkov, createDetector)
