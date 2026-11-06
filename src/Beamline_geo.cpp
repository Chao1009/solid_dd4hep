// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_Beamline: generic beam-pipe builder (entrance and exit pipes of every SoLID configuration).
//
// Construction (follows GEMC geometry/beamline_moved): every <pipe> is a solid wall polycone
// (planes z, rout) with a fill (vacuum) polycone daughter of the same z planes (planes z, rvac),
// so the pipe ends are open like in GEMC. <window> and <kill> disks are daughters of the fill;
// their radius follows the bore profile (rvac) clipped to [zmin, zmax] unless r= is given.
//
//   <detector id="Beamline_ID" name="Beamline" type="SoLID_Beamline">
//     <pipe name="Exit" material="Aluminum" fill="Vacuum" vis="..." fill_vis="...">
//       <zplane z="-325*cm" rout="1.05*cm" rvac="1.0*cm"/> ...           (z increasing, >= 2 planes)
//       <window name="ExitWindow" material="Beryllium" zmin=".." zmax=".." [r=".."] vis=".."/>
//       <kill name="KillDisk" zmin=".." zmax=".." [material="Vacuum"] limits="Beamline_KillLimits"/>
//     </pipe>
//   </detector>
//
// <kill> replaces the GEMC "Kryptonite" disk (gemc kills every track entering a Kryptonite volume):
// here it is a fill-material disk with a DD4hep LimitSet (track_length_max = 0) that makes
// G4UserSpecialCuts (registered by G4StepLimiterPhysics, which ddsim always adds) stop every track
// that enters it.  All numbers come from the compact file; there are no readouts.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;

namespace {

/// linear interpolation of a piecewise-linear profile r(z) (z increasing)
double profile_at(const std::vector<double>& z, const std::vector<double>& r, double zz)
{
  if (zz <= z.front()) {
    return r.front();
  }
  if (zz >= z.back()) {
    return r.back();
  }
  for (size_t i = 1; i < z.size(); ++i) {
    if (zz <= z[i]) {
      const double dz = z[i] - z[i - 1];
      if (dz <= 0.) {
        return r[i];
      }
      return r[i - 1] + (r[i] - r[i - 1]) * (zz - z[i - 1]) / dz;
    }
  }
  return r.back();
}

/// planes of the profile r(z) restricted to [z0, z1] (end points interpolated, inner break points kept)
void clip_profile(const std::vector<double>& z, const std::vector<double>& r, double z0, double z1,
                  std::vector<double>& zo, std::vector<double>& ro)
{
  zo.clear();
  ro.clear();
  zo.push_back(z0);
  ro.push_back(profile_at(z, r, z0));
  for (size_t i = 0; i < z.size(); ++i) {
    if (z[i] > z0 && z[i] < z1) {
      zo.push_back(z[i]);
      ro.push_back(r[i]);
    }
  }
  zo.push_back(z1);
  ro.push_back(profile_at(z, r, z1));
}

} // namespace

static Ref_t create_beamline(Detector& description, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  DetElement det(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");
  assembly.setVisAttributes(description.invisible());

  int npipes = 0;
  for (xml_coll_t ip(x_det, _Unicode(pipe)); ip; ++ip, ++npipes) {
    xml_comp_t x_pipe = ip;
    const std::string pname = det_name + "_" + x_pipe.nameStr();

    std::vector<double> z, rout, rvac;
    for (xml_coll_t iz(x_pipe, _Unicode(zplane)); iz; ++iz) {
      xml_comp_t x_z = iz;
      z.push_back(x_z.z());
      rout.push_back(x_z.attr<double>(_Unicode(rout)));
      rvac.push_back(x_z.attr<double>(_Unicode(rvac)));
    }
    if (z.size() < 2) {
      throw std::runtime_error("SoLID_Beamline: pipe " + pname + " needs at least two <zplane>s");
    }
    for (size_t i = 0; i < z.size(); ++i) {
      if ((i > 0 && z[i] < z[i - 1]) || rvac[i] >= rout[i] || rvac[i] <= 0.) {
        throw std::runtime_error("SoLID_Beamline: pipe " + pname +
                                 ": z planes must increase and 0 < rvac < rout at every plane");
      }
    }
    const std::vector<double> zeros(z.size(), 0.);

    const std::string fill = getAttrOrDefault<std::string>(x_pipe, _Unicode(fill), "Vacuum");
    Volume wall_vol(pname, Polycone(0., 2. * M_PI, zeros, rout, z), description.material(x_pipe.materialStr()));
    wall_vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_pipe, _Unicode(vis), "Beamline_PipeVis"));
    Volume fill_vol(pname + "_fill", Polycone(0., 2. * M_PI, zeros, rvac, z), description.material(fill));
    fill_vol.setVisAttributes(description,
                              getAttrOrDefault<std::string>(x_pipe, _Unicode(fill_vis), "InvisibleWithDaughters"));

    // windows and kill disks: disks filling the bore between zmin and zmax
    auto make_disk = [&](xml_comp_t x_d, const std::string& dname, const std::string& mat) {
      const double zmin = x_d.attr<double>(_Unicode(zmin));
      const double zmax = x_d.attr<double>(_Unicode(zmax));
      if (zmin >= zmax || zmin < z.front() || zmax > z.back()) {
        throw std::runtime_error("SoLID_Beamline: " + dname + " must satisfy zmin < zmax inside the pipe");
      }
      std::vector<double> dz, dr;
      if (x_d.hasAttr(_Unicode(r))) {
        const double r = x_d.attr<double>(_Unicode(r));
        dz = {zmin, zmax};
        dr = {r, r};
        if (r > profile_at(z, rvac, zmin) + 1e-9 || r > profile_at(z, rvac, zmax) + 1e-9) {
          throw std::runtime_error("SoLID_Beamline: " + dname + " radius larger than the bore");
        }
      } else {
        clip_profile(z, rvac, zmin, zmax, dz, dr);
      }
      const std::vector<double> dzero(dz.size(), 0.);
      Volume vol(dname, Polycone(0., 2. * M_PI, dzero, dr, dz), description.material(mat));
      fill_vol.placeVolume(vol);
      printout(DEBUG, "SoLID_Beamline", "%s: %s z %.4f..%.4f cm, r %.4f..%.4f cm", dname.c_str(), mat.c_str(),
               zmin / cm, zmax / cm, dr.front() / cm, dr.back() / cm);
      return vol;
    };

    int iw = 0;
    for (xml_coll_t iwin(x_pipe, _Unicode(window)); iwin; ++iwin, ++iw) {
      xml_comp_t x_w = iwin;
      const std::string wname = pname + "_" + getAttrOrDefault<std::string>(x_w, _U(name), "Window" + std::to_string(iw));
      Volume v = make_disk(x_w, wname, x_w.materialStr());
      v.setVisAttributes(description, getAttrOrDefault<std::string>(x_w, _Unicode(vis), "Beamline_WindowVis"));
    }
    int ik = 0;
    for (xml_coll_t ikill(x_pipe, _Unicode(kill)); ikill; ++ikill, ++ik) {
      xml_comp_t x_k = ikill;
      const std::string kname = pname + "_" + getAttrOrDefault<std::string>(x_k, _U(name), "Kill" + std::to_string(ik));
      Volume v = make_disk(x_k, kname, getAttrOrDefault<std::string>(x_k, _U(material), fill));
      v.setVisAttributes(description, getAttrOrDefault<std::string>(x_k, _Unicode(vis), "Beamline_KillVis"));
      const std::string limits = getAttrOrDefault<std::string>(x_k, _Unicode(limits), "");
      if (!limits.empty()) {
        v.setLimitSet(description, limits);
      } else {
        printout(WARNING, "SoLID_Beamline", "%s has no limits= attribute: it will NOT stop tracks", kname.c_str());
      }
    }

    wall_vol.placeVolume(fill_vol);
    assembly.placeVolume(wall_vol);
    printout(INFO, "SoLID_Beamline", "%s: %s/%s, z %.3f..%.3f cm, rout %.4f..%.4f cm, %d window(s), %d kill disk(s)",
             pname.c_str(), x_pipe.materialStr().c_str(), fill.c_str(), z.front() / cm, z.back() / cm,
             rout.front() / cm, rout.back() / cm, iw, ik);
  }
  if (npipes == 0) {
    printout(WARNING, "SoLID_Beamline", "%s: no <pipe> defined", det_name.c_str());
  }

  PlacedVolume pv = description.pickMotherVolume(det).placeVolume(assembly);
  det.setPlacement(pv);
  return det;
}

DECLARE_DETELEMENT(SoLID_Beamline, create_beamline)
