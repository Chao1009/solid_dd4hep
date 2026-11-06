// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Whitney Armstrong, Zhongling Ji (ePIC-derived TrapEndcapTracker: ACTS DetElement structure)
// Copyright (C) 2026 SoLID collaboration
//
// SoLID_GEMTracker: planar GEM tracking planes of SoLID. One plugin for
//   GEMTracker (full SoLID SIDIS/J/psi/DDVCS 6 planes, PVDIS 5 planes), GEMTrackerLA / GEMTrackerFA (SoLID-Lite)
//   and MuonGEMTracker (DDVCS muon GEMs). All numbers come from the compact XML.
//
// Geometry (solid_gemc gem_moved / gem_separated / muon_moved; docs/geometry_comparison.md):
//   * every plane has `nsectors` azimuthal sectors (SoLID: 30 x 12 deg). Sector s (1..N) is centred at
//         phi_s = phi_sector1 + (s-1) * 360deg/N + phi_offset(layer)        (GEMC: 96 + 12(s-1) + offrot deg)
//     and the module id is the GEMC sector number s;
//   * <layer z=...> is the CENTRE of the module stack (GEMC PlateZ);
//   * the stack is the ordered list of <slice> elements, the first one is the most upstream (-z, facing the target);
//     slice k (1-based) = GEMC slice k; slices with sensitive="true" are sensitive (default: the 3 mm drift gap);
//   * module shape, <module_layout shape="...">:
//       "trapezoid" (default) planar trapezoid that COVERS the GEMC annular sector [rmin, rmax] x dphi:
//                   inner flat at d1 = rmin cos(dphi/2), outer flat at d2 = rmax, sides on the radial lines
//                   phi_s +- dphi/2 (neighbouring modules touch, no overlap).
//                   Module local frame: x = phi-hat, y = +z (beam, stack axis), z = r-hat (narrow end inward).
//                   ACTS: PlaneSurface + TrapezoidBounds (axis_definitions "XZY"); readout CartesianGridXZ.
//       "sector"    tube segment identical to GEMC (rmin..rmax, +-dphi/2), for 1:1 GEMC rate/occupancy checks.
//                   Local frame: x = radial along the sector centre line, y = phi-hat, z = beam;
//                   readout CartesianGridXY. (ACTS would see a disc segment; not meant for the ACTS chain.)
//   * optional <frame> (trapezoid only): border of `material` with perpendicular width side_width along the two
//     radial sides and inner_width / outer_width at the two ends; the slices then fill the inner active trapezoid.
//
// DetElement / volume hierarchy (ePIC TrapEndcapTracker structure, what the ACTS DD4hep conversion expects):
//   <det>                             Assembly (physvolid system), DetType flags from <type_flags>
//    +- <det>_layer<n>                Tube envelope (Air) centred at z_n (physvolid layer = n), layer_material params
//        +- <det>_layer<n>_module<s>  module container (Air; physvolid module = s)
//            +- ..._slice<k>          sensitive slice (physvolid slice = k): DDRec VolPlane + axis_definitions
// The ACTS legacy DD4hep conversion (EICrecon ActsGeometryProvider, acts.examples.dd4hep) only translates endcaps
// that sit in a DD4hep_SubdetectorAssembly ("compound") holding exactly one negative and one positive endcap,
// classified by the sign of the endcap DetElement placement z. <acts side="positive|negative" z_offset="1*um"/>
// places the detector assembly at +-z_offset and compensates the layer z, so hall-frame positions stay exact.
// A detector with no <layer> builds an empty endcap (placeholder partner for a single-endcap tracker).
// The layer_pattern parameter is set with DD4hep_ParametersPlugin in the compact fragment.
//
// Readout id fields used: system, layer, module, slice (+ the segmentation fields).
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"

#include "SoLIDTrackerHelpers.h"

#include <cmath>
#include <string>
#include <vector>

using namespace dd4hep;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::VariantParameters;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;
using SoLIDTrackerHelpers::ensureExtension;
using SoLIDTrackerHelpers::parseIndexList;
using SoLIDTrackerHelpers::xmlToProtoSurfaceMaterial;

namespace {

struct SliceDef {
  int index = 0;       // 1-based (= GEMC slice number)
  std::string name;
  Material material;
  double thickness = 0.;
  double zc = 0.;      // centre along the stack axis, relative to the stack centre (upstream < 0)
  bool sensitive = false;
  std::string vis;
  std::string limits;
};

std::string attrStr(xml_comp_t x, const char* name, const std::string& def) {
  const dd4hep::xml::Strng_t tag(name);
  return x.hasAttr(tag) ? x.attr<std::string>(tag) : def;
}

void setVis(Detector& desc, Volume vol, const std::string& vis) {
  if (vis.empty())
    vol.setVisAttributes(desc.invisible());
  else
    vol.setVisAttributes(desc.visAttributes(vis));
}

} // namespace

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  const int det_id = x_det.id();
  const std::string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& det_params = ensureExtension<VariantParameters>(sdet);
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat)
    xmlToProtoSurfaceMaterial(xml_comp_t(bmat), det_params, "boundary_material");

  Material air = desc.material("Air");
  Assembly assembly(det_name);
  setVis(desc, assembly, x_det.hasAttr(_U(vis)) ? x_det.visStr() : std::string());

  // ------------------------------------------------------------------ ACTS endcap side
  double z_offset = 0.;
  if (x_det.hasChild(_Unicode(acts))) {
    xml_comp_t x_acts = x_det.child(_Unicode(acts));
    const std::string side = attrStr(x_acts, "side", "positive");
    if (side != "positive" && side != "negative")
      throw std::runtime_error(det_name + ": <acts side> must be 'positive' or 'negative', got '" + side + "'");
    const double eps = std::abs(getAttrOrDefault<double>(x_acts, _Unicode(z_offset), 1.0 * dd4hep::um));
    z_offset = (side == "negative" ? -eps : eps);
  }

  // ------------------------------------------------------------------ layout
  std::string shape = "trapezoid";
  int nsectors = 30;
  double phi_sector1 = 96. * dd4hep::deg;
  std::string default_sectors;
  std::string module_vis, layer_vis, slice_vis;
  double pad_r = 0., pad_z = 0.;
  if (x_det.hasChild(_Unicode(module_layout))) {
    xml_comp_t x_lay = x_det.child(_Unicode(module_layout));
    shape = attrStr(x_lay, "shape", shape);
    nsectors = getAttrOrDefault<int>(x_lay, _Unicode(nsectors), nsectors);
    phi_sector1 = getAttrOrDefault<double>(x_lay, _Unicode(phi_sector1), phi_sector1);
    default_sectors = attrStr(x_lay, "sectors", "");
    module_vis = attrStr(x_lay, "module_vis", "");
    layer_vis = attrStr(x_lay, "layer_vis", "");
    slice_vis = attrStr(x_lay, "slice_vis", "");
    pad_r = getAttrOrDefault<double>(x_lay, _Unicode(envelope_pad_r), 0.);
    pad_z = getAttrOrDefault<double>(x_lay, _Unicode(envelope_pad_z), 0.);
  }
  if (shape != "trapezoid" && shape != "sector")
    throw std::runtime_error(det_name + ": module_layout shape must be 'trapezoid' or 'sector', got '" + shape + "'");
  if (nsectors < 3)
    throw std::runtime_error(det_name + ": nsectors must be >= 3");
  const bool trap = (shape == "trapezoid");
  const double dphi = 2. * M_PI / nsectors;
  const double half_phi = 0.5 * dphi;
  if (default_sectors.empty())
    default_sectors = "1-" + std::to_string(nsectors);

  // ------------------------------------------------------------------ module stack
  std::vector<SliceDef> slices;
  double total_thickness = 0.;
  {
    int k = 0;
    for (xml_coll_t xs(x_det, _Unicode(slice)); xs; ++xs) {
      xml_comp_t x_slice = xs;
      SliceDef sl;
      sl.index = ++k;
      sl.name = attrStr(x_slice, "name", "slice" + std::to_string(k));
      sl.material = desc.material(x_slice.materialStr());
      sl.thickness = x_slice.thickness();
      sl.sensitive = x_slice.hasAttr(_U(sensitive)) ? x_slice.isSensitive() : false;
      sl.vis = attrStr(x_slice, "vis", slice_vis);
      sl.limits = attrStr(x_slice, "limits", "");
      if (sl.thickness <= 0.)
        throw std::runtime_error(det_name + ": slice " + sl.name + " has non-positive thickness");
      total_thickness += sl.thickness;
      slices.push_back(sl);
    }
  }
  const double half_thickness = 0.5 * total_thickness;
  {
    double z = -half_thickness;
    for (auto& sl : slices) {
      sl.zc = z + 0.5 * sl.thickness;
      z += sl.thickness;
    }
  }

  const bool have_sens = sens.isValid();
  if (have_sens)
    sens.setType("tracker");

  // ------------------------------------------------------------------ optional frame (trapezoid only)
  bool use_frame = false;
  double frame_side = 0., frame_in = 0., frame_out = 0.;
  Material frame_mat;
  std::string frame_vis;
  if (x_det.hasChild(_Unicode(frame))) {
    xml_comp_t x_fr = x_det.child(_Unicode(frame));
    use_frame = getAttrOrDefault<bool>(x_fr, _Unicode(enable), true);
    if (use_frame && !trap) {
      printout(WARNING, det_name, "<frame> is only implemented for shape=\"trapezoid\"; ignored");
      use_frame = false;
    }
    if (use_frame) {
      frame_side = x_fr.attr<double>(_Unicode(side_width));
      frame_in = getAttrOrDefault<double>(x_fr, _Unicode(inner_width), frame_side);
      frame_out = getAttrOrDefault<double>(x_fr, _Unicode(outer_width), frame_side);
      frame_mat = desc.material(x_fr.materialStr());
      frame_vis = attrStr(x_fr, "vis", "");
    }
  }

  // ------------------------------------------------------------------ layers
  int n_layers = 0;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer = li;
    ++n_layers;
    const int l_id = x_layer.id();
    const double z_layer = x_layer.z();
    const double rmin = x_layer.rmin();
    const double rmax = x_layer.rmax();
    const double phi_off = getAttrOrDefault<double>(x_layer, _Unicode(phi_offset), 0.);
    const std::vector<int> sectors = parseIndexList(attrStr(x_layer, "sectors", default_sectors));
    const std::string layer_name = det_name + "_layer" + std::to_string(l_id);

    if (rmin >= rmax || rmin < 0.)
      throw std::runtime_error(layer_name + ": invalid rmin/rmax");
    if (slices.empty())
      throw std::runtime_error(det_name + ": no <slice> elements");
    for (int sec : sectors)
      if (sec < 1 || sec > nsectors)
        throw std::runtime_error(layer_name + ": sector " + std::to_string(sec) + " outside 1.." +
                                 std::to_string(nsectors));

    // --- module container and slice shapes
    double env_rmin = 0., env_rmax = 0., r_c = 0.;
    double x1 = 0., x2 = 0., hh = 0.;           // trapezoid half widths (inner, outer) and half height
    double ax1 = 0., ax2 = 0., ah = 0., a_dz = 0.; // active trapezoid (= module without frame), shift along r-hat
    Solid module_solid;
    if (trap) {
      const double d1 = rmin * std::cos(half_phi);
      const double d2 = rmax;
      x1 = d1 * std::tan(half_phi);
      x2 = d2 * std::tan(half_phi);
      hh = 0.5 * (d2 - d1);
      r_c = 0.5 * (d1 + d2);
      env_rmin = d1 - pad_r;
      env_rmax = d2 / std::cos(half_phi) + pad_r;
      module_solid = Trapezoid(x1, x2, half_thickness, half_thickness, hh);
      ax1 = x1, ax2 = x2, ah = hh;
      if (use_frame) {
        const double dx = frame_side / std::cos(half_phi);
        ax1 = x1 + frame_in * std::tan(half_phi) - dx;
        ax2 = x2 - frame_out * std::tan(half_phi) - dx;
        ah = hh - 0.5 * (frame_in + frame_out);
        a_dz = 0.5 * (frame_in - frame_out);
        if (ax1 <= 0. || ah <= 0.)
          throw std::runtime_error(layer_name + ": frame larger than the module");
      }
    } else {
      env_rmin = std::max(0., rmin - pad_r);
      env_rmax = rmax + pad_r;
      module_solid = Tube(rmin, rmax, half_thickness, -half_phi, half_phi);
    }
    if (env_rmin < 0.)
      throw std::runtime_error(layer_name + ": envelope rmin < 0");

    Volume module_vol(layer_name + "_module", module_solid, air);
    setVis(desc, module_vol, module_vis);

    if (use_frame) {
      Trapezoid outer(x1, x2, half_thickness, half_thickness, hh);
      // the cut-out is 0.1 mm thicker than the module so that the subtraction has no coincident faces
      Trapezoid inner(ax1, ax2, half_thickness + 0.1 * dd4hep::mm, half_thickness + 0.1 * dd4hep::mm, ah);
      SubtractionSolid frame_solid(outer, inner, Position(0., 0., a_dz));
      Volume frame_vol(layer_name + "_frame", frame_solid, frame_mat);
      setVis(desc, frame_vol, frame_vis);
      module_vol.placeVolume(frame_vol);
    }

    // --- slices; remember the sensitive placements and their measurement surfaces
    std::vector<PlacedVolume> sens_pv;
    std::vector<VolPlane> sens_surf;
    std::vector<int> sens_index;
    for (const auto& sl : slices) {
      const double hs = 0.5 * sl.thickness;
      Solid slice_solid = trap ? Solid(Trapezoid(ax1, ax2, hs, hs, ah)) : Solid(Tube(rmin, rmax, hs, -half_phi, half_phi));
      Volume slice_vol(layer_name + "_slice" + std::to_string(sl.index) + "_" + sl.name, slice_solid, sl.material);
      setVis(desc, slice_vol, sl.vis);
      if (!sl.limits.empty())
        slice_vol.setLimitSet(desc, sl.limits);
      const Position pos = trap ? Position(0., sl.zc, a_dz) : Position(0., 0., sl.zc);
      PlacedVolume spv = module_vol.placeVolume(slice_vol, pos);
      if (!sl.sensitive)
        continue;
      if (!have_sens)
        throw std::runtime_error(det_name + ": sensitive slice " + sl.name + " but no readout given");
      slice_vol.setSensitiveDetector(sens);
      spv.addPhysVolID("slice", sl.index);
      // measurement surface at the slice centre (DDRec: "inner" thickness along -n, "outer" along +n)
      const double up_t = sl.zc + half_thickness;    // material upstream of the surface (target side)
      const double down_t = half_thickness - sl.zc;  // material downstream of the surface
      // trapezoid: u = local x (phi-hat), v = local z (r-hat), n = u x v = -local y (upstream, towards the target);
      //            the same local frame as ACTS'sl "XZY" conversion of the TGeoTrd2 (loc0 = phi-hat, loc1 = r-hat)
      // sector:    u = local x (radial), v = local y (phi-hat), n = local z (+z, downstream)
      const Vector3D u = Vector3D(1., 0., 0.);
      const Vector3D v = trap ? Vector3D(0., 0., 1.) : Vector3D(0., 1., 0.);
      const Vector3D n = trap ? Vector3D(0., -1., 0.) : Vector3D(0., 0., 1.);
      const double inner_t = trap ? down_t : up_t;
      const double outer_t = trap ? up_t : down_t;
      sens_pv.push_back(spv);
      sens_surf.emplace_back(slice_vol, SurfaceType(SurfaceType::Sensitive), inner_t, outer_t, u, v, n);
      sens_index.push_back(sl.index);
    }

    // --- layer envelope
    Tube layer_solid(env_rmin, env_rmax, half_thickness + pad_z);
    Volume layer_vol(layer_name, layer_solid, air);
    setVis(desc, layer_vol, layer_vis);
    PlacedVolume lpv = assembly.placeVolume(layer_vol, Position(0., 0., z_layer - z_offset));
    lpv.addPhysVolID("layer", l_id);
    DetElement layer_de(sdet, layer_name, l_id);
    layer_de.setPlacement(lpv);
    auto& layer_params = ensureExtension<VariantParameters>(layer_de);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat)
      xmlToProtoSurfaceMaterial(xml_comp_t(lmat), layer_params, "layer_material");

    // --- modules (one DetElement per module and per sensitive slice: the VolumeManager resolves every cell)
    for (int sec : sectors) {
      const double phi_s = phi_sector1 + (sec - 1) * dphi + phi_off;
      // trapezoid: local z -> r-hat(phi_s), local y -> +z, local x -> phi-hat(phi_s)
      const Transform3D tr = trap ? Transform3D(RotationZ(phi_s + 0.5 * M_PI) * RotationX(0.5 * M_PI),
                                                Position(r_c * std::cos(phi_s), r_c * std::sin(phi_s), 0.))
                                  : Transform3D(RotationZ(phi_s), Position(0., 0., 0.));
      PlacedVolume mpv = layer_vol.placeVolume(module_vol, tr);
      mpv.addPhysVolID("module", sec);
      const std::string mod_name = layer_name + "_module" + std::to_string(sec);
      DetElement mod_de(layer_de, mod_name, sec);
      mod_de.setPlacement(mpv);
      for (std::size_t i = 0; i < sens_pv.size(); ++i) {
        DetElement sd(mod_de, mod_name + "_slice" + std::to_string(sens_index[i]), sens_index[i]);
        sd.setPlacement(sens_pv[i]);
        auto& sp = ensureExtension<VariantParameters>(sd);
        sp.set<std::string>("axis_definitions", trap ? "XZY" : "XYZ");
        volSurfaceList(sd)->push_back(sens_surf[i]);
      }
    }

    printout(DEBUG, det_name,
             "layer %d: z = %.4f cm, R = %.3f-%.3f cm (%s, envelope %.4f-%.4f cm), %zu of %d sectors, "
             "phi(sector 1) = %.3f deg, stack %zu slices / %.4f mm, %zu sensitive",
             l_id, z_layer / dd4hep::cm, rmin / dd4hep::cm, rmax / dd4hep::cm, shape.c_str(), env_rmin / dd4hep::cm,
             env_rmax / dd4hep::cm, sectors.size(), nsectors, (phi_sector1 + phi_off) / dd4hep::deg, slices.size(),
             total_thickness / dd4hep::mm, sens_pv.size());
  }

  Volume mother = desc.pickMotherVolume(sdet);
  PlacedVolume pv = mother.placeVolume(assembly, Position(0., 0., z_offset));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  printout(INFO, det_name, "%d planes, %s modules, %zu-slice stack of %.4f mm%s", n_layers, shape.c_str(),
           slices.size(), total_thickness / dd4hep::mm,
           z_offset != 0. ? (z_offset > 0. ? ", ACTS positive endcap" : ", ACTS negative endcap") : "");
  return sdet;
}

DECLARE_DETELEMENT(SoLID_GEMTracker, create_detector)
