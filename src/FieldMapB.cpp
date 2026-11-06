// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck, Dhevan Gangadharan, Aranya Giri; 2026 SoLID collaboration
//
// SoLID_FieldMapB: field map on a regular grid, bilinear (2D, BrBz) or trilinear (3D, BxByBz) interpolation.
//
//   <field type="SoLID_FieldMapB" name="SolenoidField" field_type="magnetic" coord_type="BrBz"
//          field_map="fieldmaps/solenoid_CLEOv9.dat"            path read if it exists (env vars expanded), else
//          cache="${DETECTOR_PATH}/fieldmaps"                    looked up (by file name or URL hash) in the cache
//          url="https://raw.githubusercontent.com/..."           list, else downloaded into the first writable cache
//          scale="1.0">                                          dimensionless multiplier (polarity, current scaling)
//     <dimensions>                                               optional if the file has a GEMC <mfield> header
//       <R step="1*cm" min="0*cm" max="500*cm"/>  <Z step="1*cm" min="-600*cm" max="600*cm"/>   (BrBz)
//       <X .../> <Y .../> <Z .../>                                                             (BxByBz)
//       <translationCoord x="0*cm" y="0*cm" z="0*cm"/>           position of the map origin in the hall
//       <rotationCoord x="0*deg" y="0*deg" z="0*deg"/>           ACTIVE rotation of the map frame, Rz(z)Ry(y)Rx(x)
//       <rotationField x="0" y="0" z="0"/>                       legacy: extra rotation of the B vector only (bare degrees)
//     </dimensions>
//   </field>
//
// File format: whitespace-separated rows "c1 c2 [c3] B1 B2 [B3]" (any row order), optionally preceded by the GEMC
// ASCII header (<mfield> ... </mfield>): <first|second|third npoints min max units>, <field unit>, <symmetry type>.
//   * grid: from <dimensions> and/or the header; if both are given they must agree (else an exception is thrown);
//     every grid point must be present in the file (checked).
//   * coordinate units: header "units" (default cm, the legacy behaviour); field unit: header <field unit>, or the
//     XML attribute field_unit (default tesla, i.e. legacy maps in gauss need scale="1e-4" and no header).
//     B = file value * unit * scale.
//   * symmetry (BrBz): header type cylindrical-z|x|y or XML attribute axis="z|x|y" gives the map axis in the hall
//     (GEMC: cylindrical-x maps "longitudinal" to hall x, e.g. the NH3 target field oxford_ptarget.dat).
//     The <shift> line of the header is NOT applied (gemc 2.9 ignores it too): use translationCoord.
//   * BxByBz with mirror_xy="true": map stored for x, y >= 0 only (first-quadrant TOSCA maps); (x, y) -> (|x|, |y|),
//     Bx -> sign(x) Bx, By -> sign(y) By.
// Position transformation: p_map = (T * R)^-1 p_hall, B_hall = R_field * R * B_map, R = R_user * R_axis.
// Outside the grid the map contributes nothing. fieldComponents() only reads shared data (thread safe).

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace fs = std::filesystem;

#include "FileLoaderHelper.h"

using namespace dd4hep;

namespace {

  struct GridAxis {
    double min{0}, max{0}, step{0};
    int    n{0};
    bool   valid() const { return n >= 2 && step > 0; }
  };

  double unit_value(std::string u, const std::string& what)
  {
    std::transform(u.begin(), u.end(), u.begin(), ::tolower);
    static const std::map<std::string, double> units = {
        {"mm", dd4hep::mm},       {"cm", dd4hep::cm},          {"m", dd4hep::m},
        {"gauss", dd4hep::gauss}, {"g", dd4hep::gauss},        {"kilogauss", dd4hep::kilogauss},
        {"kg", dd4hep::kilogauss}, {"tesla", dd4hep::tesla},   {"t", dd4hep::tesla},
        {"mt", 1e-3 * dd4hep::tesla}};
    auto it = units.find(u);
    if (it == units.end()) {
      throw std::runtime_error("FieldMapB: unknown " + what + " unit '" + u + "'");
    }
    return it->second;
  }

  std::string tag_attr(const std::string& tag, const std::string& attr)
  {
    std::smatch match;
    std::regex  re(attr + R"(\s*=\s*\"([^\"]*)\")");
    if (std::regex_search(tag, match, re))
      return match[1];
    return "";
  }

  // GEMC <mfield> header
  struct MapHeader {
    bool                     present{false};
    std::vector<GridAxis>    axes;       // first, second, third (in dd4hep units)
    std::vector<std::string> axis_names; // e.g. transverse, longitudinal
    std::string              field_unit;
    std::string              symmetry;
    std::string              shift;
  };

  MapHeader parse_header(const std::string& text)
  {
    MapHeader h;
    h.present = true;
    for (const char* which : {"first", "second", "third"}) {
      std::smatch match;
      std::regex  re(std::string("<") + which + R"(\b[^>]*>)");
      if (!std::regex_search(text, match, re))
        break;
      const std::string tag = match[0];
      GridAxis          a;
      a.n              = std::atoi(tag_attr(tag, "npoints").c_str());
      std::string un   = tag_attr(tag, "units");
      const double  cu = unit_value(un.empty() ? "cm" : un, "coordinate");
      a.min            = std::atof(tag_attr(tag, "min").c_str()) * cu;
      a.max            = std::atof(tag_attr(tag, "max").c_str()) * cu;
      a.step           = (a.n > 1) ? (a.max - a.min) / (a.n - 1) : 0.0;
      h.axes.push_back(a);
      h.axis_names.push_back(tag_attr(tag, "name"));
    }
    std::smatch match;
    if (std::regex_search(text, match, std::regex(R"(<field\b[^>]*>)")))
      h.field_unit = tag_attr(match[0], "unit");
    if (std::regex_search(text, match, std::regex(R"(<symmetry\b[^>]*>)")))
      h.symmetry = tag_attr(match[0], "type");
    if (std::regex_search(text, match, std::regex(R"(<shift\b[^>]*>)")))
      h.shift = match[0];
    return h;
  }

} // namespace

// implementation of the field map
class FieldMapB : public dd4hep::CartesianField::Object {
public:
  enum class Coord { BrBz, BxByBz };

  FieldMapB(const std::string& field_type_str, Coord coord) : m_coord(coord)
  {
    std::string ftype = field_type_str;
    std::transform(ftype.begin(), ftype.end(), ftype.begin(), ::tolower);
    if (ftype == "magnetic") {
      field_type = CartesianField::MAGNETIC;
    } else if (ftype == "electric") {
      field_type = CartesianField::ELECTRIC;
    } else {
      field_type = CartesianField::UNKNOWN;
      printout(ERROR, "FieldMapB", "Unknown field type " + ftype);
    }
  }

  // read the file: header (optional) + data; xml_axes may be empty (grid from header)
  void Load(const std::string& map_file, std::vector<GridAxis> xml_axes, const std::string& xml_field_unit,
            const std::string& xml_coord_unit, double scale, const std::string& xml_axis);

  void SetTransforms(const Transform3D& placement_user, const RotationZYX& field_rot)
  {
    m_place_user = placement_user;
    m_field_rot  = Rotation3D(field_rot);
  }
  void SetMirrorXY(bool on) { m_mirror_xy = on; }

  void fieldComponents(const double* pos, double* field) override;

  std::string Summary() const;

private:
  void finalize_transforms();

  Coord                 m_coord;
  std::vector<GridAxis> m_axes;
  std::vector<float>    m_B; // flattened [i0][i1]([i2])[component]
  int                   m_ncomp{2};
  bool                  m_mirror_xy{false};
  Transform3D           m_place_user;          // T * R_user (from XML)
  Rotation3D            m_axis_rot;            // map axis -> hall (symmetry)
  Rotation3D            m_field_rot;           // legacy rotationField
  Transform3D           m_hall_to_map;         // inverse of the full placement
  Rotation3D            m_map_to_hall_B;       // rotation of B vectors map -> hall
  bool                  m_identity{true};      // fast path: no transform at all
  bool                  m_translate_only{true};
  double                m_tx{0}, m_ty{0}, m_tz{0};
  std::string           m_info;
};

void FieldMapB::finalize_transforms()
{
  Transform3D placement = m_place_user * Transform3D(m_axis_rot);
  m_hall_to_map         = placement.Inverse();
  Rotation3D rot;
  placement.GetRotation(rot);
  m_map_to_hall_B = m_field_rot * rot;
  Rotation3D ident;
  m_translate_only = (rot == ident) && (m_field_rot == ident);
  const Translation3D t = placement.Translation();
  m_tx       = t.Vect().X();
  m_ty       = t.Vect().Y();
  m_tz       = t.Vect().Z();
  m_identity = m_translate_only && m_tx == 0 && m_ty == 0 && m_tz == 0;
}

void FieldMapB::Load(const std::string& map_file, std::vector<GridAxis> xml_axes, const std::string& xml_field_unit,
                     const std::string& xml_coord_unit, double scale, const std::string& xml_axis)
{
  const int ndim = (m_coord == Coord::BrBz) ? 2 : 3;
  m_ncomp        = ndim;

  std::ifstream input(map_file);
  if (!input) {
    throw std::runtime_error("FieldMapB: file " + map_file + " cannot be read");
  }

  // ---- header
  std::string line, header_text;
  MapHeader   hdr;
  std::streampos data_start = input.tellg();
  while (std::getline(input, line)) {
    auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) {
      data_start = input.tellg();
      continue;
    }
    if (line[first] == '<' || !header_text.empty()) {
      header_text += line + "\n";
      data_start = input.tellg();
      if (line.find("</mfield>") != std::string::npos)
        break;
      continue;
    }
    break; // first data line
  }
  if (!header_text.empty()) {
    hdr = parse_header(header_text);
  }
  input.clear();
  input.seekg(data_start);

  // ---- grid
  std::vector<GridAxis> axes;
  if (!xml_axes.empty()) {
    axes = xml_axes;
    if (hdr.present && static_cast<int>(hdr.axes.size()) == ndim) {
      for (int i = 0; i < ndim; ++i) {
        const auto& a = axes[i];
        const auto& b = hdr.axes[i];
        const double tol = 1e-3 * dd4hep::cm;
        if (std::fabs(a.min - b.min) > tol || std::fabs(a.max - b.max) > tol || a.n != b.n) {
          throw std::runtime_error(fmt::format(
              "FieldMapB: <dimensions> axis {} (min {} max {} n {}) does not match the header of {} (min {} max {} n {}) [cm]",
              i, a.min / dd4hep::cm, a.max / dd4hep::cm, a.n, map_file, b.min / dd4hep::cm, b.max / dd4hep::cm, b.n));
        }
      }
    }
  } else if (hdr.present && static_cast<int>(hdr.axes.size()) == ndim) {
    axes = hdr.axes;
  } else {
    throw std::runtime_error("FieldMapB: no grid for " + map_file + " (no <dimensions> and no usable <mfield> header)");
  }
  for (auto& a : axes) {
    if (!a.valid())
      throw std::runtime_error("FieldMapB: invalid grid axis (need step > 0 and >= 2 points)");
  }
  m_axes = axes;

  // ---- units
  std::string funit = hdr.field_unit;
  if (!xml_field_unit.empty()) {
    if (!funit.empty() && unit_value(funit, "field") != unit_value(xml_field_unit, "field")) {
      throw std::runtime_error("FieldMapB: field_unit=\"" + xml_field_unit + "\" contradicts the header unit \"" +
                               funit + "\" of " + map_file);
    }
    funit = xml_field_unit;
  }
  if (funit.empty())
    funit = "tesla"; // legacy: headerless maps, unit conversion through `scale`
  const double bunit = unit_value(funit, "field");
  if (hdr.present && !hdr.field_unit.empty() && std::fabs(scale) < 1e-2) {
    printout(WARNING, "FieldMapB",
             "scale = %g with a map whose header gives the field unit (%s): the unit conversion is automatic, "
             "scale is a dimensionless factor. Is scale=\"1e-4\" a leftover?",
             scale, hdr.field_unit.c_str());
  }
  double cunit = dd4hep::cm;
  if (!xml_coord_unit.empty())
    cunit = unit_value(xml_coord_unit, "coordinate");
  else if (hdr.present) {
    // header coordinate units were already applied to the grid; data rows use the same units
    std::smatch match;
    if (std::regex_search(header_text, match, std::regex(R"(<first\b[^>]*>)"))) {
      std::string un = tag_attr(match[0], "units");
      if (!un.empty())
        cunit = unit_value(un, "coordinate");
    }
  }

  // ---- map axis (2D symmetry)
  std::string axis = xml_axis;
  if (axis.empty() && !hdr.symmetry.empty()) {
    if (hdr.symmetry == "cylindrical-x")
      axis = "x";
    else if (hdr.symmetry == "cylindrical-y")
      axis = "y";
    else if (hdr.symmetry == "cylindrical-z")
      axis = "z";
  }
  if (axis.empty())
    axis = "z";
  if (m_coord == Coord::BrBz) {
    if (axis == "z")
      m_axis_rot = Rotation3D(); // identity
    else if (axis == "x")        // map z -> hall x, map x -> hall y, map y -> hall z
      m_axis_rot = Rotation3D(0, 0, 1, 1, 0, 0, 0, 1, 0);
    else if (axis == "y") // map z -> hall y, map x -> hall z, map y -> hall x
      m_axis_rot = Rotation3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
    else
      throw std::runtime_error("FieldMapB: axis must be x, y or z, got " + axis);
  }
  finalize_transforms();

  // ---- data
  std::size_t npts = 1;
  for (auto& a : m_axes)
    npts *= a.n;
  m_B.assign(npts * m_ncomp, 0.f);
  std::vector<char> filled(npts, 0);
  const double      bfac = bunit * scale;

  std::size_t nrows = 0, nout = 0, noffgrid = 0, ndup = 0, nbad = 0;
  double      v[6];
  while (std::getline(input, line)) {
    const char* p = line.c_str();
    int         k = 0;
    char*       end;
    for (; k < 2 * ndim; ++k) {
      v[k] = std::strtod(p, &end);
      if (end == p)
        break;
      p = end;
    }
    if (k == 0) {
      if (line.find_first_not_of(" \t\r") != std::string::npos)
        ++nbad;
      continue;
    }
    if (k < 2 * ndim) {
      ++nbad;
      continue;
    }
    ++nrows;
    std::size_t idx  = 0;
    bool        skip = false;
    for (int d = 0; d < ndim; ++d) {
      const auto&  a = m_axes[d];
      const double c = v[d] * cunit;
      const double u = (c - a.min) / a.step;
      const long   i = std::lround(u);
      if (i < 0 || i >= a.n) {
        skip = true;
        ++nout;
        break;
      }
      if (std::fabs(u - i) > 1e-3) {
        ++noffgrid;
        skip = true;
        break;
      }
      idx = idx * a.n + i;
    }
    if (skip)
      continue;
    if (filled[idx])
      ++ndup;
    filled[idx] = 1;
    for (int c = 0; c < m_ncomp; ++c)
      m_B[idx * m_ncomp + c] = static_cast<float>(v[ndim + c] * bfac);
  }
  const std::size_t nfilled = std::count(filled.begin(), filled.end(), 1);
  if (nbad > 0)
    printout(WARNING, "FieldMapB", "%s: %zu unparsable lines skipped", map_file.c_str(), nbad);
  if (nout > 0)
    printout(INFO, "FieldMapB", "%s: %zu rows outside the requested grid skipped", map_file.c_str(), nout);
  if (noffgrid > 0)
    printout(WARNING, "FieldMapB", "%s: %zu rows not on the grid nodes skipped", map_file.c_str(), noffgrid);
  if (ndup > 0)
    printout(WARNING, "FieldMapB", "%s: %zu duplicated grid points", map_file.c_str(), ndup);
  if (nfilled != npts) {
    throw std::runtime_error(fmt::format("FieldMapB: {}: only {} of {} grid points found in the file ({} rows read); "
                                         "check the <dimensions> / header grid",
                                         map_file, nfilled, npts, nrows));
  }

  std::ostringstream os;
  os << map_file << ": " << (m_coord == Coord::BrBz ? "BrBz" : "BxByBz") << " grid";
  for (auto& a : m_axes)
    os << " [" << a.min / dd4hep::cm << ", " << a.max / dd4hep::cm << "] cm x" << a.n;
  os << ", field unit " << funit << ", scale " << scale << ", axis " << axis << (hdr.present ? ", GEMC header" : "")
     << ", " << nfilled << " points";
  if (!hdr.shift.empty())
    os << " (header " << hdr.shift << " ignored, use translationCoord)";
  m_info = os.str();
}

std::string FieldMapB::Summary() const { return m_info; }

// get field components (thread safe: only local state is modified)
void FieldMapB::fieldComponents(const double* pos, double* field)
{
  double px = pos[0], py = pos[1], pz = pos[2];
  if (!m_identity) {
    if (m_translate_only) {
      px -= m_tx;
      py -= m_ty;
      pz -= m_tz;
    } else {
      auto p = m_hall_to_map * ROOT::Math::XYZPoint(px, py, pz);
      px     = p.x();
      py     = p.y();
      pz     = p.z();
    }
  }

  // locate a coordinate on an axis: index i (0..n-2) and fraction d; false if outside
  auto locate = [](const GridAxis& a, double c, int& i, double& d) {
    const double u = (c - a.min) / a.step;
    if (u < 0.0 || u > a.n - 1)
      return false;
    i = static_cast<int>(u);
    if (i >= a.n - 1)
      i = a.n - 2;
    d = u - i;
    return true;
  };

  double bx = 0, by = 0, bz = 0;
  if (m_coord == Coord::BrBz) {
    const double r = std::sqrt(px * px + py * py);
    int          ir, iz;
    double       dr, dz;
    if (!locate(m_axes[0], r, ir, dr) || !locate(m_axes[1], pz, iz, dz))
      return;
    const int    nz  = m_axes[1].n;
    const float* p00 = &m_B[(static_cast<std::size_t>(ir) * nz + iz) * 2];
    const float* p01 = p00 + 2;          // iz + 1
    const float* p10 = p00 + 2 * nz;     // ir + 1
    const float* p11 = p10 + 2;
    const double w00 = (1 - dr) * (1 - dz), w01 = (1 - dr) * dz, w10 = dr * (1 - dz), w11 = dr * dz;
    const double Br  = p00[0] * w00 + p01[0] * w01 + p10[0] * w10 + p11[0] * w11;
    const double Bz  = p00[1] * w00 + p01[1] * w01 + p10[1] * w10 + p11[1] * w11;
    if (r > 0) {
      bx = Br * px / r;
      by = Br * py / r;
    }
    bz = Bz;
  } else {
    double sx = 1, sy = 1;
    if (m_mirror_xy) {
      if (px < 0) {
        px = -px;
        sx = -1;
      }
      if (py < 0) {
        py = -py;
        sy = -1;
      }
    }
    int    ix, iy, iz;
    double dx, dy, dz;
    if (!locate(m_axes[0], px, ix, dx) || !locate(m_axes[1], py, iy, dy) || !locate(m_axes[2], pz, iz, dz))
      return;
    const std::size_t ny = m_axes[1].n, nz = m_axes[2].n;
    auto at = [&](int a, int b, int c, int comp) {
      return m_B[(((static_cast<std::size_t>(ix + a) * ny) + (iy + b)) * nz + (iz + c)) * 3 + comp];
    };
    double b[3];
    for (int comp = 0; comp < 3; ++comp) {
      const double b00 = at(0, 0, 0, comp) * (1 - dx) + at(1, 0, 0, comp) * dx;
      const double b01 = at(0, 0, 1, comp) * (1 - dx) + at(1, 0, 1, comp) * dx;
      const double b10 = at(0, 1, 0, comp) * (1 - dx) + at(1, 1, 0, comp) * dx;
      const double b11 = at(0, 1, 1, comp) * (1 - dx) + at(1, 1, 1, comp) * dx;
      const double b0  = b00 * (1 - dy) + b10 * dy;
      const double b1  = b01 * (1 - dy) + b11 * dy;
      b[comp]          = b0 * (1 - dz) + b1 * dz;
    }
    bx = sx * b[0];
    by = sy * b[1];
    bz = b[2];
  }

  if (!m_translate_only) {
    auto B = m_map_to_hall_B * ROOT::Math::XYZVector(bx, by, bz);
    bx     = B.x();
    by     = B.y();
    bz     = B.z();
  }
  field[0] += bx;
  field[1] += by;
  field[2] += bz;
}

// assign the field map to CartesianField
static Ref_t create_field_map_b(Detector& /*lcdd*/, xml::Handle_t handle)
{
  xml_comp_t x_par(handle);

  if (!x_par.hasAttr(_Unicode(field_map))) {
    throw std::runtime_error("FieldMapB Error: must have an xml attribute \"field_map\" for the field map.");
  }

  const std::string field_type = getAttrOrDefault<std::string>(x_par, _Unicode(field_type), "magnetic");
  const std::string coord_type = getAttrOrDefault<std::string>(x_par, _Unicode(coord_type), "BrBz");
  FieldMapB::Coord  coord;
  if (coord_type == "BrBz")
    coord = FieldMapB::Coord::BrBz;
  else if (coord_type == "BxByBz")
    coord = FieldMapB::Coord::BxByBz;
  else
    throw std::runtime_error("FieldMapB: coord_type " + coord_type + " is neither BrBz nor BxByBz");

  // grid from <dimensions> (optional when the file has a GEMC header)
  std::vector<GridAxis> axes;
  Transform3D           placement; // T * R_user
  RotationZYX           field_rot(0., 0., 0.);
  if (x_par.hasChild(_Unicode(dimensions))) {
    xml_comp_t x_dim = x_par.child(_Unicode(dimensions));
    std::vector<const char*> names =
        (coord == FieldMapB::Coord::BrBz) ? std::vector<const char*>{"R", "Z"} : std::vector<const char*>{"X", "Y", "Z"};
    int nfound = 0;
    for (auto axis_name : names) {
      if (x_dim.hasChild(xml::Strng_t(axis_name)))
        ++nfound;
    }
    if (nfound == static_cast<int>(names.size())) {
      for (auto axis_name : names) {
        xml_comp_t el = x_dim.child(xml::Strng_t(axis_name));
        GridAxis   a;
        a.step = el.step();
        a.min  = getAttrOrDefault<double>(el, _Unicode(min), 0.);
        a.max  = getAttrOrDefault<double>(el, _Unicode(max), 0.);
        a.n    = (a.step > 0) ? static_cast<int>(std::lround((a.max - a.min) / a.step)) + 1 : 0;
        axes.push_back(a);
      }
    } else if (nfound > 0) {
      throw std::runtime_error("FieldMapB: incomplete <dimensions> for coord_type " + coord_type);
    }
    Position t(0, 0, 0);
    if (x_dim.hasChild(_Unicode(translationCoord))) {
      xml_comp_t tr = x_dim.child(_Unicode(translationCoord));
      t = Position(getAttrOrDefault<double>(tr, _Unicode(x), 0.), getAttrOrDefault<double>(tr, _Unicode(y), 0.),
                   getAttrOrDefault<double>(tr, _Unicode(z), 0.));
    }
    RotationZYX rot_user(0., 0., 0.);
    if (x_dim.hasChild(_Unicode(rotationCoord))) {
      xml_comp_t rc = x_dim.child(_Unicode(rotationCoord));
      rot_user = RotationZYX(getAttrOrDefault<double>(rc, _Unicode(z), 0.), getAttrOrDefault<double>(rc, _Unicode(y), 0.),
                             getAttrOrDefault<double>(rc, _Unicode(x), 0.));
    }
    placement = Transform3D(Rotation3D(rot_user), t);
    if (x_dim.hasChild(_Unicode(rotationField))) { // legacy: bare numbers in degrees
      xml_comp_t   rf    = x_dim.child(_Unicode(rotationField));
      const double deg2r = M_PI / 180.;
      field_rot = RotationZYX(getAttrOrDefault<double>(rf, _Unicode(z), 0.) * deg2r,
                              getAttrOrDefault<double>(rf, _Unicode(y), 0.) * deg2r,
                              getAttrOrDefault<double>(rf, _Unicode(x), 0.) * deg2r);
    }
  }

  const std::string field_map_file  = x_par.attr<std::string>(_Unicode(field_map));
  const std::string field_map_url   = getAttrOrDefault<std::string>(x_par, _Unicode(url), "");
  const std::string field_map_cache = getAttrOrDefault<std::string>(x_par, _Unicode(cache), "");
  const std::string path = EnsureFileFromURLExists(field_map_url, field_map_file, field_map_cache, false);

  const double      scale      = getAttrOrDefault<double>(x_par, _Unicode(scale), 1.0);
  const std::string field_unit = getAttrOrDefault<std::string>(x_par, _Unicode(field_unit), "");
  const std::string coord_unit = getAttrOrDefault<std::string>(x_par, _Unicode(coord_unit), "");
  const std::string axis       = getAttrOrDefault<std::string>(x_par, _Unicode(axis), "");

  auto map = new FieldMapB(field_type, coord);
  map->SetTransforms(placement, field_rot);
  map->SetMirrorXY(getAttrOrDefault<bool>(x_par, _Unicode(mirror_xy), false));
  map->Load(path, axes, field_unit, coord_unit, scale, axis);
  printout(INFO, "FieldMapB", "%s: %s", x_par.nameStr().c_str(), map->Summary().c_str());

  CartesianField field;
  field.assign(map, x_par.nameStr(), "FieldMapB");
  return field;
}

DECLARE_XMLELEMENT(SoLID_FieldMapB, create_field_map_b)
