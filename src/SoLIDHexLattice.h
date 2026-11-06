// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 SoLID collaboration
//
// Framework-neutral helpers for the SoLID EC module layout:
//   * the ANL 2013 hexagonal lattice ("SOLID Tube Coordinates", T. O'Connor 2013-06-28; GEMC
//     id<i>_{x,y,sector} block of every ec_segmented*/ *__parameters_Original.txt), encoded as a row table;
//   * module selection by centre radius (strict Rmin < R < Rmax, as GEMC make_ec_module);
//   * an optional azimuthal mask (N segments of width dphi, explicit or evenly spaced centres) for SoLID-Lite
//     ("4seg24degFA");
//   * the ANL sector rule (sector = floor(phi/12 deg) + 1, phi = atan2(y,x) in [0, 360) deg).
// No DD4hep dependency: the same code can drive a future block-array (sPHENIX/STAR FCS) EC plugin or a
// reconstruction-side module map. All lengths/angles are in whatever unit system the caller uses
// (the DD4hep plugin passes DD4hep internal units).
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace solid::ec {

/// One lattice position. `id` is the ANL index (1-based, stable across configurations).
struct LatticeSite {
  int id = 0;    ///< ANL module id (1..N)
  int row = 0;   ///< ANL row index k (1..54), "idy"
  int col = 0;   ///< position in row, 1-based ("idx")
  double x = 0;  ///< centre x
  double y = 0;  ///< centre y
  int sector = 0;
};

/// A run of column indices n = a, a+2, ..., b in one row (x = n * pitch / 2).
struct RowRun { int a = 0, b = 0; };
struct Row { int k = 0; std::vector<RowRun> runs; };

/// Parse "k:a:b[,a:b...];k:a:b;..." (whitespace ignored).
inline std::vector<Row> parse_rows(const std::string& spec) {
  std::string s;
  for (char c : spec) if (!std::isspace(static_cast<unsigned char>(c))) s += c;
  std::vector<Row> rows;
  std::stringstream ss(s);
  std::string rowtok;
  while (std::getline(ss, rowtok, ';')) {
    if (rowtok.empty()) continue;
    auto p = rowtok.find(':');
    if (p == std::string::npos) throw std::runtime_error("SoLIDHexLattice: bad row token '" + rowtok + "'");
    Row r;
    r.k = std::stoi(rowtok.substr(0, p));
    std::stringstream rs(rowtok.substr(p + 1));
    std::string runtok;
    while (std::getline(rs, runtok, ',')) {
      auto q = runtok.find(':');
      if (q == std::string::npos) throw std::runtime_error("SoLIDHexLattice: bad run token '" + runtok + "'");
      RowRun run{std::stoi(runtok.substr(0, q)), std::stoi(runtok.substr(q + 1))};
      if (run.b < run.a || ((run.b - run.a) % 2) != 0)
        throw std::runtime_error("SoLIDHexLattice: run must be a..b in steps of 2: '" + runtok + "'");
      r.runs.push_back(run);
    }
    rows.push_back(r);
  }
  return rows;
}

inline double phi_0_2pi(double x, double y) {
  double phi = std::atan2(y, x);
  if (phi < 0) phi += 2. * M_PI;
  return phi;
}

/// Build all sites of the row table: x = n * pitch/2, y = (k - row_center) * row_pitch,
/// ids increase along +x within a row, rows from k = 1 (most negative y) upwards.
inline std::vector<LatticeSite> build_sites(const std::vector<Row>& rows, double pitch, double row_pitch,
                                            double row_center, double sector_dphi) {
  std::vector<LatticeSite> out;
  int id = 0;
  for (const auto& r : rows) {
    int col = 0;
    for (const auto& run : r.runs) {
      for (int n = run.a; n <= run.b; n += 2) {
        LatticeSite s;
        s.id = ++id;
        s.row = r.k;
        s.col = ++col;
        s.x = n * pitch / 2.;
        s.y = (r.k - row_center) * row_pitch;
        s.sector = static_cast<int>(std::floor(phi_0_2pi(s.x, s.y) / sector_dphi + 1e-12)) + 1;
        out.push_back(s);
      }
    }
  }
  return out;
}

/// Vertices of a regular polygon (apothem `a`, first vertex at angle `phi0`) centred at (x, y).
inline std::vector<std::array<double, 2>> polygon(double x, double y, int nsides, double apothem, double phi0) {
  std::vector<std::array<double, 2>> v;
  const double rc = apothem / std::cos(M_PI / nsides);
  for (int i = 0; i < nsides; ++i) {
    double p = phi0 + 2. * M_PI * i / nsides;
    v.push_back({x + rc * std::cos(p), y + rc * std::sin(p)});
  }
  return v;
}

/// Minimum and maximum distance of a convex polygon from the origin (origin assumed outside the polygon).
inline std::pair<double, double> polygon_rminmax(const std::vector<std::array<double, 2>>& v) {
  double rmin = 1e300, rmax = 0;
  const size_t n = v.size();
  for (size_t i = 0; i < n; ++i) {
    const auto& p = v[i];
    const auto& q = v[(i + 1) % n];
    rmax = std::max(rmax, std::hypot(p[0], p[1]));
    double dx = q[0] - p[0], dy = q[1] - p[1];
    double t = -(p[0] * dx + p[1] * dy) / (dx * dx + dy * dy);
    t = std::clamp(t, 0., 1.);
    rmin = std::min(rmin, std::hypot(p[0] + t * dx, p[1] + t * dy));
  }
  return {rmin, rmax};
}

/// Optional azimuthal mask: keep a module when its centre (mode Centre) or any part of its footprint
/// (mode Overlap) lies inside one of the windows [c_i - dphi/2, c_i + dphi/2].
struct PhiMask {
  enum class Mode { Centre = 0, Overlap = 1 };
  bool enabled = false;
  Mode mode = Mode::Centre;
  double dphi = 0;
  std::vector<double> centres;

  static double wrap_pm_pi(double a) {
    a = std::fmod(a + M_PI, 2. * M_PI);
    if (a < 0) a += 2. * M_PI;
    return a - M_PI;
  }
  bool phi_inside(double phi) const {
    for (double c : centres)
      if (std::abs(wrap_pm_pi(phi - c)) <= dphi / 2.) return true;
    return false;
  }
  /// does the line through the origin at angle `e` cross the (far-away, convex) polygon?
  static bool edge_crosses(double e, const std::vector<std::array<double, 2>>& v) {
    const double ux = std::cos(e), uy = std::sin(e);
    bool pos = false, neg = false, front = false;
    for (const auto& p : v) {
      double cr = ux * p[1] - uy * p[0];
      (cr >= 0 ? pos : neg) = true;
      if (ux * p[0] + uy * p[1] > 0) front = true;
    }
    return pos && neg && front;
  }
  bool accept(double x, double y, const std::vector<std::array<double, 2>>& footprint) const {
    if (!enabled) return true;
    if (phi_inside(phi_0_2pi(x, y))) return true;
    if (mode == Mode::Centre) return false;
    for (const auto& p : footprint)
      if (phi_inside(phi_0_2pi(p[0], p[1]))) return true;
    for (double c : centres)
      if (edge_crosses(c - dphi / 2., footprint) || edge_crosses(c + dphi / 2., footprint)) return true;
    return false;
  }
};

} // namespace solid::ec
