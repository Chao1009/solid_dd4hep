// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sakib Rahman (ensureExtension / xmlToProtoSurfaceMaterial, from ePIC DD4hepDetectorHelper.h)
// Copyright (C) 2026 SoLID collaboration
//
// Small helpers shared by the SoLID tracker plugins:
//   * ensureExtension<T>()            attach (or fetch) a DDRec extension such as VariantParameters;
//   * xmlToProtoSurfaceMaterial()     translate <layer_material surface=... binning=... bins0=... bins1=.../>
//                                     into the VariantParameters keys read by the ACTS DD4hep plugin
//                                     (layer_material, layer_material_<surface>, layer_material_<surface>_<bin>);
//   * parseIndexList()                "1-30", "4,5,11,12", "1-3,7,20-22" -> sorted unique list of integers.
#pragma once

#include "DD4hep/DetElement.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace SoLIDTrackerHelpers {

template <typename T> T& ensureExtension(dd4hep::DetElement& elt) {
  T* ext = elt.extension<T>(false);
  if (ext == nullptr) {
    ext = new T();
    elt.addExtension<T>(ext);
  }
  return *ext;
}

/// <layer_material surface="inner|outer|representing" binning="binPhi,binR" bins0="30" bins1="20"/>
inline void xmlToProtoSurfaceMaterial(const xml_comp_t& x_material, dd4hep::rec::VariantParameters& params,
                                      const std::string& baseTag) {
  params.set<bool>(baseTag, true);
  const std::string surface = x_material.attr<std::string>(_Unicode(surface));
  const std::string binning = x_material.attr<std::string>(_Unicode(binning));
  std::vector<std::string> bins;
  std::stringstream ss(binning);
  for (std::string tok; std::getline(ss, tok, ',');) {
    tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
    if (!tok.empty())
      bins.push_back(tok);
  }
  if (bins.size() != 2)
    throw std::runtime_error("layer_material: binning must name two bin values, got '" + binning + "'");
  const std::string tag = baseTag + "_" + surface;
  params.set<bool>(tag, true);
  params.set<int>(tag + "_" + bins[0], x_material.attr<int>(_Unicode(bins0)));
  params.set<int>(tag + "_" + bins[1], x_material.attr<int>(_Unicode(bins1)));
}

/// Parse "1-30", "4,5,11,12", "1-3, 7, 20-22" (whitespace ignored) into a sorted unique list.
inline std::vector<int> parseIndexList(const std::string& text) {
  std::vector<int> out;
  std::stringstream ss(text);
  for (std::string tok; std::getline(ss, tok, ',');) {
    tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
    if (tok.empty())
      continue;
    const auto dash = tok.find('-', 1);
    if (dash == std::string::npos) {
      out.push_back(std::stoi(tok));
    } else {
      const int a = std::stoi(tok.substr(0, dash));
      const int b = std::stoi(tok.substr(dash + 1));
      for (int i = std::min(a, b); i <= std::max(a, b); ++i)
        out.push_back(i);
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

} // namespace SoLIDTrackerHelpers
