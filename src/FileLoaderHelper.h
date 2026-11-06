// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck; 2026 SoLID collaboration (cache by file name, ${VAR}, no symlink mode)
//
// EnsureFileFromURLExists(url, file, cache, link): make sure a data file (e.g. a field map) is available and
// return the path to read it from. Resolution order:
//   1. `file` itself (environment variables expanded) if it exists and is non-empty (e.g. a manual copy);
//   2. every directory of the colon-separated `cache` list (environment variables expanded, e.g.
//      "${DETECTOR_PATH}/fieldmaps"): a file named like `file` (basename) or named by the 16-hex-digit hash of the
//      URL (ePIC convention), also searched recursively below the cache directory;
//   3. download `url` with curl into the FIRST WRITABLE cache directory (as <basename>, atomically via a temporary
//      file + rename; the URL is recorded in <basename>.url), or, if no cache directory is writable, next to `file`.
// With link = true (FileLoader plugin, ePIC behaviour) a symlink `file` -> resolved path is created when they differ
// and `file` is returned; with link = false the resolved path is returned and nothing is created next to `file`.
// Offline use: copy the file (same basename) into the cache directory, or to the `file` path itself.
// Errors throw std::runtime_error.

#pragma once

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Primitives.h>
#include <DD4hep/Printout.h>

#include <fmt/core.h>

#include <unistd.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

using dd4hep::printout;
using dd4hep::ERROR, dd4hep::WARNING, dd4hep::INFO, dd4hep::DEBUG;

namespace FileLoaderHelper {
  static constexpr const char* const kCommand = "curl --silent --show-error --retry 5 --location --fail {0} --output {1}";

  // expand $VAR and ${VAR} (undefined variables expand to "")
  inline std::string ExpandEnv(std::string s)
  {
    std::string::size_type pos = 0;
    while ((pos = s.find('$', pos)) != std::string::npos) {
      std::string::size_type name_begin = pos + 1, name_end, erase_end;
      if (name_begin < s.size() && s[name_begin] == '{') {
        name_end = s.find('}', name_begin + 1);
        if (name_end == std::string::npos)
          break;
        erase_end  = name_end + 1;
        name_begin = name_begin + 1;
      } else {
        name_end  = s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_", name_begin);
        if (name_end == std::string::npos)
          name_end = s.size();
        erase_end = name_end;
      }
      const std::string name(s.substr(name_begin, name_end - name_begin));
      const char*       env = std::getenv(name.c_str());
      const std::string value(env != nullptr ? env : "");
      if (env == nullptr)
        printout(WARNING, "FileLoader", "environment variable $" + name + " is not set");
      s.replace(pos, erase_end - pos, value);
      pos += value.size();
    }
    return s;
  }

  inline std::vector<std::string> SplitPaths(const std::string& s)
  {
    std::vector<std::string> out;
    std::string::size_type   start = 0;
    while (start <= s.size()) {
      auto end = s.find(':', start);
      if (end == std::string::npos)
        end = s.size();
      if (end > start)
        out.emplace_back(s.substr(start, end - start));
      start = end + 1;
    }
    return out;
  }

  inline bool IsUsableFile(const fs::path& p)
  {
    std::error_code ec;
    return fs::exists(p, ec) && fs::is_regular_file(fs::status(p, ec)) && fs::file_size(p, ec) > 0 && !ec;
  }

  inline bool IsWritableDir(const fs::path& d)
  {
    std::error_code ec;
    if (!fs::exists(d, ec)) {
      fs::create_directories(d, ec);
      if (ec)
        return false;
    }
    return fs::is_directory(d, ec) && ::access(d.c_str(), W_OK) == 0;
  }

  // look for <name> or <hash> in dir, then recursively below it
  inline fs::path FindInCache(const fs::path& dir, const std::string& name, const std::string& hash)
  {
    std::error_code ec;
    if (!fs::is_directory(dir, ec))
      return {};
    for (const auto& cand : {dir / name, dir / hash}) {
      if (IsUsableFile(cand))
        return cand;
    }
    for (auto it = fs::recursive_directory_iterator(dir, fs::directory_options::skip_permission_denied, ec);
         !ec && it != fs::recursive_directory_iterator(); it.increment(ec)) {
      if (!it->is_directory(ec))
        continue;
      for (const auto& cand : {it->path() / name, it->path() / hash}) {
        if (IsUsableFile(cand))
          return cand;
      }
    }
    return {};
  }

  inline void CheckRecordedURL(const fs::path& found, const std::string& url)
  {
    fs::path      rec(found.string() + ".url");
    std::ifstream in(rec);
    std::string   recorded;
    if (in && std::getline(in, recorded) && !recorded.empty() && recorded != url) {
      printout(WARNING, "FileLoader", "%s was downloaded from %s, the compact file asks for %s; delete it to refetch",
               found.c_str(), recorded.c_str(), url.c_str());
    }
  }

  inline void MakeLink(const fs::path& file_path, const fs::path& target)
  {
    std::error_code ec;
    if (!file_path.parent_path().empty())
      fs::create_directories(file_path.parent_path(), ec);
    if (fs::is_symlink(fs::symlink_status(file_path, ec)))
      fs::remove(file_path, ec);
    fs::path link_target = target.is_absolute() ? target : fs::proximate(target, file_path.parent_path());
    fs::create_symlink(link_target, file_path, ec);
    if (ec) {
      throw std::runtime_error("FileLoader: unable to link " + file_path.string() + " -> " + link_target.string() +
                               ": " + ec.message());
    }
  }
} // namespace FileLoaderHelper

// Function to make sure a file is available; returns the path to read it from
inline std::string EnsureFileFromURLExists(std::string url, std::string file, std::string cache_str = "",
                                           bool link = true)
{
  using namespace FileLoaderHelper;
  file      = ExpandEnv(file);
  cache_str = ExpandEnv(cache_str);
  url       = ExpandEnv(url);
  if (file.empty()) {
    throw std::runtime_error("FileLoader: no file specified");
  }

  const fs::path    file_path(file);
  const std::string name = file_path.filename().string();
  const std::string hash = fmt::format("{:016x}", dd4hep::detail::hash64(url));

  // 1. the file itself (manual copy, or a valid link from a previous run)
  if (IsUsableFile(file_path)) {
    printout(INFO, "FileLoader", "using " + file);
    return file;
  }
  std::error_code ec;
  if (fs::is_symlink(fs::symlink_status(file_path, ec))) {
    printout(INFO, "FileLoader", "removing broken symlink " + file);
    fs::remove(file_path, ec);
  }

  // 2. cache directories
  const auto caches = SplitPaths(cache_str);
  for (const auto& c : caches) {
    fs::path found = FindInCache(fs::path(c), name, hash);
    if (!found.empty()) {
      printout(INFO, "FileLoader", "found " + name + " in cache: " + found.string());
      CheckRecordedURL(found, url);
      if (link) {
        MakeLink(file_path, found);
        return file;
      }
      return found.string();
    }
  }

  // 3. download into the first writable cache directory (or next to `file`)
  if (url.empty()) {
    throw std::runtime_error("FileLoader: " + file + " not found (cache: '" + cache_str + "') and no url given");
  }
  fs::path dest_dir;
  for (const auto& c : caches) {
    if (IsWritableDir(fs::path(c))) {
      dest_dir = fs::path(c);
      break;
    }
  }
  if (dest_dir.empty()) {
    dest_dir = file_path.parent_path().empty() ? fs::path(".") : file_path.parent_path();
    if (!IsWritableDir(dest_dir)) {
      throw std::runtime_error("FileLoader: no writable directory to download " + name + " (cache: '" + cache_str +
                               "', file: " + file + "); download " + url + " manually into the cache directory");
    }
  }
  const fs::path dest = dest_dir / name;
  const fs::path tmp  = dest_dir / fmt::format(".{}.download.{}", name, ::getpid());
  const std::string cmd = fmt::format(kCommand, url, tmp.string());
  printout(INFO, "FileLoader", "downloading " + url + " -> " + dest.string());
  const int ret = std::system(cmd.c_str());
  if (ret != 0 || !IsUsableFile(tmp)) {
    fs::remove(tmp, ec);
    throw std::runtime_error("FileLoader: download failed (" + cmd + ", status " + std::to_string(ret) +
                             "). Offline? Copy " + name + " into " + dest_dir.string());
  }
  fs::rename(tmp, dest, ec);
  if (ec) {
    throw std::runtime_error("FileLoader: cannot rename " + tmp.string() + " to " + dest.string());
  }
  {
    std::ofstream rec(dest.string() + ".url");
    rec << url << "\n";
  }
  if (link && fs::absolute(dest) != fs::absolute(file_path)) {
    MakeLink(file_path, dest);
    return file;
  }
  return dest.string();
}
