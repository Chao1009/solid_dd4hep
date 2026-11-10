// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck; 2026 SoLID collaboration

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Primitives.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <fmt/core.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "FileLoaderHelper.h"

using namespace dd4hep;

void usage(int argc, char** argv)
{
  std::cerr << "Usage: -plugin <name> -arg [-arg]                                                  \n"
               "     cache:<string>           cache location (may be read-only)                    \n"
               "     file:<string>            file location                                        \n"
               "     url:<string>             url location                                         \n"
               "     sha256:<string>          optional SHA-256 (64 hex digits) the file must have  \n"
               "\tArguments given: "
            << arguments(argc, argv) << std::endl;
  std::exit(EINVAL);
}

// Plugin to download files (see FileLoaderHelper.h for the cache / download logic).
// Example: <plugin name="epic_FileLoader"><arg value="file:fieldmaps/x.dat"/><arg value="url:https://..."/>
//          <arg value="cache:${DETECTOR_PATH}/fieldmaps"/></plugin>
long load_file(Detector& /* desc */, int argc, char** argv)
{
  // argument parsing
  std::string cache, file, url, sha256;
  for (int i = 0; i < argc && argv[i]; ++i) {
    if (0 == std::strncmp("cache:", argv[i], 6))
      cache = (argv[i] + 6);
    else if (0 == std::strncmp("file:", argv[i], 5))
      file = (argv[i] + 5);
    else if (0 == std::strncmp("url:", argv[i], 4))
      url = (argv[i] + 4);
    else if (0 == std::strncmp("sha256:", argv[i], 7))
      sha256 = (argv[i] + 7);
    else {
      std::cerr << "Unexpected argument \"" << argv[i] << "\"" << std::endl;
      usage(argc, argv);
    }
  }
  printout(DEBUG, "FileLoader", "arg cache: " + cache);
  printout(DEBUG, "FileLoader", "arg file: " + file);
  printout(DEBUG, "FileLoader", "arg url: " + url);

  // if file or url is empty, do nothing
  if (file.empty()) {
    printout(WARNING, "FileLoader", "no file specified");
  }
  if (url.empty()) {
    printout(WARNING, "FileLoader", "no url specified");
  }

  const std::string path = EnsureFileFromURLExists(url, file, cache, true, sha256);
  printout(INFO, "FileLoader", "file " + file + " available as " + path);

  return 1;
}

DECLARE_APPLY(epic_FileLoader, load_file)
DECLARE_APPLY(SoLID_FileLoader, load_file)
