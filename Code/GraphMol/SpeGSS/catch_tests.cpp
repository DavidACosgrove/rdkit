//
// Copyright (C) David Cosgrove 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <cstdio>
#include <filesystem>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SpeGSS/SpeGSSSurface.H>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SpeGSS;

TEST_CASE("Generate Surface") {
  {
    const auto ethane =
        "CC |(0.8151,-0.5383,0.4928;2.1548,-0.0516,-0.0117)|"_smiles;
    REQUIRE(ethane);
    SpeGSSParams params{.contourValue = 0.5};
    auto surf = std::make_unique<SpeGSSSurface>(*ethane, params);
    REQUIRE(surf);
    CHECK(surf->getNumVertices() == 980);
    CHECK(surf->getNumTriangles() == 1956);
    CHECK(surf->isManifold());
  }
}

TEST_CASE("Write and read surfaces") {
  const auto benzene =
      "c1ccccc1 |(0.0737885,1.41082,0.0189081;1.21777,0.629842,0.0620298;1.1384,-0.751516,0.0431948;-0.0705508,-1.41453,-0.0187994;-1.20555,-0.63279,-0.061502;-1.1425,0.763757,-0.0432364)|"_smiles;
  REQUIRE(benzene);
  SpeGSSParams params{.contourValue = 0.5};
  auto surf = std::make_unique<SpeGSSSurface>(*benzene, params);
  REQUIRE(surf);
  CHECK(surf->getNumVertices() == 1816);
  CHECK(surf->getNumTriangles() == 3628);
  CHECK(surf->isManifold());

  std::string tmpFile = std::tmpnam(nullptr);
  surf->writeFile(tmpFile);
  std::cout << tmpFile << std::endl;

  auto newSurf = std::make_unique<SpeGSSSurface>();
  REQUIRE(newSurf);
  newSurf->readFile(tmpFile);
  CHECK(newSurf->getNumVertices() == 1816);
  CHECK(newSurf->getNumTriangles() == 3628);
  CHECK(newSurf->isManifold());
  CHECK_THAT(newSurf->getVertices().front(),
             Catch::Matchers::WithinAbs(surf->getVertices().front(), 1e-5));
  CHECK_THAT(newSurf->getVertices().back(),
             Catch::Matchers::WithinAbs(surf->getVertices().back(), 1e-5));
  CHECK(newSurf->getTriangles().front() == surf->getTriangles().front());
  CHECK(newSurf->getTriangles().back() == surf->getTriangles().back());
  CHECK_THAT(newSurf->getNormals().front(),
             Catch::Matchers::WithinAbs(surf->getNormals().front(), 1e-5));
  CHECK_THAT(newSurf->getNormals().back(),
             Catch::Matchers::WithinAbs(surf->getNormals().back(), 1e-5));

  std::filesystem::remove(tmpFile);
}
