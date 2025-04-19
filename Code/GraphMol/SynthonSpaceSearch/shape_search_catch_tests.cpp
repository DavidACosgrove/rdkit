//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <cstdio>

#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;
using namespace RDKit::RascalMCES;

const char *rdbase = getenv("RDBASE");

TEST_CASE("Shape Small tests") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::vector<std::string> libNames{
      fullRoot + "amide_space.txt",
      fullRoot + "triazole_space.txt",
      fullRoot + "urea_space.txt",
  };
  std::vector<std::string> enumLibNames{
      fullRoot + "amide_space_enum.smi",
      fullRoot + "triazole_space_enum.smi",
      fullRoot + "urea_space_enum.smi",
  };

  std::vector<std::string> querySmis{
      "c1ccccc1C(=O)N1CCCC1",
      "CC1CCN(c2nnc(CO)n2C2CCCC2)C1",
      "C[C@@H]1CC(NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1",
  };

  std::vector<size_t> expNumHits{6, 4, 1};

  for (size_t i = 0; i < libNames.size(); i++) {
    if (i != 0) {
      continue;
    }
    SynthonSpace synthonspace;
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    synthonspace.buildSynthonConformers();

    SynthonSpaceSearchParams params;
    params.similarityCutoff = 1.4;
    params.numConformers = 100;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.shapeSearch(*queryMol, params);
    std::cout << "Num hits : " << results.getHitMolecules().size() << " : "
              << results.getMaxNumResults() << std::endl;
  }
}

TEST_CASE("Shape DB Writer") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumReactions() == 1);
  synthonspace.buildSynthonConformers();

  auto spaceName = std::tmpnam(nullptr);

  synthonspace.writeDBFile(spaceName);

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile(spaceName);
  CHECK(newsynthonspace.getNumReactions() == 1);
  std::shared_ptr<SynthonSet> irxn;
  CHECK_NOTHROW(irxn = newsynthonspace.getReaction("doebner-miller-quinoline"));

  const auto &orxn = synthonspace.getReaction("doebner-miller-quinoline");
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    REQUIRE(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons().size(); ++j) {
      REQUIRE(irxn->getSynthons()[i][j].second->getShapes().size() ==
              orxn->getSynthons()[i][j].second->getShapes().size());
      for (size_t k = 0;
           k < irxn->getSynthons()[i][j].second->getShapes().size(); ++k) {
        const auto ishape =
            irxn->getSynthons()[i][j].second->getShapes()[k].get();
        const auto oshape =
            orxn->getSynthons()[i][j].second->getShapes()[k].get();
        CHECK(ishape->sov == Catch::Approx(oshape->sov));
      }
    }
  }
}