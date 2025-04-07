//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

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
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
}
}