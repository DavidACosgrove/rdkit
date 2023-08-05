//
//  Copyright (C) 2023 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <chrono>
#include <random>
#include <vector>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#define CATCH_CONFIG_MAIN
#include "../../../External/catch/catch/single_include/catch2/catch.hpp"

#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
#include <GraphMol/RascalMCES/RascalDetails.h>

TEST_CASE("Small test", "[basics]") {
  std::string fName = getenv("RDBASE");
  fName += "/Contrib/Fastcluster/cdk2.smi";
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  std::cout << "Read " << mols.size() << " mols" << std::endl;
  RDKit::RascalMCES::RascalOptions opts;
  opts.similarityThreshold = 0.7;
  RDKit::RascalMCES::rascalCluster(mols, opts);
}

TEST_CASE("Medium test", "[basics]") {
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/RascalMCES/BLSets_selected_actives_0.05.smi";
  std::cout << fName << std::endl;
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  std::cout << "Read " << mols.size() << " mols" << std::endl;
  RDKit::RascalMCES::RascalOptions opts;
  opts.similarityThreshold = 0.7;
  RDKit::RascalMCES::rascalCluster(mols, opts);
}

TEST_CASE("Monster test", "[basics]") {
  std::string fName = "/Users/david/Projects/Moonshot/activity_data.csv";
  RDKit::SmilesMolSupplier suppl(fName, ",", 0, 1, true);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  std::cout << "Read " << mols.size() << " mols" << std::endl;
  RDKit::RascalMCES::RascalOptions opts;
  opts.similarityThreshold = 0.7;
  RDKit::RascalMCES::rascalCluster(mols, opts);
}
