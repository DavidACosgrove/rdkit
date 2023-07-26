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

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#define CATCH_CONFIG_MAIN
#include "../../../External/catch/catch/single_include/catch2/catch.hpp"

#include <GraphMol/RascalMCES/rascal_mces.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

// 'The paper' referenced below is :
// RASCAL: Calculation of Graph Similarity using Maximum Common
// Edge Subgraphs, John W. Raymond, Eleanor J. Gardiner, Peter Willett
// 'The Computer Journal', 45, 631-644 (2002).
// https://eprints.whiterose.ac.uk/3568/1/willets3.pdf

namespace RDKit {
namespace RascalMCES {
double tier_1_sim(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
                  std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                  std::map<int, std::vector<std::pair<int, int>>> &degSeqs2);

double tier2Sim(const ROMol &mol1, const ROMol &mol2,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs2,
                const std::vector<unsigned int> &bondLabels1,
                const std::vector<unsigned int> &bondLabels2);

void get_bond_labels(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
                     const RascalOptions &opts,
                     std::vector<unsigned int> &bondLabels1,
                     std::vector<unsigned int> &bondLabels2);
}  // namespace RascalMCES
}  // namespace RDKit

using namespace RDKit;
using namespace RDKit::RascalMCES;

void check_smarts_ok(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
                     const RascalResult &res) {
  std::unique_ptr<RDKit::ROMol> qmol(RDKit::SmartsToMol(res.smarts()));
  RDKit::MatchVectType dont_care;
  REQUIRE(RDKit::SubstructMatch(mol1, *qmol, dont_care));
  REQUIRE(RDKit::SubstructMatch(mol2, *qmol, dont_care));
}

void check_smarts_found(const RascalResult &res,
                        const std::vector<std::string> &equiv_smarts) {
  auto found_it =
      std::find(equiv_smarts.begin(), equiv_smarts.end(), res.smarts());
  REQUIRE(found_it != equiv_smarts.end());
}

void check_expected_bonds(
    const RascalResult &res,
    const std::vector<std::vector<std::pair<int, int>>> &exp_bond_matches) {
  auto found_it = std::find(exp_bond_matches.begin(), exp_bond_matches.end(),
                            res.bondMatches());
  REQUIRE(found_it != exp_bond_matches.end());
}

TEST_CASE("Very small test", "[basics]") {
  auto m1 = "OCC(=O)N"_smiles;
  REQUIRE(m1);
  auto m2 = "NC(=O)C=O"_smiles;
  REQUIRE(m2);

  RascalOptions opts;
  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  auto tier1_sim = tier_1_sim(*m1, *m2, degSeqs1, degSeqs2);
  REQUIRE_THAT(tier1_sim, Catch::Matchers::WithinAbs(1.0000, 0.0001));

  std::vector<unsigned int> bondLabels1, bondLabels2;
  get_bond_labels(*m1, *m2, opts, bondLabels1, bondLabels2);
  auto tier2_sim =
      tier2Sim(*m1, *m2, degSeqs1, degSeqs2, bondLabels1, bondLabels2);
  REQUIRE_THAT(tier2_sim, Catch::Matchers::WithinAbs(0.7901, 0.0001));

  auto res = rascalMces(*m1, *m2, opts);

  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 3);
  std::vector<std::pair<int, int>> exp_bond_matches{{1, 2}, {2, 1}, {3, 0}};
  REQUIRE(exp_bond_matches == res.front().bondMatches());
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.6049, 0.0001));

  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("Default options", "[basics]") {
  auto m1 = "OCC(=O)N"_smiles;
  REQUIRE(m1);
  auto m2 = "NC(=O)C=O"_smiles;
  REQUIRE(m2);

  auto res = rascalMces(*m1, *m2);

  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 3);
  std::vector<std::pair<int, int>> exp_bond_matches{{1, 2}, {2, 1}, {3, 0}};
  REQUIRE(exp_bond_matches == res.front().bondMatches());
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.6049, 0.0001));

  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("Juglone vs Scopoletin test", "[basics]") {
  auto m1 = "Oc1cccc2C(=O)C=CC(=O)c12"_smiles;
  REQUIRE(m1);
  auto m2 = "O1C(=O)C=Cc2cc(OC)c(O)cc12"_smiles;
  REQUIRE(m2);

  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  auto tier1_sim = tier_1_sim(*m1, *m2, degSeqs1, degSeqs2);
  REQUIRE_THAT(tier1_sim, Catch::Matchers::WithinAbs(0.8633, 0.0001));

  RascalOptions opts;
  std::vector<unsigned int> bondLabels1, bondLabels2;
  get_bond_labels(*m1, *m2, opts, bondLabels1, bondLabels2);
  auto tier2_sim =
      tier2Sim(*m1, *m2, degSeqs1, degSeqs2, bondLabels1, bondLabels2);
  REQUIRE_THAT(tier2_sim, Catch::Matchers::WithinAbs(0.5632, 0.0001));

  // Note that this differs from the paper, because RDKit thinks both
  // rings in scopoletin are aromatic.
  opts.similarityThreshold = 0.5;
  auto res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 8);
  // symmetry means there is more than one solution.  Re-factoring may
  // produce different, equivalent, results.
  std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
  exp_bond_matches.push_back(
      {{0, 10}, {1, 11}, {2, 12}, {3, 14}, {4, 5}, {10, 1}, {12, 9}, {13, 6}});
  exp_bond_matches.push_back(
      {{0, 7}, {1, 6}, {2, 5}, {3, 14}, {4, 12}, {10, 1}, {12, 9}, {13, 11}});
  exp_bond_matches.push_back(
      {{0, 7}, {1, 9}, {2, 11}, {3, 12}, {4, 14}, {6, 1}, {12, 6}, {13, 5}});
  exp_bond_matches.push_back(
      {{0, 10}, {1, 11}, {2, 12}, {3, 14}, {4, 5}, {6, 1}, {12, 9}, {13, 6}});
  check_expected_bonds(res.front(), exp_bond_matches);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.3691, 0.0001));

  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("Methadone vs mepiridine test", "[basics]") {
  auto m1 = "c1ccccc1C(C(=O)CC)(c1ccccc1)CC(C)N(C)C"_smiles;
  REQUIRE(m1);
  auto m2 = "c1ccccc1C1(CCN(C)CC1)C(=O)OCC"_smiles;
  REQUIRE(m2);

  // It's a general requirement that the first mol is smaller than the second.
  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  auto tier1_sim = tier_1_sim(*m1, *m2, degSeqs1, degSeqs2);
  REQUIRE_THAT(tier1_sim, Catch::Matchers::WithinAbs(0.7044, 0.0001));

  std::vector<unsigned int> bondLabels1, bondLabels2;
  RascalOptions opts;
  get_bond_labels(*m1, *m2, opts, bondLabels1, bondLabels2);
  auto tier2_sim =
      tier2Sim(*m1, *m2, degSeqs1, degSeqs2, bondLabels1, bondLabels2);
  REQUIRE_THAT(tier2_sim, Catch::Matchers::WithinAbs(0.6262, 0.0001));

  opts.similarityThreshold = 0.6;
  auto res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 16);
  std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
  // doing this with a full initialiser confuses the auto-formatter in a
  // very unsavoury way.
  exp_bond_matches.push_back({{0, 3},
                              {1, 2},
                              {2, 1},
                              {3, 0},
                              {4, 17},
                              {5, 5},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {10, 6},
                              {16, 18},
                              {17, 11},
                              {19, 10},
                              {20, 8},
                              {21, 9},
                              {22, 4}});
  exp_bond_matches.push_back({{5, 18},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {10, 5},
                              {11, 4},
                              {12, 3},
                              {13, 2},
                              {14, 1},
                              {15, 0},
                              {16, 6},
                              {17, 7},
                              {19, 8},
                              {20, 9},
                              {21, 10},
                              {23, 17}});
  exp_bond_matches.push_back({{5, 18},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {10, 5},
                              {11, 17},
                              {12, 0},
                              {13, 1},
                              {14, 2},
                              {15, 3},
                              {16, 6},
                              {17, 7},
                              {19, 8},
                              {20, 10},
                              {21, 9},
                              {23, 4}});
  exp_bond_matches.push_back({{5, 6},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {10, 5},
                              {11, 17},
                              {12, 0},
                              {13, 1},
                              {14, 2},
                              {15, 3},
                              {16, 18},
                              {17, 11},
                              {19, 10},
                              {20, 8},
                              {21, 9},
                              {23, 4}});
  exp_bond_matches.push_back({{0, 3},
                              {1, 2},
                              {2, 1},
                              {3, 0},
                              {4, 17},
                              {5, 5},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {10, 18},
                              {16, 6},
                              {17, 7},
                              {19, 8},
                              {20, 10},
                              {21, 9},
                              {22, 4}});
  REQUIRE(res.front().smarts() == "c1ccccc1C(-C=O)(-[#6])-CCN(-C)-C.CC");

  check_expected_bonds(res.front(), exp_bond_matches);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.6262, 0.0001));

  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("testosterone vs estradiol", "[basics]") {
  auto m1 = "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"_smiles;
  REQUIRE(m1);
  auto m2 = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"_smiles;
  REQUIRE(m2);

  RascalOptions opts;
  opts.similarityThreshold = 0.6;
  auto res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 16);
  std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
  exp_bond_matches.push_back({{0, 0},
                              {1, 1},
                              {2, 2},
                              {3, 3},
                              {4, 4},
                              {5, 5},
                              {6, 6},
                              {7, 7},
                              {8, 8},
                              {9, 9},
                              {10, 10},
                              {11, 11},
                              {12, 12},
                              {20, 19},
                              {21, 20},
                              {22, 21}});
  std::sort(exp_bond_matches.begin(), exp_bond_matches.end());
  bool matched = false;
  for (const auto &ebm : exp_bond_matches) {
    if (ebm == res.front().bondMatches()) {
      matched = true;
      break;
    }
  }
  REQUIRE(matched);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.4966, 0.0001));
  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("Symmetrical esters test", "[basics]") {
  auto m1 = "c1c(OC)c(OC)c(OC)cc1C(=O)OCCCOC(=O)c1cc(OC)c(OC)c(OC)c1"_smiles;
  REQUIRE(m1);
  auto m2 = "c1c(OC)c(OC)c(OC)cc1C(=O)OCCOC(=O)c1cc(OC)c(OC)c(OC)c1"_smiles;
  REQUIRE(m2);

  // It's a general requirement that the first mol is smaller than the second.
  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  auto tier1_sim = tier_1_sim(*m1, *m2, degSeqs1, degSeqs2);
  REQUIRE_THAT(tier1_sim, Catch::Matchers::WithinAbs(0.9701, 0.0001));

  std::vector<unsigned int> bondLabels1, bondLabels2;
  RascalOptions opts;
  get_bond_labels(*m1, *m2, opts, bondLabels1, bondLabels2);
  auto tier2_sim =
      tier2Sim(*m1, *m2, degSeqs1, degSeqs2, bondLabels1, bondLabels2);
  REQUIRE_THAT(tier2_sim, Catch::Matchers::WithinAbs(0.9701, 0.0001));

  opts.similarityThreshold = 0.7;
  auto res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 32);
  std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
  exp_bond_matches.push_back(
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},   {6, 6},
       {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11}, {12, 12}, {13, 13},
       {14, 14}, {15, 15}, {18, 17}, {19, 18}, {20, 19}, {21, 32}, {22, 30},
       {23, 28}, {24, 29}, {25, 27}, {26, 25}, {27, 26}, {28, 24}, {29, 22},
       {30, 23}, {31, 21}, {32, 31}, {33, 20}});
  exp_bond_matches.push_back(
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},   {6, 6},
       {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11}, {12, 12}, {13, 13},
       {14, 14}, {17, 16}, {18, 17}, {19, 18}, {20, 19}, {21, 20}, {22, 21},
       {23, 22}, {24, 23}, {25, 24}, {26, 25}, {27, 26}, {28, 27}, {29, 28},
       {30, 29}, {31, 30}, {32, 31}, {33, 32}});
  exp_bond_matches.push_back(
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},   {6, 6},
       {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11}, {12, 12}, {13, 13},
       {14, 14}, {15, 15}, {18, 17}, {19, 18}, {20, 19}, {21, 20}, {22, 21},
       {23, 22}, {24, 23}, {25, 24}, {26, 25}, {27, 26}, {28, 27}, {29, 28},
       {30, 29}, {31, 30}, {32, 31}, {33, 32}});
  exp_bond_matches.push_back(
      {{0, 30},  {1, 28},  {2, 29},  {3, 27},  {4, 25},  {5, 26},  {6, 24},
       {7, 22},  {8, 23},  {9, 21},  {10, 20}, {11, 19}, {12, 18}, {13, 17},
       {16, 15}, {17, 14}, {18, 13}, {19, 12}, {20, 11}, {21, 10}, {22, 9},
       {23, 7},  {24, 8},  {25, 6},  {26, 4},  {27, 5},  {28, 3},  {29, 1},
       {30, 2},  {31, 0},  {32, 32}, {33, 31}});
  exp_bond_matches.push_back(
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},   {6, 6},
       {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11}, {12, 12}, {13, 13},
       {14, 14}, {17, 16}, {18, 17}, {19, 18}, {20, 19}, {21, 32}, {22, 30},
       {23, 28}, {24, 29}, {25, 27}, {26, 25}, {27, 26}, {28, 24}, {29, 22},
       {30, 23}, {31, 21}, {32, 31}, {33, 20}});
  exp_bond_matches.push_back(
      {{0, 9},   {1, 7},   {2, 8},   {3, 6},   {4, 4},   {5, 5},   {6, 3},
       {7, 1},   {8, 2},   {9, 0},   {10, 31}, {11, 11}, {12, 12}, {13, 13},
       {14, 14}, {17, 16}, {18, 17}, {19, 18}, {20, 19}, {21, 32}, {22, 30},
       {23, 28}, {24, 29}, {25, 27}, {26, 25}, {27, 26}, {28, 24}, {29, 22},
       {30, 23}, {31, 21}, {32, 10}, {33, 20}});
  exp_bond_matches.push_back(
      {{0, 30},  {1, 28},  {2, 29},  {3, 27},  {4, 25},  {5, 26},  {6, 24},
       {7, 22},  {8, 23},  {9, 21},  {10, 20}, {11, 19}, {12, 18}, {13, 17},
       {14, 16}, {17, 14}, {18, 13}, {19, 12}, {20, 11}, {21, 10}, {22, 9},
       {23, 7},  {24, 8},  {25, 6},  {26, 4},  {27, 5},  {28, 3},  {29, 1},
       {30, 2},  {31, 0},  {32, 32}, {33, 31}});
  exp_bond_matches.push_back(
      {{0, 9},   {1, 7},   {2, 8},   {3, 6},   {4, 4},   {5, 5},   {6, 3},
       {7, 1},   {8, 2},   {9, 0},   {10, 31}, {11, 11}, {12, 12}, {13, 13},
       {16, 15}, {17, 16}, {18, 17}, {19, 18}, {20, 19}, {21, 32}, {22, 30},
       {23, 28}, {24, 29}, {25, 27}, {26, 25}, {27, 26}, {28, 24}, {29, 22},
       {30, 23}, {31, 21}, {32, 10}, {33, 20}});
  REQUIRE(
      res.front().smarts() ==
      "c1c(-OC):c(-OC):c(-OC):cc1C(=O)-O.CCOC(=O)-c1cc(-OC):c(-OC):c(-OC):c1");
  check_expected_bonds(res.front(), exp_bond_matches);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.9405, 0.0001));
  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("dyphylline similarities") {
  auto dyphylline = "OCC(O)CN1C=NC2=C1C(=O)N(C)C(=O)N(C)2"_smiles;
  // The paper has the similarity scores as 0.78, 0.57, 0.51, 0.63 and 0.35
  // respectively. This implementation makes Viagra 0.11 because of the default
  // completeAromaticRings. Without that, it comes to 0.26.  The 0.63 for
  // enprofylline could be a typo, as this gets 0.73.  The 0.51 for captogon is
  // a mystery - with 15 atoms and 16 bonds in the MCES, 0.48 is correct.  To
  // get 0.51 there would need to be 16 atoms in the MCES which doesn't look
  // right.
  std::vector<std::tuple<std::shared_ptr<RDKit::RWMol>, std::string, double>> mols{
      {"CN1C=NC2=C1C(=O)N(C)C(=O)N(C)2"_smiles, "caffeine", 0.78},
      {"C12C(=O)NC(=O)NC(NC(=O)N2)=1"_smiles, "uric acid", 0.57},
      {"c1ccccc1CC(C)N(C)CCN1C=NC2=C1C(=O)N(C)C(=O)N(C)2"_smiles, "captagon",
       0.48},
      {"CCCN1C(=O)NC(=O)C2N=CNC1=2"_smiles, "enprofylline", 0.73},
      {"CN1CCN(CC1)S(=O)(=O)c1ccc(OCC)c(c1)C1=NC(=O)C2N(C)N=C(CCC)C(N1)=2"_smiles,
       "viagra", 0.11},
      {"OCCN(C)C1=NC2N(C)C(=O)N(C)C(=O)C(N1C)=2"_smiles, "cafaminol", 0.80}};
  RascalOptions opts;
  opts.similarityThreshold = 0.3;
  for (size_t i = 0; i < 6; ++i) {
    auto m = mols[i];
    auto res = rascalMces(*dyphylline, *std::get<0>(m), opts);
    REQUIRE(res.size() == 1);
    REQUIRE_THAT(res.front().similarity(),
                 Catch::Matchers::WithinAbs(std::get<2>(m), 0.01));
    check_smarts_ok(*dyphylline, *std::get<0>(m), res.front());
  }
}

TEST_CASE("compare chirality") {
  auto m1 = "OC(C)Cl"_smiles;
  REQUIRE(m1);
  auto m2 = "O[C@H](C)Cl"_smiles;
  REQUIRE(m2);
  auto m3 = "O[C@@H](C)Cl"_smiles;
  REQUIRE(m3);
  RascalOptions opts;
  {
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().smarts() == "OC(-C)-Cl");
    check_smarts_ok(*m1, *m2, res.front());
  }
  opts.exactChirality = true;
  opts.similarityThreshold = 0.3;
  {
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.empty());
  }
  {
    auto res = rascalMces(*m2, *m3, opts);
    REQUIRE(res.empty());
  }
  {
    auto res = rascalMces(*m2, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().smarts() == "O-[#6@&A](-C)-Cl");
    check_smarts_ok(*m2, *m2, res.front());
  }
  {
    auto res = rascalMces(*m3, *m3, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().smarts() == "O-[#6@@&A](-C)-Cl");
    check_smarts_ok(*m3, *m3, res.front());
  }
}

TEST_CASE("delta-y exchange", "[basics]") {
  RascalOptions opts;
  opts.similarityThreshold = 0.1;
  {
    // A delta-y exchange causes a false match, which in this
    // case will have the cyclopropyl ring matching part of the
    // t-butyl, a 3-bond match.
    auto m1 = "CC1CC1"_smiles;
    REQUIRE(m1);
    auto m2 = "CC(C)(C)C"_smiles;
    REQUIRE(m2);

    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == 3);
    std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
    exp_bond_matches.push_back({{0, 1}, {1, 3}, {3, 0}});
    exp_bond_matches.push_back({{0, 2}, {1, 3}, {3, 1}});
    exp_bond_matches.push_back({{0, 1}, {1, 2}, {3, 3}});
    check_expected_bonds(res.front(), exp_bond_matches);
    REQUIRE_THAT(res.front().similarity(),
                 Catch::Matchers::WithinAbs(0.6806, 0.0001));
    check_smarts_ok(*m1, *m2, res.front());
  }
  {
    // A delta-y exchange causes a false match, which in this
    // case will have the cyclopropyl ring matching part of the
    // t-butyl, a 3-bond match.
    auto m1 = "C1CCCCC12CC2"_smiles;
    REQUIRE(m1);
    auto m2 = "C1CCCCC1(C)(C)"_smiles;
    REQUIRE(m2);

    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == 8);
    std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
    exp_bond_matches.push_back(
        {{0, 3}, {1, 2}, {2, 1}, {3, 0}, {4, 7}, {5, 6}, {7, 4}, {8, 5}});
    exp_bond_matches.push_back(
        {{0, 3}, {1, 2}, {2, 1}, {3, 0}, {4, 7}, {5, 5}, {7, 4}, {8, 6}});
    check_expected_bonds(res.front(), exp_bond_matches);
    REQUIRE_THAT(res.front().similarity(),
                 Catch::Matchers::WithinAbs(0.9412, 0.0001));
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("bad aromatics 1") {
  std::vector<std::tuple<std::string, std::string, unsigned int, unsigned int>>
      tests = {
          {"c1ccccc1C(=O)c1ccncc1", "c1ccccc1C(=O)c1ccccc1", 9, 13},
          {"c1ccccc1C(=O)c1cc[nH]c1", "c1ccccc1C(=O)c1cnccc1", 9, 13},
          {"c1ccccc1C(=O)c1ccncc1", "c1ccccc1C(=O)c1cnccc1", 14, 14},
          {"c1ccccc1C(=O)c1ccc2cnccc2c1", "c1ccccc1C(=O)c1cnccc1", 14, 14}};

  RascalOptions opts;
  opts.similarityThreshold = 0.7;
  for (auto &test : tests) {
    std::unique_ptr<RDKit::RWMol> m1(RDKit::SmilesToMol(std::get<0>(test)));
    std::unique_ptr<RDKit::RWMol> m2(RDKit::SmilesToMol(std::get<1>(test)));
    opts.completeAromaticRings = true;
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == std::get<2>(test));
    check_smarts_ok(*m1, *m2, res.front());

    opts.completeAromaticRings = false;
    res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == std::get<3>(test));
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("timeout") {
  auto m1 =
      "O[C@@H]1CC[C@H](C[C@H]1OC)C[C@@H](C)[C@@H]4CC(=O)[C@H](C)/C=C(\\C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\\C=C\\C=C\\C=C(/C)[C@@H](OC)C[C@@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N3CCCC[C@H]3C(=O)O4"_smiles;
  REQUIRE(m1);
  auto m2 =
      "CCC1C=C(CC(CC(C2C(CC(C(O2)(C(=O)C(=O)N3CCCCC3C(=O)OC(C(C(CC1=O)O)C)C(=CC4CCC(C(C4)OC)O)C)O)C)OC)OC)C)C"_smiles;
  REQUIRE(m2);

  {
    RascalOptions opts;
    opts.timeout = 10;
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() >= 39);
    check_smarts_ok(*m1, *m2, res.front());
  }
  {
    RascalOptions opts;
    opts.timeout = 70;
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() >= 44);
    check_smarts_ok(*m1, *m2, res.front());
  }
  {
    RascalOptions opts;
    opts.timeout = 120;
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() >= 44);
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("single fragment") {
  RascalOptions opts;
  std::vector<std::tuple<std::string, std::string, unsigned int, unsigned int>>
      tests = {
          {"C1CC1CCC1NC1", "C1CC1CCCCC1NC1", 8, 6},
          {"c1cnccc1CCc1ncccc1", "c1cnccc1CCCCCCc1ncccc1", 14, 9},
          {"c1ccccc1c1cccc(c1)CCc1ccccc1", "c1ccccc1c1cccc(c1)CCCCCc1ccccc1",
           21, 16},
          {"Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008",
           "COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972",
           27, 21}};
  opts.similarityThreshold = 0.7;
  for (auto &test : tests) {
    opts.ringMatchesRingOnly = true;
    opts.singleLargestFrag = false;
    std::unique_ptr<RDKit::RWMol> m1(RDKit::SmilesToMol(std::get<0>(test)));
    std::unique_ptr<RDKit::RWMol> m2(RDKit::SmilesToMol(std::get<1>(test)));
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.front().numFrags() == 2);
    REQUIRE(res.front().bondMatches().size() == std::get<2>(test));
    check_smarts_ok(*m1, *m2, res.front());
    opts.singleLargestFrag = true;
    res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.front().numFrags() == 1);
    REQUIRE(res.front().bondMatches().size() == std::get<3>(test));
    REQUIRE(res.front().largestFragSize() == res.front().atomMatches().size());
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("minimum fragment sizes") {
  RascalOptions opts;

  std::vector<std::tuple<std::string, std::string, unsigned int, unsigned int>>
      tests = {
          {"Oc1cccc2C(=O)C=CC(=O)c12", "O1C(=O)C=Cc2cc(OC)c(O)cc12", 8, 7},
          {"c1ccccc1c1cccc(c1)CCc1ccccc1", "c1ccccc1c1cccc(c1)CCCCCc1ccccc1",
           21, 21},
      };
  opts.similarityThreshold = 0.5;
  for (auto &test : tests) {
    std::unique_ptr<RDKit::RWMol> m1(RDKit::SmilesToMol(std::get<0>(test)));
    std::unique_ptr<RDKit::RWMol> m2(RDKit::SmilesToMol(std::get<1>(test)));
    opts.minFragSize = 1;
    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == std::get<2>(test));
    check_smarts_ok(*m1, *m2, res.front());

    opts.minFragSize = 3;
    res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().bondMatches().size() == std::get<3>(test));
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("maximum fragment separation") {
  RascalOptions opts;

  std::vector<std::tuple<std::string, std::string, unsigned int, unsigned int>>
      tests = {
          {"c1ccccc1CCc1ccccc1", "c1ccccc1CCCCCCCCc1ccccc1", 14, 9},
          {"c1ccccc1c1cccc(c1)CCc1ccccc1", "c1ccccc1c1cccc(c1)CCCCCCCc1ccccc1",
           21, 16},
      };
  opts.similarityThreshold = 0.5;
  for (auto &test : tests) {
    std::unique_ptr<RDKit::RWMol> m1(RDKit::SmilesToMol(std::get<0>(test)));
    std::unique_ptr<RDKit::RWMol> m2(RDKit::SmilesToMol(std::get<1>(test)));
    {
      opts.maxFragSeparation = -1;
      auto res = rascalMces(*m1, *m2, opts);
      REQUIRE(res.size() == 1);
      REQUIRE(res.front().bondMatches().size() == std::get<2>(test));
      check_smarts_ok(*m1, *m2, res.front());
    }
    {
      opts.maxFragSeparation = 3;
      opts.timeout = 300;
      auto res = rascalMces(*m1, *m2, opts);
      REQUIRE(res.size() == 1);
      REQUIRE(res.front().bondMatches().size() == std::get<3>(test));
      check_smarts_ok(*m1, *m2, res.front());
    }
  }
}

TEST_CASE("ring matches ring", "[basics]") {
  // This is methadone vs mepiridine.  The default result in this case is
  // very unsavoury due to 2 single bond matches.
  auto m1 = "c1ccccc1C(C(=O)CC)(c1ccccc1)CC(C)N(C)C"_smiles;
  REQUIRE(m1);
  auto m2 = "c1ccccc1C1(CCN(C)CC1)C(=O)OCC"_smiles;
  REQUIRE(m2);

  RascalOptions opts;
  opts.similarityThreshold = 0.6;
  opts.ringMatchesRingOnly = true;
  auto res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 11);
  std::vector<std::vector<std::pair<int, int>>> exp_bond_matches;
  // doing this with a full initialiser confuses the auto-formatter in a
  // very unsightly way.
  exp_bond_matches.push_back({{0, 0},
                              {1, 1},
                              {2, 2},
                              {3, 3},
                              {4, 4},
                              {5, 5},
                              {6, 12},
                              {7, 13},
                              {9, 16},
                              {19, 9},
                              {22, 17}});
  exp_bond_matches.push_back({{6, 12},
                              {7, 13},
                              {10, 5},
                              {11, 4},
                              {12, 3},
                              {13, 2},
                              {14, 1},
                              {15, 0},
                              {17, 16},
                              {21, 9},
                              {23, 17}});
  exp_bond_matches.push_back({{6, 12},
                              {7, 13},
                              {10, 5},
                              {11, 17},
                              {12, 0},
                              {13, 1},
                              {14, 2},
                              {15, 3},
                              {17, 16},
                              {21, 9},
                              {23, 4}});
  exp_bond_matches.push_back({{6, 12},
                              {7, 13},
                              {10, 5},
                              {11, 17},
                              {12, 0},
                              {13, 1},
                              {14, 2},
                              {15, 3},
                              {17, 16},
                              {20, 9},
                              {23, 4}});
  RascalMCES::printBondMatches(res.front(), std::cout);

  check_expected_bonds(res.front(), exp_bond_matches);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.3312, 0.0001));
  check_smarts_ok(*m1, *m2, res.front());

  // This reduces it to the single largest fragment
  opts.singleLargestFrag = true;
  res = rascalMces(*m1, *m2, opts);
  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 9);
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.1863, 0.0001));
  REQUIRE(res.front().smarts() == "C(-C=O)-c1ccccc1");
  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("multiple cliques returned") {
  std::vector<std::tuple<std::string, std::string, unsigned int, unsigned int,
                         std::vector<std::pair<int, int>>>>
      tests = {
          {"Oc1cccc2C(=O)C=CC(=O)c12",
           "O1C(=O)C=Cc2cc(OC)c(O)cc12",
           8,
           8,
           {{0, 10},
            {1, 9},
            {2, 6},
            {3, 5},
            {4, 14},
            {6, 1},
            {12, 11},
            {13, 12}}},
          {"Oc1cccc2C(=O)C(C)=CC(=O)c12",
           "O1C(=O)C(C)=Cc2cc(OC)c(O)cc12",
           16,
           9,
           {{0, 11},
            {1, 10},
            {2, 7},
            {3, 6},
            {4, 15},
            {8, 3},
            {11, 1},
            {13, 12},
            {14, 13}}},
      };

  RascalOptions opts;
  opts.similarityThreshold = 0.5;
  opts.allBestMCESs = true;

  for (auto &test : tests) {
    opts.allBestMCESs = true;
    std::unique_ptr<RDKit::RWMol> m1(RDKit::SmilesToMol(std::get<0>(test)));
    std::unique_ptr<RDKit::RWMol> m2(RDKit::SmilesToMol(std::get<1>(test)));

    auto res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == std::get<2>(test));
    REQUIRE(res.front().bondMatches().size() == std::get<3>(test));
    REQUIRE(res.front().bondMatches() == std::get<4>(test));
    check_smarts_ok(*m1, *m2, res.front());

    opts.allBestMCESs = false;
    res = rascalMces(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    // The clique won't necessarily be the best-scoring one, but it should be
    // the same size.
    REQUIRE(res.front().bondMatches().size() == std::get<4>(test).size());
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("benchmarks") {
  // As well as timings, this also tests that the same result is produced from
  // random ordering of the input molecules i.e. the results aren't input order
  // dependent.
  std::vector<std::tuple<std::string, std::string, std::string, double, int,
                         unsigned int, int>>
      tests = {{"juglone", "Oc1cccc2C(=O)C=CC(=O)c12",
                "O1C(=O)C=Cc2cc(OC)c(O)cc12", 0.5, 100, 8, 10},
               {"methadone", "c1ccccc1C(C(=O)CC)(c1ccccc1)CC(C)N(C)C",
                "c1ccccc1C1(CCN(C)CC1)C(=O)OCC", 0.6, 100, 16, 20},
               {"symmetrical",
                "c1c(OC)c(OC)c(OC)cc1C(=O)OCCCOC(=O)c1cc(OC)c(OC)c(OC)c1",
                "c1c(OC)c(OC)c(OC)cc1C(=O)OCCOC(=O)c1cc(OC)c(OC)c(OC)c1", 0.7,
                100, 32, 50},
               {"testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
                "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O", 0.6, 100, 16, 40}};
  std::vector<double> timings;
  std::random_device rd;
  std::mt19937 g(rd());
  // It's convenient if the same test is run each time.
  g.seed(1);
  RascalOptions opts;
  for (const auto &test : tests) {
    auto m1 = std::unique_ptr<ROMol>(RDKit::SmilesToMol(std::get<1>(test)));
    auto m2 = std::unique_ptr<ROMol>(RDKit::SmilesToMol(std::get<2>(test)));
    opts.similarityThreshold = std::get<3>(test);
    std::vector<unsigned int> atom1Inds(m1->getNumAtoms(), 0);
    std::iota(atom1Inds.begin(), atom1Inds.end(), 0);
    std::vector<unsigned int> atom2Inds(m2->getNumAtoms(), 0);
    std::iota(atom2Inds.begin(), atom2Inds.end(), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < std::get<4>(test); ++i) {
      std::shuffle(atom1Inds.begin(), atom1Inds.end(), g);
      std::unique_ptr<ROMol> randmol1(MolOps::renumberAtoms(*m1, atom1Inds));
      std::shuffle(atom2Inds.begin(), atom2Inds.end(), g);
      std::unique_ptr<ROMol> randmol2(MolOps::renumberAtoms(*m2, atom2Inds));

      auto res = rascalMces(*randmol1, *randmol2, opts);
      REQUIRE(res.size() == 1);
      REQUIRE(std::get<5>(test) == res.front().bondMatches().size());
      check_smarts_ok(*randmol1, *randmol2, res.front());
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    timings.push_back(
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() /
        std::get<4>(test));
  }

  // Make sure we haven't slowed it down a lot.
  //  for (size_t i = 0; i < tests.size(); ++i) {
  //    std::cout << "TIMING : " << std::get<0>(tests[i]) << " took " <<
  //    timings[i]
  //              << " milliseconds on average of " << std::get<4>(tests[i])
  //              << " runs\n";
  //  }
  for (size_t i = 0; i < tests.size(); ++i) {
    REQUIRE(timings[i] < std::get<6>(tests[i]));
  }
}

// Anything starting FMCS is taken from the test data for RDKit's FMCS tests.
TEST_CASE("FMCS test1Basics") {
  auto m1 = "CC1CCC(N)CC1"_smiles;
  REQUIRE(m1);
  auto m2 = "CC1CC(C)CC(C)C1"_smiles;
  REQUIRE(m2);

  RascalOptions opts;
  opts.similarityThreshold = 0.6;
  auto res = rascalMces(*m1, *m2, opts);

  REQUIRE(res.size() == 1);
  REQUIRE(res.front().bondMatches().size() == 7);
  REQUIRE(res.front().smarts() == "CC1CCCCC1");
  REQUIRE_THAT(res.front().similarity(),
               Catch::Matchers::WithinAbs(0.6806, 0.0001));

  check_smarts_ok(*m1, *m2, res.front());
}

TEST_CASE("FMCS test32") {
  std::vector<std::shared_ptr<RDKit::ROMol>> mols{
      {"O=C(Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCOCC1)C=Cc1ccc(Cl)cc1 CHEMBL1515359"_smiles},
      {"c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)=O)cc1 CHEMBL1590658"_smiles},
      {"Cc1ccc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)cc1 CHEMBL1447567"_smiles},
      {"c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCCC2)=O)cc1 CHEMBL1384017"_smiles},
      {"O=C(C=Cc1ccc(F)cc1)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1 CHEMBL1456416"_smiles},
      {"c1cc(F)cc(C=CC(=O)Nc2c(N3CCCC3)ccc(S(N3CCOCC3)(=O)=O)c2)c1 CHEMBL1308819"_smiles},
      {"CCN1CCN(c2ccc(S(N3CCOCC3)(=O)=O)cc2NC(=O)C=Cc2ccc(C)cc2)CC1 CHEMBL1703007"_smiles},
      {"c1cc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)c([N+]([O-])=O)cc1 CHEMBL1707819"_smiles},
      {"N#CC(=Cc1ccccc1)C(=O)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1 CHEMBL1500793"_smiles},
      {"C(=Cc1ccc2c(c1)OCO2)C(Nc1cc(S(=O)(=O)N2CCOCC2)ccc1N1CCOCC1)=O CHEMBL1334715"_smiles}};
  RascalOptions opts;
  std::vector<std::tuple<unsigned int, unsigned int, std::string>> exp_res{
      {32, 35, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCOCC1)-C=Cc1ccccc1"},
      {32, 35, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCOCC1)-C=Cc1ccccc1"},
      {31, 33, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC)-C=Cc1ccccc1"},
      {31, 33, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC)-C=Cc1ccccc1"},
      {31, 33, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC)-C=Cc1ccccc1"},
      {31, 33, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC)-C=Cc1ccccc1"},
      {32, 35, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCOCC1)-C=Cc1ccccc1"},
      {31, 33, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC)-C=Cc1ccccc1"},
      {32, 35, "O=C(-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCOCC1)-C=Cc1ccccc1"},
      {32, 35, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {32, 35, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {32, 35, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):cc1"},
      {31, 33, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):cc1"},
      {31, 33, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):cc1"},
      {32, 34,
       "Cc1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):cc1"},
      {32, 35, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2):cc1"},
      {31, 33, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):cc1"},
      {32, 35, "c1ccc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2):cc1"},
      {31, 34, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCCC2)=O):cc1"},
      {31, 34, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCCC2)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {31, 34, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCCC2)=O):cc1"},
      {31, 33, "c1ccc(-C=CC(-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC)=O):cc1"},
      {32, 34,
       "O=C(-C=C)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCCC1.c1ccc(-F):cc1"},
      {31, 33, "O=C(-C=Cc1ccccc1)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC"},
      {31, 33, "O=C(-C=Cc1ccccc1)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC"},
      {31, 34, "O=C(-C=Cc1ccccc1)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N1CCCC1"},
      {31, 33, "O=C(-C=Cc1ccccc1)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC"},
      {31, 33,
       "c1cccc(-C=CC(=O)-Nc2c(-N(-CC)-CC):ccc(-S(-N3CCOCC3)(=O)=O):c2):c1"},
      {31, 33,
       "c1cccc(-C=CC(=O)-Nc2c(-N(-CC)-CC):ccc(-S(-N3CCOCC3)(=O)=O):c2):c1"},
      {31, 34,
       "c1cccc(-C=CC(=O)-Nc2c(-N3CCCC3):ccc(-S(-N3CCOCC3)(=O)=O):c2):c1"},
      {31, 33,
       "c1cccc(-C=CC(=O)-Nc2c(-N(-CC)-CC):ccc(-S(-N3CCOCC3)(=O)=O):c2):c1"},
      {31, 33, "CCN(-c1ccc(-S(-N2CCOCC2)(=O)=O):cc1NC(=O)-C=Cc1ccccc1)-CC"},
      {31, 33, "CCN(-c1ccc(-S(-N2CCOCC2)(=O)=O):cc1NC(=O)-C=Cc1ccccc1)-CC"},
      {31, 33, "CCN(-c1ccc(-S(-N2CCOCC2)(=O)=O):cc1NC(=O)-C=Cc1ccccc1)-CC"},
      {31, 33, "c1cc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N(-CC)-CC):ccc1"},
      {32, 35, "c1cc(-C=CC(=O)-Nc2cc(-S(-N3CCOCC3)(=O)=O):ccc2N2CCOCC2):ccc1"},
      {31, 33, "C(=Cc1ccccc1)-C(=O)-Nc1cc(-S(-N2CCOCC2)(=O)=O):ccc1N(-CC)-CC"}};
  size_t k = 0;
  for (size_t i = 0; i < mols.size() - 1; ++i) {
    for (size_t j = i + 1; j < mols.size(); ++j, ++k) {
      auto res = rascalMces(*mols[i], *mols[j], opts);
      check_smarts_ok(*mols[i], *mols[j], res.front());
      REQUIRE(res.front().atomMatches().size() == std::get<0>(exp_res[k]));
      REQUIRE(res.front().bondMatches().size() == std::get<1>(exp_res[k]));
      REQUIRE(res.front().smarts() == std::get<2>(exp_res[k]));
    }
  }
}

TEST_CASE("FMCS test190") {
  std::vector<std::shared_ptr<RDKit::ROMol>> mols{
      {"COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)oc2cc1  CHEMBL1479679"_smiles},
      {"COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1333382"_smiles},
      {"Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)nc2cc1  CHEMBL1437584"_smiles},
      {"COc1c(NC(=O)CSc2ccc(Cl)cc2)cc(-c2nc3ccccc3o2)cc1  CHEMBL1601350"_smiles},
      {"Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008"_smiles},
      {"Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)nc2cc1  CHEMBL1612903"_smiles},
      {"COc1cc2nc(-c3cc(NC(=O)Cc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1316483"_smiles},
      {"Cc1c(NC(=O)CSc2ccc(Cl)cc2)cccc1-c1nc2cc(Cl)ccc2o1  CHEMBL1568754"_smiles},
      {"COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972"_smiles},
      {"Cc1ccc(SCC(=O)Nc2cc(-c3nc4cc(C)ccc4o3)c(O)cc2)cc1  CHEMBL1611932"_smiles},
  };
  RascalOptions opts;
  opts.allBestMCESs = true;
  std::vector<std::tuple<unsigned int, unsigned int, std::string>> exp_res{
      {29, 32, "COc1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {29, 31, "CO.c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {28, 30, "COc1cc2nc(-c3cc(-NC(=O)-C):ccc3):oc2cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {27, 29, "COc1cc2nc(-c3cc(-NC(=O)-C):ccc3):oc2cc1.c1ccccc1"},
      {26, 29, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccccc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {29, 31, "CO.c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {28, 31, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):c(-C):cc3):oc2cc1"},
      {29, 31, "COc1cc2nc(-c3cc(-NC(=O)-C):c(-C):cc3):oc2cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):oc2cc1"},
      {28, 30, "COc1cc2nc(-c3cc(-NC(=O)-C):c(-C):cc3):oc2cc1.c1ccccc1"},
      {26, 29, "c1cc2nc(-c3cc(-NC(=O)-CSc4ccccc4):ccc3):oc2cc1"},
      {27, 30, "c1cc2oc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):nc2cc1"},
      {28, 31, "Cc1cc2oc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):nc2cc1"},
      {28, 31, "Cc1cc2oc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):nc2cc1"},
      {26, 28, "c1cc2oc(-c3cc(-NC(=O)-C):ccc3):nc2cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1cc2oc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):nc2cc1"},
      {25, 27, "c1cc2oc(-c3cc(-NC(=O)-C):ccc3):nc2cc1.c1ccccc1"},
      {26, 29, "c1cc2oc(-c3cc(-NC(=O)-CSc4ccccc4):ccc3):nc2cc1"},
      {27, 30, "c1c(-NC(=O)-CSc2ccc(-Cl):cc2):cc(-c2nc3ccccc3o2):cc1"},
      {27, 30, "c1c(-NC(=O)-CSc2ccc(-Cl):cc2):cc(-c2nc3ccccc3o2):cc1"},
      {28, 29, "CO.c1c(-NC(=O)-C):cc(-c2nc3ccccc3o2):cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1c(-NC(=O)-CSc2ccc(-Cl):cc2):cc(-c2nc3ccccc3o2):cc1"},
      {27, 28, "CO.c1c(-NC(=O)-C):cc(-c2nc3ccccc3o2):cc1.c1ccccc1"},
      {26, 29, "c1c(-NC(=O)-CSc2ccccc2):cc(-c2nc3ccccc3o2):cc1"},
      {28, 31, "c1cc2nc(-c3cccc(-NC(=O)-CSc4ccc(-Cl):cc4):c3):oc2cc1C"},
      {26, 28, "c1cc2nc(-c3cccc(-NC(=O)-C):c3):oc2cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1cc2nc(-c3cccc(-NC(=O)-CSc4ccc(-Cl):cc4):c3):oc2cc1"},
      {25, 27, "c1cc2nc(-c3cccc(-NC(=O)-C):c3):oc2cc1.c1ccccc1"},
      {27, 30, "Cc1cc2nc(-c3cccc(-NC(=O)-CSc4ccccc4):c3):oc2cc1"},
      {27, 29, "c1cc2oc(-c3cc(-NC(=O)-C):c(-C):cc3):nc2cc1.c1ccc(-Cl):cc1"},
      {27, 30, "c1cc2oc(-c3cc(-NC(=O)-CSc4ccc(-Cl):cc4):ccc3):nc2cc1"},
      {26, 28, "c1cc2oc(-c3cc(-NC(=O)-C):c(-C):cc3):nc2cc1.c1ccccc1"},
      {26, 29, "c1cc2oc(-c3cc(-NC(=O)-CSc4ccccc4):ccc3):nc2cc1"},
      {26, 28, "c1cc2nc(-c3cc(-NC(=O)-C):ccc3):oc2cc1.c1ccc(-Cl):cc1"},
      {28, 30, "COc1cc2nc(-c3cc(-NC(=O)-C):c(-C):cc3):oc2cc1.c1ccccc1"},
      {25, 27, "c1cc2nc(-c3cc(-NC(=O)-C):ccc3):oc2cc1.c1ccccc1"},
      {25, 27, "c1c(-NC(=O)-C):cccc1-c1nc2ccccc2o1.c1ccccc1"},
      {26, 29, "c1c(-NC(=O)-CSc2ccccc2):cccc1-c1nc2ccccc2o1"},
      {26, 28, "c1ccc2oc(-c3cccc(-NC(=O)-C):c3):nc2c1.c1cc(-C):ccc1"}};
  size_t k = 0;
  for (size_t i = 0; i < mols.size() - 1; ++i) {
    for (size_t j = i + 1; j < mols.size(); ++j, ++k) {
      auto res = rascalMces(*mols[i], *mols[j], opts);
      check_smarts_ok(*mols[i], *mols[j], res.front());
      REQUIRE(res.front().atomMatches().size() == std::get<0>(exp_res[k]));
      REQUIRE(res.front().bondMatches().size() == std::get<1>(exp_res[k]));
      REQUIRE(res.front().smarts() == std::get<2>(exp_res[k]));
    }
  }
}

TEST_CASE("FMCS test3") {
  std::vector<std::shared_ptr<RDKit::ROMol>> mols{
      {"CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934"_smiles},
      {"CN(C)c1ccc(CC(=O)NCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL152361"_smiles},
      {"CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336"_smiles},
      {"CN(C)c1ccc(CC(=O)NCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157429"_smiles},
      {"CN(C)c1ccc(CC(=O)NCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL357551"_smiles},
      {"CN(C)c1ccc(CC(=O)NCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL421974"_smiles},
      {"CN(C)c1ccc(CC(NCCCCCC(NO)=O)=O)cc1 CHEMBL484488"_smiles},
      {"CC(C)Cc1ccc(C(C)C(=O)NC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL564780"_smiles},
      {"c1cc([N+]([O-])=O)ccc1CC(=O)NC1CCCCCC1 CHEMBL1553142"_smiles},
      {"CC1(C)NC(C)(C)CC(NC(=O)Cc2ccccc2)C1 CHEMBL1703640"_smiles},
  };
  RascalOptions opts;
  opts.similarityThreshold = 0.1;
  // Because different compilers can give different equivalent results,
  // get all MCESs, which should then be sorted into a consistent order
  // so the first one should always be the same.  This results in extra
  // run time, but means the tests should work on all platforms.
  //  opts.allBestMCESs = true;
  std::vector<std::tuple<unsigned int, unsigned int, std::vector<std::string>>>
      exp_res{
          {31,
           33,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {34,
           36,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {33,
           35,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {32,
           34,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {34,
           36,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {31,
           33,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {31,
           33,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {31,
           33,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {31,
           33,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {33,
           35,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {32,
           34,
           {{"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}}},
          {35,
           37,
           {{"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {32,
           34,
           {{"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}}},
          {33,
           35,
           {"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {32,
           34,
           {{"CN(-C)-c1ccc(-CC(=O)-NCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3"}}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {19, 19, {"CN(-C)-c1ccc(-CC(=O)-NCCCCCC):cc1"}},
          {24, 25, {"c1ccc(-C):cc1.CCC.CCCNC12CC3CC(-C1)-CC(-C2)-C3"}},
          {18, 18, {"Nc1ccc(-CC(=O)-NCCCCCCC):cc1"}},
          {19, 18, {{"c1ccc(-CC(=O)-NCCCC):cc1.NC(-C)(-C)-C"}}},
          {16, 16, {"c1ccc(-CC(-NCCCCCC)=O):cc1"}},
          {18, 17, {{"c1ccc(-CC(-NCCCCCC)=O):cc1.NO"}}},
          {17, 16, {"c1ccc(-CC(-NCCCC)=O):cc1.CCN"}},
          {17, 17, {{"c1ccc(-CC(=O)-NC(-CCCCC)-C):cc1"}}},
          {18, 18, {"c1ccc(-CC(=O)-NC(-CC(-C)-C)-CCC):cc1"}},
          {17, 17, {"c1ccccc1CC(=O)-NC(-CCC)-CCC"}}};
  size_t k = 0;
  for (size_t i = 0; i < mols.size() - 1; ++i) {
    for (size_t j = i + 1; j < mols.size(); ++j, ++k) {
      auto res = rascalMces(*mols[i], *mols[j], opts);
      check_smarts_ok(*mols[i], *mols[j], res.front());
      check_smarts_found(res.front(), get<2>(exp_res[k]));
      REQUIRE(res.front().atomMatches().size() == std::get<0>(exp_res[k]));
      REQUIRE(res.front().bondMatches().size() == std::get<1>(exp_res[k]));
    }
  }
}
