//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <chrono>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>

#include <GraphMol/RascalMCES/PartitionSet.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

namespace RDKit {
namespace RascalMCES {

class TimedOutException : public std::exception {
 public:
  TimedOutException(long long int run_time,
                    std::vector<std::vector<unsigned int>> &best_cliques)
      : d_cliques(best_cliques) {
    d_message = "Timed out after " + std::to_string(run_time) + " second";
    if (run_time == 1) {
      d_message += ".";
    } else {
      d_message += "s.";
    }
  }

  const char *what() const throw() override { return d_message.c_str(); }

  std::vector<std::vector<unsigned int>> d_cliques;

 private:
  std::string d_message;
};

// This is in lap_a_la_scipy.cpp and solves the linear assignment problem.
int lap_maximize(const std::vector<std::vector<int>> &costsMat,
                 std::vector<size_t> &a, std::vector<size_t> &b);

// Contains the information used to start off a Rascal job.
struct RascalStartPoint {
  // The tier one and tier two similarities.  The initial PartitionSet
  // will only be constructed if they are both above the RascalOptions
  // threshold.
  double d_tier1Sim{-1.0};
  double d_tier2Sim{-1.0};
  // The lower bound on the maximum clique size.  Depends, amongst other things,
  // on opts.similarityThreshold.
  int d_lowerBound{0};
  // a Delta-Y exchange requires extra treatment.  They're rare, though.
  bool d_deltaYPoss{false};

  // Some parts require mol2 to be the larger molecule.  This records if they
  // have been swapped with respect to the input molecules.
  bool d_swapped{false};
  // Pointers the input molecules, swapped if necessary, so that d_mol1 is
  // always the smaller molecule.  These should never be deleted.
  const ROMol *d_mol1;
  const ROMol *d_mol2;

  std::vector<std::vector<int>> d_adjMatrix1, d_adjMatrix2;
  std::vector<std::pair<int, int>> d_vtxPairs;
  std::vector<std::vector<char>> d_modProd;

  // We might need to know which bonds are symmetrical equivalent.
  std::vector<int> d_equivBonds1, d_equivBonds2;

  // The initial partition set.  May be empty if the thresholds weren't met.
  std::shared_ptr<PartitionSet> d_partSet;
};

// Get the sorted degree sequences for the molecule, one sequence for each
// atomic number in the molecule.  Each element in the degree sequence is
// the degree of the atom and its index.
void sorted_degree_seqs(
    const ROMol &mol,
    std::map<int, std::vector<std::pair<int, int>>> &deg_seqs) {
  for (const auto &a : mol.atoms()) {
    auto apos = deg_seqs.find(a->getAtomicNum());
    if (apos == deg_seqs.end()) {
      auto res = deg_seqs.insert(std::make_pair(
          a->getAtomicNum(), std::vector<std::pair<int, int>>()));
      apos = res.first;
    }
    apos->second.push_back(std::make_pair(a->getDegree(), a->getIdx()));
  }
  for (auto &it : deg_seqs) {
    std::sort(it.second.begin(), it.second.end(),
              [](const std::pair<int, int> &p1, const std::pair<int, int> &p2)
                  -> bool { return p1.first > p2.first; });
  }
}

// calculate the tier 1 similarity between the 2 mols.
double tier_1_sim(const ROMol &mol1, const ROMol &mol2,
                  std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                  std::map<int, std::vector<std::pair<int, int>>> &degSeqs2) {
  sorted_degree_seqs(mol1, degSeqs1);
  sorted_degree_seqs(mol2, degSeqs2);
  int vg1g2 = 0;
  int eg1g2 = 0;
  for (const auto &it1 : degSeqs1) {
    const auto &seq2 = degSeqs2.find(it1.first);
    if (seq2 != degSeqs2.end()) {
      vg1g2 += std::min(it1.second.size(), seq2->second.size());
      auto numToDo = std::min(it1.second.size(), seq2->second.size());
      for (auto i = 0; i < numToDo; ++i) {
        eg1g2 += std::min(it1.second[i].first, seq2->second[i].first);
      }
    }
  }
  eg1g2 /= 2;
  double sim = double((vg1g2 + eg1g2) * (vg1g2 + eg1g2)) /
               double((mol1.getNumAtoms() + mol1.getNumBonds()) *
                      (mol2.getNumAtoms() + mol2.getNumBonds()));
  // std::cout << "vg1g2 : " << vg1g2 << " and eg1g2 : " << eg1g2 << " gives
  // tier 1 sim = " << sim << std::endl;
  return sim;
}

// Make labels for the atoms - by default the atomic symbol.
// If includeChirality and the atom is chiral, adds (R) or (S).
void get_atom_labels(const ROMol &mol, const RascalOptions &opts,
                     std::vector<std::string> &atomLabels) {
  atomLabels = std::vector<std::string>(mol.getNumAtoms(), "");
  std::string cip;
  for (const auto &a : mol.atoms()) {
    std::string label = a->getSymbol();
    if (opts.exactChirality) {
      cip = "";
      a->getPropIfPresent(common_properties::_CIPCode, cip);
      if (!cip.empty()) {
        label += "(" + cip + ")";
      }
    }
    atomLabels[a->getIdx()] = label;
  }
}

// For each bond in the molecule, encode it with its type and the labels of the
// two end atoms, returning the results as strings.
// Note that the molecule should not be in kekulized form.
void get_bond_labels(const ROMol &mol, const RascalOptions &opts,
                     std::vector<std::string> &bondLabels) {
  std::vector<std::string> atomLabels;
  get_atom_labels(mol, opts, atomLabels);
  bondLabels = std::vector<std::string>(mol.getNumBonds(), "");
  for (const auto &b : mol.bonds()) {
    if (b->getBeginAtom()->getAtomicNum() < b->getEndAtom()->getAtomicNum()) {
      bondLabels[b->getIdx()] = atomLabels[b->getBeginAtomIdx()] +
                                std::to_string(b->getBondType()) +
                                atomLabels[b->getEndAtomIdx()];
    } else {
      bondLabels[b->getIdx()] = atomLabels[b->getEndAtomIdx()] +
                                std::to_string(b->getBondType()) +
                                atomLabels[b->getBeginAtomIdx()];
    }
    //        std::cout << b->getIdx() << " : " << bondLabels[b->getIdx()] <<
    //        std::endl;
  }
}

// Fills bondLabels[12] with a small integer denoting the type of the bond and
// the types of the atoms at each end.  Both molecules need to be done at the
// same time so that the labels are consistent across both.
void get_bond_labels(const ROMol &mol1, const ROMol &mol2,
                     const RascalOptions &opts,
                     std::vector<unsigned int> &bondLabels1,
                     std::vector<unsigned int> &bondLabels2) {
  std::vector<std::string> tmpBondLabels1, tmpBondLabels2;

  get_bond_labels(mol1, opts, tmpBondLabels1);
  get_bond_labels(mol2, opts, tmpBondLabels2);

  // convert the bond labels, which are currently encoding the atoms and
  // bond type to a small set of sequential integers for ease of use later.
  // This results in loss of information, but that information is not currently
  // used anywhere.
  std::set<std::string> allLabels;
  for (auto bl : tmpBondLabels1) {
    allLabels.insert(bl);
  }
  for (auto bl : tmpBondLabels2) {
    allLabels.insert(bl);
  }
  auto recodeBondLabels = [&](std::vector<std::string> &strBondLabels,
                              std::vector<unsigned int> &bondLabels) {
    for (auto &bl : strBondLabels) {
      auto it = allLabels.find(bl);
      bondLabels.push_back(std::distance(allLabels.begin(), it));
    }
  };
  recodeBondLabels(tmpBondLabels1, bondLabels1);
  recodeBondLabels(tmpBondLabels2, bondLabels2);
}

int calc_cost(const std::vector<unsigned int> &atomiBLs,
              const std::vector<unsigned int> &atomjBLs,
              std::vector<unsigned int> &uniqAtomiBLs) {
  uniqAtomiBLs.clear();
  uniqAtomiBLs = atomiBLs;
  std::sort(uniqAtomiBLs.begin(), uniqAtomiBLs.end());
  uniqAtomiBLs.erase(std::unique(uniqAtomiBLs.begin(), uniqAtomiBLs.end()),
                     uniqAtomiBLs.end());
  int cost = 0;
  for (const auto &uai : uniqAtomiBLs) {
    int numAtomi = std::count(atomiBLs.begin(), atomiBLs.end(), uai);
    int numAtomj = std::count(atomjBLs.begin(), atomjBLs.end(), uai);
    cost += std::min(numAtomi, numAtomj);
  }
  return cost;
}

// assign the costs of matching each atom in atomDegrees1 to each atom in
// atomDegrees2.
void assignCosts(const std::vector<std::pair<int, int>> &atomDegrees1,
                 const std::vector<std::pair<int, int>> &atomDegrees2,
                 const std::vector<unsigned int> &bondLabels1,
                 const std::vector<unsigned int> &bondLabels2,
                 const ROMol &mol1, const ROMol &mol2,
                 std::vector<std::vector<int>> &costsMat) {
  //    std::cout << "assign costs" << std::endl;
  std::vector<unsigned int> atomiBLs, atomjBLs, uniqAtomiBLs;
  for (auto i = 0u; i < atomDegrees1.size(); ++i) {
    atomiBLs.clear();
    const auto atomi = mol1.getAtomWithIdx(atomDegrees1[i].second);
    for (const auto b : mol1.atomBonds(atomi)) {
      atomiBLs.push_back(bondLabels1[b->getIdx()]);
    }
    for (auto j = 0u; j < atomDegrees2.size(); ++j) {
      const auto atomj = mol2.getAtomWithIdx(atomDegrees2[j].second);
      atomjBLs.clear();
      for (const auto b : mol2.atomBonds(atomj)) {
        atomjBLs.push_back(bondLabels2[b->getIdx()]);
      }
      costsMat[i][j] = calc_cost(atomiBLs, atomjBLs, uniqAtomiBLs);
    }
  }
}

// Return the assignment score for the best match of the atoms and bonds in mol1
// to the atoms and bonds in mol2.
int getAssignmentScore(const std::vector<std::pair<int, int>> &atomDegrees1,
                       const std::vector<std::pair<int, int>> &atomDegrees2,
                       const std::vector<unsigned int> &bondLabels1,
                       const std::vector<unsigned int> &bondLabels2,
                       const ROMol &mol1, const ROMol &mol2) {
  std::vector<std::vector<int>> costsMat(
      atomDegrees1.size(), std::vector<int>(atomDegrees2.size(), 9999));
  assignCosts(atomDegrees1, atomDegrees2, bondLabels1, bondLabels2, mol1, mol2,
              costsMat);
  std::vector<size_t> a(std::min(atomDegrees1.size(), atomDegrees2.size()),
                        99999999);
  std::vector<size_t> b(std::min(atomDegrees1.size(), atomDegrees2.size()),
                        99999999);
  int retVal = lap_maximize(costsMat, a, b);
  if (retVal < 0) {
    // no solution for the LAP was possible.
    return 0;
  }
  int totalCost = 0;
  for (auto i = 0u; i < a.size(); ++i) {
    totalCost += costsMat[a[i]][b[i]];
  }
  return totalCost;
}

// Calculate the tier 2 similarity between the 2 mols.
double tier2Sim(const ROMol &mol1, const ROMol &mol2,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs2,
                const std::vector<unsigned int> &bondLabels1,
                const std::vector<unsigned int> &bondLabels2) {
  int vg1g2 = 0;
  int eg1g2 = 0;
  for (const auto &it1 : degSeqs1) {
    const auto &seq2 = degSeqs2.find(it1.first);
    if (seq2 != degSeqs2.end()) {
      vg1g2 += std::min(it1.second.size(), seq2->second.size());
      eg1g2 += getAssignmentScore(it1.second, seq2->second, bondLabels1,
                                  bondLabels2, mol1, mol2);
    }
  }
  eg1g2 /= 2;
  double sim = double((vg1g2 + eg1g2) * (vg1g2 + eg1g2)) /
               double((mol1.getNumAtoms() + mol1.getNumBonds()) *
                      (mol2.getNumAtoms() + mol2.getNumBonds()));
  return sim;
}

// make the line graph for the molecule, as an adjacency matrix.  Each
// row/column is a bond, with a connection between 2 bonds if they share an
// atom.  The adjacency matrix is 0 for no bond, the atomic number of the
// connecting atom otherwise.
void makeLineGraph(const ROMol &mol, std::vector<std::vector<int>> &adjMatrix) {
  adjMatrix = std::vector<std::vector<int>>(
      mol.getNumBonds(), std::vector<int>(mol.getNumBonds(), 0));
  for (const auto &a : mol.atoms()) {
    for (const auto &b1 : mol.atomBonds(a)) {
      for (const auto &b2 : mol.atomBonds(a)) {
        if (b1 != b2) {
          adjMatrix[b1->getIdx()][b2->getIdx()] = a->getAtomicNum();
        }
      }
    }
  }
}

// make sure that mol1_bond in mol1 and mol2_bond in mol2 are, if aromatic, in
// at least one ring that is the same.
bool checkAromaticRings(const ROMol &mol1,
                        std::vector<std::string> &mol1RingSmiles,
                        int mol1BondIdx, const ROMol &mol2,
                        std::vector<std::string> &mol2RingSmiles,
                        int mol2BondIdx) {
  auto mol1Bond = mol1.getBondWithIdx(mol1BondIdx);
  auto mol2Bond = mol2.getBondWithIdx(mol2BondIdx);
  if (!mol1Bond->getIsAromatic() || !mol2Bond->getIsAromatic()) {
    return true;
  }
  auto mol1BondRings = mol1.getRingInfo()->bondRings();
  auto mol2BondRings = mol2.getRingInfo()->bondRings();
  for (size_t i = 0u; i < mol1BondRings.size(); ++i) {
    if (std::find(mol1BondRings[i].begin(), mol1BondRings[i].end(),
                  mol1BondIdx) == mol1BondRings[i].end()) {
      continue;
    }
    for (size_t j = 0u; j < mol2BondRings.size(); ++j) {
      if (std::find(mol2BondRings[j].begin(), mol2BondRings[j].end(),
                    mol2BondIdx) == mol2BondRings[j].end()) {
        continue;
      }
      if (mol1RingSmiles[i] == mol2RingSmiles[j]) {
        return true;
      }
    }
  }

  return false;
}

// Extract the rings from the given molecule, both as mol objects and SMILES
// strings. The mol objects will have the original bond indices stored in the
// property ORIG_INDEX.
void extractRings(const ROMol &mol,
                  std::vector<std::unique_ptr<ROMol>> &molRings,
                  std::vector<std::string> &molRingSmiles) {
  auto molBondRings = mol.getRingInfo()->bondRings();
  for (size_t i = 0u; i < molBondRings.size(); ++i) {
    std::unique_ptr<RWMol> ringMol(new RWMol(mol));
    auto molAtomRings = mol.getRingInfo()->atomRings();
    std::vector<char> atomsInRing(mol.getNumAtoms(), 0);
    for (auto a : molAtomRings[i]) {
      atomsInRing[a] = 1;
    }
    for (auto ring_bond_idx : molBondRings[i]) {
      auto ringBond = ringMol->getBondWithIdx(ring_bond_idx);
      ringBond->setProp<int>("ORIG_INDEX", ringBond->getIdx());
    }
    ringMol->beginBatchEdit();
    for (auto b : ringMol->bonds()) {
      if (!b->hasProp("ORIG_INDEX)")) {
        if (!atomsInRing[b->getBeginAtomIdx()]) {
          ringMol->removeAtom(b->getBeginAtom());
        }
        if (!atomsInRing[b->getEndAtomIdx()]) {
          ringMol->removeAtom(b->getEndAtom());
        }
      }
    }
    ringMol->commitBatchEdit();
    molRings.push_back(std::unique_ptr<ROMol>(new ROMol(*ringMol)));
    molRingSmiles.push_back(MolToSmiles(*ringMol));
  }
}

bool checkRingMatchesRing(const ROMol &mol1, int mol1BondIdx, const ROMol &mol2,
                          int mol2BondIdx) {
  if (mol1.getRingInfo()->numBondRings(mol1BondIdx) &&
      !mol2.getRingInfo()->numBondRings(mol2BondIdx)) {
    return false;
  }
  if (!mol1.getRingInfo()->numBondRings(mol1BondIdx) &&
      mol2.getRingInfo()->numBondRings(mol2BondIdx)) {
    return false;
  }
  return true;
}

// Make the set of pairs of vertices, where they're a pair if the labels match.
void buildPairs(const ROMol &mol1, const std::vector<unsigned int> &vtxLabels1,
                const ROMol &mol2, const std::vector<unsigned int> &vtxLabels2,
                RascalOptions opts,
                std::vector<std::pair<int, int>> &vtxPairs) {
  std::vector<std::string> mol1RingSmiles, mol2RingSmiles;
  std::vector<std::unique_ptr<ROMol>> mol1Rings, mol2Rings;
  // For these purposes, it is correct that n1cccc1 and [nH]1cccc1 match - the
  // former would be from an N-substituted pyrrole, the latter from a plain one.
  static const std::regex reg(R"(\[([np])H\])");
  if (opts.completeAromaticRings) {
    extractRings(mol1, mol1Rings, mol1RingSmiles);
    for (auto &mrs : mol1RingSmiles) {
      mrs = std::regex_replace(mrs, reg, "$1");
    }
    extractRings(mol2, mol2Rings, mol2RingSmiles);
    for (auto &mrs : mol1RingSmiles) {
      mrs = std::regex_replace(mrs, reg, "$1");
    }
  }

  for (auto i = 0u; i < vtxLabels1.size(); ++i) {
    for (auto j = 0u; j < vtxLabels2.size(); ++j) {
      if (vtxLabels1[i] == vtxLabels2[j]) {
        if (opts.completeAromaticRings &&
            !checkAromaticRings(mol1, mol1RingSmiles, i, mol2, mol2RingSmiles,
                                j)) {
          continue;
        }
        if (opts.ringMatchesRingOnly &&
            !checkRingMatchesRing(mol1, i, mol2, j)) {
          continue;
        }
        vtxPairs.push_back(std::make_pair(i, j));
      }
    }
  }
}

// Make the modular product between the 2 graphs passed in.  Each node in the
// graph is a pair of vertices, one from the first graph, the other from the
// second, whose labels match.  Two vertices are connected in the modular
// product if either the 2 matching vertices in the 2 input vertices are
// connected by edges with the same label, or neither is connected.
void makeModularProduct(const ROMol &mol1,
                        const std::vector<std::vector<int>> &adjMatrix1,
                        const std::vector<unsigned int> &vtxLabels1,
                        const std::vector<std::vector<int>> &distMatrix1,
                        const ROMol &mol2,
                        const std::vector<std::vector<int>> &adjMatrix2,
                        const std::vector<unsigned int> &vtxLabels2,
                        const std::vector<std::vector<int>> &distMatrix2,
                        RascalOptions opts,
                        std::vector<std::pair<int, int>> &vtxPairs,
                        std::vector<std::vector<char>> &modProd) {
  buildPairs(mol1, vtxLabels1, mol2, vtxLabels2, opts, vtxPairs);
  if (vtxPairs.empty()) {
    // There was nothing in common at all.  But, what was the screening doing?
    modProd.clear();
    return;
  }

  modProd = std::vector<std::vector<char>>(
      vtxPairs.size(), std::vector<char>(vtxPairs.size(), 0));
  for (auto i = 0u; i < vtxPairs.size() - 1; ++i) {
    for (auto j = i + 1; j < vtxPairs.size(); ++j) {
      if (vtxPairs[i].first == vtxPairs[j].first ||
          vtxPairs[i].second == vtxPairs[j].second) {
        continue;
      }
      bool distsOk = true;
      if (opts.maxFragSeparation != -1) {
        if (std::abs(distMatrix1[vtxPairs[i].first][vtxPairs[j].first] -
                     distMatrix2[vtxPairs[i].second][vtxPairs[j].second]) >
            opts.maxFragSeparation) {
          distsOk = false;
        }
      }
      if (distsOk && adjMatrix1[vtxPairs[i].first][vtxPairs[j].first] ==
                         adjMatrix2[vtxPairs[i].second][vtxPairs[j].second]) {
        modProd[i][j] = modProd[j][i] = 1;
      }
    }
  }
}

// Calculate the lower bound on the size of the MCES.  This requires that mol1
// has more atoms than mol2 which is not checked.  Returns a minimum of 1.
int calcLowerBound(const ROMol &mol1, const ROMol &mol2, double simThresh) {
  std::set<int> mol1Atnos, mol2Atnos;
  for (const auto &a : mol1.atoms()) {
    mol1Atnos.insert(a->getAtomicNum());
  }
  for (const auto &a : mol2.atoms()) {
    mol2Atnos.insert(a->getAtomicNum());
  }
  int deltaVg1 = 0;
  for (auto mol1_atno : mol1Atnos) {
    if (mol2Atnos.find(mol1_atno) == mol2Atnos.end()) {
      ++deltaVg1;
    }
  }
  double lb = sqrt((mol1.getNumAtoms() + mol1.getNumBonds()) *
                   (mol2.getNumAtoms() + mol2.getNumBonds()));
  lb = lb * simThresh - mol1.getNumAtoms() + deltaVg1;
  int ilb = int(lb);
  if (ilb < 1) {
    ilb = 1;
  }
  return ilb;
}

void printClique(const std::vector<unsigned int> &clique,
                 const std::vector<std::pair<int, int>> &vtxPairs, bool swapped,
                 std::ostream &os) {
  os << "Clique : " << clique.size() << " :";
  for (auto mem : clique) {
    os << " " << mem;
  }
  os << std::endl;
  for (auto mem : clique) {
    if (swapped) {
      os << "{" << vtxPairs[mem].second << ", " << vtxPairs[mem].first << "},";
    } else {
      os << "{" << vtxPairs[mem].first << ", " << vtxPairs[mem].second << "},";
    }
  }
  std::cout << std::endl;
  std::cout << "mol 1 bonds : [";
  for (auto mem : clique) {
    if (swapped) {
      os << vtxPairs[mem].second << ", ";
    } else {
      os << vtxPairs[mem].first << ", ";
    }
  }
  std::cout << "]" << std::endl;
  std::cout << "mol 2 bonds : [";
  for (auto mem : clique) {
    if (swapped) {
      os << vtxPairs[mem].first << ", ";
    } else {
      os << vtxPairs[mem].second << ", ";
    }
  }
  std::cout << "]" << std::endl;
}

// if the clique involves a delta-y exchange, returns true.  Should only be
// called if it's a possibility.
bool deltaYInClique(const std::vector<unsigned int> &clique, const ROMol &mol1,
                    const ROMol &mol2,
                    const std::vector<std::pair<int, int>> &vtxPairs) {
  if (clique.size() < 3) {
    // there must be 3 bonds for a delta-y exchange, obs.
    return false;
  }
  // Map the clique onto the 2 molecules, counting the degrees of the atoms
  // if they are involved in the clique.  When sorted, they will be the same
  // if no delta-y exchange has occurred.
  std::vector<std::pair<int, int>> bondMatches;
  for (auto mem : clique) {
    bondMatches.push_back(
        std::make_pair(vtxPairs[mem].first, vtxPairs[mem].second));
  }
  std::vector<int> cliqueDegs1(mol1.getNumAtoms(), 0);
  std::vector<int> cliqueDegs2(mol2.getNumAtoms(), 0);
  for (const auto &bm : bondMatches) {
    const auto b1 = mol1.getBondWithIdx(bm.first);
    cliqueDegs1[b1->getBeginAtomIdx()]++;
    cliqueDegs1[b1->getEndAtomIdx()]++;
    const auto b2 = mol2.getBondWithIdx(bm.second);
    cliqueDegs2[b2->getBeginAtomIdx()]++;
    cliqueDegs2[b2->getEndAtomIdx()]++;
  }
  cliqueDegs1.erase(std::remove(cliqueDegs1.begin(), cliqueDegs1.end(), 0),
                    cliqueDegs1.end());
  std::sort(cliqueDegs1.begin(), cliqueDegs1.end());
  cliqueDegs2.erase(std::remove(cliqueDegs2.begin(), cliqueDegs2.end(), 0),
                    cliqueDegs2.end());
  std::sort(cliqueDegs2.begin(), cliqueDegs2.end());
  return cliqueDegs1 != cliqueDegs2;
}

// Return a molecule with the clique in it.  Each atom will have the property
// ORIG_INDEX giving its index in the original molecule.
RWMol *makeCliqueFrags(const ROMol &mol,
                       const std::vector<unsigned int> &clique,
                       const std::vector<std::pair<int, int>> &vtxPairs,
                       int pairNum) {
  auto *mol_frags = new RWMol(mol);
  std::vector<char> aInClique(mol.getNumAtoms(), 0);
  std::vector<char> bInClique(mol.getNumBonds(), 0);
  for (auto mem : clique) {
    const Bond *bond = nullptr;
    if (pairNum == 1) {
      bond = mol_frags->getBondWithIdx(vtxPairs[mem].first);
    } else {
      bond = mol_frags->getBondWithIdx(vtxPairs[mem].second);
    }
    bInClique[bond->getIdx()] = 1;
    aInClique[bond->getBeginAtomIdx()] = 1;
    bond->getBeginAtom()->setProp<int>("ORIG_INDEX", bond->getBeginAtomIdx());
    aInClique[bond->getEndAtomIdx()] = 1;
    bond->getEndAtom()->setProp<int>("ORIG_INDEX", bond->getEndAtomIdx());
  }
  mol_frags->beginBatchEdit();
  for (auto &a : mol_frags->atoms()) {
    if (!aInClique[a->getIdx()]) {
      mol_frags->removeAtom(a);
    }
  }
  for (auto &b : mol_frags->bonds()) {
    if (!bInClique[b->getIdx()]) {
      mol_frags->removeBond(b->getBeginAtomIdx(), b->getEndAtomIdx());
    }
  }
  mol_frags->commitBatchEdit();
  return mol_frags;
}

// Calculate the shortest bond distance between the 2 fragments in the molecule.
int minFragSeparation(const ROMol &mol, const ROMol &molFrags,
                      std::vector<int> &fragMapping, int frag1, int frag2) {
  auto extractFragAtoms = [&](int frag_num, std::vector<int> &frag_atoms) {
    for (size_t i = 0u; i < fragMapping.size(); ++i) {
      if (fragMapping[i] == frag_num) {
        int orig_idx = molFrags.getAtomWithIdx(i)->getProp<int>("ORIG_INDEX");
        frag_atoms.push_back(orig_idx);
      }
    }
  };
  std::vector<int> frag1Atoms, frag2Atoms;
  extractFragAtoms(frag1, frag1Atoms);
  extractFragAtoms(frag2, frag2Atoms);
  auto pathMatrix = MolOps::getDistanceMat(mol);
  double minDist = std::numeric_limits<double>::max();
  for (const auto &at1 : frag1Atoms) {
    for (const auto &at2 : frag2Atoms) {
      auto dist = pathMatrix[mol.getNumAtoms() * at1 + at2];
      if (dist < minDist) {
        minDist = dist;
      }
    }
  }
  return std::nearbyint(minDist);
}

// Assess the clique in terms of opts, returning true if it satisfies them all
bool cliqueOk(const std::vector<unsigned int> clique, const RascalOptions &opts,
              const ROMol &mol1, const ROMol &mol2,
              const std::vector<std::pair<int, int>> &vtxPairs) {
  std::unique_ptr<RWMol> mol1Frags, mol2Frags;
  std::vector<int> mol1FragMapping, mol2FragMapping;
  int numMol1Frags = 0, numMol2Frags = 0;

  auto buildFrags = [&]() -> void {
    if (mol1Frags) {
      return;
    }
    mol1Frags.reset(makeCliqueFrags(mol1, clique, vtxPairs, 1));
    mol2Frags.reset(makeCliqueFrags(mol2, clique, vtxPairs, 2));
    numMol1Frags = MolOps::getMolFrags(*mol1Frags, mol1FragMapping);
    numMol2Frags = MolOps::getMolFrags(*mol2Frags, mol2FragMapping);
  };

  if (opts.minFragSize > 0) {
    buildFrags();
    // only need to do it for mol1, as the fragments should match.
    for (int i = 0; i < numMol1Frags; ++i) {
      auto fragSize =
          std::count(mol1FragMapping.begin(), mol1FragMapping.end(), i);
      if (fragSize < opts.minFragSize) {
        return false;
      }
    }
  }

  return true;
}

// If this clique warrants it, update maxCliques.
void updateMaxClique(const std::vector<unsigned int> &clique, bool deltaYPoss,
                     const RascalOptions &opts, const ROMol &mol1,
                     const ROMol &mol2,
                     const std::vector<std::pair<int, int>> &vtxPairs,
                     std::vector<std::vector<unsigned int>> &maxCliques,
                     int &lower_bound) {
  if (!maxCliques.empty() && clique.size() < maxCliques.front().size()) {
    return;
  }
  bool didDeltaY =
      !deltaYPoss ? false : deltaYInClique(clique, mol1, mol2, vtxPairs);
  if (!didDeltaY) {
    if (maxCliques.empty()) {
      if (cliqueOk(clique, opts, mol1, mol2, vtxPairs)) {
        maxCliques.push_back((clique));
      }
    } else {
      bool goodClique = false, didCliqueOk = false;
      if (clique.size() > maxCliques.front().size()) {
        goodClique = cliqueOk(clique, opts, mol1, mol2, vtxPairs);
        didCliqueOk = true;
        if (goodClique) {
          maxCliques.clear();
        }
      }
      if (!didCliqueOk) {
        goodClique = cliqueOk(clique, opts, mol1, mol2, vtxPairs);
      }
      if (goodClique &&
          (maxCliques.empty() || clique.size() == maxCliques.front().size())) {
        maxCliques.push_back(clique);
      }
    }
    if (!maxCliques.empty() && maxCliques.front().size() > lower_bound) {
      lower_bound = maxCliques.front().size();
    }
  }
}

// If the current time is beyond the timeout limit, throws a TimedOutException.
void checkTimeout(
    std::chrono::time_point<std::chrono::high_resolution_clock> &startTime,
    const RascalOptions &opts, const std::vector<unsigned int> &clique,
    std::vector<std::vector<unsigned int>> &maxCliques,
    unsigned long long &numSteps) {
  ++numSteps;
  if (numSteps == 100) {
    // This clock is very convenient, but seems quite expensive.  Calling it
    // every step added 10% to the runtime.
    auto currTime = std::chrono::high_resolution_clock::now();
    auto runTime =
        std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime)
            .count();
    if (runTime > opts.timeout) {
      if (maxCliques.empty()) {
        maxCliques.push_back(clique);
      } else {
        if (clique.size() > maxCliques.front().size()) {
          maxCliques.clear();
        }
        if (clique.size() >= maxCliques.front().size()) {
          maxCliques.push_back(clique);
        }
      }
      throw TimedOutException(runTime, maxCliques);
    }
    numSteps = 0ULL;
  }
}

bool equivalentRootAlreadyDone(unsigned int rootVtx,
                               const std::vector<std::pair<int, int>> &vtxPairs,
                               const std::vector<int> &equivBonds1,
                               const std::vector<int> &equivBonds2,
                               std::set<std::pair<int, int>> &rootClasses) {
  std::pair<int, int> newClasses{equivBonds1[vtxPairs[rootVtx].first],
                                 equivBonds2[vtxPairs[rootVtx].second]};
  if (newClasses.first == -1) {
    return false;
  }
  if (!rootClasses.empty() &&
      rootClasses.find(newClasses) != rootClasses.end()) {
    return true;
  }
  rootClasses.insert(newClasses);
  return false;
}

// There are some simple substructures for which equivalent bond pruning isn't
// allowed.
bool checkEquivalentsAllowed(const ROMol &mol) {
  const static std::vector<ROMol *> notStructs{
      SmartsToMol("*~*"), SmartsToMol("*~*1~*~*~1"),
      SmartsToMol("*12~*~*~2~*~1"), SmartsToMol("*14~*(~*~2~3~4)~*~2~*~3~1")};
  const static std::vector<std::pair<int, int>> notStats{
      {2, 1}, {4, 4}, {4, 5}, {5, 8}};
  MatchVectType dontCare;
  for (size_t i = 0; i < notStructs.size(); ++i) {
    if (mol.getNumAtoms() == notStats[i].first &&
        mol.getNumBonds() == notStats[i].second &&
        SubstructMatch(mol, *notStructs[i], dontCare)) {
      return false;
    }
  }
  return true;
}

void explore_partitions(
    RascalStartPoint &starter,
    std::chrono::time_point<std::chrono::high_resolution_clock> &startTime,
    const RascalOptions &opts,
    std::vector<std::vector<unsigned int>> &maxCliques) {
  unsigned long long numSteps = 0ULL;
  std::vector<std::shared_ptr<PartitionSet>> parts(1, starter.d_partSet);
  std::vector<unsigned int> clique;
  std::set<std::pair<int, int>> rootClasses;
  bool canDoEquivs = false;
  if (opts.doEquivBondPruning) {
    canDoEquivs = checkEquivalentsAllowed(*starter.d_mol1) &&
                  checkEquivalentsAllowed(*starter.d_mol2);
  }
  while (!parts.empty()) {
    if (opts.timeout != -1) {
      checkTimeout(startTime, opts, clique, maxCliques, numSteps);
    }
    auto part = parts.back();
    bool goDeeper = false;
    bool backtrack = false;
    if (opts.allBestMCESs) {
      if (clique.size() + part->num_parts() < starter.d_lowerBound) {
        backtrack = true;
      }
    } else {
      if (clique.size() + part->num_parts() <= starter.d_lowerBound) {
        backtrack = true;
      }
    }
    if (!backtrack) {
      if (opts.allBestMCESs) {
        goDeeper = clique.size() + part->upper_bound() >= starter.d_lowerBound;
      } else {
        goDeeper = clique.size() + part->upper_bound() > starter.d_lowerBound;
      }
      if (goDeeper) {
        if (!part->is_empty()) {
          std::shared_ptr<PartitionSet> nextPart(new PartitionSet(*part));
          clique.push_back(nextPart->pop_last_vertex());
          if (clique.size() == 1 && canDoEquivs &&
              equivalentRootAlreadyDone(clique.front(), starter.d_vtxPairs,
                                        starter.d_equivBonds1,
                                        starter.d_equivBonds2, rootClasses)) {
            clique.pop_back();
            backtrack = true;
          } else {
            nextPart->prune_vertices(clique.back());
            updateMaxClique(clique, starter.d_deltaYPoss, opts, *starter.d_mol1,
                            *starter.d_mol2, starter.d_vtxPairs, maxCliques,
                            starter.d_lowerBound);
            parts.push_back(nextPart);
          }
        } else {
          backtrack = true;
        }
      } else {
        backtrack = true;
      }
    }
    if (backtrack || (!parts.empty() && parts.back()->is_empty())) {
      while (!parts.empty()) {
        if (parts.back()->is_empty()) {
          parts.pop_back();
          if (!clique.empty()) {
            clique.pop_back();
          }
        } else {
          parts.back()->pop_last_vertex();
          if (!parts.back()->is_empty()) {
            break;
          }
        }
      }
    }
    if (parts.empty()) {
      break;
    }
  }
}

bool deltaYExchangePossible(const ROMol &mol1, const ROMol &mol2) {
  // A Delta-y exchange is an incorrect match when a cyclopropyl ring (the
  // delta) is matched to a C(C)(C) group (the y) because they both have
  // isomorphic line graphs.  This checks to see if that's something we need to
  // worry about for these molecules.
  const static std::unique_ptr<ROMol> delta(SmartsToMol("C1CC1"));
  const static std::unique_ptr<ROMol> y(SmartsToMol("C(C)C"));
  MatchVectType dontCare;
  return (SubstructMatch(mol1, *delta, dontCare) &&
          SubstructMatch(mol2, *y, dontCare)) ||
         (SubstructMatch(mol2, *delta, dontCare) &&
          SubstructMatch(mol1, *y, dontCare));
}

void findEquivalentBonds(const ROMol &mol, std::vector<int> &equivBonds) {
  equivBonds = std::vector<int>(mol.getNumBonds(), -1);
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  bool breakTies = false;
  Canon::rankMolAtoms(mol, ranks, breakTies);
  int nextClass = 0;
  for (const auto &b1 : mol.bonds()) {
    for (const auto &b2 : mol.bonds()) {
      if (b1->getIdx() != b2->getIdx()) {
        if ((ranks[b1->getBeginAtomIdx()] == ranks[b2->getBeginAtomIdx()] &&
             ranks[b1->getEndAtomIdx()] == ranks[b2->getEndAtomIdx()]) ||
            (ranks[b1->getBeginAtomIdx()] == ranks[b2->getEndAtomIdx()] &&
             ranks[b1->getEndAtomIdx()] == ranks[b2->getBeginAtomIdx()])) {
          if (equivBonds[b1->getIdx()] == -1 &&
              equivBonds[b2->getIdx()] == -1) {
            equivBonds[b1->getIdx()] = nextClass;
            equivBonds[b2->getIdx()] = nextClass;
            ++nextClass;
          } else if (equivBonds[b1->getIdx()] == -1) {
            equivBonds[b1->getIdx()] = equivBonds[b2->getIdx()];
          } else if (equivBonds[b2->getIdx()] == -1) {
            equivBonds[b2->getIdx()] = equivBonds[b1->getIdx()];
          }
        }
      }
    }
  }
}

// Use the Floyd-Warshall algorithm to compute the distance matrix from the
// adjacency matrix.
void calcDistMatrix(const std::vector<std::vector<int>> &adjMatrix,
                    std::vector<std::vector<int>> &distMatrix) {
  distMatrix = std::vector<std::vector<int>>(
      adjMatrix.size(),
      std::vector<int>(adjMatrix.size(), adjMatrix.size() + 1));
  for (size_t i = 0u; i < adjMatrix.size(); ++i) {
    distMatrix[i][i] = 0;
    for (size_t j = 0; j < adjMatrix.size(); ++j) {
      if (adjMatrix[i][j]) {
        distMatrix[i][j] = 1;
        distMatrix[j][i] = 1;
      }
    }
  }
  for (size_t i = 0u; i < adjMatrix.size(); ++i) {
    for (size_t j = 0u; j < adjMatrix.size(); ++j) {
      for (size_t k = 0u; k < adjMatrix.size(); ++k) {
        if (distMatrix[i][j] > distMatrix[i][k] + distMatrix[k][j]) {
          distMatrix[i][j] = distMatrix[i][k] + distMatrix[k][j];
        }
      }
    }
  }
}

RascalStartPoint makeInitialPartitionSet(const ROMol *mol1, const ROMol *mol2,
                                         const RascalOptions &opts) {
  RascalStartPoint starter;
  if (mol1->getNumAtoms() <= mol2->getNumAtoms()) {
    starter.d_swapped = false;
    starter.d_mol1 = mol1;
    starter.d_mol2 = mol2;
  } else {
    starter.d_swapped = true;
    starter.d_mol1 = mol2;
    starter.d_mol2 = mol1;
  }
  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  starter.d_tier1Sim =
      tier_1_sim(*starter.d_mol1, *starter.d_mol2, degSeqs1, degSeqs2);
  if (starter.d_tier1Sim < opts.similarityThreshold) {
    return starter;
  }
  std::vector<unsigned int> bondLabels1, bondLabels2;
  get_bond_labels(*starter.d_mol1, *starter.d_mol2, opts, bondLabels1,
                  bondLabels2);
  starter.d_tier2Sim = tier2Sim(*starter.d_mol1, *starter.d_mol2, degSeqs1,
                                degSeqs2, bondLabels1, bondLabels2);
  if (starter.d_tier2Sim < opts.similarityThreshold) {
    return starter;
  }

  // Get the line graphs for the two molecules as adjacency matrices.
  makeLineGraph(*starter.d_mol1, starter.d_adjMatrix1);
  makeLineGraph(*starter.d_mol2, starter.d_adjMatrix2);

  std::vector<std::vector<int>> dist_mat1, dist_mat2;
  if (opts.maxFragSeparation > -1) {
    calcDistMatrix(starter.d_adjMatrix1, dist_mat1);
    calcDistMatrix(starter.d_adjMatrix2, dist_mat2);
  }

  // pairs are vertices in the 2 line graphs that are the same type.
  // mod_prod is the modular product/correspondence graph of the two
  // line graphs.
  makeModularProduct(*starter.d_mol1, starter.d_adjMatrix1, bondLabels1,
                     dist_mat1, *starter.d_mol2, starter.d_adjMatrix2,
                     bondLabels2, dist_mat2, opts, starter.d_vtxPairs,
                     starter.d_modProd);
  if (starter.d_vtxPairs.empty()) {
    return starter;
  }
  starter.d_lowerBound = calcLowerBound(*starter.d_mol1, *starter.d_mol2,
                                        opts.similarityThreshold);

  starter.d_partSet.reset(new PartitionSet(starter.d_modProd,
                                           starter.d_vtxPairs, bondLabels1,
                                           bondLabels2, starter.d_lowerBound));
  starter.d_deltaYPoss =
      deltaYExchangePossible(*starter.d_mol1, *starter.d_mol2);

  if (opts.doEquivBondPruning) {
    // if equiv_bonds1[i] and equiv_bonds1[j] are equal, the bonds are
    // equivalent.
    findEquivalentBonds(*starter.d_mol1, starter.d_equivBonds1);
    findEquivalentBonds(*starter.d_mol2, starter.d_equivBonds2);
  } else {
    starter.d_equivBonds1 = std::vector<int>(starter.d_mol1->getNumBonds(), -1);
    starter.d_equivBonds2 = std::vector<int>(starter.d_mol2->getNumBonds(), -1);
  }

  return starter;
}

std::vector<RascalResult> findMces(RascalStartPoint &starter,
                                   RascalOptions opts) {
  // opts.singleLargestFrag requires opts.allBestMCESs
  bool orig_allBestMCESs = opts.allBestMCESs;
  if (opts.singleLargestFrag) {
    opts.allBestMCESs = true;
  }
  std::vector<unsigned int> clique;
  std::vector<std::vector<unsigned int>> maxCliques;
  auto start_time = std::chrono::high_resolution_clock::now();
  bool timed_out = false;
  try {
    explore_partitions(starter, start_time, opts, maxCliques);
    auto curr_time = std::chrono::high_resolution_clock::now();
    auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                        curr_time - start_time)
                        .count();
  } catch (TimedOutException &e) {
    std::cout << e.what() << std::endl;
    maxCliques = e.d_cliques;
    timed_out = true;
  }
  opts.allBestMCESs = orig_allBestMCESs;

  std::vector<RascalResult> results;
  for (const auto &c : maxCliques) {
    results.push_back(RascalResult(
        *starter.d_mol1, *starter.d_mol2, starter.d_adjMatrix1,
        starter.d_adjMatrix2, c, starter.d_vtxPairs, timed_out,
        starter.d_swapped, opts.exactChirality, opts.maxFragSeparation));
    if (opts.singleLargestFrag) {
      results.back().largestFragOnly();
    }
  }
  std::sort(results.begin(), results.end(), resultSort);
  return results;
}

// calculate the RASCAL MCES between the 2 molecules, provided it is within the
// similarity threshold given.
std::vector<RascalResult> rascalMces(const ROMol &mol1, const ROMol &mol2,
                                     RascalOptions opts) {
  auto starter = makeInitialPartitionSet(&mol1, &mol2, opts);
  //  std::cout << "tier 1 and 2 similarities : " << starter.d_tier1Sim << " and
  //  "
  //            << starter.d_tier2Sim << std::endl;
  if (!starter.d_partSet) {
    return std::vector<RascalResult>();
  }
  //  std::cout << "clique size bounds : " << starter.d_lowerBound << " to "
  //            << starter.d_partSet->upper_bound() << std::endl;
  //  std::cout << "tier 1 and 2 similarities : " << starter.d_tier1Sim << " and
  //  "
  //            << starter.d_tier2Sim << std::endl;

  auto results = findMces(starter, opts);
  if (!opts.allBestMCESs && results.size() > 1) {
    results.erase(results.begin() + 1, results.end());
  }
  return results;
}

}  // namespace RascalMCES
}  // namespace RDKit