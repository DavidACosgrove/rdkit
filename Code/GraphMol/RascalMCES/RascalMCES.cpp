//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains the driving functions for the determination of
// the Maximum Common Edge Substructure (MCES) between 2 molecules.
// It uses the RASCAL algorithm of John Raymond:
// RASCAL: Calculation of Graph Similarity using Maximum Common
// Edge Subgraphs, John W. Raymond, Eleanor J. Gardiner, Peter Willett
// 'The Computer Journal', 45, 631-644 (2002).
// https://eprints.whiterose.ac.uk/3568/1/willets3.pdf

#include <chrono>
#include <iostream>
#include <map>
#include <regex>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include <GraphMol/atomic_data.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/PartitionSet.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
#include <GraphMol/RascalMCES/RascalDetails.h>

namespace RDKit {
namespace RascalMCES {

class TimedOutException : public std::exception {
 public:
  TimedOutException(long long int run_time,
                    std::vector<std::vector<unsigned int>> &bestCliques)
      : d_cliques(bestCliques) {
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
int lapMaximize(const std::vector<std::vector<int>> &costsMat,
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
  unsigned int d_lowerBound{0};
  // a Delta-Y exchange requires extra treatment.  They're rare, though.
  bool d_deltaYPoss{false};

  // Some parts require mol2 to be the larger molecule.  This records if they
  // have been swapped with respect to the input molecules.
  bool d_swapped{false};
  // Pointers to copies of the input molecules, swapped if necessary,
  // so that d_mol1 is always the smaller molecule.
  std::unique_ptr<ROMol> d_mol1;
  std::unique_ptr<ROMol> d_mol2;

  std::vector<std::vector<const boost::dynamic_bitset<> *>> d_adjMatrix1,
      d_adjMatrix2;
  std::vector<std::pair<int, int>> d_vtxPairs;
  std::vector<boost::dynamic_bitset<>> d_modProd;

  // We might need to know which bonds are symmetrical equivalent.
  std::vector<int> d_equivBonds1, d_equivBonds2;

  // The initial partition set.  May be empty if the thresholds weren't met.
  std::shared_ptr<PartitionSet> d_partSet;
};

// Get the sorted degree sequences for the molecule, one sequence for each
// atomic number in the molecule.  Each element in the degree sequence is
// the degree of the atom and its index.
void sortedDegreeSeqs(
    const ROMol &mol,
    std::map<int, std::vector<std::pair<int, int>>> &degSeqs) {
  for (const auto &a : mol.atoms()) {
    degSeqs[a->getAtomicNum()].push_back(
        std::make_pair(a->getDegree(), a->getIdx()));
  }
  for (auto &it : degSeqs) {
    std::sort(it.second.begin(), it.second.end(),
              [](const std::pair<int, int> &p1, const std::pair<int, int> &p2)
                  -> bool { return p1.first > p2.first; });
  }
}

void sortedDegreeSeqs(
    const ROMol &mol, const std::vector<boost::dynamic_bitset<>> &atomLabels,
    std::map<int, std::vector<std::pair<int, int>>> &degSeqs) {
  for (const auto &a : mol.atoms()) {
    const auto &al = atomLabels[a->getIdx()];
    for (unsigned int i = 0; i < al.size(); i++) {
      if (al[i]) {
        degSeqs[i].push_back(std::make_pair(a->getDegree(), a->getIdx()));
      }
    }
  }
  for (auto &it : degSeqs) {
    std::sort(it.second.begin(), it.second.end(),
              [](const std::pair<int, int> &p1, const std::pair<int, int> &p2)
                  -> bool { return p1.first > p2.first; });
  }
}

bool bitsetLess(const boost::dynamic_bitset<> &bs1,
                const boost::dynamic_bitset<> &bs2) {
  // bs1 is less than bs2 if the first bit set in bs1 for which the
  // corresponding bit in bs2 isn't set has a lower index than the first
  // bit set in bs2.
  PRECONDITION(bs1.size() == bs2.size(), "bitsets different sizes.");
  for (size_t i = 0; i < bs1.size(); ++i) {
    if (bs1[i] && !bs2[i]) {
      return true;
    }
    if (bs2[i] && !bs1[i]) {
      return false;
    }
  }
  return false;
}

// Find the number of bonds incident on atom i that match a bond incident
// on atom j.
unsigned int calcCost(
    const std::vector<
        std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>>> &atomiBLs,
    const std::vector<
        std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>>> &atomjBLs,
    unsigned int atomiDegree, unsigned int atomjDegree) {
  std::vector<std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>>>
      uniqAtomiBLs(atomiBLs);
  std::sort(uniqAtomiBLs.begin(), uniqAtomiBLs.end(),
            [](const auto &p1, const auto &p2) -> bool {
              // The bitsets should be the same size.
              if (p1.first == p2.first) {
                return bitsetLess(p1.second, p2.second);
              }
              return bitsetLess(p1.first, p2.first);
            });
  uniqAtomiBLs.erase(std::unique(uniqAtomiBLs.begin(), uniqAtomiBLs.end()),
                     uniqAtomiBLs.end());

  auto countMatches =
      [](const std::vector<
             std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>>> &vec,
         const std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>> &q)
      -> unsigned int {
    unsigned int count = 0;
    for (const auto &p : vec) {
      if ((p.first & q.first).count() && (p.second & q.second).count()) {
        ++count;
      }
    }
    return count;
  };

  unsigned int cost = 0;
  for (const auto &uai : uniqAtomiBLs) {
    int numAtomi = countMatches(atomiBLs, uai);
    int numAtomj = countMatches(atomjBLs, uai);
    cost += std::min(numAtomi, numAtomj);
  }
  return std::min({cost, atomiDegree, atomjDegree});
}

void assignCosts(const std::vector<std::pair<int, int>> &atomDegrees1,
                 const std::vector<std::pair<int, int>> &atomDegrees2,
                 const std::vector<boost::dynamic_bitset<>> atomLabels1,
                 const std::vector<boost::dynamic_bitset<>> &atomLabels2,
                 const std::vector<boost::dynamic_bitset<>> &bondLabels1,
                 const std::vector<boost::dynamic_bitset<>> &bondLabels2,
                 const ROMol &mol1, const ROMol &mol2,
                 std::vector<std::vector<int>> &costsMat) {
  // For each pair of atoms in mol1 and mol2, find the number of incident
  // matching bonds with each atom that match.
  std::vector<std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<>>>
      atomiBLs, atomjBLs;
  for (auto i = 0u; i < atomDegrees1.size(); ++i) {
    atomiBLs.clear();
    const auto atomi = mol1.getAtomWithIdx(atomDegrees1[i].second);
    for (const auto b : mol1.atomBonds(atomi)) {
      atomiBLs.push_back(
          std::make_pair(bondLabels1[b->getIdx()],
                         atomLabels1[b->getOtherAtomIdx(atomi->getIdx())]));
    }
    for (auto j = 0u; j < atomDegrees2.size(); ++j) {
      atomjBLs.clear();
      const auto atomj = mol2.getAtomWithIdx(atomDegrees2[j].second);
      for (const auto b : mol2.atomBonds(atomj)) {
        atomjBLs.push_back(
            std::make_pair(bondLabels2[b->getIdx()],
                           atomLabels2[b->getOtherAtomIdx(atomj->getIdx())]));
      }
      costsMat[i][j] =
          calcCost(atomiBLs, atomjBLs, atomi->getDegree(), atomj->getDegree());
    }
  }
}

// Return the assignment score for the best match of the atoms and bonds in mol1
// to the atoms and bonds in mol2 for this vertex label (the atomDegrees[12]).
int getAssignmentScore(const std::vector<std::pair<int, int>> &atomDegrees1,
                       const std::vector<std::pair<int, int>> &atomDegrees2,
                       const std::vector<boost::dynamic_bitset<>> atomLabels1,
                       const std::vector<boost::dynamic_bitset<>> &atomLabels2,
                       const std::vector<boost::dynamic_bitset<>> &bondLabels1,
                       const std::vector<boost::dynamic_bitset<>> &bondLabels2,
                       const ROMol &mol1, const ROMol &mol2) {
  constexpr int bigScore(9999);
  constexpr size_t unassignedValue(99999999);
  std::vector<std::vector<int>> costsMat(
      atomDegrees1.size(), std::vector<int>(atomDegrees2.size(), bigScore));
  assignCosts(atomDegrees1, atomDegrees2, atomLabels1, atomLabels2, bondLabels1,
              bondLabels2, mol1, mol2, costsMat);
  std::vector<size_t> a(std::min(atomDegrees1.size(), atomDegrees2.size()),
                        unassignedValue);
  std::vector<size_t> b(std::min(atomDegrees1.size(), atomDegrees2.size()),
                        unassignedValue);
  int retVal = lapMaximize(costsMat, a, b);
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

namespace details {

double tier1Sim(const ROMol &mol1, const ROMol &mol2,
                const std::vector<boost::dynamic_bitset<>> &atomLabels1,
                const std::vector<boost::dynamic_bitset<>> &atomLabels2,
                std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                std::map<int, std::vector<std::pair<int, int>>> &degSeqs2) {
  sortedDegreeSeqs(mol1, atomLabels1, degSeqs1);
  sortedDegreeSeqs(mol2, atomLabels2, degSeqs2);
  unsigned int vg1g2 = 0;
  unsigned int eg1g2 = 0;
  for (const auto &it1 : degSeqs1) {
    const auto &seq2 = degSeqs2.find(it1.first);
    if (seq2 != degSeqs2.end()) {
      vg1g2 += std::min(it1.second.size(), seq2->second.size());
      auto numToDo = std::min(it1.second.size(), seq2->second.size());
      for (auto i = 0U; i < numToDo; ++i) {
        eg1g2 += std::min(it1.second[i].first, seq2->second[i].first);
      }
    }
  }
  eg1g2 /= 2;

  // An atom can be in more than 1 label class, but the number of atoms
  // and bonds in the MCES clearly can't be more than the
  // number of atoms and bonds in the smaller molecule.
  vg1g2 = std::min({vg1g2, mol1.getNumAtoms(), mol2.getNumAtoms()});
  eg1g2 = std::min({eg1g2, mol1.getNumBonds(), mol2.getNumBonds()});
  double sim = double((vg1g2 + eg1g2) * (vg1g2 + eg1g2)) /
               double((mol1.getNumAtoms() + mol1.getNumBonds()) *
                      (mol2.getNumAtoms() + mol2.getNumBonds()));

  return sim;
}

double tier2Sim(const ROMol &mol1, const ROMol &mol2,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
                const std::map<int, std::vector<std::pair<int, int>>> &degSeqs2,
                const std::vector<boost::dynamic_bitset<>> &atomLabels1,
                const std::vector<boost::dynamic_bitset<>> &atomLabels2,
                const std::vector<boost::dynamic_bitset<>> &bondLabels1,
                const std::vector<boost::dynamic_bitset<>> &bondLabels2) {
  unsigned int vg1g2 = 0;
  unsigned int eg1g2 = 0;
  for (const auto &seq1 : degSeqs1) {
    if (const auto &seq2 = degSeqs2.find(seq1.first); seq2 != degSeqs2.end()) {
      vg1g2 += std::min(seq1.second.size(), seq2->second.size());
      eg1g2 +=
          getAssignmentScore(seq1.second, seq2->second, atomLabels1,
                             atomLabels2, bondLabels1, bondLabels2, mol1, mol2);
    }
  }
  eg1g2 /= 2;
  vg1g2 = std::min({vg1g2, mol1.getNumAtoms(), mol2.getNumAtoms()});
  eg1g2 = std::min({eg1g2, mol1.getNumBonds(), mol2.getNumBonds()});
  double sim = double((vg1g2 + eg1g2) * (vg1g2 + eg1g2)) /
               double((mol1.getNumAtoms() + mol1.getNumBonds()) *
                      (mol2.getNumAtoms() + mol2.getNumBonds()));
  return sim;
}

}  // namespace details

// make the line graph for the molecule, as an adjacency matrix.  Each
// row/column is a bond, with a connection between 2 bonds if they share an
// atom.  The adjacency matrix is 0 for no bond, the atomic number of the
// connecting atom otherwise.
void makeLineGraph(
    const ROMol &mol, const std::vector<boost::dynamic_bitset<>> &atomLabels,
    std::vector<std::vector<const boost::dynamic_bitset<> *>> &adjMatrix) {
  adjMatrix = std::vector<std::vector<const boost::dynamic_bitset<> *>>(
      mol.getNumBonds(),
      std::vector<const boost::dynamic_bitset<> *>(mol.getNumBonds(), 0));
  for (const auto &a : mol.atoms()) {
    for (const auto &b1 : mol.atomBonds(a)) {
      for (const auto &b2 : mol.atomBonds(a)) {
        if (b1 != b2) {
          adjMatrix[b1->getIdx()][b2->getIdx()] = &atomLabels[a->getIdx()];
        }
      }
    }
  }
}

namespace {
// Take the string holding the _smilesAtomOutputOrder and convert it to
// a vector of unsigned ints.
std::vector<unsigned int> orderStringToInts(const std::string &str) {
  std::string sht(str.substr(1, str.size() - 2));
  std::vector<unsigned int> res;
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(",");
  tokenizer tokens{sht, sep};
  for (const auto &token : tokens) {
    res.push_back(std::stoul(token));
  }
  return res;
}

}  // namespace
// make sure that mol1_bond in mol1 and mol2_bond in mol2 are, if aromatic,
// in at least one ring that is the same.
bool checkAromaticRings(const ROMol &mol1,
                        std::vector<std::string> &mol1RingSmiles,
                        const std::vector<std::unique_ptr<ROMol>> &mol1Rings,
                        const std::vector<boost::dynamic_bitset<>> &vtxLabels1,
                        const std::vector<boost::dynamic_bitset<>> &edgeLabels1,
                        int mol1BondIdx, const ROMol &mol2,
                        std::vector<std::string> &mol2RingSmiles,
                        const std::vector<std::unique_ptr<ROMol>> &mol2Rings,
                        const std::vector<boost::dynamic_bitset<>> &vtxLabels2,
                        const std::vector<boost::dynamic_bitset<>> &edgeLabels2,
                        int mol2BondIdx) {
  auto mol1Bond = mol1.getBondWithIdx(mol1BondIdx);
  auto mol2Bond = mol2.getBondWithIdx(mol2BondIdx);
  if (!mol1Bond->getIsAromatic() || !mol2Bond->getIsAromatic()) {
    return true;
  }

  // If neither bond was in a ring, but they were marked aromatic, then
  // the two mols are fragments so it's ok to match these bonds.
  const auto &mol1BondRings = mol1.getRingInfo()->bondRings();
  const auto &mol2BondRings = mol2.getRingInfo()->bondRings();
  bool mol1BondInRing = mol1.getRingInfo()->numBondRings(mol1BondIdx);
  bool mol2BondInRing = mol2.getRingInfo()->numBondRings(mol2BondIdx);
  if (!mol1BondInRing && !mol2BondInRing) {
    return true;
  }

  for (size_t i = 0u; i < mol1BondRings.size(); ++i) {
    if (std::find(mol1BondRings[i].begin(), mol1BondRings[i].end(),
                  mol1BondIdx) == mol1BondRings[i].end()) {
      continue;
    }
    for (size_t j = 0u; j < mol2BondRings.size(); ++j) {
      if (mol1Rings[i]->getNumAtoms() != mol2Rings[j]->getNumAtoms()) {
        continue;
      }
      if (std::find(mol2BondRings[j].begin(), mol2BondRings[j].end(),
                    mol2BondIdx) == mol2BondRings[j].end()) {
        continue;
      }
      if (mol1RingSmiles[i] == mol2RingSmiles[j]) {
        return true;
      }
      if (vtxLabels1.front().size() > elementNames.size() + 1) {
        bool match = true;
        // If the atom bitstring is bigger than the number of elements
        // RDKit knows about, check the case where 2 atoms have different
        // atomic numbers but have otherwise been deemed equivalent.
        // If the input contained equivalentAtoms or equivalentBonds the
        // SMILES won't necessarily match (e.g. with c1ccccc1 and c1cccnc1
        // and [c,n] for equivalentAtoms, so we need to check the vertex
        // labels.
        auto mol1AtomOrder =
            orderStringToInts(mol1Rings[i]->getProp<std::string>(
                common_properties::_smilesAtomOutputOrder));
        auto mol2AtomOrder =
            orderStringToInts(mol2Rings[j]->getProp<std::string>(
                common_properties::_smilesAtomOutputOrder));
        for (size_t k = 0; k < mol1AtomOrder.size(); ++k) {
          auto origMol1Atom = mol1Rings[i]
                                  ->getAtomWithIdx(mol1AtomOrder[k])
                                  ->getProp<int>("ORIG_INDEX");
          auto origMol2Atom = mol2Rings[j]
                                  ->getAtomWithIdx(mol2AtomOrder[k])
                                  ->getProp<int>("ORIG_INDEX");
          if (!(vtxLabels1[origMol1Atom] & vtxLabels2[origMol2Atom]).count()) {
            match = false;
            break;
          }
        }
        if (match) {
          auto mol1BondOrder =
              orderStringToInts(mol1Rings[i]->getProp<std::string>(
                  common_properties::_smilesBondOutputOrder));
          auto mol2BondOrder =
              orderStringToInts(mol2Rings[j]->getProp<std::string>(
                  common_properties::_smilesBondOutputOrder));
          for (size_t k = 0; k < mol1BondOrder.size(); ++k) {
            auto origMol1Bond = mol1Rings[i]
                                    ->getBondWithIdx(mol1BondOrder[k])
                                    ->getProp<int>("ORIG_INDEX");
            auto origMol2bond = mol2Rings[j]
                                    ->getBondWithIdx(mol2BondOrder[k])
                                    ->getProp<int>("ORIG_INDEX");
            if (!(edgeLabels1[origMol1Bond] & edgeLabels2[origMol2bond])
                     .count()) {
              match = false;
              break;
            }
          }
        }
        if (match) {
          return true;
        }
      }
    }
  }
  return false;
}

// Extract the rings from the given molecule, both as mol objects and SMILES
// strings. The mol objects will have the original atom and bond indices
// stored in the property ORIG_INDEX.
void extractRings(const ROMol &mol,
                  std::vector<std::unique_ptr<ROMol>> &molRings,
                  std::vector<std::string> &molRingSmiles) {
  const auto &molBondRings = mol.getRingInfo()->bondRings();
  for (size_t i = 0u; i < molBondRings.size(); ++i) {
    std::unique_ptr<RWMol> ringMol(new RWMol(mol));
    const auto &molAtomRings = mol.getRingInfo()->atomRings();
    boost::dynamic_bitset<> atomsInRing(mol.getNumAtoms());
    for (auto a : molAtomRings[i]) {
      atomsInRing.set(a);
    }
    for (auto ringBondIdx : molBondRings[i]) {
      auto ringBond = ringMol->getBondWithIdx(ringBondIdx);
      ringBond->setProp<int>("ORIG_INDEX", ringBond->getIdx());
    }
    ringMol->beginBatchEdit();
    for (auto b : ringMol->bonds()) {
      if (!b->hasProp("ORIG_INDEX)")) {
        if (!atomsInRing[b->getBeginAtomIdx()]) {
          ringMol->removeAtom(b->getBeginAtom());
        } else {
          b->getBeginAtom()->setProp<int>("ORIG_INDEX",
                                          b->getBeginAtom()->getIdx());
        }
        if (!atomsInRing[b->getEndAtomIdx()]) {
          ringMol->removeAtom(b->getEndAtom());
        } else {
          b->getEndAtom()->setProp<int>("ORIG_INDEX",
                                        b->getEndAtom()->getIdx());
        }
      }
    }
    ringMol->commitBatchEdit();
    molRingSmiles.push_back(MolToSmiles(*ringMol));
    molRings.push_back(std::move(ringMol));
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

bool bondsMatch(const Bond &bond1, const Bond &bond2,
                const std::vector<boost::dynamic_bitset<>> &vtxLabels1,
                const std::vector<boost::dynamic_bitset<>> &edgeLabels1,
                const std::vector<boost::dynamic_bitset<>> &vtxLabels2,
                const std::vector<boost::dynamic_bitset<>> &edgeLabels2,
                bool ignoreAtomAromaticity) {
  if ((edgeLabels1[bond1.getIdx()] & edgeLabels2[bond2.getIdx()]).count()) {
    if (!ignoreAtomAromaticity) {
      if (!((bond1.getBeginAtom()->getIsAromatic() ==
                 bond2.getBeginAtom()->getIsAromatic() &&
             bond1.getEndAtom()->getIsAromatic() ==
                 bond2.getEndAtom()->getIsAromatic()) ||
            (bond1.getBeginAtom()->getIsAromatic() ==
                 bond2.getEndAtom()->getIsAromatic() &&
             bond1.getEndAtom()->getIsAromatic() ==
                 bond2.getBeginAtom()->getIsAromatic()))) {
        return false;
      }
    }
    if ((vtxLabels1[bond1.getBeginAtomIdx()] &
         vtxLabels2[bond2.getBeginAtomIdx()])
            .count() &&
        (vtxLabels1[bond1.getEndAtomIdx()] & vtxLabels2[bond2.getEndAtomIdx()])
            .count()) {
      return true;
    }
    if ((vtxLabels1[bond1.getBeginAtomIdx()] &
         vtxLabels2[bond2.getEndAtomIdx()])
            .count() &&
        (vtxLabels1[bond1.getEndAtomIdx()] &
         vtxLabels2[bond2.getBeginAtomIdx()])
            .count()) {
      return true;
    }
  }
  return false;
}

// Make the set of pairs of vertices, where they're a pair if the labels
// match.
void buildPairs(const ROMol &mol1,
                const std::vector<boost::dynamic_bitset<>> &vtxLabels1,
                const std::vector<boost::dynamic_bitset<>> &edgeLabels1,
                const ROMol &mol2,
                const std::vector<boost::dynamic_bitset<>> &vtxLabels2,
                const std::vector<boost::dynamic_bitset<>> &edgeLabels2,
                const RascalOptions &opts,
                std::vector<std::pair<int, int>> &vtxPairs) {
  std::vector<std::string> mol1RingSmiles, mol2RingSmiles;
  std::vector<std::unique_ptr<ROMol>> mol1Rings, mol2Rings;
  // For these purposes, it is correct that n1cccc1 and [nH]1cccc1 match - the
  // former would be from an N-substituted pyrrole, the latter from a plain
  // one.
  static const std::regex reg(R"(\[([np])H\])");
  if (opts.completeAromaticRings) {
    extractRings(mol1, mol1Rings, mol1RingSmiles);
    for (auto &mrs : mol1RingSmiles) {
      mrs = std::regex_replace(mrs, reg, "$1");
    }
    extractRings(mol2, mol2Rings, mol2RingSmiles);
    for (auto &mrs : mol2RingSmiles) {
      mrs = std::regex_replace(mrs, reg, "$1");
    }
  }

  for (const auto &bond1 : mol1.bonds()) {
    for (const auto &bond2 : mol2.bonds()) {
      if (bondsMatch(*bond1, *bond2, vtxLabels1, edgeLabels1, vtxLabels2,
                     edgeLabels2, opts.ignoreAtomAromaticity)) {
        if (opts.completeAromaticRings &&
            !checkAromaticRings(mol1, mol1RingSmiles, mol1Rings, vtxLabels1,
                                edgeLabels1, bond1->getIdx(), mol2,
                                mol2RingSmiles, mol2Rings, vtxLabels2,
                                edgeLabels2, bond2->getIdx())) {
          continue;
        }
        if (opts.ringMatchesRingOnly &&
            !checkRingMatchesRing(mol1, bond1->getIdx(), mol2,
                                  bond2->getIdx())) {
          continue;
        }
        vtxPairs.push_back(std::make_pair(bond1->getIdx(), bond2->getIdx()));
      }
    }
  }
}

// Make the modular product between the 2 graphs passed in.  Each node in the
// graph is a pair of vertices, one from the first graph, the other from the
// second, whose labels match.  Two vertices are connected in the modular
// product if either the 2 matching vertices in the 2 input vertices are
// connected by edges with the same label, or neither is connected.
void makeModularProduct(
    const ROMol &mol1,
    const std::vector<std::vector<const boost::dynamic_bitset<> *>> &adjMatrix1,
    const std::vector<boost::dynamic_bitset<>> &vtxLabels1,
    const std::vector<boost::dynamic_bitset<>> &edgeLabels1,
    const std::vector<std::vector<int>> &distMatrix1, const ROMol &mol2,
    const std::vector<std::vector<const boost::dynamic_bitset<> *>> &adjMatrix2,
    const std::vector<boost::dynamic_bitset<>> &vtxLabels2,
    const std::vector<boost::dynamic_bitset<>> &edgeLabels2,
    const std::vector<std::vector<int>> &distMatrix2, const RascalOptions &opts,
    std::vector<std::pair<int, int>> &vtxPairs,
    std::vector<boost::dynamic_bitset<>> &modProd) {
  buildPairs(mol1, vtxLabels1, edgeLabels1, mol2, vtxLabels2, edgeLabels2, opts,
             vtxPairs);
  if (vtxPairs.empty()) {
    // There was nothing in common at all.  But, what was the screening doing?
    modProd.clear();
    return;
  }
  if (vtxPairs.size() > opts.maxBondMatchPairs) {
    BOOST_LOG(rdErrorLog) << "Too many matching bond pairs (" << vtxPairs.size()
                          << ") so can't continue." << std::endl;
    modProd.clear();
    return;
  }
  modProd = std::vector<boost::dynamic_bitset<>>(
      vtxPairs.size(), boost::dynamic_bitset<>(vtxPairs.size()));
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
      if (opts.singleLargestFrag &&
          distMatrix1[vtxPairs[i].first][vtxPairs[j].first] !=
              distMatrix2[vtxPairs[i].second][vtxPairs[j].second]) {
        distsOk = false;
      }
      if (distsOk && (!adjMatrix1[vtxPairs[i].first][vtxPairs[j].first] &&
                      !adjMatrix2[vtxPairs[i].second][vtxPairs[j].second]) ||
          (adjMatrix1[vtxPairs[i].first][vtxPairs[j].first] &&
           adjMatrix2[vtxPairs[i].second][vtxPairs[j].second] &&
           (*adjMatrix1[vtxPairs[i].first][vtxPairs[j].first] &
            *adjMatrix2[vtxPairs[i].second][vtxPairs[j].second])
               .count())) {
        modProd[i][j] = modProd[j][i] = 1;
      }
    }
  }
}

// Calculate the lower bound on the size of the MCES.  This requires that mol1
// has more atoms than mol2 which is not checked.  Returns a minimum of 1.
unsigned int calcLowerBound(const ROMol &mol1, const ROMol &mol2,
                            double simThresh) {
  std::unordered_set<int> mol1AtNos;
  int maxAtNo = 0;
  for (const auto &a : mol1.atoms()) {
    mol1AtNos.insert(a->getAtomicNum());
    maxAtNo = std::max(a->getAtomicNum(), maxAtNo);
  }
  boost::dynamic_bitset<> mol2AtNos(maxAtNo + 1);
  for (const auto &a : mol2.atoms()) {
    // since we're interested in the atoms that match in the 2 molecules,
    // it doesn't matter if mol2 has an atomic number higher than anything
    // in mol1 - that can't be a match.
    if (a->getAtomicNum() < maxAtNo) {
      mol2AtNos.set(a->getAtomicNum());
    }
  }
  int deltaVg1 = 0;
  for (auto mol1AtNo : mol1AtNos) {
    if (!mol2AtNos[mol1AtNo]) {
      ++deltaVg1;
    }
  }
  double lb = sqrt((mol1.getNumAtoms() + mol1.getNumBonds()) *
                   (mol2.getNumAtoms() + mol2.getNumBonds()));
  lb = lb * simThresh - mol1.getNumAtoms() + deltaVg1;
  lb = lb < 0 ? 0 : lb;
  unsigned int ilb(lb);
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
  os << std::endl;
  os << "mol 1 bonds : [";
  for (auto mem : clique) {
    if (swapped) {
      os << vtxPairs[mem].second << ", ";
    } else {
      os << vtxPairs[mem].first << ", ";
    }
  }
  os << "]" << std::endl;
  os << "mol 2 bonds : [";
  for (auto mem : clique) {
    if (swapped) {
      os << vtxPairs[mem].first << ", ";
    } else {
      os << vtxPairs[mem].second << ", ";
    }
  }
  os << "]" << std::endl;
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
  auto *molFrags = new RWMol(mol);
  boost::dynamic_bitset<> aInClique(mol.getNumAtoms());
  boost::dynamic_bitset<> bInClique(mol.getNumBonds());
  for (auto mem : clique) {
    const Bond *bond = nullptr;
    if (pairNum == 1) {
      bond = molFrags->getBondWithIdx(vtxPairs[mem].first);
    } else {
      bond = molFrags->getBondWithIdx(vtxPairs[mem].second);
    }
    bInClique[bond->getIdx()] = 1;
    aInClique.set(bond->getBeginAtomIdx());
    bond->getBeginAtom()->setProp<int>("ORIG_INDEX", bond->getBeginAtomIdx());
    aInClique.set(bond->getEndAtomIdx());
    bond->getEndAtom()->setProp<int>("ORIG_INDEX", bond->getEndAtomIdx());
  }
  molFrags->beginBatchEdit();
  for (auto &a : molFrags->atoms()) {
    if (!aInClique[a->getIdx()]) {
      molFrags->removeAtom(a);
    }
  }
  for (auto &b : molFrags->bonds()) {
    if (!bInClique[b->getIdx()]) {
      molFrags->removeBond(b->getBeginAtomIdx(), b->getEndAtomIdx());
    }
  }
  molFrags->commitBatchEdit();
  return molFrags;
}

// Calculate the shortest bond distance between the 2 fragments in the
// molecule.
int minFragSeparation(const ROMol &mol, const ROMol &molFrags,
                      std::vector<int> &fragMapping, int frag1, int frag2) {
  auto extractFragAtoms = [&](int fragNum, std::vector<int> &fragAtoms) {
    for (size_t i = 0u; i < fragMapping.size(); ++i) {
      if (fragMapping[i] == fragNum) {
        int origIdx = molFrags.getAtomWithIdx(i)->getProp<int>("ORIG_INDEX");
        fragAtoms.push_back(origIdx);
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
      minDist = std::min(dist, minDist);
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
                     unsigned int &lowerBound) {
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
    if (!maxCliques.empty() && maxCliques.front().size() > lowerBound) {
      lowerBound = maxCliques.front().size();
    }
  }
}

// If the current time is beyond the timeout limit, throws a
// TimedOutException.
void checkTimeout(
    const std::chrono::time_point<std::chrono::high_resolution_clock>
        &startTime,
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

namespace {
bool hasSubstructMatch(const ROMol &mol, const ROMol &query) {
  SubstructMatchParameters ps;
  ps.maxMatches = 1;
  return !SubstructMatch(mol, query, ps).empty();
}
}  // namespace

// There are some simple substructures for which equivalent bond pruning isn't
// allowed.
bool checkEquivalentsAllowed(const ROMol &mol) {
  const static std::vector<std::string> notSmarts{
      "*~*", "*~*1~*~*~1", "*12~*~*~2~*~1", "*14~*(~*~2~3~4)~*~2~*~3~1"};
  static std::vector<std::unique_ptr<ROMol>> notStructs;
  if (notStructs.empty()) {
    for (const auto &smt : notSmarts) {
      notStructs.emplace_back(SmartsToMol(smt));
    }
  }
  const static std::vector<std::pair<unsigned int, unsigned int>> notStats{
      {2, 1}, {4, 4}, {4, 5}, {5, 8}};
  for (size_t i = 0; i < notStructs.size(); ++i) {
    if (mol.getNumAtoms() == notStats[i].first &&
        mol.getNumBonds() == notStats[i].second &&
        hasSubstructMatch(mol, *notStructs[i])) {
      return false;
    }
  }
  return true;
}

void explorePartitions(
    RascalStartPoint &starter,
    const std::chrono::time_point<std::chrono::high_resolution_clock>
        &startTime,
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
      if (clique.size() + part->numParts() < starter.d_lowerBound) {
        backtrack = true;
      }
    } else {
      if (clique.size() + part->numParts() <= starter.d_lowerBound) {
        backtrack = true;
      }
    }
    if (!backtrack) {
      if (opts.allBestMCESs) {
        goDeeper = clique.size() + part->upperBound() >= starter.d_lowerBound;
      } else {
        goDeeper = clique.size() + part->upperBound() > starter.d_lowerBound;
      }
      if (goDeeper) {
        if (!part->isEmpty()) {
          std::shared_ptr<PartitionSet> nextPart(new PartitionSet(*part));
          clique.push_back(nextPart->popLastVertex());
          if (clique.size() == 1 && canDoEquivs &&
              equivalentRootAlreadyDone(clique.front(), starter.d_vtxPairs,
                                        starter.d_equivBonds1,
                                        starter.d_equivBonds2, rootClasses)) {
            clique.pop_back();
            backtrack = true;
          } else {
            nextPart->pruneVertices(clique.back());
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
    if (backtrack || (!parts.empty() && parts.back()->isEmpty())) {
      while (!parts.empty()) {
        if (parts.back()->isEmpty()) {
          parts.pop_back();
          if (!clique.empty()) {
            clique.pop_back();
          }
        } else {
          parts.back()->popLastVertex();
          if (!parts.back()->isEmpty()) {
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
  // isomorphic line graphs.  This checks to see if that's something we need
  // to worry about for these molecules.
  const static std::unique_ptr<ROMol> delta(SmartsToMol("C1CC1"));
  const static std::unique_ptr<ROMol> y(SmartsToMol("C(C)C"));
  return (hasSubstructMatch(mol1, *delta) && hasSubstructMatch(mol2, *y)) ||
         (hasSubstructMatch(mol2, *delta) && hasSubstructMatch(mol1, *y));
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
// Adapted from https://en.wikipedia.org/wiki/Floyd–Warshall_algorithm
void calcDistMatrix(
    const std::vector<std::vector<const boost::dynamic_bitset<> *>> &adjMatrix,
    std::vector<std::vector<int>> &distMatrix) {
  distMatrix = std::vector<std::vector<int>>(
      adjMatrix.size(),
      std::vector<int>(adjMatrix.size(), adjMatrix.size() + 1));
  for (size_t i = 0u; i < adjMatrix.size(); ++i) {
    distMatrix[i][i] = 0;
    for (size_t j = 0u; j < adjMatrix.size(); ++j) {
      if (i != j && adjMatrix[i][j]) {
        distMatrix[i][j] = 1;
      }
    }
  }
  for (size_t k = 0u; k < adjMatrix.size(); ++k) {
    for (size_t i = 0u; i < adjMatrix.size(); ++i) {
      for (size_t j = 0u; j < adjMatrix.size(); ++j) {
        if (distMatrix[i][j] > distMatrix[i][k] + distMatrix[k][j]) {
          distMatrix[i][j] = distMatrix[i][k] + distMatrix[k][j];
        }
      }
    }
  }
}

namespace details {
// Set a bit in the bitset for each atom for its atomic number OR
// all the equivalent atom classes it is in.  If an O atom in one
// molecule is in a class "H bond donor", for example, and the only
// O atom in the other molecule is not, we don't want them to match.
// Add any matching SMARTS to the atom as the property "EQUIV_SMARTS".
void makeAtomLabels(const ROMol &mol, const std::string &equivalentAtoms,
                    std::vector<boost::dynamic_bitset<>> &atomLabels) {
  std::vector<std::string> classSmarts;
  if (!equivalentAtoms.empty()) {
    boost::split(classSmarts, equivalentAtoms, boost::is_any_of(" "),
                 boost::token_compress_on);
    classSmarts.erase(
        std::remove_if(classSmarts.begin(), classSmarts.end(),
                       [](const auto &s) -> bool { return s.empty(); }),
        classSmarts.end());
  }
  const unsigned int numBits = elementNames.size() + 1 + classSmarts.size();
  atomLabels.resize(mol.getNumAtoms());
  for (const auto &atom : mol.atoms()) {
    atomLabels[atom->getIdx()] = boost::dynamic_bitset<>(numBits);
    atomLabels[atom->getIdx()][atom->getAtomicNum()] = true;
  }

  for (size_t i = 0; i < classSmarts.size(); ++i) {
    if (classSmarts[i].empty()) {
      continue;
    }
    auto qmol = v2::SmilesParse::MolFromSmarts(classSmarts[i]);
    std::vector<RDKit::MatchVectType> hits_vect;
    if (RDKit::SubstructMatch(mol, *qmol, hits_vect)) {
      for (const auto &hv : hits_vect) {
        for (const auto &h : hv) {
          auto a = mol.getAtomWithIdx(h.second);
          atomLabels[h.second][i + elementNames.size()] = true;
          atomLabels[h.second][a->getAtomicNum()] = false;
          if (a->hasProp("EQUIV_SMARTS")) {
            auto currProp = a->getProp<std::string>("EQUIV_SMARTS");
            currProp += " " + classSmarts[i];
            a->setProp("EQUIV_SMARTS", currProp);
          } else {
            a->setProp<std::string>("EQUIV_SMARTS", classSmarts[i]);
          }
        }
      }
    }
  }
}

void makeBondBitstrings(const ROMol &mol, const bool ignoreBondOrders,
                        std::vector<boost::dynamic_bitset<>> &bondLabels) {
  bondLabels.resize(mol.getNumBonds());
  for (const auto &bond : mol.bonds()) {
    bondLabels[bond->getIdx()] =
        boost::dynamic_bitset<>(Bond::BondType::ZERO + 1);
    if (ignoreBondOrders) {
      bondLabels[bond->getIdx()][0] = true;
    } else {
      bondLabels[bond->getIdx()][bond->getBondType()] = true;
    }
  }
}
}  // namespace details

// Add the unique entries in bitsets to uniqueBitSets.
void uniqueBitsets(const std::vector<boost::dynamic_bitset<>> &bitsets,
                   std::vector<boost::dynamic_bitset<>> &uniqueBitsets) {
  uniqueBitsets.insert(uniqueBitsets.end(), bitsets.begin(), bitsets.end());
  std::sort(uniqueBitsets.begin(), uniqueBitsets.end(), bitsetLess);
  uniqueBitsets.erase(std::unique(uniqueBitsets.begin(), uniqueBitsets.end()),
                      uniqueBitsets.end());
}

// Make the bond labels as indices into the unique bond and atom labels in
// order bond, atom1, atom2 where atom1 and atom2 are in their bitset
// order, i.e. uniqueAtomLabels[atom1] < uniqueAtomLabels[atom2].
void makeBondLabels(
    const ROMol &mol, const std::vector<boost::dynamic_bitset<>> &bondLabels,
    const std::vector<boost::dynamic_bitset<>> &uniqueBondLabels,
    const std::vector<boost::dynamic_bitset<>> &atomLabels,
    const std::vector<boost::dynamic_bitset<>> &uniqueAtomLabels,
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        &bondTypeLabels) {
  auto getType =
      [](const boost::dynamic_bitset<> &bs,
         const std::vector<boost::dynamic_bitset<>> &types) -> unsigned int {
    const auto it = std::find(types.begin(), types.end(), bs);
    return std::distance(types.begin(), it);
  };

  bondTypeLabels.resize(mol.getNumBonds());
  for (const auto bond : mol.bonds()) {
    unsigned int atom1 = bond->getBeginAtomIdx();
    unsigned int atom2 = bond->getEndAtomIdx();
    if (atomLabels[atom2] < atomLabels[atom1]) {
      std::swap(atom1, atom2);
    }
    bondTypeLabels[bond->getIdx()] =
        std::make_tuple(getType(bondLabels[bond->getIdx()], uniqueBondLabels),
                        getType(atomLabels[atom1], uniqueAtomLabels),
                        getType(atomLabels[atom2], uniqueAtomLabels));
  }
}

RascalStartPoint makeInitialPartitionSet(const ROMol *mol1, const ROMol *mol2,
                                         const RascalOptions &opts) {
  RascalStartPoint starter;
  if (mol1->getNumAtoms() <= mol2->getNumAtoms()) {
    starter.d_swapped = false;
    starter.d_mol1.reset(new ROMol(*mol1));
    starter.d_mol2.reset(new ROMol(*mol2));
  } else {
    starter.d_swapped = true;
    starter.d_mol1.reset(new ROMol(*mol2));
    starter.d_mol2.reset(new ROMol(*mol1));
  }
  std::vector<boost::dynamic_bitset<>> atomLabels1, atomLabels2;
  details::makeAtomLabels(*starter.d_mol1, opts.equivalentAtoms, atomLabels1);
  details::makeAtomLabels(*starter.d_mol2, opts.equivalentAtoms, atomLabels2);

  std::map<int, std::vector<std::pair<int, int>>> degSeqs1, degSeqs2;
  starter.d_tier1Sim =
      details::tier1Sim(*starter.d_mol1, *starter.d_mol2, atomLabels1,
                        atomLabels2, degSeqs1, degSeqs2);
  if (starter.d_tier1Sim < opts.similarityThreshold) {
    return starter;
  }

  std::vector<boost::dynamic_bitset<>> bondStrings1, bondStrings2;
  details::makeBondBitstrings(*starter.d_mol1, opts.ignoreBondOrders,
                              bondStrings1);
  details::makeBondBitstrings(*starter.d_mol2, opts.ignoreBondOrders,
                              bondStrings2);

  starter.d_tier2Sim =
      details::tier2Sim(*starter.d_mol1, *starter.d_mol2, degSeqs1, degSeqs2,
                        atomLabels1, atomLabels2, bondStrings1, bondStrings2);
  if (starter.d_tier2Sim < opts.similarityThreshold) {
    return starter;
  }

  // Get the line graphs for the two molecules as adjacency matrices.
  makeLineGraph(*starter.d_mol1, atomLabels1, starter.d_adjMatrix1);
  makeLineGraph(*starter.d_mol2, atomLabels2, starter.d_adjMatrix2);

  std::vector<std::vector<int>> distMat1, distMat2;
  if (opts.maxFragSeparation > -1 || opts.singleLargestFrag) {
    calcDistMatrix(starter.d_adjMatrix1, distMat1);
    calcDistMatrix(starter.d_adjMatrix2, distMat2);
  }

  // pairs are vertices in the 2 line graphs that are the same type.
  // d_modProd is the modular product/correspondence graph of the two
  // line graphs.
  // makeModularProduct(*starter.d_mol1, starter.d_adjMatrix1, bondLabels1,
  // distMat1, *starter.d_mol2, starter.d_adjMatrix2,
  // bondLabels2, distMat2, opts, starter.d_vtxPairs,
  // starter.d_modProd);
  makeModularProduct(*starter.d_mol1, starter.d_adjMatrix1, atomLabels1,
                     bondStrings1, distMat1, *starter.d_mol2,
                     starter.d_adjMatrix2, atomLabels2, bondStrings2, distMat2,
                     opts, starter.d_vtxPairs, starter.d_modProd);
  if (starter.d_modProd.empty()) {
    return starter;
  }
  if (opts.minCliqueSize > 0) {
    starter.d_lowerBound = opts.minCliqueSize;
  } else {
    starter.d_lowerBound = calcLowerBound(*starter.d_mol1, *starter.d_mol2,
                                          opts.similarityThreshold);
  }
  std::vector<boost::dynamic_bitset<>> uniqueAtomLabels, uniqueBondLabels;
  uniqueBitsets(bondStrings1, uniqueBondLabels);
  uniqueBitsets(bondStrings2, uniqueBondLabels);
  uniqueBitsets(atomLabels1, uniqueAtomLabels);
  uniqueBitsets(atomLabels2, uniqueAtomLabels);
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> bondLabels1,
      bondLabels2;
  makeBondLabels(*starter.d_mol1, bondStrings1, uniqueBondLabels, atomLabels1,
                 uniqueAtomLabels, bondLabels1);
  makeBondLabels(*starter.d_mol2, bondStrings2, uniqueBondLabels, atomLabels2,
                 uniqueAtomLabels, bondLabels2);

  starter.d_partSet.reset(new PartitionSet(
      starter.d_modProd, starter.d_vtxPairs, bondLabels1, bondLabels2,
      uniqueAtomLabels, uniqueBondLabels, starter.d_lowerBound));

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

std::vector<RascalResult> findMCES(RascalStartPoint &starter,
                                   const RascalOptions &opts) {
  std::vector<unsigned int> clique;
  std::vector<std::vector<unsigned int>> maxCliques;
  auto startTime = std::chrono::high_resolution_clock::now();
  bool timedOut = false;
  RascalOptions tmpOpts{opts};
  if (opts.singleLargestFrag) {
    tmpOpts.allBestMCESs = true;
  }
  try {
    explorePartitions(starter, startTime, tmpOpts, maxCliques);
  } catch (TimedOutException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
    maxCliques = e.d_cliques;
    timedOut = true;
  }
  std::vector<RascalResult> results;
  for (const auto &c : maxCliques) {
    results.push_back(
        RascalResult(*starter.d_mol1, *starter.d_mol2, starter.d_adjMatrix1,
                     starter.d_adjMatrix2, c, starter.d_vtxPairs, timedOut,
                     starter.d_swapped, starter.d_tier1Sim, starter.d_tier2Sim,
                     opts.ringMatchesRingOnly, opts.singleLargestFrag,
                     opts.maxFragSeparation, opts.exactConnectionsMatch,
                     opts.equivalentAtoms, opts.ignoreBondOrders));
  }
  if (opts.singleLargestFrag) {
    std::sort(
        results.begin(), results.end(),
        [](const RascalResult &r1, const RascalResult &r2) -> bool {
          if (r1.getAtomMatches().size() == r2.getAtomMatches().size()) {
            if (r1.getBondMatches().size() == r2.getBondMatches().size()) {
              if (r1.getAtomMatches() == r2.getAtomMatches()) {
                return (r1.getBondMatches() < r2.getBondMatches());
              }
              return r1.getAtomMatches() < r2.getAtomMatches();
            }
            return r1.getBondMatches().size() > r2.getBondMatches().size();
          }
          return (r1.getAtomMatches().size() > r2.getAtomMatches().size());
        });

    // the singleLargestFrag method throws bits of solutions out, so there may
    // now be duplicates and results that are different sizes.
    results.erase(
        std::unique(results.begin(), results.end(),
                    [](const RascalResult &r1, const RascalResult &r2) -> bool {
                      return (r1.getAtomMatches() == r2.getAtomMatches() &&
                              r1.getBondMatches() == r2.getBondMatches());
                    }),
        results.end());
    boost::dynamic_bitset<> want(results.size());
    want.set();
    for (size_t i = 1; i < results.size(); ++i) {
      if (results[i].getAtomMatches().size() <
              results[0].getAtomMatches().size() ||
          results[i].getBondMatches().size() <
              results[0].getBondMatches().size()) {
        want[i] = false;
      }
    }
    if (want.count() < results.size()) {
      size_t j = 0;
      for (size_t i = 0; i < results.size(); ++i) {
        if (want[i]) {
          results[j++] = results[i];
        }
      }
      results.erase(results.begin() + j, results.end());
    }
  } else {
    // If 2 cliques are the same size, this sort puts the one with the smaller
    // number of fragments first, which may have fewer atoms.
    std::sort(results.begin(), results.end(), details::resultCompare);
  }
  return results;
}

// calculate the RASCAL MCES between the 2 molecules, provided it is within
// the similarity threshold given.
std::vector<RascalResult> rascalMCES(const ROMol &mol1, const ROMol &mol2,
                                     const RascalOptions &opts) {
  auto starter = makeInitialPartitionSet(&mol1, &mol2, opts);
  if (!starter.d_partSet) {
    if (opts.returnEmptyMCES) {
      return std::vector<RascalResult>(
          1, RascalResult(starter.d_tier1Sim, starter.d_tier2Sim));
    }
    return std::vector<RascalResult>();
  }
  auto results = findMCES(starter, opts);
  if (results.empty() && opts.returnEmptyMCES) {
    return std::vector<RascalResult>(
        1, RascalResult(starter.d_tier1Sim, starter.d_tier2Sim));
  }
  if (!opts.allBestMCESs && results.size() > 1) {
    results.erase(results.begin() + 1, results.end());
  }
  return results;
}

}  // namespace RascalMCES
}  // namespace RDKit
