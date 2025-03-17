//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <memory>

#include "PartitionSet.h"

namespace RDKit {

namespace RascalMCES {
// This is in lap_a_la_scipy.cpp and solves the linear assignment problem.
int lapMaximize(const std::vector<std::vector<int>> &costsMat,
                std::vector<size_t> &a, std::vector<size_t> &b);

PartitionSet::PartitionSet(
    const std::vector<boost::dynamic_bitset<>> &modProd,
    const std::vector<std::pair<int, int>> &vtxPairs,
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        &vtx1Labels,
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        &vtx2Labels,
    const std::vector<boost::dynamic_bitset<>> &uniqueAtomLabels,
    const std::vector<boost::dynamic_bitset<>> &uniqueBondLabels,
    unsigned int lowerBound)
    : d_ModProd(new std::vector<boost::dynamic_bitset<>>(modProd)),
      d_VtxPairs(new std::vector<std::pair<int, int>>(vtxPairs)),
      d_vtx1Labels(
          new std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>(
              vtx1Labels)),
      d_vtx2Labels(
          new std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>(
              vtx2Labels)),
      d_uniqueAtomLabels(
          new std::vector<boost::dynamic_bitset<>>(uniqueAtomLabels)),
      d_uniqueBondLabels(
          new std::vector<boost::dynamic_bitset<>>(uniqueBondLabels)) {
  d_vtx1Counts = std::vector<int>(d_vtx1Labels->size(), 0);
  d_vtx2Counts = std::vector<int>(d_vtx2Labels->size(), 0);
  int firstVtx = -1;
  // Clearly, a vertex in one of the line graphs can only match one vertex
  // in the other.  Thus, the initial partitions can be set up so that
  // all vertices in a partition have the same vertex in the first
  // line graph.
  for (size_t i = 0; i < vtxPairs.size(); ++i) {
    auto &vp = vtxPairs[i];
    if (vp.first != firstVtx) {
      d_parts.push_back(std::vector<unsigned int>());
      d_parts.back().push_back(i);
      firstVtx = vp.first;
    } else {
      d_parts.back().push_back(i);
    }
    d_vtx1Counts[vp.first]++;
    d_vtx2Counts[vp.second]++;
  }
  if (d_parts.empty()) {
    return;
  }
  // Now sort the partitions by size.  This means that the vertices at the
  // top of the partition set, above the lowerBound (or Pex as Raymond
  // calls it in the paper), are the ones that match the least number of
  // vertices in the other line graph.  This has a dramatic effect on the
  // speed compared with other things tried.  I think it is what Raymond
  // means when he says "Perform an initial partitioning of the vertices...
  // using the labeled edge projection procedure."
  sortPartitions();
  // Now reassign vertices from above Pex to below it if possible.
  // This also improves the speed of finding a large clique early.
  // A vertex is moved to a partition where it isn't connected to a vertex
  // in the modular product graph that is in the partition.
  for (size_t i = d_parts.size() - 1; i > lowerBound; --i) {
    bool reassigned = false;
    for (auto &iv : d_parts[i]) {
      for (size_t k = 0; k <= lowerBound; ++k) {
        bool conn = false;
        for (auto kv : d_parts[k]) {
          if (modProd[iv][kv]) {
            conn = true;
            break;
          }
        }
        if (!conn) {
          d_parts[k].push_back(iv);
          iv = std::numeric_limits<unsigned int>::max();
          reassigned = true;
          break;
        }
      }
    }
    if (reassigned) {
      d_parts[i].erase(std::remove(d_parts[i].begin(), d_parts[i].end(),
                                   std::numeric_limits<unsigned int>::max()),
                       d_parts[i].end());
    }
  }
  d_parts.erase(std::remove_if(d_parts.begin(), d_parts.end(),
                               [](const std::vector<unsigned int> &v) {
                                 return v.empty();
                               }),
                d_parts.end());
  // Sort again, to make sure the large partitions are dealt with as late as
  // possible.
  sortPartitions();

  // Get the info together for the upper bound calculation.
  calcVtxTypeCounts();
}

int PartitionSet::upperBound() {
  int upperBound = 0;
  if (d_vtxInMoreThanOneClass) {
    // It's complicated because we don't want to over-count the matching
    // vertices.  Use the LAP algorithm to match a single vertex in mol1
    // with a single vertex in mol2.  It'll be a lot slower.
    constexpr size_t unassignedValue(99999999);
    std::vector<std::vector<int>> costsMat(
        d_vtx1Labels->size(), std::vector<int>(d_vtx2Labels->size(), 0));
    for (size_t i = 0; i < d_vtx1TypeCountVtxes.size(); ++i) {
      if (!d_vtx1TypeCounts[i]) {
        continue;
      }
      for (size_t j = 0; j < d_vtx1TypeCountVtxes[i].size(); ++j) {
        if (d_vtx1TypeCountVtxes[i][j]) {
          for (size_t k = 0; k < d_vtx2TypeCountVtxes[i].size(); ++k) {
            if (d_vtx2TypeCountVtxes[i][k]) {
              costsMat[j][k] = 1;
            }
          }
        }
      }
    }
    std::vector<size_t> a(std::min(d_vtx1Labels->size(), d_vtx2Labels->size()),
                          unassignedValue);
    std::vector<size_t> b(std::min(d_vtx1Labels->size(), d_vtx2Labels->size()),
                          unassignedValue);
    int retVal = lapMaximize(costsMat, a, b);
    if (retVal >= 0) {
      for (auto i = 0u; i < a.size(); ++i) {
        upperBound += costsMat[a[i]][b[i]];
      }
      return upperBound;
    }
  }

  // If we're here, it's either a straightforward calculation or the LAP failed.
  // In the latter case, we'll be getting an over-estimate which is ok but means
  // unnecessary levels in the search tree will be attempted.
  for (size_t i = 0; i < d_vtx1TypeCounts.size(); ++i) {
    if (!d_vtx1TypeCounts[i]) {
      continue;
    }
    upperBound += std::min(d_vtx1TypeCounts[i], d_vtx2TypeCounts[i]);
  }
  return upperBound;
}

unsigned int PartitionSet::popLastVertex() {
  if (d_parts.empty()) {
    throw std::runtime_error("PartitionSet set is empty.");
  }
  unsigned int ret_val = d_parts.back().back();
  d_parts.back().pop_back();
  if (d_parts.back().empty()) {
    d_parts.pop_back();
  }
  decrementVertexCounts(ret_val);
  return ret_val;
}

void PartitionSet::pruneVertices(unsigned int vtx_num) {
  for (auto &part : d_parts) {
    size_t i = 0;
    while (i < part.size()) {
      if (!(*d_ModProd)[part[i]][vtx_num]) {
        decrementVertexCounts(part[i]);
        part[i] = part.back();
        part.pop_back();
      } else {
        ++i;
      }
    }
  }
  d_parts.erase(std::remove_if(d_parts.begin(), d_parts.end(),
                               [](const std::vector<unsigned int> &v) {
                                 return v.empty();
                               }),
                d_parts.end());
  sortPartitions();
}

void PartitionSet::sortPartitions() {
  // When sorting lists with duplicate values, the order of the
  // duplicates isn't defined.  Different compilers do it differently.
  // This can affect the results in the case where more than 1 MCES is
  // possible, because the partition orders and hence the search tree
  // traversal will be different.  The results should be equivalent,
  // though.  To make things consistent, the sort is done with a
  // tie-breaker on the first value in vectors of the same size.  It
  // doesn't slow things down very much on average, and it makes things
  // tidier.
  std::sort(d_parts.begin(), d_parts.end(),
            [](const std::vector<unsigned int> &v1,
               const std::vector<unsigned int> &v2) {
              if (v1.size() == v2.size() && !v1.empty()) {
                return v1.front() < v2.front();
              } else {
                return v1.size() > v2.size();
              }
            });
}

namespace {
void makeIndVertexTypes(
    const std::vector<boost::dynamic_bitset<>> &uniqueLabels,
    std::vector<boost::dynamic_bitset<>> &indVtxTypes) {
  boost::dynamic_bitset<> allAtomTypes(uniqueLabels.front().size());
  for (const auto &ual : uniqueLabels) {
    allAtomTypes |= ual;
  }
  indVtxTypes.resize(allAtomTypes.count());
  for (size_t i = 0; i < allAtomTypes.count(); ++i) {
    indVtxTypes[i] = boost::dynamic_bitset<>(uniqueLabels.front().size());
  }
  size_t j = 0;
  for (size_t k = 0; k < allAtomTypes.size(); ++k) {
    if (allAtomTypes[k]) {
      indVtxTypes[j++][k] = true;
    }
  }
}
}  // namespace

void PartitionSet::calcVtxTypeCounts() {
  makeIndVertexTypes(*d_uniqueAtomLabels, d_indAtomTypes);
  makeIndVertexTypes(*d_uniqueBondLabels, d_indBondTypes);

  size_t numVtxTypes =
      d_indBondTypes.size() * d_indAtomTypes.size() * d_indAtomTypes.size();
  calcVtxTypeCounts(numVtxTypes, d_vtx1Counts, *d_vtx1Labels, d_vtx1TypeCounts,
                    d_vtx1TypeCountVtxes);
  calcVtxTypeCounts(numVtxTypes, d_vtx2Counts, *d_vtx2Labels, d_vtx2TypeCounts,
                    d_vtx2TypeCountVtxes);
}

void PartitionSet::calcVtxTypeCounts(
    size_t numVtxTypes, const std::vector<int> &vtxCounts,
    const std::vector<std::tuple<unsigned, unsigned, unsigned>> &vtxLabels,
    std::vector<int> &vtxTypeCounts,
    std::vector<boost::dynamic_bitset<>> &vtxTypeCountVtxes) {
  vtxTypeCounts.resize(numVtxTypes);
  vtxTypeCountVtxes = std::vector<boost::dynamic_bitset<>>(
      numVtxTypes, boost::dynamic_bitset<>(vtxCounts.size()));
  boost::dynamic_bitset<> seenVtx(vtxCounts.size());
  for (size_t i = 0; i < vtxCounts.size(); ++i) {
    const auto &bondLabel = (*d_uniqueBondLabels)[std::get<0>(vtxLabels[i])];
    const auto &atom1Label = (*d_uniqueAtomLabels)[std::get<1>(vtxLabels[i])];
    const auto &atom2Label = (*d_uniqueAtomLabels)[std::get<2>(vtxLabels[i])];
    if (vtxCounts[i]) {
      size_t num = 0;
      for (size_t bond = 0; bond < d_indBondTypes.size(); ++bond) {
        for (size_t atom1 = 0; atom1 < d_indAtomTypes.size(); ++atom1) {
          for (size_t atom2 = 0; atom2 < d_indAtomTypes.size();
               ++atom2, ++num) {
            if ((bondLabel & d_indBondTypes[bond]).count() &&
                (atom1Label & d_indAtomTypes[atom1]).count() &&
                (atom2Label & d_indAtomTypes[atom2]).count()) {
              vtxTypeCounts[num]++;
              vtxTypeCountVtxes[num][i] = true;
              if (seenVtx[i]) {
                d_vtxInMoreThanOneClass = true;
              } else {
                seenVtx[i] = true;
              }
            }
          }
        }
      }
    }
  }
}

namespace {
void reassessVertexInMoreThanOneClass(
    const std::vector<boost::dynamic_bitset<>> &vtxTypeCountVtxes,
    unsigned int numVtxes, bool &vtxInMoreThanOneClass) {
  boost::dynamic_bitset<> seenVtx(numVtxes);
  for (const auto &vtxTypeCount : vtxTypeCountVtxes) {
    for (size_t i = 0; i < vtxTypeCount.size(); ++i) {
      if (vtxTypeCount[i]) {
        if (seenVtx[i]) {
          vtxInMoreThanOneClass = true;
          return;
        } else {
          seenVtx[i] = true;
        }
      }
    }
  }
}
}  // namespace

void PartitionSet::decrementVertexCounts(int vtxNum) {
  --d_vtx1Counts[(*d_VtxPairs)[vtxNum].first];
  std::unique_ptr<boost::dynamic_bitset<>> seenVtx1s, seenVtx2s;
  bool reassess = false;
  // If this vertex count has gone to zero, see if it has emptied
  // a type.  If it has, the upperBound will have changed.
  if (!d_vtx1Counts[(*d_VtxPairs)[vtxNum].first]) {
    for (size_t i = 0; i < d_vtx1TypeCountVtxes.size(); i++) {
      if (d_vtx1TypeCountVtxes[i].size()) {
        d_vtx1TypeCountVtxes[i][(*d_VtxPairs)[vtxNum].first] = false;
        d_vtx1TypeCounts[i] = static_cast<int>(d_vtx1TypeCountVtxes[i].count());
        reassess = true;
      }
    }
  }
  --d_vtx2Counts[(*d_VtxPairs)[vtxNum].second];
  if (!d_vtx2Counts[(*d_VtxPairs)[vtxNum].second]) {
    for (size_t i = 0; i < d_vtx2TypeCountVtxes.size(); i++) {
      if (d_vtx2TypeCountVtxes[i].size()) {
        d_vtx2TypeCountVtxes[i][(*d_VtxPairs)[vtxNum].second] = false;
        d_vtx2TypeCounts[i] = static_cast<int>(d_vtx2TypeCountVtxes[i].count());
        reassess = true;
      }
    }
  }
  // Re-assess d_vtxInMoreThanOneClass, because if we can avoid the LAP
  // in upperBound it will save a lot of time.
  if (reassess && d_vtxInMoreThanOneClass) {
    d_vtxInMoreThanOneClass = false;
    reassessVertexInMoreThanOneClass(d_vtx1TypeCountVtxes, d_vtx1Counts.size(),
                                     d_vtxInMoreThanOneClass);
    if (!d_vtxInMoreThanOneClass) {
      reassessVertexInMoreThanOneClass(
          d_vtx2TypeCountVtxes, d_vtx2Counts.size(), d_vtxInMoreThanOneClass);
    }
  }
}

std::ostream &operator<<(std::ostream &os, const PartitionSet &pt) {
  for (size_t i = 0; i < pt.d_parts.size(); ++i) {
    os << i << " :: " << pt.d_parts[i].size() << " ::";
    for (auto &mem : pt.d_parts[i]) {
      os << " " << mem << " (" << (*pt.d_VtxPairs)[mem].first << ","
         << (*pt.d_VtxPairs)[mem].second << ")";
    }
    os << std::endl;
  }
  os << "vtx1_counts :";
  for (auto vc : pt.d_vtx1Counts) {
    os << " " << vc;
  }
  os << std::endl;
  os << "vtx2_counts :";
  for (auto vc : pt.d_vtx2Counts) {
    os << " " << vc;
  }
  os << std::endl;
  return os;
}

}  // namespace RascalMCES
}  // namespace RDKit
