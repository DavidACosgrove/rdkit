//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RASCALMCES_PARTITION_SET_H
#define RASCALMCES_PARTITION_SET_H

#include <map>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace RascalMCES {

class PartitionSet {
 public:
  // Make a partition set from the modular product and the labels
  // of the vertices from the first graph.  Each element in vtxPairs
  // has a row/column in modProd.  The partitions are sorted
  // into descending order of sizes.
  PartitionSet(
      const std::vector<boost::dynamic_bitset<>> &modProd,
      const std::vector<std::pair<int, int>> &vtxPairs,
      const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
          &vtx1Labels,
      const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
          &vtx2Labels,
      const std::vector<boost::dynamic_bitset<>> &uniqueAtomLabels,
      const std::vector<boost::dynamic_bitset<>> &uniqueBondLabels,
      unsigned int lowerBound);

  bool isEmpty() const { return d_parts.empty(); }

  size_t numParts() const { return d_parts.size(); }

  // Compute the upper bound on the clique that can be extracted from
  // the current partition.
  int upperBound();

  friend std::ostream &operator<<(std::ostream &os, const PartitionSet &pt);

  // removes the last element of the last partition and returns
  // its value. Throws a runtime_error if empty.
  unsigned int popLastVertex();

  // remove from the partitions any vertex not connected to the given
  // vertex
  void pruneVertices(unsigned int vtx_num);

 private:
  std::shared_ptr<const std::vector<boost::dynamic_bitset<>>> d_ModProd;
  std::shared_ptr<const std::vector<std::pair<int, int>>> d_VtxPairs;
  std::shared_ptr<
      const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>>
      d_vtx1Labels;
  std::shared_ptr<
      const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>>
      d_vtx2Labels;
  std::shared_ptr<const std::vector<boost::dynamic_bitset<>>>
      d_uniqueAtomLabels;
  std::shared_ptr<const std::vector<boost::dynamic_bitset<>>>
      d_uniqueBondLabels;
  std::vector<std::vector<unsigned int>> d_parts;
  // counts of the number of times each vertex appears in the partitions
  std::vector<int> d_vtx1Counts, d_vtx2Counts;
  // counts of the number of times the d_vtx[12]_labels appear in the partitions
  std::vector<int> d_vtx1TypeCounts, d_vtx2TypeCounts;
  // These are the vertices of the type, so
  // d_vtx1TypeCountVtxes[i].count() == d_vtx1TypeCounts[i].
  std::vector<boost::dynamic_bitset<>> d_vtx1TypeCountVtxes,
      d_vtx2TypeCountVtxes;
  bool d_vtxInMoreThanOneClass{false};
  // Bitsets for the individual atom and bond types present in the vertices
  std::vector<boost::dynamic_bitset<>> d_indAtomTypes;
  std::vector<boost::dynamic_bitset<>> d_indBondTypes;

  void sortPartitions();

  void calcVtxTypeCounts();
  void calcVtxTypeCounts(
      size_t numVtxTypes, const std::vector<int> &vtxCounts,
      const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
          &vtxLabels,
      std::vector<int> &vtxTypeCounts,
      std::vector<boost::dynamic_bitset<>> &vtxTypeCountVtxes);

  void decrementVertexCounts(int vtxNum);
};
}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALMCES_PARTITION_SET_H
