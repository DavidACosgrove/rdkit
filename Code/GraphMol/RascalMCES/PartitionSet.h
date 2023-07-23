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

namespace RDKit {

namespace RascalMCES {

class PartitionSet {
 public:
  // Make a partition set from the modular product and the labels
  // of the vertices from the first graph.  Each element in vtxPairs
  // has a row/column in modProd.  The partitions are sorted
  // into descending order of sizes. Because that uses std::sort,
  // the final order may differ from compiler to compiler, because
  // the standard doesn't define what should happen in the case of
  // duplicate values.  It is quite common to have partitions of the
  // same size.
  PartitionSet(const std::vector<std::vector<char>> &modProd,
               const std::vector<std::pair<int, int>> &vtxPairs,
               const std::vector<unsigned int> &vtx1Labels,
               const std::vector<unsigned int> &vtx2Labels,
               unsigned int lowerBound);

  bool isEmpty() const { return d_parts.empty(); }

  size_t numParts() const { return d_parts.size(); }

  // Compute the upper bound on the clique that can be extracted from
  // the current partition.
  int upperBound();

  // Print the partitions to os.  Very useful for debugging, but not
  // used otherwise.
  void printPartitions(std::ostream &os) const;

  // removes the last element of the last partition and returns
  // its value. Throws a runtime_error if empty.
  unsigned int popLastVertex();

  // remove from the partitions any vertex not connected to the given
  // vertex
  void pruneVertices(unsigned int vtx_num);

 private:
  std::shared_ptr<const std::vector<std::vector<char>>> d_ModProd;
  std::shared_ptr<const std::vector<std::pair<int, int>>> d_VtxPairs;
  std::shared_ptr<const std::vector<unsigned int>> d_vtx1Labels;
  std::shared_ptr<const std::vector<unsigned int>> d_vtx2Labels;
  std::vector<std::vector<unsigned int>> d_parts;
  // counts of the number of times each vertex appears in the partitions
  std::vector<int> d_vtx1Counts, d_vtx2Counts;
  // counts of the number of times the d_vtx[12]_labels appear in the partitions
  std::vector<int> d_vtx1TypeCounts, d_vtx2TypeCounts;

  void addVertex(unsigned int vtxNum);

  void sortPartitions();

  void calcVtxTypeCounts();

  void decrementVertexCounts(int vtxNum);
};
}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALMCES_PARTITION_SET_H
