//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This file declares a concrete class derived from SynthonSpaceSearcher
// that does shape similarity searching of the SynthonSpace using the
// PubChem code.

#ifndef SYNTHONSPACESHAPESEARCHER_H
#define SYNTHONSPACESHAPESEARCHER_H

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>

namespace RDKit::SynthonSpaceSearch {

// Concrete class that does the search by fingerprint similarity.
class SynthonSpaceShapeSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceShapeSearcher() = delete;
  SynthonSpaceShapeSearcher(const ROMol &query,
                            const SynthonSpaceSearchParams &params,
                            SynthonSpace &space);

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      const std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const override;

 private:
  ShapeSet d_queryShapes;
  // These are the fragment shapes for this search, derived from
  // d_query.  The shapes in d_fragShapes are keyed on the address
  // of the corresponding fragment.  d_fragShapesPool is never read,
  // it is just used a repository of the shapes for the duration of
  // the search.
  std::vector<ShapeSet> d_fragShapesPool;
  std::vector<std::pair<void *, ShapeSet *>> d_fragShapes;

  void extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  bool quickVerify(const SynthonSpaceHitSet *hitset,
                   const std::vector<size_t> &synthNums) const override;
  bool verifyHit(ROMol &hit) const override;
};
}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACESHAPESEARCHER_H
