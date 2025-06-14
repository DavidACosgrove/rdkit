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
#include <GraphMol/SynthonSpaceSearch/SearchShapeInput.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>

namespace RDKit::SynthonSpaceSearch {

// A hash function used to hash a pair of any kind.  Taken from
// https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
// But with the hash combining from
// https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes/27952689#27952689
// Which is apparently effectively how boost::hash_combine does it.
struct hash_address_pair {
  size_t operator()(const std::pair<const void *, const void *> &p) const {
    // Hash the first element
    size_t hash1 = std::hash<const void *>{}(p.first);
    // Hash the second element
    size_t hash2 = std::hash<const void *>{}(p.second);
    // Combine the two hash values
    return hash1 ^ (hash2 + 0x517cc1b727220a95 + (hash1 << 6) + (hash1 >> 2));
  }
};

// Used for storing the pre-computed matches between fragments and
// synthons.  The actual type of the first address in the pair will vary,
// the second should always be Synthon *.  There's an entry only if
// the fragment->synthon similarity exceeded the threshold.
using FragSynthonSims =
    std::unordered_map<std::pair<const void *, const void *>, float,
                       hash_address_pair>;

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

  // Use d_fragSynthonSims to decide if the fragment matched the
  // Synthon.
  bool fragMatchedSynthon(const void *frag, const void *synthon,
                          float &sim) const;
  bool hasPrecomputedSims() const { return !d_fragSynthonSims.empty(); }

 private:
  // Shapes for all the conformers of the query.
  std::unique_ptr<SearchShapeInput> dp_queryShapes;
  // If a conformational expansion was done, keep it here, otherwise
  // just copy the query.
  std::unique_ptr<RWMol> dp_queryConfs;
  // These are the fragment shapes for this search, derived from
  // d_query.  The shapes in d_fragShapes are keyed on the address
  // of the corresponding fragment.  d_fragShapesPool is never read,
  // it is just used a repository of the shapes for the duration of
  // the search.
  std::vector<std::unique_ptr<SearchShapeInput>> d_fragShapesPool;
  std::vector<std::pair<void *, SearchShapeInput *>> d_fragShapes;

  // Precomputed similarities between fragments and synthons.  This
  // speeds things up because synthons are re-used and for the shape
  // search we have to calculate the similarity between every query
  // fragment and all the synthons in each SynthonSet.  Fingerprint
  // and RASCAL searching can dismiss a lot of fragment/synthon pairs
  // as not matching without calculating the full similarity so doing
  // the full set of similarity calculations up front slows things
  // down for them.
  FragSynthonSims d_fragSynthonSims;

  bool extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  bool quickVerify(const SynthonSpaceHitSet *hitset,
                   const std::vector<size_t> &synthNums) const override;
  bool verifyHit(ROMol &hit) const override;

  // Fill in the d_fragSynthonSims map.
  bool computeFragSynthonSims();
};
}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACESHAPESEARCHER_H
