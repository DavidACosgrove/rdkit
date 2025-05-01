//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef SYNTHONSPACESEARCHHELPERS_H
#define SYNTHONSPACESEARCHHELPERS_H

#include <../External/pubchem_shape/PubChemShape.hpp>
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif
namespace RDKit::SynthonSpaceSearch {

// This the maximum number of connectors that we can deal with at the moment.
// In reality, there may be fewer than this.  However, the key limit is in
// The symbols used for the connectors in Enamine REAL etc.
const std::vector<std::string> CONNECTOR_SYMBOLS{"[U]", "[Np]", "[Pu]", "[Am]"};
constexpr unsigned int MAX_CONNECTOR_NUM{4};

struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceSearchParams {
  std::int64_t maxHits{1000};  // The maximum number of hits to return.
                               // Use -1 for no maximum.
  std::uint64_t maxNumFragSets{
      100000};  // The maximum number of fragment sets the query can
                // be broken into.  Big molecules will create huge
                // numbers of fragment sets that may cause excessive
                // memory use.  If the number of fragment sets hits this
                // number, fragmentation stops and the search results
                // will likely be incomplete.
  std::int64_t toTryChunkSize{2500000};  // For similarity searching, especially
  // fingerprint similarity, there can be a
  // very large number of possible hits to
  // screen which can use a lot of memory and
  // crash the program.  It will also be very
  // slow.  To alleviate the memory use, the
  // possible hits are processed in chunks.
  // This parameter sets the chunk size.

  std::int64_t hitStart{0};  // Sequence number of hit to start from.  So that
                             // you can return the next N hits of a search
                             // having already obtained N-1.
  bool randomSample{false};  // If true, returns a random sample of the hit
                             // hits, up to maxHits in number.
  int randomSeed{-1};        // Seed for random-number generator.  -1 means use
                             // a random seed (std::random_device).
  bool buildHits{true};  // If false, reports the maximum number of hits that
                         // the search could produce, but doesn't return them.
  int numRandomSweeps{10};  // The random sampling doesn't always produce the
                            // required number of hits in 1 go.  This parameter
                            // controls how many loops it makes to try and get
                            // the hits before giving up.
  double similarityCutoff{0.5};  // Similarity cutoff for returning hits by
                                 // fingerprint similarity.  The default is
                                 // appropriate for a Morgan fingerprint of
                                 // radius=2, it may need changing for other
                                 // fingerprint types.
  double fragSimilarityAdjuster{
      0.1};  // Similarity values for fragments are generally low
             // due to low bit densities.  For the fragment
             // matching, reduce the similarity cutoff
             // by this amount.  A higher number will give slower search
             // times, a lower number will give faster searches at the
             // risk of missing some hits.  The value you use should have
             // a positive correlation with your FOMO.
  double approxSimilarityAdjuster{
      0.1};  // The fingerprint search uses an approximate similarity method
             // before building a product and doing a final check.  The
             // similarityCutoff is reduced by this value for the approximate
             // check.  A lower value will give faster run times at the
             // risk of missing some hits.  The value you use should have a
             // positive correlation with your FOMO.  The default is
             // appropriate for Morgan fingerprints.  With RDKit fingerprints,
             // 0.05 is adequate, and higher than that has been seen to
             // produce long run times.
  unsigned int minHitHeavyAtoms{0};  // Minimum number of heavy atoms in a hit.
  unsigned int maxHitHeavyAtoms{0};  // Maximum number of heavy atoms in a hit.
  // 0 means no maximum.
  double minHitMolWt{0};  // Minimum molecular weight for a hit.
  double maxHitMolWt{0};  // Maximum molecular weight for a hit.  0.0 means
  // no maximum.
  unsigned int minHitChiralAtoms{
      0};  // Minimum number of chiral atoms in a hit.
  unsigned int maxHitChiralAtoms{0};  // Maximum number of chiral atoms in a
  // hit. 0 means no maximum.
  unsigned int numConformers{100};  // When doing a shape search, the number of
                                    // conformers to use for each molecule.
  double confRMSThreshold{1.0};  // When doing a shape search, the RMS threshold
                                 // to use when pruning conformers.  Passed
                                 // directly EmbedMultipleConfs.
  std::uint64_t timeOut{600};  // Maximum number of seconds to spend on a single
                               // search.  0 means no maximum.
  int numThreads = 1;  // The number of threads to use.  If > 0, will use that
  // number.  If <= 0, will use the number of hardware
  // threads plus this number.  So if the number of
  // hardware threads is 8, and numThreads is -1, it will
  // use 7 threads.
};

// Make a subclass of ShapeInput with some extra info
struct SearchShapeInput : ShapeInput {
  SearchShapeInput() = default;
  SearchShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss(str);
    boost::archive::text_iarchive ia(ss);
    ia &*this;
#endif
  }

  SearchShapeInput(const ShapeInput &other) : ShapeInput(other) {}
  SearchShapeInput(const SearchShapeInput &other) = default;
  SearchShapeInput(SearchShapeInput &&other) = default;
  SearchShapeInput &operator=(const SearchShapeInput &other) = default;
  SearchShapeInput &operator=(SearchShapeInput &&other) = default;
  ~SearchShapeInput() override = default;

  std::string toString() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa &*this;
    return ss.str();
#endif
  }

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    ar &boost::serialization::base_object<ShapeInput>(*this);
    ar & numDummies;
    ar & dummyVol;
  }
#endif

  unsigned int numDummies{0};
  double dummyVol{0.0};
};

using ShapeSet = std::vector<std::unique_ptr<SearchShapeInput>>;

}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACESEARCHHELPERS_H
