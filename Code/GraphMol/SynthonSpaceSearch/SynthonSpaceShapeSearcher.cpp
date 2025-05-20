//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <ranges>

#include <../External/pubchem_shape/PubChemShape.hpp>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearchHelpers.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <Geometry/point.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceShapeSearcher::SynthonSpaceShapeSearcher(
    const ROMol &query, const SynthonSpaceSearchParams &params,
    SynthonSpace &space)
    : SynthonSpaceSearcher(query, params, space) {
  if (space.getNumConformers() == 0) {
    throw std::runtime_error("No conformers found in SynthonSpaceSearch");
  }
  // For the fragmentation, we need to be able to keep track of the
  // original atoms indices.
  for (auto atom : query.atoms()) {
    atom->setProp<unsigned int>("ORIG_IDX", atom->getIdx());
  }
}

namespace {

// Take the fragged mol ShapeSets and flag all those synthons that have a
// fragment as a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<SearchShapeInput *> &fragShapes,
    const double similarityCutoff, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonSetOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction.getSynthons().size());

  // It makes sense to match fragments against synthon sets in order of
  // smallest synthon set first because if a fragment doesn't have a match
  // in a synthon set, the whole thing's a bust.  So if fragShapes[0] is matched
  // against 1000 synthons and then fragShapes[1] is matched against 10 synthons
  // and doesn't match any of them, the first set of matches was wasted time.
  std::vector<std::pair<unsigned int, size_t>> fragOrders(
      synthonSetOrder.size());
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    fragOrders[i].first = i;
    fragOrders[i].second = reaction.getSynthons()[synthonSetOrder[i]].size();
  }
  std::ranges::sort(fragOrders, [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  std::vector<float> matrix(12, 0.0);
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    const auto fragNum = fragOrders[i].first;
    const auto &synthons = reaction.getSynthons()[synthonSetOrder[fragNum]];
    // std::cout << fragNum << " vs " << synthonOrder[fragNum] << std::endl;
    // Get the smallest fragment volume.
    bool fragMatched = false;
    // Because the combination score is the sum of 2 tanimotos, it's not
    // possible to use the threshold to set upper and lower bonds on the
    // search space as is done with fingerprints and Rascal similarity.
    // We just need to plough through them in order.
    for (size_t j = 0; j < synthons.size(); j++) {
      // std::cout << "synthon " << j << " : " << synthons[j].first << " : "
      //           << synthons[j].second->getSmiles() << " vs frag "
      //           << synthonOrder[i] << std::endl;
      if (!synthons[j].second->getShapes() ||
          synthons[j].second->getShapes()->hasNoShapes()) {
        continue;
      }
      if (const auto sim = BestSimilarity(*fragShapes[fragNum],
                                          *synthons[j].second->getShapes(),
                                          matrix, similarityCutoff);
          sim.first + sim.second >= similarityCutoff) {
        // std::cout << "sim : " << sim.first + sim.second << " : " << sim.first
        // << ", " << sim.second << std::endl;
        synthonsToUse[synthonSetOrder[fragNum]][j] = true;
        fragSims[synthonSetOrder[fragNum]].emplace_back(j,
                                                        sim.first + sim.second);
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      // std::cout << "No matches" << std::endl;
      // No synthons matched this fragment, so the whole fragment set is a
      // bust.
      return retSynthons;
    }
  }

  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now order the synthons in descending order of their similarity to
  // the corresponding fragment.
  for (size_t i = 0; i < fragShapes.size(); i++) {
    if (fragSims[i].empty()) {
      // This one will have been filled in by expandBitSet so we need to use
      // all the synthons and a dummy similarity.
      fragSims[i].resize(synthonsToUse[i].size());
      for (size_t j = 0; j < fragSims[i].size(); j++) {
        fragSims[i][j] = std::make_pair(j, 0.0);
      }
    } else {
      std::ranges::sort(
          fragSims[i].begin(), fragSims[i].end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
    }
    retSynthons[i].clear();
    std::ranges::transform(fragSims[i], std::back_inserter(retSynthons[i]),
                           [](const auto &fs) { return fs.first; });
  }

  return retSynthons;
}

}  // namespace
std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceShapeSearcher::searchFragSet(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;

  // if (fragSet.size() != 3) {
  //   return results;
  // }
  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }
  // if (MolToSmiles(*fragSet[0]) != "[1*]NC1C[C@@H](C)CN(C(=O)c([2*])[3*])C1")
  // {
  //   return results;
  // }
  // std::cout << "\nSearchFragSet : ";
  // for (const auto &f : fragSet) {
  //   std::cout << MolToSmiles(*f) << " ";
  // }
  // std::cout << std::endl;
  // Collect the ShapeSets for the fragSet
  std::vector<SearchShapeInput *> fragShapes;
  fragShapes.reserve(fragSet.size());
  for (auto &frag : fragSet) {
    std::pair<void *, ShapeSet *> tmp{frag.get(), nullptr};
    const auto it = std::ranges::lower_bound(
        d_fragShapes, tmp, [](const auto &p1, const auto &p2) -> bool {
          return p1.first > p2.first;
        });
    fragShapes.push_back(it->second);
  }

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  const auto synthConnPatts = reaction.getSynthonConnectorPatterns();

  // Get all the possible permutations of connector numbers compatible with
  // the number of synthon sets in this reaction.  So if the
  // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
  // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
  // that might be how they're labelled in the reaction database.
  const auto connCombConnPatterns =
      details::getConnectorPermutations(connPatterns, reaction.getConnectors());

  // Need to try all combinations of synthon orders.
  const auto synthonOrders =
      details::permMFromN(fragSet.size(), reaction.getSynthons().size());

  for (const auto &synthonOrder : synthonOrders) {
    for (auto &connCombPatt : connCombConnPatterns) {
      // Make sure that for this connector combination, the synthons in this
      // order have something similar.  All query fragment connectors must
      // match something in the corresponding synthon.  The synthon can
      // have unused connectors.
      bool skip = false;
      for (size_t i = 0; i < connCombPatt.size(); ++i) {
        if ((connCombPatt[i] & synthConnPatts[synthonOrder[i]]).count() <
            connCombPatt[i].count()) {
          skip = true;
          break;
        }
      }
      if (skip) {
        continue;
      }
      auto theseSynthons = getHitSynthons(
          fragShapes,
          getParams().similarityCutoff - getParams().fragSimilarityAdjuster,
          reaction, synthonOrder);
      // std::cout << "Num hit synthons : " << theseSynthons.size() <<
      // std::endl;
      if (!theseSynthons.empty()) {
        std::unique_ptr<SynthonSpaceHitSet> hs(
            new SynthonSpaceHitSet(reaction, theseSynthons, fragSet));
        if (hs->numHits) {
          results.push_back(std::move(hs));
        }
      }
    }
  }

  return results;
}

namespace {
std::unique_ptr<SearchShapeInput> generateShapes(const ROMol &queryConfs,
                                                 const ROMol &frag,
                                                 double pruneThreshold) {
  // The fragSets molecules will have their atoms labelled with
  // ORIG_IDX apart from the dummy atoms, but we need coords
  // for them, too.  They are normally copied from the atom at
  // the other end of the broken bond so find that atom too.
  std::vector<unsigned int> fragAtoms;
  fragAtoms.reserve(frag.getNumAtoms());
  boost::dynamic_bitset<> inFrag(queryConfs.getNumAtoms());
  // std::cout << "\nGenerate shapes for " << MolToSmiles(frag) << " of "
  //           << queryConfs.getNumConformers() << " conformers" << std::endl;
  for (auto atom : frag.atoms()) {
    unsigned int origIdx;
    if (atom->getPropIfPresent<unsigned int>("ORIG_IDX", origIdx)) {
      fragAtoms.emplace_back(origIdx);
      inFrag[origIdx] = true;
    }
  }
  std::ranges::sort(fragAtoms);
  fragAtoms.erase(std::unique(fragAtoms.begin(), fragAtoms.end()),
                  fragAtoms.end());
  std::vector<std::pair<unsigned int, double>> dummyRadii;
  std::vector<unsigned int> notColorAtoms;
  for (auto atom : frag.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->getIsotope() >= 1 &&
        atom->getIsotope() <= MAX_CONNECTOR_NUM) {
      auto nbr = *frag.atomNeighbors(atom).begin();
      auto origNbr =
          queryConfs.getAtomWithIdx(nbr->getProp<unsigned int>("ORIG_IDX"));
      for (auto nbrNbr : queryConfs.atomNeighbors(origNbr)) {
        if (!inFrag[nbrNbr->getIdx()]) {
          dummyRadii.emplace_back(nbrNbr->getIdx(), 2.16);
          notColorAtoms.emplace_back(nbrNbr->getIdx());
        }
      }
    }
  }
  std::ranges::sort(dummyRadii);
  dummyRadii.erase(std::unique(dummyRadii.begin(), dummyRadii.end()),
                   dummyRadii.end());
  std::transform(dummyRadii.begin(), dummyRadii.end(),
                 std::back_inserter(fragAtoms),
                 [](const auto &p) -> unsigned int { return p.first; });

  // Build shapes with and without dummies to get a value for the dummy
  // atom volume in each conformation.
  ShapeInputOptions opts;
  opts.atomSubset = fragAtoms;
  opts.atomRadii = dummyRadii;
  opts.notColorAtoms = notColorAtoms;
  auto shapes = PrepareConformers(queryConfs, opts, pruneThreshold);
  return shapes;
}

void generateSomeShapes(
    const std::vector<ROMol *> &fragsForShape, unsigned int beginFrag,
    unsigned int endFrag, const ROMol &queryMolHs, double pruneThreshold,
    std::vector<std::unique_ptr<SearchShapeInput>> &fragShapes) {
  if (beginFrag >= fragsForShape.size()) {
    return;
  }
  if (endFrag >= fragsForShape.size()) {
    endFrag = fragsForShape.size();
  }
  for (unsigned int fragIdx = beginFrag; fragIdx < endFrag; ++fragIdx) {
    fragShapes[fragIdx] =
        generateShapes(queryMolHs, *fragsForShape[fragIdx], pruneThreshold);
    if (ControlCHandler::getGotSignal()) {
      return;
    }
  }
}
}  // namespace

void SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  // Use the given conformers if there are some unless it looks like a
  // 2D molecule.
  ShapeInputOptions opts;
  auto queryMol = std::unique_ptr<RWMol>(new RWMol(getQuery()));
  if (!queryMol->getNumConformers() || !queryMol->getConformer().is3D()) {
    std::cout << "Making query conformers" << std::endl;
    MolOps::addHs(*queryMol);
    auto dgParams = DGeomHelpers::ETKDGv3;
    dgParams.numThreads = getParams().numThreads;
    dgParams.pruneRmsThresh = getParams().confRMSThreshold;
    dgParams.randomSeed = getParams().randomSeed;
    // Build shapes for multiple conformations of the query molecule.
    auto cids = DGeomHelpers::EmbedMultipleConfs(
        *queryMol, getParams().numConformers, dgParams);
    if (cids.empty()) {
      BOOST_LOG(rdWarningLog)
          << "Couldn't generate conformers for query molecule." << std::endl;
      return;
    }
    MolOps::removeHs(*static_cast<RWMol *>(queryMol.get()));
  }
  std::cout << "Generating  query shapes" << std::endl;
  dp_queryShapes = PrepareConformers(*queryMol, opts, 1.9);
  // std::cout << "Query mol : " << MolToSmiles(getQuery())
  //           << " num confs = " << queryMolHs->getNumConformers() << " :: ";
  // for (auto a : getQuery().atoms()) {
  //   std::cout << a->getIdx() << ", " << a->getProp<unsigned int>("ORIG_IDX")
  //             << ", " << a->getAtomicNum() << " : ";
  // }
  // std::cout << "  num shapes : " << dp_queryShapes->confCoords.size()
  //           << std::endl;
  // std::cout << MolToCXSmiles(*queryMolHs) << std::endl;

  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return;
  }

  // Compute ShapeSets for the fragments
  std::cout << "Making shapes for fragments" << std::endl;
  d_fragShapesPool.resize(fragSmiToFrag.size());
  std::vector<ROMol *> fragsForShape;
  fragsForShape.reserve(fragSmiToFrag.size());
  std::transform(fragSmiToFrag.begin(), fragSmiToFrag.end(),
                 back_inserter(fragsForShape),
                 [](const auto &p) -> ROMol * { return p.second.front(); });

  unsigned int fragNum = 0;
  if (const auto numThreads = getNumThreadsToUse(getParams().numThreads);
      numThreads > 1) {
    std::cout << "numThreads: " << numThreads << std::endl;
    const size_t eachThread = 1 + fragsForShape.size() / numThreads;
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(generateSomeShapes, std::ref(fragsForShape),
                                    start, start + eachThread,
                                    std::ref(*queryMol), 1.9,
                                    std::ref(d_fragShapesPool)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    generateSomeShapes(fragsForShape, 0, fragsForShape.size(), *queryMol, 1.9,
                       d_fragShapesPool);
  }
  // Use the pooled ShapeSets to populate the vectors for each fragSet
  fragNum = 0;
  d_fragShapes.reserve(fragSmiToFrag.size());
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    for (auto &frag : frags) {
      d_fragShapes.emplace_back(frag, d_fragShapesPool[fragNum].get());
    }
    ++fragNum;
  }
  std::ranges::sort(d_fragShapes, [](const auto &p1, const auto &p2) -> bool {
    return p1.first > p2.first;
  });
  std::cout << "Done extra setup" << std::endl;
}

bool SynthonSpaceShapeSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  double maxVol = 0.0;
  double featureVol = 0.0;
  // The synthon shapes are sorted in descending order of sov + sof.
  // Assume therefore that the maximum volume of the hit is the sum
  // of the sov's of the first shape in each synthon, minus the volume
  // of their dummy atoms.
  for (size_t i = 0; i < synthNums.size(); i++) {
    const auto &synth = hitset->synthonsToUse[i][synthNums[i]].second;
    const auto &shapes = synth->getShapes();
    maxVol += shapes->sovs.front() - shapes->dummyVols.front();
    featureVol += shapes->sofs.front();
    // std::cout << "quickVerify synth " << i << " : " << synth->getSmiles()
    //           << " : " << shapes->sovs.front() << " and "
    //           << shapes->dummyVols.front() << std::endl;
  }
  double maxSt = std::min(maxVol, dp_queryShapes->sovs.front()) /
                 std::max(maxVol, dp_queryShapes->sovs.front());
  double maxCt = std::min(featureVol, dp_queryShapes->sofs.front()) /
                 std::max(featureVol, dp_queryShapes->sofs.front());
  // std::cout << maxSt << " " << maxCt << " : " << maxSt + maxCt << " : "
  //           << (maxSt + maxCt >= getParams().similarityCutoff -
  //                                    getParams().approxSimilarityAdjuster)
  //           << std::endl;
  return maxSt + maxCt >=
         getParams().similarityCutoff - getParams().approxSimilarityAdjuster;
}

bool SynthonSpaceShapeSearcher::verifyHit(ROMol &hit) const {
  auto dgParams = DGeomHelpers::ETKDGv3;
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = getParams().confRMSThreshold;
  dgParams.randomSeed = getParams().randomSeed;
  std::unique_ptr<ROMol> hitMolHs(MolOps::addHs(hit));
  DGeomHelpers::EmbedMultipleConfs(*hitMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*static_cast<RWMol *>(hitMolHs.get()));
  std::vector<float> matrix(12, 0.0);
  // std::cout << "Verifying hit for " << MolToSmiles(hit) << std::endl;
  for (size_t i = 0U; i < dp_queryShapes->confCoords.size(); ++i) {
    dp_queryShapes->setActiveConformer(i);
    for (unsigned int j = 0u; j < hitMolHs->getNumConformers(); ++j) {
      // std::cout << "Checking conf " << j << " against query conf " << i <<
      // std::endl;
      auto [st, ct] = AlignMolecule(*dp_queryShapes, *hitMolHs, matrix, j);
      if (st + ct >= getParams().similarityCutoff) {
        // std::cout << "sims : " << st << ", " << ct << " : " << st + ct
        // << std::endl;
        hit.setProp<double>("Similarity", st + ct);
        hit.setProp<unsigned int>("Query_Conformer", i);
        hit.addConformer(new Conformer(hitMolHs->getConformer(j)));
        return true;
      }
    }
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch