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
#include <Geometry/Transform3D.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SynthonSpaceSearch/ProgressBar.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearchHelpers.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

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
  for (auto atom : getQuery().atoms()) {
    atom->setProp<unsigned int>("ORIG_IDX", atom->getIdx());
  }
}

bool SynthonSpaceShapeSearcher::fragMatchedSynthon(const void *frag,
                                                   const void *synthon,
                                                   float &sim) const {
  auto it = d_fragSynthonSims.find(std::make_pair(frag, synthon));
  if (it == d_fragSynthonSims.end()) {
    sim = -1.0;
    return false;
  }
  sim = it->second;
  return true;
}

namespace {

// Take the fragged mol ShapeSets and flag all those synthons that have a
// fragment as a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<SearchShapeInput *> &fragShapes,
    const double similarityCutoff, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonSetOrder,
    const SynthonSpaceShapeSearcher &shapeSearcher) {
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
  float sim;
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    const auto fragNum = fragOrders[i].first;
    const auto &synthons = reaction.getSynthons()[synthonSetOrder[fragNum]];
    // Get the smallest fragment volume.
    bool fragMatched = false;
    // Because the combination score is the sum of 2 tanimotos, it's not
    // possible to use the threshold to set upper and lower bonds on the
    // search space as is done with fingerprints and Rascal similarity.
    // We just need to plough through them in order.
    for (size_t j = 0; j < synthons.size(); j++) {
      if (!synthons[j].second->getShapes() ||
          synthons[j].second->getShapes()->hasNoShapes()) {
        continue;
      }
      if (shapeSearcher.hasPrecomputedSims()) {
        if (shapeSearcher.fragMatchedSynthon(fragShapes[fragNum],
                                             synthons[j].second, sim)) {
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(j, sim);
          fragMatched = true;
        }
      } else {
        if (const auto sim = bestSimilarity(*fragShapes[fragNum],
                                            *synthons[j].second->getShapes(),
                                            matrix, similarityCutoff);
            sim.first + sim.second >= similarityCutoff) {
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(
              j, sim.first + sim.second);
          fragMatched = true;
        }
      }
    }
    if (!fragMatched) {
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
    retSynthons[i].reserve(fragSims[i].size());
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
  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }
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
          reaction, synthonOrder, *this);
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
    const TimePoint *endTime,
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
    if (!(fragIdx % 100) && details::checkTimeOut(endTime)) {
      return;
    }
  }
}
}  // namespace

bool SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime) {
  // Use the given conformers if there are some unless it looks like a
  // 2D molecule.  We assume that the steroisomer is defined.
  auto queryMol = std::unique_ptr<RWMol>(new RWMol(getQuery()));
  if (!queryMol->getNumConformers() || !queryMol->getConformer().is3D()) {
    if (details::hasUnspecifiedStereo(*queryMol) &&
        !getParams().enumerateUnspecifiedStereo) {
      BOOST_LOG(rdErrorLog)
          << "The query molecule has unspecified stereochemistry."
             "  You need either to correct this or set the option "
             "'enumerateUnspecifiedStereo' to true."
          << std::endl;
      return false;
    }
    auto dgParams = DGeomHelpers::ETKDGv3;
    dgParams.numThreads = getParams().numThreads;
    dgParams.pruneRmsThresh = getParams().confRMSThreshold;
    dgParams.randomSeed = getParams().randomSeed;
    dgParams.timeout = getParams().timeOut;
    // Make conformers for this molecule, but without generating
    // isomers.  If that was needed, it will already have been done.
    auto queryMols = details::generateIsomerConformers(
        *queryMol, getParams().numConformers, false, getParams().stereoEnumOpts,
        dgParams);
    if (queryMols.empty() || details::checkTimeOut(endTime)) {
      return false;
    }
    dp_queryConfs = std::move(queryMols.front());
  } else {
    dp_queryConfs = std::make_unique<RWMol>(getQuery());
  }
  BOOST_LOG(rdInfoLog) << "Generating query shapes for "
                       << dp_queryConfs->getNumConformers() << " conformers of "
                       << MolToSmiles(*dp_queryConfs) << std::endl;
  ShapeInputOptions opts;
  dp_queryShapes = PrepareConformers(*dp_queryConfs, opts, 1.9);
  dp_queryShapes->setActiveShape(0);

  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return false;
  }
  // Compute ShapeSets for the fragments
  d_fragShapesPool.resize(fragSmiToFrag.size());
  std::vector<ROMol *> fragsForShape;
  fragsForShape.reserve(fragSmiToFrag.size());
  std::transform(fragSmiToFrag.begin(), fragSmiToFrag.end(),
                 back_inserter(fragsForShape),
                 [](const auto &p) -> ROMol * { return p.second.front(); });

  unsigned int fragNum = 0;
  if (const auto numThreads = getNumThreadsToUse(getParams().numThreads);
      numThreads > 1) {
    const size_t eachThread = 1 + fragsForShape.size() / numThreads;
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(generateSomeShapes, std::ref(fragsForShape),
                                    start, start + eachThread,
                                    std::ref(*dp_queryConfs), 1.9, endTime,
                                    std::ref(d_fragShapesPool)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    generateSomeShapes(fragsForShape, 0, fragsForShape.size(), *dp_queryConfs,
                       1.9, endTime, d_fragShapesPool);
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

  // The search will require that all fragments are compared against all
  // synthons at least once.  It's likely that a fragment occurs more than once
  // in the fragment sets, so compute them all up front.  This makes maximum use
  // of the parallel environment, doesn't do anything that wouldn't be done at
  // some point and may prevent duplicated comparisons.
  if (!computeFragSynthonSims(endTime)) {
    return false;
  }

  return true;
}

namespace {
// Process a chunk of the shape->synthon similarity calculations.
void computeSomeFragSynthonSims(
    const std::vector<std::pair<SearchShapeInput *, Synthon *>> &toDo,
    std::atomic<std::int64_t> &nextToDo, float threshold, std::mutex &mtx,
    const TimePoint *endTime, FragSynthonSims fragSynthonSims,
    std::unique_ptr<ProgressBar> &pbar) {
  SearchShapeInput *fragShape;
  Synthon *synthon;
  std::vector<float> matrix(12, 0.0);
  int numProgs = 1000;
  int numTimes = 1;
  while (true) {
    std::int64_t thisPair = ++nextToDo;
    if (thisPair >= static_cast<std::int64_t>(toDo.size())) {
      return;
    }
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    --numTimes;
    if (!numTimes) {
      numTimes = 100;
      if (details::checkTimeOut(endTime)) {
        return;
      }
    }
    fragShape = toDo[thisPair].first;
    synthon = toDo[thisPair].second;
    if (synthon->getShapes() && !synthon->getShapes()->hasNoShapes()) {
      const auto sim = bestSimilarity(*fragShape, *synthon->getShapes().get(),
                                      matrix, threshold);

      if (sim.first + sim.second >= threshold) {
        std::pair<const void *, const void *> p{fragShape, synthon};
        std::unique_lock lock1{mtx};
        fragSynthonSims.insert(std::make_pair(p, sim.first + sim.second));
      }
      if (pbar) {
        --numProgs;
        if (!numProgs) {
          numProgs = 1000;
          pbar->increment(1000);
        }
      }
    }
  }
}

// Calculate shape->synthon similarities for this block of the full search
// space, in parallel if required.
void processShapeSynthonList(
    const std::vector<std::pair<SearchShapeInput *, Synthon *>> &toDo,
    float threshold, const TimePoint *endTime, FragSynthonSims fragSynthonSims,
    std::unique_ptr<ProgressBar> &pbar, unsigned int numThreadsToUse) {
  std::atomic<std::int64_t> nextToDo = -1;
  std::vector<std::thread> threads;
  std::mutex mtx;
  if (numThreadsToUse > 1) {
    for (unsigned int i = 0u; i < numThreadsToUse; ++i) {
      threads.emplace_back(computeSomeFragSynthonSims, std::ref(toDo),
                           std::ref(nextToDo), threshold, std::ref(mtx),
                           endTime, std::ref(fragSynthonSims), std::ref(pbar));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    computeSomeFragSynthonSims(toDo, nextToDo, threshold, std::ref(mtx),
                               endTime, fragSynthonSims, pbar);
  }
}
}  // namespace

bool SynthonSpaceShapeSearcher::computeFragSynthonSims(
    const TimePoint *endTime) {
  float threshold =
      getParams().similarityCutoff - getParams().fragSimilarityAdjuster;

  std::unique_ptr<ProgressBar> pbar;
  if (getParams().useProgressBar) {
    pbar.reset(new ProgressBar(
        getParams().useProgressBar,
        d_fragShapesPool.size() * getSpace().d_synthonPool.size()));
    std::cout << "Computing fragment/synthon shape similarities for "
              << d_fragShapesPool.size() << " fragments against "
              << getSpace().d_synthonPool.size() << " synthons with "
              << getNumThreadsToUse(getParams().numThreads) << " threads."
              << std::endl;
  }

  if (const auto numThreadsToUse = getNumThreadsToUse(getParams().numThreads);
      numThreadsToUse > 1) {
    std::vector<std::pair<SearchShapeInput *, Synthon *>> toDo;
    toDo.reserve(2500000);
    for (size_t i = 0U; i < d_fragShapesPool.size(); ++i) {
      for (size_t j = 0U; j < getSpace().d_synthonPool.size(); ++j) {
        toDo.push_back(
            std::make_pair(d_fragShapesPool[i].get(),
                           getSpace().d_synthonPool[j].second.get()));
        if (toDo.size() == 2500000) {
          processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims,
                                  pbar, numThreadsToUse);
          toDo.clear();
        }
      }
    }
    processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims, pbar,
                            numThreadsToUse);
  }
  return !ControlCHandler::getGotSignal();
}  // namespace

bool SynthonSpaceShapeSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  auto approxSim = approxSimilarity(hitset, synthNums);
  return approxSim >=
         getParams().similarityCutoff - getParams().approxSimilarityAdjuster;
}

double SynthonSpaceShapeSearcher::approxSimilarity(
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
    if (!shapes || shapes->hasNoShapes()) {
      return false;
    }
    maxVol += shapes->sovs.front() - shapes->dummyVols.front();
    featureVol += shapes->sofs.front();
  }
  double maxSt = std::min(maxVol, dp_queryShapes->sovs.front()) /
                 std::max(maxVol, dp_queryShapes->sovs.front());
  double maxCt = std::min(featureVol, dp_queryShapes->sofs.front()) /
                 std::max(featureVol, dp_queryShapes->sofs.front());
  return maxSt + maxCt;
}

bool SynthonSpaceShapeSearcher::verifyHit(ROMol &hit) const {
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  std::cout << "verifyHit : " << MolToSmiles(hit) << std::endl;
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = getParams().confRMSThreshold;
  dgParams.randomSeed = getParams().randomSeed;
  dgParams.timeout = getParams().timeOut;
  auto hitConfs =
      details::generateIsomerConformers(hit, getParams().numConformers, true,
                                        getParams().stereoEnumOpts, dgParams);
  bool foundHit = false;
  double bestSim = getParams().similarityCutoff;
  ShapeInputOptions opts;

  std::vector<float> matrix(12, 0.0);
  for (auto &isomer : hitConfs) {
    auto hitShapes = PrepareConformers(*isomer, opts, 1.9);
    for (size_t i = 0U; i < dp_queryShapes->getNumShapes(); ++i) {
      dp_queryShapes->setActiveShape(i);
      for (unsigned int j = 0u; j < hitShapes->getNumShapes(); ++j) {
        hitShapes->setActiveShape(j);
        auto [st, ct] = AlignShape(*dp_queryShapes, *hitShapes, matrix);
        if (st + ct >= bestSim) {
          hit.setProp<double>("Similarity", st + ct);
          hit.setProp<unsigned int>("Query_Conformer",
                                    dp_queryShapes->molConfs[i]);
          // Make a molecule of this query conformer, and add it to the hit.
          ROMol thisConf(*dp_queryConfs, false, dp_queryShapes->molConfs[i]);
          hit.setProp<std::string>("Query_CXSmiles", MolToCXSmiles(thisConf));
          // Copy the conformer into the hit.
          hit.clearConformers();
          auto hitConf =
              new Conformer(isomer->getConformer(hitShapes->molConfs[j]));
          auto hitShape = hitShapes->makeSingleShape(j);
          TransformConformer(dp_queryShapes->shift, matrix, hitShape, *hitConf);
          hit.addConformer(hitConf, true);
          MolOps::assignStereochemistryFrom3D(hit);
          // If we're only interested in whether there's a shape match, and
          // not in finding the best shape, we're done.
          if (!getParams().bestHit) {
            return true;
          }
          foundHit = true;
          bestSim = st + ct;
        }
      }
    }
  }
  return foundHit;
}

}  // namespace RDKit::SynthonSpaceSearch
