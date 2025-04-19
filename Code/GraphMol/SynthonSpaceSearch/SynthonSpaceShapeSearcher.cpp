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
  for (auto atom : query.atoms()) {
    atom->setProp<unsigned int>("ORIG_IDX", atom->getIdx());
  }
}

namespace {

// Returns -1.0 if no shapes match within the threshold.
double bestShapeMatch(const ShapeSet &shapeSet1, const ShapeSet &shapeSet2,
                      double similarityCutoff) {
  double bestSim = -1.0;
  std::vector<float> matrix(12, 0.0);
  for (const auto &shape1 : shapeSet1) {
    for (const auto &shape2 : shapeSet2) {
      double maxSim = std::min(shape1->sov, shape2->sov) /
                          std::max(shape1->sov, shape2->sov) +
                      std::min(shape1->sof, shape2->sof) /
                          std::max(shape1->sof, shape2->sof);
      // The best similarity achievable is the lower of the two self scores
      std::cout << shape1->sov << ", " << shape1->sof << ", "
                << shape1->coord.size() / 3 << " vs " << shape2->sov << ", "
                << shape2->sof << ", " << shape2->coord.size() / 3
                << " :: " << maxSim << std::endl;

      if (maxSim > similarityCutoff) {
        auto [sov, sof] = AlignShapes(*shape1, *shape2, matrix);
        std::cout << shape1->sov << ", " << shape1->sof << " vs " << shape2->sov
                  << ", " << shape2->sof << " and " << sov << ", " << sof
                  << std::endl;
        if (sov + sof > bestSim) {
          bestSim = sov + sof;
        }
      }
    }
  }
  std::cout << "bestSim : " << bestSim << std::endl;
  return bestSim;
}

// Take the fragged mol ShapeSets and flag all those synthons that have a
// fragment as a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<ShapeSet *> &fragShapes, const double similarityCutoff,
    const SynthonSet &reaction, const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction.getSynthons().size());

  std::cout << "getHitSynthons : " << std::endl;
  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthons = reaction.getSynthons()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthons.size(); j++) {
      std::cout << "synthon " << j << " : " << synthons[j].first << " : "
                << synthons[j].second->getSmiles() << std::endl;
      if (const auto sim =
              bestShapeMatch(*fragShapes[i], synthons[j].second->getShapes(),
                             similarityCutoff);
          sim >= similarityCutoff) {
        synthonsToUse[synthonOrder[i]][j] = true;
        fragSims[synthonOrder[i]].emplace_back(j, sim);
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      std::cout << "No matching fragment found" << std::endl;
      // No synthons matched this fragment, so the whole fragment set is a
      // bust.
      return retSynthons;
    }
  }
  // Now order the synthons in descending order of their similarity to
  // the corresponding fragFP.
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

  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }
  std::cout << "\nSearchFragSet : ";
  for (const auto &f : fragSet) {
    std::cout << MolToSmiles(*f) << " ";
  }
  std::cout << std::endl;
  // Collect the ShapeSets for the fragSet
  std::vector<ShapeSet *> fragShapes;
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
      std::cout << "Num hit synthons : " << theseSynthons.size() << std::endl;
      if (!theseSynthons.empty()) {
        std::unique_ptr<SynthonSpaceHitSet> hs(
            new SynthonSpaceFPHitSet(reaction, theseSynthons, fragSet));
        if (hs->numHits) {
          results.push_back(std::move(hs));
        }
      }
    }
  }

  return results;
}

namespace {
ShapeSet generateShapes(const ROMol &queryConfs, const ROMol &frag) {
  // The fragSets molecules will have their atoms labelled with
  // ORIG_IDX apart from the dummy atoms, but we need coords
  // for them, too.  They are normally copied from the atom at
  // the other end of the broken bond so find that atom too.
  std::vector<unsigned int> fragAtoms;
  ShapeSet shapes;
  fragAtoms.reserve(frag.getNumAtoms());
  boost::dynamic_bitset<> inFrag(queryConfs.getNumAtoms());
  std::cout << "\nGenerate shape for " << MolToSmiles(frag) << " of "
            << MolToSmiles(queryConfs) << std::endl;
  std::cout << "frag atoms : ";
  for (auto atom : frag.atoms()) {
    unsigned int origIdx;
    if (atom->getPropIfPresent<unsigned int>("ORIG_IDX", origIdx)) {
      fragAtoms.emplace_back(origIdx);
      std::cout << fragAtoms.back() << " ";
      inFrag[origIdx] = true;
    }
  }
  std::cout << std::endl;
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
          std::cout << "Dummy atom " << nbrNbr->getIdx() << std::endl;
          fragAtoms.emplace_back(nbrNbr->getIdx());
          dummyRadii.emplace_back(nbrNbr->getIdx(), 2.16);
        }
      }
    }
  }
  std::ranges::sort(fragAtoms);
  fragAtoms.erase(std::unique(fragAtoms.begin(), fragAtoms.end()),
                  fragAtoms.end());
  std::ranges::sort(dummyRadii);
  dummyRadii.erase(std::unique(dummyRadii.begin(), dummyRadii.end()),
                   dummyRadii.end());

  ShapeInputOptions opts;
  opts.atomSubset = fragAtoms;
  opts.atomRadii = dummyRadii;
  shapes.resize(queryConfs.getNumConformers());
  for (unsigned int k = 0u; k < queryConfs.getNumConformers(); ++k) {
    auto shape = PrepareConformer(queryConfs, k, opts);
    shapes[k] = std::make_unique<ShapeInput>(shape);
    std::cout << k << " : " << shape.sov << " : " << shape.sof << " : "
              << shape.coord.size() / 3 << std::endl;
  }
  return shapes;
}
}  // namespace

void SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  auto queryMolHs = std::unique_ptr<ROMol>(MolOps::addHs(getQuery()));
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = getParams().numThreads;
  dgParams.pruneRmsThresh = 1.0;
  // Build shapes for multiple conformations of the query molecule.
  DGeomHelpers::EmbedMultipleConfs(*queryMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*static_cast<RWMol *>(queryMolHs.get()));
  for (unsigned int i = 0u; i < queryMolHs->getNumConformers(); ++i) {
    std::unique_ptr<ShapeInput> shape(
        new ShapeInput(PrepareConformer(*queryMolHs, i)));
    d_queryShapes.emplace_back(std::move(shape));
  }
  std::cout << "Query mol : " << MolToSmiles(getQuery())
            << " num confs = " << queryMolHs->getNumConformers() << " :: ";
  for (auto a : getQuery().atoms()) {
    std::cout << a->getIdx() << ", " << a->getProp<unsigned int>("ORIG_IDX")
              << ", " << a->getAtomicNum() << " : ";
  }
  std::cout << std::endl;

  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return;
  }

  // Compute ShapeSets for the fragments
  d_fragShapesPool.resize(fragSmiToFrag.size());
  unsigned int fragNum = 0;
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    d_fragShapesPool[fragNum++] = generateShapes(*queryMolHs, *frags.front());
  }

  // Use the pooled ShapeSets to populate the vectors for each fragSet
  fragNum = 0;
  d_fragShapes.reserve(fragSmiToFrag.size());
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    for (auto &frag : frags) {
      d_fragShapes.emplace_back(frag, &d_fragShapesPool[fragNum]);
    }
    ++fragNum;
  }
  std::ranges::sort(d_fragShapes, [](const auto &p1, const auto &p2) -> bool {
    return p1.first > p2.first;
  });
}

bool SynthonSpaceShapeSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  return true;
}

bool SynthonSpaceShapeSearcher::verifyHit(ROMol &hit) const {
  auto dgParams = DGeomHelpers::ETKDGv3;
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = 1.0;
  auto hitMolHs = MolOps::addHs(hit);
  DGeomHelpers::EmbedMultipleConfs(*hitMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*hitMolHs);
  std::vector<float> matrix;
  for (const auto &qshape : d_queryShapes) {
    for (unsigned int i = 0u; i < hitMolHs->getNumConformers(); ++i) {
      auto sims = AlignMolecule(*qshape, *hitMolHs, matrix, i);
      if (sims.first + sims.second >= getParams().similarityCutoff) {
        hit.setProp<double>("Similarity", sims.first + sims.second);
        auto coords = hitMolHs->getConformer(i).getPositions();
        for (unsigned int j = 0u; j < coords.size(); ++j) {
          hit.getConformer().setAtomPos(j, coords[j]);
        }
        return true;
      }
    }
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch