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

#include <Geometry/point.h>

namespace RDKit::SynthonSpaceSearch {

// Make a subclass of ShapeInput with some extra info
struct SearchShapeInput : ShapeInput {
  SearchShapeInput() = default;
  SearchShapeInput(const ShapeInput &other) : ShapeInput(other) {}
  SearchShapeInput(const SearchShapeInput &other) = default;
  SearchShapeInput(SearchShapeInput &&other) = default;
  SearchShapeInput &operator=(const SearchShapeInput &other) = default;
  SearchShapeInput &operator=(SearchShapeInput &&other) = default;
  ~SearchShapeInput() = default;

  unsigned int numDummies{0};
  double dummyVol{0.0};
};

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
      // The best score achievable is when the smaller volume is entirely inside
      // the larger volume.  The Shape tanimoto  is the fraction of volume in
      // common.
      double maxSt = std::min(shape1->sov, shape2->sov) /
                     std::max(shape1->sov, shape2->sov);
      double maxCt = std::min(shape1->sof, shape2->sof) /
                     std::max(shape1->sof, shape2->sof);
      double maxSim = maxSt + maxCt;
      if (maxSim > similarityCutoff) {
        // We want the score, but want the shape in the same place afterwards.
        auto keepCoord = shape2->coord;
        auto [sov, sof] = AlignShape(*shape1, *shape2, matrix);
        shape2->coord = keepCoord;
        if (sov + sof > bestSim) {
          bestSim = sov + sof;
        }
      }
    }
  }
  std::cout << "bestSim = " << bestSim << std::endl << std::endl;
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

  std::cout << "getHitSynthons : ";
  for (auto so : synthonOrder) {
    std::cout << so << " ";
  }
  std::cout << " : " << similarityCutoff << std::endl;
  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthons = reaction.getSynthons()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthons.size(); j++) {
      std::cout << "synthon " << j << " : " << synthons[j].first << " : "
                << synthons[j].second->getSmiles() << " vs frag "
                << synthonOrder[i] << std::endl;
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

  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now order the synthons in descending order of their similarity to
  // the corresponding fragFP.
  std::cout << "fragSims : " << fragSims.size() << std::endl;
  for (const auto &frs : fragSims) {
    for (const auto &fr : frs) {
      std::cout << "fr " << fr.first << " : " << fr.second << " :: ";
    }
    std::cout << std::endl;
  }
  std::cout << "retSynthons : " << retSynthons.size() << std::endl;
  for (auto so : retSynthons) {
    for (auto r : so) {
      std::cout << r << " ";
    }
    std::cout << std::endl;
  }
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
    std::cout << MolToSmiles(*f) << " "
              << f->getProp<std::string>("_smilesAtomOutputOrder") << " ";
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
ShapeSet generateShapes(const ROMol &queryConfs, const ROMol &frag,
                        double pruneThreshold) {
  // The fragSets molecules will have their atoms labelled with
  // ORIG_IDX apart from the dummy atoms, but we need coords
  // for them, too.  They are normally copied from the atom at
  // the other end of the broken bond so find that atom too.
  std::vector<unsigned int> fragAtoms;
  fragAtoms.reserve(frag.getNumAtoms());
  boost::dynamic_bitset<> inFrag(queryConfs.getNumAtoms());
  std::cout << "\nGenerate shapes for " << MolToSmiles(frag) << " of "
            << queryConfs.getNumConformers() << " conformers" << std::endl;
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
          // fragAtoms.emplace_back(nbrNbr->getIdx());
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

  ShapeInputOptions opts;
  opts.atomSubset = fragAtoms;
  opts.atomRadii = dummyRadii;
  opts.notColorAtoms = notColorAtoms;
  ShapeSet shapes(queryConfs.getNumConformers());
  ShapeInputOptions noDummyOpts;
  noDummyOpts = opts;
  noDummyOpts.atomRadii.clear();
  noDummyOpts.includeDummies = false;
  for (unsigned int k = 0u; k < queryConfs.getNumConformers(); ++k) {
    auto shape = PrepareConformer(queryConfs, k, opts);
    auto noDummyShape = PrepareConformer(queryConfs, k, noDummyOpts);
    SearchShapeInput *ss = new SearchShapeInput(shape);
    ss->numDummies = dummyRadii.size();
    ss->dummyVol = shape.sov - noDummyShape.sov;
    shapes[k].reset(ss);
  }

  details::pruneShapes(shapes, pruneThreshold);
  std::ranges::sort(shapes, [](const auto &s1, const auto &s2) -> bool {
    return s1->sov + s1->sof > s2->sov + s2->sof;
  });
  return shapes;
}
}  // namespace

void SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  auto queryMolHs = std::unique_ptr<ROMol>(MolOps::addHs(getQuery()));
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = getParams().numThreads;
  dgParams.pruneRmsThresh = getParams().confRMSThreshold;
  // Build shapes for multiple conformations of the query molecule.
  DGeomHelpers::EmbedMultipleConfs(*queryMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*static_cast<RWMol *>(queryMolHs.get()));
  for (unsigned int i = 0u; i < queryMolHs->getNumConformers(); ++i) {
    std::unique_ptr<ShapeInput> shape(
        new ShapeInput(PrepareConformer(*queryMolHs, i)));
    d_queryShapes.emplace_back(std::move(shape));
  }
  details::pruneShapes(d_queryShapes, 1.9);
  std::cout << "Query mol : " << MolToSmiles(getQuery())
            << " num confs = " << queryMolHs->getNumConformers() << " :: ";
  for (auto a : getQuery().atoms()) {
    std::cout << a->getIdx() << ", " << a->getProp<unsigned int>("ORIG_IDX")
              << ", " << a->getAtomicNum() << " : ";
  }
  std::cout << "  num shapes : " << d_queryShapes.size() << std::endl;
  std::cout << MolToCXSmiles(*queryMolHs) << std::endl;

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
    d_fragShapesPool[fragNum++] =
        generateShapes(*queryMolHs, *frags.front(), 1.9);
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
  double maxVol = 0.0;
  // The synthon shapes are sorted in descending order of sov + sof.
  // Assume therefore that the maximum volume of the hit is the sum
  // of the sov's of the first shape in each synthon, minus the volume
  // of their dummy atoms.
  for (unsigned int i = 0u; i < synthNums.size(); ++i) {
  }
  return true;
}

bool SynthonSpaceShapeSearcher::verifyHit(ROMol &hit) const {
  auto dgParams = DGeomHelpers::ETKDGv3;
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  dgParams.numThreads = 1;
  // dgParams.pruneRmsThresh = 1.0;
  std::unique_ptr<ROMol> hitMolHs(MolOps::addHs(hit));
  DGeomHelpers::EmbedMultipleConfs(*hitMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*static_cast<RWMol *>(hitMolHs.get()));
  std::vector<float> matrix(12, 0.0);
  std::cout << "Verifying hit for " << MolToSmiles(hit) << std::endl;
  for (const auto &qshape : d_queryShapes) {
    for (unsigned int i = 0u; i < hitMolHs->getNumConformers(); ++i) {
      // std::cout << "Checking conf " << i << std::endl;
      auto [st, ct] = AlignMolecule(*qshape, *hitMolHs, matrix, i);
      // std::cout << "sims : " << st << ", " << ct << std::endl;
      if (st + ct >= getParams().similarityCutoff) {
        hit.setProp<double>("Similarity", st + ct);
        hit.addConformer(new Conformer(hitMolHs->getConformer(i)));
        // std::cout << "transferred" << std::endl;
        return true;
      }
    }
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch