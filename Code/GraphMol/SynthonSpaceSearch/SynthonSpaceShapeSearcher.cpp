//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <../External/pubchem_shape/PubChemShape.hpp>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>

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

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceShapeSearcher::searchFragSet(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;

  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }

  return results;
}

void SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  auto queryMolHs = std::unique_ptr<ROMol>(MolOps::addHs(getQuery()));
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = getParams().numThreads;

  // Build shapes for multiple conformations of the query molecule.
  DGeomHelpers::EmbedMultipleConfs(*queryMolHs, getParams().numConformers,
                                   dgParams);
  MolOps::removeHs(*static_cast<RWMol *>(queryMolHs.get()));
  for (unsigned int i = 0u; i < queryMolHs->getNumConformers(); ++i) {
    std::unique_ptr<ShapeInput> shape(
        new ShapeInput(PrepareConformer(*queryMolHs, i, true)));
    d_queryShapes.emplace_back(std::move(shape));
  }

  // Now compute ShapeInput objects for all the fragments in
  // fragSets.
  boost::dynamic_bitset<> fragAtoms(getQuery().getNumAtoms());
  d_fragShapes.resize(fragSets.size());
  for (size_t i = 0u; i < fragSets.size(); ++i) {
    d_fragShapes[i].resize(fragSets[i].size());
    for (size_t j = 0u; j < fragSets[i].size(); ++j) {
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(fragSets[i][j].get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);

      std::cout << i << ", " << j << " : " << fragSets[i][j]->getNumAtoms()
                << " : " << fragSets[i][j]->getNumConformers() << " : "
                << MolToSmiles(*fragSets[i][j]) << " :: " << std::endl;
      // The fragSets molecules will have their atoms labelled with
      // ORIG_IDX apart from the dummy atoms, but we need coords
      // for them, too.  They are normally copied from the atom at
      // the other end of the broken bond so find that atom too.
      std::vector<unsigned int> fragAtoms;
      fragAtoms.reserve(fragSets[i][j]->getNumAtoms());
      boost::dynamic_bitset<> inFrag(queryMolHs->getNumAtoms());
      for (auto atom : fragSets[i][j]->atoms()) {
        unsigned int origIdx;
        if (atom->getPropIfPresent<unsigned int>("ORIG_IDX", origIdx)) {
          fragAtoms.emplace_back(origIdx);
          std::cout << fragAtoms.back() << " ";
          inFrag[origIdx] = true;
        }
      }
      std::cout << std::endl;
      for (auto atom : fragSets[i][j]->atoms()) {
        if (atom->getAtomicNum() == 0 && atom->getIsotope() >= 1 &&
            atom->getIsotope() <= MAX_CONNECTOR_NUM) {
          auto nbr = *fragSets[i][j]->atomNeighbors(atom).begin();
          auto origNbr = queryMolHs->getAtomWithIdx(
              nbr->getProp<unsigned int>("ORIG_IDX"));
          for (auto nbrNbr : queryMolHs->atomNeighbors(origNbr)) {
            if (!inFrag[nbrNbr->getIdx()]) {
              fragAtoms.emplace_back(nbrNbr->getIdx());
            }
          }
        }
      }
      std::sort(fragAtoms.begin(), fragAtoms.end());
      fragAtoms.erase(std::unique(fragAtoms.begin(), fragAtoms.end()),
                      fragAtoms.end());

      // Now copy the fragment, make the dummies Fr atoms, copy the
      // coords from the conformers and create the shapes.
      ROMol fragCp(*fragSets[i][j]);
      std::cout << "fragCP num atoms : " << fragCp.getNumAtoms() << " : "
                << fragAtoms.size() << " :: ";
      for (auto fa : fragAtoms) {
        std::cout << fa << " ";
      }
      std::cout << std::endl;
      for (auto atom : fragCp.atoms()) {
        if (atom->getAtomicNum() == 0 && atom->getIsotope() >= 1 &&
            atom->getIsotope() <= MAX_CONNECTOR_NUM) {
          atom->setAtomicNum(87);
        }
      }
      d_fragShapes[i][j].resize(queryMolHs->getNumConformers());
      for (unsigned int k = 0u; k < queryMolHs->getNumConformers(); ++k) {
        const auto &wholeConf = queryMolHs->getConformer(k);
        Conformer *newConformer = new Conformer(fragAtoms.size());
        newConformer->set3D(true);
        for (size_t l = 0u; l < fragAtoms.size(); ++l) {
          newConformer->setAtomPos(j, wholeConf.getAtomPos(fragAtoms[l]));
        }
        std::cout << "conformer num atoms : " << newConformer->getNumAtoms()
                  << std::endl;
        fragCp.addConformer(newConformer, true);
        // d_fragShapes[i][j][k] = std::move(std::unique_ptr<ShapeInput>(
        //     new ShapeInput(PrepareConformer(fragCp, k, true))));
      }
    }
  }
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