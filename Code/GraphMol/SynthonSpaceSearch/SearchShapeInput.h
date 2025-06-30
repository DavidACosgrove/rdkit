//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// A class derived from the pubchem-shape ShapeInput object with more
// stuff.

#ifndef SEARCHSHAPEINPUT_H
#define SEARCHSHAPEINPUT_H

#include <../External/pubchem_shape/PubChemShape.hpp>

// Make a subclass of ShapeInput with some extra info, including allowing
// for multiple conformations of the same atoms.  ShapeInput is in the
// global namespace, so this must be too.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SearchShapeInput : ShapeInput {
  SearchShapeInput() = default;
  SearchShapeInput(const std::string &str);

  SearchShapeInput(const ShapeInput &other);
  SearchShapeInput(const SearchShapeInput &other) = default;
  SearchShapeInput(SearchShapeInput &&other) = default;
  SearchShapeInput &operator=(const SearchShapeInput &other) = default;
  SearchShapeInput &operator=(SearchShapeInput &&other) = default;
  ~SearchShapeInput() override = default;

  // Merge the other ShapeInputSet, assuming it has the correct number
  // of atoms etc.  Unless this is empty, in which case it copies it
  // straight over.  Empties other, unless they can't be merged
  // in which case it returns unscathed.
  void merge(SearchShapeInput &other);

  bool hasNoShapes() const { return confCoords.empty(); }
  size_t getNumShapes() const { return confCoords.size(); }

  // Make a single ShapeInput from the given shape number.
  // If shapeNum is out of range, use the first shape.
  ShapeInput makeSingleShape(unsigned int shapeNum) const;

  void setActiveShape(unsigned int shapeNum);

  std::string toString() const;

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int);
#endif

  unsigned int numDummies{0};
  double dummyVol{0.0};
  unsigned int actConf{0};
  std::vector<std::vector<float>> confCoords;
  std::vector<unsigned int> molConfs;  // the conformer from the input molecule
                                       // that this shape refers to.  The shapes
                                       // are pruned and sorted so they may end
                                       // up not in the original order.
  std::vector<double> dummyVols;
  std::vector<double> sovs;
  std::vector<double> sofs;
  std::vector<std::vector<double>> shifts;
};

// Cut the shapes down so that none of them are more than the simThreshold
// in similarity with each other.
RDKIT_SYNTHONSPACESEARCH_EXPORT void pruneShapes(SearchShapeInput &shapes,
                                                 double simThreshold);

// Make a SearchShapeInput from all conformations of a molecule and then
// prune them at the given threshold. so that all selected shapes have
// a similarity score with each other that's less than the threshold.
// Assumes that shapeOpts.atomRadii only includes the dummy atoms, and
// doesn't use non-standard radii for other atoms.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<SearchShapeInput>
PrepareConformers(const RDKit::ROMol &mol,
                  const ShapeInputOptions &shapeOpts = ShapeInputOptions(),
                  double pruneThreshold = 1.9);

// Find the best similarity score between all the shapes in the 2
// SearchShapeInputs.  Stops as soon as it gets something above the
// threshold.  The threshold is applied to the sum of the shape
// tanimoto and color tanimoto, but the two values are returned separately.
// The similarity is between 0.0 and 2.0 so the default threshold of 3.0
// effectively means no threshold.
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> bestSimilarity(
    SearchShapeInput &refShape, SearchShapeInput &fitShape,
    std::vector<float> &matrix, double threshold = 3.0, double opt_param = 1.0,
    unsigned int max_preiters = 10u, unsigned int max_postiters = 30u);

// Mock a molecule up from the shape for visual inspection.  No bonds.
// Atoms are C, features are N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<RDKit::RWMol> shapeToMol(
    const ShapeInput &shape);

#endif  // SEARCHSHAPEINPUT_H
