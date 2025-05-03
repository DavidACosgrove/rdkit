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

  // Make a single ShapeInput from the given conformer number.
  // If confNum is out of range, use the first conformer.
  ShapeInput makeSingleShape(unsigned int confNum) const;

  void setActiveConformer(unsigned int confNum);

  std::string toString() const;

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int);
#endif

  unsigned int numDummies{0};
  double dummyVol{0.0};
  unsigned int actConf{0};
  std::vector<std::vector<float>> confCoords;
  std::vector<double> dummyVols;
  std::vector<double> sovs;
  std::vector<double> sofs;
};

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
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> BestSimilarity(
    SearchShapeInput &refShape, SearchShapeInput &fitShape,
    std::vector<float> &matrix, double threshold = 3.0, double opt_param = 1.0,
    unsigned int max_preiters = 10u, unsigned int max_postiters = 30u);

// Mock a molecule up from the shape for visual inspection.  No bonds.
// Atoms are C, features are N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<RDKit::RWMol> shapeToMol(
    const ShapeInput &shape);

#endif  // SEARCHSHAPEINPUT_H
