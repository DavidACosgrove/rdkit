//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#ifndef RDKIT_SUBSHAPEOPTIONS_H
#define RDKIT_SUBSHAPEOPTIONS_H

namespace RDKit {
namespace GaussianShape {

struct SubshapeOptions {
  double gridSpacing = 0.5;
  double r_b = 0.25;  // half the thickness of the atom shell
  double n_max = 7;   // max neighbours for terminal points
  double r_w = 3.0;   // The window radius for the terminal points.  The paper
                      // says 4.0, Greg uses 3.0.
  bool sortBeforeClustering =
      false;  // The original algorithm clusters the terminal points
              // in the order they are produced.  The results thus
              // on the input atom order.  If this is true, they
              // are first sorted in descending number of neighbours
              // which gives a less prominent dependence on input order.
  unsigned int numTerminalPoints = 5;  // This is the ideal number.  We may get
                                       // fewer, but we need at least 3.
  double pointRadScale = 0.75;         // For the clustering
  double stepSize =
      1.0;  // The minimum distance of a skeleton point from a terminal point
  double maxDistC =
      15.0;  // A mysterious distance used when pruning skeleton points
  double symFactor =
      1.5;  // An another weird constant used when pruning skeleton points
};

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SUBSHAPEOPTIONS_H
