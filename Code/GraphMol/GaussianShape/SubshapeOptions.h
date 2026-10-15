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
  double r_w = 3.0;  // The window radius for finding the neighbours of atoms
                     // when building points.  The paper says 4.0, Greg
                     // uses 3.0.
  double clusterStartsThreshold =
      2.0;  // Clusters the start transformations at
            // this RMSD threshold so that similar start points
            // are not used.  < 0.0 means no clustering.
  double triangleSideTolerance = 1.0;  // When deciding if 2 triangles match so
                                       // are to be used as a start point, this
                                       // is the threshold in the difference in
                                       // side lengths.
  double rmsdTriangleTolerance = 0.5;  // When assessing whether two triangles
                                       // are overlaid well enough to be used
                                       // as start point, this is the RMSD
                                       // tolerance.
};

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SUBSHAPEOPTIONS_H
