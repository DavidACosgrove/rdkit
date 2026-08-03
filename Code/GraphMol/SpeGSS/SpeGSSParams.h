// Copyright (C) David Cosgrove 2026
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#ifndef RDKIT_SPEGSSPARAMS_H
#define RDKIT_SPEGSSPARAMS_H

namespace RDKit {
namespace SpeGSS {

struct SpeGSSParams {
  double moleculePadding =
      4.5;  // Distance that grid extends beyond molecule.
            // Same as default distance cutoff for GaussianShape
  double gridSpacing =
      0.25;  // Spacing between grid points - x, y and z the same.
  double contourValue = 0.05;  // Contour value for surface.
};

}  // namespace SpeGSS
}  // namespace RDKit
#endif  // RDKIT_SPEGSSPARAMS_H
