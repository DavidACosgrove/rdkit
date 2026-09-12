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

#include "SubshapePoint.h"

#include <GraphMol/GaussianShape/SubshapeGrid.h>

namespace RDKit {
namespace GaussianShape {

SubshapePoint::SubshapePoint(
    const std::vector<std::unique_ptr<SubshapePoint>> &points,
    const boost::dynamic_bitset<> &pointsToUse) {
  for (size_t i = 0; i < points.size(); ++i) {
    if (pointsToUse.test(i)) {
      d_pos += points[i]->getPosition();
    }
  }
  d_pos /= pointsToUse.count();
}

}  // namespace GaussianShape
}  // namespace RDKit
