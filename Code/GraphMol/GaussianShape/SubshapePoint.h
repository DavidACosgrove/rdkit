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

// Points for use in the SubspaceStarts code.

#ifndef RDKIT_SUBSHAPEPOINT_H
#define RDKIT_SUBSHAPEPOINT_H

#include <boost/dynamic_bitset.hpp>

#include "Geometry/point.h"

namespace RDKit {
class ROMol;
namespace GaussianShape {

class SubshapePoint {
 public:
  explicit SubshapePoint(const RDGeom::Point3D &pos) : d_pos(pos) {}
  // Make a new point at the average position of the points given.
  SubshapePoint(const std::vector<std::unique_ptr<SubshapePoint>> &points,
                const boost::dynamic_bitset<> &pointsToUse);
  SubshapePoint(const SubshapePoint &other) = default;
  SubshapePoint(SubshapePoint &&other) = default;
  SubshapePoint &operator=(const SubshapePoint &other) = default;
  SubshapePoint &operator=(SubshapePoint &&other) = default;
  ~SubshapePoint() = default;

  const RDGeom::Point3D &getPosition() const { return d_pos; }
  void setPosition(const RDGeom::Point3D &pos) { d_pos = pos; }
  const RDGeom::Point3D &getDirection() const { return d_dir; }
  void setDirection(const RDGeom::Point3D &dir) { d_dir = dir; }

 private:
  RDGeom::Point3D d_pos;  // position
  RDGeom::Point3D d_dir;  // direction, as defined in the paper
};

}  // namespace GaussianShape
}  // namespace RDKit
#endif  // RDKIT_SUBSHAPEPOINT_H
