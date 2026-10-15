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

// This is an implementation of the method of Putta et al,
// JCICS, 2003, 43, 1623-1635,
// 'A Novel Subshape Molecular Descriptor' for generating
// start points for passing into the shape overlay.  It
// generates the triangles for the 2 molecules and then
// transformations that orient the fit molecule or fit
// ShapeInput onto the reference based on the matching
// triangles.

#ifndef RDKIT_SUBSHAPESTARTS_H
#define RDKIT_SUBSHAPESTARTS_H

#include "Geometry/Transform3D.h"
#include "Geometry/UniformGrid3D.h"
#include "GraphMol/GaussianShape/SubshapeOptions.h"
#include "GraphMol/GaussianShape/SubshapePoint.h"

#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;

namespace GaussianShape {

class RDKIT_GAUSSIANSHAPE_EXPORT SubshapeStarts {
 public:
  SubshapeStarts(const ROMol &ref, const ROMol &fit, int refConfId = -1,
                 int fitConfId = -1,
                 const SubshapeOptions &options = SubshapeOptions());
  SubshapeStarts(const SubshapeStarts &other);
  SubshapeStarts(SubshapeStarts &&other) = default;
  SubshapeStarts &operator=(const SubshapeStarts &other);
  SubshapeStarts &operator=(SubshapeStarts &&other) = default;
  ~SubshapeStarts() = default;

  std::uint64_t getNumPossibleStarts() const { return d_possStarts.size(); }

  // To initialise the possible starts (valid triangle matches) iterations.
  void gotoFirstPossibleStart() { d_nextStart = 0; }
  // Cycle through the valid triangle matches, returning the transform
  // for the next match.  Returns an empty one at the end.  If ssd is
  // not null returns the sum of the squared distances in it, as calculated
  // by RDNumerics::Alignments::AlignPoints.  Since there are only 3 points
  // in the alignment, it is easy to change the tolerance to use this.
  std::unique_ptr<RDGeom::Transform3D> getNextStartTransform(
      double *ssd = nullptr);
  // Likewise, returning the transformation as a quaternion and translation
  // such as the SingleConformerAlignment likes - first 4
  // elements the quaternion, last 3 the translation.
  std::array<double, 7> getNextStartQuaternion(double *ssd = nullptr);

 private:
  std::vector<std::unique_ptr<SubshapePoint>> d_refPoints;
  std::vector<std::unique_ptr<SubshapePoint>> d_fitPoints;
  std::vector<std::vector<double>> d_refPointsDists;
  std::vector<std::vector<double>> d_fitPointsDists;
  std::uint64_t d_nextStart = 0;
  std::vector<std::array<std::uint64_t, 6>> d_possStarts;

  // Make the squared distance matrices for the points.
  void buildDistMatrices();

  // Cluster the start points at the threshold - only keep start
  // points that put the target start overlays more than the
  // threshold apart so we don't do an optimisation from a very
  // similar position (by RMSD) to one that has already been attempted.
  void clusterStartPoints(const ROMol &fit, int fitConfId,
                          const SubshapeOptions &options);
  // Sort the start points in ascending ssd, as returns by
  // getNextStartTransform so the ones that have the best
  // match between SubShapePoint triangles are at the start.
  void sortStartPoints();
};

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SUBSHAPESTARTS_H
