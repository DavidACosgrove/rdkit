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

// #include <GraphMol/GaussianShape/SubshapeGrid.h>
#include <Geometry/UniformGrid3D.h>
#include <GraphMol/GaussianShape/SubshapeOptions.h>
#include <GraphMol/GaussianShape/SubshapePoint.h>

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

 private:
  std::unique_ptr<RDGeom::UniformGrid3D> d_refGrid;
  std::unique_ptr<RDGeom::UniformGrid3D> d_fitGrid;
  std::vector<std::unique_ptr<SubshapePoint>> d_refPoints;
  std::vector<std::unique_ptr<SubshapePoint>> d_fitPoints;
};

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SUBSHAPESTARTS_H
