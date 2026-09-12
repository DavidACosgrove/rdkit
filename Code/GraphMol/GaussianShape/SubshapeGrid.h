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

// This is an implementation of the occupancy grid used in
// Putta et al., JCICS, 2003, 43, 1623-1635
// A Novel Subshape Molecular Descriptor

#ifndef RDKIT_SUBSHAPEGRID_H
#define RDKIT_SUBSHAPEGRID_H

#include <ostream>
#include <vector>

#include <GraphMol/GaussianShape/SubshapeOptions.h>
#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;

namespace GaussianShape {

enum class SubshapeGridType : unsigned char {
  CORE = 3,
  INNER_SKIN = 2,
  OUTER_SKIN = 1,
  VOID = 0,
};

inline std::ostream &operator<<(std::ostream &os, const SubshapeGridType &t) {
  switch (t) {
    case SubshapeGridType::CORE:
      os << "CORE";
      break;
    case SubshapeGridType::INNER_SKIN:
      os << "INNER_SKIN";
      break;
    case SubshapeGridType::OUTER_SKIN:
      os << "OUTER_SKIN";
      break;
    case SubshapeGridType::VOID:
      os << "VOID";
  }
  return os;
}

class RDKIT_GAUSSIANSHAPE_EXPORT SubshapeGrid {
 public:
  explicit SubshapeGrid(const ROMol &mol, int confId = -1,
                        const SubshapeOptions &options = SubshapeOptions());
  SubshapeGrid(const SubshapeGrid &other) = default;
  SubshapeGrid(SubshapeGrid &&other) = default;
  ~SubshapeGrid() = default;

  SubshapeGrid &operator=(SubshapeGrid &&other) = default;
  SubshapeGrid &operator=(const SubshapeGrid &other) = default;

 private:
  double d_gridSpacing;
  std::array<double, 3> d_gridOrigin{0.0, 0.0, 0.0};
  std::array<std::uint32_t, 3> d_gridSizes{0, 0, 0};
  std::vector<SubshapeGridType> d_occ_grid;

  void calcGridDimensions(const ROMol &mol, int confId,
                          const SubshapeOptions &options);
  void fillGrid(const ROMol &mol, int confId, const SubshapeOptions &options);

  std::uint32_t getGridIndex(const std::uint32_t i, const std::uint32_t j,
                             const std::uint32_t k) const;
};

inline std::uint32_t SubshapeGrid::getGridIndex(const std::uint32_t i,
                                                const std::uint32_t j,
                                                const std::uint32_t k) const {
  return i * (d_gridSizes[1] * d_gridSizes[2]) + j * d_gridSizes[2] + k;
}

}  // namespace GaussianShape
}  // namespace RDKit
#endif  // RDKIT_SUBSPACEGRID_H
