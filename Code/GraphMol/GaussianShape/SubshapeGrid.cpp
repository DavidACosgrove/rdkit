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

#include <limits>

#include <GraphMol/GaussianShape/AtomRadii.h>
#include <GraphMol/GaussianShape/SubshapeGrid.h>

#include <GraphMol/ROMol.h>

namespace RDKit {
namespace GaussianShape {

SubshapeGrid::SubshapeGrid(const ROMol &mol, int confId,
                           const SubshapeOptions &options)
    : d_gridSpacing(options.gridSpacing) {
  PRECONDITION(mol.getNumConformers() > 0,
               "Molecule needs coordinates for SubshapeGrid.");
  calcGridDimensions(mol, confId, options);
  std::cout << "Grid origin : " << d_gridOrigin[0] << ", " << d_gridOrigin[1]
            << ", " << d_gridOrigin[2] << std::endl;
  std::cout << "Grid size : " << d_gridSizes[0] << ", " << d_gridSizes[1]
            << ", " << d_gridSizes[2] << std::endl;
  d_occ_grid.resize(d_gridSizes[0] * d_gridSizes[1] * d_gridSizes[2],
                    SubshapeGridType::VOID);
  std::cout << "Occ grid size : " << d_occ_grid.size() << std::endl;
  fillGrid(mol, confId, options);
}

void SubshapeGrid::calcGridDimensions(const ROMol &mol, int confId,
                                      const SubshapeOptions &options) {
  double padding = d_gridSpacing + 2.0 * options.r_b;
  double xMin, xMax, yMin, yMax, zMin, zMax;
  xMin = yMin = zMin = std::numeric_limits<double>::max();
  xMax = yMax = zMax = std::numeric_limits<double>::lowest();
  const auto &conf = mol.getConformer(confId);
  for (const auto atom : mol.atoms()) {
    const auto pos = conf.getAtomPos(atom->getIdx());
    if (pos.x < xMin) {
      xMin = pos.x;
    }
    if (pos.x > xMax) {
      xMax = pos.x;
    }
    if (pos.y < yMin) {
      yMin = pos.y;
    }
    if (pos.y > yMax) {
      yMax = pos.y;
    }
    if (pos.z < zMin) {
      zMin = pos.z;
    }
    if (pos.z > zMax) {
      zMax = pos.z;
    }
  }
  d_gridOrigin[0] = xMin - padding;
  d_gridOrigin[1] = yMin - padding;
  d_gridOrigin[2] = zMin - padding;

  double xSide = xMax - xMin + 2.0 * padding;
  double ySide = yMax - yMin + 2.0 * padding;
  double zSide = zMax - zMin + 2.0 * padding;
  d_gridSizes[0] = int(xSide / d_gridSpacing) + 1;
  d_gridSizes[1] = int(ySide / d_gridSpacing) + 1;
  d_gridSizes[2] = int(zSide / d_gridSpacing) + 1;
}

void SubshapeGrid::fillGrid(const ROMol &mol, int confId,
                            const SubshapeOptions &options) {
  const auto &conf = mol.getConformer(confId);
  for (const auto atom : mol.atoms()) {
    const auto pos = conf.getAtomPos(atom->getIdx());
    // std::cout << pos.x << ' ' << pos.y << ' ' << pos.z << std::endl;
    auto xDisp = pos.x - d_gridOrigin[0];
    std::uint32_t xNum((d_gridSpacing * 0.5 + xDisp) / d_gridSpacing);
    auto yDisp = pos.y - d_gridOrigin[1];
    std::uint32_t yNum((d_gridSpacing * 0.5 + yDisp) / d_gridSpacing);
    auto zDisp = pos.z - d_gridOrigin[2];
    std::uint32_t zNum((d_gridSpacing * 0.5 + zDisp) / d_gridSpacing);
    // std::cout << "nearest point " << xNum << " " << yNum << " " << zNum
    // << std::endl;
    // std::cout << "at " << d_gridOrigin[0] + xNum * d_gridSpacing << ", "
    // << d_gridOrigin[1] + yNum * d_gridSpacing << ", "
    // << d_gridOrigin[2] + zNum * d_gridSpacing << std::endl;
#if 0
    size_t nearXNum, nearYNum, nearZNum;
    double nearX = std::numeric_limits<double>::max();
    double nearY = std::numeric_limits<double>::max();
    double nearZ = std::numeric_limits<double>::lowest();
    double closestDist = std::numeric_limits<double>::max();
    for (std::uint32_t i = 0; i < d_gridSizes[0]; ++i) {
      for (std::uint32_t j = 0; j < d_gridSizes[1]; ++j) {
        for (std::uint32_t k = 0; k < d_gridSizes[2]; ++k) {
          double x = d_gridOrigin[0] + i * d_gridSpacing;
          double y = d_gridOrigin[1] + j * d_gridSpacing;
          double z = d_gridOrigin[2] + k * d_gridSpacing;
          double dist_sq = (pos.x - x) * (pos.x - x) +
                           (pos.y - y) * (pos.y - y) +
                           (pos.z - z) * (pos.z - z);
          if (dist_sq < closestDist) {
            closestDist = dist_sq;
            nearX = x;
            nearY = y;
            nearZ = z;
            nearXNum = i;
            nearYNum = j;
            nearZNum = k;
          }
        }
      }
    }
    if (nearXNum != xNum || nearYNum != yNum || nearZNum != zNum) {
      std::cout << "nearest point " << xNum << " " << yNum << " " << zNum
                << std::endl;
      std::cout << "BF nearest point " << nearXNum << " " << nearYNum << " "
                << nearZNum << " at " << nearX << ", " << nearY << ", " << nearZ
                << " :: " << closestDist << " vs "
                << (d_gridOrigin[0] + xNum * d_gridSpacing - pos.x) *
                           (d_gridOrigin[0] + xNum * d_gridSpacing - pos.x) +
                       (d_gridOrigin[1] + yNum * d_gridSpacing - pos.y) *
                           (d_gridOrigin[1] + yNum * d_gridSpacing - pos.y) +
                       (d_gridOrigin[2] + zNum * d_gridSpacing - pos.z) *
                           (d_gridOrigin[2] + zNum * d_gridSpacing - pos.z)
                << std::endl
                << std::endl;
    }
#endif
    double rad = getStandardAtomRadius(atom->getAtomicNum());
    double extRad = rad + 2 * options.r_b;
    const std::uint32_t disp = 1 + extRad / d_gridSpacing;
    const std::uint32_t xStart = xNum - disp < 0 ? 0 : xNum - disp;
    const std::uint32_t yStart = yNum - disp < 0 ? 0 : yNum - disp;
    const std::uint32_t zStart = zNum - disp < 0 ? 0 : zNum - disp;
    double radSq = rad * rad;
    double innerRad = (rad + options.r_b) * (rad + options.r_b);
    double outerRad = (rad + 2.0 * options.r_b) * (rad + 2.0 * options.r_b);
    for (std::uint32_t i = xStart; i < xStart + 2 * disp; ++i) {
      for (std::uint32_t j = yStart; j < yStart + 2 * disp; ++j) {
        for (std::uint32_t k = zStart; k < zStart + 2 * disp; ++k) {
          double x = d_gridOrigin[0] + i * d_gridSpacing;
          double y = d_gridOrigin[1] + j * d_gridSpacing;
          double z = d_gridOrigin[2] + k * d_gridSpacing;
          double sqDist = (pos.x - x) * (pos.x - x) +
                          (pos.y - y) * (pos.y - y) + (pos.z - z) * (pos.z - z);
          auto index = getGridIndex(i, j, k);
          SubshapeGridType occ = d_occ_grid[index];
          if (sqDist <= radSq) {
            d_occ_grid[index] = SubshapeGridType::CORE;
          } else if (sqDist <= innerRad && occ < SubshapeGridType::INNER_SKIN) {
            d_occ_grid[index] = SubshapeGridType::INNER_SKIN;
          } else if (sqDist <= outerRad && occ < SubshapeGridType::OUTER_SKIN) {
            d_occ_grid[index] = SubshapeGridType::OUTER_SKIN;
          }
        }
      }
    }
  }
  int numCore = 0;
  int numInner = 0;
  int numOuter = 0;
  int numVoid = 0;
  for (size_t i = 0; i < d_occ_grid.size(); ++i) {
    switch (d_occ_grid[i]) {
      case SubshapeGridType::CORE:
        numCore++;
        break;
      case SubshapeGridType::INNER_SKIN:
        numInner++;
        break;
      case SubshapeGridType::OUTER_SKIN:
        numOuter++;
        break;
      case SubshapeGridType::VOID:
        numVoid++;
    }
  }
  std::cout << "numCore : " << numCore << " numInner : " << numInner
            << " numOuter : " << numOuter << " numVoid : " << numVoid
            << " total = " << numCore + numInner + numOuter + numVoid
            << "  total poss = "
            << d_gridSizes[0] * d_gridSizes[1] * d_gridSizes[2] << std::endl;
}
}  // namespace GaussianShape
}  // namespace RDKit