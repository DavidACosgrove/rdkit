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

#include <algorithm>
#include <numbers>

#include <boost/dynamic_bitset.hpp>

#include "Geometry/GridUtils.h"
#include "GraphMol/ROMol.h"
#include <GraphMol/GaussianShape/SubshapeStarts.h>
#include <GraphMol/ShapeHelpers/ShapeEncoder.h>

namespace RDKit {
namespace GaussianShape {

namespace {

// This is all largely a transliteration of Greg's Python implementation
// in $RDBASE/rdkit/Chem/Subshape/BuilderUtils.py
void calcGridDimensions(const ROMol &mol, int confId,
                        const SubshapeOptions &options, RDGeom::Point3D &origin,
                        double dims[3]) {
  double padding = options.gridSpacing + 2.0 * options.r_w;
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
  origin.x = xMin - padding;
  origin.y = yMin - padding;
  origin.z = zMin - padding;
  dims[0] = xMax - xMin + 2.0 * padding;
  dims[1] = yMax - yMin + 2.0 * padding;
  dims[2] = zMax - zMin + 2.0 * padding;
}

std::unique_ptr<RDGeom::UniformGrid3D> buildGrid(
    const ROMol &mol, int confId, const SubshapeOptions &options) {
  RDGeom::Point3D origin;
  double dims[3];
  calcGridDimensions(mol, confId, options, origin, dims);
  std::cout << "grid dims : " << dims[0] << ", " << dims[1] << ", " << dims[2]
            << std::endl;
  std::cout << "origin : " << origin[0] << ", " << origin[1] << ", "
            << origin[2] << std::endl;
  auto grid = std::make_unique<RDGeom::UniformGrid3D>(
      dims[0], dims[1], dims[2], options.gridSpacing,
      RDKit::DiscreteValueVect::TWOBITVALUE, &origin);
  MolShapes::EncodeShape(mol, *grid, confId, nullptr, 0.8, options.r_b, -1,
                         true);
  int num[4] = {0, 0, 0, 0};
  const auto shapeVect = grid->getOccupancyVect();
  for (unsigned int i = 0; i < shapeVect->size(); ++i) {
    if ((*shapeVect)[i] == 3.0) {
      ++num[3];
    } else if ((*shapeVect)[i] == 2.0) {
      ++num[2];
    } else if ((*shapeVect)[i] == 1.0) {
      ++num[1];
    } else {
      ++num[0];
    }
  }
  std::cout << "nums : " << num[0] << ", " << num[1] << ", " << num[2] << ", "
            << num[3] << "  tot = " << num[0] + num[1] + num[2] + num[3]
            << std::endl;
  return grid;
}

void makeInitialTerminalPoints(
    const ROMol &mol, int confId, const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  const auto &conf = mol.getConformer(confId);
  double rWSq = options.r_w * options.r_w;
  std::cout << "makeInitialTerminalPoints at " << options.r_w << std::endl;
  std::vector<std::vector<unsigned int>> nbours(mol.getNumAtoms());

  for (const auto atom1 : mol.atoms()) {
    unsigned int at1Idx = atom1->getIdx();
    const auto pos1 = conf.getAtomPos(at1Idx);
    for (const auto atom2 : mol.atoms()) {
      unsigned int at2Idx = atom2->getIdx();
      if (at1Idx <= at2Idx) {
        continue;
      }
      const auto pos2 = conf.getAtomPos(at2Idx);
      double sqDist = (pos1 - pos2).lengthSq();
      if (sqDist < rWSq) {
        nbours[at1Idx].push_back(at2Idx);
        nbours[at2Idx].push_back(at1Idx);
      }
    }
  }
  for (size_t i = 0; i < nbours.size(); ++i) {
    // The n_max includes i itself, which isn't in the list.
    if (nbours[i].size() < options.n_max - 1) {
      // This isn't as described in the paper, because that didn't make any
      // sense, and indeed the algorithm in Figure 6 is different from
      // equation 3.  This follows what Greg did in his Python implementation
      // in $RDBASE/rdkit/Chem/Subshape/BuilderUtils.py
      auto pnt = conf.getAtomPos(i);
      double totWt = 1.0;  // For atom i
      for (const auto n : nbours[i]) {
        // Add one to the weights because Greg put i in the nbours lists, which
        // I haven't done.
        double wt =
            static_cast<double>(nbours[i].size() + 1) / (nbours[n].size() + 1);
        pnt += conf.getAtomPos(n) * wt;
        totWt += wt;
      }
      pnt /= totWt;
      termPoints.push_back(std::make_unique<SubshapePoint>(pnt));
      // std::cout << pnt.x << "," << pnt.y << "," << pnt.z << ";";
    }
  }
  // std::cout << " :: " << termPoints.size() << std::endl;
}

void clusterTerminalPoints(
    const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  const double radCutoffSq =
      options.pointRadScale * options.pointRadScale * options.r_w * options.r_w;
  std::vector<boost::dynamic_bitset<>> nbours(
      termPoints.size(), boost::dynamic_bitset<>(termPoints.size()));
  for (size_t i = 0; i < termPoints.size(); ++i) {
    nbours[i].set(i);
  }
  for (size_t i = 1; i < nbours.size(); ++i) {
    const auto &p1 = termPoints[i];
    for (size_t j = 0; j < i; ++j) {
      if ((p1->getPosition() - termPoints[j]->getPosition()).lengthSq() <
          radCutoffSq) {
        // std::cout
        //     << i << " -> " << j << " : "
        //     << (p1->getPosition() - termPoints[j]->getPosition()).length()
        //     << " : "
        //     << (p1->getPosition() - termPoints[j]->getPosition()).lengthSq()
        //     << " vs " << radCutoffSq << std::endl;
        nbours[i].set(j);
        nbours[j].set(i);
      }
    }
  }
  // The paper's clustering algorithm always starts with the first
  // terminal point which clearly depends on the order of atoms in
  // the molecule.  Try and make it a bit better by sorting in
  // descending order of number of neighbours.  There will still
  // be ties, though, which means some atom order dependence seems
  // unavoidable.
  if (options.sortBeforeClustering) {
    std::ranges::sort(nbours, [](const auto &a, const auto &b) -> bool {
      return a.count() > b.count();
    });
  }
  std::vector<std::unique_ptr<SubshapePoint>> clusteredPoints;
  while (!nbours.empty()) {
    clusteredPoints.emplace_back(
        std::make_unique<SubshapePoint>(termPoints, nbours.front()));
    // knock members of nbours.front() out of the other lists
    for (size_t i = 1; i < nbours.size(); ++i) {
      for (size_t j = 0; j < nbours.front().size(); ++j) {
        if (nbours.front()[j]) {
          nbours[i][j] = 0;
        }
      }
    }
    nbours.front().clear();
    std::erase_if(nbours, [](const auto &a) -> bool { return !a.count(); });
  }
  termPoints = std::move(clusteredPoints);
}

RDGeom::Point3D findFurthestGridPoint(const RDGeom::UniformGrid3D &grid,
                                      const RDGeom::Point3D &pt,
                                      const SubshapeOptions &options) {
  const auto &shapeVect = grid.getOccupancyVect();
  double distMaxSq = -1.0;
  RDGeom::Point3D furthestPt;
  for (size_t i = 0; i < shapeVect->size(); ++i) {
    if ((*shapeVect)[i] < 3.0) {
      continue;
    }
    const auto pos = grid.getGridPointLoc(i);
    double distSq = (pos - pt).lengthSq();
    if (distSq > distMaxSq) {
      distMaxSq = distSq;
      furthestPt = pos;
    }
  }
  double weightSum = 0.0;
  const auto res =
      RDGeom::computeGridCentroid(grid, furthestPt, options.r_w, weightSum);
  return res;
}

RDGeom::Point3D findGridPointBetweenPoints(const RDGeom::Point3D &pt1,
                                           const RDGeom::Point3D &pt2,
                                           const RDGeom::UniformGrid3D &grid,
                                           double windowRadius) {
  auto mid = (pt1 + pt2) / 2.0;
  double distSq = 1.0e16;
  double weightSum = 0.0;
  while (distSq > grid.getSpacing()) {
    const auto centroid =
        RDGeom::computeGridCentroid(grid, mid, windowRadius, weightSum);
    distSq = (mid - centroid).lengthSq();
    mid = centroid;
  }
  return mid;
}

void getMoreTerminalPoints(
    const RDGeom::UniformGrid3D &grid, const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  // Use a max-min algorithm to get enough points
  const auto &shapeVect = grid.getOccupancyVect();

  while (termPoints.size() < options.numTerminalPoints) {
    double maxMinSq = -1.0;
    RDGeom::Point3D bestPt;
    for (size_t i = 0; i < shapeVect->size(); ++i) {
      if ((*shapeVect)[i] < 3.0) {
        continue;
      }
      double minDistSq = 1.0e16;
      const auto posI = grid.getGridPointLoc(i);
      for (const auto &pt : termPoints) {
        double distSq = (pt->getPosition() - posI).lengthSq();
        if (distSq < minDistSq) {
          minDistSq = distSq;
        }
      }
      if (minDistSq > maxMinSq) {
        maxMinSq = minDistSq;
        bestPt = posI;
      }
    }
    double weightSum = 0.0;
    const auto centroid =
        RDGeom::computeGridCentroid(grid, bestPt, options.r_w, weightSum);
    std::cout << "centroid : " << centroid << "  weightSum = " << weightSum
              << std::endl;
    termPoints.emplace_back(std::make_unique<SubshapePoint>(centroid));
  }
}

void addExtraTerminalPoints(
    const RDGeom::UniformGrid3D &grid, const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  // This is taken from Greg's Python implementation
  if (termPoints.size() == 1) {
    // If there's only 1 point, add the one with max value that is
    // furthest from it.
    const auto pt =
        findFurthestGridPoint(grid, termPoints[0]->getPosition(), options);
    termPoints.emplace_back(std::make_unique<SubshapePoint>(pt));
  }
  if (termPoints.size() == 2) {
    // Add a point roughly in the middle.
    const auto pt1 = termPoints[0]->getPosition();
    const auto pt2 = termPoints[1]->getPosition();
    const auto mid = findGridPointBetweenPoints(pt1, pt2, grid, options.r_w);
    termPoints.emplace_back(std::make_unique<SubshapePoint>(mid));
  }
  if (termPoints.size() < options.numTerminalPoints) {
    getMoreTerminalPoints(grid, options, termPoints);
  }
}

void buildTerminalPoints(
    const ROMol &mol, int confId, const SubshapeOptions &options,
    const RDGeom::UniformGrid3D &grid,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  makeInitialTerminalPoints(mol, confId, options, termPoints);
  std::cout << "Initial number of terminal points: " << termPoints.size()
            << std::endl;
  clusterTerminalPoints(options, termPoints);
  std::cout << "Clustered number of terminal points: " << termPoints.size()
            << std::endl;
  if (termPoints.size() < options.numTerminalPoints) {
    addExtraTerminalPoints(grid, options, termPoints);
  }
}

void addSkeletonPoints(
    const RDGeom::UniformGrid3D &grid, const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  const auto shapeVect = grid.getOccupancyVect();
  const double stepSq = options.stepSize * options.stepSize;
  std::cout << "stepDist : " << options.stepSize << std::endl;
  // Build the initial set of skeleton points and add them to end
  // end of termpPoints, keeping track of where they started
  size_t numOrigTermPoints = termPoints.size();
  int num3 = 0;
  for (size_t i = 1; i < shapeVect->size(); ++i) {
    if ((*shapeVect)[i] < 3.0) {
      continue;
    }
    ++num3;
    const auto posI = grid.getGridPointLoc(i);
    bool ok = true;
    for (size_t j = 0; j < numOrigTermPoints; ++j) {
      double distSq = (termPoints[j]->getPosition() - posI).lengthSq();
      if (distSq < stepSq) {
        ok = false;
        break;
      }
    }
    if (ok) {
      termPoints.emplace_back(new SubshapePoint(posI));
    }
  }
  std::cout << "num3 : " << num3 << std::endl;
  std::cout << "Orig number of terminal points: " << numOrigTermPoints
            << " num added points : " << termPoints.size() - numOrigTermPoints
            << std::endl;
  // Now trim the skeleton points.  This came from Greg's code and I don't
  // entirely understand the maxDistC bit.
  std::vector<double> fracVols(termPoints.size());
  double weightSum = 0.0;
  const double gridBoxVolume =
      grid.getSpacing() * grid.getSpacing() * grid.getSpacing();
  // The 3.0 in the numerator is the grid value for fully occupied grid points.
  const double sphereVol =
      4.0 * std::numbers::pi * options.r_w * options.r_w * options.r_w / 3.0;
  const double maxVol = sphereVol * 3.0 / gridBoxVolume;
  const double maxDistCSq = options.maxDistC * options.maxDistC;

  for (size_t i = numOrigTermPoints; i < termPoints.size(); ++i) {
    const auto centroid = RDGeom::computeGridCentroid(
        grid, termPoints[i]->getPosition(), options.r_w, weightSum);
    double centroidPointDistSq =
        (termPoints[i]->getPosition() - centroid).lengthSq();
    if (centroidPointDistSq > maxDistCSq) {
      termPoints[i].reset();
    } else {
      fracVols[i] = weightSum / maxVol;
      termPoints[i]->setPosition(centroid);
    }
  }
  // Continue trimming
  const double sqDistCutoff = options.symFactor * options.symFactor *
                              options.stepSize * options.stepSize;
  for (size_t i = 0; i < termPoints.size(); ++i) {
    if (!termPoints[i]) {
      continue;
    }
    int p = -1;
    double mFrac = 0.0;
    const auto ptI = termPoints[i]->getPosition();

    int startJ = std::max(i + 1, numOrigTermPoints);
    for (int j = startJ; j < static_cast<int>(termPoints.size()); ++j) {
      if (!termPoints[j]) {
        continue;
      }
      const auto ptJ = termPoints[j]->getPosition();
      const double distCSq = (ptI - ptJ).lengthSq();
      if (distCSq < sqDistCutoff && fracVols[j] > mFrac) {
        p = j;
        mFrac = fracVols[j];
      }
    }
    if (p > -1) {
      for (int j = startJ; j < static_cast<int>(termPoints.size()); ++j) {
        if (j == p || !termPoints[j]) {
          continue;
        }
        const auto ptJ = termPoints[j]->getPosition();
        double distCSq = (ptI - ptJ).lengthSq();
        if (distCSq < sqDistCutoff) {
          termPoints[j].reset();
        }
      }
    }
  }
  std::erase_if(termPoints, [](const auto &p) -> bool { return !p; });
}

}  // namespace

SubshapeStarts::SubshapeStarts(const ROMol &ref, const ROMol &fit,
                               int refConfId, int fitConfId,
                               const SubshapeOptions &options) {
  d_refGrid = buildGrid(ref, refConfId, options);
  d_fitGrid = buildGrid(fit, fitConfId, options);

  buildTerminalPoints(ref, refConfId, options, *d_refGrid, d_refPoints);
  std::cout << "Ref points : " << d_refPoints.size() << " " << std::endl;
  for (const auto &p : d_refPoints) {
    std::cout << p->getPosition().x << "," << p->getPosition().y << ","
              << p->getPosition().z << ";";
  }
  std::cout << std::endl;
  buildTerminalPoints(fit, fitConfId, options, *d_fitGrid, d_fitPoints);
  std::cout << "Fit points : " << d_fitPoints.size() << " " << std::endl;
  for (const auto &p : d_fitPoints) {
    std::cout << p->getPosition().x << "," << p->getPosition().y << ","
              << p->getPosition().z << ";";
  }
  std::cout << std::endl;
  addSkeletonPoints(*d_fitGrid, options, d_fitPoints);
  std::cout << "Fit points + skeleton points : " << d_fitPoints.size() << " "
            << std::endl;
  for (const auto &p : d_fitPoints) {
    std::cout << p->getPosition().x << "," << p->getPosition().y << ","
              << p->getPosition().z << ";";
  }
  std::cout << std::endl;
  double minDistSq = std::numeric_limits<double>::max();
  size_t minI, minJ;
  for (size_t i = 1; i < d_fitPoints.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      double distSq =
          (d_fitPoints[i]->getPosition() - d_fitPoints[j]->getPosition())
              .lengthSq();
      if (distSq < minDistSq) {
        minDistSq = distSq;
        minI = i;
        minJ = j;
      }
      if (distSq < 1.0) {
        std::cout << sqrt(distSq) << " for " << i << " -> " << j << std::endl;
      }
    }
  }
  std::cout << "Closest : " << sqrt(minDistSq) << " for " << minI << " -> "
            << minJ << std::endl;
}

namespace {
void copyPoints(const std::vector<std::unique_ptr<SubshapePoint>> &inPoints,
                std::vector<std::unique_ptr<SubshapePoint>> &outPoints) {
  outPoints.clear();
  outPoints.reserve(inPoints.size());
  for (const auto &p : inPoints) {
    outPoints.emplace_back(new SubshapePoint(*p));
  }
}
}  // namespace

SubshapeStarts::SubshapeStarts(const SubshapeStarts &other) {
  d_refGrid = std::make_unique<RDGeom::UniformGrid3D>(*other.d_refGrid);
  d_fitGrid = std::make_unique<RDGeom::UniformGrid3D>(*other.d_fitGrid);
  copyPoints(other.d_refPoints, d_refPoints);
  copyPoints(other.d_fitPoints, d_fitPoints);
}

SubshapeStarts &SubshapeStarts::operator=(const SubshapeStarts &other) {
  if (this == &other) {
    return *this;
  }
  d_refGrid = std::make_unique<RDGeom::UniformGrid3D>(*other.d_refGrid);
  d_fitGrid = std::make_unique<RDGeom::UniformGrid3D>(*other.d_fitGrid);
  copyPoints(other.d_refPoints, d_refPoints);
  copyPoints(other.d_fitPoints, d_fitPoints);
  return *this;
}
}  // namespace GaussianShape
}  // namespace RDKit
