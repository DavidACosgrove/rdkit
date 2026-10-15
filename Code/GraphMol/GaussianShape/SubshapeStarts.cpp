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

#include <Eigen/Geometry>

#include "Geometry/Transform3D.h"
#include "Numerics/Alignment/AlignPoints.h"
#include "GraphMol/GaussianShape/ShapeInput.h"
#include "GraphMol/GaussianShape/SubshapeStarts.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"
#include "GraphMol/MolAlign/AlignMolecules.h"
#include "GraphMol/MolTransforms/MolTransforms.h"

#include <SimDivPickers/LeaderPicker.h>

namespace RDKit {
namespace GaussianShape {

namespace {

void makeInitialPoints(
    const ShapeInput &shape, const SubshapeOptions &options,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  auto atomPoints = shape.getAtomPoints();
  double rWSq = options.r_w * options.r_w;
  std::vector<std::vector<unsigned int>> nbours(atomPoints.size());

  for (size_t i = 1; i < atomPoints.size(); ++i) {
    for (size_t j = 1; j < i; ++j) {
      double sqDist = (atomPoints[i] - atomPoints[j]).lengthSq();
      if (sqDist < rWSq) {
        nbours[i].push_back(j);
        nbours[j].push_back(i);
      }
    }
  }

  // The paper only does this for atoms with a small number of neighbours
  // by distance.  We're doing it for all as we are not just interested
  // in terminal atoms.
  for (size_t i = 0; i < nbours.size(); ++i) {
    // This isn't as described in the paper, because that didn't make any
    // sense, and indeed the algorithm in Figure 6 is different from
    // equation 3.  This follows what Greg did in his Python implementation
    // in $RDBASE/rdkit/Chem/Subshape/BuilderUtils.py
    auto pnt = atomPoints[i];
    double totWt = 1.0;  // For atom i
    for (const auto n : nbours[i]) {
      // Add one to the weights because Greg put i in the nbours lists, which
      // I haven't done.
      double wt =
          static_cast<double>(nbours[i].size() + 1) / (nbours[n].size() + 1);
      pnt += atomPoints[n] * wt;
      totWt += wt;
    }
    pnt /= totWt;
    termPoints.push_back(
        std::make_unique<SubshapePoint>(pnt, nbours[i].size() + 1));
  }
}

void prunePoints(std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  class DistFunctor {
   public:
    explicit DistFunctor(
        const std::vector<std::unique_ptr<SubshapePoint>> &termPoints)
        : d_points(&termPoints) {}
    ~DistFunctor() = default;
    double operator()(const unsigned int i, const unsigned int j) {
      return ((*d_points)[i]->getPosition() - (*d_points)[j]->getPosition())
          .lengthSq();
    }
    const std::vector<std::unique_ptr<SubshapePoint>> *d_points;
  };
  const RDPickers::LeaderPicker leaderPicker;
  DistFunctor distFunctor(termPoints);
  // The LeaderPicker starts with the first entry for the things to be
  // picked.  Make that the point with the fewest neighbours (the "most
  // terminal point" you might say).
  std::ranges::sort(termPoints, [](const auto &pt1, const auto &pt2) -> bool {
    return pt1->getNumNbours() < pt2->getNumNbours();
  });
  auto picks = leaderPicker.lazyPick(distFunctor, termPoints.size(), 0, 1.0);
  // Allow for mysterious LeaderPicker behaviour where it returns a
  // vector of 0s when it should only pick 1 shape.
  std::ranges::sort(picks);
  auto [first, last] = std::ranges::unique(picks);
  picks.erase(first, last);
  std::vector<std::unique_ptr<SubshapePoint>> tmpPoints;
  for (const auto &pick : picks) {
    tmpPoints.emplace_back(std::move(termPoints[pick]));
  }
  termPoints = std::move(tmpPoints);
}

}  // namespace

SubshapeStarts::SubshapeStarts(const ShapeInput &ref, const ShapeInput &fit,
                               const SubshapeOptions &options) {
  std::cout << "RRRRRRRRRRRRRef" << std::endl;
  makeInitialPoints(ref, options, d_refPoints);
  prunePoints(d_refPoints);
  std::cout << "Ref points : " << d_refPoints.size() << " " << std::endl;
  for (const auto &p : d_refPoints) {
    std::cout << p->getPosition().x << "," << p->getPosition().y << ","
              << p->getPosition().z << ";";
  }
  std::cout << std::endl;
  std::cout << "FFFFFFFFFFFit" << std::endl;
  makeInitialPoints(fit, options, d_fitPoints);
  prunePoints(d_fitPoints);
  std::cout << "Orig Fit points : " << d_fitPoints.size() << " " << std::endl;
  for (const auto &p : d_fitPoints) {
    std::cout << p->getPosition().x << "," << p->getPosition().y << ","
              << p->getPosition().z << ";";
  }
  std::cout << std::endl;

  buildDistMatrices();
  std::array<double, 3> refTriDists;
  std::array<double, 3> fitTriDists;

  for (size_t ref_i = 2; ref_i < d_refPoints.size(); ++ref_i) {
    for (size_t ref_j = 1; ref_j < ref_i; ++ref_j) {
      for (size_t ref_k = 0; ref_k < ref_j; ++ref_k) {
        refTriDists[0] = d_refPointsDists[ref_i][ref_j];
        refTriDists[1] = d_refPointsDists[ref_i][ref_k];
        refTriDists[2] = d_refPointsDists[ref_j][ref_k];
        for (size_t tgt_i = 2; tgt_i < d_fitPoints.size(); ++tgt_i) {
          for (size_t tgt_j = 1; tgt_j < tgt_i; ++tgt_j) {
            for (size_t tgt_k = 0; tgt_k < tgt_j; ++tgt_k) {
              fitTriDists[0] = d_fitPointsDists[tgt_i][tgt_j];
              fitTriDists[1] = d_fitPointsDists[tgt_i][tgt_k];
              fitTriDists[2] = d_fitPointsDists[tgt_j][tgt_k];
              if (fabs(refTriDists[0] - fitTriDists[0]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[1]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[2]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_i, ref_j, ref_k, tgt_i, tgt_j, tgt_k});
              }
              if (fabs(refTriDists[0] - fitTriDists[1]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[0]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[2]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_i, ref_k, ref_j, tgt_i, tgt_j, tgt_k});
              }
              if (fabs(refTriDists[0] - fitTriDists[0]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[2]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[1]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_j, ref_i, ref_k, tgt_i, tgt_j, tgt_k});
              }
              if (fabs(refTriDists[0] - fitTriDists[2]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[0]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[1]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_j, ref_k, ref_i, tgt_i, tgt_j, tgt_k});
              }
              if (fabs(refTriDists[0] - fitTriDists[1]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[2]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[0]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_k, ref_i, ref_j, tgt_i, tgt_j, tgt_k});
              }
              if (fabs(refTriDists[0] - fitTriDists[2]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[1] - fitTriDists[1]) >
                      options.triangleSideTolerance ||
                  fabs(refTriDists[2] - fitTriDists[0]) >
                      options.triangleSideTolerance) {
                d_possStarts.emplace_back(std::array<std::uint64_t, 6>{
                    ref_k, ref_j, ref_i, tgt_i, tgt_j, tgt_k});
              }
            }
          }
        }
      }
    }
  }
  std::cout << "Total starts :: " << d_possStarts.size() << std::endl;
  if (options.clusterStartsThreshold > 0.0) {
    clusterStartPoints(fit, options);
  }
  std::cout << "Clustered starts : " << d_possStarts.size() << std::endl;
  sortStartPoints();
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
  copyPoints(other.d_refPoints, d_refPoints);
  copyPoints(other.d_fitPoints, d_fitPoints);
}

SubshapeStarts &SubshapeStarts::operator=(const SubshapeStarts &other) {
  if (this == &other) {
    return *this;
  }
  copyPoints(other.d_refPoints, d_refPoints);
  copyPoints(other.d_fitPoints, d_fitPoints);
  return *this;
}

std::unique_ptr<RDGeom::Transform3D> SubshapeStarts::getNextStartTransform(
    double *ssd) {
  if (d_nextStart >= d_possStarts.size()) {
    return std::unique_ptr<RDGeom::Transform3D>();
  }
  auto trans = std::make_unique<RDGeom::Transform3D>();
  std::vector<const RDGeom::Point3D *> refs(3);
  refs[0] = &d_refPoints[d_possStarts[d_nextStart][0]]->getPosition();
  refs[1] = &d_refPoints[d_possStarts[d_nextStart][1]]->getPosition();
  refs[2] = &d_refPoints[d_possStarts[d_nextStart][2]]->getPosition();
  std::vector<const RDGeom::Point3D *> fits(3);
  fits[0] = &d_fitPoints[d_possStarts[d_nextStart][3]]->getPosition();
  fits[1] = &d_fitPoints[d_possStarts[d_nextStart][4]]->getPosition();
  fits[2] = &d_fitPoints[d_possStarts[d_nextStart][5]]->getPosition();
  auto rssd = RDNumeric::Alignments::AlignPoints(refs, fits, *trans);
  ++d_nextStart;
  if (ssd) {
    *ssd = rssd;
  }
  return trans;
}

std::unique_ptr<std::array<double, 7>> SubshapeStarts::getNextStartQuaternion(
    double *ssd) {
  auto trans = getNextStartTransform(ssd);
  std::unique_ptr<std::array<double, 7>> retVal;
  if (!trans) {
    return retVal;
  }

  Eigen::Matrix3d r;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r(i, j) = trans->getValUnchecked(i, j);
    }
  }
  auto q = Eigen::Quaterniond(r);
  retVal.reset(new std::array<double, 7>());
  (*retVal)[0] = q.x();
  (*retVal)[1] = q.y();
  (*retVal)[2] = q.z();
  (*retVal)[3] = q.w();
  (*retVal)[4] = trans->getValUnchecked(0, 3);
  (*retVal)[5] = trans->getValUnchecked(1, 3);
  (*retVal)[6] = trans->getValUnchecked(2, 3);

  return retVal;
}

namespace {
std::vector<std::vector<double>> buildDistMatrix(
    const std::vector<std::unique_ptr<SubshapePoint>> &points) {
  std::vector<std::vector<double>> dists(points.size(),
                                         std::vector<double>(points.size()));
  for (size_t i = 1; i < points.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      dists[i][j] = dists[j][i] =
          (points[i]->getPosition() - points[j]->getPosition()).length();
    }
  }
  return dists;
}
}  // namespace

void SubshapeStarts::buildDistMatrices() {
  d_refPointsDists = buildDistMatrix(d_refPoints);
  d_fitPointsDists = buildDistMatrix(d_fitPoints);
}

void SubshapeStarts::clusterStartPoints(const ShapeInput &fit,
                                        const SubshapeOptions &options) {
  std::vector<size_t> selTrans;
  std::vector<std::vector<RDGeom::Point3D>> selShapes;
  selTrans.reserve(getNumPossibleStarts());
  double ssdCutoff =
      options.rmsdTriangleTolerance * options.rmsdTriangleTolerance * 3.0;
  double shapeCutoff =
      options.clusterStartsThreshold * options.clusterStartsThreshold;

  for (size_t i = 0; i < getNumPossibleStarts(); ++i) {
    double ssd = 0.0;
    const auto trans = getNextStartTransform(&ssd);
    // reject any where the triangle points end up too far apart.
    if (ssd > ssdCutoff) {
      continue;
    }
    auto fitPoints = fit.getAtomPoints();
    std::for_each(fitPoints.begin(), fitPoints.end(),
                  [&trans](RDGeom::Point3D &p) { trans->TransformPoint(p); });

    if (selTrans.empty()) {
      selTrans.push_back(i);
      selShapes.emplace_back(std::move(fitPoints));
    } else {
      bool foundSimilar = false;
      for (const auto &sp : selShapes) {
        double msd = 0.0;
        for (size_t j = 0; j < sp.size(); ++j) {
          msd += (sp[j] - fitPoints[j]).lengthSq();
        }
        msd /= sp.size();
        if (msd < shapeCutoff) {
          foundSimilar = true;
          break;
        }
      }
      if (!foundSimilar) {
        selTrans.push_back(i);
        selShapes.emplace_back(std::move(fitPoints));
      }
    }
  }
  std::vector<std::array<std::uint64_t, 6>> tmpStarts;
  tmpStarts.reserve(selTrans.size());
  for (size_t i = 0; i < selTrans.size(); ++i) {
    tmpStarts.push_back(d_possStarts[selTrans[i]]);
  }
  d_nextStart = 0;
  d_possStarts = std::move(tmpStarts);
}

void SubshapeStarts::sortStartPoints() {
  std::vector<std::pair<size_t, double>> ssd_indices;
  ssd_indices.reserve(getNumPossibleStarts());
  for (size_t i = 0; i < getNumPossibleStarts(); ++i) {
    double ssd = 0.0;
    const auto trans = getNextStartTransform(&ssd);
    ssd_indices.push_back(std::make_pair(i, ssd));
  }
  std::ranges::sort(ssd_indices, [](const auto &p1, const auto &p2) -> bool {
    return p1.second < p2.second;
  });
  std::vector<std::array<std::uint64_t, 6>> tmpStarts;
  tmpStarts.reserve(getNumPossibleStarts());
  for (const auto &c : ssd_indices) {
    tmpStarts.push_back(std::move(d_possStarts[c.first]));
  }
  d_nextStart = 0;
  d_possStarts = std::move(tmpStarts);
}

}  // namespace GaussianShape
}  // namespace RDKit
