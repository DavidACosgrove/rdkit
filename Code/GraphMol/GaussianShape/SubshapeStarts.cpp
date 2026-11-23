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

#include <Eigen/Geometry>

#include "Geometry/Transform3D.h"
#include "Numerics/Alignment/AlignPoints.h"
#include "GraphMol/GaussianShape/ShapeInput.h"
#include "GraphMol/GaussianShape/SubshapeStarts.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"
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
    for (size_t j = 0; j < i; ++j) {
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

void makeFeaturePoints(
    const ShapeInput &shape,
    std::vector<std::unique_ptr<SubshapePoint>> &termPoints) {
  const auto &coords = shape.getCoords();
  const auto &types = shape.getFeatureTypes();
  for (unsigned int i = 0, j = shape.getNumAtoms(), k = 3 * shape.getNumAtoms();
       i < shape.getNumFeatures(); ++i, ++j, k += 3) {
    termPoints.push_back(std::make_unique<SubshapePoint>(
        RDGeom::Point3D(coords[k], coords[k + 1], coords[k + 2]), types[j], 0));
  }
}

}  // namespace

SubshapeStarts::SubshapeStarts(const ShapeInput &ref, const ShapeInput &fit,
                               const SubshapeOptions &options) {
  PRECONDITION(!ref.getCoords().empty() && !fit.getCoords().empty(),
               "Both ShapeInput objects must have at least 1 point.")
  if (options.byPharmacophoreFeature) {
    makeFeaturePoints(ref, d_refPoints);
  } else {
    makeInitialPoints(ref, options, d_refPoints);
    prunePoints(d_refPoints);
  }
#if 0
  std::cout << "Ref shape : " << MolToCXSmiles(*ref.shapeToMol(true))
            << std::endl;
  std::cout << "Ref types : ";
  for (unsigned int i = 0; i < ref.getNumFeatures(); ++i) {
    std::cout << ref.getFeatureTypes()[i + ref.getNumAtoms()] << " ";
  }
  std::cout << std::endl;
  std::cout << "Ref points : ";
  for (unsigned int i = 0; i < d_refPoints.size(); i++) {
    std::cout << "[Cl]";
    if (i < d_refPoints.size() - 1) {
      std::cout << ".";
    } else {
      std::cout << " |(";
    }
  }
  for (const auto &p : d_refPoints) {
    std::cout << p->getPosition()[0] << "," << p->getPosition()[1] << ","
              << p->getPosition()[2] << ";";
  }
  std::cout << ")|" << std::endl;
#endif
  if (options.byPharmacophoreFeature) {
    makeFeaturePoints(fit, d_fitPoints);
  } else {
    makeInitialPoints(fit, options, d_fitPoints);
    prunePoints(d_fitPoints);
  }
#if 0
  std::cout << "Fit shape : " << MolToCXSmiles(*fit.shapeToMol(true))
            << std::endl;
  for (unsigned int i = 0; i < fit.getNumFeatures(); ++i) {
    std::cout << fit.getFeatureTypes()[i + fit.getNumAtoms()] << " ";
  }
  std::cout << std::endl;
  std::cout << "Fit points : ";
  for (unsigned int i = 0; i < d_fitPoints.size(); i++) {
    std::cout << "[Cl]";
    if (i < d_fitPoints.size() - 1) {
      std::cout << ".";
    } else {
      std::cout << " |(";
    }
  }
  for (const auto &p : d_fitPoints) {
    std::cout << p->getPosition()[0] << "," << p->getPosition()[1] << ","
              << p->getPosition()[2] << ";";
  }
  std::cout << ")|" << std::endl;
#endif
  buildDistMatrices();
  buildStarts(options);
  // std::cout << "InitialNumber of starts : " << getNumPossibleStarts()
  // << std::endl;

  if (getNumPossibleStarts() > 1 && options.clusterStartsThreshold > 0.0) {
    clusterStartPoints(fit, options);
  }
  sortStartPoints();
  // std::cout << "Number of starts : " << getNumPossibleStarts() << std::endl;
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
  auto xform = getTransform(d_nextStart, ssd);
  ++d_nextStart;
  return xform;
}

std::unique_ptr<std::array<double, 7>> SubshapeStarts::getNextStartQuaternion(
    double *ssd) {
  auto xform = getNextStartTransform(ssd);
  if (!xform) {
    return std::unique_ptr<std::array<double, 7>>();
  }

  // Le Chat provided this.  Apparently Eigen and the RDKit have different
  // conventions.
  std::unique_ptr<std::array<double, 7>> retVal;
  Eigen::Matrix3d m;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m(i, j) = xform->getValUnchecked(i, j);
    }
  }
  auto q = Eigen::Quaterniond(m);
  q.normalize();
  retVal.reset(new std::array<double, 7>());
  (*retVal)[0] = q.w();
  (*retVal)[1] = -q.x();
  (*retVal)[2] = -q.y();
  (*retVal)[3] = -q.z();
  (*retVal)[4] = xform->getValUnchecked(0, 3);
  (*retVal)[5] = xform->getValUnchecked(1, 3);
  (*retVal)[6] = xform->getValUnchecked(2, 3);
  return retVal;
}

std::unique_ptr<RDGeom::Transform3D> SubshapeStarts::getTransform(
    std::uint64_t num, double *ssd) {
  if (num >= d_possStarts.size()) {
    return std::unique_ptr<RDGeom::Transform3D>();
  }
  std::unique_ptr<RDGeom::Transform3D> trans =
      std::make_unique<RDGeom::Transform3D>(*d_possStarts[num].first);
  if (ssd) {
    *ssd = d_possStarts[num].second;
  }
  return trans;
}

std::unique_ptr<RDGeom::Transform3D> SubshapeStarts::calcTransform(
    size_t ref_i, size_t ref_j, size_t ref_k, size_t fit_i, size_t fit_j,
    size_t fit_k, double &ssd) {
  auto trans = std::make_unique<RDGeom::Transform3D>();
  std::vector<const RDGeom::Point3D *> refs(3);
  refs[0] = &d_refPoints[ref_i]->getPosition();
  refs[1] = &d_refPoints[ref_j]->getPosition();
  refs[2] = &d_refPoints[ref_k]->getPosition();
  std::vector<const RDGeom::Point3D *> fits(3);
  fits[0] = &d_fitPoints[fit_i]->getPosition();
  fits[1] = &d_fitPoints[fit_j]->getPosition();
  fits[2] = &d_fitPoints[fit_k]->getPosition();
  ssd = RDNumeric::Alignments::AlignPoints(refs, fits, *trans);
  return trans;
}

std::unique_ptr<RDGeom::Transform3D> SubshapeStarts::calcTransform(
    size_t ref_i, size_t ref_j, size_t fit_i, size_t fit_j, double &ssd) {
  auto trans = std::make_unique<RDGeom::Transform3D>();
  std::vector<const RDGeom::Point3D *> refs(2);
  refs[0] = &d_refPoints[ref_i]->getPosition();
  refs[1] = &d_refPoints[ref_j]->getPosition();
  std::vector<const RDGeom::Point3D *> fits(2);
  fits[0] = &d_fitPoints[fit_i]->getPosition();
  fits[1] = &d_fitPoints[fit_j]->getPosition();
  ssd = RDNumeric::Alignments::AlignPoints(refs, fits, *trans);
  return trans;
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

void SubshapeStarts::buildStarts(const SubshapeOptions &options) {
  if (d_refPoints.size() == 1 || d_fitPoints.size() == 1) {
    // If one of the shapes only has 1 point, we can only do point-onto-point
    // and rotate starts.
    buildPointStarts(options);
  } else if (d_refPoints.size() == 2 || d_fitPoints.size() == 2) {
    // If one of the shapes only has 2 points, we can only do line-onto-line
    // and rotate starts.
    buildLineStarts(options);
  } else {
    // Both shapes have at least 3 points,
    buildTriangleStarts(options);
  }
}
void SubshapeStarts::buildPointStarts(const SubshapeOptions &options) {
  const static std::vector<std::array<double, 4>> quats{{1.0, 0.0, 0.0, 0.0},
                                                        {0.0, 1.0, 0.0, 0.0},
                                                        {0.0, 0.0, 1.0, 0.0},
                                                        {0.0, 0.0, 0.0, 1.0}};
  static std::vector<RDGeom::Transform3D> rots;
  if (rots.empty()) {
    for (const auto &q : quats) {
      rots.push_back(RDGeom::Transform3D());
      rots.back().SetRotationFromQuaternion(q.data());
    }
  }

  for (size_t ref_i = 0; ref_i < d_refPoints.size(); ++ref_i) {
    RDGeom::Transform3D fromOrigin;
    fromOrigin.SetTranslation(d_refPoints[ref_i]->getPosition());
    for (size_t tgt_i = 0; tgt_i < d_fitPoints.size(); ++tgt_i) {
      // Move each target point onto each reference point and do
      // a ROTATE_180-style initialisation.
      RDGeom::Transform3D toOrigin;
      toOrigin.SetTranslation(d_fitPoints[tgt_i]->getPosition());
      for (size_t i = 0; i < 4; ++i) {
        RDGeom::Transform3D xform = fromOrigin * rots[i] * toOrigin;
        d_possStarts.emplace_back(
            std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                std::make_unique<RDGeom::Transform3D>(xform), 0.0));
      }
    }
  }
}

void SubshapeStarts::buildLineStarts(const SubshapeOptions &options) {
  double ssd;
  // Rotations of 90, 180, 270 degrees.
  static const std::array<double, 4> sinT{0.0, sin(std::numbers::pi / 2.0),
                                          sin(std::numbers::pi),
                                          sin(std::numbers::pi * 1.5)};
  static const std::array<double, 4> cosT{1.0, cos(std::numbers::pi / 2.0),
                                          cos(std::numbers::pi),
                                          cos(std::numbers::pi * 1.5)};

  for (size_t ref_i = 1; ref_i < d_refPoints.size(); ++ref_i) {
    for (size_t ref_j = 0; ref_j < ref_i; ++ref_j) {
      RDGeom::Point3D rotAxis =
          d_refPoints[ref_i]->getPosition().directionVector(
              d_refPoints[ref_j]->getPosition());
      RDGeom::Point3D midRef = (d_refPoints[ref_i]->getPosition() +
                                d_refPoints[ref_j]->getPosition()) /
                               2.0;
      RDGeom::Transform3D toOrigin, fromOrigin;
      toOrigin.SetTranslation(
          RDGeom::Point3D(-midRef[0], -midRef[1], -midRef[2]));
      fromOrigin.SetTranslation(
          RDGeom::Point3D(midRef[0], midRef[1], midRef[2]));

      for (size_t tgt_i = 1; tgt_i < d_fitPoints.size(); ++tgt_i) {
        for (size_t tgt_j = 0; tgt_j < tgt_i; ++tgt_j) {
          if (fabs(d_refPointsDists[ref_i][ref_j] -
                   d_fitPointsDists[tgt_i][tgt_j]) <
              options.triangleSideTolerance) {
            auto trans1 = calcTransform(ref_i, ref_j, tgt_i, tgt_j, ssd);
            for (int k = 0; k < 4; ++k) {
              RDGeom::Transform3D rot;
              rot.SetRotation(sinT[k], cosT[k], rotAxis);
              RDGeom::Transform3D fullXform =
                  fromOrigin * rot * toOrigin * (*trans1);
              d_possStarts.emplace_back(
                  std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                      std::make_unique<RDGeom::Transform3D>(fullXform), ssd));
            }
            auto trans2 = calcTransform(ref_j, ref_i, tgt_i, tgt_j, ssd);
            for (int k = 0; k < 4; ++k) {
              RDGeom::Transform3D rot;
              rot.SetRotation(sinT[k], cosT[k], rotAxis);
              RDGeom::Transform3D fullXform =
                  fromOrigin * rot * toOrigin * (*trans2);
              d_possStarts.emplace_back(
                  std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                      std::make_unique<RDGeom::Transform3D>(fullXform), ssd));
            }
          }
        }
      }
    }
  }
}

namespace {
bool checkTriangleMatch(const std::array<double, 3> &refTriDists,
                        const std::array<double, 3> &fitTriDists,
                        const std::array<int, 3> &refTypes,
                        std::array<int, 3> &fitTypes, const unsigned int i,
                        const unsigned int j, const unsigned int k,
                        const double triangleSideTolerance) {
  if (refTypes[0] != fitTypes[i] || refTypes[1] != fitTypes[j] ||
      refTypes[2] != fitTypes[k]) {
    return false;
  }
  if (fabs(refTriDists[0] - fitTriDists[i]) >= triangleSideTolerance ||
      fabs(refTriDists[1] - fitTriDists[j]) >= triangleSideTolerance ||
      fabs(refTriDists[2] - fitTriDists[k]) >= triangleSideTolerance) {
    return true;
  }
  return true;
}
}  // namespace
void SubshapeStarts::buildTriangleStarts(const SubshapeOptions &options) {
  std::array<double, 3> refTriDists;
  std::array<double, 3> fitTriDists;
  std::array<int, 3> refTypes;
  std::array<int, 3> fitTypes;
  double ssd;

  for (size_t ref_i = 2; ref_i < d_refPoints.size(); ++ref_i) {
    refTypes[0] = d_refPoints[ref_i]->getType();
    for (size_t ref_j = 1; ref_j < ref_i; ++ref_j) {
      refTypes[1] = d_refPoints[ref_j]->getType();
      for (size_t ref_k = 0; ref_k < ref_j; ++ref_k) {
        refTypes[2] = d_refPoints[ref_k]->getType();
        refTriDists[0] = d_refPointsDists[ref_i][ref_j];
        refTriDists[1] = d_refPointsDists[ref_i][ref_k];
        refTriDists[2] = d_refPointsDists[ref_j][ref_k];
        for (size_t tgt_i = 2; tgt_i < d_fitPoints.size(); ++tgt_i) {
          fitTypes[0] = d_fitPoints[tgt_i]->getType();
          for (size_t tgt_j = 1; tgt_j < tgt_i; ++tgt_j) {
            fitTypes[1] = d_fitPoints[tgt_j]->getType();
            for (size_t tgt_k = 0; tgt_k < tgt_j; ++tgt_k) {
              fitTypes[2] = d_fitPoints[tgt_k]->getType();
              fitTriDists[0] = d_fitPointsDists[tgt_i][tgt_j];
              fitTriDists[1] = d_fitPointsDists[tgt_i][tgt_k];
              fitTriDists[2] = d_fitPointsDists[tgt_j][tgt_k];
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 0, 1, 2,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_i, ref_j, ref_k, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 1, 0, 2,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_i, ref_k, ref_j, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 0, 2, 1,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_j, ref_i, ref_k, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 2, 0, 1,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_j, ref_k, ref_i, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 1, 2, 0,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_k, ref_i, ref_j, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
              if (checkTriangleMatch(refTriDists, fitTriDists, fitTypes,
                                     refTypes, 2, 1, 0,
                                     options.triangleSideTolerance)) {
                auto trans = calcTransform(ref_k, ref_j, ref_i, tgt_i, tgt_j,
                                           tgt_k, ssd);
                d_possStarts.emplace_back(
                    std::pair<std::unique_ptr<RDGeom::Transform3D>, double>(
                        std::move(trans), ssd));
              }
            }
          }
        }
      }
    }
  }
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

  for (size_t i = 0; i < d_possStarts.size(); ++i) {
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
  std::vector<std::pair<std::unique_ptr<RDGeom::Transform3D>, double>>
      tmpStarts;
  tmpStarts.reserve(selTrans.size());
  for (size_t i = 0; i < selTrans.size(); ++i) {
    tmpStarts.push_back(std::move(d_possStarts[selTrans[i]]));
  }
  d_nextStart = 0;
  d_possStarts = std::move(tmpStarts);
}

void SubshapeStarts::sortStartPoints() {
  std::vector<std::pair<size_t, double>> ssd_indices;
  ssd_indices.reserve(d_possStarts.size());
  for (size_t i = 0; i < d_possStarts.size(); ++i) {
    double ssd = 0.0;
    const auto trans = getNextStartTransform(&ssd);
    ssd_indices.push_back(std::make_pair(i, ssd));
  }
  std::ranges::sort(ssd_indices, [](const auto &p1, const auto &p2) -> bool {
    return p1.second < p2.second;
  });
  std::vector<std::pair<std::unique_ptr<RDGeom::Transform3D>, double>>
      tmpStarts;
  tmpStarts.reserve(getNumPossibleStarts());
  for (const auto &c : ssd_indices) {
    tmpStarts.push_back(std::move(d_possStarts[c.first]));
  }
  d_nextStart = 0;
  d_possStarts = std::move(tmpStarts);
}

}  // namespace GaussianShape
}  // namespace RDKit
