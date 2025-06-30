//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SearchShapeInput.h>
#include <SimDivPickers/LeaderPicker.h>

SearchShapeInput::SearchShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  std::stringstream ss(str);
  boost::archive::text_iarchive ia(ss);
  ia &*this;
#endif
}

SearchShapeInput::SearchShapeInput(const ShapeInput &other)
    : ShapeInput(other) {
  confCoords = std::vector<std::vector<float>>(1, coord);
  // We don't actually know what conformer it came from
  molConfs = std::vector<unsigned int>(1, 0);
  // Keep the dummyVols in synch with the other data, but
  // clearly flag it as not calculated.
  dummyVols = std::vector<double>(1, -1);
  sovs = std::vector<double>(1, sov);
  sofs = std::vector<double>(1, sof);
  shifts = std::vector<std::vector<double>>(1, shift);
}

void SearchShapeInput::merge(SearchShapeInput &other) {
  if (!confCoords.empty() &&
      confCoords.front().size() != other.confCoords.front().size()) {
    BOOST_LOG(rdWarningLog) << "Can't merge shapes as different sizes.\n";
    return;
  }
  if (confCoords.empty()) {
    coord = other.coord;
    alpha_vector = other.alpha_vector;
    atom_type_vector = other.atom_type_vector;
    volumeAtomIndexVector = other.volumeAtomIndexVector;
    colorAtomType2IndexVectorMap = other.colorAtomType2IndexVectorMap;
    shift = other.shift;
    sov = other.sov;
    sof = other.sof;
    numDummies = other.numDummies;
    dummyVol = other.dummyVol;
    actConf = other.actConf;
    confCoords = std::move(other.confCoords);
    molConfs = std::move(other.molConfs);
    dummyVols = std::move(other.dummyVols);
    sovs = std::move(other.sovs);
    sofs = std::move(other.sofs);
    shifts = std::move(other.shifts);
  } else {
    confCoords.reserve(confCoords.size() + other.confCoords.size());
    confCoords.insert(confCoords.end(),
                      std::make_move_iterator(other.confCoords.begin()),
                      std::make_move_iterator(other.confCoords.end()));
    molConfs.reserve(molConfs.size() + other.molConfs.size());
    molConfs.insert(molConfs.end(),
                    std::make_move_iterator(other.molConfs.begin()),
                    std::make_move_iterator(other.molConfs.end()));
    dummyVols.reserve(dummyVols.size() + other.dummyVols.size());
    dummyVols.insert(dummyVols.end(),
                     std::make_move_iterator(other.dummyVols.begin()),
                     std::make_move_iterator(other.dummyVols.end()));
    sovs.reserve(sovs.size() + other.sovs.size());
    sovs.insert(sovs.end(), std::make_move_iterator(other.sovs.begin()),
                std::make_move_iterator(other.sovs.end()));
    sofs.reserve(sofs.size() + other.sofs.size());
    sofs.insert(sofs.end(), std::make_move_iterator(other.sofs.begin()),
                std::make_move_iterator(other.sofs.end()));
    shifts.reserve(shifts.size() + other.shifts.size());
    shifts.insert(shifts.end(), std::make_move_iterator(other.shifts.begin()),
                  std::make_move_iterator(other.shifts.end()));
  }
  other.confCoords.clear();
  other.molConfs.clear();
  other.dummyVols.clear();
  other.sovs.clear();
  other.sofs.clear();
  other.shifts.clear();
}

ShapeInput SearchShapeInput::makeSingleShape(unsigned int shapeNum) const {
  if (shapeNum >= confCoords.size()) {
    shapeNum = 0;
  }
  ShapeInput shape;
  shape.coord = confCoords[shapeNum];
  shape.alpha_vector = alpha_vector;
  shape.atom_type_vector = atom_type_vector;
  shape.volumeAtomIndexVector = volumeAtomIndexVector;
  shape.colorAtomType2IndexVectorMap = colorAtomType2IndexVectorMap;
  shape.shift = shifts[shapeNum];
  shape.sov = sov;
  shape.sof = sof;
  return shape;
}

void SearchShapeInput::setActiveShape(unsigned int shapeNum) {
  PRECONDITION(shapeNum < confCoords.size(), "confNum is out of bounds");
  actConf = shapeNum;
  coord = confCoords[shapeNum];
  dummyVol = dummyVols[shapeNum];
  sov = sovs[shapeNum];
  sof = sofs[shapeNum];
  shift = shifts[shapeNum];
}

std::string SearchShapeInput::toString() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  std::stringstream ss;
  boost::archive::text_oarchive oa(ss);
  oa &*this;
  return ss.str();
#endif
}

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void SearchShapeInput::serialize(Archive &ar, const unsigned int) {
  ar &boost::serialization::base_object<ShapeInput>(*this);
  ar & numDummies;
  ar & dummyVol;
  ar & actConf;
  ar & confCoords;
  ar & molConfs;
  ar & dummyVols;
  ar & sovs;
  ar & sofs;
  ar & shifts;
}
#endif

void pruneShapes(SearchShapeInput &shapes, double simThreshold) {
  if (shapes.confCoords.size() < 2) {
    return;
  }
  class DistFunctor {
   public:
    DistFunctor(const SearchShapeInput &shapes) : d_shapes(shapes) {
      d_shapei = d_shapes.makeSingleShape(0);
      d_shapej = d_shapes.makeSingleShape(1);
      d_matrix = std::vector<float>(12, 0.0);
    }
    ~DistFunctor() = default;
    double operator()(unsigned int i, unsigned int j) {
      d_shapei.coord = d_shapes.confCoords[i];
      d_shapej.coord = d_shapes.confCoords[j];
      auto [st, ct] = AlignShape(d_shapei, d_shapej, d_matrix);
      double res = 2.0 - (st + ct);
      return res;
    }
    const SearchShapeInput &d_shapes;
    ShapeInput d_shapei, d_shapej;
    std::vector<float> d_matrix;
  };
  RDPickers::LeaderPicker leaderPicker;
  DistFunctor distFunctor(shapes);
  auto picks = leaderPicker.lazyPick(distFunctor, shapes.confCoords.size(), 0,
                                     2.0 - simThreshold);
  // Allow for mysterious LeaderPicker behaviour where it returns a full vector
  // of 0s when it should only pick 1 shape.
  if (picks.size() == shapes.confCoords.size()) {
    std::ranges::sort(picks);
    auto [first, last] = std::ranges::unique(picks);
    picks.erase(first, last);
  }
  std::vector<std::vector<float>> newCoords;
  newCoords.reserve(picks.size());
  std::vector<unsigned int> newMolConfs;
  newMolConfs.reserve(picks.size());
  std::vector<double> newDummyVols;
  newDummyVols.reserve(picks.size());
  std::vector<double> newSovs;
  newSovs.reserve(picks.size());
  std::vector<double> newSofs;
  newSofs.reserve(picks.size());
  std::vector<std::vector<double>> newShifts;
  newShifts.reserve(picks.size());
  for (auto p : picks) {
    newCoords.push_back(std::move(shapes.confCoords[p]));
    newMolConfs.push_back(shapes.molConfs[p]);
    newDummyVols.push_back(shapes.dummyVols[p]);
    newSovs.push_back(shapes.sovs[p]);
    newSofs.push_back(shapes.sofs[p]);
    newShifts.push_back(shapes.shifts[p]);
  }
  shapes.confCoords = std::move(newCoords);
  shapes.molConfs = std::move(newMolConfs);
  shapes.dummyVols = std::move(newDummyVols);
  shapes.sovs = std::move(newSovs);
  shapes.sofs = std::move(newSofs);
  shapes.shifts = std::move(newShifts);
}

namespace {
// Sort the shapes in descending order of sov + sof;
void sortShapes(SearchShapeInput &shapes) {
  std::vector<std::pair<double, size_t>> vals;
  vals.reserve(shapes.confCoords.size());
  for (size_t i = 0; i < shapes.confCoords.size(); i++) {
    vals.push_back(std::make_pair(shapes.sofs[i] + shapes.sovs[i], i));
  }
  std::ranges::sort(vals,
                    [](const std::pair<double, size_t> &a,
                       const std::pair<double, size_t> &b) -> bool {
                      return a.first > b.first;
                    });
  std::vector<std::vector<float>> newCoords;
  newCoords.reserve(vals.size());
  std::vector<unsigned int> newMolConfs;
  newMolConfs.reserve(vals.size());
  std::vector<double> newDummyVols;
  newDummyVols.reserve(vals.size());
  std::vector<double> newSovs;
  newSovs.reserve(vals.size());
  std::vector<double> newSofs;
  newSofs.reserve(vals.size());
  std::vector<std::vector<double>> newShifts;
  newShifts.reserve(vals.size());
  for (auto v : vals) {
    newCoords.push_back(std::move(shapes.confCoords[v.second]));
    newMolConfs.push_back(shapes.molConfs[v.second]);
    newDummyVols.push_back(shapes.dummyVols[v.second]);
    newSovs.push_back(shapes.sovs[v.second]);
    newSofs.push_back(shapes.sofs[v.second]);
    newShifts.push_back(shapes.shifts[v.second]);
  }
  shapes.confCoords = std::move(newCoords);
  shapes.molConfs = std::move(newMolConfs);
  shapes.dummyVols = std::move(newDummyVols);
  shapes.sovs = std::move(newSovs);
  shapes.sofs = std::move(newSofs);
  shapes.shifts = std::move(newShifts);
}
}  // namespace

std::unique_ptr<SearchShapeInput> PrepareConformers(
    const RDKit::ROMol &mol, const ShapeInputOptions &shapeOpts,
    double pruneThreshold) {
  PRECONDITION(
      mol.getNumConformers() > 0,
      "SearchShapeInput object needs the molecule to have conformers.  " +
          mol.getProp<std::string>("_Name") + "  " + MolToSmiles(mol));
  auto first = PrepareConformer(mol, 0, shapeOpts);
  auto result = std::make_unique<SearchShapeInput>(first);

  result->numDummies = shapeOpts.atomRadii.size();
  result->confCoords.reserve(mol.getNumConformers());
  result->molConfs.reserve(mol.getNumConformers());
  result->dummyVols.reserve(mol.getNumConformers());
  result->sovs.reserve(mol.getNumConformers());
  result->sofs.reserve(mol.getNumConformers());
  result->shifts.reserve(mol.getNumConformers());

  ShapeInputOptions noDummyOpts(shapeOpts);
  noDummyOpts.includeDummies = false;
  noDummyOpts.atomRadii.clear();
  auto otherFirst = PrepareConformer(mol, 0, noDummyOpts);
  // A value of -1.0 will have been placed in dummyVols[0] as a placekeeper.
  result->dummyVols[0] = first.sov - otherFirst.sov;

  for (unsigned int i = 1; i < mol.getNumConformers(); i++) {
    auto shape = PrepareConformer(mol, i, shapeOpts);
    result->confCoords.push_back(shape.coord);
    result->molConfs.push_back(i);
    result->sovs.push_back(shape.sov);
    result->sofs.push_back(shape.sof);
    auto otherShape = PrepareConformer(mol, i, noDummyOpts);
    result->dummyVols.push_back(shape.sov - otherShape.sov);
    result->shifts.push_back(shape.shift);
  }
  pruneShapes(*result, pruneThreshold);
  sortShapes(*result);
  return result;
}

std::pair<double, double> bestSimilarity(SearchShapeInput &refShape,
                                         SearchShapeInput &fitShape,
                                         std::vector<float> &matrix,
                                         double threshold, double opt_param,
                                         unsigned int max_preiters,
                                         unsigned int max_postiters) {
  // The best score achievable is when the smaller volume is entirely inside
  // the larger volume.  The Shape tanimoto is the fraction of volume in
  // common.  The scores for the different conformations are sorted in
  // descending order.  So try the smallest refShape score against the
  // largest fitShape score and vice versa.
  auto maxScore = [](double v1, double v2, double f1, double f2) -> double {
    double maxSt = std::min(v1, v2) / std::max(v1, v2);
    double maxCt = std::min(f1, f2) / std::max(f1, f2);
    return maxSt + maxCt;
  };
  if (maxScore(fitShape.sovs.front(), refShape.sovs.back(),
               fitShape.sofs.front(), refShape.sofs.back()) < threshold &&
      maxScore(fitShape.sovs.back(), refShape.sovs.front(),
               fitShape.sofs.back(), refShape.sofs.front()) < threshold) {
    return std::make_pair(-1.0, -1.0);
  }
  double best_st = -1.0;
  double best_ct = -1.0;
  double best_combo_t = -1.0;

  for (size_t i = 0; i < refShape.confCoords.size(); i++) {
    refShape.setActiveShape(i);
    for (size_t j = 0; j < fitShape.confCoords.size(); j++) {
      auto maxSim = maxScore(refShape.sovs[i], fitShape.sovs[j],
                             fitShape.sofs[i], fitShape.sofs[j]);
      if (maxSim > threshold) {
        fitShape.setActiveShape(j);
        auto [st, ct] = AlignShape(refShape, fitShape, matrix, opt_param,
                                   max_preiters, max_postiters);
        double combo_t = st + ct;
        if (combo_t > threshold) {
          return std::make_pair(st, ct);
        }
        if (combo_t > best_combo_t) {
          best_combo_t = combo_t;
          best_st = st;
          best_ct = ct;
        }
      }
    }
  }

  return std::make_pair(best_st, best_ct);
}

std::unique_ptr<RDKit::RWMol> shapeToMol(const ShapeInput &shape) {
  auto mol = std::make_unique<RDKit::RWMol>();
  unsigned int numAtoms = shape.coord.size() / 3;
  for (unsigned int i = 0; i < numAtoms; i++) {
    RDKit::Atom *atom = nullptr;
    if (shape.atom_type_vector[i]) {
      atom = new RDKit::Atom(7);
    } else {
      atom = new RDKit::Atom(6);
    }
    mol->addAtom(atom, true, true);
  }
  RDKit::Conformer *conf = new RDKit::Conformer(numAtoms);
  RDGeom::Point3D ave;
  unsigned int nAves = 0;
  for (unsigned int i = 0; i < numAtoms; i++) {
    auto &pos = conf->getAtomPos(i);
    pos.x = shape.coord[3 * i];
    pos.y = shape.coord[3 * i + 1];
    pos.z = shape.coord[3 * i + 2];
    if (!shape.atom_type_vector[i]) {
      ave += pos;
      ++nAves;
    }
  }
  ave /= nAves;
  mol->addConformer(conf, true);
  return mol;
}
