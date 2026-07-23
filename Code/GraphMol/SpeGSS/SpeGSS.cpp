//
// Copyright (C) David Cosgrove 2026
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is an implementation of the Spectral Geometry Shape Similarity
// measure described in
// Alignment-Free Molecular Shape Comparison Using Spectral Geometry:
// The Framework.
// Matthew P Seddon, David A Cosgrove, Martin J Packer and Valerie J Gillet
// JCIM 59, 98-116 (2019)
// https://pubs.acs.org/doi/10.1021/acs.jcim.8b00676

#include <limits>
#include <numbers>

#include <GraphMol/ROMol.h>
#include <Geometry/point.h>
#include <GraphMol/SpeGSS/SpeGSS.h>

namespace RDKit {
namespace SpeGSS {
// Bondi radii
// You can find more of these in Table 12 of this publication:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
const std::map<unsigned int, double> vdw_radii = {
    {0, 2.16},   // Dummy, same as Xe.
    {1, 1.10},   // H
    {2, 1.40},   // He
    {3, 1.81},   // Li
    {4, 1.53},   // Be
    {5, 1.92},   // B
    {6, 1.70},   // C
    {7, 1.55},   // N
    {8, 1.52},   // O
    {9, 1.47},   // F
    {10, 1.54},  // Ne
    {11, 2.27},  // Na
    {12, 1.73},  // Mg
    {13, 1.84},  // Al
    {14, 2.10},  // Si
    {15, 1.80},  // P
    {16, 1.80},  // S
    {17, 1.75},  // Cl
    {18, 1.88},  // Ar
    {19, 2.75},  // K
    {20, 2.31},  // Ca
    {31, 1.87},  // Ga
    {32, 2.11},  // Ge
    {33, 1.85},  // As
    {34, 1.90},  // Se
    {35, 1.83},  // Br
    {36, 2.02},  // Kr
    {37, 3.03},  // Rb
    {38, 2.49},  // Sr
    {49, 1.93},  // In
    {50, 2.17},  // Sn
    {51, 2.06},  // Sb
    {52, 2.06},  // Te
    {53, 1.98},  // I
    {54, 2.16},  // Xe
    {55, 3.43},  // Cs
    {56, 2.68},  // Ba
    {81, 1.96},  // Tl
    {82, 2.02},  // Pb
    {83, 2.07},  // Bi
    {84, 1.97},  // Po
    {85, 2.02},  // At
    {86, 2.20},  // Rn
    {87, 3.48},  // Fr
    {88, 2.83},  // Ra
};

namespace {
void getMoleculeExtremes(const ROMol &mol, int confId, double &xmin,
                         double &xmax, double &ymin, double &ymax, double &zmin,
                         double &zmax, double molPadding) {
  xmin = std::numeric_limits<double>::max();
  xmax = std::numeric_limits<double>::lowest();
  ymin = std::numeric_limits<double>::max();
  ymax = std::numeric_limits<double>::lowest();
  zmin = std::numeric_limits<double>::max();
  zmax = std::numeric_limits<double>::lowest();
  const auto &conf = mol.getConformer(confId);
  for (const auto &atom : mol.atoms()) {
    const auto &pos = conf.getAtomPos(atom->getIdx());
    std::cout << atom->getIdx() << " : " << pos.x << " " << pos.y << " "
              << pos.z << std::endl;
    if (pos.x < xmin) {
      xmin = pos.x;
    }
    if (pos.x > xmax) {
      xmax = pos.x;
    }
    if (pos.y < ymin) {
      ymin = pos.y;
    }
    if (pos.y > ymax) {
      ymax = pos.y;
    }
    if (pos.z < zmin) {
      zmin = pos.z;
    }
    if (pos.z > zmax) {
      zmax = pos.z;
    }
  }
  xmin -= molPadding;
  xmax += molPadding;
  ymin -= molPadding;
  ymax += molPadding;
  zmin -= molPadding;
  zmax += molPadding;
}

void calcNumGridPoints(double xrange, double yrange, double zrange,
                       double gridSpacing, unsigned int &numX,
                       unsigned int &numY, unsigned int &numZ) {
  numX = 1 + int((xrange) / gridSpacing);
  numY = 1 + int((yrange) / gridSpacing);
  numZ = 1 + int((zrange) / gridSpacing);
}

double getStandardAtomRadius(const unsigned int atomicNum) {
  // Mostly they will be carbons, so just return that without lookup.
  if (atomicNum == 6) {
    return 1.7;
  }
  if (const auto rad = vdw_radii.find(static_cast<unsigned int>(atomicNum));
      rad != vdw_radii.end()) {
    return rad->second;
  }
  throw ValueErrorException("No VdW radius for atom with Z=" +
                            std::to_string(atomicNum));
}

void extractAlphasAndCoords(const ROMol &mol, int confId,
                            std::vector<double> &alphas,
                            std::vector<RDGeom::Point3D> &coords) {
  static double constexpr alpha_bit = std::numbers::pi / 1.3401;
  alphas.resize(mol.getNumAtoms());
  coords.resize(3 * mol.getNumAtoms());
  const auto &conf = mol.getConformer(confId);
  for (const auto &a : mol.atoms()) {
    auto rad = getStandardAtomRadius(a->getAtomicNum());
    alphas[a->getIdx()] = alpha_bit / (rad * rad);
    coords[a->getIdx()] = conf.getAtomPos(a->getIdx());
    std::cout << a->getIdx() << " : " << alphas[a->getIdx()] << " : "
              << coords[a->getIdx()] << std::endl;
  }
}

double calcDensity(const std::vector<double> &alphas,
                   const std::vector<RDGeom::Point3D> &coords, const double x,
                   const double y, const double z) {
  static double constexpr p =
      4.0 * std::numbers::pi / (3.0 * 1.5514);  // 4 * PI / 3 * lambda
  double dens = 0.0;
  RDGeom::Point3D gp{x, y, z};
  for (size_t i = 0; i < alphas.size(); ++i) {
    double alpha_dist = -alphas[i] * (coords[i] - gp).lengthSq();
    dens += p * exp(alpha_dist);
  }
  return dens;
}

void fillGrid(grid3d &grid, const ROMol &mol, const SpeGSSParams &params,
              const int confId = -1) {
  double xmin, xmax, ymin, ymax, zmin, zmax;
  getMoleculeExtremes(mol, confId, xmin, xmax, ymin, ymax, zmin, zmax,
                      params.moleculePadding + params.gridSpacing);
  std::cout << "x : " << xmin << " -> " << xmax << std::endl;
  std::cout << "y : " << ymin << " -> " << ymax << std::endl;
  std::cout << "z : " << zmin << " -> " << zmax << std::endl;
  unsigned int numX, numY, numZ;
  calcNumGridPoints(xmax - xmin, ymax - ymin, zmax - zmin, params.gridSpacing,
                    numX, numY, numZ);
  std::cout << "Num grid points : " << numX << " x " << numY << " x " << numZ
            << std::endl;
  grid.set_grid_dimensions(numX, numY, numZ);
  grid.set_r0(xmin, ymin, zmin);
  grid.set_ratio_aspect(params.gridSpacing, params.gridSpacing,
                        params.gridSpacing);

  std::vector<double> alphas;
  std::vector<RDGeom::Point3D> coords;
  double maxDens = 0.0;
  double minDens = std::numeric_limits<double>::max();

  extractAlphasAndCoords(mol, confId, alphas, coords);
  for (unsigned int i = 0; i < numX; i++) {
    double x = xmin + params.gridSpacing * i;
    for (unsigned int j = 0; j < numY; j++) {
      double y = ymin + params.gridSpacing * j;
      for (unsigned int k = 0; k < numZ; k++) {
        double z = zmin + params.gridSpacing * k;
        const auto dens = calcDensity(alphas, coords, x, y, z);
        if (dens > maxDens) {
          maxDens = dens;
        }
        if (dens < minDens) {
          minDens = dens;
        }
        grid.set_grid_value(i, j, k, dens);
      }
    }
  }
  std::cout << "Density range : " << minDens << " -> " << maxDens << std::endl;
}
}  // namespace

std::unique_ptr<surface> createSurface(const ROMol &mol,
                                       const SpeGSSParams &params) {
  PRECONDITION(mol.getNumConformers(),
               "molecule must have at least 1 conformer");
  grid3d grid;
  fillGrid(grid, mol, params);
  grid.save_raw_file("filename.dat");
  MC33 mc;
  mc.set_grid3d(grid);
  auto surf = std::make_unique<surface>();
  mc.calculate_isosurface(*surf, params.isoValue);
  return surf;
}

}  // namespace SpeGSS
}  // namespace RDKit