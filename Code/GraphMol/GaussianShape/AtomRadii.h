//
// Created by David Cosgrove on 03/09/2026.
//

#ifndef RDKIT_ATOMRADII_H
#define RDKIT_ATOMRADII_H

#include <map>

#include <RDGeneral/Exceptions.h>

constexpr double CARBON_RAD = 1.70;
constexpr double DUMMY_RAD = 2.16;  // same as Xe

// Bondi radii
// You can find more of these in Table 12 of this publication:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
const static std::map<unsigned int, double> vdw_radii = {
    {0, DUMMY_RAD},   // Dummy, same as Xe.
    {1, 1.10},        // H
    {2, 1.40},        // He
    {3, 1.81},        // Li
    {4, 1.53},        // Be
    {5, 1.92},        // B
    {6, CARBON_RAD},  // C
    {7, 1.55},        // N
    {8, 1.52},        // O
    {9, 1.47},        // F
    {10, 1.54},       // Ne
    {11, 2.27},       // Na
    {12, 1.73},       // Mg
    {13, 1.84},       // Al
    {14, 2.10},       // Si
    {15, 1.80},       // P
    {16, 1.80},       // S
    {17, 1.75},       // Cl
    {18, 1.88},       // Ar
    {19, 2.75},       // K
    {20, 2.31},       // Ca
    {31, 1.87},       // Ga
    {32, 2.11},       // Ge
    {33, 1.85},       // As
    {34, 1.90},       // Se
    {35, 1.83},       // Br
    {36, 2.02},       // Kr
    {37, 3.03},       // Rb
    {38, 2.49},       // Sr
    {49, 1.93},       // In
    {50, 2.17},       // Sn
    {51, 2.06},       // Sb
    {52, 2.06},       // Te
    {53, 1.98},       // I
    {54, 2.16},       // Xe
    {55, 3.43},       // Cs
    {56, 2.68},       // Ba
    {81, 1.96},       // Tl
    {82, 2.02},       // Pb
    {83, 2.07},       // Bi
    {84, 1.97},       // Po
    {85, 2.02},       // At
    {86, 2.20},       // Rn
    {87, 3.48},       // Fr
    {88, 2.83},       // Ra
};

static double getStandardAtomRadius(const unsigned int atomicNum) {
  // Mostly they will be carbons, so just return that without lookup.
  if (atomicNum == 6) {
    return CARBON_RAD;
  }
  if (const auto rad = vdw_radii.find(static_cast<unsigned int>(atomicNum));
      rad != vdw_radii.end()) {
    return rad->second;
  }
  throw ValueErrorException("No VdW radius for atom with Z=" +
                            std::to_string(atomicNum));
}

#endif  // RDKIT_ATOMRADII_H
