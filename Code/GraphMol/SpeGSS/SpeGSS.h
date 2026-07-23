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

#include <memory>

#include <RDGeneral/export.h>

#ifndef RDKIT_SPEGSS_H
#define RDKIT_SPEGSS_H

#include <GraphMol/SpeGSS/SpeGSSParams.h>

// This is the Marching Cubes 33 implementation for generating the
// surface from the Gaussian density of a molecule.
#define mc33cpp_implementation
#include <MC33/mc33/mc33cpp.h>

namespace RDKit {
class ROMol;

namespace SpeGSS {

std::unique_ptr<surface> createSurface(
    const ROMol &mol, const SpeGSSParams &params = SpeGSSParams());

}
}  // namespace RDKit

#endif  // RDKIT_SPEGSS_H
