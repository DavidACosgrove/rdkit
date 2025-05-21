//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <cstdio>

#include <../External/pubchem_shape/PubChemShape.hpp>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;
using namespace RDKit::RascalMCES;

const char *rdbase = getenv("RDBASE");

void prepareMolecule(RWMol *mol) {
  MolOps::addHs(*mol);
  auto dgParams = DGeomHelpers::ETKDGv3;
  // dgParams.pruneRmsThresh = 1.0;
  dgParams.randomSeed = 1;
  DGeomHelpers::EmbedMultipleConfs(*mol, 100, dgParams);
  MolOps::removeHs(*mol);
}

std::map<std::string, std::unique_ptr<ROMol>> loadLibrary(
    const std::string inFilename) {
  v2::FileParsers::SmilesMolSupplierParams params;
  params.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(inFilename, params);
  std::map<std::string, std::unique_ptr<ROMol>> mols;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    if (mol) {
      prepareMolecule(mol.get());
      std::string molName = mol->getProp<std::string>(common_properties::_Name);
      mols.insert(std::make_pair(
          molName,
          std::unique_ptr<ROMol>(static_cast<ROMol *>(mol.release()))));
    }
  }
  return mols;
};

TEST_CASE("Shape Small tests") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::vector<std::string> libNames{
      fullRoot + "amide_space.txt",
      fullRoot + "triazole_space.txt",
      fullRoot + "urea_space.txt",
  };
  std::vector<std::string> enumLibNames{
      fullRoot + "amide_space_enum.smi",
      fullRoot + "triazole_space_enum.smi",
      fullRoot + "urea_space_enum.smi",
  };
  std::vector<std::string> enumOutputNames{
      "amide_space_enum_out.sdf",
      "triazole_space_enum_out.sdf",
      "urea_space_enum_out.sdf",
  };
  std::vector<std::string> searchOutputNames{
      "amide_space_search_out.sdf",
      "triazole_search_out.sdf",
      "urea_space_search_out.sdf",
  };

  // The search of the enumerated libraries give 4, 8, 4 hits
  // respectively.
  std::vector<std::string> querySmis{
      "c1ccccc1C(=O)N1CCCC1",
      "CC1CCN(c2nnc(CO)n2C2CCCC2)C1",
      "C[C@@H]1CC(NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1",
  };

  // The synthon search gives 1 hit for the urea space, where the
  // brute-force search gives 4 because the fragment similarities fall
  // below the threshold.  For example, comparing [2*]c1nccnc1F from
  // the query with synthon N#CCc(cncc1)c1[2*] (689988332-107515102)
  // when the dummy atoms are aligned, which they should be for a
  // good synthon match, the feature score is low because the nitrogen
  // acceptors don't align.  In the full molecule overlay, that is
  // compensated for by other things.
  std::vector<size_t> expNumHits{3, 5, 1};
  unsigned int numConfs = 100;
  double rmsThreshold = 1.0;
  int numThreads = 1;

  for (size_t i = 0; i < libNames.size(); i++) {
    // if (i != 0) {
    // continue;
    // }
    SynthonSpace synthonspace;
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    synthonspace.buildSynthonShapes(numConfs, rmsThreshold, numThreads);

    SynthonSpaceSearchParams params;
    params.similarityCutoff = 1.6;
    params.numConformers = numConfs;
    params.numThreads = numThreads;
    params.confRMSThreshold = rmsThreshold;
    params.timeOut = 0;
    params.randomSeed = 1;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.shapeSearch(*queryMol, params);
    std::cout << "Num hits : " << results.getHitMolecules().size() << " : "
              << results.getMaxNumResults() << std::endl;
    for (const auto &hit : results.getHitMolecules()) {
      std::cout << hit->getProp<std::string>(common_properties::_Name) << " : "
                << hit->getProp<double>("Similarity") << std::endl;
    }
    CHECK(expNumHits[i] == results.getHitMolecules().size());
    RDKit::SDWriter sdw(searchOutputNames[i]);
    for (const auto &hit : results.getHitMolecules()) {
      sdw.write(*hit);
    }
#if 0
    auto mols = loadLibrary(enumLibNames[i]);
    prepareMolecule(queryMol.get());
    RDKit::SDWriter sdw2(enumOutputNames[i]);
    std::vector<float> matrix(12, 0.0);
    unsigned int numHits = 0;
    for (auto &[smiles, mol] : mols) {
      bool foundHit = false;
      for (unsigned int i = 0; i < queryMol->getNumConformers(); ++i) {
        for (unsigned int j = 0; j < mol->getNumConformers(); ++j) {
          auto [st, ct] = AlignMolecule(*queryMol, *mol, matrix, i, j);
          if (st + ct > params.similarityCutoff) {
            std::cout << mol->getProp<std::string>(common_properties::_Name)
                      << " hit at " << st << ", " << ct << " : " << st + ct
                      << " for " << i << ", " << j << std::endl;
            ++numHits;
            foundHit = true;
            sdw2.write(*mol);
            break;
          }
        }
        if (foundHit) {
          break;
        }
      }
    }
#endif
  }
}

TEST_CASE("Shape DB Writer") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumReactions() == 1);
  synthonspace.buildSynthonShapes();

  auto spaceName = std::tmpnam(nullptr);

  synthonspace.writeDBFile(spaceName);

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile(spaceName);
  CHECK(newsynthonspace.getNumReactions() == 1);
  std::shared_ptr<SynthonSet> irxn;
  CHECK_NOTHROW(irxn = newsynthonspace.getReaction("doebner-miller-quinoline"));

  const auto &orxn = synthonspace.getReaction("doebner-miller-quinoline");
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    REQUIRE(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons().size(); ++j) {
      REQUIRE(
          irxn->getSynthons()[i][j].second->getShapes()->confCoords.size() ==
          orxn->getSynthons()[i][j].second->getShapes()->confCoords.size());
      for (size_t k = 0;
           k < irxn->getSynthons()[i][j].second->getShapes()->confCoords.size();
           ++k) {
        const auto ishape = irxn->getSynthons()[i][j].second->getShapes().get();
        const auto oshape = orxn->getSynthons()[i][j].second->getShapes().get();
        CHECK(ishape->sovs[k] == Catch::Approx(oshape->sovs[k]));
      }
    }
  }
}

TEST_CASE("Build conformer DB") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  libName =
      "/Users/david/Projects/SynthonSpaceTests/FreedomSpace/2024-09_Freedom_synthons.txt";
  bool cancelled = false;
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName, cancelled);
  std::cout << "Number of reactions " << synthonspace.getNumReactions()
            << std::endl;
  std::cout << "Number of products : " << synthonspace.getNumProducts()
            << std::endl;
  // synthonspace.buildSynthonShapes(100, 1.0, -1);
  auto spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_confs.spc";
  synthonspace.writeDBFile(spaceName);
}

TEST_CASE("Hits back onto query") {
  // Make sure the hits are properly translated to the reference
  // frame of the query.
  SynthonSpace synthonspace;
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
[2H]C([2H])([2H])Oc1ccc(N[1*:1])cc1	iz40kScoVtWPpF4jD0a9CQ	1	m_27bcb	3
c1cc(-c2ccsc2)nc([1*:1])n1	t17csEp9XHjGy__M0m_BaA	2	m_27bcb	3
FC1(F)CCC(N[1*:1])CC1	Fq0QBDWKFd1IEAFgT9fo9Q	1	m_282030abb	3
COCc1nc([1*:1])cc([2*:2])n1	G5GZo2pyGFPFUrga0tLhmQ	2	m_282030abb	3
CCCC1CN([2*:2])CCO1	AB2bmNAkx_loJm9IO9xV4w	3	m_282030abb	3
CC1CC(N[1*:1])CCS1	AtgNFHa8gpi1jiwP9pX30g	1	m_27bbd	3
c1cnc(-c2nccc([1*:1])n2)nc1	1V-_7VgP1WANh7JxGLtcEg	2	m_27bbd	3
C#CC1(CN[1*:1])CCCC1	MfUUSSfGb-glCBjJJ8Pk6A	1	m_27bbh	3
Brc1ccc(-c2csc([1*:1])n2)cc1	LO-9o3amJr7aNrTbxWz2xA	2	m_27bbh	3
FC1(F)CCC(N[1*:1])CC1	Fq0QBDWKFd1IEAFgT9fo9Q	1	m_282030abb	3
COCc1nc([1*:1])cc([2*:2])n1	G5GZo2pyGFPFUrga0tLhmQ	2	m_282030abb	3
CC1(C(F)(F)F)CCCN([2*:2])C1	xJny03i9orYjBwsnkBp1Sg	3	m_282030abb	3)");
  bool cancelled = false;
  synthonspace.readStream(iss, cancelled);
  std::cout << "Number of reactions " << synthonspace.getNumReactions()
            << std::endl;
  std::cout << "Number of products : " << synthonspace.getNumProducts()
            << std::endl;

  unsigned int numConfs = 100;
  double rmsThreshold = 1.0;
  int numThreads = -1;
  int seed = 1;
  synthonspace.buildSynthonShapes(numConfs, rmsThreshold, numThreads, seed);

  auto tagrisso_pdb_core =
      "c1cc(Nc2nccc(c3cn(C)c4ccccc34)n2)ccc1 |(-30.966,18.467,-10.003;-29.741,18.8,-10.881;-29.776,18.58,-12.402;-28.626,18.878,-13.264;-27.858,20.11,-13.139;-26.809,20.446,-14.135;-26.039,21.676,-14.006;-26.301,22.606,-12.864;-27.356,22.266,-11.866;-27.643,23.19,-10.674;-26.776,24.159,-10.172;-27.396,24.761,-9.099;-26.842,25.83,-8.286;-28.633,24.178,-8.929;-29.782,24.445,-7.884;-31.052,23.635,-7.939;-31.218,22.587,-8.984;-30.11,22.344,-9.979;-28.784,23.198,-9.912;-28.114,21.037,-12.005;-31.044,18.019,-13.045;-32.253,17.694,-12.176;-32.227,17.912,-10.676)|"_smiles;
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = 1;
  params.similarityCutoff = 1.2;
  params.numConformers = 100;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = 1;

  RDGeom::Point3D tag_centre;
  for (const auto atom : tagrisso_pdb_core->atoms()) {
    tag_centre += tagrisso_pdb_core->getConformer().getAtomPos(atom->getIdx());
  }
  tag_centre /= tagrisso_pdb_core->getNumAtoms();
  // The random nature of the conformation generation etc means that we don't
  // always get a hit.
  for (int i = 0; i < 5; ++i) {
    auto results = synthonspace.shapeSearch(*tagrisso_pdb_core, params);
    std::string outFileName = "tagrisso_core_hits.sdf";
    if (!results.getHitMolecules().empty()) {
      std::cout << "Writing " << results.getHitMolecules().size() << " to "
                << outFileName << std::endl;
      SDWriter writer(outFileName);
      for (const auto &m : results.getHitMolecules()) {
        writer.write(*m);
        RDGeom::Point3D hit_centre;
        for (const auto atom : m->atoms()) {
          hit_centre += m->getConformer().getAtomPos(atom->getIdx());
        }
        hit_centre /= m->getNumAtoms();
        std::cout << "hit_centre : " << hit_centre << " : "
                  << m->getProp<std::string>("Similarity") << " : "
                  << m->getProp<std::string>("_Name") << std::endl;
        CHECK((hit_centre - tag_centre).length() < 2.0);
      }
      break;
    }
  }
}
