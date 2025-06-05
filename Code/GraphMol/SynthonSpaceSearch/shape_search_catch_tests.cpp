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
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
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
      "C[C@H]1CCN(c2nnc(CO)n2C2CCCC2)C1",
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
  std::vector<size_t> expNumHits{3, 4, 1};
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 1.0;
  shapeBuildOptions.numThreads = 1;

  for (size_t i = 0; i < libNames.size(); i++) {
    // if (i != 1) {
    // continue;
    // }
    SynthonSpace synthonspace;
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);

    SynthonSpaceSearchParams params;
    params.similarityCutoff = 1.6;
    params.numConformers = shapeBuildOptions.numConfs;
    params.numThreads = shapeBuildOptions.numThreads;
    params.confRMSThreshold = shapeBuildOptions.rmsThreshold;
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
  synthonspace.buildSynthonShapes(cancelled);

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
      "/Users/david/Projects/SynthonSpaceTests/REAL/2024-09_RID-4-Cozchemix/random_real_1.txt";
  bool cancelled = false;
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName, cancelled);
  std::cout << "Number of reactions " << synthonspace.getNumReactions()
            << std::endl;
  std::cout << "Number of products : " << synthonspace.getNumProducts()
            << std::endl;
  std::cout << "Number of unique synthons : " << synthonspace.getNumSynthons()
            << std::endl;

  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 1.0;
  shapeBuildOptions.numThreads = -1;

  synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);
  auto spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_confs.spc";
  spaceName =
      "/Users/david/Projects/SynthonSpaceTests/REAL/2024-09_RID-4-Cozchemix/random_real_1_shapes.spc";
  std::cout << "writing to " << spaceName << std::endl;
  synthonspace.writeDBFile(spaceName);
}

TEST_CASE("Hits back onto query") {
  // Make sure the hits are properly translated to the reference
  // frame of the query.
  SynthonSpace synthonspace;
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
O=C(c1coc2cc(Cl)ccc12)[1*:1]	IncOd2mbQnsC2WFE9PMsTA	2	m_22dbk	3
OC(CC1CCCN1[1*:1])c1ccco1	EOreohGYwWx7O_cjcanAew	1	m_22dbk	3
Cc1cnc([1*:1])nc1C	vFcPIiVuZ0Al-L1KW7ldug	2	m_27bbd	3
Cn1cc(N2CCCC(N[1*:1])C2)cn1	N0LeCStfWMnRJMSRZOYwtA	1	m_27bbd	3
Cc1cc(N[1*:1])cc(C(F)(F)F)c1	rh9OClWoshYUVeH-i89h_Q	1	m_27bcb	3
Cc1nc(-c2ccccn2)cc([1*:1])n1	nzBmAGnPKQoiQvq1zBbDqA	2	m_27bcb	3
CC(C)(C)C1(C)CN([1*:1])CCO1	j0k-TeXWtaeUeUx8eT_pow	1	m_282030abb	3
CCCC1CN([2*:2])CCO1	AB2bmNAkx_loJm9IO9xV4w	3	m_282030abb	3
COC1C2CCCC2C1N[2*:2]	yrGYN7qlqhFPPC4BUMEUSQ	3	m_282030abb	3
COCc1nc([1*:1])cc([2*:2])n1	G5GZo2pyGFPFUrga0tLhmQ	2	m_282030abb	3
FC1(F)CCC(N[1*:1])CC1	Fq0QBDWKFd1IEAFgT9fo9Q	1	m_282030abb	3)");
  bool cancelled = false;
  synthonspace.readStream(iss, cancelled);
  std::cout << "Number of reactions " << synthonspace.getNumReactions()
            << std::endl;
  std::cout << "Number of products : " << synthonspace.getNumProducts()
            << std::endl;
  synthonspace.writeEnumeratedFile("tagrisso_hits_space_enum.smi");
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 1.0;
  shapeBuildOptions.numThreads = -1;
  shapeBuildOptions.randomSeed = -1;
  synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);

  auto tagrisso_pdb_core =
      "c1cc(Nc2nccc(c3cn(C)c4ccccc34)n2)ccc1 |(-30.966,18.467,-10.003;-29.741,18.8,-10.881;-29.776,18.58,-12.402;-28.626,18.878,-13.264;-27.858,20.11,-13.139;-26.809,20.446,-14.135;-26.039,21.676,-14.006;-26.301,22.606,-12.864;-27.356,22.266,-11.866;-27.643,23.19,-10.674;-26.776,24.159,-10.172;-27.396,24.761,-9.099;-26.842,25.83,-8.286;-28.633,24.178,-8.929;-29.782,24.445,-7.884;-31.052,23.635,-7.939;-31.218,22.587,-8.984;-30.11,22.344,-9.979;-28.784,23.198,-9.912;-28.114,21.037,-12.005;-31.044,18.019,-13.045;-32.253,17.694,-12.176;-32.227,17.912,-10.676)|"_smiles;
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = 1;
  params.similarityCutoff = 1.0;
  params.numConformers = 100;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = -1;

  RDGeom::Point3D tag_centre;
  for (const auto atom : tagrisso_pdb_core->atoms()) {
    tag_centre += tagrisso_pdb_core->getConformer().getAtomPos(atom->getIdx());
  }
  tag_centre /= tagrisso_pdb_core->getNumAtoms();
  // The random nature of the conformation generation etc means that we don't
  // always get a hit.
  bool foundHit = false;
  for (int i = 0; i < 1; ++i) {
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
        foundHit = true;
      }
      break;
    } else {
      std::cout << "No hits in run " << i << std::endl;
    }
  }
  CHECK(foundHit);
}

TEST_CASE("Unspecified Stereo") {
  auto m1 = "C[C@H](Cl)CCOOC(Cl)F"_smiles;
  REQUIRE(m1);
  CHECK(details::hasUnspecifiedStereo(*m1) == true);
  CHECK(details::countChiralAtoms(*m1) == 2);

  auto m2 = "C[C@H](Cl)CCOO[C@@H](Cl)F"_smiles;
  REQUIRE(m2);
  CHECK(details::hasUnspecifiedStereo(*m2) == false);
  CHECK(details::countChiralAtoms(*m2) == 2);

  auto m3 = "C[C@H](Cl)CCOO[C@@](Cl)(F)CC=CC"_smiles;
  REQUIRE(m3);
  CHECK(details::hasUnspecifiedStereo(*m3) == true);

  auto m4 = R"(C[C@H](Cl)CCOO[C@@](Cl)(F)C\C=C/C)"_smiles;
  REQUIRE(m4);
  CHECK(details::hasUnspecifiedStereo(*m4) == false);

  SynthonSpace space;
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
O=C(c1coc2cc(Cl)ccc12)[1*:1]	IncOd2mbQnsC2WFE9PMsTA	2	m_22dbk	3
OC(CC1CCCN1[1*:1])c1ccco1	EOreohGYwWx7O_cjcanAew	1	m_22dbk	3
Cc1cnc([1*:1])nc1C	vFcPIiVuZ0Al-L1KW7ldug	2	m_27bbd	3
Cn1cc(N2CCCC(N[1*:1])C2)cn1	N0LeCStfWMnRJMSRZOYwtA	1	m_27bbd	3
Cc1cc(N[1*:1])cc(C(F)(F)F)c1	rh9OClWoshYUVeH-i89h_Q	1	m_27bcb	3
Cc1nc(-c2ccccn2)cc([1*:1])n1	nzBmAGnPKQoiQvq1zBbDqA	2	m_27bcb	3
CC(C)(C)C1(C)CN([1*:1])CCO1	j0k-TeXWtaeUeUx8eT_pow	1	m_282030abb	3
CCCC1CN([2*:2])CCO1	AB2bmNAkx_loJm9IO9xV4w	3	m_282030abb	3
COC1C2CCCC2C1N[2*:2]	yrGYN7qlqhFPPC4BUMEUSQ	3	m_282030abb	3
COCc1nc([1*:1])cc([2*:2])n1	G5GZo2pyGFPFUrga0tLhmQ	2	m_282030abb	3
FC1(F)CCC(N[1*:1])CC1	Fq0QBDWKFd1IEAFgT9fo9Q	1	m_282030abb	3)");
  bool cancelled = false;
  space.readStream(iss, cancelled);

  std::ostringstream oss;
  space.enumerateToStream(oss);
  std::cout << oss.str() << std::endl;
  ShapeBuildParams shapeOptions;
  shapeOptions.randomSeed = 1;
  space.buildSynthonShapes(cancelled, shapeOptions);

  SynthonSpaceSearchParams params;
  params.similarityCutoff = 1.6;
  params.enumerateUnspecifiedStereo = false;

  // This should bale with no results because there's unspecified
  // stereochem.
  auto results = space.shapeSearch(*m1, params);
  CHECK(results.getHitMolecules().empty());

  // This is one of the molecules in the library, so should always
  // be a hit.
  auto m5 = R"(Cc1cnc(NC2CCCN(c3cnn(C)c3)C2)nc1C)"_smiles;
  REQUIRE(m5);
  CHECK(details::hasUnspecifiedStereo(*m5) == true);

  params.enumerateUnspecifiedStereo = true;
  params.randomSeed = 1;
  results = space.shapeSearch(*m5, params);
  REQUIRE(results.getHitMolecules().size() == 1);
  auto &hitMol1 = results.getHitMolecules().front();
  std::cout << hitMol1->getProp<std::string>(common_properties::_Name) << " : "
            << hitMol1->getProp<double>("Similarity") << " : "
            << hitMol1->getProp<std::string>("Query_CXSmiles") << std::endl;
  double firstSim = hitMol1->getProp<double>("Similarity");

  params.bestHit = true;
  results = space.shapeSearch(*m5, params);
  REQUIRE(results.getHitMolecules().size() == 1);
  auto &hitMol2 = results.getHitMolecules().front();
  std::cout << hitMol2->getProp<std::string>(common_properties::_Name) << " : "
            << hitMol2->getProp<double>("Similarity") << " : "
            << hitMol2->getProp<std::string>("Query_CXSmiles") << std::endl;
  double bestSim = hitMol2->getProp<double>("Similarity");
  CHECK(bestSim > firstSim);
}
