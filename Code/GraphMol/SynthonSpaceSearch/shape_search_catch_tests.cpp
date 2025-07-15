//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <cstdio>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
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
  // respectively.  The first of these queries will have
  // conformers generated, the second and third will just use
  // the one given.
  std::vector<std::string> querySmis{
      "c1ccccc1C(=O)N1CCCC1",
      "C[C@H]1CCN(c2nnc(CO)n2C2CCCC2)C1 |(3.88187,-1.9608,1.02401;3.40947,-0.556473,0.685633;3.49763,-0.278493,-0.787477;2.18424,0.313597,-1.21035;1.39217,0.480256,0.0110759;0.36454,1.42337,0.278193;0.656593,2.66561,0.766052;-0.477075,3.33073,0.923413;-1.48168,2.52297,0.540741;-2.93641,2.8664,0.551747;-3.33354,3.43671,-0.657732;-0.965924,1.32594,0.134935;-1.71735,0.224133,-0.333621;-1.11549,-0.643316,-1.37025;-2.26431,-1.65079,-1.54776;-2.58926,-1.96617,-0.0975393;-2.00229,-0.836476,0.740437;1.90383,-0.551243,0.906815),wD:1.0|",
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(1.04478,-0.224363,-2.91518;1.46111,-1.07435,-1.76352;0.306282,-1.55748,-0.916018;-0.365201,-0.475593,-0.123586;-1.76672,-0.342998,-0.464794;-2.78181,-0.454019,0.51802;-2.5057,-0.664137,1.73574;-4.16163,-0.338949,0.199908;-5.17788,-0.4597,1.20931;-6.35254,0.434956,1.04892;-7.08067,-0.58036,0.437714;-6.19369,-1.54598,0.869789;0.333956,0.862517,-0.269047;1.75934,0.639384,-0.124627;2.53261,1.42843,0.743485;1.99398,2.35859,1.42205;3.96629,1.24952,0.929582;4.47592,0.048496,1.32782;5.80644,-0.103658,1.49611;6.61864,0.989857,1.25231;6.12678,2.18712,0.856797;4.80783,2.3177,0.695962;4.34378,3.52947,0.29874;2.46143,-0.410331,-0.868145),wD:1.0,3.3|",
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
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 1;

  for (size_t i = 0; i < libNames.size(); i++) {
    // if (i != 0) {
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
    params.bestHit = true;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    std::cout << "Number of query confs from " << querySmis[i]
              << " :: " << queryMol->getNumConformers() << std::endl;
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

  auto spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_confs.spc";

  bool cancelled = false;
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName, cancelled);

  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 1.0;
  shapeBuildOptions.numThreads = -1;
  shapeBuildOptions.maxSynthonAtoms = 0;
  shapeBuildOptions.maxEmbedAttempts = 10;
  shapeBuildOptions.interimFile = spaceName;
  shapeBuildOptions.interimWrites = 100;

  synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);
  synthonspace.writeDBFile(spaceName);
}

TEST_CASE("Hits back onto query") {
  // Make sure the hits are properly translated to the reference
  // frame of the query.
  SynthonSpace synthonspace;
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
O=C(c1coc2cc(Cl)ccc12)[1*:1]	A	2	r1	3
OC(CC1CCCN1[1*:1])c1ccco1	B	1	r1	3
Cc1cnc([1*:1])nc1C	C	2	r2	3
Cn1cc(N2CCCC(N[1*:1])C2)cn1	D	1	r2	3
Cc1cc(N[1*:1])cc(C(F)(F)F)c1	E	1	r3	3
Cc1nc(-c2ccccn2)cc([1*:1])n1	F	2	r3	3
CC(C)(C)C1(C)CN([1*:1])CCO1	G	1	r4	3
CCCC1CN([2*:2])CCO1	H	3	r4	3
COC1C2CCCC2C1N[2*:2]	I	3	r4	3
COCc1nc([1*:1])cc([2*:2])n1	J	2	r4	3
FC1(F)CCC(N[1*:1])CC1	K	1	r4	3)");
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
  std::cout << sizeof(void *) << std::endl;

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
O=C(c1coc2cc(Cl)ccc12)[1*:1]	A	2	r1	3
OC(CC1CCCN1[1*:1])c1ccco1	B	1	r1	3
Cc1cnc([1*:1])nc1C	C	2	r2	3
Cn1cc(N2CCCC(N[1*:1])C2)cn1	D	1	r2	3
Cc1cc(N[1*:1])cc(C(F)(F)F)c1	E	1	r3	3
Cc1nc(-c2ccccn2)cc([1*:1])n1	F	2	r3	3
CC(C)(C)C1(C)CN([1*:1])CCO1	G	1	r4	3
CCCC1CN([2*:2])CCO1	H	3	r4	3
COC1C2CCCC2C1N[2*:2]	I	3	r4	3
COCc1nc([1*:1])cc([2*:2])n1	J	2	r4	3
FC1(F)CCC(N[1*:1])CC1	K	1	r4	3)");
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

TEST_CASE("Bad elements") {
  SynthonSpace space;
  // This bit of FreedomSpace 2024-09 threw an exception originally
  // due to the Pd atom not being recognised by the shape builder.
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
Cc1sc2ncnc(NN[U])c2c1C	100000182719	1	20a	2024-09
ClC(Cl)=C(Cl)c1cccc(N[U])c1	100016989384	1	20a	2024-09
Cl[Pd]c1cccc(-c2ccccc2)c1N[U]	114813040699	1	20a	2024-09
C=CCC1(S(=O)(=O)[U])CC1	100000101377	2	20a	2024-09
C=Cc1ccc(S(=O)(=O)[U])cc1	100000034458	2	20a	2024-09
)");
  bool cancelled = false;
  space.readStream(iss, cancelled);
  ShapeBuildParams shapeOptions;
  shapeOptions.randomSeed = 1;
  CHECK_NOTHROW(space.buildSynthonShapes(cancelled, shapeOptions));
}

TEST_CASE("Shape Best Hit Found") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/small_freedom_shapes.spc";

  SynthonSpace space;
  space.readDBFile(spaceName);
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = 1;
  params.numConformers = 200;
  params.confRMSThreshold = 0.25;
  params.randomSeed = 0xdac;

  params.fragSimilarityAdjuster = 0.1;
  params.approxSimilarityAdjuster = 0.1;
  params.similarityCutoff = 1.0;
  params.timeOut = 0;

  auto mol =
      "CCC(C(=O)NCc1ccco1)N(Cc1sccc1C)C(C)C |(1.19967,-2.26511,-1.8853;-0.0674677,-1.53728,-1.54329;-0.0395195,-0.735921,-0.2957;0.978688,0.316536,-0.356493;0.610079,1.54481,-0.391485;2.37979,0.114038,-0.377756;3.31574,1.24811,-0.443092;4.723,0.774447,-0.458657;5.58509,0.534718,0.592748;6.76773,0.108395,-0.00740365;6.56938,0.109237,-1.38929;5.34763,0.509833,-1.60401;-1.33892,-0.260032,0.0922388;-2.08764,0.476264,-0.832263;-3.39845,0.92709,-0.383151;-4.99177,0.223473,-0.828611;-5.94415,1.34626,0.133161;-5.00918,2.1798,0.729651;-3.7235,1.96002,0.462253;-2.60368,2.81164,1.05297;-1.99138,-1.01276,1.09008;-1.10607,-0.965532,2.34165;-2.26411,-2.45544,0.780694)|"_smiles;
  REQUIRE(mol);
  auto results = space.shapeSearch(*mol, params);
  std::cout << "Number of hits : " << results.getHitMolecules().size()
            << std::endl;
  CHECK(results.getHitMolecules().empty());
  auto &bestHit = results.getBestHit();
  std::cout << "Best hit found : " << bestHit->getProp<double>("Similarity")
            << std::endl;
  CHECK(bestHit);
  CHECK(bestHit->getProp<double>("Similarity") < 1.0);
}
