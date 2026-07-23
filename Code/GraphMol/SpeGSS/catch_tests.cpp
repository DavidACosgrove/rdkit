//
// Copyright (C) David Cosgrove 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SpeGSS/SpeGSS.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SpeGSS;

TEST_CASE("Generate Surface") {
  {
    const auto ethane =
        "CC |(0.8151,-0.5383,0.4928;2.1548,-0.0516,-0.0117)|"_smiles;
    REQUIRE(ethane);
    SpeGSSParams params{.isoValue = 0.5};
    auto surf = createSurface(*ethane, params);
    REQUIRE(surf);
    std::cout << "Num vertices : " << surf->get_num_vertices() << std::endl;
    std::cout << "Num triangles : " << surf->get_num_triangles() << std::endl;
    CHECK(surf->get_num_vertices() == 980);
    CHECK(surf->get_num_triangles() == 1956);
    surf->save_txt("ethane.surf");
  }
#if 1
  {
    const auto benzene =
        "c1ccccc1 |(-1.39928,0.010791,-0.0463459;-0.701518,-1.17621,-0.047681;0.676993,-1.18937,-0.00208102;1.38771,-0.0146145,0.0458731;0.719034,1.19151,0.048576;-0.646501,1.17791,0.00285385)|"_smiles;
    REQUIRE(benzene);
    SpeGSSParams params{.isoValue = 1.5};
    auto surf = createSurface(*benzene, params);
    REQUIRE(surf);
    std::cout << "Num vertices : " << surf->get_num_vertices() << std::endl;
    std::cout << "Num triangles : " << surf->get_num_triangles() << std::endl;
    CHECK(surf->get_num_vertices() == 1118);
    CHECK(surf->get_num_triangles() == 2232);
    surf->save_txt("benzene.surf");
  }
#endif
}