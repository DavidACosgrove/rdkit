
rdkit_library(SpeGSS
    SpeGSS.cpp
    LINK_LIBRARIES SmilesParse GraphMol)
target_compile_definitions(SpeGSS PRIVATE RDKIT_SPEGGS_BUILD)

rdkit_headers(SpeGSS.h
    DEST GraphMol/SpeGSS)

rdkit_catch_test(testSpeGSS catch_tests.cpp LINK_LIBRARIES SpeGSS)