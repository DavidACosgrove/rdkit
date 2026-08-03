
rdkit_library(SpeGSS
    SpeGSSSurface.cpp
    LINK_LIBRARIES SmilesParse GraphMol)
target_compile_definitions(SpeGSS PRIVATE RDKIT_SPEGGS_BUILD)

rdkit_headers(SpeGSSSurface.h SpeGSSParams.h
    DEST GraphMol/SpeGSS)

rdkit_catch_test(testSpeGSS catch_tests.cpp LINK_LIBRARIES SpeGSS)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
