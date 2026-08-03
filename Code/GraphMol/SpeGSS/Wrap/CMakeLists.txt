rdkit_python_extension(rdSpeGSS
        rdSpeGSS.cpp
        DEST Chem
        LINK_LIBRARIES SpeGSS)

#add_pytest(pySynthonSpaceSearch ${CMAKE_CURRENT_SOURCE_DIR}/testSynthonSpaceSearch.py)
