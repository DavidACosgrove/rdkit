rdkit_headers(AndQuery.h
              EqualityQuery.h
              GreaterEqualQuery.h
              GreaterQuery.h
              LessEqualQuery.h
              LessQuery.h
              NullQueryAlgebra.h
              OrQuery.h
              Query.h
              QueryObjects.h
              RangeQuery.h
              SetQuery.h
              XOrQuery.h DEST Query)

rdkit_catch_test(testQuery test.cpp LINK_LIBRARIES RDGeneral)

