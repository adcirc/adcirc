SET( METIS_SOURCES ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/coarsen.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/fm.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/initpart.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/match.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/ccgraph.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/memory.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/pmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/pqueue.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/refine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/util.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/timing.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/debug.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/bucketsort.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/graph.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/stat.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kwayrefine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kwayfm.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/balance.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/ometis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/srefine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/sfm.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/separator.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mincover.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mmd.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mesh.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/meshpart.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/frename.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/fortran.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/myqsort.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/compress.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/parmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/estmem.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mpmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mcoarsen.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mmatch.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/minitpart.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mbalance.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mrefine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mutil.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mfm.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mkmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mkwayrefine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mkwayfmh.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mrefine2.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/minitpart2.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mbalance2.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/mfm2.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kvmetis.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolrefine.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolfm.c 
                   ${CMAKE_SOURCE_DIR}/thirdparty/metis/Lib/subdomains.c )
ADD_LIBRARY(metis STATIC ${METIS_SOURCES})
TARGET_INCLUDE_DIRECTORIES(metis PRIVATE ${CMAKE_SOURCE_DIR}/metis/Lib)
SET_TARGET_PROPERTIES(metis PROPERTIES EXCLUDE_FROM_ALL TRUE)

#...In gcc-10+, we need to ignore some things that were elevated to errors. This is a thirdparty
#   library, so adcirc devs are in nofix mode
GET_FILENAME_COMPONENT(C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
IF(${C_COMPILER_NAME} MATCHES "gcc.*" OR ${C_COMPILER_NAME} MATCHES "cc" )
    IF( ${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_C_COMPILER_VERSION} VERSION_EQUAL 10 )
        SET_TARGET_PROPERTIES(metis PROPERTIES COMPILE_FLAGS "-Wno-implicit-function-declaration" )
    ENDIF()
ENDIF()




