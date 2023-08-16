set(METIS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/coarsen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/fm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/initpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/match.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/ccgraph.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/memory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/pmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/pqueue.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/refine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/util.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/timing.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/debug.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/bucketsort.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/graph.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/stat.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/balance.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/ometis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/srefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/sfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/separator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mincover.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mmd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mesh.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/meshpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/frename.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/fortran.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/myqsort.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/compress.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/parmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/estmem.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mpmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mcoarsen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mmatch.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/minitpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mbalance.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkwayrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkwayfmh.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mrefine2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/minitpart2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mbalance2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mfm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kvmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/subdomains.c)
add_library(metis STATIC ${METIS_SOURCES})
target_include_directories(metis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/metis/Lib)
set_target_properties(metis PROPERTIES EXCLUDE_FROM_ALL TRUE)

# ...In gcc-10+, we need to ignore some things that were elevated to errors. This is a thirdparty library, so adcirc
# devs are in nofix mode
get_filename_component(C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
if(${C_COMPILER_NAME} MATCHES "gcc.*" OR ${C_COMPILER_NAME} MATCHES "cc")
  if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_C_COMPILER_VERSION} VERSION_EQUAL 10)
    set_target_properties(metis PROPERTIES COMPILE_FLAGS "-Wno-implicit-function-declaration")
  endif()
endif()
