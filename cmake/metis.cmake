SET( METIS_SOURCES metis/Lib/coarsen.c metis/Lib/fm.c metis/Lib/initpart.c metis/Lib/match.c 
                  metis/Lib/ccgraph.c metis/Lib/memory.c metis/Lib/pmetis.c metis/Lib/pqueue.c 
                  metis/Lib/refine.c metis/Lib/util.c metis/Lib/timing.c metis/Lib/debug.c 
                  metis/Lib/bucketsort.c metis/Lib/graph.c metis/Lib/stat.c metis/Lib/kmetis.c 
                  metis/Lib/kwayrefine.c metis/Lib/kwayfm.c metis/Lib/balance.c metis/Lib/ometis.c 
                  metis/Lib/srefine.c metis/Lib/sfm.c metis/Lib/separator.c metis/Lib/mincover.c 
                  metis/Lib/mmd.c metis/Lib/mesh.c metis/Lib/meshpart.c metis/Lib/frename.c 
                  metis/Lib/fortran.c metis/Lib/myqsort.c metis/Lib/compress.c metis/Lib/parmetis.c 
                  metis/Lib/estmem.c metis/Lib/mpmetis.c metis/Lib/mcoarsen.c metis/Lib/mmatch.c 
                  metis/Lib/minitpart.c metis/Lib/mbalance.c metis/Lib/mrefine.c metis/Lib/mutil.c 
                  metis/Lib/mfm.c metis/Lib/mkmetis.c metis/Lib/mkwayrefine.c metis/Lib/mkwayfmh.c 
                  metis/Lib/mrefine2.c metis/Lib/minitpart2.c metis/Lib/mbalance2.c metis/Lib/mfm2.c 
                  metis/Lib/kvmetis.c metis/Lib/kwayvolrefine.c metis/Lib/kwayvolfm.c 
                  metis/Lib/subdomains.c )
ADD_LIBRARY(metis STATIC ${METIS_SOURCES})
TARGET_INCLUDE_DIRECTORIES(metis PRIVATE ${CMAKE_SOURCE_DIR}/metis/Lib)
SET_TARGET_PROPERTIES(metis PROPERTIES EXCLUDE_FROM_ALL TRUE)
