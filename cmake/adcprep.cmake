IF(BUILD_ADCPREP)
    
    SET( METIS_SOURCES    metis/Lib/coarsen.c metis/Lib/fm.c metis/Lib/initpart.c metis/Lib/match.c 
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
    
    SET( ADCPREP_SOURCES  src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F src/hashtable.F
                          src/mesh.F src/global_3dvs.F wind/vortex.F src/owiwind.F src/rs2.F 
                          src/owi_ice.F src/wind.F prep/presizes.F prep/pre_global.F prep/metis.F 
                          prep/subprep.F prep/adcprep.F prep/decomp.F prep/prep_weir.F 
                          src/nodalattr.F src/harm.F prep/read_global.F src/netcdfio.F prep/prep.F  
                          prep/interp.F prep/machdep.F )
    
    ADD_LIBRARY(metis STATIC ${METIS_SOURCES})
    
    ADD_EXECUTABLE(adcprep ${ADCPREP_SOURCES})
    
    SET_TARGET_PROPERTIES(adcprep PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcprep_mod)
    
    SET(ADCPREP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${MACHINE_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCPREP}")
    
    IF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${PREP_SWAN_FLAG}")
    ENDIF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
   
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcprep netcdf netcdff)
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCPREP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${ADCPREP_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcprep netcdf netcdff XdmfCore XdmfUtils Xdmf)
        TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
   
    SET_TARGET_PROPERTIES(adcprep PROPERTIES COMPILE_FLAGS ${ADCPREP_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_LINK_LIBRARIES(adcprep metis version mkdir)
    
    ADD_DEPENDENCIES(adcprep metis)
    ADD_DEPENDENCIES(adcprep mkdir)
    ADD_DEPENDENCIES(adcprep version)
    
    INSTALL(TARGETS adcprep RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCPREP)
