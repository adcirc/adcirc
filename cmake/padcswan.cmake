
IF(BUILD_PADCSWAN AND PERL_FOUND)

    #...Note that we need to build PADCSWAN in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    ADD_EXECUTABLE(padcswan src/driver.F)
   
    SET( SWAN1PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swmod1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/couple2adcirc.f90 )
    
    SET( SWAN2PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanmain.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanpre2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom4.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanout1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanser.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpids.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpmix.f )
     
    SET( PADCSWAN1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                            src/global.F src/boundaries.F src/global_3dvs.F
                            src/messenger.F )

    SET( PADCSWAN2_SOURCES  src/mesh.F src/harm.F wind/vortex.F src/wind.F 
                            src/owiwind.F src/rs2.F src/owi_ice.F 
                            src/itpackv.F src/nodalattr.F src/globalio.F 
                            src/netcdfio.F src/control.F src/xdmfio.F )

    SET( PADCSWAN3_SOURCES  src/writer.F )
    
    SET( PADCSWAN4_SOURCES  src/write_output.F src/couple2swan.F )
    
    SET( PADCSWAN5_SOURCES  src/adcirc.F src/subdomain.F 
                            src/weir_boundary.F src/read_input.F src/cstart.F 
                            src/hstart.F src/timestep.F src/vsmy.F 
                            src/transport.F src/driver.F )
    
    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel SWAN Sources..."
        )
    ENDIF(WIN32)

    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source)

    ADD_LIBRARY(templib_swan1parallel STATIC ${SWAN1PARALLEL_SOURCES})
    ADD_LIBRARY(templib_swan2parallel STATIC ${SWAN2PARALLEL_SOURCES})
    ADD_LIBRARY(templib_padcswan1     STATIC ${PADCSWAN1_SOURCES})
    ADD_LIBRARY(templib_padcswan2     STATIC ${PADCSWAN2_SOURCES})
    ADD_LIBRARY(templib_padcswan3     STATIC ${PADCSWAN3_SOURCES})
    ADD_LIBRARY(templib_padcswan4     STATIC ${PADCSWAN4_SOURCES})
    ADD_LIBRARY(templib_padcswan5     STATIC ${PADCSWAN5_SOURCES})
   
    SET(PADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel
                                       mkdir netcdf netcdff ${MPI_Fortran_LIBRARIES}) 

    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
   
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${PADCSWAN_LINKER_FLAGS} )
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel
                                       mkdir netcdf netcdff XdmfCore XdmfUtils Xdmf
                                       ${MPI_Fortran_LIBRARIES})
        
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel mkdir
                                       ${MPI_Fortran_LIBRARIES})

    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
    SET_TARGET_PROPERTIES(templib_swan1parallel PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_swan2parallel PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan1     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan2     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan3     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan4     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan5     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan              PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan1     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan2     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan3     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan4     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan5     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan              PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    
    TARGET_INCLUDE_DIRECTORIES(padcswan              PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan1     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan2     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan3     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan4     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan5     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_INCLUDE_DIRECTORIES(padcswan              PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan1     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan2     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan3     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan4     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan5     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(swan1parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(swan2parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF()
    
    ADD_DEPENDENCIES(padcswan              templib_swan2parallel)
    ADD_DEPENDENCIES(padcswan              templib_padcswan3)
    ADD_DEPENDENCIES(templib_swan2parallel templib_padcswan5)
    ADD_DEPENDENCIES(templib_padcswan5     templib_padcswan4)
    ADD_DEPENDENCIES(templib_padcswan4     templib_padcswan2)
    ADD_DEPENDENCIES(templib_padcswan4     templib_padcswan3)
    ADD_DEPENDENCIES(templib_padcswan3     templib_padcswan2)
    ADD_DEPENDENCIES(templib_padcswan2     templib_padcswan1)
    ADD_DEPENDENCIES(templib_padcswan1     templib_swan1parallel)
    ADD_DEPENDENCIES(templib_padcswan1     mkdir)
    
    INSTALL(TARGETS padcswan RUNTIME DESTINATION bin)

ENDIF(BUILD_PADCSWAN AND PERL_FOUND)
