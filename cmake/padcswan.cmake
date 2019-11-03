
IF(BUILD_PADCSWAN AND PERL_FOUND)

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
     
    SET( PADCSWAN1_SOURCES  src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F src/global_3dvs.F
                            src/messenger.F src/mesh.F src/harm.F wind/vortex.F src/wind.F src/hashtable.F 
                            src/owiwind.F src/rs2.F src/owi_ice.F src/itpackv.F src/nodalattr.F src/globalio.F 
                            src/write_output.F src/writer.F src/couple2swan.F src/netcdfio.F src/subdomain.F 
                            src/wetdry.F src/gwce.F src/momentum.F src/xdmfio.F src/control.F )

    SET( PADCSWAN_SOURCES   src/couple2swan.F src/adcirc.F src/subdomain.F 
                            src/weir_boundary.F src/read_input.F src/cstart.F src/hstart.F 
                            src/timestep.F src/vsmy.F src/transport.F src/driver.F )
    
    #...SWAN Configuration
    swanConfigurePadcswan()

    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source)

    ADD_LIBRARY(templib_swan1parallel STATIC ${SWAN1PARALLEL_SOURCES})
    ADD_LIBRARY(templib_swan2parallel STATIC ${SWAN2PARALLEL_SOURCES})
    ADD_LIBRARY(templib_padcswan1     STATIC ${PADCSWAN1_SOURCES})
    ADD_EXECUTABLE(padcswan                     ${PADCSWAN_SOURCES})

    addCompilerFlags(templib_padcswan1)
    addCompilerFlags(padcswan)
    addCompilerFlagsSwan(templib_swan1parallel ${ADDITIONAL_FLAGS_SWAN})
    addCompilerFlagsSwan(templib_swan2parallel ${ADDITIONAL_FLAGS_SWAN})

    addLibMkdir(templib_padcswan1)
    addLibMkdir(padcswan)

    addLibVersion(padcswan)
    
    addMPI(templib_padcswan1)
    addMPI(padcswan)
    addMPI(templib_swan1parallel)
    addMPI(templib_swan2parallel)

    TARGET_COMPILE_DEFINITIONS(templib_padcswan1 PRIVATE CSWAN)
    TARGET_COMPILE_DEFINITIONS(padcswan          PRIVATE CSWAN)

    TARGET_INCLUDE_DIRECTORIES(templib_padcswan1     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcswan1)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
    TARGET_INCLUDE_DIRECTORIES(padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcswan1)
    TARGET_INCLUDE_DIRECTORIES(padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
    TARGET_INCLUDE_DIRECTORIES(padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan2parallel)

    TARGET_LINK_LIBRARIES(padcswan templib_swan2parallel templib_padcswan1 templib_swan1parallel) 
   
    ADD_DEPENDENCIES(padcswan              templib_swan2parallel templib_padcswan1)
    ADD_DEPENDENCIES(templib_swan2parallel templib_padcswan1 templib_swan1parallel)
    ADD_DEPENDENCIES(templib_padcswan1     mkdir version templib_swan1parallel )
    
    INSTALL(TARGETS padcswan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_PADCSWAN AND PERL_FOUND)
