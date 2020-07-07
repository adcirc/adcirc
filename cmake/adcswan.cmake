
IF(BUILD_ADCSWAN AND PERL_FOUND)

    SET( SWAN1SERIAL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swmod1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swmod2.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/m_constants.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/serv_xnl4v5.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGriddata.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanCompdata.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/couple2adcirc.f90 )
    
    SET( SWAN2SERIAL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanmain.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanpre1.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanpre2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swancom1.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swancom2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swancom3.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swancom4.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swancom5.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanout1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanout2.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanser.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/swanparll.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanReadADCGrid.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanReadEasymeshGrid.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanCheckGrid.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGridTopology.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGridCell.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPrintGridInfo.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPointinMesh.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPrepComp.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanCompUnstruc.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPropvelX.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPropvelS.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanTranspX.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGSECorr.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanInterpolatePoint.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanInterpolateOutput.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanConvStopc.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanCrossObstacle.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanIntgratSpc.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanReadfort18.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanSumOverNodes.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanMaxOverNodes.f90 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/ocpids.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/ocpcre.f 
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/ocpmix.f ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SdsBabanin.f90  
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/SwanGradVel.f90 )
    
    SET( ADCSWAN1_SOURCES    ${CMAKE_SOURCE_DIR}/src/sizes.F 
                             ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F 
                             ${CMAKE_SOURCE_DIR}/src/global.F 
                             ${CMAKE_SOURCE_DIR}/src/boundaries.F 
                             ${CMAKE_SOURCE_DIR}/src/mesh.F 
                             ${CMAKE_SOURCE_DIR}/src/hashtable.F
                             ${CMAKE_SOURCE_DIR}/src/global_3dvs.F 
                             ${CMAKE_SOURCE_DIR}/src/harm.F 
                             ${CMAKE_SOURCE_DIR}/wind/vortex.F 
                             ${CMAKE_SOURCE_DIR}/src/wind.F 
                             ${CMAKE_SOURCE_DIR}/src/owiwind.F 
                             ${CMAKE_SOURCE_DIR}/src/rs2.F
                             ${CMAKE_SOURCE_DIR}/src/owi_ice.F 
                             ${CMAKE_SOURCE_DIR}/src/itpackv.F 
                             ${CMAKE_SOURCE_DIR}/src/nodalattr.F 
                             ${CMAKE_SOURCE_DIR}/src/globalio.F 
                             ${CMAKE_SOURCE_DIR}/src/netcdfio.F 
                             ${CMAKE_SOURCE_DIR}/src/subdomain.F 
                             ${CMAKE_SOURCE_DIR}/src/gwce.F 
                             ${CMAKE_SOURCE_DIR}/src/wetdry.F 
                             ${CMAKE_SOURCE_DIR}/src/momentum.F 
                             ${CMAKE_SOURCE_DIR}/src/control.F 
                             ${CMAKE_SOURCE_DIR}/src/xdmfio.F 
                             ${CMAKE_SOURCE_DIR}/src/write_output.F 
                             ${CMAKE_SOURCE_DIR}/src/couple2swan.F 
                             ${CMAKE_SOURCE_DIR}/src/sponge_layer.F 
                             ${CMAKE_SOURCE_DIR}/src/quadtrature.F 
                             ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F)

    SET( ADCSWAN_SOURCES     ${CMAKE_SOURCE_DIR}/src/adcirc.F 
                             ${CMAKE_SOURCE_DIR}/src/weir_boundary.F 
                             ${CMAKE_SOURCE_DIR}/src/read_input.F 
                             ${CMAKE_SOURCE_DIR}/src/cstart.F 
                             ${CMAKE_SOURCE_DIR}/src/hstart.F 
                             ${CMAKE_SOURCE_DIR}/src/timestep.F 
                             ${CMAKE_SOURCE_DIR}/src/vsmy.F 
                             ${CMAKE_SOURCE_DIR}/src/transport.F 
                             ${CMAKE_SOURCE_DIR}/src/driver.F )
    
    #...SWAN Configuration
    swanConfigureAdcswan()

    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source)

    ADD_LIBRARY(templib_swan1serial STATIC ${SWAN1SERIAL_SOURCES})
    ADD_LIBRARY(templib_swan2serial STATIC ${SWAN2SERIAL_SOURCES})
    ADD_LIBRARY(templib_adcswan1 STATIC ${ADCSWAN1_SOURCES})
    ADD_EXECUTABLE(adcswan ${ADCSWAN_SOURCES})

    addCompilerFlagsSwan(templib_swan1serial ${ADDITIONAL_FLAGS_SWAN})
    addCompilerFlagsSwan(templib_swan2serial ${ADDITIONAL_FLAGS_SWAN})
    addCompilerFlags(templib_adcswan1 ${ADDITIONAL_FLAGS_ADCIRC})
    addCompilerFlags(adcswan ${ADDITIONAL_FLAGS_ADCIRC})

    TARGET_COMPILE_DEFINITIONS(templib_adcswan1 PRIVATE CSWAN)
    TARGET_COMPILE_DEFINITIONS(adcswan          PRIVATE CSWAN)

    TARGET_INCLUDE_DIRECTORIES(templib_adcswan1    PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1serial)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2serial PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1serial)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2serial PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_adcswan1)
    TARGET_INCLUDE_DIRECTORIES(adcswan             PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_adcswan1)   
    TARGET_INCLUDE_DIRECTORIES(adcswan             PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1serial)
    TARGET_INCLUDE_DIRECTORIES(adcswan             PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan2serial)

    TARGET_LINK_LIBRARIES(adcswan templib_adcswan1 templib_swan2serial templib_swan1serial )
    
    ADD_DEPENDENCIES(adcswan             templib_adcswan1 templib_swan2serial templib_swan1serial )
    ADD_DEPENDENCIES(templib_swan2serial templib_adcswan1 templib_swan1serial )
    ADD_DEPENDENCIES(templib_adcswan1    templib_swan1serial)
    
    INSTALL(TARGETS adcswan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_ADCSWAN AND PERL_FOUND)
