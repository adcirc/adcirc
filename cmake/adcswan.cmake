
IF(BUILD_ADCSWAN AND PERL_FOUND)

    #...Note that we need to build ADCSWAN in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    ADD_EXECUTABLE(adcswan src/driver.F)
   
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
                             ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/ocpmix.f )
    
    SET( ADCSWAN1_SOURCES    src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F src/mesh.F 
                             src/global_3dvs.F src/harm.F wind/vortex.F src/wind.F src/owiwind.F src/rs2.F
                             src/owi_ice.F src/itpackv.F src/nodalattr.F src/globalio.F src/netcdfio.F 
                             src/control.F src/xdmfio.F )

    SET( ADCSWAN2_SOURCES    src/write_output.F src/couple2swan.F )

    SET( ADCSWAN3_SOURCES    src/adcirc.F src/subdomain.F src/weir_boundary.F src/read_input.F src/cstart.F 
                             src/hstart.F src/timestep.F src/vsmy.F src/transport.F )
    
    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ENDIF(WIN32)

    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source)

    ADD_LIBRARY(templib_swan1serial STATIC ${SWAN1SERIAL_SOURCES})
    ADD_LIBRARY(templib_swan2serial STATIC ${SWAN2SERIAL_SOURCES})
    ADD_LIBRARY(templib_adcswan1 STATIC ${ADCSWAN1_SOURCES})
    ADD_LIBRARY(templib_adcswan2 STATIC ${ADCSWAN2_SOURCES})
    ADD_LIBRARY(templib_adcswan3 STATIC ${ADCSWAN3_SOURCES})
   
    TARGET_LINK_LIBRARIES(adcswan templib_swan1serial)
    TARGET_LINK_LIBRARIES(adcswan templib_swan2serial)
    
    SET(ADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_OPTION_FLAGS}")
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${ADCSWAN_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(adcswan PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    SET_TARGET_PROPERTIES(templib_swan1serial PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_swan2serial PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan1    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan2    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan3    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan             PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan1    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_adcswan2    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_adcswan3    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(adcswan             PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(templib_swan1serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(templib_swan2serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
    
    TARGET_INCLUDE_DIRECTORIES(adcswan  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan1 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan2 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan3 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2serial PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)

    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version
                                      netcdf netcdff )
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version
                                      netcdf netcdff XdmfCore XdmfUtils Xdmf )
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version )
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    ADD_DEPENDENCIES(adcswan             templib_adcswan3)
    ADD_DEPENDENCIES(templib_adcswan3    templib_swan2serial)
    ADD_DEPENDENCIES(templib_swan2serial templib_adcswan2)
    ADD_DEPENDENCIES(templib_adcswan3    templib_adcswan2)
    ADD_DEPENDENCIES(templib_adcswan2    templib_adcswan1)
    ADD_DEPENDENCIES(templib_adcswan1    templib_swan1serial)
    ADD_DEPENDENCIES(templib_adcswan1    version)
    
    INSTALL(TARGETS adcswan RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCSWAN AND PERL_FOUND)
