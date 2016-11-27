
IF(BUILD_SWAN AND PERL_FOUND)
    
    SET( SWANONLY_SERIAL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swmod1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/couple2adcirc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanmain.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanpre2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom4.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanout1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanser.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpids.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpmix.f )
    
    ADD_EXECUTABLE(swan ${SWANONLY_SERIAL_SOURCES})

    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES} 
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )
    ENDIF(WIN32)
    
    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source)

    SET_TARGET_PROPERTIES(swan PROPERTIES COMPILE_FLAGS "${Fortran_COMPILER_SPECIFIC_FLAG} 
                                                         ${ADDITIONAL_FLAGS_SWAN}")
    SET_TARGET_PROPERTIES(swan PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/swan_mod)
    
    TARGET_INCLUDE_DIRECTORIES(swan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/swan_mod)
    
    INSTALL(TARGETS swan RUNTIME DESTINATION bin)

ENDIF(BUILD_SWAN AND PERL_FOUND)
