
IF(BUILD_PUNSWAN AND PERL_FOUND)
    
    SET( SWANONLY1_PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/couple2adcirc.f90 ) 
                               
    SET( SWANONLY2_PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanmain.f  ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom2.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom4.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout1.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanser.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpids.f ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpmix.f )

    SET( MSGLIB_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                         src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )
    
    ADD_LIBRARY(templib_punmsglib ${MSGLIB_SOURCES})
    ADD_LIBRARY(templib_punswan1  ${SWANONLY1_PARALLEL_SOURCES})
    ADD_EXECUTABLE(punswan        ${SWANONLY2_PARALLEL_SOURCES})

    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel unSWAN Sources..."
        )
    ENDIF(WIN32)

    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source)
    
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(punswan1       PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(punswan        PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF()
    
    SET_TARGET_PROPERTIES(templib_punswan1  PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    SET_TARGET_PROPERTIES(templib_punmsglib PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    SET_TARGET_PROPERTIES(punswan           PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    
    SET(MSGLIB_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    SET_TARGET_PROPERTIES(templib_punmsglib PROPERTIES COMPILE_FLAGS ${MSGLIB_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/punswan_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_punswan1  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_punswan1  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_punmsglib PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_punmsglib PRIVATE ${MPI_Fortran_INCLUDE_PATH})

    TARGET_LINK_LIBRARIES(punswan templib_punmsglib templib_punswan1 mkdir ${MPI_Fortran_LIBRARIES})

    ADD_DEPENDENCIES(templib_punmsglib templib_punswan1)
    ADD_DEPENDENCIES(punswan           templib_punswan1)
    ADD_DEPENDENCIES(punswan           templib_punmsglib)
    ADD_DEPENDENCIES(templib_punmsglib version)
    ADD_DEPENDENCIES(templib_punmsglib mkdir)
    
    INSTALL(TARGETS punswan RUNTIME DESTINATION bin)

ENDIF(BUILD_PUNSWAN AND PERL_FOUND)
