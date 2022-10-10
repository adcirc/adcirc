if(BUILD_PUNSWAN AND PERL_FOUND)

  set(SWANONLY1_PARALLEL_SOURCES
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSpectPart.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_constants.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_fileio.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/serv_xnl4v5.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/mod_xnl4v5.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGriddata.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridobjects.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompdata.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/couple2adcirc.f90)

  set(SWANONLY2_PARALLEL_SOURCES
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanmain.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom3.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom4.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom5.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanser.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanparll.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadADCGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadTriangleGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadEasymeshGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInitCompGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCheckGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCreateEdges.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridTopology.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridVert.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridCell.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridFace.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrintGridInfo.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindPoint.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPointinMesh.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBpntlist.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrepComp.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVertlist.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompUnstruc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDispParm.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelX.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSweepSel.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelS.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspAc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspX.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDiffPar.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGSECorr.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGradDepthorK.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolatePoint.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateAc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateOutput.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvAccur.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvStopc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindObstacles.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCrossObstacle.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanComputeForce.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanIntgratSpc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBndStruc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadfort18.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPunCollect.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSumOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMinOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMaxOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpids.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpcre.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpmix.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SdsBabanin.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGradVel.f90)

  set(MSGLIB_SOURCES
      ${CMAKE_SOURCE_DIR}/src/sizes.F
      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_SOURCE_DIR}/src/global.F
      ${CMAKE_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_SOURCE_DIR}/src/messenger.F)

  add_library(templib_punmsglib ${MSGLIB_SOURCES})
  add_library(templib_punswan1 ${SWANONLY1_PARALLEL_SOURCES})
  add_executable(punswan ${SWANONLY2_PARALLEL_SOURCES})

  # ...SWAN Configuration
  swanconfigureparallel()

  addcompilerflags(templib_punmsglib)
  addcompilerflagsswan(templib_punswan1 ${ADDITIONAL_FLAGS_SWAN})
  addcompilerflagsswan(punswan ${ADDITIONAL_FLAGS_SWAN})

  addlibversion(templib_punmsglib)
  addlibversion(punswan)
  addmpi(templib_punmsglib)
  addmpi(templib_punswan1)
  addmpi(punswan)
  addnetcdf(punswan)
  addxdmf(punswan)

  set_directory_properties(
    PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
               ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source)

  target_include_directories(
    templib_punswan1
    PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_punmsglib)
  target_include_directories(
    punswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_punswan1)
  target_include_directories(
    punswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_punmsglib)

  target_link_libraries(punswan templib_punmsglib templib_punswan1 mkdir
                        ${MPI_Fortran_LIBRARIES})

  add_dependencies(templib_punmsglib templib_punswan1)
  add_dependencies(punswan templib_punswan1)
  add_dependencies(punswan templib_punmsglib)
  add_dependencies(templib_punmsglib version)
  add_dependencies(templib_punmsglib mkdir)

  install(TARGETS punswan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_PUNSWAN AND PERL_FOUND)
