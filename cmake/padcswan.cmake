if(BUILD_PADCSWAN AND PERL_FOUND)

  set(SWAN1PARALLEL_SOURCES
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swmod1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swmod2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSpectPart.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/m_constants.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/m_fileio.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/serv_xnl4v5.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/mod_xnl4v5.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGriddata.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridobjects.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCompdata.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/couple2adcirc.f90)

  set(SWAN2PARALLEL_SOURCES
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanmain.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanpre1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanpre2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom3.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom4.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swancom5.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanout1.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanout2.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanser.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/swanparll.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadADCGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadTriangleGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadEasymeshGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInitCompGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCheckGrid.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCreateEdges.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridTopology.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridVert.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridCell.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGridFace.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPrintGridInfo.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanFindPoint.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPointinMesh.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanBpntlist.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPrepComp.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanVertlist.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCompUnstruc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanDispParm.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPropvelX.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSweepSel.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPropvelS.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanTranspAc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanTranspX.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanDiffPar.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGSECorr.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGradDepthorK.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolatePoint.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolateAc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanInterpolateOutput.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanConvAccur.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanConvStopc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanFindObstacles.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanCrossObstacle.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanComputeForce.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanIntgratSpc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanBndStruc.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanReadfort18.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanPunCollect.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanSumOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanMinOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanMaxOverNodes.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpids.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpcre.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/ocpmix.f
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SdsBabanin.f90
      ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/SwanGradVel.f90)

  set(PADCSWAN1_SOURCES
      ${CMAKE_SOURCE_DIR}/src/sizes.F
      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_SOURCE_DIR}/src/global.F
      ${CMAKE_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_SOURCE_DIR}/src/messenger.F
      ${CMAKE_SOURCE_DIR}/src/mesh.F
      ${CMAKE_SOURCE_DIR}/src/harm.F
      ${CMAKE_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_SOURCE_DIR}/src/wind.F
      ${CMAKE_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_SOURCE_DIR}/src/owiwind_netcdf.F
      ${CMAKE_SOURCE_DIR}/src/rs2.F
      ${CMAKE_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_SOURCE_DIR}/src/itpackv.F
      ${CMAKE_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_SOURCE_DIR}/src/globalio.F
      ${CMAKE_SOURCE_DIR}/src/write_output.F
      ${CMAKE_SOURCE_DIR}/src/writer.F
      ${CMAKE_SOURCE_DIR}/src/couple2swan.F
      ${CMAKE_SOURCE_DIR}/src/netcdfio.F
      ${CMAKE_SOURCE_DIR}/src/subdomain.F
      ${CMAKE_SOURCE_DIR}/src/sponge_layer.F
      ${CMAKE_SOURCE_DIR}/src/quadrature.F
      ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F
      ${CMAKE_SOURCE_DIR}/src/wetdry.F
      ${CMAKE_SOURCE_DIR}/src/gwce.F
      ${CMAKE_SOURCE_DIR}/src/momentum.F
      ${CMAKE_SOURCE_DIR}/src/xdmfio.F
      ${CMAKE_SOURCE_DIR}/src/control.F)

  set(PADCSWAN_SOURCES
      ${CMAKE_SOURCE_DIR}/src/couple2swan.F
      ${CMAKE_SOURCE_DIR}/src/adcirc.F
      src/subdomain.F
      ${CMAKE_SOURCE_DIR}/src/weir_boundary.F
      ${CMAKE_SOURCE_DIR}/src/read_input.F
      ${CMAKE_SOURCE_DIR}/src/cstart.F
      ${CMAKE_SOURCE_DIR}/src/hstart.F
      ${CMAKE_SOURCE_DIR}/src/timestep.F
      ${CMAKE_SOURCE_DIR}/src/vsmy.F
      ${CMAKE_SOURCE_DIR}/src/transport.F
      ${CMAKE_SOURCE_DIR}/src/driver.F)

  # ...SWAN Configuration
  swanconfigurepadcswan()

  set_directory_properties(
    PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
               ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source)

  add_library(templib_swan1parallel STATIC ${SWAN1PARALLEL_SOURCES})
  add_library(templib_swan2parallel STATIC ${SWAN2PARALLEL_SOURCES})
  add_library(templib_padcswan1 STATIC ${PADCSWAN1_SOURCES})
  add_executable(padcswan ${PADCSWAN_SOURCES})

  addcompilerflags(templib_padcswan1)
  addcompilerflags(padcswan)
  addcompilerflagsswan(templib_swan1parallel ${ADDITIONAL_FLAGS_SWAN})
  addcompilerflagsswan(templib_swan2parallel ${ADDITIONAL_FLAGS_SWAN})

  addlibmkdir(templib_padcswan1)
  addlibmkdir(padcswan)

  addlibversion(padcswan)

  addmpi(templib_padcswan1)
  addmpi(padcswan)
  addmpi(templib_swan1parallel)
  addmpi(templib_swan2parallel)

  target_compile_definitions(templib_padcswan1 PRIVATE CSWAN)
  target_compile_definitions(padcswan PRIVATE CSWAN)

  target_include_directories(
    templib_padcswan1
    PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
  target_include_directories(
    templib_swan2parallel
    PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcswan1)
  target_include_directories(
    templib_swan2parallel
    PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
  target_include_directories(
    padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcswan1)
  target_include_directories(
    padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan1parallel)
  target_include_directories(
    padcswan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_swan2parallel)

  target_link_libraries(padcswan templib_swan2parallel templib_padcswan1
                        templib_swan1parallel)

  add_dependencies(padcswan templib_swan2parallel templib_padcswan1)
  add_dependencies(templib_swan2parallel templib_padcswan1
                   templib_swan1parallel)
  add_dependencies(templib_padcswan1 mkdir version templib_swan1parallel)

  install(TARGETS padcswan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_PADCSWAN AND PERL_FOUND)
