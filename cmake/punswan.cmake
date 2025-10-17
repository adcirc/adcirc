# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# ######################################################################################################################
if(BUILD_PUNSWAN AND PERL_FOUND)

  set(SWANONLY1_PARALLEL_SOURCES
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swmod2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSpectPart.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_constants.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/m_fileio.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/serv_xnl4v5.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/mod_xnl4v5.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGriddata.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridobjects.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompdata.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/couple2adcirc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanIEM.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBraggScat.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanQCM.f90)

  set(SWANONLY2_PARALLEL_SOURCES
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanmain.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanpre2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom3.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom4.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swancom5.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanout2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanser.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/swanparll.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadADCGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadTriangleGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadEasymeshGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInitCompGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCheckGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCreateEdges.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridTopology.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridVert.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridCell.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGridFace.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrintGridInfo.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindPoint.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPointinMesh.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBpntlist.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPrepComp.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVertlist.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCompUnstruc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDispParm.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelX.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSweepSel.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPropvelS.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspAc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanTranspX.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanDiffPar.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGSECorr.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGradDepthorK.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolatePoint.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateAc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanInterpolateOutput.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvAccur.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanConvStopc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanFindObstacles.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanCrossObstacle.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanComputeForce.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanIntgratSpc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanBndStruc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanReadfort18.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanPunCollect.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanSumOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMinOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanMaxOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpids.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpcre.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/ocpmix.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SdsBabanin.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanGradVel.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVTKWriteHeader.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVTKWriteData.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/SwanVTKPDataSets.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/fftpack51.f90)

  set(MSGLIB_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/messenger.F90)

  add_library(templib_punmsglib OBJECT ${MSGLIB_SOURCES})
  add_library(templib_punswan1 OBJECT ${SWANONLY1_PARALLEL_SOURCES})
  add_executable(punswan ${SWANONLY2_PARALLEL_SOURCES})

  # Set the linker language to Fortran for the executable
  set_target_properties(punswan PROPERTIES LINKER_LANGUAGE Fortran)

  # ...SWAN Configuration
  adcirc_swan_configure_parallel()

  adcirc_add_compiler_flags(templib_punmsglib)
  adcirc_add_compiler_flags_swan(templib_punswan1 ${ADDITIONAL_FLAGS_SWAN})
  adcirc_add_compiler_flags_swan(punswan ${ADDITIONAL_FLAGS_SWAN})
  adcirc_add_mpi(templib_punmsglib)
  adcirc_add_mpi(templib_punswan1)
  adcirc_add_mpi(punswan)
  adcirc_add_datetime_definitions(templib_punmsglib)
  adcirc_add_datetime_definitions(templib_punswan1)
  adcirc_add_datetime_definitions(punswan)
  adcirc_add_datetime_libraries(punswan)

  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source)

  target_include_directories(templib_punswan1 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/templib_punmsglib)
  target_include_directories(punswan PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/templib_punswan1)
  target_include_directories(punswan PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/templib_punmsglib)

  target_link_libraries(
    punswan
    templib_punmsglib
    templib_punswan1
    mkdir
    ${MPI_Fortran_LIBRARIES})

  add_dependencies(templib_punmsglib templib_punswan1)
  add_dependencies(punswan templib_punswan1)
  add_dependencies(punswan templib_punmsglib)
  add_dependencies(templib_punmsglib version)
  add_dependencies(templib_punmsglib mkdir)

  install(TARGETS punswan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # Conditionally enable strict compiler flags for developers
  enable_developer_mode(${MSGLIB_SOURCES})

endif(BUILD_PUNSWAN AND PERL_FOUND)
