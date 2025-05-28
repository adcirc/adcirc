# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
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
if(BUILD_SWAN AND PERL_FOUND)

  set(SWANONLY_SERIAL_SOURCES
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swmod1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swmod2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSpectPart.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/m_constants.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/m_fileio.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/serv_xnl4v5.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/mod_xnl4v5.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGriddata.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridobjects.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCompdata.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/couple2adcirc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanmain.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanpre1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanpre2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom3.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom4.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swancom5.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanout1.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanout2.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanser.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/swanparll.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadADCGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadTriangleGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadEasymeshGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInitCompGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCheckGrid.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCreateEdges.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridTopology.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridVert.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridCell.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGridFace.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPrintGridInfo.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanFindPoint.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPointinMesh.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanBpntlist.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPrepComp.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanVertlist.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCompUnstruc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanDispParm.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPropvelX.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSweepSel.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPropvelS.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanTranspAc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanTranspX.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanDiffPar.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGSECorr.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGradDepthorK.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolatePoint.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolateAc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanInterpolateOutput.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanConvAccur.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanConvStopc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanFindObstacles.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanCrossObstacle.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanComputeForce.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanIntgratSpc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanBndStruc.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanReadfort18.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanPunCollect.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanSumOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanMinOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanMaxOverNodes.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpids.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpcre.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/ocpmix.f
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SdsBabanin.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanGradVel.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanIEM.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanBraggScat.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanQCM.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanVTKWriteHeader.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanVTKWriteData.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/SwanVTKPDataSets.f90
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/fftpack51.f90)

  add_executable(swan ${SWANONLY_SERIAL_SOURCES})

  # ...SWAN Configuration
  adcirc_swan_configure_serial()

  adcirc_add_compiler_flags_swan(swan ${ADDITIONAL_FLAGS_SWAN})

  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source)

  install(TARGETS swan RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_SWAN AND PERL_FOUND)
