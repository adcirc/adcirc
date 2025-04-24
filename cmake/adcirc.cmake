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
if(BUILD_ADCIRC)

  set(ADCIRC_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/logging.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F
      ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/vew1d.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/harm.F
      ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/rs2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/itpackv.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/globalio.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subdomain.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gwce.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wetdry.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/momentum.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/control.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/write_output.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2swan.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/adcirc.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/weir_boundary.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/read_input.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/cstart.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hstart.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/timestep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/vsmy.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/transport.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/driver.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sponge_layer.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/quadrature.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2loc_mapping.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/internaltide.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/astronomic.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sun.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/moon.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sun_moon_system.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ephemerides.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/tidalpotential.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subgridLookup.F)

  if(NETCDF_WORKING)
    set(ADCIRC_SOURCES
        ${ADCIRC_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind_netcdf.F
        ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfio.F
        ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdf_error.F90)
  endif()

  if(XDMF_WORKING)
    set(ADCIRC_SOURCES ${ADCIRC_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/xdmfio.F)
  endif()

  add_executable(adcirc ${ADCIRC_SOURCES})
  set(ADCIRC_COMPILER_FLAGS "${ADDITIONAL_FLAGS_ADCIRC} ${ADCIRC_OPTION_FLAGS}")
  addcompilerflags(adcirc ${ADDITIONAL_FLAGS_ADCIRC})
  addnetcdflibraries(adcirc)
  addgrib2libraries(adcirc)
  adddatetimelibraries(adcirc)
  addxdmflibraries(adcirc)
  addversionlibrary(adcirc)
  addmkdirlibrary(adcirc)
  install(TARGETS adcirc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_ADCIRC)
