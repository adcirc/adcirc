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
if(BUILD_ASWIP)

  set(ASWIP_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/internaltide.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2loc_mapping.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/grid_search.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws08.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subgridLookup.F
      ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_CURRENT_SOURCE_DIR}/wind/aswip.F)

  if(NETCDF_WORKING)
    set(ASWIP_SOURCES
        ${ASWIP_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws13.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws13data.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdf_error.F90)
  endif()

  add_executable(aswip ${ASWIP_SOURCES})

  adcirc_add_compiler_flags(aswip ${ADDITIONAL_FLAGS_ASWIP})
  adcirc_add_libraries(aswip)

  install(TARGETS aswip RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # Conditionally enable strict compiler flags for developers
  enable_developer_mode(${ASWIP_SOURCES})

  # Set the linker language to Fortran for the executable
  set_target_properties(aswip PROPERTIES LINKER_LANGUAGE Fortran)

endif(BUILD_ASWIP)
