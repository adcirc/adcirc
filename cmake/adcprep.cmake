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
if(BUILD_ADCPREP)

  set(ADCPREP_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/vew1d.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/rs2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/nws08.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/normal_flow_boundary.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/presizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/pre_global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/metis.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/subprep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/adcprep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/decomp.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/prep_weir.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/itpackv.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/harm.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/read_global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subdomain.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gwce.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wetdry.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/momentum.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/prep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/interp.F
      ${CMAKE_CURRENT_SOURCE_DIR}/prep/machdep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sponge_layer.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/quadrature.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2loc_mapping.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/internaltide.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subgridLookup.F)

  if(NETCDF_WORKING)
    set(ADCPREP_SOURCES
        ${ADCPREP_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind_netcdf.F
        ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfio.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdf_error.F90)
  endif()

  add_executable(adcprep ${ADCPREP_SOURCES})

  adcirc_add_compiler_flags(adcprep ${ADDITIONAL_FLAGS_ADCPREP})
  adcirc_add_libraries(adcprep)
  adcirc_add_metis_library(adcprep)

  if(BUILD_PADCSWAN OR BUILD_PUNSWAN)
    target_compile_definitions(adcprep PRIVATE ${PREP_SWAN_FLAG})
  endif(BUILD_PADCSWAN OR BUILD_PUNSWAN)

  target_include_directories(adcprep PRIVATE prep)

  install(TARGETS adcprep RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # Conditionally enable strict compiler flags for developers
  enable_developer_mode(${ADCPREP_SOURCES})

endif(BUILD_ADCPREP)
