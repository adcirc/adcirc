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

if(WIN32)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/version.F)
    add_custom_target(
      generate_version_file
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/version.F
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/version_default.F)
    add_custom_target(
      generate_version_file
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/version_default.F
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  endif()
else()
  add_custom_target(
    generate_version_file
    COMMAND ./adcirc_version.py --create-version-file --directory ${CMAKE_CURRENT_SOURCE_DIR} >/dev/null
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/version.F
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    COMMENT "Generating ADCIRC version...")
endif()

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F PROPERTIES GENERATED TRUE)

add_library(adcirc_version OBJECT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F)
add_library(adcirc::version ALIAS adcirc_version)

add_dependencies(adcirc_version generate_version_file)

adcirc_set_module_directory(adcirc_version)
set_target_properties(adcirc_version PROPERTIES EXCLUDE_FROM_ALL TRUE)

if(Fortran_LINELENGTH_FLAG OR Fortran_COMPILER_SPECIFIC_FLAG)
  set(_version_flags "")
  if(Fortran_LINELENGTH_FLAG)
    list(APPEND _version_flags ${Fortran_LINELENGTH_FLAG})
  endif()
  if(Fortran_COMPILER_SPECIFIC_FLAG)
    list(APPEND _version_flags ${Fortran_COMPILER_SPECIFIC_FLAG})
  endif()
  target_compile_options(adcirc_version PRIVATE ${_version_flags})
endif()
