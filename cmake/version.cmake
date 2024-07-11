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
add_library(version ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F)
if(WIN32)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/version.F)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/version.F
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/version_default.F)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/version_default.F
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/version.F)
else(WIN32)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
    COMMAND ./adcirc_version.py --create-version-file --directory ${CMAKE_CURRENT_SOURCE_DIR} >/dev/null
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/version.F ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    COMMENT "Generating ADCIRC version...")
  add_custom_target(
    version_generate
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F)
endif(WIN32)
set_target_properties(version PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/version_mod)
set_target_properties(version PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
set_target_properties(version PROPERTIES EXCLUDE_FROM_ALL TRUE)
add_dependencies(version version_generate)
