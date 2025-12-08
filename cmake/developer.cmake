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
function(enable_developer_mode)
  if(ADCIRC_DEVELOPER_MODE)
    add_strict_compiler_flags(${ARGN})
  endif()
endfunction()

function(add_strict_compiler_flags)
  if(${CMAKE_Fortran_COMPILER_ID} MATCHES "IntelLLVM")
    set(STRICT_FLAGS "-warn all -diag-enable remark -implicit-none")
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set(STRICT_FLAGS
        "-Werror -Wall -Wextra -Wconversion -pedantic -fimplicit-none -Wuninitialized -Wsurprising -Wuse-without-only -Wimplicit-procedure -Winteger-division -Wconversion-extra"
    )
  else()
    message(WARNING "No developer compiler flags defined for ${CMAKE_Fortran_COMPILER_ID}. No action taken.")
  endif()

  foreach(SOURCE_FILE ${ARGN})
    if(${SOURCE_FILE} MATCHES ".*\\.F90$")
      set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS "${STRICT_FLAGS}")
    endif()
  endforeach()
endfunction()
