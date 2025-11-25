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

include(CheckFortranCompilerFlag)

message(STATUS "System Architecture: ${CMAKE_SYSTEM_PROCESSOR}")

# Include MPI detection and interface libraries
include(${CMAKE_CURRENT_LIST_DIR}/mpi_check.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/interface_libraries.cmake)

# ######################################################################################################################
# Compiler-specific options
# ######################################################################################################################

# Initialize flag variables
set(Fortran_LINELENGTH_FLAG "")
set(Fortran_COMPILER_SPECIFIC_FLAG "")

# ######################################################################################################################
# Detect Fortran line length flag using feature detection
# ######################################################################################################################
# Test common Fortran line length flags in order of preference
set(LINE_LENGTH_FLAGS
    "-ffixed-line-length-none" # GNU (unlimited)
    "-132" # Intel
    "-Mextend" # PGI/NVHPC
    "-ffixed-line-length=132" # LLVM Flang
    "-extend-source" # Cray
    "-qfixed=132" # IBM XL
)

foreach(flag ${LINE_LENGTH_FLAGS})
  string(MAKE_C_IDENTIFIER "COMPILER_SUPPORTS${flag}" flag_var)
  CHECK_Fortran_COMPILER_FLAG("${flag}" ${flag_var})
  if(${flag_var})
    set(Fortran_LINELENGTH_FLAG "${flag}")
    break()
  endif()
endforeach()

if(NOT Fortran_LINELENGTH_FLAG)
  message(WARNING "Could not determine Fortran line length flag. Long fixed-form lines may fail to compile.")
endif()

# ######################################################################################################################
# Intel-specific compiler flags
# ######################################################################################################################
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(Fortran_COMPILER_SPECIFIC_FLAG "")
  list(
    APPEND
    Fortran_COMPILER_SPECIFIC_FLAG
    "-assume"
    "byterecl")

  if(UNIX)
    execute_process(
      COMMAND sh -c "ulimit -s"
      OUTPUT_VARIABLE STACKSIZE
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

    if(STACKSIZE STREQUAL "unlimited")
      message(STATUS "Unlimited stack size detected. No heap array flag added.")
    elseif(STACKSIZE)
      message(STATUS "Stack size: ${STACKSIZE}KB - adding heap-arrays flag")
      message(STATUS "Note: netcdf-fortran should also be compiled with -heap-arrays ${STACKSIZE}")
      list(
        APPEND
        Fortran_COMPILER_SPECIFIC_FLAG
        "-heap-arrays"
        "${STACKSIZE}")
    endif()
  endif()
endif()

# Display final flags
message(STATUS "Fortran line length flag: ${Fortran_LINELENGTH_FLAG}")
if(Fortran_COMPILER_SPECIFIC_FLAG)
  message(STATUS "Fortran compiler flags: ${Fortran_COMPILER_SPECIFIC_FLAG}")
endif()
