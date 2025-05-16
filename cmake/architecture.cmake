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
message(STATUS "System Architecture Detected: ${CMAKE_SYSTEM_PROCESSOR}")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mpi_check.cmake)

# ######################################################################################################################
# ...Compiler specific options
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set(Fortran_LINELENGTH_FLAG
      "-ffixed-line-length-none"
      CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(Fortran_COMPILER_SPECIFIC_FLAG
        "-mcmodel=medium"
        CACHE STRING "Compiler specific flags")
  else()
    set(Fortran_COMPILER_SPECIFIC_FLAG
        ""
        CACHE STRING "Compiler specific flag")
  endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
  # ifort/ifx
  set(Fortran_LINELENGTH_FLAG
      "-132"
      CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # Heap array allocation
  execute_process(COMMAND sh -c "ulimit -s" OUTPUT_VARIABLE STACKSIZE)
  string(STRIP "${STACKSIZE}" STACKSIZE_TRIMMED)
  if("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(STATUS "Unlimited stack size detected. No heap array flag added.")
    set(heaparray_FLAG " ")
  else("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(
      STATUS
        "The compiler flag -heap-arrays ${STACKSIZE_TRIMMED} is being added. This should also be used to compile netcdf-fortran."
    )
    set(heaparray_FLAG "-heap-arrays ${STACKSIZE_TRIMMED}")
  endif("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")

  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl -mcmodel=medium")
  else()
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl")
  endif()

  string(STRIP ${ifort_FLAG} ifort_FLAG_TRIMMED)
  set(Fortran_COMPILER_SPECIFIC_FLAG
      ${ifort_FLAG_TRIMMED}
      CACHE STRING "Compiler specific flags")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  # pgf90
  set(Fortran_LINELENGTH_FLAG
      "-Mextend"
      CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(Fortran_COMPILER_SPECIFIC_FLAG
        "-Mlarge_arrays"
        CACHE STRING "Compiler specific flags")
  endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang")  
  # flang
  set(Fortran_LINELENGTH_FLAG 
    "-ffixed-line-length=132"
    CACHE STRING "Compiler specific flag to enable extended Fortran line length")

else()
  message(WARNING "Unknown Fortran Compiler. Fortran Compiler ID detected as ${CMAKE_Fortran_COMPILER_ID}")
  message(
    WARNING
      "No known predefined Fortran extended line length flag known. Please manually set the Fortran_LINELENGTH_FLAG")
  set(Fortran_LINELENGTH_FLAG
      ""
      CACHE STRING "Compiler specific flag to enable extended Fortran line length")
  set(Fortran_COMPILER_SPECIFIC_FLAG
      ""
      CACHE STRING "Compiler specific flags")
endif()
