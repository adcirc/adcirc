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

find_package(MPI COMPONENTS Fortran)

if(MPI_FOUND)
  # Check if FindMPI successfully found the Fortran module interface If both F90 and F08 module interfaces are
  # unavailable, it means FindMPI fell back to the legacy mpif.h header interface, which often indicates a compiler
  # mismatch (e.g., gfortran-compiled MPI with Intel compiler). This won't work, so we bail out on building anything
  # requiring MPI.
  if(NOT MPI_Fortran_HAVE_F90_MODULE AND NOT MPI_Fortran_HAVE_F08_MODULE)
    set(MPI_FOUND FALSE)
    message(
      WARNING "MPI was detected but the Fortran module interface is not available.\n"
              "This typically indicates a mismatch between the MPI installation and Fortran compiler.\n"
              "For example, this occurs when using an MPI compiled with gfortran with the Intel compiler.\n"
              "MPI include path: ${MPI_Fortran_INCLUDE_DIRS}\n"
              "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}\n"
              "Parallel ADCIRC/SWAN builds (PADCIRC, PADCSWAN, etc.) will not be available.\n"
              "Serial builds (ADCIRC, ADCSWAN) will still work.\n"
              "To fix: Install or specify an MPI compiled with ${CMAKE_Fortran_COMPILER_ID}.")
  endif()
endif()

if((BUILD_PADCIRC
    OR BUILD_PADCSWAN
    OR BUILD_LIBADCIRC_STATIC
    OR BUILD_LIBADCIRC_SHARED
   )
   AND NOT MPI_FOUND)
  message(
    SEND_ERROR "PADCIRC, PADCSWAN, and parallel ADCIRC libraries require MPI, but MPI was not found or is unusable.\n"
               "Please install an MPI implementation compatible with your Fortran compiler and re-run CMake.")
endif()
