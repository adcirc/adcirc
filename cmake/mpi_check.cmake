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
find_package(MPI)
if(MPI_FOUND)
  # ...Check if we have mpi.mod for Fortran or need to use mpif.h, which is discouraged
  set(mpi_f90mod_check
      "       PROGRAM MPIMOD_CHECK
            USE MPI
            IMPLICIT NONE
            INTEGER :: IERR
            CALL MPI_INIT(IERR)
            END PROGRAM
    ")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90" "${mpi_f90mod_check}")
  try_compile(
    MPI_COMPILE "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}" "-DLINK_LIBRARIES=${MPI_Fortran_LIBRARIES}"
    OUTPUT_VARIABLE LOG)
  if(NOT MPI_COMPILE)
    set(MPI_FOUND FALSE)
    message(
      WARNING
        "The MPI library specified does not function with the specified compilers. Parallel ADCIRC/SWAN compilation will not be enabled until this is corrected."
    )
  endif(NOT MPI_COMPILE)
endif(MPI_FOUND)
