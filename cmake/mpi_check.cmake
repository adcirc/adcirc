find_package(MPI)
if (MPI_FOUND)
  # ...Check if we have mpi.mod for Fortran or need to use mpif.h, which is
  # discouraged
  set(mpi_f90mod_check
          "       PROGRAM MPIMOD_CHECK
            USE MPI
            IMPLICIT NONE
            INTEGER :: IERR
            CALL MPI_INIT(IERR)
            END PROGRAM
    ")
  set(mpi_f90inc_check
          "       PROGRAM MPIMOD_CHECK
            INCLUDE 'mpif.h'
            IMPLICIT NONE
            INTEGER :: IERR
            CALL MPI_INIT(IERR)
            END PROGRAM
    ")
  file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90"
          "${mpi_f90mod_check}")
  file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/mpif90_inc_check.f90"
          "${mpi_f90inc_check}")
  try_compile(
          MPIMOD_COMPILE "${CMAKE_CURRENT_BINARY_DIR}"
          "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90"
          CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}"
          "-DLINK_LIBRARIES=${MPI_Fortran_LIBRARIES}"
          OUTPUT_VARIABLE LOG)
  if (MPIMOD_COMPILE)
    set(MPIMOD_FLAG "-DHAVE_MPI_MOD")
  else (MPIMOD_COMPILE)
    try_compile(
            MPIINC_COMPILE "${CMAKE_CURRENT_BINARY_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpif90_inc_check.f90"
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}"
            "-DLINK_LIBRARIES=${MPI_Fortran_LIBRARIES}"
            OUTPUT_VARIABLE LOG)
    if (MPIINC_COMPILE)
      set(MPIMOD_FLAG "")
    else (MPIINC_COMPILE)
      set(MPI_FOUND FALSE)
      message(
              WARNING
              "The MPI library specified does not function with the specified compilers. Parallel ADCIRC/SWAN compilation will not be enabled until this is corrected."
      )
    endif (MPIINC_COMPILE)
  endif (MPIMOD_COMPILE)
endif (MPI_FOUND)
