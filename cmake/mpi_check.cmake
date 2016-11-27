
FIND_PACKAGE(MPI)
IF(MPI_FOUND)
    #...Check if we have mpi.mod for Fortran or need
    #   to use mpif.h, which is discouraged
    SET(mpi_f90mod_check 
    "       PROGRAM MPIMOD_CHECK
            USE MPI
            IMPLICIT NONE
            INTEGER :: IERR
            CALL MPI_INIT(IERR)
            END PROGRAM
    ")
    FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90" "${mpi_f90mod_check}")
    TRY_COMPILE(MPIMOD_COMPILE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpif90_mod_check.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}" "-DLINK_LIBRARIES=${MPI_Fortran_LIBRARIES}" OUTPUT_VARIABLE LOG)
    IF(MPIMOD_COMPILE)
        SET(MPIMOD_FLAG "-DHAVE_MPI_MOD")
    ELSE(MPIMOD_COMPILE)
        SET(MPI_FOUND FALSE)
        MESSAGE(WARNING "The MPI library specified does not function with the specified compilers. Parallel ADCIRC/SWAN compilation will not be enabled until this is corrected.")
    ENDIF(MPIMOD_COMPILE)
ENDIF(MPI_FOUND)
