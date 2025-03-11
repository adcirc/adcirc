# CMake Toolchain file for Apple Silicon (arm64) using Homebrew
SET(CMAKE_SYSTEM_NAME Darwin)
SET(CMAKE_SYSTEM_PROCESSOR arm64)

# Set the build type (Release, Debug, RelWithDebInfo, MinSizeRel)
SET(CMAKE_BUILD_TYPE RelWithDebInfo)

# Compiler Paths
SET(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-14)
SET(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-14)
SET(CMAKE_Fortran_COMPILER /opt/homebrew/bin/gfortran-14)
SET(MPI_C_COMPILER /opt/homebrew/bin/mpicc)
SET(MPI_CXX_COMPILER /opt/homebrew/bin/mpicxx)
SET(MPI_Fortran_COMPILER /opt/homebrew/bin/mpif90)

# Compiler Flags for the four build types
SET(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3 -funroll-loops -mtune=native")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-DNDEBUG -g -O2 -funroll-loops -mtune=native")
SET(CMAKE_C_FLAGS_MINSIZEREL "-DNDEBUG -Os -mtune=native")
SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")

SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -funroll-loops -mtune=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -g -O2 -funroll-loops -mtune=native")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "-DNDEBUG -Os -mtune=native")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops -mtune=native")
SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -fbacktrace -funroll-loops -mtune=native")
SET(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -mtune=native")
SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fbounds-check")

# Options
SET(ADCIRC_ENABLE_OUTPUT_NETCDF ON)
SET(ADCIRC_ENABLE_GRIB2 ON)
SET(ADCIRC_ENABLE_DATETIME ON)

# Library Paths
SET(NETCDFHOME /opt/homebrew)
SET(NETCDF_F90_ROOT /opt/homebrew)

# Executables
SET(ADCIRC_BUILD_ADCIRC ON)
SET(ADCIRC_BUILD_PADCIRC ON)
SET(ADCIRC_BUILD_ADCPREP ON)
SET(ADCIRC_BUILD_PADCSWAN ON)
SET(ADCIRC_BUILD_ADCSWAN ON)
SET(ADCIRC_BUILD_ASWIP ON)

# Adjust default behaviors
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)