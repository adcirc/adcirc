set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR "x86_64")

# Compiler Paths
set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries_2020.4.304/linux/bin/intel64/icc)
set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2020.4.304/linux/bin/intel64/icpc)
set(CMAKE_Fortran_COMPILER /opt/intel/compilers_and_libraries_2020.4.304/linux/bin/intel64/ifort)
set(MPI_C_COMPILER /opt/packages/openmpi/intel/4.0.2-intel20.4/bin/mpicc)
set(MPI_CXX_COMPILER /opt/packages/openmpi/intel/4.0.2-intel20.4/bin/mpicxx)
set(MPI_Fortran_COMPILER /opt/packages/openmpi/intel/4.0.2-intel20.4/bin/mpif90)

# Compiler Flags
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3 -mavx2 -funroll-loops")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-DNDEBUG -g -O2 -mavx2")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -mavx2 -funroll-loops")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -g -O2 -mavx2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mavx2 -funroll-loops")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -traceback -mavx2")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -check bounds")

# Library Paths SET(ENABLE_OUTPUT_NETCDF ON) SET(NETCDFHOME /path/to/netcdf-root) SET(NETCDF_F90_ROOT
# /path/to/netcdf-fortran-root)

# Adjust default behaviors
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
