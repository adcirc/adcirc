# Cmake version, project name, language
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(fortran-metis-interface Fortran)

# Set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release
# If not specified default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

# Handle integer size
if (INT)
    message("Configuring build for ${INT}-bit integers")
    add_definitions(-DINT${INT})
else()
    message("Configuring build for 32-bit integers")
    add_definitions(-DINT32)
endif()

# Handle real size
if (REAL)
    message("Configuring build for ${INT}-bit reals")
    add_definitions(-DREAL${REAL})
else()
    message("Configuring build for 32-bit reals")
    add_definitions(-DREAL32)
endif()


# Compiler flags for gfortran.
if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(common "-Wall -cpp")
    set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -O0 -g -C -fbacktrace")
endif()

# Compiler flags for ifort.
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(common "-warn all -fpp")
    set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -O0 -g -C -traceback")
endif ()

# Find METIS library
if (NOT METIS_LIB)
    find_library(METIS_LIB metis)
endif ()

# Output messages
message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")


# Library to archive (libfmetis.a)
add_library(fmetis src/lib/metis_interface.f90)
target_link_libraries(fmetis ${METIS_LIB})

install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")
install(DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")
