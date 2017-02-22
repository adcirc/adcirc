###########################################################################
#                 CMake Build File for ADCIRC(+SWAN)
#
#    Written By: Zach Cobell
#
###########################################################################
#
# The CMake build system enable ADCIRC (and SWAN) to be deployed and built
# in a cross platform environment. 
#
###########################################################################

#...Set the default build type
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
        build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
        Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
        options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
        RelWithDebInfo MinSizeRel.")
ENDIF()

###########################################################################
#  GENERAL OPTIONS
###########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(adcirc_model)

###########################################################################


###########################################################################
#  COMPILERS/LANGUAGE
###########################################################################


###########################################################################
#...Perl
#   Perl is required to use SWAN since it generates the source files. 
#   All SWAN options will be disabled if perl cannot be found. The user
#   is given a chance to specify its location
FIND_PROGRAM(PERL perl)
IF( ${PERL} STREQUAL "PERL-NOTFOUND")
    SET(PERL_FOUND FALSE)
ELSE()
    SET(PERL_FOUND TRUE)
ENDIF()
###########################################################################



###########################################################################
#...Language Specifications
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
###########################################################################


###########################################################################
#...Put the static libraries in the CMakeFiles folder so they don't
#   contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
###########################################################################

#...Include Macros to build components
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/macros.cmake)

#...Determine architecture specific parameters
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/architecture.cmake)

#...Generate the selection list of options
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/options_select.cmake)

#...Translate the selected options list into compiler flags
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/options_flags.cmake)

#...Ensure that the netCDF libraries are working
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/netcdf_check.cmake)

#...Ensure that the XDMF libraries are working
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/xdmf_check.cmake)

#...Build the version module
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/version.cmake)

#...Build the mkdir library
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/mkdir.cmake)

#...Build ADCIRC
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/adcirc.cmake)

#...Build METIS
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/metis.cmake)

#...Build ADCPREP
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/adcprep.cmake)

#...Build ADCSWAN (Serial ADCIRC+SWAN)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/adcswan.cmake)

#...Build SWAN (Serial structured SWAN)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/swan.cmake)

#...Build PADCIRC (Parallel ADCIRC)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/padcirc.cmake)

#...Build PADCSWAN (Parallel ADCIRC+SWAN)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/padcswan.cmake)

#...Build PUNSWAN (Parallel Unstructured SWAN)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/punswan.cmake)

#...Build ASWIP
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/aswip.cmake)

#...Build ADCIRC Utility Codes
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)
