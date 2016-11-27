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


###########################################################################
#...MPI Compilers
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
###########################################################################



###########################################################################
#...Architecture Specifications. Determine if 
#   the system is 32 bit or 64 bit. If the 64 bit
#   integer pointer is not detected, an error is 
#   thrown during compile
SET(archdetect_c_code 
"#include <stdint.h>
int main()
{
#if INTPTR_MAX == INT64_MAX
        return 0;
#elif INTPTR_MAX == INT32_MAX
#error 32-bit max integer pointer
#else
#error Unknown pointer size or missing size macros
#endif
}"
)

FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/architecture_check.c" "${archdetect_c_code}")

TRY_COMPILE(ARCH_TEST_COMPILE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/architecture_check.c")
IF(ARCH_TEST_COMPILE)
    SET(ARCH 64)
ELSE(ARCH_TEST_COMPILE)
    SET(ARCH 32)
ENDIF(ARCH_TEST_COMPILE)
###########################################################################



###########################################################################
#...Compiler specific options
GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  SET(Fortran_LINELENGTH_FLAG "-ffixed-line-length-none" CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-mcmodel=medium" CACHE STRING "Compiler specific flags")  
  ENDIF(ARCH EQUAL 64)

ELSEIF(Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  SET(Fortran_LINELENGTH_FLAG "-132" CACHE STRING "Compiler specific flag to enable extended Fortran line length")
  
  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-mcmodel=medium" CACHE STRING "Compiler specific flags")  
  ENDIF(ARCH EQUAL 64)

ELSEIF(Fortran_COMPILER_NAME MATCHES "pgf90.*")
  # pgf90
  SET(Fortran_LINELENGTH_FLAG "-Mextend" CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-Mlarge_arrays" CACHE STRING "Compiler specific flags")  
  ENDIF(ARCH EQUAL 64)

ELSE(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  MESSAGE("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  MESSAGE("Fortran compiler: " ${Fortran_COMPILER_NAME})
  MESSAGE("No known predefined Fortran extended line length flag known. Please manually set the Fortran_LINELENGTH_FLAG")
  SET(Fortran_LINELENGTH_FLAG "" CACHE STRING "Compiler specific flag to enable extended Fortran line length")
  SET(Fortran_COMPILER_SPECIFIC_FLAG "" CACHE STRING "Compiler specific flags")
ENDIF(Fortran_COMPILER_NAME MATCHES "gfortran.*")
###########################################################################



###########################################################################
# OPTIONS
###########################################################################

#...Option for 8-byte real numbers (or 4-byte if turned off)
OPTION(PRECISION_8BYTE "Use 8-byte real numbers instead of 4-byte." ON)
MARK_AS_ADVANCED(PRECISION_8BYTE)


###########################################################################
#...Output format options
OPTION(ENABLE_OUTPUT_NETCDF "Turn on netCDF output options" OFF)
IF(ENABLE_OUTPUT_NETCDF)
    OPTION(ENABLE_OUTPUT_XDMF "Turn on XDMF output options" OFF)
ENDIF(ENABLE_OUTPUT_NETCDF)    
###########################################################################


###########################################################################
#...Executables
OPTION(BUILD_ADCIRC "Build the serial ADCIRC executable" OFF)

IF(PERL_FOUND)
    OPTION(BUILD_ADCSWAN  "Build the serial SWAN+ADCIRC executable" OFF)
    OPTION(BUILD_SWAN   "Build the serial SWAN executable" OFF)
ENDIF(PERL_FOUND)

IF(MPI_FOUND)
    OPTION(BUILD_ADCPREP  "Build the MPI parallel ADCIRC preprocessor" OFF)
    OPTION(BUILD_PADCIRC  "Build the MPI parallel ADCIRC executable" OFF)
    IF(PERL_FOUND)
        OPTION(BUILD_PADCSWAN "Build the MPI parallel SWAN+ADCIRC executable" OFF)
        OPTION(BUILD_PUNSWAN  "Build the MPI parallel unstructured SWAN executable" OFF)
    ENDIF(PERL_FOUND)
ENDIF(MPI_FOUND)

OPTION(BUILD_ASWIP "Build ASWIP (ASymmetric Wind Input Preprocessor)")
OPTION(BUILD_UTILITIES "Build the ADCIRC utility programs" OFF)
###########################################################################


###########################################################################
#...No need to show CXX compilers on main screen
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER CMAKE_C_COMPILER CMAKE_Fortran_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS_RELEASE CMAKE_Fortran_FLAGS_RELEASE)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_DEBUG CMAKE_C_FLAGS_DEBUG CMAKE_Fortran_FLAGS_DEBUG)
###########################################################################


###########################################################################
#...Library paths
IF(ENABLE_OUTPUT_NETCDF)
    SET(NETCDFHOME "NETCDF-NOTFOUND" CACHE STRING "netCDF home path containing lib and include")
ELSE(ENABLE_OUTPUT_NETCDF)
    UNSET(NETCDFHOME CACHE)
ENDIF(ENABLE_OUTPUT_NETCDF)

IF(ENABLE_OUTPUT_XDMF)
    SET(XDMFHOME "XDMF-NOTFOUND" CACHE STRING "XDMF home path containing lib and include")
ELSE(ENABLE_OUTPUT_XDMF)
    UNSET(XDMFHOME CACHE)
ENDIF(ENABLE_OUTPUT_XDMF)
###########################################################################


###########################################################################
#...Additional flags
SET(ADDITIONAL_FLAGS_ADCIRC "" CACHE STRING "Additional flags to compile ADCIRC with")
IF(PERL_FOUND)
    SET(ADDITIONAL_FLAGS_SWAN "" CACHE STRING "Additional flags to compile SWAN with")
ENDIF(PERL_FOUND)
SET(ADDITIONAL_FLAGS_ADCPREP "" CACHE STRING "Additional flags to compile ADCPREP with")
SET(ADDITIONAL_FLAGS_ASWIP "" CACHE STRING "Additional flags to compile ASWIP with")
SET(ADDITIONAL_FLAGS_UTLIITIES "" CACHE STRING "Additional flags for utility programs")
###########################################################################



###########################################################################
#...Options enabled via compiler flags within the code
OPTION(ENABLE_WARN_ELEV_DEBUG "Enable writing of the fort.69 debug file" OFF)

IF(BUILD_ADCSWAN OR BUILD_PADCSWAN)
    OPTION(ENABLE_SWAN_MODIFIED_FRICTION "Enable Ethan/Joannes' modified friction" OFF)
ENDIF(BUILD_ADCSWAN OR BUILD_PADCSWAN)

OPTION(IBM    "Format code for IBM based architectures"    OFF)
OPTION(SGI    "Format code for SGI based architectures"    OFF)
OPTION(SUN    "Format code for SUN based architectures"    OFF)
OPTION(CRAY   "Format code for CRAY based architectures"   OFF)
OPTION(CRAYX1 "Format code for CRAYX1 based architectures" OFF)
MARK_AS_ADVANCED(IBM SGI SUN CRAY CRAYX1)

OPTION(DEBUG_FULL_STACK "Write the detailed stack trace during debugging" OFF)
OPTION(DEBUG_FLUSH_MESSAGES "Do not allow caching of screen printed messages" OFF)
OPTION(DEBUG_LOG_LEVEL "Force debug log level for screen messages" OFF)
OPTION(DEBUG_ALL_TRACE "Write all tracing debug information to screen output during debugging" OFF)
OPTION(DEBUG_GLOBALIO_TRACE "Write tracing debug information from the GLOBALIO module" OFF)
OPTION(DEBUG_WRITER_TRACE "Write tracing debug information for writer processors" OFF)
OPTION(DEBUG_WRITE_OUTPUT_TRACE "Write tracing debug information for the write output module" OFF)
OPTION(DEBUG_WIND_TRACE "Write tracing debug information for the wind module" OFF)
OPTION(DEBUG_WEIR_TRACE "Write tracing debug information for the weir module" OFF)
OPTION(DEBUG_TVW_TRACE "Write tracing debug information for the time varying weir module" OFF)
OPTION(DEBUG_VSMY_TRACE "Write tracing debug information for the 3D momentum equation module" OFF)
OPTION(DEBUG_TIMESTEP_TRACE "Write the tracing debug information for the timestepping module" OFF)
OPTION(DEBUG_SUBPREP_TRACE "Write the tracing debug information for the subdomain prep module" OFF)
OPTION(DEBUG_SUBDOMAIN_TRACE "Write the tracing debug information for the subdomain module" OFF)
OPTION(DEBUG_READ_INPUT_TRACE "Write the tracing debug information for the read_input module" OFF)
OPTION(DEBUG_OWIWIND_TRACE "Write the tracing debug information for the OWI wind module" OFF)
OPTION(DEBUG_NODALATTR_TRACE "Write the tracing debug information for the nodal attributes module" OFF)
OPTION(DEBUG_NETCDF_TRACE "Write the tracing debug information for the netCDF modle" OFF)
OPTION(DEBUG_MESSENGER_TRACE "Write the tracing debug information for the message passing (MPI) module" OFF)
OPTION(DEBUG_MESH_TRACE "Write the tracing debug information for the mesh module" OFF)
OPTION(DEBUG_HOTSTART_TRACE "Write the tracing debug information for the hotstart module" OFF)
OPTION(DEBUG_GLOBAL_TRACE "Write the tracing debug information for the global module" OFF)
OPTION(DEBUG_HARM_TRACE "Write the tracing debug information for the harmonics module" OFF)
OPTION(DEBUG_COLDSTART_TRACE "Write the tracing debug information for the coldstart module" OFF)
OPTION(DEBUG_COUPLE2SWAN_TRACE "Write the tracing debug information for the couple2swan module" OFF)
OPTION(DEBUG_ADCIRC_TRACE "Write the tracing debug information for the main ADCIRC module" OFF)
MARK_AS_ADVANCED(DEBUG_GLOBALIO_TRACE DEBUG_WRITER_TRACE DEBUG_WRITE_OUTPUT_TRACE DEBUG_WIND_TRACE 
                 DEBUG_WEIR_TRACE DEBUG_TVW_TRACE DEBUG_VSMY_TRACE DEBUG_TIMESTEP_TRACE DEBUG_SUBPREP_TRACE 
                 DEBUG_SUBDOMAIN_TRACE DEBUG_READ_INPUT_TRACE DEBUG_OWIWIND_TRACE DEBUG_NODALATTR_TRACE 
                 DEBUG_NETCDF_TRACE DEBUG_MESSENGER_TRACE DEBUG_MESH_TRACE DEBUG_HOTSTART_TRACE DEBUG_GLOBAL_TRACE
                 DEBUG_HARM_TRACE DEBUG_COLDSTART_TRACE DEBUG_COUPLE2SWAN_TRACE DEBUG_ADCIRC_TRACE IBM )

OPTION(ENABLE_POWELL "Force Powell wind drag to be enabled. Warning: Overrides any other options specified at run time." OFF)
MARK_AS_ADVANCED(ENABLE_POWELL)

OPTION(VECTOR_COMPUTER "Assume the system is a vector computer" OFF)
MARK_AS_ADVANCED(VECTOR_COMPUTER)
###########################################################################



###########################################################################
#                       COMPILER FLAGS
###########################################################################
IF(PRECISION_8BYTE)
    SET(PRECISION_FLAG "-DREAL8")
ELSE(PRECISION_8BYTE)
    SET(PRECISION_FLAG "-DREAL4")
ENDIF(PRECISION_8BYTE)

IF(SUN)
    SET(MACHINE_FLAG "-DCMACHSUN")
ELSEIF(SGI)
    SET(MACHINE_FLAG "-DSGI")
ELSEIF(CRAY)
    SET(MACHINE_FLAG "-DCRAY")
ELSEIF(CRAYX1)
    SET(MACHINE_FLAG "-DCRAYX1")
ELSEIF(UNIX)
    SET(MACHINE_FLAG "-DLINUX")
ELSEIF(CYGWIN)
    SET(MACHINE_FLAG "-DLINUX")
ELSEIF(WIN32)
    SET(MACHINE_FLAG "-DWINDOWS")
ELSEIF(APPLE)
    SET(MACHINE_FLAG "-DLINUX")
ENDIF(SUN)

IF(SUN_MACHINE)
    SET(MACHINE_FLAG "${MACHINE_FLAG} -CMACHSUN")
ENDIF(SUN_MACHINE)

IF(ENABLE_WARN_ELEV_DEBUG)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DDEBUG_WARN_ELEV")
ENDIF(ENABLE_WARN_ELEV_DEBUG)

IF(ENABLE_POWELL)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DPOWELL")
ENDIF(ENABLE_POWELL)

IF(DEBUG_FULL_STACK)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DFULL_STACK")
ENDIF(DEBUG_FULL_STACK)

IF(DEBUG_ALL_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DALL_TRACE")
ENDIF(DEBUG_ALL_TRACE)

IF(DEBUG_FLUSH_MESSAGES)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DFLUSH_MESSAGES")
ENDIF(DEBUG_FLUSH_MESSAGES)

IF(DEBUG_LOG_LEVEL)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DDEBUG")
ENDIF(DEBUG_LOG_LEVEL)

IF(DEBUG_GLOBALIO_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DGLOBALIO_TRACE")
ENDIF(DEBUG_GLOBALIO_TRACE)

IF(DEBUG_WRITER_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DWRITER_TRACE")
ENDIF(DEBUG_WRITER_TRACE)

IF(DEBUG_WRITE_OUTPUT_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DWRITE_OUTPUT_TRACE")
ENDIF(DEBUG_WRITE_OUTPUT_TRACE)

IF(DEBUG_WIND_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DWIND_TRACE")
ENDIF(DEBUG_WIND_TRACE)

IF(DEBUG_WEIR_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DWEIR_TRACE")
ENDIF(DEBUG_WEIR_TRACE)

IF(DEBUG_TVW_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DTVW_TRACE")
ENDIF(DEBUG_TVW_TRACE)

IF(DEBUG_VSMY_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DVSMY_TRACE")
ENDIF(DEBUG_VSMY_TRACE)

IF(DEBUG_TIMESTEP_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DTIMESTEP_TRACE")
ENDIF(DEBUG_TIMESTEP_TRACE)

IF(DEBUG_SUBPREP_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DSUBPREP_TRACE")
ENDIF(DEBUG_SUBPREP_TRACE)

IF(DEBUG_SUBDOMAIN_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DSUBDOMAIN_TRACE")
ENDIF(DEBUG_SUBDOMAIN_TRACE)

IF(DEBUG_READ_INPUT_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DSUBDOMAIN_TRACE")
ENDIF(DEBUG_READ_INPUT_TRACE)

IF(DEBUG_OWIWIND_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DOWIWIND_TRACE")
ENDIF(DEBUG_OWIWIND_TRACE)

IF(DEBUG_NODALATTR_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DNODALATTR_TRACE")
ENDIF(DEBUG_NODALATTR_TRACE)
    
IF(DEBUG_NETCDF_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DNETCDF_TRACE")
ENDIF(DEBUG_NETCDF_TRACE)
    
IF(DEBUG_MESSENGER_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DMESSENGER_TRACE")
ENDIF(DEBUG_MESSENGER_TRACE)
    
IF(DEBUG_MESH_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DMESH_TRACE")
ENDIF(DEBUG_MESH_TRACE)

IF(DEBUG_HOTSTART_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DHOTSTART_TRACE")
ENDIF(DEBUG_HOTSTART_TRACE)

IF(DEBUG_GLOBAL_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DGLOBAL_TRACE")
ENDIF(DEBUG_GLOBAL_TRACE)

IF(DEBUG_HARM_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DHARM_TRACE")
ENDIF(DEBUG_HARM_TRACE)

IF(DEBUG_COLDSTART_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCOLDSTART_TRACE")
ENDIF(DEBUG_COLDSTART_TRACE)

IF(DEBUG_COUPLE2SWAN_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCOUPLE2SWAN_TRACE")
ENDIF(DEBUG_COUPLE2SWAN_TRACE)

IF(DEBUG_ADCIRC_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DADCIRC_TRACE")
ENDIF(DEBUG_ADCIRC_TRACE)

IF(IBM)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DIBM")
ENDIF(IBM) 

IF(VECTOR_COMPUTER)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCVEC")
ELSE(VECTOR_COMPUTER)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCSCA")
ENDIF(VECTOR_COMPUTER)

SET(SWAN_FLAG "-DCSWAN")
SET(PREP_SWAN_FLAG "-DADCSWAN")
SET(ADCIRC_MPI_FLAG "-DCMPI ${MPIMOD_FLAG}")


###########################################################################
#...Check that the netCDF library specified works with the compilers
IF(ENABLE_OUTPUT_NETCDF)
    
    SET(netcdf3_f90_code 
"
        PROGRAM netCDF3Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID

            IERR = NF90_OPEN('test.nc',NF90_NOWRITE,NCID)

        END PROGRAM
"
    )
    SET(netcdf4_f90_code
"
        PROGRAM netCDF4Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID
            INTEGER :: VARID

            IERR = NF90_DEF_VAR_DEFLATE(NCID,VARID,1,1,2)

        END PROGRAM
"
    )

    IF(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
        MESSAGE(SEND_ERROR "Specify the netCDF path on the following screen")
    ELSE(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
        
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf3check.f90" "${netcdf3_f90_code}")
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf4check.f90" "${netcdf4_f90_code}")
        TRY_COMPILE(NETCDF_TEST1 "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf3check.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDFHOME}/include" "-DLINK_DIRECTORIES=${NETCDFHOME}/lib" LINK_LIBRARIES netcdf LINK_LIBRARIES netcdff OUTPUT_VARIABLE LOG1)
        TRY_COMPILE(NETCDF_TEST2 "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf4check.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDFHOME}/include" "-DLINK_DIRECTORIES=${NETCDFHOME}/lib" LINK_LIBRARIES netcdf LINK_LIBRARIES netcdff OUTPUT_VARIABLE LOG2)

        IF(NETCDF_TEST1)
            SET(NETCDF_FLAG "-DADCNETCDF -I${NETCDFHOME}/include")
            SET(NETCDF_LINKER_FLAG "-L${NETCDFHOME}/lib")
            SET(NETCDF_WORKING TRUE)
            IF(NETCDF_TEST2)
                SET(NETCDF_COMPRESSION_FLAG "-DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE")
            ELSE(NETCDF_TEST2)
                SET(NETCDF_COMPRESSION_FLAG "")
            ENDIF(NETCDF_TEST2)
        ELSE(NETCDF_TEST1)
            MESSAGE(SEND_ERROR "The netCDF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable netCDF. Ensure that you specify the same compilers to build ADCIRC as were used to build the netCDF library.")
            SET(NETCDF_WORKING FALSE)
        ENDIF(NETCDF_TEST1)

    ENDIF(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
ENDIF(ENABLE_OUTPUT_NETCDF)
###########################################################################


###########################################################################
#...Check that the XDMF library works with the specified compilers
IF(ENABLE_OUTPUT_XDMF)

    SET(xdmf_f90_code
"
        PROGRAM XDMFCHECK
            IMPLICIT NONE
            INCLUDE 'Xdmf.f'
            INTEGER :: xdmfunit
            CALL xdmfInit(xdmfunit)
        END PROGRAM
"
)
    IF(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
        MESSAGE(SEND_ERROR "Specify the XDMF path on the following screen")
    ELSE(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" "${xdmf_f90_code}")
        TRY_COMPILE(XDMF_TEST "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${XDMFHOME}/include" "-DLINK_DIRECTORIES=${XDMFHOME}/lib" LINK_LIBRARIES XdmfCore LINK_LIBRARIES XdmfUtils LINK_LIBRARIES Xdmf OUTPUT_VARIABLE XDMFLOG)

        IF(XDMF_TEST)
            SET(XDMF_FLAG "-DADCXDMF -I${XDMFHOME}/include -I${CMAKE_SOURCE_DIR}/src")
            SET(XDMF_LINKER_FLAG "-L${XDMFHOME}/lib")
            SET(XDMF_WORKING TRUE)
        ELSE(XDMF_TEST)
            MESSAGE(SEND_ERROR "The XDMF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable XDMF. Ensure that you specify the same compilers to build ADCIRC as were used to build the XDMF library.")
            SET(XDMF_WORKING FALSE)
        ENDIF(XDMF_TEST)
    ENDIF(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
ENDIF(ENABLE_OUTPUT_XDMF)
###########################################################################



###########################################################################
#                         ADCIRC VERSION
###########################################################################
ADD_LIBRARY(version version.F)
#ADD_LIBRARY(version ${CMAKE_BINARY_DIR}/version.F)
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_BINARY_DIR}/version_cmake.F
#    COMMAND ./generateVersion.sh ${CMAKE_BINARY_DIR}/version_cmake.F
#    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
#    COMMENT "Generating ADCIRC version...")
SET_TARGET_PROPERTIES(version PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/version_mod)
###########################################################################


###########################################################################
#                            MKDIR
###########################################################################
ADD_LIBRARY(mkdir prep/mkdir.c)
TARGET_INCLUDE_DIRECTORIES(mkdir PRIVATE ${CMAKE_SOURCE_DIR}/src)
###########################################################################


###########################################################################
#                            ADCIRC
###########################################################################
IF(BUILD_ADCIRC)
    
    #...Note that we need to build ADCIRC in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    SET( ADCIRC1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                          src/global.F src/boundaries.F src/mesh.F 
                          src/global_3dvs.F src/harm.F wind/vortex.F 
                          src/wind.F src/owiwind.F src/rs2.F 
                          src/owi_ice.F src/itpackv.F src/nodalattr.F 
                          src/globalio.F src/netcdfio.F src/control.F
                          src/xdmfio.F )

    SET( ADCIRC2_SOURCES  src/write_output.F 
                          src/couple2swan.F src/adcirc.F src/subdomain.F 
                          src/weir_boundary.F src/read_input.F src/cstart.F 
                          src/hstart.F src/timestep.F src/vsmy.F 
                          src/transport.F src/driver.F )

    ADD_LIBRARY(templib_adcirc1 ${ADCIRC1_SOURCES})
    
    ADD_EXECUTABLE(adcirc ${ADCIRC2_SOURCES})
    
    SET(ADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1 netcdf netcdff)
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${ADCIRC_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1 netcdf netcdff XdmfCore XdmfUtils Xdmf)
        TARGET_INCLUDE_DIRECTORIES(adcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    TARGET_INCLUDE_DIRECTORIES(adcirc          PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcirc1 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)

    SET_TARGET_PROPERTIES(templib_adcirc1 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcirc_mod)
    SET_TARGET_PROPERTIES(templib_adcirc1 PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    SET_TARGET_PROPERTIES(adcirc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcirc_mod)
    SET_TARGET_PROPERTIES(adcirc PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    
    TARGET_LINK_LIBRARIES(adcirc  version)
    TARGET_LINK_LIBRARIES(templib_adcirc1 version)
    
    ADD_DEPENDENCIES(adcirc          templib_adcirc1)
    ADD_DEPENDENCIES(templib_adcirc1 version)
    
    INSTALL(TARGETS adcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCIRC)
###########################################################################



###########################################################################
#                             ADCPREP
###########################################################################
IF(BUILD_ADCPREP)
    
    SET( METIS_SOURCES    metis/Lib/coarsen.c metis/Lib/fm.c metis/Lib/initpart.c metis/Lib/match.c 
                          metis/Lib/ccgraph.c metis/Lib/memory.c metis/Lib/pmetis.c metis/Lib/pqueue.c 
                          metis/Lib/refine.c metis/Lib/util.c metis/Lib/timing.c metis/Lib/debug.c 
                          metis/Lib/bucketsort.c metis/Lib/graph.c metis/Lib/stat.c metis/Lib/kmetis.c 
                          metis/Lib/kwayrefine.c metis/Lib/kwayfm.c metis/Lib/balance.c metis/Lib/ometis.c 
                          metis/Lib/srefine.c metis/Lib/sfm.c metis/Lib/separator.c metis/Lib/mincover.c 
                          metis/Lib/mmd.c metis/Lib/mesh.c metis/Lib/meshpart.c metis/Lib/frename.c 
                          metis/Lib/fortran.c metis/Lib/myqsort.c metis/Lib/compress.c metis/Lib/parmetis.c 
                          metis/Lib/estmem.c metis/Lib/mpmetis.c metis/Lib/mcoarsen.c metis/Lib/mmatch.c 
                          metis/Lib/minitpart.c metis/Lib/mbalance.c metis/Lib/mrefine.c metis/Lib/mutil.c 
                          metis/Lib/mfm.c metis/Lib/mkmetis.c metis/Lib/mkwayrefine.c metis/Lib/mkwayfmh.c 
                          metis/Lib/mrefine2.c metis/Lib/minitpart2.c metis/Lib/mbalance2.c metis/Lib/mfm2.c 
                          metis/Lib/kvmetis.c metis/Lib/kwayvolrefine.c metis/Lib/kwayvolfm.c 
                          metis/Lib/subdomains.c )
    
    SET( ADCPREP_SOURCES  src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F 
                          src/mesh.F src/global_3dvs.F wind/vortex.F src/owiwind.F src/rs2.F 
                          src/owi_ice.F src/wind.F prep/presizes.F prep/pre_global.F prep/metis.F 
                          prep/subprep.F prep/adcprep.F prep/decomp.F prep/prep_weir.F 
                          src/nodalattr.F src/harm.F prep/read_global.F src/netcdfio.F prep/prep.F  
                          prep/interp.F prep/machdep.F )
    
    ADD_LIBRARY(metis STATIC ${METIS_SOURCES})
    
    ADD_EXECUTABLE(adcprep ${ADCPREP_SOURCES})
    
    SET_TARGET_PROPERTIES(adcprep PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcprep_mod)
    
    SET(ADCPREP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${MACHINE_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCPREP}")
    
    IF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${PREP_SWAN_FLAG}")
    ENDIF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
   
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcprep netcdf netcdff)
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCPREP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${ADCPREP_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcprep netcdf netcdff XdmfCore XdmfUtils Xdmf)
        TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
   
    SET_TARGET_PROPERTIES(adcprep PROPERTIES COMPILE_FLAGS ${ADCPREP_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_LINK_LIBRARIES(adcprep metis version mkdir)
    
    ADD_DEPENDENCIES(adcprep metis)
    ADD_DEPENDENCIES(adcprep mkdir)
    ADD_DEPENDENCIES(adcprep version)
    
    INSTALL(TARGETS adcprep RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCPREP)
###########################################################################


###########################################################################
#                             ADCSWAN
###########################################################################
IF(BUILD_ADCSWAN AND PERL_FOUND)

    #...Note that we need to build ADCSWAN in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    ADD_EXECUTABLE(adcswan src/driver.F)
   
    SET( SWAN1SERIAL_SOURCES ${CMAKE_BINARY_DIR}/swan_serial_source/swmod1.f ${CMAKE_BINARY_DIR}/swan_serial_source/swmod2.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/swan_serial_source/m_constants.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/serv_xnl4v5.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGriddata.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanCompdata.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/couple2adcirc.f90 )
    
    SET( SWAN2SERIAL_SOURCES ${CMAKE_BINARY_DIR}/swan_serial_source/swanmain.f ${CMAKE_BINARY_DIR}/swan_serial_source/swanpre1.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/swanpre2.f ${CMAKE_BINARY_DIR}/swan_serial_source/swancom1.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/swancom2.f ${CMAKE_BINARY_DIR}/swan_serial_source/swancom3.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/swancom4.f ${CMAKE_BINARY_DIR}/swan_serial_source/swancom5.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/swanout1.f ${CMAKE_BINARY_DIR}/swan_serial_source/swanout2.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/swanser.f ${CMAKE_BINARY_DIR}/swan_serial_source/swanparll.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanReadADCGrid.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanReadEasymeshGrid.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanCheckGrid.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGridTopology.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGridCell.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPrintGridInfo.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPointinMesh.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPrepComp.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanCompUnstruc.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPropvelX.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPropvelS.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanTranspX.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGSECorr.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanInterpolatePoint.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanInterpolateOutput.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanConvStopc.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanCrossObstacle.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanIntgratSpc.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanReadfort18.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanSumOverNodes.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/SwanMaxOverNodes.f90 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/ocpids.f ${CMAKE_BINARY_DIR}/swan_serial_source/ocpcre.f 
                             ${CMAKE_BINARY_DIR}/swan_serial_source/ocpmix.f )
    
    SET( ADCSWAN1_SOURCES    src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F src/mesh.F 
                             src/global_3dvs.F src/harm.F wind/vortex.F src/wind.F src/owiwind.F src/rs2.F
                             src/owi_ice.F src/itpackv.F src/nodalattr.F src/globalio.F src/netcdfio.F 
                             src/control.F src/xdmfio.F )

    SET( ADCSWAN2_SOURCES    src/write_output.F src/couple2swan.F )

    SET( ADCSWAN3_SOURCES    src/adcirc.F src/subdomain.F src/weir_boundary.F src/read_input.F src/cstart.F 
                             src/hstart.F src/timestep.F src/vsmy.F src/transport.F )
    
    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/swan_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/swan_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/swan_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/swan_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/swan_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/swan_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ENDIF(WIN32)

    ADD_LIBRARY(templib_swan1serial STATIC ${SWAN1SERIAL_SOURCES})
    ADD_LIBRARY(templib_swan2serial STATIC ${SWAN2SERIAL_SOURCES})
    ADD_LIBRARY(templib_adcswan1 STATIC ${ADCSWAN1_SOURCES})
    ADD_LIBRARY(templib_adcswan2 STATIC ${ADCSWAN2_SOURCES})
    ADD_LIBRARY(templib_adcswan3 STATIC ${ADCSWAN3_SOURCES})
   
    TARGET_LINK_LIBRARIES(adcswan templib_swan1serial)
    TARGET_LINK_LIBRARIES(adcswan templib_swan2serial)
    
    SET(ADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_OPTION_FLAGS}")
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${ADCSWAN_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(adcswan PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    SET_TARGET_PROPERTIES(templib_swan1serial PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_swan2serial PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan1    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan2    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan3    PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan             PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcswan_mod)
    SET_TARGET_PROPERTIES(templib_adcswan1    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_adcswan2    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_adcswan3    PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(adcswan             PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(templib_swan1serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(templib_swan2serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
    
    TARGET_INCLUDE_DIRECTORIES(adcswan  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan1 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan2 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcswan3 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2serial PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)

    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version
                                      netcdf netcdff )
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version
                                      netcdf netcdff XdmfCore XdmfUtils Xdmf )
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcswan templib_adcswan3 templib_adcswan2
                                      templib_swan2serial templib_adcswan1
                                      templib_swan1serial version )
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    ADD_DEPENDENCIES(adcswan             templib_adcswan3)
    ADD_DEPENDENCIES(templib_adcswan3    templib_swan2serial)
    ADD_DEPENDENCIES(templib_swan2serial templib_adcswan2)
    ADD_DEPENDENCIES(templib_adcswan3    templib_adcswan2)
    ADD_DEPENDENCIES(templib_adcswan2    templib_adcswan1)
    ADD_DEPENDENCIES(templib_adcswan1    templib_swan1serial)
    ADD_DEPENDENCIES(templib_adcswan1    version)
    
    INSTALL(TARGETS adcswan RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCSWAN AND PERL_FOUND)
###########################################################################



###########################################################################
#                             SWAN 
###########################################################################
IF(BUILD_SWAN AND PERL_FOUND)
    
    SET( SWANONLY_SERIAL_SOURCES ${CMAKE_BINARY_DIR}/swanonly_serial_source/swmod1.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/couple2adcirc.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanmain.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanpre2.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swancom2.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swancom4.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanout1.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanser.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/ocpids.f ${CMAKE_BINARY_DIR}/swanonly_serial_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/swanonly_serial_source/ocpmix.f )
    
    ADD_EXECUTABLE(swan ${SWANONLY_SERIAL_SOURCES})

    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/swanonly_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/swanonly_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/swanonly_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/swanonly_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES} 
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/swanonly_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/swanonly_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )
    ENDIF(WIN32)
    
    SET_TARGET_PROPERTIES(swan PROPERTIES COMPILE_FLAGS "${Fortran_COMPILER_SPECIFIC_FLAG} 
                                                         ${ADDITIONAL_FLAGS_SWAN}")
    SET_TARGET_PROPERTIES(swan PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/swan_mod)
    
    TARGET_INCLUDE_DIRECTORIES(swan PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/swan_mod)
    
    INSTALL(TARGETS swan RUNTIME DESTINATION bin)

ENDIF(BUILD_SWAN AND PERL_FOUND)
###########################################################################



###########################################################################
#                           PADCIRC
###########################################################################
IF(BUILD_PADCIRC)
    
    #...Note that we need to build PADCIRC in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    SET( PADCIRC1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                           src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )

    SET( PADCIRC2_SOURCES  src/mesh.F src/harm.F wind/vortex.F src/wind.F 
                           src/owiwind.F src/rs2.F src/owi_ice.F 
                           src/itpackv.F src/nodalattr.F src/globalio.F 
                           src/netcdfio.F src/control.F src/xdmfio.F )

    SET( PADCIRC3_SOURCES  src/writer.F )

    SET( PADCIRC4_SOURCES  src/write_output.F src/couple2swan.F src/adcirc.F src/subdomain.F 
                           src/weir_boundary.F src/read_input.F src/cstart.F 
                           src/hstart.F src/timestep.F src/vsmy.F 
                           src/transport.F src/driver.F )

    ADD_LIBRARY(templib_padcirc1 ${PADCIRC1_SOURCES})
    ADD_LIBRARY(templib_padcirc2 ${PADCIRC2_SOURCES})
    ADD_LIBRARY(templib_padcirc3 ${PADCIRC3_SOURCES})
    
    ADD_EXECUTABLE(padcirc ${PADCIRC4_SOURCES})

    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1 PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(padcirc          PRIVATE prep)
    
    SET(PADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${PADCIRC_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(padcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
    
    SET_TARGET_PROPERTIES(templib_padcirc1 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc1 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(templib_padcirc2 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc2 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    SET_TARGET_PROPERTIES(templib_padcirc3 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc3 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(padcirc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(padcirc PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc2  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcirc           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc2  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(padcirc           PRIVATE ${MPI_Fortran_INCLUDE_PATH})

    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir netcdf netcdff ${MPI_Fortran_LIBRARIES})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir netcdf netcdff XdmfCore XdmfUtils Xdmf
                                      ${MPI_Fortran_LIBRARIES})
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir ${MPI_Fortran_LIBRARIES})
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
   
    ADD_DEPENDENCIES(padcirc          templib_padcirc3)
    ADD_DEPENDENCIES(templib_padcirc2 templib_padcirc1)
    ADD_DEPENDENCIES(templib_padcirc3 templib_padcirc2)
    ADD_DEPENDENCIES(templib_padcirc1 version)
    ADD_DEPENDENCIES(templib_padcirc1 mkdir)
    
    INSTALL(TARGETS padcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_PADCIRC)
###########################################################################


###########################################################################
#                            PADCSWAN
###########################################################################
IF(BUILD_PADCSWAN AND PERL_FOUND)

    #...Note that we need to build PADCSWAN in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    ADD_EXECUTABLE(padcswan src/driver.F)
   
    SET( SWAN1PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/swan_parallel_source/swmod1.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/swan_parallel_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/couple2adcirc.f90 )
    
    SET( SWAN2PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/swan_parallel_source/swanmain.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/swanpre2.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/swancom2.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/swancom4.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/swanout1.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/swanser.f ${CMAKE_BINARY_DIR}/swan_parallel_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/ocpids.f ${CMAKE_BINARY_DIR}/swan_parallel_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/swan_parallel_source/ocpmix.f )
     
    SET( PADCSWAN1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                            src/global.F src/boundaries.F src/global_3dvs.F
                            src/messenger.F )

    SET( PADCSWAN2_SOURCES  src/mesh.F src/harm.F wind/vortex.F src/wind.F 
                            src/owiwind.F src/rs2.F src/owi_ice.F 
                            src/itpackv.F src/nodalattr.F src/globalio.F 
                            src/netcdfio.F src/control.F src/xdmfio.F )

    SET( PADCSWAN3_SOURCES  src/writer.F )
    
    SET( PADCSWAN4_SOURCES  src/write_output.F src/couple2swan.F )
    
    SET( PADCSWAN5_SOURCES  src/adcirc.F src/subdomain.F 
                            src/weir_boundary.F src/read_input.F src/cstart.F 
                            src/hstart.F src/timestep.F src/vsmy.F 
                            src/transport.F src/driver.F )
    
    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/swan_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/swan_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/swan_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/swan_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/swan_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/swan_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel SWAN Sources..."
        )
    ENDIF(WIN32)

    ADD_LIBRARY(templib_swan1parallel STATIC ${SWAN1PARALLEL_SOURCES})
    ADD_LIBRARY(templib_swan2parallel STATIC ${SWAN2PARALLEL_SOURCES})
    ADD_LIBRARY(templib_padcswan1     STATIC ${PADCSWAN1_SOURCES})
    ADD_LIBRARY(templib_padcswan2     STATIC ${PADCSWAN2_SOURCES})
    ADD_LIBRARY(templib_padcswan3     STATIC ${PADCSWAN3_SOURCES})
    ADD_LIBRARY(templib_padcswan4     STATIC ${PADCSWAN4_SOURCES})
    ADD_LIBRARY(templib_padcswan5     STATIC ${PADCSWAN5_SOURCES})
   
    SET(PADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel
                                       mkdir netcdf netcdff ${MPI_Fortran_LIBRARIES}) 

    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
   
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${PADCSWAN_LINKER_FLAGS} )
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel
                                       mkdir netcdf netcdff XdmfCore XdmfUtils Xdmf
                                       ${MPI_Fortran_LIBRARIES})
        
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        TARGET_LINK_LIBRARIES(padcswan templib_padcswan5 templib_padcswan4
                                       templib_padcswan3 templib_swan2parallel
                                       templib_padcswan2 templib_padcswan1 
                                       version templib_swan1parallel mkdir
                                       ${MPI_Fortran_LIBRARIES})

    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
    SET_TARGET_PROPERTIES(templib_swan1parallel PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_swan2parallel PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan1     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan2     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan3     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan4     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan5     PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan              PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcswan_mod)
    SET_TARGET_PROPERTIES(templib_padcswan1     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan2     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan3     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan4     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(templib_padcswan5     PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan              PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    
    TARGET_INCLUDE_DIRECTORIES(padcswan              PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan1     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan2     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan3     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan4     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan5     PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_INCLUDE_DIRECTORIES(padcswan              PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan1     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan2     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan3     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan4     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcswan5     PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_swan2parallel PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(swan1parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(swan2parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF()
    
    ADD_DEPENDENCIES(padcswan              templib_swan2parallel)
    ADD_DEPENDENCIES(padcswan              templib_padcswan3)
    ADD_DEPENDENCIES(templib_swan2parallel templib_padcswan5)
    ADD_DEPENDENCIES(templib_padcswan5     templib_padcswan4)
    ADD_DEPENDENCIES(templib_padcswan4     templib_padcswan2)
    ADD_DEPENDENCIES(templib_padcswan4     templib_padcswan3)
    ADD_DEPENDENCIES(templib_padcswan3     templib_padcswan2)
    ADD_DEPENDENCIES(templib_padcswan2     templib_padcswan1)
    ADD_DEPENDENCIES(templib_padcswan1     templib_swan1parallel)
    ADD_DEPENDENCIES(templib_padcswan1     mkdir)
    
    INSTALL(TARGETS padcswan RUNTIME DESTINATION bin)

ENDIF(BUILD_PADCSWAN AND PERL_FOUND)
###########################################################################



###########################################################################
#                          PUNSWAN 
###########################################################################
IF(BUILD_PUNSWAN AND PERL_FOUND)
    
    SET( SWANONLY1_PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swmod1.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swmod2.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanSpectPart.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/m_constants.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/m_fileio.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/serv_xnl4v5.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/mod_xnl4v5.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGriddata.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGridobjects.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanCompdata.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/couple2adcirc.f90 ) 
                               
    SET( SWANONLY2_PARALLEL_SOURCES ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanmain.f  ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanpre1.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanpre2.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swancom1.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swancom2.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swancom3.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swancom4.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swancom5.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanout1.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanout2.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanser.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/swanparll.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanReadGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanReadADCGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanReadTriangleGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanReadEasymeshGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanInitCompGrid.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanCheckGrid.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanCreateEdges.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGridTopology.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGridVert.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGridCell.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGridFace.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPrintGridInfo.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanFindPoint.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPointinMesh.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanBpntlist.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPrepComp.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanVertlist.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanCompUnstruc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanDispParm.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPropvelX.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanSweepSel.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPropvelS.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanTranspAc.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanTranspX.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanDiffPar.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGSECorr.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanGradDepthorK.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanInterpolatePoint.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanInterpolateAc.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanInterpolateOutput.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanConvAccur.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanConvStopc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanFindObstacles.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanCrossObstacle.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanComputeForce.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanIntgratSpc.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanBndStruc.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanReadfort18.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanPunCollect.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanSumOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanMinOverNodes.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/SwanMaxOverNodes.f90 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/ocpids.f ${CMAKE_BINARY_DIR}/swanonly_parallel_source/ocpcre.f 
                               ${CMAKE_BINARY_DIR}/swanonly_parallel_source/ocpmix.f )

    SET( MSGLIB_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                         src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )
    
    ADD_LIBRARY(templib_punmsglib ${MSGLIB_SOURCES})
    ADD_LIBRARY(templib_punswan1  ${SWANONLY1_PARALLEL_SOURCES})
    ADD_EXECUTABLE(punswan        ${SWANONLY2_PARALLEL_SOURCES})

    #...SWAN Configuration
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/swanonly_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/swanonly_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/swanonly_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/swanonly_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/swanonly_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/swanonly_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel unSWAN Sources..."
        )
    ENDIF(WIN32)
    
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(punswan1       PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(punswan        PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF()
    
    SET_TARGET_PROPERTIES(templib_punswan1  PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    SET_TARGET_PROPERTIES(templib_punmsglib PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    SET_TARGET_PROPERTIES(punswan           PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/punswan_mod)
    
    SET(MSGLIB_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    SET_TARGET_PROPERTIES(templib_punmsglib PROPERTIES COMPILE_FLAGS ${MSGLIB_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/punswan_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan           PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_punswan1  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_punswan1  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_punmsglib PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_punmsglib PRIVATE ${MPI_Fortran_INCLUDE_PATH})

    TARGET_LINK_LIBRARIES(punswan templib_punmsglib templib_punswan1 mkdir ${MPI_Fortran_LIBRARIES})

    ADD_DEPENDENCIES(templib_punmsglib templib_punswan1)
    ADD_DEPENDENCIES(punswan           templib_punswan1)
    ADD_DEPENDENCIES(punswan           templib_punmsglib)
    ADD_DEPENDENCIES(templib_punmsglib version)
    ADD_DEPENDENCIES(templib_punmsglib mkdir)
    
    INSTALL(TARGETS punswan RUNTIME DESTINATION bin)

ENDIF(BUILD_PUNSWAN AND PERL_FOUND)
###########################################################################


###########################################################################
#                              ASWIP
###########################################################################
IF(BUILD_ASWIP)
    
    SET(ASWIP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ASWIP}")

    SET(ASWIP_SOURCES src/sizes.F src/global.F src/global_3dvs.F src/boundaries.F 
                      src/mesh.F src/wind.F src/owiwind.F KDTREE2/kdtree2.F src/owi_ice.F 
                      wind/vortex.F wind/aswip_1.0.3.F )


    ADD_EXECUTABLE(aswip ${ASWIP_SOURCES})


    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff )

    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        SET(ASWIP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${ASWIP_LINKER_FLAGS} )
        
        TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_SOURCE_DIR}/src)
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff XdmfCore XdmfUtils Xdmf )

    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING) 

    TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_LINK_LIBRARIES(aswip version)

    SET_TARGET_PROPERTIES(aswip PROPERTIES COMPILE_FLAGS ${ASWIP_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(aswip PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/aswip_mod)

    INSTALL(TARGETS aswip RUNTIME DESTINATION bin)

ENDIF(BUILD_ASWIP)
###########################################################################


###########################################################################
# UTILITY PROGRAMS
###########################################################################
IF(BUILD_UTILITIES)
    SET(UTILITY_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS}")
    
    ADD_EXECUTABLE(adccmp util/adccmp.F)
    ADD_EXECUTABLE(p15 wind/p15.F)
    ADD_EXECUTABLE(owi22 wind/owi22.F)
    ADD_EXECUTABLE(build13 util/build13.F)
    ADD_EXECUTABLE(buildstwave23 util/buildstwave23.F)
    ADD_EXECUTABLE(hot2asc util/hot2asc.F)
    ADD_EXECUTABLE(inflate util/inflate.F)
    ADD_EXECUTABLE(hstime util/hstime.F)
    
    SET_TARGET_PROPERTIES(adccmp PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(p15 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(owi22 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(build13 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(buildstwave23 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(hot2asc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(inflate PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(hstime PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)

    SET_TARGET_PROPERTIES(adccmp PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(p15 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(owi22 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(build13 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(buildstwave23 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(hot2asc PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(inflate PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(hstime PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})

    INSTALL(TARGETS adccmp p15 owi22 build13 buildstwave23 hot2asc inflate hstime 
            RUNTIME DESTINATION bin)

ENDIF(BUILD_UTILITIES)
###########################################################################
