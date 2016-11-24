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

#...MPI Compilers
FIND_PACKAGE(MPI)

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

#...Language Specifications
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)

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
  MESSAGE("No known predefined Fortran extended line length flag known. Please set manually in Fortran_F77F90_FLAG")
  SET(Fortran_LINELENGTH_FLAG "" CACHE STRING "Compiler specific flag to enable extended Fortran line length")
  SET(Fortran_COMPILER_SPECIFIC_FLAG "" CACHE STRING "Compiler specific flags")
ENDIF(Fortran_COMPILER_NAME MATCHES "gfortran.*")

#...Substitute in the MPI compiler 
IF(MPI_FOUND)
    FIND_PROGRAM(MPI_Fortran_COMPILER
        NAMES mpif90
        HINTS "${MPI_BASE_DIR}"
        PATH_SUFFIXES bin
        DOC "MPI Fortran compiler. Warning: Overrides standard Fortran compiler")
    MARK_AS_ADVANCED(MPI_Fortran_COMPILER)
    
    IF(BUILD_PADCIRC OR BUILD_PADCSWAN OR BUILD_PUNSWAN)
        SET( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
    ENDIF(BUILD_PADCIRC OR BUILD_PADCSWAN OR BUILD_PUNSWAN)
ENDIF(MPI_FOUND)
###########################################################################



###########################################################################
# OPTIONS
###########################################################################
OPTION(PRECISION_8BYTE "Use 8-byte real numbers. This should be enabled." ON)
MARK_AS_ADVANCED(PRECISION_8BYTE)

#...Output options
OPTION(ENABLE_OUTPUT_NETCDF "Turn on netCDF output options" OFF)
IF(ENABLE_OUTPUT_NETCDF)
    OPTION(ENABLE_OUTPUT_XDMF "Turn on XDMF output options" OFF)
ENDIF(ENABLE_OUTPUT_NETCDF)    

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

#...Compilers
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER CMAKE_C_COMPILER CMAKE_Fortran_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS_RELEASE CMAKE_Fortran_FLAGS_RELEASE)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_DEBUG CMAKE_C_FLAGS_DEBUG CMAKE_Fortran_FLAGS_DEBUG)

#...Library paths
IF(ENABLE_OUTPUT_NETCDF)
    SET(NETCDFHOME "NETCDF-NOTFOUND" CACHE STRING "netCDF home path containing lib and include")
ENDIF(ENABLE_OUTPUT_NETCDF)

IF(ENABLE_OUTPUT_XDMF)
    SET(XDMFHOME "XDMF-NOTFOUND" CACHE STRING "XDMF home path containing lib and include")
ENDIF(ENABLE_OUTPUT_XDMF)


#...Additional flags
SET(ADDITIONAL_FLAGS_ADCIRC "" CACHE STRING "Additional flags to compile ADCIRC with")
IF(PERL_FOUND)
    SET(ADDITIONAL_FLAGS_SWAN "" CACHE STRING "Additional flags to compile SWAN with")
ENDIF(PERL_FOUND)
SET(ADDITIONAL_FLAGS_ADCPREP "" CACHE STRING "Additional flags to compile ADCPREP with")
SET(ADDITIONAL_FLAGS_ASWIP "" CACHE STRING "Additional flags to compile ASWIP with")
SET(ADDITIONAL_FLAGS_UTLIITIES "" CACHE STRING "Additional flags for utility programs")


#...Options enabled via compiler flags within the code
OPTION(ENABLE_WARN_ELEV_DEBUG "Enable writing of the fort.69 debug file" OFF)

IF(BUILD_ADCSWAN OR BUILD_PADCSWAN)
    OPTION(ENABLE_SWAN_MODIFIED_FRICTION "Enable Ethan/Joannes' modified friction" OFF)
ENDIF(BUILD_ADCSWAN OR BUILD_PADCSWAN)

OPTION(DEBUG_FULL_STACK "Write the detailed stack trace during debugging" OFF)
OPTION(DEBUG_TRACE "Write tracing information to screen output during debugging" OFF)
OPTION(ENABLE_POWELL "Force Powell wind drag to be enabled. Warning: Overrides any other options specified at run time." OFF)
MARK_AS_ADVANCED(ENABLE_POWELL)

option(VECTOR_COMPUTER "Assume the system is a vector computer" OFF)
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

IF(UNIX)
    SET(MACHINE_FLAG "-DLINUX")
ELSEIF(CYGWIN)
    SET(MACHINE_FLAG "-DLINUX")
ELSEIF(WIN32)
    SET(MACHINE_FLAG "-DWINDOWS")
ELSEIF(APPLE)
    SET(MACHINE_FLAG "-DLINUX")
ENDIF(UNIX)

IF(ENABLE_WARN_ELEV_DEBUG)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DDEBUG_WARN_ELEV")
ENDIF(ENABLE_WARN_ELEV_DEBUG)

IF(ENABLE_POWELL)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DPOWELL")
ENDIF(ENABLE_POWELL)

IF(DEBUG_FULL_STACK)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DFULL_STACK")
ENDIF(DEBUG_FULL_STACK)

IF(DEBUG_TRACE)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DALL_TRACE -DFLUSH_MESSAGES")
ENDIF(DEBUG_TRACE)

IF(VECTOR_COMPUTER)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCVEC")
ELSE(VECTOR_COMPUTER)
    SET(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} -DCSCA")
ENDIF(VECTOR_COMPUTER)

SET(SWAN_FLAG "-DCSWAN")
SET(PREP_SWAN_FLAG "-DADCSWAN")
SET(ADCIRC_MPI_FLAG "-DCMPI")

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
        MESSAGE("WARNING: Specify the netCDF path on the following screen")
    ELSE(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
        
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf3check.f90" "${netcdf3_f90_code}")
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf4check.f90" "${netcdf4_f90_code}")
        TRY_COMPILE(NETCDF_TEST1 "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf3check.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDFHOME}/include" "-DLINK_DIRECTORIES=${NETCDFHOME}/lib" LINK_LIBRARIES netcdf LINK_LIBRARIES netcdff OUTPUT_VARIABLE LOG1)
        TRY_COMPILE(NETCDF_TEST2 "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf4check.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDFHOME}/include" "-DLINK_DIRECTORIES=${NETCDFHOME}/lib" LINK_LIBRARIES netcdf LINK_LIBRARIES netcdff OUTPUT_VARIABLE LOG2)

        IF(NETCDF_TEST1)
            SET(NETCDF_FLAG "-DADCNETCDF -I${NETCDFHOME}/include")
            SET(NETCDF_LINKER_FLAG "-L${NETCDFHOME}/lib")
            IF(NETCDF_TEST2)
                SET(NETCDF_COMPRESSION_FLAG "-DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE")
            ELSE(NETCDF_TEST2)
                SET(NETCDF_COMPRESSION_FLAG "")
            ENDIF(NETCDF_TEST2)
        ELSE(NETCDF_TEST1)
            MESSAGE("ERROR: The netCDF library specified is not compatible with the specified compilers. It will not be enabled.")
        ENDIF(NETCDF_TEST1)

    ENDIF(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
ENDIF(ENABLE_OUTPUT_NETCDF)

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
        MESSAGE("WARNING: Specify the XDMF path on the following screen")
    ELSE(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" "${xdmf_f90_code}")
        TRY_COMPILE(XDMF_TEST "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDFHOME}/include" "-DINCLUDE_DIRECTORIES=${XDMFHOME}/include" "-DLINK_DIRECTORIES=${NETCDFHOME}/lib" "-DLINK_DIRECTORIES=${XDMFHOME}/lib" LINK_LIBRARIES netcdf LINK_LIBRARIES netcdff LINK_LIBRARIES XdmfCore LINK_LIBRARIES XdmfUtils LINK_LIBRARIES Xdmf OUTPUT_VARIABLE XDMFLOG)
    
        IF(XDMF_TEST)
            SET(XDMF_FLAG "-DADCXDMF -I${XDMFHOME}/include -I${CMAKE_SOURCE_DIR}/src")
            SET(XDMF_LINKER_FLAG "-L${XDMFHOME}/lib")
        ELSE(XDMF_TEST)
            MESSAGE("ERROR: The XDMF library specified is not compatible with the specified compilers. It will not be enabled.")
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
SET_TARGET_PROPERTIES(version PROPERTIES Fortran_MODULE_DIRECTORY version_mod)
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

    ADD_LIBRARY(adcirc1 ${ADCIRC1_SOURCES})
    
    ADD_EXECUTABLE(adcirc ${ADCIRC2_SOURCES})
    
    SET(ADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcirc ${CMAKE_BINARY_DIR}/libadcirc1.a -lnetcdf -lnetcdff)
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${ADCIRC_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcirc ${CMAKE_BINARY_DIR}/libadcirc1.a -lnetcdf -lnetcdff -lXdmfCore -lXdmfUtils -lXdmf)
        TARGET_INCLUDE_DIRECTORIES(adcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(adcirc ${CMAKE_BINARY_DIR}/libadcirc1.a)
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)

    TARGET_INCLUDE_DIRECTORIES(adcirc  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(adcirc1 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)

    SET_TARGET_PROPERTIES(adcirc1 PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(adcirc1 PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    SET_TARGET_PROPERTIES(adcirc PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(adcirc PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    
    TARGET_LINK_LIBRARIES(adcirc  version)
    TARGET_LINK_LIBRARIES(adcirc1 version)
    
    ADD_DEPENDENCIES(adcirc adcirc1)
    ADD_DEPENDENCIES(adcirc1 version)
    
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
    
    SET_TARGET_PROPERTIES(adcprep PROPERTIES Fortran_MODULE_DIRECTORY adcprep_mod)
    
    SET(ADCPREP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${MACHINE_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCPREP}")
    
    IF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${PREP_SWAN_FLAG}")
    ENDIF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
   
    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcprep ${CMAKE_BINARY_DIR}/libmetis.a -lnetcdf -lnetcdff)
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        SET(ADCPREP_COMPILER_FLAGS "${ADCPREP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCPREP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcprep PROPERTIES LINK_FLAGS ${ADCPREP_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcprep ${CMAKE_BINARY_DIR}/libmetis.a -lnetcdf -lnetcdff -lXdmfCore -lXdmfUtils -lXdmf)
        TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(adcprep ${CMAKE_BINARY_DIR}/libmetis.a)
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
   
    SET_TARGET_PROPERTIES(adcprep PROPERTIES COMPILE_FLAGS ${ADCPREP_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    
    TARGET_LINK_LIBRARIES(adcprep version mkdir)
    
    ADD_DEPENDENCIES(adcprep metis)
    
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

    ADD_LIBRARY(swan1serial STATIC ${SWAN1SERIAL_SOURCES})
    ADD_LIBRARY(swan2serial STATIC ${SWAN2SERIAL_SOURCES})
    ADD_LIBRARY(adcswan1 STATIC ${ADCSWAN1_SOURCES})
    ADD_LIBRARY(adcswan2 STATIC ${ADCSWAN2_SOURCES})
    ADD_LIBRARY(adcswan3 STATIC ${ADCSWAN3_SOURCES})
   
    #...Note - libraries appear twice since the objects in libraries need to
    #          see the other libraries. Do not remove the doubles
    TARGET_LINK_LIBRARIES(adcswan ${CMAKE_BINARY_DIR}/libswan1serial.a)
    TARGET_LINK_LIBRARIES(adcswan ${CMAKE_BINARY_DIR}/libswan2serial.a)
    
    SET(ADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_OPTION_FLAGS}")
    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        SET(ADCSWAN_COMPILER_FLAGS "${ADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcswan PROPERTIES LINK_FLAGS ${ADCSWAN_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(adcswan PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)

    SET_TARGET_PROPERTIES(swan1serial PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(swan2serial PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan1 PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan2 PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan3 PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan PROPERTIES Fortran_MODULE_DIRECTORY adcswan_mod)
    SET_TARGET_PROPERTIES(adcswan1 PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(adcswan2 PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(adcswan3 PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(adcswan PROPERTIES COMPILE_FLAGS ${ADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(swan1serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(swan2serial PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
    
    TARGET_INCLUDE_DIRECTORIES(adcswan  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(adcswan1 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(adcswan2 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(adcswan3 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(swan2serial PRIVATE ${CMAKE_BINARY_DIR}/version_mod)

    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(adcswan ${CMAKE_BINARY_DIR}/libadcswan3.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan2.a 
                                      ${CMAKE_BINARY_DIR}/libswan2serial.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan1.a 
                                      ${CMAKE_BINARY_DIR}/libswan1serial.a 
                                      ${CMAKE_BINARY_DIR}/libversion.a
                                      netcdf netcdff )
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(adcswan ${CMAKE_BINARY_DIR}/libadcswan3.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan2.a 
                                      ${CMAKE_BINARY_DIR}/libswan2serial.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan1.a 
                                      ${CMAKE_BINARY_DIR}/libswan1serial.a  
                                      ${CMAKE_BINARY_DIR}/libversion.a 
                                      netcdf netcdff XdmfCore XdmfUtils Xdmf )
    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(adcswan ${CMAKE_BINARY_DIR}/libadcswan3.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan2.a 
                                      ${CMAKE_BINARY_DIR}/libswan2serial.a 
                                      ${CMAKE_BINARY_DIR}/libadcswan1.a 
                                      ${CMAKE_BINARY_DIR}/libversion.a
                                      ${CMAKE_BINARY_DIR}/libswan1serial.a )
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)

    ADD_DEPENDENCIES(adcswan adcswan3)
    ADD_DEPENDENCIES(adcswan3 swan2serial)
    ADD_DEPENDENCIES(swan2serial adcswan2)
    ADD_DEPENDENCIES(adcswan3 adcswan2)
    ADD_DEPENDENCIES(adcswan2 adcswan1)
    ADD_DEPENDENCIES(adcswan1 swan1serial)
    ADD_DEPENDENCIES(adcswan1 version)
    
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
    SET_TARGET_PROPERTIES(swan PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/swan_mod)
    
    TARGET_INCLUDE_DIRECTORIES(swan PRIVATE ${CMAKE_BINARY_DIR}/swan_mod)
    
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

    ADD_LIBRARY(padcirc1 ${PADCIRC1_SOURCES})
    ADD_LIBRARY(padcirc2 ${PADCIRC2_SOURCES})
    ADD_LIBRARY(padcirc3 ${PADCIRC3_SOURCES})
    
    ADD_EXECUTABLE(padcirc ${PADCIRC4_SOURCES})

    TARGET_INCLUDE_DIRECTORIES(padcirc1 PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(padcirc  PRIVATE prep)
    
    SET(PADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${PADCIRC_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(padcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(padcirc ${CMAKE_BINARY_DIR}/libpadcirc2.a)
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
    
    SET_TARGET_PROPERTIES(padcirc1 PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(padcirc1 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(padcirc2 PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(padcirc2 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    SET_TARGET_PROPERTIES(padcirc3 PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(padcirc3 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(padcirc PROPERTIES Fortran_MODULE_DIRECTORY adcirc_mod)
    SET_TARGET_PROPERTIES(padcirc PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(padcirc1  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcirc2  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcirc3  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcirc  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)

    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(padcirc ${CMAKE_BINARY_DIR}/libpadcirc3.a  
                                      ${CMAKE_BINARY_DIR}/libpadcirc2.a 
                                      ${CMAKE_BINARY_DIR}/libpadcirc1.a 
                                      ${CMAKE_BINARY_DIR}/libversion.a 
                                      ${CMAKE_BINARY_DIR}/libmkdir.a
                                      netcdf netcdff )
    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(padcirc ${CMAKE_BINARY_DIR}/libpadcirc3.a  
                                      ${CMAKE_BINARY_DIR}/libpadcirc2.a 
                                      ${CMAKE_BINARY_DIR}/libpadcirc1.a 
                                      ${CMAKE_BINARY_DIR}/libversion.a 
                                      ${CMAKE_BINARY_DIR}/libmkdir.a
                                      netcdf netcdff XdmfCore XdmfUtils Xdmf )
    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        TARGET_LINK_LIBRARIES(padcirc ${CMAKE_BINARY_DIR}/libpadcirc3.a  
                                      ${CMAKE_BINARY_DIR}/libpadcirc2.a 
                                      ${CMAKE_BINARY_DIR}/libpadcirc1.a 
                                      ${CMAKE_BINARY_DIR}/libversion.a
                                      ${CMAKE_BINARY_DIR}/libmkdir.a )
    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
   

    ADD_DEPENDENCIES(padcirc padcirc3)
    ADD_DEPENDENCIES(padcirc2 padcirc1)
    ADD_DEPENDENCIES(padcirc3 padcirc2)
    ADD_DEPENDENCIES(padcirc1 version)
    ADD_DEPENDENCIES(padcirc1 mkdir)
    
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

    ADD_LIBRARY(swan1parallel STATIC ${SWAN1PARALLEL_SOURCES})
    ADD_LIBRARY(swan2parallel STATIC ${SWAN2PARALLEL_SOURCES})
    ADD_LIBRARY(padcswan1 STATIC ${PADCSWAN1_SOURCES})
    ADD_LIBRARY(padcswan2 STATIC ${PADCSWAN2_SOURCES})
    ADD_LIBRARY(padcswan3 STATIC ${PADCSWAN3_SOURCES})
    ADD_LIBRARY(padcswan4 STATIC ${PADCSWAN4_SOURCES})
    ADD_LIBRARY(padcswan5 STATIC ${PADCSWAN5_SOURCES})
   
    SET(PADCSWAN_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    
    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(padcswan ${CMAKE_BINARY_DIR}/libpadcswan5.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan4.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan3.a
                                   ${CMAKE_BINARY_DIR}/libswan2parallel.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan2.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan1.a 
                                   ${CMAKE_BINARY_DIR}/libversion.a 
                                   ${CMAKE_BINARY_DIR}/libswan1parallel.a 
                                   ${CMAKE_BINARY_DIR}/libmkdir.a
                                   netcdf netcdff )

    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
   
        SET(PADCSWAN_COMPILER_FLAGS "${PADCSWAN_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCSWAN_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcswan PROPERTIES LINK_FLAGS ${PADCSWAN_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(padcswan PRIVATE ${CMAKE_SOURCE_DIR}/src)
        
        TARGET_LINK_LIBRARIES(padcswan ${CMAKE_BINARY_DIR}/libpadcswan5.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan4.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan3.a
                                   ${CMAKE_BINARY_DIR}/libswan2parallel.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan2.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan1.a 
                                   ${CMAKE_BINARY_DIR}/libversion.a 
                                   ${CMAKE_BINARY_DIR}/libswan1parallel.a 
                                   ${CMAKE_BINARY_DIR}/libmkdir.a
                                   netcdf netcdff XdmfCore XdmfUtils Xdmf )

    ELSE(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_Xswan)

        TARGET_LINK_LIBRARIES(padcswan ${CMAKE_BINARY_DIR}/libpadcswan5.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan4.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan3.a
                                   ${CMAKE_BINARY_DIR}/libswan2parallel.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan2.a 
                                   ${CMAKE_BINARY_DIR}/libpadcswan1.a 
                                   ${CMAKE_BINARY_DIR}/libversion.a 
                                   ${CMAKE_BINARY_DIR}/libswan1parallel.a
                                   ${CMAKE_BINARY_DIR}/libmkdir.a )

    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
    
    SET_TARGET_PROPERTIES(swan1parallel PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(swan2parallel PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan1 PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan2 PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan3 PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan4 PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan5 PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan  PROPERTIES Fortran_MODULE_DIRECTORY padcswan_mod)
    SET_TARGET_PROPERTIES(padcswan1 PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan2 PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan3 PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan4 PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan5 PROPERTIES COMPILE_FLAGS ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    SET_TARGET_PROPERTIES(padcswan  PROPERTIES COMPILE_FLAGS  ${PADCSWAN_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC})
    
    TARGET_INCLUDE_DIRECTORIES(padcswan  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcswan1 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcswan2 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcswan3 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcswan4 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcswan5 PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(swan2parallel PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    
    IF( NOT ${ADDITIONAL_FLAGS_SWAN} STREQUAL "" )
        SET_TARGET_PROPERTIES(swan1parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
        SET_TARGET_PROPERTIES(swan2parallel PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS_SWAN})
    ENDIF()
    
    ADD_DEPENDENCIES(padcswan swan2parallel)
    ADD_DEPENDENCIES(padcswan padcswan3)
    ADD_DEPENDENCIES(swan2parallel padcswan5)
    ADD_DEPENDENCIES(padcswan5 padcswan4)
    ADD_DEPENDENCIES(padcswan4 padcswan2)
    ADD_DEPENDENCIES(padcswan4 padcswan3)
    ADD_DEPENDENCIES(padcswan3 padcswan2)
    ADD_DEPENDENCIES(padcswan2 padcswan1)
    ADD_DEPENDENCIES(padcswan1 swan1parallel)
    ADD_DEPENDENCIES(padcswan1 mkdir)
    
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
    
    ADD_LIBRARY(msglib ${MSGLIB_SOURCES})
    ADD_LIBRARY(punswan1 ${SWANONLY1_PARALLEL_SOURCES})
    ADD_EXECUTABLE(punswan ${SWANONLY2_PARALLEL_SOURCES})

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
    
    SET_TARGET_PROPERTIES(punswan1  PROPERTIES Fortran_MODULE_DIRECTORY punswan_mod)
    SET_TARGET_PROPERTIES(msglib    PROPERTIES Fortran_MODULE_DIRECTORY punswan_mod)
    SET_TARGET_PROPERTIES(punswan   PROPERTIES Fortran_MODULE_DIRECTORY punswan_mod)
    
    SET(MSGLIB_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${SWAN_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS}")
    SET_TARGET_PROPERTIES(msglib PROPERTIES COMPILE_FLAGS ${MSGLIB_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(punswan   PRIVATE ${CMAKE_BINARY_DIR}/punswan_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan   PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(punswan1  PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_INCLUDE_DIRECTORIES(msglib    PRIVATE ${CMAKE_BINARY_DIR}/version_mod)

    TARGET_LINK_LIBRARIES(punswan msglib punswan1 mkdir)

    ADD_DEPENDENCIES(msglib  punswan1)
    ADD_DEPENDENCIES(punswan punswan1)
    ADD_DEPENDENCIES(punswan msglib)
    ADD_DEPENDENCIES(msglib  version)
    ADD_DEPENDENCIES(msglib  mkdir)
    
    INSTALL(TARGETS punswan RUNTIME DESTINATION bin)

ENDIF(BUILD_PUNSWAN AND PERL_FOUND)
###########################################################################


###########################################################################
# ASWIP
###########################################################################
IF(BUILD_ASWIP)
    
    SET(ASWIP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ASWIP}")

    SET(ASWIP_SOURCES src/sizes.F src/global.F src/global_3dvs.F src/boundaries.F 
                      src/mesh.F src/wind.F src/owiwind.F KDTREE2/kdtree2.F src/owi_ice.F 
                      wind/vortex.F wind/aswip_1.0.3.F )


    ADD_EXECUTABLE(aswip ${ASWIP_SOURCES})


    IF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff )

    ELSEIF(ENABLE_OUTPUT_NETCDF AND ENABLE_OUTPUT_XDMF)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        SET(ASWIP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${ASWIP_LINKER_FLAGS} )
        
        TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_SOURCE_DIR}/src)
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff XdmfCore XdmfUtils Xdmf )

    ENDIF(ENABLE_OUTPUT_NETCDF AND NOT ENABLE_OUTPUT_XDMF) 

    TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_BINARY_DIR}/version_mod)
    TARGET_LINK_LIBRARIES(aswip version)

    SET_TARGET_PROPERTIES(aswip PROPERTIES COMPILE_FLAGS ${ASWIP_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(aswip PROPERTIES Fortran_MODULE_DIRECTORY aswip_mod)

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
    
    SET_TARGET_PROPERTIES(adccmp PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(p15 PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(owi22 PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(build13 PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(buildstwave23 PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(hot2asc PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(inflate PROPERTIES Fortran_MODULE_DIRECTORY util_mod)
    SET_TARGET_PROPERTIES(hstime PROPERTIES Fortran_MODULE_DIRECTORY util_mod)

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
