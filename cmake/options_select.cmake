
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
    OPTION(BUILD_LIBADC "Build the library version of ADCIRC" OFF)
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
IF(ENABLE_OUTPUT_XDMF)
    IF(NOT ${XDMFHOME} STREQUAL "")
        SET(XDMFHOME ${XDMFHOME} CACHE STRING "XDMF home path containing lib and include")
    ELSEIF(NOT $ENV{XDMFHOME} STREQUAL "") 
        SET(XDMFHOME $ENV{XDMFHOME} CACHE STRING "XDMF home path containing lib and include")
    ELSE(NOT ${XDMFHOME} STREQUAL "")
        SET(XDMFHOME "XDMF-NOTFOUND" CACHE STRING "XDMF home path containing lib and include")
    ENDIF(NOT ${XDMFHOME} STREQUAL "")
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
