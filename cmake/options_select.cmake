# ...Option for 8-byte real numbers (or 4-byte if turned off)
option(PRECISION_8BYTE "Use 8-byte real numbers instead of 4-byte." ON)
mark_as_advanced(PRECISION_8BYTE)

# ##############################################################################
# ...Output format options
option(ENABLE_OUTPUT_NETCDF "Turn on netCDF output options" OFF)
if(ENABLE_OUTPUT_NETCDF)
  option(ENABLE_OUTPUT_XDMF "Turn on XDMF output options" OFF)
endif(ENABLE_OUTPUT_NETCDF)
# ##############################################################################

# ##############################################################################
# ...Executables
option(BUILD_ADCIRC "Build the serial ADCIRC executable" OFF)

if(PERL_FOUND)
  option(BUILD_ADCSWAN "Build the serial SWAN+ADCIRC executable" OFF)
  option(BUILD_SWAN "Build the serial SWAN executable" OFF)
  mark_as_advanced(BUILD_SWAN)
endif(PERL_FOUND)

if(MPI_FOUND)
  option(BUILD_ADCPREP "Build the MPI parallel ADCIRC preprocessor" OFF)
  option(BUILD_PADCIRC "Build the MPI parallel ADCIRC executable" OFF)
  option(BUILD_LIBADCIRC_STATIC "Build the static library version of ADCIRC"
         OFF)
  option(BUILD_LIBADCIRC_SHARED "Build the shared library version of ADCIRC"
         OFF)
  if(PERL_FOUND)
    option(BUILD_PADCSWAN "Build the MPI parallel SWAN+ADCIRC executable" OFF)
    option(BUILD_PUNSWAN "Build the MPI parallel unstructured SWAN executable"
           OFF)
    mark_as_advanced(BUILD_PUNSWAN)
  endif(PERL_FOUND)
endif(MPI_FOUND)

option(BUILD_ASWIP "Build ASWIP (ASymmetric Wind Input Preprocessor)")
option(BUILD_UTILITIES "Build the ADCIRC utility programs" OFF)
option(ENABLE_GRIB2 "Use GRIB2API static libraries." OFF)
option(ENABLE_DATETIME "Use DATETIME static libraries." OFF)
# ##############################################################################

# ##############################################################################
# ...No need to show CXX compilers on main screen
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER CMAKE_C_COMPILER
                 CMAKE_Fortran_COMPILER)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS_RELEASE
                 CMAKE_Fortran_FLAGS_RELEASE)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_DEBUG CMAKE_C_FLAGS_DEBUG
                 CMAKE_Fortran_FLAGS_DEBUG)
# ##############################################################################

# ##############################################################################
# ...Library paths
if(ENABLE_OUTPUT_XDMF)
  if(NOT ${XDMFHOME} STREQUAL "")
    set(XDMFHOME
        ${XDMFHOME}
        CACHE STRING "XDMF home path containing lib and include")
  elseif(NOT $ENV{XDMFHOME} STREQUAL "")
    set(XDMFHOME
        $ENV{XDMFHOME}
        CACHE STRING "XDMF home path containing lib and include")
  else(NOT ${XDMFHOME} STREQUAL "")
    set(XDMFHOME
        "XDMF-NOTFOUND"
        CACHE STRING "XDMF home path containing lib and include")
  endif(NOT ${XDMFHOME} STREQUAL "")
else(ENABLE_OUTPUT_XDMF)
  unset(XDMFHOME CACHE)
endif(ENABLE_OUTPUT_XDMF)
# ##############################################################################

# ##############################################################################
# ...Additional flags
set(ADDITIONAL_FLAGS_ADCIRC
    ""
    CACHE STRING "Additional flags to compile ADCIRC with")
if(PERL_FOUND)
  set(ADDITIONAL_FLAGS_SWAN
      ""
      CACHE STRING "Additional flags to compile SWAN with")
endif(PERL_FOUND)
set(ADDITIONAL_FLAGS_ADCPREP
    ""
    CACHE STRING "Additional flags to compile ADCPREP with")
set(ADDITIONAL_FLAGS_ASWIP
    ""
    CACHE STRING "Additional flags to compile ASWIP with")
set(ADDITIONAL_FLAGS_UTLIITIES
    ""
    CACHE STRING "Additional flags for utility programs")
# ##############################################################################

# ##############################################################################
# ...Options enabled via compiler flags within the code
option(ENABLE_WARN_ELEV_DEBUG "Enable writing of the fort.69 debug file" OFF)

option(IBM "Format code for IBM based architectures" OFF)
option(SGI "Format code for SGI based architectures" OFF)
option(SUN "Format code for SUN based architectures" OFF)
option(CRAY "Format code for CRAY based architectures" OFF)
option(CRAYX1 "Format code for CRAYX1 based architectures" OFF)
mark_as_advanced(IBM SGI SUN CRAY CRAYX1)

option(DEBUG_FULL_STACK "Write the detailed stack trace during debugging" OFF)
option(DEBUG_FLUSH_MESSAGES "Do not allow caching of screen printed messages"
       OFF)
option(DEBUG_LOG_LEVEL "Force debug log level for screen messages" OFF)
option(DEBUG_ALL_TRACE
       "Write all tracing debug information to screen output during debugging"
       OFF)
option(DEBUG_GLOBALIO_TRACE
       "Write tracing debug information from the GLOBALIO module" OFF)
option(DEBUG_WRITER_TRACE
       "Write tracing debug information for writer processors" OFF)
option(DEBUG_WRITE_OUTPUT_TRACE
       "Write tracing debug information for the write output module" OFF)
option(DEBUG_WIND_TRACE "Write tracing debug information for the wind module"
       OFF)
option(DEBUG_WEIR_TRACE "Write tracing debug information for the weir module"
       OFF)
option(DEBUG_TVW_TRACE
       "Write tracing debug information for the time varying weir module" OFF)
option(DEBUG_VSMY_TRACE
       "Write tracing debug information for the 3D momentum equation module"
       OFF)
option(DEBUG_TIMESTEP_TRACE
       "Write the tracing debug information for the timestepping module" OFF)
option(DEBUG_SUBPREP_TRACE
       "Write the tracing debug information for the subdomain prep module" OFF)
option(DEBUG_SUBDOMAIN_TRACE
       "Write the tracing debug information for the subdomain module" OFF)
option(DEBUG_READ_INPUT_TRACE
       "Write the tracing debug information for the read_input module" OFF)
option(DEBUG_OWIWIND_TRACE
       "Write the tracing debug information for the OWI wind module" OFF)
option(DEBUG_NODALATTR_TRACE
       "Write the tracing debug information for the nodal attributes module"
       OFF)
option(DEBUG_NETCDF_TRACE
       "Write the tracing debug information for the netCDF modle" OFF)
option(
  DEBUG_MESSENGER_TRACE
  "Write the tracing debug information for the message passing (MPI) module"
  OFF)
option(DEBUG_MESH_TRACE
       "Write the tracing debug information for the mesh module" OFF)
option(DEBUG_HOTSTART_TRACE
       "Write the tracing debug information for the hotstart module" OFF)
option(DEBUG_GLOBAL_TRACE
       "Write the tracing debug information for the global module" OFF)
option(DEBUG_HARM_TRACE
       "Write the tracing debug information for the harmonics module" OFF)
option(DEBUG_COLDSTART_TRACE
       "Write the tracing debug information for the coldstart module" OFF)
option(DEBUG_COUPLE2SWAN_TRACE
       "Write the tracing debug information for the couple2swan module" OFF)
option(DEBUG_ADCIRC_TRACE
       "Write the tracing debug information for the main ADCIRC module" OFF)
option(DEBUG_HOLLAND
       "Write the debugging information for the symmetric Holland model" OFF)
option(DEBUG_NWS14 "Write the debugging information for NWS=14 interpolation"
       OFF)
mark_as_advanced(
  DEBUG_GLOBALIO_TRACE
  DEBUG_WRITER_TRACE
  DEBUG_WRITE_OUTPUT_TRACE
  DEBUG_WIND_TRACE
  DEBUG_WEIR_TRACE
  DEBUG_TVW_TRACE
  DEBUG_VSMY_TRACE
  DEBUG_TIMESTEP_TRACE
  DEBUG_SUBPREP_TRACE
  DEBUG_SUBDOMAIN_TRACE
  DEBUG_READ_INPUT_TRACE
  DEBUG_OWIWIND_TRACE
  DEBUG_NODALATTR_TRACE
  DEBUG_NETCDF_TRACE
  DEBUG_MESSENGER_TRACE
  DEBUG_MESH_TRACE
  DEBUG_HOTSTART_TRACE
  DEBUG_GLOBAL_TRACE
  DEBUG_HARM_TRACE
  DEBUG_COLDSTART_TRACE
  DEBUG_COUPLE2SWAN_TRACE
  DEBUG_ADCIRC_TRACE
  DEBUG_HOLLAND
  DEBUG_NWS14
  IBM)

option(
  ENABLE_POWELL
  "Force Powell wind drag to be enabled. Warning: Overrides any other options specified at run time."
  OFF)
mark_as_advanced(ENABLE_POWELL)

option(VECTOR_COMPUTER "Assume the system is a vector computer" OFF)
mark_as_advanced(VECTOR_COMPUTER)

option(ADCIRC_NOF2008
       "The fortran compiler being used does not have F2008 intrinsics" OFF)
mark_as_advanced(ADCIRC_NOF2008)
