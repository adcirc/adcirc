MACRO(addCompilerFlags TARGET)

    SET(LOCAL_COMPILER_FLAGS "${ARGN} ${ADCIRC_OPTION_FLAGS} ${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${MACHINE_FLAG} ${ARGN}")

    addLibVersion(${TARGET})
    addNetCDF(${TARGET})
    addXDMF(${TARGET})

    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/${TARGET})
    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})

    SET(LOCAL_COMPILER_FLAGS "")
    
ENDMACRO(addCompilerFlags)


MACRO(addCompilerFlagsSwan TARGET)

    SET(LOCAL_COMPILER_FLAGS "${ARGN}")

    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/${TARGET})
    IF(NOT "${LOCAL_COMPILER_FLAGS}" STREQUAL "")
        SET_TARGET_PROPERTIES(${TARGET} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})
    ENDIF(NOT "${LOCAL_COMPILER_FLAGS}" STREQUAL "")

    SET(LOCAL_COMPILER_FLAGS "")
    
ENDMACRO(addCompilerFlagsSwan)


MACRO(addNetCDF TARGET)
    IF(NETCDF_WORKING)
        TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE ADCNETCDF)
        IF(NETCDF4_WORKING)
            TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE HAVE_NETCDF4)
            TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE NETCDF_CAN_DEFLATE)
        ENDIF(NETCDF4_WORKING)
        TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${NETCDF_INCLUDE_DIRS})
        TARGET_LINK_LIBRARIES(${TARGET} ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})
    ENDIF(NETCDF_WORKING)
ENDMACRO(addNetCDF)


MACRO(addXDMF TARGET)
    IF(XDMF_WORKING)
        TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/src)
        TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${XDMFHOME}/include)
        TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE ADCXDMF)
        TARGET_LINK_LIBRARIES(${TARGET} ${XDMF_LibXdmfCore} ${XDMF_LibXdmfUtils} ${XDMF_LibXdmf} ${XDMF_AdditionalLibs})
    ENDIF(XDMF_WORKING)
ENDMACRO(addXDMF)


MACRO(addLibVersion TARGET)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_LINK_LIBRARIES(${TARGET} version)
    ADD_DEPENDENCIES(${TARGET} version)
ENDMACRO(addLibVersion)


MACRO(addLibMkdir TARGET)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/prep)
    TARGET_LINK_LIBRARIES(${TARGET} mkdir)
    ADD_DEPENDENCIES(${TARGET} mkdir)
ENDMACRO(addLibMkdir)


MACRO(addLibMetis TARGET)
    TARGET_LINK_LIBRARIES(${TARGET} metis)
    ADD_DEPENDENCIES(${TARGET} metis)
ENDMACRO(addLibMetis)


MACRO(addMPI TARGET)
   TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE CMPI)
   TARGET_LINK_LIBRARIES(${TARGET} ${MPI_Fortran_LIBRARIES})
ENDMACRO(addMPI)


MACRO(swanConfigureAdcswan)
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )
    ENDIF(WIN32)
ENDMACRO(swanConfigureAdcswan)


MACRO(swanConfigurePadcswan)
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swan_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel SWAN Sources..."
        )
    ENDIF(WIN32)
ENDMACRO(swanConfigurePadcswan)


MACRO(swanConfigureSerial)
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES}
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY_SERIAL_SOURCES} 
            COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN stand alone Sources..."
        )
    ENDIF(WIN32)
ENDMACRO(swanConfigureSerial)


MACRO(swanConfigureParallel)
    IF(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND if not exist \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\" mkdir \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\"
            COMMAND move /y *.f \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
            COMMAND move /y *.f90 \"${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Serial SWAN Sources..."
        )  
    ELSE(WIN32)
        ADD_CUSTOM_COMMAND( OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
            COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source
            COMMAND mv *.f *.f90 ${CMAKE_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/swan
            COMMENT "Generating Parallel unSWAN Sources..."
        )
    ENDIF(WIN32)
ENDMACRO(swanConfigureParallel)
