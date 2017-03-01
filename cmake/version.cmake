
ADD_LIBRARY(version ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F)
IF(WIN32)
    IF(EXISTS ${CMAKE_SOURCE_DIR}/version.F)
        ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F 
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/version.F ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
            COMMENT "Generating ADCIRC version...")
    ELSEIF(EXISTS ${CMAKE_SOURCE_DIR}/version_default.F)
        ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/version_default.F ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
            COMMENT "Generating ADCIRC version...")
    ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/version.F)
ELSE(WIN32)
    ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
        COMMAND ./adcircVersion.sh ${CMAKE_SOURCE_DIR} >/dev/null
        COMMAND cp ${CMAKE_SOURCE_DIR}/version.F ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/scripts
        COMMENT "Generating ADCIRC version...")
ENDIF(WIN32)
SET_TARGET_PROPERTIES(version PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/version_mod)
SET_TARGET_PROPERTIES(version PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}") 
