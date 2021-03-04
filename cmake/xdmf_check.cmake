
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

    FIND_LIBRARY(XDMF_LibXdmfCore  NAMES XdmfCore  HINTS "${XDMFHOME}/lib" HINTS "${XDMFHOME}/lib/x86_64-linux-gnu")
    FIND_LIBRARY(XDMF_LibXdmfUtils NAMES XdmfUtils HINTS "${XDMFHOME}/lib" HINTS "${XDMFHOME}/lib/x86_64-linux-gnu")
    FIND_LIBRARY(XDMF_LibXdmf      NAMES Xdmf      HINTS "${XDMFHOME}/lib" HINTS "${XDMFHOME}/lib/x86_64-linux-gnu") 
    SET(XDMF_AdditionalLibs "" CACHE STRING "Additional libraries that may be required for XDMF")

    IF(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
        MESSAGE(SEND_ERROR "Specify the XDMF path on the following screen")
    ELSE(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
        FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" "${xdmf_f90_code}")

        TRY_COMPILE(XDMF_TEST "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/xdmfcheck.f90" CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${XDMFHOME}/include" LINK_LIBRARIES ${XDMF_LibXdmfCore} LINK_LIBRARIES ${XDMF_LibXdmfUtils} LINK_LIBRARIES ${XDMF_LibXdmf} LINK_LIBRARIES ${XDMF_AdditionalLibs} OUTPUT_VARIABLE XDMFLOG)

        IF(XDMF_TEST)
            SET(XDMF_WORKING TRUE)
            LINK_DIRECTORIES(${XDMFHOME}/lib)
        ELSE(XDMF_TEST)
            MESSAGE(SEND_ERROR "The XDMF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable XDMF. Ensure that you specify the same compilers to build ADCIRC as were used to build the XDMF library.")
            SET(XDMF_WORKING FALSE)
        ENDIF(XDMF_TEST)
    ENDIF(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
ENDIF(ENABLE_OUTPUT_XDMF)
