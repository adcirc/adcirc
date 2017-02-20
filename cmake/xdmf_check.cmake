
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
