if(ENABLE_OUTPUT_XDMF)

  find_package( HDF5 COMPONENTS C HL )

  set(xdmf_f90_code
      "
        PROGRAM XDMFCHECK
            IMPLICIT NONE
            INCLUDE 'Xdmf.f'
            INTEGER :: xdmfunit
            CALL xdmfInit(xdmfunit)
        END PROGRAM
")

  find_library(
    XDMF_LibXdmfCore
    NAMES XdmfCore
    HINTS "${XDMFHOME}/lib"
    HINTS "${XDMFHOME}/lib/x86_64-linux-gnu")
  find_library(
    XDMF_LibXdmfUtils
    NAMES XdmfUtils
    HINTS "${XDMFHOME}/lib"
    HINTS "${XDMFHOME}/lib/x86_64-linux-gnu")
  find_library(
    XDMF_LibXdmf
    NAMES Xdmf
    HINTS "${XDMFHOME}/lib"
    HINTS "${XDMFHOME}/lib/x86_64-linux-gnu")
  set(XDMF_AdditionalLibs
      ""
      CACHE STRING "Additional libraries that may be required for XDMF")

  if(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
    message(SEND_ERROR "Specify the XDMF path on the following screen")
  else(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/xdmfcheck.f90"
         "${xdmf_f90_code}")

    try_compile(
      XDMF_TEST "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/xdmfcheck.f90"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${XDMFHOME}/include"
      LINK_LIBRARIES ${XDMF_LibXdmfCore}
      LINK_LIBRARIES ${XDMF_LibXdmfUtils}
      LINK_LIBRARIES ${XDMF_LibXdmf}
      LINK_LIBRARIES ${XDMF_AdditionalLibs}
      LINK_LIBRARIES ${HDF5_LIBRARIES}
      OUTPUT_VARIABLE XDMFLOG)

    if(XDMF_TEST)
      set(XDMF_WORKING TRUE)
    else(XDMF_TEST)
      message(
        SEND_ERROR
          "The XDMF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable XDMF. Ensure that you specify the same compilers to build ADCIRC as were used to build the XDMF library."
      )
      set(XDMF_WORKING FALSE)
    endif(XDMF_TEST)
  endif(${XDMFHOME} STREQUAL "XDMF-NOTFOUND")
endif(ENABLE_OUTPUT_XDMF)
