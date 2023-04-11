if(ENABLE_OUTPUT_NETCDF)

  if(NOT
     "${NETCDFHOME}"
     STREQUAL
     "")
    set(NETCDF_DIR "${NETCDFHOME}")
  elseif(
    NOT
    $ENV{NETCDFHOME}
    STREQUAL
    "")
    set(NETCDF_DIR $ENV{NETCDFHOME})
  endif(
    NOT
    "${NETCDFHOME}"
    STREQUAL
    "")

  set(NETCDF_F90 "YES")
  find_package(NetCDF)

  set(NETCDF_AdditionalLibs
      ""
      CACHE STRING "Additional libraries that may be required for netCDF")

  set(netcdf3_f90_code
      "
        PROGRAM netCDF3Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID

            IERR = NF90_OPEN('test.nc',NF90_NOWRITE,NCID)

        END PROGRAM
")
  set(netcdf4_f90_code
      "
        PROGRAM netCDF4Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID
            INTEGER :: VARID

            IERR = NF90_DEF_VAR_DEFLATE(NCID,VARID,1,1,2)

        END PROGRAM
")

  if(NOT NETCDF_FOUND)
    message(SEND_ERROR "Specify the netCDF path on the following screen")
  else(NOT NETCDF_FOUND)

    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf3check.f90" "${netcdf3_f90_code}")
    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/netcdf4check.f90" "${netcdf4_f90_code}")
    try_compile(
      NETCDF_TEST1 "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf3check.f90"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDF_INCLUDE_DIRS}"
      LINK_LIBRARIES "${NETCDF_LIBRARIES}"
      LINK_LIBRARIES "${NETCDF_AdditionalLibs}"
      OUTPUT_VARIABLE LOG1)
    try_compile(
      NETCDF_TEST2 "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/netcdf4check.f90"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NETCDF_INCLUDE_DIRS}"
      LINK_LIBRARIES "${NETCDF_LIBRARIES}"
      LINK_LIBRARIES "${NETCDF_AdditionalLibs}"
      OUTPUT_VARIABLE LOG2)

    if(NETCDF_TEST1)
      set(NETCDF_WORKING TRUE)
      if(NETCDF_TEST2)
        set(NETCDF4_WORKING TRUE)
      else(NETCDF_TEST2)
        set(NETCDF4_WORKING FALSE)
      endif(NETCDF_TEST2)
    else(NETCDF_TEST1)
      message(
        SEND_ERROR
          "The netCDF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable netCDF. Ensure that you specify the same compilers to build ADCIRC as were used to build the netCDF library."
      )
      set(NETCDF_WORKING FALSE)
    endif(NETCDF_TEST1)

  endif(NOT NETCDF_FOUND)
endif(ENABLE_OUTPUT_NETCDF)
