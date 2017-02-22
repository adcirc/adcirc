
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
            SET(NETCDF_FLAG "-DADCNETCDF")
            SET(NETCDF_LINKER_FLAG "-L${NETCDFHOME}/lib")
            SET(NETCDF_WORKING TRUE)
            LINK_DIRECTORIES(${NETCDFHOME}/lib)
            IF(NETCDF_TEST2)
                SET(NETCDF4_WORKING TRUE)
                SET(NETCDF_COMPRESSION_FLAG "-DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE")
            ELSE(NETCDF_TEST2)
                SET(NETCDF4_WORKING FALSE)
                SET(NETCDF_COMPRESSION_FLAG "")
            ENDIF(NETCDF_TEST2)
        ELSE(NETCDF_TEST1)
            MESSAGE(SEND_ERROR "The netCDF library specified is not compatible with the specified compilers. It will not be enabled. Specify a different path or disable netCDF. Ensure that you specify the same compilers to build ADCIRC as were used to build the netCDF library.")
            SET(NETCDF_WORKING FALSE)
        ENDIF(NETCDF_TEST1)

    ENDIF(${NETCDFHOME} STREQUAL "NETCDF-NOTFOUND")
ENDIF(ENABLE_OUTPUT_NETCDF)
