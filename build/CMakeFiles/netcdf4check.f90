
        PROGRAM netCDF4Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID
            INTEGER :: VARID

            IERR = NF90_DEF_VAR_DEFLATE(NCID,VARID,1,1,2)

        END PROGRAM
