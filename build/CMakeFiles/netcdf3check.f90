
        PROGRAM netCDF3Test
            USE NETCDF
            IMPLICIT NONE

            INTEGER :: IERR
            INTEGER :: NCID

            IERR = NF90_OPEN('test.nc',NF90_NOWRITE,NCID)

        END PROGRAM
