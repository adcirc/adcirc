!-----------------------------------------------------------------------
! adcircResultCompare.F90
!   Written by Zachary Cobell
!              ARCADIS, INC.
!              Zachary.Cobell@Arcadis.com
!
! This code will comapre two sets of ADCIRC files to see if the
! solution is different between them. The solution can be compared
! across formats as well. A return code of 0 is returned if the files
! are considered to match. A return code of 1 is returned otherwise.
!-----------------------------------------------------------------------
! @jasonfleming 20181129: Sample command line for compiling this program:
! 
! gfortran -o adcircResultsComparison -I/usr/include -L/usr/lib adcircResultsComparison.F90 -lnetcdf -lnetcdff
!
! RENCI hatteras:
! ifort -o adcircResultsComparison -I/usr/share/Modules/software/CentOS-7/netcdf-Fortran/4.4.0_intel-18.0.0/include -L/usr/share/Modules/software/CentOS-7/netcdf-Fortran/4.4.0_intel-18.0.0/lib adcircResultsComparison.F90 -lnetcdf -lnetcdff
!-----------------------------------------------------------------------

        MODULE adcircCompare_module
            USE NETCDF

            INTEGER,PARAMETER   :: nNetCDFVariables = 38
            CHARACTER(200),SAVE :: netcdf_types(nNetCDFVariables)
            CHARACTER(200),SAVE :: nc_longname(nNetCDFVariables)
            CHARACTER(200),SAVE :: nc_stdname(nNetCDFVariables)
            REAL(8),PARAMETER   :: eps = EPSILON(1D0)

            PRIVATE nNetCDFVariables,netcdf_types,nc_longname,initializeNetcdf

            CONTAINS
            
            LOGICAL FUNCTION findFile(filename,returnfalse)
                CHARACTER(*),INTENT(IN) :: Filename
                LOGICAL,INTENT(IN),OPTIONAL :: returnfalse
                LOGICAL                 :: exists

                INQUIRE(FILE=TRIM(filename),EXIST=exists)

                FindFile=.TRUE.
                IF(.NOT.exists)THEN
                    IF(PRESENT(RETURNFALSE))THEN
                        IF(returnfalse)THEN
                            FindFile=.FALSE.
                            RETURN
                        ELSE
                            WRITE(*,'(3A)') "Specified file ",TRIM(filename),&
                                " does not exist."
                            STOP 1
                        ENDIF
                    ELSE
                        WRITE(*,'(3A)') "Specified file ",TRIM(filename),&
                            " does not exist."
                        STOP 1
                    ENDIF
                ENDIF
            END FUNCTION findFile
        

            INTEGER FUNCTION getFreeUnit()
                INTEGER :: I
                LOGICAL :: ISOPEN
                I = 0
                DO
                    I = I + 1
                    INQUIRE(UNIT=I,OPENED=ISOPEN)
                    IF(ISOPEN)CYCLE
                    GETFREEUNIT = I
                    EXIT
                ENDDO
                RETURN
            END FUNCTION getFreeUnit
           
            SUBROUTINE CHECK(stat,error,fatal)
                USE NETCDF
                IMPLICIT NONE
                INTEGER,INTENT(IN)             :: stat
                LOGICAL,INTENT(OUT),OPTIONAL   :: error
                LOGICAL,INTENT(IN),OPTIONAL    :: fatal
                LOGICAL                        :: fatal_local
                INTEGER,ALLOCATABLE            :: Dmy(:)

                IF(.NOT.PRESENT(FATAL))THEN
                    FATAL_LOCAL = .TRUE.
                ELSE
                    FATAL_LOCAL = FATAL
                ENDIF
                
                IF(fatal_local)THEN
                    IF(stat.NE.NF90_NOERR)THEN
                        WRITE(*,'(A,A)') "FATAL ERROR from ",TRIM(NF90_STRERROR(stat))
#ifdef EBUG
                        !...Intentional segfault to trigger stack trace
                        Dmy(1) = 1.0D0
#endif
                        STOP 1
                    ELSE
                        IF(PRESENT(error))error = .FALSE.
                    ENDIF
                ELSE
                    IF(stat.NE.NF90_NOERR)THEN
                        IF(PRESENT(error))error = .TRUE.
                    ELSE
                        IF(PRESENT(error))error = .FALSE.
                    ENDIF
                ENDIF
            END SUBROUTINE CHECK
           

            SUBROUTINE initializeNetcdf()
                IMPLICIT NONE

                NETCDF_TYPES(:)  = ""
                NETCDF_TYPES(1)  = "sigmat"       !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(2)  = "salinity"     !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(3)  = "temperature"  !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(4)  = "u-vel3D"      !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(5)  = "v-vel3D"      !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(6)  = "w-vel3D"      !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(7)  = "q20"          !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(8)  = "l"            !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(9)  = "ev"           !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(10) = "qsurfkp1"     !...Not implemented, placed here just to follow ADCIRC
                NETCDF_TYPES(11) = "zeta"
                NETCDF_TYPES(12) = "zeta_max"
                NETCDF_TYPES(13) = "u-vel"
                NETCDF_TYPES(14) = "v-vel"
                NETCDF_TYPES(15) = "vel_max"
                NETCDF_TYPES(16) = "pressure"
                NETCDF_TYPES(17) = "pressure_min"
                NETCDF_TYPES(18) = "windx"
                NETCDF_TYPES(19) = "windy"
                NETCDF_TYPES(20) = "wind_max"
                NETCDF_TYPES(21) = "radstress_x"
                NETCDF_TYPES(22) = "radstress_y"
                NETCDF_TYPES(23) = "radstress_max"
                NETCDF_TYPES(24) = "swan_HS"
                NETCDF_TYPES(25) = "swan_HS_max"
                NETCDF_TYPES(26) = "swan_DIR"
                NETCDF_TYPES(27) = "swan_DIR_max"
                NETCDF_TYPES(28) = "swan_TM01"
                NETCDF_TYPES(29) = "swan_TM01_max"
                NETCDF_TYPES(30) = "swan_TPS"
                NETCDF_TYPES(31) = "swan_TPS_max"
                NETCDF_TYPES(32) = "swan_windx"
                NETCDF_TYPES(33) = "swan_windy"
                NETCDF_TYPES(34) = "swan_wind_max"
                NETCDF_TYPES(35) = "swan_TM02"
                NETCDF_TYPES(36) = "swan_TM02_max"
                NETCDF_TYPES(37) = "swan_TMM10"
                NETCDF_TYPES(38) = "swan_TMM10_max"

                NC_LONGNAME(:)  = ""
                NC_LONGNAME(1)  = "water column vertically varying density"
                NC_LONGNAME(2)  = "water column vertically varying salinity"
                NC_LONGNAME(3)  = "water column vertically varying temperature"
                NC_LONGNAME(4)  = "water column vertically varying east/west velocity"
                NC_LONGNAME(5)  = "water column vertically varying north/south velocity"
                NC_LONGNAME(6)  = "water column vertically varying up/down velocity"
                NC_LONGNAME(7)  = "water column vertically varying turbulent kinetic energy"
                NC_LONGNAME(8)  = "water column vertically varying mixing length"
                NC_LONGNAME(9)  = "water column vertically varying eddy viscosity"
                NC_LONGNAME(10) = "sea surface temperature at the k+1 time level"
                NC_LONGNAME(11) = "water surface elevation above geoid"
                NC_LONGNAME(12) = "maximum water surface elevation above geoid"
                NC_LONGNAME(13) = "water column vertically averaged east/west velocity"
                NC_LONGNAME(14) = "water column vertically averaged north/south velocity"
                NC_LONGNAME(15) = "maximum water column vertically averaged velocity"
                NC_LONGNAME(16) = "air pressure at sea level"
                NC_LONGNAME(17) = "minimum air pressure at sea level"
                NC_LONGNAME(18) = "wind velocity in x-direction"
                NC_LONGNAME(19) = "wind velocity in y-direction"
                NC_LONGNAME(20) = "maximum wind velocity"
                NC_LONGNAME(21) = "radiation stress gradient x component"
                NC_LONGNAME(22) = "radiation stress gradient y component"
                NC_LONGNAME(23) = "maximum radiation stress gradient"
                NC_LONGNAME(24) = "significant wave height"
                NC_LONGNAME(25) = "maximum significant wave height"
                NC_LONGNAME(26) = "mean wave direction"
                NC_LONGNAME(27) = "maximum mean wave direction"
                NC_LONGNAME(28) = "mean absolute wave period"
                NC_LONGNAME(29) = "maximum TM01 mean wave period"
                NC_LONGNAME(30) = "smoothed peak period"
                NC_LONGNAME(31) = "maximum smoothed peak period"
                NC_LONGNAME(32) = "wind velocity in x-direction"
                NC_LONGNAME(33) = "wind velocity in y-direction"
                NC_LONGNAME(34) = "maximum wind stress"
                NC_LONGNAME(35) = "mean absoloute zero crossing period"
                NC_LONGNAME(36) = "maximum TM02 mean wave period"
                NC_LONGNAME(37) = "mean absolute wave period"
                NC_LONGNAME(38) = "maximum TMM10 mean wave period"

                NC_STDNAME(:)  = ""
                NC_STDNAME(1)  = "water_density_vertically_varying"
                NC_STDNAME(2)  = "water_salinity_vertically_varying"
                NC_STDNAME(3)  = "water_temperature_vertically_varying"
                NC_STDNAME(4)  = "eastward_water_velocity_vertically_varying"
                NC_STDNAME(5)  = "northward_water_velocity_vertically_varying"
                NC_STDNAME(6)  = "upward_water_velocity_vertically_varying"
                NC_STDNAME(7)  = "turbulent_kinetic_energy_vertically_varying"
                NC_STDNAME(8)  = "water_mixing_length_vertically_varying"
                NC_STDNAME(9)  = "water_eddy_viscosity_vertically_varying"
                NC_STDNAME(10) = "future sea surface temperature"
                NC_STDNAME(11) = "sea_surface_height_above_geoid"
                NC_STDNAME(12) = "maximum_sea_surface_height_above_geoid"
                NC_STDNAME(13) = "x_water_velocity_depth_averaged"
                NC_STDNAME(14) = "y_water_velocity_depth_averaged"
                NC_STDNAME(15) = "maximum_water_velocity_depth_averaged"
                NC_STDNAME(16) = "air_pressure_at_sea_level"
                NC_STDNAME(17) = "minimum_air_pressure_at_sea_level"
                NC_STDNAME(18) = "x_wind"
                NC_STDNAME(19) = "y_wind"
                NC_STDNAME(20) = "maximum_wind"
                NC_STDNAME(21) = "radiation_stress_gradient_x"
                NC_STDNAME(22) = "radiation_stress_gradient_y"
                NC_STDNAME(23) = "maximum_radiation_stress"
                NC_STDNAME(24) = "sea_surface_wave_significant_height"
                NC_STDNAME(25) = "maximum_sea_surface_wave_significant_height"
                NC_STDNAME(26) = "sea_surface_wave_to_direction"
                NC_STDNAME(27) = "maximum_sea_surface_wave_to_direction"
                NC_STDNAME(28) = "sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"
                NC_STDNAME(29) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"
                NC_STDNAME(30) = "sea_surface_wave_period_at_variance_spectral_density_maximum"
                NC_STDNAME(31) = "maximum_sea_surface_wave_period_at_variance_spectral_density_maximum"
                NC_STDNAME(32) = "x_wind"
                NC_STDNAME(33) = "y_wind"
                NC_STDNAME(34) = "maximum_wind"
                NC_STDNAME(35) = "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment"
                NC_STDNAME(36) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment"
                NC_STDNAME(37) = "sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"
                NC_STDNAME(38) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"

                RETURN

            END SUBROUTINE initializeNetcdf


            SUBROUTINE findMyNetCDFVariable(NCID)
                IMPLICIT NONE
                INTEGER,INTENT(IN)  :: NCID
                INTEGER :: I
                INTEGER :: J
                INTEGER :: NVAR
                CHARACTER(200) :: NC_NAME

                CALL CHECK(NF90_INQUIRE(NCID,NVARIABLES=NVAR))

                CALL initializeNetcdf()
                DO I = 1,NVAR
                    CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,I,&
                        NAME=NC_NAME))
                    DO J = 1,SIZE(NETCDF_TYPES)
                        IF(TRIM(ADJUSTL(NC_NAME)).EQ.&
                                TRIM(ADJUSTL(NETCDF_TYPES(J))))THEN
                            RETURN
                        ENDIF
                    ENDDO
                    IF(I.EQ.NVAR)THEN
                        WRITE(*,'(A)') &
                            "ADCIRC NetCDF Variable not found in file."
                        STOP 1
                    ENDIF
                ENDDO

            END SUBROUTINE findMyNetCDFVariable


            SUBROUTINE getNetCDFVarId(NCID,VARID1,VARID2,NCOLS,VarName1,VarName2)
                USE netcdf
                IMPLICIT NONE
                INTEGER,INTENT(IN)  :: NCID
                INTEGER,INTENT(OUT) :: VARID1
                INTEGER,INTENT(OUT) :: VARID2
                INTEGER,INTENT(OUT) :: NCOLS
                CHARACTER(*),INTENT(OUT),OPTIONAL :: VarName1,VarName2
                INTEGER             :: I
                INTEGER             :: J
                INTEGER             :: NVAR
                CHARACTER(200)      :: NC_NAME

                CALL CHECK(NF90_INQUIRE(NCID,NVARIABLES=NVAR))

                DO I = 1,NVAR
                    CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,I,&
                        NAME=NC_NAME))
                    DO J = 1,SIZE(NETCDF_TYPES)
                        IF(TRIM(NC_NAME).EQ.TRIM(NETCDF_TYPES(J)))THEN
                            VARID1 = I
                            IF(PRESENT(VarName1))THEN
                                VarName1 = NC_NAME
                            ENDIF
                            SELECT CASE(J)
                                CASE(13,18,21,33)
                                    NCOLS=2
                                    CALL CHECK(NF90_INQ_VARID(NCID,&
                                        TRIM(NETCDF_TYPES(J+1)),VARID2))
                                    IF(PRESENT(VarName2))THEN
                                        VarName2 = NETCDF_TYPES(J+1)
                                    ENDIF
                                CASE DEFAULT
                                    NCOLS=1
                                    VARID2=-1
                            END SELECT
                            RETURN
                        ENDIF
                    ENDDO
                ENDDO
                STOP 1
            END SUBROUTINE getNetCDFVarId

            SUBROUTINE showHelp()
                IMPLICIT NONE

                WRITE(*,'(A)') "Usage: adcircCompare [OPTIONS]"
                WRITE(*,'(A)') ""
                WRITE(*,'(A)') "Options:"
                WRITE(*,'(A)') "  -h,--help          Displays this help"
                WRITE(*,'(A)') "  -f1,--file1        File 1 for comparison"
                WRITE(*,'(A)') "  -f2,--file2        File 2 for comparison"
                WRITE(*,'(A)') "  -t,--tolerance     Maximum allowable difference in results"
                WRITE(*,'(A)') "  -w,--wetdry        Ignore wet/dry differences"
                WRITE(*,'(A)') "  -v,--verbose       Verbose output"
                WRITE(*,'(A)') "  -c,--continue      Continue even when output is considered to have differences"
            END SUBROUTINE showHelp


            SUBROUTINE processCommandLineArgs(file1,file2,tolerance,wetdry,verbose,cont)
                IMPLICIT NONE
                
                CHARACTER(*),INTENT(OUT)   :: file1,file2
                REAL(8),INTENT(OUT)        :: tolerance
                LOGICAL,INTENT(OUT)        :: wetdry
                LOGICAL,INTENT(OUT)        :: verbose
                LOGICAL,INTENT(OUT)        :: cont
                INTEGER                    :: iargc
                INTEGER                    :: I
                CHARACTER(2000)            :: CMD
                LOGICAL                    :: foundFile1,foundFile2,exists

                tolerance = 0.00001D0
                wetDry = .FALSE.
                verbose = .FALSE.
                cont = .FALSE.
                foundFile1 = .FALSE.
                foundFile2 = .FALSE.

                IF(iargc().EQ.0)THEN
                    WRITE(*,'(A)') "ERROR: No command line arguments specified."
                    CALL showHelp()
                    STOP 1
                ENDIF

                I = 0
                DO WHILE(I.LT.iargc())
                
                    I = I + 1
                    CALL GETARG(I,CMD)

                    SELECT CASE(TRIM(CMD))
                        CASE("-f1","--file1")
                            I = I + 1
                            CALL GETARG(I,file1)
                            foundFile1 = .TRUE.
                        CASE("-f2","--file2")
                            I = I + 1
                            CALL GETARG(I,file2)
                            foundFile2 = .TRUE.
                        CASE("-t","--tolerance")
                            I = I + 1
                            CALL GETARG(I,CMD)
                            READ(CMD,*) tolerance
                        CASE("-w","--wetdry")
                            wetdry = .TRUE.
                        CASE("-h","--help")
                            CALL showHelp()
                            STOP 0
                        CASE("-v","--verbose")
                            verbose = .TRUE.
                        CASE("-c","--continue")
                            cont = .TRUE.
                        CASE DEFAULT
                            WRITE(*,'(2A)') "ERROR: Unrecognized argument ",TRIM(CMD)
                            CALL showHelp()
                            STOP 1
                    END SELECT
                ENDDO

                IF(.NOT.foundFile1)THEN
                    WRITE(*,'(A)') "ERROR: File 1 not specified."
                    CALL showHelp()
                    STOP 1
                ENDIF

                IF(.NOT.foundFile2)THEN
                    WRITE(*,'(A)') "ERROR: File 2 not specified."
                    CALL showHelp()
                    STOP 1
                ENDIF

                exists = FindFile(file1)

                exists = FindFile(file2)

                RETURN

            END SUBROUTINE processCommandLineArgs
            
            INTEGER FUNCTION determineFileType(filename)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename

                !...Check the various output formats to determine 
                !   the type of file that is specified

                !...Check 1: netcdf
                IF(checkFormat_netcdf(filename))THEN
                    determineFileType = 3
                    RETURN
                ENDIF

                !...Check 2: xdmf
                IF(checkFormat_xdmf(filename))THEN
                    determineFileType = 7
                    RETURN
                ENDIF
                
                !...Check 3: full format ascii
                IF(checkFormat_fullFormatAscii(filename))THEN
                    determineFileType = 1
                    RETURN
                ENDIF

                !...Check 4: sparse format ascii
                IF(checkFormat_sparseFormatAscii(filename))THEN
                    determineFileType = 4
                    RETURN
                ENDIF

                determineFileType = -9999
                RETURN

            END FUNCTION determineFiletype

            LOGICAL FUNCTION checkFormat_netcdf(filename)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER                 :: ncid
                INTEGER                 :: ierr

                ierr = NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid)
                IF(ierr.EQ.NF90_NOERR)THEN
                    checkFormat_netCDF = .TRUE.
                    RETURN
                ELSE
                    checkFormat_netCDF = .FALSE.
                ENDIF
                RETURN

            END FUNCTION checkFormat_netcdf

            LOGICAL FUNCTION checkFormat_xdmf(filename)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename

                checkFormat_xdmf = .FALSE.

            END FUNCTION checkFormat_xdmf

            LOGICAL FUNCTION checkFormat_fullFormatAscii(filename)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER                 :: io
                CHARACTER(400)          :: header
                INTEGER                 :: JunkI
                REAL(8)                 :: JunkR

                io = getFreeUnit()
                OPEN(FILE=TRIM(filename),UNIT=io,ACTION="READ")
                READ(io,'(A)') header
                READ(io,'(A)') header
                READ(io,'(A)') header
                CLOSE(io)

                READ(header,*,ERR=100,END=100) JunkR,JunkI,JunkI,JunkR
                checkFormat_fullFormatAscii = .FALSE.
                RETURN
100             CONTINUE
                checkFormat_fullFormatAscii = .TRUE.
                RETURN
            END FUNCTION checkFormat_fullFormatAscii
            
            LOGICAL FUNCTION checkFormat_sparseFormatAscii(filename)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER                 :: io
                CHARACTER(400)          :: header
                INTEGER                 :: JunkI
                REAL(8)                 :: JunkR

                io = getFreeUnit()
                OPEN(FILE=TRIM(filename),UNIT=io,ACTION="READ")
                READ(io,'(A)') header
                READ(io,'(A)') header
                READ(io,'(A)') header
                CLOSE(io)

                READ(header,*,ERR=100,END=100) JunkR,JunkI,JunkI,JunkR
                checkFormat_sparseFormatAscii = .TRUE.
                RETURN
100             CONTINUE
                checkFormat_sparseFormatAscii = .FALSE.
                RETURN
            END FUNCTION checkFormat_sparseFormatAscii

            SUBROUTINE getMetadata(filename,noutfile,numsnaps,numnodes,ncol)
                IMPLICIT NONE
                INTEGER,INTENT(IN)      :: noutfile
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER,INTENT(OUT)     :: numsnaps,numnodes,ncol

                ncol = 0

                IF(noutfile.EQ.1.OR.noutfile.EQ.4)THEN
                    CALL getMetadataAscii(filename,numsnaps,numnodes,ncol)
                ELSEIF(noutfile.EQ.3)THEN
                    CALL getMetadataNetCDF(filename,numsnaps,numnodes)
                ELSE
                    numsnaps = -9999
                    numnodes = -9999
                ENDIF

                RETURN

            END SUBROUTINE getMetadata

            SUBROUTINE getMetadataAscii(filename,numsnaps,numnodes,numvalues)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER,INTENT(OUT)     :: numsnaps,numnodes,numvalues
                INTEGER                 :: io
                INTEGER                 :: junki
                REAL(8)                 :: junkr
                CHARACTER(200)          :: header
                io = getFreeUnit()
                OPEN(FILE=TRIM(filename),UNIT=io,ACTION="READ")
                READ(io,'(A)') header
                READ(io,'(A)') header
                CLOSE(io)
                READ(header,*) numsnaps,numnodes,junkR,junkI,numvalues
                RETURN
            END SUBROUTINE getMetadataAscii

            SUBROUTINE getMetadataNetCDF(filename,numsnaps,numnodes)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER,INTENT(OUT)     :: numsnaps,numnodes
                INTEGER                 :: ncid
                INTEGER                 :: dimid_time
                INTEGER                 :: dimid_node
                INTEGER                 :: ierr

                CALL CHECK(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
                CALL CHECK(NF90_INQ_DIMID(ncid,"time",dimid_time))

                ierr = NF90_INQ_DIMID(ncid,"node",dimid_node)
                IF(ierr.NE.NF90_NOERR)THEN
                    ierr = NF90_INQ_DIMID(ncid,"station",dimid_node)
                ENDIF
                CALL CHECK(NF90_INQUIRE_DIMENSION(ncid,dimid_time,LEN=numsnaps))
                CALL CHECK(NF90_INQUIRE_DIMENSION(ncid,dimid_node,LEN=numnodes))
                CALL FindMyNetCDFVariable(ncid)
                CALL CHECK(NF90_CLOSE(ncid))

                RETURN
            END SUBROUTINE getMetadataNetCDF

            SUBROUTINE openFile(filename,filetype,fileunit)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: filename
                INTEGER,INTENT(IN)      :: filetype
                INTEGER,INTENT(OUT)     :: fileunit
                CHARACTER(200)          :: header

                IF(filetype.EQ.1.OR.filetype.EQ.4)THEN
                    fileunit = getFreeUnit()
                    OPEN(FILE=TRIM(filename),UNIT=fileunit,ACTION="READ")
                    READ(fileunit,*) header
                    READ(fileunit,*) header
                ELSEIF(filetype.EQ.3)THEN
                    CALL CHECK(NF90_OPEN(TRIM(filename),NF90_NOWRITE,fileunit))
                ELSE
                    fileunit = -9999
                ENDIF

                RETURN
            END SUBROUTINE openFile

            SUBROUTINE closeFile(fileunit,filetype)
                IMPLICIT NONE
                INTEGER,INTENT(IN) :: fileunit
                INTEGER,INTENT(IN) :: filetype

                IF(filetype.EQ.1.OR.filetype.EQ.4)THEN
                    CLOSE(fileunit)
                ELSEIF(filetype.EQ.3)THEN
                    !CALL CHECK(NF90_CLOSE(fileunit))
                ELSE
                    STOP 1
                ENDIF

                RETURN
            END SUBROUTINE closeFile

            SUBROUTINE readNextSnapFullAscii(io_unit,nnodes,nvalues,nodaldata,time)
                IMPLICIT NONE
                INTEGER,INTENT(IN)    :: io_unit
                INTEGER,INTENT(IN)    :: nnodes
                INTEGER,INTENT(IN)    :: nvalues
                REAL(8),INTENT(INOUT) :: nodaldata(:,:)
                REAL(8),INTENT(OUT)   :: time
                INTEGER               :: timestep
                INTEGER               :: I,J
                INTEGER               :: junki

                READ(io_unit,*) time,timestep
                DO I = 1,nnodes
                   READ(io_unit,*) junkI,(nodalData(I,J),J=1,nvalues)
                ENDDO

                RETURN
            END SUBROUTINE readNextSnapFullAscii

            SUBROUTINE readNextSnapSparseAscii(io_unit,nnodes,nvalues,nodaldata,time)
                IMPLICIT NONE
                INTEGER,INTENT(IN)    :: io_unit
                INTEGER,INTENT(IN)    :: nnodes
                INTEGER,INTENT(IN)    :: nvalues
                REAL(8),INTENT(INOUT) :: nodaldata(:,:)
                REAL(8),INTENT(OUT)   :: time
                INTEGER               :: node
                INTEGER               :: timestep
                INTEGER               :: nnondefault
                INTEGER               :: I,J
                REAL(8)               :: defaultvalue

                READ(io_unit,*) time,timestep,nnondefault,defaultvalue
                nodaldata(:,:) = defaultvalue
                DO I = 1,nnondefault
                    READ(io_unit,*) node,(nodalData(node,J),J=1,nvalues)
                ENDDO

                RETURN
            END SUBROUTINE readNextSnapSparseAscii

            SUBROUTINE readNextSnapNetCDF(io_unit,snap,nnodes,nvalues,varid1,varid2,nodaldata)
                IMPLICIT NONE
                INTEGER,INTENT(IN)    :: io_unit
                INTEGER,INTENT(IN)    :: nnodes
                INTEGER,INTENT(IN)    :: nvalues
                INTEGER,INTENT(IN)    :: varid1,varid2
                INTEGER,INTENT(IN)    :: snap
                REAL(8),INTENT(INOUT) :: nodaldata(:,:)

                CALL CHECK(NF90_GET_VAR(io_unit,varid1,nodaldata(:,1),START=(/1,snap/),COUNT=(/nnodes,1/)))
                IF(nvalues.EQ.2)THEN
                    CALL CHECK(NF90_GET_VAR(io_unit,varid2,nodaldata(:,2),START=(/1,snap/),COUNT=(/nnodes,1/)))
                ENDIF
                
                RETURN
            END SUBROUTINE readNextSnapNetCDF

            SUBROUTINE compareData(nnodes,data1,data2,tolerance,wetdry,ndiff,avgdiff,maxdiff,nerror)
                IMPLICIT NONE
                INTEGER,INTENT(IN)  :: nnodes
                REAL(8),INTENT(IN)  :: data1(:)
                REAL(8),INTENT(IN)  :: data2(:)
                REAL(8),INTENT(IN)  :: tolerance
                LOGICAL,INTENT(IN)  :: wetdry
                INTEGER,INTENT(OUT) :: ndiff
                REAL(8),INTENT(OUT) :: avgdiff
                REAL(8),INTENT(OUT) :: maxdiff
                INTEGER,INTENT(OUT) :: nerror
                INTEGER             :: I
                REAL(8)             :: d

                ndiff = 0
                nerror = 0
                avgdiff = 0D0
                maxdiff = 0D0

                DO I = 1,nnodes
                    IF(wetdry)THEN
                        IF(data1(I).LE.-9999D0.OR.data2(I).LE.-9999D0)THEN
                            CYCLE
                        ENDIF
                    ENDIF
                    d = data1(I)-data2(I)
                    IF(ABS(d).GT.eps)THEN
                        ndiff = ndiff + 1
                    ENDIF
                    IF(ABS(d).GT.tolerance)THEN
                        nerror = nerror + 1
                    ENDIF
                    IF(ABS(d).GT.ABS(maxdiff).AND.ABS(d).GT.eps)THEN
                        maxdiff = d
                    ENDIF
                    avgdiff = avgdiff + d
                ENDDO
                avgdiff = avgdiff / nnodes

                RETURN
            END SUBROUTINE compareData

            SUBROUTINE displaySummary(nsnaps,nvalues,tolerance,wetdry,ndiff,avgdiff,maxdiff,nerror)
                IMPLICIT NONE
                INTEGER,INTENT(IN) :: nsnaps
                INTEGER,INTENT(IN) :: nvalues
                REAL(8),INTENT(IN) :: tolerance
                LOGICAL,INTENT(IN) :: wetdry
                REAL(8),INTENT(IN) :: avgdiff(:,:)
                REAL(8),INTENT(IN) :: maxdiff(:,:)
                INTEGER,INTENT(IN) :: nerror(:,:)
                INTEGER,INTENT(IN) :: ndiff(:,:)
                INTEGER            ::I

                WRITE(*,'(A)') "Summary of results file comparison"
                WRITE(*,'(A)') ""
                WRITE(*,'(A)') "|-----------------------------------------------------------------------|"
                WRITE(*,'(A)') "|  TIMESNAP  |  NUMDIFF  |  AVERAGE DIFF  |  MAXIMUM DIFF  |  NUMERROR  |"
                WRITE(*,'(A)') "|-----------------------------------------------------------------------|"
                DO I = 1,nsnaps
                    IF(nvalues.EQ.1)THEN
                        IF(ndiff(I,1).GT.0)THEN
                            WRITE(*,'(A,I9,A,I9,A,1PE12.4E3,A,1PE12.4E3,A,I9,A)') &
                               "|",I,"   | ",ndiff(I,1)," |  ",avgdiff(I,1),"  |  ",maxdiff(I,1),"  | ",&
                                nerror(I,1),"  |"
                        ENDIF
                    ELSEIF(nvalues.EQ.2)THEN
                        IF(ndiff(I,1).GT.0.OR.ndiff(I,2).GT.0)THEN
                            WRITE(*,'(A,I9,A,I9,A,1PE12.4E3,A,1PE12.4E3,A,I9,A)') &
                                "|",I,"   | ",ndiff(I,1)+ndiff(I,2)," |  ",(avgdiff(I,1)+avgdiff(I,2))/2D0,&
                                "  |  ",MAX(ABS(maxdiff(I,1)),ABS(maxdiff(I,2))),"  | ",&
                                nerror(I,1),"  |"
                        ENDIF
                    ENDIF
                ENDDO
                WRITE(*,'(A)') "|-----------------------------------------------------------------------|"
                WRITE(*,'(A)') ""
                WRITE(*,'(A,1PE12.4E3)') "  -->       Error Tolerance: ",tolerance
                IF(wetdry)THEN
                    WRITE(*,'(A)') "  --> Wet/Dry Diffs Skipped: Yes"
                ELSE
                    WRITE(*,'(A)') "  --> Wet/Dry Diffs Skipped: No"
                ENDIF

                RETURN
            END SUBROUTINE displaySummary

        END MODULE adcircCompare_module

        PROGRAM adcircResultCompare
            USE adcircCompare_module

            IMPLICIT NONE

            CHARACTER(2000)     :: file1,file2
            INTEGER             :: noutfile1,noutfile2
            INTEGER             :: nsnapfile1,nsnapfile2
            INTEGER             :: nnodefile1,nnodefile2
            INTEGER             :: nvalues1,nvalues2
            INTEGER             :: ncol1,ncol2
            INTEGER             :: io_unit1,io_unit2
            INTEGER             :: varid11,varid12,varid21,varid22
            INTEGER             :: I,J
            INTEGER,ALLOCATABLE :: nerror(:,:),ndiff(:,:)
            REAL(8),ALLOCATABLE :: maxdiff(:,:),avgdiff(:,:)
            REAL(8)             :: tolerance
            REAL(8)             :: time1,time2
            REAL(8),ALLOCATABLE :: nodaldata1(:,:),nodaldata2(:,:)
            LOGICAL             :: wetdry,verbose,cont

            !...Process the command line arguments from the user
            CALL processCommandLineArgs(file1,file2,tolerance,wetdry,verbose,cont)

            !...Determine the file type of each file
            noutfile1 = determineFileType(file1)
            noutfile2 = determineFileType(file2)

            !...Pull the number of time snaps and nodes from each file
            CALL getMetadata(file1,noutfile1,nsnapfile1,nnodefile1,ncol1)
            CALL getMetadata(file2,noutfile2,nsnapfile2,nnodefile2,ncol2)

            !...Check for the number of snaps in each file
            IF(nsnapfile1.NE.nsnapfile2)THEN
                WRITE(*,'(A)') "ERROR: There are different numbers of output snaps in the files."
                STOP 1
            ENDIF

            !...Check for the number of nodes in each file
            IF(nnodefile1.NE.nnodefile2)THEN
                WRITE(*,'(A)') "ERROR: There are different numbers of nodes in the files."
                STOP 1
            ENDIF
            
            io_unit1 = 0
            io_unit2 = 0
            CALL openFile(file1,noutfile1,io_unit1)
            CALL openFile(file2,noutfile2,io_unit2)
            
            !...For netCDF files, locate the netCDF variable ID's for later use
            IF(noutfile1.EQ.3)THEN
                CALL getNetCDFVarId(io_unit1,varid11,varid12,ncol1)
            ENDIF
            IF(noutfile2.EQ.3)THEN
                CALL getNetCDFVarId(io_unit2,varid21,varid22,ncol2)
            ENDIF

            !...Check if comparing like files
            IF(ncol1.NE.ncol2)THEN
                WRITE(*,'(A)') "ERROR: The number of values in the files is different."
                STOP 1
            ENDIF

            ALLOCATE(nodaldata1(1:nnodefile1,1:ncol1))
            ALLOCATE(nodaldata2(1:nnodefile2,1:ncol2))
            ALLOCATE(ndiff(1:nsnapfile1,1:ncol1))
            ALLOCATE(nerror(1:nsnapfile1,1:ncol1))
            ALLOCATE(avgdiff(1:nsnapfile1,1:ncol1))
            ALLOCATE(maxdiff(1:nsnapfile1,1:ncol1))

            !...Time comparison loop
            DO I = 1,nsnapfile1
                
                IF(verbose)THEN
                    WRITE(*,'(A,I6,A,I6)') "Processing snap ",I," of ",nsnapfile1
                ENDIF

                !...Read data from file 1
                IF(noutfile1.EQ.1)THEN
                    CALL readNextSnapFullAscii(io_unit1,nnodefile1,ncol1,nodaldata1,time1)
                ELSEIF(noutfile1.EQ.4)THEN
                    CALL readNextSnapSparseAscii(io_unit1,nnodefile1,ncol1,nodaldata1,time1)
                ELSEIF(noutfile1.EQ.3)THEN
                    CALL readNextSnapNetCDF(io_unit1,I,nnodefile1,ncol1,varid11,varid12,nodaldata1)
                ELSE
                    STOP 1
                ENDIF
                
                !...Read data from file 2
                IF(noutfile2.EQ.1)THEN
                    CALL readNextSnapFullAscii(io_unit2,nnodefile2,ncol2,nodaldata2,time2)
                ELSEIF(noutfile2.EQ.4)THEN
                    CALL readNextSnapSparseAscii(io_unit2,nnodefile2,ncol2,nodaldata2,time2)
                ELSEIF(noutfile1.EQ.3)THEN
                    CALL readNextSnapNetCDF(io_unit2,I,nnodefile2,ncol2,varid21,varid22,nodaldata2)
                ELSE
                    STOP 1
                ENDIF

                !...Perform comparison
                CALL compareData(nnodefile1,nodalData1(:,1),nodalData2(:,1),tolerance,wetdry,&
                                 ndiff(I,1),avgdiff(I,1),maxdiff(I,1),nerror(I,1))
                IF(ncol1.EQ.2)THEN
                    CALL compareData(nnodefile1,nodalData1(:,2),nodalData2(:,2),tolerance,wetdry,&
                                ndiff(I,2),avgdiff(I,2),maxdiff(I,2),nerror(I,2))
                ENDIF

                !...Check to stop code on any error
                IF(.NOT.cont)THEN
                    IF(ncol1.EQ.2)THEN
                        IF(nerror(I,1).GT.0.OR.nerror(I,2).GT.0)THEN
                            WRITE(*,'(A,I0,A)') "Code stopped at snap ",I,"."
                            CALL displaySummary(I,ncol1,tolerance,wetdry,ndiff,avgdiff,maxdiff,nerror)
                            CALL closeFile(io_unit1,noutfile1)
                            CALL closeFile(io_unit2,noutfile2)
                            WRITE(*,'(A)') "Results do not match within specified tolerance."
                            STOP 1
                        ENDIF
                    ELSE
                        IF(nerror(I,1).GT.0)THEN
                            WRITE(*,'(A,I0,A)') "Code stopped at snap ",I,"."
                            CALL displaySummary(I,ncol1,tolerance,wetdry,ndiff,avgdiff,maxdiff,nerror)
                            CALL closeFile(io_unit1,noutfile1)
                            CALL closeFile(io_unit2,noutfile2)
                            WRITE(*,'(A)') "Results do not match within specified tolerance."
                            STOP 1
                        ENDIF
                    ENDIF
                ENDIF

            ENDDO

            CALL closeFile(io_unit1,noutfile1)
            CALL closeFile(io_unit2,noutfile2)

            CALL displaySummary(nsnapfile1,ncol1,tolerance,wetdry,ndiff,avgdiff,maxdiff,nerror)
            DO I = 1,nsnapfile1
                IF(nerror(I,1).GT.0)THEN
                    WRITE(*,'(A)') "Results do not match within specified tolerance."
                    STOP 1
                ENDIF
                IF(ncol1.EQ.2)THEN
                    IF(nerror(I,2).GT.0)THEN
                        WRITE(*,'(A)') "Results do not match within specified tolerance."
                        STOP 1
                    ENDIF
                ENDIF
            ENDDO
            WRITE(*,'(A)') "Results match within specified tolerance."

        END PROGRAM adcircResultCompare
