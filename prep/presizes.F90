    MODULE MEMORY_USAGE

! This module allows monitoring of dynamic memory usage
! vjp 10/6/2006

    type MemoryDescript_t
    integer :: highmem
    integer :: currmem
    end type MemoryDescript_t

    type (MemoryDescript_t),save :: mem_descript

    CONTAINS

    subroutine memory_init( )
    mem_descript % highmem = 0
    mem_descript % currmem = 0
    end subroutine memory_init


    subroutine memory_status( )
    print *, " "
    print *, "memory currently allocated = ", mem_descript % currmem, " bytes"
    print *, "memory high water mark     = ", mem_descript % highmem, " bytes"
    print *, " "
    end subroutine memory_status


    subroutine memory_alloc( nbytes)
    integer :: nbytes
    mem_descript % currmem =  mem_descript % currmem + nbytes
    if (mem_descript % currmem > mem_descript % highmem) then
        mem_descript % highmem =  mem_descript % currmem + nbytes
    endif
    end subroutine memory_alloc


    subroutine memory_dealloc( nbytes)
    integer :: nbytes
    mem_descript % currmem =  mem_descript % currmem - nbytes
    end subroutine memory_dealloc


    END MODULE MEMORY_USAGE



!----------------------------------------------------------------------------

!                           MODULE PRESIZES

!----------------------------------------------------------------------------

!                  For use with ADCPREP Version 2.2 (  09/19/2006 )

!                     current for ADCIRC v44.17d   5/24/2004
!                 jgf Updated for ADCIRC v45.06   10/07/2005
!                 jgf Updated for ADCIRC v45.07   11/17/2005
!                 jgf Updated for ADCIRC v45.10   01/12/2006
!                 jgf Updated for ADCIRC v45.11   02/02/2006
!                 jgf Updated for ADCIRC v45.12   03/17/2006
!                 vjp Updated for ADCIRC v46.34   09/20/2006
!----------------------------------------------------------------------------

!  Program Development History
!  ---------------------------

!  45.06 jgf Added contribution dated 4/09/04 written by M. Brown:
!   - Modified SIZEUP() subroutine to pass a parameter, USE_DEFAULT, that
!     instructs SIZEUP to use default filenames. If default filename does
!     not exist, ADCPREP resorts to user-specified action.

!  45.07 jgf Changed REAL variable declarations
!   - Needed to make the precision of REALs consistent with ADCIRC to
!     prevent any possible loss of precision during pre/post processing.

!  45.10 jgf Partially updated input routines
!   - File format for 3D input files had changed in ADCIRC and the process
!     of updating ADCPREP for these changes is partially complete.

!  45.11 jgf Updated for 3D recording stations
!   - 3D recording stations are now defined by coordinates rather than
!     node numbers; routines for handling this are now complete but not
!     extensively tested.

!  45.12 jgf Added features to skip over user specified vertical grid
!     spacing as well as user specified eddy viscosity profile.

!  46.34 vjp added code for partition-only command and Shintaro's winds

!  46.44 vjp added code for prep13-only command
!----------------------------------------------------------------------------

    MODULE PRESIZES
#ifdef ADCSWAN
    USE GLOBAL, ONLY : &
    SWAN_OutputHS,SWAN_OutputDIR,SWAN_OutputTM01, &
    SWAN_OutputTPS,SWAN_OutputWIND,SWAN_OutputTM02, &
    SWAN_OutputTMM10,USE_TVW,TVW_FILE,NOUT_TVW, &
    TOUTS_TVW,TOUTF_TVW,NSPOOL_TVW,WaveWindMultiplier
#else
    USE GLOBAL, ONLY : USE_TVW,TVW_FILE,NOUT_TVW, &
    TOUTS_TVW,TOUTF_TVW,NSPOOL_TVW,WaveWindMultiplier
! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    USE GLOBAL, ONLY : WaveWindMultiplier
#endif
    USE WIND, ONLY :  WindDragLimit,DragLawString,rhoAir, &
    invertedBarometerOnElevationBoundary, &
    initWindModule

    IMPLICIT NONE

!...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
!...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef REAL4
    INTEGER, PARAMETER :: SZ = 4
    INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL8
    INTEGER, PARAMETER :: SZ = 8
    INTEGER, PARAMETER :: NBYTE=8
#endif

!...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
!...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
!   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
!   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

    INTEGER, PARAMETER ::  NPREC=7

!     variables which describe the full domain
    INTEGER :: MNPROC ! number of subdomains
    INTEGER :: MNE    ! number of elements
    INTEGER :: MNP    ! number of nodes
    INTEGER :: MNEI   ! maxnodes connected to any node
    INTEGER :: MNOPE  ! number of elevation boundary segments
    INTEGER :: MNETA  ! number of Elevation Boundary Nodes
    INTEGER :: MNBOU  ! number of Land Boundary Segments
    INTEGER :: MNVEL  ! number of Land Boundary Nodes
    INTEGER :: MNTIF  ! number of Tidal Potential Constituents
    INTEGER :: MNBFR  ! number of Periodic Elevation Bndry Forcing Constituents
    INTEGER :: MNFFR  ! number of Periodic Norm. Flow Bndry Forcing Constituents
    INTEGER :: MNSTAE ! number of Elevation Recording Stations
    INTEGER :: MNSTAV ! number of Velocity Recording Stations
    INTEGER :: MNSTAC ! number of Concentration Recording Stations
    INTEGER :: MNSTAM ! number of meteorological Recording Stations
    INTEGER :: MNWP   ! 1 If No Meteorlogical Forcing else = MNP
    INTEGER :: MNWLAT ! number of Latitudes
    INTEGER :: MNWLON ! number of Longitudes

    INTEGER :: NSTA3DD! number of 3D density recording stations
    INTEGER :: NSTA3DV! number of 3D velocity recording stations
    INTEGER :: NSTA3DT! number of 3D turbulence recording stations

! kmd48.33bc - added these variables for boundary conditions
    INTEGER :: RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP

    INTEGER :: MNPP   ! number of max nodes of any subdomain
    INTEGER :: MNEP   ! number of max elements of any subdomain
!     Model Type:
    INTEGER :: IDEN        ! determines whether prognostic or diagnostic
    LOGICAL :: C2DDI       ! .TRUE. if 2D depth integrated
    LOGICAL :: C2D_BTrans  ! .TRUE. if 2D with baroclinic transport
    LOGICAL :: C2D_PTrans  ! .TRUE. if 2D with passive transport
    LOGICAL :: CBaroclinic ! .TRUE. if baroclinic (density forcing)
    LOGICAL :: C3D         ! .TRUE. if 3D
    LOGICAL :: C3DDSS      ! .TRUE. if 3D stress formulation
    LOGICAL :: C3DVS       ! .TRUE. if 3D velocity formulation
    LOGICAL :: C3D_BTrans  ! .TRUE. if 3D with baroclinic transport
    LOGICAL :: C3D_PTrans  ! .TRUE. if 3D with passive transport
!     GWCE Lumping:
    LOGICAL :: CLUMP  ! .TRUE. only if ILUMP > 0
!     Tidal Forcing:
    LOGICAL :: CTIP   ! .TRUE. only if NTIP <> 0
!     Solver Type:
    LOGICAL :: CSOLIT ! .TRUE. only if ITITER > 0
    LOGICAL :: CSOLDIR! .TRUE. only if ITITER = 0
    LOGICAL :: CSOLDIA! .TRUE. only if ITITER < 0
! kmd - added for baroclinic portion of rivers
    LOGICAL :: BndBCRiver ! .TRUE. if 3D, baroclinc and includes rivers in the simulation
    INTEGER :: NRIVBCN    ! number of river nodes to be evaluate for salinity and temperature
    INTEGER :: NDISC
    INTEGER :: skipNETA ! avoid collision with NETA in pre_global module
    LOGICAL :: fluxBoundary = .FALSE. ! avoid collision with NFLUXF in pre_global module
! tcm v51.20.03 additions for station files outside the fort.15
    LOGICAL :: USE_ELEV_STAT_FILE    ! .TRUE. if an elevation station file exists
    LOGICAL :: USE_VEL_STAT_FILE     ! .TRUE. if a velocity station file exists
    LOGICAL :: USE_CONC_STAT_FILE    ! .TRUE. if a concentration station file exists
    LOGICAL :: USE_MET_STAT_FILE     ! .TRUE. if a met station file exists


! For Definition of Working Directory

    INTEGER,SAVE :: LNAME = 6
    CHARACTER*6,SAVE :: DIRNAME = 'PE0000'

! Logicals added for adcprep paths
    LOGICAL :: PARTITION   ! .TRUE. if only mesh partition is to be performed
    LOGICAL :: USE_DEFAULT ! .TRUE. iff fort.x to be used as input
    LOGICAL :: PREP_ALL    ! .TRUE. if all input files should be written
    LOGICAL :: PREP_15     ! .TRUE. if only RunInfo file is to be localized
    LOGICAL :: PREP_13     ! .TRUE. if only nodal attributes file is to be localized
    LOGICAL :: HOT_LOCAL    ! .TRUE. if only hotstart file is to be localized
    LOGICAL :: HOT_GLOBAL   ! .TRUE. if only hotstart files is to be globalized
    LOGICAL :: PREP_20  ! .TRUE. if user specified --prep20 cmd line option
! OGICAL :: PREP_88  ! .true. if user specified --prep88 cmd line option !commented out by tcm v51.27

    LOGICAL :: useNetCDF = .FALSE.  !jgf48.07
    LOGICAL :: useXDMF = .FALSE.    !jgf51.27 trigger metadata read from fort.15
    INTEGER, ALLOCATABLE :: NM(:,:) !jgf48.08


#if defined CSWAN || defined ADCSWAN
! jgf50.60.09: Namelist for turning SWAN output files on or off
    NAMELIST /SWANOutputControl/ SWAN_OutputHS, SWAN_OutputDIR, &
    SWAN_OutputTM01, SWAN_OutputTPS, SWAN_OutputWIND, &
    SWAN_OutputTM02, SWAN_OutputTMM10
#endif
    NAMELIST /metControl/ WindDragLimit,DragLawString,rhoAir, &
    invertedBarometerOnElevationBoundary
    NAMELIST /TVWControl/ USE_TVW,TVW_FILE,NOUT_TVW,TOUTS_TVW, &
    TOUTF_TVW,NSPOOL_TVW
    LOGICAL :: NETCDF_AVAIL ! true if netCDF was compiled in
    LOGICAL :: XDMF_AVAIL ! true if XDMF was compiled in
! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    NAMELIST /waveCoupling/ WaveWindMultiplier

!---------------------end of data declarations--------------------------------C

    CONTAINS

!-----------------------------------------------------------------------------C
!  SUBROUTINE SIZEUP
!  Determine sizes of full domain arrays for ADCIRC.                          C
!                                                                             C
!  Determine runtime option logicals
!                                                                             C
! vjp 12/08/99                                                                C
! rl  10/10/01                                                                C

!-----------------------------------------------------------------------------C

!     jgf51.21.27: Broke the sizeup routine into two pieces, one for
!     reading the fort.14 and one for reading the fort.15. This enables
!     the reading of mesh files in different formats without affecting
!     the reading of the fort.15 file.
    SUBROUTINE SIZEUP14()
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER, ALLOCATABLE :: NNEIGH(:)
    INTEGER ::  I,IL,IDUM,N,N1,N2,N3,NDIF1,NDIF2,NDIF3,NBN,NVEL, &
    IK,IBN,NBBN,NFLUXB,NFLUXI,NIBP,IBTYPE,IBN1,NETA
    CHARACTER(95) :: LINEI
    CHARACTER(1) :: CHARI(95)
    EQUIVALENCE (LINEI,CHARI(1))
    LOGICAL :: FOUND
    CHARACTER(60) :: GRID,RUNINFO
    CHARACTER(len=80) skipped !jgf46.00 data that we want to skip
    INTEGER :: ios  ! return value of i/o operation

! bell - Added for found time varying weir namelist
    LOGICAL :: found_tvw_nml
!.... tcm v51.20.04 additions for station location file
    integer ::  ios_stations
    integer :: MNSTAE2,MNSTAV2,MNSTAC2,MNSTAM2

#ifdef ADCNETCDF
    NETCDF_AVAIL = .TRUE. 
#else
    NETCDF_AVAIL = .FALSE. 
#endif

!...OPEN AND PROCESS THE UNIT 14 ADCIRC GRID FILE TO DETERMINE SIZES

    FOUND = .FALSE. 
    IF (USE_DEFAULT) THEN
        GRID='fort.14'
        INQUIRE(FILE=GRID,EXIST=FOUND)
        IF(FOUND) THEN
            GOTO 32
        ELSE
            print *, GRID, " not found"
            STOP
        ENDIF
    ELSE
        31 WRITE(*,*) 'Enter the name of the ADCIRC UNIT 14 (Grid) file:'
        READ(*,'(a60)') GRID
        INQUIRE(FILE=trim(GRID),EXIST=FOUND)
        IF(FOUND) THEN
            GOTO 32
        ELSE
            GOTO 31
        ENDIF
    ENDIF
    32 WRITE(*,1011) trim(GRID)
    OPEN(14,FILE=trim(GRID))

    33 READ(14,'(a80)') LINEI                       !SKIP OVER AGRID
    DO I=1,95
        IF(CHARI(I) /= ' ') GOTO 34
    END DO
    GOTO 33
    34 READ(14,*) MNE,MNP                            !PROCESS MNE,MNP

    ALLOCATE (NNEIGH(MNP))                        !Allocate Neighbor Table
    nbytes = 4*mnp
    call memory_alloc(nbytes)
!     jgf48.08 For netCDF
    ALLOCATE( NM(MNE,3) )
    nbytes = 4*(3*mne)
    call memory_alloc(nbytes)

    DO IL=1,MNP                                   !SKIP OVER NODES
        READ(14,*) IDUM
        NNEIGH(IL)=0
    END DO

    DO IL=1,MNE                                   !READ IN THE ELEMENT
        READ(14,*) N,IDUM,N1,N2,N3                 !CONNECTIVITY TABLE
        NNEIGH(N1)=NNEIGH(N1)+1                    !DETERMINE THE# OF NEIGHBORS
        NNEIGH(N2)=NNEIGH(N2)+1
        NNEIGH(N3)=NNEIGH(N3)+1
    ! jgf48.08 populate element table for netCDF
        NM(N,1)=N1
        NM(N,2)=N2
        NM(N,3)=N3
    ENDDO

    NETA=0                                        !PROCESS OPEN BOUNDARIES
    READ(14,*) MNOPE
    READ(14,*) MNETA

    MNEI=0                                        !PROCESS MAX# NEIGHBORS
    DO IL=1,MNOPE
        READ(14,*) NBN
        NETA=NETA+NBN
        DO IK=1,NBN
            READ(14,*) IBN
            IF (NNEIGH(IBN) /= 0) THEN
                NNEIGH(IBN)=NNEIGH(IBN)+1
                IF (NNEIGH(IBN) > MNEI) MNEI=NNEIGH(IBN)
                NNEIGH(IBN) = 0
            ENDIF
        ENDDO
    ENDDO
    NETA = MNETA
    skipNETA = MNETA
    IF(MNOPE == 0) MNOPE=1
    IF(MNETA == 0) MNETA=1

    NVEL=0                            !PROCESS LAND BOUNDARIES
    NDISC=0                           !non-zero normal discharge
    NBBN=0                            !NO. OF MAINLAND BARRIER BOUNDARY NODES
    NFLUXB=0                          !SPECIFIED MAINLAND BARRIER BC
    NIBP=0                            !NO. OF INTERNAL BARRIER BOUNDARY PAIRS
    NFLUXI=0                          !SPECIFIED INTERNAL BARRIER BC
    fluxBoundary= .FALSE.             !SPECIFIED FLUX BC

    READ(14,*) MNBOU                  !Land Boundary Segments
    READ(14,*) MNVEL                  !Land Boundary Nodes

!     jgf46.21 Added support for IBTYPE=52.
    DO IL=1,MNBOU
        READ(14,*) NBN,IBTYPE
    ! kmd - added in rivers for baroclinic simulations
        IF (ABS(IBTYPE/100) == 1) THEN
            BndBCRiver= .TRUE. 
            NRIVBCN=NRIVBCN+NBN
            IBTYPE=(ABS(IBTYPE)-100)*(IBTYPE/ABS(IBTYPE))
        END IF
    ! kmd - continue on with the normal dividing of the land boundaries
        IF((IBTYPE == 2) .OR. (IBTYPE == 12) .OR. (IBTYPE == 22) .OR. &
        (IBTYPE == 32) .OR. (IBTYPE == 52)) THEN
            fluxBoundary = .TRUE. 
            NDISC=NDISC+NBN
        ENDIF
        IF((IBTYPE == 3) .OR. (IBTYPE == 13) .OR. (IBTYPE == 23)) THEN
            NFLUXB=1
            NBBN=NBBN+NBN
        ENDIF
        IF((IBTYPE == 4) .OR. (IBTYPE == 24)) THEN
            NFLUXI=1
            NIBP=NIBP+NBN
        ENDIF
        IF((IBTYPE /= 2) .AND. (IBTYPE /= 12) .AND. (IBTYPE /= 22) .AND. &
        (IBTYPE /= 52) .AND. &
        (IBTYPE /= 3) .AND. (IBTYPE /= 13) .AND. (IBTYPE /= 23) .AND. &
        (IBTYPE /= 4) .AND. (IBTYPE /= 24)) THEN
            NVEL=NVEL+NBN
        ENDIF
        IBN1=0
        DO IK=1,NBN
            READ(14,*) IBN
            IF (NNEIGH(IBN) /= 0) THEN
                NNEIGH(IBN)=NNEIGH(IBN)+1
                IF (NNEIGH(IBN) > MNEI) MNEI=NNEIGH(IBN)
                NNEIGH(IBN) = 0
            ENDIF
            IF ((IBTYPE == 1) .OR. (IBTYPE == 11) .OR. (IBTYPE == 21)) THEN
                IF ((IK == NBN) .AND. (IBN /= IBN1)) NVEL=NVEL+1
            ENDIF
            IF (IK == 1) IBN1=IBN
        ENDDO
    ENDDO

    MNVEL=NVEL+NDISC+NBBN+2*NIBP
    IF(MNBOU == 0) MNBOU=1
    MNVEL=MNVEL+1

    DO IL=1,MNP   ! FINISH DET. MAX# NEIGHBORS
        IF(NNEIGH(IL) > MNEI) MNEI=NNEIGH(IL)
    END DO
    MNEI=MNEI+1

! if ONLY partitioning mesh using metis or a prephot
    IF (PARTITION .OR. HOT_LOCAL .OR. HOT_GLOBAL) THEN
        CLOSE(14)
        return
    ELSE
        REWIND(14)         ! performing full prep
    ENDIF

! jgf51.21.27: This subroutine rewinds unit 14 but leaves the
! file open as it will be read again in read_global.F
! (subroutine read14)
    return
    1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
!-----------------------------------------------------------------------------C
    end subroutine sizeup14
!-----------------------------------------------------------------------------C




!-----------------------------------------------------------------------------C
!  SUBROUTINE SIZEUP
!  Determine sizes of full domain arrays for ADCIRC.                          C
!                                                                             C
!  Determine runtime option logicals
!                                                                             C
! vjp 12/08/99                                                                C
! rl  10/10/01                                                                C

!-----------------------------------------------------------------------------C

! jgf51.21.27: This subroutine was created out of the second half
! of the original sizeup() subroutine to separate the reading of
! the control file from the reading of the mesh file. This
! enables different subroutines to be called for reading the mesh
! file, depending on the format.
    SUBROUTINE SIZEUP15()
    use memory_usage

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER ::  I,IL,IDUM,N,N1,N2,N3,NDIF1,NDIF2,NDIF3, &
    IM,NWP,NCOR,NTIP,NWS,NRAMP,IREFYR,NTIF,IL1,Il2,NHARF,NHAINC,IDUM2, &
    ITITER,ILUMP,IHOT &
    , NHSTAR                                                        !for netCDF
    INTEGER :: IGC  ! 3DVS vertical grid code
    INTEGER :: NFEN ! 3DVS number of nodes in the vertical grid
    INTEGER :: NRS
    INTEGER :: NCICE !tcm v49.64.01 -- added for ice
    REAL(8) RSTIMINC,DT,STATIM,REFTIM,WTIMINC
    REAL(8) CICE_TIMINC !tcm v49.64.01 -- added for ice
    REAL(SZ) GRAVITY,TAU0,RDUM,A00,B00,C00,THAS,THAF,FMV
    INTEGER :: IEVC,I3DSD,I3DSV,I3DST
    CHARACTER(95) :: LINEI
    CHARACTER(1) :: CHARI(95)
    EQUIVALENCE (LINEI,CHARI(1))
    LOGICAL :: FOUND
    CHARACTER(60) :: GRID,RUNINFO
    CHARACTER(len=80) skipped !jgf46.00 data that we want to skip

!....tcm v50.66.02 additions for time varying bathymetry
    INTEGER :: NDDT
    REAL(8) BTIMINC,BCHGTIMINC
    INTEGER :: ios_nddt
    logical :: found_tbc_nml   !flag to determine if the timebathycontrol namelist was present
    NAMELIST /TimeBathyControl/ nddt,BTIMINC,BCHGTIMINC
    INTEGER :: ios  ! return value of i/o operation
!...
!... tcm v50.79 addition for metControl namelist
    logical :: found_metCon_nml  !flag to determine if the metControl namelist was present
! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    LOGICAL :: FOUND_WC_NML
    LOGICAL :: found_tvw_nml
!.... tcm v51.20.04 additions for station location file
    integer ::  ios_stations
    integer :: MNSTAE2,MNSTAV2,MNSTAC2,MNSTAM2

#ifdef ADCNETCDF
    NETCDF_AVAIL = .TRUE. 
#else
    NETCDF_AVAIL = .FALSE. 
#endif
#ifdef ADCXDMF
    XDMF_AVAIL = .TRUE. 
#else
    XDMF_AVAIL = .FALSE. 
#endif

!     tcm v50.66.02 -- additions for time varying bathymetry
    NDDT = 0
    IOS_NDDT = 0
    FOUND_TBC_NML = .FALSE. 

!     tcm v50.79 -- addition for metControl namelist
    FOUND_metCon_NML = .FALSE. 

!...OPEN AND PROCESS THE UNIT 15 ADCIRC EXTERNAL MODE (2DDI)
!...HORIZONTAL RUN INFORMATION FILE

!--Enter, Locate, Open, and Read the ADCIRC UNIT 15 (Run Info) File

    FOUND = .FALSE. 
    IF (USE_DEFAULT) THEN
        RUNINFO='fort.15'
        INQUIRE(FILE=RUNINFO,EXIST=FOUND)
        IF(FOUND) THEN
            GOTO 132
        ELSE
            print *, RUNINFO, " not found"
            STOP
        ENDIF
    ELSE
        131 WRITE(*,*) 'Enter the name of ADCIRC UNIT 15 (Run Info) file:'
        READ(*,'(A)') RUNINFO
        INQUIRE(FILE=trim(RUNINFO),EXIST=FOUND)
        IF(FOUND) THEN
            GOTO 132
        ELSE
            GOTO 131
        ENDIF
    ENDIF
    132 WRITE(*,1011) RUNINFO
    OPEN(15,FILE=RUNINFO)

!...  tcm v50.66.02 Addtions for Time Varying Bathymetry
!...  read through the fort.15 file for the namelist (TimeBathyControl) for
!...  the time varying bathymetry.  This namelist must be at the bottom of the
!...  fort.15 file. If found, then set the appropriate values (btiminc,bchgtiminc,
!...  and nddt).  If the namelist is not there, then the time varying bathymetry
!...  will not be used.
!...
!...  After this search and read, we will close the file and then reopen it
!...  for further processing the traditional non-namelist components.
!...
    READ(UNIT = 15,NML = TimeBathyControl,IOSTAT = IOS_NDDT)

    IF (IOS_NDDT < 0) THEN
    !.....   it is possible for the namelist to be present in the file and occuring at the end
    !.....   of the file with no line breaks after the ending "\" which causes the iostat to
    !.....   return a negative value.  By checking to be sure a namelist variable was set to
    !....    a non-default value we can determine this was the case.
        IF (NDDT /= 0) THEN
        !            WRITE(*,*) 'NAMELIST PRESENT, BUT AT THE END OF FILE WITH',
        !     &                 ' NO ADVANCING CHARACTER'
            found_tbc_nml = .TRUE. 
        ELSE
        !            WRITE(*,*) 'NAMELIST NOT PRESENT'
        ENDIF
    ELSEIF (IOS_NDDT == 0) THEN
    !         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
        found_tbc_nml = .TRUE. 
    ELSE
        found_tbc_nml = .TRUE. 
    !         WRITE(*,*) 'THERE WAS A PROBLEM PROCESSING THE TimeBathyControl NAMELIST'
    !         WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCIRC NOW.'
        stop
    ENDIF
!.....REWIND THE FORT.15 FILE IN ORDER TO PROCESS THE REMAINING PIECES
    REWIND(15)

#if defined CSWAN || defined ADCSWAN
! jgf50.60.08: Add a namelist for the user to turn SWAN output
! on and off. Similar to tcm's timevaryingbathy namelist.
    READ(UNIT=15,NML=SWANOutputControl,IOSTAT=IOS)
    IF (IOS > 0) THEN
        write(*,*)"INFO: The SWANOutputControl namelist was not found."
    ENDIF
    REWIND(15)
#endif
! jgf50.60.13: Add a namelist for the user to control met forcing.
! Similar to tcm's timevaryingbathy namelist.
    call initWindModule()
    READ(UNIT=15,NML=metControl,IOSTAT=IOS)
    IF (IOS > 0) THEN
        write(*,*) "INFO: The metControl namelist was not found."
    ELSE
        found_metCon_nml = .TRUE.  !tcm v50.79 added
    ENDIF
    REWIND(15)
          
! zc: Add namelist for user to control time varying weirs
    READ(UNIT=15,NML=TVWControl,IOSTAT=IOS)
    IF(IOS > 0)THEN
        WRITE(*,'(2A)') "INFO: The time varying weir control namelist", &
        " was not found."
        found_tvw_nml = .FALSE. 
    ELSE
        WRITE(*,'(2A)') "INFO: The time varying weir control namelist", &
        " was found."
        found_tvw_nml = .TRUE. 
    ENDIF
    REWIND(15)

! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    READ(UNIT=15,NML=waveCoupling,IOSTAT=IOS)
    IF(IOS > 0)THEN
        WRITE(*,*) "INFO: The waveCoupling namelist was not found."
        FOUND_WC_NML = .FALSE. 
    ELSE
        FOUND_WC_NML = .TRUE. 
    ENDIF
    REWIND(15)



!  Initialize all runtime option logicals to false

    C2DDI  = .FALSE. 
    C3D    = .FALSE. 
    C3DDSS = .FALSE. 
    C3DVS  = .FALSE. 
    C2D_BTrans  = .FALSE. !jgf46.28 Added transport/baroclinic support
    C2D_PTrans  = .FALSE. 
    CBaroclinic = .FALSE. 
    C3D_BTrans  = .FALSE. 
    C3D_PTrans  = .FALSE. 
    CLUMP  = .FALSE. 
    CTIP   = .FALSE. 
    CSOLIT = .FALSE. 

    23 READ(15,80) LINEI                             !SKIP OVER RUNDES
    DO I=1,95
        IF(CHARI(I) /= ' ') GOTO 24
    END DO
    GOTO 23
    24 READ(15,80) LINEI                             !SKIP OVER RUNID
    DO I=1,95
        IF(CHARI(I) /= ' ') GOTO 25
    ENDDO
    GOTO 24

    25 READ(15,*) IDUM  ! SKIP NFOVER
    READ(15,*) IDUM  ! SKIP NABOUT
    READ(15,*) IDUM  ! SKIP NSCREEN
    READ(15,*) IHOT
    IF ((IHOT == 367) .OR. (IHOT == 368) .OR. &
    (IHOT == 567) .OR. (IHOT == 568)) THEN
        useNetCDF = .TRUE. 
    ENDIF
    READ(15,*) IDUM  ! SKIP ICS

    READ(15,*) IM                                 !READ IM (model type)

    SELECT CASE (IM) ! jgf46.28 added transport/baroclinic support
    CASE(0)
    C2DDI = .TRUE. 
    CASE(1)
    C3D  = .TRUE. 
    C3DVS  = .TRUE. 
    CASE(2)
!        C3D  = .TRUE.
!        C3DDSS = .TRUE.
    print *, "DSS model type not presently supported"
    stop
    CASE(10)
    C2DDI = .TRUE. 
    C2D_PTrans    = .TRUE. 
    CASE(11)
    C3D           = .TRUE. 
    C3DVS         = .TRUE. 
    C3D_PTrans    = .TRUE. 
    CASE(20)
    C2DDI         = .TRUE. 
    CBaroclinic   = .TRUE. 
    CASE(21)
    C3D           = .TRUE. 
    C3DVS         = .TRUE. 
    CBaroclinic   = .TRUE. 
    CASE(30)
    C2DDI         = .TRUE. 
    C2D_PTrans    = .TRUE. 
    CBaroclinic   = .TRUE. 
    CASE(31)
    C3D           = .TRUE. 
    C3DVS         = .TRUE. 
    C3D_PTrans    = .TRUE. 
    CBaroclinic   = .TRUE. 
    CASE DEFAULT
    IF ((IM >= 111111) .AND. (IM <= 534322)) THEN
        C2DDI = .TRUE. 
    ELSEIF ((IM >= 611111) .AND. (IM <= 634322)) THEN
        C3D  = .TRUE. 
        C3DVS  = .TRUE. 
    ELSE
        print *, "model type not supported"
        stop
    ENDIF
    END SELECT

    IDEN=0
    IF (CBaroclinic) READ(15,*) IDEN
    SELECT CASE(IDEN)
    CASE(0) ! Barotropic
! do nothing, this is valid when IDEN is read in 3D section
    CASE(1) ! 2DDI Prognostic Baroclinic ADCIRC run with SigmaT forcing
    C2D_BTrans = .TRUE. 
    CASE(-1)! 2DDI Diagnostic Baroclinic ADCIRC run with SigmaT forcing
! do nothing
    CASE(2) ! 2DDI Prognostic Baroclinic ADCIRC run with Salinity forcing
    C2D_BTrans = .TRUE. 
    CASE(-2)! 2DDI Diagnostic Baroclinic ADCIRC run with Salinity forcing
! do nothing
    CASE(3) ! 2DDI Prognostic Baroclinic ADCIRC run with Temperature forcing
    C2D_BTrans = .TRUE. 
    CASE(-3)! 2DDI Diagnostic Baroclinic ADCIRC run with Temperature forcing
! do nothing
    CASE(4) ! 2DDI Prognostic Baroclinic ADCIRC run with Salinity
!               and Temperature forcing
    C2D_BTrans = .TRUE. 
    CASE(-4)! 2DDI Diagnostic Baroclinic ADCIRC run with Salinity
!               and Temperature forcing'
! do nothing
    CASE DEFAULT
    print *, "IDEN=",IDEN," not supported"
    stop
    END SELECT

    DO IL=1,4                                     !SKIP OVER NOLIBF-NOLICAT
        READ(15,*) IDUM
    ENDDO
    READ(15,*) NWP                                !SKIP NWP
    DO IL=1,NWP
        READ(15,*) skipped                !jgf46.00 skip over nodal attributes
    ENDDO
    READ(15,*) NCOR                               !SKIP OVER NCOR
    READ(15,*) NTIP                               !READ NTIP
    IF (NTIP /= 0) CTIP = .TRUE. 
    READ(15,*) NWS                                !READ NWS
    READ(15,*) NRAMP                              !SKIP OVER NRAMP
    READ(15,*) GRAVITY                            !SKIP OVER GRAVITY
    READ(15,*) TAU0                               !SKIP OVER TAU0
    IF ( (TAU0 <= -5.d0) .AND. (TAU0 > -6.d0) ) THEN
        READ(15,*) RDUM, RDUM                      !SKIP OVER TAU0 MIN AND MAX
    ENDIF
    READ(15,*) DT                                 !SKIP OVER DT
    READ(15,*) STATIM                             !SKIP OVER STATIM
    READ(15,*) REFTIM                             !SKIP OVER REFTIM
    MNWLAT = 1
    MNWLON = 1
    MNWP=1

    NRS=0
    NCICE=0   !tcm v49.64.01 -- added for ice
!       IF(ABS(NWS/100).EQ.1) THEN ! sb46.28sb03
!          NRS=1
!          NWS=(ABS(NWS)-100)*(NWS/ABS(NWS))
!       ENDIF
! C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
!       IF(ABS(NWS/100).EQ.2) THEN
!          NRS=2
!          NWS=(ABS(NWS)-200)*(NWS/ABS(NWS))
!       ENDIF
! Casey 090302: Added NWS=3xx for coupling to SWAN.
!       IF(ABS(NWS/100).EQ.3) THEN
!          NRS=3
!          NWS=(ABS(NWS)-300)*(NWS/ABS(NWS))
!       ENDIF
!. tcm v49.46 unified the special cases above
!.....tcm v49.64.01 Additions for ice
    IF(NWS == 0) THEN
        NWS = 0
        NRS = 0
        NCICE = 0
    ELSE
        NCICE = INT(ABS(NWS)/1000)
        NRS=INT((ABS(NWS) - NCICE*1000)/100)
        NWS=INT((ABS(NWS)- NCICE*1000 - NRS*100))*INT(NWS/ABS(NWS))
    ENDIF

    IF((NWS == 0) .AND. (NRS >= 1)) READ(15,*) RSTIMINC ! sb46.28sb03
    IF((NWS == 1) .AND. (NRS >= 1)) READ(15,*) RSTIMINC ! sb46.28sb03
!     jgf46.02 added NWS=8
!     jgfdebug46.02 added NWS=45
!     jgf46.16 merged:
!     cf & cm added NWS=9: asymmetric hurricane wind model
!     rjw added NWS=19: asymmetric hurricane wind model v2.0
!     jie added NWS=20: generalized asymmetric vortex model
!     sb46.28sb01 added NWS=12: OWI format

!     tcm v49.64.01 broke nws = 2 out of the group below for additions to ice
    IF(ABS(NWS) == 2) then
        IF(NRS == 0) READ(15,*) WTIMINC
        IF(NRS >= 1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
        MNWP=MNP
    ENDIF

!.... tcm v48.4641 additions for ice
    IF((ABS(NWS) == 4) .OR. (ABS(NWS) == 5) &
     .OR. (ABS(NWS) == 45) .OR. (NWS == 8) &
     .OR. (ABS(NWS) == 15)   & !jgf50.38.05: Added NWS=15 for HWind.
     .OR. (ABS(NWS) == 20)   & !jie NWS=20
     .OR. (ABS(NWS) == 16)   & !tcm v51.06.02 added for nws=16 GFDL Met
     .OR. (ABS(NWS) == 12) .OR. (NWS == 19)) THEN
    IF((NCICE == 0) .AND. (NRS == 0)) READ(15,*) WTIMINC
    IF((NCICE == 0) .AND. ((NRS == 1) .OR. (NRS == 2) .OR. (NRS == 4))) &
    READ(15,*) WTIMINC,RSTIMINC
    IF((NCICE == 0) .AND. (NRS == 3)) THEN  !Casey 090825: Fix for NRS = 3.
        IF(ABS(NWS) == 8) THEN
            READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
        ELSE
            READ(15,*) WTIMINC,RSTIMINC
        ENDIF
    ENDIF
    IF((NCICE >= 1) .AND. (NRS == 0)) READ(15,*) WTIMINC,CICE_TIMINC
    IF((NCICE >= 1) .AND. ((NRS == 1) .OR. (NRS == 2) .OR. (NRS == 4))) &
    READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
    IF((NCICE >= 1) .AND. (NRS == 3)) THEN  !Casey 090825: Fix for NRS = 3.
        IF(ABS(NWS) == 8) THEN
            READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM, &
            RSTIMINC,CICE_TIMINC
        ELSE
            READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
        ENDIF
    ENDIF
    MNWP=MNP
ENDIF

! jgf49.0804: Read WTIMINC line if NWS is 29 (VortexOWI)
    IF ((ABS(NWS)) == 29) THEN
        SELECT CASE(NRS)
        CASE(0) ! no wave radiation stress
    ! YYYY MM DD HH24 StNum BLAdj WTIMINC pureVortex pureBackgrnd
        READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,WTIMINC,RDUM,RDUM
        CASE DEFAULT ! some kind of wave radiation stress
    ! YYYY MM DD HH24 StNum BLAdj WTIMINC RSTIMINC pureVortex pureBack
        READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM, &
        WTIMINC,RSTIMINC,RDUM,RDUM
        END SELECT
    ENDIF

    IF(NWS == 3) THEN
        READ(15,*) IREFYR                          !SKIP THE REST OF THIS LINE
        READ(15,*) MNWLAT,MNWLON                   !SKIP THE REST OF THIS LINE
        MNWP=MNP
    ENDIF
!.... tcm v49.64.01 additions for ice
    IF(ABS(NWS) == 6) THEN
        IF((NCICE == 0) .AND. (NRS == 0)) READ(15,*) WTIMINC            !SKIP OVER WTIMINC
        IF((NCICE == 0) .AND. (NRS >= 1)) READ(15,*) WTIMINC,RSTIMINC
        IF((NCICE >= 1) .AND. (NRS == 0)) READ(15,*) WTIMINC,CICE_TIMINC            !SKIP OVER WTIMINC
        IF((NCICE >= 1) .AND. (NRS >= 1)) READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
        MNWP=MNP
    ENDIF

    DO IL=1,2                                    !SKIP OVER RNDAY,DRAMP
        READ(15,*) RDUM
    ENDDO
    READ(15,*) A00,B00,C00                       !READ IN GWCE TIME WEIGHTING COEFFS
    DO IL=1,5                                    !SKIP OVER H0 - CORI
        READ(15,*) RDUM
    ENDDO
    READ(15,*) NTIF                              !PROCESS NTIF
    DO IL=1,NTIF                                 !SKIP OVER TIPOTAG & TPK, AMIGT,etc.
        26 READ(15,80) LINEI
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 27
        ENDDO
        GOTO 26
        27 READ(15,*) RDUM
    ENDDO
    MNTIF=NTIF
    IF(NTIF == 0) MNTIF=1
    READ(15,*) MNBFR                             !PROCESS MNBFR
    DO IL=1,MNBFR                                !SKIP OVER BOUNTAG, & AMIG, FF,etc.
        28 READ(15,80) LINEI
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 29
        ENDDO
        GOTO 28
        29 READ(15,*) RDUM

    ENDDO
    DO IL1=1,MNBFR
        40 READ(15,80) LINEI                          !SKIP OVER ALPHA
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 41
        END DO
        GOTO 40
        41 DO IL2=1,skipNETA
            READ(15,*) RDUM                          !SKIP OVER BOUNDARY FORCINGS
        ENDDO
    ENDDO
    IF(MNBFR == 0) MNBFR=1
    READ(15,*)  RDUM                             !SKIP OVER ANGIN
    MNFFR=0
    IF (fluxBoundary.eqv. .TRUE. ) READ(15,*) MNFFR   !# FREQ IN NORMAL FLUX B.C.

    DO IL=1,MNFFR                                !SKIP OVER BOUNTAG, & AMIG, FF,etc.
        42 READ(15,80) LINEI
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 43
        ENDDO
        GOTO 42
        43 READ(15,*) RDUM
    ENDDO
    DO IL1=1,MNFFR
        44 READ(15,80) LINEI                          !SKIP OVER ALPHA
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 45
        ENDDO
        GOTO 44
        45 DO IL2=1,NDISC
            READ(15,*) RDUM                          !SKIP OVER BOUNDARY FORCINGS
        ENDDO
    ENDDO
    IF((MNFFR == 0) .OR. (MNFFR == -1)) MNFFR=1

!..... ELEVATION STATIONS

    MNSTAE=0
    MNSTAE2=0
    READ(15,*)  IDUM                              !SKIP OVER NOUTE, TOUTSE...
    call setFormatUsage(idum)

    READ(15,*) MNSTAE                             !PROCESS MNSTAE
! tcm v51.20.03 -- added negative mnstae to signal use of external station file elev_stat.151
    if (MNSTAE < 0) then
        write(*,*) "External File Used for Elevation Station Locations"
        USE_ELEV_STAT_FILE = .TRUE. 
        MNSTAE = abs(MNSTAE)  !reset to a positive value
    ! need to open the elev station file and read the number of points
        open(unit=151,file='elev_stat.151',status='old',err=7690,iostat=ios_stations)
        read(151,*) MNSTAE2
        if (abs(MNSTAE2) /= abs(MNSTAE)) then
            MNSTAE = abs(MNSTAE2)  !reset the value to what's in the file
        endif
        7690 if (ios_stations /= 0) then
            write(*,*) "Error in Reading Elevation Station File: elev_stat.151"
            stop  ! there is a stop here
        endif
        close(151)
    ELSE
        write(*,*) "Elevation Station Locations contained in fort.15"
        USE_ELEV_STAT_FILE = .FALSE. 
        DO IL=1,MNSTAE
            READ(15,*) RDUM                            !SKIP OVER STA COORDS
        ENDDO
    endif
    IF(MNSTAE == 0) MNSTAE=1

!..... VELOCITY STATIONS

    MNSTAV=0
    MNSTAV2=0
    READ(15,*) IDUM                               !SKIP OVER NOUTV, TOUTSV...
    call setFormatUsage(idum)

    READ(15,*) MNSTAV                             !PROCESS MNSTAV
! tcm v51.20.03 -- added negative mnstav to signal use of external station file vel_stat.151
    if (MNSTAV < 0 ) then
        write(*,*) "External File Used for Velocity Station Locations"
        USE_VEL_STAT_FILE = .TRUE. 
        MNSTAV = abs(MNSTAV) !reset to a positive value
        ios_stations = 0
        open(unit=151,file='vel_stat.151',status='old',err=7691,iostat=ios_stations)
        read(151,*) MNSTAV2
        if (abs(MNSTAV2) /= abs(MNSTAV)) then
            MNSTAV = abs(MNSTAV2)  !reset the value to what's in the file
        endif
        7691 if (ios_stations /= 0) then
            write(*,*) "Error in Reading Velocity Station File: vel_stat.151"
            stop  ! there is a stop here
        endif
        close(151)
    else
        write(*,*) "Velocity Station Locations Contained in fort.15"
        USE_VEL_STAT_FILE = .FALSE. 
        DO IL=1,MNSTAV
            READ(15,*) RDUM                            !SKIP OVER STA COORDS
        ENDDO
    endif
    IF(MNSTAV == 0) MNSTAV=1

!..... CONCENTRATION STATIONS

    MNSTAC=0
    MNSTAC2=0
    IF(IM == 10) THEN
        READ(15,*) IDUM                            !SKIP OVER NOUTC, TOUTSC...
        call setFormatUsage(idum)

        READ(15,*) MNSTAC                          !PROCESS MNSTAC
    ! cm v51.20.03 -- added negative mnstac to signal use of external station file conc_stat.151
        if (MNSTAC < 0) then
            write(*,*) "External File Used for Concentration Station Locations"
            USE_CONC_STAT_FILE = .TRUE. 
            MNSTAC = abs(MNSTAC) !reset to a positive value
            ios_stations = 0
            open(unit=151,file='conc_stat.151',status='old',err=7692,iostat=ios_stations)
            read(151,*) MNSTAC2
            if (abs(MNSTAC2) /= abs(MNSTAC)) then
                MNSTAC = abs(MNSTAC2)  !reset the value to what's in the file
            endif
            7692 if (ios_stations /= 0) then
                write(*,*) "Error in Reading Concentration Station File: conc_stat.151"
                stop  ! there is a stop here
            endif
            close(151)
        else
            write(*,*) "Concentration Station Locations Contained in fort.15"
            USE_CONC_STAT_FILE = .FALSE. 
            DO IL=1,MNSTAC
                READ(15,*) RDUM                         !SKIP OVER STA COORDS
            ENDDO
        endif
    ENDIF
    IF(MNSTAC == 0) MNSTAC=1

!..... MET STATIONS

    MNSTAM=0
    MNSTAM2=0
    IF(NWS /= 0) THEN
        READ(15,*) IDUM                           !SKIP OVER NOUTM, TOUTSM...
        call setFormatUsage(idum)

        READ(15,*) MNSTAM                         !PROCESS MNSTAM
    ! cm v51.20.03 -- added negative mnstam to signal use of external station file met_stat.151
        if (MNSTAM < 0 ) then
            write(*,*) "External File Used for MET Station Locations"
            USE_MET_STAT_FILE = .TRUE. 
            MNSTAM = ABS(MNSTAM)  !reset to a positive number
            ios_stations = 0
            open(unit=151,file='met_stat.151',status='old', &
            err=7693,iostat=ios_stations)
            read(151,*) MNSTAM2
            if (abs(MNSTAM2) /= abs(MNSTAM)) then
                MNSTAM = abs(MNSTAM2)  !reset the value to what's in the file
            endif
            7693 if (ios_stations /= 0) then
                write(*,*) "Error in Reading MET Station File: met_stat.151"
                stop  ! there is a stop here
            endif
            close(151)
        else
            write(*,*) "MET Station Locations Contained in fort.15"
            USE_MET_STAT_FILE = .FALSE. 
            DO IL=1,MNSTAM
                READ(15,*) RDUM                        !SKIP OVER STA COORDS
            END DO
        endif
    ENDIF
    IF(MNSTAM == 0) MNSTAM=1

    READ(15,*) IDUM                            !SKIP OVER NOUTGE, TOUTSGE...
    call setFormatUsage(idum)

    READ(15,*) IDUM                            !SKIP OVER NOUTGV, TOUTSGV...
    call setFormatUsage(idum)
          
    IF (IM == 10) THEN
        READ(15,*) IDUM               !SKIP OVER NOUTGC, TOUTSGC...
        call setFormatUsage(idum)
    ENDIF
    IF (NWS /= 0) THEN
        READ(15,*) IDUM               !SKIP OVER NOUTGW, TOUTSGW...
        call setFormatUsage(idum)
    ENDIF
    READ(15,*) NHARF                           !PROCESS MNHARF
    DO IL1=1,NHARF
        47 READ(15,80) LINEI                       !SKIP OVER HAFNAM
        DO I=1,95
            IF(CHARI(I) /= ' ') GOTO 48
        ENDDO
        GOTO 47
        48 READ(15,*) RDUM,RDUM,RDUM             !SKIP OVER HAFREQ,HAFF,HAFACE
    ENDDO

    READ(15,*) THAS,THAF,NHAINC,FMV               !READ THAS,...FMV
    READ(15,*) IDUM,IDUM,IDUM,IDUM                !SKIP OVER NHASE,NHASV,...
    READ(15,*) NHSTAR,IDUM                          !SKIP OVER NHSTAR,NHSINC
    PRINT *, "NHSTAR = ", NHSTAR
    IF ((ABS(NHSTAR) == 3) .OR. (ABS(NHSTAR) == 367) &
     .OR. (ABS(NHSTAR) == 368) .OR. &
    (ABS(NHSTAR) == 5) .OR. (ABS(NHSTAR) == 567) &
     .OR. (ABS(NHSTAR) == 568) ) THEN
        useNetCDF = .TRUE.                          !for netCDF
    ENDIF
    IF ( (useNetCDF.eqv. .TRUE. ) .AND. (NETCDF_AVAIL.eqv. .FALSE. ) ) THEN
        WRITE(*,*)'ERROR: NetCDF input or output files were specified.'
        WRITE(*,*) 'but adcprep was not compiled with NetCDF support.'
        WRITE(*,*) 'Please recompile adcprep with NetCDF libraries.'
        STOP
    ENDIF
    IF ( (useXDMF.eqv. .TRUE. ) .AND. (XDMF_AVAIL.eqv. .FALSE. ) ) THEN
        WRITE(*,*)'ERROR: XDMF output files were specified.'
        WRITE(*,*) 'but adcprep was not compiled with XDMF support.'
        WRITE(*,*) 'Please recompile adcprep with XDMF libraries.'
        STOP
    ENDIF

!...THIS SECTION TO LUMP THE GWCE MATRIX
! jp 11/30/99 made lumping a compile time option

!      PRINT *, "Made it to LUMP definition"
#ifdef LUMP
    CLUMP = .TRUE. 
    ILUMP=1
#else
    CLUMP = .FALSE. 
    ILUMP=0
#endif

    READ(15,*) ITITER,IDUM,RDUM,IDUM2             !READ SOLVER TYPE
    CSOLIT = .TRUE. 


!--Read in 3D info

    IF(C3DVS) THEN  !3DVS
    !     jgf45.10 removed IDIAG
        READ(15,*) IDUM                               !Skip IDEN
        READ(15,*) IDUM,RDUM                          !Skip ISLIP,KP
        READ(15,*) RDUM,RDUM                          !Skip Z0S,Z0B
        READ(15,*) RDUM,RDUM,RDUM                     !Skip ALP1,ALP2,ALP3
    !     jgf45.12 add code to handle user specified vertical spacing.
        READ(15,*) IGC,NFEN
        IF (IGC == 0) THEN
            DO I=1,NFEN
                READ(15,*) RDUM  !skip over vertical spacing
            ENDDO
        ENDIF
    !     jgf45.12 add code to handle user specified vertical eddy viscosity.
        READ(15,*) IEVC,RDUM,RDUM !Process IEVC
        IF (IEVC == 0) THEN
            DO I=1,NFEN
                READ(15,*) RDUM  !skip over vertical eddy viscosity
            ENDDO
        ENDIF
        IF((IEVC == 50) .OR. (IEVC == 51)) READ(15,*) RDUM,RDUM  !Skip THETA1,THETA2
        READ(15,*) I3DSD,RDUM,RDUM,IDUM               !Process I3DSD
        READ(15,*) NSTA3DD                            !Process NSTA3DD
        IF(NSTA3DD /= 0) THEN !kmd : changed to match 2D options
            DO I=1,NSTA3DD
                READ(15,*) RDUM,RDUM !Skip density stations
            END DO
        ENDIF
        READ(15,*) I3DSV,RDUM,RDUM,IDUM               !Process I3DSV
        READ(15,*) NSTA3DV                            !Process NSTA3DV
        IF (NSTA3DV /= 0) THEN  !kmd : changed to match 2D options
            DO I=1,NSTA3DV
                READ(15,*) RDUM,RDUM !Skip velocity stations
            END DO
        ENDIF
        READ(15,*) I3DST,RDUM,RDUM,IDUM !Process I3DST
        READ(15,*) NSTA3DT                            !Process NSTA3DT
        IF (NSTA3DT /= 0) THEN  !kmd : changed to match 2D options
            DO I=1,NSTA3DT
                READ(15,*) RDUM,RDUM !Skip turbulence stations
            ENDDO
        ENDIF
        READ(15,*) IDUM,RDUM,RDUM,IDUM !Skip 3D global density output
        READ(15,*) IDUM,RDUM,RDUM,IDUM                !Skip 3D global velocity output
        READ(15,*) IDUM,RDUM,RDUM,IDUM                !Skip 3D global turbulence output

    ! md48.33bc - added this information for the boundary conditions
        IF (CBAROCLINIC) THEN
            READ(15,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP    !Process RES_BC_FLAG
            IF ((RES_BC_FLAG < 0) .OR. (RES_BC_FLAG == 1)) THEN
                READ(15,*) RDUM
                READ(15,*) RDUM
            ELSE IF ((RES_BC_FLAG == 2) .OR. (RES_BC_FLAG == 3)) THEN
                READ(15,*) RDUM,RDUM
                READ(15,*) RDUM,RDUM
                IF (BCFLAG_TEMP /= 0) THEN
                    READ(15,*) RDUM, RDUM
                END IF
            ELSE IF (RES_BC_FLAG == 4) THEN
                READ(15,*) RDUM,RDUM,RDUM
                READ(15,*) RDUM,RDUM,RDUM
                IF (BCFLAG_TEMP /= 0) THEN
                    READ(15,*) RDUM, RDUM
                END IF
            END IF
        END IF

    !         PRINT *, "Made it to sponge information"
        IF (CBAROCLINIC) THEN
            READ(15,*) RDUM
        END IF


    !     ELSEIF(C3DDSS) THEN  !3DDSS
    
    ENDIF

    REWIND(15)                                    !FINISHED WITH UNIT 15 FILE


!...
    WRITE(*,3000) MNPROC,MNE,MNP,MNEI,MNOPE,MNETA, &
    MNBOU,MNVEL,MNTIF,MNBFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNWP, &
    MNWLAT,MNWLON,MNFFR

    IF(NWS == 0) WRITE(*,3011)
    IF(NWS == 1) WRITE(*,3012)
    IF(ABS(NWS) == 2) WRITE(*,3013)
    IF(NWS == 3) WRITE(*,3014)
    IF(ABS(NWS) == 4) WRITE(*,3015)
    IF(ABS(NWS) == 45) WRITE(*,3015) !jgfdebug46.02
    IF(ABS(NWS) == 5) WRITE(*,3016)
    IF(ABS(NWS) == 6) WRITE(*,3017)
    IF(ABS(NWS) == 8) WRITE(*,3217)
    IF(NRS >= 1) WRITE(*,3018) ! sb46.28sb03
    IF(NWS == 10) WRITE(*,3019)
    IF(NWS == 11) WRITE(*,3020)
!     sb46.28sb01 added NWS=12: OWI format
    IF(ABS(NWS) == 12) WRITE(*,3023)
    IF(ABS(NWS) == 19) WRITE(*,3219)
    IF(ABS(NWS) == 20) WRITE(*,3220)
    IF(ABS(NWS) == 16) WRITE(*,3026)
    IF(NCICE >= 1) WRITE(*,3221)  !tcm v49.64.01
    IF(ABS(NWS) == 29) WRITE(*,3223)
    IF((NHARF == 0) .OR. (FMV == 0.)) WRITE(*,3021)
    IF((NHARF >= 1) .AND. (FMV /= 0.)) WRITE(*,3022)
    IF(ILUMP == 0) WRITE(*,3031)
    IF(ILUMP == 1) WRITE(*,3032)
    IF(IM == 0) WRITE(*,3101)
    IF(IM == 10) WRITE(*,3109)
    IF(IM == 1) WRITE(*,3102)
    IF(IM == 2) WRITE(*,3103)
    IF(ITITER == 0) WRITE(*,3104)
    IF(ITITER > 0) WRITE(*,3105)
    IF(ITITER < 0) WRITE(*,3106)
    IF(USE_ELEV_STAT_FILE) WRITE(*,3180)  !tcm v51.20.03
    IF(USE_VEL_STAT_FILE) WRITE(*,3181)   !tcm v51.20.03
    IF(USE_MET_STAT_FILE) WRITE(*,3182)   !tcm v51.20.03
    IF(USE_CONC_STAT_FILE) WRITE(*,3183)  !tcm v51.20.03
    WRITE(*,3108)

    3000 FORMAT(' *****************************************************',/, &
    ' *   Based on input and information extracted from   *',/, &
    ' *   the ADCIRC UNIT 14 and 15 (grid and horiz run   *',/, &
    ' *   info) files the following paramter values will  *',/, &
    ' *   be set:                                         *',/, &
    ' *                                                   *',/, &
    ' *       MNPROC = ',I5,'                             *',/, &
    ' *       MNE = ',I8,1X,'     MNP = ',I8,1X,'         *',/, &
    ' *       MNEI = ',I7,'                               *',/, &
    ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/, &
    ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/, &
    ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/, &
    ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/, &
    ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/, &
    ' *       MNWP = ',I7,'                               *',/, &
    ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/, &
    ' *       MNFFR = ',I6,3X,'                           *',/, &
    ' *                                                   *')
    3011 FORMAT(' *   Also, NO wind forcing will be used,             *')
    3012 FORMAT(' *   Also, NWS=1 meteorological forcing is used,     *')
    3013 FORMAT(' *   Also, NWS=+-2 meteorological forcing is used,   *')
    3014 FORMAT(' *   Also, NWS=3 meteorological forcing is used,     *')
    3015 FORMAT(' *   Also, NWS=+-4 meteorological forcing is used,   *')
    3016 FORMAT(' *   Also, NWS=+-5 meteorological forcing is used,   *')
    3017 FORMAT(' *   Also, NWS=+-6 meteorological forcing is used,   *')
    3217 FORMAT(' *   Also, NWS=+-8 Holland wind forcing is used,     *')
    3219 FORMAT(' *   Also, NWS=+-19 Asymmetric Holland wind v2.0     *',/, &
    ' *                 forcing is used,                  *')
    3220 FORMAT(' *   Also, NWS=+-20 Generalized Asym. Model is used  *')
    3026 FORMAT(' *   Also, NWS=+-16 GFDL Met Data is used            *')
    3221 FORMAT(' *   Also, ABS(NWS)>=1000 ice concentration are used,*')
    3223 FORMAT(' *   Also, NWS=+-29 VortexOWI met forcing is used,   *')
    3018 FORMAT(' *   Also, ABS(NWS)>=100 wave stress forcing is used,*')
    3019 FORMAT(' *   Also, AVN wind & pressure forcing will be used, *')
    3020 FORMAT(' *   Also, ETA wind & pressure forcing will be used, *')
    3023 FORMAT(' *   Also, NWS=+-12 meteorological forcing is used,  *')
    3021 FORMAT(' *   means and variance calculation will NOT be made,*')
    3022 FORMAT(' *   means and variance calculation will be made,    *')
    3031 FORMAT(' *   the GWCE matrix will be left in consistent form *')
    3032 FORMAT(' *   the GWCE matrix will be LUMPED                  *')
    3101 FORMAT(' *   the model will be set up for a 2DDI run,        *')
    3109 FORMAT(' *   the model will be set up for a 2DDI run + transp*')
    3102 FORMAT(' *   the model will be set up for a 3D-VS run,       *')
    3103 FORMAT(' *   the model will be set up for a 3D-DSS run,      *')
    3104 FORMAT(' *   and the direct band solver will be used.        *')
    3105 FORMAT(' *   and an iterative solver will be used            *')
    3106 FORMAT(' *   and no external solver will be used             *')
    3180 FORMAT(' *   An external elevation station file is used      *')  !tcm v51.20.03
    3181 FORMAT(' *   An external velocity station file is used       *')  !tcm v51.20.03
    3182 FORMAT(' *   An external met. station file is used           *')  !tcm v51.20.03
    3183 FORMAT(' *   An external concentration station file is used  *') !tcm v51.20.03
    3108 FORMAT(' *****************************************************',/)

    60 FORMAT(A60)
    80 FORMAT(A95)
    180 FORMAT(95A1)
    1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
    1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)


    RETURN
!-----------------------------------------------------------------
    END SUBROUTINE SIZEUP15
!-----------------------------------------------------------------


!-----------------------------------------------------------------
!  S U B R O U T I N E   S E T   F O R M A T   U S A G E
!-----------------------------------------------------------------
! If netCDF or XDMF file formats are used, the reading of extra
! metadata from the fort.15 file is triggered.
!-----------------------------------------------------------------
    subroutine setFormatUsage(idum)
    use sizes, only : NETCDF3, NETCDF4, XDMF
    implicit none
    integer, intent(in) :: idum
          
    select case(abs(idum))
    case(NETCDF3,NETCDF4)
    useNetCDF = .TRUE. 
    case(XDMF)
    useXDMF = .TRUE. 
    case default
! don't trigger the reading of metadata in fort.15
    end select

!-----------------------------------------------------------------
    end subroutine setFormatUsage
!-----------------------------------------------------------------


!-----------------------------------------------------------------
!-----------------------------------------------------------------
    END MODULE PRESIZES
!-----------------------------------------------------------------      !-----------------------------------------------------------------
