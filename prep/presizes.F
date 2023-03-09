      MODULE MEMORY_USAGE

      ! This module allows monitoring of dynamic memory usage
      ! vjp 10/6/2006

      type MemoryDescript_t
        integer(8) :: highmem
        integer(8) :: currmem
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


      subroutine memory_alloc(nbytes)
      integer(8),intent(in) :: nbytes
      mem_descript % currmem =  mem_descript % currmem + nbytes
      if (mem_descript % currmem > mem_descript % highmem) then
         mem_descript % highmem =  mem_descript % currmem + nbytes
      endif
      end subroutine memory_alloc


      subroutine memory_dealloc( nbytes)
      integer(8),intent(in) :: nbytes
      mem_descript % currmem =  mem_descript % currmem - nbytes
      end subroutine memory_dealloc


      END MODULE MEMORY_USAGE



C----------------------------------------------------------------------------
C
C                           MODULE PRESIZES
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 2.2 (  09/19/2006 )
C
C                     current for ADCIRC v44.17d   5/24/2004
C                 jgf Updated for ADCIRC v45.06   10/07/2005
C                 jgf Updated for ADCIRC v45.07   11/17/2005
C                 jgf Updated for ADCIRC v45.10   01/12/2006
C                 jgf Updated for ADCIRC v45.11   02/02/2006
C                 jgf Updated for ADCIRC v45.12   03/17/2006
C                 vjp Updated for ADCIRC v46.34   09/20/2006
C----------------------------------------------------------------------------
C
C  Program Development History
C  ---------------------------
C
C  45.06 jgf Added contribution dated 4/09/04 written by M. Brown:
C   - Modified SIZEUP() subroutine to pass a parameter, USE_DEFAULT, that
C     instructs SIZEUP to use default filenames. If default filename does
C     not exist, ADCPREP resorts to user-specified action.
C
C  45.07 jgf Changed REAL variable declarations
C   - Needed to make the precision of REALs consistent with ADCIRC to
C     prevent any possible loss of precision during pre/post processing.
C
C  45.10 jgf Partially updated input routines
C   - File format for 3D input files had changed in ADCIRC and the process
C     of updating ADCPREP for these changes is partially complete.
C
C  45.11 jgf Updated for 3D recording stations
C   - 3D recording stations are now defined by coordinates rather than
C     node numbers; routines for handling this are now complete but not
C     extensively tested.
C
C  45.12 jgf Added features to skip over user specified vertical grid
C     spacing as well as user specified eddy viscosity profile.
C
C  46.34 vjp added code for partition-only command and Shintaro's winds
C
C  46.44 vjp added code for prep13-only command
C----------------------------------------------------------------------------

      MODULE PRESIZES
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      USE CONSTANTS, ONLY: waveWindMultiplier
      USE GLOBAL, ONLY :
     &    USE_TVW,TVW_FILE,NOUT_TVW,TOUTS_TVW,TOUTF_TVW,NSPOOL_TVW,
     &    outputWindDrag,slim,windlim,directvelWD,useHF
#if defined CSWAN || defined ADCSWAN
     &   ,SWAN_OutputHS,SWAN_OutputDIR,SWAN_OutputTM01,
     &    SWAN_OutputTPS,SWAN_OutputWIND,SWAN_OutputTM02,
     &    SWAN_OutputTMM10
#endif
      USE WIND, ONLY :  WindDragLimit,DragLawString,rhoAir,
     &                  initWindModule
      USE OWIWIND, ONLY: nPowellSearchDomains
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      USE OWIWIND_NETCDF, ONLY: NWS13ColdStartString,
     &                          NWS13File,
     &                          NWS13GroupForPowell,
     &                          NWS13WindMultiplier
#endif
C
      IMPLICIT NONE
C
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
C...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef REAL4
      INTEGER, PARAMETER :: SZ = 4
      INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL8
      INTEGER, PARAMETER :: SZ = 8
      INTEGER, PARAMETER :: NBYTE=8
#endif

C...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
C...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
C   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
C   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

      INTEGER, PARAMETER ::  NPREC=7
C
C     variables which describe the full domain
      INTEGER MNPROC ! number of subdomains
      INTEGER MNE    ! number of elements
      INTEGER MNP    ! number of nodes
      INTEGER MNEI   ! maxnodes connected to any node
      INTEGER MNOPE  ! number of elevation boundary segments
      INTEGER MNETA  ! number of Elevation Boundary Nodes
      INTEGER MNBOU  ! number of Land Boundary Segments
      INTEGER MNVEL  ! number of Land Boundary Nodes
      INTEGER MNTIF  ! number of Tidal Potential Constituents
      INTEGER MNBFR  ! number of Periodic Elevation Bndry Forcing Constituents
      INTEGER MNFFR  ! number of Periodic Norm. Flow Bndry Forcing Constituents
      INTEGER MNSTAE ! number of Elevation Recording Stations
      INTEGER MNSTAV ! number of Velocity Recording Stations
      INTEGER MNSTAC ! number of Concentration Recording Stations
      INTEGER MNSTAM ! number of meteorological Recording Stations
      INTEGER MNWP   ! 1 If No Meteorlogical Forcing else = MNP
      INTEGER MNWLAT ! number of Latitudes
      INTEGER MNWLON ! number of Longitudes

      INTEGER NSTA3DD! number of 3D density recording stations
      INTEGER NSTA3DV! number of 3D velocity recording stations
      INTEGER NSTA3DT! number of 3D turbulence recording stations

! kmd48.33bc - added these variables for boundary conditions
      INTEGER RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP

      INTEGER MNPP   ! number of max nodes of any subdomain
      INTEGER MNEP   ! number of max elements of any subdomain
C     Model Type:
      INTEGER IDEN        ! determines whether prognostic or diagnostic
      LOGICAL C2DDI       ! .TRUE. if 2D depth integrated
      LOGICAL C2D_BTrans  ! .True. if 2D with baroclinic transport
      LOGICAL C2D_PTrans  ! .True. if 2D with passive transport
      LOGICAL CBaroclinic ! .True. if baroclinic (density forcing)
      LOGICAL C3D         ! .TRUE. if 3D
      LOGICAL C3DDSS      ! .TRUE. if 3D stress formulation
      LOGICAL C3DVS       ! .TRUE. if 3D velocity formulation
      LOGICAL C3D_BTrans  ! .true. if 3D with baroclinic transport
      LOGICAL C3D_PTrans  ! .true. if 3D with passive transport
C     GWCE Lumping:
      LOGICAL CLUMP  ! .TRUE. only if ILUMP > 0
C     Tidal Forcing:
      LOGICAL CTIP   ! .TRUE. only if NTIP <> 0
C     Solver Type:
      LOGICAL CSOLIT ! .TRUE. only if ITITER > 0
      LOGICAL CSOLDIR! .TRUE. only if ITITER = 0
      LOGICAL CSOLDIA! .TRUE. only if ITITER < 0
! kmd - added for baroclinic portion of rivers
      LOGICAL BndBCRiver ! .TRUE. if 3D, baroclinc and includes rivers in the simulation
      INTEGER NRIVBCN    ! number of river nodes to be evaluate for salinity and temperature
      INTEGER :: NDISC
      INTEGER :: skipNETA ! avoid collision with NETA in pre_global module
      LOGICAL :: fluxBoundary = .false. ! avoid collision with NFLUXF in pre_global module
! tcm v51.20.03 additions for station files outside the fort.15
      LOGICAL USE_ELEV_STAT_FILE    !.true. if an elevation station file exists
      LOGICAL USE_VEL_STAT_FILE     !.true. if a velocity station file exists
      LOGICAL USE_CONC_STAT_FILE    !.true. if a concentration station file exists
      LOGICAL USE_MET_STAT_FILE     !.true. if a met station file exists
      LOGICAL :: invertedBarometerOnElevationBoundary
C
C For Definition of Working Directory
C
      INTEGER,SAVE :: LNAME = 6
      CHARACTER*6,SAVE :: DIRNAME = 'PE0000'

C Logicals added for adcprep paths
      LOGICAL PARTITION   ! .true. if only mesh partition is to be performed
      LOGICAL USE_DEFAULT ! .true. iff fort.x to be used as input
      LOGICAL PREP_ALL    ! .true. if all input files should be written
      LOGICAL PREP_15     ! .true. if only RunInfo file is to be localized
      LOGICAL PREP_13     ! .true. if only nodal attributes file is to be localized
      LOGICAL HOT_LOCAL    ! .true. if only hotstart file is to be localized
      LOGICAL HOT_GLOBAL   ! .true. if only hotstart files is to be globalized
      LOGICAL :: PREP_20  ! .true. if user specified --prep20 cmd line option
      !LOGICAL :: PREP_88  ! .true. if user specified --prep88 cmd line option !commented out by tcm v51.27

      LOGICAL :: useNetCDF = .false.  !jgf48.07
      LOGICAL :: useXDMF = .false.    !jgf51.27 trigger metadata read from fort.15
      INTEGER, ALLOCATABLE :: NM(:,:) !jgf48.08


#if defined CSWAN || defined ADCSWAN
      ! jgf50.60.09: Namelist for turning SWAN output files on or off
      NAMELIST /SWANOutputControl/ SWAN_OutputHS, SWAN_OutputDIR,
     &    SWAN_OutputTM01, SWAN_OutputTPS, SWAN_OutputWIND,
     &    SWAN_OutputTM02, SWAN_OutputTMM10
#endif
      NAMELIST /metControl/ WindDragLimit,DragLawString,outputWindDrag,rhoAir,
     &                       invertedBarometerOnElevationBoundary,nPowellSearchDomains

      NAMELIST /TVWControl/ USE_TVW,TVW_FILE,NOUT_TVW,TOUTS_TVW,
     &    TOUTF_TVW,NSPOOL_TVW
      LOGICAL :: NETCDF_AVAIL ! true if netCDF was compiled in
      LOGICAL :: XDMF_AVAIL ! true if XDMF was compiled in
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      NAMELIST /waveCoupling/ WaveWindMultiplier
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      NAMELIST /owiWindNetcdf/ NWS13File,NWS13ColdStartString,
     &    NWS13WindMultiplier,NWS13GroupForPowell
#endif
      !
      ! @jasonfleming: reduce memory usage if there are no weir boundaries
      ! with cross barrier pipes
      logical :: meshHasWeirWithPipes
      ! Warn Elev control namelist
c.... !wjp added for warnElev namelist
      REAL(SZ) :: WarnElev = 20.0  !
      REAL(SZ) :: ErrorElev = 1000.0  !
      LOGICAL :: WarnElevDump = .FALSE.  
      INTEGER :: WarnElevDumpLimit = 50  !
      LOGICAL :: WARN_ELEV_FLAG = .FALSE.   
      LOGICAL :: FOUND_WARNELEV_NML = .false.  !
      NAMELIST /WarnElevControl/ WarnElev, WarnElevDump,
     &  WarnElevDumpLimit, ErrorElev
      ! Ali Dispersion control namelist
      REAL(SZ) :: Bd = 0.23394_SZ  ! exponent of depth (default)
      REAL(SZ) :: Ad = 0.0050189_SZ ! coefficient of depth (default)
      REAL(SZ) :: Cs = 1500.0_SZ ! Speed of sound in water
      LOGICAL :: CAliDisp= .false.  !
      LOGICAL :: FOUND_ALIDISP_NML = .false.  !
      NAMELIST /AliDispersionControl/ CAliDisp, Cs, Ad, Bd
      ! density control namelist
      LOGICAL :: FOUND_DENSITY_NML = .false.  !
      character(len=100) :: densityRunType = 'none' ! diagnostic or prognostic
      character(len=100) :: densityForcingType = 'SigmaT' ! sigmaT, salinity, etc
      character(len=100) :: densityFileName = 'fort.11.nc' ! filename
      integer :: densityTimeIterator = 1 ! BC3D time step iteration skipper
      NAMELIST /densityControl/ densityRunType, densityFileName, 
     &  densityTimeIterator, densityForcingType
      ! itfric control namelist
      LOGICAL :: FOUND_ITFRIC_NML = .FALSE. !

C---------------------end of data declarations--------------------------------C

      CONTAINS
C
C-----------------------------------------------------------------------------C
C  SUBROUTINE SIZEUP
C  Determine sizes of full domain arrays for ADCIRC.                          C
C                                                                             C
C  Determine runtime option logicals
C                                                                             C
C vjp 12/08/99                                                                C
C rl  10/10/01                                                                C
C
C-----------------------------------------------------------------------------C
C
C     jgf51.21.27: Broke the sizeup routine into two pieces, one for
C     reading the fort.14 and one for reading the fort.15. This enables
C     the reading of mesh files in different formats without affecting
C     the reading of the fort.15 file.
      SUBROUTINE SIZEUP14()
      use memory_usage
      IMPLICIT NONE
      integer(8) :: nbytes = 0
      INTEGER, ALLOCATABLE :: NNEIGH(:)
      INTEGER  I,IL,IDUM,N,N1,N2,N3,NDIF1,NDIF2,NDIF3,NBN,NVEL,
     &  IK,IBN,NBBN,NFLUXB,NFLUXI,NFLUXIPipe,NIBP,NIBPPipe,IBTYPE,IBN1,NETA
      CHARACTER*95 LINEI
      CHARACTER*1 CHARI(95)
      EQUIVALENCE (LINEI,CHARI(1))
      LOGICAL FOUND
      CHARACTER*60 GRID,RUNINFO
      CHARACTER(len=80) skipped !jgf46.00 data that we want to skip
      INTEGER :: ios  ! return value of i/o operation

Cobell - Added for found time varying weir namelist
      LOGICAL :: found_tvw_nml
C.... tcm v51.20.04 additions for station location file
      integer ::  ios_stations
      integer :: MNSTAE2,MNSTAV2,MNSTAC2,MNSTAM2

#ifdef ADCNETCDF
      NETCDF_AVAIL = .true.
#else
      NETCDF_AVAIL = .false.
#endif
      meshHasWeirWithPipes = .false.
C
C...OPEN AND PROCESS THE UNIT 14 ADCIRC GRID FILE TO DETERMINE SIZES
C
      FOUND = .FALSE.
      IF (USE_DEFAULT) THEN
         GRID='fort.14'
         INQUIRE(FILE=GRID,EXIST=FOUND)
         IF(.NOT.FOUND) THEN
            print *, GRID, " not found"
            CALL EXIT(1)
         ENDIF
      ELSE
         DO WHILE(.NOT.FOUND)
            WRITE(*,*) 
     &        'Enter the name of the ADCIRC UNIT 14 (Grid) file:'
            READ(*,'(a60)') GRID
            INQUIRE(FILE=trim(GRID),EXIST=FOUND)
         ENDDO
      ENDIF
      WRITE(*,1011) trim(GRID)
      OPEN(14,FILE=trim(GRID))

      L14: DO WHILE(.TRUE.)
         READ(14,'(a80)') LINEI                       !SKIP OVER AGRID
         DO I=1,95
            IF(CHARI(I).NE.' ') EXIT L14
         END DO
      ENDDO L14
      READ(14,*) MNE,MNP                            !PROCESS MNE,MNP
C
      ALLOCATE (NNEIGH(MNP))                        !Allocate Neighbor Table
      nbytes = 4*mnp
      call memory_alloc(nbytes)
C     jgf48.08 For netCDF
      ALLOCATE( NM(MNE,3) )
      nbytes = 4*(3*mne)
      call memory_alloc(nbytes)
C
      DO IL=1,MNP                                   !SKIP OVER NODES
         READ(14,*) IDUM
         NNEIGH(IL)=0
      END DO
C
      DO IL=1,MNE                                   !READ IN THE ELEMENT
         READ(14,*) N,IDUM,N1,N2,N3                 !CONNECTIVITY TABLE
         NNEIGH(N1)=NNEIGH(N1)+1                    !DETERMINE THE # OF NEIGHBORS
         NNEIGH(N2)=NNEIGH(N2)+1
         NNEIGH(N3)=NNEIGH(N3)+1
         !jgf48.08 populate element table for netCDF
         NM(N,1)=N1
         NM(N,2)=N2
         NM(N,3)=N3
      ENDDO
C
      NETA=0                                        !PROCESS OPEN BOUNDARIES
      READ(14,*) MNOPE
      READ(14,*) MNETA
C
      MNEI=0                                        !PROCESS MAX # NEIGHBORS
      DO IL=1,MNOPE
         READ(14,*) NBN
         NETA=NETA+NBN
         DO IK=1,NBN
            READ(14,*) IBN
            IF (NNEIGH(IBN).NE.0) THEN
              NNEIGH(IBN)=NNEIGH(IBN)+1
              IF (NNEIGH(IBN).GT.MNEI) MNEI=NNEIGH(IBN)
                NNEIGH(IBN) = 0
              ENDIF
         ENDDO
      ENDDO
      NETA = MNETA
      skipNETA = MNETA
      IF(MNOPE.EQ.0) MNOPE=1
      IF(MNETA.EQ.0) MNETA=1
C
      NVEL=0                            !PROCESS LAND BOUNDARIES
      NDISC=0                           !non-zero normal discharge
      NBBN=0                            !NO. OF MAINLAND BARRIER BOUNDARY NODES
      NFLUXB=0                          !SPECIFIED MAINLAND BARRIER BC
      NIBP=0                            !NO. OF INTERNAL BARRIER BOUNDARY PAIRS
      NIBPPipe=0                        !number of internal barrier bounary (with pipes) pairs
      NFLUXI=0                          !SPECIFIED INTERNAL BARRIER BC
      NFLUXIPipe=0                      !internal barrier boundary bc with pipes
      fluxBoundary=.false.             !SPECIFIED FLUX BC
C
      READ(14,*) MNBOU                  !Land Boundary Segments
      READ(14,*) MNVEL                  !Land Boundary Nodes
C
C     jgf46.21 Added support for IBTYPE=52.
      DO IL=1,MNBOU
         READ(14,*) NBN,IBTYPE
! kmd - added in rivers for baroclinic simulations
         IF (ABS(IBTYPE/100).EQ.1) THEN
            BndBCRiver=.TRUE.
            NRIVBCN=NRIVBCN+NBN
            IBTYPE=(ABS(IBTYPE)-100)*(IBTYPE/ABS(IBTYPE))
         END IF
! kmd - continue on with the normal dividing of the land boundaries
         select case(ibtype)
         case(2,12,22,32,52)
            fluxBoundary = .true.
            NDISC=NDISC+NBN
         case(3,13,23)
            NFLUXB=1
            NBBN=NBBN+NBN
         case(4,24,64)
            NFLUXI=1
            NIBP=NIBP+NBN
         case(5,25)
            meshHasWeirWithPipes = .false.
            NFLUXIPipe=1
            NIBPPipe=NIBPPipe+NBN
         case default  ! mainland and island
            NVEL=NVEL+NBN
         end select
         IBN1=0
         DO IK=1,NBN
            READ(14,*) IBN
            IF (NNEIGH(IBN).NE.0) THEN
              NNEIGH(IBN)=NNEIGH(IBN)+1
              IF (NNEIGH(IBN).GT.MNEI) MNEI=NNEIGH(IBN)
              NNEIGH(IBN) = 0
            ENDIF
            IF ((IBTYPE.EQ.1).OR.(IBTYPE.EQ.11).OR.(IBTYPE.EQ.21)) THEN
               IF ((IK.EQ.NBN).AND.(IBN.NE.IBN1)) NVEL=NVEL+1
            ENDIF
            IF (IK.EQ.1) IBN1=IBN
         ENDDO
      ENDDO
C
      MNVEL=NVEL+NDISC+NBBN+2*NIBP+2*NIBPPipe
      IF(MNBOU.EQ.0) MNBOU=1
      MNVEL=MNVEL+1
C
      DO IL=1,MNP   ! FINISH DET. MAX # NEIGHBORS
         IF(NNEIGH(IL).GT.MNEI) MNEI=NNEIGH(IL)
      END DO
      MNEI=MNEI+1
C
C if ONLY partitioning mesh using metis or a prephot
      IF (PARTITION .or. HOT_LOCAL .or. HOT_GLOBAL) THEN
         CLOSE(14)
         return
      ELSE
         REWIND(14)         ! performing full prep
      ENDIF
      !
      ! jgf51.21.27: This subroutine rewinds unit 14 but leaves the
      ! file open as it will be read again in read_global.F
      ! (subroutine read14)
      return
1011  FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
C-----------------------------------------------------------------------------C
      end subroutine sizeup14
C-----------------------------------------------------------------------------C




C-----------------------------------------------------------------------------C
C  SUBROUTINE SIZEUP
C  Determine sizes of full domain arrays for ADCIRC.                          C
C                                                                             C
C  Determine runtime option logicals
C                                                                             C
C vjp 12/08/99                                                                C
C rl  10/10/01                                                                C
C
C-----------------------------------------------------------------------------C
C
      ! jgf51.21.27: This subroutine was created out of the second half
      ! of the original sizeup() subroutine to separate the reading of
      ! the control file from the reading of the mesh file. This
      ! enables different subroutines to be called for reading the mesh
      ! file, depending on the format.
      SUBROUTINE SIZEUP15()
      use memory_usage

      IMPLICIT NONE
      integer(8) :: nbytes = 0
      INTEGER  I,IL,IDUM,N,N1,N2,N3,NDIF1,NDIF2,NDIF3,
     &  IM,NWP,NCOR,NTIP,NWS,NRAMP,IREFYR,NTIF,IL1,Il2,NHARF,NHAINC,IDUM2,
     &  ITITER,ILUMP,IHOT
     &  , NHSTAR                                                        !for netCDF
      INTEGER IGC  ! 3DVS vertical grid code
      INTEGER NFEN ! 3DVS number of nodes in the vertical grid
      INTEGER NRS
      INTEGER NCICE !tcm v49.64.01 -- added for ice
      INTEGER :: NHASE,NHASV,NHAGE,NHAGV
      REAL(8) RSTIMINC,DT,STATIM,REFTIM,WTIMINC
      REAL(8) CICE_TIMINC !tcm v49.64.01 -- added for ice
      REAL(SZ) GRAVITY,TAU0,RDUM,A00,B00,C00,THAS,THAF,FMV
      INTEGER IEVC,I3DSD,I3DSV,I3DST
      CHARACTER*95 LINEI
      CHARACTER*1 CHARI(95)
      EQUIVALENCE (LINEI,CHARI(1))
      LOGICAL FOUND
      CHARACTER*60 GRID,RUNINFO
      CHARACTER(len=80) skipped !jgf46.00 data that we want to skip

c....tcm v50.66.02 additions for time varying bathymetry
      INTEGER NDDT
      REAL(8) BTIMINC,BCHGTIMINC
      INTEGER :: ios_nddt
      logical :: found_tbc_nml   !flag to determine if the timebathycontrol namelist was present
      NAMELIST /TimeBathyControl/ nddt,BTIMINC,BCHGTIMINC
      INTEGER :: ios  ! return value of i/o operation
c...
c... tcm v50.79 addition for metControl namelist
      logical :: found_metCon_nml  !flag to determine if the metControl namelist was present
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      LOGICAL :: FOUND_WC_NML
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      LOGICAL :: FOUND_OWINC_NML
#endif
      LOGICAL :: found_tvw_nml
      LOGICAL :: found_WarnElev_nml = .false.
      LOGICAL :: found_AliDisp_nml = .false.
      LOGICAL :: found_density_nml = .false.
C.... tcm v51.20.04 additions for station location file
      integer ::  ios_stations
      integer :: MNSTAE2,MNSTAV2,MNSTAC2,MNSTAM2

C.... tcm v52.30.01 additions for Smagorinski Eddy Viscosity Namelist
      integer :: ios_smag = 0
      REAL(8) :: SMAG_UPPER_LIM = 100.d0  ! m^2/s (default upper limit)
      REAL(8) :: SMAG_LOWER_LIM = 1.0d-8  ! m^2/s (default lower limit)
      LOGICAL :: SMAG_COMP_FLAG = .FALSE. ! BY DEFAULT DO NOT USE SMAG
      LOGICAL :: FOUND_Smag_NML = .false. ! flag to determine if the Smag_Control namelist was present
      NAMELIST /Smag_Control/ SMAG_UPPER_LIM,SMAG_LOWER_LIM,SMAG_COMP_FLAG

#ifdef ADCNETCDF
      NETCDF_AVAIL = .true.
#else
      NETCDF_AVAIL = .false.
#endif
#ifdef ADCXDMF
      XDMF_AVAIL = .true.
#else
      XDMF_AVAIL = .false.
#endif
C
C     tcm v50.66.02 -- additions for time varying bathymetry
      NDDT = 0
      IOS_NDDT = 0
      FOUND_TBC_NML = .FALSE.
C
C     tcm v50.79 -- addition for metControl namelist
      FOUND_metCon_NML = .FALSE.
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      FOUND_OWINC_NML = .FALSE.
#endif
C
C...OPEN AND PROCESS THE UNIT 15 ADCIRC EXTERNAL MODE (2DDI)
C...HORIZONTAL RUN INFORMATION FILE
C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 15 (Run Info) File
C
      FOUND = .FALSE.
      IF (USE_DEFAULT) THEN
         RUNINFO='fort.15'
         INQUIRE(FILE=RUNINFO,EXIST=FOUND)
         IF(.NOT.FOUND) THEN
            print *, RUNINFO, " not found"
            CALL EXIT(1)
         ENDIF
      ELSE
         DO WHILE(.NOT.FOUND)
            WRITE(*,*) 
     &        'Enter the name of ADCIRC UNIT 15 (Run Info) file:'
            READ(*,'(A)') RUNINFO
            INQUIRE(FILE=trim(RUNINFO),EXIST=FOUND)
         ENDDO
      ENDIF
      WRITE(*,1011) RUNINFO
      OPEN(15,FILE=RUNINFO)

c...  tcm v50.66.02 Addtions for Time Varying Bathymetry
c...  read through the fort.15 file for the namelist (TimeBathyControl) for
c...  the time varying bathymetry.  This namelist must be at the bottom of the
c...  fort.15 file. If found, then set the appropriate values (btiminc,bchgtiminc,
c...  and nddt).  If the namelist is not there, then the time varying bathymetry
c...  will not be used.
C...
C...  After this search and read, we will close the file and then reopen it
c...  for further processing the traditional non-namelist components.
C...
      READ(UNIT = 15,NML = TimeBathyControl,IOSTAT = IOS_NDDT)

      IF (IOS_NDDT < 0) THEN
c.....   it is possible for the namelist to be present in the file and occuring at the end
c.....   of the file with no line breaks after the ending "\" which causes the iostat to
c.....   return a negative value.  By checking to be sure a namelist variable was set to
c....    a non-default value we can determine this was the case.
         IF (NDDT.NE.0) THEN
!            WRITE(*,*) 'NAMELIST PRESENT, BUT AT THE END OF FILE WITH',
!     &                 ' NO ADVANCING CHARACTER'
            found_tbc_nml = .true.
         ELSE
!            WRITE(*,*) 'NAMELIST NOT PRESENT'
         ENDIF
      ELSEIF (IOS_NDDT == 0) THEN
!         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
         found_tbc_nml = .true.
      ELSE
         found_tbc_nml = .true.
!         WRITE(*,*) 'THERE WAS A PROBLEM PROCESSING THE TimeBathyControl NAMELIST'
!         WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCIRC NOW.'
         CALL EXIT(1)
      ENDIF
C.....REWIND THE FORT.15 FILE IN ORDER TO PROCESS THE REMAINING PIECES
      REWIND(15)

#if defined CSWAN || defined ADCSWAN
      ! jgf50.60.08: Add a namelist for the user to turn SWAN output
      ! on and off. Similar to tcm's timevaryingbathy namelist.
      IOS = 0
      READ(UNIT=15,NML=SWANOutputControl,IOSTAT=IOS)
      IF (IOS.ne.0) THEN
         write(*,*)"INFO: The SWANOutputControl namelist was not found."
         SWAN_OutputHS=.TRUE.
         SWAN_OutputTPS=.TRUE.
         SWAN_OutputDIR=.TRUE.
      ENDIF
      REWIND(15)
#endif
      ! jgf50.60.13: Add a namelist for the user to control met forcing.
      ! Similar to tcm's timevaryingbathy namelist.
      call initWindModule()
      READ(UNIT=15,NML=metControl,IOSTAT=IOS)
      IF (IOS.gt.0) THEN
         write(*,*) "INFO: The metControl namelist was not found."
      ELSE
         found_metCon_nml = .true.  !tcm v50.79 added
      ENDIF
      REWIND(15)

      ! zc: Add namelist for user to control time varying weirs
      IOS = 0
      READ(UNIT=15,NML=TVWControl,IOSTAT=IOS)
      IF(IOS.GT.0)THEN
         IF(NOUT_TVW.NE.0)THEN
             WRITE(*,'(2A)') "INFO: The time varying weir control namelist",
     &          " was found."
             found_tvw_nml = .TRUE.
         ELSE
             found_tvw_nml = .FALSE.
         ENDIF
      ELSEIF(IOS.EQ.0) then
         WRITE(*,'(2A)') "INFO: The time varying weir control namelist",
     &      " was found."
         found_tvw_nml = .TRUE.
      ENDIF
      REWIND(15)

Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      IOS = 0
      READ(UNIT=15,NML=waveCoupling,IOSTAT=IOS)
      IF(IOS.GT.0)THEN
         WRITE(*,*) "INFO: The waveCoupling namelist was not found."
         FOUND_WC_NML = .FALSE.
      ELSE
         FOUND_WC_NML = .TRUE.
      ENDIF
      REWIND(15)

#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      READ(UNIT=15,NML=owiWindNetcdf,IOSTAT=IOS)
      IF(IOS.GT.0)THEN
         WRITE(*,*) "INFO: The owiWindNetcdf namelist was not found."
         FOUND_OWINC_NML = .FALSE.
      ELSE
         FOUND_OWINC_NML = .TRUE.
      ENDIF
      REWIND(15)
#endif

c...  tcm v52.30.01 Addtions for using the Smagorinski Control Namelist
c...  read through the fort.15 file for the namelist (Smag_Control) for
c...  using upper and lower limits on Smagorinski.  This namelist must be at the bottom
c...  of the fort.15 file. If found, then set the appropriate values (Smag_Comp_Flag,
c...  SMAG_UPPER_LIM, and Smag_Lower_Lim).  If the namelist is not there, then the
c...  Smagorinski limit adjustments will not be used.
C...
C...  After this search and read, we will close the file and then reopen it
c...  for further processing the traditional non-namelist components.
C...
      IOS_SMAG = 0
      READ(UNIT = 15,NML = Smag_Control,IOSTAT = IOS_SMAG)

      IF (IOS_SMAG < 0) THEN
c.....   it is possible for the namelist to be present in the file and occuring at the end
c.....   of the file with no line breaks after the ending "\" which causes the iostat to
c.....   return a negative value.  By checking to be sure a namelist variable was set to
c....    a non-default value we can determine this was the case.
         IF (SMAG_COMP_FLAG .eqv. .true. ) THEN
!            WRITE(*,*) 'Smag_Control NAMELIST PRESENT, BUT AT THE',
!     &                 ' END OF THE FILE WITH NO ADVANCING',
!     &                 ' CHARACTER'
            found_smag_nml = .true.
         ELSE
!            WRITE(*,*) 'Smag_Control NAMELIST NOT PRESENT'
         ENDIF
      ELSEIF (IOS_SMAG == 0) THEN
!         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
         found_smag_nml = .true.
      ELSE
         found_smag_nml = .true.
!         WRITE(*,*) 'THERE WAS A PROBLEM PROCESSING THE Smag_Control NAMELIST'
!         WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCPREP NOW.'
         CALL EXIT(1)
      ENDIF
C.....REWIND THE FORT.15 FILE IN ORDER TO PROCESS THE REMAINING PIECES
      REWIND(15)

      READ(UNIT=15,NML=WarnElevControl,IOSTAT=IOS)
      IF (IOS.gt.0) THEN
         write(*,*) "INFO: The WarnElevControl namelist was not found."
      ELSE
         found_WarnElev_nml = .true.  
      ENDIF
      REWIND(15)
      
      READ(UNIT=15,NML=AliDispersionControl,IOSTAT=IOS)
      IF (IOS.gt.0) THEN
         write(*,*) "INFO: The AliDisp Control namelist was not found."
      ELSE
         found_AliDisp_nml = .true.  
      ENDIF
      REWIND(15)
      
      READ(UNIT=15,NML=densityControl,IOSTAT=IOS)
      IF (IOS.gt.0) THEN
         write(*,*) "INFO: The density Control namelist was not found."
      ELSE
         found_density_nml = .true.  
      ENDIF
      REWIND(15)


C
C
C  Initialize all runtime option logicals to false
C
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
C
      L15A: DO WHILE(.TRUE.)
         READ(15,80) LINEI                             !SKIP OVER RUNDES
         DO I=1,95
            IF(CHARI(I).NE.' ') EXIT L15A
         END DO
      ENDDO L15A
      L15B: DO WHILE(.TRUE.)
         READ(15,80) LINEI                             !SKIP OVER RUNID
         DO I=1,95
            IF(CHARI(I).NE.' ') EXIT L15B
         ENDDO
      ENDDO L15B
C
 25   READ(15,*) IDUM  ! SKIP NFOVER
      READ(15,*) IDUM  ! SKIP NABOUT
      READ(15,*) IDUM  ! SKIP NSCREEN
      READ(15,*) IHOT
      IF ((IHOT.EQ.367).OR.(IHOT.EQ.368).OR.
     &    (IHOT.EQ.567).OR.(IHOT.EQ.568)) THEN
         useNetCDF = .true.
      ENDIF
      READ(15,*) IDUM  ! SKIP ICS
C
      READ(15,*) IM                                 !READ IM (model type)
C
      SELECT CASE (IM) ! jgf46.28 added transport/baroclinic support
      CASE(0)
        C2DDI = .TRUE.
      CASE(1)
         C3D  = .TRUE.
         C3DVS  = .TRUE.
      CASE(2)
c        C3D  = .TRUE.
c        C3DDSS = .TRUE.
         print *, "DSS model type not presently supported"
         CALL EXIT(1)
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
         IF ((IM.GE.111111).AND.(IM.LE.534328)) THEN
            C2DDI = .TRUE.
            IF (IM-int(IM*0.1)*10 > 4) THEN
                ! Setting Baroclinic true 
                CBaroclinic = .TRUE.
            ENDIF
         ELSEIF ((IM.GE.611111).AND.(IM.LE.634322)) THEN
           C3D  = .TRUE.
           C3DVS  = .TRUE.
         ELSE
            print *, "model type not supported"
            CALL EXIT(1)
         ENDIF
      END SELECT
C
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
C               and Temperature forcing
         C2D_BTrans = .TRUE.
      CASE(-4)! 2DDI Diagnostic Baroclinic ADCIRC run with Salinity
C               and Temperature forcing'
         ! do nothing
      CASE(-5)! 2DDI Diagnostic Baroclinic ADCIRC run reading from fort.11
C               with pre-computed BPG
      CASE(5) ! 2DDI prognostic Baroclinic ADCIRC run reading from fort.11.nc
C               with pre-computed BPG on structured grid
      CASE(-6)! 2DDI Diagnostic Baroclinic ADCIRC run reading from fort.11
C               with pre-computed momnetum dispersion and BPG
      CASE(6) ! 2DDI prognostic Baroclinic ADCIRC run reading from fort.11.nc
C             ! with pre-computed momentum dispersion and BPG on structured grid
      CASE DEFAULT
         print *, "IDEN=",IDEN," not supported"
         CALL EXIT(1)
      END SELECT
C
      DO IL=1,4                                     !SKIP OVER NOLIBF-NOLICAT
         READ(15,*) IDUM
      ENDDO
      READ(15,*) NWP                                !SKIP NWP
      DO IL=1,NWP
         READ(15,*) skipped                !jgf46.00 skip over nodal attributes
      ENDDO
      READ(15,*) NCOR                               !SKIP OVER NCOR
      READ(15,*) NTIP                               !READ NTIP
      IF (NTIP.NE.0) CTIP = .TRUE.
      READ(15,*) NWS                                !READ NWS
      READ(15,*) NRAMP                              !SKIP OVER NRAMP
      READ(15,*) GRAVITY                            !SKIP OVER GRAVITY
      READ(15,*) TAU0                               !SKIP OVER TAU0
      IF ( (TAU0.LE.-5.d0).AND.(TAU0.GT.-6.d0) ) THEN
         READ(15,*) RDUM, RDUM                      !SKIP OVER TAU0 MIN AND MAX
      ENDIF
      READ(15,*) DT                                 !SKIP OVER DT
      READ(15,*) STATIM                             !SKIP OVER STATIM
      READ(15,*) REFTIM                             !SKIP OVER REFTIM
      MNWLAT = 1
      MNWLON = 1
      MNWP=1
C
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
c. tcm v49.46 unified the special cases above
C.....tcm v49.64.01 Additions for ice
      IF(NWS.EQ.0) THEN
         NWS = 0
         NRS = 0
         NCICE = 0
      ELSE
         NCICE = INT(ABS(NWS)/1000)
         NRS=INT((ABS(NWS) - NCICE*1000)/100)
         NWS=INT((ABS(NWS)- NCICE*1000 - NRS*100))*INT(NWS/ABS(NWS))
      ENDIF

      IF((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
      IF((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
C     jgf46.02 added NWS=8
C     jgfdebug46.02 added NWS=45
C     jgf46.16 merged:
C     cf & cm added NWS=9: asymmetric hurricane wind model
C     rjw added NWS=19: asymmetric hurricane wind model v2.0
C     jie added NWS=20: generalized asymmetric vortex model
C     sb46.28sb01 added NWS=12: OWI format

C     tcm v49.64.01 broke nws = 2 out of the group below for additions to ice
      IF(ABS(NWS).eq.2) then
        IF(NRS.EQ.0) READ(15,*) WTIMINC
        IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
        MNWP=MNP
      ENDIF

C.... tcm v48.4641 additions for ice
      IF((ABS(NWS).EQ.4).OR.(ABS(NWS).eq.5)
     &     .OR.(ABS(NWS).EQ.45).OR.(NWS.EQ.8)
     &     .OR.(ABS(NWS).EQ.15)  !jgf50.38.05: Added NWS=15 for HWind.
     &     .OR.(ABS(NWS).EQ.20)  !jie NWS=20
     &     .OR.(ABS(NWS).EQ.16)  !tcm v51.06.02 added for nws=16 GFDL Met
     &     .OR.(NWS.EQ.10)       !yf added for nws=10 GFS Met
     &     .OR.(ABS(NWS).EQ.12).OR.(ABS(NWS).EQ.13).OR.(ABS(NWS).EQ.14)
     &     .OR.(NWS.EQ.19))THEN
        IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
C++! GML added 5
        IF((NCICE.EQ.0).AND.((NRS.EQ.1).OR.(NRS.EQ.2).OR.(NRS.EQ.4))) 
     &         READ(15,*) WTIMINC,RSTIMINC
        IF((NCICE.EQ.0).AND.((NRS.EQ.3).OR.(NRS.EQ.5))) THEN  !Casey 090825: Fix for NRS = 3.
          IF(ABS(NWS).EQ.8) THEN
              READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
           ELSEIF (ABS(NWS).EQ.20) THEN
              READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
           ELSE
              READ(15,*) WTIMINC,RSTIMINC
           ENDIF
        ENDIF
        IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,CICE_TIMINC
C++! GML added 5
        IF((NCICE.GE.1).AND.((NRS.EQ.1).OR.(NRS.EQ.2).OR.(NRS.EQ.4).OR.(NRS.EQ.5)))
     &                READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC ! NEED TO BE FIXED FOR NWS=20
        IF((NCICE.GE.1).AND.(NRS.EQ.3)) THEN  !Casey 090825: Fix for NRS = 3.
          IF(ABS(NWS).EQ.8) THEN
              READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,
     &                                     RSTIMINC,CICE_TIMINC
           ELSE
              READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
           ENDIF
        ENDIF
        MNWP=MNP
      ENDIF

      !jgf49.0804: Read WTIMINC line if NWS is 29 (VortexOWI)
      IF ((ABS(NWS)).EQ.29) THEN
         SELECT CASE(NRS)
         CASE(0) ! no wave radiation stress
            ! YYYY MM DD HH24 StNum BLAdj WTIMINC pureVortex pureBackgrnd
            READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,WTIMINC,RDUM,RDUM
         CASE DEFAULT ! some kind of wave radiation stress
            ! YYYY MM DD HH24 StNum BLAdj WTIMINC RSTIMINC pureVortex pureBack
            READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,
     &                 WTIMINC,RSTIMINC,RDUM,RDUM
         END SELECT
      ENDIF

      !jgf: Read WTIMINC line if NWS is 30 (GAHM + OWI)
      IF ((ABS(NWS)).EQ.30) THEN
         SELECT CASE(NRS)
         CASE(0) ! no wave radiation stress
                      ! YYYY  MM   DD HH24  StNum BLAdj GEOFACTOR
                      ! WTIMINC pureVortex pureBackgrnd
            READ(15,*) IDUM,IDUM,IDUM,IDUM, IDUM, RDUM,   IDUM,
     &                  WTIMINC,   RDUM,      RDUM
         CASE DEFAULT ! some kind of wave radiation stress
                      ! YYYY MM DD HH24 StNum BLAdj GEOFACTOR WTIMINC RSTIMINC pureVortex pureBack
            READ(15,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,IDUM,
     &                 WTIMINC,RSTIMINC,RDUM,RDUM
         END SELECT
      ENDIF

      IF(NWS.EQ.3) THEN
        READ(15,*) IREFYR                          !SKIP THE REST OF THIS LINE
        READ(15,*) MNWLAT,MNWLON                   !SKIP THE REST OF THIS LINE
        MNWP=MNP
      ENDIF
C.... tcm v49.64.01 additions for ice
C.... xyc v52.30 NWS=7 follow the same format as NWS=6
      IF((ABS(NWS).EQ.6).OR.(ABS(NWS).EQ.7)) THEN
        IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC            !SKIP OVER WTIMINC
        IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC
        IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,CICE_TIMINC            !SKIP OVER WTIMINC
        IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
        MNWP=MNP
      ENDIF

C..... saeed added,
      IF(ABS(NWS).EQ.17)  THEN
        IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC                        !SKIP OVER WTIMINC
        IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC
        IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,CICE_TIMINC            !SKIP OVER WTIMINC
        IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
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
         L15C: DO WHILE(.TRUE.)
            READ(15,80) LINEI
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15C
            ENDDO
         ENDDO L15C
         READ(15,*) RDUM
      ENDDO
      MNTIF=NTIF
      IF(NTIF.EQ.0) MNTIF=1
      READ(15,*) MNBFR                             !PROCESS MNBFR
      DO IL=1,MNBFR                                !SKIP OVER BOUNTAG, & AMIG, FF,etc.
         L15D: DO WHILE(.TRUE.)
            READ(15,80) LINEI
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15D
            ENDDO
         ENDDO L15D
         READ(15,*) RDUM
      ENDDO
      DO IL1=1,MNBFR
         L15E: DO WHILE(.TRUE.)
            READ(15,80) LINEI                          !SKIP OVER ALPHA
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15E
            END DO
         ENDDO L15E
         DO IL2=1,skipNETA
            READ(15,*) RDUM                          !SKIP OVER BOUNDARY FORCINGS
         ENDDO
      ENDDO
      IF(MNBFR.EQ.0) MNBFR=1
      READ(15,*)  RDUM                             !SKIP OVER ANGIN
      MNFFR=0
      IF (fluxBoundary.eqv..true.) READ(15,*) MNFFR   !# FREQ IN NORMAL FLUX B.C.

      DO IL=1,MNFFR                                !SKIP OVER BOUNTAG, & AMIG, FF,etc.
         L15F: DO WHILE(.TRUE.)
            READ(15,80) LINEI
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15F
            ENDDO
         ENDDO L15F
         READ(15,*) RDUM
      ENDDO
      DO IL1=1,MNFFR
         L15G: DO WHILE(.TRUE.)
            READ(15,80) LINEI                          !SKIP OVER ALPHA
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15G
            ENDDO
         ENDDO L15G
         DO IL2=1,NDISC
            READ(15,*) RDUM                          !SKIP OVER BOUNDARY FORCINGS
         ENDDO
      ENDDO
      IF((MNFFR.EQ.0).OR.(MNFFR.EQ.-1)) MNFFR=1
C
C..... ELEVATION STATIONS
C
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
         if (abs(MNSTAE2).ne. abs(MNSTAE)) then
            MNSTAE = abs(MNSTAE2)  !reset the value to what's in the file
         endif
7690     if (ios_stations .ne. 0) then
            write(*,*) "Error in Reading Elevation Station File: elev_stat.151"
            CALL EXIT(1)  ! there is a stop here
         endif
         close(151)
      ELSE
         write(*,*) "Elevation Station Locations contained in fort.15"
         USE_ELEV_STAT_FILE = .FALSE.
         DO IL=1,MNSTAE
            READ(15,*) RDUM                            !SKIP OVER STA COORDS
         ENDDO
      endif
      IF(MNSTAE.EQ.0) MNSTAE=1
C
C..... VELOCITY STATIONS
C
      MNSTAV=0
      MNSTAV2=0
      READ(15,*) IDUM                               !SKIP OVER NOUTV, TOUTSV...
      call setFormatUsage(idum)

      READ(15,*) MNSTAV                             !PROCESS MNSTAV
      ! tcm v51.20.03 -- added negative mnstav to signal use of external station file vel_stat.151
      if (MNSTAV < 0 ) then
         write(*,*) "External File Used for Velocity Station Locations"
         USE_VEL_STAT_FILE = .true.
         MNSTAV = abs(MNSTAV) !reset to a positive value
         ios_stations = 0
         open(unit=151,file='vel_stat.151',status='old',err=7691,iostat=ios_stations)
         read(151,*) MNSTAV2
         if (abs(MNSTAV2).ne. abs(MNSTAV)) then
            MNSTAV = abs(MNSTAV2)  !reset the value to what's in the file
         endif
7691     if (ios_stations .ne. 0) then
            write(*,*) "Error in Reading Velocity Station File: vel_stat.151"
            CALL EXIT(1)  ! there is a stop here
         endif
         close(151)
      else
         write(*,*) "Velocity Station Locations Contained in fort.15"
         USE_VEL_STAT_FILE = .false.
         DO IL=1,MNSTAV
            READ(15,*) RDUM                            !SKIP OVER STA COORDS
         ENDDO
      endif
      IF(MNSTAV.EQ.0) MNSTAV=1
C
C..... CONCENTRATION STATIONS
C
      MNSTAC=0
      MNSTAC2=0
      IF(IM.EQ.10) THEN
         READ(15,*) IDUM                            !SKIP OVER NOUTC, TOUTSC...
         call setFormatUsage(idum)

         READ(15,*) MNSTAC                          !PROCESS MNSTAC
         !tcm v51.20.03 -- added negative mnstac to signal use of external station file conc_stat.151
         if (MNSTAC < 0) then
            write(*,*) "External File Used for Concentration Station Locations"
            USE_CONC_STAT_FILE = .true.
            MNSTAC = abs(MNSTAC) !reset to a positive value
            ios_stations = 0
            open(unit=151,file='conc_stat.151',status='old',err=7692,iostat=ios_stations)
            read(151,*) MNSTAC2
            if (abs(MNSTAC2).ne. abs(MNSTAC)) then
               MNSTAC = abs(MNSTAC2)  !reset the value to what's in the file
            endif
7692        if (ios_stations .ne. 0) then
               write(*,*) "Error in Reading Concentration Station File: conc_stat.151"
               CALL EXIT(1)  ! there is a stop here
            endif
            close(151)
         else
            write(*,*) "Concentration Station Locations Contained in fort.15"
            USE_CONC_STAT_FILE = .false.
            DO IL=1,MNSTAC
               READ(15,*) RDUM                         !SKIP OVER STA COORDS
            ENDDO
         endif
      ENDIF
      IF(MNSTAC.EQ.0) MNSTAC=1
C
C..... MET STATIONS
C
      MNSTAM=0
      MNSTAM2=0
      IF(NWS.NE.0) THEN
        READ(15,*) IDUM                           !SKIP OVER NOUTM, TOUTSM...
        call setFormatUsage(idum)

        READ(15,*) MNSTAM                         !PROCESS MNSTAM
        !tcm v51.20.03 -- added negative mnstam to signal use of external station file met_stat.151
        if (MNSTAM < 0 ) then
           write(*,*) "External File Used for MET Station Locations"
           USE_MET_STAT_FILE = .true.
           MNSTAM = ABS(MNSTAM)  !reset to a positive number
           ios_stations = 0
           open(unit=151,file='met_stat.151',status='old',
     &           err=7693,iostat=ios_stations)
           read(151,*) MNSTAM2
           if (abs(MNSTAM2).ne. abs(MNSTAM)) then
              MNSTAM = abs(MNSTAM2)  !reset the value to what's in the file
           endif
7693       if (ios_stations .ne. 0) then
              write(*,*) "Error in Reading MET Station File: met_stat.151"
              CALL EXIT(1)  ! there is a stop here
           endif
           close(151)
        else
           write(*,*) "MET Station Locations Contained in fort.15"
           USE_MET_STAT_FILE = .false.
           DO IL=1,MNSTAM
              READ(15,*) RDUM                        !SKIP OVER STA COORDS
           END DO
        endif
      ENDIF
      IF(MNSTAM.EQ.0) MNSTAM=1
C
      READ(15,*) IDUM                            !SKIP OVER NOUTGE, TOUTSGE...
      call setFormatUsage(idum)

      READ(15,*) IDUM                            !SKIP OVER NOUTGV, TOUTSGV...
      call setFormatUsage(idum)

      IF (IM.EQ.10) THEN
         READ(15,*) IDUM               !SKIP OVER NOUTGC, TOUTSGC...
         call setFormatUsage(idum)
      ENDIF
      IF (NWS.NE.0) THEN
         READ(15,*) IDUM               !SKIP OVER NOUTGW, TOUTSGW...
         call setFormatUsage(idum)
      ENDIF
      READ(15,*) NHARF                           !PROCESS MNHARF
      DO IL1=1,NHARF
         L15H: DO WHILE(.TRUE.)
            READ(15,80) LINEI                       !SKIP OVER HAFNAM
            DO I=1,95
               IF(CHARI(I).NE.' ') EXIT L15H
            ENDDO
         ENDDO L15H
         READ(15,*) RDUM,RDUM,RDUM             !SKIP OVER HAFREQ,HAFF,HAFACE
      ENDDO
C
      READ(15,*) THAS,THAF,NHAINC,FMV               !READ THAS,...FMV
      READ(15,*) NHASE,NHASV,NHAGE,NHAGV            !READ NHASE,NHASV,...
      CALL setFormatUsage(NHASE)
      CALL setFormatUsage(NHASV)
      CALL setFormatUsage(NHAGE)
      CALL setFormatUsage(NHAGV)
      READ(15,*) NHSTAR,IDUM                          !SKIP OVER NHSTAR,NHSINC
      PRINT *, "NHSTAR = ", NHSTAR
      IF ((ABS(NHSTAR).EQ.3).OR.(ABS(NHSTAR).EQ.367)
     &    .OR.(ABS(NHSTAR).EQ.368).OR.
     &    (ABS(NHSTAR).EQ.5).OR.(ABS(NHSTAR).EQ.567)
     &    .OR.(ABS(NHSTAR).EQ.568) ) THEN
          useNetCDF = .true.                          !for netCDF
      ENDIF
      IF ( (useNetCDF.eqv..true.).AND.(NETCDF_AVAIL.eqv..false.) ) THEN
         WRITE(*,*)'ERROR: NetCDF input or output files were specified.'
         WRITE(*,*) 'but adcprep was not compiled with NetCDF support.'
         WRITE(*,*) 'Please recompile adcprep with NetCDF libraries.'
         CALL EXIT(1)
      ENDIF
      IF ( (useXDMF.eqv..true.).AND.(XDMF_AVAIL.eqv..false.) ) THEN
         WRITE(*,*)'ERROR: XDMF output files were specified.'
         WRITE(*,*) 'but adcprep was not compiled with XDMF support.'
         WRITE(*,*) 'Please recompile adcprep with XDMF libraries.'
         CALL EXIT(1)
      ENDIF
C
C...THIS SECTION TO LUMP THE GWCE MATRIX
Cvjp 11/30/99 made lumping a compile time option

C      PRINT *, "Made it to LUMP definition"
#ifdef LUMP
       CLUMP = .TRUE.
       ILUMP=1
#else
       CLUMP = .FALSE.
       ILUMP=0
#endif

      READ(15,*) ITITER,IDUM,RDUM,IDUM2             !READ SOLVER TYPE
      CSOLIT = .TRUE.

C
C--Read in 3D info
C
      IF(C3DVS) THEN  !3DVS
c     jgf45.10 removed IDIAG
         READ(15,*) IDUM                               !Skip IDEN
         READ(15,*) IDUM,RDUM                          !Skip ISLIP,KP
         READ(15,*) RDUM,RDUM                          !Skip Z0S,Z0B
         READ(15,*) RDUM,RDUM,RDUM                     !Skip ALP1,ALP2,ALP3
C     jgf45.12 add code to handle user specified vertical spacing.
         READ(15,*) IGC,NFEN
         IF (IGC.EQ.0) THEN
            DO I=1,NFEN
               READ(15,*) RDUM  !skip over vertical spacing
            ENDDO
         ENDIF
C     jgf45.12 add code to handle user specified vertical eddy viscosity.
         READ(15,*) IEVC,RDUM,RDUM !Process IEVC
         IF (IEVC.EQ.0) THEN
            DO I=1,NFEN
               READ(15,*) RDUM  !skip over vertical eddy viscosity
            ENDDO
         ENDIF
         IF((IEVC.EQ.50).OR.(IEVC.EQ.51)) READ(15,*) RDUM,RDUM  !Skip THETA1,THETA2
         READ(15,*) I3DSD,RDUM,RDUM,IDUM               !Process I3DSD
         READ(15,*) NSTA3DD                            !Process NSTA3DD
         IF(NSTA3DD.NE.0) THEN !kmd : changed to match 2D options
            DO I=1,NSTA3DD
               READ(15,*) RDUM,RDUM !Skip density stations
            END DO
         ENDIF
         READ(15,*) I3DSV,RDUM,RDUM,IDUM               !Process I3DSV
         READ(15,*) NSTA3DV                            !Process NSTA3DV
         IF (NSTA3DV.NE.0) THEN  !kmd : changed to match 2D options
            DO I=1,NSTA3DV
               READ(15,*) RDUM,RDUM !Skip velocity stations
            END DO
         ENDIF
         READ(15,*) I3DST,RDUM,RDUM,IDUM !Process I3DST
         READ(15,*) NSTA3DT                            !Process NSTA3DT
         IF (NSTA3DT.NE.0) THEN  !kmd : changed to match 2D options
            DO I=1,NSTA3DT
               READ(15,*) RDUM,RDUM !Skip turbulence stations
            ENDDO
         ENDIF
         READ(15,*) IDUM,RDUM,RDUM,IDUM !Skip 3D global density output
         READ(15,*) IDUM,RDUM,RDUM,IDUM                !Skip 3D global velocity output
         READ(15,*) IDUM,RDUM,RDUM,IDUM                !Skip 3D global turbulence output

!kmd48.33bc - added this information for the boundary conditions
         IF (CBAROCLINIC) THEN
           READ(15,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP    !Process RES_BC_FLAG
           IF ((RES_BC_FLAG.LT.0).OR.(RES_BC_FLAG.EQ.1)) THEN
              READ(15,*) RDUM
              READ(15,*) RDUM
           ELSE IF ((RES_BC_FLAG.EQ.2).OR.(RES_BC_FLAG.EQ.3)) THEN
              READ(15,*) RDUM,RDUM
              READ(15,*) RDUM,RDUM
              IF (BCFLAG_TEMP.NE.0) THEN
                 READ(15,*) RDUM, RDUM
              END IF
           ELSE IF (RES_BC_FLAG.EQ.4) THEN
              READ(15,*) RDUM,RDUM,RDUM
              READ(15,*) RDUM,RDUM,RDUM
              IF (BCFLAG_TEMP.NE.0) THEN
                 READ(15,*) RDUM, RDUM
              END IF
           END IF
         END IF

C         PRINT *, "Made it to sponge information"
         IF (CBAROCLINIC) THEN
            READ(15,*) RDUM
         END IF


c     ELSEIF(C3DDSS) THEN  !3DDSS
c
      ENDIF
C
      REWIND(15)                                    !FINISHED WITH UNIT 15 FILE


C...
      WRITE(*,3000) MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,
     &  MNBOU,MNVEL,MNTIF,MNBFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNWP,
     &  MNWLAT,MNWLON,MNFFR
C
      IF(NWS.EQ.0) WRITE(*,3011)
      IF(NWS.EQ.1) WRITE(*,3012)
      IF(ABS(NWS).EQ.2) WRITE(*,3013)
      IF(NWS.EQ.3) WRITE(*,3014)
      IF(ABS(NWS).EQ.4) WRITE(*,3015)
      IF(ABS(NWS).EQ.45) WRITE(*,3015) !jgfdebug46.02
      IF(ABS(NWS).EQ.5) WRITE(*,3016)
      IF(ABS(NWS).EQ.6) WRITE(*,3017)
      IF(ABS(NWS).EQ.7) WRITE(*,3033)  !xyc52.30 added
      IF(ABS(NWS).EQ.8) WRITE(*,3217)
      IF(NRS.GE.1) WRITE(*,3018) ! sb46.28sb03
      IF(NWS.EQ.10) WRITE(*,3019)
      IF(NWS.EQ.11) WRITE(*,3020)
C     sb46.28sb01 added NWS=12: OWI format 
      IF(ABS(NWS).EQ.12) WRITE(*,3023)
      IF(ABS(NWS).EQ.19) WRITE(*,3219)
      IF(ABS(NWS).EQ.20) WRITE(*,3220)
      IF(ABS(NWS).EQ.16) WRITE(*,3026)
      IF(NCICE.GE.1) WRITE(*,3221)  !tcm v49.64.01
      IF(ABS(NWS).EQ.29) WRITE(*,3223)
      IF(ABS(NWS).EQ.30) WRITE(*,3224)
      IF((NHARF.EQ.0).OR.(FMV.EQ.0.)) WRITE(*,3021)
      IF((NHARF.GE.1).AND.(FMV.NE.0.)) WRITE(*,3022)
      IF(ILUMP.EQ.0) WRITE(*,3031)
      IF(ILUMP.EQ.1) WRITE(*,3032)
      IF(IM.EQ.0) WRITE(*,3101)
      IF(IM.EQ.10) WRITE(*,3109)
      IF(IM.EQ.1) WRITE(*,3102)
      IF(IM.EQ.2) WRITE(*,3103)
      IF(ITITER.EQ.0) WRITE(*,3104)
      IF(ITITER.GT.0) WRITE(*,3105)
      IF(ITITER.LT.0) WRITE(*,3106)
      IF(USE_ELEV_STAT_FILE) WRITE(*,3180)  !tcm v51.20.03
      IF(USE_VEL_STAT_FILE) WRITE(*,3181)   !tcm v51.20.03
      IF(USE_MET_STAT_FILE) WRITE(*,3182)   !tcm v51.20.03
      IF(USE_CONC_STAT_FILE) WRITE(*,3183)  !tcm v51.20.03
      WRITE(*,3108)
C
 3000 FORMAT(' *****************************************************',/,
     &       ' *   Based on input and information extracted from   *',/,
     &       ' *   the ADCIRC UNIT 14 and 15 (grid and horiz run   *',/,
     &       ' *   info) files the following paramter values will  *',/,
     &       ' *   be set:                                         *',/,
     &       ' *                                                   *',/,
     &       ' *       MNPROC = ',I5,'                             *',/,
     &       ' *       MNE = ',I8,1X,'     MNP = ',I8,1X,'         *',/,
     &       ' *       MNEI = ',I7,'                               *',/,
     &       ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/,
     &       ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/,
     &       ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/,
     &       ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/,
     &       ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/,
     &       ' *       MNWP = ',I7,'                               *',/,
     &       ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/,
     &       ' *       MNFFR = ',I6,3X,'                           *',/,
     &       ' *                                                   *')
 3011 FORMAT(' *   Also, NO wind forcing will be used,             *')
 3012 FORMAT(' *   Also, NWS=1 meteorological forcing is used,     *')
 3013 FORMAT(' *   Also, NWS=+-2 meteorological forcing is used,   *')
 3014 FORMAT(' *   Also, NWS=3 meteorological forcing is used,     *')
 3015 FORMAT(' *   Also, NWS=+-4 meteorological forcing is used,   *')
 3016 FORMAT(' *   Also, NWS=+-5 meteorological forcing is used,   *')
 3017 FORMAT(' *   Also, NWS=+-6 meteorological forcing is used,   *')
 3033 FORMAT(' *   Also, NWS=+-7 meteorological forcing is used,   *')
 3217 FORMAT(' *   Also, NWS=+-8 Holland wind forcing is used,     *')
 3219 FORMAT(' *   Also, NWS=+-19 Asymmetric Holland wind v2.0     *',/,
     &       ' *                 forcing is used,                  *')
 3220 FORMAT(' *   Also, NWS=+-20 Generalized Asym. Model is used  *')
 3026 FORMAT(' *   Also, NWS=+-16 GFDL Met Data is used            *')
 3221 FORMAT(' *   Also, ABS(NWS)>=1000 ice concentration are used,*')
 3223 FORMAT(' *   Also, NWS=+-29 VortexOWI met forcing is used,   *')
 3224 FORMAT(' *   Also, NWS=+-30 GAHM+OWI met forcing is used,    *')
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
C
  60  FORMAT(A60)
  80  FORMAT(A95)
 180  FORMAT(95A1)
1010  FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
1011  FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
C
C
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
            useNetCDF = .true.
         case(XDMF)
            useXDMF = .true.
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
