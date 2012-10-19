C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.07                                    *
C                                                                             *
C  mod history                                                                *
C  v45.07     - 11/07/05 - jgf- from 45.06 - removed EP from hot start file   *
C  v44.02     - 09/30/03 - rl - from 44.01 - additional changes associated    *
C                                          with removing local variables from *
C                                          global module.                     *
C  v44.01     - 10/06/03 - rl - from 43.03 - changed names of load vectors in *
C                                          2D solution                        *
C  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
C                                          output buffer flush (m.cobb)       *
C                                          3D fixes (k.dresback)              *
C                                          drop MNPROC in fort.15 (t.campbell)*
C                                          various bug fixes in RBCs          *
C                                          ZSURFBUOY/BCPG calc                *
C  v43.00a    - sum  /02 - tc - from 36.01 (3D) & 41.12? (2D), create F90/    *
C                                               parallel unified 2D/3D source *
C  v41.11     - 09/14/01 - rl - from 41.10 - added NWS = -2 capability        *
C  v41.09     - 06/30/01 - jw - from 41.08 - minor mods per vp version 41.05  *
C  v40.02m002 - 12/22 - jjw/vjp - Vic suggested to avoid compiler conflicts   *
C  v40.02m001 - 12/21 - jjw - add cross barrier pipes cjjwm001                *
C                                                                             *
C******************************************************************************
C
      SUBROUTINE HOTSTART()
C
C**************************************************************************
C
C     HOT START PROGRAM SETUP ROUTINE
C
C**************************************************************************
C
      USE GLOBAL
      USE GLOBAL_IO
      USE GLOBAL_3DVS, ONLY: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT,
     &    qsurf, NFEN, DUU, DUV, DVV, UU, VV, BSX, BSY, Q, WZ, q20, l,
     &    SigT, Sal, Temp, N3DSD, I3DSDRec, N3DSV, I3DSVRec, N3DST,
     &    I3DSTRec, N3DGD, I3DGDRec, N3DGV, I3DGVRec, N3DGT, I3DGTRec,
     &    DUU_g, DUV_g, DVV_g, UU_g, VV_g, BSX_g, BSY_g, Q_g, WZ_g,
     &    q20_g, l_g, SigT_g, Sal_g, Temp_g, iy
C
      USE SIZES
#ifdef CMPI
      USE MESSENGER,ONLY:COMM,IERR
#endif      
      USE HARM, ONLY : IHARIND, IHABEG, ITHAS, NHAINC, ICHA, CHARMV,
     &                 NHAGE, NHAGV, NHASE, NHASV, readBinaryHAHotstart,
     &                 checkHAHotstart
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE WIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET   ! sb46.28sb01 added 09/xx/2006
      USE OWI_ICE, ONLY : NCICE1_INIT, NCICE1_GET  !tcm v49.64.01
      USE RS2,ONLY : RS2INIT,RS2GET           ! sb46.28sb01 added 09/xx/2006
C     Addition for netCDF by MCF 5/18/08 modified by jgf49.28
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : readNetCDFHotstart,
     &                     readNetCDFHotstart3D,
     &                     readNetCDFHotstartHarmonic,
     &                     readNetCDFHotstartHarmonicMeansVariances
#endif
      USE NodalAttributes, ONLY :
     &     nolibf, nwp, tau0, cf, eslm,
     &     LoadDirEffRLen, LoadCanopyCoef,
     &     ApplyDirectionalWindReduction, ApplyCanopyCoefficient,
     &     OutputTau0, LoadEleSlopeLim, elemental_slope_limiter_active
      IMPLICIT NONE
#ifdef CMPI       
#ifndef HAVE_MPI_MOD
                  include 'mpif.h'
#endif
#endif
      INTEGER IT
      INTEGER I, J, NH, N              !local loop counters
      INTEGER NCyc, NA, InputFileFmtVn
      INTEGER sd_node_number,fd_node_number

C   kmd48.33bc - add variable for timestep change
      INTEGER NTHS

      REAL(SZ) HollandTime
      REAL(SZ) ArgT, ArgTP, ArgSAlt
      REAL(SZ) H2
      REAL(SZ) CCSFEA
      REAL(SZ) QTRatio
      REAL(SZ) SAltMul, S2SFEA
      REAL(SZ) TPMul
      REAL(SZ) WTRatio, WindX, WindY, WindMag, WDragCo
      REAL(8) TimeLoc
      REAL(SZ) PRDIFF,PRBCKGRND_MH2O     ! tcm v49.16 20100617 added

C     jgf46.08 Fine grained ramp functions
      REAL(SZ) RampExtFlux1   ! Ramp for external flux b.c.s @ITHS-1
      REAL(SZ) RampExtFlux2   ! Ramp for external flux b.c.s @ITHS
      REAL(SZ) RampIntFlux1   ! Ramp for internal flux b.c.s @ITHS-1
      REAL(SZ) RampIntFlux2   ! Ramp for internal flux b.c.s @ITHS
      REAL(SZ) RampTip2    ! Ramp for tidal potential @ITHS
      REAL(SZ) RampMete2   ! Ramp for wind and atmospheric pressure @ITHS
      REAL(SZ) RampWRad2   ! Ramp for wave radiation stress @ITHS
C     jgf46.14 Check for file existence before attempting to open.
      LOGICAL fileFound    ! .True. if the the file exists
      INTEGER ErrorIO      ! zero if file opened successfully
      INTEGER NP_G_IN, NE_G_IN  ! Global and
      INTEGER NP_A_IN, NE_A_IN  ! active number nodes, elements.
      CHARACTER(len=7) :: hsfn  ! name of binary hotstart file
C
C     jgf48.4628 Add capability to turn off solution if only met output
C     was requested.
      LOGICAL metOutputOn     ! .true. if met output was requested
      LOGICAL nonMetOutputOff ! .true. no other output was requested
C
C     jgf48.4636 Add capability to read min/max files upon hotstart,
C     if they are available
      CHARACTER(len=12) minMaxFN ! name of the min or max file to read in
      INTEGER node               ! node number of the min/max data
      CHARACTER(len=80) skipline ! dummy variable for min/max header data
      INTEGER numLinesInMinMaxFile ! count the lines to report to log file
      LOGICAL tooFewMinMaxLines  ! true if couldn't read enough lines from file
C
C     kmd48.33bc - add variables for heat flux
      INTEGER :: NOD
      REAL(SZ),ALLOCATABLE :: TMP(:,:)
      REAL(SZ) :: CD, CDQ, QWIND
      CHARACTER(80) :: CDUM80
      REAL(SZ) :: RNDAYHS
      REAL(SZ),ALLOCATABLE :: RealQ_Tmp(:,:)
      REAL(SZ),ALLOCATABLE :: ImQ_Tmp(:,:)
      INTEGER :: IDenHS
      REAL(8) :: dtChange ! difference in time step size (sec) between 
                          ! current run and the one we're hotstarting from
C
C     jgf49.44: Add variables for reading fulldomain binary hotstart file
C     in parallel.
      INTEGER subdomainNode  ! subdomain node number to map fulldomain data to
      INTEGER fulldomainNode ! fulldomain node number from hotstart file
      INTEGER subdomainElement ! subdomain element num to map fulldoman data to
      INTEGER fulldomainElement ! fulldomain element number from hotstart file
C
C     tcm v49.64.01 -- ADDED variables for Ice concentration fields
      REAL(SZ) PIC,CICE_TRatio  !ICE VARIABLES
C
C     tcm v50.66.01 added variables for time varying bathymetry
      real(sz) DPTMP,BTRatio
      integer IJK
C
C     jgf49.1001 Added for blending NWS19 and NWS12 (from NAM) met data
      REAL(SZ) bf ! blending factor for the two wind data types

C     kmd -- Added variables for initial conditions
      CHARACTER(80) :: gridinfo

      INTEGER :: moderbin ! determines whether in binary or real
      INTEGER :: ncHundreds ! hundreds place of IHOT, used for netcdf hs file

      call setMessageSource("hotstart")
#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      metOutputOn = .FALSE.
      nonMetOutputOff = .FALSE.
C
C     -----------------------------------
C
C     R E A D   H O T S T A R T   F I L E
C
C     -----------------------------------
C     The value of IHOT indicates the type of
C     the file as well as the file name.
      SELECT CASE(IHOT)

      CASE(17) ! initial conditions or ascii hot start

C      kmd49.49 - redo case(67,68) for initial conditions or ascii
C                 format for hot start.
C...
C...     READ IN 2DDI HOT START INITIAL CONDITION OVER WHOLE DOMAIN
C...
C        Determine if the hot start file exists.
         fileFound=.False.
         IF (IHOT.eq.17) THEN
            hsfn = "fort.17"
         ENDIF
         IF ((MNPROC.eq.1).or.
     &       (READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
            INQUIRE(FILE=TRIM(INPUTDIR)//'/'//hsfn,
     &         EXIST=fileFound)
         ELSE
            ! look for fulldomain hotstart file
            INQUIRE(FILE=TRIM(GBLINPUTDIR)//'/'//hsfn,
     &         EXIST=fileFound)
         ENDIF
C
         IF (fileFound.eqv..false.) THEN
            WRITE(16,1701)         ! hot start file
            WRITE(16,1711)         ! was not found.
            WRITE(16,9773)         ! execution terminated
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1701)
               WRITE(ScreenUnit,1711)
               WRITE(16,9773)      ! execution terminated
            ENDIF
            CALL ADCIRC_Terminate()
         ENDIF
C        start reading ascii initial conditions file
         WRITE(16,1740)
 1740    FORMAT(/,9X,
     &    'Initial condition or ascii hot start file was found.',
     &     ' Opening file.')

         IF ((MNPROC.eq.1).or.
     &       (READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
            OPEN(IHOT,FILE=TRIM(INPUTDIR)//'/'//hsfn,ACTION='READ',
     &           ACCESS='SEQUENTIAL',RECL=8,IOSTAT=ErrorIO,STATUS='OLD')
         ELSE
            ! fulldomain hotstart file
            OPEN(IHOT,FILE=TRIM(GBLINPUTDIR)//'/'//hsfn,ACTION='READ',
     &           ACCESS='SEQUENTIAL',RECL=8,IOSTAT=ErrorIO,STATUS='OLD')
         ENDIF

         IF (ErrorIO.GT.0) THEN
            WRITE(16,1701)         ! hot start file
            WRITE(16,1705)         ! exists but can't be opened
            WRITE(16,9773)         ! execution terminated
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1701)
               WRITE(ScreenUnit,1705)
               WRITE(ScreenUnit,9773)
            ENDIF
#ifdef CMPI
            CALL MSG_FINI()
#endif
            STOP                   ! We're toast.
         ENDIF
C
 1701    FORMAT('ERROR: The initial condition or ascii hot start file')
 1711    FORMAT('was not found.')
 1712    FORMAT('was a nonmatching version')
 1705    FORMAT('exists but cannot be opened.')
 9773    FORMAT(/,1X,'!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!',//)
C
C        Now read the readable initial condition or hot start file.
         moderbin=2
         IHOTSTP=1
         READ(IHOT,*) gridinfo
         READ(IHOT,*) IMHS
         READ(IHOT,*) TimeLoc
         READ(IHOT,*) ITHS

         READ(IHOT,*) NP_G_IN
         READ(IHOT,*) NE_G_IN
         READ(IHOT,*) NP_A_IN
         READ(IHOT,*) NE_A_IN
C
C        jgf49.44: Add capability to read fulldomain binary hotstart
C        file in parallel.
         CALL readAndMapToSubdomain2D(ETA1, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(ETA2, IHOT, IHOTSTP, moderbin)
C        jgf46.34 Added support for IBTYPE=52.
         CALL readAndMapToSubdomain2D(ETADisc, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(UU2, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(VV2, IHOT, IHOTSTP, moderbin)
         IF(IMHS.EQ.10) THEN
            CALL readAndMapToSubdomain2D(CH1, IHOT, IHOTSTP, moderbin)
         ENDIF
         CALL readAndMapToSubdomainInt2D(NNODECODE, IHOT, IHOTSTP, moderbin)
         CALL readAndMapEleToSubdomainInt2D(NOFF, IHOT, IHOTSTP, moderbin)
C
         READ(IHOT,*) IESTP
         READ(IHOT,*) NSCOUE
         READ(IHOT,*) IVSTP
         READ(IHOT,*) NSCOUV
         READ(IHOT,*) ICSTP
         READ(IHOT,*) NSCOUC
         READ(IHOT,*) IPSTP
         READ(IHOT,*) IWSTP
         READ(IHOT,*) NSCOUM
         READ(IHOT,*) IGEP
         READ(IHOT,*) NSCOUGE
         READ(IHOT,*) IGVP
         READ(IHOT,*) NSCOUGV
         READ(IHOT,*) IGCP
         READ(IHOT,*) NSCOUGC
         READ(IHOT,*) IGPP
         READ(IHOT,*) IGWP
         READ(IHOT,*) NSCOUGW
         !
         !     R E A D   3 D   H O T S T A R T   D A T A
         !
         IF (C3D) THEN
            READ(IHOT,*) IDenHS
            !3D station density output time counter
            READ(IHOT,*) N3DSD
            READ(IHOT,*) I3DSDRec
            !3D station velocity output time counter
            READ(IHOT,*) N3DSV
            READ(IHOT,*) I3DSVRec
            !3D station turbulence output time counter
            READ(IHOT,*) N3DST
            READ(IHOT,*) I3DSTRec
            READ(IHOT,*) N3DGD
            READ(IHOT,*) I3DGDRec
            READ(IHOT,*) N3DGV
            READ(IHOT,*) I3DGVRec
            READ(IHOT,*) N3DGT
            READ(IHOT,*) I3DGTRec
C           kmd48.33bc - changed to the format used with the global IO
C           jgf49.44: Add capability to read fulldomain 3D binary hotstart
C           file in parallel.
            CALL readAndMapToSubdomain2D(DUU, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(DUV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(DVV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(UU, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(VV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(BSX, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(BSY, IHOT, IHOTSTP, moderbin)
C           kmd48.33bc - allocate realQ_tmp
            ALLOCATE(RealQ_Tmp(NP,NFEN))
            ALLOCATE(ImQ_Tmp(NP,NFEN))
C           kmd48.33bc - changed to format used in global IO, layer output
            CALL readAndMapToSubdomain3D(RealQ_Tmp, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(ImQ_Tmp, NFEN, IHOT, IHOTSTP, moderbin)
            Q = RealQ_Tmp + iy*ImQ_Tmp ! reconstruct Q matrix from parts
            DEALLOCATE(RealQ_Tmp, ImQ_Tmp)
            CALL readAndMapToSubdomain3D(WZ, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(q20, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(l, NFEN, IHOT, IHOTSTP, moderbin)
C           kmd48.33bc - need to determine if IDEN parameter in hotstart
C           is different than IDEN from fort.15 file. If so, determine which
C           one should be used in the simulation.
            IF (IDenHS.NE.IDEN) THEN
               IF ((CBAROCLINIC).AND.(RES_BC_FLAG.GE.1)) THEN
                  ! IDEN value from fort.15 should be used
               ELSE IF ((CBAROCLINIC).AND.(RES_BC_FLAG.LE.0)) THEN
                  IDEN=IDenHS
               END IF
            END IF
            SELECT CASE(ABS(IDen))
            CASE(0)
               CALL logMessage(INFO,
     &          "IDEN=0; This is a barotropic model run.")
            CASE(1)
               CALL readAndMapToSubdomain3D(SigT, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(2)
               CALL readAndMapToSubdomain3D(Sal, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(3)
               CALL readAndMapToSubdomain3D(Temp, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(4)
               CALL readAndMapToSubdomain3D(Sal, NFEN, IHOT, IHOTSTP, moderbin)
               CALL readAndMapToSubdomain3D(Temp, NFEN, IHOT, IHOTSTP, moderbin)
            CASE DEFAULT
               ! should not happen
               WRITE(scratchMessage,1722) IDEN
               CALL allMessage(ERROR,scratchMessage)
1722             FORMAT("IDEN is ",I2,
     &             " but the valid range is +/- 0 to 4.")
            END SELECT
         ENDIF  ! C3D
C...
C......  HOT START INFORMATION FOR HARMONIC ANALYSIS
C...
!         IF(IHARIND.EQ.1) THEN
!            IHABEG=ITHAS+NHAINC
C...
C........   IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........   HARMONIC ANALYSIS, MEAN AND SQUARE INFO
C...
!            IF(ITHS.GT.ITHAS) THEN
!               READ(IHOT,REC=IHOTSTP) ICHA
!            ENDIF
!            IF(ITHS.GE.IHABEG) THEN
!               CALL readBinaryHAHotstart(IHOT, IHOTSTP)
!            ENDIF
!         ENDIF    !IHARIND
         CLOSE(IHOT)  ! close portable ascii hotstart file
C
C

      CASE(67,68) ! non-portable binary file format

C...
C...     READ IN 2DDI HOT START INITIAL CONDITION OVER WHOLE DOMAIN
C...     THIS FILE ALWAYS HAS A RECL=8 BECAUSE IT IS ASSUMED THAT THE HARMONIC
C...     ANALYSIS IS ALWAYS DONE IN 64 BITS, EVEN ON A WORKSTATION
C...
C        Determine if the hot start file exists.
         fileFound=.False.
         IF (IHOT.eq.67) THEN
            hsfn = "fort.67"
         ELSE
            hsfn = "fort.68"
         ENDIF
         IF ((MNPROC.eq.1).or.
     &       (READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
            INQUIRE(FILE=TRIM(INPUTDIR)//'/'//hsfn,
     &         EXIST=fileFound)
         ELSE
            ! look for fulldomain hotstart file
            INQUIRE(FILE=TRIM(GBLINPUTDIR)//'/'//hsfn,
     &         EXIST=fileFound)
         ENDIF
C
         IF (fileFound.eqv..false.) THEN
            WRITE(16,1001)         ! hot start file
            WRITE(16,1011)         ! was not found.
            WRITE(16,9973)         ! execution terminated
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)
               WRITE(ScreenUnit,1011)
               WRITE(16,9973)      ! execution terminated
            ENDIF
#ifdef CMPI
            CALL MSG_FINI()
#endif
            STOP
         ENDIF
C        start reading binary hot start file
         WRITE(16,240)
 240     FORMAT(/,9X,'Hot start file was found.',
     &     ' Opening file.')

         IF ((MNPROC.eq.1).or.
     &       (READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
            OPEN(IHOT,FILE=TRIM(INPUTDIR)//'/'//hsfn,ACTION='READ',
     &            ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO,STATUS='OLD')
         ELSE
            ! fulldomain hotstart file
            OPEN(IHOT,FILE=TRIM(GBLINPUTDIR)//'/'//hsfn,ACTION='READ',
     &            ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO,STATUS='OLD')
         ENDIF

         IF (ErrorIO.GT.0) THEN
            WRITE(16,1001)         ! hot start file
            WRITE(16,1005)         ! exists but can't be opened
            WRITE(16,9973)         ! execution terminated
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)
               WRITE(ScreenUnit,1005)
               WRITE(ScreenUnit,9973)
            ENDIF
#ifdef CMPI
            CALL MSG_FINI()
#endif
            STOP                   ! We're toast.
         ENDIF
C
 1001    FORMAT('ERROR: The hot start file')
 1011    FORMAT('was not found.')
 1012    FORMAT('was a nonmatching version')
 1005    FORMAT('exists but cannot be opened.')
 9973    FORMAT(/,1X,'!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!',//)
C
C        Now read the non-portable binary hot start file.
         moderbin=1 ! kmd - sets the flag for the readAndMapToSubdomain subroutines
         IHOTSTP=1
         READ(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
         if (.not.
     &     CMP_VERSION_NUMBERS(InputFileFmtVn,FileFmtVersion))then
            WRITE(16, 1001)
            write(16, 1012)
            write(16, 9973)
            if (NScreen /= 0 .and. myproc == 0) then
               WRITE(ScreenUnit, 1001)
               write(ScreenUnit, 1012)
               write(ScreenUnit, 9973)
            endif
#ifdef CMPI
            call msg_fini()
#endif
            stop
         endif
         READ(IHOT,REC=IHOTSTP) IMHS         ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) TimeLoc         ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) ITHS         ; IHOTSTP = IHOTSTP + 1

         READ(IHOT,REC=IHOTSTP) NP_G_IN      ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NE_G_IN      ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NP_A_IN      ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NE_A_IN      ; IHOTSTP = IHOTSTP + 1
C
C        jgf49.44: Add capability to read fulldomain binary hotstart
C        file in parallel.
         CALL readAndMapToSubdomain2D(ETA1, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(ETA2, IHOT, IHOTSTP, moderbin)
C        jgf46.34 Added support for IBTYPE=52.
         CALL readAndMapToSubdomain2D(ETADisc, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(UU2, IHOT, IHOTSTP, moderbin)
         CALL readAndMapToSubdomain2D(VV2, IHOT, IHOTSTP, moderbin)
         IF(IMHS.EQ.10) THEN
            CALL readAndMapToSubdomain2D(CH1, IHOT, IHOTSTP, moderbin)
         ENDIF
         CALL readAndMapToSubdomainInt2D(NNODECODE, IHOT, IHOTSTP, moderbin)
         CALL readAndMapEleToSubdomainInt2D(NOFF, IHOT, IHOTSTP, moderbin)
C
         READ(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
         READ(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
         !
         !     R E A D   3 D   H O T S T A R T   D A T A
         !
         IF (C3D) THEN
            READ(IHOT,REC=IHOTSTP) IDenHS
            IHOTSTP=IHOTSTP+1
            !3D station density output time counter
            READ(IHOT,REC=IHOTSTP) N3DSD
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DSDRec
            IHOTSTP=IHOTSTP+1
            !3D station velocity output time counter
            READ(IHOT,REC=IHOTSTP) N3DSV
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DSVRec
            IHOTSTP=IHOTSTP+1
            !3D station turbulence output time counter
            READ(IHOT,REC=IHOTSTP) N3DST
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DSTRec
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) N3DGD
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DGDRec
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) N3DGV
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DGVRec
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) N3DGT
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) I3DGTRec
            IHOTSTP=IHOTSTP+1
C           kmd48.33bc - changed to the format used with the global IO
C           jgf49.44: Add capability to read fulldomain 3D binary hotstart
C           file in parallel.
            CALL readAndMapToSubdomain2D(DUU, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(DUV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(DVV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(UU, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(VV, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(BSX, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain2D(BSY, IHOT, IHOTSTP, moderbin)
C           kmd48.33bc - allocate realQ_tmp
            ALLOCATE(RealQ_Tmp(NP,NFEN))
            ALLOCATE(ImQ_Tmp(NP,NFEN))
C           kmd48.33bc - changed to format used in global IO, layer output
            CALL readAndMapToSubdomain3D(RealQ_Tmp, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(ImQ_Tmp, NFEN, IHOT, IHOTSTP, moderbin)
            Q = RealQ_Tmp + iy*ImQ_Tmp ! reconstruct Q matrix from parts
            DEALLOCATE(RealQ_Tmp, ImQ_Tmp)
            CALL readAndMapToSubdomain3D(WZ, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(q20, NFEN, IHOT, IHOTSTP, moderbin)
            CALL readAndMapToSubdomain3D(l, NFEN, IHOT, IHOTSTP, moderbin)
C           kmd48.33bc - need to determine if IDEN parameter in hotstart
C           is different than IDEN from fort.15 file. If so, determine which
C           one should be used in the simulation.
            IF (MYPROC.EQ.0) THEN
               PRINT *, "IDEN = ", IDEN
               PRINT *, "IDENHS = ", IDENHS
            END IF
            IF (IDenHS.NE.IDEN) THEN
               IF ((CBAROCLINIC).AND.(RES_BC_FLAG.GE.1)) THEN
                  ! IDEN value from fort.15 should be used
               ELSE IF ((CBAROCLINIC).AND.(RES_BC_FLAG.LE.0)) THEN
                  IDEN=IDenHS
               END IF
            END IF
            SELECT CASE(ABS(IDen))
            CASE(0)
               CALL logMessage(INFO,
     &          "IDEN=0; This is a barotropic model run.")
            CASE(1)
               CALL readAndMapToSubdomain3D(SigT, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(2)
               CALL readAndMapToSubdomain3D(Sal, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(3)
               CALL readAndMapToSubdomain3D(Temp, NFEN, IHOT, IHOTSTP, moderbin)
            CASE(4)
               CALL readAndMapToSubdomain3D(Sal, NFEN, IHOT, IHOTSTP, moderbin)
               CALL readAndMapToSubdomain3D(Temp, NFEN, IHOT, IHOTSTP, moderbin)
            CASE DEFAULT
               ! should not happen
               WRITE(scratchMessage,12) IDEN
               CALL allMessage(ERROR,scratchMessage)
12             FORMAT("IDEN is ",I2,
     &             " but the valid range is +/- 0 to 4.")
            END SELECT
         ENDIF  ! C3D
C...
C......  HOT START INFORMATION FOR HARMONIC ANALYSIS
C...
         IF(IHARIND.EQ.1) THEN
            IHABEG=ITHAS+NHAINC
C...
C........   IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........   HARMONIC ANALYSIS, MEAN AND SQUARE INFO
C...
            IF(ITHS.GT.ITHAS) THEN
               READ(IHOT,REC=IHOTSTP) ICHA
            ENDIF
            IF(ITHS.GE.IHABEG) THEN
               CALL readBinaryHAHotstart(IHOT, IHOTSTP)
            ENDIF
         ENDIF    !IHARIND
         CLOSE(IHOT)  ! close non-portable binary hotstart file
C
C
#ifdef ADCNETCDF
      CASE(367,368,567,568)  ! netcdf hotstart file
         ncHundreds = 300
         IF ((IHOT.eq.567).OR.(IHOT.eq.568)) THEN
            ncHundreds = 500
         ENDIF
         CALL readNetCDFHotstart(IHOT-ncHundreds, timeLoc)
         IF (C3D) THEN
            CALL readNetCDFHotstart3D(IHOT-ncHundreds)
         ENDIF
         IF ((IHARIND.eq.1).AND.(ITHS.gt.ITHAS)) THEN
            CALL readNetCDFHotstartHarmonic(IHOT-ncHundreds)
            IF (CHARMV.eqv..true.) THEN
               CALL readNetCDFHotstartHarmonicMeansVariances(IHOT
     &             -ncHundreds)
            ENDIF
         ENDIF
#endif
C
C
      CASE DEFAULT
         ! This code should be unreachable since IHOT was already
         ! checked to be sure that it is valid.
         write(scratchMessage,'(A,I3,A)') "IHOT was set to '",IHOT,
     &    "' but this value of IHOT is not valid."
         call allMessage(ERROR, scratchMessage)
         call ADCIRC_Terminate()
      END SELECT
C
C     Compare the harmonic analysis data read in from the hotstart file
C     with data from the fort.15 file to be sure that they are consistent.
      IF((IHARIND.EQ.1).AND.(ITHS.GE.IHABEG)) THEN
         CALL checkHAHotstart()
      ENDIF

C     ---------------------------------------------
C
C     I N I T I A L I Z E   H O T S T A R T   R U N
C
C     ---------------------------------------------
      !jgf51.01: Determine if there is a difference between the time 
      ! step size in the hotstart file, and the time step size for this
      ! run (i.e, in the fort.15 file) ... if so, this affects any 
      ! calculation of time in seconds since cold start, and any 
      ! feature that counts time steps since cold start
      DTDPHS=TimeLoc/ITHS  ! Set up time step information for hotstart file.
      dtChange = DTDPHS - DTDP
      ! we don't use .ne. for dtChange b/c of the limitations in machine precision
      IF ( (abs(dtChange).gt.1.0e-10).OR.
     &     ( (IDEN.GE.1).AND.
     &       (C3D.eqv..true.).AND.
     &       (CBAROCLINIC.eqv..true.) ) ) THEN
         CHOTHS=.TRUE.
         IF (DTDPHS.NE.DTDP) THEN
            ! must recalculate the ending time step number
            RNDAYHS=TimeLoc/86400.d0
            NTHS=INT((RNDAY-RNDAYHS)*(86400.d0/DTDP)+0.5d0)
            NT=ITHS+NTHS
            call logMessage(INFO,"The time step size is different "//
     &        "from the time step size listed in the hotstart file. "//
     &        "The total number of time steps in this "//
     &        "run has been recalculated accordingly.")
            WRITE(scratchMessage,
     &          '("The ending time step number is now",I9,".")') NT
            call logMessage(INFO,scratchMessage)
            DO I=1,NP
               ETA1(I)=ETA2(I) + DTDP*((ETA1(I)-ETA2(I))/DTDPHS)
            END DO
         END IF
      END IF

C     kmd48.33bc - added time and time step to be passed to the hotstart file
      IF(C3D) THEN
        CALL HOTSTART_3D(TimeLoc,ITHS)
      ENDIF
C
      DO I=1, NP
         ETAS(I)=ETA2(I)-ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         QX2(I)=UU2(I)*H2
         QY2(I)=VV2(I)*H2
      END DO
C     jgf46.08 Fine grained ramp functions
C     jgf46.21 Split flux into internal and external, added support
C     for IBTYPE=52.
!Kendra: Changed these to DTDPHS because if timestep is changed then
!       the value of DTDP will not match with the ITHS.
!
      IF (NRamp.eq.0) THEN
         RampExtFlux1=1.0d0
         RampExtFlux2=1.0d0
         RampIntFlux1=1.0d0
         RampIntFlux2=1.0d0
         RampTip2=1.0d0
         RampMete2=1.0d0
         RampWRad2=1.0d0
      ELSE
         RampExtFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampExtFlux)
         RampExtFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampExtFlux)
         RampIntFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampIntFlux)
         RampIntFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampIntFlux)
         RampTip2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampTip)
Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
         RampMete2=TANH((2.D0*((ITHS)*DTDPHS/86400.D0-DUnRampMete))/DRampMete)
         RampWRad2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampWRad)
      ENDIF
#ifdef IBM
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
#else
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
#endif
      IF(ITHS.LT.(FluxSettlingIT+10)) THEN
         RampIntFlux1=0.0d0
         RampIntFlux2=0.0d0
         RampTip2=0.0d0
         RampMete2=0.0d0
         RampWRad2=0.0d0
      ENDIF
C
C     jgf48.4627 If user only requests meteorological output, then
C     set a flag that turns off the calculation of the GWCE and
C     momentum equations
      IF (NWS.NE.0) THEN
         IF ((NOUTE.EQ.0) .AND. (NOUTV.EQ.0) .AND. (NOUTC.EQ.0) .AND.
     &       (NOUTGE.EQ.0) .AND. (NOUTGV.EQ.0) .AND. (NOUTGC.EQ.0 ).AND.
     &       (I3DSD.EQ.0) .AND. (I3DSV.EQ.0) .AND. (I3DST.EQ.0) .AND.
     &       (I3DGD.EQ.0) .AND. (I3DGV.EQ.0) .AND. (I3DGT.EQ.0) .AND.
     &       (NHASE.EQ.0) .AND. (NHASV.EQ.0) .AND.
     &       (NHAGE.EQ.0) .AND. (NHAGV.EQ.0) )
     &       THEN
            nonMetOutputOff = .TRUE.
         ENDIF
         IF ( (NOUTM.NE.0) .OR. (NOUTGW.NE.0) ) THEN
            metOutputOn = .TRUE.
         ENDIF
         IF ( (metOutputOn.eqv..true.) .AND. 
     &        (nonMetOutputOff.eqv..true.) ) THEN
            call allMessage(INFO,'Only meterological output was '//
     &       'requested. ADCIRC will not solve the GWCE or '//
     &       'momentum equations.')
             METONLY = .TRUE. ! set flag
         ENDIF
      ENDIF
C
C.... SET POSITIONS IN BOUNDARY CONDITION, WIND AND OUTPUT FILES
C
      WRITE(16,1112)
      WRITE(16,1794)
 1794 FORMAT(//,' INFORMATION ABOUT RE-STARTING THE TIME SERIES',
     &     ' OUTPUT FILES (UNITS 61-64,71-74,81,83),',
     &     /,' WIND/PRESSURE FILE (UNIT 22) AND FLOW BOUNDARY',
     &     ' CONDITION FILE (UNIT 20)',//)

C......INITIALLY, ZERO OUT THE NORMAL FLOW ON ALL BOUNDARIES

      DO I=1,NVEL
         QN2(I)=0.D0
         QN1(I)=0.D0
         QN0(I)=0.D0
      END DO

C.... FIND THE PROPER PLACE IN THE APERIODIC ELEVATION SPECIFIED
C.... BOUNDARY CONDITION FILE IF IT IS REQURIED.

      IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1976)
         WRITE(16,1976)
 1976    FORMAT(/,1X,'LOCATING ELEVATION SPECIFIED INFORMATION IN ',
     &        'UNIT 19',/)
         OPEN(19,FILE=TRIM(INPUTDIR)//'/'//'fort.19')
         READ(19,*) ETIMINC
         ETIME1=STATIM*86400.D0
         ETIME2=ETIME1+ETIMINC
         DO J=1,NETA
            READ(19,*) ESBIN1(J)
         END DO
         DO J=1,NETA
            READ(19,*) ESBIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
            IF(TIMEIT.GT.ETIME2) THEN
               ETIME1=ETIME2
               ETIME2=ETIME1+ETIMINC
               DO J=1,NETA
                  ESBIN1(J)=ESBIN2(J)
                  READ(19,*) ESBIN2(J)
               END DO
            ENDIF
         END DO
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TIMEIT-ETIME1)/ETIMINC
      ENDIF

C......FIND PROPER PLACE IN THE APERIODIC NORMAL FLOW BOUNDARY CONDITION
C......FILE IF IT IS REQUIRED

      IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1978)
         WRITE(16,1978)
 1978    FORMAT(/,1X,'LOCATING NORMAL FLOW INFORMATION IN UNIT 20',/)
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC
         QTIME1=STATIM*86400.D0
         QTIME2=QTIME1+FTIMINC
         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN1(J)
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
            IF(TIMEIT.GT.QTIME2) THEN
               QTIME1=QTIME2
               QTIME2=QTIME2+FTIMINC
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     QNIN1(J)=QNIN2(J)
                     READ(20,*) QNIN2(J)
                  ENDIF
               END DO
            ENDIF
         END DO
         QTRATIO=(TIMEIT-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN1(I)=RampExtFlux1*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
         IF(TimeLoc.GT.QTIME2) THEN
            QTIME1=QTIME2
            QTIME2=QTIME1+FTIMINC
            DO J=1,NVEL
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &              .OR.(LBCODEI(J).EQ.22)) THEN
                  QNIN1(J)=QNIN2(J)
                  READ(20,*) QNIN2(J)
               ENDIF
            END DO
         ENDIF
         QTRATIO=(TimeLoc-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN2(I)=RampExtFlux2*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
      ENDIF

! kmd - added in a option to start the river information at the hot start time
      IF((NFLUXF.EQ.1).AND.(NFFR.EQ.-1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1978)
         WRITE(16,1978)
! 1978    FORMAT(/,1X,'LOCATING NORMAL FLOW INFORMATION IN UNIT 20',/)
C         IF (MYPROC.EQ.7) THEN
C            PRINT *, "Made it to the reading in of the river file"
C         END IF
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC
         QTIME1=TIMELOC
         QTIME2=QTIME1+FTIMINC
         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN1(J)
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN2(J)
         END DO
         QTRATIO=(TIMEIT-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN1(I)=RampExtFlux1*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
         IF(TIMELOC.GT.QTIME2) THEN
            QTIME1=QTIME2
            QTIME2=QTIME1+FTIMINC
            DO J=1,NVEL
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &              .OR.(LBCODEI(J).EQ.22)) THEN
                  QNIN1(J)=QNIN2(J)
                  READ(20,*) QNIN2(J)
               ENDIF
            END DO
         ENDIF
         QTRATIO=(TIMELOC-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN2(I)=RampExtFlux2*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
      ENDIF

C......RESTART THE PERIODIC NORMAL FLOW BOUNDARY CONDITION

      IF((NFLUXF.EQ.1).AND.(NFFR.GT.0)) THEN
         DO J=1,NFFR
            IF(FPER(J).EQ.0.) THEN
               NCYC=0.
            ELSE
#ifdef IBM
               NCYC=INT(TimeLoc/FPER(J),KIND(0.0d0))
#else
               NCYC=INT(TimeLoc/FPER(J))
#endif
            ENDIF
            ARGJ1=FAMIG(J)*(TimeLoc-DTDPHS-NCYC*FPER(J))+FFACE(J) ! kmd48.33bc - changed from DTDP to DTDPHS
            ARGJ2=FAMIG(J)*(TimeLoc-NCYC*FPER(J))+FFACE(J)
            RFF1=FFF(J)*RampExtFlux1
            RFF2=FFF(J)*RampExtFlux2
            DO I=1,NVEL
               ARG1=ARGJ1-QNPH(J,I)
               ARG2=ARGJ2-QNPH(J,I)
               QN1(I)=QN1(I)+QNAM(J,I)*RFF1*COS(ARG1)
               QN2(I)=QN2(I)+QNAM(J,I)*RFF2*COS(ARG2)
            END DO
         END DO
      ENDIF

C...
C...  RESTART SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C.... EXTERNAL BARRIER BOUNDARY NODES
C...
      IF(NFLUXB.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.3).OR.(LBCODEI(I).EQ.13)
     &           .OR.(LBCODEI(I).EQ.23)) THEN
               NNBB=NBV(I)
               RBARWL=2.D0*(ETA1(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN1(I)=-RampIntFlux1
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
               RBARWL=2.D0*(ETA2(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN2(I)=-RampIntFlux2
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
            ENDIF
         END DO
      ENDIF

C...
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
cjj   wm001 - modified/added the following 3 lines
C.... INTERNAL BARRIERS AND FOR INTERNAL BARRIER BOUNDARIES
C.... WITH CROSS BARRIER PIPES
C.... THIS SECTION ONLY RESTARTS THE OVER BARRIER FLOW COMPONENT
C...
cjj   wm001 - modified following line
      IF(NFLUXIB.EQ.1) THEN
         DO I=1,NVEL
cjj   wm001 - modified following 2 lines
            IF((LBCODEI(I).EQ.4).OR.(LBCODEI(I).EQ.24)
     &           .OR.(LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-BARINHT(I)
               RBARWL2=ETA1(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN1(I)=-RampIntFlux1*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1998
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN1(I)=-RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN1(I)=RampIntFlux1*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1998
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN1(I)=RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
 1998          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-BARINHT(I)
               RBARWL2=ETA2(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN2(I)=-RampIntFlux2*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1999
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN2(I)=-RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN2(I)=RampIntFlux2*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1999
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN2(I)=RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
 1999          CONTINUE
            ENDIF
         END DO
      ENDIF

cjj   wm001 - start add
C...
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
C.... INTERNAL BARRIERS WITH CROSS BARRIER PIPES
C.... THIS SECTION RESTARTS THE PIPE FLOW COMPONENT
C.... NOTE THAT PIPE FLOW COMPONENT IS ADDED INTO BARRIER FLOW COMPONENT
C.... THAT WAS PREVIOUSLY SET
C...
      IF(NFLUXIBP.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-PIPEHT(I)
               RBARWL2=ETA1(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
 2002          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-PIPEHT(I)
               RBARWL2=ETA2(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
 2003          CONTINUE
            ENDIF
         END DO
      ENDIF
cjjwm001 - end add

C-------------------TIME VARYING BATHYMETRY FIELDS-------------------------
C... update the time varying bathymetry if used
C... tcm v50.66.01 added this section
C

C     time varying bathymetry values are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the 
C     bathymetry information with the model time step.
      IF(Nddt.EQ.1) THEN
         OPEN(141,FILE=TRIM(INPUTDIR)//'/'//'fort.141')
         bTIME1 = STATIM*86400.D0
         bTIME2 = bTIME1 + bTIMINC
         BTIME_END = BTIME1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
         DP1(:) = DP(:)
         DP2(:) = DP(:)
         READ(141,*) (IJK,dp2(IJK),I=1,NP)
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
         IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
            DO I=1,NP
               IF (DP2(I).LT.H0) DP2(I) = H0
            ENDDO
         ENDIF
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.bTIME2) THEN
               bTIME1=bTIME2
               bTIME2=bTIME2+bTIMINC
               BTIME_END = BTIME1 + BCHGTIMINC
               dp1(:) = dp2(:)
               DO I=1,NP
                  READ(141,*) IJK,DP2(IJK)
               END DO
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
               IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
                  DO I=1,NP
                     IF (DP2(I).LT.H0) DP2(I) = H0
                  ENDDO
               ENDIF
            ENDIF  !timeit
         END DO  !IT

C.......If time is during the bathymetry change interval, then update bathymetry
         IF(TIMEloc.LE.BTIME_END) THEN
            bTRatio=(TIMEloc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I))  !Determine incremental amount to adjust bathymetry
               DP(I) = DP1(I) + DPTMP           !updating bathymetry
            ENDDO !I
         ENDIF
         IF(TIMEloc.EQ.BTIME_END) THEN
            DO I=1,NP
               DP(I) = DP2(I)         !updating bathymetry to final value
            ENDDO !I
         ENDIF
      ENDIF

      IF(NDDT.EQ.-1) THEN
         OPEN(141,FILE=TRIM(INPUTDIR)//'/'//'fort.141')
         bTIME1 = timeLoc
         bTIME2 = bTIME1 + bTIMINC
         BTIME_END = BTIME1 + BCHGTIMINC
         DP1(:) = DP(:)
         DP2(:) = DP(:)
         READ(141,*) (IJK,DP2(IJK),I=1,NP)
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
         IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
            DO I=1,NP
               IF (DP2(I).LT.H0) DP2(I) = H0
            ENDDO
         ENDIF
C.......If time is during the bathymetry change interval, then update bathymetry 
         IF(TIMEloc.LE.BTIME_END) THEN
            bTRatio=(Timeloc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I))  !Determine incremental amount to adjust bathymetry
               DP(I) = DP1(I) + DPTMP           !updating bathymetry
            ENDDO !I
         endif
         IF(TIMEloc.EQ.BTIME_END) THEN
            DO I=1,NP
               DP(I) = DP2(I)         !updating bathymetry to final value
            ENDDO !I
         ENDIF
      ENDIF

      IF(NDDT.EQ.2) THEN
         OPEN(141,FILE=TRIM(INPUTDIR)//'/'//'fort.141')
         bTIME1 = STATIM*86400.D0
         bTIME2 = bTIME1 + bTIMINC
         BTIME_END = BTIME1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
         DP1(:) = DP(:)
         DP2(:) = DP(:)

!!!     go get first record for only some nodes, all
!!!     other nodes keep their current value
        CALL NDDT2GET(141,DP2(:),-99999.d0 )
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
        IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
           DO I=1,NP
              IF (DP2(I).LT.H0) DP2(I) = H0
           ENDDO
        ENDIF

         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.bTIME2) THEN
               bTIME1=bTIME2
               bTIME2=bTIME2+bTIMINC
               BTIME_END = BTIME1 + BCHGTIMINC
               DP1(:) = DP2(:)
!!!            go get next record for only some nodes, all
!!!            other nodes keep their current value
               CALL NDDT2GET(141,DP2(:),-99999.d0 )
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
               IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
                  DO I=1,NP
                     IF (DP2(I).LT.H0) DP2(I) = H0
                  ENDDO
               ENDIF
            ENDIF
         END DO  !IT

C.......If time is during the bathymetry change interval, then update bathymetry
         IF(TIMEloc.LE.BTIME_END) THEN
         bTRatio=(Timeloc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I))  !Determine incremental amount to adjust bathymetry
               DP(I) = DP1(I) + DPTMP           !updating bathymetry
            ENDDO !I
         ENDIF
         IF(timeloc.EQ.BTIME_END) THEN
            DO I=1,NP
               DP(I) = DP2(I)         !updating bathymetry to final value
            ENDDO !I
         ENDIF
      ENDIF

      IF(NDDT.EQ.-2) THEN
         OPEN(141,FILE=TRIM(INPUTDIR)//'/'//'fort.141')
         bTIME1 = timeLoc
         bTIME2 = bTIME1 + bTIMINC
         BTIME_END = BTIME1 + BCHGTIMINC
         DP1(:) = DP(:)
         DP2(:) = DP(:)

!!!     go get first record for only some nodes, all
!!!     other nodes keep their current value
        CALL NDDT2GET(141,DP2(:),-99999.d0 )
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
         IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
            DO I=1,NP
               IF (DP2(I).LT.H0) DP2(I) = H0
            ENDDO
         ENDIF
C.......If time is during the bathymetry change interval, then update bathymetry
         IF(TIMEloc.LE.BTIME_END) THEN
            bTRatio=(Timeloc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I))  !Determine incremental amount to adjust bathymetry
               DP(I) = DP1(I) + DPTMP           !updating bathymetry
            ENDDO !I
         endif
         IF(TIMEloc.EQ.BTIME_END) THEN
            DO I=1,NP
               DP(I) = DP2(I)         !updating bathymetry to final value
            ENDDO !I
         ENDIF
      ENDIF


C------------------------ICE FIELDS----------------------------------------
C...  UPDATE THE ICE CONCENTRATION FIELDS FROM UNIT 25, 225, 227
c...  TCM V49.64.01 ADDED THE ICE FIELDS SECTION
C
      IF (NCICE.EQ.12) THEN
          CICE_TIME1 = STATIM*86400.D0
          CICE_TIME2 = CICE_TIME1 + CICE_TIMINC
          ! This just initializes some variables and might set CICE,etc...
          ! if skipping ahead is used
          CALL NCICE1_INIT(CICE1,NP)
          CALL NCICE1_GET(CICE1,NP)
          CALL NCICE1_GET(CICE2,NP)
          DO IT=1,ITHS
             TIMEIT=IT*DTDP + STATIM*86400.D0
             IF(TIMEIT.GT.CICE_TIME2) THEN
                CICE_TIME1 = CICE_TIME2
                CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
                DO I=1,NP
                   CICE1(I)=CICE2(I)
                END DO
                CALL NCICE1_GET(CICE2,NP)
             ENDIF
          END DO
      ENDIF
C
C

C------------------------MET FORCING---------------------------------------
C
C......RESTART WIND AND PRESSURE INFORMATION
C
C.... tcm v49.16 20100617
C.... convert background pressure from millibars to meters of water
      PRBCKGRND_MH2O = 100.0D0*PRBCKGRND/(RHOWAT0*G)

C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         DO J=1,ITHS
            DO I=1,NP
               READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I)
            END DO
         END DO
         DO I=1,NP
            WSX2(I)=RampMete2*WSX2(I)
            WSY2(I)=RampMete2*WSY2(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PR2(I)
            PRDIFF = RampMete2*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(NWS.EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WINDX
            WSY2(I) = RampMete2*WINDY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO I=1,NP
            WSX2(I) = RampMete2*WVNX1(I)
            WSY2(I) = RampMete2*WVNY1(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
 2223    CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED,NP,
     &        NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &        NScreen,ScreenUnit)
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIMEP
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIME
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIME
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,
     &              IWYR,WTIMED,NP,NWLON,NWLAT,WLATMAX,WLONMIN,
     &              WLATINC,WLONINC,ICS,NScreen,ScreenUnit)
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &              WRITE(ScreenUnit,*)
     &              'WIND FILE ADVANCED TO TIME',' = ', IWTIME
               WRITE(16,*) 'WIND FILE ADVANCED TO TIME = ',IWTIME
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP              !INTERPOLATE IN TIME
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG=SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I)=RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)=RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(NWS.EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  !TCM v49.02 Added read for rec. 1
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
               END DO
               CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(NWS.EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
C     jjw-42.06j wrote, jgf46.00 added following two lines:
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF

      IF(NWS.EQ.-5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC !jgf46.00 added
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WINDX !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
            WVNYOUT(I)=RampMete2*WINDY !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
C     The following 3 lines are a hardwire to allow a non standard met
C     file to be read in at time zero in a hot start.  They should be
C     eliminated or commented out for normal operation
c     OPEN(199,FILE=TRIM(INPUTDIR)//'/'//'fort.199')
c     READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c     CLOSE(199)
C     The following CALL statement should be uncommented for normal operation
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
      IF(NWS.EQ.7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)  !TCMv49.02 Changed (WVNX2,WVNY2,PRN2) to (WVNX1,WVNY1,PRN1)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

      IF(NWS.EQ.-7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
C     jgfdebug46.01 How do we convert these marine wind stresses into
C     directional land surface stresses?
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         HollandTime = TimeLoc
         CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF

!     rjw added nws = 19: asymmetric hurricane winds v2.0

      IF(ABS(NWS).EQ.19) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
             CALL NWS19GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TimeLoc, ICS)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF
C
C     jgf49.1001 Added NWS29, asymmetric vortex winds (NWS19)
C     embedded in an OWI basin scale gridded met field (NWS12)
C     from NAM data
      IF(NWS.EQ.29) THEN
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
C        This just initializes some variables and might
C        set WVNX1, etc ... if skipping ahead.
         CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
C        !TCM v49.02 Added read for rec. 1
         CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO IT=1,ITHS
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                END DO
                CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
             ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
            ! jgf49.1001 NAM winds already contain wind reduction
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*PRN1(I)
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF
      IF(NWS.EQ.-29) THEN
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
            ! jgf49.1001 NAM winds already contain wind reduction
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*PRN1(I)
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
         END DO
      ENDIF
C     Now blend in the vortex winds into the background OWI wind field.
      IF(ABS(NWS).EQ.29) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &        STATUS='OLD')
         ALLOCATE(vortexWVNX2(NP),vortexWVNY2(NP),vortexPRN2(NP))
         CALL NWS19GET(SLAM,SFEA,vortexWVNX2,vortexWVNY2,vortexPRN2,
     &                NP,TimeLoc, ICS)
         DO I=1,NP
            WindX = vortexWVNX2(I)
            WindY = vortexWVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            CALL getBlendFactor(I, SLAM, SFEA, bf)
            WSX2(I) = bf*RampMete2*0.001293d0*WDragCo*WindX*WindMag
     &               +(1.0d0-bf)*WSX2(I)
            WSY2(I) = bf*RampMete2*0.001293d0*WDragCo*WindY*WindMag
     &               +(1.0d0-bf)*WSY2(I)
            PR2(I) = bf*RampMete2*vortexPRN2(I)
     &               +(1.0d0-bf)*PRN2(I)
            WVNXOUT(I)=bf*RampMete2*WindX+(1.0d0-bf)*WVNXOUT(I)
            WVNYOUT(I)=bf*RampMete2*WindY+(1.0d0-bf)*WVNYOUT(I)
         ENDDO
      ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         NWSGGWI=-1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=0
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1197)
         WRITE(16,1197)
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !JUST COMPUTE INTERPOLATING FACTORS
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1198)
         WRITE(16,1198)
         NWSEGWI=0
         IDSETFLG=8
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &        RHOWAT0,G)        !NOW COMPUTE HOTSTART WIND FILED
         NWSEGWI=1
         IDSETFLG=1
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !NOW COMPUTE NEXT WIND FIELD
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C.....sb46.28sb01 NWS=12,-12 were added to deal with raw OWI files.  09/xx/2006
      IF(NWS.EQ.12) THEN
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  ! This just initializes some variables and might set WVNX1,etc... if skipping ahead.
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  !TCM v49.02 Added read for rec. 1
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO IT=1,ITHS
             TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
             IF(TIMEIT.GT.WTIME2) THEN
                WTIME1=WTIME2
                WTIME2=WTIME2+WTIMINC
                DO I=1,NP
                   WVNX1(I)=WVNX2(I)
                   WVNY1(I)=WVNY2(I)
                   PRN1(I)=PRN2(I)
                 END DO
                 CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
              ENDIF
          END DO
          WTRATIO=(TimeLoc-WTIME1)/WTIMINC
          DO I=1,NP
             WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
             WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, I)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I) = RampMete2*PRN1(I)
             PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF

       IF(NWS.EQ.-12) THEN
          WTIME1 = TimeLoc
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO I=1,NP
             WINDX = WVNX1(I)
             WINDY = WVNY1(I)
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, I)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
             PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX  ! tcm 20100618 v49.16 changed RampMete to RampMete2
             WVNYOUT(I)=RampMete2*WINDY  ! tcm 20100618 v49.16 changed RampMete to RampMete2
          END DO
       ENDIF

      !
      ! jgf50.38.05: Added NWS=15 for reading HWind data
       IF(ABS(NWS).EQ.15) THEN
          CALL NWS15INIT(timeloc)
          CALL NWS15GET(WVNX2,WVNY2,PRN2,timeloc)
          DO I=1,NP
             WDragCo = WindDrag(windSpeeds(i), i)
             windx = wvnx2(i)
             windy = wvny2(i)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(i, WDragCo,
     &              windSpeeds(i), DP(I), ETA2(I), H0, G,
     &              windx, windy)
                WindMag = SQRT(windx**2+windy**2)
                WDragCo = WindDrag(WindMag, i)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(i,windx,windy)
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
             PR2(I) = PRBCKGRND_MH2O + RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF


C
C--------------------END MET FORCING---------------------------------------

C......RESTART THE WAVE RADIATION STRESS

        IF(NRS.GE.1) THEN ! sb46.28sb03
          RSTIME1 = TimeLoc
          RSTIME2 = RSTIME1 + RSTIMINC
C          IF (MYPROC.EQ.0) THEN
C             PRINT *, "Made it to the radiation stress"
C             PRINT *, TIMELOC, RSTIME2, NRS
C          END IF
          ! NRS=2 was added. sb46.28sb03 09/xx/2006
          IF(NRS.EQ.1) THEN
            OPEN(23,FILE=TRIM(INPUTDIR)//'/'//'fort.23')
            CALL RSGET(RSNX1,RSNY1,NP)
            CALL RSGET(RSNX2,RSNY2,NP)
          ENDIF
          IF(NRS.EQ.2) THEN
            CALL RS2INIT(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX2,RSNY2,NP)
          ENDIF
          DO I=1,NP
            WSX2(I) = WSX2(I)+RampWRad2*RSNX1(I)
            WSY2(I) = WSY2(I)+RampWRad2*RSNY1(I)
          END DO
        ENDIF

       if (CTIP) then
Cjromo 11-01-00  Initialize TIP2 for HOTSTART
          DO I=1,NP
             TIP2(I)=0.0
          END DO
CTIP  LINES TO USE TIDAL POTENTIAL FORCING
       IF(NTIP.GE.1) THEN
         DO J=1,NTIF
           IF(PERT(J).EQ.0.) THEN
             NCYC=0
             ELSE
#ifdef IBM
             NCYC=INT(TimeLoc/PERT(J),KIND(0.0d0))
#else
             NCYC=INT(TimeLoc/PERT(J))
#endif
             ENDIF
           ARGT=AMIGT(J)*(TimeLoc-NCYC*PERT(J))+FACET(J)
           TPMUL=RampTip2*ETRF(J)*TPK(J)*FFT(J)
           SALTMUL=RampTip2*FFT(J)
#ifdef IBM
           NA=NINT(0.00014/AMIGT(J),KIND(0.0d0))
#else
           NA=NINT(0.00014/AMIGT(J))
#endif
           IF(NA.EQ.1) THEN                        !SEMI-DIURNAL SPECIES
             DO I=1,NP
               ARGTP=ARGT+2.*SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
               CCSFEA=COS(SFEA(I))
               CCSFEA=CCSFEA*CCSFEA
               TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           IF(NA.EQ.2) THEN
             DO I=1,NP
               ARGTP=ARGT+SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
cjjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
               S2SFEA=SIN(2.d0*SFEA(I))
#else
               S2SFEA=SIN(2.e0*SFEA(I))
#endif
               TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP)
     &                +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           END DO
         ENDIF
      endif     !   CTIP


C...
C....SET UP TO RESTART TIMESERIES OUTPUT FILES
C....
C...
      IF(NBYTE.EQ.4) ITEMPSTP=20
      IF(NBYTE.EQ.8) ITEMPSTP=10

C.... TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
C.... SET THE BATHY VARIABLES EQUAL TO THE ELVATION VALUES
C.... FOR USE LATER ON
      NTCYSB = NTCYSE
      NTCYFB = NTCYFE
      NSCOUB = NSCOUE
      IBSTP = IESTP

C...
C....IF RESTARTING THE ELEVATION STATION OUTPUT FILE, GO TO THE PROPER PLACE
C....IN THE FILE.  OTHERWISE ZERO OUT NSCOUE.
C...
      WRITE(16,1040) IESTP,NSCOUE
 1040 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN ELEVATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)

      IF(NOUTE.LT.0) THEN
        IESTP=0
        NSCOUE=0
        IF((NTCYSE.LT.ITHS).AND.(NSPOOLE.GT.0)) THEN
          NTCYSE=NTCYSE+((REAL(ITHS)-REAL(NTCYSE))/
     &           REAL(NSPOOLE))*REAL(NSPOOLE)
          IF(NTCYSE.LT.ITHS) NTCYSE=NTCYSE+NSPOOLE
          IF(NSPOOLE.NE.0) NTRSPE=(NTCYFE-NTCYSE)/NSPOOLE
        ENDIF
        WRITE(16,1041)
 1041   FORMAT(//,' A NEW ELEVATION STATION FILE WILL BE STARTED')
      ENDIF
      IF(NOUTE.EQ.-2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(61,REC=IESTP+I) RDES4(I)
          ENDDO
          IESTP=IESTP+8
          DO I=1,6
            WRITE(61,REC=IESTP+I) RID4(I)
          ENDDO
          IESTP=IESTP+6
          DO I=1,6
            WRITE(61,REC=IESTP+I) AID4(I)
          ENDDO
          IESTP=IESTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(61,REC=IESTP+I) RDES8(I)
          ENDDO
          IESTP=IESTP+4
          DO I=1,3
            WRITE(61,REC=IESTP+I) RID8(I)
          ENDDO
          IESTP=IESTP+3
          DO I=1,3
            WRITE(61,REC=IESTP+I) AID8(I)
          ENDDO
          IESTP=IESTP+3
        ENDIF
        WRITE(61,REC=IESTP+1) NTRSPE
        WRITE(61,REC=IESTP+2) NSTAE
        WRITE(61,REC=IESTP+3) DT*NSPOOLE
        WRITE(61,REC=IESTP+4) NSPOOLE
        WRITE(61,REC=IESTP+5) 1
        IESTP=IESTP+5
        CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTE.EQ.-1) THEN
        CALL OPEN_GBL_FILE(61, TRIM(GLOBALDIR)//'/'//'fort.61',
     $    NSTAE_G, NSTAE, HEADER61)
        IESTP=2
      ENDIF

      IF(NOUTE.EQ.2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(61,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(61)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE VELOCITY STATION OUTPUT FILE
C...
      WRITE(16,1042) IVSTP,NSCOUV
 1042 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN VELOCITY ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER =',I6)

      IF(NOUTV.LT.0) THEN
        IVSTP=0
        NSCOUV=0
        IF((NTCYSV.LT.ITHS).AND.(NSPOOLV.GT.0)) THEN
          NTCYSV=NTCYSV+((REAL(ITHS)-REAL(NTCYSV))/
     &           REAL(NSPOOLV))*REAL(NSPOOLV)
          IF(NTCYSV.LT.ITHS) NTCYSV=NTCYSV+NSPOOLV
          NTRSPV=(NTCYFV-NTCYSV)/NSPOOLV
        ENDIF
        WRITE(16,1043)
 1043   FORMAT(//,' A NEW VELOCITY STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTV.EQ.-2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(62,REC=IVSTP+I) RDES4(I)
          ENDDO
          IVSTP=IVSTP+8
          DO I=1,6
            WRITE(62,REC=IVSTP+I) RID4(I)
          ENDDO
          IVSTP=IVSTP+6
          DO I=1,6
            WRITE(62,REC=IVSTP+I) AID4(I)
          ENDDO
          IVSTP=IVSTP+6
          ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(62,REC=IVSTP+I) RDES8(I)
          ENDDO
          IVSTP=IVSTP+4
          DO I=1,3
            WRITE(62,REC=IVSTP+I) RID8(I)
          ENDDO
          IVSTP=IVSTP+3
          DO I=1,3
            WRITE(62,REC=IVSTP+I) AID8(I)
          ENDDO
          IVSTP=IVSTP+3
        ENDIF
        WRITE(62,REC=IVSTP+1) NTRSPV
        WRITE(62,REC=IVSTP+2) NSTAV
        WRITE(62,REC=IVSTP+3) DT*NSPOOLV
        WRITE(62,REC=IVSTP+4) NSPOOLV
        WRITE(62,REC=IVSTP+5) 2
        IVSTP=IVSTP+5
        CLOSE(62)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTV.EQ.-1) THEN
        CALL OPEN_GBL_FILE(62, TRIM(GLOBALDIR)//'/'//'fort.62',
     $    NSTAV_G, NSTAV, HEADER62)
        IVSTP=2
      ENDIF
      IF(NOUTV.EQ.2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(62,REC=ITEMPSTP+1) NTRSPV   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(62)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE CONCENTRATION STATION OUTPUT FILE
C...
      WRITE(16,1044) ICSTP,NSCOUC
 1044 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN CONCENTRATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
      IF(NOUTC.LT.0) THEN
        ICSTP=0
        NSCOUC=0
        IF((NTCYSC.LT.ITHS).AND.(NSPOOLC.GT.0)) THEN
          NTCYSC=NTCYSC+((REAL(ITHS)-REAL(NTCYSC))/
     &           REAL(NSPOOLC))*REAL(NSPOOLC)
          IF(NTCYSC.LT.ITHS) NTCYSC=NTCYSC+NSPOOLC
          NTRSPC=(NTCYFC-NTCYSC)/NSPOOLC
        ENDIF
        WRITE(16,1045)
 1045   FORMAT(//,' A NEW CONCENTRATION STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTC.EQ.-2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(81,REC=ICSTP+I) RDES4(I)
          ENDDO
          ICSTP=ICSTP+8
          DO I=1,6
            WRITE(81,REC=ICSTP+I) RID4(I)
          ENDDO
          ICSTP=ICSTP+6
          DO I=1,6
            WRITE(81,REC=ICSTP+I) AID4(I)
          ENDDO
          ICSTP=ICSTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(81,REC=ICSTP+I) RDES8(I)
          ENDDO
          ICSTP=ICSTP+4
          DO I=1,3
            WRITE(81,REC=ICSTP+I) RID8(I)
          ENDDO
          ICSTP=ICSTP+3
          DO I=1,3
            WRITE(81,REC=ICSTP+I) AID8(I)
          ENDDO
          ICSTP=ICSTP+3
        ENDIF
        WRITE(81,REC=ICSTP+1) NTRSPC
        WRITE(81,REC=ICSTP+2) NSTAC
        WRITE(81,REC=ICSTP+3) DT*NSPOOLC
        WRITE(81,REC=ICSTP+4) NSPOOLC
        WRITE(81,REC=ICSTP+5) 1
        ICSTP=ICSTP+5
        CLOSE(81)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTC.EQ.-1) THEN
        CALL OPEN_GBL_FILE(81, TRIM(GLOBALDIR)//'/'//'fort.81',
     $    NSTAC_G, NSTAC, HEADER81)
        ICSTP=2
        CLOSE(81)
      ENDIF
      IF(NOUTC.EQ.2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(81,REC=ITEMPSTP+1) NTRSPC   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(81)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE METEOROLOGICAL STATION OUTPUT FILE
C...
c... v49.64.01 tcm -added support for ice stations
      IF(NCICE.EQ.0)  THEN
         WRITE(16,1038) IWSTP,IPSTP,NSCOUM
      ELSE
         IICESTP = IPSTP  !SET TO THE SAME PLACE AS THE PRESSURE FILE
         WRITE(16,2038) IWSTP,IPSTP,IICESTP,NSCOUM
      ENDIF
      WRITE(16,1038) IWSTP,IPSTP,NSCOUM
 1038 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTIN IN THE PRES STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
 2038 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTEN IN THE PRES STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTEN IN THE ICE STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
C
      IF(NOUTM.LT.0) THEN
        IPSTP=0
        IWSTP=0
        IICESTP = 0  !v49.64.01 tcm added ice station
        NSCOUM=0
        IF((NTCYSM.LT.ITHS).AND.(NSPOOLM.GT.0)) THEN
          NTCYSM=NTCYSM+((REAL(ITHS)-REAL(NTCYSM))/
     &           REAL(NSPOOLM))*REAL(NSPOOLM)
          IF(NTCYSM.LT.ITHS) NTCYSM=NTCYSM+NSPOOLM
          NTRSPM=(NTCYFM-NTCYSM)/NSPOOLM
        ENDIF
        WRITE(16,1039)
 1039   FORMAT(//,' A NEW METEOROLOGICAL STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTM.EQ.-2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
!       v49.64.01 tcm added ice station
        IF(NCICE.NE.0) THEN
           OPEN(91,FILE=TRIM(LOCALDIR)//'/'//'fort.91',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        ENDIF
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(71,REC=IPSTP+I) RDES4(I)
            WRITE(72,REC=IWSTP+I) RDES4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RDES4(I)
          ENDDO
          IPSTP=IPSTP+8
          IWSTP=IWSTP+8
          IF (NCICE.NE.0) IICESTP = IICESTP+8
          DO I=1,6
            WRITE(71,REC=IPSTP+I) RID4(I)
            WRITE(72,REC=IWSTP+I) RID4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
          IF (NCICE.NE.0) IICESTP = IICESTP+6
          DO I=1,6
            WRITE(71,REC=IPSTP+I) AID4(I)
            WRITE(72,REC=IWSTP+I) AID4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) AID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
          IF (NCICE.NE.0) IICESTP = IICESTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(71,REC=IPSTP+I) RDES8(I)
            WRITE(72,REC=IWSTP+I) RDES8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RDES8(I)
          ENDDO
          IPSTP=IPSTP+4
          IWSTP=IWSTP+4
          IF (NCICE.NE.0) IICESTP = IICESTP+4
          DO I=1,3
            WRITE(71,REC=IPSTP+I) RID8(I)
            WRITE(72,REC=IWSTP+I) RID8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RID8(I)
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
          IF (NCICE.NE.0) IICESTP = IICESTP+3
          DO I=1,3
            WRITE(71,REC=IPSTP+I) AID8(I)
            WRITE(72,REC=IWSTP+I) AID8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) AID8(I)
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
          IF (NCICE.NE.0) IICESTP = IICESTP+3
        ENDIF
        WRITE(71,REC=IPSTP+1) NTRSPM
        WRITE(71,REC=IPSTP+2) NSTAM
        WRITE(71,REC=IPSTP+3) DT*NSPOOLM
        WRITE(71,REC=IPSTP+4) NSPOOLM
        WRITE(71,REC=IPSTP+5) 1
        WRITE(72,REC=IWSTP+1) NTRSPM
        WRITE(72,REC=IWSTP+2) NSTAM
        WRITE(72,REC=IWSTP+3) DT*NSPOOLM
        WRITE(72,REC=IWSTP+4) NSPOOLM
        WRITE(72,REC=IWSTP+5) 2
        IPSTP=IPSTP+5
        IWSTP=IWSTP+5
        CLOSE(71)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
        ! v49.64.01 tcm -- added ice stations
        IF(NCICE.NE.0) THEN
           WRITE(91,REC=IICESTP+1) NTRSPM
           WRITE(91,REC=IICESTP+2) NSTAM
           WRITE(91,REC=IICESTP+3) DT*NSPOOLM
           WRITE(91,REC=IICESTP+4) NSPOOLM
           WRITE(91,REC=IICESTP+5) 1
           IICESTP = IICESTP + 5
           CLOSE(91)
        ENDIF
      ENDIF
      IF(NOUTM.EQ.-1) THEN
        CALL OPEN_GBL_FILE(71, TRIM(GLOBALDIR)//'/'//'fort.71',
     $    NSTAM_G, NSTAM, HEADER71)
        IPSTP=2
        CALL OPEN_GBL_FILE(72, TRIM(GLOBALDIR)//'/'//'fort.72',
     $    NSTAM_G, NSTAM, HEADER72)
        IWSTP=2
!       v49.64.01 tcm added ice station support
        IF(NCICE.NE.0) THEN
           CALL OPEN_GBL_FILE(91, TRIM(GLOBALDIR)//'/'//'fort.91',
     $       NSTAM_G, NSTAM, HEADER91)
           IICESTP=2
        ENDIF
      ENDIF
      IF(NOUTM.EQ.2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(71,REC=ITEMPSTP+1) NTRSPM   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
c         WRITE(72,REC=ITMEPSTP+1) NTRSPM
        CLOSE(71)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
!       v49.64.01 tcm added ice station support
        IF(NCICE.NE.0) THEN
           OPEN(91,FILE=TRIM(LOCALDIR)//'/'//'fort.91',
     &                      ACCESS='DIRECT',RECL=NBYTE)
           CLOSE(91)
        ENDIF
      ENDIF


C... TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
C... THIS FILE IS CONTROLLED BY THE ELEVATION STATION DATA
C...
C....IF RESTARTING THE BATHYMETRY STATION OUTPUT FILE, GO TO THE PROPER PLACE
C....IN THE FILE.  OTHERWISE ZERO OUT NSCOUE.
C...

      IF (NDDT.NE.0) WRITE(16,1035) IBSTP,NSCOUB
 1035 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN BATHYMETRY ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)

      IF ((NDDT.NE.0).AND.(NOUTE.LT.0)) THEN
        IBSTP=0
        NSCOUB=0
        IF((NTCYSB.LT.ITHS).AND.(NSPOOLE.GT.0)) THEN
          NTCYSB=NTCYSB+((REAL(ITHS)-REAL(NTCYSB))/
     &           REAL(NSPOOLE))*REAL(NSPOOLE)
          IF(NTCYSB.LT.ITHS) NTCYSB=NTCYSB+NSPOOLE
          IF(NSPOOLE.NE.0) NTRSPB=(NTCYFB-NTCYSB)/NSPOOLE
        ENDIF
        WRITE(16,1036)
 1036   FORMAT(//,' A NEW BATHYMETRY STATION FILE WILL BE STARTED')
      ENDIF
      IF(NOUTE.EQ.-2) THEN
        OPEN(75,FILE=TRIM(LOCALDIR)//'/'//'fort.75',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(75,REC=IBSTP+I) RDES4(I)
          ENDDO
          IBSTP=IBSTP+8
          DO I=1,6
            WRITE(75,REC=IBSTP+I) RID4(I)
          ENDDO
          IBSTP=IBSTP+6
          DO I=1,6
            WRITE(75,REC=IBSTP+I) AID4(I)
          ENDDO
          IBSTP=IBSTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(75,REC=IBSTP+I) RDES8(I)
          ENDDO
          IBSTP=IBSTP+4
          DO I=1,3
            WRITE(75,REC=IBSTP+I) RID8(I)
          ENDDO
          IBSTP=IBSTP+3
          DO I=1,3
            WRITE(75,REC=IBSTP+I) AID8(I)
          ENDDO
          IBSTP=IBSTP+3
        ENDIF
        WRITE(75,REC=IBSTP+1) NTRSPB
        WRITE(75,REC=IBSTP+2) NSTAE
        WRITE(75,REC=IBSTP+3) DT*NSPOOLE
        WRITE(75,REC=IBSTP+4) NSPOOLE
        WRITE(75,REC=IBSTP+5) 1
        IBSTP=IBSTP+5
        CLOSE(75)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTE.EQ.-1) THEN
        CALL OPEN_GBL_FILE(75, TRIM(GLOBALDIR)//'/'//'fort.75',
     $    NSTAE_G, NSTAE, HEADER61)
        IESTP=2
      ENDIF

      IF ((NDDT.NE.0).AND.(NOUTE.EQ.2)) THEN
        OPEN(75,FILE=TRIM(LOCALDIR)//'/'//'fort.75',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(75,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(75)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL ELEVATION OUTPUT FILE
C...
        WRITE(16,1046) IGEP,NSCOUGE
 1046   FORMAT(//,1X,I8,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ELEVATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I8)

C....   TCM V50.66.01 ADDITIONS FOR TIME VARYING BATYMETRY
C....   SET THESE VARIABLES FOR USE LATER 
        IGBP = IGEP
        NTCYSGB = NTCYSGE
        NTCYFGB = NTCYFGE
        NSCOUGB = NSCOUGE

        IF(NOUTGE.LT.0) THEN
          IGEP=0
          NSCOUGE=0
          IF((NTCYSGE.LT.ITHS).AND.(NSPOOLGE.GT.0)) THEN
            NTCYSGE=NTCYSGE+((REAL(ITHS)-REAL(NTCYSGE))/
     &              REAL(NSPOOLGE))*REAL(NSPOOLGE)
            IF(NTCYSGE.LT.ITHS) NTCYSGE=NTCYSGE+NSPOOLGE
C     kmd48.33bc - changed NDSETSE for use of new DTDPHS
            IF(DTDP.NE.DTDPHS) THEN
              NDSETSE=((NTCYFGE*DTDP)-(NTCYSGE*DTDPHS))/
     &                (NSPOOLGE*DTDP)
            ELSE
              NDSETSE=(NTCYFGE-NTCYSGE)/NSPOOLGE
            END IF
            ENDIF
          WRITE(16,1047)
 1047     FORMAT(//,' A NEW GLOBAL ELEVATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGE.EQ.-2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=IGEP+I) RDES4(I)
              ENDDO
            IGEP=IGEP+8
            DO I=1,6
              WRITE(63,REC=IGEP+I) RID4(I)
              ENDDO
            IGEP=IGEP+6
            DO I=1,6
              WRITE(63,REC=IGEP+I) AID4(I)
              ENDDO
            IGEP=IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=IGEP+I) RDES8(I)
              ENDDO
            IGEP=IGEP+4
            DO I=1,3
              WRITE(63,REC=IGEP+I) RID8(I)
              ENDDO
            IGEP=IGEP+3
            DO I=1,3
              WRITE(63,REC=IGEP+I) AID8(I)
              ENDDO
            IGEP=IGEP+3
            ENDIF
          WRITE(63,REC=IGEP+1) NDSETSE
          WRITE(63,REC=IGEP+2) NP
          WRITE(63,REC=IGEP+3) DT*NSPOOLGE
          WRITE(63,REC=IGEP+4) NSPOOLGE
          WRITE(63,REC=IGEP+5) 1
          IGEP=IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

Casey 090310: Added the possibility for NOUTGE = -4.
        IF((NOUTGE.EQ.-1).OR.(NOUTGE.EQ.-4)) THEN
           CALL OPEN_GBL_FILE(63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $         NP_G, NP, HEADER63)
C
C        jgf47.06 Tau0 output is produced on the same schedule as
C        global elevation
           IF (OUTPUTTAU0) THEN
Casey 090310: Corrected the unit number from 10 to 90.
               CALL writeDomainHeader(90,
     &              TRIM(GLOBALDIR)//'/'//'fort.tau0',
     $               NP_G, NP, 'Tau0      ')
           ENDIF

C kmd48.33bc - information for writing out the sponge layer
           IF(OUTPUTSPONGE) THEN
              CALL writeDomainHeader(92,
     &             TRIM(GLOBALDIR)//'/'//'fort.92',
     $             NP_G, NP, 'sponge      ')
           END IF

           IGEP=2
        ENDIF
        IF(NOUTGE.EQ.2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(63,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(63)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL VELOCITY OUTPUT FILE
C...
        WRITE(16,1048) IGVP,NSCOUGV
 1048   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'VELOCITY FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGV.LT.0) THEN
          IGVP=0
          NSCOUGV=0
          IF((NTCYSGV.LT.ITHS).AND.(NSPOOLGV.GT.0)) THEN
            NTCYSGV=NTCYSGV+((REAL(ITHS)-REAL(NTCYSGV))/
     &              REAL(NSPOOLGV))*REAL(NSPOOLGV)
            IF(NTCYSGV.LT.ITHS) NTCYSGV=NTCYSGV+NSPOOLGV
C  kmd48.33bc - change the NDSETSV for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN
              NDSETSV=((NTCYFGV*DTDP)-(NTCYSGV*DTDPHS))/
     &                (NSPOOLGV*DTDP)
            ELSE
              NDSETSV=(NTCYFGV-NTCYSGV)/NSPOOLGV
            END IF
            ENDIF
          WRITE(16,1049)
 1049     FORMAT(//,' A NEW GLOBAL VELOCITY FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGV.EQ.-2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=IGVP+I) RDES4(I)
              ENDDO
            IGVP=IGVP+8
            DO I=1,6
              WRITE(64,REC=IGVP+I) RID4(I)
              ENDDO
            IGVP=IGVP+6
            DO I=1,6
              WRITE(64,REC=IGVP+I) AID4(I)
              ENDDO
            IGVP=IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=IGVP+I) RDES8(I)
              ENDDO
            IGVP=IGVP+4
            DO I=1,3
              WRITE(64,REC=IGVP+I) RID8(I)
              ENDDO
            IGVP=IGVP+3
            DO I=1,3
              WRITE(64,REC=IGVP+I) AID8(I)
              ENDDO
            IGVP=IGVP+3
            ENDIF
          WRITE(64,REC=IGVP+1) NDSETSV
          WRITE(64,REC=IGVP+2) NP
          WRITE(64,REC=IGVP+3) DT*NSPOOLGV
          WRITE(64,REC=IGVP+4) NSPOOLGV
          WRITE(64,REC=IGVP+5) 2
          IGVP=IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF((NOUTGV.EQ.-1).OR.(NOUTGV.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $      NP_G, NP, HEADER64)
          IGVP=2
        ENDIF
        IF(NOUTGV.EQ.2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(64,REC=ITEMPSTP+1) NDSETSV  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(64)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL CONCENTRATION OUTPUT FILE
C...
        WRITE(16,1053) IGCP,NSCOUGC
 1053   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                'CONCENTRATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)

        IF(NOUTGC.LT.0) THEN
          IGCP=0
          NSCOUGC=0
          IF((NTCYSGC.LT.ITHS).AND.(NSPOOLGC.GT.0)) THEN
            NTCYSGC=NTCYSGC+((REAL(ITHS)-REAL(NTCYSGC))/
     &              REAL(NSPOOLGC))*REAL(NSPOOLGC)
            IF(NTCYSGC.LT.ITHS) NTCYSGC=NTCYSGC+NSPOOLGC
C  kmd48.33bc - change the NDSETSC for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN
              NDSETSC=((NTCYFGC*DTDP)-(NTCYSGC*DTDPHS))/
     &                (NSPOOLGC*DTDP)
            ELSE
              NDSETSC=(NTCYFGC-NTCYSGC)/NSPOOLGC
            END IF
            ENDIF
          WRITE(16,1054)
 1054     FORMAT(//,' A NEW GLOBAL CONCENTRATION FILE WILL BE STARTED')
        ENDIF

        IF(NOUTGC.EQ.-2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &           ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(83,REC=IGCP+I) RDES4(I)
              ENDDO
            IGCP=IGCP+8
            DO I=1,6
              WRITE(83,REC=IGCP+I) RID4(I)
              ENDDO
            IGCP=IGCP+6
            DO I=1,6
              WRITE(83,REC=IGCP+I) AID4(I)
              ENDDO
            IGCP=IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(83,REC=IGCP+I) RDES8(I)
              ENDDO
            IGCP=IGCP+4
            DO I=1,3
              WRITE(83,REC=IGCP+I) RID8(I)
              ENDDO
            IGCP=IGCP+3
            DO I=1,3
              WRITE(83,REC=IGCP+I) AID8(I)
              ENDDO
            IGCP=IGCP+3
            ENDIF
          WRITE(83,REC=IGCP+1) NDSETSC
          WRITE(83,REC=IGCP+2) NP
          WRITE(83,REC=IGCP+3) DT*NSPOOLGC
          WRITE(83,REC=IGCP+4) NSPOOLGC
          WRITE(83,REC=IGCP+5) 1
          IGCP=IGCP+5
          CLOSE(83)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF(NOUTGC.EQ.-1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          WRITE(83,3220) RUNDES,RUNID,AGRID
          WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
          IGCP=2
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          DO I=1,IGCP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSC
            READ(83,1050)
 1050       FORMAT(1X)
          ENDDO
          ENDFILE(83)
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(83,REC=ITEMPSTP+1) NDSETSC  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(83)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL METEOROLOGICAL OUTPUT FILES
C...
!       V49.64.01 TCM ADDED GLOBAL ICE CONCENTRATION FILES
        IF(NCICE.NE.0) THEN
           WRITE(16,2055) IGWP,IGPP,IGIP,NSCOUGW
        ELSE
           WRITE(16,1055) IGWP,IGPP,NSCOUGW
        ENDIF
 1055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
 2055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ICE FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
c
        IF(NOUTGW.LT.0) THEN
          IGPP=0
          IGWP=0
          NSCOUGW=0
          IGIP=0
          IF((NTCYSGW.LT.ITHS).AND.(NSPOOLGW.GT.0)) THEN
            NTCYSGW=NTCYSGW+((REAL(ITHS)-REAL(NTCYSGW))/
     &              REAL(NSPOOLGW))*REAL(NSPOOLGW)
            IF(NTCYSGW.LT.ITHS) NTCYSGW=NTCYSGW+NSPOOLGW
C  kmd48.33bc - change the NDSETSC for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN
              NDSETSW=((NTCYFGW*DTDP)-(NTCYSGW*DTDPHS))/
     &                (NSPOOLGW*DTDP)
            ELSE
              NDSETSW=(NTCYFGW-NTCYSGW)/NSPOOLGW
            END IF
          ENDIF
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             WRITE(16,2056)
          ELSE
             WRITE(16,1056)
          ENDIF
          WRITE(16,1056)
 1056     FORMAT(//,' NEW GLOBAL WIND & pressure FILEs WILL BE STARTED')
 2056   FORMAT(//,' NEW GLOBAL WIND & pressure & Ice FILEs WILL BE STARTED')
        ENDIF

        IF(NOUTGW.EQ.-2) THEN
          OPEN(73,file=trim(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             OPEN(93,FILE=TRIM(LOCALDIR)//'/'//'fort.93',
     &          ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(73,rec=igpp+i) rdes4(i)
              WRITE(74,REC=IGWP+I) RDES4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rdes4(i)
              ENDDO
            igpp=igpp+8
            IGWP=IGWP+8
            if(ncice.ne.0) igip = igip+8
            DO I=1,6
              write(73,rec=igpp+i) rid4(i)
              WRITE(74,REC=IGWP+I) RID4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
            DO I=1,6
              write(73,rec=igpp+i) aid4(i)
              WRITE(74,REC=IGWP+I) AID4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) aid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
          ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(73,rec=igpp+i) rdes8(i)
              WRITE(74,REC=IGWP+I) RDES8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rdes8(i)
            ENDDO
            igpp=igpp+4
            IGWP=IGWP+4
            if(ncice.ne.0) igip = igip+4
            DO I=1,3
              write(73,rec=igpp+i) rid8(i)
              WRITE(74,REC=IGWP+I) RID8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
            DO I=1,3
              write(73,rec=igpp+i) aid8(i)
              WRITE(74,REC=IGWP+I) AID8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) aid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
          ENDIF
          WRITE(73,rec=igpp+1) ndsetsw
          WRITE(73,rec=igpp+2) np
          WRITE(73,rec=igpp+3) dt*nspoolgw
          WRITE(73,rec=igpp+4) nspoolgw
          WRITE(73,rec=igpp+5) 2
          IGPP=IGPP+5
          CLOSE(73)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          WRITE(74,REC=IGWP+1) NDSETSW
          WRITE(74,REC=IGWP+2) NP
          WRITE(74,REC=IGWP+3) DT*NSPOOLGW
          WRITE(74,REC=IGWP+4) NSPOOLGW
          WRITE(74,REC=IGWP+5) 2
          IGWP=IGWP+5
          CLOSE(74)                    ! DO THIS TO FLUSH THE WRITE BUFFER
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             WRITE(93,REC=IGIP+1) NDSETSW
             WRITE(93,REC=IGIP+2) NP
             WRITE(93,REC=IGIP+3) DT*NSPOOLGW
             WRITE(93,REC=IGIP+4) NSPOOLGW
             WRITE(93,REC=IGIP+5) 2
             IGIP=IGIP+5
             CLOSE(93)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          ENDIF
        ENDIF

Casey 090310: Added the possibility for NOUTGW = -4.
        IF((NOUTGW.EQ.-1).OR.(NOUTGW.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $      NP_G, NP, HEADER73)
          IGPP=2
          CALL OPEN_GBL_FILE(74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $      NP_G, NP, HEADER74)
          IGWP=2
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF (NCICE.NE.0) THEN
            CALL OPEN_GBL_FILE(93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $        NP_G, NP, HEADER93)
            IGIP=2
          ENDIF
        ENDIF

! tcm v50.75 removed ifdef cswan and added test for nrs=3 or 4
        IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
Casey 090310: Added the output of radiation stress gradients.
           IF(NOUTGW.LT.0)THEN
              IGRadS=0
           ENDIF
           IF(NOUTGW.EQ.-2)THEN
              OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &            ACCESS='DIRECT',RECL=NBYTE)
              IF(NBYTE.EQ.4) THEN
                 DO I=1,8
                    WRITE(164,REC=IGRadS+I) RDES4(I)
                 ENDDO
                 IGRadS=IGRadS+8
                 DO I=1,6
                    WRITE(164,REC=IGRadS+I) RID4(I)
                 ENDDO
                 IGRadS=IGRadS+6
                 DO I=1,6
                    WRITE(164,REC=IGRadS+I) AID4(I)
                 ENDDO
                 IGRadS=IGRadS+6
              ENDIF
              IF(NBYTE.EQ.8) THEN
                 DO I=1,4
                    WRITE(164,REC=IGRadS+I) RDES8(I)
                 ENDDO
                 IGRadS=IGRadS+4
                 DO I=1,3
                    WRITE(164,REC=IGRadS+I) RID8(I)
                 ENDDO
                 IGRadS=IGRadS+3
                 DO I=1,3
                    WRITE(164,REC=IGRadS+I) AID8(I)
                 ENDDO
                 IGRadS=IGRadS+3
              ENDIF
              WRITE(164,REC=IGRadS+1) NDSETSW
              WRITE(164,REC=IGRadS+2) NP
              WRITE(164,REC=IGRadS+3) DT*NSPOOLGW
              WRITE(164,REC=IGRadS+4) NSPOOLGW
              WRITE(164,REC=IGRadS+5) 2
              IGRadS=IGRadS+5
              CLOSE(164)                    ! DO THIS TO FLUSH THE WRITE BUFFER
           ENDIF
           IF((NOUTGW.EQ.-1).OR.(NOUTGW.EQ.-4)) THEN
              CALL OPEN_GBL_FILE(164, TRIM(GLOBALDIR)//'/'//'rads.64',
     &            NP_G, NP, HEADER74)
              IGRadS=2
           ENDIF
        ENDIF
!#endif

        IF(NOUTGW.EQ.2) THEN
          open(73,file=TRIM(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
c         WRITE(73,REC=ITempStp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          close(73)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(74,REC=ITEMPSTP+1) NDSETSW  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(74)                         ! DO THIS TO FLUSH THE WRITE BUFFER
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             OPEN(93,FILE=TRIM(LOCALDIR)//'/'//'fort.93',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(93,REC=ITempStp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
             CLOSE(93)                      ! DO THIS TO FLUSH THE WRITE BUFFER
          ENDIF
        ENDIF

C... TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
C... THIS FILE IS CONTROLED BY THE ELEVATION OUTPUT OPTIONS
C....GO TO THE PROPER PLACE IN THE GLOBAL BATHYMETRY OUTPUT FILE
C...
        IF (NDDT.NE.0) WRITE(16,1051) IGBP,NSCOUGB
 1051   FORMAT(//,1X,I8,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'BATHYMETRY FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I8)

        IF ((NDDT.NE.0).AND.(NOUTGE.LT.0)) THEN
           IGBP=0
           NSCOUGB=0
           IF((NTCYSGB.LT.ITHS).AND.(NSPOOLGE.GT.0)) THEN
              NTCYSGB=NTCYSGB+((REAL(ITHS)-REAL(NTCYSGB))/
     &               REAL(NSPOOLGE))*REAL(NSPOOLGE)
              IF(NTCYSGB.LT.ITHS) NTCYSGB=NTCYSGB+NSPOOLGE
              IF(DTDP.NE.DTDPHS) THEN
                 NDSETSB=((NTCYFGB*DTDP)-(NTCYSGB*DTDPHS))/
     &                (NSPOOLGE*DTDP)
              ELSE
                 NDSETSB=(NTCYFGB-NTCYSGB)/NSPOOLGE
              ENDIF
           ENDIF
           WRITE(16,1052)
 1052      FORMAT(//,' A NEW GLOBAL BATHYMETRY FILE WILL BE STARTED')
        ENDIF

        IF ((NDDT.NE.0).AND.(NOUTGE.EQ.-2)) THEN
           OPEN(76,FILE=TRIM(LOCALDIR)//'/'//'fort.76',
     &          ACCESS='DIRECT',RECL=NBYTE)
           IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(76,REC=IGBP+I) RDES4(I)
              ENDDO
              IGBP=IGBP+8
              DO I=1,6
                 WRITE(76,REC=IGBP+I) RID4(I)
              ENDDO
              IGBP=IGBP+6
              DO I=1,6
                 WRITE(76,REC=IGBP+I) AID4(I)
              ENDDO
              IGBP=IGBP+6
           ENDIF
           IF(NBYTE.EQ.8) THEN
              DO I=1,4
                 WRITE(76,REC=IGBP+I) RDES8(I)
              ENDDO
              IGBP=IGBP+4
              DO I=1,3
                 WRITE(76,REC=IGBP+I) RID8(I)
              ENDDO
              IGBP=IGBP+3
              DO I=1,3
                 WRITE(76,REC=IGBP+I) AID8(I)
              ENDDO
              IGBP=IGBP+3
           ENDIF
           WRITE(76,REC=IGBP+1) NDSETSE
           WRITE(76,REC=IGBP+2) NP
           WRITE(76,REC=IGBP+3) DT*NSPOOLGE
           WRITE(76,REC=IGBP+4) NSPOOLGE
           WRITE(76,REC=IGBP+5) 1
           IGBP=IGBP+5
           CLOSE(76)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

        IF ((NDDT.NE.0).AND.
     &            ((NOUTGE.EQ.-1).OR.(NOUTGE.EQ.-4))) THEN
           CALL OPEN_GBL_FILE(76, TRIM(GLOBALDIR)//'/'//'fort.76',
     &         NP_G, NP, HEADER63)
           IGBP=2
        ENDIF
        IF ((NDDT.NE.0).AND.(NOUTGE.EQ.2)) THEN
           OPEN(76,FILE=TRIM(LOCALDIR)//'/'//'fort.76',
     &          ACCESS='DIRECT',RECL=NBYTE)
c          WRITE(76,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
           CLOSE(76)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF


C
C     jgf48.4635 If the maxele.63 (etc) files are present, we need to
C     load them up, to preserve high water marks across hotstarts.
C     If they are not present, thats okay too.
      
      IF(LoadEleSlopeLim)THEN
        IF(MNPROC.GT.1)ALLOCATE(HOT_ESLONOFF_G(NP_G))
        ALLOCATE(HOT_ESLONOFF(NP))
      ENDIF  
      
      if ((MNPROC.gt.1).and.(MyProc.eq.0)
     &     .and.(.not.WRITE_LOCAL_FILES)) then
         ALLOCATE(HOT_ETAMAX_G(NP_G))
         HOT_ETAMAX_G(:) = -99999.d0  ! initialize the value
         CALL readMinMax('maxele.63',HOT_ETAMAX_G)
         ALLOCATE(HOT_UMAX_G(NP_G))
         HOT_UMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxvel.63',HOT_UMAX_G)
         ALLOCATE(HOT_WVNOUTMAX_G(NP_G))
         HOT_WVNOUTMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxwvel.63',HOT_WVNOUTMAX_G)
         ALLOCATE(HOT_PrMin_G(NP_G))
         HOT_PrMin_G(:) = 99999.d0  ! initialize the value
         CALL readMinMax('minpr.63',HOT_PrMin_G)
         ALLOCATE(HOT_RSNMAX_G(NP_G))
         HOT_RSNMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxrs.63',HOT_RSNMAX_G)
         IF(LoadEleSlopeLim)THEN
            HOT_ESLONOFF_G(:) = 0.d0 ! initialize the value
            CALL readMinMax('ESLNodes.63',HOT_ESLONOFF_G)
         ENDIF   
      endif
      if ( (MNPROC.eq.1).or.(WRITE_LOCAL_FILES) ) then
         ALLOCATE(HOT_ETAMAX(NP))
         HOT_ETAMAX(:) = -99999.d0    ! initialize the value
         CALL readMinMax('maxele.63',HOT_ETAMAX)
         ALLOCATE(HOT_UMAX(NP))
         HOT_UMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxvel.63',HOT_UMAX)
         ALLOCATE(HOT_WVNOUTMAX(NP))
         HOT_WVNOUTMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxwvel.63',HOT_WVNOUTMAX)
         ALLOCATE(HOT_PrMin(NP))
         HOT_PrMin(:) = 99999.d0  ! initialize the value
         CALL readMinMax('minpr.63',HOT_PrMin)
         ALLOCATE(HOT_RSNMAX(NP))
         HOT_RSNMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxrs.63',HOT_RSNMAX)
         IF(LoadEleSlopeLim)THEN
            HOT_ESLONOFF(:) = 0.d0 ! initialize the value
            CALL readMinMax('ESLNodes.63',HOT_ESLONOFF)
         ENDIF   
      endif

C...SEND THE VALUES FROM THE PREVIOUS ELEMENTAL SLOPE LIMITER RUN
C   TO THE OTHER PROCESSORS AND PROCESS TO THE SUBDOMAINS
      IF(LoadEleSlopeLim)THEN
        IF(MNPROC.GT.1)THEN
#ifdef CMPI
            CALL MPI_BCAST( HOT_ESLONOFF_G, NP_G, REALTYPE, 0, COMM, 
     &          IERR )
#endif
            DO sd_node_number = 1,NP
                fd_node_number = ABS(NODES_LG(sd_node_number))
                HOT_ESLONOFF(sd_node_number) = 
     &              HOT_ESLONOFF_G(fd_node_number)
                IF(HOT_ESLONOFF(sd_node_number).EQ.1D0)THEN
                    elemental_slope_limiter_active(sd_node_number) = 
     &                  .TRUE.
                ELSE
                    elemental_slope_limiter_active(sd_node_number) = 
     &                  .FALSE.
                ENDIF  
            ENDDO
        ELSEIF(MNPROC.EQ.0)THEN
            DO I = 1,NP
                IF(HOT_ESLONOFF(sd_node_number).EQ.1D0)THEN
                    elemental_slope_limiter_active(sd_node_number) = 
     &                  .TRUE.
                ELSE
                    elemental_slope_limiter_active(sd_node_number) =
     &                  .FALSE.               
                ENDIF  
            ENDDO    
        ENDIF
      ENDIF  

C
 1112 FORMAT(/,1X,79('_'))
 1197 FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &                'BEING COMPUTED ')
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END SUBROUTINE HOTSTART

C******************************************************************************
C   Subroutine to initialize the 3D routines for a hot start including        *
C   finding proper places in 3D output files.                                 *
C                                                                             *
C   jgf49.43: Had to move the reading of the 3D portion of the hotstart file
C   to the hotstart() subroutine, to facilitate separation of the reading of
C   the file from initialization. This is required if the file is to be read
C   in various formats (e.g., netcdf).
C                                                                             *
c******************************************************************************

      SUBROUTINE HOTSTART_3D(TimeLoc,ITHS)
      USE GLOBAL_3DVS

C   kmd48.33bc - variable information needed from the global.F file
! kmd - added variables for river information
      USE GLOBAL, ONLY : RES_BC_FLAG, RBCTIME1, RBCTIME2,
     &                   RBCTIMEINC, SBCTIME1, SBCTIME2,
     &                   SBCTIMEINC, TBCTIME1, TBCTIME2,
     &                   TBCTIMEINC, BCSTATIM, RBCRATIO,
     &                   LNM_BC1, LNM_BC2, TIMEIT, STATIM,
     &                   SBCSTATIM, TBCSTATIM, SBCRATIO,
     &                   TBCRATIO, DTDPHS, LNM_BC,
     &                   BCFLAG_LNM, BCFLAG_TEMP, TTBCTIME1,
     &                   TTBCTIME2, TTBCSTATIM, TTBCTIMEINC,
     &                   q_heat1, q_heat2, HFLUX, TTBCRATIO,
     &                   setMessageSource, unsetMessageSource,
     &                   DEBUG, ECHO, INFO, WARNING, ERROR, 
     &                   NP_G, allMessage, screenMessage,
     &                   NBD, RIVBCTIMINC, RIVBCTIME1,RIVBCTIME2,
     &                   totalbcrivernodes, BCRivRATIO, NOPE,
     &                   RIVBCSTATIM

      IMPLICIT NONE

      INTEGER :: IRType                  !number of fields in output file
      INTEGER :: NH, N                   !horizontal & vertical loop counters
!      REAL(SZ) :: RealPartOfQ
C  kmd48.33bc - new variables
      INTEGER :: NumofBCNodes, J, NOD, IT
      INTEGER, INTENT(IN) :: ITHS
      REAL(SZ), INTENT(IN) :: TimeLoc

      CHARACTER(80) :: CDUM80

      INTEGER :: IINDX
      INTEGER :: NumofNodes, K
      REAL(SZ),ALLOCATABLE :: TMP(:,:)

      REAL(SZ) :: tempp, temp1, temp2

      call setMessageSource("hotstart_3D")
#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
C...  Define format statements used to initialize 3D output files
C...
 499  FORMAT(1X,A32,2X,A24,2X,A24)
C.. RJW bug fix in format 498 (kendra found this)
 498  FORMAT(1X,I10,1X,I10,1X,E15.7,I10,1X,I10,1X,I3)
 497  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE ASCII')
 496  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE BINARY')

C  kmd48.33bc - restart boundary condtions for baroclinic simulations
       IF((RES_BC_FLAG.GT.0).AND.(CBAROCLINIC)) THEN
         IF ((ABS(RES_BC_FLAG).GE.1).AND.(NOPE.GT.0)) THEN
            IF (BCFLAG_LNM.EQ.1) THEN
              OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
              RBCTIME1=BCSTATIM*86400.d0
              RBCTIME2=RBCTIME1+RBCTIMEINC
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC1(NumofBCNodes)
              END DO
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC2(NumofBCNodes)
              END DO
              DO IT=1,ITHS
                 TIMEIT=IT*DTDPHS + STATIM*86400.D0
                 IF(TIMEIT.GT.RBCTIME2) THEN
                    RBCTIME1=RBCTIME2
                    RBCTIME2=RBCTIME1+RBCTIMEINC
                    READ(35,'(A)') CDUM80
                    DO NumofBCNodes=1,NETA
                       LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                       READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                    END DO
                 END IF
              END DO
              IF(TimeLoc.GT.RBCTIME2) THEN
                 RBCTIME1=RBCTIME2
                 RBCTIME2=RBCTIME2+RBCTIMEINC
                 READ(35,'(A)') CDUM80
                 DO NumofBCNodes=1,NETA
                   LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                   READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                 END DO
              END IF
              RBCRATIO=(TIMEIT-RBCTIME1)/RBCTIMEINC
            ELSE IF (BCFLAG_LNM.EQ.2) THEN
              ! to be added later but currently information
              ! is taken in two ways from the HYCOM results
            ELSE IF (BCFLAG_LNM.EQ.3) THEN
              DO NumofBCNodes=1,NETA
                 IINDX=NBD(NumofBCNodes)
                 LNM_BC1(NumofBCNodes)=ETA2(IINDX)
                 LNM_BC2(NumofBCNodes)=LNM_BC1(NumofBCNodes)
              END DO
            END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.2).OR.(ABS(RES_BC_FLAG).EQ.4)) THEN
            IF (NOPE.GT.0) THEN
              OPEN(36,FILE=TRIM(INPUTDIR)//'/'//'fort.36')
              SBCTIME1=SBCSTATIM*86400.D0
              SBCTIME2=SBCTIME1+SBCTIMEINC
              READ(36,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(36,*) NOD,(RESSAL1(NumofBCNodes,J),J=1,NFEN)
              END DO
              READ(36,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
              END DO
              DO IT=1,ITHS
                 TIMEIT=IT*DTDPHS + STATIM*86400.D0
                 IF(TIMEIT.GT.SBCTIME2) THEN
                    SBCTIME1=SBCTIME2
                    SBCTIME2=SBCTIME1+SBCTIMEINC
                    READ(36,'(A)') CDUM80
                    DO NumofBCNodes=1,NETA
                       DO J=1,NFEN
                         RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                       END DO
                       READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
                    END DO
                 END IF
              END DO
              IF(TIMELOC.GT.SBCTIME2) THEN
                 SBCTIME1=SBCTIME2
                 SBCTIME2=SBCTIME2+SBCTIMEINC
                 READ(36,'(A)') CDUM80
                 DO NumofBCNodes=1,NETA
                   DO J=1,NFEN
                     RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                   END DO
                   READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
                 END DO
              END IF
              SBCRATIO=(TIMEIT-SBCTIME1)/SBCTIMEINC
            END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &            (ABS(RES_BC_FLAG).EQ.4)) THEN
            IF (NOPE.GT.0) THEN
              OPEN(37,FILE=TRIM(INPUTDIR)//'/'//'fort.37')
              TBCTIME1=TBCSTATIM*86400.D0
              TBCTIME2=TBCTIME1+TBCTIMEINC
              READ(37,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(37,*) NOD,(RESTEMP1(NumofBCNodes,J),J=1,NFEN)
              END DO
              READ(37,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
              END DO
              DO IT=1,ITHS
                 TIMEIT=IT*DTDPHS + STATIM*86400.D0
                 IF(TIMEIT.GT.TBCTIME2) THEN
                    TBCTIME1=TBCTIME2
                    TBCTIME2=TBCTIME1+TBCTIMEINC
                    READ(37,'(A)') CDUM80
                    DO NumofBCNodes=1,NETA
                       DO J=1,NFEN
                         RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                       END DO
                       READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
                    END DO
                 END IF
              END DO
              IF(TIMELOC.GT.TBCTIME2) THEN
                 TBCTIME1=TBCTIME2
                 TBCTIME2=TBCTIME2+TBCTIMEINC
                 READ(37,'(A)') CDUM80
                 DO NumofBCNodes=1,NETA
                   DO J=1,NFEN
                     RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                   END DO
                   READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
                 END DO
              END IF
              TBCRATIO=(TIMEIT-TBCTIME1)/TBCTIMEINC
            END IF
            IF (BCFLAG_TEMP.EQ.1) THEN ! read in file with multiple values but no calcs
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               DO NumofNodes=1,NP
                  READ(38,*) NOD, q_heat1(Numofnodes)
               END DO
               DO NumofNodes=1,NP
                  READ(38,*) NOD, q_heat2(NumofNodes)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NH
                       q_heat1(NumofNodes)=q_heat2(NumofNodes)
                       READ(38,*) NOD, q_heat2(NumofNodes)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
            ELSE IF (BCFLAG_TEMP.EQ.2) THEN ! read in from file with 6 components
               ALLOCATE (TMP(NP,6))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &                     TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &                     TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &                     TMP(NumofNodes,6)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            ELSE IF (BCFLAG_TEMP.EQ.3) THEN ! read in from file with 4 components
               ALLOCATE (TMP(NP,4))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=TMP(NumofNodes,4)+
     &               TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &               TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &                     TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &                     TMP(NumofNodes,2)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            END IF ! end if for heat flux file
         END IF ! end if for temperature boundary and heat flux
      ELSE IF ((RES_BC_FLAG.LT.0).AND.(CBAROCLINIC).AND.(NOPE.GT.0)) THEN
         IF (BCFLAG_LNM.EQ.1) THEN
            OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
            READ(35,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC(NumofBCNodes)
            END DO
         ELSE IF (BCFLAG_LNM.EQ.2) THEN
            ! add later
         ELSE IF (BCFLAG_LNM.EQ.3) THEN
            DO NumofBCNodes=1,NETA
               IINDX=NBD(NumofBCNodes)
               LNM_BC(NumofBCNodes)=ETA2(IINDX)
             END DO
         END IF
      END IF

C kmd - River information for boundary conditions
      IF (BndBCRiver) THEN
        OPEN(39,FILE=TRIM(INPUTDIR)//'/'//'fort.39')
        READ(39,*) RIVBCTIMINC, RIVBCSTATIM
        RIVBCTIME1=RIVBCSTATIM*86400.d0
        RIVBCTIME2=RIVBCTIME1+RIVBCTIMINC
        DO J=1,totalbcrivernodes
         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN1(J,K), BCRivTempN1(J,K),K=1,NFEN)
         END IF
        END DO
        DO J=1,totalbcrivernodes
         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN2(J,K), BCRivTempN2(J,K),K=1,NFEN)
         END IF
        END DO
        DO IT=1,ITHS-1
           TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
           IF(TIMEIT.GT.RIVBCTIME2) THEN
              RIVBCTIME1=RIVBCTIME2
              RIVBCTIME2=RIVBCTIME2+RIVBCTIMINC
              DO J=1,totalbcrivernodes
                IF (IDEN.EQ.2) THEN
                 DO K=1,NFEN
                   BCRivSalN1(J,K)=BCRivSalN2(J,K)
                 END DO
                 READ(39,*) (BCRivSalN2(J,K),K=1,NFEN)
                ELSE IF (IDEN.EQ.3) THEN
                 DO K=1,NFEN
                   BCRivTempN1(J,K)=BCRivTempN2(J,K)
                 END DO
                 READ(39,*) (BCRivTempN2(J,K),K=1,NFEN)
                ELSE IF (IDEN.EQ.4) THEN
                 DO K=1,NFEN
                   BCRivSalN1(J,K)=BCRivSalN2(J,K)
                   BCRivTempN1(J,K)=BCRivTempN2(J,K)
                 END DO
                 READ(39,*) (BCRivSalN2(J,K), BCRivTempN2(J,K),K=1,NFEN)
                END IF
              END DO
           END IF
        END DO
        BCRivRATIO=(TIMEIT-RIVBCTIME1)/RIVBCTIMINC
        DO J=1,totalbcrivernodes
           IF (IDEN.EQ.2) THEN
             DO K=1,NFEN
               BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO*
     &             (BCRivSalN2(J,K)-BCRivSalN1(J,K))
             END DO
           ELSE IF (IDEN.EQ.3) THEN
             DO K=1,NFEN
               BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO*
     &             (BCRivTempN2(J,K)-BCRivTempN1(J,K))
             END DO
           ELSE IF (IDEN.EQ.4) THEN
             DO K=1,NFEN
               BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO*
     &             (BCRivSalN2(J,K)-BCRivSalN1(J,K))
               BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO*
     &             (BCRivTempN2(J,K)-BCRivTempN1(J,K))
             END DO
           END IF
        END DO
      END IF

C.... Initialize the station 3D density output file
!kmd48.33bc - added in the information in front of all
!             the other output files for resetting the
!             information.
      IF (I3DSD.LT.0) THEN
         CALL setOutputTiming(N3DSD, I3DSDRec,
     &      NTO3DSDS, NTO3DSDF, NSpo3DSD, NDSet3DSD, iths, dtdphs)
      END IF


      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3

      CALL initOutput3D(41, I3DSD, NDSET3DSD, NSTA3DD, NSTA3DD_g,
     &                        NSPO3DSD, IRType, I3DSDRec)

C.... Initialize the 3D velocity station output file (Unit 42)
!kmd48.33bc - added in the information in front of all
!             the other output files for resetting the
!             information.
      IF (I3DSV.LT.0) THEN
         CALL setOutputTiming(N3DSV, I3DSVRec,
     &      NTO3DSVS, NTO3DSVF, NSpo3DSV, NDSet3DSV, iths, dtdphs)
      END IF

      IRType=3

      CALL initOutput3D(42, I3DSV, NDSET3DSV, NSTA3DV, NSTA3DV_g,
     &                        NSPO3DSV, IRType, I3DSVRec)


C.... Initialize the 3D turbulence station output file (Unit 43)
!kmd48.33bc - added in the information in front of all
!             the other output files for resetting the
!             information.
      IF (I3DST.LT.0) THEN
         CALL setOutputTiming(N3DST, I3DSTRec,
     &      NTO3DSTS, NTO3DSTF, NSpo3DST, NDSet3DST, iths, dtdphs)
      END IF

      IRType=3

C  kmd48.33bc - changed NDSet3DST to use DTDPHS
      IF (DTDP.NE.DTDPHS) THEN
         NDSet3DST=((NTO3DSTF*DTDP)-(NTO3DSTS*DTDPHS))/
     &              (NSpo3DST*DTDP)
      END IF

      CALL initOutput3D(43, I3DST, NDSET3DST, NSTA3DT, NSTA3DT_g,
     &                        NSPO3DST, IRType, I3DSTRec)

C.... Initialize the global 3D density output file (Unit 44)

C  kmd48.33bc - change NDSet3DGD to use DTDPHS
      IF(I3DGD.LT.0) THEN
         CALL setOutputTiming(N3DGD, I3DGDRec,
     &      NTO3DGDS, NTO3DGDF, NSpo3DGD, NDSet3DGD, iths, dtdphs)
      END IF

      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
C     jgf46.27 Replaced IRType with IDen
!kmd48.33bc - added in the information to make a file with all the
!             top temperature boundary condition information in

      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_g,
     &                        NSPO3DGD, IRType, I3DGDRec)
      IF(I3DGD.EQ.-1) THEN         !start a new ASCII file
        IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
          IF (BCFLAG_TEMP.NE.0) THEN
            OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47')
            WRITE(47,499) RUNDES,RUNID,AGRID
            WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
            CLOSE(47)
          END IF
        END IF
      ENDIF

c.... Initialize the global 3D velocity output file (Unit 45)

C  kmd48.33bc - changed NDSet3DGV to use DTDPHS and added REALs to
C              calculation of the start time due to timestep change
      IF(I3DGV.LT.0) THEN
         CALL setOutputTiming(N3DGV, I3DGVRec,
     &      NTO3DGVS, NTO3DGVF, NSpo3DGV, NDSet3DGV, iths, dtdphs)
      ENDIF

      CALL initOutput3D(45, I3DGV, NDSET3DGV, NP, NP_g,
     &                        NSPO3DGV, IRType, I3DGVRec)

c.... Initialize the global 3D turbulence output file (Unit 46)

C  kmd48.33bc - changed NDSet3DGT to use DTDPHS
      IF(I3DGT.LT.0) THEN
         CALL setOutputTiming(N3DGT, I3DGTRec,
     &      NTO3DGTS, NTO3DGTF, NSpo3DGT, NDSet3DGT, iths, dtdphs)
      ENDIF

      CALL initOutput3D(46, I3DGT, NDSET3DGT, NP, NP_g,
     &                        NSPO3DGT, IRType, I3DGTRec)

C.... Set up a few final odds and ends for a 3D run

      CALL VSSTUP ()

#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE HOTSTART_3D
C----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E T   O U T P U T   T I M I N G
C-----------------------------------------------------------------------
C     jgf49.48.04 For output files that start anew on hotstart, calculate
C     the timing of outputs, taking into account that the run that wrote
C     the hotstart file may have had a different time step than this run.
C-----------------------------------------------------------------------
      SUBROUTINE setOutputTiming(tsCounter, recCounter,
     &   startTS, endTS, tsPeriod, nSets, iths, dtdphs)
      USE GLOBAL, ONLY : DTDP, setMessageSource, unsetMessageSource,
     &   screenMessage, DEBUG, allMessage, CHOTHS
      IMPLICIT NONE
      INTEGER, intent(out) :: tsCounter  ! time step counter btw outputs
      INTEGER, intent(out) :: recCounter ! counts lines in the output file
      INTEGER, intent(inout) :: startTS  ! time step for output to start
      INTEGER, intent(inout) :: endTS    ! time step for output to end
      INTEGER, intent(inout) :: tsPeriod ! period of time steps between outputs
      INTEGER, intent(out) :: nSets      ! num data sets in output file
      INTEGER, intent(in) :: iths
      REAL(8), intent(in) :: dtdphs
C
      call setMessageSource("setOutputTiming")
#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      tsCounter=0
      recCounter=0
      IF ((startTS.LT.ITHS).AND.(tsPeriod.GT.0)) THEN
         startTS=startTS+((REAL(ITHS)-REAL(startTS))/
     &           REAL(tsPeriod))*REAL(tsPeriod)
         IF (startTS.LT.ITHS) startTS=startTS+tsPeriod
         IF (CHOTHS.eqv..true.) THEN
            nSets=((endTS*DTDP)-(startTS*DTDPHS))/
     &                (tsPeriod*DTDP)
         ELSE
           nSets=(endTS-startTS)/tsPeriod
         END IF
      END IF
#if defined(HOTSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C----------------------------------------------------------------------
      END SUBROUTINE setOutputTiming
C----------------------------------------------------------------------
