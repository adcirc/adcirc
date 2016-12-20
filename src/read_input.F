C******************************************************************************
C PADCIRC VERSION 46.00 xx/xx/2006                                            *
C  last changes in this file VERSION 46.00                                    *
C                                                                             *
C This module handles most of the model input.  The primary 2d input is read  *
C in subroutine READ_INPUT.  The primary 3D input is read in subroutine       *
C READ_INPUT_3D.  Initial conditions that are read in for a cold start are    *
C handled in the cold start subroutines.                                      *
C                                                                             *
C******************************************************************************

C-----------------------------------------------------------------------
C     S U B R O U T I N E     R E A D _ I N P U T
C-----------------------------------------------------------------------
C
C     READS INPUT FILES
C
C-----------------------------------------------------------------------
      SUBROUTINE READ_INPUT()
      USE SIZES               !since GLOBAL uses SIZES, is this necessary?
      USE GLOBAL
      USE MESH, ONLY : NE, SLAM0, SFEA0, BCXY, RMAX, ICS, X, Y, DP,
     &   SLAM, SFEA, NM, readMesh, initializeMesh, initializeBoundaries, 
     &   freeMesh
      USE BOUNDARIES, ONLY : ANGINN, COSTSET, NETA, NVEL, NVELME, 
     &   NFLUXF, NBD, ME2GW, LBCODEI, NBV
      USE GLOBAL_IO
      USE KDTREE2_MODULE     ! V49.48.02 tcm -- added for fast searching
      USE GLOBAL_3DVS, ONLY: C3D_PTrans, NSTA3DD, NSTA3DV, NSTA3DT,
     &    I3DSD,I3DSV,I3DST,I3DGD,I3DGV,I3DGT, NSTA3DD_G, NSTA3DV_G,
     &    NSTA3DT_G
      USE SUBDOMAIN, ONLY : subdomainOn
      USE WRITE_OUTPUT, ONLY : outputNodeCode, outputNOFF
      USE HARM
      USE WEIR
      USE TIME_VARYING_WEIR_BOUNDARY
      USE WIND
      USE ITPACKV
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate    
      USE VERSION
!     jgf46.00
      USE NodalAttributes, ONLY :
     &     NoLiBF, NWP, Tau0, HBreak, FTheta, FGamma, Tau, CF,
     &     InitNAModule, ReadNodalAttr, InitNodalAttr, ESLM, ESLC,
     &     Tau0FullDomainMin, Tau0FullDomainMax, Z0b_var
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
#ifdef CSWAN
Casey 100205: Enable hot-start file generation by SWAN.
      USE Couple2Swan, ONLY: SwanHotStartUnit
#endif
      IMPLICIT NONE
      INTEGER NIBP,IBN1,IK,NDISC,NBBN,NVEL2
      INTEGER NTimeVaryingWeir
      INTEGER IDUM80
      INTEGER IMDig1,IMDig2,IMDig3,IMDig4,IMDig5,IMDig6
      REAL(SZ) RampVal
      CHARACTER(len=80) :: CDUM80
      LOGICAL fileFound
      INTEGER :: I, J, JKI
C
#if defined CSWAN || defined ADCSWAN
      ! jgf50.60.09: Namelist for turning SWAN output files on or off
      NAMELIST /SWANOutputControl/ SWAN_OutputHS, SWAN_OutputDIR,
     &    SWAN_OutputTM01, SWAN_OutputTPS, SWAN_OutputWIND,
     &    SWAN_OutputTM02, SWAN_OutputTMM10
#endif
      CHARACTER(len=1000) namelistSpecifier ! used in log messages
      NAMELIST /metControl/ WindDragLimit,DragLawString,rhoAir,
     &  invertedBarometerOnElevationBoundary
      NAMELIST /subdomainModeling/ subdomainOn
      NAMELIST /wetDryControl/ outputNodeCode, outputNOFF, noffActive,
     &   outputUU1VV1, outputEta1, outputEta2
      NAMELIST /inundationOutputControl/ inundationOutput, inunThresh
      NAMELIST /TVWControl/ USE_TVW,TVW_FILE,NOUT_TVW,TOUTS_TVW,
     &                      TOUTF_TVW,NSPOOL_TVW
      
      CHARACTER(len=160) :: origLine     ! raw line of input from fort.15
      CHARACTER(len=200) :: modifiedLine ! converted to namelist
      CHARACTER(len=2000) :: origLineTVW
      INTEGER :: commentStart            ! location of fort.15 comment start
      INTEGER :: ios = 0                 ! i/o status of read operation

      !jgf51.52.35: for log messages about IDEN
      character(len=100) :: densityRunType ! diagnostic or prognostic
      character(len=100) :: densityDimensions ! 2DDI or 3D
      character(len=100) :: densityForcingType ! sigmaT, salinity, etc

      !tcm v51.20.04 Additions for External Station Location File
      INTEGER :: IOS_STATIONS = 0
      INTEGER :: STAT_LUN = 15  !station unit number fort.15 file be default
      INTEGER :: NSTAE2,NSTAV2,NSTAC2,NSTAM2      
      LOGICAL :: USE_ELEV_STAT_FILE = .false.    !.true. only if an elevation station file exists
      LOGICAL :: USE_VEL_STAT_FILE = .false.     !.true. only if a velocity station file exists
      LOGICAL :: USE_CONC_STAT_FILE = .false.    !.true. only if a concentration station file exists
      LOGICAL :: USE_MET_STAT_FILE = .false.     !.true. only  if a met station file exists

!   tcm v50.66.02 -- added timebathycontrol namelist related variables
      INTEGER :: ios_nddt
      logical :: found_tbc_nml = .false.   !flag to determine if the timebathycontrol namelist was present
      NAMELIST /timebathycontrol/ NDDT,BTIMINC,BCHGTIMINC

Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      NAMELIST /waveCoupling/ WaveWindMultiplier

      call setMessageSource("read_input")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
C...PRINT OUT HEADER FOR OUTPUT INCLUDING VERSION NUMBER AND COPYRIGHT
C...
      WRITE(16,1112)
      WRITE(16,1112)
      WRITE(16,1114) TRIM(ADC_VERSION)
      WRITE(16,1112)

 1114 FORMAT(//,19X,'PROGRAM ADCIRC   VERSION ',A,
     &     //,5X,'AN ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL ',
     &     'SEAS AND ESTUARIES',
     &     ///,7X,'-  DEVELOPED BY',
     &     //,10X,'R.A. LUETTICH, JR',
     &     /,12X,'UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL',
     &     /,12X,'INSTITUTE OF MARINE SCIENCES',
     &     //,10X,'J.J. WESTERINK ',
     &     /,12X,
     &     'DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES',
     &     /,12X,'UNIVERSITY OF NOTRE DAME',
     &     ///,7X,'-  THE ADCIRC SOURCE CODE IS COPYRIGHTED BY',
     &     //,10X,'R.A. LUETTICH, JR. AND J.J. WESTERINK, 1994-2006',
     &     //,7X,
     &     'NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED',
     &     /,10X,'WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS',//)
C...
C...  WRITE OUT HEADER INFORMATION DESCRIBING HOW THE CODE HAS BE SET UP
C...
      WRITE(16,1210)
 1210 FORMAT(//,1X,'THE ADCIRC SOURCE CODE HAS BEEN CONFIGURED ',
     &             'BY THE PREPROCESSOR AS FOLLOWS:',/)

#ifdef CMACHSUN
      WRITE(16,*) '      - CODE SETUP TO RUN ON SUN 4 OR SPARC ',
     &            'COMPUTERS'
#endif

#ifdef REAL4
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 4 byte REALS'
      WRITE(16,4312)
4312  FORMAT(/,'**** W A R N I N G ****',/
            ' You are running ADCIRC in SINGLE PRECISION! ',/
            ' It is always recommended that you run ADCIRC in ',/
            ' DOUBLE PRECISION ONLY.'
#else
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 8 byte REALS'
#endif

#ifdef CVEC
      WRITE(16,*) '      - CODE OPTIMIZED FOR A VECTOR COMPUTER'
#endif

#ifdef CSCA
      WRITE(16,*) '      - CODE OPTIMIZED FOR A SCALAR COMPUTER'
#endif
      write(16,1112) 
C
C     jgf46.00 Zero out all the variables in the Nodal Attributes
C     Module.
      CALL initNAModule()
C
C     jgf50.32: Initialize default values in wind module.
      CALL initWindModule()

c     tcm v50.66.01 addition for time varying bathymetry
      NDDT = 0  !set the default value for time varying bathymetry e.g. no changes
      ios_nddt = 0
C
      fileFound = .false.
C
C     jgf51.12.11: moved the reading of the mesh and boundaries file
C     (called fort.14 by default) to the mesh module. Also added the
C     capability to read different formats and use nondefault filenames
C     via command line options in sizes.F.
      call readMesh()
C
C     ALLOCATE ARRAYS Dimensioned by MNP and MNE
      call alloc_main1a()
      ! jgf51.21.13: Allocate boundary condition arrays dimensioned by
      ! mnope and mneta
      call alloc_main2()
      ! jgf51.21.13: Allocate boundary condition arrays dimensioned by
      ! mnvel 
      call alloc_main3()
C...
C...  OPEN UNIT 15 INPUT FILE (control parameters and periodic boundary conditions)
C...
      CALL openFileForRead(15,TRIM(INPUTDIR)//'/'//'fort.15',ios)
      IF (ios.gt.0) THEN
         CALL ADCIRC_Terminate()
      ENDIF
C
C...GENERAL PURPOSE FORMAT STATEMENTS for subtly expressed error messages
C...
 1112 FORMAT(/,1X,79('_'))
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
 9974 FORMAT(/,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!!!',//)


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
      namelistSpecifier = 'TimeBatyhControl'
      READ(UNIT = 15,NML = TimeBathyControl,IOSTAT = IOS_NDDT) 
      call logNamelistReadStatus(namelistSpecifier, ios_nddt)     
c     it is possible for the namelist to be present in the file and 
C     occuring at the end of the file with no line breaks after the 
C     ending "\" which causes the iostat to return a negative value.
C     By checking to be sure a namelist variable was set to
c     a non-default value we can determine this was the case.
      select case(ios_nddt)
      case(:-1) 
         if (nddt.ne.0) then
            call allMessage(WARNING,
     &        'TimeBathyControl NAMELIST PRESENT, BUT AT THE '
     &     // 'END OF FILE WITH NO ADVANCING CHARACTER.')
            found_tbc_nml = .true.
         endif
      case(0)
         found_tbc_nml = .true.
      case(1:)
         found_tbc_nml = .true.
         call allMessage(ERROR,
     &   'THERE WAS A PROBLEM PROCESSING THE TimeBathyControl NAMELIST '
     &    //'IN THE FORT.15 FILE.  SHUTTING DOWN ADCIRC NOW.')
         call adcirc_terminate()
      end select     
      write(scratchMessage,'(a,l,a,i0,a,e15.8,a,e15.8)') 
     &  'found_tbc_nml=',found_tbc_nml,
     &  ' nddt=',nddt,' btiminc=',btiminc,' bchgtiminc=',bchgtiminc
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)

#if defined CSWAN || defined ADCSWAN
      ! jgf50.60.08: Add a namelist for the user to turn SWAN output
      ! on and off. Similar to tcm's timevaryingbathy namelist.
      namelistSpecifier = 'SWANOutputControl'
      read(unit=15,nml=SWANOutputControl,iostat=ios)
      call logNamelistReadStatus(namelistSpecifier, ios)     
      call logMessage(ECHO,
     &    "The values of SWANOutputControl are as follows:")
      write(scratchMessage,*) "SWAN_OutputHS=",SWAN_OutputHS,
     &    " SWAN_OutputDIR=",SWAN_OutputDIR,
     &    " SWAN_OutputTM01=",SWAN_OutputTM01, 
     &    " SWAN_OutputTPS=",SWAN_OutputTPS, 
     &    " SWAN_OutputWIND=",SWAN_OutputWIND,
     &    " SWAN_OutputTM02=",SWAN_OutputTM02, 
     &    " SWAN_OutputTMM10=",SWAN_OutputTMM10
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)
#endif
      !
      ! jgf51.42: Add a namelist for the user to control subdomain
      ! modeling. 
      namelistSpecifier = 'subdomainModeling'
      read(unit=15,nml=subdomainModeling,iostat=ios)
      call logNamelistReadStatus(namelistSpecifier, ios)
      write(scratchMessage,*) "subdomainOn=",subdomainOn
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)
      !
      ! jgf50.60.13: Add a namelist for the user to control met forcing.
      ! Similar to tcm's timevaryingbathy namelist.
      namelistSpecifier = 'metControl'
      read(unit=15,nml=metControl,iostat=ios)
      call logNamelistReadStatus(namelistSpecifier, ios)
      write(scratchMessage,*) "WindDragLimit=",WindDragLimit,
     & " DragLawString=",DragLawString, " rhoAir=",rhoAir,
     & ' invertedBarometerOnElevationBoundary=', 
     &   invertedBarometerOnElevationBoundary 
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)

      !...Read in time varying weir control from fort.15
      NOUT_TVW = -99999
      TOUTS_TVW = -99999
      TOUTF_TVW = -99999
      NSPOOL_TVW = -99999
      READ(UNIT=15,NML=TVWControl,IOSTAT=IOS)
      IF(IOS.NE.0)THEN
        CALL logMessage(INFO,"The tvwControl namelist was not found.")
        found_tvw_nml = .FALSE.
      ELSE
        CALL logMessage(INFO,"The tvwControl namelist was found.")
        found_tvw_nml = .TRUE.
        IF(USE_TVW)THEN
            IF(NOUT_TVW.EQ.-99999D0)THEN
                NOUT_TVW=NOUTGE
                IF(TOUTS_TVW.EQ.-99999D0)TOUTS_TVW=TOUTSGE
                IF(TOUTF_TVW.EQ.-99999D0)TOUTF_TVW=TOUTFGE
                IF(NSPOOL_TVW.EQ.-99999D0)NSPOOL_TVW=NSPOOLGE
            ELSE
                IF(TOUTS_TVW.EQ.-99999D0)TOUTS_TVW=TOUTSGE
                IF(TOUTF_TVW.EQ.-99999D0)TOUTF_TVW=TOUTFGE
                IF(NSPOOL_TVW.EQ.-99999D0)NSPOOL_TVW=NSPOOLGE
            ENDIF    
        ELSE
            NOUT_TVW = 0
        ENDIF   
      ENDIF
      REWIND(15)

Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      namelistSpecifier = 'waveCoupling'
      read(unit=15,nml=waveCoupling,iostat=ios)
      call logNamelistReadStatus(namelistSpecifier,ios)
      write(scratchMessage,*) "WaveWindMultiplier=",WaveWindMultiplier
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)
      !
      ! jgf52.05: Add a namelist for the user to control
      ! output of NODECODE and NOFF.
      namelistSpecifier = 'wetDryControl'
      read(unit=15,nml=wetDryControl,iostat=IOS)
      call logNameListReadStatus(namelistSpecifier,ios)
      write(scratchMessage,'(a,l)') "outputNodeCode=",outputNodeCode
      call logMessage(ECHO,trim(scratchMessage))
      write(scratchMessage,'(a,l)') "outputNOFF=",outputNOFF
      call logMessage(ECHO,trim(scratchMessage))
      write(scratchMessage,'(a,l)') "noffActive=",noffActive
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)

      ! jgf52.08.01: Add a namelist for the user to control
      ! output of detailed inundation data. 
      namelistSpecifier = 'inundationOutputControl'
      READ(UNIT=15,NML=inundationOutputControl,IOSTAT=IOS)
      call logNamelistReadStatus(namelistSpecifier,ios)
      write(scratchMessage,'(a,l)') "inundationOutput=",inundationOutput
      call logMessage(ECHO,trim(scratchMessage))
      write(scratchMessage,'(a,e15.8)') "inunThresh=",inunThresh
      call logMessage(ECHO,trim(scratchMessage))
      rewind(15)
C...
C...  INPUT FROM UNIT 15 AND OUTPUT TO UNIT 16 RUN DESCRIPTION AND RUN
C...  IDENTIFICATION
C...
      READ(15,'(A80)') RUNDES
      READ(15,'(A80)') RUNID
      WRITE(16,1) RUNDES
  1   FORMAT(//,1X,'RUN DESCRIPTION : ',A80)
      WRITE(16,209) RUNID
209   FORMAT(/,1X,'RUN IDENTIFICATION : ',A80)

      do I=1,20
         J=(I-1)*4+1
         RDES4(I)=RUNDES(J:J+3)
         RID4(I) =RUNID (J:J+3)
      end do
      do I=1,10
         J=(I-1)*8+1
         RDES8(I)=RUNDES(J:J+7)
         RID8(I) =RUNID (J:J+7)
      end do

C...
C... READ AND PROCESS NFOVER - NONFATAL ERROR OVERRIDE OPTION
C...
C
C     jgf46.10 Add user-controllable warning, output, and stop criteria
C     for elevations. Initialize default values.
      WarnElev = 20.0         ! default
      iWarnElevDump = 0       ! init
      WarnElevDump = .False.  ! default
      WarnElevDumpLimit = 50  ! default
      WarnElevDumpCounter = 0 ! init
      ErrorElev = 1000.0      ! default
C
#ifndef DEBUG_WARN_ELEV
      READ(15,*) NFOVER
#else
      READ(15,*) NFOVER, WarnElev, iWarnElevDump, WarnElevDumpLimit,
     &  ErrorElev
#endif
C
      WRITE(16,1112)
      WRITE(16,1250)
 1250 FORMAT(//,1X,'GENERAL RUN INFORMATION',/)
      IF(NFOVER.EQ.1) THEN
         WRITE(16,1951) NFOVER
 1951    FORMAT(5X,'NFOVER = ',I2,
     &        /,9X,'IF NON-FATAL INPUT ERRORS ARE DETECTED, THEY WILL ',
     &        'BE CORRECTED AND EXECUTION CONTINUED')
      ELSE
         WRITE(16,1952) NFOVER
 1952    FORMAT(/,5X,'NFOVER = ',I3,
     &        /,9X,'NON-FATAL INPUT ERRORS WILL STOP EXECUTION ',/)
      ENDIF
#ifdef DEBUG_WARN_ELEV
      IF (iWarnElevDump .ne. 0) WarnElevDump = .True.
      WRITE(16,1953) WarnElev,WarnElevDump,WarnElevDumpLimit,ErrorElev
 1953 FORMAT(//,5X,
     &     'A warning will be issued if elevation exceeds WarnElev = ',
     &     e16.8,
     &     /,5X,'A global elevation file (fort.69) will be written if '
     &     /,5X,'WarnElev is exceeded and WarnElevDump is true: ',L2,
     &     /,5X,'Execution will be terminated if ',
     &     '(WarnElevDumpLimit = 'I3,') ',
     &     /,5X,'global elevation files have been written as warning.'
     &     /,5X,'Execution will be terminated if elevation exceeds'
     &     ' ErrorElev =',e16.8)
#endif
C...
C...  READ AND PROCESS NABOUT - ABBREVIATED UNIT 16 OUTPUT OPTION
C...
      READ(15,*) NABOUT
      IF(NABOUT.EQ.1) THEN
         WRITE(16,3501) NABOUT
 3501    FORMAT(5X,'NABOUT = ',I2,
     &        /,9X,'ABREVIATED OUTPUT WILL BE PROVIDED TO UNIT 16',/,9X,
     &        'UNIT 14, 21, 22 INPUT DATA WILL NOT BE ECHO PRINTED',/)
      ELSE
         WRITE(16,3502) NABOUT
 3502    FORMAT(/,5X,'NABOUT = ',I3,
     &        /,9X,'DETAILED OUTPUT WILL BE PROVIDED TO UNIT 16',/,9X,
     &        'UNIT 14, 15, 21, 22 INPUT DATA WILL BE ECHO PRINTED',/)
      ENDIF

C...
C...  READ AND PROCESS NSCREEN - SCREEN OUTPUT OPTION
C...
C
C     jgf46.00 Added option to output data to the screen every NSCREEN
C     time steps, rather than on every time step.

C     jgf46.19 Added option to output "screen" data to fort.999 file
C     rather than the screen. This can be superior to a shell redirect
C     because on some platforms, the redirected log file is not
C     available until the run is complete.
C
      READ(15,*) NSCREEN
      ScreenUnit=0
      IF(NSCREEN.GT.0.AND.MYPROC.EQ.0) THEN
         ScreenUnit=6
         WRITE(16,3561) NSCREEN
 3561    FORMAT(5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO UNIT 6',
     &        /,9X,'EVERY NSCREEN TIME STEPS.',/)
      ELSEIF((NSCREEN.LT.0).AND.(MYPROC.EQ.0)) THEN
         ScreenUnit=999
         WRITE(16,3562) NSCREEN
 3562    FORMAT(/,5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO adcirc.log',/)
         OPEN(ScreenUnit,FILE=TRIM(GLOBALDIR)//'/'//'adcirc.log',
     &        STATUS='REPLACE')
      ELSE
         WRITE(16,3563) NSCREEN
 3563    FORMAT(/,5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL NOT BE PROVIDED',/)
      ENDIF

      IF (MYPROC.EQ.0.and.NScreen.ne.0) THEN
         WRITE(ScreenUnit,1112)
         WRITE(ScreenUnit,1114) TRIM(ADC_VERSION)
         WRITE(ScreenUnit,1112)
      ENDIF
C...
C...  READ AND PROCESS IHOT - HOT START OPTION
C...
      READ(15,*) IHOT
      WRITE(scratchMessage,'("IHOT=",I3,".")') IHOT
      CALL logMessage(ECHO,scratchMessage)

C     Logic to set the output unit number and output file name of the next
C     hotstart file to be written (so that the LUNs and file names alternate).
C
C     kmd48.33bc - added in the hot start option for the fort.17
C     jgf49.39: netcdf hotstart is 367 or 368. Changed if/then structure
C     to SELECT CASE to avoid confusing myself. The netcdf module makes
C     its own file name based on the unit number, and ignores hss%filename.
      SELECT CASE(IHOT)
      CASE(0,17,68,368,568)
         hss%lun = 67
         hss%filename = 'fort.67'
         IF ((IHOT.EQ.368).OR.(IHOT.EQ.568)) THEN
            useNetCDF = .true.
         ENDIF
      CASE(67,367,567)
         hss%lun = 68
         hss%filename = 'fort.68'
         IF ((IHOT.EQ.367).OR.(IHOT.EQ.567)) THEN
            useNetCDF = .true.
         ENDIF
      CASE DEFAULT
         write(scratchMessage,'("IHOT=",i0," is not valid.")') IHOT
         call allMessage(ERROR,scratchMessage)
         call ADCIRC_Terminate()
      END SELECT
#ifdef CSWAN
Casey 100205: Enable hot-start file generation by SWAN.
      SwanHotStartUnit = hss%lun
#endif

      IF(IHOT.NE.0) THEN
         WRITE(16,9733) IHOT
 9733    FORMAT(/,5X,'ADCIRC will be hot started using information ',
     &        'on UNIT ',I3,'.')
      ELSE
         CALL logMessage(INFO,'ADCIRC will be cold started.')
      ENDIF
C...
C...  READ AND PROCESS ICS - CARTESIAN/SPHERICAL COORDINATE OPTION
C...
      READ(15,*) ICS
      IF((ICS.NE.1).AND.(ICS.NE.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'ICS =',ICS
            WRITE(ScreenUnit,9735)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'ICS =',ICS
         WRITE(16,9735)
         WRITE(16,9973)
 9735    FORMAT(/,1X,'Your selection of ICS (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      IF(ICS.EQ.1) THEN
         WRITE(16,9736) ICS
 9736    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Cartesian coordinates')
      ELSE
         WRITE(16,9737) ICS
 9737    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Spherical coordinates',
     &        /,9X,'mapped using a CPP projection')
      ENDIF

C...
C...  READ AND PROCESS IM - 2D/3D MODEL FORMULATION OPTION
C...
      READ(15,*) IM
      WRITE(16,*) ' '
      WRITE(16,*) '    IM = ',IM
      WRITE(16,*) ' '
c     - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IM.LT.100) THEN
c     jgf Set defaults for model type (IM). All LOGICAL variables are
c     initialized to .FALSE. when declared in global.F.
         CGWCE_New        = .TRUE.
         CGWCE_LS_KGQ     = .TRUE.
         CGWCE_Advec_NC   = .TRUE.
c     jgf To use original momentum equations, uncomment the following
c     line and comment out the following two lines.
c     CME_Orig         = .TRUE. !uncomment for original momentum eqs
         CME_New_NC       = .TRUE. !comment out for original momentum eqs
         CME_LS_IBPV      = .TRUE. !comment out for original momentum eqs
         CME_AreaInt_Corr = .TRUE.
      ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IM.EQ.0) THEN
         C2DDI         = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 2DDI ',
     &        'run using: New GWCE and Momentum Eq formulations'
      ELSEIF (IM.EQ.10) THEN
         C2DDI         = .TRUE.
         C2D_PTrans    = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 2DDI ',
     &        'run using: New GWCE and Momentum Eq '
         WRITE(16,*) '          formulations + Passive Scalar Transport'
      ELSEIF (IM.EQ.20) THEN
         C2DDI         = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 2DDI ',
     &        'run using: New GWCE and Momentum Eq formulations'
      ELSEIF (IM.EQ.30) THEN
         C2DDI         = .TRUE.
         C2D_PTrans    = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 2DDI ',
     &        'run using: New GWCE and Momentum Eq '
         WRITE(16,*) '          formulations + Passive Scalar Transport'
      ELSEIF (IM.EQ.1) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 3D ',
     &        'run using: New GWCE and velocity based ',
     &        'Momentum Eqs.'
      ELSEIF (IM.EQ.11) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         C3D_PTrans    = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 3D ',
     &        'run using: New GWCE and velocity based '
         WRITE(16,*) '          Momentum Eqs + Passive Scalar Transport'
      ELSEIF (IM.EQ.21) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 3D ',
     &        'run using: New GWCE and velocity based ',
     &        'Momentum Eqs.'
      ELSEIF (IM.EQ.31) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         C3D_PTrans    = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 3D ',
     &        'run using: New GWCE and velocity based '
         WRITE(16,*) '          Momentum Eqs + Passive Scalar Transport'
      ELSEIF (IM.EQ.2) THEN
         C3D           = .TRUE.
         C3DDSS        = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a 3D run using',
     &        ': New GWCE and stress based Momentum Eqs.'
c     - - - - - - - - - - - - - - - - - - - - - - - -
c     f i n e   g r a i n e d   o p t i o n s ( i m )
c     - - - - - - - - - - - - - - - - - - - - - - - -
      ELSEIF ((IM.GE.111111).AND.(IM.LE.634322)) THEN
         IMDig1 = IM/100000
         IMDig2 = (IM - 100000*IMDig1)/10000
         IMDig3 = (IM - 100000*IMDig1 - 10000*IMDig2)/1000
         IMDig4 = (IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3)/100
         IMDig5 = (IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3
     &        -  100*IMDig4)/10
         IMDig6 =  IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3
     &        - 100*IMDig4 -   10*IMDig5

         C2DDI     = .TRUE.
         CGWCE_New = .TRUE.
         WRITE(16,*) '    ADCIRC is configured for a 2DDI run using'
         WRITE(16,*) '    the new GWCE routine and:'
         IF(IMDig1.EQ.1) THEN
            CGWCE_LS_KGQ     = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Kolar-Gray, flux based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.2) THEN
            CGWCE_LS_2PartQ  = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.3) THEN
            CGWCE_LS_2PartV  = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.4) THEN
            CGWCE_LS_2PartSQ  = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'symmetric stress in GWCE'
         ELSEIF(IMDig1.EQ.5) THEN
            CGWCE_LS_2PartSV  = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'symmetric stress in GWCE'
         ELSEIF(IMDig1.EQ.6) THEN
            C2DDI        = .FALSE.
            CGWCE_LS_KGQ  = .TRUE.
            C3D           = .TRUE.
            C3DVS         = .TRUE.
            ILump=0
          WRITE(16,*) '    ADCIRC is configured for a Barotropic 3D ',
     &        'run using: New GWCE and velocity based ',
     &        'Momentum Eqs.'
           WRITE(16,*) '        Kolar-Gray, flux based lateral ',
     &           'stress in GWCE'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig2.EQ.1) THEN
            CGWCE_Advec_NC   = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Non conservative advection in GWCE'
         ELSEIF(IMDig2.EQ.2) THEN
            CGWCE_Advec_C1   = .TRUE.
            WRITE(16,*) '        Conservative form 1 advection in GWCE'
         ELSEIF(IMDig2.EQ.3) THEN
            CGWCE_Advec_C2   = .TRUE.
            WRITE(16,*) '        Conservative form 2 advection in GWCE'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig3.EQ.1) THEN
            CME_LS_IBPV      = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Integration by parts, velocity based ',
     &           'lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.2) THEN
            CME_LS_IBPQ      = .TRUE.
            WRITE(16,*) '        Integration by parts, flux based ',
     &           'lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.3) THEN
            CME_LS_IBPSV      = .TRUE.
            WRITE(16,*) '        Integration by parts, velocity based ',
     &           'symmetric lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.4) THEN
            CME_LS_IBPSQ      = .TRUE.
            WRITE(16,*) '        Integration by parts, flux based ',
     &           'symmetric lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.5) THEN
            CME_LS_2PartV    = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.6) THEN
            CME_LS_2PartQ    = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'stress in Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig4.EQ.1) THEN
            CME_New_NC     = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Non conservative advection in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig4.EQ.2) THEN
            CME_New_C1       = .TRUE.
            WRITE(16,*) '        Conservative form 1 advection in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig4.EQ.3) THEN
            CME_New_C2     = .TRUE.
            WRITE(16,*) '        Conservative form 2 advection in ',
     &           'Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig5.EQ.1) THEN
            CME_AreaInt_Corr = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Corrected Area Integration in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig5.EQ.2) THEN
            CME_AreaInt_Orig = .TRUE.
            WRITE(16,*) '        Original Area Integration in ',
     &           'Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig6.EQ.2) THEN
            ILump=1
            CGWCE_Lump = .TRUE.
            WRITE(16,*) '        Lumped GWCE mass matrix'
         ELSE
            ILump=0
            WRITE(16,*) '        Consistent GWCE mass matrix'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'IM =',IM
            WRITE(ScreenUnit,9721)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'IM =',IM
         WRITE(16,9721)
         WRITE(16,9973)
 9721    FORMAT(/,1X,'Your selection of IM (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF

      IDEN=0
      if (CBaroclinic.eqv..true.) then
         READ(15,*) IDEN
         ! jgf51.52.35: Removed 0 as a valid choice for IDEN because the
         ! choice between baroclinic and barotropic is specified by IM; 
         ! also removed "2DDI" from the log messages because the choice 
         ! between 2DDI and 3D is specified by IM.  
         write(scratchMessage,'("IDEN is set to ",i0,".")') IDEN
         call logMessage(ECHO,scratchMessage)
         !
         densityRunType = ' diagnostic '
         if (iden.lt.0) then
            densityRunType = ' prognostic '
         endif
         !        
         densityDimensions = ' 2DDI '
         if (C3D.eqv..true.) then
            densityDimensions = ' 3D '
         endif
         !
         select case(abs(iden))
         case(1)
            densityForcingType = ' sigmaT '       
         case(2)
            densityForcingType = ' salinity '
         case(3)
            densityForcingType = ' temperature '
         case(4)
            densityForcingType = ' salinity and temperature '         
         case default
            call allMessage(ERROR,
     &       'The absolute value of IDEN must be 1, 2, 3, or 4.')
            call adcirc_terminate()
         end select
        call logMessage(INFO,'This run will be'//trim(densityDimensions)
     &    //trim(densityRunType)//'baroclinic with'
     &    //trim(densityForcingType)//'forcing.')
         !
         if ( (C2DDI.eqv..true.).and.(CBaroclinic.eqv..true.) ) then
            C2D_BTrans = .true.
         endif 
      endif


      WRITE(16,*) ' '

      WRITE(16,*) '     The ADCIRC logical variables are set to:'
      WRITE(16,*) '         C2DDI            = ',C2DDI
      WRITE(16,*) '         C3D              = ',C3D
      WRITE(16,*) '         C3DDSS           = ',C3DDSS
      WRITE(16,*) '         C3DVS            = ',C3DVS
      WRITE(16,*) '         C2D_BTrans       = ',C2D_BTrans
      WRITE(16,*) '         C2D_PTrans       = ',C2D_PTrans
c     WRITE(16,*) '         C3D_BTrans       = ',C3D_BTrans            !haven't yet read 3D input
      WRITE(16,*) '         C3D_PTrans       = ',C3D_PTrans
      WRITE(16,*) '         CBaroclinic      = ',CBaroclinic
      WRITE(16,*) '         CGWCE_Lump       = ',CGWCE_Lump
      WRITE(16,*) '         CGWCE_LS_KGQ     = ',CGWCE_LS_KGQ
      WRITE(16,*) '         CGWCE_LS_2PartQ  = ',CGWCE_LS_2PartQ
      WRITE(16,*) '         CGWCE_LS_2PartV  = ',CGWCE_LS_2PartV
      WRITE(16,*) '         CGWCE_LS_2PartSQ = ',CGWCE_LS_2PartSQ
      WRITE(16,*) '         CGWCE_LS_2PartSV = ',CGWCE_LS_2PartSV
      WRITE(16,*) '         CGWCE_Advec_NC   = ',CGWCE_Advec_NC
      WRITE(16,*) '         CGWCE_Advec_C1   = ',CGWCE_Advec_C1
      WRITE(16,*) '         CGWCE_Advec_C2   = ',CGWCE_Advec_C2
      WRITE(16,*) '         CME_Orig         = ',CME_Orig
      WRITE(16,*) '         CME_New_NC       = ',CME_New_NC
      WRITE(16,*) '         CME_New_C1       = ',CME_New_C1
      WRITE(16,*) '         CME_New_C2       = ',CME_New_C2
      WRITE(16,*) '         CME_LS_IBPQ      = ',CME_LS_IBPQ
      WRITE(16,*) '         CME_LS_IBPV      = ',CME_LS_IBPV
      WRITE(16,*) '         CME_LS_IBPSQ     = ',CME_LS_IBPSQ
      WRITE(16,*) '         CME_LS_IBPSV     = ',CME_LS_IBPSV
      WRITE(16,*) '         CME_LS_2PartQ    = ',CME_LS_2PartQ
      WRITE(16,*) '         CME_LS_2PartV    = ',CME_LS_2PartV
      WRITE(16,*) '         CME_AreaInt_Orig = ',CME_AreaInt_Orig
      WRITE(16,*) '         CME_AreaInt_Corr = ',CME_AreaInt_Corr
      WRITE(16,*) '         CTIP             = ',CTIP
      WRITE(16,*) '         CHARMV           = ',CHARMV

C...
C...  READ AND PROCESS NOLIBF - NONLINEAR BOTTOM FRICTION OPTION
C...
      READ(15,*) NOLIBF
      IF((NOLIBF.LT.0).OR.(NOLIBF.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLIBF =',NOLIBF
            WRITE(ScreenUnit,9722)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLIBF =',NOLIBF
         WRITE(16,9722)
         WRITE(16,9973)
 9722    FORMAT(/,1X,'Your selection of NOLIBF (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      WRITE(16,9845) NOLIBF
 9845 FORMAT(/,5X,'NOLIBF = ',I3)
      IF (NOLIBF.EQ.0) WRITE(16,2050)
 2050 FORMAT(9X,'THE MODEL WILL USE LINEAR BOTTOM FRICTION')
      IF (NOLIBF.EQ.1) WRITE(16,2051)
 2051 FORMAT(9X,'THE MODEL WILL USE NONLINEAR BOTTOM FRICTION')
      IF (NOLIBF.EQ.2) WRITE(16,2052)
 2052 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION',
     &     'IN DEEP WATER ',
     &     /,9X,'AND A FRICTION FACTOR THAT INCREASES AS THE DEPTH ',
     &     'DECREASES IN SHALLOW WATER')
C...
C... READ AND PROCESS NOLIFA - NONLINEAR FINITE AMPLITUDE OPTION
C...
      READ(15,*) NOLIFA
      IF((NOLIFA.LT.0).OR.(NOLIFA.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)       ! input error
            WRITE(ScreenUnit,*) 'NOLIFA =',NOLIFA
            WRITE(ScreenUnit,9723)       ! not allowable
         ENDIF
         WRITE(16,9972)         ! input error
         WRITE(16,*) 'NOLIFA =',NOLIFA
         WRITE(16,9723)         ! not allowable
 9723    FORMAT(/,1X,'Your selection of NOLIFA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         IF (NoLiFA.eq.3.and.NFOver.eq.1) THEN
            WRITE(16,8735)
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,8735)
 8735       FORMAT(/,1X,'WARNING: The StartDry file was replaced ',
     &           'by surface_submergence_state in ',
     &           /,1X,'the Nodal Attributes file (unit 13).'
     &           //,1X,'ACTION: NOLIFA will be corrected to 2; the ',
     &           'loading of StartDry data will not ',
     &           /,1X,'be triggered now, although ',
     &           'it may be triggered later in the NWP section.',/)
            NoLiFA = 2
         ELSE
            IF (NoLiFA.eq.3.and.NFOver.eq.0) THEN
               WRITE(16,7624)   ! startdry replaced
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,7624)
            ENDIF
 7624       FORMAT(/,1X,'ERROR: NOLIFA=3 formerly triggered the ',
     &           'loading of StartDry data.'
     &           /,1X,'However, the StartDry file was replaced ',
     &           'by surface_submergence_state in ',
     &           /,1X,'the Nodal Attributes file (unit 13). Please ',
     &           'use NWP to load this data.',/)
            WRITE(16,9973)      ! execution will terminate
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            CALL ADCIRC_Terminate()
         ENDIF
      ENDIF
      WRITE(16,9846) NOLIFA
 9846 FORMAT(/,5X,'NOLIFA = ',I3)
      IF(NOLIFA.EQ.0) WRITE(16,2053)
 2053 FORMAT(9X,'THE MODEL WILL NOT USE FINITE AMPLITUDE TERMS OR ',
     &     'WETTING AND DRYING')
      IF(NOLIFA.EQ.1) WRITE(16,2054)
 2054 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS BUT NO ',
     &     'WETTING AND DRYING')
      IF(NOLIFA.EQ.2) WRITE(16,2049)
 2049 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',
     &     'WETTING AND DRYING')
C...
C...  READ AND PROCESS NOLICA - ADVECTIVE TERM SPATIAL GRADIENT
C...
      READ(15,*) NOLICA
      IF((NOLICA.LT.0).OR.(NOLICA.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICA =',NOLICA
            WRITE(ScreenUnit,9724)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICA =',NOLICA
         WRITE(16,9724)
         WRITE(16,9973)
 9724    FORMAT(/,1X,'Your selection of NOLICA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      WRITE(16,9847) NOLICA
 9847 FORMAT(/,5X,'NOLICA = ',I3)
      IF(NOLICA.EQ.0) WRITE(16,2055)
 2055 FORMAT(9X,'THE MODEL WILL NOT USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')
      IF(NOLICA.EQ.1) WRITE(16,2056)
 2056 FORMAT(9X,'THE MODEL WILL USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')

C...
C...  READ AND PROCESS NOLICAT - GWCE ADVECTIVE TERM TIME DERIVATIVE
C...
      READ(15,*) NOLICAT
      IF((NOLICAT.LT.0).OR.(NOLICAT.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9725)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9725)
         WRITE(16,9973)
 9725    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      !
      ! jgf52.05: Removed this error message because it is alarming
      ! and because it contradicts common practice. 
C      IF((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
C         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
C            WRITE(ScreenUnit,9972)
C            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
C            WRITE(ScreenUnit,9726)
C            IF(NFOVER.EQ.1) THEN
C               WRITE(ScreenUnit,9974)
C            ELSE
C               WRITE(ScreenUnit,9973)
C            ENDIF
C         ENDIF
C         WRITE(16,9972)
C         WRITE(16,*) 'NOLICAT =',NOLICAT
C         WRITE(16,9726)
C         WRITE(16,9974)
 9726    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLIFA and may lead to mass ',
     &        'balance problems')
C         IF(NFOVER.EQ.1) THEN
C            if (myproc == 0) WRITE(ScreenUnit,9974)
C         ELSE
C            if (myproc == 0) WRITE(ScreenUnit,9973)
C            call ADCIRC_Terminate()
C         ENDIF
C      ENDIF

      IF((NOLIFA.EQ.0).AND.(NOLICAT.EQ.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9726)
            IF(NFOVER.EQ.1) THEN
               WRITE(ScreenUnit,9974)
            ELSE
               WRITE(ScreenUnit,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9727)
         WRITE(16,9974)
         IF(NFOVER.EQ.1) THEN
            if (myproc == 0) WRITE(ScreenUnit,9974)
         ELSE
            if (myproc == 0) WRITE(ScreenUnit,9973)
            CALL ADCIRC_Terminate()
         ENDIF
      ENDIF

      IF(NOLICA.NE.NOLICAT) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9726)
            IF(NFOVER.EQ.1) THEN
               WRITE(ScreenUnit,9974)
            ELSE
               WRITE(ScreenUnit,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9727)
         WRITE(16,9974)
 9727    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLICA and may lead to mass ',
     &        'balance problems')
         IF(NFOVER.EQ.1) THEN
            if (myproc == 0) WRITE(ScreenUnit,9974)
         ELSE
            if (myproc == 0) WRITE(ScreenUnit,9973)
            call ADCIRC_Terminate()
         ENDIF
      ENDIF

      WRITE(16,9848) NOLICAT
 9848 FORMAT(/,5X,'NOLICAT = ',I3)
      IF(NOLICAT.EQ.0) WRITE(16,2057)
 2057 FORMAT(9X,'THE MODEL WILL NOT USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
      IF(NOLICAT.EQ.1) WRITE(16,2058)
 2058 FORMAT(9X,'THE MODEL WILL USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
C
C     READ AND PROCESS NWP jgf46.00 Read in nodal attributes such as
C     tau0, bottom friction, directional wind speed reduction factor,
C     startdry, etc. The full initialization and error checking of these
C     data must wait until the grid has been read in from unit 14.
      READ(15,*) NWP
      CALL ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
C...
C...  READ AND PROCESS NCOR - SPATIALLY VARYING CORIOLIS PARAMETER
C...
      READ(15,*) NCOR
      IF((NCOR.NE.0).AND.(NCOR.NE.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCOR =',NCOR
            WRITE(ScreenUnit,9729)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9729)
         WRITE(16,9973)
 9729    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      IF((ICS.EQ.1).AND.(NCOR.EQ.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCOR =',NCOR
            WRITE(ScreenUnit,9730)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9730)
         WRITE(16,9973)
 9730    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Spatially ',
     &        'variable Coriolis should be used only with ',
     &        /,1X,'Spherical coordinates')
         CALL ADCIRC_Terminate()
      ENDIF
      IF(NCOR.EQ.0) THEN
         WRITE(16,233) NCOR
 233     FORMAT(/,5X,'NCOR = ',I2,
     &    /,9X,'A CONSTANT VALUE OF THE CORIOLIS PARAMETER WILL BE ',
     &    /,9X,'USED THROUGHOUT THE DOMAIN')
      ELSE
         WRITE(16,234) NCOR
 234     FORMAT(/,5X,'NCOR = ',I2,
     &    /,9X,'SPATIALLY VARYING CORIOLIS VALUES WILL BE COMPUTED ',
     &    'FROM INPUT LATITUDES')
      ENDIF

C...
C...  READ AND PROCESS NTIP - TIDAL POTENTIAL FORCING
C...
      READ(15,*) NTIP
      IF((NTIP.LT.0).OR.(NTIP.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NTIP =',NTIP
            WRITE(ScreenUnit,9710)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9710)
         WRITE(16,9973)
 9710    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      IF((ICS.EQ.1).AND.(NTIP.GE.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NTIP =',NTIP
            WRITE(ScreenUnit,9711)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9711)
         WRITE(16,9973)
 9711    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Tidal',
     &        'potential forcing should be used only with ',
     &        /,1X,'Spherical coordinates')
         CALL ADCIRC_Terminate()
      ENDIF
      IF (NTIP.NE.0) CTIP = .TRUE.
      IF(NTIP.EQ.0) THEN
         WRITE(16,235) NTIP
 235     FORMAT(/,5X,'NTIP = ',I2,/,9X,
     &        'TIDAL POTENTIAL FORCING IS NOT USED IN THE COMPUTATION')
      ENDIF
      IF(NTIP.GE.1) THEN
         WRITE(16,236) NTIP
 236     FORMAT(/,5X,'NTIP = ',I2,
     &    /,9X,'TIDAL POTENTIAL FORCING IS USED IN THE COMPUTATION ',
     &   'BASED ON INPUT LONGITUDES/LATITUDES')
      ENDIF
      IF(NTIP.EQ.2) THEN
         WRITE(16,239)
 239     FORMAT(9X,'SELF ATTRACTION/LOAD TIDE FORCING IS ALSO USED ',
     &        'IN THE COMPUTATION')
      ENDIF
C...
C...  READ AND PROCESS NWS - WIND AND PRESSURE FORCING & WAVE RADIATION
C...  STRESS FORCING
C...
C     jgf46.00 Added NWS=7 (direct surface stress).
C     jgf46.03 Added NWS=8 (Holland wind model)
C     jgf46.16 Merged:
C     cf & cm  Added NWS=9 (asymmetric hurricane winds)
C     rjw      Added NWS=19(asymmetric hurricane winds v2.0)
C     jie      ADDed NWS=20(generalized asymmetric vortex winds)
C     sb46.28sb01 Added NWS=12 (OWI format) 09/xx/2006
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
C     tcm v49.46 Added NWS = 4xx for tight coupling with STWAVE
C     tcm v49.64.01 Added NWS = 1xxx's for ice concentration
      NCICE = 0  !set ice type to be 0
C
      READ(15,FMT=*,END=99998,ERR=99999) NWS
C
C.... BREAK OUT THE ICE CONCENTRATION FLAG FROM NWS/NRS
      IF (NWS.EQ.0) THEN
         NCICE = 0
      ELSE
         NCICE = INT(ABS(NWS)/1000)
         NWS = INT(ABS(NWS)-NCICE*1000)*INT(NWS/ABS(NWS))  !RESETTING NWS/NRS
                                                           !BY REMOVING THE 1000'S
                                                           !PLACE FOR ICE
      ENDIF
      nwsOK = .false.
      do i=0,4
         if ( any(nws.eq.(allowable_nws + i*100)).or.
     &        any(nws.eq.(-1*(allowable_nws + i*100))) ) then
            nwsOK = .true.
            exit
         endif
      end do
      if (nwsOK.eqv..false.) then
         write(scratchMessage,
     &        '("NWS=",I3," is not an allowable value.")') NWS
         call allMessage(ERROR,scratchMessage)
         call ADCIRC_Terminate()
      endif

C.... SET WAVE RADIATION STRESS FLAG AND ADJUST NWS ACCORDINGLY

!       NRS=0
!       IF(ABS(NWS/100).EQ.1) THEN ! sb46.28sb03
!          NRS=1
!          NWS=(ABS(NWS)-100)*(NWS/ABS(NWS))
!       ENDIF
! C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
!       IF(ABS(NWS/100).EQ.2) THEN
!          NRS=2
!          NWS=(ABS(NWS)-200)*(NWS/ABS(NWS))
!       ENDIF
! #ifdef CSWAN
! Casey 090302: Added the option for coupling directly to SWAN.
!       IF(ABS(NWS/100).EQ.3) THEN
!          NRS=3
!          NWS=(ABS(NWS)-300)*(NWS/ABS(NWS))
!       ENDIF
! #endif
C.... tcm v49.46 rewrote to combine different 100's power NRS
      NRS = 0
      IF(NWS.EQ.0) THEN
         NWS = 0
         NRS = 0
      ELSE
         NRS=ABS(NWS/100)
         NWS=(ABS(NWS)-NRS*100)*(NWS/ABS(NWS))
      ENDIF

      IF(NWS.EQ.0) THEN
         WRITE(16,237) NWS
 237     FORMAT(/,5X,'NWS = ',I3,/,9X,
     &        'WIND STRESS OR SURFACE PRESSURE ARE NOT USED TO FORCE',
     &        'THE COMPUTATION')
      ENDIF
      IF(NWS.EQ.1) THEN
         WRITE(16,238) NWS
 238     FORMAT(/,5X,'NWS = ',I3,
     &        /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &        /,9X,' THE COMPUTATION',
     &        /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &        /,9X,' EVERY MODEL TIME STEP')
      ENDIF
      IF(NWS.EQ.2) THEN
         WRITE(16,2381) NWS
 2381    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-2) THEN
         WRITE(16,2380) NWS
 2380    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.3) THEN
         WRITE(16,2382) NWS
 2382    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS ONLY IS USED TO FORCE THE COMPUTATION.',
     &    /,9X,'WIND SPEEDS AND DIRECTIONS ARE READ FROM A FLEET ',
     &    /,9X,'NUMERIC FORMAT FILE AT UNIT 22 AND INTERPOLATED TO',
     &    /,9X,'THE ADCIRC GRID. ',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &     'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.4) THEN
         WRITE(16,2383) NWS
 2383    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-4) THEN
         WRITE(16,2388) NWS
 2388    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.5) THEN
         WRITE(16,2384) NWS
 2384    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-5) THEN
         WRITE(16,2389) NWS
 2389    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.6) THEN
         WRITE(16,2385) NWS
 2385    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM A ',
     &    /,9X,'REGULARLY SPACED GRID FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A REGULAR GRID TO THE ADCIRC GRID.'
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
C     jgf46.00 Added NWS=7 (direct surface stress).
      IF(NWS.EQ.7) THEN
         WRITE(16,1234) NWS
 1234    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'SURFACE STRESS AND PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STRESS DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-7) THEN
         WRITE(16,1235) NWS
 1235    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'SURFACE STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STRESS DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
C     jgf46.03 Added NWS=8 (Holland wind model).
      IF(NWS.EQ.8) THEN
         WRITE(16,1237) NWS
 1237    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'HURRICANE PARAMETERS AND THE HOLLAND WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-8) THEN
         WRITE(16,1238) NWS
 1238    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'HURRICANE PARAMETERS AND THE HOLLAND WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.10) THEN
         WRITE(16,2386) NWS
 2386    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY N',
     &    /,9X,' HOURS FROM A DIFFERENT FILE AT UNITS 200, 200+N,',
     &   ' 200+2N, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A GAUSSIAN GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.11) THEN
         WRITE(16,2387) NWS
 2387    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY 3 ',
     &    /,9X,'HOURS FROM ETA-29 FILES AT UNITS 200, 201, 202, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'WIND DATA FROM THE 29 KM E GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &   'DRAG LAW.')
      ENDIF
C.....sb46_28sb01 added for NWS=-12,12 09/xx/2006
      IF(NWS.EQ.12) THEN
         WRITE(16,12384) NWS
12384    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-12) THEN
         WRITE(16,12389) NWS
12389    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
      IF (ABS(NWS).EQ.15) THEN
         WRITE(16,10389) NWS
10389    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY VALUES ARE READ FROM RAW ',
     &    /,9X,'HWIND FILES (UNIT 220). DATA WILL BE INTERPOLATED IN',
     &    /,9X,'SPACE BY REPROJECTING THE ADCIRC MESH IN MERCATOR.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
C
      IF(NWS.EQ.16) THEN
         WRITE(16,12374) NWS
12374    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'GFDL MET DATA FILES (UNIT 220).',
     &    /,9X,'META DATA IS READ FROM UNIT 22.',
     &    /,9X,'INTERPOLATION IN SPACE IS DONE WITH A WEIGHTED ',
     &    /,9x,'DISTANCE NEAREST NEIGHBOR METHOD.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-16) THEN
         WRITE(16,12379) NWS
12379    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'GFDL MET DATA FILES (UNIT 220).',
     &    /,9X,'META DATA IS READ FROM UNIT 22.',
     &    /,9X,'INTERPOLATION IN SPACE IS DONE WITH A WEIGHTED ',
     &    /,9x,'DISTANCE NEAREST NEIGHBOR METHOD.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
C
C     rjw added nws = 19: asymmetric hurricane winds v2.0
      IF(NWS.EQ.19) THEN
         WRITE(16,2400) NWS
 2400    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM TEH ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-19) THEN
         WRITE(16,2401) NWS
 2401    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM TEH ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

C     jie added nws = 20: generalized asymmetric wind model
      IF(NWS.EQ.20) THEN
         WRITE(16,2404) NWS
 2404    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'THE GENERALIZED ASYMMETRIC VORTEX MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM TEH ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-20) THEN
         WRITE(16,2405) NWS
 2405    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'HTHE GENERALIZED ASYMMETRIC VORTEX MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM TEH ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

C....Surface Wave Stresses
C     jgf49.1001 Added option to embed NWS19 in a NAM background wind field
      IF(NWS.EQ.29) THEN
         WRITE(16,2402) NWS
 2402    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM THE ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.'
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-29) THEN
         WRITE(16,2403) NWS
 2403    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ IN FOR THE STORM FROM UNIT 22',
     &    /,9X,
     &    'WHICH IS CREATED FROM THE ATCF FILE USING THE ASWIP PROGRAM',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.'
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NRS.EQ.0) THEN
         WRITE(16,2390) NRS
 2390    FORMAT(/,5X,'NRS = ',I3,
     &        /,9X,'WAVE RADIATION STRESS IS NOT USED TO FORCE THE ',
     &        'COMPUTATION')
      ENDIF
      IF(NRS.EQ.1) THEN
         WRITE(16,2391) NRS
 2391    FORMAT(/,5X,'NRS = ',I3,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',
     &    /,9X,'PBL TYPE FILE AT UNIT 23.  INTERPOLATION IN TIME IS ',
     &    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
      IF(NRS.EQ.2) THEN
         WRITE(16,2392) NRS
 2392    FORMAT(/,5X,'NRS = ',I3,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',
     &    /,9X,'UNIT 23 FILE GENERATED BY UTIL/BUILDSTWAVE23. ',
     &    /,9X,'NO DECOMPOSITION IS NEEDED FOR THIS UNIT 23 FILE EVEN',
     &    /,9X,'IN A PARALLEL EXCECUTION.  INTERPOLATION IN TIME IS ',
     &    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
#ifdef CSWAN
Casey 090302: Added the following lines.
      IF(NRS.EQ.3) THEN
         WRITE(16,2393) NRS
 2393    FORMAT(/,5X,'NRS = ',I3,
     &    /,9X,'WAVES WILL BE COUPLED TO SWAN!')
      ENDIF
#endif

C     NWS=4xx added for STWAVE tighly coupled run  (v49.46 tcm)
      IF(NRS.EQ.4) THEN
         CPL2STWAVE = .TRUE.  !flag indicating coupling with STWAVE 
         WRITE(16,2394) NRS
 2394    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE COMPUTED ON THE FLY BY STWAVE COMPUTE    ',
     &    /,9X,'PROCS.')
      ENDIF
C
C------------------ BEGIN ICE CONCENTRATION ---------------------------C
C
C     TCM V49.64.01 ADDITION FOR ICE CONCENTRATION FIELDS

!     TEST TO BE SURE NCICE HAS AN ALLOWABLE VALUE
      IF((NCICE.NE.0).AND.(NCICE.NE.12)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCICE =',NCICE
            WRITE(ScreenUnit,9812)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCICE =',NCICE
         WRITE(16,9812)
         WRITE(16,9973)
 9812    FORMAT(/,1X,'Your selection of NCICE (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
C...  TCM v49.64.02 -- added
C...  BE SURE NWS AND NCICE ARE COMPATABLE
      IF((NCICE.GT.0).AND.
     &       ((NWS.EQ.1).OR.(NWS.EQ.2).OR.(NWS.EQ.7))) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCICE =',NCICE
            WRITE(ScreenUnit,*) 'NWS =',NWS
            WRITE(ScreenUnit,9813)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCICE =',NCICE
         WRITE(16,*) 'NWS =',NWS
         WRITE(16,9813)
         WRITE(16,9973)
 9813    FORMAT(/,1X,'Your selection of NCICE (a UNIT 15 input ',
     &    'parameter) is not compatable',
     &    /,1x,'with your NWS value. ',
     &    'NCICE is not allowed for use with abs(NWS)=1,2, or 7.')
         CALL ADCIRC_Terminate()
      ENDIF

      IF(NCICE.EQ.0) THEN
         WRITE(16,3294) NCICE
 3294    FORMAT(/,5X,'NCICE = ',I2,
     &        /,9X,'ICE CONCENTRATION FIELDS ARE NOT USED TO ',
     &        'ADJUST WIND STRESS COMPUTATIONS')
      ENDIF
      IF(NCICE.EQ.12) THEN
         WRITE(16,3295) NCICE
 3295    FORMAT(/,5X,'NCICE = ',I2,
     &    /,9X,'ICE CONCENTRATION FIELDS ARE USED TO ',
     &        'ADJUST WIND STRESS COMPUTATIONS'
     &    /,9X,'ICE CONCENTRATION FIELDS ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 225).',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE ICE DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 225 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 225 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
C
C------------------ END ICE CONCENTRATION ------------------------------C
C

C... tcm v50.66.02 addition for time varying bathymetry
C
C------------------ BEGIN TIME VARYING BATHY---------------------------C
C
C   The value of NDDT is set in the TimeBathyControl Namelist if present
C
      IF ( (ABS(NDDT).NE.0).AND.(abs(NDDT).NE.1)
     &                     .AND.(abs(NDDT).NE.2) ) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NDDT = ',NDDT
            WRITE(ScreenUnit,9814)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NDDT =',NDDT
         WRITE(16,9812)
         WRITE(16,9973)
 9814    FORMAT(/,1X,'Your selection of NDDT (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL ADCIRC_Terminate()
      ENDIF
      IF(NDDT.EQ.0) THEN
         WRITE(16,837) NDDT
 837     FORMAT(/,5X,'NDDT = ',I2,/,9X,
     &        'A TIME VARYING BATHYMETRY IS NOT USED DURING',
     &        ' THE COMPUTATION')
      ENDIF
      IF(NDDT.EQ.1) THEN
         WRITE(16,838) NDDT
 838     FORMAT(/,5X,'NDDT = ',I2,
     &        /,9X,'A TIME VARYING BATHYMETRY IS USED DURING',
     &             ' THE COMPUTATION',
     &    /,9X,'NEW BATHYMETRY VALUES ARE READ AT ALL',
     &    /,9X,'ADCIRC GRID NODES FROM A UNIT 141.',
     &    /,9X,'THE UNIT 141 FILE BEGINS AT TIME=STATIM+BTIMINC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE BATHYMETRY ',
     &    /,9X,'DATA WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NDDT.EQ.-1) THEN
         WRITE(16,839) NDDT
 839     FORMAT(/,5X,'NDDT = ',I2,
     &        /,9X,'A TIME VARYING BATHYMETRY IS USED DURING',
     &             ' THE COMPUTATION',
     &    /,9X,'NEW BATHYMETRY VALUES ARE READ AT ALL',
     &    /,9X,'ADCIRC GRID NODES FROM A UNIT 141.',
     &    /,9X,'THE UNIT 141 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'PLUS BTIMINC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE BATHYMETRY ',
     &    /,9X,'DATA WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NDDT.EQ.2) THEN
         WRITE(16,840) NDDT
 840     FORMAT(/,5X,'NDDT = ',I2,
     &        /,9X,'A TIME VARYING BATHYMETRY IS USED DURING',
     &             ' THE COMPUTATION',
     &    /,9X,'NEW BATHYMETRY VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A UNIT 141.',
     &    /,9X,'THE UNIT 141 FILE BEGINS AT TIME=STATIM+BTIMINC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE BATHYMETRY ',
     &    /,9X,'DATA WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NDDT.EQ.-2) THEN
         WRITE(16,841) NDDT
 841     FORMAT(/,5X,'NDDT = ',I2,
     &        /,9X,'A TIME VARYING BATHYMETRY IS USED DURING',
     &             ' THE COMPUTATION',
     &    /,9X,'NEW BATHYMETRY VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A UNIT 141.',
     &    /,9X,'THE UNIT 141 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'PLUS BTIMINC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE BATHYMETRY ',
     &    /,9X,'DATA WITH THE MODEL TIME STEP.')
      ENDIF
C...
C...
C
C------------------ END TIME VARYING BATHY ------------------------------C
C

C...
C...  READ AND PROCESS NRAMP - WHETHER A RAMP FUNCTION WILL BE USED
C...
C     jgf46.08 Change to the number of ramp functions that will be used.
      READ(15,*) NRAMP
      IF((NRAMP.NE.0).AND.(NRAMP.gt.8)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NRAMP =',NRAMP
            WRITE(ScreenUnit,9713)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
 9713    FORMAT(/,1X,'Your selection of NRAMP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      IF(NRAMP.EQ.0) THEN
         WRITE(16,240) NRAMP
 240     FORMAT(/,5X,'NRAMP = ',I2,
     &        /,9X,'NO RAMP FUNCTION IS USED IN THE COMPUTATION')
      ELSE
         WRITE(16,241) NRAMP
 241     FORMAT(/,5X,'NRAMP = ',I2,
     &    /,9X,'HYPERBOLIC TANGENT RAMP(S) WILL BE APPLIED TO THE ',
     &    'FORCING FUNCTIONS')
      ENDIF
C...
C...  PROCESS G - GRAVITY
C...
      READ(15,*) G
      IF((ICS.EQ.2).AND.(abs(G-9.81).gt.0.01)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'G =',G
            WRITE(ScreenUnit,9714)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'G =',G
         WRITE(16,9714)
         WRITE(16,9973)
 9714    FORMAT(/,1X,'Your specification of the gravitational ',
     &        'constant, G, (a UNIT 15 input) is not ',
     &        /,1X,'consistant with the use of spherical coordinates.',
     &        '  G must be in units of m/s^2')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      WRITE(16,5) G
    5 FORMAT(///,5X,'GRAVITATIONAL CONSTANT G =',F10.5,/)

C...
C...  READ AND PROCESS TAU0 - WEIGHTING COEFFICIENT IN THE GWCE
C
C...  jgf45.12 Added three-tier tau0 scheme.
C
C     jgf46.00 Added user specified spatially varying tau0 scheme using
C     the Nodal Attributes File (unit 13) and removed three tier tau0
C     scheme. Since Tau0 is a nodal attribute, it is initialized later
C     in this subroutine using a CALL to the InitNodalAttr subroutine in
C     the NodalAttributes module.
      READ(15,*) Tau0
C
C     jgf47.11 Added a line to the fort.15 file to read in the min
C     and max tau0 values, if the user has chosen to use the time
C     varying tau0.
      IF ( (Tau0.le.-5.d0).and.(Tau0.gt.-6.d0) ) THEN
        READ(15,*) Tau0FullDomainMin, Tau0FullDomainMax
      ENDIF
C...
C...  INPUT FROM UNIT 15 AND OUTPUT TO UNIT 16 TIME INTEGRATION
C...  INFORMATION INCLUDING DT,STATIM,REFTIM,AND RNDAY
C...
      WRITE(16,1112)
      WRITE(16,245)
  245 FORMAT(//,1X,'TIME INTEGRATION INFORMATION',//)

C...
C...  READ AND PROCESS DT - MODEL TIME STEP
C...
ckmd   Changed the time step to allow for negative values in
ckmd   order to turn on the predictor-corrector algorithm.
      READ(15,*) DTDP
      IF(DTDP.LT.0.d0) THEN
         CPRECOR    = .TRUE.
         CGWCE_New  = .FALSE. !jgf Turn off the default.
         CME_New_NC = .FALSE. !jgf Turn off the default.
         CME_New_C1 = .FALSE.
         CME_New_C2 = .FALSE.
         DT=-DTDP
         DTDP=DT
         WRITE(16,*) ' ADCIRC is configured for a 2DDI run using'
         WRITE(16,*) ' the predictor-corrector algorithm and'
         WRITE(16,*) ' the ADCIRC logical variable is set to:    '
         WRITE(16,*) '         CPRECOR           = ',CPRECOR
      ELSE IF(DTDP.GT.0.d0) THEN
         DT=DTDP
         WRITE(16,*) ' ADCIRC is configured for a 2DDI run '
         WRITE(16,*) ' without the predictor-corrector algorithm and'
         WRITE(16,*) ' the ADCIRC logical variable is set to:  '
         WRITE(16,*) '         CPRECOR           = ',CPRECOR
      END IF
      WRITE(16,9) DTDP
    9 FORMAT(5X,'TIME STEP =',F12.6,5X,'SECONDS',/)
      DTDPHS = DTDP ! jgf51.14: Initialize; may be changed in hstart()     

C...
C...  READ AND PROCESS STATIM - SIMULATION STARTING TIME
C...
      READ(15,*) STATIM
      WRITE(16,1113) STATIM
 1113 FORMAT(5X,'STARTING TIME FOR SIMULATION = ',F14.6,' DAYS',/)

C...
C...  READ AND PROCESS REFTIM - Harmonic REFERNCE TIME
C...
      READ(15,*) REFTIM
      WRITE(16,1115) REFTIM
 1115 FORMAT(5X,'Harmonic REFERENCE TIME = ',F14.6,' DAYS',/)

C...
C...  Read in and process additional timing information for wind.
C...
      IF((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
      IF((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
      IF(ABS(NWS).EQ.2) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
      IF(NWS.EQ.3) THEN
         READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC
         WRITE(16,1116) IREFMO,IREFDAY,IREFYR,IREFHR,IREFMIN,REFSEC
 1116    FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &        I2,'/',I2,'/',I2,'  ',I2,':',I2,':',f7.4,/)
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC,
     &        WREFTIM, MyProc, NScreen, ScreenUnit)
C...     TCM V49.64.01 CHANGES FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) NWLAT,NWLON,WLATMAX,
     &        WLONMIN,WLATINC,WLONINC,WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) NWLAT,NWLON,WLATMAX,  ! sb46.28sb03
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,RSTIMINC
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) NWLAT,NWLON,WLATMAX,
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) NWLAT,NWLON,WLATMAX,  ! sb46.28sb03
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
      ENDIF
      IF(ABS(NWS).EQ.4) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
      ENDIF
      IF(ABS(NWS).EQ.5) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,CICE_TIMINC
      ENDIF
      IF(NWS.EQ.6) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) NWLAT,NWLON,WLATMAX,
     &        WLONMIN,WLATINC,WLONINC,WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) NWLAT,NWLON,WLATMAX, ! sb46.28sb03
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,RSTIMINC
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) NWLAT,NWLON,WLATMAX,
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) NWLAT,NWLON,WLATMAX, ! sb46.28sb03
     &        WLONMIN,WLATINC,WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
      ENDIF
C     jgf46.00 Added NWS=7 (direct surface stress).
      IF(ABS(NWS).EQ.7) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
C     jgf46.05 Added NWS=8 (Holland Wind Model).
C     jgf46.28 Changed WTIMINC to StormNumber for activating
C     wind multiplier to final wind speeds from Holland model.
      IF(ABS(NWS).EQ.8) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION AND FIXED BUG
C...     WHERE NRS=2, NRS=4 WERE NOT BEING INCLUDED
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC,CICE_TIMINC
         ENDIF
         WRITE(16,6111) IREFMO,IREFDAY,IREFYR,IREFHR
 6111    FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &           I2,'/',I2,'/',I2,'  ',I2,'H',/)
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,0,0.0d0,
     &        WindRefTime, MyProc, NScreen, ScreenUnit)
      ENDIF
C     jgf46.16 Merged:
C     cf & cm added nws = 9: asymmetric hurricane winds
      IF(ABS(NWS).EQ.9) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,
     &                                                CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,
     &                                                CICE_TIMINC
         WRITE(16,1117) WTIMINC
      ENDIF
      IF(NWS.EQ.10) THEN
         NWLAT=190
         NWLON=384
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,
     &                                                 CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC, ! sb46.28sb03
     &                                                 CICE_TIMINC
      ENDIF
      IF(NWS.EQ.11) THEN
         NWLAT=271
         NWLON=181
         WTIMINC=10800.
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMINC
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) RSTIMINC,
     &                                                 CICE_TIMINC
         !READ(15,*) NWLAT,NWLON,WTIMINC
      ENDIF
C     sb46.28sb01 Added NWS=12 (OWI format)
      IF(ABS(NWS).EQ.12) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,
     &                                                 CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC, ! sb46.28sb03
     &                                                 CICE_TIMINC
      ENDIF
C     rjw added nws = 19: asymmetric hurricane winds
      IF(ABS(NWS).EQ.19) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC,CICE_TIMINC
         ELSE
         ENDIF
C LETS not use this now but have the option to i nthe future
c            WRITE(16,6111) IREFMO,IREFDAY,IREFYR,IREFHR
c 6112       FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
c     &           I2,'/',I2,'/',I2,'  ',I2,'H',/)
c         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,0,0.0d0,
c     &        WindRefTime, MyProc, NScreen, ScreenUnit)
      ENDIF

C     jie added nws = 20: generalized asymmetric vortex model
C     2014.07 read in the geofactor, which controls the on or off
C     of the Coriolis term in the geostrophic balance

      IF(ABS(NWS).EQ.20) THEN
C...     TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &          GEOFACTOR       
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           GEOFACTOR,RSTIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           GEOFACTOR,CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           GEOFACTOR,RSTIMINC,CICE_TIMINC
         ELSE
         ENDIF
C LETS not use this now but have the option to i nthe future
c            WRITE(16,6111) IREFMO,IREFDAY,IREFYR,IREFHR
c 6112       FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
c     &           I2,'/',I2,'/',I2,'  ',I2,'H',/)
c         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,0,0.0d0,
c     &        WindRefTime, MyProc, NScreen, ScreenUnit)
      ENDIF

C     jgf49.1001 Added NWS29 format for embedding an asymmetric vortex inside
C     an OWI basin scale met field derived from NAM winds.
      IF(ABS(NWS).EQ.29) THEN
         IF(NRS.EQ.0) THEN
             READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &                  WTIMINC, pureVortex, pureBackground
         ENDIF
         IF(NRS.GE.1) THEN
             READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &                  WTIMINC, RSTIMINC, pureVortex, pureBackground
         ENDIF
      ENDIF

      !
      ! jgf50.38.05: Added capability to use hwind data.
      IF(ABS(NWS).EQ.15) THEN
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,
     &                                                 CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC, ! sb46.28sb03
     &                                                 CICE_TIMINC
      ENDIF

      ! TCM 51.06.02: Added capability to use gfdl met data.
      IF(ABS(NWS).EQ.16) THEN
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(15,*) WTIMINC
         IF((NCICE.EQ.0).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(15,*) WTIMINC,
     &                                                 CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.GE.1)) READ(15,*) WTIMINC,RSTIMINC,
     &                                                 CICE_TIMINC
      ENDIF

      IF(NWS.NE.0) WRITE(16,1117) WTIMINC
 1117 FORMAT(5X,'WIND TIME INCREMENT (SEC) = ',F10.2,/)
      IF(NRS.NE.0) WRITE(16,1118) RSTIMINC
 1118 FORMAT(5X,'RADIATION STRESS TIME INCREMENT (SEC) = ',F10.2,/)
C... TCM V49.64.01 ADDITIONS FOR ICE CONCENTRATION
      IF(NCICE.NE.0) WRITE(16,1119) CICE_TIMINC
 1119 FORMAT(5X,'ICE CONCENTRATION FIELD TIME INCREMENT (SEC) = ',
     &         F10.2,/)


C... tcm v50.66.02 -- addition for time varying bathymetry
C...
C... PROCESS BTIMINC, BCHGTIMINC if NDDT .NE. 0
C...
C   The value of NDDT,BTIMINC, and BCHGTIMINC are all
C   set in the TimeBathyControl Namelist if present

      IF (ABS(NDDT).NE.0) THEN
c...     !READ(15,*) BTIMINC,BCHGTIMINC  !read from namelist
         if (btiminc.lt.abs(dtdp)) then
            write(16,*) 'BATHYMETRY FIELD RECORD TIME LENGTH ',
     &      btiminc,'MUST BE GREATER THAN TIME STEP SIZE ',abs(dtdp) 
            CALL ADCIRC_Terminate()
         endif
         if (BCHGTIMINC.gt.btiminc) then
            write(16,*) 'BATHYMETRY TRANSITION LENGTH MUST',
     &         ' BE NO GREATER THAN BATHYMETRY FIELD RECORD',
     &         ' LENGTHS.  RESETTING.'
            BCHGTIMINC = btiminc
         endif
         if (BCHGTIMINC.lt.abs(dtdp)) then
            write(16,*) 'BATHYMETRY TRANSITION LENGTH MUST',
     &         ' BE NO SMALLER THAN THE TIME STEP SIZE.  RESETTING.'
            BCHGTIMINC = ABS(DTDP)
         ENDIF

C....... MAKE SURE THAT BATHYMETRY TRANSITION LENGTH IS AN INTEGER MULTIPLE OF THE TIMESTEP SIZE
         IF ( BCHGTIMINC/ABS(DTDP) - REAL(INT(BCHGTIMINC/ABS(DTDP))) .NE. 0.D0 ) THEN
            WRITE(16,*) 'BATHYMETRY TRANSITION LENGTH MUST',
     &         ' BE AN INTEGER MULTIPLE OF THE TIME STEP SIZE.  RESETTING.'
            BCHGTIMINC = ABS(DTDP)*REAL(INT(BCHGTIMINC/ABS(DTDP) + 0.5D0))
            DO WHILE ((BCHGTIMINC.GT.BTIMINC).AND.(BCHGTIMINC.GE.ABS(DTDP)) )
               BCHGTIMINC = ABS(DTDP)*REAL(INT(BCHGTIMINC/ABS(DTDP) + 0.5D0)-1)
            ENDDO
         ENDIF
      ELSE
         BTIMINC = 0.d0
         BCHGTIMINC = abs(dtdp)
      ENDIF
      IF (ABS(NDDT).NE.0) WRITE(16,1120) BTIMINC,BCHGTIMINC
 1120 FORMAT(5X,'BATHYMETRY FIELD TIME INCREMENT (SEC) = ',
     &         F10.2,/,
     &       5X,'BATHYMETRY TRANSITIONING LENGTH (SEC) = ',
     &         F10.2/)


C...
C...  READ AND PROCESS RNDAY - SIMULATION DURATION IN DAYS
C...
      READ(15,*) RNDAY
      WRITE(16,10) RNDAY
 10   FORMAT(5X,'TOTAL LENGTH OF NUMERICAL SIMULATION =',F12.4,
     &       5X,'DAYS',/)
C     NWS=4xx added for STWAVE tighly coupled run 01/03/2008
      IF(NRS.EQ.4) THEN
         WRITE(16,2395) NRS
 2395    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE COMPUTED ON THE FLY BY STWAVE COMPUTE    ',
     &    /,9X,'PROCS.')
      ENDIF
C...
C...  COMPUTE TOTAL NUMBER OF TIME STEPS NT
C...
#ifdef IBM
      NT=INT(RNDAY*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
#else
      NT=INT(RNDAY*(86400.D0/DTDP)+0.5d0)
#endif
      WRITE(16,1920) NT
 1920 FORMAT(5X,'NUMBER OF TIME STEPS  =',I8,/)
C...
C...  READ AND PROCESS EFFECTIVE LENGTH OF THE HYPERBOLIC TANGENT RAMP(S)
C...  IN DAYS
C...
C     jgf46.08 Add fine-grained ramp functions.
C     jgf46.21 Add FluxSettlingTime for IBTYPE=52 to accomodate
C     MS river during Katrina, split ramps for flux b.c.s into internal
C     and external.
      FluxSettlingTime = 0.0d0
      DRamp = 1.0d0
Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
      DUnRampMete=0.D0
      SELECT CASE(NRamp)
C     ---------
      CASE(0,1)! Either no ramp, or same ramp for all forcings
C     ---------
         READ(15,*) DRamp
         DRampIntFlux=DRamp
         DRampExtFlux=DRamp
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(2) ! Ramp for external flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime
         DRampIntFlux=DRamp
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(3) ! Ramp for internal flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(4) ! Ramp for surface elevation specified boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(5) ! Ramp for tidal potential
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(6) ! Ramp for wind and atmospheric pressure
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete
         DRampWRad=DRamp
C     -------
      CASE(7) ! Ramp for wave radiation stress
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete,DRampWRad
C     -------
      CASE(8) ! Start Time for Mete Ramping
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete,DRampWRad,DUnRampMete
C     ------------
      CASE DEFAULT ! fall-through
C     ------------
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NRAMP =',NRAMP
            WRITE(ScreenUnit,9713)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
#ifdef CMPI
         call msg_fini()
#endif
         STOP
C
      END SELECT
C
      IF(NRAMP.NE.0) THEN
         WRITE(16,8763) DRAMP
 8763    FORMAT(/,5X,'VALUE FOR DRAMP USED IN RAMP EVALUATION =',F12.4,
     &        5X,'DAYS',/)

         WRITE(16,5841)
 5841    FORMAT(11X,' DAYS OF SIMULATION',2X,' TIME  ',6X,'  RAMP VALUE'
     &                                                               ,/)
         IF (DRAMP.lt.1.0d-6) THEN !jgf49.44: cover the case where DRAMP is zero
            WRITE(16,*)
     &         "WARNING: DRAMP=",DRAMP,". It will be rounded to zero."
            WRITE(16,*)
     &         "All forcing and boundary conditions will be at full"
            WRITE(16,*) "strength from the start of the simulation."
         ELSE
            Day=0.0d0
 999        RampVal=TANH(Day*2.d0/DRAMP)
            WRITE(16,5845) Day,Day+StaTim,RampVal
 5845       FORMAT(15X,F8.2,6X,F8.2,2X,F15.7)
            DAY=DAY+0.5d0
            IF(Day.LT.DRAMP*1.25) GOTO 999
         ENDIF
      ENDIF
C...
C...  READ GWCE TIME WEIGHTING FACTORS
C...
      READ(15,*) A00,B00,C00
      WRITE(16,14)
14    FORMAT(//,5X,'TIME WEIGHTING FACTORS IN THE WAVE EQUATION :'/)
      WRITE(16,15) A00,B00,C00
15    FORMAT(9X,'AT TIME LEVEL K+1 : ',F8.5,
     &  /,9X,'AT TIME LEVEL K   : ',F8.5,
     &  /,9X,'AT TIME LEVEL K-1 : ',F8.5,/)
C...
C...  READ MINIMUM DEPTH OR WET/DRY PARAMETERS FROM UNIT 15
C...
      IF(NOLIFA.NE.2) THEN
         READ(15,*) H0
         WRITE(16,16) H0
 16      FORMAT(//,5X,'THE BATHYMETRIC DEPTH AT ALL NODES WILL BE ',
     &               'INCREASED TO H0= ',F12.4,' IF NECESSARY'/)
      ENDIF
      IF(NOLIFA.EQ.2) THEN
         READ(15,*) H0,NODEDRYMIN,NODEWETMIN,VELMIN
         WRITE(16,17) H0,NODEWETMIN,VELMIN,NODEDRYMIN
 17      FORMAT(//,5X,'DRYING WILL OCCUR WHEN THE WATER DEPTH < H0',
     &          /,5X,'H0 = ',F10.6,
     &          /,5X,'AND NODEREP > NODEWETMIN = ',I6,' TIME STEPS',
     &          /,5X,'NODEREP = NUMBER OF TIME STEPS SINCE A NODE ',
     &               'CHANGED STATE (EITHER WETTED OR DRIED)',
     &         //,5X,'WETTING WILL OCCUR WHEN THERE IS A FAVORABLE ',
     &               'PRESSURE GRADIENT THAT',
     &          /,5X,'WOULD DRIVE A STEADY VELOCITY TOWARDS A DRY NODE',
     &          /,5X,'THAT IS GREATER THAN VELMIN = ',F10.5,
     &          /,5X,'AND NODEREP > NODEDRYMIN = ',I6,' TIME STEPS',/)
      ENDIF
C
C     jgf46.00 Read longitude and latitude on which the CPP coordinate
C     projection is centered (in degrees) if ICS = 2. (The reading of
C     the top of the grid file, including NE and NP, was moved nearer to
C     the beginning of this subroutine.)
      READ(15,*) SLAM0,SFEA0
      WRITE(16,1112)
      WRITE(16,246)
 246  FORMAT(//,1X,'GRID INFORMATION',//)

C     jgf51.12.13: Adjust node table based on ics, nolifa, h0,
C     slam0, sfea0. Check for sufficient precision. Compute neighbor table.
      call initializeMesh()
            
C...  v49.48.02 tcm -- Allocate space for kdtree search
C...       Be sure the maximum search depth is not larger than
C...       the number of elements being kept
      IF (NE.LT.SRCHDP) SRCHDP = NE

C...  Create the search tree
      tree => kdtree2_create(bcxy,rearrange=.true.,sort=.true.)

C...  allocate space for the search results from the tree
C...  this space will be deallocated later in the subroutine
      ALLOCATE(KDRESULTS(SRCHDP))

C...
C...IF A 2DDI BAROCLINIC RUN, READ IN INITIAL CONDITION DENSITY FIELDS
C...
      IF((C2DDI).AND.(CBaroclinic)) THEN
         OPEN(11,FILE=TRIM(INPUTDIR)//'/'//'fort.11')
         READ(11,*)
         READ(11,*)
         READ(11,*) NP2
         IF(NP2.NE.NP) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9943)
            WRITE(16,9943)
 9943       FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &               //,' THE NUMBER OF NODES IN THE BAROCLINIC',
     &                  ' INITIAL CONDITION FILE (UNIT 11) ',
     &                /,' MUST EQUAL THE NUMBER OF NODES (NP) IN ',
     &                  'THE ADCIRC GRID FILE (UNIT 14)'
     &               //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
            CALL ADCIRC_Terminate()
         ENDIF

         IF    (ABS(IDEN).EQ.1) THEN
            DO I=1,NP
               READ(11,*) JKI,DASigT(JKI)
            END DO
         ELSEIF(ABS(IDEN).EQ.2) THEN
            DO I=1,NP
               READ(11,*) JKI,DASalt(JKI)
            END DO
c           CALL CALC_SIGMAT_2D()      !need to activate
         ELSEIF(ABS(IDEN).EQ.3) THEN
            DO I=1,NP
               READ(11,*) JKI,DATemp(JKI)
c           CALL CALC_SIGMAT_2D()      !need to activate
            END DO
         ELSEIF(ABS(IDEN).EQ.4) THEN
            DO I=1,NP
               READ(11,*) JKI,DATemp(JKI),DASalt(JKI)
            END DO
c           CALL CALC_SIGMAT_2D()      !need to activate
         ENDIF
         CLOSE(11)
      ENDIF

C...
C...READ IN 2DDI PASSIVE SCALAR TRANSPORT INITIAL CONDITIONS
C...
      IF((C2DDI).AND.(C2D_PTrans)) THEN
         OPEN(10,FILE=TRIM(INPUTDIR)//'/'//'fort.10')
         READ(10,*)
         READ(10,*)
         READ(10,*) NP2
         IF(NP2.NE.NP) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9944)
            WRITE(16,9943)
 9944       FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &               //,' THE NUMBER OF NODES IN THE SCALAR ',
     &                  ' INITIAL CONDITION FILE (UNIT 10) ',
     &                /,' MUST EQUAL THE NUMBER OF NODES (NP) IN ',
     &                  'THE ADCIRC GRID FILE (UNIT 14)'
     &               //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
            CALL ADCIRC_Terminate()
         ENDIF

         DO I=1,NP
            READ(10,*) JKI,DAConc(JKI)
         END DO
      ENDIF
C...
C...READ INFORMATION CONCERNING BOTTOM FRICTION COEFFICIENT
C
C     jgf46.00 If some type of spatially varying bottom friction is
C     specified in the NWP section, these inputs are ignored, and the
C     friction coefficients that are read in from the nodal attributes
C     file will take precedence.
C
C     jgf47.04 If ManningsN is loaded from the nodal attributes (fort.13)
C     file, the value of BFCdLLimit is set to CF (see nodal attributes
C     module).

      WRITE(16,1112)

      IF(NOLIBF.EQ.0) READ(15,*) TAU
      IF(NOLIBF.EQ.1) READ(15,*) CF
      IF(NOLIBF.EQ.2) READ(15,*) CF,HBREAK,FTHETA,FGAMMA
C
      WRITE(16,2045)
 2045 FORMAT(//,' BOTTOM FRICTION INFORMATION',//)
      IF(NOLIBF.EQ.0) THEN
         WRITE(16,106) TAU
 106     FORMAT(5X,'LINEAR BOTTOM FRICTION TAU =',F12.8,5X,'1/sec'/)
         IF(TAU.NE.TAU0) THEN   !CHECK TAU VALUE AGAINST TAU0
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9951)
            WRITE(16,9951)
 9951       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',
     &           //,1X,'TYPICALLY YOUR INPUT VALUE FOR ',
     &           'TAU0 SHOULD BE SET EQUAL TO TAU')
         ENDIF
      ENDIF
      IF(NOLIBF.EQ.1) THEN
         WRITE(16,8) CF
 8       FORMAT(5X,'NONLINEAR FRICTION COEFFICIENT CF =',F12.8,/)
      ENDIF
      IF(NOLIBF.EQ.2) THEN
         WRITE(16,101) CF,HBREAK,FTHETA,FGAMMA
 101     FORMAT(5X,'HYBRID FRICTION RELATIONSHIP PARAMTERS, CFMIN =',
     &        F12.8,'  HBREAK = ',F8.2,
     &        /,5X,'FTHETA = ',F8.2,'  FGAMMA = ',F10.4,//)
      ENDIF
C
C     jgf46.00 Bottom friction coefficients are initialized along with other
C     nodal attributes in the InitNodalAttr subroutine of the
C     NodalAttributes module.
C...
C...OUTPUT TO UNIT 16 GRID INFORMATION INCLUDING AGRID,NE,NP
C....H0 AND NODAL COORDINATES AND BATHYMETRY
C
C     jgf46.00 Modified this output routine so that it does not print a
C     STARTDRY column, whether the STARTDRY array has been loaded from a
C     file or not.
C
      WRITE(16,2039) AGRID
 2039 FORMAT(/,5X,'GRID IDENTIFICATION : ',A80,/)
      WRITE(16,3) NP
    3 FORMAT(5X,'TOTAL NUMBER OF NODES =',I8,/)
      WRITE(16,4) NE
    4 FORMAT(5X,'TOTAL NUMBER OF ELEMENTS =',I8,/)
      IF(ICS.EQ.2) WRITE(16,13) SLAM0*RAD2DEG,SFEA0*RAD2DEG
 13   FORMAT(5X,'LONGITUDE ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SLAM0 = ',F9.4,' DEGREES',
     &     /,5X,'LATITUDE  ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SFEA0 = ',F9.4,' DEGREES',/)
      IF(NABOUT.LT.1) THEN
         WRITE(16,24)
 24      FORMAT(/,1X,'NODAL COORDINATES AND BATHYMETRY :')
         IF(ICS.EQ.1) THEN
            IF((NTIP.EQ.0).AND.(NCOR.EQ.0)) THEN
               WRITE(16,25)
 25            FORMAT(/,10X,'NODE NO.',10X,'X',20X,'Y',15X,'DP',/)
               DO I=1,NP
                  WRITE (16,2008) I,X(I),Y(I),DP(I)
 2008             FORMAT(5X,I6,2(2X,F20.2),2X,F12.2)
               END DO
            ELSE
               WRITE(16,9195)
 9195          FORMAT(/,1X,'   NODE ',7X,'X',14X,'Y',9X,
     &              'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'DP',/)
               DO I=1,NP
                  WRITE (16,9197) I,X(I),Y(I),SLAM(I)*RAD2DEG,
     &                 SFEA(I)*RAD2DEG,DP(I)
 9197             FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2)
               END DO
            ENDIF
         ELSE
            WRITE(16,9225)
 9225       FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,
     &           'XCP',14X,'YCP',11X,'DP',/)
            DO I=1,NP
               WRITE (16,9228) I,SLAM(I)*RAD2DEG,SFEA(I)*RAD2DEG,
     &              X(I),Y(I),DP(I)
 9228          FORMAT(1X,i0,2(1X,F14.8),2(1X,F15.1),1X,F10.2)
            END DO
         ENDIF
      ELSE
         WRITE(16,3511)
 3511    FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',
     &        ' INFORMATION IS AVAILABLE IN THE',
     &        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
C...
C...OUTPUT TO UNIT 16 THE GLOBAL CONNECTIVITY TABLE (NODE NUMBERS FOR ELEMENTS)
C...
      IF(NABOUT.NE.1) THEN
         WRITE(16,26)
 26      FORMAT(//,5X,'GLOBAL NODE NUMBERS FOR EACH ELEMENT :')
         WRITE(16,27)
 27      FORMAT(/,9X,'ELEMENT',8X,'N1',9X,'N2',10X,'N3',/)
         DO I=1,NE
            WRITE(16,2009) I,NM(I,1),NM(I,2),NM(I,3)
 2009       FORMAT(8X,4(I7,4X))
         END DO
      ELSE
         WRITE(16,3512)
 3512    FORMAT(/,5X,'THE GLOBAL CONNECTIVITY TABLE',
     &        ' INFORMATION IS AVAILABLE IN THE',
     &        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
C...
C...READ IN AND WRITE OUT EDDY VISCOSITY/DIFFUSIVITY COEFFICIENTS
C...
C     jgf46.18 Made EVM and EVC nodal attributes. Their values are
C     initialized in the call to InitNodalAttr.
      IF (IM.EQ.10) THEN
        READ(15,*) ESLM,ESLC
        WRITE(16,111) ESLM,ESLC
 111    FORMAT(5X,'EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,/,
     &         5X,'EVC, EDDY DIFFUSIVITY COEFFICIENT =',E15.8,//)
      ELSE
        READ(15,*) ESLM
        IF(ESLM.LT.0.) THEN
          CSmag_Eh= .TRUE.
          ESLM=ABS(ESLM)
          WRITE(16,1111) ESLM
 1111     FORMAT(5X,'Smagorinski lateral stress coefficient with ',
     &              'constant =',E15.8,//)

          IF((CGWCE_LS_KGQ).OR.(CME_Orig)) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
              WRITE(screenunit,1112)
              WRITE(screenunit,9972)
              WRITE(screenunit,9973)
              WRITE(screenunit,1110)
              ENDIF
            WRITE(16,1112)
            WRITE(16,9972)
            WRITE(16,9973)
            WRITE(16,1110)
1110        FORMAT(5X,'The Smagorinski lateral stress coefficient is ',
     &                'not compatible with the',//
     &             5X,'Momentum Eqn formulations or with the ',
     &                'Kolar & Gray GWCE lateral stress formulation',//)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
          ENDIF
        ELSE
          WRITE(16,11) ESLM
 11       FORMAT(5X,'Constant lateral stress coefficient =',E15.8,//)
        ENDIF
      ENDIF
C...
C...  READ CORIOLIS INFORMATION AND COMPUTE THE CORIOLIS VECTOR
C...  OUTPUT RESULTING CORIOLIS INFORMATION
C...
      WRITE(16,1112)
      WRITE(16,2090)
 2090 FORMAT(//,1X,'CORIOLIS INFORMATION ',//)

      READ(15,*) CORI
      IF(NCOR.EQ.0) THEN
         DO I=1,NP
            CORIF(I)=CORI
         END DO
      ENDIF
      IF(NCOR.EQ.1) THEN
         DO I=1,NP
            CORIF(I)=2.0d0*7.29212d-5*SIN(SFEA(I))
         END DO
      ENDIF

      IF(NCOR.EQ.0) THEN
         WRITE(16,12) CORI
 12      FORMAT(5X,'CONSTANT CORIOLIS COEFFICIENT =',E15.8,5X,'1/SEC',/)
      ENDIF
      IF(NCOR.EQ.1) THEN
         WRITE(16,3604)
 3604    FORMAT(/,5X,'LATITUDES ARE USED TO COMPUTE VARIABLE CORIOLIS',
     &        /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES',/)
         IF(NABOUT.NE.1) THEN
            WRITE(16,2092)
 2092       FORMAT(/,10X,' NODE ',5X,'NODAL CORIOLIS CORIF',/)
            DO I=1,NP
               WRITE(16,2096) I,CORIF(I)
 2096          FORMAT(7X,I6,10X,E16.9)
            END DO
         ENDIF
      ENDIF
C...
C...  READ AND PROCESS INFORMATION ABOUT THE TIDAL POTENTIAL CONSTITUENTS
C...
      READ(15,*) NTIF
      mntif = ntif
      if (ntif .eq. 0) mntif = 1

C...  allocate tidal potential arrays
      call alloc_main4a()
C...  READ TIDAL POTENTIAL AMPLITUDE, FREQUENCIES, NODAL FACTORS,
C...  EQUILIBRIUM ARGUMENTS AND ALPHANUMERIC LABEL
C....
      DO I=1,NTIF
         READ(15,'(A5)')  TIPOTAG(I)
         READ(15,*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
         IF(AMIGT(I).EQ.0.) THEN
            PERT(I)=0.
         ELSE
            PERT(I)=2.D0*PI/AMIGT(I)
         ENDIF
      END DO

C...  LINES TO USE EARTH LOAD/SELF-ATTRACTION PART OF TIDAL POTENTIAL FORCING

      CALL ALLOC_MAIN4b()
      IF(NTIP.EQ.2) THEN
         OPEN(24,FILE='fort.24')
         DO I=1,NTIF
            READ(24,9930)
 9930       FORMAT(///)
            DO J=1,NP
               READ(24,*) JJ,SALTAMP(I,JJ),SALTPHA(I,JJ)
               SALTPHA(I,JJ)=SALTPHA(I,JJ)*DEG2RAD
            END DO
         END DO
      ELSE
         DO I=1,NTIF
            DO J=1,NP
               SALTAMP(I,J)=0.d0
               SALTPHA(I,J)=0.d0
            END DO
         END DO
         CLOSE(24)
      ENDIF
C...
C...  OUTPUT TO UNIT 16 INFORMATION ABOUT TIDAL POTENTIAL FORCING
C...  OUTPUT WILL VARY DEPENDING ON VALUES OF NTIP,NTIF AND NCOR
C...
      WRITE(16,1112)
      WRITE(16,2102)
 2102 FORMAT(//,1X,'TIDAL POTENTIAL FORCING INFORMATION ',//)
      WRITE(16,22) NTIF
 22   FORMAT(/,1X,'TIDAL POTENTIAL IS FORCED FOR ',I5,
     &     ' CONSTITUENT(S) ')
      IF(NTIF.GT.0) WRITE(16,23)
 23   FORMAT(/,1X,'AMPLITUDE',4X,'FREQUENCY',5X,
     &     '    ETRF      ','NODAL FACTOR',2X,
     &     'EQU.ARG(DEG)',1X,'CONSTITUENT',/)
      DO I=1,NTIF
         WRITE(16,2107) TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I),
     &        TIPOTAG(I)
 2107    FORMAT(1X,F10.7,1X,F15.12,2X,F10.7,5X,F10.7,1X,F10.3,7X,A5)
      END DO
C...
C...  CONVERT FACET(I) VALUES FROM DEGREES TO RADIANS
C...
      DO I=1,NTIF
         FACET(I)=FACET(I)*DEG2RAD
      END DO
C...
C...  CHECK CONSISTENCY OF INPUT PARAMETERS NTIF AND NTIP
C...
      IF(((NTIP.EQ.0).AND.(NTIF.NE.0)).OR.((NTIP.NE.0).AND.
     &     (NTIF.EQ.0))) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9961)
         WRITE(16,9961)
 9961    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NTIF AND NTIP (UNIT 15 INPUT ',
     &        'PARAMETERS) IS INCONSISTENT',
     &        /,1X,'PLEASE CHECK THESE VALUES')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9987)
            WRITE(16,9987)
 9987       FORMAT(/,1X,'PROGRAM WILL OVERRIDE THE SPECIFIED ',
     &           'INPUT AND NEGLECT TIDAL POTENTIAL TERMS',
     &           /,1X,' AND/OR RESET NTIP = 0',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NTIP=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
            CALL ADCIRC_Terminate()
         ENDIF
         GOTO 1893
      ENDIF
C...
C...  PRINT OUT LAT/LON VALUES TO BE USED IN COMPUTING TIDAL POTENTIAL
C...  IF NOT ALREADY DONE SO IN CORIOLIS SECTION AND TIDAL POTENTIAL IS
C...  ACTIVATED WITH NTIP=1
C...
      IF(NTIP.GE.1) THEN
         IF(ICS.EQ.1) THEN
            WRITE(16,3605)
 3605      FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &             ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &             /,7X,'AND ARE BASED ON AN INVERSE CPP PROJECTION ',
     &             'OF THE INPUT COORDINATES',/)
         ELSE
            WRITE(16,2109)
 2109       FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &             ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &             /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES ',/)
         ENDIF
      ENDIF
C...
C...  INPUT FROM UNIT 15 THE TIDAL FORCING FREQUENCIES ON THE ELEVATION
C...  SPECIFIED BOUNDARIES: INCLUDING NBFR, FREQUENCIES, NODAL FACTORS,
C...  EQUILIBRIUM ARGUMENTS AND AN ELEVATION BOUNDARY CONDITION
C...  ALPHANUMERIC DESCRIPTOR
C...
1893  READ(15,*) NBFR
      MNBFR = NBFR
      IF (NBFR.EQ.0) MNBFR = 1

C     Allocate arrays dimensioned by MNBFR
      call alloc_main5()

      WRITE(16,1112)
      WRITE(16,2106)
 2106 FORMAT(//,1X,'ELEVATION SPECIFIED BOUNDARY FORCING INFORMATION '
     &       ,//)
      WRITE(16,20) NBFR
 20   FORMAT(/,5X,'NUMBER OF PERIODIC, ELEVATION SPECIFIED ',
     &       'CONSTITUENTS =',I5)
      IF(NBFR.GE.1) WRITE(16,21)
 21   FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &       3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
      DO I=1,NBFR
         READ(15,'(A5)') BOUNTAG(I)
         READ(15,*) AMIG(I),FF(I),FACE(I)
         WRITE(16,1850) I,AMIG(I),FF(I),FACE(I),BOUNTAG(I)
 1850    FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
         FACE(I)=FACE(I)*DEG2RAD
         IF(AMIG(I).EQ.0.) THEN
            PER(I)=0.
         ELSE
            PER(I)=2.D0*PI/AMIG(I)
         ENDIF
      END DO

C...  INPUT FORCING CONDITIONS ON PERIODIC ELEVATION SPECIFIED
C...  BOUNDARIES FOR EACH OF THE ELEVATION FORCING FREQUENCIES FROM UNIT
C...  15
C...
      ALLOCATE(ELEVALPHA(NBFR))
      DO I=1,NBFR
         READ(15,'(A10)') ELEVALPHA(i)
         DO J=1,NETA
            READ(15,*) EMO(I,J),EFA(I,J)
         END DO
      END DO

C.....READ THE MINIMUM INNER ANGLE FOR WHICH VELOCITY AT FLOW BOUNDARY NODES
C.....WILL BE ZEROED IN THE TANGENTIAL DIRECTIONS WHEN NORMAL FLOW IS AN
C.....ESSENTIAL B.C.

      READ(15,*) ANGINN
      WRITE(16,1112)
      WRITE(16,7654) ANGINN
7654  FORMAT(//,5X,'ANGINN = ',F8.2,' DEGREES',
     &  /,5X,'ALL FLOW BOUNDARY NODES WITH NORMAL FLOW AS AN ',
     &       'ESSENTIAL B.C. AND ',
     &  /,9X,'INNER ANGLES LESS THAN ANGINN WILL HAVE BOTH NORMAL ',
     &  /,9X,'AND TANGENTIAL VELOCITY COMPONENTS ZEROED',/)
      COSTSET=COS(ANGINN*DEG2RAD)
      ! 
      ! jgf51.21.12: Now that the value of anginn has been read, 
      ! the boundaries can be checked and the boundary arrays can
      ! be constructed. This initialization also determines if there
      ! are flux boundaries in the mesh, and sets the value of 
      ! nfluxf accordingly. The nfluxf value is required for further
      ! parsing of the fort.15 control file; it determines whether
      ! the value of NFFR should be read below. 
      call initializeBoundaries()
      !
      ! jgf51.21.12: Now write out log messages that include NBD (the
      ! NBD array was also computed in the call to initializeBoundaries().
      DO I=1,NBFR
         WRITE(16,29) I,BOUNTAG(I)
 29      FORMAT(////,5X,'ELEVATION BOUNDARY TIDAL FORCING FOR',
     &        ' CONSTITUENT NUMBER',i0,1X,'DESIGNATED : ',A5)
         WRITE(16,31) ELEVALPHA(i)
 31      FORMAT(9X,'VERIFICATION OF CONSTITUENT : ',A10,/)
         WRITE(16,30)
 30      FORMAT(14X,'NODE',11X,'AMPL.',9X,'PHASE(DEG)',/)
         DO J=1,NETA
            WRITE(16,1870) NBD(J),EMO(I,J),EFA(I,J)
 1870       FORMAT(10X,I8,4X,F14.5,4X,F12.3)
            EFA(I,J)=EFA(I,J)*DEG2RAD
         END DO
      END DO

   
C...IF ANY NON ZERO NORMAL FLOW BOUNDARIES WERE SPECIFIED, (NFLUXF=1)
C.....READ FORCING INFORMATION FROM UNIT 15 FILE

      NFFR = 0
      IF(NFLUXF.EQ.1) THEN

C.....INPUT FROM THE NUMBER OF FREQUENCIES PRESENT IN NORMAL FLOW FORCING
C......DATA.  IF THIS = 0, NORMAL FLOW DATA IS READ IN FROM THE FORT.20 FILE.
Ckmd  Made it to where if NFFR = -1, then the normal flow data is to be read
C     from the fort.20 file, too. However, it will not be read from time=0 but
C     from time=hotstart time

         READ(15,*) NFFR
         MNFFR = NFFR
         IF (NFFR.EQ.0) MNFFR = 1
         IF (NFFR.EQ.-1) MNFFR = 1

C.....Allocate space for periodic normal flow boundary conditions
         call alloc_main6()
C
         DO I=1,NVELME
            J=ME2GW(I)
            QNAM(1,J)=0.
            QNPH(1,J)=0.
         END DO

C.....READ IN AND WRITE OUT INFO ON SPECIFIED NORMAL FLOW BOUNDARIES
Ckmd  Added in a NFFR=-1 for reading information from the fort.20 files
C     for hot starting the run, so file does not have to include all
C     river data from time=0.

         WRITE(16,1112)
         WRITE(16,2200)
 2200   FORMAT(//,1X,'NORMAL FLOW BOUNDARY FORCING INFORMATION ',//)
         IF((NFFR.EQ.0).OR.(NFFR.EQ.-1)) THEN
            WRITE(16,2201)
 2201      FORMAT(/,5X,'NORMAL FLOW VALUES WILL BE READ FROM UNIT 20 ',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE FLOW DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
         ENDIF
         IF((NFFR.NE.0).AND.(NFFR.NE.-1)) THEN
            WRITE(16,2202) NFFR
 2202      FORMAT(/,5X,'NUMBER OF PERIODIC NORMAL FLOW CONSTITUENTS =',
     &                                                               I5)
            WRITE(16,2203)
 2203      FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &         3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
            DO I=1,NFFR
               READ(15,'(A5)') FBOUNTAG(I)
               READ(15,*) FAMIG(I),FFF(I),FFACE(I)
               WRITE(16,2204) I,FAMIG(I),FFF(I),FFACE(I),FBOUNTAG(I)
 2204          FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
               FFACE(I)=FFACE(I)*DEG2RAD
               IF(FAMIG(I).EQ.0.) THEN
                  FPER(I)=0.
               ELSE
                  FPER(I)=2.D0*PI/FAMIG(I)
               ENDIF
            END DO

C.......INPUT PERIODIC NORMAL FLOW FORCING CONDITIONS ON DESIGNATED FLOW BOUNDARIES
C........FOR EACH OF THE FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16
            DO I=1,NFFR
               WRITE(16,2206) I,FBOUNTAG(I)
 2206          FORMAT(////,5X,'PERIODIC NORMAL FLOW CONSTITUENT ',
     &                     'NUMBER',I4,1X,'DESIGNATED : ',A5)
               READ(15,'(A10)') ALPHA
               WRITE(16,31) ALPHA
               WRITE(16,30)
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12).OR.
     &                 (LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.52)) THEN
                     READ(15,*) QNAM(I,J),QNPH(I,J)
                     WRITE(16,2205) NBV(J),QNAM(I,J),QNPH(I,J)
 2205                FORMAT(10X,I8,4X,F14.5,4X,F12.3)
                     QNPH(I,J)=QNPH(I,J)*DEG2RAD
                  ENDIF
                  IF(LBCODEI(J).EQ.32) THEN
                     READ(15,*) QNAM(I,J),QNPH(I,J),ENAM(I,J),ENPH(I,J)
                     WRITE(16,2207) NBV(J),QNAM(I,J),QNPH(I,J),
     &                  ENAM(I,J),ENPH(I,J)
 2207                FORMAT(10X,I8,4X,F14.5,4X,F12.3,4X,F14.5,4X,F12.3)
                     QNPH(I,J)=QNPH(I,J)*DEG2RAD
                     ENPH(I,J)=ENPH(I,J)*DEG2RAD
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
C...
C...READ IN INFORMATION CONCERNING OUTPUT REQUIREMENTS FROM UNIT 15 AND
C...OUTPUT THIS TO UNIT 16
C...
      WRITE(16,1112)
      WRITE(16,3000)
3000  FORMAT(//,1X,'OUTPUT INFORMATION WILL BE PROVIDED AS'
     &  ,' FOLLOWS :')

C...
C...INPUT INFORMATION FOR ELEVATION RECORDING STATIONS
C...

C....READ IN NOUTE,TOUTSE,TOUTFE,NSPOOLE : IF ABS(NOUTE)>0, INTERPOLATED
C....ELEVATIONS AT ELEVATION STATIONS ARE SPOOLED TO UNIT 61 EVERY NSPOOLE
C....TIME STEPS BETWEEN TIMES TOUTSE AND TOUTFE

      READ(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      WRITE(16,3001) NOUTE
 3001 FORMAT(///,1X,'ELEVATION RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTE = ',I2)

C....CHECK INPUT PARAMETER NOUTE
      SELECT CASE(ABS(NOUTE))
         CASE(0)
            ! IF STATION ELEVATION OUTPUT WILL NOT BE GENERATED
            CALL logMessage(INFO,
     &     'NO OUTPUT WILL BE SPOOLED AT ELEVATION RECORDING STATIONS.')
         CASE(1)
            CALL logMessage(INFO,'UNIT 61 FORMAT WILL BE ASCII.')
         CASE(2)
            CALL logMessage(INFO,'UNIT 61 FORMAT WILL BE BINARY.')
         CASE(3)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &       'UNIT 61 WILL BE NETCDF CLASSIC MODEL / NETCDF3 FORMAT.')
         CASE(5)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &  'UNIT 61 WILL BE NETCDF CLASSIC MODEL / NETCDF4 (HDF5) FORMAT.')
         CASE(4,6:)
            call allMessage(ERROR,"This NOUTE value is invalid.")
            call ADCIRC_Terminate()
         CASE DEFAULT
            ! do nothing, the other cases handled below
      END SELECT

C....IF STATION ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTE.NE.0) THEN

C......COMPUTE NTCYSE, NTCYFE, WHICH = TOUTSE AND TOUTFE IN TIMESTEPS

#ifdef IBM
         NTCYSE=INT((TOUTSE-STATIM)*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
         NTCYFE=INT((TOUTFE-STATIM)*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
#else
         NTCYSE=INT((TOUTSE-STATIM)*(86400.D0/DTDP)+0.5d0)
         NTCYFE=INT((TOUTFE-STATIM)*(86400.D0/DTDP)+0.5d0)
#endif
         IF(NTCYFE.GT.NT) NTCYFE=NT

C......COMPUTE NTRSPE = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 61

         IF(NSPOOLE.EQ.0) NTRSPE=0
#ifdef IBM
         IF(NSPOOLE.NE.0) NTRSPE=INT((NTCYFE-NTCYSE)/NSPOOLE,KIND(0.0d0))
#else
         IF(NSPOOLE.NE.0) NTRSPE=INT((NTCYFE-NTCYSE)/NSPOOLE)
#endif
C......WRITE TOUTSE,TOUTFE,NTCYSE,NTCYFE,NSPOOLE TO UNIT 16

         WRITE(16,3004) TOUTSE,NTCYSE,TOUTFE,NTCYFE,NSPOOLE
 3004    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSE =',F8.3,
     &              ' DAY(S) RELATIVE',
     &         /,9X,'TO THE STARTING TIME OR',I9,
     &              ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 61 EVERY',
     &              ' NSPOOLE =',I8,' TIME STEPS')

      ENDIF
C....REGARDLESS OF WHETHER NOUTE=0, READ IN THE NUMBER OF ELEVATION
C....RECORDING STATIONS

C.... tcm v51.20.04 - additions for external specification of 
C       elevation station locations
      READ(15,*) NSTAE
      STAT_LUN = 15
      IF (NSTAE < 0) THEN
         USE_ELEV_STAT_FILE = .TRUE.
         IOS_STATIONS = 0
         STAT_LUN = 151
         NSTAE = ABS(NSTAE) !SET TO POSITIVE
         WRITE(16,*) '    ELEVATION RECORDING STATIONS WILL BE READ FROM',
     &               ' AN EXTERNAL FILE.'
         NSTAE2 = 0
         OPEN(unit=stat_lun,file=TRIM(INPUTDIR)//'/'//'elev_stat.151',
     &        status='old',err=7690,iostat=ios_stations)
         READ(151,*) NSTAE2
         IF (ABS(NSTAE2).NE. ABS(NSTAE)) THEN
            NSTAE = ABS(NSTAE2)  !RESET THE VALUE TO WHAT'S IN THE FILE
         ENDIF
         WRITE(16,3007) NSTAE
7690     IF (IOS_STATIONS .NE. 0) THEN
            WRITE(16,*) "ERROR IN READING ELEVATION STATION FILE: elev_stat.151"
            WRITE(16,*) " Stopping Execution"
            call allMessage(ERROR,"Problem Reading Elevation Station File.")
            call ADCIRC_Terminate()
#ifdef CMPI
            call msg_fini()
#endif
            stop  ! there is a stop here
         ENDIF                  
      ELSE         
         WRITE(16,3007) NSTAE
         IF (NSTAE.NE.0) THEN
            WRITE(16,*) "ELEVATION STATION LOCATIONS WILL BE READ FROM FORT.15"
         ENDIF
      ENDIF
 3007 FORMAT(///,5X,'NUMBER OF INPUT ELEVATION RECORDING STATIONS = ',
     &              I5)
C
C
      IF(NSTAE.GT.0) THEN
         IF(ICS.EQ.1) WRITE(16,3008)
 3008     FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
         IF(ICS.EQ.2) WRITE(16,3009)
 3009      FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         MNSTAE = NSTAE
      ENDIF
      IF (NSTAE.EQ.0) MNSTAE = 1

C  Allocate arrays dimensioned by MNSTAE
      call alloc_main7()
C
C....INPUT COORDINATES OF ELEVATION RECORDING STATIONS THEN COMPUTE
C....THE ELEMENT NO. THE STATION LIES IN
      ALLOCATE(STATNAME(MNSTAE))
      CALL readStations(STATNAME, NSTAE, NNE, XEL, YEL, SLEL, SFEL,
     &                  STAIE1, STAIE2, STAIE3,STAT_LUN,
     &                  'ELEVATION RECORDING STATION   ')
      !tcm v51.20.04 addition for external station file
      IF ((USE_ELEV_STAT_FILE).AND.(STAT_LUN ==151)) CLOSE(STAT_LUN)
 9911 FORMAT(F12.3,2X,F12.3, 6X, A50)
 9111 FORMAT(F8.3,2X,F8.3, 6X, A50)
C...
C...INPUT INFORMATION FOR VELOCITY RECORDING STATIONS
C...

C....READ IN NOUTV,TOUTSV,TOUTFV,NSPOOLV : IF NOUTV<>0,INTERPOLATED VELOCITIES AT
C....VELOCITY STATIONS ARE SPOOLED TO UNIT 62 EVERY NSPOOLV TIME STEPS BETWEEN
C....TIMES TOUTSV AND TOUTFV; IF ABS(NOUTV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      WRITE(16,3101) NOUTV
 3101 FORMAT(////,1X,'VELOCITY RECORDING STATION OUTPUT : ',
     &         //,5X,'NOUTV = ',I2)

C....CHECK INPUT PARAMETER NOUTV
      SELECT CASE(ABS(NOUTV))
         CASE(0)
            ! IF STATION OUTPUT WILL NOT BE GENERATED
            CALL logMessage(INFO,
     &      'NO OUTPUT WILL BE SPOOLED AT VELOCITY RECORDING STATIONS')
         CASE(1)
            CALL logMessage(INFO,'UNIT 62 FORMAT WILL BE ASCII.')
         CASE(2)
            CALL logMessage(INFO,'UNIT 62 FORMAT WILL BE BINARY.')
         CASE(3)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &       'UNIT 62 WILL BE NETCDF CLASSIC MODEL / NETCDF3 FORMAT.')
         CASE(5)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &  'UNIT 62 WILL BE NETCDF CLASSIC MODEL / NETCDF4 (HDF5) FORMAT.')
         CASE(4,6:)
            call allMessage(ERROR,"This NOUTV value is invalid.")
            call ADCIRC_Terminate()
         CASE DEFAULT
            ! do nothing, the other cases handled below
      END SELECT

C....IF STATION VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTV.NE.0) THEN

C......  COMPUTE NTCYSV, NTCYFV, WHICH = TOUTSV AND TOUTFV IN TIME STEPS
#ifdef IBM
         NTCYSV=INT((TOUTSV-STATIM)*(86400.D0/DTDP) + 0.5d0,KIND(0.0d0))
         NTCYFV=INT((TOUTFV-STATIM)*(86400.D0/DTDP) + 0.5d0,KIND(0.0d0))
#else
         NTCYSV=INT((TOUTSV-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFV=INT((TOUTFV-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
         IF(NTCYFV.GT.NT) NTCYFV=NT

C......CALCULATE NTRSPV = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 62

         IF(NSPOOLV.EQ.0) NTRSPV=0
#ifdef IBM
         IF(NSPOOLV.NE.0) NTRSPV=INT((NTCYFV-NTCYSV)/NSPOOLV,KIND(0.0d0))
#else
         IF(NSPOOLV.NE.0) NTRSPV=INT((NTCYFV-NTCYSV)/NSPOOLV)
#endif

C......WRITE NOUTV,TOUTSV,TOUTFV,NTCYSV,NTCYFV,NSPOOLV TO UNIT 16

         WRITE(16,3104) TOUTSV,NTCYSV,TOUTFV,NTCYFV,NSPOOLV
 3104    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSV =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFV =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 62 EVERY ',
     &              ' NSPOOLV =',I8,' TIME STEPS')

      ENDIF
C....REGARDLESS OF WHETHER NOUTV=0, READ IN THE NUMBER OF VELOCITY
C....RECORDING STATIONS

      READ(15,*) NSTAV
      STAT_LUN = 15
      IF (NSTAV < 0) THEN
         USE_VEL_STAT_FILE = .TRUE.
         ios_stations = 0
         stat_lun = 151
         NSTAV = ABS(NSTAV) !SET TO POSITIVE
         WRITE(16,*) '   VELOCITY RECORDING STATIONS WILL BE READ FROM',
     &               ' AN EXTERNAL FILE.'
         NSTAV2 = 0
         OPEN(unit=stat_lun,file=TRIM(INPUTDIR)//'/'//'vel_stat.151',
     &        status='old', err=7691,iostat=ios_stations)
         READ(151,*) NSTAV2
         IF (ABS(NSTAV2).NE. ABS(NSTAV)) THEN
            NSTAV = ABS(NSTAV2)  !RESET THE VALUE TO WHAT'S IN THE FILE
         ENDIF
         WRITE(16,3107) NSTAV
7691     IF (IOS_STATIONS .NE. 0) THEN
            WRITE(16,*) "ERROR IN READING VELOCITY STATION FILE: vel_stat.151"
            WRITE(16,*) " STOPPING EXECUTION"
            call allMessage(ERROR,"Problem Reading Velocity Station File.")
            call ADCIRC_Terminate()
#ifdef CMPI
            call msg_fini()
#endif
            stop  ! there is a stop here
         ENDIF                  
      ELSE         
         WRITE(16,3107) NSTAV
         IF (NSTAV.NE.0) THEN
            WRITE(16,*) "VELOCITY STATION LOCATIONS WILL BE READ FROM FORT.15"
         ENDIF
      ENDIF
 3107 FORMAT(////,5X,'NUMBER OF INPUT VELOCITY RECORDING STATIONS = ',
     &            I5)

      IF(NSTAV.GT.0) THEN
         IF(ICS.EQ.1) WRITE(16,3108)
 3108     FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
         IF(ICS.EQ.2) WRITE(16,3109)
 3109     FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &             4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         MNSTAV = NSTAV
      ENDIF
      IF (NSTAV.EQ.0) MNSTAV = 1

C      Allocate arrays dimensioned by MNSTAV
      call alloc_main8()

C....INPUT COORDINATES OF VELOCITY RECORDING STATIONS
C....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES
      ALLOCATE(STATNAMEV(MNSTAV))
      CALL readStations(STATNAMEV, NSTAV, NNV, XEV, YEV, SLEV, SFEV,
     &                  STAIV1, STAIV2, STAIV3,STAT_LUN,
     &                 'VELOCITY RECORDING STATION    ' )
      !tcm v51.20.04 addition for external station file
      IF ((USE_VEL_STAT_FILE).AND.(STAT_LUN ==151)) CLOSE(STAT_LUN)

C...
C...
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION FOR CONCENTRATION
C...  RECORDING STATIONS
C...
      NOUTC=0
      IF(IM.EQ.10) THEN

C...  READ IN NOUTC,TOUTSC,TOUTFC,NSPOOLC : IF NOUTC<>0,INTERPOLATED
C...  CONCENTRATIONS ARE SPOOLED TO UNIT 81 EVERY NSPOOLC TIME STEPS
C...  BETWEEN TIMES TOUTSC AND TOUTFC; IF ABS(NOUTC)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
         WRITE(16,3201) NOUTC
 3201    FORMAT(///,1X,'CONCENTRATION RECORDING STATION OUTPUT : ',
     &          //,5X,'NOUTC = ',I2)

C...     CHECK INPUT PARAMETER NOUTC
         SELECT CASE(ABS(NOUTC))
            CASE(0)
               ! IF STATION OUTPUT WILL NOT BE GENERATED
               CALL logMessage(INFO,
     &         'NO OUTPUT WILL BE SPOOLED AT CONC. RECORDING STATIONS')
            CASE(1)
               CALL logMessage(INFO,'UNIT 81 FORMAT WILL BE ASCII.')
            CASE(2)
               CALL logMessage(INFO,'UNIT 81 FORMAT WILL BE BINARY.')
            CASE(3)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &       'UNIT 81 WILL BE NETCDF CLASSIC MODEL / NETCDF3 FORMAT.')
            CASE(5)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &  'UNIT 81 WILL BE NETCDF CLASSIC MODEL / NETCDF4 (HDF5) FORMAT.')
            CASE(4,6:)
               call allMessage(ERROR,"This NOUTC value is invalid.")
               call ADCIRC_Terminate()
            CASE DEFAULT
               ! do nothing, the other cases handled below
         END SELECT

C...  IF STATION CONCENTRATION OUTPUT WILL BE GENERATED

         NSTAC = 0
         IF(NOUTC.NE.0) THEN

C...  COMPUTE NTCYSC, NTCYFC, WHICH = TOUTSC AND TOUTFC IN TIMESTEPS
#ifdef IBM
            NTCYSC=INT((TOUTSC-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
            NTCYFC=INT((TOUTFC-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
#else
            NTCYSC=INT((TOUTSC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFC=INT((TOUTFC-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
            IF(NTCYFC.GT.NT) NTCYFC=NT

C...  COMPUTE NTRSPC = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 81

            IF(NSPOOLC.EQ.0) NTRSPC=0
#ifdef IBM
            IF(NSPOOLC.NE.0) NTRSPC=INT((NTCYFC-NTCYSC)/NSPOOLC,
     &           KIND(0.0d0))
#else
            IF(NSPOOLC.NE.0) NTRSPC=INT((NTCYFC-NTCYSC)/NSPOOLC)
#endif

C...  WRITE TOUTSC,TOUTFC,NTCYSC,NTCYFC,NSPOOLC TO UNIT 16

            WRITE(16,3204) TOUTSC,NTCYSC,TOUTFC,NTCYFC,NSPOOLC
 3204       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 81 EVERY',
     &                ' NSPOOLC =',I8,' TIME STEPS')
         ENDIF

C...  REGARDLESS OF WHETHER NOUTC=0, READ IN THE NUMBER OF CONCENTRATION
C...  RECORDING STATIONS

         READ(15,*) NSTAC
         STAT_LUN = 15
         IF (NSTAC < 0) THEN
            USE_CONC_STAT_FILE = .TRUE.
            IOS_STATIONS = 0
            STAT_LUN = 151
            NSTAC = ABS(NSTAC) !SET TO POSITIVE
            WRITE(16,*) 'CONCENTRATION RECORDING STATIONS WILL BE READ FROM',
     &               ' AN EXTERNAL FILE.'
            NSTAC2 = 0
            OPEN(unit=stat_lun,file=TRIM(INPUTDIR)//'/'//'conc_stat.151',
     &           status='old',err=7692,iostat=ios_stations)
            READ(151,*) NSTAC2
            IF (ABS(NSTAC2).NE. ABS(NSTAC)) THEN
               NSTAC = ABS(NSTAC2)  !RESET THE VALUE TO WHAT'S IN THE FILE
            ENDIF
            WRITE(16,3207) NSTAC
7692        IF (IOS_STATIONS .NE. 0) THEN
               WRITE(16,*) "ERROR IN READING CONCENTRATION STATION FILE: conc_stat.151"
               WRITE(16,*) " Stopping Execution"
               call allMessage(ERROR,"Problem Reading Concentration Station File.")
               call ADCIRC_Terminate()
#ifdef CMPI
               call msg_fini()
#endif
               stop  ! there is a stop here
            ENDIF                  
         ELSE         
            WRITE(16,3207) NSTAC
            IF (NSTAC.NE.0) THEN
               WRITE(16,*) "CONC. STATION LOCATIONS WILL BE READ FROM FORT.15"
            ENDIF
         ENDIF
 3207    FORMAT(///,5X,'NUMBER OF INPUT CONCENTRATION RECORDING ',
     &        'STATIONS = ',I5)

         IF(NSTAC.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3208)
 3208       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3209)
 3209       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAC = NSTAC
         ENDIF
         IF (NSTAC.EQ.0) MNSTAC = 1

C  Allocate arrays dimensioned by MNSTAC
         call alloc_main9()


C...  INPUT COORDINATES OF CONCENTRATION RECORDING STATIONS
C...  THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES
         ALLOCATE(STATNAMEC(NSTAC))
         CALL readStations(STATNAMEC, NSTAC, NNC, XEC, YEC, SLEC, SFEC,
     &                  STAIC1, STAIC2, STAIC3,STAT_LUN,
     &                  'CONCENTRATION REC. STATION    ')
         !tcm v51.20.04 addition for external station file
         IF ((USE_CONC_STAT_FILE).AND.(STAT_LUN ==151)) CLOSE(STAT_LUN)
      ENDIF

C...  IF METEOROLOICAL FORCING IS INCLUDED IN THE RUN, INPUT
C...  INFORMATION FOR MET RECORDING STATIONS - OUTPUT
C...
      NOUTM=0
      NSTAM=0
C
      IF(NWS.NE.0) THEN

C...  READ IN NOUTM,TOUTSM,TOUTFM,NSPOOLM : IF NOUTM<>0,INTERPOLATED
C...  MET DATA ARE SPOOLED TO UNITS 71&72 EVERY NSPOOLM TIME STEPS
C...  BETWEEN TIMES TOUTSM AND TOUTFM; IF ABS(NOUTM)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         WRITE(16,3211) NOUTM
 3211    FORMAT(///,1X,'METEOROLOGICAL RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTM = ',I2)

C...     CHECK INPUT PARAMETER NOUTM
         SELECT CASE(ABS(NOUTM))
            CASE(0)
               ! IF STATION OUTPUT WILL NOT BE GENERATED
               CALL logMessage(INFO,
     &         'NO OUTPUT WILL BE SPOOLED AT MET. RECORDING STATIONS')
            CASE(1)
               CALL logMessage(INFO,
     &         'UNIT 71 AND 72 FORMATS WILL BE ASCII.')
            CASE(2)
               CALL logMessage(INFO,
     &          'UNIT 71 AND 72 FORMATS WILL BE BINARY.')
            CASE(3)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &          'UNIT 71 AND 72 FORMATS WILL BE NETCDF CLASSIC MODEL'
     &           //' / NETCDF3 FORMAT.')
            CASE(5)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &         'UNIT 71 AND 72 FORMATS WILL BE NETCDF CLASSIC MODEL'
     &         //' / NETCDF4 (HDF5) FORMAT.')
            CASE(4,6:)
               call allMessage(ERROR,"This NOUTM value is invalid.")
               call ADCIRC_Terminate()
            CASE DEFAULT
               ! do nothing, the other cases handled below
         END SELECT

C...  IF STATION MET OUTPUT WILL BE GENERATED

         IF(NOUTM.NE.0) THEN

C...  COMPUTE NTCYSM, NTCYFM, WHICH = TOUTSM AND TOUTFM IN TIMESTEPS
#ifdef IBM
            NTCYSM=INT((TOUTSM-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
            NTCYFM=INT((TOUTFM-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
#else
            NTCYSM=INT((TOUTSM-STATIM)*(86400.D0/DTDP))
            NTCYFM=INT((TOUTFM-STATIM)*(86400.D0/DTDP))
#endif

            IF(NTCYFM.GT.NT) NTCYFM=NT

C...  COMPUTE NTRSPM = THE NO. OF DATA SETS TO BE SPOOLED TO UNITS 71&72

            IF(NSPOOLM.EQ.0) NTRSPM=0
#ifdef IBM
            IF(NSPOOLM.NE.0) NTRSPM=INT((NTCYFM-NTCYSM)/NSPOOLM,
     &           KIND(0.0d0))
#else
            IF(NSPOOLM.NE.0) NTRSPM=INT((NTCYFM-NTCYSM)/NSPOOLM)
#endif

C...  WRITE TOUTSM,TOUTFM,NTCYSM,NTCYFM,NSPOOLM TO UNIT 16

            WRITE(16,3214) TOUTSM,NTCYSM,TOUTFM,NTCYFM,NSPOOLM
 3214       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNITS 71&72',
     &           ' EVERY NSPOOLM =',I8,' TIME STEPS')
         ENDIF

C...  REGARDLESS OF WHETHER NOUTM=0, READ IN THE NUMBER OF METEOROLOGICAL
C...  RECORDING STATIONS

         READ(15,*) NSTAM
         STAT_LUN = 15
         IF (NSTAM < 0) THEN
            USE_MET_STAT_FILE = .TRUE.
            IOS_STATIONS = 0
            STAT_LUN = 151
            NSTAM = ABS(NSTAM) !SET TO POSITIVE
            WRITE(16,*) 'MET RECORDING STATIONS WILL BE READ FROM',
     &               ' AN EXTERNAL FILE.'
            NSTAM2 = 0
            OPEN(unit=stat_lun,file=TRIM(INPUTDIR)//'/'//'met_stat.151',
     &       status='old',err=7693,iostat=ios_stations)
            READ(151,*) NSTAM2
            IF (ABS(NSTAM2).NE. ABS(NSTAM)) THEN
               NSTAM = ABS(NSTAM2)  !RESET THE VALUE TO WHAT'S IN THE FILE
            ENDIF
            WRITE(16,3217) NSTAM
7693        IF (IOS_STATIONS .NE. 0) THEN
               WRITE(16,*) "ERROR IN READING MET STATION FILE: elev_stat.151"
               WRITE(16,*) " Stopping Execution"
               call allMessage(ERROR,"Problem Reading Met Station File.")
               call ADCIRC_Terminate()
#ifdef CMPI
               call msg_fini()
#endif
               stop  ! there is a stop here
            ENDIF                  
         ELSE         
            WRITE(16,3217) NSTAM
            IF (NSTAM.NE.0) THEN
               WRITE(16,*) "MET STATION LOCATIONS WILL BE READ FROM FORT.15"
            ENDIF
         ENDIF
 3217    FORMAT(///,5X,'NUMBER OF INPUT METEOROLOGICAL RECORDING ',
     &        'STATIONS = ',I5)

         IF(NSTAM.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3218)
 3218       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3219)
 3219       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAM = NSTAM
         ENDIF
         IF (NSTAM.EQ.0) MNSTAM = 1

C  Allocate arrays dimensioned by MNSTAM
         call alloc_main10()

C...  INPUT COORDINATES OF METEOROLOGICAL RECORDING STATIONS
C...  THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES
         ALLOCATE(STATNAMEM(MNSTAM))
         CALL readStations(STATNAMEM, NSTAM, NNM, XEM, YEM, SLEM, SFEM,
     &                     STAIM1, STAIM2, STAIM3,STAT_LUN,
     &                      'METEOROLOGICAL REC. STATION   ')
         !tcm v51.20.04 addition for external station file
         IF ((USE_MET_STAT_FILE).AND.(STAT_LUN ==151)) CLOSE(STAT_LUN)
      ENDIF

C...
C...  INPUT INFORMATION ABOUT GLOBAL ELEVATION DATA OUTPUT
C...

C...  READ IN NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : IF NOUTGE<>0, GLOBAL ELEV.
C...  OUTPUT IS SPOOLED TO UNIT 63 EVERY NSPOOLGE TIME STEPS BETWEEN
C...  TIMES TOUTSGE AND TOUTFGE; IF ABS(NOUTGE)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      WRITE(16,3301) NOUTGE
 3301 FORMAT(////,1X,'GLOBAL NODAL ELEVATION INFORMATION OUTPUT: ',
     &     //,5X,'NOUTGE = ',I2)

C...  CHECK INPUT PARAMETER NOUTGE
      SELECT CASE(ABS(NOUTGE))
         CASE(0)
            ! IF STATION OUTPUT WILL NOT BE GENERATED
            CALL logMessage(INFO,
     &         'NO GLOBAL ELEVATION OUTPUT WILL BE SPOOLED.')
         CASE(1)
            CALL logMessage(INFO,'UNIT 63 FORMAT WILL BE ASCII.')
         CASE(2)
            CALL logMessage(INFO,'UNIT 63 FORMAT WILL BE BINARY.')
         CASE(3)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &          'UNIT 63 FORMAT WILL BE NETCDF CLASSIC MODEL'
     &           //' / NETCDF3 FORMAT.')
         CASE(4)
            CALL logMessage(INFO,
     &       'UNIT 63 FORMAT WILL BE COMPACT ASCII.')
         CASE(5)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &         'UNIT 63 FORMAT WILL BE NETCDF CLASSIC MODEL'
     &         //' / NETCDF4 (HDF5) FORMAT.')
         CASE(7)
            useXDMF = .true.
            call logMessage(INFO,'UNIT 63 FORMAT WILL BE XDMF.')              
         CASE(6,8:)
            call allMessage(ERROR,"This NOUTGE value is invalid.")
            call ADCIRC_Terminate()
         CASE DEFAULT
            ! do nothing, the other cases handled below
      END SELECT

C...  IF GLOBAL ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTGE.NE.0) THEN

C...  COMPUTE NTCYSGE, NTCYFGE, WHICH = TOUTSGE AND TOUTFGE IN TIMESTEPS
#ifdef IBM
         NTCYSGE=INT((TOUTSGE-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &        KIND(0.0d0))
         NTCYFGE=INT((TOUTFGE-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &        KIND(0.0d0))
#else
         NTCYSGE=INT((TOUTSGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFGE=INT((TOUTFGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
         IF(NTCYFGE.GT.NT) NTCYFGE=NT

C...  CALCULATE NDSETSE = THE # OF DATA SETS TO BE SPOOLED TO UNIT 63

         IF(NSPOOLGE.EQ.0) NDSETSE=0
#ifdef IBM
         IF(NSPOOLGE.NE.0) NDSETSE=INT((NTCYFGE-NTCYSGE)/NSPOOLGE,
     &        KIND(0.0d0))
#else
         IF(NSPOOLGE.NE.0) NDSETSE=INT((NTCYFGE-NTCYSGE)/NSPOOLGE)
#endif

C...  WRITE NOUTGE,TOUTSGE,TOUTFGE,NTCYSGE,NTCYFGE,NSPOOLGE TO UNIT 16

         WRITE(16,3304) TOUTSGE,NTCYSGE,TOUTFGE,NTCYFGE,NSPOOLGE
 3304    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 63 EVERY ',
     &              'NSPOOLGE =',I8,' TIME STEPS')
      ENDIF
C...
C...  INPUT INFORMATION ABOUT GLOBAL VELOCITY DATA OUTPUT
C...

C...  READ IN NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : IF NOUTGV<>0, GLOBAL VEL.
C...  OUTPUT IS SPOOLED TO UNIT 64 EVERY NSPOOLGV TIME STEPS BETWEEN
C...  TIMES TOUTSGV AND TOUTFGV; IF ABS(NOUTGV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      WRITE(16,3351) NOUTGV
 3351 FORMAT(////,1X,'GLOBAL NODAL VELOCITY INFORMATION OUTPUT : ',
     &     //,5X,'NOUTGV = ',I2)

C...  CHECK INPUT PARAMETER NOUTGV
      SELECT CASE(ABS(NOUTGV))
         CASE(0)
            ! IF STATION OUTPUT WILL NOT BE GENERATED
            CALL logMessage(INFO,
     &         'NO GLOBAL VELOCITY OUTPUT WILL BE SPOOLED.')
         CASE(1)
            CALL logMessage(INFO,'UNIT 64 FORMAT WILL BE ASCII.')
         CASE(2)
            CALL logMessage(INFO,'UNIT 64 FORMAT WILL BE BINARY.')
         CASE(3)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &          'UNIT 64 FORMAT WILL BE NETCDF CLASSIC MODEL'
     &           //' / NETCDF3 FORMAT.')
         CASE(4)
            CALL logMessage(INFO,
     &       'UNIT 64 FORMAT WILL BE COMPACT ASCII.')
         CASE(5)
            useNetCDF = .true.
            useNetCDFOutput = .true.
            CALL logMessage(INFO,
     &         'UNIT 64 FORMAT WILL BE NETCDF CLASSIC MODEL'
     &         //' / NETCDF4 (HDF5) FORMAT.')
         CASE(7)
            useXDMF = .true.
            CALL logMessage(INFO,'UNIT 64 FORMAT WILL BE XDMF.')
         CASE(6,8:)
            call allMessage(ERROR,"This NOUTGV value is invalid.")
            call ADCIRC_Terminate()
         CASE DEFAULT
            ! do nothing, the other cases handled below
      END SELECT

C...  IF GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTGV.NE.0) THEN

C...  COMPUTE NTCYSGV, NTCYFGV, WHICH = TOUTSGV AND TOUTFGV IN TIMESTEPS
#ifdef IBM
          NTCYSGV=INT((TOUTSGV-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &        KIND(0.0d0))
          NTCYFGV=INT((TOUTFGV-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &         KIND(0.0d0))
#else
          NTCYSGV=INT((TOUTSGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
          NTCYFGV=INT((TOUTFGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
          IF(NTCYFGV.GT.NT) NTCYFGV=NT

C...  CALCULATE NDSETSV = THE # OF DATA SETS TO BE SPOOLED TO UNIT 64

          IF(NSPOOLGV.EQ.0) NDSETSV=0
#ifdef IBM
          IF(NSPOOLGV.NE.0) NDSETSV=INT((NTCYFGV-NTCYSGV)/NSPOOLGV,
     &         KIND(0.0d0))
#else
          IF(NSPOOLGV.NE.0) NDSETSV=INT((NTCYFGV-NTCYSGV)/NSPOOLGV)
#endif
C...  WRITE NOUTGV,TOUTSGV,TOUTFGV,NTCYSGV,NTCYFGV,NSPOOLGV TO UNIT 16

          WRITE(16,3354) TOUTSGV,NTCYSGV,TOUTFGV,NTCYFGV,NSPOOLGV
 3354     FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGV =',F8.3,
     &         ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &         I9,' TIME STEPS INTO THE SIMULATION',
     &         //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGV =',F8.3,
     &         ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &         I9,' TIME STEPS INTO THE SIMULATION',
     &         //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 64 EVERY ',
     &         'NSPOOLGV =',I8,' TIME STEPS')

       ENDIF



C...  COMPUTE PARAMETERS FOR TIME VARIABLE WEIR OUTPUT
      IF(USE_TVW.AND.NOUT_TVW.NE.0)THEN
        ALLOCATE(TVW(1:MNP))
        TVW(:) = 0D0
#ifdef IBM      
        NTCYS_TVW=INT((TOUTS_TVW-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &      KIND(0.0d0))
        NTCYF_TVW=INT((TOUTF_TVW-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &      KIND(0.0d0))   
#else
        NTCYS_TVW=INT((TOUTS_TVW-STATIM)*(86400.D0/DTDP) + 0.5d0)
        NTCYF_TVW=INT((TOUTF_TVW-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif        
        IF(NTCYF_TVW.GT.NT)NTCYF_TVW=NT
        IF(NSPOOL_TVW.EQ.0)THEN
            NDSETS_TVW = 0
        ELSE
#ifdef IBM       
             NDSETS_TVW=INT((NTCYF_TVW-NTCYS_TVW)/NSPOOL_TVW,KIND(0D0))
#else
             NDSETS_TVW=INT((NTCYF_TVW-NTCYS_TVW)/NSPOOL_TVW)
#endif             
        ENDIF
      ENDIF


             

C...
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION ABOUT GLOBAL
C...  CONCENTRATION DATA OUTPUT
C...
      NOUTGC=0
      IF(IM.EQ.10) THEN

C...  READ IN NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC : IF NOUTGC<>0, GLOBAL
C...  CONCENTRATION OUTPUT IS SPOOLED TO UNIT 73 EVERY NSPOOLGC TIME
C...  STEPS BETWEEN TIMES TOUTSGC AND TOUTFGC; IF ABS(NOUTGC)=2, OUTPUT
C...  WILL BE BINARY

         READ(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
         WRITE(16,3401) NOUTGC
 3401    FORMAT(////,1X,'GLOBAL NODAL CONCENTRATION INFORMATION OUTPUT:',
     &         //,5X,'NOUTGC = ',I2)

C...  CHECK INPUT PARAMETER NOUTGC
         SELECT CASE(ABS(NOUTGC))
            CASE(0)
               ! IF STATION OUTPUT WILL NOT BE GENERATED
               CALL logMessage(INFO,
     &         'NO GLOBAL CONCENTRATION OUTPUT WILL BE SPOOLED.')
            CASE(1)
               CALL logMessage(INFO,'UNIT 83 FORMAT WILL BE ASCII.')
            CASE(2)
               CALL logMessage(INFO,'UNIT 83 FORMAT WILL BE BINARY.')
            CASE(3:)
               call allMessage(ERROR,"This NOUTGV value is invalid.")
               call ADCIRC_Terminate()
            CASE DEFAULT
               ! do nothing, the other cases handled below
         END SELECT

C...  IF GLOBAL CONCENTRATION OUTPUT WILL BE GENERATED

         IF(NOUTGC.NE.0) THEN

C...  COMPUTE NTCYSGC, NTCYFGC, WHICH = TOUTSGC AND TOUTFGC IN TIMESTEPS
#ifdef IBM
           NTCYSGC=INT((TOUTSGC-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &          KIND(0.0d0))
           NTCYFGC=INT((TOUTFGC-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &          KIND(0.0d0))
#else
            NTCYSGC=INT((TOUTSGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFGC=INT((TOUTFGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
            IF(NTCYFGC.GT.NT) NTCYFGC=NT

C...  CALCULATE NDSETSC = THE # OF DATA SETS TO BE SPOOLED TO UNIT 73

            IF(NSPOOLGC.EQ.0) NDSETSC=0
#ifdef IBM
           IF(NSPOOLGC.NE.0) NDSETSC=INT((NTCYFGC-NTCYSGC)/NSPOOLGC,
     &          KIND(0.0d0))
#else
            IF(NSPOOLGC.NE.0) NDSETSC=INT((NTCYFGC-NTCYSGC)/NSPOOLGC)
#endif

C...  WRITE NOUTGC,TOUTSGC,TOUTFGC,NTCYSGC,NTCYFGC,NSPOOLGC TO UNIT 16

            WRITE(16,3404) TOUTSGC,NTCYSGC,TOUTFGC,NTCYFGC,NSPOOLGC
 3404       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 73 EVERY ',
     &                'NSPOOLGC =',I8,' TIME STEPS')
         ENDIF

      ENDIF

C...
C...  IF NWS<>0   INPUT INFORMATION ABOUT GLOBAL WIND DATA OUTPUT
C...
      IF(NWS.NE.0) THEN

C...  READ IN NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW : IF NOUTGW<>0, GLOBAL
C...  WIND OUTPUT IS SPOOLED TO UNIT 74 EVERY NSPOOLGW TIME STEPS
C...  BETWEEN TIMES TOUTSGW AND TOUTFGW; IF ABS(NOUTGW)=2, OUTPUT WILL
C...  BE BINARY

         READ(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
         WRITE(16,3451) NOUTGW
 3451    FORMAT(////,1X,'GLOBAL WIND STRESS INFORMATION OUTPUT : ',
     &           //,5X,'NOUTGW = ',I2)

C...  CHECK INPUT PARAMETER NOUTGW
         SELECT CASE(ABS(NOUTGW))
            CASE(0)
               CALL logMessage(INFO,
     &         'NO GLOBAL METEOROLOGICAL OUTPUT WILL BE SPOOLED.')
            CASE(1)
               CALL logMessage(INFO,
     &         'UNIT 73 AND 74 FORMATS WILL BE ASCII.')
            CASE(2)
               CALL logMessage(INFO,
     &          'UNIT 73 AND 74 FORMATS WILL BE BINARY.')
            CASE(3)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &          'UNIT 73 AND 74 FORMATS WILL BE NETCDF CLASSIC MODEL'
     &           //' / NETCDF3 FORMAT.')
            CASE(4)
               CALL logMessage(INFO,
     &          'UNIT 73 AND 74 FORMATS WILL BE COMPACT ASCII.')
            CASE(5)
               useNetCDF = .true.
               useNetCDFOutput = .true.
               CALL logMessage(INFO,
     &         'UNIT 73 AND 74 FORMATS WILL BE NETCDF CLASSIC MODEL'
     &         //' / NETCDF4 (HDF5) FORMAT.')
            CASE(7)
               CALL logMessage(INFO,
     &            'UNIT 73 AND 74 FORMATS WILL BE XDMF.')
            CASE(6,8:)
               call allMessage(ERROR,"This NOUTGW value is invalid.")
               call ADCIRC_Terminate()
            CASE DEFAULT
               ! do nothing, the other cases handled below
         END SELECT

C...  IF GLOBAL WIND STRESS OUTPUT WILL NOT BE GENERATED
         IF(NOUTGW.NE.0) THEN

C...  COMPUTE NTCYSGW, NTCYFGW, WHICH = TOUTSGW AND TOUTFGW IN TIMESTEPS
#ifdef IBM
            NTCYSGW=INT((TOUTSGW-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &          KIND(0.0d0))

            NTCYFGW=INT((TOUTFGW-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &          KIND(0.0d0))
#else
            NTCYSGW=INT((TOUTSGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFGW=INT((TOUTFGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
            IF(NTCYFGW.GT.NT) NTCYFGW=NT

C...  CALCULATE NDSETSW = THE # OF DATA SETS TO BE SPOOLED TO UNIT 74
            IF(NSPOOLGW.EQ.0) NDSETSW=0
#ifdef IBM
            IF(NSPOOLGW.NE.0) NDSETSW=INT((NTCYFGW-NTCYSGW)/NSPOOLGW,
     &          KIND(0.0d0))
#else
            IF(NSPOOLGW.NE.0) NDSETSW=INT((NTCYFGW-NTCYSGW)/NSPOOLGW)
#endif

C...  WRITE NOUTGW,TOUTSGW,TOUTFGW,NTCYSGW,NTCYFGW,NSPOOLGW TO UNIT 16

            WRITE(16,3454) TOUTSGW,NTCYSGW,TOUTFGW,NTCYFGW,NSPOOLGW
 3454       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGW =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGW =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 74 EVERY ',
     &                'NSPOOLGW =',I8,' TIME STEPS')

         ENDIF
      ENDIF

C...
C...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
C...
      READ(15,*) NFREQ
      WRITE(16,99392) NFREQ
99392 FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',
     &  //,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
      MNHARF = NFREQ

      IF (NFREQ.EQ.0) MNHARF = 1

C     allocate harmonic analysis arrays

      IF (NFREQ.GT.0) THEN
         CALL ALLOC_HA()
         CALL ALLOC_MAIN14()
      ENDIF

      IF(NFREQ.LT.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99391)
         WRITE(16,99391)
99391    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NFREQ (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL ADCIRC_Terminate()
C#ifdef CMPI
C         CALL ADCIRC_LOCALTERMINATE()
Cendif
         STOP
      ENDIF
      IF(NFREQ.GT.0) WRITE(16,2330)
 2330 FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.ARG(DEG)',
     +     1X,'CONSTITUENT',/)
      DO 1201 I=1,NFREQ
         READ(15,'(A10)') NAMEFR(I)
         READ(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331    FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
 1201 CONTINUE

C     read in interval information for harmonic analysis
C     compute thas and thaf in terms of the number of time steps

      READ(15,*) THAS,THAF,NHAINC,FMV
#ifdef IBM
      ITHAS=INT((THAS-STATIM)*(86400.D0/DTDP) + 0.5d0, KIND(0.0d0))
#else
      ITHAS=INT((THAS-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
      THAS=ITHAS*DTDP/86400.D0 + STATIM
#ifdef IBM
      ITHAF=INT((THAF-STATIM)*(86400.D0/DTDP) + 0.5d0, KIND(0.0d0))
#else
      ITHAF=INT((THAF-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
      THAF=ITHAF*DTDP/86400.D0 + STATIM
      ITMV = ITHAF - (ITHAF-ITHAS)*FMV
      IF(NFREQ.GT.0) THEN
         WRITE(16,34634) THAS,ITHAS,THAF,ITHAF,NHAINC
34634    FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER THAS =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'HARMONIC ANALYSIS WILL STOP AFTER THAF =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION'
     &        ,//,5X,'INFORMATION WILL BE ANALYZED EVERY ',
     &        'NHAINC =',I8,' TIME STEPS.')
         WRITE(16,34639) FMV*100.,ITMV
34639    FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE ',
     &        'FINAL ',F10.5,' %',/9X,'OF THE HARMONIC ANALYSIS ',
     &        'PERIOD OR AFTER ',I9,' TIME STEPS INTO THE ',
     &        'SIMULATION.',/9X,' RESULTS ARE WRITTEN TO UNIT 55.')

      ELSE
         WRITE(16,34645)
34645    FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
      ENDIF

      IF ((FMV.GT.0.).AND.(NFREQ.GT.0).AND.(C2DDI)) CHARMV = .TRUE.

C     read in and write out information on where harmonic analysis will
C     be done

      READ(15,*) NHASE,NHASV,NHAGE,NHAGV
      WRITE(scratchMessage,
     &   '("NHASE=",I1," NHASV=",I1," NHAGE=",I1," NHAGV=",I1,".")')
     &    NHASE, NHASV, NHAGE, NHAGV
      CALL logMessage(ECHO,scratchMessage)
      CALL checkHarmonicParameters()

C     compute flag telling whether any harmonic analysis will be done

      IHARIND=NFREQ*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C...
C...  Input information about hot start output
C...
C     jgf45.07 added undocumented option to STOP after writing hot start file.
C     This option will be used in testing ADCIRC's hot start capabilities.
      READ(15,*) NHSTAR,NHSINC

      CALL logMessage(INFO,'HOT START OUTPUT INFORMATION : ')
      WRITE(scratchMessage,'("NHSTAR=",I3," NHSINC=",I8,".")')
     &      NHSTAR, NHSINC
      CALL logMessage(ECHO,scratchMessage)
C
      SELECT CASE(NHSTAR)
      CASE(0)
         CALL logMessage(INFO,'HOT START OUTPUT WILL NOT BE GENERATED.')
      CASE(1,67,68)
         CALL logMessage(INFO,'HOT START OUTPUT WILL BE GENERATED'
     &        //' IN NON-PORTABLE BINARY FORMAT.')
      CASE(-1)  !tcm v51.26 added for time-stamped hotstart files
         CALL logMessage(INFO,'HOT START OUTPUT WILL BE GENERATED'
     &        //' IN NON-PORTABLE BINARY FORMAT IN TIME-STAMPED FILES.')
      CASE(3,367,368)
         useNetCDF = .true.
         CALL logMessage(INFO,'HOT START OUTPUT WILL BE GENERATED'
     &        //' IN PORTABLE NETCDF CLASSIC / NETCDF3 FORMAT.')
      CASE(5,567,568)
         useNetCDF = .true.
         CALL logMessage(INFO,'HOT START OUTPUT WILL BE GENERATED'
     &        //' IN PORTABLE NETCDF CLASSIC / NETCDF4 (HDF5) FORMAT.')
      CASE DEFAULT
         CALL allMessage(ERROR,"Input value of NHSTAR is invalid.")
         CALL ADCIRC_Terminate()
      END SELECT

      IF((NHSINC.EQ.0).AND.(NHSTAR.NE.0)) THEN
         CALL allMessage(ERROR,"Input value of NHSINC is 0"
     &     //" but the input value of NHSTAR is nonzero.")
         CALL allMessage(ERROR,"Please specify a time step increment"
     &     //" for writing hotstart files.")
         CALL ADCIRC_Terminate()
      ENDIF
C
      WRITE(16,34636) NHSINC
34636 FORMAT(/,5X,'HOT START OUTPUT WILL BE WRITTEN TO UNIT',
     &        ' 67 OR 68 EVERY ',I10,' TIME STEPS')
C
      IF((NHSTAR.EQ.67).OR.(NHSTAR.EQ.68).OR.
     &      (NHSTAR.EQ.367).OR.(NHSTAR.EQ.368).OR.
     &      (NHSTAR.EQ.567).OR.(NHSTAR.EQ.568) ) THEN
         WRITE(16,34626) NHSTAR
34626    FORMAT(/,5X,'ADCIRC will stop after writing to unit ',I3)
      ENDIF
C
      if (NHSINC <= 0) NHSINC = 1  ! rtm 46.xx NHSINC must have a
                                   ! reasonable value even when not
                                   ! generating hot start files.
C...
C...  Input information about GWCE solver
C...

C     read in and check matrix solver parameters

      READ(15,*) ITITER,ISLDIA,CONVCR,ITMAX

      WRITE(16,99656)
99656 FORMAT(//,1X,'SOLVER INFORMATION OUTPUT : ')
      IF((ISLDIA.LT.0).OR.(ISLDIA.GT.5)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9920)
         WRITE(16,9920)
 9920     FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',
     &                 '!!!!!!!!!',
     &           //,1X,'ISLDIA (A UNIT 15 INPUT PARAMETER) MUST BE 0-5',
     &            /,1X,'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9921)
            WRITE(16,9921)
 9921       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &                ' AND SET ISLDIA EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            ISLDIA=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            CALL ADCIRC_Terminate()
#endif
            STOP
         ENDIF
      ENDIF

C     jgf48.4619: Accommodate Seizos changes for explicit solve
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
C        allocate arrays needed by GWCE matrix and iterative solver
         call alloc_main11()
C        initialize parameter arrays needed by iterative solver
         CALL DFAULT(IPARM,RPARM)
         IPARM(1)=ITMAX
         IPARM(2)=ISLDIA
         !cms51.06: moved opening of fort.33 to openLogFile sub in global.F
         IPARM(4)=33
         RPARM(1)=CONVCR
         NW = 4*NP + 4*ITMAX
      ELSE ! lumped LHS
         call alloc_main11_lumped()
      ENDIF

C...
C...  Read input for 3D run
C...
      IF(C3D) THEN
        CALL READ_INPUT_3D(StaTim,NT)
c     ELSEIF(C3DDSS) THEN
c     CALL READ_INPUT_3DDSS(STATIM,NT)
      ENDIF
C     RJWW jgf46.00 Initialize nodal attributes, now that grid has been read
C     in from unit 14 file.
C     RJW IOOS move until after ZOB definition
      CALL InitNodalAttr(DP, NP, G, NScreen, ScreenUnit, MyProc, NAbOut)
CS

C...  INITIALIZE NIBNODECODE(I)
      DO I=1,NP
         NIBNODECODE(I)=0
      END DO

C - tcm v50.66.02 -- additions for time varying bathymetry
C - allocate global arrays for time varying bathymetry
      if (abs(nddt).ne.0) call ALLOC_MAIN13()

C - allocate arrays dealing with wind forcing
      call alloc_main12()

      if (mnproc == 1) then
         NP_G    = NP               !
         NE_G    = NE               !
         NSTAE_G = NSTAE
         NSTAV_G = NSTAV
         NSTAM_G = NSTAM
         NSTAC_G = NSTAC
         IF (C3D.eqv..true.) THEN
            NSTA3DD_G = NSTA3DD
            NSTA3DV_G = NSTA3DV
            NSTA3DT_G = NSTA3DT
         ENDIF
      endif

C     write table of ADCIRC parameter sizes
C tcm v51.09 added output for MNWPROH -- Number of Hot Start Writer Procs
      WRITE(16,4010)MNPROC,MNWPROC,MNWPROH,MNE,MNP,MNei,MNOPE,MNETA,
     & MNBOU,MNVEL,MNTIF,MNBFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,NWLAT,NWLON,
     & MNHARF,MNFFR
      IF(NWS.EQ.0) WRITE(16,4011)
      IF(NWS.EQ.1) WRITE(16,4012)
      IF(ABS(NWS).EQ.2) WRITE(16,4013)
      IF(NWS.EQ.3) WRITE(16,4014)
      IF(ABS(NWS).EQ.4) WRITE(16,4015)
      IF(ABS(NWS).EQ.5) WRITE(16,4115)
      IF(ABS(NWS).EQ.7) WRITE(16,4013) !46.00 Added NWS=7 (direct stress)
      IF(NWS.EQ.9) WRITE(16,4018)      !cf & cm added nws = 9 (jgf46.16 merged)
      IF(NWS.EQ.9) WRITE(16,4019)      !cf & cm added nws = 9 (jgf46.16 merged)
      IF(NWS.EQ.10) WRITE(16,4016)
      IF(NWS.EQ.11) WRITE(16,4017)
      IF(ABS(NWS).EQ.12) WRITE(16,4032) ! sb46.28sb01 Added NWS=12 (OWI format)
      IF(ABS(NWS).EQ.16) WRITE(16,4026) ! tcm v51.06.02  !gfdl met data
      IF(NWS.EQ.19) WRITE(16,4018)      !rjw added nws = 19
      IF(NWS.EQ.19) WRITE(16,4019)      !rjw added nws = 19
      IF(NWS.EQ.20) WRITE(16,4020)      !jie added nws = 20
      IF((NFREQ.EQ.0).OR.(FMV.EQ.0.)) WRITE(16,4021)
      IF((NFREQ.GE.1).AND.(FMV.NE.0.)) WRITE(16,4022)
      IF(ILUMP.EQ.0) WRITE(16,4031)
      IF(ILUMP.EQ.1) WRITE(16,4032)
      IF(IM.EQ.0) WRITE(16,4101)
      IF(IM.EQ.10) WRITE(16,4109)
      IF(IM.EQ.1) WRITE(16,4102)
      IF(IM.EQ.2) WRITE(16,4103)
      WRITE(16,4105)
      IF(USE_ELEV_STAT_FILE) WRITE(16,3180)  !tcm v51.20.05
      IF(USE_VEL_STAT_FILE) WRITE(16,3181)   !tcm v51.20.05 
      IF(USE_MET_STAT_FILE) WRITE(16,3182)   !tcm v51.20.05 
      IF(USE_CONC_STAT_FILE) WRITE(16,3183)  !tcm v51.20.05 
      WRITE(16,4108)
C
C tcm v51.09 added output for MNWPROH -- Number of Hot Start Writer Procs
 4010 FORMAT(' *****************************************************',/,
     &       ' *   Based on information extracted from the ADCIRC  *',/,
     &       ' *   UNIT 14 and 15 (grid and horiz run info) files  *',/,
     &       ' *   the following paramter values will be set:      *',/,
     &       ' *                                                   *',/,
     &       ' *       MNPROC = ',I5,1x,'     MWPROC = ',I5,7x,'   *',/,
     &       ' *       MWPROH = ',I5,4x,'                          *',/,
     &       ' *       MNE = ',I8,1X,'     MNP  = ',I8,1X,'        *',/,
     &       ' *       MNEI = ',I7,2X,'                            *',/,
     &       ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/,
     &       ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/,
     &       ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/,
     &       ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/,
     &       ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/,
     &       ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/,
     &       ' *       MNHARF = ',I5,4X,'  MNFFR = ',I6,3X,'       *',/,
     &       ' *                                                   *')
 4011 FORMAT(' *   Also, NO wind forcing will be used,             *')
 4012 FORMAT(' *   Also, Standard wind stress and pres will be used,*')
 4013 FORMAT(' *   Also, Semi-standard wind forcing will be used,  *')
 4014 FORMAT(' *   Also, Fleet numeric wind forcing will be used,  *')
 4015 FORMAT(' *   Also, PBL/JAG wind forcing will be used,        *')
 4026 FORMAT(' *   Also, GFDL Met Data wind and pres will be used, *')
 4115 FORMAT(' *   Also, Standard wind vel and pres will be used,  *')
 1236 FORMAT(' *   Also, surface stress forcing will be used,      *')
 4016 FORMAT(' *   Also, AVN wind & pressure forcing will be used, *')
 4017 FORMAT(' *   Also, ETA wind & pressure forcing will be used, *')
 4033 FORMAT(' *   Also, OWI format wind vel and pres will be used,*')
 4018 FORMAT(' *   Asymmetric hurricane wind and pressure forcing  *')
 4019 FORMAT(' *              will be used,                        *')
 4020 FORMAT(' *   Generalized Asymmetric Vortex Model forcing     *')
 4021 FORMAT(' *   means and variance calculation will NOT be made,*')
 4022 FORMAT(' *   means and variance calculation will be made,    *')
 4031 FORMAT(' *   the GWCE matrix will be left in consistent form *')
 4032 FORMAT(' *   the GWCE matrix will be LUMPED                  *')
 4101 FORMAT(' *   the model will be set up for a 2DDI run,        *')
 4109 FORMAT(' *   the model will be set up for a 2DDI run + transp*')
 4102 FORMAT(' *   the model will be set up for a 3D-VS run,       *')
 4103 FORMAT(' *   the model will be set up for a 3D-DSS run,      *')
 4105 FORMAT(' *   and an iterative solver will be used            *')
 3180 FORMAT(' *   An external elevation station file is used      *')  !tcm v51.20.05
 3181 FORMAT(' *   An external velocity station file is used       *')  !tcm v51.20.05
 3182 FORMAT(' *   An external met. station file is used           *')  !tcm v51.20.05
 3183 FORMAT(' *   An external concentration station file is used  *')  !tcm v51.20.05
 4108 FORMAT(' *****************************************************',/)
C

      IF ((useNetCDF.eqv..true.).and.(NETCDF_AVAIL.eqv..false.)) THEN
         CALL allMessage(ERROR,"NetCDF input and/or output was"
     &     //" indicated in the control parameters of the fort.15 file"
     &     //" but it is not supported by this ADCIRC executable"
     &     //" program. This program must be recompiled with NetCDF"
     &     //" libraries in order to enable NetCDF input or output.")
         CALL ADCIRC_Terminate()
      ENDIF
      
      if ((useNetCDFOutput.eqv..true.).and.
     &    (WRITE_LOCAL_FILES.eqv..true.)) then
         call allMessage(ERROR,'Some of the output file format '
     &      //'specifications in the fort.15 were for NetCDF format. '
     &      //'However, the command line option -L was also used to '
     &      //'specify local (subdomain) output files. The problem is '
     &      //'that ADCIRC cannot produce subdomain output files in '
     &      //'NetCDF format. '
     &      //'Please change output file formats to ASCII in the ' 
     &      //'fort.15 so the output files can be written '
     &      //'as subdomain output files. ' 
     &      //'Alternatively, remove the -L command line '
     &      //'option so that fulldomain output files are produced '
     &      //'in NetCDF format.')
         call adcirc_terminate()
      endif
         

      IF (useNetCDF.eqv..true.) THEN
         CALL logMessage(INFO,"Now reading metadata from fort.15 file"
     &     //" for use in NetCDF files.")
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) title
         CALL logMessage(ECHO,"metadata, title: "//trim(title))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) institution
         CALL logMessage(ECHO,"metadata, institution: "
     &     //trim(institution))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) source
         CALL logMessage(ECHO,"metadata, source: "//trim(source))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) history
         CALL logMessage(ECHO,"metadata, history: "//trim(history))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) references
         CALL logMessage(ECHO,"metadata, references: "
     &     //trim(references))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) comments
         CALL logMessage(ECHO,"metadata, comments: "//trim(comments))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) host
         CALL logMessage(ECHO,"metadata, host: "//trim(host))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) convention
         CALL logMessage(ECHO,"metadata, convention: "//trim(convention))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) contact
         CALL logMessage(ECHO,"metadata, contact: "//trim(contact))
         READ(15,'(A80)',err=99999,end=99998,iostat=ios) base_date
         CALL logMessage(ECHO,"metadata, base_date: "//trim(base_date))
      ENDIF
C...
C...  CLOSE FORT.15
C...
      CLOSE(15)

C
C... v49.48.02 tcm -- Deallocating rmax, bcxy and kdresults
C...                 used for searching and finding points in elements
      IF(ALLOCATED(KDRESULTS)) DEALLOCATE(KDRESULTS)
      call freeMesh()

C... zc50.80 - All compute processors pass through here and check if anyone
C              encountered an error while reading the input files.
C#ifdef CMPI      
C      CALL ADCIRC_CHECKLOCALTERMINATE()
C#endif

C...
c.RJW merged 09/02/2008 Casey 071219: Added the following subroutine call to compute the RESELEM array.
C             The subroutine is located at the begining of the file 'massbal.F.'
C commented out until can fix for 3D only
C#ifdef CMPI
C      CALL COMPUTE_RESELEM
C#endif
c.
C
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      RETURN
C
C     jgf50.41: This section is where we jump if there was an error
C     reading a file.
99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
99999 call allMessage(ERROR,"I/O error during file access.") !  ERR jumps here
      if (ios.gt.0) then
         write(scratchMessage,'(A,I3,A)')
     &       'The value of the i/o error flag was ',ios,'.'
         call allMessage(ERROR,scratchMessage)
      endif
C#ifdef CMPI
C      CALL ADCIRC_LOCALTERMINATE()
C#else      
      CALL ADCIRC_Terminate()
C#endif      
C
c******************************************************************************
      END SUBROUTINE READ_INPUT
c******************************************************************************



c******************************************************************************
C   Subroutine to read in 3D portion of fort.15 file                          *
C                                                                             *
C  Note, initial conditions on density, temperature and/or salinity are read  *
C        in for a cold start in subroutine COLD_START_3D                      *
C                                                                             *
c******************************************************************************

      SUBROUTINE READ_INPUT_3D(StaTime,NT)
      USE SIZES
C    kmd48.33bc - added the variables from Global
      USE GLOBAL, ONLY: RES_BC_FLAG, BCSTATIM, RBCTIME1,
     &                  RBCTIME2, SBCTIME1, SBCTIME2, TBCTIME1,
     &                  TBCTIME2, RBCTIMEINC, SBCTIMEINC,
     &                  TBCTIMEINC, SBCSTATIM, TBCSTATIM,
     &                  BCFLAG_LNM, BCFLAG_TEMP, TTBCSTATIM,
     &                  TTBCTIMEINC, SPONGEDIST, Sponge,
     &                  scratchMessage, INFO, DEBUG, ERROR, allMessage,
     &                  setMessageSource, logMessage, unsetMessageSource,
     &                  OUTPUTSPONGE, screenMessage, ECHO, sec2day
      USE MESH, ONLY : NP, DP, X, Y, CPP
      USE BOUNDARIES, ONLY : NETA, NOPE, NBD
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE GLOBAL_3DVS
      USE NodalAttributes, ONLY : Z0b_var

#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE !jgf45.09 added
C...  Declaration and definition of local variables used in this subroutine

      REAL(8), intent(in)  :: StaTime          !Model start time
      INTEGER, intent(in)  :: NT   !Total number of time steps in model run
C
      REAL(SZ) :: HH1  !domain averaged depth used for some vertical FE grids
      INTEGER  :: N,K,NN,J         !loop counters
      INTEGER  :: NH               !horizontal node loop counter
      CHARACTER(len=80) :: CDUM80

C    kmd48.33bc - additional variables for sponge layers
      REAL(SZ) :: Xloc1, Xloc2, Yloc1, Yloc2, Xloc, Yloc
      REAL(SZ) :: distloc, distpoint, slope1, slope2, xpoint
      REAL(SZ) :: ypoint, xpart, ypart
      REAL(SZ),ALLOCATABLE :: compdist(:)
      INTEGER :: BCnode, counter
      INTEGER, ALLOCATABLE :: nodedist(:)
      INTEGER :: IDen3D ! second instance of IDEN in the fort.15 file
      INTEGER :: ios   ! i/o status of read operation
C
      call setMessageSource("read_input_3D")
#if defined(READ_INPUT_TRACE) || defined (ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      ios = 0
C
C...  General invalid entry format statement used throughout
C
 350  FORMAT(//,2X,'***** INVALID ENTRY IN THE 3D INPUT SECTION OF ',
     &     ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

C...
C...  BEGIN READING VERTICAL PARAMETER INFORMATION
C...
      WRITE(16,300)
 300  FORMAT(//,1X,'3D SOLUTION INFORMATION',/)

C...  SPECIFY WHETHER A BAROTROPIC OR BAROCLINIC RUN
      !
      ! jgf51.52.35: Made this consistent with the previous read of IDEN.
      ! Not sure why we read IDEN twice, but need to keep this for backward
      ! compatibility. 
      READ(15,*) IDen3D
      if (IDen3D.ne.IDEN) then
         call allMessage(ERROR,'Both IDEN values must be the same.')
         call adcirc_terminate()
      endif
      if (IDen.GT.0) then
         C3D_BTrans  = .TRUE.
      endif

C...  READ IN THE TYPE OF BOTTOM BOUNDARY CONDITION AND THE SLIP COEFFICIENTS

      READ(15,*) ISlip,KP
      WRITE(16,355) ISlip,KP
 355  FORMAT(/,5X,'ISlip = ',I3,' KP = ',E12.5)
      IF((ISlip.LT.0).OR.(ISlip.GT.3)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,360)
        ENDIF
        WRITE(16,350)
        WRITE(16,360)
 360    FORMAT(/,2X,'    The Bottom Slip Code Must = 0,1,2,OR 3')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

C...  READ IN THE SURFACE AND BOTTOM ROUGHNESSES
C
C     Made Z0B bottom roughness a nodal attribute read in from fort.13 file if desired
C
      READ(15,*) Z0S, Z0B
      WRITE(16,380) Z0S,Z0B
 380  FORMAT(/,5X,'Z0S = ',E12.5,' Z0B = ',E12.5)


C...  READ IN THE TIME STEPPING COEFFICIENTS AND COMPUTE ASSOCIATED VARIABLES

      READ(15,*) Alp1,Alp2,Alp3
      WRITE(16,390) Alp1,Alp2,Alp3
 390  FORMAT(/,5X,'3D TIME STEPPING COEFFS Alp1 = ',E9.2,' Alp2 = ',
     &     E9.2,' Alp3 = ',E9.2)

      IDTAlp1 = iy*DelT*Alp1
      IDT1MAlp1 = iy*DelT*(1.-Alp1)
      DTAlp3 = DelT*Alp3
      DT1MAlp3 = DelT*(1-Alp3)
      DTAlp2 = DelT*Alp2
      DT1MAlp2 = DelT*(1.-Alp2)

C...  READ IN IGC & NFEN: F.E. GRID CODE & # NODES IN F.E. GRID

      READ(15,*,err=99999,end=99998,iostat=ios) IGC,NFEN
      WRITE(16,400) IGC,NFEN
 400  FORMAT(/,5X,'Vertical grid code IGC = ',I3,
     &            '  Number of vertical nodes (NFEN) = ',I5)
      IF((IGC.LT.0).OR.(IGC.GT.6)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,408)
        ENDIF
        WRITE(16,350)
        WRITE(16,408)
 408    FORMAT(/,2X,'    IGC MUST BE 0, 1, 2, 3, 4, 5, 6')
        CALL ADCIRC_Terminate()
      ENDIF
      IF(NFEN.LT.0) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,409)
        ENDIF
        WRITE(16,350)
        WRITE(16,409)
 409    FORMAT(/,2X,'    NFEN MUST BE > 0')
        CALL ADCIRC_Terminate()
      ENDIF

C...  SET MNFEN = NFEN

      MNFEN = NFEN

C...  ALLOCATE GENERAL 3D ARRAYS

      CALL ALLOC_3DVS()

C...  READ IN OR SET UP Vertical F.E. GRID

c     IGC = 0 - Read in grid from UNIT 15

      IF(IGC.EQ.0) then
        DO N=1,NFEN
          READ(15,*) Sigma(N)
        ENDDO
        IF(Sigma(1).NE.B) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(screenunit,350)
            WRITE(screenunit,1011)
            WRITE(screenunit,1012)
          ENDIF
          WRITE(16,350)
          WRITE(16,1011)
          WRITE(16,1012)
 1011     FORMAT(' Error reading in the vertical finite element grid')
 1012     FORMAT(' The first point in the finite element grid ',
     +           'must = b (-1) : run terminated'/)
          CALL ADCIRC_Terminate()
        ENDIF
        IF(Sigma(NFEN).NE.A) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(screenunit,350)
            WRITE(screenunit,1011)
            WRITE(screenunit,1013)
          ENDIF
          WRITE(16,350)
          WRITE(16,1011)
          WRITE(16,1013)
 1013     FORMAT(' The last point in the finite element grid ',
     +           'must = a (1) : run terminated'/)
          CALL ADCIRC_Terminate()
        ENDIF
      ENDIF

c     IGC <> 0 - Set up grid in subroutine FEGRIDS

      IF(IGC.NE.0) THEN
         HH1=0.d0
         DO NH=1,NP
            HH1=HH1+DP(NH)
         ENDDO
         HH1=HH1/NP                        !domain averaged depth
! tcm v49.74 -- removed NH from call to FEGRIDS which does not need it
         !CALL FEGRIDS(IGC,HH1,NH)
         CALL FEGRIDS(IGC,HH1)
      ENDIF

c     write out the vertical grid in fort.16 file

      WRITE(16,1000)
 1000 FORMAT(//,5X,'Vertical Grid Information')
      WRITE(16,1001)
 1001 FORMAT(/,5X,'V. Node #',5X,'V. Position',/)
      DO N = 1,NFEN
         WRITE(16,*) N,Sigma(N)
      ENDDO


C...  SPECIFY TYPE OF EDDY VISCOSITY PROFILE

      READ(15,*) IEVC,EVMin,EVCon
      WRITE(16,410) IEVC,EVMin,EVCon
 410  FORMAT(/,5X,'IEVC = ',I3,2X,'EVMin = ',E15.8,2X,'EVCon = ',E15.8)
      IF((IEVC.NE.0 ).and.(IEVC.NE.1 ).and.
     &   (IEVC.NE.10).and.(IEVC.NE.11).and.
     &   (IEVC.NE.20).and.(IEVC.NE.21).and.
     &   (IEVC.NE.22).and.(IEVC.NE.23).and.
     &   (IEVC.NE.30).and.(IEVC.NE.31).and.(IEVC.NE.32).and.
     &   (IEVC.NE.33).and.
     &   (IEVC.NE.40).and.(IEVC.NE.41).and.(IEVC.NE.42).and.
     &   (IEVC.NE.43).and.
     &   (IEVC.NE.50).and.
     &   (IEVC.NE.51)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,411)
        ENDIF
        WRITE(16,350)
        WRITE(16,411)
 411    FORMAT(/,2X,'    IEVC MUST BE 0,1,10,11,20,21,22,23,',
     &        '30,31,32,33,40,41,42,43,50,51')
        CALL ADCIRC_Terminate()
      ENDIF
      IF((IEVC.EQ.50).or.(IEVC.EQ.51)) THEN
         READ(15,*,err=99999,end=99998,iostat=ios) Theta1,Theta2
         WRITE(scratchMessage,'("theta1=",E15.8," theta2=",E15.8,".")')
     &      theta1, theta2   !tcm v50.85 20120829 changed E15.9 to E15.8 
      CALL logMessage(ECHO,scratchMessage)
      ENDIF

C...  FOR IEVC=0, CONSTANT EDDY VISCOSITY, READ IN PROFILE

      IF(IEVC.EQ.0) THEN
        DO N=1,NFEN
          READ(15,*) EVTot(N)
        ENDDO
        WRITE(16,*) ' Vertical E.V. read in from UNIT 15'
      ENDIF

C
C...  READ IN 3D OUTPUT CONTROLS, COMPUTE NEEDED ANCILLARY PARAMETERS
C

C...  Format statements used for 3D Station output diagnostic information

 3108 FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
 3109 FORMAT(/,5X,'STATION #   ELEMENT',3X,'LAMBDA(DEG)',
     &             4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
 1880 FORMAT(8X,I6,5X,I9,2(2X,F14.2))
 1883 FORMAT(6X,I6,5X,I9,2(2X,F13.8),2X,2(1X,F13.2))
 9790 FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &       /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &      //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)

C.... General variable definitions

C     NE       - total number of elements in grid
C     Areas(K) - 2x area of element K

C...3D Station Density Output (fort.41)
C
C     TO3DSDS  - starting time in days for 3D station density output
C     TO3DSDF  - ending time in days for 3D station density output
C
      READ(15,*) I3DSD,TO3DSDS,TO3DSDF,NSpo3DSD

C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D density station"
      CALL checkAndFillIOParameters(I3DSD, CDUM80, staTime,
     &   NSpo3DSD, TO3DSDS, TO3DSDF, NTO3DSDS, NTO3DSDF, NDSET3DSD,
     &   N3DSD, I3DSDRec)
C
C.....WRITE Density station output information to UNIT 16
      IF (I3DSD.ne.0) THEN
         WRITE(scratchMessage,571) trim(CDUM80), NSpo3DSD,TO3DSDS,
     &                (TO3DSDS+NSpo3DSD*DTDP*Sec2Day),
     &                (NTO3DSDS+NSpo3DSD),TO3DSDF,NTO3DSDF,NSpo3DSD,
     &                NDSet3DSD
         CALL logMessage(INFO,scratchMessage)
 571     FORMAT(5X,A,' data output will start',I9,
     &              ' timesteps after day = ',F9.3,
     &         9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &             ' relative to the simulation Start Time.',
     &         9X,'Data will be written every ',I7,
     &              ' time steps at total of ',I7,' times.')
      ENDIF

C.... REGARDLESS OF WHETHER I3DSD=0, READ IN THE NUMBER OF 3D DENSITY
C.... RECORDING STATIONS

      READ(15,*) NSta3DD
      IF(I3DSD.NE.0) WRITE(16,581) NSta3DD
 581  FORMAT(/,5X,' Output will be written at ',I7,' stations.')

      IF (NSta3DD.EQ.0) THEN
         MNSta3DD=1
      ELSE
         MNSta3DD = NSta3DD
      ENDIF
C
C     Allocate arrays for station density output
      call alloc_3DSD()
C
C.... Input the coordinates of the density output stations
C.... and then compute the element # containing each station
      IF(NSta3DD.GT.0) THEN
         CALL readStations(STATNAMED, NSta3DD, NE3DD, XED, YED,
     &         SLED, SFED, StaI3DD1, StaI3DD2, StaI3DD3,15,
     &                      '3D DENSITY RECORDING STATION  ')
      ENDIF
C
C...3D Station Velocity Output (fort.42)
C
C     TO3DSVS  - starting time in days for 3D station velocity output
C     TO3DSVF  - ending time in days for 3D station velocity output
C
      READ(15,*) I3DSV,TO3DSVS,TO3DSVF,NSpo3DSV
C
C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D velocity station"
      CALL checkAndFillIOParameters(I3DSV, CDUM80, staTime,
     &   NSpo3DSV, TO3DSVS, TO3DSVF, NTO3DSVS, NTO3DSVF, NDSET3DSV,
     &   N3DSV, I3DSVRec)
C
C.....WRITE Velocity station output information to UNIT 16
      IF (I3DSV.ne.0) THEN
         WRITE(scratchMessage,571) CDUM80,NSpo3DSV,TO3DSVS,
     &                (TO3DSVS+NSpo3DSV*DTDP*Sec2Day),
     &                (NTO3DSVS+NSpo3DSV),TO3DSVF,NTO3DSVF,NSpo3DSV,
     &                NDSet3DSV
         CALL logMessage(INFO,scratchMessage)
      ENDIF
C
C.... REGARDLESS OF WHETHER I3DSV=0, READ IN THE NUMBER OF 3D VELOCITY
C.... RECORDING STATIONS
      READ(15,*) NSta3DV
      IF(I3DSV.NE.0) WRITE(16,582) NSta3DV
 582  FORMAT(/,5X,' Output will be written at ',I7,' stations.')
C
      IF (NSta3DV.EQ.0) THEN
        MNSta3DV=1
      ELSE
        MNSta3DV = NSta3DV
      ENDIF
C
C  Allocate arrays for station velocity output
      call alloc_3DSV()
C
C....Input the coordinates of the velocity output stations
C....and then compute the element # containing each station
      IF(NSta3DV.GT.0) THEN
         CALL readStations(STATNAMEV3D, NSta3DV, NE3DV, XE3DV, YE3DV,
     &         SLE3DV, SFE3DV, StaI3DV1, StaI3DV2, StaI3DV3,15,
     &                      '3D VELOCITY RECORDING STATION ')
      ENDIF
C
C...3D Station Turbulence Output  (fort.43)
C     TO3DSTS  - starting time in days for 3D station turbulence output
C     TO3DSTF  - ending time in days for 3D station turbulence output
      READ(15,*) I3DST,TO3DSTS,TO3DSTF,NSpo3DST
C
C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D turbulence station"
      CALL checkAndFillIOParameters(I3DST, CDUM80, staTime,
     &   NSpo3DST, TO3DSTS, TO3DSTF, NTO3DSTS, NTO3DSTF, NDSET3DST,
     &   N3DST, I3DSTRec)
C
C.....  Write turbulence station output information to UNIT 16
      IF (I3DST.ne.0) THEN
         WRITE(scratchMessage,571) CDUM80, NSpo3DST,TO3DSTS,
     &                (TO3DSTS+NSpo3DST*DTDP*Sec2Day),
     &                (NTO3DSTS+NSpo3DST),TO3DSTF,NTO3DSTF,NSpo3DST,
     &                NDSet3DST
         CALL logMessage(INFO,scratchMessage)
      ENDIF
C
C.... REGARDLESS OF WHETHER I3DST=0, READ IN THE NUMBER OF 3D Turbulence
C.... STATIONS
      READ(15,*) NSta3DT
      IF(I3DST.NE.0) WRITE(16,583) NSta3DT
 583  FORMAT(/,5X,'Output will be written at ',I7,' stations')
C
      IF (NSta3DT.EQ.0) THEN
         MNSta3DT=1
      ELSE
         MNSta3DT = NSta3DT
      ENDIF
C
C  Allocate arrays for station turbulence output
      call alloc_3DST()
C
C....Input the coordinates of the turbulence output stations
C....and then compute the element # containing each station
      IF(NSta3DT.GT.0) THEN
         CALL readStations(STATNAMET, NSta3DT, NE3DT, XET, YET,
     &         SLET, SFET, StaI3DT1, StaI3DT2, StaI3DT3,15,
     &                      '3D TURBULENCE REC. STATION    ')
      ENDIF
C
C...3D Global Density Output (fort.44)
C
C     TO3DGDS  - starting time in days for 3D global density output
C     TO3DGDF  - ending time in days for 3D global density output
C
      READ(15,*) I3DGD,TO3DGDS,TO3DGDF,NSpo3DGD
C
C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D fulldomain density"
      CALL checkAndFillIOParameters(I3DGD, CDUM80, staTime,
     &   NSpo3DGD, TO3DGDS, TO3DGDF, NTO3DGDS, NTO3DGDF, NDSET3DGD,
     &   N3DGD, I3DGDRec)
C
C.....  Write global 3D Density output information to UNIT 16
      IF (I3DGD.ne.0) THEN
         WRITE(scratchMessage,571) CDUM80,
     &                NSpo3DGD,TO3DGDS,(TO3DGDS+NSpo3DGD*DTDP*Sec2Day),
     &                (NTO3DGDS+NSpo3DGD),TO3DGDF,NTO3DGDF,NSpo3DGD,
     &                NDSet3DGD
         CALL logMessage(INFO,scratchMessage)
      ENDIF
C
C...3D Global Velocity Output  (fort.45)
C
C     TO3DGVS  - starting time in days for 3D global velocity output
C     TO3DGVF  - ending time in days for 3D global velocity output
C
      READ(15,*) I3DGV,TO3DGVS,TO3DGVF,NSpo3DGV
C
C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D fulldomain velocity"
      CALL checkAndFillIOParameters(I3DGV, CDUM80, staTime,
     &   NSpo3DGV, TO3DGVS, TO3DGVF, NTO3DGVS, NTO3DGVF, NDSET3DGV,
     &   N3DGV, I3DGVRec)
C
C.....  Write global velocity output information to UNIT 16
      IF (I3DGV.ne.0) THEN
         WRITE(scratchMessage,571) CDUM80, NSpo3DGV,TO3DGVS,
     &                (TO3DGVS+NSpo3DGV*DTDP*Sec2Day),
     &                (NTO3DGVS+NSpo3DGV),TO3DGVF,NTO3DGVF,NSpo3DGV,
     &                NDSet3DGV
         CALL logMessage(INFO,scratchMessage)
      ENDIF
C
C...3D Global Turbulence Output  (fort.46)
C
C     TO3DGTS  - starting time in days for 3D global turbulence output
C     TO3DGTF  - ending time in days for 3D global turbulence output
C
      READ(15,*) I3DGT,TO3DGTS,TO3DGTF,NSpo3DGT
C
C.... CHECK INPUT PARAMETERS
      CDUM80 = "3D fulldomain turbulence"
      CALL checkAndFillIOParameters(I3DGT, CDUM80, staTime,
     &   NSpo3DGT, TO3DGTS, TO3DGTF, NTO3DGTS, NTO3DGTF, NDSET3DGT,
     &   N3DGT, I3DGTRec)
C
C.....  Write global turbulence output information to UNIT 16
      IF (I3DGT.ne.0) THEN
         WRITE(scratchMessage,571) CDUM80, NSpo3DGT,TO3DGTS,
     &                (TO3DGTS+NSpo3DGT*DTDP*Sec2Day),
     &                (NTO3DGTS+NSpo3DGT),TO3DGTF,NTO3DGTF,NSpo3DGT,
     &                NDSet3DGT
         CALL logMessage(INFO,scratchMessage)
      ENDIF
C
C    kmd48.33bc - added in the information of the 3D boundary conditions
C                 these boundary conditions are the level of no motion,
C                 salinity and temperature forcings.
      IF (CBAROCLINIC) THEN
         READ(15,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP !bc flags for elevation and temperature
         WRITE(16,429) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP
 429     FORMAT(/,5x,'RES_BC_FLAG = ',I3,5X,'BCFLAG_LNM = ',I3,5X,
     &    'BCFLAG_TEMP = ',I3)
         !
         !jgf51.52.35: Make sure the RES_BC_FLAG value matches the IDEN value.
         if (RES_BC_FLAG.ne.IDEN) then
            call allMessage(ERROR,
     & 'The value of RES_BC_FLAG must be the same as the value of IDEN.')
            call adcirc_terminate()
         endif         
         ! 
         IF ((RES_BC_FLAG.LT.0)) THEN  ! Diagnostic
            IF (ABS(RES_BC_FLAG).GE.1) THEN ! only have one set of values to read
               IF (NOPE.GT.0) THEN
                  READ(15,*) RBCTIMEINC
                  READ(15,*) BCSTATIM
                  WRITE(16,430) RBCTIMEINC
               END IF
            END IF
         ELSE IF ((RES_BC_FLAG.GT.0)) THEN
            IF (ABS(RES_BC_FLAG).EQ.1) THEN ! only read in the elevation changes
               IF (NOPE.GT.0) THEN
                  READ(15,*) RBCTIMEINC
                  READ(15,*) BCSTATIM
                  WRITE(16,430) RBCTIMEINC
 430             FORMAT(/,5X,'Read in elevation boundary conditions
     &                every', E9.2 ,'seconds')
               END IF
            ELSE IF ((ABS(RES_BC_FLAG).EQ.2)) THEN
               IF (NOPE.GT.0) THEN
                  READ(15,*) RBCTIMEINC, SBCTIMEINC
                  READ(15,*) BCSTATIM, SBCSTATIM
                  WRITE(16,431) RBCTIMEINC, SBCTIMEINC
 431              FORMAT(/,5X,'Read in elevation boundary conditions
     &                every', E9.2 ,'seconds and salinity
     &                boundary conditions every', E9.2,
     &                'seconds')
               END IF
            ELSE IF ((ABS(RES_BC_FLAG).EQ.3)) THEN
               IF (NOPE.GT.0) THEN
                  READ(15,*) RBCTIMEINC, TBCTIMEINC
                  READ(15,*) BCSTATIM, TBCSTATIM
                  WRITE(16,432) RBCTIMEINC, TBCTIMEINC
 432              FORMAT(/,5X,'Read in elevation boundary conditions
     &                every', E9.2 ,'seconds and temperature
     &                boundary conditions every', E9.2,
     &                'seconds')
                  IF (BCFLAG_TEMP.NE.0) THEN
                     READ(15,*) TTBCTIMEINC, TTBCSTATIM
                     WRITE(16,434) TTBCTIMEINC
 434                 FORMAT(/,5X,'Read in the top temperature
     &                   boundary condition every', E9.2 ,'
     &                   seconds.')
                  END IF
               END IF
            ELSE IF ((ABS(RES_BC_FLAG).EQ.4)) THEN
               IF (NOPE.GT.0) THEN
                  READ(15,*) RBCTIMEINC, SBCTIMEINC, TBCTIMEINC
                  READ(15,*) BCSTATIM, SBCSTATIM, TBCSTATIM
                  WRITE(16,433) RBCTIMEINC, SBCTIMEINC, TBCTIMEINC
 433              FORMAT(/,5X,'Read in elevation boundary conditions
     &                every', E9.2 ,'seconds and salinity
     &                boundary conditions every', E9.2,
     &                'seconds and temperature boundary conditions
     &                every', E9.2, 'seconds')
                  IF (BCFLAG_TEMP.NE.0) THEN
                     READ(15,*) TTBCTIMEINC, TTBCSTATIM
                     WRITE(16,434) TTBCTIMEINC
                 END IF
               END IF
            ELSE
               WRITE(16,350)
               WRITE(16,*) 'RES_BC_FLAG = ',RES_BC_FLAG
               WRITE(16,9722)
 9722         FORMAT(/,1X,'Your selection of RES_BC_FLAG (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
               CALL ADCIRC_Terminate()
            END IF
         END IF
      END IF

C    kmd48.33bc - added distance information for the sponge layer.
C                 Note that this is only used in the wind and advective terms

      IF (CBAROCLINIC) THEN
        READ(15,*) SPONGEDIST
        WRITE(16,435) SPONGEDIST
 435    FORMAT(/,5x,'SPONGEDIST = ',E9.2)
        IF (SPONGEDIST.GE.0.d0) THEN
           OUTPUTSPONGE=.TRUE.
        END IF
      END IF

      IF ((CBAROCLINIC).AND.(IDEN.NE.0)) THEN
         ALLOCATE (compdist(2),nodedist(2))
         DO J=1,NP
            counter=0
            DO K=1,2
              compdist(K)=0.d0
              nodedist(K)=0
            END DO
            DO N=1,NETA
              Xloc1=X(J)
              Yloc1=Y(J)
              BCnode=NBD(N)
              Xloc2=X(BCnode)
              Yloc2=Y(BCnode)
              Xloc=(Xloc2-Xloc1)
              Yloc=(Yloc2-Yloc1)
              Distloc=SQRT((Xloc*Xloc)+(Yloc*Yloc))
              IF (Distloc.GT.(2.d0*spongedist)) THEN
                CYCLE  ! don't need to store any information
              ELSE IF (Distloc.LE.(2.d0*spongedist)) THEN
                counter=counter+1
                IF (counter.LE.2) THEN
                   compdist(counter)=distloc
                   nodedist(counter)=NBD(N)
                ELSE
                 IF (distloc.LT.compdist(1)) THEN
                  compdist(1)=Distloc ! replace with new distance
                  nodedist(1)=NBD(N)
                  CYCLE
                 END IF
                 IF (distloc.LT.compdist(2)) THEN
                  compdist(2)=Distloc ! replace with new distance
                  nodedist(2)=NBD(N)
                  CYCLE
                 END IF
                END IF
              END IF
            END DO
            IF ((nodedist(1).EQ.J).OR.(nodedist(2).EQ.J)) THEN
               sponge(J)=0.d0   ! sponge layer will start with these nodes
            ELSE IF (counter.LT.2) THEN
               sponge(J)=1.d0   !  no sponge layer needed on these nodes
            ELSE
               IF (X(nodedist(1)).EQ.X(nodedist(2))) THEN
                  distpoint=ABS(X(J)-X(nodedist(1)))
               ELSE IF (Y(nodedist(1)).EQ.Y(nodedist(2))) THEN
                  distpoint=ABS(Y(J)-Y(nodedist(1)))
               ELSE
                  slope1=((Y(nodedist(2))-Y(nodedist(1)))/
     &                   (X(nodedist(2))-X(nodedist(1))))
                  slope2=-1/slope1
                  xpoint=((slope1*X(nodedist(1)))-(slope2*X(J))+
     &                   Y(J)-Y(nodedist(1)))/((slope1-slope2))
                  ypoint=slope1*(X(J)-X(nodedist(1)))+Y(nodedist(1))
                  xpart=(X(J)-xpoint)
                  ypart=(Y(J)-ypoint)
                  distpoint=SQRT((xpart*xpart)+(ypart*ypart))
               END IF
               IF (distpoint.GT.spongedist) THEN
                  sponge(J)=1.d0
               ELSE IF (distpoint.LE.spongedist) THEN
                  sponge(J)=distpoint/spongedist
               END IF
            END IF
         END DO
         DEALLOCATE(nodedist)
         DEALLOCATE(compdist)
      ELSE
         DO J=1,NP
           sponge(J)=1.d0
         END DO
      END IF
!kmd - end of additions

!Kendra: Add in information for equation of state
      IF (CBAROCLINIC) THEN
         READ(15,*) Eqnstate
         IF ((Eqnstate.EQ.2).OR.(Eqnstate.EQ.3)) THEN
            IF (ABS(IDEN).NE.4) THEN
              WRITE(16,424)
 424          FORMAT(/,1X,
     &        'Your selection of Eqnstate is not allowed with
     &               your choice of IDEN')
            END IF
         END IF
         IF (Eqnstate.eq.1) THEN
           WRITE(16,418)
 418       FORMAT(/,5X,
     &        'Equation of state uses the simple equation from Mellor')
         ELSE IF (Eqnstate.eq.2) THEN
           WRITE(16,419)
 419       FORMAT(/,5X,
     &'Equation of state uses the equation from McDougall et al(2003)')
         ELSE IF (Eqnstate.eq.3) THEN
           WRITE(16,420)
 420       FORMAT(/,5X,
     &'Equation of state uses the equation from UNESCO(1980)')
         ELSE
           WRITE(16,422)
 422       FORMAT(/,1X,'Your selection of Eqnstate (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
           CALL ADCIRC_Terminate()
         END IF
      END IF
C
!     Kendra45.12: Read in the new input values for the transport equation
C     jgf45.12: Made READs conditional on value of C3D_BTrans.

      if (C3D_Btrans) then

!     Kendra45.12: Must read in new values for lateral and vertical
!     diffusion
C...  READ IN NLSD, NLTD, NVTD & NVSD: Lateral and vertical diffusion
!     coefficients.
         READ(15,*) NLSD, NVSD
         WRITE(16,416) NLSD, NVSD
 416     FORMAT(/,5X,'Salinity Lateral Diffusion Coefficient = ',E9.2,
     &        'Salinity Vertical Diffusion Coefficient = ',E9.2)

         READ(15,*) NLTD, NVTD
         WRITE(16,417) NLTD, NVTD
 417     FORMAT(/,5X,'Temperature Lateral Diffusion Coefficient = ',
     &        E9.2,'Temperature Vertical Diffusion Coefficient = ',E9.2)

!     Kendra45.12: Read in the time stepping coefficient associated with the
!     transport equation terms.
         READ(15,*) ALP4
         WRITE(16,445) ALP4
 445     FORMAT(/,5X,'3D TIME STEPPING COEFFS ALP4 = ',E9.2)

         DTAlp4 = DelT*Alp4
         DT1MAlp4 = DelT*(1-Alp4)

C   kmd48.33bc - remove this boundary condition information due to the
C                new information being used.
!     Kendra45.12: Read in the temperature boundary condition file type
!c     jgf45.12: Made READ conditional on dynamic temperature forcing.
!         if ( IDEN .eq. 3 .or. IDEN .eq. 4 ) then
!            READ(15,*) NTF
!            WRITE(16,444) NTF
! 444        FORMAT(/,5X,'Temperature flux conditions are ', I7)
!         endif

      endif
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C
C
C     jgf50.41: This section is where we jump if there was an error
C     reading a file.
99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
99999 call allMessage(ERROR,"I/O error during file access.") !  ERR jumps here
      call allMessage(ERROR,
     &   "Check the fort.16 file for more information."
     &  //" Also, reducing the value of NABOUT to 0"
     &  //" will maximize the information written to the fort.16 file,"
     &  //" which may aid in troubleshooting this issue.")
      if (ios.gt.0) then
         write(scratchMessage,'(A,I3,A)')
     &       'The value of the i/o error flag was ',ios,'.'
         call allMessage(ERROR,scratchMessage)
      endif
      CALL ADCIRC_Terminate()

C-----------------------------------------------------------------------
      END SUBROUTINE READ_INPUT_3D
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C          C H E C K   A N D   F I L L   I O   P A R A M E T E R S
C-----------------------------------------------------------------------
C     jgf49.48.01 Checks i/o parameters and fills in parameters that
C     must be calculated.
C-----------------------------------------------------------------------
      SUBROUTINE checkAndFillIOParameters(specifier, description,
     &   staTime, tsPeriod, startTime, endTime, startTS, endTS, nSets,
     &   tsCounter, recCounter)
      USE GLOBAL, ONLY : scratchMessage, ECHO, INFO, WARNING, ERROR,
     &   DTDP, Day2Sec, NT, setMessageSource, logMessage, allMessage,
     &   unsetMessageSource, DEBUG, screenMessage, useNetCDF
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: specifier   ! format, 1=ascii, 2=binary, etc
      CHARACTER(len=80), intent(in) :: description ! type of output data
      REAL(8), intent(in) :: staTime
      INTEGER, intent(inout) :: tsPeriod ! period of time steps between outputs
      REAL(8), intent(inout) :: startTime! time for output to start (days)
      REAL(8), intent(inout) :: endTime  ! time for output to end (days)
      INTEGER, intent(inout) :: startTS  ! time step for output to start
      INTEGER, intent(inout) :: endTS    ! time step for output to end
      INTEGER, intent(out) :: nSets      ! num data sets in output file
      INTEGER, intent(out) :: tsCounter  ! time step counter btw outputs
      INTEGER, intent(out) :: recCounter ! counts lines in the output file
C
      call setMessageSource("checkAndFillIOParameters")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      tsCounter = 0  ! init
      recCounter = 0 ! init
C
      WRITE(scratchMessage,501) trim(description), specifier
      CALL logMessage(ECHO,scratchMessage)
 501  FORMAT(A,' output specifier = ',I2)
C
C     jgf49.48.01: Check to make sure that we can write in the
C     specified format.
      SELECT CASE(ABS(specifier))
         CASE(0)
            CALL logMessage(INFO,
     &         'OUTPUT WILL NOT BE SPOOLED FOR THIS FILE TYPE.')
         CASE(1)
            CALL logMessage(INFO,'OUTPUT FORMAT WILL BE ASCII.')
         CASE(2)
            CALL logMessage(INFO,'OUTPUT FORMAT WILL BE BINARY.')
         CASE(3)
            useNetCDF = .true.
            CALL logMessage(INFO,
     &          'OUTPUT FORMAT WILL BE NETCDF CLASSIC MODEL'
     &           //' / NETCDF3 FORMAT.')
         CASE(5)
            useNetCDF = .true.
            CALL logMessage(INFO,
     &         'OUTPUT FORMAT WILL BE NETCDF CLASSIC MODEL'
     &         //' / NETCDF4 (HDF5) FORMAT.')
         CASE(4,6:)
            WRITE(scratchMessage,350)
            CALL allMessage(ERROR, scratchMessage)
            WRITE(scratchMessage,511) specifier, trim(description)
            CALL allMessage(ERROR, scratchMessage)
 511        FORMAT('YOUR SETTING OF ',I2,' FOR THE ',A,
     &   ' OUTPUT PARAMETER IS NOT VALID. ','CHECK YOUR INPUT!!')
            call ADCIRC_Terminate()
         CASE DEFAULT
            ! do nothing, the other cases handled below
      END SELECT
C
C     jgf49.48.01: Check to make sure we have a valid output period
C     (if output was requested).
      IF ((specifier.NE.0).AND.(tsPeriod.EQ.0)) THEN
         WRITE(scratchMessage,350)
         CALL allMessage(ERROR,scratchMessage)
         WRITE(scratchMessage,561) description
         CALL allMessage(ERROR,scratchMessage)
 350     FORMAT('***** INVALID ENTRY IN THE 3D INPUT SECTION OF',
     &     ' FILE (UNIT 15) ****')
 561     FORMAT(' Time step increment for ',A,
     &      ' output data was 0, but it must be greater than zero.')
         CALL ADCIRC_Terminate()
      ENDIF
C
C     jgf49.48.01: Calculate output parameters and check their values.
      IF (specifier.NE.0) THEN
C....    COMPUTE startTS, endTS, WHICH = startTime AND endTime IN TIME STEPS
#ifdef IBM
         startTS=INT((startTime-StaTime)*Day2Sec/DTDP+0.5d0,
     &       KIND(0.0d0))      !jgf45.11 was NINT
         endTS=INT((endTime-StaTime)*Day2Sec/DTDP+0.5d0,
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
         startTS=INT((startTime-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
         endTS=INT((endTime-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif
C
C        jgf49.48.01: Check to make sure the start time step for output
C        is later than the actual start of the simulation.
         IF (startTS.LT.0) THEN
            WRITE(scratchMessage,531) description, startTime
            CALL allMessage(WARNING,scratchMessage)
 531        FORMAT('Start time for output of ',A,' data = ',E14.6,
     &         ' which is before the start time of the simulation. ',
     &         'It has been reset to coincide with the start time.')
            startTime=StaTime
            startTS=0
         ENDIF
C
C        jgf49.48.01: Check to make sure the end time step for output
C        is later than the start time for output.
         IF(endTS.LT.startTS) THEN
            WRITE(scratchMessage,541) description, endTime
            CALL allMessage(WARNING,scratchMessage)
 541        FORMAT('End time for output of ',A,' data = ',E14.6,
     &       ' which is before the start time for output of these data.'
     &       ' It has been reset to coincide with the start time.')
            endTime=startTime
            endTS=startTS
         ENDIF
C
C        jgf49.48.01: Check to see if the end time step for output
C        is later than the end of the simulation.
         IF(endTS.GT.NT) THEN
            WRITE(scratchMessage,551) description, endTime
            CALL logMessage(INFO,scratchMessage)
 551        FORMAT('End time for output of ',A,' data = ',E14.6,
     &         ' is later than the end of the simulation (RNDAY). ',
     &  'It has been reset to coincide with the end of the simulation.')
            endTS=NT
         ENDIF
         nSets = (endTS-startTS)/tsPeriod
      ENDIF
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call screenMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE checkAndFillIOParameters
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   R E A D   S T A T I O N S
C-----------------------------------------------------------------------
C     jgf45.12 Subroutine to read in station names and coordinates,
C     determine the element that the station falls in, and compute
C     interpolating factors.
C     tcm v49.48.01 -- added the description variable for correct
C     referencing of the type of station.
C     tcm v51.20.04 -- added stat_lun to specify the unit number
C     to read the station information from.  The default stat_lun = 15
C     but for external files it is stat_lun = 151.
C-----------------------------------------------------------------------
      SUBROUTINE readStations(names, num_stations, nnv, xcoord, ycoord,
     &                        lat, lon, sta1, sta2, sta3, stat_lun,
     &                        Description)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : DEG2RAD, RAD2DEG, parse,
     &   a2f, DEBUG, allMessage, setMessageSource, unsetMessageSource
      USE MESH, ONLY : ICS, CPP, SLAM0, SFEA0
      IMPLICIT NONE
      CHARACTER(50) :: names(num_stations)
      INTEGER num_stations, stat_lun
      INTEGER, dimension(num_stations) :: nnv
      REAL(SZ), dimension(num_stations) :: xcoord
      REAL(SZ), dimension(num_stations) :: ycoord
      REAL(SZ), dimension(num_stations) :: lat
      REAL(SZ), dimension(num_stations) :: lon
      REAL(SZ), dimension(num_stations) :: sta1
      REAL(SZ), dimension(num_stations) :: sta2
      REAL(SZ), dimension(num_stations) :: sta3
      INTEGER I
      CHARACTER(132) STATLINE
      CHARACTER(50) LVAR(3)
      CHARACTER(30), INTENT(IN) :: Description
C
      call setMessageSource("readStations")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      DO I=1,num_stations
         NNV(I)=0
         READ(stat_lun,'(A132)') STATLINE
         call parse(STATLINE, LVAR)
         names(I)=LVAR(3)
         IF(ICS.EQ.1) THEN
            xcoord(I)=a2f(LVAR(1))
            ycoord(I)=a2f(LVAR(2))
         ELSE
            lat(I)=a2f(LVAR(1))*DEG2RAD
            lon(I)=a2f(LVAR(2))*DEG2RAD
            CALL CPP(xcoord(I),ycoord(I),lat(I),lon(I),SLAM0,SFEA0)
         ENDIF
C
C... v49.48.02 -- tcm replaced with call to kdtsearch
C         CALL CoordinateToElement(xcoord(I), ycoord(I),
C     &         NNV(I), I, Description)
         CALL KDTSEARCH(xcoord(I), ycoord(I),
     &         NNV(I), I, Description)

         IF(ICS.EQ.1) THEN
            WRITE(16,1880) I,NNV(I),xcoord(I),xcoord(I)
         ELSE
            WRITE(16,1883) I,NNV(I),lat(I)*RAD2DEG,lon(I)*RAD2DEG,
     &            xcoord(I),ycoord(I)
         ENDIF

C....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT VEL. RECORDING STATIONS
         CALL ComputeInterpolatingFactors(xcoord(I), ycoord(I), NNV(I),
     &         sta1(I), sta2(I), sta3(I))

      END DO
1880  FORMAT(8X,I3,6X,I7,2(2X,F14.2))
1883  FORMAT(6X,I3,4X,I7,2(2X,F13.8),2X,2(1X,F13.2))
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readStations
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   K D T S E A R C H
C-----------------------------------------------------------------------
C
C  Subroutine that uses the KDTREE2 algorithm for finding
C      which element a point lies in.
C
C  Written by:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39056
C  Added in v49.48.02
C
C-----------------------------------------------------------------------

      SUBROUTINE kdtsearch(InputXCoordinate, InputYCoordinate,
     &     OutputElement, StationNumber, Description)
      use sizes, only : sz, MyProc
      use global, only : NFOver, NScreen, ScreenUnit, srchdp, tree, 
     &   kdresults, DEBUG, screenMessage, allMessage, setMessageSource, 
     &   unsetMessageSource
      use mesh, only : ne, nm, x, y, areas, rmax, bcxy 
      use adcirc_mod, only : adcirc_terminate
      use kdtree2_module
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      implicit none
      REAL(sz), intent(in) :: InputXCoordinate                  ! cartesian
      REAL(sz), intent(in) :: InputYCoordinate                  ! cartesian
      INTEGER, intent(out) :: OutputElement
      INTEGER, intent(in) :: StationNumber                     ! for err. mesg.
      CHARACTER(len=30), intent(in) :: Description             ! for err. mesg.

      INTEGER Element         ! element loop counter
      INTEGER ClosestElement  ! element with closest match
      INTEGER :: ielm(3),itc,iek
      REAL(sz) X1, X2, X3, X4, Y1, Y2, Y3, Y4,Xsta,Ysta       ! geometry
      REAL(sz) A1, A2, A3, AE, AREASK, AA            ! area
      real(sz) :: elmmin(2),xelm(3),yelm(3),dist
      LOGICAL ElementFound  ! .true. when a corresponding element is found

      REAL(sz), PARAMETER :: Tolerance = 1.0d-5     ! area difference for match
C
      call setMessageSource("kdtsearch")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ElementFound = .false.

      Xsta = InputXCoordinate
      Ysta = InputYCoordinate

      call kdtree2_n_nearest(tp=tree,qv=(/Xsta,Ysta/),
     &                  nn=srchdp,results=KDRESULTS)
C...    Check to see if the point lies with rmax of any of these elements

      ITC = 1
      ClosestElement = KDRESULTS(itc)%idx

      elmmin = minval(sqrt(KDRESULTS(1:srchdp)%dis)
     &                              - rmax(KDRESULTS(1:srchdp)%idx) )

      if(elmmin(1).le.0.0D0) then  ! Point lies within search radius of an element
C...        loop through the elements in the search list
         do while ((ElementFound.eqv..false.).and.(itc.le.srchdp))
            iek = KDRESULTS(itc)%idx  !Current search element number
C...           Get the distance from this point to the barycenter of the
C...           current element
            dist = sqrt(KDRESULTS(itc)%dis)
C...           If the distance is less than or equal to rmax (rmax=1.5*element radius)
C...           Then the point is near the element and might be in it
C...           Proceed with the weights test
            if(dist-rmax(iek).le.0.0d0) then
               !get the shape function for this element
               ielm(:) = NM(iek,(/1,2,3/))  !element's node numbers
               xelm(:) = X(ielm(:))      !element's vertex x-values
               yelm(:) = Y(ielm(:))      !element's vertex y-values
               X1=xelm(1)
               X2=xelm(2)
               X3=xelm(3)
               Y1=yelm(1)
               Y2=yelm(2)
               Y3=yelm(3)
               A1=(Xsta-X3)*(Y2-Y3)+(X2-X3)*(Y3-Ysta)
               A2=(Xsta-X1)*(Y3-Y1)-(Ysta-Y1)*(X3-X1)
               A3=(Ysta-Y1)*(X2-X1)-(Xsta-X1)*(Y2-Y1)
               AA=ABS(A1)+ABS(A2)+ABS(A3)
               AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
               AE=ABS(AA-AREASK)/AREASK
               IF (AE.LT.Tolerance) THEN
                  ElementFound = .true.
                  ClosestElement = iek
                  OutputElement = ClosestElement
               else !not in this element keep looking
                  itc = itc + 1
               endif !End area ratio test
            else !
C...             point is too far away from the barycenter of the
C...             element to possibly be in the element, so move to
C...             the next element
              itc = itc + 1
            endif !end Radius test
         enddo !end the while loop
      endif
      IF (.not. ElementFound ) THEN
         IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(ScreenUnit,9892) Description, StationNumber
         ENDIF
         WRITE(16,9892) Description, StationNumber
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9890) sqrt(KDRESULTS(1)%dis)
            ENDIF
            WRITE(16,9890) sqrt(KDRESULTS(1)%dis)
            OutputElement = ClosestElement
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9891) sqrt(KDRESULTS(1)%dis)
            ENDIF
            WRITE(16,9891) sqrt(KDRESULTS(1)%dis)
            call ADCIRC_Terminate()
         ENDIF
      ENDIF

 9892 FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &     'INPUT ERROR  !!!!!!!!!',//
     &     ,1X,A30,1X,I6,' DOES ',
     &     'NOT LIE WITHIN ANY ELEMENT IN THE DEFINED',
     &     /,1X,'COMPUTATIONAL DOMAIN.   PLEASE CHECK THE ',
     &     'INPUT COORDINATES FOR THIS STATION')

 9890 FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &     /,1X,'DISTANCE TO NEAREST ELEMENT IS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)

 9891 FORMAT(/,1X,'PROGRAM WILL NOT CORRECT ERROR ',
     &     'SINCE NON-FATAL ERROR OVERIDE OPTION, NFOVER,',
     &     /,1X,'HAS BEEN SELECTED EQUAL TO 0',
     &     /,1X,'DISTANCE TO NEAREST ELEMENT IS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',
     &     //)
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      END SUBROUTINE KDTSEARCH


C-----------------------------------------------------------------------
C     S U B R O U T I N E   C O O R D I N A T E  T O  E L E M E N T
C-----------------------------------------------------------------------
C
C     jgf45.12 Subroutine to take an X and Y cartesian coordinate and
C     find the corresponding element.
C
C-----------------------------------------------------------------------
      SUBROUTINE CoordinateToElement(InputXCoordinate, InputYCoordinate,
     &     OutputElement, StationNumber, Description)
      USE SIZES, ONLY : SZ, MyProc
      USE GLOBAL, ONLY: NFOver, NScreen, screenMessage, allMessage,
     & DEBUG, ScreenUnit, setMessageSource, unsetMessageSource
      USE MESH, ONLY : NE, NM, X, Y, Areas
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE

      REAL(8), intent(in) :: InputXCoordinate                  ! cartesian
      REAL(8), intent(in) :: InputYCoordinate                  ! cartesian
      INTEGER, intent(out) :: OutputElement
      INTEGER, intent(in) :: StationNumber                     ! for err. mesg.
      CHARACTER(len=30), intent(in) :: Description             ! for err. mesg.

      INTEGER Element         ! element loop counter
      INTEGER ClosestElement  ! element with closest match
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4       ! geometry
      REAL(8) A1, A2, A3, AE, AEMIN, AA            ! area
      LOGICAL ElementFound  ! .true. when a corresponding element is found
      REAL(8), PARAMETER :: Tolerance = 1.0E-5     ! area difference for match
C
      call setMessageSource("CoordinateToElement")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ElementFound = .false.
      AEMIN=1.0E+25
      ClosestElement=0
      DO Element=1,NE
         X1=X(NM(Element,1))
         X2=X(NM(Element,2))
         X3=X(NM(Element,3))
         X4=InputXCoordinate
         Y1=Y(NM(Element,1))
         Y2=Y(NM(Element,2))
         Y3=Y(NM(Element,3))
         Y4=InputYCoordinate
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-Areas(Element))/Areas(Element)
         IF (AE.LT.AEMIN) THEN
            AEMIN=AE
            ClosestElement=Element
         ENDIF
         IF (AE.LT.Tolerance) THEN
            ElementFound = .true.
            OutputElement=Element
         ENDIF
      ENDDO

      IF (.not. ElementFound ) THEN
         IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(ScreenUnit,593) Description, StationNumber
         ENDIF
         WRITE(16,593) Description, StationNumber
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9790) AEMIN
            ENDIF
            WRITE(16,9790) AEMIN
            OutputElement = ClosestElement
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9791) AEMIN
            ENDIF
            WRITE(16,9791) AEMIN
            call ADCIRC_Terminate()
         ENDIF
      ENDIF

 593  FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &     'INPUT ERROR  !!!!!!!!!',//
     &     ,1X,A30,1X,I6,' does ',
     &     'not lie within any element in the defined',
     &     /,1X,'computational domain.   PLEASE CHECK THE ',
     &     'INPUT COORDINATES FOR THIS STATION')
 9790 FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &     /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
 9791 FORMAT(/,1X,'PROGRAM WILL NOT CORRECT ERROR ',
     &     'SINCE NON-FATAL ERROR OVERIDE OPTION, NFOVER,',
     &     /,1X,'HAS BEEN SELECTED EQUAL TO 0',
     &     /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',
     &     //)
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE CoordinateToElement
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C     C O M P U T E  I N T E R P O L A T I N G  F A C T O R S
C-----------------------------------------------------------------------
C
C     jgf45.12 Subroutine to pre-compute the interpolating factors for a
C     recording station.
C
C-----------------------------------------------------------------------
      SUBROUTINE ComputeInterpolatingFactors(InputXCoordinate,
     &     InputYCoordinate, InputElement, Factor1, Factor2, Factor3)

      USE GLOBAL, ONLY : DEBUG, screenMessage,
     &    setMessageSource, unsetMessageSource, allMessage
      USE MESH, ONLY : NM, X, Y, Areas
      IMPLICIT NONE
      REAL(8), intent(in) :: InputXCoordinate                  ! cartesian
      REAL(8), intent(in) :: InputYCoordinate                  ! cartesian
      INTEGER, intent(in) :: InputElement
      REAL(8), intent(out):: Factor1, Factor2, Factor3

      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4                   ! geometry
C
      call setMessageSource("ComputeInterpolatingFactors")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      X1=X(NM(InputElement,1))
      X2=X(NM(InputElement,2))
      X3=X(NM(InputElement,3))
      X4=InputXCoordinate
      Y1=Y(NM(InputElement,1))
      Y2=Y(NM(InputElement,2))
      Y3=Y(NM(InputElement,3))
      Y4=InputYCoordinate

      Factor1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/Areas(InputElement)
      Factor2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/Areas(InputElement)
      Factor3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/Areas(InputElement)
C
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeInterpolatingFactors
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C                    S U B R O U T I N E   
C          L O G   N A M E L I S T   R E A D   S T A T U S
C-----------------------------------------------------------------------
C     jgf52.08.02: Record the i/o status associated with reading the
C     namelist using the negative, zero, or positive value of ios. 
C-----------------------------------------------------------------------
      subroutine logNamelistReadStatus(nmlname, ios)
      use global, only : logMessage, scratchMessage, DEBUG, ECHO, INFO,
     &   WARNING, ERROR, setMessageSource, unsetMessageSource,
     &   allMessage
      implicit none
      character(len=1000), intent(in) :: nmlname
      integer, intent(in) :: ios

      call setMessageSource("logNameListReadStatus")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      
      select case(ios)
      ! negative value indicates that we reached the end of the file
      ! before reading the namelist (not an error, since namelists
      ! are generally used for optional input)
      case(:-1) 
         call logMessage(INFO,
     &      'End-of-file when searching for '//trim(nmlName)//'.')
      ! zero indicates success
      case(0)
         call logMessage(INFO,
     &      'The '//trim(nmlName)//' namelist was found.')
      ! positive values indicate some sort of i/o error, other than
      ! reaching the end-file before finding the namelist
      case(1:)
         write(scratchMessage,'(a,i0,a)')
     &      'Could not read '//trim(nmlName)//
     &     ' namelist. The Fortran i/o error code was ',ios,'.'
         call allMessage(ERROR,scratchMessage)
      end select    
      
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine logNamelistReadStatus
C-----------------------------------------------------------------------

