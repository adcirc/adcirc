C******************************************************************************
C (P)ADCIRC VERSION 45.12 03/17/2006                                          *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C                                                                             *
C******************************************************************************
C                                                                             *
C                           (P)ADCIRC                                         *
C                                                                             *
C    A (PARALLEL) ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL SEAS        *
C                         AND ESTUARIES                                       *
C                                                                             *
C                                                                             *
C                          DEVELOPED BY:                                      *
C                                                                             *
C                      DR. R.A. LUETTICH, JR                                  *
C                                                                             *
C             UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL                     *
C                   INSTITUTE OF MARINE SCIENCES                              *
C                       3431 ARENDELL ST.                                     *
C                   MOREHEAD CITY, NC, 28557                                  *
C                         252-726-6841 EXT. 137                               *
C                   EMAIL  RICK_LUETTICH@UNC.EDU                              *
C                                                                             *
C                        DR. J.J. WESTERINK                                   *
C                                                                             *
C          DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES            *
C                     UNIVERSITY OF NOTRE DAME                                *
C                       NOTRE DAME, IN 46556                                  *
C                         219-631-6475                                        *
C                    EMAIL JJW@PHOTIUS.CE.ND.EDU                              *
C                                                                             *
C                                                                             *
C        MAJOR FUNDING FOR THE DEVELOPMENT OF ADCIRC WAS PROVIDED BY          *
C                                                                             *
C                       DEPARTMENT OF THE ARMY                                *
C                    WATERWAYS EXPERIMENT STATION                             *
C                 COASTAL ENGINEERING RESEARCH CENTER                         *
C                        3909 HALLS FERRY RD                                  *
C                      VICKSBURG, MI 39180-6199                               *
C                                                                             *
C******************************************************************************
C                                                                             *
C          THE ADCIRC SOURCE CODE IS COPYRIGHTED, 1994-2006 BY:               *
C                                                                             *
C                 R.A. LUETTICH, JR AND J.J. WESTERINK                        *
C                                                                             *
C         NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED             *
C               WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS                 *
C                                                                             *
C******************************************************************************
C                                                                             *
C                  Parallelization of ADCIRC 2D was done by                   *
C                    the center for subsurface modeling                       *
C                         The University of Texas                             *
C                             Austin, TX 78712                                *
C                            03/06/98 - 5/19/99                               *
C                                                                             *
C                   Translation of 2D code to Fortran90                       *
C                             Victor J. Parr                                  *
C                              John B. Romo                                   *
C                               8/31/99                                       *
C                                                                             *
C               Parallelization, consolidation of 2D & 3D codes,              *
C             translation of 3D routines to Fortran 90 was done by            *
C                               Tim Campbell                                  *
C                   Naval Research Lab, Stennis Space Center                  *
C                             summer 2002                                     *
C                                                                             *
C******************************************************************************


C******************************************************************************
C******************************************************************************
      MODULE ADCIRC_Mod
C******************************************************************************
C******************************************************************************
C
      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM
      USE GLOBAL_3DVS, ONLY : GORho, GORhoOAMB, AMB
      USE VERSION
#ifdef CMPI
      USE MESSENGER
      USE WRITER, ONLY: WRITER_MAIN, WRITER_INIT, WRITER_PAUSE !tcm added pause
      USE HSWRITER, ONLY: HSWRITER_MAIN, HSWRITER_INIT,!st3 for hstart
     &                    HSWRITER_PAUSE !tcm v49.52.01 added pause
#endif
      IMPLICIT NONE
      PUBLIC

      INTEGER, SAVE :: ITIME_BGN, ITIME_END


      CONTAINS


C******************************************************************************
C*    ADCIRC initialize routine
C******************************************************************************
      SUBROUTINE ADCIRC_Init(COMM, ROOTD)
      !USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      USE SIZES, ONLY : rootDir, mnproc, mnwproc, mnallproc,
     &  read_local_hot_start_files, make_dirname, write_local_harm_files
      USE MESH, ONLY : Areas, TotalArea, totalAreaCalc
      USE HARM, ONLY : allocateFullDomainHAIOArrays,
     &                 initHarmonicParameters
      USE GLOBAL_IO, ONLY : allocateFullDomainIOArrays
      USE WRITE_OUTPUT, ONLY : initOutput2D
      IMPLICIT NONE
C
      INTEGER, OPTIONAL :: COMM
      CHARACTER(*), OPTIONAL :: ROOTD
      REAL(SZ) :: StaTimHS
      REAL(SZ) :: RefTimHS
      REAL(SZ) :: TimeLoc
      INTEGER :: iargc    ! function to return the number of command line arguments
      INTEGER :: argcount ! number of command line arguments
      CHARACTER(2048) :: cmdlinearg ! a single command line argument
      INTEGER :: i        ! command line argument counter
      !
      ! jgf51.52.07: Added command line option to just write out
      ! the ADCIRC version number and exit. This is to allow 
      ! external applications to query the version number so that
      ! they can behave in version-appropriate ways. 
      argcount  = iargc()  ! get number of command line arguments
      if ((argcount.gt.0).and.(myproc.eq.0)) then
         i = 0
         do while (i.lt.argcount)
            i = i + 1
            call getarg(i, cmdlinearg)
            select case(cmdlinearg(1:2))
            case("-V","-v") 
               ! write the version number string ... actually it is
               ! a series of two or three numbers with two digits 
               ! separated by dots, e.g., 51.02 or 51.52.07
               write(6,'(a)') trim(adc_version)                 
               stop
            case default
               ! do nothing, this is some other command line argument
            end select
         end do
      endif
      !
      call initLogging()
      nabout = 0 !jgf52.08.01: Initialize log level to ECHO.
#ifdef DEBUG
      nabout = -1 !jgf52.08.01: Initialize log level to DEBUG if so compiled
#endif
      call setMessageSource("ADCIRC_Init")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call screenMessage(DEBUG,"Enter.") ! log to screen; don't have log dirname
#endif
      IF (PRESENT(ROOTD)) THEN
        ROOTDIR = TRIM(ROOTD)
      ELSE
        ROOTDIR = '.'
      ENDIF
C
#ifdef CMPI
      CALL GET_NUMWRITERS()  ! Get number of writer processors from command line arguments
      IF (PRESENT(COMM)) THEN
        CALL MSG_INIT(COMM)  ! Init MPI and get MPI-rank of this cpu
      ELSE
        CALL MSG_INIT()      ! Init MPI and get MPI-rank of this cpu
      ENDIF
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
      call openLogFile()     ! jgf50.65: open subdomain log files (fort.16)
      CALL WRITER_INIT()     ! Initialize WRITER module
      CALL HSWRITER_INIT()   ! Initialize HSWRITER module  !st3 for hstart  100711
      if(((mnwproh > 0).or.(mnwproc > 0)) .and. MYPROC >= MNPROC) then
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()     
         return
      endif
      CALL READ_INPUT()      ! Establish sizes by reading fort.14 and fort.15
      CALL MSG_TABLE()       ! Read Message-Passing Tables
      CALL MSG_START()       ! Startup message passing
#else
      MNPROC = 1             ! Init number of compute procs for serial
      MNWPROC = 0            ! Init number of writer procs for serial
      MNALLPROC = 1          ! Init number of all procs for serial
      MYPROC = 0             ! Init task id for serial
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
      call openLogFile()     ! jgf50.65: open fort.16 file
      CALL READ_INPUT()      ! Establish sizes by reading fort.14 and fort.15
#endif

      ! File Fmt Version  must match Major and Minor numbers.
      FileFmtVersion =
     &   VERSION_NUMBER(FileFmtMajor, FileFmtMinor, FileFmtRev)
      IF ((NSCREEN.NE.0).AND.(MYPROC.EQ.0)) THEN
         WRITE(ScreenUnit,'(a)')"ADCIRC Version is "//ADC_VERSION
      ENDIF

      !jgf49.44: Allocate memory for reading and writing full domain arrays
      IF ( (MNPROC.gt.1).and.
     &        ( (myProc.eq.0)
     &          .or.(READ_LOCAL_HOT_START_FILES.eqv..false.) ) ) THEN
          CALL allocateFullDomainIOArrays()
      ENDIF
      IF ( (MNPROC.gt.1).and.(myProc.eq.0).and.
     &     (WRITE_LOCAL_HARM_FILES.eqv..false.) ) THEN
          CALL allocateFullDomainHAIOArrays()
      ENDIF

      !jgf49.44: Initialize harmonic analysis parameters.
      CALL initHarmonicParameters()

#ifdef CMACHSUN
      CALL ABRUPT_UNDERFLOW() ! Abrupt underflow for Sun computers
#endif

C...Compute the reciprocal of the number of nodes in the entire domain

#ifdef CMPI
      CALL ALLNODES(NP_GLOBAL)
      RNP_GLOBAL = 1.0D0/NP_GLOBAL
#endif

C     jgf45.08 put this here so it can be used in hot start
      IF(NOLIFA.EQ.0) THEN
         IFNLFA=0
      ELSE
         IFNLFA=1
      ENDIF
C     kmd 48.33bc - include a new initial condition file for baroclinic
C...
C...******************** PROGRAM SETUP SECTION ************************
C...

      SELECT CASE(IHOT)
         CASE(0)
            CALL COLDSTART()
            ITHS = 0  !TCM20091029 v49.02 Initialized to zero for coldstarts
         ! kmd : added in the initial condition read
         CASE(17,67,68,367,368,567,568) ! non-portable binary or netcdf hotstart file
            CALL HOTSTART()
         CASE DEFAULT
            ! the IHOT value was vetted already, so this code should
            ! not be reachable
            IF ((NSCREEN.ne.0).and.(myProc.eq.0)) THEN
               WRITE(screenUnit,*) "ERROR: IHOT=",IHOT," is not valid."
            ENDIF
      END SELECT

C...Default beginning and ending time steps

      ITIME_BGN = ITHS+1
      ITIME_END = NT  !NT is set in read_input.F

C...  Determine the number of active elements (MJU), the total number of
C...  elements (NODELE) and the total area of elements (TotalArea) 
C...  attached to each node
      call totalAreaCalc()
      !
      !
      ! jgf52.08.20: Compute master time which is referenced to the beginning 
      ! time of the model run. This is needed here only because we need to 
      ! determine whether to load min/max files in subroutine initOutput2D().
      ! If the time in the min/max file is later than the current time
      ! in ADCIRC, we don't load the min/max files. TODO: This is computed
      ! again in subroutine timestep().
      TimeLoc=ITIME_BGN*DTDP + StaTim*86400.D0
      IF (CHOTHS.eqv..true.) THEN
         StaTimHS=((ITIME_BGN-1)*DTDPHS)/86400.D0
         RefTimHS=((ITIME_BGN-1)*DTDP)/86400.D0
         TimeLoc=ITIME_BGN*DTDP + (StaTimHS - RefTimHS)*86400.D0
      END IF
C
C     jgf51.21.24: Initialize output.
      call initOutput2D(timeloc)
C...
C...************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING ***********
C...

C...Set flags for nonlinear terms

      IF(NOLICA.EQ.0) THEN
         IFNLCT=0
      ELSE
         IFNLCT=1
      ENDIF
      IF(NOLICAT.EQ.0) THEN
         IFNLCAT=0
      ELSE
         IFNLCAT=1
      ENDIF

Corbitt 1203022: Local Advection Flags
      IFNLCTE  = IFNLCT
      IFNLCATE = IFNLCAT

      IFWIND=1
      IF(IM.EQ.1) IFWIND=0

C...Compute several constant coefficients

      GA00=G*A00
      GC00=G*C00
Corbitt 120322: Handled in Timestep.F elementally
C      TADVODT=IFNLCAT/Dt
      GB00A00=G*(B00+A00)
      GFAO2=G*IFNLFA/2.D0
      GO3=G/3.D0
      GO2=G/2.D0
      DTO2=DT/2.D0
      DT2=DT*2.D0
      GDTO2=G*DT/2.D0
Corbitt 120322: Handled in Timestep.F elementally
C      SADVDTO3=IFNLCT*DT/3.D0
      GORho=G/RhoWat0
      GORhoOAMB=GORho/AMB

C...
C.....Compute necesary wet/dry coefficients
C...
      HABSMIN=0.8d0*H0
      HOFF=1.2d0*H0

      WRITE(16,1112)
      WRITE(16,17931)
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,17931)
 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)

      NCChange=1 !jgf45.06 set flag to set up GWCE LHS on 1st iteration
C
 1999 format(1X, a," ADCIRC Version: ",a)
C
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      
      RETURN
      END SUBROUTINE ADCIRC_Init


C******************************************************************************
C*    ADCIRC run routine
C******************************************************************************
      SUBROUTINE ADCIRC_Run(NTIME_STP)
C
#ifdef CSWAN
Casey 090302: We need the following information to couple to unstructured SWAN.
      USE Couple2Swan, ONLY: CouplingInterval,
     &                       PADCSWAN_RUN
#endif
      USE WRITE_OUTPUT,ONLY: writeOutput2D
      IMPLICIT NONE
C
      INTEGER, OPTIONAL :: NTIME_STP
C
      call setMessageSource("ADCIRC_Run")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!TCM v49.52.01 added this ifdef (brought in from adcirc_init()
#ifdef CMPI
      if (((mnwproh > 0).or.(mnwproc > 0)) .and. myproc >= mnproc) then
         if( myproc < mnproc+mnwproc ) then   !st3
            call writer_main()
            call unsetMessageSource()  !tcm 20150924 -- added
            return
         else
            call hswriter_main()
            call unsetMessageSource()  !tcm 20150924 -- added
            return
         endif
      endif
#endif

      IF (PRESENT(NTIME_STP)) THEN
        ITIME_END = ITIME_BGN+NTIME_STP-1
      ENDIF
C
      call allMessage(INFO,"Begin timestepping.")
      ! print *, ' In adcirc.F  ...  line 383', 'ITIME_BGN,ITIME_END > ',ITIME_BGN,ITIME_END
      DO 100 ITIME=ITIME_BGN,ITIME_END
      ! print *, ' In adcirc.F  ...  line 384', 'ITIME_BGN,ITIME_END > ',ITIME_BGN,ITIME_END
        CALL TIMESTEP(ITIME)
C
#ifdef CSWAN
Casey 090302: If it is time, then call the following subroutine
C             to then call the SWAN time-stepping subroutine.
        IF(MOD(ITIME,CouplingInterval).EQ.0)THEN
          CALL PADCSWAN_RUN(ITIME)
        ENDIF
#endif
        !...WRITE OUTPUT
        !   zc - moved this here so that output occurs after
        !        the SWAN time step has been completed.
        CALL writeOutput2D(ITIME,ITIME*DTDP+StaTim*86400D0) ! =>see write_output.F
100   CONTINUE
C
      IF (PRESENT(NTIME_STP)) THEN
        ITIME_BGN = ITIME_END+1
      ENDIF
C

!TCM v49.52.01 added the following ifdef to pause writer processors
! before exiting adcirc_run phase.  This is to allow them to be
! used again if adcirc_run is called more than once.
#ifdef CMPI
      IF (MYPROC == 0 .and. MNWPROC > 0) CALL WRITER_PAUSE()
      IF (MYPROC == 0 .and. MNWPROH > 0) CALL HSWRITER_PAUSE()
#endif

#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      
      RETURN
      END SUBROUTINE ADCIRC_Run


C******************************************************************************
C*    ADCIRC finalize routine
C******************************************************************************
      SUBROUTINE ADCIRC_Final(NO_MPI_FINALIZE)
      USE HARM, ONLY : solveHarmonicAnalysis
      USE WRITE_OUTPUT, ONLY : writeHarmonicAnalysisOutput
      IMPLICIT NONE
C
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      INTEGER :: I
C
      call setMessageSource("ADCIRC_Final")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C...
C...*************** SOLVE THE HARMONIC ANALYSIS PROBLEM ****************
C...
      CALL solveHarmonicAnalysis(ITIME)
      CALL writeHarmonicAnalysisOutput(ITIME)
C
#ifdef CMPI
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
C
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END SUBROUTINE ADCIRC_Final


C******************************************************************************
C*    ADCIRC terminate routine
C*    jgf50.44: Mostly just needed so that we can clean up mpi when
C*    bombing out.
C******************************************************************************
      SUBROUTINE ADCIRC_Terminate(NO_MPI_FINALIZE)
      IMPLICIT NONE
C
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      call setMessageSource("ADCIRC_Terminate")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      STOP 
C
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      END SUBROUTINE ADCIRC_Terminate

C*****************************************************************************
C  ADCIRC Local Terminate Routine and corresponding collection routine
C  zc50.81 In the case where there is an error on a single processor, allow
C          everyone to check if we need to bomb out. For example, mesh errors
C          might only occur on a single processor and hang the code. Care
C          must be taken not to inadvertantly hang the code by allowing for
C          cases where not all compute processors pass through the same call
C          to this routine.
C*****************************************************************************
#ifdef CMPI
      SUBROUTINE ADCIRC_LocalTerminate()
        USE GLOBAL,ONLY:TERMINATE_LOCALPROC
        IMPLICIT NONE
        !...Set the local terminator and then go wait for the rest of
        !   the compute processors to check in
        TERMINATE_LOCALPROC = 1
        CALL ADCIRC_CheckLocalTerminate()
      END SUBROUTINE

      SUBROUTINE ADCIRC_CheckLocalTerminate()
        USE MESSENGER
        USE GLOBAL,ONLY:TERMINATE_LOCALPROC
#ifdef HAVE_MPI_MOD
        USE MPI
#endif        
        IMPLICIT NONE
        INTEGER :: TERM_LOCAL(2)
        INTEGER :: TERM_GLOBAL(2)
#ifndef HAVE_MPI_MOD
        INCLUDE 'mpif.h'
#endif
        TERM_LOCAL(1) = TERMINATE_LOCALPROC
        TERM_LOCAL(2) = MYPROC
        CALL MPI_ALLREDUCE(TERM_LOCAL,TERM_GLOBAL,1,
     &      MPI_2INTEGER,MPI_MAXLOC,COMM,IERR)
        IF(TERM_GLOBAL(2).NE.0)THEN
            IF(MYPROC.EQ.0)THEN
                WRITE(ScreenUnit,'(A,I5)') "ERROR: ADCIRC Terminating "//
     &          "from error originating from MYPROC = ",TERM_GLOBAL(2)
                WRITE(ScreenUnit,'(A)')    "       Check the local PE "//
     &          "fort.16 output file for details."
            ENDIF
            CALL ADCIRC_TERMINATE()
        ENDIF
      END SUBROUTINE   
#endif

C******************************************************************************
C******************************************************************************
      END MODULE ADCIRC_Mod
C******************************************************************************
C******************************************************************************

