!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
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
#ifdef CMPI
      USE mpi
#endif

      IMPLICIT NONE
      PUBLIC

      INTEGER, SAVE :: ITIME_BGN, ITIME_END

      REAL (8), private:: init_wtime = 0.D0, run_wtime = 0.D0, 
     &   final_wtime = 0.D0, total_wtime = 0.D0
      REAL (8), private:: start_total_wtime = 0.D0, 
     &   end_total_wtime = 0.D0

      PRIVATE::  collect_run_wall_time

      CONTAINS


C******************************************************************************
C*    ADCIRC initialize routine
C******************************************************************************
      SUBROUTINE ADCIRC_Init(COMM, ROOTD)
      !USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      USE VERSION, ONLY : adc_version, adc_hash, fileFmtMajor, fileFmtMinor, 
     &   FileFmtRev
      USE SIZES, ONLY : rootDir, mnproc, mnwproc, mnallproc, myproc,
     &   read_local_hot_start_files, make_dirname, write_local_harm_files,
     &   get_numwriters, mnstae
      USE GLOBAL, ONLY : usingDynamicWaterLevelCorrection,
     &   dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, 
     &   dynamicWaterLevelCorrectionSta, DEBUG, INFO, allMessage, itime, dtdp, statim,
     &   CHotHS, dtdpHS, FileFmtVersion, ifnlcat, nabout, ncchange,
     &   ifnlcate, ifnlct, ifnlcte, ifnlfa, ifwind, ihot, im, iths,
     &   nolica, nolicat, nolifa, nscreen, nstae, nt, screenunit,
     &   version_number, initLogging, setMessageSource, screenMessage,
     &   openLogFile, unsetMessageSource, np_global, rnp_global,
     &   scratchMessage, CTIP, L_N, writer_id
      USE NodalAttributes, ONLY: activateVEW1DChannelWetPerimeter, 
     &   initVEW1DChannelWetPerimeter
      USE ADC_CONSTANTS, ONLY: g, rhowat0
      USE GLOBAL_3DVS, ONLY : AMB, goRho, goRhoOAMB 
      USE MESH, ONLY : Areas, TotalArea, totalAreaCalc, ne, np, SFEA
      USE GWCE, ONLY : initGWCE
      USE MOMENTUM, ONLY : initMomentumEq
      USE HARM, ONLY : allocateFullDomainHAIOArrays,
     &                 initHarmonicParameters
      USE GLOBAL_IO, ONLY : allocateFullDomainIOArrays
      USE WRITE_OUTPUT, ONLY : initOutput2D
      USE WIND, ONLY : initWindModule, initDynamicWaterLevelCorrection,
     &                 getDynamicWaterLevelCorrections
#ifdef CMPI
      USE MESSENGER, ONLY : mnwproh, msg_init, msg_table, msg_start,
     &   allnodes
      USE WRITER, ONLY : writer_init, sendLabelInfoToWriters
      USE HSWRITER, ONLY : hswriter_init
#endif
      IMPLICIT NONE
C
      INTEGER, OPTIONAL :: COMM
      CHARACTER(*), OPTIONAL :: ROOTD
      REAL(8) :: StaTimHS
      REAL(8) :: RefTimHS
      REAL(8) :: TimeLoc
      INTEGER :: iargc    ! function to return the number of command line arguments
      INTEGER :: argcount ! number of command line arguments
      CHARACTER(2048) :: cmdlinearg ! a single command line argument
      INTEGER :: i        ! command line argument counter

      !
      REAL (8):: start_time = 0.D0, end_time = 0.D0

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
            select case(trim(cmdlinearg))
            case("-V","-v","--version") 
               ! write the version number string ... actually it is
               ! a series of two or three numbers with two digits 
               ! separated by dots, e.g., 51.02 or 51.52.07
               write(6,'(a)') trim(adc_version)                 
               CALL EXIT(0) 
            case("--hash")
               write(6,'(a)') trim(adc_hash) 
               CALL EXIT(0) 
            case default
               ! do nothing, this is some other command line argument
            end select
         end do
      endif
      !
      call initLogging()
      nabout = 0 !jgf52.08.01: Initialize log level to ECHO.
#ifdef EBUG
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
      if ( (mnwproh > 0.or.mnwproc > 0) .and. myProc >= mnproc) then
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
         WRITE(ScreenUnit,'(a)')"ADCIRC Commit hash is "//ADC_HASH
      ENDIF

#ifdef WTIMER
#ifdef CMPI
      start_time = MPI_Wtime()
      start_total_wtime = MPI_Wtime() 
#endif
#endif
      
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

      ! Initialize variables for bottom friction evaluation considering wet perimeter 11/06/2023
      if ( activateVEW1DChannelWetPerimeter.eqv..true. ) THEN
            CALL initVEW1DChannelWetPerimeter()
      ENDIF

      !jgf49.44: Initialize harmonic analysis parameters.
      CALL initHarmonicParameters()
C
C     jgf: Initialize default values in wind module and allocate
C     memory for internal use.
      call initWindModule()

      !jgf53.dev: allocate memory for continuity and momentum eqns
      call initGWCE()
      call initMomentumEq()

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
      
      !...Read the self attraction and loading file (fort.24)
      call readSelfAttractionAndLoadingFile()

C...************* Make the L_N Tidal potential multiple function *************
C... WJP: Feb 21 2018; moved here before hstart Feb 27 2019
C..  Corrected long term species (unfortunately there is a mistake in old ADCIRC 
C..  theory report that also propagated into some subsequent publications)
      IF (CTIP) THEN

         DO I = 1,NP
C           LONG-TERM SPECIES
            L_N(0,I) = 1.5d0*COS(SFEA(I))*COS(SFEA(I)) - 1d0
C           DIURNAL SPECIES
            L_N(1,I) = SIN(2d0*SFEA(I))
C           SEMIDIURNAL SPECIES
            L_N(2,I) = COS(SFEA(I))*COS(SFEA(I))
         END DO

      ENDIF

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

      !jgf: Initialize water level offset if needed
      if (usingDynamicWaterLevelCorrection.eqv..true.) then
         allocate(dynamicWaterLevelCorrection1(np))
         allocate(dynamicWaterLevelCorrection2(np)) ! used in timestepping
         allocate(dynamicWaterLevelCorrectionSta(mnstae))
         call initDynamicWaterLevelCorrection(timeloc) ! initializes the i/o 
         call getDynamicWaterLevelCorrections(dynamicWaterLevelCorrection2, timeloc) ! populate 
      endif
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

      GORho=G/RhoWat0
      GORhoOAMB=GORho/AMB

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

#ifdef WTIMER      
#ifdef CMPI
      end_time = MPI_Wtime()
#endif
      init_wtime = init_wtime + (end_time - start_time)  
#endif

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
      USE SIZES, ONLY : mnproc, mnwproc, mnwproh, myproc
      USE GLOBAL, ONLY : iTime, DEBUG, INFO, setMessageSource, allMessage,
     &   unsetMessageSource, statim, dtdp, screenunit, logMessage, earlyterminate
#ifdef CMPI
      USE WRITER, ONLY : writer_main, writer_pause, sendLabelInfoToWriters
      USE HSWRITER, ONLY : hswriter_main, hswriter_pause
#endif
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
      REAL(8) :: TimeLoc
C
 
      REAL(8):: start_time, end_time

      call setMessageSource("ADCIRC_Run")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

#ifdef WTIMER      
#ifdef CMPI
      start_time = MPI_Wtime() ; 
#endif      
#endif

!TCM v49.52.01 added this ifdef (brought in from adcirc_init()
#ifdef CMPI
      if (((mnwproh > 0).or.(mnwproc > 0)) .and. myproc >= mnproc) then
         if ( myproc < mnproc+mnwproc ) then   !st3
            !call sendLabelInfoToWriters()
            call logMessage(INFO,"Enter writer main.")
            call writer_main()
            call unsetMessageSource()  !tcm 20150924 -- added
            return
         else
            call hswriter_main()
            call unsetMessageSource()  !tcm 20150924 -- added
            return
         endif
      elseif ( mnwproc > 0 .and. myproc == 0 ) then
         call sendLabelInfoToWriters()
      endif
#endif

      IF (PRESENT(NTIME_STP)) THEN
        ITIME_END = ITIME_BGN+NTIME_STP-1
      ENDIF
C
      call allMessage(INFO,"Begin timestepping.")
      DO 100 ITIME=ITIME_BGN,ITIME_END
        CALL TIMESTEP(ITIME,TimeLoc)
C
#ifdef CSWAN
Casey 090302: If it is time, then call the following subroutine
C             to then call the SWAN time-stepping subroutine.
        IF(MOD(ITIME,CouplingInterval).EQ.0)THEN
          CALL PADCSWAN_RUN(ITIME)
        ENDIF
#endif
        if ( earlyterminate.ne.0 ) then 
           call writeOutput2D(ITIME_END,TimeLoc) 
           exit
        else
        !...WRITE OUTPUT
        !   zc - moved this here so that output occurs after
        !        the SWAN time step has been completed.
            CALL writeOutput2D(ITIME,TimeLoc) ! =>see write_output.F
        endif
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

#ifdef WTIMER
#ifdef CMPI
      end_time = MPI_Wtime() ; 
#endif 
      run_wtime = run_wtime + (end_time - start_time) 
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
      USE GLOBAL, ONLY : iTime, DEBUG, setMessageSource, allMessage,
     &   unsetMessageSource
      USE HARM, ONLY : solveHarmonicAnalysis
      USE WRITE_OUTPUT, ONLY : writeHarmonicAnalysisOutput
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      USE SIZES, ONLY : mnproc, mnwproc, mnwproh, myproc

      IMPLICIT NONE
C
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      INTEGER :: I

      REAL (8):: start_time, end_time
C
      call setMessageSource("ADCIRC_Final")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

#ifdef WTIMER
#ifdef CMPI
      start_time = MPI_Wtime() ; 
#endif 
#endif

C
C...
C...*************** SOLVE THE HARMONIC ANALYSIS PROBLEM ****************
C...
      CALL solveHarmonicAnalysis(ITIME)
      CALL writeHarmonicAnalysisOutput(ITIME)
C
#ifdef WTIMER      
#ifdef CMPI
      end_time = MPI_Wtime() ; 
      end_total_wtime = MPI_Wtime()
#endif     
      final_wtime = final_wtime + (end_time - start_time) 
      total_wtime = total_Wtime + 
     &              (end_total_wtime - start_total_wtime) 

      if ( myproc < mnproc) then
        !
        CALL collect_run_wall_time()
      end if  
#endif

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
C-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Final
C-----------------------------------------------------------------------

C******************************************************************************
C*    ADCIRC terminate routine
C*    jgf50.44: Mostly just needed so that we can clean up mpi when
C*    bombing out.
C******************************************************************************
      SUBROUTINE ADCIRC_Terminate(NO_MPI_FINALIZE)
      USE GLOBAL, ONLY : DEBUG, INFO, setMessageSource, allMessage, 
     &   unsetMessageSource
#ifdef CMPI
      USE MESSENGER, ONLY : subdomainFatalError, msg_fini
#endif

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
      CALL EXIT(1) 
C
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Terminate
C-----------------------------------------------------------------------

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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ADCIRC_LocalTerminate()
        USE GLOBAL,ONLY: terminate_localproc, screenUnit
        IMPLICIT NONE
        !...Set the local terminator and then go wait for the rest of
        !   the compute processors to check in
        TERMINATE_LOCALPROC = 1
        CALL ADCIRC_CheckLocalTerminate()
C-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_LocalTerminate
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ADCIRC_CheckLocalTerminate()
        USE MESSENGER
        USE GLOBAL,ONLY: terminate_localproc, screenUnit
        USE MPI
        IMPLICIT NONE
        INTEGER :: TERM_LOCAL(2)
        INTEGER :: TERM_GLOBAL(2)
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
C-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_CheckLocalTerminate
C-----------------------------------------------------------------------
#endif

C******************************************************************************
C******************************************************************************

C......................................................................
      SUBROUTINE collect_run_wall_time( )
         USE SIZES, ONLY: MYPROC, MNPROC

#ifdef CMPI         
         USE MESSENGER, ONLY : MSG_RScalar_Reduce
#endif

         IMPLICIT NONE

#ifdef WTIMER         
         REAL (8):: init_max_rank, init_min_rank
         REAL (8):: init_max_wall_time, init_min_wall_time, init_sum_wall_time

         REAL (8):: run_max_rank, run_min_rank
         REAL (8):: run_max_wall_time, run_min_wall_time, run_sum_wall_time

         REAL (8):: final_max_rank, final_min_rank
         REAL (8):: final_max_wall_time, final_min_wall_time, final_sum_wall_time

         REAL (8):: total_max_rank, total_min_rank
         REAL (8):: total_max_wall_time, total_min_wall_time, total_sum_wall_time

#ifdef CMPI
         init_max_rank = DBLE(MYPROC)
         init_min_rank = DBLE(MYPROC) 
         
         run_max_rank = DBLE(MYPROC)
         run_min_rank = DBLE(MYPROC)

         final_max_rank = DBLE(MYPROC)
         final_min_rank = DBLE(MYPROC)

         total_max_rank = DBLE(MYPROC)
         total_min_rank = DBLE(MYPROC)

         ! init
         CALL MSG_RScalar_Reduce( 'MAX', init_wtime, 
     &        init_max_wall_time, init_max_rank )

         CALL MSG_RScalar_Reduce( 'MIN', init_wtime, 
     &        init_min_wall_time, init_min_rank )

         CALL MSG_RScalar_Reduce( 'SUM', init_wtime,
     &        init_sum_wall_time ) ;

         ! run 
         CALL MSG_RScalar_Reduce( 'MAX', run_wtime, 
     &        run_max_wall_time, run_max_rank )

         CALL MSG_RScalar_Reduce( 'MIN', run_wtime, 
     &        run_min_wall_time, run_min_rank )

         CALL MSG_RScalar_Reduce( 'SUM', run_wtime,
     &        run_sum_wall_time ) ;

         ! finalize
         CALL MSG_RScalar_Reduce( 'MAX', final_wtime, 
     &        final_max_wall_time, final_max_rank )

         CALL MSG_RScalar_Reduce( 'MIN', final_wtime, 
     &        final_min_wall_time, final_min_rank )

         CALL MSG_RScalar_Reduce( 'SUM', final_wtime,
     &        final_sum_wall_time ) ;


         ! Total
         CALL MSG_RScalar_Reduce( 'MAX', total_wtime, 
     &        total_max_wall_time, total_max_rank )

         CALL MSG_RScalar_Reduce( 'MIN', total_wtime, 
     &        total_min_wall_time, total_min_rank )

         CALL MSG_RScalar_Reduce( 'SUM', total_wtime,
     &        total_sum_wall_time ) ;


         IF ( MYPROC == 0 ) THEN
          OPEN( UNIT = 830, FILE = 'runtime_report.out', 
     &          ACTION="write") ;

          !
          WRITE( 830, '(A)' ) '% Time (in second) spent in the driver'//
     &            ' subroutines.'
          WRITE( 830, '(A)') ' max_time, max_rank, min_time,'// 
     &                ' min_rank, average_time'

          WRITE( 830, '(A)') '%% ADCIRC_Init : '  
          WRITE( 830, '(F19.5,1X,I6,1X,F19.5,I6,1X,F19.5)' )
     &       init_max_wall_time, INT(init_max_rank), init_min_wall_time,
     &       INT(init_min_rank), init_sum_wall_time/DBLE(MNPROC) 

          WRITE( 830, '(A)') '%% ADCIRC_run : ' 
          WRITE( 830, '(F19.5,1X,I6,1X,F19.5,I6,1X,F19.5)' )
     &       run_max_wall_time, INT(run_max_rank), run_min_wall_time,
     &       INT(run_min_rank), run_sum_wall_time/DBLE(MNPROC)

          WRITE( 830, '(A)') '%% ADCIRC_final : ' 
          WRITE( 830, '(F19.5,1X,I6,1X,F19.5,I6,1X,F19.5)' )
     &       final_max_wall_time, INT(final_max_rank), final_min_wall_time,
     &       INT(final_min_rank), final_sum_wall_time/DBLE(MNPROC)

          WRITE( 830, '(A)') '%% Total :' 
          WRITE( 830, '(F19.5,1X,I6,1X,F19.5,I6,1X,F19.5)' )
     &       total_max_wall_time, INT(total_max_rank), total_min_wall_time,
     &       INT(total_min_rank), total_sum_wall_time/DBLE(MNPROC)


          CLOSE(830) ;
         END IF     

#endif
#endif

      END SUBROUTINE collect_run_wall_time
C----------------------------------------------------------------------

      END MODULE ADCIRC_Mod
C******************************************************************************
C******************************************************************************

