!-------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation, either version 3 of the License, or
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
<<<<<<< HEAD
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
#include "logging_macros.h"
=======
!-----------------------------------------------------------------------
! (P)ADCIRC VERSION 45.12 03/17/2006
!  last changes in this file VERSION 45.12
!
!-----------------------------------------------------------------------
!
!                           (P)ADCIRC
!
!    A (PARALLEL) ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL SEAS
!                         AND ESTUARIES
!
!
!                          DEVELOPED BY:
!
!                      DR. R.A. LUETTICH, JR
!
!             UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL
!                   INSTITUTE OF MARINE SCIENCES
!                       3431 ARENDELL ST.
!                   MOREHEAD CITY, NC, 28557
!                         252-726-6841 EXT. 137
!                   EMAIL  RICK_LUETTICH@UNC.EDU
!
!                        DR. J.J. WESTERINK
!
!          DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES
!                     UNIVERSITY OF NOTRE DAME
!                       NOTRE DAME, IN 46556
!                         219-631-6475
!                    EMAIL JJW@PHOTIUS.CE.ND.EDU
!
!
!        MAJOR FUNDING FOR THE DEVELOPMENT OF ADCIRC WAS PROVIDED BY
!
!                       DEPARTMENT OF THE ARMY
!                    WATERWAYS EXPERIMENT STATION
!                 COASTAL ENGINEERING RESEARCH CENTER
!                        3909 HALLS FERRY RD
!                      VICKSBURG, MI 39180-6199
!
!-----------------------------------------------------------------------
!
!          THE ADCIRC SOURCE CODE IS COPYRIGHTED, 1994-2006 BY:
!
!                 R.A. LUETTICH, JR AND J.J. WESTERINK
!
!         NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED
!               WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS
!
!-----------------------------------------------------------------------
!
!                  Parallelization of ADCIRC 2D was done by
!                    the center for subsurface modeling
!                         The University of Texas
!                             Austin, TX 78712
!                            03/06/98 - 5/19/99
!
!                   Translation of 2D code to Fortran90
!                             Victor J. Parr
!                              John B. Romo
!                               8/31/99
!
!               Parallelization, consolidation of 2D & 3D codes,
!             translation of 3D routines to Fortran 90 was done by
!                               Tim Campbell
!                   Naval Research Lab, Stennis Space Center
!                             summer 2002
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      MODULE ADCIRC_Mod
!-----------------------------------------------------------------------
#ifdef CMPI
      use mpi_f08, only: MPI_Comm
#endif

      IMPLICIT NONE
      PUBLIC

      INTEGER, SAVE :: ITIME_BGN, ITIME_END


      CONTAINS

<<<<<<< HEAD
C******************************************************************************
C*    ADCIRC initialize routine
C******************************************************************************
=======

!------------------------------------------------------------------------
!> @brief ADCIRC initialization routine.
!------------------------------------------------------------------------
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      SUBROUTINE ADCIRC_Init(COMM, ROOTD)
!------------------------------------------------------------------------
      !USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      USE VERSION, ONLY : adc_version, adc_hash, fileFmtMajor,
     &   fileFmtMinor, FileFmtRev
      USE SIZES, ONLY : rootDir, mnproc, mnwproc, mnallproc, myproc,
<<<<<<< HEAD
     &   read_local_hot_start_files, make_dirname, write_local_harm_files,
     &   get_numwriters, mnstae, mnwproh, mnwproc, localdir
      USE GLOBAL, ONLY : usingDynamicWaterLevelCorrection,
     &   dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2,
     &   dynamicWaterLevelCorrectionSta, itime, dtdp, statim,
     &   CHotHS, dtdpHS, FileFmtVersion, ifnlcat, ncchange,
     &   ifnlcate, ifnlct, ifnlcte, ifnlfa, ifwind, ihot, im, iths,
     &   nolica, nolicat, nolifa, nstae, nt,
     &   version_number, np_global, rnp_global,
     &   CTIP, L_N, writer_id, dynamicWaterLevelCorrectionDelta
      USE mod_logging, ONLY : allMessage, screenMessage,
     &   openLogFile, DEBUG, INFO, screenUnit,
     &   nabout, nscreen, ECHO, t_log_level, t_log_scope, init_log_scope
=======
     &   read_local_hot_start_files, make_dirname, get_numwriters,
     &   write_local_harm_files, mnstae, mnwproh, mnwproc
      USE GLOBAL, ONLY : usingDynamicWaterLevelCorrection,
     &   dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, 
     &   dynamicWaterLevelCorrectionSta, DEBUG, INFO, allMessage, itime,
     &   dtdp, statim, CHotHS, dtdpHS, FileFmtVersion, ifnlcat, nabout,
     &   ncchange, ifnlcate, ifnlct, ifnlcte, ifnlfa, ifwind, ihot, im,
     &   iths, nolica, nolicat, nolifa, nscreen, nstae, nt, screenunit,
     &   version_number, initLogging, setMessageSource, screenMessage,
     &   scratchMessage, openLogFile, unsetMessageSource, np_global,
     &   rnp_global, CTIP, L_N, writer_id,
     &   dynamicWaterLevelCorrectionDelta
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      USE NodalAttributes, ONLY: activateVEW1DChannelWetPerimeter, 
     &   initVEW1DChannelWetPerimeter
      use mod_terminate, only: terminate, ADCIRC_EXIT_SUCCESS
      USE ADC_CONSTANTS, ONLY: g, rhowat0
<<<<<<< HEAD
      USE GLOBAL_3DVS, ONLY : AMB, goRho, goRhoOAMB,
     &   NSTA3DD, NSTA3DV, NSTA3DT, NSTA3DD_G, NSTA3DV_G, NSTA3DT_G,
     &   IMAP_STA3DD_LG, IMAP_STA3DV_LG, IMAP_STA3DT_LG 
      USE MESH, ONLY : Areas, TotalArea, totalAreaCalc, ne, np, SFEA
=======
      USE GLOBAL_3DVS, ONLY : AMB, goRho, goRhoOAMB 
      USE MESH, ONLY : totalAreaCalc, ne, np, SFEA
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      USE GWCE, ONLY : initGWCE
      USE MOMENTUM, ONLY : initMomentumEq
      USE HARM, ONLY : allocateFullDomainHAIOArrays,
     &                 initHarmonicParameters
      USE GLOBAL_IO, ONLY : allocateFullDomainIOArrays
      USE WRITE_OUTPUT, ONLY : initOutput2D
      USE WIND, ONLY : initWindModule, initDynamicWaterLevelCorrection,
     &                 getDynamicWaterLevelCorrections
#ifdef CMPI
      USE MESSENGER, ONLY : msg_init, msg_table, msg_start,
     &   allnodes
      USE WRITER, ONLY : writer_init, sendLabelInfoToWriters
      USE HSWRITER, ONLY : hswriter_init
#endif
      IMPLICIT NONE

<<<<<<< HEAD
#ifdef CMPI
      type(MPI_Comm), OPTIONAL :: COMM
#else
=======
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      INTEGER, OPTIONAL :: COMM
#endif

      CHARACTER(*), OPTIONAL :: ROOTD
      REAL(8) :: StaTimHS
      REAL(8) :: RefTimHS
      REAL(8) :: TimeLoc
      INTEGER :: iargc    ! function to get number command line args
      INTEGER :: argcount ! number of command line arguments
      CHARACTER(2048) :: cmdlinearg ! a single command line argument
      character(1024) :: scratchMessage
      INTEGER :: i        ! command line argument counter
      REAL (8) :: start_time = 0.D0, end_time = 0.D0

      LOG_SCOPE_TRACED("ADCIRC_Init", ADCIRC_TRACING)

      ! jgf51.52.07: Added command line option to just write out
      ! the ADCIRC version number and exit. This is to allow 
      ! external applications to query the version number so that
      ! they can behave in version-appropriate ways. 
      argcount  = COMMAND_ARGUMENT_COUNT()  ! get number of command line arguments
      if ((argcount.gt.0).and.(myproc.eq.0)) then
         i = 0
         do while (i.lt.argcount)
            i = i + 1
            call GET_COMMAND_ARGUMENT(i, cmdlinearg)
            select case(trim(cmdlinearg))
            case("-h","--help")
               call PRINT_HELP()
               call exit(ADCIRC_EXIT_SUCCESS)
            case("-V","-v","--version")
               ! write the version number string ... actually it is
               ! a series of two or three numbers with two digits
               ! separated by dots, e.g., 51.02 or 51.52.07
               write(*,'(a)') trim(adc_version)
               call exit(ADCIRC_EXIT_SUCCESS)
            case("--hash")
               write(*,'(a)') trim(adc_hash)
               call exit(ADCIRC_EXIT_SUCCESS)
            case default
               ! do nothing, this is some other command line argument
            end select
         end do
      endif
<<<<<<< HEAD

      nabout = ECHO
#ifdef EBUG
      nabout = DEBUG
=======

      call initLogging()
      nabout = 0 !jgf52.08.01: Initialize log level to ECHO.
#ifdef DEBUG
      nabout = -1 !jgf52.08.01: Initialize log level to DEBUG
#endif
      call setMessageSource("ADCIRC_Init")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call screenMessage(DEBUG,"Enter.") ! log to screen
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
#endif
      IF (PRESENT(ROOTD)) THEN
        ROOTDIR = TRIM(ROOTD)
      ELSE
        ROOTDIR = '.'
      ENDIF

#ifdef CMPI
      CALL GET_NUMWRITERS()  ! Get number of writer processors
      IF (PRESENT(COMM)) THEN
        CALL MSG_INIT(COMM)  ! Init MPI and get MPI-rank of this cpu
      ELSE
        CALL MSG_INIT()      ! Init MPI and get MPI-rank of this cpu
      ENDIF
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
<<<<<<< HEAD
      call openLogFile(localdir)  ! jgf50.65: open subdomain log files (fort.16)
=======
      call openLogFile()     ! jgf50.65: open subdomain log (fort.16)
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      CALL WRITER_INIT()     ! Initialize WRITER module
      CALL HSWRITER_INIT()   ! Initialize HSWRITER module: hstart 100711
      if ( (mnwproh > 0.or.mnwproc > 0) .and. myProc >= mnproc) then
         return
      endif
<<<<<<< HEAD
      CALL READ_INPUT()      ! Establish sizes by reading fort.14 and fort.15
      CALL MSG_TABLE(NSTA3DD, NSTA3DV, NSTA3DT,       ! Read Message-Passing Tables
     &               NSTA3DD_G, NSTA3DV_G, NSTA3DT_G,
     &               IMAP_STA3DD_LG, IMAP_STA3DV_LG, IMAP_STA3DT_LG)
=======
      CALL READ_INPUT()      ! Establish sizes: read fort.14 and fort.15
      CALL MSG_TABLE()       ! Read Message-Passing Tables
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      CALL MSG_START()       ! Startup message passing
#else
      MNPROC = 1             ! Init number of compute procs for serial
      MNWPROC = 0            ! Init number of writer procs for serial
      MNALLPROC = 1          ! Init number of all procs for serial
      MYPROC = 0             ! Init task id for serial
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
      call openLogFile()     ! jgf50.65: open fort.16 file
      CALL READ_INPUT()      ! Establish sizes: read fort.14 and fort.15
#endif

      ! File Fmt Version  must match Major and Minor numbers.
      FileFmtVersion =
     &   VERSION_NUMBER(FileFmtMajor, FileFmtMinor, FileFmtRev)
      IF ((NSCREEN.NE.0).AND.(MYPROC.EQ.0)) THEN
         WRITE(ScreenUnit,'(a)')"ADCIRC Version is "//ADC_VERSION
         WRITE(ScreenUnit,'(a)')"ADCIRC Commit hash is "//ADC_HASH
      ENDIF

      
      !jgf49.44: Allocate memory for full domain array I/O operations
      IF ( (MNPROC.gt.1).and.
     &        ( (myProc.eq.0)
     &          .or.(READ_LOCAL_HOT_START_FILES.eqv..false.) ) ) THEN
          CALL allocateFullDomainIOArrays()
      ENDIF
      IF ( (MNPROC.gt.1).and.(myProc.eq.0).and.
     &     (WRITE_LOCAL_HARM_FILES.eqv..false.) ) THEN
          CALL allocateFullDomainHAIOArrays()
      ENDIF

      ! Initialize variables for bottom friction (wet perimeter) 11/06/23
      if ( activateVEW1DChannelWetPerimeter.eqv..true. ) THEN
            CALL initVEW1DChannelWetPerimeter()
      ENDIF

      !jgf49.44: Initialize harmonic analysis parameters.
      CALL initHarmonicParameters()

<<<<<<< HEAD
C     jgf: Initialize default values in wind module and allocate
C     memory for internal use.
=======
      !jgf: Initialize default values in wind module and allocate
      ! memory for internal use.
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      call initWindModule()

      !jgf53.dev: allocate memory for continuity and momentum eqns
      call initGWCE()
      call initMomentumEq()

<<<<<<< HEAD
C...Compute the reciprocal of the number of nodes in the entire domain

=======
#ifdef CMACHSUN
      CALL ABRUPT_UNDERFLOW() ! Abrupt underflow for Sun computers
#endif

      ! Compute the reciprocal of the number of nodes in the mesh
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
#ifdef CMPI
      CALL ALLNODES(NP_GLOBAL)
      RNP_GLOBAL = 1.0D0/NP_GLOBAL
#endif

      !jgf45.08 put this here so it can be used in hot start
      IF(NOLIFA.EQ.0) THEN
         IFNLFA=0
      ELSE
         IFNLFA=1
      ENDIF
      
      !...Read the self attraction and loading file (fort.24)
      call readSelfAttractionAndLoadingFile()

      !--------- Make the L_N Tidal potential multiple function ---------
      !WJP: Feb 21 2018; moved here before hstart Feb 27 2019
      ! Corrected long term species (unfortunately there is a mistake in
      ! old ADCIRC theory report that also propagated into some
      ! publications)
      IF (CTIP) THEN

         DO I = 1,NP
            ! LONG-TERM SPECIES
            L_N(0,I) = 1.5d0*COS(SFEA(I))*COS(SFEA(I)) - 1d0
            ! DIURNAL SPECIES
            L_N(1,I) = SIN(2d0*SFEA(I))
            ! SEMIDIURNAL SPECIES
            L_N(2,I) = COS(SFEA(I))*COS(SFEA(I))
         END DO

      ENDIF

      !-----------------------------------------------------------------
      !                  PROGRAM SETUP SECTION
      !-----------------------------------------------------------------

      SELECT CASE(IHOT)
         CASE(0)
            CALL COLDSTART()
            ITHS = 0  !TCM20091029 v49.02 Initialize to zero for cstart
         ! kmd : added in the initial condition read
         CASE(17,67,68,367,368,567,568) ! read the hotstart file
            CALL HOTSTART()
         CASE DEFAULT
            ! the IHOT value was vetted already, so this code should
            ! not be reachable
            IF ((NSCREEN.ne.0).and.(myProc.eq.0)) THEN
               WRITE(screenUnit,*) "ERROR: IHOT=",IHOT," is not valid."
            ENDIF
      END SELECT

      ! Default beginning and ending time steps
      ITIME_BGN = ITHS+1
      ITIME_END = NT  !NT is set in read_input.F

      ! Determine the number of active elements (MJU), the total number
      ! of elements (NODELE) and the total area of elements (TotalArea)
      ! attached to each node.
      ! CMS 11/20/25: moved inside cstart/hstart modules after nnodecode
      ! is updated, since it is needed for the hydrology initial setup.
      ! Removed areas and totalarea from use mesh statement above.

      ! call totalAreaCalc()
 
      ! jgf52.08.20: Compute master time which is referenced to the
      ! beginning time of the model run. This is needed here only
      ! because we need to determine whether to load min/max files in
      ! subroutine initOutput2D().
      !
      ! If the time in the min/max file is later than the current time
      ! in ADCIRC, we don't load the min/max files. TODO: This is
      ! computed again in subroutine timestep().
      TimeLoc=ITIME_BGN*DTDP + StaTim*86400.D0
      IF (CHOTHS.eqv..true.) THEN
         StaTimHS=((ITIME_BGN-1)*DTDPHS)/86400.D0
         RefTimHS=((ITIME_BGN-1)*DTDP)/86400.D0
         TimeLoc=ITIME_BGN*DTDP + (StaTimHS - RefTimHS)*86400.D0
      END IF

      !jgf: Initialize water level offset if needed
      ! time difference in offset: tgaContOffset
      allocate(dynamicWaterLevelCorrectionDelta(np))
      dynamicWaterLevelCorrectionDelta=0.D0
      if (usingDynamicWaterLevelCorrection.eqv..true.) then
         allocate(dynamicWaterLevelCorrection1(np))
         allocate(dynamicWaterLevelCorrection2(np)) ! for timestepping
         allocate(dynamicWaterLevelCorrectionSta(mnstae))
         call initDynamicWaterLevelCorrection(timeloc) ! init the i/o
         call getDynamicWaterLevelCorrections(
     &        dynamicWaterLevelCorrection2, timeloc) ! populate
      endif


<<<<<<< HEAD
C     jgf51.21.24: Initialize output.
=======
      !jgf51.21.24: Initialize output.
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      call initOutput2D(timeloc)


      !-----------------------------------------------------------------
      !        SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING
      !-----------------------------------------------------------------

      ! Set flags for nonlinear terms
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

      ! Corbitt 1203022: Local Advection Flags
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

 1999 format(1X, a," ADCIRC Version: ",a)



      RETURN
!------------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Init
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!> @brief ADCIRC run routine.
!------------------------------------------------------------------------
      SUBROUTINE ADCIRC_Run(NTIME_STP)
!------------------------------------------------------------------------
      USE SIZES, ONLY : mnproc, mnwproc, mnwproh, myproc
<<<<<<< HEAD
      USE GLOBAL, ONLY : iTime, statim, dtdp,
     &   earlyterminate
      USE mod_logging, ONLY : allMessage,
     &   logMessage, DEBUG, INFO, t_log_scope, init_log_scope
=======
      USE GLOBAL, ONLY : iTime, DEBUG, INFO, setMessageSource,
     &    allMessage, unsetMessageSource, statim, dtdp, screenunit,
     &    logMessage, earlyterminate
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
#ifdef CMPI
      USE WRITER, ONLY : writer_main, writer_pause,
     &    sendLabelInfoToWriters
      USE HSWRITER, ONLY : hswriter_main, hswriter_pause
#endif

#ifdef CSWAN
      !Casey 090302: for coupling to unstructured SWAN
      USE Couple2Swan, ONLY: CouplingInterval, PADCSWAN_RUN
#endif
      USE WRITE_OUTPUT,ONLY: writeOutput2D
#ifdef CMPI      
      USE MESSENGER, ONLY: writers_active, hs_writers_active
#endif      
      IMPLICIT NONE

      INTEGER, OPTIONAL :: NTIME_STP
      REAL(8) :: TimeLoc


      REAL(8):: start_time, end_time

      LOG_SCOPE_TRACED("ADCIRC_Run", ADCIRC_TRACING)

      !TCM v49.52.01 added #ifdef test from adcirc_init()
#ifdef CMPI
      if (((mnwproh > 0).or.(mnwproc > 0)) .and. myproc >= mnproc) then
         if ( myproc < mnproc+mnwproc ) then   !st3
            !call sendLabelInfoToWriters()
            call logMessage(INFO,"Enter writer main.")
            call writer_main()
            return
         else
            call hswriter_main()
            return
         endif
      elseif ( (mnwproc > 0 .or. mnwproh > 0) .and. myproc == 0 ) then
         if(mnwproc > 0)then
             call sendLabelInfoToWriters()
             writers_active = .true.
         endif
         if(mnwproh > 0)then
             hs_writers_active = .true.
         endif
      endif
#endif

      IF (PRESENT(NTIME_STP)) THEN
        ITIME_END = ITIME_BGN+NTIME_STP-1
      ENDIF

      call allMessage(INFO,"Begin timestepping.")
      DO ITIME=ITIME_BGN,ITIME_END
        CALL TIMESTEP(ITIME,TimeLoc)

#ifdef CSWAN
        !Casey 090302: For coupled swan if it is time, call SWAN
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
<<<<<<< HEAD
      END DO
=======
100   CONTINUE
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart

      IF (PRESENT(NTIME_STP)) THEN
        ITIME_BGN = ITIME_END+1
      ENDIF
<<<<<<< HEAD

=======
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart


      !TCM v49.52.01 added #ifdef to pause writer processors before
      ! exiting adcirc_run phase.  This is to allow them to be used
      ! again if adcirc_run is called more than once.
#ifdef CMPI
      IF (MYPROC == 0 .and. MNWPROC > 0) CALL WRITER_PAUSE()
      IF (MYPROC == 0 .and. MNWPROH > 0) CALL HSWRITER_PAUSE()
#endif

      RETURN
!------------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Run
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!> @brief ADCIRC finalize routine.
!------------------------------------------------------------------------
      SUBROUTINE ADCIRC_Final(NO_MPI_FINALIZE)
<<<<<<< HEAD
      USE GLOBAL, ONLY : iTime
      USE mod_logging, ONLY : allMessage,
     &   DEBUG, t_log_scope, init_log_scope
=======
!------------------------------------------------------------------------
      USE GLOBAL, ONLY : iTime, DEBUG, setMessageSource, allMessage,
     &   unsetMessageSource
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      USE HARM, ONLY : solveHarmonicAnalysis
      USE WRITE_OUTPUT, ONLY : writeHarmonicAnalysisOutput
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      USE SIZES, ONLY : mnproc, mnwproc, mnwproh, myproc
      use mod_terminate, only: terminate, ADCIRC_EXIT_SUCCESS

      IMPLICIT NONE

      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
<<<<<<< HEAD

=======
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      INTEGER :: I
      REAL (8):: start_time, end_time
<<<<<<< HEAD

      LOG_SCOPE_TRACED("ADCIRC_Final", ADCIRC_TRACING)
C...
C...*************** SOLVE THE HARMONIC ANALYSIS PROBLEM ****************
C...
=======

      call setMessageSource("ADCIRC_Final")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


      ! Solve and write output for harmonic analysis
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      CALL solveHarmonicAnalysis(ITIME)
      CALL writeHarmonicAnalysisOutput(ITIME)

      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL terminate(exit_code=ADCIRC_EXIT_SUCCESS,
     &                 do_finalize_mpi=.not.NO_MPI_FINALIZE,
     &                 message="ADCIRC terminating normally")
      ELSE
        CALL terminate(exit_code=ADCIRC_EXIT_SUCCESS,
     &                 do_finalize_mpi=.true.,
     &                 message="ADCIRC terminating normally")
      ENDIF
<<<<<<< HEAD

=======
#endif

#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Final
!-----------------------------------------------------------------------


<<<<<<< HEAD
C-----------------------------------------------------------------------
C     SUBROUTINE PRINT_HELP
C-----------------------------------------------------------------------
C     Print detailed help information for adcirc/padcirc command-line
C     options
C-----------------------------------------------------------------------
      SUBROUTINE PRINT_HELP()
      USE VERSION, ONLY: adc_version
      IMPLICIT NONE
=======
!-----------------------------------------------------------------------
!> @brief jgf50.44: Exit mpi cleanly when terminate unexpectedly.
!-----------------------------------------------------------------------
      SUBROUTINE ADCIRC_Terminate(NO_MPI_FINALIZE)
!-----------------------------------------------------------------------
      USE GLOBAL, ONLY : DEBUG, INFO, setMessageSource, allMessage, 
     &   unsetMessageSource
#ifdef CMPI
      USE MESSENGER, ONLY : subdomainFatalError, msg_fini
#endif

      IMPLICIT NONE

      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE

      call setMessageSource("ADCIRC_Terminate")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart

#ifdef CSWAN
#ifdef CMPI
      WRITE(*,'(A)') 'Usage: padcswan [OPTION]...'
#else
      WRITE(*,'(A)') 'Usage: adcswan [OPTION]...'
#endif
<<<<<<< HEAD
      WRITE(*,'(A)') 'ADvanced CIRCulation model coupled with SWAN for shelves, coastal seas, and estuaries.'
#else
#ifdef CMPI
      WRITE(*,'(A)') 'Usage: padcirc [OPTION]...'
#else
      WRITE(*,'(A)') 'Usage: adcirc [OPTION]...'
#endif
      WRITE(*,'(A)') 'ADvanced CIRCulation model for shelves, coastal seas, and estuaries.'
#endif

      WRITE(*,'(A)') ''
      WRITE(*,'(A)') 'General options:'
      WRITE(*,'(A)') '  -h, --help             display this help and exit'
      WRITE(*,'(A)') '  -v, --version          output version information and exit'
      WRITE(*,'(A)') '      --hash             output git commit hash and exit'
      WRITE(*,'(A)') ''
#ifdef CMPI
      WRITE(*,'(A)') 'Parallel execution options:'
      WRITE(*,'(A)') '  -W NUM                 number of dedicated writer processors for output files'
      WRITE(*,'(A)') '  -Ws NUM                number of dedicated writer processors with split file output'
      WRITE(*,'(A)') '  -H NUM                 number of dedicated writer processors for hotstart files'
      WRITE(*,'(A)') ''
#endif
      WRITE(*,'(A)') 'Report bugs to <https://github.com/adcirc/adcirc/issues>.'
      WRITE(*,'(A)') 'ADCIRC home page: <https://adcirc.org>.'
      WRITE(*,'(A)') 'User documentation: <https://adcirc.github.io/adcirc/>.'

      RETURN
      END SUBROUTINE PRINT_HELP
C-----------------------------------------------------------------------
=======
      CALL EXIT(1) 

#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_Terminate
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief subroutines to check for errors on other processors.
!>
!> @details
!>  ADCIRC Local Terminate Routine and corresponding collection routine
!>  zc50.81 In the case where there is an error on a single processor,
!>  allow everyone to check if we need to bomb out. For example, mesh
!>  errors might only occur on a single processor and hang the code.
!>  Care must be taken not to inadvertantly hang the code by allowing
!>  for cases where not all compute processors pass through the same
!>  call to this routine.
!-----------------------------------------------------------------------
#ifdef CMPI
!-----------------------------------------------------------------------
      SUBROUTINE ADCIRC_LocalTerminate()
!-----------------------------------------------------------------------
        USE GLOBAL,ONLY: terminate_localproc, screenUnit
        IMPLICIT NONE
        !...Set the local terminator and then go wait for the rest of
        !   the compute processors to check in
        TERMINATE_LOCALPROC = 1
        CALL ADCIRC_CheckLocalTerminate()
!-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_LocalTerminate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE ADCIRC_CheckLocalTerminate()
!-----------------------------------------------------------------------
        USE MESSENGER
        USE GLOBAL,ONLY: terminate_localproc, screenUnit, comm
        USE SIZES,ONLY: MYPROC
        USE MPI
        IMPLICIT NONE

        INTEGER :: TERM_LOCAL(2)
        INTEGER :: TERM_GLOBAL(2)
        INTEGER :: ierr

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
!-----------------------------------------------------------------------
      END SUBROUTINE ADCIRC_CheckLocalTerminate
!-----------------------------------------------------------------------
#endif
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart


<<<<<<< HEAD
=======

!-----------------------------------------------------------------------
>>>>>>> Subject: Edit adcirc_init to move totalareacalc call to c/hstart
      END MODULE ADCIRC_Mod
!-----------------------------------------------------------------------
