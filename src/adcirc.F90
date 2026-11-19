!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
!******************************************************************************
! (P)ADCIRC VERSION 45.12 03/17/2006                                          *
!  last changes in this file VERSION 45.12                                    *
!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!                           (P)ADCIRC                                         *
!                                                                             *
!    A (PARALLEL) ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL SEAS        *
!                         AND ESTUARIES                                       *
!                                                                             *
!                                                                             *
!                          DEVELOPED BY:                                      *
!                                                                             *
!                      DR. R.A. LUETTICH, JR                                  *
!                                                                             *
!             UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL                     *
!                   INSTITUTE OF MARINE SCIENCES                              *
!                       3431 ARENDELL ST.                                     *
!                   MOREHEAD CITY, NC, 28557                                  *
!                         252-726-6841 EXT. 137                               *
!                   EMAIL  RICK_LUETTICH@UNC.EDU                              *
!                                                                             *
!                        DR. J.J. WESTERINK                                   *
!                                                                             *
!          DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES            *
!                     UNIVERSITY OF NOTRE DAME                                *
!                       NOTRE DAME, IN 46556                                  *
!                         219-631-6475                                        *
!                    EMAIL JJW@PHOTIUS.CE.ND.EDU                              *
!                                                                             *
!                                                                             *
!        MAJOR FUNDING FOR THE DEVELOPMENT OF ADCIRC WAS PROVIDED BY          *
!                                                                             *
!                       DEPARTMENT OF THE ARMY                                *
!                    WATERWAYS EXPERIMENT STATION                             *
!                 COASTAL ENGINEERING RESEARCH CENTER                         *
!                        3909 HALLS FERRY RD                                  *
!                      VICKSBURG, MI 39180-6199                               *
!                                                                             *
!******************************************************************************
!                                                                             *
!          THE ADCIRC SOURCE CODE IS COPYRIGHTED, 1994-2006 BY:               *
!                                                                             *
!                 R.A. LUETTICH, JR AND J.J. WESTERINK                        *
!                                                                             *
!         NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED             *
!               WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS                 *
!                                                                             *
!******************************************************************************
!                                                                             *
!                  Parallelization of ADCIRC 2D was done by                   *
!                    the center for subsurface modeling                       *
!                         The University of Texas                             *
!                             Austin, TX 78712                                *
!                            03/06/98 - 5/19/99                               *
!                                                                             *
!                   Translation of 2D code to Fortran90                       *
!                             Victor J. Parr                                  *
!                              John B. Romo                                   *
!                               8/31/99                                       *
!                                                                             *
!               Parallelization, consolidation of 2D & 3D codes,              *
!             translation of 3D routines to Fortran 90 was done by            *
!                               Tim Campbell                                  *
!                   Naval Research Lab, Stennis Space Center                  *
!                             summer 2002                                     *
!                                                                             *
!******************************************************************************

!******************************************************************************
!******************************************************************************
#include "logging_macros.h"
module ADCIRC_Mod
   !******************************************************************************
   !******************************************************************************
#ifdef CMPI
   use mpi_f08, only: MPI_Comm
#endif

   implicit none

   private

   public :: ADCIRC_Init, ADCIRC_Run, ADCIRC_Final

   ! Setup routines that live outside any module (read_input.F, cstart.F,
   ! hstart.F)
   interface
      subroutine READ_INPUT()
         implicit none
      end subroutine READ_INPUT
      subroutine readSelfAttractionAndLoadingFile()
         implicit none
      end subroutine readSelfAttractionAndLoadingFile
      subroutine COLDSTART()
         implicit none
      end subroutine COLDSTART
      subroutine HOTSTART()
         implicit none
      end subroutine HOTSTART
   end interface

   integer :: ITIME_BGN, ITIME_END

contains

   !******************************************************************************
   !*    ADCIRC initialize routine
   !******************************************************************************
   subroutine ADCIRC_Init(COMM, ROOTD)
      !USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      use VERSION, only: adc_version, adc_hash, fileFmtMajor, fileFmtMinor, &
                         FileFmtRev
      use MOD_SIZES, only: rootDir, mnproc, mnwproc, myproc, &
                           read_local_hot_start_files, make_dirname, write_local_harm_files, &
                           get_numwriters, mnstae
      use GLOBAL, only: usingDynamicWaterLevelCorrection, &
                        dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, &
                        dynamicWaterLevelCorrectionSta, dtdp, statim, &
                        CHotHS, dtdpHS, FileFmtVersion, ifnlcat, ncchange, &
                        ifnlcate, ifnlct, ifnlcte, ifnltg, ifnlfa, ifwind, ihot, im, &
                        iths, &
                        nolica, nolicat, nolifa, nt, &
                        version_number, &
                        CTIP, L_N, dynamicWaterLevelCorrectionDelta
      use mod_logging, only: allMessage, screenMessage, &
                             openLogFile, screenUnit, &
                             nabout, nscreen, ECHO, t_log_level, t_log_scope, init_log_scope
      use NodalAttributes, only: activateVEW1DChannelWetPerimeter, &
                                 initVEW1DChannelWetPerimeter
      use mod_terminate, only: terminate, ADCIRC_EXIT_SUCCESS
      use ADC_CONSTANTS, only: g, rhowat0
      use GLOBAL_3DVS, only: AMB, goRho, goRhoOAMB
      use MESH, only: totalAreaCalc, np, SFEA
      use GWCE, only: initGWCE
      use MOMENTUM2D, only: init2DMomentumEq
      use MOD_HARMONIC_ANALYSIS, only: allocateFullDomainHAIOArrays, &
                                       initHarmonicParameters
      use GLOBAL_IO, only: allocateFullDomainIOArrays
      use WRITE_OUTPUT, only: initOutput2D
      use WIND, only: initWindModule, initDynamicWaterLevelCorrection, &
                      getDynamicWaterLevelCorrections
#ifdef CMPI
      use MESSENGER, only: msg_init, msg_table, msg_start, &
                           allnodes
      use WRITER, only: writer_init, sendLabelInfoToWriters
      use MOD_SIZES, only: localdir, hotstartWritersRequested
      use GLOBAL, only: np_global, rnp_global
      use mod_terminate, only: ADCIRC_EXIT_FAILURE
      use GLOBAL_3DVS, only: NSTA3DD, NSTA3DV, NSTA3DT, NSTA3DD_G, NSTA3DV_G, &
                             NSTA3DT_G, IMAP_STA3DD_LG, IMAP_STA3DV_LG, IMAP_STA3DT_LG
#else
      use MOD_SIZES, only: mnallproc
#endif
      implicit none

#ifdef CMPI
      type(MPI_Comm), optional, intent(in) :: COMM
#else
      integer, optional, intent(in) :: COMM
#endif

      character(*), optional, intent(in) :: ROOTD
      real(8) :: StaTimHS
      real(8) :: RefTimHS
      real(8) :: TimeLoc
      integer :: argcount ! number of command line arguments
      character(2048) :: cmdlinearg ! a single command line argument
      integer :: i ! command line argument counter

      LOG_SCOPE_TRACED("ADCIRC_Init", ADCIRC_TRACING)

#ifndef CMPI
      ! COMM is only used in parallel builds; serial builds keep the argument
      ! so that ADCIRC_Init has the same interface in both
      if (present(COMM)) continue
#endif

      ! jgf51.52.07: Added command line option to just write out
      ! the ADCIRC version number and exit. This is to allow
      ! external applications to query the version number so that
      ! they can behave in version-appropriate ways.
      argcount = command_argument_count() ! get number of command line arguments
      if ((argcount > 0) .and. (myproc == 0)) then
         i = 0
         do while (i < argcount)
            i = i + 1
            call get_command_argument(i, cmdlinearg)
            select case (trim(cmdlinearg))
            case ("-h", "--help")
               call PRINT_HELP()
               call exit(ADCIRC_EXIT_SUCCESS)
            case ("-V", "-v", "--version")
               ! write the version number string ... actually it is
               ! a series of two or three numbers with two digits
               ! separated by dots, e.g., 51.02 or 51.52.07
               write (*, '(a)') trim(adc_version)
               call exit(ADCIRC_EXIT_SUCCESS)
            case ("--hash")
               write (*, '(a)') trim(adc_hash)
               call exit(ADCIRC_EXIT_SUCCESS)
            case default
               ! do nothing, this is some other command line argument
            end select
         end do
      end if

      nabout = ECHO
#ifdef EBUG
      nabout = DEBUG
#endif
      if (present(ROOTD)) then
         ROOTDIR = trim(ROOTD)
      else
         ROOTDIR = '.'
      end if

#ifdef CMPI
      call GET_NUMWRITERS() ! Get number of writer processors from command line arguments
      if (present(COMM)) then
         call MSG_INIT(COMM) ! Init MPI and get MPI-rank of this cpu
      else
         call MSG_INIT() ! Init MPI and get MPI-rank of this cpu
      end if
! Dedicated hotstart writer processors (-H) are no longer supported
      if (hotstartWritersRequested) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="The -H option (dedicated hotstart writer " &
                        //"processors) is no longer supported. Run without -H.")
      end if
      call MAKE_DIRNAME() ! Establish Working Directory Name
      call openLogFile(localdir) ! jgf50.65: open subdomain log files (fort.16)
      call WRITER_INIT() ! Initialize WRITER module
      if (mnwproc > 0 .and. myProc >= mnproc) then
         return
      end if
      call READ_INPUT() ! Establish sizes by reading fort.14 and fort.15
      call MSG_TABLE(NSTA3DD, NSTA3DV, NSTA3DT, & ! Read Message-Passing Tables
                     NSTA3DD_G, NSTA3DV_G, NSTA3DT_G, &
                     IMAP_STA3DD_LG, IMAP_STA3DV_LG, IMAP_STA3DT_LG)
      call MSG_START() ! Startup message passing
#else
      MNPROC = 1 ! Init number of compute procs for serial
      MNWPROC = 0 ! Init number of writer procs for serial
      MNALLPROC = 1 ! Init number of all procs for serial
      MYPROC = 0 ! Init task id for serial
      call MAKE_DIRNAME() ! Establish Working Directory Name
      call openLogFile() ! jgf50.65: open fort.16 file
      call READ_INPUT() ! Establish sizes by reading fort.14 and fort.15
#endif

      ! File Fmt Version  must match Major and Minor numbers.
      FileFmtVersion = &
         VERSION_NUMBER(FileFmtMajor, FileFmtMinor, FileFmtRev)
      if ((NSCREEN /= 0) .and. (MYPROC == 0)) then
         write (ScreenUnit, '(a)') "ADCIRC Version is "//ADC_VERSION
         write (ScreenUnit, '(a)') "ADCIRC Commit hash is "//ADC_HASH
      end if

      !jgf49.44: Allocate memory for reading and writing full domain arrays
      if ((MNPROC > 1) .and. &
          ((myProc == 0) &
           .or. (READ_LOCAL_HOT_START_FILES .eqv. .false.))) then
         call allocateFullDomainIOArrays()
      end if
      if ((MNPROC > 1) .and. (myProc == 0) .and. &
          (WRITE_LOCAL_HARM_FILES .eqv. .false.)) then
         call allocateFullDomainHAIOArrays()
      end if

      ! Initialize variables for bottom friction evaluation considering wet perimeter 11/06/2023
      if (activateVEW1DChannelWetPerimeter .eqv. .true.) then
         call initVEW1DChannelWetPerimeter()
      end if

      !jgf49.44: Initialize harmonic analysis parameters.
      call initHarmonicParameters()

      !     jgf: Initialize default values in wind module and allocate
      !     memory for internal use.
      call initWindModule()

      !jgf53.dev: allocate memory for continuity and momentum eqns
      call initGWCE()
      call init2DMomentumEq()

      !...Compute the reciprocal of the number of nodes in the entire domain

#ifdef CMPI
      call ALLNODES(NP_GLOBAL)
      RNP_GLOBAL = 1.0d0/dble(NP_GLOBAL)
#endif

      !     jgf45.08 put this here so it can be used in hot start
      if (NOLIFA == 0) then
         IFNLFA = 0
      else
         IFNLFA = 1
      end if

      !...Read the self attraction and loading file (fort.24)
      call readSelfAttractionAndLoadingFile()

      !...************* Make the L_N Tidal potential multiple function *************
      !... WJP: Feb 21 2018; moved here before hstart Feb 27 2019
      !..  Corrected long term species (unfortunately there is a mistake in old ADCIRC
      !..  theory report that also propagated into some subsequent publications)
      if (CTIP) then

         do I = 1, NP
            !           LONG-TERM SPECIES
            L_N(0, I) = 1.5d0*cos(SFEA(I))*cos(SFEA(I)) - 1d0
            !           DIURNAL SPECIES
            L_N(1, I) = sin(2d0*SFEA(I))
            !           SEMIDIURNAL SPECIES
            L_N(2, I) = cos(SFEA(I))*cos(SFEA(I))
         end do

      end if

      !...
      !...******************** PROGRAM SETUP SECTION ************************
      !...

      select case (IHOT)
      case (0)
         call COLDSTART()
         ITHS = 0 !TCM20091029 v49.02 Initialized to zero for coldstarts
         ! kmd : added in the initial condition read
      case (17, 67, 68, 367, 368, 567, 568) ! non-portable binary or netcdf hotstart file
         call HOTSTART()
      case DEFAULT
         ! the IHOT value was vetted already, so this code should
         ! not be reachable
         if ((NSCREEN /= 0) .and. (myProc == 0)) then
            write (screenUnit, *) "ERROR: IHOT=", IHOT, " is not valid."
         end if
      end select

      !...Default beginning and ending time steps

      ITIME_BGN = ITHS + 1
      ITIME_END = NT !NT is set in read_input.F

      !...  Determine the number of active elements (MJU), the total number of
      !...  elements (NODELE) and the total area of elements (TotalArea)
      !...  attached to each node
      call totalAreaCalc()
      !
      !
      ! jgf52.08.20: Compute master time which is referenced to the beginning
      ! time of the model run. This is needed here only because we need to
      ! determine whether to load min/max files in subroutine initOutput2D().
      ! If the time in the min/max file is later than the current time
      ! in ADCIRC, we don't load the min/max files. TODO: This is computed
      ! again in subroutine timestep().
      TimeLoc = dble(ITIME_BGN)*DTDP + StaTim*86400.d0
      if (CHOTHS .eqv. .true.) then
         StaTimHS = (dble(ITIME_BGN - 1)*DTDPHS)/86400.d0
         RefTimHS = (dble(ITIME_BGN - 1)*DTDP)/86400.d0
         TimeLoc = dble(ITIME_BGN)*DTDP + (StaTimHS - RefTimHS)*86400.d0
      end if

      !jgf: Initialize water level offset if needed
      allocate (dynamicWaterLevelCorrectionDelta(np)) !time difference in offset tgaContOffset
      dynamicWaterLevelCorrectionDelta = 0.d0
      if (usingDynamicWaterLevelCorrection .eqv. .true.) then
         allocate (dynamicWaterLevelCorrection1(np))
         allocate (dynamicWaterLevelCorrection2(np)) ! used in timestepping
         allocate (dynamicWaterLevelCorrectionSta(mnstae))
         call initDynamicWaterLevelCorrection(timeloc) ! initializes the i/o
         call getDynamicWaterLevelCorrections(dynamicWaterLevelCorrection2, timeloc) ! populate
      end if

      !     jgf51.21.24: Initialize output.
      call initOutput2D(timeloc)

      !...
      !...************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING ***********
      !...

      !...Set flags for nonlinear terms

      if (NOLICA == 0) then
         IFNLCT = 0
      else
         IFNLCT = 1
      end if
      if (NOLICA == 2) then
         IFNLTG = 1
      else
         IFNLTG = 0
      end if
      if (NOLICAT == 0) then
         IFNLCAT = 0
      else
         IFNLCAT = 1
      end if

      ! Corbitt 1203022: Local Advection Flags
      IFNLCTE = IFNLCT
      IFNLCATE = IFNLCAT

      IFWIND = 1
      if (IM == 1) IFWIND = 0

      GORho = G/RhoWat0
      GORhoOAMB = GORho/AMB

      write (16, 1112)
      write (16, 17931)
      if (NSCREEN /= 0 .and. MYPROC == 0) write (ScreenUnit, 1112)
      if (NSCREEN /= 0 .and. MYPROC == 0) write (ScreenUnit, 17931)
1112  format(/, 1x, 79('_'))
17931 format(//, 1x, 'LIMITED RUNTIME INFORMATION SECTION ', //)

      NCChange = 1 !jgf45.06 set flag to set up GWCE LHS on 1st iteration

      return
   end subroutine ADCIRC_Init

   !******************************************************************************
   !*    ADCIRC run routine
   !******************************************************************************
   subroutine ADCIRC_Run(NTIME_STP)
      use GLOBAL, only: iTime, earlyterminate
      use mod_logging, only: allMessage, &
                             logMessage, INFO, t_log_scope, init_log_scope
#ifdef CMPI
      use MOD_SIZES, only: mnproc, mnwproc, myproc
      use WRITER, only: writer_main, writer_pause, sendLabelInfoToWriters
#endif

#ifdef CSWAN
! Casey 090302: We need the following information to couple to unstructured SWAN.
      use Couple2Swan, only: CouplingInterval, &
                             PADCSWAN_RUN
#endif
      use WRITE_OUTPUT, only: writeOutput2D
      use mod_timestep, only: TIMESTEP
#ifdef CMPI
      use MESSENGER, only: writers_active
#endif
      implicit none

      integer, optional, intent(in) :: NTIME_STP
      real(8) :: TimeLoc

      LOG_SCOPE_TRACED("ADCIRC_Run", ADCIRC_TRACING)

      !TCM v49.52.01 added this ifdef (brought in from adcirc_init()
#ifdef CMPI
      if (mnwproc > 0 .and. myproc >= mnproc) then
         call logMessage(INFO, "Enter writer main.")
         call writer_main()
         return
      elseif (mnwproc > 0 .and. myproc == 0) then
         call sendLabelInfoToWriters()
         writers_active = .true.
      end if
#endif

      if (present(NTIME_STP)) then
         ITIME_END = ITIME_BGN + NTIME_STP - 1
      end if

      call allMessage(INFO, "Begin timestepping.")
      do ITIME = ITIME_BGN, ITIME_END
         call TIMESTEP(ITIME, ITIME_BGN, TimeLoc)

#ifdef CSWAN
! Casey 090302: If it is time, then call the following subroutine
!             to then call the SWAN time-stepping subroutine.
         if (mod(ITIME, CouplingInterval) == 0) then
            call PADCSWAN_RUN(ITIME)
         end if
#endif
         if (earlyterminate /= 0) then
            call writeOutput2D(ITIME_END, TimeLoc)
            exit
         else
            !...WRITE OUTPUT
            !   zc - moved this here so that output occurs after
            !        the SWAN time step has been completed.
            call writeOutput2D(ITIME, TimeLoc) ! =>see write_output.F
         end if
      end do

      if (present(NTIME_STP)) then
         ITIME_BGN = ITIME_END + 1
      end if

      !TCM v49.52.01 added the following ifdef to pause writer processors
      ! before exiting adcirc_run phase.  This is to allow them to be
      ! used again if adcirc_run is called more than once.
#ifdef CMPI
      if (MYPROC == 0 .and. MNWPROC > 0) call WRITER_PAUSE()
#endif

      return
   end subroutine ADCIRC_Run

   !******************************************************************************
   !*    ADCIRC finalize routine
   !******************************************************************************
   subroutine ADCIRC_Final(NO_MPI_FINALIZE)
      use GLOBAL, only: iTime
      use mod_logging, only: allMessage, &
                             t_log_scope, init_log_scope
      use MOD_HARMONIC_ANALYSIS, only: solveHarmonicAnalysis
      use WRITE_OUTPUT, only: writeHarmonicAnalysisOutput
#ifdef CMPI
      use MESSENGER, only: msg_fini
#endif
      use MOD_SIZES, only: mnproc, myproc
      use mod_terminate, only: terminate, ADCIRC_EXIT_SUCCESS

      implicit none

      logical, optional, intent(in) :: NO_MPI_FINALIZE

      LOG_SCOPE_TRACED("ADCIRC_Final", ADCIRC_TRACING)
      !...
      !...*************** SOLVE THE HARMONIC ANALYSIS PROBLEM ****************
      !...
      !...  Not on the writer ranks (MYPROC >= MNPROC): they leave ADCIRC_Init
      !...  before the input is read, so the harmonic analysis parameters are
      !...  not set there
      if (MYPROC < MNPROC) then
         call solveHarmonicAnalysis(ITIME)
         call writeHarmonicAnalysisOutput(ITIME)
      end if

      if (present(NO_MPI_FINALIZE)) then
         call terminate(exit_code=ADCIRC_EXIT_SUCCESS, &
                        do_finalize_mpi=.not. NO_MPI_FINALIZE, &
                        message="ADCIRC terminating normally")
      else
         call terminate(exit_code=ADCIRC_EXIT_SUCCESS, &
                        do_finalize_mpi=.true., &
                        message="ADCIRC terminating normally")
      end if

      return
      !-----------------------------------------------------------------------
   end subroutine ADCIRC_Final
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     SUBROUTINE PRINT_HELP
   !-----------------------------------------------------------------------
   !     Print detailed help information for adcirc/padcirc command-line
   !     options
   !-----------------------------------------------------------------------
   subroutine PRINT_HELP()
      implicit none

#ifdef CSWAN
#ifdef CMPI
      write (*, '(A)') 'Usage: padcswan [OPTION]...'
#else
      write (*, '(A)') 'Usage: adcswan [OPTION]...'
#endif
      write (*, '(A)') 'ADvanced CIRCulation model coupled with SWAN for shelves, coastal seas, and estuaries.'
#else
#ifdef CMPI
      write (*, '(A)') 'Usage: padcirc [OPTION]...'
#else
      write (*, '(A)') 'Usage: adcirc [OPTION]...'
#endif
      write (*, '(A)') 'ADvanced CIRCulation model for shelves, coastal seas, and estuaries.'
#endif

      write (*, '(A)') ''
      write (*, '(A)') 'General options:'
      write (*, '(A)') '  -h, --help             display this help and exit'
      write (*, '(A)') '  -v, --version          output version information and exit'
      write (*, '(A)') '      --hash             output git commit hash and exit'
      write (*, '(A)') ''
#ifdef CMPI
      write (*, '(A)') 'Parallel execution options:'
      write (*, '(A)') '  -W NUM                 number of dedicated writer processors for output files'
      write (*, '(A)') '  -Ws NUM                number of dedicated writer processors with split file output'
      write (*, '(A)') '  -H NUM                 number of dedicated writer processors for hotstart files'
      write (*, '(A)') ''
#endif
      write (*, '(A)') 'Report bugs to <https://github.com/adcirc/adcirc/issues>.'
      write (*, '(A)') 'ADCIRC home page: <https://adcirc.org>.'
      write (*, '(A)') 'User documentation: <https://adcirc.github.io/adcirc/>.'

      return
   end subroutine PRINT_HELP
   !-----------------------------------------------------------------------

   !******************************************************************************
   !******************************************************************************

end module ADCIRC_Mod
!******************************************************************************
!******************************************************************************
