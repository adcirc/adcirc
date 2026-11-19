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
#include "logging_macros.h"
!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!    last changes in this file prior to VERSION 44.15                         *
!                                                                             *
!                                                                             *
!   added local LOGICAL CHARMV declaration 04/21/2004                         *
!******************************************************************************
!                                                                      *
!   PADCIRC MODULE  ( HARM )                                           *
!                                                                      *
!   HA_SUBS.FOR     V3.01        11/9/95                               *
!                                                                      *
!   Least Square harmonic analysis of timeseries from ADCIRC2DDI_v27   *
!                                                                      *
!    Notes:                                                            *
!    1.)  Both the left hand side matrix and the right hand side       *
!         forcing vectors are continuously updated in time.  This      *
!         eliminates the need to store time series outputs for later   *
!         harmonic analysis.                                           *
!    2.)  The left hand side matrix and the right hand side forcing    *
!         vectors are output in the hotstart file and can be used to   *
!         perform harmonic analysis on an incomplete run.              *
!    3.)  Frequencies should be in rad/sec,times should be in sec.     *
!                                                                      *
!***********************************************************************
!                                                                      *
!    Program Written by:                                               *
!          R.A. Luettich, IMS UNC                                      *
!          J.J. Westerink, CE ND                                       *
!                                                                      *
!    Program Development History:                                      *
!    1.) lsq_stations_v004 by JJW                                      *
!    2.) LSQEX by RL used in 2D binary extr program                    *
!    3.) LSQRL by RL used in 1D test codes                             *
!    4.) LSQ2D v1.00-v2.26 by RL real time Harmonic Analysis for ADCIRC*
!    5.) HA_SUBS v3.01 by RL real time HA for ADCIRC separate          *
!        subroutines for elevation station, velocity station,          *
!        global elevation and global velocity harmonic analysis        *
!                                                                      *
!***********************************************************************
!                                                                      *
! SUBROUTINE LSQUPDLHS updates the LHS matrix                          *
! SUBROUTINE LSQUPDES updates the RHS load vector for elev stations    *
! SUBROUTINE LSQUPDVS updates the RHS load vector for velocity stations*
! SUBROUTINE LSQUPDEG updates the RHS load vector for elevation global *
! SUBROUTINE LSQUPDVG updates the RHS load vector for velocity global  *
! SUBROUTINE FULSOL fills out, decomposes and solves the matricies     *
! SUBROUTINE LSQSOLES solves & writes output for elevation stations    *
! SUBROUTINE LSQSOLVS solves & writes output for velocity stations     *
! SUBROUTINE LSQSOLEG solves & writes output for elevation global      *
! SUBROUTINE LSQSOLVG solves & writes output for velocity global       *
! SUBROUTINE HAHOUT writes HA parameters & LHS matrix to hotstart file *
! SUBROUTINE HAHOUTES writes elev sta RHS load vector to hotstart file *
! SUBROUTINE HAHOUTVS writes vel sta RHS load vector to hotstart file  *
! SUBROUTINE HAHOUTEG writes glob elev RHS load vector to hotstart file*
! SUBROUTINE HAHOUTVG writes glob vel RHS load vector to hotstart file *
! SUBROUTINE HACOLDS initializes HA param & LHS matrix for cold start  *
! SUBROUTINE HACOLDSES initializes elev sta RHS load vec for cold start*
! SUBROUTINE HACOLDSVS initializes vel sta RHS load vec for cold start *
! SUBROUTINE HACOLDSEG initializes glob ele RHS load vec for cold start*
! SUBROUTINE HACOLDSVG initializes glob vel RHS load vec for cold start*
! SUBROUTINE HAHOTS initializes HA params & LHS matrix for a hot start *
! SUBROUTINE HAHOTSES initializes elev sta RHS load vec for hot start  *
! SUBROUTINE HAHOTSVS initializes vel sta RHS load vec for hot start   *
! SUBROUTINE HAHOTSEG initializes glob elev RHS load vec for hot start *
! SUBROUTINE HAHOTSVG initializes glob vel RHS load vec for hot start  *
!                                                                      *
!***********************************************************************
!                                                                      *
!    INPUT FILES:                                                      *
!      - Frequency information is read in by ADCIRC from unit 15.      *
!                                                                      *
!      - If the model is hot start, input is read from UNIT 67 or 68   *
!                                                                      *
!    OUTPUT FILES:                                                     *
!      UNIT 51 : HARMONIC CONSTITUENT ELEVATION VALUES AT SPECIFIED    *
!                  ELEVATION RECORDING STATION COORDINATES (ASCII)     *
!      UNIT 52 : HARMONIC CONSTITUENT VELOCITY VALUES AT SPECIFIED     *
!                  VELOCITY RECORDING STATION COORDINATES  (ASCII)     *
!      UNIT 53 : HARMONIC CONSTITUENT ELEVATIONS AT ALL NODES (ASCII)  *
!      UNIT 54 : HARMONIC CONSTITUENT VELOCITIES AT ALL NODES (ASCII)  *
!      UNIT 55 : COMPARISON BETWEEN THE MEAN AND VARIANCE OF THE TIME  *
!                  SERIES GENERATED BY THE MODEL AND THE MEAN AND      *
!                  VARIANCE OF A TIME SERIES RESYNTHESIZED FROM THE    *
!                  COMPUTED HARMONIC CONSTITUENTS.  THIS GIVES AN      *
!                  INDICATION OF HOW COMPLETE THE HARMONIC ANALYSIS    *
!                  WAS. (ASCII)                                        *
!      UNIT 67 or 68 : HOT START FILES (BINARY)                        *
!                                                                      *
!***********************************************************************

module mod_harmonic_analysis

   use MOD_SIZES, only: MNP, MNHARF, MNSTAE, MNSTAV, LOCALDIR, &
                        WRITE_LOCAL_HARM_FILES
   use mod_logging, only: screenMessage, logMessage, &
                          t_log_scope, init_log_scope, &
                          allMessage, &
                          DEBUG, ECHO, INFO, WARNING, ERROR
   use GLOBAL, only: IFSPROTS
   use ADC_CONSTANTS, only: pi
   implicit none

   private

   public :: ALLOC_HA, ALLOC_MAIN14, allocateFullDomainHAIOArrays, binaryRead3D, CHARMV, &
             checkHAHotstart, checkHarmonicParameters, EAV, EAVDIF, ELAV, ELAV_g, ELVA, &
             ELVA_g, EMAG, EMAGT, ESQ, EVADIF, FMV, FNAM8, FNAME, GLOELV, GLOELV_g, &
             GLOULV, GLOULV_g, GLOVLV, GLOVLV_g, HA, HAFACE, HAFF, HAFREQ, ICALL, ICHA, &
             IFACE, IFF, IFREQ, IHABEG, IHARIND, IMM, INAMEFR, INF, INFREQ, INHAGE, &
             INHAGV, INHASE, INHASV, initHarmonicParameters, INP, INSTAE, INSTAV, INZ, &
             ITHAF, ITHAS, ITMV, ITUD, MM, NAMEFR, NF, NFREQ, NHAGE, NHAGV, NHAINC, &
             NHASE, NHASV, NTSTEPS, NZ, PHASEDE, PHASEDEN, PHASEDU, PHASEDUT, PHASEDV, &
             PHASEDVT, readBinaryHAHotstart, solveHarmonicAnalysis, STAELV, STAELV_g, &
             STAULV, STAULV_g, STAVLV, STAVLV_g, THAF, THAS, TIMEUD, UAV, UAVDIF, UMAG, &
             UMAGT, updateHarmonicAnalysis, USQ, UVADIF, VAV, VAVDIF, VMAG, VMAGT, VSQ, &
             VVADIF, XELAV, XELVA, XVELAV, XVELAV_g, XVELVA, XVELVA_g, YELAV, YELVA, &
             YVELAV, YVELAV_g, YVELVA, YVELVA_g

   ! jgf49.44: parameters describing the harmonic constituents to be
   ! included in the harmonic analysis of model results
   integer :: IHARIND ! =1 if harmonic analysis was requested, 0 otherwise
   logical :: CHARMV = .false. ! .true. for timeseries reconstruction
   integer :: NFREQ ! number of frequencies
   character(LEN=10), allocatable :: NAMEFR(:) ! "M2", "S2", etc
   real(8), allocatable :: HAFREQ(:) ! frequency (rad/s)
   real(8), allocatable :: HAFF(:) ! nodal factor
   real(8), allocatable :: HAFACE(:) ! equilibrium argument (deg)
   !
   integer :: NZ ! nz=0 if there is a steady component, nz=1 otherwise
   integer :: NF ! nf=0 if no steady constituent, nf=1 otherwise
   !
   integer :: ITHAS ! time step upon which harmonic analysis starts
   integer :: ITHAF ! time step upon which harmonic analysis finishes
   integer :: ITMV ! time step upon which means and variance calcs start
   real(8) :: TIMEBEG ! model time (sec) at which means and var calcs start
   integer :: IHABEG ! 1 h.a. time increment past the start of h.a.
   integer :: ICHA ! spool counter for HA updates
   integer :: NTSTEPS ! number of time steps into means and variance calcs
   integer :: ICALL ! number of times HA update sub has been called
   integer :: ITUD ! time step of most recent HA update
   real(8) :: TIMEUD ! model time of most recent HA update
   !
   real(8), allocatable :: HA(:, :) ! LHS matrix
   integer :: MM ! 2nd dim of HA matrix; 2x num freqs (+ steady)
   real(8), allocatable :: HAP(:) ! used in matrix solution
   real(8), allocatable :: HAX(:) ! used in matrix solution
   !
   !     Load Vectors
   real(8), allocatable, target :: GLOELV(:, :) ! fulldomain elevation
   real(8), allocatable, target :: GLOULV(:, :) ! fulldomain u velocity
   real(8), allocatable, target :: GLOVLV(:, :) ! fulldomain v velocity
   real(8), allocatable, target :: STAELV(:, :) ! station elevation
   real(8), allocatable, target :: STAULV(:, :) ! station u velocity
   real(8), allocatable, target :: STAVLV(:, :) ! station v velocity
   !
   !     Arrays for time series reconstruction.
   real(8), allocatable, target :: XVELAV(:)
   real(8), allocatable, target :: YVELAV(:)
   real(8), allocatable, target :: XVELVA(:)
   real(8), allocatable, target :: YVELVA(:)
   real(8), allocatable, target :: ELAV(:)
   real(8), allocatable, target :: ELVA(:)
   !
   !     jgf49.44: Parameters that control the spatial locations where
   !     harmonic analysis is performed:
   integer :: NHASE ! /= 0 to perform HA at elevation recording stations
   integer :: NHASV ! /= 0 to perform HA at velocity recording stations
   integer :: NHAGE ! /= 0 to do elevation HA at all nodes in the mesh
   integer :: NHAGV ! /= 0 to do velocity HA at all nodes in the mesh
   !
   !     jgf49.44: Parameters that control the calculation of harmonic
   !     constituents:
   real(8) :: THAS ! number of days b/f harmonic analysis starts
   real(8) :: THAF ! number of days at which harmonic analysis ends
   integer :: NHAINC ! time step increment for including data in HA
   real(8) :: FMV ! fraction of HA period for means and var. calcs (0to1)
   !
   !     jgf49.44: Harmonic analysis hotstarting:
   real(8), allocatable, target :: GLOELV_g(:, :)
   real(8), allocatable, target :: STAELV_g(:, :)
   real(8), allocatable, target :: GLOULV_g(:, :)
   real(8), allocatable, target :: GLOVLV_g(:, :)
   real(8), allocatable, target :: STAULV_g(:, :)
   real(8), allocatable, target :: STAVLV_g(:, :)
   real(8), allocatable, target :: ELAV_g(:)
   real(8), allocatable, target :: ELVA_g(:)
   real(8), allocatable, target :: XVELAV_g(:)
   real(8), allocatable, target :: YVELAV_g(:)
   real(8), allocatable, target :: XVELVA_g(:)
   real(8), allocatable, target :: YVELVA_g(:)
   !
   !     jgf49.44: Equivalence variables for reading in hotstart data, either
   !     in serial or in parallel.
   real(8), pointer :: STAELV_pt(:, :), STAULV_pt(:, :), STAVLV_pt(:, :)
   real(8), pointer :: GLOELV_pt(:, :), GLOULV_pt(:, :), GLOVLV_pt(:, :)
   real(8), pointer :: ELAV_pt(:), ELVA_pt(:)
   real(8), pointer :: XVELAV_pt(:), YVELAV_pt(:)
   real(8), pointer :: XVELVA_pt(:), YVELVA_pt(:)
   !
   !     jgf49.44: The following variables are read from the hotstart file and
   !     later compared with the values from the fort.15 to ensure that they
   !     match.
   integer :: INHASE, INHASV, INHAGE, INHAGV, IFLAG
   integer :: INFREQ, INSTAE, INSTAV, INP, INZ, INF, IMM, IICALL
   real(8), allocatable :: IFREQ(:), IFF(:), IFACE(:)
   character(LEN=10), allocatable :: INAMEFR(:)
   !
   character(LEN=16) :: FNAME
   character(LEN=8) :: FNAM8(2)
   !
   !     jgf49.44: Raised output variables to module level and added
   !     full domain counterparts for globalio. These arrays represent
   !     magnitudes and phases for each frequency at each station or node.
   !     Stations:
   real(8), allocatable, target :: EMAG(:, :) ! elevation magnitudes
   real(8), allocatable, target :: PHASEDE(:, :) ! elevation phases
   real(8), allocatable, target :: UMAG(:, :) ! u velocity magnitudes
   real(8), allocatable, target :: VMAG(:, :) ! v velocity magnitudes
   real(8), allocatable, target :: PHASEDU(:, :) ! u velocity phases
   real(8), allocatable, target :: PHASEDV(:, :) ! v velocity phases
   !     Nodes:
   real(8), allocatable, target :: EMAGT(:, :) ! elevation magnitudes
   real(8), allocatable, target :: PHASEDEN(:, :) ! elevation phases
   real(8), allocatable, target :: UMAGT(:, :) ! u velocity magnitudes
   real(8), allocatable, target :: VMAGT(:, :) ! v velocity magnitudes
   real(8), allocatable, target :: PHASEDUT(:, :) ! u velocity phases
   real(8), allocatable, target :: PHASEDVT(:, :) ! v velocity phases
   !     Time series reconstruction:
   real(8), allocatable, target :: EAV(:)
   real(8), allocatable, target :: ESQ(:)
   real(8), allocatable, target :: EAVDIF(:)
   real(8), allocatable, target :: EVADIF(:)
   !
   real(8), allocatable, target :: UAV(:)
   real(8), allocatable, target :: USQ(:)
   real(8), allocatable, target :: UAVDIF(:)
   real(8), allocatable, target :: UVADIF(:)
   !
   real(8), allocatable, target :: VAV(:)
   real(8), allocatable, target :: VSQ(:)
   real(8), allocatable, target :: VAVDIF(:)
   real(8), allocatable, target :: VVADIF(:)

   !jgf51.52.01: Variables used in adcpost.
   real(8), allocatable, target :: XELAV(:)
   real(8), allocatable, target :: YELAV(:)
   real(8), allocatable, target :: XELVA(:)
   real(8), allocatable, target :: YELVA(:)

   !-----------------END OF DECLARATIONS---------------------------------------

contains
   !
   !***********************************************************************
   !  Allocate arays used by LSQ_HARM.
   !
   !  vjp 8/99
   !***********************************************************************
   subroutine ALLOC_HA()
      implicit none
      !
      LOG_SCOPE_TRACED("ALLOC_HA", HARM_TRACING)
      allocate (HAFREQ(MNHARF), HAFF(MNHARF), HAFACE(MNHARF))
      allocate (NAMEFR(MNHARF))
      !
      allocate (HA(2*MNHARF, 2*MNHARF))
      allocate (HAP(2*MNHARF), HAX(2*MNHARF))
      allocate (GLOELV(2*MNHARF, MNP))
      allocate (GLOULV(2*MNHARF, MNP), GLOVLV(2*MNHARF, MNP))
      allocate (STAELV(2*MNHARF, MNSTAE))
      allocate (STAULV(2*MNHARF, MNSTAV), STAVLV(2*MNHARF, MNSTAV))
      !
      allocate (IFREQ(MNHARF), IFF(MNHARF), IFACE(MNHARF))
      allocate (INAMEFR(MNHARF))
      !
      allocate (EMAG(MNHARF, MNSTAE), PHASEDE(MNHARF, MNSTAE))
      allocate (UMAG(MNHARF, MNSTAV), PHASEDU(MNHARF, MNSTAV))
      allocate (VMAG(MNHARF, MNSTAV), PHASEDV(MNHARF, MNSTAV))
      allocate (EMAGT(MNHARF, MNP), PHASEDEN(MNHARF, MNP))
      allocate (UMAGT(MNHARF, MNP), PHASEDUT(MNHARF, MNP))
      allocate (VMAGT(MNHARF, MNP), PHASEDVT(MNHARF, MNP))
      !
      return
      !----------------------------------------------------------------------
   end subroutine ALLOC_HA
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   !...
   !...Allocate space for harmonic analysis means and variance
   !...calculations
   !...
   !----------------------------------------------------------------------
   subroutine ALLOC_MAIN14()
      implicit none
      !
      LOG_SCOPE_TRACED("ALLOC_MAIN14", HARM_TRACING)
      !
      allocate (XVELAV(MNP), YVELAV(MNP), XVELVA(MNP), YVELVA(MNP))
      allocate (ELAV(MNP), ELVA(MNP))
      !
      allocate (EAV(MNP), ESQ(MNP), EAVDIF(MNP), EVADIF(MNP))
      allocate (UAV(MNP), USQ(MNP), UAVDIF(MNP), UVADIF(MNP))
      allocate (VAV(MNP), VSQ(MNP), VAVDIF(MNP), VVADIF(MNP))
      !

      return
      !----------------------------------------------------------------------
   end subroutine ALLOC_MAIN14
   !----------------------------------------------------------------------

   !--------------------------------------------------------------------
   !                    S U B R O U T I N E
   !     A L L O C A T E   F U L L   D O M A I N   H A   I O   A R R A Y S
   !--------------------------------------------------------------------
   !     jgf49.44: Allocates memory for the full domain arrays for i/o
   !     purposes when executing in parallel.
   !--------------------------------------------------------------------
   subroutine allocateFullDomainHAIOArrays()
      use GLOBAL, only: NSTAE_G, NSTAV_G, NP_G
      implicit none
      !
      LOG_SCOPE_TRACED("allocateFullDomainHAIOArrays", HARM_TRACING)
      if (IHARIND == 1) then
         allocate (GLOELV_g(2*MNHARF, NP_G))
         allocate (STAELV_g(2*MNHARF, NSTAE_G))
         allocate (GLOULV_g(2*MNHARF, NP_G))
         allocate (GLOVLV_g(2*MNHARF, NP_G))
         allocate (STAULV_g(2*MNHARF, NSTAV_G))
         allocate (STAVLV_g(2*MNHARF, NSTAV_G))
      end if
      if (CHARMV) then
         allocate (ELAV_g(NP_G))
         allocate (ELVA_g(NP_G))
         allocate (XVELAV_g(NP_G))
         allocate (YVELAV_g(NP_G))
         allocate (XVELVA_g(NP_G))
         allocate (YVELVA_g(NP_G))
      end if
      !--------------------------------------------------------------------
   end subroutine allocateFullDomainHAIOArrays
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E   W R I T E   H A
   !--------------------------------------------------------------------
   !     jgf50.97: Writes LHS matrix to a file for use in debugging.
   !--------------------------------------------------------------------
   subroutine writeHA(timesec, it)
      use MOD_SIZES, only: globaldir, myproc
      implicit none
      real(8), intent(in) :: timesec ! time in seconds since cold start
      integer, intent(in) :: it ! current adcirc time step number
      integer :: i, j ! loop counters
      !
      LOG_SCOPE_TRACED("writeHA", HARM_TRACING)
      if (myproc == 0) then
         open (56, FILE=trim(GLOBALDIR)//'/'//'fort.56', &
               STATUS='UNKNOWN', ACCESS='SEQUENTIAL', ACTION='WRITE', &
               POSITION='APPEND')
         write (56, 2120) timesec, it
         do i = 1, 2*MNHARF
            do j = 1, 2*MNHARF
               write (56, 2130) i, j, ha(i, j)
            end do
         end do
         close (56)
      end if

2120  format(2x, 1pe20.10e3, 5x, I10)
2130  format('ha(', I2, ',', I2, ')=', 1pe20.10e3)

      !--------------------------------------------------------------------
   end subroutine writeHA
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  C H E C K   H A R M O N I C   P A R A M E T E R S
   !--------------------------------------------------------------------
   !     jgf50.41: Checks harmonic analysis settings.
   !     WJP 02.20.2018: Allowing for > 1 options for netcdf output
   !--------------------------------------------------------------------
   subroutine checkHarmonicParameters()
      use GLOBAL, only: NFOVER, useNetCDF, useNetCDFOutput
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      !
      LOG_SCOPE_TRACED("checkHarmonicParameters", HARM_TRACING)

      if ((NHASE < 0) .or. (NHASE > 5)) then
         if (NFOVER == 1) then
            call allMessage(WARNING, "Input value of NHASE is not valid.")
            call allMessage(WARNING, "NHASE will be set to zero.")
         else
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Input value of NHASE is not valid')
         end if
      else if (NHASE == 1) then
         call logMessage(INFO, &
                         'STATION ELEVATION HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 51 IN ASCII FORMAT.')
      elseif (NHASE > 1) then
         call logMessage(INFO, &
                         'STATION ELEVATION HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 51 IN NETCDF FORMAT.')
         useNetCDF = .true.
         useNetCDFOutput = .true.
      end if

      if ((NHASV < 0) .or. (NHASV > 5)) then
         if (NFOVER == 1) then
            call allMessage(WARNING, "Input value of NHASV is not valid.")
            call allMessage(WARNING, "NHASV will be set to zero.")
         else
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Input value of NHASV is not valid')
         end if
      else if (NHASV == 1) then
         call logMessage(INFO, &
                         'STATION VELOCITY HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 52 IN ASCII FORMAT.')
      else if (NHASV > 1) then
         call logMessage(INFO, &
                         'STATION VELOCITY HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 52 IN NETCDF FORMAT.')
         useNetCDF = .true.
         useNetCDFOutput = .true.
      end if

      if ((NHAGE < 0) .or. (NHAGE > 5)) then
         if (NFOVER == 1) then
            call allMessage(WARNING, "Input value of NHAGE is not valid.")
            call allMessage(WARNING, "NHAGE will be set to zero.")
         else
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Input value of NHAGE is not valid')
         end if
      else if (NHAGE == 1) then
         call logMessage(INFO, &
                         'FULL DOMAIN ELEVATION HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 53 IN ASCII FORMAT.')
      else if (NHAGE > 1) then
         call logMessage(INFO, &
                         'FULL DOMAIN ELEVATION HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 53 IN NETCDF FORMAT.')
         useNetCDF = .true.
         useNetCDFOutput = .true.
      end if

      if ((NHAGV < 0) .or. (NHAGV > 5)) then
         if (NFOVER == 1) then
            call allMessage(WARNING, "Input value of NHAGV is not valid.")
            call allMessage(WARNING, "NHAGV will be set to zero.")
         else
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Input value of NHAGV is not valid')
         end if
      else if (NHAGV == 1) then
         call logMessage(INFO, &
                         'FULL DOMAIN VELOCITY HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 54 ASCII FORMAT.')
      else if (NHAGV > 1) then
         call logMessage(INFO, &
                         'FULL DOMAIN VELOCITY HARMONIC ANALYSIS' &
                         //' WILL BE WRITTEN TO UNIT 54 NETCDF FORMAT.')
         useNetCDF = .true.
         useNetCDFOutput = .true.
      end if

      !--------------------------------------------------------------------
   end subroutine checkHarmonicParameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  I N I T   H A R M O N I C   P A R A M E T E R S
   !--------------------------------------------------------------------
   !     jgf49.44: Initializes harmonic analysis arrays and parameters.
   !     Based on HACOLDS, written by rl.
   !     WJP 02.20.2018: Allowing for > 1 options for netcdf output
   !--------------------------------------------------------------------
   subroutine initHarmonicParameters()
      use MOD_SIZES, only: READ_LOCAL_HOT_START_FILES, MNPROC, myProc, &
                           WRITE_LOCAL_HARM_FILES, &
                           WRITE_LOCAL_HOT_START_FILES
      use GLOBAL, only: NHSTAR
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      LOG_SCOPE_TRACED("initHarmonicParameters", HARM_TRACING)
      if (IHARIND == 0) then

         return ! EARLY RETURN if harmonic analysis is not part of this run
      end if
      !
      ! jgf51.44: If subdomain harmonic analysis files will be written
      ! in parallel, subdomain hotstart files are also required for
      ! saving the load vectors and rhs of the harmonic analysis problem.
      if ((mnproc > 1) .and. &
          (write_local_harm_files .eqv. .true.)) then
         call logMessage(INFO, 'Subdomain harmonic analysis files ' &
                         //'were specified with -m on the command line.')
         if ((NHSTAR /= 0) .and. &
             (write_local_hot_start_files .eqv. .false.)) then
            call logMessage(INFO, 'Fulldomain hotstart files will be ' &
                            //'created by default.')
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Fulldomain hotstart output files ' &
                           //'cannot be used in combination with subdomain ' &
                           //'harmonic analysis output files. ' &
                           //'The -s argument can be used on the ADCIRC ' &
                           //'command line to specify that hotstart files ' &
                           //'should be written for each subdomain.')
         end if
         if ((NHSTAR /= 0) .and. &
             (write_local_hot_start_files .eqv. .true.)) then
            call logMessage(INFO, 'Subdomain hotstart files will ' &
                            //'contain subdomain harmonic analysis data.')
         end if
      end if

      IHABEG = ITHAS + NHAINC
      ICHA = 0
      icall = 0
      !
      if (abs(hafreq(1)) <= 0.d0) then
         nz = 0
         nf = 1
      else
         nz = 1
         nf = 0
      end if
      !
      nfreq = nfreq - nf
      mm = 2*nfreq + nf
      !
      ha(:, :) = 0.0d0
      !
      if (NHASE /= 0) then
         STAELV(:, :) = 0.d0
      end if
      if (NHASV /= 0) then
         STAULV(:, :) = 0.d0
         STAVLV(:, :) = 0.d0
      end if
      if (NHAGE /= 0) then
         GLOELV(:, :) = 0.d0
      end if
      if (NHAGV /= 0) then
         GLOULV(:, :) = 0.d0
         GLOVLV(:, :) = 0.d0
      end if
      if (CHARMV .eqv. .true.) then
         ELAV(:) = 0.d0
         XVELAV(:) = 0.d0
         YVELAV(:) = 0.d0
         ELVA(:) = 0.d0
         XVELVA(:) = 0.d0
         YVELVA(:) = 0.d0
      end if !  charmv
      !
      !     jgf49.44: Initialize all the full domain i/o arrays to zero, if
      !     they will be used (i.e., if we are in serial or if we might be
      !     reading or writing full domain arrays in parallel).
      if ((WRITE_LOCAL_HARM_FILES .eqv. .false.) .and. &
          (READ_LOCAL_HOT_START_FILES .eqv. .false.) .and. &
          (MNPROC > 1) .and. &
          (myProc == 0)) then
         GLOELV_g(:, :) = 0.d0
         STAELV_g(:, :) = 0.d0
         GLOULV_g(:, :) = 0.d0
         GLOVLV_g(:, :) = 0.d0
         STAULV_g(:, :) = 0.d0
         STAVLV_g(:, :) = 0.d0
         if (CHARMV .eqv. .true.) then
            ELAV_g(:) = 0.d0
            ELVA_g(:) = 0.d0
            XVELAV_g(:) = 0.d0
            YVELAV_g(:) = 0.d0
            XVELVA_g(:) = 0.d0
            YVELVA_g(:) = 0.d0
         end if
      end if
      !--------------------------------------------------------------------
   end subroutine initHarmonicParameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  U P D A T E   H A R M O N I C  A N A L Y S I S
   !--------------------------------------------------------------------
   !     jgf49.44: Removed from timestep.F and placed here.
   !     WJP 02.20.2018: Allowing for > 1 options for netcdf output
   !...
   !...  IF IHARIND=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
   !...  AND ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE
   !...  HARMONIC ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD
   !...  SHOULD BE USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN
   !...  ON 32 BIT WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS
   !...  IS DONE IN DOUBLE PRECISION.
   !...
   !--------------------------------------------------------------------
   subroutine updateHarmonicAnalysis(it, timeh)
      use GLOBAL, only: ET00, UU00, VV00, ETA2, UU2, VV2, &
                        STAIE1, STAIE2, STAIE3, STAIV1, STAIV2, STAIV3, &
                        NSTAE, NSTAV, NNE, NNV
      use MESH, only: NP, NM, UVECTMP, VVECTMP, DRVMAP2DSPVEC
      implicit none
      integer, intent(in) :: it ! current ADCIRC timestep
      real(8), intent(in) :: timeh ! time (sec) incl. STATIM and REFTIM
      !
      integer :: i ! loop counters
      real(8) :: EE1, EE2, EE3 ! nodal elevations around an element
      real(8) :: U11, U22, U33 ! nodal u velocities around an element
      real(8) :: V11, V22, V33 ! nodal v velocities around an element
      !
      LOG_SCOPE_TRACED("updateHarmonicAnalysis", HARM_TRACING)
      !...
      if (IHARIND == 1) then
         if ((IT > ITHAS) .and. (IT <= ITHAF)) then
            if (ICHA == NHAINC) ICHA = 0
            ICHA = ICHA + 1
            if (ICHA == NHAINC) then
               !...
               !.....UPDATE THE LHS MATRIX
               !...
               !.... BEG: DW 2019
               if ((NHASV == 1) .or. (NHAGV == 1)) then
                  if (IFSPROTS == 1) then
                     UVECTMP = UU2; 
                     VVECTMP = VV2; 
                     call DRVMAP2DSPVEC(UU2, VV2, &
                                        UVECTMP, VVECTMP, NP, FWD=.false.)
                  end if
               end if
               !.... DW

               call LSQUPDLHS(timeh, IT)

               !...  IF DESIRED COMPUTE ELEVATION STATION INFORMATION AND UPDATE LOAD
               !.....VECTOR
               !...
               if (NHASE /= 0) then
                  do I = 1, NSTAE
                     EE1 = ETA2(NM(NNE(I), 1))
                     EE2 = ETA2(NM(NNE(I), 2))
                     EE3 = ETA2(NM(NNE(I), 3))
                     ET00(I) = EE1*STAIE1(I) + EE2*STAIE2(I) + EE3*STAIE3(I)
                  end do
                  call LSQUPDES(ET00, NSTAE)
               end if
               !...  IF DESIRED COMPUTE VELOCITY STATION INFORMATION AND UPDATE LOAD
               !.....VECTOR
               !...
               if (NHASV /= 0) then
                  do I = 1, NSTAV
                     U11 = UU2(NM(NNV(I), 1))
                     U22 = UU2(NM(NNV(I), 2))
                     U33 = UU2(NM(NNV(I), 3))
                     V11 = VV2(NM(NNV(I), 1))
                     V22 = VV2(NM(NNV(I), 2))
                     V33 = VV2(NM(NNV(I), 3))
                     UU00(I) = U11*STAIV1(I) + U22*STAIV2(I) + U33*STAIV3(I)
                     VV00(I) = V11*STAIV1(I) + V22*STAIV2(I) + V33*STAIV3(I)
                  end do
                  call LSQUPDVS(UU00, VV00, NSTAV)
               end if
               !...
               !.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
               !...
               if (NHAGE /= 0) call LSQUPDEG(ETA2, NP)
               !...
               !.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
               !...
               if (NHAGV /= 0) call LSQUPDVG(UU2, VV2, NP)

               !.... BEG: DW 2019
               if ((NHASV == 1) .or. (NHAGV == 1)) then
                  if (IFSPROTS == 1) then
                     UU2 = UVECTMP; 
                     VV2 = VVECTMP; 
                  end if
               end if
               !.... DW

            end if
         end if

         !...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            if (IT > ITMV) then
               NTSTEPS = NTSTEPS + 1
               do I = 1, NP
                  ELAV(I) = ELAV(I) + ETA2(I)
                  XVELAV(I) = XVELAV(I) + UU2(I)
                  YVELAV(I) = YVELAV(I) + VV2(I)
                  ELVA(I) = ELVA(I) + ETA2(I)*ETA2(I)
                  XVELVA(I) = XVELVA(I) + UU2(I)*UU2(I)
                  YVELVA(I) = YVELVA(I) + VV2(I)*VV2(I)
               end do
            end if
         end if !   charmv

      end if

      !--------------------------------------------------------------------
   end subroutine updateHarmonicAnalysis
   !--------------------------------------------------------------------

   !***********************************************************************
   !   Subroutine to update the Left Hand Side Matrix                     *
   !                                                                      *
   !  TIMELOC  - ABSOLUTE MODEL TIME (SEC)                                   *
   !  IT    - MODEL TIME STEP                                             *
   !  icall - number of times the subroutine has been called              *
   !  a     - Left Hand Side Matrix                                       *
   !                                                                      *
   !                        RL 11/7/95                                    *
   !***********************************************************************
   !
   subroutine LSQUPDLHS(TIMELOC, IT)
      implicit none
      integer, intent(in) :: IT
      integer :: I, J, I1, I2, J1, J2
      real(8) :: TF1, TF2
      real(8), intent(in) :: TIMELOC
      !
      LOG_SCOPE_TRACED("LSQUPDLHS", HARM_TRACING)
      !
      icall = icall + 1
      !
      !***** Update the Left Hand Side Matrix
      !     Note: this is a symmetric matrix and therefore only store the
      !     upper triangular part.  The lower part will be filled out in
      !     SUBROUTINE FULSOL prior to the matrix's decomposition

      !     Take care of the steady constituent if included in the analysis

      if (nf == 1) then
         ha(1, 1) = dble(icall)
         do j = 1, nfreq
            tf1 = hafreq(j + nf)*TIMELOC
            ha(1, 2*j) = ha(1, 2*j) + cos(tf1)
            ha(1, 2*j + 1) = ha(1, 2*j + 1) + sin(tf1)
         end do
      end if

      !   Take care of the other constituents

      do i = 1, nfreq
         do j = i, nfreq
            i1 = 2*i - (1 - nf)
            i2 = i1 + 1
            j1 = 2*j - (1 - nf)
            j2 = j1 + 1
            tf1 = hafreq(i + nf)*TIMELOC
            tf2 = hafreq(j + nf)*TIMELOC
            ha(i1, j1) = ha(i1, j1) + cos(tf1)*cos(tf2)
            ha(i1, j2) = ha(i1, j2) + cos(tf1)*sin(tf2)
            ha(i2, j2) = ha(i2, j2) + sin(tf1)*sin(tf2)
            if (i2 <= j1) ha(i2, j1) = ha(i2, j1) + sin(tf1)*cos(tf2)
         end do
      end do

      !   Record update time and time step

      TIMEUD = TIMELOC
      ITUD = IT

      return
   end subroutine LSQUPDLHS

   !***********************************************************************
   !   Subroutine to update the Right Hand Side Load Vectors for the      *
   !   elevation station harmonic analysis.                               *
   !                                                                      *
   !  STAE  - STATION ELEVATION VALUES USED TO UPDATE LOAD VECTORS        *
   !  NSTAE - NUMBER OF TIDAL ELEVATION RECORDING STATIONS                *
   !                                                                      *
   !  STAELV - station elevation load vector                              *
   !                                                                      *
   !                        RL 11/8/95                                    *
   !***********************************************************************
   !
   subroutine LSQUPDES(STAE, NSTAE)
      implicit none
      integer, intent(in) :: NSTAE
      integer :: N, I, I1, I2
      real(8) :: TF1, CTF1, STF1
      real(8), intent(in) :: STAE(MNSTAE)
      !
      LOG_SCOPE_TRACED("LSQUPDES", HARM_TRACING)
      !
      !***** Update the Right Hand Side Load Vectors
      !
      !   Take care of the steady constituent if included in the analysis

      if (nz == 0) then
         do n = 1, NSTAE
            STAELV(1, N) = STAELV(1, N) + STAE(N)
         end do
      end if

      !   Take care of the other constituents

      do i = 1, nfreq
         i1 = 2*i - nz
         i2 = i1 + 1
         tf1 = hafreq(i + nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n = 1, NSTAE
            STAELV(I1, N) = STAELV(I1, N) + STAE(N)*CTF1
            STAELV(I2, N) = STAELV(I2, N) + STAE(N)*STF1
         end do
      end do
      !
      return
   end subroutine LSQUPDES

   !***********************************************************************
   !   Subroutine to update the Right Hand Side Load Vectors for the      *
   !   velocity station harmonic analysis.                                *
   !                                                                      *
   !  STAU  - STATION U VELOCITY VALUES USED TO UPDATE LOAD VECTORS       *
   !  STAV  - STATION V VELOCITY VALUES USED TO UPDATE LOAD VECTORS       *
   !  NSTAV - NUMBER OF TIDAL CURRENT RECORDING STATIONS                  *
   !                                                                      *
   !  STAULV - station u velocity load vector                             *
   !  STAVLV - station v velocity load vector                             *
   !                                                                      *
   !                        RL 11/8/95                                    *
   !***********************************************************************
   !
   subroutine LSQUPDVS(STAU, STAV, NSTAV)
      implicit none
      integer, intent(in) :: NSTAV
      integer :: N, I, I1, I2
      real(8) :: TF1, CTF1, STF1
      real(8), intent(in) :: STAU(MNSTAV), STAV(MNSTAV)
      !
      LOG_SCOPE_TRACED("LSQUPDVS", HARM_TRACING)
      !
      !***** Update the Right Hand Side Load Vectors
      !
      !     Take care of the steady constituent if included in the analysis

      if (nz == 0) then
         do n = 1, NSTAV
            STAULV(1, N) = STAULV(1, N) + STAU(N)
            STAVLV(1, N) = STAVLV(1, N) + STAV(N)
         end do
      end if

      !     Take care of the other constituents

      do i = 1, nfreq
         i1 = 2*i - nz
         i2 = i1 + 1
         tf1 = hafreq(i + nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n = 1, NSTAV
            STAULV(I1, N) = STAULV(I1, N) + STAU(N)*CTF1
            STAVLV(I1, N) = STAVLV(I1, N) + STAV(N)*CTF1
            STAULV(I2, N) = STAULV(I2, N) + STAU(N)*STF1
            STAVLV(I2, N) = STAVLV(I2, N) + STAV(N)*STF1
         end do
      end do
      !
      return
   end subroutine LSQUPDVS

   !***********************************************************************
   !   Subroutine to update the Right Hand Side Load Vectors for the      *
   !   global elevation harmonic analysis.                                *
   !                                                                      *
   !  GLOE  - GLOBAL ELEVATION VALUES USED TO UPDATE LOAD VECTORS         *
   !  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
   !                                                                      *
   !  GLOELV - global elevation load vector                               *
   !                                                                      *
   !                        RL 11/8/95                                    *
   !***********************************************************************
   !
   subroutine LSQUPDEG(GLOE, NP)
      implicit none
      integer, intent(in) :: NP
      integer :: I, N, I1, I2
      real(8) :: TF1, CTF1, STF1
      real(8), intent(in) :: GLOE(MNP)
      !
      LOG_SCOPE_TRACED("LSQUPDEG", HARM_TRACING)
      !
      !*****Update the Right Hand Side Load Vectors
      !
      !     Take care of the steady constituent if included in the analysis

      if (nz == 0) then
         do n = 1, np
            GLOELV(1, N) = GLOELV(1, N) + GLOE(N)
         end do
      end if

      !     Take care of the other constituents

      do i = 1, nfreq
         i1 = 2*i - nz
         i2 = i1 + 1
         tf1 = hafreq(i + nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n = 1, np
            GLOELV(I1, N) = GLOELV(I1, N) + GLOE(N)*CTF1
            GLOELV(I2, N) = GLOELV(I2, N) + GLOE(N)*STF1
         end do
      end do
      !

      return
   end subroutine LSQUPDEG

   !***********************************************************************
   !   Subroutine to update the Right Hand Side Load Vectors for the      *
   !   global velocity harmonic analysis.                                 *
   !                                                                      *
   !  GLOU  - GLOBAL U VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
   !  GLOV  - GLOBAL V VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
   !  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
   !                                                                      *
   !  GLOULV - global u velocity load vector                              *
   !  GLOVLV - global v velocity load vector                              *
   !                                                                      *
   !                        RL 11/8/95                                    *
   !***********************************************************************
   !
   subroutine LSQUPDVG(GLOU, GLOV, NP)
      implicit none
      integer, intent(in) :: NP
      integer :: I1, I2, N, I
      real(8) :: TF1, CTF1, STF1
      real(8), intent(in) :: GLOU(MNP), GLOV(MNP)
      !
      LOG_SCOPE_TRACED("LSQUPDVG", HARM_TRACING)
      !
      !*****Update the Right Hand Side Load Vectors
      !
      !     Take care of the steady constituent if included in the analysis

      if (nz == 0) then
         do n = 1, np
            GLOULV(1, N) = GLOULV(1, N) + GLOU(N)
            GLOVLV(1, N) = GLOVLV(1, N) + GLOV(N)
         end do
      end if

      !     Take care of the other constituents

      do i = 1, nfreq
         i1 = 2*i - nz
         i2 = i1 + 1
         tf1 = hafreq(i + nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n = 1, np
            GLOULV(I1, N) = GLOULV(I1, N) + GLOU(N)*CTF1
            GLOVLV(I1, N) = GLOVLV(I1, N) + GLOV(N)*CTF1
            GLOULV(I2, N) = GLOULV(I2, N) + GLOU(N)*STF1
            GLOVLV(I2, N) = GLOVLV(I2, N) + GLOV(N)*STF1
         end do
      end do
      !

      return
   end subroutine LSQUPDVG

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  S O L V E   H A R M O N I C   A N A L Y S I S
   !--------------------------------------------------------------------
   !     jgf49.44: Solves for amplitudes and phases at the nodes and/or
   !     stations as specified. Also finishes time series reconstruction
   !     if necessary.
   !     WJP 02.20.2018: Allowing for > 1 options for netcdf output
   !--------------------------------------------------------------------
   subroutine solveHarmonicAnalysis(ITIME)
      use GLOBAL, only: STATIM, REFTIM, DTDP
      use MESH, only: NP
      implicit none
      integer, intent(in) :: ITIME
      integer :: I
      !
      LOG_SCOPE_TRACED("solveHarmonicAnalysis", HARM_TRACING)
      !
      if ((IHARIND /= 1) .or. (ITIME <= ITHAS)) then
         return ! EARLY RETURN if we weren't supposed to do harmonic analysis
         ! or it hasn't started yet.
      end if
      !
      !...Compute means and variances for checking the harmonic analysis results
      !...Accumulate mean and variance at each node.
      if ((CHARMV .eqv. .true.) .and. (FMV > 1.0d-3)) then
         do I = 1, NP
            ELAV(I) = ELAV(I)/dble(NTSTEPS)
            XVELAV(I) = XVELAV(I)/dble(NTSTEPS)
            YVELAV(I) = YVELAV(I)/dble(NTSTEPS)
            ELVA(I) = ELVA(I)/dble(NTSTEPS) - ELAV(I)*ELAV(I)
            XVELVA(I) = XVELVA(I)/dble(NTSTEPS) - XVELAV(I)*XVELAV(I)
            YVELVA(I) = YVELVA(I)/dble(NTSTEPS) - YVELAV(I)*YVELAV(I)
         end do
         TIMEBEG = dble(ITMV)*DTDP + (STATIM - REFTIM)*86400.d0
      end if

      !......Fill out and decompose the LHS harmonic analaysis matrix

      call FULSOL(0)

      !......Solve the harmonic analysis problems

      if (NHAGE /= 0) call LSQSOLEG()
      if (NHAGV /= 0) call LSQSOLVG()
      if (NHASE /= 0) call LSQSOLES()
      if (NHASV /= 0) call LSQSOLVS()

      !--------------------------------------------------------------------
   end subroutine solveHarmonicAnalysis
   !--------------------------------------------------------------------

   !***********************************************************************
   !   Subroutine to fill out, decompose and solve the lsq system         *
   !   Solves system a*x=b by l*d*l(tr) decomp in full storage mode       *
   !                                                                      *
   !   NOTE: This routine has been modified so that the filling out and   *
   !         decomposition (and only those operations) are done if        *
   !         idecom=0.                                                    *
   !                                                                      *
   !   mm  -  actual dimension of a matrix                                *
   !                                                                      *
   !                        rl 11/7/95                                    *
   !***********************************************************************
   !
   subroutine fulsol(idecom)
      implicit none
      integer, intent(in) :: idecom
      integer :: i, j, ir, ire, k, jr
      real(8), allocatable :: c(:), y(:)
      !
      LOG_SCOPE_TRACED("fulsol", HARM_TRACING)
      !
      !**** If only want to fill out matrix and decompose
      !
      if (idecom == 0) then

         !     Set up the lower triangular part of the LHS a matrix

         do j = 1, mm
            do i = j, mm
               ha(i, j) = ha(j, i)
            end do
         end do

         !     Decomposition of matrix a

         outer: do ir = 1, mm
            ire = ir + 1
            do j = ire, mm
               ha(ir, j) = ha(ir, j)/ha(ir, ir)
            end do
            if (ire > mm) exit outer
            do j = ire, mm
               do k = ire, mm
                  ha(k, j) = ha(k, j) - ha(k, ir)*ha(ir, j)
               end do
            end do
            do j = ire, mm
               ha(j, ir) = 0.0d0
            end do
         end do outer

         return ! EARLY RETURN if idecom.eq.0
      end if

      !...  solve for y by forward substitution for l*y=p

      allocate (c(2*MNHARF), y(2*MNHARF))
      !
      do ir = 1, mm
         y(ir) = hap(ir)
         do jr = 1, ir - 1
            y(ir) = y(ir) - ha(jr, ir)*y(jr)
         end do
      end do

      !...  calculate c=d**(-1)*y

      do ir = 1, mm
         c(ir) = y(ir)/ha(ir, ir)
      end do

      !...  solve for x by back-substituting for l(tr)*x=c

      ir = mm
      do while (.true.)
         hax(ir) = c(ir)
         do jr = ir + 1, mm
            hax(ir) = hax(ir) - ha(ir, jr)*hax(jr)
         end do
         ir = ir - 1
         if (ir >= 1) cycle
         exit
      end do

   end subroutine fulsol

   !***********************************************************************
   !   Subroutine to solve the system and write output for elevation      *
   !   stations.                                                          *
   !                                                                      *
   !   nf=0  if no steady constituent                                     *
   !   nf=1 if steady constituent                                         *
   !                                                                      *
   !  closed unit 51            08/23/05                                  *
   !***********************************************************************
   !
   subroutine LSQSOLES()
      use GLOBAL, only: NSTAE
      implicit none
      integer :: N, I, K, I1, I2
      real(8) :: CONVRD
      real(8) :: PHASEE
      !
      LOG_SCOPE_TRACED("LSQSOLES", HARM_TRACING)
      !
      convrd = 180.d0/pi
      !
      !**** AT each STATION TRANSFER each load vector to p and solve the system
      !
      do N = 1, NSTAE
         do k = 1, mm
            hap(k) = STAELV(k, n)
         end do
         call fulsol(n)
         !
         !        Compute amplitude and phase for each frequency making sure that the
         !        phase is between 0 and 360 deg.
         do i = 1, nfreq + nf
            if ((nf == 1) .and. (i == 1)) then
               emag(i, n) = hax(i)/haff(i)
               phasee = 0.d0
            else
               i1 = 2*i - 1 - nf
               i2 = i1 + 1
               emag(i, n) = sqrt(hax(i1)*hax(i1) + hax(i2)*hax(i2))/haff(i)
               if ((abs(hax(i1)) <= 0.d0) .and. (abs(hax(i2)) <= 0.d0)) then
                  phasee = 0.d0
               else
                  phasee = atan2(hax(i2), hax(i1))
               end if
            end if
            phasede(i, n) = convrd*phasee + haface(i)
            if (phasede(i, n) < 0.d0) phasede(i, n) = phasede(i, n) + 360.d0
            if (phasede(i, n) >= 360.d0) phasede(i, n) = phasede(i, n) - 360.d0
         end do
      end do
      !

      return
   end subroutine lsqsoles

   !***********************************************************************
   !   Subroutine to solve the system and write output for velocity       *
   !   stations.                                                          *
   !                                                                      *
   !   nf=0  if no steady constituent                                     *
   !   nf=1  if steady constituent                                        *
   !                                                                      *
   !                        R.L. 11/8/95                                  *
   !   closed unit 52            08/23/05                                 *
   !***********************************************************************
   !
   subroutine LSQSOLVS()
      use GLOBAL, only: NSTAV
      implicit none
      integer :: I, N, K, I1, I2
      real(8) :: CONVRD
      real(8) :: PHASEU, PHASEV
      real(8), allocatable :: Y(:)
      !
      LOG_SCOPE_TRACED("LSQSOLVS", HARM_TRACING)
      !
      convrd = 180.d0/pi
      !
      allocate (Y(2*MNHARF))
      !
      !**** AT each STATION, transfer each load vector to p, solve system,
      !**** then write results
      !
      do N = 1, NSTAV
         do k = 1, mm
            hap(k) = STAVLV(k, n)
         end do
         call fulsol(n)
         do k = 1, mm
            y(k) = hax(k)
         end do
         do k = 1, mm
            hap(k) = STAULV(k, n)
         end do
         call fulsol(n)

         !        Compute amplitude and phase for each frequency making sure that the
         !        phase is between 0 and 360 deg.

         do i = 1, nfreq + nf
            if ((nf == 1) .and. (i == 1)) then
               umag(i, n) = hax(i)/haff(i)
               vmag(i, n) = y(i)/haff(i)
               phaseu = 0.d0
               phasev = 0.d0
            else
               i1 = 2*i - 1 - nf
               i2 = i1 + 1
               umag(i, n) = sqrt(hax(i1)*hax(i1) + hax(i2)*hax(i2))/haff(i)
               vmag(i, n) = sqrt(y(i1)*y(i1) + y(i2)*y(i2))/haff(i)
               if ((abs(hax(i1)) <= 0.d0) .and. (abs(hax(i2)) <= 0.d0)) then
                  phaseu = 0.d0
               else
                  phaseu = atan2(hax(i2), hax(i1))
               end if
               if ((abs(y(i1)) <= 0.d0) .and. (abs(y(i2)) <= 0.d0)) then
                  phasev = 0.d0
               else
                  phasev = atan2(y(i2), y(i1))
               end if
            end if
            phasedu(i, n) = convrd*phaseu + haface(i)
            if (phasedu(i, n) < 0.d0) phasedu(i, n) = phasedu(i, n) + 360.d0
            if (phasedu(i, n) >= 360.d0) phasedu(i, n) = phasedu(i, n) - 360.d0
            phasedv(i, n) = convrd*phasev + haface(i)
            if (phasedv(i, n) < 0.d0) phasedv(i, n) = phasedv(i, n) + 360.d0
            if (phasedv(i, n) >= 360.d0) phasedv(i, n) = phasedv(i, n) - 360.d0
         end do
      end do

      return
   end subroutine lsqsolvs

   !***********************************************************************
   !   Subroutine to solve the system and write output for elevation      *
   !   globally.                                                          *
   !                                                                      *
   !   nf=0  if no steady constituent                                     *
   !   nf=1  if steady constituent                                        *
   !                                                                      *
   !                        R.L. 11/8/95                                  *
   !                                                                      *
   !   added local LOGICAL CHARMV declaration 04/09/2004                  *
   !   closed unit 53 08/23/05                                            *
   !***********************************************************************
   !
   subroutine LSQSOLEG()
      use GLOBAL, only: DTDP
      implicit none
      integer :: N, K, I, I1, I2, IT, IFR, NEAVMAX, NEAVMIN, &
                 NEVAMAX, NEVAMIN
      real(8)  :: CONVRD
      real(8) :: EAVMAX, EVAMAX, EAVMIN, EVAMIN
      real(8) :: TIMELOC, RSE, FTIME
      real(8), allocatable :: PHASEE(:), EMAG(:)
      !
      LOG_SCOPE_TRACED("LSQSOLEG", HARM_TRACING)
      !
      convrd = 180.d0/pi
      !
      allocate (PHASEE(MNHARF), EMAG(MNHARF))
      !
      if (CHARMV) then
         EAVMAX = -999.d0
         EVAMAX = -999.d0
         EAVMIN = 999.d0
         EVAMIN = 999.d0
      end if
      !
      !***** AT each node transfer each load vector to p, solve and write output
      !
      do N = 1, MNP
         do k = 1, mm
            hap(k) = GLOELV(k, n)
         end do
         call fulsol(n)
         !
         !        Compute amplitude and phase for each frequency making sure that the
         !        phase is between 0 and 360 deg.  Then write output.
         !
         do i = 1, nfreq + nf
            if ((nf == 1) .and. (i == 1)) then
               emag(i) = hax(i)
               emagt(i, n) = emag(i)/haff(i)
               phasee(i) = 0.d0
            else
               i1 = 2*i - 1 - nf
               i2 = i1 + 1
               emag(i) = sqrt(hax(i1)*hax(i1) + hax(i2)*hax(i2))
               emagt(i, n) = emag(i)/haff(i)
               if ((abs(hax(i1)) <= 0.d0) .and. (abs(hax(i2)) <= 0.d0)) then
                  phasee(i) = 0.d0
               else
                  phasee(i) = atan2(hax(i2), hax(i1))
               end if
            end if
            phaseden(i, n) = convrd*phasee(i) + haface(i)
            if (phaseden(i, n) < 0.d0) phaseden(i, n) = phaseden(i, n) + 360.d0
            if (phaseden(i, n) >= 360.d0) phaseden(i, n) = phaseden(i, n) - 360.d0
         end do

         if (CHARMV) then
            eav(n) = 0.d0
            esq(n) = 0.d0
            do it = 1, ntsteps
               TIMELOC = TIMEBEG + DTDP*dble(IT)
               rse = 0.d0
               do ifr = 1, nfreq + nf
                  ftime = hafreq(ifr)*TIMELOC
                  rse = rse + emag(ifr)*cos(ftime - phasee(ifr))
               end do
               eav(n) = eav(n) + rse
               esq(n) = esq(n) + rse*rse
            end do
            !
            eav(n) = eav(n)/dble(ntsteps)
            esq(n) = esq(n)/dble(ntsteps) - eav(n)*eav(n)
            if (abs(elav(n)) <= 0.d0) then
               if (abs(eav(n)) <= 0.d0) eavdif(n) = 1.0d0
               if (abs(eav(n)) > 0.d0) eavdif(n) = 99d19
            else
               eavdif(n) = eav(n)/elav(n)
            end if
            if (abs(elva(n)) <= 0.d0) then
               if (abs(esq(n)) <= 0.d0) evadif(n) = 1.0d0
               if (abs(esq(n)) > 0.d0) evadif(n) = 99d19
            else
               evadif(n) = esq(n)/elva(n)
            end if
            !
            if (EAVDIF(n) > EAVMAX) then
               EAVMAX = EAVDIF(n)
               NEAVMAX = n
            end if
            if (EAVDIF(n) < EAVMIN) then
               EAVMIN = EAVDIF(n)
               NEAVMIN = n
            end if
            if (EVADIF(n) > EVAMAX) then
               EVAMAX = EVADIF(n)
               NEVAMAX = n
            end if
            if (EVADIF(n) < EVAMIN) then
               EVAMIN = EVADIF(n)
               NEVAMIN = n
            end if
         end if ! charmv
      end do
      !
      if (charmv) then
         !
         write (16, 7740)
7740     format(///, 5x, 'THE LARGEST VALUES OF THE RATIO ', &
                 'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
         write (16, 7741) EAVMAX, NEAVMAX
         write (16, 7742) EVAMAX, NEVAMAX
         write (16, 7747)
7747     format(/, 5x, 'THE LOWEST VALUES OF THE RATIO ', &
                 'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
         write (16, 7741) EAVMIN, NEAVMIN
         write (16, 7742) EVAMIN, NEVAMIN
7741     format(9x, '  AVERAGE ELEVATION RATIO = ', E15.7, ' AT NODE ', I8)
7742     format(9x, ' VARIANCE ELEVATION RATIO = ', E15.7, ' AT NODE ', I8)
         !
      end if ! charmv
      !

      return
   end subroutine lsqsoleg

   !***********************************************************************
   !   Subroutine to solve the system and write output for velocity       *
   !   globally.                                                          *
   !                                                                      *
   !   nf=0  if no steady constituent                                     *
   !   nf=1  if steady constituent                                        *
   !                                                                      *
   !                        R.L. 11/10/95                                 *
   !                                                                      *
   !   added local LOGICAL CHARMV declaration 04/09/2004                  *
   !   closed unit 54 08/23/05                                            *
   !***********************************************************************
   !
   subroutine LSQSOLVG()
      use GLOBAL, only: DTDP
      implicit none
      integer :: I, N, K, I1, I2, IT, IFR
      integer :: NUAVMAX, NUAVMIN, NVAVMAX, NVAVMIN, NUVAMAX, NUVAMIN, &
                 NVVAMAX, NVVAMIN
      real(8) :: TIMELOC, FTIME, RSU, RSV
      real(8) :: UAVMAX, VAVMAX, UVAMAX, VVAMAX, UAVMIN, VAVMIN, &
                 UVAMIN, VVAMIN
      real(8) :: CONVRD
      real(8), allocatable :: UMAGG(:), VMAGG(:), PHASEUG(:), PHASEVG(:)
      real(8), allocatable :: Y(:)
      !
      LOG_SCOPE_TRACED("LSQSOLVG", HARM_TRACING)
      !
      convrd = 180.d0/pi
      !
      allocate (UMAGG(MNHARF), VMAGG(MNHARF))
      allocate (PHASEUG(MNHARF), PHASEVG(MNHARF))
      allocate (Y(2*MNHARF))
      !
      !
      if (charmv) then
         UAVMAX = -999.d0
         VAVMAX = -999.d0
         UVAMAX = -999.d0
         VVAMAX = -999.d0
         UAVMIN = 999.d0
         VAVMIN = 999.d0
         UVAMIN = 999.d0
         VVAMIN = 999.d0
      end if ! charmv
      !
      !***** AT each node transfer each load vector to p, solve and write output
      !
      do N = 1, MNP
         do k = 1, mm
            hap(k) = GLOVLV(k, n)
         end do
         call fulsol(n)
         do k = 1, mm
            y(k) = hax(k)
         end do
         do k = 1, mm
            hap(k) = GLOULV(k, n)
         end do
         call fulsol(n)
         do i = 1, nfreq + nf
            if ((nf == 1) .and. (i == 1)) then
               umagg(i) = hax(i)
               umagt(i, n) = umagg(i)/haff(i)
               vmagg(i) = y(i)
               vmagt(i, n) = vmagg(i)/haff(i)
               phaseug(i) = 0.d0
               phasevg(i) = 0.d0
            else
               i1 = 2*i - 1 - nf
               i2 = i1 + 1
               umagg(i) = sqrt(hax(i1)*hax(i1) + hax(i2)*hax(i2))
               umagt(i, n) = umagg(i)/haff(i)
               vmagg(i) = sqrt(y(i1)*y(i1) + y(i2)*y(i2))
               vmagt(i, n) = vmagg(i)/haff(i)
               if ((abs(hax(i1)) <= 0.d0) .and. (abs(hax(i2)) <= 0.d0)) then
                  phaseug(i) = 0.d0
               else
                  phaseug(i) = atan2(hax(i2), hax(i1))
               end if
               if ((abs(y(i1)) <= 0.d0) .and. (abs(y(i2)) <= 0.d0)) then
                  phasevg(i) = 0.d0
               else
                  phasevg(i) = atan2(y(i2), y(i1))
               end if
            end if
            phasedut(i, n) = convrd*phaseug(i) + haface(i)
            if (phasedut(i, n) < 0.d0) phasedut(i, n) = phasedut(i, n) + 360.d0
            if (phasedut(i, n) >= 360.d0) phasedut(i, n) = phasedut(i, n) - 360.d0
            phasedvt(i, n) = convrd*phasevg(i) + haface(i)
            if (phasedvt(i, n) < 0.d0) phasedvt(i, n) = phasedvt(i, n) + 360.d0
            if (phasedvt(i, n) >= 360.d0) phasedvt(i, n) = phasedvt(i, n) - 360.d0
         end do

         !HARMV...UNCOMMENT THE FOLLOWING LINES TO COMPUTE MEANS AND VARIANCES
         !HARMV...FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
         !HARMV...Resynthesize the time series to compute the average and variances.
         !HARMV...Compare resynthesized values with those computed during time stepping.
         if (charmv) then
            uav(n) = 0.d0
            vav(n) = 0.d0
            usq(n) = 0.d0
            vsq(n) = 0.d0
            do it = 1, ntsteps
               TIMELOC = TIMEBEG + DTDP*dble(IT)
               rsu = 0.d0
               rsv = 0.d0
               do ifr = 1, nfreq + nf
                  ftime = hafreq(ifr)*TIMELOC
                  rsu = rsu + umagg(ifr)*cos(ftime - phaseug(ifr))
                  rsv = rsv + vmagg(ifr)*cos(ftime - phasevg(ifr))
               end do
               uav(n) = uav(n) + rsu
               vav(n) = vav(n) + rsv
               usq(n) = usq(n) + rsu*rsu
               vsq(n) = vsq(n) + rsv*rsv
            end do

            uav(n) = uav(n)/dble(ntsteps)
            vav(n) = vav(n)/dble(ntsteps)
            usq(n) = usq(n)/dble(ntsteps) - uav(n)*uav(n)
            vsq(n) = vsq(n)/dble(ntsteps) - vav(n)*vav(n)
            if (abs(xvelav(n)) <= 0.d0) then
               if (abs(uav(n)) <= 0.d0) uavdif(n) = 1.0d0
               if (abs(uav(n)) > 0.d0) uavdif(n) = 99d19
            else
               uavdif(n) = uav(n)/xvelav(n)
            end if
            if (abs(yvelav(n)) <= 0.d0) then
               if (abs(vav(n)) <= 0.d0) vavdif(n) = 1.0d0
               if (abs(vav(n)) > 0.d0) vavdif(n) = 99d19
            else
               vavdif(n) = vav(n)/yvelav(n)
            end if
            if (abs(xvelva(n)) <= 0.d0) then
               if (abs(usq(n)) <= 0.d0) uvadif(n) = 1.0d0
               if (abs(usq(n)) > 0.d0) uvadif(n) = 99d19
            else
               uvadif(n) = usq(n)/xvelva(n)
            end if
            if (abs(yvelva(n)) <= 0.d0) then
               if (abs(vsq(n)) <= 0.d0) vvadif(n) = 1.0d0
               if (abs(vsq(n)) > 0.d0) vvadif(n) = 99d19
            else
               vvadif(n) = vsq(n)/yvelva(n)
            end if

            if (UAVDIF(n) > UAVMAX) then
               UAVMAX = UAVDIF(n)
               NUAVMAX = n
            end if
            if (UAVDIF(n) < UAVMIN) then
               UAVMIN = UAVDIF(n)
               NUAVMIN = n
            end if
            if (VAVDIF(n) > VAVMAX) then
               VAVMAX = VAVDIF(n)
               NVAVMAX = n
            end if
            if (VAVDIF(n) < VAVMIN) then
               VAVMIN = VAVDIF(n)
               NVAVMIN = n
            end if
            if (UVADIF(n) > UVAMAX) then
               UVAMAX = UVADIF(n)
               NUVAMAX = n
            end if
            if (UVADIF(n) < UVAMIN) then
               UVAMIN = UVADIF(n)
               NUVAMIN = n
            end if
            if (VVADIF(n) > VVAMAX) then
               VVAMAX = VVADIF(n)
               NVVAMAX = n
            end if
            if (VVADIF(n) < VVAMIN) then
               VVAMIN = VVADIF(n)
               NVVAMIN = n
            end if

         end if !  charmv

      end do

      if (charmv) then
         !
         write (16, 7740)
7740     format(///, 5x, 'THE LARGEST VALUES OF THE RATIO ', &
                 'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         write (16, 7743) UAVMAX, NUAVMAX
         write (16, 7744) UVAMAX, NUVAMAX
         write (16, 7745) VAVMAX, NVAVMAX
         write (16, 7746) VVAMAX, NVVAMAX
         write (16, 7747)
7747     format(//, 5x, 'THE LOWEST VALUES OF THE RATIO ', &
                 'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         write (16, 7743) UAVMIN, NUAVMIN
         write (16, 7744) UVAMIN, NUVAMIN
         write (16, 7745) VAVMIN, NVAVMIN
         write (16, 7746) VVAMIN, NVVAMIN
7743     format(9x, ' AVERAGE U VELOCITY RATIO = ', E15.7, ' AT NODE ', I8)
7744     format(9x, 'VARIANCE U VELOCITY RATIO = ', E15.7, ' AT NODE ', I8)
7745     format(9x, ' AVERAGE V VELOCITY RATIO = ', E15.7, ' AT NODE ', I8)
7746     format(9x, 'VARIANCE V VELOCITY RATIO = ', E15.7, ' AT NODE ', I8)
         !
      end if ! charmv
      !

      return
   end subroutine lsqsolvg

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  R E A D  B I N A R Y  H A  H O T S T A R T
   !--------------------------------------------------------------------
   !     jgf49.44: Reads harmonic analysis data from a binary hotstart
   !     file. Based on HAHOTS, written by rl.
   !     WJP 02.20.2018: Allowing for > 1 options for netcdf output
   !--------------------------------------------------------------------
   subroutine readBinaryHAHotstart(lun, counter)
      use MOD_SIZES, only: MNPROC, READ_LOCAL_HOT_START_FILES
      use GLOBAL, only: NSTAE, NSTAE_G, NSTAV, NSTAV_G, &
                        ITHS, IMAP_STAE_LG, NP_G, NODES_LG, &
                        IMAP_STAV_LG
      use MESH, only: NP
      implicit none
      integer, intent(in) :: lun ! i/o logical unit number
      integer, intent(inout) :: counter ! i/o record
      integer :: i, j, n ! loop counters
      integer :: num_elev_sta ! number of elevation stations to read
      integer :: num_vel_sta ! number of velocity stations to read
      integer :: num_nodes ! number of nodes to read
      integer :: subdomainStation ! number of station to map data to
      integer :: fulldomainStation ! number of station to map data from
      integer :: subdomainNode ! number of node to map data to
      integer :: fulldomainNode ! number of node to map data from
      !
      LOG_SCOPE_TRACED("readBinaryHAHotstart", HARM_TRACING)
      !
      !     Settings for reading in hotstart file in serial, or for reading
      !     subdomain hotstart files in parallel.
      if ((MNPROC == 1) &
          .or. (READ_LOCAL_HOT_START_FILES .eqv. .true.)) then
         STAELV_pt => STAELV
         num_elev_sta = NSTAE
         STAULV_pt => STAULV
         STAVLV_pt => STAVLV
         num_vel_sta = NSTAV
         GLOELV_pt => GLOELV
         num_nodes = NP
         GLOULV_pt => GLOULV
         GLOVLV_pt => GLOVLV
         ELAV_pt => ELAV
         ELVA_pt => ELVA
         XVELAV_pt => XVELAV
         YVELAV_pt => YVELAV
         XVELVA_pt => XVELVA
         YVELVA_pt => YVELVA
      else
         ! read full domain hotstart file in parallel
         STAELV_pt => STAELV_G
         num_elev_sta = NSTAE_G
         STAULV_pt => STAULV_G
         STAVLV_pt => STAVLV_G
         num_vel_sta = NSTAV_G
         GLOELV_pt => GLOELV_G
         num_nodes = NP_G
         GLOULV_pt => GLOULV_G
         GLOVLV_pt => GLOVLV_G
         ELAV_pt => ELAV_g
         ELVA_pt => ELVA_g
         XVELAV_pt => XVELAV_g
         YVELAV_pt => YVELAV_g
         XVELVA_pt => XVELVA_g
         YVELVA_pt => YVELVA_g
      end if
      !
      !***** Read in parameter values
      !
      read (lun, REC=counter + 1) inz
      read (lun, REC=counter + 2) inf
      read (lun, REC=counter + 3) imm
      read (lun, REC=counter + 4) inp
      read (lun, REC=counter + 5) instae
      read (lun, REC=counter + 6) instav
      read (lun, REC=counter + 7) inhase
      read (lun, REC=counter + 8) inhasv
      read (lun, REC=counter + 9) inhage
      read (lun, REC=counter + 10) inhagv
      read (lun, REC=counter + 11) iicall
      read (lun, REC=counter + 12) infreq
      counter = counter + 12
      !
      do i = 1, nfreq + nf
         read (lun, REC=counter + 1) FNAM8(1)
         read (lun, REC=counter + 2) FNAM8(2)
         counter = counter + 2
         FNAME = FNAM8(1)//FNAM8(2)
         INAMEFR(I) = FNAME(1:10)
         read (lun, REC=counter + 1) ifreq(i)
         read (lun, REC=counter + 2) iff(i)
         read (lun, REC=counter + 3) iface(i)
         counter = counter + 3
      end do
      !
      !***** Read in time of most recent H.A. update
      !
      read (lun, REC=counter + 1) TIMEUD
      read (lun, REC=counter + 2) ITUD
      counter = counter + 2
      !
      !     Read in LHS Matrix
      counter = counter + 1

      ! jgf50.96: do not use binaryRead3D(HA, mm, mm, lun, counter)
      ! b/c HA is written with an inner j loop while binaryRead3D is
      ! written with an inner i loop
      do i = 1, mm
         do j = 1, mm
            read (lun, REC=counter) HA(I, J)
            counter = counter + 1
         end do
      end do
      !
      !
      !     Read in Station Elevation LHS load vector
      if (NHASE /= 0) then
         call binaryRead3D(STAELV_pt, mm, num_elev_sta, lun, counter)
      end if
      !     Read in Station Velocity LHS load vector
      if (NHASV /= 0) then
         do n = 1, num_vel_sta
            do i = 1, mm
               read (lun, REC=counter) STAULV_pt(I, N)
               read (lun, REC=counter + 1) STAVLV_pt(I, N)
               counter = counter + 2
            end do
         end do
      end if
      !     Read in Global Elevation LHS load vector
      if (NHAGE /= 0) then
         call binaryRead3D(GLOELV_pt, mm, num_nodes, lun, counter)
      end if
      !     Read in Global Velocity LHS load vector
      if (NHAGV /= 0) then
         do n = 1, num_nodes
            do i = 1, mm
               read (lun, REC=counter) GLOULV_pt(I, N)
               read (lun, REC=counter + 1) GLOVLV_pt(I, N)
               counter = counter + 2
            end do
         end do
      end if
      !..   Read in Means and Squares
      if (CHARMV .eqv. .true.) then
         if ((abs(FMV) > 0.d0) .and. (ITHS > ITMV)) then
            read (lun, REC=counter) NTSTEPS; counter = counter + 1
            if (NHAGE /= 0) then
               do I = 1, num_nodes
                  read (lun, REC=counter) ELAV_pt(I)
                  counter = counter + 1
                  read (lun, REC=counter) ELVA_pt(I)
                  counter = counter + 1
               end do
            end if
            if (NHAGV /= 0) then
               do I = 1, num_nodes
                  read (lun, REC=counter) XVELAV_pt(I)
                  counter = counter + 1
                  read (lun, REC=counter) YVELAV_pt(I)
                  counter = counter + 1
                  read (lun, REC=counter) XVELVA_pt(I)
                  counter = counter + 1
                  read (lun, REC=counter) YVELVA_pt(I)
                  counter = counter + 1
               end do
            end if
         end if
      end if !  charmv
      !
      !     Map data from fulldomain to subdomain in parallel
      if ((MNPROC > 1) &
          .and. (READ_LOCAL_HOT_START_FILES .eqv. .false.)) then
         if (NHASE /= 0) then
            do subdomainStation = 1, NSTAE
               fulldomainStation = abs(IMAP_STAE_LG(subdomainStation))
               do i = 1, mm
                  STAELV(i, subdomainStation) = &
                     STAELV_g(i, fulldomainStation)
               end do
            end do
         end if
         if (NHASV /= 0) then
            do subdomainStation = 1, NSTAV
               fulldomainStation = abs(IMAP_STAV_LG(subdomainStation))
               do i = 1, mm
                  STAULV(i, subdomainStation) = &
                     STAULV_g(i, fulldomainStation)
                  STAVLV(i, subdomainStation) = &
                     STAVLV_g(i, fulldomainStation)
               end do
            end do
         end if
         if (NHAGE /= 0) then
            do subdomainNode = 1, NP
               fulldomainNode = abs(NODES_LG(subdomainNode))
               do i = 1, mm
                  GLOELV(i, subdomainNode) = GLOELV_g(i, fulldomainNode)
               end do
            end do
         end if
         if (NHAGV /= 0) then
            do subdomainNode = 1, NP
               fulldomainNode = abs(NODES_LG(subdomainNode))
               do i = 1, mm
                  GLOULV(i, subdomainNode) = GLOULV_g(i, fulldomainNode)
                  GLOVLV(i, subdomainNode) = GLOVLV_g(i, fulldomainNode)
               end do
            end do
         end if

         !..      Map Means and Squares
         if (CHARMV .eqv. .true.) then
            if ((abs(FMV) > 0.d0) .and. (ITHS > ITMV)) then
               if (NHAGE /= 0) then
                  do subdomainNode = 1, NP
                     fulldomainNode = abs(NODES_LG(subdomainNode))
                     ELAV(subdomainNode) = ELAV_g(fulldomainNode)
                     ELVA(subdomainNode) = ELVA_g(fulldomainNode)
                  end do
               end if
               if (NHAGV /= 0) then
                  do subdomainNode = 1, NP
                     fulldomainNode = abs(NODES_LG(subdomainNode))
                     XVELAV(subdomainNode) = XVELAV_g(fulldomainNode)
                     YVELAV(subdomainNode) = YVELAV_g(fulldomainNode)
                     XVELVA(subdomainNode) = XVELVA_g(fulldomainNode)
                     YVELVA(subdomainNode) = YVELVA_g(fulldomainNode)
                  end do
               end if
            end if
         end if
      end if ! MNPROC.gt.1 etc

      !--------------------------------------------------------------------
   end subroutine readBinaryHAHotstart
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E  C H E C K   H A  H O T S T A R T
   !--------------------------------------------------------------------
   !     jgf49.44: Checks harmonic analysis data from a hotstart file
   !     by comparing with parameters read from the fort.15.
   !     Based on HAHOTS, written by rl.
   !--------------------------------------------------------------------
   subroutine checkHAHotstart()
      use MOD_SIZES, only: READ_LOCAL_HOT_START_FILES, MYPROC
      use GLOBAL, only: NSTAE_G, NSTAV_G, NSTAE, NSTAV, &
                        NP_G
      use MESH, only: NP
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      integer :: i ! loop counter
      real(8) :: FDIFF ! difference between frequencies btw hotstart and fort.15
      character(len=1024) :: message
      !
      LOG_SCOPE_TRACED("checkHAHotstart", HARM_TRACING)
      call logMessage(INFO, &
                      "Comparing harmonic parameters in hotstart file to fort.15...")
      !
      iflag = 0
      if (nz /= inz) iflag = 1
      if (nf /= inf) iflag = 1
      if (mm /= imm) iflag = 1
      if ((myProc == 0) .and. (read_local_hot_start_files .eqv. .false.)) then
         if (np_g /= inp) iflag = 1
         if (nstae_g /= instae) iflag = 1
         if (nstav_g /= instav) iflag = 1
      else if (READ_LOCAL_HOT_START_FILES .eqv. .true.) then
         if (np /= inp) iflag = 1
         if (nstae /= instae) iflag = 1
         if (nstav /= instav) iflag = 1
      end if
      if (nhase /= inhase) iflag = 1
      if (nhasv /= inhasv) iflag = 1
      if (nhage /= inhage) iflag = 1
      if (nhagv /= inhagv) iflag = 1
      if (nfreq /= infreq) iflag = 1
      !
      do i = 1, nfreq + nf
         if (namefr(i) /= inamefr(i)) iflag = 1
         if (abs(hafreq(i) + ifreq(i)) < 1.0d-30) then
            fdiff = 0.d0
         else
            fdiff = abs(hafreq(i) - ifreq(i))/abs(hafreq(i) + ifreq(i))
         end if
         if (fdiff >= 1.d-6) iflag = 1
         if (abs(HAFF(i) + iFF(i)) < 1d-30) then
            fdiff = 0.d0
         else
            fdiff = abs(HAFF(i) - iFF(i))/abs(HAFF(i) + iFF(i))
         end if
         if (fdiff >= 1.d-6) iflag = 1
         if (abs(HAFACE(i) + iFACE(i)) < 1d-30) then
            fdiff = 0.d0
         else
            fdiff = abs(HAFACE(i) - iFACE(i))/abs(HAFACE(i) + iFACE(i))
         end if
         if (fdiff >= 1.d-6) iflag = 1
      end do
      !
      if (iflag == 0) then
         call logMessage(INFO, &
                         "Harmonic data in hotstart file matches fort.15. PASSED.")
      else
         !
         !***** FATAL Error Messages
         !
         write (message, 1000)
         call allMessage(ERROR, message)
1000     format('***** DISCREPANCY IN HARMONIC ANALYSIS HOT ', &
                'START FILE *****')
         if (nz /= inz) then
            write (message, 2010) nz, inz
            call allMessage(ERROR, message)
2010        format('NZ COMPUTED FROM UNIT 14 INPUT = ', I2, &
                   ', NZ READ IN FROM HOT START FILE = ', I2)
         end if
         if (nf /= inf) then
            write (message, 2020) nf, inf
            call allMessage(ERROR, message)
2020        format('NF COMPUTED FROM UNIT 14 INPUT = ', I2, &
                   ', NF READ IN FROM HOT START FILE = ', I2)
         end if
         if (mm /= imm) then
            write (message, 2030) mm, imm
            call allMessage(ERROR, message)
2030        format('MM COMPUTED FROM UNIT 14 INPUT = ', I2, &
                   ', MM READ IN FROM HOT START FILE = ', I2)
         end if
         if (myproc == 0) then
            if (np_g /= inp) then
               write (message, 2040) np_g, inp
               call allMessage(ERROR, message)
2040           format('NP READ IN FROM UNIT 14 = ', I2, &
                      ', NP READ IN FROM HOT START FILE = ', I2)
            end if
            if (nstae_g /= instae) then
               write (message, 2050) nstae_g, instae
               call allMessage(ERROR, message)
2050           format('NSTAE READ IN FROM UNIT 15 = ', I2, &
                      ', NSTAE READ IN FROM HOT START FILE = ', I2)
            end if
            if (nstav_g /= instav) then
               write (message, 2060) nstav_g, instav
               call allMessage(ERROR, message)
2060           format('NSTAV READ IN FROM UNIT 15 = ', I2, &
                      ', NSTAV READ IN FROM HOT START FILE = ', I2)
            end if
         end if
         if (nhase /= inhase) then
            write (message, 2070) NHASE, INHASE
            call allMessage(ERROR, message)
2070        format('NHASE READ IN FROM UNIT 15 = ', I2, &
                   ', NHASE READ IN FROM HOT START FILE = ', I2)
         end if
         if (nhasv /= inhasv) then
            write (message, 2080) NHASV, INHASV
            call allMessage(ERROR, message)
2080        format('NHASV READ IN FROM UNIT 15 = ', I2, &
                   ', NHASV READ IN FROM HOT START FILE = ', I2)
         end if
         if (nhage /= inhage) then
            write (message, 2090) NHAGE, INHAGE
            call allMessage(ERROR, message)
2090        format('NHAGE READ IN FROM UNIT 15 = ', I2, &
                   ', NHAGE READ IN FROM HOT START FILE = ', I2)
         end if
         if (nhagv /= inhagv) then
            write (message, 2100) NHAGV, INHAGV
            call allMessage(ERROR, message)
2100        format('NHAGV READ IN FROM UNIT 15 = ', I2, &
                   ', NHAGV READ IN FROM HOT START FILE = ', I2)
         end if
         if (nfreq /= infreq) then
            write (message, 2110) NFREQ, INFREQ
            call allMessage(ERROR, message)
2110        format('NFREQ COMPUTED FROM UNIT 15 INPUT = ', I2, &
                   ', NFREQ READ IN FROM HOT START FILE = ', I2)
         end if
         do i = 1, nfreq + nf
            if (namefr(i) /= inamefr(i)) then
               write (message, 2120) i, namefr(i), inamefr(i)
               call allMessage(ERROR, message)
2120           format('FOR CONSTITUENT # ', I3, &
                      ', NAMEFR READ IN FROM UNIT 15 = ', A10, &
                      ', NAMEFR READ IN FROM HOT START FILE = ', A10)
            end if
            if (abs(hafreq(i) - ifreq(i)) > 0.d0) then
               write (message, 2130) i, hafreq(i), ifreq(i)
               call allMessage(ERROR, message)
2130           format('FOR CONSTITUENT # ', I3, &
                      ', FREQ READ IN FROM UNIT 15 = ', D20.10, &
                      ', FREQ READ IN FROM HOT START FILE = ', D20.10)
            end if
            if (abs(HAFF(i) - iFF(i)) > 0.d0) then
               write (message, 2140) i, haff(i), iff(i)
               call allMessage(ERROR, message)
2140           format('FOR CONSTITUENT # ', I3, &
                      ', FF READ IN FROM UNIT 15 = ', F10.5, &
                      ', FF READ IN FROM HOT START FILE = ', F10.5)
            end if
            if (abs(HAFACE(i) - iFACE(i)) > 0.d0) then
               write (message, 2150) i, haface(i), iface(i)
               call allMessage(ERROR, message)
2150           format('FOR CONSTITUENT # ', I3, &
                      ', FACE READ IN FROM UNIT 15 = ', F10.5, &
                      ', FACE READ IN FROM HOT START FILE = ', F10.5)
            end if
         end do
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='Discrepancy in harmonic analysis ' &
                        //'hotstart file parameters')
      end if

      !--------------------------------------------------------------------
   end subroutine checkHAHotstart
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     S U B R O U T I N E      B I N A R Y   R E A D   3 D
   !--------------------------------------------------------------------
   !     jgf49.44: Same as binaryRead2D except it reads data with
   !     both horizontal and vertical lengths. Also used for reading
   !     harmonic analysis arrays from binary hotstart file.
   !--------------------------------------------------------------------
   subroutine binaryRead3D(array, iLength, jLength, lun, &
                           counter)
      use mod_logging, only: t_log_scope, &
                             allMessage
      implicit none
      integer, intent(in) :: iLength ! the number of horiz values to read
      integer, intent(in) :: jLength ! the number of layers
      real(8), intent(out), dimension(ilength, jlength) :: array ! data to read from the hotstart file
      integer, intent(in) :: lun ! fortran logical unit number to read from
      integer, intent(inout) :: counter ! i/o position in the file
      !
      integer :: i, j ! array indices
      !
      LOG_SCOPE_TRACED("binaryRead3D", HARM_TRACING)
      do j = 1, jLength
         do i = 1, iLength
            read (lun, REC=counter) array(i, j)
            counter = counter + 1
         end do
      end do
      !
      !--------------------------------------------------------------------
   end subroutine binaryRead3D
   !--------------------------------------------------------------------

   !--------------------------------------------------------------
end module mod_harmonic_analysis
!--------------------------------------------------------------
