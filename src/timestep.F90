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
!-----------------------------------------------------------------------
!  MODULE TIMESTEP
!-----------------------------------------------------------------------
!> @brief Main timestepping module for ADCIRC hydrodynamic calculations
!>
!> This module provides the primary timestepping routine for ADCIRC,
!> orchestrating the solution of the GWCE and momentum equations. It handles
!> time level shifts, ramp functions, boundary conditions, meteorological
!> forcing, wetting/drying, baroclinic pressure gradients, and output.
!>
!> The timestepping module supports multiple algorithm configurations:
!> - C2DDI: 2D Depth Integrated model
!> - C3D: 3D model
!> - Baroclinic transport calculations
!> - Predictor-corrector
!-----------------------------------------------------------------------
module mod_timestep

   implicit none

   private

   public :: TIMESTEP

   real(8) :: hotstart_start_time = 0.0d0
   real(8) :: hotstart_reference_time = 0.0d0

   real(8), parameter :: eps = epsilon(1.0d0)

contains

!-----------------------------------------------------------------------
!> @brief Main timestepping subroutine for ADCIRC hydrodynamic model
!>
!> This subroutine performs a single timestep of the ADCIRC model, including:
!> - Computing time variables and ramp functions
!> - Shifting solution arrays to previous time levels
!> - Applying meteorological forcing (wind, pressure)
!> - Computing boundary conditions (elevation, flux, barriers)
!> - Solving the GWCE for water surface elevation
!> - Applying wetting and drying algorithm
!> - Solving 2D or 3D momentum equations
!> - Applying absorbing sponge layers
!> - Collecting output statistics and writing hot starts
!>
!> Algorithm selection is controlled by logical flags set in READ_INPUT
!> based on the IM parameter:
!> - C2DDI: 2D Depth Integrated model
!> - C3D: 3D model (C3DDSS for stress form, C3DVS for velocity form)
!> - C2D_BTrans/C3D_BTrans: Baroclinic transport
!> - C2D_PTrans/C3D_PTrans: Passive transport
!> - CBaroclinic: Include baroclinic terms
!> - CPRECOR: Predictor-corrector algorithm
!>
!> @param[in] IT        Current timestep number
!> @param[in] ITIME_BGN Starting timestep number (for hot starts)
!> @param[out] TimeLoc  Simulation time at current timestep (seconds)
!-----------------------------------------------------------------------
   subroutine TIMESTEP(IT, ITIME_BGN, TimeLoc)

      use GLOBAL, only: METONLY, NDDT, CTIP, C2DDI, C3DVS, C2D_PTrans, NT, &
                        NWS, WSX2, WSY2, PR2, WVNXOUT, WVNYOUT, &
                        usingDynamicWaterLevelCorrection
      use GWCE, only: solveGWCE
      use MOMENTUM2D, only: solve2DMomentumEq
      use MOMENTUM3D, only: solve3DVSMomentumEq
      use WIND, only: CLOSE_MET_FILES, getMeteorologicalForcing
      use mod_logging, only: t_log_scope, init_log_scope
      use SPONGELAYER, only: updateAbsorbingLayerBoundary, applyAbsorbingSpongeLayer
      use mod_baroclinic, only: computeBaroclinicPressureGradient, updateBaroclinicBoundaryConditions
      use mod_flux_boundaries, only: applyFluxBoundaryConditions
      use mod_transport_2d, only: SCALAR_TRANS_2D

      implicit none
      integer, intent(in) :: IT
      integer, intent(in) :: ITIME_BGN
      real(8), intent(out) :: TimeLoc
      real(8) :: TimeH

      LOG_SCOPE_TRACED("timestep", TIMESTEP_TRACING)

      call computeTimeVariables(IT, TimeLoc, TimeH)
      call computeRampFunctions(IT, TimeLoc)
      call shiftTimeLevels(IT)
      call updateAbsorbingLayerBoundary(TimeLoc, TimeH)
      if (NDDT /= 0) call updateTimeDependentBathymetry(TimeLoc)
      call computeBottomFriction(IT)
      call updateIceFields(TimeLoc)
      if (usingDynamicWaterLevelCorrection) call updateWaterLevelCorrections(TimeLoc)
      call getMeteorologicalForcing(nws, timeloc, wsx2, wsy2, pr2, wvnxout, wvnyout)
      call updateWaveRadiationStress(TimeLoc)

      ! Skip the GWCE and momentum solutions if only meteorological output
      ! was requested
      if (.not. METONLY) then

         call updateBaroclinicBoundaryConditions(TimeLoc)
         if (CTIP) call applyTidePotentialForcing(TimeLoc, TimeH)
         call computeBaroclinicPressureGradient(TimeLoc)

         call applyFluxBoundaryConditions(IT, TimeLoc, TimeH)

         call solveGWCE(it, ITIME_BGN, timeloc, timeh)

         call applyWettingAndDrying(IT)

         call updateFlowDepth()

         if (C2DDI) call solve2DMomentumEq(IT, TimeLoc)

         call applyAbsorbingSpongeLayer()

         if (C3DVS) call solve3DVSMomentumEq(IT, TimeLoc)

         if (C2D_PTrans) call SCALAR_TRANS_2D()

      end if

      call collectOutputStatistics(IT, TimeLoc, TimeH)
      call writeHotstartIfNeeded(IT, TimeLoc)
      call shiftBarrierHeights()

      call progressScreen(IT, TimeLoc)

      if (IT == NT) call CLOSE_MET_FILES()

   end subroutine TIMESTEP
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute tidal potential forcing at mesh nodes
!>
!> Computes the gravitational tide-generating potential at each node,
!> including Earth tide reduction factors. Supports both the full formula
!> (tidePotential%active) and the traditional constituent-based approach.
!>
!> For the traditional approach, tidal potential is computed as a sum of
!> constituents with amplitudes, frequencies, and phase corrections. The
!> self-attraction and loading (SAL) terms are also included.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine applyTidePotentialForcing(TimeLoc, TimeH)
      use mod_tidepotential, only: tidePotential
      use global, only: NTIP, TIP2, NTIF, RampTip, PERT, FACET, &
                        AMIGT, FFT, TPK, ETRF, SALTPHA, SALTAMP, L_N
      use mesh, only: np, SLAM
      implicit none
      real(8), intent(in) :: TimeLoc, TimeH
      real(8) :: ArgT, ArgTP, ArgSAlt, SALTMUL, TPMUL
      integer :: i, j, na, ncyc

      ! use the full formula
      if (tidePotential%active()) then
         TIP2 = tidePotential%compute(TimeLoc, NP, SLAM)*RampTip
      end if

      ! The traditional approach: a sum of different constituents
      ! TIP & SAL - for tidePotential%active() = .TRUE.
      ! SAL       - for tidePotential%active() = .FALSE.
      if (NTIP == 2 .or. .not. tidePotential%active()) then
         do J = 1, NTIF
            if (abs(PERT(J)) <= eps) then
               NCYC = 0
            else
               NCYC = int(timeh/PERT(J))
            end if

            ARGT = AMIGT(J)*(timeh - dble(NCYC)*PERT(J)) + FACET(J)
            TPMUL = RampTip*ETRF(J)*TPK(J)*FFT(J)
            SALTMUL = RampTip*FFT(J)

            !! WJP: We actually want to compare against diurnal and get 0,1,2
            !!      or higher species (was opposite beforehand)
            NA = min(nint(AMIGT(J)/7d-5), 2)

            !! WJP: Rewritten so that we have a general formula for all
            !!      species
            if (tidePotential%active()) then
               do I = 1, NP
                  ARGSALT = ARGT - SALTPHA(J, I)
                  TIP2(I) = TIP2(I) + SALTMUL*SALTAMP(J, I)*cos(ARGSALT)
               end do
            else
               do I = 1, NP
                  ARGTP = ARGT + dble(NA)*SLAM(I)
                  ARGSALT = ARGT - SALTPHA(J, I)
                  TIP2(I) = TIP2(I) + TPMUL*L_N(NA, I)*cos(ARGTP) &
                            + SALTMUL*SALTAMP(J, I)*cos(ARGSALT)
               end do
            end if

         end do
      end if
   end subroutine applyTidePotentialForcing
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Shift solution arrays to previous time levels
!>
!> Shifts flux per unit width, depth averaged velocities, bottom stress,
!> wind stress, surface pressure and tidal potentials to previous time step.
!> Zeros out new forcing terms and load vectors. On the first timestep,
!> initializes H1 and H2 from bathymetry and initial water surface elevation.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine shiftTimeLevels(IT)
      use global, only: MOM_LV_X, MOM_LV_Y, QX1, QY1, QX2, QY2, &
                        UU1, VV1, UU2, VV2, CPRECOR, UU0, VV0, QX0, QY0, &
                        IM, TRANS_LV_B, TRANS_LV_A, NWS, NRS, &
                        WSX1, WSY1, WSX2, WSY2, PR1, PR2, CTIP, TIP1, TIP2, &
                        usingDynamicWaterLevelCorrection, &
                        dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, &
                        H1, H2, ETA2, IFNLFA, &
                        QN0, QN1, QN2, EN0, EN1, EN2
      use gwce, only: gwce_lv_reset
      use wind, only: PRBCKGRND_MH2O
      use mesh, only: DP
      use boundaries, only: NVEL
      use mod_subgrid, only: subgrid_level0, gridDepthVertETA2
      implicit none

      integer, intent(in) :: IT

      ! On first timestep of a cold start, initialize H1 and H2
      ! (the grid-averaged depth H_G under subgrid level0)
      if (IT == 1) then
         if (subgrid_level0) then
            H2(:) = gridDepthVertETA2(:)
         else
            H2(:) = DP(:) + dble(IFNLFA)*ETA2(:)
         end if
         H1(:) = H2(:)
      end if

      ! Shift values in time for predictor-corrector algorithm. This must
      ! happen before level 1 is overwritten by level 2 below.
      if (CPRECOR) then
         UU0(:) = UU1(:)
         VV0(:) = VV1(:)
         QX0(:) = QX1(:)
         QY0(:) = QY1(:)
      end if

      MOM_LV_X(:) = 0d0
      MOM_LV_Y(:) = 0d0
      QX1(:) = QX2(:)
      QY1(:) = QY2(:)
      UU1(:) = UU2(:)
      VV1(:) = VV2(:)
      call gwce_lv_reset(0d0)

      ! Transport
      if (IM == 10) then
         TRANS_LV_B(:) = 0.d0
         TRANS_LV_A(:) = 0.d0
      end if

      ! Wind (& wave radiation stress if used)
      if ((NWS /= 0) .or. (NRS /= 0)) then
         WSX1(:) = WSX2(:)
         WSX2(:) = 0.d0
         WSY1(:) = WSY2(:)
         WSY2(:) = 0.d0
         PR1(:) = PR2(:)
         PR2(:) = PRBCKGRND_MH2O
      else
         PR2(:) = PRBCKGRND_MH2O
      end if

      ! Tidal potential forcing
      if (CTIP) then
         TIP1(:) = TIP2(:)
         TIP2(:) = 0.d0
      end if

      if (usingDynamicWaterLevelCorrection) then
         dynamicWaterLevelCorrection1(:) = dynamicWaterLevelCorrection2(:)
         dynamicWaterLevelCorrection2(:) = 0.d0
      end if

      ! Normal flow boundary conditions
      QN0(1:NVEL) = QN1(1:NVEL)
      QN1(1:NVEL) = QN2(1:NVEL)
      QN2(1:NVEL) = 0.d0
      EN0(1:NVEL) = EN1(1:NVEL)
      EN1(1:NVEL) = EN2(1:NVEL)
      EN2(1:NVEL) = 0.d0

   end subroutine shiftTimeLevels
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Save the current barrier heights into the previous time level
!>
!> Called at the end of the timestep, after the hot start file has been
!> written.
!-----------------------------------------------------------------------
   subroutine shiftBarrierHeights()
      use boundaries, only: NFLUXIB, NFLUXB
      use weir, only: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2
      implicit none

      if (NFLUXIB == 1) BARINHT1(:) = BARINHT2(:)
      if (NFLUXB == 1) BARLANHT1(:) = BARLANHT2(:)

   end subroutine shiftBarrierHeights
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Calculate temporal ramp functions for forcing terms
!>
!> Computes hyperbolic tangent ramp functions for gradual introduction of
!> boundary elevation, wind/pressure, tidal potential, and wave radiation
!> stress forcing. Handles special cases for flux settling periods.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeRampFunctions(IT, TimeLoc)
      use global, only: NRamp, Ramp, RampExtFlux, RampIntFlux, &
                        RampElev, RampTip, RampMete, RampWRad, &
                        DRamp, DRampExtFlux, DRampIntFlux, DRampElev, &
                        DRampTip, DRampMete, DRampWRad, DUnRampMete, &
                        FluxSettlingIT, DTDP, STATIM, CHOTHS
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      real(8) :: TimeLocFluxSettling

      !     jgf46.08 Calculate ramp functions.
      !     jgf46.21 Modify to match behavior of 46.02
      if (NRamp == 0) then
         Ramp = 1.0d0
         RampExtFlux = 1.0d0
         RampIntFlux = 1.0d0
         RampElev = 1.0d0
         RampTip = 1.0d0
         RampMete = 1.0d0
         RampWRad = 1.0d0
      else
         Ramp = tanh((2.d0*TimeLoc/86400.d0)/DRamp)
         RampExtFlux = tanh((2.d0*TimeLoc/86400.d0)/DRampExtFlux)
         RampIntFlux = tanh((2.d0*TimeLoc/86400.d0)/DRampIntFlux)
         RampElev = tanh((2.d0*TimeLoc/86400.d0)/DRampElev)
         RampTip = tanh((2.d0*TimeLoc/86400.d0)/DRampTip)
         RampMete = tanh((2.d0*(TimeLoc/86400.d0 - DUnRampMete))/DRampMete)
         RampWRad = tanh((2.d0*TimeLoc/86400.d0)/DRampWRad)
      end if

      !     jgf46.21 If there is an external flux (i.e. river) boundary, turn
      !     off all forcings except the river flux forcing for the duration of
      !     the FluxSettlingTime. When the FluxSettlingTime has ended, turn
      !     all forcings back on.
      if (NRamp > 1 .and. FluxSettlingIT > 0) then
         if (IT < (FluxSettlingIT + 10)) then
            Ramp = 0.0d0
            RampIntFlux = 0.0d0
            RampElev = 0.0d0
            RampTip = 0.0d0
            RampMete = 0.0d0
            RampWRad = 0.0d0
         else
            TimeLocFluxSettling = dble(IT - FluxSettlingIT - 10)*DTDP &
                                  + STATIM*86400.d0
            if (CHOTHS) then
               TimeLocFluxSettling = dble(IT - FluxSettlingIT - 10)*DTDP &
                                     + (hotstart_start_time - hotstart_reference_time)*86400.d0
            end if

            Ramp = tanh((2.d0*TimeLocFluxSettling/86400.d0)/DRamp)
            RampIntFlux = tanh((2.d0 &
                                *TimeLocFluxSettling/86400.d0)/DRampIntFlux)
            RampElev = tanh((2.d0 &
                             *TimeLocFluxSettling/86400.d0)/DRampElev)
            RampTip = tanh((2.d0 &
                            *TimeLocFluxSettling/86400.d0)/DRampTip)
            if (NRamp /= 8) then
               RampMete = tanh((2.d0 &
                                *TimeLocFluxSettling/86400.d0)/DRampMete)
            end if
            RampWRad = tanh((2.d0 &
                             *TimeLocFluxSettling/86400.d0)/DRampWRad)
         end if
      end if

      ! jgf49.44: Cover the case where the ramp length is zero.
      if (DRamp < 1.0d-6) Ramp = 1.0d0
      if (DRampExtFlux < 1.0d-6) RampExtFlux = 1.0d0
      if (DRampIntFlux < 1.0d-6) RampIntFlux = 1.0d0
      if (DRampElev < 1.0d-6) RampElev = 1.0d0
      if (DRampTip < 1.0d-6) RampTip = 1.0d0
      if (DRampMete < 1.0d-6) RampMete = 1.0d0
      if (DRampWRad < 1.0d-6) RampWRad = 1.0d0

   end subroutine computeRampFunctions
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute time variables for current timestep
!>
!> Computes the master simulation time (TimeLoc), harmonic analysis
!> reference time (TimeH), and flux settling iteration count. Handles
!> timestep changes when hot starting from a different timestep size.
!>
!> @param[in]  IT      Current timestep number
!> @param[out] TimeLoc Simulation time at current timestep (seconds)
!> @param[out] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine computeTimeVariables(IT, TimeLoc, TimeH)
      use global, only: CHOTHS, FluxSettlingTime, FluxSettlingIT, &
                        DTDPHS, DTDP, STATIM, REFTIM, ITHS
      implicit none

      integer, intent(in) :: IT
      real(8), intent(out) :: TimeLoc, TimeH

      !     kmd48.33bc - changed for timestep changes in the hot start files
      !     jgf46.21 Combined flux/radiation b.c. for rivers
      if (CHOTHS) then
         FluxSettlingIT = int(FluxSettlingTime*86400.d0/DTDPHS)
      else
         FluxSettlingIT = int(FluxSettlingTime*86400.d0/DTDP)
      end if

      !     Compute master time referenced to beginning of model run
      !     TimeLoc: master simulation time
      !     TimeH: time for harmonic calculations (includes REFTIM offset)
      TimeLoc = dble(IT)*DTDP + STATIM*86400.d0
      TimeH = dble(IT)*DTDP + (STATIM - REFTIM)*86400.d0

      if (CHOTHS) then
         if ((ITHS + 1) == IT) then
            hotstart_start_time = (dble(IT - 1)*DTDPHS)/86400.d0
            hotstart_reference_time = (dble(IT - 1)*DTDP)/86400.d0
         end if
         TimeLoc = dble(IT)*DTDP + (hotstart_start_time - hotstart_reference_time)*86400.d0
         TimeH = dble(IT)*DTDP + ((hotstart_start_time - STATIM) - REFTIM)*86400.d0
      end if

      !     Harmonic time recalculation (includes hot start timestep changes)
      TimeH = dble(IT)*DTDP + (STATIM - REFTIM)*86400.d0
      if (CHOTHS) then
         if ((ITHS + 1) == IT) then
            hotstart_start_time = (dble(IT - 1)*DTDPHS)/86400.d0
            STATIM = (dble(IT - 1)*DTDP)/86400.d0
         end if
         TimeH = dble(IT)*DTDP + ((hotstart_start_time - STATIM) - REFTIM)*86400.d0
      end if

   end subroutine computeTimeVariables
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update bathymetry for time-dependent bathymetry runs
!>
!> Updates mesh bathymetry (DP) from input file for simulations with
!> time-varying bathymetry. Supports two modes:
!> - NDDT=1: Read bathymetry for all nodes
!> - NDDT=2: Read bathymetry for sparse node subset
!>
!> Bathymetry is linearly interpolated between snapshots:
!> @verbatim
!>    DP1                         DP2
!>  btime1                      btime2      btime2 = btime1 + btiminc
!>     |---------x-----------------|        btime_end = btime1 + bchgtiminc
!>            btime_end
!> @endverbatim
!>
!> Water surface elevation is adjusted to maintain total water column depth.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateTimeDependentBathymetry(TimeLoc)
      use global, only: NDDT, BTIME1, BTIME2, BTIMINC, BCHGTIMINC, &
                        BTIME_END, dp1, dp2, ETA2, NOLIFA, H0
      use mesh, only: NP, DP
      use mod_logging, only: logMessage, allMessage, INFO
      use global_io, only: NDDT2GET
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: bTRatio, DPTMP, DPTMP2
      integer :: I, IJK
      character(256) :: scratchMessage

      if (abs(NDDT) == 1) then
         ! Mode 1: Read bathymetry for all nodes
         if (TimeLoc > BTIME2) then
            BTIME1 = BTIME2
            BTIME2 = BTIME2 + BTIMINC
            BTIME_END = BTIME1 + BCHGTIMINC

            dp1(:) = dp2(:)
            ! Read new record for all nodes
            do I = 1, NP
               read (141, *) IJK, DP2(IJK)
            end do
            ! Ensure depths >= H0 if no wetting/drying
            if (NOLIFA == 0 .or. NOLIFA == 1) then
               where (DP2(1:NP) < H0) DP2(1:NP) = H0
            end if
            write (scratchMessage, '(A36,1X,E15.8,1X,A5)') &
               'BATHYMETRY RECORDS UPDATED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
         end if

         ! Interpolate during bathymetry change interval
         if (TimeLoc < BTIME_END) then
            bTRatio = (TimeLoc - BTIME1)/BCHGTIMINC
            do I = 1, NP
               DPTMP = bTRatio*(DP2(I) - DP1(I))
               DPTMP2 = DP1(I) + DPTMP
               DPTMP = DP(I) - DPTMP2
               DP(I) = DPTMP2
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

         ! Finalize bathymetry at end of change interval
         if (abs(TimeLoc - BTIME_END) <= eps) then
            write (scratchMessage, '(A42,1X,E15.8,1X,A5)') &
               'BATHYMETRY VALUES ARE NOW FIXED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
            do I = 1, NP
               DPTMP = DP2(I) - DP(I)
               DP(I) = DP2(I)
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

      else if (abs(NDDT) == 2) then
         ! Mode 2: Read bathymetry for sparse nodes
         if (TimeLoc > BTIME2) then
            BTIME1 = BTIME2
            BTIME2 = BTIME2 + BTIMINC
            BTIME_END = BTIME1 + BCHGTIMINC
            dp1(:) = dp2(:)
            ! Read new record for selected nodes only
            call NDDT2GET(141, DP2(:), -99999.d0)
            ! Ensure depths >= H0 if no wetting/drying
            if (NOLIFA == 0 .or. NOLIFA == 1) then
               where (DP2(1:NP) < H0) DP2(1:NP) = H0
            end if
            write (scratchMessage, '(A36,1X,E15.8,1X,A5)') &
               'BATHYMETRY RECORDS UPDATED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
         end if

         ! Interpolate during bathymetry change interval
         if (TimeLoc < BTIME_END) then
            bTRatio = (TimeLoc - BTIME1)/BCHGTIMINC
            do I = 1, NP
               DPTMP = bTRatio*(DP2(I) - DP1(I))
               DPTMP2 = DP1(I) + DPTMP
               DPTMP = DP(I) - DPTMP2
               DP(I) = DPTMP2
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

         ! Finalize bathymetry at end of change interval
         if (abs(TimeLoc - BTIME_END) <= eps) then
            write (scratchMessage, '(A42,1X,E15.8,1X,A5)') &
               'Bathymetry values are now fixed at time =', TimeLoc, '(sec).'
            call allMessage(INFO, scratchMessage)
            do I = 1, NP
               DPTMP = DP2(I) - DP(I)
               DP(I) = DP2(I)
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if
      end if

   end subroutine updateTimeDependentBathymetry
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute bottom friction coefficients
!>
!> Computes 2D or 3D bottom friction coefficients (TK) based on the
!> selected friction formulation. For baroclinic runs, also computes
!> internal wave drag and momentum dispersion coefficients (TKM).
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine computeBottomFriction(IT)
      use global, only: C2DDI, C3D, CBaroclinic, UU1, VV1, ETA2, &
                        IFNLFA, TK, TKM
      use global_3dvs, only: Q, SIGMA, NFEN, Z0B, KP
      use mesh, only: NP, DP
      use ADC_CONSTANTS, only: G
      use NodalAttributes, only: Apply2DBottomFriction, &
                                 Apply2DInternalWaveDrag, Apply3DBottomFriction, &
                                 Apply2DMomentumDisp
      implicit none

      integer, intent(in) :: IT

      ! 2DDI: Set up the 2D friction coefficient
      if (C2DDI) then
         call Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, IFNLFA, &
                                    NP, TK)
         ! Get wave drag coefficients into TKM (includes TK)
         if (CBaroclinic) then
            call Apply2DInternalWaveDrag(NP, TK, TKM, UU1, VV1, DP, IT)
            ! Vertically averaged momentum dispersion (power law)
            call Apply2DMomentumDisp(UU1, VV1, DP, ETA2, IFNLFA, &
                                     NP, TKM)
         else
            call Apply2DInternalWaveDrag(NP, TK, TKM)
         end if
      elseif (C3D) then
         ! 3D: Set up the 3D friction coefficient
         call Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, IFNLFA, &
                                    NP, TK, NFEN, Z0B, KP)
      end if

   end subroutine computeBottomFriction
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update ice concentration fields
!>
!> Reads ice concentration data from input file when needed and linearly
!> interpolates ice values (CICEOUT) for the current timestep.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateIceFields(TimeLoc)
      use global, only: NCICE, CICE_TIME1, CICE_TIME2, CICE_TIMINC, &
                        CICE1, CICE2, CICEOUT
      use mesh, only: NP
      use owi_ice, only: NCICE1_GET
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: CICE_TRatio
      integer :: I

      ! Read new ice data if needed
      if (NCICE == 12) then
         if (TimeLoc > CICE_TIME2) then
            CICE_TIME1 = CICE_TIME2
            CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
            CICE1(:) = CICE2(:)
            call NCICE1_GET(CICE2, NP)
         end if
      end if

      ! Interpolate ice concentration for current time
      if (NCICE /= 0 .and. NCICE /= 14 .and. NCICE /= 17) then
         CICE_TRatio = (TimeLoc - CICE_TIME1)/CICE_TIMINC
         do I = 1, NP
            CICEOUT(I) = CICE1(I) + CICE_TRatio*(CICE2(I) - CICE1(I))
         end do
      end if

   end subroutine updateIceFields
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update dynamic water level corrections
!>
!> Retrieves dynamic water level corrections for the current timestep
!> and computes the delta for mass adjustment if enabled.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateWaterLevelCorrections(TimeLoc)
      use global, only: usingDynamicWaterLevelCorrection, &
                        dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, &
                        dynamicWaterLevelCorrectionDelta
      use wind, only: getDynamicWaterLevelCorrections, &
                      dynamicWaterLevelCorrectionMassAdjust
      implicit none

      real(8), intent(in) :: TimeLoc

      if (usingDynamicWaterLevelCorrection) then
         call getDynamicWaterLevelCorrections( &
            dynamicWaterLevelCorrection2, TimeLoc)
         if (dynamicWaterLevelCorrectionMassAdjust) then
            dynamicWaterLevelCorrectionDelta = &
               dynamicWaterLevelCorrection2 - &
               dynamicWaterLevelCorrection1
         end if
      end if

   end subroutine updateWaterLevelCorrections
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update wave radiation stress forcing
!>
!> Reads wave radiation stress data from file or coupled wave model and
!> adds to wind stress arrays (WSX2, WSY2). Supports multiple modes:
!> - NRS=1: Read from fort.23
!> - NRS=2: Read from RS2 format
!> - NRS=3: SWAN coupling
!> - NRS=4: Wave time extrapolation
!> - NRS=5: Linear extrapolation
!>
!> Optionally limits wave stress gradients to prevent numerical instability.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateWaveRadiationStress(TimeLoc)
      use global, only: NRS, RSTIME1, RSTIME2, RSTIMINC, &
                        RSNX1, RSNY1, RSNX2, RSNY2, WSX2, WSY2, &
                        RSNXOUT, RSNYOUT, ENDWAVE, Limit_WaveStressGrad, &
                        RampWRad
      use ADC_CONSTANTS, only: WaveStressGrad_Cap
      use mesh, only: NP
      use wind, only: RSGET
      use rs2, only: RS2GET
#ifdef CSWAN
      use global, only: NWS, wvnxout, wvnyout, RampMete, &
                        SWAN_RSNX1, SWAN_RSNX2, SWAN_RSNY1, SWAN_RSNY2
      use ADC_CONSTANTS, only: waveWindMultiplier
      use owiwind, only: windMultiplier
      use couple2swan, only: coupwind, swan_wx2, swan_wy2, &
                             InterpoWeight, ComputeWaveDrivenForces
#endif
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: RStRatio, RSX, RSY, RSXYMAG, RSSCALE, RSXLIM, RSYLIM
      integer :: I

#ifdef CSWAN
      ! Coupling winds to SWAN
      if (coupwind) then
         select case (abs(nws))
         case (12, 16)
            do i = 1, np
               swan_wx2(i, 2) = wvnxout(i)*waveWindMultiplier &
                                /(rampMete*windMultiplier)
               swan_wy2(i, 2) = wvnyout(i)*waveWindMultiplier &
                                /(rampMete*windMultiplier)
            end do
         case default
            do i = 1, np
               swan_wx2(i, 2) = wvnxout(i)*waveWindMultiplier/rampMete
               swan_wy2(i, 2) = wvnyout(i)*waveWindMultiplier/rampMete
            end do
         end select
      end if
#endif

      if (NRS /= 0) then
         if ((NRS == 1) .or. (nrs == 2) .or. (nrs == 3) .or. (nrs == 5)) then
            if (TimeLoc > RSTIME2) then
               RSTIME1 = RSTIME2
               RSTIME2 = RSTIME2 + RSTIMINC
               RSNX1 = RSNX2
               RSNY1 = RSNY2
               if (NRS == 1) then
                  call RSGET(RSNX2, RSNY2)
               elseif (NRS == 2) then
                  call RS2GET(RSNX2, RSNY2, NP)
               elseif (NRS == 3) then
#ifdef CSWAN
                  InterpoWeight = 1.0d0
                  call ComputeWaveDrivenForces()
                  RSNX2 = 2.0d0*SWAN_RSNX2 - SWAN_RSNX1
                  RSNY2 = 2.0d0*SWAN_RSNY2 - SWAN_RSNY1
#endif
               elseif (NRS == 5) then
                  do I = 1, NP
                     RSX = RSNX1(I)
                     RSY = RSNY1(I)
                     RSNX1(I) = RSNX2(I)
                     RSNY1(I) = RSNY2(I)
                     RSNX2(I) = RSNX2(I) + (RSNX2(I) - RSX)
                     RSNY2(I) = RSNY2(I) + (RSNY2(I) - RSY)
                  end do
               end if
            end if

            RStRatio = (TimeLoc - RSTIME1)/RSTIMINC

            if (.not. Limit_WaveStressGrad) then
               do I = 1, NP
                  RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I) - RSNX1(I)))
                  RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I) - RSNY1(I)))
                  WSX2(I) = WSX2(I) + RSX
                  WSY2(I) = WSY2(I) + RSY
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            else
               do I = 1, NP
                  RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I) - RSNX1(I)))
                  RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I) - RSNY1(I)))
                  RSXYMAG = sqrt(RSX*RSX + RSY*RSY)
                  if ((RSXYMAG > WaveStressGrad_Cap) .and. &
                      (RSXYMAG > 0.0d0)) then
                     RSSCALE = WaveStressGrad_Cap/RSXYMAG
                     RSXLIM = RSX*RSSCALE
                     RSYLIM = RSY*RSSCALE
                  else
                     RSXLIM = RSX
                     RSYLIM = RSY
                  end if
                  WSX2(I) = WSX2(I) + RSXLIM
                  WSY2(I) = WSY2(I) + RSYLIM
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            end if

         end if

         if (NRS == 4) then
            if (TimeLoc >= RSTIME2) then
               RSTIME1 = RSTIME2
               RSTIME2 = RSTIME2 + RSTIMINC
            end if
            if (TimeLoc > ENDWAVE + RSTIMINC) then
               RSNX2(:) = 0.0d0; RSNY2(:) = 0.0d0
            end if

            if (.not. Limit_WaveStressGrad) then
               do I = 1, NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)
                  WSX2(I) = WSX2(I) + RSX
                  WSY2(I) = WSY2(I) + RSY
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            else
               do I = 1, NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)
                  RSXYMAG = sqrt(RSX*RSX + RSY*RSY)
                  if ((RSXYMAG > WaveStressGrad_Cap) .and. &
                      (RSXYMAG > 0.0d0)) then
                     RSSCALE = WaveStressGrad_Cap/RSXYMAG
                     RSXLIM = RSX*RSSCALE
                     RSYLIM = RSY*RSSCALE
                  else
                     RSXLIM = RSX
                     RSYLIM = RSY
                  end if
                  WSX2(I) = WSX2(I) + RSXLIM
                  WSY2(I) = WSY2(I) + RSYLIM
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            end if
         end if
      end if

   end subroutine updateWaveRadiationStress
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply wetting and drying algorithm
!>
!> Applies the wetting and drying algorithm to handle tidal flats and
!> inundation. Under subgrid level0 the subgrid lookups are first updated
!> to the new water surface elevation.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine applyWettingAndDrying(IT)
      use global, only: ETA2
      use mod_subgrid, only: subgrid_level0, subgrid_update
      use wetdry, only: computeWettingAndDrying
      implicit none

      integer, intent(in) :: IT

      if (subgrid_level0) call subgrid_update(ETA2)
      call computeWettingAndDrying(IT)

   end subroutine applyWettingAndDrying
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update total water depth
!>
!> Updates the total water depth (H2) from bathymetry and water surface
!> elevation, or sets it to the grid-averaged depth H_G under subgrid
!> level0. Called after wetting/drying since effective bathymetry
!> could change with subgrid barriers.
!-----------------------------------------------------------------------
   subroutine updateFlowDepth()
      use global, only: H2, ETA2, IFNLFA
      use mesh, only: DP
      use mod_subgrid, only: subgrid_level0, gridDepthVertETA2
      implicit none

      if (subgrid_level0) then
         H2 = gridDepthVertETA2
      else
         H2 = DP + dble(IFNLFA)*ETA2
      end if

   end subroutine updateFlowDepth
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Collect output statistics and write subdomain output
!>
!> Collects min/max elevation and velocity data, inundation statistics
!> (if enabled), updates harmonic analysis, and writes subdomain output
!> files (fort.065/066/067) if subdomain modeling is enabled.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TIMEH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine collectOutputStatistics(IT, TimeLoc, TIMEH)
      use global, only: inundationOutput
      use write_output, only: collectMinMaxData, collectInundationData
      use subdomain, only: subdomainOn, NOutGS, &
                           writeFort065, writeFort066, writeFort067
      use HARM, only: updateHarmonicAnalysis
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: TIMEH

      call collectMinMaxData(TimeLoc)

      if (inundationOutput) call collectInundationData(TimeLoc, IT)

      if (subdomainOn) then
         if (NOutGS == 1) call writeFort065(IT)
         if (NOutGS == 2) call writeFort066(IT)
         if (NOutGS == 2) call writeFort067(IT)
      end if

      !     If harmonic analysis was requested and the current time
      !     is within the harmonic analysis period, update the left hand side
      !     of the harmonic analysis matrix. Also update the load vectors for
      !     each type of analysis. If timeseries reconstruction was specified,
      !     also update the timeseries.
      call updateHarmonicAnalysis(IT, TIMEH)

   end subroutine collectOutputStatistics
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Write hot start file if needed
!>
!> Writes hot start (restart) files at specified intervals for simulation
!> checkpointing. Supports both direct file writing and the HS writer for
!> parallel runs. Uses 8-byte record length for double precision data.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine writeHotstartIfNeeded(IT, TimeLoc)
      use sizes, only: MNWPROH
      use global, only: NHSINC, NHSTAR, IHOT, C3D, NHOUTONCE
      use write_output, only: writeHotstart
#ifdef CMPI
      use hswriter, only: writeHotstart_through_hswriter
#endif
#ifdef CSWAN
      use couple2swan, only: WriteSwanHotStart
#endif
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      integer :: ITEST

      ITEST = (IT/NHSINC)*NHSINC

      if ((abs(NHSTAR) > 0 .and. ITEST == IT) .or. (-IHOT == IT)) then
         if (MNWPROH > 0) then
            if (.not. C3D) then
#ifdef CMPI
               call writeHotstart_through_hswriter(TimeLoc, IT)
#endif
            else
               write (6, *) 'HS writer does not support C3D'
            end if
         else
            if (NHOUTONCE) then
               NHSINC = huge(NHSINC)
            end if
            call writeHotstart(TimeLoc, IT)
         end if
#ifdef CSWAN
         WriteSwanHotStart = .true.
#endif
      end if

   end subroutine writeHotstartIfNeeded
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Monitor and report simulation progress
!>
!> Monitors solution progress by finding and reporting maximum elevation
!> and velocity magnitudes. Checks values against warning and error
!> thresholds (WarnElev, ErrorElev, WarnVel). Detects NaN values in
!> velocity fields. Terminates simulation if error conditions are met.
!>
!> @param[in] IT        Current timestep number
!> @param[in] TimeLoc   Current simulation time (seconds)
!> @param[in] ITIME_BGN Starting timestep number (for progress percentage)
!-----------------------------------------------------------------------
   subroutine progressScreen(IT, TimeLoc)
      use global, only: ETA2, NODECODE, UU2, VV2, NT, ITHS, &
                        WarnElev, ErrorElev, Flag_ElevError, WarnElevDump
      use mod_logging, only: screenUnit, nscreen
      use mesh, only: NP
      use gwce, only: numitr
      use write_output, only: WriteWarnElev
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
#if defined(CMPI) || defined(VELCHECK)
      use sizes, only: MYPROC
      use global, only: nodes_lg
#endif
#ifdef VELCHECK
      use global, only: earlyterminate, IT_foundNaN, WarnVel
      use mesh, only: DP, NNeigh, NeiTab
#endif
#ifdef CMPI
      use messenger, only: WarnElevSum, EarlyTermSum
#endif
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc

      integer :: I
      integer :: KEMAX, KVMAX
      integer :: ITEST
#ifdef CMPI
      integer :: WarnElevExceeded, ErrorElevExceeded
#endif
#ifdef VELCHECK
      integer :: J, JJ
      integer :: WarnVelExceeded
      logical :: foundNaN
#endif
      real(8) :: ELMAX, VELMAX, velabs
      real(8) :: runperccomplete

      if (NSCREEN /= 0) then
         ELMAX = 0.0d0
         VELMAX = 0.0d0
         KEMAX = 0
         KVMAX = 0
         do I = 1, NP
            if ((NODECODE(I) == 1) .and. (abs(ETA2(I)) > ELMAX)) then
               ELMAX = abs(ETA2(I))
               KEMAX = I
            end if
            VELABS = UU2(I)*UU2(I) + VV2(I)*VV2(I)
            if (VELABS > VELMAX) then
               VELMAX = VELABS
               KVMAX = I
            end if
         end do
         VELMAX = VELMAX**0.5d0
         ITEST = (IT/NSCREEN)*NSCREEN

         runperccomplete = real(IT - ITHS, 8)/real(NT - ITHS, 8)*100d0
#ifdef CMPI
         if (MYPROC == 0 .and. ELMAX < WarnElev .and. ITEST == IT) then
            if (KEMAX > 0) then
               write (ScreenUnit, 1991) &
                  IT, runperccomplete, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
                  VELMAX, KVMAX, MYPROC
            else
               write (ScreenUnit, 1991) &
                  IT, runperccomplete, NUMITR, TimeLoc, 0., KEMAX, VELMAX, &
                  KVMAX, MYPROC
            end if
1991        format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', &
                   5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I8, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I8, &
                   2x, 'ON MYPROC = ', I4)
         end if
         WarnElevExceeded = 0
         if (ELMAX > WarnElev) then
            write (ScreenUnit, 1993) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 1993) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
1993        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, &
                   3x, '** WARNING: Elevation.gt.WarnElev **')
            if (WarnElevDump) WarnElevExceeded = 1
         end if

         if (WarnElevDump) then
            call WarnElevSum(WarnElevExceeded)
            if (WarnElevExceeded /= 0) then
               call WriteWarnElev(TimeLoc, IT)
            end if
         end if
         ErrorElevExceeded = 0
         if (ELMAX > ErrorElev) then
            write (ScreenUnit, 1995) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 1995) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
1995        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, /, &
                   2x, '** ERROR: Elevation.gt.ErrorElev,', &
                   ' ADCIRC stopping. **')
            ErrorElevExceeded = 1
         end if
         call WarnElevSum(ErrorElevExceeded)
         if (ErrorElevExceeded /= 0) then
            Flag_ElevError = .true.
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Elevation exceeded ErrorElev threshold')
         end if

#ifdef VELCHECK
         WarnVelExceeded = 0
         if (VELMAX > WarnVel) then
            write (ScreenUnit, 3995) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 3995) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
3995        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, /, &
                   2x, '** WARNING: Velocity.gt.WarnVel **')
            WarnVelExceeded = 1
         end if
#endif

#else
         if (ELMAX < WarnElev .and. ITEST == IT) then
            if (KEMAX > 0) then
               write (ScreenUnit, 1992) &
                  IT, runperccomplete, NUMITR, TimeLoc, ETA2(KEMAX), &
                  KEMAX, VELMAX, KVMAX
            else
               write (ScreenUnit, 1992) &
                  IT, runperccomplete, NUMITR, TimeLoc, 0., KEMAX, &
                  VELMAX, KVMAX
            end if
1992        format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', &
                   5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9)
         end if
         if (ELMAX > WarnElev) then
            write (ScreenUnit, 1994) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 1994) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
1994        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, '** WARNING: Elevation.gt.WarnElev **')
            if (WarnElevDump) call WriteWarnElev(TimeLoc, IT)
         end if
         if (ELMAX > ErrorElev) then
            write (6, *) 'elmax, errorelev', elmax, errorelev
            Flag_ElevError = .true.
            write (ScreenUnit, 1996) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 1996) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
1996        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, /, &
                   2x, '** ERROR: Elevation.gt.ErrorElev, ', &
                   'ADCIRC stopping. **')
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Elevation exceeded ErrorElev threshold')
         end if
#ifdef VELCHECK
         if (VELMAX > WarnVel) then
            write (ScreenUnit, 3996) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 3996) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
3996        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, /, &
                   2x, '** WARNING: Velocity.gt.WarnVel **')
         end if
#endif

#endif

      end if

#ifdef VELCHECK
      foundNaN = .false.
      do I = 1, NP
         if (ISNAN(UU2(I)) .or. ISNAN(VV2(I))) then
            foundNaN = .true.
            if (NSCREEN /= 0) then
               write (ScreenUnit, 2993) &
                  IT, NUMITR, TimeLoc, NODES_LG(I), &
                  ETA2(I), DP(I), MYPROC
               write (16, 2993) IT, NUMITR, TimeLoc, &
                  NODES_LG(I), ETA2(I), DP(I), MYPROC
2993           format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                      5x, 'TIME = ', E15.8, &
                      2x, 'SPEED IS NaN AT NODE I =', I9, &
                      /, 2x, 'ETA2(I) = ', 1pe12.4e3, &
                      2x, 'DP(I) = ', 1pe12.4e3, &
                      2x, 'ON MYPROC = ', I4)
               do J = 2, NNeigh(I)
                  JJ = NeiTab(I, J)
                  write (ScreenUnit, 2995) NODES_LG(JJ), ETA2(JJ), DP(JJ), &
                     UU2(JJ), VV2(JJ), MYPROC
                  write (16, 2995) NODES_LG(JJ), ETA2(JJ), DP(JJ), &
                     UU2(JJ), VV2(JJ), MYPROC
2995              format(1x, 'AT NEIGHBOR NODE I =', I9, &
                         /, 2x, 'ETA2(I) = ', 1pe12.4e3, &
                         2x, 'DP(I) = ', 1pe12.4e3, &
                         2x, 'UU2(I) = ', 1pe12.4e3, &
                         2x, 'VV2(I) = ', 1pe12.4e3, &
                         2x, 'ON MYPROC = ', I4)
               end do
            end if
         end if
      end do
      if (foundNaN) then
         if (NSCREEN /= 0) then
            write (ScreenUnit, 2994) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 2994) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
2994        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, &
                   3x, '** WARNING: Velocity is NaN **')
         end if
      end if
      if ((foundNaN) .and. (IT_foundNaN == 0)) then
         IT_foundNaN = IT
      elseif ((IT == IT_foundNaN + 20) .and. (IT_foundNaN /= 0)) then
         earlyterminate = 1
      end if
#ifdef CMPI
      call EarlyTermSum(earlyterminate)
#endif
#endif

   end subroutine progressScreen
   !----------------------------------------------------------------------

end module mod_timestep
