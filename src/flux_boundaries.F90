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
!-----------------------------------------------------------------------
!  MODULE FLUX_BOUNDARIES
!-----------------------------------------------------------------------
!> @brief Flux-type boundary conditions applied each timestep
!>
!> This module computes the normal flux (QN1/QN2) and elevation (EN2)
!> boundary values used by the GWCE and momentum equations for:
!> - Specified normal flow boundaries (periodic or fort.20 time series)
!> - Radiation and zero normal velocity gradient boundaries
!> - External and internal barrier (weir) boundaries, including vertical
!>   element walls and cross-barrier pipes
!> - Subdomain modeling boundary data
!-----------------------------------------------------------------------
module mod_flux_boundaries

   implicit none

   private

   public :: applyBoundaryConditions

   ! Whether the river radiation b.c. elevations still need to be filled
   ! from EtaDisc after flux settling
   logical :: EtaDisc_Fill = .true.

   real(8), parameter :: eps = epsilon(1.0d0)

contains

!-----------------------------------------------------------------------
!> @brief Compute specified normal flow boundary conditions
!>
!> Computes normal flux (QN2) and elevation (EN2) boundary conditions
!> for specified flux boundaries. Supports periodic forcing with tidal
!> constituents and time-varying flux input from fort.20.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine computeSpecifiedNormalFlowBC(IT, TimeLoc, TimeH)
      use global, only: NFFR, FPER, FAMIG, FFACE, FFF, RampExtFlux, &
                        QNAM, QNPH, ENAM, ENPH, QN2, EN2, &
                        QTIME1, QTIME2, FTIMINC, QNIN1, QNIN2, ENIN1, ENIN2, &
                        FluxSettlingIT, EtaDisc, ElevDisc, Eta2
      use mesh, only: NP, ICS
      use boundaries, only: NVEL, LBCODEI, NBV
      use mod_normal_flow_boundary, only: rotate_normal_flux
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc, TimeH
      integer :: I, J, NCYC, NNBB
      real(8) :: ARGJ, RFF, QTRATIO, QNIN_TMP

      if (NFFR > 0) then
         do J = 1, NFFR
            if (abs(FPER(J)) <= eps) then
               NCYC = 0
            else
               NCYC = int(TimeH/FPER(J))
            end if
            ARGJ = FAMIG(J)*(TimeH - dble(NCYC)*FPER(J)) + FFACE(J)
            RFF = FFF(J)*RampExtFlux
            do I = 1, NVEL
               select case (LBCODEI(I))
               case (2, 12, 22, 52)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF*cos(ARGJ - QNPH(J, I)))
               case (32)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF*cos(ARGJ - QNPH(J, I)))
                  EN2(I) = EN2(I) + ENAM(J, I)*RFF*cos(ARGJ - ENPH(J, I))
               end select
            end do
         end do
      elseif ((NFFR == 0) .or. (NFFR == -1)) then
         if (TimeLoc > QTIME2) then
            QTIME1 = QTIME2
            QTIME2 = QTIME2 + FTIMINC
            do J = 1, NVEL
               if ((LBCODEI(J) == 2) .or. (LBCODEI(J) == 12) &
                   .or. (LBCODEI(J) == 22)) then
                  QNIN1(J) = QNIN2(J)
                  read (20, *) QNIN_TMP
                  QNIN2(J) = rotate_normal_flux(ICS, J, QNIN_TMP)
               elseif (LBCODEI(J) == 32) then
                  QNIN1(J) = QNIN2(J)
                  ENIN1(J) = ENIN2(J)
                  read (20, *) QNIN_TMP, ENIN2(J)
                  QNIN2(J) = rotate_normal_flux(ICS, J, QNIN_TMP)
               end if
            end do
         end if

         QTRATIO = (TimeLoc - QTIME1)/FTIMINC
         QN2 = RampExtFlux*(QNIN1 + QTRATIO*(QNIN2 - QNIN1))
         EN2 = RampExtFlux*(ENIN1 + QTRATIO*(ENIN2 - ENIN1))
      end if

      ! Collect elevation information for river radiation b.c.
      if (IT == FluxSettlingIT) then
         EtaDisc_Fill = .false.
         do I = 1, NP
            EtaDisc(I) = Eta2(I)
         end do
         do I = 1, NVEL
            if (LBCODEI(I) == 52) then
               NNBB = NBV(I)
               ElevDisc(I) = Eta2(NNBB)
            end if
         end do
      else if (EtaDisc_Fill .and. IT > FluxSettlingIT) then
         EtaDisc_Fill = .false.
         do I = 1, NVEL
            if (LBCODEI(I) == 52) then
               NNBB = NBV(I)
               ElevDisc(I) = EtaDisc(NNBB)
            end if
         end do
      end if

   end subroutine computeSpecifiedNormalFlowBC
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute radiation boundary condition discharge
!>
!> Computes discharge (QN1) at radiation boundary nodes (LBCODE=30)
!> based on previous timestep velocities and current depth.
!-----------------------------------------------------------------------
   subroutine computeRadiationDischargeBC()
      use global, only: UU1, VV1, H2, QN1
      use boundaries, only: NVEL, LBCODEI, NBV, CSII, SIII
      implicit none

      integer :: J, NNBB
      real(8) :: UN1

      do J = 1, NVEL
         if (LBCODEI(J) == 30) then
            NNBB = NBV(J)
            UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
            QN1(J) = H2(NNBB)*UN1
         end if
      end do

   end subroutine computeRadiationDischargeBC
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute discharge for zero normal velocity gradient BC
!>
!> Computes discharge (QN1) at zero normal velocity gradient boundary
!> nodes (LBCODE=40,41) based on current velocities and depth.
!-----------------------------------------------------------------------
   subroutine computeZeroNormalVelocityBCDischarge()
      use global, only: UU1, VV1, H2, QN1
      use boundaries, only: NVEL, LBCODEI, NBV, CSII, SIII
      implicit none

      integer :: J, NNBB
      real(8) :: UN1

      do J = 1, NVEL
         if ((LBCODEI(J) == 40) .or. (LBCODEI(J) == 41)) then
            NNBB = NBV(J)
            UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
            QN1(J) = H2(NNBB)*UN1
         end if
      end do

   end subroutine computeZeroNormalVelocityBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute external barrier boundary flux
!>
!> Computes supercritical outward normal flow (QN2) over external barrier
!> boundary nodes (LBCODE=3,13,23) using weir formulas.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeExternalBarrierBCDischarge(TimeLoc)
      use global, only: QN2
      use boundaries, only: NVEL, LBCODEI
      use weir_flux, only: COMPUTE_EXTERNAL_BOUNDARY_FLUX
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I

      do I = 1, NVEL
         select case (LBCODEI(I))
         case (3, 13, 23)
            call COMPUTE_EXTERNAL_BOUNDARY_FLUX(I, TimeLoc, QN2(I))
         end select
      end do

   end subroutine computeExternalBarrierBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Set submergence flags for vertical element wall boundaries
!>
!> Sets the ISSUBMERGED flag for vertical element wall (VEW) boundaries
!> (LBCODE=64) based on current water levels.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine setVewBCSubmergedFlag(TimeLoc)
      use boundaries, only: NBOU, LBCODEI, NVELL
      use WEIR_FLUX, only: SET_SUBMERGED64_AT
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I, J, K

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               call SET_SUBMERGED64_AT(I, J, K, TimeLoc)
            end do
            I = I + NVELL(K)
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine setVewBCSubmergedFlag
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute internal barrier boundary flux
!>
!> Computes inward/outward normal flow (QN2) over internal barrier
!> boundary nodes (LBCODE=4,24,5,25,64) using weir formulas. Handles
!> both permeable and impermeable barrier types.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeInternalBarrierBCDischarge(TimeLoc)
      use global, only: QN2, NIBNODECODE
      use boundaries, only: NBOU, LBCODEI, NVELL
      use WEIR_FLUX, only: COMPUTE_INTERNAL_BOUNDARY_FLUX, &
                           COMPUTE_INTERNAL_BOUNDARY64_FLUX
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I, J, K
      logical :: ISFRONT

      NIBNODECODE(:) = 0
      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (4, 24, 5, 25)
            do J = 1, NVELL(K)*2
               I = I + 1
               call COMPUTE_INTERNAL_BOUNDARY_FLUX(I, J, K, &
                                                   TimeLoc, QN2(I))
            end do
         case (64)
            do J = 1, NVELL(K)*2
               I = I + 1
               if (J <= NVELL(K)) then
                  ISFRONT = .true.
               else
                  ISFRONT = .false.
               end if
               call COMPUTE_INTERNAL_BOUNDARY64_FLUX(I, J, K, &
                                                     TimeLoc, QN2(I))
            end do
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine computeInternalBarrierBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute cross-barrier pipe flux
!>
!> Computes inward/outward normal flow (QN2) through cross-barrier pipes
!> for internal barrier boundary nodes (LBCODE=5,25).
!-----------------------------------------------------------------------
   subroutine computeInternalBarrierPipeFlux()
      use global, only: QN2
      use boundaries, only: NVEL, LBCODEI
      use WEIR_FLUX, only: COMPUTE_CROSS_BARRIER_PIPE_FLUX
      implicit none

      integer :: I
      real(8) :: PIPE_FLUX

      do I = 1, NVEL
         if ((LBCODEI(I) == 5) .or. (LBCODEI(I) == 25)) then
            PIPE_FLUX = 0d0
            call COMPUTE_CROSS_BARRIER_PIPE_FLUX(I, PIPE_FLUX)
            QN2(I) = QN2(I) + PIPE_FLUX
         end if
      end do

   end subroutine computeInternalBarrierPipeFlux
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply subdomain boundary condition data
!>
!> Reads and applies subdomain boundary condition data from fort.019,
!> fort.020, and fort.021 files for subdomain modeling.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine applySubdomainBCData(IT)
      use subdomain, only: enforceBN, &
                           readFort019, readFort020, readFort021
      implicit none

      integer, intent(in) :: IT

      if (enforceBN == 1) call readFort019(IT)
      if (enforceBN == 2) call readFort020(IT)
      if (enforceBN == 2) call readFort021(IT)

   end subroutine applySubdomainBCData
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply all boundary conditions
!>
!> Master routine that applies all boundary conditions including normal
!> flow, radiation, barrier, VEW, and subdomain boundaries.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine applyBoundaryConditions(IT, TimeLoc, TimeH)
      use boundaries, only: NFLUXIB, NFLUXF, NFLUXIBP, NFLUXGBC, &
                            NFLUXB, NFLUXIB64, NFLUXRBC
      use subdomain, only: subdomainOn
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc, TimeH

      if (NFLUXF == 1) call computeSpecifiedNormalFlowBC(IT, TimeLoc, TimeH)
      if (NFLUXRBC == 1) call computeRadiationDischargeBC()
      if (NFLUXGBC == 1) call computeZeroNormalVelocityBCDischarge()
      if (NFLUXB == 1) call computeExternalBarrierBCDischarge(TimeLoc)
      if (NFLUXIB64 == 1) call setVewBCSubmergedFlag(TimeLoc)
      if (NFLUXIB == 1) call computeInternalBarrierBCDischarge(TimeLoc)
      if (NFLUXIBP == 1) call computeInternalBarrierPipeFlux()

      if (subdomainOn) call applySubdomainBCData(IT)

   end subroutine applyBoundaryConditions
   !-----------------------------------------------------------------------
end module mod_flux_boundaries
