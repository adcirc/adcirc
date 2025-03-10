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
!
module mod_timestep_bc_forcing

   implicit none

   public :: apply_boundary_conditions

   private

contains

   subroutine apply_boundary_conditions(IT, FluxSettlingIT, NFFR, timeh, TimeLoc, FTIMINC, RampExtFlux, &
                                        NODECODE, NIBNODECODE, Eta2, DP, H2, UU1, VV1, QTIME1, QTIME2, ENIN1, &
                                        ENIN2, QNIN1, QNIN2, EtaDisc, QN1, QN2, EN2, ElevDisc, &
                                        TVW, EtaDisc_Fill, ISSUBMERGED64, ISSUBMERGED64P)

      use boundaries, only: NBOU, NFLUXF, NFLUXIB, LBCODEI, NBV, NVEL, NVELL, CSII, SIII, IBCONN
      use global, only: qnph, qnam, fper, famig, fface, fff, enph, enam, h0, ifnlfa
      use mesh, only: LBArray_Pointer, NeiTab, NNeigh, X, Y
      use global, only: qnph, qnam, fper, famig, fface, fff, enph, enam, h0, ifnlfa
      use subdomain, only: subdomainOn, enforceBN

      implicit none

      integer, intent(in) :: IT
      integer, intent(in) :: FluxSettlingIT
      integer, intent(in) :: NFFR
      real(8), intent(in) :: timeh
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: FTIMINC
      real(8), intent(in) :: RampExtFlux
      real(8), intent(in) :: Eta2(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(in) :: DP(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: QTIME1
      real(8), intent(inout) :: QTIME2
      real(8), intent(inout) :: ENIN1(:)
      real(8), intent(inout) :: ENIN2(:)
      real(8), intent(inout) :: QNIN1(:)
      real(8), intent(inout) :: QNIN2(:)
      real(8), intent(inout) :: EtaDisc(:)
      real(8), intent(inout) :: QN1(:)
      real(8), intent(inout) :: QN2(:)
      real(8), intent(inout) :: EN2(:)
      real(8), intent(inout) :: ElevDisc(:)
      real(8), intent(inout) :: TVW(:)
      logical, intent(inout) :: EtaDisc_Fill
      integer, intent(inout) :: ISSUBMERGED64(NVEL)
      integer, intent(inout) :: ISSUBMERGED64P(NVEL)

      call apply_specified_normal_flow_boundary_conditions(IT, NFLUXF, &
                                                           NFFR, NVEL, FluxSettlingIT, &
                                                           LBCODEI, NBV, timeh, TimeLoc, FTIMINC, &
                                                           RampExtFlux, QNPH, QNAM, FPER, FAMIG, &
                                                           FFACE, FFF, QNIN1, QNIN2, ENPH, ENAM, &
                                                           QTIME1, QTIME2, ENIN1, ENIN2, Eta2, &
                                                           EtaDisc, QN2, EN2, ElevDisc, &
                                                           EtaDisc_Fill)
      call apply_computed_normal_flow_boundary_conditions(NFLUXIB, NBOU, NVELL, LBCODEI, NBV, &
                                                          TIMELOC, TVW, UU1, VV1, CSII, SIII, H2, X, Y, &
                                                          DP, ETA2, H0, IBCONN, IFNLFA, LBArray_Pointer, &
                                                          NeiTab, NNeigh, NODECODE, NIBNODECODE, QN1, QN2, &
                                                          ISSUBMERGED64, ISSUBMERGED64P)

      if (subdomainOn) call enforce_subdomain_boundaries(IT, enforceBN)

   end subroutine apply_boundary_conditions

   !*******************************************************************************
   !> Apply flow boundary conditions that are computed internally
   !>
   !> @param[in] NFLUXIB If there are internal weir flow boundaries
   !> @param[in] NBOU Number of boundaries
   !> @param[in] NVELL Number of velocity components
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] NBV Array of boundary vertices
   !> @param[in] TIMELOC Current time
   !> @param[in] TVW Array of time varying weir elevations
   !> @param[in] UU1 Array of x-velocity components
   !> @param[in] VV1 Array of y-velocity components
   !> @param[in] CSII Array of cosine of the angle
   !> @param[in] SIII Array of sine of the angle
   !> @param[in] H2 Array of depths
   !> @param[in,out] NIBNODECODE Array of internal boundary node codes
   !> @param[in,out] QN1 Array of normal flow boundary fluxes
   !> @param[in,out] QN2 Array of normal flow boundary fluxes
   !*******************************************************************************
   subroutine apply_computed_normal_flow_boundary_conditions(NFLUXIB, NBOU, NVELL, LBCODEI, NBV, &
                                                             TIMELOC, TVW, UU1, VV1, CSII, SIII, H2, X, Y, &
                                                             DP, ETA2, H0, IBCONN, IFNLFA, LBArray_Pointer, &
                                                             NeiTab, NNeigh, NODECODE, NIBNODECODE, QN1, QN2, &
                                                             ISSUBMERGED64, ISSUBMERGED64P)

      implicit none

      integer, intent(in) :: NFLUXIB
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: LBArray_Pointer(:)
      integer, intent(in) :: NeiTab(:, :)
      integer, intent(in) :: NNeigh(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: IFNLFA
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(in) :: DP(:)
      real(8), intent(in) :: ETA2(:)
      real(8), intent(in) :: X(:)
      real(8), intent(in) :: Y(:)
      real(8), intent(in) :: H0
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: TVW(:)
      real(8), intent(inout) :: QN1(:)
      real(8), intent(inout) :: QN2(:)
      integer, intent(inout) :: ISSUBMERGED64(:)
      integer, intent(inout) :: ISSUBMERGED64P(:)

      integer :: I, K

      if (NFLUXIB > 0) NIBNODECODE = 0

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (30)
            call apply_radiation_boundary_flux(I, NVELL(K), NBV, UU1, VV1, CSII, SIII, H2, QN1)
            I = I + NVELL(K)
         case (40, 41)
            call apply_zero_normal_velocity_gradient_flux(I, NVELL(K), NBV, UU1, VV1, CSII, SIII, H2, QN1)
            I = I + NVELL(K)
         case (3, 13, 23)
            call apply_external_weir_boundary_flux(I, NVELL(K), TIMELOC, TVW, QN2)
            I = I + NVELL(K)
         case (4, 24)
            call apply_internal_weir_boundary_flux(K, I, NVELL(K), TIMELOC, TVW, NIBNODECODE, QN2)
            I = I + NVELL(K)*2
         case (5, 25)
            call apply_internal_weir_with_pipes_flux(K, I, NVELL(K), TIMELOC, TVW, NIBNODECODE, QN2)
            I = I + NVELL(K)*2
         case (64)
            call check_vew1d_submerged(K, I, NVELL(K), NBV, TIMELOC, X, Y, ETA2, DP, IBCONN, &
                                       LBArray_Pointer, NeiTab, NNeigh, NODECODE, H0, IFNLFA, ISSUBMERGED64, ISSUBMERGED64P)
            call apply_vew1d_boundary_flux(K, I, NVELL(K), TIMELOC, TVW, QN2)
            I = I + NVELL(K)*2
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine apply_computed_normal_flow_boundary_conditions

   subroutine apply_radiation_boundary_flux(I_in, NVELL, NBV, UU1, VV1, CSII, SIII, H2, QN1)

      implicit none

      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      integer, intent(in) :: NBV(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(inout) :: QN1(:)

      integer :: I, J, NNBB

      I = I_in
      do J = 1, NVELL
         I = I + 1
         NNBB = NBV(I)
         QN1(I) = compute_radiation_boundary_flux(UU1(NNBB), VV1(NNBB), CSII(I), SIII(I), H2(NNBB))
      end do

   end subroutine apply_radiation_boundary_flux

   subroutine apply_zero_normal_velocity_gradient_flux(I_in, NVELL, NBV, UU1, VV1, CSII, SIII, H2, QN1)

      implicit none

      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      integer, intent(in) :: NBV(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(inout) :: QN1(:)

      integer :: I, J, NNBB

      I = I_in
      do J = 1, NVELL
         I = I + 1
         NNBB = NBV(I)
         QN1(I) = compute_zero_normal_velocity_gradient_flux(UU1(NNBB), VV1(NNBB), CSII(I), SIII(I), H2(NNBB))
      end do

   end subroutine apply_zero_normal_velocity_gradient_flux

   subroutine apply_external_weir_boundary_flux(I_in, NVELL, TIMELOC, TVW, QN2)

      use mod_weir_flow, only: compute_external_boundary_flux

      implicit none

      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I, J

      I = I_in
      do J = 1, NVELL
         I = I + 1
         QN2(I) = compute_external_boundary_flux(I, TIMELOC, TVW)
      end do

   end subroutine apply_external_weir_boundary_flux

   subroutine apply_internal_weir_boundary_flux(K, I_in, NVELL, TIMELOC, TVW, NIBNODECODE, QN2)

      use mod_weir_flow, only: compute_internal_boundary_flux

      implicit none

      integer, intent(in) :: K
      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I, J

      I = I_in
      do J = 1, NVELL*2
         I = I + 1
         QN2(I) = compute_internal_boundary_flux(I, J, K, TIMELOC, NIBNODECODE, TVW)
      end do

   end subroutine apply_internal_weir_boundary_flux

   subroutine apply_internal_weir_with_pipes_flux(K, I_in, NVELL, TIMELOC, TVW, NIBNODECODE, QN2)

      use mod_weir_flow, only: compute_internal_boundary_flux, compute_cross_barrier_pipe_flux

      implicit none

      integer, intent(in) :: K
      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I, J

      I = I_in
      do J = 1, NVELL*2
         I = I + 1
         QN2(I) = compute_internal_boundary_flux(I, J, K, TIMELOC, NIBNODECODE, TVW) + &
                  compute_cross_barrier_pipe_flux(I, NIBNODECODE)
      end do

   end subroutine apply_internal_weir_with_pipes_flux

   !*******************************************************************************
   !> Check if nodes along vew1d boundaries are submerged
   !*******************************************************************************
   subroutine check_vew1d_submerged(K, I_in, NVELL, NBV, TIMELOC, X, Y, ETA2, DP, &
                                    IBCONN, LBArray_Pointer, NeiTab, NNeigh, &
                                    NODECODE, H0, IFNLFA, ISSUBMERGED64, ISSUBMERGED64P)
      use mod_weir_flow, only: SET_SUBMERGED64_AT

      implicit none

      integer, intent(in) :: K
      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: LBArray_Pointer(:)
      integer, intent(in) :: NeiTab(:, :)
      integer, intent(in) :: NNeigh(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: IFNLFA
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: X(:)
      real(8), intent(in) :: Y(:)
      real(8), intent(in) :: ETA2(:)
      real(8), intent(in) :: DP(:)
      real(8), intent(in) :: H0
      integer, intent(inout) :: ISSUBMERGED64(:)
      integer, intent(inout) :: ISSUBMERGED64P(:)

      integer :: I, J
      integer :: NNBB1, NNBBNEI, NNBB2, NWETNEI
      real(8) :: HTOT, X1, Y1, ET1, X2, Y2, ET2, LEN, SLP
      logical :: UNSUBMERGE
      integer :: JJ, II

      I = I_in

      do J = 1, NVELL
         I = I + 1
         call SET_SUBMERGED64_AT(I, J, K, TIMELOC)
         ISSUBMERGED64P(I) = ISSUBMERGED64(I)
      end do

      ! Cancelling ISSUBMERGED64 flag if any of the adjacent nodes along the boundary is not submerged. sb 8/9/2024
      do J = 1, NVELL
         I = I + 1
         NNBB1 = NBV(I)
         HTOT = ETA2(NNBB1) + IFNLFA*DP(NNBB1)
         if (ISSUBMERGED64(I) /= 0 .and. &
             NODECODE(NNBB1) /= 0 .and. &
             HTOT < 4.d0*H0) then
            UNSUBMERGE = .false.
            X1 = X(NNBB1)
            Y1 = Y(NNBB1)
            ET1 = ETA2(NNBB1)
            if (J > 1) then
               NNBBNEI = NBV(I - 1)
               X2 = X(NNBBNEI)
               Y2 = Y(NNBBNEI)
               ET2 = ETA2(NNBBNEI)
               LEN = sqrt((X1 - X2)**2 + (Y1 - Y2)**2)
               SLP = abs(ET1 - ET2)/LEN
               if (SLP > 0.0001) then
                  UNSUBMERGE = .true.
               end if
            elseif (J < NVELL) then
               NNBBNEI = NBV(I + 1)
               X2 = X(NNBBNEI)
               Y2 = Y(NNBBNEI)
               ET2 = ETA2(NNBBNEI)
               LEN = sqrt((X1 - X2)**2 + (Y1 - Y2)**2)
               SLP = abs(ET1 - ET2)/LEN
               if (SLP > 0.0001) then
                  UNSUBMERGE = .true.
               end if
            end if
            if (.not. UNSUBMERGE) then
               NWETNEI = 0
               do JJ = 2, NNeigh(NNBB1)
                  II = NeiTab(NNBB1, JJ)
                  if (NODECODE(II) /= 0) then
                     NWETNEI = NWETNEI + 1
                  end if
               end do
               if (NWETNEI < 3) then
                  UNSUBMERGE = .true.
               end if
            else
               ISSUBMERGED64(I) = 0
               NNBB2 = IBCONN(I)
               if (NNBB2 > 0) then
                  ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
               end if
            end if
         end if
      end do

   end subroutine check_vew1d_submerged

   subroutine apply_vew1d_boundary_flux(K, I_in, NVELL, TIMELOC, TVW, QN2)

      use mod_weir_flow, only: compute_internal_boundary64_flux

      implicit none

      integer, intent(in) :: K
      integer, intent(in) :: I_in
      integer, intent(in) :: NVELL
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I, J

      I = I_in
      do J = 1, NVELL*2
         I = I + 1
         QN2(I) = compute_internal_boundary64_flux(I, J, K, TIMELOC, TVW)
      end do

   end subroutine apply_vew1d_boundary_flux

   !*******************************************************************************
   !> Compute radiation boundary flux
   !>
   !> @param[in] UU1 Velocity component in x-direction
   !> @param[in] VV1 Velocity component in y-direction
   !> @param[in] CS Cosine of the angle
   !> @param[in] SI Sine of the angle
   !> @param[in] H2 Depth
   !> @return Q Radiation boundary flux
   !*******************************************************************************
   real(8) pure function compute_radiation_boundary_flux(UU1, VV1, CS, SI, H2) result(Q)
      implicit none

      real(8), intent(in) :: UU1
      real(8), intent(in) :: VV1
      real(8), intent(in) :: CS
      real(8), intent(in) :: SI
      real(8), intent(in) :: H2

      Q = H2*(UU1*CS + VV1*SI)

   end function compute_radiation_boundary_flux

   !*******************************************************************************
   !> Compute zero normal velocity gradient flux
   !>
   !> @param[in] UU1 Velocity component in x-direction
   !> @param[in] VV1 Velocity component in y-direction
   !> @param[in] CS Cosine of the angle
   !> @param[in] SI Sine of the angle
   !> @param[in] H2 Depth
   !> @return Q Zero normal velocity gradient flux
   !*******************************************************************************
   real(8) pure function compute_zero_normal_velocity_gradient_flux(UU1, VV1, CS, SI, H2) result(Q)

      real(8), intent(in) :: UU1
      real(8), intent(in) :: VV1
      real(8), intent(in) :: CS
      real(8), intent(in) :: SI
      real(8), intent(in) :: H2

      Q = H2*(UU1*CS + VV1*SI)

   end function compute_zero_normal_velocity_gradient_flux

   !*******************************************************************************
   !> Compute normal flow boundary flux
   !>
   !> @param[in] IT Current time step
   !> @param[in] NFLUXF Number of normal flow boundary fluxes
   !> @param[in] NFFR Number of harmonic constituents
   !> @param[in] NVEL Number of velocity components
   !> @param[in] FluxSettlingIT Time step for flux settling
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] NBV Array of boundary vertices
   !> @param[in] timeh Current time
   !> @param[in] TimeLoc Current time
   !> @param[in] FTIMINC Time increment
   !> @param[in] RampExtFlux Ramp factor for external forcing
   !> @param[in] QNPH Array of phase angles for normal flow boundary fluxes
   !> @param[in] QNAM Array of amplitudes for normal flow boundary fluxes
   !> @param[in] FPER Array of periods for normal flow boundary fluxes
   !> @param[in] FAMIG Array of phase angles for normal flow boundary fluxes
   !> @param[in] FFACE Array of phase angles for normal flow boundary fluxes
   !> @param[in] FFF Array of phase angles for normal flow boundary fluxes
   !> @param[in] QNIN1 Array of normal flow boundary fluxes at time QTIME1
   !> @param[in] QNIN2 Array of normal flow boundary fluxes at time QTIME2
   !> @param[in] ENPH Array of phase angles for normal flow boundary fluxes
   !> @param[in] ENAM Array of amplitudes for normal flow boundary fluxes
   !> @param[in,out] QTIME1 Time of the first normal flow boundary flux
   !> @param[in,out] QTIME2 Time of the second normal flow boundary flux
   !> @param[in,out] ENIN1 Array of normal flow boundary fluxes at time QTIME1
   !> @param[in,out] ENIN2 Array of normal flow boundary fluxes at time QTIME2
   !> @param[in,out] Eta2 Array of water surface elevations
   !> @param[in,out] EtaDisc Array of water surface elevations at the discontinuity
   !> @param[in,out] QN2 Array of normal flow boundary fluxes
   !> @param[in,out] EN2 Array of normal flow boundary fluxes
   !> @param[in,out] ElevDisc Array of water surface elevations at the discontinuity
   !> @param[in,out] EtaDisc_Fill Flag to indicate if the water surface elevation at the discontinuity is filled
   !*******************************************************************************
   subroutine apply_specified_normal_flow_boundary_conditions(IT, NFLUXF, &
                                                              NFFR, NVEL, FluxSettlingIT, &
                                                              LBCODEI, NBV, timeh, TimeLoc, FTIMINC, &
                                                              RampExtFlux, QNPH, QNAM, FPER, FAMIG, &
                                                              FFACE, FFF, QNIN1, QNIN2, ENPH, ENAM, &
                                                              QTIME1, QTIME2, ENIN1, ENIN2, Eta2, &
                                                              EtaDisc, QN2, EN2, ElevDisc, &
                                                              EtaDisc_Fill)
      implicit none

      integer, intent(in) :: IT
      integer, intent(in) :: NFLUXF
      integer, intent(in) :: NFFR
      integer, intent(in) :: NVEL
      integer, intent(in) :: FluxSettlingIT
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      real(8), intent(in) :: timeh
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: FTIMINC
      real(8), intent(in) :: RampExtFlux
      real(8), intent(in) :: QNPH(:, :)
      real(8), intent(in) :: QNAM(:, :)
      real(8), intent(in) :: FPER(:)
      real(8), intent(in) :: FAMIG(:)
      real(8), intent(in) :: FFACE(:)
      real(8), intent(in) :: FFF(:)
      real(8), intent(in) :: ENPH(:, :)
      real(8), intent(in) :: ENAM(:, :)
      real(8), intent(in) :: Eta2(:)
      real(8), intent(inout) :: QTIME1
      real(8), intent(inout) :: QTIME2
      real(8), intent(inout) :: ENIN1(:)
      real(8), intent(inout) :: ENIN2(:)
      real(8), intent(inout) :: QNIN1(:)
      real(8), intent(inout) :: QNIN2(:)
      real(8), intent(inout) :: EtaDisc(:)
      real(8), intent(inout) :: QN2(:)
      real(8), intent(inout) :: EN2(:)
      real(8), intent(inout) :: ElevDisc(:)
      logical, intent(inout) :: EtaDisc_Fill

      integer :: I
      integer :: NNBB

      if (NFLUXF == 1) then
         if (NFFR > 0) then
            call apply_periodic_normal_flow_boundary_conditions(NFFR, NVEL, LBCODEI, timeh, RampExtFlux, &
                                                                FPER, FAMIG, FFACE, FFF, QNPH, QNAM, ENPH, ENAM, &
                                                                QN2, EN2)
         elseif (NFFR == 0 .or. NFFR == -1) then
            call apply_normal_flow_boundary_conditions_from_file(NVEL, TimeLoc, FTIMINC, &
                                                                 RampExtFlux, LBCODEI, QTIME1, QTIME2, &
                                                                 QNIN1, QNIN2, ENIN1, ENIN2, QN2, EN2)
         end if

         if (IT == FluxSettlingIT) then
            EtaDisc_Fill = .false. ! sb v46.48 11/06/2006
            EtaDisc = Eta2
            do I = 1, NVEL
               if (LBCODEI(I) == 52) then
                  ElevDisc(I) = Eta2(NBV(I))
               end if
            end do
         else if (EtaDisc_Fill .and. IT > FluxSettlingIT) then
            EtaDisc_Fill = .false.
            do I = 1, NVEL
               if (LBCODEI(I) == 52) then
                  NNBB = NBV(I)
                  ElevDisc(I) = EtaDisc(NNBB) ! sb v46.48 11/06/2006
               end if
            end do
         end if

      end if

   end subroutine apply_specified_normal_flow_boundary_conditions

   !*******************************************************************************
   !> Apply normal flow boundary conditions from file
   !>
   !> @param[in] NVEL Number of velocity components
   !> @param[in] TimeLoc Current time
   !> @param[in] FTIMINC Time increment
   !> @param[in] RampExtFlux Ramp factor for external forcing
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in,out] QTIME1
   !> @param[in,out] QTIME2
   !> @param[in,out] QNIN1
   !> @param[in,out] QNIN2
   !> @param[in,out] ENIN1
   !> @param[in,out] ENIN2
   !> @param[in,out] QN2
   !> @param[in,out] EN2
   !*******************************************************************************
   subroutine apply_normal_flow_boundary_conditions_from_file(NVEL, TimeLoc, FTIMINC, &
                                                              RampExtFlux, LBCODEI, QTIME1, QTIME2, &
                                                              QNIN1, QNIN2, ENIN1, ENIN2, QN2, EN2)

      implicit none

      integer, intent(in) :: NVEL
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: FTIMINC
      real(8), intent(in) :: RampExtFlux
      integer, intent(in) :: LBCODEI(:)
      real(8), intent(inout) :: QTIME1
      real(8), intent(inout) :: QTIME2
      real(8), intent(inout) :: QNIN1(:)
      real(8), intent(inout) :: QNIN2(:)
      real(8), intent(inout) :: ENIN1(:)
      real(8), intent(inout) :: ENIN2(:)
      real(8), intent(inout) :: QN2(:)
      real(8), intent(inout) :: EN2(:)

      integer :: J
      real(8) :: QTRATIO

      if (TimeLoc > QTIME2) then
         QTIME1 = QTIME2
         QTIME2 = QTIME2 + FTIMINC
         QNIN1 = QNIN2
         ENIN1 = ENIN2

         do J = 1, NVEL
            select case (LBCODEI(J))
            case (2, 12, 22)
               read (20, *) QNIN2(J)
            case (32)
               read (20, *) QNIN2(J), ENIN2(J)
            end select
         end do
      end if

      QTRATIO = (TimeLoc - QTIME1)/FTIMINC
      QN2 = RampExtFlux*(QNIN1 + QTRATIO*(QNIN2 - QNIN1))
      EN2 = RampExtFlux*(ENIN1 + QTRATIO*(ENIN2 - ENIN1))
   end subroutine apply_normal_flow_boundary_conditions_from_file

   !*******************************************************************************
   !> Apply periodic normal flow boundary conditions from harmonic constituents
   !>
   !> @param[in] NFFR Number of harmonic constituents
   !> @param[in] NVEL Number of velocity components
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] timeh Current time
   !> @param[in] RampExtFlux Ramp factor for external forcing
   !> @param[in] FPER
   !> @param[in] FAMIG
   !> @param[in] FFACE
   !> @param[in] FFF
   !> @param[in] QNPH
   !> @param[in] QNAM
   !> @param[in] ENPH
   !> @param[in] ENAM
   !> @param[in,out] QN2
   !> @param[in,out] EN2
   !*******************************************************************************
   subroutine apply_periodic_normal_flow_boundary_conditions(NFFR, NVEL, LBCODEI, timeh, RampExtFlux, &
                                                             FPER, FAMIG, FFACE, FFF, QNPH, QNAM, ENPH, ENAM, &
                                                             QN2, EN2)

      implicit none

      real(8), parameter  :: eps = epsilon(1.d0)

      integer, intent(in) :: NFFR
      integer, intent(in) :: NVEL
      integer, intent(in) :: LBCODEI(:)
      real(8), intent(in) :: timeh
      real(8), intent(in) :: RampExtFlux
      real(8), intent(in) :: FPER(:)
      real(8), intent(in) :: FAMIG(:)
      real(8), intent(in) :: FFACE(:)
      real(8), intent(in) :: FFF(:)
      real(8), intent(in) :: QNPH(:, :)
      real(8), intent(in) :: QNAM(:, :)
      real(8), intent(in) :: ENPH(:, :)
      real(8), intent(in) :: ENAM(:, :)
      real(8), intent(inout) :: QN2(:)
      real(8), intent(inout) :: EN2(:)

      integer :: I, J
      integer :: NCYC
      real(8) :: ARG, ARGJ
      real(8) :: RFF

      do J = 1, NFFR
         ! Check if the period is zero
         if (abs(FPER(J)) <= eps) then
            NCYC = 0
         else
            NCYC = int(timeh/FPER(J))
         end if

         ARGJ = FAMIG(J)*(timeh - NCYC*FPER(J)) + FFACE(J)
         RFF = FFF(J)*RampExtFlux

         do I = 1, NVEL
            ARG = ARGJ - QNPH(J, I)
            QN2(I) = QN2(I) + QNAM(J, I)*RFF*cos(ARG)
            if (LBCODEI(I) == 32) then
               ARG = ARGJ - ENPH(J, I)
               EN2(I) = EN2(I) + ENAM(J, I)*RFF*cos(ARG)
            end if
         end do
      end do

   end subroutine apply_periodic_normal_flow_boundary_conditions

   !*******************************************************************************
   !> Enforce subdomain boundaries
   !>
   !> @param[in] it Current time step
   !> @param[in] enforceBN Flag to indicate which boundaries to enforce
   !*******************************************************************************
   subroutine enforce_subdomain_boundaries(it, enforceBN)
      use subdomain, only: readFort019, readFort020, readFort021

      implicit none

      integer, intent(in) :: it
      integer, intent(in) :: enforceBN

      if (enforceBN == 1) then
         call readFort019(it)
      elseif (enforceBN == 2) then
         call readFort020(it)
         call readFort021(it)
      end if

   end subroutine enforce_subdomain_boundaries

end module mod_timestep_bc_forcing
