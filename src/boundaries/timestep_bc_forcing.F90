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

      use mesh, only: NP, LBArray_Pointer, NeiTab, NNeigh, X, Y
      use boundaries, only: NBOU, NFLUXB, NFLUXRBC, NFLUXF, NFLUXGBC, NFLUXIB, NFLUXIBP, LBCODEI, &
                            NBV, NVEL, NVELL, CSII, SIII, IBCONN, NFLUXIB64
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
      call apply_radiation_boundary_discharge(NFLUXRBC, NVEL, LBCODEI, &
                                              NBV, UU1, VV1, CSII, SIII, H2, &
                                              QN1)
      call apply_zero_normal_velocity_gradient_boundary_discharge(NFLUXGBC, NVEL, LBCODEI, &
                                                                  NBV, UU1, VV1, CSII, SIII, H2, &
                                                                  QN1)
      call apply_exterior_weir_boundary_discharge(NFLUXB, NVEL, LBCODEI, TIMELOC, TVW, QN2)
      call check_vew1d_submerged(NFLUXIB64, NFLUXIB, NFLUXIBP, NP, NBOU, NVELL, LBCODEI, NBV, &
              TIMELOC, X, Y, ETA2, DP, IBCONN, LBArray_Pointer, NeiTab, NNeigh, &
              NODECODE, H0, IFNLFA, ISSUBMERGED64, ISSUBMERGED64P)
      call apply_interior_weir_boundary_discharge(NFLUXIB, NBOU, NVELL, LBCODEI, TIMELOC, NIBNODECODE, TVW, QN2)
      call apply_cross_barrier_pipe_discharge(NFLUXIBP, NVEL, LBCODEI, NIBNODECODE, QN2)
      call enforce_subdomain_boundaries(IT, subdomainOn, enforceBN)

   end subroutine apply_boundary_conditions

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
         elseif ((NFFR == 0) .or. (NFFR == -1)) then
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
   !> Apply radiation boundary discharge
   !>
   !> @param[in] NFLUXRBC Flag to indicate if there are radiation boundaries
   !> @param[in] NP Number of nodes in the mesh
   !> @param[in] NVEL Number of boundary nodes
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] NBV Array of boundary nodes
   !> @param[in] UU1 Array of x velocities
   !> @param[in] VV1 Array of y velocities
   !> @param[in] CSII Array of x components of the normal vector
   !> @param[in] SIII Array of y components of the normal vector
   !> @param[in] H2 Array of water depths
   !> @param[in,out] QN1 Array of boundary discharges
   !*******************************************************************************
   subroutine apply_radiation_boundary_discharge(NFLUXRBC, NVEL, LBCODEI, &
                                                 NBV, UU1, VV1, CSII, SIII, H2, &
                                                 QN1)

      implicit none

      integer, intent(in) :: NFLUXRBC
      integer, intent(in) :: NVEL
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(inout) :: QN1(:)

      integer :: J
      integer :: NNBB
      real(8) :: UN1

      if (NFLUXRBC == 1) then
         do J = 1, NVEL
            if (LBCODEI(J) == 30) then
               NNBB = NBV(J)
               UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
               QN1(J) = H2(NNBB)*UN1
            end if
         end do
      end if

   end subroutine apply_radiation_boundary_discharge

   !*******************************************************************************
   !> Apply zero normal velocity gradient boundary discharge
   !>
   !> @param[in] NFLUXGBC Flag to indicate if there are zero normal velocity gradient boundaries
   !> @param[in] NVEL Number of boundary nodes
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] NBV Array of boundary nodes
   !> @param[in] UU1 Array of x velocities
   !> @param[in] VV1 Array of y velocities
   !> @param[in] CSII Array of x components of the normal vector
   !> @param[in] SIII Array of y components of the normal vector
   !> @param[in] H2 Array of water depths
   !> @param[in,out] QN1 Array of boundary discharges
   !*******************************************************************************
   subroutine apply_zero_normal_velocity_gradient_boundary_discharge(NFLUXGBC, NVEL, LBCODEI, &
                                                                     NBV, UU1, VV1, CSII, SIII, H2, &
                                                                     QN1)
      implicit none

      integer, intent(in) :: NFLUXGBC
      integer, intent(in) :: NVEL
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(inout) :: QN1(:)

      integer :: J
      integer :: NNBB
      real(8) :: UN1

      if (NFLUXGBC == 1) then
         do J = 1, NVEL
            if ((LBCODEI(J) == 40) .or. (LBCODEI(J) == 41)) then
               NNBB = NBV(J)
               UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
               QN1(J) = H2(NNBB)*UN1
            end if
         end do
      end if

   end subroutine apply_zero_normal_velocity_gradient_boundary_discharge

   !*******************************************************************************
   !> Apply exterior weir boundary discharge
   !>
   !> @param[in] NFLUXB Flag to indicate if there are exterior weir boundaries
   !> @param[in] NVEL Number of boundary nodes
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] TIMELOC Current time
   !> @param[in,out] QN2 Array of boundary discharges
   !*******************************************************************************
   subroutine apply_exterior_weir_boundary_discharge(NFLUXB, NVEL, LBCODEI, TIMELOC, TVW, QN2)
      use mod_weir_flow, only: COMPUTE_EXTERNAL_BOUNDARY_FLUX

      implicit none

      integer, intent(in) :: NFLUXB
      integer, intent(in) :: NVEL
      integer, intent(in) :: LBCODEI(:)
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I

      if (NFLUXB == 1) then
         do I = 1, NVEL
            select case (LBCODEI(I))
            case (3, 13, 23)
               QN2(I) = COMPUTE_EXTERNAL_BOUNDARY_FLUX(I, TIMELOC, TVW)
            end select
         end do
      end if

   end subroutine apply_exterior_weir_boundary_discharge

   !*******************************************************************************
   !>  COMPUTE INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED INTERNAL BARRIER
   !>  BOUNDARY (PERMEABLE OR NOT) NODES
   !>
   !>  NFLUXIB is set to 1 in read_input.F if there are internal barrier
   !>  boundaries in the fort.14 (mesh) file.
   !>
   !>  IBSTART is a flag that indicates the first time through the time
   !>  stepping loop; set to 0 in read_input.F and set to 1 here.
   !>
   !>  NIBNODECODE seems to be set to 1 for nodes receiving water across
   !>  the barrier
   !>
   !>  BARMIN is used in several places, mainly as the minimum elevation
   !>  above the levee for flow to occur. It is a parameter and is set to
   !>  0.04 in global.F.
   !>
   !> @param[in] NFLUXIB Flag to indicate if there are internal barrier boundaries
   !> @param[in] NBOU Number of boundaries
   !> @param[in] NVELL Number of boundary nodes per boundary
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] TIMELOC Current time
   !> @param[in] NIBNODECODE Array of internal barrier boundary codes
   !> @param[in,out] QN2 Array of boundary discharges
   !*******************************************************************************
   subroutine apply_interior_weir_boundary_discharge(NFLUXIB, NBOU, NVELL, LBCODEI, TIMELOC, NIBNODECODE, TVW, QN2)
      use mod_weir_flow, only: COMPUTE_INTERNAL_BOUNDARY_FLUX, COMPUTE_INTERNAL_BOUNDARY64_FLUX

      implicit none

      integer, intent(in) :: NFLUXIB
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      real(8), intent(in) :: TIMELOC
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: TVW(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I, J, K
      logical :: ISFRONT

      if (NFLUXIB == 1) then
         NIBNODECODE(:) = 0
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (4, 24, 5, 25)
               do J = 1, NVELL(K)*2
                  I = I + 1
                  QN2(I) = COMPUTE_INTERNAL_BOUNDARY_FLUX(I, J, K, TIMELOC, NIBNODECODE, TVW)
               end do
            case (64)
               do J = 1, NVELL(K)*2
                  I = I + 1
                  if (J <= NVELL(K)) then
                     ISFRONT = .true.
                  else
                     ISFRONT = .false.
                  end if
                  QN2(I) = COMPUTE_INTERNAL_BOUNDARY64_FLUX(I, J, K, TIMELOC, TVW)
               end do
            case DEFAULT
               I = I + NVELL(K)
            end select
         end do
      end if

   end subroutine apply_interior_weir_boundary_discharge

   !*******************************************************************************
   !> COMPUTE INWARD/OUTWARD NORMAL FLOW FOR INTERNAL BARRIER
   !> BOUNDARY NODES THROUGH CROSS BARRIER PIPES
   !> NOTE THAT THIS ADDS AN ADDITIONAL FLOW COMPONENT INTO QN2
   !>
   !> @param[in] NFLUXIBP Flag to indicate if there are internal barrier boundaries
   !> @param[in] NVEL Number of boundary nodes
   !> @param[in] LBCODEI Array of boundary condition codes
   !> @param[in] TIMELOC Current time
   !> @param[in,out] QN2 Array of boundary discharges
   !*******************************************************************************
   subroutine apply_cross_barrier_pipe_discharge(NFLUXIBP, NVEL, LBCODEI, NIBNODECODE, QN2)
      use mod_weir_flow, only: COMPUTE_CROSS_BARRIER_PIPE_FLUX

      implicit none

      integer, intent(in) :: NFLUXIBP
      integer, intent(in) :: NVEL
      integer, intent(in) :: LBCODEI(:)
      integer, intent(inout) :: NIBNODECODE(:)
      real(8), intent(inout) :: QN2(:)

      integer :: I

      if (NFLUXIBP == 1) then
         do I = 1, NVEL
            if ((LBCODEI(I) == 5) .or. (LBCODEI(I) == 25)) then
               QN2(I) = QN2(I) + COMPUTE_CROSS_BARRIER_PIPE_FLUX(I, NIBNODECODE)
            end if
         end do
      end if

   end subroutine apply_cross_barrier_pipe_discharge

   !*******************************************************************************
   !> Check if nodes along vew1d boundaries are
   subroutine check_vew1d_submerged(NFLUXIB64, NFLUXIB, NFLUXIBP, NP, NBOU, NVELL, LBCODEI, NBV, &
                                    TIMELOC, X, Y, ETA2, DP, IBCONN, LBArray_Pointer, NeiTab, NNeigh, &
                                    NODECODE, H0, IFNLFA, ISSUBMERGED64, ISSUBMERGED64P)
      use mod_weir_flow, only: SET_SUBMERGED64_AT
      implicit none

      integer, intent(in) :: NFLUXIB
      integer, intent(in) :: NFLUXIBP
      integer, intent(in) :: NFLUXIB64
      integer, intent(in) :: NP
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: LBArray_Pointer(:)
      integer, intent(in) :: NeiTab(:,:)
      integer, intent(in) :: NNeigh(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: IFNLFA
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: X(:)
      real(8), intent(in) :: Y(:)
      real(8), intent(in) :: ETA2(:)
      real(8), intent(in) :: DP(:)
      real(8), intent(in) :: H0
      integer, intent(inout) :: ISSUBMERGED64(:)
      integer, intent(inout) :: ISSUBMERGED64P(:)

      integer :: I, J, K
      integer :: NNBB1, NNBBNEI, NNBB2, NWETNEI
      real(8) :: HTOT, X1, Y1, ET1, X2, Y2, ET2, LEN, SLP
      logical :: UNSUBMERGE
      integer :: JJ, II

      if (NFLUXIB64 == 1) then

         ! Setting ISSUBMERGED64 flag
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case (64)
               do J = 1, NVELL(K)
                  I = I + 1
                  call SET_SUBMERGED64_AT(I, J, K, TIMELOC)
                  ISSUBMERGED64P(I) = ISSUBMERGED64(I)
               end do
               I = I + NVELL(K)
            case DEFAULT
               I = I + NVELL(K)
            end select
         end do

         ! Cancelling ISSUBMERGED64 flag if any of the adjacent nodes along the boundary is not submerged. sb 8/9/2024
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case (64)
               do J = 1, NVELL(K)
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
                     end if
                     if (J < NVELL(K)) then
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
                     end if
                     if (UNSUBMERGE) then
                        ISSUBMERGED64(I) = 0
                        NNBB2 = IBCONN(I)
                        if (NNBB2 > 0) then
                           ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
                        end if
                     end if
                  end if
               end do
               I = I + NVELL(K)
            case DEFAULT
               I = I + NVELL(K)
            end select
         end do
      end if

   end subroutine check_vew1d_submerged

   !*******************************************************************************
   !> Enforce subdomain boundaries
   !>
   !> @param[in] it Current time step
   !> @param[in] subdomainOn Flag to indicate if subdomain boundaries are on
   !> @param[in] enforceBN Flag to indicate which boundaries to enforce
   !*******************************************************************************
   subroutine enforce_subdomain_boundaries(it, subdomainOn, enforceBN)
      use subdomain, only: readFort019, readFort020, readFort021

      implicit none

      integer, intent(in) :: it
      logical, intent(in) :: subdomainOn
      integer, intent(in) :: enforceBN

      if (subdomainOn) then
         if (enforceBN == 1) then
            call readFort019(it)
         elseif (enforceBN == 2) then
            call readFort020(it)
            call readFort021(it)
         end if
      end if

   end subroutine enforce_subdomain_boundaries

end module mod_timestep_bc_forcing
