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
module mod_momentum_bc_forcing

   implicit none

   integer, parameter :: conservative = 0
   integer, parameter :: non_conservative = 1
   integer, parameter :: predictor_corrector = 2

   private

   public :: apply_velocity_boundary_conditions, apply_zero_normal_velocity_gradient, &
             apply_subdomain_boundaries, UPDATE_U_PERSLNODES, conservative, &
             non_conservative, predictor_corrector

contains

   !************************************************************************
   !> Apply boundary conditions to the momentum equations.
   !>
   !> Modify the momentum equations to impose velocity boundary conditions.
   !> In each case the equations are manipulated to maintain the LHS matrix
   !> structure of: AUV11=AUV22; AUV12=-AUV21.
   !>
   !************************************************************************
   subroutine apply_velocity_boundary_conditions(formulation, NODECODE, QN2, H2, &
                                                 TKM, TK, UU1, VV1, QX1, QY1, &
                                                 MOM_LV_X, MOM_LV_Y, &
                                                 AUV, AUV11, AUV12)

      use mesh, only: NeiTab, NeiTabEle, NNeigh, FDXE, FDYE, SFMXEle, &
                      SFMYEle, SFacEle, MJU, TotalArea
      use global, only: NOFF, IFSFM, ILUMP, flgNodesMultipliedByTotalArea
      use boundaries, only: NBV, SIII, CSII, LBCODEI, ME2GW, NVELME, NFLUXIB64_GBL, IBCONN, &
                            ISSUBMERGED64, NBOU, NVELL
      use nodalattributes, only: ListCondensedNodes, NListCondensedNodes, NNodesListCondensedNodes, &
                                 LoadCondensedNodes

      implicit none

      integer, intent(in) :: formulation
      integer, intent(in) :: NODECODE(:)
      real(8), intent(in) :: QN2(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(in) :: TKM(:, :)
      real(8), intent(in) :: TK(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: QX1(:)
      real(8), intent(in) :: QY1(:)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout), optional :: AUV11(:), AUV12(:)
      integer :: J

      do J = 1, NVELME
         select case (LBCODEI(ME2GW(J)))
         case (0:9)
            call essential_normal_flow_free_tangential_slip(formulation, J, ME2GW, NBV, NODECODE, QN2, &
                                                            H2, TKM, TK, SIII, CSII, MOM_LV_X, &
                                                            MOM_LV_Y, AUV, AUV11, AUV12)
         case (10:19)
            call essential_normal_flow_no_tangential_slip(formulation, J, ME2GW, NBV, NODECODE, QN2, &
                                                          H2, SIII, CSII, MOM_LV_X, MOM_LV_Y, AUV, &
                                                          AUV11, AUV12)
         case (41)
            call zero_normal_velocity_gradient(formulation, J, IFSFM, ME2GW, NBV, NODECODE, NOFF, &
                                               NeiTab, NeiTabEle, NNeigh, SIII, CSII, UU1, VV1, QX1, QY1, &
                                               SFacEle, SFMXEle, SFMYEle, FDXE, FDYE, MOM_LV_X, &
                                               MOM_LV_Y, AUV, AUV11, AUV12)
         end select
      end do

      if (formulation == non_conservative .and. ILUMP /= 0) then
         call apply_vew1d_and_condensed_nodes(NFLUXIB64_GBL, NBOU, &
                                              NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, NODECODE, &
                                              NListCondensedNodes, NNodesListCondensedNodes, &
                                              ListCondensedNodes, LoadCondensedNodes, MJU, TotalArea, &
                                              flgNodesMultipliedByTotalArea, AUV, MOM_LV_X, MOM_LV_Y)
      end if

   end subroutine apply_velocity_boundary_conditions

   !************************************************************************
   !> Essential normal flow with free tangential slip.
   !>
   !> This subroutine modifies the momentum equations to impose an essential
   !> normal flow boundary condition with free tangential slip.
   !>
   !> @param[in] formulation A flag to indicate whether the conservative or non-conservative formulation should be used.
   !> @param[in] J The index of the boundary node.
   !> @param[in] ME2GW The mapping from the boundary node from momentum equations to gwce
   !> @param[in] NBV The mapping from the global node index to the boundary node index.
   !> @param[in] NODECODE The node code of the boundary node (i.e. wet/dry state).
   !> @param[in] QN2 The flux at the boundary node.
   !> @param[in] H2 The depth at the boundary node.
   !> @param[in] TKM
   !> @param[in] TK
   !> @param[in] SIII The sine of the angle between the normal and the x-axis.
   !> @param[in] CSII The cosine of the angle between the normal and the x-axis.
   !> @param[in,out] MOM_LV_X The x-component of the momentum at the boundary node.
   !> @param[in,out] MOM_LV_Y The y-component of the momentum at the boundary node.
   !> @param[in,out] AUV
   !> @param[in,out] AUV11 Predictor-corrector AUV values
   !> @param[in,out] AUV12 Predictor-corrector AUV values
   !>
   !************************************************************************
   subroutine essential_normal_flow_free_tangential_slip(formulation, J, ME2GW, NBV, NODECODE, &
                                                         QN2, H2, TKM, TK, SIII, CSII, MOM_LV_X, &
                                                         MOM_LV_Y, AUV, AUV11, AUV12)

      implicit none

      real(8), parameter :: eps = epsilon(1.d0)
      real(8), parameter :: sqrt2 = sqrt(2.d0)

      integer, intent(in) :: formulation
      integer, intent(in) :: J
      integer, intent(in) :: ME2GW(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: NODECODE(:)
      real(8), intent(in) :: QN2(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(in) :: TKM(:, :)
      real(8), intent(in) :: TK(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout), optional :: AUV11(:), AUV12(:)

      integer :: I
      integer :: NBDI
      real(8) :: NCI
      real(8) :: VelNorm

      I = ME2GW(J)
      NBDI = NBV(I)
      NCI = dble(NODECODE(NBDI))

      if (formulation == conservative) then
         MOM_LV_X(NBDI) = (SIII(I)*MOM_LV_X(NBDI) - CSII(I)*MOM_LV_Y(NBDI))*NCI !Tangetial Eqn RHS
         MOM_LV_Y(NBDI) = -QN2(I)*NCI !Normal Eqn RHS
         AUV(1, NBDI) = AUV(1, NBDI)*SIII(I) - AUV(4, NBDI)*CSII(I)
         AUV(3, NBDI) = AUV(3, NBDI)*SIII(I) - AUV(2, NBDI)*CSII(I)
         AUV(4, NBDI) = CSII(I)
         AUV(2, NBDI) = SIII(I)
      elseif (formulation == non_conservative) then
         VelNorm = -QN2(I)/H2(NBDI)
         if (abs(adcirc_norm2(TKM(:, NBDI)) - sqrt2*TK(NBDI)) < eps) then
            ! WJP 03.6.2018 In the case of the symmetric matrix..
            ! Should be equivalent to the non-symmetric formula below
            ! but for consistency when testing using this formula helps
            ! to keep the same solution
            MOM_LV_X(NBDI) = (SIII(I)*MOM_LV_X(NBDI) - CSII(I)*MOM_LV_Y(NBDI) - VelNorm*AUV(3, NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI) = VelNorm*AUV(1, NBDI)*NCI !Normal Eqn RHS
            AUV(3, NBDI) = -CSII(I)*AUV(1, NBDI)
            AUV(4, NBDI) = -AUV(3, NBDI)
            AUV(1, NBDI) = SIII(I)*AUV(1, NBDI)
            AUV(2, NBDI) = AUV(1, NBDI)
         else
            ! WJP 02.24.2018 in the case of the non-symmetric matrix
            MOM_LV_X(NBDI) = (SIII(I)*MOM_LV_X(NBDI) - CSII(I)*MOM_LV_Y(NBDI))*NCI !Tangetial Eqn RHS
            MOM_LV_Y(NBDI) = VelNorm*NCI !Normal Eqn RHS
            AUV(1, NBDI) = AUV(1, NBDI)*SIII(I) - AUV(4, NBDI)*CSII(I)
            AUV(3, NBDI) = AUV(3, NBDI)*SIII(I) - AUV(2, NBDI)*CSII(I)
            AUV(4, NBDI) = CSII(I)
            AUV(2, NBDI) = SIII(I)
         end if
      elseif (formulation == predictor_corrector) then
         VelNorm = -QN2(I)/H2(NBDI)
         MOM_LV_X(NBDI) = (SIII(I)*MOM_LV_X(NBDI) - CSII(I)*MOM_LV_Y(NBDI) - VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
         MOM_LV_Y(NBDI) = VelNorm*AUV11(NBDI)*NCI !Normal Eqn RHS
         AUV12(NBDI) = -CSII(I)*AUV11(NBDI)
         AUV11(NBDI) = SIII(I)*AUV11(NBDI)
      end if

   end subroutine essential_normal_flow_free_tangential_slip

   !************************************************************************
   !> Essential normal flow with no tangential slip.
   !>
   !> This subroutine modifies the momentum equations to impose an essential
   !> normal flow boundary condition with no tangential slip.
   !>
   !> @param[in] formulation A flag to indicate whether the conservative or non-conservative formulation should be used.
   !> @param[in] J The index of the boundary node.
   !> @param[in] NVEL The number of boundary nodes.
   !> @param[in] ME2GW The mapping from the boundary node from momentum equations to gwce
   !> @param[in] NBV The mapping from the global node index to the boundary node index.
   !> @param[in] NODECODE The node code of the boundary node (i.e. wet/dry state).
   !> @param[in] QN2 The flux at the boundary node.
   !> @param[in] H2 The depth at the boundary node.
   !> @param[in] SIII The sine of the angle between the normal and the x-axis.
   !> @param[in] CSII The cosine of the angle between the normal and the x-axis.
   !> @param[in,out] MOM_LV_X The x-component of the momentum at the boundary node.
   !> @param[in,out] MOM_LV_Y The y-component of the momentum at the boundary node.
   !> @param[in,out] AUV
   !> @param[in,out] AUV11 Predictor-corrector AUV values
   !> @param[in,out] AUV12 Predictor-corrector AUV values
   !************************************************************************
   subroutine essential_normal_flow_no_tangential_slip(formulation, J, ME2GW, NBV, NODECODE, &
                                                       QN2, H2, SIII, CSII, MOM_LV_X, MOM_LV_Y, &
                                                       AUV, AUV11, AUV12)

      implicit none

      integer, intent(in) :: formulation
      integer, intent(in) :: J
      integer, intent(in) :: ME2GW(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: NODECODE(:)
      real(8), intent(in) :: QN2(:)
      real(8), intent(in) :: H2(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout), optional :: AUV11(:), AUV12(:)

      integer :: I
      integer :: NBDI
      real(8) :: NCI

      I = ME2GW(J)
      NBDI = NBV(I)
      NCI = dble(NODECODE(NBDI))

      if (int(NCI) == 0) then
         MOM_LV_X(NBDI) = 0.d0
         MOM_LV_Y(NBDI) = 0.d0
      else
         if (formulation == conservative) then
            MOM_LV_X(NBDI) = 0.0d0 !Tangential Eqn RHS
            MOM_LV_Y(NBDI) = -QN2(I)*NCI !Normal Eqn RHS
            AUV(1, NBDI) = SIII(I)
            AUV(2, NBDI) = SIII(I)
            AUV(3, NBDI) = -CSII(I)
            AUV(4, NBDI) = CSII(I)
         elseif (formulation == non_conservative) then
            MOM_LV_X(NBDI) = 0.0d0 !Tangential Eqn RHS
            MOM_LV_Y(NBDI) = -QN2(I)/H2(NBDI)*NCI !Normal Eqn RHS
            AUV(1, NBDI) = SIII(I)
            AUV(2, NBDI) = SIII(I)
            AUV(3, NBDI) = -CSII(I)
            AUV(4, NBDI) = CSII(I)
         elseif (formulation == predictor_corrector) then
            MOM_LV_X(NBDI) = 0.0d0 !Tangential Eqn RHS
            MOM_LV_Y(NBDI) = -QN2(I)/H2(NBDI)*NCI !Normal Eqn RHS
            AUV11(NBDI) = SIII(I)
            AUV12(NBDI) = -CSII(I)
         end if
      end if

   end subroutine essential_normal_flow_no_tangential_slip

   !************************************************************************
   !> Zero normal velocity gradient boundary condition using a Galerkin approximation to
   !> the normal derivatives. Note: this is fully explicit and therefore
   !> the velocity at the boundary is computed entirely from surrounding
   !> velocities at the previous time step.
   !>
   !> @param[in] formulation A flag to indicate whether the conservative or non-conservative formulation should be used.
   !> @param[in] J The index of the boundary node.
   !> @param[in] IFSFM The flag to indicate whether to use the SFM or SFMx/SFMy.
   !> @param[in] ME2GW The mapping from the boundary node from momentum equations to gwce
   !> @param[in] NBV The mapping from the global node index to the boundary node index.
   !> @param[in] NODECODE The node code of the boundary node (i.e. wet/dry state).
   !> @param[in] NOFF The offset for the element.
   !> @param[in] NeiTab The neighbor table.
   !> @param[in] NeiTabEle The element table.
   !> @param[in] NNeigh The number of neighbors.
   !> @param[in] SIII The sine of the angle between the normal and the x-axis.
   !> @param[in] CSII The cosine of the angle between the normal and the x-axis.
   !> @param[in] UU1 The x-component of the velocity.
   !> @param[in] VV1 The y-component of the velocity.
   !> @param[in] QX1 The x-component of the velocity at the previous time step.
   !> @param[in] QY1 The y-component of the velocity at the previous time step.
   !> @param[in,out] MOM_LV_X The x-component of the momentum at the boundary node.
   !> @param[in,out] MOM_LV_Y The y-component of the momentum at the boundary node.
   !> @param[in,out] AUV
   !> @param[in,out] AUV11 Predictor-corrector AUV values
   !> @param[in,out] AUV12 Predictor-corrector AUV values
   !************************************************************************
   subroutine zero_normal_velocity_gradient(formulation, J, IFSFM, ME2GW, NBV, &
                                            NODECODE, NOFF, NeiTab, NeiTabEle, &
                                            NNeigh, SIII, CSII, UU1, VV1, QX1, QY1, SFacEle, &
                                            SFMXEle, SFMYEle, FDXE, FDYE, &
                                            MOM_LV_X, MOM_LV_Y, AUV, AUV11, AUV12)
      implicit none

      integer, intent(in) :: formulation
      integer, intent(in) :: J
      integer, intent(in) :: IFSFM
      integer, intent(in) :: ME2GW(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: NOFF(:)
      integer, intent(in) :: NeiTab(:, :)
      integer, intent(in) :: NeiTabEle(:, :)
      integer, intent(in) :: NNeigh(:)
      real(8), intent(in) :: SIII(:)
      real(8), intent(in) :: CSII(:)
      real(8), intent(in) :: UU1(:)
      real(8), intent(in) :: VV1(:)
      real(8), intent(in) :: SFacEle(:)
      real(8), intent(in) :: SFMXEle(:)
      real(8), intent(in) :: SFMYEle(:)
      real(8), intent(in) :: FDXE(:, :)
      real(8), intent(in) :: FDYE(:, :)
      real(8), intent(in) :: QX1(:)
      real(8), intent(in) :: QY1(:)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout), optional :: AUV11(:), AUV12(:)

      integer :: I
      integer :: NBDI
      integer :: NCI
      integer :: NM1, NM2, NM3
      integer :: N, NEle, NCEle, NNFirst
      real(8) :: SFacAvg, SFmxAvg, SFmyAvg
      real(8) :: sfdxfac, sfdyfac
      real(8) :: FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      real(8) :: ZNGRHS1, ZNGRHS2, ZNGLHS

      I = ME2GW(J)
      NBDI = NBV(I)
      NCI = NODECODE(NBDI)

      NM1 = NBDI
      ZNGRHS1 = 0.d0 !Zero Norm Grad of U Eqn
      ZNGRHS2 = 0.d0 !Zero Norm Grad of V Eqn
      ZNGLHS = 0.d0
      NM2 = NeiTab(NBDI, 2) !operate on 1st neighbor
      NNFirst = NM2 !save these values until end

      do N = 3, NNeigh(NBDI) !operate on rest of neighbors
         NM3 = NM2 !shift previously computed values
         NM2 = NEITAB(NBDI, N) !select new neighbor to work on
         NEle = NeiTabEle(NBDI, N - 2) !element # defined by nodes NM1,NM2,NM3
         NCEle = NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)

         if ((NEle /= 0) .and. (NCEle /= 0)) then !if element is active, compute contribution
            SFacAvg = SFacEle(NEle)
            SFmxAvg = SFMXEle(NEle)
            SFmyAvg = SFMYEle(NEle)
            sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg
            sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg
            FDX1 = FDXE(1, NEle)*sfdxfac !c FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
            FDX2 = FDXE(2, NEle)*sfdxfac !c FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
            FDX3 = FDXE(3, NEle)*sfdxfac !c FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
            FDY1 = FDYE(1, NEle)*sfdyfac !c  FDY1=X(NM3)-X(NM2) !a1
            FDY2 = FDYE(2, NEle)*sfdyfac !c  FDY2=X(NM1)-X(NM3) !a2
            FDY3 = FDYE(3, NEle)*sfdyfac !c  FDY3=X(NM2)-X(NM1) !a3

            if (formulation == conservative) then
               ZNGRHS1 = ZNGRHS1 - (CSII(I)*FDX2 + SIII(I)*FDY2)*QX1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*QX1(NM3)
               ZNGRHS2 = ZNGRHS2 - (CSII(I)*FDX2 + SIII(I)*FDY2)*QY1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*QY1(NM3)
            elseif (formulation == non_conservative .or. formulation == predictor_corrector) then
               ZNGRHS1 = ZNGRHS1 - (CSII(I)*FDX2 + SIII(I)*FDY2)*UU1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*UU1(NM3)
               ZNGRHS2 = ZNGRHS2 - (CSII(I)*FDX2 + SIII(I)*FDY2)*VV1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*VV1(NM3)
            end if

            ZNGLHS = ZNGLHS + CSII(I)*FDX1 + SIII(I)*FDY1
         end if

      end do

      NM3 = NM2 !wrap back to beginning to get final contribution
      NM2 = NNFirst
      NEle = NeiTabEle(NBDI, NNeigh(NBDI) - 1)
      NCEle = NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)

      if ((NEle /= 0) .and. (NCEle /= 0)) then
         SFacAvg = SFacEle(NEle)
         SFmxAvg = SFMXEle(NEle)
         SFmyAvg = SFMYEle(NEle)
         sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg
         sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg
         FDX1 = FDXE(1, NEle)*sfdxfac
         FDX2 = FDXE(2, NEle)*sfdxfac
         FDX3 = FDXE(3, NEle)*sfdxfac
         FDY1 = FDYE(1, NEle)*sfdyfac
         FDY2 = FDYE(2, NEle)*sfdyfac
         FDY3 = FDYE(3, NEle)*sfdyfac

         if (formulation == conservative) then
            ZNGRHS1 = ZNGRHS1 - (CSII(I)*FDX2 + SIII(I)*FDY2)*QX1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*QX1(NM3)
            ZNGRHS2 = ZNGRHS2 - (CSII(I)*FDX2 + SIII(I)*FDY2)*QY1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*QY1(NM3)
         elseif (formulation == non_conservative .or. formulation == predictor_corrector) then
            ZNGRHS1 = ZNGRHS1 - (CSII(I)*FDX2 + SIII(I)*FDY2)*UU1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*UU1(NM3)
            ZNGRHS2 = ZNGRHS2 - (CSII(I)*FDX2 + SIII(I)*FDY2)*VV1(NM2) - (CSII(I)*FDX3 + SIII(I)*FDY3)*VV1(NM3)
         end if
         ZNGLHS = ZNGLHS + CSII(I)*FDX1 + SIII(I)*FDY1
      end if

      if (NCI == 0) then
         MOM_LV_X(NBDI) = 0.d0
         MOM_LV_Y(NBDI) = 0.d0
      else
         MOM_LV_X(NBDI) = ZNGRHS1/ZNGLHS
         MOM_LV_Y(NBDI) = ZNGRHS2/ZNGLHS
      end if

      if (formulation == conservative .or. formulation == non_conservative) then
         AUV(1, NBDI) = 1.d0
         AUV(2, NBDI) = 1.d0
         AUV(3, NBDI) = 0.d0
         AUV(4, NBDI) = 0.d0
      elseif (formulation == predictor_corrector) then
         AUV11(NBDI) = 1.d0
         AUV12(NBDI) = 0.d0
      end if

   end subroutine zero_normal_velocity_gradient

   !************************************************************************
   !> VEW1D (=1D channel) 08-11-2022 SB
   !> Nodal equations at IBTYPE=64 VEW boundary nodes and at condensed nodes are
   !> summed up together in the following part. First, the value on the front side,
   !> i.e., the floodplain side, is summed to the back side, i.e., the channel side.
   !> Secondly, the values at the condensed nodes which are most likely on channel
   !> beds are summed together. Notice that by this second step the nodes on channel
   !> bed hold the sum of all the values on the floodplain side and the grouped
   !> condensed nodes. And then finally, the value on the backside, i.e. on the
   !> channel bed, is copied to the front side, i.e., the floodplain side. Through
   !> this procedure, the values at the floodplain nodes and (condensed) channel nodes
   !> have the same values on both sides of the equations.
   !>
   !> @param[in] NFLUXIB64_GBL Number of 1D channel boundary nodes
   !> @param[in] NBOU Number of boundaries
   !> @param[in] NVELL Number of vertices on each boundary
   !> @param[in] LBCODEI Boundary code for each vertex
   !> @param[in] NBV Global node number on this side of the barrier
   !> @param[in] IBCONN Global node number on opposite side of the barrier
   !> @param[in] ISSUBMERGED64 Submerged flag for each vertex
   !> @param[in] NODECODE Node code for each vertex
   !> @param[in] NListCondensedNodes Number of condensed nodes
   !> @param[in] NNodesListCondensedNodes Number of nodes in each condensed node group
   !> @param[in] ListCondensedNodes List of condensed nodes
   !> @param[in] LoadCondensedNodes Flag to indicate whether to load condensed nodes
   !> @param[in] MJU Node code for each vertex
   !> @param[in] TotalArea Total area of the 1D channel boundary
   !> @param[inout] flgNodesMultipliedByTotalArea Flag to indicate whether the nodes are multiplied by the total area
   !> @param[inout] AUV AUV values
   !> @param[inout] MOM_LV_X X-momentum values
   !> @param[inout] MOM_LV_Y Y-momentum values
   !************************************************************************
   subroutine apply_vew1d_and_condensed_nodes(NFLUXIB64_GBL, NBOU, NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, &
                                              NodeCode, NListCondensedNodes, NNodesListCondensedNodes, &
                                              ListCondensedNodes, LoadCondensedNodes, MJU, TotalArea, &
                                              flgNodesMultipliedByTotalArea, AUV, MOM_LV_X, MOM_LV_Y)

#ifdef CMPI
      use messenger, only: UPDATEM4R, UPDATER
#endif

      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: ISSUBMERGED64(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(:)
      integer, intent(in) :: ListCondensedNodes(:, :)
      logical, intent(in) :: LoadCondensedNodes
      integer, intent(in) :: MJU(:)
      real(8), intent(in) :: TotalArea(:)
      integer, intent(inout) :: flgNodesMultipliedByTotalArea(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)

#ifdef CMPI
      real(8), allocatable :: DUMY1(:)
#endif

      if (NFLUXIB64_GBL > 0) then
         call vew1d_sum_front_to_back_side(NBOU, NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, &
                                           NODECODE, TotalArea, MJU, flgNodesMultipliedByTotalArea, AUV, MOM_LV_X, MOM_LV_Y)
      end if

      if (LoadCondensedNodes) then
         call condensed_nodes_sum(NListCondensedNodes, NNodesListCondensedNodes, ListCondensedNodes, &
                                  NodeCode, MJU, TotalArea, flgNodesMultipliedByTotalArea, AUV, MOM_LV_X, MOM_LV_Y)
      end if

#ifdef CMPI
      if (NFLUXIB64_GBL > 0 .or. LoadCondensedNodes) then
         call UPDATEM4R(AUV)
         call UPDATER(MOM_LV_X, MOM_LV_Y, DUMY1, 2)
      end if
#endif

      if (NFLUXIB64_GBL > 0) then
         call vew1d_copy_back_to_front_side(NBOU, NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, NODECODE, &
                                            AUV, MOM_LV_X, MOM_LV_Y)
      end if

   end subroutine apply_vew1d_and_condensed_nodes

   !************************************************************************
   !> Copy values from back side to front side of the 1D channel boundary
   !>
   !> @param[in] NBOU Number of boundaries
   !> @param[in] NVELL Number of vertices on each boundary
   !> @param[in] LBCODEI Boundary code for each vertex
   !> @param[in] NBV Global node number on this side of the barrier
   !> @param[in] IBCONN Global node number on opposite side of the barrier
   !> @param[in] ISSUBMERGED64 Submerged flag for each vertex
   !> @param[in] NODECODE Node code for each vertex
   !> @param[inout] AUV AUV values
   !> @param[inout] MOM_LV_X X-momentum values
   !> @param[inout] MOM_LV_Y Y-momentum values
   !************************************************************************
   subroutine vew1d_copy_back_to_front_side(NBOU, NVELL, LBCODEI, NBV, IBCONN, &
                                            ISSUBMERGED64, NODECODE, AUV, MOM_LV_X, MOM_LV_Y)

      implicit none

      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: ISSUBMERGED64(:)
      integer, intent(in) :: NODECODE(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)

      integer :: I, J, K
      integer :: NNBB1, NNBB2

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               if (ISSUBMERGED64(I) /= 0) then
                  NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                  NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
                  if (NODECODE(NNBB1) /= 0 .and. NODECODE(NNBB2) /= 0) then
                     AUV(:, NNBB1) = AUV(:, NNBB2)
                     MOM_LV_X(NNBB1) = MOM_LV_X(NNBB2)
                     MOM_LV_Y(NNBB1) = MOM_LV_Y(NNBB2)
                  end if
               end if
            end do
            I = I + NVELL(K)
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine vew1d_copy_back_to_front_side

   !************************************************************************
   !> Sum the condensed nodes
   !>
   !> @param[in] NListCondensedNodes Number of sets of condensed nodes
   !> @param[in] NNodesListCondensedNodes Number of nodes for each condensed node
   !> @param[in] ListCondensedNodes List of condensed nodes
   !> @param[in] NodeCode Node code for each vertex
   !> @param[in] MJU MJU values
   !> @param[in] TotalArea Total area for each vertex
   !> @param[inout] flgNodesMultipliedByTotalArea Flags for nodes multiplied by total areas
   !> @param[inout] AUV AUV values
   !> @param[inout] MOM_LV_X X-momentum values
   !> @param[inout] MOM_LV_Y Y-momentum values
   !************************************************************************
   subroutine condensed_nodes_sum(NListCondensedNodes, NNodesListCondensedNodes, &
                                  ListCondensedNodes, NodeCode, MJU, TotalArea, flgNodesMultipliedByTotalArea, &
                                  AUV, MOM_LV_X, MOM_LV_Y)

      implicit none

      real(8), parameter :: eps = epsilon(1.d0)
      logical, parameter :: CME_AreaInt_Corr = .true.
      logical, parameter :: CME_AreaInt_Orig = .false.

      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(:)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: NodeCode(:)
      integer, intent(in) :: MJU(:)
      real(8), intent(in) :: TotalArea(:)
      integer, intent(inout) :: flgNodesMultipliedByTotalArea(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)

      integer :: I, J, K, L
      real(8) :: TotalArea1

      do K = 1, NListCondensedNodes
         I = ListCondensedNodes(K, 1)
         if (I == 0) cycle
         if ((NODECODE(I) /= 0)) then
            ! 1) Mutiply LHS & RHS by total area at Node I
            if ((flgNodesMultipliedByTotalArea(I) == 0) .and. (abs(TotalArea(I)) > eps)) then
               if (CME_AreaInt_Corr) then !Correct area integration
                  TotalArea1 = TotalArea(I)
               elseif (CME_AreaInt_Orig) then !Original (incorrect) area integration
                  TotalArea1 = dble(MJU(I))
               end if
            else
               TotalArea1 = 1.d0
            end if
            AUV(:, I) = TotalArea1*AUV(:, I)
            MOM_LV_X(I) = TotalArea1*MOM_LV_X(I)
            MOM_LV_Y(I) = TotalArea1*MOM_LV_Y(I)

            ! 2) Sum them up
            do L = 2, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               if ((flgNodesMultipliedByTotalArea(J) == 0) .and. (abs(TotalArea(J)) > eps)) then
                  if (CME_AreaInt_Corr) then !Correct area integration
                     TotalArea1 = TotalArea(J)
                  elseif (CME_AreaInt_Orig) then !Original (incorrect) area integration
                     TotalArea1 = dble(MJU(J))
                  end if
               else
                  TotalArea1 = 1.d0
               end if
               AUV(:, I) = AUV(:, I) + TotalArea1*AUV(:, J)
               MOM_LV_X(I) = MOM_LV_X(I) + TotalArea1*MOM_LV_X(J)
               MOM_LV_Y(I) = MOM_LV_Y(I) + TotalArea1*MOM_LV_Y(J)
            end do

            ! 3) Distribute them
            do L = 2, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               AUV(:, J) = AUV(:, I)
               MOM_LV_X(J) = MOM_LV_X(I)
               MOM_LV_Y(J) = MOM_LV_Y(I)
            end do
         end if
      end do

   end subroutine condensed_nodes_sum

   !************************************************************************
   !> Sum front side values to back side of the 1D channel boundary
   !>
   !> @param[in] NBOU Number of boundaries
   !> @param[in] NVELL Number of vertices on each boundary
   !> @param[in] LBCODEI Boundary code for each vertex
   !> @param[in] NBV Global node number on this side of the barrier
   !> @param[in] IBCONN Global node number on opposite side of the barrier
   !> @param[in] ISSUBMERGED64 Submerged flag for each vertex
   !> @param[in] NODECODE Node code for each vertex
   !> @param[in] TotalArea Total area for each vertex
   !> @param[inout] flgNodesMultipliedByTotalArea Flags for nodes multiplied by total areas
   !> @param[inout] AUV AUV values
   !> @param[inout] MOM_LV_X X-momentum values
   !> @param[inout] MOM_LV_Y Y-momentum values
   !************************************************************************
   subroutine vew1d_sum_front_to_back_side(NBOU, NVELL, LBCODEI, NBV, IBCONN, &
                                           ISSUBMERGED64, NODECODE, TotalArea, MJU, &
                                           flgNodesMultipliedByTotalArea, &
                                           AUV, MOM_LV_X, MOM_LV_Y)

      implicit none

      real(8), parameter :: eps = epsilon(1.d0)
      logical, parameter :: CME_AreaInt_Corr = .true.
      logical, parameter :: CME_AreaInt_Orig = .false.

      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: ISSUBMERGED64(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: MJU(:)
      real(8), intent(in) :: TotalArea(:)
      integer, intent(inout) :: flgNodesMultipliedByTotalArea(:)
      real(8), intent(inout) :: AUV(:, :)
      real(8), intent(inout) :: MOM_LV_X(:)
      real(8), intent(inout) :: MOM_LV_Y(:)

      integer :: I, J, K
      integer :: NNBB1, NNBB2
      real(8) :: TotalArea1, TotalArea2

      flgNodesMultipliedByTotalArea(:) = 0 ! Initialize flags for nodes multiplied by total areas

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               if (ISSUBMERGED64(I) /= 0) then
                  NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                  NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER

                  if (NODECODE(NNBB1) /= 0 .and. NODECODE(NNBB2) /= 0) then
                     if ((flgNodesMultipliedByTotalArea(NNBB1) == 0) .and. (abs(TotalArea(NNBB1)) > eps)) then
                        if (CME_AreaInt_Corr) then !Correct area integration
                           TotalArea1 = TotalArea(NNBB1)
                        elseif (CME_AreaInt_Orig) then !Original (incorrect) area integration
                           TotalArea1 = dble(MJU(NNBB1))
                        end if
                     else
                        TotalArea1 = 1.d0
                     end if
                     if ((flgNodesMultipliedByTotalArea(NNBB2) == 0) .and. (abs(TotalArea(NNBB2)) > eps)) then
                        if (CME_AreaInt_Corr) then !Correct area integration
                           TotalArea2 = TotalArea(NNBB2)
                        elseif (CME_AreaInt_Orig) then !Original (incorrect) area integration
                           TotalArea2 = dble(MJU(NNBB2))
                        end if
                     else
                        TotalArea2 = 1.d0
                     end if

                     AUV(:, NNBB2) = TotalArea1*AUV(:, NNBB1) + TotalArea2*AUV(:, NNBB2)
                     MOM_LV_X(NNBB2) = TotalArea1*MOM_LV_X(NNBB1) + TotalArea2*MOM_LV_X(NNBB2)
                     MOM_LV_Y(NNBB2) = TotalArea1*MOM_LV_Y(NNBB1) + TotalArea2*MOM_LV_Y(NNBB2)

                     AUV(:, NNBB1) = 0.d0 ! Set zero to avoid duplicated additions
                     MOM_LV_X(NNBB1) = 0.d0 !
                     MOM_LV_Y(NNBB1) = 0.d0 !

                     flgNodesMultipliedByTotalArea(NNBB1) = 1
                     flgNodesMultipliedByTotalArea(NNBB2) = 1
                  end if
               end if
            end do
            I = I + NVELL(K)
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine vew1d_sum_front_to_back_side

   !************************************************************************
   !> Apply zero normal velocity gradient boundary condition
   !>
   !> @param[in] formulation flag to indicate whether the conservative or nonconservative form of the boundary condition is to be applied
   !> @param[in] NFLUXGBC Integer flag to indicate whether the boundary condition is to be applied
   !> @param[in] NVEL Integer number of velocity components
   !> @param[in] NVELME Integer number of velocity components in the mixed element
   !> @param[in] ME2GW Integer array mapping mixed element velocity components to global velocity components
   !> @param[in] NBV Integer array mapping velocity components to global nodes
   !> @param[in] LBCODEI Integer array of boundary condition codes
   !> @param[in] NEleZNG Integer array mapping mixed elements to global elements
   !> @param[in] NODECODE Integer array of node codes
   !> @param[in] NOFF Integer array of node offsets
   !> @param[in] NM Integer array of node mappings
   !> @param[in] ZNGIF1 Real array of fictitious point interpolation coefficients
   !> @param[in] ZNGIF2 Real array of fictitious point interpolation coefficients
   !> @param[in] ZNGIF3 Real array of fictitious point interpolation coefficients
   !> @param[in,out] QX2 Real array of x-momentum values
   !> @param[in,out] QY2 Real array of y-momentum values
   !> @param[in,out] UU2 Real array of x-velocity values
   !> @param[in,out] VV2 Real array of y-velocity values
   !************************************************************************
   subroutine apply_zero_normal_velocity_gradient(formulation, NODECODE, NOFF, &
                                                  QX2, QY2, UU2, VV2)

      use mesh, only: NM
      use boundaries, only: NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, NEleZNG, &
                            ZNGIF1, ZNGIF2, ZNGIF3

      implicit none

      integer, intent(in) :: formulation
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: NOFF(:)
      real(8), intent(inout) :: QX2(:)
      real(8), intent(inout) :: QY2(:)
      real(8), intent(inout) :: UU2(:)
      real(8), intent(inout) :: VV2(:)

      integer :: I, J
      integer :: NBDI, NM1, NM2, NM3, NC1, NC2, NC3
      real(8) :: NCEle

      if (NFLUXGBC == 1) then
         do J = 1, NVELME
            I = ME2GW(J)
            NBDI = NBV(I)
            if (LBCODEI(I) == 40) then
               NM1 = NM(NEleZNG(I), 1)
               NM2 = NM(NEleZNG(I), 2)
               NM3 = NM(NEleZNG(I), 3)
               NC1 = NODECODE(NM1)
               NC2 = NODECODE(NM2)
               NC3 = NODECODE(NM3)
               NCEle = dble(NC1*NC2*NC3*NOFF(NEleZNG(I)))
               if (formulation == conservative) then
                  QX2(NBDI) = NCEle*(QX2(NM1)*ZNGIF1(I) + QX2(NM2)*ZNGIF2(I) + QX2(NM3)*ZNGIF3(I))
                  QY2(NBDI) = NCEle*(QY2(NM1)*ZNGIF1(I) + QY2(NM2)*ZNGIF2(I) + QY2(NM3)*ZNGIF3(I))
               else
                  UU2(NBDI) = NCEle*(UU2(NM1)*ZNGIF1(I) + UU2(NM2)*ZNGIF2(I) + UU2(NM3)*ZNGIF3(I))
                  VV2(NBDI) = NCEle*(VV2(NM1)*ZNGIF1(I) + VV2(NM2)*ZNGIF2(I) + VV2(NM3)*ZNGIF3(I))
               end if
            end if
         end do
      end if

   end subroutine apply_zero_normal_velocity_gradient

   !************************************************************************
   !> @brief Internal norm2 function
   !> If the compiler supports Fortran 2008 then the intrinsic is used.
   !> Otherwise, the implementation is provided here. Use the compiler
   !> flag -DNOF2008 to enable the internal implementation. norm2
   !> calculates the Euclidean vector norm (L2 norm) of an array
   !>
   !> @param[in] x array to calculate the L2 norm of
   !> @result euclidean vector norm of x
   !************************************************************************
   real(8) pure function adcirc_norm2(x) result(v)
      implicit none
      real(8), intent(in) :: x(:)
#ifdef NOF2008
      intrinsic         :: dot_product, sqrt
      v = sqrt(dot_product(x, x))
#else
      intrinsic            :: norm2
      v = norm2(x)
#endif
   end function adcirc_norm2

   !************************************************************************
   !> Apply the subdomain boundaries for velocity components
   !>
   !> @param[in] subdomainOn Logical flag to indicate whether the subdomain boundaries are to be applied
   !> @param[in] enforceBN Integer flag to indicate whether the boundary nodes are to be enforced
   !************************************************************************
   subroutine apply_subdomain_boundaries(subdomainOn, enforceBN)

      use subdomain, only: enforceUVcb, enforceUVob

      implicit none

      logical, intent(in) :: subdomainOn
      integer, intent(in) :: enforceBN

      if (subdomainOn) then
         if (enforceBN == 1) then
            call enforceUVcb()
         elseif (enforceBN == 2) then
            call enforceUVob()
         end if
      end if

   end subroutine apply_subdomain_boundaries

   !************************************************************************
   !> Update the U and V values at the nodes on the periodic sponge layer
   !> boundary
   !>
   !> @param[in] NNPERBC Number of nodes on the periodic boundary
   !> @param[in] IPERCONN Integer array of periodic node connections
   !> @param[in,out] UU2 Real array of x-velocity values
   !> @param[in,out] VV2 Real array of y-velocity values
   !************************************************************************
   subroutine UPDATE_U_PERSLNODES(NNPERBC, IPERCONN, UU2, VV2)
      implicit none

      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      real(8), intent(inout) :: UU2(:)
      real(8), intent(inout) :: VV2(:)

      integer :: I, I1, I2

      ! Update values of secondary node !
      do I = 1, NNPERBC
         I1 = IPERCONN(I, 1)
         I2 = IPERCONN(I, 2)

         UU2(I2) = UU2(I1)
         VV2(I2) = VV2(I1)
      end do
   end subroutine UPDATE_U_PERSLNODES

end module mod_momentum_bc_forcing
