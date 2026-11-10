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
!*******************************************************************
!
!  MODULE VEW1D
!
!  This module declares subroutines for vertical element walls and
!  1D condensation
!
!  2022-09-21: Created by SB
!
!*******************************************************************
module vew1d
   implicit none

   private

   public :: ROTATE_AT_CONDENSEDNODES_ALL, ROTATEBACK_AT_CONDENSEDNODES_ALL, &
             REMOVE_NORMAL_AT_CONDENSEDNODES, GET_STREAMDIRECTION_AT_ELEMENT, &
             COMPUTE_VEW_ADVECT_VELOCITY, UADV, VADV, &
             VEW_ADVECT_FR_ATTACH, VEW_ADVECT_FR_FULL

   !...Recovered crest velocity at coupled (submerged) IBTYPE=64 VEW node pairs, on with
   !   VEW_CREST_HYDRAULICS (WEIR module, VEW1DChannelControl namelist). A coupled pair shares one
   !   velocity, so the shallower (wall-top) side carries less than the discharge reaching the wall.
   !   The momentum and GWCE advective terms then use, at the nodes of each coupled pair, a velocity
   !   whose wall-normal component (boundary normal CSII/SIII of the pair) is
   !      v'n = q / min(H1, H2),
   !   q the wall-normal discharge averaged over the wet neighbours of the pair's foot (deeper) node
   !   that are not on an IBTYPE=64 boundary; the tangential component is kept as is. Where the foot
   !   has no such neighbour, or both nodes have the same depth, v'n = vn*0.5*(H1 + H2)/min(H1, H2).
   !   UADV/VADV hold the velocities those terms use. Further:
   !   - Lee loss: the foot node of a pair whose bed side is downstream (flow leaving the wall, judged
   !     from the neighbours' flow) keeps its own velocity. The approach side and the whole crest carry
   !     v', the lee-side bed elements see the v' -> v drop, so the velocity head over the crest is lost
   !     as in the broad-crested weir formulas (expansion loss) instead of recovered.
   !   - Froude taper: with Fr = |v'n|/sqrt(g*min(H)), the recovery has weight 0 (flow attached to the
   !     wall, no recovery and no lee loss) below VEW_ADVECT_FR_ATTACH, 1 above VEW_ADVECT_FR_FULL,
   !     linear between.
   !   - Froude cap: |v'n| <= sqrt(g*min(H1, H2)), so the recovered flow stays subcritical.
   !   - Condensed wall tops: both faces of a wall share one discharge, the one of the face
   !     whose foot is on the upstream side, stored per condensed group (QGROUP, keyed by the group's
   !     first node).
   real(8), allocatable, target :: UADV(:), VADV(:)
   real(8) :: VEW_ADVECT_FR_ATTACH = 0.05d0
   real(8) :: VEW_ADVECT_FR_FULL = 0.1d0
   real(8), allocatable :: QGROUP(:)
   logical, allocatable :: QGROUPSET(:)
   logical, allocatable :: ISVEWNODE(:)

contains

   !-------------------------------------------------------------
   subroutine GET_STREAMDIRECTION_AT_ELEMENT &
      (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
       CAX, CAY, CBX, CBY, SX, SY, STATUS)
      use MESH, only: X, Y, findIBTYPEAtNode
      use NodalAttributes, only: NCondensedNodes, CondensedNodes

      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      integer, intent(OUT) :: IA1, IA2, IB1, IB2
      real(8), intent(OUT) :: CAX, CAY, CBX, CBY, SX, SY
      integer, intent(OUT) :: STATUS
      integer :: IBTYPE_A1, IBTYPE_A2, IBTYPE_B1
      real(8) :: LEN
      logical :: ON_BOUNDARY

      STATUS = 0

      if ((NCondensedNodes(nm1) < 2 .or. & ! if any of the nodes are not condensed nodes
           NCondensedNodes(nm2) < 2 .or. &
           NCondensedNodes(nm3) < 2) .or. & ! or
          (NCondensedNodes(nm1) > 2 .and. & ! if this is a joint element
           NCondensedNodes(nm2) > 2 .and. &
           NCondensedNodes(nm3) > 2)) then
         status = -1
         return
      end if

      if (CondensedNodePairFound(NM1, NM2)) then
         IA1 = NM1
         IA2 = NM2
         IB1 = NM3
      else if (CondensedNodePairFound(NM2, NM1)) then
         IA1 = NM2
         IA2 = NM1
         IB1 = NM3
      else if (CondensedNodePairFound(NM2, NM3)) then
         IA1 = NM2
         IA2 = NM3
         IB1 = NM1
      else if (CondensedNodePairFound(NM3, NM2)) then
         IA1 = NM3
         IA2 = NM2
         IB1 = NM1
      else if (CondensedNodePairFound(NM3, NM1)) then
         IA1 = NM3
         IA2 = NM1
         IB1 = NM2
      else if (CondensedNodePairFound(NM1, NM3)) then
         IA1 = NM1
         IA2 = NM3
         IB1 = NM2
      else
         status = -1
         return !when condensed node pair not found
      end if
      IB2 = FIND_IB2(IB1, IA1, IA2)
      if (IB2 < 1) then
         status = -1
         return !when condensed node pair not found
      end if
      if (CondensedNodes(IB2, 1) == 0) then
         status = -1
         return !when condensed node pair not found
      end if

      ibtype_a1 = findIBTYPEAtNode(ia1)
      ibtype_a2 = findIBTYPEAtNode(ia2)
      ibtype_b1 = findIBTYPEAtNode(ib1)

      on_boundary = .false.
      if ((ibtype_a1 == 20 .and. ibtype_b1 == 20) .or. &
          (ibtype_a2 == 20 .and. ibtype_b1 == 20) .or. &
          (ibtype_a1 == 64 .and. ibtype_b1 == 64) .or. &
          (ibtype_a2 == 64 .and. ibtype_b1 == 64)) then
         on_boundary = .true.
      end if

      if (.not. on_boundary) then
         status = -1
         return
      end if

      CAX = 0.5d0*(X(IA1) + X(IA2))
      CAY = 0.5d0*(Y(IA1) + Y(IA2))
      CBX = 0.5d0*(X(IB1) + X(IB2))
      CBY = 0.5d0*(Y(IB1) + Y(IB2))
      SX = CAX - CBX
      SY = CAY - CBY
      LEN = sqrt(SX*SX + SY*SY)
      SX = SX/LEN
      SY = SY/LEN

      return
   end subroutine GET_STREAMDIRECTION_AT_ELEMENT

   !-------------------------------------------------------------
   subroutine ROTATE_AT_CONDENSEDNODES(VECX, VECY, IND, SX, SY)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      real(8), intent(IN) :: SX, SY
      integer, intent(IN) :: IND
      real(8) :: CS, SI, VLEN, P
      integer :: IBND
      logical :: DO_ROTATE

      if (.not. LoadCondensedNodes) return
      if (.not. NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Rotate the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      CS = SX
      SI = SY

      VLEN = sqrt(VECX*VECX + VECY*VECY)

      P = VECX*CS + VECY*SI

      VECX = VLEN*CS*sign(1.d0, P)
      VECY = VLEN*SI*sign(1.d0, P)

   end subroutine ROTATE_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   subroutine ROTATEBACK_AT_CONDENSEDNODES(VECX, VECY, IND)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      use BOUNDARIES, only: CSIICN, SIIICN
      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      integer, intent(IN) :: IND
      real(8) :: CS, SI, VLEN, P
      integer :: IBND
      logical :: DO_ROTATE

      if (.not. LoadCondensedNodes) return
      if (NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Rotate back the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      VLEN = sqrt(VECX*VECX + VECY*VECY)
      CS = -SIIICN(IBND)
      SI = CSIICN(IBND)
      P = VECX*CS + VECY*SI
      VECX = VLEN*CS*sign(1.d0, P)
      VECY = VLEN*SI*sign(1.d0, P)

      return
   end subroutine ROTATEBACK_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   subroutine ROTATE_AT_CONDENSEDNODES_ALL &
      (NM1, NM2, NM3, &
       U1N1, V1N1, U1N2, V1N2, U1N3, V1N3, &
       QX1N1, QY1N1, QX1N2, QY1N2, QX1N3, QY1N3, &
       STATUS)
      !-------------------------------------------------------------
      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      real(8), intent(INOUT) :: &
         U1N1, V1N1, U1N2, V1N2, U1N3, V1N3, &
         QX1N1, QY1N1, QX1N2, QY1N2, QX1N3, QY1N3
      integer, intent(OUT) :: STATUS
      real(8) :: CAX, CAY, CBX, CBY, SX, SY
      integer :: IA1, IA2, IB1, IB2

      call GET_STREAMDIRECTION_AT_ELEMENT &
         (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
          CAX, CAY, CBX, CBY, SX, SY, STATUS)

      if (STATUS /= 0) return

      call ROTATE_AT_CONDENSEDNODES(U1N1, V1N1, NM1, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N1, QY1N1, NM1, SX, SY)

      call ROTATE_AT_CONDENSEDNODES(U1N2, V1N2, NM2, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N2, QY1N2, NM2, SX, SY)

      call ROTATE_AT_CONDENSEDNODES(U1N3, V1N3, NM3, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N3, QY1N3, NM3, SX, SY)

      return
   end subroutine ROTATE_AT_CONDENSEDNODES_ALL

   !-------------------------------------------------------------
   subroutine ROTATEBACK_AT_CONDENSEDNODES_ALL &
      (NM1, NM2, NM3, TEMP_LV_A1, TEMP_LV_B1, &
       TEMP_LV_A2, TEMP_LV_B2, TEMP_LV_A3, TEMP_LV_B3, &
       STATUS)
      !-------------------------------------------------------------
      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      real(8), intent(INOUT) :: &
         TEMP_LV_A1, TEMP_LV_B1, TEMP_LV_A2, TEMP_LV_B2, &
         TEMP_LV_A3, TEMP_LV_B3
      integer, intent(OUT) :: STATUS
      real(8) :: CAX, CAY, CBX, CBY, SX, SY
      integer :: IA1, IA2, IB1, IB2

      call GET_STREAMDIRECTION_AT_ELEMENT &
         (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
          CAX, CAY, CBX, CBY, SX, SY, STATUS)

      if (STATUS /= 0) return

      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A1, TEMP_LV_B1, NM1)
      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A2, TEMP_LV_B2, NM2)
      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A3, TEMP_LV_B3, NM3)

      return
   end subroutine ROTATEBACK_AT_CONDENSEDNODES_ALL

   !-------------------------------------------------------------
   subroutine REMOVE_NORMAL_AT_CONDENSEDNODES(VECX, VECY, IND, STATUS)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      use BOUNDARIES, only: CSIICN, SIIICN

      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      integer, intent(IN) :: IND
      integer, intent(OUT) :: STATUS
      real(8) :: CS, SI, CSCS, CSSI, SISI, fBUF1, fBUF2
      integer :: IBND
      logical :: DO_ROTATE

      STATUS = -1

      if (.not. LoadCondensedNodes) return
      if (NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Remove the normal component if there is no water on the
      ! bank, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      CS = -SIIICN(IBND)
      SI = CSIICN(IBND)
      CSCS = CS*CS
      CSSI = CS*SI
      SISI = SI*SI

      fBUF1 = VECX*CSCS + VECY*CSSI
      fBUF2 = VECX*CSSI + VECY*SISI

      VECX = fBUF1
      VECY = fBUF2

      STATUS = 0

      return
   end subroutine REMOVE_NORMAL_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   ! FUNCTION ROTATE_OR_NOT
   !
   ! Judges if velocity solutions should be rotated at a node
   ! along a boundary. They should be rotated if 1) a node is
   ! along a IBTYPE=64 vertical element wall boundary and the
   ! boundary is not submerged, or if 2) a node is along any
   ! zero normal flow boundary. This function is applied to a
   ! node where the condensed_nodes nodal attribute is given.
   !-------------------------------------------------------------
   function ROTATE_OR_NOT(IBND)
      use GLOBAL, only: NODECODE
      use BOUNDARIES, only: IBCONN, LBCODEI, NBV, ISSUBMERGED64

      implicit none

      integer, intent(IN) :: IBND
      integer :: NNBB1, NNBB2
      logical :: ROTATE_OR_NOT

      ROTATE_OR_NOT = .false.
      if (LBCODEI(IBND) == 64) then
         NNBB1 = NBV(IBND)
         NNBB2 = IBCONN(IBND)
         if ((ISSUBMERGED64(IBND) == 0) .or. &
             (NODECODE(NNBB1) == 0 .or. NODECODE(NNBB2) == 0)) then
            ROTATE_OR_NOT = .true.
         end if
      else if (LBCODEI(IBND) == 0 .or. LBCODEI(IBND) == 1 .or. &
               LBCODEI(IBND) == 10 .or. LBCODEI(IBND) == 1 .or. &
               LBCODEI(IBND) == 20 .or. LBCODEI(IBND) == 21) then
         ROTATE_OR_NOT = .true.
      end if
   end function ROTATE_OR_NOT

   function CondensedNodePairFound(N1, N2)
      use NodalAttributes, only: CondensedNodes, CondensedNodesNoOfVals
      implicit none
      integer, intent(in) :: N1, N2
      logical :: CondensedNodePairFound
      integer :: I, N

      CondensedNodePairFound = .false.
      if (CondensedNodes(N1, 1) < 0) then
         N = abs(CondensedNodes(N1, 1))
      else
         N = N1
      end if
      do I = 1, CondensedNodesNoOfVals
         if (CondensedNodes(N, I) == N2) then
            CondensedNodePairFound = .true.
            exit
         end if
      end do

   end function CondensedNodePairFound

   function FIND_IB2(N1, N2, N3)
      use MESH, only: NNEIGHELE, NEITABELE, NM
      implicit none
      integer, intent(IN) :: N1, N2, N3
      integer :: FIND_IB2
      integer :: nnm1, nnm2, nnm3, tmp
      integer :: I, IE
      logical :: FOUND

      find_ib2 = -1
      found = .false.

      do i = 1, nneighele(n1)
         ie = neitabele(n1, i)
         if (ie <= 0) cycle ! sb251003: Added this check to avoid out-of-bounds access.
         nnm1 = nm(ie, 1)
         nnm2 = nm(ie, 2)
         nnm3 = nm(ie, 3)
         if (nnm1 == n1) then
            continue ! do nothing
         else if (nnm2 == n1) then
            tmp = nnm1
            nnm1 = nnm2
            nnm2 = nnm3
            nnm3 = tmp
         else if (nnm3 == n1) then
            tmp = nnm1
            nnm1 = nnm3
            nnm3 = nnm2
            nnm2 = tmp
         else
            stop 'Error in FIND_IB2: node not found in neighbor element'
         end if
         if ((nnm2 == n2 .and. nnm3 == n3) .or. &
             (nnm2 == n3 .and. nnm3 == n2)) cycle
         if (nnm2 == n2 .or. nnm2 == n3) then
            find_ib2 = nnm3
            found = .true.
         else if (nnm3 == n2 .or. nnm3 == n3) then
            find_ib2 = nnm2
            found = .true.
         end if
         if (found) exit
      end do
   end function FIND_IB2

   !-------------------------------------------------------------
   !  COMPUTE_VEW_ADVECT_VELOCITY
   !  Fill UADV/VADV with UU1/VV1 and, if VEW_CREST_HYDRAULICS is set, overwrite them
   !  at both nodes of each coupled IBTYPE=64 pair (both nodes wet, min depth > H0)
   !  with the recovered velocity above. Called once per time step before the GWCE.
   !-------------------------------------------------------------
   subroutine COMPUTE_VEW_ADVECT_VELOCITY()
      use GLOBAL, only: UU1, VV1, ETA2, NODECODE, IFNLFA, H0
      use ADC_CONSTANTS, only: G
      use MESH, only: NP, DP, X, Y, NNeigh, NeiTab
      use BOUNDARIES, only: NBOU, NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, CSII, SIII
      use NodalAttributes, only: LoadCondensedNodes, CondensedNodes
      use WEIR, only: VEW_CREST_HYDRAULICS
      implicit none
      integer :: I, J, K, N1, N2, NF, JJ, NT, KG
      real(8) :: HA, HB, HMIN, FAC, UAVG, VAVG, UN, U0, V0, DX, DY, QN, VNP, VCRIT
      real(8) :: UNB, VNB, WREC, FRC
      integer :: NB, CNT, II
      logical :: LEEFOOT

      if (.not. allocated(UADV)) allocate (UADV(NP), VADV(NP))
      UADV(:) = UU1(:)
      VADV(:) = VV1(:)
      if (.not. VEW_CREST_HYDRAULICS) return

      if (.not. allocated(ISVEWNODE)) then
         allocate (ISVEWNODE(NP))
         ISVEWNODE(:) = .false.
         II = 0
         do K = 1, NBOU
            select case (LBCODEI(II + 1))
            case (64)
               do J = 1, NVELL(K)
                  II = II + 1
                  ISVEWNODE(NBV(II)) = .true.
                  ISVEWNODE(IBCONN(II)) = .true.
               end do
               II = II + NVELL(K)
            case (4, 24, 5, 25)
               II = II + NVELL(K)*2
            case default
               II = II + NVELL(K)
            end select
         end do
      end if

      !...Pass 1 (condensed wall tops): discharge of the upstream face of each condensed wall-top group
      if (LoadCondensedNodes) then
         if (.not. allocated(QGROUP)) allocate (QGROUP(NP), QGROUPSET(NP))
         QGROUPSET(:) = .false.
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (64)
               do J = 1, NVELL(K)
                  I = I + 1
                  if (ISSUBMERGED64(I) == 0) cycle
                  N1 = NBV(I)
                  N2 = IBCONN(I)
                  if (DP(N1) == DP(N2)) cycle
                  if (NODECODE(N1) == 0 .or. NODECODE(N2) == 0) cycle
                  NF = merge(N2, N1, DP(N2) > DP(N1))
                  NT = merge(N1, N2, NF == N2)
                  KG = CONDENSED_GROUP(NT)
                  if (KG == 0) cycle
                  call NEIGHBOR_FLOW(NF, I, QN, UNB, VNB, DX, DY, CNT)
                  if (CNT == 0) cycle
                  !...upstream face: the flow goes from the foot's neighbours towards the wall
                  if (UNB*DX + VNB*DY < 0.d0) then
                     QGROUP(KG) = QN
                     QGROUPSET(KG) = .true.
                  end if
               end do
               I = I + NVELL(K)
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case default
               I = I + NVELL(K)
            end select
         end do
      end if

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               if (ISSUBMERGED64(I) == 0) cycle
               N1 = NBV(I)
               N2 = IBCONN(I)
               if (NODECODE(N1) == 0 .or. NODECODE(N2) == 0) cycle
               HA = DP(N1) + IFNLFA*ETA2(N1)
               HB = DP(N2) + IFNLFA*ETA2(N2)
               HMIN = min(HA, HB)
               if (HMIN <= H0) cycle
               FAC = 0.5d0*(HA + HB)/HMIN
               UAVG = 0.5d0*(UU1(N1) + UU1(N2))
               VAVG = 0.5d0*(VV1(N1) + VV1(N2))
               U0 = UAVG
               V0 = VAVG
               !...Lee-side foot: the deeper node of the pair, if the flow goes from it into the
               !   bed-side domain. Its neighbours all lie on the bed side (the coincident wall-top
               !   node is not connected to it), so their mean position minus its own points into it.
               LEEFOOT = .false.
               if (DP(N1) /= DP(N2)) then
                  NF = N1
                  if (DP(N2) > DP(N1)) NF = N2
                  DX = 0.d0
                  DY = 0.d0
                  do JJ = 2, NNeigh(NF)
                     DX = DX + X(NeiTab(NF, JJ)) - X(NF)
                     DY = DY + Y(NeiTab(NF, JJ)) - Y(NF)
                  end do
                  LEEFOOT = (U0*DX + V0*DY) > 0.d0
               end if
               !...scale only the normal component: u' = u + (FAC - 1)*(u.n)n
               UN = UAVG*CSII(I) + VAVG*SIII(I)
               UAVG = UAVG + (FAC - 1.d0)*UN*CSII(I)
               VAVG = VAVG + (FAC - 1.d0)*UN*SIII(I)
               if (DP(N1) /= DP(N2)) then
                  NF = N1
                  if (DP(N2) > DP(N1)) NF = N2
                  call NEIGHBOR_FLOW(NF, I, QN, UNB, VNB, DX, DY, CNT)
                  if (CNT > 0) then
                     !...one discharge per wall with condensed tops: the upstream face's (magnitude), with
                     !   the sign of this face's normal flux
                     if (LoadCondensedNodes) then
                        KG = CONDENSED_GROUP(merge(N1, N2, NF == N2))
                        if (KG > 0) then
                           if (QGROUPSET(KG)) QN = sign(abs(QGROUP(KG)), QN)
                        end if
                     end if
                     VNP = QN/HMIN
                     UN = U0*CSII(I) + V0*SIII(I)
                     UAVG = U0 + (VNP - UN)*CSII(I)
                     VAVG = V0 + (VNP - UN)*SIII(I)
                     !...lee side from the neighbours' flow, not the pair's own velocity
                     LEEFOOT = (UNB*DX + VNB*DY) > 0.d0
                     !...taper with the crest Froude number
                     FRC = abs(VNP)/sqrt(G*HMIN)
                     if (FRC <= VEW_ADVECT_FR_ATTACH) then
                        WREC = 0.d0
                     elseif (FRC >= VEW_ADVECT_FR_FULL) then
                        WREC = 1.d0
                     else
                        WREC = (FRC - VEW_ADVECT_FR_ATTACH)/(VEW_ADVECT_FR_FULL - VEW_ADVECT_FR_ATTACH)
                     end if
                     if (WREC <= 0.d0) cycle ! attached: plain velocities at both nodes, no lee loss
                     UAVG = U0 + WREC*(UAVG - U0)
                     VAVG = V0 + WREC*(VAVG - V0)
                  end if
               end if
               !...Froude cap
               VCRIT = sqrt(G*HMIN)
               UN = UAVG*CSII(I) + VAVG*SIII(I)
               if (abs(UN) > VCRIT) then
                  VNP = sign(VCRIT, UN)
                  UAVG = UAVG + (VNP - UN)*CSII(I)
                  VAVG = VAVG + (VNP - UN)*SIII(I)
               end if
               if (LEEFOOT) then
                  !...recovered velocity at the wall-top node only; the lee foot keeps its own
                  N1 = merge(N2, N1, NF == N1) ! N1 := the wall-top node
                  UADV(N1) = UAVG
                  VADV(N1) = VAVG
               else
                  UADV(N1) = UAVG
                  UADV(N2) = UAVG
                  VADV(N1) = VAVG
                  VADV(N2) = VAVG
               end if
            end do
            I = I + NVELL(K)
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case default
            I = I + NVELL(K)
         end select
      end do
   end subroutine COMPUTE_VEW_ADVECT_VELOCITY

   !-------------------------------------------------------------
   !  NEIGHBOR_FLOW: mean wall-normal discharge QN (boundary normal of entry IB) and mean velocity
   !  (UNB, VNB) over the wet neighbours of the foot node NF that are not on a VEW (CNT of them),
   !  and the sum (DX, DY) of the offsets of all its neighbours (points from NF into its bed side).
   !-------------------------------------------------------------
   subroutine NEIGHBOR_FLOW(NF, IB, QN, UNB, VNB, DX, DY, CNT)
      use GLOBAL, only: UU1, VV1, ETA2, NODECODE, IFNLFA
      use MESH, only: DP, X, Y, NNeigh, NeiTab
      use BOUNDARIES, only: CSII, SIII
      implicit none
      integer, intent(IN) :: NF, IB
      real(8), intent(OUT) :: QN, UNB, VNB, DX, DY
      integer, intent(OUT) :: CNT
      integer :: JJ, NB

      QN = 0.d0
      UNB = 0.d0
      VNB = 0.d0
      DX = 0.d0
      DY = 0.d0
      CNT = 0
      do JJ = 2, NNeigh(NF)
         NB = NeiTab(NF, JJ)
         DX = DX + X(NB) - X(NF)
         DY = DY + Y(NB) - Y(NF)
         if (NB == NF .or. ISVEWNODE(NB) .or. NODECODE(NB) == 0) cycle
         QN = QN + (UU1(NB)*CSII(IB) + VV1(NB)*SIII(IB))*(DP(NB) + IFNLFA*ETA2(NB))
         UNB = UNB + UU1(NB)
         VNB = VNB + VV1(NB)
         CNT = CNT + 1
      end do
      if (CNT > 0) then
         QN = QN/CNT
         UNB = UNB/CNT
         VNB = VNB/CNT
      end if
   end subroutine NEIGHBOR_FLOW

   !-------------------------------------------------------------
   !  CONDENSED_GROUP: first node of the condensed node group of node N (0 if not condensed)
   !-------------------------------------------------------------

   integer function CONDENSED_GROUP(N)
      use NodalAttributes, only: CondensedNodes
      implicit none
      integer, intent(IN) :: N
      CONDENSED_GROUP = 0
      if (CondensedNodes(N, 1) > 0) then
         CONDENSED_GROUP = N
      elseif (CondensedNodes(N, 1) < 0) then
         CONDENSED_GROUP = -CondensedNodes(N, 1)
      end if
   end function CONDENSED_GROUP

end module vew1d
