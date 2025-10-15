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
             REMOVE_NORMAL_AT_CONDENSEDNODES, GET_STREAMDIRECTION_AT_ELEMENT

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

end module vew1d
