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
C*******************************************************************
C
C  MODULE VEW1D
C
C  This module declares subroutines for vertical element walls and 
C  1D condensation
C
C  2022-09-21: Created by SB
C
C*******************************************************************
      module vew1d
      implicit none

      contains

      !-------------------------------------------------------------
      SUBROUTINE GET_STREAMDIRECTION_AT_ELEMENT
     &           (NM1,NM2,NM3,IA1,IA2,IB1,IB2,
     &            CAX,CAY,CBX,CBY,SX,SY,STATUS)
      !-------------------------------------------------------------
      USE MESH,ONLY : X,Y,findIBTYPEAtNode
      USE NodalAttributes,ONLY : NCondensedNodes,CondensedNodes

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: NM1,NM2,NM3
      INTEGER,INTENT(OUT) :: IA1,IA2,IB1,IB2
      REAL(8),INTENT(OUT) :: CAX,CAY,CBX,CBY,SX,SY
      INTEGER,INTENT(OUT) :: STATUS
      INTEGER :: IBTYPE_A1,IBTYPE_A2,IBTYPE_B1
      REAL(8) :: LEN
      LOGICAL :: ON_BOUNDARY

      STATUS = 0

      if ((NCondensedNodes(nm1) < 2.or.  ! if any of the nodes are not condensed nodes
     &     NCondensedNodes(nm2) < 2.or.
     &     NCondensedNodes(nm3) < 2).or. ! or
     &    (NCondensedNodes(nm1) > 2.and. ! if this is a joint element
     &     NCondensedNodes(nm2) > 2.and.
     &     NCondensedNodes(nm3) > 2)) then
         status = -1
         return
      end if

      IF(CondensedNodePairFound(NM1,NM2)) THEN
         IA1 = NM1
         IA2 = NM2
         IB1 = NM3 
      ELSE IF(CondensedNodePairFound(NM2,NM1)) THEN
         IA1 = NM2
         IA2 = NM1
         IB1 = NM3 
      ELSE IF(CondensedNodePairFound(NM2,NM3)) THEN
         IA1 = NM2
         IA2 = NM3
         IB1 = NM1
      ELSE IF(CondensedNodePairFound(NM3,NM2)) THEN
         IA1 = NM3
         IA2 = NM2
         IB1 = NM1
      ELSE IF(CondensedNodePairFound(NM3,NM1)) THEN
         IA1 = NM3
         IA2 = NM1
         IB1 = NM2
      ELSE IF(CondensedNodePairFound(NM1,NM3)) THEN
         IA1 = NM1
         IA2 = NM3
         IB1 = NM2
      ELSE
         status = -1
         return !when condensed node pair not found
      ENDIF
      IB2 = FIND_IB2(IB1,IA1,IA2)
      IF(IB2 < 1) THEN
         status = -1
         return !when condensed node pair not found
      ENDIF
      IF(CondensedNodes(IB2,1).EQ.0) THEN
         status = -1
         return !when condensed node pair not found
      ENDIF

      ibtype_a1 = findIBTYPEAtNode(ia1)
      ibtype_a2 = findIBTYPEAtNode(ia2)
      ibtype_b1 = findIBTYPEAtNode(ib1)

      on_boundary = .false.
      if ((ibtype_a1.eq.20.and.ibtype_b1.eq.20).or.
     &    (ibtype_a2.eq.20.and.ibtype_b1.eq.20).or.
     &    (ibtype_a1.eq.64.and.ibtype_b1.eq.64).or.
     &    (ibtype_a2.eq.64.and.ibtype_b1.eq.64)) then
         on_boundary = .true.
      end if

      if (.not.on_boundary) then
         status = -1
         return
      end if

      CAX = 0.5*(X(IA1)+X(IA2))
      CAY = 0.5*(Y(IA1)+Y(IA2)) 
      CBX = 0.5*(X(IB1)+X(IB2))
      CBY = 0.5*(Y(IB1)+Y(IB2)) 
      SX = CAX - CBX
      SY = CAY - CBY
      LEN = SQRT(SX*SX+SY*SY)
      SX = SX/LEN
      SY = SY/LEN

      RETURN
      END SUBROUTINE GET_STREAMDIRECTION_AT_ELEMENT

      !-------------------------------------------------------------
      SUBROUTINE ROTATE_AT_CONDENSEDNODES(VECX,VECY,IND,SX,SY)
      !-------------------------------------------------------------
      USE GLOBAL,ONLY : NODECODE
      USE NodalAttributes,ONLY : LoadCondensedNodes,NCondensedNodes
      USE Mesh,ONLY : LBArray_Pointer
      USE BOUNDARIES,ONLY : LBCODEI,NBV,IBCONN,ISSUBMERGED64
      IMPLICIT NONE

      REAL(8),INTENT(INOUT) :: VECX, VECY
      REAL(8),INTENT(IN) :: SX,SY
      INTEGER,INTENT(IN) :: IND
      REAL(8) :: CS,SI,CSCS,CSSI,SISI,fBUF1,fBUF2,VLEN,P
      INTEGER :: IBND,NNBB1,NNBB2
      LOGICAL :: DO_ROTATE
      
      IF(.NOT.LoadCondensedNodes) RETURN
      IF(.NOT.NCondensedNodes(IND).NE.2) RETURN
      
      IBND = LBArray_Pointer(IND)
      IF (IBND.LE.0) RETURN

      ! Rotate the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls 
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      IF(.NOT.DO_ROTATE) RETURN

      CS = SX
      SI = SY

      VLEN = SQRT(VECX*VECX + VECY*VECY)

      P = VECX*CS + VECY*SI

      VECX = VLEN*CS*SIGN(1.D0,P)
      VECY = VLEN*SI*SIGN(1.D0,P)

      RETURN
      END SUBROUTINE ROTATE_AT_CONDENSEDNODES

      !-------------------------------------------------------------
      SUBROUTINE ROTATEBACK_AT_CONDENSEDNODES(VECX,VECY,IND)
      !-------------------------------------------------------------
      USE GLOBAL,ONLY : NODECODE
      USE NodalAttributes,ONLY : LoadCondensedNodes,NCondensedNodes
      USE Mesh,ONLY : LBArray_Pointer
      USE BOUNDARIES,ONLY : CSIICN,SIIICN,LBCODEI,NBV,IBCONN,
     &                      ISSUBMERGED64
      IMPLICIT NONE

      REAL(8),INTENT(INOUT) :: VECX, VECY
      INTEGER,INTENT(IN) :: IND
      REAL(8) CS,SI,CSCS,CSSI,SISI,fBUF1,fBUF2,VLEN,P
      INTEGER :: IBND,NNBB1,NNBB2
      LOGICAL :: DO_ROTATE
      
      IF(.NOT.LoadCondensedNodes) RETURN
      IF(NCondensedNodes(IND).NE.2) RETURN

      IBND = LBArray_Pointer(IND)
      IF (IBND.LE.0) RETURN

      ! Rotate back the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      IF(.NOT.DO_ROTATE) RETURN

      VLEN = SQRT(VECX*VECX+VECY*VECY)
      CS = -SIIICN(IBND)
      SI = CSIICN(IBND)
      P = VECX*CS + VECY*SI
      VECX = VLEN*CS*SIGN(1.D0,P)
      VECY = VLEN*SI*SIGN(1.D0,P)

      RETURN
      END SUBROUTINE ROTATEBACK_AT_CONDENSEDNODES
C
      !-------------------------------------------------------------
      SUBROUTINE ROTATE_AT_CONDENSEDNODES_ALL
     &           (NM1,NM2,NM3,
     &            U1N1,V1N1,U1N2,V1N2,U1N3,V1N3,
     &            QX1N1,QY1N1,QX1N2,QY1N2,QX1N3,QY1N3,
     &            STATUS)
      !-------------------------------------------------------------
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: NM1,NM2,NM3
      REAL(8),INTENT(INOUT) :: 
     &     U1N1,V1N1,U1N2,V1N2,U1N3,V1N3,
     &     QX1N1,QY1N1,QX1N2,QY1N2,QX1N3,QY1N3
      INTEGER,INTENT(OUT) :: STATUS
      REAL(8) :: CAX,CAY,CBX,CBY,SX,SY
      INTEGER :: IA1,IA2,IB1,IB2

      CALL GET_STREAMDIRECTION_AT_ELEMENT
     &     (NM1,NM2,NM3,IA1,IA2,IB1,IB2,
     &      CAX,CAY,CBX,CBY,SX,SY,STATUS)

      IF (STATUS.NE.0) RETURN

      CALL ROTATE_AT_CONDENSEDNODES(U1N1,V1N1,NM1,SX,SY)
      CALL ROTATE_AT_CONDENSEDNODES(QX1N1,QY1N1,NM1,SX,SY)

      CALL ROTATE_AT_CONDENSEDNODES(U1N2,V1N2,NM2,SX,SY)
      CALL ROTATE_AT_CONDENSEDNODES(QX1N2,QY1N2,NM2,SX,SY)

      CALL ROTATE_AT_CONDENSEDNODES(U1N3,V1N3,NM3,SX,SY)
      CALL ROTATE_AT_CONDENSEDNODES(QX1N3,QY1N3,NM3,SX,SY)

      RETURN
      END SUBROUTINE ROTATE_AT_CONDENSEDNODES_ALL
C
      !-------------------------------------------------------------
      SUBROUTINE ROTATEBACK_AT_CONDENSEDNODES_ALL
     &           (NM1,NM2,NM3,TEMP_LV_A1,TEMP_LV_B1,
     &            TEMP_LV_A2,TEMP_LV_B2,TEMP_LV_A3,TEMP_LV_B3,
     &            STATUS)
      !-------------------------------------------------------------
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: NM1,NM2,NM3
      REAL(8),INTENT(INOUT) :: 
     &     TEMP_LV_A1,TEMP_LV_B1,TEMP_LV_A2,TEMP_LV_B2,
     &     TEMP_LV_A3,TEMP_LV_B3
      INTEGER,INTENT(OUT) :: STATUS
      REAL(8) :: CAX,CAY,CBX,CBY,SX,SY
      INTEGER :: IA1,IA2,IB1,IB2

      CALL GET_STREAMDIRECTION_AT_ELEMENT
     &     (NM1,NM2,NM3,IA1,IA2,IB1,IB2,
     &      CAX,CAY,CBX,CBY,SX,SY,STATUS)

      IF (STATUS.NE.0) RETURN

      CALL ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A1,TEMP_LV_B1,NM1)
      CALL ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A2,TEMP_LV_B2,NM2)
      CALL ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A3,TEMP_LV_B3,NM3)

      RETURN
      END SUBROUTINE ROTATEBACK_AT_CONDENSEDNODES_ALL
C
      !-------------------------------------------------------------
      SUBROUTINE REMOVE_NORMAL_AT_CONDENSEDNODES
     &           (VECX,VECY,IND,STATUS)
      !-------------------------------------------------------------
      USE GLOBAL,ONLY : NODECODE
      USE NodalAttributes,ONLY : LoadCondensedNodes,NCondensedNodes
      USE Mesh,ONLY : LBArray_Pointer
      USE BOUNDARIES,ONLY : CSIICN,SIIICN,LBCODEI,NBV,IBCONN,
     &                      ISSUBMERGED64
      
      IMPLICIT NONE

      REAL(8),INTENT(INOUT) :: VECX, VECY
      INTEGER,INTENT(IN) :: IND
      INTEGER,INTENT(OUT) :: STATUS
      REAL(8) CS,SI,CSCS,CSSI,SISI,fBUF1,fBUF2,LEN
      INTEGER :: IBND,NNBB1,NNBB2
      LOGICAL :: DO_ROTATE
      
      STATUS = -1

      IF(.NOT.LoadCondensedNodes) RETURN
      IF(NCondensedNodes(IND).NE.2) RETURN

      IBND = LBArray_Pointer(IND)
      IF (IBND.LE.0) RETURN

      ! Remove the normal component if there is no water on the 
      ! bank, i.e., along land boundaries, or along not-submerged
      ! vertical element walls 
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      IF(.NOT.DO_ROTATE) RETURN  

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

      RETURN
      END SUBROUTINE REMOVE_NORMAL_AT_CONDENSEDNODES

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
      FUNCTION ROTATE_OR_NOT(IBND)
      !-------------------------------------------------------------
      USE GLOBAL,ONLY : NODECODE
      USE BOUNDARIES,ONLY : IBCONN,LBCODEI,NBV,ISSUBMERGED64
      
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IBND
      INTEGER :: NNBB1,NNBB2
      LOGICAL :: ROTATE_OR_NOT

      ROTATE_OR_NOT = .false.
      IF(LBCODEI(IBND).EQ.64) THEN
         NNBB1=NBV(IBND)
         NNBB2=IBCONN(IBND)
         IF ((ISSUBMERGED64(IBND).EQ.0).OR.
     &       (NODECODE(NNBB1).EQ.0.OR.NODECODE(NNBB2).EQ.0)) THEN
            ROTATE_OR_NOT = .true.
         ENDIF
      ELSE IF(LBCODEI(IBND).EQ.0.OR.LBCODEI(IBND).EQ.1.OR.
     &        LBCODEI(IBND).EQ.10.OR.LBCODEI(IBND).EQ.1.OR.
     &        LBCODEI(IBND).EQ.20.OR.LBCODEI(IBND).EQ.21) THEN
         ROTATE_OR_NOT = .true.
      END IF
      END FUNCTION

      FUNCTION CondensedNodePairFound(N1,N2)
      USE NodalAttributes, ONLY : NCondensedNodes,CondensedNodes,
     &                            CondensedNodesNoOfVals
      IMPLICIT NONE
      INTEGER :: N1,N2
      LOGICAL :: CondensedNodePairFound
      INTEGER :: I,N

      CondensedNodePairFound = .false.
      IF (CondensedNodes(N1,1) < 0) THEN
         N = ABS(CondensedNodes(N1,1))
      ELSE
         N = N1
      END IF
      DO I=1,CondensedNodesNoOfVals
         IF (CondensedNodes(N,I).EQ.N2) THEN
            CondensedNodePairFound = .true.
            EXIT
         END IF
      END DO

      END FUNCTION CondensedNodePairFound

      FUNCTION FIND_IB2(N1,N2,N3)
      USE MESH,ONLY : NNEIGHELE,NEITABELE,NM
      IMPLICIT NONE
      INTEGER :: FIND_IB2,N1,N2,N3
      INTEGER :: nnm1,nnm2,nnm3,tmp
      INTEGER :: I,IE
      LOGICAL :: FOUND

      find_ib2 = -1
      found = .false.

      do i=1,nneighele(n1)
         ie = neitabele(n1,i)
         if (ie.le.0) cycle ! sb251003: Added this check to avoid out-of-bounds access.
         nnm1 = nm(ie,1)
         nnm2 = nm(ie,2)
         nnm3 = nm(ie,3)
         if (nnm1.eq.n1) then
            continue ! do nothing
         else if (nnm2.eq.n1) then
            tmp = nnm1
            nnm1 = nnm2
            nnm2 = nnm3
            nnm3 = tmp
         else if (nnm3.eq.n1) then
            tmp = nnm1
            nnm1 = nnm3
            nnm3 = nnm2
            nnm2 = tmp
         else
            STOP('Error in FIND_IB2: node not found in neighbor element')
         end if
         if ((nnm2.eq.n2.and.nnm3.eq.n3).or.
     &       (nnm2.eq.n3.and.nnm3.eq.n2)) cycle
         if (nnm2.eq.n2.or.nnm2.eq.n3) then
            find_ib2 = nnm3
            found = .true.
         else if (nnm3.eq.n2.or.nnm3.eq.n3) then
            find_ib2 = nnm2
            found = .true.
         end if
         if (found) exit
      end do
      END FUNCTION FIND_IB2

      end module vew1d
