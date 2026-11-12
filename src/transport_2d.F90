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
!  MODULE TRANSPORT_2D
!-----------------------------------------------------------------------
!> @brief Depth-averaged (2D) passive scalar transport
!>
!> This module solves the 2D depth-averaged scalar transport equation
!> (C2D_PTrans) for the concentration CH1 using a finite element
!> formulation, with cohesive sediment source/sink terms (Krone
!> deposition, Partheniades erosion).
!-----------------------------------------------------------------------
module mod_transport_2d

   implicit none

   private

   public :: SCALAR_TRANS_2D

contains

!-----------------------------------------------------------------------
!> @brief Compute 2D scalar transport equation
!>
!> Solves the depth-averaged scalar transport equation using a finite
!> element formulation. Computes source/sink terms using classical cohesive
!> sediment transport relations (Krone's deposition, Partheniades' erosion).
!>
!-----------------------------------------------------------------------
   subroutine SCALAR_TRANS_2D()

      use GLOBAL, only: ch1, trans_lv_b, nodecode, noff, lumpt, dtdp, &
                        trans_lv_a, soursin, uu1, vv1, bedstr, tk, H1

      use ADC_CONSTANTS, only: rhowat0
      use MESH, only: NE, NP, NM, X, Y, AREAS, SFAC
      use NodalAttributes, only: EVC
      use mod_logging, only: t_log_scope, init_log_scope

      implicit none

      integer :: IE, I !local loop counters
      integer :: NM1, NM2, NM3
      integer :: NC1, NC2, NC3, NCEle, NCI

      real(8) :: C1, CBEDSTRD, CBEDSTRE, CCRITD
      real(8) :: CH1N1, CH1N2, CH1N3, CHSUM
      real(8) :: DHDX, DHDY
      real(8) :: DXXYY11, DXXYY12, DXXYY13
      real(8) :: DXXYY22, DXXYY23
      real(8) :: DXXYY33
      real(8) :: ECONST
      real(8) :: EVC1, EVC2, EVC3, EVCEA
      real(8) :: FDDDODT, FDDODODT
      real(8) :: HEA
      real(8) :: H1N1, H1N2, H1N3
      real(8) :: HSD, HSE
      real(8) :: SFacAvg
      real(8) :: SS1N1, SS1N2, SS1N3
      real(8) :: TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      real(8) :: TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      real(8) :: U1N1, U1N2, U1N3
      real(8) :: V1N1, V1N2, V1N3
      real(8) :: UV1
      real(8) :: UEA, VEA, UPEA, VPEA
      real(8) :: WS, WSMOD

      real(8) :: AreaIE2
      real(8) :: FDDD, FDDOD
      real(8) :: FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      real(8) :: FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      real(8) :: DDX1, DDX2, DDX3, DDY1, DDY2, DDY3
      real(8) :: DXX11, DXX12, DXX13, DXX22, DXX23, DXX33
      real(8) :: DYY11, DYY12, DYY13, DYY22, DYY23, DYY33

      LOG_SCOPE_TRACED("scalar_trans_2d", TIMESTEP_TRACING)

      !...  NOTE: THE VARIABLE CH1(I) IS ACTUALLY C*H
      !.... COMPUTE SOURCE/SINK TERM AT THE NODES USING CLASSICAL COHESIVE
      !.... SEDIMENT TRANSPORT RELATIONS

      WS = 0.0001d0 ! particle fall velocity [m/s]
      CBEDSTRD = 0.15d0 ! critical shear stress for deposition [N/m^2]
      CCRITD = 0.30d0 ! critical concentration for hindered settling [kg/m^3]
      ECONST = 0.00001d0 ! erosion rate constant [kg/m^2/sec]
      CBEDSTRE = 0.4d0 ! critical shear stress for erosion [N/m^2]

      do I = 1, NP
         UV1 = sqrt(UU1(I)*UU1(I) + VV1(I)*VV1(I))
         BEDSTR = H1(I)*UV1*TK(I)*RhoWat0 ![N/m^2]
         C1 = CH1(I)/H1(I)

         !.....Calculate the deposition rate using Krone's (1962) formulation:
         !.....dC/dt = -P*WSMOD*C/D     where
         !.....WSMOD=WS          when C < Ccrit  and
         !.....WSMOD=K*C**1.33   when C > Ccrit
         !.....D is the average depth through which particles settle D = H/2,
         !.....H is the water depth
         !.....C is the depth-averaged sediment concentration,
         !.....P is the sticking probability  P = (1-BEDSTR/CBEDSTRD),
         !.....CBEDSTRD is the critical bottom stress above which no deposition occurs.
         !.....It was assumed that the constant K could be backed out by setting
         !.....WSMOD = WS when C = Ccrit.

         WSMOD = WS
         if (C1 > CCRITD) WSMOD = WS*(C1/CCRITD)**1.33d0
         HSD = 0.d0
         if (BEDSTR < CBEDSTRD) HSD = -(2.d0*WSMOD*C1)* &
                                      (1.0d0 - BEDSTR/CBEDSTRD)
         if (HSD > 0.d0) HSD = 0.d0

         !.....Calculate the surface erosion rate for cohesive sediment using
         !.....the Ariathurai et at. (1977) adaption of Partheniades' (1962) findings

         HSE = 0.d0
         if (BEDSTR > CBEDSTRE) HSE = ECONST*(BEDSTR/CBEDSTRE - 1.0d0)

         !.....Determine the total source sink term

         SOURSIN(I) = HSD + HSE
      end do

      !.... UPDATE THE TRANSPORT EQUATION ELEMENT BY ELEMENT BY FORMING
      !.... TEMPORARY VECTORS AND THEN ASSEMBLING.  NOTE: TRANS_LV_B(I), TRANS_LV_A(I) ARE
      !.... ZEROED OUT AT THE TOP OF THE TIME STEPPING LOOP.  AGAIN THESE
      !.... LOOPS HAVE BEEN UNROLLED TO OPTIMIZE VECTORIZATION

      do IE = 1, NE

         !.....SET NODAL VALUES FOR EACH ELEMENT

         NM1 = NM(IE, 1)
         NM2 = NM(IE, 2)
         NM3 = NM(IE, 3)
         NC1 = NODECODE(NM1)
         NC2 = NODECODE(NM2)
         NC3 = NODECODE(NM3)
         NCELE = NC1*NC2*NC3*NOFF(IE)
         U1N1 = UU1(NM1)
         U1N2 = UU1(NM2)
         U1N3 = UU1(NM3)
         V1N1 = VV1(NM1)
         V1N2 = VV1(NM2)
         V1N3 = VV1(NM3)
         CH1N1 = CH1(NM1)
         CH1N2 = CH1(NM2)
         CH1N3 = CH1(NM3)
         EVC1 = EVC(NM1)
         EVC2 = EVC(NM2)
         EVC3 = EVC(NM3)
         SS1N1 = SOURSIN(NM1)
         SS1N2 = SOURSIN(NM2)
         SS1N3 = SOURSIN(NM3)
         H1N1 = H1(NM1)
         H1N2 = H1(NM2)
         H1N3 = H1(NM3)
         SFacAvg = (SFAC(NM1) + SFAC(NM2) + SFAC(NM3))/3.d0

         !.....COMPUTE ELEMENTAL MATRICIES

         AREAIE2 = AREAS(IE) !2*element area
         FDX1 = (Y(NM2) - Y(NM3))*SFacAvg !b1
         FDX2 = (Y(NM3) - Y(NM1))*SFacAvg !b2
         FDX3 = (Y(NM1) - Y(NM2))*SFacAvg !b3
         FDY1 = X(NM3) - X(NM2) !a1
         FDY2 = X(NM1) - X(NM3) !a2
         FDY3 = X(NM2) - X(NM1) !a3
         FDX1O2A = FDX1/AREAIE2 !dphi1/dx
         FDY1O2A = FDY1/AREAIE2 !dphi1/dy
         FDX2O2A = FDX2/AREAIE2 !dphi2/dx
         FDY2O2A = FDY2/AREAIE2 !dphi2/dy
         FDX3O2A = FDX3/AREAIE2 !dphi3/dx
         FDY3O2A = FDY3/AREAIE2 !dphi3/dy

         DDX1 = FDX1/3.d0 !<2*(dphi1/dx)*phij> j=1,2,3
         DDY1 = FDY1/3.d0 !<2*(dphi1/dy)*phij> j=1,2,3
         DXX11 = FDX1O2A*FDX1 !<2*(dphi1/dx)*(dphi1/dx)>
         DYY11 = FDY1O2A*FDY1 !<2*(dphi1/dy)*(dphi1/dy)>
         DXXYY11 = DXX11 + DYY11
         DXX12 = FDX1O2A*FDX2 !<2*(dphi1/dx)*(dphi2/dx)>
         DYY12 = FDY1O2A*FDY2 !<2*(dphi1/dy)*(dphi2/dy)>
         DXXYY12 = DXX12 + DYY12
         DXX13 = FDX1O2A*FDX3 !<2*(dphi1/dx)*(dphi3/dx)>
         DYY13 = FDY1O2A*FDY3 !<2*(dphi1/dy)*(dphi3/dy)>
         DXXYY13 = DXX13 + DYY13

         DDX2 = FDX2/3.d0 !<2*(dphi2/dx)*phij> j=1,2,3
         DDY2 = FDY2/3.d0 !<2*(dphi2/dy)*phij> j=1,2,3
         DXX22 = FDX2O2A*FDX2 !<2*(dphi2/dx)*(dphi2/dx)>
         DYY22 = FDY2O2A*FDY2 !<2*(dphi2/dy)*(dphi2/dy)>
         DXXYY22 = DXX22 + DYY22
         DXX23 = FDX2O2A*FDX3 !<2*(dphi2/dx)*(dphi3/dx)>
         DYY23 = FDY2O2A*FDY3 !<2*(dphi2/dy)*(dphi3/dy)>
         DXXYY23 = DXX23 + DYY23

         DDX3 = FDX3/3.d0 !<2*(dphi3/dx)*phij> j=1,2,3
         DDY3 = FDY3/3.d0 !<2*(dphi3/dy)*phij> j=1,2,3
         DXX33 = FDX3O2A*FDX3 !<2*(dphi3/dx)*(dphi3/dx)>
         DYY33 = FDY3O2A*FDY3 !<2*(dphi3/dy)*(dphi3/dy)>
         DXXYY33 = DXX33 + DYY33

         LUMPT = 1 !=1/0; LUMP/DO NOT LUMP THE TRANSPORT EQN
         FDDD = dble(1 + LUMPT)*AREAIE2/6.d0 !<2*(phii*phij) i=j>
         FDDOD = dble(1 - LUMPT)*AREAIE2/12.d0 !<2*(phii*phij) i<>j>
         FDDDODT = FDDD/DTDP
         FDDODODT = FDDOD/DTDP

         !.....COMPUTE ELEMENTAL QUANTITIES

         UEA = (U1N1 + U1N2 + U1N3)/3.d0
         VEA = (V1N1 + V1N2 + V1N3)/3.d0
         HEA = (H1N1 + H1N2 + H1N3)/3.d0
         EVCEA = (EVC1 + EVC2 + EVC3)/3.d0
         DHDX = H1N1*FDX1O2A + H1N2*FDX2O2A + H1N3*FDX3O2A
         DHDY = H1N1*FDY1O2A + H1N2*FDY2O2A + H1N3*FDY3O2A
         UPEA = UEA + DHDX*EVCEA/HEA
         VPEA = VEA + DHDY*EVCEA/HEA

         !.....ASSEMBLE PARTIAL PRODUCT

         CHSUM = CH1N1 + CH1N2 + CH1N3

         !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A1 AND
         !.....TEMP_LV_B1 VECTORS FOR NODE NM1

         TEMP_LV_B1 = & !LOAD VECTOR
            FDDDODT*CH1N1 + FDDODODT*(CH1N2 + CH1N3) &
            - EVCEA*(DXXYY11*CH1N1 + DXXYY12*CH1N2 + DXXYY13*CH1N3) &
            + (UPEA*DDX1 + VPEA*DDY1)*CHSUM &
            + FDDD*SS1N1 + FDDOD*(SS1N2 + SS1N3)
         TEMP_LV_A1 = & !LHS VECTOR
            FDDDODT + 2.d0*FDDODODT

         !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A2 AND
         !.....TEMP_LV_B2 VECTOR FOR NODE NM2

         TEMP_LV_B2 = & !LOAD VECTOR
            FDDDODT*CH1N2 + FDDODODT*(CH1N1 + CH1N3) &
            - EVCEA*(DXXYY12*CH1N1 + DXXYY22*CH1N2 + DXXYY23*CH1N3) &
            + (UPEA*DDX2 + VPEA*DDY2)*CHSUM &
            + FDDD*SS1N2 + FDDOD*(SS1N1 + SS1N3)
         TEMP_LV_A2 = & !LHS VECTOR
            FDDDODT + 2.d0*FDDODODT

         !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A3 AND
         !.....TEMP_LV_B3 VECTOR FOR NODE NM3

         TEMP_LV_B3 = & !LOAD VECTOR
            FDDDODT*CH1N3 + FDDODODT*(CH1N1 + CH1N2) &
            - EVCEA*(DXXYY13*CH1N1 + DXXYY23*CH1N2 + DXXYY33*CH1N3) &
            + (UPEA*DDX3 + VPEA*DDY3)*CHSUM &
            + FDDD*SS1N3 + FDDOD*(SS1N1 + SS1N2)
         TEMP_LV_A3 = & !LHS VECTOR
            FDDDODT + 2.d0*FDDODODT

         !     FINALIZE THE ASSEMBLY PROCESS FOR QC AND TRANS_LV_A
         !     USING THE TEMPORARY VECTORS
         TRANS_LV_B(NM1) = TRANS_LV_B(NM1) + TEMP_LV_B1*dble(NCELE) !LOAD VECTOR
         TRANS_LV_B(NM2) = TRANS_LV_B(NM2) + TEMP_LV_B2*dble(NCELE) !LOAD VECTOR
         TRANS_LV_B(NM3) = TRANS_LV_B(NM3) + TEMP_LV_B3*dble(NCELE) !LOAD VECTOR
         TRANS_LV_A(NM1) = TRANS_LV_A(NM1) + TEMP_LV_A1*dble(NCELE) !LUMPED LHS MATRIX
         TRANS_LV_A(NM2) = TRANS_LV_A(NM2) + TEMP_LV_A2*dble(NCELE) !LUMPED LHS MATRIX
         TRANS_LV_A(NM3) = TRANS_LV_A(NM3) + TEMP_LV_A3*dble(NCELE) !LUMPED LHS MATRIX

      end do

      !.... SOLVE FOR C*H NODE BY NODE

      do I = 1, NP
         NCI = NODECODE(I)
         if (NCI /= 0) CH1(I) = TRANS_LV_B(I)/TRANS_LV_A(I)
         !     IF(LBArray_Pointer(I).NE.0) CH1(I)=0.d0  !ESSENTIAL C=0 BOUNDARY CONDITION
      end do

   end subroutine SCALAR_TRANS_2D
   !-----------------------------------------------------------------------

end module mod_transport_2d
