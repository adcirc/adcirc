C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C 
C                       M O D U L E   M O M E N T U M
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C
C Sets up and solves the momentum equations. 
C
C******************************************************************************
C Logical Variable List (default value .FALSE., set in global.f)              *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     C2DDI            - 2D Depth Integrated model run                        *
C     C3D              - 3D model run                                         *
C     C3DDSS           - Stress form of the 3D momentum equations             *
C     C3DVS            - Velocity form of the 3D momentum equations           *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     C2D_BTrans       - Include a 2D baroclinic transport calculation        *
C     C2D_PTrans       - Include a 2D passive transport calculation           *
C     C3D_BTrans       - Include a 3D baroclinic transport calculation        *
C                        (used in 3D subroutines only)                        *
C     C3D_PTrans       - Include a 3D passive transport calculation           *
C                        (used in 3D subroutines only)                        *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CBaroclinic      - Include baroclinic terms                             *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CME_Orig         - Original Momentum Eq. formulation                    *
C     CME_New_NC       - Non-conservative advection formulation in the        *
C                        Momentum Eqs. (same as original formulation)         *
C     CME_New_C1       - Conservative advection formulation 1 in the          *
C                        Momentum Eqs.                                        *
C     CME_New_C2       - Conservative advection formulation 2 in the          *
C                        Momentum Eqs.                                        *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     The following are not available in Mom_Eqs_Original()                   *
C     CME_LS_IBPQ      - Integration by parts, flux-based, lateral stress     *
C                        formulation in the Momentum Eqs.                     *
C     CME_LS_IBPV      - Integration by parts, velocity-based, lateral stress *
C                        formulation in the Momentum Eqs.                     *
C                        (same as original formulation)                       *
C     CME_LS_IBPSQ     - Integration by parts, flux-based, symmetric lateral  *
C                        stress formulation in the Momentum Eqs.              *
C     CME_LS_IBPSV     - Integration by parts, velocity-based, symmetric      *
C                        lateral stress formulation in the Momentum Eqs.      *
C     CME_LS_2PartQ    - 2 Part, flux-based, lateral stress formulation in    *
C                        the Momentum Eqs.  (NOT IMPLEMENTED)                 *
C     CME_LS_2PartV    - 2 Part, velocity-based, lateral stress formulation in*
C                        the Momentum Eqs.  (NOT IMPLEMENTED)                 *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CME_AreaInt_Orig - Original area integration in the Momentum Eqs.       *
C                        (incorrect, but same as original formulation)        *
C     CME_AreaInt_Corr - Corrected area integration in the Momentum Eqs.      *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CSmag_Eh         - Use Smagorinski, spatially varying, vertically       *
C                        integrated lateral viscosity coefficient             *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C   See header.f for a summary history of code modifications.                 *
C******************************************************************************
      module momentum
      
      contains

C*******************************************************************************
C                                                                              *
C   Subroutine to compute the velocity and from that the flux/unit width using *
C   a 2DDI non conservative momentum equation.                                 *
C                                                                              *
C   Options are provided for either the correct area integration or the        *
C   original incorrect area integration.                                       *
C                                                                              *
C   Options are provided to use either flux or velocity based lateral          *
C   viscosity.                                                                 *
C                                                                              *
C   For a uniform grid and velocity based lateral viscosity, this subroutine   *
C   should give the same results as the original nonconservative formulation.  *
C                                                                              *
C   This subroutine follows the naming convention and formulation in the new   *
C   ADCIRC theory report.                                                      *
C                                                                              *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE Mom_Eqs_New_NC()
C
      USE GLOBAL
      USE WIND
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, Areas, TotalArea, 
     &                  NNeigh, NeiTab, NeiTabEle, MJU, SFAC
      USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI,
     &    CSII, SIII, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3                 
      USE NodalAttributes, ONLY: EVM,LoadAdvectionState
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      USE Couple2Swan, ONLY: TKXX,
     &                       TKXY,
     &                       TKYY
#endif
#endif
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceUVOB, 
     &                       enforceUVCB

      IMPLICIT NONE

      INTEGER IE, I, J, N                           !local loop counters
      INTEGER  NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI
      INTEGER NBDI
      INTEGER NNFirst

      REAL(SZ) BTP1N1, BTP1N2, BTP1N3, BTP2N1, BTP2N2, BTP2N3
      REAL(SZ) DBTPDXA, DBTPDYA
      REAL(SZ) DDU
      REAL(SZ) DQX1DX, DQX1DY, DQY1DX, DQY1DY
      REAL(SZ) DU1DX, DU1DY, DV1DX, DV1DY
      REAL(SZ) DU1DXA, DU1DYA, DV1DXA, DV1DYA
      REAL(SZ) EVMH1N1, EVMH1N2, EVMH1N3, EVMEle, EVMSmag
      REAL(SZ) H1, H1N1, H1N2, H1N3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) LSXXN1, LSXXN2, LSXXN3
      REAL(SZ) LSXYN1, LSXYN2, LSXYN3
      REAL(SZ) LSYXN1, LSYXN2, LSYXN3
      REAL(SZ) LSYYN1, LSYYN2, LSYYN3
      REAL(SZ) QTan
      REAL(SZ) QX1N1, QX1N2, QX1N3
      REAL(SZ) QY1N1, QY1N2, QY1N3
      REAL(SZ) SFacAvg
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) U1N1, U1N2, U1N3, U1Avg
      REAL(SZ) U1AvgDU1DXA, U1AvgDV1DXA
      REAL(SZ) V1N1, V1N2, V1N3, V1Avg
      REAL(SZ) V1AvgDU1DYA, V1AvgDV1DYA
      REAL(SZ) VCoef1, VCoef2
      REAL(SZ) VelNorm, VelTan
      REAL(SZ) VIDBCPDX, VIDBCPDY
      REAL(SZ) WSX, WSY
      REAL(SZ) ZNGLHS,ZNGRHS1,ZNGRHS2

      REAL(8) AreaIE, AreaIE2
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      REAL(SZ) :: AUVXX(NP)
      REAL(SZ) :: AUVXY(NP)
      REAL(SZ) :: AUVYX(NP)
      REAL(SZ) :: AUVYY(NP)
      REAL(SZ) :: VCOEFXX
      REAL(SZ) :: VCOEFXY
      REAL(SZ) :: VCOEFYX
      REAL(SZ) :: VCOEFYY
#endif
#endif
      call setMessageSource("mom_eqs_new_nc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
C...        THE TIME STEPPING LOOP.

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT


Corbitt 120322: Localized Advection
         IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         H1N1=DP(NM1)+IFNLFA*ETA1(NM1)
         H1N2=DP(NM2)+IFNLFA*ETA1(NM2)
         H1N3=DP(NM3)+IFNLFA*ETA1(NM3)
         H2N1=DP(NM1)+IFNLFA*ETA2(NM1)
         H2N2=DP(NM2)+IFNLFA*ETA2(NM2)
         H2N3=DP(NM3)+IFNLFA*ETA2(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         AreaIE2=Areas(IE)
         AreaIE =AreaIE2/2.d0
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg  !b1
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg  !b2
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg  !b3
         FDY1=X(NM3)-X(NM2)            !a1
         FDY2=X(NM1)-X(NM3)            !a2
         FDY3=X(NM2)-X(NM1)            !a3

C...  Compute element averaged quantities

         U1Avg =(U1N1+U1N2+U1N3)/3.d0
         V1Avg =(V1N1+V1N2+V1N3)/3.d0

         DU1DXA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
         DU1DYA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
         DV1DXA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
         DV1DYA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0

         EVMEle=(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
         !If using Smagorinski vertically-integrated lateral stress coefficient
         IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
            EVMSmag=EVMEle*
     &             sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA)
     &                 +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
            !tcm v52.30.01 added test for limits
            IF(Smag_Comp_Flag) then
               IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
               IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
            ENDIF
            EVMEle=EVMSmag
         ENDIF         

C...  Compute terms associated with the barotropic pressure

         BTP1N1=ETA1(NM1)
         BTP2N1=ETA2(NM1)
         BTP1N2=ETA1(NM2)
         BTP2N2=ETA2(NM2)
         BTP1N3=ETA1(NM3)
         BTP2N3=ETA2(NM3)

C.......If using atm pressure add it into the barotropic pressure

         IF(NWS.NE.0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
            ENDIF

C.......If using tidal potential terms, add these into the barotropic pressure

         IF (CTIP) THEN
            BTP1N1=BTP1N1-TiP1(NM1)
            BTP2N1=BTP2N1-TiP2(NM1)
            BTP1N2=BTP1N2-TiP1(NM2)
            BTP2N2=BTP2N2-TiP2(NM2)
            BTP1N3=BTP1N3-TiP1(NM3)
            BTP2N3=BTP2N3-TiP2(NM3)
            ENDIF

C...  Compute the barotropic pressure gradient x area for the element

         DBTPDXA=((BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3)
     &           +(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
         DBTPDYA=((BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3)
     &           +(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

C...  Compute the advective term gradients x area for the element

         U1AvgDU1DXA=U1Avg*DU1DXA
         V1AvgDU1DYA=V1Avg*DU1DYA
         U1AvgDV1DXA=U1Avg*DV1DXA
         V1AvgDV1DYA=V1Avg*DV1DYA

C...  Compute the lateral viscous terms for the element (flux formulation)

         IF (CME_LS_IBPQ) THEN
            DQX1DX=(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)/AreaIE2
            DQX1DY=(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)/AreaIE2
            DQY1DX=(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)/AreaIE2
            DQY1DY=(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=EVMEle*DQX1DY
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=EVMEle*DQY1DX
            LSYXN2=LSYXN1
            LSYXN3=LSYXN1
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric flux formulation)

         IF (CME_LS_IBPSQ) THEN
            DQX1DX=(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)/AreaIE2
            DQX1DY=(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)/AreaIE2
            DQY1DX=(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)/AreaIE2
            DQY1DY=(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=0.5d0*EVMEle*(DQX1DY+DQY1DX)
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (velocity formulation)

         IF (CME_LS_IBPV) THEN
            DU1DX=DU1DXA/AreaIE
            DU1DY=DU1DYA/AreaIE
            DV1DX=DV1DXA/AreaIE
            DV1DY=DV1DYA/AreaIE
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=EVMH1N1*DU1DY
            LSXYN2=EVMH1N2*DU1DY
            LSXYN3=EVMH1N3*DU1DY
            LSYXN1=EVMH1N1*DV1DX
            LSYXN2=EVMH1N2*DV1DX
            LSYXN3=EVMH1N3*DV1DX
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

         IF (CME_LS_IBPSV) THEN
            DU1DX=DU1DXA/AreaIE
            DU1DY=DU1DYA/AreaIE
            DV1DX=DV1DXA/AreaIE
            DV1DY=DV1DYA/AreaIE
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=0.5d0*EVMH1N1*(DU1DY+DV1DX)
            LSXYN2=0.5d0*EVMH1N2*(DU1DY+DV1DX)
            LSXYN3=0.5d0*EVMH1N3*(DU1DY+DV1DX)
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM1

         TEMP_LV_A1=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1)/H1N1
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM2
C...
         TEMP_LV_A2=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2)/H1N2
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM3
C...
         TEMP_LV_A3=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3)/H1N3
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM1

         TEMP_LV_B1=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1)/H1N1
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM2

         TEMP_LV_B2=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2)/H1N2
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM3

         TEMP_LV_B3=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3)/H1N3
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C     Original (incorrect) area integration - for historical comparison

         IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
         ENDIF

C     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
         TEMP_LV_A(IE,1)=TEMP_LV_A1
         TEMP_LV_A(IE,2)=TEMP_LV_A2
         TEMP_LV_A(IE,3)=TEMP_LV_A3
         TEMP_LV_B(IE,1)=TEMP_LV_B1
         TEMP_LV_B(IE,2)=TEMP_LV_B2
         TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

C     LINES TO RUN ON A SCALAR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

         ENDDO


C     LINES TO RUN ON A VECTOR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND AUV
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A(IE,1)
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A(IE,2)
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A(IE,3)
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B(IE,1)
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B(IE,2)
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B(IE,3)
         END DO
#endif

C...  Update the momentum equation LHS coefficients and load vectors at each
C...  node by dividing by the area of all active elements attached to the node
C...  and adding in the lumped terms, bottom friction and boundary conditions

      WSX=0.D0
      WSY=0.D0
      VIDBCPDX=0.D0
      VIDBCPDY=0.D0
      DO I=1,NP
         NCI=NODECODE(I)
         IF(TotalArea(I).ne.0.d0) THEN
            IF (CME_AreaInt_Corr) THEN       !Correct area integration
               MOM_LV_X(I)=MOM_LV_X(I)/TotalArea(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/TotalArea(I)
            ENDIF
            IF (CME_AreaInt_Orig) THEN       !Original (incorrect) area integration
               MOM_LV_X(I)=MOM_LV_X(I)/MJU(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/MJU(I)
            ENDIF
         ENDIF
         H1=DP(I)+IFNLFA*ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
         ENDIF
         VCoef1=DTO2*TK(I)                     !TK = Kslip/H
         VCoef2=DTO2*CORIF(I)
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         VCOEFXX = DTO2*TKXX(I)
         VCOEFYY = DTO2*TKYY(I)
         VCOEFXY = DTO2*(TKXY(I)-CORIF(I))
         VCOEFYX = DTO2*(TKXY(I)+CORIF(I))
#endif
#endif
         IF(CBaroclinic) THEN
            VIDBCPDX=DT*VIDBCPDXOH(I)
            VIDBCPDY=DT*VIDBCPDYOH(I)
         ENDIF

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCOEFXX)*UU1(I)-VCOEFXY*VV1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCOEFYY)*VV1(I)-VCOEFYX*UU1(I)-VIDBCPDX)
#else
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I)
     &                                     +VCoef2*VV1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I)
     &                                     -VCoef2*UU1(I)-VIDBCPDY)
#endif
#else
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I)
     &                                     +VCoef2*VV1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I)
     &                                     -VCoef2*UU1(I)-VIDBCPDY)
#endif

         AUV11(I)=1.D0+VCoef1*NCI
         AUV12(I)=-VCoef2*NCI

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         AUVXX(I) = 1.D0 + VCOEFXX*NCI
         AUVYY(I) = 1.D0 + VCOEFYY*NCI
         AUVXY(I) = VCOEFXY*NCI
         AUVYX(I) = VCOEFXY*NCI
#endif
#endif

         END DO

C...  Modify the momentum equations to impose velocity boundary
C...  conditions In each case the equations are manipulated to
C...  maintain the LHS matrix structure of AUV11=AUV22;
C...  AUV12=-AUV21)

      DO J=1,NVELME
         I=ME2GW(J)
         NBDI=NBV(I)
         H2=DP(NBDI)+IFNLFA*ETA2(NBDI)
         NCI=NODECODE(NBDI)

C      Specified essential normal flow and free tangential slip

         IF((LBCODEI(I).GE.0).AND.(LBCODEI(I).LE.9)) THEN
            VelNorm=-QN2(I)/H2
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VelNorm*AUVXY(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUVXX(NBDI)*NCI   !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVXY(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
            ENDIF

C     Specified essential normal flow and no tangential slip

         IF((LBCODEI(I).GE.10).AND.(LBCODEI(I).LE.19)) THEN
            VelNorm=-QN2(I)/H2
            VelTan=0.D0
            MOM_LV_X(NBDI)=VelTan*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*NCI !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI) = SIII(I)
            AUVXY(NBDI) = - CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#endif
#endif
            ENDIF

C     Zero normal velocity gradient using a Galerkin approximation to
C     the normal derivatives. Note: this is fully explicit and therefore
C     the velocity at the boundary is computed entirely from surrounding
C     velocities at the previous time step.

         IF(LBCODEI(I).EQ.41) THEN
            NM1=NBDI
            ZNGRHS1=0.d0     !Zero Norm Grad of U Eqn
            ZNGRHS2=0.d0     !Zero Norm Grad of V Eqn
            ZNGLHS=0.d0
            NM2=NeiTab(NBDI,2) !operate on 1st neighbor
            NNFirst=NM2      !save these values until end
            DO N=3,NNeigh(NBDI) !operate on rest of neighbors
               NM3=NM2       !shift previously computed values
               NM2=NEITAB(NBDI,N) !select new neighbor to work on
               SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
               NEle=NeiTabEle(NBDI,N-2) !element # defined by nodes NM1,NM2,NM3
               NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
               IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN !if element is active, compute contribution
                  FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                  FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                  FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                  FDY1 = X(NM3)-X(NM2) !a1
                  FDY2 = X(NM1)-X(NM3) !a2
                  FDY3 = X(NM2)-X(NM1) !a3
                  ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                  ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                  ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                  ENDIF
               END DO
            NM3=NM2          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN
               FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1 = X(NM3)-X(NM2) !a1
               FDY2 = X(NM1)-X(NM3) !a2
               FDY3 = X(NM2)-X(NM1) !a3
               ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
               ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
               ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
               ENDIF
            IF(NCI.EQ.0) THEN
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=0.d0
               MOM_LV_Y(NBDI)=0.d0
               ELSE
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=ZNGRHS1/ZNGLHS
               MOM_LV_Y(NBDI)=ZNGRHS2/ZNGLHS
               ENDIF
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI)=1.D0
            AUVXY(NBDI)=0.D0
#endif
#endif
            ENDIF

         ENDDO
C...
C...  SOLVE FOR VELOCITY AT NEW LEVEL  (s+1)
C...

C.....Note: This includes the comparison between MJU and NODELE to
C.....determine if the node is an interface node.  If MJU < NODELE the
C.....velocity can be zeroed out to obtain an essential zero velocity at
C.....interface nodes.

      DO I=1,NP
         AUV22=AUV11(I)
         AUV21=-AUV12(I)
         DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
         UU2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
         VV2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         AUVYY(I)= AUVXX(I)
         AUVYX(I)= - AUVXY(I)
         DDU = AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
         UU2(I)=(MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
         VV2(I)=(MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
c        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
c           UBAR2(I)=0.D0    !no slip and normal flux
c           VBAR2(I)=0.D0    !on wet/dry interface nodes
c           ENDIF
      END DO

      if(subdomainOn.and.enforceBN.eq.1) call enforceUVcb() ! NCSU Subdomain
      if(subdomainOn.and.enforceBN.eq.2) call enforceUVob() ! NCSU Subdomain
C...
C...  Impose a zero normal velocity gradient based on interpolating the
C...  velocity at a fictitious point in the interior of the domain,
C...  normal to a specified boundary node and setting the boundary
C...  velocity equal to the interpolated value at the fictitious point.
C...  Provided the fictitious point does not lie in an element that
C...  contains a boundary point, this is an entirely implicit
C...  calculation.
C...
      IF(NFLUXGBC.EQ.1) THEN
         DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I).EQ.40) THEN
               NM1=NM(NEleZNG(I),1)
               NM2=NM(NEleZNG(I),2)
               NM3=NM(NEleZNG(I),3)
               NC1=NODECODE(NM1)
               NC2=NODECODE(NM2)
               NC3=NODECODE(NM3)
               NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
               UU2(NBDI)=NCEle*(UU2(NM1)*ZNGIF1(I)+UU2(NM2)*ZNGIF2(I)
     &                                            +UU2(NM3)*ZNGIF3(I))
               VV2(NBDI)=NCEle*(VV2(NM1)*ZNGIF1(I)+VV2(NM2)*ZNGIF2(I)
     &                                            +VV2(NM3)*ZNGIF3(I))
               ENDIF
            ENDDO
         ENDIF

C...  Compute fluxes

      DO I=1,NP
         H2=DP(I)+IFNLFA*ETA2(I)
         QX2(I)=UU2(I)*H2
         QY2(I)=VV2(I)*H2
      ENDDO

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C*******************************************************************************
      END SUBROUTINE MOM_EQS_NEW_NC
C*******************************************************************************



C*******************************************************************************
C                                                                              *
C   Subroutine to compute the flux/unit width and from that the velocity using *
C   2DDI conservative momentum equation formulations version 1 or 2.           *
C                                                                              *
C   Options are provided for either the correct area integration or the        *
C   original incorrect area integration.                                       *
C                                                                              *
C   Options are provided to use either flux or velocity based lateral          *
C   viscosity.                                                                 *
C                                                                              *
C   This subroutine follows the naming convention and formulation in the new   *
C   ADCIRC theory report.                                                      *
C                                                                              *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE Mom_Eqs_New_Conserv()
C
      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, TotalArea, MJU,
     &                  NEITAB, NEITABELE, NNEIGH, SFAC 
      USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, 
     &   CSII, SIII, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
      USE WIND
      USE NodalAttributes, ONLY: EVM,LoadAdvectionState
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      USE Couple2Swan, ONLY: TKXX,
     &                       TKXY,
     &                       TKYY
#endif
#endif
      IMPLICIT NONE

      INTEGER IE, I, J, N                           !local loop counters
      INTEGER  NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI
      INTEGER NBDI
      INTEGER NNFirst

      REAL(SZ) BTP1N1, BTP1N2, BTP1N3, BTP2N1, BTP2N2, BTP2N3
      REAL(SZ) DBTPDXA, DBTPDYA
      REAL(SZ) DDU
      REAL(SZ) DU1DX, DU1DY, DV1DX, DV1DY
      REAL(SZ) DU1DX2A, DU1DY2A, DV1DX2A, DV1DY2A
      REAL(SZ) DQX1DX, DQX1DY, DQY1DX, DQY1DY
      REAL(SZ) DQX1DXA, DQX1DYA, DQY1DXA, DQY1DYA
      REAL(SZ) DQX1DX2A, DQX1DY2A, DQY1DX2A, DQY1DY2A
      REAL(SZ) DU1QX1DXA, DV1QX1DYA, DU1QY1DXA, DV1QY1DYA
      REAL(SZ) EVMH1N1, EVMH1N2, EVMH1N3, EVMEle, EVMSmag
      REAL(SZ) H1, H1N1, H1N2, H1N3, H1Avg
      REAL(SZ) H2, H2N1, H2N2, H2N3, H2Avg
      REAL(SZ) LSXXN1, LSXXN2, LSXXN3
      REAL(SZ) LSXYN1, LSXYN2, LSXYN3
      REAL(SZ) LSYXN1, LSYXN2, LSYXN3
      REAL(SZ) LSYYN1, LSYYN2, LSYYN3
      REAL(SZ) QTan
      REAL(SZ) QX1N1, QX1N2, QX1N3, QX1Avg, QX1DU1DXA, QX1DV1DYA
      REAL(SZ) QY1N1, QY1N2, QY1N3, QY1Avg, QY1DU1DXA, QY1DV1DYA
      REAL(SZ) SFacAvg
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) U1N1, U1N2, U1N3, U1Avg, U1DQX1DXA, U1DQY1DXA
      REAL(SZ) V1N1, V1N2, V1N3, V1Avg, V1DQX1DYA, V1DQY1DYA
      REAL(SZ) VCoef1, VCoef2
      REAL(SZ) VIDBCPDX, VIDBCPDY
      REAL(SZ) WSX, WSY
      REAL(SZ) ZNGLHS,ZNGRHS1,ZNGRHS2

      REAL(8) AREAIE, AREAIE2
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      REAL(SZ) :: AUVXX(NP)
      REAL(SZ) :: AUVXY(NP)
      REAL(SZ) :: AUVYX(NP)
      REAL(SZ) :: AUVYY(NP)
      REAL(SZ) :: VCOEFXX
      REAL(SZ) :: VCOEFXY
      REAL(SZ) :: VCOEFYX
      REAL(SZ) :: VCOEFYY
#endif
#endif
      call setMessageSource("mom_eqs_new_conserv")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
C...        THE TIME STEPPING LOOP.

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT


Corbitt 120322: Localized Advection
         IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         H1N1=DP(NM1)+IFNLFA*ETA1(NM1)
         H1N2=DP(NM2)+IFNLFA*ETA1(NM2)
         H1N3=DP(NM3)+IFNLFA*ETA1(NM3)
         H2N1=DP(NM1)+IFNLFA*ETA2(NM1)
         H2N2=DP(NM2)+IFNLFA*ETA2(NM2)
         H2N3=DP(NM3)+IFNLFA*ETA2(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         AreaIE2=Areas(IE)
         AreaIE =AreaIE2/2.d0
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg  !b1
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg  !b2
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg  !b3
         FDY1=X(NM3)-X(NM2)            !a1
         FDY2=X(NM1)-X(NM3)            !a2
         FDY3=X(NM2)-X(NM1)            !a3

C...  Compute element averaged quantities

         H1Avg = (H1N1+H1N2+H1N3)/3.d0
         H2Avg = (H2N1+H2N2+H2N3)/3.d0

         DQX1DX2A=QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3
         DQX1DY2A=QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3
         DQY1DX2A=QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3
         DQY1DY2A=QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3

         DU1DX2A=U1N1*FDX1+U1N2*FDX2+U1N3*FDX3
         DU1DY2A=U1N1*FDY1+U1N2*FDY2+U1N3*FDY3
         DV1DX2A=V1N1*FDX1+V1N2*FDX2+V1N3*FDX3
         DV1DY2A=V1N1*FDY1+V1N2*FDY2+V1N3*FDY3

         EVMEle=(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0 
         !If using Smagorinski vertically-integrated lateral stress coefficient         
         IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
            EVMSmag=0.5d0*EVMEle*
     &             sqrt((DU1DX2A-DV1DY2A)*(DU1DX2A-DV1DY2A)
     &                 +(DU1DY2A+DV1DX2A)*(DU1DY2A+DV1DX2A))
            !tcm v52.30.01 added test for limits
            IF(Smag_Comp_Flag) then
               IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
               IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
            ENDIF
            EVMEle=EVMSmag
         ENDIF

C...  Compute terms associated with the barotropic pressure

         BTP1N1=ETA1(NM1)
         BTP2N1=ETA2(NM1)
         BTP1N2=ETA1(NM2)
         BTP2N2=ETA2(NM2)
         BTP1N3=ETA1(NM3)
         BTP2N3=ETA2(NM3)

C.......If using atm pressure add it into the barotropic pressure

         IF(NWS.NE.0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
         ENDIF

C.......If using tidal potential terms, add these into the barotropic pressure

         IF (CTIP) THEN
            BTP1N1=BTP1N1-TIP1(NM1)
            BTP2N1=BTP2N1-TIP2(NM1)
            BTP1N2=BTP1N2-TIP1(NM2)
            BTP2N2=BTP2N2-TIP2(NM2)
            BTP1N3=BTP1N3-TIP1(NM3)
            BTP2N3=BTP2N3-TIP2(NM3)
         ENDIF

C...  Compute the element avg depth x barotropic pressure gradient x area for the element

         DBTPDXA=(H1Avg*(BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3)
     &           +H2Avg*(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
         DBTPDYA=(H1Avg*(BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3)
     &           +H2Avg*(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

C...  Compute the advective term gradients for the element version 1

         IF (CME_New_C1) THEN
            DU1QX1DXA=(U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2+U1N3*QX1N3*FDX3)
     &                                                             /2.d0
            DV1QX1DYA=(V1N1*QX1N1*FDY1+V1N2*QX1N2*FDY2+V1N3*QX1N3*FDY3)
     &                                                             /2.d0
            DU1QY1DXA=(U1N1*QY1N1*FDX1+U1N2*QY1N2*FDX2+U1N3*QY1N3*FDX3)
     &                                                             /2.d0
            DV1QY1DYA=(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2+V1N3*QY1N3*FDY3)
     &                                                             /2.d0
         ENDIF

C...  Compute the advective term gradients for the element version 2

         IF (CME_New_C2) THEN
            U1Avg =(U1N1+U1N2+U1N3)/3.d0
            V1Avg =(V1N1+V1N2+V1N3)/3.d0
            QX1Avg=(QX1N1+QX1N2+QX1N3)/3.d0
            QY1Avg=(QY1N1+QY1N2+QY1N3)/3.d0
            U1DQX1DXA=U1Avg *DQX1DX2A/2.d0
            QX1DU1DXA=QX1Avg*DU1DX2A/2.d0
            V1DQX1DYA=V1Avg *DQX1DY2A/2.d0
            QX1DV1DYA=QX1Avg*DV1DY2A/2.d0
            U1DQY1DXA=U1Avg *DQY1DX2A/2.d0
            QY1DU1DXA=QY1Avg*DU1DX2A/2.d0
            V1DQY1DYA=V1Avg *DQY1DY2A/2.d0
            QY1DV1DYA=QY1Avg*DV1DY2A/2.d0
            DU1QX1DXA=U1DQX1DXA+QX1DU1DXA
            DV1QX1DYA=V1DQX1DYA+QX1DV1DYA
            DU1QY1DXA=U1DQY1DXA+QY1DU1DXA
            DV1QY1DYA=V1DQY1DYA+QY1DV1DYA
         ENDIF

C...  Compute the lateral viscous terms for the element (flux formulation)

         IF (CME_LS_IBPQ) THEN
            DQX1DX=DQX1DX2A/AreaIE2
            DQX1DY=DQX1DY2A/AreaIE2
            DQY1DX=DQY1DX2A/AreaIE2
            DQY1DY=DQY1DY2A/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=EVMEle*DQX1DY
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=EVMEle*DQY1DX
            LSYXN2=LSYXN1
            LSYXN3=LSYXN1
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric flux formulation)

         IF (CME_LS_IBPSQ) THEN
            DQX1DX=DQX1DX2A/AreaIE2
            DQX1DY=DQX1DY2A/AreaIE2
            DQY1DX=DQY1DX2A/AreaIE2
            DQY1DY=DQY1DY2A/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=0.5d0*EVMEle*(DQX1DY+DQY1DX)
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (velocity formulation)

         IF (CME_LS_IBPV) THEN
            DU1DX=DU1DX2A/AreaIE2
            DU1DY=DU1DY2A/AreaIE2
            DV1DX=DV1DX2A/AreaIE2
            DV1DY=DV1DY2A/AreaIE2
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=EVMH1N1*DU1DY
            LSXYN2=EVMH1N2*DU1DY
            LSXYN3=EVMH1N3*DU1DY
            LSYXN1=EVMH1N1*DV1DX
            LSYXN2=EVMH1N2*DV1DX
            LSYXN3=EVMH1N3*DV1DX
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

         IF (CME_LS_IBPSV) THEN
            DU1DX=DU1DX2A/AreaIE2
            DU1DY=DU1DY2A/AreaIE2
            DV1DX=DV1DX2A/AreaIE2
            DV1DY=DV1DY2A/AreaIE2
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=0.5d0*EVMH1N1*(DU1DY+DV1DX)
            LSXYN2=0.5d0*EVMH1N2*(DU1DY+DV1DX)
            LSXYN3=0.5d0*EVMH1N3*(DU1DY+DV1DX)
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM1

         TEMP_LV_A1=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QX1DXA+DV1QX1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM2
C...
         TEMP_LV_A2=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QX1DXA+DV1QX1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM3
C...
         TEMP_LV_A3=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QX1DXA+DV1QX1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM1

         TEMP_LV_B1=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QY1DXA+DV1QY1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM2

         TEMP_LV_B2=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QY1DXA+DV1QY1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM3

         TEMP_LV_B3=NCELE*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*(DU1QY1DXA+DV1QY1DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3)
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C     Original (incorrect) area integration - for historical comparison

         IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
         ENDIF

C     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
         TEMP_LV_A(IE,1)=TEMP_LV_A1
         TEMP_LV_A(IE,2)=TEMP_LV_A2
         TEMP_LV_A(IE,3)=TEMP_LV_A3
         TEMP_LV_B(IE,1)=TEMP_LV_B1
         TEMP_LV_B(IE,2)=TEMP_LV_B2
         TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

C     LINES TO RUN ON A SCALAR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

         ENDDO



C     LINES TO RUN ON A VECTOR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND AUV
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A(IE,1)
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A(IE,2)
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A(IE,3)
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B(IE,1)
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B(IE,2)
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B(IE,3)
         END DO
#endif


C...  Update the momentum equation LHS coefficients and load vectors at each
C...  node by dividing by the area of all active elements attached to the node
C...  and adding in the lumped terms, bottom friction and boundary conditions

      WSX=0.D0
      WSY=0.D0
      VIDBCPDX=0.D0
      VIDBCPDY=0.D0
      DO I=1,NP
         NCI=NODECODE(I)
         IF(TotalArea(I).ne.0.d0) THEN
            IF (CME_AreaInt_Corr) THEN !Correct area integration
               MOM_LV_X(I)=MOM_LV_X(I)/TotalArea(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/TotalArea(I)
            ENDIF
            IF (CME_AreaInt_Orig) THEN !Original (incorrect) area integration
               MOM_LV_X(I)=MOM_LV_X(I)/MJU(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/MJU(I)
            ENDIF
         ENDIF
         H1=DP(I)+IFNLFA*ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)+WSX2(I))
            WSY=DTO2*IFWIND*(WSY1(I)+WSY2(I))
         ENDIF
         VCoef1=DTO2*TK(I)
         VCoef2=DTO2*CORIF(I)
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         VCOEFXX = DTO2*TKXX(I)
         VCOEFYY = DTO2*TKYY(I)
         VCOEFXY = DTO2*(TKXY(I)-CORIF(I))
         VCOEFYX = DTO2*(TKXY(I)+CORIF(I))
#endif
#endif
         VIDBCPDX=DT*VIDBCPDXOH(I)*H2
         VIDBCPDY=DT*VIDBCPDYOH(I)*H2

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCOEFXX)*QX1(I)-VCOEFXY*QY1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCOEFYY)*QY1(I)-VCOEFYX*QX1(I)-VIDBCPDX)
#else
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*QX1(I)
     &        +VCoef2*QY1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*QY1(I)
     &        -VCoef2*QX1(I)-VIDBCPDY)
#endif
#else
         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*QX1(I)
     &        +VCoef2*QY1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*QY1(I)
     &        -VCoef2*QX1(I)-VIDBCPDY)
#endif

         AUV11(I)=1.D0+VCoef1*NCI
         AUV12(I)=-VCoef2*NCI

#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         AUVXX(I) = 1.D0 + VCOEFXX*NCI
         AUVYY(I) = 1.D0 + VCOEFYY*NCI
         AUVXY(I) = VCOEFXY*NCI
         AUVYX(I) = VCOEFXY*NCI
#endif
#endif

      END DO

C...  Modify the momentum equations to impose velocity boundary
C...  conditions In each case the equations are manipulated to
C...  maintain the LHS matrix structure of AUV11=AUV22;
C...  AUV12=-AUV21)

      DO J=1,NVELME
         I=ME2GW(J)
         NBDI=NBV(I)
         NCI=NODECODE(NBDI)

C      Specified essential normal flow and free tangential slip

         IF((LBCODEI(I).GE.0).AND.(LBCODEI(I).LE.9)) THEN
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     +QN2(I)*AUVXY(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUVXX(NBDI)*NCI           !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVYX(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     +QN2(I)*AUV12(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUV11(NBDI)*NCI           !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     +QN2(I)*AUV12(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUV11(NBDI)*NCI           !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
            ENDIF

C     Specified essential normal flow and no tangential slip

         IF((LBCODEI(I).GE.10).AND.(LBCODEI(I).LE.19)) THEN
            QTAN=0.D0
            MOM_LV_X(NBDI)=QTan*NCI                          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*NCI                       !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI) = SIII(I)
            AUVXY(NBDI) = - CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#endif
#endif
            ENDIF

C     Zero normal flux gradient using a Galerkin approximation to
C     the normal derivatives. Note: this is fully explicit and therefore
C     the flux at the boundary is computed entirely from surrounding
C     fluxes at the previous time step.

         IF(LBCODEI(I).EQ.41) THEN
            NM1=NBDI
            ZNGRHS1=0.d0                                     !Zero Norm Grad of U Eqn
            ZNGRHS2=0.d0                                     !Zero Norm Grad of V Eqn
            ZNGLHS=0.d0
            NM2=NeiTab(NBDI,2)                               !operate on 1st neighbor
            NNFirst=NM2                                      !save these values until end
            DO N=3,NNeigh(NBDI)                              !operate on rest of neighbors
               NM3=NM2                                       !shift previously computed values
               NM2=NEITAB(NBDI,N)                            !select new neighbor to work on
               SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
               NEle=NeiTabEle(NBDI,N-2)                      !element # defined by nodes NM1,NM2,NM3
               NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
               IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN         !if element is active, compute contribution
                  FDX1 = (Y(NM2)-Y(NM3))*SFacAvg             !b1
                  FDX2 = (Y(NM3)-Y(NM1))*SFacAvg             !b2
                  FDX3 = (Y(NM1)-Y(NM2))*SFacAvg             !b3
                  FDY1 = X(NM3)-X(NM2)                       !a1
                  FDY2 = X(NM1)-X(NM3)                       !a2
                  FDY3 = X(NM2)-X(NM1)                       !a3
                  ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*QX1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*QX1(NM3)
                  ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*QY1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*QY1(NM3)
                  ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                  ENDIF
               END DO
            NM3=NM2                                          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN
               FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1 = X(NM3)-X(NM2) !a1
               FDY2 = X(NM1)-X(NM3) !a2
               FDY3 = X(NM2)-X(NM1) !a3
               ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*QX1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*QX1(NM3)
               ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*QY1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*QY1(NM3)
               ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
               ENDIF
            IF(NCI.EQ.0) THEN
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=0.d0
               MOM_LV_Y(NBDI)=0.d0
               ELSE
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=ZNGRHS1/ZNGLHS
               MOM_LV_Y(NBDI)=ZNGRHS2/ZNGLHS
               ENDIF
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI)=1.d0
            AUVXY(NBDI)=0.d0
#endif
#endif
            ENDIF

         ENDDO

C...
C...  SOLVE FOR FLUX AT NEW LEVEL  (K+1)
C...

C.....Note: This includes the comparison between MJU and NODELE to
C.....determine if the node is an interface node.  If MJU < NODELE the
C.....velocity can be zeroed out to obtain an essential zero velocity at
C.....interface nodes.

      DO I=1,NP
         AUV22=AUV11(I)
         AUV21=-AUV12(I)
         DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
         QX2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
         QY2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
         AUVYY(I)= AUVXX(I)
         AUVYX(I)= - AUVXY(I)
         DDU=AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
         QX2(I)=(MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
         QY2(I)=(MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
c        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
c           QX2(I)=0.D0    !no slip and normal flux
c           QY2(I)=0.D0    !on wet/dry interface nodes
c           ENDIF
       END DO

C...
C...  Impose a zero normal flux gradient based on interpolating the
C...  flux at a fictitious point in the interior of the domain,
C...  normal to a specified boundary node and setting the boundary
C...  flux equal to the interpolated value at the fictitious point.
C...  Provided the fictitious point does not lie in an element that
C...  contains a boundary point, this is an entirely implicit
C...  calculation.
C...
      IF(NFLUXGBC.EQ.1) THEN
         DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I).EQ.40) THEN
               NM1=NM(NEleZNG(I),1)
               NM2=NM(NEleZNG(I),2)
               NM3=NM(NEleZNG(I),3)
               NC1=NODECODE(NM1)
               NC2=NODECODE(NM2)
               NC3=NODECODE(NM3)
               NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
               QX2(NBDI)=NCEle*(QX2(NM1)*ZNGIF1(I)+QX2(NM2)*ZNGIF2(I)
     &                                            +QX2(NM3)*ZNGIF3(I))
               QY2(NBDI)=NCEle*(QY2(NM1)*ZNGIF1(I)+QY2(NM2)*ZNGIF2(I)
     &                                            +QY2(NM3)*ZNGIF3(I))
            ENDIF
         ENDDO
      ENDIF

C...  Compute velocities

      DO I=1,NP
         H2=DP(I)+IFNLFA*ETA2(I)
         IF(H2.NE.0.) THEN
            UU2(I)=QX2(I)/H2
            VV2(I)=QY2(I)/H2
            ELSE
            WRITE(16,*) ''
            WRITE(16,*) ''
            WRITE(16,*) '*******************************************'
            WRITE(16,*) '*******************************************'
            WRITE(16,*) 'WARNING: Total water depth = 0 at node  ',I
            WRITE(16,*) '         Velocities set = -999.'
            WRITE(16,*) '*******************************************'
            WRITE(16,*) '*******************************************'
            WRITE(16,*) ''
            WRITE(16,*) ''
            UU2(I)=-999.
            VV2(I)=-999.
            ENDIF
         ENDDO

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE MOM_EQS_NEW_CONSERV
C**********************************************************************


C*******************************************************************************
C                                                                              *
C   Subroutine to compute the velocity and from that the flux/unit width using *
C   a 2DDI non conservative momentum equation.                                 *
C                                                                              *
C   Options are provided for either the correct area integration or the        *
C   original incorrect area integration.                                       *
C                                                                              *
C   Options are provided to use either flux or velocity based lateral          *
C   viscosity.                                                                 *
C                                                                              *
C   For a uniform grid and velocity based lateral viscosity, this subroutine   *
C   should give the same results as the original nonconservative formulation.  *
C                                                                              *
C   This subroutine follows the naming convention and formulation in the new   *
C   ADCIRC theory report.                                                      *
C                                                                              *
C   This subroutine provides the corrector part of the momentum equation for   *
C   the predictor-corrector algorithm and obtains the corrected velocities.    *
C                                                                              *
C                            k.d.  06/24/2004                                  *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE Mom_Eqs_Non_Conserv_pc()
C
      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, TotalArea, MJU,
     &                  NEITAB, NEITABELE, NNEIGH, SFAC
      USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI,
     &   CSII, SIII, ME2GW, NBV, LBCODEI, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
      USE WIND
      USE NodalAttributes, ONLY: EVM,LoadAdvectionState

      IMPLICIT NONE

      INTEGER IE, I, J, N                           !local loop counters
      INTEGER  NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI
      INTEGER NBDI
      INTEGER NNFirst

      REAL(SZ) BTP1N1, BTP1N2, BTP1N3, BTP2N1, BTP2N2, BTP2N3
      REAL(SZ) DBTPDXA, DBTPDYA
      REAL(SZ) DDU
      REAL(SZ) DQX1DX, DQX1DY, DQY1DX, DQY1DY
      REAL(SZ) DU1DX, DU1DY, DV1DX, DV1DY
      REAL(SZ) DU1DXA, DU1DYA, DV1DXA, DV1DYA
      REAL(SZ) EVMH1N1, EVMH1N2, EVMH1N3, EVMEle, EVMSmag
      REAL(SZ) H1, H1N1, H1N2, H1N3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) LSXXN1, LSXXN2, LSXXN3
      REAL(SZ) LSXYN1, LSXYN2, LSXYN3
      REAL(SZ) LSYXN1, LSYXN2, LSYXN3
      REAL(SZ) LSYYN1, LSYYN2, LSYYN3
      REAL(SZ) QTan
      REAL(SZ) QX1N1, QX1N2, QX1N3
      REAL(SZ) QY1N1, QY1N2, QY1N3
      REAL(SZ) SFacAvg
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) U1N1, U1N2, U1N3, U1Avg
      REAL(SZ) U1AvgDU1DXA, U1AvgDV1DXA
      REAL(SZ) V1N1, V1N2, V1N3, V1Avg
      REAL(SZ) V1AvgDU1DYA, V1AvgDV1DYA
      REAL(SZ) VCoef1, VCoef2
      REAL(SZ) VelNorm, VelTan
      REAL(SZ) VIDBCPDX, VIDBCPDY
      REAL(SZ) WSX, WSY
      REAL(SZ) ZNGLHS,ZNGRHS1,ZNGRHS2

ckmd   Added in parameters for the pc algorithm
      REAL(SZ) QX0N1, QX0N2, QX0N3
      REAL(SZ) QY0N1, QY0N2, QY0N3
      REAL(SZ) QX2N1, QX2N2, QX2N3
      REAL(SZ) QY2N1, QY2N2, QY2N3
      REAL(SZ) U0N1, U0N2, U0N3, U0Avg
      REAL(SZ) U0AvgDU0DXA, U0AvgDV0DXA
      REAL(SZ) V0N1, V0N2, V0N3, V0Avg
      REAL(SZ) V0AvgDU0DYA, V0AvgDV0DYA
      REAL(SZ) U2N1, U2N2, U2N3, U2Avg
      REAL(SZ) U2AvgDU2DXA, U2AvgDV2DXA
      REAL(SZ) V2N1, V2N2, V2N3, V2Avg
      REAL(SZ) V2AvgDU2DYA, V2AvgDV2DYA
      REAL(SZ) DU0DXA, DU0DYA, DV0DXA, DV0DYA
      REAL(SZ) DU2DXA, DU2DYA, DV2DXA, DV2DYA
      REAL(SZ) timewtmom0,timewtmom1,timewtmom2,timebfflag
      REAL(SZ) VCoef12


      REAL(8) AreaIE, AreaIE2
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
C
      call setMessageSource("mom_eqs_non_conserv_pc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
C...        THE TIME STEPPING LOOP.
ckmd    Must reset the result vectors to zero before recomputing
ckmd  the next time level.
      DO I=1,NP
         MOM_LV_X(I)=0.D0
         MOM_LV_Y(I)=0.D0
      END DO

ckmd  Add in the time weights for the corrector loop
         timewtmom0=0.0d0
         timewtmom1=0.5d0
         timewtmom2=0.5d0
         timebfflag=1.0d0

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT


Corbitt 120322: Localized Advection
         IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         U0N1=UU0(NM1)
         U0N2=UU0(NM2)
         U0N3=UU0(NM3)
         V0N1=VV0(NM1)
         V0N2=VV0(NM2)
         V0N3=VV0(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         U2N1=UU2(NM1)
         U2N2=UU2(NM2)
         U2N3=UU2(NM3)
         V2N1=VV2(NM1)
         V2N2=VV2(NM2)
         V2N3=VV2(NM3)
         H1N1=DP(NM1)+IFNLFA*ETA1(NM1)
         H1N2=DP(NM2)+IFNLFA*ETA1(NM2)
         H1N3=DP(NM3)+IFNLFA*ETA1(NM3)
         H2N1=DP(NM1)+IFNLFA*ETA2(NM1)
         H2N2=DP(NM2)+IFNLFA*ETA2(NM2)
         H2N3=DP(NM3)+IFNLFA*ETA2(NM3)
         QX0N1=QX0(NM1)
         QX0N2=QX0(NM2)
         QX0N3=QX0(NM3)
         QY0N1=QY0(NM1)
         QY0N2=QY0(NM2)
         QY0N3=QY0(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         QX2N1=QX2(NM1)
         QX2N2=QX2(NM2)
         QX2N3=QX2(NM3)
         QY2N1=QY2(NM1)
         QY2N2=QY2(NM2)
         QY2N3=QY2(NM3)

         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         AreaIE2=Areas(IE)
         AreaIE =AreaIE2/2.d0
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg  !b1
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg  !b2
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg  !b3
         FDY1=X(NM3)-X(NM2)            !a1
         FDY2=X(NM1)-X(NM3)            !a2
         FDY3=X(NM2)-X(NM1)            !a3

C...  Compute element averaged quantities

         U0Avg =(U0N1+U0N2+U0N3)/3.d0
         V0Avg =(V0N1+V0N2+V0N3)/3.d0
         U1Avg =(U1N1+U1N2+U1N3)/3.d0
         V1Avg =(V1N1+V1N2+V1N3)/3.d0
         U2Avg =(U2N1+U2N2+U2N3)/3.d0
         V2Avg =(V2N1+V2N2+V2N3)/3.d0

         DU0DXA=(UU0(NM1)*FDX1+UU0(NM2)*FDX2+UU0(NM3)*FDX3)/2.d0
         DU0DYA=(UU0(NM1)*FDY1+UU0(NM2)*FDY2+UU0(NM3)*FDY3)/2.d0
         DV0DXA=(VV0(NM1)*FDX1+VV0(NM2)*FDX2+VV0(NM3)*FDX3)/2.d0
         DV0DYA=(VV0(NM1)*FDY1+VV0(NM2)*FDY2+VV0(NM3)*FDY3)/2.d0
         DU1DXA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
         DU1DYA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
         DV1DXA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
         DV1DYA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0
         DU2DXA=(UU2(NM1)*FDX1+UU2(NM2)*FDX2+UU2(NM3)*FDX3)/2.d0
         DU2DYA=(UU2(NM1)*FDY1+UU2(NM2)*FDY2+UU2(NM3)*FDY3)/2.d0
         DV2DXA=(VV2(NM1)*FDX1+VV2(NM2)*FDX2+VV2(NM3)*FDX3)/2.d0
         DV2DYA=(VV2(NM1)*FDY1+VV2(NM2)*FDY2+VV2(NM3)*FDY3)/2.d0

         EVMEle=(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
         !If using Smagorinski vertically-integrated lateral stress coefficient
         IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
            EVMSmag=EVMEle*
     &             sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA)
     &                 +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
            !tcm v52.30.01 added test for limits
            IF(Smag_Comp_Flag) then
               IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
               IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
            ENDIF
            EVMEle=EVMSmag
         ENDIF

C...  Compute terms associated with the barotropic pressure

         BTP1N1=ETA1(NM1)
         BTP2N1=ETA2(NM1)
         BTP1N2=ETA1(NM2)
         BTP2N2=ETA2(NM2)
         BTP1N3=ETA1(NM3)
         BTP2N3=ETA2(NM3)

C.......If using atm pressure add it into the barotropic pressure

         IF(NWS.NE.0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
            ENDIF

C.......If using tidal potential terms, add these into the barotropic pressure

         IF (CTIP) THEN
            BTP1N1=BTP1N1-TiP1(NM1)
            BTP2N1=BTP2N1-TiP2(NM1)
            BTP1N2=BTP1N2-TiP1(NM2)
            BTP2N2=BTP2N2-TiP2(NM2)
            BTP1N3=BTP1N3-TiP1(NM3)
            BTP2N3=BTP2N3-TiP2(NM3)
            ENDIF

C...  Compute the barotropic pressure gradient x area for the element

         DBTPDXA=((BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3)
     &           +(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
         DBTPDYA=((BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3)
     &           +(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

C...  Compute the advective term gradients x area for the element

         U0AvgDU0DXA=U0Avg*DU0DXA
         V0AvgDU0DYA=V0Avg*DU0DYA
         U0AvgDV0DXA=U0Avg*DV0DXA
         V0AvgDV0DYA=V0Avg*DV0DYA
         U1AvgDU1DXA=U1Avg*DU1DXA
         V1AvgDU1DYA=V1Avg*DU1DYA
         U1AvgDV1DXA=U1Avg*DV1DXA
         V1AvgDV1DYA=V1Avg*DV1DYA
         U2AvgDU2DXA=U2Avg*DU2DXA
         V2AvgDU2DYA=V2Avg*DU2DYA
         U2AvgDV2DXA=U2Avg*DV2DXA
         V2AvgDV2DYA=V2Avg*DV2DYA

C...  Compute the lateral viscous terms for the element (flux formulation)

         IF (CME_LS_IBPQ) THEN
            DQX1DX=(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)/AreaIE2
            DQX1DY=(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)/AreaIE2
            DQY1DX=(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)/AreaIE2
            DQY1DY=(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=EVMEle*DQX1DY
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=EVMEle*DQY1DX
            LSYXN2=LSYXN1
            LSYXN3=LSYXN1
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric flux formulation)

         IF (CME_LS_IBPSQ) THEN
            DQX1DX=(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)/AreaIE2
            DQX1DY=(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)/AreaIE2
            DQY1DX=(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)/AreaIE2
            DQY1DY=(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)/AreaIE2
            LSXXN1=EVMEle*DQX1DX
            LSXXN2=LSXXN1
            LSXXN3=LSXXN1
            LSXYN1=0.5d0*EVMEle*(DQX1DY+DQY1DX)
            LSXYN2=LSXYN1
            LSXYN3=LSXYN1
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMEle*DQY1DY
            LSYYN2=LSYYN1
            LSYYN3=LSYYN1
         ENDIF

C...  Compute the lateral viscous terms for the element (velocity formulation)

         IF (CME_LS_IBPV) THEN
            DU1DX=DU1DXA/AreaIE
            DU1DY=DU1DYA/AreaIE
            DV1DX=DV1DXA/AreaIE
            DV1DY=DV1DYA/AreaIE
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=EVMH1N1*DU1DY
            LSXYN2=EVMH1N2*DU1DY
            LSXYN3=EVMH1N3*DU1DY
            LSYXN1=EVMH1N1*DV1DX
            LSYXN2=EVMH1N2*DV1DX
            LSYXN3=EVMH1N3*DV1DX
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

         IF (CME_LS_IBPSV) THEN
            DU1DX=DU1DXA/AreaIE
            DU1DY=DU1DYA/AreaIE
            DV1DX=DV1DXA/AreaIE
            DV1DY=DV1DYA/AreaIE
            EVMH1N1=EVMEle*H1N1
            EVMH1N2=EVMEle*H1N2
            EVMH1N3=EVMEle*H1N3
            LSXXN1=EVMH1N1*DU1DX
            LSXXN2=EVMH1N2*DU1DX
            LSXXN3=EVMH1N3*DU1DX
            LSXYN1=0.5d0*EVMH1N1*(DU1DY+DV1DX)
            LSXYN2=0.5d0*EVMH1N2*(DU1DY+DV1DX)
            LSXYN3=0.5d0*EVMH1N3*(DU1DY+DV1DX)
            LSYXN1=LSXYN1
            LSYXN2=LSXYN2
            LSYXN3=LSXYN3
            LSYYN1=EVMH1N1*DV1DY
            LSYYN2=EVMH1N2*DV1DY
            LSYYN3=EVMH1N3*DV1DY
         ENDIF

C
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM1

         TEMP_LV_A1=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1)/H1N1
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM2
C...
         TEMP_LV_A2=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2)/H1N2
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A VECTOR FOR NODE NM3
C...
         TEMP_LV_A3=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDXA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3)/H1N3
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM1

         TEMP_LV_B1=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1)/H1N1
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM2

         TEMP_LV_B2=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2)/H1N2
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
C...  TEMP_LV_B VECTOR FOR NODE NM3

         TEMP_LV_B3=NCEle*DT*(
C...  ADVECTIVE TERMS
     &             -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA)
     &             -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA)
     &             -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA)
C...  BAROTROPIC PRESSURE GRADIENT
     &             -GO2*DBTPDYA
C...  LATERAL VISCOUS TERMS
     &             -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3)/H1N3
C...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
     &             )

C     Original (incorrect) area integration - for historical comparison

         IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
         ENDIF

C     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
         TEMP_LV_A(IE,1)=TEMP_LV_A1
         TEMP_LV_A(IE,2)=TEMP_LV_A2
         TEMP_LV_A(IE,3)=TEMP_LV_A3
         TEMP_LV_B(IE,1)=TEMP_LV_B1
         TEMP_LV_B(IE,2)=TEMP_LV_B2
         TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

C     LINES TO RUN ON A SCALAR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

         ENDDO


C     LINES TO RUN ON A VECTOR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
C           AND AUV
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A(IE,1)
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A(IE,2)
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A(IE,3)
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B(IE,1)
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B(IE,2)
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B(IE,3)
         END DO
#endif


C...  Update the momentum equation LHS coefficients and load vectors at each
C...  node by dividing by the area of all active elements attached to the node
C...  and adding in the lumped terms, bottom friction and boundary conditions

      WSX=0.D0
      WSY=0.D0
      VIDBCPDX=0.D0
      VIDBCPDY=0.D0
      DO I=1,NP
         NCI=NODECODE(I)
         IF(TotalArea(I).ne.0.d0) THEN
            IF (CME_AreaInt_Corr) THEN     !Correct area integration
               MOM_LV_X(I)=MOM_LV_X(I)/TotalArea(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/TotalArea(I)
            ENDIF
            IF (CME_AreaInt_Orig) THEN     !Original (incorrect) area integration
               MOM_LV_X(I)=MOM_LV_X(I)/MJU(I)
               MOM_LV_Y(I)=MOM_LV_Y(I)/MJU(I)
            ENDIF
         ENDIF
         H1=DP(I)+IFNLFA*ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
         ENDIF
ckmd
ckmd  Added time weights to tau terms in the momentum equation here.
ckmd  Note the weighting should follow the weighting in the GWCE.
ckmd
         VCoef1=DTO2*TK(I)                     !TK = Kslip/H
         VCoef12=DTO2*(TK2(I)*timebfflag+TK(I)*(1.D0-timebfflag))
         VCoef2=DTO2*CORIF(I)
         IF(CBaroclinic) THEN
            VIDBCPDX=DT*VIDBCPDXOH(I)
            VIDBCPDY=DT*VIDBCPDYOH(I)
            ENDIF

         MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I)
     &                                     +VCoef2*VV1(I)-VIDBCPDX)
         MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I)
     &                                     -VCoef2*UU1(I)-VIDBCPDY)

ckmd    Change for the corrector formulation
         AUV11(I)=1.D0+VCoef12*NCI
         AUV12(I)=-VCoef2*NCI
         END DO

C...  Modify the momentum equations to impose velocity boundary
C...  conditions In each case the equations are manipulated to
C...  maintain the LHS matrix structure of AUV11=AUV22;
C...  AUV12=-AUV21)

      DO J=1,NVELME
         I=ME2GW(J)
         NBDI=NBV(I)
         H2=DP(NBDI)+IFNLFA*ETA2(NBDI)
         NCI=NODECODE(NBDI)

C      Specified essential normal flow and free tangential slip

         IF((LBCODEI(I).GE.0).AND.(LBCODEI(I).LE.9)) THEN
            VelNorm=-QN2(I)/H2
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
            ENDIF

C     Specified essential normal flow and no tangential slip

         IF((LBCODEI(I).GE.10).AND.(LBCODEI(I).LE.19)) THEN
            VelNorm=-QN2(I)/H2
            VelTan=0.D0
            MOM_LV_X(NBDI)=VelTan*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*NCI !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
            ENDIF

C     Zero normal velocity gradient using a Galerkin approximation to
C     the normal derivatives. Note: this is fully explicit and therefore
C     the velocity at the boundary is computed entirely from surrounding
C     velocities at the previous time step.

         IF(LBCODEI(I).EQ.41) THEN
            NM1=NBDI
            ZNGRHS1=0.d0     !Zero Norm Grad of U Eqn
            ZNGRHS2=0.d0     !Zero Norm Grad of V Eqn
            ZNGLHS=0.d0
            NM2=NeiTab(NBDI,2) !operate on 1st neighbor
            NNFirst=NM2      !save these values until end
            DO N=3,NNeigh(NBDI) !operate on rest of neighbors
               NM3=NM2       !shift previously computed values
               NM2=NEITAB(NBDI,N) !select new neighbor to work on
               SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
               NEle=NeiTabEle(NBDI,N-2) !element # defined by nodes NM1,NM2,NM3
               NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
               IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN !if element is active, compute contribution
                  FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                  FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                  FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                  FDY1 = X(NM3)-X(NM2) !a1
                  FDY2 = X(NM1)-X(NM3) !a2
                  FDY3 = X(NM2)-X(NM1) !a3
                  ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                  ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2)
     &                           -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                  ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                  ENDIF
               END DO
            NM3=NM2          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle.NE.0).AND.(NCEle.NE.0)) THEN
               FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1 = X(NM3)-X(NM2) !a1
               FDY2 = X(NM1)-X(NM3) !a2
               FDY3 = X(NM2)-X(NM1) !a3
               ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
               ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2)
     &                        -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
               ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
               ENDIF
            IF(NCI.EQ.0) THEN
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=0.d0
               MOM_LV_Y(NBDI)=0.d0
               ELSE
               AUV11(NBDI)=1.d0
               AUV12(NBDI)=0.d0
               MOM_LV_X(NBDI)=ZNGRHS1/ZNGLHS
               MOM_LV_Y(NBDI)=ZNGRHS2/ZNGLHS
               ENDIF
            ENDIF

         ENDDO

C...
C...  SOLVE FOR VELOCITY AT NEW LEVEL  (s+1)
C...

C.....Note: This includes the comparison between MJU and NODELE to
C.....determine if the node is an interface node.  If MJU < NODELE the
C.....velocity can be zeroed out to obtain an essential zero velocity at
C.....interface nodes.

      DO I=1,NP
         AUV22=AUV11(I)
         AUV21=-AUV12(I)
         DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
         UU2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
         VV2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU

c        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
c           UBAR2(I)=0.D0    !no slip and normal flux
c           VBAR2(I)=0.D0    !on wet/dry interface nodes
c           ENDIF
         END DO

C...
C...  Impose a zero normal velocity gradient based on interpolating the
C...  velocity at a fictitious point in the interior of the domain,
C...  normal to a specified boundary node and setting the boundary
C...  velocity equal to the interpolated value at the fictitious point.
C...  Provided the fictitious point does not lie in an element that
C...  contains a boundary point, this is an entirely implicit
C...  calculation.
C...
      IF(NFLUXGBC.EQ.1) THEN
         DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I).EQ.40) THEN
               NM1=NM(NEleZNG(I),1)
               NM2=NM(NEleZNG(I),2)
               NM3=NM(NEleZNG(I),3)
               NC1=NODECODE(NM1)
               NC2=NODECODE(NM2)
               NC3=NODECODE(NM3)
               NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
               UU2(NBDI)=NCEle*(UU2(NM1)*ZNGIF1(I)+UU2(NM2)*ZNGIF2(I)
     &                                            +UU2(NM3)*ZNGIF3(I))
               VV2(NBDI)=NCEle*(VV2(NM1)*ZNGIF1(I)+VV2(NM2)*ZNGIF2(I)
     &                                            +VV2(NM3)*ZNGIF3(I))
               ENDIF
            ENDDO
         ENDIF

C...  Compute fluxes

      DO I=1,NP
         H2=DP(I)+IFNLFA*ETA2(I)
         QX2(I)=UU2(I)*H2
         QY2(I)=VV2(I)*H2
         ENDDO

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE MOM_EQS_NON_CONSERV_PC
C**********************************************************************


C----------------------------------------------------------------------
C----------------------------------------------------------------------
      end module momentum
C----------------------------------------------------------------------
C----------------------------------------------------------------------
