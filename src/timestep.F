C#define SB_WETDRY
C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C The timestepping module is configured to allow selection of a number of     *
C alternative algorithms within the overall FE framework.  These algorithms   *
C are selected by the TRUE/FALSE state of the logical variables listed below. *
C These variables are set in READ_INPUT.F, based on the value of the fort.15  *
C input parameter IM.  The only exception is CGWCE_Lump which is set by a     *
C preprocessor flag at compile time. The variables are passed in GLOBAL.F     *
C                                                                             *
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
C     CGWCE_Orig       - Original ADCIRC GWCE formulation                     *
C     CGWCE_New        - New ADCIRC GWCE formulation (old algorithm, new code)*
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CGWCE_Lump       - Lump the GWCE matrix                                 *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     GWCE_New() and GWCE_New_pc() only:                                      *
C     CGWCE_LS_KGQ     - Kolar-Gray, flux-based, lateral stress formulation   *
C                        in the GWCE (same as original formulation)           *
C     CGWCE_LS_2PartQ  - 2 Part, flux-based, lateral stress formulation       *
C                        in the GWCE                                          *
C     CGWCE_LS_2PartV  - 2 Part, velocity-based, lateral stress formulation   *
C                        in the GWCE                                          *
C     CGWCE_LS_2PartSQ - 2 Part, flux-based, symmetric lateral stress         *
C                        formulation in the GWCE                              *
C     CGWCE_LS_2PartSV - 2 Part, velocity-based, symmetric lateral stress     *
C                        formulation in the GWCE                              *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     GWCE_New() and GWCE_New_pc() only:                                      *
C     CGWCE_Advec_NC   - Non-conservative advection formulation in the GWCE   *
C                        (same as original formulation)                       *
C     CGWCE_Advec_C1   - Use conservative advection formulation 1 in the GWCE *
C     CGWCE_Advec_C2   - Use conservative advection formulation 2 in the GWCE *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CME_Orig         - Original Momentum Eq. formulation                    *
C     CME_New_NC       - Non-conservative advection formulation in the        *
C                        Momentum Eqs. (same as original formulation)         *
C     CME_New_C1       - Conservative advection formulation 1 in the          *
C                        Momentum Eqs.                                        *
C     CME_New_C2       - Conservative advection formulation 2 in the          *
C                        Momentum Eqs.                                        *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CPRECOR          - Use the predictor-corrector algorithm for GWCE       *
C                        and momentum equations (package deal)                *
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
C 
      SUBROUTINE TIMESTEP(IT)
C     
      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE GLOBAL_IO
      USE HARM, ONLY : updateHarmonicAnalysis
      USE WIND
C.....sb46.28sb01 added 09/xx/2006
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET
C.... TCM V49.64.01 ADDITIONS FOR ICE
      USE OWI_ICE,ONLY : NCICE1_INIT,NCICE1_GET
C.....sb46.28sb03 added 09/xx/2006
      USE RS2,ONLY : RS2INIT,RS2GET
      USE NodalAttributes, ONLY : FRIC, HBREAK, FTHETA, FGAMMA, ! for wet/dry 
     &     IFLINBF, IFNLBF, IFHYBF, ! for wet/dry
     &     Apply2DBottomFriction,
     &     Apply3DBottomFriction,
     &     ApplyDirectionalWindReduction,
     &     LoadDirEffRLen,
     &     ApplyCanopyCoefficient,
     &     LoadCanopyCoef,
     &     ManningsN, LoadManningsN, !sb46.28sb02
     &     Z0b_var, LoadZ0B_var,!Rjw IOOS 2010
     &     BFCdLLimit !sb46.28sb02/jgf47.04 lower limit of Cd for bot. friction


c. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
Cjgfdebug      USE GLOBAL_3DVS, ONLY:
Cjgfdebug     &     A, B, BSX, BSY, EVTOT, ISLIP, KP, Q, SIGMA
      USE GLOBAL_3DVS, ONLY:
     &     A, B, BSX, BSY, EVTOT, ISLIP, KP, Q, SIGMA, Z0B,NFEN

#ifdef CMPI
      USE MESSENGER
      USE HSWRITER, ONLY: writeHotstart_through_HSwriter  !st3 100711 for hsfile

#endif

#ifdef CSWAN
Casey 090302: We need these values from other places.
      USE OWIWIND,     ONLY: WindMultiplier
      USE Couple2Swan, ONLY: ComputeWaveDrivenForces,
     &                       CouplingInterval,
     &                       COUPWIND,
     &                       InterpoWeight,
     &                       SWAN_WX2,
     &                       SWAN_WY2
Casey 100205: Add a variable for writing of SWAN hot-start files.
     &                      ,WriteSwanHotStart
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
     &                      ,ComputeModifiedWaveFriction,
     &                       ComputeWaveFrictionProperties,
     &                       TKXX,
     &                       TKXY,
     &                       TKYY,
     &                       WAVE_A,
     &                       WAVE_A1,
     &                       WAVE_A2,
     &                       WAVE_H,
     &                       WAVE_H1,
     &                       WAVE_H2,
     &                       WAVE_T,
     &                       WAVE_T1,
     &                       WAVE_T2
#endif
#endif


      IMPLICIT NONE 
      INTEGER, intent(in) :: IT
C
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      REAL(SZ) :: WAVA
      REAL(SZ) :: WAVH
      REAL(SZ) :: WAVT
#endif
#endif

c. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
      CHARACTER(6) :: TEMPDIRNAME
      COMPLEX(SZ) :: DUDS !jgf48.50 declare size SZ instead of plain COMPLEX
      REAL(SZ) :: KSLIP

      INTEGER IE, I, J, K                  !local loop counters
      INTEGER NM1, NM2, NM3, NM123
      INTEGER NC1, NC2, NC3, NCEle, NCI
      INTEGER NCyc, NA,NITEMS
      INTEGER NMN1,NMN2,NMN3,NWETNEI,NWETADJ  ! sb v46.28.sb05.06 11/01/2006

      logical, save ::  EtaDisc_Fill = .TRUE.

      REAL(SZ) ArgT, ArgTP, ArgSAlt
      REAL(SZ) CCSFEA
      REAL(SZ) ElMax
      REAL(SZ) Fr, FricBP               !for bridge pilings
      REAL(SZ) Z0B1                     !for varying roughness
      REAL(SZ) H1
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) HTOT
      REAL(SZ) EtaN1,EtaN2,EtaN3,HTotN1,HTotN2,HTotN3               
      REAL(SZ) QTRatio
      REAL(SZ) RStRatio, RSX, RSY
      REAL(SZ) SAltMul, S2SFEA
      REAL(SZ) SFacAvg
      REAL(SZ) TPMul
      REAL(SZ) WTRatio, WindX, WindY, WindMag, WDragCo
      REAL(SZ) DARhoMRho0N1,DARhoMRho0N2,DARhoMRho0N3
      REAL(SZ) DEta2DX,DEta2DY,DRhoDX,DRhoDY
      REAL(SZ) UN1
      REAL(SZ) UV0, UV1, UV2
      REAL(SZ) VIDBCPDXOHN1,VIDBCPDXOHN2,VIDBCPDXOHN3,VIDBCPDXOHAvgArea
      REAL(SZ) VIDBCPDYOHN1,VIDBCPDYOHN2,VIDBCPDYOHN3,VIDBCPDYOHAvgArea
      REAL(SZ) WindDragLimit    !jgf46.01 Cap on max wind drag coefficient
      REAL(SZ) DPMIN      ! sb v46.28.sb05.06 11/01/2006
      REAL(SZ) PRDIFF,PRBCKGRND_MH2O     ! tcm v49.16 20100617 added

      REAL(8) AreaIE2,AreaEle
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) TimeLoc, TimeH
      INTEGER WarnElevExceeded, ErrorElevExceeded
      REAL(SZ) HollandTime, AsymmetricTime !for parametric hurricane models
      ! jgf49.1001 factor used to blend vortex and background winds for NWS29:
      REAL(SZ) bf ! blending factor, 0.0 to 1.0
C   kmd48.33bc - added in the heat flux variables
      CHARACTER(80) :: CDUM80
      INTEGER :: NumofBCNode
      REAL(SZ), SAVE :: StaTimHS, RefTimHS
      REAL(SZ), PARAMETER :: RHOAIR=1.225
      INTEGER :: NOD
      REAL(SZ), ALLOCATABLE :: TMP(:,:)
      REAL(SZ) :: CD, CDQ, QWIND

C   TCM V49.64.01 -- ADDED FOR ICE CONCENTRATION FIELDS
      REAL(SZ) PIC,CICE_TRatio  !ICE VARIABLES


C.... tcm v49.16 20100617 added
C.... convert background pressure from millibars to meters of water
      PRBCKGRND_MH2O = 100.0D0*PRBCKGRND/(RHOWAT0*G)

C    kmd48.33bc - changed for timestep changes in the hot start files
C     jgf46.21 Combined flux/radiation b.c. for rivers
#ifdef IBM
      IF ((CHOTHS).AND.(DTDP.NE.DTDPHS)) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
      ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP,KIND(0.0D0))
      END IF
#else
      IF ((CHOTHS).AND.(DTDP.NE.DTDPHS)) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
      ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP)
      END IF
#endif

C...  COMPUTE MASTER TIME WHICH IS REFERENCED TO THE BEGINNING TIME OF
C...  THE MODEL RUN
C...  
      TimeLoc=IT*DTDP + StaTim*86400.D0
C    kmd48.33bc - added for the changes in the timestep in the hot start files
      IF ((CHOTHS).AND.(DTDP.NE.DTDPHS)) THEN
         IF ((ITHS+1).EQ.IT) THEN
           StaTimHS=((IT-1)*DTDPHS)/86400.D0
           RefTimHS=((IT-1)*DTDP)/86400.D0
         END IF
         TimeLoc=IT*DTDP + (StaTimHS - RefTimHS)*86400.D0
      END IF
      
C...  HARMONIC CALCULATIONS ARE MADE FOR TIME WHICH INCLUDES THE REFTIM
C...  TO ALLOW FOR THE POSSIBILITY THAT THE EQUILIBRIUM ARGUMENTS MAY
C...  BE FOR A TIME OTHER THAN THE MODEL STARTING TIME.
C...  
      TimeH=IT*DTDP + (StaTim - RefTim)*86400.D0
!kmd - added this for the cases where the timestep changes from hot start timestep
      IF ((CHOTHS).AND.(DTDP.NE.DTDPHS)) THEN
         IF ((ITHS+1).EQ.IT) THEN
           StaTimHS=((IT-1)*DTDPHS)/86400.D0
           StaTim=((IT-1)*DTDP)/86400.D0
         END IF
         TimeH=IT*DTDP + ((StaTimHS - StaTim) - RefTim)*86400.D0
      END IF


C...  SHIFT THE FLUX PER UNIT WIDTH, DEPTH AVERAGED VELOCITIES, BOTTOM STRESS,
C...  WIND STRESS, SURFACE PRESSURE AND TIDAL POTENTIALS TO PREVIOUS TIME STEP.
C...  ZERO OUT THE NEW FORCING TERMS AND LOAD VECTORS
C...  COMPUTE A NEW BOTTOM FRICTION COEFFICIENT
C...  
ckmd	Shift values in time for predictor-corrector algorithm 
         
      DO I=1,NP

         if(CPRECOR) THEN 
            UU0(I)=UU1(I)
            VV0(I)=VV1(I)
            QX0(I)=QX1(I)
            QY0(I)=QY1(I)
         end if
         QX1(I)=QX2(I)
         QY1(I)=QY2(I)
         UU1(I)=UU2(I)
         VV1(I)=VV2(I)
         GWCE_LV(I) =0.D0
         MOM_LV_X(I)=0.D0
         MOM_LV_Y(I)=0.D0

C...  Transport
         IF(IM.EQ.10) THEN
           TRANS_LV_B(I)=0.D0
           TRANS_LV_A(I)=0.D0
           ENDIF

C...  Wind (& wave radiation stress if used)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
            WSX1(I)=WSX2(I)
            WSX2(I)=0.D0
            WSY1(I)=WSY2(I)
            WSY2(I)=0.D0
            PR1(I)=PR2(I)
!            PR2(I)=0.D0   !tcm v49.16 20100617 
            PR2(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
            ENDIF

C     TIP..Tidal potential forcing
         if(CTIP) then
            TIP1(I)=TIP2(I)
            TIP2(I)=0.D0
         endif

      END DO

C     2DDI.For the 2DDI version of the code
      if(C2DDI) then
C     2DDI.Set up the 2D friction coefficient
         CALL Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, IFNLFA,
     &        NP, TK) 
C
      endif 
C..RJW 3D for the 3D version of the code
      if(C3D) then
C     3D.Set up the 3D friction coefficient
         CALL Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, IFNLFA,
     &        NP, TK, NFEN)
C
      endif

C
C...  SHIFT THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION TO PREVIOUS
C...  TIME STEPS.  ZERO OUT THE NEW SPECIFIED NORMAL FLOW BOUNDARY
C...  CONDITION
C...  
      DO I=1,NVEL
         QN0(I)=QN1(I)
         QN1(I)=QN2(I)
         QN2(I)=0.D0
         EN0(I)=EN1(I)
         EN1(I)=EN2(I)
         EN2(I)=0.D0 
      END DO
C...  
C...  DEFINE Ramp FUNCTION FOR BOUNDARY ELEVATION FORCING, WIND AND PRESSURE
C.... FORCING AND TIDAL POTENTIAL FORCING
C...  
C
C     jgf46.08 Calculate ramp functions.
C     jgf46.21 Modify to match behavior of 46.02
      SELECT CASE(NRamp)
      CASE(0)
         Ramp=1.0D0
         RampExtFlux=1.0D0
         RampIntFlux=1.0D0
         RampElev=1.0D0
         RampTip=1.0D0
         RampMete=1.0D0
         RampWRad=1.0D0
      CASE(1,2,3,4,5,6,7)
         Ramp=TANH((2.D0*IT*DTDP/86400.D0)/DRamp)
         RampExtFlux=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampExtFlux)
         RampIntFlux=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampIntFlux)
         RampElev=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampElev)
         RampTip=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampTip)
         RampMete=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampMete)
         RampWRad=TANH((2.D0*(IT)*DTDP/86400.D0)/DRampWRad)
      END SELECT
C    kmd48.33bc - ramp changes with baroclinic when timestep is changed
      IF ((CHOTHS).AND.(DTDP.NE.DTDPHS).AND.(NRAMP.NE.0)) THEN
         Ramp=TANH((2.D0*TimeLoc/86400.D0)/DRamp)
         RampExtFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampExtFlux)
         RampIntFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampIntFlux)
         RampElev=TANH((2.D0*TimeLoc/86400.D0)/DRampElev)
         RampTip=TANH((2.D0*TimeLoc/86400.D0)/DRampTip)
         RampMete=TANH((2.D0*TimeLoc/86400.D0)/DRampMete)
         RampWRad=TANH((2.D0*TimeLoc/86400.D0)/DRampWRad)
      END IF

C     jgf46.21 If there is an external flux (i.e. river) boundary, turn
C     off all forcings except the river flux forcing for the duration of
C     the FluxSettlingTime. When the FluxSettlingTime has ended, turn
C     all forcings back on.
      IF (NRamp.gt.1) THEN
         IF(IT.LT.(FluxSettlingIT+10)) THEN
            Ramp=0.0
            RampIntFlux=0.0
            RampElev=0.0
            RampTip=0.0
            RampMete=0.0
            RampWRad=0.0
         ELSE
            Ramp=TANH((2.D0*(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRamp)
            RampIntFlux=TANH((2.D0
     &           *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampIntFlux)
            RampElev=TANH((2.D0
     &           *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampElev)
            RampTip=TANH((2.D0
     &           *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampTip)
            RampMete=TANH((2.D0
     &           *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampMete)
            RampWRad=TANH((2.D0
     &           *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampWRad)
         ENDIF
         !jgf49.44: Cover the case where the ramp length is zero.
         IF (DRamp.lt.1.0e-6) Ramp = 1.0d0
         IF (DRampExtFlux.lt.1.0e-6) RampExtFlux = 1.0d0
         IF (DRampIntFlux.lt.1.0e-6) RampIntFlux = 1.0d0
         IF (DRampElev.lt.1.0e-6) RampElev = 1.0d0
         IF (DRampTip.lt.1.0e-6) RampTip = 1.0d0
         IF (DRampMete.lt.1.0e-6) RampMete = 1.0d0
         IF (DRampWRad.lt.1.0e-6) RampWRad = 1.0d0
      ELSE
         !jgf49.44: Cover the case where the ramp length is zero.
         IF (DRamp.lt.1.0e-6) Ramp = 1.0d0
      ENDIF

C
C------------------------ICE FIELDS----------------------------------------
C...  UPDATE THE ICE CONCENTRATION FIELDS FROM UNIT 25,225,227
c...  TCM V49.64.01 ADDED THE ICE FIELDS SECTION
C
      IF (NCICE.EQ.12) THEN 
         IF(TimeLoc.GT.CICE_TIME2) THEN
            CICE_TIME1 = CICE_TIME2
            CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
            DO I=1,NP
              CICE1(I)=CICE2(I)
            END DO
            CALL NCICE1_GET(CICE2,NP)
          ENDIF
      ENDIF


C------------------------MET FORCING---------------------------------------
C
C...  UPDATE THE WIND STRESS AND SURFACE PRESSURE AND READ IN NEW VALUES
C...  FROM UNIT 22.  APPLY Ramp FUNCTION.
C
C  tcm v49.16 20100617 -- Changed pressure ramping so that we apply the ramp
C     to the difference between the background pressure and the forced pressure,
C     then add that ramped difference back to the background pressure

Casey 091217: Lower the wind drag limit from 0.0035 to 0.0020.
      WindDragLimit = 0.0035
C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         DO I=1,NP
            READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I) !read in
            WSX2(I)=RampMete*WSX2(I) ! apply ramp function
            WSY2(I)=RampMete*WSY2(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete*PR2(I)
            PRDIFF = RampMete*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) ! for met recording station output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(ABS(NWS).EQ.2) THEN
         IF(TimeLoc.GT.WTIME2) THEN !determine if met file time incr. is exceeded
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I) ! move current data to old
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
               READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I) ! read in 
            END DO
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC ! interpolate
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete*WindX !apply ramp
            WSY2(I) = RampMete*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED,
     &           NP,NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &           NSCREEN, ScreenUnit)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(ABS(NWS).EQ.4) THEN 
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(ABS(NWS).EQ.5) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
               READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
            END DO
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            NWSGGWI=NWSGGWI+1
            CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &           WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
      IF(ABS(NWS).EQ.7) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &           WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete*WindX !apply ramp
            WSY2(I) = RampMete*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate 
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         HollandTime = TimeLoc
         CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         ENDDO
      ENDIF       

C     jgf46.16 Merged:
!     cf & cm added nws = 9: asymmetric hurricane winds
      AsymmetricTime = TimeLoc
      IF(NWS.EQ.9) THEN
         CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,AsymmetricTime,ICS)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            NWSGGWI=NWSGGWI+1
            CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,
     &           G,NWLON,NWLAT,WTIMINC)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            IDSETFLG=IDSETFLG+1
            IF(IDSETFLG.GT.8) THEN
               NWSEGWI=NWSEGWI+1
               IDSETFLG=1
            ENDIF
            CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,
     &           NP,RHOWAT0,G)
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF

C...sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
        IF(ABS(NWS).EQ.12) THEN
          IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              PRN1(I)=PRN2(I)
            END DO
            CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          ENDIF
          WTRATIO=(TimeLoc-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &             WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY) 
              WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF 
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX/WindMultiplier
               SWAN_WY2(I,2) = WindY/WindMultiplier
            ENDIF
#endif
          END DO
        ENDIF

C     RJW ffpl Merged:
!     rjw added nws = 19: asymmetric hurricane winds
      AsymmetricTime = TimeLoc
      IF(NWS.EQ.19) THEN
         CALL NWS19GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,AsymmetricTime,ICS)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF

            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
               CICEOUT(I) = PIC
            ENDIF
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
            IF(COUPWIND)THEN
               SWAN_WX2(I,2) = WindX
               SWAN_WY2(I,2) = WindY
            ENDIF
#endif
         END DO
      ENDIF
C
C     jgf49.1001 Added NWS29 for embedding an asymmetric vortex from
C     NWS19 into an OWI basin scale met field from NWS12 (derived from NAM).
      IF(ABS(NWS).EQ.29) THEN
         ! bring in next set of OWI met data 
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         ENDIF
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         ! bring in next set of asymmetric vortex met data to separate arrays
         CALL NWS19GET(SLAM,SFEA,vortexWVNX2,vortexWVNY2,vortexPRN2,
     &                 NP,AsymmetricTime, ICS)
         DO I=1,NP
            ! compute wind stress due to background met at this node
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ! jgf49.1001 NAM winds already contain wind reduction
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I))) 
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
            ! compute wind stress due to vortex met at this node
            WindX = vortexWVNX2(I)
            WindY = vortexWVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            IF (LoadCanopyCoef) 
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            CALL getBlendFactor(I, SLAM, SFEA, bf)
            ! blend the wind stresses and barometric pressures
            WSX2(I) = bf*RampMete*0.001293d0*WDragCo*WindX*WindMag
     &               +(1.0d0-bf)*WSX2(I)
            WSY2(I) = bf*RampMete*0.001293d0*WDragCo*WindY*WindMag
     &               +(1.0d0-bf)*WSY2(I)
            PR2(I)= bf*RampMete*vortexPRN2(I)
     &               +(1.0d0-bf)*PRN2(I)
            WVNXOUT(I)=bf*RampMete*WindX+(1.0d0-bf)*WVNXOUT(I)
            WVNYOUT(I)=bf*RampMete*WindY+(1.0d0-bf)*WVNYOUT(I)
         END DO
      ENDIF


C    kmd48.33bc - added in information for the sponge layer
       IF ((C3D).AND.(RES_BC_FLAG.GT.0).AND.(NWS.NE.0)) THEN
         DO I=1,NP
            WSX2(I)=sponge(I)*WSX2(I)
            WSY2(I)=sponge(I)*WSY2(I)
         END DO
       END IF

C
C--------------------END MET FORCING---------------------------------------



C...  UPDATE THE WAVE RADIATION STRESS AND READ IN NEW VALUES FROM
C.... UNIT 23.  APPLY Ramp FUNCTION.  ADD RADIATION STRESS TO WIND
C...  STRESS
C...  
C...  NRS=2 was added.  sb46.28sb03 09/xx/2006
C...  TCM v49.48 Restructured the wave stress updates in order to 
c...             include NRS=4
      IF(NRS.NE.0) THEN
         if((NRS.eq.1).or.(nrs.eq.2).or.(nrs.eq.3)) then
            IF(TimeLoc.GT.RSTIME2) THEN
               RSTIME1=RSTIME2
               RSTIME2=RSTIME2+RSTIMINC
               DO I=1,NP
                  RSNX1(I)=RSNX2(I)
                  RSNY1(I)=RSNY2(I)
               END DO
               IF(NRS.EQ.1) THEN
                 CALL RSGET(RSNX2,RSNY2,NP)
               ENDIF
               IF(NRS.EQ.2) THEN
                 CALL RS2GET(RSNX2,RSNY2,NP)
               ENDIF
#ifdef CSWAN
Casey 090302: Added for coupling to SWAN.
               IF(NRS.EQ.3) THEN
                 InterpoWeight = 1.0
                 CALL ComputeWaveDrivenForces
Casey 090707: We want to extrapolate forward in time.  Load the latest (current) forces
C             into RSNX1/RSNY1, and then load the future forces into RSNX2/RSNY2.
                 DO I=1,NP
                   RSX = RSNX1(I)
                   RSY = RSNY1(I)
                   RSNX1(I) = RSNX2(I)
                   RSNY1(I) = RSNY2(I)
                   RSNX2(I) = RSNX2(I) + (RSNX2(I)-RSX)
                   RSNY2(I) = RSNY2(I) + (RSNY2(I)-RSY)
                 ENDDO
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.  Compute the wave properties
C             for height, period, angle and dissipation.  Then extrapolate forward
C             in time using the same logic as above.
                 DO I=1,NP
                   WAVE_A1(I) = WAVE_A2(I)
                   WAVE_H1(I) = WAVE_H2(I)
                   WAVE_T1(I) = WAVE_T2(I)
                 ENDDO
                 CALL ComputeWaveFrictionProperties
                 DO I=1,NP
                   WAVA = WAVE_A1(I)
                   WAVH = WAVE_H1(I)
                   WAVT = WAVE_T1(I)
                   WAVE_A1(I) = WAVE_A2(I)
                   WAVE_H1(I) = WAVE_H2(I)
                   WAVE_T1(I) = WAVE_T2(I)
                   WAVE_A2(I) = WAVE_A2(I) + (WAVE_A2(I)-WAVA)
                   WAVE_H2(I) = WAVE_H2(I) + (WAVE_H2(I)-WAVH)
                   WAVE_T2(I) = WAVE_T2(I) + (WAVE_T2(I)-WAVT)
                 ENDDO
#endif
               ENDIF
#endif
            ENDIF
            RStRatio=(TimeLoc-RSTIME1)/RSTIMINC
            DO I=1,NP
               RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I)-RSNX1(I)))
               RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I)-RSNY1(I)))
               WSX2(I) = WSX2(I) + RSX
               WSY2(I) = WSY2(I) + RSY
#ifdef CSWAN
Casey 090302: Added these lines for output to the rads.64 file.
               RSNXOUT(I) = RSX
               RSNYOUT(I) = RSY
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
               WAVE_A(I) = WAVE_A1(I) + RSTRATIO*(WAVE_A2(I)-WAVE_A1(I))
               WAVE_H(I) = WAVE_H1(I) + RSTRATIO*(WAVE_H2(I)-WAVE_H1(I))
               WAVE_T(I) = WAVE_T1(I) + RSTRATIO*(WAVE_T2(I)-WAVE_T1(I))
#endif
#endif
            ENDDO
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            CALL ComputeModifiedWaveFriction(TK)
#endif
#endif
         ENDIF  !nrs = 1,2,or 3

         ! Tightly Coupled Code with STWAVE
         ! Apply the ramping function and add wave stress to WSX2,WSY2
         ! This cases uses a step function in time
         IF(NRS == 4) THEN    ! vjp modified Jan 22 2010
            IF(TimeLoc.GE.RSTIME2) THEN   !Get a new wave record
               RSTIME1=RSTIME2
               RSTIME2=RSTIME2+RSTIMINC
            ENDIF
            IF (TimeLoc > ENDWAVE+RSTIMINC) THEN
              RSNX2(:) = 0.0d0;  RSNY2(:) = 0.0d0
            ENDIF
            DO I=1,NP
               RSX = RampWRad*RSNX2(I)
               RSY = RampWRad*RSNY2(I)
               WSX2(I) = WSX2(I) + RSX
               WSY2(I) = WSY2(I) + RSY
            END DO
         ENDIF !(NRS = 4)
      
      ENDIF  !(End test for updating wave radiation stress)

C     jgf48.4627 Skip past GWCE and momentum calculations if only
C     meteorological output was requested.
      IF (METONLY) THEN 
         goto 9999
      ENDIF

C   kmd48.33 - added in information for the elevation boundary conditions
C              used in prognostic runs. The diagnostic information is read 
C              in once and used during the simulation. Note that there is 
C              no ramp utilized for this boundary condition. 
      IF((C3D).AND.(RES_BC_FLAG.GT.0).AND.(CBAROCLINIC)) THEN
         IF (ABS(RES_BC_FLAG).GE.1) THEN
            IF(TimeLoc.GT.RBCTIME2) THEN
               RBCTIME1=RBCTIME2
               RBCTIME2=RBCTIME2+RBCTIMEINC
               READ(35,'(A)') CDUM80
               DO I=1,NETA
                 LNM_BC1(I)=LNM_BC2(I)
                 READ(35,*) NOD,LNM_BC2(I)
               END DO
            END IF
            RBCRATIO=(TimeLoc-RBCTIME1)/RBCTIMEINC
            DO NumofBCNode=1,NETA
                  LNM_BC(NumofBCNode)=LNM_BC1(NumofBCNode)+
     &            RBCRATIO*(LNM_BC2(NumofBCNode)-
     &            LNM_BC1(NumofBCNode))
            END DO
         END IF
      END IF

C...  
C...  Tidal Potential Forcing
C...  Note, the Earth tide potential reduction factor, ETRF(J) has been
C...        incorporated into this calculation.
C... 
      IF(CTIP) THEN
         DO J=1,NTIF
            IF(PERT(J).EQ.0.) THEN
               NCYC=0
            ELSE
#ifdef IBM
               NCYC=INT(timeh/PERT(J),KIND(0.0d0))
#else
               NCYC=INT(timeh/PERT(J))
#endif
            ENDIF
            ARGT=AMIGT(J)*(timeh-NCYC*PERT(J))+FACET(J)
            TPMUL=RampTip*ETRF(J)*TPK(J)*FFT(J)
            SALTMUL=RampTip*FFT(J)
#ifdef IBM
            NA=NINT(0.00014/AMIGT(J),KIND(0.0d0))
#else
            NA=NINT(0.00014/AMIGT(J))
#endif
            IF(NA.EQ.1) THEN    !SEMI-DIURNAL SPECIES
               DO I=1,NP
                  ARGTP=ARGT+2.d0*SLAM(I)
                  ARGSALT=ARGT-SALTPHA(J,I)
                  CCSFEA=COS(SFEA(I))
                  CCSFEA=CCSFEA*CCSFEA
                  TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
            ENDIF
            IF(NA.EQ.2) THEN    !DIURNAL SPECIES
               DO I=1,NP
                  ARGTP=ARGT+SLAM(I)
                  ARGSALT=ARGT-SALTPHA(J,I)
#ifdef REAL4
                  S2SFEA=SIN(2.e0*SFEA(I))
#else 
                  S2SFEA=SIN(2.d0*SFEA(I))
#endif
                  TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
            ENDIF
         END DO
      ENDIF

C...
C...  Depth Averaged Baroclinic Forcing needed by GWCE and 2DDI Momentum
C...  Compute this (divided by H, i.e., Bx/H, By/H) as a nodally averaged
C...  quantity for smoothing
C...
      DO J=1,NP
         VIDBCPDXOH(J)=0.D0
         VIDBCPDYOH(J)=0.D0
cjgf45.06 TotalArea(J)=0.D0 ! What if CBaroclinic=F and nolifa=0?
      ENDDO
  
      IF(CBaroclinic) Then

         DO J=1,NP              !jgf45.06 try this instead
            TotalArea(J)=0.D0   !jgf45.06
         ENDDO                  !jgf45.06

!     Kendra45.12 - Test placement of BPG 
         if (C3DVS) CALL BPG3D()
         
c     jgf45.12 this algorithm only works in 2D
         if (C2DDI) then 
            DO IE=1,NE
               NM1=NM(IE,1)
               NM2=NM(IE,2)
               NM3=NM(IE,3)
               NC1=NODECODE(NM1)
               NC2=NODECODE(NM2)
               NC3=NODECODE(NM3)
               NCEle=NC1*NC2*NC3*NOFF(IE)
               H2N1=DP(NM1)+IFNLFA*ETA2(NM1) !jgf45.11 add IFNLFA (kd bug fix)
               H2N2=DP(NM2)+IFNLFA*ETA2(NM2) ! "
               H2N3=DP(NM3)+IFNLFA*ETA2(NM3) ! "
               EtaN1=IFNLFA*Eta2(NM1)
               EtaN2=IFNLFA*Eta2(NM2)
               EtaN3=IFNLFA*Eta2(NM3)
               SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
               
               AreaIE2=Areas(IE)
               AreaEle=NCEle*AreaIE2/2.D0
               FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1=X(NM3)-X(NM2) !a1
               FDY2=X(NM1)-X(NM3) !a2
               FDY3=X(NM2)-X(NM1) !a3
               FDX1O2A=FDX1/AreaIE2 !dphi1/dx
               FDY1O2A=FDY1/AreaIE2 !dphi1/dy
               FDX2O2A=FDX2/AreaIE2 !dphi2/dx
               FDY2O2A=FDY2/AreaIE2 !dphi2/dy
               FDX3O2A=FDX3/AreaIE2 !dphi3/dx
               FDY3O2A=FDY3/AreaIE2 !dphi3/dy
               
               DARhoMRho0N1=(DASigT(NM1)-SigT0)/RhoWat0
               DARhoMRho0N2=(DASigT(NM2)-SigT0)/RhoWat0
               DARhoMRho0N3=(DASigT(NM3)-SigT0)/RhoWat0
               DEta2DX=EtaN1*FDX1O2A+EtaN2*FDX2O2A+EtaN3*FDX3O2A
               DEta2DY=EtaN1*FDY1O2A+EtaN2*FDY2O2A+EtaN3*FDY3O2A
               DRhoDX=DARhoMRho0N1*FDX1O2A+DARhoMRho0N2*FDX2O2A
     &              +DARhoMRho0N3*FDX3O2A
               DRhoDY=DARhoMRho0N1*FDY1O2A+DARhoMRho0N2*FDY2O2A
     &              +DARhoMRho0N3*FDY3O2A
               VIDBCPDXOHN1=Ramp*G*
     &              (DARhoMRho0N1*DEta2DX+0.5d0*H2N1*DRhoDX)
               VIDBCPDXOHN2=Ramp*G*
     &              (DARhoMRho0N2*DEta2DX+0.5d0*H2N2*DRhoDX)
               VIDBCPDXOHN3=Ramp*G*
     &              (DARhoMRho0N3*DEta2DX+0.5d0*H2N3*DRhoDX)
               VIDBCPDYOHN1=Ramp*G*
     &              (DARhoMRho0N1*DEta2DY+0.5d0*H2N1*DRhoDY)
               VIDBCPDYOHN2=Ramp*G*
     &              (DARhoMRho0N2*DEta2DY+0.5d0*H2N2*DRhoDY)
               VIDBCPDYOHN3=Ramp*G*
     &              (DARhoMRho0N3*DEta2DY+0.5d0*H2N3*DRhoDY)
               VIDBCPDXOHAvgArea=AreaEle*(VIDBCPDXOHN1+VIDBCPDXOHN2
     &              +VIDBCPDXOHN3)/3.D0
               VIDBCPDYOHAvgArea=AreaEle*(VIDBCPDYOHN1+VIDBCPDYOHN2
     &              +VIDBCPDYOHN3)/3.D0
               VIDBCPDXOH(NM1)=VIDBCPDXOH(NM1)+VIDBCPDXOHAvgArea
               VIDBCPDXOH(NM2)=VIDBCPDXOH(NM2)+VIDBCPDXOHAvgArea
               VIDBCPDXOH(NM3)=VIDBCPDXOH(NM3)+VIDBCPDXOHAvgArea
               VIDBCPDYOH(NM1)=VIDBCPDYOH(NM1)+VIDBCPDYOHAvgArea
               VIDBCPDYOH(NM2)=VIDBCPDYOH(NM2)+VIDBCPDYOHAvgArea
               VIDBCPDYOH(NM3)=VIDBCPDYOH(NM3)+VIDBCPDYOHAvgArea
               TotalArea(NM1)=TotalArea(NM1)+AreaEle
               TotalArea(NM2)=TotalArea(NM2)+AreaEle
               TotalArea(NM3)=TotalArea(NM3)+AreaEle
            ENDDO
            
            DO J=1,NP
               IF(TotalArea(J).NE.0.) THEN
                  VIDBCPDXOH(J)=VIDBCPDXOH(J)/TotalArea(J)
                  VIDBCPDYOH(J)=VIDBCPDYOH(J)/TotalArea(J)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

C...  
C...  COMPUTE SPECIFIED NORMAL FLOW BOUNDARY CONDITION
C...  
      IF(NFLUXF.EQ.1) THEN
         DO J=1,NFFR
            IF(FPER(J).EQ.0.D0) THEN
               NCYC=0.
            ELSE
#ifdef IBM
               NCYC=INT(timeh/FPER(J),KIND(0.0d0))
#else
               NCYC=INT(timeh/FPER(J))
#endif
            ENDIF
            ARGJ=FAMIG(J)*(timeh-NCYC*FPER(J))+FFACE(J)
            RFF=FFF(J)*RampExtFlux  !jgf46.02 use river ramp for Katrina
            DO I=1,NVEL
               ARG=ARGJ-QNPH(J,I)
               QN2(I)=QN2(I)+QNAM(J,I)*RFF*COS(ARG)
               IF(LBCODEI(J).EQ.32) THEN
                  ARG=ARGJ-ENPH(J,I)
                  EN2(I)=EN2(I)+ENAM(J,I)*RFF*COS(ARG)
               ENDIF
            END DO
         END DO
         IF(NFFR.EQ.0) THEN
            IF(TimeLoc.GT.QTIME2) THEN
               QTIME1=QTIME2
               QTIME2=QTIME2+FTIMINC
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     QNIN1(J)=QNIN2(J)
                     READ(20,*) QNIN2(J)
                  ENDIF
                  IF(LBCODEI(J).EQ.32) THEN
                     QNIN1(J)=QNIN2(J)
                     ENIN1(J)=ENIN2(J)
                     READ(20,*) QNIN2(J),ENIN2(J)
                  ENDIF
               END DO
            ENDIF
            QTRATIO=(TimeLoc-QTIME1)/FTIMINC
            DO I=1,NVEL
               QN2(I)=RampExtFlux*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
               EN2(I)=RampExtFlux*(ENIN1(I)+QTRATIO*(ENIN2(I)-ENIN1(I)))
            END DO
         ENDIF

CRAL_add_42.06f
C     jgf46.21 Collect elevation information for river radiation b.c.
         IF(IT.EQ.FluxSettlingIT) THEN
            EtaDisc_Fill = .FALSE.    ! sb v46.48 11/06/2006
            DO I=1, NP
               EtaDisc(I) = Eta2(I)   ! EtaDisc written to hotstart file
            ENDDO
            DO I=1,NVEL
               IF(LBCODEI(I).EQ.52) THEN
                  NNBB=NBV(I)
                  ElevDisc(I)=Eta2(NNBB)
               ENDIF
            END DO
         ELSE IF(EtaDisc_Fill .and. IT > FluxSettlingIT) THEN
            EtaDisc_Fill = .FALSE.
            DO I=1,NVEL
               IF(LBCODEI(I).EQ.52) THEN
                  NNBB=NBV(I)
                  ElevDisc(I)=EtaDisc(NNBB)   ! sb v46.48 11/06/2006
               ENDIF
            END DO
         ENDIF

      ENDIF
C...  
C...  COMPUTE DISCHARGE CONTRIBUTION FROM RADIATION BOUNDARY CONDITION
C...  
      IF(NFLUXRBC.EQ.1) THEN
         DO J=1,NVEL
            IF(LBCODEI(J).EQ.30) THEN
               NNBB=NBV(J)
               H1=DP(NNBB)+IFNLFA*ETA2(NNBB)
               UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
               QN1(J)=H1*UN1
            ENDIF
         END DO
      ENDIF

C...  COMPUTE DISCHARGE CONTRIBUTION FROM ZERO NORMAL VELOCITY GRADIENT
C...  BOUNDARY CONDITION
C...  
        IF(NFLUXGBC.EQ.1) THEN
          DO J=1,NVEL
            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) THEN
              NNBB=NBV(J)
              H1=DP(NNBB)+IFNLFA*ETA2(NNBB)
              UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
              QN1(J)=H1*UN1
              ENDIF
            END DO
          ENDIF
C...  
C...  COMPUTE SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C.... EXTERNAL BARRIER BOUNDARY NODES
C...  
        IF(NFLUXB.EQ.1) THEN
          DO I=1,NVEL
            IF((LBCODEI(I).EQ.3).OR.(LBCODEI(I).EQ.13).OR.
     &        (LBCODEI(I).EQ.23)) THEN
              NNBB=NBV(I)
              RBARWL=2.D0*(ETA2(NNBB)-BARLANHT(I))/3.D0
              IF(RBARWL.GT.BARMIN) THEN
                QN2(I)=-RampIntFlux
     &                *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
              ELSE
                QN2(I)=0.0D0
              ENDIF
            ENDIF
          END DO
        ENDIF

C...  COMPUTE INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED INTERNAL BARRIER
C...  BOUNDARY (PERMEABLE OR NOT) NODES
C...  
C     jgf46.03 Begin block of notes for internal barrier boundaries
C     
C     NFLUXIB is set to 1 in read_input.F if there are internal barrier
C     boundaries in the fort.14 (mesh) file.
C
C     IBSTART is a flag that indicates the first time through the time
C     stepping loop; set to 0 in read_input.F and set to 1 here.
C
C     BARAVGWT was apparently intended for use in averaging internal
C     barrier water levels. It is set to 0 in read_input.F, which seems
C     to turn off any time averaging here.
C
C     NIBNODECODE seems to be set to 1 for nodes receiving water across
C     the barrier
C 
C     BARMIN is used in several places, mainly as the minimum elevation
C     above the levee for flow to occur. It is a parameter and is set to
C     0.04 in global.F.
C
C
      IF(NFLUXIB.EQ.1) THEN
         DO I=1,NP
            NIBNODECODE(I)=0
         END DO
         I = 0                          !Seizo 05/06/2008
         NBOU_LOOP_K: DO K = 1, NBOU    !Seizo 05/06/2008
            IE = 1                      !Seizo 05/06/2008
            IF ( (LBCODEI(I+1)==4) .OR. (LBCODEI(I+1)==24)        !Seizo 05/06/2008
     &      .OR. (LBCODEI(I+1)==5) .OR. (LBCODEI(I+1)==25) ) then !Seizo 05/06/2008
               IE = 2                   !Seizo 05/06/2008
            ENDIF                       !Seizo 05/06/2008
            IF( IE /= 2 ) then          !Seizo 05/06/2008
               I = I + NVELL(K)         !Seizo 05/06/2008
            ELSE                        !Seizo 05/06/2008
            NVELL_LOOP_J: DO J = 1, NVELL(K) * IE   !Seizo 05/06/2008
               I = I + 1                            !Seizo 05/06/2008
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
#if 1
               NNBB1WN = 0      ! COUNT NUMBER OF WET NEIGHBORS sb46.28.sb05.05 11/01/2006
               NNBB2WN = 0      ! COUNT NUMBER OF WET NEIGHBORS sb46.28.sb05.05 11/01/2006
               IF( (J==1).OR.(J==NVELL(K)+1) ) THEN !Seizo 05/06/2008
                            NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))     ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))  ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
               ELSE                                 !Seizo 05/06/2008
                    IF( (J==NVELL(K)).OR.(J==NVELL(K)*IE) ) THEN       !                            !Seizo 05/06/2008
                            NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))     ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))  ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                    ELSE                            !Seizo 05/06/2008
                            NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))     ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                            NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))     ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))  ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))  ! sb46.28.sb05.05 11/01/2006 !Seizo 05/06/2008
                    ENDIF                           !Seizo 05/06/2008
               ENDIF                                !Seizo 05/06/2008
#endif
               IF(IBSTART.EQ.0)THEN
                  RBARWL1AVG(I)=ETA2(NNBB1)-BARINHT(I)
                  RBARWL2AVG(I)=ETA2(NNBB2)-BARINHT(I)
                  IBSTART=1
               ELSE
                  RBARWL1AVG(I)=(ETA2(NNBB1)-BARINHT(I)+BARAVGWT
     &                 *RBARWL1AVG(I))/(1.D0+BARAVGWT)
                  RBARWL2AVG(I)=(ETA2(NNBB2)-BARINHT(I)+BARAVGWT
     &                 *RBARWL2AVG(I))/(1.D0+BARAVGWT)
               ENDIF
               RBARWL1=RBARWL1AVG(I)
               RBARWL2=RBARWL2AVG(I)
               RBARWL1F=2.D0*RBARWL1/3.D0
               RBARWL2F=2.D0*RBARWL2/3.D0
               QN2(I)=0.D0
               IF((RBARWL1.LT.0.D0).AND.(RBARWL2.LT.0.D0)) THEN
C...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
                  QN2(I)=0.D0
                  GOTO 1034
               ENDIF
               IF(ABS(RBARWL1-RBARWL2).LT.0.01D0) THEN
cjjw               IF(ABS(RBARWL1-RBARWL2).LT.BARMIN) THEN
C...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER 
C................TO WITHIN TOLERANCE BARMIN -> CASE 2
                  QN2(I)=0.D0
                  GOTO 1034
               ENDIF
               IF((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
C...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
C................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
C................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
                  IF(RBARWL2.GT.RBARWL1F) THEN
C..................OUTWARD SUBCRITICAL FLOW -> CASE 3
#if 1
                     IF(NODECODE(NNBB1).EQ.0.OR.NNBB1WN.EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#else
                     IF(NODECODE(NNBB1).EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#endif
                        QN2(I)=0.0D0
                     ELSE
                        QN2(I)=-RampIntFlux*RBARWL2*BARINCFSB(I)*
     &                     (2.D0*G*(RBARWL1-RBARWL2))**0.5D0
                     ENDIF
                  ELSE
C..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
#if 1
                     IF(NODECODE(NNBB1).EQ.0.OR.NNBB1WN.EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#else
                     IF(NODECODE(NNBB1).EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#endif
                        QN2(I)=0.0D0
                     ELSE
                        QN2(I)=-RampIntFlux*BARINCFSP(I)*RBARWL1F*
     &                    (RBARWL1F*G)**0.5D0
                     ENDIF
                  ENDIF
                  GOTO 1034
               ENDIF
               IF((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN
C...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
C................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
C................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
                  IF(RBARWL1.GT.RBARWL2F) THEN
C..................INWARD SUBCRITICAL FLOW -> CASE 5
#if 1
                     IF(NODECODE(NNBB2).EQ.0.OR.NNBB2WN.EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#else
                     IF(NODECODE(NNBB2).EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#endif
                        QN2(I)=0.0D0
                     ELSE
                        QN2(I)=RampIntFlux*RBARWL1*BARINCFSB(I)*
     &                     (2.0D0*G*(RBARWL2-RBARWL1))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ELSE
C..................INWARD SUPERCRITICAL FLOW -> CASE 6
#if 1
                     IF(NODECODE(NNBB2).EQ.0.OR.NNBB2WN.EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#else
                     IF(NODECODE(NNBB2).EQ.0) THEN ! sb46.28.sb05.05 11/01/2006
#endif
                        QN2(I)=0.0D0
                     ELSE
                        QN2(I)=RampIntFlux*BARINCFSP(I)*RBARWL2F*
     &                     (RBARWL2F*G)**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ENDIF
                  GOTO 1034
               ENDIF
 1034          CONTINUE
               END DO NVELL_LOOP_J
            ENDIF
         END DO NBOU_LOOP_K
      ENDIF
   
C...  
C...  COMPUTE INWARD/OUTWARD NORMAL FLOW FOR INTERNAL BARRIER 
C.... BOUNDARY NODES THROUGH CROSS BARRIER PIPES
C.... NOTE THAT THIS ADDS AN ADDITIONAL FLOW COMPONENT INTO QN2
C...  
      IF(NFLUXIBP.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
               IF(IBSTART.EQ.0)THEN
                  RPIPEWL1AVG(I)=ETA2(NNBB1)-PIPEHT(I)
                  RPIPEWL2AVG(I)=ETA2(NNBB2)-PIPEHT(I)
                  IBSTART=1
               ELSE
                  RPIPEWL1AVG(I)=(ETA2(NNBB1)-PIPEHT(I)+BARAVGWT
     &                 *RPIPEWL1AVG(I))/(1.D0+BARAVGWT)
                  RPIPEWL2AVG(I)=(ETA2(NNBB2)-PIPEHT(I)+BARAVGWT
     &                 *RPIPEWL2AVG(I))/(1.D0+BARAVGWT)
               ENDIF
               RBARWL1=RPIPEWL1AVG(I)
               RBARWL2=RPIPEWL2AVG(I)
               IF((RBARWL1.LT.0.D0).AND.(RBARWL2.LT.0.D0)) THEN
C...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW PIPE -> CASE 1
                  QN2(I)=QN2(I)+0.D0
                  GOTO 1036
               ENDIF
               IF(ABS(RBARWL1-RBARWL2).LT.BARMIN) THEN
C...............WATER LEVEL EQUAL ON BOTH SIDES OF PIPE -> CASE 2
                  QN2(I)=QN2(I)+0.D0
                  GOTO 1036
               ENDIF
               IF((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
C...............WATER LEVEL GREATER ON THIS SIDE OF THE PIPE AND IS SUCH
C................THAT OUTWARD DISCHARGE IS OCCURING
C................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW THROUGH THE PIPE  
C................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE PIPE TO
C................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
C................THE PIPE
                  IF(RBARWL2.LE.0.D0) THEN
C..................OUTWARD FREE DISCHARGE -> CASE 3
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        QN2(I)=QN2(I)+0.0D0
                     ELSE
                        QN2(I)=QN2(I)-RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(I))**2
     &                       *(2.D0*G*RBARWL1/(1.D0+PIPECOEF(I)))**0.5D0
                     ENDIF
                  ELSE
C..................OUTWARD SUBMERGED DISCHARGE -> CASE 4
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        QN2(I)=QN2(I)+0.0D0
                     ELSE
                        QN2(I)=QN2(I)-RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     ENDIF
                     GOTO 1036
                  ENDIF
               ENDIF
               IF((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN
C...............WATER LEVEL LOWER ON THIS SIDE OF PIPE AND IS SUCH
C................THAT INWARD DISCHARGE IS OCCURING
C................THUS THIS IS THE RECEIVING SIDE OF THE FLOW THROUGH THE  PIPE 
C................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE PIPE TO
C................REMAIN WET WHEN THERE IS FLOW ACROSS THE PIPE.
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW THROUGH 
C................THE PIPE
                  IF(RBARWL1.LE.0) THEN
C..................INWARD FREE DISCHARGE -> CASE 5
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        QN2(I)=QN2(I)+0.0D0
                     ELSE
                        QN2(I)=QN2(I)+RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(I))**2
     &                       *(2.D0*G*RBARWL2/(1.D0+PIPECOEF(I)))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ELSE
C..................INWARD SUBMERGED DISCHARGE -> CASE 6
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        QN2(I)=0.0D0
                     ELSE
                        QN2(I)=QN2(I)+RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(I))**2
     &                   *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ENDIF
                  GOTO 1036
               ENDIF
 1036          CONTINUE
            ENDIF
         END DO
      ENDIF

C...
C...  Compute the water surface elevation from the GWCE form of the 2D continuity eq.
C...

ckmd	  Changed to include the predictor-corrector algorithm
      IF(CPRECOR) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
         CALL Mom_Eqs_New_NC()
C...  If running in parallel, update velocities & fluxes on all processors
         
#ifdef CMPI
         CALL UPDATER(UU2,VV2,DUMY1,2)
         CALL UPDATER(QX2,QY2,DUMY1,2)
#endif
         CALL GWCE_New_pc(IT,TimeLoc,TimeH)
      ENDIF  
      IF(CGWCE_Orig) THEN
         CALL GWCE_Original(IT,TimeLoc,TimeH)
      ENDIF
      IF(CGWCE_New) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
      ENDIF

CWET...
CWET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
CWET...
CWET...NOTE:NNODECODE is a working variable that can change more than once
CWET...               during a time step 
CWET...     NNODECODE = 0 for a dry node
CWET...     NNODECODE = 1 for a wet node
CWET...     NODECODE  - is a more static version of NNODECODE that is reconciled
CWET...                 once and for all at the end time step
CWET...
CWET...
CWET...        (   DRYING CRITERIA   )
CWET...
CWET...A node should be dry under two conditions.
CWET...D1.) If the total water depth falls below H0. 
CWET .......Note: if the total water depth falls below 0.8*H0, the surface elevation
CWET........is lifted up so that the total water depth = 0.8*H0.
CWET......
CWET...D2.) If the node is connected to only nonfunctioning (dry) elements.  In
CWET........this case the node is dried due to becoming landlocked.
CWET........Note: this criteria is applied after all other wetting and drying criteria
CWET...
CWET...An element should be dry under the following conditions.
CWET...DE3.) This is an elemental check section designed to avoid artificial wetting of
CWET.........of control sections
CWET.........All elements where downhill flow originates from a barely wet node
CWET.........(defined as 1.2*H0) into wet nodes are forced inactive; the only exception
CWET......... is receiving overtopped barrier nodes
CWET...
CWET...        (   WETTING CRITERIA   )
CWET...
CWET...A node should be wet under two conditions.
CWET...W1.) If 2 nodes on an element are wet and one is dry, wet the dry node 
CWET........if the water level at one of the wet nodes is greater than the 
CWET........water level at the dry node and the steady state velocity that 
CWET........would result from a balance between the water level gradient and 
CWET........bottom friction would yield a velocity > VELMIN.  
CWET........Note that the criteria outlined in DE3 must also be satified before 
CWET........the node is allowed to wet           
CWET...
CWET...W2.) If an element has a node lying on a receiving internal barrier boundary or 
CWET......specified discharge boundary that is actively discharging flow into the
CWET......domain at that node, all nodes in this element must stay wet. 
CWET...
CWET...
CWET...        (  VELOCITY BOUNDARY CONDITION  )
CWET...
CWET...Either a natural or essential boundary condition can be used as a velocity 
CWET...boundary condition in the momentum equation solution along a wet/dry boudary
CWET...To use a natural boundary condition, do nothing along the wet/dry interface.
CWET...To use an essential, no velocity boundary condition, identify the nodes along
CWET...the wet/dry interface and zero out the velocity at the nodes.  Interface nodes
CWET...can easily be identified by comparing the number of active elements a node is
CWET...connected to (MJU) to the total number of elements a node is connected to (NODELE).
CWET...If MJU < NODELE for any node, it must lie along the wet/dry interface.  See
CWET...further comments at the end of the momentum equation solution section.
CWET...

CWET...
CWET...WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP    
CWET...
      IF(NOLIFA.EQ.2) THEN      !  This goes on until end of part 6
         DO I=1,NP
            NIBCNT(I) = 0        
         ENDDO
         DO I=1,NE
            NOFFOLD(I)=NOFF(I)  
            NOFF(I)=1
         ENDDO

CWET...
CWET...WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
CWET....Drying Criteria D1: this depends on NODECODE and updates NODECODE
CWET...

         DO I=1,NP
            IF(NODECODE(I).EQ.1) THEN
               HTOT=DP(I)+ETA2(I)
               IF(HTOT.LE.H0) THEN
                  IF(HTOT.LT.HABSMIN) ETA2(I)=HABSMIN-DP(I)
                     NNODECODE(I)=0
                     NODECODE(I)=0
                     NCCHANGE=NCCHANGE+1 !NCCHANGE=0 set near beginning of GWCE
C                  ENDIF
               ENDIF         
            ENDIF
         ENDDO
CWET...
CWET...END WET/DRY SECTION - PART 1
CWET...

cjjwC     Use Message-Passing to update nodecode and nnodecode at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(NODECODE,NNODECODE,2)
cjjw#endif

CWET...
CWET...WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
CWET...
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            
CWET...
CWET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
CWET...
            NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
            IF(NCTOT.EQ.2) THEN
               ETAN1=ETA2(NM1)
               ETAN2=ETA2(NM2)
               ETAN3=ETA2(NM3)
               HTOTN1=DP(NM1)+ETA2(NM1)
               HTOTN2=DP(NM2)+ETA2(NM2)
               HTOTN3=DP(NM3)+ETA2(NM3)
               IF((NODECODE(NM1).EQ.1).AND.(NODECODE(NM2).EQ.1)) THEN
                  IF((HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF)) THEN
                     NM123=NM1
                     IF(ETA2(NM2).GT.ETA2(NM1)) NM123=NM2
                     DELDIST=SQRT((y(NM3)-y(NM123))**2.D0
     &                    +(X(NM3)-X(NM123))**2.D0)         
                     DELETA=ETA2(NM123)-ETA2(NM3)
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged from Casey 071219: Added the following logic for 3D friction.
c. RJW modified the following for 3D friction
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))

                        IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                        VEL=G*(DELETA/DELDIST)/TKWET

                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                       IF(LoadZ0B_var) THEN
                           Z0B1 = Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1 = Z0B
                       ENDIF
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                    ENDIF

                     IF(VEL.GT.VELMIN) THEN
C    ....         third node met criteria and is also wet
                        NNODECODE(NM3)=1
c. RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
c                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
                           TKXX(NM123) = TK(NM123)
                           TKYY(NM123) = TK(NM123)
                           TKXY(NM123) = 0.D0
#endif
#endif
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM3)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              ENDIF
                              BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM3)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ELSEIF((NODECODE(NM2).EQ.1).AND.(NODECODE(NM3).EQ.1)) 
     &          THEN
                  IF((HTOTN2.GE.HOFF).AND.(HTOTN3.GE.HOFF)) THEN
                     NM123=NM2
                     IF(ETA2(NM3).GT.ETA2(NM2)) NM123=NM3
                     DELDIST=SQRT((Y(NM1)-Y(NM123))**2.D0
     &                    +(X(NM1)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM1)
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                        IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                        VEL=G*(DELETA/DELDIST)/TKWET

                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                        IF(LoadZ0B_var) THEN
                           Z0B1=Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1=Z0B
                        ENDIF
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     ENDIF

                     IF(VEL.GT.VELMIN) THEN
                        NNODECODE(NM1)=1
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
                           TKXX(NM123) = TK(NM123)
                           TKYY(NM123) = TK(NM123)
                           TKXY(NM123) = 0.D0
#endif
#endif
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM1)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              ENDIF
                              BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM1)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                           ENDIF
                        ENDIF

                     ENDIF
                  ENDIF 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ELSEIF((NODECODE(NM3).EQ.1).AND.(NODECODE(NM1).EQ.1)) 
     &          THEN
                  IF((HTOTN3.GE.HOFF).AND.(HTOTN1.GE.HOFF)) THEN
                     NM123=NM3
                     IF(ETA2(NM1).GT.ETA2(NM3)) NM123=NM1
                     DELDIST=SQRT((Y(NM2)-Y(NM123))**2.D0
     &                     +(X(NM2)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM2)
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                       IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                       VEL=G*(DELETA/DELDIST)/TKWET
                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                        IF(LoadZ0B_var) THEN
                           Z0B1 = Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1 = Z0B
                        ENDIF
                      VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     ENDIF

                     IF(VEL.GT.VELMIN) THEN
                        NNODECODE(NM2)=1
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
                           TKXX(NM123) = TK(NM123)
                           TKYY(NM123) = TK(NM123)
                           TKXY(NM123) = 0.D0
#endif
#endif
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM2)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              ENDIF
                              BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM2)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM2)=KSLIP*AIMAG(Q(NM123,1))
                           ENDIF
                        ENDIF

                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
CWET...
CWET...Nodal Wetting Criteria W2a
CWET...
            NBNCTOT=NIBNODECODE(NM1)+NIBNODECODE(NM2)+NIBNODECODE(NM3)
            NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
            NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
            NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

         ENDDO

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
         CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif

CWET...
CWET...Nodal Wetting Criteria W2b
CWET...Check for adjacent nodes and force nodes wet when attached
CWET...to receiving barrier nodes
CWET...
         DO I=1,NP
            IF((NIBCNT(I).GT.0).AND.(NNODECODE(I).EQ.0)) THEN
               NNODECODE(I)=1  
            ENDIF
         ENDDO

cjjwC     Use Message-Passing to update nnodecode at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(NNODECODE,IDUMY,1)
cjjw#endif

CWET...
CWET...END WET/DRY SECTION - PART 2
CWET...

CWET...
CWET...START WET/DRY SECTION  - PART 3
CWET...Elemental drying criteria DE1
CWET...This is an elemental check section designed to avoid artificial wetting of
CWET....of control sections
CWET...All elements where downhill flow originates from a barely wet node
CWET....into wet nodes are forced inactive; the only exception is receiving
CWET....overtopped barrier nodes
CWET...
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NBNCTOT=NIBCNT(NM1)*NIBCNT(NM2)*NIBCNT(NM3)
            IF(NBNCTOT.EQ.0) THEN   !No barrier/pipe receiving nodes in this elem
               ETAN1=ETA2(NM1)
               ETAN2=ETA2(NM2)
               ETAN3=ETA2(NM3)
               HTOTN1=DP(NM1)+ETA2(NM1)
               HTOTN2=DP(NM2)+ETA2(NM2)
               HTOTN3=DP(NM3)+ETA2(NM3)
#ifdef SB_WETDRY
C...Find the heighest point on the bed in the element.  sb v46.28.sb05.06 11/01/2006
               IF(DP(NM1).LE.DP(NM2).AND.DP(NM1).LE.DP(NM3)) THEN
                  DPMIN = DP(NM1)
               ELSE IF(DP(NM2).LE.DP(NM3).AND.DP(NM2).LE.DP(NM1)) THEN
                  DPMIN = DP(NM2)
               ELSE IF(DP(NM3).LE.DP(NM1).AND.DP(NM3).LE.DP(NM2)) THEN
                  DPMIN = DP(NM3)
               ENDIF
#endif
#ifndef SB_WETDRY
C...ABC pattern
               IF((ETAN1.GE.ETAN2).AND.(ETAN2.GT.ETAN3)) THEN
                 IF((HTOTN1.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
               ENDIF
               IF((ETAN2.GE.ETAN3).AND.(ETAN3.GT.ETAN1)) THEN
                 IF((HTOTN2.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
               ENDIF
               IF((ETAN3.GE.ETAN1).AND.(ETAN1.GT.ETAN2)) THEN
                 IF((HTOTN3.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
               ENDIF
C...ACB pattern
               IF((ETAN1.GE.ETAN3).AND.(ETAN3.GT.ETAN2)) THEN
                 IF((HTOTN1.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
               ENDIF
               IF((ETAN2.GE.ETAN1).AND.(ETAN1.GT.ETAN3)) THEN
                 IF((HTOTN2.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
               ENDIF
               IF((ETAN3.GE.ETAN2).AND.(ETAN2.GT.ETAN1)) THEN
                 IF((HTOTN3.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
               ENDIF
#else
C...ABC pattern
               IF((ETAN1.GE.ETAN2).AND.(ETAN2.GT.ETAN3)) THEN
                 IF((HTOTN1.LT.HOFF)) NOFF(I)=0
                 IF((ETAN1-ETAN2).LT.(ETAN2-ETAN3)) THEN
                   IF(HTOTN2.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
               IF((ETAN2.GE.ETAN3).AND.(ETAN3.GT.ETAN1)) THEN
                 IF((HTOTN2.LT.HOFF)) NOFF(I)=0
                 IF((ETAN2-ETAN3).LT.(ETAN3-ETAN1)) THEN
                   IF(HTOTN3.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
               IF((ETAN3.GE.ETAN1).AND.(ETAN1.GT.ETAN2)) THEN
                 IF((HTOTN3.LT.HOFF)) NOFF(I)=0
                 IF((ETAN3-ETAN1).LT.(ETAN1-ETAN2)) THEN
                   IF(HTOTN1.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
C...ACB pattern
               IF((ETAN1.GE.ETAN3).AND.(ETAN3.GT.ETAN2)) THEN
                 IF((HTOTN1.LT.HOFF)) NOFF(I)=0
                 IF((ETAN1-ETAN3).LT.(ETAN3-ETAN2)) THEN
                   IF(HTOTN3.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
               IF((ETAN2.GE.ETAN1).AND.(ETAN1.GT.ETAN3)) THEN
                 IF((HTOTN2.LT.HOFF)) NOFF(I)=0
                 IF((ETAN2-ETAN1).LT.(ETAN1-ETAN3)) THEN
                   IF(HTOTN1.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
               IF((ETAN3.GE.ETAN2).AND.(ETAN2.GT.ETAN1)) THEN
                 IF((HTOTN3.LT.HOFF)) NOFF(I)=0
                 IF((ETAN3-ETAN2).LT.(ETAN2-ETAN1)) THEN
                   IF(HTOTN2.LT.HOFF) NOFF(I)=0
                 ENDIF
               ENDIF
#endif    
    
#ifdef SB_WETDRY
C...An element is set to be dry if it is determined to be a flooding type
C...wetting element.  An element is a flooding type wetting element if
C...the bed elevation at the node with the biggest water column height  
C...is lower than the heighest point on the bed in the element.  
C...sb v46.28.sb05.06 11/01/2006
C...This is applied only when NOFF flag of the element at the previous time step
C...is 0, which means that this logic works to prevent an element from 
C...re-wetting.
C...sb v46.52.03
               IF(NOFFOLD(I).EQ.0) THEN
                 IF(HTOTN1.GE.HTOTN2.AND.HTOTN1.GE.HTOTN3) THEN
                    IF(ETAN1.LT.(-DPMIN+H0)) NOFF(I) = 0
                 ENDIF
                 IF(HTOTN2.GE.HTOTN3.AND.HTOTN2.GE.HTOTN1) THEN
                    IF(ETAN2.LT.(-DPMIN+H0)) NOFF(I) = 0
                 ENDIF
                 IF(HTOTN3.GE.HTOTN1.AND.HTOTN3.GE.HTOTN2) THEN
                    IF(ETAN3.LT.(-DPMIN+H0)) NOFF(I) = 0
                 ENDIF
              ENDIF
#endif

            ENDIF
         ENDDO

#ifdef SB_WETDRY
CWET......added by sb on 11/02/2006
CWET...
CWET...This section is added after we realize that it's not possible 
CWET...to compute a correct flow going through two elements if
CWET...the elements are connected just by one node. i.e., elements need
CWET...to share an edge to let the flow go through between the elements.
CWET...Therefore, in this section, a node is determined to be dry
CWET...if two elements are connected at one node, not sharing an edge.
CWET...As it seemed this procedure needed NOFF information,
CWET...although this section changes NNODECODE,
CWET...I put this section here, rather than the end of PART 2.
CWET...I tried setting NOFF(I) = 0, but it didn't shut down the flow.
CWET...
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            IF(NOFF(I).EQ.1.AND.
     &         NNODECODE(NM1).EQ.1.AND.
     &         NNODECODE(NM2).EQ.1.AND.
     &         NNODECODE(NM3).EQ.1) THEN
               DO K=1,3
                  NM1=NM(I,K)
                  NM2=NM(I,MOD(K+0,3)+1)
                  NM3=NM(I,MOD(K+1,3)+1)

                  NWETNEI = 0
                  NWETADJ = 0
                  DO J=1,MNEI
                     N=NeiTabEle(NM1,J)
                     IF(N.EQ.0) CYCLE
                     IF(N.EQ.I) CYCLE

                     NMN1=NM(N,1)
                     NMN2=NM(N,2)
                     NMN3=NM(N,3)
                     IF(NOFF(N).EQ.1.AND.
     &                  NNODECODE(NMN1).EQ.1.AND.
     &                  NNODECODE(NMN2).EQ.1.AND.
     &                  NNODECODE(NMN3).EQ.1) THEN
                        NWETNEI = NWETNEI + 1

                        IF((NMN1.EQ.NM2.OR.NMN1.EQ.NM3).OR.
     &                     (NMN2.EQ.NM2.OR.NMN2.EQ.NM3).OR.
     &                     (NMN3.EQ.NM2.OR.NMN3.EQ.NM3)) THEN
                           NWETADJ = NWETADJ + 1
                        ENDIF
                     ENDIF
                  ENDDO

                  IF(NWETNEI.GT.0.AND.NWETADJ.EQ.0.AND.
     &               NIBCNT(NM1).EQ.0) THEN
                     NNODECODE(NM1) = 0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
#endif

CWET...
CWET...END WET/DRY SECTION  - PART 3
CWET...

CWET...  
CWET...START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2 
CWET...Update number of active elements (MJU) and the total area (TotalArea) connected  
CWET...to a node. If these are zero, the node is landlocked and should be dried.
CWET...These depend on NNODECODE which varies during the time step
CWET...
         DO I=1,NP
            MJU(I)=0
            TotalArea(I)=0.d0
         ENDDO
         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NNODECODE(NM1)
            NC2=NNODECODE(NM2)
            NC3=NNODECODE(NM3)

            NCEle=NC1*NC2*NC3*NOFF(IE)
            AreaEle=NCEle*Areas(IE)/2.d0
            MJU(NM1)=MJU(NM1)+NCEle
            MJU(NM2)=MJU(NM2)+NCEle
            MJU(NM3)=MJU(NM3)+NCEle
            TotalArea(NM1)=TotalArea(NM1)+AreaEle
            TotalArea(NM2)=TotalArea(NM2)+AreaEle
            TotalArea(NM3)=TotalArea(NM3)+AreaEle
         ENDDO

cjjwnote - looks like this is used later in momentum equations
cjjwnote - this has implications on making this into a subroutine   

         DO I=1,NP
            IF((NNODECODE(I).EQ.1).AND.(MJU(I).EQ.0)) THEN
               NNODECODE(I)=0
            ENDIF
            IF(MJU(I).EQ.0) MJU(I)=1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
         ENDDO

C     WET...
C     WET...END WET/DRY SECTION - PART 4
C     WET...

cjjwnote - may have to pass TotalArea and mju as well               
 
C     Use Message-Passing to update nnodecode at ghost nodes
#ifdef CMPI
         CALL UPDATEI(NNODECODE,IDUMY,1)
#endif
CWET...  
CWET...WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE 
CWET...Check to see if any wetting occurred & update NODECODE 
CWET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
CWET...       
         DO I=1,NP
            IF(NNODECODE(I).NE.NODECODE(I)) THEN
               NODECODE(I)=NNODECODE(I)
               NCCHANGE=NCCHANGE+1             
            ENDIF
         ENDDO
CWET...
CWET...END WET/DRY SECTION - PART 5
CWET...

CWET...  
CWET...WET/DRY SECTION - PART 6 
CWET...Check to see if any NOFF changed requiring the matrix to be reset
CWET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
CWET...      
         DO I=1,NE
            IF(NOFF(I).NE.NOFFOLD(I)) NCCHANGE=NCCHANGE+1
         ENDDO
CWET...
CWET... jgf45.06 If there has been any wetting or drying in any
CWET... of the subdomains, the NCCHANGE flag will be activated on all 
CWET... of the subdomains, to prevent them from getting out of sync
CWET... with their MPI calls as some reset the GWCE and others do not.
CWET...
#ifdef CMPI
         !jgf48.4619 implementing Seizo's changes for Lumped, fully 
         ! explicit operation. In that case, the GWCE LHS matrix is 
         ! recalculated on each individual subdomain that has wetted
         ! or dried, without recourse to MPI, eliminating the need
         ! for the call to the subroutine WetDrySum.
         IF ( ILump.eq.0 ) THEN
            call WetDrySum(NCCHANGE)
         ELSE
            NCCHANGE=NCCHANGE ! jgf48.4619 do nothing
         ENDIF
#endif
CWET...
CWET...END WET/DRY SECTION - PART 6
CWET...
      ENDIF                     !  This is started in Part 1 of CWET
c.....
c. RJW merged 08/26/2008 Casey 071219: Added these output lines for wet/dry visualization.
!      IF(NOUTGE.NE.0)THEN
!         IF((IT.GT.NTCYSGE).AND.(IT.LE.NTCYFGE))THEN
!            IF(NSCOUGE.EQ.(NSPOOLGE-1))THEN
!               CLOSE(123)
!#ifdef CMPI
!               WRITE(TEMPDIRNAME,'(A2,I4.4)') 'PE',MYPROC
!               OPEN(123,FILE=TEMPDIRNAME//'/'//'fort.123',
!     &            ACCESS='SEQUENTIAL',POSITION='APPEND')
!#else
!               OPEN(123,FILE='fort.123',
!     &            ACCESS='SEQUENTIAL',POSITION='APPEND')
!#endif
!               DO I=1,NP
!                  WRITE(123,*) NODECODE(I)
!               ENDDO
!               CLOSE(124)
!#ifdef CMPI
!               OPEN(124,FILE=TEMPDIRNAME//'/'//'fort.124',
!     &            ACCESS='SEQUENTIAL',POSITION='APPEND')
!#else
!               OPEN(124,FILE='fort.124',
!     &            ACCESS='SEQUENTIAL',POSITION='APPEND')
!#endif
!               DO IE=1,NE
!                  WRITE(124,*) NOFF(IE)
!               ENDDO
!            ENDIF
!         ENDIF
!      ENDIF

c.end rjw merge vis....
C...  
C...  2DDI Momentum Equation Solution
C...  
      IF (C2DDI) THEN    
         IF (CME_Orig) THEN
            CALL Mom_Eqs_Original()
         ENDIF
         IF (CME_New_NC) THEN
            CALL Mom_Eqs_New_NC()
         ENDIF
         IF ((CME_New_C1).OR.(CME_New_C2)) THEN
            CALL Mom_Eqs_New_Conserv()
         ENDIF
         IF (CPRECOR) THEN
            CALL Mom_Eqs_Non_Conserv_pc()
         ENDIF

C...  If running in parallel, update velocities & fluxes on all processors

#ifdef CMPI
         CALL UPDATER(UU2,VV2,DUMY1,2)
         CALL UPDATER(QX2,QY2,DUMY1,2)
#endif
     
      ENDIF 

C...  
C     C2DDI....END OF 2DDI MOMENTUM EQUATION SOLUTION
C...  



C...  
C...  3DVS Momentum Equation Solution
C...
      IF (C3DVS) THEN

C... Load the vector MOM_LV_X(I) with barotropic pressure terms 
C...     including atmospheric pressure, water level and tidal potential
C...     averaged between time levels s and s+1, (time levels 1 and 2).  Note:
C...     MOM_LV_X gets renamed as BTP in global_3dvs.f

        DO I=1,NP
           MOM_LV_X(I)=ETA1(I)+ETA2(I)
           IF(NWS.NE.0) MOM_LV_X(I)=MOM_LV_X(I)+PR1(I)+PR2(I) !atmospheric pressure
           IF (CTIP) MOM_LV_X(I)=MOM_LV_X(I)-TIP1(I)-TIP2(I) !tidal potential
           MOM_LV_X(I)=G*MOM_LV_X(I)/2.d0
        ENDDO

C...  Solve for velocity at the new time level (K+1)

        CALL VSSOL(IT,TimeLoc)

      ENDIF
C...  
C...  End of 3DVS Momentum Equation Solution
C...  


C...  
C...  IF 2D TRANSPORT IS INCLUDED SOLVE FOR THE CONCENTRATION
C...  
      IF(C2D_PTrans) THEN

        CALL SCALAR_TRANS_2D (IT, TimeLoc)

      ENDIF
C...
C...  End of 2D Scalar Transport Solution
C...
C     jgf48.4627 Jump to here if METONLY is .TRUE., i.e., if only 
C     meteorological output is requested.
9999  CONTINUE  
C...
C...  Collect maximum values of variables   v46.50 sb 11/11/2006
C...
      do i=1,np
         if(nodecode(i).eq.1.and.etamax(i).lt.eta2(i)) then
            etamax(i) = eta2(i)
         endif
         if((umax(i)*umax(i)).lt.(uu2(i)*uu2(i)+vv2(i)*vv2(i))) then
            umax(i) = sqrt(uu2(i)*uu2(i)+vv2(i)*vv2(i))
         endif
         if(prmin(i).gt.pr2(i)) then
            prmin(i) = pr2(i)
         endif
         if((wvnoutmax(i)*wvnoutmax(i)).lt.
     &        (wvnxout(i)*wvnxout(i)+wvnyout(i)*wvnyout(i))) then
            wvnoutmax(i) 
     &          = sqrt(wvnxout(i)*wvnxout(i)+wvnyout(i)*wvnyout(i))
         endif
         if((rsnmax(i)*rsnmax(i)).lt.
     &        (rsnx2(i)*rsnx2(i)+rsny2(i)*rsny2(i))) then
            rsnmax(i) = sqrt(rsnx2(i)*rsnx2(i)+rsny2(i)*rsny2(i))
         endif
      enddo
C      
C...
C...  WRITE OUTPUT
C...
      CALL WRITE_OUTPUT_2D (IT,TimeLoc) ! =>see write_output.F
C
C 
C     jgf49.44: If harmonic analysis was requested and the current time
C     is within the harmonic analysis period, update the left hand side 
C     of the harmonic analysis matrix. Also update the load vectors for
C     each type of analysis. If timeseries reconstruction was specified,
C     also update the timeseries.
      CALL updateHarmonicAnalysis(IT, TIMEH)
C...
C...  WRITE OUT HOT START INFORMATION IF NHSTAR=1 AND AT CORRECT TIME
C.... STEP
C...  NOTE: THE HOT START FILES USE A RECORD LENGTH OF 8 ON BOTH 32 BIT
C.... WORKSTATIONS AND THE 64 BIT CRAY.  THIS IS BECAUSE THE HARMONIC
C.... ANALYSIS IS DONE IN DOUBLE PRECISION (64 BITS) ON WORKSTATIONS.
C...
     
      ITEST=(IT/NHSINC)*NHSINC  
      
!     IF(myproc.eq.0) PRINT *, " ITEST **********", ITEST
      
      if ((NHSTAR > 0 .and. ITEST == IT) .or. (-IHOT == IT)) then
        if( (MNWPROH > 0) ) then   !Writer for HSfile
          if( (.not.C3D) ) then                     !Writer for HSfile
#ifdef CMPI
          CALL writeHotstart_through_hswriter(TimeLoc,IT)  !st3 hsfile
#endif
          else
            write(6,*) 'HS writer does not support C3D'
          endif
        else
!     IF(myproc.eq.0) PRINT *, " writeHotstart **********"
         CALL writeHotstart(TimeLoc, IT)
        endif

#ifdef CSWAN
Casey 100205: Enable writing of SWAN hot-start file.  We need to wait
C             and do this after the next SWAN time step, so that
C             everything is up-to-date.
         WriteSwanHotStart = .TRUE.
#endif
      ENDIF


C...
c...  find and print to unit 6, the maximum elevation, the maximum
c...  velocity and the node numbers at which they occur on myproc=0 if
c...  elmax exceeds threshold, print information on all processors where
c...  this occurs
C...
C
C     jgf46.00 Added option to output data to the screen every NSCREEN
C     time steps, rather than on every time step, as long as there are
C     no high elevations. In the case of high elevations, the warning 
C     messages are sent to the screen each time they are generated.
C
      IF(NSCREEN.NE.0) THEN
         ELMAX=0.0d0
         VELMAX=0.0d0
         KEMAX = 0
         KVMAX = 0
         DO I=1,NP
            IF((NODECODE(I).EQ.1).AND.(ABS(ETA2(I)).GT.ELMAX))THEN
               ELMAX=ABS(ETA2(I))
               KEMAX=I
            ENDIF
            VELABS=UU2(I)*UU2(I)+VV2(I)*VV2(I)
            IF (VELABS.GT.VELMAX) THEN
               VELMAX=VELABS
               KVMAX=I
            ENDIF
         END DO
         VELMAX=VELMAX**0.5d0
         ITEST=(IT/NSCREEN)*NSCREEN 

C     jgf46.10 Added the ability to for the user to control the warning
C     and error elevations. Also added the ability for the user to write
C     a fort.69 (global elevation debug) file. 

C     jgf46.12 Removed the dependence on KEMAX for producing output to
C     the screen.
#ifdef CMPI
         IF(MYPROC.EQ.0.AND.ELMAX.LT.WarnElev.AND.ITEST.EQ.IT) THEN
            IF (KEMAX.GT.0) THEN
               WRITE(ScreenUnit,1991) 
     &            IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,MYPROC
            ELSE
               WRITE(ScreenUnit,1991)
     &            IT,NUMITR,TimeLoc,0.,KEMAX,VELMAX,KVMAX,MYPROC
            ENDIF
 1991       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E14.8,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &           2X,'ON MYPROC = ',I4)
         ENDIF
         WarnElevExceeded = 0
         IF(ELMAX.GT.WarnElev) THEN
            WRITE(ScreenUnit,1993) 
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,
     &           KVMAX,MYPROC
            WRITE(16,1993) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,
     &           KVMAX,MYPROC
 1993       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E14.8,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &           2X,'ON MYPROC = ',I4,
     &           3X,'** WARNING: Elevation.gt.WarnElev **')
            IF (WarnElevDump) WarnElevExceeded=1
         ENDIF
#ifdef DEBUG_WARN_ELEV
         call WarnElevSum(WarnElevExceeded)
         IF (WarnElevExceeded.ne.0) THEN
            CALL WriteWarnElev(TimeLoc, IT)
         ENDIF
#endif
         ErrorElevExceeded = 0                ! Clint's Zombie Slyaer
         IF(ELMAX.GT.ErrorElev) THEN
            WRITE(ScreenUnit,1995) 
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,
     &           KVMAX,MYPROC
            WRITE(16,1995) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,
     &           KVMAX,MYPROC
 1995       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E14.8,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &           2X,'ON MYPROC = ',I4,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev,',
     &           ' ADCIRC stopping. **')
         ErrorElevExceeded = 1                ! Clint's Zombie Slayer
         ENDIF
         call WarnElevSum(ErrorElevExceeded)  ! Clint's Zombie Slayer 2010.08.07
         IF( ErrorElevExceeded /= 0 ) THEN    !  Finalize MPI Environment,
           CALL MSG_FINI()                    !  if there are Error Elvation nodes.
           STOP                               !
         ENDIF                                !st3
#else 
         IF(ELMAX.LT.WarnElev.AND.ITEST.EQ.IT) THEN
            IF (KEMAX.GT.0) THEN
               WRITE(ScreenUnit,1992) 
     &              IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            ELSE
               WRITE(ScreenUnit,1992) 
     &              IT,NUMITR,TimeLoc,0.,KEMAX,VELMAX,KVMAX
            ENDIF
 1992       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7)
         ENDIF
         IF(ELMAX.GT.WarnElev) THEN
            WRITE(ScreenUnit,1994) 
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1994) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,
     &           2X,'** WARNING: Elevation.gt.WarnElev **')
#ifdef DEBUG_WARN_ELEV
            IF (WarnElevDump) CALL WriteWarnElev(TimeLoc, IT)
#endif
         ENDIF
         IF(ELMAX.GT.ErrorElev) THEN
            WRITE(ScreenUnit,1996) 
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1996) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1996       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev, '
     &           'ADCIRC stopping. **')
            STOP
         ENDIF
#endif
      ENDIF

C...  
C...  ****************** TIME STEPPING LOOP ENDS HERE ********************
C...  

c. RJW merged 08/26/2008 Casey 071219: Added the folowing call to the mass balance subroutine.
C      IF(C3DVS)THEN
C         CALL MASSBAL3D(IT)
C      ENDIF

      RETURN
      END


C*******************************************************************************
C                                                                              *
C    Subroutine to compute the elevation using the original GWCE formluation   *
C                                                                              *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE GWCE_Original(IT,TimeLoc,TimeH)
C     
      USE GLOBAL
      USE SIZES
      USE WIND
      USE ITPACKV
      USE NodalAttributes, ONLY : 
     &     Tau0Var, LoadGeoidOffset, GeoidOffset, EVM
#ifdef CMPI
      USE MESSENGER
#endif
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      USE Couple2Swan, ONLY: TKXX,
     &                       TKXY,
     &                       TKYY
#endif
#endif

      IMPLICIT NONE 

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI

      REAL(SZ) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3
      REAL(SZ) BSXAvg3, BSYAvg3
      REAL(SZ) BXN1, BXN2, BXN3, BYN1, BYN2, BYN3
      REAL(SZ) CorifAvg
      REAL(SZ) CoefO2
      REAL(SZ) DXN1, DXN2, DXN3, DYN1, DYN2, DYN3
      REAL(SZ) DEta1SqDX, DEta1SqDY
      REAL(SZ) DPAvg
      REAL(SZ) DPrMTiPDX, DPrMTiPDY
      REAL(SZ) DU1DX, DU1DY, DV1DX, DV1DY
      REAL(SZ) DU1QX1DX,  DV1QX1DY,  DU1QY1DX,  DV1QY1DY
      REAL(SZ) DUU1N1, DUU1N2, DUU1N3
      REAL(SZ) DUV1N1, DUV1N2, DUV1N3
      REAL(SZ) DVV1N1, DVV1N2, DVV1N3
      REAL(SZ) DXXYY11, DXXYY12, DXXYY13, DXXYY21
      REAL(SZ) DXXYY22, DXXYY23, DXXYY31
      REAL(SZ) DXXYY32, DXXYY33, DXYH11, DXYH12
      REAL(SZ) DXYH13, DXYH21, DXYH22
      REAL(SZ) DXYH23, DXYH31, DXYH32, DXYH33
      REAL(SZ) E0N1, E0N2, E0N3
      REAL(SZ) E1N1, E1N2, E1N3, E1N1SQ, E1N2SQ, E1N3SQ 
      REAL(SZ) ESN1, ESN2, ESN3
      REAL(SZ) EVMAvgODT
      REAL(SZ) GFac, GHAvg, GH1Fac
      REAL(SZ) H1, H1N1, H1N2, H1N3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) HPPO2
      REAL(SZ) PR1N1, PR1N2, PR1N3
      REAL(SZ) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3
      REAL(SZ) QX1Avg, QY1Avg, QX1Avg3, QY1Avg3
      REAL(SZ) SFacAvg
      REAL(SZ) SpaVarTau0
      REAL(SZ) T0N1,T0N2, T0N3, T0XN1, T0XN2, T0XN3, T0YN1, T0YN2, T0YN3
      REAL(SZ) Tau0Avg, T0XAvg3, T0YAvg3
      REAL(SZ) TT0L, TT0R
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) TIPN1, TIPN2, TIPN3
      REAL(SZ) U1N1,U1N2,U1N3, U1Avg
      REAL(SZ) V1N1,V1N2,V1N3, V1Avg
      REAL(SZ) VIDBCPDXAvg3, VIDBCPDYAvg3
      REAL(SZ) WSXN1, WSXN2, WSXN3, WSYN1, WSYN2, WSYN3
      REAL(SZ) WSXAvg3, WSYAvg3

      REAL(8) AREAIE, AREAIE2, AO6, AO12
      REAL(8) FDDD, FDDOD    
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) DDX1,DDX2,DDX3,DDY1,DDY2,DDY3
      REAL(8) DXX11,DXX12,DXX13,DXX21,DXX22,DXX23,DXX31,DXX32,DXX33
      REAL(8) DYY11,DYY12,DYY13,DYY21,DYY22,DYY23,DYY31,DYY32,DYY33
      REAL(8) DXY11,DXY12,DXY13,DXY21,DXY22,DXY23,DXY31,DXY32,DXY33
      REAL(8) TimeLoc, TimeH

C...  RECOMPUTE THE GWCE SYSTEM MATRIX AT THE FIRST TIME STEP OR IF ANY
C...  WETTING OR DRYING OCCURRED IN THE PREVIOUS TIME STEP.
C...  
      IF(NCCHANGE.GT.0) THEN !if any subdomain grid has changed
         NCCHANGE=0
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3806)
c     WRITE(16,3806)
 3806    FORMAT(/,1X,'RE-SETTING GWCE SYSTEM MATRIX',/)

C.....Set up the LHS matrix (for the iterative matrix solver)
         COEF(:,:)=0.0d0

         DO IE=1,NE
            NMI1=NM(IE,1)
            NMI2=NM(IE,2)
            NMI3=NM(IE,3)
            NMJ1=NMI1
            NMJ2=NMI2
            NMJ3=NMI3
            NC1=NODECODE(NMI1)
            NC2=NODECODE(NMI2)
            NC3=NODECODE(NMI3)
            NCELE=NC1*NC2*NC3*NOFF(IE)

            SFacAvg=(SFAC(NMI1)+SFAC(NMI2)+SFAC(NMI3))/3.d0
            Tau0Avg=(TAU0VAR(NMI1)+TAU0VAR(NMI2)+TAU0VAR(NMI3))/3.d0
            TT0L=((1.0d0+0.5d0*DT*Tau0Avg)/DT)/DT

            FDX1 = (Y(NMI2)-Y(NMI3))*SFacAvg         !b1
            FDX2 = (Y(NMI3)-Y(NMI1))*SFacAvg         !b2
            FDX3 = (Y(NMI1)-Y(NMI2))*SFacAvg         !b3
            FDY1 = X(NMI3)-X(NMI2)                   !a1
            FDY2 = X(NMI1)-X(NMI3)                   !a2
            FDY3 = X(NMI2)-X(NMI1)                   !a3

            DPAvg=(DP(NMI1)+DP(NMI2)+DP(NMI3))/3.d0

            AREAIE=AREAS(IE)/2.0d0
            AO6=AREAIE/6.d0
            AO12=AREAIE/12.d0
            FDDD=(1+ILUMP)*AO6             !<phi*phj> diagonal terms
            FDDOD=(1-ILUMP)*AO12           !<phi*phj> off diagonal terms

            DO JN=2,NEIMAX
               IF(NeiTab(NMI1,JN).EQ.NMJ2) J12=JN
               IF(NeiTab(NMI1,JN).EQ.NMJ3) J13=JN
               IF(NeiTab(NMI2,JN).EQ.NMJ1) J21=JN
               IF(NeiTab(NMI2,JN).EQ.NMJ3) J23=JN
               IF(NeiTab(NMI3,JN).EQ.NMJ1) J31=JN
               IF(NeiTab(NMI3,JN).EQ.NMJ2) J32=JN
            END DO

            GFac=GA00*DPAvg/(4.d0*AREAIE)
            COEF(NMI1,1)  =COEF(NMI1,1)
     &               +2.d0*(TT0L*FDDD +GFac*(FDY1*FDY1+FDX1*FDX1))*NCELE
            COEF(NMI1,J12)=COEF(NMI1,J12)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY1*FDY2+FDX1*FDX2))*NCELE
            COEF(NMI1,J13)=COEF(NMI1,J13)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY1*FDY3+FDX1*FDX3))*NCELE
            COEF(NMI2,J21)=COEF(NMI2,J21)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY2*FDY1+FDX2*FDX1))*NCELE
            COEF(NMI2,1)  =COEF(NMI2,1)
     &               +2.d0*(TT0L*FDDD +GFac*(FDY2*FDY2+FDX2*FDX2))*NCELE
            COEF(NMI2,J23)=COEF(NMI2,J23)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY2*FDY3+FDX2*FDX3))*NCELE
            COEF(NMI3,J31)=COEF(NMI3,J31)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY3*FDY1+FDX3*FDX1))*NCELE
            COEF(NMI3,J32)=COEF(NMI3,J32)
     &               +2.d0*(TT0L*FDDOD+GFac*(FDY3*FDY2+FDX3*FDX2))*NCELE
            COEF(NMI3,1)  =COEF(NMI3,1)
     &               +2.d0*(TT0L*FDDD +GFac*(FDY3*FDY3+FDX3*FDX3))*NCELE

c...  NOTE: the factor of 2 that multiplies each of the lines above was
c...  only only added for backward compatability with the RHS
c...  formluation in earlier versions of the code.  It should be removed
c...  as soon as possible.  RL 4/12/02

         END DO

C...  Modify the matrix "COEF" by imposing the elevation specified
C...  boundary conditions while maintaining the symmetry of the system
c...  This could be redone every time matrix is recomputed...............rl 2/25/04
c...  Then will no longer need to pass IT
cjjw         IF ((IT.EQ.1).OR.(IT.EQ.ITHS+1)) THEN
        IF (IT.EQ.1) THEN
#ifdef CMPI
            EP = PSDOT(NP,COEF(1,1),COEF(1,1))
            EP = SQRT(RNP_GLOBAL*EP)
#else 
            EP=0.0D0
            DO I=1,NP
               EP=EP+COEF(I,1)*COEF(I,1)
            END DO
            EP=SQRT(EP/NP)
#endif
         ENDIF

c...  for each elevation specified boundary node, zero all off diagonal
c...  terms on the row and set diagonal term to EP

         DO I=1,NETA
            COEF(NBD(I),1)=EP
            DO J=2,NNEIGH(NBD(I))
               COEF(NBD(I),J)=0.0d0
            END DO
         END DO

c...  for each elevation specified boundary node, zero all off diagnoal
c...  terms on the column but save these to be multiplied by the
c...  boundary value and subtracted from the RHS

         DO I=1,NETA
            DO J=2,NNeigh(NBD(I))
               DO IJ=2,NNeigh(NeiTab(NBD(I),J))
                  IF(NBD(I).EQ.NeiTab(NeiTab(NBD(I),J),IJ)) THEN
                     OBCCOEF(I,J-1)=COEF(NeiTab(NBD(I),J),IJ)
                     COEF(NeiTab(NBD(I),J),IJ)=0.0d0
                  ENDIF
               END DO
            END DO
         END DO

C.....Check that all the diagonal elements in "COEF" are > 0.
         DO I=1,NP
            IF(COEF(I,1).EQ.0.d0) COEF(I,1)=EP
            IF(COEF(I,1).LT.0.d0) THEN
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &          WRITE(ScreenUnit,1019) I,COEF(I,1)
               WRITE(16,1019) I,COEF(I,1)
 1019          FORMAT(/,1X,'!!!!!!!!  ERROR  !!!!!!!',
     &           /,1X,'THE DIAGONAL TERM IN THE EQUATION FOR NODE ',I10,
     &             '= ',E15.6,' AND IS <= 0',/)
               STOP
            END IF
         END DO

      ENDIF                     !END OF GWCE MATRIX SETUP
C...  
C...  COMPUTE THE THE LOAD VECTOR GWCE_LV FOR THE GWCE ELEMENT BY ELEMENT
C...  BY FORMING TEMPORARY VECTORS AND THEN ASSEMBLING AT THE END.
C...  THE FOLLOWING ASSEMBLY LOOPS HAVE BEEN UNROLLED TO OPTIMIZE
C...  VECORIZATION
C...  

C...  Initialize variables to zero if forcing is not used
C...  
      if((NWS.NE.0).OR.(NRS.NE.0)) then
      else
         WSXN1=0.d0
         WSXN2=0.d0
         WSXN3=0.d0
         WSYN1=0.d0
         WSYN2=0.d0
         WSYN3=0.d0
         PR1N1=0.d0
         PR1N2=0.d0
         PR1N3=0.d0
      endif

      if (CTIP) then
      else
         TIPN1=0.d0
         TIPN2=0.d0
         TIPN3=0.d0
      endif

      if (CBaroclinic) then
      else
         VIDBCPDXAvg3=0.D0
         VIDBCPDYAvg3=0.D0
      endif


      DO 1037 IE=1,NE
C...  
C...  SET NODAL VALUES FOR EACH ELEMENT
C...  
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         E0N1=ETA1(NM1)
         E0N2=ETA1(NM2)
         E0N3=ETA1(NM3)
         E1N1=ETA2(NM1)
         E1N2=ETA2(NM2)
         E1N3=ETA2(NM3)
         E1N1SQ=E1N1*E1N1
         E1N2SQ=E1N2*E1N2
         E1N3SQ=E1N3*E1N3
         ESN1=ETAS(NM1)
         ESN2=ETAS(NM2)
         ESN3=ETAS(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

C     jjw - m10 - define various tau0 values                   
         T0N1=Tau0Var(NM1) 
         T0N2=Tau0Var(NM2) 
         T0N3=Tau0Var(NM3) 
         Tau0Avg=(T0N1+T0N2+T0N3)/3.d0
         TT0R=((0.5d0*Tau0Avg*DT-1.0d0)/DT)/DT
         T0XN1=T0N1*QX1N1
         T0YN1=T0N1*QY1N1
         T0XN2=T0N2*QX1N2
         T0YN2=T0N2*QY1N2
         T0XN3=T0N3*QX1N3
         T0YN3=T0N3*QY1N3

         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            WSXN1=WSX1(NM1)
            WSXN2=WSX1(NM2)
            WSXN3=WSX1(NM3)
            WSYN1=WSY1(NM1)
            WSYN2=WSY1(NM2)
            WSYN3=WSY1(NM3)
            PR1N1=PR1(NM1)
            PR1N2=PR1(NM2)
            PR1N3=PR1(NM3)
         ENDIF

CTIP..If using tidal potential terms
         if (CTIP) then         !tidal potential
            TIPN1=TIP1(NM1)
            TIPN2=TIP1(NM2)
            TIPN3=TIP1(NM3)
         endif

C2DDI.For the 2DDI version of the code
C2DDI.Compute bottom friction
         if (C2DDI) then        !2D bottom friction
            BSXN1=TK(NM1)*QX1N1
            BSYN1=TK(NM1)*QY1N1
            BSXN2=TK(NM2)*QX1N2
            BSYN2=TK(NM2)*QY1N2
            BSXN3=TK(NM3)*QX1N3
            BSYN3=TK(NM3)*QY1N3
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            BSXN1 = TKXX(NM1)*QX1N1 + TKXY(NM1)*QY1N1
            BSYN1 = TKXY(NM1)*QX1N1 + TKYY(NM1)*QY1N1
            BSXN2 = TKXX(NM2)*QX1N2 + TKXY(NM2)*QY1N2
            BSYN2 = TKXY(NM2)*QX1N2 + TKYY(NM2)*QY1N2
            BSXN3 = TKXX(NM3)*QX1N3 + TKXY(NM3)*QY1N3
            BSYN3 = TKXY(NM3)*QX1N3 + TKYY(NM3)*QY1N3
#endif
#endif
         endif      

C3D...For the 3D version of the code
C3D...Setup bottom friction and velocity dispersion
         if (C3D) then          !3D bottom friction & dispersion
            BSXN1=BSX1(NM1)
            BSXN2=BSX1(NM2)
            BSXN3=BSX1(NM3)
            BSYN1=BSY1(NM1)
            BSYN2=BSY1(NM2)
            BSYN3=BSY1(NM3)
            DVV1N1=DVV1(NM1)
            DVV1N2=DVV1(NM2)
            DVV1N3=DVV1(NM3)
            DUV1N1=DUV1(NM1)
            DUV1N2=DUV1(NM2)
            DUV1N3=DUV1(NM3)
            DUU1N1=DUU1(NM1)
            DUU1N2=DUU1(NM2)
            DUU1N3=DUU1(NM3)
         endif


C...  
C...  COMPUTE ELEMENT AVERAGED QUANTITIES
C... 
         QX1Avg3=QX1N1+QX1N2+QX1N3
         QY1Avg3=QY1N1+QY1N2+QY1N3
         QX1Avg=QX1Avg3/3.d0
         QY1Avg=QY1Avg3/3.d0
         DPAvg=(DP(NM1)+DP(NM2)+DP(NM3))/3.d0
         GHAvg=GO3*(H1N1+H1N2+H1N3)
         U1Avg=(U1N1+U1N2+U1N3)/3.d0
         V1Avg=(V1N1+V1N2+V1N3)/3.d0
         EVMAvgODT=((EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0)/DT
         CorifAvg=(CORIF(NM1)+CORIF(NM2)+CORIF(NM3))/3.d0
         BSXAvg3=BSXn1+BSXn2+BSXn3
         BSYAvg3=BSYn1+BSYn2+BSYn3
         WSXAvg3=WSXN1+WSXN2+WSXN3
         WSYAvg3=WSYN1+WSYN2+WSYN3
         T0XAvg3=T0XN1+T0XN2+T0XN3
         T0YAvg3=T0YN1+T0YN2+T0YN3
         IF(CBaroclinic) THEN
            VIDBCPDXAvg3=H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                       +H1N3*VIDBCPDXOH(NM3)
            VIDBCPDYAvg3=H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                       +H1N3*VIDBCPDYOH(NM3)
         ENDIF

C...
C...  Compute elemental coefficients
C...
         AREAIE2=AREAS(IE)               !2A

         FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
         FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
         FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
         FDY1 = X(NM3)-X(NM2)           !a1
         FDY2 = X(NM1)-X(NM3)           !a2
         FDY3 = X(NM2)-X(NM1)           !a3

         FDX1O2A=FDX1/AREAIE2     !dphi1/dx
         FDY1O2A=FDY1/AREAIE2     !dphi1/dy
         FDX2O2A=FDX2/AREAIE2     !dphi2/dx
         FDY2O2A=FDY2/AREAIE2     !dphi2/dy
         FDX3O2A=FDX3/AREAIE2     !dphi3/dx
         FDY3O2A=FDY3/AREAIE2     !dphi3/dy

C...  
C...  COMPUTE ADDITIONAL ELEMENTAL COEFFICIENTS
C...  
         DDX1=FDX1/3.d0         !<2*(dphi1/dx)*phij> j=1,2,3
         DDY1=FDY1/3.d0         !<2*(dphi1/dy)*phij> j=1,2,3
         DXX11=FDX1O2A*FDX1      !<2*(dphi1/dx)*(dphi1/dx)>
         DYY11=FDY1O2A*FDY1      !<2*(dphi1/dy)*(dphi1/dy)>
         DXY11=FDX1O2A*FDY1      !<2*(dphi1/dx)*(dphi1/dy)>
         DXXYY11=DXX11+DYY11
         DXYH11=DPAvg*DXXYY11   
         DXX12=FDX1O2A*FDX2      !<2*(dphi1/dx)*(dphi2/dx)>
         DYY12=FDY1O2A*FDY2      !<2*(dphi1/dy)*(dphi2/dy)>
         DXY12=FDX1O2A*FDY2      !<2*(dphi1/dx)*(dphi2/dy)>
         DXXYY12=DXX12+DYY12
         DXYH12=DPAvg*DXXYY12
         DXX13=FDX1O2A*FDX3      !<2*(dphi1/dx)*(dphi3/dx)>
         DYY13=FDY1O2A*FDY3      !<2*(dphi1/dy)*(dphi3/dy)>
         DXY13=FDX1O2A*FDY3      !<2*(dphi1/dx)*(dphi3/dy)>
         DXXYY13=DXX13+DYY13
         DXYH13=DPAvg*DXXYY13

         DDX2=FDX2/3.d0         !<2*(dphi2/dx)*phij> j=1,2,3
         DDY2=FDY2/3.d0         !<2*(dphi2/dy)*phij> j=1,2,3
         DXX21=DXX12            !<2*(dphi2/dx)*(dphi1/dx)>
         DYY21=DYY12            !<2*(dphi2/dy)*(dphi1/dy)>
         DXY21=FDX2O2A*FDY1      !<2*(dphi2/dx)*(dphi1/dy)>
         DXXYY21=DXXYY12
         DXYH21=DXYH12
         DXX22=FDX2O2A*FDX2      !<2*(dphi2/dx)*(dphi2/dx)>
         DYY22=FDY2O2A*FDY2      !<2*(dphi2/dy)*(dphi2/dy)>
         DXY22=FDX2O2A*FDY2      !<2*(dphi2/dx)*(dphi2/dy)>
         DXXYY22=DXX22+DYY22
         DXYH22=DPAvg*DXXYY22      
         DXX23=FDX2O2A*FDX3      !<2*(dphi2/dx)*(dphi3/dx)>
         DYY23=FDY2O2A*FDY3      !<2*(dphi2/dy)*(dphi3/dy)>
         DXY23=FDX2O2A*FDY3      !<2*(dphi2/dx)*(dphi3/dy)>
         DXXYY23=DXX23+DYY23
         DXYH23=DPAvg*DXXYY23      

         DDX3=FDX3/3.d0         !<2*(dphi3/dx)*phij> j=1,2,3
         DDY3=FDY3/3.d0         !<2*(dphi3/dy)*phij> j=1,2,3
         DXX31=DXX13            !<2*(dphi3/dx)*(dphi1/dx)>
         DYY31=DYY13            !<2*(dphi3/dy)*(dphi1/dy)>
         DXY31=FDX3O2A*FDY1      !<2*(dphi3/dx)*(dphi1/dy)>
         DXXYY31=DXXYY13
         DXYH31=DXYH13
         DXX32=DXX23            !<2*(dphi3/dx)*(dphi2/dx)>
         DYY32=DYY23            !<2*(dphi3/dy)*(dphi2/dy)>
         DXY32=FDX3O2A*FDY2      !<2*(dphi3/dx)*(dphi2/dy)>
         DXXYY32=DXXYY23
         DXYH32=DXYH23
         DXX33=FDX3O2A*FDX3      !<2*(dphi3/dx)*(dphi3/dx)>
         DYY33=FDY3O2A*FDY3      !<2*(dphi3/dy)*(dphi3/dy)>
         DXY33=FDX3O2A*FDY3      !<2*(dphi3/dx)*(dphi3/dy)>
         DXXYY33=DXX33+DYY33
         DXYH33=DPAvg*DXXYY33      

         AO6=AREAIE2/6.d0
         AO12=AREAIE2/12.d0
         FDDD=(1+ILUMP)*AO6     !2*<phi*phj> diagonal terms
         FDDOD=(1-ILUMP)*AO12   !2*<phi*phj> off diagonal terms

C
C...  COMPUTE THE CONTRIBUTION OF A SPATIALLY VARIABLE Tau0
C 
        SpaVarTau0=(T0N1*DDX1+T0N2*DDX2+T0N3*DDX3)*QX1Avg
     &            +(T0N1*DDY1+T0N2*DDY2+T0N3*DDY3)*QY1Avg

C...  
C...  COMPUTE THE RHS GWCE FORCING AND PUT INTO TEMP_LV_A VECTOR FOR NODE NM1
C...  
         TEMP_LV_A1=
C...  TRANSIENT AND Tau0 TERMS FROM LHS
     &        -(FDDD*ESN1+FDDOD*ESN2+FDDOD*ESN3)*TT0R

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K-1)
     &        -(DXYH11*E0N1+DXYH12*E0N2+DXYH13*E0N3)*GC00

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K)
     &        -(DXYH11*E1N1+DXYH12*E1N2+DXYH13*E1N3)*GB00A00

C...  BOTTOM FRICTION & Tau0
     &        +(T0XAvg3-BSXAvg3)*DDX1
     &        +(T0YAvg3-BSYAvg3)*DDY1
     &        +SpaVarTau0

C...  CORIOLIS FORCE
     &        +CorifAvg*(QY1Avg3*DDX1-QX1Avg3*DDY1)

C...  WIND STRESS & ATMONPHERIC PRESSURE GRADIENTS
     &        +WSXAvg3*DDX1+WSYAvg3*DDY1
     &        -GHAvg*(PR1N1*DXXYY11+PR1N2*DXXYY12+PR1N3*DXXYY13)

C...  TIDAL POTENTIAL FORCING
     &        +GHAvg*(TIPN1*DXXYY11+TIPN2*DXXYY12+TIPN3*DXXYY13)

C...  LATERAL VISCOSITY TERM  
     &        -EVMAvgODT*(DXXYY11*ESN1+DXXYY12*ESN2+DXXYY13*ESN3)

C...  FINITE AMPLITUDE
     &        -GFAO2*(E1N1SQ*DXXYY11+E1N2SQ*DXXYY12+E1N3SQ*DXXYY13)

C... SPATIAL DERIVATIVE PORTION OF THE ADVECTIVE TERMS
     &        -IFNLCT*(QX1Avg*(U1N1*DXX11+U1N2*DXX21+U1N3*DXX31
     &        +                V1N1*DXY11+V1N2*DXY21+V1N3*DXY31)
     &        +        QY1Avg*(U1N1*DXY11+U1N2*DXY12+U1N3*DXY13
     &        +                V1N1*DYY11+V1N2*DYY21+V1N3*DYY31))

C...  Advective terms (time derivative portion in gwce) which must be
c...  bundled in with the finite amplitude terms in order to get good
c...  mass conservation when the advective terms are shut down
     &        +TADVODT*(U1Avg*DDX1+V1Avg*DDY1)*(ESN1+ESN2+ESN3)

C...  DENSITY TERMS 
     &        -VIDBCPDXAvg3*DDX1-VIDBCPDYAvg3*DDY1

C     3D.3D Velocity dispersion
         if (C3D) then
           TEMP_LV_A1=TEMP_LV_A1-IFNLCT*
     &           (DUU1N1*DXX11+DUU1N2*DXX12+DUU1N3*DXX13
     &           +DUV1N1*DXY11+DUV1N2*DXY12+DUV1N3*DXY13
     &           +DUV1N1*DXY11+DUV1N2*DXY21+DUV1N3*DXY31
     &           +DVV1N1*DYY11+DVV1N2*DYY12+DVV1N3*DYY13)
           endif

C...  
C...  COMPUTE THE RHS GWCE FORCING AND PUT INTO TEMP_LV_A VECTOR FOR NODE NM2
C...  
         TEMP_LV_A2=
C...  TRANSIENT AND Tau0 TERMS FROM LHS
     &        -(FDDOD*ESN1+FDDD*ESN2+FDDOD*ESN3)*TT0R

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K-1)
     &        -(DXYH12*E0N1+DXYH22*E0N2+DXYH23*E0N3)*GC00

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K)
     &        -(DXYH12*E1N1+DXYH22*E1N2+DXYH23*E1N3)*GB00A00

C...  BOTTOM FRICTION & Tau0
     &        +(T0XAvg3-BSXAvg3)*DDX2
     &        +(T0YAvg3-BSYAvg3)*DDY2
     &        +SpaVarTau0

C...  CORIOLIS FORCE
     &        +CorifAvg*(QY1Avg3*DDX2-QX1Avg3*DDY2)

C...  WIND AND ATMOSPHERIC PRESSURE FORCING
     &        +WSXAvg3*DDX2+WSYAvg3*DDY2
     &        -GHAvg*(PR1N1*DXXYY21+PR1N2*DXXYY22+PR1N3*DXXYY23)

C...  TIDAL POTENTIAL FORCING
     &        +GHAvg*(TIPN1*DXXYY21+TIPN2*DXXYY22+TIPN3*DXXYY23)

C...  LATERAL VISCOSITY TERM
     &        -EVMAvgODT*(DXXYY12*ESN1+DXXYY22*ESN2+DXXYY23*ESN3)

C...  FINITE AMPLITUDE
     &        -GFAO2*(E1N1SQ*DXXYY21+E1N2SQ*DXXYY22+E1N3SQ*DXXYY23)

C...  SPATIAL DERIVATIVE PORTION OF THE ADVECTIVE TERMS
     &        -IFNLCT*(QX1Avg*(U1N1*DXX12+U1N2*DXX22+U1N3*DXX32
     &        +                V1N1*DXY12+V1N2*DXY22+V1N3*DXY32)
     &        +        QY1Avg*(U1N1*DXY21+U1N2*DXY22+U1N3*DXY23
     &        +                V1N1*DYY12+V1N2*DYY22+V1N3*DYY32))

C...  Advective terms (time derivative portion in gwce) which must be
c...  bundled in with the finite amplitude terms in order to get good
c...  mass conservation when the advective terms are shut down - non-conservative

     &        +TADVODT*(U1Avg*DDX2+V1Avg*DDY2)*(ESN1+ESN2+ESN3)

C...  DENSITY TERMS 
     &        -VIDBCPDXAvg3*DDX2-VIDBCPDYAvg3*DDY2

C     3D.3D Velocity dispersion
         if (C3D) then
           TEMP_LV_A2=TEMP_LV_A2-IFNLCT*
     &           (DUU1N1*DXX12+DUU1N2*DXX22+DUU1N3*DXX23
     &           +DUV1N1*DXY21+DUV1N2*DXY22+DUV1N3*DXY23
     &           +DUV1N1*DXY12+DUV1N2*DXY22+DUV1N3*DXY32
     &           +DVV1N1*DYY12+DVV1N2*DYY22+DVV1N3*DYY23)
           endif

C...  
C...  COMPUTE THE RHS GWCE FORCING AND PUT INTO TEMP_LV_A VECTOR FOR NODE NM1
C...  
         TEMP_LV_A3=
C...  TRANSIENT AND Tau0 TERMS FROM LHS
     &        -(FDDOD*ESN1+FDDOD*ESN2+FDDD*ESN3)*TT0R

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K-1)
     &        -(DXYH13*E0N1+DXYH23*E0N2+DXYH33*E0N3)*GC00

C...  FREE SURFACE TERMS FROM LHS (TIME LEVEL K)
     &        -(DXYH13*E1N1+DXYH23*E1N2+DXYH33*E1N3)*GB00A00

C...  BOTTOM FRICTION & Tau0
     &        +(T0XAvg3-BSXAvg3)*DDX3
     &        +(T0YAvg3-BSYAvg3)*DDY3
     &        +SpaVarTau0

C...  CORIOLIS FORCE
     &        +CorifAvg*(QY1Avg3*DDX3-QX1Avg3*DDY3)

C...  WIND AND ATMOSPHERIC PRESSURE FORCING
     &        +WSXAvg3*DDX3+WSYAvg3*DDY3
     &        -GHAvg*(PR1N1*DXXYY31+PR1N2*DXXYY32+PR1N3*DXXYY33)

C...  TIDAL POTENTIAL FORCING
     &        +GHAvg*(TIPN1*DXXYY31+TIPN2*DXXYY32+TIPN3*DXXYY33)

C...  LATERAL VISCOSITY TERM
     &        -EVMAvgODT*(DXXYY13*ESN1+DXXYY23*ESN2+DXXYY33*ESN3)

C...  FINITE AMPLITUDE
     &        -GFAO2*(E1N1SQ*DXXYY31+E1N2SQ*DXXYY32+E1N3SQ*DXXYY33)

C...  SPATIAL DERIVATIVE PORTION OF THE ADVECTIVE TERMS
     &        -IFNLCT*(QX1Avg*(U1N1*DXX13+U1N2*DXX23+U1N3*DXX33
     &        +                V1N1*DXY13+V1N2*DXY23+V1N3*DXY33)
     &        +        QY1Avg*(U1N1*DXY31+U1N2*DXY32+U1N3*DXY33
     &        +                V1N1*DYY13+V1N2*DYY23+V1N3*DYY33))

C...  Advective terms (time derivative portion in gwce) which must be
c...  bundled in with the finite amplitude terms in order to get good
c...  mass conservation when the advective terms are shut down - non-conservative

     &        +TADVODT*(U1Avg*DDX3+V1Avg*DDY3)*(ESN1+ESN2+ESN3)
           
C...  DENSITY TERMS 
     &        -VIDBCPDXAvg3*DDX3-VIDBCPDYAvg3*DDY3

C     3D.3D Velocity dispersion
         if (C3D) then
           TEMP_LV_A3=TEMP_LV_A3-IFNLCT*
     &           (DUU1N1*DXX13+DUU1N2*DXX23+DUU1N3*DXX33
     &           +DUV1N1*DXY31+DUV1N2*DXY32+DUV1N3*DXY33
     &           +DUV1N1*DXY13+DUV1N2*DXY23+DUV1N3*DXY33
     &           +DVV1N1*DYY13+DVV1N2*DYY23+DVV1N3*DYY33)
           endif

C     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
         TEMP_LV_A(IE,1)=TEMP_LV_A1*NCELE
         TEMP_LV_A(IE,2)=TEMP_LV_A2*NCELE
         TEMP_LV_A(IE,3)=TEMP_LV_A3*NCELE
#endif

C     LINES TO RUN ON A SCALAR COMPUTER.
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR GWCE_LV
C     ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
         GWCE_LV(NM1)=GWCE_LV(NM1)+TEMP_LV_A1*NCELE
         GWCE_LV(NM2)=GWCE_LV(NM2)+TEMP_LV_A2*NCELE
         GWCE_LV(NM3)=GWCE_LV(NM3)+TEMP_LV_A3*NCELE
#endif

c        IF(IE.EQ.1) THEN
c           WRITE(101,*) ' '
c           WRITE(101,*) '  ************ GWCE Load Vector/2 ***********'
c           WRITE(101,*) '  Time Step = ',IT
c           ENDIF
c        WRITE(101,*) IE, ESN1, ESN2, ESN3
c        WRITE(101,*) IE, TEMP_LV_A1/2., TEMP_LV_A2/2., TEMP_LV_A3/2.
 
 1037 CONTINUE

C     LINES TO RUN ON A VECTOR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR GWCE_LV
C     ON A VECTOR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         GWCE_LV(NM1)=GWCE_LV(NM1)+TEMP_LV_A(IE,1)
         GWCE_LV(NM2)=GWCE_LV(NM2)+TEMP_LV_A(IE,2)
         GWCE_LV(NM3)=GWCE_LV(NM3)+TEMP_LV_A(IE,3)
      END DO
#endif

C...  
C...  SAVE THE ELEVATION AT THE PAST TIME STEP INTO ETA1 AND ZERO ETA2
C...  
      DO I=1,NP
         ETA1(I)=ETA2(I)
      END DO
      
C...  AT ELEVATION BOUNDARY CONDITION NODES DETERMINE ELEVATION AT NEXT
C...  TIME STEP
C...  
C...  FOR PERIODIC ELEVATION BOUNDARY CONDITION

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0.
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            ETA2(NBDI)=ETA2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO
C
C...  FOR APERIODIC ELEVATION BOUNDARY CONDITION
      IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TimeLoc-ETIME1)/ETIMINC
         DO I=1,NETA
            NBDI=NBD(I)
            ETA2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
         END DO
      ENDIF
C
C     jgf46.02 Added the ability to include geoid offset on the boundary.
      IF (LoadGeoidOffset) THEN
         DO I=1,NETA
            ETA2(NBD(I))=ETA2(NBD(I))+GeoidOffset(NBD(I))
         END DO
      ENDIF

C   kmd48.33bc add information for the levels of no motion boundary conditions
C              these are considered the steric adjustments. 
      IF ((ABS(RES_BC_FLAG).GE.1).AND.(CBaroclinic))  THEN
         DO I=1,NETA
           NBDI=NBD(I)
           ETA2(NBDI) = ETA2(NBDI) + LNM_BC(I)
         END DO
      END IF

C
C...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
C...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

C...  Note 1, these values all must be multiplied by 2 since all
C...  elemental coefficients have been

C...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
C...  assume that QN is positive into the domain.  QFORCEJ has a -1
C...  built in and the terms are not explicitly negated. Boundary
C...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
C...  flux that  is positive out of the domain.  Therefore, to match
C...  the formulation these terms must be explicitly multiplied by -1.

C...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN 
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 + 
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN 
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)

            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN 
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0

            NCI=NODECODE(NBDI)
            NCJ=NODECODE(NBDJ)
            NCBND=NCI*NCJ
            BNDLEN2O3NC=NCBND*BNDLEN2O3(J-1)
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BNDLEN2O3NC*(QFORCEI+QFORCEJ/2.D0)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BNDLEN2O3NC*(QFORCEJ+QFORCEI/2.D0)
            END DO
        ENDIF

C...
C...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
C...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
C...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
C...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
C...  EQUATIONS
C...  
      DO I=1,NETA
         NBDI=NBD(I)
         ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
         GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
         DO J=2,NNEIGH(NBDI)
            GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &           -ETAS(NBDI)*OBCCOEF(I,J-1)
         END DO
      END DO

C...  
C...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
C...  

C...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
C...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
#endif

C...  JCG ITERATIVE MATRIX SOLVER
      IPARM(1)=ITMAX
      CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &     IWKSP,NW,WKSP,IPARM,RPARM,IER)

      NUMITR=IPARM(1)
      DO I=1,NP
         ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
      END DO
     
C     UPDATE ELEVATIONS

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif

      RETURN
      END SUBROUTINE GWCE_ORIGINAL



C******************************************************************************
C                                                                             *
C    Subroutine to compute the elevation using the GWCE formluation           *
C    Re-written to conform to the ADCIRC Theory Report                        *
C                                                                             *
C                            r.l.  06/22/2005                                 *
C******************************************************************************

      SUBROUTINE GWCE_New(IT,TimeLoc,TimeH)
C     
      USE GLOBAL
      USE SIZES
      USE WIND
      USE ITPACKV
      USE NodalAttributes, ONLY : 
     &     LoadGeoidOffset, GeoidOffset, EVM, 
     &     TAU0VAR, HighResTimeVaryingTau0, FullDomainTimeVaryingTau0,
     &     CalculateTimeVaryingTau0
#ifdef CMPI
      USE MESSENGER
#endif
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      USE Couple2Swan, ONLY: TKXX,
     &                       TKXY,
     &                       TKYY
#endif
#endif

      IMPLICIT NONE 

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      LOGICAL  DIE
      REAL(SZ) A00pB00
      REAL(SZ) BCXAvg, BCYAvg
      REAL(SZ) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(SZ) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(SZ) CorifAvg
      REAL(SZ) DPAvg, GDPAvgOAreaIE4
      REAL(SZ) DispX, DispY, DispXAvg, DispYAvg
      REAL(SZ) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(SZ) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(SZ) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(SZ) ESN1, ESN2, ESN3, ESAvg
      REAL(SZ) EVMH, EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(SZ) EVMEle, EVMSmag
      REAL(SZ) GA00DPAvgOAreaIE4
      REAL(SZ) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(SZ) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg, H1, H2
      REAL(SZ) H2OTotalArea
      REAL(SZ) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(SZ) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(SZ) MsFacR, MsFacLOnDiag, MsFacLOffDiag
      REAL(SZ) MX, MY, MXAvg, MYAvg
      REAL(SZ) JXAvg, JYAvg
      REAL(SZ) Pr1N1, Pr1N2, Pr1N3
      REAL(SZ) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(SZ) SFacAvg
      REAL(SZ) T0N1,T0N2, T0N3
      REAL(SZ) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(SZ) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TiPN1, TiPN2, TiPN3
      REAL(SZ) U1N1,U1N2,U1N3, U1Avg
      REAL(SZ) V1N1,V1N2,V1N3, V1Avg
      REAL(SZ) WSXAvg, WSYAvg
      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) TimeLoc, TimeH

      REAL(SZ) HH1 !jgf46.02 Added for Katrina.
      REAL(SZ) RDIAG ! jgf48.4619 Seizo parameter for fully explicit mode
Casey 050711 : Added for averaged variable Tau0.
C      REAL(SZ) :: CaseySum
C      REAL(SZ), ALLOCATABLE :: TAU0VARTEMP(:)

C     jgf45.11 Bug fix: calculate the integers OnDiag and OffDiag here
C     instead of inside the GWCE lhs (system matrix) setup, since they
C     are also used in the calculation of the GWCE load vector gwce_lv.

C...  Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1
      OnDiag=(1+ILump)*2        !diagonal coefficient
      OffDiag=(1-ILump)         !off diagonal coefficient
C...
C...  Recompute the GWCE system matrix at the first time step or if any
C...  wetting or drying occurred in the previous time step.
C...  
      IF(NCChange.GT.0) THEN !if any subdomain grid has changed
         NCChange=0
c        IF(NScreen.GT.0.AND.MYProc.EQ.0) WRITE(screenunit,3806) 
c        WRITE(16,3806)                                 
 3806    FORMAT(/,1X,'RE-SETTING GWCE SYSTEM MATRIX',/)

C.....Set up the LHS matrix (for the iterative matrix solver)
         IF ( ILump.eq.0 ) THEN ! default, fully consistent case
            Coef(:,:)=0.0d0
         ELSE ! jgf48.4619: ILump.eq.1, only need the diagonals (Seizo)
            Coefd(:)=0.0d0    ! Only Diagnal
         ENDIF
C
C        jgf47.08 Moved time-varying tau0 subroutine to nodalattr.F
         IF(C2DDI.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
               TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         IF(C3D.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
               TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         
         
         !jgf48.4619: Added Seizo's handling of Lumping vs non-Lumping
         IF ( ILump.eq.0 ) THEN ! default case: fully consistent LHS
            DO IE=1,NE
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg=(SFac(NMI1)+SFac(NMI2)+SFac(NMI3))/3.d0

               FDX1 = (Y(NMI2)-Y(NMI3))*SFacAvg         !b1 = 2*Area*dphi1/dx
               FDX2 = (Y(NMI3)-Y(NMI1))*SFacAvg         !b2 = 2*Area*dphi2/dx
               FDX3 = (Y(NMI1)-Y(NMI2))*SFacAvg         !b3 = 2*Area*dphi3/dx
               FDY1 = X(NMI3)-X(NMI2)                   !a1 = 2*Area*dphi1/dy
               FDY2 = X(NMI1)-X(NMI3)                   !a2 = 2*Area*dphi2/dy
               FDY3 = X(NMI2)-X(NMI1)                   !a3 = 2*Area*dphi3/dy

               AreaIE2=Areas(IE)
               AreaIE =AreaIE2/2.0d0
               AreaIE4=AreaIE2*2.0d0

               DPAvg=(DP(NMI1)+DP(NMI2)+DP(NMI3))/3.d0
               GA00DPAvgOAreaIE4=G*A00*DPAvg/AreaIE4
               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0
               MsFacLOnDiag =OnDiag *AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0
               MsFacLOffDiag=OffDiag*AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0

               DO JN=2,NEIMAX
                  IF(NeiTab(NMI1,JN).EQ.NMJ2) J12=JN
                  IF(NeiTab(NMI1,JN).EQ.NMJ3) J13=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ1) J21=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ3) J23=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ1) J31=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ2) J32=JN
               END DO

               Coef(NMI1,1)  =Coef(NMI1,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX1+FDY1*FDY1))*NCELE
               Coef(NMI1,J12)=Coef(NMI1,J12) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX2+FDY1*FDY2))*NCELE
               Coef(NMI1,J13)=Coef(NMI1,J13) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX3+FDY1*FDY3))*NCELE
               Coef(NMI2,J21)=Coef(NMI2,J21) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX1+FDY2*FDY1))*NCELE
               Coef(NMI2,1)  =Coef(NMI2,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX2+FDY2*FDY2))*NCELE
               Coef(NMI2,J23)=Coef(NMI2,J23) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX3+FDY2*FDY3))*NCELE
               Coef(NMI3,J31)=Coef(NMI3,J31) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX1+FDY3*FDY1))*NCELE
               Coef(NMI3,J32)=Coef(NMI3,J32) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX2+FDY3*FDY2))*NCELE
               Coef(NMI3,1)  =Coef(NMI3,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX3+FDY3*FDY3))*NCELE

            ENDDO
         ELSE
            !jgf48.4619: Add Seizo's construction of Lumped LHS matrix
            DO IE=1,NE ! Make LHS Lumped Matrix
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg=(SFac(NMI1)+SFac(NMI2)+SFac(NMI3))/3.d0

               AreaIE2=Areas(IE)
               AreaIE =AreaIE2/2.0d0

               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0
               MsFacLOnDiag =OnDiag *AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0

               Coefd(NMI1)  =Coefd(NMI1)   + (MsFacLOnDiag)*NCELE
               Coefd(NMI2)  =Coefd(NMI2)   + (MsFacLOnDiag)*NCELE
               Coefd(NMI3)  =Coefd(NMI3)   + (MsFacLOnDiag)*NCELE
            ENDDO          
         ENDIF

C...  Modify the matrix "COEF" by imposing the elevation specified
C...  boundary conditions while maintaining the symmetry of the system

         IF (ILump.eq.0) THEN    
#ifdef CMPI
            EP = PSDOT(NP,Coef(1,1),Coef(1,1))
            EP = SQRT(RNP_GLOBAL*EP)
#else 
            EP=0.0D0
            DO I=1,NP
               EP=EP+Coef(I,1)*Coef(I,1)
            ENDDO
            EP=SQRT(EP/NP)
#endif
c...        for each elevation specified boundary node, zero all off diagonal
c...        terms on the row and set diagnoal term to EP
            DO I=1,NETA
               Coef(NBD(I),1)=EP
               DO J=2,NNEIGH(NBD(I))
                  Coef(NBD(I),J)=0.0d0
               ENDDO
            ENDDO
c...        for each elevation specified boundary node, zero all off diagonal
c...        terms on the column but save these to be multiplied by the
c...        boundary value and subtracted from the RHS
            DO I=1,NETA
               DO J=2,NNeigh(NBD(I))
                  DO IJ=2,NNeigh(NeiTab(NBD(I),J))
                     IF(NBD(I).EQ.NeiTab(NeiTab(NBD(I),J),IJ)) THEN
                        OBCCoef(I,J-1)=Coef(NeiTab(NBD(I),J),IJ)
                        Coef(NeiTab(NBD(I),J),IJ)=0.0d0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
C.....      Check that all the diagonal elements in "COEF" are > 0.
            DIE = .FALSE.
            DO I=1,NP
               IF(COEF(I,1).EQ.0.d0) COEF(I,1)=EP
               IF(COEF(I,1).LT.0.d0) THEN
                  IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              WRITE(ScreenUnit,1019) I,COEF(I,1)
                  WRITE(16,1019) I,COEF(I,1)
                  DIE = .TRUE.
               ENDIF
            ENDDO
            IF (DIE) STOP

         ELSE  ! jgf48.4619: include Seizo's changes for lumped LHS

            ! Seizo: Explicit scheme can solve localy. (the efect is small?)
            EP=0.0D0
            DO I=1,NP
               EP=EP+Coefd(I)*Coefd(I)
            ENDDO
            EP=SQRT(EP/NP)

            ! set diagonal term to EP
            DO I=1,NETA
               Coefd(NBD(I))=EP
            ENDDO

C.....      Check that all the diagonal elements in "COEFD" are > 0.
            DIE = .FALSE.
            DO I=1,NP
               IF(COEFD(I).EQ.0.d0) COEFD(I)=EP
               IF(COEFD(I).LT.0.d0) THEN
                  IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              WRITE(ScreenUnit,1019) I,COEFD(I)
                  WRITE(16,1019) I,COEFD(I)
 1019             FORMAT(/,1X,'!!!!!!!!  ERROR  !!!!!!!',
     &              /,1X,'THE DIAGONAL TERM IN THE EQUATION FOR NODE ',I10,
     &                '= ',E15.6,' AND IS <= 0',/)
                  DIE = .TRUE.
               ENDIF
            ENDDO
            IF (DIE) STOP
         ENDIF
        
      ENDIF                     !End of GWCE matrix setup


C...  
C...  Compute the GWCE load vector GWCE_LV
C...  This is done primarily element by element by forming 
C...  temporary vectors and then assembling at the end.  
C...  This has been set up to unroll loops to optimize performance
C...  on vector processors.  
C...
C...  Elevation and flux boundary conditions are imposed after the
C...  element by element assembly section.  
C...

C...  Initialize variables to zero if these forcings are not used

      IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
      ELSE
        WSXAvg=0.d0
        WSYAvg=0.d0
        Pr1N1=0.d0
        Pr1N2=0.d0
        Pr1N3=0.d0
      ENDIF

      IF (CTIP) THEN
      ELSE
        TiPN1=0.d0
        TiPN2=0.d0
        TiPN3=0.d0
      ENDIF

      IF(C3D) THEN
      ELSE
        DispXAvg=0.d0
        DispYAvg=0.d0
      ENDIF

      IF(CBaroclinic) THEN
      ELSE
        BCXAvg=0.d0
        BCYAvg=0.d0
      ENDIF

C...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartV).OR.(CGWCE_LS_2PartSV)) THEN

         DO I=1,NP
            LSXX(I)=0.d0
            LSXY(I)=0.d0
            LSYX(I)=0.d0
            LSYY(I)=0.d0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0   !A*DUDX
            LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0   !A*DUDY
            LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0   !A*DVDX
            LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0   !A*DVDY
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               EVMEle=EVMSmag
            ENDIF
            LSXXEle = LSXXGradA*EVMEle            
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
               IF (CGWCE_LS_2PartV) THEN          !nonsymmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=H2OTotalArea*LSXY(I)
                  LSYX(I)=H2OTotalArea*LSYX(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
               IF (CGWCE_LS_2PartSV) THEN         !symmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=0.5d0*H2OTotalArea*(LSXY(I)+LSYX(I))
                  LSYX(I)=LSXY(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
            ELSE
               LSXX(I)=0.d0
               LSXY(I)=0.d0
               LSYX(I)=0.d0
               LSYY(I)=0.d0
            ENDIF
         ENDDO

      ENDIF

C...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartSQ)) THEN

         DO I=1,NP
            LSXX(I)=0.d0
            LSXY(I)=0.d0
            LSYX(I)=0.d0
            LSYY(I)=0.d0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
              LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
              LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
              LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
              LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0
              EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
              EVMEle=EVMSmag
            ENDIF
            LSXXGradA=(QX1(NM1)*FDX1+QX1(NM2)*FDX2+QX1(NM3)*FDX3)/2.d0
            LSXYGradA=(QX1(NM1)*FDY1+QX1(NM2)*FDY2+QX1(NM3)*FDY3)/2.d0
            LSYXGradA=(QY1(NM1)*FDX1+QY1(NM2)*FDX2+QY1(NM3)*FDX3)/2.d0
            LSYYGradA=(QY1(NM1)*FDY1+QY1(NM2)*FDY2+QY1(NM3)*FDY3)/2.d0
            LSXXEle = LSXXGradA*EVMEle            
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               IF (CGWCE_LS_2PartQ) THEN          !nonsymmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=LSXY(I)/TotalArea(I)
                  LSYX(I)=LSYX(I)/TotalArea(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
               IF (CGWCE_LS_2PartSQ) THEN         !symmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=0.5d0*(LSXY(I)+LSYX(I))/TotalArea(I)
                  LSYX(I)=LSXY(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
            ELSE
               LSXX(I)=0.d0
               LSXY(I)=0.d0
               LSYX(I)=0.d0
               LSYY(I)=0.d0
            ENDIF
         ENDDO

      ENDIF

C...  Assemble the GWCE RHS except for the boundary integral terms

      DO 1037 IE=1,NE
  
C...     Set nodal values for each element
  
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         E0N1=ETA1(NM1)
         E0N2=ETA1(NM2)
         E0N3=ETA1(NM3)
         E1N1=ETA2(NM1)
         E1N2=ETA2(NM2)
         E1N3=ETA2(NM3)
         E1N1SQ=E1N1*E1N1
         E1N2SQ=E1N2*E1N2
         E1N3SQ=E1N3*E1N3
         ESN1=ETAS(NM1)
         ESN2=ETAS(NM2)
         ESN3=ETAS(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         T0N1=Tau0Var(NM1) 
         T0N2=Tau0Var(NM2) 
         T0N3=Tau0Var(NM3) 

         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            Pr1N1=PR1(NM1)
            Pr1N2=PR1(NM2)
            Pr1N3=PR1(NM3)
         ENDIF
         IF (CTIP) THEN                        !tidal potential
            TiPN1=TiP1(NM1)
            TiPN2=TiP1(NM2)
            TiPN3=TiP1(NM3)
         ENDIF
         IF (C2DDI) THEN                       !2D bottom friction
            BSXN1=TK(NM1)*QX1N1
            BSYN1=TK(NM1)*QY1N1
            BSXN2=TK(NM2)*QX1N2
            BSYN2=TK(NM2)*QY1N2
            BSXN3=TK(NM3)*QX1N3
            BSYN3=TK(NM3)*QY1N3
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            BSXN1 = TKXX(NM1)*QX1N1 + TKXY(NM1)*QY1N1
            BSYN1 = TKXY(NM1)*QX1N1 + TKYY(NM1)*QY1N1
            BSXN2 = TKXX(NM2)*QX1N2 + TKXY(NM2)*QY1N2
            BSYN2 = TKXY(NM2)*QX1N2 + TKYY(NM2)*QY1N2
            BSXN3 = TKXX(NM3)*QX1N3 + TKXY(NM3)*QY1N3
            BSYN3 = TKXY(NM3)*QX1N3 + TKYY(NM3)*QY1N3
#endif
#endif
         ENDIF      
         IF (C3D) THEN                         !3D bottom friction
            BSXN1=BSX1(NM1)
            BSXN2=BSX1(NM2)
            BSXN3=BSX1(NM3)
            BSYN1=BSY1(NM1)
            BSYN2=BSY1(NM2)
            BSYN3=BSY1(NM3)
         ENDIF

         AreaIE2=Areas(IE)               !2A
         AreaIE=AreaIE2/2.d0             ! A
         AreaIE4=2.d0*AreaIE2            !4A

         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1 = 2*Area*dphi1/dx
         FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2 = 2*Area*dphi2/dx
         FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3 = 2*Area*dphi3/dx
         FDY1 = X(NM3)-X(NM2)           !a1 = 2*Area*dphi1/dy
         FDY2 = X(NM1)-X(NM3)           !a2 = 2*Area*dphi2/dy
         FDY3 = X(NM2)-X(NM1)           !a3 = 2*Area*dphi3/dy

C...     Compute part of several spatial gradients for use below

         E0XGrad2A = 0.0d0
         E0YGrad2A = 0.0d0
         IF (ILump.eq.0) THEN
            E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
            E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
         ENDIF
         E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
         E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
         Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
         Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

C...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM

         IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.d0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &         +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)
     &         -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)
     &         +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)
     &         -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
         ENDIF

C...     Compute the remainder of the 2 Part lateral stress terms

         IF((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartV) .OR.
     &        (CGWCE_LS_2PartSQ).OR.(CGWCE_LS_2PartSV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
         ENDIF

C...     Compute the spatial gradients of the velocity dispersion terms if 3D

         IF (C3D) THEN                         !3D bottom friction
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
         ENDIF

C...     Compute elemental averages                   

         CorifAvg=(Corif(NM1)+Corif(NM2)+Corif(NM3))/3.d0
         Tau0Avg=(T0N1+T0N2+T0N3)/3.d0
         Tau0QXAvg=(T0N1*QX1N1+T0N2*QX1N2+T0N3*QX1N3)/3.d0
         Tau0QYAvg=(T0N1*QY1N1+T0N2*QY1N2+T0N3*QY1N3)/3.d0
         U1Avg=(U1N1+U1N2+U1N3)/3.d0
         V1Avg=(V1N1+V1N2+V1N3)/3.d0
         QX1Avg=(QX1N1+QX1N2+QX1N3)/3.d0
         QY1Avg=(QY1N1+QY1N2+QY1N3)/3.d0
         ESAvg=(ESN1+ESN2+ESN3)/3.d0
         DPAvg=(DP(NM1)+DP(NM2)+DP(NM3))/3.d0
         GDPAvgOAreaIE4=G*DPAvg/AreaIE4
         HAvg=(H1N1+H1N2+H1N3)/3.d0
         GHAvg=G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.d0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.d0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
         ENDIF
         IF (C3D) THEN                !3D velocity dispersion
            DispXAvg=IFNLCT*DispX
            DispYAvg=IFNLCT*DispY
         ENDIF
         IF(CBaroclinic) THEN
            BCXAvg=(H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                       +H1N3*VIDBCPDXOH(NM3))/3.d0
            BCYAvg=(H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                       +H1N3*VIDBCPDYOH(NM3))/3.d0
         ENDIF

C...     Compute additional partial factors

         MsFacR=AreaIE*(1.d0/DT-Tau0Avg/2.d0)/DT/12.d0
         GOAreaIE4=G/AreaIE4
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
         A00pB00=A00+B00

C...     Compute the JX, JY terms less the advection terms

         JXAvg = CorifAvg*QY1Avg
     &          -IFNLFA*GOAreaIE4*(E1N1SQ*FDX1+E1N2SQ*FDX2
     &                                        +E1N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1
     &                     +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3)
     &          +WSXAvg-BSXAvg+MXAvg-DispXAvg-BCXAvg+Tau0QXAvg 

         JYAvg =-CorifAvg*QX1Avg
     &          -IFNLFA*GOAreaIE4*(E1N1SQ*FDY1+E1N2SQ*FDY2
     &                                        +E1N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1
     &                    +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3)
     &          +WSYAvg-BSYAvg+MYAvg-DispYAvg-BCYAvg+Tau0QYAvg

C...     Complete the JX, JY terms depending on the advection formulation

         IF(CGWCE_Advec_NC) THEN        !nonconservative advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2
     &             +IFNLCAT*U1Avg*ESAvg/DT
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2
     &             +IFNLCAT*V1Avg*ESAvg/DT
         ENDIF
         IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
           JXAvg = JXAvg - IFNLCT*(
     &              (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2
     &                              +U1N3*QX1N3*FDX3)
     &             +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2
     &                              +U1N3*QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2
     &                              +V1N3*QX1N3*FDX3)
     &             +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2
     &                              +V1N3*QY1N3*FDY3))/AreaIE2
         ENDIF
         IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3)
     &             +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3)
     &             +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
         ENDIF

C...     Assemble forcing for node NM1 (local index j=1) 

         Temp_LV_A1=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN1 + OffDiag*(ESN2+ESN3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX1*E0XGrad2A+FDY1*E0YGrad2A)
     &                        +A00pB00*(FDX1*E1XGrad2A+FDY1*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX1 + JYAvg*FDY1)/2.d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar

C...     Assemble forcing for node NM2 (local index j=2) 

         Temp_LV_A2=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN2 + OffDiag*(ESN1+ESN3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX2*E0XGrad2A+FDY2*E0YGrad2A)
     &                        +A00pB00*(FDX2*E1XGrad2A+FDY2*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX2 + JYAvg*FDY2)/2.d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar


C...     Assemble forcing for node NM3 (local index j=3) 

         Temp_LV_A3=

C...     Transient and Tau0 terms from LHS
C...    (consistent mass matrix: ILump=0, lumped mass matrix: ILump=1)
     &         (OnDiag*ESN3 + OffDiag*(ESN1+ESN2))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX3*E0XGrad2A+FDY3*E0YGrad2A)
     &                        +A00pB00*(FDX3*E1XGrad2A+FDY3*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX3 + JYAvg*FDY3)/2.d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar


C...     Put these partial products into further elemental storage for a vector computer
C...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
         Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
         Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
         Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

C...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif
            
 1037 CONTINUE                  !End of elemental loop

C...  Put load vector elemental values into nodal storage for a vector computer
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A(IE,1)
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A(IE,2)
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A(IE,3)
      END DO
#endif

C     jgfdebug4608 
C      if (it.eq.16) then
C         open(unit=99,file=trim(localdir)//'/'//'fort.99')
C         do i=1,np
C            write(99,'(i5,2x,e16.8)') i, gwce_lv(i)
C         enddo
C         close(99)
C      endif

C...  Save the elevation at the past time step into Eta1 and zero out Eta2
ckmd	Need to save z(s-1) and etas(s-1) for the corrector loop
 
      DO I=1,NP
         IF(CPRECOR) THEN 
            ETAS0(I)=ETAS(I)
            Eta0(I)=Eta1(I)
         END IF 
         Eta1(I)=Eta2(I)
         Eta2(I)=0.0d0
      END DO

C...  At elevation boundary condition nodes, determine the elevation at
C...  the s+1 time step
C...  
C...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0.
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

C...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

      IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TimeLoc-ETIME1)/ETIMINC
         DO I=1,NETA
            NBDI=NBD(I)
            Eta2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
         END DO
      ENDIF
C
C     jgf46.02 Added the ability to include geoid offset on the boundary.
      IF (LoadGeoidOffset) THEN
         DO I=1,NETA
            ETA2(NBD(I))=ETA2(NBD(I))+GeoidOffset(NBD(I))
         END DO
      ENDIF
C
C     jgf48.04 Added an inverted barometer boundary condition so that
C     low pressure systems can cross the boundary without creating an
C     elevation anomaly.
C      DO I=1,NETA
C         ETA2(NBD(I))=ETA2(NBD(I))
C     &       + RampMete*(101300.d0/(RHOWAT0*G)) - PR2(NBD(I))
C      END DO

C   kmd48.33bc add information for the levels of no motion boundary conditions
C              these are considered the steric adjustments. 
      IF ((ABS(RES_BC_FLAG).GE.1).AND.(CBaroclinic))  THEN
         DO I=1,NETA
           NBDI=NBD(I)
           ETA2(NBDI) = ETA2(NBDI) + LNM_BC(I)
         END DO
      END IF

C...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
C...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

C...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
C...  assume that QN is positive into the domain.  QFORCEJ has a -1
C...  built in and the terms are not explicitly negated. Boundary
C...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
C...  flux that  is positive out of the domain.  Therefore, to match
C...  the formulation these terms must be explicitly multiplied by -1.

C...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN 
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 + 
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN 
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

C     jgf46.21 Added IBTYPE=52.
         IF(LBCODEI(1).EQ.52) THEN
            QFORCEJ=(QN2(1)-QN0(1))/DT2 + Tau0VAR(NBDJ)*QN1(1)
            IF (IT.GT.FluxSettlingIT) THEN
               HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
               Celerity=SQRT(G*HH1)
               QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &              + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(1)))
            ENDIF
         ENDIF
C
         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)
            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN 
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0


            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)
C
C     jgf46.21 Added IBTYPE=52
            IF(LBCODEI(J).EQ.52) THEN 
               QFORCEJ=(QN2(J)-QN0(J))/DT2 + Tau0VAR(NBDJ)*QN1(J)
               IF (IT.GT.FluxSettlingIT) THEN
                  HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
                  Celerity=SQRT(G*HH1)
                  QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &                 + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(J)))
               ENDIF
            ENDIF

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.d0*QForceJ+QForceI)
          ENDDO
        ENDIF

C...
C...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
C...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
C...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
C...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
C...  EQUATIONS
C...  
      IF ( ILump.eq.0 ) THEN ! default, fully consistent GWCE LHS
         DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
            DO J=2,NNEIGH(NBDI)
               GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &              -ETAS(NBDI)*OBCCOEF(I,J-1)
            END DO
         END DO
      ELSE                   ! ILump.eq.1, lumped GWCE
         DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*COEFD(NBDI)
         END DO
      ENDIF
C...  
C...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
C...  

C...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
C...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
         CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
      ELSE  ! lumped LHS
         CALL UPDATER(GWCE_LV,COEFD,DUMY1,2)
      ENDIF 
#endif
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
C...  JCG ITERATIVE MATRIX SOLVER
         IPARM(1)=ITMAX
         CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &        IWKSP,NW,WKSP,IPARM,RPARM,IER)

         NUMITR=IPARM(1)
         DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ELSE ! lumped LHS
         DO I = 1, NP
            IF (COEFD(I).eq.0.0d0) THEN
               RDIAG = 0.0d0
            ELSE
               RDIAG = 1.0d0 / COEFD(I)
            ENDIF
            ETAS(I) = GWCE_LV(I) * RDIAG
         ENDDO
         NUMITR=0
         DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ENDIF

C     UPDATE ELEVATIONS

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
      RETURN
      END SUBROUTINE GWCE_NEW



C*******************************************************************************
C                                                                              *
C    Subroutine to compute the velocity and flux/unit width using a 2DDI non-  *
C    conservative momentum equation formluation as contained in original       *
C    versions of ADCIRC (e.g., pre v45.XX).  This formulation is correct if    *
C    the element size is spatially constant but does not account properly for  *
C    spatially variable elemental areas.                                       *
C    Despite its incorrect formulation, we have maintained this subroutine for *
C    posterity and comparison with earlier runs.                               *
C                                                                              *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE Mom_Eqs_Original()

C     
      USE GLOBAL
      USE WIND
      USE NodalAttributes, ONLY: EVM
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
      INTEGER	NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI
      INTEGER NBDI
      INTEGER NNFirst

      REAL(SZ) DDU
      REAL(SZ) DXXYY11, DXXYY12, DXXYY13, DXXYY21
      REAL(SZ) DXXYY22, DXXYY23, DXXYY31
      REAL(SZ) DXXYY32, DXXYY33, DXYH11, DXYH12
      REAL(SZ) EVMN1, EVMN2, EVMN3
      REAL(SZ) EVMAvgDT
      REAL(SZ) FIIN
      REAL(SZ) H1, H1N1, H1N2, H1N3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) QTan
      REAL(SZ) SFacAvg
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) U1N1, U1N2, U1N3
      REAL(SZ) U1AvgDT, U1AvgDTDDX1, U1AvgDTDDX2, U1AvgDTDDX3
      REAL(SZ) V1N1, V1N2, V1N3
      REAL(SZ) V1AvgDT, V1AvgDTDDY1, V1AvgDTDDY2, V1AvgDTDDY3
      REAL(SZ) VCoef1, VCoef2
      REAL(SZ) VCoef3N1, VCoef3N2, VCoef3N3, VCoef3X, VCoef3Y
      REAL(SZ) VelNorm, VelTan
      REAL(SZ) VIDBCPDX, VIDBCPDY                          
      REAL(SZ) WSX, WSY
      REAL(SZ) ZNGLHS,ZNGRHS1,ZNGRHS2

      REAL(8) AREAIE, AREAIE2
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) DDX1,DDX2,DDX3,DDY1,DDY2,DDY3
      REAL(8) DXX11,DXX12,DXX13,DXX21,DXX22,DXX23,DXX31,DXX32,DXX33
      REAL(8) DYY11,DYY12,DYY13,DYY21,DYY22,DYY23,DYY31,DYY32,DYY33
      REAL(8) DXY11,DXY12,DXY13,DXY21,DXY22,DXY23,DXY31,DXY32,DXY33

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

C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF 
C...        THE TIME STEPPING LOOP.

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT

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
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         H1N1=DP(NM1)+IFNLFA*ETA1(NM1)
         H1N2=DP(NM2)+IFNLFA*ETA1(NM2)
         H1N3=DP(NM3)+IFNLFA*ETA1(NM3)
         H2N1=DP(NM1)+IFNLFA*ETA2(NM1)
         H2N2=DP(NM2)+IFNLFA*ETA2(NM2)
         H2N3=DP(NM3)+IFNLFA*ETA2(NM3)
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         AREAIE2=AREAS(IE)
         AREAIE =AREAIE2/2.d0
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
         FDY1=X(NM3)-X(NM2)  !a1
         FDY2=X(NM1)-X(NM3)  !a2
         FDY3=X(NM2)-X(NM1)  !a3
         FDX1O2A=FDX1/AREAIE2  !dphi1/dx
         FDY1O2A=FDY1/AREAIE2  !dphi1/dy
         FDX2O2A=FDX2/AREAIE2  !dphi2/dx
         FDY2O2A=FDY2/AREAIE2  !dphi2/dy
         FDX3O2A=FDX3/AREAIE2  !dphi3/dx
         FDY3O2A=FDY3/AREAIE2  !dphi3/dy
            
         DDX1=FDX1/3.d0      !<2*(dphi1/dx)*phij> j=1,2,3
         DDY1=FDY1/3.d0      !<2*(dphi1/dy)*phij> j=1,2,3
         DXX11=FDX1O2A*FDX1   !<2*(dphi1/dx)*(dphi1/dx)>
         DYY11=FDY1O2A*FDY1   !<2*(dphi1/dy)*(dphi1/dy)>
         DXXYY11=DXX11+DYY11
         DXX12=FDX1O2A*FDX2   !<2*(dphi1/dx)*(dphi2/dx)>
         DYY12=FDY1O2A*FDY2   !<2*(dphi1/dy)*(dphi2/dy)>
         DXXYY12=DXX12+DYY12
         DXX13=FDX1O2A*FDX3   !<2*(dphi1/dx)*(dphi3/dx)>
         DYY13=FDY1O2A*FDY3   !<2*(dphi1/dy)*(dphi3/dy)>
         DXXYY13=DXX13+DYY13

         DDX2=FDX2/3.d0      !<2*(dphi2/dx)*phij> j=1,2,3
         DDY2=FDY2/3.d0      !<2*(dphi2/dy)*phij> j=1,2,3
         DXXYY21=DXXYY12
         DXX22=FDX2O2A*FDX2   !<2*(dphi2/dx)*(dphi2/dx)>
         DYY22=FDY2O2A*FDY2   !<2*(dphi2/dy)*(dphi2/dy)>
         DXXYY22=DXX22+DYY22
         DXX23=FDX2O2A*FDX3   !<2*(dphi2/dx)*(dphi3/dx)>
         DYY23=FDY2O2A*FDY3   !<2*(dphi2/dy)*(dphi3/dy)>
         DXXYY23=DXX23+DYY23

         DDX3=FDX3/3.d0      !<2*(dphi3/dx)*phij> j=1,2,3
         DDY3=FDY3/3.d0      !<2*(dphi3/dy)*phij> j=1,2,3
         DXXYY31=DXXYY13
         DXXYY32=DXXYY23
         DXX33=FDX3O2A*FDX3   !<2*(dphi3/dx)*(dphi3/dx)>
         DYY33=FDY3O2A*FDY3   !<2*(dphi3/dy)*(dphi3/dy)>
         DXXYY33=DXX33+DYY33

         FIIN=AREAIE2/3.D0    !2*<phi*phj> lumped

C...  Accumulate nodal values for terms associated with the barotropic pressure gradient

         VCOEF3N1=ETA1(NM1)+ETA2(NM1)
         VCOEF3N2=ETA1(NM2)+ETA2(NM2)
         VCOEF3N3=ETA1(NM3)+ETA2(NM3)

C......If using atm pressure

         IF(NWS.NE.0) THEN
            VCOEF3N1=VCOEF3N1+PR1(NM1)+PR2(NM1)
            VCOEF3N2=VCOEF3N2+PR1(NM2)+PR2(NM2)
            VCOEF3N3=VCOEF3N3+PR1(NM3)+PR2(NM3)
            ENDIF

C......If using tidal potential terms

         if (CTIP) then
            VCOEF3N1=VCOEF3N1-TIP1(NM1)-TIP2(NM1)
            VCOEF3N2=VCOEF3N2-TIP1(NM2)-TIP2(NM2)
            VCOEF3N3=VCOEF3N3-TIP1(NM3)-TIP2(NM3)
            endif

         VCOEF3N1=VCOEF3N1*GDTO2
         VCOEF3N2=VCOEF3N2*GDTO2
         VCOEF3N3=VCOEF3N3*GDTO2

C...  COMPUTE ELEMENT AVERAGES QUANTITIES

         U1AvgDT=SADVDTO3*(U1N1+U1N2+U1N3)
         V1AvgDT=SADVDTO3*(V1N1+V1N2+V1N3)
         U1AvgDTDDX1=U1AvgDT*DDX1
         U1AvgDTDDX2=U1AvgDT*DDX2
         U1AvgDTDDX3=U1AvgDT*DDX3
         V1AvgDTDDY1=V1AvgDT*DDY1
         V1AvgDTDDY2=V1AvgDT*DDY2
         V1AvgDTDDY3=V1AvgDT*DDY3
         EVMAvgDT=((EVMN1+EVMN2+EVMN3)/3.d0)*DT

C...  ASSEMBLE PARTIAL PRODUCTS

         VCOEF3X=VCOEF3N1*DDX1+VCOEF3N2*DDX2+VCOEF3N3*DDX3
         VCOEF3Y=VCOEF3N1*DDY1+VCOEF3N2*DDY2+VCOEF3N3*DDY3
         ADVECX=(U1AvgDTDDX1+V1AvgDTDDY1)*U1N1
     &           +(U1AvgDTDDX2+V1AvgDTDDY2)*U1N2
     &           +(U1AvgDTDDX3+V1AvgDTDDY3)*U1N3
         ADVECY=(U1AvgDTDDX1+V1AvgDTDDY1)*V1N1
     &           +(U1AvgDTDDX2+V1AvgDTDDY2)*V1N2
     &           +(U1AvgDTDDX3+V1AvgDTDDY3)*V1N3
C
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
C...  TEMP_LV_A
C.... VECTOR FOR NODE NM1

         TEMP_LV_A1=NCELE*(               
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3X
C...  LATERAL VISCOUS TERMS
     &           -EVMAvgDT*(DXXYY11*U1N1+DXXYY12*U1N2+DXXYY13*U1N3)
C...  ADVECTIVE TERMS
     &           -ADVECX
C...  COMMON DIVISION OPERATION
     &           )/FIIN
C...  
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO TEMP_LV_A
C.... VECTOR FOR NODE NM2
C...  
         TEMP_LV_A2=NCELE*(   
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3X
C...  LATERAL VISCOUS TERMS]
     &           -EVMAvgDT*(DXXYY21*U1N1+DXXYY22*U1N2+DXXYY23*U1N3)
C...  ADVECTIVE TERMS
     &           -ADVECX
C...  COMMON DIVISION OPERATION
     &           )/FIIN
C...  
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO TEMP_LV_A
C.... VECTOR FOR NODE NM3
C...  
         TEMP_LV_A3=NCELE*(   
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3X
C...  LATERAL VISCOUS TERMS
     &           -EVMAvgDT*(DXXYY31*U1N1+DXXYY32*U1N2+DXXYY33*U1N3)
C...  ADVECTIVE TERMS
     &           -ADVECX
C...  COMMON DIVISION OPERATION
     &           )/FIIN
C...  
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
C.... VECTOR FOR NODE NM1
C...  
         TEMP_LV_B1=NCELE*( 
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3Y
C...  LATERAL VISCOUS TERMS
     &           -EVMAvgDT*(DXXYY11*V1N1+DXXYY12*V1N2+DXXYY13*V1N3)
C...  ADVECTIVE TERMS
     &           -ADVECY
C...  COMMON DIVISION OPERATION
     &           )/FIIN
C...  
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
C.... VECTOR FOR NODE NM2
C...  
         TEMP_LV_B2=NCELE*(  
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3Y
C...  LATERAL VISCOUS TERMS
     &           -EVMAvgDT*(DXXYY21*V1N1+DXXYY22*V1N2+DXXYY23*V1N3)
C...  ADVECTIVE TERMS
     &           -ADVECY
C...  COMMON DIVISION OPERATION
     &           )/FIIN
C...  
C...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
C.... VECTOR FOR NODE NM3
C...  
         TEMP_LV_B3=NCELE*(   
C...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
     &           -VCOEF3Y
C...  LATERAL VISCOUS TERMS
     &           -EVMAvgDT*(DXXYY31*V1N1+DXXYY32*V1N2+DXXYY33*V1N3)
C...  ADVECTIVE TERMS
     &           -ADVECY
C...  COMMON DIVISION OPERATION
     &           )/FIIN

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

         ENDDO                        !end of elemental assembly loop

  
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
         AreaIE=Areas(IE)/2.d0
         NCEle=NC1*NC2*NC3*NOFF(IE)
         MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A(IE,1)
         MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A(IE,2)
         MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A(IE,3)
         MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B(IE,1)
         MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B(IE,2)
         MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B(IE,3)
         END DO
#endif


C...  UPDATE MOMENTUM EQUATION LHS COEFFICIENTS AND LOAD VECTORS AT EACH
C...  NODE BY DIVIDING BY NUMBER OF ELEMENTS THE NODE IS ASSOCIATED WITH
C...  AND ADDING IN LUMPED TERMS AND BOTTOM FRICTION AND TAKING ACCOUNT
C...  OF THE BOUNDARY CONDITION

      WSX=0.D0
      WSY=0.D0
      VIDBCPDX=0.D0
      VIDBCPDY=0.D0
      DO I=1,NP
         NCI=NODECODE(I)
         MOM_LV_X(I)=MOM_LV_X(I)/MJU(I)
         MOM_LV_Y(I)=MOM_LV_Y(I)/MJU(I)
         H1=DP(I)+IFNLFA*ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
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
            VELNORM=-QN2(I)/H2 
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VELNORM*AUVXY(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUVXX(NBDI)*NCI !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVYX(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VELNORM*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUV11(NBDI)*NCI !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI)
     &                     -CSII(I)*MOM_LV_Y(NBDI)
     &                     -VELNORM*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUV11(NBDI)*NCI !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
            ENDIF

C     Specified essential normal flow and no tangential slip

         IF((LBCODEI(I).GE.10).AND.(LBCODEI(I).LE.19)) THEN
            VELNORM=-QN2(I)/H2
            VELTAN=0.D0
            MOM_LV_X(NBDI)=VELTAN*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*NCI !Normal Eqn RHS
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
               NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
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
C...  SOLVE FOR VELOCITY AT NEW LEVEL  (K+1)
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
         AUVYY(I) = AUVXX(I)
         AUVYX(I) = - AUVXY(I)
         DDU = AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
         UU2(I) = (MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
         VV2(I) = (MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
cjjw
c        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
c           UU2(I)=0.D0    !no slip and normal flux
c           VV2(I)=0.D0    !on wet/dry interface nodes
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

      RETURN
      END


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
      USE NodalAttributes, ONLY: EVM
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
      INTEGER	NM1, NM2, NM3
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

C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF 
C...        THE TIME STEPPING LOOP.

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT

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
         IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
            EVMSmag=EVMEle*
     &             sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA)
     &                 +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
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

      RETURN
      END


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
      USE WIND
      USE NodalAttributes, ONLY: EVM
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
      INTEGER	NM1, NM2, NM3
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

C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF 
C...        THE TIME STEPPING LOOP.

C.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

      DO IE=1,NE

C...  SET NODAL VALUES FOR EACH ELEMENT

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
         IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient

            EVMSmag=0.5d0*EVMEle*
     &             sqrt((DU1DX2A-DV1DY2A)*(DU1DX2A-DV1DY2A)
     &                 +(DU1DY2A+DV1DX2A)*(DU1DY2A+DV1DX2A))
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
         
      RETURN
      END

C*******************************************************************************
C                                                                              *
C    Subroutine to compute the elevation using the GWCE formluation            *
C    This subroutine is the corrector step for the predictor-corrector         *
C    algorithm and obtains the corrected elevations                            *
C    Re-written to conform to the ADCIRC Theory Report                         *
C                                                                              *
C                            k.d.  06/24/2004                                  *
C                            r.l.  06/22/2005                                  *
C*******************************************************************************

      SUBROUTINE GWCE_new_pc(IT,TimeLoc,TimeH)
C     
      USE GLOBAL
      USE SIZES
      USE WIND
      USE ITPACKV
      USE NodalAttributes, ONLY : FRIC, Tau0Var, HBREAK, FTHETA, FGAMMA,
     &     IFLINBF, IFNLBF, IFHYBF, EVM
#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE 

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER	NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      REAL(SZ) A00pB00
      REAL(SZ) BCXAvg, BCYAvg
      REAL(SZ) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(SZ) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(SZ) CorifAvg
      REAL(SZ) DPAvg, GDPAvgOAreaIE4
      REAL(SZ) DispX, DispY, DispXAvg, DispYAvg
      REAL(SZ) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(SZ) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(SZ) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(SZ) ESN1, ESN2, ESN3, ESAvg
      REAL(SZ) EVMH, EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(SZ) EVMEle, EVMSmag
      REAL(SZ) GA00DPAvgOAreaIE4
      REAL(SZ) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(SZ) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg, H1, H2
      REAL(SZ) H2OTotalArea
      REAL(SZ) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(SZ) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(SZ) MsFacR
      REAL(SZ) MX, MY, MXAvg, MYAvg
      REAL(SZ) JXAvg, JYAvg
      REAL(SZ) Pr1N1, Pr1N2, Pr1N3
      REAL(SZ) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(SZ) SFacAvg
      REAL(SZ) T0N1,T0N2, T0N3
      REAL(SZ) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(SZ) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TiPN1, TiPN2, TiPN3
      REAL(SZ) UV0, UV1, UV2
      REAL(SZ) U1N1,U1N2,U1N3, U1Avg
      REAL(SZ) V1N1,V1N2,V1N3, V1Avg
      REAL(SZ) WSXAvg, WSYAvg

      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) TimeLoc, TimeH

ckmd	 Added in parameters for the pc algorithm
      REAL(SZ) BSX0N1, BSX0N2, BSX0N3, BSY0N1
      REAL(SZ) BSY0N2, BSY0N3, BSX0Avg, BSY0Avg
      REAL(SZ) BSX2N1, BSX2N2, BSX2N3, BSY2N1
      REAL(SZ) BSY2N2, BSY2N3, BSX2Avg, BSY2Avg
      REAL(SZ) E2N1,E2N2,E2N3
      REAL(SZ) E0N1SQ, E0N2SQ, E0N3SQ
      REAL(SZ) E2N1SQ, E2N2SQ, E2N3SQ
      REAL(SZ) H0N1, H0N2, H0N3, H00
      REAL(SZ) QX0N1, QX0N2, QX0N3, QY0N1, QY0N2, QY0N3, QX0Avg, QY0Avg
      REAL(SZ) QX2N1, QX2N2, QX2N3, QY2N1, QY2N2, QY2N3, QX2Avg, QY2Avg
      REAL(SZ) Tau0QX0Avg, Tau0QY0Avg, Tau0QX2Avg, Tau0QY2Avg
      REAL(SZ) Tau0SpaVar0, Tau0SpaVar2
      REAL(SZ) U0N1,U0N2,U0N3, U0Avg
      REAL(SZ) V0N1,V0N2,V0N3, V0Avg
      REAL(SZ) U2N1,U2N2,U2N3, U2Avg
      REAL(SZ) V2N1,V2N2,V2N3, V2Avg
      REAL(SZ) timewtgwce0,timewtgwce1,timewtgwce2,timeagflag


ckmd    Must reset the result vector to zero before recomputing
ckmd	the next time level.
      DO I=1,NP
         GWCE_LV(I) =0.D0
      END DO

C     Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1 
C     Re-compute these local values

      OnDiag=(1+ILump)*2                         !diagonal coefficient
      OffDiag=(1-ILump)                          !off diagonal coefficient

c...  
C...  Compute the GWCE load vector GWCE_LV
C...  This is done primarily element by element by forming 
C...  temporary vectors and then assembling at the end.  
C...  This has been set up to unroll loops to optimize performance
C...  on vector processors.  
C...
C...  Elevation and flux boundary conditions are imposed after the
C...  element by element assembly section.  
C...

C...  Initialize variables to zero if these forcings are not used
 
      IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
      ELSE
        WSXAvg=0.d0
        WSYAvg=0.d0
        Pr1N1=0.d0
        Pr1N2=0.d0
        Pr1N3=0.d0
      ENDIF

      IF (CTIP) THEN
      ELSE
        TiPN1=0.d0
        TiPN2=0.d0
        TiPN3=0.d0
      ENDIF

      IF(C3D) THEN
      ELSE
        DispXAvg=0.d0
        DispYAvg=0.d0
      ENDIF

      IF(CBaroclinic) THEN
	ELSE
        BCXAvg=0.d0
        BCYAvg=0.d0
      ENDIF

C...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartV).OR.(CGWCE_LS_2PartSV)) THEN

         DO I=1,NP
            LSXX(I)=0.d0
            LSXY(I)=0.d0
            LSYX(I)=0.d0
            LSYY(I)=0.d0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0   !A*DUDX
            LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0   !A*DUDY
            LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0   !A*DVDX
            LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0   !A*DVDY
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               EVMEle=EVMSmag
            ENDIF
            LSXXEle = LSXXGradA*EVMEle            
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
               IF (CGWCE_LS_2PartV) THEN          !nonsymmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=H2OTotalArea*LSXY(I)
                  LSYX(I)=H2OTotalArea*LSYX(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
               IF (CGWCE_LS_2PartSV) THEN         !symmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=0.5d0*H2OTotalArea*(LSXY(I)+LSYX(I))
                  LSYX(I)=LSXY(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
            ELSE
               LSXX(I)=0.d0
               LSXY(I)=0.d0
               LSYX(I)=0.d0
               LSYY(I)=0.d0
            ENDIF
         ENDDO

      ENDIF

C...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartSQ)) THEN

         DO I=1,NP
            LSXX(I)=0.d0
            LSXY(I)=0.d0
            LSYX(I)=0.d0
            LSYY(I)=0.d0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
              LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
              LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
              LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
              LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0
              EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
              EVMEle=EVMSmag
            ENDIF
            LSXXGradA=(QX1(NM1)*FDX1+QX1(NM2)*FDX2+QX1(NM3)*FDX3)/2.d0
            LSXYGradA=(QX1(NM1)*FDY1+QX1(NM2)*FDY2+QX1(NM3)*FDY3)/2.d0
            LSYXGradA=(QY1(NM1)*FDX1+QY1(NM2)*FDX2+QY1(NM3)*FDX3)/2.d0
            LSYYGradA=(QY1(NM1)*FDY1+QY1(NM2)*FDY2+QY1(NM3)*FDY3)/2.d0
            LSXXEle = LSXXGradA*EVMEle            
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               IF (CGWCE_LS_2PartQ) THEN          !nonsymmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=LSXY(I)/TotalArea(I)
                  LSYX(I)=LSYX(I)/TotalArea(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
               IF (CGWCE_LS_2PartSQ) THEN         !symmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=0.5d0*(LSXY(I)+LSYX(I))/TotalArea(I)
                  LSYX(I)=LSXY(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
            ELSE
               LSXX(I)=0.d0
               LSXY(I)=0.d0
               LSYX(I)=0.d0
               LSYY(I)=0.d0
            ENDIF
         ENDDO

      ENDIF

      DO I=1,NP
ckmd
ckmd	Added in the three time levels for the tau term.
ckmd	Every term is updated for the three time levels.
ckmd
            UV0=SQRT(UU0(I)*UU0(I)+VV0(I)*VV0(I))
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            UV2=SQRT(UU2(I)*UU2(I)+VV2(I)*VV2(I))
            H00=DP(I)+IFNLFA*ETA0(I)
            H1=DP(I)+IFNLFA*ETA1(I)
            H2=DP(I)+IFNLFA*ETA2(I)
            TK0(I)=FRIC(I)*(IFLINBF + (UV0/H00)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H00)**FTHETA)**(FGAMMA/FTHETA)))
            TK(I)=FRIC(I)*(IFLINBF + (UV1/H1)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
            TK2(I)=FRIC(I)*(IFLINBF + (UV2/H2)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H2)**FTHETA)**(FGAMMA/FTHETA)))
       END DO

ckmd      Added in the time weights
C...     Time weights for the nonlinear terms in the GWCE for 
C...       the corrector step
          timewtgwce0=0.33d0
          timewtgwce1=0.34d0
          timewtgwce2=0.33d0 
          timeagflag=1.0d0

C...  Assemble the GWCE RHS except for the boundary integral terms
ckmd	Renumber the GWCE loop for the corrector step

      DO 1038 IE=1,NE
  
C...     Set nodal values for each element
ckmd
ckmd	Define the needed product terms at three time levels
ckmd

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         E0N1=ETA0(NM1)
         E0N2=ETA0(NM2)
         E0N3=ETA0(NM3)
         E1N1=ETA1(NM1)
         E1N2=ETA1(NM2)
         E1N3=ETA1(NM3)
         E2N1=ETA2(NM1)
         E2N2=ETA2(NM2)
         E2N3=ETA2(NM3)
         E0N1SQ=E0N1*E0N1
         E0N2SQ=E0N2*E0N2
         E0N3SQ=E0N3*E0N3
         E1N1SQ=E1N1*E1N1
         E1N2SQ=E1N2*E1N2
         E1N3SQ=E1N3*E1N3
         E2N1SQ=E2N1*E2N1
         E2N2SQ=E2N2*E2N2
         E2N3SQ=E2N3*E2N3
         ESN1=ETAS0(NM1)
         ESN2=ETAS0(NM2)
         ESN3=ETAS0(NM3)
         U0N1=UU0(NM1)
         U0N2=UU0(NM2)
         U0N3=UU0(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         U2N1=UU2(NM1)
         U2N2=UU2(NM2)
         U2N3=UU2(NM3)
         V0N1=VV0(NM1)
         V0N2=VV0(NM2)
         V0N3=VV0(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         V2N1=VV2(NM1)
         V2N2=VV2(NM2)
         V2N3=VV2(NM3)
         QX0N1=QX0(NM1)
         QX0N2=QX0(NM2)
         QX0N3=QX0(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QX2N1=QX2(NM1)
         QX2N2=QX2(NM2)
         QX2N3=QX2(NM3)
         QY0N1=QY0(NM1)
         QY0N2=QY0(NM2)
         QY0N3=QY0(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         QY2N1=QY2(NM1)
         QY2N2=QY2(NM2)
         QY2N3=QY2(NM3)
         H0N1=DP(NM1)+IFNLFA*E0N1
         H0N2=DP(NM2)+IFNLFA*E0N2
         H0N3=DP(NM3)+IFNLFA*E0N3
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         H2N1=DP(NM1)+IFNLFA*E2N1
         H2N2=DP(NM2)+IFNLFA*E2N2
         H2N3=DP(NM3)+IFNLFA*E2N3
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         T0N1=Tau0Var(NM1) 
         T0N2=Tau0Var(NM2) 
         T0N3=Tau0Var(NM3) 
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            Pr1N1=PR1(NM1)
            Pr1N2=PR1(NM2)
            Pr1N3=PR1(NM3)
         ENDIF
         IF (CTIP) THEN                        !tidal potential
            TiPN1=TiP1(NM1)
            TiPN2=TiP1(NM2)
            TiPN3=TiP1(NM3)
         ENDIF
         IF (C2DDI) THEN                       !2D bottom friction
            BSX0N1=TK0(NM1)*QX0N1
            BSY0N1=TK0(NM1)*QY0N1
            BSX0N2=TK0(NM2)*QX0N2
            BSY0N2=TK0(NM2)*QY0N2
            BSX0N3=TK0(NM3)*QX0N3
            BSY0N3=TK0(NM3)*QY0N3
            BSXN1=TK(NM1)*QX1N1
            BSYN1=TK(NM1)*QY1N1
            BSXN2=TK(NM2)*QX1N2
            BSYN2=TK(NM2)*QY1N2
            BSXN3=TK(NM3)*QX1N3
            BSYN3=TK(NM3)*QY1N3
            BSX2N1=TK2(NM1)*QX2N1
            BSY2N1=TK2(NM1)*QY2N1
            BSX2N2=TK2(NM2)*QX2N2
            BSY2N2=TK2(NM2)*QY2N2
            BSX2N3=TK2(NM3)*QX2N3
            BSY2N3=TK2(NM3)*QY2N3

         ENDIF      
         IF (C3D) THEN                         !3D bottom friction
            BSXN1=BSX1(NM1)
            BSXN2=BSX1(NM2)
            BSXN3=BSX1(NM3)
            BSYN1=BSY1(NM1)
            BSYN2=BSY1(NM2)
            BSYN3=BSY1(NM3)
         ENDIF

         AreaIE2=Areas(IE)               !2A
         AreaIE=AreaIE2/2.d0             ! A
         AreaIE4=2.d0*AreaIE2            !4A

         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0

         FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1 = 2*Area*dphi1/dx
         FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2 = 2*Area*dphi2/dx
         FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3 = 2*Area*dphi3/dx
         FDY1 = X(NM3)-X(NM2)           !a1 = 2*Area*dphi1/dy
         FDY2 = X(NM1)-X(NM3)           !a2 = 2*Area*dphi2/dy
         FDY3 = X(NM2)-X(NM1)           !a3 = 2*Area*dphi3/dy

C...     Compute part of several spatial gradients for use below

         E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
         E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
         E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
         E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
         Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
         Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

C...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM

         IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.d0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &         +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)
     &         -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)
     &         +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)
     &         -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
         ENDIF

C...     Compute the remainder of the 2 Part lateral stress terms

         IF((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
         ENDIF

C...     Compute the spatial gradients of the velocity dispersion terms if 3D

         IF (C3D) THEN                         !3D bottom friction
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
         ENDIF

C...     Compute elemental averages                   

         CorifAvg=(Corif(NM1)+Corif(NM2)+Corif(NM3))/3.d0
         Tau0Avg=(T0N1+T0N2+T0N3)/3.d0
         Tau0QX0Avg=(T0N1*QX0N1+T0N2*QX0N2+T0N3*QX0N3)/3.d0
         Tau0QY0Avg=(T0N1*QY0N1+T0N2*QY0N2+T0N3*QY0N3)/3.d0
         Tau0QXAvg=(T0N1*QX1N1+T0N2*QX1N2+T0N3*QX1N3)/3.d0
         Tau0QYAvg=(T0N1*QY1N1+T0N2*QY1N2+T0N3*QY1N3)/3.d0
         Tau0QX2Avg=(T0N1*QX2N1+T0N2*QX2N2+T0N3*QX2N3)/3.d0
         Tau0QY2Avg=(T0N1*QY2N1+T0N2*QY2N2+T0N3*QY2N3)/3.d0
         U0Avg=(U0N1+U0N2+U0N3)/3.d0
         V0Avg=(V0N1+V0N2+V0N3)/3.d0
         U1Avg=(U1N1+U1N2+U1N3)/3.d0
         V1Avg=(V1N1+V1N2+V1N3)/3.d0
         U2Avg=(U2N1+U2N2+U2N3)/3.d0
         V2Avg=(V2N1+V2N2+V2N3)/3.d0
         QX0Avg=(QX0N1+QX0N2+QX0N3)/3.d0
         QY0Avg=(QY0N1+QY0N2+QY0N3)/3.d0
         QX1Avg=(QX1N1+QX1N2+QX1N3)/3.d0
         QY1Avg=(QY1N1+QY1N2+QY1N3)/3.d0
         QX2Avg=(QX2N1+QX2N2+QX2N3)/3.d0
         QY2Avg=(QY2N1+QY2N2+QY2N3)/3.d0
         ESAvg=(ESN1+ESN2+ESN3)/3.d0
         DPAvg=(DP(NM1)+DP(NM2)+DP(NM3))/3.d0
         GDPAvgOAreaIE4=G*DPAvg/AreaIE4
         HAvg=(H1N1+H1N2+H1N3)/3.d0
         GHAvg=G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSX0Avg=(BSX0N1+BSX0N2+BSX0N3)/3.d0
         BSY0Avg=(BSY0N1+BSY0N2+BSY0N3)/3.d0
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.d0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.d0
         BSX2Avg=(BSX2N1+BSX2N2+BSX2N3)/3.d0
         BSY2Avg=(BSY2N1+BSY2N2+BSY2N3)/3.d0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
         ENDIF
         IF (C3D) THEN                !3D velocity dispersion
            DispXAvg=IFNLCT*DispX
            DispYAvg=IFNLCT*DispY
         ENDIF
         IF(CBaroclinic) THEN
            BCXAvg=(H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                       +H1N3*VIDBCPDXOH(NM3))/3.d0
            BCYAvg=(H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                       +H1N3*VIDBCPDYOH(NM3))/3.d0
         ENDIF

C...     Compute additional partial factors

         MsFacR=AreaIE*(1.d0/DT-Tau0Avg/2.d0)/DT/12.d0
         GOAreaIE4=G/AreaIE4
         Tau0SpaVar0=(QX0Avg*Tau0XGrad2A+QY0Avg*Tau0YGrad2A)/6.d0
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
         Tau0SpaVar2=(QX2Avg*Tau0XGrad2A+QY2Avg*Tau0YGrad2A)/6.d0
         A00pB00=A00+B00

C...     Compute the JX, JY terms less the advection terms
         JXAvg = timewtgwce0*CorifAvg*QY0Avg+timewtgwce1*CorifAvg
     &                          *QY1Avg+timewtgwce2*CorifAvg*QY2Avg
     &          -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDX1+E0N2SQ*FDX2
     &                                        +E0N3SQ*FDX3)
     &          -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDX1+E1N2SQ*FDX2
     &                                        +E1N3SQ*FDX3)
     &          -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDX1+E2N2SQ*FDX2
     &                                        +E2N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1
     &                     +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3)
     &          +WSXAvg-timewtgwce0*BSX0Avg-timewtgwce1*BSXAvg
     &                                        -timewtgwce2*BSX2Avg
     &          +MXAvg-DispXAvg-BCXAvg
     &          +timewtgwce0*Tau0QX0Avg+timewtgwce1*Tau0QXAvg
     &                                 +timewtgwce2*Tau0QX2Avg 
         JYAvg =-timewtgwce0*CorifAvg*QX0Avg-timewtgwce1*CorifAvg
     &                          *QX1Avg-timewtgwce2*CorifAvg*QX2Avg
     &          -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDY1+E0N2SQ*FDY2
     &                                        +E0N3SQ*FDY3)
     &          -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDY1+E1N2SQ*FDY2
     &                                        +E1N3SQ*FDY3)
     &          -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDY1+E2N2SQ*FDY2
     &                                        +E2N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1
     &                    +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3)
     &          +WSYAvg-timewtgwce0*BSY0Avg-timewtgwce1*BSYAvg
     &                                      -timewtgwce2*BSY2Avg
     &          +MYAvg-DispYAvg-BCYAvg
     &          +timewtgwce0*Tau0QY0Avg+timewtgwce1*Tau0QYAvg
     &                                    +timewtgwce2*Tau0QY2Avg

C...     Complete the JX, JY terms depending on the advection formulation
         IF(CGWCE_Advec_NC) THEN        !nonconservative advection
           JXAvg = JXAvg - IFNLCT*timewtgwce0*(
     &              QX0Avg*(U0N1*FDX1+U0N2*FDX2+U0N3*FDX3)
     &             +QY0Avg*(U0N1*FDY1+U0N2*FDY2+U0N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce1
     &             *(QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce2
     &             *(QX2Avg*(U2N1*FDX1+U2N2*FDX2+U2N3*FDX3)
     &             +QY2Avg*(U2N1*FDY1+U2N2*FDY2+U2N3*FDY3))/AreaIE2
     &             +IFNLCAT*(timewtgwce0*U0Avg+timewtgwce1*U1Avg
     &             +timewtgwce2*U2Avg)
     &             *(timeagflag*0.5d0*((E2N1-E0N1+E2N2-E0N2
     &             +E2N3-E0N3)/DT)
     &             +(1.D0-timeagflag)*ESAvg/DT)
           JYAvg = JYAvg - IFNLCT*timewtgwce0*(
     &              QX0Avg*(V0N1*FDX1+V0N2*FDX2+V0N3*FDX3)
     &             +QY0Avg*(V0N1*FDY1+V0N2*FDY2+V0N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce1
     &             *(QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce2
     &             *(QX2Avg*(V2N1*FDX1+V2N2*FDX2+V2N3*FDX3)
     &             +QY2Avg*(V2N1*FDY1+V2N2*FDY2+V2N3*FDY3))/AreaIE2
     &             +IFNLCAT*(timewtgwce0*V0Avg+timewtgwce1*V1Avg
     &             +timewtgwce2*V2Avg)
     &             *(timeagflag*0.5d0*((E2N1-E0N1+E2N2-E0N2
     &             +E2N3-E0N3)/DT)
     &             +(1.D0-timeagflag)*ESAvg/DT)
         ENDIF
         IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
           JXAvg = JXAvg - IFNLCT*(
     &              (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2
     &                              +U1N3*QX1N3*FDX3)
     &             +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2
     &                              +U1N3*QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2
     &                              +V1N3*QX1N3*FDX3)
     &             +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2
     &                              +V1N3*QY1N3*FDY3))/AreaIE2
         ENDIF
         IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3)
     &             +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3)
     &             +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
         ENDIF


C...     Assemble forcing for node NM1 (local index j=1) 

         Temp_LV_A1=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN1 + OffDiag*(ESN2+ESN3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX1*E0XGrad2A+FDY1*E0YGrad2A)
     &                        +A00pB00*(FDX1*E1XGrad2A+FDY1*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX1 + JYAvg*FDY1)/2.d0

C...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                    +timewtgwce2*Tau0SpaVar2


C...     Assemble forcing for node NM2 (local index j=2) 

         Temp_LV_A2=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN2 + OffDiag*(ESN1+ESN3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX2*E0XGrad2A+FDY2*E0YGrad2A)
     &                        +A00pB00*(FDX2*E1XGrad2A+FDY2*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX2 + JYAvg*FDY2)/2.d0

C...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                   +timewtgwce2*Tau0SpaVar2


C...     Assemble forcing for node NM3 (local index j=3) 

         Temp_LV_A3=

C...     Transient and Tau0 terms from LHS
C...    (consistent mass matrix: ILump=0, lumped mass matrix: ILump=1)
     &         (OnDiag*ESN3 + OffDiag*(ESN1+ESN2))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX3*E0XGrad2A+FDY3*E0YGrad2A)
     &                        +A00pB00*(FDX3*E1XGrad2A+FDY3*E1YGrad2A)) 

C...     Terms from momentum eqs.
     &        +(JXAvg*FDX3 + JYAvg*FDY3)/2.d0

C...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                    +timewtgwce2*Tau0SpaVar2


C...     Put these partial products into further elemental storage for a vector computer
C...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
         Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
         Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
         Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

C...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif

c        IF(IE.EQ.1) THEN
c           WRITE(101,*) ' '
c           WRITE(101,*) '  ************* GWCE Load Vector ************'
c           WRITE(101,*) '  Time Step = ',IT
c           ENDIF
c        WRITE(101,*) IE, ESN1, ESN2, ESN3
c        WRITE(101,*) IE, TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3

ckmd  Change the number of the loop
 1038 CONTINUE      !End of elemental loop


C...  Put load vector elemental values into nodal storage for a vector computer
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A(IE,1)
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A(IE,2)
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A(IE,3)
      END DO
#endif


C...  Save the elevation at the past time step into Eta1 and zero out Eta2
ckmd
ckmd	Already did this and haven't advanced in time yet, so just
ckmd	comment out the save elevation line. But we do need to zero
ckmd	out eta2 because it's involved in the summation below.
ckmd
      DO I=1,NP
c         Eta1(I)=Eta2(I)
         Eta2(I)=0.0d0
      END DO

C...  At elevation boundary condition nodes, determine the elevation at
C...  the s+1 time step
C...  
C...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0.
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

C...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

      IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TimeLoc-ETIME1)/ETIMINC
         DO I=1,NETA
            NBDI=NBD(I)
            Eta2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
         END DO
      ENDIF


C...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
C...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

C...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
C...  assume that QN is positive into the domain.  QFORCEJ has a -1
C...  built in and the terms are not explicitly negated. Boundary
C...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
C...  flux that  is positive out of the domain.  Therefore, to match
C...  the formulation these terms must be explicitly multiplied by -1.

C...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN 
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 + 
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN 
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)

            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN 
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.d0*QForceJ+QForceI)
          ENDDO
        ENDIF

C...
C...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
C...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
C...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
C...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
C...  EQUATIONS
C...  
      DO I=1,NETA
         NBDI=NBD(I)
         ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
         GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
         DO J=2,NNEIGH(NBDI)
            GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &           -ETAS(NBDI)*OBCCOEF(I,J-1)
         END DO
      END DO

C...  
C...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
C...  

C...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
C...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
#endif

C...  JCG ITERATIVE MATRIX SOLVER
      IPARM(1)=ITMAX
      CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &     IWKSP,NW,WKSP,IPARM,RPARM,IER)

      NUMITR=IPARM(1)
      DO I=1,NP
         ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
      END DO
     
C     UPDATE ELEVATIONS

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif

ckmd
ckmd	Eta2 values are now corrected elevations at time level s+1.
ckmd

      RETURN
      END

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
      USE WIND
      USE NodalAttributes, ONLY: EVM

      IMPLICIT NONE 
 
      INTEGER IE, I, J, N                           !local loop counters
      INTEGER	NM1, NM2, NM3
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

ckmd	 Added in parameters for the pc algorithm
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

C...
C...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
C...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF 
C...        THE TIME STEPPING LOOP.
ckmd    Must reset the result vectors to zero before recomputing
ckmd	the next time level.
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
         IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
            EVMSmag=EVMEle*
     &             sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA)
     &                 +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
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
ckmd	Added time weights to tau terms in the momentum equation here.
ckmd	Note the weighting should follow the weighting in the GWCE.
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

ckmd	  Change for the corrector formulation
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

      RETURN
      END


C****************************************************************************************
C   Subroutine to compute Scalar Transport                                              *
C                                                                                       *
C   Note, this is not set up for parallel operation.                                    *
C****************************************************************************************
C
      SUBROUTINE SCALAR_TRANS_2D (IT,TimeLoc)

      USE GLOBAL
      USE NodalAttributes, ONLY: EVC
#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE 
 
      INTEGER IE, I                         !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI

      REAL(SZ) C1, C2, C3, CBEDSTRD, CBEDSTRE, CCRITD
      REAL(SZ) CH1N1, CH1N2, CH1N3, CHSUM
      REAL(SZ) DHDX, DHDY
      REAL(SZ) DXXYY11, DXXYY12, DXXYY13
      REAL(SZ) DXXYY22, DXXYY23
      REAL(SZ) DXXYY33
      REAL(SZ) ECONST
      REAL(SZ) EVC1, EVC2, EVC3, EVCEA
      REAL(SZ) FDDDODT, FDDODODT
      REAL(SZ) HEA                
      REAL(SZ) H1, H1N1, H1N2, H1N3
      REAL(SZ) HSD, HSE
      REAL(SZ) SFacAvg
      REAL(SZ) SS1N1, SS1N2, SS1N3
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(SZ) U1N1, U1N2, U1N3
      REAL(SZ) V1N1, V1N2, V1N3
      REAL(SZ) UV1
      REAL(SZ) UEA, VEA, UPEA, VPEA
      REAL(SZ) WS, WSMOD

      REAL(8) AreaIE2
      REAL(8) FDDD, FDDOD    
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) DDX1,DDX2,DDX3,DDY1,DDY2,DDY3
      REAL(8) DXX11,DXX12,DXX13,DXX21,DXX22,DXX23,DXX31,DXX32,DXX33
      REAL(8) DYY11,DYY12,DYY13,DYY21,DYY22,DYY23,DYY31,DYY32,DYY33
      REAL(8) TimeLoc


C...  NOTE: THE VARIABLE CH1(I) IS ACTUALLY C*H
C.... COMPUTE SOURCE/SINK TERM AT THE NODES USING CLASSICAL COHESIVE
C.... SEDIMENT TRANSPORT RELATIONS

      WS = 0.0001d0          ! particle fall velocity [m/s]
      CBEDSTRD = 0.15d0      ! critical shear stress for deposition [N/m^2]
      CCRITD = 0.30d0        ! critical concentration for hindered settling [kg/m^3]
      ECONST = 0.00001d0     ! erosion rate constant [kg/m^2/sec]
      CBEDSTRE = 0.4d0       ! critical shear stress for erosion [N/m^2]

      DO I=1,NP
        UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
        H1=DP(I)+IFNLFA*ETA1(I)
        BEDSTR=H1*UV1*TK(I)*RhoWat0                           ![N/m^2]
        C1=CH1(I)/H1

C.....Calculate the deposition rate using Krone's (1962) formulation:
C.....dC/dt = -P*WSMOD*C/D     where
C.....WSMOD=WS          when C < Ccrit  and
C.....WSMOD=K*C**1.33   when C > Ccrit
C.....D is the average depth through which particles settle D = H/2,
C.....H is the water depth
C.....C is the depth-averaged sediment concentration,
C.....P is the sticking probability  P = (1-BEDSTR/CBEDSTRD),
C.....CBEDSTRD is the critical bottom stress above which no deposition occurs.
C.....It was assumed that the constant K could be backed out by setting
C.....WSMOD = WS when C = Ccrit.

        WSMOD=WS
        IF(C1.GT.CCRITD) WSMOD=WS*(C1/CCRITD)**1.33d0
        HSD=0.d0
        IF(BEDSTR.LT.CBEDSTRD) HSD=-(2.d0*WSMOD*C1)*
     &                                           (1.0d0-BEDSTR/CBEDSTRD)
        IF(HSD.GT.0.d0) HSD=0.d0

C.....Calculate the surface erosion rate for cohesive sediment using
C.....the Ariathurai et at. (1977) adaption of Partheniades' (1962) findings

        HSE=0.
        IF(BEDSTR.GT.CBEDSTRE) HSE=ECONST*(BEDSTR/CBEDSTRE-1.0)

C.....Determine the total source sink term

        SOURSIN(I)=HSD+HSE
      END DO

C.... UPDATE THE TRANSPORT EQUATION ELEMENT BY ELEMENT BY FORMING
C.... TEMPORARY VECTORS AND THEN ASSEMBLING.  NOTE: TRANS_LV_B(I), TRANS_LV_A(I) ARE
C.... ZEROED OUT AT THE TOP OF THE TIME STEPPING LOOP.  AGAIN THESE
C.... LOOPS HAVE BEEN UNROLLED TO OPTIMIZE VECTORIZATION

      DO IE=1,NE

C.....SET NODAL VALUES FOR EACH ELEMENT

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
         CH1N1=CH1(NM1)
         CH1N2=CH1(NM2)
         CH1N3=CH1(NM3)
         EVC1=EVC(NM1)
         EVC2=EVC(NM2)
         EVC3=EVC(NM3)
         SS1N1=SOURSIN(NM1)
         SS1N2=SOURSIN(NM2)
         SS1N3=SOURSIN(NM3)
         H1N1=DP(NM1)+IFNLFA*ETA1(NM1)
         H1N2=DP(NM2)+IFNLFA*ETA1(NM2)
         H1N3=DP(NM3)+IFNLFA*ETA1(NM3)
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.

C.....COMPUTE ELEMENTAL MATRICIES

         AREAIE2=AREAS(IE)    !2*element area
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
         FDY1=X(NM3)-X(NM2)  !a1
         FDY2=X(NM1)-X(NM3)  !a2
         FDY3=X(NM2)-X(NM1)  !a3
         FDX1O2A=FDX1/AREAIE2  !dphi1/dx
         FDY1O2A=FDY1/AREAIE2  !dphi1/dy
         FDX2O2A=FDX2/AREAIE2  !dphi2/dx
         FDY2O2A=FDY2/AREAIE2  !dphi2/dy
         FDX3O2A=FDX3/AREAIE2  !dphi3/dx
         FDY3O2A=FDY3/AREAIE2  !dphi3/dy

         DDX1=FDX1/3.        !<2*(dphi1/dx)*phij> j=1,2,3
         DDY1=FDY1/3.        !<2*(dphi1/dy)*phij> j=1,2,3
         DXX11=FDX1O2A*FDX1   !<2*(dphi1/dx)*(dphi1/dx)>
         DYY11=FDY1O2A*FDY1   !<2*(dphi1/dy)*(dphi1/dy)>
         DXXYY11=DXX11+DYY11
         DXX12=FDX1O2A*FDX2   !<2*(dphi1/dx)*(dphi2/dx)>
         DYY12=FDY1O2A*FDY2   !<2*(dphi1/dy)*(dphi2/dy)>
         DXXYY12=DXX12+DYY12
         DXX13=FDX1O2A*FDX3   !<2*(dphi1/dx)*(dphi3/dx)>
         DYY13=FDY1O2A*FDY3   !<2*(dphi1/dy)*(dphi3/dy)>
         DXXYY13=DXX13+DYY13

         DDX2=FDX2/3.        !<2*(dphi2/dx)*phij> j=1,2,3
         DDY2=FDY2/3.        !<2*(dphi2/dy)*phij> j=1,2,3
         DXX22=FDX2O2A*FDX2   !<2*(dphi2/dx)*(dphi2/dx)>
         DYY22=FDY2O2A*FDY2   !<2*(dphi2/dy)*(dphi2/dy)>
         DXXYY22=DXX22+DYY22
         DXX23=FDX2O2A*FDX3   !<2*(dphi2/dx)*(dphi3/dx)>
         DYY23=FDY2O2A*FDY3   !<2*(dphi2/dy)*(dphi3/dy)>
         DXXYY23=DXX23+DYY23

         DDX3=FDX3/3.        !<2*(dphi3/dx)*phij> j=1,2,3
         DDY3=FDY3/3.        !<2*(dphi3/dy)*phij> j=1,2,3
         DXX33=FDX3O2A*FDX3   !<2*(dphi3/dx)*(dphi3/dx)>
         DYY33=FDY3O2A*FDY3   !<2*(dphi3/dy)*(dphi3/dy)>
         DXXYY33=DXX33+DYY33

         LUMPT=1             !=1/0; LUMP/DO NOT LUMP THE TRANSPORT EQN
         FDDD=(1+LUMPT)*AREAIE2/6.D0 !<2*(phii*phij) i=j>
         FDDOD=(1-LUMPT)*AREAIE2/12.D0 !<2*(phii*phij) i<>j>
         FDDDODT=FDDD/DTDP
         FDDODODT=FDDOD/DTDP

C.....COMPUTE ELEMENTAL QUANTITIES

         UEA=(U1N1+U1N2+U1N3)/3.
         VEA=(V1N1+V1N2+V1N3)/3.
         HEA=(H1N1+H1N2+H1N3)/3.
         EVCEA=(EVC1+EVC2+EVC3)/3.
         DHDX=H1N1*FDX1O2A+H1N2*FDX2O2A+H1N3*FDX3O2A
         DHDY=H1N1*FDY1O2A+H1N2*FDY2O2A+H1N3*FDY3O2A
         UPEA=UEA+DHDX*EVCEA/HEA
         VPEA=VEA+DHDY*EVCEA/HEA

C.....ASSEMBLE PARTIAL PRODUCT

         CHSUM=CH1N1+CH1N2+CH1N3

C.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A1 AND
C.....TEMP_LV_B1 VECTORS FOR NODE NM1

         TEMP_LV_B1=             !LOAD VECTOR
C......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
     &           FDDDODT*CH1N1+FDDODODT*(CH1N2+CH1N3)
C......LATERAL SGS TERMS
     &           -EVCEA*(DXXYY11*CH1N1+DXXYY12*CH1N2+DXXYY13*CH1N3)
C......ADVECTIVE TERMS
     &           +(UPEA*DDX1+VPEA*DDY1)*CHSUM
C......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
     &           +FDDD*SS1N1+FDDOD*(SS1N2+SS1N3)
         TEMP_LV_A1=             !LHS VECTOR
C......TRANSIENT TERM (LUMPED)
     &           FDDDODT+2.*FDDODODT

C.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A2 AND
C.....TEMP_LV_B2 VECTOR FOR NODE NM2

         TEMP_LV_B2=             !LOAD VECTOR
C......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
     &           FDDDODT*CH1N2+FDDODODT*(CH1N1+CH1N3)
C......LATERAL SGS TERMS
     &           -EVCEA*(DXXYY12*CH1N1+DXXYY22*CH1N2+DXXYY23*CH1N3)
C......ADVECTIVE TERMS
     &           +(UPEA*DDX2+VPEA*DDY2)*CHSUM
C......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
     &           +FDDD*SS1N2+FDDOD*(SS1N1+SS1N3)
         TEMP_LV_A2=             !LHS VECTOR
C......TRANSIENT TERM (LUMPED)
     &           FDDDODT+2.*FDDODODT

C.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A3 AND
C.....TEMP_LV_B3 VECTOR FOR NODE NM3

         TEMP_LV_B3=             !LOAD VECTOR
C......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
     &           FDDDODT*CH1N3+FDDODODT*(CH1N1+CH1N2)
C......LATERAL SGS TERMS
     &           -EVCEA*(DXXYY13*CH1N1+DXXYY23*CH1N2+DXXYY33*CH1N3)
C......ADVECTIVE TERMS
     &           +(UPEA*DDX3+VPEA*DDY3)*CHSUM
C......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
     &           +FDDD*SS1N3+FDDOD*(SS1N1+SS1N2)
         TEMP_LV_A3=             !LHS VECTOR
C......TRANSIENT TERM (LUMPED)
     &           FDDDODT+2.*FDDODODT

C     VEC...LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
         TEMP_LV_B(IE,1)=TEMP_LV_B1*NCELE !LOAD VECTOR
         TEMP_LV_B(IE,2)=TEMP_LV_B2*NCELE !LOAD VECTOR
         TEMP_LV_B(IE,3)=TEMP_LV_B3*NCELE !LOAD VECTOR
         TEMP_LV_A(IE,1)=TEMP_LV_A1*NCELE !LUMPED LHS MATRIX
         TEMP_LV_A(IE,2)=TEMP_LV_A2*NCELE !LUMPED LHS MATRIX
         TEMP_LV_A(IE,3)=TEMP_LV_A3*NCELE !LUMPED LHS MATRIX
#endif

C     LINES TO RUN ON A SCALAR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR QC AND TRANS_LV_A
C     ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
         TRANS_LV_B(NM1)=TRANS_LV_B(NM1)+TEMP_LV_B1*NCELE !LOAD VECTOR
         TRANS_LV_B(NM2)=TRANS_LV_B(NM2)+TEMP_LV_B2*NCELE !LOAD VECTOR
         TRANS_LV_B(NM3)=TRANS_LV_B(NM3)+TEMP_LV_B3*NCELE !LOAD VECTOR
         TRANS_LV_A(NM1)=TRANS_LV_A(NM1)+TEMP_LV_A1*NCELE !LUMPED LHS MATRIX
         TRANS_LV_A(NM2)=TRANS_LV_A(NM2)+TEMP_LV_A2*NCELE !LUMPED LHS MATRIX
         TRANS_LV_A(NM3)=TRANS_LV_A(NM3)+TEMP_LV_A3*NCELE !LUMPED LHS MATRIX
#endif

      ENDDO

C     LINES TO RUN ON A VECTOR COMPUTER
C     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR QC, TRANS_LV_A
#ifdef CVEC
      DO IE=1,NE
        NM1=NM(IE,1)
        NM2=NM(IE,2)
        NM3=NM(IE,3)
        TRANS_LV_A(NM1)=TRANS_LV_A(NM1)+TEMP_LV_A(IE,1) !LUMPED LHS MATRIX
        TRANS_LV_A(NM2)=TRANS_LV_A(NM2)+TEMP_LV_A(IE,2) !LUMPED LHS MATRIX
        TRANS_LV_A(NM3)=TRANS_LV_A(NM3)+TEMP_LV_A(IE,3) !LUMPED LHS MATRIX
        TRANS_LV_B(NM1)=TRANS_LV_B(NM1)+TEMP_LV_B(IE,1) !LOAD VECTOR
        TRANS_LV_B(NM2)=TRANS_LV_B(NM2)+TEMP_LV_B(IE,2) !LOAD VECTOR
        TRANS_LV_B(NM3)=TRANS_LV_B(NM3)+TEMP_LV_B(IE,3) !LOAD VECTOR
      END DO
#endif

C.... SOLVE FOR C*H NODE BY NODE

      DO I=1,NP
         NCI=NODECODE(I)
         IF(NCI.NE.0) CH1(I)=TRANS_LV_B(I)/TRANS_LV_A(I)
C     IF(LBArray_Pointer(I).NE.0) CH1(I)=0.d0  !ESSENTIAL C=0 BOUNDARY CONDITION
      END DO

      RETURN
      END


C****************************************************************************************
C   Subroutine to compute 2D SigmaT fields from 3D salinity and/or temperature fields   *
C                                                                                       *
C                                    R.L.  6/22/05                                      *
C****************************************************************************************
C

      SUBROUTINE CALC_SIGMAT_2D ()

      USE GLOBAL

      INTEGER NH

      IF(ABS(IDEN).EQ.2) THEN
        DO NH=1,NP
        ENDDO
      ELSEIF(ABS(IDEN).EQ.3) THEN
        DO NH=1,NP
        ENDDO
      ELSEIF(ABS(IDEN).EQ.4) THEN
        DO NH=1,NP
        ENDDO
      ENDIF

      RETURN
      END


c******************************************************************************
!Kendra: Eliminated the bpg calculation in the vsmy subroutine and added a 
!        subroutine for the bpg calculation v45.12
c******************************************************************************
c  SUBROUTINE BPG3D                                                           *
c                                                                             *
c  Note, the following time stepping coefficients are computed in             *
C     VSSTUP and passed in a common block.                                    *
c                                                                             *
c  IDTAlp1      = I*DelT*Alp1        - weights coriolis term in LHS matrix    *
c  IDT1MAlp1    = I*DelT*(1.-Alp1)   - weights coriolis term in RHS forcing   *
c  DTAlp3       = DelT*Alp3          - weights vert diff term in LHS matrix   *
c  DT1MAlp3     = DelT*(1-Alp3)      - weights vert diff term in RHS forcing  *
c  DTAlp2       = DelT*Alp2          - weights bot stress term in LHS matrix  *
c  DT1MAlp2     = DelT*(1.-Alp2)     - weights bot stress term in RHS forcing *
c                                                                             *
c  q(MNP,MNodes) - 3D Complex Velocity field (GAMMA) from past time step.     *
c                                                                             *
c                                                                             *
c  NH - horizontal node counter                                               *
c  NP - number of nodes in horizontal grid                                    *
c  NFEN - number of nodes in the vertical grid                                *
c  BTP - total barotropic pressure (atmos press, water level, tidal potential)*
c                 at time levels s+1/2                                        *
c******************************************************************************

      SUBROUTINE BPG3D()

Casey: Added the following variable declarations from GLOBAL.
C
      USE GLOBAL, ONLY: IFNLFA, IFNLCT, NODELE, NM,
     &      NEITABELE, NODECODE, NOFF, VIDBCPDXOH, VIDBCPDYOH,
     &      SFAC

      USE GLOBAL_3DVS
#ifdef CMPI
      USE MESSENGER
      IMPLICIT NONE
      REAL(SZ) :: DUMV1(1),DUMV2(1)
#endif

Casey: Added the following local variable declarations.
C
      INTEGER :: NCELE
      INTEGER :: TEMPNCELE
      INTEGER :: TEMPSTOP

      INTEGER :: NEle           !local value of NetTabEle
      INTEGER :: k              !vertical node loop counter (1-bottom, NFEN-surf)
      INTEGER :: NH             !horizontal node loop counter
      INTEGER :: N              !neighbor node loop counter
      INTEGER :: N1,N2,N3,NNFirst !local node numbers used to compute gradients
      INTEGER :: NN             !output loop counter

      REAL(8) :: Hs             !Total water depth at time level s
      REAL(8) :: HsOAMB         !Hs/(a-b)
      REAL(8) :: HsHsOAMBAMB    !(Hs/(a-b))^2
      REAL(8) :: Hsp1           !Total water depth at time level s+1
      REAL(8) :: Hsp1OAMB       !Hsp1/(a-b)
      REAL(8) :: Hsp1Hsp1OAMBAMB !(Hsp1/(a-b))^2

      REAL(SZ) :: Zk            !z depth of any node k in the vertical
      REAL(SZ) :: DelSig        ! sigma(k+1)-sigma(k)
      REAL(SZ) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
      REAL(SZ) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
      REAL(SZ) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
      REAL(SZ) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
      REAL(SZ) :: SigmaNN       !Sigma value of a neighbor node
      REAL(SZ) :: HsN2          !Depth value of a neighbor node

      REAL(SZ) :: SFacAvg       ! kmd48.33bc add in spherical factors

      REAL(SZ) :: BCPN1,BCPN2,BCPN3,BCPNFirst !nodal values of BCP
      REAL(SZ) :: BCPDX2A,BCPDY2A !(Horiz. grads of BCP)*2*Element Area
      REAL(SZ) :: SigTAvg       !avg SigT between 2 vertical nodes
      REAL(SZ) :: HGORhoOAMB    !depth*gravity/(reference density)/(a-b)

      REAL(SZ) :: a1,a2,a3,b1,b2,b3
      REAL(SZ) :: TotalBCPGArea2

      REAL(SZ) :: DBCPDX2A  
      REAL(SZ) :: DBCPDY2A  
      COMPLEX(SZ) :: BCPG(NFEN)   !baroclinic pressure gradient
      COMPLEX(SZ) :: VIBCPG         !baroclinic pressure gradient
C     
C     INCREMENT THE TIMESTEP SINCE START COUNTER
C     

C*************************************************************************************
C     Check whether it is time to print various 3D outputs      

C     
C     If a baroclinic run, compute the 3D baroclinic pressure field
C     The buoyancy field is defined as
C     BCP(z)    =(gravity/rho ref)*          integral (SigT) from surface down to z
C     BCP(sigma)=(gravity/rho ref)*(H/(a-b))*integral (SigT) from a down to sigma
C     where 
C     SigT = Sigma T = Rho - 1000 = density - 1000
C     SigT0 = Sigma t value of reference density (typically = 0)
C     Sigma = dimensionless vertical coordinate
C     
      IF((IDEN.GE.1).OR.(IDEN.LE.-1)) THEN
         DO NH=1,NP             !loop over horizontal nodes

Casey: Changed "NolIFA" to "IFNLFA."
C
            Hs=DP(NH)+IFNLFA*Eta2(NH) !total depth at previous (s) timestep
            HGORhoOAMB=GORhoOAMB*Hs !(gravity/rho ref)*(H/(a-b))
            BCP(NH,NFEN)=0.d0
            DO k=NFEN-1,1,-1    !loop over vertical nodes, starting at top and working down
               SigTAvg=(SigT(NH,k+1)+SigT(NH,k))/2.d0
               DelSig=Sigma(k+1)-Sigma(k)
               BCP(NH,k)=BCP(NH,k+1)+HGORhoOAMB*(SigTAvg-SigT0)*DelSig
            ENDDO
         ENDDO
#ifdef CMPI
C     Update BCP on ghost nodes
C      CALL UPDATER3D(BCP) !!!!Don't know if this is needed at this time
#endif
      ENDIF


C*************************************************************************************
C     Compute 3D baroclinic pressure gradients    

C     
C     Loop over each horizontal node to compute the horizontal velocity
C     
      DO NH=1,NP                !loop over horizontal nodes

         Hs  = DP(NH)+IFNLFA*Eta2(NH) !Total depth at previous (s) timestep
         HsOAMB=Hs/AMB

c     Zero out baroclinic pressure gradient and vertically integrated
c     baroclinic pressure gradient for a barotropic run

         IF (IDEN.EQ.0) THEN
            DO k=1,NFEN
               BCPG(k)=(0.d0,0.d0)
            END DO
            VIDBCPDXOH(NH)=0.d0
            VIDBCPDYOH(NH)=0.d0
         ENDIF

c     Start computing baroclinic terms

         IF ((IDEN.GE.1).OR.(IDEN.LE.-1)) THEN

c     Start computing baroclinic pressure gradient (computed in level
c     coordinates) at each node in the vertical

            DO k=1,NFEN

               DBCPDX2A=0.d0
               DBCPDY2A=0.d0
               TotalBCPGArea2=0.d0
               N1=NH
               BCPN1=BCP(NH,k)

               Zk=HsOAMB*(Sigma(k)-B)-DP(NH) !determine z corresponding to sigma level k
               N2=NEITAB(NH,2)  !operate on 1st neighbor

Casey: Changed "NolIFA" to "IFNLFA."
C
               HsN2=DP(N2)+IFNLFA*Eta2(N2)
               SigmaNN=B+AMB*(Zk+DP(N2))/HsN2 !equivalent sigma value at neighbor
               CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
               NNFirst=N2       !save these values until end
               BCPNFirst=BCPN2  !save these values until end

               DO N=3,NNeigh(NH) !operate on rest of neighbors
                  N3=N2         !shift previously computed values
                  BCPN3=BCPN2   !shift previously computed values
                  N2=NeiTab(NH,N) !select new neighbor to work on

Casey: Changed "NolIFA" to "IFNLFA."
C
                  HsN2=DP(N2)+IFNLFA*Eta2(N2)
                  SigmaNN=B+AMB*(Zk+DP(N2))/HsN2 !equivalent sigma value at neighbor
                  CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
                  NEle=NeiTabEle(NH,N-2) !element # defined by nodes NH,NN2,NN1
                  ! jgf49.58: The NeiTabEle matrix is semi sparse; see
                  ! Casey's comments in vsmy.F. NOFF array lookups fail
                  ! if NEle comes back 0, so I will just cycle to the next
                  ! N value here if that happens ... TODO: somebody please
                  ! confirm that this is the right answer here.
                  IF (NEle.eq.0) THEN
                     CYCLE
                  ENDIF
Casey: Added the computation of "NCELE" and the last part of the IF statement.
C
                  NCELE = NODECODE(NH)*NODECODE(N2)
     &                *NODECODE(N3)*NOFF(NELE)
                  IF((BCPN2.NE.-999.).AND.(BCPN3.NE.-999.)
     &                 .AND.(NEle.NE.0).AND.(NCELE.NE.0)) THEN !if all 3 nodes are active, compute bu
                     TotalBCPGArea2=TotalBCPGArea2+Areas(NEle)
C    kmd48.33bc add in spherical factors
                     SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                     a1=X(N3)-X(N2)
                     a2=X(N1)-X(N3)
                     a3=X(N2)-X(N1)
                     b1=(Y(N2)-Y(N3))*SFacAvg
                     b2=(Y(N3)-Y(N1))*SFacAvg
                     b3=(Y(N1)-Y(N2))*SFacAvg
                     DBCPDX2A=DBCPDX2A+(BCPN1*b1+BCPN2*b2+BCPN3*b3)
                     DBCPDY2A=DBCPDY2A+(BCPN1*a1+BCPN2*a2+BCPN3*a3)
                  ENDIF
               END DO

               N3=N2            !wrap back to beginning to get final contributio
               N2=NNFirst
               BCPN3=BCPN2
               BCPN2=BCPNFirst
               NEle=NeiTabEle(NH,NNeigh(NH)-1)

Casey: Added the computation of "NCELE" and the last part of the IF statement.
C
               ! jgf49.58 NOFF lookups fail if NELE comes back 0.
               IF (NELE.ne.0) THEN
                  NCELE = NODECODE(NH)*NODECODE(N2)
     &                *NODECODE(N3)*NOFF(NELE)
               ENDIF
               IF((BCPN2.NE.-999.).AND.(BCPN3.NE.-999.)
     &              .AND.(NEle.NE.0).AND.(NCELE.NE.0)) THEN
                  TotalBCPGArea2=TotalBCPGArea2+Areas(NEle)
C    kmd48.33bc add in spherical factors
                  SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  a1=X(N3)-X(N2)
                  a2=X(N1)-X(N3)
                  a3=X(N2)-X(N1)
                  b1=(Y(N2)-Y(N3))*SFacAvg
                  b2=(Y(N3)-Y(N1))*SFacAvg
                  b3=(Y(N1)-Y(N2))*SFacAvg
                  DBCPDX2A=DBCPDX2A+(BCPN1*b1+BCPN2*b2+BCPN3*b3)
                  DBCPDY2A=DBCPDY2A+(BCPN1*a1+BCPN2*a2+BCPN3*a3)
               ENDIF

C    kmd48.33bc changed the BPG calculation for bottom boundary issues
               IF(TotalBCPGArea2.EQ.0.) THEN 
!                  IF (k.eq.NFEN) THEN !kd46.01
                     BCPG(k)=(0.d0,0.d0)
!                  ELSE
!                     BCPG(k)=BCPG(k+1)
!                  END IF 
               ELSE
                  BCPG(k)=(DBCPDX2A+iy*DBCPDY2A)/TotalBCPGArea2
               ENDIF

            ENDDO

C    kmd48.33bc added this to BPG
            DO k=NFEN,1,-1
               IF ((BCPG(k).EQ.(0.d0,0.d0)).AND.(k.NE.NFEN)) THEN
                   BCPG(k)=BCPG(k+1)
               END IF
            END DO


            DO k=1,NFEN
               BPG(NH,k) = BCPG(k) 
            END DO

c     Finished computing baroclinic pressure gradient (computed in level
c     coordinates) at each node in the vertical


c     Compute vertically integrated baroclinic pressure gradient for use
c     in the wave equation.  NOTE: For a prognostic model in which the
c     density field evolves in time, this calculation should be done
c     after the new density field is computed.  In this case one would
c     integrate over the vertical first and differentiate second.

            VIBCPG=(0.d0,0.d0)
            DO k=NFEN-1,1,-1
               VIBCPG=VIBCPG+0.5d0*(BCPG(k+1)+BCPG(k))
     &              *(Sigma(k+1)-Sigma(k))
            ENDDO
            VIDBCPDXOH(NH)=REAL(VIBCPG)/AMB
            VIDBCPDYOH(NH)=AIMAG(VIBCPG)/AMB

         ENDIF

c     Finished computing baroclinic terms

      ENDDO

C     Finish loop over horizontal nodes to compute the horizontal velocity

#ifdef CMPI
C***********************************************************************
C     Update new 3D baroclinic pressure gradient and the vertically
C     integrated baroclinic pressure gradient on ghost nodes
C     
      CALL UPDATEC3D(BPG)    
      CALL UPDATER(VIDBCPDXOH,VIDBCPDYOH,DUMV1,2) 
#endif

      RETURN 
      END SUBROUTINE BPG3D



C*************************************************************************
C     Subroutine to interpolate baroclinic pressure (BCP) to a specified
C     sigma value (SigmaNN) given an initial guess of which sigma
C     level is closest to the specified value.
C
C                                    R.L.  5/04/01
C                                    R.L.  5.19/03
C*************************************************************************
C
      SUBROUTINE ZSURFBUOY(SigmaNN,BCPressNN,NN,J)

      USE GLOBAL_3DVS
      IMPLICIT NONE
      REAL(SZ) :: BCPressNN
      REAL(SZ) :: SigmaNN     !Sigma value of a neighbor node
C     jgf46.00 Explicitly declared the following variables
      INTEGER NN
      INTEGER J
      INTEGER LBelo
      INTEGER LAbov
      INTEGER LTry
      INTEGER IDiag
      REAL(SZ) SigBelo
      REAL(SZ) SigAbov
      REAL(SZ) SigTry

      IDiag=0 

      IF(SigmaNN.LE.1.0001*b) THEN !if into ground then skip
         SigBelo=-999
         SigAbov=-999
         BCPressNN=-999.
         GOTO 100
      ENDIF
      IF((SigmaNN.GT.1.0001*b).AND.(SigmaNN.LE.b)) THEN !at bottom then use bottom
         LBelo=1
         BCPressNN=BCP(NN,LBelo)
         SigBelo=b
         SigAbov=b
         GOTO 100
      ENDIF
      IF(SigmaNN.GE.a) THEN     !into air use surface
         LAbov=NFEN
         BCPressNN=BCP(NN,LAbov)
         SigBelo=a
         SigAbov=a
         GOTO 100
      ENDIF

      LTry=J                    !start search for SIGABOV and SIGBELO
      SigTry=Sigma(LTry)                    
      IF(SigmaNN.GT.SigTry) THEN !too low
         SigBelo=SigTry         !SIGBELO may = SIGTRY
         LBelo=LTry
         LTry=LTry+1            !look at next level higher
 90      SigTry=Sigma(LTry)
         IF(SigmaNN.GT.SigTry) THEN !still too low
            SigBelo=SigTry
            LBelo=LTry
            LTry=LTry+1
            GOTO 90 
         ENDIF
         SigAbov=SigTry         !found upper bracketing sigma
         LAbov=LTry
         GOTO 99                !go interpolate
      ENDIF
      IF(SigmaNN.LE.SigTry) THEN !to high
         SigAbov=SigTry         !SIGABOV may = SIGTRY
         LAbov=LTry                          
         LTry=LTry-1            !look at next level lower
 91      SigTry=Sigma(LTry)
         IF(SigmaNN.LE.SigTry) THEN !still too high
            SigAbov=SigTry
            LAbov=LTry
            LTry=LTry-1
            GOTO 91
         ENDIF
         SigBelo=SigTry         !found lower bracketing sigma
         LBelo=LTry
      ENDIF

 99   BCPressNN=(BCP(NN,LAbov)-BCP(NN,LBelo)) !interpolation
     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + BCP(NN,LBelo)

 100  CONTINUE

      IF(IDiag.EQ.2) THEN
         WRITE(2,*) '******** ZSURFBUOY **********'
         WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov',
     &        '      BCPressNN'
         WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,BCPressNN
 777     FORMAT(I7,I5,3(F10.3),E14.5)
      ENDIF

      RETURN
      END


