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
C     CPRECOR          - Use the predictor-corrector algorithm for GWCE       *
C                        and momentum equations (package deal)                *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C   See header.f for a summary history of code modifications.                 *
C******************************************************************************
C
      SUBROUTINE TIMESTEP(IT,TimeLoc)
C
#ifdef IEEE_DEBUG
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      USE SIZES, ONLY : MNP, MNWPROH, MYPROC
      USE GLOBAL, ONLY : UU1, VV1, QX1, QY1, QX2, QY2, UU2, VV2, MOM_LV_X,
     &    MOM_LV_Y, WSX2, WSY2, PR1, PR2, TiP1, TiP2, UU0, VV0, QY0, 
     &    GWCE_LV, TRANS_LV_B, TRANS_LV_A, WSX1, WSY1, dp1, dp2, QN0, 
     &    EN0, EN1, EN2, CICE1, CICE2, CICEOUT, RSNX1, RSNX2, RSNY1, RSNY2,
     &    QX0, QN1, QN2, RSNXOUT, RSNYOUT, LNM_BC1, LNM_BC2, H1, H2,
     &    QNIN1, QNIN2, ENIN1, ENIN2, EtaDisc, NIBNODECODE,
     &    BCHGTIMINC, ETA2,ElevDisc, BTIME1, BTIME2, BTIME_END, C2D_PTrans,
     &    C2DDI, C3D, C3DVS, CBAROCLINIC, CGWCE_New, CHOTHS, CICE_TIME1,
     &    CICE_TIME2, CICE_TIMINC, BTIMINC, CME_New_C1, CME_New_C2, 
     &    CME_New_NC, CPRECOR, CTiP, DRAMP, DRampElev, DRampExtFlux, 
     &    DRampIntFlux, DRampMete, DRampTip, DRampWRad, DTDP, DTDPHS,
     &    DUnRampMete, ENDWAVE, ErrorElev, Flag_ElevError, FluxSettlingIT,
     &    FluxSettlingTime, FTIMINC, IFNLFA, H0, IHOT, IM,
     &    inundationOutput, ITHS, METONLY, NCICE, NDDT, NFFR, NHSINC,
     &    NHSTAR, NOLIFA, NRamp, NRS, NSCREEN, NTIF, QTIME1, QTIME2, Ramp,
     &    RampElev, RampExtFlux, RampMete, RampTip, RampWRad, RBCRATIO,
     &    RBCTIME1, RBCTIME2, RBCTIMEINC, REFTIM, RES_BC_FLAG, 
     &    RSTIME1, RSTIME2, NWS, RSTIMINC, STATIM, TK,
     &    WVNXOUT, WVNYOUT, WarnElev, PERT, FACET, FFT, AMIGT, SALTPHA,
     &    SALTAMP, NODECODE, scratchMessage, TPK, NOFF, DASigT, FFF, 
     &    QNPH, QNAM, ENPH, ENAM, ETRF, FPER, FFACE, ETA1, FAMIG, LNM_BC,
     &    BCFLAG_LNM, usingDynamicWaterLevelCorrection,
     &    dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, L_N, TKM, 
     &    WarnElevDump, nodes_lg, NT, NM, NB, DEBUG, NHOUTONCE,
     &    Limit_WaveStressGrad,
     &    earlyterminate, IT_foundNaN, WarnVel, Flag_VelError,
     &    dynamicWaterLevelCorrectionDelta, RampIntFlux
      use mod_logging, only: ScreenUnit, INFO, WARNING,
     &    ERROR, logMessage, allMessage, setMessageSource, unsetMessageSource
      USE ADC_CONSTANTS, ONLY: rhoWat0, sigT0, waveWindMultiplier,
     &                         WaveStressGrad_Cap, G
     
      USE GWCE, ONLY : solveGWCE, numitr
      USE MOMENTUM, ONLY : solveMomentumEq
      USE WETDRY, ONLY : computeWettingAndDrying
      USE WRITE_OUTPUT, ONLY : writeOutput2D, writeHotStart, 
     &    writeWarnElev, collectInundationData, collectMinMaxData
      USE MESH, ONLY : NE, NP, DP, SLAM, SFEA, ICS, TotalArea,
     &    MJU, Areas, SFAC, NNeigh, NeiTab, X, Y,
     &    LBARRAY_POINTER
      USE BOUNDARIES, ONLY : NOPE, NETA, NBOU, NVEL, LBCODEI, NBV, SIII,
     &    NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, NFLUXIBP, NFLUXRBC, CSII,
     &    NFLUXIB64,
     &    BARLANHT, BARLANCFSP, NVELL, IBCONN, BARINHT, BARINCFSP,
     &    BARINCFSB, PIPEHT, PIPECOEF, PIPEDIAM, ISSUBMERGED64, ISSUBMERGED64P
      USE GLOBAL_IO, ONLY : nddt2get
      USE HARM, ONLY : updateHarmonicAnalysis
      USE WIND, ONLY : getMeteorologicalForcing, PRBCKGRND_MH2O, rsget,
     &    getdynamicWaterLevelCorrections,
     &    dynamicWaterLevelCorrectionMassAdjust,
     &    CLOSE_MET_FILES

      USE mod_terminate, only: terminate
C.... TCM V49.64.01 ADDITIONS FOR ICE
      USE OWI_ICE,ONLY : NCICE1_INIT,NCICE1_GET
C.....sb46.28sb03 added 09/xx/2006
      USE RS2,ONLY : RS2INIT,RS2GET
      USE NodalAttributes, ONLY : 
     &     Apply2DBottomFriction,
     &     Apply2DInternalWaveDrag,
     &     Apply3DBottomFriction,
     &     Apply2DMomentumDisp,
     &     GeoidOffset,
     &     BFCdLLimit, fric, fgamma, ftheta, manningsn, 
     &     hbreak, loadManningsN, iflinbf, ifnlbf, ifhybf
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, NOutGS, 
     &    writeFort066, writeFort067, writeFort065, 
     &    readFort020, readFort021, readFort019 
c. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
      USE GLOBAL_3DVS, ONLY:
     &     A, B, BSX, BSY, EVTOT, ISLIP, KP, Q, SIGMA, Z0B,NFEN
! arash 1/20/2016
     &     , BPG
      USE ADCIRC_Mod, ONLY: ITIME_BGN

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
#endif
      USE WEIR_FLUX, ONLY: COMPUTE_EXTERNAL_BOUNDARY_FLUX, COMPUTE_INTERNAL_BOUNDARY_FLUX,
     & BARINHT1, BARINHT2, BARLANHT1, BARLANHT2, COMPUTE_INTERNAL_BOUNDARY64_FLUX,
     & COMPUTE_CROSS_BARRIER_PIPE_FLUX, SET_SUBMERGED64_AT
      USE SPONGELAYER, ONLY: LoadAbsLayerSigma, sponge_dis_mthd, sponge_opsplit0,
     &     adjust_sponge_sigma, getabslayerext, sponge_shift_soln

      USE subgrid, ONLY: level0, getVertLookup

      ! full luni-solar tidal potential !
      USE mod_tidepotential, only: tidePotential

      IMPLICIT NONE
      INTEGER, intent(in) :: IT
C
c. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
      INTEGER IE, I, J, K, II, JJ                  !local loop counters
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle
      INTEGER NCyc
      INTEGER :: itest
      INTEGER :: kemax
      INTEGER :: kvmax
      INTEGER :: nnbb
      REAL(8) :: velabs
      REAL(8) :: vel
      REAL(8) :: velmax
      INTEGER :: errorElevExceeded
      INTEGER :: warnVelExceeded
      INTEGER :: warnElevExceeded
      logical, save ::  EtaDisc_Fill = .TRUE.

      REAL(8) CCSFEA
      REAL(8) ElMax
      REAL(8) PIPE_FLUX
      REAL(8) H2N1, H2N2, H2N3
      REAL(8) EtaN1,EtaN2,EtaN3
      REAL(8) QTRatio
      REAL(8) RStRatio, RSX, RSY
      REAL(8) RSXLIM,RSYLIM,RSXYMAG,RSSCALE
      REAL(8) SAltMul, S2SFEA
      REAL(8) TPMul
      REAL(8) UN1

      REAL(8) AreaIE2,AreaEle
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8),intent(out) :: TimeLoc
      REAL(8) TimeLocFluxSettling
      REAL(8) TimeH
      REAL(8) runperccomplete
C   kmd48.33bc - added in the heat flux variables
      CHARACTER(80) :: CDUM80
      INTEGER :: NumofBCNode
      REAL(8), SAVE :: StaTimHS, RefTimHS
      INTEGER :: NOD

C   TCM V49.64.01 -- ADDED FOR ICE CONCENTRATION FIELDS
      REAL(8) CICE_TRatio  !ICE VARIABLES

C   TCM v50.66.02 -- Added for Time Varying Bathymetry
      INTEGER :: IJK
      REAL(8) ETA2TMP,DPTMP,DPTMP2,BTRATIO  !tcm v50.66.01 bathymetry changes

      REAL(8) :: ARG, ARGJ, RFF
C....... DW
      REAL(8):: RAMPTMP
      LOGICAL, SAVE:: FIRST_ADJ_ABS_LAYER = .TRUE. 
      LOGICAL, SAVE:: FIRST_ENTER = .TRUE. 
      INTEGER:: NBDI, istep
C.......
      LOGICAL ISFRONT
      LOGICAL foundNaN
      INTEGER NNBB1, NNBB2, NNBBNEI, NWETNEI
      LOGICAL UNSUBMERGE
      REAL(8) X1, X2, Y1, Y2, ET1, ET2, LEN, SLP, HTOT
      
      IF (IT .EQ. 1) THEN
!.......DMW202401 On first timestep of a cold start, find H2 at start of timestepper since ETA2 = ETA1
        DO I = 1,NP
           H2(I) = DP(I)+IFNLFA*ETA2(I)
           H1(I) = H2(I)
        ENDDO
      ENDIF

      call setMessageSource("timestep")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


C    kmd48.33bc - changed for timestep changes in the hot start files
C     jgf46.21 Combined flux/radiation b.c. for rivers
#ifdef IBM
      IF (CHOTHS.eqv..true.) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
      ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP,KIND(0.0D0))
      END IF
#else
      IF (CHOTHS.eqv..true.) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
      ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP)
      END IF
#endif

C...  COMPUTE MASTER TIME WHICH IS REFERENCED TO THE BEGINNING TIME OF
C...  THE MODEL RUN
C...
C...  TIMEH: HARMONIC CALCULATIONS ARE MADE FOR TIME WHICH INCLUDES THE REFTIM
C...  TO ALLOW FOR THE POSSIBILITY THAT THE EQUILIBRIUM ARGUMENTS MAY
C...  BE FOR A TIME OTHER THAN THE MODEL STARTING TIME.
C...  
C...  Computation of TimeLoc and TimeH are consolidated. SB,ZC,RL 2022-06-22
C...
      TimeLoc=IT*DTDP + StaTim*86400.D0
      TimeH=IT*DTDP + (StaTim - RefTim)*86400.D0
      IF (CHOTHS.eqv..true.) THEN
         IF ((ITHS+1).EQ.IT) THEN
           StaTimHS=((IT-1)*DTDPHS)/86400.D0
           RefTimHS=((IT-1)*DTDP)/86400.D0
         END IF
         TimeLoc=IT*DTDP + (StaTimHS - RefTimHS)*86400.D0
         TimeH=IT*DTDP + ((StaTimHS - StaTim) - RefTim)*86400.D0
      END IF

C...  HARMONIC CALCULATIONS ARE MADE FOR TIME WHICH INCLUDES THE REFTIM
C...  TO ALLOW FOR THE POSSIBILITY THAT THE EQUILIBRIUM ARGUMENTS MAY
C...  BE FOR A TIME OTHER THAN THE MODEL STARTING TIME.
C...
      TimeH=IT*DTDP + (StaTim - RefTim)*86400.D0
!kmd - added this for the cases where the timestep changes from hot start timestep
      IF (CHOTHS.eqv..true.) THEN
         IF ((ITHS+1).EQ.IT) THEN
           StaTimHS=((IT-1)*DTDPHS)/86400.D0
           StaTim=((IT-1)*DTDP)/86400.D0
         END IF
         TimeH=IT*DTDP + ((StaTimHS - StaTim) - RefTim)*86400.D0
      END IF

C-------------------TEMPORAL RAMP-------------------------------------
C...  THIS SECTION IS REWRITTEN BY SB, ZC, and RL 2022-06-22
C...
C...  NOTE: MOVED HERE BY DW/WJP : May 3rd 2018
C...  DEFINE Ramp FUNCTION FOR BOUNDARY ELEVATION FORCING, WIND AND PRESSURE
C.... FORCING AND TIDAL POTENTIAL FORCING
C...
C
C     jgf46.08 Calculate ramp functions.
C     jgf46.21 Modify to match behavior of 46.02
      IF(NRamp.EQ.0) THEN
         Ramp=1.0D0
         RampExtFlux=1.0D0
         RampIntFlux=1.0D0
         RampElev=1.0D0
         RampTip=1.0D0
         RampMete=1.0D0
         RampWRad=1.0D0
      ELSE
         Ramp=TANH((2.D0*TimeLoc/86400.D0)/DRamp)
         RampExtFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampExtFlux)
         RampIntFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampIntFlux)
         RampElev=TANH((2.D0*TimeLoc/86400.D0)/DRampElev)
         RampTip=TANH((2.D0*TimeLoc/86400.D0)/DRampTip)
         RampMete=TANH((2.D0*(TimeLoc/86400.D0-DUnRampMete))/DRampMete)
         RampWRad=TANH((2.D0*TimeLoc/86400.D0)/DRampWRad)
      ENDIF
C
C     jgf46.21 If there is an external flux (i.e. river) boundary, turn
C     off all forcings except the river flux forcing for the duration of
C     the FluxSettlingTime. When the FluxSettlingTime has ended, turn
C     all forcings back on.
      IF (NRamp.gt.1.AND.FluxSettlingIT.GT.0) THEN
         IF(IT.LT.(FluxSettlingIT+10)) THEN
            Ramp=0.0
            RampIntFlux=0.0
            RampElev=0.0
            RampTip=0.0
            RampMete=0.0
            RampWRad=0.0
         ELSE
            TimeLocFluxSettling=(IT-FluxSettlingIT-10)*DTDP
     &           + StaTim*86400.D0
            IF (CHOTHS.eqv..true.) THEN
               TimeLocFluxSettling=(IT-FluxSettlingIT-10)*DTDP
     &           + (StaTimHS - RefTimHS)*86400.D0
            END IF

            Ramp=TANH((2.D0*TimeLocFluxSettling/86400.D0)/DRamp)
            RampIntFlux=TANH((2.D0
     &           *TimeLocFluxSettling/86400.D0)/DRampIntFlux)
            RampElev=TANH((2.D0
     &           *TimeLocFluxSettling/86400.D0)/DRampElev)
            RampTip=TANH((2.D0
     &           *TimeLocFluxSettling/86400.D0)/DRampTip)
            IF(NRamp.ne.8) then
               RampMete=TANH((2.D0
     &           *TimeLocFluxSettling/86400.D0)/DRampMete)
            ENDIF
            RampWRad=TANH((2.D0
     &           *TimeLocFluxSettling/86400.D0)/DRampWRad)
         ENDIF
      ENDIF

      !jgf49.44: Cover the case where the ramp length is zero.
      IF (DRamp.lt.1.0e-6) Ramp = 1.0d0
      IF (DRampExtFlux.lt.1.0e-6) RampExtFlux = 1.0d0
      IF (DRampIntFlux.lt.1.0e-6) RampIntFlux = 1.0d0
      IF (DRampElev.lt.1.0e-6) RampElev = 1.0d0
      IF (DRampTip.lt.1.0e-6) RampTip = 1.0d0
      IF (DRampMete.lt.1.0e-6) RampMete = 1.0d0
      IF (DRampWRad.lt.1.0e-6) RampWRad = 1.0d0
C     
C------DW
      IF ( LoadAbsLayerSigma ) THEN
         IF (  FIRST_ADJ_ABS_LAYER ) THEN
            CALL Adjust_Sponge_Sigma( abs(DTDP) ) ;

            FIRST_ADJ_ABS_LAYER = .FALSE. ;
         END IF
      END IF   
C-----END DW

C-------------------TEMPORAL RAMP-------------------------------------


C...  SHIFT THE FLUX PER UNIT WIDTH, DEPTH AVERAGED VELOCITIES, BOTTOM STRESS,
C...  WIND STRESS, SURFACE PRESSURE AND TIDAL POTENTIALS TO PREVIOUS TIME STEP.
C...  ZERO OUT THE NEW FORCING TERMS AND LOAD VECTORS
C...  COMPUTE A NEW BOTTOM FRICTION COEFFICIENT
C...
ckmd  Shift values in time for predictor-corrector algorithm

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
            PR2(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
         ELSE
C WJP 02.26.2018 inserted from Dam's code    
!c dw ------------- otherwise, tide run could fail      
            PR2(I) = PRBCKGRND_MH2O 
         END IF

C     TIP..Tidal potential forcing
         if(CTIP) then
            TIP1(I)=TIP2(I)
            TIP2(I)=0.D0
         endif

      END DO
      !
      ! jgf: shift water level offset values
      if (usingDynamicWaterLevelCorrection.eqv..true.) then
         dynamicWaterLevelCorrection1(:) = dynamicWaterLevelCorrection2(:)
         dynamicWaterLevelCorrection2(:) = 0.d0
      endif
C
C.............................DW, sponge layer ..................
      ! Get 
      !   Eta2_AbsLayer, uu2_AbsLayer, vv2_AbsLayer
      IF ( LoadAbsLayerSigma ) THEN
         !
         IF ( .NOT. FIRST_ENTER ) THEN
            !c u_{ext}^{n + 1} -- > eta2, 
            CALL GETABSLAYEREXT( IT, TimeLoc, TimeH, RampElev,
     &           GeoidOffset ) 
         ELSE

            RampTMP = RampElev ;
            DO istep = 2, 1, -1
               IF ( NRamp > 0 ) THEN
                  RampTMP=TANH((2.D0*(TimeLoc - istep*DTDP)/86400.D0)/DRampElev)
               END IF
               RampTMP = max(0.D0,RampTMP) ;                   
               
               !  u_{ext}^{n + 1 - istep} --> uu2
               CALL GETABSLAYEREXT( IT - istep, TimeLoc - istep*DTDP,
     &              TimeH - istep*DTDP, RampTMP,
     &              GeoidOffset ) 

               CALL SPONGE_SHIFT_SOLN( ) ; 
               ! eta1 --> eta0, eta2 --> eta1, u2 --> u1          
            END DO
            !c u_{ext}^{n + 1} -- > eta2, 
            CALL GETABSLAYEREXT( IT, TimeLoc, TimeH, RampElev,
     &           GeoidOffset ) 
            FIRST_ENTER = .FALSE. ;
         END IF
         !
      END IF
C.... END DW..........................................................
      
C...TCM V50.66.01 -- ADDING TIME DEPENDENT BATHYMETRY
C... so that total water column height is unchanged
C
C  DP is linearly interpolated between DP1 and DP2
C  during the time interval btime1 and btime_end.
C  After btime_end, DP is equal to DP2
C
C      DP1                         DP2
C    btime1                      btime2      btime2 = btime1 + btiminc
C       |---------x-----------------|        btime_end = btime1 + btime_end < btime2
C              btime_end
C

      IF(abs(NDDT).EQ.1) THEN
!        Get a new bathymetry from file if time to do so
         IF(TimeLoc.GT.BTIME2) THEN !determine if bathy file time incr. is exceeded
            BTIME1=BTIME2  !new starting time for this record
            BTIME2=BTIME2+BTIMINC  !new ending time for this record
            BTIME_END = BTIME1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval

            DO I=1,NP
               dp1(I) = dp2(I)  ! move current data to old
            END DO
!!!        go get new record for all nodes
            DO I=1,NP
               READ(141,*) IJK,DP2(IJK)
            ENDDO
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
            IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
               DO I=1,NP
                  IF (DP2(I).LT.H0) DP2(I) = H0
               ENDDO
            ENDIF
            write(scratchMessage,'(A36,1X,E15.8,1X,A5)') 
     &         'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            call allMessage(INFO,scratchMessage)
         ENDIF  !test for updating bathymetry records
C.......If time is during the bathymetry change interval, then update bathymetry
         IF(timeloc.LT.BTIME_END) THEN  !tcm 20150728 changed .LE. to .LT.
            bTRatio=(TimeLoc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I))  !Determine incremental amount to adjust bathymetry from DP1
               DPTMP2 = DP1(I) + DPTMP  !this is what will be the new bathymetry to use
               DPTMP = DP(I)-DPTMP2  !this is now the adjustment in bathymetry for this timestep (how much to adjust eta2 by)
               DP(I) = DPTMP2           !updating bathymetry to new value
               ETA2TMP = ETA2(I)-DPTMP   !subtracting elevation by incremental amount from one step step
               ETA2(I) = ETA2TMP        !updating elevation
            ENDDO !I
         ENDIF
         IF(timeloc.EQ.BTIME_END) THEN
            write(scratchMessage,'(A42,1X,E15.8,1X,A5)')
     &         'BATHYMETRY VALUES ARE NOW FIXED AT TIME =',TIMELOC,'(SEC)'
            call allMessage(INFO,scratchMessage)
            DO I=1,NP
               DPTMP = DP2(I)-DP(I)   !figuring what the incremental amount will be to get to the final bathy value
               DP(I) = DP2(I)         !updating bathymetry to final value
               ETA2TMP = ETA2(I)-DPTMP !subtracting elevation by incremental amount
               ETA2(I) = ETA2TMP       !updating elevation
            ENDDO !I
         ENDIF

      ENDIF  !NDDT = 1

      IF(abs(NDDT).EQ.2) THEN
!        Get a new bathymetry from file if time to do so
         IF(TimeLoc.GT.BTIME2) THEN !determine if bathy file time incr. is exceeded
            BTIME1=BTIME2  !new starting time for this record
            BTIME2=BTIME2+BTIMINC  !new ending time for this record
            BTIME_END = BTIME1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
            DO I=1,NP
               dp1(I) = dp2(I)  ! move current data to old
            END DO

!!!        go get new record for only some nodes, all
!!!        other nodes keep their current value
            CALL NDDT2GET( 141,DP2(:),-99999.d0 )
C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
            IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
               DO I=1,NP
                  IF (DP2(I).LT.H0) DP2(I) = H0
               ENDDO
            ENDIF
            write(scratchMessage,'(A36,1X,E15.8,1X,A5)')
     &         'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            call allMessage(INFO,scratchMessage)
         ENDIF  !test for updating bathymetry records

C.......If time is during the bathymetry change interval, then update bathymetry
         IF(timeloc.LT.BTIME_END) THEN !tcm 20150728 changed .LE. to .LT.
            bTRatio=(TimeLoc-bTIME1)/BCHGTIMINC ! interpolate
            DO I=1,NP
               DPTMP =  btratio*(DP2(I)-DP1(I)) !Determine incremental amount to adjust bathymetry from DP1
               DPTMP2 = DP1(I) + DPTMP  !this is what will be the new bathymetry to use
               DPTMP = DP(I)-DPTMP2  !this is now the adjustment in bathymetry for this timestep (how much to adjust eta2 by)
               DP(I) = DPTMP2 !updating bathymetry to new value
               ETA2TMP = ETA2(I)-DPTMP   !subtracting elevation by incremental amount from one step step
               ETA2(I) = ETA2TMP        !updating elevation
            ENDDO !I
         ENDIF
C         
         IF (timeloc.EQ.BTIME_END) THEN
            write(scratchMessage,'(A42,1X,E15.8,1X,A5)')
     &         'Bathymetry values are now fixed at time =',timeloc,'(sec).'
            call allMessage(INFO,scratchMessage)
            DO I=1,NP
               DPTMP = DP2(I)-DP(I)   !figuring what the incremental amount will be to get to the final bathy value
               DP(I) = DP2(I)         !updating bathymetry to final value 
               ETA2TMP = ETA2(I)-DPTMP !subtracting elevation by incremental amount
               ETA2(I) = ETA2TMP       !updating elevation
            ENDDO !I
         ENDIF

      ENDIF  !NDDT = 2



C
C
C-------------------BOTTOM FRICTION-------------------------------------
C
C     2DDI.Set up the 2D friction coefficient
      if (C2DDI) then
         CALL Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, IFNLFA, NP,
     &                              TK)
         ! WJP 02.24.2018 Get the wave drag coefficients
         ! that are put into TKM (which include TK from bottom friction)
         IF (CBaroclinic) THEN
            CALL Apply2DInternalWaveDrag(NP, TK, TKM, UU1, VV1, DP, IT)
            ! CPB 3/9/2023: vertically averaged momentum dispersion as a
            ! power law function. If NOLIBF =/= 3 then this just returns
            CALL Apply2DMomentumDisp(UU1, VV1, DP, ETA2, IFNLFA, NP,
     &                               TKM)
         ELSE
            CALL Apply2DInternalWaveDrag(NP, TK, TKM)
         ENDIF         
      endif
C..RJW Set up the 3D friction coefficient
      if (C3D) then
         CALL Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, IFNLFA, NP,
     &                              TK, NFEN, Z0B)
      endif
      
C-------------------BOTTOM FRICTION-------------------------------------
C
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

C
C---------------------BEGIN ICE FIELDS----------------------------------------
C...  UPDATE THE ICE CONCENTRATION FIELDS FROM UNIT 25,225,227
c...  TCM V49.64.01 ADDED THE ICE FIELDS SECTION
C
      IF (NCICE.EQ.12) THEN
         IF(TimeLoc.GT.CICE_TIME2) THEN
            CICE_TIME1 = CICE_TIME2
            CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
            DO I=1,NP
              CICE1(I) = CICE2(I)
            END DO
            CALL NCICE1_GET(CICE2,NP)
          ENDIF
      ENDIF

C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
C++ GML added ncice .ne. 17
      if (ncice.ne.0.and.(ncice.ne.14) .and. (ncice.ne.17)) then
         CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
         do i=1,np
            CICEOUT(I) = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
         end do
      endif
C--------------------END ICE FIELDS----------------------------------------
C---------------------WATER LEVEL OFFSET-----------------------------------
C     jgf: update water level correction fields
      if (usingDynamicWaterLevelCorrection.eqv..true.) then
         call getDynamicWaterLevelCorrections(dynamicWaterLevelCorrection2, timeloc)
         if (dynamicWaterLevelCorrectionMassAdjust.eqv..true.) then
            dynamicWaterLevelCorrectionDelta=dynamicWaterLevelCorrection2-dynamicWaterLevelCorrection1
         endif
      endif
C
C---------------------BEGIN MET FORCING----------------------------------------

      call getMeteorologicalForcing(nws, wsx2, wsy2, pr2, 
     &   wvnxout, wvnyout, timeloc)
C
C----------------------END MET FORCING----------------------------------------


C...  UPDATE THE WAVE RADIATION STRESS AND READ IN NEW VALUES FROM
C.... UNIT 23.  APPLY Ramp FUNCTION.  ADD RADIATION STRESS TO WIND
C...  STRESS
C...
C...  NRS=2 was added.  sb46.28sb03 09/xx/2006
C...  TCM v49.48 Restructured the wave stress updates in order to
c...             include NRS=4
#ifdef CSWAN
Casey 090302: Added these lines for coupling winds to SWAN.
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      if (coupwind) then
         select case(abs(nws))
         case(12,16)
            do i=1,np
               swan_wx2(i,2) = wvnxout(i) * waveWindMultiplier 
     &                         / ( rampMete * windMultiplier ) 
               swan_wy2(i,2) = wvnyout(i) * waveWindMultiplier 
     &                         / ( rampMete * windMultiplier ) 
            end do
         case default
            do i=1,np
               swan_wx2(i,2) = wvnxout(i) * waveWindMultiplier / rampMete
               swan_wy2(i,2) = wvnyout(i) * waveWindMultiplier / rampMete
            end do
         end select
      endif
#endif
      IF(NRS.NE.0) THEN

        if((NRS.eq.1).or.(nrs.eq.2).or.(nrs.eq.3) .or. (nrs.eq.5)) then  !saeed added 5
            IF(TimeLoc.GT.RSTIME2) THEN
               RSTIME1=RSTIME2
               RSTIME2=RSTIME2+RSTIMINC
               DO I=1,NP
                  RSNX1(I)=RSNX2(I)
                  RSNY1(I)=RSNY2(I)
               END DO
               IF(NRS.EQ.1) THEN
                 CALL RSGET(RSNX2,RSNY2)
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
               ENDIF
#endif
C... v50.xx sm -- nrs=5 for NUOPC coupling wave data from coupler
               NRS5: IF (NRS == 5) THEN    ! saeed

                 DO I=1,NP
                   RSX = RSNX1(I)
                   RSY = RSNY1(I)
                   RSNX1(I) = RSNX2(I)
                   RSNY1(I) = RSNY2(I)
                   RSNX2(I) = RSNX2(I) + (RSNX2(I)-RSX) ! extrapolate: 2*Rc - Rp
                   RSNY2(I) = RSNY2(I) + (RSNY2(I)-RSY)
                 ENDDO
               ENDIF NRS5 !(NRS = 5)

            ENDIF  !timeloc gt rstime2

            RStRatio=(TimeLoc-RSTIME1)/RSTIMINC
            
! tcm v56 2024 -- added option for limiting wave stress gradient magnitude
            IF (Limit_WaveStressGrad.eqv..false.) then            
               DO I=1,NP
                 RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I)-RSNX1(I)))
                 RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I)-RSNY1(I)))
                 WSX2(I) = WSX2(I) + RSX
                 WSY2(I) = WSY2(I) + RSY

Casey 090302: Added these lines for output to the rads.64 file.
!TCM v56 2024: removed the test for NRS=3, to allow for outputting if any NRS>0
!                IF(ABS(NRS).EQ.3) then
                 RSNXOUT(I) = RSX
                 RSNYOUT(I) = RSY
!                ENDIF
               ENDDO
            ELSE !test and limit wave stress
               DO I=1,NP
                 RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I)-RSNX1(I)))
                 RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I)-RSNY1(I)))
                 RSXYMAG = SQRT(RSX*RSX + RSY*RSY)
                 IF ((RSXYMAG .gt. WaveStressGrad_Cap ).and.
     &               (RSXYMAG.gt.0.0d0)) then
                    RSSCALE = WaveStressGrad_Cap/RSXYMAG
                    RSXLIM = RSX*RSSCALE
                    RSYLIM = RSY*RSSCALE
                 ELSE
                    RSXLIM = RSX
                    RSYLIM = RSY
                 ENDIF
                 WSX2(I) = WSX2(I) + RSXLIM
                 WSY2(I) = WSY2(I) + RSYLIM
                 RSNXOUT(I) = RSX    !keep unlimited value for writing
                 RSNYOUT(I) = RSY    !keep unlimited value for writing
               ENDDO
            ENDIF  ! test for limiting wave stress gradient magnitudes
            
         ENDIF  !nrs = 1,2,3, or 5

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

!tcm v56 2024: Added test for limiting wave stress gradient magnitudes
            IF (Limit_WaveStressGrad.eqv..false.) then
               DO I=1,NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)
                  WSX2(I) = WSX2(I) + RSX
                  WSY2(I) = WSY2(I) + RSY
                  !  tcm v50.75 added for use whenever nrs=3 or nrs=4
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               ENDDO
               
            ELSE !test and limit wave stress gradients
             
               DO I=1,NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)

                  RSXYMAG = SQRT(RSX*RSX + RSY*RSY)
                  IF ((RSXYMAG .gt. WaveStressGrad_Cap ).and.
     &                (RSXYMAG.gt.0.0d0)) then
                     RSSCALE = WaveStressGrad_Cap/RSXYMAG
                     RSXLIM = RSX*RSSCALE
                     RSYLIM = RSY*RSSCALE
                  ELSE
                     RSXLIM = RSX
                     RSYLIM = RSY
                  ENDIF

                  WSX2(I) = WSX2(I) + RSXLIM
                  WSY2(I) = WSY2(I) + RSYLIM
                  RSNXOUT(I) = RSX !keep unlimited value for writing
                  RSNYOUT(I) = RSY !keep unlimited value for writing
               ENDDO
            ENDIF !test for limit wave stress gradients      
         ENDIF !(NRS = 4)
      ENDIF  !(End test for updating wave radiation stress)

C     jgf48.4627 Skip past GWCE and momentum calculations if only
C     meteorological output was requested.
      IF (.NOT.METONLY) THEN

C   kmd48.33 - added in information for the elevation boundary conditions
C              used in prognostic runs. The diagnostic information is read
C              in once and used during the simulation. Note that there is
C              no ramp utilized for this boundary condition.
      IF((C3D).AND.(RES_BC_FLAG.GT.0).AND.(CBAROCLINIC)) THEN
         IF ((ABS(RES_BC_FLAG).GE.1).AND.(NOPE.GT.0)) THEN
Casey 140701: Added the following IF statement.
           IF(BCFLAG_LNM.EQ.1)THEN
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
Casey 140701: Ending the IF statement.
           ENDIF
         END IF
      END IF


! Apply the tide potential forcing subroutine      
      call applyTidePotentialForcing(TimeLoc, TimeH)

C...
C...  Depth Averaged Baroclinic Forcing needed by GWCE and 2DDI Momentum
C...  Compute this (divided by H, i.e., Bx/H, By/H) as a nodally averaged
C...  quantity for smoothing
C...
      IF (CBaroclinic) THEN

!     Kendra45.12 - Test placement of BPG
         IF (C3DVS) CALL BPG3D(IT)
! arash 01/20/2016

c     WJP 03.23.18 adding in the 2D baroclinic option
         IF (C2DDI) CALL BPG2D(IT,TimeLoc)

      ENDIF
C...
C...  COMPUTE SPECIFIED NORMAL FLOW BOUNDARY CONDITION
C...
      IF(NFLUXF.EQ.1) THEN
         IF (NFFR.GT.0) THEN
          DO J=1,NFFR
            IF(FPER(J).EQ.0.D0) THEN
               NCYC=0
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
               IF (LBCODEI(I).EQ.32) THEN
                  ARG=ARGJ-ENPH(J,I)
                  EN2(I)=EN2(I)+ENAM(J,I)*RFF*COS(ARG)
               ENDIF
            END DO
          END DO
         END IF

         IF((NFFR.EQ.0).OR.(NFFR.EQ.-1)) THEN
            IF(TimeLoc.GT.QTIME2) THEN
               QTIME1=QTIME2
               QTIME2=QTIME2+FTIMINC
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     QNIN1(J)=QNIN2(J)
                     READ(20,*) QNIN2(J)
                  ELSEIF(LBCODEI(J).EQ.32) THEN
                     QNIN1(J)=QNIN2(J)
                     ENIN1(J)=ENIN2(J)
                     READ(20,*) QNIN2(J),ENIN2(J)
                  ENDIF
               END DO
            ENDIF

            QTRATIO=(TimeLoc-QTIME1)/FTIMINC
            QN2=RampExtFlux*(QNIN1+QTRATIO*(QNIN2-QNIN1))
            EN2=RampExtFlux*(ENIN1+QTRATIO*(ENIN2-ENIN1))
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
               UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
               QN1(J)=H2(NNBB)*UN1
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
              UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
              QN1(J)=H2(NNBB)*UN1
              ENDIF
            END DO
          ENDIF
C...
C...  COMPUTE SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C.... EXTERNAL BARRIER BOUNDARY NODES
C...  COBELL - MOVED TO WEIR_BOUNDARY.F
        IF(NFLUXB.EQ.1) THEN
          DO I=1,NVEL
            SELECT CASE(LBCODEI(I))
            CASE(3,13,23)
                CALL COMPUTE_EXTERNAL_BOUNDARY_FLUX(I,TIMELOC,QN2(I))
            END SELECT 
          END DO
        ENDIF

      IF(NFLUXIB64.EQ.1) THEN
         ! Setting ISSUBMERGED64 flag
         I = 0                         
         DO K = 1, NBOU  
            SELECT CASE(LBCODEI(I+1))
               CASE(4,24,5,25)
                  I = I + NVELL(K)*2
               CASE(64)
                  DO J = 1,NVELL(K)
                     I = I + 1
                     CALL SET_SUBMERGED64_AT(I,J,K,TIMELOC)
                     ISSUBMERGED64P(I) = ISSUBMERGED64(I)
                  ENDDO
                  I = I + NVELL(K)
               CASE DEFAULT
                  I = I + NVELL(K)
            END SELECT    
         ENDDO
         ! Cancelling ISSUBMERGED64 flag if any of the adjacent nodes along the boundary is not submerged. sb 8/9/2024
         I = 0                         
         DO K = 1, NBOU  
            SELECT CASE(LBCODEI(I+1))
               CASE(4,24,5,25)
                  I = I + NVELL(K)*2
               CASE(64)
                  DO J = 1,NVELL(K)
                     I = I + 1
                     NNBB1 = NBV(I)
                     HTOT = ETA2(NNBB1) + IFNLFA*DP(NNBB1)
                     IF(ISSUBMERGED64(I).NE.0.AND.
     &                  NODECODE(NNBB1).NE.0.AND.
     &                  HTOT.LT.4.D0*H0) THEN
                        UNSUBMERGE = .FALSE.
                        X1 = X(NNBB1)
                        Y1 = Y(NNBB1)
                        ET1 = ETA2(NNBB1)
                        IF(J>1) THEN
                           NNBBNEI = NBV(I-1)
                           X2 = X(NNBBNEI)
                           Y2 = Y(NNBBNEI)
                           ET2 = ETA2(NNBBNEI)
                           LEN = SQRT((X1-X2)**2+(Y1-Y2)**2)
                           SLP = ABS(ET1-ET2)/LEN
                           IF(SLP.GT.0.0001) THEN
                              UNSUBMERGE = .TRUE.
                           ENDIF
                        ENDIF
                        IF(J<NVELL(K)) THEN
                           NNBBNEI = NBV(I+1)
                           X2 = X(NNBBNEI)
                           Y2 = Y(NNBBNEI)
                           ET2 = ETA2(NNBBNEI)
                           LEN = SQRT((X1-X2)**2+(Y1-Y2)**2)
                           SLP = ABS(ET1-ET2)/LEN
                           IF(SLP.GT.0.0001) THEN
                              UNSUBMERGE = .TRUE.
                           ENDIF
                        ENDIF
                        IF(.NOT.UNSUBMERGE) THEN
                           NWETNEI = 0
                           DO JJ = 2,NNeigh(NNBB1)
                              II = NeiTab(NNBB1,JJ)
                              IF(NODECODE(II).NE.0) THEN
                                 NWETNEI = NWETNEI + 1
                              ENDIF
                           ENDDO
                           IF(NWETNEI.LT.3) THEN
                              UNSUBMERGE = .TRUE.
                           ENDIF
                        ENDIF
                        IF(UNSUBMERGE) THEN
                           ISSUBMERGED64(I) = 0
                           NNBB2 = IBCONN(I)
                           IF(NNBB2>0) THEN
                              ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
                  I = I + NVELL(K)
               CASE DEFAULT
                  I = I + NVELL(K)
            END SELECT    
         ENDDO 
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
C...ZC - SIMPLIFIED THIS SECTION, MOVED FLUX COMPUTATION TO WEIR_BOUNDARY.F
C   USE THE COMPILER FLAGS TO CHANGE IMPLEMENTATION:
C        -DORIGWEIR - JOANNES WESTERINK ET AL IMPLEMENTATION (ORIGINAL)
C                     THIS DOES NOT APPEAR TO HAVE BEEN USED IN QUITE 
C                     SOME TIME
C         DEFAULT   - SHINTARO BUNYA/SEIZO TANAKA IMPLEMENTATION FOR 
C                     CHECKING FOR A WET EDGE BEFORE PASSING FLOW
C                     ACROSS A WEIR
      IF(NFLUXIB.EQ.1) THEN
         NIBNODECODE(:) = 0
         I = 0                         
         DO K = 1, NBOU  
            SELECT CASE(LBCODEI(I+1))
                CASE(4,24,5,25)
                    DO J = 1,NVELL(K)*2
                        I = I + 1
                        CALL COMPUTE_INTERNAL_BOUNDARY_FLUX(I,J,K,
     &                      TIMELOC,QN2(I))
                    ENDDO
                CASE(64)
                    DO J = 1,NVELL(K)*2
                        I = I + 1
                        IF(J.LE.NVELL(K)) THEN
                           ISFRONT = .TRUE.
                        ELSE
                           ISFRONT = .FALSE.
                        ENDIF
                        CALL COMPUTE_INTERNAL_BOUNDARY64_FLUX(I,J,K,
     &                      TIMELOC,QN2(I),ISFRONT)
                    ENDDO
                CASE DEFAULT
                    I = I + NVELL(K)
            END SELECT    
         ENDDO
      ENDIF

C...
C...  COMPUTE INWARD/OUTWARD NORMAL FLOW FOR INTERNAL BARRIER
C.... BOUNDARY NODES THROUGH CROSS BARRIER PIPES
C.... NOTE THAT THIS ADDS AN ADDITIONAL FLOW COMPONENT INTO QN2
C...
      IF(NFLUXIBP.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
                PIPE_FLUX = 0D0
                CALL COMPUTE_CROSS_BARRIER_PIPE_FLUX(I,TIMELOC,PIPE_FLUX)
                QN2(I) = QN2(I) + PIPE_FLUX
            ENDIF
         ENDDO
      ENDIF 

      if(subdomainOn.and.enforceBN.eq.1) call readFort019(it)  ! NCSU Subdomain
      if(subdomainOn.and.enforceBN.eq.2) call readFort020(it)  ! NCSU Subdomain
      if(subdomainOn.and.enforceBN.eq.2) call readFort021(it)  ! NCSU Subdomain


!----------------------SOLVE GWCE---------------------------------------
C...
C...  Compute new water surface elevation 
      call solveGWCE(it, ITIME_BGN, timeloc, timeh)

!---------------------END SOLVE GWCE------------------------------------

!JLW: adding subgrid
      IF(level0)THEN
#ifdef ADCNETCDF
         !JLW: update vertex arrays
         call getVertLookup()
         !JLW: update nodecode, ncele, etc
         call computeWettingAndDrying(it)
#else
         call allMessage(ERROR, "Subgrid calculation requires building"
     &                          //"with netcdf enabled")
         call terminate()
#endif         
      ELSE
!--------------------WETTING AND DRYING--------------------------------

         call computeWettingAndDrying(it)

!--------------------WETTING AND DRYING--------------------------------
      ENDIF



!.....DMW202401 Add update for flow depth at current timestep
!......We do this after wetting/drying since bathy could change with subgrid barriers.
      DO I = 1,NP
         H2(I) = DP(I)+IFNLFA*ETA2(I)
      ENDDO
C...   
C---------------2DDI Momentum Equation Solution-------------------------
C...
      call solveMomentumEq()
C...
C-------------END OF 2DDI MOMENTUM EQUATION SOLUTION--------------------
C...

C.....DW - ABS SPONGE LAYER     
      IF ( C2DDI .AND. LoadAbsLayerSigma .AND. (sponge_dis_mthd == 0) ) THEN
         CALL SPONGE_OPSPLIT0( IT, DTDP/2, TimeLoc, TimeH)  
      END IF

      IF ( LoadAbsLayerSigma ) THEN
         CALL SPONGE_SHIFT_SOLN( ) 
      END IF
C     
C.....END DW
C     

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
      ENDIF
C...
C...  Collect maximum values of variables   v46.50 sb 11/11/2006
      call collectMinMaxData(timeloc)
      !
      ! jgf52.08.01: If detailed inundation output was requested, 
      ! collect detailed inundation data. 
      if (inundationOutput.eqv..true.) then
         call collectInundationData(timeloc, it)
      endif
C
C...
C...  WRITE OUTPUT
C...
C
      ! NCSU Subdomain Modeling   
      if (subdomainOn.and.NOUTGS.eq.1) call writeFort065(it) 
      if (subdomainOn.and.NOUTGS.eq.2) call writeFort066(it) 
      if (subdomainOn.and.NOUTGS.eq.2) call writeFort067(it) 
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

      if ((ABS(NHSTAR) > 0 .and. ITEST == IT) .or. (-IHOT == IT)) then  !tcm v51.26 added abs(nhstar) to handle when nhstar = -1
         if( (MNWPROH > 0) ) then   !Writer for HSfile
            if( (.not.C3D) ) then                     !Writer for HSfile
#ifdef CMPI
               CALL writeHotstart_through_hswriter(TimeLoc,IT)  !st3 hsfile
#endif
            else
               write(6,*) 'HS writer does not support C3D'
            endif
         else
            if (NHOUTONCE) then
               ! Make sure never enters again
               NHSINC = huge(NHSINC)
            endif
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

C...  SAVE THE CURRENT TIME LEVEL OF BARRIER INTO THE PREVIOUS TIME LEVEL
      IF(NFLUXIB.EQ.1)THEN
        BARINHT1(:)=BARINHT2(:)
      ENDIF
      IF(NFLUXB.EQ.1)THEN
        BARLANHT1(:)=BARLANHT2(:)
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

         runperccomplete=REAL(IT-ITHS)/REAL(NT-ITHS)*100d0    !percentage of time steps completed
#ifdef CMPI
         IF(MYPROC.EQ.0.AND.ELMAX.LT.WarnElev.AND.ITEST.EQ.IT) THEN
            IF (KEMAX.GT.0) THEN
               WRITE(ScreenUnit,1991)
     &            IT,runperccomplete,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,
     &            VELMAX,KVMAX,MYPROC
            ELSE
               WRITE(ScreenUnit,1991)
     &            IT,runperccomplete,NUMITR,TimeLoc,0.,KEMAX,VELMAX,
     &            KVMAX,MYPROC

            ENDIF
 1991       FORMAT(1X,'TIME STEP =',I8,1x,F6.2,'% COMPLETE',
     &           5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,
     &           2X,'ON MYPROC = ',I4)
         ENDIF
         WarnElevExceeded = 0
         IF(ELMAX.GT.WarnElev) THEN
C     WJP 03.27.2018 with the addition of the hash table we now need to
C     output the global node number to debug the local PE fort.14s
            WRITE(ScreenUnit,1993)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),NODES_LG(KEMAX),
     &           VELMAX,NODES_LG(KVMAX),MYPROC
            WRITE(16,1993) IT,NUMITR,TimeLoc,ETA2(KEMAX),
     &           NODES_LG(KEMAX),VELMAX,NODES_LG(KVMAX),MYPROC
 1993       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,
     &           2X,'ON MYPROC = ',I4,
     &           3X,'** WARNING: Elevation.gt.WarnElev **')
            IF (WarnElevDump) WarnElevExceeded=1
         ENDIF

         IF(WarnElevDump)THEN
             call WarnElevSum(WarnElevExceeded)
             IF (WarnElevExceeded.ne.0) THEN
                CALL WriteWarnElev(TimeLoc, IT)
             ENDIF
         ENDIF
         ErrorElevExceeded = 0                ! Clint's Zombie Slyaer
         IF(ELMAX.GT.ErrorElev) THEN

            WRITE(ScreenUnit,1995)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),NODES_LG(KEMAX),
     &           VELMAX,NODES_LG(KVMAX),MYPROC
            WRITE(16,1995) IT,NUMITR,TimeLoc,ETA2(KEMAX),
     &           NODES_LG(KEMAX),VELMAX,NODES_LG(KVMAX),MYPROC
 1995       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,
     &           2X,'ON MYPROC = ',I4,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev,',
     &           ' ADCIRC stopping. **')
         ErrorElevExceeded = 1                ! Clint's Zombie Slayer
         ENDIF
         call WarnElevSum(ErrorElevExceeded)  ! Clint's Zombie Slayer 2010.08.07
         IF( ErrorElevExceeded /= 0 ) THEN    !  Finalize MPI Environment,
           Flag_ElevError = .true. 
           CALL MSG_FINI()                    !  if there are Error Elvation nodes.
           CALL EXIT(1)                             !
         ENDIF                                !st3

#ifdef VELCHECK
         WarnVelExceeded = 0                
         IF (IT .EQ. ITIME_BGN) THEN
            WRITE(16,*) 'dmw202401 THIS CODE IS NOW MAKING ITS FIRST CHECK FOR WARNING VELS'
         ENDIF
         IF(VELMAX.GT.WarnVel) THEN
         
            WRITE(ScreenUnit,3995)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),NODES_LG(KEMAX),
     &           VELMAX,NODES_LG(KVMAX),MYPROC
            WRITE(16,3995) IT,NUMITR,TimeLoc,ETA2(KEMAX),
     &           NODES_LG(KEMAX),VELMAX,NODES_LG(KVMAX),MYPROC
 3995       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,
     &           2X,'ON MYPROC = ',I4,/,
     &           2X,'** WARNING: Velocity.gt.WarnVel **')
            WarnVelExceeded = 1                
         ENDIF   
#endif
		 
#else
         IF(ELMAX.LT.WarnElev.AND.ITEST.EQ.IT) THEN
            IF (KEMAX.GT.0) THEN
               WRITE(ScreenUnit,1992)
     &              IT,runperccomplete,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            ELSE
               WRITE(ScreenUnit,1992)
     &              IT,runperccomplete,NUMITR,TimeLoc,0.,KEMAX,VELMAX,KVMAX
            ENDIF
 1992       FORMAT(1X,'TIME STEP =',I8,1x,F6.2,'% COMPLETE',
     &           5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9)
         ENDIF
         IF(ELMAX.GT.WarnElev) THEN
            WRITE(ScreenUnit,1994)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1994) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,
     &           2X,'** WARNING: Elevation.gt.WarnElev **')
            IF (WarnElevDump) CALL WriteWarnElev(TimeLoc, IT)
         ENDIF
         IF(ELMAX.GT.ErrorElev) THEN
            write(6,*) 'elmax, errorelev',elmax, errorelev !jgf
            Flag_ElevError = .true.
            WRITE(ScreenUnit,1996)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1996) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1996       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev, '
     &           'ADCIRC stopping. **')
            CALL terminate()
         ENDIF
#ifdef VELCHECK
         IF (IT .EQ. ITIME_BGN) THEN
            WRITE(16,*) 'dmw202401 THIS CODE IS NOW MAKING ITS FIRST CHECK FOR ERROR VELS'
         ENDIF
         IF(VELMAX.GT.WarnVel) THEN
            WRITE(ScreenUnit,3996)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,3996) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 3996       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,/,
     &           2X,'** WARNING: Velocity.gt.WarnVel **')
         ENDIF
#endif

#endif

      ENDIF

#ifdef VELCHECK
!....DMW 20240109 - WARNING FOR NaNs IN VELOCITY SOLN
      foundNaN = .false.
      DO I = 1,NP
         IF ((I.EQ.1).AND.(IT.EQ.ITIME_BGN)) THEN
            WRITE(16,*) 'dmw202401 THIS CODE IS NOW MAKING ITS FIRST CHECK FOR NAN VELS'
         ENDIF
         IF(ISNAN(UU2(I)).OR.ISNAN(VV2(I))) THEN
            foundNaN = .true.
            IF(NSCREEN.NE.0) THEN
               WRITE(ScreenUnit,2993)
     &              IT,NUMITR,TimeLoc,NODES_LG(I),
     &              ETA2(I),DP(I),MYPROC
               WRITE(16,2993) IT,NUMITR,TimeLoc,
     &              NODES_LG(I),ETA2(I),DP(I),MYPROC
 2993          FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &              5X,'TIME = ',E15.8,
     &              2X,'SPEED IS NaN AT NODE I =',I9,
     &              /,2X,'ETA2(I) = ', 1pE12.4E3,
     &              2X,'DP(I) = ', 1pE12.4E3,
     &              2X,'ON MYPROC = ',I4)
               DO J = 2,NNeigh(I)
                  JJ = NeiTab(I,J)
                  WRITE(ScreenUnit,2995) NODES_LG(JJ),ETA2(JJ),DP(JJ),               
     &                 UU2(JJ),VV2(JJ),MYPROC
              
                  WRITE(16,2995) NODES_LG(JJ),ETA2(JJ),DP(JJ),
     &                 UU2(JJ),VV2(JJ),MYPROC
 2995             FORMAT(1X,'AT NEIGHBOR NODE I =',I9,
     &                 /,2X,'ETA2(I) = ', 1pE12.4E3,
     &                 2X,'DP(I) = ', 1pE12.4E3,
     &                 2X,'UU2(I) = ', 1pE12.4E3,
     &                 2X,'VV2(I) = ', 1pE12.4E3,
     &                 2X,'ON MYPROC = ',I4)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      IF (foundNaN) THEN
         IF(NSCREEN.NE.0) THEN
            WRITE(ScreenUnit,2994)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),NODES_LG(KEMAX),
     &           VELMAX,NODES_LG(KVMAX),MYPROC
            WRITE(16,2994) IT,NUMITR,TimeLoc,ETA2(KEMAX),
     &           NODES_LG(KEMAX),VELMAX,NODES_LG(KVMAX),MYPROC
 2994       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I9,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I9,
     &           2X,'ON MYPROC = ',I4,
     &           3X,'** WARNING: Velocity is NaN **')
         ENDIF
      ENDIF
      IF ((foundNaN).AND.(IT_foundNaN.EQ.0)) THEN
         IT_foundNaN = IT
      ELSEIF ((IT.EQ.IT_foundNaN + 20).AND.(IT_foundNaN.NE.0)) THEN
         earlyterminate = 1
      ENDIF
#ifdef CMPI
      call EarlyTermSum( earlyterminate )
#endif

#endif


C...
C...  ****************** TIME STEPPING LOOP ENDS HERE ********************
C...

c. RJW merged 08/26/2008 Casey 071219: Added the folowing call to the mass balance subroutine.
C      IF(C3DVS)THEN
C         CALL MASSBAL3D(IT)
C      ENDIF

      IF ( IT == NT ) THEN
        ! close NWS14 NetCdf files !
        CALL CLOSE_MET_FILES() ; 
      END IF        

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

C******************************************************************************

      END SUBROUTINE TIMESTEP
C******************************************************************************


!****************************************************************************************
! Subroutine to compute tide potential at the nodes                                     
!                                                                                       
! Tidal Potential Forcing
! Note, the Earth tide potential reduction factor, ETRF(J) has been
!       incorporated into this calculation.
!
! WJP 02.21.2018 Rewriting for general formula for all species,
! cleaner code, and correct identification of 0,1,2 or higher species
!****************************************************************************************
      subroutine applyTidePotentialForcing(TimeLoc, TimeH)
         use mod_tidepotential, only: tidePotential
         use global, only: NTIP, TIP1, TIP2, CTIP, NTIF, RampTip, PERT, FACET,
     &                     AMIGT, FFT, TPK, ETRF, SALTPHA, SALTAMP, L_N
         use mesh , only: np, SLAM
         implicit none
         real(8), intent(IN) :: TimeLoc, TimeH
         real(8) :: ArgT, ArgTP, ArgSAlt, SALTMUL, TPMUL
         integer :: i, j, na, ncyc

         IF(CTIP) THEN

            ! use the full formula 
            IF ( tidePotential%active() ) THEN
               TIP2 = tidePotential%compute(TimeLoc, NP, SLAM) * RampTip
            ENDIF  

            ! The traditional approach: a sum of different constituents
            ! TIP & SAL - for tidePotential%active() = .TRUE.  
            ! SAL       - for tidePotential%active() = .FALSE.
            IF(NTIP.EQ.2.or..not.tidePotential%active())THEN 
               DO J=1,NTIF
                  IF (PERT(J).EQ.0D0) THEN
                     NCYC=0
                  ELSE
                     NCYC=INT(timeh/PERT(J))
                  ENDIF

                  ARGT=AMIGT(J)*(timeh-NCYC*PERT(J))+FACET(J)
                  TPMUL=RampTip*ETRF(J)*TPK(J)*FFT(J)
                  SALTMUL=RampTip*FFT(J)

   !              WJP: We actually want to compare against diurnal and get 0,1,2
   !              or higher species (was opposite beforehand)
                  NA = MIN(NINT(AMIGT(J)/7d-5),2)

   !              WJP: Rewritten so that we have a general formula for all
   !              species
                  IF(tidePotential%active())THEN
                     DO I = 1, NP
                        ARGSALT = ARGT - SALTPHA(J,I) ; 
                        TIP2(I) = TIP2(I) + SALTMUL * SALTAMP(J,I) * COS(ARGSALT)
                     END DO
                  ELSE
                     DO I = 1,NP
                        ARGTP   = ARGT + NA*SLAM(I)
                        ARGSALT = ARGT - SALTPHA(J,I)
                        TIP2(I) = TIP2(I) + TPMUL * L_N(NA,I) * COS(ARGTP)
     &                            + SALTMUL * SALTAMP(J,I) * COS(ARGSALT)
                     ENDDO
                  ENDIF

               ENDDO
            ENDIF   
         ENDIF
      end subroutine applyTidePotentialForcing


C****************************************************************************************
C   Subroutine to compute Scalar Transport                                              *
C                                                                                       *
C   Note, this is not set up for parallel operation.                                    *
C****************************************************************************************
C
      SUBROUTINE SCALAR_TRANS_2D (IT,TimeLoc)

      USE GLOBAL, only: ch1, trans_lv_b, nodecode, noff, lumpt, dtdp,
     &    trans_lv_a, soursin, uu1, vv1, ifnlfa, eta1, bedstr, tk,
     &    setmessagesource, unsetmessagesource, DEBUG,
     &    allMessage, H1

      USE ADC_CONSTANTS, ONLY: rhowat0
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, SFAC
      USE NodalAttributes, ONLY: EVC

      IMPLICIT NONE

      INTEGER IE, I                         !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI

      REAL(8) C1, C2, C3, CBEDSTRD, CBEDSTRE, CCRITD
      REAL(8) CH1N1, CH1N2, CH1N3, CHSUM
      REAL(8) DHDX, DHDY
      REAL(8) DXXYY11, DXXYY12, DXXYY13
      REAL(8) DXXYY22, DXXYY23
      REAL(8) DXXYY33
      REAL(8) ECONST
      REAL(8) EVC1, EVC2, EVC3, EVCEA
      REAL(8) FDDDODT, FDDODODT
      REAL(8) HEA
      REAL(8) H1N1, H1N2, H1N3
      REAL(8) HSD, HSE
      REAL(8) SFacAvg
      REAL(8) SS1N1, SS1N2, SS1N3
      REAL(8) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(8) TEMP_LV_B1, TEMP_LV_B2, TEMP_LV_B3
      REAL(8) U1N1, U1N2, U1N3
      REAL(8) V1N1, V1N2, V1N3
      REAL(8) UV1
      REAL(8) UEA, VEA, UPEA, VPEA
      REAL(8) WS, WSMOD

      REAL(8) AreaIE2
      REAL(8) FDDD, FDDOD
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) DDX1,DDX2,DDX3,DDY1,DDY2,DDY3
      REAL(8) DXX11,DXX12,DXX13,DXX21,DXX22,DXX23,DXX31,DXX32,DXX33
      REAL(8) DYY11,DYY12,DYY13,DYY21,DYY22,DYY23,DYY31,DYY32,DYY33
      REAL(8) TimeLoc

      call setMessageSource("scalar_trans_2D")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
        BEDSTR=H1(I)*UV1*TK(I)*RhoWat0                           ![N/m^2]
        C1=CH1(I)/H1(I)

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
         H1N1=H1(NM1)
         H1N2=H1(NM2)
         H1N3=H1(NM3)
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

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE SCALAR_TRANS_2D
C**********************************************************************

C****************************************************************************************
C   Subroutine to compute 2D SigmaT fields from 3D salinity and/or temperature fields   *
C                                                                                       *
C                                    R.L.  6/22/05                                      *
C****************************************************************************************
C

      SUBROUTINE CALC_SIGMAT_2D ()

      USE GLOBAL

      INTEGER NH

      call setMessageSource("calc_sigmat_2D")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


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
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
c*********************************************************************
      END SUBROUTINE CALC_SIGMAT_2D
c*********************************************************************

c******************************************************************************
c  SUBROUTINE BPG2D                                                           *
c                                                                             *
c  WJP: FINDING VIDBCPDX(DY) FROM EITHER THE DEPTH-AVERAGED EQUATION OR       *
C  FROM THE FULL 3D INTEGRAL USING INFORMATION FROM A 3D BAROCLINIC MODEL.    *
C  IN THE LATTER CASE MOMENTUM DISPERSION IS ALSO ESTIMATED                   *
c                                                                             *
c******************************************************************************
      
      SUBROUTINE BPG2D(IT,TimeLoc)

      USE SIZES, ONLY : myproc
      USE ADC_CONSTANTS, ONLY: G, RhoWat0, SigT0
      USE GLOBAL, ONLY: IFNLFA, IFNLCT, IDEN, RAMP, NODECODE, NOFF, 
     &    VIDBCPDXOH, VIDBCPDYOH, ETA2, DTDP, DASigT, 
     &    IFSFM, H2
      USE MESH, ONLY : NE, NM, FDXE, FDYE, NP, DP, AREAS, NODELE, 
     &      NEITABELE, NEITAB, NNEIGH, TOTALAREA,
     &      SFacEle, SFMYEle, SFMXEle
      USE Couple2BC3D, ONLY : Update_BC3D_Info, FBPG_Disp_from_BC3D
      IMPLICIT NONE
    
      INTEGER,INTENT(IN) :: IT 
      REAL(8),INTENT(IN) :: TimeLoc 
      INTEGER :: IE             !element loop counter
      INTEGER :: J              !node loop counter
      INTEGER :: K              !vertical node loop counter (1-top, BCK-last z contour)
      INTEGER :: NCELE          !element code
      INTEGER :: NEle           !local value of NetTabEle
      INTEGER :: NM1,NM2,NM3    !local node numbers used to compute gradients
      INTEGER :: NC1,NC2,NC3    !local node codes

      REAL(8) :: SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac 
      REAL(8) :: AreaIE2, AreaEle
      REAL(8) :: H2N1, H2N2, H2N3
      REAL(8) :: FDX1, FDX2, FDX3
      REAL(8) :: FDY1, FDY2, FDY3
      REAL(8) :: FDX1O2A, FDX2O2A, FDX3O2A
      REAL(8) :: FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) :: EtaN1, EtaN2, EtaN3
      REAL(8) :: DEta2DX, DEta2DY
      REAL(8) :: DRhoDX, DRhoDY
      REAL(8) :: DARhoMRho0N1,DARhoMRho0N2,DARhoMRho0N3
      REAL(8) :: VIDBCPDXOHN1,VIDBCPDXOHN2,VIDBCPDXOHN3
      REAL(8) :: VIDBCPDYOHN1,VIDBCPDYOHN2,VIDBCPDYOHN3
      REAL(8) :: VIDBCPDXOHAvgArea,VIDBCPDYOHAvgArea

      ! For the Coupling to 3D Baroclinic mode 
      IF (abs(IDEN).ge.5) THEN
         ! Update the 3D BC Information for the current timestep 
         call Update_BC3D_Info(IDEN,TimeLoc) 

         if (abs(IDEN).ne.7.AND.abs(IDEN).ne.9) then
            ! Calculating the full Baroclinic pressure gradient and
            ! momentum dispersion terms required at every timestep
            ! from the 3D BC information that is on a longer time scale 
            call FBPG_Disp_from_BC3D()
         endif 
         RETURN
      ENDIF  
     
      ! For the depth-averaged baroclinic mode without coupling
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         H2N1=H2(NM1)
         H2N2=H2(NM2)
         H2N3=H2(NM3)
         EtaN1=IFNLFA*Eta2(NM1)
         EtaN2=IFNLFA*Eta2(NM2)
         EtaN3=IFNLFA*Eta2(NM3)
         SFacAvg= SFacEle(IE)
C..... BEG DW/WJP
         SFmxAvg= SFMXEle(IE) ; 
         SFmyAvg= SFMYEle(IE) ; 
         sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ;
         sfdyfac = (1 - IFSFM)*1.0D0 + IFSFM*SFmyAvg ;
C..... END DW/WJP

         AreaIE2=Areas(IE)
         AreaEle=NCEle*AreaIE2*0.5d0
         FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*MX !b1*mx
         FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*MX !b2*mx
         FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*MX !b3*mx
         FDY1 = FDYE(1,IE)*sfdyfac ; !c FDY1=(X(NM3)-X(NM2))*MY !a1*my
         FDY2 = FDYE(2,IE)*sfdyfac ; !c FDY2=(X(NM1)-X(NM3))*MY !a2*my
         FDY3 = FDYE(3,IE)*sfdyfac ; !c FDY3=(X(NM2)-X(NM1))*MY !a3*my         
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
     &         +DARhoMRho0N3*FDX3O2A
         DRhoDY=DARhoMRho0N1*FDY1O2A+DARhoMRho0N2*FDY2O2A
     &         +DARhoMRho0N3*FDY3O2A
         VIDBCPDXOHN1=Ramp*G*
     &        (DARhoMRho0N1*DEta2DX+0.5d0*H2N1*DRhoDX)
         VIDBCPDXOHN2=Ramp*G*
     &        (DARhoMRho0N2*DEta2DX+0.5d0*H2N2*DRhoDX)
         VIDBCPDXOHN3=Ramp*G*
     &        (DARhoMRho0N3*DEta2DX+0.5d0*H2N3*DRhoDX)
         VIDBCPDYOHN1=Ramp*G*
     &        (DARhoMRho0N1*DEta2DY+0.5d0*H2N1*DRhoDY)
         VIDBCPDYOHN2=Ramp*G*
     &        (DARhoMRho0N2*DEta2DY+0.5d0*H2N2*DRhoDY)
         VIDBCPDYOHN3=Ramp*G*
     &        (DARhoMRho0N3*DEta2DY+0.5d0*H2N3*DRhoDY)
         VIDBCPDXOHAvgArea=AreaEle*(VIDBCPDXOHN1+VIDBCPDXOHN2
     &        +VIDBCPDXOHN3)/3.D0
         VIDBCPDYOHAvgArea=AreaEle*(VIDBCPDYOHN1+VIDBCPDYOHN2
     &        +VIDBCPDYOHN3)/3.D0
         VIDBCPDXOH(NM1)=VIDBCPDXOH(NM1)+VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM2)=VIDBCPDXOH(NM2)+VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM3)=VIDBCPDXOH(NM3)+VIDBCPDXOHAvgArea
         VIDBCPDYOH(NM1)=VIDBCPDYOH(NM1)+VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM2)=VIDBCPDYOH(NM2)+VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM3)=VIDBCPDYOH(NM3)+VIDBCPDYOHAvgArea
      ENDDO

      DO J = 1,NP
         IF (TotalArea(J).NE.0.d0) THEN
            VIDBCPDXOH(J) = VIDBCPDXOH(J)/TotalArea(J)
            VIDBCPDYOH(J) = VIDBCPDYOH(J)/TotalArea(J)
         ENDIF
      ENDDO

C***********************************************************************
      END SUBROUTINE BPG2D
C***********************************************************************

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

      SUBROUTINE BPG3D(IT)

Casey: Added the following variable declarations from GLOBAL.
C
      USE GLOBAL, ONLY: IFNLFA, IFNLCT,
     &      NODECODE, NOFF, VIDBCPDXOH, VIDBCPDYOH,
     &      ETA2, H0, DTDP, H2
      use mod_logging, only: setMessageSource, unsetMessageSource, allMessage,
     &      logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      USE ADC_CONSTANTS, ONLY: SigT0
      USE MESH, ONLY : NM, X, Y, NP, DP, AREAS, NODELE, 
     &      NEITABELE, NEITAB, NNEIGH, SFAC
      USE GLOBAL_3DVS, ONLY : BCP, BPG, NFEN, SIGT, SIGMA, IDEN,
     &      GORhoOAMB, AMB, B, IY
! arash: this is not necessary. I just want to print this here
     &     ,Sal

#ifdef CMPI
      USE MESSENGER
      IMPLICIT NONE
      REAL(8) :: DUMV1(1),DUMV2(1)
#else
      IMPLICIT NONE
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

Casey 140701: Debug.
!     REAL(8) :: Hs             !Total water depth at time level s
!     REAL(8) :: HsOAMB         !Hs/(a-b)
!     REAL(8) :: HsHsOAMBAMB    !(Hs/(a-b))^2
!     REAL(8) :: Hsp1           !Total water depth at time level s+1
!     REAL(8) :: Hsp1OAMB       !Hsp1/(a-b)
!     REAL(8) :: Hsp1Hsp1OAMBAMB !(Hsp1/(a-b))^2
      REAL(8) :: Hs             !Total water depth at time level s
      REAL(8) :: HsOAMB         !Hs/(a-b)
      REAL(8) :: HsHsOAMBAMB    !(Hs/(a-b))^2
      REAL(8) :: Hsp1           !Total water depth at time level s+1
      REAL(8) :: Hsp1OAMB       !Hsp1/(a-b)
      REAL(8) :: Hsp1Hsp1OAMBAMB !(Hsp1/(a-b))^2

      REAL(8) :: Zk            !z depth of any node k in the vertical
      REAL(8) :: DelSig        ! sigma(k+1)-sigma(k)
      REAL(8) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
      REAL(8) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
      REAL(8) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
      REAL(8) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
      REAL(8) :: SigmaNN       !Sigma value of a neighbor node
      REAL(8) :: HsN2          !Depth value of a neighbor node

      REAL(8) :: SFacAvg       ! kmd48.33bc add in spherical factors

      REAL(8) :: BCPN1,BCPN2,BCPN3,BCPNFirst !nodal values of BCP
      REAL(8) :: BCPDX2A,BCPDY2A !(Horiz. grads of BCP)*2*Element Area
      REAL(8) :: SigTAvg       !avg SigT between 2 vertical nodes
      REAL(8) :: HGORhoOAMB    !depth*gravity/(reference density)/(a-b)

      REAL(8) :: a1,a2,a3,b1,b2,b3
      REAL(8) :: TotalBCPGArea2

      REAL(8) :: DBCPDX2A
      REAL(8) :: DBCPDY2A
      COMPLEX(8) :: BCPG(NFEN)   !baroclinic pressure gradient
      COMPLEX(8) :: VIBCPG         !baroclinic pressure gradient

! af:
      Integer      :: IT
      Real(8)     :: time_local
      Real(8)     ::  day_local
      Real(8)     ::   diff_tol = 1.0d-5

      Real(8)     :: SigT_mid_point, SigmaMP, f1, f2, f3
      COMPLEX(8)  :: BCPG_mid_point

C
      call setMessageSource("bpg3D")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
            HGORhoOAMB=GORhoOAMB*H2(NH) !(gravity/rho ref)*(H/(a-b))
            BCP(NH,NFEN)=0.d0
! arash: this integration can be improved
            DO k=NFEN-1,1,-1    !loop over vertical nodes, starting at top and working down
               SigTAvg=(SigT(NH,k+1)+SigT(NH,k))/2.d0
               DelSig=Sigma(k+1)-Sigma(k)
               BCP(NH,k)=BCP(NH,k+1)+HGORhoOAMB*(SigTAvg-SigT0)*DelSig
            ENDDO
!--------------------------------------------------------------------------------------------------
! arash 21/1/2016: Simpson integration, using quadratic function interpolation.
! arash 25/2/2016: fixed a bug.
!------------------------------------------------------------------------------
            ! 1. initialize
            BCP (NH,NFEN) = 0.d0

            ! 2. to compute the integral at node NFEN-1, we need the mid-point value, which we evaluate by using a quadratic
            !    interpolation, using nodes NFEN, NFEN-1, and NFEN-2.
            SigmaMP = ( Sigma(NFEN) + Sigma(NFEN-1) ) / 2.0d0

            f1 = ( (SigmaMP       - Sigma(NFEN-1)) * (SigmaMP       - Sigma(NFEN-2)) )
     &         / ( (Sigma(NFEN)   - Sigma(NFEN-1)) * (Sigma(NFEN)   - Sigma(NFEN-2)) )
            f2 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-2)) ) 
     &         / ( (Sigma(NFEN-1) - Sigma(NFEN))   * (Sigma(NFEN-1) - Sigma(NFEN-2)) )
            f3 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-1)) ) 
     &         / ( (Sigma(NFEN-2) - Sigma(NFEN))   * (Sigma(NFEN-2) - Sigma(NFEN-1)) )

            SigT_mid_point = f1 * SigT(NH,NFEN) + f2 * SigT(NH,NFEN-1) + f3 * SigT(NH,NFEN-2)

            ! 3. evaluate the integral at node NFEN-1
            BCP (NH,NFEN-1) = (1.0/6.0) * ( Sigma(NFEN) - Sigma(NFEN-1) ) * HGORhoOAMB *
     &                        ( (SigT(NH,NFEN)-SigT0) + 4.0d0 * (SigT_mid_point-SigT0) +  (SigT(NH,NFEN-1)-SigT0) )

            ! 4. evaluate the integral at the rest of the vertical nodes
            Do k = NFEN-2, 1, -1

                ! 4.1. 
               SigmaMP = ( Sigma(k+1) + Sigma(k) ) / 2.0d0

               f1 = ( (SigmaMP    - Sigma(k+1)) * (SigmaMP    - Sigma(k)) )
     &            / ( (Sigma(k+2) - Sigma(k+1)) * (Sigma(k+2) - Sigma(k)) )
               f2 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k)) ) 
     &            / ( (Sigma(k+1) - Sigma(k+2)) * (Sigma(k+1) - Sigma(k)) )
               f3 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k+1)) ) 
     &            / ( (Sigma(k)   - Sigma(k+2)) * (Sigma(k)   - Sigma(k+1)) )

               SigT_mid_point = f1 * SigT(NH,k+2) + f2 * SigT(NH,k+1) + f3 * SigT(NH,k)

                ! 4.2. evaluate the integral
               BCP (NH,k) = BCP (NH,k+1) + (1.0/6.0) * ( Sigma(k+1) - Sigma(k) ) * HGORhoOAMB 
     &                    * ( (SigT(NH,k+1)-SigT0) + 4.0d0 * (SigT_mid_point-SigT0) +  (SigT(NH,k)-SigT0) )

            End Do

            ! 5. let's hope we are better of now.

!--------------------------------------------------------------------------------------------------
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

         HsOAMB=H2(NH)/AMB

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

               SigmaNN=B+AMB*(Zk+DP(N2))/H2(N2) !equivalent sigma value at neighbor
               CALL ZSURFBUOY_p3(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
!!!!!!!!!         CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
               NNFirst=N2       !save these values until end
               BCPNFirst=BCPN2  !save these values until end

               DO N=3,NNeigh(NH) !operate on rest of neighbors
                  N3=N2         !shift previously computed values
                  BCPN3=BCPN2   !shift previously computed values
                  N2=NeiTab(NH,N) !select new neighbor to work on

Casey: Changed "NolIFA" to "IFNLFA."

                  SigmaNN=B+AMB*(Zk+DP(N2))/H2(N2) !equivalent sigma value at neighbor
!!!!!!!!!            CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
                  CALL ZSURFBUOY_p3(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
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

            DO k=1,NFEN
Casey 150128: Debug.
               BPG(NH,k) = BCPG(k)
C              BPG(NH,k) = (0.D0,0.D0)
            END DO

c     Finished computing baroclinic pressure gradient (computed in level
c     coordinates) at each node in the vertical


c     Compute vertically integrated baroclinic pressure gradient for use
c     in the wave equation.  NOTE: For a prognostic model in which the
c     density field evolves in time, this calculation should be done
c     after the new density field is computed.  In this case one would
c     integrate over the vertical first and differentiate second.

! this old scheme has been replaced by the new one below ...
            VIBCPG=(0.d0,0.d0)
            DO k=NFEN-1,1,-1
               VIBCPG=VIBCPG+0.5d0*(BCPG(k+1)+BCPG(k))
     &              *(Sigma(k+1)-Sigma(k))
            ENDDO
!--------------------------------------------------------------------------------------------------
! arash 25/2/2016: Simpson integration, using quadratic function interpolation.
!------------------------------------------------------------------------------
            ! 1. initialize
            VIBCPG = ( 0.d0 , 0.d0 )

            ! 2. to compute the integral at node NFEN-1, we need the mid-point value, which we evaluate by using a quadratic
            !    interpolation, using nodes NFEN, NFEN-1, and NFEN-2.
            SigmaMP = ( Sigma(NFEN) + Sigma(NFEN-1) ) / 2.0d0

            f1 = ( (SigmaMP       - Sigma(NFEN-1)) * (SigmaMP       - Sigma(NFEN-2)) )
     &         / ( (Sigma(NFEN)   - Sigma(NFEN-1)) * (Sigma(NFEN)   - Sigma(NFEN-2)) )
            f2 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-2)) ) 
     &         / ( (Sigma(NFEN-1) - Sigma(NFEN))   * (Sigma(NFEN-1) - Sigma(NFEN-2)) )
            f3 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-1)) ) 
     &         / ( (Sigma(NFEN-2) - Sigma(NFEN))   * (Sigma(NFEN-2) - Sigma(NFEN-1)) )

            BCPG_mid_point = f1 * BCPG(NFEN) + f2 * BCPG(NFEN-1) + f3 * BCPG(NFEN-2)

            ! 3. evaluate the integral at node NFEN-1
            VIBCPG = (1.0/6.0) * ( Sigma(NFEN) - Sigma(NFEN-1) ) * ( BCPG(NFEN) + 4.0d0 * BCPG_mid_point + BCPG(NFEN-1) )

            ! 4. evaluate the integral at the rest of the vertical nodes
            Do k = NFEN-2, 1, -1

                ! 4.1. 
               SigmaMP = ( Sigma(k+1) + Sigma(k) ) / 2.0d0

               f1 = ( (SigmaMP    - Sigma(k+1)) * (SigmaMP    - Sigma(k)) )
     &            / ( (Sigma(k+2) - Sigma(k+1)) * (Sigma(k+2) - Sigma(k)) )
               f2 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k)) ) 
     &            / ( (Sigma(k+1) - Sigma(k+2)) * (Sigma(k+1) - Sigma(k)) )
               f3 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k+1)) ) 
     &            / ( (Sigma(k)   - Sigma(k+2)) * (Sigma(k)   - Sigma(k+1)) )

               BCPG_mid_point = f1 * BCPG(k+2) + f2 * BCPG(k+1) + f3 * BCPG(k)

                ! 4.2. evaluate the integral
               VIBCPG = VIBCPG + (1.0/6.0) * ( Sigma(k+1) - Sigma(k) ) * ( BCPG(k+1) + 4.0d0 * BCPG_mid_point + BCPG(k) )

            End Do

            ! 5. let's hope we are better of now.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
            VIDBCPDXOH(NH)=REAL(VIBCPG)/AMB
            VIDBCPDYOH(NH)=AIMAG(VIBCPG)/AMB

         ENDIF

c     Finished computing baroclinic terms

      ENDDO

C     Finish loop over horizontal nodes to compute the horizontal velocity

#ifdef CMPI
C     Update new 3D baroclinic pressure gradient and the vertically
C     integrated baroclinic pressure gradient on ghost nodes
C
      CALL UPDATEC3D(BPG)
      CALL UPDATER(VIDBCPDXOH,VIDBCPDYOH,DUMV1,2)
#endif

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()


!----------------------------------------------------------------------------------------------------------------------------------------------------
! arash - print out a vertical slice of the BPG
!      Do NH = 1 , NP                             ! loop on all horizontal nodes
!         If ( abs (Y(NH)) < 1.0d-6 ) Then     ! we are on the slice
!             Hs = DP(NH) + IFNLFA * Eta2(NH)     !Total depth at current time step
!             HsOAMB = Hs / AMB
!             Do k = 1 , NFEN                     ! loop on vertical nodes
!                Zk = HsOAMB * (Sigma(k)-B) - DP(NH)        !determine z corresponding to sigma level k
!                Write(*,'(5E15.5E3)') X(NH), Y(NH), Zk, Real( BPG(NH,k) ), Aimag( BPG(NH,k) )
!             End Do
!         End if
!      End Do
!----------------------------------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------------------------------------------------
! arash - print out a vertical slice of Salinity 
! times for Kendra's problem: 1, 3.4, 4, 12, 20, 30
!       time_local = IT * DTDP 
!       day_local = time_local / ( 24.0d0 * 3600.0d0 )
! Kendra's example:
!      If ( (abs(time_local - 1.0d-3) < diff_tol) .or. 
!     &     (abs(time_local - 1.0d0) < diff_tol) .or.
!     &     (abs(time_local - 2.0d0) < diff_tol) .or.
!     &     (abs(time_local - 3.0d0) < diff_tol) .or.
!     &     (abs(time_local - 3.4d0) < diff_tol) .or.
!     &     (abs(time_local - 4.0d0) < diff_tol) .or.   
!     &     (abs(time_local - 8.0d0) < diff_tol) .or.
!     &     (abs(time_local -12.0d0) < diff_tol) .or.   
!     &     (abs(time_local -16.0d0) < diff_tol) .or.
!     &     (abs(time_local -20.0d0) < diff_tol) .or.   
!     &     (abs(time_local -25.0d0) < diff_tol) .or.
!     &     (abs(time_local -30.0d0) < diff_tol) .or.
!     &     (abs(time_local-100.0d0) < diff_tol) .or.
!     &     (abs(time_local-300.0d0) < diff_tol)  
!     &   ) Then

! Casey's example
!      If ( (abs(time_local  -   1.0d0) < diff_tol) .or. 
!     &     (abs(time_local  - 300.0d0) < diff_tol) .or.
!     &     (abs(time_local  -3600.0d0) < diff_tol) .or.   
!     &     (abs(time_local -21600.0d0) < diff_tol) .or. ! 6 hours
!     &     (abs(day_local -     1.0d0) < diff_tol) .or.   
!     &     (abs(day_local -     4.0d0) < diff_tol) .or.   
!     &     (abs(day_local -    10.0d0) < diff_tol) .or.   
!     &     (abs(day_local -    40.0d0) < diff_tol) .or.   
!     &     (abs(day_local -    100.0d0) < diff_tol) .or.   
!     &     (abs(day_local -    180.0d0) < diff_tol) 
!     &   ) Then

!         Do NH = 1 , NP                             ! loop on all horizontal nodes
!            If ( abs (Y(NH) - 0.05830d0 ) < 1.0d-3 ) Then     ! we are on the slice
!                Hs = DP(NH) + IFNLFA * Eta2(NH)     !Total depth at current time step
!                HsOAMB = Hs / AMB
!                Do k = 1 , NFEN                     ! loop on vertical nodes
!                   Zk = HsOAMB * (Sigma(k)-B) - DP(NH)        !determine z corresponding to sigma level k
!                   Write(*,'(5E15.5E3)') X(NH), Y(NH), Zk, Sal(NH,k)
!!                   Write(*,'(5E15.5E3)') X(NH), Y(NH), Zk, Real( BPG(NH,k) ), Aimag( BPG(NH,k) )
!                End Do
!            End if
!         End Do

!      End If
!----------------------------------------------------------------------------------------------------------------------------------------------------


      RETURN
C***********************************************************************
      END SUBROUTINE BPG3D
C***********************************************************************


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
      REAL(8) :: BCPressNN
      REAL(8) :: SigmaNN     !Sigma value of a neighbor node
C     jgf46.00 Explicitly declared the following variables
      INTEGER NN
      INTEGER J
      INTEGER LBelo
      INTEGER LAbov
      INTEGER LTry
      INTEGER IDiag
      REAL(8) SigBelo
      REAL(8) SigAbov
      REAL(8) SigTry

      call setMessageSource("zsurfbuoy")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE ZSURFBUOY
C**********************************************************************



C*************************************************************************
C     Subroutine to interpolate baroclinic pressure (BCP) to a specified
C     sigma value (SigmaNN) given an initial guess of which sigma
C     level is closest to the specified value.
C
C                                    R.L.  5/04/01
C                                    R.L.  5.19/03
C                                    a.f.  May 2015 - use cubic polynomials
C*************************************************************************
C
      SUBROUTINE ZSURFBUOY_p3 (SigmaNN,BCPressNN,NN,J)

      USE GLOBAL_3DVS
      IMPLICIT NONE
      REAL(8) :: BCPressNN
      REAL(8) :: SigmaNN     !Sigma value of a neighbor node
C     jgf46.00 Explicitly declared the following variables
      INTEGER NN
      INTEGER J
!      INTEGER :: LBelo = -999
!      INTEGER :: LAbov = -999
      INTEGER LBelo
      INTEGER LAbov
      INTEGER LTry
      INTEGER IDiag
      REAL(8) SigBelo
      REAL(8) SigAbov
      REAL(8) SigTry
! af: variables for cubic interpolation
      Integer  :: k1, k2, k3, k4                 ! nodes used for cubic polynomial
      REAL(8) :: f1, f2, f3, f4                 ! shape functions
      REAL(8) :: top, bot
      Integer  :: loop_situation = 0

      call setMessageSource("zsurfbuoy_p3")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IDiag=0

      IF(SigmaNN.LE.(1.00010d0*b)) THEN !if into ground then skip
         SigBelo=-999
         SigAbov=-999
         BCPressNN=-999.
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! try this Jan 20 2016: this is a bad idea
!         LBelo=1
!         BCPressNN=BCP(NN,LBelo)
!         SigBelo=b
!         SigAbov=b
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
         GOTO 100
      ENDIF
      IF((SigmaNN.GT.(1.0001d0*b)).AND.(SigmaNN.LE.b)) THEN !at bottom then use bottom
         LBelo=1
         BCPressNN=BCP(NN,LBelo)
         SigBelo=b
         SigAbov=b
         GOTO 100
      ENDIF

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
! arash 150513
      IF(SigmaNN.GE.(1.0001d0*a)) THEN            ! in the air - non-existent 
         BCPressNN=-999.0d0
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! try this Jan 20 2016: this is a bad idea >> instability
!         LAbov=NFEN
!         BCPressNN=BCP(NN,LAbov) ! this is zero
!         SigBelo=a
!         SigAbov=a
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
         GOTO 100
      ENDIF
      IF((SigmaNN.GE.a).and.(SigmaNN.LT.(1.0001d0*a))) THEN     ! surface
         LAbov=NFEN
         BCPressNN=BCP(NN,LAbov) ! this is zero
         SigBelo=a
         SigAbov=a
         GOTO 100
      ENDIF
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--


      LTry=J                    !start search for SIGABOV and SIGBELO
      SigTry=Sigma(LTry)
      IF(SigmaNN.GT.SigTry) THEN !too low
         loop_situation = 1
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
         loop_situation = 2
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
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
! linear interpolation
! 99   BCPressNN=(BCP(NN,LAbov)-BCP(NN,LBelo)) !interpolation
!     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + BCP(NN,LBelo)
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
 99   Continue

! arash added this because of the terminations at the SURA model. 
! this shouldn't be needed though.
! 0.
!      If ( loop_situation == 0 ) Then
!         BCPressNN = -999.0d0
!         GoTo 100
!      End If
! 1.
!      If ( (LAbov == -999) .or. (LBelo == -999) ) Then
!         BCPressNN = -999.0d0
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!         write(*,*) 'SigBelo, SigAbov = ', SigBelo, SigAbov
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         GoTo 100
!      End If
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
! arash - do some checking
! 0.
!      if ( J < 1 .or. J > NFEN ) then
!         write(*,*) 'J - mistake'
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         stop
!      end if

! 1.
!      if ( LAbov == LBelo ) then
!         write(*,*) 'LAbov == LBelo - mistake'
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         stop
!      end if
! 2.
!      if ( LAbov < LBelo ) then
!         write(*,*) 'LAbov < LBelo - mistake'
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         stop
!      end if
! 3.
!      if ( SigBelo > SigAbov ) then
!         write(*,*) 'SigBelo > SigAbov - mistake'
!         write(*,*) 'SigBelo, SigAbov = ', SigBelo, SigAbov
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         stop
!      end if

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
! use a cubic polynomial interpolation, instead of the above (99) linear interpolation
! 1- decide which nodes to pick: k1:lowest, k2, k3, k4:top-most
      If ( LAbov == NFEN ) Then
         k4 = NFEN
         k3 = NFEN - 1
         k2 = NFEN - 2
         k1 = NFEN - 3
      Else If ( LBelo == 1 ) Then
         k4 = 4
         k3 = 3
         k2 = 2
         k1 = 1
      Else
         k4 = LAbov + 1
         k3 = LAbov
         k2 = LBelo
         k1 = LBelo - 1
      End If

!      if ( (k1 == 0) .or. (k4 == nfen+1) ) then
!         write(*,*) 'wrong in ZSURFBUOY_p3'
!         write(*,*) 'Sig, BCPres, NN, J = ', SigmaNN, BCPressNN, NN, J
!         write(*,*) 'k1, k2, k3, k4 = ', k1, k2, k3, k4
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!      end if
!      if (k2 == k3) then
!!         write(*,*) 'wrong in ZSURFBUOY_p3'
!          write(*,*) 'Sig, BCPres, NN, J = ', SigmaNN, BCPressNN, NN, J
!         write(*,*) 'k1, k2, k3, k4 = ', k1, k2, k3, k4
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!      end if

! 2- compute the shape functions
      top = ((SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)) * 
     &       (SigmaNN - Sigma(k4))) 
      bot = ((Sigma(k1)-Sigma(k2)) * (Sigma(k1)-Sigma(k3)) * 
     &       (Sigma(k1)-Sigma(k4))) 
      f1  = top / bot

      top = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k3)) * 
     &       (SigmaNN - Sigma(k4)))
      bot = ((Sigma(k2)-Sigma(k1)) * (Sigma(k2)-Sigma(k3)) * 
     &       (Sigma(k2)-Sigma(k4)))
      f2  = top / bot

      top = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * 
     &       (SigmaNN - Sigma(k4)))
      bot = ((Sigma(k3)-Sigma(k1)) * (Sigma(k3)-Sigma(k2)) * 
     &       (Sigma(k3)-Sigma(k4)))
      f3  = top / bot

      top = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * 
     &       (SigmaNN - Sigma(k3)))
      bot = ((Sigma(k4)-Sigma(k1)) * (Sigma(k4)-Sigma(k2)) * 
     &       (Sigma(k4)-Sigma(k3)))
      f4  = top / bot

! 3- interpolate
      BCPressNN = f1 * BCP(NN,k1) + f2 * BCP(NN,k2) 
     &          + f3 * BCP(NN,k3) + f4 * BCP(NN,k4)

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
 100  CONTINUE

      IF(IDiag.EQ.2) THEN
         WRITE(2,*) '******** ZSURFBUOY_p3 **********'
         WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov',
     &        '      BCPressNN'
         WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,BCPressNN
 777     FORMAT(I7,I5,3(F10.3),E14.5)
      ENDIF

      RETURN
C**********************************************************************
      END SUBROUTINE ZSURFBUOY_p3
C**********************************************************************


!******************************************************************************
!    loop on all nodes, and apply a 1-layer Gaussian filter for BPG
! arash: January 20 2016
!******************************************************************************

       Subroutine Gaussian_filter_complex_BPG ( X, Y, NeiTab, NNeigh, quantity, NP, NFEN, MNEI )

       Use Global,      Only : H2
       Use Mesh,        Only : DP
       Use Global_3DVS, Only : SIGMA, AMB, B

       Implicit none
       Integer :: I, J, K
       Integer :: Neighbor_Node
       Integer :: MNEI, NP, NFEN
       Real(8) :: sigma_Gauss, length, length_max, factor, value, numerator, denominator, sum_Weights
       Real(8) :: velocity_magnitude
       Real(8) :: velocity_magnitude_max_in_xy
       Real(8) :: element_peclet_number
       Real(8) :: constant_diffusivity = 50.0d0
       Real(8) :: X(NP), Y(NP)
       INTEGER :: NeiTab(NP,MNEI), NNeigh(NP)
       REAL(8), ALLOCATABLE :: Weights(:)
       ! velocities are stored as complex quantities: u+i*v
       Complex(8), ALLOCATABLE :: quantity_smooth(:,:) 
       Complex(8), ALLOCATABLE :: inarray_vertical (:)
       Complex(8) :: quantity(NP,NFEN)
       Complex(8) :: interpolated_value

       REAL(8)  :: Hs, Hs_Neighbor_Node     !Total water depth at time level s
       REAL(8)  :: HsOAMB                   !Hs/(a-b)
       REAL(8)  :: Zk                       ! z depth of any node k in the vertical
       Real(8)  :: SigmaNN

       Allocate ( quantity_smooth(NP,NFEN) )
       Allocate ( inarray_vertical(NFEN) )

       quantity_smooth(:,:) = (0.0d0 , 0.0d0)

       ! loop on all nodes, and apply the Gaussian filter
       Do I = 1, NP
          ! 1. compute maximum distance between the node of interest, and the neighbor nodes
          !    this will be used to compute sigma2
          length_max = 0.0d0
          Do J = 2, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             length = SQRT ( ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 )
             If ( length > length_max ) length_max = length
          End Do

!----------------------------------------------------------------------------------------
          factor = 0.18d0
!----------------------------------------------------------------------------------------
          sigma_Gauss = length_max * factor

          ! 2. compute Gaussian filter weights
          IF ( ALLOCATED(Weights) ) DEALLOCATE( Weights )
          ALLOCATE( Weights( NNeigh ( I ) ) )
          Weights = 0.0d0
          Do J = 1, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             numerator = ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 
             denominator = 2.0d0 * sigma_Gauss**2
             Weights (J) = exp ( -1.0 * numerator / denominator )
          End Do
          ! 2.1 normalize weights
          sum_Weights = 0.0d0
          Do J = 1, NNeigh ( I )
             sum_Weights = sum_Weights + Weights (J)
          End Do
          If ( sum_Weights < 1.0d-6 ) Then
             Write(*,*) 'Possible error in Gaussian filter - error code = 1'
             CALL EXIT(1)
          End If
         Weights (:) = Weights (:) / sum_Weights
         ! 2.3 debug
         !Write(*,'(I5,11f6.3)') NNeigh ( I ), weights(1:NNeigh ( I ))
         !Pause
         ! 2.4. compute elevation corresponding to this node
         ! Total depth at current timestep - 'center' node
         ! arash: not sure if the next line is necessary, but probably not. we aren't using it in the bpg routine either.
         ! IF( Hs .LT. H0 ) Hs = H0
         HsOAMB = H2(I) / AMB

         ! 3. apply the Gaussian filter
          Do J = 1, NNeigh ( I )
             Neighbor_Node = NeiTab ( I , J )
             inarray_vertical (:) = quantity( Neighbor_Node , : )
             Do K = 1, NFEN
                ! 3.0. obtain the z and sigma value corresponding to the 'center' node
                Zk = HsOAMB * ( Sigma(K) - B ) - DP(I)
                SigmaNN = B + AMB * ( Zk + DP(Neighbor_Node) ) / H2(Neighbor_Node)
                ! 3.1. linearly interpolate surrounding BPG on the same horizontal plane
                Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, interpolated_value, K )
                quantity_smooth (I,K) = quantity_smooth (I,K) + Weights (J) * interpolated_value
             End Do
          End Do
       End Do

       ! transfer smooth bathymetry to the original array
       quantity (:,:) = quantity_smooth (:,:)

       End Subroutine Gaussian_filter_complex_BPG
