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
      SUBROUTINE TIMESTEP(IT)
C
#ifdef IEEE_DEBUG
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE GWCE
      USE MOMENTUM
      USE WRITE_OUTPUT, ONLY : writeOutput2D, writeHotStart, 
     &    writeWarnElev, collectInundationData, collectMinMaxData
      USE MESH, ONLY : NE, NP, DP, NM, X, Y, SLAM, SFEA, ICS, TotalArea,
     &                  MJU, Areas, SFAC
      USE BOUNDARIES, ONLY : NOPE, NETA, NBOU, NVEL, LBCODEI, NBV, SIII,
     &    NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, NFLUXIBP, NFLUXRBC, CSII,
     &    BARLANHT, BARLANCFSP, NVELL, IBCONN, BARINHT, BARINCFSP,
     &    BARINCFSB, PIPEHT, PIPECOEF, PIPEDIAM
      USE GLOBAL_IO
      USE HARM, ONLY : updateHarmonicAnalysis
      USE WIND, ONLY : getMeteorologicalForcing, PRBCKGRND_MH2O, rsget
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
C.... TCM V49.64.01 ADDITIONS FOR ICE
      USE OWI_ICE,ONLY : NCICE1_INIT,NCICE1_GET
C.....sb46.28sb03 added 09/xx/2006
      USE RS2,ONLY : RS2INIT,RS2GET
      USE NodalAttributes, ONLY : FRIC, HBREAK, FTHETA, FGAMMA, ! for wet/dry
     &     IFLINBF, IFNLBF, IFHYBF, ! for wet/dry
     &     Apply2DBottomFriction,
     &     Apply3DBottomFriction,
     &     LoadDirEffRLen,
     &     LoadEleSlopeLim,
     &     ManningsN, LoadManningsN, !sb46.28sb02
     &     Z0b_var, LoadZ0B_var,!Rjw IOOS 2010
     &     BFCdLLimit !sb46.28sb02/jgf47.04 lower limit of Cd for bot. friction
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, NOutGS, enforceWDOB,
     &    writeFort066, writeFort067, writeFort065, readFort020,
     &    readFort021, readFort019, enforceWDCB
c. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
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
      USE WEIR_FLUX

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
      REAL(SZ) PIPE_FLUX
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
      REAL(SZ) DARhoMRho0N1,DARhoMRho0N2,DARhoMRho0N3

      REAL(SZ) UN1
      REAL(SZ) UV0, UV1, UV2
      REAL(SZ) VIDBCPDXOHN1,VIDBCPDXOHN2,VIDBCPDXOHN3,VIDBCPDXOHAvgArea
      REAL(SZ) VIDBCPDYOHN1,VIDBCPDYOHN2,VIDBCPDYOHN3,VIDBCPDYOHAvgArea
      REAL(SZ) DPMIN      ! sb v46.28.sb05.06 11/01/2006
      REAL(SZ) DEta2DX,DEta2DY
      REAL(SZ) DRhoDX,DRhoDY

      REAL(8) AreaIE2,AreaEle
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
      REAL(8) TimeLoc, TimeH
      INTEGER WarnElevExceeded, ErrorElevExceeded
C   kmd48.33bc - added in the heat flux variables
      CHARACTER(80) :: CDUM80
      INTEGER :: NumofBCNode
      REAL(SZ), SAVE :: StaTimHS, RefTimHS
      INTEGER :: NOD
      REAL(SZ), ALLOCATABLE :: TMP(:,:)
      REAL(SZ) :: CD, CDQ, QWIND

C   TCM V49.64.01 -- ADDED FOR ICE CONCENTRATION FIELDS
      REAL(SZ) CICE_TRatio  !ICE VARIABLES

C   TCM v50.66.02 -- Added for Time Varying Bathymetry
      INTEGER :: IJK
      REAL(SZ) ETA2TMP,DPTMP,DPTMP2,BTRATIO  !tcm v50.66.01 bathymetry changes


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
      TimeLoc=IT*DTDP + StaTim*86400.D0
C    kmd48.33bc - added for the changes in the timestep in the hot start files
      IF (CHOTHS.eqv..true.) THEN
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
      IF (CHOTHS.eqv..true.) THEN
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
!            PR2(I)=0.D0   !tcm v49.16 20100617
            PR2(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
         ENDIF

C     TIP..Tidal potential forcing
         if(CTIP) then
            TIP1(I)=TIP2(I)
            TIP2(I)=0.D0
         endif

      END DO


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
#ifdef CMPI
            IF (MYPROC.EQ.0) then
               WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
               WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'       
            ENDIF
#else
            WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'        
#endif
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
#ifdef CMPI
            IF (MYPROC.EQ.0) then
               WRITE(ScreenUnit,'(A42,1X,E15.8,1X,A5)') 'BATHYMETRY VALUES ARE NOW FIXED AT TIME =',TIMELOC,'(SEC)'
               WRITE(16,'(A42,1X,E15.8,1X,A5)') 'BATHYMETRY VALUES ARE NOW FIXED AT TIME =',TIMELOC,'(SEC)'
            ENDIF
#else
            WRITE(ScreenUnit,'(A42,1X,E15.8,1X,A5)') 'BATHYMETRY VALUES ARE NOW FIXED AT TIME =',TIMELOC,'(SEC)'
            WRITE(16,'(A42,1X,E15.8,1X,A5)') 'BATHYMETRY VALUES ARE NOW FIXED AT TIME =',TIMELOC,'(SEC)'   
#endif
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

#ifdef CMPI
            IF (MYPROC.EQ.0) then
               WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
               WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'       
            ENDIF
#else
            WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'    
#endif               
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


C-------------------TEMPORAL RAMP-------------------------------------
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
Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
C    kmd48.33bc - ramp changes with baroclinic when timestep is changed
      CASE(1,2,3,4,5,6,7,8)
         Ramp=TANH((2.D0*TimeLoc/86400.D0)/DRamp)
         RampExtFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampExtFlux)
         RampIntFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampIntFlux)
         RampElev=TANH((2.D0*TimeLoc/86400.D0)/DRampElev)
         RampTip=TANH((2.D0*TimeLoc/86400.D0)/DRampTip)
Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
         RampMete=TANH((2.D0*(TimeLoc/86400.D0-DUnRampMete))/DRampMete)
         RampWRad=TANH((2.D0*TimeLoc/86400.D0)/DRampWRad)
      END SELECT
C
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
Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
            IF(NRamp.eq.8) then
               RampMete=TANH((2.D0*((((IT)*DTDP)/86400.D0)-DUnRampMete))/DRampMete)
            endif
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
C-------------------TEMPORAL RAMP-------------------------------------

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
              CICE1(I)=CICE2(I)
            END DO
            CALL NCICE1_GET(CICE2,NP)
          ENDIF
      ENDIF

C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
      if (ncice.ne.0) then
         CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
         do i=1,np
            CICEOUT(I) = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
         end do
      endif
C--------------------END ICE FIELDS----------------------------------------
C
C
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

         if((NRS.eq.1).or.(nrs.eq.2).or.(nrs.eq.3)) then
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
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN
Casey 090302: Added these lines for output to the rads.64 file.
               IF(ABS(NRS).EQ.3) then
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               ENDIF
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
               WAVE_A(I) = WAVE_A1(I) + RSTRATIO*(WAVE_A2(I)-WAVE_A1(I))
               WAVE_H(I) = WAVE_H1(I) + RSTRATIO*(WAVE_H2(I)-WAVE_H1(I))
               WAVE_T(I) = WAVE_T1(I) + RSTRATIO*(WAVE_T2(I)-WAVE_T1(I))
#endif
!#endif
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
               !  tcm v50.75 added for use whenever nrs=3 or nrs=4
               RSNXOUT(I) = RSX
               RSNYOUT(I) = RSY
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
         IF ((ABS(RES_BC_FLAG).GE.1).AND.(NOPE.GT.0)) THEN
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
         IF (NFFR.GT.0) THEN
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
C...  COBELL - MOVED TO WEIR_BOUNDARY.F
        IF(NFLUXB.EQ.1) THEN
          DO I=1,NVEL
            SELECT CASE(LBCODEI(I))
            CASE(3,13,23)
                CALL COMPUTE_EXTERNAL_BOUNDARY_FLUX(I,TIMELOC,QN2(I))
            END SELECT 
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


      if(subdomainOn.and.enforceBN.eq.1) call readFort019(it)  ! NCSU Subdomain
      if(subdomainOn.and.enforceBN.eq.2) call readFort020(it)  ! NCSU Subdomain
      if(subdomainOn.and.enforceBN.eq.2) call readFort021(it)  ! NCSU Subdomain


C...
C...  Compute the water surface elevation from the GWCE form of the 2D
C...  continuity eq.
C...

ckmd    Changed to include the predictor-corrector algorithm
      IF(CPRECOR) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
c...tcm added call for slope limiting         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
         CALL Mom_Eqs_New_NC()
C...  If running in parallel, update velocities & fluxes on all processors

#ifdef CMPI
         CALL UPDATER(UU2,VV2,DUMY1,2)
         CALL UPDATER(QX2,QY2,DUMY1,2)
#endif
         CALL GWCE_New_pc(IT,TimeLoc,TimeH)
c...tcm added call to slope limiter         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
      ENDIF
      IF(CGWCE_New) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
c.... tcm added call to slope limiter         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
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
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
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
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
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
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
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

         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
         CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif


Cwet...
CWET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
C*******************************************************************************************
C tcm v50.66.01 -- This is an additional test for wetting only when time varying 
C                  bathymetry is used and is only performed during the period of 
C                  bathymetry evolution.
C
         IF ((NDDT.NE.0).AND.(IT.LE.BTIME_END+1) ) THEN
            DO I=1,NE
               NM1=NM(I,1)
               NM2=NM(I,2)
               NM3=NM(I,3)
               NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
               IF(NCTOT.lt.3) THEN   !If not wet from previous time step
                  NCTOT=NNODECODE(NM1)+NNODECODE(NM2)+NNODECODE(NM3)
                  if(NCTOT.lt.3) then !if not alreay made wet for this time step
                     ETAN1=ETA2(NM1)
                     ETAN2=ETA2(NM2)
                     ETAN3=ETA2(NM3)
                     HTOTN1=DP(NM1)+ETA2(NM1)
                     HTOTN2=DP(NM2)+ETA2(NM2)
                     HTOTN3=DP(NM3)+ETA2(NM3)

!                    if all nodes have a depth greater than or equal to 
!                    hoff = 1.2*H0, then make the element wet
                     IF( (HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF).AND.
     &                   (HTOTN3.GE.HOFF) ) THEN
                     !THE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                        !  Make Node 1 Wet and set parameters
                        IF(NNODECODE(NM1).NE.1) THEN  !node 1
                           NNODECODE(NM1)=1
                           NM123 = NM1
                           IF(C2DDI)THEN
C<<                           Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                              /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL=VELMIN
                           H1 = HTOTN1
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                          (IFNLBF+IFHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
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
     &                              /(SIGMA(2)-SIGMA(1))
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
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
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
                        ENDIF  !end node 1

                        !  Make Node 2 Wet and set parameters
                        IF (NNODECODE(NM2).NE.1) THEN
                           NNODECODE(NM2) = 1
                           NM123=NM2
                           IF(C2DDI)THEN
C<<                        Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL = VELMIN
                           H1 = HTOTN2
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                          (IFNLBF+IFHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
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
     &                              /(SIGMA(2)-SIGMA(1))
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
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                 (DP(NM123)+IFNLFA*ETA2(NM123)))
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
                        ENDIF  !node 2

                        !  Make Node 3 Wet and set parameters
                        IF(NNODECODE(NM3).NE.1) THEN
                           NNODECODE(NM3)=1
                           NM123 = NM3
                           IF(C2DDI)THEN
C<<                         Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                           /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL=VELMIN
                           H1 = HTOTN3
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                           (IFNLBF+IFHYBF*
     &                           (1.D0+(HBREAK/H1)**FTHETA)**
     &                           (FGAMMA/FTHETA)))
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
     &                             /(SIGMA(2)-SIGMA(1))
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
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
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
                        ENDIF !node3

                     ENDIF  !ALL DEPTHS GREATER THAN HOFF
                  ENDIF  !IF NNODECODE SUM LESS THAN 3
               ENDIF   ! IF NODECODE SUM LESS THAN 3
            ENDDO  !LOOP OVER ELEMENTS

            if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
            CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif
C
         ENDIF !IT TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME
C
C.... END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY
CWET..
CWET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
C*******************************************************************************************

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
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
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
               endif
#else
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
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
               endif !noffActive.eqv..true.
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
               if (noffActive.eqv..true.) then 
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
              endif
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

         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain


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

C...
C...  2DDI Momentum Equation Solution
C...
      IF (C2DDI) THEN
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
      !CALL writeOutput2D(IT,TimeLoc) ! =>zc - moved to adcirc.F so that SWAN writes 
                                      !        correctly
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
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,
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
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,
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
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,
     &           2X,'ON MYPROC = ',I4,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev,',
     &           ' ADCIRC stopping. **')
         ErrorElevExceeded = 1                ! Clint's Zombie Slayer
         ENDIF
         call WarnElevSum(ErrorElevExceeded)  ! Clint's Zombie Slayer 2010.08.07
         IF( ErrorElevExceeded /= 0 ) THEN    !  Finalize MPI Environment,
           Flag_ElevError = .true. 
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
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8)
         ENDIF
         IF(ELMAX.GT.WarnElev) THEN
            WRITE(ScreenUnit,1994)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1994) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,
     &           2X,'** WARNING: Elevation.gt.WarnElev **')
#ifdef DEBUG_WARN_ELEV
            IF (WarnElevDump) CALL WriteWarnElev(TimeLoc, IT)
#endif
         ENDIF
         IF(ELMAX.GT.ErrorElev) THEN
            Flag_ElevError = .true.
            WRITE(ScreenUnit,1996)
     &           IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1996) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1996       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           5X,'TIME = ',E15.8,
     &           /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8,
     &           2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,/,
     &           2X,'** ERROR: Elevation.gt.ErrorElev, '
     &           'ADCIRC stopping. **')
            CALL ADCIRC_Terminate()
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

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

C******************************************************************************
      END SUBROUTINE TIMESTEP
C******************************************************************************


C****************************************************************************************
C   Subroutine to compute Scalar Transport                                              *
C                                                                                       *
C   Note, this is not set up for parallel operation.                                    *
C****************************************************************************************
C
      SUBROUTINE SCALAR_TRANS_2D (IT,TimeLoc)

      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, SFAC
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
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY: IFNLFA, IFNLCT,
     &      NODECODE, NOFF, VIDBCPDXOH, VIDBCPDYOH,
     &      setMessageSource, unsetMessageSource, allMessage,
     &      logMessage, DEBUG, ECHO, INFO, WARNING, ERROR,
     &      ETA2, SIGT0
      USE MESH, ONLY : NM, X, Y, NP, DP, AREAS, NODELE, 
     &      NEITABELE, NEITAB, NNEIGH, SFAC
      USE GLOBAL_3DVS, ONLY : BCP, BPG, NFEN, SIGT, SIGMA, IDEN,
     &      GORhoOAMB, AMB, B, IY
#ifdef CMPI
      USE MESSENGER
      IMPLICIT NONE
      REAL(SZ) :: DUMV1(1),DUMV2(1)
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
C     Subroutine to determine local advection state for each element
C
C                                   Corbitt 120322
C*************************************************************************

      SUBROUTINE ADVECTLOCAL(IE)

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE MESH, ONLY : NM, DP
      USE NodalAttributes, ONLY : AdvectionState

      Implicit None

      Integer  :: NM1, NM2, NM3,IE

      NM1=NM(IE,1)
      NM2=NM(IE,2)
      NM3=NM(IE,3)

      IF ((DP(NM1).GE.AdvectionState(NM1)).AND.
     &    (DP(NM2).GE.AdvectionState(NM2)).AND.
     &    (DP(NM3).GE.AdvectionState(NM3))) THEN
             IFNLCT = IFNLCTE
             IFNLCAT = IFNLCATE
      ELSE
         IFNLCT = 0
         IFNLCAT = 0
      ENDIF

C***********************************************************************
      END SUBROUTINE ADVECTLOCAL
C***********************************************************************

C*************************************************************************
C   Subroutine to check if elemental slope limiting is needed
C       Chris Massey, USACE-ERDC-CHL Dec. 8, 2014 -- Removed this
C         section from timestep subroutine and put it in its own
C         subroutine.
C*************************************************************************
      subroutine check_slopes(it,TimeLoc)
      USE SIZES, ONLY : SZ, mnproc
      use global, only : eta2,nodecode,NOFF,ESLONOFF,screenUnit,
     &      MyProc,setMessageSource, unsetMessageSource, scratchMessage,
     &      allMessage,logMessage, DEBUG, ECHO, INFO, WARNING, ERROR,
     &      nodes_lg
      use mesh, only : x, y, nm, ne, areas, sfac
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active,
     &     elemental_slope_limiter_grad_max,
     &     elemental_slope_limiter_max_exceeded
      implicit none
      INTEGER, intent(in) :: IT
      Real(8), intent(in) :: TimeLoc
      INTEGER IE,I
      INTEGER NM1, NM2, NM3, NM123
      INTEGER NC1, NC2, NC3, NCEle, NCI
      REAL(SZ) DEta2DX,DEta2DY,DEta2Mag,SFacAvg
      REAL(SZ) FDX1,FDX2,FDX3,FDY1,FDY2,FDY3
      REAL(8) :: AreaIE2
      integer nodeNumber ! fulldomain node number where wse slope is exceeded
      

      call setMessageSource("check_slopes")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! jgf51.47: Per Rick's specification, start off by deactivating
      ! any slope limitation triggered in the previous time step;
      ! we have a fresh start each time step in determining whether
      ! to activate slope limiting at each node. 
      ! jgf51.50: Added this back to the subroutine version of the
      ! slope limiter. 
      elemental_slope_limiter_active(:) = .false.

Cobell  CHECK IF THE LOCAL GRADIENT FOR ELEMENTS IS EXCEEDED AND TRIGGER THE
C...  ELEMENTAL SLOPE LIMITER ACCORDINGLY (original routine by Crystal Fulcher)

      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         IF(NCEle.EQ.0)THEN
            CYCLE  ! this element is dry, go to the next one
         ENDIF
         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
         AreaIE2=Areas(IE)
         FDX1=(Y(NM2)-Y(NM3))*SFacAvg
         FDX2=(Y(NM3)-Y(NM1))*SFacAvg
         FDX3=(Y(NM1)-Y(NM2))*SFacAvg
         FDY1=X(NM3)-X(NM2)
         FDY2=X(NM1)-X(NM3)
         FDY3=X(NM2)-X(NM1)
         dEta2Dx  = (Eta2(NM1)*FDX1+Eta2(NM2)*FDX2+Eta2(NM3)*FDX3)
     &                 /AreaIE2
         dEta2Dy  = (Eta2(NM1)*FDY1+Eta2(NM2)*FDY2+Eta2(NM3)*FDY3)
     &                 /AreaIE2
         dEta2Mag = sqrt(dEta2Dx*dEta2Dx + dEta2Dy*dEta2Dy)
         !
         ! jgf51.51: Now that the slope limiter gets reset at 
         ! every time step, I had to rewrite the logging so that
         ! a log message is only written the first time the slope
         ! limiter is activated at a node during a particular run.
         DO I=1,3
            ! If the limiter is on already, go to the next node. 
            IF (elemental_slope_limiter_active(NM(IE,I))) CYCLE
            ! Compare the elemental slope to the maximum elemental gradient.
            IF (dEta2Mag.GE.
     &           ABS(elemental_slope_limiter_grad_max(NM(IE,I)))) THEN
               ! jgf51.51: Log the fulldomain node number.
               nodeNumber = nm(ie,i)
               if (mnproc.gt.1) then
                  nodeNumber = nodes_lg(nm(ie,i))
               endif
               ! zc - If gradmax is positive or zero, activate slope
               ! limiting. 
               if (elemental_slope_limiter_grad_max(nm(ie,i))
     &               .ge.0.0d0) then
                  ! If it is the first time that the slope limiter
                  ! has been activated at this node, write a log 
                  ! message. 
                  if (eslonoff(nm(ie,i)).eq.0) then
                     write(scratchMessage,1983) nodeNumber,dEta2Mag, 
     &                     elemental_slope_limiter_grad_max(NM(IE,I)),
     &                     it, timeLoc
                     call allMessage(INFO,scratchMessage)
                     eslonoff(nm(ie,i)) = 1 ! for output file
                  endif
                  elemental_slope_limiter_active(NM(IE,I)) = .TRUE.
               else
                  ! Just print log message the first time the 
                  ! gradient is exceeded. 
                  IF (elemental_slope_limiter_max_exceeded(NM(IE,I))
     &                  .eqv..false.) THEN
                     write(scratchMessage,1984) nodeNumber,dEta2Mag, 
     &                     elemental_slope_limiter_grad_max(NM(IE,I)),
     &                     it, timeLoc
                     call allMessage(INFO,scratchMessage)
                     elemental_slope_limiter_max_exceeded(NM(IE,I)) =
     &                      .true.
                  endif
               endif
            endif
         enddo ! loop around nodes of an element
      enddo ! loop over the elements


1983  format('Elemental slope limiter turned on at fulldomain node ',i0,
     &  ' where the elemental slope is ',1pE12.4E3,
     &  ' and the maximum elemental slope is ',1pE12.4E3,
     &  ' on time step ',i0,' and time = ',e15.8,'.')

1984  format('Maximum elemental slope exceeded at fulldomain node ',i0,
     &  ' where the elemental slope is ',1pE12.4E3,
     &  ' and the maximum elemental slope is ',1pE12.4E3,
     &  ' on time step ',i0,' and time = ',e15.8,'.')
     
C...  END CHECKING ELEMENT GRADIENTS

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !-----------------------------------------------------------------
      END SUBROUTINE CHECK_SLOPES
      !-----------------------------------------------------------------   
C
C***********************************************************************

C***********************************************************************
C  Apply Elemental Slope Limiter 
C    Chris Massey, USACE-ERDC-CHL, Dec. 8, 2014 
C       Made into a subroutine
C
C***********************************************************************
C
      SUBROUTINE APPLY_SLOPE_LIMITS(ETA2Lim,LocNP)
      USE SIZES, ONLY : SZ
      use global, only : nodecode,NOFF,IFNLFA,
     &      setMessageSource, unsetMessageSource, allMessage,
     &      logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      use mesh, only : ne, nm, areas, totalArea
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active
      implicit none
      integer :: IE,NM1,NM2,NM3,NC1,NC2,NC3,NCEle
      integer, intent(in) :: LocNP
      REAL(SZ) :: EtaN1,EtaN2,EtaN3,EtaN123
      real(sz), intent(inout) :: Eta2lim(LocNP)
      REAL(8) :: AreaEle
      REAL(sz), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
      LOGICAL, SAVE :: firstCall = .true.


      call setMessageSource("check_slopes")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      IF (LoadEleSlopeLim.eqv..true.) THEN
         IF (firstCall.eqv..true.) THEN
            allocate(elevSum(LocNP))
            firstCall = .false.
            elevSum(:) = 0.d0
         ENDIF
      ENDIF


C       ELEMENTAL SLOPE LIMITER
C
C        CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
C        TIME STEP AND APPLY SMOOTHING IF DESIRED. -ZC

      IF (LoadEleSlopeLim.eqv..true.) THEN
C
         elevSum(:) = 0.d0
         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            EtaN1=IFNLFA*Eta2Lim(NM1)
            EtaN2=IFNLFA*Eta2Lim(NM2)
            EtaN3=IFNLFA*Eta2Lim(NM3)
            AreaEle=NCEle*Areas(IE)/2.d0
            EtaN123=(EtaN1+EtaN2+EtaN3)/3.d0
            elevSum(NM1)=elevSum(NM1)+AreaEle*EtaN123
            elevSum(NM2)=elevSum(NM2)+AreaEle*EtaN123
            elevSum(NM3)=elevSum(NM3)+AreaEle*EtaN123
         ENDDO

Cobell   CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
C        TIME STEP AND APPLY SMOOTHING IF DESIRED. THIS ROUTINE PARALLELS THE
C        ABOVE ROUTINE.
         WHERE ((elemental_slope_limiter_active.eqv..true.).and.
     &         (TotalArea.ne.0.d0))
               Eta2Lim = elevSum / TotalArea
         END WHERE
      ENDIF

C... Will apply the updating outside the subroutine
!#ifdef CMPI
!      CALL UPDATER(ETA2Lim,DUMY1,DUMY2,1)
!#endif


#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()


      RETURN

      END SUBROUTINE APPLY_SLOPE_LIMITS

C***********************************************************************

