! define SB_WETDRY
!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 45.12                                    *
!                                                                             *
! The timestepping module is configured to allow selection of a number of     *
! alternative algorithms within the overall FE framework.  These algorithms   *
! are selected by the TRUE/FALSE state of the logical variables listed below. *
! These variables are set in READ_INPUT.F, based on the value of the fort.15  *
! input parameter IM.  The only exception is CGWCE_Lump which is set by a     *
! preprocessor flag at compile time. The variables are passed in GLOBAL.F     *
!                                                                             *
! Logical Variable List (default value .FALSE., set in global.f)              *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     C2DDI            - 2D Depth Integrated model run                        *
!     C3D              - 3D model run                                         *
!     C3DDSS           - Stress form of the 3D momentum equations             *
!     C3DVS            - Velocity form of the 3D momentum equations           *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     C2D_BTrans       - Include a 2D baroclinic transport calculation        *
!     C2D_PTrans       - Include a 2D passive transport calculation           *
!     C3D_BTrans       - Include a 3D baroclinic transport calculation        *
!                        (used in 3D subroutines only)                        *
!     C3D_PTrans       - Include a 3D passive transport calculation           *
!                        (used in 3D subroutines only)                        *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CBaroclinic      - Include baroclinic terms                             *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CGWCE_New        - New ADCIRC GWCE formulation (old algorithm, new code)*
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CGWCE_Lump       - Lump the GWCE matrix                                 *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     GWCE_New() and GWCE_New_pc() only:                                      *
!     CGWCE_LS_KGQ     - Kolar-Gray, flux-based, lateral stress formulation   *
!                        in the GWCE (same as original formulation)           *
!     CGWCE_LS_2PartQ  - 2 Part, flux-based, lateral stress formulation       *
!                        in the GWCE                                          *
!     CGWCE_LS_2PartV  - 2 Part, velocity-based, lateral stress formulation   *
!                        in the GWCE                                          *
!     CGWCE_LS_2PartSQ - 2 Part, flux-based, symmetric lateral stress         *
!                        formulation in the GWCE                              *
!     CGWCE_LS_2PartSV - 2 Part, velocity-based, symmetric lateral stress     *
!                        formulation in the GWCE                              *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     GWCE_New() and GWCE_New_pc() only:                                      *
!     CGWCE_Advec_NC   - Non-conservative advection formulation in the GWCE   *
!                        (same as original formulation)                       *
!     CGWCE_Advec_C1   - Use conservative advection formulation 1 in the GWCE *
!     CGWCE_Advec_C2   - Use conservative advection formulation 2 in the GWCE *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CME_Orig         - Original Momentum Eq. formulation                    *
!     CME_New_NC       - Non-conservative advection formulation in the        *
!                        Momentum Eqs. (same as original formulation)         *
!     CME_New_C1       - Conservative advection formulation 1 in the          *
!                        Momentum Eqs.                                        *
!     CME_New_C2       - Conservative advection formulation 2 in the          *
!                        Momentum Eqs.                                        *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CPRECOR          - Use the predictor-corrector algorithm for GWCE       *
!                        and momentum equations (package deal)                *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     The following are not available in Mom_Eqs_Original()                   *
!     CME_LS_IBPQ      - Integration by parts, flux-based, lateral stress     *
!                        formulation in the Momentum Eqs.                     *
!     CME_LS_IBPV      - Integration by parts, velocity-based, lateral stress *
!                        formulation in the Momentum Eqs.                     *
!                        (same as original formulation)                       *
!     CME_LS_IBPSQ     - Integration by parts, flux-based, symmetric lateral  *
!                        stress formulation in the Momentum Eqs.              *
!     CME_LS_IBPSV     - Integration by parts, velocity-based, symmetric      *
!                        lateral stress formulation in the Momentum Eqs.      *
!     CME_LS_2PartQ    - 2 Part, flux-based, lateral stress formulation in    *
!                        the Momentum Eqs.  (NOT IMPLEMENTED)                 *
!     CME_LS_2PartV    - 2 Part, velocity-based, lateral stress formulation in*
!                        the Momentum Eqs.  (NOT IMPLEMENTED)                 *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CME_AreaInt_Orig - Original area integration in the Momentum Eqs.       *
!                        (incorrect, but same as original formulation)        *
!     CME_AreaInt_Corr - Corrected area integration in the Momentum Eqs.      *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CSmag_Eh         - Use Smagorinski, spatially varying, vertically       *
!                        integrated lateral viscosity coefficient             *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!   See header.f for a summary history of code modifications.                 *
!******************************************************************************

    SUBROUTINE TIMESTEP(IT)

#ifdef IEEE_DEBUG
    USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
    USE SIZES, ONLY : SZ
    USE GLOBAL
    USE WRITE_OUTPUT, ONLY : writeOutput2D, writeHotStart, &
    writeWarnElev, collectInundationData, collectMinMaxData
    USE MESH, ONLY : NE, NP, DP, NM, X, Y, SLAM, SFEA, ICS, TotalArea, &
    MJU, Areas, SFAC
    USE BOUNDARIES, ONLY : NOPE, NETA, NBOU, NVEL, LBCODEI, NBV, SIII, &
    NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, NFLUXIBP, NFLUXRBC, CSII, &
    BARLANHT, BARLANCFSP, NVELL, IBCONN, BARINHT, BARINCFSP, &
    BARINCFSB, PIPEHT, PIPECOEF, PIPEDIAM
    USE GLOBAL_IO
    USE HARM, ONLY : updateHarmonicAnalysis
    USE WIND
    USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
!.....sb46.28sb01 added 09/xx/2006
    USE OWIWIND,ONLY : NWS12INIT,NWS12GET
!.... TCM V49.64.01 ADDITIONS FOR ICE
    USE OWI_ICE,ONLY : NCICE1_INIT,NCICE1_GET
!.....sb46.28sb03 added 09/xx/2006
    USE RS2,ONLY : RS2INIT,RS2GET
    USE NodalAttributes, ONLY : FRIC, HBREAK, FTHETA, FGAMMA,  & ! for wet/dry
    IFLINBF, IFNLBF, IFHYBF,  & ! for wet/dry
    Apply2DBottomFriction, &
    Apply3DBottomFriction, &
    ApplyDirectionalWindReduction, &
    LoadDirEffRLen, &
    ApplyCanopyCoefficient, &
    LoadCanopyCoef, &
    LoadEleSlopeLim, &
    ManningsN, LoadManningsN,  & !sb46.28sb02
    Z0b_var, LoadZ0B_var, & !Rjw IOOS 2010
    BFCdLLimit !sb46.28sb02/jgf47.04 lower limit of Cd for bot. friction
    USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, NOutGS, enforceWDOB, &
    writeFort066, writeFort067, writeFort065, readFort020, &
    readFort021, readFort019, enforceWDCB
!. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
    USE GLOBAL_3DVS, ONLY: &
    A, B, BSX, BSY, EVTOT, ISLIP, KP, Q, SIGMA, Z0B,NFEN

#ifdef CMPI
    USE MESSENGER
    USE HSWRITER, ONLY: writeHotstart_through_HSwriter  !st3 100711 for hsfile
#endif
          
#ifdef CSWAN
! Casey 090302: We need these values from other places.
    USE OWIWIND,     ONLY: WindMultiplier
#ifndef CSWANFRIC
    USE Couple2Swan, ONLY: ComputeWaveDrivenForces, &
    CouplingInterval, COUPWIND, InterpoWeight, SWAN_WX2, &
    SWAN_WY2,WriteSwanHotStart
#else
    USE Couple2Swan, ONLY: ComputeWaveDrivenForces, &
    CouplingInterval, COUPWIND, InterpoWeight, SWAN_WX2, &
    SWAN_WY2,WriteSwanHotStart,ComputeModifiedWaveFriction, &
    ComputeWaveFrictionProperties,TKXX, TKXY, TKYY, WAVE_A, &
    WAVE_A1, WAVE_A2, WAVE_H, WAVE_H1, WAVE_H2, WAVE_T, &
    WAVE_T1, WAVE_T2
#endif
#endif
    USE WEIR_FLUX

    IMPLICIT NONE
    INTEGER, intent(in) :: IT

#ifdef CSWAN
#ifdef CSWANFRIC
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
    REAL(SZ) :: WAVA
    REAL(SZ) :: WAVH
    REAL(SZ) :: WAVT
#endif
#endif

!. RJW merged 08/26/2008 from Casey 071219:Added the following variables for 3D wet/dry.
    CHARACTER(6) :: TEMPDIRNAME
    COMPLEX(SZ) :: DUDS !jgf48.50 declare size SZ instead of plain COMPLEX
    REAL(SZ) :: KSLIP

    INTEGER :: IE, I, J, K                  !local loop counters
    INTEGER :: NM1, NM2, NM3, NM123
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    INTEGER :: NCyc, NA,NITEMS
    INTEGER :: NMN1,NMN2,NMN3,NWETNEI,NWETADJ  ! sb v46.28.sb05.06 11/01/2006

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
    REAL(SZ) WTRatio, WindX, WindY, WindMag, WDragCo
    REAL(SZ) DARhoMRho0N1,DARhoMRho0N2,DARhoMRho0N3

    REAL(SZ) UN1
    REAL(SZ) UV0, UV1, UV2
    REAL(SZ) VIDBCPDXOHN1,VIDBCPDXOHN2,VIDBCPDXOHN3,VIDBCPDXOHAvgArea
    REAL(SZ) VIDBCPDYOHN1,VIDBCPDYOHN2,VIDBCPDYOHN3,VIDBCPDYOHAvgArea
    REAL(SZ) DPMIN      ! sb v46.28.sb05.06 11/01/2006
    REAL(SZ) PRDIFF,PRBCKGRND_MH2O     ! tcm v49.16 20100617 added
    REAL(SZ) DEta2DX,DEta2DY
    REAL(SZ) DRhoDX,DRhoDY

    REAL(8) AreaIE2,AreaEle
    REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
    REAL(8) FDX1O2A, FDX2O2A, FDX3O2A, FDY1O2A, FDY2O2A, FDY3O2A
    REAL(8) TimeLoc, TimeH
    INTEGER :: WarnElevExceeded, ErrorElevExceeded
    REAL(SZ) HollandTime, AsymmetricTime, GeneralAsymTime !for parametric hurricane models
! jgf49.1001 factor used to blend vortex and background winds for NWS29:
    REAL(SZ) bf ! blending factor, 0.0 to 1.0
!   kmd48.33bc - added in the heat flux variables
    CHARACTER(80) :: CDUM80
    INTEGER :: NumofBCNode
    REAL(SZ), SAVE :: StaTimHS, RefTimHS
    INTEGER :: NOD
    REAL(SZ), ALLOCATABLE :: TMP(:,:)
    REAL(SZ) :: CD, CDQ, QWIND

!   TCM V49.64.01 -- ADDED FOR ICE CONCENTRATION FIELDS
    REAL(SZ) PIC,CICE_TRatio  !ICE VARIABLES

!   TCM v50.66.02 -- Added for Time Varying Bathymetry
    INTEGER :: IJK
    REAL(SZ) ETA2TMP,DPTMP,DPTMP2,BTRATIO  !tcm v50.66.01 bathymetry changes


    call setMessageSource("timestep")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!.... tcm v49.16 20100617 added
!.... convert background pressure from millibars to meters of water
    PRBCKGRND_MH2O = 100.0D0*PRBCKGRND/(RHOWAT0*G)

!    kmd48.33bc - changed for timestep changes in the hot start files
!     jgf46.21 Combined flux/radiation b.c. for rivers
#ifdef IBM
    IF (CHOTHS.eqv. .TRUE. ) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
    ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP,KIND(0.0D0))
    END IF
#else
    IF (CHOTHS.eqv. .TRUE. ) THEN
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
    ELSE
        FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDP)
    END IF
#endif

!...  COMPUTE MASTER TIME WHICH IS REFERENCED TO THE BEGINNING TIME OF
!...  THE MODEL RUN
!...
    TimeLoc=IT*DTDP + StaTim*86400.D0
!    kmd48.33bc - added for the changes in the timestep in the hot start files
    IF (CHOTHS.eqv. .TRUE. ) THEN
        IF ((ITHS+1) == IT) THEN
            StaTimHS=((IT-1)*DTDPHS)/86400.D0
            RefTimHS=((IT-1)*DTDP)/86400.D0
        END IF
        TimeLoc=IT*DTDP + (StaTimHS - RefTimHS)*86400.D0
    END IF

!...  HARMONIC CALCULATIONS ARE MADE FOR TIME WHICH INCLUDES THE REFTIM
!...  TO ALLOW FOR THE POSSIBILITY THAT THE EQUILIBRIUM ARGUMENTS MAY
!...  BE FOR A TIME OTHER THAN THE MODEL STARTING TIME.
!...
    TimeH=IT*DTDP + (StaTim - RefTim)*86400.D0
! md - added this for the cases where the timestep changes from hot start timestep
    IF (CHOTHS.eqv. .TRUE. ) THEN
        IF ((ITHS+1) == IT) THEN
            StaTimHS=((IT-1)*DTDPHS)/86400.D0
            StaTim=((IT-1)*DTDP)/86400.D0
        END IF
        TimeH=IT*DTDP + ((StaTimHS - StaTim) - RefTim)*86400.D0
    END IF


!...  SHIFT THE FLUX PER UNIT WIDTH, DEPTH AVERAGED VELOCITIES, BOTTOM STRESS,
!...  WIND STRESS, SURFACE PRESSURE AND TIDAL POTENTIALS TO PREVIOUS TIME STEP.
!...  ZERO OUT THE NEW FORCING TERMS AND LOAD VECTORS
!...  COMPUTE A NEW BOTTOM FRICTION COEFFICIENT
!...
! md  Shift values in time for predictor-corrector algorithm

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

    !...  Transport
        IF(IM == 10) THEN
            TRANS_LV_B(I)=0.D0
            TRANS_LV_A(I)=0.D0
        ENDIF

    !...  Wind (& wave radiation stress if used)
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN
            WSX1(I)=WSX2(I)
            WSX2(I)=0.D0
            WSY1(I)=WSY2(I)
            WSY2(I)=0.D0
            PR1(I)=PR2(I)
        !            PR2(I)=0.D0   !tcm v49.16 20100617
            PR2(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
        ENDIF

    !     TIP..Tidal potential forcing
        if(CTIP) then
            TIP1(I)=TIP2(I)
            TIP2(I)=0.D0
        endif

    END DO


!...TCM V50.66.01 -- ADDING TIME DEPENDENT BATHYMETRY
!... so that total water column height is unchanged

!  DP is linearly interpolated between DP1 and DP2
!  during the time interval btime1 and btime_end.
!  After btime_end, DP is equal to DP2

!      DP1                         DP2
!    btime1                      btime2      btime2 = btime1 + btiminc
!       |---------x-----------------|        btime_end = btime1 + btime_end < btime2
!              btime_end


    IF(abs(NDDT) == 1) THEN
    !        Get a new bathymetry from file if time to do so
        IF(TimeLoc > BTIME2) THEN !determine if bathy file time incr. is exceeded
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
        !...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
        !...     DEPTHS ARE > OR = TO H0.
            IF ((NOLIFA == 0) .OR. (NOLIFA == 1)) THEN
                DO I=1,NP
                    IF (DP2(I) < H0) DP2(I) = H0
                ENDDO
            ENDIF
#ifdef CMPI
            IF (MYPROC == 0) then
                WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
                WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            ENDIF
#else
            WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
#endif
        ENDIF  !test for updating bathymetry records
    !.......If time is during the bathymetry change interval, then update bathymetry
        IF(timeloc < BTIME_END) THEN  !tcm 20150728 changed <= to < 
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
        IF(timeloc == BTIME_END) THEN
#ifdef CMPI
            IF (MYPROC == 0) then
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

    IF(abs(NDDT) == 2) THEN
    !        Get a new bathymetry from file if time to do so
        IF(TimeLoc > BTIME2) THEN !determine if bathy file time incr. is exceeded
            BTIME1=BTIME2  !new starting time for this record
            BTIME2=BTIME2+BTIMINC  !new ending time for this record
            BTIME_END = BTIME1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
            DO I=1,NP
                dp1(I) = dp2(I)  ! move current data to old
            END DO

        !!!        go get new record for only some nodes, all
        !!!        other nodes keep their current value
            CALL NDDT2GET( 141,DP2(:),-99999.d0 )
        !...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
        !...     DEPTHS ARE > OR = TO H0.
            IF ((NOLIFA == 0) .OR. (NOLIFA == 1)) THEN
                DO I=1,NP
                    IF (DP2(I) < H0) DP2(I) = H0
                ENDDO
            ENDIF

#ifdef CMPI
            IF (MYPROC == 0) then
                WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
                WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            ENDIF
#else
            WRITE(ScreenUnit,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
            WRITE(16,'(A36,1X,E15.8,1X,A5)') 'BATHYMETRY RECORDS UPDATED AT TIME =',TIMELOC,'(SEC)'
#endif
        ENDIF  !test for updating bathymetry records

    !.......If time is during the bathymetry change interval, then update bathymetry
        IF(timeloc < BTIME_END) THEN !tcm 20150728 changed <= to < 
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
    
        IF(timeloc == BTIME_END) THEN
#ifdef CMPI
            IF (MYPROC == 0) then
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

    ENDIF  !NDDT = 2


!     2DDI.For the 2DDI version of the code
    if(C2DDI) then
    !     2DDI.Set up the 2D friction coefficient
        CALL Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, IFNLFA, &
        NP, TK)
    
    endif
!..RJW 3D for the 3D version of the code
    if(C3D) then
    !     3D.Set up the 3D friction coefficient
        CALL Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, IFNLFA, &
        NP, TK, NFEN)
    
    endif


!...  SHIFT THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION TO PREVIOUS
!...  TIME STEPS.  ZERO OUT THE NEW SPECIFIED NORMAL FLOW BOUNDARY
!...  CONDITION
!...
    DO I=1,NVEL
        QN0(I)=QN1(I)
        QN1(I)=QN2(I)
        QN2(I)=0.D0
        EN0(I)=EN1(I)
        EN1(I)=EN2(I)
        EN2(I)=0.D0
    END DO
!...
!...  DEFINE Ramp FUNCTION FOR BOUNDARY ELEVATION FORCING, WIND AND PRESSURE
!.... FORCING AND TIDAL POTENTIAL FORCING
!...

!     jgf46.08 Calculate ramp functions.
!     jgf46.21 Modify to match behavior of 46.02
    SELECT CASE(NRamp)
    CASE(0)
    Ramp=1.0D0
    RampExtFlux=1.0D0
    RampIntFlux=1.0D0
    RampElev=1.0D0
    RampTip=1.0D0
    RampMete=1.0D0
    RampWRad=1.0D0
! Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
!    kmd48.33bc - ramp changes with baroclinic when timestep is changed
    CASE(1,2,3,4,5,6,7,8)
    Ramp=TANH((2.D0*TimeLoc/86400.D0)/DRamp)
    RampExtFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampExtFlux)
    RampIntFlux=TANH((2.D0*TimeLoc/86400.D0)/DRampIntFlux)
    RampElev=TANH((2.D0*TimeLoc/86400.D0)/DRampElev)
    RampTip=TANH((2.D0*TimeLoc/86400.D0)/DRampTip)
! Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
    RampMete=TANH((2.D0*(TimeLoc/86400.D0-DUnRampMete))/DRampMete)
    RampWRad=TANH((2.D0*TimeLoc/86400.D0)/DRampWRad)
    END SELECT

!     jgf46.21 If there is an external flux (i.e. river) boundary, turn
!     off all forcings except the river flux forcing for the duration of
!     the FluxSettlingTime. When the FluxSettlingTime has ended, turn
!     all forcings back on.
    IF (NRamp > 1) THEN
        IF(IT < (FluxSettlingIT+10)) THEN
            Ramp=0.0
            RampIntFlux=0.0
            RampElev=0.0
            RampTip=0.0
            RampMete=0.0
            RampWRad=0.0
        ELSE
            Ramp=TANH((2.D0*(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRamp)
            RampIntFlux=TANH((2.D0 &
            *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampIntFlux)
            RampElev=TANH((2.D0 &
            *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampElev)
            RampTip=TANH((2.D0 &
            *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampTip)
            RampMete=TANH((2.D0 &
            *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampMete)
            RampWRad=TANH((2.D0 &
            *(IT-FluxSettlingIT-10)*DTDP/86400.D0)/DRampWRad)
        ! Corbitt 1203022: Added Zach's Fix for Assigning a Start Time to Mete Ramping
            IF(NRamp == 8) then
                RampMete=TANH((2.D0*((((IT)*DTDP)/86400.D0)-DUnRampMete))/DRampMete)
            endif
        ENDIF
    ! jgf49.44: Cover the case where the ramp length is zero.
        IF (DRamp < 1.0e-6) Ramp = 1.0d0
        IF (DRampExtFlux < 1.0e-6) RampExtFlux = 1.0d0
        IF (DRampIntFlux < 1.0e-6) RampIntFlux = 1.0d0
        IF (DRampElev < 1.0e-6) RampElev = 1.0d0
        IF (DRampTip < 1.0e-6) RampTip = 1.0d0
        IF (DRampMete < 1.0e-6) RampMete = 1.0d0
        IF (DRampWRad < 1.0e-6) RampWRad = 1.0d0
    ELSE
    ! jgf49.44: Cover the case where the ramp length is zero.
        IF (DRamp < 1.0e-6) Ramp = 1.0d0
    ENDIF


!------------------------ICE FIELDS----------------------------------------
!...  UPDATE THE ICE CONCENTRATION FIELDS FROM UNIT 25,225,227
!...  TCM V49.64.01 ADDED THE ICE FIELDS SECTION

    IF (NCICE == 12) THEN
        IF(TimeLoc > CICE_TIME2) THEN
            CICE_TIME1 = CICE_TIME2
            CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
            DO I=1,NP
                CICE1(I)=CICE2(I)
            END DO
            CALL NCICE1_GET(CICE2,NP)
        ENDIF
    ENDIF


!------------------------MET FORCING---------------------------------------

!...  UPDATE THE WIND STRESS AND SURFACE PRESSURE AND READ IN NEW VALUES
!...  FROM UNIT 22.  APPLY Ramp FUNCTION.

!  tcm v49.16 20100617 -- Changed pressure ramping so that we apply the ramp
!     to the difference between the background pressure and the forced pressure,
!     then add that ramped difference back to the background pressure

!     No wind, radiation stress or atmospheric pressure forcings are used.
    IF(NWS == 1) THEN
        DO I=1,NP
            READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I) !read in
            WSX2(I)=RampMete*WSX2(I) ! apply ramp function
            WSY2(I)=RampMete*WSY2(I)
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*PR2(I)
            PRDIFF = RampMete*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) ! for met recording station output
            wvnyout(i)=WSY2(i)
        END DO
    ENDIF

!     Wind stress and atmospheric pressure are read in at all grid nodes
!     at a time interval that does not equal the model time
!     step. Interpolation in time is used to synchronize the wind and
!     pressure information with the model time step.
    IF(ABS(NWS) == 2) THEN
        IF(TimeLoc > WTIME2) THEN !determine if met file time incr. is exceeded
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
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
        END DO
    ENDIF

!     Wind velocity in US Navy Fleet Numeric format interpolated in
!     space onto the ADCIRC grid and in time to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from the wind velocity.
    IF(NWS == 3) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
            END DO
            CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED, &
            NP,NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS, &
            NSCREEN, ScreenUnit)
        ENDIF
        WTRatio=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     Wind velocity and atmospheric pressure are read in (PBL/JAG
!     format) at selected ADCIRC grid nodes. Interpolation in time is
!     used to synchronize the wind and pressure information with the
!     model time step. Garratt's formula is used to compute wind stress
!     from wind velocity.
    IF(ABS(NWS) == 4) THEN
        IF(TimeLoc > WTIME2) THEN
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
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     Wind velocity and atmospheric pressure are read in at all grid
!     nodes. Interpolation in time is used to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from wind velocity.
    IF(ABS(NWS) == 5) THEN
        IF(TimeLoc > WTIME2) THEN
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
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     Wind velocity and atmospheric pressure are read in for a
!     rectangular grid (either in Longitude, Latitude or Cartesian
!     coordinates, consistent with the grid coordinates) and
!     interpolated in space onto the ADCIRC grid and in time to
!     synchronize the wind and pressure information with the model time
!     step. Garratt's formula is used to compute wind stress from the
!     wind velocity.
    IF(NWS == 6) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
            END DO
            NWSGGWI=NWSGGWI+1
            CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT, &
            WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
        ENDIF
        WTRatio=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     jgf46.01 New option to read in surface wind stress and atmospheric
!     pressure for a rectangular grid (either in Longitude, Latitude or
!     Cartesian coordinates, consistent with the grid coordinates) and
!     interpolate in space onto the ADCIRC grid. Interpolation in time
!     is used to synchronize the wind and pressure information with the
!     model time step.
    IF(ABS(NWS) == 7) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
            END DO
            CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT, &
            WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
        ENDIF
        WTRatio=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete*WindX !apply ramp
            WSY2(I) = RampMete*WindY
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
        END DO
    ENDIF

!     jgf46.02 New option to read in hurricane locations and generate
!     generate hurricane winds from the Holland Wind Model.
    IF(ABS(NWS) == 8) THEN
        HollandTime = TimeLoc
        CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP, &
        ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
        DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        ENDDO
    ENDIF

!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of National Weather Service (NWS) Aviation (AVN) model
!     output files. Each AVN file is assumed to contain data on a
!     Gaussian longitude, latitude grid at a single time. Consecutive
!     files in the sequence are separated by N hours in time. Garratt's
!     formula is used to compute wind stress from the wind velocity.
    IF(NWS == 10) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
            END DO
            NWSGGWI=NWSGGWI+1
            CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0, &
            G,NWLON,NWLAT,WTIMINC)
        ENDIF
        WTRatio=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of stripped down National Weather Service (NWS) ETA 29km
!     model output files. Each ETA file is assumed to contain data on an
!     E grid for a single day (8 data sets, one every 3 hours, beginning
!     @ 03:00 and continuing through 24:00 of the given day). The wind
!     data is converted to an east-west, north-south coordinate system
!     inside ADCIRC. Garratt's formula is used to compute wind stress
!     from the wind velocity.
    IF(NWS == 11) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
            END DO
            IDSETFLG=IDSETFLG+1
            IF(IDSETFLG > 8) THEN
                NWSEGWI=NWSEGWI+1
                IDSETFLG=1
            ENDIF
            CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2, &
            NP,RHOWAT0,G)
        ENDIF
        WTRatio=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
        !     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = &
                WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!...sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
    IF(ABS(NWS) == 12) THEN
        IF(TimeLoc > WTIME2) THEN
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
        ! asey 110518: Enable Mark Powell's sector-based wind drag.
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
            ! asey 110518: Enable Mark Powell's sector-based wind drag.
                WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
        !  tcm v49.16 20100617
        !            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
        ! RBITT (TO OUTPUT WIND STRESS UNCOMMENT BELOW)
        !            WVNXOUT(I)=WSX2(I)
        !            WVNYOUT(I)=WSY2(I)

#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX / WindMultiplier * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY / WindMultiplier * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     RJW ffpl Merged:
!     rjw added nws = 19: asymmetric hurricane winds
    AsymmetricTime = TimeLoc
    IF(NWS == 19) THEN
        CALL NWS19GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,AsymmetricTime,ICS)
        DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
        !            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     jie added nws = 20: generalized asymmetric vortex model
    GeneralAsymTime = TimeLoc
    IF(NWS == 20) THEN
        CALL NWS20GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,GeneralAsymTime,ICS)
        DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*0.001293d0*WDragCo*WindY*WindMag
        !  tcm v49.16 20100617C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
        !            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!     jgf49.1001 Added NWS29 for embedding an asymmetric vortex from
!     NWS19 into an OWI basin scale met field from NWS12 (derived from NAM).
    IF(ABS(NWS) == 29) THEN
    ! bring in next set of OWI met data
        IF(TimeLoc > WTIME2) THEN
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
        CALL NWS19GET(SLAM,SFEA,vortexWVNX2,vortexWVNY2,vortexPRN2, &
        NP,AsymmetricTime, ICS)
        DO I=1,NP
        ! compute wind stress due to background met at this node
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)

        !...zc51.06 - Make sure SWAN sees blended winds
#ifdef CSWAN
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX
                SWAN_WY2(I,2) = WindY
            ENDIF
#endif
        ! jgf49.1001 NAM winds already contain wind reduction
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
        ! compute wind stress due to vortex met at this node
            WindX = vortexWVNX2(I)
            WindY = vortexWVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
            CALL getBlendFactor(I, SLAM, SFEA, bf)
        ! blend the wind stresses and barometric pressures
            WSX2(I) = bf*RampMete*0.001293d0*WDragCo*WindX*WindMag &
            +(1.0d0-bf)*WSX2(I)
            WSY2(I) = bf*RampMete*0.001293d0*WDragCo*WindY*WindMag &
            +(1.0d0-bf)*WSY2(I)
            PR2(I)= bf*RampMete*vortexPRN2(I) &
            +(1.0d0-bf)*PRN2(I)
            WVNXOUT(I)=bf*RampMete*WindX+(1.0d0-bf)*WVNXOUT(I)
            WVNYOUT(I)=bf*RampMete*WindY+(1.0d0-bf)*WVNYOUT(I)
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
        ! c...make sure swan sees blended winds
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = ((WindX*bf)+(1d0-bf)*SWAN_WX2(I,2)) &
                * WaveWindMultiplier
                SWAN_WY2(I,2) = ((WindY*bf)+(1d0-bf)*SWAN_WY2(I,2)) &
                * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF


! jgf50.38.05: Added NWS=15 for reading HWind data
    IF(ABS(NWS) == 15) THEN
        CALL NWS15GET(WVNX2,WVNY2,PRN2,timeloc)
        DO I=1,NP
            windx = wvnx2(i)
            windy = wvny2(i)
            windMag = sqrt(windx**2 + windy**2)
            WDragCo = WindDrag(windMag, i)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(i, WDragCo, &
                windMag, DP(I), ETA2(I), H0, G, &
                windx, windy)
                WindMag = SQRT(windx**2+windy**2)
                WDragCo = WindDrag(WindMag, i)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(i,windx,windy)
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = PRBCKGRND_MH2O + RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
#ifdef CSWAN
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF
!.... tcm v51.06.02 added for GFDL Met Data
    IF(ABS(NWS) == 16) THEN
        IF(TimeLoc > WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
            END DO
            CALL NWS16GET(timeloc,WVNX2,WVNY2,PRN2)
        ENDIF
        WTRATIO=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
        ! asey 110518: Enable Mark Powell's sector-based wind drag.
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo, &
                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
            ! asey 110518: Enable Mark Powell's sector-based wind drag.
                WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef) &
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
        !    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE /= 0) THEN
                CICE_TRatio = (TimeLoc-CICE_TIME1)/CICE_TIMINC
                PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
                WDragCo = WindIceDrag(WDragCo,PIC)
                CICEOUT(I) = PIC
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I))) &
            -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
#ifdef CSWAN
        ! Casey 090302: Added these lines for coupling winds to SWAN.
        ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
            IF(COUPWIND)THEN
                SWAN_WX2(I,2) = WindX / WindMultiplier * WaveWindMultiplier
                SWAN_WY2(I,2) = WindY / WindMultiplier * WaveWindMultiplier
            ENDIF
#endif
        END DO
    ENDIF

!    kmd48.33bc - added in information for the sponge layer
    IF ((C3D) .AND. (RES_BC_FLAG > 0) .AND. (NWS /= 0)) THEN
        DO I=1,NP
            WSX2(I)=sponge(I)*WSX2(I)
            WSY2(I)=sponge(I)*WSY2(I)
        END DO
    END IF


!--------------------END MET FORCING---------------------------------------



!...  UPDATE THE WAVE RADIATION STRESS AND READ IN NEW VALUES FROM
!.... UNIT 23.  APPLY Ramp FUNCTION.  ADD RADIATION STRESS TO WIND
!...  STRESS
!...
!...  NRS=2 was added.  sb46.28sb03 09/xx/2006
!...  TCM v49.48 Restructured the wave stress updates in order to
!...             include NRS=4
    IF(NRS /= 0) THEN
        if((NRS == 1) .OR. (nrs == 2) .OR. (nrs == 3)) then
            IF(TimeLoc > RSTIME2) THEN
                RSTIME1=RSTIME2
                RSTIME2=RSTIME2+RSTIMINC
                DO I=1,NP
                    RSNX1(I)=RSNX2(I)
                    RSNY1(I)=RSNY2(I)
                END DO
                IF(NRS == 1) THEN
                    CALL RSGET(RSNX2,RSNY2)
                ENDIF
                IF(NRS == 2) THEN
                    CALL RS2GET(RSNX2,RSNY2,NP)
                ENDIF
#ifdef CSWAN
            ! Casey 090302: Added for coupling to SWAN.
                IF(NRS == 3) THEN
                    InterpoWeight = 1.0
                    CALL ComputeWaveDrivenForces
                ! Casey 090707: We want to extrapolate forward in time.  Load the latest (current) forces
                !             into RSNX1/RSNY1, and then load the future forces into RSNX2/RSNY2.
                    DO I=1,NP
                        RSX = RSNX1(I)
                        RSY = RSNY1(I)
                        RSNX1(I) = RSNX2(I)
                        RSNY1(I) = RSNY2(I)
                        RSNX2(I) = RSNX2(I) + (RSNX2(I)-RSX)
                        RSNY2(I) = RSNY2(I) + (RSNY2(I)-RSY)
                    ENDDO
#ifdef CSWANFRIC
                ! Casey 091020: Adopt Ethan's/Joannes's modified friction.  Compute the wave properties
                !             for height, period, angle and dissipation.  Then extrapolate forward
                !             in time using the same logic as above.
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
            ! ifdef CSWAN
            ! Casey 090302: Added these lines for output to the rads.64 file.
                IF(ABS(NRS) == 3) then
                    RSNXOUT(I) = RSX
                    RSNYOUT(I) = RSY
                ENDIF
#ifdef CSWANFRIC
            ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                WAVE_A(I) = WAVE_A1(I) + RSTRATIO*(WAVE_A2(I)-WAVE_A1(I))
                WAVE_H(I) = WAVE_H1(I) + RSTRATIO*(WAVE_H2(I)-WAVE_H1(I))
                WAVE_T(I) = WAVE_T1(I) + RSTRATIO*(WAVE_T2(I)-WAVE_T1(I))
#endif
            ! endif
            ENDDO
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            CALL ComputeModifiedWaveFriction(TK)
#endif
#endif
        ENDIF  !nrs = 1,2,or 3

    ! Tightly Coupled Code with STWAVE
    ! Apply the ramping function and add wave stress to WSX2,WSY2
    ! This cases uses a step function in time
        IF(NRS == 4) THEN    ! vjp modified Jan 22 2010
            IF(TimeLoc >= RSTIME2) THEN   !Get a new wave record
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

!     jgf48.4627 Skip past GWCE and momentum calculations if only
!     meteorological output was requested.
    IF (METONLY) THEN
        goto 9999
    ENDIF

!   kmd48.33 - added in information for the elevation boundary conditions
!              used in prognostic runs. The diagnostic information is read
!              in once and used during the simulation. Note that there is
!              no ramp utilized for this boundary condition.
    IF((C3D) .AND. (RES_BC_FLAG > 0) .AND. (CBAROCLINIC)) THEN
        IF ((ABS(RES_BC_FLAG) >= 1) .AND. (NOPE > 0)) THEN
            IF(TimeLoc > RBCTIME2) THEN
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
                LNM_BC(NumofBCNode)=LNM_BC1(NumofBCNode)+ &
                RBCRATIO*(LNM_BC2(NumofBCNode)- &
                LNM_BC1(NumofBCNode))
            END DO
        END IF
    END IF

!...
!...  Tidal Potential Forcing
!...  Note, the Earth tide potential reduction factor, ETRF(J) has been
!...        incorporated into this calculation.
!...
    IF(CTIP) THEN
        DO J=1,NTIF
            IF(PERT(J) == 0.) THEN
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
            IF(NA == 1) THEN    !SEMI-DIURNAL SPECIES
                DO I=1,NP
                    ARGTP=ARGT+2.d0*SLAM(I)
                    ARGSALT=ARGT-SALTPHA(J,I)
                    CCSFEA=COS(SFEA(I))
                    CCSFEA=CCSFEA*CCSFEA
                    TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP) &
                    +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
                END DO
            ENDIF
            IF(NA == 2) THEN    !DIURNAL SPECIES
                DO I=1,NP
                    ARGTP=ARGT+SLAM(I)
                    ARGSALT=ARGT-SALTPHA(J,I)
#ifdef REAL4
                    S2SFEA=SIN(2.e0*SFEA(I))
#else
                    S2SFEA=SIN(2.d0*SFEA(I))
#endif
                    TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP) &
                    +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
                END DO
            ENDIF
        END DO
    ENDIF

!...
!...  Depth Averaged Baroclinic Forcing needed by GWCE and 2DDI Momentum
!...  Compute this (divided by H, i.e., Bx/H, By/H) as a nodally averaged
!...  quantity for smoothing
!...
    DO J=1,NP
        VIDBCPDXOH(J)=0.D0
        VIDBCPDYOH(J)=0.D0
    ! jgf45.06 TotalArea(J)=0.D0 ! What if CBaroclinic=F and nolifa=0?
    ENDDO

    IF(CBaroclinic) Then

        DO J=1,NP              !jgf45.06 try this instead
            TotalArea(J)=0.D0   !jgf45.06
        ENDDO                  !jgf45.06

    !     Kendra45.12 - Test placement of BPG
        if (C3DVS) CALL BPG3D()

    !     jgf45.12 this algorithm only works in 2D
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
                DRhoDX=DARhoMRho0N1*FDX1O2A+DARhoMRho0N2*FDX2O2A &
                +DARhoMRho0N3*FDX3O2A
                DRhoDY=DARhoMRho0N1*FDY1O2A+DARhoMRho0N2*FDY2O2A &
                +DARhoMRho0N3*FDY3O2A
                VIDBCPDXOHN1=Ramp*G* &
                (DARhoMRho0N1*DEta2DX+0.5d0*H2N1*DRhoDX)
                VIDBCPDXOHN2=Ramp*G* &
                (DARhoMRho0N2*DEta2DX+0.5d0*H2N2*DRhoDX)
                VIDBCPDXOHN3=Ramp*G* &
                (DARhoMRho0N3*DEta2DX+0.5d0*H2N3*DRhoDX)
                VIDBCPDYOHN1=Ramp*G* &
                (DARhoMRho0N1*DEta2DY+0.5d0*H2N1*DRhoDY)
                VIDBCPDYOHN2=Ramp*G* &
                (DARhoMRho0N2*DEta2DY+0.5d0*H2N2*DRhoDY)
                VIDBCPDYOHN3=Ramp*G* &
                (DARhoMRho0N3*DEta2DY+0.5d0*H2N3*DRhoDY)
                VIDBCPDXOHAvgArea=AreaEle*(VIDBCPDXOHN1+VIDBCPDXOHN2 &
                +VIDBCPDXOHN3)/3.D0
                VIDBCPDYOHAvgArea=AreaEle*(VIDBCPDYOHN1+VIDBCPDYOHN2 &
                +VIDBCPDYOHN3)/3.D0
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
                IF(TotalArea(J) /= 0.) THEN
                    VIDBCPDXOH(J)=VIDBCPDXOH(J)/TotalArea(J)
                    VIDBCPDYOH(J)=VIDBCPDYOH(J)/TotalArea(J)
                ENDIF
            ENDDO
        ENDIF
    ENDIF

!...
!...  COMPUTE SPECIFIED NORMAL FLOW BOUNDARY CONDITION
!...
    IF(NFLUXF == 1) THEN
        IF (NFFR > 0) THEN
            DO J=1,NFFR
                IF(FPER(J) == 0.D0) THEN
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
                    IF(LBCODEI(J) == 32) THEN
                        ARG=ARGJ-ENPH(J,I)
                        EN2(I)=EN2(I)+ENAM(J,I)*RFF*COS(ARG)
                    ENDIF
                END DO
            END DO
        END IF
        IF((NFFR == 0) .OR. (NFFR == -1)) THEN
            IF(TimeLoc > QTIME2) THEN
                QTIME1=QTIME2
                QTIME2=QTIME2+FTIMINC
                DO J=1,NVEL
                    IF((LBCODEI(J) == 2) .OR. (LBCODEI(J) == 12) &
                     .OR. (LBCODEI(J) == 22)) THEN
                        QNIN1(J)=QNIN2(J)
                        READ(20,*) QNIN2(J)
                    ENDIF
                    IF(LBCODEI(J) == 32) THEN
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

    ! AL_add_42.06f
    !     jgf46.21 Collect elevation information for river radiation b.c.
        IF(IT == FluxSettlingIT) THEN
            EtaDisc_Fill = .FALSE.    ! sb v46.48 11/06/2006
            DO I=1, NP
                EtaDisc(I) = Eta2(I)   ! EtaDisc written to hotstart file
            ENDDO
            DO I=1,NVEL
                IF(LBCODEI(I) == 52) THEN
                    NNBB=NBV(I)
                    ElevDisc(I)=Eta2(NNBB)
                ENDIF
            END DO
        ELSE IF(EtaDisc_Fill .AND. IT > FluxSettlingIT) THEN
            EtaDisc_Fill = .FALSE. 
            DO I=1,NVEL
                IF(LBCODEI(I) == 52) THEN
                    NNBB=NBV(I)
                    ElevDisc(I)=EtaDisc(NNBB)   ! sb v46.48 11/06/2006
                ENDIF
            END DO
        ENDIF

    ENDIF
!...
!...  COMPUTE DISCHARGE CONTRIBUTION FROM RADIATION BOUNDARY CONDITION
!...
    IF(NFLUXRBC == 1) THEN
        DO J=1,NVEL
            IF(LBCODEI(J) == 30) THEN
                NNBB=NBV(J)
                H1=DP(NNBB)+IFNLFA*ETA2(NNBB)
                UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
                QN1(J)=H1*UN1
            ENDIF
        END DO
    ENDIF

!...  COMPUTE DISCHARGE CONTRIBUTION FROM ZERO NORMAL VELOCITY GRADIENT
!...  BOUNDARY CONDITION
!...
    IF(NFLUXGBC == 1) THEN
        DO J=1,NVEL
            IF((LBCODEI(J) == 40) .OR. (LBCODEI(J) == 41)) THEN
                NNBB=NBV(J)
                H1=DP(NNBB)+IFNLFA*ETA2(NNBB)
                UN1=UU1(NNBB)*CSII(J)+VV1(NNBB)*SIII(J)
                QN1(J)=H1*UN1
            ENDIF
        END DO
    ENDIF
!...
!...  COMPUTE SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
!.... EXTERNAL BARRIER BOUNDARY NODES
!...  COBELL - MOVED TO WEIR_BOUNDARY.F
    IF(NFLUXB == 1) THEN
        DO I=1,NVEL
            SELECT CASE(LBCODEI(I))
            CASE(3,13,23)
            CALL COMPUTE_EXTERNAL_BOUNDARY_FLUX(I,TIMELOC,QN2(I))
            END SELECT
        END DO
    ENDIF

!...  COMPUTE INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED INTERNAL BARRIER
!...  BOUNDARY (PERMEABLE OR NOT) NODES
!...
!     jgf46.03 Begin block of notes for internal barrier boundaries

!     NFLUXIB is set to 1 in read_input.F if there are internal barrier
!     boundaries in the fort.14 (mesh) file.

!     IBSTART is a flag that indicates the first time through the time
!     stepping loop; set to 0 in read_input.F and set to 1 here.

!     BARAVGWT was apparently intended for use in averaging internal
!     barrier water levels. It is set to 0 in read_input.F, which seems
!     to turn off any time averaging here.

!     NIBNODECODE seems to be set to 1 for nodes receiving water across
!     the barrier

!     BARMIN is used in several places, mainly as the minimum elevation
!     above the levee for flow to occur. It is a parameter and is set to
!     0.04 in global.F.


!...ZC - SIMPLIFIED THIS SECTION, MOVED FLUX COMPUTATION TO WEIR_BOUNDARY.F
!   USE THE COMPILER FLAGS TO CHANGE IMPLEMENTATION:
!        -DORIGWEIR - JOANNES WESTERINK ET AL IMPLEMENTATION (ORIGINAL)
!                     THIS DOES NOT APPEAR TO HAVE BEEN USED IN QUITE
!                     SOME TIME
!         DEFAULT   - SHINTARO BUNYA/SEIZO TANAKA IMPLEMENTATION FOR
!                     CHECKING FOR A WET EDGE BEFORE PASSING FLOW
!                     ACROSS A WEIR
    IF(NFLUXIB == 1) THEN
        NIBNODECODE(:) = 0
        I = 0
        DO K = 1, NBOU
            SELECT CASE(LBCODEI(I+1))
            CASE(4,24,5,25)
            DO J = 1,NVELL(K)*2
                I = I + 1
                CALL COMPUTE_INTERNAL_BOUNDARY_FLUX(I,J,K, &
                TIMELOC,QN2(I))
            ENDDO
            CASE DEFAULT
            I = I + NVELL(K)
            END SELECT
        ENDDO
    ENDIF

!...
!...  COMPUTE INWARD/OUTWARD NORMAL FLOW FOR INTERNAL BARRIER
!.... BOUNDARY NODES THROUGH CROSS BARRIER PIPES
!.... NOTE THAT THIS ADDS AN ADDITIONAL FLOW COMPONENT INTO QN2
!...
    IF(NFLUXIBP == 1) THEN
        DO I=1,NVEL
            IF((LBCODEI(I) == 5) .OR. (LBCODEI(I) == 25)) THEN
                PIPE_FLUX = 0D0
                CALL COMPUTE_CROSS_BARRIER_PIPE_FLUX(I,TIMELOC,PIPE_FLUX)
                QN2(I) = QN2(I) + PIPE_FLUX
            ENDIF
        ENDDO
    ENDIF

    if(subdomainOn .AND. enforceBN == 1) call readFort019(it)  ! NCSU Subdomain
    if(subdomainOn .AND. enforceBN == 2) call readFort020(it)  ! NCSU Subdomain
    if(subdomainOn .AND. enforceBN == 2) call readFort021(it)  ! NCSU Subdomain


    if(subdomainOn .AND. enforceBN == 1) call readFort019(it)  ! NCSU Subdomain
    if(subdomainOn .AND. enforceBN == 2) call readFort020(it)  ! NCSU Subdomain
    if(subdomainOn .AND. enforceBN == 2) call readFort021(it)  ! NCSU Subdomain


!...
!...  Compute the water surface elevation from the GWCE form of the 2D
!...  continuity eq.
!...

! md    Changed to include the predictor-corrector algorithm
    IF(CPRECOR) THEN
        CALL GWCE_New(IT,TimeLoc,TimeH)
    !...tcm added call for slope limiting
        IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
        ENDIF
        CALL Mom_Eqs_New_NC()
    !...  If running in parallel, update velocities & fluxes on all processors

#ifdef CMPI
        CALL UPDATER(UU2,VV2,DUMY1,2)
        CALL UPDATER(QX2,QY2,DUMY1,2)
#endif
        CALL GWCE_New_pc(IT,TimeLoc,TimeH)
    !...tcm added call to slope limiter
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
    !.... tcm added call to slope limiter
        IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
        ENDIF
    ENDIF

! ET...
! ET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
! ET...
! ET...NOTE:NNODECODE is a working variable that can change more than once
! ET...               during a time step
! ET...     NNODECODE = 0 for a dry node
! ET...     NNODECODE = 1 for a wet node
! ET...     NODECODE  - is a more static version of NNODECODE that is reconciled
! ET...                 once and for all at the end time step
! ET...
! ET...
! ET...        (   DRYING CRITERIA   )
! ET...
! ET...A node should be dry under two conditions.
! ET...D1.) If the total water depth falls below H0.
! ET .......Note: if the total water depth falls below 0.8*H0, the surface elevation
! ET........is lifted up so that the total water depth = 0.8*H0.
! ET......
! ET...D2.) If the node is connected to only nonfunctioning (dry) elements.  In
! ET........this case the node is dried due to becoming landlocked.
! ET........Note: this criteria is applied after all other wetting and drying criteria
! ET...
! ET...An element should be dry under the following conditions.
! ET...DE3.) This is an elemental check section designed to avoid artificial wetting of
! ET.........of control sections
! ET.........All elements where downhill flow originates from a barely wet node
! ET.........(defined as 1.2*H0) into wet nodes are forced inactive; the only exception
! ET......... is receiving overtopped barrier nodes
! ET...
! ET...        (   WETTING CRITERIA   )
! ET...
! ET...A node should be wet under two conditions.
! ET...W1.) If 2 nodes on an element are wet and one is dry, wet the dry node
! ET........if the water level at one of the wet nodes is greater than the
! ET........water level at the dry node and the steady state velocity that
! ET........would result from a balance between the water level gradient and
! ET........bottom friction would yield a velocity > VELMIN.
! ET........Note that the criteria outlined in DE3 must also be satified before
! ET........the node is allowed to wet
! ET...
! ET...W2.) If an element has a node lying on a receiving internal barrier boundary or
! ET......specified discharge boundary that is actively discharging flow into the
! ET......domain at that node, all nodes in this element must stay wet.
! ET...
! ET...
! ET...        (  VELOCITY BOUNDARY CONDITION  )
! ET...
! ET...Either a natural or essential boundary condition can be used as a velocity
! ET...boundary condition in the momentum equation solution along a wet/dry boudary
! ET...To use a natural boundary condition, do nothing along the wet/dry interface.
! ET...To use an essential, no velocity boundary condition, identify the nodes along
! ET...the wet/dry interface and zero out the velocity at the nodes.  Interface nodes
! ET...can easily be identified by comparing the number of active elements a node is
! ET...connected to (MJU) to the total number of elements a node is connected to (NODELE).
! ET...If MJU < NODELE for any node, it must lie along the wet/dry interface.  See
! ET...further comments at the end of the momentum equation solution section.
! ET...

! ET...
! ET...WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP
! ET...
    IF(NOLIFA == 2) THEN      !  This goes on until end of part 6
        DO I=1,NP
            NIBCNT(I) = 0
        ENDDO
        DO I=1,NE
            NOFFOLD(I)=NOFF(I)
            NOFF(I)=1
        ENDDO

    ! ET...
    ! ET...WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
    ! ET....Drying Criteria D1: this depends on NODECODE and updates NODECODE
    ! ET...
        DO I=1,NP
            IF(NODECODE(I) == 1) THEN
                HTOT=DP(I)+ETA2(I)
                IF(HTOT <= H0) THEN
                    IF(HTOT < HABSMIN) ETA2(I)=HABSMIN-DP(I)
                    NNODECODE(I)=0
                    NODECODE(I)=0
                    NCCHANGE=NCCHANGE+1 !NCCHANGE=0 set near beginning of GWCE
                !                  ENDIF
                ENDIF
            ENDIF
        ENDDO
    ! ET...
    ! ET...END WET/DRY SECTION - PART 1
    ! ET...

    ! jwC     Use Message-Passing to update nodecode and nnodecode at ghost nodes
    ! jw#ifdef CMPI
    ! jw         CALL UPDATEI(NODECODE,NNODECODE,2)
    ! jw#endif

    ! ET...
    ! ET...WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
    ! ET...
        DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)

        ! ET...
        ! ET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
        ! ET...
            NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
            IF(NCTOT == 2) THEN
                ETAN1=ETA2(NM1)
                ETAN2=ETA2(NM2)
                ETAN3=ETA2(NM3)
                HTOTN1=DP(NM1)+ETA2(NM1)
                HTOTN2=DP(NM2)+ETA2(NM2)
                HTOTN3=DP(NM3)+ETA2(NM3)
                IF((NODECODE(NM1) == 1) .AND. (NODECODE(NM2) == 1)) THEN
                    IF((HTOTN1 >= HOFF) .AND. (HTOTN2 >= HOFF)) THEN
                        NM123=NM1
                        IF(ETA2(NM2) > ETA2(NM1)) NM123=NM2
                        DELDIST=SQRT((y(NM3)-y(NM123))**2.D0 &
                        +(X(NM3)-X(NM123))**2.D0)
                        DELETA=ETA2(NM123)-ETA2(NM3)
                    ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA < 0.d0) DELETA = 0.d0
                        H1=ETA2(NM123)+DP(NM123)
                    !. RJW merged from Casey 071219: Added the following logic for 3D friction.
                    !. RJW modified the following for 3D friction
                        IF(C2DDI)THEN
                        ! b46.28sb02
                        !<<                     Convert Manning's N to Cd, if necessary.
                            IF (LoadManningsN) THEN
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                **(1.d0/3.d0) )
                                IF(FRIC(NM123) < BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                ENDIF
                            ENDIF
                        !>>
                            TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)* &
                            (IFNLBF+IFHYBF* &
                            (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))

                            IF(TKWET < 0.0001d0) TKWET=0.0001d0
                            VEL=G*(DELETA/DELDIST)/TKWET

                        ELSEIF(C3D)THEN
                        ! solve for the depth averaged velocity,U, from the relation :
                        !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                        !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                        ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                            IF(LoadZ0B_var) THEN
                                Z0B1 = Z0B_var(NM123)
                            ELSEIF (LoadManningsN) THEN
                                Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+ &
                                ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/ &
                                (ManningsN(NM123)*sqrt(g)) ) ))
                            ELSE
                                Z0B1 = Z0B
                            ENDIF
                            VEL=sqrt(g*H1*(DELETA/DELDIST)) &
                            * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                        ENDIF

                        IF(VEL > VELMIN) THEN
                        !    ....         third node met criteria and is also wet
                            NNODECODE(NM3)=1
                        !. RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                            IF(C2DDI)THEN
                            !                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
                                TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                (IFNLBF+IFHYBF* &
                                (1.D0+(HBREAK/H1)**FTHETA)** &
                                (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                            ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                TKXX(NM123) = TK(NM123)
                                TKYY(NM123) = TK(NM123)
                                TKXY(NM123) = 0.D0
#endif
#endif
                            ELSEIF(C3D)THEN
                                IF(ISLIP == 0)THEN
                                    DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                    /(SIGMA(2)-SIGMA(1))
                                    BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM3)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                                ENDIF
                                IF(ISLIP /= 0)THEN
                                    IF(ISLIP == 1)THEN
                                        KSLIP=KP
                                    ENDIF
                                    IF(ISLIP == 2)THEN
                                        KSLIP = (1.D0 / &
                                        ( (1.D0/0.41D0) * &
                                        LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                        (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                        +Z0B1) &
                                        / (Z0B1) ) ) )**2.D0 &
                                        * ABS(Q(NM123,1))
                                        IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                        KSLIP=1.D0* ABS(Q(NM123,1))
                                        IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                        KSLIP=0.0025D0*ABS(Q(NM123,1))
                                    ENDIF
                                    BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                    BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                    BSX1(NM3)=KSLIP*REAL(Q(NM123,1))
                                    BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ELSEIF((NODECODE(NM2) == 1) .AND. (NODECODE(NM3) == 1)) &
                    THEN
                    IF((HTOTN2 >= HOFF) .AND. (HTOTN3 >= HOFF)) THEN
                        NM123=NM2
                        IF(ETA2(NM3) > ETA2(NM2)) NM123=NM3
                        DELDIST=SQRT((Y(NM1)-Y(NM123))**2.D0 &
                        +(X(NM1)-X(NM123))**2.D0)
                        DELETA=ETA2(NM123)-ETA2(NM1)
                    ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA < 0.d0) DELETA = 0.d0
                        H1=ETA2(NM123)+DP(NM123)
                    !. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                        IF(C2DDI)THEN
                        ! b46.28sb02
                        !<<                     Convert Manning's N to Cd, if necessary.
                            IF (LoadManningsN) THEN
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                **(1.d0/3.d0) )
                                IF(FRIC(NM123) < BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                ENDIF
                            ENDIF
                        !>>
                            TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)* &
                            (IFNLBF+IFHYBF* &
                            (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                            IF(TKWET < 0.0001d0) TKWET=0.0001d0
                            VEL=G*(DELETA/DELDIST)/TKWET

                        ELSEIF(C3D)THEN
                        ! solve for the depth averaged velocity,U, from the relation :
                        !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                        !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                        ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                            IF(LoadZ0B_var) THEN
                                Z0B1=Z0B_var(NM123)
                            ELSEIF (LoadManningsN) THEN
                                Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+ &
                                ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/ &
                                (ManningsN(NM123)*sqrt(g)) ) ))
                            ELSE
                                Z0B1=Z0B
                            ENDIF
                            VEL=sqrt(g*H1*(DELETA/DELDIST)) &
                            * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                        ENDIF

                        IF(VEL > VELMIN) THEN
                            NNODECODE(NM1)=1
                        !. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                            IF(C2DDI)THEN
                                TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                (IFNLBF+IFHYBF* &
                                (1.D0+(HBREAK/H1)**FTHETA)** &
                                (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                            ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                TKXX(NM123) = TK(NM123)
                                TKYY(NM123) = TK(NM123)
                                TKXY(NM123) = 0.D0
#endif
#endif
                            ELSEIF(C3D)THEN
                                IF(ISLIP == 0)THEN
                                    DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                    /(SIGMA(2)-SIGMA(1))
                                    BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM1)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                                ENDIF
                                IF(ISLIP /= 0)THEN
                                    IF(ISLIP == 1)THEN
                                        KSLIP=KP
                                    ENDIF
                                    IF(ISLIP == 2)THEN
                                        KSLIP = (1.D0 / &
                                        ( (1.D0/0.41D0) * &
                                        LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                        (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                        +Z0B1) &
                                        / (Z0B1) ) ) )**2.D0 &
                                        * ABS(Q(NM123,1))
                                        IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                        KSLIP=1.D0* ABS(Q(NM123,1))
                                        IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                        KSLIP=0.0025D0*ABS(Q(NM123,1))
                                    ENDIF
                                    BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                    BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                    BSX1(NM1)=KSLIP*REAL(Q(NM123,1))
                                    BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                                ENDIF
                            ENDIF

                        ENDIF
                    ENDIF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ELSEIF((NODECODE(NM3) == 1) .AND. (NODECODE(NM1) == 1)) &
                    THEN
                    IF((HTOTN3 >= HOFF) .AND. (HTOTN1 >= HOFF)) THEN
                        NM123=NM3
                        IF(ETA2(NM1) > ETA2(NM3)) NM123=NM1
                        DELDIST=SQRT((Y(NM2)-Y(NM123))**2.D0 &
                        +(X(NM2)-X(NM123))**2.D0)
                        DELETA=ETA2(NM123)-ETA2(NM2)
                    ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA < 0.d0) DELETA = 0.d0
                        H1=ETA2(NM123)+DP(NM123)
                    !. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                        IF(C2DDI)THEN
                        ! b46.28sb02
                        !<<                     Convert Manning's N to Cd, if necessary.
                            IF (LoadManningsN) THEN
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                **(1.d0/3.d0) )
                                IF(FRIC(NM123) < BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                ENDIF
                            ENDIF
                        !>>
                            TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)* &
                            (IFNLBF+IFHYBF* &
                            (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                            IF(TKWET < 0.0001d0) TKWET=0.0001d0
                            VEL=G*(DELETA/DELDIST)/TKWET
                        ELSEIF(C3D)THEN
                        ! solve for the depth averaged velocity,U, from the relation :
                        !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                        !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                        ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                            IF(LoadZ0B_var) THEN
                                Z0B1 = Z0B_var(NM123)
                            ELSEIF (LoadManningsN) THEN
                                Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+ &
                                ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/ &
                                (ManningsN(NM123)*sqrt(g)) ) ))
                            ELSE
                                Z0B1 = Z0B
                            ENDIF
                            VEL=sqrt(g*H1*(DELETA/DELDIST)) &
                            * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                        ENDIF

                        IF(VEL > VELMIN) THEN
                            NNODECODE(NM2)=1
                        !. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                            IF(C2DDI)THEN
                                TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                (IFNLBF+IFHYBF* &
                                (1.D0+(HBREAK/H1)**FTHETA)** &
                                (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                            ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                TKXX(NM123) = TK(NM123)
                                TKYY(NM123) = TK(NM123)
                                TKXY(NM123) = 0.D0
#endif
#endif
                            ELSEIF(C3D)THEN
                                IF(ISLIP == 0)THEN
                                    DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                    /(SIGMA(2)-SIGMA(1))
                                    BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM2)=EVTOT(1)*REAL(DUDS)
                                    BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                                ENDIF
                                IF(ISLIP /= 0)THEN
                                    IF(ISLIP == 1)THEN
                                        KSLIP=KP
                                    ENDIF
                                    IF(ISLIP == 2)THEN
                                        KSLIP = (1.D0 / &
                                        ( (1.D0/0.41D0) * &
                                        LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                        (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                        +Z0B1) &
                                        / (Z0B1) ) ) )**2.D0 &
                                        * ABS(Q(NM123,1))
                                        IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                        KSLIP=1.D0* ABS(Q(NM123,1))
                                        IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                        KSLIP=0.0025D0*ABS(Q(NM123,1))
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
        ! ET...
        ! ET...Nodal Wetting Criteria W2a
        ! ET...
            NBNCTOT=NIBNODECODE(NM1)+NIBNODECODE(NM2)+NIBNODECODE(NM3)
            NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
            NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
            NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

        ENDDO

        if (subdomainOn .AND. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
        if (subdomainOn .AND. enforceBN == 2) call enforceWDob() ! NCSU Subdomain

    !     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
        CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif


    ! et...
    ! ET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
    !*******************************************************************************************
    ! tcm v50.66.01 -- This is an additional test for wetting only when time varying
    !                  bathymetry is used and is only performed during the period of
    !                  bathymetry evolution.
    
        IF ((NDDT /= 0) .AND. (IT <= BTIME_END+1) ) THEN
            DO I=1,NE
                NM1=NM(I,1)
                NM2=NM(I,2)
                NM3=NM(I,3)
                NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
                IF(NCTOT < 3) THEN   !If not wet from previous time step
                    NCTOT=NNODECODE(NM1)+NNODECODE(NM2)+NNODECODE(NM3)
                    if(NCTOT < 3) then !if not alreay made wet for this time step
                        ETAN1=ETA2(NM1)
                        ETAN2=ETA2(NM2)
                        ETAN3=ETA2(NM3)
                        HTOTN1=DP(NM1)+ETA2(NM1)
                        HTOTN2=DP(NM2)+ETA2(NM2)
                        HTOTN3=DP(NM3)+ETA2(NM3)

                    !                    if all nodes have a depth greater than or equal to
                    !                    hoff = 1.2*H0, then make the element wet
                        IF( (HTOTN1 >= HOFF) .AND. (HTOTN2 >= HOFF) .AND. &
                        (HTOTN3 >= HOFF) ) THEN
                        ! HE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                        !  Make Node 1 Wet and set parameters
                            IF(NNODECODE(NM1) /= 1) THEN  !node 1
                                NNODECODE(NM1)=1
                                NM123 = NM1
                                IF(C2DDI)THEN
                                !<<                           Convert Manning's N to Cd, if necessary.
                                    IF (LoadManningsN) THEN
                                        FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                        /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                        **(1.d0/3.d0) )
                                        IF(FRIC(NM123) < BFCdLLimit) THEN
                                            FRIC(NM123) = BFCdLLimit
                                        ENDIF
                                    ENDIF
                                ENDIF
                                VEL=VELMIN
                                H1 = HTOTN1
                                IF(C2DDI)THEN
                                    TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                    (IFNLBF+IFHYBF* &
                                    (1.D0+(HBREAK/H1)**FTHETA)** &
                                    (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                                ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                    TKXX(NM123) = TK(NM123)
                                    TKYY(NM123) = TK(NM123)
                                    TKXY(NM123) = 0.D0
#endif
#endif
                                ELSEIF(C3D)THEN
                                    IF(ISLIP == 0)THEN
                                        DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                        /(SIGMA(2)-SIGMA(1))
                                        BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                        BSX1(NM3)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                                    ENDIF
                                    IF(ISLIP /= 0)THEN
                                        IF(ISLIP == 1)THEN
                                            KSLIP=KP
                                        ENDIF
                                        IF(ISLIP == 2)THEN
                                            KSLIP = (1.D0 / &
                                            ( (1.D0/0.41D0) * &
                                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                            (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                            +Z0B1) &
                                            / (Z0B1) ) ) )**2.D0 &
                                            * ABS(Q(NM123,1))
                                            IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                            KSLIP=1.D0* ABS(Q(NM123,1))
                                            IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                            KSLIP=0.0025D0*ABS(Q(NM123,1))
                                        ENDIF
                                        BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                        BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                        BSX1(NM3)=KSLIP*REAL(Q(NM123,1))
                                        BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                                    ENDIF
                                ENDIF
                            ENDIF  !end node 1

                        !  Make Node 2 Wet and set parameters
                            IF (NNODECODE(NM2) /= 1) THEN
                                NNODECODE(NM2) = 1
                                NM123=NM2
                                IF(C2DDI)THEN
                                !<<                        Convert Manning's N to Cd, if necessary.
                                    IF (LoadManningsN) THEN
                                        FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                        /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                        **(1.d0/3.d0) )
                                        IF(FRIC(NM123) < BFCdLLimit) THEN
                                            FRIC(NM123) = BFCdLLimit
                                        ENDIF
                                    ENDIF
                                ENDIF
                                VEL = VELMIN
                                H1 = HTOTN2
                                IF(C2DDI)THEN
                                    TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                    (IFNLBF+IFHYBF* &
                                    (1.D0+(HBREAK/H1)**FTHETA)** &
                                    (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                                ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                    TKXX(NM123) = TK(NM123)
                                    TKYY(NM123) = TK(NM123)
                                    TKXY(NM123) = 0.D0
#endif
#endif
                                ELSEIF(C3D)THEN
                                    IF(ISLIP == 0)THEN
                                        DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                        /(SIGMA(2)-SIGMA(1))
                                        BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                        BSX1(NM1)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                                    ENDIF
                                    IF(ISLIP /= 0)THEN
                                        IF(ISLIP == 1)THEN
                                            KSLIP=KP
                                        ENDIF
                                        IF(ISLIP == 2)THEN
                                            KSLIP = (1.D0 / &
                                            ( (1.D0/0.41D0) * &
                                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                            (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                            +Z0B1) &
                                            / (Z0B1) ) ) )**2.D0 &
                                            * ABS(Q(NM123,1))
                                            IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                            KSLIP=1.D0* ABS(Q(NM123,1))
                                            IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                            KSLIP=0.0025D0*ABS(Q(NM123,1))
                                        ENDIF
                                        BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                        BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                        BSX1(NM1)=KSLIP*REAL(Q(NM123,1))
                                        BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                                    ENDIF
                                ENDIF
                            ENDIF  !node 2

                        !  Make Node 3 Wet and set parameters
                            IF(NNODECODE(NM3) /= 1) THEN
                                NNODECODE(NM3)=1
                                NM123 = NM3
                                IF(C2DDI)THEN
                                !<<                         Convert Manning's N to Cd, if necessary.
                                    IF (LoadManningsN) THEN
                                        FRIC(NM123)=g*ManningsN(NM123)**2.d0 &
                                        /( ( DP(NM123)+IFNLFA*ETA2(NM123) ) &
                                        **(1.d0/3.d0) )
                                        IF(FRIC(NM123) < BFCdLLimit) THEN
                                            FRIC(NM123) = BFCdLLimit
                                        ENDIF
                                    ENDIF
                                ENDIF
                                VEL=VELMIN
                                H1 = HTOTN3
                                IF(C2DDI)THEN
                                    TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)* &
                                    (IFNLBF+IFHYBF* &
                                    (1.D0+(HBREAK/H1)**FTHETA)** &
                                    (FGAMMA/FTHETA)))
#ifdef CSWAN
#ifdef CSWANFRIC
                                ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
                                    TKXX(NM123) = TK(NM123)
                                    TKYY(NM123) = TK(NM123)
                                    TKXY(NM123) = 0.D0
#endif
#endif
                                ELSEIF(C3D)THEN
                                    IF(ISLIP == 0)THEN
                                        DUDS=(Q(NM123,2)-Q(NM123,1)) &
                                        /(SIGMA(2)-SIGMA(1))
                                        BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                        BSX1(NM2)=EVTOT(1)*REAL(DUDS)
                                        BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                                    ENDIF
                                    IF(ISLIP /= 0)THEN
                                        IF(ISLIP == 1)THEN
                                            KSLIP=KP
                                        ENDIF
                                        IF(ISLIP == 2)THEN
                                            KSLIP = (1.D0 / &
                                            ( (1.D0/0.41D0) * &
                                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))* &
                                            (DP(NM123)+IFNLFA*ETA2(NM123))) &
                                            +Z0B1) &
                                            / (Z0B1) ) ) )**2.D0 &
                                            * ABS(Q(NM123,1))
                                            IF(KSLIP > 1.D0*ABS(Q(NM123,1))) &
                                            KSLIP=1.D0* ABS(Q(NM123,1))
                                            IF(KSLIP < 0.0025D0*ABS(Q(NM123,1))) &
                                            KSLIP=0.0025D0*ABS(Q(NM123,1))
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

            if (subdomainOn .AND. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn .AND. enforceBN == 2) call enforceWDob() ! NCSU Subdomain

        !     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
            CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif
        
        ENDIF !IT TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME
    
    !.... END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY
    ! ET..
    ! ET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
    !*******************************************************************************************

    ! ET...
    ! ET...Nodal Wetting Criteria W2b
    ! ET...Check for adjacent nodes and force nodes wet when attached
    ! ET...to receiving barrier nodes
    ! ET...
        DO I=1,NP
            IF((NIBCNT(I) > 0) .AND. (NNODECODE(I) == 0)) THEN
                NNODECODE(I)=1
            ENDIF
        ENDDO

    ! jwC     Use Message-Passing to update nnodecode at ghost nodes
    ! jw#ifdef CMPI
    ! jw         CALL UPDATEI(NNODECODE,IDUMY,1)
    ! jw#endif

    ! ET...
    ! ET...END WET/DRY SECTION - PART 2
    ! ET...

    ! ET...
    ! ET...START WET/DRY SECTION  - PART 3
    ! ET...Elemental drying criteria DE1
    ! ET...This is an elemental check section designed to avoid artificial wetting of
    ! ET....of control sections
    ! ET...All elements where downhill flow originates from a barely wet node
    ! ET....into wet nodes are forced inactive; the only exception is receiving
    ! ET....overtopped barrier nodes
    ! ET...
        DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NBNCTOT=NIBCNT(NM1)*NIBCNT(NM2)*NIBCNT(NM3)
            IF(NBNCTOT == 0) THEN   !No barrier/pipe receiving nodes in this elem
                ETAN1=ETA2(NM1)
                ETAN2=ETA2(NM2)
                ETAN3=ETA2(NM3)
                HTOTN1=DP(NM1)+ETA2(NM1)
                HTOTN2=DP(NM2)+ETA2(NM2)
                HTOTN3=DP(NM3)+ETA2(NM3)
#ifdef SB_WETDRY
            !...Find the heighest point on the bed in the element.  sb v46.28.sb05.06 11/01/2006
                IF(DP(NM1) <= DP(NM2) .AND. DP(NM1) <= DP(NM3)) THEN
                    DPMIN = DP(NM1)
                ELSE IF(DP(NM2) <= DP(NM3) .AND. DP(NM2) <= DP(NM1)) THEN
                    DPMIN = DP(NM2)
                ELSE IF(DP(NM3) <= DP(NM1) .AND. DP(NM3) <= DP(NM2)) THEN
                    DPMIN = DP(NM3)
                ENDIF
#endif
#ifndef SB_WETDRY
            ! jgf52.08.08: Analyst can eliminate noff from
            ! consideration in fort.15 file.
                if (noffActive.eqv. .TRUE. ) then
                !...ABC pattern
                    IF((ETAN1 >= ETAN2) .AND. (ETAN2 > ETAN3)) THEN
                        IF((HTOTN1 < HOFF) .OR. (HTOTN2 < HOFF)) NOFF(I)=0
                    ENDIF
                    IF((ETAN2 >= ETAN3) .AND. (ETAN3 > ETAN1)) THEN
                        IF((HTOTN2 < HOFF) .OR. (HTOTN3 < HOFF)) NOFF(I)=0
                    ENDIF
                    IF((ETAN3 >= ETAN1) .AND. (ETAN1 > ETAN2)) THEN
                        IF((HTOTN3 < HOFF) .OR. (HTOTN1 < HOFF)) NOFF(I)=0
                    ENDIF
                !...ACB pattern
                    IF((ETAN1 >= ETAN3) .AND. (ETAN3 > ETAN2)) THEN
                        IF((HTOTN1 < HOFF) .OR. (HTOTN3 < HOFF)) NOFF(I)=0
                    ENDIF
                    IF((ETAN2 >= ETAN1) .AND. (ETAN1 > ETAN3)) THEN
                        IF((HTOTN2 < HOFF) .OR. (HTOTN1 < HOFF)) NOFF(I)=0
                    ENDIF
                    IF((ETAN3 >= ETAN2) .AND. (ETAN2 > ETAN1)) THEN
                        IF((HTOTN3 < HOFF) .OR. (HTOTN2 < HOFF)) NOFF(I)=0
                    ENDIF
                endif
#else
            ! jgf52.08.08: Analyst can eliminate noff from
            ! consideration in fort.15 file.
                if (noffActive.eqv. .TRUE. ) then
                !...ABC pattern
                    IF((ETAN1 >= ETAN2) .AND. (ETAN2 > ETAN3)) THEN
                        IF((HTOTN1 < HOFF)) NOFF(I)=0
                        IF((ETAN1-ETAN2) < (ETAN2-ETAN3)) THEN
                            IF(HTOTN2 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                    IF((ETAN2 >= ETAN3) .AND. (ETAN3 > ETAN1)) THEN
                        IF((HTOTN2 < HOFF)) NOFF(I)=0
                        IF((ETAN2-ETAN3) < (ETAN3-ETAN1)) THEN
                            IF(HTOTN3 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                    IF((ETAN3 >= ETAN1) .AND. (ETAN1 > ETAN2)) THEN
                        IF((HTOTN3 < HOFF)) NOFF(I)=0
                        IF((ETAN3-ETAN1) < (ETAN1-ETAN2)) THEN
                            IF(HTOTN1 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                !...ACB pattern
                    IF((ETAN1 >= ETAN3) .AND. (ETAN3 > ETAN2)) THEN
                        IF((HTOTN1 < HOFF)) NOFF(I)=0
                        IF((ETAN1-ETAN3) < (ETAN3-ETAN2)) THEN
                            IF(HTOTN3 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                    IF((ETAN2 >= ETAN1) .AND. (ETAN1 > ETAN3)) THEN
                        IF((HTOTN2 < HOFF)) NOFF(I)=0
                        IF((ETAN2-ETAN1) < (ETAN1-ETAN3)) THEN
                            IF(HTOTN1 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                    IF((ETAN3 >= ETAN2) .AND. (ETAN2 > ETAN1)) THEN
                        IF((HTOTN3 < HOFF)) NOFF(I)=0
                        IF((ETAN3-ETAN2) < (ETAN2-ETAN1)) THEN
                            IF(HTOTN2 < HOFF) NOFF(I)=0
                        ENDIF
                    ENDIF
                endif !noffActive.eqv. .TRUE. 
#endif

#ifdef SB_WETDRY
            !...An element is set to be dry if it is determined to be a flooding type
            !...wetting element.  An element is a flooding type wetting element if
            !...the bed elevation at the node with the biggest water column height
            !...is lower than the heighest point on the bed in the element.
            !...sb v46.28.sb05.06 11/01/2006
            !...This is applied only when NOFF flag of the element at the previous time step
            !...is 0, which means that this logic works to prevent an element from
            !...re-wetting.
            !...sb v46.52.03
                if (noffActive.eqv. .TRUE. ) then
                    IF(NOFFOLD(I) == 0) THEN
                        IF(HTOTN1 >= HTOTN2 .AND. HTOTN1 >= HTOTN3) THEN
                            IF(ETAN1 < (-DPMIN+H0)) NOFF(I) = 0
                        ENDIF
                        IF(HTOTN2 >= HTOTN3 .AND. HTOTN2 >= HTOTN1) THEN
                            IF(ETAN2 < (-DPMIN+H0)) NOFF(I) = 0
                        ENDIF
                        IF(HTOTN3 >= HTOTN1 .AND. HTOTN3 >= HTOTN2) THEN
                            IF(ETAN3 < (-DPMIN+H0)) NOFF(I) = 0
                        ENDIF
                    ENDIF
                endif
#endif
            ENDIF
        ENDDO

#ifdef SB_WETDRY
    ! ET......added by sb on 11/02/2006
    ! ET...
    ! ET...This section is added after we realize that it's not possible
    ! ET...to compute a correct flow going through two elements if
    ! ET...the elements are connected just by one node. i.e., elements need
    ! ET...to share an edge to let the flow go through between the elements.
    ! ET...Therefore, in this section, a node is determined to be dry
    ! ET...if two elements are connected at one node, not sharing an edge.
    ! ET...As it seemed this procedure needed NOFF information,
    ! ET...although this section changes NNODECODE,
    ! ET...I put this section here, rather than the end of PART 2.
    ! ET...I tried setting NOFF(I) = 0, but it didn't shut down the flow.
    ! ET...
        DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            IF(NOFF(I) == 1 .AND. &
            NNODECODE(NM1) == 1 .AND. &
            NNODECODE(NM2) == 1 .AND. &
            NNODECODE(NM3) == 1) THEN
                DO K=1,3
                    NM1=NM(I,K)
                    NM2=NM(I,MOD(K+0,3)+1)
                    NM3=NM(I,MOD(K+1,3)+1)

                    NWETNEI = 0
                    NWETADJ = 0
                    DO J=1,MNEI
                        N=NeiTabEle(NM1,J)
                        IF(N == 0) CYCLE
                        IF(N == I) CYCLE

                        NMN1=NM(N,1)
                        NMN2=NM(N,2)
                        NMN3=NM(N,3)
                        IF(NOFF(N) == 1 .AND. &
                        NNODECODE(NMN1) == 1 .AND. &
                        NNODECODE(NMN2) == 1 .AND. &
                        NNODECODE(NMN3) == 1) THEN
                            NWETNEI = NWETNEI + 1

                            IF((NMN1 == NM2 .OR. NMN1 == NM3) .OR. &
                            (NMN2 == NM2 .OR. NMN2 == NM3) .OR. &
                            (NMN3 == NM2 .OR. NMN3 == NM3)) THEN
                                NWETADJ = NWETADJ + 1
                            ENDIF
                        ENDIF
                    ENDDO

                    IF(NWETNEI > 0 .AND. NWETADJ == 0 .AND. &
                    NIBCNT(NM1) == 0) THEN
                        NNODECODE(NM1) = 0
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
#endif

    ! ET...
    ! ET...END WET/DRY SECTION  - PART 3
    ! ET...

    ! ET...
    ! ET...START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2
    ! ET...Update number of active elements (MJU) and the total area (TotalArea) connected
    ! ET...to a node. If these are zero, the node is landlocked and should be dried.
    ! ET...These depend on NNODECODE which varies during the time step
    ! ET...
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

    ! jwnote - looks like this is used later in momentum equations
    ! jwnote - this has implications on making this into a subroutine

        DO I=1,NP
            IF((NNODECODE(I) == 1) .AND. (MJU(I) == 0)) THEN
                NNODECODE(I)=0
            ENDIF
            IF(MJU(I) == 0) MJU(I)=1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
        ENDDO

    !     WET...
    !     WET...END WET/DRY SECTION - PART 4
    !     WET...

    ! jwnote - may have to pass TotalArea and mju as well

        if (subdomainOn .AND. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
        if (subdomainOn .AND. enforceBN == 2) call enforceWDob() ! NCSU Subdomain


    !     Use Message-Passing to update nnodecode at ghost nodes
#ifdef CMPI
        CALL UPDATEI(NNODECODE,IDUMY,1)
#endif
    ! ET...
    ! ET...WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE
    ! ET...Check to see if any wetting occurred & update NODECODE
    ! ET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
    ! ET...
        DO I=1,NP
            IF(NNODECODE(I) /= NODECODE(I)) THEN
                NODECODE(I)=NNODECODE(I)
                NCCHANGE=NCCHANGE+1
            ENDIF
        ENDDO
    ! ET...
    ! ET...END WET/DRY SECTION - PART 5
    ! ET...

    ! ET...
    ! ET...WET/DRY SECTION - PART 6
    ! ET...Check to see if any NOFF changed requiring the matrix to be reset
    ! ET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
    ! ET...
        DO I=1,NE
            IF(NOFF(I) /= NOFFOLD(I)) NCCHANGE=NCCHANGE+1
        ENDDO
    ! ET...
    ! ET... jgf45.06 If there has been any wetting or drying in any
    ! ET... of the subdomains, the NCCHANGE flag will be activated on all
    ! ET... of the subdomains, to prevent them from getting out of sync
    ! ET... with their MPI calls as some reset the GWCE and others do not.
    ! ET...
#ifdef CMPI
    ! jgf48.4619 implementing Seizo's changes for Lumped, fully
    ! explicit operation. In that case, the GWCE LHS matrix is
    ! recalculated on each individual subdomain that has wetted
    ! or dried, without recourse to MPI, eliminating the need
    ! for the call to the subroutine WetDrySum.
        IF ( ILump == 0 ) THEN
            call WetDrySum(NCCHANGE)
        ELSE
            NCCHANGE=NCCHANGE ! jgf48.4619 do nothing
        ENDIF
#endif
    ! ET...
    ! ET...END WET/DRY SECTION - PART 6
    ! ET...
    ENDIF                     !  This is started in Part 1 of CWET
!.....

!...
!...  2DDI Momentum Equation Solution
!...
    IF (C2DDI) THEN
        IF (CME_Orig) THEN
            CALL Mom_Eqs_Original()
        ENDIF
        IF (CME_New_NC) THEN
            CALL Mom_Eqs_New_NC()
        ENDIF
        IF ((CME_New_C1) .OR. (CME_New_C2)) THEN
            CALL Mom_Eqs_New_Conserv()
        ENDIF
        IF (CPRECOR) THEN
            CALL Mom_Eqs_Non_Conserv_pc()
        ENDIF

    !...  If running in parallel, update velocities & fluxes on all processors

#ifdef CMPI
        CALL UPDATER(UU2,VV2,DUMY1,2)
        CALL UPDATER(QX2,QY2,DUMY1,2)
#endif

    ENDIF

!...
!     C2DDI....END OF 2DDI MOMENTUM EQUATION SOLUTION
!...

!...
!...  3DVS Momentum Equation Solution
!...
    IF (C3DVS) THEN

    !... Load the vector MOM_LV_X(I) with barotropic pressure terms
    !...     including atmospheric pressure, water level and tidal potential
    !...     averaged between time levels s and s+1, (time levels 1 and 2).  Note:
    !...     MOM_LV_X gets renamed as BTP in global_3dvs.f

        DO I=1,NP
            MOM_LV_X(I)=ETA1(I)+ETA2(I)
            IF(NWS /= 0) MOM_LV_X(I)=MOM_LV_X(I)+PR1(I)+PR2(I) !atmospheric pressure
            IF (CTIP) MOM_LV_X(I)=MOM_LV_X(I)-TIP1(I)-TIP2(I) !tidal potential
            MOM_LV_X(I)=G*MOM_LV_X(I)/2.d0
        ENDDO

    !...  Solve for velocity at the new time level (K+1)

        CALL VSSOL(IT,TimeLoc)

    ENDIF
!...
!...  End of 3DVS Momentum Equation Solution
!...

!...
!...  IF 2D TRANSPORT IS INCLUDED SOLVE FOR THE CONCENTRATION
!...
    IF(C2D_PTrans) THEN

        CALL SCALAR_TRANS_2D (IT, TimeLoc)

    ENDIF
!...
!...  End of 2D Scalar Transport Solution
!...
!     jgf48.4627 Jump to here if METONLY is .TRUE., i.e., if only
!     meteorological output is requested.
    9999 CONTINUE
!...
!...  Collect maximum values of variables   v46.50 sb 11/11/2006
    call collectMinMaxData(timeloc)

! jgf52.08.01: If detailed inundation output was requested,
! collect detailed inundation data.
    if (inundationOutput.eqv. .TRUE. ) then
        call collectInundationData(timeloc, it)
    endif

!...
!...  WRITE OUTPUT
!...
! ALL writeOutput2D(IT,TimeLoc) ! =>zc - moved to adcirc.F so that SWAN writes
!        correctly

! NCSU Subdomain Modeling
    if (subdomainOn .AND. NOUTGS == 1) call writeFort065(it)
    if (subdomainOn .AND. NOUTGS == 2) call writeFort066(it)
    if (subdomainOn .AND. NOUTGS == 2) call writeFort067(it)

!     jgf49.44: If harmonic analysis was requested and the current time
!     is within the harmonic analysis period, update the left hand side
!     of the harmonic analysis matrix. Also update the load vectors for
!     each type of analysis. If timeseries reconstruction was specified,
!     also update the timeseries.
    CALL updateHarmonicAnalysis(IT, TIMEH)
!...
!...  WRITE OUT HOT START INFORMATION IF NHSTAR=1 AND AT CORRECT TIME
!.... STEP
!...  NOTE: THE HOT START FILES USE A RECORD LENGTH OF 8 ON BOTH 32 BIT
!.... WORKSTATIONS AND THE 64 BIT CRAY.  THIS IS BECAUSE THE HARMONIC
!.... ANALYSIS IS DONE IN DOUBLE PRECISION (64 BITS) ON WORKSTATIONS.
!...
    ITEST=(IT/NHSINC)*NHSINC

!     IF(myproc.eq.0) PRINT *, " ITEST **********", ITEST

    if ((ABS(NHSTAR) > 0 .AND. ITEST == IT) .OR. (-IHOT == IT)) then  !tcm v51.26 added abs(nhstar) to handle when nhstar = -1
        if( (MNWPROH > 0) ) then   !Writer for HSfile
            if( ( .NOT. C3D) ) then                     !Writer for HSfile
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
    ! Casey 100205: Enable writing of SWAN hot-start file.  We need to wait
    !             and do this after the next SWAN time step, so that
    !             everything is up-to-date.
        WriteSwanHotStart = .TRUE. 
#endif
    ENDIF



!...  SAVE THE CURRENT TIME LEVEL OF BARRIER INTO THE PREVIOUS TIME LEVEL
    IF(NFLUXIB == 1)THEN
        BARINHT1(:)=BARINHT2(:)
    ENDIF
    IF(NFLUXB == 1)THEN
        BARLANHT1(:)=BARLANHT2(:)
    ENDIF

!...
!...  find and print to unit 6, the maximum elevation, the maximum
!...  velocity and the node numbers at which they occur on myproc=0 if
!...  elmax exceeds threshold, print information on all processors where
!...  this occurs
!...

!     jgf46.00 Added option to output data to the screen every NSCREEN
!     time steps, rather than on every time step, as long as there are
!     no high elevations. In the case of high elevations, the warning
!     messages are sent to the screen each time they are generated.

    IF(NSCREEN /= 0) THEN
        ELMAX=0.0d0
        VELMAX=0.0d0
        KEMAX = 0
        KVMAX = 0
        DO I=1,NP
            IF((NODECODE(I) == 1) .AND. (ABS(ETA2(I)) > ELMAX))THEN
                ELMAX=ABS(ETA2(I))
                KEMAX=I
            ENDIF
            VELABS=UU2(I)*UU2(I)+VV2(I)*VV2(I)
            IF (VELABS > VELMAX) THEN
                VELMAX=VELABS
                KVMAX=I
            ENDIF
        END DO
        VELMAX=VELMAX**0.5d0
        ITEST=(IT/NSCREEN)*NSCREEN

    !     jgf46.10 Added the ability to for the user to control the warning
    !     and error elevations. Also added the ability for the user to write
    !     a fort.69 (global elevation debug) file.

    !     jgf46.12 Removed the dependence on KEMAX for producing output to
    !     the screen.
#ifdef CMPI
        IF(MYPROC == 0 .AND. ELMAX < WarnElev .AND. ITEST == IT) THEN
            IF (KEMAX > 0) THEN
                WRITE(ScreenUnit,1991) &
                IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,MYPROC
            ELSE
                WRITE(ScreenUnit,1991) &
                IT,NUMITR,TimeLoc,0.,KEMAX,VELMAX,KVMAX,MYPROC
            ENDIF
            1991 FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8, &
            &            2X,'ON MYPROC = ',I4)
        ENDIF
        WarnElevExceeded = 0
        IF(ELMAX > WarnElev) THEN
            WRITE(ScreenUnit,1993) &
            IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX, &
            KVMAX,MYPROC
            WRITE(16,1993) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX, &
            KVMAX,MYPROC
            1993 FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8, &
            &            2X,'ON MYPROC = ',I4, &
            &            3X,'** WARNING: Elevation > WarnElev **')
            IF (WarnElevDump) WarnElevExceeded=1
        ENDIF
#ifdef DEBUG_WARN_ELEV
        call WarnElevSum(WarnElevExceeded)
        IF (WarnElevExceeded /= 0) THEN
            CALL WriteWarnElev(TimeLoc, IT)
        ENDIF
#endif
        ErrorElevExceeded = 0                ! Clint's Zombie Slyaer
        IF(ELMAX > ErrorElev) THEN
            WRITE(ScreenUnit,1995) &
            IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX, &
            KVMAX,MYPROC
            WRITE(16,1995) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX, &
            KVMAX,MYPROC
            1995 FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8, &
            &            2X,'ON MYPROC = ',I4,/, &
            &            2X,'** ERROR: Elevation > ErrorElev,', &
            ' ADCIRC stopping. **')
            ErrorElevExceeded = 1                ! Clint's Zombie Slayer
        ENDIF
        call WarnElevSum(ErrorElevExceeded)  ! Clint's Zombie Slayer 2010.08.07
        IF( ErrorElevExceeded /= 0 ) THEN    !  Finalize MPI Environment,
            Flag_ElevError = .TRUE. 
            CALL MSG_FINI()                    !  if there are Error Elvation nodes.
            STOP                               !
        ENDIF                                !st3
#else
        IF(ELMAX < WarnElev .AND. ITEST == IT) THEN
            IF (KEMAX > 0) THEN
                WRITE(ScreenUnit,1992) &
                IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            ELSE
                WRITE(ScreenUnit,1992) &
                IT,NUMITR,TimeLoc,0.,KEMAX,VELMAX,KVMAX
            ENDIF
            1992 FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8)
        ENDIF
        IF(ELMAX > WarnElev) THEN
            WRITE(ScreenUnit,1994) &
            IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1994) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            1994 FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8, &
            &            2X,'** WARNING: Elevation > WarnElev **')
#ifdef DEBUG_WARN_ELEV
            IF (WarnElevDump) CALL WriteWarnElev(TimeLoc, IT)
#endif
        ENDIF
        IF(ELMAX > ErrorElev) THEN
            Flag_ElevError = .TRUE. 
            WRITE(ScreenUnit,1996) &
            IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            WRITE(16,1996) IT,NUMITR,TimeLoc,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
            1996 FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5, &
            &            5X,'TIME = ',E15.8, &
            /,2X,'ELMAX = ', 1pE12.4E3,' AT NODE ',I8, &
            &            2X,'SPEEDMAX = ',1pE12.4E3,' AT NODE ',I8,/, &
            &            2X,'** ERROR: Elevation > ErrorElev, ' &
            'ADCIRC stopping. **')
            CALL ADCIRC_Terminate()
        ENDIF
#endif
    ENDIF

!...
!...  ****************** TIME STEPPING LOOP ENDS HERE ********************
!...

!. RJW merged 08/26/2008 Casey 071219: Added the folowing call to the mass balance subroutine.
!      IF(C3DVS)THEN
!         CALL MASSBAL3D(IT)
!      ENDIF

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN

!******************************************************************************
    END SUBROUTINE TIMESTEP
!******************************************************************************

!******************************************************************************
!                                                                             *
!    Subroutine to compute the elevation using the GWCE formluation           *
!    Re-written to conform to the ADCIRC Theory Report                        *
!                                                                             *
!                            r.l.  06/22/2005                                 *
!******************************************************************************

    SUBROUTINE GWCE_New(IT,TimeLoc,TimeH)

#ifdef IEEE_DEBUG
    USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
    USE SIZES
    USE GLOBAL
    USE MESH, ONLY : NE, NP, NM, X, Y, DP, NNeigh, NeiTab, TotalArea, &
    Areas, NEIMAX, SFAC, nneighele, neitabele
    USE BOUNDARIES, ONLY : NETA, NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, &
    NFLUXRBC, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3
    USE WIND
    USE ITPACKV
    USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
    USE NodalAttributes, ONLY : &
    LoadGeoidOffset, GeoidOffset, EVM, &
    TAU0VAR, HighResTimeVaryingTau0, FullDomainTimeVaryingTau0, &
    CalculateTimeVaryingTau0, &
    LoadEleSlopeLim,elemental_slope_limiter_active, &
    LoadAdvectionState, &
    elemental_slope_limiter_grad_max, &
    elemental_slope_limiter_max_exceeded

#ifdef CMPI
    USE MESSENGER
#endif
#ifdef CSWAN
#ifdef CSWANFRIC
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
    USE Couple2Swan, ONLY: TKXX, &
    TKXY, &
    TKYY
#endif
#endif
    USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceECB, &
    checkChange, enforceEOB, enforceEIB, enforceGWCELVOB

    IMPLICIT NONE

    INTEGER :: IE, JN, IJ, I, J                           !local loop counters
    INTEGER :: IT
    INTEGER :: NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI, NCJ
    INTEGER :: NCyc
    INTEGER :: NBDI
    INTEGER :: OnDiag, OffDiag

    LOGICAL ::  DIE
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

!     jgf50.44: Added water surface elevation smoothing
    REAL(sz) :: EtaN1, ETaN2, EtaN3, AreaEle, EtaN123
    INTEGER :: IT_SmoothTime
!      moved to subroutine -- tcm
!      REAL(sz), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
!      LOGICAL, SAVE :: firstCall = .true.

    REAL(SZ) HH1 !jgf46.02 Added for Katrina.
    REAL(SZ) RDIAG ! jgf48.4619 Seizo parameter for fully explicit mode
    REAL(SZ) dEta2Mag ! magnitude of the slope of the water surface elevation
    REAL(SZ) DEta2DX,DEta2DY
    REAL(SZ) DRhoDX,DRhoDY

    call setMessageSource("gwce_new")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!   tcm moved to subroutine
!      IF (LoadEleSlopeLim.eqv..true.) THEN
!         IF (firstCall.eqv..true.) THEN
!            allocate(elevSum(mnp))
!            firstCall = .false.
!         ENDIF
!      ENDIF

    if (subdomainOn .AND. enforceBN == 2) call checkChange()    ! NCSU Subdomain

! Casey 050711 : Added for averaged variable Tau0.
!      REAL(SZ) :: CaseySum
!      REAL(SZ), ALLOCATABLE :: TAU0VARTEMP(:)

!     jgf45.11 Bug fix: calculate the integers OnDiag and OffDiag here
!     instead of inside the GWCE lhs (system matrix) setup, since they
!     are also used in the calculation of the GWCE load vector gwce_lv.

!...  Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1
    OnDiag=(1+ILump)*2        !diagonal coefficient
    OffDiag=(1-ILump)         !off diagonal coefficient
!...
!...  Recompute the GWCE system matrix at the first time step or if any
!...  wetting or drying occurred in the previous time step.
!...
    IF(NCChange > 0) THEN !if any subdomain grid has changed
        NCChange=0
    !        IF(NScreen.GT.0.AND.MYProc.EQ.0) WRITE(screenunit,3806)
    !        WRITE(16,3806)
        3806 FORMAT(/,1X,'RE-SETTING GWCE SYSTEM MATRIX',/)

    !.....Set up the LHS matrix (for the iterative matrix solver)
        IF ( ILump == 0 ) THEN ! default, fully consistent case
            Coef(:,:)=0.0d0
        ELSE ! jgf48.4619: ILump == 1, only need the diagonals (Seizo)
            Coefd(:)=0.0d0    ! Only Diagnal
        ENDIF
    
    !        jgf47.08 Moved time-varying tau0 subroutine to nodalattr.F
        IF(C2DDI .AND. &
        ((FullDomainTimeVaryingTau0 .OR. HighResTimeVaryingTau0))) &
        THEN
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
                TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
        ENDIF

        IF(C3D .AND. &
        ((FullDomainTimeVaryingTau0 .OR. HighResTimeVaryingTau0))) &
        THEN
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
                TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
        ENDIF

    ! jgf48.4619: Added Seizo's handling of Lumping vs non-Lumping
        IF ( ILump == 0 ) THEN ! default case: fully consistent LHS
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
                    IF(NeiTab(NMI1,JN) == NMJ2) J12=JN
                    IF(NeiTab(NMI1,JN) == NMJ3) J13=JN
                    IF(NeiTab(NMI2,JN) == NMJ1) J21=JN
                    IF(NeiTab(NMI2,JN) == NMJ3) J23=JN
                    IF(NeiTab(NMI3,JN) == NMJ1) J31=JN
                    IF(NeiTab(NMI3,JN) == NMJ2) J32=JN
                END DO

                Coef(NMI1,1)  =Coef(NMI1,1)   + (MsFacLOnDiag &
                +GA00DPAvgOAreaIE4*(FDX1*FDX1+FDY1*FDY1))*NCELE
                Coef(NMI1,J12)=Coef(NMI1,J12) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX1*FDX2+FDY1*FDY2))*NCELE
                Coef(NMI1,J13)=Coef(NMI1,J13) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX1*FDX3+FDY1*FDY3))*NCELE
                Coef(NMI2,J21)=Coef(NMI2,J21) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX2*FDX1+FDY2*FDY1))*NCELE
                Coef(NMI2,1)  =Coef(NMI2,1)   + (MsFacLOnDiag &
                +GA00DPAvgOAreaIE4*(FDX2*FDX2+FDY2*FDY2))*NCELE
                Coef(NMI2,J23)=Coef(NMI2,J23) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX2*FDX3+FDY2*FDY3))*NCELE
                Coef(NMI3,J31)=Coef(NMI3,J31) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX3*FDX1+FDY3*FDY1))*NCELE
                Coef(NMI3,J32)=Coef(NMI3,J32) + (MsFacLOffDiag &
                +GA00DPAvgOAreaIE4*(FDX3*FDX2+FDY3*FDY2))*NCELE
                Coef(NMI3,1)  =Coef(NMI3,1)   + (MsFacLOnDiag &
                +GA00DPAvgOAreaIE4*(FDX3*FDX3+FDY3*FDY3))*NCELE

            ENDDO
        ELSE
        ! jgf48.4619: Add Seizo's construction of Lumped LHS matrix
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

    !...  Modify the matrix "COEF" by imposing the elevation specified
    !...  boundary conditions while maintaining the symmetry of the system

        IF (ILump == 0) THEN
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
        !...        for each elevation specified boundary node, zero all off diagonal
        !...        terms on the row and set diagnoal term to EP
            DO I=1,NETA
                Coef(NBD(I),1)=EP
                DO J=2,NNEIGH(NBD(I))
                    Coef(NBD(I),J)=0.0d0
                ENDDO
            ENDDO
        !...        for each elevation specified boundary node, zero all off diagonal
        !...        terms on the column but save these to be multiplied by the
        !...        boundary value and subtracted from the RHS
            DO I=1,NETA
                DO J=2,NNeigh(NBD(I))
                    DO IJ=2,NNeigh(NeiTab(NBD(I),J))
                        IF(NBD(I) == NeiTab(NeiTab(NBD(I),J),IJ)) THEN
                            OBCCoef(I,J-1)=Coef(NeiTab(NBD(I),J),IJ)
                            Coef(NeiTab(NBD(I),J),IJ)=0.0d0
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        !.....      Check that all the diagonal elements in "COEF" are > 0.
            DIE = .FALSE. 
            DO I=1,NP
                IF(COEF(I,1) == 0.d0) COEF(I,1)=EP
                IF(COEF(I,1) < 0.d0) THEN
                    IF(NSCREEN /= 0 .AND. MYPROC == 0) &
                    WRITE(ScreenUnit,1019) I,COEF(I,1)
                    WRITE(16,1019) I,COEF(I,1)
                    DIE = .TRUE. 
                ENDIF
            ENDDO
            IF (DIE) THEN
            ! jgfdebug
                open(899,file='debug.txt',status='replace',action='write')
                do i=1,np
                    if (coef(i,1) < 0.d0) then
                        do j=2,nneighele(i)
                            if (neitabele(i,j) /= 0) then
                                write(6, &
                                '("Node ",i0," element ",i0," area=",f15.7)') i, j, &
                                areas(neitabele(i,j))
                            endif
                        enddo
                    endif
                enddo
                close(899)
            ! end jgfdebug
                CALL ADCIRC_Terminate()
            ENDIF

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

        !.....      Check that all the diagonal elements in "COEFD" are > 0.
            DIE = .FALSE. 
            DO I=1,NP
                IF(COEFD(I) == 0.d0) COEFD(I)=EP
                IF(COEFD(I) < 0.d0) THEN
                    IF(NSCREEN /= 0 .AND. MYPROC == 0) &
                    WRITE(ScreenUnit,1019) I,COEFD(I)
                    WRITE(16,1019) I,COEFD(I)
                    1019 FORMAT(/,1X,'!!!!!!!!  ERROR  !!!!!!!', &
                    /,1X,'THE DIAGONAL TERM IN THE EQUATION FOR NODE ',I10, &
                    '= ',E15.6,' AND IS <= 0',/)
                    DIE = .TRUE. 
                ENDIF
            ENDDO
            IF (DIE) CALL ADCIRC_Terminate()
        ENDIF
    ENDIF                     !End of GWCE matrix setup
!...
!...  Compute the GWCE load vector GWCE_LV
!...  This is done primarily element by element by forming
!...  temporary vectors and then assembling at the end.
!...  This has been set up to unroll loops to optimize performance
!...  on vector processors.
!...
!...  Elevation and flux boundary conditions are imposed after the
!...  element by element assembly section.
!...

!...  Initialize variables to zero if these forcings are not used

    IF((NWS /= 0) .OR. (NRS /= 0)) THEN
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

!...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

    IF ((CGWCE_LS_2PartV) .OR. (CGWCE_LS_2PartSV)) THEN

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
                EVMSmag=EVMEle* &
                sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA) &
                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
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
            IF(TotalArea(I) /= 0.) THEN
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

!...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

    IF ((CGWCE_LS_2PartQ) .OR. (CGWCE_LS_2PartSQ)) THEN

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
                EVMSmag=EVMEle* &
                sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA) &
                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
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
            IF(TotalArea(I) /= 0.) THEN
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

!...  Assemble the GWCE RHS except for the boundary integral terms

    DO 1037 IE=1,NE

    !...     Set nodal values for each element

    ! Corbitt 120322: Localized Advection
        IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

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

        IF((NWS /= 0) .OR. (NRS /= 0)) THEN     !wind or radiation stress
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
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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

    !...     Compute part of several spatial gradients for use below

        E0XGrad2A = 0.0d0
        E0YGrad2A = 0.0d0
        IF (ILump == 0) THEN
            E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
            E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
        ENDIF
        E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
        E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
        Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
        Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

    !...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM

        IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.d0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3) &
            -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3) &
            +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3) &
            -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
        ENDIF

    !...     Compute the remainder of the 2 Part lateral stress terms

        IF((CGWCE_LS_2PartQ) .OR. (CGWCE_LS_2PartV) .OR. &
        (CGWCE_LS_2PartSQ) .OR. (CGWCE_LS_2PartSV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3 &
            +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3 &
            +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
        ENDIF

    !...     Compute the spatial gradients of the velocity dispersion terms if 3D

        IF (C3D) THEN                         !3D bottom friction
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3 &
            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3 &
            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
        ENDIF

    !...     Compute elemental averages

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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
        ENDIF

        IF (C3D) THEN                !3D velocity dispersion
            DispXAvg=IFNLCT*DispX
            DispYAvg=IFNLCT*DispY
        ENDIF
        IF(CBaroclinic) THEN
            BCXAvg=(H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2) &
            +H1N3*VIDBCPDXOH(NM3))/3.d0
            BCYAvg=(H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2) &
            +H1N3*VIDBCPDYOH(NM3))/3.d0
        ENDIF

    !...     Compute additional partial factors

        MsFacR=AreaIE*(1.d0/DT-Tau0Avg/2.d0)/DT/12.d0
        GOAreaIE4=G/AreaIE4
        Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
        A00pB00=A00+B00

    !...     Compute the JX, JY terms less the advection terms

        JXAvg = CorifAvg*QY1Avg &
        -IFNLFA*GOAreaIE4*(E1N1SQ*FDX1+E1N2SQ*FDX2 &
        +E1N3SQ*FDX3) &
        -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1 &
        +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3) &
        +WSXAvg-BSXAvg+MXAvg-DispXAvg-BCXAvg+Tau0QXAvg

        JYAvg =-CorifAvg*QX1Avg &
        -IFNLFA*GOAreaIE4*(E1N1SQ*FDY1+E1N2SQ*FDY2 &
        +E1N3SQ*FDY3) &
        -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1 &
        +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3) &
        +WSYAvg-BSYAvg+MYAvg-DispYAvg-BCYAvg+Tau0QYAvg

    !...     Complete the JX, JY terms depending on the advection formulation

        IF(CGWCE_Advec_NC) THEN        !nonconservative advection
            JXAvg = JXAvg - IFNLCT*( &
            QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3) &
            +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2 &
            +IFNLCAT*U1Avg*ESAvg/DT
            JYAvg = JYAvg - IFNLCT*( &
            QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3) &
            +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2 &
            +IFNLCAT*V1Avg*ESAvg/DT
        ENDIF
        IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
            JXAvg = JXAvg - IFNLCT*( &
            (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2 &
            +U1N3*QX1N3*FDX3) &
            +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2 &
            +U1N3*QY1N3*FDY3))/AreaIE2
            JYAvg = JYAvg - IFNLCT*( &
            (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2 &
            +V1N3*QX1N3*FDX3) &
            +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2 &
            +V1N3*QY1N3*FDY3))/AreaIE2
        ENDIF
        IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
            JXAvg = JXAvg - IFNLCT*( &
            QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3) &
            +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3) &
            +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
            JYAvg = JYAvg - IFNLCT*( &
            QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3) &
            +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3) &
            +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
        ENDIF

    !...     Assemble forcing for node NM1 (local index j=1)

        Temp_LV_A1= &

    !...     Transient and Tau0 terms from LHS
        (OnDiag*ESN1 + OffDiag*(ESN2+ESN3))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX1*E0XGrad2A+FDY1*E0YGrad2A) &
        +A00pB00*(FDX1*E1XGrad2A+FDY1*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX1 + JYAvg*FDY1)/2.d0 &

    !...     Spatially varying Tau0 terms
        +Tau0SpaVar

    !...     Assemble forcing for node NM2 (local index j=2)

        Temp_LV_A2= &

    !...     Transient and Tau0 terms from LHS
        (OnDiag*ESN2 + OffDiag*(ESN1+ESN3))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX2*E0XGrad2A+FDY2*E0YGrad2A) &
        +A00pB00*(FDX2*E1XGrad2A+FDY2*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX2 + JYAvg*FDY2)/2.d0 &

    !...     Spatially varying Tau0 terms
        +Tau0SpaVar


    !...     Assemble forcing for node NM3 (local index j=3)

        Temp_LV_A3= &

    !...     Transient and Tau0 terms from LHS
    !...    (consistent mass matrix: ILump=0, lumped mass matrix: ILump=1)
        (OnDiag*ESN3 + OffDiag*(ESN1+ESN2))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX3*E0XGrad2A+FDY3*E0YGrad2A) &
        +A00pB00*(FDX3*E1XGrad2A+FDY3*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX3 + JYAvg*FDY3)/2.d0 &

    !...     Spatially varying Tau0 terms
        +Tau0SpaVar

    !...     Put these partial products into further elemental storage for a vector computer
    !...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
        Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
        Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
        Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

    !...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
        GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
        GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
        GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif

        CONTINUE                  !End of elemental loop
    1037 END DO

!...  Put load vector elemental values into nodal storage for a vector computer
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

!...  Save the elevation at the past time step into Eta1 and zero out Eta2
! md  Need to save z(s-1) and etas(s-1) for the corrector loop

    DO I=1,NP
        IF(CPRECOR) THEN
            ETAS0(I)=ETAS(I)
            Eta0(I)=Eta1(I)
        END IF
        Eta1(I)=Eta2(I)
        Eta2(I)=0.0d0
    END DO

!...  At elevation boundary condition nodes, determine the elevation at
!...  the s+1 time step
!...
!...  For periodic elevation boundary conditions

    DO J=1,NBFR
        IF(PER(J) == 0.) THEN
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

!...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

    if (subdomainOn) then                     ! NCSU Subdomain
        if(enforceBN == 1) call enforceEcb()  ! NCSU Subdomain
    else                                      ! NCSU Subdomain

        IF((NBFR == 0) .AND. (NOPE > 0)) THEN
            IF(TimeLoc > ETIME2) THEN
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
                Eta2(NBDI)=RampElev &
                *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
            END DO
        ENDIF
    endif                                   ! NCSU Subdomain

! jgf46.02 Added the ability to include geoid offset on the boundary.

! aaltuntas51.48: Deactivated geoidOffset for subdomain boundary
! conditions.
    if (subdomainOn .AND. enforceBN == 1) then  ! NCSU Subdomain
        continue                                ! NCSU Subdomain
    else
        IF (LoadGeoidOffset) THEN
            DO I=1,NETA
                ETA2(NBD(I))=ETA2(NBD(I))+GeoidOffset(NBD(I))
            END DO
        ENDIF
    endif                                      ! NCSU Subdomain


!     jgf48.04 Added an inverted barometer boundary condition so that
!     low pressure systems can cross the boundary without creating an
!     elevation anomaly.
!     jgf52.30.04: Included a parameter from fort.15 to turn this
!     on and off according to analyst preference.
    if (invertedBarometerOnElevationBoundary.eqv. .TRUE. ) then
        DO I=1,NETA
            ETA2(NBD(I))=ETA2(NBD(I)) &
            + RampMete*(101300.d0/(RHOWAT0*G)) - PR2(NBD(I))
        END DO
    endif

!   kmd48.33bc add information for the levels of no motion boundary conditions
!              these are considered the steric adjustments.
    IF ((ABS(RES_BC_FLAG) >= 1) .AND. (CBaroclinic) .AND. (NOPE > 0))THEN
        DO I=1,NETA
            NBDI=NBD(I)
            ETA2(NBDI) = ETA2(NBDI) + LNM_BC(I)
        END DO
    END IF

!...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
!...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

!...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
!...  assume that QN is positive into the domain.  QFORCEJ has a -1
!...  built in and the terms are not explicitly negated. Boundary
!...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
!...  flux that  is positive out of the domain.  Therefore, to match
!...  the formulation these terms must be explicitly multiplied by -1.

!...Note 3, Eta1 is the latest computed elevation (it was updated above).

    IF((NFLUXF == 1) .OR. (NFLUXB == 1) .OR. (NFLUXIB == 1) &
     .OR. (NFLUXGBC == 1) .OR. (NFLUXRBC == 1)) THEN
        NBDJ=NBV(1)
        IF(LBCODEI(1) <= 29) QFORCEJ=(QN2(1)-QN0(1))/DT2 + &
        Tau0VAR(NBDJ)*QN1(1)

        IF(LBCODEI(1) == 30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
        ENDIF

        IF(LBCODEI(1) == 32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT &
            -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT &
            +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
        ENDIF

        IF((LBCODEI(1) == 40) .OR. (LBCODEI(1) == 41)) QFORCEJ= &
        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

    !     jgf46.21 Added IBTYPE=52.
        IF(LBCODEI(1) == 52) THEN
            QFORCEJ=(QN2(1)-QN0(1))/DT2 + Tau0VAR(NBDJ)*QN1(1)
            IF (IT > FluxSettlingIT) THEN
                HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
                Celerity=SQRT(G*HH1)
                QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT &
                + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(1)))
            ENDIF
        ENDIF
    
        DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J) <= 29) QFORCEJ=(QN2(J)-QN0(J))/DT2+ &
            Tau0VAR(NBDJ)*QN1(J)
            IF(LBCODEI(J) == 30) THEN
                H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
                CELERITY=SQRT(G*H1)
                QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
            ENDIF

            IF(LBCODEI(J) == 32) THEN
                H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
                CELERITY=SQRT(G*H1)
                QFORCEJ=(QN1(J)-QN0(J))/DT &
                -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT &
                +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
            ENDIF

            IF((LBCODEI(J) == 40) .OR. (LBCODEI(J) == 41)) QFORCEJ= &
            -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0


            IF(LBCODEI(J) <= 29) QFORCEJ=(QN2(J)-QN0(J))/DT2+ &
            Tau0VAR(NBDJ)*QN1(J)
        
        !     jgf46.21 Added IBTYPE=52
            IF(LBCODEI(J) == 52) THEN
                QFORCEJ=(QN2(J)-QN0(J))/DT2 + Tau0VAR(NBDJ)*QN1(J)
                IF (IT > FluxSettlingIT) THEN
                    HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
                    Celerity=SQRT(G*HH1)
                    QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT &
                    + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(J)))
                ENDIF
            ENDIF

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI) &
            + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ) &
            + BndLenO6NC*(2.d0*QForceJ+QForceI)
        ENDDO
    ENDIF

!...
!...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
!...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
!...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
!...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
!...  EQUATIONS
!...
    IF ( ILump == 0 ) THEN ! default, fully consistent GWCE LHS
        DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
            DO J=2,NNEIGH(NBDI)
                GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J)) &
                -ETAS(NBDI)*OBCCOEF(I,J-1)
            END DO
        END DO
    ELSE                   ! ILump == 1, lumped GWCE
        DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*COEFD(NBDI)
        END DO
    ENDIF
!...
!...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
!...

!...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
!...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
    IF (ILump == 0) THEN ! default, fully consistent LHS
        CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
    ELSE  ! lumped LHS
        CALL UPDATER(GWCE_LV,COEFD,DUMY1,2)
    ENDIF
#endif

    if (subdomainOn .AND. enforceBN == 2) call enforceEob() ! NCSU Subdomain
    if (subdomainOn .AND. enforceBN == 2) call enforceEib() ! NCSU Subdomain
    if (subdomainOn .AND. enforceBN == 2) call enforceGWCELVob() ! NCSU Subdomain

    IF (ILump == 0) THEN ! default, fully consistent LHS
    !...  JCG ITERATIVE MATRIX SOLVER
        IPARM(1)=ITMAX
        CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS, &
        IWKSP,NW,WKSP,IPARM,RPARM,IER)

        NUMITR=IPARM(1)
        DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
        END DO
    ELSE ! lumped LHS
        DO I = 1, NP
            IF (COEFD(I) == 0.0d0) THEN
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

#ifdef CMPI
    CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN

!**********************************************************************
    END SUBROUTINE GWCE_NEW
!**********************************************************************


!*******************************************************************************
!                                                                              *
!    Subroutine to compute the velocity and flux/unit width using a 2DDI non-  *
!    conservative momentum equation formluation as contained in original       *
!    versions of ADCIRC (e.g., pre v45.XX).  This formulation is correct if    *
!    the element size is spatially constant but does not account properly for  *
!    spatially variable elemental areas.                                       *
!    Despite its incorrect formulation, we have maintained this subroutine for *
!    posterity and comparison with earlier runs.                               *
!                                                                              *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

    SUBROUTINE Mom_Eqs_Original()


    USE GLOBAL
    USE WIND
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, MJU, &
    NEITAB, NEITABELE, NNEIGH, SFAC
    USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, &
    CSII, SIII, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
    USE NodalAttributes, ONLY: EVM,LoadAdvectionState
#ifdef CSWAN
#ifdef CSWANFRIC
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
    USE Couple2Swan, ONLY: TKXX, &
    TKXY, &
    TKYY
#endif
#endif

    IMPLICIT NONE

    INTEGER :: IE, I, J, N                           !local loop counters
    INTEGER ::  NM1, NM2, NM3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    INTEGER :: NBDI
    INTEGER :: NNFirst

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
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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
    call setMessageSource("mom_eqs_original")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif
!...
!...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
!...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
!...        THE TIME STEPPING LOOP.

!.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

    DO IE=1,NE

    !...  SET NODAL VALUES FOR EACH ELEMENT


    ! Corbitt 120322: Localized Advection
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

    !...  Accumulate nodal values for terms associated with the barotropic pressure gradient

        VCOEF3N1=ETA1(NM1)+ETA2(NM1)
        VCOEF3N2=ETA1(NM2)+ETA2(NM2)
        VCOEF3N3=ETA1(NM3)+ETA2(NM3)

    !......If using atm pressure

        IF(NWS /= 0) THEN
            VCOEF3N1=VCOEF3N1+PR1(NM1)+PR2(NM1)
            VCOEF3N2=VCOEF3N2+PR1(NM2)+PR2(NM2)
            VCOEF3N3=VCOEF3N3+PR1(NM3)+PR2(NM3)
        ENDIF

    !......If using tidal potential terms

        if (CTIP) then
            VCOEF3N1=VCOEF3N1-TIP1(NM1)-TIP2(NM1)
            VCOEF3N2=VCOEF3N2-TIP1(NM2)-TIP2(NM2)
            VCOEF3N3=VCOEF3N3-TIP1(NM3)-TIP2(NM3)
        endif

        VCOEF3N1=VCOEF3N1*GDTO2
        VCOEF3N2=VCOEF3N2*GDTO2
        VCOEF3N3=VCOEF3N3*GDTO2

    !...  COMPUTE ELEMENT AVERAGES QUANTITIES
    ! Corbitt 120322: Local Advection
        U1AvgDT=IFNLCT*DT/3.D0*(U1N1+U1N2+U1N3)
        V1AvgDT=IFNLCT*DT/3.D0*(V1N1+V1N2+V1N3)
        U1AvgDTDDX1=U1AvgDT*DDX1
        U1AvgDTDDX2=U1AvgDT*DDX2
        U1AvgDTDDX3=U1AvgDT*DDX3
        V1AvgDTDDY1=V1AvgDT*DDY1
        V1AvgDTDDY2=V1AvgDT*DDY2
        V1AvgDTDDY3=V1AvgDT*DDY3
        EVMAvgDT=((EVMN1+EVMN2+EVMN3)/3.d0)*DT

    !...  ASSEMBLE PARTIAL PRODUCTS

        VCOEF3X=VCOEF3N1*DDX1+VCOEF3N2*DDX2+VCOEF3N3*DDX3
        VCOEF3Y=VCOEF3N1*DDY1+VCOEF3N2*DDY2+VCOEF3N3*DDY3
        ADVECX=(U1AvgDTDDX1+V1AvgDTDDY1)*U1N1 &
        +(U1AvgDTDDX2+V1AvgDTDDY2)*U1N2 &
        +(U1AvgDTDDX3+V1AvgDTDDY3)*U1N3
        ADVECY=(U1AvgDTDDX1+V1AvgDTDDY1)*V1N1 &
        +(U1AvgDTDDX2+V1AvgDTDDY2)*V1N2 &
        +(U1AvgDTDDX3+V1AvgDTDDY3)*V1N3
    
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A
    !.... VECTOR FOR NODE NM1

        TEMP_LV_A1=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3X &
    !...  LATERAL VISCOUS TERMS
        -EVMAvgDT*(DXXYY11*U1N1+DXXYY12*U1N2+DXXYY13*U1N3) &
    !...  ADVECTIVE TERMS
        -ADVECX &
    !...  COMMON DIVISION OPERATION
        )/FIIN
    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO TEMP_LV_A
    !.... VECTOR FOR NODE NM2
    !...
        TEMP_LV_A2=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3X &
    !...  LATERAL VISCOUS TERMS]
        -EVMAvgDT*(DXXYY21*U1N1+DXXYY22*U1N2+DXXYY23*U1N3) &
    !...  ADVECTIVE TERMS
        -ADVECX &
    !...  COMMON DIVISION OPERATION
        )/FIIN
    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO TEMP_LV_A
    !.... VECTOR FOR NODE NM3
    !...
        TEMP_LV_A3=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3X &
    !...  LATERAL VISCOUS TERMS
        -EVMAvgDT*(DXXYY31*U1N1+DXXYY32*U1N2+DXXYY33*U1N3) &
    !...  ADVECTIVE TERMS
        -ADVECX &
    !...  COMMON DIVISION OPERATION
        )/FIIN
    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
    !.... VECTOR FOR NODE NM1
    !...
        TEMP_LV_B1=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3Y &
    !...  LATERAL VISCOUS TERMS
        -EVMAvgDT*(DXXYY11*V1N1+DXXYY12*V1N2+DXXYY13*V1N3) &
    !...  ADVECTIVE TERMS
        -ADVECY &
    !...  COMMON DIVISION OPERATION
        )/FIIN
    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
    !.... VECTOR FOR NODE NM2
    !...
        TEMP_LV_B2=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3Y &
    !...  LATERAL VISCOUS TERMS
        -EVMAvgDT*(DXXYY21*V1N1+DXXYY22*V1N2+DXXYY23*V1N3) &
    !...  ADVECTIVE TERMS
        -ADVECY &
    !...  COMMON DIVISION OPERATION
        )/FIIN
    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO TEMP_LV_B
    !.... VECTOR FOR NODE NM3
    !...
        TEMP_LV_B3=NCELE*( &
    !...  SURFACE GRADIENT, ATMOSPHERIC PRESSURE AND TIDAL POTENTIAL TERMS
        -VCOEF3Y &
    !...  LATERAL VISCOUS TERMS
        -EVMAvgDT*(DXXYY31*V1N1+DXXYY32*V1N2+DXXYY33*V1N3) &
    !...  ADVECTIVE TERMS
        -ADVECY &
    !...  COMMON DIVISION OPERATION
        )/FIIN

    !     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
        TEMP_LV_A(IE,1)=TEMP_LV_A1
        TEMP_LV_A(IE,2)=TEMP_LV_A2
        TEMP_LV_A(IE,3)=TEMP_LV_A3
        TEMP_LV_B(IE,1)=TEMP_LV_B1
        TEMP_LV_B(IE,2)=TEMP_LV_B2
        TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

    !     LINES TO RUN ON A SCALAR COMPUTER
    !     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
    !           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
        MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
        MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
        MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
        MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
        MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
        MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

    ENDDO                        !end of elemental assembly loop


!     LINES TO RUN ON A VECTOR COMPUTER
!     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
!           AND AUV
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


!...  UPDATE MOMENTUM EQUATION LHS COEFFICIENTS AND LOAD VECTORS AT EACH
!...  NODE BY DIVIDING BY NUMBER OF ELEMENTS THE NODE IS ASSOCIATED WITH
!...  AND ADDING IN LUMPED TERMS AND BOTTOM FRICTION AND TAKING ACCOUNT
!...  OF THE BOUNDARY CONDITION

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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
        ENDIF
        VCoef1=DTO2*TK(I)
        VCoef2=DTO2*CORIF(I)
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCOEFXX)*UU1(I)-VCOEFXY*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCOEFYY)*VV1(I)-VCOEFYX*UU1(I)-VIDBCPDX)
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I) &
        +VCoef2*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I) &
        -VCoef2*UU1(I)-VIDBCPDY)
#endif
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I) &
        +VCoef2*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I) &
        -VCoef2*UU1(I)-VIDBCPDY)
#endif

        AUV11(I)=1.D0+VCoef1*NCI
        AUV12(I)=-VCoef2*NCI

#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVXX(I) = 1.D0 + VCOEFXX*NCI
        AUVYY(I) = 1.D0 + VCOEFYY*NCI
        AUVXY(I) = VCOEFXY*NCI
        AUVYX(I) = VCOEFXY*NCI
#endif
#endif
    END DO

!...  Modify the momentum equations to impose velocity boundary
!...  conditions In each case the equations are manipulated to
!...  maintain the LHS matrix structure of AUV11=AUV22;
!...  AUV12=-AUV21)

    DO J=1,NVELME
        I=ME2GW(J)
        NBDI=NBV(I)
        H2=DP(NBDI)+IFNLFA*ETA2(NBDI)
        NCI=NODECODE(NBDI)

    !      Specified essential normal flow and free tangential slip

        IF((LBCODEI(I) >= 0) .AND. (LBCODEI(I) <= 9)) THEN
            VELNORM=-QN2(I)/H2
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VELNORM*AUVXY(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUVXX(NBDI)*NCI !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVYX(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VELNORM*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUV11(NBDI)*NCI !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VELNORM*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*AUV11(NBDI)*NCI !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
        ENDIF

    !     Specified essential normal flow and no tangential slip

        IF((LBCODEI(I) >= 10) .AND. (LBCODEI(I) <= 19)) THEN
            VELNORM=-QN2(I)/H2
            VELTAN=0.D0
            MOM_LV_X(NBDI)=VELTAN*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VELNORM*NCI !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI) = SIII(I)
            AUVXY(NBDI) = - CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#endif
#endif
        ENDIF

    !     Zero normal velocity gradient using a Galerkin approximation to
    !     the normal derivatives. Note: this is fully explicit and therefore
    !     the velocity at the boundary is computed entirely from surrounding
    !     velocities at the previous time step.

        IF(LBCODEI(I) == 41) THEN
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
                NEle=NeiTabEle(NBDI,N-2) !element# defined by nodes NM1,NM2,NM3
                NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
                IF((NEle /= 0) .AND. (NCEle /= 0)) THEN !if element is active, compute contribution
                    FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                    FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                    FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                    FDY1 = X(NM3)-X(NM2) !a1
                    FDY2 = X(NM1)-X(NM3) !a2
                    FDY3 = X(NM2)-X(NM1) !a3
                    ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                    ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                    ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                ENDIF
            END DO
            NM3=NM2          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle /= 0) .AND. (NCEle /= 0)) THEN
                FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                FDY1 = X(NM3)-X(NM2) !a1
                FDY2 = X(NM1)-X(NM3) !a2
                FDY3 = X(NM2)-X(NM1) !a3
                ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
            ENDIF
            IF(NCI == 0) THEN
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
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI)=1.D0
            AUVXY(NBDI)=0.D0
#endif
#endif
        ENDIF

    ENDDO

!...
!...  SOLVE FOR VELOCITY AT NEW LEVEL  (K+1)
!...

!.....Note: This includes the comparison between MJU and NODELE to
!.....determine if the node is an interface node.  If MJU < NODELE the
!.....velocity can be zeroed out to obtain an essential zero velocity at
!.....interface nodes.

    DO I=1,NP
        AUV22=AUV11(I)
        AUV21=-AUV12(I)
        DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
        UU2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
        VV2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVYY(I) = AUVXX(I)
        AUVYX(I) = - AUVXY(I)
        DDU = AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
        UU2(I) = (MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
        VV2(I) = (MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
    ! jw
    !        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
    !           UU2(I)=0.D0    !no slip and normal flux
    !           VV2(I)=0.D0    !on wet/dry interface nodes
    !           ENDIF
    END DO

!...
!...  Impose a zero normal velocity gradient based on interpolating the
!...  velocity at a fictitious point in the interior of the domain,
!...  normal to a specified boundary node and setting the boundary
!...  velocity equal to the interpolated value at the fictitious point.
!...  Provided the fictitious point does not lie in an element that
!...  contains a boundary point, this is an entirely implicit
!...  calculation.
!...
    IF(NFLUXGBC == 1) THEN
        DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I) == 40) THEN
                NM1=NM(NEleZNG(I),1)
                NM2=NM(NEleZNG(I),2)
                NM3=NM(NEleZNG(I),3)
                NC1=NODECODE(NM1)
                NC2=NODECODE(NM2)
                NC3=NODECODE(NM3)
                NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
                UU2(NBDI)=NCEle*(UU2(NM1)*ZNGIF1(I)+UU2(NM2)*ZNGIF2(I) &
                +UU2(NM3)*ZNGIF3(I))
                VV2(NBDI)=NCEle*(VV2(NM1)*ZNGIF1(I)+VV2(NM2)*ZNGIF2(I) &
                +VV2(NM3)*ZNGIF3(I))
            ENDIF
        ENDDO
    ENDIF

!...  Compute fluxes

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
!**********************************************************************
    END SUBROUTINE MOM_EQS_ORIGINAL
!**********************************************************************

!*******************************************************************************
!                                                                              *
!   Subroutine to compute the velocity and from that the flux/unit width using *
!   a 2DDI non conservative momentum equation.                                 *
!                                                                              *
!   Options are provided for either the correct area integration or the        *
!   original incorrect area integration.                                       *
!                                                                              *
!   Options are provided to use either flux or velocity based lateral          *
!   viscosity.                                                                 *
!                                                                              *
!   For a uniform grid and velocity based lateral viscosity, this subroutine   *
!   should give the same results as the original nonconservative formulation.  *
!                                                                              *
!   This subroutine follows the naming convention and formulation in the new   *
!   ADCIRC theory report.                                                      *
!                                                                              *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

    SUBROUTINE Mom_Eqs_New_NC()

    USE GLOBAL
    USE WIND
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, Areas, TotalArea, &
    NNeigh, NeiTab, NeiTabEle, MJU, SFAC
    USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, &
    CSII, SIII, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
    USE NodalAttributes, ONLY: EVM,LoadAdvectionState
#ifdef CSWAN
#ifdef CSWANFRIC
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
    USE Couple2Swan, ONLY: TKXX, &
    TKXY, &
    TKYY
#endif
#endif
    USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceUVOB, &
    enforceUVCB

    IMPLICIT NONE

    INTEGER :: IE, I, J, N                           !local loop counters
    INTEGER ::  NM1, NM2, NM3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    INTEGER :: NBDI
    INTEGER :: NNFirst

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
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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
!...
!...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
!...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
!...        THE TIME STEPPING LOOP.

!.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

    DO IE=1,NE

    !...  SET NODAL VALUES FOR EACH ELEMENT


    ! Corbitt 120322: Localized Advection
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

    !...  Compute element averaged quantities

        U1Avg =(U1N1+U1N2+U1N3)/3.d0
        V1Avg =(V1N1+V1N2+V1N3)/3.d0

        DU1DXA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
        DU1DYA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
        DV1DXA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
        DV1DYA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0

        EVMEle=(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
        IF(CSmag_Eh) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
            EVMSmag=EVMEle* &
            sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA) &
            +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
            EVMEle=EVMSmag
        ENDIF

    !...  Compute terms associated with the barotropic pressure

        BTP1N1=ETA1(NM1)
        BTP2N1=ETA2(NM1)
        BTP1N2=ETA1(NM2)
        BTP2N2=ETA2(NM2)
        BTP1N3=ETA1(NM3)
        BTP2N3=ETA2(NM3)

    !.......If using atm pressure add it into the barotropic pressure

        IF(NWS /= 0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
        ENDIF

    !.......If using tidal potential terms, add these into the barotropic pressure

        IF (CTIP) THEN
            BTP1N1=BTP1N1-TiP1(NM1)
            BTP2N1=BTP2N1-TiP2(NM1)
            BTP1N2=BTP1N2-TiP1(NM2)
            BTP2N2=BTP2N2-TiP2(NM2)
            BTP1N3=BTP1N3-TiP1(NM3)
            BTP2N3=BTP2N3-TiP2(NM3)
        ENDIF

    !...  Compute the barotropic pressure gradient x area for the element

        DBTPDXA=((BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3) &
        +(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
        DBTPDYA=((BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3) &
        +(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

    !...  Compute the advective term gradients x area for the element

        U1AvgDU1DXA=U1Avg*DU1DXA
        V1AvgDU1DYA=V1Avg*DU1DYA
        U1AvgDV1DXA=U1Avg*DV1DXA
        V1AvgDV1DYA=V1Avg*DV1DYA

    !...  Compute the lateral viscous terms for the element (flux formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric flux formulation)

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

    !...  Compute the lateral viscous terms for the element (velocity formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

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

    
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM1

        TEMP_LV_A1=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1)/H1N1 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM2
    !...
        TEMP_LV_A2=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2)/H1N2 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM3
    !...
        TEMP_LV_A3=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDU1DXA+V1AvgDU1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3)/H1N3 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM1

        TEMP_LV_B1=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1)/H1N1 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM2

        TEMP_LV_B2=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2)/H1N2 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM3

        TEMP_LV_B3=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(U1AvgDV1DXA+V1AvgDV1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3)/H1N3 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !     Original (incorrect) area integration - for historical comparison

        IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
        ENDIF

    !     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
        TEMP_LV_A(IE,1)=TEMP_LV_A1
        TEMP_LV_A(IE,2)=TEMP_LV_A2
        TEMP_LV_A(IE,3)=TEMP_LV_A3
        TEMP_LV_B(IE,1)=TEMP_LV_B1
        TEMP_LV_B(IE,2)=TEMP_LV_B2
        TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

    !     LINES TO RUN ON A SCALAR COMPUTER
    !     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
    !           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
        MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
        MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
        MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
        MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
        MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
        MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

    ENDDO


!     LINES TO RUN ON A VECTOR COMPUTER
!     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
!           AND AUV
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

!...  Update the momentum equation LHS coefficients and load vectors at each
!...  node by dividing by the area of all active elements attached to the node
!...  and adding in the lumped terms, bottom friction and boundary conditions

    WSX=0.D0
    WSY=0.D0
    VIDBCPDX=0.D0
    VIDBCPDY=0.D0
    DO I=1,NP
        NCI=NODECODE(I)
        IF(TotalArea(I) /= 0.d0) THEN
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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
        ENDIF
        VCoef1=DTO2*TK(I)                     !TK = Kslip/H
        VCoef2=DTO2*CORIF(I)
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCOEFXX)*UU1(I)-VCOEFXY*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCOEFYY)*VV1(I)-VCOEFYX*UU1(I)-VIDBCPDX)
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I) &
        +VCoef2*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I) &
        -VCoef2*UU1(I)-VIDBCPDY)
#endif
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I) &
        +VCoef2*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I) &
        -VCoef2*UU1(I)-VIDBCPDY)
#endif

        AUV11(I)=1.D0+VCoef1*NCI
        AUV12(I)=-VCoef2*NCI

#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVXX(I) = 1.D0 + VCOEFXX*NCI
        AUVYY(I) = 1.D0 + VCOEFYY*NCI
        AUVXY(I) = VCOEFXY*NCI
        AUVYX(I) = VCOEFXY*NCI
#endif
#endif

    END DO

!...  Modify the momentum equations to impose velocity boundary
!...  conditions In each case the equations are manipulated to
!...  maintain the LHS matrix structure of AUV11=AUV22;
!...  AUV12=-AUV21)

    DO J=1,NVELME
        I=ME2GW(J)
        NBDI=NBV(I)
        H2=DP(NBDI)+IFNLFA*ETA2(NBDI)
        NCI=NODECODE(NBDI)

    !      Specified essential normal flow and free tangential slip

        IF((LBCODEI(I) >= 0) .AND. (LBCODEI(I) <= 9)) THEN
            VelNorm=-QN2(I)/H2
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VelNorm*AUVXY(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUVXX(NBDI)*NCI   !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVXY(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
        ENDIF

    !     Specified essential normal flow and no tangential slip

        IF((LBCODEI(I) >= 10) .AND. (LBCODEI(I) <= 19)) THEN
            VelNorm=-QN2(I)/H2
            VelTan=0.D0
            MOM_LV_X(NBDI)=VelTan*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*NCI !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI) = SIII(I)
            AUVXY(NBDI) = - CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#endif
#endif
        ENDIF

    !     Zero normal velocity gradient using a Galerkin approximation to
    !     the normal derivatives. Note: this is fully explicit and therefore
    !     the velocity at the boundary is computed entirely from surrounding
    !     velocities at the previous time step.

        IF(LBCODEI(I) == 41) THEN
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
                NEle=NeiTabEle(NBDI,N-2) !element# defined by nodes NM1,NM2,NM3
                NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
                IF((NEle /= 0) .AND. (NCEle /= 0)) THEN !if element is active, compute contribution
                    FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                    FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                    FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                    FDY1 = X(NM3)-X(NM2) !a1
                    FDY2 = X(NM1)-X(NM3) !a2
                    FDY3 = X(NM2)-X(NM1) !a3
                    ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                    ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                    ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                ENDIF
            END DO
            NM3=NM2          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle /= 0) .AND. (NCEle /= 0)) THEN
                FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                FDY1 = X(NM3)-X(NM2) !a1
                FDY2 = X(NM1)-X(NM3) !a2
                FDY3 = X(NM2)-X(NM1) !a3
                ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
            ENDIF
            IF(NCI == 0) THEN
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
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI)=1.D0
            AUVXY(NBDI)=0.D0
#endif
#endif
        ENDIF

    ENDDO
!...
!...  SOLVE FOR VELOCITY AT NEW LEVEL  (s+1)
!...

!.....Note: This includes the comparison between MJU and NODELE to
!.....determine if the node is an interface node.  If MJU < NODELE the
!.....velocity can be zeroed out to obtain an essential zero velocity at
!.....interface nodes.

    DO I=1,NP
        AUV22=AUV11(I)
        AUV21=-AUV12(I)
        DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
        UU2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
        VV2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVYY(I)= AUVXX(I)
        AUVYX(I)= - AUVXY(I)
        DDU = AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
        UU2(I)=(MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
        VV2(I)=(MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
    !        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
    !           UBAR2(I)=0.D0    !no slip and normal flux
    !           VBAR2(I)=0.D0    !on wet/dry interface nodes
    !           ENDIF
    END DO

    if(subdomainOn .AND. enforceBN == 1) call enforceUVcb() ! NCSU Subdomain
    if(subdomainOn .AND. enforceBN == 2) call enforceUVob() ! NCSU Subdomain
!...
!...  Impose a zero normal velocity gradient based on interpolating the
!...  velocity at a fictitious point in the interior of the domain,
!...  normal to a specified boundary node and setting the boundary
!...  velocity equal to the interpolated value at the fictitious point.
!...  Provided the fictitious point does not lie in an element that
!...  contains a boundary point, this is an entirely implicit
!...  calculation.
!...
    IF(NFLUXGBC == 1) THEN
        DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I) == 40) THEN
                NM1=NM(NEleZNG(I),1)
                NM2=NM(NEleZNG(I),2)
                NM3=NM(NEleZNG(I),3)
                NC1=NODECODE(NM1)
                NC2=NODECODE(NM2)
                NC3=NODECODE(NM3)
                NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
                UU2(NBDI)=NCEle*(UU2(NM1)*ZNGIF1(I)+UU2(NM2)*ZNGIF2(I) &
                +UU2(NM3)*ZNGIF3(I))
                VV2(NBDI)=NCEle*(VV2(NM1)*ZNGIF1(I)+VV2(NM2)*ZNGIF2(I) &
                +VV2(NM3)*ZNGIF3(I))
            ENDIF
        ENDDO
    ENDIF

!...  Compute fluxes

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
!*******************************************************************************
    END SUBROUTINE MOM_EQS_NEW_NC
!*******************************************************************************



!*******************************************************************************
!                                                                              *
!   Subroutine to compute the flux/unit width and from that the velocity using *
!   2DDI conservative momentum equation formulations version 1 or 2.           *
!                                                                              *
!   Options are provided for either the correct area integration or the        *
!   original incorrect area integration.                                       *
!                                                                              *
!   Options are provided to use either flux or velocity based lateral          *
!   viscosity.                                                                 *
!                                                                              *
!   This subroutine follows the naming convention and formulation in the new   *
!   ADCIRC theory report.                                                      *
!                                                                              *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

    SUBROUTINE Mom_Eqs_New_Conserv()

    USE GLOBAL
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, TotalArea, MJU, &
    NEITAB, NEITABELE, NNEIGH, SFAC
    USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, &
    CSII, SIII, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
    USE WIND
    USE NodalAttributes, ONLY: EVM,LoadAdvectionState
#ifdef CSWAN
#ifdef CSWANFRIC
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
    USE Couple2Swan, ONLY: TKXX, &
    TKXY, &
    TKYY
#endif
#endif
    IMPLICIT NONE

    INTEGER :: IE, I, J, N                           !local loop counters
    INTEGER ::  NM1, NM2, NM3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    INTEGER :: NBDI
    INTEGER :: NNFirst

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
! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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

!...
!...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
!...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
!...        THE TIME STEPPING LOOP.

!.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

    DO IE=1,NE

    !...  SET NODAL VALUES FOR EACH ELEMENT


    ! Corbitt 120322: Localized Advection
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

    !...  Compute element averaged quantities

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

            EVMSmag=0.5d0*EVMEle* &
            sqrt((DU1DX2A-DV1DY2A)*(DU1DX2A-DV1DY2A) &
            +(DU1DY2A+DV1DX2A)*(DU1DY2A+DV1DX2A))
            EVMEle=EVMSmag
        ENDIF

    !...  Compute terms associated with the barotropic pressure

        BTP1N1=ETA1(NM1)
        BTP2N1=ETA2(NM1)
        BTP1N2=ETA1(NM2)
        BTP2N2=ETA2(NM2)
        BTP1N3=ETA1(NM3)
        BTP2N3=ETA2(NM3)

    !.......If using atm pressure add it into the barotropic pressure

        IF(NWS /= 0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
        ENDIF

    !.......If using tidal potential terms, add these into the barotropic pressure

        IF (CTIP) THEN
            BTP1N1=BTP1N1-TIP1(NM1)
            BTP2N1=BTP2N1-TIP2(NM1)
            BTP1N2=BTP1N2-TIP1(NM2)
            BTP2N2=BTP2N2-TIP2(NM2)
            BTP1N3=BTP1N3-TIP1(NM3)
            BTP2N3=BTP2N3-TIP2(NM3)
        ENDIF

    !...  Compute the element avg depth x barotropic pressure gradient x area for the element

        DBTPDXA=(H1Avg*(BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3) &
        +H2Avg*(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
        DBTPDYA=(H1Avg*(BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3) &
        +H2Avg*(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

    !...  Compute the advective term gradients for the element version 1

        IF (CME_New_C1) THEN
            DU1QX1DXA=(U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2+U1N3*QX1N3*FDX3) &
            /2.d0
            DV1QX1DYA=(V1N1*QX1N1*FDY1+V1N2*QX1N2*FDY2+V1N3*QX1N3*FDY3) &
            /2.d0
            DU1QY1DXA=(U1N1*QY1N1*FDX1+U1N2*QY1N2*FDX2+U1N3*QY1N3*FDX3) &
            /2.d0
            DV1QY1DYA=(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2+V1N3*QY1N3*FDY3) &
            /2.d0
        ENDIF

    !...  Compute the advective term gradients for the element version 2

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

    !...  Compute the lateral viscous terms for the element (flux formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric flux formulation)

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

    !...  Compute the lateral viscous terms for the element (velocity formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

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

    
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM1

        TEMP_LV_A1=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QX1DXA+DV1QX1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM2
    !...
        TEMP_LV_A2=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QX1DXA+DV1QX1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM3
    !...
        TEMP_LV_A3=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QX1DXA+DV1QX1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM1

        TEMP_LV_B1=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QY1DXA+DV1QY1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM2

        TEMP_LV_B2=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QY1DXA+DV1QY1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM3

        TEMP_LV_B3=NCELE*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*(DU1QY1DXA+DV1QY1DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3) &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !     Original (incorrect) area integration - for historical comparison

        IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
        ENDIF

    !     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
        TEMP_LV_A(IE,1)=TEMP_LV_A1
        TEMP_LV_A(IE,2)=TEMP_LV_A2
        TEMP_LV_A(IE,3)=TEMP_LV_A3
        TEMP_LV_B(IE,1)=TEMP_LV_B1
        TEMP_LV_B(IE,2)=TEMP_LV_B2
        TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

    !     LINES TO RUN ON A SCALAR COMPUTER
    !     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
    !           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
        MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
        MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
        MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
        MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
        MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
        MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

    ENDDO



!     LINES TO RUN ON A VECTOR COMPUTER
!     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
!           AND AUV
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


!...  Update the momentum equation LHS coefficients and load vectors at each
!...  node by dividing by the area of all active elements attached to the node
!...  and adding in the lumped terms, bottom friction and boundary conditions

    WSX=0.D0
    WSY=0.D0
    VIDBCPDX=0.D0
    VIDBCPDY=0.D0
    DO I=1,NP
        NCI=NODECODE(I)
        IF(TotalArea(I) /= 0.d0) THEN
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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)+WSX2(I))
            WSY=DTO2*IFWIND*(WSY1(I)+WSY2(I))
        ENDIF
        VCoef1=DTO2*TK(I)
        VCoef2=DTO2*CORIF(I)
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
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
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCOEFXX)*QX1(I)-VCOEFXY*QY1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCOEFYY)*QY1(I)-VCOEFYX*QX1(I)-VIDBCPDX)
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*QX1(I) &
        +VCoef2*QY1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*QY1(I) &
        -VCoef2*QX1(I)-VIDBCPDY)
#endif
#else
        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*QX1(I) &
        +VCoef2*QY1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*QY1(I) &
        -VCoef2*QX1(I)-VIDBCPDY)
#endif

        AUV11(I)=1.D0+VCoef1*NCI
        AUV12(I)=-VCoef2*NCI

#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVXX(I) = 1.D0 + VCOEFXX*NCI
        AUVYY(I) = 1.D0 + VCOEFYY*NCI
        AUVXY(I) = VCOEFXY*NCI
        AUVYX(I) = VCOEFXY*NCI
#endif
#endif

    END DO

!...  Modify the momentum equations to impose velocity boundary
!...  conditions In each case the equations are manipulated to
!...  maintain the LHS matrix structure of AUV11=AUV22;
!...  AUV12=-AUV21)

    DO J=1,NVELME
        I=ME2GW(J)
        NBDI=NBV(I)
        NCI=NODECODE(NBDI)

    !      Specified essential normal flow and free tangential slip

        IF((LBCODEI(I) >= 0) .AND. (LBCODEI(I) <= 9)) THEN
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            +QN2(I)*AUVXY(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUVXX(NBDI)*NCI           !Normal Eqn RHS
            AUVXX(NBDI) = AUVXX(NBDI)*SIII(I) - AUVYX(NBDI)*CSII(I)
            AUVXY(NBDI) = AUVXY(NBDI)*SIII(I) - AUVYY(NBDI)*CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            +QN2(I)*AUV12(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUV11(NBDI)*NCI           !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
#else
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            +QN2(I)*AUV12(NBDI))*NCI          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*AUV11(NBDI)*NCI           !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
#endif
        ENDIF

    !     Specified essential normal flow and no tangential slip

        IF((LBCODEI(I) >= 10) .AND. (LBCODEI(I) <= 19)) THEN
            QTAN=0.D0
            MOM_LV_X(NBDI)=QTan*NCI                          !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=-QN2(I)*NCI                       !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
#ifdef CSWAN
#ifdef CSWANFRIC
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI) = SIII(I)
            AUVXY(NBDI) = - CSII(I)
            AUVYX(NBDI) = CSII(I)
            AUVYY(NBDI) = SIII(I)
#endif
#endif
        ENDIF

    !     Zero normal flux gradient using a Galerkin approximation to
    !     the normal derivatives. Note: this is fully explicit and therefore
    !     the flux at the boundary is computed entirely from surrounding
    !     fluxes at the previous time step.

        IF(LBCODEI(I) == 41) THEN
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
                NEle=NeiTabEle(NBDI,N-2)                      !element# defined by nodes NM1,NM2,NM3
                NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
                IF((NEle /= 0) .AND. (NCEle /= 0)) THEN         !if element is active, compute contribution
                    FDX1 = (Y(NM2)-Y(NM3))*SFacAvg             !b1
                    FDX2 = (Y(NM3)-Y(NM1))*SFacAvg             !b2
                    FDX3 = (Y(NM1)-Y(NM2))*SFacAvg             !b3
                    FDY1 = X(NM3)-X(NM2)                       !a1
                    FDY2 = X(NM1)-X(NM3)                       !a2
                    FDY3 = X(NM2)-X(NM1)                       !a3
                    ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*QX1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*QX1(NM3)
                    ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*QY1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*QY1(NM3)
                    ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                ENDIF
            END DO
            NM3=NM2                                          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle /= 0) .AND. (NCEle /= 0)) THEN
                FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                FDY1 = X(NM3)-X(NM2) !a1
                FDY2 = X(NM1)-X(NM3) !a2
                FDY3 = X(NM2)-X(NM1) !a3
                ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*QX1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*QX1(NM3)
                ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*QY1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*QY1(NM3)
                ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
            ENDIF
            IF(NCI == 0) THEN
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
        ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
            AUVXX(NBDI)=1.d0
            AUVXY(NBDI)=0.d0
#endif
#endif
        ENDIF

    ENDDO

!...
!...  SOLVE FOR FLUX AT NEW LEVEL  (K+1)
!...

!.....Note: This includes the comparison between MJU and NODELE to
!.....determine if the node is an interface node.  If MJU < NODELE the
!.....velocity can be zeroed out to obtain an essential zero velocity at
!.....interface nodes.

    DO I=1,NP
        AUV22=AUV11(I)
        AUV21=-AUV12(I)
        DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
        QX2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
        QY2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU
#ifdef CSWAN
#ifdef CSWANFRIC
    ! Casey 091020: Adopt Ethan's/Joannes's modified friction.
        AUVYY(I)= AUVXX(I)
        AUVYX(I)= - AUVXY(I)
        DDU=AUVXX(I)*AUVYY(I)-AUVXY(I)*AUVYX(I)
        QX2(I)=(MOM_LV_X(I)*AUVYY(I)-MOM_LV_Y(I)*AUVXY(I))/DDU
        QY2(I)=(MOM_LV_Y(I)*AUVXX(I)-MOM_LV_X(I)*AUVYX(I))/DDU
#endif
#endif
    !        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
    !           QX2(I)=0.D0    !no slip and normal flux
    !           QY2(I)=0.D0    !on wet/dry interface nodes
    !           ENDIF
    END DO

!...
!...  Impose a zero normal flux gradient based on interpolating the
!...  flux at a fictitious point in the interior of the domain,
!...  normal to a specified boundary node and setting the boundary
!...  flux equal to the interpolated value at the fictitious point.
!...  Provided the fictitious point does not lie in an element that
!...  contains a boundary point, this is an entirely implicit
!...  calculation.
!...
    IF(NFLUXGBC == 1) THEN
        DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I) == 40) THEN
                NM1=NM(NEleZNG(I),1)
                NM2=NM(NEleZNG(I),2)
                NM3=NM(NEleZNG(I),3)
                NC1=NODECODE(NM1)
                NC2=NODECODE(NM2)
                NC3=NODECODE(NM3)
                NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
                QX2(NBDI)=NCEle*(QX2(NM1)*ZNGIF1(I)+QX2(NM2)*ZNGIF2(I) &
                +QX2(NM3)*ZNGIF3(I))
                QY2(NBDI)=NCEle*(QY2(NM1)*ZNGIF1(I)+QY2(NM2)*ZNGIF2(I) &
                +QY2(NM3)*ZNGIF3(I))
            ENDIF
        ENDDO
    ENDIF

!...  Compute velocities

    DO I=1,NP
        H2=DP(I)+IFNLFA*ETA2(I)
        IF(H2 /= 0.) THEN
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
!**********************************************************************
    END SUBROUTINE MOM_EQS_NEW_CONSERV
!**********************************************************************


!*******************************************************************************
!                                                                              *
!    Subroutine to compute the elevation using the GWCE formluation            *
!    This subroutine is the corrector step for the predictor-corrector         *
!    algorithm and obtains the corrected elevations                            *
!    Re-written to conform to the ADCIRC Theory Report                         *
!                                                                              *
!                            k.d.  06/24/2004                                  *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

    SUBROUTINE GWCE_new_pc(IT,TimeLoc,TimeH)

    USE SIZES
    USE GLOBAL
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, TotalArea, Areas, NeiTab, &
    NNEIGH, SFAC
    USE BOUNDARIES, ONLY : NETA, NFLUXB, NFLUXF, NFLUXGBC, NFLUXRBC, &
    NFLUXIB, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3
    USE WIND
    USE ITPACKV
    USE NodalAttributes, ONLY : FRIC, Tau0Var, HBREAK, FTHETA, FGAMMA, &
    IFLINBF, IFNLBF, IFHYBF, EVM, LoadAdvectionState
#ifdef CMPI
    USE MESSENGER
#endif

    IMPLICIT NONE

    INTEGER :: IE, JN, IJ, I, J                           !local loop counters
    INTEGER :: IT
    INTEGER ::  NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI, NCJ
    INTEGER :: NCyc
    INTEGER :: NBDI
    INTEGER :: OnDiag, OffDiag

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

! md   Added in parameters for the pc algorithm
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

    call setMessageSource("gwce_new_pc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

! md    Must reset the result vector to zero before recomputing
! md  the next time level.
    DO I=1,NP
        GWCE_LV(I) =0.D0
    END DO

!     Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1
!     Re-compute these local values

    OnDiag=(1+ILump)*2                         !diagonal coefficient
    OffDiag=(1-ILump)                          !off diagonal coefficient

!...
!...  Compute the GWCE load vector GWCE_LV
!...  This is done primarily element by element by forming
!...  temporary vectors and then assembling at the end.
!...  This has been set up to unroll loops to optimize performance
!...  on vector processors.
!...
!...  Elevation and flux boundary conditions are imposed after the
!...  element by element assembly section.
!...

!...  Initialize variables to zero if these forcings are not used

    IF((NWS /= 0) .OR. (NRS /= 0)) THEN
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

!...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

    IF ((CGWCE_LS_2PartV) .OR. (CGWCE_LS_2PartSV)) THEN

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
                EVMSmag=EVMEle* &
                sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA) &
                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
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
            IF(TotalArea(I) /= 0.) THEN
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

!...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

    IF ((CGWCE_LS_2PartQ) .OR. (CGWCE_LS_2PartSQ)) THEN

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
                EVMSmag=EVMEle* &
                sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA) &
                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
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
            IF(TotalArea(I) /= 0.) THEN
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
    ! md
    ! md  Added in the three time levels for the tau term.
    ! md  Every term is updated for the three time levels.
    ! md
        UV0=SQRT(UU0(I)*UU0(I)+VV0(I)*VV0(I))
        UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
        UV2=SQRT(UU2(I)*UU2(I)+VV2(I)*VV2(I))
        H00=DP(I)+IFNLFA*ETA0(I)
        H1=DP(I)+IFNLFA*ETA1(I)
        H2=DP(I)+IFNLFA*ETA2(I)
        TK0(I)=FRIC(I)*(IFLINBF + (UV0/H00)*(IFNLBF + IFHYBF* &
        (1+(HBREAK/H00)**FTHETA)**(FGAMMA/FTHETA)))
        TK(I)=FRIC(I)*(IFLINBF + (UV1/H1)*(IFNLBF + IFHYBF* &
        (1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
        TK2(I)=FRIC(I)*(IFLINBF + (UV2/H2)*(IFNLBF + IFHYBF* &
        (1+(HBREAK/H2)**FTHETA)**(FGAMMA/FTHETA)))
    END DO

! md      Added in the time weights
!...     Time weights for the nonlinear terms in the GWCE for
!...       the corrector step
    timewtgwce0=0.33d0
    timewtgwce1=0.34d0
    timewtgwce2=0.33d0
    timeagflag=1.0d0

!...  Assemble the GWCE RHS except for the boundary integral terms
! md  Renumber the GWCE loop for the corrector step

    DO 1038 IE=1,NE

    !...     Set nodal values for each element
    ! md
    ! md  Define the needed product terms at three time levels
    ! md


    ! Corbitt 120322: Localized Advection
        IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN     !wind or radiation stress
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

    !...     Compute part of several spatial gradients for use below

        E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
        E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
        E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
        E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
        Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
        Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

    !...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM

        IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.d0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3) &
            -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3) &
            +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3) &
            -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
        ENDIF

    !...     Compute the remainder of the 2 Part lateral stress terms

        IF((CGWCE_LS_2PartQ) .OR. (CGWCE_LS_2PartV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3 &
            +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3 &
            +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
        ENDIF

    !...     Compute the spatial gradients of the velocity dispersion terms if 3D

        IF (C3D) THEN                         !3D bottom friction
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3 &
            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3 &
            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
        ENDIF

    !...     Compute elemental averages

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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
        ENDIF
        IF (C3D) THEN                !3D velocity dispersion
            DispXAvg=IFNLCT*DispX
            DispYAvg=IFNLCT*DispY
        ENDIF
        IF(CBaroclinic) THEN
            BCXAvg=(H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2) &
            +H1N3*VIDBCPDXOH(NM3))/3.d0
            BCYAvg=(H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2) &
            +H1N3*VIDBCPDYOH(NM3))/3.d0
        ENDIF

    !...     Compute additional partial factors

        MsFacR=AreaIE*(1.d0/DT-Tau0Avg/2.d0)/DT/12.d0
        GOAreaIE4=G/AreaIE4
        Tau0SpaVar0=(QX0Avg*Tau0XGrad2A+QY0Avg*Tau0YGrad2A)/6.d0
        Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
        Tau0SpaVar2=(QX2Avg*Tau0XGrad2A+QY2Avg*Tau0YGrad2A)/6.d0
        A00pB00=A00+B00

    !...     Compute the JX, JY terms less the advection terms
        JXAvg = timewtgwce0*CorifAvg*QY0Avg+timewtgwce1*CorifAvg &
        *QY1Avg+timewtgwce2*CorifAvg*QY2Avg &
        -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDX1+E0N2SQ*FDX2 &
        +E0N3SQ*FDX3) &
        -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDX1+E1N2SQ*FDX2 &
        +E1N3SQ*FDX3) &
        -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDX1+E2N2SQ*FDX2 &
        +E2N3SQ*FDX3) &
        -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1 &
        +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3) &
        +WSXAvg-timewtgwce0*BSX0Avg-timewtgwce1*BSXAvg &
        -timewtgwce2*BSX2Avg &
        +MXAvg-DispXAvg-BCXAvg &
        +timewtgwce0*Tau0QX0Avg+timewtgwce1*Tau0QXAvg &
        +timewtgwce2*Tau0QX2Avg
        JYAvg =-timewtgwce0*CorifAvg*QX0Avg-timewtgwce1*CorifAvg &
        *QX1Avg-timewtgwce2*CorifAvg*QX2Avg &
        -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDY1+E0N2SQ*FDY2 &
        +E0N3SQ*FDY3) &
        -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDY1+E1N2SQ*FDY2 &
        +E1N3SQ*FDY3) &
        -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDY1+E2N2SQ*FDY2 &
        +E2N3SQ*FDY3) &
        -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1 &
        +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3) &
        +WSYAvg-timewtgwce0*BSY0Avg-timewtgwce1*BSYAvg &
        -timewtgwce2*BSY2Avg &
        +MYAvg-DispYAvg-BCYAvg &
        +timewtgwce0*Tau0QY0Avg+timewtgwce1*Tau0QYAvg &
        +timewtgwce2*Tau0QY2Avg

    !...     Complete the JX, JY terms depending on the advection formulation
        IF(CGWCE_Advec_NC) THEN        !nonconservative advection
            JXAvg = JXAvg - IFNLCT*timewtgwce0*( &
            QX0Avg*(U0N1*FDX1+U0N2*FDX2+U0N3*FDX3) &
            +QY0Avg*(U0N1*FDY1+U0N2*FDY2+U0N3*FDY3))/AreaIE2 &
            - IFNLCT*timewtgwce1 &
            *(QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3) &
            +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2 &
            - IFNLCT*timewtgwce2 &
            *(QX2Avg*(U2N1*FDX1+U2N2*FDX2+U2N3*FDX3) &
            +QY2Avg*(U2N1*FDY1+U2N2*FDY2+U2N3*FDY3))/AreaIE2 &
            +IFNLCAT*(timewtgwce0*U0Avg+timewtgwce1*U1Avg &
            +timewtgwce2*U2Avg) &
            *(timeagflag*0.5d0*((E2N1-E0N1+E2N2-E0N2 &
            +E2N3-E0N3)/DT) &
            +(1.D0-timeagflag)*ESAvg/DT)
            JYAvg = JYAvg - IFNLCT*timewtgwce0*( &
            QX0Avg*(V0N1*FDX1+V0N2*FDX2+V0N3*FDX3) &
            +QY0Avg*(V0N1*FDY1+V0N2*FDY2+V0N3*FDY3))/AreaIE2 &
            - IFNLCT*timewtgwce1 &
            *(QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3) &
            +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2 &
            - IFNLCT*timewtgwce2 &
            *(QX2Avg*(V2N1*FDX1+V2N2*FDX2+V2N3*FDX3) &
            +QY2Avg*(V2N1*FDY1+V2N2*FDY2+V2N3*FDY3))/AreaIE2 &
            +IFNLCAT*(timewtgwce0*V0Avg+timewtgwce1*V1Avg &
            +timewtgwce2*V2Avg) &
            *(timeagflag*0.5d0*((E2N1-E0N1+E2N2-E0N2 &
            +E2N3-E0N3)/DT) &
            +(1.D0-timeagflag)*ESAvg/DT)
        ENDIF
        IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
            JXAvg = JXAvg - IFNLCT*( &
            (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2 &
            +U1N3*QX1N3*FDX3) &
            +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2 &
            +U1N3*QY1N3*FDY3))/AreaIE2
            JYAvg = JYAvg - IFNLCT*( &
            (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2 &
            +V1N3*QX1N3*FDX3) &
            +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2 &
            +V1N3*QY1N3*FDY3))/AreaIE2
        ENDIF
        IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
            JXAvg = JXAvg - IFNLCT*( &
            QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3) &
            +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3) &
            +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
            JYAvg = JYAvg - IFNLCT*( &
            QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3) &
            +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3) &
            +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3) &
            +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
        ENDIF


    !...     Assemble forcing for node NM1 (local index j=1)

        Temp_LV_A1= &

    !...     Transient and Tau0 terms from LHS
        (OnDiag*ESN1 + OffDiag*(ESN2+ESN3))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX1*E0XGrad2A+FDY1*E0YGrad2A) &
        +A00pB00*(FDX1*E1XGrad2A+FDY1*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX1 + JYAvg*FDY1)/2.d0 &

    !...     Spatially varying Tau0 terms
        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar &
        +timewtgwce2*Tau0SpaVar2


    !...     Assemble forcing for node NM2 (local index j=2)

        Temp_LV_A2= &

    !...     Transient and Tau0 terms from LHS
        (OnDiag*ESN2 + OffDiag*(ESN1+ESN3))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX2*E0XGrad2A+FDY2*E0YGrad2A) &
        +A00pB00*(FDX2*E1XGrad2A+FDY2*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX2 + JYAvg*FDY2)/2.d0 &

    !...     Spatially varying Tau0 terms
        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar &
        +timewtgwce2*Tau0SpaVar2


    !...     Assemble forcing for node NM3 (local index j=3)

        Temp_LV_A3= &

    !...     Transient and Tau0 terms from LHS
    !...    (consistent mass matrix: ILump=0, lumped mass matrix: ILump=1)
        (OnDiag*ESN3 + OffDiag*(ESN1+ESN2))*MsFacR &

    !...     Free surface terms from LHS (time levels s-1 & s)
        -GDPAvgOAreaIE4*(  C00  *(FDX3*E0XGrad2A+FDY3*E0YGrad2A) &
        +A00pB00*(FDX3*E1XGrad2A+FDY3*E1YGrad2A)) &

    !...     Terms from momentum eqs.
        +(JXAvg*FDX3 + JYAvg*FDY3)/2.d0 &

    !...     Spatially varying Tau0 terms
        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar &
        +timewtgwce2*Tau0SpaVar2


    !...     Put these partial products into further elemental storage for a vector computer
    !...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
        Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
        Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
        Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

    !...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
        GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
        GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
        GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif

    !        IF(IE.EQ.1) THEN
    !           WRITE(101,*) ' '
    !           WRITE(101,*) '  ************* GWCE Load Vector ************'
    !           WRITE(101,*) '  Time Step = ',IT
    !           ENDIF
    !        WRITE(101,*) IE, ESN1, ESN2, ESN3
    !        WRITE(101,*) IE, TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3

    ! md  Change the number of the loop
        CONTINUE      !End of elemental loop
    1038 END DO


!...  Put load vector elemental values into nodal storage for a vector computer
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


!...  Save the elevation at the past time step into Eta1 and zero out Eta2
! md
! md  Already did this and haven't advanced in time yet, so just
! md  comment out the save elevation line. But we do need to zero
! md  out eta2 because it's involved in the summation below.
! md
    DO I=1,NP
    !         Eta1(I)=Eta2(I)
        Eta2(I)=0.0d0
    END DO

!...  At elevation boundary condition nodes, determine the elevation at
!...  the s+1 time step
!...
!...  For periodic elevation boundary conditions

    DO J=1,NBFR
        IF(PER(J) == 0.) THEN
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

!...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

    IF((NBFR == 0) .AND. (NOPE > 0)) THEN
        IF(TimeLoc > ETIME2) THEN
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
            Eta2(NBDI)=RampElev &
            *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
        END DO
    ENDIF


!...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
!...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

!...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
!...  assume that QN is positive into the domain.  QFORCEJ has a -1
!...  built in and the terms are not explicitly negated. Boundary
!...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
!...  flux that  is positive out of the domain.  Therefore, to match
!...  the formulation these terms must be explicitly multiplied by -1.

!...Note 3, Eta1 is the latest computed elevation (it was updated above).

    IF((NFLUXF == 1) .OR. (NFLUXB == 1) .OR. (NFLUXIB == 1) &
     .OR. (NFLUXGBC == 1) .OR. (NFLUXRBC == 1)) THEN
        NBDJ=NBV(1)
        IF(LBCODEI(1) <= 29) QFORCEJ=(QN2(1)-QN0(1))/DT2 + &
        Tau0VAR(NBDJ)*QN1(1)

        IF(LBCODEI(1) == 30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
        ENDIF

        IF(LBCODEI(1) == 32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT &
            -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT &
            +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
        ENDIF

        IF((LBCODEI(1) == 40) .OR. (LBCODEI(1) == 41)) QFORCEJ= &
        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

        DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J) <= 29) QFORCEJ=(QN2(J)-QN0(J))/DT2+ &
            Tau0VAR(NBDJ)*QN1(J)

            IF(LBCODEI(J) == 30) THEN
                H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
                CELERITY=SQRT(G*H1)
                QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
            ENDIF

            IF(LBCODEI(J) == 32) THEN
                H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
                CELERITY=SQRT(G*H1)
                QFORCEJ=(QN1(J)-QN0(J))/DT &
                -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT &
                +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
            ENDIF

            IF((LBCODEI(J) == 40) .OR. (LBCODEI(J) == 41)) QFORCEJ= &
            -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI) &
            + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ) &
            + BndLenO6NC*(2.d0*QForceJ+QForceI)
        ENDDO
    ENDIF

!...
!...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
!...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
!...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
!...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
!...  EQUATIONS
!...
    DO I=1,NETA
        NBDI=NBD(I)
        ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
        GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
        DO J=2,NNEIGH(NBDI)
            GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J)) &
            -ETAS(NBDI)*OBCCOEF(I,J-1)
        END DO
    END DO

!...
!...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
!...

!...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
!...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
    CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
#endif

!...  JCG ITERATIVE MATRIX SOLVER
    IPARM(1)=ITMAX
    CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS, &
    IWKSP,NW,WKSP,IPARM,RPARM,IER)

    NUMITR=IPARM(1)
    DO I=1,NP
        ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
    END DO

!     UPDATE ELEVATIONS

#ifdef CMPI
    CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif

! md
! md  Eta2 values are now corrected elevations at time level s+1.
! md
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!**********************************************************************
    END SUBROUTINE GWCE_NEW_PC
!**********************************************************************

!*******************************************************************************
!                                                                              *
!   Subroutine to compute the velocity and from that the flux/unit width using *
!   a 2DDI non conservative momentum equation.                                 *
!                                                                              *
!   Options are provided for either the correct area integration or the        *
!   original incorrect area integration.                                       *
!                                                                              *
!   Options are provided to use either flux or velocity based lateral          *
!   viscosity.                                                                 *
!                                                                              *
!   For a uniform grid and velocity based lateral viscosity, this subroutine   *
!   should give the same results as the original nonconservative formulation.  *
!                                                                              *
!   This subroutine follows the naming convention and formulation in the new   *
!   ADCIRC theory report.                                                      *
!                                                                              *
!   This subroutine provides the corrector part of the momentum equation for   *
!   the predictor-corrector algorithm and obtains the corrected velocities.    *
!                                                                              *
!                            k.d.  06/24/2004                                  *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

    SUBROUTINE Mom_Eqs_Non_Conserv_pc()

    USE GLOBAL
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, TotalArea, MJU, &
    NEITAB, NEITABELE, NNEIGH, SFAC
    USE BOUNDARIES, ONLY : NFLUXGBC, NVELME, ME2GW, NBV, LBCODEI, &
    CSII, SIII, ME2GW, NBV, LBCODEI, NEleZNG, ZNGIF1, ZNGIF2, ZNGIF3
    USE WIND
    USE NodalAttributes, ONLY: EVM,LoadAdvectionState

    IMPLICIT NONE

    INTEGER :: IE, I, J, N                           !local loop counters
    INTEGER ::  NM1, NM2, NM3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    INTEGER :: NBDI
    INTEGER :: NNFirst

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

! md   Added in parameters for the pc algorithm
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

    call setMessageSource("mom_eqs_non_conserv_pc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif
!...
!...  UPDATE LOAD VECTOR MOM_LV_X(I) AND MOM_LV_Y(I)
!...  NOTE: MOM_LV_X, MOM_LV_Y AND AUV ARE ZEROED OUT AT THE TOP OF
!...        THE TIME STEPPING LOOP.
! md    Must reset the result vectors to zero before recomputing
! md  the next time level.
    DO I=1,NP
        MOM_LV_X(I)=0.D0
        MOM_LV_Y(I)=0.D0
    END DO

! md  Add in the time weights for the corrector loop
    timewtmom0=0.0d0
    timewtmom1=0.5d0
    timewtmom2=0.5d0
    timebfflag=1.0d0

!.....FIRST TREAT THE NON-LUMPED PART OF THE EQUATIONS.

    DO IE=1,NE

    !...  SET NODAL VALUES FOR EACH ELEMENT


    ! Corbitt 120322: Localized Advection
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

    !...  Compute element averaged quantities

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
            EVMSmag=EVMEle* &
            sqrt((DU1DXA-DV1DYA)*(DU1DXA-DV1DYA) &
            +(DU1DYA+DV1DXA)*(DU1DYA+DV1DXA))
            EVMEle=EVMSmag
        ENDIF

    !...  Compute terms associated with the barotropic pressure

        BTP1N1=ETA1(NM1)
        BTP2N1=ETA2(NM1)
        BTP1N2=ETA1(NM2)
        BTP2N2=ETA2(NM2)
        BTP1N3=ETA1(NM3)
        BTP2N3=ETA2(NM3)

    !.......If using atm pressure add it into the barotropic pressure

        IF(NWS /= 0) THEN
            BTP1N1=BTP1N1+PR1(NM1)
            BTP2N1=BTP2N1+PR2(NM1)
            BTP1N2=BTP1N2+PR1(NM2)
            BTP2N2=BTP2N2+PR2(NM2)
            BTP1N3=BTP1N3+PR1(NM3)
            BTP2N3=BTP2N3+PR2(NM3)
        ENDIF

    !.......If using tidal potential terms, add these into the barotropic pressure

        IF (CTIP) THEN
            BTP1N1=BTP1N1-TiP1(NM1)
            BTP2N1=BTP2N1-TiP2(NM1)
            BTP1N2=BTP1N2-TiP1(NM2)
            BTP2N2=BTP2N2-TiP2(NM2)
            BTP1N3=BTP1N3-TiP1(NM3)
            BTP2N3=BTP2N3-TiP2(NM3)
        ENDIF

    !...  Compute the barotropic pressure gradient x area for the element

        DBTPDXA=((BTP1N1*FDX1+BTP1N2*FDX2+BTP1N3*FDX3) &
        +(BTP2N1*FDX1+BTP2N2*FDX2+BTP2N3*FDX3))/2.d0
        DBTPDYA=((BTP1N1*FDY1+BTP1N2*FDY2+BTP1N3*FDY3) &
        +(BTP2N1*FDY1+BTP2N2*FDY2+BTP2N3*FDY3))/2.d0

    !...  Compute the advective term gradients x area for the element

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

    !...  Compute the lateral viscous terms for the element (flux formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric flux formulation)

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

    !...  Compute the lateral viscous terms for the element (velocity formulation)

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

    !...  Compute the lateral viscous terms for the element (symmetric velocity formulation)

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

    
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM1

        TEMP_LV_A1=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN1*FDX1+LSXYN1*FDY1)/H1N1 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM2
    !...
        TEMP_LV_A2=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN2*FDX2+LSXYN2*FDY2)/H1N2 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...
    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR X-MOMENTUM EQUATION INTO
    !...  TEMP_LV_A VECTOR FOR NODE NM3
    !...
        TEMP_LV_A3=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDU0DXA+V0AvgDU0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDU1DXA+V1AvgDU1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDU2DXA+V2AvgDU2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDXA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSXXN3*FDX3+LSXYN3*FDY3)/H1N3 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM1

        TEMP_LV_B1=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN1*FDX1+LSYYN1*FDY1)/H1N1 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM2

        TEMP_LV_B2=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN2*FDX2+LSYYN2*FDY2)/H1N2 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !...  LOAD NON-LUMPED ELEMENTAL COMPONENTS FOR Y-MOMENTUM EQUATION INTO
    !...  TEMP_LV_B VECTOR FOR NODE NM3

        TEMP_LV_B3=NCEle*DT*( &
    !...  ADVECTIVE TERMS
        -IFNLCT*timewtmom0*(U0AvgDV0DXA+V0AvgDV0DYA) &
        -IFNLCT*timewtmom1*(U1AvgDV1DXA+V1AvgDV1DYA) &
        -IFNLCT*timewtmom2*(U2AvgDV2DXA+V2AvgDV2DYA) &
    !...  BAROTROPIC PRESSURE GRADIENT
        -GO2*DBTPDYA &
    !...  LATERAL VISCOUS TERMS
        -1.5d0*(LSYXN3*FDX3+LSYYN3*FDY3)/H1N3 &
    !...  STILL NEED TO DIVIDE BY TOTAL AREA AROUND A NODE
        )

    !     Original (incorrect) area integration - for historical comparison

        IF (CME_AreaInt_Orig) THEN
            TEMP_LV_A1=TEMP_LV_A1/AreaIE
            TEMP_LV_A2=TEMP_LV_A2/AreaIE
            TEMP_LV_A3=TEMP_LV_A3/AreaIE
            TEMP_LV_B1=TEMP_LV_B1/AreaIE
            TEMP_LV_B2=TEMP_LV_B2/AreaIE
            TEMP_LV_B3=TEMP_LV_B3/AreaIE
        ENDIF

    !     LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
        TEMP_LV_A(IE,1)=TEMP_LV_A1
        TEMP_LV_A(IE,2)=TEMP_LV_A2
        TEMP_LV_A(IE,3)=TEMP_LV_A3
        TEMP_LV_B(IE,1)=TEMP_LV_B1
        TEMP_LV_B(IE,2)=TEMP_LV_B2
        TEMP_LV_B(IE,3)=TEMP_LV_B3
#endif

    !     LINES TO RUN ON A SCALAR COMPUTER
    !     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
    !           AND QUV ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
        MOM_LV_X(NM1)=MOM_LV_X(NM1)+TEMP_LV_A1
        MOM_LV_X(NM2)=MOM_LV_X(NM2)+TEMP_LV_A2
        MOM_LV_X(NM3)=MOM_LV_X(NM3)+TEMP_LV_A3
        MOM_LV_Y(NM1)=MOM_LV_Y(NM1)+TEMP_LV_B1
        MOM_LV_Y(NM2)=MOM_LV_Y(NM2)+TEMP_LV_B2
        MOM_LV_Y(NM3)=MOM_LV_Y(NM3)+TEMP_LV_B3
#endif

    ENDDO


!     LINES TO RUN ON A VECTOR COMPUTER
!     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR MOM_LV_X, MOM_LV_Y
!           AND AUV
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


!...  Update the momentum equation LHS coefficients and load vectors at each
!...  node by dividing by the area of all active elements attached to the node
!...  and adding in the lumped terms, bottom friction and boundary conditions

    WSX=0.D0
    WSY=0.D0
    VIDBCPDX=0.D0
    VIDBCPDY=0.D0
    DO I=1,NP
        NCI=NODECODE(I)
        IF(TotalArea(I) /= 0.d0) THEN
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
        IF((NWS /= 0) .OR. (NRS /= 0)) THEN
            WSX=DTO2*IFWIND*(WSX1(I)/H1+WSX2(I)/H2)
            WSY=DTO2*IFWIND*(WSY1(I)/H1+WSY2(I)/H2)
        ENDIF
    ! md
    ! md  Added time weights to tau terms in the momentum equation here.
    ! md  Note the weighting should follow the weighting in the GWCE.
    ! md
        VCoef1=DTO2*TK(I)                     !TK = Kslip/H
        VCoef12=DTO2*(TK2(I)*timebfflag+TK(I)*(1.D0-timebfflag))
        VCoef2=DTO2*CORIF(I)
        IF(CBaroclinic) THEN
            VIDBCPDX=DT*VIDBCPDXOH(I)
            VIDBCPDY=DT*VIDBCPDYOH(I)
        ENDIF

        MOM_LV_X(I)=NCI*(MOM_LV_X(I)+WSX+(1.D0-VCoef1)*UU1(I) &
        +VCoef2*VV1(I)-VIDBCPDX)
        MOM_LV_Y(I)=NCI*(MOM_LV_Y(I)+WSY+(1.D0-VCoef1)*VV1(I) &
        -VCoef2*UU1(I)-VIDBCPDY)

    ! md    Change for the corrector formulation
        AUV11(I)=1.D0+VCoef12*NCI
        AUV12(I)=-VCoef2*NCI
    END DO

!...  Modify the momentum equations to impose velocity boundary
!...  conditions In each case the equations are manipulated to
!...  maintain the LHS matrix structure of AUV11=AUV22;
!...  AUV12=-AUV21)

    DO J=1,NVELME
        I=ME2GW(J)
        NBDI=NBV(I)
        H2=DP(NBDI)+IFNLFA*ETA2(NBDI)
        NCI=NODECODE(NBDI)

    !      Specified essential normal flow and free tangential slip

        IF((LBCODEI(I) >= 0) .AND. (LBCODEI(I) <= 9)) THEN
            VelNorm=-QN2(I)/H2
            MOM_LV_X(NBDI)=(SIII(I)*MOM_LV_X(NBDI) &
            -CSII(I)*MOM_LV_Y(NBDI) &
            -VelNorm*AUV12(NBDI))*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*AUV11(NBDI)*NCI   !Normal Eqn RHS
            AUV12(NBDI)=-CSII(I)*AUV11(NBDI)
            AUV11(NBDI)=SIII(I)*AUV11(NBDI)
        ENDIF

    !     Specified essential normal flow and no tangential slip

        IF((LBCODEI(I) >= 10) .AND. (LBCODEI(I) <= 19)) THEN
            VelNorm=-QN2(I)/H2
            VelTan=0.D0
            MOM_LV_X(NBDI)=VelTan*NCI !Tangential Eqn RHS
            MOM_LV_Y(NBDI)=VelNorm*NCI !Normal Eqn RHS
            AUV11(NBDI)=SIII(I)
            AUV12(NBDI)=-CSII(I)
        ENDIF

    !     Zero normal velocity gradient using a Galerkin approximation to
    !     the normal derivatives. Note: this is fully explicit and therefore
    !     the velocity at the boundary is computed entirely from surrounding
    !     velocities at the previous time step.

        IF(LBCODEI(I) == 41) THEN
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
                NEle=NeiTabEle(NBDI,N-2) !element# defined by nodes NM1,NM2,NM3
                NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NEle)
                IF((NEle /= 0) .AND. (NCEle /= 0)) THEN !if element is active, compute contribution
                    FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                    FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                    FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                    FDY1 = X(NM3)-X(NM2) !a1
                    FDY2 = X(NM1)-X(NM3) !a2
                    FDY3 = X(NM2)-X(NM1) !a3
                    ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                    ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                    -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                    ZNGLHS =ZNGLHS  +  CSII(I)*FDX1+SIII(I)*FDY1
                ENDIF
            END DO
            NM3=NM2          !wrap back to beginning to get final contribution
            NM2=NNFirst
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.d0
            NEle=NeiTabEle(NBDI,NNeigh(NBDI)-1)
            NCEle=NCI*NodeCode(NM2)*NodeCode(NM3)*NOFF(NELE)
            IF((NEle /= 0) .AND. (NCEle /= 0)) THEN
                FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1
                FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2
                FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3
                FDY1 = X(NM3)-X(NM2) !a1
                FDY2 = X(NM1)-X(NM3) !a2
                FDY3 = X(NM2)-X(NM1) !a3
                ZNGRHS1=ZNGRHS1-(CSII(I)*FDX2+SIII(I)*FDY2)*UU1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*UU1(NM3)
                ZNGRHS2=ZNGRHS2-(CSII(I)*FDX2+SIII(I)*FDY2)*VV1(NM2) &
                -(CSII(I)*FDX3+SIII(I)*FDY3)*VV1(NM3)
                ZNGLHS =ZNGLHS + CSII(I)*FDX1+SIII(I)*FDY1
            ENDIF
            IF(NCI == 0) THEN
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

!...
!...  SOLVE FOR VELOCITY AT NEW LEVEL  (s+1)
!...

!.....Note: This includes the comparison between MJU and NODELE to
!.....determine if the node is an interface node.  If MJU < NODELE the
!.....velocity can be zeroed out to obtain an essential zero velocity at
!.....interface nodes.

    DO I=1,NP
        AUV22=AUV11(I)
        AUV21=-AUV12(I)
        DDU=AUV11(I)*AUV22-AUV12(I)*AUV21
        UU2(I)=(MOM_LV_X(I)*AUV22-MOM_LV_Y(I)*AUV12(I))/DDU
        VV2(I)=(MOM_LV_Y(I)*AUV11(I)-MOM_LV_X(I)*AUV21)/DDU

    !        IF(MJU(I).NE.NODELE(I)) THEN !uncomment for essential
    !           UBAR2(I)=0.D0    !no slip and normal flux
    !           VBAR2(I)=0.D0    !on wet/dry interface nodes
    !           ENDIF
    END DO

!...
!...  Impose a zero normal velocity gradient based on interpolating the
!...  velocity at a fictitious point in the interior of the domain,
!...  normal to a specified boundary node and setting the boundary
!...  velocity equal to the interpolated value at the fictitious point.
!...  Provided the fictitious point does not lie in an element that
!...  contains a boundary point, this is an entirely implicit
!...  calculation.
!...
    IF(NFLUXGBC == 1) THEN
        DO J=1,NVELME
            I=ME2GW(J)
            NBDI=NBV(I)
            IF(LBCODEI(I) == 40) THEN
                NM1=NM(NEleZNG(I),1)
                NM2=NM(NEleZNG(I),2)
                NM3=NM(NEleZNG(I),3)
                NC1=NODECODE(NM1)
                NC2=NODECODE(NM2)
                NC3=NODECODE(NM3)
                NCEle=NC1*NC2*NC3*NOFF(NEleZNG(I))
                UU2(NBDI)=NCEle*(UU2(NM1)*ZNGIF1(I)+UU2(NM2)*ZNGIF2(I) &
                +UU2(NM3)*ZNGIF3(I))
                VV2(NBDI)=NCEle*(VV2(NM1)*ZNGIF1(I)+VV2(NM2)*ZNGIF2(I) &
                +VV2(NM3)*ZNGIF3(I))
            ENDIF
        ENDDO
    ENDIF

!...  Compute fluxes

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
!**********************************************************************
    END SUBROUTINE MOM_EQS_NON_CONSERV_PC
!**********************************************************************

!****************************************************************************************
!   Subroutine to compute Scalar Transport                                              *
!                                                                                       *
!   Note, this is not set up for parallel operation.                                    *
!****************************************************************************************

    SUBROUTINE SCALAR_TRANS_2D (IT,TimeLoc)

    USE GLOBAL
    USE MESH, ONLY : NE, NP, NM, DP, X, Y, AREAS, SFAC
    USE NodalAttributes, ONLY: EVC
#ifdef CMPI
    USE MESSENGER
#endif

    IMPLICIT NONE

    INTEGER :: IE, I                         !local loop counters
    INTEGER :: IT
    INTEGER :: NM1, NM2, NM3
    INTEGER :: NC1, NC2, NC3, NCEle, NCI

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
!...  NOTE: THE VARIABLE CH1(I) IS ACTUALLY C*H
!.... COMPUTE SOURCE/SINK TERM AT THE NODES USING CLASSICAL COHESIVE
!.... SEDIMENT TRANSPORT RELATIONS

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

        WSMOD=WS
        IF(C1 > CCRITD) WSMOD=WS*(C1/CCRITD)**1.33d0
        HSD=0.d0
        IF(BEDSTR < CBEDSTRD) HSD=-(2.d0*WSMOD*C1)* &
        (1.0d0-BEDSTR/CBEDSTRD)
        IF(HSD > 0.d0) HSD=0.d0

    !.....Calculate the surface erosion rate for cohesive sediment using
    !.....the Ariathurai et at. (1977) adaption of Partheniades' (1962) findings

        HSE=0.
        IF(BEDSTR > CBEDSTRE) HSE=ECONST*(BEDSTR/CBEDSTRE-1.0)

    !.....Determine the total source sink term

        SOURSIN(I)=HSD+HSE
    END DO

!.... UPDATE THE TRANSPORT EQUATION ELEMENT BY ELEMENT BY FORMING
!.... TEMPORARY VECTORS AND THEN ASSEMBLING.  NOTE: TRANS_LV_B(I), TRANS_LV_A(I) ARE
!.... ZEROED OUT AT THE TOP OF THE TIME STEPPING LOOP.  AGAIN THESE
!.... LOOPS HAVE BEEN UNROLLED TO OPTIMIZE VECTORIZATION

    DO IE=1,NE

    !.....SET NODAL VALUES FOR EACH ELEMENT

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

    !.....COMPUTE ELEMENTAL MATRICIES

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

    !.....COMPUTE ELEMENTAL QUANTITIES

        UEA=(U1N1+U1N2+U1N3)/3.
        VEA=(V1N1+V1N2+V1N3)/3.
        HEA=(H1N1+H1N2+H1N3)/3.
        EVCEA=(EVC1+EVC2+EVC3)/3.
        DHDX=H1N1*FDX1O2A+H1N2*FDX2O2A+H1N3*FDX3O2A
        DHDY=H1N1*FDY1O2A+H1N2*FDY2O2A+H1N3*FDY3O2A
        UPEA=UEA+DHDX*EVCEA/HEA
        VPEA=VEA+DHDY*EVCEA/HEA

    !.....ASSEMBLE PARTIAL PRODUCT

        CHSUM=CH1N1+CH1N2+CH1N3

    !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A1 AND
    !.....TEMP_LV_B1 VECTORS FOR NODE NM1

        TEMP_LV_B1=              & !LOAD VECTOR
    !......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
        FDDDODT*CH1N1+FDDODODT*(CH1N2+CH1N3) &
    !......LATERAL SGS TERMS
        -EVCEA*(DXXYY11*CH1N1+DXXYY12*CH1N2+DXXYY13*CH1N3) &
    !......ADVECTIVE TERMS
        +(UPEA*DDX1+VPEA*DDY1)*CHSUM &
    !......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
        +FDDD*SS1N1+FDDOD*(SS1N2+SS1N3)
        TEMP_LV_A1=              & !LHS VECTOR
    !......TRANSIENT TERM (LUMPED)
        FDDDODT+2.*FDDODODT

    !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A2 AND
    !.....TEMP_LV_B2 VECTOR FOR NODE NM2

        TEMP_LV_B2=              & !LOAD VECTOR
    !......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
        FDDDODT*CH1N2+FDDODODT*(CH1N1+CH1N3) &
    !......LATERAL SGS TERMS
        -EVCEA*(DXXYY12*CH1N1+DXXYY22*CH1N2+DXXYY23*CH1N3) &
    !......ADVECTIVE TERMS
        +(UPEA*DDX2+VPEA*DDY2)*CHSUM &
    !......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
        +FDDD*SS1N2+FDDOD*(SS1N1+SS1N3)
        TEMP_LV_A2=              & !LHS VECTOR
    !......TRANSIENT TERM (LUMPED)
        FDDDODT+2.*FDDODODT

    !.....LOAD ELEMENTAL COMPONENTS FOR TRANSPORT EQUATION INTO TEMP_LV_A3 AND
    !.....TEMP_LV_B3 VECTOR FOR NODE NM3

        TEMP_LV_B3=              & !LOAD VECTOR
    !......TRANSIENT TERM (EITHER LUMPED OR CONSISTENT)
        FDDDODT*CH1N3+FDDODODT*(CH1N1+CH1N2) &
    !......LATERAL SGS TERMS
        -EVCEA*(DXXYY13*CH1N1+DXXYY23*CH1N2+DXXYY33*CH1N3) &
    !......ADVECTIVE TERMS
        +(UPEA*DDX3+VPEA*DDY3)*CHSUM &
    !......SOURCE SINK TERMS (EITHER LUMPED OR CONSISTENT)
        +FDDD*SS1N3+FDDOD*(SS1N1+SS1N2)
        TEMP_LV_A3=              & !LHS VECTOR
    !......TRANSIENT TERM (LUMPED)
        FDDDODT+2.*FDDODODT

    !     VEC...LINES TO RUN ON A VECTOR COMPUTER
#ifdef CVEC
        TEMP_LV_B(IE,1)=TEMP_LV_B1*NCELE !LOAD VECTOR
        TEMP_LV_B(IE,2)=TEMP_LV_B2*NCELE !LOAD VECTOR
        TEMP_LV_B(IE,3)=TEMP_LV_B3*NCELE !LOAD VECTOR
        TEMP_LV_A(IE,1)=TEMP_LV_A1*NCELE !LUMPED LHS MATRIX
        TEMP_LV_A(IE,2)=TEMP_LV_A2*NCELE !LUMPED LHS MATRIX
        TEMP_LV_A(IE,3)=TEMP_LV_A3*NCELE !LUMPED LHS MATRIX
#endif

    !     LINES TO RUN ON A SCALAR COMPUTER
    !     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR QC AND TRANS_LV_A
    !     ON A SCALAR COMPUTER USING THE TEMPORARY VECTORS
#ifdef CSCA
        TRANS_LV_B(NM1)=TRANS_LV_B(NM1)+TEMP_LV_B1*NCELE !LOAD VECTOR
        TRANS_LV_B(NM2)=TRANS_LV_B(NM2)+TEMP_LV_B2*NCELE !LOAD VECTOR
        TRANS_LV_B(NM3)=TRANS_LV_B(NM3)+TEMP_LV_B3*NCELE !LOAD VECTOR
        TRANS_LV_A(NM1)=TRANS_LV_A(NM1)+TEMP_LV_A1*NCELE !LUMPED LHS MATRIX
        TRANS_LV_A(NM2)=TRANS_LV_A(NM2)+TEMP_LV_A2*NCELE !LUMPED LHS MATRIX
        TRANS_LV_A(NM3)=TRANS_LV_A(NM3)+TEMP_LV_A3*NCELE !LUMPED LHS MATRIX
#endif

    ENDDO

!     LINES TO RUN ON A VECTOR COMPUTER
!     NOTE: THESE LINES FINALIZE THE ASSEMBLY PROCESS FOR QC, TRANS_LV_A
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

!.... SOLVE FOR C*H NODE BY NODE

    DO I=1,NP
        NCI=NODECODE(I)
        IF(NCI /= 0) CH1(I)=TRANS_LV_B(I)/TRANS_LV_A(I)
    !     IF(LBArray_Pointer(I).NE.0) CH1(I)=0.d0  !ESSENTIAL C=0 BOUNDARY CONDITION
    END DO

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!**********************************************************************
    END SUBROUTINE SCALAR_TRANS_2D
!**********************************************************************

!****************************************************************************************
!   Subroutine to compute 2D SigmaT fields from 3D salinity and/or temperature fields   *
!                                                                                       *
!                                    R.L.  6/22/05                                      *
!****************************************************************************************


    SUBROUTINE CALC_SIGMAT_2D ()

    USE GLOBAL

    INTEGER :: NH

    call setMessageSource("calc_sigmat_2D")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif


    IF(ABS(IDEN) == 2) THEN
        DO NH=1,NP
        ENDDO
    ELSEIF(ABS(IDEN) == 3) THEN
        DO NH=1,NP
        ENDDO
    ELSEIF(ABS(IDEN) == 4) THEN
        DO NH=1,NP
        ENDDO
    ENDIF
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!*********************************************************************
    END SUBROUTINE CALC_SIGMAT_2D
!*********************************************************************




!******************************************************************************
! endra: Eliminated the bpg calculation in the vsmy subroutine and added a
!        subroutine for the bpg calculation v45.12
!******************************************************************************
!  SUBROUTINE BPG3D                                                           *
!                                                                             *
!  Note, the following time stepping coefficients are computed in             *
!     VSSTUP and passed in a common block.                                    *
!                                                                             *
!  IDTAlp1      = I*DelT*Alp1        - weights coriolis term in LHS matrix    *
!  IDT1MAlp1    = I*DelT*(1.-Alp1)   - weights coriolis term in RHS forcing   *
!  DTAlp3       = DelT*Alp3          - weights vert diff term in LHS matrix   *
!  DT1MAlp3     = DelT*(1-Alp3)      - weights vert diff term in RHS forcing  *
!  DTAlp2       = DelT*Alp2          - weights bot stress term in LHS matrix  *
!  DT1MAlp2     = DelT*(1.-Alp2)     - weights bot stress term in RHS forcing *
!                                                                             *
!  q(MNP,MNodes) - 3D Complex Velocity field (GAMMA) from past time step.     *
!                                                                             *
!                                                                             *
!  NH - horizontal node counter                                               *
!  NP - number of nodes in horizontal grid                                    *
!  NFEN - number of nodes in the vertical grid                                *
!  BTP - total barotropic pressure (atmos press, water level, tidal potential)*
!                 at time levels s+1/2                                        *
!******************************************************************************

    SUBROUTINE BPG3D()

! Casey: Added the following variable declarations from GLOBAL.

    USE SIZES, ONLY : SZ
    USE GLOBAL, ONLY: IFNLFA, IFNLCT, &
    NODECODE, NOFF, VIDBCPDXOH, VIDBCPDYOH, &
    setMessageSource, unsetMessageSource, allMessage, &
    logMessage, DEBUG, ECHO, INFO, WARNING, ERROR, &
    ETA2, SIGT0
    USE MESH, ONLY : NM, X, Y, NP, DP, AREAS, NODELE, &
    NEITABELE, NEITAB, NNEIGH, SFAC
    USE GLOBAL_3DVS, ONLY : BCP, BPG, NFEN, SIGT, SIGMA, IDEN, &
    GORhoOAMB, AMB, B, IY
#ifdef CMPI
    USE MESSENGER
    IMPLICIT NONE
    REAL(SZ) :: DUMV1(1),DUMV2(1)
#else
    IMPLICIT NONE
#endif

! Casey: Added the following local variable declarations.

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

    call setMessageSource("bpg3D")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!     INCREMENT THE TIMESTEP SINCE START COUNTER


!*************************************************************************************
!     Check whether it is time to print various 3D outputs


!     If a baroclinic run, compute the 3D baroclinic pressure field
!     The buoyancy field is defined as
!     BCP(z)    =(gravity/rho ref)*          integral (SigT) from surface down to z
!     BCP(sigma)=(gravity/rho ref)*(H/(a-b))*integral (SigT) from a down to sigma
!     where
!     SigT = Sigma T = Rho - 1000 = density - 1000
!     SigT0 = Sigma t value of reference density (typically = 0)
!     Sigma = dimensionless vertical coordinate

    IF((IDEN >= 1) .OR. (IDEN <= -1)) THEN
        DO NH=1,NP             !loop over horizontal nodes

        ! Casey: Changed "NolIFA" to "IFNLFA."
        
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
    !     Update BCP on ghost nodes
    !      CALL UPDATER3D(BCP) !!!!Don't know if this is needed at this time
#endif
    ENDIF


!*************************************************************************************
!     Compute 3D baroclinic pressure gradients


!     Loop over each horizontal node to compute the horizontal velocity

    DO NH=1,NP                !loop over horizontal nodes

        Hs  = DP(NH)+IFNLFA*Eta2(NH) !Total depth at previous (s) timestep
        HsOAMB=Hs/AMB

    !     Zero out baroclinic pressure gradient and vertically integrated
    !     baroclinic pressure gradient for a barotropic run

        IF (IDEN == 0) THEN
            DO k=1,NFEN
                BCPG(k)=(0.d0,0.d0)
            END DO
            VIDBCPDXOH(NH)=0.d0
            VIDBCPDYOH(NH)=0.d0
        ENDIF

    !     Start computing baroclinic terms

        IF ((IDEN >= 1) .OR. (IDEN <= -1)) THEN

        !     Start computing baroclinic pressure gradient (computed in level
        !     coordinates) at each node in the vertical

            DO k=1,NFEN

                DBCPDX2A=0.d0
                DBCPDY2A=0.d0
                TotalBCPGArea2=0.d0
                N1=NH
                BCPN1=BCP(NH,k)

                Zk=HsOAMB*(Sigma(k)-B)-DP(NH) !determine z corresponding to sigma level k
                N2=NEITAB(NH,2)  !operate on 1st neighbor

            ! Casey: Changed "NolIFA" to "IFNLFA."
            
                HsN2=DP(N2)+IFNLFA*Eta2(N2)
                SigmaNN=B+AMB*(Zk+DP(N2))/HsN2 !equivalent sigma value at neighbor
                CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
                NNFirst=N2       !save these values until end
                BCPNFirst=BCPN2  !save these values until end

                DO N=3,NNeigh(NH) !operate on rest of neighbors
                    N3=N2         !shift previously computed values
                    BCPN3=BCPN2   !shift previously computed values
                    N2=NeiTab(NH,N) !select new neighbor to work on

                ! Casey: Changed "NolIFA" to "IFNLFA."
                
                    HsN2=DP(N2)+IFNLFA*Eta2(N2)
                    SigmaNN=B+AMB*(Zk+DP(N2))/HsN2 !equivalent sigma value at neighbor
                    CALL ZSURFBUOY(SigmaNN,BCPN2,N2,k) !interp BCP at neighbor
                    NEle=NeiTabEle(NH,N-2) !element# defined by nodes NH,NN2,NN1
                ! jgf49.58: The NeiTabEle matrix is semi sparse; see
                ! Casey's comments in vsmy.F. NOFF array lookups fail
                ! if NEle comes back 0, so I will just cycle to the next
                ! N value here if that happens ... TODO: somebody please
                ! confirm that this is the right answer here.
                    IF (NEle == 0) THEN
                        CYCLE
                    ENDIF
                ! Casey: Added the computation of "NCELE" and the last part of the IF statement.
                
                    NCELE = NODECODE(NH)*NODECODE(N2) &
                    *NODECODE(N3)*NOFF(NELE)
                    IF((BCPN2 /= -999.) .AND. (BCPN3 /= -999.) &
                     .AND. (NEle /= 0) .AND. (NCELE /= 0)) THEN !if all 3 nodes are active, compute bu
                        TotalBCPGArea2=TotalBCPGArea2+Areas(NEle)
                    !    kmd48.33bc add in spherical factors
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

            ! Casey: Added the computation of "NCELE" and the last part of the IF statement.
            
            ! jgf49.58 NOFF lookups fail if NELE comes back 0.
                IF (NELE /= 0) THEN
                    NCELE = NODECODE(NH)*NODECODE(N2) &
                    *NODECODE(N3)*NOFF(NELE)
                ENDIF
                IF((BCPN2 /= -999.) .AND. (BCPN3 /= -999.) &
                 .AND. (NEle /= 0) .AND. (NCELE /= 0)) THEN
                    TotalBCPGArea2=TotalBCPGArea2+Areas(NEle)
                !    kmd48.33bc add in spherical factors
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

            !    kmd48.33bc changed the BPG calculation for bottom boundary issues
                IF(TotalBCPGArea2 == 0.) THEN
                !                  IF (k.eq.NFEN) THEN !kd46.01
                    BCPG(k)=(0.d0,0.d0)
                !                  ELSE
                !                     BCPG(k)=BCPG(k+1)
                !                  END IF
                ELSE
                    BCPG(k)=(DBCPDX2A+iy*DBCPDY2A)/TotalBCPGArea2
                ENDIF

            ENDDO

        !    kmd48.33bc added this to BPG
            DO k=NFEN,1,-1
                IF ((BCPG(k) == (0.d0,0.d0)) .AND. (k /= NFEN)) THEN
                    BCPG(k)=BCPG(k+1)
                END IF
            END DO


            DO k=1,NFEN
                BPG(NH,k) = BCPG(k)
            END DO

        !     Finished computing baroclinic pressure gradient (computed in level
        !     coordinates) at each node in the vertical


        !     Compute vertically integrated baroclinic pressure gradient for use
        !     in the wave equation.  NOTE: For a prognostic model in which the
        !     density field evolves in time, this calculation should be done
        !     after the new density field is computed.  In this case one would
        !     integrate over the vertical first and differentiate second.

            VIBCPG=(0.d0,0.d0)
            DO k=NFEN-1,1,-1
                VIBCPG=VIBCPG+0.5d0*(BCPG(k+1)+BCPG(k)) &
                *(Sigma(k+1)-Sigma(k))
            ENDDO
            VIDBCPDXOH(NH)=REAL(VIBCPG)/AMB
            VIDBCPDYOH(NH)=AIMAG(VIBCPG)/AMB

        ENDIF

    !     Finished computing baroclinic terms

    ENDDO

!     Finish loop over horizontal nodes to compute the horizontal velocity

#ifdef CMPI
!     Update new 3D baroclinic pressure gradient and the vertically
!     integrated baroclinic pressure gradient on ghost nodes

    CALL UPDATEC3D(BPG)
    CALL UPDATER(VIDBCPDXOH,VIDBCPDYOH,DUMV1,2)
#endif

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!***********************************************************************
    END SUBROUTINE BPG3D
!***********************************************************************


!*************************************************************************
!     Subroutine to interpolate baroclinic pressure (BCP) to a specified
!     sigma value (SigmaNN) given an initial guess of which sigma
!     level is closest to the specified value.

!                                    R.L.  5/04/01
!                                    R.L.  5.19/03
!*************************************************************************

    SUBROUTINE ZSURFBUOY(SigmaNN,BCPressNN,NN,J)

    USE GLOBAL_3DVS
    IMPLICIT NONE
    REAL(SZ) :: BCPressNN
    REAL(SZ) :: SigmaNN     !Sigma value of a neighbor node
!     jgf46.00 Explicitly declared the following variables
    INTEGER :: NN
    INTEGER :: J
    INTEGER :: LBelo
    INTEGER :: LAbov
    INTEGER :: LTry
    INTEGER :: IDiag
    REAL(SZ) SigBelo
    REAL(SZ) SigAbov
    REAL(SZ) SigTry

    call setMessageSource("zsurfbuoy")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif
    IDiag=0

    IF(SigmaNN <= 1.0001*b) THEN !if into ground then skip
        SigBelo=-999
        SigAbov=-999
        BCPressNN=-999.
        GOTO 100
    ENDIF
    IF((SigmaNN > 1.0001*b) .AND. (SigmaNN <= b)) THEN !at bottom then use bottom
        LBelo=1
        BCPressNN=BCP(NN,LBelo)
        SigBelo=b
        SigAbov=b
        GOTO 100
    ENDIF
    IF(SigmaNN >= a) THEN     !into air use surface
        LAbov=NFEN
        BCPressNN=BCP(NN,LAbov)
        SigBelo=a
        SigAbov=a
        GOTO 100
    ENDIF

    LTry=J                    !start search for SIGABOV and SIGBELO
    SigTry=Sigma(LTry)
    IF(SigmaNN > SigTry) THEN !too low
        SigBelo=SigTry         !SIGBELO may = SIGTRY
        LBelo=LTry
        LTry=LTry+1            !look at next level higher
        90 SigTry=Sigma(LTry)
        IF(SigmaNN > SigTry) THEN !still too low
            SigBelo=SigTry
            LBelo=LTry
            LTry=LTry+1
            GOTO 90
        ENDIF
        SigAbov=SigTry         !found upper bracketing sigma
        LAbov=LTry
        GOTO 99                !go interpolate
    ENDIF
    IF(SigmaNN <= SigTry) THEN !to high
        SigAbov=SigTry         !SIGABOV may = SIGTRY
        LAbov=LTry
        LTry=LTry-1            !look at next level lower
        91 SigTry=Sigma(LTry)
        IF(SigmaNN <= SigTry) THEN !still too high
            SigAbov=SigTry
            LAbov=LTry
            LTry=LTry-1
            GOTO 91
        ENDIF
        SigBelo=SigTry         !found lower bracketing sigma
        LBelo=LTry
    ENDIF

    99 BCPressNN=(BCP(NN,LAbov)-BCP(NN,LBelo))  & !interpolation
    *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + BCP(NN,LBelo)

    100 CONTINUE

    IF(IDiag == 2) THEN
        WRITE(2,*) '******** ZSURFBUOY **********'
        WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov', &
        '      BCPressNN'
        WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,BCPressNN
        777 FORMAT(I7,I5,3(F10.3),E14.5)
    ENDIF

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!**********************************************************************
    END SUBROUTINE ZSURFBUOY
!**********************************************************************


!*************************************************************************
!     Subroutine to determine local advection state for each element

!                                   Corbitt 120322
!*************************************************************************

    SUBROUTINE ADVECTLOCAL(IE)

    USE SIZES, ONLY : SZ
    USE GLOBAL
    USE MESH, ONLY : NM, DP
    USE NodalAttributes, ONLY : AdvectionState

    Implicit None

    Integer :: NM1, NM2, NM3,IE

    NM1=NM(IE,1)
    NM2=NM(IE,2)
    NM3=NM(IE,3)

    IF ((DP(NM1) >= AdvectionState(NM1)) .AND. &
    (DP(NM2) >= AdvectionState(NM2)) .AND. &
    (DP(NM3) >= AdvectionState(NM3))) THEN
        IFNLCT = IFNLCTE
        IFNLCAT = IFNLCATE
    ELSE
        IFNLCT = 0
        IFNLCAT = 0
    ENDIF

!***********************************************************************
    END SUBROUTINE ADVECTLOCAL
!***********************************************************************

!*************************************************************************
!   Subroutine to check if elemental slope limiting is needed
!       Chris Massey, USACE-ERDC-CHL Dec. 8, 2014 -- Removed this
!         section from timestep subroutine and put it in its own
!         subroutine.
!*************************************************************************
    subroutine check_slopes(it,TimeLoc)
    USE SIZES, ONLY : SZ, mnproc
    use global, only : eta2,nodecode,NOFF,ESLONOFF,screenUnit, &
    MyProc,setMessageSource, unsetMessageSource, scratchMessage, &
    allMessage,logMessage, DEBUG, ECHO, INFO, WARNING, ERROR, &
    nodes_lg
    use mesh, only : x, y, nm, ne, areas, sfac
    use NodalAttributes, ONLY : LoadEleSlopeLim, &
    elemental_slope_limiter_active, &
    elemental_slope_limiter_grad_max, &
    elemental_slope_limiter_max_exceeded
    implicit none
    INTEGER, intent(in) :: IT
    Real(8), intent(in) :: TimeLoc
    INTEGER :: IE,I
    INTEGER :: NM1, NM2, NM3, NM123
    INTEGER :: NC1, NC2, NC3, NCEle, NCI
    REAL(SZ) DEta2DX,DEta2DY,DEta2Mag,SFacAvg
    REAL(SZ) FDX1,FDX2,FDX3,FDY1,FDY2,FDY3
    REAL(8) :: AreaIE2
    integer :: nodeNumber ! fulldomain node number where wse slope is exceeded
          

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
    elemental_slope_limiter_active(:) = .FALSE. 

! bell  CHECK IF THE LOCAL GRADIENT FOR ELEMENTS IS EXCEEDED AND TRIGGER THE
!...  ELEMENTAL SLOPE LIMITER ACCORDINGLY (original routine by Crystal Fulcher)

    DO IE=1,NE
        NM1=NM(IE,1)
        NM2=NM(IE,2)
        NM3=NM(IE,3)
        NC1=NODECODE(NM1)
        NC2=NODECODE(NM2)
        NC3=NODECODE(NM3)
        NCEle=NC1*NC2*NC3*NOFF(IE)
        IF(NCEle == 0)THEN
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
        dEta2Dx  = (Eta2(NM1)*FDX1+Eta2(NM2)*FDX2+Eta2(NM3)*FDX3) &
        /AreaIE2
        dEta2Dy  = (Eta2(NM1)*FDY1+Eta2(NM2)*FDY2+Eta2(NM3)*FDY3) &
        /AreaIE2
        dEta2Mag = sqrt(dEta2Dx*dEta2Dx + dEta2Dy*dEta2Dy)
    
    ! jgf51.51: Now that the slope limiter gets reset at
    ! every time step, I had to rewrite the logging so that
    ! a log message is only written the first time the slope
    ! limiter is activated at a node during a particular run.
        DO I=1,3
        ! If the limiter is on already, go to the next node.
            IF (elemental_slope_limiter_active(NM(IE,I))) CYCLE
        ! Compare the elemental slope to the maximum elemental gradient.
            IF (dEta2Mag >= &
            ABS(elemental_slope_limiter_grad_max(NM(IE,I)))) THEN
            ! jgf51.51: Log the fulldomain node number.
                nodeNumber = nm(ie,i)
                if (mnproc > 1) then
                    nodeNumber = nodes_lg(nm(ie,i))
                endif
            ! zc - If gradmax is positive or zero, activate slope
            ! limiting.
                if (elemental_slope_limiter_grad_max(nm(ie,i)) &
                 >= 0.0d0) then
                ! If it is the first time that the slope limiter
                ! has been activated at this node, write a log
                ! message.
                    if (eslonoff(nm(ie,i)) == 0) then
                        write(scratchMessage,1983) nodeNumber,dEta2Mag, &
                        elemental_slope_limiter_grad_max(NM(IE,I)), &
                        it, timeLoc
                        call allMessage(INFO,scratchMessage)
                        eslonoff(nm(ie,i)) = 1 ! for output file
                    endif
                    elemental_slope_limiter_active(NM(IE,I)) = .TRUE. 
                else
                ! Just print log message the first time the
                ! gradient is exceeded.
                    IF (elemental_slope_limiter_max_exceeded(NM(IE,I)) &
                    .eqv. .FALSE. ) THEN
                        write(scratchMessage,1984) nodeNumber,dEta2Mag, &
                        elemental_slope_limiter_grad_max(NM(IE,I)), &
                        it, timeLoc
                        call allMessage(INFO,scratchMessage)
                        elemental_slope_limiter_max_exceeded(NM(IE,I)) = &
                         .TRUE. 
                    endif
                endif
            endif
        enddo ! loop around nodes of an element
    enddo ! loop over the elements


    1983 format('Elemental slope limiter turned on at fulldomain node ',i0, &
    ' where the elemental slope is ',1pE12.4E3, &
    ' and the maximum elemental slope is ',1pE12.4E3, &
    ' on time step ',i0,' and time = ',e15.8,'.')

    1984 format('Maximum elemental slope exceeded at fulldomain node ',i0, &
    ' where the elemental slope is ',1pE12.4E3, &
    ' and the maximum elemental slope is ',1pE12.4E3, &
    ' on time step ',i0,' and time = ',e15.8,'.')
         
!...  END CHECKING ELEMENT GRADIENTS

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
!-----------------------------------------------------------------
    END SUBROUTINE CHECK_SLOPES
!-----------------------------------------------------------------

!***********************************************************************

!***********************************************************************
!  Apply Elemental Slope Limiter
!    Chris Massey, USACE-ERDC-CHL, Dec. 8, 2014
!       Made into a subroutine

!***********************************************************************

    SUBROUTINE APPLY_SLOPE_LIMITS(ETA2Lim,LocNP)
    USE SIZES, ONLY : SZ
    use global, only : nodecode,NOFF,IFNLFA, &
    setMessageSource, unsetMessageSource, allMessage, &
    logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
    use mesh, only : ne, nm, areas, totalArea
    use NodalAttributes, ONLY : LoadEleSlopeLim, &
    elemental_slope_limiter_active
    implicit none
    integer :: IE,NM1,NM2,NM3,NC1,NC2,NC3,NCEle
    integer, intent(in) :: LocNP
    REAL(SZ) :: EtaN1,EtaN2,EtaN3,EtaN123
    real(sz), intent(inout) :: Eta2lim(LocNP)
    REAL(8) :: AreaEle
    REAL(sz), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
    LOGICAL, SAVE :: firstCall = .TRUE. 


    call setMessageSource("check_slopes")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    IF (LoadEleSlopeLim.eqv. .TRUE. ) THEN
        IF (firstCall.eqv. .TRUE. ) THEN
            allocate(elevSum(LocNP))
            firstCall = .FALSE. 
            elevSum(:) = 0.d0
        ENDIF
    ENDIF


!       ELEMENTAL SLOPE LIMITER

!        CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
!        TIME STEP AND APPLY SMOOTHING IF DESIRED. -ZC

    IF (LoadEleSlopeLim.eqv. .TRUE. ) THEN
    
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

    ! bell   CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
    !        TIME STEP AND APPLY SMOOTHING IF DESIRED. THIS ROUTINE PARALLELS THE
    !        ABOVE ROUTINE.
        WHERE ((elemental_slope_limiter_active.eqv. .TRUE. ) .AND. &
            (TotalArea /= 0.d0))
            Eta2Lim = elevSum / TotalArea
        END WHERE
    ENDIF

!... Will apply the updating outside the subroutine
! ifdef CMPI
!      CALL UPDATER(ETA2Lim,DUMY1,DUMY2,1)
! endif


#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()


    RETURN

    END SUBROUTINE APPLY_SLOPE_LIMITS

!***********************************************************************

