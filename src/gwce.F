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
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C 
C                       M O D U L E   G W C E
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C
C Generalized Wave Continuity Equation (GWCE) sets up and solves the 
C continuity (conservation of mass) equation. 
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
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CPRECOR          - Use the predictor-corrector algorithm for GWCE       *
C                        and momentum equations (package deal)                *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     CSmag_Eh         - Use Smagorinski, spatially varying, vertically       *
C                        integrated lateral viscosity coefficient             *
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C******************************************************************************

      module gwce
      use mod_logging, only : DEBUG, INFO, WARNING, ERROR, logMessage,
     &   setMessageSource, unsetMessageSource, allMessage
      use mod_terminate, only: terminate
      implicit none

      real(8),allocatable :: lsxx(:),lsxy(:),lsyx(:),lsyy(:)
      real(8),allocatable :: etas0(:)      
      real(8),allocatable :: eta0(:)
      real(8) :: convcr
      integer :: isldia
      integer :: itmax
      real(8) :: a00,b00,c00
      real(8) :: dt2
      real(8) :: ga00 

C     Arrays used by JCG iterative solver
      integer :: nw
      integer, allocatable :: iwksp(:)
      real(8), allocatable :: wksp(:)
      integer, allocatable :: iparm(:)
      real(8), allocatable :: rparm(:)
      integer :: numitr ! number of iterations per timestep 

cjjw...added 2 lines
C.....for matrix conditioning,global to allow for reading in from hotstart file
      REAL(8) :: EP !jgf45.08 EP global for predictor-corrector
      !
      ! for managing datasets for aperiodic elevation specified boundaries (fort.19)
      real(8) :: etime1, etime2, etiminc 

      ! jgf52.30.04: Allow the inverted barometer condition to be 
      ! activated on elevation-specified boundaries so that low pressure
      ! systems crossing the boundary won't trigger numerical instabilities.
      ! This is normally .false. but can be set to .true. in the fort.15
      ! in the metControl namelist. 
      logical :: invertedBarometerOnElevationBoundary = .false.

      private

      public :: initGwce, solveGwce, a00, b00, c00, convcr, isldia, itmax,
     &          invertedBarometerOnElevationBoundary, etime1, etime2, etiminc,
     &          numitr

      contains

!----------------------------------------------------------------------
!           S U B R O U T I N E   I N I T  G  W  C  E
!----------------------------------------------------------------------
!     Allocates arrays and initializes constants.
!----------------------------------------------------------------------
      subroutine initGWCE()
      use sizes, only : mnp
      use global, only : dt, ilump, alloc_main11, alloc_main11_lumped
      USE ADC_CONSTANTS, only: g
      use mesh, only : np
      use itpackv, only : dfault
      implicit none

      allocate(lsxx(mnp),lsxy(mnp),lsyx(mnp),lsyy(mnp))
      allocate(etas0(mnp))
      allocate(eta0(mnp))
      allocate( iparm(12),rparm(12) )
      allocate( iwksp(3*mnp),wksp(4*mnp+400) )
      dt2=dt*2.d0
      ga00=g*a00

C     jgf48.4619: Accommodate Seizos changes for explicit solve
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
C        allocate arrays needed by GWCE matrix and iterative solver
         call alloc_main11()
C        initialize parameter arrays needed by iterative solver
         CALL DFAULT(IPARM,RPARM)
         IPARM(1)=ITMAX
         IPARM(2)=ISLDIA
         !cms51.06: moved opening of fort.33 to openLogFile sub in global.F
         IPARM(4)=33
         RPARM(1)=CONVCR
         NW = 4*NP + 4*ITMAX
      ELSE ! lumped LHS
         call alloc_main11_lumped()
      ENDIF

!----------------------------------------------------------------------
      end subroutine initGWCE
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     S U B R O U T I N E   S O L V E   G   W   C   E
!----------------------------------------------------------------------
C     Two mutually exclusive possibilities are (1)
C     the conventional GWCE or (2) the predictor-corrector GWCE.
!----------------------------------------------------------------------
      subroutine solveGWCE(it, ITIME_BGN, timeloc, timeh)
      use sizes, only : mnp
      use global, only : CGWCE_New, CPRECOR, ETA2, QX2, QY2, UU2, VV2
#ifdef CMPI
     & , dumy1, dumy2
      use messenger
#endif
      use momentum, only : momentum_eq_nonconservative
      use nodalattributes, only : loadEleSlopeLim
      implicit none
      integer, intent(in) :: it, ITIME_BGN ! time step 
      real(8), intent(in) :: timeloc ! simulation time in seconds
      real(8), intent(in) :: timeh

C       -------------PREDICTOR-CORRECTOR-GWCE---------------------------
ckmd    Changed to include the predictor-corrector algorithm
      IF(CPRECOR) THEN
         CALL GWCE_New(IT,ITIME_BGN,TimeLoc,TimeH)
c...tcm added call for slope limiting         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
         CALL momentum_eq_nonconservative()
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
C       -------------PREDICTOR-CORRECTOR-GWCE---------------------------

C                         -OR-

C      -------------------GWCE------------------------------------------
      IF(CGWCE_New) THEN
         CALL GWCE_New(IT,ITIME_BGN,TimeLoc,TimeH)
c.... tcm added call to slope limiter         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
      ENDIF
C     -------------------GWCE------------------------------------------
!----------------------------------------------------------------------
      end subroutine solveGWCE
!----------------------------------------------------------------------


C******************************************************************************
C                                                                             *
C    Subroutine to compute the elevation using the GWCE formluation           *
C    Re-written to conform to the ADCIRC Theory Report                        *
C                                                                             *
C                            r.l.  06/22/2005                                 *
C******************************************************************************
      SUBROUTINE GWCE_New(IT,ITIME_BGN,TimeLoc,TimeH)

#ifdef IEEE_DEBUG
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : GWCE_LV, COEF, COEFD, OBCCOEF, ETAS, ETA1, ETA2,
     &   ESBIN1, ESBIN2, C2DDI, C3D, CBaroclinic, CGWCE_Advec_C1, H1, H2,
     &   CGWCE_Advec_C2, CGWCE_Advec_NC, CGWCE_LS_2PartQ, CGWCE_LS_2PartSQ,
     &   CGWCE_LS_2PartV, CGWCE_LS_2PartSV, CGWCE_LS_KGQ, CPRECOR,
     &   CSmag_Eh, Smag_Comp_Flag, CTIP, DT, FluxSettlingIT, IFNLCAT,
     &   IFNLCT, IFNLFA, NBFR, NODECODE, NRS, NSCREEN, NCCHANGE, IDEN,
     &   NWS, RAMPELEV, RAMPMETE, RES_BC_FLAG,
     &   SMAG_LOWER_LIM, SMAG_UPPER_LIM, TK, NOFF, UU1, VV1, QX1, QY1,
     &   PR1, TiP1, BSX1, BSY1, DUU1, DUV1, DVV1, WSX1, WSY1, VIDBCPDXOH,
     &   VIDBCPDYOH, CORIF, PER, FACE, FF, EFA, EMO, AMIG, PR2, LNM_BC,
     &   QN1, EN1, QN0, ElevDisc, QN2, EN0, iLump, BCFLAG_LNM, EN2,
     &   TKM, NPERSEG, NTIP, NNPERBC, IPERCONN, VIDispDXOH,
     &   VIDispDYOH, IFSFM, IFNL_HDP, CGWCE_HDP, CAliDisp, H0, ALPHAL, windlim,
     &   usingDynamicWaterLevelCorrection, dynamicWaterLevelCorrection1,
     &   dynamicWaterLevelCorrectionDelta, 
     &   COEFDTemp, COEFDTempMem  ! Added for VED1D  08-11-2022 SB
#ifdef CMPI
     &    , dumy1, dumy2, rnp_global
#endif
      use mod_logging, only: screenUnit
      USE ADC_CONSTANTS, ONLY: G, RHOWAT0, Bd, Ad, Cs2
      USE MESH, ONLY : NE, NP, NM, DP, NNeigh, NeiTab, TotalArea, FDXE,
     &   Areas, NEIMAX, SFAC, nneighele, neitabele, FDYE,
     &   SWITCH_ELTAB_PERBC, SFacEle, SFMYEle, SFMXEle, SFCXEle,
     &   SFCYEle, YCSFacEle, TANPHI, SFCT, TANPHIEle, LBArray_Pointer,
     &   X, Y
      USE BOUNDARIES, ONLY : NETA, NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, 
     &   NFLUXRBC, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3, NPEBC,
     &   NBOU, NVELL, IBCONN, ISSUBMERGED64, NFLUXIB64_GBL  ! Added for VEW1D  08-11-2022 SB
      USE ITPACKV
      USE NodalAttributes, ONLY :
     &     LoadGeoidOffset, GeoidOffset, EVM, NOLIBF,
     &     TAU0VAR, HighResTimeVaryingTau0, FullDomainTimeVaryingTau0,
     &     CalculateTimeVaryingTau0, LoadAdvectionState, advectlocal,
     &     LoadCondensedNodes, ListCondensedNodes, NListCondensedNodes, ! Added for VEW1D  08-11-2022 SB
     &     NNodesListCondensedNodes, NCondensedNodes,                    !
C... DW
     &     LoadAbsLayerSigma, absorblayer_sigma_eta,
     &     absorblayer_sigma_mnx, absorblayer_sigma_mny
      USE SPONGELAYER
C... DW
C... SB
      USE VEW1D, ONLY : ROTATE_AT_CONDENSEDNODES_ALL
C... SB
#ifdef CMPI
      USE MESSENGER
#endif
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceECB, 
     &  checkChange, enforceEOB, enforceEIB, enforceGWCELVOB
      USE WIND, ONLY: PRBCKGRND_MH2O, windLimiter

!JLW: add subgrid additions
      USE subgrid, ONLY: subgrid_level0 => level0, subgrid_level1 => level1
     &   , subgridVertList
     &   , gridDepthVertETA1, gridDepthVertETA2
     &   , wetFracVertETA1, wetFracVertETA2
     &   , cadvVertETA2, dphidt
      use mod_gwce_bc_forcing, only: apply_boundary_conditions, apply_vew1d_and_condensed_nodes,
     &                           apply_subdomain_boundary_nodes, update_coef_periodic_sponge_layer_lumped,
     &                           update_coef_periodic_sponge_layer_consistent, update_periodic_sponge_layer_nodes

      IMPLICIT NONE

      integer,intent(in) :: IT
      integer,intent(in) :: ITIME_BGN
      real(8),intent(in) :: TimeLoc
      real(8),intent(in) :: TimeH

      INTEGER IE, JN, IJ, I, J, K, L                     !local loop counters
      INTEGER NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      LOGICAL  DIE
      REAL(8) A00pB00
      REAL(8) BCXAvg, BCYAvg
      REAL(8) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(8) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(8) CorifAvg
      REAL(8) DPAvg, GDPAvgOAreaIE4
      REAL(8) DPAvgS0, GDPAvgOAreaIE4_S0
      REAL(8) DispX, DispY, DispXAvg, DispYAvg
      REAL(8) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(8) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(8) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(8) ESN1, ESN2, ESN3, ESAvg
      REAL(8) EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(8) EVMEle, EVMSmag
      REAL(8) GA00DPAvgOAreaIE4
      REAL(8) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(8) H1N1, H1N2, H1N3, HAvg, H2N
      REAL(8) H2OTotalArea
      REAL(8) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(8) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(8) MsFacR, MsFacLOnDiag, MsFacLOffDiag
      REAL(8) MX, MY, MXAvg, MYAvg
      REAL(8) JXAvg, JYAvg
      REAL(8) Pr1N1, Pr1N2, Pr1N3
      REAL(8) Off1N1, Off1N2, Off1N3 ! water level offsets
      REAL(8) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(8) SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac 
      REAL(8) SFCXAvg, SFCYAvg, YCFACAvg, SFCXADJ, SFCYADJ, TANPHIAvg
      REAL(8) T0N1,T0N2, T0N3
      REAL(8) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(8) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(8) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(8) TiPN1, TiPN2, TiPN3
      REAL(8) U1N1,U1N2,U1N3, U1Avg
      REAL(8) V1N1,V1N2,V1N3, V1Avg
      REAL(8) WSXAvg, WSYAvg
      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) FDXX1, FDXX2, FDXX3, FDYY1, FDYY2, FDYY3
      REAL(8) HH1 !jgf46.02 Added for Katrina.
      REAL(8) RDIAG ! jgf48.4619 Seizo parameter for fully explicit mode
      INTEGER :: j12, j13, j21, j23, j31, j32
      REAL(8) :: arg, argj
      REAL(8) :: celerity
      REAL(8) :: etratio
      INTEGER :: ier
      INTEGER :: nbdj
      REAL(8) :: rff
      REAL(8) :: qforcei, qforcej
C...... SB
      REAL(8) :: SX, SY
      INTEGER :: ROT_STATUS
C...... SB
C...... DW
      REAL(8):: SigmaAvg, Tau0Avg_S
      REAL(8):: SgN1_eta, SgN2_eta, SgN3_eta
      REAL(8):: SgN1_mnx, SgN2_mnx, SgN3_mnx
      REAL(8):: SgN1_mny, SgN2_mny, SgN3_mny
      REAL(8):: Tau0SigmaAvg, SigmaQxAvg, SigmaQyAvg
      REAL(8):: AvgSigmaHAbsU, AvgSigmaHAbsV
      REAL(8):: MsFacR_S, MsFacRS_T, MsFacRS_ABC00
      REAL(8):: E0AN1, E0AN2, E0AN3
      REAL(8):: E1AN1, E1AN2, E1AN3
      REAL(8):: E2AN1, E2AN2, E2AN3
      REAL(8):: ESAN1, ESAN2, ESAN3 
      REAL(8):: U1AN1, U1AN2, U1AN3
      REAL(8):: V1AN1, V1AN2, V1AN3
      REAL(8):: ABC00_EAN1, ABC00_EAN2, ABC00_EAN3
      REAL(8):: fwsx1, fwsx2, fwsx3
      REAL(8):: fwsy1, fwsy2, fwsy3
      REAL(8) :: SPM1, SPM2, SPM3 
      LOGICAL :: ETA_SPONGE = .false., GWCE_SPONGE = .false.
      REAL(8) :: Ma2, Cfac = 1.0D0, CfacS0 = 1.0D0
      INTEGER :: NNBB1, NNBB2  ! Added for VEW1D  08-11-2022 SB
      REAL(8):: fBuf          !
!JLW: adding total water depth variable for previous TS
      REAL(8) :: H0N1, H0N2, H0N3
      REAL(8) :: CADV1, CADV2, CADV3
      REAL(8) :: PHIAVG1, PHIAVG2
      REAL(8) :: PHI1N1, PHI1N2, PHI1N3
      REAL(8) :: PHI2N1, PHI2N2, PHI2N3

      IF ( LoadAbsLayerSigma ) THEN
         IF (NumNodesAbsLayer(1) > 0 ) ETA_SPONGE = .true.
         IF (SUM(NumNodesAbsLayer) > 0 ) GWCE_SPONGE = .true.
      ENDIF
C..... DW

      call setMessageSource("gwce_new")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C....... DW
C     Switch to the element table to the primary node
      CALL SWITCH_ELTAB_PERBC( IDIREC = 1 ) ;
C.....

      if (subdomainOn.and.enforceBN.eq.2) call checkChange()    ! NCSU Subdomain

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
C..... DW, Dev 2019
C.....  if CGWCE_HDP = .true. use H instead of DP in the surface gradient term
      NCIFBLCK: IF(NCChange.GT.0 .or. CGWCE_HDP ) THEN !if any subdomain grid has changed
         NCChange=0
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
            DO I=1,NP
C              WJP 02.24.2018 Get the infinity norm of the matrix
               TK(I) = max(abs(TKM(1,I))+abs(TKM(3,I)),
     &                     abs(TKM(2,I))+abs(TKM(3,I)))
            ENDDO
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         IF(C3D.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         !jgf48.4619: Added Seizo's handling of Lumping vs non-Lumping
         COEFASSYM : IF ( ILump.eq.0 ) THEN ! default case: fully consistent LHS
            FULLMASS: DO IE=1,NE
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               !JLW: add subgrid additions
               IF(subgrid_level0)THEN
                  H1N1 = gridDepthVertETA2(NMI1)
                  H1N2 = gridDepthVertETA2(NMI2)
                  H1N3 = gridDepthVertETA2(NMI3)
                  PHI1N1 = wetFracVertETA1(NMI1)
                  PHI1N2 = wetFracVertETA1(NMI2)
                  PHI1N3 = wetFracVertETA1(NMI3)
                  PHI2N1 = wetFracVertETA2(NMI1)
                  PHI2N2 = wetFracVertETA2(NMI2)
                  PHI2N3 = wetFracVertETA2(NMI3)
                  PHIAVG1 = (PHI1N1+PHI1N2+PHI1N3)/3.d0
                  PHIAVG2 = (PHI2N1+PHI2N2+PHI2N3)/3.d0
               ENDIF
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg = SFacEle(IE)
C.... Beg DW/WJP:
               SFmxAvg = SFMXEle(IE) ;
               SFmyAvg = SFMYEle(IE) ;
               sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
               sfdyfac = (1 - IFSFM)*1.0D0  + IFSFM*SFmyAvg ; 
               SFCXAvg = SFCXEle(IE) ;
               SFCYAvg = SFCYEle(IE) ;
               YCFACAvg = YCSFACEle(IE) ;
               SFCXADJ = SFCXAvg/sfdxfac ; 
               SFCYADJ = YCFACAvg*SFCYAvg/sfdyfac ; 
               FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1 = FDYE(1,IE)*sfdyfac ; !c  FDY1=X(NM3)-X(NM2) !a1
               FDY2 = FDYE(2,IE)*sfdyfac ; !c  FDY2=X(NM1)-X(NM3) !a2
               FDY3 = FDYE(3,IE)*sfdyfac ; !c  FDY3=X(NM2)-X(NM1) !a3
C.... End WJP/DW
               AreaIE2=Areas(IE)
               AreaIE =AreaIE2*0.5d0
               AreaIE4=AreaIE2*2.0d0
!JLW: adding subgrid to water surface elevation gradiant term
               IF(subgrid_level0)THEN
                  DPAvg = (H1N1+H1N2+H1N3)/3.d0
               ELSE
C.... DW, Nov 2019 if IFNL_HDP = 1 use DP + ETA = H in the GWCE
                  DPAvg = (DP(NMI1)+DP(NMI2)+DP(NMI3))/3.d0 +   ! BAVG if GWCE_NL_HDP = .FALSE.
     &              IFNLFA*IFNL_HDP*(ETA2(NMI1) + ! HAVG otherwise
     &              ETA2(NMI2)+ETA2(NMI3))/3.d0 ; 
               ENDIF   
C.....         WJP adjust the dispersion equation through
C              Ali et al. formula
               if (CAliDisp) then
                  Ma2   = G*DPavg/Cs2
                  Cfac = (1.0D0 - 0.25D0*Ma2 - Ad*DPavg**Bd)**2
               endif
               ! DMW 2022/06 Apply slope limiting to gravity
               GA00DPAvgOAreaIE4 = ALPHAL(IE)*A00*G*DPAvg*Cfac/AreaIE4
               Tau0Avg = (Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0

C........... DW, Absorbing layer
C    
               SigmaAvg = 0.D0 ; Tau0SigmaAvg = 0.D0
               Tau0Avg_S = Tau0Avg
               IF ( ETA_SPONGE ) THEN
                  ! \overline( \sigma ) 
                  SigmaAvg = (absorblayer_sigma_eta(NMI1,1) 
     &                     +  absorblayer_sigma_eta(NMI2,1) 
     &                     +  absorblayer_sigma_eta(NMI3,1))/3.D0

                  ! \overline(\tau_{0} \sigma)
                  Tau0SigmaAvg = (absorblayer_sigma_eta(NMI1,1)*Tau0Var(NMI1)
     &              +  absorblayer_sigma_eta(NMI2,1)*Tau0Var(NMI2)
     &              +  absorblayer_sigma_eta(NMI3,1)*Tau0Var(NMI3))/3.D0

                  Tau0Avg_S = Tau0Avg + 
     &             (A00*2.D0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
               END IF              
!
               MsFacLOnDiag =OnDiag*AreaIE*(1.d0/DT+Tau0Avg_S/2.d0)/DT/12.d0
               MsFacLOffDiag=OffDiag*AreaIE*(1.d0/DT+Tau0Avg_S/2.d0)/DT/12.d0
!JLW: adding subgrid to the on and off diagonal terms
               IF(subgrid_level0)THEN
                  MSFacLOnDiag = MSFacLOnDiag*PHIAVG2 + OnDiag
     &               *dphidt*AreaIE*((PHIAVG2-PHIAVG1)/2.d0)/DT/12.d0
                  MSFacLOffDiag = MSFacLOffDiag*PHIAVG2 + OffDiag
     &               *dphidt*AreaIE*((PHIAVG2-PHIAVG1)/2.d0)/DT/12.d0
               ENDIF
C
C.......... DW
               DO JN=2,NEIMAX
                  IF(NeiTab(NMI1,JN).EQ.NMJ2) J12=JN
                  IF(NeiTab(NMI1,JN).EQ.NMJ3) J13=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ1) J21=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ3) J23=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ1) J31=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ2) J32=JN
               END DO
C.........     WJP multiplying by cosine(lat) if doing the 
C              spherical correction 
C              Node-by-Node
               SPM1 = SFCT(NMI1) 
               SPM2 = SFCT(NMI2)
               SPM3 = SFCT(NMI3)
C              Assemble the LHS matrix
               Coef(NMI1,1)  = Coef(NMI1,1)   + (MsFacLOnDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX1 +
     &                                        SFCYADJ*FDY1*FDY1) )*NCELE
               Coef(NMI1,J12)= Coef(NMI1,J12) + (MsFacLOffDiag*SPM2
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX2 + 
     &                                        SFCYADJ*FDY1*FDY2) )*NCELE
               Coef(NMI1,J13)= Coef(NMI1,J13) + (MsFacLOffDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX3 +
     &                                        SFCYADJ*FDY1*FDY3) )*NCELE

               Coef(NMI2,J21)= Coef(NMI2,J21) + (MsFacLOffDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX1 + 
     &                                        SFCYADJ*FDY2*FDY1) )*NCELE
               Coef(NMI2,1)  = Coef(NMI2,1)   + (MsFacLOnDiag*SPM2
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX2 + 
     &                                        SFCYADJ*FDY2*FDY2) )*NCELE
               Coef(NMI2,J23)= Coef(NMI2,J23) + (MsFacLOffDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX3 +
     &                                        SFCYADJ*FDY2*FDY3) )*NCELE

               Coef(NMI3,J31)= Coef(NMI3,J31) + (MsFacLOffDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX1 + 
     &                                        SFCYADJ*FDY3*FDY1) )*NCELE
               Coef(NMI3,J32)= Coef(NMI3,J32) + (MsFacLOffDiag*SPM2
     &                   + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX2 +
     &                                        SFCYADJ*FDY3*FDY2) )*NCELE
               Coef(NMI3,1)  = Coef(NMI3,1)   + (MsFacLOnDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX3 + 
     &                                        SFCYADJ*FDY3*FDY3) )*NCELE
            ENDDO FULLMASS
         ELSE
            !jgf48.4619: Add Seizo's construction of Lumped LHS matrix
            LUMPEDMASS: DO IE=1,NE ! Make LHS Lumped Matrix
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               AreaIE2=Areas(IE)
               AreaIE =AreaIE2*0.5d0
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0
C
C.... DW, absorbing layer
               SigmaAvg = 0.D0 ; Tau0sigmaAvg = 0.D0
               Tau0Avg_S = Tau0Avg
               IF ( ETA_SPONGE ) THEN 
                 SigmaAvg = (absorblayer_sigma_eta(NMI1,1) 
     &             +  absorblayer_sigma_eta(NMI2,1) 
     &             +  absorblayer_sigma_eta(NMI3,1))/3.D0
                 Tau0SigmaAvg = (absorblayer_sigma_eta(NMI1,1)*Tau0Var(NMI1)
     &             +  absorblayer_sigma_eta(NMI2,1)*Tau0Var(NMI2)
     &             +  absorblayer_sigma_eta(NMI3,1)*Tau0Var(NMI3))/3.D0

                 Tau0Avg_S = Tau0Avg + 
     &             (A00*2.D0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
               END IF
               MsFacLOnDiag =OnDiag*AreaIE*(1.d0/DT+Tau0Avg_S/2.d0)/DT/12.d0
C.............DW
C.........     WJP multiplying by cosine(lat) if doing the 
C              spherical correction 
C              Node-by-Node
               SPM1 = SFCT(NMI1) 
               SPM2 = SFCT(NMI2)
               SPM3 = SFCT(NMI3)

               Coefd(NMI1) = Coefd(NMI1) + MsFacLOnDiag*SPM1*NCELE
               Coefd(NMI2) = Coefd(NMI2) + MsFacLOnDiag*SPM2*NCELE
               Coefd(NMI3) = Coefd(NMI3) + MsFacLOnDiag*SPM3*NCELE
            ENDDO LUMPEDMASS
         ENDIF COEFASSYM 

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
            OBCCoef = 0.d0 !jgf53.dev: later code assumes this is initialized to zero
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

            CALL update_coef_periodic_sponge_layer_consistent( NPERSEG, NNPERBC, IPERCONN, NNEIGH, EP, COEF)

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
            IF (DIE) THEN
               ! jgfdebug
               open(899,file='debug.txt',status='replace',action='write')
                  do i=1,np
                     if (coef(i,1).lt.0.d0) then
                        do j=2,nneighele(i)
                           if (neitabele(i,j).ne.0) then
                              write(6,
     & '("Node ",i0," element ",i0," area=",f15.7)') i, j, 
     &  areas(neitabele(i,j))
                           endif
                        enddo
                     endif
                  enddo
               close(899)
               ! end jgfdebug
               CALL terminate()
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
C......  DW
            CALL update_coef_periodic_sponge_layer_lumped(NPERSEG, NNPERBC, IPERCONN, EP, COEFD)
C......  DW
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
            IF (DIE) CALL terminate()
         ENDIF
      ENDIF NCIFBLCK                    !End of GWCE matrix setup
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
      IF((NWS==0).AND.(NRS==0)) THEN
        WSXAvg=0.d0
        WSYAvg=0.d0
        Pr1N1=0.d0
        Pr1N2=0.d0
        Pr1N3=0.d0
      ENDIF

      ! jgf: Adding support for water level offset as a pseudo
      ! barometric pressure
      if (usingDynamicWaterLevelCorrection.eqv..false.) then
         Off1N1 = 0.d0
         Off1N2 = 0.d0
         Off1N3 = 0.d0
      endif

      IF (.not.CTIP) THEN
        TiPN1=0.d0
        TiPN2=0.d0
        TiPN3=0.d0
      ENDIF

      IF(.not.C3D) THEN
        DispXAvg=0.d0
        DispYAvg=0.d0
      ENDIF

      IF(.not.CBaroclinic) THEN
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
            SFacAvg = SFacEle(IE)
            SFmxAvg=SFMXEle(IE) ;
            SFmyAvg=SFMYEle(IE) ;
            sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
            sfdyfac=(1 - IFSFM)*1.0D0  + IFSFM*SFmyAvg ;            
C........ BEG DW:
            FDX1 = FDXE(1,IE)*sfdxfac ; 
            FDX2 = FDXE(2,IE)*sfdxfac ; 
            FDX3 = FDXE(3,IE)*sfdxfac ; 
            FDY1 = FDYE(1,IE)*sfdyfac ; 
            FDY2 = FDYE(2,IE)*sfdyfac ; 
            FDY3 = FDYE(3,IE)*sfdyfac ; 
C........ END DW:
            LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0   !A*DUDX
            LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0   !A*DUDY
            LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0   !A*DVDX
            LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0   !A*DVDY
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            !
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
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


!JLW: adding subgrid to velocity based lateral stress
               IF(subgrid_level0)THEN
                  H2N = gridDepthVertETA2(I)
               ELSE
C..............DMW202401 Use saved H2
                  H2N= H2(I)
               ENDIF
               H2OTotalArea=H2N/TotalArea(I)

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
!JLW: add subgrid wet/dry to lateral stress term so we can use Smag
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg = SFacEle(IE)
C..... Beg DW:
            SFmxAvg=SFMXEle(IE) ;
            SFmyAvg=SFMYEle(IE) ;
            sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
            sfdyfac=(1 - IFSFM)*1.0D0  + IFSFM*SFmyAvg ;            

            FDX1 = FDXE(1,IE)*sfdxfac ; 
            FDX2 = FDXE(2,IE)*sfdxfac ; 
            FDX3 = FDXE(3,IE)*sfdxfac ; 
            FDY1 = FDYE(1,IE)*sfdyfac ; 
            FDY2 = FDYE(2,IE)*sfdyfac ; 
            FDY3 = FDYE(3,IE)*sfdyfac ; 
C..... End DW
            EVMEle = NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.d0
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
               LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
               LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
               LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
               LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0
               EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               ! tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
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
      
#ifdef CMPI
      IF ((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartSQ).OR.(CGWCE_LS_2PartV).OR.(CGWCE_LS_2PartSV)) THEN
         CALL UPDATER(LSXX,LSXY,DUMY1,2)
         CALL UPDATER(LSYX,LSYY,DUMY1,2)
      ENDIF
#endif

C...  Assemble the GWCE RHS except for the boundary integral terms

      DO IE=1,NE

C...     Set nodal values for each element

Corbitt 120322: Localized Advection
        IF (LoadAdvectionState) CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
!JLW: moving NCEle calulation further down. 
!DW: uncomment line NCELE =  below to make the subgrid code run  
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

!JLW: similar to the LHS, look up wet area fraction, grid averaged water
!depth, and advection
        IF(subgrid_level0)THEN
            H0N1 = gridDepthVertETA1(NM1)
            H0N2 = gridDepthVertETA1(NM2)
            H0N3 = gridDepthVertETA1(NM3)
            H1N1 = gridDepthVertETA2(NM1)
            H1N2 = gridDepthVertETA2(NM2)
            H1N3 = gridDepthVertETA2(NM3)
            IF(subgrid_level1)THEN
               CADV1 = cadvVertETA2(NM1)
               CADV2 = cadvVertETA2(NM2)
               CADV3 = cadvVertETA2(NM3)
            ENDIF
            PHI1N1 = wetFracVertETA1(NM1)
            PHI1N2 = wetFracVertETA1(NM2)
            PHI1N3 = wetFracVertETA1(NM3)
            PHI2N1 = wetFracVertETA2(NM1)
            PHI2N2 = wetFracVertETA2(NM2)
            PHI2N3 = wetFracVertETA2(NM3)
            PHIAVG1 = (PHI1N1+PHI1N2+PHI1N3)/3.d0
            PHIAVG2 = (PHI2N1+PHI2N2+PHI2N3)/3.d0
         ELSE

C...........DMW202401 Use saved H2, H1
            H0N1=H1(NM1)
            H0N2=H2(NM1)
            H0N3=H1(NM3) 

C...........DMW202401 Use saved H2
            H1N1=H2(NM1)
            H1N2=H2(NM2)
            H1N3=H2(NM3)
            NCELE=NC1*NC2*NC3*NOFF(IE)
         ENDIF
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         T0N1=Tau0Var(NM1)
         T0N2=Tau0Var(NM2)
         T0N3=Tau0Var(NM3)
C..... SB 2022-09-14
C        Rorate the velocity and flux vectors at the condensed nodes
C        so that they align with the stream direction
         IF (LoadCondensedNodes) THEN
            CALL ROTATE_AT_CONDENSEDNODES_ALL
     &        (NM1,NM2,NM3,
     &         U1N1,V1N1,U1N2,V1N2,U1N3,V1N3,
     &         QX1N1,QY1N1,QX1N2,QY1N2,QX1N3,QY1N3,
     &         ROT_STATUS)
         ENDIF
C..... SB
C
C.... DW: begin ! absorbing layer
C 
         E0AN1 = 0.D0 ; E0AN2 = 0.D0 ; E0AN3 = 0.D0 ;
         E1AN1 = 0.D0 ; E1AN2 = 0.D0 ; E1AN3 = 0.D0 ;
         E2AN1 = 0.D0 ; E2AN2 = 0.D0 ; E2AN3 = 0.D0 ;
         ESAN1 = 0.D0 ; ESAN2 = 0.D0 ; ESAN3 = 0.D0 ;
 
         ABC00_EAN1 = 0.D0 ; 
         ABC00_EAN2 = 0.D0 ;
         ABC00_EAN3 = 0.D0 ;

         U1AN1 = 0.D0 ; U1AN2 = 0.D0 ; U1AN3 = 0.D0 ;
         V1AN1 = 0.D0 ; V1AN2 = 0.D0 ; V1AN3 = 0.D0 ;

         SgN1_eta = 0.D0 ; SgN2_eta = 0.D0 ; SgN3_eta = 0.D0 ;
         SgN1_mnx = 0.D0 ; SgN2_mnx = 0.D0 ; SgN3_mnx = 0.D0 ;
         SgN1_mny = 0.D0 ; SgN2_mny = 0.D0 ; SgN3_mny = 0.D0 ;
  
         IF ( GWCE_SPONGE ) THEN 
            SgN1_eta = absorblayer_sigma_eta(NM1,1) ; 
            SgN2_eta = absorblayer_sigma_eta(NM2,1) ; 
            SgN3_eta = absorblayer_sigma_eta(NM3,1) ;

            SgN1_mnx = absorblayer_sigma_mnx(NM1,1) ; 
            SgN2_mnx = absorblayer_sigma_mnx(NM2,1) ; 
            SgN3_mnx = absorblayer_sigma_mnx(NM3,1) ;

            SgN1_mny = absorblayer_sigma_mny(NM1,1) ; 
            SgN2_mny = absorblayer_sigma_mny(NM2,1) ; 
            SgN3_mny = absorblayer_sigma_mny(NM3,1) ;

            E0AN1 = eta0_abslayer(NM1) ;
            E0AN2 = eta0_abslayer(NM2) ; 
            E0AN3 = eta0_abslayer(NM3) ; 

            E1AN1 = eta1_abslayer(NM1) ; 
            E1AN2 = eta1_abslayer(NM2) ; 
            E1AN3 = eta1_abslayer(NM3) ; 

            E2AN1 = eta2_abslayer(NM1) ;
            E2AN2 = eta2_abslayer(NM2) ; 
            E2AN3 = eta2_abslayer(NM3) ;

            ESAN1 = E2AN1 - E0AN1 ;
            ESAN2 = E2AN2 - E0AN2 ;
            ESAN3 = E2AN3 - E0AN3 ;

            ABC00_EAN1 = A00*E2AN1 + B00*E1AN1 + C00*E0AN1 ; 
            ABC00_EAN2 = A00*E2AN2 + B00*E1AN2 + C00*E0AN2 ; 
            ABC00_EAN3 = A00*E2AN3 + B00*E1AN3 + C00*E0AN3 ; 

            U1AN1 = uu1_abslayer(NM1) ; 
            U1AN2 = uu1_abslayer(NM2) ; 
            U1AN3 = uu1_abslayer(NM3) ; 
 
            V1AN1 = vv1_abslayer(NM1) ;
            V1AN2 = vv1_abslayer(NM2) ;
            V1AN3 = vv1_abslayer(NM3) ;          
         END IF
C
         IF ((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            Pr1N1=PR1(NM1)
            Pr1N2=PR1(NM2)
            Pr1N3=PR1(NM3)

            !c Pressure c!
            IF ( NO_MET_IN_SPONGE ) THEN
               IF ( absorblayer_sigma_eta(NM1,1) > 1.0e-9 ) THEN
                  Pr1N1 = PRBCKGRND_MH2O
               END IF
               IF ( absorblayer_sigma_eta(NM2,1) > 1.0e-9 ) THEN
                  Pr1N2 = PRBCKGRND_MH2O
               END IF
               IF ( absorblayer_sigma_eta(NM3,1) > 1.0e-9 ) THEN
                  Pr1N3 = PRBCKGRND_MH2O
               END IF
            END IF
         ENDIF
C ....,... END DW
C
         ! jgf: Adding support for water level offset as a pseudo barometric
         ! pressure
         if (usingDynamicWaterLevelCorrection.eqv..true.) then
            Off1N1 = dynamicWaterLevelCorrection1(NM1)
            Off1N2 = dynamicWaterLevelCorrection1(NM2)
            Off1N3 = dynamicWaterLevelCorrection1(NM3)
         endif
         IF (CTIP) THEN                        !tidal potential
            TiPN1=TiP1(NM1)
            TiPN2=TiP1(NM2)
            TiPN3=TiP1(NM3)
         ENDIF
         IF (C2DDI) THEN                       !2D bottom friction
            !WJP taking into account the symmetrical friction matrix
            ! CPB: added resynthesized tidal velocity option
            BSXN1 = TKM(1,NM1)*QX1N1 + TKM(3,NM1)*QY1N1
            BSYN1 = TKM(3,NM1)*QX1N1 + TKM(2,NM1)*QY1N1
            BSXN2 = TKM(1,NM2)*QX1N2 + TKM(3,NM2)*QY1N2
            BSYN2 = TKM(3,NM2)*QX1N2 + TKM(2,NM2)*QY1N2
            BSXN3 = TKM(1,NM3)*QX1N3 + TKM(3,NM3)*QY1N3
            BSYN3 = TKM(3,NM3)*QX1N3 + TKM(2,NM3)*QY1N3
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
         SFacAvg= SFacEle(IE)
C..... DW
         ! derivative: accomodate corrected spherical coordinate
         SFmxAvg=SFMXEle(IE)  ;
         SFmyAvg=SFMYEle(IE) ;
         sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
         sfdyfac=(1 - IFSFM)*1.0D0  + IFSFM*SFmyAvg ; 
         ! continuity
         SFCXAvg=SFCXEle(IE) ;
         SFCYAvg=SFCYEle(IE) ;
         YCFACAvg=YCSFACEle(IE) ;
         TANPHIAvg=TANPHIEle(IE) ; 
         SFCXADJ = SFCXAvg/sfdxfac ; 
         SFCYADJ = YCFACAvg*SFCYAvg/sfdyfac ; 

         FDX1 = FDXE(1,IE)*sfdxfac ; !c dx/dlambda*b1*2A
         FDX2 = FDXE(2,IE)*sfdxfac ; !c dx/dlambda*b2*2A
         FDX3 = FDXE(3,IE)*sfdxfac ; !c dx/dlambda*b3*2A
         FDY1 = FDYE(1,IE)*sfdyfac ; !c dy/dphi*a1*2A
         FDY2 = FDYE(2,IE)*sfdyfac ; !c dy/dphi*a1*2A
         FDY3 = FDYE(3,IE)*sfdyfac ; !c dy/dphi*a1*2A
C..... DW

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

C...... BEG: DW, Nov 2019
         Tau0XGrad2A=SFCXADJ*Tau0XGrad2A      !2*Area*dTau0/dx
         Tau0YGrad2A=SFCYADJ*Tau0YGrad2A      !2*Area*dTau0/dy
C...... END DW

C...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM
C...     ********     Need to think about this more: DW ********
C...     ********     Kolar & Gray does not seem to be suitable for  
C...     ********     the spherical coordinates        
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
C           X-direction
            LSXYEle = (LSXY(NM1)+LSXY(NM2)+LSXY(NM3))/3.d0
            LSYXEle = (LSYX(NM1)+LSYX(NM2)+LSYX(NM3))/3.d0
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
     &         -TANPHIAvg*(LSXYEle + LSYXEle)
C           Y-direction 
            LSXXEle = (LSXX(NM1)+LSXX(NM2)+LSXX(NM3))/3.d0
            LSYYEle = (LSYY(NM1)+LSYY(NM2)+LSYY(NM3))/3.d0
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
     &         +TANPHIAvg*(LSXXEle - LSYYEle)
         ENDIF

C...     Compute the spatial gradients of the velocity dispersion terms if 3D

         IF (C3D) THEN                         !3D bottom friction
            DispX = (DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &          +  DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY = (DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &          +  DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
            DispXAvg = IFNLCT*DispX
            DispYAvg = IFNLCT*DispY
!         ELSEIF (CBaroclinic.and.abs(IDEN).eq.6) THEN
!            ! If 2D and baroclinic and getting information from a 3D BC
!            ! model (the gradient has already been calculated)
!            DispXAvg=(H1N1*VIDISPDXOH(NM1)+H1N2*VIDISPDXOH(NM2)
!     &                                       +H1N3*VIDISPDXOH(NM3))/3.d0
!            DispYAvg=(H1N1*VIDISPDYOH(NM1)+H1N2*VIDISPDYOH(NM2)
!     &                                       +H1N3*VIDISPDYOH(NM3))/3.d0
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
!JLW: add subgrid
         IF(subgrid_level0)THEN
            DPAvg = (H1N1+H1N2+H1N3)/3.d0
            DPAvgS0 = (H0N1+H0N2+H0N3)/3.d0
         ELSE
C.....   DW, Nov 2019 
            DPAvg = (DP(NM1)+DP(NM2)+DP(NM3))/3.d0 + ! BAVG if GWCE_HDP = .FALSE. 
     &           IFNLFA*IFNL_HDP*(E1N1 + E1N2 + E1N3)/3.d0 ; ! H^{s}_{avg} otherwise
            
            DPAvgS0 = (DP(NM1)+DP(NM2)+DP(NM3))/3.d0 + ! BAVG if GWCE_HDP = .FALSE. 
     &           IFNLFA*IFNL_HDP*(E0N1 + E0N2 + E0N3)/3.d0 ; ! H^{s-1}_{avg} otherwise
         ENDIF
C.....   WJP adjust the dispersion equation through
C        Ali et al. formula
         if (CAliDisp) then
            Ma2    = G*DPavg/Cs2
            Cfac   = (1.0D0 - 0.25D0*Ma2 - Ad*DPavg**Bd)**2D0
            Ma2    = G*DPavg/Cs2
            CfacS0 = (1.0D0 - 0.25D0*Ma2 - Ad*DPavgS0**Bd)**2D0
         endif
         ! DMW 202206 Apply slope limiting to gravity
         GDPAvgOAreaIE4 = ALPHAL(IE)*G*DPAvg*Cfac/AreaIE4
         ! DMW 202206 Apply slope limiting to gravity
         GDPAvgOAreaIE4_S0 = ALPHAL(IE)*G*DPAvgS0*CfacS0/AreaIE4
C....  END DW & WJP
         HAvg=(H1N1+H1N2+H1N3)/3.d0
         ! DMW 202206 Apply slope limiting to gravity
         GHAvg=ALPHAL(IE)*G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.d0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.d0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
C        WJP: Add on the contribution from spherical correction term
C        from the momentum equation
         ! (tan \phi V_{lambda})/R + f
         CorifAvg = CorifAvg + IFNLCT*TANPHIAvg*U1Avg
C
C...... Eliminate Coriols at condensed nodes  10/17/2022 sb
         IF (LoadCondensedNodes.AND.ROT_STATUS.EQ.0) CorifAvg = 0.D0
C
C...... DW absorbing layer
C
         SigmaAvg = (SgN1_eta + SgN2_eta + SgN3_eta)/3.D0
         Tau0SigmaAvg = (SgN1_eta*T0N1 + SgN2_eta*T0N2 + 
     &                   SgN3_eta*T0N3)/3.d0

         ! Use in the momentum equation parts
         SigmaQxAvg = (SgN1_mnx*QX1N1 +
     &         SgN2_mnx*QX1N2 + SgN3_mnx*QX1N3)/3.d0 ;
         SigmaQyAvg = (SgN1_mny*QY1N1 +
     &         SgN2_mny*QY1N2 + SgN3_mny*QY1N3)/3.d0 ;

         AvgSigmaHAbsU = (SgN1_mnx*H1N1*U1AN1 +
     &        SgN2_mnx*H1N2*U1AN2 + SgN3_mnx*H1N3*U1AN3)/3.D0 ; 

         AvgSigmaHAbsV = (SgN1_mny*H1N1*V1AN1 +
     &        SgN2_mny*H1N2*V1AN2 + SgN3_mny*H1N3*V1AN3)/3.D0 ; 

         IF ((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            fwsx1 = 1.D0 ; fwsx2 = 1.D0 ; fwsx3 = 1.D0 ;
            fwsy1 = 1.D0 ; fwsy2 = 1.D0 ; fwsy3 = 1.D0 ;
            
C.....      DMW 202207 tail off wind forcing in very shallow water
            IF (WINDLIM.eqv..true.) THEN
               CALL windLimiter(H1N1,fwsx1)
               fwsy1 = fwsx1
               CALL windLimiter(H1N2,fwsx2)
               fwsy2 = fwsx2
               CALL windLimiter(H1N3,fwsx3)
               fwsy3 = fwsx3
            ENDIF
C.....      DMW            
            
            IF ( NO_MET_IN_SPONGE ) THEN
               IF ( absorblayer_sigma_mnx(NM1,1) > 1e-9 ) fwsx1 = 0.D0
               IF ( absorblayer_sigma_mnx(NM2,1) > 1e-9 ) fwsx2 = 0.D0
               IF ( absorblayer_sigma_mnx(NM3,1) > 1e-9 ) fwsx3 = 0.D0
           
               IF ( absorblayer_sigma_mny(NM1,1) > 1e-9 ) fwsy1 = 0.D0
               IF ( absorblayer_sigma_mny(NM2,1) > 1e-9 ) fwsy2 = 0.D0
               IF ( absorblayer_sigma_mny(NM3,1) > 1e-9 ) fwsy3 = 0.D0
            END IF
         
            WSXAvg=(fwsx1*WSX1(NM1)+fwsx2*WSX1(NM2)+fwsx3*WSX1(NM3))/3.d0
            WSYAvg=(fwsy1*WSY1(NM1)+fwsy2*WSY1(NM2)+fwsy3*WSY1(NM3))/3.d0
!JLW: adding subgrid
            IF(subgrid_level0)THEN
               WSXAvg = PHIAVG2*WSXAvg
               WSYAvg = PHIAVG2*WSYAvg
            ENDIF
         ENDIF
C------- END DW
C 
         IF (CBaroclinic) THEN
#ifndef NOFSBPG
            BCXAvg = (H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                    +H1N3*VIDBCPDXOH(NM3))/3.d0
            BCYAvg = (H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                    +H1N3*VIDBCPDYOH(NM3))/3.d0
#else
            BCXAvg = (DP(NM1)*VIDBCPDXOH(NM1) + DP(NM2)*VIDBCPDXOH(NM2)
     &                + DP(NM3)*VIDBCPDXOH(NM3))/3.d0
            BCYAvg = (DP(NM1)*VIDBCPDYOH(NM1) + DP(NM2)*VIDBCPDYOH(NM2)
     &                + DP(NM3)*VIDBCPDYOH(NM3))/3.d0
#endif
         ENDIF

C...     Compute additional partial factors
C........DW: begin add absorbing layer term
         Tau0Avg_S = Tau0Avg + 
     &     ((A00 + B00)*2.D0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
   
         MsFacR=AreaIE*(1.d0/DT-Tau0Avg_S/2.d0)/DT/12.d0
         MsFacR_S = (-(A00 + B00 + C00)*AreaIE*
     &          Tau0SigmaAvg/12.D0)*sponge_dis_mthd 

         MsFacRS_ABC00 = (Tau0SigmaAvg*AreaIE/12.D0)*sponge_dis_mthd
         MsFacRS_T = ((SigmaAvg/(2.D0*DT))*AreaIE/12.D0)*sponge_dis_mthd
C......... end DW
C
         ! DMW 2022/06 Apply slope limiting to gravity
         GOAreaIE4=ALPHAL(IE)*G/AreaIE4
!JLW: adding the wet area fraction to the time derivative terms and the
!water surface gradient term
         IF(subgrid_level0)THEN
            MSFacR = PHIAVG2*MSFacR - 
     &         dphidt*AreaIE*((PHIAVG2-PHIAVG1)/2.d0)/DT/12.d0
!JLW: removing the following subgrid addition because incorrect
!            GOAreaIE4 = GOAreaIE4*PHIAVG2
         ENDIF
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
         A00pB00=A00+B00

C...     Compute the JX, JY terms less the advection terms
         JXAvg = CorifAvg*QY1Avg
     &          -IFNLFA*GOAreaIE4*(1 - IFNL_HDP)*(
     &               E1N1SQ*FDX1 +E1N2SQ*FDX2 + E1N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1-Off1N1)*FDX1
     &          +(PR1N2-TiPN2-Off1N2)*FDX2+(PR1N3-TiPN3-Off1N3)*FDX3)
     &          +WSXAvg-BSXAvg+MXAvg-DispXAvg-BCXAvg+Tau0QXAvg 
C--- begin: DW, absorbing layer
     &          - SigmaQxAvg*sponge_dis_mthd
     &          + AvgSigmaHAbsU*sponge_dis_mthd 
C--- end DW


         JYAvg =-CorifAvg*QX1Avg
     &          -IFNLFA*GOAreaIE4*(1 - IFNL_HDP)*(
     &                E1N1SQ*FDY1+E1N2SQ*FDY2+E1N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1-Off1N1)*FDY1
     &          +(PR1N2-TiPN2-Off1N2)*FDY2+(PR1N3-TiPN3-Off1N3)*FDY3)
     &          +WSYAvg-BSYAvg+MYAvg-DispYAvg-BCYAvg+Tau0QYAvg
C--- begin: DW, absorbing layer
     &          -SigmaQYAvg*sponge_dis_mthd 
     &          +AvgSigmaHAbsV*sponge_dis_mthd
C--- end DW

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
!JLW: adding advection corrections
           IF(subgrid_level1)THEN
              JXAvg = JXAvg - IFNLCT*(
     &                (CADV1*U1N1*QX1N1*FDX1+CADV2*U1N2*QX1N2
     &                *FDX2+CADV3*U1N3*QX1N3*FDX3)
     &                +(CADV1*U1N1*QY1N1*FDY1+CADV2*U1N2*QY1N2
     &                *FDY2+CADV3*U1N3*QY1N3*FDY3))
     &                /AreaIE2
              JYAvg = JYAvg - IFNLCT*(
     &                (CADV1*V1N1*QX1N1*FDX1+CADV2*V1N2*QX1N2
     &                *FDX2+CADV3*V1N3*QX1N3*FDX3)
     &                +(CADV1*V1N1*QY1N1*FDY1+CADV2*V1N2*QY1N2
     &                *FDY2+CADV3*V1N3*QY1N3*FDY3))
     &                /AreaIE2
           ELSE
              JXAvg = JXAvg - IFNLCT*(
     &                (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2
     &                +U1N3*QX1N3*FDX3)
     &                +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2
     &                +U1N3*QY1N3*FDY3))
     &                /AreaIE2
              JYAvg = JYAvg - IFNLCT*(
     &                (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2
     &                +V1N3*QX1N3*FDX3)
     &                +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2
     &                +V1N3*QY1N3*FDY3))
     &                /AreaIE2
            ENDIF
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

C...... BEG: DW, DEC 2019         
C...     Needs multiplying by cos(lat) in case of the correction 
C        Node by node
         SPM1 = SFCT(NM1) ; 
         SPM2 = SFCT(NM2) ; 
         SPM3 = SFCT(NM3) ;
C...... END DW         

C...     Assemble forcing for node NM1 (local index j=1)

         Temp_LV_A1=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN1*SPM1 + OffDiag*(ESN2*SPM2+ESN3*SPM3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*( C00*( SFCXADJ*FDX1*E0XGrad2A+ 
     &                                    SFCYADJ*FDY1*E0YGrad2A ) ) 
     &        - GDPAvgOAreaIE4*(A00pB00*( SFCXADJ*FDX1*E1XGrad2A+
     &                                    SFCYADJ*FDY1*E1YGrad2A ) )

C...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX1 + SFCYADJ*JYAvg*FDY1)*0.5d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar
C
C..... DW: begin
     &        +(OnDiag*E0N1*SPM1 + OffDiag*(E0N2*SPM2+E0N3*SPM3))*MsFacR_S
     &        +(onDiag*ESAN1*SPM1 + OffDiag*(ESAN2*SPM2+ESAN3*SPM3))*MsFacRS_T
     &        +(onDiag*ABC00_EAN1*SPM1 + 
     &                OffDiag*(ABC00_EAN2*SPM2 + ABC00_EAN3*SPM3))*MsFacRS_ABC00
C..... END DW
C

C...     Assemble forcing for node NM2 (local index j=2)
         Temp_LV_A2=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN2*SPM2 + OffDiag*(ESN1*SPM1+ESN3*SPM3))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*( C00*( SFCXADJ*FDX2*E0XGrad2A +
     &                                    SFCYADJ*FDY2*E0YGrad2A ) ) 
     &        - GDPAvgOAreaIE4*(A00pB00*( SFCXADJ*FDX2*E1XGrad2A +
     &                                    SFCYADJ*FDY2*E1YGrad2A ) )

C...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX2 + SFCYADJ*JYAvg*FDY2)*0.5d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar 
C
C..... DW: begin
     &        +(OnDiag*E0N2*SPM2 + OffDiag*(E0N1*SPM1+E0N3*SPM3))*MsFacR_S
     &        +(onDiag*ESAN2*SPM2 + OffDiag*(ESAN1*SPM1+ESAN3*SPM3))*MsFacRS_T
     &        +(onDiag*ABC00_EAN2*SPM2 + 
     &                 OffDiag*(ABC00_EAN1*SPM1 + ABC00_EAN3*SPM3))*MsFacRS_ABC00
C..... END DW
C


C...     Assemble forcing for node NM3 (local index j=3)

         Temp_LV_A3=

C...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN3*SPM3 + OffDiag*(ESN1*SPM1+ESN2*SPM2))*MsFacR

C...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*(  C00*( SFCXADJ*FDX3*E0XGrad2A +
     &                                     SFCYADJ*FDY3*E0YGrad2A ) )
     &        - GDPAvgOAreaIE4*( A00pB00*( SFCXADJ*FDX3*E1XGrad2A +
     &                                     SFCYADJ*FDY3*E1YGrad2A ) )

C...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX3 + SFCYADJ*JYAvg*FDY3)*0.5d0

C...     Spatially varying Tau0 terms
     &        +Tau0SpaVar
C
C..... DW: begin, absorbing layer
     &        +(OnDiag*E0N3*SPM3 + OffDiag*(E0N1*SPM1+E0N2*SPM1))*MsFacR_S
     &        +(OnDiag*ESAN3*SPM3 + OffDiag*(ESAN1*SPM1+ESAN2*SPM2))*MsFacRS_T
     &        +(OnDiag*ABC00_EAN3*SPM3 + 
     &               OffDiag*(ABC00_EAN1*SPM1 + ABC00_EAN2*SPM2))*MsFacRS_ABC00
C..... END DW
C

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

      end do                  !End of elemental loop

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
ckmd  Need to save z(s-1) and etas(s-1) for the corrector loop
      IF(CPRECOR)THEN
        ETAS0(1:NP) = ETAS(1:NP)
        ETA0(1:NP) = ETA1(1:NP)
      ENDIF
      ETA1(1:NP) = ETA2(1:NP)
      ETA2(1:NP) = 0.0D0
      H1(1:NP) = H2(1:NP)

!...  Apply boundary conditions
      call apply_boundary_conditions(IT, TimeLoc, TimeH, dt, invertedBarometerOnElevationBoundary,
     &                               esbin1, esbin2, etime1, etime2, etratio, etiminc, eta1, eta2,
     &                               etas, h1, elevdisc, en0, en1, en2, rampelev, rampmete, ep, gwce_lv)

C...
C...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
C...

#ifdef CMPI
!...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      IF (ILump==0) THEN ! default, fully consistent LHS
         IF(NPERSEG>0)THEN
             CALL UPDATER_W_PERBC(GWCE_LV,COEF(1,1),DUMY1,2)
         ELSE
             CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
         ENDIF
      ELSE  ! lumped LHS
         CALL UPDATER(GWCE_LV,COEFD,DUMY1,2)
      ENDIF
#endif

      call apply_subdomain_boundary_nodes(subdomainOn, enforceBN)

      call apply_vew1d_and_condensed_nodes(LoadCondensedNodes, ILump, NFLUXIB64_GBL, NBOU, NVELL, NBV, IBCONN,
     &                                     LBCODEI, NODECODE, NListCondensedNodes, NNodesListCondensedNodes,
     &                                     ListCondensedNodes, IsSubmerged64, COEFD,
     &                                     COEFDTemp, COEFDTempMem, GWCE_LV, ETA1)

      IF (ILump.eq.0) THEN ! default, fully consistent LHS
C...  JCG ITERATIVE MATRIX SOLVER
         IPARM(1)=ITMAX
         CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &        IWKSP,NW,WKSP,IPARM,RPARM,IER)
         NUMITR=IPARM(1)
         DO I=1,NP
            ETA2(I)=NODECODE(I)*(ETAS(I)+dynamicWaterLevelCorrectionDelta(I))+ETA1(I) !COMPUTE NEW ELEVATIONS tgaContOffset
         END DO
      ELSE ! lumped LHS
         DO I = 1, NP
            IF (COEFDTemp(I).eq.0.0d0) THEN   ! COEFD --> COEFDTemp  08-11-2022 SB
               RDIAG = 0.0d0
            ELSE
               RDIAG = 1.0d0 / COEFDTemp(I)   ! COEFD --> COEFDTemp  08-11-2022 SB
            ENDIF
            ETAS(I) = GWCE_LV(I) * RDIAG
         ENDDO
         NUMITR=0
         DO I=1,NP
            ETA2(I)=NODECODE(I)*(ETAS(I)+dynamicWaterLevelCorrectionDelta(I))+ETA1(I) !COMPUTE NEW ELEVATIONS tgaContOffset
         END DO
      ENDIF

C...... DW, periodic bcs
      call update_periodic_sponge_layer_nodes(NPERSEG, NNPERBC, IPERCONN, ETAS, ETA2)
C     Switch back to the original element table 
      CALL SWITCH_ELTAB_PERBC( IDIREC = 2 )
C......

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif      
     
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      END SUBROUTINE GWCE_NEW
C**********************************************************************


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

      USE SIZES
      USE GLOBAL, ONLY : GWCE_LV, COEF, COEFD, OBCCOEF, ETAS, ETA1, ETA2,
     &   ESBIN1, ESBIN2, C2DDI, C3D, CBaroclinic, CGWCE_Advec_C1, H1, H2,
     &   CGWCE_Advec_C2, CGWCE_Advec_NC, CGWCE_LS_2PartQ, CGWCE_LS_2PartSQ,
     &   CGWCE_LS_2PartV, CGWCE_LS_2PartSV, CGWCE_LS_KGQ, CPRECOR,
     &   CSmag_Eh, Smag_Comp_Flag, CTIP, DT, FluxSettlingIT, IFNLCAT,
     &   IFNLCT, IFNLFA, NBFR, NODECODE, NRS, NSCREEN, NCCHANGE, IDEN,
     &   NWS, RAMPELEV, RAMPMETE, RES_BC_FLAG,
     &   SMAG_LOWER_LIM, SMAG_UPPER_LIM, TK, NOFF, UU1, VV1, QX1, QY1,
     &   PR1, TiP1, BSX1, BSY1, DUU1, DUV1, DVV1, WSX1, WSY1, VIDBCPDXOH,
     &   VIDBCPDYOH, CORIF, PER, FACE, FF, EFA, EMO, AMIG, PR2, LNM_BC,
     &   QN1, EN1, QN0, ElevDisc, QN2, EN0, iLump, BCFLAG_LNM, EN2,
     &   TKM, NPERSEG, NNPERBC, IPERCONN, UU0, VV0, QX0, QY0, TK0, TK2, 
     &   QX2, QY2, UU2, VV2, windlim
#ifdef CMPI
     &   , dumy1, dumy2
#endif
      use mod_logging, only: screenUnit
      USE ADC_CONSTANTS, ONLY: G, RHOWAT0
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, TotalArea, Areas, NeiTab,
     &                  NNEIGH, SFAC
      USE BOUNDARIES, ONLY : NETA, NFLUXB, NFLUXF, NFLUXGBC, NFLUXRBC,
     &   NFLUXIB, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3, NPEBC
      USE ITPACKV
      USE NodalAttributes, ONLY : FRIC, Tau0Var, HBREAK, FTHETA, FGAMMA,
     &     IFLINBF, IFNLBF, IFHYBF, EVM, LoadAdvectionState, advectlocal
#ifdef CMPI
      USE MESSENGER
#endif
      USE WIND, ONLY: windLimiter
      USE mod_gwce_bc_forcing, only: apply_boundary_conditions

      IMPLICIT NONE

      integer,intent(in) :: IT
      real(8),intent(in) :: TimeLoc
      real(8),intent(in) :: TimeH

      INTEGER IE, I, J                           !local loop counters
      INTEGER  NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      REAL(8) A00pB00
      REAL(8) BCXAvg, BCYAvg
      REAL(8) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(8) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(8) CorifAvg
      REAL(8) DPAvg, GDPAvgOAreaIE4
      REAL(8) DispX, DispY, DispXAvg, DispYAvg
      REAL(8) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(8) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(8) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(8) ESN1, ESN2, ESN3, ESAvg
      REAL(8) EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(8) EVMEle, EVMSmag
      REAL(8) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(8) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg
      REAL(8) H2OTotalArea
      REAL(8) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(8) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(8) MsFacR
      REAL(8) MX, MY, MXAvg, MYAvg
      REAL(8) JXAvg, JYAvg
      REAL(8) Pr1N1, Pr1N2, Pr1N3
      REAL(8) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(8) SFacAvg
      REAL(8) T0N1,T0N2, T0N3
      REAL(8) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(8) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(8) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(8) TiPN1, TiPN2, TiPN3
      REAL(8) UV0, UV1, UV2
      REAL(8) U1N1,U1N2,U1N3, U1Avg
      REAL(8) V1N1,V1N2,V1N3, V1Avg
      REAL(8) WSXAvg, WSYAvg

      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3

ckmd   Added in parameters for the pc algorithm
      REAL(8) BSX0N1, BSX0N2, BSX0N3, BSY0N1
      REAL(8) BSY0N2, BSY0N3, BSX0Avg, BSY0Avg
      REAL(8) BSX2N1, BSX2N2, BSX2N3, BSY2N1
      REAL(8) BSY2N2, BSY2N3, BSX2Avg, BSY2Avg
      REAL(8) E2N1,E2N2,E2N3
      REAL(8) E0N1SQ, E0N2SQ, E0N3SQ
      REAL(8) E2N1SQ, E2N2SQ, E2N3SQ
      REAL(8) H0N1, H0N2, H0N3, H00
      REAL(8) QX0N1, QX0N2, QX0N3, QY0N1, QY0N2, QY0N3, QX0Avg, QY0Avg
      REAL(8) QX2N1, QX2N2, QX2N3, QY2N1, QY2N2, QY2N3, QX2Avg, QY2Avg
      REAL(8) Tau0QX0Avg, Tau0QY0Avg, Tau0QX2Avg, Tau0QY2Avg
      REAL(8) Tau0SpaVar0, Tau0SpaVar2
      REAL(8) U0N1,U0N2,U0N3, U0Avg
      REAL(8) V0N1,V0N2,V0N3, V0Avg
      REAL(8) U2N1,U2N2,U2N3, U2Avg
      REAL(8) V2N1,V2N2,V2N3, V2Avg
      REAL(8) timewtgwce0,timewtgwce1,timewtgwce2,timeagflag
      REAL(8) :: arg, argj
      REAL(8) :: celerity
      REAL(8) :: etratio
      INTEGER :: ier
      INTEGER :: nbdj
      REAL(8) :: rff
      REAL(8) :: qforcei, qforcej
      REAL(8):: fwsx1, fwsx2, fwsx3
      REAL(8):: fwsy1, fwsy2, fwsy3

      call setMessageSource("gwce_new_pc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

ckmd    Must reset the result vector to zero before recomputing
ckmd  the next time level.
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
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
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
C..............DMW202401 Use saved H2
               H2OTotalArea=H2(I)/TotalArea(I)
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
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN
               LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.d0
               LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.d0
               LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.d0
               LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.d0
               EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
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
ckmd  Added in the three time levels for the tau term.
ckmd  Every term is updated for the three time levels.
ckmd
            UV0=SQRT(UU0(I)*UU0(I)+VV0(I)*VV0(I))
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            UV2=SQRT(UU2(I)*UU2(I)+VV2(I)*VV2(I))
            H00=DP(I)+IFNLFA*ETA0(I)

C...........DMW202401 Use saved H1 and H2  
            TK0(I)=FRIC(I)*(IFLINBF + (UV0/H00)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H00)**FTHETA)**(FGAMMA/FTHETA)))
            TK(I)=FRIC(I)*(IFLINBF + (UV1/H1(I))*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H1(I))**FTHETA)**(FGAMMA/FTHETA)))
            TK2(I)=FRIC(I)*(IFLINBF + (UV2/H2(I))*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H2(I))**FTHETA)**(FGAMMA/FTHETA)))
       END DO

ckmd      Added in the time weights
C...     Time weights for the nonlinear terms in the GWCE for
C...       the corrector step
          timewtgwce0=0.33d0
          timewtgwce1=0.34d0
          timewtgwce2=0.33d0
          timeagflag=1.0d0

C...  Assemble the GWCE RHS except for the boundary integral terms
ckmd  Renumber the GWCE loop for the corrector step

      DO 1038 IE=1,NE

C...     Set nodal values for each element
ckmd
ckmd  Define the needed product terms at three time levels
ckmd


Corbitt 120322: Localized Advection
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
C........DMW202401 Use saved H1 and H2
         H1N1=H1(NM1)
         H1N2=H1(NM2)
         H1N3=H1(NM3)
         H2N1=H2(NM1)
         H2N2=H2(NM2)
         H2N3=H2(NM3)
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
         IF (C3D) THEN !3D velocity dispersion
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
            fwsx1 = 1.D0 ; fwsx2 = 1.D0 ; fwsx3 = 1.D0 ;
            fwsy1 = 1.D0 ; fwsy2 = 1.D0 ; fwsy3 = 1.D0 ;
            
C.....      DMW 202207 tail off wind forcing in very shallow water
            IF (WINDLIM.eqv..true.) THEN
               CALL windLimiter(H1N1,fwsx1)
               fwsy1 = fwsx1
               CALL windLimiter(H1N2,fwsx2)
               fwsy2 = fwsx2
               CALL windLimiter(H1N3,fwsx3)
               fwsy3 = fwsx3
            ENDIF
C.....      DMW            
         
            WSXAvg=(fwsx1*WSX1(NM1)+fwsx2*WSX1(NM2)+fwsx3*WSX1(NM3))/3.d0
            WSYAvg=(fwsy1*WSY1(NM1)+fwsy2*WSY1(NM2)+fwsy3*WSY1(NM3))/3.d0
         ENDIF
         IF (C3D.or.(CBaroclinic.and.abs(IDEN).eq.5)) THEN !3D velocity dispersion
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
ckmd  Already did this and haven't advanced in time yet, so just
ckmd  comment out the save elevation line. But we do need to zero
ckmd  out eta2 because it's involved in the summation below.
      Eta2(:) = 0.0d0

C...  At elevation boundary condition nodes, determine the elevation at
C...  the s+1 time step
      call apply_boundary_conditions(IT, TimeLoc, TimeH, dt, invertedBarometerOnElevationBoundary,
     &                               esbin1, esbin2, etime1, etime2, etratio, etiminc, eta1, eta2,
     &                               etas, h1, elevdisc, en0, en1, en2, rampelev, rampmete, ep, gwce_lv)

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
ckmd  Eta2 values are now corrected elevations at time level s+1.
ckmd
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE GWCE_NEW_PC
C**********************************************************************

C*************************************************************************
C   Subroutine to check if elemental slope limiting is needed
C       Chris Massey, USACE-ERDC-CHL Dec. 8, 2014 -- Removed this
C         section from timestep subroutine and put it in its own
C         subroutine.
C     jgf53.dev: Moved to GWCE module. 
C*************************************************************************
      subroutine check_slopes(it,TimeLoc)
      USE SIZES, ONLY : mnproc, myproc
      use global, only : eta2,nodecode,NOFF,ESLCOUNT,nodes_lg, IFSFM
      use mod_logging, only: screenUnit,
     &      setMessageSource, unsetMessageSource, scratchMessage,
     &      allMessage,logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      use mesh, only : nm, ne, areas, FDXE, FDYE, 
     &      SFacEle, SFMYEle, SFMXEle
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active,
     &     elemental_slope_limiter_grad_max,
     &     elemental_slope_limiter_max_exceeded,
     &     LoadCondensedNodes, NCondensedNodes

      implicit none
      INTEGER, intent(in) :: IT
      Real(8), intent(in) :: TimeLoc
      INTEGER IE,I
      INTEGER NM1, NM2, NM3, NM123
      INTEGER NC1, NC2, NC3, NCEle, NCI
      REAL(8) DEta2DX, DEta2DY, DEta2Mag
      REAL(8) SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac 
      REAL(8) FDX1,FDX2,FDX3,FDY1,FDY2,FDY3
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
         AreaIE2 = Areas(IE)
C.... Beg WP:
         SFacAvg = SFacEle(IE)
C..... BEG DW/WJP
         SFmxAvg = SFMXEle(IE) ; 
         SFmyAvg = SFMYEle(IE) ; 
         sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ;
         sfdyfac = (1 - IFSFM)*1.0D0 + IFSFM*SFmyAvg ;
C..... END DW/WJP
         FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*MX !b1*mx
         FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*MX !b2*mx
         FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*MX !b3*mx
         FDY1 = FDYE(1,IE)*sfdyfac ; !c FDY1=(X(NM3)-X(NM2))*MY !a1*my
         FDY2 = FDYE(2,IE)*sfdyfac ; !c FDY2=(X(NM1)-X(NM3))*MY !a2*my
         FDY3 = FDYE(3,IE)*sfdyfac ; !c FDY3=(X(NM2)-X(NM1))*MY !a3*my         
C.... end WP
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
            ! If this is one of the condensed nodes, go to the next node. 10/12/2023 sb
            IF (LoadCondensedNodes) THEN
                IF (NCondensedNodes(NM(IE,I)) > 0) CYCLE
            ENDIF
            ! Compare the elemental slope to the maximum elemental gradient.
            grad_check: IF (dEta2Mag.GE.
     &           ABS(elemental_slope_limiter_grad_max(NM(IE,I)))) THEN
               ! jgf51.51: Log the fulldomain node number.
               nodeNumber = nm(ie,i)
               if (mnproc.gt.1) then
                  nodeNumber = nodes_lg(nm(ie,i))
               endif
               ! zc - If gradmax is positive or zero, activate slope
               ! limiting. 
               grad_case: if (elemental_slope_limiter_grad_max(nm(ie,i))
     &               .ge.0.0d0) then
                  ! If it is the first time that the slope limiter
                  ! has been activated at this node, write a log 
                  ! message. 
                  if (eslcount(nm(ie,i)).eq.0) then
                     write(scratchMessage,1983) nodeNumber,dEta2Mag, 
     &                     elemental_slope_limiter_grad_max(NM(IE,I)),
     &                     it, timeLoc
                     call allMessage(INFO,scratchMessage)
                     eslcount(nm(ie,i)) = 1 ! for output file
                  endif
                  elemental_slope_limiter_active(NM(IE,I)) = .true.
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
               endif grad_case
            endif grad_check
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
C     jgf53.dev: Moved to GWCE module. 
C***********************************************************************
C
      SUBROUTINE APPLY_SLOPE_LIMITS(ETA2Lim,LocNP)
      use global, only : nodecode,NOFF,IFNLFA
      use mod_logging, only:
     &  setMessageSource, unsetMessageSource, allMessage,
     &  logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      use mesh, only : ne, nm, areas, totalArea
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active
      implicit none
      integer :: IE,NM1,NM2,NM3,NC1,NC2,NC3,NCEle
      integer, intent(in) :: LocNP
      REAL(8) :: EtaN1,EtaN2,EtaN3,EtaN123
      real(8), intent(inout) :: Eta2lim(LocNP)
      REAL(8) :: AreaEle
      REAL(8), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
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


C----------------------------------------------------------------------
C----------------------------------------------------------------------
      end module gwce
C----------------------------------------------------------------------
C----------------------------------------------------------------------
