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
!  MODULE TIMESTEP
!-----------------------------------------------------------------------
!> @brief Main timestepping module for ADCIRC hydrodynamic calculations
!>
!> This module provides the primary timestepping routine for ADCIRC,
!> orchestrating the solution of the GWCE and momentum equations. It handles
!> time level shifts, ramp functions, boundary conditions, meteorological
!> forcing, wetting/drying, baroclinic pressure gradients, and output.
!>
!> The timestepping module supports multiple algorithm configurations:
!> - C2DDI: 2D Depth Integrated model
!> - C3D: 3D model
!> - Baroclinic transport calculations
!> - Predictor-corrector
!-----------------------------------------------------------------------
module mod_timestep

   implicit none

   private

   public :: TIMESTEP

   logical :: EtaDisc_Fill = .true.
   real(8) :: hotstart_start_time = 0.0d0
   real(8) :: hotstart_reference_time = 0.0d0

   real(8), parameter :: eps = epsilon(1.0d0)

contains

!-----------------------------------------------------------------------
!> @brief Main timestepping subroutine for ADCIRC hydrodynamic model
!>
!> This subroutine performs a single timestep of the ADCIRC model, including:
!> - Computing time variables and ramp functions
!> - Shifting solution arrays to previous time levels
!> - Applying meteorological forcing (wind, pressure)
!> - Computing boundary conditions (elevation, flux, barriers)
!> - Solving the GWCE for water surface elevation
!> - Applying wetting and drying algorithm
!> - Solving 2D or 3D momentum equations
!> - Applying absorbing sponge layers
!> - Collecting output statistics and writing hot starts
!>
!> Algorithm selection is controlled by logical flags set in READ_INPUT
!> based on the IM parameter:
!> - C2DDI: 2D Depth Integrated model
!> - C3D: 3D model (C3DDSS for stress form, C3DVS for velocity form)
!> - C2D_BTrans/C3D_BTrans: Baroclinic transport
!> - C2D_PTrans/C3D_PTrans: Passive transport
!> - CBaroclinic: Include baroclinic terms
!> - CPRECOR: Predictor-corrector algorithm
!>
!> @param[in] IT        Current timestep number
!> @param[in] ITIME_BGN Starting timestep number (for hot starts)
!> @param[out] TimeLoc  Simulation time at current timestep (seconds)
!-----------------------------------------------------------------------
   subroutine TIMESTEP(IT, ITIME_BGN, TimeLoc)

      use GLOBAL, only: METONLY, NDDT, CTIP, C2DDI, C3DVS, C2D_PTrans, NT, &
                        NWS, WSX2, WSY2, PR2, WVNXOUT, WVNYOUT, &
                        usingDynamicWaterLevelCorrection
      use GWCE, only: solveGWCE
      use MOMENTUM2D, only: solve2DMomentumEq
      use MOMENTUM3D, only: solve3DVSMomentumEq
      use WETDRY, only: computeWettingAndDrying
      use WRITE_OUTPUT, only: writeOutput2D, writeHotStart, &
                              writeWarnElev, collectInundationData, collectMinMaxData
      use WIND, only: CLOSE_MET_FILES, getMeteorologicalForcing
      use mod_terminate, only: terminate
      use mod_logging, only: t_log_scope, init_log_scope

#ifdef CMPI
      use HSWRITER, only: writeHotstart_through_HSwriter !st3 100711 for hsfile
#endif

      use SPONGELAYER, only: sponge_opsplit0, &
                             adjust_sponge_sigma, getabslayerext, sponge_shift_soln, &
                             LoadAbsLayerSigma
      use subgrid, only: getVertLookup
      use mod_baroclinic, only: computeBaroclinicPressureGradient, updateBaroclinicBoundaryConditions

      implicit none
      integer, intent(in) :: IT
      integer, intent(in) :: ITIME_BGN
      real(8), intent(out) :: TimeLoc
      real(8) :: TimeH

      LOG_SCOPE_TRACED("timestep", TIMESTEP_TRACING)

      call computeTimeVariables(IT, TimeLoc, TimeH)
      call computeRampFunctions(IT, TimeLoc)
      call shiftTimeLevels(IT)
      if (LoadAbsLayerSigma) call initializeAbsorbingLayer(TimeLoc, TimeH)
      if (NDDT /= 0) call updateTimeDependentBathymetry(TimeLoc)
      call computeBottomFriction(IT)
      call updateIceFields(TimeLoc)
      if (usingDynamicWaterLevelCorrection) call updateWaterLevelCorrections(TimeLoc)
      call getMeteorologicalForcing(nws, timeloc, wsx2, wsy2, pr2, wvnxout, wvnyout)
      call updateWaveRadiationStress(TimeLoc)

      ! Skip the GWCE and momentum solutions if only meteorological output
      ! was requested
      if (.not. METONLY) then

         call updateBaroclinicBoundaryConditions(TimeLoc)
         if (CTIP) call applyTidePotentialForcing(TimeLoc, TimeH)
         call computeBaroclinicPressureGradient(TimeLoc)

         call applyBoundaryConditions(IT, TimeLoc, TimeH)

         call solveGWCE(it, ITIME_BGN, timeloc, timeh)

         call applyWettingAndDrying(IT)

         call updateFlowDepth()

         if (C2DDI) call solve2DMomentumEq(IT, TimeLoc)

         call applyAbsorbingSpongeLayer()

         if (C3DVS) call solve3DVSMomentumEq(IT, TimeLoc)

         if (C2D_PTrans) call SCALAR_TRANS_2D()

      end if

      call collectOutputStatistics(IT, TimeLoc, TimeH)
      call writeHotstartIfNeeded(IT, TimeLoc)
      call shiftBarrierHeights()

      call progressScreen(IT, TimeLoc)

      if (IT == NT) call CLOSE_MET_FILES()

   end subroutine TIMESTEP
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute tidal potential forcing at mesh nodes
!>
!> Computes the gravitational tide-generating potential at each node,
!> including Earth tide reduction factors. Supports both the full formula
!> (tidePotential%active) and the traditional constituent-based approach.
!>
!> For the traditional approach, tidal potential is computed as a sum of
!> constituents with amplitudes, frequencies, and phase corrections. The
!> self-attraction and loading (SAL) terms are also included.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine applyTidePotentialForcing(TimeLoc, TimeH)
      use mod_tidepotential, only: tidePotential
      use global, only: NTIP, TIP2, NTIF, RampTip, PERT, FACET, &
                        AMIGT, FFT, TPK, ETRF, SALTPHA, SALTAMP, L_N
      use mesh, only: np, SLAM
      implicit none
      real(8), intent(in) :: TimeLoc, TimeH
      real(8) :: ArgT, ArgTP, ArgSAlt, SALTMUL, TPMUL
      integer :: i, j, na, ncyc

      ! use the full formula
      if (tidePotential%active()) then
         TIP2 = tidePotential%compute(TimeLoc, NP, SLAM)*RampTip
      end if

      ! The traditional approach: a sum of different constituents
      ! TIP & SAL - for tidePotential%active() = .TRUE.
      ! SAL       - for tidePotential%active() = .FALSE.
      if (NTIP == 2 .or. .not. tidePotential%active()) then
         do J = 1, NTIF
            if (abs(PERT(J)) <= eps) then
               NCYC = 0
            else
               NCYC = int(timeh/PERT(J))
            end if

            ARGT = AMIGT(J)*(timeh - dble(NCYC)*PERT(J)) + FACET(J)
            TPMUL = RampTip*ETRF(J)*TPK(J)*FFT(J)
            SALTMUL = RampTip*FFT(J)

            !! WJP: We actually want to compare against diurnal and get 0,1,2
            !!      or higher species (was opposite beforehand)
            NA = min(nint(AMIGT(J)/7d-5), 2)

            !! WJP: Rewritten so that we have a general formula for all
            !!      species
            if (tidePotential%active()) then
               do I = 1, NP
                  ARGSALT = ARGT - SALTPHA(J, I); 
                  TIP2(I) = TIP2(I) + SALTMUL*SALTAMP(J, I)*cos(ARGSALT)
               end do
            else
               do I = 1, NP
                  ARGTP = ARGT + dble(NA)*SLAM(I)
                  ARGSALT = ARGT - SALTPHA(J, I)
                  TIP2(I) = TIP2(I) + TPMUL*L_N(NA, I)*cos(ARGTP) &
                            + SALTMUL*SALTAMP(J, I)*cos(ARGSALT)
               end do
            end if

         end do
      end if
   end subroutine applyTidePotentialForcing
   !-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!> @brief Shift solution arrays to previous time levels
!>
!> Shifts flux per unit width, depth averaged velocities, bottom stress,
!> wind stress, surface pressure and tidal potentials to previous time step.
!> Zeros out new forcing terms and load vectors. On the first timestep,
!> initializes H1 and H2 from bathymetry and initial water surface elevation.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine shiftTimeLevels(IT)
      use global, only: MOM_LV_X, MOM_LV_Y, QX1, QY1, QX2, QY2, &
                        UU1, VV1, UU2, VV2, CPRECOR, UU0, VV0, QX0, QY0, &
                        IM, TRANS_LV_B, TRANS_LV_A, NWS, NRS, &
                        WSX1, WSY1, WSX2, WSY2, PR1, PR2, CTIP, TIP1, TIP2, &
                        usingDynamicWaterLevelCorrection, &
                        dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, &
                        H1, H2, ETA2, IFNLFA, &
                        QN0, QN1, QN2, EN0, EN1, EN2
      use gwce, only: gwce_lv_reset
      use wind, only: PRBCKGRND_MH2O
      use mesh, only: DP
      use boundaries, only: NVEL
      implicit none

      integer, intent(in) :: IT

      ! On first timestep of a cold start, initialize H1 and H2
      if (IT == 1) then
         H2(:) = DP(:) + dble(IFNLFA)*ETA2(:)
         H1(:) = H2(:)
      end if

      ! Shift values in time for predictor-corrector algorithm. This must
      ! happen before level 1 is overwritten by level 2 below.
      if (CPRECOR) then
         UU0(:) = UU1(:)
         VV0(:) = VV1(:)
         QX0(:) = QX1(:)
         QY0(:) = QY1(:)
      end if

      MOM_LV_X(:) = 0d0
      MOM_LV_Y(:) = 0d0
      QX1(:) = QX2(:)
      QY1(:) = QY2(:)
      UU1(:) = UU2(:)
      VV1(:) = VV2(:)
      call gwce_lv_reset(0d0)

      ! Transport
      if (IM == 10) then
         TRANS_LV_B(:) = 0.d0
         TRANS_LV_A(:) = 0.d0
      end if

      ! Wind (& wave radiation stress if used)
      if ((NWS /= 0) .or. (NRS /= 0)) then
         WSX1(:) = WSX2(:)
         WSX2(:) = 0.d0
         WSY1(:) = WSY2(:)
         WSY2(:) = 0.d0
         PR1(:) = PR2(:)
         PR2(:) = PRBCKGRND_MH2O
      else
         PR2(:) = PRBCKGRND_MH2O
      end if

      ! Tidal potential forcing
      if (CTIP) then
         TIP1(:) = TIP2(:)
         TIP2(:) = 0.d0
      end if

      if (usingDynamicWaterLevelCorrection) then
         dynamicWaterLevelCorrection1(:) = dynamicWaterLevelCorrection2(:)
         dynamicWaterLevelCorrection2(:) = 0.d0
      end if

      ! Normal flow boundary conditions
      QN0(1:NVEL) = QN1(1:NVEL)
      QN1(1:NVEL) = QN2(1:NVEL)
      QN2(1:NVEL) = 0.d0
      EN0(1:NVEL) = EN1(1:NVEL)
      EN1(1:NVEL) = EN2(1:NVEL)
      EN2(1:NVEL) = 0.d0

   end subroutine shiftTimeLevels
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Save the current barrier heights into the previous time level
!>
!> Called at the end of the timestep, after the hot start file has been
!> written.
!-----------------------------------------------------------------------
   subroutine shiftBarrierHeights()
      use boundaries, only: NFLUXIB, NFLUXB
      use weir, only: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2
      implicit none

      if (NFLUXIB == 1) BARINHT1(:) = BARINHT2(:)
      if (NFLUXB == 1) BARLANHT1(:) = BARLANHT2(:)

   end subroutine shiftBarrierHeights
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Calculate temporal ramp functions for forcing terms
!>
!> Computes hyperbolic tangent ramp functions for gradual introduction of
!> boundary elevation, wind/pressure, tidal potential, and wave radiation
!> stress forcing. Handles special cases for flux settling periods.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeRampFunctions(IT, TimeLoc)
      use global, only: NRamp, Ramp, RampExtFlux, RampIntFlux, &
                        RampElev, RampTip, RampMete, RampWRad, &
                        DRamp, DRampExtFlux, DRampIntFlux, DRampElev, &
                        DRampTip, DRampMete, DRampWRad, DUnRampMete, &
                        FluxSettlingIT, DTDP, STATIM, CHOTHS
      use spongelayer, only: LoadAbsLayerSigma, Adjust_Sponge_Sigma, &
                             sponge_sigma_adjusted
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      real(8) :: TimeLocFluxSettling

      !     jgf46.08 Calculate ramp functions.
      !     jgf46.21 Modify to match behavior of 46.02
      if (NRamp == 0) then
         Ramp = 1.0d0
         RampExtFlux = 1.0d0
         RampIntFlux = 1.0d0
         RampElev = 1.0d0
         RampTip = 1.0d0
         RampMete = 1.0d0
         RampWRad = 1.0d0
      else
         Ramp = tanh((2.d0*TimeLoc/86400.d0)/DRamp)
         RampExtFlux = tanh((2.d0*TimeLoc/86400.d0)/DRampExtFlux)
         RampIntFlux = tanh((2.d0*TimeLoc/86400.d0)/DRampIntFlux)
         RampElev = tanh((2.d0*TimeLoc/86400.d0)/DRampElev)
         RampTip = tanh((2.d0*TimeLoc/86400.d0)/DRampTip)
         RampMete = tanh((2.d0*(TimeLoc/86400.d0 - DUnRampMete))/DRampMete)
         RampWRad = tanh((2.d0*TimeLoc/86400.d0)/DRampWRad)
      end if

      !     jgf46.21 If there is an external flux (i.e. river) boundary, turn
      !     off all forcings except the river flux forcing for the duration of
      !     the FluxSettlingTime. When the FluxSettlingTime has ended, turn
      !     all forcings back on.
      if (NRamp > 1 .and. FluxSettlingIT > 0) then
         if (IT < (FluxSettlingIT + 10)) then
            Ramp = 0.0d0
            RampIntFlux = 0.0d0
            RampElev = 0.0d0
            RampTip = 0.0d0
            RampMete = 0.0d0
            RampWRad = 0.0d0
         else
            TimeLocFluxSettling = dble(IT - FluxSettlingIT - 10)*DTDP &
                                  + STATIM*86400.d0
            if (CHOTHS) then
               TimeLocFluxSettling = dble(IT - FluxSettlingIT - 10)*DTDP &
                                     + (hotstart_start_time - hotstart_reference_time)*86400.d0
            end if

            Ramp = tanh((2.d0*TimeLocFluxSettling/86400.d0)/DRamp)
            RampIntFlux = tanh((2.d0 &
                                *TimeLocFluxSettling/86400.d0)/DRampIntFlux)
            RampElev = tanh((2.d0 &
                             *TimeLocFluxSettling/86400.d0)/DRampElev)
            RampTip = tanh((2.d0 &
                            *TimeLocFluxSettling/86400.d0)/DRampTip)
            if (NRamp /= 8) then
               RampMete = tanh((2.d0 &
                                *TimeLocFluxSettling/86400.d0)/DRampMete)
            end if
            RampWRad = tanh((2.d0 &
                             *TimeLocFluxSettling/86400.d0)/DRampWRad)
         end if
      end if

      ! jgf49.44: Cover the case where the ramp length is zero.
      if (DRamp < 1.0d-6) Ramp = 1.0d0
      if (DRampExtFlux < 1.0d-6) RampExtFlux = 1.0d0
      if (DRampIntFlux < 1.0d-6) RampIntFlux = 1.0d0
      if (DRampElev < 1.0d-6) RampElev = 1.0d0
      if (DRampTip < 1.0d-6) RampTip = 1.0d0
      if (DRampMete < 1.0d-6) RampMete = 1.0d0
      if (DRampWRad < 1.0d-6) RampWRad = 1.0d0

      !     Sponge layer sigma adjustment (DW)
      if (LoadAbsLayerSigma) then
         if (.not. sponge_sigma_adjusted) then
            call Adjust_Sponge_Sigma(abs(DTDP))
            sponge_sigma_adjusted = .true.
         end if
      end if

   end subroutine computeRampFunctions
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute time variables for current timestep
!>
!> Computes the master simulation time (TimeLoc), harmonic analysis
!> reference time (TimeH), and flux settling iteration count. Handles
!> timestep changes when hot starting from a different timestep size.
!>
!> @param[in]  IT      Current timestep number
!> @param[out] TimeLoc Simulation time at current timestep (seconds)
!> @param[out] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine computeTimeVariables(IT, TimeLoc, TimeH)
      use global, only: CHOTHS, FluxSettlingTime, FluxSettlingIT, &
                        DTDPHS, DTDP, STATIM, REFTIM, ITHS
      implicit none

      integer, intent(in) :: IT
      real(8), intent(out) :: TimeLoc, TimeH

      !     kmd48.33bc - changed for timestep changes in the hot start files
      !     jgf46.21 Combined flux/radiation b.c. for rivers
      if (CHOTHS) then
         FluxSettlingIT = int(FluxSettlingTime*86400.d0/DTDPHS)
      else
         FluxSettlingIT = int(FluxSettlingTime*86400.d0/DTDP)
      end if

      !     Compute master time referenced to beginning of model run
      !     TimeLoc: master simulation time
      !     TimeH: time for harmonic calculations (includes REFTIM offset)
      TimeLoc = dble(IT)*DTDP + STATIM*86400.d0
      TimeH = dble(IT)*DTDP + (STATIM - REFTIM)*86400.d0

      if (CHOTHS) then
         if ((ITHS + 1) == IT) then
            hotstart_start_time = (dble(IT - 1)*DTDPHS)/86400.d0
            hotstart_reference_time = (dble(IT - 1)*DTDP)/86400.d0
         end if
         TimeLoc = dble(IT)*DTDP + (hotstart_start_time - hotstart_reference_time)*86400.d0
         TimeH = dble(IT)*DTDP + ((hotstart_start_time - STATIM) - REFTIM)*86400.d0
      end if

      !     Harmonic time recalculation (includes hot start timestep changes)
      TimeH = dble(IT)*DTDP + (STATIM - REFTIM)*86400.d0
      if (CHOTHS) then
         if ((ITHS + 1) == IT) then
            hotstart_start_time = (dble(IT - 1)*DTDPHS)/86400.d0
            STATIM = (dble(IT - 1)*DTDP)/86400.d0
         end if
         TimeH = dble(IT)*DTDP + ((hotstart_start_time - STATIM) - REFTIM)*86400.d0
      end if

   end subroutine computeTimeVariables
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Initialize absorbing sponge layer boundary conditions
!>
!> Initializes the absorbing/sponge layer external boundary values. On
!> first entry, steps back in time to properly initialize the previous
!> time level states needed for the time integration scheme.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine initializeAbsorbingLayer(TimeLoc, TimeH)
      use global, only: NRamp, RampElev, DRampElev, DTDP
      use spongelayer, only: GETABSLAYEREXT, SPONGE_SHIFT_SOLN, sponge_initialized
      use NodalAttributes, only: GeoidOffset
      implicit none

      real(8), intent(in) :: TimeLoc, TimeH
      real(8) :: RampTMP
      integer :: istep

      if (sponge_initialized) then
         ! Normal timestep: get external boundary values
         call GETABSLAYEREXT(TimeLoc, TimeH, RampElev, GeoidOffset)
      else
         ! First entry: step back in time to initialize previous states
         RampTMP = RampElev
         do istep = 2, 1, -1
            if (NRamp > 0) then
               RampTMP = tanh((2.d0*(TimeLoc - dble(istep)*DTDP)/86400.d0) &
                              /DRampElev)
            end if
            RampTMP = max(0.d0, RampTMP)

            ! Get external values at previous time
            call GETABSLAYEREXT(TimeLoc - dble(istep)*DTDP, &
                                TimeH - dble(istep)*DTDP, RampTMP, GeoidOffset)
            call SPONGE_SHIFT_SOLN()
         end do
         ! Get current external values
         call GETABSLAYEREXT(TimeLoc, TimeH, RampElev, GeoidOffset)
         sponge_initialized = .true.
      end if

   end subroutine initializeAbsorbingLayer
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update bathymetry for time-dependent bathymetry runs
!>
!> Updates mesh bathymetry (DP) from input file for simulations with
!> time-varying bathymetry. Supports two modes:
!> - NDDT=1: Read bathymetry for all nodes
!> - NDDT=2: Read bathymetry for sparse node subset
!>
!> Bathymetry is linearly interpolated between snapshots:
!> @verbatim
!>    DP1                         DP2
!>  btime1                      btime2      btime2 = btime1 + btiminc
!>     |---------x-----------------|        btime_end = btime1 + bchgtiminc
!>            btime_end
!> @endverbatim
!>
!> Water surface elevation is adjusted to maintain total water column depth.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateTimeDependentBathymetry(TimeLoc)
      use global, only: NDDT, BTIME1, BTIME2, BTIMINC, BCHGTIMINC, &
                        BTIME_END, dp1, dp2, ETA2, NOLIFA, H0
      use mesh, only: NP, DP
      use mod_logging, only: logMessage, allMessage, INFO
      use global_io, only: NDDT2GET
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: bTRatio, DPTMP, DPTMP2
      integer :: I, IJK
      character(256) :: scratchMessage

      if (abs(NDDT) == 1) then
         ! Mode 1: Read bathymetry for all nodes
         if (TimeLoc > BTIME2) then
            BTIME1 = BTIME2
            BTIME2 = BTIME2 + BTIMINC
            BTIME_END = BTIME1 + BCHGTIMINC

            dp1(:) = dp2(:)
            ! Read new record for all nodes
            do I = 1, NP
               read (141, *) IJK, DP2(IJK)
            end do
            ! Ensure depths >= H0 if no wetting/drying
            if (NOLIFA == 0 .or. NOLIFA == 1) then
               where (DP2(1:NP) < H0) DP2(1:NP) = H0
            end if
            write (scratchMessage, '(A36,1X,E15.8,1X,A5)') &
               'BATHYMETRY RECORDS UPDATED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
         end if

         ! Interpolate during bathymetry change interval
         if (TimeLoc < BTIME_END) then
            bTRatio = (TimeLoc - BTIME1)/BCHGTIMINC
            do I = 1, NP
               DPTMP = bTRatio*(DP2(I) - DP1(I))
               DPTMP2 = DP1(I) + DPTMP
               DPTMP = DP(I) - DPTMP2
               DP(I) = DPTMP2
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

         ! Finalize bathymetry at end of change interval
         if (abs(TimeLoc - BTIME_END) <= eps) then
            write (scratchMessage, '(A42,1X,E15.8,1X,A5)') &
               'BATHYMETRY VALUES ARE NOW FIXED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
            do I = 1, NP
               DPTMP = DP2(I) - DP(I)
               DP(I) = DP2(I)
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

      else if (abs(NDDT) == 2) then
         ! Mode 2: Read bathymetry for sparse nodes
         if (TimeLoc > BTIME2) then
            BTIME1 = BTIME2
            BTIME2 = BTIME2 + BTIMINC
            BTIME_END = BTIME1 + BCHGTIMINC
            dp1(:) = dp2(:)
            ! Read new record for selected nodes only
            call NDDT2GET(141, DP2(:), -99999.d0)
            ! Ensure depths >= H0 if no wetting/drying
            if (NOLIFA == 0 .or. NOLIFA == 1) then
               where (DP2(1:NP) < H0) DP2(1:NP) = H0
            end if
            write (scratchMessage, '(A36,1X,E15.8,1X,A5)') &
               'BATHYMETRY RECORDS UPDATED AT TIME =', TimeLoc, '(SEC)'
            call allMessage(INFO, scratchMessage)
         end if

         ! Interpolate during bathymetry change interval
         if (TimeLoc < BTIME_END) then
            bTRatio = (TimeLoc - BTIME1)/BCHGTIMINC
            do I = 1, NP
               DPTMP = bTRatio*(DP2(I) - DP1(I))
               DPTMP2 = DP1(I) + DPTMP
               DPTMP = DP(I) - DPTMP2
               DP(I) = DPTMP2
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if

         ! Finalize bathymetry at end of change interval
         if (abs(TimeLoc - BTIME_END) <= eps) then
            write (scratchMessage, '(A42,1X,E15.8,1X,A5)') &
               'Bathymetry values are now fixed at time =', TimeLoc, '(sec).'
            call allMessage(INFO, scratchMessage)
            do I = 1, NP
               DPTMP = DP2(I) - DP(I)
               DP(I) = DP2(I)
               ETA2(I) = ETA2(I) - DPTMP
            end do
         end if
      end if

   end subroutine updateTimeDependentBathymetry
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute bottom friction coefficients
!>
!> Computes 2D or 3D bottom friction coefficients (TK) based on the
!> selected friction formulation. For baroclinic runs, also computes
!> internal wave drag and momentum dispersion coefficients (TKM).
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine computeBottomFriction(IT)
      use global, only: C2DDI, C3D, CBaroclinic, UU1, VV1, ETA2, &
                        IFNLFA, TK, TKM
      use global_3dvs, only: Q, SIGMA, NFEN, Z0B, KP
      use mesh, only: NP, DP
      use ADC_CONSTANTS, only: G
      use NodalAttributes, only: Apply2DBottomFriction, &
                                 Apply2DInternalWaveDrag, Apply3DBottomFriction, &
                                 Apply2DMomentumDisp
      implicit none

      integer, intent(in) :: IT

      ! 2DDI: Set up the 2D friction coefficient
      if (C2DDI) then
         call Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, IFNLFA, &
                                    NP, TK)
         ! Get wave drag coefficients into TKM (includes TK)
         if (CBaroclinic) then
            call Apply2DInternalWaveDrag(NP, TK, TKM, UU1, VV1, DP, IT)
            ! Vertically averaged momentum dispersion (power law)
            call Apply2DMomentumDisp(UU1, VV1, DP, ETA2, IFNLFA, &
                                     NP, TKM)
         else
            call Apply2DInternalWaveDrag(NP, TK, TKM)
         end if
      elseif (C3D) then
         ! 3D: Set up the 3D friction coefficient
         call Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, IFNLFA, &
                                    NP, TK, NFEN, Z0B, KP)
      end if

   end subroutine computeBottomFriction
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update ice concentration fields
!>
!> Reads ice concentration data from input file when needed and linearly
!> interpolates ice values (CICEOUT) for the current timestep.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateIceFields(TimeLoc)
      use global, only: NCICE, CICE_TIME1, CICE_TIME2, CICE_TIMINC, &
                        CICE1, CICE2, CICEOUT
      use mesh, only: NP
      use owi_ice, only: NCICE1_GET
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: CICE_TRatio
      integer :: I

      ! Read new ice data if needed
      if (NCICE == 12) then
         if (TimeLoc > CICE_TIME2) then
            CICE_TIME1 = CICE_TIME2
            CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
            CICE1(:) = CICE2(:)
            call NCICE1_GET(CICE2, NP)
         end if
      end if

      ! Interpolate ice concentration for current time
      if (NCICE /= 0 .and. NCICE /= 14 .and. NCICE /= 17) then
         CICE_TRatio = (TimeLoc - CICE_TIME1)/CICE_TIMINC
         do I = 1, NP
            CICEOUT(I) = CICE1(I) + CICE_TRatio*(CICE2(I) - CICE1(I))
         end do
      end if

   end subroutine updateIceFields
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update dynamic water level corrections
!>
!> Retrieves dynamic water level corrections for the current timestep
!> and computes the delta for mass adjustment if enabled.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateWaterLevelCorrections(TimeLoc)
      use global, only: usingDynamicWaterLevelCorrection, &
                        dynamicWaterLevelCorrection1, dynamicWaterLevelCorrection2, &
                        dynamicWaterLevelCorrectionDelta
      use wind, only: getDynamicWaterLevelCorrections, &
                      dynamicWaterLevelCorrectionMassAdjust
      implicit none

      real(8), intent(in) :: TimeLoc

      if (usingDynamicWaterLevelCorrection) then
         call getDynamicWaterLevelCorrections( &
            dynamicWaterLevelCorrection2, TimeLoc)
         if (dynamicWaterLevelCorrectionMassAdjust) then
            dynamicWaterLevelCorrectionDelta = &
               dynamicWaterLevelCorrection2 - &
               dynamicWaterLevelCorrection1
         end if
      end if

   end subroutine updateWaterLevelCorrections
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update wave radiation stress forcing
!>
!> Reads wave radiation stress data from file or coupled wave model and
!> adds to wind stress arrays (WSX2, WSY2). Supports multiple modes:
!> - NRS=1: Read from fort.23
!> - NRS=2: Read from RS2 format
!> - NRS=3: SWAN coupling
!> - NRS=4: Wave time extrapolation
!> - NRS=5: Linear extrapolation
!>
!> Optionally limits wave stress gradients to prevent numerical instability.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine updateWaveRadiationStress(TimeLoc)
      use global, only: NRS, RSTIME1, RSTIME2, RSTIMINC, &
                        RSNX1, RSNY1, RSNX2, RSNY2, WSX2, WSY2, &
                        RSNXOUT, RSNYOUT, ENDWAVE, Limit_WaveStressGrad, &
                        RampWRad
      use ADC_CONSTANTS, only: WaveStressGrad_Cap
      use mesh, only: NP
      use wind, only: RSGET
      use rs2, only: RS2GET
#ifdef CSWAN
      use global, only: NWS, wvnxout, wvnyout, RampMete, &
                        SWAN_RSNX1, SWAN_RSNX2, SWAN_RSNY1, SWAN_RSNY2
      use ADC_CONSTANTS, only: waveWindMultiplier
      use owiwind, only: windMultiplier
      use couple2swan, only: coupwind, swan_wx2, swan_wy2, &
                             InterpoWeight, ComputeWaveDrivenForces
#endif
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: RStRatio, RSX, RSY, RSXYMAG, RSSCALE, RSXLIM, RSYLIM
      integer :: I

#ifdef CSWAN
      ! Coupling winds to SWAN
      if (coupwind) then
         select case (abs(nws))
         case (12, 16)
            do i = 1, np
               swan_wx2(i, 2) = wvnxout(i)*waveWindMultiplier &
                                /(rampMete*windMultiplier)
               swan_wy2(i, 2) = wvnyout(i)*waveWindMultiplier &
                                /(rampMete*windMultiplier)
            end do
         case default
            do i = 1, np
               swan_wx2(i, 2) = wvnxout(i)*waveWindMultiplier/rampMete
               swan_wy2(i, 2) = wvnyout(i)*waveWindMultiplier/rampMete
            end do
         end select
      end if
#endif

      if (NRS /= 0) then
         if ((NRS == 1) .or. (nrs == 2) .or. (nrs == 3) .or. (nrs == 5)) then
            if (TimeLoc > RSTIME2) then
               RSTIME1 = RSTIME2
               RSTIME2 = RSTIME2 + RSTIMINC
               RSNX1 = RSNX2
               RSNY1 = RSNY2
               if (NRS == 1) then
                  call RSGET(RSNX2, RSNY2)
               elseif (NRS == 2) then
                  call RS2GET(RSNX2, RSNY2, NP)
               elseif (NRS == 3) then
#ifdef CSWAN
                  InterpoWeight = 1.0d0
                  call ComputeWaveDrivenForces()
                  RSNX2 = 2.0d0*SWAN_RSNX2 - SWAN_RSNX1
                  RSNY2 = 2.0d0*SWAN_RSNY2 - SWAN_RSNY1
#endif
               elseif (NRS == 5) then
                  do I = 1, NP
                     RSX = RSNX1(I)
                     RSY = RSNY1(I)
                     RSNX1(I) = RSNX2(I)
                     RSNY1(I) = RSNY2(I)
                     RSNX2(I) = RSNX2(I) + (RSNX2(I) - RSX)
                     RSNY2(I) = RSNY2(I) + (RSNY2(I) - RSY)
                  end do
               end if
            end if

            RStRatio = (TimeLoc - RSTIME1)/RSTIMINC

            if (.not. Limit_WaveStressGrad) then
               do I = 1, NP
                  RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I) - RSNX1(I)))
                  RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I) - RSNY1(I)))
                  WSX2(I) = WSX2(I) + RSX
                  WSY2(I) = WSY2(I) + RSY
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            else
               do I = 1, NP
                  RSX = RampWRad*(RSNX1(I) + RStRatio*(RSNX2(I) - RSNX1(I)))
                  RSY = RampWRad*(RSNY1(I) + RStRatio*(RSNY2(I) - RSNY1(I)))
                  RSXYMAG = sqrt(RSX*RSX + RSY*RSY)
                  if ((RSXYMAG > WaveStressGrad_Cap) .and. &
                      (RSXYMAG > 0.0d0)) then
                     RSSCALE = WaveStressGrad_Cap/RSXYMAG
                     RSXLIM = RSX*RSSCALE
                     RSYLIM = RSY*RSSCALE
                  else
                     RSXLIM = RSX
                     RSYLIM = RSY
                  end if
                  WSX2(I) = WSX2(I) + RSXLIM
                  WSY2(I) = WSY2(I) + RSYLIM
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            end if

         end if

         if (NRS == 4) then
            if (TimeLoc >= RSTIME2) then
               RSTIME1 = RSTIME2
               RSTIME2 = RSTIME2 + RSTIMINC
            end if
            if (TimeLoc > ENDWAVE + RSTIMINC) then
               RSNX2(:) = 0.0d0; RSNY2(:) = 0.0d0
            end if

            if (.not. Limit_WaveStressGrad) then
               do I = 1, NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)
                  WSX2(I) = WSX2(I) + RSX
                  WSY2(I) = WSY2(I) + RSY
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            else
               do I = 1, NP
                  RSX = RampWRad*RSNX2(I)
                  RSY = RampWRad*RSNY2(I)
                  RSXYMAG = sqrt(RSX*RSX + RSY*RSY)
                  if ((RSXYMAG > WaveStressGrad_Cap) .and. &
                      (RSXYMAG > 0.0d0)) then
                     RSSCALE = WaveStressGrad_Cap/RSXYMAG
                     RSXLIM = RSX*RSSCALE
                     RSYLIM = RSY*RSSCALE
                  else
                     RSXLIM = RSX
                     RSYLIM = RSY
                  end if
                  WSX2(I) = WSX2(I) + RSXLIM
                  WSY2(I) = WSY2(I) + RSYLIM
                  RSNXOUT(I) = RSX
                  RSNYOUT(I) = RSY
               end do
            end if
         end if
      end if

   end subroutine updateWaveRadiationStress
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute specified normal flow boundary conditions
!>
!> Computes normal flux (QN2) and elevation (EN2) boundary conditions
!> for specified flux boundaries. Supports periodic forcing with tidal
!> constituents and time-varying flux input from fort.20.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine computeSpecifiedNormalFlowBC(IT, TimeLoc, TimeH)
      use global, only: NFFR, FPER, FAMIG, FFACE, FFF, RampExtFlux, &
                        QNAM, QNPH, ENAM, ENPH, QN2, EN2, &
                        QTIME1, QTIME2, FTIMINC, QNIN1, QNIN2, ENIN1, ENIN2, &
                        FluxSettlingIT, EtaDisc, ElevDisc, Eta2
      use mesh, only: NP, ICS
      use boundaries, only: NVEL, LBCODEI, NBV
      use mod_normal_flow_boundary, only: rotate_normal_flux
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc, TimeH
      integer :: I, J, NCYC, NNBB
      real(8) :: ARGJ, RFF, QTRATIO, QNIN_TMP

      if (NFFR > 0) then
         do J = 1, NFFR
            if (abs(FPER(J)) <= eps) then
               NCYC = 0
            else
               NCYC = int(TimeH/FPER(J))
            end if
            ARGJ = FAMIG(J)*(TimeH - dble(NCYC)*FPER(J)) + FFACE(J)
            RFF = FFF(J)*RampExtFlux
            do I = 1, NVEL
               select case (LBCODEI(I))
               case (2, 12, 22, 52)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF*cos(ARGJ - QNPH(J, I)))
               case (32)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF*cos(ARGJ - QNPH(J, I)))
                  EN2(I) = EN2(I) + ENAM(J, I)*RFF*cos(ARGJ - ENPH(J, I))
               end select
            end do
         end do
      elseif ((NFFR == 0) .or. (NFFR == -1)) then
         if (TimeLoc > QTIME2) then
            QTIME1 = QTIME2
            QTIME2 = QTIME2 + FTIMINC
            do J = 1, NVEL
               if ((LBCODEI(J) == 2) .or. (LBCODEI(J) == 12) &
                   .or. (LBCODEI(J) == 22)) then
                  QNIN1(J) = QNIN2(J)
                  read (20, *) QNIN_TMP
                  QNIN2(J) = rotate_normal_flux(ICS, J, QNIN_TMP)
               elseif (LBCODEI(J) == 32) then
                  QNIN1(J) = QNIN2(J)
                  ENIN1(J) = ENIN2(J)
                  read (20, *) QNIN_TMP, ENIN2(J)
                  QNIN2(J) = rotate_normal_flux(ICS, J, QNIN_TMP)
               end if
            end do
         end if

         QTRATIO = (TimeLoc - QTIME1)/FTIMINC
         QN2 = RampExtFlux*(QNIN1 + QTRATIO*(QNIN2 - QNIN1))
         EN2 = RampExtFlux*(ENIN1 + QTRATIO*(ENIN2 - ENIN1))
      end if

      ! Collect elevation information for river radiation b.c.
      if (IT == FluxSettlingIT) then
         EtaDisc_Fill = .false.
         do I = 1, NP
            EtaDisc(I) = Eta2(I)
         end do
         do I = 1, NVEL
            if (LBCODEI(I) == 52) then
               NNBB = NBV(I)
               ElevDisc(I) = Eta2(NNBB)
            end if
         end do
      else if (EtaDisc_Fill .and. IT > FluxSettlingIT) then
         EtaDisc_Fill = .false.
         do I = 1, NVEL
            if (LBCODEI(I) == 52) then
               NNBB = NBV(I)
               ElevDisc(I) = EtaDisc(NNBB)
            end if
         end do
      end if

   end subroutine computeSpecifiedNormalFlowBC
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute radiation boundary condition discharge
!>
!> Computes discharge (QN1) at radiation boundary nodes (LBCODE=30)
!> based on previous timestep velocities and current depth.
!-----------------------------------------------------------------------
   subroutine computeRadiationDischargeBC()
      use global, only: UU1, VV1, H2, QN1
      use boundaries, only: NVEL, LBCODEI, NBV, CSII, SIII
      implicit none

      integer :: J, NNBB
      real(8) :: UN1

      do J = 1, NVEL
         if (LBCODEI(J) == 30) then
            NNBB = NBV(J)
            UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
            QN1(J) = H2(NNBB)*UN1
         end if
      end do

   end subroutine computeRadiationDischargeBC
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute discharge for zero normal velocity gradient BC
!>
!> Computes discharge (QN1) at zero normal velocity gradient boundary
!> nodes (LBCODE=40,41) based on current velocities and depth.
!-----------------------------------------------------------------------
   subroutine computeZeroNormalVelocityBCDischarge()
      use global, only: UU1, VV1, H2, QN1
      use boundaries, only: NVEL, LBCODEI, NBV, CSII, SIII
      implicit none

      integer :: J, NNBB
      real(8) :: UN1

      do J = 1, NVEL
         if ((LBCODEI(J) == 40) .or. (LBCODEI(J) == 41)) then
            NNBB = NBV(J)
            UN1 = UU1(NNBB)*CSII(J) + VV1(NNBB)*SIII(J)
            QN1(J) = H2(NNBB)*UN1
         end if
      end do

   end subroutine computeZeroNormalVelocityBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute external barrier boundary flux
!>
!> Computes supercritical outward normal flow (QN2) over external barrier
!> boundary nodes (LBCODE=3,13,23) using weir formulas.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeExternalBarrierBCDischarge(TimeLoc)
      use global, only: QN2
      use boundaries, only: NVEL, LBCODEI
      use weir_flux, only: COMPUTE_EXTERNAL_BOUNDARY_FLUX
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I

      do I = 1, NVEL
         select case (LBCODEI(I))
         case (3, 13, 23)
            call COMPUTE_EXTERNAL_BOUNDARY_FLUX(I, TimeLoc, QN2(I))
         end select
      end do

   end subroutine computeExternalBarrierBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Set submergence flags for vertical element wall boundaries
!>
!> Sets the ISSUBMERGED flag for vertical element wall (VEW) boundaries
!> (LBCODE=64) based on current water levels.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine setVewBCSubmergedFlag(TimeLoc)
      use boundaries, only: NBOU, LBCODEI, NVELL
      use WEIR_FLUX, only: SET_SUBMERGED64_AT
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I, J, K

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               call SET_SUBMERGED64_AT(I, J, K, TimeLoc)
            end do
            I = I + NVELL(K)
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine setVewBCSubmergedFlag
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute internal barrier boundary flux
!>
!> Computes inward/outward normal flow (QN2) over internal barrier
!> boundary nodes (LBCODE=4,24,5,25,64) using weir formulas. Handles
!> both permeable and impermeable barrier types.
!>
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine computeInternalBarrierBCDischarge(TimeLoc)
      use global, only: QN2, NIBNODECODE
      use boundaries, only: NBOU, LBCODEI, NVELL
      use WEIR_FLUX, only: COMPUTE_INTERNAL_BOUNDARY_FLUX, &
                           COMPUTE_INTERNAL_BOUNDARY64_FLUX
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: I, J, K
      logical :: ISFRONT

      NIBNODECODE(:) = 0
      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (4, 24, 5, 25)
            do J = 1, NVELL(K)*2
               I = I + 1
               call COMPUTE_INTERNAL_BOUNDARY_FLUX(I, J, K, &
                                                   TimeLoc, QN2(I))
            end do
         case (64)
            do J = 1, NVELL(K)*2
               I = I + 1
               if (J <= NVELL(K)) then
                  ISFRONT = .true.
               else
                  ISFRONT = .false.
               end if
               call COMPUTE_INTERNAL_BOUNDARY64_FLUX(I, J, K, &
                                                     TimeLoc, QN2(I))
            end do
         case DEFAULT
            I = I + NVELL(K)
         end select
      end do

   end subroutine computeInternalBarrierBCDischarge
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Compute cross-barrier pipe flux
!>
!> Computes inward/outward normal flow (QN2) through cross-barrier pipes
!> for internal barrier boundary nodes (LBCODE=5,25).
!-----------------------------------------------------------------------
   subroutine computeInternalBarrierPipeFlux()
      use global, only: QN2
      use boundaries, only: NVEL, LBCODEI
      use WEIR_FLUX, only: COMPUTE_CROSS_BARRIER_PIPE_FLUX
      implicit none

      integer :: I
      real(8) :: PIPE_FLUX

      do I = 1, NVEL
         if ((LBCODEI(I) == 5) .or. (LBCODEI(I) == 25)) then
            PIPE_FLUX = 0d0
            call COMPUTE_CROSS_BARRIER_PIPE_FLUX(I, PIPE_FLUX)
            QN2(I) = QN2(I) + PIPE_FLUX
         end if
      end do

   end subroutine computeInternalBarrierPipeFlux
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply subdomain boundary condition data
!>
!> Reads and applies subdomain boundary condition data from fort.019,
!> fort.020, and fort.021 files for subdomain modeling.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine applySubdomainBCData(IT)
      use subdomain, only: enforceBN, &
                           readFort019, readFort020, readFort021
      implicit none

      integer, intent(in) :: IT

      if (enforceBN == 1) call readFort019(IT)
      if (enforceBN == 2) call readFort020(IT)
      if (enforceBN == 2) call readFort021(IT)

   end subroutine applySubdomainBCData
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply all boundary conditions
!>
!> Master routine that applies all boundary conditions including normal
!> flow, radiation, barrier, VEW, baroclinic, and subdomain boundaries.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine applyBoundaryConditions(IT, TimeLoc, TimeH)
      use boundaries, only: NFLUXIB, NFLUXF, NFLUXIBP, NFLUXGBC, &
                            NFLUXB, NFLUXIB64, NFLUXRBC
      use subdomain, only: subdomainOn
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc, TimeH

      if (NFLUXF == 1) call computeSpecifiedNormalFlowBC(IT, TimeLoc, TimeH)
      if (NFLUXRBC == 1) call computeRadiationDischargeBC()
      if (NFLUXGBC == 1) call computeZeroNormalVelocityBCDischarge()
      if (NFLUXB == 1) call computeExternalBarrierBCDischarge(TimeLoc)
      if (NFLUXIB64 == 1) call setVewBCSubmergedFlag(TimeLoc)
      if (NFLUXIB == 1) call computeInternalBarrierBCDischarge(TimeLoc)
      if (NFLUXIBP == 1) call computeInternalBarrierPipeFlux()

      if (subdomainOn) call applySubdomainBCData(IT)

   end subroutine applyBoundaryConditions
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply wetting and drying algorithm
!>
!> Applies the wetting and drying algorithm to handle tidal flats and
!> inundation. Supports subgrid barrier lookup tables when enabled.
!>
!> @param[in] IT Current timestep number
!-----------------------------------------------------------------------
   subroutine applyWettingAndDrying(IT)
      use subgrid, only: subgrid_level0, getVertLookup
      use wetdry, only: computeWettingAndDrying
#ifndef ADCNETCDF
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
#endif
      implicit none

      integer, intent(in) :: IT

      if (subgrid_level0) then
#ifdef ADCNETCDF
         call getVertLookup()
         call computeWettingAndDrying(IT)
#else
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Subgrid calculation requires building " &
                        //"with netcdf enabled")
#endif
      else
         call computeWettingAndDrying(IT)
      end if

   end subroutine applyWettingAndDrying
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Update total water depth
!>
!> Updates the total water depth (H2) from bathymetry and water surface
!> elevation. Called after wetting/drying since effective bathymetry
!> could change with subgrid barriers.
!-----------------------------------------------------------------------
   subroutine updateFlowDepth()
      use global, only: H2, ETA2, IFNLFA
      use mesh, only: DP
      implicit none

      H2 = DP + dble(IFNLFA)*ETA2

   end subroutine updateFlowDepth
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply absorbing sponge layer operations
!>
!> Applies absorbing sponge layer boundary treatment after the momentum
!> solve. Uses operator splitting for the sponge relaxation terms and
!> shifts the sponge solution arrays for the next timestep.
!-----------------------------------------------------------------------
   subroutine applyAbsorbingSpongeLayer()
      use global, only: C2DDI, DTDP
      use spongelayer, only: LoadAbsLayerSigma, sponge_dis_mthd, &
                             SPONGE_OPSPLIT0, SPONGE_SHIFT_SOLN
      implicit none

      if (C2DDI .and. LoadAbsLayerSigma .and. sponge_dis_mthd == 0) &
         call SPONGE_OPSPLIT0(DTDP/2)

      if (LoadAbsLayerSigma) call SPONGE_SHIFT_SOLN()

   end subroutine applyAbsorbingSpongeLayer
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Collect output statistics and write subdomain output
!>
!> Collects min/max elevation and velocity data, inundation statistics
!> (if enabled), updates harmonic analysis, and writes subdomain output
!> files (fort.065/066/067) if subdomain modeling is enabled.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TIMEH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine collectOutputStatistics(IT, TimeLoc, TIMEH)
      use global, only: inundationOutput
      use write_output, only: collectMinMaxData, collectInundationData
      use subdomain, only: subdomainOn, NOutGS, &
                           writeFort065, writeFort066, writeFort067
      use HARM, only: updateHarmonicAnalysis
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: TIMEH

      call collectMinMaxData(TimeLoc)

      if (inundationOutput) call collectInundationData(TimeLoc, IT)

      if (subdomainOn) then
         if (NOutGS == 1) call writeFort065(IT)
         if (NOutGS == 2) call writeFort066(IT)
         if (NOutGS == 2) call writeFort067(IT)
      end if

      !     If harmonic analysis was requested and the current time
      !     is within the harmonic analysis period, update the left hand side
      !     of the harmonic analysis matrix. Also update the load vectors for
      !     each type of analysis. If timeseries reconstruction was specified,
      !     also update the timeseries.
      call updateHarmonicAnalysis(IT, TIMEH)

   end subroutine collectOutputStatistics
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Write hot start file if needed
!>
!> Writes hot start (restart) files at specified intervals for simulation
!> checkpointing. Supports both direct file writing and the HS writer for
!> parallel runs. Uses 8-byte record length for double precision data.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine writeHotstartIfNeeded(IT, TimeLoc)
      use sizes, only: MNWPROH
      use global, only: NHSINC, NHSTAR, IHOT, C3D, NHOUTONCE
      use write_output, only: writeHotstart
#ifdef CMPI
      use hswriter, only: writeHotstart_through_hswriter
#endif
#ifdef CSWAN
      use couple2swan, only: WriteSwanHotStart
#endif
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      integer :: ITEST

      ITEST = (IT/NHSINC)*NHSINC

      if ((abs(NHSTAR) > 0 .and. ITEST == IT) .or. (-IHOT == IT)) then
         if (MNWPROH > 0) then
            if (.not. C3D) then
#ifdef CMPI
               call writeHotstart_through_hswriter(TimeLoc, IT)
#endif
            else
               write (6, *) 'HS writer does not support C3D'
            end if
         else
            if (NHOUTONCE) then
               NHSINC = huge(NHSINC)
            end if
            call writeHotstart(TimeLoc, IT)
         end if
#ifdef CSWAN
         WriteSwanHotStart = .true.
#endif
      end if

   end subroutine writeHotstartIfNeeded
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Monitor and report simulation progress
!>
!> Monitors solution progress by finding and reporting maximum elevation
!> and velocity magnitudes. Checks values against warning and error
!> thresholds (WarnElev, ErrorElev, WarnVel). Detects NaN values in
!> velocity fields. Terminates simulation if error conditions are met.
!>
!> @param[in] IT        Current timestep number
!> @param[in] TimeLoc   Current simulation time (seconds)
!> @param[in] ITIME_BGN Starting timestep number (for progress percentage)
!-----------------------------------------------------------------------
   subroutine progressScreen(IT, TimeLoc)
      use global, only: ETA2, NODECODE, UU2, VV2, NT, ITHS, &
                        WarnElev, ErrorElev, Flag_ElevError, WarnElevDump
      use mod_logging, only: screenUnit, nscreen
      use mesh, only: NP
      use gwce, only: numitr
      use write_output, only: WriteWarnElev
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
#if defined(CMPI) || defined(VELCHECK)
      use sizes, only: MYPROC
      use global, only: nodes_lg
#endif
#ifdef VELCHECK
      use global, only: earlyterminate, IT_foundNaN, WarnVel
      use mesh, only: DP, NNeigh, NeiTab
#endif
#ifdef CMPI
      use messenger, only: WarnElevSum, EarlyTermSum
#endif
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc

      integer :: I
      integer :: KEMAX, KVMAX
      integer :: ITEST
#ifdef CMPI
      integer :: WarnElevExceeded, ErrorElevExceeded
#endif
#ifdef VELCHECK
      integer :: J, JJ
      integer :: WarnVelExceeded
      logical :: foundNaN
#endif
      real(8) :: ELMAX, VELMAX, velabs
      real(8) :: runperccomplete

      if (NSCREEN /= 0) then
         ELMAX = 0.0d0
         VELMAX = 0.0d0
         KEMAX = 0
         KVMAX = 0
         do I = 1, NP
            if ((NODECODE(I) == 1) .and. (abs(ETA2(I)) > ELMAX)) then
               ELMAX = abs(ETA2(I))
               KEMAX = I
            end if
            VELABS = UU2(I)*UU2(I) + VV2(I)*VV2(I)
            if (VELABS > VELMAX) then
               VELMAX = VELABS
               KVMAX = I
            end if
         end do
         VELMAX = VELMAX**0.5d0
         ITEST = (IT/NSCREEN)*NSCREEN

         runperccomplete = real(IT - ITHS, 8)/real(NT - ITHS, 8)*100d0
#ifdef CMPI
         if (MYPROC == 0 .and. ELMAX < WarnElev .and. ITEST == IT) then
            if (KEMAX > 0) then
               write (ScreenUnit, 1991) &
                  IT, runperccomplete, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
                  VELMAX, KVMAX, MYPROC
            else
               write (ScreenUnit, 1991) &
                  IT, runperccomplete, NUMITR, TimeLoc, 0., KEMAX, VELMAX, &
                  KVMAX, MYPROC
            end if
1991        format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', &
                   5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I8, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I8, &
                   2x, 'ON MYPROC = ', I4)
         end if
         WarnElevExceeded = 0
         if (ELMAX > WarnElev) then
            write (ScreenUnit, 1993) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 1993) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
1993        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, &
                   3x, '** WARNING: Elevation.gt.WarnElev **')
            if (WarnElevDump) WarnElevExceeded = 1
         end if

         if (WarnElevDump) then
            call WarnElevSum(WarnElevExceeded)
            if (WarnElevExceeded /= 0) then
               call WriteWarnElev(TimeLoc, IT)
            end if
         end if
         ErrorElevExceeded = 0
         if (ELMAX > ErrorElev) then
            write (ScreenUnit, 1995) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 1995) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
1995        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, /, &
                   2x, '** ERROR: Elevation.gt.ErrorElev,', &
                   ' ADCIRC stopping. **')
            ErrorElevExceeded = 1
         end if
         call WarnElevSum(ErrorElevExceeded)
         if (ErrorElevExceeded /= 0) then
            Flag_ElevError = .true.
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Elevation exceeded ErrorElev threshold')
         end if

#ifdef VELCHECK
         WarnVelExceeded = 0
         if (VELMAX > WarnVel) then
            write (ScreenUnit, 3995) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 3995) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
3995        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, /, &
                   2x, '** WARNING: Velocity.gt.WarnVel **')
            WarnVelExceeded = 1
         end if
#endif

#else
         if (ELMAX < WarnElev .and. ITEST == IT) then
            if (KEMAX > 0) then
               write (ScreenUnit, 1992) &
                  IT, runperccomplete, NUMITR, TimeLoc, ETA2(KEMAX), &
                  KEMAX, VELMAX, KVMAX
            else
               write (ScreenUnit, 1992) &
                  IT, runperccomplete, NUMITR, TimeLoc, 0., KEMAX, &
                  VELMAX, KVMAX
            end if
1992        format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', &
                   5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9)
         end if
         if (ELMAX > WarnElev) then
            write (ScreenUnit, 1994) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 1994) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
1994        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, '** WARNING: Elevation.gt.WarnElev **')
            if (WarnElevDump) call WriteWarnElev(TimeLoc, IT)
         end if
         if (ELMAX > ErrorElev) then
            write (6, *) 'elmax, errorelev', elmax, errorelev
            Flag_ElevError = .true.
            write (ScreenUnit, 1996) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 1996) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
1996        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, /, &
                   2x, '** ERROR: Elevation.gt.ErrorElev, ', &
                   'ADCIRC stopping. **')
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Elevation exceeded ErrorElev threshold')
         end if
#ifdef VELCHECK
         if (VELMAX > WarnVel) then
            write (ScreenUnit, 3996) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, VELMAX, KVMAX
            write (16, 3996) IT, NUMITR, TimeLoc, ETA2(KEMAX), KEMAX, &
               VELMAX, KVMAX
3996        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, /, &
                   2x, '** WARNING: Velocity.gt.WarnVel **')
         end if
#endif

#endif

      end if

#ifdef VELCHECK
      foundNaN = .false.
      do I = 1, NP
         if (ISNAN(UU2(I)) .or. ISNAN(VV2(I))) then
            foundNaN = .true.
            if (NSCREEN /= 0) then
               write (ScreenUnit, 2993) &
                  IT, NUMITR, TimeLoc, NODES_LG(I), &
                  ETA2(I), DP(I), MYPROC
               write (16, 2993) IT, NUMITR, TimeLoc, &
                  NODES_LG(I), ETA2(I), DP(I), MYPROC
2993           format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                      5x, 'TIME = ', E15.8, &
                      2x, 'SPEED IS NaN AT NODE I =', I9, &
                      /, 2x, 'ETA2(I) = ', 1pe12.4e3, &
                      2x, 'DP(I) = ', 1pe12.4e3, &
                      2x, 'ON MYPROC = ', I4)
               do J = 2, NNeigh(I)
                  JJ = NeiTab(I, J)
                  write (ScreenUnit, 2995) NODES_LG(JJ), ETA2(JJ), DP(JJ), &
                     UU2(JJ), VV2(JJ), MYPROC
                  write (16, 2995) NODES_LG(JJ), ETA2(JJ), DP(JJ), &
                     UU2(JJ), VV2(JJ), MYPROC
2995              format(1x, 'AT NEIGHBOR NODE I =', I9, &
                         /, 2x, 'ETA2(I) = ', 1pe12.4e3, &
                         2x, 'DP(I) = ', 1pe12.4e3, &
                         2x, 'UU2(I) = ', 1pe12.4e3, &
                         2x, 'VV2(I) = ', 1pe12.4e3, &
                         2x, 'ON MYPROC = ', I4)
               end do
            end if
         end if
      end do
      if (foundNaN) then
         if (NSCREEN /= 0) then
            write (ScreenUnit, 2994) &
               IT, NUMITR, TimeLoc, ETA2(KEMAX), NODES_LG(KEMAX), &
               VELMAX, NODES_LG(KVMAX), MYPROC
            write (16, 2994) IT, NUMITR, TimeLoc, ETA2(KEMAX), &
               NODES_LG(KEMAX), VELMAX, NODES_LG(KVMAX), MYPROC
2994        format(1x, 'TIME STEP =', I8, 5x, 'ITERATIONS =', I5, &
                   5x, 'TIME = ', E15.8, &
                   /, 2x, 'ELMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I9, &
                   2x, 'ON MYPROC = ', I4, &
                   3x, '** WARNING: Velocity is NaN **')
         end if
      end if
      if ((foundNaN) .and. (IT_foundNaN == 0)) then
         IT_foundNaN = IT
      elseif ((IT == IT_foundNaN + 20) .and. (IT_foundNaN /= 0)) then
         earlyterminate = 1
      end if
#ifdef CMPI
      call EarlyTermSum(earlyterminate)
#endif
#endif

   end subroutine progressScreen
   !----------------------------------------------------------------------

end module mod_timestep
