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
!  MODULE BAROCLINIC
!-----------------------------------------------------------------------
!> @brief Baroclinic forcing for the ADCIRC timestep
!>
!> This module computes the baroclinic pressure gradient terms used by the
!> GWCE and momentum equations, and manages the level of no motion (LNM)
!> elevation boundary condition for baroclinic runs.
!>
!> The baroclinic pressure gradient supports:
!> - C2DDI: depth-averaged density (BPG2D), optionally coupled to 3D
!>   baroclinic boundary information (IDEN >= 5)
!> - C3D: 3D baroclinic pressure field (BPG3D) with trapezoidal, Simpson,
!>   or Gaussian quadrature vertical integration and linear or cubic
!>   vertical interpolation
!>
!> The LNM boundary condition (fort.35, BCFLAG_LNM) supports:
!> - Prognostic runs (RES_BC_FLAG > 0): LNM_BC is linearly interpolated in
!>   time between two fort.35 snapshots, including catch-up on hotstart
!> - Diagnostic runs (RES_BC_FLAG < 0): LNM_BC is read once and held fixed
!-----------------------------------------------------------------------
module mod_baroclinic
   implicit none

   private

   public :: computeBaroclinicPressureGradient, updateBaroclinicBoundaryConditions, &
             allocateLNMBoundaryConditions, initializeLNMBoundaryConditions, &
             hotstartLNMBoundaryConditions, initializeDiagnosticLNMBoundaryConditions

   !> Level of no motion (LNM) elevation adjustment at each elevation
   !> specified boundary node, applied to the boundary elevation in the GWCE
   real(8), allocatable, public, protected :: LNM_BC(:)

   ! The two fort.35 LNM snapshots bracketing the current time and the
   ! times (seconds) they are valid at; LNM_BC is interpolated between them
   real(8), allocatable :: LNM_BC1(:), LNM_BC2(:)
   real(8) :: RBCTIME1 = 0.0d0
   real(8) :: RBCTIME2 = 0.0d0

   integer, parameter :: LNM_UNIT = 35

contains

   !-----------------------------------------------------------------------
   !> @brief Allocate the level of no motion boundary condition arrays
   !>
   !> @param[in] n_boundary_nodes Number of elevation specified boundary nodes (MNETA)
   !-----------------------------------------------------------------------
   subroutine allocateLNMBoundaryConditions(n_boundary_nodes)
      implicit none
      integer, intent(in) :: n_boundary_nodes

      if (allocated(LNM_BC)) return
      allocate (LNM_BC(n_boundary_nodes), LNM_BC1(n_boundary_nodes), LNM_BC2(n_boundary_nodes))
      LNM_BC = 0.0d0
      LNM_BC1 = 0.0d0
      LNM_BC2 = 0.0d0

   end subroutine allocateLNMBoundaryConditions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Initialize prognostic level of no motion boundary conditions
   !>
   !> BCFLAG_LNM = 1: opens fort.35 and reads the first two snapshots
   !> BCFLAG_LNM = 3: holds the boundary at the initial elevation
   !-----------------------------------------------------------------------
   subroutine initializeLNMBoundaryConditions()
      use global, only: BCFLAG_LNM, BCSTATIM, RBCTIMEINC, ETA2
      use boundaries, only: NETA, NBD
      implicit none

      if (BCFLAG_LNM == 1) then
         call openLNMFile()
         RBCTIME1 = BCSTATIM*86400.d0
         RBCTIME2 = RBCTIME1 + RBCTIMEINC
         call readLNMSnapshot(LNM_BC1)
         call readLNMSnapshot(LNM_BC2)
      else if (BCFLAG_LNM == 3) then
         LNM_BC1(1:NETA) = ETA2(NBD(1:NETA))
         LNM_BC2(1:NETA) = LNM_BC1(1:NETA)
      end if
      ! BCFLAG_LNM = 2 (HYCOM) is not yet implemented

   end subroutine initializeLNMBoundaryConditions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Initialize prognostic LNM boundary conditions on hotstart
   !>
   !> Initializes as on a cold start, then advances through fort.35 so the
   !> bracketing snapshots match the hotstart time.
   !>
   !> @param[in] ITHS    Timestep number of the hotstart
   !> @param[in] TimeLoc Simulation time of the hotstart (seconds)
   !-----------------------------------------------------------------------
   subroutine hotstartLNMBoundaryConditions(ITHS, TimeLoc)
      use global, only: BCFLAG_LNM, DTDPHS, STATIM
      implicit none
      integer, intent(in) :: ITHS
      real(8), intent(in) :: TimeLoc
      integer :: IT

      call initializeLNMBoundaryConditions()

      if (BCFLAG_LNM == 1) then
         do IT = 1, ITHS
            call advanceLNMSnapshots(dble(IT)*DTDPHS + STATIM*86400.d0)
         end do
         call advanceLNMSnapshots(TimeLoc)
      end if

   end subroutine hotstartLNMBoundaryConditions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Initialize diagnostic (RES_BC_FLAG < 0) LNM boundary conditions
   !>
   !> The diagnostic LNM is constant in time and is set directly in LNM_BC.
   !-----------------------------------------------------------------------
   subroutine initializeDiagnosticLNMBoundaryConditions()
      use global, only: BCFLAG_LNM, ETA2
      use boundaries, only: NETA, NBD
      implicit none

      if (BCFLAG_LNM == 1) then
         call openLNMFile()
         call readLNMSnapshot(LNM_BC)
      else if (BCFLAG_LNM == 3) then
         LNM_BC(1:NETA) = ETA2(NBD(1:NETA))
      end if
      ! BCFLAG_LNM = 2 (HYCOM) is not yet implemented

   end subroutine initializeDiagnosticLNMBoundaryConditions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Open the fort.35 level of no motion file for reading
   !>
   !> Terminates the run with an error message if the file is missing or
   !> cannot be opened.
   !-----------------------------------------------------------------------
   subroutine openLNMFile()
      use sizes, only: INPUTDIR
      use mod_io, only: openFileForRead
      implicit none
      integer :: errorIO

      call openFileForRead(LNM_UNIT, trim(INPUTDIR)//'/'//'fort.35', errorIO)

   end subroutine openLNMFile
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Read one LNM snapshot (header line + one value per boundary node)
   !-----------------------------------------------------------------------
   subroutine readLNMSnapshot(values)
      use boundaries, only: NETA
      implicit none
      real(8), intent(inout) :: values(:)
      integer :: I, NOD
      character(80) :: header

      read (LNM_UNIT, '(A)') header
      do I = 1, NETA
         read (LNM_UNIT, *) NOD, values(I)
      end do

   end subroutine readLNMSnapshot
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Shift to the next fort.35 snapshot if TimeLoc is past RBCTIME2
   !-----------------------------------------------------------------------
   subroutine advanceLNMSnapshots(TimeLoc)
      use global, only: RBCTIMEINC
      implicit none
      real(8), intent(in) :: TimeLoc

      if (TimeLoc > RBCTIME2) then
         RBCTIME1 = RBCTIME2
         RBCTIME2 = RBCTIME2 + RBCTIMEINC
         LNM_BC1 = LNM_BC2
         call readLNMSnapshot(LNM_BC2)
      end if

   end subroutine advanceLNMSnapshots

   !-----------------------------------------------------------------------
   !> @brief Compute baroclinic pressure gradient terms
   !>
   !> Computes the depth-averaged baroclinic forcing (Bx/H, By/H) for use
   !> in the GWCE and 2D momentum equations. Calls BPG3D for 3D velocity
   !> form or BPG2D for 2D depth-integrated runs.
   !>
   !> @param[in] IT      Current timestep number
   !> @param[in] TimeLoc Current simulation time (seconds)
   !-----------------------------------------------------------------------
   subroutine computeBaroclinicPressureGradient(TimeLoc)
      use global, only: CBaroclinic, C3DVS, C2DDI
      implicit none

      real(8), intent(in) :: TimeLoc

      if (CBaroclinic) then
         if (C3DVS) then
            call BPG3D()
         else if (C2DDI) then
            call BPG2D(TimeLoc)
         end if
      end if

   end subroutine computeBaroclinicPressureGradient
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Update baroclinic boundary conditions
   !>
   !> Updates elevation boundary conditions for 3D baroclinic runs by reading
   !> large-scale mean (LNM) boundary conditions from fort.35 input file.
   !> LNM values are linearly interpolated between input snapshots.
   !>
   !> @param[in] TimeLoc Current simulation time (seconds)
   !-----------------------------------------------------------------------
   subroutine updateBaroclinicBoundaryConditions(TimeLoc)
      use global, only: C3D, RES_BC_FLAG, CBAROCLINIC, &
                        BCFLAG_LNM, RBCTIMEINC
      use boundaries, only: NOPE, NETA
      implicit none

      real(8), intent(in) :: TimeLoc
      real(8) :: RBCRATIO

      if ((C3D) .and. (RES_BC_FLAG > 0) .and. (CBAROCLINIC)) then
         if ((abs(RES_BC_FLAG) >= 1) .and. (NOPE > 0)) then
            if (BCFLAG_LNM == 1) then
               call advanceLNMSnapshots(TimeLoc)
               RBCRATIO = (TimeLoc - RBCTIME1)/RBCTIMEINC
               LNM_BC(1:NETA) = LNM_BC1(1:NETA) + &
                                RBCRATIO*(LNM_BC2(1:NETA) - LNM_BC1(1:NETA))
            end if
         end if
      end if

   end subroutine updateBaroclinicBoundaryConditions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Compute 2D baroclinic pressure gradient
   !>
   !> Computes the vertically integrated baroclinic pressure gradient divided
   !> by depth (VIDBCPDXOH, VIDBCPDYOH) for use in the 2D momentum equations.
   !> Supports both the depth-averaged equation and coupling to 3D baroclinic
   !> models. When coupled to 3D, momentum dispersion is also estimated.
   !>
   !> @param[in] IT_unused Unused timestep parameter (for interface compatibility)
   !> @param[in] TimeLoc   Current simulation time (seconds)
   !-----------------------------------------------------------------------
   subroutine BPG2D(TimeLoc)

      use ADC_CONSTANTS, only: G, RhoWat0, SigT0
      use GLOBAL, only: IFNLFA, IDEN, RAMP, NODECODE, NOFF, &
                        VIDBCPDXOH, VIDBCPDYOH, ETA2, DASigT, &
                        IFSFM, H2
      use MESH, only: NE, NM, FDXE, FDYE, NP, AREAS, TOTALAREA, &
                      SFacEle, SFMYEle, SFMXEle
      use Couple2BC3D, only: Update_BC3D_Info, FBPG_Disp_from_BC3D
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: IE !element loop counter
      integer :: J !node loop counter
      integer :: NCELE !element code
      integer :: NM1, NM2, NM3 !local node numbers used to compute gradients
      integer :: NC1, NC2, NC3 !local node codes

      real(8) :: SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac
      real(8) :: AreaIE2, AreaEle
      real(8) :: H2N1, H2N2, H2N3
      real(8) :: FDX1, FDX2, FDX3
      real(8) :: FDY1, FDY2, FDY3
      real(8) :: FDX1O2A, FDX2O2A, FDX3O2A
      real(8) :: FDY1O2A, FDY2O2A, FDY3O2A
      real(8) :: EtaN1, EtaN2, EtaN3
      real(8) :: DEta2DX, DEta2DY
      real(8) :: DRhoDX, DRhoDY
      real(8) :: DARhoMRho0N1, DARhoMRho0N2, DARhoMRho0N3
      real(8) :: VIDBCPDXOHN1, VIDBCPDXOHN2, VIDBCPDXOHN3
      real(8) :: VIDBCPDYOHN1, VIDBCPDYOHN2, VIDBCPDYOHN3
      real(8) :: VIDBCPDXOHAvgArea, VIDBCPDYOHAvgArea

      real(8), parameter :: eps = epsilon(1.0d0)

      ! For the Coupling to 3D Baroclinic mode
      if (abs(IDEN) >= 5) then
         ! Update the 3D BC Information for the current timestep
         call Update_BC3D_Info(IDEN, TimeLoc)

         if (abs(IDEN) /= 7 .and. abs(IDEN) /= 9) then
            ! Calculating the full Baroclinic pressure gradient and
            ! momentum dispersion terms required at every timestep
            ! from the 3D BC information that is on a longer time scale
            call FBPG_Disp_from_BC3D()
         end if
         return
      end if

      ! For the depth-averaged baroclinic mode without coupling
      do IE = 1, NE
         NM1 = NM(IE, 1)
         NM2 = NM(IE, 2)
         NM3 = NM(IE, 3)
         NC1 = NODECODE(NM1)
         NC2 = NODECODE(NM2)
         NC3 = NODECODE(NM3)
         NCEle = NC1*NC2*NC3*NOFF(IE)
         H2N1 = H2(NM1)
         H2N2 = H2(NM2)
         H2N3 = H2(NM3)
         EtaN1 = dble(IFNLFA)*Eta2(NM1)
         EtaN2 = dble(IFNLFA)*Eta2(NM2)
         EtaN3 = dble(IFNLFA)*Eta2(NM3)
         SFacAvg = SFacEle(IE)
         !..... BEG DW/WJP
         SFmxAvg = SFMXEle(IE); 
         SFmyAvg = SFMYEle(IE); 
         sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg; 
         sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg; 
         !..... END DW/WJP

         AreaIE2 = Areas(IE)
         AreaEle = dble(NCEle)*AreaIE2*0.5d0
         FDX1 = FDXE(1, IE)*sfdxfac; 
         !c FDX1=(Y(NM2)-Y(NM3))*MX !b1*mx
         FDX2 = FDXE(2, IE)*sfdxfac; 
         !c FDX2=(Y(NM3)-Y(NM1))*MX !b2*mx
         FDX3 = FDXE(3, IE)*sfdxfac; 
         !c FDX3=(Y(NM1)-Y(NM2))*MX !b3*mx
         FDY1 = FDYE(1, IE)*sfdyfac; 
         !c FDY1=(X(NM3)-X(NM2))*MY !a1*my
         FDY2 = FDYE(2, IE)*sfdyfac; 
         !c FDY2=(X(NM1)-X(NM3))*MY !a2*my
         FDY3 = FDYE(3, IE)*sfdyfac; 
         !c FDY3=(X(NM2)-X(NM1))*MY !a3*my
         FDX1O2A = FDX1/AreaIE2 !dphi1/dx
         FDY1O2A = FDY1/AreaIE2 !dphi1/dy
         FDX2O2A = FDX2/AreaIE2 !dphi2/dx
         FDY2O2A = FDY2/AreaIE2 !dphi2/dy
         FDX3O2A = FDX3/AreaIE2 !dphi3/dx
         FDY3O2A = FDY3/AreaIE2 !dphi3/dy

         DARhoMRho0N1 = (DASigT(NM1) - SigT0)/RhoWat0
         DARhoMRho0N2 = (DASigT(NM2) - SigT0)/RhoWat0
         DARhoMRho0N3 = (DASigT(NM3) - SigT0)/RhoWat0
         DEta2DX = EtaN1*FDX1O2A + EtaN2*FDX2O2A + EtaN3*FDX3O2A
         DEta2DY = EtaN1*FDY1O2A + EtaN2*FDY2O2A + EtaN3*FDY3O2A
         DRhoDX = DARhoMRho0N1*FDX1O2A + DARhoMRho0N2*FDX2O2A &
                  + DARhoMRho0N3*FDX3O2A
         DRhoDY = DARhoMRho0N1*FDY1O2A + DARhoMRho0N2*FDY2O2A &
                  + DARhoMRho0N3*FDY3O2A
         VIDBCPDXOHN1 = Ramp*G* &
                        (DARhoMRho0N1*DEta2DX + 0.5d0*H2N1*DRhoDX)
         VIDBCPDXOHN2 = Ramp*G* &
                        (DARhoMRho0N2*DEta2DX + 0.5d0*H2N2*DRhoDX)
         VIDBCPDXOHN3 = Ramp*G* &
                        (DARhoMRho0N3*DEta2DX + 0.5d0*H2N3*DRhoDX)
         VIDBCPDYOHN1 = Ramp*G* &
                        (DARhoMRho0N1*DEta2DY + 0.5d0*H2N1*DRhoDY)
         VIDBCPDYOHN2 = Ramp*G* &
                        (DARhoMRho0N2*DEta2DY + 0.5d0*H2N2*DRhoDY)
         VIDBCPDYOHN3 = Ramp*G* &
                        (DARhoMRho0N3*DEta2DY + 0.5d0*H2N3*DRhoDY)
         VIDBCPDXOHAvgArea = AreaEle*(VIDBCPDXOHN1 + VIDBCPDXOHN2 &
                                      + VIDBCPDXOHN3)/3.d0
         VIDBCPDYOHAvgArea = AreaEle*(VIDBCPDYOHN1 + VIDBCPDYOHN2 &
                                      + VIDBCPDYOHN3)/3.d0
         VIDBCPDXOH(NM1) = VIDBCPDXOH(NM1) + VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM2) = VIDBCPDXOH(NM2) + VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM3) = VIDBCPDXOH(NM3) + VIDBCPDXOHAvgArea
         VIDBCPDYOH(NM1) = VIDBCPDYOH(NM1) + VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM2) = VIDBCPDYOH(NM2) + VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM3) = VIDBCPDYOH(NM3) + VIDBCPDYOHAvgArea
      end do

      do J = 1, NP
         if (abs(TotalArea(J)) > eps) then
            VIDBCPDXOH(J) = VIDBCPDXOH(J)/TotalArea(J)
            VIDBCPDYOH(J) = VIDBCPDYOH(J)/TotalArea(J)
         end if
      end do

   end subroutine BPG2D
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Compute 3D baroclinic pressure gradient
   !>
   !> Computes the 3D baroclinic pressure field (BCP) and its horizontal
   !> gradient (BPG) at each vertical level. Also computes the vertically
   !> integrated baroclinic pressure gradient for the wave equation.
   !>
   !> Supports multiple density integration schemes:
   !> - Trapezoidal: Standard trapezoidal rule integration
   !> - Simpson: Higher-order Simpson integration with quadratic interpolation
   !> - Quadrature: Gaussian quadrature integration
   !>
   !> The baroclinic pressure is computed from the density field as:
   !> BCP(sigma) = (g/rho_ref) * (H/(a-b)) * integral(SigT) from surface to sigma
   !>
   !> Time stepping coefficients (IDTAlp1, DTAlp3, etc.) are computed in VSSTUP.
   !-----------------------------------------------------------------------
   subroutine BPG3D()

      !Casey: Added the following variable declarations from GLOBAL.
      use GLOBAL, only: IFNLFA, &
                        NODECODE, NOFF, VIDBCPDXOH, VIDBCPDYOH, &
                        ETA2, H0, H2
      use mod_logging, only: t_log_scope, init_log_scope
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use ADC_CONSTANTS, only: SigT0
      use MESH, only: NM, X, Y, NP, AREAS, &
                      NEITABELE, SFAC, DP
      use GLOBAL_3DVS, only: BCP, BPG, NFEN, SIGT, SIGMA, IDEN, &
                             GORhoOAMB, AMB, B, IY, &
                             densityIntegrationID, &
                             verticalInterpolationSchemeID, &
                             NumIP, IIP, NumGQP, GQP, GQW, SIGTM, &
                             DENSITY_INTEGRATION_TRAPEZOIDAL, DENSITY_INTEGRATION_SIMPSON, &
                             VERTICAL_INTERPOLATION_LINEAR

#ifdef CMPI
      use MESSENGER, only: UPDATER, UPDATEC3D
#endif

      implicit none

#ifdef CMPI
      real(8) :: DUMV1(1)
#endif

      !asey: Added the following local variable declarations.
      integer :: NCELE

      integer :: NEle !local value of NetTabEle
      integer :: k !vertical node loop counter (1-bottom, NFEN-surf)
      integer :: NH !horizontal node loop counter
      integer :: N !neighbor node loop counter
      integer :: N1, N2, N3 !local node numbers used to compute gradients
      integer :: LocalVertex !position of NH within the neighbor element

      real(8) :: Hs !Total water depth at time level s
      real(8) :: HsOAMB !Hs/(a-b)

      real(8) :: Zk !z depth of any node k in the vertical
      real(8) :: DelSig ! sigma(k+1)-sigma(k)
      real(8) :: SigmaNN !Sigma value of a neighbor node
      real(8) :: HsN2 !Depth value of a neighbor node

      real(8) :: SFacAvg ! kmd48.33bc add in spherical factors

      real(8) :: BCPN1, BCPN2, BCPN3 !nodal values of BCP
      real(8) :: SigTAvg !avg SigT between 2 vertical nodes
      real(8) :: HGORhoOAMB !depth*gravity/(reference density)/(a-b)

      real(8) :: a1, a2, a3, b1, b2, b3
      real(8) :: TotalBCPGArea2

      real(8) :: DBCPDX2A
      real(8) :: DBCPDY2A
      complex(8) :: BCPG(NFEN) !baroclinic pressure gradient
      complex(8) :: VIBCPG !baroclinic pressure gradient

      real(8) :: SigT_mid_point, SigmaMP, f1, f2, f3
      complex(8) :: BCPG_mid_point

      !Casey 230120 : Add to compute SigT at higher resolution.
      integer :: IGQ, GQI
      real(8) :: halfdist, midpoint
      integer :: IIP1, IIP2
      real(8) :: w(NumIP)
      complex(8) :: BCPGM(NumGQP*(NFEN - 1))
      character(len=256) :: bpg_error_msg

      real(8), parameter :: one_sixth = 1.0d0/6.0d0

      LOG_SCOPE_TRACED("bpg3d", TIMESTEP_TRACING)

      !     INCREMENT THE TIMESTEP SINCE START COUNTER
      !

      !*************************************************************************************
      !     Check whether it is time to print various 3D outputs

      !
      !     If a baroclinic run, compute the 3D baroclinic pressure field
      !     The buoyancy field is defined as
      !     BCP(z)    =(gravity/rho ref)*          integral (SigT) from surface down to z
      !     BCP(sigma)=(gravity/rho ref)*(H/(a-b))*integral (SigT) from a down to sigma
      !     where
      !     SigT = Sigma T = Rho - 1000 = density - 1000
      !     SigT0 = Sigma t value of reference density (typically = 0)
      !     Sigma = dimensionless vertical coordinate
      !
      if ((IDEN >= 1) .or. (IDEN <= -1)) then
         do NH = 1, NP !loop over horizontal nodes

            Hs = DP(NH) + dble(IFNLFA)*Eta2(NH) !total depth at previous (s) timestep
            HGORhoOAMB = GORhoOAMB*Hs !(gravity/rho ref)*(H/(a-b))

            if (densityIntegrationID == DENSITY_INTEGRATION_TRAPEZOIDAL) then
               BCP(NH, NFEN) = 0.d0
               do k = NFEN - 1, 1, -1 !loop over vertical nodes, starting at top and working down
                  SigTAvg = (SigT(NH, k + 1) + SigT(NH, k))/2.d0
                  DelSig = Sigma(k + 1) - Sigma(k)
                  BCP(NH, k) = BCP(NH, k + 1) + HGORhoOAMB*(SigTAvg - SigT0)*DelSig
               end do
            elseif (densityIntegrationID == DENSITY_INTEGRATION_SIMPSON) then
               !--------------------------------------------------------------------------------------------------
               ! arash 21/1/2016: Simpson integration, using quadratic function interpolation.
               ! arash 25/2/2016: fixed a bug.
               !------------------------------------------------------------------------------
               ! 1. initialize
               BCP(NH, NFEN) = 0.d0

               ! 2. to compute the integral at node NFEN-1, we need the mid-point value, which we evaluate by using a quadratic
               !    interpolation, using nodes NFEN, NFEN-1, and NFEN-2.
               !              SigmaMP = ( Sigma(NFEN) + Sigma(NFEN-1) ) / 2.0d0
               !
               !              f1 = ( (SigmaMP       - Sigma(NFEN-1)) * (SigmaMP       - Sigma(NFEN-2)) )
               !     &           / ( (Sigma(NFEN)   - Sigma(NFEN-1)) * (Sigma(NFEN)   - Sigma(NFEN-2)) )
               !              f2 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-2)) )
               !     &           / ( (Sigma(NFEN-1) - Sigma(NFEN))   * (Sigma(NFEN-1) - Sigma(NFEN-2)) )
               !              f3 = ( (SigmaMP       - Sigma(NFEN))   * (SigmaMP       - Sigma(NFEN-1)) )
               !     &           / ( (Sigma(NFEN-2) - Sigma(NFEN))   * (Sigma(NFEN-2) - Sigma(NFEN-1)) )
               !
               !              SigT_mid_point = f1 * SigT(NH,NFEN) + f2 * SigT(NH,NFEN-1) + f3 * SigT(NH,NFEN-2)
               SigT_mid_point = SIGTM(NH, NFEN - 1)

               ! 3. evaluate the integral at node NFEN-1
               BCP(NH, NFEN - 1) = one_sixth*(Sigma(NFEN) - Sigma(NFEN - 1))*HGORhoOAMB* &
                                   ((SigT(NH, NFEN) - SigT0) + 4.0d0*(SigT_mid_point - SigT0) + (SigT(NH, NFEN - 1) - SigT0))

               ! 4. evaluate the integral at the rest of the vertical nodes
               do k = NFEN - 2, 1, -1

                  ! 4.1.
                  !                 SigmaMP = ( Sigma(k+1) + Sigma(k) ) / 2.0d0
                  !
                  !                 f1 = ( (SigmaMP    - Sigma(k+1)) * (SigmaMP    - Sigma(k)) )
                  !     &              / ( (Sigma(k+2) - Sigma(k+1)) * (Sigma(k+2) - Sigma(k)) )
                  !                 f2 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k)) )
                  !     &              / ( (Sigma(k+1) - Sigma(k+2)) * (Sigma(k+1) - Sigma(k)) )
                  !                 f3 = ( (SigmaMP    - Sigma(k+2)) * (SigmaMP    - Sigma(k+1)) )
                  !     &              / ( (Sigma(k)   - Sigma(k+2)) * (Sigma(k)   - Sigma(k+1)) )
                  !
                  !                 SigT_mid_point = f1 * SigT(NH,k+2) + f2 * SigT(NH,k+1) + f3 * SigT(NH,k)
                  SigT_mid_point = SIGTM(NH, k)

                  ! 4.2. evaluate the integral
                  BCP(NH, k) = BCP(NH, k + 1) + one_sixth*(Sigma(k + 1) - Sigma(k))*HGORhoOAMB &
                               *((SigT(NH, k + 1) - SigT0) + 4.0d0*(SigT_mid_point - SigT0) + (SigT(NH, k) - SigT0))

               end do

               ! 5. let's hope we are better of now.

               !Casey 230106 : Try Arash's higher-order integration, but with SigT actually
               !               computed at the quadrature points of each vertical cell.
            else ! default is densityIntegration = 'quadrature'
               BCP(NH, NFEN) = 0.d0
               do k = NFEN - 1, 1, -1
                  BCP(NH, k) = 0.d0
                  do IGQ = NumGQP, 1, -1
                     GQI = NumGQP*(k - 1) + IGQ
                     BCP(NH, k) = BCP(NH, k) + (Sigma(k + 1) - Sigma(k))*HGORhooAMB &
                                  *GQW(IGQ)*SIGTM(NH, GQI)
                  end do
                  BCP(NH, k) = BCP(NH, k + 1) + BCP(NH, k)/sum(GQW)
               end do
            end if
            !--------------------------------------------------------------------------------------------------
         end do
#ifdef CMPI
!     Update BCP on ghost nodes
!      CALL UPDATER3D(BCP) !!!!Don't know if this is needed at this time
#endif
      end if

      !*************************************************************************************
      !     Compute 3D baroclinic pressure gradients

      !
      !     Loop over each horizontal node to compute the horizontal velocity
      !
      do NH = 1, NP !loop over horizontal nodes

         HsOAMB = H2(NH)/AMB

         !     Zero out baroclinic pressure gradient and vertically integrated
         !     baroclinic pressure gradient for a barotropic run

         do k = 1, NFEN
            BCPG(k) = (0.d0, 0.d0)
         end do
         VIDBCPDXOH(NH) = 0.d0
         VIDBCPDYOH(NH) = 0.d0

         !     Start computing baroclinic terms

         if ((IDEN >= 1) .or. (IDEN <= -1)) then

            !     Start computing baroclinic pressure gradient (computed in level
            !     coordinates) at each node in the vertical

            do k = 1, NFEN
               ! Start from zero
               DBCPDX2A = 0.d0
               DBCPDY2A = 0.d0
               TotalBCPGArea2 = 0.d0

               BCPG(k) = (0.d0, 0.d0)
               do N = 1, size(NeiTabEle, 2)
                  ! Find neighbor element
                  NEle = NeiTabEle(NH, N)
                  ! Skip blank entries in NeiTabEle
                  if (NEle == 0) cycle
                  ! Find vertices in counterclockwise order, starting at NH
                  LocalVertex = findloc(NM(NEle, :), NH, dim=1)
                  if (LocalVertex == 0) then
                     write (bpg_error_msg, '(A,I0,A,I0,A,3(1X,I0))') &
                        "BPG3D: node ", NH, " is not a vertex of neighbor element ", &
                        NEle, "; element vertices:", NM(NEle, :)
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(bpg_error_msg))
                  end if
                  N1 = NM(NEle, LocalVertex)
                  N2 = NM(NEle, mod(LocalVertex, 3) + 1)
                  N3 = NM(NEle, mod(LocalVertex + 1, 3) + 1)
                  ! Skip dry elements
                  NCEle = NODECODE(N1)*NODECODE(N2) &
                          *NODECODE(N3)*NOFF(NEle)
                  if (NCEle == 0) cycle
                  ! For vertex of interest
                  BCPN1 = BCP(N1, k)
                  Zk = HsOAMB*(Sigma(k) - B) - DP(N1)
                  ! For first neighbor vertex
                  HsN2 = DP(N2) + dble(IFNLFA)*Eta2(N2)
                  if (HsN2 < H0) HsN2 = H0
                  SigmaNN = B + AMB*(Zk + DP(N2))/HsN2
                  !Casey 221103 : User control over cubic interpolation
                  if (verticalInterpolationSchemeID == VERTICAL_INTERPOLATION_LINEAR) then
                     call ZSURFBUOY(SigmaNN, BCPN2, N2, k)
                  else ! Default
                     call ZSURFBUOY_p3(SigmaNN, BCPN2, N2, k)
                  end if
                  ! For second neighbor vertex
                  HsN2 = DP(N3) + dble(IFNLFA)*Eta2(N3)
                  if (HsN2 < H0) HsN2 = H0
                  SigmaNN = B + AMB*(Zk + DP(N3))/HsN2
                  !Casey 221103 : User control over cubic interpolation
                  if (verticalInterpolationSchemeID == VERTICAL_INTERPOLATION_LINEAR) then
                     call ZSURFBUOY(SigmaNN, BCPN3, N3, k)
                  else ! Default
                     call ZSURFBUOY_p3(SigmaNN, BCPN3, N3, k)
                  end if
                  ! Compute gradients
                  if ((BCPN2 > -990.d0) .and. (BCPN3 > -990.d0)) then
                     TotalBCPGArea2 = TotalBCPGArea2 + Areas(NEle)
                     SFacAvg = (SFAC(N1) + SFAC(N2) + SFAC(N3))/3.d0
                     a1 = X(N3) - X(N2)
                     a2 = X(N1) - X(N3)
                     a3 = X(N2) - X(N1)
                     b1 = (Y(N2) - Y(N3))*SFacAvg
                     b2 = (Y(N3) - Y(N1))*SFacAvg
                     b3 = (Y(N1) - Y(N2))*SFacAvg
                     DBCPDX2A = DBCPDX2A + (BCPN1*b1 + BCPN2*b2 + BCPN3*b3)
                     DBCPDY2A = DBCPDY2A + (BCPN1*a1 + BCPN2*a2 + BCPN3*a3)
                  end if
               end do
               if (TotalBCPGArea2 > 0.0001d0) then
                  BCPG(k) = (cmplx(DBCPDX2A, kind=8) + iy*cmplx(DBCPDY2A, kind=8))/cmplx(TotalBCPGArea2, kind=8)
                  !Casey 230317 Debug
                  ! BCPG(k)=(0.D0,0.D0)
               end if
            end do

            do k = 1, NFEN
               BPG(NH, k) = BCPG(k)
            end do

            !     Finished computing baroclinic pressure gradient (computed in level
            !     coordinates) at each node in the vertical

            !     Compute vertically integrated baroclinic pressure gradient for use
            !     in the wave equation.  NOTE: For a prognostic model in which the
            !     density field evolves in time, this calculation should be done
            !     after the new density field is computed.  In this case one would
            !     integrate over the vertical first and differentiate second.

            ! this old scheme has been replaced by the new one below ...
            !Casey 230120 : This integration has been replaced below.
            if (densityIntegrationID == DENSITY_INTEGRATION_TRAPEZOIDAL) then
               VIBCPG = (0.d0, 0.d0)
               do k = NFEN - 1, 1, -1
                  VIBCPG = VIBCPG + cmplx(0.5d0, kind=8)*(BCPG(k + 1) + BCPG(k)) &
                           *cmplx(Sigma(k + 1) - Sigma(k), kind=8)
               end do
               VIDBCPDXOH(NH) = real(VIBCPG)/AMB
               VIDBCPDYOH(NH) = aimag(VIBCPG)/AMB
               !Casey 230120 : This integration has been replaced below.
            elseif (densityIntegrationID == DENSITY_INTEGRATION_SIMPSON) then
               !--------------------------------------------------------------------------------------------------
               ! arash 25/2/2016: Simpson integration, using quadratic function interpolation.
               !------------------------------------------------------------------------------
               ! 1. initialize
               VIBCPG = (0.d0, 0.d0)

               ! 2. to compute the integral at node NFEN-1, we need the mid-point value, which we evaluate by using a quadratic
               !    interpolation, using nodes NFEN, NFEN-1, and NFEN-2.
               SigmaMP = (Sigma(NFEN) + Sigma(NFEN - 1))/2.0d0

               f1 = ((SigmaMP - Sigma(NFEN - 1))*(SigmaMP - Sigma(NFEN - 2))) &
                    /((Sigma(NFEN) - Sigma(NFEN - 1))*(Sigma(NFEN) - Sigma(NFEN - 2)))
               f2 = ((SigmaMP - Sigma(NFEN))*(SigmaMP - Sigma(NFEN - 2))) &
                    /((Sigma(NFEN - 1) - Sigma(NFEN))*(Sigma(NFEN - 1) - Sigma(NFEN - 2)))
               f3 = ((SigmaMP - Sigma(NFEN))*(SigmaMP - Sigma(NFEN - 1))) &
                    /((Sigma(NFEN - 2) - Sigma(NFEN))*(Sigma(NFEN - 2) - Sigma(NFEN - 1)))

               BCPG_mid_point = cmplx(f1, kind=8)*BCPG(NFEN) + cmplx(f2, kind=8)*BCPG(NFEN - 1) + &
                                cmplx(f3, kind=8)*BCPG(NFEN - 2)

               ! 3. evaluate the integral at node NFEN-1
               VIBCPG = cmplx(one_sixth*(Sigma(NFEN) - Sigma(NFEN - 1)), kind=8) &
                        *(BCPG(NFEN) + cmplx(4.0d0, kind=8)*BCPG_mid_point + BCPG(NFEN - 1))

               ! 4. evaluate the integral at the rest of the vertical nodes
               do k = NFEN - 2, 1, -1

                  ! 4.1.
                  SigmaMP = (Sigma(k + 1) + Sigma(k))/2.0d0

                  f1 = ((SigmaMP - Sigma(k + 1))*(SigmaMP - Sigma(k))) &
                       /((Sigma(k + 2) - Sigma(k + 1))*(Sigma(k + 2) - Sigma(k)))
                  f2 = ((SigmaMP - Sigma(k + 2))*(SigmaMP - Sigma(k))) &
                       /((Sigma(k + 1) - Sigma(k + 2))*(Sigma(k + 1) - Sigma(k)))
                  f3 = ((SigmaMP - Sigma(k + 2))*(SigmaMP - Sigma(k + 1))) &
                       /((Sigma(k) - Sigma(k + 2))*(Sigma(k) - Sigma(k + 1)))

                  BCPG_mid_point = cmplx(f1, kind=8)*BCPG(k + 2) + cmplx(f2, kind=8)*BCPG(k + 1) + &
                                   cmplx(f3, kind=8)*BCPG(k)

                  ! 4.2. evaluate the integral
                  VIBCPG = VIBCPG + cmplx(one_sixth*(Sigma(k + 1) - Sigma(k)), kind=8) &
                           *(BCPG(k + 1) + cmplx(4.0d0, kind=8)*BCPG_mid_point + BCPG(k))

               end do

               ! 5. let's hope we are better of now.
               !-----------  -------------------------------------------------------------------
               !-----------  -------------------------------------------------------------------
               VIDBCPDXOH(NH) = real(VIBCPG)/AMB
               VIDBCPDYOH(NH) = aimag(VIBCPG)/AMB
               !Casey 230120 : Try to improve this?
            else ! default is densityIntegration = 'quadrature'
               VIBCPG = (0.d0, 0.d0)
               do k = NFEN - 1, 1, -1
                  do IGQ = NumGQP, 1, -1
                     GQI = NumGQP*(k - 1) + IGQ
                     ! Determine sigma at quadrature point
                     midpoint = 0.5d0*(Sigma(k) + Sigma(k + 1))
                     halfdist = midpoint - Sigma(k)
                     SigmaMP = midpoint + GQP(IGQ)*halfdist
                     ! Interpolation for BCPG ...
                     ! ... quadratic.
                     !                 IIP = (/ k-1 , k , k+1 /)
                     !                 IF(k.EQ.1)      IIP = IIP + 1
                     ! ... cubic.
                     IIP = [k - 1, k, k + 1, k + 2]
                     if (k == 1) IIP = IIP + 1
                     if (k == NFEN - 1) IIP = IIP - 1
                     ! ... quartic.
                     !                 IIP = (/ k-2 , k-1 , k , k+1 , k+2 /)
                     !                 IF(k.EQ.1)      IIP = IIP + 2
                     !                 IF(k.EQ.2)      IIP = IIP + 1
                     !                 IF(k.EQ.NFEN-1) IIP = IIP - 1
                     w = 1.d0
                     BCPGM(GQI) = (0.d0, 0.d0)
                     do IIP1 = 1, NumIP
                        ! Multiply up the weight for this term.
                        do IIP2 = 1, NumIP
                           if (IIP1 /= IIP2) then
                              w(IIP1) = w(IIP1) &
                                        *(SigmaMP - Sigma(IIP(IIP2))) &
                                        /(Sigma(IIP(IIP1)) - Sigma(IIP(IIP2)))
                           end if
                        end do
                        ! Add to the weighted interpolation.
                        BCPGM(GQI) = BCPGM(GQI) + cmplx(w(IIP1), kind=8)*BCPG(IIP(IIP1))
                     end do
                     VIBCPG = VIBCPG + cmplx((Sigma(k + 1) - Sigma(k)) &
                                             *GQW(IGQ)/sum(GQW), kind=8)*BCPGM(GQI)
                  end do
               end do
               VIDBCPDXOH(NH) = real(VIBCPG)/AMB
               VIDBCPDYOH(NH) = aimag(VIBCPG)/AMB
            end if

         end if

         !     Finished computing baroclinic terms

      end do

      !     Finish loop over horizontal nodes to compute the horizontal velocity

#ifdef CMPI
!     Update new 3D baroclinic pressure gradient and the vertically
!     integrated baroclinic pressure gradient on ghost nodes
!
      call UPDATEC3D(BPG)
      call UPDATER(VIDBCPDXOH, VIDBCPDYOH, DUMV1, 2)
#endif

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

      return
   end subroutine BPG3D
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Interpolate baroclinic pressure to specified sigma level (linear)
   !>
   !> Linearly interpolates the baroclinic pressure (BCP) to a specified sigma
   !> coordinate value. Uses an initial guess for the closest sigma level and
   !> searches up or down to find bracketing levels.
   !>
   !> @param[in]  SigmaNN   Target sigma coordinate value
   !> @param[out] BCPressNN Interpolated baroclinic pressure (or -999 if invalid)
   !> @param[in]  NN        Horizontal node index
   !> @param[in]  J         Initial guess for vertical level index
   !-----------------------------------------------------------------------
   subroutine ZSURFBUOY(SigmaNN, BCPressNN, NN, J)

      use GLOBAL_3DVS, only: A, B, BCP, NFEN, SIGMA
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      real(8), intent(in) :: SigmaNN !Sigma value of a neighbor node
      real(8), intent(out) :: BCPressNN
      integer, intent(in) :: NN
      integer, intent(in) :: J
      integer :: LBelo
      integer :: LAbov
      integer :: LTry
      integer :: IDiag
      real(8) :: SigBelo
      real(8) :: SigAbov
      real(8) :: SigTry
      logical :: needs_interpolation

      LOG_SCOPE_TRACED("zsurfbuoy", TIMESTEP_TRACING)

      IDiag = 0
      LBelo = 1
      LAbov = 1
      SigBelo = 0.d0
      SigAbov = 0.d0
      needs_interpolation = .true.

      if (SigmaNN <= 1.0001d0*B) then
         ! Into ground - skip
         SigBelo = -999.d0
         SigAbov = -999.d0
         BCPressNN = -999.d0
         needs_interpolation = .false.
      else if ((SigmaNN > 1.0001d0*B) .and. (SigmaNN <= B)) then
         ! At bottom - use bottom value
         LBelo = 1
         BCPressNN = BCP(NN, LBelo)
         SigBelo = B
         SigAbov = B
         needs_interpolation = .false.
      else if (SigmaNN >= A) then
         ! Into air - use surface value
         LAbov = NFEN
         BCPressNN = BCP(NN, LAbov)
         SigBelo = A
         SigAbov = A
         needs_interpolation = .false.
      else
         ! Search for bracketing sigma levels
         LTry = J
         SigTry = Sigma(LTry)
         if (SigmaNN > SigTry) then
            ! Too low - search upward
            SigBelo = SigTry
            LBelo = LTry
            LTry = LTry + 1
            SigTry = Sigma(LTry)
            do while (SigmaNN > SigTry)
               SigBelo = SigTry
               LBelo = LTry
               LTry = LTry + 1
               SigTry = Sigma(LTry)
            end do
            SigAbov = SigTry
            LAbov = LTry
         else
            ! Too high - search downward
            SigAbov = SigTry
            LAbov = LTry
            LTry = LTry - 1
            SigTry = Sigma(LTry)
            do while (SigmaNN <= SigTry)
               SigAbov = SigTry
               LAbov = LTry
               LTry = LTry - 1
               SigTry = Sigma(LTry)
            end do
            SigBelo = SigTry
            LBelo = LTry
         end if
      end if

      if (needs_interpolation) then
         BCPressNN = (BCP(NN, LAbov) - BCP(NN, LBelo)) &
                     *(SigmaNN - SigBelo)/(SigAbov - SigBelo) + BCP(NN, LBelo)
      end if

      if (IDiag == 2) then
         write (2, *) '******** ZSURFBUOY **********'
         write (2, *) '     NH  NV  SigmaNN   SigBelo   SigAbov', &
            '      BCPressNN'
         write (2, 777) NN, J, SigmaNN, SigBelo, SigAbov, BCPressNN
777      format(I7, I5, 3(F10.3), E14.5)
      end if

   end subroutine ZSURFBUOY
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Interpolate baroclinic pressure to specified sigma level (cubic)
   !>
   !> Interpolates the baroclinic pressure (BCP) to a specified sigma coordinate
   !> using cubic polynomial interpolation for improved accuracy. Uses Lagrangian
   !> interpolation with four surrounding sigma levels.
   !>
   !> @param[in]  SigmaNN   Target sigma coordinate value
   !> @param[out] BCPressNN Interpolated baroclinic pressure (or -999 if invalid)
   !> @param[in]  NN        Horizontal node index
   !> @param[in]  J         Initial guess for vertical level index
   !-----------------------------------------------------------------------
   subroutine ZSURFBUOY_p3(SigmaNN, BCPressNN, NN, J)

      use GLOBAL_3DVS, only: A, B, BCP, NFEN, SIGMA
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      real(8), intent(in) :: SigmaNN !Sigma value of a neighbor node
      real(8), intent(out) :: BCPressNN
      integer, intent(in) :: NN
      integer, intent(in) :: J
      integer :: LBelo
      integer :: LAbov
      integer :: LTry
      integer :: IDiag
      real(8) :: SigBelo
      real(8) :: SigAbov
      real(8) :: SigTry
      ! Variables for cubic interpolation
      integer :: k1, k2, k3, k4 ! nodes used for cubic polynomial
      real(8) :: f1, f2, f3, f4 ! shape functions
      real(8) :: top, bot
      logical :: needs_interpolation

      LOG_SCOPE_TRACED("zsurfbuoy_p3", TIMESTEP_TRACING)

      IDiag = 0
      SigBelo = 0.d0
      SigAbov = 0.d0
      LBelo = 1
      LAbov = 1
      needs_interpolation = .true.

      if (SigmaNN <= (1.000000010d0*B)) then
         ! Into ground - skip (Casey 221122: more restrictive check)
         SigBelo = -999.d0
         SigAbov = -999.d0
         BCPressNN = -999.d0
         needs_interpolation = .false.
      else if ((SigmaNN > (1.00000001d0*B)) .and. (SigmaNN <= B)) then
         ! At bottom - use bottom value
         LBelo = 1
         BCPressNN = BCP(NN, LBelo)
         SigBelo = B
         SigAbov = B
         needs_interpolation = .false.
      else if (SigmaNN >= (1.0001d0*A)) then
         ! In the air - non-existent (arash 150513)
         BCPressNN = -999.0d0
         needs_interpolation = .false.
      else if ((SigmaNN >= A) .and. (SigmaNN < (1.0001d0*A))) then
         ! At surface
         LAbov = NFEN
         BCPressNN = BCP(NN, LAbov)
         SigBelo = A
         SigAbov = A
         needs_interpolation = .false.
      else
         ! Search for bracketing sigma levels
         LTry = J
         SigTry = Sigma(LTry)
         if (SigmaNN > SigTry) then
            ! Too low - search upward
            SigBelo = SigTry
            LBelo = LTry
            LTry = LTry + 1
            SigTry = Sigma(LTry)
            do while (SigmaNN > SigTry)
               SigBelo = SigTry
               LBelo = LTry
               LTry = LTry + 1
               SigTry = Sigma(LTry)
            end do
            SigAbov = SigTry
            LAbov = LTry
         else
            ! Too high - search downward
            SigAbov = SigTry
            LAbov = LTry
            LTry = LTry - 1
            SigTry = Sigma(LTry)
            do while (SigmaNN <= SigTry)
               SigAbov = SigTry
               LAbov = LTry
               LTry = LTry - 1
               SigTry = Sigma(LTry)
            end do
            SigBelo = SigTry
            LBelo = LTry
         end if
      end if

      if (needs_interpolation) then
         ! Cubic polynomial interpolation
         ! 1- decide which nodes to pick: k1:lowest, k2, k3, k4:top-most
         if (LAbov == NFEN) then
            k4 = NFEN
            k3 = NFEN - 1
            k2 = NFEN - 2
            k1 = NFEN - 3
         else if (LBelo == 1) then
            k4 = 4
            k3 = 3
            k2 = 2
            k1 = 1
         else
            k4 = LAbov + 1
            k3 = LAbov
            k2 = LBelo
            k1 = LBelo - 1
         end if

         ! 2- compute the shape functions
         top = ((SigmaNN - Sigma(k2))*(SigmaNN - Sigma(k3))* &
                (SigmaNN - Sigma(k4)))
         bot = ((Sigma(k1) - Sigma(k2))*(Sigma(k1) - Sigma(k3))* &
                (Sigma(k1) - Sigma(k4)))
         f1 = top/bot

         top = ((SigmaNN - Sigma(k1))*(SigmaNN - Sigma(k3))* &
                (SigmaNN - Sigma(k4)))
         bot = ((Sigma(k2) - Sigma(k1))*(Sigma(k2) - Sigma(k3))* &
                (Sigma(k2) - Sigma(k4)))
         f2 = top/bot

         top = ((SigmaNN - Sigma(k1))*(SigmaNN - Sigma(k2))* &
                (SigmaNN - Sigma(k4)))
         bot = ((Sigma(k3) - Sigma(k1))*(Sigma(k3) - Sigma(k2))* &
                (Sigma(k3) - Sigma(k4)))
         f3 = top/bot

         top = ((SigmaNN - Sigma(k1))*(SigmaNN - Sigma(k2))* &
                (SigmaNN - Sigma(k3)))
         bot = ((Sigma(k4) - Sigma(k1))*(Sigma(k4) - Sigma(k2))* &
                (Sigma(k4) - Sigma(k3)))
         f4 = top/bot

         ! 3- interpolate
         BCPressNN = f1*BCP(NN, k1) + f2*BCP(NN, k2) &
                     + f3*BCP(NN, k3) + f4*BCP(NN, k4)
      end if

      if (IDiag == 2) then
         write (2, *) '******** ZSURFBUOY_p3 **********'
         write (2, *) '     NH  NV  SigmaNN   SigBelo   SigAbov', &
            '      BCPressNN'
         write (2, 777) NN, J, SigmaNN, SigBelo, SigAbov, BCPressNN
777      format(I7, I5, 3(F10.3), E14.5)
      end if

   end subroutine ZSURFBUOY_p3
   !-----------------------------------------------------------------------

end module mod_baroclinic
