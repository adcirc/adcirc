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
!-----------------------------------------------------------------------
!  MODULE FLUX_BOUNDARIES
!-----------------------------------------------------------------------
!> @brief Flux-type boundary conditions applied each timestep
!>
!> This module computes the normal flux (QN1/QN2) and elevation (EN2)
!> boundary values used by the GWCE and momentum equations for:
!> - Specified normal flow boundaries (periodic or fort.20 time series),
!>   including cold start/hotstart setup and conversion to rotated
!>   coordinate systems
!> - Radiation and zero normal velocity gradient boundaries
!> - External and internal barrier (weir) boundaries, including vertical
!>   element walls and cross-barrier pipes
!> - Subdomain modeling boundary data
!-----------------------------------------------------------------------
module mod_flux_boundaries

   implicit none

   private

   public :: applyBoundaryConditions, initializeNormalFlowBoundary, &
             hotstartNormalFlowBoundary

   ! Whether the river radiation b.c. elevations still need to be filled
   ! from EtaDisc after flux settling
   logical :: EtaDisc_Fill = .true.

   real(8), parameter :: eps = epsilon(1.0d0)

   integer, parameter :: NORMAL_FLOW_UNIT = 20

   ! Where ElevDisc for river radiation (LBCODE=52) nodes is captured from
   ! on the current timestep
   integer, parameter :: ELEVDISC_NONE = 0
   integer, parameter :: ELEVDISC_FROM_ETA2 = 1
   integer, parameter :: ELEVDISC_FROM_ETADISC = 2

contains

!-----------------------------------------------------------------------
!> @brief Initialize aperiodic normal flow boundary conditions (cold start)
!>
!> Opens fort.20 and reads the snapshot interval and the first two
!> snapshots of normal flux (and elevation for LBCODE=32).
!-----------------------------------------------------------------------
   subroutine initializeNormalFlowBoundary()
      use sizes, only: INPUTDIR
      use global, only: FTIMINC, QTIME1, QTIME2, STATIM, &
                        QNIN1, QNIN2, ENIN1, ENIN2
      implicit none

      open (NORMAL_FLOW_UNIT, file=trim(INPUTDIR)//'/'//'fort.20')
      read (NORMAL_FLOW_UNIT, *) FTIMINC

      QNIN1 = 0d0
      QNIN2 = 0d0
      call readNormalFlowSnapshot(QNIN1, ENIN1)
      call readNormalFlowSnapshot(QNIN2, ENIN2)
      QTIME1 = STATIM*86400.d0
      QTIME2 = QTIME1 + FTIMINC

   end subroutine initializeNormalFlowBoundary
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Restore specified normal flow boundary conditions on hotstart
!>
!> Computes the normal flux (QN1, QN2) and elevation (EN1, EN2) boundary
!> values at the previous and current hotstart time levels.
!>
!> - NFFR = 0:  advances through fort.20 from the start of the run
!> - NFFR = -1: starts fort.20 at the hotstart time
!> - NFFR > 0:  evaluates the periodic forcing
!>
!> @param[in] TimeLoc       Simulation time of the hotstart (seconds)
!> @param[in] ITHS          Timestep number of the hotstart
!> @param[in] RampExtFlux1  External flux ramp at timestep ITHS-1
!> @param[in] RampExtFlux2  External flux ramp at timestep ITHS
!-----------------------------------------------------------------------
   subroutine hotstartNormalFlowBoundary(TimeLoc, ITHS, RampExtFlux1, RampExtFlux2)
      use sizes, only: INPUTDIR
      use global, only: NFFR, FTIMINC, QTIME1, QTIME2, STATIM, DTDPHS, &
                        QNIN1, QNIN2, ENIN1, ENIN2, QN1, QN2, EN1, EN2, &
                        FPER, FAMIG, FFACE, FFF, QNAM, QNPH
      use mesh, only: ICS
      use boundaries, only: NVEL
      implicit none

      real(8), intent(in) :: TimeLoc
      integer, intent(in) :: ITHS
      real(8), intent(in) :: RampExtFlux1, RampExtFlux2
      real(8) :: TimePrevious, QTRATIO
      real(8) :: ARGJ1, ARGJ2, ARG1, ARG2, RFF1, RFF2
      integer :: I, J, IT, NCYC

      if ((NFFR == 0) .or. (NFFR == -1)) then
         open (NORMAL_FLOW_UNIT, file=trim(INPUTDIR)//'/'//'fort.20')
         read (NORMAL_FLOW_UNIT, *) FTIMINC
         if (NFFR == 0) then
            QTIME1 = STATIM*86400.d0
         else
            QTIME1 = TimeLoc
         end if
         QTIME2 = QTIME1 + FTIMINC

         QNIN1(1:NVEL) = 0d0
         QNIN2(1:NVEL) = 0d0
         call readNormalFlowSnapshot(QNIN1, ENIN1)
         call readNormalFlowSnapshot(QNIN2, ENIN2)

         if (NFFR == 0) then
            do IT = 1, ITHS - 1
               call advanceNormalFlowSnapshots(dble(IT)*DTDPHS + STATIM*86400.d0)
            end do
            TimePrevious = dble(ITHS - 1)*DTDPHS + STATIM*86400.d0
         else
            TimePrevious = TimeLoc
         end if

         QTRATIO = (TimePrevious - QTIME1)/FTIMINC
         QN1 = RampExtFlux1*(QNIN1 + QTRATIO*(QNIN2 - QNIN1))
         EN1 = RampExtFlux1*(ENIN1 + QTRATIO*(ENIN2 - ENIN1))

         call advanceNormalFlowSnapshots(TimeLoc)

         QTRATIO = (TimeLoc - QTIME1)/FTIMINC
         QN2 = RampExtFlux2*(QNIN1 + QTRATIO*(QNIN2 - QNIN1))
         EN2 = RampExtFlux2*(ENIN1 + QTRATIO*(ENIN2 - ENIN1))

      else if (NFFR > 0) then
         do J = 1, NFFR
            if (abs(FPER(J)) <= eps) then
               NCYC = 0
            else
               NCYC = int(TimeLoc/FPER(J))
            end if
            ARGJ1 = FAMIG(J)*(TimeLoc - DTDPHS - dble(NCYC)*FPER(J)) + FFACE(J)
            ARGJ2 = FAMIG(J)*(TimeLoc - dble(NCYC)*FPER(J)) + FFACE(J)
            RFF1 = FFF(J)*RampExtFlux1
            RFF2 = FFF(J)*RampExtFlux2
            do I = 1, NVEL
               ARG1 = ARGJ1 - QNPH(J, I)
               ARG2 = ARGJ2 - QNPH(J, I)
               QN1(I) = QN1(I) + rotate_normal_flux(ICS, I, QNAM(J, I)*RFF1*cos(ARG1))
               QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, QNAM(J, I)*RFF2*cos(ARG2))
            end do
         end do
      end if

   end subroutine hotstartNormalFlowBoundary
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Read one fort.20 snapshot of normal flux (and elevation)
!>
!> Reads a normal flux for each LBCODE=2,12,22 node and a normal flux and
!> elevation for each LBCODE=32 node. Entries for other nodes are left
!> unchanged.
!>
!> @param[inout] QNIN   Normal flux at each flow boundary node
!> @param[inout] ENIN   Elevation at each flow boundary node
!-----------------------------------------------------------------------
   subroutine readNormalFlowSnapshot(QNIN, ENIN)
      use mesh, only: ICS
      use boundaries, only: NVEL, LBCODEI
      implicit none

      real(8), intent(inout) :: QNIN(:), ENIN(:)
      real(8) :: QNIN_TMP
      integer :: J

      do J = 1, NVEL
         select case (LBCODEI(J))
         case (2, 12, 22)
            read (NORMAL_FLOW_UNIT, *) QNIN_TMP
         case (32)
            read (NORMAL_FLOW_UNIT, *) QNIN_TMP, ENIN(J)
         case default
            cycle
         end select
         QNIN(J) = rotate_normal_flux(ICS, J, QNIN_TMP)
      end do

   end subroutine readNormalFlowSnapshot
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Shift to the next fort.20 snapshot if TimeLoc is past QTIME2
!>
!> @param[in] TimeLoc Simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine advanceNormalFlowSnapshots(TimeLoc)
      use global, only: QTIME1, QTIME2, FTIMINC, QNIN1, QNIN2, ENIN1, ENIN2
      use boundaries, only: NVEL, LBCODEI
      implicit none

      real(8), intent(in) :: TimeLoc
      integer :: J

      if (TimeLoc > QTIME2) then
         QTIME1 = QTIME2
         QTIME2 = QTIME2 + FTIMINC
         do J = 1, NVEL
            select case (LBCODEI(J))
            case (2, 12, 22)
               QNIN1(J) = QNIN2(J)
            case (32)
               QNIN1(J) = QNIN2(J)
               ENIN1(J) = ENIN2(J)
            end select
         end do
         call readNormalFlowSnapshot(QNIN2, ENIN2)
      end if

   end subroutine advanceNormalFlowSnapshots
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply all flux boundary conditions
!>
!> Computes the per-timestep forcing state (fort.20 snapshot advance,
!> periodic constituent arguments, river radiation elevation capture), then
!> visits each flow boundary segment once and applies the treatment for its
!> boundary type. Subdomain boundary data is applied last.
!>
!> @param[in] IT      Current timestep number
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] TimeH   Harmonic analysis reference time (seconds)
!-----------------------------------------------------------------------
   subroutine applyBoundaryConditions(IT, TimeLoc, TimeH)
      use global, only: NFFR, FPER, FAMIG, FFACE, FFF, RampExtFlux, &
                        QTIME1, FTIMINC, FluxSettlingIT, EtaDisc, Eta2, &
                        NIBNODECODE
      use mesh, only: NP
      use boundaries, only: NBOU, NVELL, IBTYPE, NFLUXF, NFLUXIB
      use subdomain, only: subdomainOn
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc, TimeH

      real(8) :: ARGJ(max(NFFR, 0)), RFF(max(NFFR, 0))
      real(8) :: QTRATIO
      integer :: elevDiscSource
      integer :: I, J, K, NCYC

      ! Per-timestep state for specified normal flow boundaries
      QTRATIO = 0.d0
      elevDiscSource = ELEVDISC_NONE
      if (NFLUXF == 1) then
         if (NFFR > 0) then
            do J = 1, NFFR
               if (abs(FPER(J)) <= eps) then
                  NCYC = 0
               else
                  NCYC = int(TimeH/FPER(J))
               end if
               ARGJ(J) = FAMIG(J)*(TimeH - dble(NCYC)*FPER(J)) + FFACE(J)
               RFF(J) = FFF(J)*RampExtFlux
            end do
         elseif ((NFFR == 0) .or. (NFFR == -1)) then
            call advanceNormalFlowSnapshots(TimeLoc)
            QTRATIO = (TimeLoc - QTIME1)/FTIMINC
         end if

         ! Collect elevation information for river radiation b.c.
         if (IT == FluxSettlingIT) then
            EtaDisc_Fill = .false.
            EtaDisc(1:NP) = Eta2(1:NP)
            elevDiscSource = ELEVDISC_FROM_ETA2
         else if (EtaDisc_Fill .and. IT > FluxSettlingIT) then
            EtaDisc_Fill = .false.
            elevDiscSource = ELEVDISC_FROM_ETADISC
         end if
      end if

      if (NFLUXIB == 1) NIBNODECODE(:) = 0

      I = 0
      do K = 1, NBOU
         select case (IBTYPE(K))
         case (2, 12, 22, 32, 52)
            if (NFLUXF == 1) call specifiedNormalFlowSegment(I, NVELL(K), ARGJ, RFF, QTRATIO, elevDiscSource)
         case (30, 40, 41)
            call velocityDischargeSegment(I, NVELL(K))
         case (3, 13, 23)
            call externalBarrierSegment(I, NVELL(K), TimeLoc)
         case (4, 24)
            call internalBarrierSegment(K, I, NVELL(K), TimeLoc, pipes=.false.)
         case (5, 25)
            call internalBarrierSegment(K, I, NVELL(K), TimeLoc, pipes=.true.)
         case (64)
            call vewBarrierSegment(K, I, NVELL(K), TimeLoc)
         end select
         I = I + segmentNodeCount(IBTYPE(K), NVELL(K))
      end do

      if (subdomainOn) call applySubdomainBCData(IT)

   end subroutine applyBoundaryConditions
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Number of flow boundary node entries used by a boundary segment
!>
!> Internal barrier segments list each node pair twice (front side, then
!> back side), so they occupy 2*NVELL entries in the flow boundary arrays.
!-----------------------------------------------------------------------
   pure integer function segmentNodeCount(boundary_type, nvell_k) result(n)
      implicit none
      integer, intent(in) :: boundary_type, nvell_k

      select case (boundary_type)
      case (4, 24, 5, 25, 64)
         n = 2*nvell_k
      case default
         n = nvell_k
      end select

   end function segmentNodeCount
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Specified normal flow on one boundary segment
!>
!> Sets the normal flux (QN2) and, for LBCODE=32, elevation (EN2) from
!> either the periodic constituents (NFFR > 0) or the fort.20 snapshots
!> (NFFR = 0, -1). For LBCODE=52 river radiation nodes, captures the
!> discharge elevation (ElevDisc) once flux settling has ended.
!>
!> @param[in] I0             Flow boundary index before the segment
!> @param[in] N              Number of nodes in the segment
!> @param[in] ARGJ           Periodic constituent arguments
!> @param[in] RFF            Ramped periodic constituent nodal factors
!> @param[in] QTRATIO        fort.20 interpolation weight
!> @param[in] elevDiscSource Source for ElevDisc this timestep, if any
!-----------------------------------------------------------------------
   subroutine specifiedNormalFlowSegment(I0, N, ARGJ, RFF, QTRATIO, elevDiscSource)
      use global, only: NFFR, RampExtFlux, QNAM, QNPH, ENAM, ENPH, QN2, EN2, &
                        QNIN1, QNIN2, ENIN1, ENIN2, EtaDisc, ElevDisc, Eta2
      use mesh, only: ICS
      use boundaries, only: LBCODEI, NBV
      implicit none

      integer, intent(in) :: I0, N
      real(8), intent(in) :: ARGJ(:), RFF(:)
      real(8), intent(in) :: QTRATIO
      integer, intent(in) :: elevDiscSource
      integer :: I, J

      do I = I0 + 1, I0 + N
         if (NFFR > 0) then
            do J = 1, NFFR
               select case (LBCODEI(I))
               case (2, 12, 22, 52)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF(J)*cos(ARGJ(J) - QNPH(J, I)))
               case (32)
                  QN2(I) = QN2(I) + rotate_normal_flux(ICS, I, &
                                                       QNAM(J, I)*RFF(J)*cos(ARGJ(J) - QNPH(J, I)))
                  EN2(I) = EN2(I) + ENAM(J, I)*RFF(J)*cos(ARGJ(J) - ENPH(J, I))
               end select
            end do
         elseif ((NFFR == 0) .or. (NFFR == -1)) then
            QN2(I) = RampExtFlux*(QNIN1(I) + QTRATIO*(QNIN2(I) - QNIN1(I)))
            EN2(I) = RampExtFlux*(ENIN1(I) + QTRATIO*(ENIN2(I) - ENIN1(I)))
         end if

         if (LBCODEI(I) == 52) then
            select case (elevDiscSource)
            case (ELEVDISC_FROM_ETA2)
               ElevDisc(I) = Eta2(NBV(I))
            case (ELEVDISC_FROM_ETADISC)
               ElevDisc(I) = EtaDisc(NBV(I))
            end select
         end if
      end do

   end subroutine specifiedNormalFlowSegment
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Velocity-based discharge on one boundary segment
!>
!> Computes the outward discharge (QN1) from the previous timestep
!> velocities and current depth for radiation (LBCODE=30) and zero normal
!> velocity gradient (LBCODE=40,41) boundaries.
!>
!> @param[in] I0 Flow boundary index before the segment
!> @param[in] N  Number of nodes in the segment
!-----------------------------------------------------------------------
   subroutine velocityDischargeSegment(I0, N)
      use global, only: UU1, VV1, H2, QN1
      use boundaries, only: NBV, CSII, SIII
      implicit none

      integer, intent(in) :: I0, N
      integer :: I, NNBB
      real(8) :: UN1

      do I = I0 + 1, I0 + N
         NNBB = NBV(I)
         UN1 = UU1(NNBB)*CSII(I) + VV1(NNBB)*SIII(I)
         QN1(I) = H2(NNBB)*UN1
      end do

   end subroutine velocityDischargeSegment
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief External barrier flux on one boundary segment
!>
!> Computes supercritical outward normal flow (QN2) over external barrier
!> boundary nodes (LBCODE=3,13,23) using weir formulas.
!>
!> @param[in] I0      Flow boundary index before the segment
!> @param[in] N       Number of nodes in the segment
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine externalBarrierSegment(I0, N, TimeLoc)
      use global, only: QN2
      use weir_flux, only: COMPUTE_EXTERNAL_BOUNDARY_FLUX
      implicit none

      integer, intent(in) :: I0, N
      real(8), intent(in) :: TimeLoc
      integer :: I

      do I = I0 + 1, I0 + N
         call COMPUTE_EXTERNAL_BOUNDARY_FLUX(I, TimeLoc, QN2(I))
      end do

   end subroutine externalBarrierSegment
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Internal barrier flux on one boundary segment
!>
!> Computes inward/outward normal flow (QN2) over both sides of an
!> internal barrier (LBCODE=4,24,5,25) using weir formulas, adding the
!> cross-barrier pipe flux for LBCODE=5,25.
!>
!> @param[in] K       Boundary segment number
!> @param[in] I0      Flow boundary index before the segment
!> @param[in] N       Number of node pairs in the segment
!> @param[in] TimeLoc Current simulation time (seconds)
!> @param[in] pipes   Whether the barrier has cross-barrier pipes
!-----------------------------------------------------------------------
   subroutine internalBarrierSegment(K, I0, N, TimeLoc, pipes)
      use global, only: QN2
      use WEIR_FLUX, only: COMPUTE_INTERNAL_BOUNDARY_FLUX, &
                           COMPUTE_CROSS_BARRIER_PIPE_FLUX
      implicit none

      integer, intent(in) :: K, I0, N
      real(8), intent(in) :: TimeLoc
      logical, intent(in) :: pipes
      integer :: I, J
      real(8) :: PIPE_FLUX

      do J = 1, 2*N
         I = I0 + J
         call COMPUTE_INTERNAL_BOUNDARY_FLUX(I, J, K, TimeLoc, QN2(I))
         if (pipes) then
            PIPE_FLUX = 0d0
            call COMPUTE_CROSS_BARRIER_PIPE_FLUX(I, PIPE_FLUX)
            QN2(I) = QN2(I) + PIPE_FLUX
         end if
      end do

   end subroutine internalBarrierSegment
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Vertical element wall flux on one boundary segment
!>
!> Sets the submergence flags (ISSUBMERGED64) from the front side of a
!> vertical element wall barrier (LBCODE=64), which also sets the paired
!> back side flags, then computes the normal flow (QN2) on both sides.
!>
!> @param[in] K       Boundary segment number
!> @param[in] I0      Flow boundary index before the segment
!> @param[in] N       Number of node pairs in the segment
!> @param[in] TimeLoc Current simulation time (seconds)
!-----------------------------------------------------------------------
   subroutine vewBarrierSegment(K, I0, N, TimeLoc)
      use global, only: QN2
      use WEIR_FLUX, only: SET_SUBMERGED64_AT, COMPUTE_INTERNAL_BOUNDARY64_FLUX
      implicit none

      integer, intent(in) :: K, I0, N
      real(8), intent(in) :: TimeLoc
      integer :: I, J

      do J = 1, N
         call SET_SUBMERGED64_AT(I0 + J, J, K, TimeLoc)
      end do

      do J = 1, 2*N
         I = I0 + J
         call COMPUTE_INTERNAL_BOUNDARY64_FLUX(I, J, K, TimeLoc, QN2(I))
      end do

   end subroutine vewBarrierSegment
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
!> @brief Convert a boundary normal flux to the rotated coordinate system
!>
!> For rotated coordinate systems (ICS = 20-24), the normal flux Q at flow
!> boundary node IDX is decomposed into flux components along the rotated
!> coordinate axes using the boundary normal in rotated coordinates,
!> optionally rotated again for IFSPROTS = 1, and projected back onto the
!> boundary normal used by the model. Otherwise Q is returned unchanged.
!>
!> @param[in] ICS Coordinate system flag
!> @param[in] IDX Flow boundary node index
!> @param[in] Q   Normal flux in the input coordinate system
!> @return        Normal flux in the model coordinate system
!-----------------------------------------------------------------------
   real(8) function rotate_normal_flux(ICS, IDX, Q) result(QROT)
      use boundaries, only: NBV, CSII, SIII, CSII_OLD, SIII_OLD
      use mesh, only: SFMX, SFCX, SFCY, YCSFAC, RVELF
      use global, only: IFSPROTS

      implicit none

      integer, intent(in) :: ICS
      integer, intent(in) :: IDX
      real(8), intent(in) :: Q

      real(8) :: QX_R, QY_R, NMAG2
      real(8) :: QX_ROT, QY_ROT
      real(8) :: VO(2), VR(2)
      real(8) :: NX_R, NY_R

      select case (ICS)
      case (20:24)
         ! Boundary normal in rotated coordinates
         NX_R = CSII_OLD(IDX)/SFMX(NBV(IDX))
         NY_R = SIII_OLD(IDX)

         ! Flux vector with normal component Q and no tangential component:
         ! solves [NX_R NY_R; -NY_R NX_R] [QX_R; QY_R] = [Q; 0]. Written in
         ! closed form so a boundary segment that is exactly vertical in
         ! rotated coordinates (NY_R = 0) does not divide by zero.
         NMAG2 = NX_R*NX_R + NY_R*NY_R
         QX_R = Q*NX_R/NMAG2
         QY_R = Q*NY_R/NMAG2

         if (IFSPROTS == 1) then
            VO = [QX_R, QY_R]
            VR = matmul(RVELF(1:2, 1:2, IDX), VO)
            QX_ROT = VR(1)
            QY_ROT = VR(2)
         else
            QX_ROT = QX_R
            QY_ROT = QY_R
         end if

         ! Generate the new flux forcing terms
         QROT = SFCX(NBV(IDX))*QX_ROT*CSII(IDX) + SFCY(NBV(IDX))*QY_ROT*SIII(IDX)*YCSFAC(NBV(IDX))
      case default
         QROT = Q
      end select

   end function rotate_normal_flux
   !-----------------------------------------------------------------------

end module mod_flux_boundaries
