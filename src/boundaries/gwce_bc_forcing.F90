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
!
module gwce_bc_forcing

   implicit none

   public :: apply_boundary_conditions, apply_vew1d_and_condensed_nodes, &
             apply_subdomain_boundary_nodes, update_periodic_sponge_layer_nodes, &
             update_coef_periodic_sponge_layer_lumped, update_coef_periodic_sponge_layer_consistent

   private

contains

   ! *********************************************************************
   !> Apply the boundary conditions to the gwce fields
   !>
   !> @param IT The current time step
   !> @param TimeLoc The current time
   !> @param TimeH
   !> @param dt The time step increment
   !> @param esbin1
   !> @param esbin2
   !> @param etime1 Elevation boundary time 1
   !> @param etime2 Elevation boundary time 2
   !> @param etratio Ratio between etime1/etime2
   !> @param etiminc Time increment for the elevation boundary
   !> @param eta1 Elevation at previous time step
   !> @param eta2 ELevation at current time step
   !> @param etas Elevation solution from previous time step
   !> @param h1 Total depth at the previous time step
   !> @param ElevDisc Elevation at discharge boundaries
   !> @param en0
   !> @param en1
   !> @param en2
   !> @param Current elevation forcing ramp value
   !> @param Current meteorological ramp value
   !> @param gwce_lv The gwce load vector
   !>
   ! *********************************************************************
   subroutine apply_boundary_conditions(IT, TimeLoc, TimeH, dt, invertedBarometerOnElevationBoundary, &
                                        esbin1, esbin2, etime1, etime2, etratio, etiminc, eta1, eta2, etas, h1, &
                                        elevdisc, en0, en1, en2, rampelev, rampmete, ep, gwce_lv)

      use mesh, only: np, nneigh, neitab
      use sizes, only: mnei
      use boundaries, only: neta, nbv, nbd, npebc, nope, nvel, bndlen2o3, lbcodei, &
                            nfluxb, nfluxf, nfluxib, nfluxgbc, nfluxrbc
      use global, only: per, amig, face, ff, efa, emo, nbfr, pr2, ilump, obccoef, coefd, &
                        usingDynamicWaterLevelCorrection, dynamicWaterLevelCorrection1, &
                        dynamicWaterLevelCorrection2, RES_BC_FLAG, BCFLAG_LNM, CBaroclinic, &
                        LNM_BC, FluxSettlingIT, nodecode, qn0, qn1, qn2, nperseg, iperconn, &
                        nnperbc
      use subdomain, only: subdomainOn, enforceBN, enforceECB
      use nodalattributes, only: geoidoffset, loadGeoidOffset, tau0var
      use spongelayer, only: no_met_in_sponge
      use gwce_bc_forcing_impl, only: apply_periodic_boundaries, apply_aperiodic_boundaries, &
                                      apply_geoid_offset, apply_inverted_barometer_on_elevation_boundaries, &
                                      apply_dynamic_water_level_correction, apply_levels_of_no_motion_boundary_conditions, &
                                      apply_normal_flow_boundary_conditions, update_rhs_sponge_periodic_boundary_nodes, &
                                      update_load_vector

      implicit none

      real(8), intent(in)    :: h1(np)
      real(8), intent(in)    :: elevdisc(neta)
      real(8), intent(in)    :: en0(neta), en1(neta), en2(neta)
      real(8), intent(in)    :: etiminc, TimeLoc, TimeH
      real(8), intent(in)    :: rampelev, rampmete, dt, ep
      real(8), intent(inout) :: esbin1(neta)
      real(8), intent(inout) :: esbin2(neta)
      real(8), intent(inout) :: etime1, etime2, etratio
      real(8), intent(inout) :: eta1(np), eta2(np), etas(np)
      real(8), intent(inout) :: gwce_lv(np)

      integer, intent(in)    :: IT
      logical, intent(in)    :: invertedBarometerOnElevationBoundary

      call apply_periodic_boundaries(NBFR, NETA, NP, NBD, RampElev, TimeLoc, TimeH, PER, &
                                     AMIG, FACE, FF, EFA, EMO, Eta2)
      call apply_aperiodic_boundaries(subdomainOn, enforceBN, NETA, NPEBC, NP, NBD, RampElev, &
                                      TimeLoc, TimeH, ETIME1, ETIME2, ETIMINC, ESBIN1, ESBIN2, &
                                      ETRATIO, Eta2)
      call apply_geoid_offset(subdomainOn, LoadGeoidOffset, enforceBN, NETA, NP, NBD, GeoidOffset, ETA2)
      call apply_inverted_barometer_on_elevation_boundaries(NO_MET_IN_SPONGE, &
                                                            invertedBarometerOnElevationBoundary, NETA, &
                                                            NP, NBD, RampMete, PR2, ETA2)
      call apply_dynamic_water_level_correction(usingDynamicWaterLevelCorrection, np, neta, nbd, &
                                                dynamicWaterLevelCorrection1, eta2)
      call apply_levels_of_no_motion_boundary_conditions(RES_BC_FLAG, NP, NOPE, BCFLAG_LNM, CBaroclinic, NETA, NBD, &
                                                         LNM_BC, ETA2)
      call apply_normal_flow_boundary_conditions(IT, NFLUXF, NFLUXB, NFLUXIB, NFLUXGBC, NFLUXRBC, &
                                                 NVEL, NP, NBV, LBCODEI, NodeCode, FluxSettlingIT, &
                                                 QN0, QN1, QN2, ETAS, EN0, EN1, EN2, ETA1, ElevDisc, &
                                                 H1, TAU0Var, DT, BndLen2O3, GWCE_LV)
      call update_rhs_sponge_periodic_boundary_nodes(NP, NPERSEG, NNPERBC, IPERCONN, NODECODE, EP, ETAS, GWCE_LV)
      call update_load_vector(ILump, NETA, NP, MNEI, NBD, NNEIGH, NEITAB, NODECODE, ETA1, ETA2, &
                              OBCCOEF, COEFD, EP, ETAS, GWCE_LV)

   end subroutine apply_boundary_conditions

   ! *********************************************************************
   !> Enforce the subdomain boundary nodes in the gwce
   !>
   !> @param subdomainOn Flag to indicate if the subdomain model is active
   !> @param enforceBN Flag to enforce the boundary nodes
   ! *********************************************************************
   subroutine apply_subdomain_boundary_nodes(subdomainOn, enforceBN)
      use subdomain, only: enforceEob, enforceEib, enforceGWCELVob
      implicit none

      logical, intent(in) :: subdomainOn
      integer, intent(in) :: enforceBN

      if (subdomainOn .and. enforceBN == 2) then
         call enforceEob()
         call enforceEib()
         call enforceGWCELVob()
      end if
   end subroutine apply_subdomain_boundary_nodes

   ! *********************************************************************
   !> Sum the values at the front and back side of the vertical element wall (vew)
   !> 1d boundary nodes
   !>
   !> @param NFLUXIB64_GBL The number of global flux boundary nodes
   !> @param ILump The lumping flag
   !> @param NP The number of nodes
   !> @param NBOU The number of boundary nodes
   !> @param NVEL The number of velocity components
   !> @param NVELL The number of velocity components for each boundary node
   !> @param NBV The number of boundary nodes for each velocity component
   !> @param IBCONN The connectivity table for the boundary nodes
   !> @param LBCODEI The boundary code for each boundary node
   !> @param NODECODE The node code
   !> @param NListCondensedNodes The number of condensed nodes
   !> @param NNodesListCondensedNodes The number of nodes in the condensed nodes list
   !> @param ListCondensedNodes The list of condensed nodes
   !> @param ISSUBMERGED64 The submerged flag for the vew boundary nodes
   !> @param COEFD The coefficient for the elevation specified boundary nodes
   !> @param COEFDTempMem The temporary memory for the coefficient
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary coefficient
   ! *********************************************************************
   subroutine vew1d_sum_front_and_back_side(NFLUXIB64_GBL, ILump, NP, NBOU, NVEL, NVELL, NBV, IBCONN, LBCODEI, &
                                            NODECODE, NListCondensedNodes, NNodesListCondensedNodes, &
                                            ListCondensedNodes, ISSUBMERGED64, COEFD, COEFDTempMem, GWCE_LV, ETA1, COEFDTemp)
      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: ILump
      integer, intent(in) :: NP
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVEL
      integer, intent(in) :: NVELL(NBOU)
      integer, intent(in) :: NBV(NVEL)
      integer, intent(in) :: IBCONN(NVEL)
      integer, intent(in) :: LBCODEI(NVEL)
      integer, intent(in) :: NODECODE(NP)
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(NListCondensedNodes)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: ISSUBMERGED64(NVEL)
      real(8), intent(inout) :: COEFD(NP)
      real(8), intent(inout) :: COEFDTempMem(NP)
      real(8), intent(inout) :: GWCE_LV(NP)
      real(8), intent(inout) :: ETA1(NP)
      real(8), pointer, intent(inout) :: COEFDTemp(:)

      integer :: I, J, K, L, NNBB1, NNBB2

      if ((NFLUXIB64_GBL > 0) .and. (ILump /= 0)) then
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (64)
               do J = 1, NVELL(K)
                  I = I + 1
                  if (ISSUBMERGED64(I) /= 0) then
                     NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                     NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
                     if (NODECODE(NNBB1) /= 0 .and. NODECODE(NNBB2) /= 0) then
                        COEFDTemp(NNBB2) = COEFD(NNBB1) + COEFD(NNBB2)
                        GWCE_LV(NNBB2) = GWCE_LV(NNBB1) + GWCE_LV(NNBB2)
                        ETA1(NNBB2) = (ETA1(NNBB1) + ETA1(NNBB2))*0.5d0
                        GWCE_LV(NNBB1) = 0.d0 ! Set zero to avoid duplicated additions
                        ETA1(NNBB1) = ETA1(NNBB2) ! Substitute here to make duplicated averaging ineffective
                     end if
                  end if
               end do
               I = I + NVELL(K)
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case DEFAULT
               I = I + NVELL(K)
            end select
         end do
      end if

   end subroutine vew1d_sum_front_and_back_side

   ! *********************************************************************
   !> Sum the values at the condensed nodes
   !>
   !> @param LoadCondensedNodes Flag to indicate if the condensed nodes are loaded
   !> @param ILump The lumping flag
   !> @param NP The number of nodes
   !> @param NListCondensedNodes The number of condensed nodes
   !> @param NNodesListCondensedNodes The number of nodes in the condensed nodes list
   !> @param ListCondensedNodes The list of condensed nodes
   !> @param NODECODE The node code
   !> @param COEFD The gwce matrix diagonal entries
   !> @param COEFDTempMem The memory for the temporary diagonal entries
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary matrix diagonal entries
   !> @param LoadCondensedNodes Flag to indicate if the condensed nodes are loaded
   ! *********************************************************************
   subroutine condensed_nodes_sum_values(LoadCondensedNodes, ILump, NP, NListCondensedNodes, NNodesListCondensedNodes, &
                                         ListCondensedNodes, NODECODE, COEFDTempMem, GWCE_LV, ETA1, COEFDTemp)
      implicit none

      logical, intent(in) :: LoadCondensedNodes
      integer, intent(in) :: ILump
      integer, intent(in) :: NP
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(NListCondensedNodes)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: NODECODE(NP)
      real(8), intent(inout) :: COEFDTempMem(NP)
      real(8), intent(inout) :: GWCE_LV(NP)
      real(8), intent(inout) :: ETA1(NP)
      real(8), pointer, intent(inout) :: COEFDTemp(:)

      integer :: I, J, K, L

      if ((LoadCondensedNodes) .and. (ILump /= 0)) then
         do K = 1, NListCondensedNodes
            I = ListCondensedNodes(K, 1)
            if (I == 0) cycle
            if ((NODECODE(I) /= 0)) then
               ! 1) Sum them up
               do L = 2, NNodesListCondensedNodes(K)
                  J = ListCondensedNodes(K, L)
                  COEFDTemp(I) = COEFDTemp(I) + COEFDTemp(J)
                  GWCE_LV(I) = GWCE_LV(I) + GWCE_LV(J)
                  ETA1(I) = ETA1(I) + ETA1(J)
               end do
               ! 2) Distribute them
               ETA1(I) = ETA1(I)/NNodesListCondensedNodes(K)
               do L = 2, NNodesListCondensedNodes(K)
                  J = ListCondensedNodes(K, L)
                  COEFDTemp(J) = COEFDTemp(I)
                  GWCE_LV(J) = GWCE_LV(I)
                  ETA1(J) = ETA1(I)
               end do
            end if
         end do
      end if

   end subroutine condensed_nodes_sum_values

   ! *********************************************************************
   !> Copy the values from the back side to the front side of the vertical element wall (vew)
   !> 1d boundary nodes
   !>
   !> @param NFLUXIB64_GBL The number of global flux boundary nodes
   !> @param ILump The lumping flag
   !> @param NP The number of nodes
   !> @param NBOU The number of boundary nodes
   !> @param NVEL The number of velocity components
   !> @param NVELL The number of velocity components for each boundary node
   !> @param NBV The number of boundary nodes for each velocity component
   !> @param IBCONN The connectivity table for the boundary nodes
   !> @param LBCODEI The boundary code for each boundary node
   !> @param NODECODE The node code
   !> @param NListCondensedNodes The number of condensed nodes
   !> @param NNodesListCondensedNodes The number of nodes in the condensed nodes list
   !> @param ListCondensedNodes The list of condensed nodes
   !> @param ISSUBMERGED64 The submerged flag for the vew boundary nodes
   !> @param COEFD The gwce matrix diagonal entries
   !> @param COEFDTempMem The temporary memory for temporary diagonal entries
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary matrix diagonal entries
   ! *********************************************************************
   subroutine vew1d_copy_values_front_to_back(NFLUXIB64_GBL, ILump, NP, NBOU, NVEL, NVELL, NBV, IBCONN, LBCODEI, &
                                              NODECODE, NListCondensedNodes, NNodesListCondensedNodes, &
                                              ListCondensedNodes, ISSUBMERGED64, COEFDTempMem, GWCE_LV, &
                                              ETA1, COEFDTemp)
      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: ILump
      integer, intent(in) :: NP
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVEL
      integer, intent(in) :: NVELL(NBOU)
      integer, intent(in) :: NBV(NVEL)
      integer, intent(in) :: IBCONN(NVEL)
      integer, intent(in) :: LBCODEI(NVEL)
      integer, intent(in) :: NODECODE(NP)
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(NListCondensedNodes)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: ISSUBMERGED64(NVEL)
      real(8), intent(inout) :: COEFDTempMem(NP)
      real(8), intent(inout) :: GWCE_LV(NP)
      real(8), intent(inout) :: ETA1(NP)
      real(8), pointer, intent(inout) :: COEFDTemp(:)

      integer :: I, J, K, NNBB1, NNBB2

      if ((NFLUXIB64_GBL > 0) .and. (ILump /= 0)) then
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (64)
               do J = 1, NVELL(K)
                  I = I + 1
                  if (ISSUBMERGED64(I) /= 0) then
                     NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                     NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
                     if (NODECODE(NNBB1) /= 0 .and. NODECODE(NNBB2) /= 0) then
                        COEFDTemp(NNBB1) = COEFDTemp(NNBB2)
                        GWCE_LV(NNBB1) = GWCE_LV(NNBB2)
                        ETA1(NNBB1) = ETA1(NNBB2)
                     end if
                  end if
               end do
               I = I + NVELL(K)
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case DEFAULT
               I = I + NVELL(K)
            end select
         end do
      end if

   end subroutine vew1d_copy_values_front_to_back

   ! *********************************************************************
   !> Apply the verical element wall (vew) and condensed nodes boundary conditions
   !> to the load vector and eta1
   !>
   !> Nodal equations at IBTYPE=64 VEW boundary nodes and at condensed nodes are
   !> summed up together in the following part. First, the value on the front side,
   !> i.e., the floodplain side, is summed to the back side, i.e., the channel side.
   !> Secondly, the values at the condensed nodes which are most likely on channel
   !> beds are summed together. Notice that by this second step the nodes on channel
   !> bed hold the sum of all the values on the floodplain side and the grouped
   !> condensed nodes. And then finally, the value on the backside, i.e. on the
   !> channel bed, is copied to the front side, i.e., the floodplain side. Through
   !> this procedure, the values at the floodplain nodes and (condensed) channel nodes
   !> have the same values on both sides of the equations.
   !>
   !> @param NFLUXIB64_GBL The number of global flux boundary nodes
   !> @param ILump The lumping flag
   !> @param NP The number of nodes
   !> @param NBOU The number of boundary nodes
   !> @param NVEL The number of velocity components
   !> @param NVELL The number of velocity components for each boundary node
   !> @param NBV The number of boundary nodes for each velocity component
   !> @param IBCONN The connectivity table for the boundary nodes
   !> @param LBCODEI The boundary code for each boundary node
   !> @param NODECODE The node code
   !> @param NListCondensedNodes The number of condensed nodes
   !> @param NNodesListCondensedNodes The number of nodes in the condensed nodes list
   !> @param ListCondensedNodes The list of condensed nodes
   !> @param ISSUBMERGED64 The submerged flag for the vew boundary nodes
   !> @param LoadCondensedNodes Flag to indicate if the condensed nodes are loaded
   !> @param COEFD The gwce matrix diagonal entries
   !> @param COEFDTempMem The temporary memory for temporary diagonal entries
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary matrix diagonal entries
   ! *********************************************************************
   subroutine apply_vew1d_and_condensed_nodes(LoadCondensedNodes, ILump, NP, NFLUXIB64_GBL, NBOU, NVEL, &
                                              NVELL, NBV, IBCONN, LBCODEI, &
                                              NODECODE, NListCondensedNodes, NNodesListCondensedNodes, &
                                              ListCondensedNodes, IsSubmerged64, COEFD, &
                                              COEFDTemp, COEFDTempMem, GWCE_LV, ETA1)
#ifdef CMPI
      use global, only: dumy1
      use messenger, only: updater
#endif

      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: ILump
      integer, intent(in) :: NP
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVEL
      integer, intent(in) :: NVELL(NBOU)
      integer, intent(in) :: NBV(NVEL)
      integer, intent(in) :: IBCONN(NVEL)
      integer, intent(in) :: LBCODEI(NVEL)
      integer, intent(in) :: NODECODE(NP)
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(NListCondensedNodes)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: Issubmerged64(nvel)
      logical, intent(in) :: LoadCondensedNodes
      real(8), target, intent(inout) :: COEFD(NP)
      real(8), target, intent(inout) :: COEFDTempMem(NP)
      real(8), intent(inout) :: GWCE_LV(NP)
      real(8), intent(inout) :: ETA1(NP)
      real(8), pointer, intent(inout) :: COEFDTemp(:)

      integer :: I, J, K, L, NNBB1, NNBB2

      !.... Prep for the temporary LHS lumped array
      if (((NFLUXIB64_GBL > 0) .and. (ILump /= 0)) .or. LoadCondensedNodes) then
         COEFDTemp => COEFDTempMem
         COEFDTemp(:) = COEFD(:)
      else
         COEFDTemp => COEFD
      end if

      ! VEW: Sum front side values to back side
      call vew1d_sum_front_and_back_side(NFLUXIB64_GBL, ILump, NP, NBOU, NVEL, NVELL, NBV, IBCONN, LBCODEI, &
                                         NODECODE, NListCondensedNodes, NNodesListCondensedNodes, &
                                         ListCondensedNodes, ISSUBMERGED64, COEFD, COEFDTempMem, GWCE_LV, ETA1, COEFDTemp)

      !.... CONDENSED NODES: Summing up the values at condensed nodes
      call condensed_nodes_sum_values(LoadCondensedNodes, ILump, NP, NListCondensedNodes, NNodesListCondensedNodes, &
                                      ListCondensedNodes, NODECODE, COEFDTempMem, GWCE_LV, ETA1, COEFDTemp)

#ifdef CMPI
      if ((NFLUXIB64_GBL > 0 .or. LoadCondensedNodes) .and. (ILump /= 0)) then
         call UPDATER(COEFDTemp, GWCE_LV, DUMY1, 2)
      end if
#endif

      ! VEW: Copy values from back side to front side
      call vew1d_copy_values_front_to_back(NFLUXIB64_GBL, ILump, NP, NBOU, NVEL, NVELL, NBV, IBCONN, LBCODEI, &
                                           NODECODE, NListCondensedNodes, NNodesListCondensedNodes, &
                                           ListCondensedNodes, ISSUBMERGED64, COEFDTempMem, GWCE_LV, ETA1, COEFDTemp)

   end subroutine apply_vew1d_and_condensed_nodes

   ! *********************************************************************
   !> Update the right hand side sponge layer periodic boundary condition nodes
   !> for the consistent formulation
   !>
   !> @param NP The number of nodes
   !> @param MNEI The maximum number of neighbors
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param NNEIGH The number of neighbors
   !> @param EP The rms of the diagonal members of the gwce
   !> @param COEF The gwce matrix entries
   ! *********************************************************************
   subroutine update_coef_periodic_sponge_layer_consistent(NP, MNEI, NPERSEG, NNPERBC, IPERCONN, NNEIGH, EP, COEF)
      implicit none

      integer, intent(in) :: np
      integer, intent(in) :: mnei
      integer, intent(in) :: nperseg
      integer, intent(in) :: nnperbc
      integer, intent(in) :: iperconn(:, :)
      integer, intent(in) :: nneigh(np)
      real(8), intent(in) :: ep
      real(8), intent(inout) :: coef(np, mnei)

      integer :: I, J, I2

      !...  Modify equations associated with the slave nodes.
      !...  For each slave node, zero all off diagonal
      !......terms on the row and set diagnoal term to EP
      if (NPERSEG > 0) then
         do I = 1, NNPERBC
            I2 = IPERCONN(I, 2)
            Coef(I2, 1) = EP
            do J = 2, NNEIGH(I2)
               Coef(I2, J) = 0.0d0
            end do
         end do
      end if

   end subroutine update_coef_periodic_sponge_layer_consistent

   ! *********************************************************************
   !> Update the right hand side sponge layer periodic boundary condition nodes
   !> for the lumped formulation
   !>
   !> @param NP The number of nodes
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param EP The rms of the diagonal members of the gwce
   !> @param COEFD The coefficient for the elevation specified boundary nodes
   ! *********************************************************************
   subroutine update_coef_periodic_sponge_layer_lumped(NP, NPERSEG, NNPERBC, IPERCONN, EP, COEFD)
      implicit none

      integer, intent(in) :: np
      integer, intent(in) :: NPERSEG
      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      real(8), intent(in) :: EP
      real(8), intent(inout) :: COEFD(np)

      integer :: I, I2

      if (NPERSEG > 0) then
         do I = 1, NNPERBC
            I2 = IPERCONN(I, 2)
            Coefd(I2) = EP
         end do
      end if

   end subroutine update_coef_periodic_sponge_layer_lumped

   ! *********************************************************************
   !> Update eta for the sponge layer periodic boundary condition nodes
   !>
   !> @param NP The number of nodes
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param ETAS The eta solution
   ! *********************************************************************
   subroutine update_periodic_sponge_layer_nodes(NP, NPERSEG, NNPERBC, IPERCONN, ETAS, ETA2)
      implicit none

      integer, intent(in) :: np
      integer, intent(in) :: NPERSEG
      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      real(8), intent(inout) :: ETAS(np)
      real(8), intent(inout) :: ETA2(np)

      integer :: I, I1, I2

      if (NPERSEG > 0) then
         do I = 1, NNPERBC
            I1 = IPERCONN(I, 1)
            I2 = IPERCONN(I, 2)
            ETAS(I2) = ETAS(I1)
            ETA2(I2) = ETA2(I1)
         end do
      end if

   end subroutine update_periodic_sponge_layer_nodes

end module gwce_bc_forcing
