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
module mod_gwce_bc_forcing

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

      use mesh, only: nneigh, neitab
      use boundaries, only: neta, nbv, nbd, npebc, nope, nvel, bndlen2o3, lbcodei, &
                            nfluxb, nfluxf, nfluxib, nfluxgbc, nfluxrbc
      use global, only: per, amig, face, ff, efa, emo, nbfr, pr2, ilump, obccoef, coefd, &
                        usingDynamicWaterLevelCorrection, dynamicWaterLevelCorrection1, &
                        RES_BC_FLAG, BCFLAG_LNM, CBaroclinic, LNM_BC, FluxSettlingIT, &
                        nodecode, qn0, qn1, qn2, nperseg, iperconn, nnperbc
      use subdomain, only: subdomainOn, enforceBN, enforceECB
      use nodalattributes, only: geoidoffset, loadGeoidOffset, tau0var
      use spongelayer, only: no_met_in_sponge

      implicit none

      real(8), intent(in)    :: h1(:)
      real(8), intent(in)    :: elevdisc(:)
      real(8), intent(in)    :: en0(:), en1(:), en2(:)
      real(8), intent(in)    :: etiminc, TimeLoc, TimeH
      real(8), intent(in)    :: rampelev, rampmete, dt, ep
      real(8), intent(inout) :: esbin1(:)
      real(8), intent(inout) :: esbin2(:)
      real(8), intent(inout) :: etime1, etime2, etratio
      real(8), intent(inout) :: eta1(:), eta2(:), etas(:)
      real(8), intent(inout) :: gwce_lv(:)

      integer, intent(in)    :: IT
      logical, intent(in)    :: invertedBarometerOnElevationBoundary

      call apply_periodic_boundaries(NBFR, NETA, NBD, RampElev, TimeH, PER, &
                                     AMIG, FACE, FF, EFA, EMO, Eta2)
      call apply_aperiodic_boundaries(subdomainOn, enforceBN, NETA, NPEBC, NBD, RampElev, TimeLoc, &
                                      ETIME1, ETIME2, ETIMINC, ESBIN1, ESBIN2, ETRATIO, Eta2)
      call apply_geoid_offset(subdomainOn, LoadGeoidOffset, enforceBN, NETA, NBD, GeoidOffset, ETA2)
      call apply_inverted_barometer_on_elevation_boundaries(NO_MET_IN_SPONGE, &
                                                            invertedBarometerOnElevationBoundary, NETA, &
                                                            NBD, RampMete, PR2, ETA2)
      call apply_dynamic_water_level_correction(usingDynamicWaterLevelCorrection, neta, nbd, &
                                                dynamicWaterLevelCorrection1, eta2)
      call apply_levels_of_no_motion_boundary_conditions(RES_BC_FLAG, NOPE, BCFLAG_LNM, CBaroclinic, NETA, NBD, &
                                                         LNM_BC, ETA2)
      call apply_normal_flow_boundary_conditions(IT, NFLUXF, NFLUXB, NFLUXIB, NFLUXGBC, NFLUXRBC, &
                                                 NVEL, NBV, LBCODEI, NodeCode, FluxSettlingIT, &
                                                 QN0, QN1, QN2, ETAS, EN0, EN1, EN2, ETA1, ElevDisc, &
                                                 H1, TAU0Var, DT, BndLen2O3, GWCE_LV)
      call update_rhs_sponge_periodic_boundary_nodes(NPERSEG, NNPERBC, IPERCONN, NODECODE, EP, ETAS, GWCE_LV)
      call update_load_vector(ILump, NETA, NBD, NNEIGH, NEITAB, NODECODE, ETA1, ETA2, &
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
   !> @param NBOU The number of boundary nodes
   !> @param NVEL The number of velocity components
   !> @param NVELL The number of velocity components for each boundary node
   !> @param NBV The number of boundary nodes for each velocity component
   !> @param IBCONN The connectivity table for the boundary nodes
   !> @param LBCODEI The boundary code for each boundary node
   !> @param NODECODE The node code
   !> @param ISSUBMERGED64 The submerged flag for the vew boundary nodes
   !> @param COEFD The coefficient for the elevation specified boundary nodes
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary coefficient
   ! *********************************************************************
   subroutine vew1d_sum_front_and_back_side(NFLUXIB64_GBL, ILump, NBOU, NVELL, NBV, IBCONN, LBCODEI, &
                                            NODECODE, ISSUBMERGED64, COEFD, GWCE_LV, ETA1, &
                                            COEFDTemp)
      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: ILump
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: ISSUBMERGED64(:)
      real(8), intent(inout) :: COEFD(:)
      real(8), intent(inout) :: GWCE_LV(:)
      real(8), intent(inout) :: ETA1(:)
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
   subroutine condensed_nodes_sum_values(LoadCondensedNodes, ILump, NListCondensedNodes, NNodesListCondensedNodes, &
                                         ListCondensedNodes, NODECODE, GWCE_LV, ETA1, COEFDTemp)
      implicit none

      logical, intent(in) :: LoadCondensedNodes
      integer, intent(in) :: ILump
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(:)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: NODECODE(:)
      real(8), intent(inout) :: GWCE_LV(:)
      real(8), intent(inout) :: ETA1(:)
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
   !> @param NBOU The number of boundary nodes
   !> @param NVEL The number of velocity components
   !> @param NVELL The number of velocity components for each boundary node
   !> @param NBV The number of boundary nodes for each velocity component
   !> @param IBCONN The connectivity table for the boundary nodes
   !> @param LBCODEI The boundary code for each boundary node
   !> @param NODECODE The node code
   !> @param ISSUBMERGED64 The submerged flag for the vew boundary nodes
   !> @param COEFD The gwce matrix diagonal entries
   !> @param GWCE_LV The gwce load vector
   !> @param ETA1 The eta1 solution
   !> @param COEFDTemp The temporary matrix diagonal entries
   ! *********************************************************************
   subroutine vew1d_copy_values_front_to_back(NFLUXIB64_GBL, ILump, NBOU, NVELL, NBV, IBCONN, LBCODEI, &
                                              NODECODE, ISSUBMERGED64, GWCE_LV, ETA1, COEFDTemp)
      implicit none

      integer, intent(in) :: NFLUXIB64_GBL
      integer, intent(in) :: ILump
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: ISSUBMERGED64(:)
      real(8), intent(inout) :: GWCE_LV(:)
      real(8), intent(inout) :: ETA1(:)
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
   subroutine apply_vew1d_and_condensed_nodes(LoadCondensedNodes, ILump, NFLUXIB64_GBL, NBOU, &
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
      integer, intent(in) :: NBOU
      integer, intent(in) :: NVELL(:)
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: IBCONN(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NODECODE(:)
      integer, intent(in) :: NListCondensedNodes
      integer, intent(in) :: NNodesListCondensedNodes(:)
      integer, intent(in) :: ListCondensedNodes(:, :)
      integer, intent(in) :: Issubmerged64(:)
      logical, intent(in) :: LoadCondensedNodes
      real(8), target, intent(inout) :: COEFD(:)
      real(8), target, intent(inout) :: COEFDTempMem(:)
      real(8), intent(inout) :: GWCE_LV(:)
      real(8), intent(inout) :: ETA1(:)
      real(8), pointer, intent(inout) :: COEFDTemp(:)

      !.... Prep for the temporary LHS lumped array
      if (((NFLUXIB64_GBL > 0) .and. (ILump /= 0)) .or. LoadCondensedNodes) then
         COEFDTemp => COEFDTempMem
         COEFDTemp(:) = COEFD(:)
      else
         COEFDTemp => COEFD
      end if

      ! VEW: Sum front side values to back side
      call vew1d_sum_front_and_back_side(NFLUXIB64_GBL, ILump, NBOU, NVELL, NBV, IBCONN, LBCODEI, &
                                         NODECODE, ISSUBMERGED64, COEFD, GWCE_LV, ETA1, COEFDTemp)

      !.... CONDENSED NODES: Summing up the values at condensed nodes
      call condensed_nodes_sum_values(LoadCondensedNodes, ILump, NListCondensedNodes, NNodesListCondensedNodes, &
                                      ListCondensedNodes, NODECODE, GWCE_LV, ETA1, COEFDTemp)

#ifdef CMPI
      if ((NFLUXIB64_GBL > 0 .or. LoadCondensedNodes) .and. (ILump /= 0)) then
         call UPDATER(COEFDTemp, GWCE_LV, DUMY1, 2)
      end if
#endif

      ! VEW: Copy values from back side to front side
      call vew1d_copy_values_front_to_back(NFLUXIB64_GBL, ILump, NBOU, NVELL, NBV, IBCONN, LBCODEI, &
                                           NODECODE, ISSUBMERGED64, GWCE_LV, ETA1, COEFDTemp)

   end subroutine apply_vew1d_and_condensed_nodes

   ! *********************************************************************
   !> Update the right hand side sponge layer periodic boundary condition nodes
   !> for the consistent formulation
   !>
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param NNEIGH The number of neighbors
   !> @param EP The rms of the diagonal members of the gwce
   !> @param COEF The gwce matrix entries
   ! *********************************************************************
   subroutine update_coef_periodic_sponge_layer_consistent(NPERSEG, NNPERBC, IPERCONN, NNEIGH, EP, COEF)
      implicit none

      integer, intent(in) :: nperseg
      integer, intent(in) :: nnperbc
      integer, intent(in) :: iperconn(:, :)
      integer, intent(in) :: nneigh(:)
      real(8), intent(in) :: ep
      real(8), intent(inout) :: coef(:, :)

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
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param EP The rms of the diagonal members of the gwce
   !> @param COEFD The coefficient for the elevation specified boundary nodes
   ! *********************************************************************
   subroutine update_coef_periodic_sponge_layer_lumped(NPERSEG, NNPERBC, IPERCONN, EP, COEFD)
      implicit none

      integer, intent(in) :: NPERSEG
      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      real(8), intent(in) :: EP
      real(8), intent(inout) :: COEFD(:)

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
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The sponge layer periodic boundary connectivity table
   !> @param ETAS The eta solution
   ! *********************************************************************
   subroutine update_periodic_sponge_layer_nodes(NPERSEG, NNPERBC, IPERCONN, ETAS, ETA2)
      implicit none

      integer, intent(in) :: NPERSEG
      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      real(8), intent(inout) :: ETAS(:)
      real(8), intent(inout) :: ETA2(:)

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

   ! *********************************************************************
   !> Apply the boundary conditions to the elevation field
   !>
   !> @param NBFR The number of flow boundaries
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indices
   !> @param RampElev The ramp elevation
   !> @param TimeLoc The current time
   !> @param TimeH
   !> @param PER
   !> @param AMIG
   !> @param FACE
   !> @param FF
   !> @param EFA
   !> @param EMO
   !> @param Eta2 The elevation field
   ! *********************************************************************
   subroutine apply_periodic_boundaries(NBFR, NETA, NBD, RampElev, TimeH, PER, AMIG, FACE, FF, &
                                        EFA, EMO, Eta2)
      implicit none

      real(8), parameter  :: eps = epsilon(1.d0)

      integer, intent(in) :: NBFR
      integer, intent(in) :: NETA
      integer, intent(in) :: NBD(:)
      real(8), intent(in) :: RampElev
      real(8), intent(in) :: TimeH
      real(8), intent(in) :: PER(:)
      real(8), intent(in) :: AMIG(:)
      real(8), intent(in) :: FACE(:)
      real(8), intent(in) :: FF(:)
      real(8), intent(in) :: EFA(:, :)
      real(8), intent(in) :: EMO(:, :)
      real(8), intent(inout) :: Eta2(:)

      integer :: I, J, NCYC, NBDI
      real(8) :: ARGJ, ARG, RFF

      do J = 1, NBFR
         if (abs(PER(J)) <= eps) then
            NCYC = 0
         else
            NCYC = int(timeh/PER(J))
         end if
         ARGJ = AMIG(J)*(timeh - NCYC*PER(J)) + FACE(J)
         RFF = FF(J)*RampElev
         do I = 1, NETA
            ARG = ARGJ - EFA(J, I)
            NBDI = NBD(I)
            Eta2(NBDI) = Eta2(NBDI) + EMO(J, I)*RFF*cos(ARG)
         end do
      end do

   end subroutine apply_periodic_boundaries

   ! *********************************************************************
   !> Apply the boundary conditions to the elevation field
   !>
   !> @param subdomainOn Flag to indicate if the subdomain model is active
   !> @param enforceBN Flag to enforce the boundary nodes
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indices
   !> @param RampElev The ramp elevation
   !> @param TimeLoc The current time
   !> @param ETIME1 The start time for the elevation boundary condition
   !> @param ETIME2 The end time for the elevation boundary condition
   !> @param ETIMINC The time increment for the elevation boundary condition
   !> @param ESBIN1 The start elevation boundary condition
   !> @param ESBIN2 The end elevation boundary condition
   !> @param ETRATIO The ratio of the current time to the end time
   !> @param Eta2 The elevation field
   ! *********************************************************************
   subroutine apply_aperiodic_boundaries(subdomainOn, enforceBN, NETA, NPEBC, NBD, RampElev, TimeLoc, &
                                         ETIME1, ETIME2, ETIMINC, ESBIN1, ESBIN2, ETRATIO, Eta2)
      use subdomain, only: enforceEcb
      implicit none

      logical, intent(in) :: subdomainOn
      integer, intent(in) :: enforceBN
      integer, intent(in) :: NETA
      logical, intent(in) :: NPEBC
      integer, intent(in) :: NBD(:)
      real(8), intent(in) :: RampElev
      real(8), intent(in) :: TimeLoc
      real(8), intent(in) :: ETIMINC
      real(8), intent(inout) :: ETRATIO
      real(8), intent(inout) :: ESBIN1(:)
      real(8), intent(inout) :: ESBIN2(:)
      real(8), intent(inout) :: ETIME1
      real(8), intent(inout) :: ETIME2
      real(8), intent(inout) :: Eta2(:)
      integer :: I, J

      if (subdomainOn) then
         if (enforceBN == 1) call enforceEcb()
      else
         if (NPEBC) then
            if (TimeLoc > ETIME2) then
               ETIME1 = ETIME2
               ETIME2 = ETIME1 + ETIMINC
               ESBIN1 = ESBIN2
               do J = 1, NETA
                  read (19, *) ESBIN2(J)
               end do
            end if
            ETRATIO = (TimeLoc - ETIME1)/ETIMINC
            do I = 1, NETA
               Eta2(NBD(I)) = Eta2(NBD(I)) + RampElev*(ESBIN1(I) + ETRATIO*(ESBIN2(I) - ESBIN1(I)))
            end do
         end if
      end if
   end subroutine apply_aperiodic_boundaries

   ! *********************************************************************
   !> Apply the geoid offset to the boundary conditions
   !>
   !> Updates eta2 for the geoid offset along boundaries except when
   !> the subdomain model is active.
   !>
   !> @param subdomainOn Flag to indicate if the subdomain model is active
   !> @param LoadGeoidOffset Flag to indicate if the geoid offset is loaded
   !> @param enforceBN Flag to enforce the boundary nodes
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indicesq
   !> @param GeoidOffset The geoid offset
   !> @param ETA2 The elevation field
   ! *********************************************************************
   subroutine apply_geoid_offset(subdomainOn, LoadGeoidOffset, enforceBN, NETA, NBD, GeoidOffset, ETA2)
      implicit none
      logical, intent(in) :: subdomainOn
      logical, intent(in) :: LoadGeoidOffset
      integer, intent(in) :: enforceBN
      integer, intent(in) :: NETA
      integer, intent(in) :: NBD(:)
      real(8), intent(in) :: GeoidOffset(:)
      real(8), intent(inout) :: ETA2(:)
      integer :: i

      if (.not. subdomainOn .and. enforceBN /= 1 .and. LoadGeoidOffset) then
         do I = 1, NETA
            ETA2(NBD(I)) = ETA2(NBD(I)) + GeoidOffset(NBD(I))
         end do
      end if

   end subroutine apply_geoid_offset

   ! *********************************************************************
   !> Apply the inverted barometer boundary condition to the elevation
   !> boundaries so that low pressure systems can cross the boundary
   !> without creating an elevation anomaly.
   !>
   !> @param NO_MET_IN_SPONGE Flag to indicate if meteorological forcing is applied in the sponge layer
   !> @param invertedBarometerOnElevationBoundary Flag to apply the inverted barometer boundary condition
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indices
   !> @param RampMete The ramp meteorological value
   !> @param PR2 The pressure field at the current time step
   !> @param ETA2 The elevation field at the current time step
   !>
   ! *********************************************************************
   subroutine apply_inverted_barometer_on_elevation_boundaries(NO_MET_IN_SPONGE, invertedBarometerOnElevationBoundary, &
                                                               NETA, NBD, RampMete, PR2, ETA2)
      use adc_constants, only: G, RHOWAT0
      implicit none

      logical, intent(in) :: NO_MET_IN_SPONGE
      logical, intent(in) :: invertedBarometerOnElevationBoundary
      integer, intent(in) :: NETA
      integer, intent(in) :: NBD(:)
      real(8), intent(in) :: RampMete
      real(8), intent(in) :: PR2(:)
      real(8), intent(inout) :: ETA2(:)
      integer :: i

      if (.not. NO_MET_IN_SPONGE .and. invertedBarometerOnElevationBoundary) then
         do I = 1, NETA
            ETA2(NBD(I)) = ETA2(NBD(I)) + RampMete*(101300.d0/(RHOWAT0*G) - PR2(NBD(I)))
         end do
      end if

   end subroutine apply_inverted_barometer_on_elevation_boundaries

   ! *********************************************************************
   !> Apply the dynamic water level correction to the elevation boundaries
   !>
   !> @param usingDynamicWaterLevelCorrection Flag to apply the dynamic water level correction
   !> @param neta The number of elevation boundary nodes
   !> @param nbd The node boundary indices
   !> @param dynamicWaterLevelCorrection1 The dynamic water level correction at the current time
   !> @param eta2 The elevation field
   ! *********************************************************************
   subroutine apply_dynamic_water_level_correction(usingDynamicWaterLevelCorrection, neta, nbd, &
                                                   dynamicWaterLevelCorrection1, eta2)
      implicit none

      logical, intent(in) :: usingDynamicWaterLevelCorrection
      integer, intent(in) :: neta
      integer, intent(in) :: nbd(:)
      real(8), intent(in) :: dynamicWaterLevelCorrection1(:)
      real(8), intent(inout) :: eta2(:)
      integer :: i

      if (usingDynamicWaterLevelCorrection) then
         do I = 1, NETA
            ETA2(NBD(I)) = ETA2(NBD(I)) + dynamicWaterLevelCorrection1(nbd(i))
         end do
      end if
   end subroutine apply_dynamic_water_level_correction

   ! *********************************************************************
   !> Apply the levels of no motion boundary conditions
   !> These are considered the steric adjustments.
   !>
   !> @param RES_BC_FLAG The flag boundary condition
   !> @param NOPE The number of open boundaries
   !> @param BCFLAG_LNM The flag for the levels of no motion boundary conditions
   !> @param CBaroclinic The flag for the baroclinic model
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indices
   !> @param LNM_BC The levels of no motion boundary condition
   !> @param ETA2 The elevation field
   ! *********************************************************************
   subroutine apply_levels_of_no_motion_boundary_conditions(RES_BC_FLAG, NOPE, BCFLAG_LNM, CBaroclinic, NETA, NBD, &
                                                            LNM_BC, ETA2)
      implicit none

      integer, intent(in) :: res_bc_flag
      integer, intent(in) :: nope
      integer, intent(in) :: bcflag_lnm
      logical, intent(in) :: cbaroclinic
      integer, intent(in) :: neta
      integer, intent(in) :: nbd(:)
      real(8), intent(in) :: lnm_bc(:)
      real(8), intent(inout) :: eta2(:)
      integer :: i

      if (abs(RES_BC_FLAG) >= 1 .and. CBaroclinic .and. NOPE > 0 .and. BCFLAG_LNM > 0) then
         do I = 1, NETA
            ETA2(NBD(I)) = ETA2(NBD(I)) + LNM_BC(I)
         end do
      end if
   end subroutine apply_levels_of_no_motion_boundary_conditions

   ! *********************************************************************
   !> Compute the default inward flux boundary condition forcing
   !>
   !> @param qn0 The flow at the previous time step
   !> @param qn1 The flow at the current time step
   !> @param qn2 The flow at the next time step
   !> @param tau0 The tau0 value for the boundary node
   !> @param dt The time step
   ! *********************************************************************
   real(8) pure function compute_qforce_default(qn0, qn1, qn2, tau0, dt) result(qforce)

      implicit none

      real(8), intent(in) :: QN0
      real(8), intent(in) :: QN1
      real(8), intent(in) :: QN2
      real(8), intent(in) :: Tau0
      real(8), intent(in) :: DT

      real(8) :: dt2

      ! Specified flux boundary condition
      dt2 = dt*2.d0
      qforce = (QN2 - QN0)/dt2 + Tau0*QN1

   end function compute_qforce_default

   ! *********************************************************************
   !> Compute the qforcing for the time varying flux boundary type 30
   !>
   !> @param etas The solution to dE/dt from the prior solution
   !> @param qn1 The flow at the current time step
   !> @param h1 The total depth at the previous time step
   !> @param tau0 The tau0 value for the boundary node
   !> @param dt The time step
   ! *********************************************************************
   real(8) pure function compute_qforce_ibtype30(etas, qn1, h1, tau0, dt) result(qforce)
      use adc_constants, only: G
      implicit none

      real(8), intent(in) :: ETAS
      real(8), intent(in) :: QN1
      real(8), intent(in) :: H1
      real(8), intent(in) :: Tau0
      real(8), intent(in) :: DT

      real(8) :: celerity

      celerity = sqrt(G*H1)
      qforce = -CELERITY*ETAS/DT - Tau0*QN1
   end function compute_qforce_ibtype30

   ! *********************************************************************
   !> Compute the qforcing for time varying flux/elevation boundary type 32
   !>
   !> 2nd order extrapolation formula: \eta^{n+1} = 2*\eta^{n} - \eta^{n-1} is used
   !> and second order in time derivative is used
   !>
   !> @param etas The solution to dE/dt from the prior solution
   !> @param qn0 The flow at the previous time step
   !> @param qn1 The flow at the current time step
   !> @param qn2 The flow at the next time step
   !> @param en0 The elevation at the previous time step
   !> @param en1 The elevation at the current time step
   !> @param en2 The elevation at the next time step
   !> @param eta1 The elevation at the current time step
   !> @param h1 The total depth at the previous time step
   !> @param tau0 The tau0 value for the boundary node
   !> @param dt The time step
   ! *********************************************************************
   real(8) pure function compute_qforce_ibtype32(etas, qn0, qn1, qn2, en0, en1, en2, eta1, h1, tau0, dt) result(qforce)
      use adc_constants, only: G
      implicit none

      real(8), intent(in) :: ETAS
      real(8), intent(in) :: QN0
      real(8), intent(in) :: QN1
      real(8), intent(in) :: QN2
      real(8), intent(in) :: EN0
      real(8), intent(in) :: EN1
      real(8), intent(in) :: EN2
      real(8), intent(in) :: ETA1
      real(8), intent(in) :: H1
      real(8), intent(in) :: TAU0
      real(8), intent(in) :: DT

      real(8) :: celerity
      real(8) :: dt2

      dt2 = dt*2.d0
      celerity = sqrt(G*H1)

      qforce = (QN2 - QN0)/dt2 - CELERITY*(2.0d0*ETAS - (EN2 - EN0))/dt2 + &
               TAU0*(QN1 - CELERITY*(ETA1 - EN1))
   end function compute_qforce_ibtype32

   ! *********************************************************************
   !> Compute the qforcing for the outward flux boundary type 40, 41
   !>
   !> @param qn0 The flow at the previous time step
   !> @param qn1 The flow at the current time step
   !> @param tau0 The tau0 value for the boundary node
   !> @param dt The time step
   ! *********************************************************************
   real(8) pure function compute_qforce_ibtype4041(qn0, qn1, tau0, dt) result(qforce)
      implicit none

      real(8), intent(in) :: QN0
      real(8), intent(in) :: QN1
      real(8), intent(in) :: Tau0
      real(8), intent(in) :: DT

      ! Specified flux with flow out of the domain
      qforce = -(QN1 - QN0)/DT - Tau0*(QN1 + QN0)/2.0d0
   end function compute_qforce_ibtype4041

   ! *********************************************************************
   !> Compute the qforcing for the combined radiation and specified flux
   !> boundary condition
   !>
   !> @param qn0 The flow at the previous time step
   !> @param qn1 The flow at the current time step
   !> @param qn2 The flow at the next time step
   !> @param etas The solution to dE/dt from the prior solution
   !> @param en0 The elevation at the previous time step
   !> @param en1 The elevation at the current time step
   !> @param en2 The elevation at the next time step
   !> @param eta1 The elevation at the current time step
   !> @param elevdisc The elevation at the discharge boundaries
   !> @param h1 The total depth at the previous time step
   !> @param tau0 The tau0 value for the boundary node
   !> @param it The current time step
   !> @param dt The time step
   !> @param FluxSettlingIT The number of iterations of flux settling
   ! *********************************************************************
   real(8) pure function compute_qforce_ibtype52(qn0, qn1, qn2, etas, eta1, &
                                                 elevdisc, h1, tau0, &
                                                 it, dt, FluxSettlingIT) result(qforce)
      use adc_constants, only: G
      implicit none

      real(8), intent(in) :: QN0
      real(8), intent(in) :: QN1
      real(8), intent(in) :: QN2
      real(8), intent(in) :: EtaS
      real(8), intent(in) :: Eta1
      real(8), intent(in) :: ElevDisc
      real(8), intent(in) :: H1
      real(8), intent(in) :: Tau0
      real(8), intent(in) :: DT
      integer, intent(in) :: IT
      integer, intent(in) :: FluxSettlingIT

      real(8) :: celerity

      ! Specified flux boundary condition with radiation boundary condition
      qforce = compute_qforce_default(QN0, QN1, QN2, Tau0, DT)
      if (IT > FluxSettlingIT) then
         celerity = sqrt(G*H1)
         qforce = qforce - Celerity*(EtaS/DT + Tau0*(Eta1 - ElevDisc))
      end if
   end function compute_qforce_ibtype52

   ! *********************************************************************
   !> Compute the normal flow boundary condition force at the given boundary index
   !>
   !> Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
   !> assume that QN is positive into the domain.  QFORCEJ has a -1
   !> built in and the terms are not explicitly negated. Boundary
   !> conditions using computed fluxes (LBCODEI 30, 40) compute a normal
   !> flux that  is positive out of the domain.  Therefore, to match
   !> the formulation these terms must be explicitly multiplied by -1.
   !>
   !> Note 3, Eta1 is the latest computed elevation (it was updated previously).
   !>
   !> @param LBCode The boundary code for the boundary
   !> @param IT The current time step
   !> @param QN0 The flow at the previous time step
   !> @param QN1 The flow at the current time step
   !> @param QN2 The flow at the next time step
   !> @param ETAS THe solution to dE/dt from the prior solution
   !> @param EN0
   !> @param EN1
   !> @param EN2
   !> @param ETA1 Elevation at the current time step
   !> @param Tau0 Tau0 for the boundary node
   !> @param dt The time step
   ! *********************************************************************
   real(8) pure function compute_qforce_normal_flow_boundary(LBCode, DT, IT, FluxSettlingIT, QN0, QN1, QN2, ETAS, EN0, &
                                                             EN1, EN2, ETA1, ElevDisc, H1, TAU0) result(qforce)
      implicit none

      integer, intent(in) :: IT
      integer, intent(in) :: LBCode
      integer, intent(in) :: FluxSettlingIT
      real(8), intent(in) :: QN0
      real(8), intent(in) :: QN1
      real(8), intent(in) :: QN2
      real(8), intent(in) :: ETAS
      real(8), intent(in) :: EN0
      real(8), intent(in) :: EN1
      real(8), intent(in) :: EN2
      real(8), intent(in) :: ETA1
      real(8), intent(in) :: ElevDisc
      real(8), intent(in) :: h1
      real(8), intent(in) :: Tau0
      real(8), intent(in) :: dt

      select case (LBCode)
      case (30)
         qforce = compute_qforce_ibtype30(ETAS, QN1, H1, TAU0, DT)
      case (32)
         qforce = compute_qforce_ibtype32(ETAS, QN0, QN1, QN2, EN0, &
                                          EN1, EN2, ETA1, H1, TAU0, DT)
      case (40, 41)
         qforce = compute_qforce_ibtype4041(QN0, QN1, TAU0, DT)
      case (52)
         qforce = compute_qforce_ibtype52(QN0, QN1, QN2, ETAS, ETA1, ElevDisc, &
                                          H1, TAU0, IT, DT, FluxSettlingIT)
      case default
         qforce = compute_qforce_default(QN0, QN1, QN2, TAU0, DT)
      end select

   end function compute_qforce_normal_flow_boundary

   ! *********************************************************************
   !> Apply the normal flow boundary conditions to the load vector
   !> Impose normal flow, radiation or gradient boundary conditions
   !> along flow boundary to load vector GWCE_LV(I)
   !>
   !> @param IT The current time step
   !> @param NFLUXF The number of flow boundaries
   !> @param NFLUXB The number of boundaries
   !> @param NFLUXIB The number of internal boundaries
   !> @param NFLUXGBC The number of flux
   !> @param NFLUXRBC The number of radiation boundaries
   !> @param NVEL The number of elevation boundary nodes
   !> @param NBV The node boundary indices
   !> @param LBCODEI The boundary code
   !> @param NodeCode The wet/dry node code
   !> @param FluxSettlingIT The number of iterations of flux settling
   !> @param QN0 The flow at the previous time step
   !> @param QN1 The flow at the current time step
   !> @param QN2 The flow at the next time step
   !> @param ETAS The solution to dE/dt from the prior solution
   !> @param EN0
   !> @param EN1
   !> @param EN2
   !> @param ETA1 Elevation at the current time step
   !> @param ElevDisc Elevation at discharge boundaries
   !>
   ! *********************************************************************
   subroutine apply_normal_flow_boundary_conditions(IT, NFLUXF, NFLUXB, NFLUXIB, NFLUXGBC, NFLUXRBC, &
                                                    NVEL, NBV, LBCODEI, NodeCode, FluxSettlingIT, &
                                                    QN0, QN1, QN2, ETAS, EN0, EN1, EN2, ETA1, ElevDisc, &
                                                    H1, TAU0, DT, BndLen2O3, GWCE_LV)
      implicit none

      integer, intent(in) :: IT
      integer, intent(in) :: NFLUXF
      integer, intent(in) :: NFLUXB
      integer, intent(in) :: NFLUXIB
      integer, intent(in) :: NFLUXGBC
      integer, intent(in) :: NFLUXRBC
      integer, intent(in) :: NVEL
      integer, intent(in) :: NBV(:)
      integer, intent(in) :: LBCODEI(:)
      integer, intent(in) :: NodeCode(:)
      integer, intent(in) :: FluxSettlingIT
      real(8), intent(in) :: QN0(:)
      real(8), intent(in) :: QN1(:)
      real(8), intent(in) :: QN2(:)
      real(8), intent(in) :: ETAS(:)
      real(8), intent(in) :: EN0(:)
      real(8), intent(in) :: EN1(:)
      real(8), intent(in) :: EN2(:)
      real(8), intent(in) :: ETA1(:)
      real(8), intent(in) :: ElevDisc(:)
      real(8), intent(in) :: H1(:)
      real(8), intent(in) :: TAU0(:)
      real(8), intent(in) :: DT
      real(8), intent(in) :: BndLen2O3(:)
      real(8), intent(inout) :: GWCE_LV(:)

      real(8) :: qforcei, qforcej
      real(8) :: bndleno6nc
      integer :: j, current_node, previous_node
      integer :: nodecode_prev_node, nodecode_current_node

      if ((NFLUXF == 0) .and. (NFLUXB == 0) .and. (NFLUXIB == 0) .and. (NFLUXGBC == 0) .and. (NFLUXRBC == 0)) then
         return
      end if

      qforcej = compute_qforce_normal_flow_boundary(LBCodeI(1), DT, IT, FluxSettlingIT, QN0(1), QN1(1), QN2(1), &
                                                    ETAS(NBV(1)), EN0(1), EN1(1), EN2(1), ETA1(NBV(1)), &
                                                    ElevDisc(1), H1(NBV(1)), TAU0(NBV(1)))

      do J = 2, NVEL
         ! Get the current and previous node id
         current_node = NBV(J)
         previous_node = NBV(J - 1)

         ! Swap the qforce values and compute the new qforce value at the next node
         qforcei = qforcej
         qforcej = compute_qforce_normal_flow_boundary(LBCODEI(J), DT, IT, FluxSettlingIT, QN0(J), QN1(J), QN2(J), &
                                                       ETAS(current_node), EN0(J), EN1(J), EN2(J), &
                                                       ETA1(current_node), ElevDisc(J), &
                                                       H1(current_node), TAU0(current_node))

         ! Compute the nodecode for the two nodes in the boundary
         nodecode_prev_node = NodeCode(previous_node)
         nodecode_current_node = NodeCode(current_node)

         ! Compute the length of the boundary segment and zero for non-wet nodes
         BndLenO6NC = nodecode_prev_node*nodecode_current_node*BndLen2O3(J - 1)/4.0d0

         ! Apply the boundary conditions to the load vector
         GWCE_LV(previous_node) = GWCE_LV(previous_node) + BndLenO6NC*(2.0d0*QForceI + QForceJ)
         GWCE_LV(current_node) = GWCE_LV(current_node) + BndLenO6NC*(2.d0*QForceJ + QForceI)
      end do
   end subroutine apply_normal_flow_boundary_conditions

   ! *********************************************************************
   !> Update the right hand side sponge layer periodic boundary condition nodes
   !>
   !> @param NPERSEG The number of periodic segments
   !> @param NNPERBC The number of nodes per boundary condition
   !> @param IPERCONN The connectivity table
   !> @param NODECODE The wet/dry node code
   !> @param EP The rms of the diagonal members of the gwce
   !> @param ETAS The eta solution
   !> @param GWCE_LV The gwce load vector
   subroutine update_rhs_sponge_periodic_boundary_nodes(NPERSEG, NNPERBC, &
                                                        IPERCONN, NODECODE, EP, &
                                                        ETAS, GWCE_LV)
      implicit none

      integer, intent(in) :: NPERSEG
      integer, intent(in) :: NNPERBC
      integer, intent(in) :: IPERCONN(:, :)
      integer, intent(in) :: NODECODE(:)
      real(8), intent(in) :: EP
      real(8), intent(inout) :: ETAS(:)
      real(8), intent(inout) :: GWCE_LV(:)

      integer :: I, I2

      if (NPERSEG > 0) then
         do I = 1, NNPERBC
            I2 = IPERCONN(I, 2)
            ETAS(I2) = 0.0d0 !  ETA1 ==> ETA2
            GWCE_LV(I2) = ETAS(I2)*NODECODE(I2)*EP
         end do
      end if

   end subroutine update_rhs_sponge_periodic_boundary_nodes

   ! *********************************************************************
   !> Impost the elevation boundary condition to the gwce load vector
   !> Note: EP is the rms of the diagonal members of the gwce. It is used
   !> to scale the diagonal element for the elevation specified boundary
   !> nodes and threfore must also be used to scale the RHS of the equations
   !>
   !> @param ILump Flag to indicate if the gwce is lumped
   !> @param NETA The number of elevation boundary nodes
   !> @param NBD The node boundary indices
   !> @param NNEIGH The number of neighbors
   !> @param NEITAB The neighbor table
   !> @param NODECODE The wet/dry node code
   !> @param ETA1 Elevation at the previous time step
   !> @param ETA2 Elevation at the current time step
   !> @param OBCCOEF The obc coefficient
   !> @param COEFD The coefficient for the elevation specified boundary nodes
   !> @param EP The rms of the diagonal members of the gwce
   !> @param ETAS The eta solution
   !> @param GWCE_LV The gwce load vector
   !>
   ! *********************************************************************
   subroutine update_load_vector(ILump, NETA, NBD, NNEIGH, NEITAB, NODECODE, ETA1, ETA2, &
                                 OBCCOEF, COEFD, EP, ETAS, GWCE_LV)
      implicit none

      integer, intent(in) :: ilump
      integer, intent(in) :: neta
      integer, intent(in) :: nbd(:)
      integer, intent(in) :: nneigh(:)
      integer, intent(in) :: nodecode(:)
      integer, intent(in) :: neitab(:, :)
      real(8), intent(in) :: eta1(:)
      real(8), intent(in) :: eta2(:)
      real(8), intent(in) :: obccoef(:, :)
      real(8), intent(in) :: coefd(:)
      real(8), intent(in) :: ep
      real(8), intent(inout) :: etas(:)
      real(8), intent(inout) :: gwce_lv(:)

      integer :: i, j, nbdi

      if (ILump == 0) then
         do I = 1, NETA
            NBDI = NBD(I)
            ETAS(NBDI) = ETA2(NBDI) - ETA1(NBDI)
            GWCE_LV(NBDI) = ETAS(NBDI)*NODECODE(NBDI)*EP
            do J = 2, NNEIGH(NBDI)
               GWCE_LV(NEITAB(NBDI, J)) = GWCE_LV(NEITAB(NBDI, J)) - ETAS(NBDI)*OBCCOEF(I, J - 1)
            end do
         end do
      else
         do I = 1, NETA
            NBDI = NBD(I)
            ETAS(NBDI) = ETA2(NBDI) - ETA1(NBDI)
            GWCE_LV(NBDI) = ETAS(NBDI)*NODECODE(NBDI)*COEFD(NBDI)
         end do
      end if
   end subroutine update_load_vector

end module mod_gwce_bc_forcing
