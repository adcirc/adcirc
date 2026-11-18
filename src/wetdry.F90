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
!***********************************************************************
!                     M O D U L E   W E T   D R Y
!***********************************************************************
!     Executes wetting and drying algorithm.
!***********************************************************************
! Logical Variable List (default value .FALSE., set in global.f)       *
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     C2DDI            - 2D Depth Integrated model run                 *
!     C3D              - 3D model run                                  *
!     C3DDSS           - Stress form of the 3D momentum equations      *
!     C3DVS            - Velocity form of the 3D momentum equations    *
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!   See header.f for a summary history of code modifications.          *
!***********************************************************************
!WET...
!WET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
!WET...
!WET...NOTE:NNODECODE is a working variable that can change more than once
!WET...               during a time step
!WET...     NNODECODE = 0 for a dry node
!WET...     NNODECODE = 1 for a wet node
!WET...     NODECODE  - is a more static version of NNODECODE that is reconciled
!WET...                 once and for all at the end time step
!WET...
!WET...
!WET...        (   DRYING CRITERIA   )
!WET...
!WET...A node should be dry under two conditions.
!WET...D1.) If the total water depth falls below H0.
!WET .......Note: if the total water depth falls below 0.8*H0, the surface elevation
!WET........is lifted up so that the total water depth = 0.8*H0.
!WET......
!WET...D2.) If the node is connected to only nonfunctioning (dry) elements.  In
!WET........this case the node is dried due to becoming landlocked.
!WET........Note: this criteria is applied after all other wetting and drying criteria
!WET...
!WET...An element should be dry under the following conditions.
!WET...DE3.) This is an elemental check section designed to avoid artificial wetting of
!WET.........of control sections
!WET.........All elements where downhill flow originates from a barely wet node
!WET.........(defined as 1.2*H0) into wet nodes are forced inactive; the only exception
!WET......... is receiving overtopped barrier nodes
!WET...
!WET...        (   WETTING CRITERIA   )
!WET...
!WET...A node should be wet under two conditions.
!WET...W1.) If 2 nodes on an element are wet and one is dry, wet the dry node
!WET........if the water level at one of the wet nodes is greater than the
!WET........water level at the dry node and the steady state velocity that
!WET........would result from a balance between the water level gradient and
!WET........bottom friction would yield a velocity > VELMIN.
!WET........Note that the criteria outlined in DE3 must also be satified before
!WET........the node is allowed to wet
!WET...
!WET...W2.) If an element has a node lying on a receiving internal barrier boundary or
!WET......specified discharge boundary that is actively discharging flow into the
!WET......domain at that node, all nodes in this element must stay wet.
!WET...
!WET...
!WET...        (  VELOCITY BOUNDARY CONDITION  )
!WET...
!WET...Either a natural or essential boundary condition can be used as a velocity
!WET...boundary condition in the momentum equation solution along a wet/dry boudary
!WET...To use a natural boundary condition, do nothing along the wet/dry interface.
!WET...To use an essential, no velocity boundary condition, identify the nodes along
!WET...the wet/dry interface and zero out the velocity at the nodes.  Interface nodes
!WET...can easily be identified by comparing the number of active elements a node is
!WET...connected to (MJU) to the total number of elements a node is connected to (NODELE).
!WET...If MJU < NODELE for any node, it must lie along the wet/dry interface.  See
!WET...further comments at the end of the momentum equation solution section.
!WET...
module mod_wetdry

   implicit none

   private

   public :: initializeWettingAndDrying, computeWettingAndDrying, &
             noffActive, velmin

   real(8) :: habsmin
   real(8) :: hoff
   real(8) :: velmin

   integer, allocatable :: nibcnt(:)
   integer, allocatable :: noffold(:)
   integer, allocatable :: noff_sg(:)

   ! jgf52.08.08: Enable analyst to remove NOFF from wetting and
   ! drying algorithm; .true. by default to conform with prior
   ! ADCIRC versions.
   logical :: noffActive = .true.

contains

   !----------------------------------------------------------------------
   !                  S U B R O U T I N E
   !     I N I T I A L I Z E   W E T T I N G   A N D   D R Y I N G
   !----------------------------------------------------------------------
   !     Allocates arrays and initializes constants.
   !----------------------------------------------------------------------
   subroutine initializeWettingAndDrying()
      use sizes, only: mne, mnp
      use global, only: h0, nodecode, nnodecode
      use nodalattributes, only: loadSubgridBarrier
      implicit none

      allocate (nibcnt(mnp))
      allocate (noffold(mne))
      if (loadSubgridBarrier) allocate (noff_sg(mne))
      nnodecode = 1
      nodecode = 1
      noffold(:) = 1
      habsmin = 0.8d0*h0
      hoff = 1.2d0*h0

      !----------------------------------------------------------------------
   end subroutine initializeWettingAndDrying
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   !                   S U B R O U T I N E
   !     C O M P U T E   W E T T I N G   A N D   D R Y I N G
   !----------------------------------------------------------------------
   !     Determines which nodes should be wet and which should be dry.
   !----------------------------------------------------------------------
   subroutine computeWettingAndDrying(it)
      use ADC_CONSTANTS, only: g
      use global, only: noff, nodecode, nnodecode, eta2, tk, &
                        nolifa, &
                        bsx1, bsy1, btime_end, C2DDI, C3D, h0, ifnlfa, nddt, dtdp, &
                        nibnodecode, ncchange, tkm, directvelWD, useHF, &
                        cbaroclinic, iden
#ifdef CMPI
      use global, only: ilump, idumy, dumy1, dumy2
      use messenger, only: updatei, updater, wetDrySum
#endif
      use mesh, only: ne, np, nm, dp, mju, totalArea, nm, x, y, areas, &
                      LBArray_Pointer, nneighele, neitabele, onIBTYP64Floodplain
      use nodalattributes, only: BFCdLLimit, fgamma, ftheta, fric, &
                                 manningsn, hbreak, ifhybf, ifnlbf, iflinbf, loadManningsN, &
                                 loadSubgridBarrier, &
                                 loadCondensedNodes, &
                                 ListCondensedNodes, NListCondensedNodes, &
                                 NNodesListCondensedNodes
      use BOUNDARIES, only: IBCONN, ISSUBMERGED64, NFLUXIB64_GBL, &
                            LBCODEI
      use global_3dvs, only: islip, sigma, evtot, q, nfen, sal, temp
      use subdomain, only: subdomainOn, enforceBN, enforceWDcb, enforceWDob
      use subgrid, only: subgrid_level0, subgridVertList, gridDepthVertETA2, &
                         subgrid_level1, cmfVertETA2, cfVertETA2, wetDepthVertETA2, &
                         wetFracVertTab, gridDepthVertTab, numPhi
      implicit none
      integer, intent(in) :: it ! time step number
      complex(8) :: duds
      integer :: nc1, nc2, nc3
      integer :: nm1, nm2, nm3, nmi
      integer :: nm123
      integer :: ncele
      integer :: nctot
      real(8) :: kslip
      real(8) :: vel
      real(8) :: z0b1
      real(8) :: areaEle
      real(8) :: tkWet
      real(8) :: etaN1, etaN2, etaN3
      real(8) :: hTotN1, hTotN2, hTotN3
      real(8) :: deldist, deleta
      real(8) :: htot
      real(8) :: h1
      real(8) :: tk_tmp(np)
#ifdef CMPI
      real(8) :: tmp1(np), tmp2(np), tmp3(np)
#endif
      real(8) :: totArea1, totArea2, eta2Avg
      integer :: tk_cnt(np)
      integer :: nbnctot
      integer :: i, j, k, l
      integer :: ie
      real(8) :: tknf !dmw 202207, added for "new" formula for vel
      integer :: nnc
      integer :: nnbb2
      integer :: nwetele
      integer :: numGreater

      if (nolifa /= 2) then
         return ! wetting and drying is not active
      end if

      ! VEL is only assigned under C2DDI or C3D; the compiler cannot see that
      ! exactly one of them is always true.
      vel = 0.d0

      !WET...
      !WET...WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP
      !WET...
      TK_TMP(:) = 0.d0
      TK_CNT(:) = 0
      if (directvelWD .and. useHF) then
         TKNF = 100.d0
      else
         TKNF = 1.d0
      end if

      do I = 1, NP
         NIBCNT(I) = 0
      end do
      do I = 1, NE
         NOFFOLD(I) = NOFF(I)
         NOFF(I) = 1
      end do

      if (loadSubgridBarrier) call applySubgridBarrier()

      !WET...
      !WET...WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
      !WET....Drying Criteria D1: this depends on NODECODE and updates NODECODE
      !WET...
      !WET...
      do I = 1, NP
         if (NODECODE(I) == 1) then
            if (subgrid_level0) then
               HTOT = gridDepthVertETA2(i)
            else
               HTOT = DP(I) + ETA2(I)
            end if
            if (HTOT <= H0) then
               if (subgrid_level0) then
                  !JLW: need to interpolate here to set eta2
                  if (HTOT < HABSMIN) then
                     gridDepthVertETA2(i) = HABSMIN
                     wetDepthVertETA2(i) = HABSMIN
                     if (subgridVertList(i) == 1) then
                        numGreater = 0
                        do j = 1, numPhi
                           if (gridDepthVertETA2(i) < &
                               gridDepthVertTab(i, j)) then
                              numGreater = numGreater + 1
                           end if
                        end do
                        if (numGreater == numPhi) then
                           eta2(i) = wetFracVertTab(i, 1)
                        elseif (numGreater == 0) then
                           eta2(i) = wetFracVertTab(i, numPhi) &
                                     + (HABSMIN - gridDepthVertTab(i, numPhi))
                        else
                           !JLW: now back calculate eta2 based on this change
                           eta2(i) = ((gridDepthVertETA2(i) &
                                       - gridDepthVertTab(i, numPhi - numGreater)) &
                                      /(gridDepthVertTab(i, numPhi - numGreater + 1) &
                                        - gridDepthVertTab(i, numPhi - numGreater)) &
                                      *(wetFracVertTab(i, numPhi - numGreater + 1) &
                                        - wetFracVertTab(i, numPhi - numGreater)) &
                                      + wetFracVertTab(i, numPhi - numGreater))
                        end if
                     else
                        eta2(i) = HABSMIN - DP(I)
                     end if
                  end if
                  NNODECODE(I) = 0
                  NODECODE(I) = 0
                  NCCHANGE = NCCHANGE + 1
               else
                  if (HTOT < HABSMIN) ETA2(I) = HABSMIN - DP(I)
                  NNODECODE(I) = 0
                  NODECODE(I) = 0
                  NCCHANGE = NCCHANGE + 1 !NCCHANGE=0 set near beginning of GWCE
               end if
            end if
         end if
      end do
      !WET...
      !WET...END WET/DRY SECTION - PART 1
      !WET...

      !.....
      !.....Nodal Wetting Criteria for the ibtype=64 weir boundary
      !.....Check the node if it is along the ibtype=64 weir boundary and if the weir is
      !.....submerged. If it is the case, set the nnodecode to 1.   SB
      !.....
      if (NFLUXIB64_GBL > 0) then
         do I = 1, NP
            J = LBArray_Pointer(I)
            if (J > 0) then
               if ((LBCODEI(J)) == 64 .and. (NODECODE(I) == 0)) then
                  if ((ISSUBMERGED64(J) > 0)) then
                     NNBB2 = IBCONN(J)
                     if (NNBB2 > 0) then ! This check is needed ocasionally but not always when the barrier is splitted for paralell computation
                        ! ETA2(I) = ETA2(NNBB2)
                        ! HTOT=DP(I)+ETA2(I)
                        HTOT = DP(I) + max(ETA2(I), ETA2(NNBB2))
                        if (HTOT > HOFF) then
                           ! ETA2(I) = MAX(ETA2(I), ETA2(NNBB2))
                           TotArea1 = computeTotAreaAtNode(I)
                           TotArea2 = computeTotAreaAtNode(NNBB2)
                           ETA2AVG = (ETA2(I)*TotArea1 + ETA2(NNBB2)*TotArea2)/(TotArea1 + TotArea2)
                           ETA2(I) = ETA2AVG
                           ETA2(NNBB2) = ETA2AVG
                           NNODECODE(I) = 1
                           NODECODE(I) = 1
                           NCCHANGE = NCCHANGE + 1 !NCCHANGE=0 set near beginning of GWCE
                        end if
                     end if
                  end if
               end if
            end if
         end do
      end if

      !jjwC     Use Message-Passing to update nodecode and nnodecode at ghost nodes
      !jjw#ifdef CMPI
      !jjw         CALL UPDATEI(NODECODE,NNODECODE,2)
      !jjw#endif

      !WET...
      !WET...WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
      !WET...
      do I = 1, NE
         NM1 = NM(I, 1)
         NM2 = NM(I, 2)
         NM3 = NM(I, 3)

         !WET...
         !WET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
         !WET...
         NCTOT = NODECODE(NM1) + NODECODE(NM2) + NODECODE(NM3)
         if (NCTOT == 2) then
            ETAN1 = ETA2(NM1)
            ETAN2 = ETA2(NM2)
            ETAN3 = ETA2(NM3)
            if (subgrid_level0) then
               HTOTN1 = gridDepthVertETA2(NM1)
               HTOTN2 = gridDepthVertETA2(NM2)
               HTOTN3 = gridDepthVertETA2(NM3)
            else
               HTOTN1 = DP(NM1) + ETA2(NM1)
               HTOTN2 = DP(NM2) + ETA2(NM2)
               HTOTN3 = DP(NM3) + ETA2(NM3)
            end if
            NMI = 0
            if ((NODECODE(NM1) == 1) .and. (NODECODE(NM2) == 1)) then
               if ((HTOTN1 >= HOFF) .and. (HTOTN2 >= HOFF)) then
                  NM123 = NM1
                  if (ETA2(NM2) > ETA2(NM1)) NM123 = NM2
                  DELDIST = sqrt((y(NM3) - y(NM123))**2.d0 &
                                 + (X(NM3) - X(NM123))**2.d0)
                  DELETA = ETA2(NM123) - ETA2(NM3)
                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                  if (DELETA < 0.d0) DELETA = 0.d0
                  if (subgrid_level0) then
                     H1 = wetDepthVertETA2(NM123)
                  else
                     H1 = ETA2(NM123) + DP(NM123)
                  end if
                  !. RJW merged from Casey 071219: Added the following logic for 3D friction.
                  !. RJW modified the following for 3D friction
                  if (C2DDI) then
                     !sb46.28sb02
                     !<<                     Convert Manning's N to Cd, if necessary.
                     if (LoadManningsN) then
                        if (subgrid_level0) then
                           if (subgridVertList(NM123) == 1) then
                              if (subgrid_level1) then
                                 FRIC(NM123) = cmfVertETA2(NM123)
                              else
                                 FRIC(NM123) = cfVertETA2(NM123)
                              end if
                           else
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                           end if
                        else
                           FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                         /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                           **(1.d0/3.d0))
                        end if
                        if (FRIC(NM123) < BFCdLLimit) then
                           FRIC(NM123) = BFCdLLimit
                        end if
                     end if
                     !-- DMW Reformulate Velocity by taking sqrt of u instead of using velmin 202207
                     if (directvelWD .and. useHF) then
                        VEL = sqrt(g*H1*DELETA/DELDIST)
                     elseif (directvelWD) then
                        VEL = sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                     else
                        TKWET = max(FRIC(NM123)*(dble(IFLINBF) + (VELMIN/H1)* &
                                                 (dble(IFNLBF) + dble(IFHYBF)* &
                                                  (1.d0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA))), &
                                    0.0001d0)
                        VEL = G*(DELETA/DELDIST)/TKWET
                     end if
                  elseif (C3D) then
                     ! solve for the depth averaged velocity,U, from the relation :
                     !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                     !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                     ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                     Z0B1 = bottomRoughnessLength(NM123)
                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                           *((H1 + Z0B1)*log((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41d0)
                  end if

                  if (VEL > VELMIN) then
                     !    ....         third node met criteria and is also wet
                     NNODECODE(NM3) = 1
                     NMI = NM3
                     !. RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                     if (C2DDI) then

                        TK_TMP(NM123) = min( &
                                        FRIC(NM123)*(dble(IFLINBF) + (TKNF*VEL/H1)* &
                                                     (dble(IFNLBF) + dble(IFHYBF)* &
                                                      (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                      (FGAMMA/FTHETA))), &
                                        2.0d0/DTDP)
                        TK(NM123) = max(TK_TMP(NM123), TK(NM123))
                        !WJP needed for internal_tide stuff
                        TKM(1:2, NM123) = TK(NM123)
                        if (NM123 == NM2) then
                           TK_TMP(NM1) = TK(NM123)
                           TK(NM1) = max(TK_TMP(NM1), TK(NM1))
                           TKM(1:2, NM1) = TK(NM1)
                           if (TK_TMP(NM1) < TK(NM1)) then
                              TK(NM123) = TK(NM1)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        else
                           TK_TMP(NM2) = TK(NM123)
                           TK(NM2) = max(TK_TMP(NM2), TK(NM2))
                           TKM(1:2, NM2) = TK(NM2)
                           if (TK_TMP(NM2) < TK(NM2)) then
                              TK(NM123) = TK(NM2)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        end if
                     elseif (C3D) then
                        if (ISLIP == 0) then
                           DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                  /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                           BSX1(NM123) = EVTOT(1)*real(DUDS)
                           BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                           BSX1(NM3) = EVTOT(1)*real(DUDS)
                           BSY1(NM3) = EVTOT(1)*aimag(DUDS)
                        end if
                        if (ISLIP /= 0) then
                           KSLIP = slipCoefficient(NM123, Z0B1)
                           BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                           BSX1(NM3) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM3) = KSLIP*aimag(Q(NM123, 1))
                        end if
                        !Casey 220207 from Rosemary : Update sal/temp at newly wet node
                        if ((CBAROCLINIC) .and. (IDEN == 4)) then
                           do k = 1, NFEN
                              SAL(NM3, k) = 0.5d0*(SAL(NM1, k) + SAL(NM2, k))
                              TEMP(NM3, k) = 0.5d0*(TEMP(NM1, k) + TEMP(NM2, k))
                           end do
                        end if
                     end if
                  end if
               end if
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif ((NODECODE(NM2) == 1) .and. (NODECODE(NM3) == 1)) &
               then
               if ((HTOTN2 >= HOFF) .and. (HTOTN3 >= HOFF)) then
                  NM123 = NM2
                  if (ETA2(NM3) > ETA2(NM2)) NM123 = NM3
                  DELDIST = sqrt((Y(NM1) - Y(NM123))**2.d0 &
                                 + (X(NM1) - X(NM123))**2.d0)
                  DELETA = ETA2(NM123) - ETA2(NM1)
                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                  if (DELETA < 0.d0) DELETA = 0.d0
                  if (subgrid_level0) then
                     H1 = wetDepthVertETA2(NM123)
                  else
                     H1 = ETA2(NM123) + DP(NM123)
                  end if
                  !. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                  if (C2DDI) then
                     !sb46.28sb02
                     !<<                     Convert Manning's N to Cd, if necessary.
                     if (LoadManningsN) then
                        if (subgrid_level0) then
                           if (subgridVertList(NM123) == 1) then
                              if (subgrid_level1) then
                                 FRIC(NM123) = cmfVertETA2(NM123)
                              else
                                 FRIC(NM123) = cfVertETA2(NM123)
                              end if
                           else
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                           end if
                        else
                           FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                         /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                           **(1.d0/3.d0))
                        end if
                        if (FRIC(NM123) < BFCdLLimit) then
                           FRIC(NM123) = BFCdLLimit
                        end if
                     end if

                     !-- DMW Reformulate Velocity by taking sqrt of u instead of using velmin 202207
                     if (directvelWD .and. useHF) then
                        VEL = sqrt(g*H1*DELETA/DELDIST)
                     elseif (directvelWD) then
                        VEL = sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                     else
                        TKWET = max(FRIC(NM123)*(dble(IFLINBF) + (VELMIN/H1)* &
                                                 (dble(IFNLBF) + dble(IFHYBF)* &
                                                  (1.d0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA))), &
                                    0.0001d0)
                        VEL = G*(DELETA/DELDIST)/TKWET
                     end if

                  elseif (C3D) then
                     ! solve for the depth averaged velocity,U, from the relation :
                     !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                     !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                     ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                     Z0B1 = bottomRoughnessLength(NM123)
                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                           *((H1 + Z0B1)*log((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41d0)
                  end if

                  if (VEL > VELMIN) then
                     NNODECODE(NM1) = 1
                     NMI = NM1
                     !. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                     if (C2DDI) then
                        TK_TMP(NM123) = min( &
                                        FRIC(NM123)*(dble(IFLINBF) + (TKNF*VEL/H1)* &
                                                     (dble(IFNLBF) + dble(IFHYBF)* &
                                                      (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                      (FGAMMA/FTHETA))), &
                                        2.0d0/DTDP)
                        TK(NM123) = max(TK_TMP(NM123), TK(NM123))
                        !WJP needed for internal_tide stuff
                        TKM(1:2, NM123) = TK(NM123)
                        if (NM123 == NM2) then
                           TK_TMP(NM3) = TK(NM123)
                           TK(NM3) = max(TK_TMP(NM3), TK(NM3))
                           TKM(1:2, NM3) = TK(NM3)
                           if (TK_TMP(NM3) < TK(NM3)) then
                              TK(NM123) = TK(NM3)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        else
                           TK_TMP(NM2) = TK(NM123)
                           TK(NM2) = max(TK_TMP(NM2), TK(NM2))
                           TKM(1:2, NM2) = TK(NM2)
                           if (TK_TMP(NM2) < TK(NM2)) then
                              TK(NM123) = TK(NM2)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        end if
                     elseif (C3D) then
                        if (ISLIP == 0) then
                           DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                  /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                           BSX1(NM123) = EVTOT(1)*real(DUDS)
                           BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                           BSX1(NM1) = EVTOT(1)*real(DUDS)
                           BSY1(NM1) = EVTOT(1)*aimag(DUDS)
                        end if
                        if (ISLIP /= 0) then
                           KSLIP = slipCoefficient(NM123, Z0B1)
                           BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                           BSX1(NM1) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM1) = KSLIP*aimag(Q(NM123, 1))
                        end if
                        !Casey 220207 from Rosemary : Again in the 3 wetting criteria
                        if ((CBAROCLINIC) .and. (IDEN == 4)) then
                           do k = 1, NFEN
                              SAL(NM1, k) = 0.5d0*(SAL(NM2, k) + SAL(NM3, k))
                              TEMP(NM1, k) = 0.5d0*(TEMP(NM2, k) + TEMP(NM3, k))
                           end do
                        end if
                     end if

                  end if
               end if
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif ((NODECODE(NM3) == 1) .and. (NODECODE(NM1) == 1)) &
               then
               if ((HTOTN3 >= HOFF) .and. (HTOTN1 >= HOFF)) then
                  NM123 = NM3
                  if (ETA2(NM1) > ETA2(NM3)) NM123 = NM1
                  DELDIST = sqrt((Y(NM2) - Y(NM123))**2.d0 &
                                 + (X(NM2) - X(NM123))**2.d0)
                  DELETA = ETA2(NM123) - ETA2(NM2)
                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                  if (DELETA < 0.d0) DELETA = 0.d0
                  if (subgrid_level0) then
                     H1 = wetDepthVertETA2(NM123)
                  else
                     H1 = ETA2(NM123) + DP(NM123)
                  end if
                  !. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                  if (C2DDI) then
                     !sb46.28sb02
                     !<<                     Convert Manning's N to Cd, if necessary.
                     if (LoadManningsN) then
                        if (subgrid_level0) then
                           if (subgridVertList(NM123) == 1) then
                              if (subgrid_level1) then
                                 FRIC(NM123) = cmfVertETA2(NM123)
                              else
                                 FRIC(NM123) = cfVertETA2(NM123)
                              end if
                           else
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                           end if
                        else
                           FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                         /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                           **(1.d0/3.d0))
                        end if
                        if (FRIC(NM123) < BFCdLLimit) then
                           FRIC(NM123) = BFCdLLimit
                        end if
                     end if

                     !-- DMW Reformulate Velocity by taking sqrt of u instead of using velmin 202207
                     if (directvelWD .and. useHF) then
                        VEL = sqrt(g*H1*DELETA/DELDIST)
                     elseif (directvelWD) then
                        VEL = sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                     else
                        TKWET = max(FRIC(NM123)*(dble(IFLINBF) + (VELMIN/H1)* &
                                                 (dble(IFNLBF) + dble(IFHYBF)* &
                                                  (1.d0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA))), &
                                    0.0001d0)
                        VEL = G*(DELETA/DELDIST)/TKWET
                     end if
                  elseif (C3D) then
                     ! solve for the depth averaged velocity,U, from the relation :
                     !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                     !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                     ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                     Z0B1 = bottomRoughnessLength(NM123)
                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                           *((H1 + Z0B1)*log((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41d0)
                  end if

                  if (VEL > VELMIN) then
                     NNODECODE(NM2) = 1
                     NMI = NM2
                     !. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                     if (C2DDI) then

                        TK_TMP(NM123) = min( &
                                        FRIC(NM123)*(dble(IFLINBF) + (TKNF*VEL/H1)* &
                                                     (dble(IFNLBF) + dble(IFHYBF)* &
                                                      (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                      (FGAMMA/FTHETA))), &
                                        2.0d0/DTDP)
                        TK(NM123) = max(TK_TMP(NM123), TK(NM123))
                        !WJP needed for internal_tide stuff
                        TKM(1:2, NM123) = TK(NM123)
                        if (NM123 == NM1) then
                           TK_TMP(NM3) = TK(NM123)
                           TK(NM3) = max(TK_TMP(NM3), TK(NM3))
                           TKM(1:2, NM3) = TK(NM3)
                           if (TK_TMP(NM3) < TK(NM3)) then
                              TK(NM123) = TK(NM3)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        else
                           TK_TMP(NM1) = TK(NM123)
                           TK(NM1) = max(TK_TMP(NM1), TK(NM1))
                           TKM(1:2, NM1) = TK(NM1)
                           if (TK_TMP(NM1) < TK(NM1)) then
                              TK(NM123) = TK(NM1)
                              TKM(1:2, NM123) = TK(NM123)
                           end if
                        end if
                     elseif (C3D) then
                        if (ISLIP == 0) then
                           DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                  /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                           BSX1(NM123) = EVTOT(1)*real(DUDS)
                           BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                           BSX1(NM2) = EVTOT(1)*real(DUDS)
                           BSY1(NM2) = EVTOT(1)*aimag(DUDS)
                        end if
                        if (ISLIP /= 0) then
                           KSLIP = slipCoefficient(NM123, Z0B1)
                           BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                           BSX1(NM2) = KSLIP*real(Q(NM123, 1))
                           BSY1(NM2) = KSLIP*aimag(Q(NM123, 1))
                        end if
                        !Casey 220207 from Rosemary : Again in all wetting criteria
                        if ((CBAROCLINIC) .and. (IDEN == 4)) then
                           do k = 1, NFEN
                              SAL(NM2, k) = 0.5d0*(SAL(NM1, k) + SAL(NM3, k))
                              TEMP(NM2, k) = 0.5d0*(TEMP(NM1, k) + TEMP(NM3, k))
                           end do
                        end if
                     end if

                  end if
               end if
            end if
         end if
         !WET...
         !WET...Nodal Wetting Criteria W2a
         !WET...
         NBNCTOT = NIBNODECODE(NM1) + NIBNODECODE(NM2) + NIBNODECODE(NM3)
         NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
         NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
         NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

      end do

      if (subdomainOn .and. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
      if (subdomainOn .and. enforceBN == 2) call enforceWDob() ! NCSU Subdomain

      !     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
      call UPDATEI(NNODECODE, NIBCNT, 2)
!-- DMW 202204 Update FRIC, TK, and ETA2 at ghost nodes as well
      call UPDATER(FRIC, TK, ETA2, 3)
      tmp1 = TKM(1, :)
      tmp2 = TKM(2, :)
      tmp3 = TKM(3, :)
      call UPDATER(tmp1, tmp2, tmp3, 3)
      TKM(1, :) = tmp1
      TKM(2, :) = tmp2
      TKM(3, :) = tmp3
#endif

      !wet...
      !WET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
      !*******************************************************************************************
      ! tcm v50.66.01 -- This is an additional test for wetting only when time varying
      !                  bathymetry is used and is only performed during the period of
      !                  bathymetry evolution.
      !
      if ((NDDT /= 0) .and. (dble(IT) <= BTIME_END + 1.0d0)) then
         do I = 1, NE
            NM1 = NM(I, 1)
            NM2 = NM(I, 2)
            NM3 = NM(I, 3)
            NCTOT = NODECODE(NM1) + NODECODE(NM2) + NODECODE(NM3)
            if (NCTOT < 3) then !If not wet from previous time step
               NCTOT = NNODECODE(NM1) + NNODECODE(NM2) + NNODECODE(NM3)
               if (NCTOT < 3) then !if not alreay made wet for this time step
                  ETAN1 = ETA2(NM1)
                  ETAN2 = ETA2(NM2)
                  ETAN3 = ETA2(NM3)
                  HTOTN1 = DP(NM1) + ETA2(NM1)
                  HTOTN2 = DP(NM2) + ETA2(NM2)
                  HTOTN3 = DP(NM3) + ETA2(NM3)

                  !                    if all nodes have a depth greater than or equal to
                  !                    hoff = 1.2*H0, then make the element wet
                  if ((HTOTN1 >= HOFF) .and. (HTOTN2 >= HOFF) .and. &
                      (HTOTN3 >= HOFF)) then
                     !THE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                     !  Make Node 1 Wet and set parameters
                     if (NNODECODE(NM1) /= 1) then !node 1
                        NNODECODE(NM1) = 1
                        NM123 = NM1
                        if (C2DDI) then
                           !<<                           Convert Manning's N to Cd, if necessary.
                           if (LoadManningsN) then
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              if (FRIC(NM123) < BFCdLLimit) then
                                 FRIC(NM123) = BFCdLLimit
                              end if
                           end if
                        end if
                        VEL = VELMIN
                        H1 = HTOTN1
                        if (C2DDI) then
                           TK(NM123) = FRIC(NM123)*(dble(IFLINBF) + (VEL/H1)* &
                                                    (dble(IFNLBF) + dble(IFHYBF)* &
                                                     (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                     (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2, NM123) = TK(NM123)
                        elseif (C3D) then
                           if (ISLIP == 0) then
                              DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                     /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                              BSX1(NM123) = EVTOT(1)*real(DUDS)
                              BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                              BSX1(NM3) = EVTOT(1)*real(DUDS)
                              BSY1(NM3) = EVTOT(1)*aimag(DUDS)
                           end if
                           if (ISLIP /= 0) then
                              KSLIP = slipCoefficient(NM123, bottomRoughnessLength(NM123))
                              BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                              BSX1(NM3) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM3) = KSLIP*aimag(Q(NM123, 1))
                           end if
                           !Casey 220207 from Rosemary : In ALL wetting criteria
                           if ((CBAROCLINIC) .and. (IDEN == 4)) then
                              do k = 1, NFEN
                                 SAL(NM3, k) = 0.5d0*(SAL(NM1, k) + SAL(NM2, k))
                                 TEMP(NM3, k) = 0.5d0*(TEMP(NM1, k) + TEMP(NM2, k))
                              end do
                           end if
                        end if
                     end if !end node 1

                     !  Make Node 2 Wet and set parameters
                     if (NNODECODE(NM2) /= 1) then
                        NNODECODE(NM2) = 1
                        NM123 = NM2
                        if (C2DDI) then
                           !<<                        Convert Manning's N to Cd, if necessary.
                           if (LoadManningsN) then
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              if (FRIC(NM123) < BFCdLLimit) then
                                 FRIC(NM123) = BFCdLLimit
                              end if
                           end if
                        end if
                        VEL = VELMIN
                        H1 = HTOTN2
                        if (C2DDI) then
                           TK(NM123) = FRIC(NM123)*(dble(IFLINBF) + (VEL/H1)* &
                                                    (dble(IFNLBF) + dble(IFHYBF)* &
                                                     (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                     (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2, NM123) = TK(NM123)
                        elseif (C3D) then
                           if (ISLIP == 0) then
                              DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                     /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                              BSX1(NM123) = EVTOT(1)*real(DUDS)
                              BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                              BSX1(NM1) = EVTOT(1)*real(DUDS)
                              BSY1(NM1) = EVTOT(1)*aimag(DUDS)
                           end if
                           if (ISLIP /= 0) then
                              KSLIP = slipCoefficient(NM123, bottomRoughnessLength(NM123))
                              BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                              BSX1(NM1) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM1) = KSLIP*aimag(Q(NM123, 1))
                           end if
                           !Casey 220207 from Rosemary : And again in the next two wetting criteria
                           if ((CBAROCLINIC) .and. (IDEN == 4)) then
                              do k = 1, NFEN
                                 SAL(NM1, k) = 0.5d0*(SAL(NM2, k) + SAL(NM3, k))
                                 TEMP(NM1, k) = 0.5d0*(TEMP(NM2, k) + TEMP(NM3, k))
                              end do
                           end if
                        end if
                     end if !node 2

                     !  Make Node 3 Wet and set parameters
                     if (NNODECODE(NM3) /= 1) then
                        NNODECODE(NM3) = 1
                        NM123 = NM3
                        if (C2DDI) then
                           !<<                         Convert Manning's N to Cd, if necessary.
                           if (LoadManningsN) then
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + dble(IFNLFA)*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              if (FRIC(NM123) < BFCdLLimit) then
                                 FRIC(NM123) = BFCdLLimit
                              end if
                           end if
                        end if
                        VEL = VELMIN
                        H1 = HTOTN3
                        if (C2DDI) then
                           TK(NM123) = FRIC(NM123)*(dble(IFLINBF) + (VEL/H1)* &
                                                    (dble(IFNLBF) + dble(IFHYBF)* &
                                                     (1.d0 + (HBREAK/H1)**FTHETA)** &
                                                     (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2, NM123) = TK(NM123)
                        elseif (C3D) then
                           if (ISLIP == 0) then
                              DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                     /cmplx(SIGMA(2) - SIGMA(1), 0.0d0, kind=8)
                              BSX1(NM123) = EVTOT(1)*real(DUDS)
                              BSY1(NM123) = EVTOT(1)*aimag(DUDS)
                              BSX1(NM2) = EVTOT(1)*real(DUDS)
                              BSY1(NM2) = EVTOT(1)*aimag(DUDS)
                           end if
                           if (ISLIP /= 0) then
                              KSLIP = slipCoefficient(NM123, bottomRoughnessLength(NM123))
                              BSX1(NM123) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM123) = KSLIP*aimag(Q(NM123, 1))
                              BSX1(NM2) = KSLIP*real(Q(NM123, 1))
                              BSY1(NM2) = KSLIP*aimag(Q(NM123, 1))
                           end if
                           !Casey 220207 from Rosemary : In ALL wetting criteria
                           if ((CBAROCLINIC) .and. (IDEN == 4)) then
                              do k = 1, NFEN
                                 SAL(NM2, k) = 0.5d0*(SAL(NM1, k) + SAL(NM3, k))
                                 TEMP(NM2, k) = 0.5d0*(TEMP(NM1, k) + TEMP(NM3, k))
                              end do
                           end if
                        end if
                     end if !node3

                  end if !ALL DEPTHS GREATER THAN HOFF
               end if !IF NNODECODE SUM LESS THAN 3
            end if ! IF NODECODE SUM LESS THAN 3
         end do !LOOP OVER ELEMENTS

         if (subdomainOn .and. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn .and. enforceBN == 2) call enforceWDob() ! NCSU Subdomain

         !     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
         call UPDATEI(NNODECODE, NIBCNT, 2)
#endif
         !
      end if !IT TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME
      !
      !.... END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY
      !WET..
      !WET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
      !*******************************************************************************************

      !WET...
      !WET...Nodal Wetting Criteria W2b
      !WET...Check for adjacent nodes and force nodes wet when attached
      !WET...to receiving barrier nodes
      !WET...
      do I = 1, NP
         if ((NIBCNT(I) > 0) .and. (NNODECODE(I) == 0)) then
            NNODECODE(I) = 1
         end if
      end do

      !jjwC     Use Message-Passing to update nnodecode at ghost nodes
      !jjw#ifdef CMPI
      !jjw         CALL UPDATEI(NNODECODE,IDUMY,1)
      !jjw#endif

      !WET...
      !WET...END WET/DRY SECTION - PART 2
      !WET...

      !WET...
      !WET...START WET/DRY SECTION  - PART 3
      !WET...Elemental drying criteria DE1
      !WET...This is an elemental check section designed to avoid artificial wetting of
      !WET....of control sections
      !WET...All elements where downhill flow originates from a barely wet node
      !WET....into wet nodes are forced inactive; the only exception is receiving
      !WET....overtopped barrier nodes
      !WET...
      do I = 1, NE
         NM1 = NM(I, 1)
         NM2 = NM(I, 2)
         NM3 = NM(I, 3)
         NBNCTOT = NIBCNT(NM1)*NIBCNT(NM2)*NIBCNT(NM3)
         if (NBNCTOT == 0) then !No barrier/pipe receiving nodes in this elem
            ETAN1 = ETA2(NM1)
            ETAN2 = ETA2(NM2)
            ETAN3 = ETA2(NM3)
            if (subgrid_level0) then
               HTOTN1 = gridDepthVertETA2(NM1)
               HTOTN2 = gridDepthVertETA2(NM2)
               HTOTN3 = gridDepthVertETA2(NM3)
            else
               HTOTN1 = DP(NM1) + ETA2(NM1)
               HTOTN2 = DP(NM2) + ETA2(NM2)
               HTOTN3 = DP(NM3) + ETA2(NM3)
            end if
#ifdef SB_WETDRY
!...Find the heighest point on the bed in the element.  sb v46.28.sb05.06 11/01/2006
            if (DP(NM1) <= DP(NM2) .and. DP(NM1) <= DP(NM3)) then
               DPMIN = DP(NM1)
            else if (DP(NM2) <= DP(NM3) .and. DP(NM2) <= DP(NM1)) then
               DPMIN = DP(NM2)
            else if (DP(NM3) <= DP(NM1) .and. DP(NM3) <= DP(NM2)) then
               DPMIN = DP(NM3)
            end if
#endif
#ifndef SB_WETDRY
            ! jgf52.08.08: Analyst can eliminate noff from
            ! consideration in fort.15 file.
            if (noffActive .eqv. .true.) then
               !...ABC pattern
               if ((ETAN1 >= ETAN2) .and. (ETAN2 > ETAN3)) then
                  if ((HTOTN1 < HOFF) .or. (HTOTN2 < HOFF)) NOFF(I) = 0
               end if
               if ((ETAN2 >= ETAN3) .and. (ETAN3 > ETAN1)) then
                  if ((HTOTN2 < HOFF) .or. (HTOTN3 < HOFF)) NOFF(I) = 0
               end if
               if ((ETAN3 >= ETAN1) .and. (ETAN1 > ETAN2)) then
                  if ((HTOTN3 < HOFF) .or. (HTOTN1 < HOFF)) NOFF(I) = 0
               end if
               !...ACB pattern
               if ((ETAN1 >= ETAN3) .and. (ETAN3 > ETAN2)) then
                  if ((HTOTN1 < HOFF) .or. (HTOTN3 < HOFF)) NOFF(I) = 0
               end if
               if ((ETAN2 >= ETAN1) .and. (ETAN1 > ETAN3)) then
                  if ((HTOTN2 < HOFF) .or. (HTOTN1 < HOFF)) NOFF(I) = 0
               end if
               if ((ETAN3 >= ETAN2) .and. (ETAN2 > ETAN1)) then
                  if ((HTOTN3 < HOFF) .or. (HTOTN2 < HOFF)) NOFF(I) = 0
               end if
            end if
#else
! jgf52.08.08: Analyst can eliminate noff from
! consideration in fort.15 file.
            if (noffActive .eqv. .true.) then
!...ABC pattern
               if ((ETAN1 >= ETAN2) .and. (ETAN2 > ETAN3)) then
                  if ((HTOTN1 < HOFF)) NOFF(I) = 0
                  if ((ETAN1 - ETAN2) < (ETAN2 - ETAN3)) then
                     if (HTOTN2 < HOFF) NOFF(I) = 0
                  end if
               end if
               if ((ETAN2 >= ETAN3) .and. (ETAN3 > ETAN1)) then
                  if ((HTOTN2 < HOFF)) NOFF(I) = 0
                  if ((ETAN2 - ETAN3) < (ETAN3 - ETAN1)) then
                     if (HTOTN3 < HOFF) NOFF(I) = 0
                  end if
               end if
               if ((ETAN3 >= ETAN1) .and. (ETAN1 > ETAN2)) then
                  if ((HTOTN3 < HOFF)) NOFF(I) = 0
                  if ((ETAN3 - ETAN1) < (ETAN1 - ETAN2)) then
                     if (HTOTN1 < HOFF) NOFF(I) = 0
                  end if
               end if
!...ACB pattern
               if ((ETAN1 >= ETAN3) .and. (ETAN3 > ETAN2)) then
                  if ((HTOTN1 < HOFF)) NOFF(I) = 0
                  if ((ETAN1 - ETAN3) < (ETAN3 - ETAN2)) then
                     if (HTOTN3 < HOFF) NOFF(I) = 0
                  end if
               end if
               if ((ETAN2 >= ETAN1) .and. (ETAN1 > ETAN3)) then
                  if ((HTOTN2 < HOFF)) NOFF(I) = 0
                  if ((ETAN2 - ETAN1) < (ETAN1 - ETAN3)) then
                     if (HTOTN1 < HOFF) NOFF(I) = 0
                  end if
               end if
               if ((ETAN3 >= ETAN2) .and. (ETAN2 > ETAN1)) then
                  if ((HTOTN3 < HOFF)) NOFF(I) = 0
                  if ((ETAN3 - ETAN2) < (ETAN2 - ETAN1)) then
                     if (HTOTN2 < HOFF) NOFF(I) = 0
                  end if
               end if
            end if !noffActive.eqv..true.
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
            if (noffActive .eqv. .true.) then
               if (NOFFOLD(I) == 0) then
                  if (HTOTN1 >= HTOTN2 .and. HTOTN1 >= HTOTN3) then
                     if (ETAN1 < (-DPMIN + H0)) NOFF(I) = 0
                  end if
                  if (HTOTN2 >= HTOTN3 .and. HTOTN2 >= HTOTN1) then
                     if (ETAN2 < (-DPMIN + H0)) NOFF(I) = 0
                  end if
                  if (HTOTN3 >= HTOTN1 .and. HTOTN3 >= HTOTN2) then
                     if (ETAN3 < (-DPMIN + H0)) NOFF(I) = 0
                  end if
               end if
            end if
#endif
         end if
      end do

#ifdef SB_WETDRY
!WET......added by sb on 11/02/2006
!WET...
!WET...This section is added after we realize that it's not possible
!WET...to compute a correct flow going through two elements if
!WET...the elements are connected just by one node. i.e., elements need
!WET...to share an edge to let the flow go through between the elements.
!WET...Therefore, in this section, a node is determined to be dry
!WET...if two elements are connected at one node, not sharing an edge.
!WET...As it seemed this procedure needed NOFF information,
!WET...although this section changes NNODECODE,
!WET...I put this section here, rather than the end of PART 2.
!WET...I tried setting NOFF(I) = 0, but it didn't shut down the flow.
!WET...
      do I = 1, NE
         NM1 = NM(I, 1)
         NM2 = NM(I, 2)
         NM3 = NM(I, 3)
         if (NOFF(I) == 1 .and. &
             NNODECODE(NM1) == 1 .and. &
             NNODECODE(NM2) == 1 .and. &
             NNODECODE(NM3) == 1) then
            do K = 1, 3
               NM1 = NM(I, K)
               NM2 = NM(I, mod(K + 0, 3) + 1)
               NM3 = NM(I, mod(K + 1, 3) + 1)

               NWETNEI = 0
               NWETADJ = 0
               do J = 1, MNEI
                  N = NeiTabEle(NM1, J)
                  if (N == 0) cycle
                  if (N == I) cycle

                  NMN1 = NM(N, 1)
                  NMN2 = NM(N, 2)
                  NMN3 = NM(N, 3)
                  if (NOFF(N) == 1 .and. &
                      NNODECODE(NMN1) == 1 .and. &
                      NNODECODE(NMN2) == 1 .and. &
                      NNODECODE(NMN3) == 1) then
                     NWETNEI = NWETNEI + 1

                     if ((NMN1 == NM2 .or. NMN1 == NM3) .or. &
                         (NMN2 == NM2 .or. NMN2 == NM3) .or. &
                         (NMN3 == NM2 .or. NMN3 == NM3)) then
                        NWETADJ = NWETADJ + 1
                     end if
                  end if
               end do

               if (NWETNEI > 0 .and. NWETADJ == 0 .and. &
                   NIBCNT(NM1) == 0) then
                  NNODECODE(NM1) = 0
               end if
            end do
         end if
      end do
#endif

      !WET...
      !WET...END WET/DRY SECTION  - PART 3
      !WET...

      !.....
      !.....Nodal Drying Criteria for the ibtype=64 weir boundary
      !.....Check the node if it is along the ibtype=64 weir boundary and if the weir is
      !.....submerged. If it is the case and the node does not belong to at least
      !.....one wet element, set the nnodecode back to 0.   SB
      !.....
      if (NFLUXIB64_GBL > 0) then
         do I = 1, NP
            if (LBArray_Pointer(I) > 0) then
               J = LBArray_Pointer(I)
               if ((ISSUBMERGED64(J) > 0) .and. (NNODECODE(I) == 1)) then
                  NNBB2 = IBCONN(J)
                  NWetEle = 0
                  do K = 1, NNeighEle(I)
                     IE = NeiTabEle(I, K)
                     NCTOT = NNODECODE(NM(IE, 1)) + &
                             NNODECODE(NM(IE, 2)) + NNODECODE(NM(IE, 3))
                     if (NCTOT == 3) then
                        NWetEle = NWetEle + 1
                     end if
                  end do
                  if (NWetEle == 0) then
                     NNODECODE(I) = 0
                  end if
               end if
            end if
         end do
      end if

      !.....
      !..... Section for condensed nodes
      !..... Setting the same nodecode value at grouped condensed nodes   SB
      !.....
      if (LoadCondensedNodes) then
         do K = 1, NListCondensedNodes
            I = ListCondensedNodes(K, 1)
            if (I == 0) cycle
            ! 0) DEFAULT VALUE
            NNC = NNODECODE(I)
            ! 1) FIND UPDATED VALUE
            do L = 2, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               if (NNODECODE(J) /= NODECODE(J)) then
                  NNC = NNODECODE(J)
                  exit ! EXiT THE LOOP IF AN UPDATED NODECODE IS FOUND
               end if
            end do
            ! 2) DISTRIBUTE IT
            NNODECODE(I) = NNC
            do L = 2, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               NNODECODE(J) = NNC
            end do
         end do
      end if

      !WET...
      !WET...START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2
      !WET...Update number of active elements (MJU) and the total area (TotalArea) connected
      !WET...to a node. If these are zero, the node is landlocked and should be dried.
      !WET...These depend on NNODECODE which varies during the time step
      !WET...
      do I = 1, NP
         MJU(I) = 0
         TotalArea(I) = 0.d0
      end do
      do IE = 1, NE
         NM1 = NM(IE, 1)
         NM2 = NM(IE, 2)
         NM3 = NM(IE, 3)
         NC1 = NNODECODE(NM1)
         NC2 = NNODECODE(NM2)
         NC3 = NNODECODE(NM3)

         NCEle = NC1*NC2*NC3*NOFF(IE)
         AreaEle = dble(NCEle)*Areas(IE)/2.d0
         MJU(NM1) = MJU(NM1) + NCEle
         MJU(NM2) = MJU(NM2) + NCEle
         MJU(NM3) = MJU(NM3) + NCEle
         TotalArea(NM1) = TotalArea(NM1) + AreaEle
         TotalArea(NM2) = TotalArea(NM2) + AreaEle
         TotalArea(NM3) = TotalArea(NM3) + AreaEle
      end do

      !jjwnote - looks like this is used later in momentum equations
      !jjwnote - this has implications on making this into a subroutine

      do I = 1, NP
         if ((NNODECODE(I) == 1) .and. (MJU(I) == 0)) then
            NNODECODE(I) = 0
         end if
         if (MJU(I) == 0) MJU(I) = 1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
      end do

      !     WET...
      !     WET...END WET/DRY SECTION - PART 4
      !     WET...

      !jjwnote - may have to pass TotalArea and mju as well

      if (subdomainOn .and. enforceBN == 1) call enforceWDcb() ! NCSU Subdomain
      if (subdomainOn .and. enforceBN == 2) call enforceWDob() ! NCSU Subdomain

      !     Use Message-Passing to update nnodecode at ghost nodes
#ifdef CMPI
      call UPDATEI(NNODECODE, IDUMY, 1)
      call UPDATER(TotalArea, DUMY1, DUMY2, 1)
#endif
      !WET...
      !WET...WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE
      !WET...Check to see if any wetting occurred & update NODECODE
      !WET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
      !WET...
      do I = 1, NP
         if (NNODECODE(I) /= NODECODE(I)) then
            NODECODE(I) = NNODECODE(I)
            NCCHANGE = NCCHANGE + 1
         end if
      end do
      !WET...
      !WET...END WET/DRY SECTION - PART 5
      !WET...

      !WET...
      !WET...WET/DRY SECTION - PART 6
      !WET...Check to see if any NOFF changed requiring the matrix to be reset
      !WET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
      !WET...
      do I = 1, NE
         if (NOFF(I) /= NOFFOLD(I)) NCCHANGE = NCCHANGE + 1
      end do
      !WET...
      !WET... jgf45.06 If there has been any wetting or drying in any
      !WET... of the subdomains, the NCCHANGE flag will be activated on all
      !WET... of the subdomains, to prevent them from getting out of sync
      !WET... with their MPI calls as some reset the GWCE and others do not.
      !WET...

      if (LoadSubgridBarrier) then
         do I = 1, NE
            if (NOFF_SG(I) == 0) then
               if (NOFF(I) /= 0) then
                  NOFF(I) = 0
                  NCCHANGE = NCCHANGE + 1
               end if
            end if
         end do
      end if

      !.....   DMW 202207 Compute a slope limiting factor if needed
      call computeAlpha()

#ifdef CMPI
!jgf48.4619 implementing Seizo's changes for Lumped, fully
! explicit operation. In that case, the GWCE LHS matrix is
! recalculated on each individual subdomain that has wetted
! or dried, without recourse to MPI, eliminating the need
! for the call to the subroutine WetDrySum.
      if (ILump == 0) then
         call WetDrySum(NCCHANGE)
      else
         NCCHANGE = NCCHANGE ! jgf48.4619 do nothing
      end if
#endif
      !WET...
      !WET...END WET/DRY SECTION - PART 6
      !WET...

      !.....
      !-----------------------------------------------------------------------
   end subroutine computeWettingAndDrying
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !   Apply a subgrid barrier.
   !
   !   Assumptions: When the subgrid barrier is not being overtopped,
   !                the cross-barrier flow is disabled and the nodes
   !                are placed at the elevation found in the mesh
   !                geometry. When the barrier becomes overtopped,
   !                the nodes are elevated temporarily during the
   !                overtopping process. When the overtopping
   !                concludes, the nodes are returned to their mesh
   !                specified elevations and the cross-barrier flow
   !                is turned off again
   !-----------------------------------------------------------------------
   subroutine applySubgridBarrier()
      !-----------------------------------------------------------------------
      use MESH, only: DP, DP_ORIG, NM, NE
      use GLOBAL, only: ETA2, NOFF, H0, NODECODE
      use NODALATTRIBUTES, only: SubgridBarrier, &
                                 SubgridBarrierOvertopping
      implicit none
      integer :: I, NM1, NM2, NM3, NCELE
      real(8) :: SB1, SB2, SB3
      real(8) :: DOT, DOT1, DOT2, DOT3
      real(8) :: E1, E2, E3, ETA
      real(8), parameter :: SBDT = 200d0

      NOFF_SG(:) = 1

      do I = 1, NE
         NM1 = NM(I, 1)
         NM2 = NM(I, 2)
         NM3 = NM(I, 3)
         SB1 = subgridBarrier(NM1, 1)
         SB2 = subgridBarrier(NM2, 1)
         SB3 = subgridBarrier(NM3, 1)

         NCELE = NODECODE(NM1)*NODECODE(NM2)*NODECODE(NM3)*NOFF(I)
         if (NCELE == 0) cycle

         if (SB1 < DP_ORIG(NM1) .and. &
             SB2 < DP_ORIG(NM2) .and. &
             SB3 < DP_ORIG(NM3)) then

            E1 = ETA2(NM1)
            E2 = ETA2(NM2)
            E3 = ETA2(NM3)
            ETA = max(E1, E2, E3)
            DOT1 = E1*dble(NODECODE(NM1)) + SB1
            DOT2 = E2*dble(NODECODE(NM2)) + SB2
            DOT3 = E3*dble(NODECODE(NM3)) + SB3
            DOT = max(DOT1, DOT2, DOT3)

            if (DOT > -1.2d0*H0) then
               DP(NM1) = max(SB1, DP(NM1) + SB1/SBDT)
               DP(NM2) = max(SB2, DP(NM2) + SB2/SBDT)
               DP(NM3) = max(SB3, DP(NM3) + SB3/SBDT)
               SubgridBarrierOvertopping(NM1) = .true.
               SubgridBarrierOvertopping(NM2) = .true.
               SubgridBarrierOvertopping(NM3) = .true.
            else
               DP(NM1) = min(DP_ORIG(NM1), DP(NM1) - SB1/SBDT)
               DP(NM2) = min(DP_ORIG(NM2), DP(NM2) - SB2/SBDT)
               DP(NM3) = min(DP_ORIG(NM3), DP(NM3) - SB3/SBDT)
               NOFF_SG(I) = 0
               SubgridBarrierOvertopping(NM1) = .false.
               SubgridBarrierOvertopping(NM2) = .false.
               SubgridBarrierOvertopping(NM3) = .false.
            end if
         else
            SubgridBarrierOvertopping(NM1) = .false.
            SubgridBarrierOvertopping(NM2) = .false.
            SubgridBarrierOvertopping(NM3) = .false.
         end if
      end do

      !-----------------------------------------------------------------------
   end subroutine applySubgridBarrier
   !-----------------------------------------------------------------------

   !..... DMW 202207 Compute a slope limiting factor for scaling down the
   !..... elevation gradient if the slope of the water surface aligned with the
   !..... bathy slope becomes too large, i.e., when it exceeds slim.
   !-----------------------------------------------------------------------
   subroutine computeAlpha()
      !-----------------------------------------------------------------------
      use MESH, only: DP, NM, NE, SFacEle, SFMYEle, SFMXEle, FDXE, FDYE, &
                      areas
      use global, only: eta2, noff, nodecode, alphal, IFSFM, slim
      use subgrid, only: subgrid_level0, gridDepthVertETA2
      implicit none
      integer :: I, NM1, NM2, NM3, NC1, NC2, NC3, NCELE
      real(8) :: sfacavg, SFmxAvg, SFmyAvg, DPN1, DPN2, DPN3
      real(8) :: DETAM, DPM, sfdxfac, sfdyfac, hn1, hn2, hn3
      real(8) :: detadxa, detadya, dpdxa, dpdya
      real(8) :: fdx1, fdx2, fdx3, fdy1, fdy2, fdy3
      real(8) :: etaN1, etaN2, etaN3
      real(8) :: DDTMP !tcm 20250503

      ALPHAL(:) = 1.d0

      do I = 1, NE
         NM1 = NM(I, 1)
         NM2 = NM(I, 2)
         NM3 = NM(I, 3)
         NC1 = NODECODE(NM1)
         NC2 = NODECODE(NM2)
         NC3 = NODECODE(NM3)
         NCEle = NC1*NC2*NC3*NOFF(I)
         if (NCEle == 1) then

            if (subgrid_level0) then
               HN1 = gridDepthVertETA2(NM1)
               HN2 = gridDepthVertETA2(NM2)
               HN3 = gridDepthVertETA2(NM3)
            else
               HN1 = ETA2(NM1) + DP(NM1)
               HN2 = ETA2(NM2) + DP(NM2)
               HN3 = ETA2(NM3) + DP(NM3)
            end if
            if ((HN1 <= 5.d-1) .or. (HN2 <= 5.d-1) .or. (HN3 <= 5.d-1)) then

               SFacAvg = SFacEle(I)
               SFmxAvg = SFMXEle(I)
               SFmyAvg = SFMYEle(I)
               sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg; 
               sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg; 
               FDX1 = FDXE(1, I)*sfdxfac !b1
               FDX2 = FDXE(2, I)*sfdxfac !b2
               FDX3 = FDXE(3, I)*sfdxfac !b3
               FDY1 = FDYE(1, I)*sfdyfac !a1
               FDY2 = FDYE(2, I)*sfdyfac !a2
               FDY3 = FDYE(3, I)*sfdyfac !a3

               ETAN1 = ETA2(NM1)
               ETAN2 = ETA2(NM2)
               ETAN3 = ETA2(NM3)

               DPN1 = DP(NM1)
               DPN2 = DP(NM2)
               DPN3 = DP(NM3)

               DPDXA = DPN1*FDX1 + DPN2*FDX2 + DPN3*FDX3
               DPDYA = DPN1*FDY1 + DPN2*FDY2 + DPN3*FDY3

               DETADXA = ETAN1*FDX1 + ETAN2*FDX2 + ETAN3*FDX3
               DETADYA = ETAN1*FDY1 + ETAN2*FDY2 + ETAN3*FDY3
               DPM = sqrt(DPDXA*DPDXA + DPDYA*DPDYA)/AREAS(I)

               if (DPM > 1d-5) then
                  DDTMP = DETADXA*(-DPDXA) + DETADYA*(-DPDYA) !tcm 20250503
                  if (abs(DDTMP) > 1d-12) then !tcm 20250503 added test and else condition
                     ALPHAL(I) = SLIM*(sqrt(DPDXA*DPDXA + DPDYA*DPDYA)/ &
                                       ((1.d0/AREAS(I))*(DDTMP))) !tcm changed to DDTMP
                  else
                     ALPHAL(I) = 1.d0
                  end if
               else
                  DETAM = sqrt(DETADXA*DETADXA + DETADYA*DETADYA)/AREAS(I)
                  if (abs(DETAM) > 1d-12) then !tcm 20250503 added test and else condition
                     ALPHAL(I) = SLIM/DETAM
                  else
                     ALPHAL(I) = 1.d0
                  end if
               end if

               if (ALPHAL(I) <= 0.d0) ALPHAL(I) = 1.d0
               if (ALPHAL(I) >= 1.d0) ALPHAL(I) = 1.d0
            end if

         end if
      end do

      !-----------------------------------------------------------------------
   end subroutine computeAlpha
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   real(8) function computeTotAreaAtNode(INODE)
      !-----------------------------------------------------------------------
      use MESH, only: NNeighEle, NeiTabEle, Areas
      implicit none
      integer, intent(in) :: INODE
      integer :: I, IE
      real(8) :: TotArea

      TotArea = 0.d0
      do I = 1, NNeighEle(INODE)
         IE = NeiTabEle(INODE, I)
         TotArea = TotArea + Areas(IE)
      end do
      computeTotAreaAtNode = 0.5d0*TotArea
      return
      !-----------------------------------------------------------------------
   end function computeTotAreaAtNode

   !-----------------------------------------------------------------------
   !> Bottom roughness length used for 3D wetting: the nodal z0b attribute
   !> if loaded, otherwise derived from Manning's n, otherwise the global z0b.
   !-----------------------------------------------------------------------
   real(8) function bottomRoughnessLength(node) result(z0)
      use ADC_CONSTANTS, only: g
      use global, only: eta2, ifnlfa
      use global_3dvs, only: z0b
      use mesh, only: dp
      use nodalattributes, only: loadZ0B_var, z0b_var, loadManningsN, manningsn
      integer, intent(in) :: node

      if (loadZ0B_var) then
         z0 = z0b_var(node)
      elseif (loadManningsN) then
         z0 = (dp(node) + dble(ifnlfa)*eta2(node))*exp(-(1.0d0 + &
                                                         ((0.41d0*(dp(node) + dble(ifnlfa)*eta2(node))**(1.0d0/6.0d0))/ &
                                                          (manningsn(node)*sqrt(g)))))
      else
         z0 = z0b
      end if
   end function bottomRoughnessLength

   !-----------------------------------------------------------------------
   !> Slip coefficient for the 3D bottom stress at a newly wetted node.
   !>   ISLIP=1: KP is the linear slip coefficient
   !>   ISLIP=2: log-law quadratic coefficient, bounded to [0.0025, 1]*|U|
   !>   ISLIP=3: KP is the quadratic slip coefficient, KP*|U| bounded to [0, 0.1]
   !-----------------------------------------------------------------------
   real(8) function slipCoefficient(node, z0) result(kslip)
      use global, only: eta2, ifnlfa
      use global_3dvs, only: a, b, islip, kp, sigma, q
      use mesh, only: dp
      integer, intent(in) :: node
      real(8), intent(in) :: z0 !< bottom roughness length at node

      select case (islip)
      case (1)
         kslip = kp
      case (2)
         kslip = (1.d0/ &
                  ((1.d0/0.41d0)* &
                   log((abs(((sigma(2) - sigma(1))/(a - b))* &
                            (dp(node) + dble(ifnlfa)*eta2(node))) &
                        + z0) &
                       /(z0))))**2.d0 &
                 *abs(q(node, 1))
         if (kslip > 1.d0*abs(q(node, 1))) kslip = 1.d0*abs(q(node, 1))
         if (kslip < 0.0025d0*abs(q(node, 1))) kslip = 0.0025d0*abs(q(node, 1))
      case (3)
         kslip = kp*abs(q(node, 1))
         if (kslip > 0.1d0) kslip = 0.1d0
         if (kslip < 0.0d0) kslip = 0.0d0
      case default
         kslip = 0.d0
      end select
   end function slipCoefficient
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
end module mod_wetdry
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
