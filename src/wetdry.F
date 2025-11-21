!------------------------------------------------------------------------
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the 
! License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------
!                     M O D U L E   W E T   D R Y      
!------------------------------------------------------------------------
!> @file wetdry.F
!>
!> @author Various
!>
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module executes wetting and drying algorithm.
!>
!> @details
!>
!> Logical Variable List (default value .FALSE., set in global.f)
!>    C2DDI            - 2D Depth Integrated model run
!>    C3D              - 3D model run
!>    C3DDSS           - Stress form of the 3D momentum equations
!>    C3DVS            - Velocity form of the 3D momentum equations
!>
!> -- Important Variable Notes --
!>          nnodecode is a working variable that can change more than once
!>                    during a time step
!>          nnodecode = 0 for a dry node
!>          nnodecode = 1 for a wet node
!>          nodecode  is a more static version of nnodecode that is 
!>                      reconciled once and for all at the end time step
!>     
!> -- Algorithm Description --
!>
!>             (   DRYING CRITERIA   )
!>     
!>     A node should be dry under two conditions:
!>     D1.) If the total water depth falls below H0.
!>          Note: if the total water depth falls below 0.8*H0, the surface
!>          elevation is lifted up so that the total water depth = 0.8*H0
!>        
!>     D2.) If the node is connected to only nonfunctioning (dry) elements
!>          In this case the node is dried due to becoming landlocked.
!>          Note: this criteria is applied after all other wetting and
!>          drying criteria
!>     
!>     An element should be dry under the following conditions:
!>     DE3.) This is an elemental check section designed to avoid
!>           artificial wetting of control sections
!>           All elements where downhill flow originates from a barely
!>           wet node (defined as 1.2*H0) into wet nodes are forced
!>           inactive; the only exception is receiving overtopped
!>           barrier nodes
!>     
!>             (   WETTING CRITERIA   )
!>     
!>     A node should be wet under two conditions:
!>     W1.) If 2 nodes on an element are wet and one is dry, wet the dry
!>          node if the water level at one of the wet nodes is greater
!>          than the water level at the dry node and the steady-state
!>          velocity that would result from a balance between the water
!>          level gradient and bottom friction would yield a velocity
!>          > velmin.
!>          Note that the criteria outlined in DE3 must also be satified
!>          before the node is allowed to wet
!>     
!>     W2.) If an element has a node lying on a receiving internal barrier
!>          boundary or specified discharge boundary that is actively
!>          discharging flow into the domain at that node, all nodes in
!>          this element must stay wet.
!>     
!>             (  VELOCITY BOUNDARY CONDITION  )
!>     
!>     Either a natural or essential boundary condition can be used as a
!>     velocity boundary condition in the momentum equation solution
!>     along a wet/dry boundary. 
!>     To use a natural boundary condition, do nothing along the wet/dry
!>     interface.
!>     To use an essential, no velocity boundary condition, identify the
!>     nodes along the wet/dry interface and zero out the velocity at
!>     the nodes.  Interface nodes can easily be identified by comparing
!>     the number of active elements a node is connected to (MJU) to the
!>     total number of elements a node is connected to (NODELE). 
!>     If MJU < NODELE for any node, it must lie along the wet/dry
!>     interface.  See further comments at the end of the momentum
!>     equation solution section.
!>
!> @todo Need to add documentation above for new features:
!>       weir boundary (ibtype=64) and add (un)setMessageSource sections
!>       to all subroutines. Also, some subroutines within wetdry as
!>       well as other modules are setting their own habsmin and velmin
!>       values or hard-coding values instead of using the values here.
!>       Do a search and fix all these.
!------------------------------------------------------------------------

      module wetdry
      
      real(8) :: habsmin !local versions created in multiple modules??
      real(8) :: hoff
      real(8) :: velmin
      
      integer,allocatable ::    nibcnt(:)
      integer,allocatable ::    noffold(:)
      integer,allocatable ::    noff_sg(:)

      ! jgf52.08.08: Enable analyst to remove noff from wetting and
      ! drying algorithm; .true. by default to conform with prior
      ! ADCIRC versions. 
      logical :: noffActive = .true.

      ! cms 11/20/25: added public/private declaration
      ! Variables
      private :: habsmin, hoff, nibcnt, noffold, noff_sg
      public :: noffActive, velmin !habsmin could be reused elsewhere

      !Subroutines and functions
      private :: computeAlpha, computeTotAreaAtNode, applySubgridBarrier
      public :: initializeWettingAndDrying, computeWettingAndDrying

      contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                             P U B L I C
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief Allocates arrays and initializes wetdry states/constants.
!>
!> @details cms 11/20/25 updated to move all initializations here -mvd
!>          noff=1 from setColdWetDryStateAndWaterSurfaceElevation()
!>          in cstart to here so that all initialized at the same time.
!>          Also changed to explicit array bounds.
!----------------------------------------------------------------------
      subroutine initializeWettingAndDrying()
!----------------------------------------------------------------------
      use sizes, only : mne, mnp
      use global, only : h0, nodecode, nnodecode, noff
      use nodalattributes,only: loadSubgridBarrier
      implicit none

      allocate ( nibcnt(mnp) )
      allocate ( noffold(mne))
      if(loadSubgridBarrier)allocate(noff_sg(mne))

      ! Initialize all wetdry states (noff mvd from cstart)
      nnodecode(1:mnp) = 1
      nodecode(1:mnp) = 1
      noff(1:mne) = 1
      noffold(1:mne) = 1

      ! Compute wetdry threshold constants
      habsmin=0.8d0*h0
      hoff=1.2d0*h0

      return
!----------------------------------------------------------------------
      end subroutine initializeWettingAndDrying
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> @brief Determines which nodes should be wet or dry.
!>
!> @details cms 11/20/25 updated to use explicit array structure for
!>          all array initializations and replaced code in WET/DRY Part4
!>          with a call to updated totalAreaCalc() in mesh module.
!----------------------------------------------------------------------
      subroutine computeWettingAndDrying(it)
!----------------------------------------------------------------------
      use sizes, only : mne, mnp, mnei
      use nodalattributes,only: loadSubgridBarrier
      use ADC_CONSTANTS, only: g
      use global, only : noff, nodecode, nnodecode, eta1, eta2, tk,
     &    nolifa, bsx1, bsy1, btime_end, C2DDI, C3D, h0, ifnlfa, nddt,
     &    dtdp, nibnodecode, ilump, ncchange, tkm, directvelWD, useHF,
     &    cbaroclinic, iden 
#ifdef CMPI
     &    ,idumy,dumy1,dumy2
      use messenger
#endif
      use mesh, only : ne, np, nm, dp, totalArea, nm, x, y,
     &   LBArray_Pointer, nneighele, neitabele, onIBTYP64Floodplain,
     &   neitab, nneigh, totalAreaCalc
      use nodalattributes, only : BFCdLLimit, fgamma, ftheta, fric,
     &   manningsn, hbreak, ifhybf, ifnlbf, iflinbf, loadManningsN,
     &   loadZ0B_var, z0b_var, loadSubgridBarrier,
     &   loadCondensedNodes, NCondensedNodes,
     &   ListCondensedNodes, NListCondensedNodes,
     &   NNodesListCondensedNodes
      use boundaries, only : IBCONN, isSubmerged64, NFLUXIB64_GBL,
     &    barinht, lbcodeI
      use global_3dvs, only : a, b, islip, kp, z0b, sigma, evtot, q
!Casey 220207 from Rosemary
     &    , nfen, sal, temp
      use subdomain, only : subdomainOn, enforceBN, enforceWDcb, 
     &    enforceWDob
      use subgrid, only: level0, subgridVertList, gridDepthVertETA2,
     &   level1, cmfVertETA2, cfVertETA2, wetDepthVertETA2,
     &   wetFracVertTab,gridDepthVertTab, numPhi

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
      real(8) :: etaN1,etaN2,etaN3
      real(8) :: hTotN1,hTotN2,hTotN3      
      real(8) :: deldist,deleta
      real(8) :: htot
      real(8) :: h1
      real(8) :: tk_tmp(np),tmp1(np),tmp2(np),tmp3(np)
      real(8) :: totArea1,totArea2,eta2Avg
      integer :: tk_cnt(np)
      integer :: nbnctot
      integer :: i,j,k,kk,l
      integer :: ie          
      real(8) :: tknf !dmw 202207, added for "new" formula for vel 
      integer :: nnc
      integer :: nnbb2
      integer :: nwetele
      integer :: ibtype64floodplaincnt = 0
      integer :: numGreater
      integer :: cntwet, cntdry, fstwet, curwet
  
      if (nolifa.ne.2) then
         return ! wetting and drying is not active 
      endif

         ! WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP
         TK_TMP(:) = 0.D0
         TK_CNT(:) = 0
         if (directvelWD.AND.useHF) then
            TKNF = 100.D0
         else
            TKNF = 1.D0
         endif
 
         !cms: replaced loops with explicit array bounds
         nibcnt(1:np) = 0
         noffold(1:ne) = noff(1:ne)
         noff((1:ne) = 1
         !do I=1,NP
         !   NIBCNT(I) = 0
         !end do
         !do I=1,NE
         !   noffold(I)=noff(I)
         !   noff(I)=1
         !end do

         if(loadSubgridBarrier)call applySubgridBarrier()


         !  WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
         !Drying Criteria D1: depends on nodecode and updates nodecode
         do I=1,NP
            if(nodecode(I).eq.1) then
               if(level0)then
                  HTOT = gridDepthVertETA2(i)
               else
                  HTOT=DP(I)+ETA2(I)
               endif
               if(HTOT.LE.H0) then
                  if(level0)then
                     !JLW: need to interpolate here to set eta2
                     if(HTOT.lt.habsmin)then
                        gridDepthVertETA2(i) = habsmin
                        wetDepthVertETA2(i) = habsmin
                        if(subgridVertList(i).eq.1)then
                           numGreater = 0
                           do j=1,numPhi
                              if(gridDepthVertETA2(i).lt.
     &                           gridDepthVertTab(i,j))then
                                 numGreater = numGreater + 1
                              endif
                           end do
                           if(numGreater.eq.numPhi)then
                              eta2(i) = wetFracVertTab(i,1)
                           elseif(numGreater.eq.0)then
                              eta2(i) = wetFracVertTab(i,numPhi)
     &                           + (habsmin-gridDepthVertTab(i,numPhi))
                           else
               !JLW: now back calculate eta2 based on this change
                              eta2(i) = ((gridDepthVertETA2(i)
     &                        - gridDepthVertTab(i,numPhi-numGreater))
     &                        / (gridDepthVertTab(i,numPhi-numGreater+1)
     &                        - gridDepthVertTab(i,numPhi-numGreater))
     &                        * (wetFracVertTab(i,numPhi-numGreater+1)
     &                        - wetFracVertTab(i,numPhi-numGreater))
     &                        + wetFracVertTab(i,numPhi-numGreater))
                           endif
                        else
                           eta2(i) = habsmin - DP(I)
                        endif
                     endif
                     nnodecode(I) = 0
                     nodecode(I) = 0
                     ncchange = ncchange + 1
                  else
                     if(HTOT.lt.habsmin) ETA2(I)=habsmin-DP(I)
                        nnodecode(I)=0
                        nodecode(I)=0
                        ncchange=ncchange+1 ! set to 0 in gwce
                  endif
               endif
            endif
         end do
         ! END WET/DRY SECTION - PART 1


         ! Nodal Wetting Criteria for the ibtype=64 weir boundary
         ! Check if node is along the ibtype=64 weir boundary and
         ! if the weir is submerged. If true, set the nnodecode to 1. SB
         if (NFLUXIB64_GBL.gt.0) then
            do I=1,NP
               J = LBArray_Pointer(I)
               if(J.gt.0) then
                  if((LBCODEI(J)).eq.64.AND.(nodecode(I).eq.0)) then
                     if((ISSUBMERGED64(J).gt.0)) then
                        NNBB2=IBCONN(J)
                        if(NNBB2.gt.0) then ! This check is needed 
                        !occasionally but not always when the barrier is
                        ! split for parallel computation
                           ! ETA2(I) = ETA2(NNBB2)
                           !HTOT=DP(I)+ETA2(I)
                           HTOT = DP(I) + MAX(ETA2(I), ETA2(NNBB2))
                           if(HTOT.gt.HOFF) then
                              !ETA2(I) = MAX(ETA2(I), ETA2(NNBB2))
                              TotArea1 = computeTotAreaAtNode(I)
                              TotArea2 = computeTotAreaAtNode(NNBB2)
                              ETA2AVG = (ETA2(I)*TotArea1 + ETA2(NNBB2)*
     &                                  TotArea2)/(TotArea1 + TotArea2)
                              ETA2(I) = ETA2AVG
                              ETA2(NNBB2) = ETA2AVG
                              nnodecode(I)=1
                              nodecode(I)=1
                              ncchange=ncchange+1 !set to 0 in gwce
                           endif
                        endif
                     endif
                  endif
               endif
            end do
         endif

        ! Use Message-Passing to update at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(nodecode,nnodecode,2)
cjjw#endif


        ! WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
         do I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)

        ! Nodal Wetting Criteria W1:
        ! This depends on changes that occurred in D1
       
            ! Count number of nodes along ibtype=64 boundary 10/16/2023 sb
            if(NFLUXIB64_GBL.gt.0) then
               IBTYPE64FLOODPLAINCNT = 0
               if(onIBTYP64Floodplain(NM1)) then
                  IBTYPE64FLOODPLAINCNT = IBTYPE64FLOODPLAINCNT + 1
               endif
               if(onIBTYP64Floodplain(NM2)) then
                  IBTYPE64FLOODPLAINCNT = IBTYPE64FLOODPLAINCNT + 1
               endif
               if(onIBTYP64Floodplain(NM3)) then
                  IBTYPE64FLOODPLAINCNT = IBTYPE64FLOODPLAINCNT + 1
               endif
            endif

            NCTOT=nodecode(NM1)+nodecode(NM2)+nodecode(NM3)
            if((NCTOT.eq.2) .AND. (IBTYPE64FLOODPLAINCNT.NE.3)) then  
            ! (IBTYPE64FLOODPLAINCNT.NE.3) is added to avoid excessive
            !  wetting when all three nodes of an element on a floodplain
            !  are on a IBTYPE=64 boundary. 10/16/2023 sb
               ETAN1=ETA2(NM1)
               ETAN2=ETA2(NM2)
               ETAN3=ETA2(NM3)
               if(level0)then
                  HTOTN1 = gridDepthVertETA2(NM1)
                  HTOTN2 = gridDepthVertETA2(NM2)
                  HTOTN3 = gridDepthVertETA2(NM3)
               else
                  HTOTN1=DP(NM1)+ETA2(NM1)
                  HTOTN2=DP(NM2)+ETA2(NM2)
                  HTOTN3=DP(NM3)+ETA2(NM3)
               endif
               NMI = 0
               if((nodecode(NM1).eq.1).AND.(nodecode(NM2).eq.1)) then
                  if((HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF)) then
                     NM123=NM1
                     if(ETA2(NM2).gt.ETA2(NM1)) NM123=NM2
                     DELDIST=SQRT((y(NM3)-y(NM123))**2.D0
     &                    +(X(NM3)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM3)
                     ! jgf50.60.18: Prevent problems if DELETA is negative
                     if (DELETA.lt.0.d0) DELETA = 0.d0
                     if(level0)then
                         H1 = wetDepthVertETA2(NM123)
                      else
                         H1=ETA2(NM123)+DP(NM123)
                      endif
c  RJW merged from Casey 071219: Added the following logic for 3D friction
c  RJW modified the following for 3D friction
                     if(C2DDI)then
                        ! sb46.28sb02
                        !Convert Manning's N to Cd, if necessary.
                        if (LoadManningsN) then
                           if(level0)then
                             if(subgridVertList(NM123).eq.1)then
                                if(level1)then
                                   FRIC(NM123) = cmfVertETA2(NM123)
                                else
                                   FRIC(NM123) = cfVertETA2(NM123)
                                endif
                             else
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                   **(1.d0/3.d0) )
                             endif
                          else
                             FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                         /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                  **(1.d0/3.d0) )
                          endif
                          if(FRIC(NM123).lt.BFCdLLimit) then
                            FRIC(NM123) = BFCdLLimit
                          endif
                        endif
! DMW: 202207 Reformulate Velocity by taking sqrt of u instead of using velmin
                        if (directvelWD.AND.useHF) then
                           VEL=sqrt(g*H1*DELETA/DELDIST)
                        elseif (directvelWD) then
                           VEL=sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                        else
                           TKWET = MAX(FRIC(NM123)*(ifLINBF+(velmin/H1)*
     &                                  (ifNLBF+ifHYBF*
     &                                  (1.D0+(HBREAK/H1)**FTHETA)**
     &                                  (FGAMMA/FTHETA))),
     &                                  0.0001d0)
                           VEL=G*(DELETA/DELDIST)/TKWET
                        endif   
                     elseif(C3D)then
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coeff
                       if(LoadZ0B_var) then
                           Z0B1 = Z0B_var(NM123)
                        elseif (LoadManningsN) then
                     Z0B1 = ( DP(NM123)+ifNLFA*ETA2(NM123) )*
     &                      exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+ifNLFA*ETA2(NM123))**
     &               (1.0D0/6.0D0) )/(ManningsN(NM123)*sqrt(g)) ) ))
                        else
                           Z0B1 = Z0B
                       endif
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                    endif

                     if(VEL.gt.velmin) then
                     ! third node met criteria and is also wet
                        nnodecode(NM3)=1
                        NMI = NM3
c RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                        if(C2DDI)then

                           TK_TMP(NM123) = MIN(
     &                     FRIC(NM123)*(ifLINBF+(TKNF*VEL/H1)*
     &                        (ifNLBF+ifHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA))), 
     &                     2.0D0/DTDP)
                           TK(NM123) = MAX(TK_TMP(NM123),TK(NM123))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123)
                           if (NM123 .eq. NM2) then
                             TK_TMP(NM1) = TK(NM123)
                             TK(NM1) = MAX(TK_TMP(NM1),TK(NM1))
                             TKM(1:2,NM1) = TK(NM1)
                             if (TK_TMP(NM1).lt.TK(NM1))then
                               TK(NM123) = TK(NM1)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           else
                             TK_TMP(NM2) = TK(NM123)
                             TK(NM2) = MAX(TK_TMP(NM2),TK(NM2))
                             TKM(1:2,NM2) = TK(NM2)
                             if (TK_TMP(NM2).lt.TK(NM2))then
                               TK(NM123) = TK(NM2)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           endif
                        elseif(C3D)then
                           if(ISLIP.eq.0)then
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*real(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM3)=EVTOT(1)*real(DUDS)
                              BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                           endif
                           if(ISLIP.NE.0)then
                              if(ISLIP.eq.1)then
                                 KSLIP=KP
                              endif
                              if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              endif
                              BSX1(NM123)=KSLIP*real(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM3)=KSLIP*real(Q(NM123,1))
                              BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                           endif
!Casey 220207 from Rosemary : Update sal/temp at newly wet node
                           if((CBAROCLINIC).and.(IDEN.eq.4))then
                              do k = 1,NFEN
                               SAL(NM3,k) =0.5D0*(SAL(NM1,k)+SAL(NM2,k))
                               TEMP(NM3,k)=0.5D0*(TEMP(NM1,k)+TEMP(NM2,k))
                              end do
                           endif
                        endif
                     endif
                  endif

               elseif((nodecode(NM2).eq.1).AND.(nodecode(NM3).eq.1))
     &          then
                  if((HTOTN2.GE.HOFF).AND.(HTOTN3.GE.HOFF)) then
                     NM123=NM2
                     if(ETA2(NM3).gt.ETA2(NM2)) NM123=NM3
                     DELDIST=SQRT((Y(NM1)-Y(NM123))**2.D0
     &                    +(X(NM1)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM1)
                     ! jgf50.60.18: Prevent problems if DELETA is neg
                     if (DELETA.lt.0.d0) DELETA = 0.d0
                     if(level0)then
                        H1 = wetDepthVertETA2(NM123)
                     else
                        H1=ETA2(NM123)+DP(NM123)
                     endif
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     if(C2DDI)then
Csb46.28sb02
                        !Convert Manning's N to Cd, if necessary.
                        if (LoadManningsN) then
                          if(level0)then
                              if(subgridVertList(NM123).eq.1)then
                                if(level1)then
                                   FRIC(NM123) = cmfVertETA2(NM123)
                                else
                                   FRIC(NM123) = cfVertETA2(NM123)
                                endif
                             else
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                   **(1.d0/3.d0) )
                             endif
                          else
                             FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                         /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                  **(1.d0/3.d0) )
                          endif
                          if(FRIC(NM123).lt.BFCdLLimit) then
                            FRIC(NM123) = BFCdLLimit
                          endif
                        endif

!-- DMW Reformulate Velocity by taking sqrt of u instead of using velmin 202207                  
                        if (directvelWD.AND.useHF) then
                           VEL=sqrt(g*H1*DELETA/DELDIST)
                        elseif (directvelWD) then
                           VEL=sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                        else
                           TKWET = MAX(FRIC(NM123)*(ifLINBF+(velmin/H1)*
     &                                  (ifNLBF+ifHYBF*
     &                                  (1.D0+(HBREAK/H1)**FTHETA)**
     &                                  (FGAMMA/FTHETA))),0.0001d0)
                           VEL=G*(DELETA/DELDIST)/TKWET
                        endif

                     elseif(C3D)then
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                        if(LoadZ0B_var) then
                           Z0B1=Z0B_var(NM123)
                        elseif (LoadManningsN) then
                     Z0B1 = ( DP(NM123)+ifNLFA*ETA2(NM123) )*
     &                      exp(-(1.0D0+
     &                      ( (0.41D0*( DP(NM123)+ifNLFA*ETA2(NM123))**
     &                      (1.0D0/6.0D0))/(ManningsN(NM123)*sqrt(g)))))
                        else
                           Z0B1=Z0B
                        endif
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     endif

                     if(VEL.gt.velmin) then
                        nnodecode(NM1)=1
                        NMI = NM1
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        if(C2DDI)then
                           TK_TMP(NM123) = MIN(
     &                     FRIC(NM123)*(ifLINBF+(TKNF*VEL/H1)*
     &                        (ifNLBF+ifHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA))), 
     &                     2.0D0/DTDP)
                           TK(NM123) = MAX(TK_TMP(NM123),TK(NM123))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123) 
                           if (NM123 .eq. NM2) then
                             TK_TMP(NM3) = TK(NM123)
                             TK(NM3) = MAX(TK_TMP(NM3),TK(NM3))
                             TKM(1:2,NM3) = TK(NM3)
                             if (TK_TMP(NM3).lt.TK(NM3))then
                               TK(NM123) = TK(NM3)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           else
                             TK_TMP(NM2) = TK(NM123)
                             TK(NM2) = MAX(TK_TMP(NM2),TK(NM2))
                             TKM(1:2,NM2) = TK(NM2)
                             if (TK_TMP(NM2).lt.TK(NM2))then
                               TK(NM123) = TK(NM2)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           endif
                        elseif(C3D)then
                           if(ISLIP.eq.0)then
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*real(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM1)=EVTOT(1)*real(DUDS)
                              BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                           endif
                           if(ISLIP.NE.0)then
                              if(ISLIP.eq.1)then
                                 KSLIP=KP
                              endif
                              if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              endif
                              BSX1(NM123)=KSLIP*real(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM1)=KSLIP*real(Q(NM123,1))
                              BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                           endif
!Casey 220207 from Rosemary : Again in the 3 wetting criteria
                           if((CBAROCLINIC).and.(IDEN.eq.4))then
                              do k = 1,NFEN
                               SAL(NM1,k) =0.5D0*(SAL(NM2,k)+SAL(NM3,k))
                               TEMP(NM1,k)=0.5D0*(TEMP(NM2,k)+TEMP(NM3,k))
                              end do
                        endif
                        endif

                     endif
                  endif

               elseif((nodecode(NM3).eq.1).AND.(nodecode(NM1).eq.1))
     &          then
                  if((HTOTN3.GE.HOFF).AND.(HTOTN1.GE.HOFF)) then
                     NM123=NM3
                     if(ETA2(NM1).gt.ETA2(NM3)) NM123=NM1
                     DELDIST=SQRT((Y(NM2)-Y(NM123))**2.D0
     &                     +(X(NM2)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM2)
                     ! jgf50.60.18: Prevent problems if DELETA is neg
                     if (DELETA.lt.0.d0) DELETA = 0.d0
                     if(level0)then
                        H1 = wetDepthVertETA2(NM123)
                     else
                        H1=ETA2(NM123)+DP(NM123)
                     endif
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     if(C2DDI)then
Csb46.28sb02
                        !Convert Manning's N to Cd, if necessary.
                        if (LoadManningsN) then
                          if(level0)then
                             if(subgridVertList(NM123).eq.1)then
                                if(level1)then
                                   FRIC(NM123) = cmfVertETA2(NM123)
                                else
                                   FRIC(NM123) = cfVertETA2(NM123)
                                endif
                             else
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                   **(1.d0/3.d0) )
                             endif
                          else
                             FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                         /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                                  **(1.d0/3.d0) )
                          endif
                          if(FRIC(NM123).lt.BFCdLLimit) then
                            FRIC(NM123) = BFCdLLimit
                          endif
                        endif

!-- DMW Reformulate Velocity by taking sqrt of u instead of using velmin 202207                  
                        if (directvelWD.AND.useHF) then
                           VEL=sqrt(g*H1*DELETA/DELDIST)
                        elseif (directvelWD) then
                           VEL=sqrt((g*H1/FRIC(NM123))*(DELETA/DELDIST))
                        else
                           TKWET = MAX(FRIC(NM123)*(ifLINBF+(velmin/H1)*
     &                                  (ifNLBF+ifHYBF*
     &                                  (1.D0+(HBREAK/H1)**FTHETA)**
     &                                  (FGAMMA/FTHETA))),0.0001d0)
                           VEL=G*(DELETA/DELDIST)/TKWET
                        endif
                     elseif(C3D)then
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coeff
                        if(LoadZ0B_var) then
                           Z0B1 = Z0B_var(NM123)
                        elseif (LoadManningsN) then
                     Z0B1 = ( DP(NM123)+ifNLFA*ETA2(NM123) )*
     &                       exp(-(1.0D0+
     &                      ((0.41D0*( DP(NM123)+ifNLFA*ETA2(NM123))**
     &                      (1.0D0/6.0D0))/(ManningsN(NM123)*sqrt(g)))))
                        else
                           Z0B1 = Z0B
                        endif
                      VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     endif

                     if(VEL.gt.velmin) then
                        nnodecode(NM2)=1
                        NMI = NM2
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        if(C2DDI)then                   
                        
                           TK_TMP(NM123) = MIN(
     &                     FRIC(NM123)*(ifLINBF+(TKNF*VEL/H1)*
     &                        (ifNLBF+ifHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA))), 
     &                     2.0D0/DTDP)
                           TK(NM123) = MAX(TK_TMP(NM123),TK(NM123))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123)
                           if (NM123 .eq. NM1) then
                             TK_TMP(NM3) = TK(NM123)
                             TK(NM3) = MAX(TK_TMP(NM3),TK(NM3))
                             TKM(1:2,NM3) = TK(NM3)
                             if (TK_TMP(NM3).lt.TK(NM3))then
                               TK(NM123) = TK(NM3)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           else
                             TK_TMP(NM1) = TK(NM123)
                             TK(NM1) = MAX(TK_TMP(NM1),TK(NM1))
                             TKM(1:2,NM1) = TK(NM1)
                             if (TK_TMP(NM1).lt.TK(NM1))then
                               TK(NM123) = TK(NM1)
                               TKM(1:2,NM123) = TK(NM123)
                             endif
                           endif
                        elseif(C3D)then
                           if(ISLIP.eq.0)then
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*real(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM2)=EVTOT(1)*real(DUDS)
                              BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                           endif
                           if(ISLIP.NE.0)then
                              if(ISLIP.eq.1)then
                                 KSLIP=KP
                              endif
                              if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              endif
                              BSX1(NM123)=KSLIP*real(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM2)=KSLIP*real(Q(NM123,1))
                              BSY1(NM2)=KSLIP*AIMAG(Q(NM123,1))
                           endif
!Casey 220207 from Rosemary : Again in all wetting criteria
                           if((CBAROCLINIC).and.(IDEN.eq.4))then
                              do k = 1,NFEN
                               SAL(NM2,k) =0.5D0*(SAL(NM1,k)+SAL(NM3,k))
                               TEMP(NM2,k)=0.5D0*(TEMP(NM1,k)+TEMP(NM3,k))
                              end do
                           endif
                        endif

                     endif
                  endif
               endif
            endif
       
       !---------------------------
       ! Nodal Wetting Criteria W2a
       !---------------------------
       
            NBNCTOT=NIBnodecode(NM1)+NIBnodecode(NM2)+NIBnodecode(NM3)
            NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
            NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
            NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

         end do
    
         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
         CALL UPDATEI(nnodecode,NIBCNT,2)
!-- DMW 202204 Update FRIC, TK, and ETA2 at ghost nodes as well 
         CALL UPDATER(FRIC,TK,ETA2,3)
         tmp1 = TKM(1,:)
         tmp2 = TKM(2,:)
         tmp3 = TKM(3,:)
         CALL UPDATER(tmp1,tmp2,tmp3,3)
         TKM(1,:) = tmp1
         TKM(2,:) = tmp2
         TKM(3,:) = tmp3 
#endif


        ! ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
        ! tcm v50.66.01 -- This is an additional test for wetting only
        ! when time varying bathymetry is used and is only performed
        ! during the period of bathymetry evolution. 
         if ((NDDT.NE.0).AND.(IT.LE.BTIME_END+1) ) then
            do I=1,NE
               NM1=NM(I,1)
               NM2=NM(I,2)
               NM3=NM(I,3)
               NCTOT=nodecode(NM1)+nodecode(NM2)+nodecode(NM3)
               if(NCTOT.lt.3) then   !If not wet from previous time step
                  NCTOT=nnodecode(NM1)+nnodecode(NM2)+nnodecode(NM3)
                  if(NCTOT.lt.3) then !if not yet wetted this time step
                     ETAN1=ETA2(NM1)
                     ETAN2=ETA2(NM2)
                     ETAN3=ETA2(NM3)
                     HTOTN1=DP(NM1)+ETA2(NM1)
                     HTOTN2=DP(NM2)+ETA2(NM2)
                     HTOTN3=DP(NM3)+ETA2(NM3)

                     ! if all nodes have a depth greater than or equal to 
                     ! hoff = 1.2*H0, then make the element wet
                     if( (HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF).AND.
     &                   (HTOTN3.GE.HOFF) ) then
                     !THE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                        !  Make Node 1 Wet and set parameters
                        if(nnodecode(NM1).NE.1) then  !node 1
                           nnodecode(NM1)=1
                           NM123 = NM1
                           if(C2DDI)then
                              !onvert Manning's N to Cd, if necessary.
                              if (LoadManningsN) then
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                              /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 if(FRIC(NM123).lt.BFCdLLimit) then
                                    FRIC(NM123) = BFCdLLimit
                                 endif
                              endif
                           endif
                           VEL=velmin
                           H1 = HTOTN1
                           if(C2DDI)then
                              TK(NM123)=FRIC(NM123)*(ifLINBF+(VEL/H1)*
     &                          (ifNLBF+ifHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           elseif(C3D)then
                              if(ISLIP.eq.0)then
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                              /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*real(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM3)=EVTOT(1)*real(DUDS)
                                 BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                              endif
                              if(ISLIP.NE.0)then
                                 if(ISLIP.eq.1)then
                                    KSLIP=KP
                                 endif
                                 if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                     if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                     if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                                 endif
                                 BSX1(NM123)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                 BSX1(NM3)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                              endif
!Casey 220207 from Rosemary : In ALL wetting criteria
                              if((CBAROCLINIC).and.(IDEN.eq.4))then
                                 do k = 1,NFEN
                                    SAL(NM3,k)=0.5D0*(SAL(NM1,k)+
     &                                         SAL(NM2,k))
                                    TEMP(NM3,k)=0.5D0*(TEMP(NM1,k)+
     &                                          TEMP(NM2,k))
                                 end do
                           endif
                           endif
                        endif  !end node 1

                        !  Make Node 2 Wet and set parameters
                        if (nnodecode(NM2).NE.1) then
                           nnodecode(NM2) = 1
                           NM123=NM2
                           if(C2DDI)then
                           !Convert Manning's N to Cd, if necessary.
                              if (LoadManningsN) then
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 if(FRIC(NM123).lt.BFCdLLimit) then
                                    FRIC(NM123) = BFCdLLimit
                                 endif
                              endif
                           endif
                           VEL = velmin
                           H1 = HTOTN2
                           if(C2DDI)then
                              TK(NM123)=FRIC(NM123)*(ifLINBF+(VEL/H1)*
     &                          (ifNLBF+ifHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           elseif(C3D)then
                              if(ISLIP.eq.0)then
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                              /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*real(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM1)=EVTOT(1)*real(DUDS)
                                 BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                              endif
                              if(ISLIP.NE.0)then
                                 if(ISLIP.eq.1)then
                                    KSLIP=KP
                                 endif
                                 if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                 (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                    if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                    if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                                 endif
                                 BSX1(NM123)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                 BSX1(NM1)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                              endif
!Casey 220207 from Rosemary : And again in the next two wetting criteria
                              if((CBAROCLINIC).and.(IDEN.eq.4))then
                                 do k = 1,NFEN
                                    SAL(NM1,k)=0.5D0*(SAL(NM2,k)+
     &                                         SAL(NM3,k))
                                    TEMP(NM1,k)=0.5D0*(TEMP(NM2,k)+
     &                                          TEMP(NM3,k))
                                 end do
                           endif
                           endif
                        endif  !node 2

                        !  Make Node 3 Wet and set parameters
                        if(nnodecode(NM3).NE.1) then
                           nnodecode(NM3)=1
                           NM123 = NM3
                           if(C2DDI)then
                            !Convert Manning's N to Cd, if necessary.
                              if (LoadManningsN) then
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                           /( ( DP(NM123)+ifNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 if(FRIC(NM123).lt.BFCdLLimit) then
                                    FRIC(NM123) = BFCdLLimit
                                 endif
                              endif
                           endif
                           VEL=velmin
                           H1 = HTOTN3
                           if(C2DDI)then
                              TK(NM123)=FRIC(NM123)*(ifLINBF+(VEL/H1)*
     &                           (ifNLBF+ifHYBF*
     &                           (1.D0+(HBREAK/H1)**FTHETA)**
     &                           (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           elseif(C3D)then
                              if(ISLIP.eq.0)then
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                             /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*real(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM2)=EVTOT(1)*real(DUDS)
                                 BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                              endif
                              if(ISLIP.NE.0)then
                                 if(ISLIP.eq.1)then
                                    KSLIP=KP
                                 endif
                                 if(ISLIP.eq.2)then
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+ifNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                     if(KSLIP.gt.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                     if(KSLIP.lt.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                                 endif
                                 BSX1(NM123)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                 BSX1(NM2)=KSLIP*real(Q(NM123,1))
                                 BSY1(NM2)=KSLIP*AIMAG(Q(NM123,1))
                              endif
!Casey 220207 from Rosemary : In ALL wetting criteria
                              if((CBAROCLINIC).and.(IDEN.eq.4))then
                                 do k = 1,NFEN
                                    SAL(NM2,k)=0.5D0*(SAL(NM1,k)+
     &                                          SAL(NM3,k))
                                    TEMP(NM2,k)=0.5D0*(TEMP(NM1,k)+
     &                                          TEMP(NM3,k))
                                 end do
                              endif
                           endif
                        endif !node3

                     endif  !ALL DEPTHS GREATER THAN HOFF
                  endif  !if nnodecode SUM LESS THAN 3
               endif   ! if nodecode SUM LESS THAN 3
            end do  !LOOP OVER ELEMENTS

            if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
            CALL UPDATEI(nnodecode,NIBCNT,2)
#endif
 
         endif !if TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME
         !END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY

         !---------------------------
         ! Nodal Wetting Criteria W2b
         !---------------------------
         ! Check for adjacent nodes and force nodes wet when attached
         ! to receiving barrier nodes
         do I=1,NP
            if((NIBCNT(I).gt.0).AND.(nnodecode(I).eq.0)) then
               nnodecode(I)=1
            endif
         end do

cjjwC     Use Message-Passing to update nnodecode at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(nnodecode,IDUMY,1)
cjjw#endif

       ! END WET/DRY SECTION - PART 2

       ! START WET/DRY SECTION  - PART 3

         ! Elemental drying criteria DE1
         ! This is an elemental check section designed to avoid 
         ! artificial wetting of control sections. All elements where
         ! downhill flow originates from a barely wet node into wet nodes
         ! are forced inactive; the only exception is receiving overtopped
         ! barrier nodes
         do I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NBNCTOT=NIBCNT(NM1)*NIBCNT(NM2)*NIBCNT(NM3)
            if(NBNCTOT.eq.0) then   !No barrier/pipe receiving nodes in this elem
               ETAN1=ETA2(NM1)
               ETAN2=ETA2(NM2)
               ETAN3=ETA2(NM3)
               if(level0)then
                  HTOTN1 = gridDepthVertETA2(NM1)
                  HTOTN2 = gridDepthVertETA2(NM2)
                  HTOTN3 = gridDepthVertETA2(NM3)
               else
                  HTOTN1=DP(NM1)+ETA2(NM1)
                  HTOTN2=DP(NM2)+ETA2(NM2)
                  HTOTN3=DP(NM3)+ETA2(NM3)
               endif
#ifdef SB_WETDRY
C...Find the heighest point on the bed in the element.  sb v46.28.sb05.06 11/01/2006
               if(DP(NM1).LE.DP(NM2).AND.DP(NM1).LE.DP(NM3)) then
                  DPMIN = DP(NM1)
               else if(DP(NM2).LE.DP(NM3).AND.DP(NM2).LE.DP(NM1)) then
                  DPMIN = DP(NM2)
               else if(DP(NM3).LE.DP(NM1).AND.DP(NM3).LE.DP(NM2)) then
                  DPMIN = DP(NM3)
               endif
#endif
#ifndef SB_WETDRY
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
C...ABC pattern
                  if((ETAN1.GE.ETAN2).AND.(ETAN2.gt.ETAN3)) then
                    if((HTOTN1.lt.HOFF).OR.(HTOTN2.lt.HOFF)) noff(I)=0
                  endif
                  if((ETAN2.GE.ETAN3).AND.(ETAN3.gt.ETAN1)) then
                    if((HTOTN2.lt.HOFF).OR.(HTOTN3.lt.HOFF)) noff(I)=0
                  endif
                  if((ETAN3.GE.ETAN1).AND.(ETAN1.gt.ETAN2)) then
                     if((HTOTN3.lt.HOFF).OR.(HTOTN1.lt.HOFF)) noff(I)=0
                  endif
C...ACB pattern
                  if((ETAN1.GE.ETAN3).AND.(ETAN3.gt.ETAN2)) then
                    if((HTOTN1.lt.HOFF).OR.(HTOTN3.lt.HOFF)) noff(I)=0
                  endif
                  if((ETAN2.GE.ETAN1).AND.(ETAN1.gt.ETAN3)) then
                    if((HTOTN2.lt.HOFF).OR.(HTOTN1.lt.HOFF)) noff(I)=0
                  endif
                  if((ETAN3.GE.ETAN2).AND.(ETAN2.gt.ETAN1)) then
                    if((HTOTN3.lt.HOFF).OR.(HTOTN2.lt.HOFF)) noff(I)=0
                  endif
               endif
#else
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
C...ABC pattern
                  if((ETAN1.GE.ETAN2).AND.(ETAN2.gt.ETAN3)) then
                     if((HTOTN1.lt.HOFF)) noff(I)=0
                     if((ETAN1-ETAN2).lt.(ETAN2-ETAN3)) then
                        if(HTOTN2.lt.HOFF) noff(I)=0
                     endif
                  endif
                  if((ETAN2.GE.ETAN3).AND.(ETAN3.gt.ETAN1)) then
                     if((HTOTN2.lt.HOFF)) noff(I)=0
                     if((ETAN2-ETAN3).lt.(ETAN3-ETAN1)) then
                        if(HTOTN3.lt.HOFF) noff(I)=0
                     endif
                  endif
                  if((ETAN3.GE.ETAN1).AND.(ETAN1.gt.ETAN2)) then
                     if((HTOTN3.lt.HOFF)) noff(I)=0
                     if((ETAN3-ETAN1).lt.(ETAN1-ETAN2)) then
                        if(HTOTN1.lt.HOFF) noff(I)=0
                     endif
                  endif
C...ACB pattern
                  if((ETAN1.GE.ETAN3).AND.(ETAN3.gt.ETAN2)) then
                     if((HTOTN1.lt.HOFF)) noff(I)=0
                     if((ETAN1-ETAN3).lt.(ETAN3-ETAN2)) then
                        if(HTOTN3.lt.HOFF) noff(I)=0
                     endif
                  endif
                  if((ETAN2.GE.ETAN1).AND.(ETAN1.gt.ETAN3)) then
                     if((HTOTN2.lt.HOFF)) noff(I)=0
                     if((ETAN2-ETAN1).lt.(ETAN1-ETAN3)) then
                        if(HTOTN1.lt.HOFF) noff(I)=0
                     endif
                  endif
                  if((ETAN3.GE.ETAN2).AND.(ETAN2.gt.ETAN1)) then
                     if((HTOTN3.lt.HOFF)) noff(I)=0
                     if((ETAN3-ETAN2).lt.(ETAN2-ETAN1)) then
                        if(HTOTN2.lt.HOFF) noff(I)=0
                     endif
                  endif
               endif !noffActive.eqv..true.
#endif

#ifdef SB_WETDRY
!>An element is set to be dry if it is determined to be a flooding type
!>wetting element.  An element is a flooding type wetting element if
!>the bed elevation at the node with the biggest water column height
!>is lower than the heighest point on the bed in the element.
!>sb v46.28.sb05.06 11/01/2006
!>This is applied only when noff flag of the element at the previous time
!>step is 0, which means that this logic works to prevent an element from
!>re-wetting. sb v46.52.03
               if (noffActive.eqv..true.) then 
                  if(noffold(I).eq.0) then
                     if(HTOTN1.GE.HTOTN2.AND.HTOTN1.GE.HTOTN3) then
                        if(ETAN1.lt.(-DPMIN+H0)) noff(I) = 0
                     endif
                     if(HTOTN2.GE.HTOTN3.AND.HTOTN2.GE.HTOTN1) then
                        if(ETAN2.lt.(-DPMIN+H0)) noff(I) = 0
                     endif
                     if(HTOTN3.GE.HTOTN1.AND.HTOTN3.GE.HTOTN2) then
                        if(ETAN3.lt.(-DPMIN+H0)) noff(I) = 0
                     endif
                 endif
              endif
#endif
            endif
         end do

#ifdef SB_WETDRY
!>     ...added by sb on 11/02/2006
!>     
!>     This section is added after we realize that it's not possible
!>     to compute a correct flow going through two elements if
!>     the elements are connected just by one node. i.e., elements need
!>     to share an edge to let the flow go through between the elements.
!>     Therefore, in this section, a node is determined to be dry
!>     if two elements are connected at one node, not sharing an edge.
!>     As it seemed this procedure needed noff information,
!>     although this section changes nnodecode,
!>     I put this section here, rather than the end of PART 2.
!>     I tried setting noff(I) = 0, but it didn't shut down the flow.
!>     
         do I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            if(noff(I).eq.1.AND.
     &         nnodecode(NM1).eq.1.AND.
     &         nnodecode(NM2).eq.1.AND.
     &         nnodecode(NM3).eq.1) then
               do K=1,3
                  NM1=NM(I,K)
                  NM2=NM(I,MOD(K+0,3)+1)
                  NM3=NM(I,MOD(K+1,3)+1)

                  NWETNEI = 0
                  NWETADJ = 0
                  do J=1,MNEI
                     N=NeiTabEle(NM1,J)
                     if(N.eq.0) CYCLE
                     if(N.eq.I) CYCLE

                     NMN1=NM(N,1)
                     NMN2=NM(N,2)
                     NMN3=NM(N,3)
                     if(noff(N).eq.1.AND.
     &                  nnodecode(NMN1).eq.1.AND.
     &                  nnodecode(NMN2).eq.1.AND.
     &                  nnodecode(NMN3).eq.1) then
                        NWETNEI = NWETNEI + 1

                        if((NMN1.eq.NM2.OR.NMN1.eq.NM3).OR.
     &                     (NMN2.eq.NM2.OR.NMN2.eq.NM3).OR.
     &                     (NMN3.eq.NM2.OR.NMN3.eq.NM3)) then
                           NWETADJ = NWETADJ + 1
                        endif
                     endif
                  end do

                  if(NWETNEI.gt.0.AND.NWETADJ.eq.0.AND.
     &               NIBCNT(NM1).eq.0) then
                     nnodecode(NM1) = 0
                  endif
               end do
            endif
         end do
#endif


       ! END WET/DRY SECTION  - PART 3

       ! Nodal Drying Criteria for the ibtype=64 weir boundary
       ! Check the node if it is along the ibtype=64 weir boundary and if 
       ! the weir is submerged. If true and the node does not belong to
       ! at least one wet element, set the nnodecode back to 0.  SB
         if (NFLUXIB64_GBL.gt.0) then
            do I=1,NP
               if(LBArray_Pointer(I).gt.0) then
                  J = LBArray_Pointer(I)
                  if((ISSUBMERGED64(J).gt.0).AND.(nnodecode(I).eq.1)) then
                     NNBB2=IBCONN(J)
                     NWetEle = 0
                     do K=1,NNeighEle(I)
                        IE = NeiTabEle(I,K)
                        NCTOT = nnodecode(NM(IE,1)) + 
     &                          nnodecode(NM(IE,2)) + nnodecode(NM(IE,3))
                        if (NCTOT.eq.3) then
                           NWetEle = NWetEle + 1
                        endif
                     end do
                     if (NWetEle.eq.0) then
                        nnodecode(I) = 0
                     endif
                  endif
               endif
            end do
         endif

      ! Section for condensed nodes
      ! Setting the same nodecode value at grouped condensed nodes   SB
      if(LoadCondensedNodes) then
         do K=1,NListCondensedNodes
            I = ListCondensedNodes(K,1)
            if(I==0) CYCLE
            ! 0) DEFAUlt VALUE
            NNC = nnodecode(I)
            ! 1) FIND UPDATED VALUE 
            do L=2,NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K,L)
               if(nnodecode(J).NE.nodecode(J)) then
                  NNC = nnodecode(J)
                  EXIT ! EXiT THE LOOP if AN UPDATED nodecode IS FOUND
               endif
            end do
            ! 2) DISTRIBUTE IT
            nnodecode(I) = NNC
            do L=2,NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K,L)
               nnodecode(J) = NNC
            end do
         end do
      endif


!>       START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2
!>
!>       Update number of active elements (MJU) and the total area
!>       (TotalArea) connected to a node. If these are zero, the node is
!>       landlocked and should be dried. These depend on nnodecode which
!>       varies during the time step. Also add computation of total area
!>       reciprocal used in various modules. This code was a repeat of
!>       totalAreaCalc in the mesh module, so that subroutine was
!>       modified to use nnodecode instead of nodecode and reused here.
!>       See additional nodes in mesh module. 
!>       Of special note, the necessity of having non-zero mju values is
!>       only required when trying to match incorrect verions of adcirc
!>       where mju was used in place of the totalarea in the momentum
!>       equations; at some point this could be removed. CMS 11/18/25

         call totalAreaCalc()


!>       Start update variables
!>          jjwnote - may have to pass TotalArea and mju as well

         !NCSU subdomain
         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb()
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob()

         ! Message-Passing to update nnodecode at ghost nodes
#ifdef CMPI
         CALL UPDATEI(nnodecode,IDUMY,1)
         CALL UPDATER(TotalArea,DUMY1,DUMY2,1)
#endif


!>     WET/DRY SECTION - PART 5 - RESET nodecode USING nnodecode
!>     Check to see if any wetting occurred & update nodecode
!>     Note, ncchange=0 set near the beginning of GWCE subroutine
         do I=1,NP
            if(nnodecode(I).NE.nodecode(I)) then
               nodecode(I)=nnodecode(I)
               ncchange=ncchange+1
            endif
         end do

     
!>     WET/DRY SECTION - PART 6
!>     Check to see if any noff changed requiring the matrix to be reset
!>     Note, ncchange=0 set near the beginning of GWCE subroutine
     
         do I=1,NE
            if(noff(I).NE.noffold(I)) ncchange=ncchange+1
         end do
     
        ! jgf45.06 If there has been any wetting or drying in any
        ! of the subdomains, the ncchange flag will be activated on all
        ! of the subdomains, to prevent them from getting out of sync
        ! with their MPI calls as some reset the GWCE and others do not.

         if(LoadSubgridBarrier)then
             do I=1,NE
                 if(noff_SG(I).eq.0)then
                     if(noff(I).NE.0)then
                         noff(I) = 0
                         ncchange=ncchange+1
                     endif
                 endif
             end do
         endif
         
         ! DMW 202207 Compute a slope limiting factor if needed
         call computeAlpha()
         
#ifdef CMPI
         !jgf48.4619 implementing Seizo's changes for Lumped, fully
         ! explicit operation. In that case, the GWCE LHS matrix is
         ! recalculated on each individual subdomain that has wetted
         ! or dried, without recourse to MPI, eliminating the need
         ! for the call to the subroutine WetDrySum.
         if ( ILump.eq.0 ) then
            call WetDrySum(ncchange)
         else
            ncchange=ncchange ! jgf48.4619 do nothing
         endif
#endif
     
        ! END WET/DRY SECTION - PART 6
     

!-----------------------------------------------------------------------
      end subroutine computeWettingAndDrying
!-----------------------------------------------------------------------



!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                          P R I V A T E
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Apply a subgrid barrier.
!>
!> @details Assumptions:
!>          * When the subgrid barrier is not being overtopped, the
!>            cross-barrier flow is disabled and the nodes are placed at
!>            the elevation found in the mesh geometry. 
!>          * When the barrier becomes overtopped, the nodes are 
!>            elevated temporarily during the overtopping process.
!>          * When the overtopping concludes, the nodes are returned to
!>            their mesh specified elevations and the cross-barrier flow
!>            is turned off again.
!-----------------------------------------------------------------------
      subroutine applySubgridBarrier()
!-----------------------------------------------------------------------
          use mesh, only : DP,DP_orig,nm,ne
          USE global, only : ETA2,noff,H0,nodecode,ncChange,ilump
          USE nodalattributes, only : SubgridBarrier,
     &                              SubgridBarrierOvertopping
#ifdef CMPI
          use messenger, only : wetDrySum
#endif          
          implicit none
          integer  :: I,NM1,NM2,NM3,NCELE
          real(8) :: SB1,SB2,SB3
          real(8) :: ETABAR,BARHEIGHT
          real(8) :: doT,doT1,doT2,doT3
          real(8) :: E1,E2,E3,ETA
          real(8),parameter :: SBDT = 200D0

          noff_SG(:) = 1

          do I = 1,NE
             NM1=NM(I,1)
             NM2=NM(I,2)
             NM3=NM(I,3)
             SB1=subgridBarrier(NM1,1)
             SB2=subgridBarrier(NM2,1)
             SB3=subgridBarrier(NM3,1)

             NCELE = nodecode(NM1)*nodecode(NM2)*nodecode(NM3)*noff(I)
             if(NCELE.eq.0) cycle

             if(SB1.lt.DP_ORIG(NM1).and.
     &          SB2.lt.DP_ORIG(NM2).and.
     &          SB3.lt.DP_ORIG(NM3))then

                E1 = ETA2(NM1)
                E2 = ETA2(NM2)
                E3 = ETA2(NM3)
                ETA = MAX(E1,E2,E3)
                doT1 = E1*nodecode(NM1)+SB1
                doT2 = E2*nodecode(NM2)+SB2
                doT3 = E3*nodecode(NM3)+SB3
                doT  = MAX(doT1,doT2,doT3)

                if(doT.gt.-1.2D0*H0)then
                    DP(NM1) = MAX(SB1,DP(NM1)+SB1/SBDT)
                    DP(NM2) = MAX(SB2,DP(NM2)+SB2/SBDT)
                    DP(NM3) = MAX(SB3,DP(NM3)+SB3/SBDT)
                    SubgridBarrierOvertopping(NM1) = .true.
                    SubgridBarrierOvertopping(NM2) = .true.
                    SubgridBarrierOvertopping(NM3) = .true.
                else
                    DP(NM1) = MIN(DP_ORIG(NM1),DP(NM1)-SB1/SBDT)
                    DP(NM2) = MIN(DP_ORIG(NM2),DP(NM2)-SB2/SBDT)
                    DP(NM3) = MIN(DP_ORIG(NM3),DP(NM3)-SB3/SBDT)
                    noff_SG(I) = 0
                    SubgridBarrierOvertopping(NM1) = .false.
                    SubgridBarrierOvertopping(NM2) = .false.
                    SubgridBarrierOvertopping(NM3) = .false.
                endif
             else
                SubgridBarrierOvertopping(NM1) = .false.
                SubgridBarrierOvertopping(NM2) = .false.
                SubgridBarrierOvertopping(NM3) = .false.
             endif
          enddo

!-----------------------------------------------------------------------
      end subroutine applySubgridBarrier
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Compute slope limiting factor for elevation gradients.
!>
!> @author DWM 202207  need institution and email
!>
!> @details Compute a slope limiting factor for scaling down the 
!>          elevation gradient if the slope of the water surface aligned
!>          with the bathy slope becomes too large ( slope > slim).
!-----------------------------------------------------------------------
      subroutine computeAlpha()
!-----------------------------------------------------------------------
         use mesh, only : DP,NM,NE,sfac,SFacEle,SFMYEle,SFMXEle,FDXE,FDYE,
     &    areas         
         use global, only : eta2, noff, nodecode, alphal, ifSFM, slim
         use subgrid, only : level0, gridDepthVertETA2
         implicit none
         integer  :: I,NM1,NM2,NM3,NC1,NC2,NC3,NCELE
         real(8) :: sfacavg, SFmxAvg, SFmyAvg, DPN1, DPN2, DPN3
         real(8) :: DETAM, DPM, sfdxfac, sfdyfac, hn1, hn2, hn3
         real(8) :: detadxa, detadya, dpdxa, dpdya
         real(8) :: fdx1, fdx2, fdx3, fdy1, fdy2, fdy3
         real(8) :: etaN1, etaN2, etaN3
         real(8) :: DDTMP  !tcm 20250503
         
         ALPHAL(:) = 1.D0
         
         do I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NC1=nodecode(NM1)
            NC2=nodecode(NM2)
            NC3=nodecode(NM3)
            NCEle=NC1*NC2*NC3*noff(I)
            if (NCEle.eq.1) then
               if(level0) then
                 HN1 = gridDepthVertETA2(NM1)
                 HN2 = gridDepthVertETA2(NM2)
                 HN3 = gridDepthVertETA2(NM3)
               else
                 HN1 = ETA2(NM1) + DP(NM1)
                 HN2 = ETA2(NM2) + DP(NM2)
                 HN3 = ETA2(NM3) + DP(NM3)
               end if

               if((HN1.LE.5.D-1).OR.(HN2.LE.5.D-1).OR.(HN3.LE.5.D-1) ) then
                  SFacAvg = SFacEle(I)
                  SFmxAvg = SFMXEle(I)  
                  SFmyAvg = SFMYEle(I)  
                  sfdxfac = (1 - ifSFM)*SFacAvg + ifSFM*SFmxAvg ;
                  sfdyfac = (1 - ifSFM)*1.0D0 + ifSFM*SFmyAvg ;
               
                  FDX1 = FDXE(1,I)*sfdxfac  !b1
                  FDX2 = FDXE(2,I)*sfdxfac  !b2
                  FDX3 = FDXE(3,I)*sfdxfac  !b3
                  FDY1 = FDYE(1,I)*sfdyfac  !a1
                  FDY2 = FDYE(2,I)*sfdyfac  !a2
                  FDY3 = FDYE(3,I)*sfdyfac  !a3
               
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
                  DPM = SQRT(DPDXA*DPDXA + DPDYA*DPDYA)/AREAS(I)

                  if (DPM .gt. 1D-5) then
                     DDTMP = DETADXA*(-DPDXA) + DETADYA*(-DPDYA)  !tcm 20250503
                     if (ABS(DDTMP) .gt. 1D-12) then   !tcm 20250503 added test and else condition
                        ALPHAL(I) = SLIM*(SQRT(DPDXA*DPDXA + DPDYA*DPDYA)/
     &                     ((1.D0/AREAS(I))*( DDTMP )))  !tcm changed to DDTMP
                      else
                         ALPHAL(I) = 1.D0
                      endif
                  else
                     DETAM = SQRT(DETADXA*DETADXA + DETADYA*DETADYA)/AREAS(I)
                     if (ABS(DETAM) .gt. 1D-12) then   !tcm 20250503 added test and else condition
                        ALPHAL(I) = SLIM/DETAM
                     else
                        ALPHAL(I) = 1.D0
                     endif
                  endif
                  
                  if (ALPHAL(I) .LE. 0.D0) ALPHAL(I) = 1.D0
                  if (ALPHAL(I) .GE. 1.D0) ALPHAL(I) = 1.D0
               endif
                  
            endif 
         end do

!-----------------------------------------------------------------------
      end subroutine computeAlpha
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      real(8) function computeTotAreaAtNode(INODE)
!-----------------------------------------------------------------------
      use mesh, only : NE,NNeighEle,NeiTabEle,Areas
      implicit none
      integer, intent(in) :: INODE
      integer :: I,IE
      real(8) :: TotArea

      TotArea = 0.d0
      do I=1,NNeighEle(INODE)
         IE = NeiTabEle(INODE,I)
         TotArea = TotArea + Areas(IE)
      enddo
      computeTotAreaAtNode = 0.5d0*TotArea
      return
!-----------------------------------------------------------------------
      end function computeTotAreaAtNode
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end module wetDry
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

