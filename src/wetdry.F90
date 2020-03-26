!C***********************************************************************
!C                     M O D U L E   W E T   D R Y      
!C***********************************************************************
!C     Executes wetting and drying algorithm.
!C***********************************************************************
!C Logical Variable List (default value .FALSE., set in global.f)       *
!C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!C     C2DDI            - 2D Depth Integrated model run                 *
!C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!C   See header.f for a summary history of code modifications.          *
!C***********************************************************************
!CWET...
!CWET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
!CWET...
!CWET...NOTE:NNODECODE is a working variable that can change more than once
!CWET...               during a time step
!CWET...     NNODECODE = 0 for a dry node
!CWET...     NNODECODE = 1 for a wet node
!CWET...     NODECODE  - is a more static version of NNODECODE that is reconciled
!CWET...                 once and for all at the end time step
!CWET...
!CWET...
!CWET...        (   DRYING CRITERIA   )
!CWET...
!CWET...A node should be dry under two conditions.
!CWET...D1.) If the total water depth falls below H0.
!CWET .......Note: if the total water depth falls below 0.8*H0, the surface elevation
!CWET........is lifted up so that the total water depth = 0.8*H0.
!CWET......
!CWET...D2.) If the node is connected to only nonfunctioning (dry) elements.  In
!CWET........this case the node is dried due to becoming landlocked.
!CWET........Note: this criteria is applied after all other wetting and drying criteria
!CWET...
!CWET...An element should be dry under the following conditions.
!CWET...DE3.) This is an elemental check section designed to avoid artificial wetting of
!CWET.........of control sections
!CWET.........All elements where downhill flow originates from a barely wet node
!CWET.........(defined as 1.2*H0) into wet nodes are forced inactive; the only exception
!CWET......... is receiving overtopped barrier nodes
!CWET...
!CWET...        (   WETTING CRITERIA   )
!CWET...
!CWET...A node should be wet under two conditions.
!CWET...W1.) If 2 nodes on an element are wet and one is dry, wet the dry node
!CWET........if the water level at one of the wet nodes is greater than the
!CWET........water level at the dry node and the steady state velocity that
!CWET........would result from a balance between the water level gradient and
!CWET........bottom friction would yield a velocity > VELMIN.
!CWET........Note that the criteria outlined in DE3 must also be satified before
!CWET........the node is allowed to wet
!CWET...
!CWET...W2.) If an element has a node lying on a receiving internal barrier boundary or
!CWET......specified discharge boundary that is actively discharging flow into the
!CWET......domain at that node, all nodes in this element must stay wet.
!CWET...
!CWET...
!CWET...        (  VELOCITY BOUNDARY CONDITION  )
!CWET...
!CWET...Either a natural or essential boundary condition can be used as a velocity
!CWET...boundary condition in the momentum equation solution along a wet/dry boudary
!CWET...To use a natural boundary condition, do nothing along the wet/dry interface.
!CWET...To use an essential, no velocity boundary condition, identify the nodes along
!CWET...the wet/dry interface and zero out the velocity at the nodes.  Interface nodes
!CWET...can easily be identified by comparing the number of active elements a node is
!CWET...connected to (MJU) to the total number of elements a node is connected to (NODELE).
!CWET...If MJU < NODELE for any node, it must lie along the wet/dry interface.  See
!CWET...further comments at the end of the momentum equation solution section.
!CWET...
      MODULE WetDry
      use sizes, only : sz, mne, mnp
      use global, only : noff, nodecode, nnodecode, eta2, tk, nolifa,&
          bsx1, bsy1, btime_end, C2DDI, C3D, g, h0, ifnlfa, nddt, &
          nibnodecode, ilump, ncchange, tkm
#ifdef CMPI 
      use global, only : idumy
      use messenger
#endif
      use mesh, only : ne, np, dp, mju, totalArea, nm, x, y, areas
      use nodalattributes, only : BFCdLLimit, fgamma, ftheta, fric,&
           manningsn, hbreak, ifhybf, ifnlbf, iflinbf, loadManningsN,&
           loadZ0B_var, z0b_var
      use global_3dvs, only : a, b, islip, kp, z0b, sigma, evtot, q
      use subdomain, only : subdomainOn, enforceBN, enforceWDcb, enforceWDob

      implicit none

      complex(sz) :: duds !jgf48.50 declare size SZ instead of plain COMPLEX
           
      real(sz) :: habsmin 
      real(sz) :: hoff
      real(sz) :: velmin
     
      integer,allocatable ::    nibcnt(:)
      integer,allocatable ::    noffold(:)

      ! jgf52.08.08: Enable analyst to remove NOFF from wetting and
      ! drying algorithm; .true. by default to conform with prior
      ! ADCIRC versions. 
      LOGICAL :: noffActive = .true.

      integer,allocatable :: temp_NM(:,:)


      contains
!----------------------------------------------------------------------
!                   S U B R O U T I N E
!     C O M P U T E   W E T T I N G   A N D   D R Y I N G
!----------------------------------------------------------------------
!     Determines which nodes should be wet and which should be dry.
!----------------------------------------------------------------------
      subroutine computeWettingAndDrying(it)
      implicit none

      integer, intent(in) :: it ! time step number


      if (nolifa.ne.2) then
         return ! wetting and drying is not active 
      endif

      CALL updateNOFF()

      CALL dryingCriteriaD1()

      CALL wettingCriteriaW1_W2a()
#ifdef CMPI
      CALL updateI(NNODECODE,NIBCNT,2)
#endif
      CALL wettingCriteriaW2b()

      CALL dryingCriteriaDE1()

      CALL dryingCriteriaD2() 
#ifdef CMPI
      CALL updateI(NNODECODE,IDUMY,1)
#endif
      CALL updateNnodecodeAtGhosts() 

      CALL finalizeWettingAndDrying() 

!----------------------------------------------------------------------
      END SUBROUTINE ComputeWettingAndDrying
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!                  S U B R O U T I N E
!     I N I T I A L I Z E   W E T T I N G   A N D   D R Y I N G
!----------------------------------------------------------------------
!     Allocates arrays and initializes constants.
!----------------------------------------------------------------------
      subroutine initializeWettingAndDrying()
      use sizes, only : mne, mnp
      use global, only : h0, nodecode, nnodecode

      IMPLICIT NONE 

      allocate ( nibcnt(mnp) )
      allocate ( noffold(mne))
      ! temp element connectivity table, used for wetting only 
      allocate(temp_NM(mne,3))
      temp_NM = -1 
      nnodecode = 1
      nodecode = 1
      noffold(:) = 1
      habsmin=0.8d0*h0
      hoff=1.2d0*h0

!----------------------------------------------------------------------
      end subroutine initializeWettingAndDrying
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!                       S U B R O U T I N E 
!                        U P D A T E  N O F F
!----------------------------------------------------------------------
    subroutine updateNOFF() 
!----------------------------------------------------------------------
!   Set current values of NOFF to previous ones. 
!----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: I 

    DO I=1,NP
       NIBCNT(I) = 0
    ENDDO
    DO I=1,NE
       NOFFOLD(I)=NOFF(I)
       NOFF(I)=1
    ENDDO
!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!                   S U B R O U T I N E
!            D R Y I N G  C R I T E R I A  D1
!----------------------------------------------------------------------
    subroutine dryingCriteriaD1() 
!----------------------------------------------------------------------
!  D1) Dry the node if the total water depth falls below H0.
!      Note: if the total water depth falls below 0.8*H0, the surface elevation
!      is lifted up so that the total water depth = 0.8*H0 (HABSMIN)
!----------------------------------------------------------------------
    implicit none

    integer :: I
    real(sz) :: HTOT 

    DO I=1,NP
       IF(NODECODE(I).EQ.1) THEN
          HTOT=DP(I)+ETA2(I)
          IF(HTOT.LE.H0) THEN
             IF(HTOT.LT.HABSMIN) ETA2(I)=HABSMIN-DP(I)
             NNODECODE(I)=0
             NODECODE(I)=0
             NCCHANGE=NCCHANGE+1 !NCCHANGE=0 set near beginning of GWCE
          ENDIF
       ENDIF
    ENDDO
!----------------------------------------------------------------------
    end subroutine
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!                   S U B R O U T I N E
!            W E T T I N G  C R I T E R I A  W1
!----------------------------------------------------------------------
    subroutine wettingCriteriaW1_W2a
!----------------------------------------------------------------------
! W1.) If 2 nodes on an element are wet and one is dry, wet the dry node
!      if the water level at one of the wet nodes is greater than the
!      water level at the dry node and the steady state velocity that
!      would result from a balance between the water level gradient and
!      bottom friction would yield a velocity > VELMIN.
!      Note that the criteria outlined in DE3 must also be satified before
!      the node is allowed to wet
!
! W2a.)  If an element has a node lying on a receiving internal barrier boundary or
!      specified discharge boundary that is actively discharging flow into the
!      domain at that node, all nodes in this element must stay wet.
!
!
!      DE3: An element should be dry under the following conditions.
!      This is an elemental check section designed to avoid artificial wetting of
!      of control sections
!      All elements where downhill flow originates from a barely wet node
!      (defined as 1.2*H0) into wet nodes are forced inactive; the only exception
!       is receiving overtopped barrier nodes
!----------------------------------------------------------------------
        IMPLICIT NONE 

        INTEGER :: NM1,NM2,NM3
        INTEGER :: NM123
        INTEGER :: NCELE 
        INTEGER :: NCTOT 
        INTEGER :: NBNCTOT 
        REAL(SZ) :: ETAN1,ETAN2,ETAN3
        REAL(SZ) :: hTotN1,hTotN2,hTotN3      
        REAL(SZ) :: deldist,deleta
        REAL(sz) :: h1
        REAL(SZ) :: TKWET
        REAL(sz) :: vel
        
        INTEGER :: I
        INTEGER :: J
        INTEGER :: JJ
        INTEGER :: NEW


        ! temp element table to be rotated ccw
        DO I=1,NE
            DO J=1,3
                temp_NM(I,J) = NM(I,J)
            ENDDO
        ENDDO

        DO I=1,NE

            ! Rotate element table CCW here 
            DO J=1,3
                
                ! Do not rotate the first time!
                IF(J.GT.1) THEN 
                  temp_NM(I,1:3) = CSHIFT(NM(I,1:3),J-1)
                ENDIF

                NM1=temp_NM(I,1)
                NM2=temp_NM(I,2)
                NM3=temp_NM(I,3)

                NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)

                ! If two nodes of an element are wet
                IF(NCTOT.EQ.2) THEN

                    ! Free surface at the nodes
                    ETAN1=ETA2(NM1)
                    ETAN2=ETA2(NM2)
                    ETAN3=ETA2(NM3)

                    ! Total water depth at the nodes 
                    HTOTN1=DP(NM1)+ETA2(NM1)
                    HTOTN2=DP(NM2)+ETA2(NM2)
                    HTOTN3=DP(NM3)+ETA2(NM3)

                    IF((NODECODE(NM1).EQ.1).AND.(NODECODE(NM2).EQ.1)) THEN

                        IF((HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF)) THEN

                            ! Used to determine the maximum pressure gradient force. 
                            NM123=NM1
                            IF(ETA2(NM2).GT.ETA2(NM1)) THEN 
                                NM123=NM2
                            ENDIF
                            DELDIST=SQRT((y(NM3)-y(NM123))**2.D0 + (X(NM3)-X(NM123))**2.D0)

                            DELETA=ETA2(NM123)-ETA2(NM3)
                            IF (DELETA.lt.0.d0) THEN 
                                 DELETA = 0.d0
                            ENDIF

                            H1=ETA2(NM123)+DP(NM123)

                            ! If using Mannings N friction law 
                            IF (LoadManningsN) THEN
                             
                                FRIC(NM123)=g*ManningsN(NM123)**2.d0/ (DP(NM123)+IFNLFA*ETA2(NM123))**(1.d0/3.d0)

                                IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                ENDIF

                            ENDIF

                            TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*(IFNLBF+IFHYBF*(1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))

                            IF(TKWET.LT.0.0001d0) THEN 
                                TKWET=0.0001d0
                            ENDIF 

                            VEL=G*(DELETA/DELDIST)/TKWET

                            ! Third node of element met wetting criteria met?
                            IF(VEL.GT.VELMIN) THEN
                                ! Yup...wet it!
                                NNODECODE(NM3)=1

                                TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*(IFNLBF+IFHYBF*&
                                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
    
                                !WJP needed for internal_tide parameterization
                                TKM(1:2,NM123) = TK(NM123) 
                            ENDIF ! 
                        ENDIF ! IF THE TWL AT NODES 1 AND 2 ARE > HOFF
                    ENDIF ! IF NODES 1 AND 2 ARE WET 
                ENDIF ! IF TWO NODES ARE WET 

                ! Nodal Wetting Criteria W2a
                NBNCTOT=NIBNODECODE(NM1)+NIBNODECODE(NM2)+NIBNODECODE(NM3)
                NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
                NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
                NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT
            ENDDO ! ELEMENTAL ROTATION LOOP
        ENDDO ! ELEMENTAL LOOP
!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!                   S U B R O U T I N E 
!         W E T T I N G  C R I T E R I A  W2b 
!----------------------------------------------------------------------
    subroutine wettingCriteriaW2B() 
!----------------------------------------------------------------------
!W2b) Nodal Wetting Criteria W2b
!     Check for adjacent nodes and force nodes wet when attached
!     to receiving barrier nodes
!----------------------------------------------------------------------
         IMPLICIT NONE 

         INTEGER :: I 

         DO I=1,NP
            IF((NIBCNT(I).GT.0).AND.(NNODECODE(I).EQ.0)) THEN
               NNODECODE(I)=1
            ENDIF
         ENDDO
!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
    subroutine dryingCriteriaDE1() 
!----------------------------------------------------------------------
! DE1) Elemental drying criteria DE1
!      This is an elemental check section designed to avoid artificial wetting
!      of control sections
!
!      All elements where downhill flow originates from a barely wet node
!      into wet nodes are forced inactive; the only exception is receiving
!      overtopped barrier nodes
!----------------------------------------------------------------------
    IMPLICIT NONE 

    INTEGER :: NM1,NM2,NM3
    INTEGER :: NBNCTOT 
    REAL(SZ) :: ETAN1,ETAN2,ETAN3
    REAL(SZ) :: HTOTN1,HTOTN2,HTOTN3

    INTEGER :: I 


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
            ! jgf52.08.08: Analyst can eliminate noff from
            IF (noffActive.eqv..true.) then
                ! consideration in fort.15 file. 
                IF((ETAN1.GE.ETAN2).AND.(ETAN2.GT.ETAN3)) THEN
                  IF((HTOTN1.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
                ENDIF
                IF((ETAN2.GE.ETAN3).AND.(ETAN3.GT.ETAN1)) THEN
                  IF((HTOTN2.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
                ENDIF
                IF((ETAN3.GE.ETAN1).AND.(ETAN1.GT.ETAN2)) THEN
                   IF((HTOTN3.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
                ENDIF
                !...ACB pattern
                IF((ETAN1.GE.ETAN3).AND.(ETAN3.GT.ETAN2)) THEN
                  IF((HTOTN1.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
                ENDIF
                IF((ETAN2.GE.ETAN1).AND.(ETAN1.GT.ETAN3)) THEN
                  IF((HTOTN2.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
                ENDIF
                IF((ETAN3.GE.ETAN2).AND.(ETAN2.GT.ETAN1)) THEN
                  IF((HTOTN3.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
                ENDIF
            ENDIF  ! END NOFF ACTIVE FLAG 
        ENDIF ! END NO BARRIERS OR PIPES 
    ENDDO ! END ELEMENTAL LOOP 

!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
    subroutine dryingCriteriaD2
!----------------------------------------------------------------------
!   D2) Update number of active elements (MJU) and the total area (TotalArea) connected
!       to a node. If these are zero, the node is landlocked and should be dried.
!       These depend on NNODECODE which varies during the time step
!----------------------------------------------------------------------
        IMPLICIT NONE 

        INTEGER :: NM1,NM2,NM3 
        INTEGER :: NC1,NC2,NC3
        INTEGER :: NCEle
        REAL(SZ) :: AreaEle

        INTEGER I
        INTEGER IE 

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


         DO I=1,NP
            IF((NNODECODE(I).EQ.1).AND.(MJU(I).EQ.0)) THEN
               NNODECODE(I)=0
            ENDIF
            IF(MJU(I).EQ.0) THEN 
                MJU(I)=1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
            ENDIF 
         ENDDO

!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
    subroutine updateNnodecodeAtGhosts() 
!----------------------------------------------------------------------
!   Use Message-Passing to update nnodecode at ghost nodes
!   WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE
!   Check to see if any wetting occurred & update NODECODE
!   Note, NCCHANGE=0 set near the beginning of GWCE subroutine
!----------------------------------------------------------------------
    IMPLICIT NONE 

    INTEGER :: I 


    DO I=1,NP
       IF(NNODECODE(I).NE.NODECODE(I)) THEN
          NODECODE(I)=NNODECODE(I)
          NCCHANGE=NCCHANGE+1
       ENDIF
    ENDDO

!----------------------------------------------------------------------
    end subroutine 
!----------------------------------------------------------------------


!----------------------------------------------------------------------
    subroutine finalizeWettingAndDrying() 
!----------------------------------------------------------------------
!   WET/DRY SECTION - PART 6
!   Check to see if any NOFF changed requiring the matrix to be reset
!   Note, NCCHANGE=0 set near the beginning of GWCE subroutine
!----------------------------------------------------------------------
  IMPLICIT NONE 

  INTEGER :: I 

  DO I=1,NE
     IF(NOFF(I).NE.NOFFOLD(I)) NCCHANGE=NCCHANGE+1
  ENDDO

 !jgf45.06 If there has been any wetting or drying in any
 !of the subdomains, the NCCHANGE flag will be activated on all
 !of the subdomains, to prevent them from getting out of sync
 !with their MPI calls as some reset the GWCE and others do not.

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
!----------------------------------------------------------------------
    end subroutine
!----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END MODULE WetDry
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

