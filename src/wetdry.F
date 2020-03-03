C***********************************************************************
C                     M O D U L E   W E T   D R Y      
C***********************************************************************
C     Executes wetting and drying algorithm.
C***********************************************************************
C Logical Variable List (default value .FALSE., set in global.f)       *
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C     C2DDI            - 2D Depth Integrated model run                 *
C     C3D              - 3D model run                                  *
C     C3DDSS           - Stress form of the 3D momentum equations      *
C     C3DVS            - Velocity form of the 3D momentum equations    *
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
C   See header.f for a summary history of code modifications.          *
C***********************************************************************
CWET...
CWET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
CWET...
CWET...NOTE:NNODECODE is a working variable that can change more than once
CWET...               during a time step
CWET...     NNODECODE = 0 for a dry node
CWET...     NNODECODE = 1 for a wet node
CWET...     NODECODE  - is a more static version of NNODECODE that is reconciled
CWET...                 once and for all at the end time step
CWET...
CWET...
CWET...        (   DRYING CRITERIA   )
CWET...
CWET...A node should be dry under two conditions.
CWET...D1.) If the total water depth falls below H0.
CWET .......Note: if the total water depth falls below 0.8*H0, the surface elevation
CWET........is lifted up so that the total water depth = 0.8*H0.
CWET......
CWET...D2.) If the node is connected to only nonfunctioning (dry) elements.  In
CWET........this case the node is dried due to becoming landlocked.
CWET........Note: this criteria is applied after all other wetting and drying criteria
CWET...
CWET...An element should be dry under the following conditions.
CWET...DE3.) This is an elemental check section designed to avoid artificial wetting of
CWET.........of control sections
CWET.........All elements where downhill flow originates from a barely wet node
CWET.........(defined as 1.2*H0) into wet nodes are forced inactive; the only exception
CWET......... is receiving overtopped barrier nodes
CWET...
CWET...        (   WETTING CRITERIA   )
CWET...
CWET...A node should be wet under two conditions.
CWET...W1.) If 2 nodes on an element are wet and one is dry, wet the dry node
CWET........if the water level at one of the wet nodes is greater than the
CWET........water level at the dry node and the steady state velocity that
CWET........would result from a balance between the water level gradient and
CWET........bottom friction would yield a velocity > VELMIN.
CWET........Note that the criteria outlined in DE3 must also be satified before
CWET........the node is allowed to wet
CWET...
CWET...W2.) If an element has a node lying on a receiving internal barrier boundary or
CWET......specified discharge boundary that is actively discharging flow into the
CWET......domain at that node, all nodes in this element must stay wet.
CWET...
CWET...
CWET...        (  VELOCITY BOUNDARY CONDITION  )
CWET...
CWET...Either a natural or essential boundary condition can be used as a velocity
CWET...boundary condition in the momentum equation solution along a wet/dry boudary
CWET...To use a natural boundary condition, do nothing along the wet/dry interface.
CWET...To use an essential, no velocity boundary condition, identify the nodes along
CWET...the wet/dry interface and zero out the velocity at the nodes.  Interface nodes
CWET...can easily be identified by comparing the number of active elements a node is
CWET...connected to (MJU) to the total number of elements a node is connected to (NODELE).
CWET...If MJU < NODELE for any node, it must lie along the wet/dry interface.  See
CWET...further comments at the end of the momentum equation solution section.
CWET...
      module wetdry
      use sizes, only : sz
      
      real(sz) :: habsmin
!$omp threadprivate(habsmin)
      real(sz) :: hoff
!$omp threadprivate(hoff)
      real(sz) :: velmin
!$omp threadprivate(velmin)
      
      integer,allocatable ::    nibcnt(:)
!$omp threadprivate(nibcnt)
      integer,allocatable ::    noffold(:)
!$omp threadprivate(noffold)

      ! jgf52.08.08: Enable analyst to remove NOFF from wetting and
      ! drying algorithm; .true. by default to conform with prior
      ! ADCIRC versions. 
      LOGICAL :: noffActive = .true.
!$omp threadprivate(noffActive)

      contains

!----------------------------------------------------------------------
!                  S U B R O U T I N E
!     I N I T I A L I Z E   W E T T I N G   A N D   D R Y I N G
!----------------------------------------------------------------------
!     Allocates arrays and initializes constants.
!----------------------------------------------------------------------
      subroutine initializeWettingAndDrying()
      use sizes, only : mne, mnp
      use global, only : h0, nodecode, nnodecode
      implicit none

      allocate ( nibcnt(mnp) )
      allocate ( noffold(mne))
      nnodecode = 1
      nodecode = 1
      noffold(:) = 1
      habsmin=0.8d0*h0
      hoff=1.2d0*h0

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
      use sizes, only : sz, mne, mnp
      use global, only : noff, nodecode, nnodecode, eta2, tk, nolifa,
     &    bsx1, bsy1, btime_end, C2DDI, C3D, g, h0, ifnlfa, nddt,
     &    nibnodecode, ilump, ncchange, tkm
#ifdef CMPI
     &    ,idumy
      use messenger
#endif
      use mesh, only : ne, np, dp, mju, totalArea, nm, x, y, areas
      use nodalattributes, only : BFCdLLimit, fgamma, ftheta, fric,
     &   manningsn, hbreak, ifhybf, ifnlbf, iflinbf, loadManningsN,
     &   loadZ0B_var, z0b_var
      use global_3dvs, only : a, b, islip, kp, z0b, sigma, evtot, q
      use subdomain, only : subdomainOn, enforceBN, enforceWDcb, enforceWDob
      implicit none
      integer, intent(in) :: it ! time step number
      complex(sz) :: duds !jgf48.50 declare size SZ instead of plain COMPLEX
      integer :: nc1, nc2, nc3
      integer :: nm1, nm2, nm3
      integer :: nm123
      integer :: ncele
      integer :: nctot
      real(sz) :: kslip
      real(sz) :: vel
      real(sz) :: z0b1
      real(sz) :: areaEle
      real(sz) :: tkWet
      real(sz) :: etaN1,etaN2,etaN3
      real(sz) :: hTotN1,hTotN2,hTotN3      
      real(sz) :: deldist,deleta
      real(sz) :: htot
      real(sz) :: h1
      integer :: nbnctot
      integer :: i
      integer :: ie

      if (nolifa.ne.2) then
         return ! wetting and drying is not active 
      endif

CWET...
CWET...WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP
CWET...
         DO I=1,NP
            NIBCNT(I) = 0
         ENDDO
         DO I=1,NE
            NOFFOLD(I)=NOFF(I)
            NOFF(I)=1
         ENDDO
CWET...
CWET...WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
CWET....Drying Criteria D1: this depends on NODECODE and updates NODECODE
CWET...
         DO I=1,NP
            IF(NODECODE(I).EQ.1) THEN
               HTOT=DP(I)+ETA2(I)
               IF(HTOT.LE.H0) THEN
                  IF(HTOT.LT.HABSMIN) ETA2(I)=HABSMIN-DP(I)
                  NNODECODE(I)=0
                  NODECODE(I)=0
                  NCCHANGE=NCCHANGE+1 !NCCHANGE=0 set near beginning of GWCE
C                  ENDIF
               ENDIF
            ENDIF
         ENDDO
CWET...
CWET...END WET/DRY SECTION - PART 1
CWET...

cjjwC     Use Message-Passing to update nodecode and nnodecode at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(NODECODE,NNODECODE,2)
cjjw#endif

CWET...
CWET...WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
CWET...
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)

CWET...
CWET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
CWET...
            NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
            IF(NCTOT.EQ.2) THEN
               ETAN1=ETA2(NM1)
               ETAN2=ETA2(NM2)
               ETAN3=ETA2(NM3)
               HTOTN1=DP(NM1)+ETA2(NM1)
               HTOTN2=DP(NM2)+ETA2(NM2)
               HTOTN3=DP(NM3)+ETA2(NM3)
               IF((NODECODE(NM1).EQ.1).AND.(NODECODE(NM2).EQ.1)) THEN
                  IF((HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF)) THEN
                     NM123=NM1
                     IF(ETA2(NM2).GT.ETA2(NM1)) NM123=NM2
                     DELDIST=SQRT((y(NM3)-y(NM123))**2.D0
     &                    +(X(NM3)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM3)
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged from Casey 071219: Added the following logic for 3D friction.
c. RJW modified the following for 3D friction
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))

                        IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                        VEL=G*(DELETA/DELDIST)/TKWET

                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                       IF(LoadZ0B_var) THEN
                           Z0B1 = Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1 = Z0B
                       ENDIF
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                    ENDIF

                     IF(VEL.GT.VELMIN) THEN
C    ....         third node met criteria and is also wet
                        NNODECODE(NM3)=1
c. RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
c                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123) 
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM3)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              ENDIF
                              BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM3)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ELSEIF((NODECODE(NM2).EQ.1).AND.(NODECODE(NM3).EQ.1))
     &          THEN
                  IF((HTOTN2.GE.HOFF).AND.(HTOTN3.GE.HOFF)) THEN
                     NM123=NM2
                     IF(ETA2(NM3).GT.ETA2(NM2)) NM123=NM3
                     DELDIST=SQRT((Y(NM1)-Y(NM123))**2.D0
     &                    +(X(NM1)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM1)
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                        IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                        VEL=G*(DELETA/DELDIST)/TKWET

                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                        IF(LoadZ0B_var) THEN
                           Z0B1=Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1=Z0B
                        ENDIF
                        VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     ENDIF

                     IF(VEL.GT.VELMIN) THEN
                        NNODECODE(NM1)=1
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123) 
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM1)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                              ENDIF
                              BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                              BSX1(NM1)=KSLIP*REAL(Q(NM123,1))
                              BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                           ENDIF
                        ENDIF

                     ENDIF
                  ENDIF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ELSEIF((NODECODE(NM3).EQ.1).AND.(NODECODE(NM1).EQ.1))
     &          THEN
                  IF((HTOTN3.GE.HOFF).AND.(HTOTN1.GE.HOFF)) THEN
                     NM123=NM3
                     IF(ETA2(NM1).GT.ETA2(NM3)) NM123=NM1
                     DELDIST=SQRT((Y(NM2)-Y(NM123))**2.D0
     &                     +(X(NM2)-X(NM123))**2.D0)
                     DELETA=ETA2(NM123)-ETA2(NM2)
                     ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                     IF (DELETA.lt.0.d0) DELETA = 0.d0
                     H1=ETA2(NM123)+DP(NM123)
c. RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                     IF(C2DDI)THEN
Csb46.28sb02
C<<                     Convert Manning's N to Cd, if necessary.
                        IF (LoadManningsN) THEN
                          FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                      /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                          IF(FRIC(NM123).LT.BFCdLLimit) THEN
                            FRIC(NM123) = BFCdLLimit
                          ENDIF
                        ENDIF
C>>
                        TKWET=FRIC(NM123)*(IFLINBF+(VELMIN/H1)*
     &                    (IFNLBF+IFHYBF*
     &                    (1.D0+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                       IF(TKWET.LT.0.0001d0) TKWET=0.0001d0
                       VEL=G*(DELETA/DELDIST)/TKWET
                     ELSEIF(C3D)THEN
C solve for the depth averaged velocity,U, from the relation :
C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                        IF(LoadZ0B_var) THEN
                           Z0B1 = Z0B_var(NM123)
                        ELSEIF (LoadManningsN) THEN
                     Z0B1 = ( DP(NM123)+IFNLFA*ETA2(NM123) )* exp(-(1.0D0+
     &               ( (0.41D0*( DP(NM123)+IFNLFA*ETA2(NM123))**(1.0D0/6.0D0) )/
     &                                 (ManningsN(NM123)*sqrt(g)) ) ))
                        ELSE
                           Z0B1 = Z0B
                        ENDIF
                      VEL=sqrt(g*H1*(DELETA/DELDIST))
     &                 * ((H1+Z0B1)*LOG((H1+Z0B1)/Z0B1)-H1)/(H1*0.41D0)
                     ENDIF

                     IF(VEL.GT.VELMIN) THEN
                        NNODECODE(NM2)=1
c. RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                        IF(C2DDI)THEN
                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                        (IFNLBF+IFHYBF*
     &                        (1.D0+(HBREAK/H1)**FTHETA)**
     &                        (FGAMMA/FTHETA)))
                           !WJP needed for internal_tide stuff
                           TKM(1:2,NM123) = TK(NM123) 
                       ELSEIF(C3D)THEN
                           IF(ISLIP.EQ.0)THEN
                              DUDS=(Q(NM123,2)-Q(NM123,1))
     &                           /(SIGMA(2)-SIGMA(1))
                              BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                              BSX1(NM2)=EVTOT(1)*REAL(DUDS)
                              BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                           ENDIF
                           IF(ISLIP.NE.0)THEN
                              IF(ISLIP.EQ.1)THEN
                                 KSLIP=KP
                              ENDIF
                              IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                            LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
        IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
        IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
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
CWET...
CWET...Nodal Wetting Criteria W2a
CWET...
            NBNCTOT=NIBNODECODE(NM1)+NIBNODECODE(NM2)+NIBNODECODE(NM3)
            NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
            NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
            NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

         ENDDO

         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
         CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif


Cwet...
CWET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
C*******************************************************************************************
C tcm v50.66.01 -- This is an additional test for wetting only when time varying 
C                  bathymetry is used and is only performed during the period of 
C                  bathymetry evolution.
C
         IF ((NDDT.NE.0).AND.(IT.LE.BTIME_END+1) ) THEN
            DO I=1,NE
               NM1=NM(I,1)
               NM2=NM(I,2)
               NM3=NM(I,3)
               NCTOT=NODECODE(NM1)+NODECODE(NM2)+NODECODE(NM3)
               IF(NCTOT.lt.3) THEN   !If not wet from previous time step
                  NCTOT=NNODECODE(NM1)+NNODECODE(NM2)+NNODECODE(NM3)
                  if(NCTOT.lt.3) then !if not alreay made wet for this time step
                     ETAN1=ETA2(NM1)
                     ETAN2=ETA2(NM2)
                     ETAN3=ETA2(NM3)
                     HTOTN1=DP(NM1)+ETA2(NM1)
                     HTOTN2=DP(NM2)+ETA2(NM2)
                     HTOTN3=DP(NM3)+ETA2(NM3)

!                    if all nodes have a depth greater than or equal to 
!                    hoff = 1.2*H0, then make the element wet
                     IF( (HTOTN1.GE.HOFF).AND.(HTOTN2.GE.HOFF).AND.
     &                   (HTOTN3.GE.HOFF) ) THEN
                     !THE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                        !  Make Node 1 Wet and set parameters
                        IF(NNODECODE(NM1).NE.1) THEN  !node 1
                           NNODECODE(NM1)=1
                           NM123 = NM1
                           IF(C2DDI)THEN
C<<                           Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                              /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL=VELMIN
                           H1 = HTOTN1
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                          (IFNLBF+IFHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           ELSEIF(C3D)THEN
                              IF(ISLIP.EQ.0)THEN
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                              /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM3)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM3)=EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF(ISLIP.NE.0)THEN
                                 IF(ISLIP.EQ.1)THEN
                                    KSLIP=KP
                                 ENDIF
                                 IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                     IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                     IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                                 ENDIF
                                 BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                 BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                 BSX1(NM3)=KSLIP*REAL(Q(NM123,1))
                                 BSY1(NM3)=KSLIP*AIMAG(Q(NM123,1))
                              ENDIF
                           ENDIF
                        ENDIF  !end node 1

                        !  Make Node 2 Wet and set parameters
                        IF (NNODECODE(NM2).NE.1) THEN
                           NNODECODE(NM2) = 1
                           NM123=NM2
                           IF(C2DDI)THEN
C<<                        Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                            /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL = VELMIN
                           H1 = HTOTN2
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                          (IFNLBF+IFHYBF*
     &                          (1.D0+(HBREAK/H1)**FTHETA)**
     &                          (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           ELSEIF(C3D)THEN
                              IF(ISLIP.EQ.0)THEN
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                              /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM1)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM1)=EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF(ISLIP.NE.0)THEN
                                 IF(ISLIP.EQ.1)THEN
                                    KSLIP=KP
                                 ENDIF
                                 IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                 (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                    IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                    IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
                                 ENDIF
                                 BSX1(NM123)=KSLIP*REAL(Q(NM123,1))
                                 BSY1(NM123)=KSLIP*AIMAG(Q(NM123,1))
                                 BSX1(NM1)=KSLIP*REAL(Q(NM123,1))
                                 BSY1(NM1)=KSLIP*AIMAG(Q(NM123,1))
                              ENDIF
                           ENDIF
                        ENDIF  !node 2

                        !  Make Node 3 Wet and set parameters
                        IF(NNODECODE(NM3).NE.1) THEN
                           NNODECODE(NM3)=1
                           NM123 = NM3
                           IF(C2DDI)THEN
C<<                         Convert Manning's N to Cd, if necessary.
                              IF (LoadManningsN) THEN
                                 FRIC(NM123)=g*ManningsN(NM123)**2.d0
     &                           /( ( DP(NM123)+IFNLFA*ETA2(NM123) )
     &                                               **(1.d0/3.d0) )
                                 IF(FRIC(NM123).LT.BFCdLLimit) THEN
                                    FRIC(NM123) = BFCdLLimit
                                 ENDIF
                              ENDIF
                           ENDIF
                           VEL=VELMIN
                           H1 = HTOTN3
                           IF(C2DDI)THEN
                              TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
     &                           (IFNLBF+IFHYBF*
     &                           (1.D0+(HBREAK/H1)**FTHETA)**
     &                           (FGAMMA/FTHETA)))
                              !WJP needed for internal_tide stuff
                              TKM(1:2,NM123) = TK(NM123) 
                           ELSEIF(C3D)THEN
                              IF(ISLIP.EQ.0)THEN
                                 DUDS=(Q(NM123,2)-Q(NM123,1))
     &                             /(SIGMA(2)-SIGMA(1))
                                 BSX1(NM123)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123)=EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM2)=EVTOT(1)*REAL(DUDS)
                                 BSY1(NM2)=EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF(ISLIP.NE.0)THEN
                                 IF(ISLIP.EQ.1)THEN
                                    KSLIP=KP
                                 ENDIF
                                 IF(ISLIP.EQ.2)THEN
                                     KSLIP = (1.D0 /
     &                                     ( (1.D0/0.41D0) *
     &                               LOG((ABS(((SIGMA(2)-SIGMA(1))/(A-B))*
     &                                  (DP(NM123)+IFNLFA*ETA2(NM123)))
     &                                       +Z0B1)
     &                                       / (Z0B1) ) ) )**2.D0
     &                                     * ABS(Q(NM123,1))
                                     IF(KSLIP.GT.1.D0*ABS(Q(NM123,1)))
     &                                  KSLIP=1.D0* ABS(Q(NM123,1))
                                     IF(KSLIP.LT.0.0025D0*ABS(Q(NM123,1)))
     &                                  KSLIP=0.0025D0*ABS(Q(NM123,1))
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

            if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain

C     Use Message-Passing to update nnodecode and nibcnt at ghost nodes
#ifdef CMPI
            CALL UPDATEI(NNODECODE,NIBCNT,2)
#endif
C
         ENDIF !IT TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME
C
C.... END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY
CWET..
CWET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
C*******************************************************************************************

CWET...
CWET...Nodal Wetting Criteria W2b
CWET...Check for adjacent nodes and force nodes wet when attached
CWET...to receiving barrier nodes
CWET...
         DO I=1,NP
            IF((NIBCNT(I).GT.0).AND.(NNODECODE(I).EQ.0)) THEN
               NNODECODE(I)=1
            ENDIF
         ENDDO

cjjwC     Use Message-Passing to update nnodecode at ghost nodes
cjjw#ifdef CMPI
cjjw         CALL UPDATEI(NNODECODE,IDUMY,1)
cjjw#endif

CWET...
CWET...END WET/DRY SECTION - PART 2
CWET...

CWET...
CWET...START WET/DRY SECTION  - PART 3
CWET...Elemental drying criteria DE1
CWET...This is an elemental check section designed to avoid artificial wetting of
CWET....of control sections
CWET...All elements where downhill flow originates from a barely wet node
CWET....into wet nodes are forced inactive; the only exception is receiving
CWET....overtopped barrier nodes
CWET...
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
#ifdef SB_WETDRY
C...Find the heighest point on the bed in the element.  sb v46.28.sb05.06 11/01/2006
               IF(DP(NM1).LE.DP(NM2).AND.DP(NM1).LE.DP(NM3)) THEN
                  DPMIN = DP(NM1)
               ELSE IF(DP(NM2).LE.DP(NM3).AND.DP(NM2).LE.DP(NM1)) THEN
                  DPMIN = DP(NM2)
               ELSE IF(DP(NM3).LE.DP(NM1).AND.DP(NM3).LE.DP(NM2)) THEN
                  DPMIN = DP(NM3)
               ENDIF
#endif
#ifndef SB_WETDRY
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
C...ABC pattern
                  IF((ETAN1.GE.ETAN2).AND.(ETAN2.GT.ETAN3)) THEN
                    IF((HTOTN1.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
                  ENDIF
                  IF((ETAN2.GE.ETAN3).AND.(ETAN3.GT.ETAN1)) THEN
                    IF((HTOTN2.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
                  ENDIF
                  IF((ETAN3.GE.ETAN1).AND.(ETAN1.GT.ETAN2)) THEN
                     IF((HTOTN3.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
                  ENDIF
C...ACB pattern
                  IF((ETAN1.GE.ETAN3).AND.(ETAN3.GT.ETAN2)) THEN
                    IF((HTOTN1.LT.HOFF).OR.(HTOTN3.LT.HOFF)) NOFF(I)=0
                  ENDIF
                  IF((ETAN2.GE.ETAN1).AND.(ETAN1.GT.ETAN3)) THEN
                    IF((HTOTN2.LT.HOFF).OR.(HTOTN1.LT.HOFF)) NOFF(I)=0
                  ENDIF
                  IF((ETAN3.GE.ETAN2).AND.(ETAN2.GT.ETAN1)) THEN
                    IF((HTOTN3.LT.HOFF).OR.(HTOTN2.LT.HOFF)) NOFF(I)=0
                  ENDIF
               endif
#else
               ! jgf52.08.08: Analyst can eliminate noff from
               ! consideration in fort.15 file. 
               if (noffActive.eqv..true.) then
C...ABC pattern
                  IF((ETAN1.GE.ETAN2).AND.(ETAN2.GT.ETAN3)) THEN
                     IF((HTOTN1.LT.HOFF)) NOFF(I)=0
                     IF((ETAN1-ETAN2).LT.(ETAN2-ETAN3)) THEN
                        IF(HTOTN2.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
                  IF((ETAN2.GE.ETAN3).AND.(ETAN3.GT.ETAN1)) THEN
                     IF((HTOTN2.LT.HOFF)) NOFF(I)=0
                     IF((ETAN2-ETAN3).LT.(ETAN3-ETAN1)) THEN
                        IF(HTOTN3.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
                  IF((ETAN3.GE.ETAN1).AND.(ETAN1.GT.ETAN2)) THEN
                     IF((HTOTN3.LT.HOFF)) NOFF(I)=0
                     IF((ETAN3-ETAN1).LT.(ETAN1-ETAN2)) THEN
                        IF(HTOTN1.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
C...ACB pattern
                  IF((ETAN1.GE.ETAN3).AND.(ETAN3.GT.ETAN2)) THEN
                     IF((HTOTN1.LT.HOFF)) NOFF(I)=0
                     IF((ETAN1-ETAN3).LT.(ETAN3-ETAN2)) THEN
                        IF(HTOTN3.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
                  IF((ETAN2.GE.ETAN1).AND.(ETAN1.GT.ETAN3)) THEN
                     IF((HTOTN2.LT.HOFF)) NOFF(I)=0
                     IF((ETAN2-ETAN1).LT.(ETAN1-ETAN3)) THEN
                        IF(HTOTN1.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
                  IF((ETAN3.GE.ETAN2).AND.(ETAN2.GT.ETAN1)) THEN
                     IF((HTOTN3.LT.HOFF)) NOFF(I)=0
                     IF((ETAN3-ETAN2).LT.(ETAN2-ETAN1)) THEN
                        IF(HTOTN2.LT.HOFF) NOFF(I)=0
                     ENDIF
                  ENDIF
               endif !noffActive.eqv..true.
#endif

#ifdef SB_WETDRY
C...An element is set to be dry if it is determined to be a flooding type
C...wetting element.  An element is a flooding type wetting element if
C...the bed elevation at the node with the biggest water column height
C...is lower than the heighest point on the bed in the element.
C...sb v46.28.sb05.06 11/01/2006
C...This is applied only when NOFF flag of the element at the previous time step
C...is 0, which means that this logic works to prevent an element from
C...re-wetting.
C...sb v46.52.03
               if (noffActive.eqv..true.) then 
                  IF(NOFFOLD(I).EQ.0) THEN
                     IF(HTOTN1.GE.HTOTN2.AND.HTOTN1.GE.HTOTN3) THEN
                        IF(ETAN1.LT.(-DPMIN+H0)) NOFF(I) = 0
                     ENDIF
                     IF(HTOTN2.GE.HTOTN3.AND.HTOTN2.GE.HTOTN1) THEN
                        IF(ETAN2.LT.(-DPMIN+H0)) NOFF(I) = 0
                     ENDIF
                     IF(HTOTN3.GE.HTOTN1.AND.HTOTN3.GE.HTOTN2) THEN
                        IF(ETAN3.LT.(-DPMIN+H0)) NOFF(I) = 0
                     ENDIF
                 ENDIF
              endif
#endif
            ENDIF
         ENDDO

#ifdef SB_WETDRY
CWET......added by sb on 11/02/2006
CWET...
CWET...This section is added after we realize that it's not possible
CWET...to compute a correct flow going through two elements if
CWET...the elements are connected just by one node. i.e., elements need
CWET...to share an edge to let the flow go through between the elements.
CWET...Therefore, in this section, a node is determined to be dry
CWET...if two elements are connected at one node, not sharing an edge.
CWET...As it seemed this procedure needed NOFF information,
CWET...although this section changes NNODECODE,
CWET...I put this section here, rather than the end of PART 2.
CWET...I tried setting NOFF(I) = 0, but it didn't shut down the flow.
CWET...
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            IF(NOFF(I).EQ.1.AND.
     &         NNODECODE(NM1).EQ.1.AND.
     &         NNODECODE(NM2).EQ.1.AND.
     &         NNODECODE(NM3).EQ.1) THEN
               DO K=1,3
                  NM1=NM(I,K)
                  NM2=NM(I,MOD(K+0,3)+1)
                  NM3=NM(I,MOD(K+1,3)+1)

                  NWETNEI = 0
                  NWETADJ = 0
                  DO J=1,MNEI
                     N=NeiTabEle(NM1,J)
                     IF(N.EQ.0) CYCLE
                     IF(N.EQ.I) CYCLE

                     NMN1=NM(N,1)
                     NMN2=NM(N,2)
                     NMN3=NM(N,3)
                     IF(NOFF(N).EQ.1.AND.
     &                  NNODECODE(NMN1).EQ.1.AND.
     &                  NNODECODE(NMN2).EQ.1.AND.
     &                  NNODECODE(NMN3).EQ.1) THEN
                        NWETNEI = NWETNEI + 1

                        IF((NMN1.EQ.NM2.OR.NMN1.EQ.NM3).OR.
     &                     (NMN2.EQ.NM2.OR.NMN2.EQ.NM3).OR.
     &                     (NMN3.EQ.NM2.OR.NMN3.EQ.NM3)) THEN
                           NWETADJ = NWETADJ + 1
                        ENDIF
                     ENDIF
                  ENDDO

                  IF(NWETNEI.GT.0.AND.NWETADJ.EQ.0.AND.
     &               NIBCNT(NM1).EQ.0) THEN
                     NNODECODE(NM1) = 0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
#endif

CWET...
CWET...END WET/DRY SECTION  - PART 3
CWET...

CWET...
CWET...START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2
CWET...Update number of active elements (MJU) and the total area (TotalArea) connected
CWET...to a node. If these are zero, the node is landlocked and should be dried.
CWET...These depend on NNODECODE which varies during the time step
CWET...
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

cjjwnote - looks like this is used later in momentum equations
cjjwnote - this has implications on making this into a subroutine

         DO I=1,NP
            IF((NNODECODE(I).EQ.1).AND.(MJU(I).EQ.0)) THEN
               NNODECODE(I)=0
            ENDIF
            IF(MJU(I).EQ.0) MJU(I)=1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
         ENDDO

C     WET...
C     WET...END WET/DRY SECTION - PART 4
C     WET...

cjjwnote - may have to pass TotalArea and mju as well

         if (subdomainOn.and.enforceBN.eq.1) call enforceWDcb() ! NCSU Subdomain
         if (subdomainOn.and.enforceBN.eq.2) call enforceWDob() ! NCSU Subdomain


C     Use Message-Passing to update nnodecode at ghost nodes
#ifdef CMPI
         CALL UPDATEI(NNODECODE,IDUMY,1)
#endif
CWET...
CWET...WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE
CWET...Check to see if any wetting occurred & update NODECODE
CWET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
CWET...
         DO I=1,NP
            IF(NNODECODE(I).NE.NODECODE(I)) THEN
               NODECODE(I)=NNODECODE(I)
               NCCHANGE=NCCHANGE+1
            ENDIF
         ENDDO
CWET...
CWET...END WET/DRY SECTION - PART 5
CWET...

CWET...
CWET...WET/DRY SECTION - PART 6
CWET...Check to see if any NOFF changed requiring the matrix to be reset
CWET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
CWET...
         DO I=1,NE
            IF(NOFF(I).NE.NOFFOLD(I)) NCCHANGE=NCCHANGE+1
         ENDDO
CWET...
CWET... jgf45.06 If there has been any wetting or drying in any
CWET... of the subdomains, the NCCHANGE flag will be activated on all
CWET... of the subdomains, to prevent them from getting out of sync
CWET... with their MPI calls as some reset the GWCE and others do not.
CWET...
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
CWET...
CWET...END WET/DRY SECTION - PART 6
CWET...

c.....
!-----------------------------------------------------------------------
      end subroutine computeWettingAndDrying
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end module wetDry
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

