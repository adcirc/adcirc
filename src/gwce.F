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
      
      contains


C******************************************************************************
C                                                                             *
C    Subroutine to compute the elevation using the GWCE formluation           *
C    Re-written to conform to the ADCIRC Theory Report                        *
C                                                                             *
C                            r.l.  06/22/2005                                 *
C******************************************************************************
      SUBROUTINE GWCE_New(IT,TimeLoc,TimeH)
C
#ifdef IEEE_DEBUG
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      USE SIZES
      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM, X, Y, DP, NNeigh, NeiTab, TotalArea,
     &                  Areas, NEIMAX, SFAC, nneighele, neitabele
      USE BOUNDARIES, ONLY : NETA, NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, 
     &   NFLUXRBC, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3 
      USE WIND
      USE ITPACKV
      USE ADCIRC_MOD, only : adcirc_terminate
      USE NodalAttributes, ONLY :
     &     LoadGeoidOffset, GeoidOffset, EVM,
     &     TAU0VAR, HighResTimeVaryingTau0, FullDomainTimeVaryingTau0,
     &     CalculateTimeVaryingTau0,
     &     LoadEleSlopeLim,elemental_slope_limiter_active,
     &     LoadAdvectionState,
     &     elemental_slope_limiter_grad_max,
     &     elemental_slope_limiter_max_exceeded

#ifdef CMPI
      USE MESSENGER
#endif
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      USE Couple2Swan, ONLY: TKXX,
     &                       TKXY,
     &                       TKYY
#endif
#endif
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceECB, 
     &  checkChange, enforceEOB, enforceEIB, enforceGWCELVOB

      IMPLICIT NONE

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      LOGICAL  DIE
      REAL(SZ) A00pB00
      REAL(SZ) BCXAvg, BCYAvg
      REAL(SZ) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(SZ) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(SZ) CorifAvg
      REAL(SZ) DPAvg, GDPAvgOAreaIE4
      REAL(SZ) DispX, DispY, DispXAvg, DispYAvg
      REAL(SZ) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(SZ) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(SZ) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(SZ) ESN1, ESN2, ESN3, ESAvg
      REAL(SZ) EVMH, EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(SZ) EVMEle, EVMSmag
      REAL(SZ) GA00DPAvgOAreaIE4
      REAL(SZ) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(SZ) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg, H1, H2
      REAL(SZ) H2OTotalArea
      REAL(SZ) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(SZ) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(SZ) MsFacR, MsFacLOnDiag, MsFacLOffDiag
      REAL(SZ) MX, MY, MXAvg, MYAvg
      REAL(SZ) JXAvg, JYAvg
      REAL(SZ) Pr1N1, Pr1N2, Pr1N3
      REAL(SZ) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(SZ) SFacAvg
      REAL(SZ) T0N1,T0N2, T0N3
      REAL(SZ) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(SZ) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TiPN1, TiPN2, TiPN3
      REAL(SZ) U1N1,U1N2,U1N3, U1Avg
      REAL(SZ) V1N1,V1N2,V1N3, V1Avg
      REAL(SZ) WSXAvg, WSYAvg
      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) TimeLoc, TimeH
C
C     jgf50.44: Added water surface elevation smoothing
      REAL(sz) :: EtaN1, ETaN2, EtaN3, AreaEle, EtaN123
      INTEGER :: IT_SmoothTime
!      moved to subroutine -- tcm
!      REAL(sz), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
!      LOGICAL, SAVE :: firstCall = .true.
C
      REAL(SZ) HH1 !jgf46.02 Added for Katrina.
      REAL(SZ) RDIAG ! jgf48.4619 Seizo parameter for fully explicit mode
      REAL(SZ) dEta2Mag ! magnitude of the slope of the water surface elevation
      REAL(SZ) DEta2DX,DEta2DY
      REAL(SZ) DRhoDX,DRhoDY

      call setMessageSource("gwce_new")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
!   tcm moved to subroutine
!      IF (LoadEleSlopeLim.eqv..true.) THEN
!         IF (firstCall.eqv..true.) THEN
!            allocate(elevSum(mnp))
!            firstCall = .false.
!         ENDIF
!      ENDIF

      if (subdomainOn.and.enforceBN.eq.2) call checkChange()    ! NCSU Subdomain

Casey 050711 : Added for averaged variable Tau0.
C      REAL(SZ) :: CaseySum
C      REAL(SZ), ALLOCATABLE :: TAU0VARTEMP(:)

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
      IF(NCChange.GT.0) THEN !if any subdomain grid has changed
         NCChange=0
c        IF(NScreen.GT.0.AND.MYProc.EQ.0) WRITE(screenunit,3806)
c        WRITE(16,3806)
 3806    FORMAT(/,1X,'RE-SETTING GWCE SYSTEM MATRIX',/)

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
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
               TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         IF(C3D.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091022: Adopt Ethan's/Joannes's modified friction.
            DO I=1,NP
               TK(I)=0.25D0*(TKXX(I)+2.D0*TKXY(I)+TKYY(I))
            ENDDO
#endif
#endif
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         !jgf48.4619: Added Seizo's handling of Lumping vs non-Lumping
         IF ( ILump.eq.0 ) THEN ! default case: fully consistent LHS
            DO IE=1,NE
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg=(SFac(NMI1)+SFac(NMI2)+SFac(NMI3))/3.d0

               FDX1 = (Y(NMI2)-Y(NMI3))*SFacAvg         !b1 = 2*Area*dphi1/dx
               FDX2 = (Y(NMI3)-Y(NMI1))*SFacAvg         !b2 = 2*Area*dphi2/dx
               FDX3 = (Y(NMI1)-Y(NMI2))*SFacAvg         !b3 = 2*Area*dphi3/dx
               FDY1 = X(NMI3)-X(NMI2)                   !a1 = 2*Area*dphi1/dy
               FDY2 = X(NMI1)-X(NMI3)                   !a2 = 2*Area*dphi2/dy
               FDY3 = X(NMI2)-X(NMI1)                   !a3 = 2*Area*dphi3/dy

               AreaIE2=Areas(IE)
               AreaIE =AreaIE2/2.0d0
               AreaIE4=AreaIE2*2.0d0

               DPAvg=(DP(NMI1)+DP(NMI2)+DP(NMI3))/3.d0
               GA00DPAvgOAreaIE4=G*A00*DPAvg/AreaIE4
               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0
               MsFacLOnDiag =OnDiag *AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0
               MsFacLOffDiag=OffDiag*AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0

               DO JN=2,NEIMAX
                  IF(NeiTab(NMI1,JN).EQ.NMJ2) J12=JN
                  IF(NeiTab(NMI1,JN).EQ.NMJ3) J13=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ1) J21=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ3) J23=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ1) J31=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ2) J32=JN
               END DO

               Coef(NMI1,1)  =Coef(NMI1,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX1+FDY1*FDY1))*NCELE
               Coef(NMI1,J12)=Coef(NMI1,J12) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX2+FDY1*FDY2))*NCELE
               Coef(NMI1,J13)=Coef(NMI1,J13) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX1*FDX3+FDY1*FDY3))*NCELE
               Coef(NMI2,J21)=Coef(NMI2,J21) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX1+FDY2*FDY1))*NCELE
               Coef(NMI2,1)  =Coef(NMI2,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX2+FDY2*FDY2))*NCELE
               Coef(NMI2,J23)=Coef(NMI2,J23) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX2*FDX3+FDY2*FDY3))*NCELE
               Coef(NMI3,J31)=Coef(NMI3,J31) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX1+FDY3*FDY1))*NCELE
               Coef(NMI3,J32)=Coef(NMI3,J32) + (MsFacLOffDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX2+FDY3*FDY2))*NCELE
               Coef(NMI3,1)  =Coef(NMI3,1)   + (MsFacLOnDiag
     &                   +GA00DPAvgOAreaIE4*(FDX3*FDX3+FDY3*FDY3))*NCELE

            ENDDO
         ELSE
            !jgf48.4619: Add Seizo's construction of Lumped LHS matrix
            DO IE=1,NE ! Make LHS Lumped Matrix
               NMI1=NM(IE,1)
               NMI2=NM(IE,2)
               NMI3=NM(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg=(SFac(NMI1)+SFac(NMI2)+SFac(NMI3))/3.d0

               AreaIE2=Areas(IE)
               AreaIE =AreaIE2/2.0d0

               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.d0
               MsFacLOnDiag =OnDiag *AreaIE*(1.d0/DT+Tau0Avg/2.d0)/DT/12.d0

               Coefd(NMI1)  =Coefd(NMI1)   + (MsFacLOnDiag)*NCELE
               Coefd(NMI2)  =Coefd(NMI2)   + (MsFacLOnDiag)*NCELE
               Coefd(NMI3)  =Coefd(NMI3)   + (MsFacLOnDiag)*NCELE
            ENDDO
         ENDIF

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
               CALL ADCIRC_Terminate()
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
            IF (DIE) CALL ADCIRC_Terminate()
         ENDIF
      ENDIF                     !End of GWCE matrix setup
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
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
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

C...  Assemble the GWCE RHS except for the boundary integral terms

      DO 1037 IE=1,NE

C...     Set nodal values for each element

Corbitt 120322: Localized Advection
         IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
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
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
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
            BSXN1=TK(NM1)*QX1N1
            BSYN1=TK(NM1)*QY1N1
            BSXN2=TK(NM2)*QX1N2
            BSYN2=TK(NM2)*QY1N2
            BSXN3=TK(NM3)*QX1N3
            BSYN3=TK(NM3)*QY1N3
#ifdef CSWAN
#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
            BSXN1 = TKXX(NM1)*QX1N1 + TKXY(NM1)*QY1N1
            BSYN1 = TKXY(NM1)*QX1N1 + TKYY(NM1)*QY1N1
            BSXN2 = TKXX(NM2)*QX1N2 + TKXY(NM2)*QY1N2
            BSYN2 = TKXY(NM2)*QX1N2 + TKYY(NM2)*QY1N2
            BSXN3 = TKXX(NM3)*QX1N3 + TKXY(NM3)*QY1N3
            BSYN3 = TKXY(NM3)*QX1N3 + TKYY(NM3)*QY1N3
#endif
#endif
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

         IF((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartV) .OR.
     &        (CGWCE_LS_2PartSQ).OR.(CGWCE_LS_2PartSV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
         ENDIF

C...     Compute the spatial gradients of the velocity dispersion terms if 3D

         IF (C3D) THEN                         !3D bottom friction
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
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
         DPAvg=(DP(NM1)+DP(NM2)+DP(NM3))/3.d0
         GDPAvgOAreaIE4=G*DPAvg/AreaIE4
         HAvg=(H1N1+H1N2+H1N3)/3.d0
         GHAvg=G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.d0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.d0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
         ENDIF

         IF (C3D) THEN                !3D velocity dispersion
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
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.d0
         A00pB00=A00+B00

C...     Compute the JX, JY terms less the advection terms

         JXAvg = CorifAvg*QY1Avg
     &          -IFNLFA*GOAreaIE4*(E1N1SQ*FDX1+E1N2SQ*FDX2
     &                                        +E1N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1
     &                     +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3)
     &          +WSXAvg-BSXAvg+MXAvg-DispXAvg-BCXAvg+Tau0QXAvg

         JYAvg =-CorifAvg*QX1Avg
     &          -IFNLFA*GOAreaIE4*(E1N1SQ*FDY1+E1N2SQ*FDY2
     &                                        +E1N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1
     &                    +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3)
     &          +WSYAvg-BSYAvg+MYAvg-DispYAvg-BCYAvg+Tau0QYAvg

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
     &        +Tau0SpaVar

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
     &        +Tau0SpaVar


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
     &        +Tau0SpaVar

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

 1037 CONTINUE                  !End of elemental loop

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

      DO I=1,NP
         IF(CPRECOR) THEN
            ETAS0(I)=ETAS(I)
            Eta0(I)=Eta1(I)
         END IF
         Eta1(I)=Eta2(I)
         Eta2(I)=0.0d0
      END DO

C...  At elevation boundary condition nodes, determine the elevation at
C...  the s+1 time step
C...
C...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0.
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

C...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

      if (subdomainOn) then                     ! NCSU Subdomain
          if(enforceBN.eq.1) call enforceEcb()  ! NCSU Subdomain
      else                                      ! NCSU Subdomain

         IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
            IF(TimeLoc.GT.ETIME2) THEN
               ETIME1=ETIME2
               ETIME2=ETIME1+ETIMINC
               DO J=1,NETA
                  ESBIN1(J)=ESBIN2(J)
                  READ(19,*) ESBIN2(J)
               END DO
            ENDIF
            ETRATIO=(TimeLoc-ETIME1)/ETIMINC
            DO I=1,NETA
               NBDI=NBD(I)
               Eta2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
            END DO
         ENDIF
      endif                                   ! NCSU Subdomain
      ! 
      ! jgf46.02 Added the ability to include geoid offset on the boundary.
      !
      ! aaltuntas51.48: Deactivated geoidOffset for subdomain boundary 
      ! conditions.
      if (subdomainOn.and.enforceBN.eq.1) then  ! NCSU Subdomain
         continue                                ! NCSU Subdomain
      else 
         IF (LoadGeoidOffset) THEN
            DO I=1,NETA
               ETA2(NBD(I))=ETA2(NBD(I))+GeoidOffset(NBD(I))
            END DO
         ENDIF
      endif                                      ! NCSU Subdomain
C
C
C     jgf48.04 Added an inverted barometer boundary condition so that
C     low pressure systems can cross the boundary without creating an
C     elevation anomaly.
C     jgf52.30.04: Included a parameter from fort.15 to turn this 
C     on and off according to analyst preference. 
      if (invertedBarometerOnElevationBoundary.eqv..true.) then
         DO I=1,NETA
            ETA2(NBD(I))=ETA2(NBD(I))
     &          + RampMete*(101300.d0/(RHOWAT0*G)) - PR2(NBD(I))
         END DO
      endif

C   kmd48.33bc add information for the levels of no motion boundary conditions
C              these are considered the steric adjustments.
      IF ((ABS(RES_BC_FLAG).GE.1).AND.(CBaroclinic).AND.(NOPE.GT.0))THEN
         DO I=1,NETA
           NBDI=NBD(I)
           ETA2(NBDI) = ETA2(NBDI) + LNM_BC(I)
         END DO
      END IF

C...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
C...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

C...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
C...  assume that QN is positive into the domain.  QFORCEJ has a -1
C...  built in and the terms are not explicitly negated. Boundary
C...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
C...  flux that  is positive out of the domain.  Therefore, to match
C...  the formulation these terms must be explicitly multiplied by -1.

C...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 +
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

C     jgf46.21 Added IBTYPE=52.
         IF(LBCODEI(1).EQ.52) THEN
            QFORCEJ=(QN2(1)-QN0(1))/DT2 + Tau0VAR(NBDJ)*QN1(1)
            IF (IT.GT.FluxSettlingIT) THEN
               HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
               Celerity=SQRT(G*HH1)
               QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &              + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(1)))
            ENDIF
         ENDIF
C
         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)
            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0


            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)
C
C     jgf46.21 Added IBTYPE=52
            IF(LBCODEI(J).EQ.52) THEN
               QFORCEJ=(QN2(J)-QN0(J))/DT2 + Tau0VAR(NBDJ)*QN1(J)
               IF (IT.GT.FluxSettlingIT) THEN
                  HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
                  Celerity=SQRT(G*HH1)
                  QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &                 + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(J)))
               ENDIF
            ENDIF

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.d0*QForceJ+QForceI)
          ENDDO
        ENDIF

C...
C...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
C...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
C...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
C...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
C...  EQUATIONS
C...
      IF ( ILump.eq.0 ) THEN ! default, fully consistent GWCE LHS
         DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
            DO J=2,NNEIGH(NBDI)
               GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &              -ETAS(NBDI)*OBCCOEF(I,J-1)
            END DO
         END DO
      ELSE                   ! ILump.eq.1, lumped GWCE
         DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*COEFD(NBDI)
         END DO
      ENDIF
C...
C...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
C...

C...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
C...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
         CALL UPDATER(GWCE_LV,COEF(1,1),DUMY1,2)
      ELSE  ! lumped LHS
         CALL UPDATER(GWCE_LV,COEFD,DUMY1,2)
      ENDIF
#endif

      if (subdomainOn.and.enforceBN.eq.2) call enforceEob() ! NCSU Subdomain
      if (subdomainOn.and.enforceBN.eq.2) call enforceEib() ! NCSU Subdomain
      if (subdomainOn.and.enforceBN.eq.2) call enforceGWCELVob() ! NCSU Subdomain

      IF (ILump.eq.0) THEN ! default, fully consistent LHS
C...  JCG ITERATIVE MATRIX SOLVER
         IPARM(1)=ITMAX
         CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &        IWKSP,NW,WKSP,IPARM,RPARM,IER)

         NUMITR=IPARM(1)
         DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ELSE ! lumped LHS
         DO I = 1, NP
            IF (COEFD(I).eq.0.0d0) THEN
               RDIAG = 0.0d0
            ELSE
               RDIAG = 1.0d0 / COEFD(I)
            ENDIF
            ETAS(I) = GWCE_LV(I) * RDIAG
         ENDDO
         NUMITR=0
         DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ENDIF

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif      
     
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN

C**********************************************************************
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
C
      USE SIZES
      USE GLOBAL
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, TotalArea, Areas, NeiTab,
     &                  NNEIGH, SFAC
      USE BOUNDARIES, ONLY : NETA, NFLUXB, NFLUXF, NFLUXGBC, NFLUXRBC,
     &   NFLUXIB, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3
      USE WIND
      USE ITPACKV
      USE NodalAttributes, ONLY : FRIC, Tau0Var, HBREAK, FTHETA, FGAMMA,
     &     IFLINBF, IFNLBF, IFHYBF, EVM, LoadAdvectionState
#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER  NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      REAL(SZ) A00pB00
      REAL(SZ) BCXAvg, BCYAvg
      REAL(SZ) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(SZ) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(SZ) CorifAvg
      REAL(SZ) DPAvg, GDPAvgOAreaIE4
      REAL(SZ) DispX, DispY, DispXAvg, DispYAvg
      REAL(SZ) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(SZ) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(SZ) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(SZ) ESN1, ESN2, ESN3, ESAvg
      REAL(SZ) EVMH, EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(SZ) EVMEle, EVMSmag
      REAL(SZ) GA00DPAvgOAreaIE4
      REAL(SZ) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(SZ) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg, H1, H2
      REAL(SZ) H2OTotalArea
      REAL(SZ) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(SZ) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(SZ) MsFacR
      REAL(SZ) MX, MY, MXAvg, MYAvg
      REAL(SZ) JXAvg, JYAvg
      REAL(SZ) Pr1N1, Pr1N2, Pr1N3
      REAL(SZ) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(SZ) SFacAvg
      REAL(SZ) T0N1,T0N2, T0N3
      REAL(SZ) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(SZ) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(SZ) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(SZ) TiPN1, TiPN2, TiPN3
      REAL(SZ) UV0, UV1, UV2
      REAL(SZ) U1N1,U1N2,U1N3, U1Avg
      REAL(SZ) V1N1,V1N2,V1N3, V1Avg
      REAL(SZ) WSXAvg, WSYAvg

      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) TimeLoc, TimeH

ckmd   Added in parameters for the pc algorithm
      REAL(SZ) BSX0N1, BSX0N2, BSX0N3, BSY0N1
      REAL(SZ) BSY0N2, BSY0N3, BSX0Avg, BSY0Avg
      REAL(SZ) BSX2N1, BSX2N2, BSX2N3, BSY2N1
      REAL(SZ) BSY2N2, BSY2N3, BSX2Avg, BSY2Avg
      REAL(SZ) E2N1,E2N2,E2N3
      REAL(SZ) E0N1SQ, E0N2SQ, E0N3SQ
      REAL(SZ) E2N1SQ, E2N2SQ, E2N3SQ
      REAL(SZ) H0N1, H0N2, H0N3, H00
      REAL(SZ) QX0N1, QX0N2, QX0N3, QY0N1, QY0N2, QY0N3, QX0Avg, QY0Avg
      REAL(SZ) QX2N1, QX2N2, QX2N3, QY2N1, QY2N2, QY2N3, QX2Avg, QY2Avg
      REAL(SZ) Tau0QX0Avg, Tau0QY0Avg, Tau0QX2Avg, Tau0QY2Avg
      REAL(SZ) Tau0SpaVar0, Tau0SpaVar2
      REAL(SZ) U0N1,U0N2,U0N3, U0Avg
      REAL(SZ) V0N1,V0N2,V0N3, V0Avg
      REAL(SZ) U2N1,U2N2,U2N3, U2Avg
      REAL(SZ) V2N1,V2N2,V2N3, V2Avg
      REAL(SZ) timewtgwce0,timewtgwce1,timewtgwce2,timeagflag

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
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
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
            H1=DP(I)+IFNLFA*ETA1(I)
            H2=DP(I)+IFNLFA*ETA2(I)
            TK0(I)=FRIC(I)*(IFLINBF + (UV0/H00)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H00)**FTHETA)**(FGAMMA/FTHETA)))
            TK(I)=FRIC(I)*(IFLINBF + (UV1/H1)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
            TK2(I)=FRIC(I)*(IFLINBF + (UV2/H2)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H2)**FTHETA)**(FGAMMA/FTHETA)))
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
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         H2N1=DP(NM1)+IFNLFA*E2N1
         H2N2=DP(NM2)+IFNLFA*E2N2
         H2N3=DP(NM3)+IFNLFA*E2N3
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

         IF (C3D) THEN                         !3D bottom friction
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
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.d0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.d0
         ENDIF
         IF (C3D) THEN                !3D velocity dispersion
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
ckmd
      DO I=1,NP
c         Eta1(I)=Eta2(I)
         Eta2(I)=0.0d0
      END DO

C...  At elevation boundary condition nodes, determine the elevation at
C...  the s+1 time step
C...
C...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0.
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

C...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

      IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TimeLoc-ETIME1)/ETIMINC
         DO I=1,NETA
            NBDI=NBD(I)
            Eta2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
         END DO
      ENDIF


C...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
C...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

C...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
C...  assume that QN is positive into the domain.  QFORCEJ has a -1
C...  built in and the terms are not explicitly negated. Boundary
C...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
C...  flux that  is positive out of the domain.  Therefore, to match
C...  the formulation these terms must be explicitly multiplied by -1.

C...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 +
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.d0

         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)

            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.d0

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.d0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.d0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.d0*QForceJ+QForceI)
          ENDDO
        ENDIF

C...
C...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
C...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
C...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
C...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
C...  EQUATIONS
C...
      DO I=1,NETA
         NBDI=NBD(I)
         ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
         GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
         DO J=2,NNEIGH(NBDI)
            GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &           -ETAS(NBDI)*OBCCOEF(I,J-1)
         END DO
      END DO

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


C----------------------------------------------------------------------
C----------------------------------------------------------------------
      end module gwce
C----------------------------------------------------------------------
C----------------------------------------------------------------------
