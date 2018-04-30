C     
C--------DW for sponge layer
C     
      MODULE SPONGELAYER 
      
      USE GLOBAL, ONLY: STATIM, IFNLFA, ETA2, UU2, VV2, G, QX2, QY2
      USE SIZES, ONLY: SZ
      USE MESH, ONLY: NP, NM, NE, DP, X
!      USE WIND
!      USE ITPACKV
!      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE NodalAttributes, ONLY: LoadAbsLayerSigma, 
     &     SSIGMA_ETA => absorblayer_sigma_eta,
     &     SSIGMA_MNX => absorblayer_sigma_mnx,
     &     SSIGMA_MNY => absorblayer_sigma_mny,
     &     NumNodesAbsLayer, 
!     &     AbsLayerNodesID,  
     &     AbsLayer_Eta_NodesID, 
     &     AbsLayer_MnX_NodesID,  
     &     AbsLayer_MnY_NodesID,
     &     AbsLayerType, AbsLayerNBF,
     &     AbsLayerEtaAMIG, AbsLayerEtaFF, AbsLayerEtaFACE,
     &     AbsLayerEtaEMO, AbsLayerEtaEFA,
     &     AbsLayerQAMIG, AbsLayerQFF, AbsLayerQFACE,
     &     AbsLayerQxEMO, AbsLayerQxEFA, 
     &     AbsLayerQyEMO, AbsLayerQyEFA !, GWCE_spg_alt 

      USE QUADRATURETRI

#ifdef CMPI
      USE MESSENGER
#endif
      IMPLICIT NONE

      !c For AbsLayerType == 1  
      REAL(SZ), ALLOCATABLE:: pert_eta_abslayer(:), pert_vel_abslayer(:)
      
      REAL(SZ), ALLOCATABLE:: eta2_abslayer(:), eta1_abslayer(:),
     &                     eta0_abslayer(:)  ! (s + 1, s, s-1) !      

      REAL(SZ), ALLOCATABLE:: uu2_abslayer(:), uu1_abslayer(:) ! (s+1,s) !
      REAL(SZ), ALLOCATABLE:: vv2_abslayer(:), vv1_abslayer(:) ! (s+1,s) !

      REAL(SZ), private, ALLOCATABLE:: WorkVecEta(:)
      REAL(SZ), private, ALLOCATABLE:: WorkVecUU(:)
      REAL(SZ), private, ALLOCATABLE:: WorkVecVV(:)
      REAL(SZ), private, ALLOCATABLE:: ArbEta0(:), ArbEta(:)
      REAL(SZ), private, ALLOCATABLE:: ArbUU0(:), ArbUU(:)
      REAL(SZ), private, ALLOCATABLE:: ArbVV0(:), ArbVV(:)

      REAL(SZ) :: SPTIMINC, SPTIME1, SPTIME2

      LOGICAL, allocatable:: spgflag(:)

C------------- accurate calculation of \int_{} \sigma \phi_{i} \phi_{j} dx
      TYPE (MATARR), allocatable:: MsElmFh(:)
      TYPE (VECARR), allocatable:: LumpMsElmFh(:)
C-------------
 
C-------------
C-------------DW: sponge layers (2015)
      REAL(SZ), ALLOCATABLE:: MCOEFSPG(:,:) !c mass matrix 
      REAL(SZ), ALLOCATABLE:: MCOEFSPGD(:) !c lumped mass matrix
C------------
C------------

      INTEGER, parameter:: sponge_dis_mthd = 1  !c 1 -- solve GWCE c!
!c 0 -- Operator Spliling
      LOGICAL, parameter:: adjsigval = .FALSE. 
      LOGICAL :: NO_MET_IN_SPONGE = .FALSE.  
      LOGICAL :: NO_BPG_IN_SPONGE = .FALSE.  
      
     
      CONTAINS

C     
C     Allocate space for array used in the sponge layer through an
C     operator spliting technique: 
C     - precomputed mass matrix and lumped-mass matrix
C     
      SUBROUTINE ALLOC_MAINSPG_LUMPED()
      IMPLICIT NONE

      ALLOCATE( MCOEFSPGD(MNP) ) 

      RETURN
      END SUBROUTINE  ALLOC_MAINSPG_LUMPED

      SUBROUTINE ALLOC_MAINSPG()
      IMPLICIT NONE

      ALLOCATE( MCOEFSPG(MNP,MNEI) ) 

      RETURN
      END SUBROUTINE  ALLOC_MAINSPG
C--------DW
      
      SUBROUTINE Adjust_Sponge_Sigma( DT )
      IMPLICIT NONE

      REAL(SZ):: DT

      REAL(8):: SIGMAX
      LOGICAL, SAVE:: first = .TRUE. 

      SIGMAX = 1.D0 ;
      
      IF ( adjsigval )THEN
         SIGMAX = 0.9D0/DT ;
      END IF

      IF ( first ) THEN
!!$        SSIGMA = SIGMAX*SSIGMA ;
        SSIGMA_ETA = SIGMAX*SSIGMA_ETA ;

        first = .FALSE. ;
      
        SELECT CASE( sponge_dis_mthd )
        CASE (0)
          write(16,*) "Sponge Layer: use an operator 
     & splitting in discretization" 
        CASE (1)
          write(16,*) "Sponge Layer: implement in GWCE
     & and Non conservative Momentum equation" 
        END SELECT

      END IF
      
      RETURN       
      END SUBROUTINE Adjust_Sponge_Sigma
      
C     
C.....Simplest implementation: lumed both lhs and rhs matrix
C     
C     
      SUBROUTINE SPONGE_OPSPLIT0( IT, DT, TimeLoc, TimeH)
C     
      IMPLICIT NONE

C     
C     Dummy variables     
C     
      INTEGER:: IT
      REAL(SZ):: DT, TimeLoc, TimeH
C     
C     Local variables
C     
      INTEGER:: I
      REAL(SZ):: DTD2, SG(3), CNLHS(3), CNRHS, SIGMAX, H2    
      REAL(SZ), parameter:: smalleps = 1.D-12 ; 

      DTD2 = DT/2.D0 ;
!!$ SIGMAX = 0.8/DT ;
!!$ SIGMAX = 1.D0/DT ;
      SIGMAX = 1.D0 ; 

      DO I = 1, NP
         SG(1) = SIGMAX*SSIGMA_ETA(I,1);
         SG(2) = SIGMAX*SSIGMA_MNX(I,1);
         SG(3) = SIGMAX*SSIGMA_MNY(I,1);

         IF ( sum(abs(SG)) > smalleps ) THEN
C.....Crank-Nicolson 
C     CNLHS = 1.D0/(1.D0 + DTD2*SG) ; 
C     CNRHS = (1.D0 - DTD2*SG) ;

C.....Backward Euler
            CNLHS = 1.D0/(1.D0 + DT*SG) ;
            CNRHS = 1.D0 ; 
C     
C     Update \zeta
C     d(\zeta)/dt = - sigma*\zeta  
            ETA2(I) = CNLHS(1)*CNRHS*ETA2(I) ;

C     
C     Update U,V
C     dU/dt = -sigma*U
C     dV/dt = -sigma*U
C     
            UU2(I) = CNLHS(2)*CNRHS*UU2(I) ; 
            VV2(I) = CNLHS(3)*CNRHS*VV2(I) ; 


C     update QX
            H2=DP(I)+IFNLFA*ETA2(I)
            QX2(I)=UU2(I)*H2
            QY2(I)=VV2(I)*H2

         END IF
      END DO

      RETURN ;
      END SUBROUTINE SPONGE_OPSPLIT0


      SUBROUTINE FLAGSPONGEELEM()
        IMPLICIT NONE
        
        INTEGER:: IE, idx(3)
        REAL(8), parameter:: SMALL = 1.0D-10 ;

        IF ( NumNodesAbsLayer(1) > 0 ) THEN
           OPEN(unit=550,FILE='spgelment.dat') ;
 
           !     Flag element inside a sponge layer
           allocate( spgflag(NE) ) ; spgflag = .FALSE. ;
           DO IE = 1, NE
              idx = NM(IE,1:3) ; 

              IF ( sum(abs(SSIGMA_ETA(idx,1))) > SMALL ) THEN
                 spgflag(IE) = .TRUE. ;

                 WRITE(550,*) IE ; 
              END IF  
           END DO
           CLOSE(550) ; 

        END IF

        RETURN ;
      END SUBROUTINE FLAGSPONGEELEM
 
      SUBROUTINE ComputeSpongeMsElmMat( pc )
        IMPLICIT NONE

        !c 
        INTEGER:: pc

        !c local variables c!
        INTEGER:: II, IE, idx(3)
        REAL (SZ):: sigmaval(3)

        LOGICAL, SAVE:: FirstEnter = .TRUE. 


        IF ( NumNodesAbsLayer(1) > 0 ) THEN
           IF ( FirstEnter ) THEN
              ALLOCATE(MsElmFh(NE)) ;
              ALLOCATE(LumpMsElmFh(NE)) ;  
           END IF
           CALL GetDefaultCub2D( pc ) ;

           DO IE = 1, NE
              IF ( spgflag(IE) ) THEN

                 idx = NM(IE,1:3) ; 
                 sigmaval = SSIGMA_ETA(idx,1) ;

                 IF ( FirstEnter ) THEN
                    nullify(MsElmFh(IE)%ARRVAL) ;
                    nullify(LumpMsElmFh(IE)%VECVAL) ;  
                 END IF

                 CALL AllocMatArr( MsElmFh(IE), 3 ) ;
                 CALL AllocVecArr( LumpMsElmFh(IE), 3) ; 

                 CALL CompElmMsfh( MsElmFh(IE)%ARRVAL, sigmaval ) ;
                 DO II = 1, 3
                    LumpMsElmFh(IE)%VECVAL(II) = 
     &                   sum(MsElmFh(IE)%ARRVAL(II,1:3)) ;
                 END DO
              END IF
           END DO
        END IF
        FirstEnter = .FALSE. ;  

        RETURN ; 
      END SUBROUTINE ComputeSpongeMsElmMat

      SUBROUTINE SpongeLayerRelatedPrep( ) 
        IMPLICIT NONE

        INTEGER:: I, J, K, IE, NODEID
        INTEGER:: NAL, IMIN(1)
     
        INTEGER:: SSDUM(4)
        INTEGER:: allowabstype(4) = (/ 0, 1, -1, 2 /) ; 
        REAL(SZ), parameter:: SMALL = 1.D-10  ;

        REAL(SZ), pointer:: SSIGM_PT(:,:)
        INTEGER, pointer:: AbsLayerid_PTS(:)
        !
        IF ( NumNodesAbsLayer(1) > 0 ) THEN
          
           !IF (  GWCE_spg_alt ) THEN
           !   !c GWCE_spg_alt == .TRUE.                    
           !   !c   Use a more accurate integration for the 
           !   !c   sponge_layer term                       
           !   CGWCE_New = .FALSE. ;
           !   CGWCE_New_SPG_ALT_SOL = .TRUE. ; 
           !END IF

           NAL = UBOUND(allowabstype,1) ; 

           SSDUM = 0 ;

           DO K = 1, 3
              SELECT CASE(K)
              CASE (1)
                SSIGM_PT => SSIGMA_ETA ;
                AbsLayerid_PTS => AbsLayer_Eta_NodesID ;
              CASE (2)
                SSIGM_PT => SSIGMA_MNX ;
                AbsLayerid_PTS => AbsLayer_MnX_NodesID ;
              CASE (3)
                SSIGM_PT => SSIGMA_MNY ;
                AbsLayerid_PTS => AbsLayer_MnY_NodesID ;
              END SELECT

              DO I = 1, NumNodesAbsLayer(K)
                 NODEID = AbsLayerid_PTS(I) ;

                 DO J = 1, NAL
                    SSDUM(J) = SSDUM(J) + 
     &                   NINT(SSIGM_PT(NODEID,2) - allowabstype(J)) ; 
                 END DO   

              END DO
           END DO

           ! allocate memormy for external forcing function !
           !  MNP - in global.F, maximum no. of nodes       !
           allocate( eta2_abslayer(MNP) ) ; eta2_abslayer = 0.D0 ;
           allocate( eta1_abslayer(MNP) ) ; eta1_abslayer = 0.D0 ;
           allocate( eta0_abslayer(MNP) ) ; eta0_abslayer = 0.D0 ;

           allocate( uu2_abslayer(MNP) ) ;  uu2_abslayer = 0.D0 ;
           allocate( uu1_abslayer(MNP) ) ;  uu1_abslayer = 0.D0 ;
           allocate( vv2_abslayer(MNP) ) ;  vv2_abslayer = 0.D0 ;
           allocate( vv1_abslayer(MNP) ) ;  vv1_abslayer = 0.D0 ;

           allocate(WorkVecEta(NumNodesAbsLayer(1))) ; 
           allocate(WorkVecUU(NumNodesAbsLayer(2))) ; 
           allocate(WorkVecVV(NumNodesAbsLayer(3))) ; 


           IMIN = MINLOC(abs(SSDUM)) ;
           IF ( abs(SSDUM(IMIN(1))) > 0 ) THEN
              PRINT*, "Error: Unknown sponge-layer type" ; 
              STOP ; 
           ELSE
              AbsLayerType = allowabstype(IMIN(1)) ;
           END IF

           SELECT CASE(AbsLayerType) 
           CASE (0)
             write(16,*) "Sponge Layer: Fully Absorbing"
             AbsLayerNBF = 0 ;
             NO_MET_IN_SPONGE = .true.
             NO_BPG_IN_SPONGE = .true.
             write(16,*) "MET and BPG turned off in sponge!"
           CASE (1)
             write(16,*) "Sponge Layer: Generating(Tides)-Absorbing"
             NO_MET_IN_SPONGE = .true.
             NO_BPG_IN_SPONGE = .true.
             write(16,*) "MET and BPG turned off in sponge!"
             CALL READFORT53001()
             CALL READFORT54001()
           CASE (2)
             write(16,*) "Sponge Layer: Generating(Tides+Arb)-Absorbing"
             NO_MET_IN_SPONGE = .true.
             NO_BPG_IN_SPONGE = .true.
             write(16,*) "MET and BPG turned off in sponge!"
             CALL READFORT53001()
             CALL READFORT54001()
             CALL READFORT2001(.true.)
           CASE (-1)
             write(16,*) "Sponge Layer: Generating(Arbitrary)-Absorbing"
             write(16,*) "MET and BPG kept on in sponge!"
             CALL READFORT2001(.true.)
             AbsLayerNBF = 0 ;
           END SELECT
        END IF
        !
        RETURN ;
      END SUBROUTINE SpongeLayerRelatedPrep

      !
      ! NOTE:
      !   fort.53001 and fort.54001 
      !  have the same format as fort.53 and fort.54, respectively.
      !
      SUBROUTINE READFORT53001()
        IMPLICIT NONE
         
        INTEGER:: FUNIT, IOS

        INTEGER:: KK, JJ, NPTMP, KTMP
        CHARACTER (LEN=300):: IGNOREMSG 
        REAL(SZ):: PIL, DEGRAD  

        PIL = ACOS(-1.D0) ;
        DEGRAD = PIL/180.D0 ; 

        FUNIT = 531 ; 
        OPEN(UNIT = FUNIT, 
     &     FILE=TRIM(INPUTDIR)//'/'//'fort.53001', IOSTAT = IOS ) ;


        READ(FUNIT,*) AbsLayerNBF ; 
        IF ( .NOT. ALLOCATED(AbsLayerEtaAMIG) ) THEN
           ALLOCATE(AbsLayerEtaAMIG(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerEtaFF) ) THEN
           ALLOCATE(AbsLayerEtaFF(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerEtaFACE) ) THEN
           ALLOCATE(AbsLayerEtaFACE(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(pert_eta_abslayer) ) THEN
           ALLOCATE(pert_eta_abslayer(AbsLayerNBF)) ;
        END IF

        DO KK = 1, AbsLayerNBF
           READ(FUNIT,*)  AbsLayerEtaAMIG(KK), AbsLayerEtaFF(KK), 
     &            AbsLayerEtaFACE(KK), IGNOREMSG ;

           AbsLayerEtaFACE(KK) = AbsLayerEtaFACE(KK)*DEGRAD ; 

           IF ( AbsLayerEtaAMIG(KK) .EQ. 0 ) THEN
              pert_eta_abslayer(kk) = 0.D0 ;
           ELSE
              pert_eta_abslayer(kk) = 2.D0*PIL/AbsLayerEtaAMIG(KK) ;
           END IF
        END DO

        READ(FUNIT,*) NPTMP
        IF ( NPTMP /= NumNodesAbsLayer(1) ) THEN
           PRINT*, "Error in READFORT53001(): 1" ; 
           STOP ; 
        END IF
        
        IF ( .NOT. ALLOCATED(AbsLayerEtaEMO) ) THEN
           ALLOCATE(AbsLayerEtaEMO(NPTMP,AbsLayerNBF)) ;
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerEtaEFA) ) THEN
           ALLOCATE(AbsLayerEtaEFA(NPTMP,AbsLayerNBF)) ;
        END IF       

        DO KK = 1, NPTMP
           READ(FUNIT,*) KTMP
           
           ! PRINT*, KTMP, AbsLayer_Eta_NodesID(KK)
           IF ( KTMP /= AbsLayer_Eta_NodesID(KK) ) THEN
             PRINT*, "Error in READFORT53001(): 2" ; 
             STOP ; 
           END IF
           
           DO JJ = 1, AbsLayerNBF
              READ(FUNIT,*) AbsLayerEtaEMO(KK,JJ), AbsLayerEtaEFA(KK,JJ) ;

              !
              ! AbsLayerEtaEFA(KK,JJ) = AbsLayerEtaEFA(KK,JJ)*DEGRAD ;
           END DO

        END DO
        AbsLayerEtaEFA = AbsLayerEtaEFA*DEGRAD ;


        CLOSE(FUNIT) ;         
        WRITE(16,*) "DONE READING FORT.53001" 

        RETURN ;
      END SUBROUTINE READFORT53001

      SUBROUTINE READFORT54001()
        IMPLICIT NONE
         
        INTEGER:: FUNIT, IOS

        INTEGER:: KK, JJ, NPTMP, KTMP, NBF
        CHARACTER (LEN=300):: IGNOREMSG 
        REAL(SZ):: PIL, DEGRAD  

        PIL = ACOS(-1.D0) ;
        DEGRAD = PIL/180.D0 ; 

        FUNIT = 541 ; 
        OPEN(UNIT = FUNIT, 
     &     FILE=TRIM(INPUTDIR)//'/'//'fort.54001', IOSTAT = IOS ) ;
                

        READ(FUNIT,*) NBF ;
        IF ( NBF /= AbsLayerNBF ) THEN
           PRINT*, "Error: READFORT54001, inconsistant o. of freq" ;
           STOP ; 
        END IF

        IF ( .NOT. ALLOCATED(AbsLayerQAMIG) ) THEN
           ALLOCATE(AbsLayerQAMIG(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerQFF) ) THEN
           ALLOCATE(AbsLayerQFF(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerQFACE) ) THEN
           ALLOCATE(AbsLayerQFACE(AbsLayerNBF)) ; 
        END IF
        IF ( .NOT. ALLOCATED(pert_vel_abslayer) ) THEN
           ALLOCATE(pert_vel_abslayer(AbsLayerNBF)) ;
        END IF

        DO KK = 1, AbsLayerNBF
           READ(FUNIT,*)  AbsLayerQAMIG(KK), AbsLayerQFF(KK), 
     &            AbsLayerQFACE(KK), IGNOREMSG ;

           AbsLayerQFACE(KK)=AbsLayerQFACE(KK)*DEGRAD ;

           IF ( AbsLayerQAMIG(KK) .EQ. 0 ) THEN
!!           IF ( AbsLayerEtaAMIG(KK) .EQ. 0 ) THEN
!!              pert_eta_abslayer(kk) = 0.D0 ;
              pert_vel_abslayer(kk) = 0.D0 ;
           ELSE
              pert_vel_abslayer(kk) = 2.D0*PIL/AbsLayerQAMIG(KK) ;
           END IF

        END DO

        READ(FUNIT,*) NPTMP
        IF ( NPTMP /= NumNodesAbsLayer(2) .AND. 
     &                NPTMP /= NumNodesAbsLayer(3) ) THEN
           PRINT*, "Error in READFORT54001(): 1" ; 
           STOP ; 
        END IF
        
        IF ( .NOT. ALLOCATED(AbsLayerQxEMO) ) THEN
           ALLOCATE(AbsLayerQxEMO(NPTMP,AbsLayerNBF)) ;
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerQxEFA) ) THEN
           ALLOCATE(AbsLayerQxEFA(NPTMP,AbsLayerNBF)) ;
        END IF       
        IF ( .NOT. ALLOCATED(AbsLayerQyEMO) ) THEN
           ALLOCATE(AbsLayerQyEMO(NPTMP,AbsLayerNBF)) ;
        END IF
        IF ( .NOT. ALLOCATED(AbsLayerQyEFA) ) THEN
           ALLOCATE(AbsLayerQyEFA(NPTMP,AbsLayerNBF)) ;
        END IF       

        DO KK = 1, NPTMP
           READ(FUNIT,*) KTMP
           
           IF ( KTMP /= AbsLayer_MnX_NodesID(KK) ) THEN
             PRINT*, "Error in READFORT53001(): 2" ; 
             STOP ; 
           END IF

           DO JJ = 1, AbsLayerNBF
              READ(FUNIT,*) AbsLayerQxEMO(KK,JJ), 
     &          AbsLayerQxEFA(KK,JJ),
     &          AbsLayerQyEMO(KK,JJ), 
     &          AbsLayerQyEFA(KK,JJ) ;  

!              AbsLayerQxEFA(KK,JJ) = AbsLayerQxEFA(KK,JJ)*DEGRAD ;
!              AbsLayerQyEFA(KK,JJ) = AbsLayerQyEFA(KK,JJ)*DEGRAD ;
           END DO

        END DO
        AbsLayerQxEFA = AbsLayerQxEFA*DEGRAD ;
        AbsLayerQyEFA = AbsLayerQyEFA*DEGRAD ;    

        CLOSE(FUNIT) ;         

        WRITE(16,*) "DONE READING FORT.54001"

        RETURN ;
      END SUBROUTINE READFORT54001
      
      SUBROUTINE READFORT2001(FIRST_CALL)
        IMPLICIT NONE
         
        INTEGER:: FUNIT, IOS
        LOGICAL:: FIRST_CALL

        INTEGER:: KK, NPTMP

        FUNIT = 201 ; 
        NPTMP = NumNodesAbsLayer(1)
        IF (FIRST_CALL) THEN
           OPEN(UNIT = FUNIT, 
     &         FILE=TRIM(INPUTDIR)//'/'//'fort.2001', IOSTAT = IOS )
           READ(FUNIT,*) SPTIMINC
           SPTIME1 = STATIM*86400.d0
           SPTIME2 = SPTIME1 + SPTIMINC
           allocate(ArbEta0(NumNodesAbsLayer(1))); 
           allocate(ArbUU0(NumNodesAbsLayer(2))) ; 
           allocate(ArbVV0(NumNodesAbsLayer(3))) ; 
           allocate(ArbEta(NumNodesAbsLayer(1))); 
           allocate(ArbUU(NumNodesAbsLayer(2))) ; 
           allocate(ArbVV(NumNodesAbsLayer(3))) ; 
           ArbEta0 = 0.D0 ; ArbUU0 = 0.D0 ; ArbVV0 = 0.D0 ;
           ! reading in new value
           DO KK = 1,NPTMP
              READ(FUNIT,*) ArbEta0(KK), ArbUU0(KK), ArbVV0(KK)
           END DO
           WRITE(16,*) "READ IN FIRST TIMESTEP DATA FOR FORT.2001" 
        ELSE
           ! Updating "0" Arb vectors
           ArbEta0 = ArbEta 
           ArbUU0  = ArbUU 
           ArbVV0  = ArbVV 
        ENDIF
        ! Initialzing new Arb vectors
        ArbEta = 0.D0 ; ArbUU = 0.D0 ; ArbVV = 0.D0 ;
        ! reading in new value
        DO KK = 1,NPTMP
           READ(FUNIT,*) ArbEta(KK), ArbUU(KK), ArbVV(KK)
        END DO

        WRITE(16,*) "READ IN NEW TIMESTEP DATA FOR FORT.2001" 
        RETURN ;
      END SUBROUTINE READFORT2001

      SUBROUTINE GETABSLAYEREXT(IT, TimeLoc, TimeH, RampVal )
        IMPLICIT NONE

        INTEGER:: IT
        REAL(SZ):: TimeLoc, TimeH
        REAL(SZ):: RampVal

        INTEGER:: I, J
        REAL(SZ):: ARG, ARGT, RFF, NCYC, DTRATIO

C.............................DW, generating layer ..............
        IF ( LoadAbsLayerSigma .AND. SUM(NumNodesAbsLayer) > 0 ) THEN
          SELECT CASE(AbsLayerType)
          CASE (0)
            !  Absoribing BC !
          CASE (1,2)
            !
            !  Harmonic forcing
            !
            !  \zeta in the sponge Layer
            WorkVecEta = 0.D0 ;
            DO J=1, AbsLayerNBF
              IF(pert_eta_abslayer(J).EQ.0.) THEN
                  NCYC=0
              ELSE
#ifdef IBM
                 NCYC=INT(timeh/PERT_ETA_ABSLAYER(J),KIND(0.0d0))
#else 
                 NCYC=INT(timeh/PERT_ETA_ABSLAYER(J))
#endif
              ENDIF

              ARGT=AbsLayerEtaAMIG(J)*(timeh-NCYC*pert_eta_abslayer(J))+
     &             AbsLayerEtaFACE(J) ;
              RFF=AbsLayerEtaFF(J)*RampVal ;

              DO I=1, NumNodesAbsLayer(1) 
                  ARG = ARGT - AbsLayerEtaEFA(I,J) ;
                  WorkVecEta(I) = WorkVecEta(I) + 
     &                AbsLayerEtaEMO(I,J)*RFF*COS(ARG) ;

              END DO
            END DO 
            Eta2_AbsLayer(AbsLayer_Eta_NodesID) = WorkVecEta ;

            !  (u,v) in the sponge layer 
            WorkVecUU = 0.D0 ;
            WorkVecVV = 0.D0 ;
            DO J=1, AbsLayerNBF
              IF(pert_vel_abslayer(J).EQ.0.) THEN
                 NCYC=0
              ELSE
#ifdef IBM
                 NCYC=INT(timeh/pert_vel_abslayer(J),KIND(0.0d0))
#else
                 NCYC=INT(timeh/pert_vel_abslayer(J))
#endif
              ENDIF

              ARGT=AbsLayerQAMIG(J)*(timeh-NCYC*pert_vel_abslayer(J))+
     &             AbsLayerQFACE(J) ;
              RFF=AbsLayerQFF(J)*RampVal ;

              DO I=1, NumNodesAbsLayer(2)
                 !
                 ARG = ARGT - AbsLayerQxEFA(I,J) ;
                 WorkVecUU(I) = WorkVecUU(I) + 
     &             AbsLayerQxEMO(I,J)*RFF*COS(ARG) ;    
                 
              END DO

              DO I = 1, NumNodesAbsLayer(3)
                 !
                 ARG = ARGT - AbsLayerQyEFA(I,J) ;
                 WorkVecVV(I) = WorkVecVV(I) +
     &             AbsLayerQyEMO(I,J)*RFF*COS(ARG) ;
                 !  
              END DO
            END DO
            uu2_AbsLayer(AbsLayer_MnX_NodesID) = WorkVecUU ;            
            vv2_AbsLayer(AbsLayer_MnY_NodesID) = WorkVecVV ; 

          CASE (-1)
             ! read new values in
             IF (TimeLoc.gt.SPTIME2) THEN
                SPTIME1 = SPTIME2
                SPTIME2 = SPTIME1 + SPTIMINC
                CALL READFORT2001(.false.)
             ENDIF
             ! do the linear intepolation
             DTRATIO = (TimeLoc-SPTIME1)/SPTIMINC
             eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*(ArbEta0 +
     &                                    DTRATIO*(ArbEta - ArbEta0))
             uu2_AbsLayer(AbsLayer_MnX_NodesID) = RampVal*(ArbUU0 + 
     &                                    DTRATIO*(ArbUU - ArbUU0))
             vv2_AbsLayer(AbsLayer_MnY_NodesID) = RampVal*(ArbVV0 +
     &                                    DTRATIO*(ArbVV - ArbVV0))
          CASE DEFAULT
           !
          END SELECT
          ! Add arbitrary to the tide if case is 2
          IF (AbsLayerType.eq.2) THEN
             ! read new values in
             IF (TimeLoc.gt.SPTIME2) THEN
                SPTIME1 = SPTIME2
                SPTIME2 = SPTIME1 + SPTIMINC
                CALL READFORT2001(.false.)
             ENDIF
             ! do the linear intepolation 
             DTRATIO = (TimeLoc-SPTIME1)/SPTIMINC
             ! adding onto the tide
             eta2_AbsLayer(AbsLayer_Eta_NodesID) = 
     &                            eta2_AbsLayer(AbsLayer_Eta_NodesID) + 
     &                    RampVal*(ArbEta0 + DTRATIO*(ArbEta - ArbEta0))
             uu2_AbsLayer(AbsLayer_MnX_NodesID) = 
     &                             uu2_AbsLayer(AbsLayer_MnX_NodesID) +
     &                       RampVal*(ArbUU0 + DTRATIO*(ArbUU - ArbUU0))
             vv2_AbsLayer(AbsLayer_MnY_NodesID) = 
     &                             vv2_AbsLayer(AbsLayer_MnY_NodesID) +
     &                       RampVal*(ArbVV0 + DTRATIO*(ArbVV - ArbVV0))
          ENDIF

        END IF
C...............................................................

        RETURN 
      END SUBROUTINE GETABSLAYEREXT

      SUBROUTINE SPONGE_SHIFT_SOLN( )
        IMPLICIT NONE

        IF ( NumNodesAbsLayer(1) > 0 ) THEN
           eta0_AbsLayer(AbsLayer_Eta_NodesID) = eta1_AbsLayer(AbsLayer_Eta_NodesID) ;
           eta1_AbsLayer(AbsLayer_Eta_NodesID) = eta2_AbsLayer(AbsLayer_Eta_NodesID) ; 
          
!!$           uu1_AbsLayer(AbsLayer_MnX_NodesID) = uu2_AbsLayer(AbsLayer_MnX_NodesID) ;
!!$           vv1_AbsLayer(AbsLayer_MnY_NodesID) = vv2_AbsLayer(AbsLayer_MnY_NodesID) ;
        END IF

        IF ( NumNodesAbsLayer(2) > 0 ) THEN
           uu1_AbsLayer(AbsLayer_MnX_NodesID) = uu2_AbsLayer(AbsLayer_MnX_NodesID) ;
        END IF

        IF ( NumNodesAbsLayer(3) > 0 ) THEN
           vv1_AbsLayer(AbsLayer_MnY_NodesID) = vv2_AbsLayer(AbsLayer_MnY_NodesID) ;
        END IF

        
        RETURN 
      END SUBROUTINE SPONGE_SHIFT_SOLN

c$$$      
c$$$      !
c$$$      ! Propagting a Gaussin 
c$$$      !
c$$$      SUBROUTINE LAYERSOLTEST1( IT, TimeLoc, TimeH, RampVal )
c$$$        IMPLICIT NONE 
c$$$
c$$$        INTEGER:: IT
c$$$        REAL(SZ):: TimeLoc, TimeH
c$$$        REAL(SZ):: RampVal
c$$$        
c$$$        REAL(SZ):: bb, R2, L, xc 
c$$$        REAL(SZ):: CVEL, LS, DY, XX, YY, CVV
c$$$        INTEGER:: I, NIDX
c$$$
c$$$        bb = 74.D0 ;
c$$$        R2 = 250*100.D0 ;
c$$$        R2 = R2*R2 ; 
c$$$        
c$$$        L = 500*1000.D0 ; 
c$$$        DY = L/25.D0 ; 
c$$$        XC = 2*L + 4*DY ; 
c$$$
c$$$        CVEL = sqrt(G*bb) ;  
c$$$        CVV = sqrt(G/bb) ; 
c$$$
c$$$        LS = XC - CVEL*TimeLoc ;
c$$$        DO I=1, NumNodesAbsLayer(1) 
c$$$          NIDX = AbsLayer_Eta_NodesID(I) ; 
c$$$
c$$$          XX = X(NIDX) ;
c$$$          
c$$$          WorkVecEta(I) = exp(-(XX - LS)**2.D0/R2) ;
c$$$        END DO
c$$$        Eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*WorkVecEta ;
c$$$        !
c$$$
c$$$        DO I=1, NumNodesAbsLayer(2) 
c$$$          NIDX = AbsLayer_Mnx_NodesID(I) ; 
c$$$
c$$$          XX = X(NIDX) ;
c$$$          
c$$$          WorkVecUU(I) = exp(-(XX - LS)**2.D0/R2) ;
c$$$        END DO
c$$$        uu2_AbsLayer(AbsLayer_MnX_NodesID) = -RampVal*CVV*WorkVecUU ;        
c$$$        vv2_AbsLayer = 0.D0 ;
c$$$
c$$$        RETURN ;  
c$$$      END SUBROUTINE LAYERSOLTEST1 
C$$$


      !
      ! Propagting a Gaussin 
      !
      SUBROUTINE LAYERSOLTEST1( IT, TimeLoc, TimeH, RampVal )
        IMPLICIT NONE 

        INTEGER:: IT
        REAL(SZ):: TimeLoc, TimeH
        REAL(SZ):: RampVal
        
        REAL(SZ):: bb, R2, L, xc, LDA, PI 
        REAL(SZ):: CVEL, LS, DY, XX, YY, CVV
        INTEGER:: I, NIDX

        bb = 74.D0 ;
        R2 = 250*100.D0 ;
        R2 = R2*R2 ; 
        
        L = 500*1000.D0 ; 
        DY = L/25.D0 ; 

        XC = 2*L + L/4 ; 
        LDA = 4.0/L ; 

        CVEL = sqrt(G*bb) ;  
        CVV = sqrt(G/bb) ; 
        !c
        PI = acos(-1.D0) ; 

        LS = XC - CVEL*TimeLoc ;
        DO I=1, NumNodesAbsLayer(1) 
          NIDX = AbsLayer_Eta_NodesID(I) ; 
          XX = X(NIDX) - LS
          
          IF ( XX > -1.D0-10 ) THEN
             WorkVecEta(I) = SIN( 2.D0*PI*XX*LDA ) ; 
          END IF
        END DO
        Eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*WorkVecEta ;
        !

        DO I=1, NumNodesAbsLayer(2) 
          NIDX = AbsLayer_Mnx_NodesID(I) ; 

          XX = X(NIDX) - LS ;
          
          IF ( XX > -1.D-10 ) THEN
              WorkVecUU(I) =  SIN( 2.D0*PI*XX*LDA ) ; 
          END IF
        END DO

        uu2_AbsLayer(AbsLayer_MnX_NodesID) = -RampVal*CVV*WorkVecUU ;        
        vv2_AbsLayer = 0.D0 ;

        !!
        RETURN ;  
      END SUBROUTINE LAYERSOLTEST1 

      END MODULE SPONGELAYER 
