c******************************************************************************
c PADCIRC RELEASE VERSION 45.12 03/17/2006                                    *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
c******************************************************************************
c  SUBROUTINE TRANS_3D                                                        *
c                                                                             *
c  Note, the following time stepping coefficients are computed in             *
C     SUBROUTINE READ_INPUT_3DVS and passed in GLOBAL_3DVS.                   *
c                                                                             *
c  DTAlp4      = DelT*Alp4        - weights vert disp term in LHS matrix      *
c  DT1MAlp4    = DelT*(1.-Alp4)   - weights vert disp term in RHS forcing     *
c                                                                             *
c  Note in this subroutine, we utilize two arrays to pass in and out the      *
c  arrays associated with the salinity and temperature fields so this         * 
c  subroutine can be utilized for both fields.                                *
c                                                                             *
c                                                                             *
c  NP - number of nodes in horizontal grid                                    *
c  NFEN - number of nodes in the vertical grid                                *
c******************************************************************************

C   kmd48.33bc added time to information passed into subroutine
      SUBROUTINE TRANS_3D(inarray,diffhoriz,diffvert,outarray,field,time)

C     Include the variables needed from the GLOBAL subroutine. 
C         NODELE, NODECODE, NOFF : used in Casey's wet/dry
C                                  additions (left in for additions later) 
C         NBD, NETA : used in the weak implementation of the 
C                     boundary conditions on the ocean boundary.
C
C  kmd48.33bc add variables from global for boundary conditions and
C             removed ScreenUnit
      USE GLOBAL, ONLY: NODELE, NODECODE, NOFF, NBD, NETA,
     &                  HFLUX, sponge, TTBCTIME1, TTBCTIME2,TTBCRATIO,
     &                  q_heat1, q_heat2, TTBCTIMEINC, BCFLAG_TEMP, SFAC,
     &                  NOPE, NVDLL, NBDV
      USE GLOBAL_3DVS
C  kmd48.33bc add implicit none if not in parallel
#ifdef CMPI
      USE MESSENGER
      IMPLICIT NONE
      REAL(SZ) :: DUMV1(1),DUMV2(1)
#else
      IMPLICIT NONE
#endif

C     Local variables added in from Casey's wet/dry additions
C
      INTEGER :: NCELE
      INTEGER :: TEMPNCELE
      INTEGER :: TEMPSTOP

C     Additional local variables that are passed into or out of the subroutine
      REAL(SZ), DIMENSION(MNP,MNFEN), INTENT(IN) :: inarray
      REAL(SZ), DIMENSION(MNP,MNFEN), INTENT(OUT) :: outarray
      REAL(SZ), INTENT(IN) :: diffhoriz, diffvert, time
      INTEGER, INTENT(IN) :: field 

C     Local variables 
      INTEGER :: IT
      INTEGER :: NEle           !local value of NetTabEle
      INTEGER :: k              !vertical node loop counter (1-bottom, NFEN-surf)
      INTEGER :: NH             !horizontal node loop counter
      INTEGER :: N              !neighbor node loop counter
      INTEGER :: N1,N2,N3,NNFirst !local node numbers used to compute gradients
      INTEGER :: LBP            !value of LBArray_Pointer at present horizontal node
      INTEGER :: NN             !output loop counter
      INTEGER :: tbc            !top boundary condition

C     Local variables used again from the vsmy.F subroutine
      REAL(SZ) :: WSigma(MNFEN) !"sigma" vertical velocity
      REAL(SZ) :: Zk            !z depth of any node k in the vertical
      REAL(SZ) :: DelSig        ! sigma(k+1)-sigma(k)
      REAL(SZ) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
      REAL(SZ) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
      REAL(SZ) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
      REAL(SZ) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
      REAL(SZ) :: SigmaNN       !Sigma value of a neighbor node
      REAL(SZ) :: EtaN1,EtaN2,EtaN3,EtaNFirst !nodal values of NolIFA(Eta1+Eta2)/2
      REAL(SZ) :: hN1,hN2,hN3,hNFirst !nodal values of DP
      REAL(SZ) :: Un,Vn         !real,imaginary components of qn
      REAL(SZ) :: RCLtrans       

      REAL(SZ) :: SFacAvg  !  kmd48.33bc add in spherical factors


      REAL(8) :: KNVnm(MNFEN,3) !integral used in vertical diffusion term for the salinity field
      REAL(8) :: DEtaDT         !time derivative of water surface elev
      REAL(8) :: DEtaDX,DEtaDY  !horizontal derivatives of water surface elev at time level s
      REAL(8) :: DEtaDX2A,DEtaDY2A !(DEtaDX,DEtaDY)*2*Element Area 
      REAL(8) :: DhDX,DhDY      !horizontal derivatives of DP
      REAL(8) :: DhDX2A,DhDY2A  !(DhDX,DhDY)*2*Element Area
      REAL(8) :: TotalArea2     !2*Area of all elements around a node
      REAL(8) :: a1,a2,a3,b1,b2,b3 !elemental coefficients used in horizontal FE method

C     Local variables 
      REAL(8) :: Hstrans             !Total water depth at time level s
      REAL(8) :: HsN2trans           !Total water depth at time level s at local node N2
      REAL(8) :: HsOAMBtrans         !Hs/(a-b)
      REAL(8) :: HsHsOAMBAMBtrans    !(Hs/(a-b))^2
      REAL(8) :: Hsp1trans           !Total water depth at time level s+1
      REAL(8) :: Hsp1OAMBtrans       !Hsp1/(a-b)
      REAL(8) :: Hsp1Hsp1OAMBAMBtrans !(Hsp1/(a-b))^2
!      REAL(8) :: NTVTot(MNFEN)      ! transport diffusion coefficients
      REAL(8) :: tempbcskp1,tempbcsk,topbcs ! top boundary conditions variables for temperature
      REAL(SZ) :: cpwater=4000.0d0          ! variable utilized in the top boundary condition for temperature
      REAL(8) :: saltbcskp1,saltbcsk        ! top boundary conditions variables for salinity (later implementation) 

      REAL(SZ) :: Frtrans(MNFEN)     !right side forcing vector
      REAL(SZ) :: Mkm1trans(MNFEN)   !1st column (k-1) in left side compact storage matrix
      REAL(SZ) :: Mktrans(MNFEN)     !2nd column (k) in left side compact storage matrix
      REAL(SZ) :: Mkp1trans(MNFEN)   !3rd column (k+1) in left side cpmpact storage matrix
      REAL(SZ) :: LAdvectrans(MNFEN) !lateral advection term in transport eqn
      REAL(SZ) :: LDiffusion(MNFEN) !lateral diffusion term in transport eqn
      REAL(SZ) :: VAdvectrans(MNFEN) !vertical advection term in transport eqn
      REAL(SZ) :: VDiffusion(MNFEN) !vertical diffusion term in transport eqn
      !jgf48.50 Declare COMPLEX as size SZ
      COMPLEX(SZ) :: qn,qN1,qN2,qN3,qNFirst !nodal values of q
      REAL(SZ) :: transN1,transN2,transN3,transNFirst !nodal values of transport equation
      REAL(SZ) :: UnDtransDX,VnDtransDY  !derivatives used in lateral advection
      REAL(SZ) :: UnDtransDX2A,VnDtransDY2A !UnDqDX,UnDqDY)*2*Element Area
      REAL(SZ) :: DtransDXDPhiDX2A,DtransDYDPhiDY2A !derivatives used in lateral diffusion calc.
      REAL(SZ) :: DtransDSigmakm1,DtransDSigmakp1 !vertical deriv. of q from k-1,k and k,k+1
      !jgf48.50 Declare COMPLEX as size SZ
      COMPLEX(SZ) :: QSave(MNP,MNFEN) !Used in the testing of the transport equation (remove if works correctly)
      REAL(SZ) :: WSave(MNP,MNFEN) !Used in the testing of the transport equation (remove if works correctly)

C     Variables for the weak implementation of the open ocean boundary condtion
      REAL(SZ) :: WBC(MNP,MNFEN) 
      REAL(SZ) :: RESWBC(MNP,MNFEN)
      REAL(SZ) :: NORM1(3)
      REAL(SZ) :: TOTRES, RESWBCN
      REAL(SZ) :: VELDOTNORM, AVERESN, AVEUN, AVEVN, AVEWN
      REAL(SZ) :: Hsp1OAMBN1, Hsp1OAMBN2
      REAL(SZ) :: Hsp1N1, Hsp1N2
      REAL(SZ) :: XDIST1, YDIST1, ZDIST1, XDIST2, YDIST2, ZDIST2
      REAL(SZ) :: ZVN1, ZVN2, ZVN3, ZVN4
      REAL(SZ) :: XSQUARE, YSQUARE, LENGTHEDGE, ref 
      INTEGER :: NM1, NM2, VM1, VM2
      INTEGER :: NBC, BCSEC             !output loop counter

C  kmd48.33bc added in variables used in the heat flux boundary condition
      INTEGER :: NumofNodes, NOD, J
      REAL(SZ), ALLOCATABLE :: TMP(:,:)

C*************************************************************************************
C     Compute 3D transport equation     

C     Initial step for transport is to determine the weak implementation 
C     of the boundary condition for the ocean. This portion determines the
C     values utilized in the equations based on the previous time step values
C     for salinity or temperature. In the implementation of these boundary
C     conditions, we must determine if we have a influx or outflux boundary
C     conditions. For outflux boundary conditions, we do nothing to the 
C     boundary however, on the influx boundary we evaluate both the boundary
C     condition of the ocean and that of the equations. 
      
C     Reset boundary conditions.
      DO NH=1,NP
         DO k=1,NFEN
            RESWBC(NH,k)=0.d0
         END DO
      END DO 
      
C     Determine the weak boundary conditions that will be used in the
C     assembly process later. Note that the boundary conditions will
C     be utilized in both the salinity and temperature fields.
ckmd49 changed to allow for multiple boundary segments
      DO BCSEC=1,NOPE
        DO NBC=1,NVDLL(BCSEC)-1
C     First obtain the two nodes that with be top two nodes of
C     the rectangle and define two vertical nodes for the depths.
         NM1=NBDV(BCSEC,NBC)
         NM2=NBDV(BCSEC,NBC+1)
         VM1=NFEN
         VM2=NFEN-1
C     Compute the depths at the vertical nodes associated with
C     these two nodes.
         Hsp1N1=DP(NM1)+IFNLFA*ETA2(NM1)
         Hsp1N2=DP(NM2)+IFNLFA*ETA2(NM2)
         Hsp1OAMBN1=Hsp1N1/AMB
         Hsp1OAMBN2=Hsp1N2/AMB
         ZVN1=Hsp1OAMBN1*(Sigma(VM1)-B)-DP(NM1)
         ZVN2=Hsp1OAMBN2*(Sigma(VM1)-B)-DP(NM2)
         ZVN3=Hsp1OAMBN2*(Sigma(VM2)-B)-DP(NM2)
         ZVN4=Hsp1OAMBN1*(Sigma(VM2)-B)-DP(NM1)
         
C     Compute the distance used in the calculations of the normal
C     vector.
         XDIST1=X(NM1)-X(NM2)
         YDIST1=Y(NM1)-Y(NM2)
         ZDIST1=ZVN1-ZVN2
         ZVN1=Hsp1OAMBN1*(Sigma(VM1)-B)-DP(NM1)
         ZVN2=Hsp1OAMBN2*(Sigma(VM1)-B)-DP(NM2)
         ZVN3=Hsp1OAMBN2*(Sigma(VM2)-B)-DP(NM2)
         ZVN4=Hsp1OAMBN1*(Sigma(VM2)-B)-DP(NM1)
         
C     Compute the distance used in the calculations of the normal
C     vector.
         XDIST1=X(NM1)-X(NM2)
         YDIST1=Y(NM1)-Y(NM2)
         ZDIST1=ZVN1-ZVN2
         XDIST2=X(NM2)-X(NM2)
         YDIST2=Y(NM2)-Y(NM2)
         ZDIST2=ZVN3-ZVN2
         
C     Determine normal vector
         CALL CROSSPRODUCT(XDIST1,YDIST1,ZDIST1,XDIST2,
     &        YDIST2,ZDIST2,NORM1)
         
C     Now develop the average velocity values for the dot product
C     to be determined later.
         
         DO k=1,NFEN
            IF (k.NE.NFEN) THEN
               AVEUN=0.25d0*(REAL(Q(NM1,k))+REAL(Q(NM2,k))+
     &              REAL(Q(NM2,k+1))+REAL(Q(NM1,k+1)))
               AVEVN=0.25d0*(AIMAG(Q(NM1,k))+AIMAG(Q(NM2,k))+
     &              AIMAG(Q(NM2,k+1))+AIMAG(Q(NM1,k+1)))
               AVEWN=0.25d0*(WZ(NM1,k)+WZ(NM2,k)+
     &              WZ(NM2,k+1)+WZ(NM1,k+1))
            ELSE IF (k.EQ.NFEN) THEN
               AVEUN=0.25d0*(REAL(Q(NM1,k-1))+REAL(Q(NM2,k-1))+
     &              REAL(Q(NM2,k))+REAL(Q(NM1,k)))
               AVEVN=0.25d0*(AIMAG(Q(NM1,k-1))+AIMAG(Q(NM2,k-1))+
     &              AIMAG(Q(NM2,k))+AIMAG(Q(NM1,k)))
               AVEWN=0.25d0*(WZ(NM1,k-1)+WZ(NM2,k-1)+
     &              WZ(NM2,k)+WZ(NM1,k))
            END IF
            VELDOTNORM=AVEUN*NORM1(1)+AVEVN*NORM1(2)+
     &           AVEWN*NORM1(3)
            XSQUARE=(X(NM1)-X(NM2))*(X(NM1)-X(NM2))
            YSQUARE=(Y(NM1)-Y(NM2))*(Y(NM1)-Y(NM2))
            LENGTHEDGE=SQRT(XSQUARE+YSQUARE)
            IF (field.EQ.2) THEN 
               IF (k.NE.NFEN) THEN
                   ref=0.25d0*(RESSAL(NBC,k)+RESSAL(NBC,k+1)+
     &                RESSAL(NBC+1,k)+RESSAL(NBC+1,k+1))
               ELSE IF (k.EQ.NFEN) THEN
                   ref=0.25d0*(RESSAL(NBC,k-1)+RESSAL(NBC,k)+
     &                RESSAL(NBC+1,k-1)+RESSAL(NBC+1,k))
               END IF
               !ref=32.d0   ! kmd48.33bc changed to 32 as the salinity reference
            ELSE IF (field.EQ.3) THEN 
               IF (k.NE.NFEN) THEN
                   ref=0.25d0*(RESTEMP(NBC,k)+RESTEMP(NBC,k+1)+
     &                RESTEMP(NBC+1,k)+RESTEMP(NBC+1,k+1))
               ELSE IF (k.EQ.NFEN) THEN
                   ref=0.25d0*(RESTEMP(NBC,k-1)+RESTEMP(NBC,k)+
     &                RESTEMP(NBC+1,k-1)+RESTEMP(NBC+1,k))
               END IF
               !ref=10.d0
            END IF  
            IF (VELDOTNORM.LE.0.d0) THEN
               IF (k.NE.NFEN) THEN
                  AVERESN=0.25d0*(inarray(NM1,k)+inarray(NM2,k)+
     &                 inarray(NM2,k+1)+inarray(NM1,k+1))
               ELSE IF (k.EQ.NFEN) THEN
                  AVERESN=0.25d0*(inarray(NM1,k-1)+inarray(NM2,k-1)+
     &                 inarray(NM2,k)+inarray(NM1,k))
               END IF
               TOTRES=AVERESN-ref
C   kmd48.33bc changed this from previous version
               RESWBCN=(VELDOTNORM*(TOTRES)*(1.d0/2.d0)*LVn(k))
               RESWBC(NM1,k)=RESWBC(NM1,k)+RESWBCN
               RESWBC(NM2,k)=RESWBC(NM2,k)+RESWBCN
            ELSE IF (VELDOTNORM.GT.0.d0) THEN
               RESWBC(NM1,k)=RESWBC(NM1,k)+0.d0
               RESWBC(NM2,k)=RESWBC(NM2,k)+0.d0
            END IF
         END DO
        END DO 
      END DO

C   kmd48.33 add information for the top temperature boundary condition 
      IF (field.EQ.3) THEN
         IF (BCFLAG_TEMP.EQ.1) THEN
           IF (TIME.GT.TTBCTIME2) THEN
              TTBCTIME1=TTBCTIME2
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=q_heat2(NumofNodes)
                 READ(38,*) NOD, q_heat2(NumofNodes)
              END DO
           END IF ! ends time if loop
           TTBCRATIO=(TIME-TTBCTIME1)/TTBCTIMEINC
           DO NumofNodes=1,NP
              qsurf(NumofNodes)=qsurfkp1(NumofNodes)
              HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
           END DO
         ELSE IF (BCFLAG_TEMP.EQ.2) THEN
            ALLOCATE (TMP(NP,6))
            IF (TIME.GT.TTBCTIME2) THEN
              TTBCTIME1=TTBCTIME2
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=q_heat2(NumofNodes)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
              DO NumofNodes=1,NP
                q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
              END DO
            END IF ! ends time if loop
            TTBCRATIO=(TIME-TTBCTIME1)/TTBCTIMEINC
            DO  NumofNodes=1,NP
               qsurf(NumofNodes)=qsurfkp1(NumofNodes)
               HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
            END DO
            DEALLOCATE(TMP)
         ELSE IF (BCFLAG_TEMP.EQ.3) THEN
            ALLOCATE (TMP(NP,4))
            IF (TIME.GT.TTBCTIME2) THEN
              TTBCTIME1=TTBCTIME2
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=q_heat2(NumofNodes)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
              DO NumofNodes=1,NP
                q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &               TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
              END DO
            END IF ! ends time if loop
            TTBCRATIO=(TIME-TTBCTIME1)/TTBCTIMEINC
            DO NumofNodes=1,NP
               qsurf(NumofNodes)=qsurfkp1(NumofNodes)
               HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
            END DO
            DEALLOCATE(TMP)
         END IF ! ends the BCFLAG_TEMP flag if loop
      END IF  ! ends the field if loop

      
C     
C     Loop over each horizontal and vertical nodes to compute the transport concentrations
C

      DO NH=1,NP                !loop over horizontal nodes
         
C     Set up some values at the node being worked on
         
         Hstrans  = DP(NH)+IFNLFA*Eta1(NH) !Total depth at previous (s) timestep
         HsOAMBtrans=Hstrans/AMB
         HsHsOAMBAMBtrans=HsOAMBtrans*HsOAMBtrans
         
         Hsp1trans= DP(NH)+IFNLFA*Eta2(NH) !Total depth at present (s+1) timestep
         Hsp1OAMBtrans=Hsp1trans/AMB
         Hsp1Hsp1OAMBAMBtrans=Hsp1OAMBtrans*Hsp1OAMBtrans
         
Casey: Solve for "TEMPNCELE," which is the sum of the "NCELE" values
C     for the elements attached to the current horizontal node.  Note that,
C     for some reason, the k neighbor elements are not contained in the first
C     k registers of "NEITABELE."  (For instance, if a node is connected to four
C     elements, those element numbers might be contained in registers
C     1, 2, 4, and 5.)  So we have to be cute about pulling element numbers
C     from "NEITABELE."
C     
         TEMPNCELE = 0
         TEMPSTOP = 0
         DO K=1,NODELE(NH)
            IF(NEITABELE(NH,K).EQ.0) THEN
               TEMPSTOP = TEMPSTOP + 1
            END IF
            TEMPNCELE = TEMPNCELE +
     &           NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),1))*
     &           NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),2))*
     &           NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),3))
     &           *NOFF(NEITABELE(NH,K+TEMPSTOP))
         END DO
         
Casey: If "TEMPNCELE" is equal to zero, then the current horizontal node is
C     not attached to any active elements.  Thus, the horizontal velocities
C     must be zero, and we can set up the matrix immediately and skip to the
C     tri-diagonal solver.
C
         IF(TEMPNCELE.EQ.0) THEN
            DO K=1,NFEN
               MKM1trans(K) = 0.D0
               MKtrans(K) = 1.D0
               MKP1trans(K) = 0.D0
               FRtrans(K) = 0.D0
            END DO
            GOTO 989
         END IF
         
C        Changed this to the vertical matrix for the vertical diffusion term
C        for both salinity and temperature replaces the EV values.
C        First develop the NTVTot values - Currently only using a constant value

C  kmd48.33bc changed to use either a constant value or the value obtained
C             from the Mellor-Yamada turbulence subroutine
         DO k=1,NFEN
           IF (IEVC.EQ.50) THEN 
            NTVTot(k)=DV(NH,k)
           ELSE
            NTVTot(k)=diffvert
           END IF 
         END DO

C        Integrals for the vertical diffusion terms for both the temperature
C        and salinity values 

         KNVnm(1,1)=0.d0
         KNVnm(1,3)=-0.5d0*(NTVTot(2)+NTVTot(1))/(Sigma(2)-Sigma(1))
         KNVnm(1,2)=-(KNVnm(1,1)+KNVnm(1,3))
         DO k=2,NFEN-1
            KNVnm(k,1)=KNVnm(k-1,3)
            KNVnm(k,3)=-0.5d0*(NTVTot(k+1)+NTVTot(k))/
     &           (Sigma(k+1)-Sigma(k))
            KNVnm(k,2)=-(KNVnm(k,1)+KNVnm(k,3))
         ENDDO
         KNVnm(NFEN,1)=KNVnm(NFEN-1,3)
         KNVnm(NFEN,3)=0.d0      
         KNVnm(NFEN,2)=-(KNVnm(NFEN,1)+KNVnm(NFEN,3))
         
C     Compute time derivative of water surface position

         DEtaDT=(Eta2(NH)-Eta1(NH))/DelT
         
C     Start computing horizontal derivatives of water level, bathymetric
C     depth. Note: TotalArea2 = 2X total elemental area surrounding a node
         
         DEtaDX=0.d0
         DEtaDY=0.d0
         DhDX=0.d0
         DhDY=0.d0
         DEtaDX2A=0.d0
         DhDX2A=0.d0
         DEtaDY2A=0.d0
         DhDY2A=0.d0
         TotalArea2=0.d0

         N1=NH
         EtaN1=IFNLFA*(Eta1(N1)+Eta2(N1))/2.d0
         hN1=DP(N1)

         N2=NeiTab(NH,2)        !operate on 1st neighbor
         EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
         hN2=DP(N2)

         NNFirst=N2             !save these values until end
         EtaNFirst=EtaN2
         hNFirst=hN2

         DO N=3,NNeigh(NH)      !operate on rest of neighbors
            N3=N2               !shift previously computed values
            hN3=hN2             !shift previously computed values
            EtaN3=EtaN2
            N2=NeiTab(NH,N)     !select new neighbor to work on
            EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
            hN2=DP(N2)
            NEle=NeiTabEle(NH,N-2) !element # defined by nodes NH,NN2,NN1

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
            NCELE = NODECODE(NH)*NODECODE(N2)
     &           *NODECODE(N3)*NOFF(NELE)
            IF((NEle.NE.0).AND.(NCELE.NE.0)) THEN !if element is active, compute velocity grads
               TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C   kmd48.33bc add in spherical factors
               SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
               a1=X(N3)-X(N2)
               a2=X(N1)-X(N3)
               a3=X(N2)-X(N1)
               b1=(Y(N2)-Y(N3))*SFacAvg
               b2=(Y(N3)-Y(N1))*SFacAvg
               b3=(Y(N1)-Y(N2))*SFacAvg
               DhDX2A=DhDX2A+(hN1*b1+hN2*b2+hN3*b3)
               DhDY2A=DhDY2A+(hN1*a1+hN2*a2+hN3*a3)
               DEtaDX2A=DEtaDX2A+(EtaN1*b1+EtaN2*b2+EtaN3*b3)
               DEtaDY2A=DEtaDY2A+(EtaN1*a1+EtaN2*a2+EtaN3*a3)
            ENDIF
         END DO

         N3=N2                  !wrap back to beginning to get final contribution
         hN3=hN2
         EtaN3=EtaN2
         N2=NNFirst
         hN2=hNFirst
         EtaN2=EtaNFirst
         NEle=NeiTabEle(NH,NNeigh(NH)-1)

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
         NCELE = NODECODE(NH)*NODECODE(N2)
     &       *NODECODE(N3)*NOFF(NELE)
         IF((NEle.NE.0).AND.(NCELE.NE.0)) THEN
            TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C   kmd48.33bc add in spherical factors
            SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
            a1=X(N3)-X(N2)
            a2=X(N1)-X(N3)
            a3=X(N2)-X(N1)
            b1=(Y(N2)-Y(N3))*SFacAvg
            b2=(Y(N3)-Y(N1))*SFacAvg
            b3=(Y(N1)-Y(N2))*SFacAvg
            DhDX2A=DhDX2A+(hN1*b1+hN2*b2+hN3*b3)
            DhDY2A=DhDY2A+(hN1*a1+hN2*a2+hN3*a3)
            DEtaDX2A=DEtaDX2A+(EtaN1*b1+EtaN2*b2+EtaN3*b3)
            DEtaDY2A=DEtaDY2A+(EtaN1*a1+EtaN2*a2+EtaN3*a3)
         ENDIF

         IF(TotalArea2.NE.0.) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
         ENDIF

C     Finished computing horizontal derivatives of water level,
C     bathymetric depth.

C     Compute the "sigma" vertical velocity from the "z" vertical velocity 

         DO k=1,NFEN
            SigmaMAOAMB=(Sigma(k)-A)/AMB
            SigmaMBOAMB=(Sigma(k)-B)/AMB
            WSigma(k) = WZ(NH,k) - SigmaMBOAMB*DEtaDT
     &           - REAL(q(NH,k))*(SigmaMBOAMB*DEtaDX+SigmaMAOAMB*DhDX)
     &           - AIMAG(q(NH,k))*(SigmaMBOAMB*DEtaDY+SigmaMAOAMB*DhDY)
         ENDDO


C     Start computing advection and diffusion terms in the transport
C     equation at each level in the vertical. The inarray can be either
C     temperature or salinity. 

         DO k=1,NFEN

C     Compute the vertical advection and vertical diffusion terms

            IF(k.EQ.1) THEN
               DtransDSigmakp1=(inarray(NH,k+1)-inarray(NH,k))/
     &              (Sigma(k+1)-Sigma(k))
               VAdvectrans(k)=DtransDsigmakp1*(2.d0*WSigma(k)+
     &              WSigma(k+1))*Inm(k,3)/HsOAMBtrans
               VDiffusion(k)=(inarray(NH,k)*KNVnm(k,2)+inarray(NH,k+1)*
     &              KNVnm(k,3))/HsHsOAMBAMBtrans
            ENDIF
            IF((k.GT.1).AND.(k.LT.NFEN)) THEN
               DtransDSigmakm1=DtransDSigmakp1
               DtransDSigmakp1=(inarray(NH,k+1)-inarray(NH,k))/
     &              (Sigma(k+1)-Sigma(k))
               VAdvectrans(k)=(DtransDSigmakm1*(WSigma(k-1)+
     &              2.d0*WSigma(k))*Inm(k,1)+
     &              DtransDSigmakp1*(2.d0*WSigma(k)+
     &              WSigma(k+1))*Inm(k,3))/HsOAMBtrans
               VDiffusion(k)=(inarray(NH,k-1)*KNVnm(k,1)+inarray(NH,k)*
     &              KNVnm(k,2)+inarray(NH,k+1)*KNVnm(k,3))/
     &              HsHsOAMBAMBtrans
            ENDIF
            IF(k.EQ.NFEN) THEN
               DtransDSigmakm1=DtransDSigmakp1
               VAdvectrans(k)=DtransDSigmakm1*(WSigma(k-1)+2.d0*
     &              WSigma(k))*Inm(k,1)/HsOAMBtrans
               VDiffusion(k)=(inarray(NH,k-1)*KNVnm(k,1)+inarray(NH,k)*
     &              KNVnm(k,2))/HsHsOAMBAMBtrans
            ENDIF
            
C      kmd48.33bc added sponge layer for the vertical advective term
            VAdvectrans(k)=sponge(NH)*IFNLCT*VAdvectrans(k)

C     Compute lateral advection and lateral diffusion terms

            UnDtransDX2A=0.d0
            VnDtransDY2A=0.d0
            DtransDXDPhiDX2A=0.d0
            DtransDYDPhiDY2A=0.d0

            N1=NH               !node 1 is always the central node
            qN1=q(N1,k)       
            transN1=inarray(N1,k)       

            N2=NEITAB(NH,2)     !operate on 1st neighbor
            qN2=q(N2,k)
            transN2=inarray(N2,k)
            
            NNFirst=N2          !save these values until end
            qNFirst=qN2
            transNFirst=transN2

            DO N=3,NNEIGH(NH)   !operate on rest of neighbors
               N3=N2            !shift previously computed values
               qN3=qN2
               transN3=transN2
               N2=NEITAB(NH,N)  !select new neighbor to work on
               qN2=q(N2,k)
               transN2=inarray(N2,k)
               NEle=NeiTabEle(NH,N-2) !element # defined by nodes N1,N2,N3

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
               NCELE = NODECODE(NH)*NODECODE(N2)
     &              *NODECODE(N3)*NOFF(NELE)
               IF((NEle.NE.0).AND.(NCELE.NE.0)) THEN !if element exists, compute terms
                  qn=(qN1+qN2+qN3)/3.d0
                  Un=REAL(qn)
                  Vn=AIMAG(qn)
C   kmd48.33bc add in spherical factors
                  SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  a1=X(N3)-X(N2)
                  a2=X(N1)-X(N3)
                  a3=X(N2)-X(N1)
                  b1=(Y(N2)-Y(N3))*SFacAvg
                  b2=(Y(N3)-Y(N1))*SFacAvg
                  b3=(Y(N1)-Y(N2))*SFacAvg
                  UnDtransDX2A=UnDtransDX2A+Un*(transN1*b1+transN2*b2+
     &                         transN3*b3)
                  VnDtransDY2A=VnDtransDY2A+Vn*(transN1*a1+transN2*a2+
     &                         transN3*a3)
                  DtransDXDPhiDX2A=DtransDXDPhiDX2A
     &                         +(transN1*b1+transN2*b2+transN3*b3)*
     &                         b1/Areas(NEle)
                  DtransDYDPhiDY2A=DtransDYDPhiDY2A
     &                         +(transN1*a1+transN2*a2+transN3*a3)*
     &                         a1/Areas(NEle)
               ENDIF
            END DO

            N3=N2               !wrap back to beginning to get final contribution
            qN3=qN2
            transN3=transN2
            N2=NNFIRST
            qN2=qNFirst
            transN2=transNFirst
            NEle=NeiTabEle(NH,NNeigh(NH)-1)

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
            NCELE = NODECODE(NH)*NODECODE(N2)
     &          *NODECODE(N3)*NOFF(NELE)
            IF((NEle.NE.0).AND.(NCELE.NE.0)) THEN
               qn=(qN1+qN2+qN3)/3.d0
               Un=real(qn)
               Vn=aimag(qn)
C   kmd48.33bc add in spherical factors
               SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
               a1=X(N3)-X(N2)
               a2=X(N1)-X(N3)
               a3=X(N2)-X(N1)
               b1=(Y(N2)-Y(N3))*SFacAvg
               b2=(Y(N3)-Y(N1))*SFacAvg
               b3=(Y(N1)-Y(N2))*SFacAvg
               UnDtransDX2A=UnDtransDX2A+Un*(transN1*b1+transN2*b2+
     &                      transN3*b3)
               VnDtransDY2A=VnDtransDY2A+Vn*(transN1*a1+transN2*a2+
     &                      transN3*a3)
               DtransDXDPhiDX2A=DtransDXDPhiDX2A
     &                      +(transN1*b1+transN2*b2+transN3*b3)*
     &                      b1/Areas(NEle)
               DtransDYDPhiDY2A=DtransDYDPhiDY2A
     &                      +(transN1*a1+transN2*a2+transN3*a3)*
     &                      a1/Areas(NEle)
            ENDIF

            IF(TotalArea2.EQ.0.) THEN
               LAdvectrans(k)=(0.d0,0.d0)
               LDiffusion(k)=(0.d0,0.d0)
            ELSE

C  kmd48.33bc added sponge layer for the lateral advective terms
               LAdvectrans(k)=sponge(NH)*IFNLCT*(UnDtransDX2A+VnDtransDY2A)/
     &                        TotalArea2   
               LDiffusion(k)=3.d0*diffhoriz*(DtransDXDPhiDX2A+
     &                       DtransDYDPhiDY2A)/TotalArea2
            ENDIF

         ENDDO

c     Finished computing advection and diffusion terms in the transport
c     equation at each level in the vertical

C     Need to add the boundary condition for the top layer for
C     temperature and salinity. Currently, I have only included the 
C     temperature term. 
C  kmd48.33bc changed the way the top temperature is evaluated.
      IF (field.EQ.3) THEN
         IF (BCFLAG_TEMP.NE.0) THEN
           qsurfkp1(NH)=HFLUX(NH)
           tempbcskp1=DTAlp4*1.d0/Hsp1OAMBtrans*
     &                (-qsurfkp1(NH)/(cpwater*rhowat0))
           tempbcsk=DT1MAlp4*1.d0/HsOAMBtrans*
     &              (-qsurf(NH)/(cpwater*rhowat0))
           topbcs=tempbcskp1+tempbcsk
         ELSE IF (BCFLAG_TEMP.EQ.0) THEN
           topbcs=0.d0
         END IF
      ELSE
         topbcs = 0.d0
      END IF

C     Need to add in the boundary conditions. Note: We loop through
C     all the nodes and the only ones which have results should be
C     the boundary conditions where vdotn is greater than 0.d0.

      DO k=1,NFEN
         WBC(NH,k)=RESWBC(NH,k)*((3.d0*DelT)/TotalArea2)
      END DO
      
C     Changed the RHS forcing vector and the LHS matrix to use the 
C     information from the transport equations
C     Set up the RHS forcing vector Fr and LHS matrix in compact storage
C     (Mkm1,Mk,Mkp1)

         RCLtrans = DTAlp4/Hsp1Hsp1OAMBAMBtrans

            Frtrans(1) = (1.d0*inarray(NH,1)
     &           -DelT*(LAdvectrans(1)  +LDiffusion(1)))*Inm(1,2)
     &           + (1.d0*inarray(NH,2)
     &           -DelT*(LAdvectrans(2)  +LDiffusion(2)))*Inm(1,3)
     &           - DelT*(VAdvectrans(1))-DT1MAlp4*VDiffusion(1)     
     &           +WBC(NH,1)
            Mkm1trans(1) = 0.d0
            Mktrans(1)   = 1.d0*Inm(1,2) + RCLtrans*KNVnm(1,2)
            Mkp1trans(1) = 1.d0*Inm(1,3) + RCLtrans*KNVnm(1,3)

         DO k=2,NFEN-1
            Frtrans(k) = (1.d0*inarray(NH,k-1)
     &           -DelT*(LAdvectrans(k-1)+LDiffusion(k-1)))*Inm(k,1)
     &           + (1.d0*inarray(NH,k)
     &           -DelT*(LAdvectrans(k)  +LDiffusion(k)))*Inm(k,2)
     &           + (1.d0*inarray(NH,k+1)
     &           -DelT*(LAdvectrans(k+1)+LDiffusion(k+1)))*Inm(k,3)
     &           - DelT*(VAdvectrans(k))-DT1MAlp4*VDiffusion(k)
     &           + WBC(NH,k)
            Mkm1trans(k) = 1.d0*Inm(k,1) + RCLtrans*KNVnm(k,1)
            Mktrans(k)   = 1.d0*Inm(k,2) + RCLtrans*KNVnm(k,2)
            Mkp1trans(k) = 1.d0*Inm(k,3) + RCLtrans*KNVnm(k,3)
         END DO

         Frtrans(NFEN) = (1.d0*inarray(NH,k-1)
     &        -DelT*(LAdvectrans(k-1)+LDiffusion(k-1)))*Inm(k,1)
     &        + (1.d0*inarray(NH,k)
     &        -DelT*(LAdvectrans(k)  +LDiffusion(k)))*Inm(k,2)
     &        - DelT*(VAdvectrans(k))-DT1MAlp4*VDiffusion(k)
     &        + topbcs + WBC(NH,NFEN)
         Mkm1trans(NFEN) = 1.d0*Inm(NFEN,1) + RCLtrans*KNVnm(NFEN,1)
         Mktrans(NFEN)   = 1.d0*Inm(NFEN,2) + RCLtrans*KNVnm(NFEN,2)
         Mkp1trans(NFEN) = 0.d0

C     Set up the boundary conditions for land boundary - currently it 
C     only does nothing to the equations. However, we will have to add
C     other boundary condtions for the land boundary later. 
          LBP=LBArray_Pointer(NH) 
          IF(LBP.GT.0) THEN 
            IF((LBCODEI(LBP).GE. 0)) THEN
!                DO NOTHING
            END IF
          END IF 

C     Decompose and solve the system
 989     CONTINUE ! This addresses the Goto statement associated with the 
                  ! wetting and drying operation
  
       CALL TRIDAG2(Mkm1trans,Mktrans,Mkp1trans,Frtrans,Gammatrans,NFEN)

C      Need to save the transport results by putting them into 
C      outgoing variables 
     
         DO k=1,NFEN
            outarray(NH,k) = Gammatrans(k)
         END DO

      ENDDO

C     Finish loop over horizontal nodes to compute the horizontal velocity

      RETURN
      END SUBROUTINE TRANS_3D

C  This subroutine utilized to solve the matrix - just changed the name
C  to tridag2
C***********************************************************************
C Subroutine tridag2                                                    *
C                                                                      *
C     -----------------------------------------------------------------*
C     |SOLVER FOR A VECTOR U OF LENGTH N FROM A SET OF LINEAR          *
C     |EQUATIONS THAT CONTAINS A TRIDIAGONAL MATRIX                    *
C     |THE FORM IS                                                     *
C     |                                                                *
C     |   * B1 C1  0 ...               *     * U1 *     * R1 *         *
C     |  *                              *   *      *   *      *        *
C     |  *  A2 B2 C2 ...                *   *  U2  *   *  R2  *        *
C     |  *  ...                         * * * ...  * = * ...  *        *
C     |  *           ... An-1 Bn-1 Cn-1 *   * Un-1 *   * Rn-1 *        *
C     |  *                              *   *      *   *      *        *
C     |   *                0   An  Bn  *     * Un *     * Rn *         *
C     |                                                                *
C     |A, B, C, U ARE ARRAYS.                                          *
C     -----------------------------------------------------------------*
C                                                                      *
C    Adapted from Numerical Recipes chapter 2                          *
C***********************************************************************
c
      SUBROUTINE TRIDAG2(A,B,C,R,U,N)
      USE GLOBAL, ONLY: ScreenUnit
      USE GLOBAL_3DVS, ONLY : SZ
      IMPLICIT NONE
      INTEGER :: J, N

      REAL(SZ) :: A(N),B(N),C(N),R(N),U(N)
      REAL(SZ) :: BET,GAM(N)

c     
c     check for zero elements on diagonal
c     
      DO J=1,N
         if(B(j).EQ.0.) then
            write(screenunit,*) 'Problem in Tridag Solver.  ',
     +           'B array value in row ',j,' = 0'
            stop
         endif
      end do
      BET = B(1)
      U(1) = R(1)/BET
      DO J = 2,N
         GAM(J) = C(J-1)/BET
         BET=B(J)-A(J)*GAM(J)
         if (BET.EQ.0) then
            write(screenunit,*) ' Problem in Tridag Solver.  ',
     +           ' BET  = 0.  Solver failed.'
            stop
         endif
         U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J = N-1,1,-1
         U(J) = U(J) - GAM(J+1)*U(J+1)
      END DO
      RETURN
      END

C************************************************************************
Casey: The following subroutine gives the unit outward normal vector    *
C      to two vectors in three-dimensional space.                       *
C************************************************************************

      SUBROUTINE CROSSPRODUCT(X1,Y1,Z1,X2,Y2,Z2,NORM)

      IMPLICIT NONE

      INTRINSIC             :: SQRT

      REAL(8), DIMENSION(3), INTENT(INOUT) :: NORM
      REAL(8)               :: DOTPROD
      REAL(8)               :: MAGN1, MAGN2
      REAL(8), INTENT(IN)   :: X1, X2
      REAL(8), INTENT(IN)   :: Y1, Y2
      REAL(8), INTENT(IN)   :: Z1, Z2

      NORM(1) = Y1 * Z2 - Z1 * Y2
      NORM(2) = Z1 * X2 - X1 * Z2
      NORM(3) = X1 * Y2 - Y1 * X2

      MAGN1        = SQRT( X1 * X1 + Y1 * Y1 + Z1 * Z1 )
      MAGN2        = SQRT( X2 * X2 + Y2 * Y2 + Z2 * Z2 )

      DOTPROD      = (X1/MAGN1) * (X2/MAGN2)
     &             + (Y1/MAGN1) * (Y2/MAGN2)
     &             + (Z1/MAGN1) * (Z2/MAGN2)

      NORM(1) = NORM(1) / ( MAGN1 * MAGN2 * SQRT( 1.D0 - 
     &          DOTPROD * DOTPROD ) )
      NORM(2) = NORM(2) / ( MAGN1 * MAGN2 * SQRT( 1.D0 - 
     &          DOTPROD * DOTPROD ) )
      NORM(3) = NORM(3) / ( MAGN1 * MAGN2 * SQRT( 1.D0 - 
     &          DOTPROD * DOTPROD ) )

      RETURN
      END SUBROUTINE CROSSPRODUCT


