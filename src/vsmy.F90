
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 45.12                                    *
!                                                                             *
!                                                                             *
!  Note, input is read in READ_INPUT_3D located in module READ_INPUT          *
!        initial conditions on density, temperature and/or salinity are read  *
!        in for a cold start in subroutine COLDSTART_3D in module COLDSTART.  *
!                                                                             *
!******************************************************************************


    SUBROUTINE VSSTUP()
    USE GLOBAL_3DVS
    USE MESH, ONLY : NP
    IMPLICIT NONE
    INTEGER :: NH, N

    call setMessageSource("vsstup")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!...  FOR BOTH COLD START AND HOT START, NEED TO INITIALIZE Sigma T FROM
!...  TEMPERATURE AND SALINITY DEPENDING ON VALUE OF IDEN
!...

    IF(ABS(IDen) >= 1) CALL CALC_SIGMAT_3D ()


!...
!...  COMPUTE THE DEPTH-AVERAGED Sigma T FIELD
!...

    DO NH=1,NP
        DASigT(NH)=0.d0
        IF(CBaroclinic) THEN
            DO N=1,NFEN-1
                DASigT(NH)=DASigT(NH)+(Sigma(N+1)-Sigma(N)) &
                *(SigT(NH,N+1)+SigT(NH,N))/2.d0
            ENDDO
            DASigT(NH)=DASigT(NH)/AMB
        ENDIF
    ENDDO

!...
!...  ADDITIONAL RUN SETUP
!...

!...  COMPUTE THE INTEGRALS Inm and LVn (INDEPENDENT OF HORIZONTAL NODE)

    CALL INMInt()

    CALL LVnInt()

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!******************************************************************************
    END SUBROUTINE VSSTUP
!******************************************************************************



!******************************************************************************
!  SUBROUTINE VSSOL                                                           *
!                                                                             *
!  Note, the following time stepping coefficients are computed in             *
!     SUBROUTINE READ_INPUT_3DVS and passed in GLOBAL_3DVS.                   *
!                                                                             *
!  IDTAlp1      = I*DelT*Alp1        - weights coriolis term in LHS matrix    *
!  IDT1MAlp1    = I*DelT*(1.-Alp1)   - weights coriolis term in RHS forcing   *
!  DTAlp3       = DelT*Alp3          - weights vert diff term in LHS matrix   *
!  DT1MAlp3     = DelT*(1-Alp3)      - weights vert diff term in RHS forcing  *
!  DTAlp2       = DelT*Alp2          - weights bot stress term in LHS matrix  *
!  DT1MAlp2     = DelT*(1.-Alp2)     - weights bot stress term in RHS forcing *
!                                                                             *
!  Note, the following are some of the global variables are used in this      *
!  subroutine and passed in GLOBAL_3DVS                                       *
!                                                                             *
!  q(MNP,NFEN) - 3D Complex Velocity field (GAMMA) from past time step.       *
!  UU, VV - depth-averaged horizontal velocity                                *
!  DAFluxX, DAFluxY - depth-averaged horizontal flux                          *
!                                                                             *
!  NP - number of nodes in horizontal grid                                    *
!  NFEN - number of nodes in the vertical grid                                *
!  BTP - total barotropic pressure (atmos press, water level, tidal potential)*
!                 at time levels s+1/2                                        *
!******************************************************************************
!                                                                             *
!   THIS  PROGRAM MAPS THE TRUE Z-COORDINATE SYSTEM INTO A DIMENSIONLESS      *
!   VERTICAL COORDINATE FROM [b,a] (BOTTOM TO TOP).  VALUES OF a=1 AND b=-1   *
!   ARE SET IN THE CODE.                                                      *
!                                                                             *
!******************************************************************************

    SUBROUTINE VSSOL(IT,TimeLoc)
!   kmd48.33bc added variable declarations for baroclinic boundary conditions
! RJW merged 08282008 Casey 071219: Added the following variable declarations from GLOBAL.
    USE GLOBAL, ONLY : NODECODE, NOFF, &
    RES_BC_FLAG, SBCTIME1, SBCTIME2, TBCTIME1, TBCTIME2, &
    SBCTIMEINC, TBCTIMEINC, SBCRATIO, TBCRATIO, sponge, &
    BCFLAG_TEMP, BCRivRATIO, RIVBCTIMINC, RIVBCTIME1, RIVBCTIME2
    USE GLOBAL_3DVS
    USE WRITE_OUTPUT, ONLY : writeOutput3D
    USE MESH, ONLY : NP, DP, NM, X, Y, NODELE, AREAS, LBArray_Pointer, &
    NEITABELE, NEITAB, NNEIGH, SFAC
    USE BOUNDARIES, ONLY : NOPE, NETA, NBD, LBCODEI, &
    BndBCRiver, totalbcrivernodes, SIII, CSII
    USE NodalAttributes, ONLY: EVM, &
    ManningsN, LoadManningsN, &
! Corbitt 120328: Local Advection
    AdvectionState, LoadAdvectionState, &
    Z0b_var, LoadZ0B_var
#ifdef CMPI
    USE MESSENGER
#endif
    IMPLICIT NONE

! RJW merged 08282008 Casey 071219: Added the following local variable declarations.
    INTEGER :: NCELE
    INTEGER :: TEMPNCELE
    INTEGER :: TEMPSTOP

    INTEGER :: IT
    INTEGER :: NEle           !local value of NetTabEle or element number
    INTEGER :: k              !vertical node loop counter (1-bottom, NFEN-surf)
    INTEGER :: NH             !horizontal node loop counter
    INTEGER :: N              !neighbor node loop counter
    INTEGER :: N1,N2,N3,NNFirst !local node numbers used to compute gradients, interpolate
    INTEGER :: LBP            !value of LBArray_Pointer at present horizontal node
! jgfdebug      INTEGER :: NN             !output loop counter

    REAL(SZ) :: KSlip         !equavalent linear slip coeff
    REAL(SZ) :: Z0B1          !Bottom roughness length (from mannings, const or from Z0Bvar)
    REAL(SZ) :: WSXsNH,WSYsNH !Wind stress components at time level s at node NH

    REAL(SZ) :: WSigma(MNFEN) !"sigma" vertical velocity
    REAL(SZ) :: Wf            !weighting coefficient in adjoint correction to w
    REAL(SZ) :: WfOHH         !Wf/(H time level s+1)^2
    REAL(SZ) :: WZSurfBC      !surface boundary condition value of w
    REAL(SZ) :: WZSurf        !computed value of w at surface
    REAL(SZ) :: WZCorrection  !adjoint correction compute for w

!   kmd48.33bc add in spherical factors
    REAL(SZ) :: SFacAvg
! jgfdebug      REAL(SZ) :: Zk            !z depth of any node k in the vertical
! jgfdebug      REAL(SZ) :: DelSig        ! sigma(k+1)-sigma(k)
    REAL(SZ) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
    REAL(SZ) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
    REAL(SZ) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
    REAL(SZ) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
! jgfdebug      REAL(SZ) :: SigmaNN       !Sigma value of a neighbor node

    REAL(SZ) :: VelNorm,VelTan !-QNormsp1(NH)/Hsp1 at flux boundary node
    REAL(SZ) :: CLBP,SLBP     !local values of CSII, SIII at boundary node LBP
    REAL(SZ) :: Auv1km1,Auv2km1 !initial real,imaginary parts of Mkm1 at flux boundary node
    REAL(SZ) :: Auv1km1star,Auv2km1star !rotated real,imaginary parts of Mkm1 at flux boundary node
! jgfdebug      REAL(SZ) :: Auv1k1,Auv2k1 !initial real,imaginary parts of Mk at flux boundary node
! jgfdebug      REAL(SZ) :: Auv1k1star,Auv2k1star !rotated real,imaginary parts of Mk at flux boundary node
    REAL(SZ) :: Auv1kp1,Auv2kp1 !initial real,imaginary parts of Mkp1 at flux boundary node
    REAL(SZ) :: Auv1kp1star,Auv2kp1star !rotated real,imaginary parts of Mkp1 at flux boundary node

    REAL(SZ) :: EtaN1,EtaN2,EtaN3,EtaNFirst !nodal values of IFNLFA(Eta1+Eta2)/2
    REAL(SZ) :: hN1,hN2,hN3,hNFirst !nodal values of DP
    REAL(SZ) :: DUDX(MNFEN),DVDY(MNFEN) !horizontal derivatives of velocity used to compute w
    REAL(SZ) :: Un,Vn         !real,imaginary components of qn
    REAL(SZ) :: DelU,DelV     !real, imaginary parts of q(k)-q(k-1)

    REAL(SZ) :: BTPN1,BTPN2,BTPN3,BTPNFirst !nodal values of BTP
    REAL(SZ) :: BTPDX2A,BTPDY2A !(Horiz. grads of BTP)*2*Element Area

! jgfdebug      REAL(SZ) :: BCPN1,BCPN2,BCPN3,BCPNFirst !nodal values of BCP
! jgfdebug      REAL(SZ) :: BCPDX2A,BCPDY2A !(Horiz. grads of BCP)*2*Element Area
! jgfdebug      REAL(SZ) :: SigTAvg       !avg SigT between 2 vertical nodes
! jgfdebug      REAL(SZ) :: HGORhoOAMB    !depth*gravity/(reference density)/(a-b)

! jgfdebug      REAL(SZ) :: DBCPDX2A
! jgfdebug      REAL(SZ) :: DBCPDY2A
    REAL(SZ) :: RCL
    REAL(SZ) :: Auv1k
    REAL(SZ) :: Auv2k
    REAL(SZ) :: Auv1kstar
    REAL(SZ) :: Auv2kstar

    REAL(8) :: KVnm(MNFEN,3)  !integral used in vertical stress term
    REAL(8) :: TimeLoc           !model time at time level s+1
    REAL(8) :: DEtaDT         !time derivative of water surface elev
    REAL(8) :: DEtaDX,DEtaDY  !horizontal derivatives of water surface elev at time level s
    REAL(8) :: DEtaDX2A,DEtaDY2A !(DEtaDX,DEtaDY)*2*Element Area
    REAL(8) :: DhDX,DhDY      !horizontal derivatives of DP
    REAL(8) :: DhDX2A,DhDY2A  !(DhDX,DhDY)*2*Element Area
    REAL(8) :: TotalArea2     !2*Area of all elements around a node
! jgfdebug      REAL(8) :: TotalBCPGArea2 !2*Area of all elements around a node used to compute the BCPG
    REAL(8) :: a1,a2,a3,b1,b2,b3 !elemental coefficients used in horizontal FE method

    REAL(8) :: Hs             !Total water depth at time level s
! jgfdebug      REAL(8) :: HsN2           !Total water depth at time level s at local node N2
    REAL(8) :: HsOAMB         !Hs/(a-b)
    REAL(8) :: HsHsOAMBAMB    !(Hs/(a-b))^2
    REAL(8) :: Hsp1           !Total water depth at time level s+1
    REAL(8) :: Hsp1OAMB       !Hsp1/(a-b)
    REAL(8) :: Hsp1Hsp1OAMBAMB !(Hsp1/(a-b))^2

    COMPLEX(SZ) :: Fr(MNFEN)      !right side forcing vector
    COMPLEX(SZ) :: Frstar         !rotated right side forcing vector at flux boundary node
    COMPLEX(SZ) :: Mkm1(MNFEN)    !1st column (k-1) in left side compact storage matrix
    COMPLEX(SZ) :: Mk(MNFEN)      !2nd column (k) in left side compact storage matrix
    COMPLEX(SZ) :: Mkp1(MNFEN)    !3rd column (k+1) in left side cpmpact storage matrix
    COMPLEX(SZ) :: LAdvec(MNFEN)  !lateral advection term in momentum eqn
    COMPLEX(SZ) :: LStress(MNFEN) !lateral stress term in momentum eqn
    COMPLEX(SZ) :: VAdvec(MNFEN)  !vertical advection term in momentum eqn
    COMPLEX(SZ) :: VStress(MNFEN) !vertical stress term in momentum eqn
! jgfdebug      COMPLEX(SZ) :: BCPG(MNFEN)    !baroclinic pressure gradient
    COMPLEX(SZ) :: BTPG           !total barotropic pressure gradient (incl TP & water level)
    COMPLEX(SZ) :: CCR,CCL        !coeffs used on right,left side of momentum eqn
    COMPLEX(SZ) :: VIVel          !vertically integrated velocity
    COMPLEX(SZ) :: DUDS           !complex vertical velocity gradient between bottom two nodes
    COMPLEX(SZ) :: qn,qN1,qN2,qN3,qNFirst !nodal values of q
! jgfdebug      COMPLEX(SZ) :: UnDqDX,VnDqDY  !derivatives used in lateral advection
    COMPLEX(SZ) :: UnDqDX2A,VnDqDY2A !UnDqDX,UnDqDY)*2*Element Area
    COMPLEX(SZ) :: DqDXDPhiDX2A,DqDYDPhiDY2A !derivatives used in lateral stress calc.
    COMPLEX(SZ) :: DqDSigmakm1,DqDSigmakp1 !vertical deriv. of q from k-1,k and k,k+1
    COMPLEX(SZ) :: DqDX2A(MNFEN)  !horizontal derivatives of complex
    COMPLEX(SZ) :: DqDY2A(MNFEN)  !velocity used in w calc

!   kmd48.33bc added the following variables for boundary conditions
    INTEGER :: J, NumofBCNode
    INTEGER :: NOD
    CHARACTER(80) :: CDUM80

    call setMessageSource("vssol")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!***************************************************************************
!     Set up several variables that are needed for 3D run


!     Zero out surface stress and stress forcing at top for case of no wind

    WSXsNH=0.d0
    WSYsNH=0.d0



!     If a baroclinic run, compute the 3D baroclinic pressure field
!     The buoyancy field is defined as
!     BCP(z)    =(gravity/rho ref)*          integral (SigT) from surface down to z
!     BCP(sigma)=(gravity/rho ref)*(H/(a-b))*integral (SigT) from a down to sigma
!     where
!     SigT = Sigma T = Rho - 1000 = density - 1000
!     SigT0 = Sigma t value of reference density (typically = 0)
!     Sigma = dimensionless vertical coordinate

!!      IF(CBaroclinic) THEN
!!        DO NH=1,NP             !loop over horizontal nodes
!!            Hs=DP(NH)+IFNLFA*Eta1(NH) !total depth at previous (s) timestep
!!            HGORhoOAMB=GORhoOAMB*Hs !(gravity/rho ref)*(H/(a-b))
!!            BCP(NH,NFEN)=0.d0
!!            DO k=NFEN-1,1,-1    !loop over vertical nodes, starting at top and working down
!!               SigTAvg=(SigT(NH,k+1)+SigT(NH,k))/2.d0
!!               DelSig=Sigma(k+1)-Sigma(k)
!!               BCP(NH,k)=BCP(NH,k+1)+HGORhoOAMB*(SigTAvg-SigT0)*DelSig
!!            ENDDO
!!         ENDDO
#ifdef CMPI
!     Update BCP on ghost nodes
!      CALL UPDATER3D(BCP)  !!!!Don't know if this is needed at this time
#endif
!!      ENDIF


!***************************************************************************
!     Compute 3D horizontal velocities

!     Loop over each horizontal node to compute the horizontal velocity

    DO NH=1,NP                !loop over horizontal nodes

    ! Corbitt 120328: Applies Advection Locally
        IF (LoadAdvectionState) Then
            IF (DP(NH) >= AdvectionState(NH)) THEN
                IFNLCT = IFNLCTE
            ELSE
                IFNLCT = 0
            ENDIF
        ENDIF

    !     Set up some values at the node being worked on

        Hs  = DP(NH)+IFNLFA*Eta1(NH) !Total depth at previous (s) timestep
        HsOAMB=Hs/AMB
        HsHsOAMBAMB=HsOAMB*HsOAMB
        Hsp1= DP(NH)+IFNLFA*Eta2(NH) !Total depth at present (s+1) timestep
        Hsp1OAMB=Hsp1/AMB
        Hsp1Hsp1OAMBAMB=Hsp1OAMB*Hsp1OAMB

        IF(NWS /= 0) THEN      !wind stress
            WSXsNH=WSX1(NH)
            WSYsNH=WSY1(NH)
        ENDIF

    ! RJW merged 08282008 Casey 071219: Solve for TEMPNCELE, which is the sum of the NCELE values
    !     for the elements attached to the current horizontal node.  Note that,
    !     for some reason, the k neighbor elements are not contained in the first
    !     k registers of NEITABELE.  (For instance, if a node is connected to four
    !     elements, those element numbers might be contained in registers
    !     1, 2, 4, and 5.)  So we have to be cute about pulling element numbers
    !     from NEITABELE.
    
        TEMPNCELE = 0
        TEMPSTOP = 0
        DO K=1,NODELE(NH)
            IF(NEITABELE(NH,K) == 0)THEN
                TEMPSTOP = TEMPSTOP + 1
            ENDIF
            TEMPNCELE = TEMPNCELE + &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),1))* &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),2))* &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),3))* &
            NOFF(NEITABELE(NH,K+TEMPSTOP))
        ENDDO

    ! Casey 071219: If TEMPNCELE is zero, then the current horizontal node is
    !     not attached to any active elements.  Thus, the horizontal velocities
    !     must be zero, and we can set up the matrix immediately and skip to the
    !     tri-diagonal solver.
    
        IF(TEMPNCELE == 0)THEN
            DO K=1,NFEN
                MKM1(K) = (0.D0,0.D0)
                MK(K)   = (1.D0,-1.D0)
                MKP1(K) = (0.D0,0.D0)
                FR(K)   = (0.D0,0.D0)
            ENDDO
            GOTO 999
        ENDIF
    !.....
    !     If specified normal flow boundary node with no tangential slip is
    !     lateral boundary condition, set up the matrix immediately and skip
    !     to the solution

        LBP=LBArray_Pointer(NH)
        IF(LBP > 0) THEN
        ! ssential normal flow
            IF((LBCodeI(LBP) >= 10) .AND. (LBCodeI(LBP) <= 19)) THEN
                VelNorm=-QNormsp1(LBP)/Hsp1 !with essential no tan.
                VelTan=0.D0
                SLBP=SIII(LBP)
                CLBP=CSII(LBP)
                DO k=1,NFEN
                    Mkm1(k)=(0.D0,0.D0)
                    Mk(k)=SLBP+iy*CLBP
                    Mkp1(k)=(0.D0,0.D0)
                    Fr(k)=VelTan+iy*VelNorm
                !                  EVTot(k)=0.D0 !f77diff??
                ENDDO
                GOTO 999
            ENDIF
        ENDIF

    ! Determine the bottom roughness length either from fort.15, from Manning's n
    ! or as read in from nodal attributes
        IF(LoadZ0B_var) THEN
            Z0B1 = Z0B_var(NH)
        ELSEIF (LoadManningsN) THEN
            Z0B1 = ( DP(NH)+IFNLFA*ETA2(NH) )* exp(-(1.0D0+ &
            ( (0.41D0*( DP(NH)+IFNLFA*ETA2(NH))**(1.0D0/6.0D0) )/ &
            (ManningsN(NH)*sqrt(g)) ) ))
        ELSE
            Z0B1 = Z0B
        ENDIF

    !     Compute the vertical eddy viscosity
    !... RJW changed Hs to Hsp1
    ! Hs is the previous time step (if node is recently wet then this is no good
    ! Hsp1 is current time step, if node is newly wet then we now have agreement
    ! and Hsp1 has a value that can be used
        CALL EDDYVIS(Hsp1,UU(NH),VV(NH), &
        WSXsNH,WSYsNH,BSX(NH),BSY(NH),NH,IT,Z0B1)

    !     Compute the integral KVnm

        KVnm(1,1)=0.d0
        KVnm(1,3)=-0.5d0*(EVTot(2)+EVTot(1))/(Sigma(2)-Sigma(1))
        KVnm(1,2)=-(KVnm(1,1)+KVnm(1,3))
        DO k=2,NFEN-1
            KVnm(k,1)=KVnm(k-1,3)
            KVnm(k,3)=-0.5d0*(EVTot(k+1)+EVTot(k))/(Sigma(k+1)-Sigma(k))
            KVnm(k,2)=-(KVnm(k,1)+KVnm(k,3))
        ENDDO
        KVnm(NFEN,1)=KVnm(NFEN-1,3)
        KVnm(NFEN,3)=0.d0
        KVnm(NFEN,2)=-(KVnm(NFEN,1)+KVnm(NFEN,3))

    !     Compute time derivative of water surface position

        DEtaDT=(Eta2(NH)-Eta1(NH))/DelT

    !     Start computing horizontal derivatives of water level, bathymetric
    !     depth and total barotropic pressure (atmos pres, water level,
    !     tidal potential) Note: TotalArea2 = 2X total elemental area
    !     surrounding a node

        DEtaDX=0.d0
        DEtaDY=0.d0
        DhDX=0.d0
        DhDY=0.d0
        BTPG=(0.d0,0.d0)
        DEtaDX2A=0.d0
        DhDX2A=0.d0
        DEtaDY2A=0.d0
        DhDY2A=0.d0
        BTPDX2A=0.d0
        BTPDY2A=0.d0
        TotalArea2=0.d0

        N1=NH
        EtaN1=IFNLFA*(Eta1(N1)+Eta2(N1))/2.d0
        hN1=DP(N1)
        BTPN1=BTP(N1)

        N2=NeiTab(NH,2)        !operate on 1st neighbor
        EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
        hN2=DP(N2)
        BTPN2=BTP(N2)

        NNFirst=N2             !save these values until end
        EtaNFirst=EtaN2
        hNFirst=hN2
        BTPNFirst=BTPN2

        DO N=3,NNeigh(NH)      !operate on rest of neighbors
            N3=N2               !shift previously computed values
            hN3=hN2             !shift previously computed values
            EtaN3=EtaN2
            BTPN3=BTPN2
            N2=NeiTab(NH,N)     !select new neighbor to work on
            EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
            hN2=DP(N2)
            BTPN2=BTP(N2)
            NEle=NeiTabEle(NH,N-2) !element# defined by nodes NH,NN2,NN1
        ! RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
        ! jgf48.50 If NEle is zero, the NOFF check would be out of
        ! bounds, making it look like NOFF is always zero (depending
        ! on how the compiler handles array out-of-bounds)
            IF (NEle == 0) THEN
                CYCLE
            ENDIF
        ! RJW 06092009  modified to elimiate zero index NOFF
            IF((NEle /= 0)) THEN  !if element is active, compute velocity grads
                NCELE = NODECODE(NH)*NODECODE(N2)* &
                NODECODE(N3)*NOFF(NELE)
                IF((NCELE /= 0)) THEN  !if element is active, compute velocity grads

                    TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
                !               kmd48.33bc added in spherical factors for 3D
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
                    BTPDX2A=BTPDX2A+(BTPN1*b1+BTPN2*b2+BTPN3*b3)
                    BTPDY2A=BTPDY2A+(BTPN1*a1+BTPN2*a2+BTPN3*a3)
                ENDIF
            ENDIF
        ENDDO

        N3=N2                  !wrap back to beginning to get final contribution
        hN3=hN2
        EtaN3=EtaN2
        BTPN3=BTPN2
        N2=NNFirst
        hN2=hNFirst
        EtaN2=EtaNFirst
        BTPN2=BTPNFirst
        NEle=NeiTabEle(NH,NNeigh(NH)-1)
    !. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
    ! jgf48.50 Avoid array out of bounds
        IF((NEle /= 0)) THEN  !if element is active, compute velocity grads
            NCELE = NODECODE(NH)*NODECODE(N2)* &
            NODECODE(N3)*NOFF(NELE)
            IF((NCELE /= 0)) THEN  !if element is active, compute velocity grads
                TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
            !    kmd48.33bc added in spherical factors for 3D
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
                BTPDX2A=BTPDX2A+(BTPN1*b1+BTPN2*b2+BTPN3*b3)
                BTPDY2A=BTPDY2A+(BTPN1*a1+BTPN2*a2+BTPN3*a3)
            ENDIF
        ENDIF

        IF(TotalArea2 /= 0.) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
            BTPG=(BTPDX2A+iy*BTPDY2A)/TotalArea2
        ENDIF

    !     Finished computing horizontal derivatives of water level,
    !     bathymetric depth and total barotropic pressure (atmos pres, water
    !     level, tidal potential)

    !     Compute the "sigma" vertical velocity from the "z" vertical velocity

        DO k=1,NFEN
            SigmaMAOAMB=(Sigma(k)-A)/AMB
            SigmaMBOAMB=(Sigma(k)-B)/AMB
            WSigma(k) = WZ(NH,k) - SigmaMBOAMB*DEtaDT &
            - REAL(q(NH,k))*(SigmaMBOAMB*DEtaDX+SigmaMAOAMB*DhDX) &
            - AIMAG(q(NH,k))*(SigmaMBOAMB*DEtaDY+SigmaMAOAMB*DhDY)
        ENDDO


    !     Start computing advection and stress terms in the momentum
    !     equation at each level in the vertical

        DO k=1,NFEN

        !     Compute the vertical advection and vertical stress terms

            IF(k == 1) THEN
                DqDSigmakp1=(q(NH,k+1)-q(NH,k))/(Sigma(k+1)-Sigma(k))
                VAdvec(k)=DqDsigmakp1* &
                (2.d0*WSigma(k)+WSigma(k+1))*Inm(k,3)/HsOAMB
                VStress(k)=(q(NH,k)*KVnm(k,2)+q(NH,k+1)*KVnm(k,3)) &
                /HsHsOAMBAMB
            ENDIF
            IF((k > 1) .AND. (k < NFEN)) THEN
                DqDSigmakm1=DqDSigmakp1
                DqDSigmakp1=(q(NH,k+1)-q(NH,k))/(Sigma(k+1)-Sigma(k))
                VAdvec(k)= &
                ( DqDSigmakm1*(WSigma(k-1)+2.d0*WSigma(k))*Inm(k,1) &
                +DqDSigmakp1*(2.d0*WSigma(k)+WSigma(k+1))*Inm(k,3) ) &
                /HsOAMB
                VStress(k)=(q(NH,k-1)*KVnm(k,1)+q(NH,k)*KVnm(k,2) &
                +q(NH,k+1)*KVnm(k,3))/HsHsOAMBAMB
            ENDIF
            IF(k == NFEN) THEN
                DqDSigmakm1=DqDSigmakp1
                VAdvec(k)= &
                DqDSigmakm1*(WSigma(k-1)+2.d0*WSigma(k))*Inm(k,1) &
                /HsOAMB
                VStress(k)=(q(NH,k-1)*KVnm(k,1)+q(NH,k)*KVnm(k,2)) &
                /HsHsOAMBAMB
            ENDIF
            VAdvec(k)=sponge(NH)*IFNLCT*VAdvec(k)

        !     Compute lateral advection and lateral stress terms

            UnDqDX2A=0.d0
            VnDqDY2A=0.d0
            DqDXDPhiDX2A=0.d0
            DqDYDPhiDY2A=0.d0

            N1=NH               !node 1 is always the central node
            qN1=q(N1,k)

            N2=NEITAB(NH,2)     !operate on 1st neighbor
            qN2=q(N2,k)

            NNFirst=N2          !save these values until end
            qNFirst=qN2

            DO N=3,NNEIGH(NH)   !operate on rest of neighbors
                N3=N2            !shift previously computed values
                qN3=qN2
                N2=NEITAB(NH,N)  !select new neighbor to work on
                qN2=q(N2,k)
                NEle=NeiTabEle(NH,N-2) !element# defined by nodes N1,N2,N3
            !. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
            ! jgf48.50 Avoid array out of bounds
                IF (NEle == 0) THEN
                    CYCLE
                ENDIF
                IF((NEle /= 0)) THEN  !if element is active, compute velocity grads
                    NCELE = NODECODE(NH)*NODECODE(N2)* &
                    NODECODE(N3)*NOFF(NELE)
                    IF (NCELE /= 0) THEN !if element exists, compute terms
                    !    kmd48.33bc added in spherical factors for 3D
                        SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                        qn=(qN1+qN2+qN3)/3.d0
                        Un=REAL(qn)
                        Vn=AIMAG(qn)
                        a1=X(N3)-X(N2)
                        a2=X(N1)-X(N3)
                        a3=X(N2)-X(N1)
                        b1=(Y(N2)-Y(N3))*SFacAvg
                        b2=(Y(N3)-Y(N1))*SFacAvg
                        b3=(Y(N1)-Y(N2))*SFacAvg
                        UnDqDX2A=UnDqDX2A+Un*(qN1*b1+qN2*b2+qN3*b3)
                        VnDqDY2A=VnDqDY2A+Vn*(qN1*a1+qN2*a2+qN3*a3)
                        DqDXDPhiDX2A=DqDXDPhiDX2A &
                        +(qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
                        DqDYDPhiDY2A=DqDYDPhiDY2A &
                        +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
                    ENDIF
                ENDIF
            ENDDO

            N3=N2               !wrap back to beginning to get final contribution
            qN3=qN2
            N2=NNFIRST
            qN2=qNFirst
            NEle=NeiTabEle(NH,NNeigh(NH)-1)
        !. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
        ! RJW 06092009  modified to elimiate zero index NOFF
            IF((NEle /= 0)) THEN  !if element is active, compute velocity grads
                NCELE = NODECODE(NH)*NODECODE(N2)* &
                NODECODE(N3)*NOFF(NELE)
                IF((NCELE /= 0)) THEN  !if element is active, compute velocity grads
                !    kmd48.33bc added in spherical factors for 3D
                    SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                    qn=(qN1+qN2+qN3)/3.d0
                    Un=real(qn)
                    Vn=aimag(qn)
                    a1=X(N3)-X(N2)
                    a2=X(N1)-X(N3)
                    a3=X(N2)-X(N1)
                    b1=(Y(N2)-Y(N3))*SFacAvg
                    b2=(Y(N3)-Y(N1))*SFacAvg
                    b3=(Y(N1)-Y(N2))*SFacAvg
                    UnDqDX2A=UnDqDX2A+Un*(qN1*b1+qN2*b2+qN3*b3)
                    VnDqDY2A=VnDqDY2A+Vn*(qN1*a1+qN2*a2+qN3*a3)
                    DqDXDPhiDX2A=DqDXDPhiDX2A &
                    +(qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
                    DqDYDPhiDY2A=DqDYDPhiDY2A &
                    +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
                ENDIF
            ENDIF

            IF(TotalArea2 == 0.) THEN
                LAdvec(k)=(0.d0,0.d0)
                LStress(k)=(0.d0,0.d0)
            ELSE
                LAdvec(k)=sponge(NH)*IFNLCT*(UnDqDX2A+VnDqDY2A)/TotalArea2
                LStress(k)=3.d0*EVM(NH)*(DqDXDPhiDX2A+DqDYDPhiDY2A) &
                /TotalArea2
            ENDIF

        ENDDO

    !     Finished computing advection and stress terms in the momentum
    !     equation at each level in the vertical

    !     Compute the equivalent linear slip coefficient
    ! RJW 02252009 changed KP for ISlip =2
    !     For ISlip = 1, KP is the equivalent linear slip coefficient
    !     For ISlip = 2, KP is the minimum equivalent quadratic slip coefficient
    !     For ISlip = 3, KP is the quadratic slip coefficient
        IF(ISlip == 1) KSlip=KP
        IF(ISlip == 2) THEN
        ! RJW merged 08282008 Casey 071219: Let's get funky.
        ! RJW enable this modification
        !     and add the following to base the drag coefficient on the actual
        !     physics of the boundary layer.
        ! ---

            IF (ABS( ( ( SIGMA(2)-SIGMA(1) )/2 ) &
            *(DP(NH)+IFNLFA*ETA2(NH)) ) &
             < 1.D0 )THEN


                KSlip = (1.D0 / ( (1.D0/0.41D0) * &
                LOG( (ABS( ( ( SIGMA(2)-SIGMA(1) )/2 ) &
                *( DP(NH)+IFNLFA*ETA2(NH) )) + Z0B1 ) &
                / Z0B1 ) ) )**2.D0 &
                * ABS(Q(NH,1))
            ELSE
                KSlip = (1.D0 / ( (1.D0/0.41D0) * &
                LOG( (1.D0+Z0B1) &
                / (Z0B1) ) ) )**2.D0 &
                * ABS(Q(NH,1))
            ENDIF
            IF(KSlip > 1.0d0* ABS(Q(NH,1))) KSlip=1.0d0*ABS(Q(NH,1))
        !          IF(KSlip.LT.0.0010d0*ABS(Q(NH,1))) KSlip=0.0010d0*ABS(Q(NH,1))
            IF(KSlip < KP*ABS(Q(NH,1))) KSlip=KP*ABS(Q(NH,1))

        ENDIF
        IF(ISlip == 3) THEN
            KSlip=KP*ABS(Q(NH,1))
            IF(KSlip > 0.1d0) KSlip=0.1d0
            IF(KSlip < 0.0d0) KSlip=0.0d0
        ENDIF
    !     Set up the RHS forcing vector Fr and LHS matrix in compact storage
    !     (Mkm1,Mk,Mkp1)

        CCL = 1.d0+Corif(NH)*IDTAlp1
        CCR = 1.d0-Corif(NH)*IDT1MAlp1
        RCL = DTAlp3/Hsp1Hsp1OAMBAMB

        IF(ISlip == 0) THEN    !no slip bottom boundary condition
            Fr(1)   = (0.d0,0.d0)
            Mkm1(1) = (0.d0,0.d0)
            Mk(1)   = (1.d0,-1.d0) !note: -I*IV=V
            Mkp1(1) = (0.d0,0.d0)
        ELSE                   ! slip bottom boundary condition
            Fr(1) = (CCR*q(NH,1) &
            -DelT*(LAdvec(1)  +LStress(1)  +BPG(NH,1)))*Inm(1,2) &
            + (CCR*q(NH,2) &
            -DelT*(LAdvec(2)  +LStress(2)  +BPG(NH,2)))*Inm(1,3) &
            - DelT*(VAdvec(1)+BTPG*LVn(1))-DT1MAlp3*VStress(1) &
            - q(NH,1)*DT1MAlp2*KSlip/HsOAMB
            Mkm1(1) = (0.d0,0.d0)
            Mk(1)   = CCL*Inm(1,2) + RCL*KVnm(1,2) &
            + DTAlp2*KSlip/Hsp1OAMB
            Mkp1(1) = CCL*Inm(1,3) + RCL*KVnm(1,3)
        ENDIF

        DO k=2,NFEN-1
            Fr(k) = (CCR*q(NH,k-1) &
            -DelT*(LAdvec(k-1)+LStress(k-1)+BPG(NH,k-1)))*Inm(k,1) &
            + (CCR*q(NH,k) &
            -DelT*(LAdvec(k)  +LStress(k)  +BPG(NH,k)  ))*Inm(k,2) &
            + (CCR*q(NH,k+1) &
            -DelT*(LAdvec(k+1)+LStress(k+1)+BPG(NH,k+1)))*Inm(k,3) &
            - DelT*(VAdvec(k)+BTPG*LVn(k))-DT1MAlp3*VStress(k)
            Mkm1(k) = CCL*Inm(k,1) + RCL*KVnm(k,1)
            Mk(k)   = CCL*Inm(k,2) + RCL*KVnm(k,2)
            Mkp1(k) = CCL*Inm(k,3) + RCL*KVnm(k,3)
        ENDDO
    ! RJW 09/01/2009 RJW k = NFEN-1
    !            try to fix this by setting k=NFEN
    !         k=NFEN
        Fr(NFEN) = (CCR*q(NH,k-1) &
        -DelT*(LAdvec(k-1)+LStress(k-1)+BPG(NH,k-1)))*Inm(k,1) &
        + (CCR*q(NH,k) &
        -DelT*(LAdvec(k)  +LStress(k)  +BPG(NH,k)  ))*Inm(k,2) &
        - DelT*(VAdvec(k)+BTPG*LVn(k))-DT1MAlp3*VStress(k) &
        + DelT*0.5d0*((WSX2(NH)+iy*WSY2(NH))/Hsp1OAMB &
        +(WSX1(NH)+iy*WSY1(NH))/HsOAMB)
        Mkm1(NFEN) = CCL*Inm(NFEN,1) + RCL*KVnm(NFEN,1)
        Mk(NFEN)   = CCL*Inm(NFEN,2) + RCL*KVnm(NFEN,2)
        Mkp1(NFEN) = (0.d0,0.d0)

    !     Start section to modify equations depending on normal flux
    !     boundary condition
    !     0 <= LBcodeI <= 10, essential normal flux and free tangential slip
    !     this b.c. is taken care of in the code section below
    !     10 <= LBcodeI <= 19, essential normal flux and zero tangential slip
    !     this b.c. is taken care of above
    !     20 <= LBcodeI <= 29, natural normal flux and free tangential slip
    !     this b.c. requires on manipulation of momentum eqns.  Do nothing!

        LBP=LBArray_Pointer(NH)
        IF(LBP > 0) THEN      !flux boundary
            IF((LBCODEI(LBP) >= 0) .AND. (LBCODEI(LBP) <= 9)) THEN
                SLBP=SIII(LBP)
                CLBP=CSII(LBP)
                VelNorm=-QNormsp1(LBP)/Hsp1

                Mkm1(1)    =(0.d0,0.d0)
                Auv1k      =Real(Mk(1))
                Auv2k      =AImag(Mk(1))
                Auv1kstar  =Auv1k*SLBP
                Auv2kstar  =Auv1k*CLBP
                Mk(1)      =Auv1kstar+iy*Auv2kstar
                Auv1kp1    =Real(Mkp1(1))
                Auv2kp1    =AImag(Mkp1(1))
                Auv1kp1star=Auv1kp1*SLBP
                Auv2kp1star=Auv1kp1*CLBP
                Mkp1(1)    =Auv1kp1star+iy*Auv2kp1star
                Frstar     =Real(Fr(1))*SLBP-AImag(Fr(1))*CLBP &
                +  (Auv2k+Auv2kp1)*VelNorm &
                +iy*(Auv1k+Auv1kp1)*VelNorm
                Fr(1)      =Frstar
                DO k=2,NFEN-1
                    Auv1km1    =Real(Mkm1(k))
                    Auv2km1    =AImag(Mkm1(k))
                    Auv1km1star=Auv1km1*SLBP
                    Auv2km1star=Auv1km1*CLBP
                    Mkm1(k)    =Auv1km1star+iy*Auv2km1star
                    Auv1k      =Real(Mk(k))
                    Auv2k      =AImag(Mk(k))
                    Auv1kstar  =Auv1k*SLBP
                    Auv2kstar  =Auv1k*CLBP
                    Mk(k)      =Auv1kstar+iy*Auv2kstar
                    Auv1kp1    =Real(Mkp1(k))
                    Auv2kp1    =AImag(Mkp1(k))
                    Auv1kp1star=Auv1kp1*SLBP
                    Auv2kp1star=Auv1kp1*CLBP
                    Mkp1(k)    =Auv1kp1star+iy*Auv2kp1star
                    Frstar     =Real(Fr(k))*SLBP-AImag(Fr(k))*CLBP &
                    +  (Auv2km1+Auv2k+Auv2kp1)*VelNorm &
                    +iy*(Auv1km1+Auv1k+Auv1kp1)*VelNorm
                    Fr(k)      =Frstar
                ENDDO
                Auv1km1    =Real(Mkm1(NFEN))
                Auv2km1    =AImag(Mkm1(NFEN))
                Auv1km1star=Auv1km1*SLBP
                Auv2km1star=Auv1km1*CLBP
                Mkm1(NFEN) =Auv1km1star+iy*Auv2km1star
                Auv1k      =Real(Mk(NFEN))
                Auv2k      =AImag(Mk(NFEN))
                Auv1kstar  =Auv1k*SLBP
                Auv2kstar  =Auv1k*CLBP
                Mk(NFEN)   =Auv1kstar+iy*Auv2kstar
                Mkp1(NFEN) =(0.d0,0.d0)
                Frstar     =Real(Fr(NFEN))*SLBP-AImag(Fr(NFEN))*CLBP &
                +  (Auv2km1+Auv2k)*VelNorm &
                +iy*(Auv1km1+Auv1k)*VelNorm
                Fr(NFEN)   =Frstar
            ENDIF
        ENDIF

    !     Finished section to modify equations depending on normal flux
    !     boundary condition
    !     Decompose and solve the system


        999 CALL TRIDIAG(Mkm1,Mk,Mkp1,Fr,Gamma,NFEN)

    !     Changed this q value to qkp1 since we need to save q and qkp1 for
    !     transport calculations
    !     Save the horizontal velocity and Eddy Viscosity solutions
    !   kmd48.33bc added the vertical diffusion term for the transport equation

        DO k=1,NFEN
            qkp1(NH,k) = Gamma(k)
            EV(NH,k)=EVTot(k)
            DV(NH,k)=NTVTot(k)
        END DO

    ENDDO

!     Finished computing horizontal velocities


!****************************************************************************
!     Update horizontal velocity solution on ghost nodes

#ifdef CMPI
    CALL UPDATEC3D(Qkp1)
#endif


!****************************************************************************
!     Compute the depth averaged velocity, depth averaged flux, bottom stress and
!     velocity dispersion

    DO NH=1,NP                                !loop over horizontal nodes

        Hsp1= DP(NH)+IFNLFA*Eta2(NH)           !Total depth at present (s+1) timestep
        VIVel = (0.d0,0.d0)
        DO k=1,NFEN
            VIVel = VIVel + Qkp1(NH,k)*LVn(k)
        END DO
        VIVel = VIVel/amb
        UU(NH) = REAL(VIVel)
        VV(NH) = AIMAG(VIVel)

        DAFluxX(NH)=UU(NH)*Hsp1
        DAFluxY(NH)=VV(NH)*Hsp1

        IF(ISlip == 0) THEN
            DUDS = (Qkp1(NH,2)-Qkp1(NH,1))/(Sigma(2)-Sigma(1))
            BSX(NH) = EV(NH,1)*REAL(DUDS)
            BSY(NH) = EV(NH,1)*AIMAG(DUDS)
        ENDIF
        IF(ISlip /= 0) THEN
            BSX(NH) = KSlip*REAL(Qkp1(NH,1))
            BSY(NH) = KSlip*AIMAG(Qkp1(NH,1))
        ENDIF

        DUU(NH) = 0.d0
        DUV(NH) = 0.d0
        DVV(NH) = 0.d0

    ! Corbitt 120328: Applies Advection Locally
        IF (LoadAdvectionState) Then
            IF (DP(NH) >= AdvectionState(NH)) THEN
                IFNLCT = IFNLCTE
            ELSE
                IFNLCT = 0
            ENDIF
        ENDIF

        IF(IFNLCT == 1) THEN
            CALL VSDISP (IT,NH,Hsp1,UU(NH),VV(NH), &
            DUU(NH),DUV(NH),DVV(NH))
        ENDIF

    ENDDO

!***************************************************************************
!     Compute "z" vertical velocity


!     Loop over each horizontal node to compute the "z" version of the
!     vertical velocity

    DO NH=1,NP

    !. RJW merged 08282008 Casey 071219: Solve for TEMPNCELE, same as above.
    
        TEMPNCELE = 0
        TEMPSTOP = 0
        DO K=1,NODELE(NH)
            IF(NEITABELE(NH,K) == 0)THEN
                TEMPSTOP = TEMPSTOP + 1
            ENDIF
            TEMPNCELE = TEMPNCELE + &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),1))* &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),2))* &
            NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),3))* &
            NOFF(NEITABELE(NH,K+TEMPSTOP))
        ENDDO

        IF(TEMPNCELE == NODELE(NH))THEN
            TEMPNCELE = 1
        ELSEIF(TEMPNCELE == 0)THEN
            TEMPNCELE = 0
        ELSE
            TEMPNCELE = 1
        ENDIF

    !     Set up some values at the node being worked on

        Hsp1=DP(NH)+IFNLFA*Eta2(NH)
        Hsp1OAMB=Hsp1/AMB

    !     Compute time derivative of water surface position

        DEtaDT=(Eta2(NH)-Eta1(NH))/DelT

    !     Compute horizontal derivatives of water surface position,
    !     bathymetric depth and horizontal velocity Note: TotalArea2 = 2X
    !     total elemental area surrounding a node

        DO k=1,NFEN
            DqDX2A(k)=(0.d0,0.d0)
            DqDY2A(k)=(0.d0,0.d0)
            DUDX(k)=0.d0
            DVDY(k)=0.d0
        ENDDO
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
        EtaN1=IFNLFA*Eta2(N1)
        hN1=DP(N1)

        N2=NeiTab(NH,2)        !operate on 1st neighbor
        EtaN2=IFNLFA*Eta2(N2)
        hN2=DP(N2)

        NNFirst=N2             !save these values until end
        EtaNFirst=EtaN2
        hNFirst=hN2

        DO N=3,NNeigh(NH)      !operate on rest of neighbors
            N3=N2               !shift previously computed values
            hN3=hN2             !shift previously computed values
            EtaN3=EtaN2
            N2=NeiTab(NH,N)     !select new neighbor to work on
            EtaN2=IFNLFA*Eta2(N2)
            hN2=DP(N2)
            NEle=NeiTabEle(NH,N-2) !element# defined by nodes NH,NN2,NN1
        !. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of this IF statement.
        ! jgf48.50 Avoid array out of bounds
            IF (NEle == 0) THEN
                CYCLE
            ENDIF
            NCELE = NODECODE(NH)*NODECODE(N2)* &
            NODECODE(N3)*NOFF(NELE)
            IF(NCELE /= 0) THEN  !if element is active, compute velocity grads
                TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
            !    kmd48.33bc added in spherical factors for 3D
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
                DO k=1,NFEN
                    qN1=qkp1(N1,k)
                    qN2=qkp1(N2,k)
                    qN3=qkp1(N3,k)
                    DqDX2A(k)=DqDX2A(k)+(qN1*b1+qN2*b2+qN3*b3)
                    DqDY2A(k)=DqDY2A(k)+(qN1*a1+qN2*a2+qN3*a3)
                ENDDO
            ENDIF
        ENDDO

        N3=N2                  !wrap back to beginning to get final contribution
        hN3=hN2
        EtaN3=EtaN2
        N2=NNFirst
        hN2=hNFirst
        EtaN2=EtaNFirst
        NEle=NeiTabEle(NH,NNeigh(NH)-1)
    !. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of this IF statement.
    ! jgf48.50 Avoid array out-of-bounds
        IF (NEle /= 0 ) THEN
            NCELE = NODECODE(NH)*NODECODE(N2)* &
            NODECODE(N3)*NOFF(NELE)
            IF (NCELE /= 0) THEN  !if element is active, compute velocity grads
                TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
            !    kmd48.33bc added in spherical factors for 3D
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
                DO k=1,NFEN
                    qN1=qkp1(N1,k)
                    qN2=qkp1(N2,k)
                    qN3=qkp1(N3,k)
                    DqDX2A(k)=DqDX2A(k)+(qN1*b1+qN2*b2+qN3*b3)
                    DqDY2A(k)=DqDY2A(k)+(qN1*a1+qN2*a2+qN3*a3)
                ENDDO
            ENDIF
        ENDIF

        IF(TotalArea2 /= 0.) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
            DO k=1,NFEN
                DUDX(k)=REAL(DqDX2A(k))/TotalArea2
                DVDY(k)=AIMAG(DqDY2A(k))/TotalArea2
            ENDDO
        ENDIF

    !     Evaluate the "z" vertical velocity

    !. RJW merged 08282008 Casey 071219: Added an IF statement so that the vertical velocities are only computed
    !     and corrected at nodes that are connected to at least one active element.
    
        IF(TEMPNCELE /= 0)THEN
        !.
            WZkp1(NH,1)=-REAL(qkp1(NH,1))*DhDX-AIMAG(qkp1(NH,1))*DhDY

            DO k=2,NFEN
                DelSigO2=(Sigma(k)-Sigma(k-1))/2.d0
                SigAvgMAOAMB=((Sigma(k)+Sigma(k-1))/2.d0 - A)/AMB
                DelU=REAL(qkp1(NH,k)-qkp1(NH,k-1))
                DelV=AIMAG(qkp1(NH,k)-qkp1(NH,k-1))
                WZkp1(NH,k)=WZkp1(NH,k-1) &
                - DelSigO2*Hsp1OAMB*(DUDX(k)+DVDY(k)+DUDX(k-1)+DVDY(k-1)) &
                + (DEtaDX+SigAvgMAOAMB*(DhDX+DEtaDX))*DelU &
                + (DEtaDY+SigAvgMAOAMB*(DhDY+DEtaDY))*DelV
            ENDDO

        !     Correct this using Adjoint method

            Wf=0.d0                !This value should match surface B.C. exactly
            Hsp1= DP(NH)+IFNLFA*Eta2(NH) !Total depth at present (s+1) timestep
            WfOHH=Wf/(Hsp1*Hsp1)
            WZSurfBC=DEtaDT+REAL(qkp1(NH,NFEN))*DEtaDX &
            +AIMAG(qkp1(NH,NFEN))*DEtaDY
            WZSurf=WZkp1(NH,NFEN)

            DO k=1,NFEN
                WZCorrection=(WZSurfBC-WZSurf)*(WfOHH+(Sigma(k)-b)/AMB) &
                /(2.d0*WfOHH+1)
                WZkp1(NH,k)=WZkp1(NH,k)+WZCorrection
            ENDDO

        ! RJW merged 08282008 Casey 071219: Added an ELSEIF statement to zero out the vertical velocities at nodes
        !     that are not connected to any active elements.
        
        ELSEIF(TEMPNCELE == 0)THEN
            DO K=1,NFEN
                WZkp1(NH,K) = (0.D0,0.D0)
            ENDDO
        ENDIF
    !.
    
    !     End loop over horizontal nodes to compute vertical velocity
    
    ENDDO

!***************************************************************************
!     Update vertical velocity on ghost nodes - this is necessary

#ifdef CMPI
    CALL UPDATER3D(WZkp1)
#endif


!**************************************************************************
!     For a prognostic baroclinic run, compute new density, temperature, salinity fields

!   kmd48.33bc need to update the boundary condition information for both
!              salinity and temperature in prognostic runs. The boundary
!              conditions are read in similar to the aperiodic elevation
!              boundary conditions
    IF ((CBAROCLINIC) .AND. (RES_BC_FLAG > 0) .AND. (NOPE > 0)) THEN
        IF ((ABS(RES_BC_FLAG) == 2) .OR. (ABS(RES_BC_FLAG) == 4)) THEN
            IF(TimeLoc > SBCTIME2) THEN
                SBCTIME1=SBCTIME2
                SBCTIME2=SBCTIME1+SBCTIMEINC
                READ(36,'(A)') CDUM80
                DO NumofBCNode=1,NETA
                    DO J=1,NFEN
                        RESSAL1(NumofBCNode,J)=RESSAL2(NumofBCNode,J)
                    END DO
                    READ(36,*) NOD,(RESSAL2(NumofBCNode,J),J=1,NFEN)
                END DO
            END IF
            SBCRATIO=(TimeLoc-SBCTIME1)/SBCTIMEINC
            DO NumofBCNode=1,NETA
                DO J=1,NFEN
                    RESSAL(NumofBCNode,J)=RESSAL1(NumofBCNode,J)+ &
                    SBCRATIO*(RESSAL2(NumofBCNode,J)- &
                    RESSAL1(NumofBCNode,J))
                END DO
            END DO
        END IF
        IF ((ABS(RES_BC_FLAG) == 3) .OR. &
        (ABS(RES_BC_FLAG) == 4)) THEN
            IF(TimeLoc > TBCTIME2) THEN
                TBCTIME1=TBCTIME2
                TBCTIME2=TBCTIME1+TBCTIMEINC
                READ(37,'(A)') CDUM80
                DO NumofBCNode=1,NETA
                    DO J=1,NFEN
                        RESTEMP1(NumofBCNode,J)=RESTEMP2(NumofBCNode,J)
                    END DO
                    READ(37,*) NOD,(RESTEMP2(NumofBCNode,J),J=1,NFEN)
                END DO
            END IF
            TBCRATIO=(TimeLoc-TBCTIME1)/TBCTIMEINC
            DO NumofBCNode=1,NETA
                DO J=1,NFEN
                    RESTEMP(NumofBCNode,J)=RESTEMP1(NumofBCNode,J)+ &
                    TBCRATIO*(RESTEMP2(NumofBCNode,J)- &
                    RESTEMP1(NumofBCNode,J))
                END DO
            END DO
        END IF
    END IF

! md - adding information for rivers - baroclinic

    IF ((CBaroclinic) .AND. (BndBCRiver)) THEN
        IF (TIMELOC > RIVBCTIME2) THEN
            RIVBCTIME1=RIVBCTIME2
            RIVBCTIME2=RIVBCTIME1+RIVBCTIMINC
            DO J=1,totalbcrivernodes
                IF (IDEN == 2) THEN
                    DO K=1,NFEN
                        BCRivSalN1(J,K)=BCRivSalN2(J,K)
                    END DO
                    READ(39,*) NOD, (BCRivSalN2(J,K),K=1,NFEN)
                ELSE IF (IDEN == 3) THEN
                    DO K=1,NFEN
                        BCRivTempN1(J,K)=BCRivTempN2(J,K)
                    END DO
                    READ(39,*) NOD, (BCRivTempN2(J,K),K=1,NFEN)
                ELSE IF (IDEN == 4) THEN
                    DO K=1,NFEN
                        BCRivSalN1(J,K)=BCRivSalN2(J,K)
                        BCRivTempN1(J,K)=BCRivTempN2(J,K)
                    END DO
                    READ(39,*) NOD, (BCRivSalN2(J,K),BCRivTempN2(J,K),K=1,NFEN)
                END IF
            END DO
        END IF
        BCRivRATIO=(TIMELOC-RIVBCTIME1)/RIVBCTIMINC
        DO J=1,totalbcrivernodes
            IF (IDEN == 2) THEN
                DO K=1,NFEN
                    BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO* &
                    (BCRivSalN2(J,K)-BCRivSalN1(J,K))
                END DO
            ELSE IF (IDEN == 3) THEN
                DO K=1,NFEN
                    BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO* &
                    (BCRivTempN2(J,K)-BCRivTempN1(J,K))
                END DO
            ELSE IF (IDEN == 4) THEN
                DO K=1,NFEN
                    BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO* &
                    (BCRivSalN2(J,K)-BCRivSalN1(J,K))
                    BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO* &
                    (BCRivTempN2(J,K)-BCRivTempN1(J,K))
                END DO
            END IF
        END DO
    END IF

!   Start computing the new density, temperature and salinity values
    IF(C3D_BTrans) THEN
        IF(IDEN == 1) THEN
        !           CALL 3D_SIGMAT_TRANS()
        ELSE IF(IDEN == 2) THEN
            CALL TRANS_3D(SAL,NLSD,NVSD,SALkp1,IDEN,TimeLoc)
        ELSE IF(IDEN == 3) THEN
            CALL TRANS_3D(TEMP,NLTD,NVTD,TEMPkp1,IDEN,TimeLoc)
        ELSE IF(IDEN == 4) THEN
            CALL TRANS_3D(SAL,NLSD,NVSD,SALkp1,IDEN-2,TimeLoc)
            CALL TRANS_3D(TEMP,NLTD,NVTD,TEMPkp1,IDEN-1,TimeLoc)
        ENDIF
    
    !     Update new temperature and salinity fields on ghost nodes
    
#ifdef CMPI
    !      CALL UPDATER3D(SIGT) !jgf45.12
        IF((IDEN == 3) .OR. (IDEN == 4)) CALL UPDATER3D(Tempkp1)
        IF((IDEN == 2) .OR. (IDEN == 4)) CALL UPDATER3D(Salkp1)
#endif

    !     Need to update and move results from the future time levels to the
    !     present time levels so we can advance in time and before obtaining
    !     new sigma-t values. We will update velocities later.
        DO NH = 1,NP
            DO k = 1,NFEN
                IF ((IDEN == 2) .OR. (IDEN == -2)) THEN
                    SAL(NH,k)=SALkp1(NH,k)
                ELSE IF ((IDEN == 3) .OR. (IDEN == -3)) THEN
                    TEMP(NH,k)=TEMPkp1(NH,k)
                ELSE IF ((IDEN == 4) .OR. (IDEN == -4)) THEN
                    SAL(NH,k)=SALkp1(NH,k)
                    TEMP(NH,k)=TEMPkp1(NH,k)
                END IF
            END DO
        END DO
        IF(IDEN > 1) CALL CALC_SIGMAT_3D ()

    !     Compute Depth averaged SigmaT

        DO NH=1,NP
            DASigT(NH)=0.d0
            DO k=1,NFEN-1
                DASigT(NH)=DASigT(NH)+(Sigma(k+1)-Sigma(k)) &
                *(SigT(NH,k+1)+SigT(NH,k))/2.d0
            ENDDO
            DASigT(NH)=DASigT(NH)/AMB
        ENDDO
    ENDIF

!     Need to update and move results from the future time levels to the
!     present time levels so we can advance in time.
    DO NH = 1,NP
        DO k = 1,NFEN
            q(NH,k) = qkp1(NH,k)
            WZ(NH,k) = WZkp1(NH,k)
        END DO
    END DO

!     WRITE 3D MODEL OUTPUT
    CALL writeOutput3D(TimeLoc,IT) !jgf48_11 moved to write_output.F

!     RETURN TO THE 2-D MODEL TO COMPUTE H AT NEXT TIME STEP

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!***************************************************************************
    END SUBROUTINE VSSOL
!***************************************************************************




!******************************************************************************
!     Subroutine to compute the eddy viscosity profile.                       *
!                                                                             *
!  ievc, evmin, evcon - E.V. code, E.V. minimum value and E.V. constant       *
!                                                                             *
!        NOTE: evcon only used for some of the E.V. formulations as           *
!                  discussed below.                                           *
!        NOTE: In cases where EV is specified to vary linearly over the       *
!              lower 20% of the water column, it actually varies linearly     *
!              with a constant slope up to the vertical FE grid node that is  *
!              less than or equal to the 20% location.  The value is constant *
!              as specified at all FE grid nodes above the 20% location.      *
!              The E.V. above and below the 20% level is joined by one        *
!              additional linearly varying segment.                           *
!        NOTE: The E.V. is constrained to always be greater than or equal to  *
!              EVMIN as specified in the UNIT 15 file.                        *
!                                                                             *
!        ievc=0-9, EV constant in time & horizontal space                     *
!             0 - EV read in from UNIT 15 (may vary vertically) - EVCON is    *
!                    not used                                                 *
!             1 - EV = EVCON                                                  *
!                                                                             *
!        ievc=10-19 EV proportional to omega*h*h  (Lynch and Officer (1986)   *
!                                              Lynch and Werner (1987, 1991)) *
!             10 - EV = omega*h*h/10 over the entire water column             *
!             11 - EV = omega*h*h/1000 at bottom                              *
!                       varies linear over lower 20% of wc                    *
!                     = omega*h*h/10 in upper 80% of w.c.                     *
!            NOTE:For this EV formulation, evcon is not used and omega is     *
!                  hardwired for a 12.42 hour tide.                           *
!                                                                             *
!        ievc=20-29 EV proportional to kappa U* z                             *
!             20 - EV = 0.41U*Zo at bottom                                    *
!                     = 0.41U*Z over entire water column                      *
!             21 - EV = 0.41U*Zo at bottom                                    *
!                     = 0.41U*Z in lower 20% of water col                     *
!                     = 0.082U*h in upper 80% of water col                    *
!            WHERE: U* is the friction velocity                               *
!            NOTE: For this EV formulation, evcon is not used.                *
!                                                                             *
!        ievc=30-39, EV proportional to Uh (Davies 1990)                      *
!             30 - EV = 0.025|U|h/9.001 over entire water column              *
!             31 - EV = evcon|U|h over entire water column                    *
!             32 - EV = 0.025|U|h/9.001 in upper 80% of wc                    *
!                     = 0.000025h|U|/9.001 at bottom                          *
!                       varies linear over lower 20% of wc                    *
!             33 - EV = evcon|U|h in upper 80% of wc                          *
!                     = evcon|U|h/1000. at bottom                             *
!                       varies linear over lower 20% of wc                    *
!            WHERE: U is depth averaged velocity                              *
!            NOTE: For this EV formluation, evcon is used only for ievc=31,33 *
!                                                                             *
!        ievc=40-49, EV proportional to U*U (Davies 1990)                     *
!             40 - EV = 2|UU|/9.001 over entire water column                  *
!             41 - EV = evcon|UU| over entire water column                    *
!             42 - EV = 2|UU|/9.001 in upper 80% of wc                        *
!                     = 0.002|UU|/9.001 at bottom                             *
!                       varies linear over lower 20% of wc                    *
!             43 - EV = evcon|UU| in upper 80% of wc                          *
!                     = evcon|UU|/1000. at bottom                             *
!                       varies linear over lower 20% of wc                    *
!            WHERE: U is depth averaged velocity                              *
!            NOTE: For this EV formluation, evcon is used only for ievc=41,43 *
!                                                                             *
!        ievc=50, EV computed from Mellor-Yamada L2.5 closure                 *
!            NOTE: For this EV formulation, evcon is not used.                *
!            =51, EV computed from Mellor-Yamada L2.5 closure                 *
!                 with enhanced mixing in the surface layer                   *
!            NOTE: For this EV formulation, evcon is the value                *
!                  alpha coefficient to the surface roughness equation        *
!                  The TKE production parameter (a) is hardwired to           *
!                  a value of 150                                             *
!                                                                             *
!                               04/06/11                                      *
!******************************************************************************

    SUBROUTINE EDDYVIS(H,UU,VV,WSX,WSY,BSX,BSY,NODE,IT,Z0B1)
! jgfdebug USE GLOBAL, ONLY : ScreenUnit
! jgfdebug  USE GLOBAL_3DVS, ONLY : SZ,LOCALDIR,IEVC,EVMIN,EVCON,EVTOT,SIGMA,
! jgfdebug     &     NFEN,Z0B,Z0S,A,B,AMB,NSCREEN,NWS,DELT
    USE GLOBAL_3DVS, ONLY : SZ,IEVC,EVMIN,EVCON,EVTOT,SIGMA, &
    NFEN,Z0B,Z0S,A,B,AMB,setMessageSource, unsetMessageSource, &
    logMessage, allMessage, DEBUG, ECHO, INFO, WARNING, ERROR
! jwdebug  USE NodalAttributes, ONLY : Z0B_var
    IMPLICIT NONE
    REAL(8) :: H
    REAL(SZ) :: UU,VV,WSX,WSY,BSX,BSY
    REAL(SZ) :: Z0B1

!     jgf46.00 explicitly declared the following variables
    INTEGER :: NODE
    INTEGER, SAVE :: istart = 0
    INTEGER :: IT
    INTEGER :: J
    REAL(SZ) RKAPPA
    REAL(SZ) OMEGA
    REAL(SZ) EVBASE
    REAL(SZ) SLOPE
    REAL(SZ) BREAK
    REAL(SZ) EVBEGIN
    REAL(SZ) USTARB
    REAL(SZ) USTARS
    REAL(SZ) UAVMAG
    REAL(SZ) UAVMAGS

    call setMessageSource("eddyvis")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

! jgfdebug 350  FORMAT(//,2X,'***** INVALID INPUT IN THE PRIMARY VERTICAL INPUT',
! jgfdebug     &     ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

    RKAPPA = 0.41

!...
!...  IEVC=0 READ IN FROM UNIT 15 in subroutine READ_INPUT_3DVS
!...

!...
!...  IEVC=1 SET EV = EVCON ON THE FIRST TIME STEP
!     .

    IF((IEVC == 1) .AND. (istart == 0)) THEN
        istart = 1
        DO J=1,NFEN
            EVTOT(J)=EVCON
        ENDDO
    ENDIF

    IF(IEVC == 1) GOTO 100

!...
!...  OMEGA*H*H formulation FOLLOWING
!...  LYNCH AND OFFICER (1986), LYNCH AND WERNER (1987, 1991)
!...
    IF(IEVC == 10) THEN
        OMEGA=1.40525D-4
        EVBASE=OMEGA*H*H/10.D0
        IF(EVBASE < EVMIN) EVBASE=EVMIN
        DO J=1,NFEN
            EVTOT(J) = EVBASE
        END DO
        GOTO 100
    ENDIF
    IF(IEVC == 11) THEN
        OMEGA=1.40525D-4
        EVBASE=OMEGA*H*H/10.D0
        EVTOT(1)=EVBASE/100.D0
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-b*SLOPE
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF

!...
!...  KAPPA USTAR Z FORMULATION
!...
    IF(IEVC == 20) THEN
        USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))
        EVTOT(1)=RKAPPA*USTARB*Z0B1
        EVBASE=RKAPPA*USTARB*H
        IF(EVBASE < EVTOT(1)) EVBASE=EVTOT(1)
        SLOPE=(EVBASE-EVTOT(1))/amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J) !=EVTOT(1)+SLOPE*(SIGMA(J)-b)
        ENDDO
        GOTO 100
    ENDIF
    IF(IEVC == 21) THEN
        USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))
        EVTOT(1)=RKAPPA*USTARB*Z0B1
        EVBASE=RKAPPA*USTARB*H*0.2D0
        IF(EVBASE < EVTOT(1)) EVBASE=EVTOT(1)
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF
    IF(IEVC == 22) THEN
        USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))

        USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
        EVTOT(1)=RKAPPA*0.5*(USTARB+USTARS)*Z0B1
        EVBASE=RKAPPA*0.5*(USTARB+USTARS)*H*0.2D0
        IF(EVBASE < EVTOT(1)) EVBASE=EVTOT(1)
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF
    IF(IEVC == 23) THEN
        USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))

        USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
        EVTOT(1)=0.1*0.5*(USTARB+USTARS)*H
        DO J=2,NFEN
            EVTOT(J)=0.1*0.5*(USTARB+USTARS)*H
        ENDDO
        GOTO 100
    ENDIF

!...
!...  H UAVG FORMULATION FOLLOWING DAVIES (1990) + [EQ. (33)]
!...
    IF(IEVC == 30) THEN
        UAVMAG=SQRT(UU*UU+VV*VV)
        EVBASE=0.025*H*UAVMAG/9.001D0
        DO J=1,NFEN
            EVTOT(J) = EVBASE
        END DO
        GOTO 100
    ENDIF
    IF(IEVC == 31) THEN
        UAVMAG=SQRT(UU*UU+VV*VV)
        EVBASE = EVCON*H*UAVMAG
        DO J=1,NFEN
            EVTOT(J) = EVBASE
        END DO
        GOTO 100
    ENDIF
    IF(IEVC == 32) THEN
        UAVMAG=SQRT(UU*UU+VV*VV)
        EVBASE = 0.025*H*UAVMAG/9.001D0
        EVTOT(1) = EVBASE/1000.D0
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF
    IF(IEVC == 33) THEN
        UAVMAG=SQRT(UU*UU+VV*VV)
        EVBASE=EVCON*H*UAVMAG
        EVTOT(1)=EVBASE/1000.D0
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF

!...
!...  UAVG SQUARED FORMULATION FOLLOWING DAVIES (1990) + [EQ. (34)]
!...
    IF(IEVC == 40) THEN
        UAVMAGS=UU*UU+VV*VV
        EVBASE=2.D0*UAVMAGS/9.001D0
        DO J=1,NFEN
            EVTOT(J) = EVBASE
        END DO
        GOTO 100
    ENDIF
    IF(IEVC == 41) THEN
        UAVMAGS=UU*UU+VV*VV
        EVBASE=EVCON*UAVMAGS
        DO J=1,NFEN
            EVTOT(J) = EVBASE
        END DO
        GOTO 100
    ENDIF
    IF(IEVC == 42) THEN
        UAVMAGS=UU*UU+VV*VV
        EVBASE=2.D0*H*UAVMAGS/9.001D0
        EVTOT(1)=EVBASE/1000.D0
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF
    IF(IEVC == 43) THEN
        UAVMAGS=UU*UU+VV*VV
        EVBASE=EVCON*UAVMAGS
        EVTOT(1)=EVBASE/1000.D0
        SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
        BREAK=b+0.2D0*amb
        EVBEGIN=EVTOT(1)-SLOPE*b
        DO J=2,NFEN
            IF(SIGMA(J) <= BREAK) THEN
                EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
                EVTOT(J)=EVBASE
            ENDIF
        ENDDO
        GOTO 100
    ENDIF

!...
!...  MELLOR-YAMADA LEVEL 2.5 CLOSURE
!...
    IF((IEVC == 50) .OR. (IEVC == 51)) THEN
        CALL TURB(NODE,H,BSX,BSY,WSX,WSY,IStart,Z0B1)

        GOTO 100
    ENDIF

!...
!...  ONCE EDDY VISCOSITY IS COMPUTED
!...


    100 CONTINUE

!...
!...  CHECK SO THAT EDDY VISCOSITY NEVER GETS BELOW MINIMUM VALUE
!...
    DO J=1,NFEN
        IF(EVTOT(J) < EVMIN) THEN
            EVTOT(J)=EVMIN
        ENDIF
    ENDDO
!...
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!***********************************************************************
    END SUBROUTINE EDDYVIS
!***********************************************************************



!***********************************************************************
!                                                                      *
!   Solver for a vector U of length nfen for a tridiagonal system of   *
!          equations with the form                                     *
!                                                                      *
!    **                               **   **    **   **    **         *
!    * Bn   An                         *   *  Un  *   *  Rn  *         *
!    *                                 *   *      *   *      *         *
!    * Cn-1 Bn-1 An-1                  *   * Un-1 *   * Rn-1 *         *
!    *                                 *   *      *   *      *         *
!    *      Cn-2 Bn-2 An-2             *   * Un-2 *   * Rn-2 *         *
!    *                                 *   *      *   *      *         *
!    *             .......             *   * .... *   * .... *         *
!    *                                 *   *      *   *      *         *
!    *                  C3 B3 A3       *   *  U3  * = *  R3  *         *
!    *                                 *   *      *   *      *         *
!    *                     C2 B2 A2    *   *  U2  *   *  R2  *         *
!    *                                 *   *      *   *      *         *
!    *                        C1 B1    *   *  U1  *   *  R1  *         *
!    **                               **   **    **   **    **         *
!                                                                      *
!      A, B, C, U, R are adjustable size arrays                        *
!                                                                      *
!      U(1 - nfen) are the complex velocities from bottom to top       *
!                                                                      *
!                                                                      *
!                         R.L.  11/18/94                               *
!***********************************************************************

!  kmd48.33bc changed to use Implicit None
    SUBROUTINE TRIDIAG(A,B,C,R,U,nfen)
    USE GLOBAL, ONLY : SZ, ScreenUnit, setMessageSource, &
    unsetMessageSource, logMessage, allMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE
    INTEGER :: J,nfen

    COMPLEX(SZ) :: A(nfen), B(nfen), C(nfen), R(nfen), U(nfen)
    COMPLEX(SZ) :: P1

    call setMessageSource("tridiag")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    DO 11 J=nfen-1,1,-1
        IF(ABS(B(J+1)) == 0.) then
            write(screenunit,*) 'Diagonal term in the VS matrix is zero'
            write(screenunit,*) '*********** Fatal error *************'
            write(1,*) 'Diagonal term in the VS matrix is zero'
            write(1,*) '*********** Fatal error *************'
            stop
        endif
        P1=C(J)/B(J+1)
        B(J)=B(J)-A(J+1)*P1
        R(J)=R(J)-R(J+1)*P1
    11 END DO

    IF(ABS(B(1)) == 0.) then
        write(screenunit,*) 'B1 term in the VS matrix is zero'
        write(screenunit,*) '*********** Fatal error *************'
        write(1,*) 'B1 term in the VS matrix is zero'
        write(1,*) '*********** Fatal error *************'
        stop
    endif

    U(1)=R(1)/B(1)

    DO 12 J=2,nfen
        U(J)=(R(J)-A(J)*U(J-1))/B(J)
    12 END DO

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
!***********************************************************************
    END SUBROUTINE TRIDIAG
!***********************************************************************


!***********************************************************************
!                                                                      *
!                    VSDISP.FOR - VERSION                              *
!                                                                      *
!  This subroutine computes the dispersion terms: Duu, Duv, Dvv        *
!  for the FE VS method as derived in "Dispersion Terms 1/27/92"       *
!                                                                      *
!                                                                      *
!                          R.L. 05/25/00                               *
!                          R.L. 06/22/05                               *
!***********************************************************************

    Subroutine VSDISP(it,node,H,UU,VV,Duu,Duv,Dvv)
    USE GLOBAL, ONLY : ScreenUnit, allMessage, DEBUG, scratchMessage
! jgfdebug      USE GLOBAL_3DVS, ONLY : Qkp1,SZ,NFEN,SIGMA,GAMMA,AMB,I,NSCREEN
    USE GLOBAL_3DVS, ONLY : Qkp1,SZ,NFEN,SIGMA,GAMMA,AMB,IY, &
    setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE

    REAL(SZ) :: UU,VV,Duu,Duv,Dvv
    REAL(8) :: delst,D1P1,Txx,Txy,Tyy,Ump,Vmp,Um,Vm,H
!     jgf46.00 explicitly declared the following variables
    INTEGER :: it
    INTEGER :: node
! jgfdebug      INTEGER ierr
    INTEGER :: j
! jgfdebug      INTEGER k
    INTEGER :: m
! jgfdebug      INTEGER n
    INTEGER :: mp
! jgfdebug      INTEGER istart
! jgfdebug      INTEGER il
! jgfdebug      INTEGER ibc
! jgfdebug      REAL(SZ) rkap
! jgfdebug      REAL(SZ) B1
! jgfdebug      REAL(SZ) B123
! jgfdebug      REAL(SZ) g2
! jgfdebug      REAL(SZ) g3
! jgfdebug      REAL(SZ) g4
! jgfdebug      REAL(SZ) g5
! jgfdebug      REAL(SZ) g6
! jgfdebug      REAL(SZ) E1

!     parameters

! jgfdebug      ierr = 0

!     don't waste time if entire profile is zero

    call setMessageSource("vsdisp")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    do j=1,nfen
        gamma(j)=Qkp1(node,j)
    end do

    do j=1,nfen
        if(gamma(j) /= (0.d0,0.d0)) goto 1
    end do
    Duu = 0.d0
    Duv = 0.d0
    Dvv = 0.d0
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    return


!     Domain is split up into intervals according to sigma grid and
!     contribution to total integral is computed on each interval.


    1 Txx = 0.d0
    Txy = 0.d0
    Tyy = 0.d0
    Ump = Real(gamma(1))
    Vmp = -Real(iy*gamma(1))
    Do 100 m=1,nfen-1
        mp = m+1
        Um = Ump
        Vm = Vmp
        Ump = Real(gamma(mp))
        Vmp = -Real(iy*gamma(mp))
        delst = sigma(mp) - sigma(m)
        D1P1 = delst/3.d0
        Txx = Txx + (Um*Um + Ump*Ump + Um*Ump)*D1P1
        Txy = Txy + ((Um*Vm + Ump*Vmp) + (Um*Vmp + Vm*Ump)/2.d0)*D1P1
        Tyy = Tyy + (Vm*Vm + Vmp*Vmp + Vm*Vmp)*D1P1
    100 END DO

    if(Txx < 0.) then
        write(screenunit,1001)
        1001 format(/'**** Serious Error Detected in SUBROUTINE VSDISP ***'/ &
        '          the partial uu dispersion term < 0        '/ &
        '        Diagnostic information written to unit 1    ')
        write(screenunit,1101) node,it
        write(1,1002)
        1002 format(/'**** Serious Error Detected in SUBROUTINE VSDISP ***'/ &
        '          the partial uu dispersion term < 0        ')
        write(1,1101) node,it
        write(1,1021)
        do m=1,NFEN            !ntotn -> NFEN
            write(1,1022) sigMA(m), gamma(m) !sigtot -> sigMA
        enddo
        write(1,2101) Txx
        2101 format('  Par uu = ',e14.6,' It will be set to 0.')
        Txx = 0.d0
    endif

    if(Tyy < 0.) then
        write(screenunit,1003)
        1003 format(/'**** Serious Error Detected in SUBROUTINE VSSOL ****'/ &
        '          the partial vv dispersion term < 0        '/ &
        '       Diagnostic information written to unit 1    ')
        write(screenunit,1101) node,it
        write(1,1004)
        1004 format(/'**** Serious Error Detected in SUBROUTINE VSSOL ****'/ &
        '          the partial vv dispersion term < 0        ')
        write(1,1101) node,it
        write(1,1021)
        do m=1,NFEN            !ntotn -> NFEN
            write(1,1022) sigMA(m),gamma(m) !sigtot -> sigMA
        enddo
        write(1,2102) Tyy
        2102 format('  Par vv = ',e14.6,' It will be set to 0.')
        Tyy = 0.d0
    endif

    1101 format(/' Error occurred at node ',I6,' time step ',I8/)
    1021 format(8x,'sigMA',11x,'  u ',12x,'  v ')
    1022 format(1x,3e16.8)

    Duu = Txx*H/amb - H*UU*UU
    Duv = Txy*H/amb - H*UU*VV
    Dvv = Tyy*H/amb - H*VV*VV

    if(Duu < 0.) Duu = 0.d0
    if(Dvv < 0.) Dvv = 0.d0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    return
!***********************************************************************
    end subroutine vsdisp
!***********************************************************************



!***********************************************************************
!     Subroutine to compute the Inm integral                           *
!                                                                      *
!     Note, Inm is based only on the f.e. grid and therefore is the    *
!     same for all horizontal nodes.                                   *
!                                                                      *
!     11/24/01                                                         *
!                                                                      *
!***********************************************************************

    subroutine InmINT()
! jgfdebug      USE GLOBAL, ONLY : ScreenUnit
    USE GLOBAL_3DVS, ONLY : Sigma, NFEN, INM, setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE
    INTEGER :: k

    call setMessageSource("InmINT")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!     integral over lower element of psi(k-1)*psi(k)
    Inm(1,1) = 0.d0
!     integral over upper element of psi(k+1)*psi(k)
    Inm(1,3) = (Sigma(2)-Sigma(1))/6.d0
!     integral over both elements of psi(k)*psi(k)
    Inm(1,2) = 2.d0*Inm(1,3)

    do k=2,NFEN-1
    !     integral over lower element of psi(k-1)*psi(k)
        Inm(k,1) = Inm(k-1,3)
    !     integral over upper element of psi(k+1)*psi(k)
        Inm(k,3) = (Sigma(k+1)-Sigma(k))/6.d0
    !     integral over both elements of psi(k)*psi(k)
        Inm(k,2) = 2.d0*(Inm(k,1)+Inm(k,3))
    enddo

!     integral over lower element of psi(k-1)*psi(k)
    Inm(NFEN,1) = Inm(NFEN-1,3)
!     integral over both elements of psi(k)*psi(k)
    Inm(NFEN,2) = 2.d0*Inm(NFEN,1)
!     integral over upper element of psi(k+1)*psi(k)
    Inm(NFEN,3) = 0.d0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    return
!***********************************************************************
    end subroutine InmINT
!***********************************************************************


!***********************************************************************
!      Subroutine to compute the LVn integral (used only in VS)        *
!                                                                      *
!                            11/14/94                                  *
!***********************************************************************

    subroutine LVnInt()
! jgfdebug      USE GLOBAL, ONLY : ScreenUnit
    USE GLOBAL_3DVS, ONLY : Sigma, NFEN, LVN, setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE
    INTEGER :: n

    call setMessageSource("LVnInt")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif


    LVn(1) = (sigma(2) - sigma(1))/2.d0
    DO n=2,nfen-1
        LVn(n) = (sigma(n+1) - sigma(n-1))/2.d0
    enddo
    LVn(nfen) = (sigma(nfen) - sigma(nfen-1))/2.d0


#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!***********************************************************************
    END subroutine LVnInt
!***********************************************************************

!***********************************************************************
!      Subroutine to compute the KQnm integral                         *
!                                                                      *
!                            01/26/00                                  *
!***********************************************************************

    subroutine KQnmInt(KQnm,Kq)
! jgfdebug      USE GLOBAL, ONLY : ScreenUnit
    USE GLOBAL_3DVS, ONLY : SZ,SIGMA,NFEN,MNFEN, setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE
    REAL(SZ) :: KQnm(MNFEN,3),Kq(MNFEN)
    REAL(8) :: EM,EP
    INTEGER :: n

    call setMessageSource("KQnmInt")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif



!     do 1st element by hand

    EP=(Kq(2)+Kq(1))/(sigma(2)-sigma(1))/2.d0
    KQnm(1,1)=0.D0
    KQnm(1,2)=EP
    KQnm(1,3)=-EP

!     loop through interior elements
!     NOTE: the integrals from sigma(n-1) to sigma(n) are simply integrals
!     the integrals from sigma(n) to sigma(n+1) from previous element

    DO n=2,nfen-1
        EM=EP
        EP=(Kq(n+1)+Kq(n))/(sigma(n+1)-sigma(n))/2.d0
        KQnm(n,1)=-EM
        KQnm(n,2)=EM+EP
        KQnm(n,3)=-EP
    ENDDO

!     do last element by hand

    EM=EP
    KQnm(nfen,1)=-EM
    KQnm(nfen,2)=EM
    KQnm(nfen,3)=0.D0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    return
!************************************************************************
    end subroutine KQnmInt
!************************************************************************

!************************************************************************
! MY2.5 TURBULENCE MODEL PROGRAM
! VELOCITY SOLUTION VERSION

! written by R. Luettich based on earlier subroutines by R. Grenier

! OVERVIEW

!     This code uses the quasi-equilibrium version of the Mellor-Yamada
!     turbulence scheme (Mellor and Yamada, 1982, Blumberg and Mellor,
!     1987 and Galperin et al., 1988) to solve transport equations for
!     q**2 and q**2l.  The parameters q and l are used to compute the
!     eddy viscosity according to the relation:

!             Km = Smql

!     where Km is the eddy viscosity and Sm is a stability parameter.

!     This routine is called from the EDDYVIS subroutine during the
!     internal mode solution for ievc = 50.  There are separate DSS and
!     VS versions of the code.

!    *** SEE SECTION TITLED "USERS GUIDE" BELOW ***

!     This routine should be linked with the 3D Code subroutines (ADCIRC
!     + VS) after running the setup program

! PARAMETER DEFINITIONS:

! argument list:

! nh          : horizontal node counter
! H           : total depth at current time step (includes finite amplitude)
! it          : current time step value
! delt        : simulation incremental time step (uses same step as
!               internal and external modes)
! nws         : switch for wind stress (nws=1 for wind on,
!               nws=0 for wind off). This affects the surface
!               boundary condition. Set in the external mode and passed in.
! bsx,bsy     : x,y bottom stresses
! wsx,wsy     : x,y surface stresses
! Z0B         : bottom roughness/mixing length (see USER's notes
!               below)
! Z0S         : surface roughness/mixing length (see USER's notes
!               below)
! Z0Sw        : surface roughness/mixing length (see USER's notes
!               below) if Z0S is greater than or equal to water depth

! Coeff block:

! Sq          : stability function used in definition of the eddy
!               diffusivity of the turbulence parameters (Sq = 0.2)
! Sm          : stability function used in definition of the eddy
!               viscosity (see Galperin, et al. for exact form)
! Sh          : stability function used in definition of the eddy
!               diffusivity of the density (see Galperin, et al.
!               for exact form)
! B1,E1,E2,E3 : empirical constants (values given below; see Mellor and
!               Yamada, 1982 & Blumberg et al, 1992 for discussion)

! Turbmod block:

! Q (real, imag): x,y direction dependent variables from the vertical
!               solution

! Vgrid block, Setup Block

! parameters as defined in the internal mode subroutines

! Other calculation parameters:

! q2l         : 2xTKE times the master length scale
! q20         : previous time step value of q2
! l           : master length scale
! Km          : momentum eddy viscosity
! Kq          : turbulence eddy viscosity
! Kh          :
! Mqa,Mqb,Mqc : LHS matrix diagonals for the q2 and q2l solutions
! LVq         : RHS load vector for the q2 and q2l solutions
! q2          : solution vector returned from the tridiagonal solver
! w           : wall function (see Blumberg at al., 1992)
! BVflux2     : Brunt Vaislai frequency squared
! Gh          : dimensionless density function
! SIGT        : profile of density (sigma T)

! Flags, run control and miscellaneous parameters

! iden        : density flag (iden <> 0 density included,iden = 0 no density)
! im          : run type flag - vs (im = 1) or dss (im = 2)
! il          : length scale flag - il = 2 or 3 for algebraic length scale (see below),
!               il = 1 for length scale computed from the q2l equation (see user notes below)
! ibc         : boundary condition flag for q2 (0=constant/zero,1=no-flux)
!               User should set ibc = 0 for nws = 1


! SUBROUTINES

! turb        : main module - handles input, run control and output.
! TRIDAG      : tridiagonal matrix solver


! USER'S GUIDE

!     The user must set a number of parameter statements and flags prior
!     to operation:

! parameter statements **Check all subroutines**
!  MNFEN = maximum number of vertical nodes
!  mnp = maximum number of horizontal nodes

! flags (SEE DISCUSSION ABOVE UNDER flags, run control and miscellaneous parameters)

!     The value of Z0B should be set consistent with the bottom boundary
!     condition used in the internal mode solution.  For a no-slip BBC,
!     choose Z0B as a physical roughness height, e.g. 0.005m.  For a
!     slip case, set Z0B to a value consistent with the slip coefficient
!     used via the log profile, (e.g., 1 m).  Tests suggest that when a
!     no-slip condition is used with the VS model (for which the bottom
!     nodes are very tightly spaced) it is best to set the time
!     weighting parameter for the momentum diffusion term (alpha3 in
!     input unit 15) to 1.0 to avoid instability problems. When a slip
!     condition is used (or whenever the bottom grid spacing is about
!     1m) a Crank-Nicholson approach (alpha3 = 0.5) is acceptable.

!     The number of vertical nodes used in the solution of the
!     turbulence equations is the same as the number of nodes used in
!     the solution of the dependent variable (velocity or stress) in the
!     internal mode.  This is unlike the other forms of eddy viscosity,
!     for which the two grids are different and the number of nodes used
!     to define the eddy viscosity is generally less.

!     This code assumes that any density field is passed into the
!     routine via a common block called "DENSITY3D", which includes a
!     density profile on the internal mode solution grid at each point
!     in the horizontal.  Any updating of this profile must be done
!     externally and passed into this routine.  The density be passed as
!     sigma t units, and the background density is RHOWAT0.

!     Model output is limited to printing results for a single
!     horizontal node.  Additional coding would be required to create
!     full output files.


! REFERENCES

!     Blumberg, A.F. and G.L. Mellor, A Description of a
!     Three-Dimensional Coastal Ocean Circulation Model, In:
!     Three-Dimensional Coastal Ocean Models, edited by N.S. Heaps,
!     pp. 1-16, American Geophysical Union, Washington, D.C., 1987.

!     Blumberg, A.F., B. Galperin and D.J. O'Connor, Modeling vertical
!     structure of open channel flows, Journal of Hydraulic Engineering,
!     118, 1119-1134., 1992.

!     Galperin, B., L.H. Kantha, S. Hassid and A. Rosati, A
!     quasi-equilibrium turbulent energy model for geophysical flows,
!     Journal of the Atmospheric Sciences, 45, 55-62, 1988.

!     Mellor, G.L. and T. Yamada, Development of a turbulence closure
!     model for geophysical fluid problems, Reviews of Geophysics and
!     Space Physics, 20, 851-875, 1982.

!************************************************************************
    subroutine turb(nh,H,bsx_loc,bsy_loc,wsx,wsy,istart,Z0B1)
    USE SIZES
! jgfdebug      USE GLOBAL, ONLY : ScreenUnit
! jgfdebug      USE GLOBAL_3DVS, EXCEPT_BSX => BSX ,EXCEPT_BSY => BSY
    USE GLOBAL_3DVS
! jwdebug     USE NodalAttributes, ONLY : Z0B_var
!>>>>>>> IOOS_w_3D/src/vsmy.F
    IMPLICIT NONE

    COMPLEX(SZ) :: dQdz,dQdz1,dQdz2
    REAL(8) :: H,H_crit

    REAL(SZ) :: bsx_loc,bsy_loc,WSX,WSY
    REAL(SZ),SAVE :: H2
!  kmd48.33bc need to change Sh to an array for transport
    REAL(SZ),SAVE,ALLOCATABLE :: Sh(:)
    REAL(SZ) :: Z0Sw
    REAL(SZ) :: Z0B1
    REAL(SZ),SAVE,ALLOCATABLE :: KQnm(:,:)
    REAL(SZ),SAVE,ALLOCATABLE :: Mqa(:),Mqb(:),Mqc(:)
    REAL(SZ),SAVE,ALLOCATABLE :: LVq(:),Sm(:)
    REAL(SZ),SAVE,ALLOCATABLE :: q2(:),q2prev(:)
    REAL(SZ),SAVE,ALLOCATABLE :: q2l(:),q2lprev(:)
    REAL(SZ),SAVE,ALLOCATABLE :: wall(:),rmlen(:),rmlen2(:)
    REAL(SZ),SAVE,ALLOCATABLE :: BVfreq2(:),spgrad2(:)
    REAL(SZ),SAVE,ALLOCATABLE :: Kq(:),Km(:),Kh(:)
    REAL(SZ),SAVE,ALLOCATABLE :: prod(:),diss(:)
!     jgf46.00 explicitly declared the following variables
!     jgf46.12 gave them SAVE and PARAMETER properties (from jgf45.18)
    INTEGER :: nh, istart, n
    INTEGER,SAVE :: il, ibc
    REAL(SZ), PARAMETER :: rkap = 0.41d0
    REAL(SZ), PARAMETER :: B1 = 16.6d0
    REAL(SZ), PARAMETER :: g2 = 0.39327d0
    REAL(SZ), PARAMETER :: g3 = 3.0858d0
    REAL(SZ), PARAMETER :: g4 = 34.676d0
    REAL(SZ), PARAMETER :: g5 = 6.1272d0
    REAL(SZ), PARAMETER :: g6 = 0.49393d0
    REAL(SZ), PARAMETER :: E1 = 1.8d0
    REAL(SZ), PARAMETER :: E2 = 1.33d0
    REAL(SZ), PARAMETER :: E3 = 0.25d0
    REAL(SZ), PARAMETER :: q2min=1.d-8
    REAL(SZ), SAVE :: B123, Sq, sig, rl1, rl2, HOamb, zval, HOamb2
    REAL(SZ), SAVE :: dsig, drhodz, dudz, dvdz, dsig1, dsig2
    REAL(SZ), SAVE :: drhodz1, drhodz2, BSlay, SSlay, db, db2, ds, ds2
    REAL(SZ), SAVE :: qprev, qlprev, Gh, tdiss, elmax
    REAL(SZ), SAVE :: coef1, coef2, coef3, coef4, coef5
    REAL(SZ) :: USTARS

    call setMessageSource("turb")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))

!     allocate local arrays
    if( .NOT. turb_allocated) then
        allocate(KQnm(MNFEN,3))
        allocate(Mqa(MNFEN),Mqb(MNFEN),Mqc(MNFEN))
        allocate(LVq(MNFEN),Sm(MNFEN))
        allocate(q2(MNFEN),q2prev(MNFEN))
        allocate(q2l(MNFEN),q2lprev(MNFEN))
        allocate(wall(MNFEN),rmlen(MNFEN),rmlen2(MNFEN))
        allocate(BVfreq2(MNFEN),spgrad2(MNFEN))
        allocate(Kq(MNFEN),Km(MNFEN),Kh(MNFEN))
        allocate(prod(MNFEN),diss(MNFEN))
        allocate(Sh(MNFEN))  ! kmd48.33bc added Sh as an array
        turb_allocated = .TRUE. 
    endif


    Z0Sw=Z0S
    IF (IEVC == 51) THEN
        Z0Sw=EVCON*USTARS*USTARS/g
    ELSE
        Z0Sw=Z0S
    ENDIF
!  RJW taken from Jones monismith 2008
!  modify surface roughness length to account for wave breaking
!         Z0Sw = 1400.0d0*USTARS*USTARS/g   ! ( Bye (1988)  , Craig and Banner (1994))
!         Z0Sw = 14000.0d0*USTARS*USTARS/g  ! (Stipps et al. 2005 extended to wave breaking)
!         Z0Sw = 200000.0d0*USTARS*USTARS/g  ! ( O(10^5)  Stacy et al. (1999) for wave ages > 25)
!         Z0Sw = 1.3*Wave Height            ! (sig wave height) (Jones and Monosmith (2007) k-omega)
!         Z0Sw = (0.25-1.0)*Wave Height     ! (sig wave height) (Burchard (2001) k-epsilon )

    if (Z0Sw > H-Z0B1) then
        Z0Sw = H-Z0B1 ! bound Z0S to water depth
    else
        Z0Sw=Z0Sw
    endif
    if (Z0Sw < Z0S ) Z0Sw = Z0S



!     At time step 1, set flags and initialize variables

!.  RJW moved up and outside of if statement to allow variable computation in case of hotstart
    il = 1  ! length scale flag, = 1 for length scale from q2l eqn.
!         il = 2  ! length scale flag, = 2 for length scale from algebraic eq. (see below)
!         il = 3  ! length scale flag, = 3 for length scale from algebraic eq. (see below)
    ibc = 1
! bc = 0 ! surface b.c. flag, = 0 for specified stress, = 1 for no flux
! bc = 1 ! surface b.c. flag, = 0 for specified stress, = 1 for no flux
! f((nws.eq.0).or.(nws.eq.100)) ibc =1
!     set constants and other parameters
    B123 = B1**(2.d0/3.d0)
    Sq = 0.2d0 ! initialize the stability constant stability function
!     initialize the Brunt-Vaisala freq squared = 0 if density not considered
    if(iden == 0) then
        do n=1,nfen
            BVfreq2(n)=0.d0
        end do
    endif

    if(istart == 1) then
    !     initialization for a cold start only
        if(ihot == 0) then
            do n = 1,nfen
                q20(nh,n) = q2min ! initialize q2 to a minimal value
            !     initilize l to a minimal value if computed from q2l equation
                if(il == 1) then
                !           l(nh,n)=rkap*Z0B1                        !rog way
                    l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
                endif
            !     set l to an exponential type length scale (davies and xing)
                if(il == 2)then
                    sig = (sigma(n)+1.d0)/amb
                    rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
                    rl2 = 1.d0/(rkap*(H-sig*H+Z0Sw))
                    l(nh,n) = 1.d0/(rl1+rl2)
                endif
            !     set l to a linear variation with kz over lower 15% with constant above
                if(il == 3)then
                    HOamb=H/amb
                    zval = (sigma(n)+1.d0)*HOamb-H
                    if(sigma(n) <= -0.7d0) l(nh,n)=rkap*(H+zval+Z0B1)
                    if(sigma(n) > -0.7d0) l(nh,n)=rkap*(0.15d0*H+Z0B1)
                endif
            end do
        endif
    !     end cold start initialization
    endif                     !end of 1st time step section

!     Begin calculations for each time step

! RJW..
! identify if node was dry last timestep but is now wet
! the initialized value of L can be recomputed to a more meaningful
! value before used to compute q2

    if (L(nh,1) == -9999d0 ) Then
        do n = 1,nfen
            l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
        enddo
    endif
! now check q20 it may also be zero
    if (q20(nh,1) == -9999d0 ) Then
        do n = 1,nfen
            q20(nh,n) = q2min
        end do
    Endif

!         if(il.eq.2)then
!                  sig = (sigma(n)+1.d0)/amb
!                  rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
!                  rl2 = 1.d0/(rkap*(H-sig*H+Z0S))
!                  l(nh,n) = 1.d0/(rl1+rl2)
!               endif

! jgfdebug      H2 = H*H
    HOamb=H/amb
    HOamb2=(H/amb)*(H/amb)

!     Compute the speed gradient squared, density gradient, BV freq
!     and split out the mixing length

    rmlen(1)=l(nh,1)
    rmlen2(1)=rmlen(1)*rmlen(1)
    dsig=sigma(2)-sigma(1)
    if(iden /= 0) then
        drhodz=((SIGT(nh,2)-SIGT(nh,1))/dsig)/HOamb
        BVfreq2(1)=-GORHO*drhodz
    endif
    dQdz=((Q(nh,2)-Q(nh,1))/dsig)/HOamb
    dudz=real(dQdz)
    dvdz=aimag(dQdz)
    spgrad2(1)=dudz*dudz+dvdz*dvdz

    do n=2,nfen-1
        rmlen(n)=l(nh,n)
        rmlen2(n)=rmlen(n)*rmlen(n)
        dsig1=sigma(n+1)-sigma(n)
        dsig2=sigma(n)-sigma(n-1)
        if(iden /= 0) then
            drhodz1=((SIGT(nh,n+1)-SIGT(nh,n))/dsig1)/HOamb
            drhodz2=((SIGT(nh,n)-SIGT(nh,n-1))/dsig2)/HOamb
            BVfreq2(n)=-GORHO*(drhodz1+drhodz2)/2.d0
        endif
        dQdz1=((Q(nh,n+1)-Q(nh,n))/dsig1)/HOamb
        dQdz2=((Q(nh,n)-Q(nh,n-1))/dsig2)/HOamb
        dQdz=(dQdz1+dQdz2)/2.d0
        dudz=real(dQdz)
        dvdz=aimag(dQdz)
        spgrad2(n)=dudz*dudz+dvdz*dvdz
    enddo

    rmlen(nfen)=l(nh,nfen)
    rmlen2(nfen)=rmlen(nfen)*rmlen(nfen)
    dsig=sigma(nfen)-sigma(nfen-1)
    if(iden /= 0) then
        drhodz=((SIGT(nh,nfen)-SIGT(nh,nfen-1))/dsig)/HOamb
        BVfreq2(nfen)=-GORHO*drhodz
    endif
    dQdz=((Q(nh,nfen)-Q(nh,nfen-1))/dsig)/HOamb
    dudz=real(dQdz)
    dvdz=aimag(dQdz)
    spgrad2(nfen)=dudz*dudz+dvdz*dvdz

!     Compute the wall function if the mixing length is determined from
!     q2l eqn

    if(il == 1)then
        BSlay=Z0B1
        SSlay=Z0Sw
        do n = 1,nfen
            db=(HOamb*(sigma(n)-b)+BSlay)*rkap
            db2=db*db
            ds=(HOamb*(a-sigma(n))+SSlay)*rkap
            ds2=ds*ds
            wall(n) = 1.d0 + E2*rmlen2(n)/db2 + E3*rmlen2(n)/ds2
        enddo
    endif

!     Compute the stability functions, eddy viscosity and partial
!     turbulence production & dissipation terms using information from
!     the previous time step
! RJW rweaver add wave production of TKE
! alpha * USTARS*USTARS*USTARS, USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
! and alpha = 60 - 250

    do n = 1,nfen
        q2prev(n)=q20(nh,n)
        q2lprev(n)=q2prev(n)*rmlen(n)
        qprev=sqrt(q2prev(n))
        qlprev=qprev*rmlen(n)
        Gh=-BVfreq2(n)*rmlen2(n)/q2prev(n)
        if(Gh > 0.0233) Gh=0.0233
        Sm(n)=(g2-g3*Gh)/((1.d0-g4*Gh)*(1.d0-g5*Gh))
        Sh(n)=g6/(1.d0-g4*Gh) ! kmd48.33bc changed Sh to an array
        Km(n)=Sm(n)*qlprev
        if(Km(n) < EVMIN) Km(n)=EVMIN
    !     if(istart.eq.1) Km(n)=EVMIN                        !rog way
        Kq(n)=Sq*qlprev
        if(Kq(n) < EVMIN) Kq(n)=EVMIN
    !     Kq(n)=Km(n)*Sq/Sm(n)                           !rog way
        Kh(n)=Sh(n)*qlprev ! kmd48.33bc changed Sh to an array
        if(Kh(n) < EVMIN) Kh(n)=EVMIN
        prod(n)=Km(n)*spgrad2(n)-Kh(n)*BVfreq2(n)
        diss(n)=qprev/(B1*rmlen(n))
    !     diss(n)=q2prev(n)*Sm(n)/(Km(n)*B1)             !rog way
        if ((n == nfen) .AND. (IEVC == 51)) then
            prod(n)=prod(n)+(150.0*USTARS*USTARS*USTARS)
        endif
    enddo

!     Compute the q2 LHS Matrix and RHS Load Vector

    call KQnmInt(KQnm,Kq)

    Mqa(1) = 0.d0
    Mqb(1) = 1.d0
    Mqc(1) = 0.d0
    LVq(1) = B123*sqrt(bsx_loc*bsx_loc+bsy_loc*bsy_loc)

    coef2 = delt*theta1/HOamb2
!     coef2 = theta1/HOamb2                !rog way
    coef4 = 2.d0*delt
!     coef4 = 2.d0                         !rog way
    coef5 = delt*(1.d0-theta1)/HOamb2
!     coef5 = (1.d0-theta1)/HOamb2         !rog way
    do n=2,nfen-1
        tdiss = 2.d0*diss(n)
        coef1 = 1.d0 + delt*theta2*tdiss
    !     coef1 = 1.d0/delt + theta2*tdiss       !rog way
        Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
    !     Mqa(n) = KQnm(n,1)*coef2               !lumping
        Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
    !     Mqb(n) = (Inm(n,1)+Inm(n,2)+Inm(n,3))*coef1+KQnm(n,2)*coef2  !lumping
        Mqc(n) = Inm(n,3)*coef1+KQnm(n,3)*coef2
    !     Mqc(n) = KQnm(n,3)*coef2                                     !lumping

        coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
    !     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss                      !rog way
        LVq(n) = Inm(n,1)*(coef3*q2prev(n-1)+coef4*prod(n-1)) &
        -KQnm(n,1)*coef5*q2prev(n-1) &
        + Inm(n,2)*(coef3*q2prev(n  )+coef4*prod(n  )) &
        -KQnm(n,2)*coef5*q2prev(n) &
        + Inm(n,3)*(coef3*q2prev(n+1)+coef4*prod(n+1)) &
        -KQnm(n,3)*coef5*q2prev(n+1)
    !     LVq(n) =                 -KQnm(n,1)*coef5*q2prev(n-1)     !lumping
    !     &         + (Inm(n,1)+Inm(n,2)+Inm(n,3))                  !lumping
    !     &         *(coef3*q2prev(n)+coef4*prod(n))     !lumping
    !     &         -KQnm(n,2)*coef5*q2prev(n)       !lumping
    !     &         -KQnm(n,3)*coef5*q2prev(n+1)     !lumping

    enddo

    if(ibc == 0) then
        Mqa(nfen) = 0.d0
        Mqb(nfen) = 1.d0
        Mqc(nfen) = 0.d0
        LVq(nfen) = B123*sqrt(wsx*wsx+wsy*wsy)
    endif

    if(ibc == 1)then
        n=nfen
        tdiss = 2.d0*diss(n)
        coef1 = 1.d0 + delt*theta2*tdiss
    !     coef1 = 1.d0/delt + theta2*tdiss                   !rog way
        Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
    !     Mqa(n) = KQnm(n,1)*coef2                           !lumping
        Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
    !     Mqb(n) = (Inm(n,1)+Inm(n,2))*coef1+KQnm(n,2)*coef2 !lumping
        Mqc(n) = 0.d0

        coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
    !     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss            !rog way
        LVq(n) = Inm(n,1)*(coef3*q2prev(n-1)+coef4*prod(n-1)) &
        -KQnm(n,1)*coef5*q2prev(n-1) &
        + Inm(n,2)*(coef3*q2prev(n  )+coef4*prod(n  )) &
        -KQnm(n,2)*coef5*q2prev(n)
    !     LVq(n) =    -KQnm(n,1)*coef5*q2prev(n-1)   !lumping
    !     &    + (Inm(n,1)+Inm(n,2))*(coef3*q2prev(n)+coef4*prod(n))       !lumping
    !     &         -KQnm(n,2)*coef5*q2prev(n)     !lumping
    endif

!     Solve the system for q2

    CALL TRIDAG(Mqa,Mqb,Mqc,LVq,q2,nfen)

!     Transfer to global array and check for zero or negative values
!     (generally for startup)

    do n = 1,nfen
        if(q2(n) <= 0.)then
            q2(n) = q2min
        endif
        q20(nh,n) = q2(n)
    enddo

!     Compute length scale

    if(il == 2)then
        do n=1,nfen
            sig = (sigma(n)+1.d0)/amb
            rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
            rl2 = 1.d0/(rkap*(H-sig*H+Z0Sw))
            l(nh,n) = 1.d0/(rl1+rl2)

            if(l(nh,n) < 0.001) l(nh,n)=0.001
            if(l(nh,n) > H) l(nh,n)=H
        enddo
    endif

    if(il == 3)then
        do n=1,nfen
            HOamb=H/amb
            zval = (sigma(n)+1.d0)*HOamb-H
            if(sigma(n) <= -0.7d0) l(nh,n)=rkap*(H+zval+Z0B1)
            if(sigma(n) > -0.7d0) l(nh,n)=rkap*(0.15d0*H+Z0B1)
        enddo
    endif

    if(il == 1)then

    !     Compute the q2l LHS Matrix and RHS Load Vector

        Mqa(1) = 0.d0
        Mqb(1) = 1.d0
        Mqc(1) = 0.d0
        LVq(1) = rkap*BSlay*q2(1)

        coef2 = delt*theta1/HOamb2
    !     coef2 = theta1/HOamb2                                     !rog way
        coef5 = delt*(1.d0-theta1)/HOamb2
    !     coef5 = (1.d0-theta1)/HOamb2                              !rog way
        do n=2,nfen-1
            tdiss = wall(n)*diss(n)
            coef1 = 1.d0 + delt*theta2*tdiss
        !     coef1 = 1.d0/delt + theta2*tdiss                        !rog way
            Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
        !     Mqa(n) = KQnm(n,1)*coef2                                !lumping
            Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
        !     Mqb(n) = (Inm(n,1)+Inm(n,2)+Inm(n,3))*coef1+KQnm(n,2)*coef2   !lumping
            Mqc(n) = Inm(n,3)*coef1+KQnm(n,3)*coef2
        !     Mqc(n) = KQnm(n,3)*coef2                                !lumping

            coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
        !     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss                 !rog way
            coef4 = E1*rmlen(n)*delt
        !     coef4 = E1*rmlen(n)                                     !rog way
            LVq(n) = Inm(n,1)*(coef3*q2lprev(n-1)+coef4*prod(n-1)) &
            -KQnm(n,1)*coef5*q2lprev(n-1) &
            + Inm(n,2)*(coef3*q2lprev(n  )+coef4*prod(n  )) &
            -KQnm(n,2)*coef5*q2lprev(n) &
            + Inm(n,3)*(coef3*q2lprev(n+1)+coef4*prod(n+1)) &
            -KQnm(n,3)*coef5*q2lprev(n+1)
        !     LVq(n) =      -KQnm(n,1)*coef5*q2lprev(n-1)  !lumping
        !     &    + (Inm(n,1)+Inm(n,2)+Inm(n,3))                           !lumping
        !     &     *(coef3*q2lprev(n)+coef4*prod(n))                       !lumping
        !     &       -KQnm(n,2)*coef5*q2lprev(n)    !lumping
        !     &     -KQnm(n,3)*coef5*q2lprev(n+1)  !lumping
        
        enddo

        Mqa(nfen) = 0.d0
        Mqb(nfen) = 1.d0
        Mqc(nfen) = 0.d0
        LVq(nfen) = SSlay*rkap*q2(nfen)

    !     Solve the system for q2l

        CALL TRIDAG(Mqa,Mqb,Mqc,LVq,q2l,nfen)

    !     Transfer to global array and check for stability limit

        l(nh,1) = q2l(1)/q2(1)
        do n = 2,nfen-1
            l(nh,n) = q2l(n)/q2(n)
            if(l(nh,n) < 0.) &
            l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
            if(l(nh,n) > H) l(nh,n)=H
            if((iden /= 0) .AND. (BVfreq2(n) > 0.0)) then
                elmax = 0.53D0*sqrt(q2(n))/sqrt(BVfreq2(n))
                if(l(nh,n) > elmax) l(nh,n)=elmax
            endif
        end do
        l(nh,nfen) = q2l(nfen)/q2(nfen)

    ! RJW 02252009 ...
    ! Now we do a little bookeeping, to be sure things stay in
    ! reasonable bounds.
    ! 1) set a lower limit on the turbulent length scale of 10 cm
    ! 2)in shallow water we want to modify the calculation of l
    ! set Hlimit = 1.0 and l_2 = H
    ! If the total water depth is less than 1.0 meter
    ! then we want to set the turbulent length scale to
    ! equal the water depth
        H_crit = 1.0

        do n = 1, nfen

            if (l(nh,n) < 0.10 ) l(nh,n) = 0.10

            if ( H <= H_crit ) then

                l(nh,n) = H/H_crit*l(nh,n) + (1-H/H_crit)*H

            endif

        enddo
    !...

    endif


!     compute eddy viscosity and store variables for output and next step

!   kmd48.33bc need to add the vertical diffusion term for transport
!              the lower limit will be the same as that for the eddy
!              viscosity
    do n = 1,nfen
        EVTOT(n) = Sm(n)*sqrt(q2(n))*l(nh,n)
        NTVTOT(n) = Sh(n)*sqrt(q2(n))*l(nh,n)
        if(EVTOT(n) < EVMIN) EVTOT(n)=EVMIN
        if(EVTOT(n) > 100.0d0) EVTOT(n)=100.0d0
        if(NTVTOT(n) < EVMIN) NTVTOT(n)=EVMIN
    end do

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    return
!***********************************************************************
    END subroutine turb
!***********************************************************************



!***********************************************************************
! Subroutine tridag                                                    *
!                                                                      *
!     -----------------------------------------------------------------*
!     |SOLVER FOR A VECTOR U OF LENGTH N FROM A SET OF LINEAR          *
!     |EQUATIONS THAT CONTAINS A TRIDIAGONAL MATRIX                    *
!     |THE FORM IS                                                     *
!     |                                                                *
!     |   * B1 C1  0 ...               *     * U1 *     * R1 *         *
!     |  *                              *   *      *   *      *        *
!     |  *  A2 B2 C2 ...                *   *  U2  *   *  R2  *        *
!     |  *  ...                         * * * ...  * = * ...  *        *
!     |  *           ... An-1 Bn-1 Cn-1 *   * Un-1 *   * Rn-1 *        *
!     |  *                              *   *      *   *      *        *
!     |   *                0   An  Bn  *     * Un *     * Rn *         *
!     |                                                                *
!     |A, B, C, U ARE ARRAYS.                                          *
!     -----------------------------------------------------------------*
!                                                                      *
!    Adapted from Numerical Recipes chapter 2                          *
!***********************************************************************

    SUBROUTINE TRIDAG(A,B,C,R,U,N)
    USE GLOBAL, ONLY : ScreenUnit
    USE GLOBAL_3DVS, ONLY : SZ, setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
    IMPLICIT NONE
    INTEGER :: J, N

    REAL(SZ) :: A(N),B(N),C(N),R(N),U(N)
    REAL(SZ) :: BET,GAM(N)

    call setMessageSource("tridag")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!     check for zero elements on diagonal

    DO J=1,N
        if(B(j) == 0.) then
            write(screenunit,*) 'Problem in Tridag Solver.  ', &
            'B array value in row ',j,' = 0'
            stop
        endif
    end do
    BET = B(1)
    U(1) = R(1)/BET
    DO J = 2,N
        GAM(J) = C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        if (BET == 0) then
            write(screenunit,*) ' Problem in Tridag Solver.  ', &
            ' BET  = 0.  Solver failed.'
            stop
        endif
        U(J)=(R(J)-A(J)*U(J-1))/BET
    END DO
    DO J = N-1,1,-1
        U(J) = U(J) - GAM(J+1)*U(J+1)
    END DO

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!***********************************************************************
    END SUBROUTINE TRIDAG
!***********************************************************************



!****************************************************************************************
!   Subroutine to interpolate baroclinic pressure (BCP) to a specified sigma value      *
!   (SigmaNN) given an initial guess of which sigma level is closest to the specified   *
!   value.                                                                              *
!                                                                                       *
!                                    R.L.  5/04/01                                      *
!                                    R.L.  5.19/03                                      *
!                                    R.L.  5/08/05                                      *
!****************************************************************************************

!      SUBROUTINE ZSURFBUOY(SigmaNN,BCPressNN,NN,J)
!      USE GLOBAL, ONLY : ScreenUnit
!      USE GLOBAL_3DVS
!      REAL(SZ) :: BCPressNN
!      REAL(SZ) :: SigmaNN     !Sigma value of a neighbor node

!      IF(SigmaNN.LE.1.0001*b) THEN !if into ground then skip
!         SigBelo=-999
!         SigAbov=-999
!         BCPressNN=-999.
!         GOTO 100
!      ENDIF
!      IF((SigmaNN.GT.1.0001*b).AND.(SigmaNN.LE.b)) THEN !at bottom then use bottom
!         LBelo=1
!         BCPressNN=BCP(NN,LBelo)
!         SigBelo=b
!         SigAbov=b
!         GOTO 100
!      ENDIF
!      IF(SigmaNN.GE.a) THEN     !into air use surface
!         LAbov=NFEN
!         BCPressNN=BCP(NN,LAbov)
!         SigBelo=a
!         SigAbov=a
!         GOTO 100
!      ENDIF

!      LTry=J                    !start search for SIGABOV and SIGBELO
!      SigTry=Sigma(LTry)
!      IF(SigmaNN.GT.SigTry) THEN !too low
!         SigBelo=SigTry         !SIGBELO may = SIGTRY
!         LBelo=LTry
!         LTry=LTry+1            !look at next level higher
! 90      SigTry=Sigma(LTry)
!         IF(SigmaNN.GT.SigTry) THEN !still too low
!            SigBelo=SigTry
!            LBelo=LTry
!            LTry=LTry+1
!            GOTO 90
!         ENDIF
!         SigAbov=SigTry         !found upper bracketing sigma
!         LAbov=LTry
!         GOTO 99                !go interpolate
!      ENDIF
!      IF(SigmaNN.LE.SigTry) THEN !to high
!         SigAbov=SigTry         !SIGABOV may = SIGTRY
!         LAbov=LTry
!         LTry=LTry-1            !look at next level lower
! 91      SigTry=Sigma(LTry)
!         IF(SigmaNN.LE.SigTry) THEN !still too high
!            SigAbov=SigTry
!            LAbov=LTry
!            LTry=LTry-1
!            GOTO 91
!         ENDIF
!         SigBelo=SigTry         !found lower bracketing sigma
!         LBelo=LTry
!      ENDIF

! 99   BCPressNN=(BCP(NN,LAbov)-BCP(NN,LBelo)) !interpolation
!     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + BCP(NN,LBelo)

! 100  CONTINUE

!      RETURN
!      END


!***************************************************************************
!   Subroutine to write out 3D Hot Start info

!                                    R.L.  8/16/05

!     jgf49.44: This subroutine is only called if a binary hotstart
!     file was specified.
!***************************************************************************

!  kmd48.33bc changed to follow the global IO format
    SUBROUTINE HSTART3D_OUT(IT)
    USE SIZES
    USE GLOBAL
    USE GLOBAL_3DVS
    USE GLOBAL_IO
#ifdef CMPI
    USE MESSENGER
#endif
#ifdef ADCNETCDF
    USE NETCDFIO, ONLY : initNetCDFHotstart3D, writeNetCDFHotstart3D, &
    writeNetCDFHotstart3DVar
#endif

    IMPLICIT NONE
    INTEGER, intent(in) :: IT ! current time step
    INTEGER :: NH, I, J ! loop counters
    INTEGER :: N,M

    type(OutputDataDescript_t),SAVE :: Duudescript
    type(OutputDataDescript_t),SAVE :: Duvdescript
    type(OutputDataDescript_t),SAVE :: Dvvdescript
    type(OutputDataDescript_t),SAVE :: Uudescript
    type(OutputDataDescript_t),SAVE :: Vvdescript
    type(OutputDataDescript_t),SAVE :: Bsxdescript
    type(OutputDataDescript_t),SAVE :: Bsydescript
!     jgf49.49.02: Added variables for 3D fields
    type(OutputDataDescript_t), save :: SigTDescript
    type(OutputDataDescript_t), save :: SalDescript
    type(OutputDataDescript_t), save :: TempDescript
    type(OutputDataDescript_t), save :: RealQDescript
    type(OutputDataDescript_t), save :: ImaginaryQDescript
    type(OutputDataDescript_t), save :: WZDescript
    type(OutputDataDescript_t), save :: Q20Descript
    type(OutputDataDescript_t), save :: LDescript
    REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rp(:,:) !real part subdomain data
    REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ip(:,:) !imaginary part subdom dat
    REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rp_g(:,:) !real part fulldomain
    REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ip_g(:,:) !imaginary part fulldom
    REAL(SZ), ALLOCATABLE, TARGET, SAVE :: layer(:) ! horiz layer in 3D
    INTEGER :: numHotstartWrites
    INTEGER :: nextLun

    LOGICAL, SAVE :: FirstCall3D = .TRUE. 
    call setMessageSource("hstart3D_out")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

!     kmd48.33bc write out results in global IO format
    IF (FirstCall3D.eqv. .TRUE. ) THEN

        ALLOCATE(rp(NP,NFEN),ip(NP,NFEN))

        IF ((MNPROC > 1) .AND. (myProc == 0) .AND. &
        (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
            ALLOCATE(layer(np_g))
        ENDIF
        IF ((MNPROC == 1) &
         .OR. (WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. )) THEN
            ALLOCATE(layer(np))
        ENDIF

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(DUU_g(NP_G))
        END IF
        Duudescript % specifier            =  NHSTAR
        Duudescript % initial_value        =  0.0
        Duudescript % num_items_per_record =  1
        Duudescript % num_fd_records       =  NP_G
        Duudescript % num_records_this     =  NP
        Duudescript % imap                 => nodes_lg
        Duudescript % array                => DUU1
        Duudescript % array_g              => DUU_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(DUV_g(NP_G))
        END IF
        Duvdescript % specifier            =  NHSTAR
        Duvdescript % initial_value        =  0.0
        Duvdescript % num_items_per_record =  1
        Duvdescript % num_fd_records       =  NP_G
        Duvdescript % num_records_this     =  NP
        Duvdescript % imap                 => nodes_lg
        Duvdescript % array                => DUV1
        Duvdescript % array_g              => DUV_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(DVV_g(NP_G))
        END IF
        Dvvdescript % specifier            =  NHSTAR
        Dvvdescript % initial_value        =  0.0
        Dvvdescript % num_items_per_record =  1
        Dvvdescript % num_fd_records       =  NP_G
        Dvvdescript % num_records_this     =  NP
        Dvvdescript % imap                 => nodes_lg
        Dvvdescript % array                => DVV1
        Dvvdescript % array_g              => DVV_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(UU_g(NP_G))
        END IF
        Uudescript % specifier            =  NHSTAR
        Uudescript % initial_value        =  0.0
        Uudescript % num_items_per_record =  1
        Uudescript % num_fd_records       =  NP_G
        Uudescript % num_records_this     =  NP
        Uudescript % imap                 => nodes_lg
        Uudescript % array                => UU2
        Uudescript % array_g              => UU_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(VV_g(NP_G))
        END IF
        Vvdescript % specifier            =  NHSTAR
        Vvdescript % initial_value        =  0.0
        Vvdescript % num_items_per_record =  1
        Vvdescript % num_fd_records       =  NP_G
        Vvdescript % num_records_this     =  NP
        Vvdescript % imap                 => nodes_lg
        Vvdescript % array                => VV2
        Vvdescript % array_g              => VV_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(BSX_g(NP_G))
        END IF
        Bsxdescript % specifier            =  NHSTAR
        Bsxdescript % initial_value        =  0.0
        Bsxdescript % num_items_per_record =  1
        Bsxdescript % num_fd_records       =  NP_G
        Bsxdescript % num_records_this     =  NP
        Bsxdescript % imap                 => nodes_lg
        Bsxdescript % array                => BSX1
        Bsxdescript % array_g              => BSX_g

        IF ((MNPROC > 1) .AND. (MyProc == 0) ) THEN
            Allocate(BSY_g(NP_G))
        END IF
        Bsydescript % specifier            =  NHSTAR
        Bsydescript % initial_value        =  0.0
        Bsydescript % num_items_per_record =  1
        Bsydescript % num_fd_records       =  NP_G
        Bsydescript % num_records_this     =  NP
        Bsydescript % imap                 => nodes_lg
        Bsydescript % array                => BSY1
        Bsydescript % array_g              => BSY_g

        SigTDescript % specifier            =  I3DGD
        SigTDescript % initial_value        =  0.0
        SigTDescript % num_items_per_record =  NFEN
        SigTDescript % num_fd_records       =  NP_G
        SigTDescript % num_records_this     =  NP
        SigTDescript % imap                 => nodes_lg
        SigTDescript % array2D              => SigT
        SigTDescript % array2D_g            => SigT_g
        SigTDescript % considerWetDry       = .FALSE. 
        SigTDescript % alternate_value      = -99999.0
        SigTDescript % field_name           = 'SigmaT'

        SalDescript % specifier            =  I3DGD
        SalDescript % initial_value        =  0.0
        SalDescript % num_items_per_record =  NFEN
        SalDescript % num_fd_records       =  NP_G
        SalDescript % num_records_this     =  NP
        SalDescript % imap                 => nodes_lg
        SalDescript % array2D              => Sal
        SalDescript % array2D_g            => Sal_g
        SalDescript % considerWetDry       = .FALSE. 
        SalDescript % alternate_value      = -99999.0
        SalDescript % field_name           = 'Salinity'

        TempDescript % specifier            =  I3DGD
        TempDescript % initial_value        =  0.0
        TempDescript % num_items_per_record =  NFEN
        TempDescript % num_fd_records       =  NP_G
        TempDescript % num_records_this     =  NP
        TempDescript % imap                 => nodes_lg
        TempDescript % array2D              => Temp
        TempDescript % array2D_g            => Temp_g
        TempDescript % considerWetDry       = .FALSE. 
        TempDescript % alternate_value      = -99999.0
        TempDescript % field_name           = 'Temperature'

        RealQDescript % specifier            =  I3DGV
        RealQdescript % initial_value        =  0.0
        RealQDescript % num_items_per_record =  NFEN
        RealQDescript % num_fd_records       =  NP_G
        RealQDescript % num_records_this     =  NP
        RealQDescript % imap                 => nodes_lg
        RealQDescript % array2D              => rp
        RealQDescript % array2D_g            => rp_g
        RealQDescript % considerWetDry       = .FALSE. 
        RealQDescript % alternate_value      = -99999.0
        RealQDescript % field_name           = 'u-vel3D'

        ImaginaryQDescript % specifier            =  I3DGV
        ImaginaryQDescript % initial_value        =  0.0
        ImaginaryQDescript % num_items_per_record =  NFEN
        ImaginaryQDescript % num_fd_records       =  NP_G
        ImaginaryQDescript % num_records_this     =  NP
        ImaginaryQDescript % imap                 => nodes_lg
        ImaginaryQDescript % array2D              => ip
        ImaginaryQDescript % array2D_g            => ip_g
        ImaginaryQDescript % considerWetDry       = .FALSE. 
        ImaginaryQDescript % alternate_value      = -99999.0
        ImaginaryQDescript % field_name           = 'v-vel3D'

        WZDescript % specifier            =  I3DGV
        WZDescript % initial_value        =  0.0
        WZDescript % num_items_per_record =  NFEN
        WZDescript % num_fd_records       =  NP_G
        WZDescript % num_records_this     =  NP
        WZDescript % imap                 => nodes_lg
        WZDescript % array2D              => WZ
        WZDescript % array2D_g            => WZ_g
        WZDescript % considerWetDry       = .FALSE. 
        WZDescript % alternate_value      = -99999.0
        WZDescript % field_name           = 'w-vel3D'

        Q20Descript % specifier            =  I3DGT
        Q20Descript % initial_value        =  0.0
        Q20Descript % num_items_per_record =  NFEN
        Q20Descript % num_fd_records       =  NP_G
        Q20Descript % num_records_this     =  NP
        Q20Descript % imap                 => nodes_lg
        Q20Descript % array2D              => q20
        Q20Descript % array2D_g            => q20_g
        Q20Descript % considerWetDry       = .FALSE. 
        Q20Descript % alternate_value      = -99999.0
        Q20Descript % field_name           = 'q20'

        LDescript % specifier            =  I3DGT
        LDescript % initial_value        =  0.0
        LDescript % num_items_per_record =  NFEN
        LDescript % num_fd_records       =  NP_G
        LDescript % num_records_this     =  NP
        LDescript % imap                 => nodes_lg
        LDescript % array2D              => l
        LDescript % array2D_g            => l_g
        LDescript % considerWetDry       = .FALSE. 
        LDescript % alternate_value      = -99999.0
        LDescript % field_name           = 'l'


        IF ((NHSTAR == 3) .OR. (NHSTAR == 367) .OR. (NHSTAR == 368) .OR. &
        (NHSTAR == 5) .OR. (NHSTAR == 567) .OR. (NHSTAR == 568)) THEN
#ifdef ADCNETCDF
        ! in parallel, we don't need to create the hotstart file, it was
        ! created by adcprep.
            IF (MNPROC == 1) THEN
                CALL initNetCDFHotstart3D(hss%lun, NHSTAR)
                numHotstartWrites = (NT-IT)/NHSINC
                IF (numHotstartWrites >= 1) THEN
                    IF (hss%lun == 67) THEN
                        nextLun = 68
                    ELSE
                        nextLun = 67
                    ENDIF
                    CALL initNetCDFHotstart3D(nextLun, NHSTAR)
                ENDIF
            ENDIF
#endif
        ENDIF

        FirstCall3d = .FALSE. 
    END IF

#ifdef CMPI
    CALL collectFullDomainArray(Duudescript,packOne,unpackOne)
    CALL collectFullDomainArray(Duvdescript,packOne,unpackOne)
    CALL collectFullDomainArray(Dvvdescript,packOne,unpackOne)
    CALL collectFullDomainArray(Uudescript,packOne,unpackOne)
    CALL collectFullDomainArray(Vvdescript,packOne,unpackOne)
    CALL collectFullDomainArray(Bsxdescript,packOne,unpackOne)
    CALL collectFullDomainArray(Bsydescript,packOne,unpackOne)
#endif

    SELECT CASE (NHSTAR)

    CASE(-1,1,67,68) ! non-portable binary  !tcm v51.26 mod for time-stamped

    IF((MYPROC == 0) .OR. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. )) THEN
        WRITE(hss % lun,REC=IHOTSTP) IDen  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DSD  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DSDRec  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DSV  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DSVRec  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DST  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DSTRec  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DGD  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DGDRec  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DGV  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DGVRec  ;  IHOTSTP = IHOTSTP + 1

        WRITE(hss % lun,REC=IHOTSTP) N3DGT  ;  IHOTSTP = IHOTSTP + 1
        WRITE(hss % lun,REC=IHOTSTP) I3DGTRec  ;  IHOTSTP = IHOTSTP + 1
    ENDIF

    IF ((MNPROC > 1) .AND. (myProc == 0) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        CALL binaryWrite2D(DUU_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(DUV_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(DVV_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(UU_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(VV_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(BSX_g, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(BSY_g, np_g, hss%lun, ihotstp)
    ENDIF
    IF ((MNPROC == 1) .OR. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. )) THEN
        CALL binaryWrite2D(DUU, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(DUV, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(DVV, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(UU, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(VV, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(BSX, np_g, hss%lun, ihotstp)
        CALL binaryWrite2D(BSY, np_g, hss%lun, ihotstp)
    ENDIF

    CASE(3,367,368,5,567,568) ! netcdf classic; netcdf3 or netcdf4
#ifdef ADCNETCDF
    IF (myProc == 0) THEN
        CALL writeNetCDFHotstart3D(hss%lun,DUUDescript, &
        DUVDescript, DVVDescript, UUDescript, VVDescript, &
        BSXDescript, BSYDescript)
    ENDIF
#endif
    CASE DEFAULT

    END SELECT

!     3 D  V E L O C I T Y
    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        IF (myProc == 0) THEN
            ALLOCATE(rp_g(NP_G,NFEN),ip_g(NP_G,NFEN),WZ_g(NP_G,NFEN))
            RealQDescript % array2D_g => rp_g
            ImaginaryQDescript % array2D_g => ip_g
            WZDescript % array2D_g => WZ_g
        ENDIF
        rp = real(q)
        CALL collectFullDomainArray(RealQDescript, &
        packNPbyM, unpackNPbyM)
        ip = aimag(q)
        CALL collectFullDomainArray(ImaginaryQDescript, &
        packNPbyM, unpackNPbyM)
        CALL collectFullDomainArray(WZDescript, &
        packNPbyM, unpackNPbyM)
    ENDIF

    SELECT CASE(NHSTAR)
    CASE(-1,1,67,68) ! non-portable binary   !tcm v51.26 mod for time-stamped
    IF ((MNPROC == 1) .OR. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. )) THEN
        DO N=1,NFEN
            layer(:) = real(q(:,N))
            CALL binaryWrite2D(layer, np, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            layer(:) = aimag(q(:,N))
            CALL binaryWrite2D(layer, np, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            CALL binaryWrite2D(WZ(:,N), np, hss%lun, ihotstp)
        END DO
    ENDIF
    IF ((MNPROC > 1) .AND. (myProc == 0) &
     .AND. (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        DO N=1,NFEN
            CALL binaryWrite2D(rp_g(:,N), np_g, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            CALL binaryWrite2D(ip_g(:,N), np_g, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            CALL binaryWrite2D(WZ_g(:,N), np_g, hss%lun, ihotstp)
        END DO
    ENDIF
    CASE(3,367,368,5,567,568)
#ifdef ADCNETCDF
    IF ((MNPROC == 1) .OR. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. )) THEN
        rp = real(q)
        ip = aimag(q)
    ENDIF
    IF (myProc == 0) THEN
        CALL writeNetCDFHotstart3DVar(hss%lun, RealQDescript)
        CALL writeNetCDFHotstart3DVar(hss%lun,ImaginaryQDescript)
        CALL writeNetCDFHotstart3DVar(hss%lun, WZDescript)
    ENDIF
#endif
    CASE DEFAULT

    END SELECT

    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. ) .AND. &
    (myProc == 0)) THEN
        DEALLOCATE(rp_g,ip_g,WZ_g)
    ENDIF

!     3 D  T U R B U L E N C E
    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        IF (myProc == 0) THEN
            ALLOCATE(q20_g(NP_G,NFEN),l_g(NP_G,NFEN))
            Q20Descript % array2D_g => q20_g
            LDescript % array2D_g => l_g
        ENDIF
        CALL collectFullDomainArray(Q20Descript, &
        packNPbyM, unpackNPbyM)
        CALL collectFullDomainArray(LDescript, &
        packNPbyM, unpackNPbyM)
    ENDIF
    SELECT CASE(NHSTAR)
    CASE(-1,1,67,68) ! non-portable binary  !tcm v51.26 mod for time-stamped hotstarts (-1)
    IF ((MNPROC == 1) .OR. (WRITE_LOCAL_FILES.eqv. .TRUE. )) THEN
        DO N=1,NFEN
            CALL binaryWrite2D(q20(:,N), np, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            CALL binaryWrite2D(l(:,N), np, hss%lun, ihotstp)
        END DO
    ENDIF
    IF ((MNPROC > 1) .AND. (myProc == 0) &
     .AND. (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        DO N=1,NFEN
            CALL binaryWrite2D(q20_g(:,N), np_g, hss%lun, ihotstp)
        END DO
        DO N=1,NFEN
            CALL binaryWrite2D(l_g(:,N), np_g, hss%lun, ihotstp)
        END DO
    ENDIF
    CASE(3,367,368,5,567,568)
#ifdef ADCNETCDF
    IF (myProc == 0) THEN
        CALL writeNetCDFHotstart3DVar(hss%lun, Q20Descript)
        CALL writeNetCDFHotstart3DVar(hss%lun, LDescript)
    ENDIF
#endif
    CASE DEFAULT

    END SELECT
    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. ) .AND. &
    (myProc == 0)) THEN
        DEALLOCATE(q20_g,l_g)
    ENDIF

!     3 D  D E N S I T Y
    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        IF (myProc == 0) THEN
            IF (ABS(IDEN) == 1) THEN
                ALLOCATE(SigT_g(NP_G,NFEN))
                SigTDescript % array2D_g => SigT_g
            ENDIF
            IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
                ALLOCATE(Sal_g(NP_G,NFEN))
                SalDescript % array2D_g => Sal_g
            ENDIF
            IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
                ALLOCATE(Temp_g(NP_G,NFEN))
                TempDescript % array2D_g => Temp_g
            ENDIF
        ENDIF
        IF (ABS(IDEN) == 1) THEN
            CALL collectFullDomainArray(SigTDescript, &
            packNPbyM, unpackNPbyM)
        ENDIF
        IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
            CALL collectFullDomainArray(SalDescript, &
            packNPbyM, unpackNPbyM)
        ENDIF
        IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
            CALL collectFullDomainArray(TempDescript, &
            packNPbyM, unpackNPbyM)
        ENDIF
    ENDIF

    SELECT CASE(NHSTAR)
    CASE(-1,1,67,68) ! non-portable binary  !tcm v51.26 mod for time-stamped
    IF ((MNPROC == 1) .OR. (WRITE_LOCAL_FILES.eqv. .TRUE. )) THEN
        IF (ABS(IDEN) == 1) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(SigT(:,N), np, hss%lun, ihotstp)
            END DO
        ENDIF
        IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(Sal(:,N), np, hss%lun, ihotstp)
            END DO
        ENDIF
        IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(Temp(:,N), np, hss%lun, ihotstp)
            END DO
        ENDIF
    ENDIF
    IF ((MNPROC > 1) .AND. (myProc == 0) &
     .AND. (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        IF (ABS(IDEN) == 1) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(SigT_g(:,N), np_g, hss%lun, &
                ihotstp)
            END DO
        ENDIF
        IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(Sal_g(:,N), np_g, hss%lun, &
                ihotstp)
            END DO
        ENDIF
        IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
            DO N=1,NFEN
                CALL binaryWrite2D(Temp_g(:,N), np_g, hss%lun, &
                ihotstp)
            END DO
        ENDIF
    ENDIF
    CASE(3,367,368,5,567,568)
#ifdef ADCNETCDF
    IF (myProc == 0) THEN
        IF (ABS(IDEN) == 1) THEN
            CALL writeNetCDFHotstart3DVar(hss%lun, SigTDescript)
        ENDIF
        IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
            CALL writeNetCDFHotstart3DVar(hss%lun, SalDescript)
        ENDIF
        IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
            CALL writeNetCDFHotstart3DVar(hss%lun, TempDescript)
        ENDIF
    ENDIF
#endif
    END SELECT

    IF ((MNPROC > 1) .AND. &
    (WRITE_LOCAL_HOT_START_FILES.eqv. .FALSE. )) THEN
        IF (myProc == 0) THEN
            IF (ABS(IDEN) == 1) THEN
                DEALLOCATE(SigT_g)
            ENDIF
            IF ((ABS(IDEN) == 2) .OR. (ABS(IDEN) == 4)) THEN
                DEALLOCATE(Sal_g)
            ENDIF
            IF ((ABS(IDEN) == 3) .OR. (ABS(IDEN) == 4)) THEN
                DEALLOCATE(Temp_g)
            ENDIF
        ENDIF
    ENDIF
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!-----------------------------------------------------------------------
    END SUBROUTINE HSTART3D_OUT
!-----------------------------------------------------------------------
!      DO NH=1,NP
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) DUU(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) DUV(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) DVV(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) UU(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) VV(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) BSX(NH)
!         IHOTSTP=IHOTSTP+1
!         WRITE(hss % lun,REC=IHOTSTP) BSY(NH)
!      ENDDO

!      DO NH=1,NP
!         DO N=1,NFEN
!            IHOTSTP=IHOTSTP+1
!            WRITE(hss % lun,REC=IHOTSTP) REAL(Q(NH,N))
!            IHOTSTP=IHOTSTP+1
!            WRITE(hss % lun,REC=IHOTSTP) AIMAG(Q(NH,N))
!            IHOTSTP=IHOTSTP+1
!            WRITE(hss % lun,REC=IHOTSTP) WZ(NH,N)
!            IHOTSTP=IHOTSTP+1
!            WRITE(hss % lun,REC=IHOTSTP) q20(NH,N)
!            IHOTSTP=IHOTSTP+1
!            WRITE(hss % lun,REC=IHOTSTP) l(NH,N)
!            IF(ABS(IDEN).EQ.1) THEN
!               IHOTSTP=IHOTSTP+1
!               WRITE(hss % lun,REC=IHOTSTP) SigT(NH,N)
!            ENDIF
!            IF(ABS(IDEN).EQ.2) THEN
!               IHOTSTP=IHOTSTP+1
!               WRITE(hss % lun,REC=IHOTSTP) Sal(NH,N)
!            ENDIF
!            IF(ABS(IDEN).EQ.3) THEN
!               IHOTSTP=IHOTSTP+1
!               WRITE(hss % lun,REC=IHOTSTP) Temp(NH,N)
!            ENDIF
!            IF(ABS(IDEN).EQ.4) THEN
!               IHOTSTP=IHOTSTP+1
!               WRITE(hss % lun,REC=IHOTSTP) Sal(NH,N)
!               IHOTSTP=IHOTSTP+1
!               WRITE(hss % lun,REC=IHOTSTP) Temp(NH,N)
!            ENDIF
!         ENDDO
!      ENDDO



!****************************************************************************
!   Subroutine to compute 3D SigmaT fields from 3D salinity
!   and/or temperature fields

!                                    R.L.  6/22/05
!                                    K.D.  2/03/06 updated
!*****************************************************************************


    SUBROUTINE CALC_SIGMAT_3D ()

! endra: Added in information needed for the calculation of sigma-t.
!        Currently, I have added only the equation of state from
!        Cushman-Roisin book entitled "Introduction to Geophysical
!        Fluid Dynamics". Later, we can update this to include more
!        complex equations of state.
!   kmd48.33bc need to add variables for new equation of state
    USE GLOBAL_3DVS, ONLY : SAL, TEMP, SIGT, IDEN, NFEN, SZ, &
    SIGMA, EQNSTATE, A, B, setMessageSource, &
    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
! jgfdebug      USE GLOBAL, ONLY : NP,RHOWAT0, ScreenUnit
    USE GLOBAL, ONLY : RHOWAT0, ETA2, IFNLFA, G
    USE MESH, ONLY : NP, DP
    IMPLICIT NONE
    INTEGER :: NH,N
    REAL(SZ) :: T_alpha=0.00017d0, S_beta=0.00076d0
!      REAL(SZ) :: T_ref=10.0d0, S_ref=35.0d0
    REAL(SZ) :: T_ref, S_ref ! changed to new values
! endra: Parameters for the 2nd equation of state
    REAL(SZ), PARAMETER :: Q1 = 196637339.d0
    REAL(SZ), PARAMETER :: Q2 = 196668.928d0
    REAL(SZ), PARAMETER :: a1 = 1446045.44d0
    REAL(SZ), PARAMETER :: a2 = -10769.1980d0
    REAL(SZ), PARAMETER :: a3 = 85.2955498d0
    REAL(SZ), PARAMETER :: b1 = 579854.265d0
    REAL(SZ), PARAMETER :: b2 = -1477.31528d0
    REAL(SZ), PARAMETER :: b3 = 419.489109d0
    REAL(SZ), PARAMETER :: c1 = 2001.22349d0
    REAL(SZ), PARAMETER :: c2 = 0.0213851025d0
    REAL(SZ), PARAMETER :: c3 = 0.968784141d0
    REAL(SZ), PARAMETER :: c4 = -0.00654785602d0
    REAL(SZ), PARAMETER :: c5 = -0.00000250468726d0
    REAL(SZ), PARAMETER :: c6 = 0.0000000628902345d0
    REAL(SZ), PARAMETER :: c7 = 0.00000000282414763d0
    REAL(SZ), PARAMETER :: d1 = 1433.02205d0
    REAL(SZ), PARAMETER :: d2 = -9.09231525d0
    REAL(SZ), PARAMETER :: d3 = 0.0791654429d0
    REAL(SZ), PARAMETER :: d4 = 0.0000398630534d0
    REAL(SZ), PARAMETER :: e1 = 417.831720d0
    REAL(SZ), PARAMETER :: e2 = -1.87581316d0
    REAL(SZ), PARAMETER :: e3 = -0.0000387902837d0
    REAL(SZ), PARAMETER :: f1 = 1.00765828d0
    REAL(SZ), PARAMETER :: f2 = 0.000312912597d0
    REAL(SZ), PARAMETER :: T_max=40.d0, T_min=0.d0
    REAL(SZ), PARAMETER :: S_max=42.d0, S_min=0.d0
! endra: Variables for the 2nd equation of state
    REAL(SZ) :: P1, P2
    REAL(SZ) :: p(NP), RHO(NP,NFEN)
    REAL(SZ) :: R_ref
    REAL(SZ) :: Hs, HsOAMB, Zk
    REAL(SZ) :: Z(NP,NFEN)
! endra: Parameters for the 3rd equation of state
    REAL(SZ), PARAMETER :: c00=999.842594d0, c11=0.06793952d0
    REAL(SZ), PARAMETER :: c12=-0.00909529d0, c13=0.0001001685d0
    REAL(SZ), PARAMETER :: c14=-0.000001120083d0, c15=0.000000006536332d0
    REAL(SZ), PARAMETER :: c21=0.824493d0, c22=-0.0040899d0
    REAL(SZ), PARAMETER :: c23=0.000076438d0, c24=0.00000082467d0
    REAL(SZ), PARAMETER :: c25=0.0000000053875d0, c31=-0.00572466d0
    REAL(SZ), PARAMETER :: c32=0.00010227d0, c33=-0.0000016546d0
    REAL(SZ), PARAMETER :: c41=0.00048314d0
    REAL(SZ), PARAMETER :: d00=190925.6d0, d11=2098.925d0
    REAL(SZ), PARAMETER :: d12=-30.41638d0, d13=-0.01852732d0
    REAL(SZ), PARAMETER :: d14=-0.0001361629d0, d21=1044.077d0
    REAL(SZ), PARAMETER :: d22=-65.00517d0, d23=1.553190d0
    REAL(SZ), PARAMETER :: d24=0.002326469d0, d31=-55.87545d0
    REAL(SZ), PARAMETER :: d32=7.390729d0, d33=-0.1909078d0
    REAL(SZ), PARAMETER :: d41=-4.721788d0, d42=-0.1028859d0
    REAL(SZ), PARAMETER :: d43=0.002512549d0, d44=0.000005939910d0
    REAL(SZ), PARAMETER :: d51=0.1571896d0, d52=0.002598241d0
    REAL(SZ), PARAMETER :: d53=0.00007267926d0, d61=-0.02042967d0
    REAL(SZ), PARAMETER :: d71=0.0001045941d0, d72=-0.000000005782165d0
    REAL(SZ), PARAMETER :: d73=0.000001296821d0, d81=-0.000002595994d0
    REAL(SZ), PARAMETER :: d82=-0.00000001248266d0, d83=-0.00000003508914d0
! endra: Variables for the 3rd equation of state
    REAL(SZ) :: Kappa

    call setMessageSource("calc_sigmat_3D")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif


! endra: Equation of state from Cushman-Roisin book entitled "Introduction to
!        Geophysical Fluid Dynamics".
    IF (Eqnstate == 1) THEN
        T_ref=10.d0
        S_ref=35.d0
        IF(ABS(IDEN) == 2) THEN
            DO NH=1,NP
                DO N=1,NFEN
                    SIGT(NH,N) = 1028.d0*(1.d0+S_beta*(SAL(NH,N)- &
                    S_ref))-RHOWAT0
                ENDDO
            ENDDO
        ELSEIF(ABS(IDEN) == 3) THEN
            DO NH=1,NP
                DO N=1,NFEN
                    SIGT(NH,N) = 1028.d0*(1.d0-T_alpha*(TEMP(NH,N)- &
                    T_ref))-RHOWAT0
                ENDDO
            ENDDO
        ELSEIF(ABS(IDEN) == 4) THEN
            DO NH=1,NP
                DO N=1,NFEN
                    SIGT(NH,N) = 1028.d0*(1.d0-T_alpha*(TEMP(NH,N)- &
                    T_ref)+S_beta*(SAL(NH,N)-S_ref))- &
                    RHOWAT0
                ENDDO
            ENDDO
        ENDIF
    ELSE IF (Eqnstate == 2) THEN
        T_ref=10.d0
        S_ref=32.d0
    !  kmd48.33bc Equation of state from McDougall et al. (2003) - only used
    !        with both salinity and temperature. This is one of the equation
    !        of states used in QUODDY
        IF (ABS(IDEN) == 4) THEN
        ! Look at the values for Temperature and Salinity to see if still
        ! inside the range of the equation of state
            DO NH=1,NP
                DO N=1,NFEN
                    IF (TEMP(NH,N) < T_min) THEN
                        TEMP(NH,N) = T_min
                    END IF
                    IF (TEMP(NH,N) > T_max) THEN
                        TEMP(NH,N) = T_max
                    END IF
                    IF (SAL(NH,N) < S_min) THEN
                        SAL(NH,N) = S_min
                    END IF
                    IF (SAL(NH,N) > S_max) THEN
                        SAL(NH,N) = S_max
                    END IF
                END DO
            END DO
        ! First, determine the reference density
            P1 = Q1 + (a1 + (a2 + a3*T_ref)*T_ref)* &
            T_ref + (b1 + b2*T_ref + &
            b3*S_ref)*S_ref
            P2 = Q2 + (d1 + (d2 + (d3 + d4*T_ref)* &
            T_ref)*T_ref)*T_ref + &
            (e1 + (f1 + f2*T_ref*T_ref)* &
            SQRT(S_ref) + (e2 + e3*T_ref* &
            T_ref)*T_ref)*S_ref
            R_ref = P1/P2
        ! Second, determine the depths of the layers
            DO NH=1,NP
                Hs = DP(NH) + IFNLFA*ETA2(NH)
                HsOAMB = Hs/(A-B)
                DO N=1, NFEN
                    Zk = HsOAMB * (Sigma(n)-B)-DP(NH)
                    Z(NH,N) = Zk
                END DO
            END DO
        ! Third, calculate the new density
            DO NH=1,NP
                DO N=NFEN,1,-1
                    IF (N == NFEN) THEN
                        p(NH) = 0.d0
                        P1 = Q1 + (a1 + (a2 + a3*TEMP(NH,N))*TEMP(NH,N))* &
                        TEMP(NH,N) + (b1 + b2*TEMP(NH,N) + &
                        b3*SAL(NH,N))*SAL(NH,N)
                        P2 = Q2 + (d1 + (d2 + (d3 + d4*TEMP(NH,N))* &
                        TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) + &
                        (e1 + (f1 + f2*TEMP(NH,N)*TEMP(NH,N))* &
                        SQRT(SAL(NH,N)) + (e2 + e3*TEMP(NH,N)* &
                        TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)
                        RHO(NH,N) = P1/P2
                    ELSE IF (N /= NFEN) THEN
                        p(NH) = p(NH) + g*RHO(NH,N+1)*(Z(NH,N+1)- &
                        Z(NH,N))*0.0001d0
                        P1 = Q1 + (a1 + (a2 + a3*TEMP(NH,N))*TEMP(NH,N))* &
                        TEMP(NH,N) + (b1 + b2*TEMP(NH,N) + &
                        b3*SAL(NH,N))*SAL(NH,N) + (c1 + c2* &
                        TEMP(NH,N)*TEMP(NH,N) + c3*SAL(NH,N) + &
                        (c4 + ((c5+c6*TEMP(NH,N))*TEMP(NH,N) + &
                        c7*p(NH))*TEMP(NH,N))*p(NH))*p(NH)
                        P2 = Q2 + (d1 + (d2 + (d3 + d4*TEMP(NH,N))* &
                        TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) + &
                        (e1 + (f1 + f2*TEMP(NH,N)*TEMP(NH,N))* &
                        SQRT(SAL(NH,N)) + (e2 + e3*TEMP(NH,N)* &
                        TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)
                        RHO(NH,N) = P1/(P2 + p(NH))
                    END IF
                END DO
            ! Fourth, calculate the new values of sigma-t for calculation
            ! of the baroclinic pressure gradient terms.
                DO N=1,NFEN
                    SIGT(NH,N) = RHO(NH,N)-RHOWAT0
                END DO
            END DO
        END IF
    ELSE IF(Eqnstate == 3) THEN
        T_ref=10.d0
        S_ref=32.d0
    !  kmd48.33bc Equation of state from UNESCO (1980) and Haidvogel and Beckman (1999) - only used
    !        with both salinity and temperature.
        IF (ABS(IDEN) == 4) THEN
        ! Look at the values for Temperature and Salinity to see if still
        ! inside the range of the equation of state
            DO NH=1,NP
                DO N=1,NFEN
                    IF (TEMP(NH,N) < T_min) THEN
                        TEMP(NH,N) = T_min
                    END IF
                    IF (TEMP(NH,N) > T_max) THEN
                        TEMP(NH,N) = T_max
                    END IF
                    IF (SAL(NH,N) < S_min) THEN
                        SAL(NH,N) = S_min
                    END IF
                    IF (SAL(NH,N) > S_max) THEN
                        SAL(NH,N) = S_max
                    END IF
                END DO
            END DO
        ! Kendra: the initial step is to define the density for the surface
            DO NH=1,NP
                N=NFEN
                RHO(NH,N) = c00 + (c11 + (c12 +(c13 +(c14 + c15* &
                TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))* &
                TEMP(NH,N) + (c21 + (c22 + (c23 + (c24 + c25* &
                TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))* &
                SAL(NH,N) + (c31 + (c32 + c33*TEMP(NH,N))* &
                TEMP(NH,N))*SAL(NH,N)*SQRT(SAL(NH,N)) + &
                c41*SAL(NH,N)*SAL(NH,N)
            END DO
        ! Second, determine the depths of the layers
            DO NH=1,NP
                Hs = DP(NH) + IFNLFA*ETA2(NH)
                HsOAMB = Hs/(A-B)
                DO N=1, NFEN
                    Zk = HsOAMB * (Sigma(n)-B)-DP(NH)
                    Z(NH,N) = Zk
                END DO
            END DO
        ! Third, calculate the new density
            DO NH=1,NP
                DO N=NFEN-1,1,-1
                    Kappa = d00 + (d11 + (d12 + (d13 + d14*TEMP(NH,N))* &
                    TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) + (d21 + &
                    (d22 + (d23 + d24*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))* &
                    SAL(NH,N) + (d31 + (d32 + d33*TEMP(NH,N))*TEMP(NH,N))* &
                    SAL(NH,N)*SQRT(SAL(NH,N)) + (d41 + (d42 + (d43 + d44* &
                    TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N) + (d51 + &
                    (d52 + d53*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N)*SAL(NH,N) + &
                    d61*Z(NH,N)*SAL(NH,N)*SQRT(SAL(NH,N)) + (d71 + (d72 + &
                    d73*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N)*Z(NH,N) + (d81 + &
                    (d82 + d83*TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)*Z(NH,N)* &
                    Z(NH,N)
                    RHO(NH,N) = RHO(NH,NFEN)/(1.d0+(Z(NH,N)/Kappa))
                END DO
                DO N=1,NFEN
                    SIGT(NH,N) = RHO(NH,N)-RHOWAT0
                END DO
            END DO
        END IF
    END IF

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
    RETURN
!*****************************************************************************
    END SUBROUTINE CALC_SIGMAT_3D
!*****************************************************************************

!--------------------------------------------------------------------
!     S U B R O U T I N E      B I N A R Y   W R I T E   2 D
!--------------------------------------------------------------------
!     jgf49.49.01: Writes a real(sz) array with the specified length
!     to a binary hotstart file; takes the file counter as an
!     input argument and returns the new value of the file counter
!     at the end of the write.
!--------------------------------------------------------------------
    SUBROUTINE binaryWrite2D(array, length, lun, counter)
    USE SIZES, ONLY : SZ
    USE GLOBAL, ONLY : setMessageSource, unsetMessageSource, DEBUG, &
    allMessage
    IMPLICIT NONE
    REAL(SZ), intent(in), dimension(length) :: array ! data to write
    INTEGER, intent(in) :: length  ! the number of values to write
    INTEGER, intent(in) :: lun     ! fortran logical unit number to write to
    INTEGER, intent(inout) :: counter ! i/o position in the file

    INTEGER :: i      ! array index

    call setMessageSource("binaryWrite2D")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif
    DO i=1, length
        WRITE(lun,REC=counter) array(i)
        counter = counter + 1
    END DO

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
!--------------------------------------------------------------------
    END SUBROUTINE binaryWrite2D
!--------------------------------------------------------------------

