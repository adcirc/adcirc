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

c PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C                                                                             *
C  Note, input is read in READ_INPUT_3D located in module READ_INPUT          *
C        initial conditions on density, temperature and/or salinity are read  *
C        in for a cold start in subroutine COLDSTART_3D in module COLDSTART.  *
C                                                                             *
C******************************************************************************


      SUBROUTINE VSSTUP()
      USE GLOBAL_3DVS
      USE MESH, ONLY : NP      
      IMPLICIT NONE
      INTEGER :: NH, N
C
      call setMessageSource("vsstup")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C...  FOR BOTH COLD START AND HOT START, NEED TO INITIALIZE Sigma T FROM
C...  TEMPERATURE AND SALINITY DEPENDING ON VALUE OF IDEN
C...

      IF(ABS(IDen).GE.1) CALL CALC_SIGMAT_3D ()

C...
C...  COMPUTE THE DEPTH-AVERAGED Sigma T FIELD
C...

      DO NH=1,NP
         DASigT(NH)=0.d0
         IF(CBaroclinic) THEN
            DO N=1,NFEN-1
               DASigT(NH)=DASigT(NH)+(Sigma(N+1)-Sigma(N))
     &                                   *(SigT(NH,N+1)+SigT(NH,N))/2.d0
             ENDDO
            DASigT(NH)=DASigT(NH)/AMB
          ENDIF
       ENDDO

C...
C...  ADDITIONAL RUN SETUP
C...

C...  COMPUTE THE INTEGRALS Inm and LVn (INDEPENDENT OF HORIZONTAL NODE)

      CALL INMInt()

      CALL LVnInt()

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
c******************************************************************************
      END SUBROUTINE VSSTUP
c******************************************************************************



c******************************************************************************
c  SUBROUTINE VSSOL                                                           *
c                                                                             *
c  Note, the following time stepping coefficients are computed in             *
C     SUBROUTINE READ_INPUT_3DVS and passed in GLOBAL_3DVS.                   *
c                                                                             *
c  IDTAlp1      = I*DelT*Alp1        - weights coriolis term in LHS matrix    *
c  IDT1MAlp1    = I*DelT*(1.-Alp1)   - weights coriolis term in RHS forcing   *
c  DTAlp3       = DelT*Alp3          - weights vert diff term in LHS matrix   *
c  DT1MAlp3     = DelT*(1-Alp3)      - weights vert diff term in RHS forcing  *
c  DTAlp2       = DelT*Alp2          - weights bot stress term in LHS matrix  *
c  DT1MAlp2     = DelT*(1.-Alp2)     - weights bot stress term in RHS forcing *
c                                                                             *
c  Note, the following are some of the global variables are used in this      *
c  subroutine and passed in GLOBAL_3DVS                                       *
c                                                                             *
c  q(MNP,NFEN) - 3D Complex Velocity field (GAMMA) from past time step.       *
c  UU, VV - depth-averaged horizontal velocity                                *
c  DAFluxX, DAFluxY - depth-averaged horizontal flux                          *
c                                                                             *
c  NP - number of nodes in horizontal grid                                    *
c  NFEN - number of nodes in the vertical grid                                *
c  BTP - total barotropic pressure (atmos press, water level, tidal potential)*
c                 at time levels s+1/2            
! arash:
! Biharmonic_auxiliary_var(MNP,NFEN): since we need to compute the Laplacian twice, we need an auxiliary variable to store the first Laplacian for U and V.
c******************************************************************************
c                                                                             *
c   THIS  PROGRAM MAPS THE TRUE Z-COORDINATE SYSTEM INTO A DIMENSIONLESS      *
C   VERTICAL COORDINATE FROM [b,a] (BOTTOM TO TOP).  VALUES OF a=1 AND b=-1   *
C   ARE SET IN THE CODE.                                                      *
c                                                                             *
c******************************************************************************

      SUBROUTINE VSSOL(IT,TimeLoc)
C   kmd48.33bc added variable declarations for baroclinic boundary conditions
c.RJW merged 08282008 Casey 071219: Added the following variable declarations from GLOBAL.
      USE GLOBAL, ONLY : NODECODE, NOFF,
     &                   RES_BC_FLAG,
     &                   SBCTIME1, SBCTIME2, TBCTIME1,
     &                   TBCTIME2, SBCTIMEINC, TBCTIMEINC,
     &                   SBCRATIO, TBCRATIO, sponge,
     &                   BCFLAG_TEMP, BCRivRATIO,
     &                   RIVBCTIMINC, RIVBCTIME1, RIVBCTIME2,TK,
Casey 140701: Prevent depths from becoming negative in bottom roughness calculations.
     &                   H0, H1, H2
      USE GLOBAL_3DVS
      USE ADC_CONSTANTS, ONLY: G
      USE WRITE_OUTPUT, ONLY : writeOutput3D
      USE MESH, ONLY : NP, DP, NM, X, Y, NODELE, AREAS, LBArray_Pointer,  
     &                  NEITABELE, NEITAB, NNEIGH, SFAC
      USE BOUNDARIES, ONLY : NOPE, NETA, NBD, LBCODEI, 
     &   BndBCRiver, totalbcrivernodes, SIII, CSII
      USE NodalAttributes, ONLY: EVM,
     &     ManningsN, LoadManningsN,
Corbitt 120328: Local Advection
     &     AdvectionState, LoadAdvectionState,
     &     Z0b_var, LoadZ0B_var
#ifdef CMPI
      USE MESSENGER
#endif      
      IMPLICIT NONE

c.RJW merged 08282008 Casey 071219: Added the following local variable declarations.
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
Cjgfdebug      INTEGER :: NN             !output loop counter

      REAL(8) :: KSlip         !equavalent linear slip coeff
      REAL(8) :: Z0B1          !Bottom roughness length (from mannings, const or from Z0Bvar)
      REAL(8) :: WSXsNH,WSYsNH !Wind stress components at time level s at node NH

      REAL(8) :: WSigma(MNFEN) !"sigma" vertical velocity
      REAL(8) :: Wf            !weighting coefficient in adjoint correction to w
      REAL(8) :: WfOHH         !Wf/(H time level s+1)^2
      REAL(8) :: WZSurfBC      !surface boundary condition value of w
      REAL(8) :: WZSurf        !computed value of w at surface
      REAL(8) :: WZCorrection  !adjoint correction compute for w

C   kmd48.33bc add in spherical factors
      REAL(8) :: SFacAvg
Cjgfdebug      REAL(8) :: Zk            !z depth of any node k in the vertical
Cjgfdebug      REAL(8) :: DelSig        ! sigma(k+1)-sigma(k)
      REAL(8) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
      REAL(8) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
      REAL(8) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
      REAL(8) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
Cjgfdebug      REAL(8) :: SigmaNN       !Sigma value of a neighbor node

      REAL(8) :: VelNorm,VelTan !-QNormsp1(NH)/Hsp1 at flux boundary node
      REAL(8) :: CLBP,SLBP     !local values of CSII, SIII at boundary node LBP
      REAL(8) :: Auv1km1,Auv2km1 !initial real,imaginary parts of Mkm1 at flux boundary node
      REAL(8) :: Auv1km1star,Auv2km1star !rotated real,imaginary parts of Mkm1 at flux boundary node
Cjgfdebug      REAL(8) :: Auv1k1,Auv2k1 !initial real,imaginary parts of Mk at flux boundary node
Cjgfdebug      REAL(8) :: Auv1k1star,Auv2k1star !rotated real,imaginary parts of Mk at flux boundary node
      REAL(8) :: Auv1kp1,Auv2kp1 !initial real,imaginary parts of Mkp1 at flux boundary node
      REAL(8) :: Auv1kp1star,Auv2kp1star !rotated real,imaginary parts of Mkp1 at flux boundary node

      REAL(8) :: EtaN1,EtaN2,EtaN3,EtaNFirst !nodal values of IFNLFA(Eta1+Eta2)/2
      REAL(8) :: hN1,hN2,hN3,hNFirst !nodal values of DP
      REAL(8) :: DUDX(MNFEN),DVDY(MNFEN) !horizontal derivatives of velocity used to compute w
      REAL(8) :: Un,Vn         !real,imaginary components of qn
      REAL(8) :: DelU,DelV     !real, imaginary parts of q(k)-q(k-1)

      REAL(8) :: BTPN1,BTPN2,BTPN3,BTPNFirst !nodal values of BTP
      REAL(8) :: BTPDX2A,BTPDY2A !(Horiz. grads of BTP)*2*Element Area

Cjgfdebug      REAL(8) :: BCPN1,BCPN2,BCPN3,BCPNFirst !nodal values of BCP
Cjgfdebug      REAL(8) :: BCPDX2A,BCPDY2A !(Horiz. grads of BCP)*2*Element Area
Cjgfdebug      REAL(8) :: SigTAvg       !avg SigT between 2 vertical nodes
Cjgfdebug      REAL(8) :: HGORhoOAMB    !depth*gravity/(reference density)/(a-b)

Cjgfdebug      REAL(8) :: DBCPDX2A
Cjgfdebug      REAL(8) :: DBCPDY2A
      REAL(8) :: RCL
      REAL(8) :: Auv1k
      REAL(8) :: Auv2k
      REAL(8) :: Auv1kstar
      REAL(8) :: Auv2kstar

      REAL(8) :: KVnm(MNFEN,3)  !integral used in vertical stress term
      REAL(8) :: TimeLoc           !model time at time level s+1
      REAL(8) :: DEtaDT         !time derivative of water surface elev
      REAL(8) :: DEtaDX,DEtaDY  !horizontal derivatives of water surface elev at time level s
      REAL(8) :: DEtaDX2A,DEtaDY2A !(DEtaDX,DEtaDY)*2*Element Area
      REAL(8) :: DhDX,DhDY      !horizontal derivatives of DP
      REAL(8) :: DhDX2A,DhDY2A  !(DhDX,DhDY)*2*Element Area
      REAL(8) :: TotalArea2     !2*Area of all elements around a node
Cjgfdebug      REAL(8) :: TotalBCPGArea2 !2*Area of all elements around a node used to compute the BCPG
      REAL(8) :: a1,a2,a3,b1,b2,b3 !elemental coefficients used in horizontal FE method

      REAL(8) :: Hs             !Total water depth at time level s
Cjgfdebug      REAL(8) :: HsN2           !Total water depth at time level s at local node N2
      REAL(8) :: HsOAMB         !Hs/(a-b)
      REAL(8) :: HsHsOAMBAMB    !(Hs/(a-b))^2
      REAL(8) :: Hsp1           !Total water depth at time level s+1
      REAL(8) :: Hsp1OAMB       !Hsp1/(a-b)
      REAL(8) :: Hsp1Hsp1OAMBAMB !(Hsp1/(a-b))^2


Casey 140701: Prevent depths from becoming negative in bottom roughness calculations.
! arash 12/14/2015
      Real(8) :: Biharmonic_viscosity

      COMPLEX(8) :: Fr(MNFEN)      !right side forcing vector
      COMPLEX(8) :: Frstar         !rotated right side forcing vector at flux boundary node
      COMPLEX(8) :: Mkm1(MNFEN)    !1st column (k-1) in left side compact storage matrix
      COMPLEX(8) :: Mk(MNFEN)      !2nd column (k) in left side compact storage matrix
      COMPLEX(8) :: Mkp1(MNFEN)    !3rd column (k+1) in left side cpmpact storage matrix
      COMPLEX(8) :: LAdvec(MNFEN)  !lateral advection term in momentum eqn
      COMPLEX(8) :: LStress(MNFEN) !lateral stress term in momentum eqn
! arash 12/14/2015
      COMPLEX(8) :: Biharmonic_LStress(MNFEN) !lateral stress term in momentum eqn - biharmonic viscosity

      COMPLEX(8) :: VAdvec(MNFEN)  !vertical advection term in momentum eqn
      COMPLEX(8) :: VStress(MNFEN) !vertical stress term in momentum eqn
Cjgfdebug      COMPLEX(8) :: BCPG(MNFEN)    !baroclinic pressure gradient
      COMPLEX(8) :: BTPG           !total barotropic pressure gradient (incl TP & water level)
      COMPLEX(8) :: CCR,CCL        !coeffs used on right,left side of momentum eqn
      COMPLEX(8) :: VIVel          !vertically integrated velocity
      COMPLEX(8) :: DUDS           !complex vertical velocity gradient between bottom two nodes
      COMPLEX(8) :: qn,qN1,qN2,qN3,qNFirst !nodal values of q
Cjgfdebug      COMPLEX(8) :: UnDqDX,VnDqDY  !derivatives used in lateral advection
      COMPLEX(8) :: UnDqDX2A,VnDqDY2A !UnDqDX,UnDqDY)*2*Element Area
      COMPLEX(8) :: DqDXDPhiDX2A,DqDYDPhiDY2A !derivatives used in lateral stress calc.
! arash 12/14/2015
      COMPLEX(8) :: Biharmonic_DqDXDPhiDX2A, Biharmonic_DqDYDPhiDY2A !derivatives of Biharmonic auxiliary variables.
      COMPLEX(8) :: Biharmonic_qn, Biharmonic_qN1, Biharmonic_qN2, Biharmonic_qN3, Biharmonic_qNFirst
! arash 12/21/2015: modified Leith viscosity
      COMPLEX(8) :: d_dx_q_element, d_dy_q_element
      COMPLEX(8) :: d_dx_qNFirst, d_dy_qNFirst
      COMPLEX(8) :: d2_dx2_q_element, d2_dy2_q_element, d2_dxy_q_element
      COMPLEX(8) :: d_dx_qN1, d_dx_qN2, d_dx_qN3, d_dy_qN1, d_dy_qN2, d_dy_qN3
      Real(8)     :: modified_Leith_delta_x, modified_Leith_Lambda_vorticity, modified_Leith_Lambda_divergence
      Real(8)     :: modified_Leith_grad_vorticity, modified_Leith_grad_divergence
! arash 01/04/2016: compute diffusion at z-levels
      REAL(8)    :: Zk            !z depth of any node k in the vertical
      REAL(8)    :: SigmaNN       !Sigma value of a neighbor node
      REAL(8)     :: HsN2          !Total water depth at time level s at local node N2
      COMPLEX(8) :: qN_z, qN1_z, qN2_z, qN3_z, qNFirst_z !nodal values of q at z-levels
      COMPLEX(8) :: Biharmonic_qn_z, Biharmonic_qN1_z, Biharmonic_qN2_z, Biharmonic_qN3_z, Biharmonic_qNFirst_z
      COMPLEX(8) :: d_dx_qN1_z, d_dx_qN2_z, d_dx_qN3_z, d_dy_qN1_z, d_dy_qN2_z, d_dy_qN3_z
      COMPLEX(8) :: d_dx_qNFirst_z, d_dy_qNFirst_z
      Complex(8) :: inarray_vertical(NFEN) 

      COMPLEX(8) :: DqDSigmakm1,DqDSigmakp1 !vertical deriv. of q from k-1,k and k,k+1
      COMPLEX(8) :: DqDX2A(MNFEN)  !horizontal derivatives of complex
      COMPLEX(8) :: DqDY2A(MNFEN)  !velocity used in w calc
! arash Feb/18/2016
      INTEGER     :: Neighbor_Node
      REAL(8)    :: length, length_max
      REAL(8)    :: velocity_magnitude, element_peclet_number

C   kmd48.33bc added the following variables for boundary conditions
      INTEGER :: J, NumofBCNode
      INTEGER :: NOD
      CHARACTER(80) :: CDUM80

      call setMessageSource("vssol")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C***************************************************************************
C     Set up several variables that are needed for 3D run
C
C     Zero out surface stress and stress forcing at top for case of no wind
C
      WSXsNH=0.d0
      WSYsNH=0.d0
C
C     If a baroclinic run, compute the 3D baroclinic pressure field
C     The buoyancy field is defined as
C     BCP(z)    =(gravity/rho ref)*          integral (SigT) from surface down to z
C     BCP(sigma)=(gravity/rho ref)*(H/(a-b))*integral (SigT) from a down to sigma
C     where
C     SigT = Sigma T = Rho - 1000 = density - 1000
C     SigT0 = Sigma t value of reference density (typically = 0)
C     Sigma = dimensionless vertical coordinate
C
C!      IF(CBaroclinic) THEN
C!        DO NH=1,NP             !loop over horizontal nodes
C!            Hs=DP(NH)+IFNLFA*Eta1(NH) !total depth at previous (s) timestep
C!            HGORhoOAMB=GORhoOAMB*Hs !(gravity/rho ref)*(H/(a-b))
C!            BCP(NH,NFEN)=0.d0
C!            DO k=NFEN-1,1,-1    !loop over vertical nodes, starting at top and working down
C!               SigTAvg=(SigT(NH,k+1)+SigT(NH,k))/2.d0
C!               DelSig=Sigma(k+1)-Sigma(k)
C!               BCP(NH,k)=BCP(NH,k+1)+HGORhoOAMB*(SigTAvg-SigT0)*DelSig
C!            ENDDO
C!         ENDDO
#ifdef CMPI
C     Update BCP on ghost nodes
C      CALL UPDATER3D(BCP)  !!!!Don't know if this is needed at this time
#endif
C!      ENDIF

! -------------------------------------------------------------------------------------------------
! arash 12/14/2015: biharmonic diffusion for the momentum equation
! before computing 3D horizontal velocities, we need to compute the inner laplacian of the biharmonic operator.
! Del Del q = Del Biharmonic_auxiliary_var
      Biharmonic_auxiliary_var (:,:) = ( 0.0d0 , 0.0d0 )
      d_dx_q_nodes (:,:) = ( 0.0d0 , 0.0d0 )   ! spatial derivarievs to compute modified Leith viscosity
      d_dy_q_nodes (:,:) = ( 0.0d0 , 0.0d0 )

      ! save modified Leith viscosity parameter and use it later in the transport equations
      Biharmonic_viscosity_modified_Leith (:,:) = 0.0d0

      DO NH = 1, NP                      !loop over horizontal nodes

C........DMW202401 Use saved H1 !Total depth at current timestep - 'center' node
         HsOAMB = H1(NH) / AMB
         ! begin the vertical loop
         DO k = 1, NFEN
            DqDXDPhiDX2A = 0.d0    ! variables for computing the biharmonic operator
            DqDYDPhiDY2A = 0.d0
            TotalArea2   = 0.d0
            d_dx_q_element = ( 0.0d0 , 0.0d0 ) ! to compute the viscosity constant
            d_dy_q_element = ( 0.0d0 , 0.0d0 )
            Zk = HsOAMB * ( Sigma(k) - B ) - DP(NH)

            N1  = NH               !node 1 is always the central node
            qN1 = q(N1,k)
            qN1_z = qN1
            N2  = NEITAB(NH,2)     !operate on 1st neighbor
            qN2 = q(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
C...........DMW202401 Use saved H1
            SigmaNN = B + AMB * ( Zk + DP(N2) ) / H1(N2) !equivalent sigma value at neighbor
            inarray_vertical (:) = q(N2,:)
            Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, qN2_z, k )

            NNFirst = N2           !save these values until end
            qNFirst = qN2
! af 01/04/2016
            qNFirst_z = qN2_z

            DO N = 3, NNEIGH(NH)   !operate on rest of neighbors
               N3  = N2            !shift previously computed values
               qN3 = qN2
! af 01/04/2016
               qN3_z = qN2_z

               N2  = NEITAB(NH,N)  !select new neighbor to work on
               qN2 = q(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
C..............DMW202401 Use saved H1
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / H1(N2) !equivalent sigma value at neighbor
               inarray_vertical (:) = q(N2,:)
               Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, qN2_z, k )

               NEle= NeiTabEle(NH,N-2) !element # defined by nodes N1,N2,N3
               IF (NEle.eq.0) THEN
                  CYCLE
               ENDIF
               IF((NEle.NE.0)) THEN    !if element is active, compute velocity grads
                  NCELE = NODECODE(NH) * NODECODE(N2) *
     &                 NODECODE(N3) * NOFF(NELE)
                  IF (NCELE.NE.0) THEN !if element exists, compute terms
                     TotalArea2 = TotalArea2 + Areas(NEle) !accumulate 2X total areas to complete calc.
                     SFacAvg = (SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                     !qn = (qN1+qN2+qN3)/3.d0
                     !Un = REAL(qn)
                     !Vn = AIMAG(qn)
                     a1 = X(N3)-X(N2)
                     a2 = X(N1)-X(N3)
                     a3 = X(N2)-X(N1)
                     b1 = (Y(N2)-Y(N3))*SFacAvg
                     b2 = (Y(N3)-Y(N1))*SFacAvg
                     b3 = (Y(N1)-Y(N2))*SFacAvg
! iso-sigma:
!                     DqDXDPhiDX2A = DqDXDPhiDX2A
!     &                 + (qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
!                     DqDYDPhiDY2A=DqDYDPhiDY2A
!     &                 +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
                     ! variables for modified Leith viscosity parameter
!                     d_dx_q_element = d_dx_q_element + qN1*b1+qN2*b2+qN3*b3
!                     d_dy_q_element = d_dy_q_element + qN1*a1+qN2*a2+qN3*a3
! z-level:
                     DqDXDPhiDX2A = DqDXDPhiDX2A
     &                 + (qN1_z*b1+qN2_z*b2+qN3_z*b3)*b1/Areas(NEle)
                     DqDYDPhiDY2A=DqDYDPhiDY2A
     &                 +(qN1_z*a1+qN2_z*a2+qN3_z*a3)*a1/Areas(NEle)
                     ! variables for modified Leith viscosity parameter
                     d_dx_q_element = d_dx_q_element + qN1_z*b1+qN2_z*b2+qN3_z*b3
                     d_dy_q_element = d_dy_q_element + qN1_z*a1+qN2_Z*a2+qN3_z*a3
                  ENDIF
               ENDIF
            ENDDO
            N3  = N2               !wrap back to beginning to get final contribution
            qN3 = qN2
! af 01/04/2016
            qN3_z = qN2_z

            N2  = NNFIRST
            qN2 = qNFirst
! af 01/04/2016
            qN2_z = qNFirst_z

            NEle=NeiTabEle(NH,NNeigh(NH)-1)
            IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
               NCELE = NODECODE(NH) * NODECODE(N2) *
     &              NODECODE(N3) * NOFF(NELE)
               IF ((NCELE.NE.0)) THEN  !if element is active, compute velocity grads
                  TotalArea2 = TotalArea2 + Areas(NEle) !accumulate 2X total areas to complete calc.
                  SFacAvg = (SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  !qn = (qN1+qN2+qN3)/3.d0
                  !Un = real(qn)
                  !Vn = aimag(qn)
                  a1 = X(N3)-X(N2)
                  a2 = X(N1)-X(N3)
                  a3 = X(N2)-X(N1)
                  b1 = (Y(N2)-Y(N3))*SFacAvg
                  b2 = (Y(N3)-Y(N1))*SFacAvg
                  b3 = (Y(N1)-Y(N2))*SFacAvg
! iso-sigma:
!                  DqDXDPhiDX2A = DqDXDPhiDX2A
!     &                 +(qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
!                  DqDYDPhiDY2A = DqDYDPhiDY2A
!     &                 +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
                  ! variables for modified Leith viscosity parameter - velocity derivatives at nodes
!                  d_dx_q_element = d_dx_q_element + qN1*b1+qN2*b2+qN3*b3
!                  d_dy_q_element = d_dy_q_element + qN1*a1+qN2*a2+qN3*a3
! z-level:
                     DqDXDPhiDX2A = DqDXDPhiDX2A
     &                 + (qN1_z*b1+qN2_z*b2+qN3_z*b3)*b1/Areas(NEle)
                     DqDYDPhiDY2A=DqDYDPhiDY2A
     &                 +(qN1_z*a1+qN2_z*a2+qN3_z*a3)*a1/Areas(NEle)
                     ! variables for modified Leith viscosity parameter
                     d_dx_q_element = d_dx_q_element + qN1_z*b1+qN2_z*b2+qN3_z*b3
                     d_dy_q_element = d_dy_q_element + qN1_z*a1+qN2_Z*a2+qN3_z*a3
               ENDIF
            ENDIF

            IF(TotalArea2.EQ.0.) THEN
               ! do nothing
            ELSE
               Biharmonic_auxiliary_var(NH,k) = -3.d0 * (DqDXDPhiDX2A+DqDYDPhiDY2A) / TotalArea2
               ! variables for modified Leith viscosity parameter
               d_dx_q_nodes (NH,k) = d_dx_q_element / TotalArea2
               d_dy_q_nodes (NH,k) = d_dy_q_element / TotalArea2
            ENDIF

         END DO ! end vertical loop
      End Do ! end horizontal loop

! end computing Biharmonic_auxiliary_var(:,:)
! -------------------------------------------------------------------------------------------------

C***************************************************************************
C     Compute 3D horizontal velocities
C
C     Loop over each horizontal node to compute the horizontal velocity
C
      DO NH=1,NP                !loop over horizontal nodes

Corbitt 120328: Applies Advection Locally
         IF (LoadAdvectionState) Then
            IF (DP(NH).GE.AdvectionState(NH)) THEN
               IFNLCT = IFNLCTE
            ELSE
               IFNLCT = 0
            ENDIF
         ENDIF

c     Set up some values at the node being worked on
C........DMW202401 Use saved H1
         Hs  = H1(NH) !Total depth at previous (s) timestep
Casey 140701: Check the total water depth.
         IF(Hs.LT.H0) Hs = H0
         HsOAMB=Hs/AMB
         HsHsOAMBAMB=HsOAMB*HsOAMB
C........DMW202401 Use saved H2
         Hsp1= H2(NH) !Total depth at present (s+1) timestep
Casey 140701: Check the total water depth.
         IF(Hsp1.LT.H0) Hsp1 = H0
         Hsp1OAMB=Hsp1/AMB
         Hsp1Hsp1OAMBAMB=Hsp1OAMB*Hsp1OAMB

         IF(NWS.NE.0) THEN      !wind stress
            WSXsNH=WSX1(NH)
            WSYsNH=WSY1(NH)
         ENDIF

c.RJW merged 08282008 Casey 071219: Solve for TEMPNCELE, which is the sum of the NCELE values
C     for the elements attached to the current horizontal node.  Note that,
C     for some reason, the k neighbor elements are not contained in the first
C     k registers of NEITABELE.  (For instance, if a node is connected to four
C     elements, those element numbers might be contained in registers
C     1, 2, 4, and 5.)  So we have to be cute about pulling element numbers
C     from NEITABELE.
C
         TEMPNCELE = 0
         TEMPSTOP = 0
Casey 140701: Debug.
C        DO K=1,NODELE(NH)
         DO K=1,SIZE(NEITABELE,2)
            IF(NEITABELE(NH,K).EQ.0)THEN
C              TEMPSTOP = TEMPSTOP + 1
               CYCLE
            ENDIF
            TEMPNCELE = TEMPNCELE +
C    &      NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),1))*
C    &      NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),2))*
C    &      NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),3))*
C    &      NOFF(NEITABELE(NH,K+TEMPSTOP))
     &      NODECODE(NM(NEITABELE(NH,K),1))*
     &      NODECODE(NM(NEITABELE(NH,K),2))*
     &      NODECODE(NM(NEITABELE(NH,K),3))*
     &      NOFF(NEITABELE(NH,K))
         ENDDO

Casey 071219: If TEMPNCELE is zero, then the current horizontal node is
C     not attached to any active elements.  Thus, the horizontal velocities
C     must be zero, and we can set up the matrix immediately and skip to the
C     tri-diagonal solver.
C

! arash: this means that we assign zero velocities to a dry node; right?

         IF(TEMPNCELE.EQ.0)THEN
            DO K=1,NFEN
               MKM1(K) = (0.D0,0.D0)
               MK(K)   = (1.D0,-1.D0)
               MKP1(K) = (0.D0,0.D0)
               FR(K)   = (0.D0,0.D0)
            ENDDO
Casey 140701: Need to set some variables for use later.
            KSlip = 0.D0
            !write(6,*) 'temncele.eq.0 skip to 999 nh=',nh !jgfdebug
            GOTO 999
         ENDIF
c.....
c     If specified normal flow boundary node with no tangential slip is
c     lateral boundary condition, set up the matrix immediately and skip
c     to the solution

         LBP=LBArray_Pointer(NH)
         IF(LBP.GT.0) THEN
                                !essential normal flow
            IF((LBCodeI(LBP).GE.10).AND.(LBCodeI(LBP).LE.19)) THEN
               VelNorm=-QNormsp1(LBP)/Hsp1 !with essential no tan.
               VelTan=0.D0
               SLBP=SIII(LBP)
               CLBP=CSII(LBP)
               DO k=1,NFEN
                  Mkm1(k)=(0.D0,0.D0)
                  Mk(k)=SLBP+iy*CLBP
                  Mkp1(k)=(0.D0,0.D0)
                  Fr(k)=VelTan+iy*VelNorm
C                  EVTot(k)=0.D0 !f77diff??
               ENDDO
               !write(6,*) 'lbp=',lbp,' lbcodei(lbp)=',lbcodei(lbp),' skip to 999 nh=',nh !jgfdebug               
               GOTO 999
            ENDIF
         ENDIF

C Determine the bottom roughness length either from fort.15, from Manning's n
C or as read in from nodal attributes
         IF(LoadZ0B_var) THEN
            Z0B1 = Z0B_var(NH)
         ELSEIF (LoadManningsN) THEN
Casey 140701: Prevent depths from becoming negative.  Cannot raise a negative number
C             to the 1/6 power.
C...........DMW202401 Use saved H2
            Hs = H2(NH)
            IF(Hs.LT.H0) Hs=H0
!            Z0B1 = ( DP(NH)+IFNLFA*ETA2(NH) )* exp(-(1.0D0+
!     &             ( (0.41D0*( DP(NH)+IFNLFA*ETA2(NH))**(1.0D0/6.0D0) )/
!     &                              (ManningsN(NH)*sqrt(g)) ) ))
! arash Sep 7 2015:
            Z0B1 = ( Hs )* exp(-(1.0D0+
     &             ( (0.41D0*( Hs )**(1.0D0/6.0D0) )/
     &                              (ManningsN(NH)*sqrt(g)) ) ))
         ELSE
            Z0B1 = Z0B
         ENDIF

c     Compute the vertical eddy viscosity
C... RJW changed Hs to Hsp1
C Hs is the previous time step (if node is recently wet then this is no good
C Hsp1 is current time step, if node is newly wet then we now have agreement
C and Hsp1 has a value that can be used
         CALL EDDYVIS(Hsp1,UU(NH),VV(NH),
     &        WSXsNH,WSYsNH,BSX(NH),BSY(NH),NH,IT,Z0B1)

c     Compute the integral KVnm

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

c     Compute time derivative of water surface position

         DEtaDT=(Eta2(NH)-Eta1(NH))/DelT

c     Start computing horizontal derivatives of water level, bathymetric
c     depth and total barotropic pressure (atmos pres, water level,
c     tidal potential) Note: TotalArea2 = 2X total elemental area
c     surrounding a node

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
            NEle=NeiTabEle(NH,N-2) !element # defined by nodes NH,NN2,NN1
c.RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
            !jgf48.50 If NEle is zero, the NOFF check would be out of
            ! bounds, making it look like NOFF is always zero (depending
            ! on how the compiler handles array out-of-bounds)
            IF (NEle.eq.0) THEN
               CYCLE
            ENDIF
c.RJW 06092009  modified to elimiate zero index NOFF
            IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
               NCELE = NODECODE(NH)*NODECODE(N2)*
     &              NODECODE(N3)*NOFF(NELE)
               IF((NCELE.NE.0)) THEN  !if element is active, compute velocity grads

                  TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C               kmd48.33bc added in spherical factors for 3D
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
c. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
         !jgf48.50 Avoid array out of bounds
          IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
            NCELE = NODECODE(NH)*NODECODE(N2)*
     &              NODECODE(N3)*NOFF(NELE)
            IF((NCELE.NE.0)) THEN  !if element is active, compute velocity grads
               TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C    kmd48.33bc added in spherical factors for 3D
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

         IF(TotalArea2.NE.0.) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
            BTPG=(BTPDX2A+iy*BTPDY2A)/TotalArea2
         ENDIF

c     Finished computing horizontal derivatives of water level,
c     bathymetric depth and total barotropic pressure (atmos pres, water
c     level, tidal potential)

c     Compute the "sigma" vertical velocity from the "z" vertical velocity

         DO k=1,NFEN
            SigmaMAOAMB=(Sigma(k)-A)/AMB
            SigmaMBOAMB=(Sigma(k)-B)/AMB
            WSigma(k) = WZ(NH,k) - SigmaMBOAMB*DEtaDT
     &           - REAL(q(NH,k))*(SigmaMBOAMB*DEtaDX+SigmaMAOAMB*DhDX)
     &           - AIMAG(q(NH,k))*(SigmaMBOAMB*DEtaDY+SigmaMAOAMB*DhDY)
         ENDDO


c     Start computing advection and stress terms in the momentum
c     equation at each level in the vertical

!--------------------------------------------------------------------------------------------------
! arash: Feb/18/2016
!     beginning of vertical integration
         length_max = 0.0d0
         Do J = 2, NNeigh ( NH )
             Neighbor_Node = NeiTab ( NH , J )
             length = SQRT ( ( X(Neighbor_Node) - X(NH) )**2 + ( Y(Neighbor_Node) - Y(NH) )**2 )
             If ( length > length_max ) length_max = length
         End Do
! -------------------------------------------------------------------------------------------------

         DO k=1,NFEN

c     Compute the vertical advection and vertical stress terms

            IF(k.EQ.1) THEN
               DqDSigmakp1=(q(NH,k+1)-q(NH,k))/(Sigma(k+1)-Sigma(k))
               VAdvec(k)=DqDsigmakp1*
     &              (2.d0*WSigma(k)+WSigma(k+1))*Inm(k,3)/HsOAMB
               VStress(k)=(q(NH,k)*KVnm(k,2)+q(NH,k+1)*KVnm(k,3))
     &              /HsHsOAMBAMB
            ENDIF
            IF((k.GT.1).AND.(k.LT.NFEN)) THEN
               DqDSigmakm1=DqDSigmakp1
               DqDSigmakp1=(q(NH,k+1)-q(NH,k))/(Sigma(k+1)-Sigma(k))
               VAdvec(k)=
     &              ( DqDSigmakm1*(WSigma(k-1)+2.d0*WSigma(k))*Inm(k,1)
     &              +DqDSigmakp1*(2.d0*WSigma(k)+WSigma(k+1))*Inm(k,3) )
     &              /HsOAMB
               VStress(k)=(q(NH,k-1)*KVnm(k,1)+q(NH,k)*KVnm(k,2)
     &              +q(NH,k+1)*KVnm(k,3))/HsHsOAMBAMB
            ENDIF
            IF(k.EQ.NFEN) THEN
               DqDSigmakm1=DqDSigmakp1
               VAdvec(k)=
     &              DqDSigmakm1*(WSigma(k-1)+2.d0*WSigma(k))*Inm(k,1)
     &              /HsOAMB
               VStress(k)=(q(NH,k-1)*KVnm(k,1)+q(NH,k)*KVnm(k,2))
     &              /HsHsOAMBAMB
            ENDIF
            VAdvec(k)=sponge(NH)*IFNLCT*VAdvec(k)

!--------------------------------------------------------------------------------------------------
c     Compute lateral advection and lateral stress terms
! -------------------------------------------------------------------------------------------------

            UnDqDX2A=0.d0
            VnDqDY2A=0.d0
            DqDXDPhiDX2A=0.d0
            DqDYDPhiDY2A=0.d0
! arash 12/14/2015: Biharmonic auxiliary variable derivatives
            Biharmonic_DqDXDPhiDX2A =  ( 0.d0 , 0.0d0 )
            Biharmonic_DqDYDPhiDY2A =  ( 0.d0 , 0.0d0 )
! arash 12/21/2015: for modified Leith viscosity
            d2_dx2_q_element = ( 0.d0 , 0.0d0 )
            d2_dy2_q_element = ( 0.d0 , 0.0d0 )
            d2_dxy_q_element = ( 0.d0 , 0.0d0 )

            Zk = HsOAMB * ( Sigma(k) - B ) - DP(NH)

            N1=NH               !node 1 is always the central node
            qN1=q(N1,k)
! af 01/04/2016
            qN1_z = qN1
! arash 12/14/2015
            Biharmonic_qN1 = Biharmonic_auxiliary_var(N1,k)
! af 01/04/2016
            Biharmonic_qN1_z = Biharmonic_qN1
! arash 12/21/2015
            d_dx_qN1 = d_dx_q_nodes (N1,k)
            d_dy_qN1 = d_dy_q_nodes (N1,k)
! af 01/04/2016
            d_dx_qN1_z = d_dx_qN1
            d_dy_qN1_z = d_dy_qN1

            N2=NEITAB(NH,2)     !operate on 1st neighbor
            qN2=q(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
C...........DMW202401 Use saved H1
            SigmaNN = B + AMB * ( Zk + DP(N2) ) / H1(N2) !equivalent sigma value at neighbor
            inarray_vertical (:) = q(N2,:)
            Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, qN2_z, k )

! arash 12/14/2015
            Biharmonic_qN2 = Biharmonic_auxiliary_var(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
            inarray_vertical (:) = Biharmonic_auxiliary_var(N2,:)
            Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, Biharmonic_qN2_z, k )

! arash 12/21/2015
            d_dx_qN2 = d_dx_q_nodes (N2,k)
            d_dy_qN2 = d_dy_q_nodes (N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
            inarray_vertical (:) = d_dx_q_nodes(N2,:)
            Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, d_dx_qN2_z, k )
            inarray_vertical (:) = d_dy_q_nodes(N2,:)
            Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, d_dy_qN2_z, k )

            NNFirst=N2          !save these values until end
            qNFirst=qN2
! arash 12/14/2015
            Biharmonic_qNFirst = Biharmonic_qN2
! arash 12/21/2015
            d_dx_qNFirst = d_dx_qN2
            d_dy_qNFirst = d_dy_qN2
! af 01/04/2016
            qNFirst_z = qN2_z
            Biharmonic_qNFirst_z = Biharmonic_qN2_z
            d_dx_qNFirst_z = d_dx_qN2_z
            d_dy_qNFirst_z = d_dy_qN2_z

            DO N = 3, NNEIGH(NH)   !operate on rest of neighbors
               N3  = N2            !shift previously computed values
               qN3 = qN2
! arash 12/22/2015
               Biharmonic_qN3 = Biharmonic_qN2
               d_dx_qN3 = d_dx_qN2
               d_dy_qN3 = d_dy_qN2
! af 01/04/2016
               qN3_z = qN2_z
               Biharmonic_qN3_z = Biharmonic_qN2_z
               d_dx_qN3_z = d_dx_qN2_z
               d_dy_qN3_z = d_dy_qN2_z

               N2  = NEITAB(NH,N)  !select new neighbor to work on
               qN2 = q(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
C..............DMW202401 Use saved H1
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / H1(N2) !equivalent sigma value at neighbor
               inarray_vertical (:) = q(N2,:)
               Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, qN2_z, k )

! arash 12/14/2015
               Biharmonic_qN2 = Biharmonic_auxiliary_var(N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
               inarray_vertical (:) = Biharmonic_auxiliary_var(N2,:)
               Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, Biharmonic_qN2_z, k )
! arash 12/21/2015
               d_dx_qN2 = d_dx_q_nodes (N2,k)
               d_dy_qN2 = d_dy_q_nodes (N2,k)
! af 01/04/2016: interpolate it at the same z level as N1
               inarray_vertical (:) = d_dx_q_nodes(N2,:)
               Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, d_dx_qN2_z, k )
               inarray_vertical (:) = d_dy_q_nodes(N2,:)
               Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, d_dy_qN2_z, k )

               NEle=NeiTabEle(NH,N-2) !element # defined by nodes N1,N2,N3
c. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
               !jgf48.50 Avoid array out of bounds
               IF (NEle.eq.0) THEN
                  CYCLE
               ENDIF
               IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
                  NCELE = NODECODE(NH)*NODECODE(N2)*
     &                 NODECODE(N3)*NOFF(NELE)
                  IF (NCELE.NE.0) THEN !if element exists, compute terms
C    kmd48.33bc added in spherical factors for 3D
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
! iso-sigma:
!                     DqDXDPhiDX2A=DqDXDPhiDX2A
!     &                 +(qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
!                     DqDYDPhiDY2A=DqDYDPhiDY2A
!     &                 +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
! arash 12/14/2015
!                     Biharmonic_DqDXDPhiDX2A = Biharmonic_DqDXDPhiDX2A
!     &                 +(Biharmonic_qN1*b1+Biharmonic_qN2*b2+Biharmonic_qN3*b3)*b1/Areas(NEle)
!                     Biharmonic_DqDYDPhiDY2A = Biharmonic_DqDYDPhiDY2A
!     &                 +(Biharmonic_qN1*a1+Biharmonic_qN2*a2+Biharmonic_qN3*a3)*a1/Areas(NEle)
! arash 12/21/2015
!                     d2_dx2_q_element = d2_dx2_q_element + d_dx_qN1*b1+d_dx_qN2*b2+d_dx_qN3*b3
!                     d2_dy2_q_element = d2_dy2_q_element + d_dy_qN1*a1+d_dy_qN2*a2+d_dy_qN3*a3
!                     d2_dxy_q_element = d2_dxy_q_element + d_dy_qN1*b1+d_dy_qN2*b2+d_dy_qN3*b3
! z-level
                     DqDXDPhiDX2A=DqDXDPhiDX2A
     &                 +(qN1_z*b1+qN2_z*b2+qN3_z*b3)*b1/Areas(NEle)
                     DqDYDPhiDY2A=DqDYDPhiDY2A
     &                 +(qN1_z*a1+qN2_z*a2+qN3_z*a3)*a1/Areas(NEle)

                     Biharmonic_DqDXDPhiDX2A = Biharmonic_DqDXDPhiDX2A
     &                 +(Biharmonic_qN1_z*b1+Biharmonic_qN2_z*b2+Biharmonic_qN3_z*b3)*b1/Areas(NEle)
                     Biharmonic_DqDYDPhiDY2A = Biharmonic_DqDYDPhiDY2A
     &                 +(Biharmonic_qN1_z*a1+Biharmonic_qN2_z*a2+Biharmonic_qN3_z*a3)*a1/Areas(NEle)

                     d2_dx2_q_element = d2_dx2_q_element + d_dx_qN1_z*b1+d_dx_qN2_z*b2+d_dx_qN3_z*b3
                     d2_dy2_q_element = d2_dy2_q_element + d_dy_qN1_z*a1+d_dy_qN2_z*a2+d_dy_qN3_z*a3
                     d2_dxy_q_element = d2_dxy_q_element + d_dy_qN1_z*b1+d_dy_qN2_z*b2+d_dy_qN3_z*b3
                  ENDIF
               ENDIF
            ENDDO

            N3 = N2               !wrap back to beginning to get final contribution
            qN3 = qN2
! af 01/04/2016
            qN3_z = qN2_z
! arash 12/14/2015
            Biharmonic_qN3 = Biharmonic_qN2
! af 01/04/2016
            Biharmonic_qN3_z = Biharmonic_qN2_z
! arash 12/21/2015
            d_dx_qN3 = d_dx_qN2
            d_dy_qN3 = d_dy_qN2
! af 01/04/2016
            d_dx_qN3_z = d_dx_qN2_z
            d_dy_qN3_z = d_dy_qN2_z

            N2=NNFIRST
            qN2=qNFirst
! arash 12/14/2015
            Biharmonic_qN2 = Biharmonic_qNFirst
! arash 12/21/2015
            d_dx_qN2 = d_dx_qNFirst
            d_dy_qN2 = d_dy_qNFirst
! af 01/04/2016
            qN2_z = qNFirst_z
            Biharmonic_qN2_z = Biharmonic_qNFirst_z
            d_dx_qN2_z = d_dx_qNFirst_z
            d_dy_qN2_z = d_dy_qNFirst_z

            NEle=NeiTabEle(NH,NNeigh(NH)-1)

! arash Feb/18/2016
            velocity_magnitude = 0.0d0
c. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of the IF statement.
c.RJW 06092009  modified to elimiate zero index NOFF
        IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
            NCELE = NODECODE(NH)*NODECODE(N2)*
     &              NODECODE(N3)*NOFF(NELE)
            IF((NCELE.NE.0)) THEN  !if element is active, compute velocity grads
C    kmd48.33bc added in spherical factors for 3D
                  SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  qn=(qN1+qN2+qN3)/3.d0
                  Un=real(qn)
                  Vn=aimag(qn)

                  velocity_magnitude = sqrt( Un**2 + Vn**2 )

                  a1=X(N3)-X(N2)
                  a2=X(N1)-X(N3)
                  a3=X(N2)-X(N1)
                  b1=(Y(N2)-Y(N3))*SFacAvg
                  b2=(Y(N3)-Y(N1))*SFacAvg
                  b3=(Y(N1)-Y(N2))*SFacAvg
                  UnDqDX2A=UnDqDX2A+Un*(qN1*b1+qN2*b2+qN3*b3)
                  VnDqDY2A=VnDqDY2A+Vn*(qN1*a1+qN2*a2+qN3*a3)
! iso-sigma
!                  DqDXDPhiDX2A=DqDXDPhiDX2A
!     &                 +(qN1*b1+qN2*b2+qN3*b3)*b1/Areas(NEle)
!                  DqDYDPhiDY2A=DqDYDPhiDY2A
!     &                 +(qN1*a1+qN2*a2+qN3*a3)*a1/Areas(NEle)
! arash 12/14/2015
!                  Biharmonic_DqDXDPhiDX2A = Biharmonic_DqDXDPhiDX2A
!     &              +(Biharmonic_qN1*b1+Biharmonic_qN2*b2+Biharmonic_qN3*b3)*b1/Areas(NEle)
!                  Biharmonic_DqDYDPhiDY2A = Biharmonic_DqDYDPhiDY2A
!     &              +(Biharmonic_qN1*a1+Biharmonic_qN2*a2+Biharmonic_qN3*a3)*a1/Areas(NEle)
! arash 12/21/2015
!                  d2_dx2_q_element = d2_dx2_q_element + d_dx_qN1*b1+d_dx_qN2*b2+d_dx_qN3*b3
!                  d2_dy2_q_element = d2_dy2_q_element + d_dy_qN1*a1+d_dy_qN2*a2+d_dy_qN3*a3
!                  d2_dxy_q_element = d2_dxy_q_element + d_dy_qN1*b1+d_dy_qN2*b2+d_dy_qN3*b3
! z-level
                  DqDXDPhiDX2A=DqDXDPhiDX2A
     &               +(qN1_z*b1+qN2_z*b2+qN3_z*b3)*b1/Areas(NEle)
                  DqDYDPhiDY2A=DqDYDPhiDY2A
     &              +(qN1_z*a1+qN2_z*a2+qN3_z*a3)*a1/Areas(NEle)

                  Biharmonic_DqDXDPhiDX2A = Biharmonic_DqDXDPhiDX2A
     &              +(Biharmonic_qN1_z*b1+Biharmonic_qN2_z*b2+Biharmonic_qN3_z*b3)*b1/Areas(NEle)
                  Biharmonic_DqDYDPhiDY2A = Biharmonic_DqDYDPhiDY2A
     &              +(Biharmonic_qN1_z*a1+Biharmonic_qN2_z*a2+Biharmonic_qN3_z*a3)*a1/Areas(NEle)

                  d2_dx2_q_element = d2_dx2_q_element + d_dx_qN1_z*b1+d_dx_qN2_z*b2+d_dx_qN3_z*b3
                  d2_dy2_q_element = d2_dy2_q_element + d_dy_qN1_z*a1+d_dy_qN2_z*a2+d_dy_qN3_z*a3
                  d2_dxy_q_element = d2_dxy_q_element + d_dy_qN1_z*b1+d_dy_qN2_z*b2+d_dy_qN3_z*b3
               ENDIF
            ENDIF

            IF(TotalArea2.EQ.0.) THEN
               LAdvec(k)=(0.d0,0.d0)
               LStress(k)=(0.d0,0.d0)
! arash 12/14/2015
               Biharmonic_LStress(k)=(0.d0,0.d0)
            ELSE
               LAdvec(k)=sponge(NH)*IFNLCT*(UnDqDX2A+VnDqDY2A)/TotalArea2
               LStress(k)=3.d0*EVM(NH)*(DqDXDPhiDX2A+DqDYDPhiDY2A)
     &              /TotalArea2
! arash 12/21/2015: modified Leith viscosity
               ! 0. preliminaries
               modified_Leith_delta_x = sqrt ( TotalArea2 * 0.5d0 ) ! grid length
               modified_Leith_Lambda_vorticity  = 2.5 !1.4d0 !             ! tuning constant O(1) ========================================================
               modified_Leith_Lambda_divergence = 2.5 !1.4d0 !             ! the modification - see MITgcm manual ========================================
               ! 1. nodal second derivatives: during the computations, we multiplied them by 2A.
               d2_dx2_q_element = d2_dx2_q_element / TotalArea2
               d2_dy2_q_element = d2_dy2_q_element / TotalArea2
               d2_dxy_q_element = d2_dxy_q_element / TotalArea2
               ! 2. magnitude of the gradient of divergence and vorticity
               modified_Leith_grad_vorticity = sqrt( ( aimag(d2_dx2_q_element) - real(d2_dxy_q_element) )**2
     &                                       +       ( aimag(d2_dxy_q_element) - real(d2_dy2_q_element) )**2 )
               modified_Leith_grad_divergence= sqrt( (  real(d2_dx2_q_element) +aimag(d2_dxy_q_element) )**2
     &                                       +       (  real(d2_dxy_q_element) +aimag(d2_dy2_q_element) )**2 )
               ! 3. compute the modified Leith viscosity
               Biharmonic_viscosity =   modified_Leith_delta_x**5 / ( 8.0d0 * 3.141592654d0**3 ) *
     &                           sqrt ( modified_Leith_Lambda_vorticity**6  * modified_Leith_grad_vorticity**2
     &                                + modified_Leith_Lambda_divergence**6 * modified_Leith_grad_divergence**2 )
               ! 3.1. engineering
               ! 3.1.1. minimum viscosity
               If ( Biharmonic_viscosity < 1.0d-2 ) Biharmonic_viscosity = 1.0d-2

               ! 3.1.2. maximum viscosity 
               element_peclet_number = velocity_magnitude * length / Biharmonic_viscosity
!               If ( element_peclet_number > 2.0 ) Then
!                  Biharmonic_viscosity = velocity_magnitude * length / 2.0d0
!               End If
               ! 3.9. save the viscosity parameter for using in the transport equation
                Biharmonic_viscosity_modified_Leith ( NH, k ) = Biharmonic_viscosity
               !write(6,'("nh, k, nele, ncele, biharmonic_viscosity_modified_leith",i0,1x,i0,1x,i0,1x,i0,f8.2)'),nh, k, nele, ncele, Biharmonic_viscosity_modified_Leith(nh,k) !jgfdebug
               ! 4. congratulations! we have computed modified Leith viscosity!
               Biharmonic_LStress(k) = -3.d0 * Biharmonic_viscosity * (Biharmonic_DqDXDPhiDX2A+Biharmonic_DqDYDPhiDY2A)
     &              /TotalArea2
            ENDIF

         ENDDO

! |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
! arash 12/14/2015: Use Biharmonic viscosity instead of a harmonic one
! We may want to make this cleaner later by having an option in the fort.15 file
         LStress(:) = Biharmonic_LStress(:)
! |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

!--------------------------------------------------------------------------------------------------
c     Finished computing advection and stress terms in the momentum
c     equation at each level in the vertical

! arash May 31 2016 commented (!!) this based on Rosemary's work
!---------------------------------------------------------------------------------------------------
c     Compute the equivalent linear slip coefficient
C RJW 02252009 changed KP for ISlip =2
C     For ISlip = 1, KP is the equivalent linear slip coefficient
C     For ISlip = 2, KP is the minimum equivalent quadratic slip coefficient
C     For ISlip = 3, KP is the quadratic slip coefficient
!!         IF(ISlip.EQ.1) KSlip=KP
!!         IF(ISlip.EQ.2) THEN
c.RJW merged 08282008 Casey 071219: Let's get funky.
c.RJW enable this modification
C     and add the following to base the drag coefficient on the actual
C     physics of the boundary layer.
C ---

Casey 140701: Adjust for small depths.
!!            Hsp1=DP(NH)+IFNLFA*ETA2(NH)
!!            IF(Hsp1.LT.H0) Hsp1=H0
!            IF (ABS( ( ( SIGMA(2)-SIGMA(1) )/2 )
!     &                *(DP(NH)+IFNLFA*ETA2(NH)) )
! arash
!!            IF ( (ABS( ( ( SIGMA(2)-SIGMA(1) )/2.0d0 ) *Hsp1 ))
!!     &                        .LT. 1.D0 )THEN

!!                KSlip = (1.D0 / ( (1.D0/0.41D0) *
!!     &                  LOG( (ABS( ( ( SIGMA(2)-SIGMA(1) )/2.0d0 )
!     &                     *( DP(NH)+IFNLFA*ETA2(NH) )) + Z0B1 )
! arash
!!     &                     * Hsp1 ) + Z0B1 )
!!     &                  / Z0B1 ) ) )**2.D0
!!     &                * ABS(Q(NH,1))
!!            ELSE
!!                KSlip = (1.D0 / ( (1.D0/0.41D0) *
!!     &                  LOG( (1.D0+Z0B1)
!!     &                  / (Z0B1) ) ) )**2.D0
!!     &                * ABS(Q(NH,1))
!!            ENDIF
!!            IF(KSlip.GT.1.0d0* ABS(Q(NH,1))) KSlip=1.0d0*ABS(Q(NH,1))
C          IF(KSlip.LT.0.0010d0*ABS(Q(NH,1))) KSlip=0.0010d0*ABS(Q(NH,1))
!!            IF(KSlip.LT.KP*ABS(Q(NH,1))) KSlip=KP*ABS(Q(NH,1))

!!         ENDIF
!!         IF(ISlip.EQ.3) THEN
!!            KSlip=KP*ABS(Q(NH,1))
!!            IF(KSlip.GT.0.1d0) KSlip=0.1d0
!!            IF(KSlip.LT.0.0d0) KSlip=0.0d0
!!         ENDIF
!---------------------------------------------------------------------------------------------------
c     Set up the RHS forcing vector Fr and LHS matrix in compact storage
c     (Mkm1,Mk,Mkp1)

         CCL = 1.d0+Corif(NH)*IDTAlp1
         CCR = 1.d0-Corif(NH)*IDT1MAlp1
         RCL = DTAlp3/Hsp1Hsp1OAMBAMB

         IF(ISlip.EQ.0) THEN    !no slip bottom boundary condition
            Fr(1)   = (0.d0,0.d0)
            Mkm1(1) = (0.d0,0.d0)
            Mk(1)   = (1.d0,-1.d0) !note: -I*IV=V
            Mkp1(1) = (0.d0,0.d0)
         ELSE                   ! slip bottom boundary condition
            Fr(1) = (CCR*q(NH,1)
     &           -DelT*(LAdvec(1)  +LStress(1)  +BPG(NH,1)))*Inm(1,2)
     &           + (CCR*q(NH,2)
     &           -DelT*(LAdvec(2)  +LStress(2)  +BPG(NH,2)))*Inm(1,3)
     &           - DelT*(VAdvec(1)+BTPG*LVn(1))-DT1MAlp3*VStress(1)
! arash May 31 2016 did this based on Rosemary's work 
!     &           - q(NH,1)*DT1MAlp2*KSlip/HsOAMB
     &           - q(NH,1)*DT1MAlp2*TK(NH)/HsOAMB

            Mkm1(1) = (0.d0,0.d0)
            Mk(1)   = CCL*Inm(1,2) + RCL*KVnm(1,2)
! arash May 31 2016 did this based on Rosemary's work 
!     &           + DTAlp2*KSlip/Hsp1OAMB
     &           + DTAlp2*TK(NH)/Hsp1OAMB
            Mkp1(1) = CCL*Inm(1,3) + RCL*KVnm(1,3)
         ENDIF

         DO k=2,NFEN-1
            Fr(k) = (CCR*q(NH,k-1)
     &           -DelT*(LAdvec(k-1)+LStress(k-1)+BPG(NH,k-1)))*Inm(k,1)
     &           + (CCR*q(NH,k)
     &           -DelT*(LAdvec(k)  +LStress(k)  +BPG(NH,k)  ))*Inm(k,2)
     &           + (CCR*q(NH,k+1)
     &           -DelT*(LAdvec(k+1)+LStress(k+1)+BPG(NH,k+1)))*Inm(k,3)
     &           - DelT*(VAdvec(k)+BTPG*LVn(k))-DT1MAlp3*VStress(k)
            Mkm1(k) = CCL*Inm(k,1) + RCL*KVnm(k,1)
            Mk(k)   = CCL*Inm(k,2) + RCL*KVnm(k,2)
            Mkp1(k) = CCL*Inm(k,3) + RCL*KVnm(k,3)
         ENDDO
C RJW 09/01/2009 RJW k = NFEN-1
C            try to fix this by setting k=NFEN
C         k=NFEN
         Fr(NFEN) = (CCR*q(NH,k-1)
     &        -DelT*(LAdvec(k-1)+LStress(k-1)+BPG(NH,k-1)))*Inm(k,1)
     &        + (CCR*q(NH,k)
     &        -DelT*(LAdvec(k)  +LStress(k)  +BPG(NH,k)  ))*Inm(k,2)
     &        - DelT*(VAdvec(k)+BTPG*LVn(k))-DT1MAlp3*VStress(k)
     &        + DelT*0.5d0*((WSX2(NH)+iy*WSY2(NH))/Hsp1OAMB
     &        +(WSX1(NH)+iy*WSY1(NH))/HsOAMB)
         Mkm1(NFEN) = CCL*Inm(NFEN,1) + RCL*KVnm(NFEN,1)
         Mk(NFEN)   = CCL*Inm(NFEN,2) + RCL*KVnm(NFEN,2)
         Mkp1(NFEN) = (0.d0,0.d0)

c     Start section to modify equations depending on normal flux
c     boundary condition
c     0 <= LBcodeI <= 10, essential normal flux and free tangential slip
c     this b.c. is taken care of in the code section below
c     10 <= LBcodeI <= 19, essential normal flux and zero tangential slip
c     this b.c. is taken care of above
c     20 <= LBcodeI <= 29, natural normal flux and free tangential slip
c     this b.c. requires on manipulation of momentum eqns.  Do nothing!

         LBP=LBArray_Pointer(NH)
         IF(LBP.GT.0) THEN      !flux boundary
            IF((LBCODEI(LBP).GE. 0).AND.(LBCODEI(LBP).LE. 9)) THEN
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
               Frstar     =Real(Fr(1))*SLBP-AImag(Fr(1))*CLBP
     &              +  (Auv2k+Auv2kp1)*VelNorm
     &              +iy*(Auv1k+Auv1kp1)*VelNorm
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
                  Frstar     =Real(Fr(k))*SLBP-AImag(Fr(k))*CLBP
     &                 +  (Auv2km1+Auv2k+Auv2kp1)*VelNorm
     &                 +iy*(Auv1km1+Auv1k+Auv1kp1)*VelNorm
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
               Frstar     =Real(Fr(NFEN))*SLBP-AImag(Fr(NFEN))*CLBP
     &              +  (Auv2km1+Auv2k)*VelNorm
     &              +iy*(Auv1km1+Auv1k)*VelNorm
               Fr(NFEN)   =Frstar
            ENDIF
         ENDIF

c     Finished section to modify equations depending on normal flux
c     boundary condition
c     Decompose and solve the system


 999     CALL TRIDIAG(Mkm1,Mk,Mkp1,Fr,Gamma,NFEN)

C     Changed this q value to qkp1 since we need to save q and qkp1 for
C     transport calculations
C     Save the horizontal velocity and Eddy Viscosity solutions
C   kmd48.33bc added the vertical diffusion term for the transport equation

         DO k=1,NFEN
            qkp1(NH,k) = Gamma(k)    ! arash: this seems to be the velocity at time step k+1
            EV(NH,k)=EVTot(k)
            DV(NH,k)=NTVTot(k)
         END DO

      ENDDO

c     Finished computing horizontal velocities

! -------------------------------------------------------------------------------------------------
!                  G a u s s i a n    F i l t e r
! -------------------------------------------------------------------------------------------------
#ifdef CMPI
      CALL UPDATEC3D(Qkp1)
#endif
! arash: apply a 1-layer Gaussian filter to the horizontal velocities
       Call Gaussian_filter_complex ( X, Y, NeiTab, NNeigh, qkp1, NP, NFEN, MNEI )
! -------------------------------------------------------------------------------------------------


C****************************************************************************
C     Update horizontal velocity solution on ghost nodes

#ifdef CMPI
      CALL UPDATEC3D(Qkp1)
#endif


C****************************************************************************
C     Compute the depth averaged velocity, depth averaged flux, bottom stress and
C     velocity dispersion

      DO NH=1,NP                                !loop over horizontal nodes
C........DMW202401 Use saved H2
         Hsp1= H2(NH)           !Total depth at present (s+1) timestep
         VIVel = (0.d0,0.d0)
         DO k=1,NFEN
            VIVel = VIVel + Qkp1(NH,k)*LVn(k)
         END DO
         VIVel = VIVel/amb
         UU(NH) = REAL(VIVel)
         VV(NH) = AIMAG(VIVel)

         DAFluxX(NH)=UU(NH)*Hsp1
         DAFluxY(NH)=VV(NH)*Hsp1

         IF(ISlip.EQ.0) THEN
            DUDS = (Qkp1(NH,2)-Qkp1(NH,1))/(Sigma(2)-Sigma(1))
            BSX(NH) = EV(NH,1)*REAL(DUDS)
            BSY(NH) = EV(NH,1)*AIMAG(DUDS)
         ENDIF
         IF(ISlip.NE.0) THEN
! arash May 31 2016 did this based on Rosemary's work 
!            BSX(NH) = KSlip*REAL(Qkp1(NH,1))
!            BSY(NH) = KSlip*AIMAG(Qkp1(NH,1))
            BSX(NH) = TK(NH)*REAL(Qkp1(NH,1))
            BSY(NH) = TK(NH)*AIMAG(Qkp1(NH,1))
         ENDIF

         DUU(NH) = 0.d0
         DUV(NH) = 0.d0
         DVV(NH) = 0.d0

Corbitt 120328: Applies Advection Locally
         IF (LoadAdvectionState) Then
            IF (DP(NH).GE.AdvectionState(NH)) THEN
               IFNLCT = IFNLCTE
            ELSE
               IFNLCT = 0
            ENDIF
         ENDIF

         IF(IFNLCT.EQ.1) THEN
            CALL VSDISP (IT,NH,Hsp1,UU(NH),VV(NH),
     &           DUU(NH),DUV(NH),DVV(NH))
         ENDIF

      ENDDO

C***************************************************************************
C     Compute "z" vertical velocity

C
C     Loop over each horizontal node to compute the "z" version of the
C     vertical velocity
C
      DO NH=1,NP

c. RJW merged 08282008 Casey 071219: Solve for TEMPNCELE, same as above.
C
         TEMPNCELE = 0
         TEMPSTOP = 0
Casey 140701: Debug.
C        DO K=1,NODELE(NH)
         DO K=1,SIZE(NEITABELE,2)
            IF(NEITABELE(NH,K).EQ.0)THEN
C              TEMPSTOP = TEMPSTOP + 1
               CYCLE
            ENDIF
            TEMPNCELE = TEMPNCELE +
C    &         NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),1))*
C    &         NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),2))*
C    &         NODECODE(NM(NEITABELE(NH,K+TEMPSTOP),3))*
C    &         NOFF(NEITABELE(NH,K+TEMPSTOP))
     &         NODECODE(NM(NEITABELE(NH,K),1))*
     &         NODECODE(NM(NEITABELE(NH,K),2))*
     &         NODECODE(NM(NEITABELE(NH,K),3))*
     &         NOFF(NEITABELE(NH,K))
         ENDDO

         IF(TEMPNCELE.EQ.NODELE(NH))THEN
            TEMPNCELE = 1
         ELSEIF(TEMPNCELE.EQ.0)THEN
            TEMPNCELE = 0
         ELSE
            TEMPNCELE = 1
         ENDIF

c     Set up some values at the node being worked on
C........DMW202401 Use saved H2
         Hsp1=H2(NH)
Casey 140701: Adjust for small depths.
         IF(Hsp1.LT.H0) Hsp1=H0
         Hsp1OAMB=Hsp1/AMB

c     Compute time derivative of water surface position

         DEtaDT=(Eta2(NH)-Eta1(NH))/DelT

c     Compute horizontal derivatives of water surface position,
c     bathymetric depth and horizontal velocity Note: TotalArea2 = 2X
c     total elemental area surrounding a node

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
            NEle=NeiTabEle(NH,N-2) !element # defined by nodes NH,NN2,NN1
c. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of this IF statement.
            !jgf48.50 Avoid array out of bounds
            IF (NEle.eq.0) THEN
               CYCLE
            ENDIF
            NCELE = NODECODE(NH)*NODECODE(N2)*
     &              NODECODE(N3)*NOFF(NELE)
            IF(NCELE.NE.0) THEN  !if element is active, compute velocity grads
               TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C    kmd48.33bc added in spherical factors for 3D
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
c. RJW merged 08282008 Casey 071219: Added the computation of NCELE and the second half of this IF statement.
         !jgf48.50 Avoid array out-of-bounds
         IF (NEle.ne.0 ) THEN
            NCELE = NODECODE(NH)*NODECODE(N2)*
     &           NODECODE(N3)*NOFF(NELE)
            IF (NCELE.NE.0) THEN  !if element is active, compute velocity grads
               TotalArea2=TotalArea2+Areas(NEle) !accumulate 2X total areas to complete calc.
C    kmd48.33bc added in spherical factors for 3D
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

         IF(TotalArea2.NE.0.) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
            DO k=1,NFEN
               DUDX(k)=REAL(DqDX2A(k))/TotalArea2
               DVDY(k)=AIMAG(DqDY2A(k))/TotalArea2
            ENDDO
         ENDIF

c     Evaluate the "z" vertical velocity

c. RJW merged 08282008 Casey 071219: Added an IF statement so that the vertical velocities are only computed
C     and corrected at nodes that are connected to at least one active element.
C
         IF(TEMPNCELE.NE.0)THEN
c.
         WZkp1(NH,1)=-REAL(qkp1(NH,1))*DhDX-AIMAG(qkp1(NH,1))*DhDY

         DO k=2,NFEN
            DelSigO2=(Sigma(k)-Sigma(k-1))/2.d0
            SigAvgMAOAMB=((Sigma(k)+Sigma(k-1))/2.d0 - A)/AMB
            DelU=REAL(qkp1(NH,k)-qkp1(NH,k-1))
            DelV=AIMAG(qkp1(NH,k)-qkp1(NH,k-1))
            WZkp1(NH,k)=WZkp1(NH,k-1)
     &      - DelSigO2*Hsp1OAMB*(DUDX(k)+DVDY(k)+DUDX(k-1)+DVDY(k-1))
     &           + (DEtaDX+SigAvgMAOAMB*(DhDX+DEtaDX))*DelU
     &           + (DEtaDY+SigAvgMAOAMB*(DhDY+DEtaDY))*DelV
         ENDDO

c     Correct this using Adjoint method

         Wf=0.d0                !This value should match surface B.C. exactly
C........DMW202401 Use saved H2
         Hsp1= H2(NH) !Total depth at present (s+1) timestep
         WfOHH=Wf/(Hsp1*Hsp1)
         WZSurfBC=DEtaDT+REAL(qkp1(NH,NFEN))*DEtaDX
     &                  +AIMAG(qkp1(NH,NFEN))*DEtaDY
         WZSurf=WZkp1(NH,NFEN)

         DO k=1,NFEN
            WZCorrection=(WZSurfBC-WZSurf)*(WfOHH+(Sigma(k)-b)/AMB)
     &           /(2.d0*WfOHH+1)
            WZkp1(NH,k)=WZkp1(NH,k)+WZCorrection
         ENDDO

c.RJW merged 08282008 Casey 071219: Added an ELSEIF statement to zero out the vertical velocities at nodes
C     that are not connected to any active elements.
C
         ELSEIF(TEMPNCELE.EQ.0)THEN
            DO K=1,NFEN
               WZkp1(NH,K) = (0.D0,0.D0)
            ENDDO
         ENDIF
c.
C
C     End loop over horizontal nodes to compute vertical velocity
C
      ENDDO

C***************************************************************************
C     Update vertical velocity on ghost nodes - this is necessary

#ifdef CMPI
      CALL UPDATER3D(WZkp1)
#endif


C**************************************************************************
C     For a prognostic baroclinic run, compute new density, temperature, salinity fields
C
c   kmd48.33bc need to update the boundary condition information for both
C              salinity and temperature in prognostic runs. The boundary
C              conditions are read in similar to the aperiodic elevation
C              boundary conditions
      IF ((CBAROCLINIC).AND.(RES_BC_FLAG.GT.0).AND.(NOPE.GT.0)) THEN
         IF ((ABS(RES_BC_FLAG).EQ.2).OR.(ABS(RES_BC_FLAG).EQ.4)) THEN
            IF(TimeLoc.GT.SBCTIME2) THEN
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
                  RESSAL(NumofBCNode,J)=RESSAL1(NumofBCNode,J)+
     &            SBCRATIO*(RESSAL2(NumofBCNode,J)-
     &            RESSAL1(NumofBCNode,J))
               END DO
            END DO
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &           (ABS(RES_BC_FLAG).EQ.4)) THEN
            IF(TimeLoc.GT.TBCTIME2) THEN
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
                  RESTEMP(NumofBCNode,J)=RESTEMP1(NumofBCNode,J)+
     &            TBCRATIO*(RESTEMP2(NumofBCNode,J)-
     &            RESTEMP1(NumofBCNode,J))
               END DO
            END DO
         END IF
      END IF

!kmd - adding information for rivers - baroclinic

      IF ((CBaroclinic).AND.(BndBCRiver)) THEN
         IF (TIMELOC.GT.RIVBCTIME2) THEN
            RIVBCTIME1=RIVBCTIME2
            RIVBCTIME2=RIVBCTIME1+RIVBCTIMINC
            DO J=1,totalbcrivernodes
               IF (IDEN.EQ.2) THEN
                  DO K=1,NFEN
                     BCRivSalN1(J,K)=BCRivSalN2(J,K)
                  END DO
Casey 140701: Shouldn't read the vertex number.
C                 READ(39,*) NOD, (BCRivSalN2(J,K),K=1,NFEN)
                  READ(39,*) (BCRivSalN2(J,K),K=1,NFEN)
               ELSE IF (IDEN.EQ.3) THEN
                  DO K=1,NFEN
                     BCRivTempN1(J,K)=BCRivTempN2(J,K)
                  END DO
Casey 140701: Shouldn't read the vertex number.
C                 READ(39,*) NOD, (BCRivTempN2(J,K),K=1,NFEN)
                  READ(39,*) (BCRivTempN2(J,K),K=1,NFEN)
               ELSE IF (IDEN.EQ.4) THEN
                  DO K=1,NFEN
                     BCRivSalN1(J,K)=BCRivSalN2(J,K)
                     BCRivTempN1(J,K)=BCRivTempN2(J,K)
                  END DO
Casey 140701: Shouldn't read the vertex number.
C                 READ(39,*) NOD, (BCRivSalN2(J,K),BCRivTempN2(J,K),K=1,NFEN)
                  READ(39,*) (BCRivSalN2(J,K),BCRivTempN2(J,K),K=1,NFEN)

                  ! arash debug -----------------------------------------------------------------------------
!                  BCRivSalN1(:,:) = 0.0
!                  BCRivTempN1(:,:) = 0.0
!                  BCRivSalN2(:,:) = 0.0
!                  BCRivTempN2(:,:) = 0.0

               END IF
            END DO
         END IF
         BCRivRATIO=(TIMELOC-RIVBCTIME1)/RIVBCTIMINC
         DO J=1,totalbcrivernodes
             IF (IDEN.EQ.2) THEN
                 DO K=1,NFEN
                    BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO*
     &                 (BCRivSalN2(J,K)-BCRivSalN1(J,K))
                 END DO
             ELSE IF (IDEN.EQ.3) THEN
                 DO K=1,NFEN
                    BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO*
     &                 (BCRivTempN2(J,K)-BCRivTempN1(J,K))
                 END DO
             ELSE IF (IDEN.EQ.4) THEN
                 DO K=1,NFEN
                    BCRivSal(J,K)=BCRivSalN1(J,K)+BCRivRATIO*
     &                 (BCRivSalN2(J,K)-BCRivSalN1(J,K))
                    BCRivTemp(J,K)=BCRivTempN1(J,K)+BCRivRATIO*
     &                 (BCRivTempN2(J,K)-BCRivTempN1(J,K))
                 END DO
             END IF
         END DO
      END IF

C   Start computing the new density, temperature and salinity values
      IF(C3D_BTrans) THEN
        IF(IDEN.EQ.1) THEN
C           CALL 3D_SIGMAT_TRANS()
        ELSE IF(IDEN.EQ.2) THEN
            CALL TRANS_3D(SAL,NLSD,NVSD,SALkp1,IDEN,TimeLoc)
        ELSE IF(IDEN.EQ.3) THEN
            CALL TRANS_3D(TEMP,NLTD,NVTD,TEMPkp1,IDEN,TimeLoc)
        ELSE IF(IDEN.EQ.4) THEN
            CALL TRANS_3D(SAL,NLSD,NVSD,SALkp1,IDEN-2,TimeLoc)
            CALL TRANS_3D(TEMP,NLTD,NVTD,TEMPkp1,IDEN-1,TimeLoc)
        ENDIF

! arash: clip negative salinity values and 04/01/2016 apply a Gaussian filter
!-----------------------------------------------------------------------
! clip
! salinity at the Dead Sea is 40 psu, 35.5 psu at global ocean, and 15 psu at the mouth of rivers.
        DO NH=1,NP
           DO N=1,NFEN
              IF (SALkp1(NH,N) .LT. 0.0d0) THEN
                  SALkp1(NH,N) = 0.0d0
              END IF
              IF (SALkp1(NH,N) .GT. 40.0d0) THEN  
                  SALkp1(NH,N) = 40.0d0
              END IF
           END DO
        END DO

! limit temperature values; very basic lower and upper bounds for now.
        DO NH=1,NP
           DO N=1,NFEN
              IF (TEMPkp1(NH,N) .LT. 0.0d0) THEN
                  TEMPkp1(NH,N) = 0.0d0
              END IF
              IF (TEMPkp1(NH,N) .GT. 100.0d0) THEN  
                  TEMPkp1(NH,N) = 100.0d0
              END IF
           END DO
        END DO


! arash added this - not sure if necessary
#ifdef CMPI
C      CALL UPDATER3D(SIGT) !jgf45.12
      IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) CALL UPDATER3D(Tempkp1)
      IF((IDEN.EQ.2).OR.(IDEN.EQ.4)) CALL UPDATER3D(Salkp1)
#endif

! Gausian filter
      IF(IDEN.EQ.1) THEN
         ELSE IF(IDEN.EQ.2) THEN
            ! salinity:
!            Call Gaussian_filter_real ( X, Y, NeiTab, NNeigh, SALkp1, NP, NFEN, MNEI )
         ELSE IF(IDEN.EQ.3) THEN
            ! temperature:
!            Call Gaussian_filter_real ( X, Y, NeiTab, NNeigh, TEMPkp1, NP, NFEN, MNEI )
         ELSE IF(IDEN.EQ.4) THEN
            ! salinity & temperature:
!            Call Gaussian_filter_real ( X, Y, NeiTab, NNeigh, SALkp1,  NP, NFEN, MNEI )
!            Call Gaussian_filter_real ( X, Y, NeiTab, NNeigh, TEMPkp1, NP, NFEN, MNEI )
      END IF
!-----------------------------------------------------------------------
C
C     Update new temperature and salinity fields on ghost nodes
C
#ifdef CMPI
C      CALL UPDATER3D(SIGT) !jgf45.12
      IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) CALL UPDATER3D(Tempkp1)
      IF((IDEN.EQ.2).OR.(IDEN.EQ.4)) CALL UPDATER3D(Salkp1)
#endif

C     Need to update and move results from the future time levels to the
C     present time levels so we can advance in time and before obtaining
C     new sigma-t values. We will update velocities later.
        DO NH = 1,NP
           DO k = 1,NFEN
              IF ((IDEN.EQ.2).OR.(IDEN.EQ.-2)) THEN
                  SAL(NH,k)=SALkp1(NH,k)
              ELSE IF ((IDEN.EQ.3).OR.(IDEN.EQ.-3)) THEN
                  TEMP(NH,k)=TEMPkp1(NH,k)
              ELSE IF ((IDEN.EQ.4).OR.(IDEN.EQ.-4)) THEN
                  SAL(NH,k)=SALkp1(NH,k)
                  TEMP(NH,k)=TEMPkp1(NH,k)
              END IF
           END DO
        END DO
        IF(IDEN.GT.1) CALL CALC_SIGMAT_3D ()

C     Compute Depth averaged SigmaT

        DO NH=1,NP
         DASigT(NH)=0.d0
          DO k=1,NFEN-1
            DASigT(NH)=DASigT(NH)+(Sigma(k+1)-Sigma(k))
     &                                   *(SigT(NH,k+1)+SigT(NH,k))/2.d0
          ENDDO
          DASigT(NH)=DASigT(NH)/AMB
        ENDDO
      ENDIF

C     Need to update and move results from the future time levels to the
C     present time levels so we can advance in time.
      DO NH = 1,NP
         DO k = 1,NFEN
            ! arash 2016 Mar 22
 !           qkm1(NH,k)  = q(NH,k)
            q(NH,k)     = qkp1(NH,k)
!            WZkm1(NH,k) = WZ(NH,k)
            WZ(NH,k)    = WZkp1(NH,k)
         END DO
      END DO

! arash 2016 Mar 22: apply an Asselin filter (see NCOM manual page 18)
!      If ( time_step_counter == 0 )
         ! do not apply the filter if you are in the first time step
!      Else
!         time_step_counter = time_step_counter + 1
         ! apply the Asselin filter with nu = 0.05
!         q = 0.05d0 * ( qkp1 + qkm1 )

!      End If

C     WRITE 3D MODEL OUTPUT
!--------------------------------------------------------------------------------------------------
      CALL writeOutput3D(TimeLoc,IT) !jgf48_11 moved to write_output.F
!--------------------------------------------------------------------------------------------------

C     RETURN TO THE 2-D MODEL TO COMPUTE H AT NEXT TIME STEP

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C***************************************************************************
      END SUBROUTINE VSSOL
C***************************************************************************




c******************************************************************************
c     Subroutine to compute the eddy viscosity profile.                       *
c                                                                             *
c  ievc, evmin, evcon - E.V. code, E.V. minimum value and E.V. constant       *
c                                                                             *
c        NOTE: evcon only used for some of the E.V. formulations as           *
c                  discussed below.                                           *
c        NOTE: In cases where EV is specified to vary linearly over the       *
c              lower 20% of the water column, it actually varies linearly     *
c              with a constant slope up to the vertical FE grid node that is  *
c              less than or equal to the 20% location.  The value is constant *
c              as specified at all FE grid nodes above the 20% location.      *
c              The E.V. above and below the 20% level is joined by one        *
c              additional linearly varying segment.                           *
c        NOTE: The E.V. is constrained to always be greater than or equal to  *
c              EVMIN as specified in the UNIT 15 file.                        *
c                                                                             *
c        ievc=0-9, EV constant in time & horizontal space                     *
c             0 - EV read in from UNIT 15 (may vary vertically) - EVCON is    *
c                    not used                                                 *
c             1 - EV = EVCON                                                  *
c                                                                             *
c        ievc=10-19 EV proportional to omega*h*h  (Lynch and Officer (1986)   *
c                                              Lynch and Werner (1987, 1991)) *
c             10 - EV = omega*h*h/10 over the entire water column             *
c             11 - EV = omega*h*h/1000 at bottom                              *
c                       varies linear over lower 20% of wc                    *
c                     = omega*h*h/10 in upper 80% of w.c.                     *
c            NOTE:For this EV formulation, evcon is not used and omega is     *
c                  hardwired for a 12.42 hour tide.                           *
c                                                                             *
c        ievc=20-29 EV proportional to kappa U* z                             *
c             20 - EV = 0.41U*Zo at bottom                                    *
c                     = 0.41U*Z over entire water column                      *
c             21 - EV = 0.41U*Zo at bottom                                    *
c                     = 0.41U*Z in lower 20% of water col                     *
c                     = 0.082U*h in upper 80% of water col                    *
c            WHERE: U* is the friction velocity                               *
c            NOTE: For this EV formulation, evcon is not used.                *
c                                                                             *
c        ievc=30-39, EV proportional to Uh (Davies 1990)                      *
c             30 - EV = 0.025|U|h/9.001 over entire water column              *
c             31 - EV = evcon|U|h over entire water column                    *
c             32 - EV = 0.025|U|h/9.001 in upper 80% of wc                    *
c                     = 0.000025h|U|/9.001 at bottom                          *
c                       varies linear over lower 20% of wc                    *
c             33 - EV = evcon|U|h in upper 80% of wc                          *
c                     = evcon|U|h/1000. at bottom                             *
c                       varies linear over lower 20% of wc                    *
c            WHERE: U is depth averaged velocity                              *
c            NOTE: For this EV formluation, evcon is used only for ievc=31,33 *
c                                                                             *
c        ievc=40-49, EV proportional to U*U (Davies 1990)                     *
c             40 - EV = 2|UU|/9.001 over entire water column                  *
c             41 - EV = evcon|UU| over entire water column                    *
c             42 - EV = 2|UU|/9.001 in upper 80% of wc                        *
c                     = 0.002|UU|/9.001 at bottom                             *
c                       varies linear over lower 20% of wc                    *
c             43 - EV = evcon|UU| in upper 80% of wc                          *
c                     = evcon|UU|/1000. at bottom                             *
c                       varies linear over lower 20% of wc                    *
c            WHERE: U is depth averaged velocity                              *
c            NOTE: For this EV formluation, evcon is used only for ievc=41,43 *
c                                                                             *
c        ievc=50, EV computed from Mellor-Yamada L2.5 closure                 *
c            NOTE: For this EV formulation, evcon is not used.                *
c            =51, EV computed from Mellor-Yamada L2.5 closure                 *
C                 with enhanced mixing in the surface layer                   *
c            NOTE: For this EV formulation, evcon is the value                *
C                  alpha coefficient to the surface roughness equation        *
C                  The TKE production parameter (a) is hardwired to           *
C                  a value of 150                                             *
C                                                                             *
C                               04/06/11                                      *
C******************************************************************************
C
       SUBROUTINE EDDYVIS(H,UU,VV,WSX,WSY,BSX,BSY,NODE,IT,Z0B1)
Cjgfdebug USE GLOBAL, ONLY : ScreenUnit
Cjgfdebug  USE GLOBAL_3DVS, ONLY : LOCALDIR,IEVC,EVMIN,EVCON,EVTOT,SIGMA,
Cjgfdebug     &     NFEN,Z0B,Z0S,A,B,AMB,NSCREEN,NWS,DELT
      USE GLOBAL_3DVS, ONLY : IEVC,EVMIN,EVCON,EVTOT,SIGMA,
     &     NFEN,Z0B,Z0S,A,B,AMB
      use mod_logging, only: setMessageSource, unsetMessageSource,
     &     logMessage, allMessage, DEBUG, ECHO, INFO, WARNING, ERROR
Crjwdebug  USE NodalAttributes, ONLY : Z0B_var
      IMPLICIT NONE
      REAL(8) :: H
      REAL(8) :: UU,VV,WSX,WSY,BSX,BSY
      REAL(8) :: Z0B1

C     jgf46.00 explicitly declared the following variables
      INTEGER NODE
      INTEGER, SAVE :: istart = 0
      INTEGER IT
      INTEGER J
      REAL(8) RKAPPA
      REAL(8) OMEGA
      REAL(8) EVBASE
      REAL(8) SLOPE
      REAL(8) BREAK
      REAL(8) EVBEGIN
      REAL(8) USTARB
      REAL(8) USTARS
      REAL(8) UAVMAG
      REAL(8) UAVMAGS

      call setMessageSource("eddyvis")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

Cjgfdebug 350  FORMAT(//,2X,'***** INVALID INPUT IN THE PRIMARY VERTICAL INPUT',
Cjgfdebug     &     ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

      RKAPPA = 0.41

C...
C...  IEVC=0 READ IN FROM UNIT 15 in subroutine READ_INPUT_3DVS
C...

C...
C...  IEVC=1 SET EV = EVCON ON THE FIRST TIME STEP
C     .

      IF((IEVC.eq.1).and.(istart.eq.0)) THEN
         istart = 1
         DO J=1,NFEN
            EVTOT(J)=EVCON
         ENDDO
      ENDIF

      IF(IEVC.eq.1) GOTO 100

C...
C...  OMEGA*H*H formulation FOLLOWING
C...  LYNCH AND OFFICER (1986), LYNCH AND WERNER (1987, 1991)
C...
      IF(IEVC.EQ.10) THEN
         OMEGA=1.40525D-4
         EVBASE=OMEGA*H*H/10.D0
         IF(EVBASE.LT.EVMIN) EVBASE=EVMIN
         DO J=1,NFEN
            EVTOT(J) = EVBASE
         END DO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.11) THEN
         OMEGA=1.40525D-4
         EVBASE=OMEGA*H*H/10.D0
         EVTOT(1)=EVBASE/100.D0
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-b*SLOPE
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF

C...
C...  KAPPA USTAR Z FORMULATION
C...
      IF(IEVC.EQ.20) THEN
         USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))
         EVTOT(1)=RKAPPA*USTARB*Z0B1
         EVBASE=RKAPPA*USTARB*H
         IF(EVBASE.LT.EVTOT(1)) EVBASE=EVTOT(1)
         SLOPE=(EVBASE-EVTOT(1))/amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J) !=EVTOT(1)+SLOPE*(SIGMA(J)-b)
         ENDDO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.21) THEN
         USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))
         EVTOT(1)=RKAPPA*USTARB*Z0B1
         EVBASE=RKAPPA*USTARB*H*0.2D0
         IF(EVBASE.LT.EVTOT(1)) EVBASE=EVTOT(1)
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.22) THEN
         USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))

         USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
         EVTOT(1)=RKAPPA*0.5*(USTARB+USTARS)*Z0B1
         EVBASE=RKAPPA*0.5*(USTARB+USTARS)*H*0.2D0
         IF(EVBASE.LT.EVTOT(1)) EVBASE=EVTOT(1)
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.23) THEN
          USTARB=SQRT(SQRT(BSX*BSX+BSY*BSY))

          USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
          EVTOT(1)=0.1*0.5*(USTARB+USTARS)*H
         DO J=2,NFEN
           EVTOT(J)=0.1*0.5*(USTARB+USTARS)*H
         ENDDO
         GOTO 100
      ENDIF

C...
C...  H UAVG FORMULATION FOLLOWING DAVIES (1990) + [EQ. (33)]
C...
      IF(IEVC.EQ.30) THEN
         UAVMAG=SQRT(UU*UU+VV*VV)
         EVBASE=0.025*H*UAVMAG/9.001D0
         DO J=1,NFEN
            EVTOT(J) = EVBASE
         END DO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.31) THEN
         UAVMAG=SQRT(UU*UU+VV*VV)
         EVBASE = EVCON*H*UAVMAG
         DO J=1,NFEN
            EVTOT(J) = EVBASE
         END DO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.32) THEN
         UAVMAG=SQRT(UU*UU+VV*VV)
         EVBASE = 0.025*H*UAVMAG/9.001D0
         EVTOT(1) = EVBASE/1000.D0
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.33) THEN
         UAVMAG=SQRT(UU*UU+VV*VV)
         EVBASE=EVCON*H*UAVMAG
         EVTOT(1)=EVBASE/1000.D0
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF

C...
C...  UAVG SQUARED FORMULATION FOLLOWING DAVIES (1990) + [EQ. (34)]
C...
      IF(IEVC.EQ.40) THEN
         UAVMAGS=UU*UU+VV*VV
         EVBASE=2.D0*UAVMAGS/9.001D0
         DO J=1,NFEN
            EVTOT(J) = EVBASE
         END DO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.41) THEN
         UAVMAGS=UU*UU+VV*VV
         EVBASE=EVCON*UAVMAGS
         DO J=1,NFEN
            EVTOT(J) = EVBASE
         END DO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.42) THEN
         UAVMAGS=UU*UU+VV*VV
         EVBASE=2.D0*H*UAVMAGS/9.001D0
         EVTOT(1)=EVBASE/1000.D0
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF
      IF(IEVC.EQ.43) THEN
         UAVMAGS=UU*UU+VV*VV
         EVBASE=EVCON*UAVMAGS
         EVTOT(1)=EVBASE/1000.D0
         SLOPE=(EVBASE-EVTOT(1))/(0.2D0*amb)
         BREAK=b+0.2D0*amb
         EVBEGIN=EVTOT(1)-SLOPE*b
         DO J=2,NFEN
            IF(SIGMA(J).LE.BREAK) THEN
               EVTOT(J)=EVBEGIN+SLOPE*SIGMA(J)
            ELSE
               EVTOT(J)=EVBASE
            ENDIF
         ENDDO
         GOTO 100
      ENDIF

C...
C...  MELLOR-YAMADA LEVEL 2.5 CLOSURE
C...
      IF((IEVC.EQ.50).or.(IEVC.EQ.51)) THEN
         CALL TURB(NODE,H,BSX,BSY,WSX,WSY,IStart,Z0B1)

         GOTO 100
      ENDIF

C...
C...  ONCE EDDY VISCOSITY IS COMPUTED
C...


 100  CONTINUE

C...
C...  CHECK SO THAT EDDY VISCOSITY NEVER GETS BELOW MINIMUM VALUE
C...
      DO J=1,NFEN
         IF(EVTOT(J).LT.EVMIN) THEN
            EVTOT(J)=EVMIN
         ENDIF
      ENDDO
C...
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C***********************************************************************
      END SUBROUTINE EDDYVIS
C***********************************************************************



C***********************************************************************
C                                                                      *
C   Solver for a vector U of length nfen for a tridiagonal system of   *
C          equations with the form                                     *
C                                                                      *
C    **                               **   **    **   **    **         *
C    * Bn   An                         *   *  Un  *   *  Rn  *         *
C    *                                 *   *      *   *      *         *
C    * Cn-1 Bn-1 An-1                  *   * Un-1 *   * Rn-1 *         *
C    *                                 *   *      *   *      *         *
C    *      Cn-2 Bn-2 An-2             *   * Un-2 *   * Rn-2 *         *
C    *                                 *   *      *   *      *         *
C    *             .......             *   * .... *   * .... *         *
C    *                                 *   *      *   *      *         *
C    *                  C3 B3 A3       *   *  U3  * = *  R3  *         *
C    *                                 *   *      *   *      *         *
C    *                     C2 B2 A2    *   *  U2  *   *  R2  *         *
C    *                                 *   *      *   *      *         *
C    *                        C1 B1    *   *  U1  *   *  R1  *         *
C    **                               **   **    **   **    **         *
C                                                                      *
C      A, B, C, U, R are adjustable size arrays                        *
C                                                                      *
C      U(1 - nfen) are the complex velocities from bottom to top       *
C                                                                      *
C                                                                      *
C                         R.L.  11/18/94                               *
C***********************************************************************
C
C  kmd48.33bc changed to use Implicit None
      SUBROUTINE TRIDIAG(A,B,C,R,U,nfen)
      use mod_terminate, only: terminate
      use sizes, only: myproc
      USE mod_logging, ONLY : ScreenUnit, setMessageSource,
     &   unsetMessageSource, logMessage, allMessage, DEBUG, ECHO, INFO
      IMPLICIT NONE
      INTEGER :: J,nfen

      COMPLEX(8) :: A(nfen), B(nfen), C(nfen), R(nfen), U(nfen)
      COMPLEX(8) :: P1
C
      call setMessageSource("tridiag")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      DO J=nfen-1,1,-1
         IF(ABS(B(J+1)).EQ.0.) then
            write(screenunit,*) 'Diagonal term in the VS matrix is zero'
            write(screenunit,*) '*********** Fatal error *************'
            write(1,*) 'Diagonal term in the VS matrix is zero'
            write(1,*) '*********** Fatal error *************'
            CALL terminate(myproc)
         endif
         P1=C(J)/B(J+1)
         B(J)=B(J)-A(J+1)*P1
         R(J)=R(J)-R(J+1)*P1
      END DO

      IF(ABS(B(1)).EQ.0.) then
         write(screenunit,*) 'B1 term in the VS matrix is zero'
         write(screenunit,*) '*********** Fatal error *************'
         write(1,*) 'B1 term in the VS matrix is zero'
         write(1,*) '*********** Fatal error *************'
         CALL terminate(myproc)
      endif

      U(1)=R(1)/B(1)

      DO J=2,nfen
         U(J)=(R(J)-A(J)*U(J-1))/B(J)
      END DO

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C***********************************************************************
      END SUBROUTINE TRIDIAG
C***********************************************************************


C***********************************************************************
C                                                                      *
C                    VSDISP.FOR - VERSION                              *
C                                                                      *
C  This subroutine computes the dispersion terms: Duu, Duv, Dvv        *
C  for the FE VS method as derived in "Dispersion Terms 1/27/92"       *
C                                                                      *
C                                                                      *
C                          R.L. 05/25/00                               *
C                          R.L. 06/22/05                               *
C***********************************************************************

      Subroutine VSDISP(it,node,H,UU,VV,Duu,Duv,Dvv)
      USE mod_logging, ONLY : ScreenUnit, allMessage, DEBUG, scratchMessage,
     &    setMessageSource, unsetMessageSource, allMessage, logMessage,
     &    DEBUG, ECHO, INFO

      USE GLOBAL_3DVS, ONLY : Qkp1,NFEN,SIGMA,GAMMA,AMB,IY

      IMPLICIT NONE

      REAL(8) :: UU,VV,Duu,Duv,Dvv
      REAL(8) :: delst,D1P1,Txx,Txy,Tyy,Ump,Vmp,Um,Vm,H
C     jgf46.00 explicitly declared the following variables
      INTEGER it
      INTEGER node
Cjgfdebug      INTEGER ierr
      INTEGER j
Cjgfdebug      INTEGER k
      INTEGER m
Cjgfdebug      INTEGER n
      INTEGER mp
Cjgfdebug      INTEGER istart
Cjgfdebug      INTEGER il
Cjgfdebug      INTEGER ibc
Cjgfdebug      REAL(8) rkap
Cjgfdebug      REAL(8) B1
Cjgfdebug      REAL(8) B123
Cjgfdebug      REAL(8) g2
Cjgfdebug      REAL(8) g3
Cjgfdebug      REAL(8) g4
Cjgfdebug      REAL(8) g5
Cjgfdebug      REAL(8) g6
Cjgfdebug      REAL(8) E1
c
c     parameters
c
Cjgfdebug      ierr = 0
c
c     don't waste time if entire profile is zero
c
      call setMessageSource("vsdisp")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      do j=1,nfen
         gamma(j)=Qkp1(node,j)
      end do

      do j=1,nfen
         if(gamma(j).ne.(0.d0,0.d0)) goto 1
      end do
      Duu = 0.d0
      Duv = 0.d0
      Dvv = 0.d0
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      return

c
c     Domain is split up into intervals according to sigma grid and
c     contribution to total integral is computed on each interval.
c

 1    Txx = 0.d0
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
 100  continue

      if(Txx.lt.0.) then
         write(screenunit,1001)
 1001    format(/'**** Serious Error Detected in SUBROUTINE VSDISP ***'/
     +        '          the partial uu dispersion term < 0        '/
     +        '        Diagnostic information written to unit 1    ')
         write(screenunit,1101) node,it
         write(1,1002)
 1002    format(/'**** Serious Error Detected in SUBROUTINE VSDISP ***'/
     +        '          the partial uu dispersion term < 0        ')
         write(1,1101) node,it
         write(1,1021)
         do m=1,NFEN            !ntotn -> NFEN
            write(1,1022) sigMA(m), gamma(m) !sigtot -> sigMA
         enddo
         write(1,2101) Txx
 2101    format('  Par uu = ',e14.6,' It will be set to 0.')
         Txx = 0.d0
      endif

      if(Tyy.lt.0.) then
         write(screenunit,1003)
 1003    format(/'**** Serious Error Detected in SUBROUTINE VSSOL ****'/
     +        '          the partial vv dispersion term < 0        '/
     +        '       Diagnostic information written to unit 1    ')
         write(screenunit,1101) node,it
         write(1,1004)
 1004    format(/'**** Serious Error Detected in SUBROUTINE VSSOL ****'/
     +        '          the partial vv dispersion term < 0        ')
         write(1,1101) node,it
         write(1,1021)
         do m=1,NFEN            !ntotn -> NFEN
            write(1,1022) sigMA(m),gamma(m) !sigtot -> sigMA
         enddo
         write(1,2102) Tyy
 2102    format('  Par vv = ',e14.6,' It will be set to 0.')
         Tyy = 0.d0
      endif

 1101 format(/' Error occurred at node ',I6,' time step ',I8/)
 1021 format(8x,'sigMA',11x,'  u ',12x,'  v ')
 1022 format(1x,3e16.8)

      Duu = Txx*H/amb - H*UU*UU
      Duv = Txy*H/amb - H*UU*VV
      Dvv = Tyy*H/amb - H*VV*VV

      if(Duu.lt.0.) Duu = 0.d0
      if(Dvv.lt.0.) Dvv = 0.d0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      return
c***********************************************************************
      end subroutine vsdisp
c***********************************************************************



c***********************************************************************
c     Subroutine to compute the Inm integral                           *
c                                                                      *
c     Note, Inm is based only on the f.e. grid and therefore is the    *
c     same for all horizontal nodes.                                   *
c                                                                      *
c     11/24/01                                                         *
c                                                                      *
c***********************************************************************

      subroutine InmINT()
Cjgfdebug      USE GLOBAL, ONLY : ScreenUnit
      USE GLOBAL_3DVS, ONLY : Sigma, NFEN, INM
      USE mod_logging, only: setMessageSource,
     &    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
      IMPLICIT NONE
      INTEGER k
C
      call setMessageSource("InmINT")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C     integral over lower element of psi(k-1)*psi(k)
      Inm(1,1) = 0.d0
C     integral over upper element of psi(k+1)*psi(k)
      Inm(1,3) = (Sigma(2)-Sigma(1))/6.d0
C     integral over both elements of psi(k)*psi(k)
      Inm(1,2) = 2.d0*Inm(1,3)

      do k=2,NFEN-1
C     integral over lower element of psi(k-1)*psi(k)
         Inm(k,1) = Inm(k-1,3)
C     integral over upper element of psi(k+1)*psi(k)
         Inm(k,3) = (Sigma(k+1)-Sigma(k))/6.d0
C     integral over both elements of psi(k)*psi(k)
         Inm(k,2) = 2.d0*(Inm(k,1)+Inm(k,3))
      enddo

C     integral over lower element of psi(k-1)*psi(k)
      Inm(NFEN,1) = Inm(NFEN-1,3)
C     integral over both elements of psi(k)*psi(k)
      Inm(NFEN,2) = 2.d0*Inm(NFEN,1)
C     integral over upper element of psi(k+1)*psi(k)
      Inm(NFEN,3) = 0.d0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
c***********************************************************************
      end subroutine InmINT
c***********************************************************************


c***********************************************************************
c      Subroutine to compute the LVn integral (used only in VS)        *
c                                                                      *
c                            11/14/94                                  *
c***********************************************************************

      subroutine LVnInt()
Cjgfdebug      USE GLOBAL, ONLY : ScreenUnit
      USE GLOBAL_3DVS, ONLY: Sigma, NFEN, LVN
      use mod_logging, only: setMessageSource,
     &    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
      IMPLICIT NONE
      INTEGER n
C
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
c***********************************************************************
      END subroutine LVnInt
c***********************************************************************

c***********************************************************************
c      Subroutine to compute the KQnm integral                         *
c                                                                      *
c                            01/26/00                                  *
c***********************************************************************

      subroutine KQnmInt(KQnm,Kq)
Cjgfdebug      USE GLOBAL, ONLY : ScreenUnit
      USE GLOBAL_3DVS, ONLY : SIGMA,NFEN,MNFEN
      use mod_logging, only: setMessageSource,
     &    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
      IMPLICIT NONE
      REAL(8) :: KQnm(MNFEN,3),Kq(MNFEN)
      REAL(8) :: EM,EP
      INTEGER n

      call setMessageSource("KQnmInt")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif



c     do 1st element by hand

      EP=(Kq(2)+Kq(1))/(sigma(2)-sigma(1))/2.d0
      KQnm(1,1)=0.D0
      KQnm(1,2)=EP
      KQnm(1,3)=-EP

c     loop through interior elements
c     NOTE: the integrals from sigma(n-1) to sigma(n) are simply integrals
c     the integrals from sigma(n) to sigma(n+1) from previous element

      DO n=2,nfen-1
         EM=EP
         EP=(Kq(n+1)+Kq(n))/(sigma(n+1)-sigma(n))/2.d0
         KQnm(n,1)=-EM
         KQnm(n,2)=EM+EP
         KQnm(n,3)=-EP
      ENDDO

c     do last element by hand

      EM=EP
      KQnm(nfen,1)=-EM
      KQnm(nfen,2)=EM
      KQnm(nfen,3)=0.D0

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
c************************************************************************
      end subroutine KQnmInt
c************************************************************************

c************************************************************************
c MY2.5 TURBULENCE MODEL PROGRAM
c VELOCITY SOLUTION VERSION
c
c written by R. Luettich based on earlier subroutines by R. Grenier
c
c OVERVIEW
c
c     This code uses the quasi-equilibrium version of the Mellor-Yamada
c     turbulence scheme (Mellor and Yamada, 1982, Blumberg and Mellor,
c     1987 and Galperin et al., 1988) to solve transport equations for
c     q**2 and q**2l.  The parameters q and l are used to compute the
c     eddy viscosity according to the relation:
c
c             Km = Smql
c
c     where Km is the eddy viscosity and Sm is a stability parameter.
c
c     This routine is called from the EDDYVIS subroutine during the
c     internal mode solution for ievc = 50.  There are separate DSS and
c     VS versions of the code.
c
c    *** SEE SECTION TITLED "USERS GUIDE" BELOW ***
c
c     This routine should be linked with the 3D Code subroutines (ADCIRC
c     + VS) after running the setup program
c
c PARAMETER DEFINITIONS:
c
c argument list:
c
c nh          : horizontal node counter
c H           : total depth at current time step (includes finite amplitude)
c it          : current time step value
c delt        : simulation incremental time step (uses same step as
c               internal and external modes)
c nws         : switch for wind stress (nws=1 for wind on,
c               nws=0 for wind off). This affects the surface
c               boundary condition. Set in the external mode and passed in.
c bsx,bsy     : x,y bottom stresses
c wsx,wsy     : x,y surface stresses
c Z0B         : bottom roughness/mixing length (see USER's notes
c               below)
c Z0S         : surface roughness/mixing length (see USER's notes
c               below)
c Z0Sw        : surface roughness/mixing length (see USER's notes
c               below) if Z0S is greater than or equal to water depth
c
c Coeff block:
c
c Sq          : stability function used in definition of the eddy
c               diffusivity of the turbulence parameters (Sq = 0.2)
c Sm          : stability function used in definition of the eddy
c               viscosity (see Galperin, et al. for exact form)
c Sh          : stability function used in definition of the eddy
c               diffusivity of the density (see Galperin, et al.
c               for exact form)
c B1,E1,E2,E3 : empirical constants (values given below; see Mellor and
c               Yamada, 1982 & Blumberg et al, 1992 for discussion)
c
c Turbmod block:
c
c Q (real, imag): x,y direction dependent variables from the vertical
c               solution
c
c Vgrid block, Setup Block
c
c parameters as defined in the internal mode subroutines
c
c Other calculation parameters:
c
c q2l         : 2xTKE times the master length scale
c q20         : previous time step value of q2
c l           : master length scale
c Km          : momentum eddy viscosity
c Kq          : turbulence eddy viscosity
c Kh          :
c Mqa,Mqb,Mqc : LHS matrix diagonals for the q2 and q2l solutions
c LVq         : RHS load vector for the q2 and q2l solutions
c q2          : solution vector returned from the tridiagonal solver
c w           : wall function (see Blumberg at al., 1992)
c BVflux2     : Brunt Vaislai frequency squared
c Gh          : dimensionless density function
c SIGT        : profile of density (sigma T)
c
c Flags, run control and miscellaneous parameters
c
c iden        : density flag (iden <> 0 density included,iden = 0 no density)
c im          : run type flag - vs (im = 1) or dss (im = 2)
c il          : length scale flag - il = 2 or 3 for algebraic length scale (see below),
c               il = 1 for length scale computed from the q2l equation (see user notes below)
c ibc         : boundary condition flag for q2 (0=constant/zero,1=no-flux)
c               User should set ibc = 0 for nws = 1
c
c
c SUBROUTINES
c
c turb        : main module - handles input, run control and output.
c ADC_TRIDAG  : tridiagonal matrix solver
c
c
c USER'S GUIDE
c
c     The user must set a number of parameter statements and flags prior
c     to operation:
c
c parameter statements **Check all subroutines**
c  MNFEN = maximum number of vertical nodes
c  mnp = maximum number of horizontal nodes
c
c flags (SEE DISCUSSION ABOVE UNDER flags, run control and miscellaneous parameters)
c
c     The value of Z0B should be set consistent with the bottom boundary
c     condition used in the internal mode solution.  For a no-slip BBC,
c     choose Z0B as a physical roughness height, e.g. 0.005m.  For a
c     slip case, set Z0B to a value consistent with the slip coefficient
c     used via the log profile, (e.g., 1 m).  Tests suggest that when a
c     no-slip condition is used with the VS model (for which the bottom
c     nodes are very tightly spaced) it is best to set the time
c     weighting parameter for the momentum diffusion term (alpha3 in
c     input unit 15) to 1.0 to avoid instability problems. When a slip
c     condition is used (or whenever the bottom grid spacing is about
c     1m) a Crank-Nicholson approach (alpha3 = 0.5) is acceptable.
c
c     The number of vertical nodes used in the solution of the
c     turbulence equations is the same as the number of nodes used in
c     the solution of the dependent variable (velocity or stress) in the
c     internal mode.  This is unlike the other forms of eddy viscosity,
c     for which the two grids are different and the number of nodes used
c     to define the eddy viscosity is generally less.
c
c     This code assumes that any density field is passed into the
c     routine via a common block called "DENSITY3D", which includes a
c     density profile on the internal mode solution grid at each point
c     in the horizontal.  Any updating of this profile must be done
c     externally and passed into this routine.  The density be passed as
c     sigma t units, and the background density is RHOWAT0.
c
c     Model output is limited to printing results for a single
c     horizontal node.  Additional coding would be required to create
c     full output files.
c
c
c REFERENCES
c
c     Blumberg, A.F. and G.L. Mellor, A Description of a
c     Three-Dimensional Coastal Ocean Circulation Model, In:
c     Three-Dimensional Coastal Ocean Models, edited by N.S. Heaps,
c     pp. 1-16, American Geophysical Union, Washington, D.C., 1987.
c
c     Blumberg, A.F., B. Galperin and D.J. O'Connor, Modeling vertical
c     structure of open channel flows, Journal of Hydraulic Engineering,
c     118, 1119-1134., 1992.
c
c     Galperin, B., L.H. Kantha, S. Hassid and A. Rosati, A
c     quasi-equilibrium turbulent energy model for geophysical flows,
c     Journal of the Atmospheric Sciences, 45, 55-62, 1988.
c
c     Mellor, G.L. and T. Yamada, Development of a turbulence closure
c     model for geophysical fluid problems, Reviews of Geophysics and
c     Space Physics, 20, 851-875, 1982.
c
c************************************************************************
      subroutine turb(nh,H,bsx_loc,bsy_loc,wsx,wsy,istart,Z0B1)
      USE SIZES
Cjgfdebug      USE GLOBAL, ONLY : ScreenUnit
Cjgfdebug      USE GLOBAL_3DVS, EXCEPT_BSX => BSX ,EXCEPT_BSY => BSY
      USE GLOBAL_3DVS
      USE ADC_CONSTANTS, ONLY: G
Crjwdebug     USE NodalAttributes, ONLY : Z0B_var
      IMPLICIT NONE

      COMPLEX(8) :: dQdz,dQdz1,dQdz2
      REAL(8) :: H,H_crit

      REAL(8) :: bsx_loc,bsy_loc,WSX,WSY
      REAL(8),SAVE :: H2
C  kmd48.33bc need to change Sh to an array for transport
      REAL(8),SAVE,ALLOCATABLE :: Sh(:)
      REAL(8) :: Z0Sw
      REAL(8) :: Z0B1
      REAL(8),SAVE,ALLOCATABLE :: KQnm(:,:)
      REAL(8),SAVE,ALLOCATABLE :: Mqa(:),Mqb(:),Mqc(:)
      REAL(8),SAVE,ALLOCATABLE :: LVq(:),Sm(:)
      REAL(8),SAVE,ALLOCATABLE :: q2(:),q2prev(:)
      REAL(8),SAVE,ALLOCATABLE :: q2l(:),q2lprev(:)
      REAL(8),SAVE,ALLOCATABLE :: wall(:),rmlen(:),rmlen2(:)
      REAL(8),SAVE,ALLOCATABLE :: BVfreq2(:),spgrad2(:)
      REAL(8),SAVE,ALLOCATABLE :: Kq(:),Km(:),Kh(:)
      REAL(8),SAVE,ALLOCATABLE :: prod(:),diss(:)
C     jgf46.00 explicitly declared the following variables
C     jgf46.12 gave them SAVE and PARAMETER properties (from jgf45.18)
      INTEGER nh, istart, n
      INTEGER,SAVE :: il, ibc
      REAL(8), PARAMETER :: rkap = 0.41d0
      REAL(8), PARAMETER :: B1 = 16.6d0
      REAL(8), PARAMETER :: g2 = 0.39327d0
      REAL(8), PARAMETER :: g3 = 3.0858d0
      REAL(8), PARAMETER :: g4 = 34.676d0
      REAL(8), PARAMETER :: g5 = 6.1272d0
      REAL(8), PARAMETER :: g6 = 0.49393d0
      REAL(8), PARAMETER :: E1 = 1.8d0
      REAL(8), PARAMETER :: E2 = 1.33d0
      REAL(8), PARAMETER :: E3 = 0.25d0
      REAL(8), PARAMETER :: q2min=1.d-8
      REAL(8), SAVE :: B123, Sq, sig, rl1, rl2, HOamb, zval, HOamb2
      REAL(8), SAVE :: dsig, drhodz, dudz, dvdz, dsig1, dsig2
      REAL(8), SAVE :: drhodz1, drhodz2, BSlay, SSlay, db, db2, ds, ds2
      REAL(8), SAVE :: qprev, qlprev, Gh, tdiss, elmax
      REAL(8), SAVE :: coef1, coef2, coef3, coef4, coef5
      REAL(8) :: USTARS

      call setMessageSource("turb")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))

c     allocate local arrays
      if(.not. turb_allocated) then
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
         turb_allocated = .true.
      endif


          Z0Sw=Z0S
          IF (IEVC.EQ.51) THEN
           Z0Sw=EVCON*USTARS*USTARS/g
          ELSE
            Z0Sw=Z0S
          ENDIF
C  RJW taken from Jones monismith 2008
C  modify surface roughness length to account for wave breaking
C         Z0Sw = 1400.0d0*USTARS*USTARS/g   ! ( Bye (1988)  , Craig and Banner (1994))
C         Z0Sw = 14000.0d0*USTARS*USTARS/g  ! (Stipps et al. 2005 extended to wave breaking)
C         Z0Sw = 200000.0d0*USTARS*USTARS/g  ! ( O(10^5)  Stacy et al. (1999) for wave ages > 25)
C         Z0Sw = 1.3*Wave Height            ! (sig wave height) (Jones and Monosmith (2007) k-omega)
C         Z0Sw = (0.25-1.0)*Wave Height     ! (sig wave height) (Burchard (2001) k-epsilon )

       if (Z0Sw .gt. H-Z0B1) then
          Z0Sw = H-Z0B1 ! bound Z0S to water depth
       else
         Z0Sw=Z0Sw
       endif
       if (Z0Sw .lt. Z0S ) Z0Sw = Z0S


c
c     At time step 1, set flags and initialize variables
c
C.  RJW moved up and outside of if statement to allow variable computation in case of hotstart
         il = 1  ! length scale flag, = 1 for length scale from q2l eqn.
c         il = 2  ! length scale flag, = 2 for length scale from algebraic eq. (see below)
c         il = 3  ! length scale flag, = 3 for length scale from algebraic eq. (see below)
        ibc = 1
        !ibc = 0 ! surface b.c. flag, = 0 for specified stress, = 1 for no flux
        !ibc = 1 ! surface b.c. flag, = 0 for specified stress, = 1 for no flux
        !if((nws.eq.0).or.(nws.eq.100)) ibc =1
c     set constants and other parameters
        B123 = B1**(2.d0/3.d0)
        Sq = 0.2d0 ! initialize the stability constant stability function
c     initialize the Brunt-Vaisala freq squared = 0 if density not considered
         if(iden.eq.0) then
            do n=1,nfen
               BVfreq2(n)=0.d0
            end do
         endif

      if(istart.eq.1) then
c     initialization for a cold start only
         if(ihot.eq.0) then
            do n = 1,nfen
               q20(nh,n) = q2min ! initialize q2 to a minimal value
c     initilize l to a minimal value if computed from q2l equation
               if(il.eq.1) then
C           l(nh,n)=rkap*Z0B1                        !rog way
          l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
               endif
c     set l to an exponential type length scale (davies and xing)
               if(il.eq.2)then
                  sig = (sigma(n)+1.d0)/amb
                 rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
                 rl2 = 1.d0/(rkap*(H-sig*H+Z0Sw))
                  l(nh,n) = 1.d0/(rl1+rl2)
               endif
c     set l to a linear variation with kz over lower 15% with constant above
               if(il.eq.3)then
                  HOamb=H/amb
                  zval = (sigma(n)+1.d0)*HOamb-H
              if(sigma(n).le.-0.7d0) l(nh,n)=rkap*(H+zval+Z0B1)
              if(sigma(n).gt.-0.7d0) l(nh,n)=rkap*(0.15d0*H+Z0B1)
               endif
            end do
         endif
c     end cold start initialization
      endif                     !end of 1st time step section
c
c     Begin calculations for each time step
c
C.RJW..
C identify if node was dry last timestep but is now wet
C the initialized value of L can be recomputed to a more meaningful
C value before used to compute q2

      if (L(nh,1) .eq. -9999d0 ) Then
         do n = 1,nfen
         l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
         enddo
       endif
C now check q20 it may also be zero
      if (q20(nh,1) .eq. -9999d0 ) Then
            do n = 1,nfen
               q20(nh,n) = q2min
            end do
       Endif

c         if(il.eq.2)then
c                  sig = (sigma(n)+1.d0)/amb
c                  rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
c                  rl2 = 1.d0/(rkap*(H-sig*H+Z0S))
c                  l(nh,n) = 1.d0/(rl1+rl2)
c               endif

Cjgfdebug      H2 = H*H
      HOamb=H/amb
      HOamb2=(H/amb)*(H/amb)

c     Compute the speed gradient squared, density gradient, BV freq
c     and split out the mixing length

      rmlen(1)=l(nh,1)
      rmlen2(1)=rmlen(1)*rmlen(1)
      dsig=sigma(2)-sigma(1)
      if(iden.ne.0) then
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
         if(iden.ne.0) then
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
      if(iden.ne.0) then
         drhodz=((SIGT(nh,nfen)-SIGT(nh,nfen-1))/dsig)/HOamb
         BVfreq2(nfen)=-GORHO*drhodz
      endif
      dQdz=((Q(nh,nfen)-Q(nh,nfen-1))/dsig)/HOamb
      dudz=real(dQdz)
      dvdz=aimag(dQdz)
      spgrad2(nfen)=dudz*dudz+dvdz*dvdz

c     Compute the wall function if the mixing length is determined from
c     q2l eqn

      if(il.eq.1)then
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

c     Compute the stability functions, eddy viscosity and partial
c     turbulence production & dissipation terms using information from
c     the previous time step
C RJW rweaver add wave production of TKE
C alpha * USTARS*USTARS*USTARS, USTARS=SQRT(SQRT(WSX*WSX+WSY*WSY))
C and alpha = 60 - 250

      do n = 1,nfen
         q2prev(n)=q20(nh,n)
         q2lprev(n)=q2prev(n)*rmlen(n)
         qprev=sqrt(q2prev(n))
         qlprev=qprev*rmlen(n)
         Gh=-BVfreq2(n)*rmlen2(n)/q2prev(n)
         if(Gh.gt.0.0233) Gh=0.0233
         Sm(n)=(g2-g3*Gh)/((1.d0-g4*Gh)*(1.d0-g5*Gh))
         Sh(n)=g6/(1.d0-g4*Gh) ! kmd48.33bc changed Sh to an array
         Km(n)=Sm(n)*qlprev
         if(Km(n).lt.EVMIN) Km(n)=EVMIN
c     if(istart.eq.1) Km(n)=EVMIN                        !rog way
         Kq(n)=Sq*qlprev
         if(Kq(n).lt.EVMIN) Kq(n)=EVMIN
c     Kq(n)=Km(n)*Sq/Sm(n)                           !rog way
         Kh(n)=Sh(n)*qlprev ! kmd48.33bc changed Sh to an array
         if(Kh(n).lt.EVMIN) Kh(n)=EVMIN
         prod(n)=Km(n)*spgrad2(n)-Kh(n)*BVfreq2(n)
         IF(B1*rmlen(n).EQ.0.d0) PRINT*, B1,rmlen(n)
         diss(n)=qprev/(B1*rmlen(n))
c     diss(n)=q2prev(n)*Sm(n)/(Km(n)*B1)             !rog way
           if ((n.eq.nfen).and.(IEVC.eq.51)) then
             prod(n)=prod(n)+(150.0*USTARS*USTARS*USTARS)
           endif
      enddo

c     Compute the q2 LHS Matrix and RHS Load Vector

      call KQnmInt(KQnm,Kq)

      Mqa(1) = 0.d0
      Mqb(1) = 1.d0
      Mqc(1) = 0.d0
      LVq(1) = B123*sqrt(bsx_loc*bsx_loc+bsy_loc*bsy_loc)

      coef2 = delt*theta1/HOamb2
c     coef2 = theta1/HOamb2                !rog way
      coef4 = 2.d0*delt
c     coef4 = 2.d0                         !rog way
      coef5 = delt*(1.d0-theta1)/HOamb2
c     coef5 = (1.d0-theta1)/HOamb2         !rog way
      do n=2,nfen-1
         tdiss = 2.d0*diss(n)
         coef1 = 1.d0 + delt*theta2*tdiss
c     coef1 = 1.d0/delt + theta2*tdiss       !rog way
         Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
c     Mqa(n) = KQnm(n,1)*coef2               !lumping
         Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
c     Mqb(n) = (Inm(n,1)+Inm(n,2)+Inm(n,3))*coef1+KQnm(n,2)*coef2  !lumping
         Mqc(n) = Inm(n,3)*coef1+KQnm(n,3)*coef2
c     Mqc(n) = KQnm(n,3)*coef2                                     !lumping

         coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
c     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss                      !rog way
         LVq(n) = Inm(n,1)*(coef3*q2prev(n-1)+coef4*prod(n-1))
     &        -KQnm(n,1)*coef5*q2prev(n-1)
     &        + Inm(n,2)*(coef3*q2prev(n  )+coef4*prod(n  ))
     &        -KQnm(n,2)*coef5*q2prev(n)
     &        + Inm(n,3)*(coef3*q2prev(n+1)+coef4*prod(n+1))
     &        -KQnm(n,3)*coef5*q2prev(n+1)
c     LVq(n) =                 -KQnm(n,1)*coef5*q2prev(n-1)     !lumping
c     &         + (Inm(n,1)+Inm(n,2)+Inm(n,3))                  !lumping
c     &         *(coef3*q2prev(n)+coef4*prod(n))     !lumping
c     &         -KQnm(n,2)*coef5*q2prev(n)       !lumping
c     &         -KQnm(n,3)*coef5*q2prev(n+1)     !lumping

      enddo

      if(ibc.eq.0) then
         Mqa(nfen) = 0.d0
         Mqb(nfen) = 1.d0
         Mqc(nfen) = 0.d0
         LVq(nfen) = B123*sqrt(wsx*wsx+wsy*wsy)
      endif

      if(ibc.eq.1)then
         n=nfen
         tdiss = 2.d0*diss(n)
         coef1 = 1.d0 + delt*theta2*tdiss
c     coef1 = 1.d0/delt + theta2*tdiss                   !rog way
         Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
c     Mqa(n) = KQnm(n,1)*coef2                           !lumping
         Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
c     Mqb(n) = (Inm(n,1)+Inm(n,2))*coef1+KQnm(n,2)*coef2 !lumping
         Mqc(n) = 0.d0

         coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
c     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss            !rog way
         LVq(n) = Inm(n,1)*(coef3*q2prev(n-1)+coef4*prod(n-1))
     &        -KQnm(n,1)*coef5*q2prev(n-1)
     &        + Inm(n,2)*(coef3*q2prev(n  )+coef4*prod(n  ))
     &        -KQnm(n,2)*coef5*q2prev(n)
c     LVq(n) =    -KQnm(n,1)*coef5*q2prev(n-1)   !lumping
c     &    + (Inm(n,1)+Inm(n,2))*(coef3*q2prev(n)+coef4*prod(n))       !lumping
c     &         -KQnm(n,2)*coef5*q2prev(n)     !lumping
      endif

c     Solve the system for q2

      CALL ADC_TRIDAG(Mqa,Mqb,Mqc,LVq,q2,nfen)

c     Transfer to global array and check for zero or negative values
c     (generally for startup)

      do n = 1,nfen
         if(q2(n).le.0.)then
            q2(n) = q2min
         endif
         q20(nh,n) = q2(n)
      enddo

c     Compute length scale

      if(il.eq.2)then
         do n=1,nfen
            sig = (sigma(n)+1.d0)/amb
            rl1 = 1.d0/(rkap*(sig*H+Z0B1)*exp(-amb*sig))
            rl2 = 1.d0/(rkap*(H-sig*H+Z0Sw))
            l(nh,n) = 1.d0/(rl1+rl2)

            if(l(nh,n).lt.0.001) l(nh,n)=0.001
            if(l(nh,n).gt.H) l(nh,n)=H
         enddo
      endif

      if(il.eq.3)then
         do n=1,nfen
            HOamb=H/amb
            zval = (sigma(n)+1.d0)*HOamb-H
            if(sigma(n).le.-0.7d0) l(nh,n)=rkap*(H+zval+Z0B1)
            if(sigma(n).gt.-0.7d0) l(nh,n)=rkap*(0.15d0*H+Z0B1)
         enddo
      endif

      if(il.eq.1)then

c     Compute the q2l LHS Matrix and RHS Load Vector

         Mqa(1) = 0.d0
         Mqb(1) = 1.d0
         Mqc(1) = 0.d0
         LVq(1) = rkap*BSlay*q2(1)

         coef2 = delt*theta1/HOamb2
c     coef2 = theta1/HOamb2                                     !rog way
         coef5 = delt*(1.d0-theta1)/HOamb2
c     coef5 = (1.d0-theta1)/HOamb2                              !rog way
         do n=2,nfen-1
            tdiss = wall(n)*diss(n)
            coef1 = 1.d0 + delt*theta2*tdiss
c     coef1 = 1.d0/delt + theta2*tdiss                        !rog way
            Mqa(n) = Inm(n,1)*coef1+KQnm(n,1)*coef2
c     Mqa(n) = KQnm(n,1)*coef2                                !lumping
            Mqb(n) = Inm(n,2)*coef1+KQnm(n,2)*coef2
c     Mqb(n) = (Inm(n,1)+Inm(n,2)+Inm(n,3))*coef1+KQnm(n,2)*coef2   !lumping
            Mqc(n) = Inm(n,3)*coef1+KQnm(n,3)*coef2
c     Mqc(n) = KQnm(n,3)*coef2                                !lumping

          coef3 = 1.d0 - delt*(1.d0-theta2)*tdiss
c     coef3 = 1.d0/delt - (1.d0-theta2)*tdiss                 !rog way
            coef4 = E1*rmlen(n)*delt
c     coef4 = E1*rmlen(n)                                     !rog way
            LVq(n) = Inm(n,1)*(coef3*q2lprev(n-1)+coef4*prod(n-1))
     &           -KQnm(n,1)*coef5*q2lprev(n-1)
     &           + Inm(n,2)*(coef3*q2lprev(n  )+coef4*prod(n  ))
     &           -KQnm(n,2)*coef5*q2lprev(n)
     &           + Inm(n,3)*(coef3*q2lprev(n+1)+coef4*prod(n+1))
     &           -KQnm(n,3)*coef5*q2lprev(n+1)
c     LVq(n) =      -KQnm(n,1)*coef5*q2lprev(n-1)  !lumping
c     &    + (Inm(n,1)+Inm(n,2)+Inm(n,3))                           !lumping
c     &     *(coef3*q2lprev(n)+coef4*prod(n))                       !lumping
c     &       -KQnm(n,2)*coef5*q2lprev(n)    !lumping
c     &     -KQnm(n,3)*coef5*q2lprev(n+1)  !lumping
c
         enddo

         Mqa(nfen) = 0.d0
         Mqb(nfen) = 1.d0
         Mqc(nfen) = 0.d0
         LVq(nfen) = SSlay*rkap*q2(nfen)

c     Solve the system for q2l

         CALL ADC_TRIDAG(Mqa,Mqb,Mqc,LVq,q2l,nfen)

c     Transfer to global array and check for stability limit

         l(nh,1) = q2l(1)/q2(1)
         do n = 2,nfen-1
            l(nh,n) = q2l(n)/q2(n)
            if(l(nh,n).lt.0.)
     &     l(nh,n)=rkap*(Z0B1*(a-sigma(n))-Z0Sw*(b-sigma(n)))/amb
            if(l(nh,n).gt.H) l(nh,n)=H
            if((iden.ne.0).and.(BVfreq2(n).gt.0.0)) then
               elmax = 0.53D0*sqrt(q2(n))/sqrt(BVfreq2(n))
               if(l(nh,n).gt.elmax) l(nh,n)=elmax
            endif
         end do
         l(nh,nfen) = q2l(nfen)/q2(nfen)

C RJW 02252009 ...
c Now we do a little bookeeping, to be sure things stay in
c reasonable bounds.
c 1) set a lower limit on the turbulent length scale of 10 cm
c 2)in shallow water we want to modify the calculation of l
c set Hlimit = 1.0 and l_2 = H
c If the total water depth is less than 1.0 meter
c then we want to set the turbulent length scale to
c equal the water depth
        H_crit = 1.0

        do n = 1, nfen

          if (l(nh,n) .lt. 0.10 ) l(nh,n) = 0.10

          if ( H .le. H_crit ) then

           l(nh,n) = H/H_crit*l(nh,n) + (1-H/H_crit)*H

          endif

        enddo
c...

      endif

c
c     compute eddy viscosity and store variables for output and next step
c
C   kmd48.33bc need to add the vertical diffusion term for transport
C              the lower limit will be the same as that for the eddy
C              viscosity
      do n = 1,nfen
         EVTOT(n) = Sm(n)*sqrt(q2(n))*l(nh,n)
         NTVTOT(n) = Sh(n)*sqrt(q2(n))*l(nh,n)
         if(EVTOT(n).lt.EVMIN) EVTOT(n)=EVMIN
         if(EVTOT(n).gt.100.0d0) EVTOT(n)=100.0d0
         if(NTVTOT(n).lt.EVMIN) NTVTOT(n)=EVMIN
      end do

#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C***********************************************************************
      END subroutine turb
C***********************************************************************



C***********************************************************************
C Subroutine ADC_TRIDAG                                                *
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
      SUBROUTINE ADC_TRIDAG(A,B,C,R,U,N)
      USE mod_logging, ONLY : setMessageSource, ScreenUnit,
     &    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
      use mod_terminate, only: terminate
      use sizes, only: myproc
      IMPLICIT NONE
      INTEGER :: J, N

      REAL(8) :: A(N),B(N),C(N),R(N),U(N)
      REAL(8) :: BET,GAM(N)

      call setMessageSource("ADC_tridag")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
c
c     check for zero elements on diagonal
c
      DO J=1,N
         if(B(j).EQ.0.) then
            write(screenunit,*) 'Problem in Tridag Solver.  ',
     +           'B array value in row ',j,' = 0'
            CALL terminate(myproc)
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
            CALL terminate(myproc)
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
C***********************************************************************
      END SUBROUTINE ADC_TRIDAG
C***********************************************************************



C****************************************************************************************
C   Subroutine to interpolate baroclinic pressure (BCP) to a specified sigma value      *
C   (SigmaNN) given an initial guess of which sigma level is closest to the specified   *
C   value.                                                                              *
C                                                                                       *
C                                    R.L.  5/04/01                                      *
C                                    R.L.  5.19/03                                      *
C                                    R.L.  5/08/05                                      *
C****************************************************************************************
C
!      SUBROUTINE ZSURFBUOY(SigmaNN,BCPressNN,NN,J)
!      USE GLOBAL, ONLY : ScreenUnit
!      USE GLOBAL_3DVS
!      REAL(8) :: BCPressNN
!      REAL(8) :: SigmaNN     !Sigma value of a neighbor node
!
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
!
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
!
! 99   BCPressNN=(BCP(NN,LAbov)-BCP(NN,LBelo)) !interpolation
!     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + BCP(NN,LBelo)
!
! 100  CONTINUE
!
!      RETURN
!      END
!

C***************************************************************************
C   Subroutine to write out 3D Hot Start info
C
C                                    R.L.  8/16/05
C
C     jgf49.44: This subroutine is only called if a binary hotstart
C     file was specified.
C***************************************************************************
C
C  kmd48.33bc changed to follow the global IO format
      SUBROUTINE HSTART3D_OUT(IT)
      USE SIZES
      USE GLOBAL
      USE GLOBAL_3DVS
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL_IO
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : initNetCDFHotstart3D, writeNetCDFHotstart3D,
     &   writeNetCDFHotstart3DVar
#endif

      IMPLICIT NONE
      INTEGER, intent(in) :: IT ! current time step
      INTEGER NH, I, J ! loop counters
      INTEGER N,M

      type(OutputDataDescript_t),SAVE :: Duudescript
      type(OutputDataDescript_t),SAVE :: Duvdescript
      type(OutputDataDescript_t),SAVE :: Dvvdescript
      type(OutputDataDescript_t),SAVE :: Uudescript
      type(OutputDataDescript_t),SAVE :: Vvdescript
      type(OutputDataDescript_t),SAVE :: Bsxdescript
      type(OutputDataDescript_t),SAVE :: Bsydescript
C     jgf49.49.02: Added variables for 3D fields
      type(OutputDataDescript_t), save :: SigTDescript
      type(OutputDataDescript_t), save :: SalDescript
      type(OutputDataDescript_t), save :: TempDescript
      type(OutputDataDescript_t), save :: RealQDescript
      type(OutputDataDescript_t), save :: ImaginaryQDescript
      type(OutputDataDescript_t), save :: WZDescript
      type(OutputDataDescript_t), save :: Q20Descript
      type(OutputDataDescript_t), save :: LDescript
      REAL(8), ALLOCATABLE, TARGET, SAVE :: rp(:,:) !real part subdomain data
      REAL(8), ALLOCATABLE, TARGET, SAVE :: ip(:,:) !imaginary part subdom dat
      REAL(8), ALLOCATABLE, TARGET, SAVE :: rp_g(:,:) !real part fulldomain
      REAL(8), ALLOCATABLE, TARGET, SAVE :: ip_g(:,:) !imaginary part fulldom
      REAL(8), ALLOCATABLE, TARGET, SAVE :: layer(:) ! horiz layer in 3D
      INTEGER :: numHotstartWrites
      INTEGER :: nextLun

      LOGICAL, SAVE :: FirstCall3D = .true.
      call setMessageSource("hstart3D_out")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     kmd48.33bc write out results in global IO format
      IF (FirstCall3D.eqv..true.) THEN

        ALLOCATE(rp(NP,NFEN),ip(NP,NFEN))

        IF ((MNPROC.GT.1).AND.(myProc.eq.0).AND.
     &     (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
           ALLOCATE(layer(np_g))
        ENDIF
        IF ((MNPROC.eq.1)
     &     .OR.(WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
           ALLOCATE(layer(np))
        ENDIF

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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

        IF ((MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
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


         IF ((NHSTAR.eq.3).or.(NHSTAR.eq.367).or.(NHSTAR.eq.368).or.
     &        (NHSTAR.eq.5).or.(NHSTAR.eq.567).or.(NHSTAR.eq.568)) THEN
#ifdef ADCNETCDF
            ! in parallel, we don't need to create the hotstart file, it was
            ! created by adcprep.
            IF (MNPROC.eq.1) THEN
               CALL initNetCDFHotstart3D(hss%lun, NHSTAR)
               numHotstartWrites = (NT-IT)/NHSINC
               IF (numHotstartWrites.ge.1) THEN
                  IF (hss%lun.eq.67) THEN
                     nextLun = 68
                  ELSE
                     nextLun = 67
                  ENDIF
                  CALL initNetCDFHotstart3D(nextLun, NHSTAR)
               ENDIF
            ENDIF
#endif
         ENDIF

         FirstCall3d = .false.
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

         IF((MYPROC.eq.0).or.
     &      (WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
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
C
         IF ((MNPROC.gt.1).and.(myProc.eq.0).and.
     &        (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
            CALL binaryWrite2D(DUU_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(DUV_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(DVV_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(UU_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(VV_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(BSX_g, np_g, hss%lun, ihotstp)
            CALL binaryWrite2D(BSY_g, np_g, hss%lun, ihotstp)
         ENDIF
         IF ((MNPROC.eq.1).or.
     &      (WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
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
            IF (myProc.eq.0) THEN
               CALL writeNetCDFHotstart3D(hss%lun,DUUDescript,
     &            DUVDescript, DVVDescript, UUDescript, VVDescript,
     &            BSXDescript, BSYDescript)
            ENDIF
#endif
         CASE DEFAULT

      END SELECT
C
C     3 D  V E L O C I T Y
      IF ((MNPROC.GT.1).and.
     &   (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
         IF (myProc.eq.0) THEN
            ALLOCATE(rp_g(NP_G,NFEN),ip_g(NP_G,NFEN),WZ_g(NP_G,NFEN))
            RealQDescript % array2D_g => rp_g
            ImaginaryQDescript % array2D_g => ip_g
            WZDescript % array2D_g => WZ_g
         ENDIF
         rp = real(q)
         CALL collectFullDomainArray(RealQDescript,
     &      packNPbyM, unpackNPbyM)
         ip = aimag(q)
         CALL collectFullDomainArray(ImaginaryQDescript,
     &      packNPbyM, unpackNPbyM)
         CALL collectFullDomainArray(WZDescript,
     &      packNPbyM, unpackNPbyM)
      ENDIF

      SELECT CASE(NHSTAR)
         CASE(-1,1,67,68) ! non-portable binary   !tcm v51.26 mod for time-stamped
            IF ((MNPROC.eq.1).or.
     &         (WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
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
            IF ((MNPROC.GT.1).AND.(myProc.eq.0)
     &         .AND.(WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
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
            IF ((MNPROC.eq.1).or.
     &         (WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
               rp = real(q)
               ip = aimag(q)
            ENDIF
            IF (myProc.eq.0) THEN
               CALL writeNetCDFHotstart3DVar(hss%lun, RealQDescript)
               CALL writeNetCDFHotstart3DVar(hss%lun,ImaginaryQDescript)
               CALL writeNetCDFHotstart3DVar(hss%lun, WZDescript)
            ENDIF
#endif
         CASE DEFAULT

      END SELECT

      IF ((MNPROC.GT.1).and.
     &    (WRITE_LOCAL_HOT_START_FILES.eqv..false.).AND.
     &     (myProc.eq.0)) THEN
            DEALLOCATE(rp_g,ip_g,WZ_g)
      ENDIF
C
C     3 D  T U R B U L E N C E
      IF ((MNPROC.GT.1).and.
     &   (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
         IF (myProc.eq.0) THEN
            ALLOCATE(q20_g(NP_G,NFEN),l_g(NP_G,NFEN))
            Q20Descript % array2D_g => q20_g
            LDescript % array2D_g => l_g
         ENDIF
         CALL collectFullDomainArray(Q20Descript,
     &      packNPbyM, unpackNPbyM)
         CALL collectFullDomainArray(LDescript,
     &      packNPbyM, unpackNPbyM)
      ENDIF
      SELECT CASE(NHSTAR)
         CASE(-1,1,67,68) ! non-portable binary  !tcm v51.26 mod for time-stamped hotstarts (-1)
            IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
               DO N=1,NFEN
                  CALL binaryWrite2D(q20(:,N), np, hss%lun, ihotstp)
               END DO
               DO N=1,NFEN
                  CALL binaryWrite2D(l(:,N), np, hss%lun, ihotstp)
               END DO
            ENDIF
            IF ((MNPROC.GT.1).AND.(myProc.eq.0)
     &         .AND.(WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
               DO N=1,NFEN
                  CALL binaryWrite2D(q20_g(:,N), np_g, hss%lun, ihotstp)
               END DO
               DO N=1,NFEN
                  CALL binaryWrite2D(l_g(:,N), np_g, hss%lun, ihotstp)
               END DO
            ENDIF
         CASE(3,367,368,5,567,568)
#ifdef ADCNETCDF
            IF (myProc.eq.0) THEN
               CALL writeNetCDFHotstart3DVar(hss%lun, Q20Descript)
               CALL writeNetCDFHotstart3DVar(hss%lun, LDescript)
            ENDIF
#endif
         CASE DEFAULT

      END SELECT
      IF ((MNPROC.GT.1).and.
     &   (WRITE_LOCAL_HOT_START_FILES.eqv..false.).AND.
     &    (myProc.eq.0)) THEN
         DEALLOCATE(q20_g,l_g)
      ENDIF
C
C     3 D  D E N S I T Y
      IF ((MNPROC.GT.1).and.
     &   (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
         IF (myProc.eq.0) THEN
            IF (ABS(IDEN).EQ.1) THEN
               ALLOCATE(SigT_g(NP_G,NFEN))
               SigTDescript % array2D_g => SigT_g
            ENDIF
            IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
               ALLOCATE(Sal_g(NP_G,NFEN))
               SalDescript % array2D_g => Sal_g
            ENDIF
            IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
               ALLOCATE(Temp_g(NP_G,NFEN))
               TempDescript % array2D_g => Temp_g
            ENDIF
         ENDIF
         IF (ABS(IDEN).EQ.1) THEN
            CALL collectFullDomainArray(SigTDescript,
     &         packNPbyM, unpackNPbyM)
         ENDIF
         IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
            CALL collectFullDomainArray(SalDescript,
     &         packNPbyM, unpackNPbyM)
         ENDIF
         IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
            CALL collectFullDomainArray(TempDescript,
     &         packNPbyM, unpackNPbyM)
         ENDIF
      ENDIF
C
      SELECT CASE(NHSTAR)
         CASE(-1,1,67,68) ! non-portable binary  !tcm v51.26 mod for time-stamped 
            IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
               IF (ABS(IDEN).EQ.1) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(SigT(:,N), np, hss%lun, ihotstp)
                  END DO
               ENDIF
               IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(Sal(:,N), np, hss%lun, ihotstp)
                  END DO
               ENDIF
               IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(Temp(:,N), np, hss%lun, ihotstp)
                  END DO
               ENDIF
            ENDIF
            IF ((MNPROC.GT.1).AND.(myProc.eq.0)
     &         .AND.(WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
               IF (ABS(IDEN).EQ.1) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(SigT_g(:,N), np_g, hss%lun,
     &                  ihotstp)
                  END DO
               ENDIF
               IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(Sal_g(:,N), np_g, hss%lun,
     &                  ihotstp)
                  END DO
               ENDIF
               IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
                  DO N=1,NFEN
                     CALL binaryWrite2D(Temp_g(:,N), np_g, hss%lun,
     &                  ihotstp)
                  END DO
               ENDIF
            ENDIF
         CASE(3,367,368,5,567,568)
#ifdef ADCNETCDF
         IF (myProc.eq.0) THEN
            IF (ABS(IDEN).EQ.1) THEN
               CALL writeNetCDFHotstart3DVar(hss%lun, SigTDescript)
            ENDIF
            IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
               CALL writeNetCDFHotstart3DVar(hss%lun, SalDescript)
            ENDIF
            IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
               CALL writeNetCDFHotstart3DVar(hss%lun, TempDescript)
            ENDIF
         ENDIF
#endif
      END SELECT
C
      IF ((MNPROC.GT.1).and.
     &   (WRITE_LOCAL_HOT_START_FILES.eqv..false.)) THEN
         IF (myProc.eq.0) THEN
            IF (ABS(IDEN).EQ.1) THEN
               DEALLOCATE(SigT_g)
            ENDIF
            IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
               DEALLOCATE(Sal_g)
            ENDIF
            IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
               DEALLOCATE(Temp_g)
            ENDIF
         ENDIF
      ENDIF
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE HSTART3D_OUT
C-----------------------------------------------------------------------
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
!
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



C****************************************************************************
C   Subroutine to compute 3D SigmaT fields from 3D salinity
C   and/or temperature fields
C
C                                    R.L.  6/22/05
C                                    K.D.  2/03/06 updated
C*****************************************************************************
C

      SUBROUTINE CALC_SIGMAT_3D ()

!Kendra: Added in information needed for the calculation of sigma-t.
!        Currently, I have added only the equation of state from
!        Cushman-Roisin book entitled "Introduction to Geophysical
!        Fluid Dynamics". Later, we can update this to include more
!        complex equations of state.
C   kmd48.33bc need to add variables for new equation of state
      USE GLOBAL_3DVS, ONLY : SAL, TEMP, SIGT, IDEN, NFEN, 
     &                        SIGMA, EQNSTATE, A, B
      use mod_logging, only: setMessageSource,
     &    unsetMessageSource, allMessage, logMessage, DEBUG, ECHO, INFO
Cjgfdebug      USE GLOBAL, ONLY : NP,RHOWAT0, ScreenUnit
      USE GLOBAL, ONLY : ETA2, IFNLFA
      USE ADC_CONSTANTS, ONLY: RHOWAT0, G
      USE MESH, ONLY : NP, DP
      IMPLICIT NONE
      INTEGER NH,N
      REAL(8) :: T_alpha=0.00017d0, S_beta=0.00076d0
!      REAL(8) :: T_ref=10.0d0, S_ref=35.0d0
      REAL(8) :: T_ref, S_ref ! changed to new values
!Kendra: Parameters for the 2nd equation of state
      REAL(8), PARAMETER :: Q1 = 196637339.d0
      REAL(8), PARAMETER :: Q2 = 196668.928d0
      REAL(8), PARAMETER :: a1 = 1446045.44d0
      REAL(8), PARAMETER :: a2 = -10769.1980d0
      REAL(8), PARAMETER :: a3 = 85.2955498d0
      REAL(8), PARAMETER :: b1 = 579854.265d0
      REAL(8), PARAMETER :: b2 = -1477.31528d0
      REAL(8), PARAMETER :: b3 = 419.489109d0
      REAL(8), PARAMETER :: c1 = 2001.22349d0
      REAL(8), PARAMETER :: c2 = 0.0213851025d0
      REAL(8), PARAMETER :: c3 = 0.968784141d0
      REAL(8), PARAMETER :: c4 = -0.00654785602d0
      REAL(8), PARAMETER :: c5 = -0.00000250468726d0
      REAL(8), PARAMETER :: c6 = 0.0000000628902345d0
      REAL(8), PARAMETER :: c7 = 0.00000000282414763d0
      REAL(8), PARAMETER :: d1 = 1433.02205d0
      REAL(8), PARAMETER :: d2 = -9.09231525d0
      REAL(8), PARAMETER :: d3 = 0.0791654429d0
      REAL(8), PARAMETER :: d4 = 0.0000398630534d0
      REAL(8), PARAMETER :: e1 = 417.831720d0
      REAL(8), PARAMETER :: e2 = -1.87581316d0
      REAL(8), PARAMETER :: e3 = -0.0000387902837d0
      REAL(8), PARAMETER :: f1 = 1.00765828d0
      REAL(8), PARAMETER :: f2 = 0.000312912597d0
      REAL(8), PARAMETER :: T_max=40.d0, T_min=0.d0
      REAL(8), PARAMETER :: S_max=42.d0, S_min=0.d0
!Kendra: Variables for the 2nd equation of state
      REAL(8) :: P1, P2
      REAL(8) :: p(NP), RHO(NP,NFEN)
      REAL(8) :: R_ref
      REAL(8) :: Hs, HsOAMB, Zk
      REAL(8) :: Z(NP,NFEN)
!Kendra: Parameters for the 3rd equation of state
      REAL(8), PARAMETER :: c00=999.842594d0, c11=0.06793952d0
      REAL(8), PARAMETER :: c12=-0.00909529d0, c13=0.0001001685d0
      REAL(8), PARAMETER :: c14=-0.000001120083d0, c15=0.000000006536332d0
      REAL(8), PARAMETER :: c21=0.824493d0, c22=-0.0040899d0
      REAL(8), PARAMETER :: c23=0.000076438d0, c24=0.00000082467d0
      REAL(8), PARAMETER :: c25=0.0000000053875d0, c31=-0.00572466d0
      REAL(8), PARAMETER :: c32=0.00010227d0, c33=-0.0000016546d0
      REAL(8), PARAMETER :: c41=0.00048314d0
      REAL(8), PARAMETER :: d00=190925.6d0, d11=2098.925d0
      REAL(8), PARAMETER :: d12=-30.41638d0, d13=-0.01852732d0
      REAL(8), PARAMETER :: d14=-0.0001361629d0, d21=1044.077d0
      REAL(8), PARAMETER :: d22=-65.00517d0, d23=1.553190d0
      REAL(8), PARAMETER :: d24=0.002326469d0, d31=-55.87545d0
      REAL(8), PARAMETER :: d32=7.390729d0, d33=-0.1909078d0
      REAL(8), PARAMETER :: d41=-4.721788d0, d42=-0.1028859d0
      REAL(8), PARAMETER :: d43=0.002512549d0, d44=0.000005939910d0
      REAL(8), PARAMETER :: d51=0.1571896d0, d52=0.002598241d0
      REAL(8), PARAMETER :: d53=0.00007267926d0, d61=-0.02042967d0
      REAL(8), PARAMETER :: d71=0.0001045941d0, d72=-0.000000005782165d0
      REAL(8), PARAMETER :: d73=0.000001296821d0, d81=-0.000002595994d0
      REAL(8), PARAMETER :: d82=-0.00000001248266d0, d83=-0.00000003508914d0
!Kendra: Variables for the 3rd equation of state
      REAL(8) :: Kappa

      call setMessageSource("calc_sigmat_3D")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


!Kendra: Equation of state from Cushman-Roisin book entitled "Introduction to
!        Geophysical Fluid Dynamics".
      IF (Eqnstate.EQ.1) THEN
       T_ref=10.d0
       S_ref=35.d0
       IF(ABS(IDEN).EQ.2) THEN
        DO NH=1,NP
          DO N=1,NFEN
             SIGT(NH,N) = 1028.d0*(1.d0+S_beta*(SAL(NH,N)-
     &                    S_ref))-RHOWAT0
          ENDDO
        ENDDO
       ELSEIF(ABS(IDEN).EQ.3) THEN
        DO NH=1,NP
          DO N=1,NFEN
             SIGT(NH,N) = 1028.d0*(1.d0-T_alpha*(TEMP(NH,N)-
     &                    T_ref))-RHOWAT0
          ENDDO
        ENDDO
       ELSEIF(ABS(IDEN).EQ.4) THEN
        DO NH=1,NP
          DO N=1,NFEN
             SIGT(NH,N) = 1028.d0*(1.d0-T_alpha*(TEMP(NH,N)-
     &                    T_ref)+S_beta*(SAL(NH,N)-S_ref))-
     &                    RHOWAT0
          ENDDO
        ENDDO
       ENDIF
      ELSE IF (Eqnstate.EQ.2) THEN
        T_ref=10.d0
        S_ref=32.d0
C  kmd48.33bc Equation of state from McDougall et al. (2003) - only used
C        with both salinity and temperature. This is one of the equation
C        of states used in QUODDY
        IF (ABS(IDEN).EQ.4) THEN
! Look at the values for Temperature and Salinity to see if still
! inside the range of the equation of state
           DO NH=1,NP
            DO N=1,NFEN
               IF (TEMP(NH,N).LT.T_min) THEN
                   TEMP(NH,N) = T_min
               END IF
               IF (TEMP(NH,N).GT.T_max) THEN
                   TEMP(NH,N) = T_max
               END IF
               IF (SAL(NH,N).LT.S_min) THEN
                   SAL(NH,N) = S_min
               END IF
               IF (SAL(NH,N).GT.S_max) THEN
                   SAL(NH,N) = S_max
               END IF
            END DO
           END DO
! First, determine the reference density
           P1 = Q1 + (a1 + (a2 + a3*T_ref)*T_ref)*
     &               T_ref + (b1 + b2*T_ref +
     &               b3*S_ref)*S_ref
           P2 = Q2 + (d1 + (d2 + (d3 + d4*T_ref)*
     &               T_ref)*T_ref)*T_ref +
     &               (e1 + (f1 + f2*T_ref*T_ref)*
     &               SQRT(S_ref) + (e2 + e3*T_ref*
     &               T_ref)*T_ref)*S_ref
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
              IF (N.EQ.NFEN) THEN
                p(NH) = 0.d0
                P1 = Q1 + (a1 + (a2 + a3*TEMP(NH,N))*TEMP(NH,N))*
     &               TEMP(NH,N) + (b1 + b2*TEMP(NH,N) +
     &               b3*SAL(NH,N))*SAL(NH,N)
                P2 = Q2 + (d1 + (d2 + (d3 + d4*TEMP(NH,N))*
     &               TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) +
     &               (e1 + (f1 + f2*TEMP(NH,N)*TEMP(NH,N))*
     &               SQRT(SAL(NH,N)) + (e2 + e3*TEMP(NH,N)*
     &               TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)
                RHO(NH,N) = P1/P2
              ELSE IF (N.NE.NFEN) THEN
                p(NH) = p(NH) + g*RHO(NH,N+1)*(Z(NH,N+1)-
     &                  Z(NH,N))*0.0001d0
                P1 = Q1 + (a1 + (a2 + a3*TEMP(NH,N))*TEMP(NH,N))*
     &               TEMP(NH,N) + (b1 + b2*TEMP(NH,N) +
     &               b3*SAL(NH,N))*SAL(NH,N) + (c1 + c2*
     &               TEMP(NH,N)*TEMP(NH,N) + c3*SAL(NH,N) +
     &               (c4 + ((c5+c6*TEMP(NH,N))*TEMP(NH,N) +
     &               c7*p(NH))*TEMP(NH,N))*p(NH))*p(NH)
                P2 = Q2 + (d1 + (d2 + (d3 + d4*TEMP(NH,N))*
     &               TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) +
     &               (e1 + (f1 + f2*TEMP(NH,N)*TEMP(NH,N))*
     &               SQRT(SAL(NH,N)) + (e2 + e3*TEMP(NH,N)*
     &               TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)
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
      ELSE IF(Eqnstate.EQ.3) THEN
        T_ref=10.d0
        S_ref=32.d0
C  kmd48.33bc Equation of state from UNESCO (1980) and Haidvogel and Beckman (1999) - only used
C        with both salinity and temperature.
        IF (ABS(IDEN).EQ.4) THEN
! Look at the values for Temperature and Salinity to see if still
! inside the range of the equation of state
           DO NH=1,NP
            DO N=1,NFEN
               IF (TEMP(NH,N).LT.T_min) THEN
                   TEMP(NH,N) = T_min
               END IF
               IF (TEMP(NH,N).GT.T_max) THEN
                   TEMP(NH,N) = T_max
               END IF
               IF (SAL(NH,N).LT.S_min) THEN
                   SAL(NH,N) = S_min
               END IF
               IF (SAL(NH,N).GT.S_max) THEN
                   SAL(NH,N) = S_max
               END IF
            END DO
           END DO
! Kendra: the initial step is to define the density for the surface
           DO NH=1,NP
              N=NFEN
              RHO(NH,N) = c00 + (c11 + (c12 +(c13 +(c14 + c15*
     &         TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*
     &         TEMP(NH,N) + (c21 + (c22 + (c23 + (c24 + c25*
     &         TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*
     &         SAL(NH,N) + (c31 + (c32 + c33*TEMP(NH,N))*
     &         TEMP(NH,N))*SAL(NH,N)*SQRT(SAL(NH,N)) +
     &         c41*SAL(NH,N)*SAL(NH,N)
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
              Kappa = d00 + (d11 + (d12 + (d13 + d14*TEMP(NH,N))*
     &          TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N) + (d21 +
     &          (d22 + (d23 + d24*TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*
     &          SAL(NH,N) + (d31 + (d32 + d33*TEMP(NH,N))*TEMP(NH,N))*
     &          SAL(NH,N)*SQRT(SAL(NH,N)) + (d41 + (d42 + (d43 + d44*
     &          TEMP(NH,N))*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N) + (d51 +
     &          (d52 + d53*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N)*SAL(NH,N) +
     &          d61*Z(NH,N)*SAL(NH,N)*SQRT(SAL(NH,N)) + (d71 + (d72 +
     &          d73*TEMP(NH,N))*TEMP(NH,N))*Z(NH,N)*Z(NH,N) + (d81 +
     &          (d82 + d83*TEMP(NH,N))*TEMP(NH,N))*SAL(NH,N)*Z(NH,N)*
     &          Z(NH,N)
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
C*****************************************************************************
      END SUBROUTINE CALC_SIGMAT_3D
C*****************************************************************************

!******************************************************************************
!    loop on all nodes, and apply a 1-layer Gaussian filter for velocities
! arash: December 9 2015
!******************************************************************************

       Subroutine Gaussian_filter_complex ( X, Y, NeiTab, NNeigh, quantity, NP, NFEN, MNEI )

       Use Global,      Only : IFNLFA, ETA2, nodecode
       Use Mesh,        Only : DP
       Use Global_3DVS, Only : SIGMA, AMB, B, Biharmonic_viscosity_modified_Leith
       Use mod_terminate, Only: terminate
       Use sizes, Only: Myproc

       Implicit none
       Integer :: I, J, K
       Integer :: Neighbor_Node
       Integer :: MNEI, NP, NFEN
       Real(8) :: sigma_Gauss, length, length_max, factor, value, numerator, denominator, sum_Weights
       Real(8) :: velocity_magnitude
       Real(8) :: velocity_magnitude_max_in_xy
       Real(8) :: element_peclet_number
       Real(8) :: constant_diffusivity = 50.0d0
       Real(8) :: average_viscosity, average_velocity, average_element_peclet_number
       Real(8) :: low_peclet_threshold, high_peclet_threshold, low_factor, high_factor
       Real(8) :: X(NP), Y(NP)
       INTEGER :: NeiTab(NP,MNEI), NNeigh(NP)
       REAL(8), ALLOCATABLE :: Weights(:)
       ! velocities are stored as complex quantities: u+i*v
       Complex(8), ALLOCATABLE :: quantity_smooth(:,:) 
       Complex(8), ALLOCATABLE :: inarray_vertical (:)
       Complex(8) :: quantity(NP,NFEN)
       Complex(8) :: interpolated_value

       REAL(8)  :: Hs, Hs_Neighbor_Node     !Total water depth at time level s
       REAL(8)  :: HsOAMB                   !Hs/(a-b)
       REAL(8)  :: Zk                       ! z depth of any node k in the vertical
       Real(8)  :: SigmaNN

       Allocate ( quantity_smooth(NP,NFEN) )
       Allocate ( inarray_vertical(NFEN) )

!---------------------------------------------------------------------------------------------------
! arash 21 July 2016: adaptive filtering, depending on local Peclet number
! I tuned these parameters for the SURA/SL16 mesh based on trial and error, pain, and engineering judgment.
       low_peclet_threshold  =  1.0d0;         low_factor = 0.15
       high_peclet_threshold =  16.0d0;        high_factor = 0.185 ! filter 1
!       high_peclet_threshold =  18.0d0;        high_factor = 0.185 ! filter 2
!       high_peclet_threshold =  4.0d0;        high_factor = 0.185 ! filter 3
!       high_peclet_threshold =  20.0d0;        high_factor = 0.185 ! filter 1
!---------------------------------------------------------------------------------------------------

       quantity_smooth(:,:) = (0.0d0 , 0.0d0)

       ! loop on all nodes, and apply the Gaussian filter
       Do I = 1, NP
          ! @jasonfleming : this loop does not take into account that a 
          ! node may be dry so I've added a statement to skip dry nodes
          if (nodecode(i).eq.0) cycle
          !
          ! 1. compute maximum distance between the node of interest, and the neighbor nodes
          !    this will be used to compute sigma2
          length_max = 0.0d0
          Do J = 2, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             length = SQRT ( ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 )
             If ( length > length_max ) length_max = length
          End Do
          ! 1.2. compute the maximum velocity in the xy plane
          velocity_magnitude_max_in_xy = 0.0d0
          average_velocity = 0.0d0
          Do K = 1, NFEN
             velocity_magnitude = sqrt( (real ( quantity( I , K ) ))**2 + (aimag ( quantity( I , K ) ))**2 )
             If ( velocity_magnitude > velocity_magnitude_max_in_xy ) velocity_magnitude_max_in_xy = velocity_magnitude
             average_velocity = average_velocity + velocity_magnitude
          End Do
          average_velocity = average_velocity / NFEN
          ! 1.3. compute the average viscosity in the xy plane
          average_viscosity = 0.0d0
          Do K = 1, NFEN
             !write(6,*) 'i, k, biharmonic_viscosity_modified_leith ',i, k, Biharmonic_viscosity_modified_Leith(i,k) !jgfdebug
             average_viscosity = average_viscosity + Biharmonic_viscosity_modified_Leith ( I, K )
          End Do
          average_viscosity = average_viscosity / NFEN

          ! 1.3.1 compute the element Peclet number
          !write(6,'("i=",i0," nodecode=",i0," average_viscosity=",f8.2)') i, nodecode(i), average_viscosity !jgfdebug
          ! @jasonfleming : the average viscosity will be zero if the node
          ! is dry, producing a division by zero error here. So I've added
          ! a statement to skip dry nodes at the top of this loop.
          element_peclet_number = length_max * velocity_magnitude_max_in_xy / average_viscosity !constant_diffusivity !average_viscosity
          !average_element_peclet_number = length_max * average_velocity / average_viscosity 
          average_element_peclet_number = length_max * velocity_magnitude_max_in_xy / average_viscosity 

          ! 1.4. specify sigma - we want to do this adaptively.
          ! I tuned this parameter for the SL16 mesh based on trial and error, pain, and engineering judgment.

             If ( average_element_peclet_number < low_peclet_threshold ) Then
                factor = low_factor
             Else If ( average_element_peclet_number > high_peclet_threshold ) Then
                factor = high_factor
             Else
                factor = low_factor + ( high_factor - low_factor ) * 
     &           ( ( average_element_peclet_number - low_peclet_threshold ) / ( high_peclet_threshold - low_peclet_threshold ) )**0.5
             End If
   
          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0
!----------------------------------------------------------------------------------------
! 43_SL16
!          If ( element_peclet_number >= 80. ) Then
!             factor = 0.19d0
!          Else If ( (element_peclet_number < 80.).And.(element_peclet_number >= 30.) ) Then
!             factor = 0.18d0
!          Else If ( (element_peclet_number < 30.).And.(element_peclet_number >= 10.) ) Then
!             factor = 0.17d0
!          Else If ( (element_peclet_number < 10.).And.(element_peclet_number >= 1.) ) Then
!             factor = 0.16d0
!          Else
!             factor = 0.15d0
!          End If
!         ! arash: this may not be necessary
!          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
! 10-SURA6
!          If ( element_peclet_number >= 40. ) Then
!             factor = 0.19d0
!          Else If ( (element_peclet_number < 40.).And.(element_peclet_number >= 20.) ) Then
!             factor = 0.18d0
!          Else If ( (element_peclet_number < 20.).And.(element_peclet_number >= 10.) ) Then
!!             factor = 0.17d0
!          Else If ( (element_peclet_number < 10.).And.(element_peclet_number >= 5.) ) Then
!             factor = 0.16d0
!          Else
!             factor = 0.15d0
!          End If
!          ! arash: this may not be necessary
!          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
! 5,6-SURA6 poly
!          If ( element_peclet_number >= 20. ) Then
!             factor = 0.19d0
!          Else If ( (element_peclet_number < 20.).And.(element_peclet_number >= 10.) ) Then
!             factor = 0.18d0
!          Else If ( (element_peclet_number < 10.).And.(element_peclet_number >= 5.) ) Then
!             factor = 0.17d0
!          Else If ( (element_peclet_number < 5.).And.(element_peclet_number >= 1.) ) Then
!             factor = 0.16d0
!          Else
!             factor = 0.15d0
!          End If
         ! arash: this may not be necessary
!          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0
!----------------------------------------------------------------------------------------
          sigma_Gauss = length_max * factor

          ! 2. compute Gaussian filter weights
          IF ( ALLOCATED(Weights) ) DEALLOCATE( Weights )
          ALLOCATE( Weights( NNeigh ( I ) ) )
          Weights = 0.0d0
          Do J = 1, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             numerator = ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 
             denominator = 2.0d0 * sigma_Gauss**2
             Weights (J) = exp ( -1.0 * numerator / denominator )
          End Do
          ! 2.1 normalize weights
          sum_Weights = 0.0d0
          Do J = 1, NNeigh ( I )
             sum_Weights = sum_Weights + Weights (J)
          End Do
          If ( sum_Weights < 1.0d-6 ) Then
             Write(*,*) 'Possible error in Gaussian filter - error code = 1'
             call terminate(myproc)
          End If
         Weights (:) = Weights (:) / sum_Weights
         ! 2.3 debug
         !Write(*,'(I5,11f6.3)') NNeigh ( I ), weights(1:NNeigh ( I ))
         !Pause
         ! 2.4. compute elevation corresponding to this node
         Hs = DP(I) + IFNLFA * Eta2(I)  !Total depth at current timestep - 'center' node
         ! arash: not sure if the next line is necessary, but probably not. we aren't using it in the bpg routine either.
         ! IF( Hs .LT. H0 ) Hs = H0
         HsOAMB = Hs / AMB

         ! 3. apply the Gaussian filter
          Do J = 1, NNeigh ( I )
             Neighbor_Node = NeiTab ( I , J )
             inarray_vertical (:) = quantity( Neighbor_Node , : )
             Do K = 1, NFEN
                ! 3.0. obtain the z and sigma value corresponding to the 'center' node
                Zk = HsOAMB * ( Sigma(K) - B ) - DP(I)
                Hs_Neighbor_Node = DP(Neighbor_Node) + IFNLFA * Eta2(Neighbor_Node)
                SigmaNN = B + AMB * ( Zk + DP(Neighbor_Node) ) / Hs_Neighbor_Node
                ! 3.1. linearly interpolate surrounding velocities on the same horizontal plane
                Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, interpolated_value, K )
                quantity_smooth (I,K) = quantity_smooth (I,K) + Weights (J) * interpolated_value
             End Do
          End Do
       End Do

       ! transfer smooth bathymetry to the original array
       quantity (:,:) = quantity_smooth (:,:)

       End Subroutine


!******************************************************************************
!    loop on all nodes, and apply a 1-layer Gaussian filter for velocities
! arash: December 9 2015
!******************************************************************************

       Subroutine Gaussian_filter_complex_new ( X, Y, NeiTab, NNeigh, quantity, NP, NFEN, MNEI )

       Use Global,      Only : IFNLFA, ETA2
       USE Mesh,        Only : DP
       Use Global_3DVS, Only : SIGMA, AMB, B, Biharmonic_viscosity_modified_Leith
       Use mod_terminate, Only: terminate
       Use sizes, only: myproc

       Implicit none
       Integer :: I, J, K
       Integer :: Neighbor_Node
       Integer :: MNEI, NP, NFEN
       !Real(8) :: sigma_Gauss, length, length_max, factor, value, numerator, denominator, sum_Weights
       Real(8) :: length, length_max, numerator
       Real(8) :: velocity_magnitude
       Real(8) :: velocity_magnitude_max_in_xy
!       Real(8) :: element_peclet_number
       Real(8) :: average_viscosity
       Real(8) :: low_peclet_threshold, high_peclet_threshold, low_factor, high_factor
       Real(8) :: X(NP), Y(NP)
       INTEGER :: NeiTab(NP,MNEI), NNeigh(NP)
       REAL(8), ALLOCATABLE :: Weights(:,:)
       REAL(8), ALLOCATABLE :: sigma_Gauss(:)
       REAL(8), ALLOCATABLE :: factor(:)
       REAL(8), ALLOCATABLE :: denominator(:)
       REAL(8), ALLOCATABLE :: sum_Weights(:)
       REAL(8), ALLOCATABLE :: element_peclet_number(:)
       ! velocities are stored as complex quantities: u+i*v
       Complex(8), ALLOCATABLE :: quantity_smooth(:,:) 
       Complex(8), ALLOCATABLE :: inarray_vertical (:)
       Complex(8) :: quantity(NP,NFEN)
       Complex(8) :: interpolated_value

       REAL(8)  :: Hs, Hs_Neighbor_Node     !Total water depth at time level s
       REAL(8)  :: HsOAMB                   !Hs/(a-b)
       REAL(8)  :: Zk                       ! z depth of any node k in the vertical
       Real(8)  :: SigmaNN

!---------------------------------------------------------------------------------------------------
! arash 21 July 2016: adaptive filtering, depending on local Peclet number
! I tuned these parameters for the SURA/SL16 mesh based on trial and error, pain, and engineering judgment.
       low_peclet_threshold  =  1.0d0;         low_factor = 0.15
       high_peclet_threshold =  10.0d0;        high_factor = 0.185
!---------------------------------------------------------------------------------------------------

       Allocate ( quantity_smooth(NP,NFEN) )
       Allocate ( inarray_vertical(NFEN) )
       Allocate ( sigma_Gauss(NFEN) )
       Allocate ( factor(NFEN) )
       Allocate ( denominator(NFEN) )
       Allocate ( sum_Weights(NFEN) )
       Allocate ( element_peclet_number(NFEN) )

       quantity_smooth(:,:) = (0.0d0 , 0.0d0)

       ! loop on all nodes, and apply the Gaussian filter
       Do I = 1, NP
          ! 1. compute maximum distance between the node of interest, and the neighbor nodes
          !    this will be used to compute sigma2
          length_max = 0.0d0
          Do J = 2, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             length = SQRT ( ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 )
             If ( length > length_max ) length_max = length
          End Do
          ! 1.2. compute the following: maximum velocity in the xy plane
          velocity_magnitude_max_in_xy = 0.0d0
          average_viscosity = 0.0d0
          Do K = 1, NFEN
             velocity_magnitude = sqrt( (real ( quantity( I , K ) ))**2 + (aimag ( quantity( I , K ) ))**2 )
             element_peclet_number(K) = length_max * velocity_magnitude / Biharmonic_viscosity_modified_Leith ( I, K )
!             If ( velocity_magnitude > velocity_magnitude_max_in_xy ) velocity_magnitude_max_in_xy = velocity_magnitude
!             average_viscosity = average_viscosity + Biharmonic_viscosity_modified_Leith ( I, K )
!          End Do
!          average_viscosity = average_viscosity / NFEN
          ! 1.3. compute the avrage element Peclet number
          !element_peclet_number = length_max * velocity_magnitude_max_in_xy / average_viscosity
          ! 1.4. specify sigma - we want to do this adaptively.
!          Do K = 1, NFEN
             If ( element_peclet_number(K) < low_peclet_threshold ) Then
                factor(K) = low_factor
             Else If ( element_peclet_number(K) > high_peclet_threshold ) Then
                factor(K) = high_factor
             Else
                factor(K) = low_factor + ( high_factor - low_factor ) * 
     &           ( ( element_peclet_number(K) - low_peclet_threshold ) / ( high_peclet_threshold - low_peclet_threshold ) )**0.5

             End If
             ! limit factor for high velocities -------------------------------
             If ( velocity_magnitude > 1.8d0 ) factor(K) = high_factor
             !-----------------------------------------------------------------
             sigma_Gauss(K) = length_max * factor(K)
          End Do

          ! 2. compute Gaussian filter weights
          IF ( ALLOCATED(Weights) ) DEALLOCATE( Weights )
          ALLOCATE( Weights( NNeigh ( I ), NFEN ) )
          Weights(:,:) = 0.0d0
          Do K = 1, NFEN
             Do J = 1, NNeigh ( I ) ! also includes node I
                Neighbor_Node = NeiTab ( I , J )
                numerator = ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 
                denominator(K) = 2.0d0 * sigma_Gauss(K)**2
                Weights (J,K) = exp ( -1.0 * numerator / denominator(K) )
             End Do
             ! 2.1 normalize weights
             sum_Weights(:) = 0.0d0
             Do J = 1, NNeigh ( I )
                sum_Weights(K) = sum_Weights(K) + Weights (J,K)
             End Do
             If ( sum_Weights(K) < 1.0d-6 ) Then
                Write(*,*) 'Possible error in Gaussian filter - error code = 1'
                call terminate(myproc)
             End If
            Weights (:,K) = Weights (:,K) / sum_Weights(K)
         End Do
         ! 2.3 debug
         !Write(*,'(I5,11f6.3)') NNeigh ( I ), weights(1:NNeigh ( I ))
         !Pause
         ! 2.4. compute elevation corresponding to this node
         Hs = DP(I) + IFNLFA * Eta2(I)  !Total depth at current timestep - 'center' node
         ! arash: not sure if the next line is necessary, but probably not. we aren't using it in the bpg routine either.
         ! IF( Hs .LT. H0 ) Hs = H0
         HsOAMB = Hs / AMB

         ! 3. apply the Gaussian filter
          Do J = 1, NNeigh ( I )
             Neighbor_Node = NeiTab ( I , J )
             inarray_vertical (:) = quantity( Neighbor_Node , : )
             Do K = 1, NFEN
                ! 3.0. obtain the z and sigma value corresponding to the 'center' node
                Zk = HsOAMB * ( Sigma(K) - B ) - DP(I)
                Hs_Neighbor_Node = DP(Neighbor_Node) + IFNLFA * Eta2(Neighbor_Node)
                SigmaNN = B + AMB * ( Zk + DP(Neighbor_Node) ) / Hs_Neighbor_Node
                ! 3.1. linearly interpolate surrounding velocities on the same horizontal plane
                Call linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, interpolated_value, K )
                quantity_smooth (I,K) = quantity_smooth (I,K) + Weights (J,K) * interpolated_value
             End Do
          End Do
       End Do

       ! transfer smooth bathymetry to the original array
       quantity (:,:) = quantity_smooth (:,:)

       DeAllocate ( quantity_smooth )

       End Subroutine

C*************************************************************************
C     Subroutine to interpolate velocities to a specified
C     sigma value (SigmaNN) given an initial guess of which sigma
C     level is closest to the specified value.
C     This subroutine is almost identical to SUBROUTINE ZSURFBUOY, written by R.L.
C
! please pay particular attention to what we are trying to interpolate, and for what purpose.
! this subroutine is used to provide velocities for filtering (smoothing). Therefore, we oly use linear
! interpolaton. Moreover, if we run into air, or, into the ground, we use the velocity value at the surface, or bottom.
! in other words, we use our engineering judgment in the smoothing procedure as well.
!
C                                    a.f.  09/12/15
! af 01/04/2016:
!I use this subroutine to compute the biharmonic viscosity operator and modified Leith viscosity as well.
!
C*************************************************************************

      Subroutine linear_vertical_interpolation_complex ( inarray_vertical, SigmaNN, interpolated_value, J )

      USE GLOBAL_3DVS
      IMPLICIT NONE
      Complex(8) :: interpolated_value
      REAL(8) :: SigmaNN     !Sigma value of a neighbor node

      INTEGER NN
      INTEGER J
!      INTEGER :: LBelo = -999
!      INTEGER :: LAbov = -999
      INTEGER LBelo
      INTEGER LAbov
      INTEGER LTry
      INTEGER IDiag
      REAL(8) SigBelo
      REAL(8) SigAbov
      REAL(8) SigTry
! af: variables for cubic interpolation
!      Integer  :: k1, k2, k3, k4                 ! nodes used for cubic polynomial
!      REAL(8) :: f1, f2, f3, f4                 ! shape functions
!      REAL(8) :: top, bot
!      Integer  :: loop_situation = 0

      Complex(8) :: inarray_vertical
      Dimension  inarray_vertical ( NFEN )

      call setMessageSource("cubic_vertical_interpolation")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IDiag=0

      IF(SigmaNN.LE.b) THEN !at bottom or lower, then use bottom
         LBelo=1
         interpolated_value = inarray_vertical(LBelo)
         SigBelo=b
         SigAbov=b
         GOTO 100
      ENDIF

      IF(SigmaNN.GE.a) THEN ! surface or above, use surface
         LAbov=NFEN
         interpolated_value=inarray_vertical(LAbov)
         SigBelo=a
         SigAbov=a
         GOTO 100
      ENDIF
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--


      LTry=J                    !start search for SIGABOV and SIGBELO
      SigTry=Sigma(LTry)
      IF(SigmaNN.GT.SigTry) THEN !too low
!         loop_situation = 1
         SigBelo=SigTry         !SIGBELO may = SIGTRY
         LBelo=LTry
         LTry=LTry+1            !look at next level higher
 90      SigTry=Sigma(LTry)
         IF(SigmaNN.GT.SigTry) THEN !still too low
            SigBelo=SigTry
            LBelo=LTry
            LTry=LTry+1
            GOTO 90
         ENDIF
         SigAbov=SigTry         !found upper bracketing sigma
         LAbov=LTry
         GOTO 99                !go interpolate
      ENDIF

      IF(SigmaNN.LE.SigTry) THEN !to high
!         loop_situation = 2
         SigAbov=SigTry         !SIGABOV may = SIGTRY
         LAbov=LTry
         LTry=LTry-1            !look at next level lower
 91      SigTry=Sigma(LTry)
         IF(SigmaNN.LE.SigTry) THEN !still too high
            SigAbov=SigTry
            LAbov=LTry
            LTry=LTry-1
            GOTO 91
         ENDIF
         SigBelo=SigTry         !found lower bracketing sigma
         LBelo=LTry
      ENDIF
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
99    Continue

! use a cubic polynomial interpolation, instead of the above (99) linear interpolation
! 1- decide which nodes to pick: k1:lowest, k2, k3, k4:top-most
!      If ( LAbov == NFEN ) Then
!         k4 = NFEN
!         k3 = NFEN - 1
!         k2 = NFEN - 2
!         k1 = NFEN - 3
!      Else If ( LBelo == 1 ) Then
!         k4 = 4
!         k3 = 3
!         k2 = 2
!         k1 = 1
!      Else
!         k4 = LAbov + 1
!         k3 = LAbov
!         k2 = LBelo
!         k1 = LBelo - 1
!      End If

! 2- compute the shape functions
!        f1 = ((SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k1)-Sigma(k2)) * (Sigma(k1)-Sigma(k3)) * (Sigma(k1)-Sigma(k4)))

!        f2 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k2)-Sigma(k1)) * (Sigma(k2)-Sigma(k3)) * (Sigma(k2)-Sigma(k4)))

!        f3 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k3)-Sigma(k1)) * (Sigma(k3)-Sigma(k2)) * (Sigma(k3)-Sigma(k4)))

!        f4 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)))
!     &     / ((Sigma(k4)-Sigma(k1)) * (Sigma(k4)-Sigma(k2)) * (Sigma(k4)-Sigma(k3)))

! 3- interpolate
!      interpolated_value = f1 * inarray_vertical(k1) + f2 * inarray_vertical(k2) 
!     &          + f3 * inarray_vertical(k3) + f4 * inarray_vertical(k4)

! linear interpolation should be enough. the goal is to smooth the velocity
      interpolated_value = (inarray_vertical(LAbov) - inarray_vertical(LBelo))
     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + inarray_vertical(LBelo)

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
 100  CONTINUE

      IF(IDiag.EQ.2) THEN
         WRITE(2,*) '******** linear_vertical_interpolation **********'
         WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov',
     &        '      interpolated_value'
         WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,interpolated_value
 777     FORMAT(I7,I5,3(F10.3),E14.5)
      ENDIF

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE linear_vertical_interpolation_complex
C**********************************************************************

!******************************************************************************
!    loop on all nodes, and apply a 1-layer Gaussian filter for salinity/temperature
! arash: 04/01/2016
! arash: March 24 2016
!******************************************************************************

       Subroutine Gaussian_filter_real ( X, Y, NeiTab, NNeigh, quantity, NP, NFEN, MNEI )

       Use Global,      Only : IFNLFA, ETA2
       Use Mesh,        Only : DP
       Use Global_3DVS, Only : SIGMA, AMB, B, qkp1, Biharmonic_viscosity_modified_Leith
       Use mod_terminate, only: terminate
       Use sizes, only: myproc

       Implicit none
       Integer :: I, J, K
       Integer :: Neighbor_Node
       Integer :: MNEI, NP, NFEN
       Real(8) :: sigma_Gauss, length, length_max, factor, value, numerator, denominator, sum_Weights
       Real(8) :: X(NP), Y(NP)
       INTEGER :: NeiTab(NP,MNEI), NNeigh(NP)
       REAL(8), ALLOCATABLE :: Weights(:)
       ! velocities are stored as complex quantities: u+i*v
       real(8), ALLOCATABLE :: quantity_smooth(:,:) 
       real(8), ALLOCATABLE :: inarray_vertical (:)
       real(8) :: quantity(NP,NFEN)
       real(8) :: interpolated_value

       Real(8) :: velocity_magnitude
       Real(8) :: velocity_magnitude_max_in_xy
       Real(8) :: element_peclet_number, average_element_peclet_number
       Real(8) :: average_viscosity
       Real(8) :: low_peclet_threshold, high_peclet_threshold, low_factor, high_factor

       REAL(8)  :: Hs, Hs_Neighbor_Node     !Total water depth at time level s
       REAL(8)  :: HsOAMB                   !Hs/(a-b)
       REAL(8)  :: Zk                       ! z depth of any node k in the vertical
       Real(8)  :: SigmaNN

       Allocate ( quantity_smooth(NP,NFEN) )
       Allocate ( inarray_vertical(NFEN) )

!---------------------------------------------------------------------------------------------------
! arash 21 July 2016: adaptive filtering, depending on local Peclet number
! I tuned these parameters for the SURA/SL16 mesh based on trial and error, pain, and engineering judgment.
       low_peclet_threshold  =  1.0d0;         low_factor = 0.15
       high_peclet_threshold =  16.0d0;        high_factor = 0.185
!---------------------------------------------------------------------------------------------------

       quantity_smooth(:,:) = 0.0d0

       ! loop on all nodes, and apply the Gaussian filter
       Do I = 1, NP
          ! 1. compute maximum distance between the node of interest, and the neighbor nodes
          !    this will be used to compute sigma2
          length_max = 0.0d0
          Do J = 2, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             length = SQRT ( ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 )
             If ( length > length_max ) length_max = length
          End Do
          ! 1.2. compute the maximum velocity in the xy plane
          velocity_magnitude_max_in_xy = 0.0d0
          Do K = 1, NFEN
             velocity_magnitude = sqrt( (real ( qkp1( I , K ) ))**2 + (aimag ( qkp1( I , K ) ))**2 )
             If ( velocity_magnitude > velocity_magnitude_max_in_xy ) velocity_magnitude_max_in_xy = velocity_magnitude
          End Do
          ! 1.3. compute the average viscosity in the xy plane
          average_viscosity = 0.0d0
          Do K = 1, NFEN
             average_viscosity = average_viscosity + Biharmonic_viscosity_modified_Leith ( I, K )
          End Do
          average_viscosity = average_viscosity / NFEN

          ! 1.3.1 compute the element Peclet number
                  element_peclet_number = length_max * velocity_magnitude_max_in_xy / average_viscosity !constant_diffusivity !average_viscosity
          average_element_peclet_number = length_max * velocity_magnitude_max_in_xy / average_viscosity 
          ! 1.4. specify sigma - we want to do this adaptively.


             If ( average_element_peclet_number < low_peclet_threshold ) Then
                factor = low_factor
             Else If ( average_element_peclet_number > high_peclet_threshold ) Then
                factor = high_factor
             Else
                factor = low_factor + ( high_factor - low_factor ) * 
     &           ( ( average_element_peclet_number - low_peclet_threshold ) / ( high_peclet_threshold - low_peclet_threshold ) )**0.5
             End If
   
          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0

!----------------------------------------------------------------------------------------
! 5,6-SURA6 poly
!          If ( element_peclet_number >= 20. ) Then
!             factor = 0.19d0
!          Else If ( (element_peclet_number < 20.).And.(element_peclet_number >= 10.) ) Then
!             factor = 0.18d0
!          Else If ( (element_peclet_number < 10.).And.(element_peclet_number >= 5.) ) Then
!             factor = 0.17d0
!          Else If ( (element_peclet_number < 5.).And.(element_peclet_number >= 1.) ) Then
!             factor = 0.16d0
!          Else
!             factor = 0.15d0
!          End If
         ! arash: this may not be necessary
!          If ( velocity_magnitude_max_in_xy > 1.80d0 ) factor = 0.19d0
!----------------------------------------------------------------------------------------


          sigma_Gauss = length_max * factor


          ! 2. compute Gaussian filter weights
          IF ( ALLOCATED(Weights) ) DEALLOCATE( Weights )
          ALLOCATE( Weights( NNeigh ( I ) ) )
          Weights = 0.0d0
          Do J = 1, NNeigh ( I ) ! also includes node I
             Neighbor_Node = NeiTab ( I , J )
             numerator = ( X(Neighbor_Node) - X(I) )**2 + ( Y(Neighbor_Node) - Y(I) )**2 
             denominator = 2.0d0 * sigma_Gauss**2
             Weights (J) = exp ( -1.0 * numerator / denominator )
          End Do
          ! 2.1 normalize weights
          sum_Weights = 0.0d0
          Do J = 1, NNeigh ( I )
             sum_Weights = sum_Weights + Weights (J)
          End Do
          If ( sum_Weights < 1.0d-6 ) Then
             Write(*,*) 'Possible error in Gaussian filter - error code = 1'
             call terminate(myproc)
          End If
         Weights (:) = Weights (:) / sum_Weights
         ! 2.3 debug
         !Write(*,'(I5,11f6.3)') NNeigh ( I ), weights(1:NNeigh ( I ))
         !Pause
         ! 2.4. compute elevation corresponding to this node
         Hs = DP(I) + IFNLFA * Eta2(I)  !Total depth at current timestep - 'center' node
         ! arash: not sure if the next line is necessary, but probably not. we aren't using it in the bpg routine either.
         ! IF( Hs .LT. H0 ) Hs = H0
         HsOAMB = Hs / AMB

         ! 3. apply the Gaussian filter
          Do J = 1, NNeigh ( I )
             Neighbor_Node = NeiTab ( I , J )
             inarray_vertical (:) = quantity( Neighbor_Node , : )
             Do K = 1, NFEN
                ! 3.0. obtain the z and sigma value corresponding to the 'center' node
                Zk = HsOAMB * ( Sigma(K) - B ) - DP(I)
                Hs_Neighbor_Node = DP(Neighbor_Node) + IFNLFA * Eta2(Neighbor_Node)
                SigmaNN = B + AMB * ( Zk + DP(Neighbor_Node) ) / Hs_Neighbor_Node
                ! 3.1. linearly interpolate surrounding sal/temp on the same horizontal plane
                Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, interpolated_value, K )
                quantity_smooth (I,K) = quantity_smooth (I,K) + Weights (J) * interpolated_value
             End Do
          End Do
       End Do

       ! transfer smooth bathymetry to the original array
       quantity (:,:) = quantity_smooth (:,:)

       End Subroutine

C*************************************************************************
C     Subroutine to interpolate sal/temp to a specified
C     sigma value (SigmaNN) given an initial guess of which sigma
C     level is closest to the specified value.
C     This subroutine is almost identical to SUBROUTINE ZSURFBUOY, written by R.L.
C
! please pay particular attention to what we are trying to interpolate, and for what purpose.
! this subroutine is used to provide sal/temp for filtering (smoothing). Therefore, we oly use linear
! interpolaton. Moreover, if we run into air, or, into the ground, we use the sal/temp value at the surface, or bottom.
! in other words, we use our engineering judgment in the smoothing procedure as well.
!
C                                    a.f.  09/12/15
! af 01/04/2016:
!I use this subroutine to compute the biharmonic viscosity operator and modified Leith viscosity as well.
!
C*************************************************************************

      Subroutine linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, interpolated_value, J )

      USE GLOBAL_3DVS
      IMPLICIT NONE
      real(8) :: interpolated_value
      REAL(8) :: SigmaNN     !Sigma value of a neighbor node

      INTEGER NN
      INTEGER J
!      INTEGER :: LBelo = -999
!      INTEGER :: LAbov = -999
      INTEGER LBelo
      INTEGER LAbov
      INTEGER LTry
      INTEGER IDiag
      REAL(8) SigBelo
      REAL(8) SigAbov
      REAL(8) SigTry
! af: variables for cubic interpolation
!      Integer  :: k1, k2, k3, k4                 ! nodes used for cubic polynomial
!      REAL(8) :: f1, f2, f3, f4                 ! shape functions
!      REAL(8) :: top, bot
!      Integer  :: loop_situation = 0

      real(8) :: inarray_vertical
      Dimension  inarray_vertical ( NFEN )

      call setMessageSource("linear_vertical_interpolation")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IDiag=0

      IF(SigmaNN.LE.b) THEN !at bottom or lower, then use bottom
         LBelo=1
         interpolated_value = inarray_vertical(LBelo)
         SigBelo=b
         SigAbov=b
         GOTO 100
      ENDIF

      IF(SigmaNN.GE.a) THEN ! surface or above, use surface
         LAbov=NFEN
         interpolated_value=inarray_vertical(LAbov)
         SigBelo=a
         SigAbov=a
         GOTO 100
      ENDIF
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--


      LTry=J                    !start search for SIGABOV and SIGBELO
      SigTry=Sigma(LTry)
      IF(SigmaNN.GT.SigTry) THEN !too low
!         loop_situation = 1
         SigBelo=SigTry         !SIGBELO may = SIGTRY
         LBelo=LTry
         LTry=LTry+1            !look at next level higher
 90      SigTry=Sigma(LTry)
         IF(SigmaNN.GT.SigTry) THEN !still too low
            SigBelo=SigTry
            LBelo=LTry
            LTry=LTry+1
            GOTO 90
         ENDIF
         SigAbov=SigTry         !found upper bracketing sigma
         LAbov=LTry
         GOTO 99                !go interpolate
      ENDIF

      IF(SigmaNN.LE.SigTry) THEN !to high
!         loop_situation = 2
         SigAbov=SigTry         !SIGABOV may = SIGTRY
         LAbov=LTry
         LTry=LTry-1            !look at next level lower
 91      SigTry=Sigma(LTry)
         IF(SigmaNN.LE.SigTry) THEN !still too high
            SigAbov=SigTry
            LAbov=LTry
            LTry=LTry-1
            GOTO 91
         ENDIF
         SigBelo=SigTry         !found lower bracketing sigma
         LBelo=LTry
      ENDIF
! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
99    Continue

! use a cubic polynomial interpolation, instead of the above (99) linear interpolation
! 1- decide which nodes to pick: k1:lowest, k2, k3, k4:top-most
!      If ( LAbov == NFEN ) Then
!         k4 = NFEN
!         k3 = NFEN - 1
!         k2 = NFEN - 2
!         k1 = NFEN - 3
!      Else If ( LBelo == 1 ) Then
!         k4 = 4
!         k3 = 3
!         k2 = 2
!         k1 = 1
!      Else
!         k4 = LAbov + 1
!         k3 = LAbov
!         k2 = LBelo
!         k1 = LBelo - 1
!      End If

! 2- compute the shape functions
!        f1 = ((SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k1)-Sigma(k2)) * (Sigma(k1)-Sigma(k3)) * (Sigma(k1)-Sigma(k4)))

!        f2 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k2)-Sigma(k1)) * (Sigma(k2)-Sigma(k3)) * (Sigma(k2)-Sigma(k4)))

!        f3 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k4)))
!     &     / ((Sigma(k3)-Sigma(k1)) * (Sigma(k3)-Sigma(k2)) * (Sigma(k3)-Sigma(k4)))

!        f4 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)))
!     &     / ((Sigma(k4)-Sigma(k1)) * (Sigma(k4)-Sigma(k2)) * (Sigma(k4)-Sigma(k3)))

! 3- interpolate
!      interpolated_value = f1 * inarray_vertical(k1) + f2 * inarray_vertical(k2) 
!     &          + f3 * inarray_vertical(k3) + f4 * inarray_vertical(k4)

! linear interpolation should be enough. the goal is to smooth the velocity
      interpolated_value = (inarray_vertical(LAbov) - inarray_vertical(LBelo))
     &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) + inarray_vertical(LBelo)

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
 100  CONTINUE

      IF(IDiag.EQ.2) THEN
         WRITE(2,*) '******** linear_vertical_interpolation **********'
         WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov',
     &        '      interpolated_value'
         WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,interpolated_value
 777     FORMAT(I7,I5,3(F10.3),E14.5)
      ENDIF

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE linear_vertical_interpolation_real
C**********************************************************************

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   W R I T E   2 D
C--------------------------------------------------------------------
C     jgf49.49.01: Writes a real(8) array with the specified length
C     to a binary hotstart file; takes the file counter as an
C     input argument and returns the new value of the file counter
C     at the end of the write.
C--------------------------------------------------------------------
      SUBROUTINE binaryWrite2D(array, length, lun, counter)
      USE mod_logging, ONLY : setMessageSource, unsetMessageSource, DEBUG,
     &                    allMessage
      IMPLICIT NONE
      REAL(8), intent(in), dimension(length) :: array ! data to write
      INTEGER, intent(in) :: length  ! the number of values to write
      INTEGER, intent(in) :: lun     ! fortran logical unit number to write to
      INTEGER, intent(inout) :: counter ! i/o position in the file
C
      INTEGER i      ! array index
C
      call setMessageSource("binaryWrite2D")
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       DO i=1, length
          WRITE(lun,REC=counter) array(i)
          counter = counter + 1
       END DO
C
#if defined(VSMY_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryWrite2D
C--------------------------------------------------------------------

