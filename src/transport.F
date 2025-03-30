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
      SUBROUTINE TRANS_3D(inarray,diffhoriz,diffvert,outarray,field,TimeLoc)

C     Include the variables needed from the GLOBAL subroutine. 
C         NODELE, NODECODE, NOFF : used in Casey's wet/dry
C                                  additions (left in for additions later) 
C         NBD, NETA : used in the weak implementation of the 
C                     boundary conditions on the ocean boundary.
C
C  kmd48.33bc add variables from global for boundary conditions and
C             removed ScreenUnit
      USE GLOBAL, ONLY: NODECODE, NOFF,  
     &                  HFLUX, sponge, TTBCTIME1, TTBCTIME2,TTBCRATIO,
     &                  q_heat1, q_heat2, TTBCTIMEINC, BCFLAG_TEMP 
Casey 140701: Avoid zero depths.
     &                 ,H0
      USE GLOBAL_3DVS
      USE MESH, ONLY : NP, X, Y, DP, NM, AREAS, SFAC, LBARRAY_POINTER,
     &                  NODELE, NEITABELE, NEITAB, NNEIGH
      USE BOUNDARIES, ONLY : NOPE, NETA, NVDLL, NBDV, NBD, LBCODEI,
     &               BndBCRiver, bcriversec, BCRNVELL, BCRNBVV
Corbitt 120328: Apply Advection Locally
      USE NodalAttributes, ONLY : LoadAdvectionState, AdvectionState
C  kmd48.33bc add implicit none if not in parallel
#ifdef CMPI
      USE MESSENGER
      IMPLICIT NONE
      REAL(8) :: DUMV1(1),DUMV2(1)
#else
      IMPLICIT NONE
#endif

C     Local variables added in from Casey's wet/dry additions
C
      INTEGER :: NCELE
      INTEGER :: TEMPNCELE
      INTEGER :: TEMPSTOP

C     Additional local variables that are passed into or out of the subroutine
      REAL(8), DIMENSION(MNP,MNFEN), INTENT(IN) :: inarray
      REAL(8), DIMENSION(MNP,MNFEN), INTENT(OUT) :: outarray
      REAL(8), INTENT(IN) :: diffhoriz, diffvert, TimeLoc
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
      REAL(8) :: WSigma(MNFEN) !"sigma" vertical velocity
      REAL(8) :: Zk            !z depth of any node k in the vertical
      REAL(8) :: DelSig        ! sigma(k+1)-sigma(k)
      REAL(8) :: DelSigO2      !(sigma(k)-sigma(k-1))/2
      REAL(8) :: SigmaMAOAMB   !(sigma(k)-A)/(a-b)
      REAL(8) :: SigmaMBOAMB   !(sigma(k)-B)/(a-b)
      REAL(8) :: SigAvgMAOAMB  !((sigma(k)+sigma(k-1))/2.d0 - A)/AMB
      REAL(8) :: SigmaNN       !Sigma value of a neighbor node
      REAL(8) :: EtaN1,EtaN2,EtaN3,EtaNFirst !nodal values of NolIFA(Eta1+Eta2)/2
      REAL(8) :: hN1,hN2,hN3,hNFirst !nodal values of DP
      REAL(8) :: Un,Vn         !real,imaginary components of qn
      REAL(8) :: RCLtrans       

      REAL(8) :: SFacAvg  !  kmd48.33bc add in spherical factors


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
      REAL(8) :: cpwater=4000.0d0          ! variable utilized in the top boundary condition for temperature
      REAL(8) :: saltbcskp1,saltbcsk        ! top boundary conditions variables for salinity (later implementation) 

      REAL(8) :: Frtrans(MNFEN)     !right side forcing vector
      REAL(8) :: Mkm1trans(MNFEN)   !1st column (k-1) in left side compact storage matrix
      REAL(8) :: Mktrans(MNFEN)     !2nd column (k) in left side compact storage matrix
      REAL(8) :: Mkp1trans(MNFEN)   !3rd column (k+1) in left side cpmpact storage matrix
      REAL(8) :: LAdvectrans(MNFEN) !lateral advection term in transport eqn
      REAL(8) :: LDiffusion(MNFEN) !lateral diffusion term in transport eqn
      REAL(8) :: VAdvectrans(MNFEN) !vertical advection term in transport eqn
      REAL(8) :: VDiffusion(MNFEN) !vertical diffusion term in transport eqn
      COMPLEX(8) :: qn,qN1,qN2,qN3,qNFirst !nodal values of q
      REAL(8) :: transN1,transN2,transN3,transNFirst !nodal values of transport equation
      REAL(8) :: UnDtransDX,VnDtransDY  !derivatives used in lateral advection
      REAL(8) :: UnDtransDX2A,VnDtransDY2A !UnDqDX,UnDqDY)*2*Element Area
      REAL(8) :: DtransDXDPhiDX2A,DtransDYDPhiDY2A !derivatives used in lateral diffusion calc.
      REAL(8) :: DtransDSigmakm1,DtransDSigmakp1 !vertical deriv. of q from k-1,k and k,k+1
      COMPLEX(8) :: QSave(MNP,MNFEN) !Used in the testing of the transport equation (remove if works correctly)
      REAL(8) :: WSave(MNP,MNFEN) !Used in the testing of the transport equation (remove if works correctly)

C     Variables for the weak implementation of the open ocean boundary condtion
      REAL(8) :: WBC(MNP,MNFEN) 
      REAL(8) :: RESWBC(MNP,MNFEN)
      REAL(8) :: NORM1(3)
      REAL(8) :: TOTRES, RESWBCN
      REAL(8) :: VELDOTNORM, AVERESN, AVEUN, AVEVN, AVEWN
      REAL(8) :: Hsp1OAMBN1, Hsp1OAMBN2
      REAL(8) :: Hsp1N1, Hsp1N2
      REAL(8) :: XDIST1, YDIST1, ZDIST1, XDIST2, YDIST2, ZDIST2
      REAL(8) :: ZVN1, ZVN2, ZVN3, ZVN4
      REAL(8) :: XSQUARE, YSQUARE, LENGTHEDGE, ref 
      INTEGER :: NM1, NM2, VM1, VM2
      INTEGER :: NBC, BCSEC             !output loop counter

C  kmd48.33bc added in variables used in the heat flux boundary condition
      INTEGER :: NumofNodes, NOD, J
      REAL(8), ALLOCATABLE :: TMP(:,:)

!----------------------------------------------------------------------------------------------------------------------------------------------------
! af: new variables to compute the horizontal diffusion at the xyz coordinate system, similar to the computation of BPG.
! to this end, we interpolate c (quantity to be transported) at the same z levels, and compute horizontal spatial derivatives.
      REAL(8) :: transN1_z, transN2_z, transN3_z, transNFirst_z  !nodal values of transport equation (at the same z)
      Real(8) :: diff_tol = 1.0d-5
      Real(8) :: minimum_diffusion
      Real(8) :: max_diffusion
      REAL(8) :: inarray_vertical(NFEN) 
      REAL(8) :: diffusion_ratio = 0.20d0

! af 01/04/2016: Biharmonic viscosity for transport
      Real(8) :: DSalOrTempDXDPhiDX2A, DSalOrTempDYDPhiDY2A
      Real(8) :: Biharmonic_DSalOrTempDXDPhiDX2A, Biharmonic_DSalOrTempDYDPhiDY2A
!      Real(8) :: d_dx_SalOrTemp_element, d_dy_SalOrTemp_element
      Real(8) :: SalOrTempN_z, SalOrTempN1_z, SalOrTempN2_z, SalOrTempN3_z, SalOrTempNFirst_z
      Real(8) :: Biharmonic_SalOrTempN_z, Biharmonic_SalOrTempN1_z, Biharmonic_SalOrTempN2_z, Biharmonic_SalOrTempN3_z
      Real(8) :: Biharmonic_SalOrTempNFirst_z
      Real(8) :: Hs, HsOAMB, HsN2
      Real(8) :: Biharmonic_viscosity
      Real(8) :: Biharmonic_LDiffusion(MNFEN) !lateral diffusion term in transport eqn - biharmonic diffusion
!Casey 220207 from Rosemary : For diffusion term in wet/dry interface
      Real(8) :: LDiffusion_WetDryInterface(MNFEN)

C*************************************************************************************
C     Compute 3D transport equation

! -------------------------------------------------------------------------------------------------
! arash 01/04/2016: biharmonic diffusion for the transport equation
! compute the inner laplacian of the biharmonic operator.
! Del Del q = Del Biharmonic_auxiliary_var
! here, q = temp/sal. in the momentum equation, q = u + i v.
      Biharmonic_auxiliary_var_SalOrTemp (:,:) = 0.0d0
!      d_dx_SalOrTemp_nodes (:,:) = 0.0d0    ! spatial derivarievs to compute modified Leith viscosity
!      d_dy_SalOrTemp_nodes (:,:) = 0.0d0

      DO NH = 1, NP                      !loop over horizontal nodes
         Hs = DP(NH) + IFNLFA * Eta1(NH) !Total depth at current timestep - 'center' node
!Casey 220207 from Rosemary
         IF(Hs.LT.H0) Hs = H0
         HsOAMB = Hs / AMB

#ifndef CASEY

!Casey 220523 : Wouldn't this be easier?

         DO k=1,NFEN
            ! Start from zero
            TotalArea2=0.d0
            DSalOrTempDXDPhiDX2A=0.d0
            DSalOrTempDYDPhiDY2A=0.d0
            DO N=1,SIZE(NeiTabEle,2)
               ! Find neighbor element
               NEle=NeiTabEle(NH,N)
               ! Skip blank entries in NeiTabEle
               IF(NEle.EQ.0) CYCLE
               ! Find vertices in counterclockwise order
               IF(NM(NEle,1).EQ.NH)THEN
                  N1=NM(NEle,1)
                  N2=NM(NEle,2)
                  N3=NM(NEle,3)
               ELSEIF(NM(NEle,2).EQ.NH)THEN
                  N1=NM(NEle,2)
                  N2=NM(NEle,3)
                  N3=NM(NEle,1)
               ELSEIF(NM(NEle,3).EQ.NH)THEN
                  N1=NM(NEle,3)
                  N2=NM(NEle,1)
                  N3=NM(NEle,2)
               ELSE
                  PRINT*, "Aack! ",NH,NEle,NM(Nele,:)
               ENDIF
               ! Skip dry elements
               NCEle=NODECODE(N1)*NODECODE(N2)
     &              *NODECODE(N3)*NOFF(NEle)
               IF(NCEle.EQ.0) CYCLE
               ! For vertex of interest
               Zk = HsOAMB * ( Sigma(k) - B ) - DP(N1)
               SalOrTempN1_z = inarray(N1,k)
               ! For first neighbor vertex 
               HsN2 = DP(N2) + IFNLFA * Eta1(N2)
               IF(HsN2.LT.H0) HsN2 = H0
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2
               inarray_vertical (:) = inarray(N2,:)
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real ( 
     &                 inarray_vertical,
     &                 SigmaNN, SalOrTempN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical, 
     &                 SigmaNN, SalOrTempN2_z, N2, k )
               ENDIF
               ! For second neighbor vertex 
               HsN2 = DP(N3) + IFNLFA * Eta1(N3)
               IF(HsN2.LT.H0) HsN2 = H0
               SigmaNN = B + AMB * ( Zk + DP(N3) ) / HsN2
               inarray_vertical (:) = inarray(N3,:)
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real ( 
     &                 inarray_vertical,
     &                 SigmaNN, SalOrTempN3_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                 SigmaNN, SalOrTempN3_z, N2, k )
               ENDIF
               ! Compute gradients
               IF(SalOrTempN2_z.GT.-990.D0.AND.
     &            SalOrTempN3_z.GT.-990.D0)THEN
                  TotalArea2 = TotalArea2 + Areas(NEle)
                  SFacAvg = (SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  a1 = X(N3)-X(N2)
                  a2 = X(N1)-X(N3)
                  a3 = X(N2)-X(N1)
                  b1 =(Y(N2)-Y(N3))*SFacAvg
                  b2 =(Y(N3)-Y(N1))*SFacAvg
                  b3 =(Y(N1)-Y(N2))*SFacAvg
                  DSalOrTempDXDPhiDX2A = DSalOrTempDXDPhiDX2A
     &               + ( SalOrTempN1_z*b1 + SalOrTempN2_z*b2
     &                 + SalOrTempN3_z*b3 )*b1/Areas(NEle)
                  DSalOrTempDYDPhiDY2A = DSalOrTempDYDPhiDY2A
     &               + ( SalOrTempN1_z*a1 + SalOrTempN2_z*a2
     &                 + SalOrTempN3_z*a3 )*a1/Areas(NEle)
               ENDIF
            ENDDO
            IF(TotalArea2.GT.0.0001D0)THEN
               Biharmonic_auxiliary_var_SalOrTemp(NH,k) = -3.d0 *
     &            (DSalOrTempDXDPhiDX2A+DSalOrTempDYDPhiDY2A)/TotalArea2
            ENDIF
         ENDDO

#else

!Casey 220523 : End of code additions, let's comment the stuff below.

         ! begin the vertical loop
         DO k = 1, NFEN
            DSalOrTempDXDPhiDX2A = 0.d0    ! variables for computing the biharmonic operator
            DSalOrTempDYDPhiDY2A = 0.d0
            TotalArea2   = 0.d0
!            d_dx_SalOrTemp_element = 0.0d0 ! to compute the viscosity constant
!            d_dy_SalOrTemp_element = 0.0d0
            Zk = HsOAMB * ( Sigma(k) - B ) - DP(NH)

            N1  = NH               !node 1 is always the central node
            SalOrTempN1_z = inarray(N1,k)

            N2  = NEITAB(NH,2)     !operate on 1st neighbor
            HsN2 = DP(N2) + IFNLFA * Eta1(N2)
!Casey 220207 from Rosemary
            IF(HsN2.LT.H0) HsN2 = H0
            SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2 !equivalent sigma value at neighbor
            inarray_vertical (:) = inarray(N2,:)
!Casey 220227 from Arash/Rosemary
!Casey 221103 : User control over cubic interpolation
            IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
               CALL linear_vertical_interpolation_real ( 
     &             inarray_vertical, SigmaNN, SalOrTempN2_z, k )
            ELSE ! Default
               CALL cubic_vertical_interpolation( inarray_vertical,
     &             SigmaNN,SalOrTempN2_z,N2,k)
            ENDIF

            NNFirst = N2           !save these values until end
            SalOrTempNFirst_z = SalOrTempN2_z

            DO N = 3, NNEIGH(NH)   !operate on rest of neighbors
               N3 = N2            !shift previously computed values
               SalOrTempN3_z = SalOrTempN2_z

               N2 = NEITAB(NH,N)  !select new neighbor to work on
               HsN2 = DP(N2) + IFNLFA * Eta1(N2)
!Casey 220207 from Rosemary
               IF(HsN2.LT.H0) HsN2 = H0
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2 !equivalent sigma value at neighbor
               inarray_vertical (:) = inarray(N2,:)
!Casey 220227 from Arash/Rosemary
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real ( 
     &                inarray_vertical, SigmaNN, SalOrTempN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                SigmaNN,SalOrTempN2_z,N2,k)
               ENDIF

               NEle= NeiTabEle(NH,N-2) !element # defined by nodes N1,N2,N3
               IF (NEle.eq.0) THEN
                  CYCLE
               ENDIF
               IF((NEle.NE.0)) THEN    !if element is active, compute velocity grads
                  NCELE = NODECODE(NH) * NODECODE(N2) *
     &                 NODECODE(N3) * NOFF(NELE)
!Casey 220227 from Arash/Rosemary
!                 IF (NCELE.NE.0) THEN !if element exists, compute terms
                  IF ((NCELE.NE.0).and.(SalOrTempN2_z.GT.-998.d0).and.
     &                    (SalOrTempN3_z.GT.-998.d0))THEN !if element exists, compute terms
                     TotalArea2 = TotalArea2 + Areas(NEle) !accumulate 2X total areas to complete calc.
                     SFacAvg = (SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                     a1 = X(N3)-X(N2)
                     a2 = X(N1)-X(N3)
                     a3 = X(N2)-X(N1)
                     b1 = (Y(N2)-Y(N3))*SFacAvg
                     b2 = (Y(N3)-Y(N1))*SFacAvg
                     b3 = (Y(N1)-Y(N2))*SFacAvg
! z-level:
                     DSalOrTempDXDPhiDX2A = DSalOrTempDXDPhiDX2A
     &                 + (SalOrTempN1_z*b1+SalOrTempN2_z*b2+SalOrTempN3_z*b3)*b1/Areas(NEle)
                     DSalOrTempDYDPhiDY2A=DSalOrTempDYDPhiDY2A
     &                 +(SalOrTempN1_z*a1+SalOrTempN2_z*a2+SalOrTempN3_z*a3)*a1/Areas(NEle)
                     ! variables for modified Leith viscosity parameter
!                     d_dx_SalOrTemp_element = d_dx_SalOrTemp_element + SalOrTempN1_z*b1+SalOrTempN2_z*b2+SalOrTempN3_z*b3
!                     d_dy_SalOrTemp_element = d_dy_SalOrTemp_element + SalOrTempN1_z*a1+SalOrTempN2_Z*a2+SalOrTempN3_z*a3
                  ENDIF
               ENDIF
            ENDDO
            N3  = N2               !wrap back to beginning to get final contribution
            SalOrTempN3_z = SalOrTempN2_z

            N2  = NNFIRST
            SalOrTempN2_z = SalOrTempNFirst_z

            NEle=NeiTabEle(NH,NNeigh(NH)-1)
            IF((NEle.NE.0)) THEN  !if element is active, compute velocity grads
               NCELE = NODECODE(NH) * NODECODE(N2) *
     &              NODECODE(N3) * NOFF(NELE)
!Casey 220227 from Arash/Rosemary
!              IF ((NCELE.NE.0)) THEN  !if element is active, compute velocity grads
               IF ((NCELE.NE.0).and.(SalOrTempN2_z.GT.-998.d0).and.
     &                   (SalOrTempN3_z.GT.-998.d0)) THEN  !if element is active, compute velocity grads
                  TotalArea2 = TotalArea2 + Areas(NEle) !accumulate 2X total areas to complete calc.
                  SFacAvg = (SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
                  a1 = X(N3)-X(N2)
                  a2 = X(N1)-X(N3)
                  a3 = X(N2)-X(N1)
                  b1 = (Y(N2)-Y(N3))*SFacAvg
                  b2 = (Y(N3)-Y(N1))*SFacAvg
                  b3 = (Y(N1)-Y(N2))*SFacAvg
! z-level:
                     DSalOrTempDXDPhiDX2A = DSalOrTempDXDPhiDX2A
     &                 + (SalOrTempN1_z*b1+SalOrTempN2_z*b2+SalOrTempN3_z*b3)*b1/Areas(NEle)
                     DSalOrTempDYDPhiDY2A=DSalOrTempDYDPhiDY2A
     &                 + (SalOrTempN1_z*a1+SalOrTempN2_z*a2+SalOrTempN3_z*a3)*a1/Areas(NEle)
                     ! variables for modified Leith viscosity parameter
!                     d_dx_SalOrTemp_element = d_dx_SalOrTemp_element + SalOrTempN1_z*b1+SalOrTempN2_z*b2+SalOrTempN3_z*b3
!                     d_dy_SalOrTemp_element = d_dy_SalOrTemp_element + SalOrTempN1_z*a1+SalOrTempN2_Z*a2+SalOrTempN3_z*a3
               ENDIF
            ENDIF

            IF(TotalArea2.EQ.0.) THEN
               ! do nothing
            ELSE
               Biharmonic_auxiliary_var_SalOrTemp (NH,k) = -3.d0 * (DSalOrTempDXDPhiDX2A+DSalOrTempDYDPhiDY2A) / TotalArea2
               ! variables for modified Leith viscosity parameter
!               d_dx_SalOrTemp_nodes (NH,k) = d_dx_SalOrTemp_element / TotalArea2
!               d_dy_SalOrTemp_nodes (NH,k) = d_dy_SalOrTemp_element / TotalArea2
            ENDIF

         END DO ! end vertical loop

!Casey 220523 : End of commenting.

#endif

      End Do ! end horizontal loop

! end computing Biharmonic_auxiliary_var(:,:)
! -------------------------------------------------------------------------------------------------

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
Casey 140701: Avoid zero depths.
         Hsp1N1=DP(NM1)+IFNLFA*ETA2(NM1)
         IF(Hsp1N1.LT.H0) Hsp1N1=H0
         Hsp1N2=DP(NM2)+IFNLFA*ETA2(NM2)
         IF(Hsp1N2.LT.H0) Hsp1N2=H0
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
           IF (TimeLoc.GT.TTBCTIME2) THEN
              TTBCTIME1=TTBCTIME2
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=q_heat2(NumofNodes)
                 READ(38,*) NOD, q_heat2(NumofNodes)
              END DO
           END IF ! ends time if loop
           TTBCRATIO=(TimeLoc-TTBCTIME1)/TTBCTIMEINC
           DO NumofNodes=1,NP
              qsurf(NumofNodes)=qsurfkp1(NumofNodes)
              HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
!--------------------------------------------------------------------------------------------------
! arash 01/08/2016: limit the surface heat flux if it is out of range - this is a reasonable range according to the literature.
!               If ( HFLUX(NumofNodes) > 350.0d0 ) HFLUX(NumofNodes) = 350.0d0
!               If ( HFLUX(NumofNodes) <-350.0d0 ) HFLUX(NumofNodes) =-350.0d0
!                HFLUX = - HFLUX
!--------------------------------------------------------------------------------------------------
           END DO
         ELSE IF (BCFLAG_TEMP.EQ.2) THEN
            ALLOCATE (TMP(NP,6))
            IF (TimeLoc.GT.TTBCTIME2) THEN
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
            TTBCRATIO=(TimeLoc-TTBCTIME1)/TTBCTIMEINC
            DO  NumofNodes=1,NP
               qsurf(NumofNodes)=qsurfkp1(NumofNodes)
               HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
!--------------------------------------------------------------------------------------------------
! arash 01/08/2016: limit the surface heat flux if it is out of range - this is a reasonable range according to the literature.
!               If ( HFLUX(NumofNodes) > 350.0d0 ) HFLUX(NumofNodes) = 350.0d0
!               If ( HFLUX(NumofNodes) <-350.0d0 ) HFLUX(NumofNodes) =-350.0d0
!                HFLUX = - HFLUX
!--------------------------------------------------------------------------------------------------
            END DO
            DEALLOCATE(TMP)
         ELSE IF (BCFLAG_TEMP.EQ.3) THEN
            ALLOCATE (TMP(NP,4))
            IF (TimeLoc.GT.TTBCTIME2) THEN
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
            TTBCRATIO=(TimeLoc-TTBCTIME1)/TTBCTIMEINC
            DO NumofNodes=1,NP
               qsurf(NumofNodes)=qsurfkp1(NumofNodes)
               HFLUX(NumofNodes)=q_heat1(NumofNodes)+TTBCRATIO*
     &            (q_heat2(NumofNodes)-q_heat1(NumofNodes))
!--------------------------------------------------------------------------------------------------
! arash 01/08/2016: limit the surface heat flux if it is out of range - this is a reasonable range according to the literature.
!               If ( HFLUX(NumofNodes) > 350.0d0 ) HFLUX(NumofNodes) = 350.0d0
!               If ( HFLUX(NumofNodes) <-350.0d0 ) HFLUX(NumofNodes) =-350.0d0
!                HFLUX = - HFLUX
!--------------------------------------------------------------------------------------------------
            END DO
            DEALLOCATE(TMP)
         END IF ! ends the BCFLAG_TEMP flag if loop
      END IF  ! ends the field if loop

!--------------------------------------------------------------------------------------------------
! Surface Heat fllux minus sign
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------


C   kmd_rivers - add information for the weak boundary condition with
C              rivers for the baroclinic code.

C     Determine the weak boundary conditions that will be used in the
C     assembly process later. Note that the boundary conditions will
C     be utilized in both the salinity and temperature fields.
      IF (BNDBCRiver) THEN
        DO BCSEC=1, bcriversec
          DO NBC=1,BCRNVELL(BCSEC)-1
C     First obtain the two nodes that with be top two nodes of
C     the rectangle and define two vertical nodes for the depths.
             NM1=BCRNBVV(BCSEC,NBC)
             NM2=BCRNBVV(BCSEC,NBC+1)
             VM1=NFEN
             VM2=NFEN-1
C     Compute the depths at the vertical nodes associated with
C     these two nodes.
Casey 140701: Avoid zero depths.
             Hsp1N1=DP(NM1)+IFNLFA*ETA2(NM1)
             IF(Hsp1N1.LT.H0) Hsp1N1=H0
             Hsp1N2=DP(NM2)+IFNLFA*ETA2(NM2)
             IF(Hsp1N2.LT.H0) Hsp1N2=H0
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
     &            YDIST2,ZDIST2,NORM1)

C     Now develop the average velocity values for the dot product
C     to be determined later.

             DO k=1,NFEN
                IF (k.NE.NFEN) THEN
                   AVEUN=0.25d0*(REAL(Q(NM1,k))+REAL(Q(NM2,k))+
     &                  REAL(Q(NM2,k+1))+REAL(Q(NM1,k+1)))
                   AVEVN=0.25d0*(AIMAG(Q(NM1,k))+AIMAG(Q(NM2,k))+
     &                  AIMAG(Q(NM2,k+1))+AIMAG(Q(NM1,k+1)))
                   AVEWN=0.25d0*(WZ(NM1,k)+WZ(NM2,k)+
     &                  WZ(NM2,k+1)+WZ(NM1,k+1))
                ELSE IF (k.EQ.NFEN) THEN
                   AVEUN=0.25d0*(REAL(Q(NM1,k-1))+REAL(Q(NM2,k-1))+
     &                  REAL(Q(NM2,k))+REAL(Q(NM1,k)))
                   AVEVN=0.25d0*(AIMAG(Q(NM1,k-1))+AIMAG(Q(NM2,k-1))+
     &                  AIMAG(Q(NM2,k))+AIMAG(Q(NM1,k)))
                   AVEWN=0.25d0*(WZ(NM1,k-1)+WZ(NM2,k-1)+
     &                  WZ(NM2,k)+WZ(NM1,k))
                END IF
                VELDOTNORM=AVEUN*NORM1(1)+AVEVN*NORM1(2)+
     &               AVEWN*NORM1(3)
                XSQUARE=(X(NM1)-X(NM2))*(X(NM1)-X(NM2))
                YSQUARE=(Y(NM1)-Y(NM2))*(Y(NM1)-Y(NM2))
                LENGTHEDGE=SQRT(XSQUARE+YSQUARE)
                IF (field.EQ.2) THEN
                   IF (k.NE.NFEN) THEN
                       ref=0.25d0*(BCRivSal(NBC,k)+BCRivSal(NBC,k+1)+
     &                    BCRivSal(NBC+1,k)+BCRivSal(NBC+1,k+1))
                   ELSE IF (k.EQ.NFEN) THEN
                       ref=0.25d0*(BCRivSal(NBC,k-1)+BCRivSal(NBC,k)+
     &                    BCRivSal(NBC+1,k-1)+BCRivSal(NBC+1,k))
                   END IF
                   !ref=32.d0   ! kmd48.33bc changed to 32 as the salinity reference
                ELSE IF (field.EQ.3) THEN
                   IF (k.NE.NFEN) THEN
                       ref=0.25d0*(BCRivTemp(NBC,k)+BCRivTemp(NBC,k+1)+
     &                    BCRivTemp(NBC+1,k)+BCRivTemp(NBC+1,k+1))
                   ELSE IF (k.EQ.NFEN) THEN
                       ref=0.25d0*(BCRivTemp(NBC,k-1)+BCRivTemp(NBC,k)+
     &                    BCRivTemp(NBC+1,k-1)+BCRivTemp(NBC+1,k))
                   END IF
                   !ref=10.d0
                END IF
                IF (VELDOTNORM.LE.0.d0) THEN
                   IF (k.NE.NFEN) THEN
                      AVERESN=0.25d0*(inarray(NM1,k)+inarray(NM2,k)+
     &                     inarray(NM2,k+1)+inarray(NM1,k+1))
                   ELSE IF (k.EQ.NFEN) THEN
                      AVERESN=0.25d0*(inarray(NM1,k-1)+inarray(NM2,k-1)+
     &                     inarray(NM2,k)+inarray(NM1,k))
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
      END IF
      
C     
C     Loop over each horizontal and vertical nodes to compute the transport concentrations
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
Casey 140829
C        IFNLCT = 1
         
C     Set up some values at the node being worked on
         
Casey 140701: Avoid zero depths.
         Hstrans  = DP(NH)+IFNLFA*Eta1(NH) !Total depth at previous (s) timestep
         IF(Hstrans.LT.H0) Hstrans=H0
         HsOAMBtrans=Hstrans/AMB
         HsHsOAMBAMBtrans=HsOAMBtrans*HsOAMBtrans
         
         Hsp1trans= DP(NH)+IFNLFA*Eta2(NH) !Total depth at present (s+1) timestep
!Casey 220207 from Rosemary
         IF(Hsp1trans.LT.H0) Hsp1trans = H0
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
           IF ((IEVC.EQ.50).or.(IEVC.EQ.51)) THEN 
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
         
#ifndef CASEY

!Casey 220523 : Wouldn't this be easier?

         ! Start from zero
         TotalArea2=0.d0
         DhDX2A=0.d0
         DhDY2A=0.d0
         DEtaDX2A=0.d0
         DEtaDY2A=0.d0
         DhDX=0.d0
         DhDY=0.d0
         DEtaDX=0.d0
         DEtaDY=0.d0
         DO N=1,SIZE(NeiTabEle,2)
            ! Find neighbor element
            NEle=NeiTabEle(NH,N)
            ! Skip blank entries in NeiTabEle
            IF(NEle.EQ.0) CYCLE
            ! Find vertices in counterclockwise order
            IF(NM(NEle,1).EQ.NH)THEN
               N1=NM(NEle,1)
               N2=NM(NEle,2)
               N3=NM(NEle,3)
            ELSEIF(NM(NEle,2).EQ.NH)THEN
               N1=NM(NEle,2)
               N2=NM(NEle,3)
               N3=NM(NEle,1)
            ELSEIF(NM(NEle,3).EQ.NH)THEN
               N1=NM(NEle,3)
               N2=NM(NEle,1)
               N3=NM(NEle,2)
            ELSE
               PRINT*, "Aack! ",NH,NEle,NM(Nele,:)
            ENDIF
            ! Skip dry elements
            NCEle=NODECODE(N1)*NODECODE(N2)
     &           *NODECODE(N3)*NOFF(NEle)
            IF(NCEle.EQ.0) CYCLE
            ! For vertex of interest
            EtaN1=IFNLFA*(Eta1(N1)+Eta2(N1))/2.d0
            hN1=DP(N1)
            IF(hN1.LT.H0) hN1 = H0
            ! For first neighbor vertex
            EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
            hN2=DP(N2)
            IF(hN2.LT.H0) hN2 = H0
            ! For second neighbor vertex
            EtaN3=IFNLFA*(Eta1(N3)+Eta2(N3))/2.d0
            hN3=DP(N3)
            IF(hN3.LT.H0) hN3 = H0
            ! Compute gradients
            TotalArea2=TotalArea2+Areas(NEle)
            SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
            a1= X(N3)-X(N2)
            a2= X(N1)-X(N3)
            a3= X(N2)-X(N1)
            b1=(Y(N2)-Y(N3))*SFacAvg
            b2=(Y(N3)-Y(N1))*SFacAvg
            b3=(Y(N1)-Y(N2))*SFacAvg
            DhDX2A=DhDX2A+(hN1*b1+hN2*b2+hN3*b3)
            DhDY2A=DhDY2A+(hN1*a1+hN2*a2+hN3*a3)
            DEtaDX2A=DEtaDX2A+(EtaN1*b1+EtaN2*b2+EtaN3*b3)
            DEtaDY2A=DEtaDY2A+(EtaN1*a1+EtaN2*a2+EtaN3*a3)
         ENDDO
         IF(TotalArea2.GT.0.0001D0) THEN
            DhDX=DhDX2A/TotalArea2
            DhDY=DhDY2A/TotalArea2
            DEtaDX=DEtaDX2A/TotalArea2
            DEtaDY=DEtaDY2A/TotalArea2
         ENDIF

#else

!Casey 220523 : End of code additions, let's comment the stuff below.

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
!Casey 220207 from Rosemary
         IF(hN1.LT.H0) hN1 = H0

         N2=NeiTab(NH,2)        !operate on 1st neighbor
         EtaN2=IFNLFA*(Eta1(N2)+Eta2(N2))/2.d0
         hN2=DP(N2)
!Casey 220207 from Rosemary
         IF(hN2.LT.H0) hN2 = H0

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
!Casey 220207 from Rosemary
            IF(hN2.LT.H0) hN2 = H0
            NEle=NeiTabEle(NH,N-2) !element # defined by nodes NH,NN2,NN1

            ! jgf49.58: NOFF array lookups fail if NEle comes back 0.
            IF (NEle.eq.0) THEN
               CYCLE
            ENDIF
 
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
         ! jgf49.58: NOFF array lookups fail if NEle comes back 0.
         IF (NEle.ne.0) THEN
            NCELE = NODECODE(NH)*NODECODE(N2)
     &          *NODECODE(N3)*NOFF(NELE)
         ENDIF
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

!Casey 220523 : End of commenting.

#endif

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

!----------------------------------------------------------------------------------------------------------------------------------------------------
! vertical integration
!----------------------------------------------------------------------------------------------------------------------------------------------------
         DO k = 1, NFEN

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

!----------------------------------------------------------------------------------------------------------------------------------------------------
C     Compute lateral advection and lateral diffusion terms - af modified computation of horizontal diffusion terms
!----------------------------------------------------------------------------------------------------------------------------------------------------

!Casey 220523 : Zero these variables.
            LAdvectrans(k) = 0.d0
            LDiffusion(k)  = 0.d0
            Biharmonic_LDiffusion(k) = 0.d0
            LDiffusion_WetDryInterface(k) = 0.d0

            Zk = HsOAMBtrans * ( Sigma(k) - B ) - DP(NH)               ! determine z corresponding to sigma level k

#ifndef CASEY

!Casey 220523 : Wouldn't this be easier?

            ! Start from zero
            TotalArea2 = 0.d0
            UnDtransDX2A=0.d0
            VnDtransDY2A=0.d0
            DtransDXDPhiDX2A=0.d0
            DtransDYDPhiDY2A=0.d0
            Biharmonic_DSalOrTempDXDPhiDX2A = 0.d0
            Biharmonic_DSalOrTempDYDPhiDY2A = 0.d0
            DO N=1,SIZE(NeiTabEle,2)
               ! Find neighbor element
               NEle=NeiTabEle(NH,N)
               ! Skip blank entries in NeiTabEle
               IF(NEle.EQ.0) CYCLE
               ! Find vertices in counterclockwise order
               IF(NM(NEle,1).EQ.NH)THEN
                  N1=NM(NEle,1)
                  N2=NM(NEle,2)
                  N3=NM(NEle,3)
               ELSEIF(NM(NEle,2).EQ.NH)THEN
                  N1=NM(NEle,2)
                  N2=NM(NEle,3)
                  N3=NM(NEle,1)
               ELSEIF(NM(NEle,3).EQ.NH)THEN
                  N1=NM(NEle,3)
                  N2=NM(NEle,1)
                  N3=NM(NEle,2)
               ELSE
                  PRINT*, "Aack! ",NH,NEle,NM(Nele,:)
               ENDIF
               ! Skip dry elements
               NCEle=NODECODE(N1)*NODECODE(N2)
     &              *NODECODE(N3)*NOFF(NEle)
               IF(NCEle.EQ.0) CYCLE
               ! For vertex of interest
               qN1=q(N1,k)
               transN1   = inarray(N1,k)
               transN1_z = transN1
               Biharmonic_SalOrTempN1_z = 
     &            Biharmonic_auxiliary_var_SalOrTemp(N1,k)
               ! For first neighbor vertex
               qN2=q(N2,k)
               transN2=inarray(N2,k)
               HsN2trans = DP(N2) + IFNLFA * Eta1(N2)
               IF(HsN2trans.LT.H0) HsN2trans = H0
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2trans
               inarray_vertical (:) = inarray(N2,:)
!Casey 220207 from Rosemary : Some confusing changes through this
!             section. In her code, the two cubic interpolations are
!             used, not sure why they would be disabled here. Let's try them.
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real( 
     &                inarray_vertical, SigmaNN, transN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                SigmaNN, transN2_z, N2, k )
               ENDIF
               inarray_vertical(:) =
     &            Biharmonic_auxiliary_var_SalOrTemp(N2,:)
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real( 
     &                inarray_vertical,
     &                SigmaNN, Biharmonic_SalOrTempN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                SigmaNN, Biharmonic_SalorTempN2_z, N2, k )
               ENDIF
               ! For second neighbor vertex
               qN3=q(N3,k)
               transN3=inarray(N3,k)
               HsN2trans = DP(N3) + IFNLFA * Eta1(N3)
               IF(HsN2trans.LT.H0) HsN2trans = H0
               SigmaNN = B + AMB * ( Zk + DP(N3) ) / HsN2trans
               inarray_vertical (:) = inarray(N3,:)
!Casey 220207 from Rosemary : Same thing here?
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real( 
     &                inarray_vertical,
     &                SigmaNN, transN3_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                SigmaNN, transN3_z, N3, k )
               ENDIF
               inarray_vertical(:) =
     &            Biharmonic_auxiliary_var_SalOrTemp(N3,:)
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real( 
     &                inarray_vertical,
     &                SigmaNN, Biharmonic_SalOrTempN3_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation( inarray_vertical,
     &                SigmaNN, Biharmonic_SalorTempN3_z, N3, k )
               ENDIF
               ! Compute gradients
               TotalArea2=TotalArea2+Areas(NEle)
               qn=(qN1+qN2+qN3)/3.d0
               Un=REAL(qn)
               Vn=AIMAG(qn)
               SFacAvg=(SFAC(N1)+SFAC(N2)+SFAC(N3))/3.d0
               a1= X(N3)-X(N2)
               a2= X(N1)-X(N3)
               a3= X(N2)-X(N1)
               b1=(Y(N2)-Y(N3))*SFacAvg
               b2=(Y(N3)-Y(N1))*SFacAvg
               b3=(Y(N1)-Y(N2))*SFacAvg
               UnDtransDX2A = UnDtransDX2A 
     &            + Un*( transN1*b1 + transN2*b2 + transN3*b3 )
               VnDtransDY2A = VnDtransDY2A
     &            + Vn*( transN1*a1 + transN2*a2 + transN3*a3 )
               IF(transN2_z.GT.-990.D0.AND.transN3_z.GT.-990.D0)THEN
                  DtransDXDPhiDX2A = DtransDXDPhiDX2A
     &               + ( transN1_z*b1 + transN2_z*b2 + transN3_z*b3 )
     &                 * b1/Areas(NEle)
                  DtransDYDPhiDY2A = DtransDYDPhiDY2A
     &               + ( transN1_z*a1 + transN2_z*a2 + transN3_z*a3 )
     &                 * a1/Areas(NEle)
               ENDIF
               IF(Biharmonic_SalOrTempN2_z.GT.-990.D0.AND.
     &            Biharmonic_SalorTempN3_z.GT.-998.D0)THEN
                  Biharmonic_DSalOrTempDXDPhiDX2A = 
     &               Biharmonic_DSalOrTempDXDPhiDX2A
     &               + ( Biharmonic_SalOrTempN1_z*b1
     &                 + Biharmonic_SalOrTempN2_z*b2
     &                 + Biharmonic_SalOrTempN3_z*b3 )
     &                 * b1/Areas(NEle)
                  Biharmonic_DSalOrTempDYDPhiDY2A =
     &               Biharmonic_DSalOrTempDYDPhiDY2A
     &               + ( Biharmonic_SalOrTempN1_z*a1
     &                 + Biharmonic_SalOrTempN2_z*a2
     &                 + Biharmonic_SalOrTempN3_z*a3 )
     &                 * a1/Areas(NEle)
               ENDIF
            ENDDO
            IF(TotalArea2.GT.0.0001D0)THEN
               LAdvectrans(k)=sponge(NH)*IFNLCT
     &            *(UnDtransDX2A+VnDtransDY2A)/TotalArea2
               LDiffusion(k)=3.d0*diffhoriz
     &            *(DtransDXDPhiDX2A+DtransDYDPhiDY2A)/TotalArea2
!Casey 221103 : User control over Leith viscosity
               IF( .NOT.useLeithViscosity )THEN
                  Biharmonic_viscosity = diffhoriz
               ELSE ! Default
                  Biharmonic_viscosity =
     &               Biharmonic_viscosity_modified_Leith ( NH, k ) * 1.00d0
               ENDIF
               Biharmonic_LDiffusion(k) = -3.d0 * Biharmonic_viscosity *
     &            ( Biharmonic_DSalOrTempDXDPhiDX2A
     &            + Biharmonic_DSalOrTempDYDPhiDY2A )/TotalArea2
            ENDIF

#else

!Casey 220523 : End of code additions, let's comment the stuff below.

            UnDtransDX2A=0.d0
            VnDtransDY2A=0.d0
            DtransDXDPhiDX2A=0.d0
            DtransDYDPhiDY2A=0.d0

! arash 01/05/2016: Biharmonic auxiliary variable derivatives and modified Leith viscosity
            Biharmonic_DSalOrTempDXDPhiDX2A = 0.d0
            Biharmonic_DSalOrTempDYDPhiDY2A = 0.d0
!Casey 220227 from Arash/Rosemary
            TotalArea2 = 0.d0
!            d2_dx2_SalOrTemp_element = 0.d0
!            d2_dy2_SalOrTemp_element = 0.d0
!            d2_dxy_SalOrTemp_element = 0.d0

            N1=NH               !node 1 is always the central node
            qN1=q(N1,k)       
            transN1=inarray(N1,k)     
! af:  
            transN1_z = transN1
! af 01/05/2016:
            Biharmonic_SalOrTempN1_z = Biharmonic_auxiliary_var_SalOrTemp (N1,k)
!            d_dx_SalOrTempN1_z = d_dx_SalOrTemp_nodes (N1,k)
!            d_dy_SalOrTempN1_z = d_dy_SalOrTemp_nodes (N1,k)

            N2=NEITAB(NH,2)     !operate on 1st neighbor
            qN2=q(N2,k)
            transN2=inarray(N2,k)
! find transN2 at the same z coordinate as transN1:
            HsN2trans = DP(N2) + IFNLFA * Eta1(N2)
!Casey 220207 from Rosemary
            IF(HsN2trans.LT.H0) HsN2trans = H0
            SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2trans !equivalent sigma value at neighbor
            inarray_vertical (:) = inarray(N2,:)
!Casey 220207 from Rosemary : Some confusing changes through this
!             section. In her code, the two cubic interpolations are used, not
!             sure why they would be disabled here. Let's try them.
!Casey 221103 : User control over cubic interpolation
            IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
               CALL linear_vertical_interpolation_real ( 
     &             inarray_vertical, SigmaNN, transN2_z, k )
            ELSE ! Default
               CALL cubic_vertical_interpolation ( inarray_vertical, 
     &             SigmaNN, transN2_z, N2, k ) !interp inarray at neighbor
            ENDIF
! af 01/05/2016:
            inarray_vertical (:) = Biharmonic_auxiliary_var_SalOrTemp (N2,:)
!Casey 221103 : User control over cubic interpolation
            IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
               CALL linear_vertical_interpolation_real ( 
     &             inarray_vertical, SigmaNN, Biharmonic_SalOrTempN2_z,
     &             k )
            ELSE ! Default
               CALL cubic_vertical_interpolation ( inarray_vertical,
     &             SigmaNN, Biharmonic_SalorTempN2_z, N2, k )
            ENDIF
!            inarray_vertical (:) = d_dx_SalOrTemp_nodes(N2,:)
!            Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, d_dx_SalOrTempN2_z, k )
!            inarray_vertical (:) = d_dy_SalOrTemp_nodes(N2,:)
!            Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, d_dy_SalOrTempN2_z, k )

            NNFirst=N2          !save these values until end
            qNFirst=qN2
            transNFirst=transN2
! af:
            transNFirst_z = transN2_z
! af 01/05/2016:
            Biharmonic_SalOrTempNFirst_z = Biharmonic_SalOrTempN2_z
!            d_dx_SalOrTempNFirst_z = d_dx_SalOrTempN2_z
!            d_dy_SalOrTempNFirst_z = d_dy_SalOrTempN2_z

            DO N=3,NNEIGH(NH)   !operate on rest of neighbors
               N3=N2            !shift previously computed values
               qN3=qN2
               transN3=transN2
! af:
               transN3_z = transN2_z
! af 01/05/2016:
               Biharmonic_SalOrTempN3_z = Biharmonic_SalOrTempN2_z
!               d_dx_SalOrTempN3_z = d_dx_SalOrTempN2_z
!               d_dy_SalOrTempN3_z = d_dy_SalOrTempN2_z

               N2=NEITAB(NH,N)  !select new neighbor to work on
               qN2=q(N2,k)
               transN2=inarray(N2,k)
! find transN2 at the same z coordinate as transN1:
               HsN2trans = DP(N2) + IFNLFA * Eta1(N2)
!Casey 220207 from Rosemary
               IF(HsN2trans.LT.H0) HsN2trans = H0
               SigmaNN = B + AMB * ( Zk + DP(N2) ) / HsN2trans !equivalent sigma value at neighbor
               inarray_vertical (:) = inarray(N2,:)
!Casey 220207 from Rosemary : Same thing here?
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real ( 
     &                inarray_vertical, SigmaNN, transN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation ( inarray_vertical,
     &                SigmaNN, transN2_z, N2, k ) !interp inarray at neighbor
               ENDIF
! af 01/05/2016:
               inarray_vertical (:) = Biharmonic_auxiliary_var_SalOrTemp (N2,:)
!Casey 221103 : User control over cubic interpolation
               IF( verticalInterpolationSchemeID==VERTICAL_INTERPOLATION_LINEAR )THEN
                  CALL linear_vertical_interpolation_real ( 
     &                inarray_vertical, SigmaNN, 
     &                Biharmonic_SalOrTempN2_z, k )
               ELSE ! Default
                  CALL cubic_vertical_interpolation ( inarray_vertical,
     &                SigmaNN, Biharmonic_SalOrTempN2_z, N2, k )
               ENDIF
!               Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, Biharmonic_SalOrTempN2_z, k )
!               inarray_vertical (:) = d_dx_SalOrTemp_nodes(N2,:)
!               Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, d_dx_SalOrTempN2_z, k )
!               inarray_vertical (:) = d_dy_SalOrTemp_nodes(N2,:)
!               Call linear_vertical_interpolation_real ( inarray_vertical, SigmaNN, d_dy_SalOrTempN2_z, k )

               NEle=NeiTabEle(NH,N-2) !element # defined by nodes N1,N2,N3

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
               ! jgf49.58: NOFF array lookups fail if NELE comes back 0.
               IF (NELE.eq.0) THEN
                  CYCLE
               ENDIF
               NCELE = NODECODE(NH)*NODECODE(N2)
     &              *NODECODE(N3)*NOFF(NELE)
               IF((NEle.NE.0).AND.(NCELE.NE.0)) THEN !if element exists, compute terms
!Casey 220227 from Arash/Rosemary 
                  TotalArea2=TotalArea2+Areas(NEle)
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
!----------------------------------------------------------------------------------------------------------------------------------------------------
! iso-sigma horizontal diffusion:
!                  DtransDXDPhiDX2A=DtransDXDPhiDX2A
!     &                         +(transN1*b1+transN2*b2+transN3*b3)*
!     &                         b1/Areas(NEle)
!                  DtransDYDPhiDY2A=DtransDYDPhiDY2A
!     &                         +(transN1*a1+transN2*a2+transN3*a3)*
!     &                         a1/Areas(NEle)
! horizontal diffusion at the same z levels:
! this if statement is not needed for the new linear vertical interpolation scheme
!                  If ( (abs( transN1_z - (-999.)) < diff_tol ) .OR. (abs( transN2_z - (-999.)) < diff_tol ) .OR. 
!     &                 (abs( transN3_z - (-999.)) < diff_tol ) ) Then
                     ! element in inactive -- add very small diffusion
! ********** **********
!      diffusion_ratio = 0.20d0
! ********** **********
!                  DtransDXDPhiDX2A=DtransDXDPhiDX2A
!     &                         +(transN1*b1+transN2*b2+transN3*b3)*
!     &                         b1/Areas(NEle) * diffusion_ratio
!                  DtransDYDPhiDY2A=DtransDYDPhiDY2A
!     &                         +(transN1*a1+transN2*a2+transN3*a3)*
!     &                         a1/Areas(NEle) * diffusion_ratio
! ********** **********
!                  Else
                     ! element is active with engineered linear interpolation
!Casey 220227 from Arash/Rosemary
                 IF((transN2_z.GT.-998.d0).and.(transN3_z.GT.-998.d0))THEN
                     DtransDXDPhiDX2A=DtransDXDPhiDX2A
     &                         +(transN1_z*b1+transN2_z*b2+transN3_z*b3)*
     &                         b1/Areas(NEle)
                     DtransDYDPhiDY2A=DtransDYDPhiDY2A
     &                         +(transN1_z*a1+transN2_z*a2+transN3_z*a3)*
     &                         a1/Areas(NEle)
!Casey 220227 from Arash/Rosemary
                 ENDIF
!                  End If
! af 01/05/2016:
!Casey 220227 from Arash/Rosemary
                 IF((Biharmonic_SalOrTempN2_z.GT.-998.d0).and.(Biharmonic_SalorTempN3_z.GT.-998.d0))THEN
                     Biharmonic_DSalOrTempDXDPhiDX2A = Biharmonic_DSalOrTempDXDPhiDX2A
     &                 +(Biharmonic_SalOrTempN1_z*b1+Biharmonic_SalOrTempN2_z*b2+Biharmonic_SalOrTempN3_z*b3)*b1/Areas(NEle)
                     Biharmonic_DSalOrTempDYDPhiDY2A = Biharmonic_DSalOrTempDYDPhiDY2A
     &                 +(Biharmonic_SalOrTempN1_z*a1+Biharmonic_SalOrTempN2_z*a2+Biharmonic_SalOrTempN3_z*a3)*a1/Areas(NEle)
!Casey 220227 from Arash/Rosemary
                 ENDIF

!                     d2_dx2_SalOrTemp_element = d2_dx2_SalOrTemp_element + d_dx_SalOrTempN1_z*b1+d_dx_SalOrTempN2_z*b2+d_dx_SalOrTempN3_z*b3
!                     d2_dy2_SalOrTemp_element = d2_dy2_SalOrTemp_element + d_dy_SalOrTempN1_z*a1+d_dy_SalOrTempN2_z*a2+d_dy_SalOrTempN3_z*a3
!                     d2_dxy_SalOrTemp_element = d2_dxy_SalOrTemp_element + d_dy_SalOrTempN1_z*b1+d_dy_SalOrTempN2_z*b2+d_dy_SalOrTempN3_z*b3
!----------------------------------------------------------------------------------------------------------------------------------------------------

               ENDIF
            END DO

            N3=N2               !wrap back to beginning to get final contribution
            qN3=qN2
            transN3=transN2
! af:
            transN3_z = transN2_z
! af 01/05/2016:
            Biharmonic_SalOrTempN3_z = Biharmonic_SalOrTempN2_z
!            d_dx_SalOrTempN3_z = d_dx_SalOrTempN2_z
!            d_dy_SalOrTempN3_z = d_dy_SalOrTempN2_z

            N2=NNFIRST
            qN2=qNFirst
            transN2=transNFirst
! af:
            transN2_z = transNFirst_z
! af 01/05/2016:
            Biharmonic_SalOrTempN2_z = Biharmonic_SalOrTempNFirst_z
!            d_dx_SalOrTempN2_z = d_dx_SalOrTempNFirst_z
!            d_dy_SalOrTempN2_z = d_dy_SalOrTempNFirst_z

            NEle=NeiTabEle(NH,NNeigh(NH)-1)

Casey: Added the computation of "NCELE" and the second half of the IF statement.
C
            ! jgf49.58: NOFF array lookups fail if NELE comes back 0.
            IF (NELE.ne.0) THEN
               NCELE = NODECODE(NH)*NODECODE(N2)
     &             *NODECODE(N3)*NOFF(NELE)
            ENDIF
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
!----------------------------------------------------------------------------------------------------------------------------------------------------
! iso-sigma horizontal diffusion:
!                  DtransDXDPhiDX2A=DtransDXDPhiDX2A
!     &                         +(transN1*b1+transN2*b2+transN3*b3)*
!     &                         b1/Areas(NEle)
!                  DtransDYDPhiDY2A=DtransDYDPhiDY2A
!     &                         +(transN1*a1+transN2*a2+transN3*a3)*
!     &                         a1/Areas(NEle)
! horizontal diffusion at the same z levels: if statement is not needed - see linear vertical interpolation above.
!                  If ( (abs( transN1_z - (-999.)) < diff_tol ) .OR. (abs( transN2_z - (-999.)) < diff_tol ) .OR. 
!     &                 (abs( transN3_z - (-999.)) < diff_tol ) ) Then
!                     ! element in inactive -- add very small diffusion
! ********** **********
!      diffusion_ratio = 0.20d0
! ********** **********
!                  DtransDXDPhiDX2A=DtransDXDPhiDX2A
!     &                         +(transN1*b1+transN2*b2+transN3*b3)*
!     &                         b1/Areas(NEle) * diffusion_ratio
!                  DtransDYDPhiDY2A=DtransDYDPhiDY2A
!     &                         +(transN1*a1+transN2*a2+transN3*a3)*
!     &                         a1/Areas(NEle) * diffusion_ratio
! ********** **********
!                  Else
                     ! element is active with engineered linear interpolation
!Casey 220227 from Arash/Rosemary
                  IF((transN2_z.GT.-998.d0).and.(transN3_z.GT.-998.d0))THEN
                     DtransDXDPhiDX2A=DtransDXDPhiDX2A
     &                         +(transN1_z*b1+transN2_z*b2+transN3_z*b3)*
     &                         b1/Areas(NEle)
                     DtransDYDPhiDY2A=DtransDYDPhiDY2A
     &                         +(transN1_z*a1+transN2_z*a2+transN3_z*a3)*
     &                         a1/Areas(NEle)
!Casey 220227 from Arash/Rosemary
                  ENDIF
!                  End If

! af 01/05/2016:
!Casey 220227 from Arash/Rosemary
                  IF((Biharmonic_SalorTempN2_z.GT.-998.d0).and.(Biharmonic_SalorTempN3_z.GT.-998.d0))THEN
                  Biharmonic_DSalOrTempDXDPhiDX2A = Biharmonic_DSalOrTempDXDPhiDX2A
     &              +(Biharmonic_SalOrTempN1_z*b1+Biharmonic_SalOrTempN2_z*b2+Biharmonic_SalOrTempN3_z*b3)*b1/Areas(NEle)
                  Biharmonic_DSalOrTempDYDPhiDY2A = Biharmonic_DSalOrTempDYDPhiDY2A
     &              +(Biharmonic_SalOrTempN1_z*a1+Biharmonic_SalOrTempN2_z*a2+Biharmonic_SalOrTempN3_z*a3)*a1/Areas(NEle)
!Casey 220227 from Arash/Rosemary
                  ENDIF

!                  d2_dx2_SalOrTemp_element = d2_dx2_SalOrTemp_element + d_dx_SalOrTempN1_z*b1+d_dx_SalOrTempN2_z*b2+d_dx_SalOrTempN3_z*b3
!                  d2_dy2_SalOrTemp_element = d2_dy2_SalOrTemp_element + d_dy_SalOrTempN1_z*a1+d_dy_SalOrTempN2_z*a2+d_dy_SalOrTempN3_z*a3
!                  d2_dxy_SalOrTemp_element = d2_dxy_SalOrTemp_element + d_dy_SalOrTempN1_z*b1+d_dy_SalOrTempN2_z*b2+d_dy_SalOrTempN3_z*b3
!----------------------------------------------------------------------------------------------------------------------------------------------------
            ENDIF

            IF(TotalArea2.EQ.0.) THEN
! af: why is this a complex value? arash changed them to real.
!               LAdvectrans(k)=(0.d0,0.d0)
!               LDiffusion(k)=(0.d0,0.d0)
               LAdvectrans(k) = 0.d0
               LDiffusion(k)  = 0.d0
! arash 12/14/2015
               Biharmonic_LDiffusion(k) = 0.d0
!Casey 220227 from Rosemary
               LDiffusion_WetDryInterface(k) = 0.d0
            ELSE

C  kmd48.33bc added sponge layer for the lateral advective terms
               LAdvectrans(k)=sponge(NH)*IFNLCT*(UnDtransDX2A+VnDtransDY2A)/
     &                        TotalArea2   
               LDiffusion(k)=3.d0*diffhoriz*(DtransDXDPhiDX2A+
     &                       DtransDYDPhiDY2A)/TotalArea2
! arash 01/05/2016: modified Leith viscosity - this is  > > > W R O N G < < < and is written for velocities
!               ! 0. preliminaries
!               modified_Leith_delta_x = sqrt ( TotalArea2 * 0.5d0 ) ! grid length
!               modified_Leith_Lambda_vorticity  = 40.0d0 !10.0d0             ! tuning constant O(1)
!               modified_Leith_Lambda_divergence = 40.0d0 !10.0d0             ! the modification - see MITgcm manual
!               ! 1. nodal second derivatives: during the computations, we multiplied them by 2A.
!               d2_dx2_SalOrTemp_element = d2_dx2_SalOrTemp_element / TotalArea2
!               d2_dy2_SalOrTemp_element = d2_dy2_SalOrTemp_element / TotalArea2
!               d2_dxy_SalOrTemp_element = d2_dxy_SalOrTemp_element / TotalArea2
!               ! 2. magnitude of the gradient of divergence and vorticity
!               modified_Leith_grad_vorticity = sqrt( ( aimag(d2_dx2_SalOrTemp_element) - real(d2_dxy__element) )**2
!     &                                       +       ( aimag(d2_dxy_SalOrTemp_element) - real(d2_dy2_q_element) )**2 )
!               modified_Leith_grad_divergence= sqrt( (  real(d2_dx2_SalOrTemp_element) +aimag(d2_dxy_q_element) )**2
!     &                                       +       (  real(d2_dxy_SalOrTemp_element) +aimag(d2_dy2_q_element) )**2 )
!               ! 3. compute the modified Leith viscosity
!               Biharmonic_viscosity = ( modified_Leith_delta_x / 3.141592654d0 )**3 *
!     &                           sqrt ( modified_Leith_Lambda_vorticity**6  * modified_Leith_grad_vorticity**2
!     &                                + modified_Leith_Lambda_divergence**6 * modified_Leith_grad_divergence**2 )
!               ! 4. congratulations! we have computed modified Leith viscosity!
!Casey 221103 : User control over Leith viscosity
                IF( .NOT.useLeithViscosity )THEN
                   Biharmonic_viscosity = diffhoriz
                ELSE ! Default
                Biharmonic_viscosity = Biharmonic_viscosity_modified_Leith ( NH, k ) * 1.00d0 ! to prevent loosing eddies
                ENDIF
                Biharmonic_LDiffusion(k) = -3.d0 * Biharmonic_viscosity * 
     &                                     (Biharmonic_DSalOrTempDXDPhiDX2A+Biharmonic_DSalOrTempDYDPhiDY2A) / TotalArea2
            ENDIF

!Casey 220523 : End of commenting.

#endif

         ENDDO

!Casey 220207 from Rosemary : To resolve diffusion problem along wet/dry
!             fronts, if any of the neigboring elements are dry, then preserve
!             the harmonic formulation (i.e. do not apply biharmonic formulation)
!             for Lateral Diffusion term.
         LDiffusion_WetDryInterface(:) = LDiffusion(:)
! |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
! arash 1/19/2016: Use Biharmonic diffusion instead of a harmonic one
! We may want to make this cleaner later by having an option in the fort.15 file
!Casey 221103 : User control over bhiarmonic operator for diffusion
         IF( diffusionOperatorID==DIFFUSION_OPERATOR_LAPLACIAN )THEN
            CONTINUE
         ELSE ! Default
         LDiffusion(:) = Biharmonic_LDiffusion(:)
         ENDIF
! |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
!Casey 220207 from Rosemary : For vertices connected to dry elements, preserve harmonic formulation
          DO N = 1,NNeigh(NH)
             NEle = NeiTabEle(NH,N)
             IF (NEle.NE.0) THEN
                IF (NOFF(NELE).EQ.0) THEN
!                  LDiffusion(:) = LDiffusion_WetDryInterface(:)
                   LDiffusion(:) = 0.d0
                   LAdvectrans(:) = 0.d0
                   EXIT
                ENDIF
             ENDIF
          ENDDO

!----------------------------------------------------------------------------------------------------------------------------------------------------

c     Finished computing advection and diffusion terms in the transport
c     equation at each level in the vertical
!----------------------------------------------------------------------------------------------------------------------------------------------------

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
!Casey 220207 from Rosemary : Heat flux was turned off at the 
!             river boundary nodes were we are specifying an incoming flux
      IF (BNDBCRiver) THEN
         DO BCSEC=1, bcriversec
            DO NBC=1,BCRNVELL(BCSEC)
               NM1=BCRNBVV(BCSEC,NBC)
               IF(NH.EQ.NM1) topbcs = 0.d0
            ENDDO
         ENDDO
      ENDIF
!Casey 220207 from Rosemary : For vertices connected to dry elements, turn off heat flux
      IF(field.EQ.3) THEN
         DO N = 1,NNeigh(NH)
            NEle = NeiTabEle(NH,N)
            IF (NEle.NE.0) THEN
               IF (NOFF(NELE).EQ.0) THEN
         topbcs = 0.d0
                  EXIT
      END IF
            ENDIF
         ENDDO
      ENDIF

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

!Casey 220227 from Rosemary
!        Frtrans(NFEN) = (1.d0*inarray(NH,k-1)
!    &        -DelT*(LAdvectrans(k-1)+LDiffusion(k-1)))*Inm(k,1)
!    &        + (1.d0*inarray(NH,k)
!    &        -DelT*(LAdvectrans(k)  +LDiffusion(k)))*Inm(k,2)
!    &        - DelT*(VAdvectrans(k))-DT1MAlp4*VDiffusion(k)
!    &        + topbcs + WBC(NH,NFEN)
         Frtrans(NFEN) = (1.d0*inarray(NH,NFEN-1)
     &        -DelT*(LAdvectrans(NFEN-1)+LDiffusion(NFEN-1)))*Inm(NFEN,1)
     &        + (1.d0*inarray(NH,NFEN)
     &        -DelT*(LAdvectrans(NFEN)  +LDiffusion(NFEN)))*Inm(NFEN,2)
     &        - DelT*(VAdvectrans(NFEN))-DT1MAlp4*VDiffusion(NFEN)
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
  
       CALL ADC_TRIDAG2(Mkm1trans,Mktrans,Mkp1trans,Frtrans,Gammatrans,NFEN)

C      Need to save the transport results by putting them into 
C      outgoing variables 
     
         DO k=1,NFEN
            outarray(NH,k) = Gammatrans(k)
         END DO
     
!Casey 220227 from Rosemary
         DO k=1,NFEN
           IF(NODECODE(NH).EQ.1)THEN
            outarray(NH,k) = Gammatrans(k)
           ELSE
             outarray(NH,k) = inarray(NH,k)
           ENDIF
         END DO

      ENDDO

C     Finish loop over horizontal nodes to compute the horizontal velocity

      RETURN
      END SUBROUTINE TRANS_3D

C  This subroutine utilized to solve the matrix - just changed the name
C  to ADC_TRIDAG2
C***********************************************************************
C Subroutine ADC_TRIDAG2                                                    *
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
      SUBROUTINE ADC_TRIDAG2(A,B,C,R,U,N)
      USE GLOBAL, ONLY: ScreenUnit
      IMPLICIT NONE
      INTEGER :: J, N

      REAL(8) :: A(N),B(N),C(N),R(N),U(N)
      REAL(8) :: BET,GAM(N)

c     
c     check for zero elements on diagonal
c     
      DO J=1,N
         if(B(j).EQ.0.) then
            write(screenunit,*) 'Problem in Tridag Solver.  ',
     +           'B array value in row ',j,' = 0'
            CALL EXIT(1)
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
            CALL EXIT(1)
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


C*************************************************************************
C     Subroutine to interpolate salinity/temperature to a specified
C     sigma value (SigmaNN) given an initial guess of which sigma
C     level is closest to the specified value.
C     This subroutine is almost identical to SUBROUTINE ZSURFBUOY, written by R.L.
C
C                                    a.f.  23/06/15
C*************************************************************************

      Subroutine cubic_vertical_interpolation ( inarray_vertical, SigmaNN, transN2_z, NN, J )

      USE GLOBAL_3DVS
      IMPLICIT NONE
      REAL(8) :: transN2_z
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
      Integer  :: k1, k2, k3, k4                 ! nodes used for cubic polynomial
      REAL(8) :: f1, f2, f3, f4                 ! shape functions
      REAL(8) :: top, bot
      Integer  :: loop_situation = 0

      REAL(8) :: inarray_vertical
      Dimension  inarray_vertical ( NFEN )

      call setMessageSource("cubic_vertical_interpolation")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IDiag=0

      IF(SigmaNN.LE.(1.00010d0*b)) THEN !if into ground then skip
         SigBelo=-999
         SigAbov=-999
         transN2_z=-999.
         GOTO 100
      ENDIF
      IF((SigmaNN.GT.(1.0001d0*b)).AND.(SigmaNN.LE.b)) THEN !at bottom then use bottom
         LBelo=1
         transN2_z=inarray_vertical(LBelo)
         SigBelo=b
         SigAbov=b
         GOTO 100
      ENDIF

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
! arash 150513
!Casey 220207 from Rosemary : Skip Arash's formulation, try different
!     IF(SigmaNN.GE.(1.0001d0*a)) THEN            ! in the air - non-existent 
!        transN2_z=-999.0d0
!        GOTO 100
!     ENDIF
!     IF((SigmaNN.GE.a).and.(SigmaNN.LT.(1.0001d0*a))) THEN     ! surface
!        LAbov=NFEN
!        transN2_z=inarray_vertical(LAbov) ! this is zero
!        SigBelo=a
!        SigAbov=a
!        GOTO 100
!     ENDIF
      IF(SigmaNN.GE.a) THEN
         LAbov=NFEN
         transN2_z=inarray_vertical(LAbov)
         SigBelo=a
         SigAbov=a
         GOTO 100
      ENDIF

! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--


      LTry=J                    !start search for SIGABOV and SIGBELO
      SigTry=Sigma(LTry)
      IF(SigmaNN.GT.SigTry) THEN !too low
         loop_situation = 1
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
         loop_situation = 2
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

! arash added this because of the terminations at the SURA model. 
! this shouldn't be needed though.
! 0.
!      If ( loop_situation == 0 ) Then
!         transN2_z = -999.0d0
!         GoTo 100
!      End If
! 1.
!      If ( (LAbov == -999) .or. (LBelo == -999) ) Then
!         transN2_z = -999.0d0
!         write(*,*) 'LAbov, LBelo = ', LAbov, LBelo
!         write(*,*) 'SigBelo, SigAbov = ', SigBelo, SigAbov
!         write(*,*) 'SigmaNN, J =', SigmaNN, J
!         GoTo 100
!      End If


! use a cubic polynomial interpolation, instead of the above (99) linear interpolation
! 1- decide which nodes to pick: k1:lowest, k2, k3, k4:top-most
      If ( LAbov == NFEN ) Then
         k4 = NFEN
         k3 = NFEN - 1
         k2 = NFEN - 2
         k1 = NFEN - 3
      Else If ( LBelo == 1 ) Then
         k4 = 4
         k3 = 3
         k2 = 2
         k1 = 1
      Else
         k4 = LAbov + 1
         k3 = LAbov
         k2 = LBelo
         k1 = LBelo - 1
      End If

! 2- compute the shape functions
        f1 = ((SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
     &     / ((Sigma(k1)-Sigma(k2)) * (Sigma(k1)-Sigma(k3)) * (Sigma(k1)-Sigma(k4)))

        f2 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k3)) * (SigmaNN - Sigma(k4)))
     &     / ((Sigma(k2)-Sigma(k1)) * (Sigma(k2)-Sigma(k3)) * (Sigma(k2)-Sigma(k4)))

        f3 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k4)))
     &     / ((Sigma(k3)-Sigma(k1)) * (Sigma(k3)-Sigma(k2)) * (Sigma(k3)-Sigma(k4)))

        f4 = ((SigmaNN - Sigma(k1)) * (SigmaNN - Sigma(k2)) * (SigmaNN - Sigma(k3)))
     &     / ((Sigma(k4)-Sigma(k1)) * (Sigma(k4)-Sigma(k2)) * (Sigma(k4)-Sigma(k3)))

! 3- interpolate
      transN2_z = f1 * inarray_vertical(k1) + f2 * inarray_vertical(k2) 
     &          + f3 * inarray_vertical(k3) + f4 * inarray_vertical(k4)
!Casey 220227 from Arash/Rosemary : They used the lines below, 
!          which look like a linear interpolation? Let's keep the cubic.
!     transN2_z = (inarray_vertical(LAbov) - inarray_vertical(LBelo))
!    &     *(SigmaNN-SigBelo)/(SigAbov-SigBelo) +inarray_vertical(LBelo)


! /--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--/--
 100  CONTINUE

      IF(IDiag.EQ.2) THEN
         WRITE(2,*) '******** cubic_vertical_interpolation **********'
         WRITE(2,*) '     NH  NV  SigmaNN   SigBelo   SigAbov',
     &        '      transN2_z'
         WRITE(2,777) NN,J,SigmaNN,SigBelo,SigAbov,transN2_z
 777     FORMAT(I7,I5,3(F10.3),E14.5)
      ENDIF

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C**********************************************************************
      END SUBROUTINE cubic_vertical_interpolation
C**********************************************************************

