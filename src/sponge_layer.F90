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
! but WITHOUT ANY WARRANTY without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
!
!--------DW for sponge layer
!
module SPONGELAYER

   use GLOBAL, only: STATIM, IFNLFA, ETA2, UU2, VV2, QX2, QY2
   use ADC_CONSTANTS, only: G
   use MESH, only: NP, NM, NE, DP, X
   use NodalAttributes, only: LoadAbsLayerSigma, &
                              SSIGMA_ETA => absorblayer_sigma_eta, &
                              SSIGMA_MNX => absorblayer_sigma_mnx, &
                              SSIGMA_MNY => absorblayer_sigma_mny, &
                              NumNodesAbsLayer, &
                              AbsLayer_Eta_NodesID, &
                              AbsLayer_MnX_NodesID, &
                              AbsLayer_MnY_NodesID, &
                              AbsLayerType, AbsLayerNBF, &
                              AbsLayerEtaAMIG, AbsLayerEtaFF, AbsLayerEtaFACE, &
                              AbsLayerEtaEMO, AbsLayerEtaEFA, &
                              AbsLayerQAMIG, AbsLayerQFF, AbsLayerQFACE, &
                              AbsLayerQxEMO, AbsLayerQxEFA, &
                              AbsLayerQyEMO, AbsLayerQyEFA !, GWCE_spg_alt

   use QUADRATURETRI, only: MATARR, VECARR, &
                            AllocMatArr, AllocVecArr, CompElmMsfh, GetDefaultCub2D

   implicit none

   private

   !c For AbsLayerType == 1
   real(8), allocatable :: pert_eta_abslayer(:), pert_vel_abslayer(:)

   real(8), allocatable :: eta2_abslayer(:), eta1_abslayer(:), &
                           eta0_abslayer(:) ! (s + 1, s, s-1) !

   real(8), allocatable :: uu2_abslayer(:), uu1_abslayer(:) ! (s+1,s) !
   real(8), allocatable :: vv2_abslayer(:), vv1_abslayer(:) ! (s+1,s) !

   real(8), private, allocatable :: WorkVecEta(:)
   real(8), private, allocatable :: WorkVecUU(:)
   real(8), private, allocatable :: WorkVecVV(:)
   real(8), private, allocatable :: ArbEta0(:), ArbEta(:)
   real(8), private, allocatable :: ArbUU0(:), ArbUU(:)
   real(8), private, allocatable :: ArbVV0(:), ArbVV(:)

   real(8) :: SPTIMINC, SPTIME1, SPTIME2

   logical, allocatable :: spgflag(:)

   !------------- accurate calculation of \int_{} \sigma \phi_{i} \phi_{j} dx
   type(MATARR), allocatable :: MsElmFh(:)
   type(VECARR), allocatable :: LumpMsElmFh(:)

   !-------------DW: sponge layers (2015)
   real(8), allocatable :: MCOEFSPG(:, :) !c mass matrix
   real(8), allocatable :: MCOEFSPGD(:) !c lumped mass matrix

   integer, parameter :: sponge_dis_mthd = 1 !c 1 -- solve GWCE c!
   !c 0 -- Operator Spliling
   logical, parameter :: adjsigval = .false.
   logical :: NO_MET_IN_SPONGE = .false.
   logical :: NO_BPG_IN_SPONGE = .false.

   public :: ALLOC_MAINSPG, ALLOC_MAINSPG_LUMPED, &
             SPONGE_OPSPLIT0, Adjust_Sponge_Sigma, &
             FLAGSPONGEELEM, ComputeSpongeMsElmMat, &
             SpongeLayerRelatedPrep, GETABSLAYEREXT, &
             SPONGE_SHIFT_SOLN, LAYERSOLTEST1, &
             NO_MET_IN_SPONGE, NO_BPG_IN_SPONGE, &
             NumNodesAbsLayer, sponge_dis_mthd, &
             uu1_AbsLayer, vv1_AbsLayer, uu2_AbsLayer, &
             vv2_AbsLayer, eta0_AbsLayer, &
             eta1_AbsLayer, eta2_AbsLayer, &
             LoadAbsLayerSigma

contains

   !
   ! Allocate space for array used in the sponge layer through an
   ! operator spliting technique:
   ! - precomputed mass matrix and lumped-mass matrix
   !
   subroutine ALLOC_MAINSPG_LUMPED()
      use SIZES, only: MNP
      implicit none

      allocate (MCOEFSPGD(MNP))

      return
   end subroutine ALLOC_MAINSPG_LUMPED

   subroutine ALLOC_MAINSPG()
      use SIZES, only: MNP, MNEI
      implicit none

      allocate (MCOEFSPG(MNP, MNEI))

      return
   end subroutine ALLOC_MAINSPG
   !--------DW

   subroutine Adjust_Sponge_Sigma(DT)
      implicit none

      real(8), intent(in) :: DT

      real(8) :: SIGMAX
      logical, save :: first = .true.

      SIGMAX = 1.d0

      if (adjsigval) then
         SIGMAX = 0.9d0/DT
      end if

      if (first) then
      !!$ SSIGMA = SIGMAX*SSIGMA
         SSIGMA_ETA = SIGMAX*SSIGMA_ETA

         first = .false.

         select case (sponge_dis_mthd)
         case (0)
            write (16, *) "Sponge Layer: use an operator"// &
               "splitting in discretization"
         case (1)
            write (16, *) "Sponge Layer: implement in GWCE"// &
               "and Non conservative Momentum equation"
         end select

      end if

      return
   end subroutine Adjust_Sponge_Sigma

   !
   !.....Simplest implementation: lumed both lhs and rhs matrix
   !
   !
   subroutine SPONGE_OPSPLIT0(DT)
      !
      implicit none

      !
      ! Dummy variables
      !
      real(8), intent(IN) :: DT
      !
      !     Local variables
      !
      integer :: I
      real(8) :: DTD2, SG(3), CNLHS(3), CNRHS, SIGMAX, H2
      real(8), parameter :: smalleps = 1.d-12

      DTD2 = DT/2.d0
    !!$ SIGMAX = 0.8/DT
    !!$ SIGMAX = 1.D0/DT
      SIGMAX = 1.d0

      do I = 1, NP
         SG(1) = SIGMAX*SSIGMA_ETA(I, 1)
         SG(2) = SIGMAX*SSIGMA_MNX(I, 1)
         SG(3) = SIGMAX*SSIGMA_MNY(I, 1)

         if (sum(abs(SG)) > smalleps) then
            !.....Crank-Nicolson
            ! CNLHS = 1.D0/(1.D0 + DTD2*SG)
            ! CNRHS = (1.D0 - DTD2*SG)

            !.....Backward Euler
            CNLHS = 1.d0/(1.d0 + DT*SG)
            CNRHS = 1.d0
            !
            ! Update \zeta
            ! d(\zeta)/dt = - sigma*\zeta
            ETA2(I) = CNLHS(1)*CNRHS*ETA2(I)

            !
            ! Update U,V
            ! dU/dt = -sigma*U
            ! dV/dt = -sigma*U
            !
            UU2(I) = CNLHS(2)*CNRHS*UU2(I)
            VV2(I) = CNLHS(3)*CNRHS*VV2(I)

            ! update QX
            H2 = DP(I) + dble(IFNLFA)*ETA2(I)
            QX2(I) = UU2(I)*H2
            QY2(I) = VV2(I)*H2

         end if
      end do

      return
   end subroutine SPONGE_OPSPLIT0

   subroutine FLAGSPONGEELEM()
      implicit none

      integer :: IE, idx(3)
      real(8), parameter :: SMALL = 1.0d-10

      if (NumNodesAbsLayer(1) > 0) then
         open (unit=550, FILE='spgelment.dat')

         !     Flag element inside a sponge layer
         allocate (spgflag(NE))
         spgflag = .false.
         do IE = 1, NE
            idx = NM(IE, 1:3)

            if (sum(abs(SSIGMA_ETA(idx, 1))) > SMALL) then
               spgflag(IE) = .true.

               write (550, *) IE
            end if
         end do
         close (550)

      end if
   end subroutine FLAGSPONGEELEM

   subroutine ComputeSpongeMsElmMat(pc)
      implicit none

      integer, intent(in) :: pc

      !c local variables c!
      integer :: II, IE, idx(3)
      real(8) :: sigmaval(3)

      logical, save :: FirstEnter = .true.

      if (NumNodesAbsLayer(1) > 0) then
         if (FirstEnter) then
            allocate (MsElmFh(NE))
            allocate (LumpMsElmFh(NE))
         end if
         call GetDefaultCub2D(pc)

         do IE = 1, NE
            if (spgflag(IE)) then

               idx = NM(IE, 1:3)
               sigmaval = SSIGMA_ETA(idx, 1)

               if (FirstEnter) then
                  nullify (MsElmFh(IE)%ARRVAL)
                  nullify (LumpMsElmFh(IE)%VECVAL)
               end if

               call AllocMatArr(MsElmFh(IE), 3)
               call AllocVecArr(LumpMsElmFh(IE), 3)

               call CompElmMsfh(MsElmFh(IE)%ARRVAL, sigmaval)
               do II = 1, 3
                  LumpMsElmFh(IE)%VECVAL(II) = &
                     sum(MsElmFh(IE)%ARRVAL(II, 1:3))
               end do
            end if
         end do
      end if
      FirstEnter = .false.

      return
   end subroutine ComputeSpongeMsElmMat

   subroutine SpongeLayerRelatedPrep()
      use SIZES, only: MNP
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      integer :: I, J, K, NODEID
      integer :: NAL, IMIN(1)

      integer :: SSDUM(4)
      integer, parameter :: allowabstype(4) = [0, 1, -1, 2]

      real(8), pointer :: SSIGM_PT(:, :)
      integer, pointer :: AbsLayerid_PTS(:)

      if (NumNodesAbsLayer(1) > 0) then

         !IF (  GWCE_spg_alt ) THEN
         !   !c GWCE_spg_alt == .TRUE.
         !   !c   Use a more accurate integration for the
         !   !c   sponge_layer term
         !   CGWCE_New = .FALSE.
         !   CGWCE_New_SPG_ALT_SOL = .TRUE.
         !END IF

         NAL = ubound(allowabstype, 1)

         SSDUM = 0

         do K = 1, 3
            select case (K)
            case (1)
               SSIGM_PT => SSIGMA_ETA
               AbsLayerid_PTS => AbsLayer_Eta_NodesID
            case (2)
               SSIGM_PT => SSIGMA_MNX
               AbsLayerid_PTS => AbsLayer_MnX_NodesID
            case (3)
               SSIGM_PT => SSIGMA_MNY
               AbsLayerid_PTS => AbsLayer_MnY_NodesID
            end select

            do I = 1, NumNodesAbsLayer(K)
               NODEID = AbsLayerid_PTS(I)

               do J = 1, NAL
                  SSDUM(J) = SSDUM(J) + &
                             nint(SSIGM_PT(NODEID, 2) - dble(allowabstype(J)))
               end do

            end do
         end do

         ! allocate memormy for external forcing function !
         !  MNP - in global.F, maximum no. of nodes       !
         allocate (eta2_abslayer(MNP))
         eta2_abslayer = 0.d0
         allocate (eta1_abslayer(MNP))
         eta1_abslayer = 0.d0
         allocate (eta0_abslayer(MNP))
         eta0_abslayer = 0.d0

         allocate (uu2_abslayer(MNP))
         uu2_abslayer = 0.d0
         allocate (uu1_abslayer(MNP))
         uu1_abslayer = 0.d0
         allocate (vv2_abslayer(MNP))
         vv2_abslayer = 0.d0
         allocate (vv1_abslayer(MNP))
         vv1_abslayer = 0.d0

         allocate (WorkVecEta(NumNodesAbsLayer(1)))
         allocate (WorkVecUU(NumNodesAbsLayer(2)))
         allocate (WorkVecVV(NumNodesAbsLayer(3)))

         IMIN = minloc(abs(SSDUM))
         if (abs(SSDUM(IMIN(1))) > 0) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="Unknown sponge-layer type")
         else
            AbsLayerType = allowabstype(IMIN(1))
         end if

         select case (AbsLayerType)
         case (0)
            write (16, *) "Sponge Layer: Fully Absorbing"
            AbsLayerNBF = 0
            NO_MET_IN_SPONGE = .true.
            NO_BPG_IN_SPONGE = .true.
            write (16, *) "MET and BPG turned off in sponge!"
         case (1)
            write (16, *) "Sponge Layer: Generating(Tides)-Absorbing"
            NO_MET_IN_SPONGE = .true.
            NO_BPG_IN_SPONGE = .true.
            write (16, *) "MET and BPG turned off in sponge!"
            call READFORT53001()
            call READFORT54001()
         case (2)
            write (16, *) "Sponge Layer: Generating(Tides+Arb)-Absorbing"
            NO_MET_IN_SPONGE = .true.
            NO_BPG_IN_SPONGE = .true.
            write (16, *) "MET and BPG turned off in sponge!"
            call READFORT53001()
            call READFORT54001()
            call READFORT2001(.true.)
         case (-1)
            write (16, *) "Sponge Layer: Generating(Arbitrary)-Absorbing"
            write (16, *) "MET and BPG kept on in sponge!"
            call READFORT2001(.true.)
            AbsLayerNBF = 0
         end select
      end if
   end subroutine SpongeLayerRelatedPrep

   !
   ! NOTE:
   !   fort.53001 and fort.54001
   !  have the same format as fort.53 and fort.54, respectively.
   !
   subroutine READFORT53001()
      use SIZES, only: INPUTDIR
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      integer :: FUNIT, IOS

      integer :: KK, JJ, NPTMP, KTMP
      character(LEN=300) :: IGNOREMSG
      real(8) :: PIL, DEGRAD
      real(8), parameter :: eps = epsilon(1d0)

      PIL = acos(-1.d0)
      DEGRAD = PIL/180.d0

      FUNIT = 531
      open (UNIT=FUNIT, &
            FILE=trim(INPUTDIR)//'/'//'fort.53001', IOSTAT=IOS)

      read (FUNIT, *) AbsLayerNBF
      if (.not. allocated(AbsLayerEtaAMIG)) then
         allocate (AbsLayerEtaAMIG(AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerEtaFF)) then
         allocate (AbsLayerEtaFF(AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerEtaFACE)) then
         allocate (AbsLayerEtaFACE(AbsLayerNBF))
      end if
      if (.not. allocated(pert_eta_abslayer)) then
         allocate (pert_eta_abslayer(AbsLayerNBF))
      end if

      do KK = 1, AbsLayerNBF
         read (FUNIT, *) AbsLayerEtaAMIG(KK), AbsLayerEtaFF(KK), &
            AbsLayerEtaFACE(KK), IGNOREMSG

         AbsLayerEtaFACE(KK) = AbsLayerEtaFACE(KK)*DEGRAD

         if (AbsLayerEtaAMIG(KK) <= eps) then
            pert_eta_abslayer(kk) = 0.d0
         else
            pert_eta_abslayer(kk) = 2.d0*PIL/AbsLayerEtaAMIG(KK)
         end if
      end do

      read (FUNIT, *) NPTMP
      if (NPTMP /= NumNodesAbsLayer(1)) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Error in READFORT53001(): 1")
      end if

      if (.not. allocated(AbsLayerEtaEMO)) then
         allocate (AbsLayerEtaEMO(NPTMP, AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerEtaEFA)) then
         allocate (AbsLayerEtaEFA(NPTMP, AbsLayerNBF))
      end if

      do KK = 1, NPTMP
         read (FUNIT, *) KTMP

         ! PRINT*, KTMP, AbsLayer_Eta_NodesID(KK)
         if (KTMP /= AbsLayer_Eta_NodesID(KK)) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="Error in READFORT53001(): 2")
         end if

         do JJ = 1, AbsLayerNBF
            read (FUNIT, *) AbsLayerEtaEMO(KK, JJ), AbsLayerEtaEFA(KK, JJ)

            !
            ! AbsLayerEtaEFA(KK,JJ) = AbsLayerEtaEFA(KK,JJ)*DEGRAD
         end do

      end do
      AbsLayerEtaEFA = AbsLayerEtaEFA*DEGRAD

      close (FUNIT)
      write (16, *) "DONE READING FORT.53001"

      return
   end subroutine READFORT53001

   subroutine READFORT54001()
      use SIZES, only: INPUTDIR
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      integer :: FUNIT, IOS

      integer :: KK, JJ, NPTMP, KTMP, NBF
      character(LEN=300) :: IGNOREMSG
      real(8) :: PIL, DEGRAD
      real(8), parameter :: eps = epsilon(1d0)

      PIL = acos(-1.d0)
      DEGRAD = PIL/180.d0

      FUNIT = 541
      open (UNIT=FUNIT, &
            FILE=trim(INPUTDIR)//'/'//'fort.54001', IOSTAT=IOS)

      read (FUNIT, *) NBF
      if (NBF /= AbsLayerNBF) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Error in READFORT54001(): 1")
      end if

      if (.not. allocated(AbsLayerQAMIG)) then
         allocate (AbsLayerQAMIG(AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerQFF)) then
         allocate (AbsLayerQFF(AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerQFACE)) then
         allocate (AbsLayerQFACE(AbsLayerNBF))
      end if
      if (.not. allocated(pert_vel_abslayer)) then
         allocate (pert_vel_abslayer(AbsLayerNBF))
      end if

      do KK = 1, AbsLayerNBF
         read (FUNIT, *) AbsLayerQAMIG(KK), AbsLayerQFF(KK), &
            AbsLayerQFACE(KK), IGNOREMSG

         AbsLayerQFACE(KK) = AbsLayerQFACE(KK)*DEGRAD

         if (AbsLayerQAMIG(KK) <= eps) then
        !! IF ( AbsLayerEtaAMIG(KK) .EQ. 0 ) THEN
        !!    pert_eta_abslayer(kk) = 0.D0
            pert_vel_abslayer(kk) = 0.d0
         else
            pert_vel_abslayer(kk) = 2.d0*PIL/AbsLayerQAMIG(KK)
         end if

      end do

      read (FUNIT, *) NPTMP
      if (NPTMP /= NumNodesAbsLayer(2) .and. &
          NPTMP /= NumNodesAbsLayer(3)) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Error in READFORT54001(): 1")
      end if

      if (.not. allocated(AbsLayerQxEMO)) then
         allocate (AbsLayerQxEMO(NPTMP, AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerQxEFA)) then
         allocate (AbsLayerQxEFA(NPTMP, AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerQyEMO)) then
         allocate (AbsLayerQyEMO(NPTMP, AbsLayerNBF))
      end if
      if (.not. allocated(AbsLayerQyEFA)) then
         allocate (AbsLayerQyEFA(NPTMP, AbsLayerNBF))
      end if

      do KK = 1, NPTMP
         read (FUNIT, *) KTMP

         if (KTMP /= AbsLayer_MnX_NodesID(KK)) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="Error in READFORT54001(): 2")
         end if

         do JJ = 1, AbsLayerNBF
            read (FUNIT, *) AbsLayerQxEMO(KK, JJ), &
               AbsLayerQxEFA(KK, JJ), &
               AbsLayerQyEMO(KK, JJ), &
               AbsLayerQyEFA(KK, JJ)

            ! AbsLayerQxEFA(KK,JJ) = AbsLayerQxEFA(KK,JJ)*DEGRAD
            ! AbsLayerQyEFA(KK,JJ) = AbsLayerQyEFA(KK,JJ)*DEGRAD
         end do

      end do
      AbsLayerQxEFA = AbsLayerQxEFA*DEGRAD
      AbsLayerQyEFA = AbsLayerQyEFA*DEGRAD

      close (FUNIT)

      write (16, *) "DONE READING FORT.54001"

   end subroutine READFORT54001

   subroutine READFORT2001(FIRST_CALL)
      use SIZES, only: INPUTDIR
      implicit none

      integer :: FUNIT, IOS
      logical, intent(in) :: FIRST_CALL

      integer :: KK, NPTMP

      FUNIT = 201
      NPTMP = NumNodesAbsLayer(1)
      if (FIRST_CALL) then
         open (UNIT=FUNIT, &
               FILE=trim(INPUTDIR)//'/'//'fort.2001', IOSTAT=IOS)
         read (FUNIT, *) SPTIMINC
         SPTIME1 = STATIM*86400.d0
         SPTIME2 = SPTIME1 + SPTIMINC
         allocate (ArbEta0(NumNodesAbsLayer(1)))
         allocate (ArbUU0(NumNodesAbsLayer(2)))
         allocate (ArbVV0(NumNodesAbsLayer(3)))
         allocate (ArbEta(NumNodesAbsLayer(1)))
         allocate (ArbUU(NumNodesAbsLayer(2)))
         allocate (ArbVV(NumNodesAbsLayer(3)))
         ArbEta0 = 0.d0
         ArbUU0 = 0.d0
         ArbVV0 = 0.d0
         !          ! reading in new value
         do KK = 1, NPTMP
            read (FUNIT, *) ArbEta0(KK), ArbUU0(KK), ArbVV0(KK)
         end do
         write (16, *) "READ IN FIRST TIMESTEP DATA FOR FORT.2001"
      else
         ! Updating "0" Arb vectors
         ArbEta0 = ArbEta
         ArbUU0 = ArbUU
         ArbVV0 = ArbVV
      end if
      ! Initialzing new Arb vectors
      ArbEta = 0.d0
      ArbUU = 0.d0
      ArbVV = 0.d0
      ! reading in new value
      do KK = 1, NPTMP
         read (FUNIT, *) ArbEta(KK), ArbUU(KK), ArbVV(KK)
      end do

      write (16, *) "READ IN NEW TIMESTEP DATA FOR FORT.2001"
      return
   end subroutine READFORT2001

   subroutine GETABSLAYEREXT(TimeLoc, TimeH, RampVal, &
                             GeoidOffset)
      implicit none

      real(8), intent(IN) :: TimeLoc, TimeH
      real(8), intent(IN) :: RampVal
      real(8), intent(IN) :: GeoidOffset(:)
      real(8), parameter :: EPS = epsilon(1d0)

      integer :: I, J
      real(8) :: ARG, ARGT, RFF, NCYC, DTRATIO

      !.............................DW, generating layer ..............
      if (LoadAbsLayerSigma .and. sum(NumNodesAbsLayer) > 0) then
         select case (AbsLayerType)
         case (0)
            !  Absoribing BC  !
         case (1, 2)
            !
            !  Harmonic forcing
            !
            !  \zeta in the sponge Layer
            WorkVecEta = 0.d0
            do J = 1, AbsLayerNBF
               if (pert_eta_abslayer(J) <= eps) then
                  NCYC = 0
               else
#ifdef IBM
                  NCYC = int(timeh/PERT_ETA_ABSLAYER(J), kind(0.0d0))
#else
                  NCYC = dble(int(timeh/PERT_ETA_ABSLAYER(J)))
#endif
               end if

               ARGT = AbsLayerEtaAMIG(J)*(timeh - NCYC*pert_eta_abslayer(J)) + AbsLayerEtaFACE(J)
               RFF = AbsLayerEtaFF(J)*RampVal

               do I = 1, NumNodesAbsLayer(1)
                  ARG = ARGT - AbsLayerEtaEFA(I, J)
                  WorkVecEta(I) = WorkVecEta(I) + AbsLayerEtaEMO(I, J)*RFF*cos(ARG)
               end do

            end do
            Eta2_AbsLayer(AbsLayer_Eta_NodesID) = WorkVecEta

            !  (u,v) in the sponge layer
            WorkVecUU = 0.d0
            WorkVecVV = 0.d0
            do J = 1, AbsLayerNBF
               if (pert_vel_abslayer(J) <= eps) then
                  NCYC = 0
               else
#ifdef IBM
                  NCYC = int(timeh/pert_vel_abslayer(J), kind(0.0d0))
#else
                  NCYC = dble(int(timeh/pert_vel_abslayer(J)))
#endif
               end if

               ARGT = AbsLayerQAMIG(J)*(timeh - NCYC*pert_vel_abslayer(J)) + AbsLayerQFACE(J)
               RFF = AbsLayerQFF(J)*RampVal

               do I = 1, NumNodesAbsLayer(2)
                  ARG = ARGT - AbsLayerQxEFA(I, J)
                  WorkVecUU(I) = WorkVecUU(I) + AbsLayerQxEMO(I, J)*RFF*cos(ARG)
               end do

               do I = 1, NumNodesAbsLayer(3)
                  ARG = ARGT - AbsLayerQyEFA(I, J)
                  WorkVecVV(I) = WorkVecVV(I) + AbsLayerQyEMO(I, J)*RFF*cos(ARG)
               end do

            end do
            uu2_AbsLayer(AbsLayer_MnX_NodesID) = WorkVecUU
            vv2_AbsLayer(AbsLayer_MnY_NodesID) = WorkVecVV

         case (-1)
            ! read new values in
            if (TimeLoc > SPTIME2) then
               SPTIME1 = SPTIME2
               SPTIME2 = SPTIME1 + SPTIMINC
               call READFORT2001(.false.)
            end if
            ! do the linear intepolation
            DTRATIO = (TimeLoc - SPTIME1)/SPTIMINC
            eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*(ArbEta0 + DTRATIO*(ArbEta - ArbEta0))
            uu2_AbsLayer(AbsLayer_MnX_NodesID) = RampVal*(ArbUU0 + DTRATIO*(ArbUU - ArbUU0))
            vv2_AbsLayer(AbsLayer_MnY_NodesID) = RampVal*(ArbVV0 + DTRATIO*(ArbVV - ArbVV0))
         case DEFAULT
         end select
         ! Add arbitrary to the tide if case is 2
         if (AbsLayerType == 2) then
            ! read new values in
            if (TimeLoc > SPTIME2) then
               SPTIME1 = SPTIME2
               SPTIME2 = SPTIME1 + SPTIMINC
               call READFORT2001(.false.)
            end if
            ! do the linear intepolation
            DTRATIO = (TimeLoc - SPTIME1)/SPTIMINC
            ! adding onto the tide
            eta2_AbsLayer(AbsLayer_Eta_NodesID) = &
               eta2_AbsLayer(AbsLayer_Eta_NodesID) + &
               RampVal*(ArbEta0 + DTRATIO*(ArbEta - ArbEta0))
            uu2_AbsLayer(AbsLayer_MnX_NodesID) = &
               uu2_AbsLayer(AbsLayer_MnX_NodesID) + &
               RampVal*(ArbUU0 + DTRATIO*(ArbUU - ArbUU0))
            vv2_AbsLayer(AbsLayer_MnY_NodesID) = &
               vv2_AbsLayer(AbsLayer_MnY_NodesID) + &
               RampVal*(ArbVV0 + DTRATIO*(ArbVV - ArbVV0))
         end if

         ! Consider GeoidOffset in sponge layer
         eta2_AbsLayer(AbsLayer_Eta_NodesID) = eta2_AbsLayer(AbsLayer_Eta_NodesID) + GeoidOffset(AbsLayer_Eta_NodesID)
      end if

   end subroutine GETABSLAYEREXT

   subroutine SPONGE_SHIFT_SOLN()
      implicit none

      if (NumNodesAbsLayer(1) > 0) then
         eta0_AbsLayer(AbsLayer_Eta_NodesID) = eta1_AbsLayer(AbsLayer_Eta_NodesID)
         eta1_AbsLayer(AbsLayer_Eta_NodesID) = eta2_AbsLayer(AbsLayer_Eta_NodesID)

      !!$           uu1_AbsLayer(AbsLayer_MnX_NodesID) = uu2_AbsLayer(AbsLayer_MnX_NodesID)
      !!$           vv1_AbsLayer(AbsLayer_MnY_NodesID) = vv2_AbsLayer(AbsLayer_MnY_NodesID)
      end if

      if (NumNodesAbsLayer(2) > 0) then
         uu1_AbsLayer(AbsLayer_MnX_NodesID) = uu2_AbsLayer(AbsLayer_MnX_NodesID)
      end if

      if (NumNodesAbsLayer(3) > 0) then
         vv1_AbsLayer(AbsLayer_MnY_NodesID) = vv2_AbsLayer(AbsLayer_MnY_NodesID)
      end if

      return
   end subroutine SPONGE_SHIFT_SOLN

   !$$$
   !$$$      !
   !$$$      ! Propagting a Gaussin
   !$$$      !
   !$$$      SUBROUTINE LAYERSOLTEST1( IT, TimeLoc, TimeH, RampVal )
   !$$$        IMPLICIT NONE
   !$$$
   !$$$        INTEGER:: IT
   !$$$        REAL(8):: TimeLoc, TimeH
   !$$$        REAL(8):: RampVal
   !$$$
   !$$$        REAL(8):: bb, R2, L, xc
   !$$$        REAL(8):: CVEL, LS, DY, XX, YY, CVV
   !$$$        INTEGER:: I, NIDX
   !$$$
   !$$$        bb = 74.D0
   !$$$        R2 = 250*100.D0
   !$$$        R2 = R2*R2
   !$$$
   !$$$        L = 500*1000.D0
   !$$$        DY = L/25.D0
   !$$$        XC = 2*L + 4*DY
   !$$$
   !$$$        CVEL = sqrt(G*bb)
   !$$$        CVV = sqrt(G/bb)
   !$$$
   !$$$        LS = XC - CVEL*TimeLoc
   !$$$        DO I=1, NumNodesAbsLayer(1)
   !$$$          NIDX = AbsLayer_Eta_NodesID(I)
   !$$$
   !$$$          XX = X(NIDX)
   !$$$
   !$$$          WorkVecEta(I) = exp(-(XX - LS)**2.D0/R2)
   !$$$        END DO
   !$$$        Eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*WorkVecEta
   !$$$        !
   !$$$
   !$$$        DO I=1, NumNodesAbsLayer(2)
   !$$$          NIDX = AbsLayer_Mnx_NodesID(I)
   !$$$
   !$$$          XX = X(NIDX)
   !$$$
   !$$$          WorkVecUU(I) = exp(-(XX - LS)**2.D0/R2)
   !$$$        END DO
   !$$$        uu2_AbsLayer(AbsLayer_MnX_NodesID) = -RampVal*CVV*WorkVecUU
   !$$$        vv2_AbsLayer = 0.D0
   !$$$
   !$$$        RETURN
   !$$$      END SUBROUTINE LAYERSOLTEST1
   !$$$

   !     !
   !     ! Propagting a Gaussin
   !     !
   subroutine LAYERSOLTEST1(TimeLoc, RampVal)
      implicit none

      real(8), intent(IN) :: TimeLoc
      real(8), intent(IN) :: RampVal

      real(8) :: bb, R2, L, xc, LDA, PI
      real(8) :: CVEL, LS, DY, XX, CVV
      integer :: I, NIDX

      bb = 74.d0
      R2 = 250d0*100.d0
      R2 = R2*R2

      L = 500d0*1000.d0
      DY = L/25.d0

      XC = 2d0*L + L/4d0
      LDA = 4.0d0/L

      CVEL = sqrt(G*bb)
      CVV = sqrt(G/bb)

      !c
      PI = acos(-1.d0)

      LS = XC - CVEL*TimeLoc
      do I = 1, NumNodesAbsLayer(1)
         NIDX = AbsLayer_Eta_NodesID(I)
         XX = X(NIDX) - LS

         if (XX > -1.d0 - 10) then
            WorkVecEta(I) = sin(2.d0*PI*XX*LDA)
         end if
      end do
      Eta2_AbsLayer(AbsLayer_Eta_NodesID) = RampVal*WorkVecEta

      do I = 1, NumNodesAbsLayer(2)
         NIDX = AbsLayer_Mnx_NodesID(I)

         XX = X(NIDX) - LS

         if (XX > -1.d-10) then
            WorkVecUU(I) = sin(2.d0*PI*XX*LDA)
         end if
      end do

      uu2_AbsLayer(AbsLayer_MnX_NodesID) = -RampVal*CVV*WorkVecUU
      vv2_AbsLayer = 0.d0

   end subroutine LAYERSOLTEST1

end module SPONGELAYER
