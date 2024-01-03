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
C******************************************************************************
C     PADCIRC VERSION 46.00 xx/xx/2006
C     last changes in this file VERSION 46.00
C
C     Written for ADCIRC v46.00 by Jason G. Fleming.
C******************************************************************************
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     M O D U L E   N O D A L  A T T R I B U T E S
C-----------------------------------------------------------------------
C
C     jgf46.00 This module manages nodal attribute data, including
C     bottom friction, tau0, startdry, directional wind speed reduction,
C     and etc. Will read the Nodal Attributes File (unit 13) and
C     initialize the nodal attribute arrays.
C
C     Handling data by label rather than an integer encoding should
C     result in increased transparency as well as ease the transition to
C     HDF5/NetCDF i/o. The labels were chosen according to the
C     guidelines of the CF Standard. Creating labels according to CF
C     Standard Guidelines should enhance interoperability with other
C     simulation frameworks.
C
C     To use a nodal attribute contained in the fort.13 file, the
C     corresponding attribute name must appear in the fort.15 file.  A
C     list of nodal attributes is read in from the fort.15 file if the
C     fort.15 parameter NWP > 0. This also signals ADCIRC to look for a
C     fort.13 file.
C
C     Summary of the file format for the Nodal Attributes File:
C
C     AGRID                       ! user's comment line - should be a cross
C                                 !    reference to the grid file
C     NumOfNodes                  ! number of nodes, must match NP
C                                 !    from grid file
C     NAttr                       ! number of attributes contained in this file
C
C     do i=1, NAttr
C        AttrName(i)              ! nodal attribute name (see
C                                 !    valid names below)
C        Units(i)                 ! physical units (ft, m/s, none)
C        ValuesPerNode(i)         ! number of values at each node for
C                                 !   a particular attribute
C        DefaultVal(i)            ! default value(s) for the nodal attribute
C     end do
C
C     do i=1, NAttr
C        AttrName(i)              ! label of the attribute, again
C        NumNodesNotDefaultVal(i) ! number of nodes with non-default values
C        do j=1, NumNodesNotDefault(i)
C           n, (Attr(n,k), k=1,ValuesPerNode(i))
C        end do
C     end do
C
C
C
C     Valid labels are as follows:
C
C     ADCIRC Variable:       CF-Style Label:
C      Tau0                  "primitive_weighting_in_continuity_equation"
C      StartDry              "surface_submergence_state"
C      Fric                  "quadratic_friction_coefficient_at_sea_floor"
C      z0Land                "surface_directional_effective_roughness_length"
C                            (z0Land has ValuesPerNode = 12)
C      VCanopy               "surface_canopy_coefficient"
C      BK,BAlpha,BDelX,POAN  "bridge_pilings_friction_parameters"
C                            (bridge_pilings... has ValuesPerNode=4)
C      ManningsN             "mannings_n_at_sea_floor"
C      Chezy                 "chezy_friction_coefficient_at_sea_floor"
C      Z0b_var               "bottom_roughness_length"
C      GeoidOffset           "sea_surface_height_above_geoid"
C      EVM                   "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth"
C      EVC                   "average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth"
C      Tao0MinMax            "min_and_max_primitive_weighting_in_continuity_equation"
C      River_et_WSE          "initial_river_elevation"
C      IT_Fric               "internal_tide_friction"
C      OverlandReductionFactor  "overland_reduction_factor"
C      subgridBarrier        "subgrid_barrier"
C      CondensedNodes        "condensed_nodes"
C
C-----------------------------------------------------------------------
      MODULE NodalAttributes
      USE SIZES
      USE GLOBAL, ONLY : scratchMessage, allMessage, logMessage,
     &                    DEBUG, ECHO, INFO, WARNING, ERROR

      USE nodal_attributes_module

      IMPLICIT NONE

      INTEGER NumOfNodes    ! number of nodes listed in unit 13 file, cf. NP
      INTEGER NAttr         ! number of nodal attributes in the unit 13 file

!     The following variables are inputs from the unit 15 model param. file
      INTEGER NWP     ! number of nodal attributes to read from file
      INTEGER NoLiBF  ! nonlinear bottom friction indicator
      REAL(8) Tau0   ! primitive continuity eqn. weight
      REAL(8) Tau    ! linear friction coefficient (1/sec)
      REAL(8) CF     ! 2DDI bottom fric. coef., effect varies based on NoLiBF
      REAL(8) HBreak ! break depth for NOLIBF .eq. 2
      REAL(8) FTheta ! dimless param. for NOLIBF .eq. 2
      REAL(8) FGamma ! dimless param. for NOLIBF .eq.
      REAL(8) ESLM   ! horizontal eddy viscosity (length^2/time)
      REAL(8) ESLC   ! horizontal eddy diffusivity (length^2/time)
      INTEGER  IFLINBF! flag to turn on linear bottom friction
      INTEGER  IFNLBF ! flag to turn on nonlinear bottom friction
      INTEGER  IFHYBF ! flag to turn on hybrid bottom friction
      REAL(8) BFCdLLimit ! lower limit of quadratic bottom friction

!     Nodal attributes.
      REAL(8), ALLOCATABLE :: FRIC(:)     ! bottom friction coefficient
      REAL(8), ALLOCATABLE :: FRIC_ORIG(:)

!------------------ For absorbing layer, DW
      REAL(8), pointer:: absorblayer_sigma(:,:), absorblayer_zero(:,:)

      REAL(8), pointer:: absorblayer_sigma_eta(:,:)
      REAL(8), pointer:: absorblayer_sigma_mnx(:,:)
      REAL(8), pointer:: absorblayer_sigma_mny(:,:)

      LOGICAL LoadAbsLayerSigma
      LOGICAL FoundAbsLayerSigma

      INTEGER AbsLayerSigmaNoOfVals  
      REAL(8) AbsLayerSigmaDefVal(2)
      
      CHARACTER (LEN=80):: AbsLayerSigmaUnits

      INTEGER:: NumNodesAbsLayer(3)
      INTEGER, pointer:: AbsLayerNodesID(:)

      INTEGER, pointer:: AbsLayer_Eta_NodesID(:)
      INTEGER, pointer:: AbsLayer_MnX_NodesID(:)
      INTEGER, pointer:: AbsLayer_MnY_NodesID(:)

      INTEGER:: AbsLayerType !  0 - Absorbing BCs,  1 - Harmonic forcing (Tides)  !
                             ! -1 - Read from file                                !

! For wave generator layer of type 1
      INTEGER:: AbsLayerNBF

      ! Water elevation
      REAL(8), ALLOCATABLE:: AbsLayerEtaAMIG(:), 
     &                       AbsLayerEtaFF(:), AbsLayerEtaFACE(:)
      REAL(8), ALLOCATABLE:: AbsLayerEtaEMO(:,:), AbsLayerEtaEFA(:,:)

      ! Tentatively for Velocity or Dischage
      REAL(8), ALLOCATABLE:: AbsLayerQAMIG(:), 
     &                AbsLayerQFF(:), AbsLayerQFACE(:)
      REAL(8), ALLOCATABLE:: AbsLayerQxEMO(:,:), AbsLayerQxEFA(:,:)
      REAL(8), ALLOCATABLE:: AbsLayerQyEMO(:,:), AbsLayerQyEFA(:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONTAINS

      subroutine read_nodal_attributes_file(filename, nnodes)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nnodes
        integer :: nnodes_file
        integer :: n_attributes
        integer :: ierr
        integer :: i
        integer, parameter :: io_unit = 13
        character(len=256) :: file_header
        character(len=256) :: attribute_name
        class(nodal_attribute_base_t), pointer :: current_attribute

        open(unit=io_unit, file=filename, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
          write(*,'(A)') '[ERROR]: Error opening nodal attributes file: ', filename
          stop 1
        end if

        read(io_unit, *) file_header
        read(io_unit, *) nnodes_file
        read(io_unit, *) n_attributes

        if (nnodes_file /= nnodes) then
          write(*,'(A)') '[ERROR]: Number of nodes in nodal attributes file does not match number of nodes in mesh'
          stop 1
        end if

        if (n_attributes < 1) then
          write(*,'(A)') '[ERROR]: Number of attributes in nodal attributes file is less than 1'
          stop 1
        end if

        do i = 1, n_attributes
          read(io_unit, *) attribute_name
          backspace(io_unit)
          current_attribute => get_nodal_attribute(attribute_name)
          call current_attribute%read_header_from_file(io_unit, nnodes)
        end do

        do i = 1, n_attributes
          read(io_unit, *) attribute_name
          backspace(io_unit)
          current_attribute => get_nodal_attribute(attribute_name)
          call current_attribute%read_body_from_file(io_unit)
        end do

      end subroutine read_nodal_attributes_file

      function get_nodal_attribute(attribute_name) result(ptr)
        implicit none
        character(len=*), intent(in) :: attribute_name
        class(nodal_attribute_base_t), pointer :: ptr

        select case (attribute_name)
        case ('primitive_weighting_in_continuity_equation')
          ptr => primitive_weighting_in_continuity_equation
        case ('surface_submergence_state')
          ptr => surface_submergence_state
        case ('quadratic_friction_coefficient')
          ptr => quadratic_friction_coefficient
        case ('surface_directional_effective_roughness_length')
          ptr => surface_directional_effective_roughness_length
        case ('surface_canopy_coefficient')
          ptr => surface_canopy_coefficient
        case ('mannings_n_at_sea_floor')
          ptr => mannings_n_at_sea_floor
        case ('chezy_friction_coefficient_at_sea_floor')
          ptr => chezy_friction_coefficient_at_sea_floor
        case ('bottom_roughness_length')
          ptr => bottom_roughness_length
        case ('sea_surface_height_above_geoid')
          ptr => sea_surface_height_above_geoid
        case ('average_horizontal_eddy_viscosity_in_sea_water_wrt_depth')
          ptr => average_horizontal_eddy_viscosity_in_sea_water_wrt_depth
        case ('average_horiztonal_eddy_diffusivity_in_sea_water_wrt_depth')
          ptr => average_horiztonal_eddy_diffusivity_in_sea_water_wrt_depth
        case ('min_and_max_primitive_weighting_in_continuity_equation')
          ptr => min_and_max_primitive_weighting_in_continuity_equation
        case ('initial_river_elevation')
          ptr => initial_river_elevation
        case ('internal_tide_friction')
          ptr => internal_tide_friction
        case ('overland_reduction_factor')
          ptr => overland_reduction_factor
        case ('wave_refraction_in_swan')
          ptr => wave_refraction_in_swan
        case ('condensed_nodes')
          ptr => condensed_nodes
        case ('subgrid_barrier')
          ptr => subgrid_barrier
        case ('elemental_slope_limiter')
          ptr => elemental_slope_limiter
        case default
          write(*,'(A)') '[ERROR]: Unknown nodal attribute: ', attribute_name
          stop 1
        end select

      end function get_nodal_attribute


      subroutine initialize_nodal_attributes_module()
        use mesh, only: np
        implicit none

        !...Some of the nodal attributes are always active, so we initialize them here
        !   and set their active flags to .true.

        call primitive_weighting_in_continuity_equation%init(
     &     np, "primitive_weighting_in_continuity_equation",
     &     "none", 1, 0.005D0)
        call average_horizontal_eddy_viscosity_in_sea_water_wrt_depth%init(
     &     np, "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth",
     &     "m^2/s", 1, 2.0D0)

      end subroutine initialize_nodal_attributes_module


!     ----------------------------------------------------------------
!     S U B R O U T I N E     R E A D  N O D A L  A T T R
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to read the nodal attributes file (unit 13).
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
      USE MESH, ONLY : NP
      USE SIZES, ONLY: INPUTDIR
      IMPLICIT NONE

      integer, intent(in) :: NScreen, ScreenUnit, MyProc, NAbOut
      character(len=256)  :: fname

      fname = TRIM(INPUTDIR)//'/'//'fort.13'
      write(fname, '(a,i0)') fname, MyProc

      call read_nodal_attributes_file(fname, NP)

!     ----------------------------------------------------------------
      END SUBROUTINE ReadNodalAttr
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  2 D  B O T T O M  F R I C T I O N
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to apply 2D bottom friction from turbulent
!     viscous effects as well as bridge pilings. This is used in the
!     time stepping loop.
!
!     ----------------------------------------------------------------
!
!     sb46.28sb02 Lower limit of Cd was added as an argument.
!     jgf47.04 Argument for lower limit of Cd was removed; this value
!     is now specified by the user in fort.15.
!
!     ----------------------------------------------------------------
      SUBROUTINE Apply2DBottomFriction(UU1, VV1, DP, ETA2, G,
     &     IFNLFA, NP, TK)
      USE GLOBAL, ONLY: MLD, DTDP
      USE subgrid, ONLY: level0, level1
     &   , cmfVertETA2, cfVertETA2, subgridVertList
     &   , wetDepthVertETA2
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                   ! number of nodes in grid
      REAL(8), intent(in), dimension(NP) :: UU1  ! x-dir velocities
      REAL(8), intent(in), dimension(NP) :: VV1  ! y-dir velocities
      REAL(8), intent(in), dimension(NP) :: DP   ! bathymetric depths
      REAL(8), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
      REAL(8), intent(in) :: G                   ! gravitational constant
      INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
      REAL(8), intent(inout), dimension(NP) :: TK! depth avg. fric.
C
      REAL(8) :: UV1   ! velocity magnitude (speed)
      REAL(8) :: H1    ! total depth
      REAL(8) :: Fr
      REAL(8) :: FricBP
      REAL(8) :: BK    ! BK(1) is pier shape factor
      REAL(8) :: BALPHA! BALPHA(2) is constriction fraction
      REAL(8) :: BDELX ! BDELX(3) is effective delx
      REAL(8) :: MLDrat ! MLD ratio to depth
      REAL(8) :: HT, WCH, ManN
      INTEGER :: I
      type(nodal_attribute_vector_real_t), pointer :: n => mannings_n_at_sea_floor
      type(nodal_attribute_subgrid_barrier_t), pointer :: sgb => subgrid_barrier
      type(nodal_attribute_vector_real_t), pointer :: chezy_coefficient => chezy_friction_coefficient_at_sea_floor
      type(nodal_attribute_matrix_real_t), pointer :: bp => bridge_pilings

      !...Reset the friction to original values
      !   when subgrid barriers might be used
      if(sgb%active.and.n%active.and.sgb%n_values==2)then
        where(n%values>0D0.and.sgb%overtopping)
          FRIC = sgb%friction
        elsewhere
          FRIC = n%values
        endwhere
      endif

!     Step 0. Convert Manning's N to Cd, if necessary.
      IF (n%active) THEN
         DO I=1, NP
            IF(level0)THEN
               IF(subgridVertList(I)==1)THEN
                  IF(level1)THEN
                     FRIC(I) = cmfVertETA2(I)
                  ELSE
                     FRIC(I) = cfVertETA2(I)
                  ENDIF
               ELSE
                  FRIC(I) = g*n%values(I)**2.d0
     &               /((DP(I)+IFNLFA*ETA2(I))**(1.d0/3.d0))
               ENDIF
            ELSE
               HT = DP(I)+IFNLFA*ETA2(I)
               IF(sgb%active)THEN
                   IF(sgb%overtopping(I).AND.sgb%n_values==2)THEN
                      FRIC(I)=MIN(BFCdLLimit, g*sgb%values(I,2)**2.d0
     &                      /( ( HT )**(1.d0/3.d0)) )
                   ELSE
                       IF (n%values(I)<=0d0) cycle ! Can use Cd where man <= 0
                       FRIC(I)=MIN(BFCdLLimit, g*n%values(I)**2.d0
     &                      /( ( HT )**(1.d0/3.d0) )) ! sb46.28sb02
                   ENDIF
               ELSE
                   IF (n%values(I)<=0d0) cycle ! Can use Cd where man <= 0
                   IF (n%values(I)>1.d0) THEN
                      WCH = floor(n%values(I))
                      ManN = n%values(I) - WCH
                      FRIC(I)=MIN(BFCdLLimit, g*((1.D0+2.D0*HT/WCH)**(2.D0/3.D0)*ManN)**2.d0
     &                  /( HT**(1.d0/3.d0) ))
                   ELSE
                      FRIC(I)=MIN(BFCdLLimit, g*n%values(I)**2.d0
     &                  /( ( HT )**(1.d0/3.d0) )) ! sb46.28sb02
                   ENDIF
               ENDIF
            ENDIF
         ENDDO
      ELSEIF (chezy_coefficient%active) THEN
         FRIC = g/chezy_coefficient%values**2.0d0
      ENDIF

      IF (NOLIBF==3) THEN
!     ... Combined Cd for bottom friction and momentum dispersion
         DO I = 1,NP
            UV1 = SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1  = DP(I) + IFNLFA*ETA2(I)
            IF (HBREAK>0d0) then
!..............CPB 3/9/2023: moved calculation of Cd as a power function
!              to a separate subroutine "apply2DMomentumDisp"
!              ... Cd as a power function [a(x^b-1)] of x = MLD
               FRIC(I) = CF ! + HBREAK*(max(4d-3,MLD(I))**FTHETA - 1d0)
            ELSE
!              ... Cd directly from the 3D model (MLD here)
               FRIC(I) = FRIC(I) + MLD(I)
            ENDIF
            TK(I) = MIN(FRIC(I) * UV1 / H1, 2.0/DTDP)
         ENDDO

      ELSE

!     .. Step 1. Apply friction arising from turbulent viscous interaction
!        with the sea floor.
         DO I=1, NP
            UV1 = SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            IF(level0)THEN
               IF(subgridVertList(I)==1)THEN
                  H1 = wetDepthVertETA2(I)
               ELSE
                  H1 = DP(I) + IFNLFA*ETA2(I)
               ENDIF
            ELSE
               H1  = DP(I) + IFNLFA*ETA2(I)
            ENDIF
            TK(I) = MIN(FRIC(I)*
     &                     (IFLINBF +                                          ! linear
     &                     (UV1/H1) * (IFNLBF +                               ! nonlinear
     &                     IFHYBF*(1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)) ), ! hybrid 
     &                  2.0/DTDP)
         END DO

      ENDIF

!    Step 2. Apply friction arising from flow interaction with bridge
!     pilings, if required.
      IF (bp%active) THEN
         DO I=1, NP
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1=DP(I)+IFNLFA*ETA2(I)
            Fr=UV1*UV1/(G*H1)
            BK = bp%values(I,1)
            BALPHA = bp%values(I,2)
            BDELX = bp%values(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0)
     &           *(BALPHA+15.d0*BALPHA**4)
            TK(I) = MIN(TK(I)+FricBP*UV1/H1, 2.0/DTDP)
         END DO
      ENDIF

!     ----------------------------------------------------------------
      END SUBROUTINE Apply2DBottomFriction
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  3 D  B O T T O M  F R I C T I O N
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to apply 3D bottom friction from turbulent
!     viscous effects as well as bridge pilings. This is used in the
!     time stepping loop.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
      SUBROUTINE Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G,
     &     IFNLFA, NP, TK, NFEN, Z0B)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP, NFEN             ! number of nodes in grid Horizontal and Vertical
      COMPLEX(8), intent(in), dimension(NP,NFEN) :: Q  ! x-dir velocities
      REAL(8), intent(in), dimension(NFEN) :: SIGMA  ! x-dir velocities
      REAL(8), intent(in), dimension(NP) :: DP   ! bathymetric depths
      REAL(8), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
      REAL(8), intent(in) :: G                   ! gravitational constant
      INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
      REAL(8), intent(inout), dimension(NP) :: TK! depth avg. fric.
      REAL(8), intent(in) :: Z0B

      INTEGER :: NH
      INTEGER :: I
      REAL(8) :: Z0B1  ! velocity magnitude (speed)
      REAL(8) :: UV1   ! velocity magnitude (speed)
      REAL(8) :: H1    ! total depth
      REAL(8) :: Fr
      REAL(8) :: FricBP
      REAL(8) :: BK    ! BK(1) is pier shape factor
      REAL(8) :: BALPHA! BALPHA(2) is constriction fraction
      REAL(8) :: BDELX ! BDELX(3) is effective delx
      type(nodal_attribute_vector_real_t), pointer :: n => mannings_n_at_sea_floor
      type(nodal_attribute_matrix_real_t), pointer :: bp => bridge_pilings

! Determine the bottom roughness length either from fort.15, from Manning's n
! or as read in from nodal attributes
      DO NH=1,NP
         H1=DP(NH)+IFNLFA*ETA2(NH)
         IF (bottom_roughness_length%active) THEN
            Z0B1 = bottom_roughness_length%values(NH)
         ELSEIF (n%active) THEN
            Z0B1 = ( H1 )* exp(-(1.0D0+
     &             ( (0.41D0*( H1 )**(1.0D0/6.0D0) )/
     &                            (n%values(NH)*sqrt(g)) ) ))
         ELSE
            Z0B1 = Z0B
         ENDIF

         FRIC(NH)= (1.D0 / ( (1.D0/0.41D0) *
     &   LOG((ABS( ( ( SIGMA(2)-SIGMA(1) )/2.d0 ) *(H1) ) + Z0B1 )/Z0B1)
     &       ) )**2.D0

         TK(NH)= FRIC(NH) * ABS(Q(NH,1))

         IF (bp%active) THEN
            Fr=ABS(Q(NH,1))*ABS(Q(NH,1))/(G*H1)
            BK = bp%values(I,1)
            BALPHA = bp%values(I,2)
            BDELX = bp%values(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0)
     &           *(BALPHA+15.d0*BALPHA**4)
            TK(I)=TK(I)+FricBP*ABS(Q(NH,1))/H1
         ENDIF
      ENDDO

!     ----------------------------------------------------------------
      END SUBROUTINE Apply3DBottomFriction
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  2 D  I N T E R N A L  W A V E  D R A G
!     ----------------------------------------------------------------
!
!     WJP March 30 2016: Adjusted for TKM input matrix and inserting
!     the scalar version in here too
!     WJP 02.24.2018 put into v53.dev
!     WJP 07.24.2018 updated for baroclinic coupling type
!     ----------------------------------------------------------------
      SUBROUTINE Apply2DInternalWaveDrag(NP, TK, TKM, U, V, DP,
     & IT)      
      USE GLOBAL, ONLY: CBaroclinic, CORIF, DTDP, NB, NM, IDEN,
     &                  IFSPROTS, RAMP
      USE MESH, ONLY: DRVMAP2DSPVEC, UVECTMP, VVECTMP, RVELF
      USE INTERNALTIDE, ONLY: UBar, VBar, MunkHPFilter
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                   ! number of nodes in grid
      REAL(8), intent(in), dimension(:) :: TK   ! scalar friction
      REAL(8), intent(inout), dimension(:,:) :: TKM
      INTEGER,  intent(in), optional :: IT          ! temporal interations
      REAL(8), intent(in), dimension(NP),optional :: U    ! x-dir velocities
      REAL(8), intent(in), dimension(NP),optional :: V    ! y-dir velocities
      REAL(8), intent(in), dimension(NP),optional :: DP   ! depth
      REAL(8), allocatable :: Hx(:), Hy(:)
      INTEGER  :: I, J 
      REAL(8) :: omega = 1.4052d-4 ![rad/s] M2 tidal frequency
      REAL(8) :: hr = 3600.0D0 ! means one hour in seconds
      REAL(8) :: SmallV = 1d-10 ! A small value
      REAL(8) :: Fit, alpha2, gamma2, slope2
      REAL(8) :: Nmdiff, Nbdiff, Fdiff, adju, adjv, adjuv
      logical, save ::  first_entry = .TRUE.
      ! CPB: for tidal velocities
      REAL(8) :: cosAMIGT, sinAMIGT, DissTidal, DissTotal, jjwgamma,
     &            Cxx, Cyy, Cxy
      INTEGER :: NCYC
      REAL(8),DIMENSION(:),ALLOCATABLE :: Udmy, Vdmy
      type(nodal_attribute_matrix_real_t), pointer :: itfric => internal_tide_friction

!     Compute TKM by adding wave drag contribution to
!     current value of the bottom friction, TK
      IF (itfric%active) THEN
         IF ( FIRST_ENTRY ) THEN
            FIRST_ENTRY = .FALSE.
            IF ( IFSPROTS==1 ) THEN
!              Rotate the IT_Fric first time step if mesh rotation used
!              (unnecessary for IT_FricNoOfVals == 1)
               IF ( itfric%n_values==3 ) THEN
                  allocate(Hx(NP),Hy(NP))
                  DO I = 1,NP
                     if (CBaroclinic.and.abs(IDEN)>=5) then
                        Hx(I) = itfric%values(I,2)
                        Hy(I) = itfric%values(I,3)
                     else
                        Hx(I) = sqrt(itfric%values(I,1))
                        Hy(I) = sign(sqrt(itfric%values(I,2)),itfric%values(I,3))
                     endif
                     UVECTMP(I) = Hx(I); 
                     VVECTMP(I) = Hy(I);
                  ENDDO
                  CALL DRVMAP2DSPVEC( Hx, Hy,
     &                UVECTMP, VVECTMP, NP, FWD = .TRUE. ) ;
                  DO I = 1,NP
                     if (CBaroclinic.and.abs(IDEN)>=5) then
                        itfric%values(I,2) = Hx(I)
                        itfric%values(I,3) = Hy(I)
                     else
                        itfric%values(I,1) = Hx(I)*Hx(I)
                        itfric%values(I,2) = Hy(I)*Hy(I)
                        itfric%values(I,3) = Hx(I)*Hy(I)
                     endif
                  ENDDO           
                  deallocate(Hx,Hy)
               ELSEIF ( itfric%n_values==5 ) THEN
                  ! need do something else
               ENDIF
            ENDIF
         ENDIF
         IF (CBaroclinic.and.abs(IDEN)>=5) THEN
!           2DDI ADCIRC with baroclinic coupling
            if (abs(IDEN)/=7) then
               CALL MunkHPFilter(U,V,IT)
            endif
            IF (itfric%n_values == 5) THEN
!           Use the "Nonlocal" Nycander (2005) type tensor wave drag
!           1st column is coefficient, 2nd and 3rd columns are the Hx and Hy slopes,
!           4th and 5th columns are the Jx and Jy slopes
               DO I = 1,NP
                  IF (omega > NB(I).or.omega < CORIF(I)) THEN
                     ! No freely propagating internal tides
                     TKM(1:2,I) = TK(I)
                     TKM(3,I)   = 0.0d0
                  ELSE 
                     ! Have freely propagating internal tides
                     ! Get coefficient
                     Fit = itfric%values(I,1)*NB(I)/DP(I)*
     &                     sqrt(1d0 - CORIF(I)**2/omega**2)
                     ! Check criticality of slope
                     alpha2 = (omega**2 - CORIF(I)**2) / 
     &                        (NB(I)**2 - omega**2)
                     gamma2 = (itfric%values(I,2)**2 + itfric%values(I,3)**2)/alpha2
                     ! Divide by criticality if > 1 
                     Fit = Fit/max(1d0,gamma2)                   
                     TKM(1,I) = TK(I) + 2*Fit*itfric%values(I,2)*itfric%values(I,4)
                     TKM(2,I) = TK(I) + 2*Fit*itfric%values(I,3)*itfric%values(I,5)
                     TKM(3,I) = Fit*(itfric%values(I,2)*itfric%values(I,5) +
     &                               itfric%values(I,3)*itfric%values(I,4))
                  ENDIF 
               ENDDO
            ELSEIF (itfric%n_values == 3.or.itfric%n_values == 1) THEN
!           For the "Local" Lyard et al. (2004) type tensor wave drag
!           1st column is coefficient, 2nd and 3rd columns are the Hx and Hy slopes
!           For the scalar Zaron & Egbert (2006) type
!           1st column is coefficient, 2nd column is the squared hypotenuse slope
               DO I = 1,NP
                  adju = 1d0; adjv = 1d0;adjuv = 1d0;
                  NBdiff = (NB(I)**2d0 - omega**2d0)
                  NMdiff = (NM(I)**2d0 - omega**2d0)
                  Fdiff  = (omega**2d0 - CORIF(I)**2d0)
                  IF (NBdiff<=SmallV.or.Nmdiff<=SmallV.or.
     &                Fdiff<=SmallV) THEN
                     ! No freely propagating internal tides
                     TKM(1:2,I) = TK(I)
                     TKM(3,I)   = 0.0d0
                  ELSE 
                     ! Have freely propagating internal tides
                     ! Get the full coefficient
                     IF (NBdiff*NMdiff<0d0) THEN
                        Fit = 0d0
                     ELSE
                        Fit = itfric%values(I,1)*sqrt(NBdiff*Nmdiff)/omega
                     ENDIF
                     ! Check criticality of slope
                     alpha2 = Fdiff / NBdiff
                     IF (alpha2<0d0) Fit = 0d0
                     if (itfric%n_values == 3) then
                        slope2 = itfric%values(I,2)**2d0 + itfric%values(I,3)**2d0
                     else
                        slope2 = itfric%values(I,2)
                     endif
                     gamma2 = slope2 / alpha2
                     ! Divide coefficient by criticality if > 1
                     ! (saturates out at high criticalities) 
                     Fit = Fit / max(1d0,gamma2)                   
                     ! Calculate adjustment based on tidal residual
                     if (abs(IDEN)/=7) then
                        ! Calculate the adjustment factors based on the delayed moving mean
                        IF (ABS(U(I))>1d-3.AND.
     &                      ABS(V(I))>1d-3) THEN
                           adju = UBar(I)/U(I)
                           adjv = VBar(I)/V(I)
                        ENDIF
                        adjuv = 0.5d0*(adju + adjv)
                     endif
                     ! Implement the friction with possible adjustment
                     IF (itfric%n_values == 3) then
                        TKM(1,I) = TK(I) + adju*Fit*itfric%values(I,2)**2d0
                        TKM(2,I) = TK(I) + adjv*Fit*itfric%values(I,3)**2d0
                        TKM(3,I) = adjuv*Fit*itfric%values(I,2)*itfric%values(I,3)
                     else
                        TKM(1:2,I) = TK(I) + adjuv*Fit*itfric%values(I,2)
                        TKM(3,I)   = 0.0d0
                     endif
                  ENDIF 
               ENDDO
            endif
         ELSE
!           Normal 2DDI ADCIRC without baroclinic coupling
            IF (itfric%n_values == 3) THEN
!              Use the "Nonlocal" Nycander (2005) or "Local" Lyard et al. (2004)
!              type tensor wave drag
!              Using three columns of data for the different segments of the matrix
               DO I=1,NP 
                  TKM(1,I) = TK(I) + itfric%values(I,1)
                  TKM(2,I) = TK(I) + itfric%values(I,2)
                  TKM(3,I) = itfric%values(I,3)
               END DO
            ELSEIF (itfric%n_values == 1) THEN
!              For the scalar type Zaron & Egbert (2006) type
               DO I=1,NP 
                  TKM(1,I) = TK(I) + itfric%values(I,1)
                  TKM(2,I) = TKM(1,I)
                  TKM(3,I) = 0.0d0
               ENDDO
            ENDIF
         ENDIF
      ELSE 
!       In case we did not load ITFric
        TKM(1,:) = TK
        TKM(2,:) = TK
        TKM(3,:) = 0.0d0
      ENDIF

!     ----------------------------------------------------------------
      END SUBROUTINE Apply2DInternalWaveDrag
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  2 D  M O M E N T U M D I S P
!     ----------------------------------------------------------------
!     CPB 3/9/2023: This subroutine applies momentum dispersion as a
!     power function of mixed-layer depth calculated from a
!     tide-omitting OGCM. It applies a Cd calculated from the power
!     function fit to the non-tidal portion of velocity calculated using
!     a 25-hour lagged average filter.
!     ----------------------------------------------------------------
      SUBROUTINE Apply2DMomentumDisp(UU1, VV1, DP, ETA2,
     &     IFNLFA, NP, TKM)
      USE GLOBAL, ONLY: MLD, DTDP
      USE INTERNALTIDE, ONLY: UBar, VBar
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                     !> number of nodes in grid
      REAL(8), intent(in), dimension(NP) :: UU1     !> x-dir velocities
      REAL(8), intent(in), dimension(NP) :: VV1     !> y-dir velocities
      REAL(8), intent(in), dimension(NP) :: DP      !> bathymetric depths
      REAL(8), intent(in), dimension(NP) :: ETA2    !> water surf. elevations
      INTEGER, intent(in) :: IFNLFA                 !> nonlin. finite amp. flag
      REAL(8), intent(inout), dimension(:,:) :: TKM !> bl/itfric matrix

      REAL(8) :: UVBar, UV1  ! non-tidal velocity magnitude
      REAL(8) :: H1    ! total depth
      REAL(8) :: Cd    ! Cd estimated as a power function
      REAL(8) :: adju, adjv ! adjustment to apply only to nontidal vel
      REAL(8) :: a, b, c ! for ensuring positive semidefinite
      INTEGER :: I
      
      ! if nolibf =/=3 then simply return
      IF (NOLIBF/=3.OR.HBREAK<=0d0) RETURN

!     ... Combined Cd for bottom friction and momentum dispersion
      DO I = 1,NP
         IF (DP(I)<1000d0) CYCLE
         UV1 = SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
         H1  = DP(I) + IFNLFA*ETA2(I)
         UVBar = SQRT(UBar(I)*UBar(I)+VBar(I)*VBar(I))
         H1  = DP(I) + IFNLFA*ETA2(I)
!        ... Cd as a power function [a(x^b-1)] of x = MLD
         Cd = HBREAK*(max(4d-3,MLD(I))**FTHETA - 1d0)
         ! adjust to put in TKM matrix
         IF ( ABS(UU1(I))>1d-3.AND.ABS(VV1(I))>1d-3 ) THEN
            adju = UBar(I)/UU1(I)
            adjv = VBar(I)/VV1(I)
         ELSE
            adju = 1d0
            adjv = 1d0
         ENDIF
         TKM(1,I) = TKM(1,I) + Cd*adju*UVBar/H1
         TKM(2,I) = TKM(2,I) + Cd*adjv*UVBar/H1
         ! ensure matrix is positive semi-definite
         a = TKM(1,I)
         b = TKM(3,I)
         c = TKM(2,I)
         TKM(1,I) = a
         TKM(3,I) = b
         TKM(2,I) = c
      ENDDO

!     ----------------------------------------------------------------
      END SUBROUTINE Apply2DMomentumDisp
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!                      S U B R O U T I N E
!         R E A D  L E G A C Y  S T A R T  D R Y  F I L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load up the legacy startdry file (unit
!     12). This is just a cut-and-paste from the section of the
!     READ_INPUT subroutine that did the same thing. This subroutine is
!     never called. It is vestigial and listed here purely as reference
!     material.
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyStartDryFile(NP, NScreen, ScreenUnit,
     &     MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

      INTEGER :: JKI          ! node number from file
      INTEGER :: NE2          ! number of elements, according to fort.12 file
      INTEGER :: NP2          ! number of nodes, according to fort.12 file
      INTEGER :: I

      CHARACTER(len=80) AGRID2 ! users comment/description line
      REAL(8) DUM1, DUM2  ! data that we want to skip

      type(nodal_attribute_vector_logical_t), pointer :: startdry => surface_submergence_state

      OPEN(12,FILE=TRIM(INPUTDIR)//'/'//'fort.12')
C
C...  READ STARTDRY INFORMATION FROM UNIT 12
      READ(12,'(A80)') AGRID2
      WRITE(16,2038) AGRID2
2038  FORMAT(5X,'STARTDRY FILE IDENTIFICATION : ',A80,/)
      READ(12,*) NE2,NP2

!...  CHECK THAT NE2 AND NP2 MATCH WITH GRID FILE
!      IF((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
       IF(NP2/=NP) THEN
         IF(NSCREEN/=0.AND.MYPROC==0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ',
     &        /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL EXIT(1)
      ENDIF

!...  READ IN STARTDRY CODE VALUES
      DO I=1,NP
         READ(12,*) JKI,DUM1,DUM2,startdry%values(JKI)
         IF(JKI/=I) THEN
            IF(NSCREEN/=0.AND.MYPROC==0) WRITE(ScreenUnit,99805)
            WRITE(16,99805)
99805       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',
     &           //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &           'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
         ENDIF
      END DO

!...  CLOSE UNIT 12 FILE
      CLOSE(12)

      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyStartDryFile
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!                       S U B R O U T I N E
!       R E A D  L E G A C Y  B O T T O M  F R I C T I O N  F I L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load up the legacy Spatially Varying
!     Friction Coefficient File (unit 21). This is just a cut-and-paste
!     from the section of the READ_INPUT subroutine that did the same
!     thing. This subroutine is never called. It is vestigial and listed
!     here purely as reference material.
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyBottomFrictionFile(NP, NScreen, ScreenUnit,
     &     MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

      CHARACTER(len=80)   :: AFRIC  ! user's comment/description line
      INTEGER             :: NHG    ! node number from file
      INTEGER             :: I

      OPEN(21,FILE=TRIM(INPUTDIR)//'/'//'fort.21')
      READ(21,'(A80)') AFRIC
      DO I=1,NP
         READ(21,*) NHG,FRIC(NHG)
         IF(NHG/=I) THEN
            IF(NSCREEN/=0.AND.MYPROC==0) WRITE(ScreenUnit,99803)
            WRITE(16,99803)
99803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',
     &           'INPUT ERROR  !!!!!!!!!',//,1X,
     &           'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',
     &           /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',//,1X,
     &           '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            CALL EXIT(1)
         ENDIF
      END DO
      WRITE(16,3601) AFRIC
 3601 FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A80,/)
      IF(NABOUT/=1) THEN
         WRITE(16,2080)
 2080    FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION FRIC',5X,/)
         DO I=1,NP
            WRITE(16,2087) I,FRIC(I)
 2087       FORMAT(7X,I6,6X,E17.10)
         END DO
      ELSE
         WRITE(16,3504)
 3504    FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',
     &        /,6X,' IN UNIT 21 INPUT FILE')
      ENDIF

      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyBottomFrictionFile
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!        S U B R O U T I N E   N A _ T E R M I N A T E
!     ----------------------------------------------------------------
!     Provide clean termination when an error occurs.
!     ----------------------------------------------------------------
      SUBROUTINE na_terminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL, ONLY : setMessageSource, unsetMessageSource, allMessage,
     &  INFO, WARNING
      IMPLICIT NONE

      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE

      call setMessageSource("terminate")
#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1) 

#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!     ----------------------------------------------------------------
      END SUBROUTINE na_terminate
!     ----------------------------------------------------------------

!-----------------------------------------------------------------------
      END MODULE NodalAttributes
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


