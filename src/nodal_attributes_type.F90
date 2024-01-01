!**********************************************************************
! MODULE: NODAL_ATTRIBUTE_TYPE_MODULE
! DESCRIPTION: Module containing types and procedures related to nodal
!              attributes in the ADCIRC model.
! WRITTEN BY:  Zachary Cobell
! CONTACT INFO: zcobell@thewaterinstitute.org
!
! This module re-implements the logic which was buried within older versions
! of nodalattr.F into a more modern and modular format. The goal is to
! use an extensible structure for nodal attributes which will allow for more
! explicit variable typing and usage. The base classes for this module are
! defined in nodal_attribute_base_type.F90.
!
!**********************************************************************

MODULE NODAL_ATTRIBUTE_TYPE_MODULE

  USE NODAL_ATTRIBUTE_BASE_TYPE_MODULE

  IMPLICIT NONE

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_SLOPE_LIMITER_T
  ! Description: Type for nodal attribute related to slope limiter
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_VECTOR_REAL_T) :: NODAL_ATTRIBUTE_SLOPE_LIMITER_T
    logical, allocatable :: elemental_slope_limiter_active(:) !< Flag indicating whether the elemental slope limiter is active
    logical, allocatable :: elemental_slope_limiter_max_exceeded(:) !< Flag indicating whether the elemental slope limiter has been exceeded
    real(8), allocatable :: elemental_slope_limiter_grad_max(:) !< Maximum gradient for the elemental slope limiter
    integer, allocatable :: count(:) !< Count for the elemental slope limiter

  contains
    final :: slope_limiter_destructor
  end type NODAL_ATTRIBUTE_SLOPE_LIMITER_T

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_SUBGRID_BARRIER_T
  ! Description: Type for nodal attribute related to subgrid barrier
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_MATRIX_REAL_T) :: NODAL_ATTRIBUTE_SUBGRID_BARRIER_T
    logical, allocatable :: overtopping(:) !< Flag indicating whether the subgrid barrier is overtopping
    real(8), allocatable :: friction(:)    !< Friction for the subgrid barrier

  contains
    final :: subgrid_barrier_destructor
  end type NODAL_ATTRIBUTE_SUBGRID_BARRIER_T

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_CONDENSED_NODES_T
  ! Description: Type for nodal attribute related to condensed nodes
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_MATRIX_INT_T) :: NODAL_ATTRIBUTE_CONDENSED_NODES_T
    integer, allocatable :: n_condensed_nodes(:) !< Number of condensed nodes for the condensed node
    integer, allocatable :: list_condensed_nodes(:,:) !< List of condensed nodes for the condensed node
    integer, allocatable :: nnodes_list_condensed_nodes(:) !< Number of condensed nodes for the condensed node
    real(8), allocatable :: CSIICN(:) !< X component of the unit vector pointing from the condensed node to the node
    real(8), allocatable :: SIIICN(:) !< Y component of the unit vector pointing from the condensed node to the node
    integer :: NListCondensedNodes = 0 !< Number of condensed nodes for the condensed node

  contains
    procedure, pass(self) :: prep => condensed_nodes_prep
    final :: condensed_nodes_destructor
  end type NODAL_ATTRIBUTE_CONDENSED_NODES_T

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_ADVECT_LOCAL_T
  ! Description: Type for nodal attribute related to advect local
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_VECTOR_REAL_T) :: NODAL_ATTRIBUTE_ADVECT_LOCAL_T
    contains
      procedure, pass(self) :: advect_local
  end type NODAL_ATTRIBUTE_ADVECT_LOCAL_T

  !--------------------------------------------------------------------
  ! Module Procedures
  !--------------------------------------------------------------------

contains

  !--------------------------------------------------------------------
  ! Subroutine: slope_limiter_destructor
  ! Description: Destructor for NODAL_ATTRIBUTE_SLOPE_LIMITER_T
  !--------------------------------------------------------------------
  subroutine slope_limiter_destructor(self)
    implicit none
    type(NODAL_ATTRIBUTE_SLOPE_LIMITER_T), intent(inout) :: self
    if(allocated(self%elemental_slope_limiter_active)) deallocate(self%elemental_slope_limiter_active)
    if(allocated(self%elemental_slope_limiter_max_exceeded)) deallocate(self%elemental_slope_limiter_max_exceeded)
    if(allocated(self%elemental_slope_limiter_grad_max)) deallocate(self%elemental_slope_limiter_grad_max)
    if(allocated(self%values)) deallocate(self%values)
  end subroutine slope_limiter_destructor

  !--------------------------------------------------------------------
  ! Subroutine: subgrid_barrier_destructor
  ! Description: Destructor for NODAL_ATTRIBUTE_SUBGRID_BARRIER_T
  !--------------------------------------------------------------------
  subroutine subgrid_barrier_destructor(self)
    implicit none
    type(NODAL_ATTRIBUTE_SUBGRID_BARRIER_T), intent(inout) :: self
    if(allocated(self%overtopping)) deallocate(self%overtopping)
    if(allocated(self%friction)) deallocate(self%friction)
    if(allocated(self%values)) deallocate(self%values)
  end subroutine subgrid_barrier_destructor

  !--------------------------------------------------------------------
  ! Subroutine: condensed_nodes_destructor
  ! Description: Destructor for NODAL_ATTRIBUTE_CONDENSED_NODES_T
  !--------------------------------------------------------------------
  subroutine condensed_nodes_destructor(self)
    implicit none
    type(NODAL_ATTRIBUTE_CONDENSED_NODES_T), intent(inout) :: self
    if(allocated(self%n_condensed_nodes)) deallocate(self%n_condensed_nodes)
    if(allocated(self%list_condensed_nodes)) deallocate(self%list_condensed_nodes)
    if(allocated(self%nnodes_list_condensed_nodes)) deallocate(self%nnodes_list_condensed_nodes)
    if(allocated(self%CSIICN)) deallocate(self%CSIICN)
    if(allocated(self%SIIICN)) deallocate(self%SIIICN)
    if(allocated(self%values)) deallocate(self%values)
  end subroutine condensed_nodes_destructor

  !--------------------------------------------------------------------
  ! Subroutine: condensed_nodes_prep
  ! Description: Prepare condensed nodes in NODAL_ATTRIBUTE_CONDENSED_NODES_T
  !--------------------------------------------------------------------
  subroutine condensed_nodes_prep(self)
    USE SIZES, ONLY : MNVEL
    USE MESH, ONLY : NP, LBArray_Pointer, X, Y
    IMPLICIT NONE
    CLASS(NODAL_ATTRIBUTE_CONDENSED_NODES_T), INTENT(INOUT) :: self
    INTEGER :: I, J, K, II, I1, I2, J1, J2
    INTEGER :: JMAX
    REAL(8) :: X1, Y1, X2, Y2, DX, DY, LEN
    LOGICAL :: ON_THE_SAME_BOUNDARY

    self%NListCondensedNodes = 0
    JMAX = 0
    DO I = 1, NP
      IF (self%values(I, 1) /= 0) THEN
        self%NListCondensedNodes = self%NListCondensedNodes + 1
        DO J = 1, self%n_values
          IF (self%values(I, J) == 0) EXIT
          IF (JMAX < J) JMAX = J
        ENDDO
      ENDIF
    ENDDO

    ALLOCATE(self%n_condensed_nodes(NP))
    ALLOCATE(self%list_condensed_nodes(self%NListCondensedNodes, JMAX + 1))
    ALLOCATE(self%nnodes_list_condensed_nodes(self%NListCondensedNodes))

    self%nnodes_list_condensed_nodes(:) = 0

    K = 0
    DO I = 1, NP
      IF (self%values(I, 1) > 0) THEN
        K = K + 1
        self%list_condensed_nodes(K, 1) = I
        DO J = 1, self%n_values
          IF (self%values(I, J) <= 0) EXIT
          II = self%values(I, J)
          self%n_condensed_nodes(I) = J + 1
          self%list_condensed_nodes(K, J + 1) = II
          self%nnodes_list_condensed_nodes(K) = J + 1
          self%values(II, 1) = -I ! To tell that node II is one of condensed nodes
        ENDDO
        DO J = 1, self%n_values
          IF (self%values(I, J) <= 0) EXIT
          II = self%values(I, J)
          self%n_condensed_nodes(II) = self%n_condensed_nodes(I)
        ENDDO
      ENDIF
    ENDDO

    ALLOCATE(self%CSIICN(MNVEL), self%SIIICN(MNVEL))
    self%CSIICN = 0.D0
    self%SIIICN = 0.D0

    DO I = 1, NP
      IF ((self%values(I, 1) > 0) .AND. (self%n_condensed_nodes(I) == 2)) THEN
        I1 = I
        I2 = self%values(I, 1)

        J1 = LBArray_Pointer(I1)
        J2 = LBArray_Pointer(I2)
        IF (J1 < 1 .OR. J2 < 1) CYCLE
        IF (J1 == (J2 + 1) .OR. J1 == (J2 - 1)) THEN ! Whether Nodes J1 and J2 are on the same boundary.
          ON_THE_SAME_BOUNDARY = .TRUE.
        ELSE
          ON_THE_SAME_BOUNDARY = .FALSE.
        END IF
        IF (ON_THE_SAME_BOUNDARY) CYCLE

        X1 = X(I1)
        Y1 = Y(I1)
        X2 = X(I2)
        Y2 = Y(I2)
        DX = X1 - X2
        DY = Y1 - Y2
        LEN = SQRT(DX * DX + DY * DY)
        DX = DX / LEN
        DY = DY / LEN

        self%CSIICN(J1) = DX
        self%SIIICN(J1) = DY

        self%CSIICN(J2) = DX
        self%SIIICN(J2) = DY
      END IF
    END DO

  END SUBROUTINE condensed_nodes_prep

  !--------------------------------------------------------------------
  ! Subroutine: advect_local
  ! Description: Apply the localized advection scheme
  !--------------------------------------------------------------------
  subroutine advect_local(self, ie)
    use GLOBAL, only : IFNLCAT, IFNLCATE, IFNLCT, IFNLCTE
    use MESH, only : NM, DP
    implicit none
    class(NODAL_ATTRIBUTE_ADVECT_LOCAL_T), intent(inout) :: self
    integer, intent(in) :: ie

    integer  :: NM1, NM2, NM3

    if(self%active.eqv..false.) return

    NM1=NM(IE,1)
    NM2=NM(IE,2)
    NM3=NM(IE,3)
    if ((DP(NM1)>=self%values(NM1)).and.&
        (DP(NM2)>=self%values(NM2)).and.&
        (DP(NM3)>=self%values(NM3))) then
      IFNLCT = IFNLCTE
      IFNLCAT = IFNLCATE
    else
      IFNLCT = 0
      IFNLCAT = 0
    endif
  end subroutine advect_local

end module NODAL_ATTRIBUTE_TYPE_MODULE
