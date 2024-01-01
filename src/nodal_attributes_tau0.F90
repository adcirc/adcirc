MODULE NODAL_ATTRIBUTE_TAU0_MODULE

  USE NODAL_ATTRIBUTE_BASE_TYPE_MODULE

  IMPLICIT NONE

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_TAU0_T
  ! Description: Type for nodal attribute related to tau0
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_VECTOR_REAL_T) :: NODAL_ATTRIBUTE_TAU0_T
    real(8), allocatable :: tau0_last(:)         !< Last value of tau0 for the nodal attribute
    real(8), allocatable :: tau0_temp(:)         !< Temporary value of tau0 for the nodal attribute
    real(8), allocatable :: tau0_base(:)         !< Base value of tau0 for the nodal attribute
    real(8), allocatable :: tau0_minmax(:,:)     !< Minimum/maximum value of tau0 for the nodal attribute
    real(8)              :: alpha_tau0 = 0.25d0  !< Alpha value for tau0 for the nodal attribute
    real(8)              :: tau0_full_domain_min
    real(8)              :: tau0_full_domain_max
    logical              :: HighResTimeVaryingTau0 = .false.
    logical              :: FullDomainTimeVaryingTau0 = .false.
    logical              :: TimeAveragedTau0 = .false.
    logical              :: BackLoadedTimeAveragedTau0 = .false.
    logical              :: outputTau0 = .false.

  contains
    final :: tau0_destructor
    procedure, pass(self) :: compute_time_varying_tau0
    procedure, pass(self) :: tau0_compute

  end type NODAL_ATTRIBUTE_TAU0_T

  contains

    !--------------------------------------------------------------------
    ! Subroutine: tau0_destructor
    ! Description: Destructor for NODAL_ATTRIBUTE_TAU0_T
    !--------------------------------------------------------------------
    subroutine tau0_destructor(self)
      implicit none
      type(NODAL_ATTRIBUTE_TAU0_T), intent(inout) :: self
      if(allocated(self%tau0_last)) deallocate(self%tau0_last)
      if(allocated(self%tau0_temp)) deallocate(self%tau0_temp)
      if(allocated(self%tau0_base)) deallocate(self%tau0_base)
      if(allocated(self%tau0_minmax)) deallocate(self%tau0_minmax)
      if(allocated(self%values)) deallocate(self%values)
    end subroutine tau0_destructor

    !--------------------------------------------------------------------
    ! Subroutine: compute_time_varying_tau0
    ! Description: Compute the time varying tau0
    !--------------------------------------------------------------------
    !
    ! jgf47.08 Subroutine to calculate a new tau0 value. Called from
    ! GWCE_New in timestep.F each time the GWCE matrix is reset (i.e.,
    ! upon startup and whenever wetting and/or drying occurs in any
    ! subdomain. Based on Casey050711.
    !
    !----------------------------------------------------------------
    SUBROUTINE compute_time_varying_tau0(self, TK, NNeigh, NeiTab, NP)
      IMPLICIT NONE
      CLASS(NODAL_ATTRIBUTE_TAU0_T), INTENT(INOUT) :: self
      REAL(8), intent(in) :: TK(:)       ! bottom friction
      INTEGER, intent(in) :: NNeigh(:)   ! number of neighbor nodes
      INTEGER, intent(in) :: NeiTab(:,:) ! table of neighbor nodes
      INTEGER, intent(in) :: NP          ! number of nodes in the domain
      REAL(8)             :: TempSum     ! sum of tau0temp values around a particular node
      INTEGER             :: I, J        ! loop counters

      ! Casey 050711 : Made changes for averaged variable Tau0.
      ! jjw46.39.sb01 :  "high/low LIMITED" variable G.
      ! jgf47.30: Distinction between fulldomain and hi res only
      IF ( self%FullDomainTimeVaryingTau0 ) THEN
        self%tau0_temp = MAX(self%Tau0_MinMax(:,2),&
                MIN(self%Tau0_MinMax(:,1), self%Tau0_MinMax(:,1) + 1.5D0*TK(:)))
      ENDIF
      IF ( self%HighResTimeVaryingTau0 ) THEN
        WHERE(self%TAU0_Base<0.025D0)
          self%Tau0_Temp = self%Tau0_Base
        ELSEWHERE
          self%Tau0_Temp = MIN(0.2D0, self%Tau0_Base + 1.5D0*TK(:))
        ENDWHERE
      ENDIF
      ! smoothing
      DO I=1, NP
        TempSum = 0.0D0
        DO J=1,NNeigh(I)
          TempSum = TempSum + self%Tau0_Temp(NeiTab(I,J))
        ENDDO
        self%values(I) = TempSum / NNeigh(I)
      ENDDO

      !     jgf47.33 Perform time averaging of tau0 if requested.
      IF (self%TimeAveragedTau0) THEN
        self%values = 0.5d0*self%values + 0.5d0*self%tau0_last
        self%tau0_Last = self%values
      ENDIF

      !     jgf48.42 Perform backloaded time averaging of tau0 if requested.
      IF (self%BackLoadedTimeAveragedTau0) THEN
        self%values = self%alpha_tau0*self%values+(1.d0-self%alpha_tau0)*self%tau0_Last
        self%tau0_last = self%values
      ENDIF
!     ----------------------------------------------------------------
    END SUBROUTINE compute_time_varying_tau0
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!      F U N C T I O N   T A U 0  N O D A L  V A L U E
!     ----------------------------------------------------------------
!
!     jgf46.27 Function to calculate tau0 based on the scheme selection
!     and the depth. This assumes that Scheme is negative.
!
!     ----------------------------------------------------------------
    REAL(8) FUNCTION tau0_compute(self, scheme, Depth)
      IMPLICIT NONE
      CLASS(NODAL_ATTRIBUTE_TAU0_T), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: scheme
      REAL(8) :: Depth

      IF (Scheme==-2) THEN
        !     Smoothly varying tau0 with depth.
        IF(Depth>=200D0) THEN
          tau0_compute=0.005D0
        ELSEIF((Depth<200D0).AND.(Depth>=1D0)) THEN
          tau0_compute=1D0/Depth
        ELSEIF(Depth<1D0) THEN
          tau0_compute=1.0D0
        ENDIF
      ELSE
        ! Abrupt variation in tau0 with depth.
        IF(Depth<=10D0) tau0_compute=0.020d0
        IF(Depth>10D0) tau0_compute=0.005d0
      ENDIF
!     ----------------------------------------------------------------
    END FUNCTION tau0_compute
!     ----------------------------------------------------------------

END MODULE NODAL_ATTRIBUTE_TAU0_MODULE