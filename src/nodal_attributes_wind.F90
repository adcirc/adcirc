MODULE NODAL_ATTRIBUTE_WIND_MODULE

  USE NODAL_ATTRIBUTE_BASE_TYPE_MODULE

  IMPLICIT NONE

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_CANOPY_T
  ! Description: Type for nodal attribute related to surface canopy
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_VECTOR_REAL_T) :: NODAL_ATTRIBUTE_CANOPY_T
  contains
    procedure, pass(self) :: apply => surface_canopy_wind_application
  end type NODAL_ATTRIBUTE_CANOPY_T

  !--------------------------------------------------------------------
  ! Type: NODAL_ATTRIBUTE_SURFACE_ROUGHNESS_T
  ! Description: Type for nodal attribute related to surface roughness
  !--------------------------------------------------------------------
  type, extends(NODAL_ATTRIBUTE_MATRIX_REAL_T) :: NODAL_ATTRIBUTE_SURFACE_ROUGHNESS_T
  contains
    procedure, pass(self) :: apply => surface_roughness_application
  end type NODAL_ATTRIBUTE_SURFACE_ROUGHNESS_T

  contains

    !--------------------------------------------------------------------
    ! Subroutine: surface_canopy_wind_application
    ! Description: Apply the surface canopy to the wind vector
    !--------------------------------------------------------------------
    subroutine surface_canopy_wind_application(self, node_number, wind_x, wind_y)
      implicit none
      class(NODAL_ATTRIBUTE_CANOPY_T), intent(inout) :: self
      integer, intent(in) :: node_number
      real(8), intent(inout) :: wind_x, wind_y

      if(self%active.eqv..false.) return

      wind_x = wind_x * self%values(node_number)
      wind_y = wind_y * self%values(node_number)

    end subroutine surface_canopy_wind_application

    !--------------------------------------------------------------------
    ! Subroutine: surface_roughness_application
    ! Description: Apply the surface roughness to the wind vector
    ! ----------------------------------------------------------------
    !>    @brief Subroutine to calculate the land wind reduction factor
    !>    based on a table of directional wind drag values. Originally
    !>    written into the hstart.F file by jjw in jjw-42.06j. This is used
    !>    in hstart.F and timestep.F.
    !>
    !>    The methology comes from:
    !>    Simiu, E. and R. Scanlan, 1986: Wind Effects on Structures. Wiley Interscience, 604 pp.
    !>
    !>    Additionally, the z0 coefficient supplied by the user is modified
    !>    by subtracting a factor of depth/30 to account for inundation
    !>
    !>    Values are linearly interpolated between heading bins to the actual
    !>    wind direction.
    !>
    !>    @param[in] NodeNumber index of node under consideration
    !>    @param[in] WindDragCo current wind drag coefficient being applied
    !>    @param[in] WindMag    current wind speed magnitude at the specified node
    !>    @param[in] BathymetricDepth bathymetric depth used to calculate total depth
    !>    @param[in] Elevation water surface elevation above datum (eta2)
    !>    @param[in,out] WindX wind velocity, x-component
    !>    @param[in,out] WindY wind velocity, y-component
    !--------------------------------------------------------------------
    subroutine surface_roughness_application(self, node_number, wind_drag_coefficient,&
            wind_magnitude, bathymetric_depth, water_surface_elevation, &
            overland_reduction_factor, wind_x, wind_y)
      use GLOBAL,only:h0
      use CONSTANTS,only:RAD2DEG,G
      use wind_drag_module, only: WindDrag
      implicit none
      class(NODAL_ATTRIBUTE_SURFACE_ROUGHNESS_T), intent(inout) :: self
      integer, intent(in)    :: node_number
      real(8), intent(inout) :: wind_drag_coefficient
      real(8), intent(inout) :: wind_magnitude
      real(8), intent(in)    :: bathymetric_depth
      real(8), intent(in)    :: water_surface_elevation
      type(nodal_attribute_vector_real_t), intent(in) :: overland_reduction_factor
      real(8), intent(inout) :: wind_x, wind_y

      integer :: i

      real(8),parameter      :: CutoffDepth = 999999.d0          !< Depth above which reduction is not applied
      real(8),parameter,dimension(12) :: z0Angles=(/(i*30,i=0,5),(i*30-180,i=0,5)/)  !< centers of directional bins

      integer                :: idir        !< code for wind direction bins to interpolate between
      integer                :: idir2       !< code for wind direction bins to interpolate between
      real(8)                :: z0m         !< marine roughness coefficient based on Garratt's formula
      real(8)                :: angle       !< direction wind is coming from
      real(8)                :: z0l         !< drag for a particular node, for particular direction
      real(8)                :: TotalDepth  !< bathymetric depth + sea surface elevation
      real(8)                :: fr          !< land wind reduction factor
      real(8), parameter     :: eps = EPSILON(1.0D0) !< machine epsilon

      if(self%active.eqv..false.) return

      ! if windspeed is zero, exit now
      if(abs(wind_x)<eps.and.abs(wind_y)<eps)return

      ! compute direction  that the wind is going to
      angle=atan2(wind_y,wind_x)*RAD2DEG
      idir=0
      idir2=0
      if(angle<-150d0)then
        idir=7
        idir2=8
      elseif(angle<-120d0)then
        idir=8
        idir2=9
      elseif(angle<-90d0)then
        idir=9
        idir2=10
      elseif(angle<-60d0)then
        idir=10
        idir2=11
      elseif(angle<-30d0)then
        idir=11
        idir2=12
      elseif(angle<0d0)then
        idir=12
        idir2=1
      elseif(angle<30d0)then
        idir=1
        idir2=2
      elseif(angle<60d0)then
        idir=2
        idir2=3
      elseif(angle<90d0)then
        idir=3
        idir2=4
      elseif(angle<120d0)then
        idir=4
        idir2=5
      elseif(angle<150d0)then
        idir=5
        idir2=6
      elseif(angle<=180d0)then
        idir=6
        idir2=7
      endif

      ! compute marine roughness coefficient based on Garratt's formula
      z0m=(0.018d0/G)*wind_drag_coefficient*wind_magnitude**2.d0

      !tga 2020-04 updated code to linearly interpolate z0 values,
      !code previously did nearest neighbor (i.e. binned) z0 values.
      !Define roughness by linearly (in angle-z0 space) interpolating
      !between binned values.
      !The 30d0 here comes from assuming the directional bins are 30
      !degrees apiece.
      z0l=(self%values(node_number,idir2)-self%values(node_number,idir))/30d0 * &
              (angle-z0Angles(idir))+self%values(node_number,idir)

      ! apply overland flooding correction
      TotalDepth = bathymetric_depth + water_surface_elevation
      if( (TotalDepth>2D0*h0).and.(bathymetric_depth<CutOffDepth)) then
        z0l=z0l-TotalDepth*overland_reduction_factor%values(node_number)
      endif

      ! compute land wind reduction factor
      ! Reduction factor is bounded to not exceed 1, i.e. assumes the
      ! land roughness is never less than the water.
      if(z0l>z0m) then
        fr=min(1.0D0,(z0l/z0m)**0.0706d0 * log(10.d0/z0l) / log(10.d0/z0m))
      else
        fr=1.000d0
      endif

      ! adjust time interpolated wind field
      wind_x = fr*wind_x
      wind_y = fr*wind_y

      ! update the input parameters
      wind_magnitude = sqrt(wind_x**2 + wind_y**2)
      wind_drag_coefficient = WindDrag(wind_magnitude, node_number)

    end subroutine surface_roughness_application

END MODULE NODAL_ATTRIBUTE_WIND_MODULE