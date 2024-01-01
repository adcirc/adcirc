
module nodal_attributes_io
  use nodal_attribute_type_module

  implicit none

  private :: get_nodal_attribute

  contains

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
    class(nodal_attribute_base_t), intent(out), pointer :: ptr

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

end module nodal_attributes_io