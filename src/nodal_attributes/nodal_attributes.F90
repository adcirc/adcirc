
module nodal_attributes_module

  use nodal_attribute_type_module
  use nodal_attribute_wind_module
  use nodal_attribute_tau0_module

  implicit none

  type(nodal_attribute_tau0_t), target :: primitive_weighting_in_continuity_equation
  type(nodal_attribute_vector_logical_t), target :: surface_submergence_state
  type(nodal_attribute_vector_real_t), target :: quadratic_friction_coefficient
  type(nodal_attribute_vector_real_t), target :: mannings_n_at_sea_floor
  type(nodal_attribute_vector_real_t), target :: chezy_friction_coefficient_at_sea_floor
  type(nodal_attribute_vector_real_t), target :: bottom_roughness_length
  type(nodal_attribute_vector_real_t), target :: sea_surface_height_above_geoid
  type(nodal_attribute_vector_real_t), target :: average_horizontal_eddy_viscosity_in_sea_water_wrt_depth
  type(nodal_attribute_vector_real_t), target :: average_horiztonal_eddy_diffusivity_in_sea_water_wrt_depth
  type(nodal_attribute_vector_real_t), target :: min_and_max_primitive_weighting_in_continuity_equation
  type(nodal_attribute_vector_real_t), target :: initial_river_elevation
  type(nodal_attribute_matrix_real_t), target :: internal_tide_friction
  type(nodal_attribute_vector_real_t), target :: overland_reduction_factor
  type(nodal_attribute_vector_logical_t), target :: wave_refraction_in_swan
  type(nodal_attribute_matrix_real_t), target :: bridge_pilings
  type(nodal_attribute_advect_local_t), target :: advection_state
  type(nodal_attribute_condensed_nodes_t), target :: condensed_nodes
  type(nodal_attribute_subgrid_barrier_t), target :: subgrid_barrier
  type(nodal_attribute_slope_limiter_t), target :: elemental_slope_limiter
  type(nodal_attribute_surface_roughness_t), target :: surface_directional_effective_roughness_length
  type(nodal_attribute_canopy_t), target :: surface_canopy_coefficient

end module nodal_attributes_module