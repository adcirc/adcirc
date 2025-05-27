.. _initial_river_elevation_description:

Initial River Elevation
=======================

The :ref:`initial_river_elevation` nodal attribute is used to set the initial water surface elevation at specified nodes. This attribute is functionally identical to the :ref:`sea_surface_height_above_geoid` attribute, except that any values assigned to elevation-specified boundary conditions in this attribute **do not** persist after the start of the simulation. As a result, this attribute is more appropriate for initializing rivers. If both attributes are specified, then the assigned elevation is their sum. Currently, this attribute is only applied if the value supplied at a given node is above zero (i.e. depth less than zero). From cstart.F:

   .. code-block:: 

      where (River_et_WSE.GT.0.d0)
         eta2 = eta2 + River_et_WSE
      end where



ADCIRC assumes by default that vertices with negative depths will be dry when the simulation starts. This is, of course, not the case for an inland river whose bed is above mean sea level. This nodal attribute is used in those cases to provide the initial water surface elevation of the river at cold start, and is typically used in conjunction with a flux or elevation :ref:`boundary condition <boundary_conditions>` at the inland boundary. See also :ref:`initial conditions <initial_conditions>`.   
