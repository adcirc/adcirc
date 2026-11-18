.. meta::
   :description: Initial conditions in ADCIRC
   :keywords: adcirc, initial conditions

.. _initial_conditions:

Initial Conditions
==================

**Initial conditions** specify the model state at start time. They are similar
to `boundary conditions <boundary_conditions>`__, which are applied along the
spatial boundaries.

Initial water elevations in ADCIRC can be modified using `nodal
attributes <nodal_attribute>`__ like
:ref:`sea_surface_height_above_geoid <sea_surface_height_above_geoid>` or
:ref:`initial_river_elevation <initial_river_elevation>`. Boundary conditions
also affect initial water elevations. If the water elevation is modified in more
than one way, their effects are additive. Note that along elevation-specified
boundary conditions (like a tidal elevation boundary),
:ref:`sea_surface_height_above_geoid <sea_surface_height_above_geoid>` persists
throughout the simulation, and so it can be both an initial and a boundary
condition.



.. raw:: html

   <style>
   .wrap-table th, .wrap-table td {
     white-space: normal !important;
     word-wrap: break-word !important;
     max-width: 100% !important;
     overflow-wrap: break-word !important;
     hyphens: auto !important;
   }
   </style>
