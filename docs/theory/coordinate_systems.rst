Coordinate Systems
==================

ADCIRC can operate in different coordinate systems depending on the application requirements and domain scale. This section describes the coordinate transformations and implementations in the model.

Spherical Coordinates
---------------------

For large domains where the Earth's curvature is significant, ADCIRC uses a spherical coordinate system:

* Longitude (:math:`\lambda`): Angular distance east or west from the Prime Meridian
* Latitude (:math:`\phi`): Angular distance north or south from the equator
* Vertical coordinate (:math:`z`): Distance above or below the reference geoid

The governing equations presented in the previous sections are formulated in spherical coordinates for global and regional applications.

Cartesian Coordinates
---------------------

For smaller domains where Earth's curvature effects are negligible, ADCIRC can operate in a Cartesian coordinate system:

* :math:`x`: Eastward distance
* :math:`y`: Northward distance
* :math:`z`: Vertical distance from reference level

The continuity equation in Cartesian coordinates is:

.. math::

    \frac{\partial \zeta}{\partial t} + \frac{\partial (UH)}{\partial x} + \frac{\partial (VH)}{\partial y} = 0

The momentum equations in Cartesian coordinates are:

.. math::

    \frac{\partial U}{\partial t} + \frac{\partial}{\partial x} \left( \frac{U^2}{H} \right) + \frac{\partial}{\partial y} \left( \frac{UV}{H} \right) - fV + g \frac{\partial \zeta}{\partial x} + \frac{\tau_{sx}}{\rho_0 H} - \frac{\tau_{bx}}{\rho_0 H} - \frac{1}{\rho_0 H} \frac{\partial p_s}{\partial x} + \frac{\partial}{\partial x} \left[ \frac{2 N_x}{H} \frac{\partial U}{\partial x} \right] + \frac{\partial}{\partial y} \left[ \frac{N_y}{H} \left( \frac{\partial U}{\partial y} + \frac{\partial V}{\partial x} \right) \right] = 0

.. math::

    \frac{\partial V}{\partial t} + \frac{\partial}{\partial x} \left( \frac{UV}{H} \right) + \frac{\partial}{\partial y} \left( \frac{V^2}{H} \right) + fU + g \frac{\partial \zeta}{\partial y} + \frac{\tau_{sy}}{\rho_0 H} - \frac{\tau_{by}}{\rho_0 H} - \frac{1}{\rho_0 H} \frac{\partial p_s}{\partial y} + \frac{\partial}{\partial x} \left[ \frac{N_x}{H} \left( \frac{\partial U}{\partial y} + \frac{\partial V}{\partial x} \right) \right] + \frac{\partial}{\partial y} \left[ \frac{2 N_y}{H} \frac{\partial V}{\partial y} \right] = 0

Coordinate Transformations
--------------------------

For regional applications, ADCIRC often employs a conformal map projection to transform between spherical and projected Cartesian coordinates. The most commonly used projections are:

1. **Mercator Projection**:
   * Preserves angles (conformal)
   * Distorts areas, especially at high latitudes
   * Transformation equations:
   
   .. math::
   
       x = R \lambda
       
   .. math::
   
       y = R \ln\left[\tan\left(\frac{\pi}{4} + \frac{\phi}{2}\right)\right]

2. **Lambert Conformal Conic**:
   * Preserves angles
   * Minimizes distortion in mid-latitudes
   * Commonly used for regional modeling in mid-latitudes

3. **Stereographic Projection**:
   * Conformal projection
   * Useful for polar regions
   * Minimal distortion near the projection center

CPP Coordinate System
---------------------

For some specialized applications, ADCIRC uses a Cartesian Coordinates with Polar Projection (CPP) system. This hybrid approach:

* Maintains the simplicity of Cartesian equations
* Accounts for Earth's curvature through carefully designed projections
* Applies correction factors to the Coriolis and other terms

The CPP coordinate system is defined by:

.. math::

    x = (R + h) \cos \phi \sin(\lambda - \lambda_0)

.. math::

    y = (R + h) [\sin\phi\cos\phi_0 - \cos\phi\sin\phi_0\cos(\lambda - \lambda_0)]

where :math:`(\lambda_0, \phi_0)` are the coordinates of the projection origin.

Vertical Coordinate Systems
---------------------------

ADCIRC employs several vertical coordinate systems:

1. **Sigma Coordinates**:
   * Terrain-following coordinates that map the water column to a uniform layer
   * :math:`\sigma = \frac{z - \zeta}{\zeta + h}` ranges from 0 (surface) to -1 (bottom)
   * Advantages include natural handling of bathymetry

2. **Z-level Coordinates**:
   * Fixed vertical levels
   * Typically used in deeper waters
   * Provides better representation of stratification

3. **Hybrid Systems**:
   * Combination of sigma and z-level approaches
   * Optimizes advantages of both systems 