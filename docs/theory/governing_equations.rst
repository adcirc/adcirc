Governing Equations
===================

The ADCIRC model solves the shallow water equations, a form of the Navier-Stokes equations with traditional hydrostatic pressure and Boussinesq approximations. This section describes the fundamental equations that form the basis of the model.

Continuity Equation
-------------------

The depth-integrated continuity equation in spherical coordinates is:

.. math::

    \frac{\partial \zeta}{\partial t} + \frac{1}{R \cos \phi} \left[ \frac{\partial}{\partial \lambda} \left( U H \right) + \frac{\partial}{\partial \phi} \left( V H \cos \phi \right) \right] = 0

where:

* :math:`\zeta` = free surface elevation relative to the geoid
* :math:`t` = time
* :math:`R` = radius of the Earth
* :math:`\phi` = latitude
* :math:`\lambda` = longitude
* :math:`H = \zeta + h` = total water column depth
* :math:`h` = bathymetric depth relative to the geoid
* :math:`U` = depth-integrated velocity in longitude direction
* :math:`V` = depth-integrated velocity in latitude direction

Momentum Equations
------------------

The depth-integrated momentum equations in spherical coordinates are:

.. math::

    \frac{\partial U}{\partial t} + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left( \frac{U^2}{H} \right) + \frac{1}{R} \frac{\partial}{\partial \phi} \left( \frac{UV}{H} \right) - \frac{UV \tan \phi}{R H} - fV \\
    + \frac{g}{R \cos \phi} \frac{\partial \zeta}{\partial \lambda} + \frac{\tau_{s\lambda}}{\rho_0 H} - \frac{\tau_{b\lambda}}{\rho_0 H} - \frac{1}{\rho_0 H} \left[ \frac{\partial}{\partial \lambda} \left( p_s \right) - \frac{\partial}{\partial \lambda} \int_{-h}^{\zeta} \int_{z'}^{\zeta} \frac{\partial \rho}{\partial \lambda'} g dz dz' \right] \\
    + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left[ \frac{2 N_\lambda}{H} \frac{\partial U}{\partial \lambda} \right] + \frac{1}{R} \frac{\partial}{\partial \phi} \left[ \frac{N_\phi}{H} \left( \frac{\partial U}{\partial \phi} + \frac{\partial V}{\partial \lambda} - V \tan \phi \right) \right] = 0

.. math::

    \frac{\partial V}{\partial t} + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left( \frac{UV}{H} \right) + \frac{1}{R} \frac{\partial}{\partial \phi} \left( \frac{V^2}{H} \right) + \frac{U^2 \tan \phi}{R H} + fU \\
    + \frac{g}{R} \frac{\partial \zeta}{\partial \phi} + \frac{\tau_{s\phi}}{\rho_0 H} - \frac{\tau_{b\phi}}{\rho_0 H} - \frac{1}{\rho_0 H} \left[ \frac{\partial}{\partial \phi} \left( p_s \right) - \frac{\partial}{\partial \phi} \int_{-h}^{\zeta} \int_{z'}^{\zeta} \frac{\partial \rho}{\partial \phi'} g dz dz' \right] \\
    + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left[ \frac{N_\lambda}{H} \left( \frac{\partial U}{\partial \phi} + \frac{\partial V}{\partial \lambda} - V \tan \phi \right) \right] + \frac{1}{R} \frac{\partial}{\partial \phi} \left[ \frac{2 N_\phi}{H} \frac{\partial V}{\partial \phi} \right] = 0

where:

* :math:`f = 2 \Omega \sin \phi` = Coriolis parameter
* :math:`\Omega` = angular speed of the Earth
* :math:`g` = gravitational acceleration
* :math:`\rho_0` = reference density of water
* :math:`\rho` = perturbation density
* :math:`p_s` = atmospheric pressure at the free surface
* :math:`\tau_{s\lambda}, \tau_{s\phi}` = surface stress components
* :math:`\tau_{b\lambda}, \tau_{b\phi}` = bottom stress components
* :math:`N_\lambda, N_\phi` = lateral eddy viscosity coefficients

Bottom Stress Formulation
-------------------------

Bottom stress is parametrized using a quadratic friction law:

.. math::

    \frac{\tau_{b\lambda}}{\rho_0} = \frac{C_f U \sqrt{U^2 + V^2}}{H}

.. math::

    \frac{\tau_{b\phi}}{\rho_0} = \frac{C_f V \sqrt{U^2 + V^2}}{H}

where :math:`C_f` is the bottom friction coefficient, which can be specified as a constant or calculated from:

.. math::

    C_f = \frac{g n^2}{H^{1/3}}

for Manning's formulation, or:

.. math::

    C_f = \max \left[ C_{fmin}, \frac{\kappa^2}{\ln^2 \left( \frac{H}{z_0} \right)} \right]

for the logarithmic formulation, where:

* :math:`n` = Manning's roughness coefficient
* :math:`\kappa` = von Karman constant (≈ 0.4)
* :math:`z_0` = bottom roughness length
* :math:`C_{fmin}` = minimum bottom friction coefficient 