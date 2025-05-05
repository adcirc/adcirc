Baroclinic Terms
================

ADCIRC can be operated in baroclinic mode, allowing it to simulate density-driven flows resulting from temperature and salinity variations. This section describes the theoretical foundation for baroclinic processes in ADCIRC.

Density Variation
-----------------

In baroclinic mode, water density varies spatially and temporally due to variations in temperature and salinity. The density is computed using an equation of state:

.. math::

    \rho = \rho(T, S, p)

where:

* :math:`\rho` = water density
* :math:`T` = temperature
* :math:`S` = salinity
* :math:`p` = pressure

ADCIRC typically employs a simplified form of the UNESCO equation of state or the linear equation of state:

.. math::

    \rho = \rho_0 [1 - \alpha_T (T - T_0) + \beta_S (S - S_0)]

where :math:`\alpha_T` is the thermal expansion coefficient and :math:`\beta_S` is the saline contraction coefficient.

Baroclinic Pressure Gradient
----------------------------

The baroclinic pressure gradient is the primary driving force for density-driven flows. In the momentum equations, this appears as:

.. math::

    \frac{1}{\rho_0} \frac{\partial}{\partial \lambda} \int_{-h}^{\zeta} \int_{z'}^{\zeta} \frac{\partial \rho}{\partial \lambda'} g dz dz'

.. math::

    \frac{1}{\rho_0} \frac{\partial}{\partial \phi} \int_{-h}^{\zeta} \int_{z'}^{\zeta} \frac{\partial \rho}{\partial \phi'} g dz dz'

These terms represent the horizontal pressure gradient due to horizontal density variations. A mathematically equivalent but numerically advantageous form is:

.. math::

    \frac{g}{\rho_0} \int_{-h}^{\zeta} \left( z - \zeta \right) \frac{\partial \rho}{\partial \lambda} dz

.. math::

    \frac{g}{\rho_0} \int_{-h}^{\zeta} \left( z - \zeta \right) \frac{\partial \rho}{\partial \phi} dz

Transport Equations
-------------------

In baroclinic mode, ADCIRC solves transport equations for temperature and salinity:

.. math::

    \frac{\partial (HT)}{\partial t} + \frac{\partial (UHT)}{\partial x} + \frac{\partial (VHT)}{\partial y} = \frac{\partial}{\partial x} \left( H K_h \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y} \left( H K_h \frac{\partial T}{\partial y} \right) + H Q_T

.. math::

    \frac{\partial (HS)}{\partial t} + \frac{\partial (UHS)}{\partial x} + \frac{\partial (VHS)}{\partial y} = \frac{\partial}{\partial x} \left( H K_h \frac{\partial S}{\partial x} \right) + \frac{\partial}{\partial y} \left( H K_h \frac{\partial S}{\partial y} \right) + H Q_S

where:

* :math:`K_h` = horizontal diffusion coefficient
* :math:`Q_T` = temperature source/sink term (e.g., heating, cooling)
* :math:`Q_S` = salinity source/sink term (e.g., evaporation, precipitation, river inflow)

Three-Dimensional Formulation
-----------------------------

For fully three-dimensional baroclinic simulations, ADCIRC uses a sigma-coordinate system with the vertical domain divided into layers. The 3D momentum equations include additional terms for vertical advection and diffusion:

.. math::

    \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + w \frac{\partial u}{\partial z} - fv + \ldots + \frac{\partial}{\partial z} \left( K_v \frac{\partial u}{\partial z} \right) = 0

.. math::

    \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + w \frac{\partial v}{\partial z} + fu + \ldots + \frac{\partial}{\partial z} \left( K_v \frac{\partial v}{\partial z} \right) = 0

where:

* :math:`u, v, w` = velocity components in 3D
* :math:`K_v` = vertical eddy viscosity

The vertical velocity (:math:`w`) is diagnosed from the continuity equation.

The 3D transport equations for temperature and salinity are:

.. math::

    \frac{\partial T}{\partial t} + u \frac{\partial T}{\partial x} + v \frac{\partial T}{\partial y} + w \frac{\partial T}{\partial z} = \frac{\partial}{\partial z} \left( K_T \frac{\partial T}{\partial z} \right) + \ldots

.. math::

    \frac{\partial S}{\partial t} + u \frac{\partial S}{\partial x} + v \frac{\partial S}{\partial y} + w \frac{\partial S}{\partial z} = \frac{\partial}{\partial z} \left( K_S \frac{\partial S}{\partial z} \right) + \ldots

where :math:`K_T` and :math:`K_S` are vertical diffusion coefficients for temperature and salinity, respectively.

Surface Heat Flux
-----------------

The heat flux at the air-sea interface is parametrized by:

.. math::

    Q_{net} = Q_{sw} - Q_{lw} - Q_{sen} - Q_{lat}

where:

* :math:`Q_{sw}` = short-wave (solar) radiation
* :math:`Q_{lw}` = long-wave (thermal) radiation
* :math:`Q_{sen}` = sensible heat flux
* :math:`Q_{lat}` = latent heat flux (evaporation)

These fluxes are computed based on meteorological data such as air temperature, humidity, wind speed, and cloud cover.

Numerical Considerations for Baroclinic Simulations
---------------------------------------------------

Baroclinic simulations introduce several numerical challenges:

1. **Stable Stratification**: In stably stratified regions, vertical mixing is suppressed, requiring adequate vertical resolution.

2. **Internal Waves**: Baroclinic modes include internal waves, which have higher frequencies than external modes, potentially requiring smaller time steps.

3. **Mode Splitting**: ADCIRC can use a mode-splitting approach, where external (barotropic) and internal (baroclinic) modes are solved with different time steps to enhance computational efficiency.

4. **Pressure Gradient Errors**: Numerical errors in computing the baroclinic pressure gradient, particularly in regions of steep bathymetry, are mitigated through specialized algorithms. 