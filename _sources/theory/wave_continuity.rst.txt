Wave Continuity Formulation
===========================

The Generalized Wave Continuity Equation (GWCE) is a key feature of the ADCIRC model. This formulation provides enhanced numerical stability compared to primitive continuity-momentum formulations.

Derivation of the GWCE
----------------------

The GWCE is derived by combining the primitive continuity equation with the time derivative of the primitive continuity equation. This process involves:

1. Taking the time derivative of the continuity equation
2. Adding it to the spatial derivative of the momentum equations, multiplied by a weighting parameter tau_0

The resulting equation is:

.. math::

    \frac{\partial^2 \zeta}{\partial t^2} + \tau_0 \frac{\partial \zeta}{\partial t} + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left( \frac{\partial U}{\partial t} + \tau_0 U + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left( \frac{U^2}{H} \right) + \frac{1}{R} \frac{\partial}{\partial \phi} \left( \frac{UV}{H} \right) - \frac{UV \tan \phi}{RH} - fV + \ldots \right) \\
    + \frac{1}{R \cos \phi} \frac{\partial}{\partial \phi} \left( \frac{\partial V}{\partial t} + \tau_0 V + \frac{1}{R \cos \phi} \frac{\partial}{\partial \lambda} \left( \frac{UV}{H} \right) + \frac{1}{R} \frac{\partial}{\partial \phi} \left( \frac{V^2}{H} \right) + \frac{U^2 \tan \phi}{RH} + fU + \ldots \right) \cos \phi = 0

In this equation, all terms from the momentum equations (except for the time derivative term) are included in the spatial derivatives, represented by ellipses above.

Advantages of the GWCE
----------------------

The GWCE formulation offers several advantages:

1. **Enhanced Stability**: The GWCE is less susceptible to spurious oscillations (known as "2Δx noise") that commonly affect primitive equation models.

2. **Better Mass Conservation**: The formulation better preserves mass conservation properties, especially important in coastal regions with complex bathymetry.

3. **Improved Numerical Performance**: The GWCE allows for larger timesteps and provides more accurate solutions, particularly for propagating waves.

4. **Prevention of Gravitational Oscillations**: The formulation damps out short-wavelength numerical noise without significantly affecting the physical solution.

Role of the Tau Parameter
-------------------------

The weighting parameter :math:`\tau_0` in the GWCE is a critical component that controls the numerical properties of the model:

* When :math:`\tau_0 = 0`, the GWCE reduces to a second-order wave equation
* As :math:`\tau_0 \rightarrow \infty`, the GWCE approaches the primitive continuity equation
* For optimal performance, ADCIRC typically uses a spatially varying :math:`\tau_0` that is:
  - Larger in deep water regions (approaching the primitive form)
  - Smaller in shallow regions (approaching the wave equation form)

The spatially varying :math:`\tau_0` is commonly defined as:

.. math::

    \tau_0 = 
    \begin{cases}
    \tau_{min} & \text{if } H \leq H_{crit} \\
    \tau_{min} + (\tau_{max} - \tau_{min}) \frac{H - H_{crit}}{H_{deep} - H_{crit}} & \text{if } H_{crit} < H < H_{deep} \\
    \tau_{max} & \text{if } H \geq H_{deep}
    \end{cases}

where :math:`H_{crit}` and :math:`H_{deep}` are user-defined depth thresholds, and :math:`\tau_{min}` and :math:`\tau_{max}` are minimum and maximum values for :math:`\tau_0`. 