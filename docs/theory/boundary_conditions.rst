Boundary Conditions
===================

ADCIRC implements various boundary conditions to represent physical constraints and forcing at domain boundaries. These boundary conditions are essential for accurately simulating coastal and ocean processes.

External Boundary Conditions
----------------------------

External boundaries are the outer boundaries of the computational domain. ADCIRC supports these types of external boundary conditions:

Elevation Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At an elevation-specified boundary, the water surface elevation is prescribed as a function of time:

.. math::

    \zeta(t) = \sum_{k=1}^{N_{tides}} f_k \alpha_k \cos(\omega_k t - \phi_k)

where:

* :math:`N_{tides}` = number of tidal constituents
* :math:`f_k` = nodal factor
* :math:`\alpha_k` = tidal amplitude
* :math:`\omega_k` = tidal frequency
* :math:`\phi_k` = tidal phase

Flow Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^

At flow-specified boundaries, the normal component of flow is prescribed:

.. math::

    U_n = Q_n(t)

where :math:`U_n` is the normal component of velocity and :math:`Q_n(t)` is the specified discharge per unit width.

Radiation Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Radiation boundaries allow waves to exit the domain with minimal reflection. ADCIRC implements the Sommerfeld radiation condition:

.. math::

    \frac{\partial \zeta}{\partial t} + c \frac{\partial \zeta}{\partial n} = 0

where :math:`c = \sqrt{gH}` is the shallow water wave speed and :math:`n` is the outward normal direction.

Combined Radiation-Flow Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For cases where both wave radiation and external flow need to be accommodated:

.. math::

    \frac{\partial \zeta}{\partial t} + c \frac{\partial \zeta}{\partial n} = \frac{Q_n(t)}{H}

Internal Boundary Conditions
----------------------------

Internal boundaries represent features within the computational domain, such as islands, barriers, and hydraulic structures.

No-Flow Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^

Applied at land boundaries where fluid cannot penetrate:

.. math::

    U_n = 0

Cross-Barrier Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For hydraulic structures like weirs and barriers, the flow across the structure is modeled as:

.. math::

    Q = C_d b H_u^{3/2} \sqrt{2g(H_u - H_d)}

where:

* :math:`Q` = flow rate across the structure
* :math:`C_d` = discharge coefficient
* :math:`b` = structure width
* :math:`H_u` = upstream water depth
* :math:`H_d` = downstream water depth

When the structure is submerged, the equation is modified to account for flow over the top of the structure.

Surface Boundary Conditions
---------------------------

Wind Stress
^^^^^^^^^^^

The wind stress at the water surface is parametrized as:

.. math::

    \frac{\tau_{s\lambda}}{\rho_0} = C_D \frac{\rho_a}{\rho_0} |W| W_\lambda

.. math::

    \frac{\tau_{s\phi}}{\rho_0} = C_D \frac{\rho_a}{\rho_0} |W| W_\phi

where:

* :math:`C_D` = wind drag coefficient
* :math:`\rho_a` = air density
* :math:`|W|` = wind speed magnitude at 10m height
* :math:`W_\lambda, W_\phi` = components of the wind velocity

The wind drag coefficient is typically modeled as a function of wind speed:

.. math::

    C_D = 
    \begin{cases}
    C_{D1} & \text{if } |W| \leq W_1 \\
    C_{D1} + (C_{D2} - C_{D1}) \frac{|W| - W_1}{W_2 - W_1} & \text{if } W_1 < |W| < W_2 \\
    C_{D2} & \text{if } |W| \geq W_2
    \end{cases}

Atmospheric Pressure
^^^^^^^^^^^^^^^^^^^^

The effect of atmospheric pressure is included in the momentum equations as the inverse barometer effect:

.. math::

    \frac{1}{\rho_0} \frac{\partial p_s}{\partial \lambda}, \frac{1}{\rho_0} \frac{\partial p_s}{\partial \phi}

where :math:`p_s` is the atmospheric pressure at the sea surface. 