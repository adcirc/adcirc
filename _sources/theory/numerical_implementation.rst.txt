Numerical Implementation
========================

ADCIRC employs a finite element method for spatial discretization and a finite difference scheme for time discretization. This section details the numerical methods used to solve the governing equations.

Finite Element Spatial Discretization
-------------------------------------

ADCIRC uses the continuous Galerkin finite element method with linear triangular elements for spatial discretization. This approach allows for:

* Flexible mesh design with high resolution in areas of interest
* Accurate representation of complex coastlines and bathymetry
* Efficient handling of varying spatial scales

The approximation for a variable :math:`\psi` (representing either :math:`\zeta`, :math:`U`, or :math:`V`) within an element is given by:

.. math::

    \psi(x,y,t) \approx \sum_{j=1}^{N_n} \psi_j(t) \phi_j(x,y)

where:

* :math:`N_n` is the number of nodes in the mesh
* :math:`\psi_j(t)` is the value of the variable at node :math:`j` at time :math:`t`
* :math:`\phi_j(x,y)` is the basis function associated with node :math:`j`

For linear triangular elements, the basis functions are:

.. math::

    \phi_j(x,y) = 
    \begin{cases}
    a_j + b_j x + c_j y & \text{within elements containing node } j \\
    0 & \text{elsewhere}
    \end{cases}

where coefficients :math:`a_j`, :math:`b_j`, and :math:`c_j` are determined by the geometry of the element.

Weak Formulation and Assembly
-----------------------------

The governing equations are converted to their weak form by:

1. Multiplying by a test function (chosen from the same space as the basis functions)
2. Integrating over the domain
3. Applying integration by parts to reduce the order of derivatives

This process yields a system of equations that can be assembled into a global system:

.. math::

    M \frac{d^2 \zeta}{dt^2} + \tau_0 M \frac{d\zeta}{dt} + K \zeta = F

for the GWCE, and:

.. math::

    M \frac{dU}{dt} = F_U

.. math::

    M \frac{dV}{dt} = F_V

for the momentum equations. Here:

* :math:`M` is the global mass matrix
* :math:`K` is the global stiffness matrix
* :math:`F`, :math:`F_U`, and :math:`F_V` contain the remaining terms

Time Discretization
-------------------

ADCIRC uses multiple time-stepping schemes:

1. **GWCE (Continuity) Equation**:
   
   A three-level explicit scheme is employed:

   .. math::

       M \frac{\zeta^{n+1} - 2\zeta^n + \zeta^{n-1}}{\Delta t^2} + \tau_0 M \frac{\zeta^{n+1} - \zeta^{n-1}}{2\Delta t} + K \zeta^n = F^n

2. **Momentum Equations**:
   
   A Crank-Nicolson implicit scheme is used for the linear terms and an explicit scheme for nonlinear terms:

   .. math::

       M \frac{U^{n+1} - U^n}{\Delta t} = \frac{1}{2} [L(U^{n+1}) + L(U^n)] + N(U^n, V^n)

   where :math:`L` represents linear terms and :math:`N` represents nonlinear terms.

Lumped Mass Matrix
------------------

To enhance computational efficiency, ADCIRC uses a lumped (diagonal) mass matrix instead of the consistent mass matrix. This is achieved by row-summing:

.. math::

    M_{ii}^{lumped} = \sum_{j=1}^{N_n} M_{ij}^{consistent}

This approximation significantly reduces computational cost while maintaining acceptable accuracy.

Stabilization Techniques
------------------------

Several numerical stabilization techniques are employed:

1. **Selective Lumping**: Different lumping ratios for different terms
2. **Spatial Filtering**: Applied periodically to remove high-frequency noise
3. **Flux-Corrected Transport**: Used for advection terms to prevent unphysical oscillations
4. **Artificial Viscosity**: Added adaptively to maintain stability in steep gradient regions

These techniques ensure the model remains stable under challenging conditions (e.g., wetting/drying, storm surge, steep bathymetry). 