Bridge Piers
============

Given the typical resolution of coastal circulation models, (e.g., in nearshore applications resolved scales are usually 10s to 100s of meters and larger), it is rarely practical to solve explicitly for the small scale flow around obstructions such as bridge pilings (diameters of meters). However, in some situations it may be desirable to include the effects of these subgrid scale obstructions on the resolved scale flow. To accomplish this, a subgrid scale obstruction parameterization has been developed and implemented in ADCIRC. Since the cross-sectional area of these obstructions is usually quite small compared to the cross-sectional area of the flow field (and thus the effect on flow continuity is quite small), their primary effect is to impart additional drag on the resolved scale flow. This section describes the theory behind the obstruction parameterization implemented in ADCIRC.

The bridge pier option is turned on by setting the parameter ``NWP`` to add for an additional Nodal Attribute in the Model Parameter and Periodic Boundary Condition File. Required input coefficients are read in at specified nodes from the Nodal Attributes file (fort.13). Additional information can be found in Luettich and Westerink (1999).

Theory
------

The two-dimensional, vertically integrated momentum equations used in ADCIRC contain bottom friction (drag) terms. To represent the extra drag caused by subgrid scale obstructions such as bridge piers, a second contribution has been added to the bottom friction terms:

* Bottom friction x + obstruction drag x 
* Bottom friction y + obstruction drag y

Options in ADCIRC exist to express the bottom friction terms as linear, quadratic or hybrid quadratic/Manning's n functions of flow velocity. The obstruction drag is assumed to be due primarily to form drag and therefore is represented as quadratic in velocity:

* obstruction drag x = :math:`C_D \frac{U\sqrt{U^2+V^2}}{H}`
* obstruction drag y = :math:`C_D \frac{V\sqrt{U^2+V^2}}{H}`

where :math:`C_D` is an obstruction (piling) drag coefficient.

Determination of Drag Coefficient
---------------------------------

To determine :math:`C_D` for a series of bridge pilings, we have utilized the extensive body of research conducted in the early to mid 1900s on flow around pilings. Yarnell (1934a,b), as summarized by Henderson (1966), fit the change in water level :math:`\Delta H` for steady, unidirectional flow past a piling to the following relation:

:math:`\frac{\Delta H}{H} = K \left( \alpha + 5\alpha^2 - 0.6 \right) F_r^2`

where:

* H is the total water depth
* K is a pier shape factor (see table below)
* :math:`F_r^2` is the square of the Froude number
* :math:`\alpha` is the fraction of the cross section obstructed by the piling
* :math:`U_s` is the velocity in the along stream (s) direction
* g is the acceleration of gravity

This equation is considered to be valid so long as :math:`\alpha < 0.5`.

.. list-table:: Recommended Pier Shape Factors
   :header-rows: 1
   :widths: 70 30

   * - Pier Shape
     - K
   * - Semicircular nose and tail
     - 0.9
   * - Lens-shaped nose and tail
     - 0.9
   * - Twin-cylinder piers with connecting diaphragm
     - 0.95
   * - Twin-cylinder piers without diaphragm
     - 1.05
   * - 90 deg triangular nose and tail
     - 1.05
   * - Square nose and tail
     - 1.25

Implementation in ADCIRC
------------------------

For steady, unidirectional flow in the vicinity of pilings, the corresponding ADCIRC momentum equations simplify to:

:math:`g \frac{\partial \zeta}{\partial s} + \tau_{bf} + \tau_{bp} = 0`

where :math:`\frac{\partial \zeta}{\partial s}` is the free surface gradient in the along stream direction. Approximating this gradient as :math:`\frac{\Delta H}{L}`, Yarnell's equation and the simplified momentum equation can be combined to obtain an expression for :math:`C_D`:

:math:`C_D = \frac{gH}{U_s^2} K \left( \alpha + 5\alpha^2 - 0.6 \right) F_r^2`

This equation is used at each time step in ADCIRC, to compute :math:`C_D` at all nodes associated with bridge pilings. This allows the obstruction drag to be computed and added to the bottom friction terms as discussed above.

Note that Yarnell's equation is based on H and :math:`U_s` values measured at a point downstream of the pilings. ADCIRC uses H and :math:`U_s` at the effective piling location to minimize complications due to changing flow direction, e.g., due to reversal of the tide. Typical elevation and velocity changes that result are small enough that the overall results are not compromised by this approximation.

References
----------

Henderson, F.M., 1966, Open Channel Flow, MacMillan Publishing Co., New York, pp 265-267.

Luettich, R.A., Jr. and J.J. Westerink, 1999, Implementation of Bridge Pilings in the ADCIRC Hydrodynamic Model: Upgrade and Documentation for ADCIRC Version 34.19, Contractors Report, department of the Army, US Army Corps of Engineers, Waterways Experiment Station, Vicksburg, MS, November 19, 1999, 8 p.

Yarnell, D.L., 1934a, Pile trestles as channel obstructions, U.S. Dept of Agriculture, Tech. Bull. 429.

Yarnell, D.L., 1934b, Bridge piers as channel obstructions, U.S. Dept of Agriculture, Tech. Bull. 442. 