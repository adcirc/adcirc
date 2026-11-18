.. meta::
   :description: ADCIRC2D+ in ADCIRC
   :keywords: adcirc, adcirc2d+

ADCIRC2D+
=========

\ **⚠ This page is still under development. Content may change or be
incomplete.**\ 

**ADCIRC2D+**, is a mode which allows for the inclusion of baroclinic
(density-driven) effects into the traditionally barotropic (depth-averaged)
system. This is achieved by coupling ADCIRC with temperature and salinity fields
from a coarser-resolution Ocean General Circulation Model (OGCM), enabling
improved representation of long-term sea level dynamics. Description of results
using this capability are available in (Pringle et al. 2019) [1]_ and (Blakely
et al. 2024) [2]_

Version
-------

Some of the capability of **ADCIRC2D+** is available in v55; however
improvements and modifications were made that are present only in v56 and above.

.. _governing_equations:

Governing Equations
-------------------

ADCIRC2D+ solves a modified form of the shallow water equations using the
Generalized Wave Continuity Equation (GWCE). The momentum equation includes
baroclinic pressure gradient and :ref:`internal wave drag
terms <internal_tide_energy_conversion>`:

.. math::

   \frac{\partial \mathbf{U}}{\partial t} + \mathbf{U} \cdot \nabla \mathbf{U} + f \mathbf{k} \times \mathbf{U} = - \nabla\left( \frac{p_s}{\rho_0} + g(\eta - \eta_{EQ} - \eta_{SAL}) \right) - \frac{\text{BPG}}{H} + \frac{\boldsymbol{\tau}_s - \boldsymbol{\tau}_b}{\rho_0 H} - \gamma_D \mathbf{C} \mathbf{U}

Where:

-  **BPG** is the depth-integrated baroclinic pressure gradient
-  **γ\ D** is a scaling parameter for internal wave drag
-  **C** is the internal wave drag tensor

The baroclinic pressure gradient is computed as:

.. math::
   \text{BPG} = \frac{g}{\rho_0} \left[ \int_{-h}^{0} \nabla \left( \int_{0}^{z} (\rho - \rho_0) \, dz' \right) dz + \eta \nabla[\eta(\rho_s - \rho_0)] \right]

The internal tide drag tensor is scaled using a dissipation ratio
:math:`\gamma_D` defined by:

.. math::
   
   \gamma_D = \frac{\text{Diss}_{\text{tidal}}}{\text{Diss}_{\text{total}}} = \frac{\mathbf{U}_{\text{tidal}} \cdot \mathbf{C} \cdot \mathbf{U}_{\text{tidal}}}{\mathbf{U} \cdot \mathbf{C} \cdot \mathbf{U}}

Tidal velocity is estimated using a lagged 25-hour filter and removing the
resultant mean signal from the total velocity. This ensures that dissipation
from internal wave drag occurs predominantly at tidal frequencies, preserving
tidal fidelity in the coupled model.

.. _densitycontrol_namelist:

``densityControl`` Namelist
----------------------------

Activating **ADCIRC2D+** is accomplished through the use of the
``densityControl`` namelist.This namelist has the following options and defaults
(denoted by ``(D)``)

-  ``densityRunType`` (string)

      ``'none' (D)``: By default **ADCIRC2D+** is not active. If
      ``densityRunType='none'`` the rest of the namelist values do not matter.
      ``'prognostic'``
      ``'diagnostic'``

-  ``densityFileName`` (string): The name of the file from which pre-computed
   baroclinic pressure gradients and stratification information are read.
-  ``densityTimeIterator`` (integer): The stride used when reading in
   ``densityFileName``. For example, if ``densityFileName`` contains data at
   hourly timesteps and ``densityTimeIterator=2`` then data will be read in
   every two hours.
-  ``densityForcingType`` (string)

      ``'SigmaT'``
      ``'Salinity'``
      ``'Temperature'``
      ``'SalinityTemperature'``
      ``'Baroclinicgradients'``
      ``'BaroclinicgradientsDispersion'``
      ``'Buoyancyfrequencies'``
      ``'BCForcingOnADCIRCGrid'``
      ``'BuoyancyFrequenciesOnGrid'``

Within the ADCIRC code, the various settings of these namelist parameters is
internally translated to an IDEN value which determines the baroclinic terms in
the shallow water equations to include.

References
----------

.. raw:: html

   </references>

.. [1]
   Pringle, W. J., Gonzalez-Lopez, J., Joyce, B. R., Westerink, J. J., van der
   Westhuysen, A. J. (2019). *Baroclinic Coupling Improves Depth-Integrated
   Modeling of Coastal Sea Level Variations Around Puerto Rico and the U.S.
   Virgin Islands*. *Journal of Geophysical Research: Oceans*, 124, 2196-2217.
   https://doi.org/10.1029/2018JC014682

.. [2]
   Blakely, C. P., Wirasaet, D., Cerrone, A. R., Pringle, W. J., Zaron, E. D.,
   Brus, S. R., Seroka, G. N., Moghimi, S., Myers, E. P., & Westerink, J. J.
   (2024). *Dissipation Scaled Internal Wave Drag in a Global Heterogeneously
   Coupled Internal/External Mode Total Water Level Model*. *Journal of Advances
   in Modeling Earth Systems*, 16, e2024MS004502.
   https://doi.org/10.1029/2024MS004502
