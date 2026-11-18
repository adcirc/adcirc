.. meta::
   :description: Horizontal eddy viscosity in ADCIRC
   :keywords: adcirc, horizontal eddy viscosity

.. _horizontal_eddy_viscosity:

Horizontal Eddy Viscosity
=========================

**Horizontal eddy viscosity** is a term in the momentum equations to cover the
well-known turbulence closure problem. The term is also useful numerically for
increased stability. For background on the turbulence closure problem, see, e.g.
the `AMS glossary
definition <http://glossary.ametsoc.org/wiki/Turbulence_closure>`__, the
`Wikipedia entry on turbulence
modeling <https://en.wikipedia.org/wiki/Turbulence_modeling>`__, or a fluid
mechanics textbook such as Kundu's [1]_. Readers unfamiliar with these concepts
should note that, contrary to the name, eddy viscosity is unrelated to true
viscosity, even though the units (length squared over time) are the same.
Typically, one should expect eddy viscosity to be many orders of magnitude
larger than the viscosity of actual water (around 1x10\ :sup:`-6`
meters\ :sup:`2`/second for seawater). In ADCIRC, the user can either specify a
constant in time horizontal eddy viscosity or employ the Smagorinsky model to
update the horizontal eddy viscosity based on local flow conditions with each
time step. In both cases the horizontal eddy viscosity or Smagorinsky
coefficient can be either spatially constant or spatially varying.

.. _user_specified_constant_in_time_eddy_viscosity:

User-Specified Constant in Time Eddy Viscosity
----------------------------------------------

A spatially constant eddy viscosity can be supplied via the ``ESLM`` parameter
in the `fort.15 file <fort.15_file>`__. Alternatively, the
:ref:`average_horizontal_eddy_viscosity_in_sea_water_wrt_depth
<average_horizontal_eddy_viscosity_in_sea_water_wrt_depth>`
nodal attribute permits a spatially variable eddy viscosity.

Usage Notes
~~~~~~~~~~~

Commonly used values range from 1 to 50 meters\ :sup:`2`/second, with 10 being a
good starting point for many coastal ocean modeling scenarios. By the nature of
turbulence closure, one might expect the value to be smaller as mesh resolution
increases. However, this may not always be the case since areas with high
resolution may also be more turbulent. In particular, some modelers use a larger
value in overland areas, which can also improve stability. However, specifying a
value that is too large can induce instabilities, as well.

.. _smagorinsky_eddy_viscosity:

Smagorinsky Eddy Viscosity
--------------------------

| ADCIRC also allows a
  `Smagorinsky-type <https://en.wikipedia.org/wiki/Turbulence_modeling#Smagorinsky_model_for_the_sub-grid_scale_eddy_viscosity>`__
  turbulence closure model [2]_. There are currently two slightly different ways
  of enabling this in the model. One is by specifying a negative value for
  :ref:`ESLM` <ESLM>` in the :ref:`fort.15 file <fort15>`. If this is done,
  the Smagorinsky turbulence closure is enabled and the absolute value of
  ``ESLM`` is used as the coefficient. Alternatively, the ``Smag_Control``
  namelist can be used, e.g. by adding the following to the end of the fort.15
  file
| ``&Smag_Control SMAG_LOWER_LIM=1.0e-8, SMAG_UPPER_LIM=100 /``
| Here again, the absolute value of ``ESLM`` is taken as the coefficient, though
  a negative ``ESLM`` value isn't necessary if the namelist is used. The
  namelist also allows the user to optionally specify lower and upper bounds on
  the actual eddy viscosity determined via the Smagorinsky formulation. The
  values given in the above example of 1E-8 m\ :sup:`2`/s and 100 m\ :sup:`2`/s
  are the default values. There is no way to supply upper and lower bounds with
  the other method.

If the
:ref:`average_horizontal_eddy_viscosity_in_sea_water_wrt_depth <average_horizontal_eddy_viscosity_in_sea_water_wrt_depth>`
nodal attribute is used when the Smagorinsky eddy viscosity is employed, then
the nodal attribute values are taken to be the Smagorinsky coefficient values.

Usage Notes
~~~~~~~~~~~

The two methods currently result in different flags being set in the code, which
have slightly inconsistent behavior. Specifically, the "negative ``ESLM``
method" sets ``CSmag_Eh=.TRUE.`` and results in a check being run at model
initialization to verify that the default (Kolar-Gray) formulation of the
lateral stresses in the GWCE is not in use, and the model errors out if it is.
That check is not done in the namelist method, which sets
``Smag_Comp_Flag=.TRUE.``. The lateral stress formulation is specified by the
first digit in a :ref:`six-digit IM code <im_code>` being set to 1.
It is also set for any regular IM value for 2D ADCIRC. Mathematically, the
Kolar-Gray formulation of the lateral stress does not take into account a
derivative of the eddy viscosity terms, meaning it is technically inconsistent
with a Smagorinsky turbulence closure model. However, it is not clear whether
the time derivative term is ever large enough to be consequential, and therefore
merit departure from ADCIRC's default GWCE formulation.

References
----------

.. raw:: html

   <references />

.. [1]
   Pijush K. Kundu, Ira M. Cohen, David R. Dowling (2012). Fluid Mechanics.

.. [2]
   Smagorinsky, J. “General Circulation Experiments with the Primitive
   Equations.” Monthly Weather Review 91, no. 3 (March 1, 1963): 99–164.
   https://doi.org/10.1175/1520-0493(1963)091\ \ <0099:GCEWTP>2.3.CO;2.
