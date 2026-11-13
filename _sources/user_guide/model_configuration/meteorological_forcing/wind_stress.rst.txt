.. meta::
   :description: Wind Stress in ADCIRC
   :keywords: adcirc, wind stress

.. _wind_stress:

Wind Stress
===========

When wind blows over the water, it exerts a shear stress at the water surface
that transfers horizontal momentum vertically downward across the air–sea
interface, driving circulation. In ADCIRC, wind stress is an input forcing term,
with several different formats provided. See the :ref:`NWS Parameter <nws_parameter>` 
for available formats. In most cases, the exact wind stress to be applied to the
model is not provided, therefore ADCIRC must determine how to convert a given
wind speed to the actual stress applied at the ocean surface. This page covers
the various aspects of this process, as well as the options available to the
user.

.. _definition_of_winds:

Definition of Winds
-------------------

The characteristics of wind forcing are often broken down in three ways:

#. Whether the winds are considered to be over-water (termed "marine exposure")
   or over-land
#. The elevation above the sea (or ground) surface of the winds
#. The time-averaging (if any) that has been applied

ADCIRC generally expects 10-meter, 10-minute winds at their actual exposure,
although the exact expectations vary depending on the input type. For instance,
when Holland-type wind inputs are provided (e.g. ``NWS=8`` or ``NWS=20``), the
wind speed is expected to be the **1-minute**\  `maximum sustained
wind <https://en.wikipedia.org/wiki/Maximum_sustained_wind>`_ at 10 meters
elevation. If marine-exposure winds are provided, then :ref:`surface
roughness <surface_directional_effective_roughness_length>` reductions may be needed  [1]_
[2]_ . If winds are provided with a different averaging time, then an
appropriate correction may be needed, though winds with averaging times of 10 to
60 minutes are generally considered to be quite similar; this is the so-called
mesoscale gap. For recommendations on wind time-scale conversions not handled
internally by ADCIRC, for tropical cyclones, see the WMO guidelines of Harper et
al.  [3]_

.. _roughness_reductions:

Roughness Reductions
--------------------

Reductions in wind speed to convert to the appropriate exposure come from a
logarithmic boundary layer formulation (see, e.g.  [1]_  [2]_) to determine a
fraction :math:`f` to reduce the winds,

.. math:: f=\left ( \frac{z_{0l}}{z_{0m}}\right ) ^{0.0706} \left ( \frac{\ln \frac{10}{z_{0l}}}{\ln \frac{10}{z_{0m}}} \right )

for marine roughness length :math:`z_{0m}` and reduced ("land") roughness length
:math:`z_{0l}`. Wind speed is then reduced as,

.. math:: \mathbf{w}'=f\mathbf{w}=f [u,v]

for x- and y- wind vector components :math:`u` and :math:`v`. The marine
roughness length is,

.. math:: z_{0m}=\frac{0.018}{g} c_d \left \Vert \mathbf{w} \right \|

for Charnock parameter :math:`0.018`, drag coefficient :math:`c_d` and
acceleration due to gravity :math:`g`. The wind drag coefficient is addressed
below in `this section <#Converting_Wind_Velocity_to_Wind_Stress>`__. As
previously noted, the reduced ("land") roughness length :math:`z_{0l}` is
specified by the user via the `surface
roughness <Fort.13_file#Surface_Roughness>`__ nodal attribute. The fraction
:math:`f` is bounded on :math:`[0,1]`, meaning the winds cannot be increased,
nor change direction.

.. _older_behavior:

Older Behavior
~~~~~~~~~~~~~~

.. _interpolating_roughness_lengths_before_v55:

Interpolating Roughness Lengths Before v55
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/wind_stress/WindTraceV2.png
   :alt: Comparison of output wind velocities using old and new methods.
   :width: 500px
   :align: center


.. raw:: mediawiki

   {{ADC version|version=55|relation=lt}}

| Before version 55, the directional wind reductions were applied by determining
  which of the 12 directional bins the wind velocity (at each time step) fell
  into, and using that roughness reduction, i.e. nearest neighbor interpolation.
  Starting in version 55, the roughness length is linearly (in angle space)
  interpolated between directional bins. In testing, this has been found to
  generally have a very small effect on water levels, but a notable effect on
  wind speeds, since time evolution of winds is smoother. It can have large
  localized effects on water levels in cases where there are large changes in
  neighboring roughness length bins coinciding with well-aligned winds, as in
  `this test case <:File:MaxeleDiffRun13MinusRun12View1.png>`__ with Hurricane
  Isaac. 
  
  |Video of wind speeds before (left) and after (right) adding linear
  interpolation of surface roughness.|

.. _roughness_reduction_bug_before_v54:

Roughness Reduction Bug Before v54
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: mediawiki

   {{ADC version|version=54|relation=lt}}

Before version 54, there was a bug in this calculation. The mistake and its
effects are addressed in this PDF document: \ **ADD A LINK TO A PDF HERE**\ .

.. _converting_wind_velocity_to_wind_stress:
.. _metControl_namelist:

Converting Wind Velocity to Wind Stress
---------------------------------------

In ADCIRC, four formulations are available to convert wind velocities to the
wind stresses applied in the momentum equations. Although there are several ways
to control this, users are generally encouraged to use the
:ref:`metControl <metControl>` namelist in the :ref:`fort.15 <fort15>` file. The
default drag formulation is the Garratt  [4]_ linear formula. An alternative for
use with tropical cyclones is the Powell formulation,  [5]_ which varies drag by
the sector of the tropical cyclone. When ice coverage is included in the model,
a wind drag formulation that accounts for this effect should be used. By
default, if ice coverage input data are supplied, ADCIRC uses a cubic function
of ice coverage, termed the "IceCube" drag formulation. Lastly, the "swell" drag
law option allows users to utilize SWAN's drag formulation when employing the
coupled model.

In all cases, the actual wind drag coefficients determined by ADCIRC can be
output to a :ref:`fort.63 <fort63>`-type file named
``winddrag.173``. Output settings (file format, output start/end
times, and output interval) match those of either the fort.63 or fort.73/74
files. Outputting of this file is enabled by setting
``outputWindDrag=.TRUE.`` in the ``&metControl`` namelist of the :ref:`fort.15 file's
namelist section <fort15>`.

.. _garratt_drag_formulation:

Garratt Drag Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~

This is the default wind drag formulation in ADCIRC. From Garratt (1977) [4]_,
the formula is,

.. math:: c_d=0.001 \left ( 0.75+0.067 \left \Vert \mathbf{w} \right \| \right )

By default, ADCIRC puts an upper bound on the drag coefficient of
:math:`c_d\le0.0035`. This upper bound ``WindDragLimit`` can be changed via the
`metControl <metControl>`__ fort.15 namelist.

.. _powell_drag_formulation:

Powell Drag Formulation
~~~~~~~~~~~~~~~~~~~~~~~

CONTRIBUTOR NEEDED Note that the code for the Powell drag law has not been
configured to reverse orientation in the southern hemisphere, and so will
produce the wrong results. For details on Powell, see `this
post <https://ccht.ccee.ncsu.edu/wind-drag-based-on-storm-sectors/>`__.

.. _icecube_drag_formulation:

IceCube Drag Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~

CONTRIBUTOR NEEDED

References
----------

.. raw:: html

   <references />

.. [1]
   Simiu, E., Scanlan, R.H., 1996. Wind effects on structures: fundamentals and
   applications to design, 3rd ed. ed. John Wiley, New York.

.. [2]
   Simiu, E., Yeo, D., 2018. Wind effects on structures: modern structural
   design for wind, Fourth edition. ed. John Wiley & Sons, Hoboken, NJ.

.. [3]
   Harper, B., Kepert, J., Ginger, J., 2010. Guidelines for converting between
   various wind averaging periods in tropical cyclone conditions (No. WMO/TD-No.
   1555). WMO, Geneva, Switzerland.

.. [4]
   Garratt, J.R., 1977. Review of Drag Coefficients over Oceans and Continents.
   Mon. Wea. Rev. 105, 915–929.
   https://doi.org/10.1175/1520-0493(1977)105\ \ <0915:RODCOO>2.0.CO;2

.. [5]
   Powell, M.D., Vickery, P.J., Reinhold, T.A., 2003. Reduced drag coefficient
   for high wind speeds in tropical cyclones. Nature 422, 279–283.
   https://doi.org/10.1038/nature01481

.. |Video of wind speeds before (left) and after (right) adding linear interpolation of surface roughness.| image:: /_static/images/user_guide/model_configuration/meteorological_forcing/wind_stress/Run03masterVsRun05interp_trimmedAndAnnotatedTry7-min.gif
   :width: 500px
