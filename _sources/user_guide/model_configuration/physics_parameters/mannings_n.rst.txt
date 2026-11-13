.. meta::
   :description: Manning's n at sea floor in ADCIRC
   :keywords: adcirc, manning's n at sea floor

.. _mannings_n:

Manning's n
===========

The **``mannings_n_at_sea_floor``** nodal attribute is one of the available
options for specifying bottom friction in ADCIRC based on
the `Manning formula <https://en.wikipedia.org/wiki/Manning_formula>`_. It is a
nodal attribute in the :ref:`fort.13 file <mannings_n_at_sea_floor>`,
meaning that it allows spatially (and temporally) varying bottom friction in the
model.

.. _attribute_summary:

Attribute Summary
-----------------

If the user elects to use this nodal attribute, which ADCIRC reads in as the
``ManningsN`` variable, :ref:`NOLIBF` must be set to 1 or the run will terminate.
During execution, the Manning’s n values specified are converted to equivalent
quadratic friction coefficients before the bottom stress is calculated. The
equivalent quadratic friction coefficient is calculated according to the
following formula at each node at each time step:

.. math:: C_d(t)=\frac{gn^2}{\sqrt[3]{h+\eta(t)}}

where,

-  *C*\ :sub:`d` drag coefficient
-  *t* time
-  *g* acceleration due to gravity
-  *n* Manning's n
-  *h* depth
-  *η* water surface elevation

The addition of the water surface elevation is conditional upon the setting of
:ref:`NOLIFA <NOLIFA>`: *η* is treated as zero if ``NOLIFA = 0`` in the
:ref:`fort.15 file <fort15>`. Finally, the value of :ref:`CF <CF>` in the
:ref:`fort.15 file <fort15>` is used to set a lower limit on the resulting
equivalent quadratic friction coefficient, under the assumption the *C*\ :sub:`d` calculated from
this formula tends to become too small in deep water.

.. figure:: /_static/images/user_guide/model_configuration/physics_parameters/manning_s_n/mnasfContours3.png
   :width: 400px
   :align: center
   :alt: Drag coefficient :math:`C_d` contours for Manning's n and total depth

   Drag coefficient :math:`C_d` contours for Manning's n and total depth

Negative *n* Values
-------------------

.. raw:: mediawiki

   {{ADC version|version=55|relation=ge}}

If a combination of Manning's and time invariant bottom
friction is desired, users can elect to set the Manning's
nodal attribute to a negative value at certain mesh nodes.
At mesh nodes where *n* is negative, the time invariant *C*\ :sub:`d`
coefficient, specified via a constant :ref:`CF <CF>` in the
:ref:`fort.15 file <fort15>` or the spatially varying
:ref:`quadratic_friction_coefficient_at_sea_floor <quadratic_friction_coefficient_at_sea_floor>`
fort.13 nodal attribute, will be used instead.

Specifying *n* Values
---------------------

Manning's *n* is often assigned using land cover datasets, when available.
Examples of commonly used land cover data in the US are the National Land Cover
Dataset and the Coastal Change Analysis Program. There is a broad literature on
specification of *n* values based on laboratory and fields studies, the most
classical example of which is Chow (1959). Ideally, field surveys or review of
on-site photography should be done to correlate land cover values to *n*.
Alternatively, review of available literature may provide some basis for
selecting values.

Utilities
---------

-  `f13builder <https://adcirc.org/home/related-software/adcirc-utility-programs/>`_
