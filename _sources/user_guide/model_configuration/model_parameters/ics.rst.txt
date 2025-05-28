.. meta::
   :description: ICS in ADCIRC
   :keywords: adcirc, ics

.. _ics_parameter:

ICS
===

**ICS** is a fundamental parameter in the :ref:`fort.15 file <fort15>`
that defines the coordinate system and the desired projection. The value of ICS
also has an important consequence for the choice of the Coriolis
:ref:`CORI <CORI>` parameter of the :ref:`fort.15 file <fort15>`.

.. _available_ics_values:

Available ICS Values
--------------------

.. list-table::
   :widths: 10 15 40
   :header-rows: 1
   :class: wrap-table

   * - ICS Value
     - Short-name
     - Description
   * - 1
     - Cartesian
     - Points in the fort.14 are already mapped onto an arbitrary Cartesian coordinate system, e.g., UTM. Also useful for idealized problems.
   * - 2
     - Geographic, CPP, no curvature
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the CPP (equidistant) cylindrical mapping by ADCIRC. The curvature of the Earth is not accounted for.
   * - 20
     - Geographic, Equal-area
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the Equal-area cylindrical mapping by ADCIRC. The curvature of the Earth is correctly accounted for.
   * - 21
     - Geographic, CPP
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the CPP (equidistant) cylindrical mapping by ADCIRC. The curvature of the Earth is correctly accounted for.
   * - 22
     - Geographic, Mercator
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the Mercator (conformal) cylindrical mapping by ADCIRC. The curvature of the Earth is correctly accounted for.
   * - 23
     - Geographic, Miller
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the Miller cylindrical mapping by ADCIRC. The curvature of the Earth is correctly accounted for.
   * - 24
     - Geographic, Gall-Stereographic
     - Points in the fort.14 are specified in geographic coordinates, which will be projected using the Gall-Stereographic cylindrical mapping by ADCIRC. The curvature of the Earth is correctly accounted for.

Implications
------------

Beginning from Version 55, ICS values equal to 20-24 will be possible. These are
intended to replace the old method of specifying ICS = 2, but this option is
retained for regression testing. When ICS = 2, the curvature of the Earth is not
correctly accounted for in the Spherical coordinate form of the governing
equations which becomes more important as the geographic size of the
computational domain size increases. In recent years global modeling using
ADCIRC has been successful [1]_, e.g.,
`GLOCOFFS <https://wpringle.github.io/GLOCOFFS/>`__ where it was found that the
old method was deficient. The new options using values ICS = 20-24 now account
for the curvature correctly, and should in general be always used on
geographical domains (ICS = 1 should still be used for Cartesian coordinate
domains). ICS = 22 is particularly attractive because it uses a conformal
mapping (Mercator) that conserve the angles on the spherical Earth, but testing
has generally found that all choices of ICS = 20-24 give effectively the same
answers.

.. _negative_ics_value_for_rotation:

Negative ICS Value for Rotation
-------------------------------

Also beginning from Version 55, the values 20-24 can also be set to a negative
value (i.e, -20, -21, ...) to implement an arbitrary rotation of the geographic
coordinates by ADCIRC. This is primarily used to ensure that Earth's poles are
rotated onto land to eliminate the singularity in the Spherical coordinate form
of the governing equations. ADCIRC's outputs will be displayed on the original
unrotated coordinate system.
The rotation is set by the :ref:`fort.rotm <fortrotm>` input file (see the link
for example formats).

References
----------

.. raw:: html

   <references />

.. [1]
   Pringle et al., Global Ocean-to-Coastal Storm Tide Modeling in ADCIRC v55:
   Unstructured Mesh Design, in preparation (2020)
