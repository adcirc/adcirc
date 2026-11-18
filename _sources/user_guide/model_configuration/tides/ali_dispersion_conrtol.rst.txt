.. meta::
   :description: AliDispersionControl in ADCIRC
   :keywords: adcirc, alidispersioncontrol

.. _ali_dispersion_control:

Ali's Dispersion Control
========================

**Ali's Dispersion Correction** is a correction to the dispersive behavior of
very long waves in a compressible ocean on an elastic Earth. It is intended to
be used in lieu of the self-attraction and loading tide prescribed through the
`fort.24 file <fort.24_file>`__. It should be used with the sixth-digit of
`IM <IM>`__ equal to 3 (fully implicit gravity wave term).

Version
-------

.. raw:: mediawiki

   {{Version support box|version=55|relation=+|support=tp}}

This is considered a technical preview in version 55. Theoretical work is still
ongoing.

.. _controlling_dispersive_behavior:

Controlling Dispersive Behavior
-------------------------------

| The feature is triggered by the presence of the &AliDispersionControl namelist
  at the bottom of the `fort.15 file <fort.15_file>`__. Here is an example of
  how this line is used (the following floats are the default values):
| ``&AliDispersionControl CAliDisp=T, Cs=1500.0, Ad = 0.0050189, Bd = 0.23394/``

-  ``CAliDisp`` logical flag to turn Ali's dispersion correction on (F=false by
   default).
-  ``Cs`` is the speed of sound in water [m/s] for the Mach number based
   dispersive correction. Set Cs to a negative value to turn the Mach number
   based correction off.
-  ``Ad`` the power law constant.
-  ``Bd`` the power law exponent.

For the following equation:
:math:`Correction = 1 - \frac{Ma^2}{4} - {A_d}H^{B_d}`

where Ma is the Mach number (:math:`Ma = \frac{\sqrt{gH}}{Cs}`) and H is the
total water depth.
