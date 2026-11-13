.. meta::
   :description: NTIP in ADCIRC
   :keywords: adcirc, ntip

.. _ntip_parameter:

NTIP
====

**NTIP** is an input in the :ref:`fort.15 file <fort15>` that selects the
astronomical forcing input type.

Parameter Summary
-----------------

The following table is a summary of possible NTIP values, description, details,
and other necessary input parameters and files.

.. list-table::
   :header-rows: 1
   :widths: 10 20 40 30
   :class: wrap-table

   * - NTIP Value
     - Description
     - Details
     - Other Required Inputs
   * - 0
     - No Astronomical Forcing
     - -
     - :ref:`NTIF <ntif_parameter>` = 0
   * - 1
     - Astronomical Tidal Potential
     - Reconstructs the tidal elevation using the analytical formulation for the equilibrium tidal potential [1]_
     - :ref:`NTIF <ntif_parameter>` > 0
   * - 2
     - Astronomical Tidal Potential plus Self-attraction and Loading (SAL) Tide [2]_
     - Reconstructs the tidal elevation by summing the contribution from the analytical formulation for the equilibrium tidal potential [1]_ with the contribution from the prescribed SAL constituent values found in the `fort.24 file <fort.24_file>`__.
     - :ref:`NTIF <ntif_parameter>` > 0, :ref:`fort.24 file <fort24>`

References
----------

.. raw:: html

   <references />

.. [1]
   Eq. (27), page 17 in Luettich, R.A., Westerink, J.J., 1992. ADCIRC: an
   advanced three-dimensional circulation model for shelves coasts and
   estuaries, report 1: theory and methodology of ADCIRC-2DDI and ADCIRC-3DL,
   Dredging Research Program. Vicksburg, MS.

.. [2]
   Ray, R.D., 1998. Ocean self-attraction and loading in numerical tidal models.
   Mar. Geod. 21, 181–192. doi:10.1080/01490419809388134
