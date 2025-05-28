.. meta::
   :description: NOLIBF in ADCIRC
   :keywords: adcirc, nolibf

.. _nolibf_parameter:

NOLIBF
======

**NOLIBF** is a parameter in the :ref:`fort.15 file <fort15>` controlling the
type of bottom stress parameterization used in a 2DDI ADCIRC run. This parameter
must be specified but is ignored in a 3D run.

Parameter Summary
-----------------

The following table is a summary of possible NOLIBF values, description,
details, and other necessary input parameters and files.

.. list-table::
   :header-rows: 1
   :widths: 10 20 50 20
   :class: wrap-table

   * - NOLIBF Value
     - Description
     - Details
     - Other Required Inputs
   * - 0
     - linear bottom friction law
     - uses a spatially constant linear bottom friction, :math:`\tau_b` of dimensions [1/T]
     - `FFACTOR <FFACTOR>`__
   * - 1
     - quadratic bottom friction law
     - can specify a spatially constant or spatially varying dimensionless coefficient, :math:`C_f`; in which :math:`\tau_b = C_f|U|/H`
     - `FFACTOR <FFACTOR>`__
   * - 2
     - hybrid nonlinear bottom friction law
     - In deep water, the friction coefficient is constant and a quadratic bottom friction law results. In shallow water the friction coefficient increases as the depth decreases (e.g. as in a Manning-type friction law). The friction coefficient is determined as: FFACTOR =FFACTORMIN*[1+(HBREAK/H)**FTHETA]**(FGAMMA/FTHETA)
     - `FFACTORMIN <FFACTORMIN>`__, `HBREAK <HBREAK>`__, `FTHETA <FTHETA>`__, `FGAMMA <FGAMMA>`__

.. _usage_notes:

Usage Notes
-----------

In the :ref:`NWP <nwp>` section, if the user selects
quadratic_friction_coefficient_at_sea_floor, mannings_n_at_sea_floor, or
chezy_friction_coefficient_at_sea_floor, then NOLIBF must be 1 (nonlinear
friction formulation) since all those formulations are nonlinear. If the NOLIBF
were anything other than 1, it is an error that will cause ADCIRC to stop.

For 3D ADCIRC runs, spatially varying bottom friction should be specified using
the :ref:`bottom_roughness_length <bottom_roughness_length>` nodal attribute.

