.. meta::
   :description: NOLIFA in ADCIRC
   :keywords: adcirc, nolifa

.. _nolifa_parameter:

NOLIFA
======

**NOLIFA** is a parameter in the :ref:`fort.15 file <fort15>` controlling
controlling the finite amplitude terms and wetting-drying in ADCIRC. The value
of NOLIFA effects the meaning of the minimum water depth parameter (:ref:`H0`)
and requires the specification of additional parameters together with
:ref:`H0`.

Parameter Summary
-----------------

The following table is a summary of possible NOLIFA values, description,
details, and other necessary or recommended input parameters.

.. list-table::
   :header-rows: 1
   :widths: 10 20 40 30
   :class: wrap-table

   * - NOLIFA Value
     - Description
     - Details
     - Required/Recommended Inputs
   * - 0
     - No finite amplitude terms and no wetting-drying
     - The depth is linearized by using the bathymetric depth, rather than the total depth, in all terms except the transient term in the continuity equation. Wetting and drying of elements is disabled. Initial water depths are assumed equal to the bathymetric water depth specified in the `fort.14 file <fort.14_file>`__.
     - :ref:`H0 <H0>`, :ref:`NOLICA <NOLICA>` = 0
   * - 1
     - Finite amplitude terms without wetting-drying
     - Finite amplitude terms are included in the model run and wetting and drying of elements is disabled. Initial water depths are assumed equal to the bathymetric water depth specified in the :ref:`fort.14 file <fort14>`.
     - :ref:`H0 <H0>`, :ref:`NOLICA <NOLICA>` = 1
   * - 2
     - Finite amplitude terms with wetting-drying
     - Finite amplitude terms are included in the model run and wetting and drying of elements is enabled. Initial water depths are assumed equal to the bathymetric water depth specified in the :ref:`fort.14 file <fort14>`.
     - :ref:`H0 <H0>`, :ref:`NOLICA <NOLICA>` = 1


Usage Notes
-----------

When the finite amplitude terms are turned on, the time derivative portion of
the advective terms should also be turned on for proper mass conservation and
consistency (i.e., when NOLIFA > 0, then `NOLICA <NOLICA>`__ = 1).


