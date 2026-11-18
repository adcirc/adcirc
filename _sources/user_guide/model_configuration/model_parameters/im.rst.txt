.. meta::
   :description: IM in ADCIRC
   :keywords: adcirc, im

.. _im_code:

IM
==

**IM** is an important parameter in the `fort.15 file <fort.15_file>`__ that
defines numerical model formulation and dimension. Among other things, ``IM``
specifies whether ADCIRC is solved in two-dimensional depth-integrated (2DDI) or
in three-dimensions (3D), solution of the governing equations is semi-implicit
or explicit in time, and whether the model formulation is barotropic or
baroclinic. Popular values for 2D barotropic ADCIRC include ``IM=0`` and
``IM=111112``; users should be aware that the
:ref:`A00 <A00>`, :ref:`B00 <B00>`, and :ref:`C00 <C00>` coefficients must be
specified differently in these two cases.

.. _default_im_values:

Default IM Values
-----------------

Default simulation option combinations can be specified through single or double
digit values, some of which are shortcuts to the six-digit codes described in
the next heading. The available ``IM`` values are specified in the table below
and in the following section on 6-digit values:

.. list-table::
   :widths: 10 10 40
   :header-rows: 1
   :class: wrap-table

   * - IM Value
     - Six-digit Equivalent
     - Description
   * - 0
     - 111111
     - Barotropic 2DDI
   * - 1
     - 611111
     - Barotropic 3D velocity-based momentum
   * - 2
     - -
     - Barotropic 3D stress-based momentum
   * - 10
     - -
     - Barotropic 2DDI with passive scalar transport
   * - 11
     - -
     - Barotropic 3D velocity-based momentum with passive scalar transport
   * - 20
     - -
     - Baroclinic 2DDI
   * - 21
     - -
     - Baroclinic 3D velocity-based momentum
   * - 30
     - -
     - Baroclinic 2DDI with passive scalar transport
   * - 31
     - -
     - Baroclinic 3D velocity-based momentum with passive scalar transport

Note that all default ``IM`` values employ the semi-implicit consistent GWCE
mass matrix solver. It has less numerical error and tends to be more stable than
the explicit mass-lumping approach at the expense of computational time and
memory.

.. _six_digit_im_codes:

Six-digit IM Codes
------------------

For fine-grained control of various options six-digit codes for ``IM`` can be
specified. Each digit represents a specific option regarding the dimension and
the formulation of certain terms or integration methods in the GWCE or momentum
equations. The available options for each digit are specified below, with the
first digit being the left-most. The internal flags that are set are listed to
help users dig through the code.

.. list-table::
   :header-rows: 1
   :widths: 10 15 15 15 15 15 15
   :class: wrap-table tighter-table

   * - Value
     - Digit 1: 2DDI/3D, Lateral Stress in GWCE [1]_
     - Digit 2: Advection in GWCE [2]_
     - Digit 3: Lateral Stress in Momentum [1]_
     - Digit 4: Advection in Momentum [2]_
     - Digit 5: Area Integration in Momentum
     - Digit 6: GWCE Mass Matrix, Barotropic/Baroclinic
   * - 1 (default)
     - 2DDI, Kolar-Gray flux-based ``CGWCE_LS_KGQ=.TRUE.``
     - Non conservative ``CGWCE_Advec_NC=.TRUE.``
     - Integration by parts, velocity-based ``CME_LS_IBPV=.TRUE.``
     - Non conservative ``CME_New_NC=.TRUE.``
     - Corrected ``CME_AreaInt_Corr=.TRUE.``
     - Consistent (implicit for linear part of gravity wave term), barotropic ``ILump=0``
   * - 2
     - 2DDI, 2-part flux-based ``CGWCE_LS_2PartQ=.TRUE.``
     - Conservative form 1 ``CGWCE_Advec_C1=.TRUE.``
     - Integration by parts, flux-based ``CME_LS_IBPQ=.TRUE.``
     - Conservative form 1 ``CME_New_C1=.TRUE.``
     - Original ``CME_AreaInt_Orig=.TRUE.``
     - Lumped (explicit), barotropic ``CGWCE_Lump=.TRUE.``, ``ILump=1``
   * - 3
     - 2DDI, 2-part velocity-based ``CGWCE_LS_2PartV=.TRUE.``
     - Conservative form 2 ``CGWCE_Advec_C2=.TRUE.``
     - Integration by parts, velocity-based symmetrical ``CME_LS_IBPSV=.TRUE.``
     - Conservative form 2 ``CME_New_C2=.TRUE.``
     - -
     - Consistent (implicit for full gravity wave term), barotropic ``CGWCE_HDP=.TRUE.``, ``IFNL_HDP=1``, ``ILump=0``
   * - 4
     - 2DDI, 2-part flux-based symmetrical ``CGWCE_LS_2PartSQ=.TRUE.``
     - -
     - Integration by parts, flux-based symmetrical ``CME_LS_IBPSQ=.TRUE.``
     - -
     - -
     - -
   * - 5
     - 2DDI, 2-part velocity-based symmetrical ``CGWCE_LS_2PartSV=.TRUE.``
     - -
     - 2 Part, velocity-based (*not implemented*) ``CME_LS_2PartV=.TRUE.``
     - -
     - -
     - -
   * - 6
     - 3D, Kolar-Gray flux-based ``C2DDI=.FALSE.``, ``CGWCE_LS_KGQ=.TRUE.``, ``C3D=.TRUE.``, ``C3DVS=.TRUE.``, ``ILump=0``
     - -
     - 2 Part, flux-based (*not implemented*) ``CME_LS_2PartQ=.TRUE.``
     - -
     - -
     - -

A common code combination is ``IM=111112``, which is identical to the default
``111111`` (same as ``IM=0``), but simulates in explicit mass-lumping mode. Note
that :ref:`A00 <A00>`, :ref:`B00 <B00>`, and :ref:`C00 <C00>` must be set to
``0.0 1.0 0.0`` when in this mode. Lumped explicit mode is a useful alternative
to the (default) semi-implicit consistent GWCE mass matrix mode, because the
latter requires a matrix solve that increases computational time and memory. By
comparison, the explicit mass-lumping mode is about twice as fast and scales to
fewer grid nodes per computational core. [3]_ Moreover, for model setups that
are sufficiently resolved in space and time, differences in the solution between
approaches should be small. Though, many users have reported somewhat lower
stability in lumped explicit mode.

The most recent version (55+) also has an option that improves the (default)
semi-implicit consistent GWCE mass matrix mode to compute the complete (total
depth) gravity wave term (free surface gradient) implicitly; toggled by setting
IMDigit-6 to 3. The default version (IMDigit-6=1), only computes the initial
still water depth component of the free surface gradient implicitly, which might
make it more susceptible to CFL violations in shallow depths and can encounter
Matrix diagonality issues overland where the initial still water depth is
negative.

References
----------

.. raw:: html

   <references />

.. [1]
   K.M. Dresback, R.L. Kolar, R.A. Luettich, Jr. (2005). On the Form of the
   Momentum Equation and Lateral Stress Closure Law in Shallow Water Modeling,
   in: Estuar. Coast. Model., American Society of Civil Engineers, Reston, VA,
   399–418. doi:10.1061/40876(209)23

.. [2]
   K.M. Dresback, R.L. Kolar, J.C. Dietrich (2005). On the Form of the Momentum Equation for Shallow Water Models Based on the Generalized Wave Continuity Equation: Conservative vs. Non-Conservative. Advances in Water Resources, 28(4), 345-358. doi:10.1016/j.advwatres.2004.11.011

.. [3]
   S. Tanaka, S. Bunya, J.J. Westerink, C. Dawson, R.A. Luettich, Scalability of an Unstructured Grid Continuous Galerkin Based Hurricane Storm Surge Model, J. Sci. Comput. 46 (2011) 329–358. doi:10.1007/s10915-010-9402-1
