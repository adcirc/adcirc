.. meta::
   :description: TAU0 in ADCIRC
   :keywords: adcirc, tau0

.. _tau0_parameter:

TAU0
====

**TAU0** is an input in the `fort.15 file <fort.15_file>`__ that influences
the degree of numerical diffusion in ADCIRC's governing equations. Specifically,
it influences the weighting factor that determines the relative contribution of
the primitive and wave portions of the the Generalized Wave-Continuity Equation
(GWCE). The weighting factor, :math:`\tau_0`, is affected by values in both the
`fort.15 file <fort.15_file>`__ and the `fort.13 file <fort.13_file>`__, if the
`primitive weighting in continuity
equation <primitive_weighting_in_continuity_equation>`__ or `min and max
primitive weighting in continuity
equation <min_and_max_primitive_weighting_in_continuity_equation>`__ nodal
attributes are specified. This page addresses both the ``TAU0`` value in the
fort.15 file and the :math:`\tau_0` parameter more generally.

.. _tau0_parameter_summary:

Parameter Summary
-----------------

Because the :ref:`TAU0` value specified in the fort.15 file can be either a flag
(indicating how ADCIRC should operate) or the actual :math:`\tau_0` value used
in solving the GWCE, it is important to distinguish between the two. All
negative :ref:`TAU0` are flags, all positive :ref:`TAU0` are :math:`\tau_0`.
Furthermore, some values are overridden by the :ref:`primitive weighting in
continuity equation <primitive_weighting_in_continuity_equation>` nodal
attribute. The following table is a summary of possible ``TAU0`` values and
their meaning. Note that for ``TAU0 = -x.1`` (where ``x`` is an integer),
behavior is the same as ``-x``, but the :math:`\tau_0` values are written to the
:ref:`fort.90 file <fort90>`. More on this below in
:ref:`Outputting <tau0_outputting>`.

.. list-table::
   :header-rows: 1
   :class: wrap-table tighter-table
   :widths: 3 1 1 1 1 1 1 1

   * - fort.15 ``TAU0``
     - >= 0
     - -1
     - -2
     - -3
     - -5
     - -6
     - -7
   * - Varies in space
     - no (mostly)
     - yes
     - yes
     - yes
     - yes
     - yes
     - yes
   * - Varies in time
     - no
     - no
     - no
     - yes
     - yes
     - yes
     - yes
   * - Space-averaged
     - no
     - no
     - no
     - yes
     - yes
     - yes
     - yes
   * - Time-averaged
     - no
     - no
     - no
     - no
     - no
     - yes
     - yes
   * - primitive weighting in continuity equation
     - yes
     - discouraged
     - discouraged
     - yes
     - discouraged
     - yes
     - yes
   * - min and max primitive weighting in continuity equation
     - no
     - no
     - no
     - no
     - yes
     - no
     - no
   * - Minimum
     - N/A
     - 0.002
     - 0.002
     - N/A
     - ``Tau0FullDomainMin``
     - N/A
     - N/A
   * - Maximum
     - N/A
     - 0.005
     - 1
     - 0.2
     - ``Tau0FullDomainMax``
     - 0.2
     - 0.2
   * - Code Flags
     -
     -
     -
     -
     -
     -
     -
   * - ``HighResTimeVaryingTau0``
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .TRUE.
     - .FALSE.
     - .TRUE.
     - .TRUE.
   * - ``FullDomainTimeVaryingTau0``
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .TRUE.
     - .FALSE.
     - .FALSE.
   * - ``TimeAveragedTau0``
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .TRUE.
     - .FALSE.
   * - ``BackLoadedTimeAveragedTau0``
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .FALSE.
     - .TRUE.

Nodal Attributes
~~~~~~~~~~~~~~~~

The :ref:`primitive weighting in continuity
equation <primitive_weighting_in_continuity_equation>` nodal attribute permits
a spatially variable :math:`\tau_0` that is (at least initially) set equal to
the nodal attribute values.

-  ``TAU0 >= 0``, :math:`\tau_0` is constant in time and equal to :ref:`primitive
   weighting in continuity
   equation <primitive_weighting_in_continuity_equation>` values
-  ``TAU0 = -3, -6, or -7``, :math:`\tau_0` varies in time and based on
   :ref:`primitive weighting in continuity
   equation <primitive_weighting_in_continuity_equation>` values

The :ref:`min and max primitive weighting in continuity
equation <min_and_max_primitive_weighting_in_continuity_equation>` only works
with ``TAU0 = -5``, details on it are below.

.. _positive_tau0:

Positive TAU0
~~~~~~~~~~~~~

A positive value is ignored if the :ref:`primitive weighting in continuity
equation <primitive_weighting_in_continuity_equation>` nodal attribute is
specified in the fort.15 file. For all positive values, the value specified in
the fort.15 file is spatially and temporally constant and applied directly, i.e.
``TAU0 =``\ :math:`\tau_0`. Practical guidance on setting a constant ``TAU0`` is
provided below in :ref:`Selecting Values <selecting_values>`. Mathematically, the
GWCE is a pure wave equation for :math:`\tau_0=0`, and approaches a pure
primitive equation as :math:`\tau_0\rightarrow\infty`, however in practice, it
behaves like a pure primitive equation once :math:`\tau_0` reaches 1 or so.

.. _negative_tau0:

Negative TAU0
~~~~~~~~~~~~~

| Unexpected behavior may arise if a mismatching nodal attribute is specified
  alongside one of these values, and this is discouraged.
| *Spatially varying, temporally constant*
| \*\ ``TAU0 = -1``, :math:`\tau_0` is computed as follows:

-  

   -  ``If depth >= 10; then``\ :math:`\tau_0`\ ``= 0.005``
   -  ``If depth < 10; then``\ :math:`\tau_0`\ ``= 0.020``
   -  This value is ignored if `primitive weighting in continuity
      equation <primitive_weighting_in_continuity_equation>`__ is specified in
      the fort.15 file.

-  ``TAU0 = -2``, :math:`\tau_0` is is computed according to depth as follows:

   -  ``If depth >= 200; then``\ :math:`\tau_0`\ ``= 0.005``
   -  ``If 1 < depth < 200; then``\ :math:`\tau_0`\ ``= 1/depth``
   -  ``If depth < 1; then``\ :math:`\tau_0 = 1`
   -  This value is ignored if `primitive weighting in continuity
      equation <primitive_weighting_in_continuity_equation>`__ is specified in
      the fort.15 file.

| *Spatially and temporally varying*
| In the unlikely event that one of the below options not paired with its
  corresponding nodal attribute (the code will exit if this is done for some
  ``TAU0``), then default values are set using the ``TAU0 = -1`` logic.

-  ``TAU0 = -3``, :math:`\tau_0` is computed from ``TAU0Base`` (read in from the
   `nodal attribute <primitive_weighting_in_continuity_equation>`__) as follows:

   -  ``If TAU0Base < 0.025; then``\ :math:`\tau_0`\ ``= TAU0Base`` (constant in
      time)
   -  ``If TAU0Base >= 0.025; then``\ :math:`\tau_0`\ ``= TAU0Base+1.5*TK``
      where ``TK=speed*Cd/depth``

-  ``TAU0 = -5``, :math:`\tau_0` is computed similar to ``TAU0 = -3`` as
   follows:

   -  :math:`\tau_0`\ ``= Tau0Min+1.5*TK``
   -  It is limited by
      ```Tau0FullDomainMin`` <Tau0FullDomainMin>`__\ ``<=``\ :math:`\tau_0`\ ``<=``\ ```Tau0FullDomainMax`` <Tau0FullDomainMax>`__,
      which are specified on the line following ``TAU0`` in the fort.15 file
      (only when ``TAU0 = -5``)
   -  If the `min and max primitive weighting in continuity
      equation <min_and_max_primitive_weighting_in_continuity_equation>`__ nodal
      attribute is used, then the nodal minimum and maximum values replace the
      full-domain values in the above calculation.

-  ``TAU0 = -6``, :math:`\tau_0` is defined using the rules for ``TAU0 = -3``,
   and then is set equal to the time-average of the current and previous
   (time-averaged) values.
-  ``TAU0 = -7``, :math:`\tau_0` is defined using the rules for ``TAU0 = -3``,
   and then is set equal to the weighted time-average of the current and
   previous (time-averaged) values as follows:

   -  :math:`\tau_0`\ ``= AlphaTau0*TAU0VAR + (1.0-AlphaTau0)*LastTau0``, where
      ``TAU0VAR`` is ``Tau0Base`` after spatial averaging, and
      ``AlphaTau0 = 0.25`` is hard-coded into the model. This means that
      :math:`\tau_0` is weighted 75% toward older values.

.. _spatial_and_temporal_updating:

Spatial and Temporal Updating
-----------------------------

For ``TAU0 = -3, -5, -6, or -7``, :math:`\tau_0` is updated (via the
``CalculateTimeVaryingTau0`` subroutine) in space and in time. An initial
"update" is done when the model starts. After this, updates are done only after
a time step in which there is a change in wet/dry state somewhere in the model
domain. For use cases that contain large number of nodes near the wet/dry
boundary, this can be the equivalent of updating every time step. However, for
use cases that have little or no wet/dry changes, **there may be little or no
updating**. The rules listed above in `TAU0 Values <TAU0#TAU0_Values>`__ are
applied during the update. Each node's :math:`\tau_0` is then spatially averaged
with all immediate neighbors. Time-averaging (for ``TAU0 = -6 or -7``) is
applied last.

.. _selecting_values:

Selecting Values
----------------

``TAU0 = -3``, paired with the `primitive weighting in continuity
equation <primitive_weighting_in_continuity_equation>`__ nodal attribute is
generally the most popular formulation. In this case, ``TAU0Base`` nodal
attribute values can be generated with the ADCIRC utility program tau0_gen.f.
The program bases generation on the following logic applied to each node
individually:

-  ``If the avg. dist. to neighboring nodes < 1750 m; then TAU0Base = 0.03``
-  Otherwise

   -  ``If depth < 10m; then TAU0Base = 0.02`` (``TAU0`` is constant in time)
   -  ``If depth > 10m; then TAU0Base = 0.005`` (``TAU0`` is constant in time)

:math:`\tau_0` needs to be smaller in deeper water where there is little
dissipation, and can also be sensitive to mesh resolution. The spatial
variation, spatial smoothing, and physics-driven time-updating typically allow
for a good balance between stability and conservation.

For manually-specified positive values (``TAU0=``\ :math:`\tau_0`), a good rule
of thumb for setting ``TAU0`` is to set it equal to the largest value of an
equivalent linear friction factor: for linear friction
``TAU0 =``\ ```TAU`` <TAU>`__; for quadratic friction
``TAU0 = max(speed*Cd/depth)``. Typical values for ``TAU0`` are in the range of
0.005 – 0.1.

.. _tau0_outputting:

Outputting
----------

For ``TAU0`` formulations that vary spatially or temporally, ADCIRC can output
the internally-calculated nodal :math:`\tau_0` values. They are written to the
:ref:`fort.90 file <fort90>`, which has the same format and output frequency
as the water surface elevation output file (:ref:`fort.63 <fort63>`). fort.90
output is activated by placing a 1 in the tenths place of ``TAU0`` in the
fort.15 file. For example, if ``TAU0 = -3.1``, the calculation of :math:`\tau_0`
is still carried out according to the description of ``TAU0 = -3`` above, and
the fort.90 output file will also be produced.
