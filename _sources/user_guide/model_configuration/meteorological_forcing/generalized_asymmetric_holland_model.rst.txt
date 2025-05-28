.. meta::
   :description: Generalized Asymmetric Holland Model in ADCIRC
   :keywords: adcirc, generalized asymmetric holland model

.. _gahm:

Generalized Asymmetric Holland Model
====================================

The Generalized Asymmetric Holland Model (hereinafter, GAHM) is a parametric
hurricane vortex model developed in ADCIRC for operational forecasting purpose,
based on the classic Holland Model  [1]_ (hereinafter, HM). The original HM was
developed to render an ideal symmetric hurricane (a.k.a., annular hurricane). To
represent hurricanes that may exhibit asymmetric structures, the Asymmetric
Holland Model (hereinafter, AHM) was developed and implemented in ADCIRC for
more practical usages [2]_, which has the same set of equations as of the HM,
but takes azimuthally-varying radius to the maximum wind to reconstruct its
spatial pressure and wind fields. Both the HM and AHM suffer from flaws that
make rendering large but weak storms difficult. As a result, efforts were made
to develop a more generalized model to fulfill the purpose. Compared to the HM
and the AHM, the GAHM removes the assumption of cyclostrophic balance at the
radius of maximum wind during derivation of its equations, and allows for a
better representation of a wide range of hurricanes. Another important feature
of the GAHM is the introduction of a composite wind method, which when activated
enables the usage of multiple storm isotaches. It makes representing complex
hurricane structures possible.

.. _the_classic_holland_model:

The Classic Holland Model
-------------------------

The HM is an analytic model that describes the radial pressure and wind profiles
of a standard hurricane. To begin with, Holland found that the normalized
pressure profiles of a number of hurricanes resemble a family of rectangular
hyperbolas and may be approximated by a hyperbolic equation, which after
antilogarithms and rearranging yields the radial pressure equation:

.. raw:: html

   <center>

:math:`P(r) = P_c + (P_n - P_c)e^{-A/r^B} \quad` (1)

.. raw:: html

   </center>

where :math:`P_c` is the central pressure, :math:`P_c` is the ambient pressure
(theoretically at infinite radius), :math:`P(r)` is the pressure at radius
:math:`r` from the center of the hurricane, and :math:`A` and :math:`B` are
shape parameters that may be empirically estimated from observations in a
hurricane.

Substituting (1) into the gradient wind equation, which describes a steady flow
balanced by the horizontal pressure gradient force, the centripetal
acceleration, and the Coriolis acceleration for a vortex above the influence of
the planetary boundary layer where the atmospheric flow decouples from surface
friction  [3]_, gives the radial wind equation of a hurricane:

.. raw:: html

   <center>

.. math::

   V_g(r) = \sqrt{AB(P_n - P_c)e^{-A/r^B}/\rho r^B + (\frac{rf}{2})^2} - \frac{rf}{2} \quad (2)

.. raw:: html

   </center>

where :math:`V_g(r)` is the gradient wind at radius :math:`r`, :math:`\rho` is
the density of air, :math:`f` is the Coriolis parameter. In the region of the
maximum winds, if we assume that the Coriolis force is negligible in comparison
to the pressure gradient and centripetal force, then the air is in cyclostrophic
balance. By removing the Coriolis term in (2) we get the cyclostrophic wind

.. raw:: html

   <center>

.. math::

   V_c(r) = \sqrt{AB(P_n - P_c)e^{-A/r^B}/\rho r^B} \quad (3)

.. raw:: html

   </center>

By setting :math:`dV_c/dr = 0` at radius to the maximum wind
:math:`r = R_{max}`, it is obtained that

.. raw:: html

   <center>

.. math::

   A = (R_{max})^B \quad (4)

.. raw:: html

   </center>

Thus the (:math:`R_{max}`) is irrelevant to the relative value of ambient and
central pressures, and is solely defined by the shape parameters :math:`A` and
:math:`B`. Substituting (4) back into (3) to get rid of :math:`A`, we get an
estimate of :math:`B` as a function of the maximum wind speed

.. raw:: html

   <center>

.. math::

   B = (V_{max})^2\rho e/(P_n - P_c) \quad (5)

.. raw:: html

   </center>

It was notable that the maximum wind speed is proportional to the square root of
:math:`B` and irrespective of the (:math:`R_{max}`), given a constant pressure
drop. It was also reasoned by Holland that a plausible range of :math:`B` would
be between 1 and 2.5 for realistic hurricanes. Substituting (4) and (5) back
into (1) and (2) yields the final radial pressure and wind profiles for the HM

.. raw:: html

   <center>

.. math::

   P(r) = P_c + (P_n - P_c)e^{-(R_{max}/r)^B} \quad (6)

.. raw:: html

   </center>

.. raw:: html

   <center>

.. math::

   V_g(r) = \sqrt{(V_{max})^2e^{1-(R_{max}/r)^B}(R_{max}/r)^B + (\frac{rf}{2})^2} - \frac{rf}{2} \quad (7)

.. raw:: html

   </center>

When sparse observations of a hurricane are given, estimates of the
:math:`R_{max}` and shape parameter :math:`B` can be estimated by fitting data
into the radial wind equation, which in turn allow us to compute :math:`P(r)`
and :math:`V_g(r)` along the radius :math:`r` of the hurricane. However,
discrepancies between wind observations and computed winds were sometimes found,
and were negatively correlated to the Rossby number at :math:`r = R_{max}`,
defined as

.. raw:: html

   <center>

.. math::

   R_o = \frac{Nonlinear Acceleration}{Coriolis force} \approx \frac{V_{max}^2/R_{max}}{V_{max}f} = \frac{V_{max}}{R_{max}f} \quad (8)

.. raw:: html

   </center>

By definition, a large :math:`R_o (\approx 10^3)` describes a system in
cyclostrophic balance that is dominated by the inertial and centrifugal force
with negligible Coriolis force, such as a tornado or the inner core of an
intense hurricane, whereas a small value :math:`(\approx 10^{-2} \sim 10^2)`
signifies a system in geostrophic balance where the Coriolis force plays an
important role, such as the outer region of a hurricane. As a result, the
assumption of cyclostrophic balance at :math:`R_{max}` made in HM is mostly
valid for describing an intense and narrow (small :math:`R_{max}`) hurricane
with a large :math:`R_o`, but not applicable for a weak and broad hurricane with
a small :math:`R_o`. This intrinsic problem with the HM calls our intention to
develop a generalized model that will work consistently for a wide range of
hurricanes, which theoretically can be accomplished by removing the above
cyclostrophic balance assumption and re-derive the radial pressure and wind
equations (6)&(7).

.. _derivation_of_the_gahm:

Derivation of the GAHM
----------------------

The GAHM also starts with the same radial pressure and wind equations (1)&(2)
with shape parameters :math:`A` and :math:`B`, as in the HM. Without assuming
cyclostrophic balance at :math:`R_{max}`, we take :math:`dV_g/dr = 0` at
:math:`r = R_{max}` to get the adjusted shape parameter :math:`B_g` as

.. raw:: html

   <center>

.. math::
   B_g = \frac{(V_{max}^2 + V_{max}R_{max}f)\rho e^\varphi}{\varphi(P_n - P_c)} = B \frac{(1+1/R_o)e^{\varphi - 1}}{\varphi} \quad (9)

.. raw:: html

   </center>

where :math:`{\varphi}` is a scaling parameter introduced to simplify the
derivation process, defined as

.. raw:: html

   <center>

.. math::
   \varphi = \frac{A}{R_{max}^B} \quad \text{or} \quad A = \varphi R_{max}^B \quad (10)

.. raw:: html

   </center>

and later derived as

.. raw:: html

   <center>

.. math::

   \varphi = 1 + \frac{V_{max}R_{max}f}{B_g(V_{max}^2+V_{max}R_{max}f)} = 1 + \frac{1/R_o}{B_g(1+1/R_o)}  \quad (11)

.. raw:: html

   </center>

Thus, the :math:`R_{max}` in the GAHM is not entirely defined by the shape
parameters :math:`A` and :math:`B` as in the HM, but also by the scaling factor
:math:`{\varphi}`, as Equation (11) indicates that :math:`{\varphi} \ge 1`.
Numerical solutions for :math:`B_g` and :math:`{\varphi}` can be solved
iteratively in the model using Equation (9)&(11). Figure 1 illustrates how
:math:`B_g/B` and :math:`\varphi` vary with :math:`\log_{10}R_o` given different
:math:`B` values. It is evident that values of both :math:`B_g/B` and
:math:`\varphi` remain close to 1 when :math:`\log_{10}R_o` is within the range
of [1,2], but increase noticeably as :math:`\log_{10}R_o` decreases below 1, and
the smaller the value of :math:`B`, the bigger the changes.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_1._GAHM.png
   :alt: Profiles of Bg⁄B (left panel) and φ (right panel)


.. raw:: html

   <center>

Figure 1. Profiles of :math:`B_g/B` (left panel) and :math:`\varphi` (right
panel) with respect to :math:`\log_{10}R_o`, given different :math:`B` values as
shown in different colors.

.. raw:: html

   </center>

Substituting (9)&(11) back into (1)&(2) yields the final radial pressure and
wind equations for the GAHM

.. raw:: html

   <center>

.. math::
   P(r) = P_c + (P_n - P_c)e^{-\varphi(R_{max}/r)^{B_g}} \quad (12)

.. raw:: html

   </center>

.. raw:: html

   <center>

.. math::
   V_g(r) = \sqrt{V_{max}^2(1+1/R_o)e^{\varphi(1-(R_{max}/r)^{B_g})}(R_{max}/r)^{B_g} + (\frac{rf}{2})^2} - \frac{rf}{2} \quad (13)

.. raw:: html

   </center>

Influence of the Coriolis force on the radial pressure and wind profiles are
evidenced by the presence of :math:`R_o` and :math:`\varphi` in (12)&(13). A
special case scenario is when we set :math:`f=0`, which corresponds to an
infinitely large :math:`R_o`, then (12)&(13) in the GAHM reduce to (6)&(7) in
the HM. However，for a hurricane with a relatively small :math:`R_o`, the
influence of the Coriolis force can only be addressed by the GAHM. It meets our
expectation that the GAHM’s solution approaches to that of the HM’s when the
influence of Coriolis force is small, but departs from it when the Coriolis
force plays an important role in the wind system.

The above reasoning can be demonstrated by the 3D plots in Figure 2, which show
the normalized gradient winds of the HM (left panel) and the GAHM (right panel)
as functions of the normalized radial distances :math:`r/R_{max}`, the Holland
:math:`B` parameter, and :math:`R_o`. In both panels, each colored surface
represents the normalized gradient winds corresponding to a unique Holland B
value. By definition, we get :math:`V_g = V_{max}` at :math:`r = R_{max}`, which
means all the surfaces in each panel should intersect with the plane of
:math:`r/R_{max} = 1` on the plane of :math:`V_g/V_{max} = 1`, no matter what
values of :math:`R_o`. However, the line of intersection (shown by the black
line) shown in the left panel deviates from the plane of :math:`V_g/V_{max} =1`
as :math:`\log_{10}R_o` decreases from 2 to close to 0 (:math:`R_o` decreases
from 100 to 1), while remains on the plane regardless of how :math:`R_o` changes
in the right panel, demonstrating that the GAHM is mathematically more coherent
than the HM.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_2._GAHM.png
   :alt: Comparison of the normalized gradient winds


.. raw:: html

   <center>

Figure 2. The normalized gradient wind profiles of the HM (left panel) and the
GAHM (right panel) as functions of the normalized radial distances and
:math:`\log_{10}R_o`, given different Holland :math:`B` values.

.. raw:: html

   </center>

To have a dissective look of the surface plots in Figure 2, we draw slices
perpendicular to the axis of :math:`\log_{10}R_o` at three different values 0,
1, 2, and plot the lines of intersection with each surface in Figure 3. It is
evident that we get :math:`V_g = V_{max}` at :math:`r = R_{max}` consistently in
the right panel for the GAHM regardless of the value of :math:`R_o`. The HM in
the left panel, however, generates distorted wind profiles with underestimated
maximum winds skewed inward towards the storm center, espeically when
:math:`\log_{10}R_o < 1`. As a results, when both models being applied to real
hurricane cases, the GAHM will perform more consistently than the HM.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_3._GAHM.png
   :alt: Normalized gradient wind profiles


.. raw:: html

   <center>

Fig 3. Slices of the normalized gradient wind profiles (as shown in Figure 2) at
:math:`\log_{10}R_o =0, 1, 2` (or correspondingly :math:`R_o =1, 10, 100`).

.. raw:: html

   </center>

.. _calculation_of_the_radius_to_the_maximum_wind:

Calculation of the Radius to the Maximum Wind
---------------------------------------------

Same with the HM and AHM, the GAHM also uses processed forecast advisories
(during active hurricanes) or best track advisories (post-hurricanes) from the
National Hurricane Center (NHC) in ATCF format as input files, which contain a
time series of storm parameters (usually at 6-hour intervals) such as storm
location, storm movement, central pressure, 1 minute averaged maximum wind,
radii to the 34-, 50-, and/or 64-kt storm isotaches in 4 storm quadrants (NE,
SE, SW, NW), etc. See meteorological input file with NWS = 20 for more details.
As a standard procedure, the :math:`B_g` and :math:`R_{max}` are pre-computed in
4 storm quadrants for all available isotaches in the ASWIP program  [4]_ and
appended to the input file prior to running an ADCIRC simulation. The following
describes the procedures to prepare the input file for the GAHM.

First, the influence of the boundary layer effect must be removed to bring the
maximum sustained wind and the 34-, 50-, and/or 64-kt isotaches from 10 meter
height to the gradient wind level. Practically, the maximum gradient wind can be
directly calculated as

.. raw:: html

   <center>

.. math::
   V_{max} = \vert \frac{\overrightarrow{V_M} - \gamma\overrightarrow{V_T}}{W_{rf}} \vert = \frac{V_M - \gamma V_T}{W_{rf}} \quad (14)

.. raw:: html

   </center>

where :math:`\overrightarrow{V_M}` is the reported maximum sustained wind at 10
meter height assuming in the same direction as :math:`\overrightarrow{V_T}`,
:math:`\overrightarrow{V_T}` is the storm translational speed calculated from
successive storm center locations, :math:`W_{rf} = 0.9` is the wind reduction
factor for reducing wind speed from the gradient wind level to the surface at 10
meter height (Powell et al., 2003), and :math:`\gamma` is the damp factor for
:math:`V_T`. The following formula of :math:`\gamma` is employed in the ASWIP
program:

.. raw:: html

   <center>

.. math::
   \gamma = \frac{V_g}{V_{max}} \quad (15)

.. raw:: html

   </center>

which is the ratio of gradient wind speed to the maximum wind speed along a
radial wind profile. Thus, :math:`\gamma` is zero at storm center, and increases
with :math:`r` until reaches a maximum value of 1 at :math:`R_{max}`, then
gradually decreases outward to zero.

In addition to the scalar reduction in wind speed, surface friction and
continuity also cause the vortex wind to flow inward across isobars, with an
inward rotation angle :math:`\beta` according to the Queensland Government's
Ocean Hazards Assessment (2001)  [5]_:

.. raw:: html

   <center>

.. math::
   \beta = \begin{cases}
   10^{\circ}, & r<R_{max} \\
   10^{\circ} + 75(r-R_{max})/R_{max}, & R_{max} \le r<1.2R_{max} \\
   25^{\circ}, & r \ge 1.2R_{max}
   \end{cases} \quad   (16)

.. raw:: html

   </center>

Thus, the gradient wind at the radii to specified storm isotaches in 4 storm
quadrants can be obtained from the observed isotaches as

.. raw:: html

   <center>

.. math::
   \begin{align}
   V_r & = \vert \overrightarrow{V_r}\vert = \vert \overrightarrow{V_{inflow}}\vert \\
   & = \frac{\vert\overrightarrow{V_{isot}} - \gamma\overrightarrow{V_T}  \vert}{W_{rf}}
   \end{align} \quad (17)

.. raw:: html

   </center>

where :math:`\overrightarrow{V_{isot}}` is the observed isotach wind speed with
an unknown angle :math:`\varepsilon`, and :math:`\overrightarrow{V_{inflow}}` is
the wind speed at radius to specified isotach before the inward rotation angle
:math:`\beta` is removed.

Rewriting (17) in x- and y-components yields:

.. raw:: html

   <center>

.. math::
   V_r\cos(quad(i)+90+\beta)=V_{isot}\cos(\varepsilon)-\gamma {\mu}_T \quad (18)

.. raw:: html

   </center>

.. raw:: html

   <center>

.. math::
   V_r\sin(quad(i)+90+\beta)=V_{isot}\sin(\varepsilon)-\gamma {\nu}_T \quad (19)

.. raw:: html

   </center>

where :math:`quad(i)` is the azimuth angle of the :math:`i-th` storm quadrant
(NE, SE, SW, NW at :math:`45^\circ, 135^\circ, 225^\circ, 315^\circ`,
respectively), :math:`V_{isot}\cos(\varepsilon)` and
:math:`V_{isot}\sin(\varepsilon)` are the zonal and meridional components of
:math:`\overrightarrow{V_{isot}}`, :math:`{\mu}_T` and :math:`{\nu}_T` are the
zonal and meridional components of :math:`\overrightarrow{V_T}`.

Given an initial guess of :math:`R_{max}`, values of :math:`B_g` and
:math:`\varphi` can be solved iteratively from (9) and (11) until both converge,
and :math:`V_r` can be estimated by combining (15), (17), (18), and (19).
Plugging :math:`V_{max}` from (14), the above calculated
:math:`B_g, \varphi, V_{max}, V_r` and the observed radius :math:`R_r` to
:math:`V_r` back into (13), a new :math:`R_{max}` can be inversely solved by a
root-finding algorithm. Since the above calculations are carried out based on an
initial guess of :math:`R_{max}`, wWe need to repeat the entire process until
:math:`R_{max}` converges.

In case where multiple isotaches are given in the forecast/best track
advisories, the :math:`R_{max}` for the highest isotach will be calculated using
the above procedure, and used as the pseudo :math:`R_{max}` for the entire storm
(physically, there is only one :math:`R_{max}` found along a radial wind profile
). For each lower isotach, :math:`R_{max}` will be calculated with the pseudo
:math:`R_{max}` set as its initial value to determine the inward rotation angle
:math:`\beta` following the above process only once. The use of the pseudo
:math:`R_{max}` across all storm isotaches ensures that the cross-isobar
frictional inflow angle changes smoothly along the radius according to (17).

Occasionally, we have to deal with situations where :math:`V_{max} < V_r`, which
violate (13) so :math:`R_{max}` couldn't be calculated. These situations mostly
happen in the right hand quadrants (in the Northern Atmosphere) of a weak storm
with a relatively high translational speed. For cases like this, we assign
:math:`V_{max} = V_r`, which is equivalent to assigning :math:`R_{max} = R_r`.

After the ASWIP program finishes processing the input file, it can be readily
used by the GAHM to construct spatial pressure and wind fields in ADCIRC for
storm surge forecast.

.. _composite_wind_generation:

Composite Wind Generation
-------------------------

Since storm parameters are only given in 4 storm quadrants (assuming at
:math:`45^\circ, 135^\circ, 225^\circ, 315^\circ` azimuthal angles,
respectively) at 3 available isotaches in the input file, spatial interpolation
of storm parameters must take place first at each ADCIRC grid node.
Traditionally, the single-isotach approach is used by the AHM, in which storm
parameters will be interpolated azimuthally from the highest isotach only. To
take advantage of the availability of multiple isotaches, a new composite wind
method is introduced in the GAHM, the multiple-isotach approach, in which storm
parameters will be interpolated both azimuthally and radially from all available
isotaches.

To begin, the relative location of a node to the storm center at given time
:math:`t` is calculated, specified by the azimuth angle :math:`\theta` and
distance :math:`d`. The angle :math:`\theta` places the node between two
adjacent quadrants :math:`i` and :math:`i+1`, where
:math:`quad(i) < \theta \le quad(i+1)`. For each storm parameter :math:`P` to be
interpolated, its value at :math:`(\theta,d)` are weighted between its values at
two pseudo nodes :math:`(quad(i),d)` and :math:`(quad(i+1),d)`:

.. raw:: html

   <center>

.. math::
   P(\theta,d)=\frac{P(quad(i),d)(90-\theta)^2+P(quad(i+1),d)\theta^2}{(90-\theta)^2+\theta^2} \quad (20)

.. raw:: html

   </center>

The distance :math:`d` then places each pseudo node between the radii of two
adjacent isotaches in its quadrant, and the value at the pseudo node is
interpolated using the inverse distance weighting (IDW) method:

.. raw:: html

   <center>

.. math::
   P(quad,d) = f_{34}P_{34}+f_{50}P_{50}+f_{64}P_{64} \quad (21)

.. raw:: html

   </center>

where :math:`P_{34}, P_{50}, P_{64}` are parameter values computed from the 34-,
50-, and 64-isotach, :math:`f_{34}, f_{50}, f_{64}` are distance weighting
factors for each isotach, calculated as

.. raw:: html

   <center>

.. math::
   \begin{array}{lll}
   \mathrm{I}.  &  r<R_{64}&  f_{64}=1,f_{50}=0,f_{34}=0 \\
   \mathrm{II}. & R_{64}\le r<R_{50} &  f_{64}=(r-R_{64})/(R_{50}-R_{64}),f_{50}=(R_{50}-r)/(R_{50}-R_{64}), f_{34}=0 \\
   \mathrm{III}. & R_{50}\le r<R_{34} &  f_{64}=0,f_{50}=(r-R_{50}/(R_{34}-R_{50})),f_{34}=(R_{34}-r/(R_{34}-R_{50})) \\
   \mathrm{IV}. & r \le R_{34} & f_{64}=0,f_{50}=0,f_{34}=1
   \end{array} \quad (22)

.. raw:: html

   </center>

and :math:`f_{34}+f_{50}+f_{64}=1`.

The above procedure is performed at each node of an ADCIRC grid. After all storm
parameters are interpolated, the pressure and gradient winds can be calculated
using (12)&(13). To bring the gradient wind down to the standard 10 meter
reference level, the same wind reduction factor :math:`W_{rf}` is applied, and
the tangential winds are rotated by an inward flow angle β according to (16).
Then, the storm translational speed is added back to the vortex winds. Last but
not least, a wind averaging factor is applied to convert resulted wind field
from 1-min to 10-min averaged winds in order to be used by ADCIRC. This new
composite wind method is simple and efficient, and more importantly, it assures
that the constructed surface winds match all observed storm isotaches provided
in NHC’s forecast or “best track” advisories.

.. _case_studies:

Case Studies
------------

Preliminary evaluation of the GAHM was carried out based on seven hurricanes
that struck the Gulf of Mexico and the Eastern United States: Katrina (2005),
Rita (2005), Gustav (2008), Ike (2008), Irene (2011), Isaac (2012), and Sandy
(2012), see Table 1. Ranging from category 1 to 5 on the Saffir-Simpson
Hurricane Wind Scale, these storms vary in storm track, forward motion, size,
intensity, and duration, but all caused severe damages to coastal states due to
destructive winds, wind-induced storm surges, and ocean waves. Their “best
track” advisories were retrieved from NHC’s ftp site
(ftp://ftp.nhc.noaa.gov/atcf; previous years’ data are located in the archive
directory) and pre-processed using the ASWIP program. The “best track” file
contains an estimate of the radius to the maximum wind for each data entry, but
will solely be used for model validation purpose as both the GAHM and AHM
calculate their own spatially-varying :math:`R_{max}`.

.. list-table:: Table 1. Seven selected hurricanes used for preliminary evaluation of the GAHM
   :class: wrap-table scroll-table
   :header-rows: 1
   :widths: 10 10 10 10 15

   * - Hurricane
     - Saffir-Simpson Wind Scale
     - Maximum Sustained Wind (knot)
     - Minimum Central Pressure (mbar)
     - Period from Formation to Dissipation
   * - **Katrina**
     - 5
     - 150
     - 902
     - 08/23-08/30, 2005
   * - **Rita**
     - 5
     - 150
     - 902
     - 09/18-09/26, 2005
   * - **Gustav**
     - 4
     - 135
     - 941
     - 08/23-09/04, 2008
   * - **Ike**
     - 4
     - 125
     - 935
     - 09/01-09/14, 2008
   * - **Irene**
     - 3
     - 105
     - 942
     - 08/21-08/30, 2011
   * - **Isaac**
     - 1
     - 70
     - 965
     - 08/21-09/03, 2012
   * - **Sandy**
     - 3
     - 95
     - 940
     - 10/22-10/01, 2012

Besides the maximum wind speed, both Holland :math:`B` and :math:`R_o` can be
used as key storm characteristics to characterize the development of the storm.
Figure 4 depicts the change of :math:`V_M`, :math:`B`, and :math:`\log_{10}R_o`
during different stages of the hurricanes along their best tracks. Typically,
both :math:`B` and :math:`R_o` increase as hurricane strengthens, and decrease
as hurricane dissipates, within the range of (0, 2.5). Previous analytical
evaluation has demonstrated that the GAHM behaves consistently better than the
HM, especially under situations where :math:`\log_{10}R_o < 1`. Here, evaluation
of model performance will be carried out by comparing the modeled winds with the
observed winds in the "best track" data, as well as the AHM, the SLOSH (Sea,
Lake, and Overland Surges from Hurricanes) winds, re-analysis H*Wind and
hindcast OWI modeled winds. The OWI winds and H*Winds are considered more mature
wind products that are able to resolve more complex structures of a hurricane
than the simple vortex models do.

.. list-table::
   :header-rows: 0

   * - .. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_4._GAHM.png

     - .. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_5._GAHM.png
     
     - .. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_6._GAHM.png

.. raw:: html

   <center>

Figure 4. The development of (a) The Maximum Wind Speed, (b) Holland MATH 0 ,
and (c) MATH 1 along the best tracks of 7 selected hurricanes

.. raw:: html

   </center>

.. _the_ahm_vs._the_gahm:

The AHM vs. The GAHM
--------------------

-  **Comparison of Radial Wind Profiles using the Single-Isotach Approach**

Since the AHM is an advanced version of the HM, here we only use model results
from the AHM for comparisons with the GAHM. First, the single-isotach approach
was evaluated using Hurricane Irene (2011) as an example. Figure 5 gives the
comparison of radial wind profiles of Hurricane Irene (2011) between the AHM and
the GAHM using the single-isotach approach at three snapshots, each representing
the developing (top panels), mature (middle panels), and dissipating (bottom
panels) stages of the hurricane.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_7._GAHM.png
   :alt: Radial wind profiles of Irene (2011) at three different stages.


.. raw:: html

   <center>

Figure 5. Comparison of radial wind profiles of Irene (2011) at three different
stages between the AHM and the GAHM.

.. raw:: html

   </center>

The cross-section radial winds from SW to NE are given in the left panels, and
NW to SE in the right panels. The observed isotaches at radii to specified
isotaches given in the "best track" file are also plotted as vertical line
segments for reference (highest isotach in black and lower isotaches in gray).
For a perfect match between the modeled winds and the isotaches, the radial
profiles should meet the tip of the line segments at the exact same height. The
:math:`B`, :math:`B_g` and :math:`\log_{10}R_o` are also computed at the same
snapshots in all 4 quadrants, given by Table 2.

.. table:: Table 2. Key storm characteristics :math:`B`, :math:`B_g` and :math:`\log_{10}R_o` at three snapshots of Irene (2011)
   :class: wrap-table scroll-table tight-table small-table
   :width: 100%

   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+
   |                                 | 2011-Aug-21 18:00                                                         | 2011-Aug-25 00:00                                                         | 2011-Aug-28 06:00                                                         |
   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+
   | **Quadrant**                    | NE               | SE               | SE               | NW               | NE               | SE               | SW               | NW               | NE               | SE               | SW               | NW               |
   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+
   | :math:`\boldsymbol{B}`          | 1.00             | 1.00             | 1.00             | 1.00             | 1.62             | 1.62             | 1.62             | 1.62             | 0.60             | 0.60             | 0.60             | 0.60             |
   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+
   | :math:`\boldsymbol{B_g}`        | 1.24             | 1.03             | 1.05             | 1.19             | 1.69             | 1.69             | 1.65             | 1.68             | 1.11             | 0.92             | 0.72             | 0.73             |
   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+
   | :math:`\boldsymbol{log_{10}R_o}`| 0.64             | 1.44             | 1.26             | 0.74             | 1.37             | 1.36             | 1.70             | 1.41             | 0.28             | 0.33             | 0.74             | 0.82             |
   +---------------------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+------------------+

It is evident that the radial wind profiles generated by the GAHM consistently
match the highest isotaches in all quadrants at different stages of Irene, no
matter how :math:`B` and :math:`\log_{10}R_o` vary. The AHM did a similarly good
job when the hurricane is strong (see middle panels), but failed to match the
highest isotaches when :math:`\log_{10}R_o < 1`. Both the AHM and the GAHM winds
died off too quickly away from the storm center, thus failed to match any lower
isotaches.

-  **Evaluation of the Maximum Winds and Radius to Maximum Winds**

Comparisons of the modeled maximum winds and radius to maximum winds to the
observed values in the input file were also carried out based on all 7 selected
hurricanes, given by the scatter plots in Figure 6. Evaluations of the maximum
winds are given in the upper panels, while the radius to maximum winds given in
lower panels, both color-coded by :math:`\log_{10}R_o`, with a simple linear
correlation given in each panel. Examination of the upper panels reveals that
the GAHM did an excellent job in estimating the maximum winds, with a few
overestimations near the lower bound of the dataset. Careful examinations of
these over estimated values revealed that they were from those "bad" dada
entries in the "best track" file that violate certain criteria in the GAHM when
solving for the :math:`R_{max}`. This phenomenon was particular common during
the dissipating stage of a hurricane. The AHM had larger discrepancies in
estimating the maximum wind compared to the GAHM, especially when
:math:`\log_{10}R_o < 1`, which was a direct consequence of the cyclostrophic
balance assumption made during the derivation of HM's equations. Examination of
the lower panels reveals that the maximum value of the modeled
azimuthally-varying :math:`R_{max}` failed to match the observed :math:`R_{max}`
values given in the input file, but the trend of the GAHM was significantly
better.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_8.1._GAHM.png
   :width: 800px


.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_8.2._GAHM.png
   :width: 800px


.. raw:: html

   <center>

Figure 6. Comparison of the modeled and “Best Track” maximum winds (upper two
panels), and the modeled and “Best Track” MATH 0 (lower two panels) between the
AHM and the GAHM based on all seven hurricanes.

.. raw:: html

   </center>

-  **Demonstration of the Multiple-Isotach Approach**

Earlier we have shown that a radial wind profile constructed by the GAHM using
the single-isotach approach would only match the highest isotach, due to
limitations of this single-fitting method. In fact, underestimations of modeled
winds at distances to isotachs other than the highest one were common, as the
radial wind profile tends to die off too quickly away from the storm center due
to the nature of GAHM’s formulas. In an effort to minimize the combined errors
mentioned above, and to improve the overall accuracy of the estimated wind
field, the multiple-isotach approach should be used whenever there is more than
one isotach present in the best track file.

The 3D plots of Irene’s radial wind profiles (left) and interpolated spatial
wind fields (right) by the GAHM using the single-isotach approach (upper panels)
versus the multiple-isotach approach (lower panels) were given by Figure 7. For
easier visualization, all available isotaches were plotted at radii to specified
isotaches in the left two panels, and as contour lines (after azimuthal
interpolation) in the right two panels. It is evident that winds generated by
the multiple-isotach approach were able to match all given isotaches in all 4
quadrants, while only the highest isotach was matched by the single-isotach
approach. Comparison of the spatial wind fields also indicated that the
multiple-isotach approach allowed the wind to die off more gradually away from
the storm center than the single-isotach approach did, demonstrated by the
smaller gradient of the contour lines in the lower panel. It is believed that
the multiple-isotach approach improves the overall accuracy and performance of
the GAHM.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_9._GAHM.png
   :width: 800px


.. raw:: html

   <center>

Figure 7. 3D plot of Irene’s radial wind profiles (left) and interpolated
spatial wind fields (right) by the single-isotach approach (upper panels) and
the multiple-isotach approach (lower panels).

.. raw:: html

   </center>

.. _evaluation_of_gahms_composite_wind_field:

Evaluation of GAHM's Composite Wind Field
-----------------------------------------

-  **Comparison of Radial Wind Profiles**

Using the multiple-isotach approach, radial wind profiles generated from the
GAHM were compared to those from the AHM, the SLOSH, H*Winds, and OWI winds,
shown in Figure 8 at three different stages of Irene (2011), with the left
panels showing the SW to NE cross-section winds, and the right showing the NW to
SE cross-section winds (same as in Figure 5). It is evident that the GAHM's
composite radial wind profiles matched all available storm isotachs at all time,
while the rest of the models failed most of the time. The SLOSH model, which is
also a parametric wind model, does not take the isotach information to construct
its wind field, so generally did a bad job matching any given isotaches. More
detailed wind structures can be observed from the radial profiles extracted from
the H*Wind and OWI winds than those from the parametric wind profiles, as
expected. However, profiles from different models in general did not match each
other due to different mechanisms involved in each model.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Picture12_GAHM.png
   :alt: Comparison of radial wind profiles of Irene (2011)


.. raw:: html

   <center>

Figure 8. Comparison of radial wind profiles of Irene (2011) at three different
stages among the AHM, the GAHM, SLOSH, H*Wind and OWI winds.

.. raw:: html

   </center>

-  **Evaluation of Model Results at Radii to Given Isotaches**

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Picture13_GAHM.png
   :alt: Comparison of specified isotachs and modeled winds


.. raw:: html

   <center>

Figure 9. Evaluation of modeled wind at radii to given isotaches based on all
seven selected storms.

.. raw:: html

   </center>

Quantitative evaluation of the GAHM's performance at radii to given isotaches
were given in Figure 9, with statistics given in Table 3, based on all seven
selected hurricanes. The GAHM had almost perfect match to each of the 34-, 50-,
and 64-kt isotaches, with a standard deviation around 0.1 kt, which was very
impressive. The AHM matched the 64-kt isotach reasonably well with a mean of
63.3 kts, but failed at the lower isotaches. The SLOSH did not take any isotach
information to construct its wind fields, thus behaved poorly in matching any
given isotaches. Its data also had the largest spread compared to other models,
with standard deviation ranging from 11 to 16 kts. For the H*Wind and the OWI
winds, the means of modeled winds were close to each specified isotach,
generally within ±3 kts, but the spreads of data were also large compared to
that of the AHM and the GAHM, with a standard deviation greater than 7 kts.

.. table:: Table 3. Statistical analysis of modeled winds at radii to
   specified isotaches based on all seven storms
   :class: wrap-table scroll-table
   :widths: 20 15 15 15 15 15 15

   +-----------+---------+--------+--------+--------+-----------------------+--------+
   |           | MEAN (kt)                          | Standard Deviation (kt)        |
   +           +---------+--------+--------+--------+-----------------------+--------+
   |           | Iso-34  | Iso-50 | Iso-64 | Iso-34 | Iso-50                | Iso-64 |
   +===========+=========+========+========+========+=======================+========+
   | **AHM**   | 26.9    | 44.7   | 63.3   | 7.96   | 7.09                  | 2.55   |
   +-----------+---------+--------+--------+--------+-----------------------+--------+
   | **GAHM**  | 34.0    | 50.0   | 64.0   | 0.10   | 0.12                  | 0.10   |
   +-----------+---------+--------+--------+--------+-----------------------+--------+
   | **SLOSH** | 30.0    | 48.9   | 69.5   | 11.00  | 14.11                 | 16.08  |
   +-----------+---------+--------+--------+--------+-----------------------+--------+
   | **OWI**   | 33.1    | 48.3   | 61.3   | 7.65   | 8.71                  | 9.42   |
   +-----------+---------+--------+--------+--------+-----------------------+--------+
   | **H*Wind**| 34.9    | 49.2   | 61.2   | 7.55   | 8.88                  | 8.89   |
   +-----------+---------+--------+--------+--------+-----------------------+--------+

-  **Comparison of Spatial Wind fields**

Comparisons of spatial wind fields among the AHM, the GAHM using the
multiple-isotach approach, the SLOSH, H*Wind and OWI winds (in rows) were given
by Figure 10 at the same snapshots (in columns) as in Figure 8. The 34-, 50-,
and 64-kt contour lines were shown in each wind snapshot. In general, the AHM
and GAHM winds shared a lot of similarities. During weaker periods of Irene, the
differences between the AHM and the GAHM's spatial wind fields were mostly
observed in the inner region of the hurricane due to differences in calculated
:math:`R_{max}`. More specifically, the calculated :math:`R_{max}` in the AHM
was under-predicted (thus is closer to the storm center) than that in the GAHM,
resulted from the faulty cyclostrophic balance assumption made in the HM and
AHM. During stronger periods of Irene, however, the differences were mostly
observed in the outer region of the hurricane due to the usage of the
multiple-isotach approach in the GAHM. Without taking into account storm
information from lower isotaches, the AHM winds tended to die off too quickly
away from the storm center. The SLOSH winds did not show much similarity to the
AHM and GAHM winds. It used an azimuthally constant :math:`R_{max}` from the
"best track" file to generate its vortex winds, and a distance-weighted
translational speed to account for storm asymmetry, which was not rendered
properly when the storm was strong. The winds fields of H*Wind and OWI winds
were much more complex than those of the simple parametric models, but the GAHM
did a relatively good job matching the spatial patterns of the hurricane.
Although it is unlikely that the parametric winds constructed over a minimum set
of data would match the level of details and complexity in the re-analysis
H*Wind and the numerical OWI winds, the ability of a the GAHM to produce
reasonable estimates of surface wind fields in a timely manner was highly
desirable for real-time forecasting.

.. figure:: /_static/images/model_configuration/meteorological_forcing/generalized_asymmetric_holland_model/Fig_10._GAHM.png
   :alt: Three snapshots (in columns) of Irene’s two-dimensional wind fields by the AHM, GAHM, SLOSH, H*Wind and OWI winds


.. raw:: html

   <center>

Figure 10. Comparison of Irene’s spatial wind fields by the AHM, GAHM, SLOSH,
H*Wind and OWI winds at three different stages

.. raw:: html

   </center>

References
----------

.. [1]
   Holland, G. J., 1980: An analytic model of the wind and pressure profiles in
   hurricanes. *Monthly Weather Review*, 108, 1212-1218.

.. [2]
   Mattocks, C., C. Forbes, and L. Ran, 2006: Design and implementation of a
   real-time storm surge and flood forecasting capability for the State of North
   Carolina. UNC-CEP Technical Report, University of North Carolina, 103pp.

.. [3]
   Powell, Mark & Uhlhorn, E. & Kepert, Jeffrey. (2009). Estimating Maximum
   Surface Winds from Hurricane Reconnaissance Measurements. *Weather and
   Forecasting*. 24. 10.1175/2008WAF2007087.1.

.. [4]
   The auxiliary preprocessing program ASWIP.F (located in the /wind directory)
   was further developed here to accommodate the GAHM.

.. [5]
   Queensland Government, 2011: Queensland Climate Change and Community
   Vulnerability to Tropical Cyclone, *Ocean Hazards Assessment – Stage 1*.




