.. meta::
   :description: Ramping in ADCIRC
   :keywords: adcirc, ramping

Ramping
=======

Ramping is a way by which terms can be steady increased over some period of time
in a simulation. This is most often done for model forcing terms like tides or
winds, in order to avoid applying a shock to the model.

Conceptual Justification
------------------------

To understand why ramping of forcing terms is needed, consider tides. Starting a
model with the water surface and velocity as zero everywhere, but with full
tidal forcing, is analogous to having a planet at rest, then instantly putting a
moon and sun in place, in motion with each other. This large, instantaneous
change in conditions (forcing) is what creates the shock, which tends to result
in spurious waves being formed. Ideally, one would initialize a tidal simulation
with the water surface and velocity at all points matching what it should be at
the point in time the model is initiated, given the phase of the tides, however
this is not generally achievable. Gradually scaling up forcing terms helps avoid
this problem.

.. _ramping_method:

Ramping Method
--------------

`NRAMP <NRAMP>`__ controls whether ramping is enabled, as well as how many
ramping terms there are. This is to permit users to apply different ramping
times for different forcing terms. All ramps in ADCIRC are applied as
(truncated) `hyperbolic
tangent <https://en.wikipedia.org/wiki/Hyperbolic_function>`__ functions over a
specified number of days. The various DRAMP\* variables control the number of
**days** over which the ramping is applied for individual forcing terms, they
are:

-  :ref:`DRAMP <DRAMP>` - any forcing terms not covered by the other DRAMP\* terms
-  :ref:`DRAMPExtFlux <DRAMPExtFlux>` - external flux boundary condition (BC)
-  :ref:`DRAMPIntFlux <DRAMPIntFlux>` - internal flux BC
-  :ref:`DRAMPElev <DRAMPElev>` - elevation specified BC
-  :ref:`DRAMPTip <DRAMPTip>` - tidal potential forcing
-  :ref:`DRAMPMete <DRAMPMete>` - meteorological forcing
-  :ref:`DRAMPWRad <DRAMPWRad>` - wave radiation stress forcing

Two other variables, :ref:`FluxSettlingTime <FluxSettlingTime>` and
:ref:`DUnRampMete <DUnRampMete>` (also in days), affect the timing of the ramping
terms.

-  If NRAMP=1, then DRAMP is relative to coldstart time.
-  If 1<NRAMP<8, then all DRAMP\* terms are relative to coldstart time plus
   FluxSettlingTime
-  If NRAMP=8, then all DRAMP\* terms are as above, except DRAMPMete is relative
   to coldstart plus FluxSettlingTime plus DUnRampMete

FluxSettlingTime is used because, for large river systems such as the
Mississippi River, it can take several days to initialize them and equilibrate.
It is a requirement when using the radiation flux boundary condition
(:ref:`IBTYPE <IBTYPE>`\ =52) because instabilities can quickly form otherwise.
DUnRampMete is useful because often a simulation is coldstarted without
meteorological forcing (such as a tide-only simulation), and meteorological
ramping is needed later when this forcing is added to a hotstarted simulation.

There are also namelist controls for some ramping terms, such as water level
offset forcing.

Typical Values
--------------

The following is a list of typical number of days for DRAMP\* variables for
different forcing terms.

-  Flux: 0.5-5 days
-  Tide: 5-30 days
-  Meteorology: 0.5-2 days

Example Usage Case
------------------

Consider the following illustrative example of how one might use all of these.
Let's say someone by the name of Arturo wants to run a hindcast simulation of
Hurricane Katrina impacting Louisiana and Mississippi with SWAN+ADCIRC, and
needs to get an accurate sense of the surge traveling up the Mississippi River.
Doing so requires specifying realistic flux boundary conditions at the upstream
end of the river. Feeling ambitious, Arturo sets NRAMP=1 and (because he has
just 5 days of met. forcing) sets DRAMP=1. His model goes unstable in the
rivers. Now feeling cautious, he first does a 30-day simulation to ramp in the
river and tide forcing, with NRAMP=3, DRAMP=15 days, DRAMPExtFlux=3 days, and
FluxSettlingTime=10 days. This means the coldstart simulation runs with only
river (external flux) forcing for the first 10 days, after which the tides are
ramped in over 15 days, and then the model is left at full-strength forcing for
another 5 days to equilibrate (10+15+5=30 days). A hotstart file is generated at
the end of the run. Arturo uses this to hotstart his 5-day
(:ref:`RNDAY <RNDAY>`\ =35) storm simulation, where DUnRAMPMete=30 ensures the
meteorological forcing ramping starts at the 30-day mark, and DRAMPMete=0.5
ensures a quick ramping of the meteorological forcing. Everything works great.
He's probably been overly cautious here, and might want to play around with
shortening some of the ramping times, but hey, at least it works.

Note that in the hotstart simulation, Arturo does not need to set any of the
other values on the DRAMP line non-zero because 30-days has already passed, and
so anything less than this is ignored by ADCIRC. Furthermore, note that
DRAMPWRad has been left at zero since ramping in the meteorological forcing
implies that the waves, and so the wave radiation stresses, should grow slowly.
