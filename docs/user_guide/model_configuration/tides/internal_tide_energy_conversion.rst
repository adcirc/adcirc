.. meta::
   :description: Internal Tide Energy Conversion in ADCIRC
   :keywords: adcirc, internal tide energy conversion

.. _internal_tide_energy_conversion:

Internal Tide Energy Conversion
===============================

Internal tide energy conversion refers to the energy conversion from barotropic
to baroclinic modes as surface tides flow over steep and rough topography in the
deep ocean generating internal tides. The "lost" barotropic tidal energy is
often accounted for through a linear friction term in large-scale numerical
models that are barotropic or not fine-scaled enough to resolve the energy
conversion. It is implemented in ADCIRC through a spatially varying nodal
attribute called
:ref:`internal_tide_friction <internal_tide_friction>`, in
the :ref:`fort.13 file <fort13>`.

Background and Theory
---------------------

The basic theory for generation of internal tides in the deep ocean was
established several decades ago. [1]_ However, it was not thought to be
incredibly important to the global energy balance of the surface tides until the
modern satellite era when it was discovered that internal tides are responsible
for approximately 30% of the global barotropic tidal dissipation. [2]_ [3]_ [4]_

Following this revelation, the past two decades have been subject to a number of
theoretical  [5]_  [6]_  [7]_  [8]_ and numerical  [9]_  [10]_  [11]_
investigations into internal tide generation and their effects on the surface
tides through parameterization of the energy conversion in large-scale numerical
tidal models.  [12]_  [13]_  [14]_  [15]_  [16]_  [17]_  [18]_

Attribute Summary
-----------------

In a computational domain covering a large portion of the deep ocean, the effect
of internal tide energy conversion may be needed to obtain more accurate tidal
solutions. The user should only elect to use the internal_tide_friction nodal
attribute when tides are included in the simulation through tidal boundary
conditions and tidal potential functions. The attribute may not be important for
domains that are small in size and/or do not cover a significant portion of the
"deep ocean" (the portion of the ocean excluding the continental shelf).

ADCIRC reads the internal_tide_friction attribute in as the *IT_Fric* variable,
which can have 1 (scalar) or 3 (tensor) dimensions. The attribute has dimensions
of [1/time], meaning that it is a linear friction term which is multiplied by
the velocity in the governing equations, and is normalized by the ocean depth
prior to simulation. Hence, it ignores the water surface elevation portion of
the total water depth, which is reasonable since the term and theory it is based
on is only applicable to deep ocean. Typically, it is only applied to ocean
depths greater than 100-500 m. [14]_ [16]_

Specifying *IT_Fric* Values
---------------------------

*IT_Fric* values are determined through analytical formulations based on Bell's
linear theory [1]_, valid in what is called sub-critical topography. [2]_

Recent publications using ADCIRC [14]_ [15]_ provide relevant formulation and
implementation details.

References
----------

.. raw:: html

   <references />

.. [1]
   T.H. Bell, Topographically generated internal waves in the open ocean, J.
   Geophys. Res. 80 (1975) 320–327. doi:10.1029/JC080i003p00320.

.. [2]
   C. Garrett, E. Kunze, Internal Tide Generation in the Deep Ocean, Annu. Rev.
   Fluid Mech. 39 (2007) 57–87. doi:10.1146/annurev.fluid.39.050905.110227.

.. [3]
   G.D. Egbert, R.D. Ray, Significant dissipation of tidal energy in the deep
   ocean inferred from satellite altimeter data, Nature. 405 (2000) 775–778.
   doi:10.1038/35015531

.. [4]
   G.D. Egbert, R.D. Ray, Estimates of M2 tidal energy dissipation from
   TOPEX/Poseidon altimeter data, J. Geophys. Res. Ocean. 106 (2001)
   22475–22502. doi:10.1029/2000JC000699.

.. [5]
   A. Melet, M. Nikurashin, C. Muller, S. Falahat, J. Nycander, P.G. Timko, B.K.
   Arbic, J.A. Goff, Internal tide generation by abyssal hills using analytical
   theory, J. Geophys. Res. Ocean. 118 (2013) 6303–6318.
   doi:10.1002/2013JC009212.

.. [6]
   F. Pétrélis, S.L. Smith, W.R. Young, F. Pétrélis, S.L. Smith, W.R. Young,
   Tidal Conversion at a Submarine Ridge, J. Phys. Oceanogr. 36 (2006)
   1053–1071. doi:10.1175/JPO2879.1.

.. [7]
   L. St. Laurent, C. Garrett, The Role of Internal Tides in Mixing the Deep
   Ocean, J. Phys. Oceanogr. 32 (2002) 2882–2899.
   doi:10.1175/1520-0485(2002)032\ \ <2882:TROITI>2.0.CO;2.

.. [8]
   L. St. Laurent, S. Stringer, C. Garrett, D. Perrault-Joncas, The generation
   of internal tides at abrupt topography, Deep Sea Res. Part I Oceanogr. Res.
   Pap. 50 (2003) 987–1003. doi:10.1016/S0967-0637(03)00096-7.

.. [9]
   J. Nycander, Generation of internal waves in the deep ocean by tides, J.
   Geophys. Res. 110 (2005) C10028. doi:10.1029/2004JC002487.

.. [10]
   A. Lefauve, C. Muller, A. Melet, A three-dimensional map of tidal dissipation
   over abyssal hills, J. Geophys. Res. C Ocean. 120 (2015) 4760–4777.
   doi:10.1002/2014JC010598.

.. [11]
   M.H. Alford, T. Peacock, J. a MacKinnon, J.D. Nash, M.C. Buijsman, L.R.
   Centuroni, S.-Y. Chao, M.-H. Chang, D.M. Farmer, O.B. Fringer, K.-H. Fu, P.C.
   Gallacher, H.C. Graber, K.R. Helfrich, S.M. Jachec, C.R. Jackson, J.M.
   Klymak, D.S. Ko, S. Jan, T.M.S. Johnston, S. Legg, I.-H. Lee, R.-C. Lien,
   M.J. Mercier, J.N. Moum, R. Musgrave, J.-H. Park, A.I. Pickering, R. Pinkel,
   L. Rainville, S.R. Ramp, D.L. Rudnick, S. Sarkar, A. Scotti, H.L. Simmons,
   L.C. St Laurent, S.K. Venayagamoorthy, Y.-H. Wang, J. Wang, Y.J. Yang, T.
   Paluszkiewicz, T.-Y.D. Tang, The formation and fate of internal waves in the
   South China Sea., Nature. 521 (2015) 65–9. doi:10.1038/nature14399.

.. [12]
   S.R. Jayne, L.C. St. Laurent, Parameterizing tidal dissipation over rough
   topography, Geophys. Res. Lett. 28 (2001) 811–814. doi:10.1029/2000GL012044.

.. [13]
   M.C. Buijsman, B.K. Arbic, J.A.M. Green, R.W. Helber, J.G. Richman, J.F.
   Shriver, P.G. Timko, A. Wallcraft, Optimizing internal wave drag in a forward
   barotropic model with semidiurnal tides, Ocean Model. 85 (2015) 42–55.
   doi:10.1016/j.ocemod.2014.11.003.

.. [14]
   W.J. Pringle, D. Wirasaet, A. Suhardjo, J. Meixner, J.J. Westerink, A.B.
   Kennedy, S. Nong, Finite-Element Barotropic Model for the Indian and Western
   Pacific Oceans: Tidal Model-Data Comparisons and Sensitivities, Ocean Model.
   129 (2018) 13–38. doi:10.1016/j.ocemod.2018.07.003.

.. [15]
   W.J. Pringle, D. Wirasaet, J.J. Westerink, Modifications to Internal Tide
   Conversion Parameterizations and Implementation into Barotropic Ocean Models,
   EarthArXiv. (2018) 9. doi:10.31223/osf.io/84w53.

.. [16]
   B.K. Arbic, A.J. Wallcraft, E.J. Metzger, Concurrent simulation of the
   eddying general circulation and tides in a global ocean model, Ocean Model.
   32 (2010) 175–187. doi:10.1016/j.ocemod.2010.01.007.

.. [17]
   A. Schmittner, G.D. Egbert, An improved parameterization of tidal mixing for
   ocean models, Geosci. Model Dev. 7 (2014) 211–224.
   doi:10.5194/gmd-7-211-2014.

.. [18]
   J.A.M. Green, J. Nycander, A Comparison of Tidal Conversion Parameterizations
   for Tidal Models, J. Phys. Oceanogr. 43 (2013) 104–119.
   doi:10.1175/JPO-D-12-023.1.
