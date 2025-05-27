.. _fort25:

Fort.25, 255/277: Ice Coverage Input Files
==========================================

ADCIRC has the ability to incorporate ice fields into the computation of wind drag coefficients. The only influence that these ice fields has on ADCIRC’s solution procedure is thaaat they are used to modify the wind drag coefficients based on wind speed and percentage of ice coverage.
When using ice fields, the wind drag coefficient is first computed by the usual Garratt formula based on wind speed, then a wind drag coefficient is computed by a new function called IceCube that is based solely on the percent of ice coverage, then the maximum of the Garratt and IceCube drag coefficients is taken but must still be less than or equal to the drag cap in ADCIRC. If there is no ice in a particular area, then the drag coefficient is the same as Garratt.
In order to use ice field data, two changes are required to an existing fort.15: the ten thousands place must be set to 12 (see the NWS Simplified Table of Values for details), and the time increment of the ice field data (CICE_TIMINC) in seconds must be supplied on the WTIMINC line (see the Supplemental Meteorology/Wave/Ice Parameters table for details). Ice coverage datasets are linearly interpolated in time by ADCIRC, as a result, at least two datasets are required.
The ASCII text fort.25 file is only used to specify two general parameters for the Ice coverage data.

File Structure
--------------

The file format is as follows:

.. parsed-literal::
   
   numIceFields
   numBlankIceSnaps

If numIceFields in the fort.25 file is set to 1, only the fort.225 file will be read; if numIceFields in the fort.25 file is set to 2 then both the fort.225 and fort.227 files will be read. The fort.225 (the basin scale, i.e., large scale coarse grid) and optional fort.227 (the region scale, i.e., smaller scale finer grid) files contain the actual ice percentage values (0.0 to 100.0%) and a value of -1.0 which means land or no data available. A value of -1.0% is treated as 0.0% ice in ADCIRC.

These files are formatted exactly like the fort.221 barometric pressure files for the OWI win/pre meteorological data format (NWS=12). The data represent the percentage of area covered by ice within a grid cell. Grid dimensions and cell sizes are contained in the header lines of the fort.225 files.

If OWI-formatted meteorological forcing and ice coverage data are both used, the grid dimensions for ice coverage and meteorological data do not have to be the same; nor do the time increments
between data sets.

Method for Specifying the Coefficient of Drag
---------------------------------------------

The ADCIRC model ordinarily uses the wind drag coefficient formulation of Garratt (1977) in the calculation of surface wind stresses. It is a widely-used formulation and it has been found to work well for storm surge applications.
An additional physical process that has been examined in ice-covered regions such as coastal Alaska and the Great Lakes is the influence of sea ice as aerodynamic roughness elements. Macklin (1983) and Pease et. al. (1983) found that measurements of wind drag coefficients over first year sea ice typically yielded values that were significantly larger and varied less with wind speed than that predicted for open water.
More recent work (Birnbaum and Lupkes (2002) and Garbrecht et. al. (2002)) has formalized the effect of form drag on the specification of wind drag coefficients within marginal ice zones. From their work, Chapman et. al. (2005) utilized an empirical fit to the range of field data for the air-ice-water wind drag coefficient, CDF, and suggested:

.. math::

   \mathrm{CDF} = [0.125 + 0.5 \mathrm{IC} (1.0 – \mathrm{IC})] 10^{-2} \text{ (RaysIce formulation, this is not the default)}

in which IC is the ice concentration varying from 0.0 to 1.0 for open water and complete ice cover conditions, respectively. To specify this quadratic formulation for the wind drag, set the DragLawString parameter in the fort.15 file to “RaysIce”.
Inspection of the air-ice-water wind drag coefficient formula shows that a maximum value of 0.0025 occurs with 50% ice coverage. This value is very close to the Macklin (1983) measurement of 0.0028 for first year ice. Furthermore, it is seen that the value of the drag coefficient is symmetrical about 50% ice coverage suggesting that drag coefficient needed to represent 75% ice coverage is close to that of 25% ice coverage.
An alternative linear fit dependence on ice concentration has been applied by Danard et. al. (1989). These notions regarding variation of wind drag coefficient with ice cover have been supported by a number of Chukchi and Beaufort Sea storm surge simulations (Henry and Heaps, 1976; Kowalik, 1984; and Schafer, 1966) in which, wind drag coefficients greater than or equal to 0.0025 where utilized.
The method adopted for inclusion in ADCIRC considers the increased wind drag due to the presence of ice as developed by Chapman et. al. (2005). The method requires reading ice field concentration files into ADCIRC and calculating the wind drag coefficient values (variable over the model domain) using the form specified above. If ice cover is present, and the increased drag coefficient exceeds the value calculated using the standard Garratt (1997) formulation, the standard Garratt drag coefficient is replaced with the increased value associated with ice cover.
The IceCube formulation is used if the DragLawString is not specified in the fort.15, or is specified as “IceCube” or “default”. IceCube is a cubic function that represents a smoother transition between Garratt and the RaysIce formulation for low ice cover percentages and low wind speeds. This function is as follows:”

.. math::

   \mathrm{CDF} = [0.075 + 0.75 \mathrm{IC} – 0.9 \mathrm{IC}^2 + 0.2 \mathrm{IC}^3 ] 10^{-2} \text{ (IceCube formulation, this is the default)}

When IC = 0 the drag coefficient is equal to the minimum value for Garratt (.000075), when IC = 50 it has a value of 0.0025 and a zero gradient, and when IC = 100 it has a value of 0.00125.

References
----------

Birnbaum, G. and Lupkes, C. (2002). “A new parameterization of surface drag in the marginal sea ice zone,” Tellus, Vol. 54 A, pp. 107-123.

Banke, E. G. and Smith, S. D. (1973), “Wind stress on Arctic Sea ice, JGR, Vol. 78, pp. 7871- 7883. Chapman, R. S., Mark, D. and A. Cialone (2005). “Regional tide and storm-induced water level prediction study for the West Coast Alaska.” Draft Report to POA,

U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS. Chapman, R. S., Kim. S-K and D. Mark (2009). “Storm-Induced Water Level Prediction Study for the Western Coast of Alaska.” Draft Report to POA, U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS.

Chapman, R. S., Kim. S-K and D. Mark (2011). “Storm-Induced Water Level Prediction Study for Shatoolik, Alaska.” Draft Report to POA, U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS.

Danard, M. B., Rasmussen, M. C., Murty, T. S., Henry, R. F., Kowalik, Z., and Venkatesh, S.: 1989, Inclusion of ice cover in a storm surge model for the Beaufort Sea, Natural Hazards 2, 153–171.

Flather, R. A. (1988). “A numerical model investigation of tides and diurnal-period continental shelf waves along Vancouver Island,” Journal of Physical Oceanography 18, 115-139.

Garratt, J.R., (1977). “Review of drag coefficients over oceans and continents,” Monthly Weather Review 105, 915-929.

Garbrecht, T., Lupkes, C., Hartmann, J. and Wolff, M. (2002), “Atmospheric drag coefficients over sea ice—validation of a parameterization concept,” Tellus, Vol. 54 A, pp. 205-219.

Kowalik, Z. (1984). “Storm surges is the Beaufort and Chukchi Seas,” JGR, Vol. 89, No. C6, pp. 10,570-10578.

Schafer, P. J. (1966). “Computation of storm surge at Barrow, Alaska,” Archiv. Meteorol., Geophys. BioKimatol. Vol. A, No. 15(3-4), pp 372-393.