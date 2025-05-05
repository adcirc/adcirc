Ice Modifications
=================

ADCIRC has been modified to incorporate the effects of ice coverage on hydrodynamic processes. These modifications allow for more accurate simulations in polar and sub-polar regions where ice coverage significantly affects water circulation patterns.

Overview
--------

The ADCIRC ice modifications enable the model to account for:

1. Ice coverage effects on surface wind stress
2. Ice friction effects on water flow
3. Changes in momentum exchange due to ice-water interaction
4. Modifications to mass conservation equations to account for ice presence

Implementation
--------------

The implementation of ice modifications in ADCIRC involves several components:

* Modification of the governing equations to include ice effects
* Addition of ice-specific parameters in input files
* New output options for ice-related variables
* Special numerical treatments for partially and fully ice-covered regions

Usage
-----

To enable ice modifications in ADCIRC, specific parameters must be set in the fort.15 file. Additionally, ice coverage data may be provided through dedicated input files.

Detailed documentation on the implementation and usage of ice modifications can be found in the ADCIRC documentation: "Modifications to ADCIRC for ICE Coverage" (available at https://adcirc.org/wp-content/uploads/sites/2255/2016/01/adcirc_modifications_for_ice.pdf).

References
----------

For complete details on the theoretical foundation, implementation, and usage of ice modifications in ADCIRC, please refer to the official documentation at the ADCIRC website: https://adcirc.org/wp-content/uploads/sites/2255/2016/01/adcirc_modifications_for_ice.pdf 