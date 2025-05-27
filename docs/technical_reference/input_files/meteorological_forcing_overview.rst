.. _meteorological_forcing_overview:

Meteorological Forcing Overview
===============================

ADCIRC supports various meteorological forcing options through the :ref:`NWS <NWS>` parameter. This overview helps you understand the different meteorological forcing options and their characteristics.

Input Files for Meteorological Forcing
--------------------------------------

Depending on the NWS value selected, ADCIRC uses different input files:

* **fort.22**: Primary meteorological forcing file for most NWS options
* **fort.200**: Used for NWS=10 (AVN model) and NWS=11 (ETA model)
* **fort.221-224**: Used for NWS=12 (OWI format)

Meteorological Forcing Types
----------------------------

ADCIRC meteorological forcing can be categorized into several types:

1. **Direct Meteorological Input**
   - NWS=1,2,-2: Wind stress and pressure specified directly
   - NWS=5,-5: Wind velocity and pressure at all nodes
   - NWS=4,-4: Wind velocity and pressure at selected nodes

2. **Gridded Meteorological Data**
   - NWS=3: US Navy Fleet Numeric format
   - NWS=6: Rectangular grid of wind/pressure
   - NWS=7,-7: Regular grid of stress/pressure
   - NWS=10: NWS Aviation (AVN) model
   - NWS=11: NWS ETA 29km model
   - NWS=12: OWI format (nested grids)

3. **Parametric Hurricane Models**
   - NWS=8: Dynamic Holland model
   - NWS=19: Asymmetric hurricane vortex
   - NWS=20: Generalized Asymmetric Holland Model (GAHM)

4. **Data-Assimilated Hurricane Models**
   - NWS=15,-15: HWind files from NOAA HRD

5. **Coupled Model Systems**
   - NWS=100-199: Wave radiation stress with meteorological forcing
   - NWS=300-399: SWAN+ADCIRC coupled model

Format Comparison
-----------------

.. list-table:: Format Comparison
   :widths: 15 15 15 15 20 20
   :header-rows: 1
   :class: tight-table

   * - Format Type
     - Input Unit
     - Interpolation
     - Coverage
     - Advantages
     - Limitations
   * - Direct Input
     - Grid nodes
     - Temporal only
     - Full domain
     - Precise control
     - Large file sizes
   * - Gridded Data
     - Regular grid
     - Spatial & temporal
     - Must cover domain
     - Standard formats
     - Interpolation errors
   * - Parametric Models
     - Track data
     - Generated on-the-fly
     - Full domain
     - Small input files
     - Idealized storm structure
   * - Data-Assimilated
     - Snapshots
     - Temporal between snapshots
     - Storm area only
     - Real observations
     - Limited spatial coverage
   * - Coupled Models
     - Multiple sources
     - Model-dependent
     - Full domain
     - Physical consistency
     - Computational cost

Selecting the Appropriate Forcing
---------------------------------

When choosing a meteorological forcing option:

1. **Consider data availability**:
   - Historical simulations often use reanalysis or observation-based options
   - Forecasts typically use forecast model output or parametric models
   - Research applications may use idealized forcing

2. **Consider model domain**:
   - Large domains benefit from gridded data or parametric models
   - Small, high-resolution domains may benefit from direct input
   - Hurricane simulations typically use parametric or data-assimilated models

3. **Consider computational resources**:
   - Parametric models reduce I/O and storage requirements
   - Direct input may require large storage for forcing files
   - Coupled models incur additional computational cost

Common Challenges
-----------------

1. **Ensuring proper coverage**: Gridded meteorological data must completely cover the ADCIRC domain
2. **Time synchronization**: Pay careful attention to start times and time intervals
3. **Format consistency**: Follow exact format specifications for each NWS option
4. **Unit compatibility**: Ensure units are consistent with ADCIRC expectations
5. **Hot start considerations**: Some NWS options have specific requirements for hot starts

For detailed descriptions of each NWS option, see the :ref:`NWS <NWS>` parameter documentation.
For file format details, see :doc:`fort22` and :doc:`fort200` documentation. 