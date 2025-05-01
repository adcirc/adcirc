Fort.26: SWAN Runtime Parameters and Coupling File
==================================================

The fort.26 file contains the SWAN runtime parameters as well as the coupling details between ADCIRC and SWAN. This file is required to run the ADCIRC+SWAN model. Note that the name "fort.26" is ADCIRC's name (and fortran unit number) for the SWAN input command file, typically named INPUT by default.

File Structure
--------------

The format is detailed in the `SWAN User's Manual, Chapter 3 <http://swanmodel.sourceforge.net/download/zip/swanuse.pdf>`_.

See Casey Dietrich's great webpage on running `ADCIRC+SWAN <https://ccht.ccee.ncsu.edu/swanadcirc/>`_.

Specifying Input
----------------

Typically, SWAN+ADCIRC runs do not force wave initial conditions. However, the SWAN portion of the model may be able to accommodate SWAN 2D spectral input as a boundary condition, though this has not yet been tested/confirmed. The BOUNDSPEC command could allow the user to identify the fort.14 boundary nodestring for applying the spectra and then define the spectral parameters. The SWAN manual includes documentation on ASCII file format for SWAN spectra, which applies to both input/output spectra: [2].

Specifying Output
-----------------

Timing of SWAN global output files (swan_HS.63, swan_TPS.63, swan_TMM10.63, and swan_DIR.63) derives from the timing specified for ADCIRC global meteorological output in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`. Note that the SWAN calculation interval in the fort.26 must coincide with the ADCIRC meteorological output interval or SWAN output will not be written.

The fort.26 allows for output of wave parameters or 1D/2D spectra at user-specified points. For SWAN+ADCIRC applications, point output specification will typically occur after the QUANTITY command and before the TEST command. First, the user defines a group of points ('sname' from the SWAN manual) by either identifying their x and y coordinates directly in the fort.26 or within a separate file ('fname' in SWAN parlance). Next, the SPECOUT command triggers spectral output, and/or TABLE triggers time series output of desired wave parameters and time steps.

Example
-------

The following is an example from a fort.26 that requests both spectral output and wave parameter point time series:

.. code-block:: none

   POINTS 'SpecPts' File 'SpecLongLat.loc'

   SPECout 'SpecPts' SPEC1D ABS '1DSpecOutABS' OUT 19640903.000000 1200 SEC

   SPECout 'SpecPts' SPEC2D ABS '2DSpecOutABS' OUT 19640903.000000 1200 SEC

   TABLE 'SpecPts' HEAD 'WavesOut.txt' HS TPS DIR OUT 19640903.000000 1200 SEC

In this example:
- Longitude/latitude coordinates are specified (space delimited without header) inside a file called SpecLongLat.loc
- The two SPECOUT lines request 1D and 2D spectra written to 1DSpecOutABS and 2DSpecOutABS
- A table of Hs, Tps, and Dir will be written to WavesOut.txt every 1200 seconds beginning 9/3/1964 at 00:00:00
- The six values after the date and decimal represent HHMMSS, not decimal day
- The output files will be written to the PE subgrid directories that contain the output location points

Example contents of SpecLongLat.loc:

.. code-block:: none

   -81.299940 29.898530
   -81.303200 29.898450
   -81.303440 29.893390
   -81.300740 29.889340
   -81.295050 29.890450 