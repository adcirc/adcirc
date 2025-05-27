.. _fort47:

Fort.47: Temperature Values at the Surface Layer
================================================

The fort.47 file records the top temperature boundary condition that is either provided via output from an atmospheric model or as an input variable into the code via a Surface Temperature Boundary Values File. This output file follows the same format as the :doc:`Water Surface Elevation <fort63>` file, as it only records the temperature values for the surface layer. The output file is only provided if the :ref:`IDEN <IDEN>` value in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>` is given as a 3 or 4, as ADCIRC evaluates the temperature changes for these two :ref:`IDEN <IDEN>` values.

Notes
-----

* Output format follows fort.63 conventions
* Only records temperature values at the surface layer
* Generated only when IDEN = 3 or 4
* Temperature boundary conditions can come from
    * Atmospheric model output
    * Surface Temperature Boundary Values File
* Output timing matches fort.63 parameters (NSPOOLGE) 