Command Line Options
====================

ADCIRC and its related tools accept various command line options that control their behavior. This page documents these options to help users automate and customize their ADCIRC runs.

ADCIRC Command Line Options
---------------------------

ADCIRC has several command line options that control file locations, I/O behavior, and more. These options can be combined as needed.

``-I INPUTDIR``
   Look for input files in the directory 'INPUTDIR' instead of looking in the directory where the ADCIRC command was issued.

``-O GLOBALDIR``
   Write fulldomain output files to the directory 'GLOBALDIR' instead of writing them to the directory where the ADCIRC command was issued.

``-L``
   When running in parallel, write subdomain output and binaryhotstart files to the PE* subdirectories instead of writing fulldomain output and binary hotstart files. This option does not affect the writing of NetCDF hotstart files, which are always written for the full domain.

``-S``
   When running in parallel, write subdomain binary hotstart files to the PE* subdirectories instead of writing fulldomain binary hotstart files. This option does not affect the writing of NetCDF hotstart file, which are always written for the full domain.

``-R``
   When running in parallel, read subdomain binary hotstart files instead of a fulldomain binary hotstart file. This option does not affect the reading of NetCDF hotstart files, which always reflect the fulldomain solution.

``-W NUM_WRITERS``
   When running in parallel, dedicate NUM_WRITERS MPI processes solely to the task of writing ascii output files (either sparse or dense ascii formats). Must be a positive integer. This option only affects the writing of the ascii formatted fort.63, fort.64, fort.73, and fort.74 files.

``-Ws NUM_SPLIT_WRITERS``
   Same as -W option above, but a new file is created each time the data are written; the integer part of the simulation time is appended to the file name.

``-H NUM_HOTSTART_WRITERS``
   Same as -W option above, except that NUM_HOTSTART_WRITERS MPI processes will be dedicated to writing binary hotstart files. This option does not affect the writing of NetCDF hotstart files. The binary hotstart files that are written using this option are only suitable for 2DDI runs without harmonic analysis.

Version 51 Additions
~~~~~~~~~~~~~~~~~~~~

``-M``
   When running in parallel, write subdomain harmonic analysis files. This option is useful for tidal analyses with many constituents (e.g., more than approximately 30) because writing fulldomain harmonic analysis output files with this many constituents can exhaust the physical memory of the host computer.

``-Mft MESH_FILE_TYPE``
   The MESH_FILE_TYPE can be specified as ascii (the default) or xdmf (introduced in ADCIRC v51).

``-Mfn MESH_FILE_NAME``
   Search for a full domain mesh file named MESH_FILE_NAME rather than the default mesh file name fort.14.

``-Nft NODAL_ATTRIBUTES_FILE_TYPE``
   The NODAL_ATTRIBUTES_FILE_TYPE can be specified as ascii (the default) or xdmf (introduced in ADCIRC v51).

``-Nfn NODAL_ATTRIBUTES_FILE_NAME``
   Search for a full domain nodal attributes file named NODAL_ATTRIBUTES_FILE_NAME rather than the default nodal attributes file name fort.13.

Examples
~~~~~~~~

Run ADCIRC with input files in a different directory::

   ./adcirc -I /path/to/input/files

Run parallel ADCIRC with input in one directory and output in another::

   mpirun -np 16 ./padcirc -I /path/to/input/files -O /path/to/output/files

Run parallel ADCIRC with subdomain output::

   mpirun -np 16 ./padcirc -L

Run parallel ADCIRC with dedicated writers::

   mpirun -np 16 ./padcirc -W 2

Run parallel ADCIRC with a custom mesh file name::

   mpirun -np 16 ./padcirc -Mfn my_custom_mesh.grd

ADCPREP Command Line Options
----------------------------

ADCPREP is used to partition the ADCIRC domain for parallel execution. If ADCPREP is run without command line options, it presents an interactive menu. Using command line options allows for automation.

``--np NUM_SUBDOMAINS``
   Decompose the domain into NUM_SUBDOMAINS subdomains. The number of subdomains should be equal to the number of compute processors that will be used in the parallel computation. This command line option is required for all other ADCPREP operations.

``--partmesh``
   Partition the mesh only; that is, decide which subdomain each of the nodes should fall into. This should be done first. The result is a file called partmesh.txt, which consists of a list of integers, one line per node in the fulldomain mesh. Each integer in the list represents the subdomain number that each node will fall into. As a result, the range of values in the file represents the range of subdomain numbers.

``--prepall``
   Use the partmesh.txt file to decompose all the ADCIRC input files into subdirectories numbered 'PExxxx' where 'xxxx' is the zero indexed, zero padded subdomain number. Expects all input files to have the default names, e.g., the mesh file must be named fort.14 instead of something like myMeshFile.grd in order for this option to work. This option requires that adcprep has already been executed with the partmesh option, and that the partmesh.txt file was written successfully.

``--prep13``
   Only redecompose the nodal attributes (fort.13) input file. Requires that prepall step has already been performed.

``--prep15``
   Only redecompose the control model parameter and periodic boundary conditions file (fort.15) input file. Requires that prepall step has already been performed.

``--prep20``
   Only redecompose the non-periodic, normal flux boundary condition file (fort.20). Requires that prepall step has already been performed.

``--prep88``
   Only redecompose the upland river initialization file (fort.88). Requires that prepall step has already been performed.

Examples
~~~~~~~~

Complete parallel preparation for 16 processors::

   # Step 1: Partition the mesh
   adcprep --np 16 --partmesh
   
   # Step 2: Prepare all input files
   adcprep --np 16 --prepall

Re-decompose just the nodal attributes file after making changes::

   adcprep --np 16 --prep13

Notes on Parallel Execution
---------------------------

When running ADCIRC in parallel, it's important to understand the workflow:

1. First run ADCPREP to partition the domain
2. Then run padcirc with the appropriate options
3. Post-processing may be needed to combine subdomain outputs

Performance considerations:

- The ``-W``, ``-Ws``, and ``-H`` options can improve performance by dedicating specific processes to I/O operations
- Using the ``-L`` option reduces memory requirements but requires post-processing to combine outputs
- For large domains with many tidal constituents, the ``-M`` option can prevent memory exhaustion 