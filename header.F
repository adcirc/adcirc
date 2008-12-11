C******************************************************************************
C           This file contains a log of changes made to (P)ADCIRC.            *
C                                                                             *
C  mod history                                                                *
C  v48.41     - 12/11/08 - rjw - fix massbal.F COMPUTE_RESELEM subroutine
C               In order to compare the global to local locations 
C               Global x and y need to be in converted units not in Lat and Lon
C               use subroutine CPP to make conversion
C  v48.40     - 12/11/08 - rjw - fix version number in version.F
C  v48.39     - 12/11/08 - rjw - applied bug fixes in cstart.F and hstart.F 
C               format statement, changed default il value to 2 in vsmy.F
C  v48.38     - 12/10/08 - jgf - Added #ifdef HAVE_MPI_MOD around mpif.h 
C               in massbal.F to fix failed compile on Cray-XT3 (sapphire).
C  v48.37     - 12/10/08 - mcf - modified netcdfio metadata
C  v48.36     - 12/03/08 - jgf - Fixed a bug in PREP22 of adcprep,
C               combined CASE(3,6) into CASE(8,9) b/c the approach of
C               using the shell to copy files in CASE(3,6) had been
C               commented out.
C  v48.35     - 12/03/08 - jgf - Fixed a bug in adcprep that was causing 
C               unprintable characters to appear in the names of
C               subdomain fort.22 files on the IBM p5 platform.
C  v48.34     - 12/02/08 - jgf - Fixed a bug in adcprep to look for
C               tau0min and tau0max values in fort.15 if tau0 is less
C               than -5.0 AND greater than -6.0. 
C  v48.33     - 11/25/08 - jgf - Added a check for negative central pressure
C               deficit in Symmetric Holland code (NWS=8). Also fixed a bug
C               introduced in the merge of Seizo's fix to internal barrier
C               boundary handling in v48.32. Reimplemented tau0base
C               formulation  when TAU0=-3.0, moved fulldomain time
C               varying TAU0 option to TAU0=-5.0 ... if tau0 output is
C               desired, the corresponding options are -3.1 or -5.1.
C  v48.32     - 09/16/08 - jgf - Merged fix from C.Mattocks to include 
C               directives in metis, e.g., <metis.h> to "metis.h".
C               Fixed HollandGet so that lat is read in as four digits, thus 
C               allowing storm locations west of 99.9W lon; integrated Seizo 
C               fixes to decomp.F (mpi communication table) and in timestep.F
C               (wet/dry checks along internal barriers). 
C  v48.31     - 09/16/08 - mcf - added hotstart file format identifiers in 
C               from version.F, modified netCDF metadata, added time counters 
C               and record numbers for hotstart files, fixed netCDF hotstart 
C               flow for preprocessor, adcirc and padcirc.
C  v48.30     - 08/29/08 - rjw - merged CD's mods for 3D wetting and Drying
C  v48.29     - 08/07/08 - mcf - modified starting character counter in parser 
C  v48.28     - 08/01/08 - mcf - included File Format Version parameters
C              in netCDF hotstart files and changed FileFmtMinor to 2 to account
C               for the removal of an extra empty record in hotstart files
C  v48.27     - 07/30/08 - mcf - corrected dt when writing initializing netCDF
C               hotstart files in prep.F
C  v48.26     - 07/30/08 - mcf - corrected noff range value when writing netCDF
C               hotstart files (range 0 to 1)
C  v48.25     - 07/30/08 - mcf - corrected memory/compiling/running issues
C               with arrays in the binary and netCDF hotstart capability
C  v48.24     - 07/25/08 - mcf - Added code in the hstart.F file to check if 
C               hotstart binary file exists in subdomains for netcdf hotstart 
C               input file option when running the parallel version.
C  v48.23     - 07/25/07 - jgf - Fixed hotstart file incrementation in 
C               write_output.F, prep.F, hstart.F and harm.F. Updated hot2asc.F
C               utility with fixed file incrementation. Restored machdep.F.
C               Updated cmplrflags.mk for apple darwin (removed -DIBM).
C  v48.22     - 07/24/07 - jgf - Fixed hotstart file incrementation in 
C               write_output.F and hstart.F.
C  v48.21     - 07/24/08 - mcf - Added code in the hstart.F file to only read
C               netCDF hotstart files when running the serial version.
C  v48.20     - 07/24/08 - mcf - Added code in the preprocessor to initialize
C               netCDF hotstart files.
C  v48.19     - 07/22/08 - mcf - Added code to assign the globalized
C               hotstart variables to netCDF variables only in processor 0.
C  v48.18     - 07/09/08 - mcf - Added code to adcprep and adcirc to 
C               initialize and write netCDF hotstart files.
C  v48.17     - 06/09/08 - jgf - Added code to adcprep to decompose UnSWAN
C               boundaries and distribute the UnSWAN input file to each of the
C               subdomain directories.
C  v48.16     - 06/05/08 - jgf - Merged additional changes from cf's latest
C               (from adc48_15_4).
C  v48.11     - 03/14/08 - jgf - Added cf's latest write_output.F, corrected
C               Chezy conversion to quadratic in nodalattr.F, added 
C               info messages to cmplrflags.mk and makefile, corrected issues
C               to allow use of g95 compiler. 
C  v48.10     - 03/13/08 - jgf - Added cf's latest prep.F, read_global.F,
C               global.F, and netcdfio.F to enable boundaries and station
C               names in the metadata and to fix bugs. 
C  v48.09     - 02/25/08 - jgf - Added cf's latest netcdfio.F and 
C               write_output.F and synchronized arguments. Also added
C               ibtypenc variable (and allocation). Fixed bug in makefile
C               that prevented compilation of adcprep without netcdf.
C  v48.08     - 02/19/08 - jgf - Removed trailing space from "pgi " in 
C               cmplrflags.mk. Removed p15 and other utilities from 'all' 
C               target in makefile. Corrected bug in FEGRIDS sub in vsmy.F 
C               where variables were declared as REAL rather than REAL(SZ). 
C               Also made all constants into doubles in FEGRIDS sub.
C               In lsq.F, changed rinv to intent(inout) on line 538 to relieve
C               compiler issue on certain platforms. Added PREPNETCDF
C               subroutine to adcprep to initialize netcdf output files for
C               parallel runs. Updated makefile to reflect the linking of 
C               the netcdf module with adcprep.
C  v48.07     - 02/07/08 - jgf - Merged changes from mcf for netCDF.
C  v48.04     - 11/19/07 - jgf - Fixed bug that was preventing NetCDF metadata
C               from being written to subdomains. Added ABS() to station
C               mappings when reading the fort.80 file in adcpost. Added
C               inverted barometer boundary condition (commented out pending
C               successful testing). 
C  v48.03     - 10/19/07 - jgf - Created write_output.F, rearranged NetCDF
C               to use globalio. Merged kmd's 3D hotstart decomposition
C               in prep.F.
C  v48.02     - 10/04/07 - jgf - Fixed bug in UPDATEC3D in messenger.F.
C  v48.01     - 10/03/07 - jgf - Changed CrayXT3 and PGI compiler flags
C               (-c89 and -Mextend). Fixed bug that caused the min and max
C               variable tau0 to be read even if not needed. Fixed bug in
C               post52 that was writing NSTAE instead of NSTAV. Fixed missing
C               comma in FORMAT statement in globalio.F. Fixed missing right
C               paren in timestep.F.
C  v48.00     - 10/02/07 - jgf - Added ONLY keyword to USE NETCDF in 
C               timestep.F to enable it to compile without conflicts.
C  v47.13     - 09/28/07 - jgf - Fixed two typo bugs in post.F for binary
C               fort.64 and fort.74 files. Fixed 3D stations bug that gave
C               incorrect output for stations inside ghost elements. Fixed
C               mixup between NOUTGE=3 and NOUTGE=4.
C  v47.12       09/05/07 - jgf - Merged code for netcdf from cf. 
C  v47.11     - 09/03/07 - jgf - Added code to read in limits of time varying
C               tau0 directly from fort.15 file.
C  v47.10     - 08/28/07 - jgf - Merged code back in to create fort.80 file 
C               for use in post processing 3D and harmonic analysis files. 
C               Merged adcpost back in.
C  v47.09     - 07/26/07 - jgf - Fixed two bugs with new nodal attribute that
C               were causing SIGSEGV.
C  v47.08     - 07/24/07 - jgf - Added new nodal attribute for floor and 
C               ceiling of variable tau0. Made "fully populated" output format 
C               standard again for fort.63, fort.64, etc. Made "sparse" output
C               format selectable by setting NOUTGE(etc)=3.
C  v47.07     - 07/23/07 - jgf - Added feature to inflate to create chunked
C               output. Removed the hardwired 2nd floor for variable tau0.
C  v47.06     - 07/16/07 - jgf - Fixed bug in timestep.F where NSTAE was
C               used to write 2D velocity station file. Added global output
C               for time varying tau0.
C  v47.05     - 07/15/07 - jgf - Enhanced inflate.F so that it count the
C               number of datasets in the file and place the correct number
C               at the top in NDSETS. Updated 3D header and store subroutines
C               in globalio. 
C  v47.04     - 07/05/07 - jgf - Merged latest asymmetric vortex wind model,
C               made CF from fort.15 into floor of equivalent quadratic
C               bottom drag when Manning's N is used, added hstime utility
C               and updated it for the latest hotstart file format.
C  v47.03     - 06/30/07 - jgf - Changed inflate to use infinite loop instead
C               of relying on NDSETS.
C  v47.02     - 06/28/07 - jgf - Rewrote the prep13 subroutine so that it
C               only uses one open file at a time to avoid exceeding the
C               Fortran compiler limitation on simultaneously open files;
C               changed to unformatted read of nodal attribute files to
C               make them more tolerant of leading spaces and end-of-line
C               characters. 
C  v47.01     - 06/21/07 - jgf - Applied the fixes from the stable version
C               of ADCIRC to v46.57. 
C  v47.00     - 06/04/07 - jgf - Applied the changes from the stable ADCIRC
C               series v46.32jgfxx to version 46.52 as follows:
C  v46.32jgf12- 02/20/07 - jgf - Fixed bug in adcpost where format string was
C               one field short.
C  v46.32jgf11- 02/19/07 - jgf - Fixed bug in 3D station post processing
C               in adcpost. Added hsconv program to convert v44.15 hot start
C               files to v46.32 hot start files for operational forecasting
C               system.
C  v46.32jgf10- 01/10/07 - jgf - Fixed bug by renaming LV_TEMP_A to TEMP_LV_A.
C               Fixed bug in Holland subroutine by using a backward difference
C               instead of forward.
C  v46.32jgf8 - 11/07/06 - jgf - Commented out section in adcprep that
C               reverses the sign of NOUTE, etc for hotstarts. Fixed bug in
C               code that found the end of output files for appending during
C               hotstart; now uses infinite loop rather than trusting the
C               number of datasets listed at the top of the file. Fixed bug
C               in logical if statement in nodal attributes by adding
C               parentheses because different compilers were treating
C               the logical statements differently. Added Shintaro's wet/dry
C               changes, including his check for infinitely narrow
C               bottlenecks. Modified compiler flags for SGI. Fixed bug in
C               makefiles to make all .Fs capital.
C  v46.32jgf6 - 10/19/06 - jgf - Updated PC makefiles.

C ----------------------------------------------------------------------------
C
C  v46.57     - 02/15/07 - sb - The change on RS2GET made at v46.54 was       *
C               removed.  The original strategy (explained in the comment on  *
C               v46.54) is now used.                                          *
C  v46.56     - 02/09/07 - sb -Dedicated writer feature was implemented.      *
C  v46.55     - xx/xx/xx - This version number was skipped.                   *
C  v46.54     - 01/03/07 - sb - RS2GET in rs2.F now takes maximum values if   *
C               stwave domains overlap. (It used to take the value of the     *
C               domain of the first appearance in a fort.23 file.)            *
C               This requires more memory,                                    *
C             - SB_WETDRY compilation flag to enable w/d fixes that had been  *
C               implemented since v46.48 was intoduced.                       * 
C               (NOT enabled by default.)                                     *
C             - Made torelant to leading white spaces in attribute names      *
C               in fort.13.                                                   *
C  v46.53     - 12/11/06 - sb - Fixed a minor bug.                            *
C  v46.52     - 11/17/06 - sb - Fixed a minor bug.                            *
C  v46.51     - 11/17/06 - sb - Compilation flags to disable some of max files*
C                have changed:                                                *
C               -DNELMAX (to suppress maxele.63)                              *
C               -DNUVMAX (to suppress maxvel.63)                              *
C               -DNPRMAX (to suppress minpr.63)                               *
C               -DNWVMAX (to suppress maxwvel.63)                             *
C               -DNRSMAX (to suppress maxrs.63)                               *
C  v46.50     - 11/03/06 - sb - Two of the three w/d modifications introduced *
C               in v46.48 are now commented out.  Only one remains. Maximum   *
C               values of several variables in time at each node are now      *
C               written out in fort.63 format at the last time step. The new  *
C               output files are as follows:                                  *
C               maxele.63, maxvel.63, maxwvel.63, minpr.63, maxrs.63          *
C               These outputs can be suppressed by adding some of the         *
C               following compiler options:                                   *
C               -DNEVMAX (to suppress maxele.63 and maxvel.63)                *
C               -DNPWMAX (to suppress minpr.63 and maxwvel.63)                *
C               -DNRSMAX (to suppress maxrs.63)                               *
C  v46.49     - 11/03/06 - sb - One of three wetting/drying algorithm         *
C               modifications introduced in v46.48 was commented out.         *
C               Sparse global output file writing was implemented.            *
C               fort.73 file output was suppressed.                           *
C  v46.48     - 11/03/06 - sb - Modified wetting and drying algorithm         *
C             - Fixed a bug in river ramping                                  *
C             - Fort.63, 64, 73 and 74 files are now in sparse format.        *
C               Use "inflate" program under "util" directory to inflate the   *
C               outputs files.                                                *
C  v46.41     - 10/04/06 - rtm removed system calls, now using                *
C               c routines.                                                   *
C  v46.40     - 10/04/06 - sb added dynamic tau0 feature.  This feature kick  *
C               in when TAU0 = -3.0 in fort.15.                               *
C  v46.39     - 10/01/06 - vjp, rtm - added global output files for units     *
C               61,62,63,64,67-68,71,72,73,74, but 51,52,53,54,55 not yet.    *
C  v46.33     - 09/21/06 - vjp, rtm, sb, jgf - Merged the following changes   *
C               from:                                                         *
C               rtm, vjp: added option 'partmesh' to preprocessor to partition*
C               mesh only; added option 'prep13' to decompose the fort.13     *
C               file only; added subroutine relocalize to relocalize selected *
C               files after a prepall; added GLOBAL_IO module to write global *
C               output files in parallel on-the-fly; improvements to          *
C               build13.F. Merged the following changes from sb: added        *
C               OWIWIND module with NWS=12 for direct input of OWI wind data  *
C               on rectangular grid (no preprocessing); added RS2 module with *
C               NWS=2xx for STWAVE output direct read; added buildstwave23.F. *
C  v46.32     - 09/15/06 - jgf - Made boundary layer adjustment in Holland    *
C               model a parameter in the fort.15 file.                        *
C  v46.31     - 09/15/06 - jgf - Fixed integer divide in Manning's N          *
C               exponent bug. Fixed bug in first read of Holland fort.22.     *
C               Now subtracting translation speed before calculating Holland  *
C               parameters.                                                   *
C  v46.30     - 09/07/06 - jgf - Fixed bug for reading new Holland reference  *
C               time and allocating memory in Holland subroutine.             *
C  v46.29     - 09/07/06 - jgf - Fixed compilation bugs in lsq.F and vortex.F.*
C               Added build12.F utility and new version of build13.F utility  *
C               from vjp. Added reference time to Holland model, subtracted   *
C               storm trans. vel. from wind speed before calcs in Holland.    *
C  v46.28     - 08/26/06 - jgf - Fixed IDENMSG bug (kd) in prep.F, updated    *
C               preprocessing code to support baroclinic and transport runs,  *
C               added PREP10 subroutine. Added call to timeconv in Holland    *
C               subroutine to calculate time between hindcasts. Added         *
C               StormNumber to replace WTIMINC for NWS=8 (Holland model) for  *
C               wind multiplier support.                                      *
C  v46.27     - 08/22/06 - jgf - Replaced IRType with IDen in header for      *
C               fort.44 files. Added scheme to handle negative default values *
C               of tau0 in nodal attributes file.                             *
C  v46.26     - 08/22/06 - jgf - Moved new wind.F into src directory where it *
C               belongs.                                                      *
C  v46.25     - 08/21/06 - jgf - Implemented choice of tau0 schemes for       *
C               negative tau0 in fort.15 file, updated asymmetric hurricane   *
C               code from cf.                                                 *
C  v46.24     - 08/16/06 - jgf - Added updated p15.F and owi22.F from sb.     *
C               widened output format in fort.63 to avoid truncation bug when *
C               elevations have negative three digit exponents. Lifted        *
C               restriction on overlap of IBTYPEs 52 and 24. Fixed bug in     *
C               prep22 for NWS=1,3,4 where outer loop was missing and only    *
C               the first set of wind data was processed.                     *
C  v46.23     - 08/13/06 - jgf - Fixed a bug in the multiple ramp functions   *
C               in timestep.F for the Katrina ramp.                           *
C  v46.22     - 08/09/06 - jgf - Added a 1/0.7 factor to max wind forecast in *
C               Holland wind model.                                           *
C  v46.21     - 08/09/06 - jgf - Added code to support new IBTYPE=52. Changes *
C               to Holland wind model: added debugging code, fixed bug that   *
C               was miscalculating wind time increment, added 0.88 factor to  *
C               convert from 1 minute to 10 minute winds. Fixed bug in        *
C               prep13 that was breaking long lines of default values onto    *
C               multiple lines.                                               *
C  v46.20     - 08/03/06 - jgf - Applied a 0.7 factor to Holland winds, fixed *
C               a bug in header for debug elevations.                         *
C  v46.19     - 07/26/06 - jgf - Discovered and fixed a bare #ifdef bug on    *
C               periodic bcs in timestep.F. Added negative NSCREEN option to  *
C               write screen output to a fort.999 file.                       *
C  v46.18     - 07/21/06 - jgf - Implemented horizontal eddy viscosity        *
C               and horizontal eddy diffusivity as nodal attributes for       *
C               research purposes, not to be used in production runs.         *
C  v46.17     - 07/21/06 - jgf - Fixed empty arguments bug in WriteWarnElev() *
C               subroutine call.                                              *
C  v46.16     - 07/20/06 - jgf - Added #ifdef IBM to INT() calls to           *
C               work around bugs in non-IBM compilers that don't like the     *
C               KIND argument. Fixed a col. 73 bug in in prep.F. Changed      *
C               METIS1 to METIS in makefile to fix CrayXT3 build bug.         *
C               Merged asymmetric hurricane wind model from cf and cm.        *
C  v46.15     - 07/10/06 - jgf - Added a FORMAT statement to prep.F for       *
C               fort.13 to keep attributes with many values on a single line. *
C  v46.14     - 07/07/06 - jgf - Fixed bug from Katrina testing in hstart.F   *
C               where too many records were read from hotstart file.          *
C  v46.13     - 07/06/06 - jgf - Fixed build bug in turb in vsmy.F where real *
C               exponentiation was used on a PARAMETER.                       *
C  v46.12     - 07/06/06 - jgf - Added debugging code for wet/dry. Removed    *
C               the check on KEMAX from screen output logic. Merged fixes     *
C               from vjp: modified two lines in src/wind.F for IBM, added a   *
C               util subdir and put build13.F in it, added wind subdir and put*
C               p15.F and owi22.F in it, modified the makefile (but did not   *
C               modify the pc oriented makefiles) to build new targets and    *
C               reduced fflags to just 3 instead of 5 types, modified the     *
C               cmplrflags.mk file. Fixed a bug in subroutine turb in vsmy.F  *
C               where constants were initialized but not saved.               *
C  v46.11     - 07/05/06 - jgf - Added debugging code for high elevations.    *
C  v46.10     - 07/03/06 - jgf - Fixed bug (ral) that prevented friction      *
C               flags from being set for spatially varying friction.          *
C  v46.09     - 06/20/06 - jgf - Fixed bug in timestep.F when NRamp=0. Fixed  *
C               handling of NOLIFA=3 in read_input.F. Fixed bugs in Holland   *
C               subroutines of wind.F (mult. lon by -1, mult RMW by 1000).    *
C  v46.08     - 06/14/06 - jgf - adcpost: Fixed bug in number of stations     *
C               reported in fort.72 files. Upgraded file handling in Holland  *
C               subroutine. Created fine-grained ramp functions by converting *
C               NRAMP in fort.15 to the number of ramp functions.             *
C  v46.07     - 06/09/06 - jgf - Fixed bug in reading Holland files.          *
C  v46.06     - 06/08/06 - jgf - Fixed bug in reading WTIMINC in adcirc       *
C               and increased line size in adcprep for Holland wind model.    *
C               Fixed compile bug on IBM p5.                                  *
C  v46.05     - 06/07/06 - jgf - Added code to calculate Cd from Mannings N   *
C               at every time step. Rewrote Holland subroutine to accept      *
C               hybrid hindcast/forecast input files. Replaced metis.F with   *
C               a new one from vjp with upgrades to weir node dual handling.  *
C               Fixed bug that prevented linking with libmetis.a on IBM p5.   *
C  v46.04     - 06/05/06 - jgf - Applied bde's fixes to compile on IBM p5.    *
C  v46.03     - 05/30/06 - jgf - Implemented Holland Wind model (NWS=8)       *
C               subroutine and added code to adcprep to copy NWS8 files       *
C               to subdirectories. Commented out some features req'd for      *
C               Katrina compliance.                                           *
C  v46.02     - 05/25/06 - jgf - Implemented GeoidOffset code to allow a      *
C               permanent water level offset from the geoid. Implemented      *
C               wind drag cap function. Implemented directional wind          *
C               roughness and canopy factor. Added a combined radiation /     *
C               flux condition for river boundaries (IBTYPE=22). Added a      *
C               hard-coded spatially varying tau0 scheme vbased on depth to   *
C               match v42. Applied kd's BPG bug fix. Fixed major bugs in      *
C               ADCPREP/ADCPOST to post-process the wind recording station    *
C               files, changing the format of the fort.80 file.               *
C  v46.00     - 04/18/06 - jgf - Adcirc: Created new nodal attributes module  *
C               to hold node-by-node data that is time invariant. Applied     *
C               vjp's fixes to itpackv.F and messenger.F, added kind parameter*
C               to int() calls, made nscreen the time step increment between  *
C               screen outputs (except for error messages) to work around     *
C               limited log file sizes on certain platforms, added option for *
C               NWS=7 (grid-independent wind stress).  Adcprep: Added code to *
C               break up new nodal attributes file (unit 13). Removed code    *
C               that decomposed the StartDry and spatially varying bottom     *
C               friction files. Added implicit none to all subroutines, added *
C               trim() function to cut down files size from adcpost, made     *
C               infinite loop over data sets to fix bug where the wrong number*
C               of data sets are listed at the top of subdomain files after   *
C               hot start. Added capability to adcprep option 2 to skip       *
C               processing of certain files. Build: added vjp's compiler      *
C               flags, meb's small fixes to makefiles, updates to PC          *
C               makefiles, new version of config.guess.                       *
C  v45.12     - 03/17/06 - kd/rk/jgf - Major additions by kd/rk:              *
C               Added 3D baroclinic transport code in transport.F. Added new  *
C               BPG subroutine in timestep.F. Updated cstart.F, global_3dvs.F *
C               and read_input.F to handle associated input data. Fixed bug   *
C               in VSSOL that was causing inadvertent mixing of current and   *
C               future timelevels. Fixed bug in VSDISP that caused an         *
C               erroneous early return. 3D BAROCLINIC NOW WORKS.              *
C               jgf: Fixed bug in read_global3dvs to handle 3D station        *
C               locations properly (SL0 and SF0). Added code to adcprep to    *
C               process new input data for 3D baroclinic cases. Added MPI     *
C               calls in timestep.F in BPG3D() to pass BPGs to ghost nodes    *
C               for 3D baroclinic in parallel (now works also). Added code    *
C               to read in vert. 3D spacing and vert. eddy vis. in adcprep.   *
C               Applied tjc's DIRNAME patch for ESMF compatibility. Rewrote   *
C               adcprep routine for fort.20, as the original did not handle   *
C               ghost nodes. Fixed another bug in fort.20 adcprep to widen    *
C               the formatting for fluxes. Added jjw's 3 tier tau0 scheme.    *
C               Added bde's bluegene cmplrflags. Increased modularity by      *
C               subroutinizing station locations and adcprep file opening.    *
C  v45.11     - 02/02/06 - jgf - Adcirc: Bugfix to remove "save" compiler     *
C               flag: OnDiag/OffDiag calculation on every entrance now to     *
C               GWCE_New. Bugfix(kd) to correct 3D elevations: flag IFNLFA    *
C               mult. by eta2 in timestep.F. Bugfix to vsmy.F to use IFNLFA   *
C               flag instead of NOLIFA in calculations. Adcprep: Added code to*
C               handle 3D recording stations defined by coordinates rather    *
C               than node number (changed format of fort.80 file), added IDEN *
C               to fort.80 file, fixed routines for post processing fort.41,  *
C               42, 43, 44, 45 and 46 files, added compiler flags for x86_64  *
C               (contrib. tjc) and to "save" subroutine values, fixed NINT    *
C               bug in read_input.F that made 3D output off by 1 time step.   *
C  v45.10     - 01/12/06 - jgf - Adcprep: partially updated preprocessing     *
C               routines to reflect changes in 3D input files.                *
C  v45.09     - 01/09/06 - jgf - Adcirc: fixed typo bugs in global_3dvs.F that*
C               prevented 3D output (NOT3DSDS->NTO3DSDS). Adcprep: added new  *
C               logical variable for aperiodic vs periodic flow boundaries to *
C               fix fort.20 breakup bug.                                      *
C  v45.08     - 12/16/05 - jgf - Adcirc: eliminated 'formatted' statements    *
C               for binary files, relying instead on compiler flags for       *
C               record length, moved setting of NOLIFA flag in adcirc.F to    *
C               fix hotstart bug, changed declarations to COMPLEX(SZ) to fix  *
C               bug in 3D hotstart.                                           *
C  v45.07     - 11/21/05 - jgf - Adcirc: added option to NHSTAR for testing   *
C               adcirc's hot start capability. Prep/post: Wrote code to       *
C               handle fort.71, fort.72, fort.73; upgraded REALs to proper    *
C               precision; modified fort.80 file format to handle NOFF array  *
C               in hotstart files; created post processing options to either  *
C               prepall or just prep14/15.                                    *
C  v45.06     - 09/07/05 - jgf - Found/fixed bug in timestep.f: TotalArea     *
C               was erroneously zeroed. Changed names of some LOGICALs for    *
C               greater consistency in naming of solution options. Selected   *
C               set of default options, including corrected area integration. *
C               Discovered and fixed bug for parallel when NCCHANGE flag was  *
C               set in one subdomain but not another. Updated pre and post    *
C               processing routines, including merging changes written by     *
C               MEB to handle fort.20 files.                                  *
C  v45.05     - 08/1- rl - from 45.04 - for NWS=6, eliminated adding 360 to   *
C               negative longitudes to match AVN model grid setup.  User is   *
C               now required to provide met and ADCIRC grid that correspond   *
C               in space without adjusting the longitude.  Also the input     *
C               variable order has been changed to U,V,P to be consistent     *
C               with other NWS input formats.                                 *
C               Modules affected: wind.f                                      *
C  v45.04     - 06/16-21/05 - rl - from 45.03 - compute baroclinic press grad *
C               from surf down, set = above value if into ground.             *
C               read in sigma grid (IGC=0) or E.V. values (IEVC=0) from       *
C               fort.15 file.                                                 *
C               setup initial condition files for baroclinic or transport     *
C               2DDI, 3D runs                                                 *
C               revised GWCE and 2D Momentum Eq. baroclinic calcs.            *
C               made new GWCE and 2D ME formulations the default              *
C               Modules affected: vssol.f, read_input.f global.f global_3D.f  *
C               timestep.f                                                    *
C  v45.03     - 06/08/05 - rl - from 45.02 - fixed 3D bugs: compute QX, QY in *
C               3D VSSOL, corrected name BcPressNN in ZSURFBuoy, added EP to  *
C               fort.67/68 hot start files & eliminated fort.99 file,         *
C               fixed roundoff error for ievc=10,11                           *
C               Modules affected: vssol.f  timestep.f  hstart.f  global_3dv.f *
C  v45.02     - 03/18/05 - rl - fixed EQUIVALENCE compile problem, changed    *
C               alpha header inputs to 80 characters long, cleaned up Kendra's*
C               pred/corr algorithm, including re-localizing OffDiag, OnDiag  *
C               Modules affected: global.f timestep.f read_input.f            *
C  v45.01     - 10/  /04 - rl - v44.19t Kendra's pred/corr algorithm          *
C               reinstated jjw's spatially variable Tau0 algorithm            *
C               Modules affected: read_input.f, global.f                      *
C  v44.19a-t  - 06/04    - jjw,rl - from 44.18 - new wetting and drying check *
C               NOFF added to hotstart file                                   *
C               Modules affected: adcirc.f cstart.f hstart.f                  *
C  v44.18     - 05/26/04 - rl - from 44.17d - small change to EP check,       *
C               overtopping barriers, baroclinic pressure gradient calc, and  *
C               IMAG to AIMAG in messg.                                       *
C               Modules affected: timestep.f  messenger.f                     *
C  v44.17d    - 05/24/04 - rl - from 44.17c - change variable declaration in  *
C               preprocessor so they are not hardwired to real*4.  Rather,    *
C               they are hardwired to real*8 or value of SZ as determined by  *
C               compiler flag.  No change to ADCIRC.                          *
C               Modules affected: presizes.f  pre_global.f                    *
C  v44.17c    - 05/21/04 - rl - from 44.17b - change output format in post    *
C               processing for fort.63, fort.64 files.  No change to ADCIRC.  *
C               Modules affected: post.f                                      *
C  v44.17b    - 05/18/04 - rl - from 44.17a - vic parrs fixes to prep22 and   *
C               prep23 to write wind info to ghost nodes. No change to ADCIRC *
C               Modules affected: prep.f                                      *
C  v44.17a    - 05/18/04 - rl - from 44.17 - islands kept as islands when they*
C               are cut by subdomains in parallel run. No change to ADCIRC    *
C               Modules affected: prep.f                                      *
C  v44.17     - 05/03/04 - rl - from 44.16 - bug fix for 2 part lateral stress*
C               calculations in GWCE.  Bug dates back to 44.15. New spatially *
C               variable Tau0 formulation.                                    *
C               Modules affected: timestep.f read_input.f                     *
C  v44.16     - 04/30/04 - rl - from 44.15 - added symmetric lateral stress   *
C               and Smagorinski lateral stress coefficient                    *
C               Modules affected: timestep.f  read_input.f  global.f          *
C  v44.15     - 04/16/04 - rl - from 44.13 - cleaned up implementation of all *
C               conservative, lateral stress, ME area integration and mass    *
C               lumping options. Added IM values to access these.             *
C               Added CBaroclinic logical & keyed baroclinic calculations to  *
C               it. Moved LOGICAL declarations from sizes.f to global.f       *
C               Added local LOGICAL CHARMV call and declaration to subroutines*
C               LSQSOLEG & LSQSOLVG in harm.f                                 *
C               Removed CLUMP logical from global_3dvs                        *
C               Fixed VSMY so that DASigT is computed for baroclinic run.     *
C               Modules affected: global.f  read_input.f  timestep.f  sizes.f *
C                                 harm.f  adcirc.f  global_3dvs  vsmy.f       *
C                                 messenger.f                                 *
C  v44.13     - 03/16/04 - rl - from 44.12 - fixed bug in ME TotalArea that   *
C               that was introduced in v44.09                                 *
C  v44.12     - 03/15/04 - rl - from 44.11 - fixed bug in GWCE lateral stress *
C               that was present in v44.10 & v44.11                           *
C  v44.11     - 03/13/04 - rl - from 44.10 - revamped spatially varying Tau0  *
C               Modules affected: read_input.f                                *
C  v44.10     - 03/09/04 - rl - from 44.09 - fixed bug in GWCE lateral stress *
C  v44.09     - 03/08/04 - rl - from 44.08 - changed lateral stresses in GWCE *
C               from Randy's method to straight forward computation           *
C               THIS VERSION HAS A BUG CORRECTED IN VERSION 44.10             *
C  v44.08     - 02/25/04 - rl - from 44.07 - moved GWCE to subroutine in      *
C               timestep.f                                                    *
C               Modules affected: global.f  read_input.f  timestep.f          *
C  v44.07     - 12/20/03 - rl - from 44.06 - major bug fix in conservative    *
C                                          GWCE formulations                  *
C                                          Modules affected                   *
C                                          timestep.f                         *
C  v44.06     - 12/17/03 - rl - from 44.05 - added 2nd nonconservative momen- *
C                                          tum Eq. (correct area integration) *
C                                          Modules affected                   *
C                                          timestep.f                         *
C  v44.05     - 10/19/03 - rl - from 44.04 - implemented 2nd conservative     *
C                                          momentum formulation using product *
C                                          rule.                              *
C  v44.04     - 10/09/03 - rl - from 44.03 - moved 2DDI momentum calculations *
C                                          to subroutines in timestep.f module*
C                                          Modules effected                   *
C                                          timestep.f                         *
C  v44.03     - 10/07/03 - rl - from 44.02 - implemented conservative momentum*
C                                          equation in main timestepping code *
C                                          continued removing local variables *
C                                          from global module.                *
C                                          Modules affected                   *
C                                          adcirc.f   cstart.f      global.f  *
C                                          hstart.f   read_input.f  timestep.f*
C  v44.02     - 09/30/03 - rl - from 44.01 - additional changes associated    *
C                                          with removing local variables from *
C                                          global module.                     *
C  v44.01     - 10/06/03 - rl - from 43.03 - changed names of load vectors in *
C                                          2D solution                        *
C  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
C                                          output buffer flush (m.cobb)       *
C                                          3D fixes (k.dresback)              *
C                                          drop MNPROC in fort.15 (t.campbell)*
C                                          various bug fixes in RBCs          *
c                                          ZSURFBUOY/BCPG calc                *
C  v43.02     - 02/06/03 - rl - from 43.01 - code clean up & documentation    *
C  v43.01     - 01/31/03 - jf - from 43.00 - reconcile with v42.07 2D code    *
C                                            changed var. names: NNEIGH->     *
C                                            NNeigh, NEIGH->NeiTab,           *
C                                            NEIGHELE->NeiTabEle,             *
C                                            UBAR1->UU1,                      *
C                                            VBAR1->VV1, UBAR2->UU2, VBAR2->  *
C                                            VV2, TAUSX1->WSX1, TAUSY1->WSY1, *
C                                            TAUSX2->WSX2, TAUSY2->WSY2,      *
C                                            TAUBX1->BSX1, TAUBY1->BSY1,      *
C                                            TAUBXN1->BSXN1, TAUBXN2->BSXN2,  *
C                                            TAUBXN3->BSXN3, TAUBYN1->BSYN1,  *
C                                            TAUBYN2->BSYN2, TAUBYN3->BSYN3,  *
C                                            TAUBXPP3->BSXPP3,TAUBYPP3->BSYPP3*
C                                               parallel unified 2D/3D source *
C  v43.00b    - 12/11/02 - jf - from 43.00a - bug fixes                       *
C  v43.00a    - sum  /02 - tc - from 36.01 (3D) & 41.12? (2D), create F90/    *
C                                            parallel unified 2D/3D source    *
C  v42.07     - 08/12/02 - rl - from 42.06 - combined discharge & wave rad bc *
C  v42.06     - 06/24/02 - rl - from 42.05 - redid no normal vel grad bc v1   *
C  v42.05     - 06/18/02 - rl - from 42.04 - combined v1 & v2 no normal vel   *
C                                                                   grad bcs  *
C  v42.04     - 06/08/02 - rl - from 42.02 - fixed bug in no normal vel grad  *
C                                                                   bc v1     *
C  v42.03     - 06/06/02 - rl - from 42.01 - no normal velocity grad bc v2    *
C  v42.02     - 06/05/02 - rl - from 42.01 - no normal velocity grad bc v1    *
C  v42.01     - 04/12/02 - rl - from 41.11 - redid GWCE formulation           *
C  v41.11     - 09/14/01 - rl - from 41.10 - added NWS = -2, eliminate MNWLON *
C                               and MNWLAT from NWS11GET                      *
C  v41.10     - 07/25/01 - rl - from 41.09 - bug fix in GWCE lateral viscosity*
C  v41.09     - 06/30/01 - jw - from 41.08 - minor mods per vp version 41.05  *
C  v41.08     - 06/22/01 - rl  - reconciled v41.07 and v41.05m009             *
C                                            changes to HABSMIN and ETA2      *
C  v41.07     - 04/09/01 - rl - from 41.06 - initialized PRN1(),PRN2() for    *
C                                                                      NRS<>0 *
C  v41.06     - 04/02/01 - rl - corrected dimensioning on most of wind reading*
C                               subroutines.  Changed MNWP to MNP in wind     *
C                               forcing ALLOCATION statements.  Note MNWLON   *
C                               and MNWLAT still need to be eliminated from   *
C                               NWS11GET subroutine.                          *
C  v41.05m009 - 06/12/01 - rl&jjw - updated NODECODE after initial drying sect*
C  v41.05m008 - 06/01/01 - jjw - changed code by adding output in wet/dry sect*
C  v41.05m007 - 05/28/01 - jjw - changed HABSMIN=0.8D0*H0                     *
C  v41.05m006 - 05/22/01 - jjw - add writes to track wetting/drying bug       *
C  v41.05m005 - 05/21/01 - jjw - modified wet/dry to allow barrier overtopping*
C                                modified from version 41.05.original         *
C  v41.05m004 - 02/15/01 - jjw - from 41.05.m003 - added logic to reset min   *
C                               depth behind overtopping barriers from wetting*
C  v41.05m003 - 02/14/01 - jjw - from 41.05.original - skipped over version   *
C                                41.05m001/m002. Fixed bug in wet/dry         *
C                                algorithm which prevented overtopping        *
C                                internal barriers from wetting               *
C  v41.03       09/15/00 - rl - added bridge piling friction, fixed several   *
C                               F90 bugs                                      *
C  v41.02       09/07/00 - rl - fixed F90 hot start bug and consolidated with *
C                               version 35.xx                                 *
C  v40_02m004 - 05/02/00 - rl - changed so that contribution to RHS forcing is*
C                               zeroed out for any element that contains a dry*
C                               node. Note: LHS is taken care of automatically*
C  v40.02m003 - 04/28/00 - jjw - Changed wet/dry interface from essential no  *
C                                normal & tangential to natural no normal     *
C  v40.02m002 - 12/22/99 - jjw/vjp - Vic change to avoid compiler problems    *
C  v40.02m001 - 12/21/99 - jjw - add cross barrier pipes cjjwm001             *
C                                                                             *
C                                                                             *
C  8/7/2000 modfied jcg_with_blas -- fixed bug in itjcg for single cpu case   *
C  8/7/2000 modfied wetting and drying section to fix bug and optimize        *
C  11/1/2000 added init of TIP2 in hstart and fixed Sun makefile for Jesse    *
C 
C******************************************************************************

