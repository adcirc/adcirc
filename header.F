C******************************************************************************
C           This file contains a log of changes made to (P)ADCIRC.            *
C                                                                             *
C  mod history                                                                *
C  v50.99.13 - 12/04/2013 - jgf - Fixed incorrect station coordinates
C              in all station netcdf files, 2D and 3D, serial and parallel,
C              where Cartesian coordinates were erroneously being written
C              to station netcdf files, coordinates are now correctly 
C              written in degrees latitude and longitude.
C  v50.99.12 - 11/11/2013 - jgf - Fixed f90 compilation rule and arguments to 
C              switch.pl (replacing -unix -f95 with swch) for compiling adcswan
C              and padcswan targets.
C  v50.99.11 - 11/01/2013 - tcm - Corrected a bug in adcprep for nws=2 where
C              wtiminc was not being correctly read from the fort.15 file.
C  v50.99.10 - 07/25/2013 - jgf - Corrected a bug in NWS8 (symmetric holland)
C              where the background pressure was being erroneously applied in 
C              mb whereas NWS8 requires Pa; fixed met-only mode by forcing
C              a recalculation of starting and ending timesteps for output
C              whenever the ADCIRC timestep has changed after a hotstart;
C              fixed issue where variable time step after hotstart was 
C              restricted to 3D baroclinic situations; fixed bug that was 
C              causing NWS8 to run out of met data one met time step too early. 
C  v50.99.09 - 05/14/2013 - jgf - Corrected 3D netcdf hotstart issue in
C              parallel where typos in netcdfio.F were causing density3D to be
C              read erroneously in place of velocity3D and turbulence3D data.
C  v50.99.08 - 05/06/2013 - jgf - Moved the logical variable writerReadMetaData
C              from global.F to netcdfio.F to avoid compilation issue on IBM
C              Power7 (AIX 6.1).
C  v50.99.07 - 05/06/2013 - jgf - Updated developers documentation to include
C              instructions for building ADCIRC (no code changes).
C  v50.99.07 - 05/01/2013 - jgf - Added a comma to a write statement in
C              hstime.F to avoid a syntax error when NetCDF is not enabled.
C  v50.99.06 - 05/01/2013 - tcm -- Added kdtree2 module to swan Makefile, which
C              is the same fix that Casey added to the trunk.
C  v50.99.05 - 02/03/2013 - zc - Corrected variable FoundEye not being set in 
C              NWS15GET, so wind drag always defaulted to Garratt formulation.
C  v50.99.04 - 01/16/2013 - jgf - Fixed issue found by JieGao in initialization
C              of quadrant angles that affected Vr calculations and caused
C              differences between calculated quadrant wind speeds and isotach 
C              values from the track file. 
C  v50.99.03 - 01/10/2013 - zc - Correcting bug that would cause error if
C              elemental_slope_limiter was not last attribute in fort.13.
C              Added check for FoundEleSlopeLim and LoadEleSlopeLim as well.
C  v50.99.02 - 11/29/2012 - jgf - Corrected a typo in HWind (NWS15) processing
C              that prevented the y-component of the wind speed from being
C              converted from 1 minute to 10 minute winds. 
C  v50.99.01 - 10/23/2012 - zc - Reversing change from 9/28/2012. This is
C              actually a leftover debugging statement in NetCDF and is corrected
C              in versions > 4.2. Also, this seemed to cause segfaults under
C              some compilers.
C  v50.99    - 10/08/2012 - jgf - Added if/then statements so that minpr.63
C              and maxwvel.63 won't be written if there was no wind forcing;
C              also, so that maxrs.63 won't be written unless there was 
C              wave forcing. 
C  v50.98    - 10/04/2012 - jgf - Added subroutine to write LHS of harmonic
C              analysis matrix for debugging purposes; fixed issue in the 
C              way that the LHS matrix is read from the binary hotstart 
C              file; fixed issue with the way that harmonic analysis
C              load vectors are written to the netcdf hotstart file; fixed
C              an issue with the on-the-fly decomposition of harmonic 
C              analysis load vectors from netcdf hotstart files; fixed 
C              issue with reading harmonic analysis load vectors for 
C              station velocity harmonics from binary hotstart file.
C  v50.97    - 09/28/2012 - zc - Corrected how the element array is written
C              to NetCDF files to prevent unnecessary screen output for some
C              builds of NetCDF. Also made sure negative values are treated
C              correctly in couple2swan max file writing.
C  v50.96    - 09/27/2012 - jgf - Removed stray debugging messages from 
C              netcdfio.F that were causing segfault upon hotstart with netcdf
C              files. 
C  v50.95    - 09/27/2012 - zc - Added SwanMaxDescript%Specifier = NOUTGW to
C              ensure NetCDF format is written for SWAN max files. 
C              Corrected SWAN output routines to always use writer processors
C              if MNWPROC > 0. Behavior is consistent with write_output.F
C  v50.94    - 09/24/2012 - jgf - Added DA compiler flags (-DREAL8 etc) to
C              compile hstime, aswip, and adccmp; fixed read format issues in 
C              nws15get; fixed cycle time dimension issue, globaldir
C              initialization issue, and lat/lon initialization in aswip
C              (these issues only affect post storm analysis for hwind). 
C  v50.93    - 09/24/2012 - zc - Fixed a bug when both writing and reading
C              netcdf hot start files where timenc was deallocated but not 
C              marked for reinitialization later.
C  v50.92    - 09/24/2012 - jcd - Fixed a bug in the NetCDF output of the sigma
C              layers in 3D.  The variable NFEN was declared in pre_global.F,
C              and thus was not available to the FEGRIDS routine within the
C              global_3dvs.F module.  Changed pre_global.F to use the NFEN
C              declared within the global_3dvs.F module.
C  v50.91    - 09/21/2012 - jcd - Merged the swan40.91 branch into the trunk.
C  v50.90    - 09/18/2012 - tcm  Fixed a typo bug in adcprep.F for the windows
C                 builds.
C  v50.89    - 09/16/20120 - zc - Merging netcdf_writer branch into trunk,
C              revisions 808:836 and resolved conflicts. Ensured that 
C              rads.64.nc could be created without CSWAN flag.
C     v50.74.06   09/12/2012 - zc - Corrected SWAN output for max wind
C                 velocities in netcdf to scalar format.
C     v50.74.05   09/12/2012 - zc - Corrected a bug in the SWAN output
C                 when not using NetCDF writers. Allow the writer
C                 processors to read metadata from ADCPREP created
C                 output files during first pass since they do not read
C                 fort.15. Tested for NETCDF3 and NETCDF4. Also, corrected
C                 SWAN max files in NetCDF which were still writing in
C                 ASCII.
C     v50.74.04   07/01/2012 - zc - Removing unused variables added for
C                 original implementation of NetCDF writer cores
C     v50.74.03   06/29/2012 - zc - Remove now unused MPI_COMM for 
C                 inter-writer only communication. Adjust makefile
C                 for correct dependencies for netcdfio.F and writer.F
C     v50.74.02   06/29/2012 - zc - NetCDF writers seem to be working
C                 properly for ADCIRC and SWAN. Continued testing needed 
C                 before merge to main trunk. 
C     v50.74.01   06/20/2012 - zc - Initial edits for NetCDF writer cores.
C                 Writer cores are recognized under ADCIRC writes only
C                 at this point. Need to ensure that writer cores have 
C                 knowledge of all NetCDF record numbers.
C  v50.88   - 09/14/2012 - zc - Corrected a single processor arriving at
C              LBCAST in couple2swan.F. Made a similar change for
C              consistency in write_output.F, but using ADCIRC_TERMINATE().
C  v50.87   - 09/14/2012 - zc - Corrected the size of the netcdf_error array
C              in couple2swan.F to prevent message passing errors.
C  v50.86   - 09/12/2012 - tcm - Adding an empty directory "windows" where
C              MS Visual Studio projects for building metis/adcprep/adcirc/
C              and padcirc can be stored.
C  v50.85   - 09/11/2012 - tcm - Added a few changes related to compiling
C              on Windows PC using MS Studio and Intel MPI.
C     v50.84.05   - 09/11/2012 - jcd - Corrected the use of MPI_COMM_WORLD
C                   when it should instead have been COMM in SwanPunCollect.
C     v50.84.04   - 09/09/2012 - jcd - Two small changes to the existing files
C                   for src/couple2swan.F and work/makefile.
C     v50.84.03   - 09/09/2012 - jcd - Added the newest release version 40.91
C                   in the swan sub-directory.  Carried forward our macros.inc
C                   files with the compiler flags for different machines.
C     v50.84.02   - 09/09/2012 - jcd - Trying to delete the existing SWAN.
C     v50.84.01   - 09/09/2012 - jcd - Fixed a small bug in write_output.F with
C                   the writing of the new NetCDF output files.
C  v50.84   - 08/16/2012 - jgf - Fixed aswip compilation dependencies to enable
C             the building of aswip on garnet at ERDC.
C  v50.83   - 08/16/2012 - jgf - Fixed (removed or changed) stray references 
C             to screen() subroutine in ADCIRC code. 
C  v50.82   - 08/15/2012 - jgf - Adjusted routines related to a single 
C             subdomain initiating code exit: if a subdomain experiences
C             a fatal error, it calls adcirc_terminate() as usual, and 
C             sets a flag that causes mpi_abort() to be called in msg_fini()
C             (after writer processors have been shut down); also added
C             execution tracing to messenger.F and fixed a conflict between
C             ADCIRC and SWAN by renaming ADCIRC's screen() subroutine to 
C             screenMessage().
C  v50.81   - 08/01/2012 - zc - Added routines so a single processor can
C             initiate code exit. ADCIRC_LOCALTERMINATE and
C             ADCIRC_CHECKLOCALTERMINATE. All processors use the check routine
C             and the processor with the error uses LOCALTERMINATE and waits
C             until the others catch up. Added this for error checking in
C             read_input.F where mesh or fort.15 errors that apply to a single
C             processor only will hang the code since other processors have
C             no knowledge of the error and hang. Readily extendable to 
C             other areas of the code as needed.
C  v50.80   - 07/31/2012 - jgf - Added underscore to mesh topology attribute
C             in netcdf files for NcML compliance.
C  v50.79   - 07/20/2012 - tcm - There was a bug in the prep.F routine in the
C             writing of the metControl namelist that was causing it to not
C             be written correctly with default values.  The writing of this
C             namelist was corrected to remove the bugs.  Then logic was 
C             added so that the metControl namelist is only written in the
C             partitioned fort.15 files if it was present in the global
C             fort.15.  Changed some float() to real(,kind=sz) to remove
C             warning messages.
C  v50.78   - 07/20/2012 - tcm - Changed real*8 to real(8) to remove
C             warning messages.
C  v50.77   - 07/20/2012 - tcm - Fixed a bug in prep.F where the default
C             filename was not initialized to blanks for certain compilers
C             and therefore could have extra characters that were not a part
C             of the acutal name, for example a file name might have been
C             fort.13(.  This error happened when default name lengths 
C             changed from 7 to 8 characters to allow for larger unit 
C             numbers.  Added the Utility Server (utils) compiler options
C              back to the cmplrflags.mk file.
C  v50.76   - 07/19/2012 - jgf - Fixed missing parenthesis on line 
C             7062 in prepNetCDF() where the code was checking if 
C             NRS eq 3 or 4; this issue was preventing compilation.
C  v50.75   - 07/16/2012 - tcm - Implemented Casey's writing of rads.64 file
C             whenever nrs = 3 or nrs = 4.  Removed the ifdef cswan and in
C             most cases replaced with testing of nrs values.
C  v50.74   - 06/11/2012 - jgf - Merged f90_netcdf branch revisions 764:806 
C             and resolved conflicts. Fixed two previously undetected 
C             minor issues: (a) SL1 and SF1 arrays were not populated in 
C             read_global.F for ICS=1, but they are required for netcdf; and
C             (b) subdomain to fulldomain station mappings were using wrong
C             dimensions during decomposition of harmonic means and variances
C             when hotstarting from a netcdf hotstart file. 
C     v50.60.21  05/24/2012 - jgf - Fixed a floating point exception when
C                writing 3D velocity stations in parallel. 
C     v50.60.20  05/24/2012 - jgf - Fixed two issues with the new metControl
C                and SWANOutputControl namelists so that (a) they are not 
C                required in the fort.15, and (b) so that the metControl 
C                namelist is written with the draglawstring properly enclosed 
C                in single quotes by adcprep upon write to subdomain fort.15s 
C                so that padcirc can read it successfully.  
C     v50.60.19  05/22/2012 - casey - Corrected a bug with the output of the
C                SWAN max files.
C     v50.60.18  05/22/2012 - jgf - Fixed dependence of presizes.F on wind
C                module; made changes to support gfortran without underscoring;
C                fixed issue with 3D bottom roughness lengths, where Z0B
C                was declared in both nodalattr.F and global_3dvs.F, resulting
C                in Z0B being used uninitialized in the apply3DBottomFriction()
C                subroutine; fixed numerical issue in wetting and drying 
C                due to DELETA sometimes being negative; fixed issue so that
C                maxrs.63 file is still produced, even without compiling in
C                SWAN coupling (purely to match past behavior); replaced some
C                stray f77-style subroutines in netcdf module; fixed
C                numerical issue with translation speed calculation in 
C                symmetric vortex (NWS8) subroutine.
C     v50.60.17  05/17/2012 - corbitt - Fixed the vector max netcdf output so
C                it was outputing scalar magnitudes instead of vectors.
C     v50.60.16  05/17/2012 - jgf - Added tracing to writer.F and 
C                couple2swan.F; removed references to netcdff library from
C                makefile in favor of adding -lnetcdff to FLIBS in
C                cmplrflags.mk for platforms that require it; fixed 
C                array out-of-bounds in netcdf_error array by increasing
C                array length to 14. 
C     v50.60.15  05/17/2012 - corbitt - Added netcdf output for maxele.63,
C                maxvel.63, minpr.63, maxrs.63, maxwvel.63. Format is dependent
C                on NOUTGE, NOUTGV, NOUTGW, so max files will be formatted same
C                as fort.63,64,73,and 74. If NOUTGE, NOUTGV, NOUTGW are set
C                to 0, standard ascii is output. Population of maxes from
C                hotstarts are still dependent on ascii format. In addition, 
C                made the 2D baroclinic and swan netcdf output CF compliant.
C                Also preliminary set up netcdf max files for swan output.
C     v50.60.14  05/15/2012 - jgf - Fixed typo that was causing array 
C                out-of-bounds when writing fort.45 in netcdf format.
C     v50.60.13  05/15/2012 - jgf - Moved metControl namelist for specification
C                of air density, drag law, and wind drag limit to end of 
C                fort.15, instead of reconstructing the NWS line; this 
C                brings this data structure into line with the 
C                TimeBathyControl and SWANOutputControl namelists. Also fixed
C                an issue in checking NWS values that was rejecting negative
C                values of NWS. 
C     v50.60.12  05/14/2012 - jgf - Fixed issue with compiler flags on 
C                blueridge at RENCI to use the path to netcdf4.1 (required
C                when using NETCDF4=enable).  
C     v50.60.11  05/14/2012 - jgf - Fixed issue in makefile to resolve
C                reference to global.F in presizes.F.
C     v50.60.10  05/14/2012 - jgf - Added code to prep.F to calculate and
C                save sigma coordinates to 3D netcdf file; moved FEGRIDS
C                subroutine from vsmy.F to global_3dvs.F so that it could
C                be used by adcprep; added an IF statement in couple2swan.F
C                so that only processor 0 would make the subroutine calls
C                to initialize the SWAN-related netcdf output files; also
C                changed SWAN_VOQR(IVTYPE+1) to SWAN_VOQR(IVTYPE)+1 in
C                couple2swan.F to remedy an array-out-of-bounds issue
C                (this appeared to be a typo). 
C     v50.60.09  05/12/2012 - jgf - Added autocompression for netcdf4 
C                formatted SWAN output files; fixed issue that was causing
C                station coordinates to be reported in CPP rather than 
C                geographic projection in netcdf output files; added optional
C                namelist to bottom of fort.15 for reading zc's logical
C                variables for turning SWAN output on and off; fixed an 
C                issue in adcprep that was trying to open the SWAN command 
C                file twice if the SWAN command file was named fort.26;
C                removed assumptions in symmetric vortex model (nws=8) about
C                northern and western hemispheres so it can now simulate
C                storms in eastern and southern hemispheres; fixed file
C                reading issue in symmetric vortex model where it was 
C                trying to read subdomain fort.22 files; removed stop
C                keyword from symmetric vortex subroutine that could
C                cause problems in parallel operation.
C     v50.60.08  05/10/2012 - zc - Added netcdf output to SWAN output files
C                and rads.64. The swan output files are now turned on and off
C                using logicals in the global.F file. Output files are set
C                to the defaults that were previously hardwired.
C     v50.60.07  05/06/2012 - jgf - Corrected issue with writing 3D fulldomain
C                velocity files in netcdf in serial; added storage of sigma
C                array to 3D netcdf files.  
C     v50.60.06  05/05/2012 - jgf - Fixed issues in netcdfio.F related to 
C                setting vertical dimension when initializing and writing 3D 
C                station and 3D fulldomain output files in netcdf format. 
C     v50.60.05  05/04/2012 - jgf - Simplified structure of netcdfio.F 
C                to reduce complexity and take advantage of commonalities
C                in the i/o of the various station and full domain files.
C                Tested netcdf output in serial and in parallel, including 
C                serial and parallel hotstarting with netcdf hotstart files
C                on 2D quarter annulus case. All smoke tests pass at this
C                point; a few more tests are required. 
C     v50.60.04  05/03/2012 - jgf - Eliminated memory allocation and copying of
C                3D data and 3D hotstart data in netcdfio.F that is no longer
C                required when using f90-style subroutine calls. Added code
C                that pulls 3D data from netcdf hotstart files directly into
C                adcirc variables in serial and in parallel. This code 
C                compiles; next step is testing.  
C     v50.60.03  05/02/2012 - jgf - Fixed an incomplete implementation of 
C                hotstarting in parallel with netcdf that was uncovered 
C                in testing previously committed changes. Applied logging
C                fix from trunk (v50.65). 
C     v50.60.02  05/01/2012 - jgf - Migration of adcirc netcdfio module 
C                to f90-style subroutine calls is complete. It now compiles too.
C                The next step is to make corresponding changes to 
C                write_output.F.
C     v50.60.01  04/29/2012 - jgf - Halfway through changing the subroutine
C                calls in netcdfio.F to f90-style. This is a work in progress.
C  v50.73   - 05/08/2012 - tcm - Removed use of reserved work (scale) from the
C             file diffmerge.F.
C  v50.72   - 05/07/2012 - tcm - Incorporated time varying bathymetry changes
C             into ADCIRC. These came from the time_varying_bathy branch.
C             If time varying bathymetry is not used, then no changes to the 
C             fort.15 occur.  To use time varying bathymetry, then add the
C             "timebathycontrol" namelist at the end of the fort.15 file.
C             Time varying bathymetry uses a fort.141 file to supply the
C             different bathymetry values and as output writes a fort.75
C             (station) and fort.76 (global) that are controled by the
C             fort.61 and fort.63 file options.
C  v50.71   - 05/01/2012 - zc - Reorganized ps2dots and ps3dots routines
C             to improve code performance. Changed ADVECTLOCAL to only
C             be called when advection is turned on, which improves code
C             speed. Moved call to updater to after slope limiter
C             routine to avoid calling it twice.
C  v50.70   - 04/22/2012 - zc - Allowed the values from ESLNodes.63 to 
C             be used in a hotstart run to reactivate slope limiter at
C             previously active nodes. Fix hot start allocation issues.
C  v50.69   - 04/16/2012 - jgf - Applied fix from Kendra Dresback for 
C             erroneous format in writeOutHarmonicArrays() subroutine
C             that was causing phase values to be off by a factor of 10. 
C  v50.68   - 04/15/2012 - jgf - Removed water surface elevation smoother
C             nodal attribute. 
C  v50.67   - 04/15/2012 - zc - Corrected the elemental slope limiter to 
C             apply smoothing if fort.13 max gradient is set to zero
C             instead of skipping the node. Also set so that a negative 
C             fort.13 value triggers gradient warnings but no slope
C             limiting
C  v50.66   - 04/09/2012 - jgf - Corrected issue in reading NWS namelist
C             that caused a compiler warning about out-of-block goto. 
C  v50.65   - 04/09/2012 - jgf - Corrected issue with logging where all 
C             subdomains were logging to an erroneous fulldomain fort.16.
C  v50.64   - 04/08/2012 - jcd - Updated a few SWAN source files to reflect
C             the latest changes by the SWAN developers.
C  v50.63   - 04/07/2012 - zc - corrected the elemental slope limiter
C             output statements so they only appear once per node instead of
C             continuously until the gradient is below the threshold
C  v50.62   - 04/07/2012 - zc - corrected my previous update to use -99999
C             when a node is dry instead of what resides in eta2
C  v50.61   - 04/05/2012 - zc - fixed a bug when writing the sparse format
C             in the fort.69 debug file. Switched to full format for
C             simplicity.
C  v50.60   - 04/05/2012 - jgf - Merged newvortex branch revisions 705 to 760
C             into the trunk (the newvortex branch was originally created 
C             from v50.38 of the ADCIRC trunk). This merge contains the new
C             NWS19 asymmetric vortex model, new aswip met preprocessing and
C             analysis code, and new NWS15 hwind input capability. Detailed 
C             changelog is shown below (v50.38.01 through v50.38.09):
C     v50.38.09  04/05/2012 - jgf - Removed extraneous use of globaldir from
C                aswip code.
C     v50.38.08  04/05/2012 - jgf - Corrected calculation of wind magnitude
C                for nws15 (hwind) in timestep.F; changed log level to echo
C                for logging drag law string from nws line in fort.15.
C     v50.38.07  03/25/2012 - jgf - Changed windMultiplier variable to
C                hWindMultiplier to avoid variable conflict between hwind and
C                owiwind modules.
C     v50.38.06  03/10/2012 - jgf - Added vmax and rmax vars to hWindData_t
C                block; added support for Powell wind drag to NWS15 (HWind);
C                added support to NWS15 for finding vmax and rmax; reorganized
C                aswip code to make it modular; added capabilities to the
C                wind analysis features of aswip to enable analysis of HWind
C                data via ADCIRC's new HWind support.
C     v50.38.05  03/05/2012 - jgf - Created NWS=15 for reading HWind data,
C                including nws15init, nws15get, loadHWindFile, and
C                hWindInterpolateVelocity subroutines, as well as supporting
C                changes in adcirc and in various routines in adcprep as needed;
C                code has been tested and appears to work correctly thus far.
C     v50.38.04  01/24/2012 - jgf - Removed dependency on precision.F, constants.F
C                fitpack.F, and lsq.F from code; removed NWS9; fixed issue
C                with percent Rmax specification in aswip.
C     v50.38.03  01/16/2012 - jgf - Changes to simplify the aswip code and
C                integrate it more closely with ADCIRC; made air density a
C                variable that can be specified by the user in the fort.15;
C                made the default value of the air density 1.15.
C     v50.38.02  01/11/2012 - jgf - Latest changes to account for storm motion
C                on Vmax and isotach wind speeds for NWS19; addition of code to
C                aswip for reading/writing additional data for visualizing
C                wind output on autogenerated unstructured mesh; addition of
C                nearest neighbor weighted interpolation to vortex.F for
C                interpolating Rmax between quadrants; some further simplification
C                of vortex.F.
C     v50.38.01  12/02/2011 - jgf - Created branch to hold ongoing changes that
C                have been made to assess the meteorological fields generated
C                by the asymmetric vortex model (NWS19), including the following:
C                very heavy modification and simplification of vortex.F,
C                reducing it from 2200 lines of code to 800 lines of code;
C                addition of code to aswip to support the generation of radial
C                and spatial profiles of wind speed; and modification of the
C                interface between NWS9/NWS19 and vortex.F.
C  v50.59   - 04/04/2012 - jgf - Added preprocessor directives so that
C             netcdfio.F can be compiled with either netcdf3 libraries or
C             netcdf4 libraries; made adcprep read the flux time increment
C             (FTIMINC) in fort.20 as a real instead of an integer; corrected
C             the swan version numbers in the changelog below for version
C             50.55.
C  v50.58   - 04/02/2012 - jcd - Two changes to adcprep: (1) to avoid an
C             endless loop when prepall does not find a default file,
C             now it will ask the user for the file name or to skip; and
C             (2) now it will allow the SWAN control file to have its default
C             name of INPUT, even if the swaninit file is not found.
C  v50.57   - 03/28/2012 - pck - Added two swan files missing from the update
C  v50.56   - 03/28/2012 - pck - Added another ramp option: NRAMP=8, where
C             the new variable DUnRampMete is the time to begin METE ramping.
C             Also added the nodal attribute (advection_state) that allows
C             advection to be initiated based on depth. If the depth is
C             greater than the advection_state depth, advection will be
C             applied (assuming that advection is turned on); otherwise
C             advection will be off.
C  v50.55   - 03/25/2012 - jgf - Merged Casey Dietrich's changes for the
C             latest version of SWAN (moving from v40.72 to v40.85 of SWAN)
C             as well as the capability to limit the spectral propagation
C             velocities in SWAN; the merged code also contains Taylor Asher's
C             new capability for modifying wave results in SWAN to account
C             for ice coverage.
C  v50.54   - 03/25/2012 - jgf - Made changes to enable compilation of netcdfio
C             on intel compilers; upgraded hstime utility to new ADCNETCDF
C             preprocessor variable; added preprocessor variables to turn off
C             automatic netcdf4 compression if the system netcdf4 library
C             does not support it (apparently support for compression was
C             added in netcdf v4.1, but is not available in netcdf 4.0.x).
C  v50.53   - 02/29/2012 - jgf - Removed massbal.F from code and makefile;
C             restored -i-dynamic compiler flag for intel to avoid
C             feupdateenv message; added swan macros.inc for g95 compiler;
C             turned on shuffle filter for compressed netcdf4 files.
C  v50.52   - 02/21/2012 - jgf - Merged the following changes from revisions
C             r728:r738 of the netcdf4 branch; this branch was developed to
C             enable compressed native netdf4 output from ADCIRC using an
C             output specifier of 5, i.e., setting NOUTGE=5 will cause ADCIRC
C             to produce fulldomain elevation output in NetCDF Classic Model,
C             NetCDF4-formatted, compressed output files:
C             v50.43.01  02/20/2011 - jgf - Changes related to netcdf4
C                  output files (output specifier = 5); changed NETCDF
C                  preprocessor definition to ADCNETCDF to avoid collision
C                  with the official NETCDF module, which was introduced
C                  to allow access to the fortran 90 interface for netcdf;
C                  automatically enabled compression for netcdf4 output format;
C                  used f90 interface for saving text attributes to netcdf
C                  files to avoid potential mismatches in character counts
C                  in nf_put_att_text(); added tracing to vsmy.F and
C                  timestep.F to track down memory issues; initialized
C                  BPG to zero in cstart.F; added a convenience subroutine
C                  ADCIRC_Terminate() to adcirc.F.
C  v50.51   - 02/13/2012 - jgf - Merged ADCIRC compiler flags and SWAN
C             macros.inc for Kraken.
C  v50.50   - 02/08/2012 - jgf - Changed organization of smoother and slope
C             limiter calculations slightly for sustainability; made negative
C             value of threshold gradient signify unlimited slope at a node
C             (setting threshold gradient to 0.0 now gives continuously limited
C             slope).
C  v50.49   - 02/08/2012 - zc - Created new nodal attribute for
C             elemental_slope_limiter which takes as a value the maximum
C             gradient to allow before applying the smoothing function which
C             is exactly the same as that used in water_surface_elevation_smooother.
C             ESLNodes.63 is written out at the end of the run showing the
C             nodes which triggered the function.
C  v50.48   - 02/08/2012 - jgf - Made stations allocate memory with MNST*
C             rather than NST* dimensions; moved statement that sets 1 station
C             as the minimum to allocate for concentration and meteorology
C             stations to avoid potential memory issues.
C  v50.47   - 02/07/2012 - jgf - Fixed tab handling issued in subroutine
C             parse() for station names that caused a segfault on pgi compiler.
C  v50.46   - 02/04/2012 - jgf - Created new nodal attribute called
C             water_surface_elevation_smoother and added code to GWCE_New
C             to apply the smoother if it has been loaded; the value of this
C             nodal attribute is the time step increment between applications
C             of the smoother; i.e., 0 represents no smoothing at that node,
C             1 represents smoothing at every time step, 2 indicates smoothing
C             every other time step, etc.
C  v50.45   - 02/04/2012 - jgf - Cleaned up stray references to CGWCE_Original.
C  v50.44   - 02/04/2012 - jgf - Set namelist scratch file to write to localdir
C             to avoid collisions when running in parallel; removed the
C             GWCE_Original() subroutine from timestep.F per Rick's request,
C             fixed 3 instances in timestep.F when STOP was called without
C             shutting down mpi first.
C  v50.43   - 02/03/2012 - jgf - Changed rewind() on scratch file to open()
C             and close() to appease g95 compiler (for nws namelist i/o).
C  v50.42   - 02/02/2012 - jgf - Changed reading of NWS line so that namelist
C             is written and then read from a scratch file rather than
C             an internal character array (fixes issue caught by pgif95).
C  v50.41   - 12/23/2011 - tcm - Added compiler options for ERDC's Utility
C              server (compiler=utils)
C  v50.40   - 12/21/2011 - jgf - Merged the r364:365 of the debug_support
C             branch (v49.06.01), as follows:
C  v49.0601   - 12/09/09 - jgf - Added capability to produce fulldomain fort.69
C               file (i.e., fulldomain elevation file on timesteps where the
C               warning elevation has been exceeded) in addition to existing
C               subdomain fort.69 files.
C  v50.39     12/08/2011 - tcm - Fixed a possible bug where eta2 was supposed
C             to be zeroed out after setting eta1 = eta2 but was not.
C  v50.38   - 11/28/2011 - jgf - Fixed CF compliance of time metadata in netcdf
C             by changing attribute text from seconds since base_date to
C             e.g., seconds since 2010-05-01 00:00:00 UTC (the date is
C             simply copied verbatim from the existing base_date metadata
C             field in the fort.15 file).
C  v50.37   - 09/23/2011 - jgf - Added extra checks to the debugging flags
C             for diamond in cmplrflags.mk; widened some integer format
C             for logging when reading fort.13, fort.14 and fort.15;
C             partial reorganization of aswip code including modification
C             to correct an issue that occasionally resulted in Rmax=0
C             used as input to NWS19 subroutine in adcirc.
C  v50.36   - 09/13/2011 - tcm - Added tabs and comma delimiters to the parse
C             subroutine to handle issues with parsing station locations
C             in the fort.15 file.  Utlitmately the station locations should
C             be read in as two real numbers followed by a ! and then a comment
C             line and the adcirc documention updated to reflect a standard
C             way of inserting comments.
C  v50.35   - 09/09/2011 - jgf - Fixed issue in adcprep when using prep13
C             option: allocated and populated imap_nod_gl in relocalize().
C  v50.34   - 09/01/2011 - jgf - Merged the following changes from r689:HEAD
C             of the nws_namelist branch:
C     v50.31.03  08/25/2011 - jgf - Checked for stationary storm when finding
C              storm direction in wind.F; added NWS value to read_input.F.
C     v50.31.02  08/19/2011 - jgf - Made different capitalizations of drag
C              law formulations valid; removed stop statement that was inserted
C              after reading of nws line for debugging.
C     v50.31.01  08/18/2011 - jgf - Reorganized reading of NWS line so that it
C              could be read as a namelist and still be backward compatible with
C              existing fort.15 files; made WindDragLimit and DragLawString
C              optionally selectable from the fort.15 with default value of
C              0.0035 for WindDragLimit and default values of Garratt and
C              IceCube (depending on whether ice drag is active) for
C              DragLawString; added storm center recording to NWS=8, 9, and 19
C              to make them compatible with Powell (sector-based) drag
C              formulation; made wind drag limit in SWAN the same as in ADCIRC.
C  v50.33   - 09/01/2011 - jgf - Fixed compiler flags for garnet at ERDC.
C  v50.32   - 09/01/2011 - jgf - Fixed compiler flags for diamond at ERDC.
C  v50.31   - 08/16/2011 - jgf - Added call to ApplyCanopyCoefficient
C             subroutine to all wind subroutines.
C  v50.30   - 08/03/2011 - jgf - Modified NWS19 so that it won't try to
C             interpolate beyond the end of the fort.22 file if the run
C             ends right at the last line of the file.
C  v50.29   - 08/01/2011 - jgf - Modified aswip so that if a storm becomes
C             stationary, its direction will be persisted and speed reduced
C             to 1 ... feeding 0 and 0 to ADCIRC NWS19 caused NaNs.
C  v50.28   - 07/26/2011 - jgf - Removed vestigial DO/ENDDO surrounding call to
C             SWMAIN in couple2swan, purely to enhance readability.
C  v50.27   - 07/12/2011 - tcm - Added the option to use dedicated writer
C             processors for full in addition to compact formats for ascii
C             global field files (fort 63,64,73,74,93,164 (cswan) ).  Renamed
C             local variable "default" to "default_val".
C  v50.26   - 07/06/2011 - jcd - Changed owiwind.F so that blank snaps are
C             inserted when the end of file is reached.
C  v50.25   - 06/28/2011 - jgf - Fixed debug flags for Jade; removed -Mextend
C             from macros.inc for Jade; corrected typo in hstime.F.
C  v50.24   - 06/27/2011 - jgf - Made NWS=19 look in the fulldomain directory
C             for fort.22 in parallel rather than the subdomain directories.
C  v50.23   - 06/27/2011 - jgf - Added NWS=319 to list of allowable NWS values
C             in read_input.F; fixed initialization in hstime.F; fixed debug
C             flags for intel compiler.
C  v50.22   - 06/27/2011 - jcd - Set the unit number swiLUN = 26 in prep.F
C  v50.21   - 06/10/2011 - tcm - Changed I8 to I12 format specifications in
C             the fort.18 files (prep.F, messenger.F, SwanReadfort18.ftn90)
C  v50.20   - 06/10/2011 - jgf - Fixed an issue in prepUnSWAN that was causing
C             the opening of the fulldomain fort.26 file to fail.
C  v50.19   - 06/07/2011 - tcm - Fixed a bug in the makefile for the
C             library build of parallel adcirc.
C  v50.18   - 06/08/2011 - rjw - Added CASE(bottom_roughness_length) to prep.F
C  v50.17   - 06/07/2011 - jgf - Added -freal-loops to g95 compiler flags
C             for padcswan; changed -module to $(IMODS) in swan compilation
C             rules in makefile; added specific format lengths for I
C             to placate g95 compiler.
C  v50.16   - 06/06/2011 - jgf - Fixed issue with netcdf hotstart where the
C             time (in seconds) was being read from the netcdf hotstart file
C             but not communicated back to hstart.F for use in initialization;
C             added enhanced logging and i/o error checking to owiwind.F, as
C             well as creating a subroutine to remove some repeated code for
C             i/o; widened the format for NSPOOLG* in globalio.F from
C             I5 to I8 and widened format for NHSINC from I5 to just I in
C             read_input.F to avoid format overflows.
C  v50.15   - 05/31/2011 - jgf - Modified prepUnSWAN routine in adcprep so that
C             it loads the name of the swan command file from the swaninit
C             file; it will then copy that command file to the subdomain
C             directories without changing the name.
C  v50.14   - 05/27/2011 - RJW - Fixed 3D hotstart bug in the turb subroutine,
C             that resulted in variables not being initialized upon return to
C             the program from hotstart.(vsmy.F)
C             Explicitly allow subroutines to access variables using the "USE"
C             statement.  This was needed to ensure MYPROC, MNPROC, SZ were
C             being properly used in subroutines (write_output.F,globalio.F,
C             global.F).
C  v50.13   - 05/27/2011 - jgf - Fixed issue caused by requesting netcdf
C             station output while specifying zero stations.
C  v50.12   - 05/26/2011 - tcm - Fixed bug that caused output to overflow the
C             formats on some compilers, for example, an output of
C             -0.1234567891E-123 would be written as -0.1234567891-123, when
C             using the Ew.d format.  The FORTRAN 90 standard is not clear
C             on how the exp part is displayed when using Ew.d so we area
C             switching to the Ew.dEe format.  Changes were made in
C             globalio.F, write_output.F and timestep.F for these.  Also,
C             made the output formats consistent between globalio.F and
C             write_output.F.
C  v50.11   - 05/26/2011 - jgf - Added capability to adcprep to prep just
C             the fort.88 file; add capability to just prep the fort.20 file;
C             added command line argument processing to adcprep to bypass
C             menu driven interface; the menu driven interface is still
C             the default, and is triggered if there aren't any command
C             line options.
C  v50.10   - 05/25/2011 - jgf - Merged rivers50 branch r643:660, which
C             contained the changes for initialization of rivers above msl.
C    v50.05.05  05/23/2011 - jgf - Removed code that had been commented out b/c
C               it was no longer being used.
C    v50.05.04  05/13/2011 - jgf - Updated hstime utility to remove hardcoding
C               of the hotstart file name; added support for netcdf hotstart
C               files.
C    v50.05.03  05/13/2011 - jgf - Commented out more leftover debugging
C               statements.
C    v50.05.02  05/12/2011 - jgf - Commented out leftover debugging statements.
C    v50.05.01  05/11/2011 - jgf - Merged and resolved conflicts for the
C             following changes from the rivers_init branch:
C    ------------------------------------------------------------------------
C    r484 | kdresback | 2010-06-28 12:19:14 -0400 (Mon, 28 Jun 2010) | 1 line
C    Changed paths:
C      M /adcirc/branches/rivers_init/prep/prep.F
C      M /adcirc/branches/rivers_init/prep/presizes.F
C      M /adcirc/branches/rivers_init/prep/read_global.F
C      M /adcirc/branches/rivers_init/src/hstart.F
C      M /adcirc/branches/rivers_init/src/read_input.F
C      M /adcirc/branches/rivers_init/src/timestep.F
C      Updates for the river initialization - adding the option for the river
C      boundary condition to not use time==0 as the start time.
C    ------------------------------------------------------------------------
C    r460 | kdresback | 2010-06-07 18:29:40 -0400 (Mon, 07 Jun 2010) | 1 line
C    Changed paths:
C      M /adcirc/branches/rivers_init/prep/pre_global.F
C      Updates for river initialization when rivers are above MSL
C    ------------------------------------------------------------------------
C    r454 | kdresback | 2010-06-02 15:52:24 -0400 (Wed, 02 Jun 2010) | 1 line
C    Changed paths:
C      M /adcirc/branches/rivers_init/prep/adcprep.F
C      M /adcirc/branches/rivers_init/prep/prep.F
C      M /adcirc/branches/rivers_init/prep/presizes.F
C      M /adcirc/branches/rivers_init/prep/read_global.F
C      M /adcirc/branches/rivers_init/src/cstart.F
C      M /adcirc/branches/rivers_init/src/global.F
C      M /adcirc/branches/rivers_init/src/global_3dvs.F
C      M /adcirc/branches/rivers_init/src/hstart.F
C      M /adcirc/branches/rivers_init/src/init_cond.F
C      M /adcirc/branches/rivers_init/src/read_input.F
C      M /adcirc/branches/rivers_init/src/transport.F
C      M /adcirc/branches/rivers_init/src/vsmy.F
C      Updates for river initialization when rivers are above MSL
C    ------------------------------------------------------------------------
C    r453 | kdresback | 2010-06-02 14:10:49 -0400 (Wed, 02 Jun 2010) | 1 line
C    Changed paths:
C      A /adcirc/branches/rivers_init (from /adcirc/trunk:447)
C      Creating a branch for river initialization when rivers are above MSL
C    ------------------------------------------------------------------------
C  v50.09   - 05/18/2011 - jcd - Incorporated Powell's sector-based wind drag
C             scheme, which is controlled by the -DPOWELL compiler flag.
C             Also corrected a few bugs.
C  v50.08   - 05/18/2011 - jgf - Changed node to rawNode in subroutine that
C             reads in min/max files ... not sure if this caused any issues
C             but the Intel compiler complained; fixed format statement for
C             log message reporting the value of IHOT.
C  v50.07   - 05/18/2011 - jgf - Added missing subroutine call to
C             readNetCDFHotstart3D in hstart.F; cleaned up typos in the
C             netcdf 3D hotstart reading subroutine.
C  v50.06   - 05/16/2011 - jgf - Fixed issue where 3D netcdf hotstart file
C             was requested but no other fulldomain 3D netcdf output was
C             requested (had to initialize num_v_nodes to NFEN in netcdfio.F).
C             Added some check_err calls to netcdfio.F  where they were
C             previously missing. Fixed an issue with turbulence length
C             scale output in netcdf (nodal_data_id was referenced rather than
C             v_nodal_data_id). Fixed format issue with NHSTAR in read_input.F.
C             Initialized E0XGrad2A and E0YGrad2A in timestep.F to zero
C             to appease Intel compiler.
C  v50.05   - 05/10/2011 - tcm - added support in adcprep for compilation
C             on Window's based PC.  Fixed a couple of stray bugs related
C             to character lengths in prep.F and added the NOFF variable to
C             hs_2_ascii.F.
C  v50.04   - 05/10/2011 - jgf - Fixed stray metadata in netcdf file due to
C             nvelnc not being defined correctly.
C  v50.03   - 05/09/2011 - jgf - Applied fixes to timestep.F and ipackv.F from
C             st to fix issue when the NeiMax node is on a boundary.
C  v50.02   - 05/06/2011 - tcm - Added a list of pre-processor definitions to
C             use in order to build Metis on a Window's based PC.
C  v50.01   - 05/06/2011 - tcm - Added a readme file about the changes to metis.
C  v50.00   - 05/06/2011 - tcm - Made a couple of modifications to the packages
C             called by metis to get it to work on a Windows based machine.
C             Added the gettimeofday.h module which is not available on a
C             Windows machine.  A person wanting to compile metis for a Windows
C             based PC would select all the routines in the Lib folder and
C             would also include the gettimeofday.h routine in the Windows
C             folder and would use a -WINDOWS pre-processor flag in addition to
C             the other pre-processor flags in the normal makefile.
C             NOTE:  Only use the gettimeofday.h for Windows machines,
C             Linux/Unix machines typically already have it installed.
C             These changes to metis do not change the usual building or
C             functionality of metis on Linux/Unix machines.
C  v49.82   - 04/26/2011 - jgf - Fixed initialization issue in netcdf when
C             netcdf hotstart file was requested in 3D but no other fulldomain
C             netcdf output was requested; removed warning message about
C             ITITER from adcprep; initialized previously uninitialized
C             variable in station parsing routine.
C  v49.81   - 04/25/2011 - rjw - Fixed bug in ./swan/Makefile now it will be able
C             to read the variables from global_3dv
C  v49.80   - 04/25/2011 - jgf - Removed stray debug msg from read_input.F.
C  v49.79   - 04/21/2011 - rjw - Fixed bug in ./prep/read_global.F for IEVC=51
C  v49.78   - 04/21/2011 - rjw - Modified files in ./prep to allow for IEVC=51
C  v49.77   - 04/17/2011 - jgf - Removed erroneous second read of WTIMINC line
C             for NWS19.
C  v49.76   - 04/14/2011 - tcm - Cleaned up more reserved words in post.F
C               and adcpost.F  (INDEX,LOC,SCALE).
C  v49.75   - 04/14/2011 - tcm - Commented out debugging flags for cray_xt4,
C               Cleaned up us of FORTRAN RESERVED WORDS (TIME,INDEX,SCALE,
C                ISIGN,LOC,IIDIM). Changed to (TIME --> TIMELOC,INDEX --> INDX,
C                ISIGN --> ISGN, LOC --> LOC2, SCALE --> SCALEVAL,
C                IIDIM --> IDIM2).  The use of reserved words can lead to
C                unintended consequences and conflicts.
C  v49.74   - 04/13/2011 - tcm - Fixed bug in call to FEGRIDS in read_input.F
C  v49.73   - 04/11/2011 - jgf - Merged vortexowi branch (r421:r536) to trunk,
C             bringing in the changes from v49.08.01 to v49.08.04 to create
C             NWS=29 option for blending vortex winds (NWS19) with
C             OWI winds (NWS12).
C  v49.0804   - 12/20/10 - jgf - Made allowance for NWS29 in read_global.F,
C               prep.F and adcprep.F; cstart.F now looks in GBLINPUTDIR for
C               the fort.22 for NWS29; corrected missing parenthesis in wind.F;
C               added reading of WTIMINC line for NWS29 in presizes.F; added
C               blend factor for blending vortex winds and OWI winds.
C  v49.0803   - 04/16/10 - jgf - Added code to NWS19 to properly export the
C               vortex center coordinates and current Rmax.
C  v49.0802   - 04/08/10 - jgf - Fixed array out of bounds in log message
C               in timestep.F. Eliminated duplicate open of fort.22 in
C               NWS19. Fixed file name to open for NWS29 in hstart.F.
C  v49.0801   - 03/19/10 - jgf - Added code for new NWS29 option that embeds an
C               asymmetric vortex (NWS19) into an OWI formatted (NWS12)
C               basin scale met field (derived from NAM). Also separated
C               the application of wind reduction factors from the application
C               of the canopy coefficient, and made several key variables
C               from NWS19 available at the wind module level for use by
C               other wind-related subroutines.
C  v49.72   - 04/11/2011 - jgf - Removed init_cond.F and related comments
C             from makefile since they are no longer needed.
C  v49.71   - 04/11/2011 - kmd - Updated the initial condition file for
C             baroclinic to be evaluated in parallel similar to the binary
C             hotstart file; updated station information for 3D.
C  v49.70   - 04/11/2011 - rjw - Updated KSlip formulation fro KSlip=2, 3
C             Added subroutine "Apply3DBottomFriction", to compute TK for 3D
C             simulations.  This fixes Tau0 bug where time varying tau0 is
C             computed as function of TK.  Added a call to "CalculateTimeVaryingTau0"
C             in timestep.F for C3D.
C  v49.69   - 04/07/2011 - jgf - Fixed typo (missing .) in netcdfio.F.
C  v49.68   - 04/07/2011 - tcm - Enlarged the field width (w) on several
C             format specifications for Ew.d fields.  Note that w should be
C             set greater than or equal to 7 more than d according to the
C             FORTRAN 90 standards.  For Fw.d, w should be greater than or
C             equal to 3 more than d.
C  v49.67   - 04/06/2011 - tcm - making additions for compiling adcirc on a
C             Windows based PC.  In vortex.F resolved a conflict with the
C             variable "error" used as an integer in global.F and as a
C             real in constants.F.
C  v49.66   - 04/06/2011 - tcm - merged ice coverage into the trunk from
C             branch versions v49.64.01 through v49.64.01.
C  v49.64.02- 04/01/2011 - tcm - Fixed bug in Harmonic Analysis output when
C             output was called to be generated but the model never reached
C             the time interval between THAS and THAF during which the data
C             is to be created.  Added a test for NWS and NCICE
C             compatability.  Included Ice support in the init_cond.F file.
C  v49.64.01- 03/31/2011 - tcm - Added ice field data for use in modifying
C             wind drag coefficients (Enabled by setting NWS >= 1000
C             uses files (fort.25 all ice types, fort.225 and fort.227
C             only for ice types (NCICE) = 12 ). Added writing of global and
C             station ice concentration files (fort.93 and fort.91).
C             Added a cubic drag coefficient function based on ice
C             concentration levels and the typical Garratt wind speed
C             computations.  Ice concentration files are only written
C             when met. files are written.
C  v49.65   - 04/06/2011 - tcm - Fixed typo bug in transport.F.
C  v49.64   - 04/06/2011 - rjw - Modified 3D code to allow user to select
C             standard Mellor-Yamada 2.5 closure (IVEC=50), or modified
C             Mellor-Yamada 2.5 with enhanced mixing (IVEC=51). (Note: for
C             IEVC=51, the alpha parameter for Z0S is read in from EVCON)
C             Added ability for 3D code use Manning's n values to compute
C             bottom roughness values based on Casey's work taken from
C             Bretschneider(1986).  Modified 3D wetting scheme to simplify the
C             calculation of BSX/Y at newly wetted node.
C  v49.63   - 03/29/2011 - rjw - Removed debugging write stament from vsmy.F
C  v49.62   - 03/28/2011 - rjw - Modified the 3D code to include a hard wired
C             enhancement to the Mellor-Yamada 2.5 Turbulence Closure scheme
C             to account for increased mixing at the surface based on magnitude
C             of surface stress.  Additionally, added nodal attribute,
C             "bottom_roughness_length" to be used in the 3D simulations.
C             included a couple of bug fixes found when implementing the
C             IM=611112 run type, and added back in the calculation of BSX and BSY
C             to the 3D wetting code in timestep.f.
C  v49.61   - 03/17/2011 - jgf - Corrected compiler flags for IBM p6 and
C             removed extra -DNETCDF for IBM p6 platform; corrected double
C             deallocation issue when writing ascii fort.44 file in parallel.
C  v49.60   - 03/09/2011 - jgf - Removed IStart3D variable from global_3dvs.F
C             to fix initialization issues with evtot; explicitly initialized
C             wsx1, wsy1, wsx2, wsy2 to avoid NaNs in vssol upon hotstart
C             of 3D runs in serial; added error handling to code that skips
C             down to the proper location in 3D output files; added flags for
C             gfortran compiler to take advantage of floating point exception
C             trapping for debugging.
C  v49.59   - 03/09/2011 - tcm - Changed "COMPILER" back to "compiler" in the
C             cmplrflags.mk file.
C  v49.58   - 03/08/2011 - jgf - Fixed NeiTabEle handling to account for its
C             semi-sparseness in subroutine bpg3D and in transport.F.
C  v49.57   - 03/08/2011 - jgf - Added private declaration to subroutines
C             in netcdf module to fix compilation on intel v9.1 compiler
C             (on kittyhawk); also tweaked intel compiler flag setup.
C  v49.56   - 03/03/2011 - jgf - Merged in netcdf3D branch revisions 562:597
C             to trunk (branch version numbers v49.48.01 through v49.48.09).
C  v49.48.09- 03/03/2011 - jgf - Added compiler flags for IBM P6; added
C             capability to activate comprehensive trace logging with a single
C             compiler flag; changed trace logging so that all processors
C             write trace messages to local fort.16, which makes debugging
C             easier during parallel execution; added conditions
C             to netcdf memory deallocation for harmonic analysis hotstart;
C             fixed issue in binary hotstart in parallel by mapping to
C             NNODECODE instead of NODECODE; removed IHOT as argument to
C             initHarmonicParameters and eliminated duplicate calculations
C             of nfreq in harm.F; added doc directory with rough draft
C             for an ADCIRC Developers Guide.
C  v49.48.08- 02/25/2011 - jgf - Fixed record accounting in netcdf module
C             so that appends place data in proper position within the netcdf
C             data file, rather than always appending at the end; fixed issue
C             where 'it' was missing from call to hstart3D_out; added
C             call to prepnetcdf() to PREP15 option in adcprep.
C  v49.48.07- 02/24/2011 - jgf - Brought over vparr's adccmp.F from ERDC_esmf
C             branch for use in testing; added experimental new autotest code
C             for use in automated testing.
C  v49.48.06- 02/24/2011 - jgf - Fixed initialization issues for 3D netcdf
C             output files and 3D netcdf hotstart files in parallel; added
C             code to adcprep to capture 3D station names; added initialization
C             of 3D station coordinates to adcprep and adcirc/write_output.F.
C  v49.48.05- 02/23/2011 - jgf - Added dedicated arrays for storing coordinates
C             of 3D stations for use in 3D netcdf station output files; fixed
C             several typos and redundancy issues in netcdf module that
C             were interfering with 3D netcdf output.
C  v49.48.04- 02/23/2011 - jgf - Added new code to adcprep to set up
C             3D netcdf output files for parallel runs and to set up netcdf
C             hotstart file file with 3D data for parallel runs; further
C             modularized mapping of full domain data to subdomains on the
C             fly during hotstart; modularized setting of output timing
C             during 3D hotstart.
C  v49.48.03- 02/11/2011 - jgf - Added new subroutines to netcdf module for
C             reading 3D data from netcdf hotstart file and for mapping
C             this 3D data from fulldomain to subdomains, combining reading
C             and mapping for each variable to conserve memory to the
C             greatest extent possible.
C  v49.48.02- 02/09/2011 - jgf - Added new subroutines to netcdf module for
C             allocating memory for 3D data in netcdf hotstart files;
C             initializing 3D data in netcdf hotstart files; and writing 3D
C             data to netcdf hotstart files; rearranged HSTART3D_OUT in
C             vsmy.F to handle both binary and netcdf output and to globalize
C             each 3D variable all at once rather than in horizontal layers
C             (although the layered output format is preserved for binary
C             output) and to allocate the memory for each 3D variable just
C             before it is needed and deallocating memory immediately after
C             its use is complete; the option of writing either subdomain or
C             fulldomain hotstart files was preserved for binary output;
C             the writing of data to binary hotstart files was also modularized
C             with the addition of a new subroutine.
C  v49.48.01- 02/03/2011 - jgf - Added new subroutines to netcdf module for
C             initializing and writing 3D data to netcdf output files;
C             modularized code for reading 3D stations in read_input.F and
C             added capability for reading 3D station names.
C  v49.55    - 03/03/2011 - tcm - Merged the last of the changes from the
C             ERDC tightly coupled with STWAVE branch.
C  v49.54.04 -03/03/2011 - tcm - Bring updates from the trunk getting
C  v49.54.03 -02/24/2011 - tcm - Added the option to make parallel adcirc
C             as a library back to the makefile.  This library option is
C             needed for the external coupler.
C  v49.54.02- 02/34/2011 - jgf - Fixed minor compilation issues in adcprep.F
C             and writer.F.
C  v49.54.01- 02/22/2011 - tcm - Modifications to adcirc.F and writer.F to
C             allow for coupling with STWAVE via external couplers.
C  v49.54   - 02/22/2011 - tcm - Fixed a bug with hot start writers related
C             to the NOFF varible in writer.F.  Cleaned up use of reserved
C             (dim, index, and scale) in the other files changed in this
C             commit.  (dim --> iidim, index --> iindx, scale --> scaleval)
C  v49.53   - 02/18/2011 - jgf - Cleaned up logging for netcdf in adcprep
C             and fixed minor bug with v-velocity metadata in netcdf hotstart
C             files.
C  v49.52   - 02/17/2011 - jgf - Declared several subroutines in netcdf
C             module as private to enable compilation on Intel v9 compiler
C             on kittyhawk at RENCI.
C  v49.51   - 02/09/2011 - tcm - Merging fast search algorithms from the
C             merge_erdc_to49 branch.  Incorporated the new search tree
C             searches into the prep routines and read_input.F. Updated the
C             readstations subroutine to accept a description variable for
C             allowing for the inputing of the correct station type.  Added
C             the writer_pause subroutine but not using it currently.
C  v49.50   - 02/09/2011 - jgf - Fixed an issue in hstart.F where harmonic
C             analysis data were checked, even though harmonic analysis was
C             not specified as part of the simulation.
C  v49.49   - 02/02/2011 - jgf - Rearranged memory allocation for 3D globalio
C             to minimize the total amount of RAM required at any given time
C             for writing 3D output files.
C  v49.48   - 01/25/2011 - tcm - Added support for NRS=4 tight coupling with
C             STWAVE in timestep.F.
C  v49.47   - 01/25/2011 - rjw - Changed format statement 2393 to 2395 on
C              cmassey STWAVE routine in Updated read_input.F
C  v49.46   - 01/25/2011 - tcm - Updated read_input.F and adcprep routines
C             to allow for NRS=4 which is tightly coupled ADCIRC+STWAVE.
C  v49.45   - 01/24/2011 - jgf - Merged harmio branch r527:556 to trunk.
C  v49.43.19  01/23/2011 - jgf - Finished code to globalize all 3D output
C             files as well as code to write fulldomain 3D files in ascii
C             format in serial and in parallel; fixed 3D output file
C             initializations in cstart.F and hstart.F so they work in
C             parallel; corrected issues with reading 3D data from binary
C             hotstart file in parallel and mapping to subdomains on the fly.
C  v49.43.18  01/21/2011 - jgf - Added code to prep18 to write 3D station
C             mappings into fort.18 for 3D runs; added code to messenger.F
C             and global_3dvs.F to read 3D station mappings from fort.18 for
C             3D runs; rearranged code to write out 3D turbulence stations
C             (fort.43) so that it is organized by output format; also added
C             code to use the new mappings to globalize and write out
C             the resulting fulldomain fort.43 files in ascii.
C  v49.43.17  01/19/2011 - jgf - Added globalization of 3D density (fort.44)
C             and 3D velocity (fort.45) to write_output.F; also added code
C             to write out these fulldomain 3D data in ascii; rearranged
C             existing 3D output code so that it is organized according to
C             output format.
C  v49.43.16  01/18/2011 - jgf - Added new globalio subroutines packNPbyM
C             and unpackNPbyM for globalizing 3D output; made changes to
C             enable globalization of fulldomain turbulence output file
C             (fort.46).
C  v49.43.15  01/18/2011 - jgf - Corrected several issues revealed by internal
C             testing of reading and writing harmonic analysis data in netcdf
C             format.
C  v49.43.14  01/15/2011 - jgf - Finished code for reading means and variance
C             calculation data for harmonic analysis from netcdf hotstart file;
C             code includes capability to read data in serial or in parallel
C             (mapping to subdomains in parallel on the fly).
C  v49.43.13  01/14/2011 - jgf - Finished code for reading harmonic analysis
C             data from netcdf hotstart file to continue harmonic analysis
C             during the hotstarted run; code includes capability to read
C             in h.a. data in in parallel as well as serial.
C  v49.43.12  01/12/2011 - jgf - Made adcprep reliant on harm.F for harmonic
C             variables; made adcprep reliant on global.F for station names;
C             completed changes to netcdfio.F and adcprep to allow harmonic
C             data to be written to hotstart files in parallel.
C  v49.43.11  01/10/2011 - jgf - Added code to netcdf module to reduce the
C             size of the 2D matrix that holds the boundary type arrays
C             so that they are dimensioned by (nope,max(nvdll)) and
C             (nbou,max(nvell)).
C  v49.43.10  01/09/2011 - jgf - Fixed initialization issues that were
C             preventing log and console messages from showing up during
C             adcprep; straightened out initialization of mesh boundaries
C             in netcdf to fix seg fault issue during adcprep of sl16 grid;
C             added a few missing calls to check_err in netcdf.
C  v49.43.09  12/24/2010 - jgf - Fixed issue in prep.F/prepnetcdf() that
C             caused compilation with netcdf to fail.
C  v49.43.08  12/23/2010 - jgf - Fixed issue where netcdf harmonic analysis
C             hotstart data were not being initialized in write_output.F;
C             added icha parameter to netcdf harmonic analysis hotstart file.
C  v49.43.07  12/22/2010 - jgf - Fixed numerous small issues in netcdf
C             module that were preventing it from writing harmonic analysis
C             load vectors to netcdf hotstart files.
C  v49.43.06  12/21/2010 - jgf - Corrected issue with initialization of netcdf
C             mesh data for hotstart file; added check for NHSTAR=3 before
C             initializing netcdf hotstart file; forced a rewrite of the
C             netcdf hotstart file each time it is called to write; moved
C             initialization of netcdf element data to initNetCDFCoord().
C  v49.43.05  12/21/2010 - jgf - Added subroutines to netcdf module to
C             free memory that was allocated for writing netcdf data; added
C             subroutine calls to prep.F to free memory from netcdf module
C             once the data had been written; added initialization to
C             screenUnit i/o parameter in adcprep which corrected a memory
C             fault issue that occurred when using netcdf; moved call to
C             prepNetCDF() to just before the call to DECOMP() in adcprep.
C  v49.43.04  12/21/2010 - jgf - Added code to reduce memory use of netcdf
C             module by only retaining a single copy of the mesh and its
C             boundaries in memory and reusing it for each netcdf output file;
C             added new data types for storing harmonic station data
C             as well as fulldomain harmonic data to netcdf files; added
C             new variables, metadata, and subroutines for storing harmonic
C             analysis load vectors to netcdf output files; added subroutine
C             calls for tracing progress of subroutine calls to netcdf
C             module.
C  v49.43.03  12/09/2010 - jgf - Added code to write out means and variance
C             calculation results for harmonic analysis to full domain output
C             file in parallel; fixed issues in globalio for new packMbyNP;
C             fixed issues in output data descriptors for harmonic data;
C             added general subroutines for writing log and screen messages
C             to aid in debugging; fixed descrepancies between reading
C             and writing of binary hotstart file.
C  v49.43.02  12/08/2010 - jgf - Added new subroutines to write_output.F
C             to write harmonic analysis data to full domain output files
C             in parallel; finished separating harmonic analysis solution
C             from writing harmonic analysis output files.
C  v49.43.01  11/30/2010 - jgf - Preserved the reading of local hotstart
C             files as an option, but made the reading of fulldomain hotstart
C             files the default. Made changes to enable reading of harmonic
C             analysis hotstart data and 3D hotstart data from fulldomain
C             hotstart file in parallel by mapping the nodal data from
C             fulldomain to subdomain on each PE on the fly. Made changes
C             to globalio by adding new subroutines that can globalize harmonic
C             analysis data and 3D data. Reorganized harmonic analysis module
C             in the following ways: (1) moved variables from global.F to
C             harm.F if they pertained only to harmonic analysis; (2) separated
C             reading from initialization when harmonic analysis hotstart
C             file is read; (3) separated solution from writing when harmonic
C             analysis data is written to output; (4) moved code that updates
C             load vectors from timestep.F to harm.F; (5) moved code that
C             calls various subroutines to solve harmonic analysis at end of
C             a run from adcirc.F to harm.F. Added new code to
C             write_output.F to write globalized harmonic analysis data to
C             fulldomain hotstart file. Added subroutines to global.F
C             to streamline reading of hotstart arrays from binary files.
C             Removed subroutines from globalio.F that were not called
C             from anywhere. Made canopy coefficient default to 1 in nodal
C             attributes. Added fix for case where NRAMP is 1 but DRAMP is 0.
C             Updated hot2asc to read harmonic analysis data.
C  v49.44   - 01/19/2011 - tcm -- For owiwind.F changed the way the header
C             line of the file was read so that if changes in the text
C             portion occur in a file, these would not cause read failures.
C             This line is an 80 character line with two 10 digit date
C             fields which are the only data used off that line.  The first
C             date field is contained in positions 56-65, and the second
C             date field is in positions 71-80.
C  v49.43   - 11/05/2010 - jgf - Removed netcdf remnants from hotlocalize
C             in prep.F and corrected a seg fault bug in hotlocalize caused
C             by six-digit IM values.
C  v49.42   - 11/05/2010 - jgf - Added compiler flags and library flags for
C             compiling with netcdf support on sapphire, jade, and diamond
C             at ERDC.
C  v49.41   - 11/04/2010 - jgf - Fixes for reading and writing netcdf hotstart
C             files in serial and in parallel; fix to allocate memory for
C             boundary arrays more efficiently in netcdf module.
C  v49.40   - 10/25/2010 - jgf - Fixes for writing netcdf metadata to trim
C             off trailing spaces when writing netcdf files from adcirc or
C             adcprep, as well as when writing subdomain fort.15 files from
C             adcprep.
C  v49.39   - 10/22/2010 - jgf - Updates to netcdf module for hotstarting,
C             including alternation of netcdf 67 and 68 files; additions to
C             support fulldomain integer arrays for nodecode and noff; updating
C             it and im in hotstart file; fix for ordering write in prep15
C             (write metadata after 3D data in fort.15).
C  v49.38   - 10/06/2010 - jgf - Significant rearrangement of hstart.F to
C             consolidate the reading of the binary hotstart file all in
C             one section of the hotstart() subroutine; this was necessary
C             b/c the subroutine mingled the reading of the file with
C             initialization and logging. In order to have a choice of
C             file formats, the reading of the hotstart file must be a
C             separate step from initialization, logging and other tasks.
C             A call to read a netcdf hotstart file if IHOT.eq.367 or
C             IHOT.eq.368 was added.
C  v49.37   - 10/05/2010 - jgf - Removed utility program called init b/c it
C             was not being used; made compilation of prep.F and read_global.F
C             dependent on adcirc modules to facilitate sharing of common
C             code (this fixes an issue with compiling adcprep w/o netcdf);
C             finished reorganization and redevelopment of netcdf hotstart
C             code in netcdf module (netcdf hotstart testing still required).
C  v49.36   - 10/04/2010 - rjw - 1) Modifications to 3D turbulence: addition of
C             enhanced surface mixed layer.  Variety of chioces hardcoded in
C             vsmy.f.  Surface roughness in fort.15 is now the minimum
C             surface roughness.  2) Modifications to 3D wetting and drying
C             in timestep.f, modified wetting velocity formulation.
C             3) Added ability to run 3D in GWCE explicit mode.  GWCE is
C             solved outside 3D routine, IM type 611122 enables all of the
C             same GWCE explicit settings as 111122 except indicates to run
C             3D instead of 2D.  To implement also change GWCE time weighting
C             factors to 0.0 1.0 0.0 and change solver flag to 0.
C  v49.35   - 10/04/2010 - jgf - Fixes for netcdf: Fixed dimensioning of
C             NVELL during adcprep (was causing segfault on ec95d.grd among
C             others); streamlined netcdf module by removing redundant
C             dimensioning and ids for multidimensional nodal and station
C             data; fixed error handling and added graceful exit when
C             netcdf file from adcprep is not present when padcirc starts;
C             added and significantly rearranged code for writing netcdf
C             hotstart files.
C  v49.34   - 09/28/2010 - jgf - Fixes for netcdf: added call to update meta
C             data when hotstarting; eliminated ncmode variable, the
C             subroutine that sets it, and all references to it throughout
C             adcirc and adcprep.
C  v49.33   - 09/27/2010 - jgf - Fixes for netcdf: boundary data now written
C             correctly by adcprep now; read_global.F now uses the parse()
C             and a2f() subroutines from global.F instead of repeating them
C             in read_global.F ... this allows adcprep to pick up the fixes
C             that were made in global.F for parsing station names and also cuts
C             down on the incidence of cut-and-pasted code in ADCIRC.
C  v49.32   - 09/23/2010 - jgf - Fixes for netcdf: SLAM0 and SFEA0
C             as well as mesh boundaries. Now works in serial. Parallel
C             should also work, but testing is required to confirm.
C  v49.31   - 09/23/2010 - jgf - Changes for writing netcdf in parallel:
C             added if/then to use array_g etc instead of array when writing
C             station data in parallel run; added a subroutine that is called
C             by adcprep to populate global variables so that they may be
C             written into the new initialized netcdf output file; moved
C             fort.15 parsing subroutines to global.F so that they can be
C             called from different places; fixed miscellaneous bug in NWS19
C             in read_input.F by reading RSTIMINC when NRS.GE.1 instead of
C             just NRS.eq.1; fixed memory allocation for station names
C             in prep/read_global.F; moved the hotstart-relevant variables
C             from local declarations in prep.F to declarations in
C             pre_global.F; created new subroutine in prep.F to initialize
C             netcdf output files for parallel runs where they will be
C             required.
C  v49.30   - 09/16/2010 - jgf - Working on aswip: added capability to set
C             the Rmax to a constant or to a percentage of the calculated
C             value based on isotach data (if any) via command line options;
C             added capability to set method of selecting quadrant Rmax
C             via command line option; removed persistence of average Rmax
C             from hindcast into the forecast; corrected issue where
C             isotach speed can sometimes be zero in source data.
C  v49.29   - 08/25/2010 - jgf - Cleaned up parsing of metadata and station
C             names such that the netcdf module no longer assumes that a "!"
C             is present on each line, but is able to use it as a metadata
C             terminator if it is present.
C  v49.28   - 08/24/2010 - jgf - Declared variables in netcdf module as
C             private if the underlying type declaration was private.
C  v49.27   - 08/15/2010 - jgf - Fixes to storage of meta data, lat/lon
C             coordinates for stations and mesh nodes, and station names in
C             netcdf.
C  v49.26   - 08/15/2010 - jgf - Fixes to netcdf so it is now capable of
C             writing all data to the eight main output files in serial.
C ------------------------------------------------------------------------
C   r501 | stanaka | 2010-08-11 10:20:22 -0400 (Wed, 11 Aug 2010) | 1 line
C             Changed paths:
C             M /adcirc/trunk/src/timestep.F
C             Implemented Clint's zombie killing process for the error
C             elevation node. If there is a error elevation node, MSG_FINI is
C             called to finalize MPI environment.
C             (jgf: let's call this v49.25.)
C  v49.24   - 08/10/2010 - jgf - Fixed issues with the construction of arrays
C             for writing in netcdf files. Data in the reorganized module were
C             not being copied correctly from ADCIRC arrays to netcdf arrays.
C  v49.23   - 08/06/2010 - jgf - Fixes to reorganized NetCDF module.
C             Files now initialize properly, but there are still issues with
C             writing actual data.
C  v49.22   - 08/06/2010 - jgf - Commented out DEBUG statements in sizes.F
C             (related to command line option processing) and read_input.F
C             (related to opening files fort.14, fort.15, and fort.16).
C  v49.21   - jgf - Merged the subdomainhotstart branch into the trunk,
C             which included a massive reorganization of the netcdf module
C             for greater flexibility. The subdomainhotstart branch included
C             the following changes:
C  v49.17.04  08/04/2010 - jgf - Fourth set of changes for netcdf.
C  v49.17.03  07/29/2010 - jgf - Third set of changes for netcdf.
C  v49.17.03  07/27/2010 - jgf - Second set of changes for netcdf.
C  v49.17.02  07/23/2010 - jgf - First set of changes for working through
C             netcdf issues.
C  v49.17.01  07/14/2010 - jgf - Added an extra command line option -S to
C             write subdomain hotstart files during a run where fulldomain
C             hotstart files are specified.
C------------------------------------------------------------------------
C   r491 | stanaka | 2010-07-19 11:18:31 -0400 (Mon, 19 Jul 2010) | 1 line
C             Changed paths:
C             M /adcirc/trunk/src/adcirc.F
C             M /adcirc/trunk/src/global.F
C             M /adcirc/trunk/src/messenger.F
C             M /adcirc/trunk/src/sizes.F
C             M /adcirc/trunk/src/timestep.F
C             M /adcirc/trunk/src/writer.F
C             M /adcirc/trunk/work/makefile
C
C             Introduced writer core approach for Hot Start file outputting.
C             This function works with '-H' option, for example
C             './padcirc -H 2'. We should use this option with less than two
C             hotstart writers, because we have only two hotstart files
C             (fort.67,68). Also, this option does not support NetCDF file
C             format.
C             (jgf: let's call this v49.20).
C ------------------------------------------------------------------------
C   r488 | stanaka | 2010-07-09 12:44:19 -0400 (Fri, 09 Jul 2010) | 1 line
C             Changed paths:
C             M /adcirc/trunk/src/couple2swan.F
C             M /adcirc/trunk/src/sizes.F
C             M /adcirc/trunk/src/writer.F
C
C             Modified the writer cores algorithm to round robin outputting,
C             and adding spliting output file option. For safety outputting,
C             this round robin outputting works only with spliting option. If
C             you want to use the writer cores with this option, you should use
C             '-Ws' option instead of original writer cores option '-W'.
C             (jgf: let's call this v49.19).
C ------------------------------------------------------------------------
C   r487 | ndill | 2010-07-07 11:10:40 -0400 (Wed, 07 Jul 2010) | 1 line
C             Changed paths:
C             M /adcirc/trunk/wind/aswip_1.0.3.F
C             fixed format for reading in atcf when abs(longitude) > 100
C             (jgf: let's call this v49.18).
C ------------------------------------------------------------------------
C  v49.17   - 06/28/2010 - jgf - Fixed a bug in hotlocalize where IM values of
C             greater than 10 triggered loading of 3D hotstart info.
C  v49.16   - 06/18/2010 - tcm -- Implemented pressure ramping differently,
C             as per Rick Luettich's suggestion:
C             We now set a background pressure at 1013 mb, then compute
C             the difference between the applied pressure and the
C             background pressure.  Then apply the meterological ramp to
C             that difference and add it back to the background pressure.
C             Fixed a couple of typos in hstart.F where RampMete was
C             used instead of RampMete2. This fix, also fixes the problem
C             with the prmin for cold-started runs.
C  v49.15   - 06/16/10 - tcm Expanded the functionality of store63 and fixed
C             a bug in that same routine.  store63 could now be used in
C             place of storeOne.
C  v49.14   - 06/04/10 - tcm Fixed a bug in write_output.F related to
C             globalizing NOFF (an elemental quantity) for parallel jobs
C             which caused hotstart runs to behave different than coldstart.
C             This bug was first worked on in v48.4638 and v48.4642.
C  v49.13   - 04/16/10 - jgf - Commented out debugging statements in adcprep.F,
C             prep.F, and read_global.F.
C  v49.12   - 04/16/10 - jgf - Commented out debugging statement in adcprep.F.
C  v49.11   - 03/23/10 - rjw - r431 -Changed KSlip formula in vsmy.F to use
C             q(n,1) for bottom velocity assuming that it occurs at height
C             delta z off the bed, where delta z is thickness of bottom
C             vertical layer
C  v49.10   - 03/17/10 - jgf - r424 -Changed WindDragLimit back to 0.0035
C             from 0.0020 in preparation for merge back to trunk.
C             r419 | cdietrich | 2010-03-09
C             Fixing an error in the passing of Manning's n values as Madsen
C             roughness lengths to SWAN.
C             r417 | cdietrich | 2010-03-02
C             Changes to adcpost to allow for meshes with 10M+ elements.
C             r415 | cdietrich | 2010-02-16
C             Added wave_refraction_in_swan as a nodal attribute in the
C             fort.13 file, so wave refraction can be enabled on a
C             vertex-by-vertex basis.
C             r412 | cdietrich | 2010-02-10
C             One last change for meshes with 10M+ elements.
C             r411 | cdietrich | 2010-02-09
C             Changing formatted writing and reading of element numbers
C             inside the fort.18 file from I8 to I12 to allow for meshes
C             with 10M+ elements.
C             r410 | cdietrich | 2010-02-08
C             Small change to SWAN to allow for hot-starting of the coupled
C             model.
C             r408 | cdietrich | 2010-02-04
C             Fixing bug with the output of the radiation stress gradients
C             to the rads.64 file.
C             r407 | cdietrich | 2010-02-02
C             Fixing bug related to writing of last record in SWAN output
C             files.  Needed to flush the writer cores from within the
C             SwanOutput routine.
C             r403 | cdietrich | 2010-01-27
C             Turning off wave refraction for further testing.
C             r391 | cdietrich | 2010-01-20
C             Implement coupling of bottom friction from ADCIRC to SWAN.
C             Our Manning's n values are used to compute a physical
C             roughness length, which is passed to SWAN and used with its
C             Madsen formulation.  Now bottom friction in SWAN varies
C             spatially and is tied to the same roughness information used
C             by ADCIRC.
C             r390 | cdietrich | 2010-01-20
C             Small change to wave-modified friction in ADCIRC.
C             r389 | cdietrich | 2010-01-20
C             Lowering the cap on the wind drag coefficient from 0.0035
C             to 0.002, to be in line with results from Mark Powell.
C             r388 | cdietrich | 2010-01-20
C             Updating SWAN and ADCIRC compiler flags for the ERDC Diamond
C             machine.
C             r379 | cdietrich | 2010-01-11
C             Updating SWAN compiler flags for the ERDC Sapphire machine.
C             r373 | cdietrich | 2010-01-06
C             Removed SWAN's check for skinny elements.
C             r366 | cdietrich | 2009-12-16
C             Committing changes related to the wave-modified friction.
C             Now the Manning n inside ADCIRC is increased in regions with
C             heavy wave breaking in SWAN.  The wave-modified friction is
C             activated with an extra CSWANFRIC compiler flag.
C             r328 | cdietrich | 2009-10-19
C             Committing SWAN half of coupling changes to repository.
C             r327 | cdietrich | 2009-10-19
C             Committing ADCIRC half of coupling changes to repository.
C             r322 | cdietrich | 2009-10-13
C             Creating new branch for coupling to the unstructured SWAN.
C  v49.09     -  03/17/10 - jgf - Removed stray cr from line endings in
C                makefile.
C  v49.08     -  03/08/10 - kmd - merged 3D_kmd branch into trunk
C             -  02/22/10 - kmd - modified the station output for the
C                3D density file (fort.41) and updated the 2D compact
C                output for fort.63, fort.64, fort.73 and fort.74
C             -  01/19/10 - kmd - updates to code for the spherical
C                factors in 3D and some of the output for 3D
C               r368  |  kdresback  |  2009-12-22
C                created a branch for 3D baroclinic changes in this
C                version of the code (v49) - merged my 3D_baroclinic
C                branch with v49.06 that contained updates by rjw.
C               r332  |  kdresback  |  2009-10-22
C                updates to the messenger.F file for bug
C               r323  |  kdresback  |  2009-10-13
C                created initial branch for baroclinic code using v48
C                of the code. This branch contains major updates
C                to the transport.F, vsmy.F, hstart.F and introduces
C                a new subroutine for initial conditions (init_cond.F).
C  v49.07     - 02/04/10 - jgf - Added entry for Cray-XT5 to
C               cmplrflags.mk.
C  v49.06     - 11/25/09 - rjw - merged 3D_rjw branch into trunk
C             - 10/30/09 - rjw - modifed turbulent lenght scale for newly wet
C               nodes and fixed bug in post processing fort.43 files and
C               commented out call to massbal3D
C             - 10/06/09 - rjw - Created a branch for 3D changes and merged
C               previous changes back into this branch to be debugged and
C               then remerged into the trunk
C               r285 | rweaver | 2009-09-01 12:31:04 -0400 (Tue, 01 Sep 2009)
C               modified vsmy.F and fixed bug in vsmy.F
C               r208 | rweaver | 2009-06-18 16:53:47 -0400 (Thu, 18 Jun 2009)
C               fix bug in Holland wind NWS=8 computation of translational velocity
C               r178 | rweaver | 2009-04-15 11:09:31 -0400 (Wed, 15 Apr 2009)
C               fixed format statements in post.f
C               r136 | rweaver | 2009-03-06 10:47:47 -0500 (Fri, 06 Mar 2009)
C               Commit changes in src and prep to 3D turbulent calculations
C  v49.05     - 11/24/09 - tcm - Renamed reserved words (count became kount,
C                  index became indx, int became in2t and sum became summ).
C                  Index in particular was being used as a user defined integer
C                  and the intrinsic function in the same routine.
C  v49.04     - 11/24/09 - tcm - Fixed a bug with type NWS=3 and NWS=6 winds
C                in prep.F by removing them from the adcprep process and
C                instead of accessing local copies of the global wind file,
C                ADCIRC/PADCIRC now access the global wind file directly which
C                must now be in the global input directory.
C  V49.03     - 11/10/09 - tcm - Adding a couplers directory where coupler
C                routines and modules can be added to allow for coupling
C                ADCIRC to other models.  Populated with an empty directory
C                to coupling with the wave model STWAVE.
C  V49.02     - 10/29/09 - tcm - Fixed bugs in hstart.F related to reading in
C                wind records for NWS = +2,+4,+5,+7,+12.  This bug caused the
C                wind records to be behind by one record when hot starting.
C                Added an initialization for ITHS = 0 in adcirc.F for cold
C                start cases. Fixed format statement in owiwind.F.
C  v49.01     - 10/28/09 - TCM - Added else to initialize numSkipSnaps to
C                be 0 in owiwind.F (This fix did not get carried over with
C                the remerge from fastfloodplain to V49.)  Also changed
C                reserved word count to kount in netcdfio.F to avoid
C                conflicts.
C  v49.00     - 09/29/09 - jgf - Incremented version number to v49.00 to
C               kick off the development of the next version of ADCIRC.
C  v48.4642   - 09/18/09 - jgf - Changed syntax for correcting NOFF in
C               write_output.F, since g95 choked on the original syntax.
C               Also fixed several bugs related to writing output
C               (particularly hotstart output) to subdomain files during
C               parallel execution.
C  v48.4641   - 09/09/09 - jgf - Fixed a potential issue in globalio.F where
C               the fortran intrinsic function 'size' was used as a variable
C               name in several older subroutines; replaced these instances
C               with the variable name 'size_this'.
C  v48.4640   - 09/04/09 - rjw - Changed aswip to persist Cent.press if a
C               zero valued for the central pressure is given
C  v48.4639   - 09/02/09 - jgf - Added feature that will read in min/max
C               files on a hotstart to preserve the full record of min and
C               max values; set NWS8 to choke if Rmax or Cp are set to zero
C               on any line of fort.22, including OFCL ("forecast").
C  v48.4638   - 09/02/09 - tcm - Fixed a bug in writing hot start values for
C               NOFF which could have the wrong values after globalizing
C               because there are no rules for which processor owns an
C               element and Noff is an elemental value. The way around is to
C               take the minimum of (Noff, 1) as NOFF is supposed to be
C               either 0 or 1.
C  v48.4637   - 08/31/09 - rjw - Fixed out of bounds array in aswip
C  v48.4636   - 08/24/09 - rjw - Fixed Bug in NWS 19 code (first lon not
C               properly converted)
C  v48.4635   - 08/14/09 - tcm - Fixed bug in call to neighb within read_input
C               updated the adcirc makefiles for a pc
C  v48.4634   - 08/13/09 - tcm - Fixed bug in globalio where multiple variables
C               where repeated in use module line.
C  v48.4633   - 08/05/09 - rjw - Fixed Bug in translation direction in ASWIP output
C  v48.4632   - 08/05/09 - rjw - Fixed formatting error in ASWIP output
C               to make it compatible with wind.F and should still be
C               compatible with utilities that parse ATCF output.
C  v48.4631   - 07/30/09 - rjw - Modified aswip set background pressure to
C               hardwired value of 1013
C  v48.4630   - 07/28/09 - jgf - Modified aswip output to include commas
C               to make it compatible with utilities that parse ATCF output.
C  v48.4629   - 07/28/09 - rjw - Added storm speed(kts) and direction(degrees)
C               information to aswip.F output fort.22 file
C               and the translation velocities (U and V) are now computed from
C               the speed and direction that is read in from
C               the new fort.22 file into NWS19 code
C  v48.4628   - 07/27/09 - jgf - Extended METONLY capability to cases where
C               3D or harmonic output is requested, as well as the case where
C               execution begins as a hotstart.
C  v48.4627   - 07/27/09 - jgf - Created METONLY logical variable: if only
C               met output is requested, ADCIRC will skip solution of the
C               GWCE and the momentum equations.
C  v48.4626   - 07/27/09 - jgf - fixed bug in presizes.F if fort.15 file is
C               not found.
C  v48.4625   - 07/22/09 - rjw - fixed bugs in NWS19 code in wind.F and aswip.F
C  v48.4624   - 07/22/09 - rjw - fixed bugs in NWS19 code in wind.F and aswip.F
C  v48.4623   - 07/22/09 - rjw - fixed confusing lines and file names in
C               ASymmetric WInd Preprocessing (ASWIP) program
C  v48.4622   - 07/22/09 - rjw - Added NWS=19 Asymmetric Holland v2.0
C               to src/ files
C  v48.4621   - 07/22/09 - rjw - Added NWS=19 Asymmetric Holland v2.0
C               to prep/ files
C  v48.4620   - 07/21/09 - jgf - Modified makefile so that aswip executable
C               may be built with 'make aswip'.
C  v48.4619   - 07/21/09 - jgf - Fixed inclusion of MNMPROC in write_output.F;
C               merged lumped/explicit mode from Seizo Tanaka; rjw
C               merged aswip.F to generate metadata for asymmetric
C               vortex model.
C  v48.4618    -07/10/09 - tcm - Fixed bug in writing out global elevation
C               values using -99999.0 that exists in serial version
C               Also added 3 fields to the OutputDataDescript_t type
C  v48.4617    -06/19/09 - jgf - Made numeric constants in NWS8 double
C               precision; in NWS9 sub added initialization of told and
C               changed time.eq.time2 to time.gt.time2.
C  v48.4616    -06/18/09 - rjw - Modified the NWS=8 Holland wind to fix
C               bug in translational velocity to account for direction.
C  v48.4615    -06/02/09 - rjw - Removed erroneous call to
C               COMPUTE_RESELEM from read_input.F.
C  v48.4614   - 05/27/09 - mcf - Fixed implementation of Rick Luettich
C               V/Vmax tapering formulation (as in NWS=8)
C               for future experimentation - not activated
C  v48.4613   - 05/27/09 - jgf - Commented out debug lines in
C               write_output.F and globalio.F that were filling up
C               fort.16.
C  v48.4612   - 05/26/09 - rjw - Fix format statement in postprocess for 3D.
C  v48.4611   - 05/22/09 - mcf - Reverted damping of translational velocity
C               back to the original formulation
C  v48.4610   - 05/22/09 - jgf - Changed air density back to 0.001293 in
C               stress calc for NWS8 and NWS9 in timestep.F.
C  v48.4609   - 05/21/09 - mcf - Fixed seg-fault when writing netcdf hotstart
C               files, incorporated fix in the define portion of netcdfio
C               for grids with no specified boundary forcing segments or nodes
C  v48.4608   - 05/19/09 - mcf - Implemented an old bug fix: units conversion
C               in asymmetric wind code (NWS=9)
C  v48.4607   - 05/18/09 - jgf - Fixes to get adcprep and adcirc to compile
C               and link with netCDF on Jade (Cray-XT4) at ERDC.
C  v48.4606   - 05/13/09 - rjw - modify makefile to compile adcprep using
C               version.F
C  v48.4605   - 05/12/09 - rjw - modify asymmetric code to use damping of
C               translational velocity based on ratio of V:Vmax (pass Vmax
C               from wind.F into vortex.F subroutine uvp), and set (in vortex.F)
C               velocities at center of vortex (dist < 1m) to zero.
C  v48.4604   - 05/08/09 - jgf - Fixed typo in version.F
C  v48.4603  -  05/08/09 - rjw-  updated the distance calculation formula
C               in Holland Model (NWS=8) to more accuarately calculate
C               distance along a shpere, and modified the coriolis formula
C               to compute a spatially varying coriolis term also for the
C               Symmetric Holland model (NWS=8)
C  v48.4602   - 05/08/09 - jgf - Changed air density for symmetric and
C               asymmetric vortex models (NWS=8 and 9) to 1.15 kg/m^3,
C               for Holland B calculation as well as conversion to wind
C               stress. Changed asymmetric far field wind damping so
C               that only the translation velocity is damped.
C  v48.4601   - 05/01/09 - jgf - Fixed improperly located check for
C               netcdfio support in presizes.F; removed FileFmt etc from USE
C               PRE_GLOBAL in prep.F (now uses VERSION instead). Merged
C               prep.F from v48.48 (trunk r161) b/c it has fastprep13.
C  v48.46     - 02/26/09 - jgf - Corrected backloaded tau0 averaging
C               formula: now just one time level.
C  v48.45     - 02/26/09 - jgf - Fixed backloaded tau0 formula by
C               removing an extra alpha factor from tau0(t-3) term.
C  v48.44     - 02/02/09 - mcf - fixed time unlimited index and time variable
C               in the fort.63.nc output file for hotstart cases in netcdfio
C               module
C  v48.43     - 01/29/09 - rjw - modified fixes to massbal.F (check for
C               ICS status) and vsmy.F (initialized turbulent length scale to
C               lmin= 0.001 in cstart.F (COLDSTART_3D)) and reformatted
C               declaration statements in global_3dvs.F.
C  v48.42     - 01/28/09 - jgf - Added two new time averaging options for tau0
C               ... activated by tau0=-6.x or 7.x in fort.15.
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

