Parameter Definitions
=====================

Input Files
-----------

fort.14, fort.15
~~~~~~~~~~~~~~~~

.. _AGRID:

AGRID
   Alpha-numeric grid identification (<=24 characters).

.. _BARINCFSB:

BARINCFSB(k,j)
   Coefficient of free surface subcritical flow at internal barrier node NBVV(k,j) and the paired node IBCONN(k,j). A typical value is BARINCFSB(k,j) = 1.0.

.. _BARINCFSP:

BARINCFSP(k,j)
   Coefficient of free surface supercritical flow at internal barrier node NBVV(k,j) and the paired node IBCONN(k,j). A typical value is BARINCFSP(k,j) = 1.0.

.. _BARINHT:

BARINHT(k,j)
   Internal barrier height (positive above the geoid and negative below the geoid) at node NBVV(k,j) and the paired node IBCONN(k,j)). The barrier height must be greater than the bathymetric depth at these nodes, i.e., recalling the different sign convention between the bathymetric depth and the barrier height, BARINHT(k,j) > –DP(NBVV(k,j) and <b >BARINHT(k,j) > –DP(IBCONN(k,j)). If this is not satisfied, the run will terminate.

.. _BARLANCFSP:

BARLANCFSP(k,j)
   Coefficient of free surface supercritical flow at external barrier node NBVV(k,j). A typical value is BARLANCFSP(k,j)= 1.0.

.. _BARLANHT:

BARLANHT(k,j)
   External barrier height (positive above the geoid and negative below the geoid) at node NBVV(k,j). The barrier height must be greater than the bathymetric depth at this node, i.e., recalling the different sign convention between the bathymetric depth and the barrier height, BARLANHT(k,j) > –DP(NBVV(k,j). If this is not satisfied, the run will terminate.

.. _BK:

BK(k)
   Bridge pier shape factor (K – see Table 1 in the section on Bridge Piers)

.. _BALPHA:

BALPHA(k)
   Fraction of the cross section occupied by all of the piers in the bridge = sum of bridge pier widths/width of cross section (corresponds to in the section on Bridge Piers)

.. _BDELX:

BDELX(k)
   Approximate nodal spacing at the location of a bridge in the upstream/downstream direction in meters or feet depending on the grid coordinate system (for lon,lat coordinates, use meters for BDELX). (Note that 2*BDELX if the bridge pier effects are distributed across 3 nodes in the alongstream direction)

.. _DP:

DP(JN)
   Bathymetric depth with respect to the geoid, positive below the geoid and negative above the geoid. Bathymetric depths above the geoid or sufficiently small that nodes will dry, require that the wetting/drying option is enabled (NOLIFA=2 in the Model Parameter and Periodic Boundary Condition File.)

.. _IBCONN:

IBCONN(k,j)
   Back face node paired with the front face node, NBVV(k,j), on an internal barrier boundary

.. _IBTYPE:

IBTYPE(k)
   Boundary type
      = 0 external boundary with no normal flow as an essential boundary condition and no constraint on tangential flow. This is applied by zeroing the normal boundary flux integral in the continuity equation and by zeroing the normal velocity in the momentum equations. This boundary condition should satisfy no normal flow in a global sense and no normal flow at each boundary node. This type of boundary represents a mainland boundary with a strong no normal flow condition and free tangential slip.

      = 1 internal boundary with no normal flow treated as an essential boundary condition and no constraint on the tangential flow. This is applied by zeroing the normal boundary flux integral in the continuity equation and by zeroing the normal velocity in the momentum equations. This boundary condition should satisfy no normal flow in a global sense and no normal flow at each boundary node. This type of boundary represents an island boundary with a strong normal flow condition and free tangential slip.

      = 2 external boundary with non-zero normal flow as an essential boundary condition and no constraint on the tangential flow. This is applied by specifying the non-zero contribution to the normal boundary flux integral in the continuity equation and by specifying the non-zero normal velocity in the momentum equations. This boundary condition should correctly satisfy the flux balance in a global sense and the normal flux at each boundary node. This type of boundary represents a river inflow or open ocean boundary with a strong specified normal flow condition and free tangential slip. Discharges are specified either in the Model Parameter and Periodic Boundary Condition File for harmonic discharge forcing or in the Non-periodic, Normal Flux Boundary Condition File for time series discharge forcing.

      = 3 external barrier boundary with either zero or non-zero normal outflow from the domain as an essential boundary condition and no constraint on the tangential flow. This is applied by specifying the contribution (zero or non-zero) to the normal boundary flux integral in the continuity equation and by specifying the (zero or non-zero) normal velocity in the momentum equations. Non-zero normal flow is computed using a supercritical, free surface weir formula if the barrier is overtopped. Zero normal flow is assumed if the barrier is not overtopped. This boundary condition should correctly satisfy the flux balance in a global sense and the normal flux at each boundary node. This type of boundary represents a mainland boundary comprised of a dike or levee with strong specified normal flow condition and free tangential slip. See External Barrier Boundary Note below for further information on exterior barrier boundaries.

      = 4 internal barrier boundary with either zero or non-zero normal flow across the barrier as an essential boundary condition and no constraint on the tangential flow. This is applied by specifying the contribution (zero or non-zero) to the normal boundary flux integral in the continuity equation and by specifying the normal velocity (zero or non-zero) in the momentum equations. Non-zero normal flow is compute using either subcritical or supercritical, free surface weir formula (based on the water level on both sides of the barrier) if the barrier is overtopped. Zero normal flow is assumed if the barrier is not overtopped. This type of boundary represents a dike or levee that lies inside the computational domain with a strong specified normal flow condition and free tangential slip. See Internal Barrier Boundary Note for further information on exterior barrier boundaries.

      = 5 internal barrier boundary with additional cross barrier pipes located under the crown. Cross barrier flow is treated as an essential normal flow boundary condition which leaves/enters the domain on one side of the barrier and enters/leaves the domain on the corresponding opposite side of the barrier flow rate and direction are based on barrier height, surface water elevation on both sides of the barrier, barrier coefficient and the appropriate barrier flow formula. In addition cross barrier pipe flow rate and direction are based on pipe crown height, surface water elevation on both sides of the barrier, pipe friction coefficient, pipe diameter and the appropriate pipe flow formula. Free tangential slip is allowed

      = 10 external boundary with no normal and no tangential flow as essential boundary conditions. This is applied by zeroing the normal boundary flux integral in the continuity equation and by setting the velocity = 0 rather than solving momentum equations along the boundary. This boundary condition should satisfy no normal flow in a global sense and zero velocity at each boundary node. This type of boundary represents a mainland boundary with strong no normal flow and no tangential slip conditions.

      = 11 internal boundary with no normal and no tangential flow as essential boundary conditions. This is applied by zeroing the normal boundary flux integral in the continuity equation and by setting the velocity = 0 rather than solving momentum equations along the boundary. This boundary condition should correctly satisfy no normal flow in a global sense and zero velocity at each boundary node. This type of boundary represents an island boundary with strong no normal flow and no tangential slip conditions.

      = 12 external boundary with non-zero normal and zero tangential flow as an essential boundary condition. This is applied by specifying the non-zero contribution to the normal boundary flux integral in the continuity equation and by setting the non-zero normal velocity and zero tangential velocity rather than solving momentum equations along the boundary. This boundary condition should correctly satisfy the flux balance in a global sense and the specified normal/zero tangential velocity at each boundary node. This type of boundary represents a river inflow or open ocean boundary in which strong normal flow is specified with no tangential slip. Discharges are specified either in the Model Parameter and Periodic Boundary Condition File for harmonic forcing or in the Non-periodic, Normal Flux Boundary Condition File for time series forcing.

      = 13 external barrier boundary with either zero or non-zero normal outflow from the domain and zero tangential flow as essential boundary conditions. This is applied by specifying the contribution (zero or non-zero) to the normal boundary flux integral in the continuity equation and by setting the (zero or non-zero) normal velocity and zero tangential velocity rather than solving momentum equations along the boundary. Non-zero normal flow is computed using a supercritical, free surface weir formula if the barrier is overtopped. Zero normal flow is assumed if the barrier is not overtopped. This boundary condition should correctly satisfy the flux balance in a global sense and the normal velocity/zero tangential velocity at each boundary node. This type of boundary represents a mainland boundary comprised of a dike or levee with strong specified normal flow and no tangential slip conditions. See External Barrier Boundary Note below for further information on exterior barrier boundaries.

      = 20 external boundary with no normal flow as a natural boundary condition and no constraint on tangential flow. This is applied by zeroing the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) applied in the momentum equations. This boundary condition should satisfy no normal flow in a global sense, but will only satisfy no normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents a mainland boundary with a weak no normal flow condition and free tangential slip.

      = 21 internal boundary with no normal flow as a natural boundary condition and no constraint on the tangential flow. This is applied by zeroing the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. This boundary condition should satisfy no normal flow in a global sense but will only satisfy no normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents an island boundary with a weak no normal flow condition and free tangential slip.

      = 22 external boundary with non-zero normal flow as a natural boundary condition and no constraint on the tangential flow. This is applied by specifying the non-zero contribution to the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. This boundary condition should correctly satisfy the flux balance in a global sense but will only satisfy the normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents a river inflow or open ocean boundary with a weak specified normal flow condition and free tangential slip. Discharges are specified either in the Model Parameter and Periodic Boundary Condition File for harmonic discharge forcing or in the Non-periodic, Normal Flux Boundary Condition File for time series discharge forcing.

      = 23 external barrier boundary with either zero or non-zero normal outflow from the domain as a natural boundary condition and no constraint on the tangential flow. This is applied by specifying the contribution (zero or non-zero) to the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. Non-zero normal flow is computed using a supercritical, free surface weir formula if the barrier is overtopped. Zero normal flow is assumed if the barrier is not overtopped. This boundary condition should correctly satisfy the flux balance in a global sense but will only satisfy the normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents a mainland boundary comprised of a dike or levee with a weak specified normal flow condition and free tangential slip. See External Barrier Boundary Note below for further information on exterior barrier boundaries.

      = 24 internal barrier boundary with either zero or non-zero normal flow across the barrier as a natural boundary condition and no constraint on the tangential flow. This is applied by specifying the contribution (zero or non-zero) to the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. Non-zero normal flow is compute using either subcritical or supercritical, free surface weir formula (based on the water level on both sides of the barrier) if the barrier is overtopped. Zero normal flow is assumed if the barrier is not overtopped. This boundary condition should correctly satisfy the flux balance in a global sense but will only satisfy the normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents a dike or levee that lies inside the computational domain with a weak specified normal flow condition and free tangential slip. See Internal Barrier Boundary Note below for further information on exterior barrier boundaries.

      = 25 internal barrier boundary with additional cross barrier pipes located under the crown. Cross barrier flow is treated as a natural normal flow boundary condition which leaves/enters the domain on one side of the barrier and enters/leaves the domain on the corresponding opposite side of the barrier. Flow rate and direction are based on barrier height, surface water elevation on both sides of the barrier, barrier coefficient and the appropriate barrier flow formula. In addition cross barrier pipe flow rate and direction are based on pipe crown height, surface water elevation on both sides of the barrier, pipe friction coefficient, pipe diameter and the appropriate pipe flow formula. Free tangential slip is allowed.

      = 30 wave radiation normal to the boundary as a natural boundary condition. This is applied by specifying the contribution to the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. The normal flow is computed using a Sommerfield radiation condition. This boundary condition should correctly satisfy the flux balance in a global sense but will only satisfy the normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents an open boundary where waves are allowed to propagate freely out of the domain.

      = 32 a combined specified normal flux and outward radiating boundary. The GWCE is forced with the total normal flux computed by adding the specified normal flux and the flux associated with the outward radiating wave. The latter is determine from a Sommerfeld type condition, flux=celerity*wave elevation. The momentum equations are used to compute the velocity field the same as for a nonboundary node.

      = 40 a zero normal velocity gradient boundary. The GWCE is forced with normal flux, the momentum eqs are sacrificed in favor of setting the velocity at a boundary node equal to the value at a fictitious point inside the domain. The fictitious point lies on the inward directed normal to the boundary a distance equal to the distance from the boundary node to its farthest 'neighbor. This should ensure that the fictitious point does not fall into an element that contains the boundary node. The velocity at the fictitious point is determined by interpolation.

      = 41 a zero normal velocity gradient boundary. The GWCE is forced with normal flux. The momentum eqs are sacrificed in favor of eqs that set the velocity gradient normal to the boundary equal to zero in the Galerkin sense.

      = 52 external boundary with periodic non-zero normal flow combined with wave radiation normal to the boundary as natural boundary conditions and no constraint on the tangential flow. This is applied by specifying the non-zero contribution to the normal boundary flux integral in the continuity equation. There is no constraint on velocity (normal or tangential) in the momentum equations. This boundary condition should correctly satisfy the flux balance in a global sense but will only satisfy the normal flow at each boundary node in the limit of infinite resolution. This type of boundary represents a periodic river inflow or open ocean boundary with a weak specified normal flow condition and free tangential slip where waves are allowed to propagate freely out of the domain. Discharges are specified in the Model Parameter and Periodic Boundary Condition File as harmonic discharge forcing. Additional parameters, including DRampExtFlux and FluxSettlingTime must also be set in the Model Parameter and Periodic Boundary Condition File in order to use this boundary type.

      = 64 vertical element wall boundary that allows the mesh to have two nodes at the same horizontal location with different depths and allows water to seemlessly flow over the boundary, by primarily consolidating the nodal equations and occasionally using the weir formula. This boundary type is useful to have steep-sided channels or other vertical features represented in the mesh while, unlike barrier boundaries that entirely depend on the weir formula, solutions along this boundary are primarily computed based on the governing equations. See :doc:`../special_features/vertical_element_walls` for more information.

      = 102, 112, or 122 flux specified baroclinic. In order to designate a river boundary as baroclinic, 100 should be added to the IBTYPE that would be appropriate in the barotropic case. For example, to convert a barotropic river boundary (IBTYPE of 22) to a baroclinic river boundary with freshwater inflow, change the IBTYPE to 122. If there is a 1 in the 100s place of the IBTYPE, ADCIRC will then try to read an input file (fort.39) for the salinity and/or temperature boundary condition. The format of the fort.39 file depends on the value of IDEN; see the documentation of the fort.39 file for more details.

.. _IBTYPEE:

IBTYPEE
   Elevation boundary type, (At present the only allowable value for IBTYPEE is 0). Elevations are specified either in the Model Parameter and Periodic Boundary Condition File for harmonic forcing or in the Non-periodic Elevation Boundary Condition File for time series forcing.

.. _PIPEHT:

PIPEHT(k,j)
   The height of the pipe above the geoid (m) at the j-th node in the k-th normal flow specified boundary segment. Only used with IBTYPE 5 or 25.

.. _PIPECOEF:

PIPECOEF(k,j)
   The bulk pipe friction factor for the internal barrier boundary cross-barrier pipe boundary node pair. This bulk friction factor is defined as: PIPECOEFR = fL/D where f is the classical pipe friction coefficient, L is the length of the pipe through the barrier (in consistent units) and D is the diameter of the cross barrier pipe (in consistent units).

.. _PIPEDIAM:

PIPEDIAM(k,j)
   Cross barrier pipe diameter in internal barrier with cross barrier pipes (diameter in units consistent with the mesh).

.. _JE:

JE
   Element number. The elements must be input in ascending order.

.. _JN:

JN
   Node number. The nodes must be input in ascending order.

.. _NBDV:

NBDV(k,j)
   Node numbers on elevation specified boundary segment k. The node numbers must be entered in order as one moves along the boundary segment, however the direction (counter clockwise or clockwise around the domain) does not matter.

.. _NBPNODES:

NBPNODES
   Total number of nodes (centerline and adjacent) in ADCIRC grid used to represent the effects of bridge piers

.. _NBNNUM:

NBNNUM(k)
   Node number in the ADCIRC grid of node k used to represent the frictional effects of bridge piers

.. _NBOU:

NBOU
   Number of normal flow (discharge) specified boundary segments. These include zero normal flow (land) boundaries, non-zero normal flow (river) boundaries, potentially overflowing external barrier boundaries, potentially overtopping internal barrier boundaries and radiation boundaries. See also General Notes for Normal Flow Boundary Conditions.

.. _NBVV:

NBVV(k,j)
   Node numbers on normal flow boundary segment k. The node numbers must be entered in order as one moves along the boundary segment with land always on the right, i.e., in a counter clockwise direction for external (e.g., mainland, external barrier) boundaries and a clockwise direction for internal (e.g., island, internal barrier) boundaries. For an internal barrier boundary (IBTYPE(k) = 4, 24) only the nodes on the front face of the boundary are specified in NBVV(k,j). The paired nodes on the back face of the boundary are specified in IBCONN(k,j).

.. _NE:

NE
   Number of elements in the horizontal grid

.. _NETA:

NETA
   Total number of elevation specified boundary nodes

.. _NHY:

NHY
   Number of nodes per element. At present the only allowable value for the number of nodes per element is 3 indicating a triangular element with linear basis functions.

.. _NM:

NM(JE,1), NM(JE,2), NM(JE,3)
   Node numbers comprising element JE. These must be specified in a counter clockwise direction moving around the element.

.. _NOPE:

NOPE
   Number of elevation specified boundary forcing segments.

.. _NP:

NP
   Number of nodes in the horizontal grid

.. _NVDLL:

NVDLL(k)
   Number of nodes in elevation boundary segment k

.. _NVEL:

NVEL
   Total number of normal flow specified boundary nodes including both the front and back nodes on internal barrier boundaries.

.. _NVELL:

NVELL(k)
   Number of nodes in normal flow specified boundary segment k. For an internal barrier boundary (IBTYPE 4 or 24), NVELL(k) includes only the nodes on the front face of the boundary as specified in NBVV(k,j) and not the paired nodes on the back face of the boundary specified in IBCONN(k,j).

.. _X:

X(JN), Y(JN)
   X and Y coordinates. If ICS=1 in the Model Parameter and Periodic Boundary Condition File then X(JN), Y(JN) are Cartesian coordinates with units of length (e.g., feet or meters) that are consistent with the definition of gravity in the Model Parameter and Periodic Boundary Condition File. If ICS=2 in the Model Parameter and Periodic Boundary Condition File then X(JN), Y(JN) are spherical coordinates in degrees of longitude (east of Greenwich is positive and west of Greenwich is negative) and degrees of latitude (north of the equator is positive and south of the equator is negative)

.. _Y:

Y(JN)
   See X(JN) above.

.. _RUNDES:

RUNDES
   Alpha-numeric run description 1 (<=32 characters)

.. _RUNID:

RUNID
   Alphanumeric run description 2 (<= 24 characters)

.. _NFOVER:

NFOVER
   Non-fatal error override option;

   = 0, inconsistent input parameters will cause program termination.

   = 1, inconsistent input parameters will, (when possible), be automatically corrected to a default or consistent value and execution continued. Be sure to read the nonfatal warning messages to see whether any parameters have been modified. Note that not all inconsistent parameters can be corrected automatically and therefore fatal error messages and program termination may still result.

   Note for NFOVER:

   Occasionally, the elevation solution becomes unphysically large due to improper numerical parameter settings, time step stability criteria violation. It can be useful for ADCIRC to terminate based on a user-specific limit to the computed water elevation. To enable this feature, ADCIRC must be compiled with the DEBUG_WARN_ELEV compiler directive. This is set in the cmplrflags.mk file. For example, for the serial model,

   DA := -DREAL8 -DLINUX -DCSCA –DDEBUG_WARN_ELEV

   and for the parallel model:

   DP := -DREAL8 -DLINUX -DCSCA -DCMPI -DDEBUG_WARN_ELEV

   This enables extra parameters in the fort.15 file, specified on the NFOVER line:

   NFOVER, WarnElev, iWarnElevDump, WarnElevDumpLimit, ErrorElev

   ADCIRC then monitors the maximum water elevation and behaves as follows:

   A warning is issued if the elevation exceeds WarnElev.

   A global elevation file (written to fort.69) will be written if WarnElev is exceeded AND iWarnElevDump == 1

   Execution will be terminated if WarnElevDumpLimit global elevation files have been written due to the above warning limits.

   Execution will be terminated if elevation exceeds ErrorElev.

   The default values are:

   WarnElev = 20.0 ! Warn at 20 meters

   iWarnElevDump = 0 ! Do not write global elevation files

   WarnElevDumpLimit = 50 ! Terminate execution if the warning level is reached 50 times

   ErrorElev = 1000.0 ! Terminate execution of the water elevation exceeds 1000 meters

   Example: to override the defaults, recompile ADCIRC as above, and modify the NFOVER line in the fort.15 file:

   1 10. 0 100 30.

   This will cause ADCIRC to warn when the water elevation exceeds 10 meters, no global elevation files will be written, and ADCIRC will terminate if 100 warnings are generated OR the elevation exceeds 30 meters.

.. _NABOUT:

NABOUT
   Logging level for output from ADCIRC to the screen or console as well as the ADCIRC log file (fort.16). ADCIRC writes log messages at 5 levels of severity: DEBUG, ECHO, INFO, WARNING, and ERROR (from lowest or least important to highest or most important). Selection of a logging level indicates that messages of that level and higher should be logged. Setting the logging level to WARNING or ERROR will reduce the size (and clutter) in the log files, but important information could be missed. ERROR messages generally result from problems that also cause the run to stop.

   NABOUT=-1:: DEBUG-level log messages and higher. This may consume a lot of disk space and slow ADCIRC down, perhaps dramatically. Generally only useful for ADCIRC developers.

   NABOUT=0:: ECHO-level log messages and higher. ECHO-level log messages include echo printing of most input files including the fort.13, fort.14 and fort.22 files.

   NABOUT=1:: INFO-level log messages and higher. These messages inform the user about something that ADCIRC has done that is important but is not the result of a problem or issue.

   NABOUT=2:: WARNING-level log messages and higher. These messages indicate a potential problem that is generally not fatal to the run.

   NABOUT=3:: ERROR-level log messages only. These messages indicate a severe problem that usually causes the run to stop.

   Prior to ADCIRC version 49, 0 (ECHO) and 1 (INFO) were the only available options. NABOUT originally stood for "Abbreviated Output" (to log files).

.. _NSCREEN:

NSCREEN
   Controls log message output to the screen (i.e., to standard output). Timestep logging will be written every abs(NSCREEN) timesteps. Output to the screen consists mainly of timestep logging.

   NSCREEN<0:: Log messages that would normally be wrtten to the screen are written to a file called adcirc.log instead.

   NSCREEN=0:: Log messages will not be written to the screen.

   NSCREEN>0:: Log messages are written to the screen (a.k.a. standard out).

.. _IHOT:

IHOT
   Parameter controlling whether the model is hot started. The hotstart facility is available for 2D and 3D runs. The hotstart file will also contain harmonic analysis data if harmonic analysis was underway, so that the harmonic analysis can be hotstarted as well.

   = 0 cold start the model

   = 17 hot start from ascii file fort.17

   = 67 hot start model using input information in hot start file fort.67

   = 68 hot start model using input information in hot start file fort.68

   = 367 hot start model using input information in netCDF hot start file fort.67.nc

   = 368 hot start model using input information in netCDF hot start file fort.68.nc

   = 567 hot start model using input information in netCDF4 hot start file fort.67.nc

   = 568 hot start model using input information in netCDF4 hot start file fort.68.nc

.. _ICS:

ICS
   Parameter controlling whether the model is run in spherical or Cartesian coordinates.

   = 1 ADCIRC governing equations are in standard Cartesian coordinates. Coordinates in the grid file (fort.14) are assumed to have units of length that are consistent with the units of gravity (as specified below). In the unlikely case that tidal potential forcing (NTIP=1 or 2) and/or a spatially variable Coriolis coefficient (NCOR=1) are desired for this type of run, an inverse map projection (Carte Parallelo-grammatique) is used to obtain longitude and latitude values for the grid. However, we strongly recommend that if the model domain is large enough for either spatially variable Coriolis or tidal potential forcing to be considered important, the model should be run with spherical governing equations (ICS=2) using a longitude, latitude grid.

   = 2 ADCIRC governing equations are in spherical coordinates transformed into Cartesian coordinates prior to discretization using a map projection (Carte Parallelo-grammatique – CPP). Coordinates in the grid file (fort.14) are in decimal degrees longitude and latitude.

.. _IM:

IM
   Model type

   = 0 Barotropic 2DDI run using New GWCE and Momentum equation formulations

   = 1 Barotropic 3D run using New GWCE and velocity based Momentum equations

   = 21 Baroclinic 3D run using New GWCE and velocity based Momentum equations

   = 111112 Barotropic 2DDI run using the lumped GWCE (instead of the default fully consistent GWCE). This option is needed to run ADCIRC in lumped + explicit mode, thereby bypassing the iterative solver. Explicit mode also requires specifying coefficients A00, B00, C00 (= 0, 1, 0) in this file. (see below)

   = 611112 Barotropic 3D run using the lumped GWCE (instead of the default fully consistent GWCE). This option is needed to run ADCIRC in lumped + explicit mode, thereby bypassing the iterative solver. Explicit mode also requires specifying coefficients A00, B00, C00 (= 0, 1, 0) in this file. (see below)

.. _IDEN:

IDEN
   Form of density forcing in a 3D run (for all Baroclinic model runs, the initial density, temperature and/or salinity field is read in from UNIT 11); ADCIRC does not support 2D baroclinic model runs.

   =-4 Diagnostic Baroclinic ADCIRC run with Salinity and Temperature forcing

   =-3 Diagnostic Baroclinic ADCIRC run with Temperature forcing

   =-2 Diagnostic Baroclinic ADCIRC run with Salinity forcing

   =-1 Diagnostic Baroclinic ADCIRC run with Sigma T forcing

   =0 Barotropic model run

   =1 Prognostic Baroclinic ADCIRC run with Sigma T forcing

   =2  Prognostic Baroclinic ADCIRC run with Salinity forcing

   =3 Prognostic Baroclinic ADCIRC run with Temperature forcing

   =4 Prognostic Baroclinic ADCIRC run with Salinity and Temperature forcing

.. _NOLIBF:

NOLIBF
   Parameter controlling the type of bottom stress parameterization used in a 2DDI ADCIRC run. This parameter must be specified but is ignored in a 3D run. Note: In the NWP section, if the user selects quadratic_friction_coefficient_at_sea_floor, mannings_n_at_sea_floor, or chezy_friction_coefficient_at_sea_floor, then NOLIBF must be 1 (nonlinear friction formulation) since all those formulations are nonlinear. If the NOLIBF were anything other than 1, it is an error that will cause ADCIRC to stop.

   = 0 linear bottom friction law. The friction coefficient (FFACTOR) is specified below.

   = 1 quadratic bottom friction law. The friction coefficient (FFACTOR) is specified below.

   = 2 hybrid nonlinear bottom friction law. In deep water, the friction coefficient is constant and a quadratic bottom friction law results. In shallow water the friction coefficient increases as the depth decreases (e.g. as in a Manning-type friction law). The friction coefficient is determined as: FFACTOR=FFACTORMIN*[1+(HBREAK/H)**FTHETA]**(FGAMMA/FTHETA). The required parameters (FFACTORMIN, HBREAK, FTHETA, FGAMMA) are specified below.

   For 3D ADCIRC runs, spatially varying bottom friction should be specified using the bottom_roughness_length nodal attribute.

.. _NOLIFA:

NOLIFA
   Parameter controlling the finite amplitude terms in ADCIRC. The value of NOLIFA effects the meaning of the minimum water depth parameter(H0) and requires the specification of additional parameters together with H0, (see below). When the finite amplitude terms are turned on, the time derivative portion of the advective terms should also be turned on for proper mass conservation and consistency (i.e. when NOLIFA>0, then NOLICAT=1).

   = 0 finite amplitude terms ARE NOT included in the model run (i.e., the depth is linearized by using the bathymetric depth, rather than the total depth, in all terms except the transient term in the continuity equation) and wetting and drying of elements is disabled. Initial water depths are assumed equal to the bathymetric water depth specified in the grid file (fort.14).

   = 1 finite amplitude terms ARE included in the model run and wetting and drying of elements is disabled. Initial water depths are assumed equal to the bathymetric water depth specified in the grid file (fort.14).

   = 2 finite amplitude terms ARE included in the model run and wetting and drying of elements is enabled. Initial water depths are assumed equal to the bathymetric water depth specified in the grid file (fort.14).

.. _NOLICA:

NOLICA
   Parameter controlling the advective terms in ADCIRC (with the exception of a time derivative portion that occurs in the GWCE form of the continuity equation and is controlled by NOLICAT). When these (spatial derivative) portions of the advective terms are included, the time derivative portion of the advective terms in the GWCE should also be included (i.e. when NOLICA=1, NOLICAT=1).

   = 0 advective terms ARE NOT included in the computations

   = 1 advective terms ARE included in the computations

   The NOLICA and NOLICAT parameters by themselves will activate or deactivate the advective terms over the full domain; however, these terms can also be activated or deactivated on an element-by-element basis. Please see the definition of the advection_state nodal attribute in the documentation for the Nodal Attributes File (fort.13) for details.

.. _NOLICAT:

NOLICAT
   Parameter controlling the time derivative portion of the advective terms that occurs in the GWCE form of the continuity equation in ADCIRC. The remainder of the advective terms in the GWCE and the entire advective terms in the momentum equation are controlled by NOLICA. These terms should be included if either the finite amplitude or the remainder of the advective terms are included to maintain mass conservation and solution consistency.

   = 0 the time derivative portion of the advective terms that occur in the GWCE continuity equation ARE NOT included in the computations.

   = 1 the time derivative portion of the advective terms that occur in the GWCE continuity equation ARE included in the computations.

   The NOLICA and NOLICAT parameters by themselves will activate or deactivate the advective terms over the full domain; however, these terms can also be activated or deactivated on an element-by-element basis. Please see the definition of the advection_state nodal attribute in the documentation for the Nodal Attributes File (fort.13) for details.

.. _NWP:

NWP
   Number of nodal attributes used in the run. Nodal attributes are properties of each node in the grid and are spatially varying but constant in time. See AttrName for examples. The nodal attribute data itself must be provided by the user in the Nodal Attributes File (fort.13). If NWP is not zero, then the names of the nodal attributes appear on the following lines, one per line:

   FOR j=1 to NWP

   AttrName(j)

   end j loop

.. _NCOR:

NCOR
   Parameter controlling whether the Coriolis parameter is constant in space and read in below or spatially varying as computed from the y-coordinates of the nodes in the grid (assumed to be in degrees Latitude). The grid coordinate system is specified by the ICS parameter (see above).

   = 0, to read in a spatially constant Coriolis parameter

   = 1, to compute a spatially variable Coriolis parameter

.. _NTIP:

NTIP
   Parameter controlling whether tidal potential and self attraction/load tide forcings will be used to force ADCIRC.

   = 0, tidal potential & self attraction/load tide forcings are not used

   = 1, tidal potential forcing is used

   = 2, tidal potential & self attraction/load tide forcings are used. In this case the self attraction/load tide information is read in for each constituent at each node in the grid from the Self Attraction/Earth Load Tide Forcing File.

.. _NWS:

NWS
   Parameter controlling whether wind velocity or stress, wave radiation stress and atmospheric pressure are used to force ADCIRC.

   = 0, no wind, radiation stress or atmospheric pressure forcings are used.

   = 1, wind stress and atmospheric pressure are read in at all grid nodes at every model time step from the Single File Meteorological Forcing Input File.

   = 2, wind stress and atmospheric pressure are read in at all grid nodes at a time interval that does not equal the model time step from the Single File Meteorological Forcing Input File. Interpolation in time is used to synchronize the wind and pressure information with the model time step. The wind time interval (WTIMINC) is specified below.

   = -2, wind stress and atmospheric pressure are read in at all grid nodes at a time interval that does not equal the model time step from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the Single File Meteorological Forcing Input File corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

   = 3, wind velocity is read in from a wind file from the Single File Meteorological Forcing Input File in US Navy Fleet Numeric format. This information is interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from the wind velocity. Several parameters (IREFYR, IREFMO, IREFDAY, IREFHR, IREFMIN, REFSEC, NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC) describing the Fleet Numeric wind file must be specified below.

   = 4, wind velocity and atmospheric pressure are read in (PBL/JAG format) at selected ADCIRC grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the fort.22 file dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

   =-4, wind velocity and atmospheric pressure are read in (PBL/JAG format) at selected ADCIRC grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the Single File Meteorological Forcing Input File corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

   = 5, wind velocity and atmospheric pressure are read in at all grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the Single File Meteorological Forcing Input File dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

   =-5, wind velocity and atmospheric pressure are read in at all grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the fort.22 file corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

   = 6, wind velocity and atmospheric pressure are read in for a rectangular grid (either in Longitude, Latitude or Cartesian coordinates, consistent with the grid coordinates) from the Single File Meteorological Forcing Input File. This information is interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from the wind velocity. Several parameters describing the rectangular grid and time increment (NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC) must be specified below. The meterological grid MUST cover the entire ADCIRC mesh; that is, the ADCIRC mesh must be ENTIRELY within the meteorological grid or an error will result.

   = 7, surface stress and pressure values are read in on a regular grid from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the Single File Meteorological Forcing Input File dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

   =-7, surface stress and pressure values are read in on a regular grid from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the fort.22 file corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

   = 8, hurricane parameters are read in from the Single File Meteorological Forcing Input File. Wind velocity and atmospheric pressure are calculated at every node on the fly by ADCIRC internally using the Dynamic Holland model. The input file is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format. Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format. Selecting NWS = 8 also requires the specification of the cold start time, storm number, and boundary layer adjustment (see YYYY MM DD HH24 StormNumber BLAdj below). Garret's formula is used to compute wind stress from the wind velocity.

   = 9, asymmetric hurricane model, no longer supported

   =10, wind velocity (10 m) and atmospheric pressure are read in from a sequence of National Weather Service (NWS) Aviation (AVN) model output files. Each AVN file is assumed to contain data on a Gaussian longitude, latitude grid at a single time. Consecutive files in the sequence are separated by N hours in time (where N=WTIMINC/3600 and WTIMINC is read in below). The files are named using the convention: fort.200 – wind & pressure at the time of a model hot start (this file is not used for a cold start); fort.XX1 (where XX1=200+1*N) – wind & pressure N hours after a cold or hot start; fort.XX2 (where XX2=200+2*N) – wind & pressure 2N hours after a cold or hot start; fort.XX3 (where XX3=200+3*N) – wind & pressure 3N hours after a cold or hot start and so on for all meteorological files. Prior to ADCIRC version 34.05 these files were in binary and created from a larger Grib form file using the program UNPKGRB1. Starting with ADCIRC version 34.05, the files are in ASCII tabular format. If ADCIRC is hot started, it must be done at an even N hour interval so that the hot start time corresponds to the time of a meteorological file. Enough meteorological files must be present to extend through the ending time of the model run. Garret's formula is used to compute wind stress from the wind velocity.

   =11, wind velocity (10 m) and atmospheric pressure are read in from a sequence of stripped down National Weather Service (NWS) ETA 29km model output files. Each ETA file is assumed to contain data on an E grid for a single day (8 data sets, one every 3 hours, beginning @ 03:00 and continuing through 24:00 of the given day). The files are named using the convention: fort.200 – wind & pressure the day before a model run is hot started. The final data in this file are used as the initial met condition for the hot start. This file is not used for a cold start. fort.201 – wind & pressure during the first day after a cold or hot start. fort.202 – wind & pressure during the second day after a cold or hot start fort.203 – wind & pressure during the third day after a cold or hot start. This sequence continues for all meteorological files. These files are in binary and have the format described below. The wind data is converted to an east-west, north-south coordinate system inside ADCIRC. If the model is hot started, it must be done at an even day interval so that the hot start time corresponds to the time of a meteorological file. Enough meteorological files must be present to extend through the ending time of the model run. Garret's formula is used to compute wind stress from the wind velocity.

   =12, wind velocity (10 minute averaged winds at 10m) and atmospheric pressure are provided in the OWI format on one or two rectangular (lat/lon) grid(s). If two grids are used, the first is designated as the large ("basin") scale grid, and the second is designated as the small ("region") scale grid. The Single File Meteorological Forcing Input File (fort.22) is only used to specify a few configuration parameters, while the actual wind fields are recorded in files named fort.221, fort.222, fort.223, and fort.224 (with fort.223 and fort.224 being optional). The time increment of the meteorological forcing is specified through WTIMINC in the fort.15 file. The wind and pressure fields are interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

   =15, HWind files are data assimilated snapshots of the wind velocity fields of tropical cyclones that are produced by the NOAA Hurricane Research Division (HRD). If the NWS value is set to +15, the hours column in the associated meterological forcing input file (fort.22) is relative to the cold start time. If NWS is set to -15, that hours column is relative to the hot start time. Please see the documentation of the Single File Meteorological Forcing Input File for complete details..

   =19, User has the ability to select which Isotach to use in each of the 4 quadrants. User also has ability to modify RMAX and Holland's B parameter using the ASWIP program. The auxiliary preprocessing program ASWIP.F (located in the /wind directory and executable is created by typing, make aswip, in the work folder after adcirc executable has been generated), will generate the fort.22 input file for NWS=19 from a NWS=9 formatted input file.

   Hurricane parameters are read in from the Single File Meteorological Forcing Input File. It is assumed that the line in the Single File Meteorological Forcing Input File with a zero as the forecast increment (i.e., column 6) corresponds to the start of the current simulation run, whether it is a hotstart or cold start. In other words, there is no option to set the NWS value negative to indicate that the file starts at the ADCIRC hotstart time. Rather, the forecast increment in hours (column 6) is used to indicate the relationship between the ADCIRC time and the data in the fort.22 file. Wind velocity and atmospheric pressure are calculated at exact finite element mesh node locations and directly coupled to ADCIRC at every time step using the asymmetric hurricane vortex formulation (Mattocks et al, 2006; Mattocks and Forbes, 2008) based on the Holland gradient wind model. The input file is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format. Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format. This option uses the radii at specific wind speeds (34, 50, 64, 100 knots) reported in the four quadrants (NE, SE, SW, NW) of the storm to calculate the radius of maximum winds as a function of the azimuthal angle. Garret's formula is used to compute wind stress from the wind velocity. The NWS=19 option allows the user to set a value for Rmax and Holland B Parameter. Additionally the user can select the isotachs to be used for each of the 4 quadrants. The utility program aswip_1.0.3.F located in the /wind folder will generate the NWS=19 fomatted file from a NWS=9 formatted fort.22 input file.

   In order to use the NWS=19 option, the file needs to be in best track format. The forecast period (column #6) needs to be edited to reflect the time of the forecast/nowcast for each track location (each line) in hours from the start of the simulation (0, 6, 12, 18, etc). There is no -19 option to indicate that the hours in column 6 are relative to the hotstart time. For the dynamic asymmetric model (NWS=19), ADCIRC always assumes that hour 0 corresponds to when the model is started, whether that is a cold start or a hot start. Therefore, ADCIRC analysts should not attempt to set NWS to -19.  The original data in that column depends on what type of best track format data is being used. The original data might have 0 or other numbers in that column. See: http://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
   It is suggested that users change the "BEST" tech type to "ASYM" in column 5 in the fort.22 file to denote that the file has been modified to accommodate the asymmetric wind formulation (the simulation time in hours in the 6th column has been added, etc.) so it will not get confused in the future with a best track file.

   The NWS=19 option requires the following variables in the fort.22 file in a best track format:

      o 1) Forecast time in hours (column 6); enter the time in hours in each record starting at 0

      o 2) Latitude of the eye (column 7)

      o 3) Longitude of the eye (column 8)

      o 4) Maximum sustained wind speed in knots (column 9)

      o 5) Minimum sea level pressure in MB (column 10)

      o 6) Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots) (column 12)

      o 7) Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM (columns 14, 15, 16, 17); ? 0

      o 8) Background pressure in MB (column 18); a standard value of 1013 can be used

      o 9) Rmax as reported in the ATCF BEST TRACK file in column 20

      o 10) Storm Name in Column 28 ATCF file format

      o 11) Time Record number in column 29. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File

      o 12) number of isotachs reported in the ATCF file for the corresponding Time record.

      o 13) Columns 31-34 indicate the selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed.

      o 14) Columns 35-38 are the designated Rmax values computed for each of the quadrants selected for each particular isotach.

      o 15) Column 39 is the Holland B parameter computed using the formulas outlines in the Holland paper, and implemented using the aswip program.

      The format of the file is fixed and users will want to use the aswip program to be sure that the input fort.22 file is properly formatted.

   =20, Generalized Asymmetric Holland Model (GAHM) Format. The track file format is similar to that of the older Dynamic Asymmetric Model (NWS = 19) but with 8 additional columns of data. See notes in the fort.22 file for more information. The theory and implementation of the GAHM was initially described at the 2013 ADCIRC Users Group Meeting.

   =100, 101, 102, -102, 103, 104, -104, 105, -105, 106, 110, 111, wave radiation stress is applied in addition to meteorological forcing. The meteorological input is specified by: SIGN(NWS)*(ABS(NWS)-100). For example, NWS=100 means include wave radiation stress with no meteorological forcing (NWS=0); NWS=101 means include wave radiation stress plus meteorological forcing corresponding to NWS=1; NWS=-104 means include wave radiation stress plus meteorological forcing corresponding to NWS=-4, etc. Wave radiation stress is read from a Wave Radiation Stress Forcing File. The format of this file is similar to the generic meteorological forcing file when NWS=-4 with the exception that no pressure values are read in. The time increment between consecutive radiation stress fields (RSTIMINC) is specified below.

   =300,301,302, -302, 303, 304, -304, 305, -305, 306, 310, 311, 312, -312, NWS values in the 300's indicate a SWAN+ADCIRC run Note padcswan or adcswan must be precompiled to use this option.

      The SWAN wave model is dynamically coupled to the ADCIRC model. Radiation stresses and currents from the SWAN model are applied in addition to meteorological forcing. The meteorological input is specified by: SIGN(NWS)*(ABS(NWS)-300). For example, NWS=300 means include wave radiation stress with no meteorological forcing (NWS=0); NWS=301 means include wave radiation stress plus meteorological forcing corresponding to NWS=1; NWS=-304 means include wave radiation stress plus meteorological forcing corresponding to NWS=-4, etc. Wave radiation stress are computed by the SWAN model every RSTIMINC seconds and passed into ADCIRC. In addition to assigning RSTIMINC the user must have a SWAN input and control file (fort.26) in the same working directory as the fort.15 ADCIRC control file.

.. _NRAMP:

NRAMP
   Ramp option parameter controlling whether a ramp is applied to ADCIRC forcing functions.

   = 0 No ramp function is used with forcing functions; full strength forcing is applied immediately upon cold start.

   = 1 A single hyperbolic tangent ramp function of specified duration (DRAMP, in days relative to the cold start time) will be applied to all forcing. See description of the DRAMP line for further information on the ramp function.

   = 2 Same as NRAMP=1, except that a second, separate hyperbolic tangent ramp of specified duration (DRAMPExtFlux, in days relative to cold start time plus FluxSettlingTime) specifically for external flux forcing (e.g., river boundary conditions) will also be read on the DRAMP line. In addition, the FluxSettlingTime parameter for IBTYPE=52 river boundaries will also be specified on the DRAMP line. If there are no IBTYPE=52 boundaries in the mesh (fort.14) file, the FluxSettlingTime will be read but ignored. See description of DRAMP for further information.

   = 3 Same as NRAMP=2, except that a third, separate hyperbolic tangent ramp of specified duration (DRAMPIntFlux, in days relative to cold start time plus FluxSettlingTime) specifically for internal flux forcing (e.g., flows over levees and under culverts) will also be read on the DRAMP line. See the description of the DRAMP line for further information.

   = 4 Same as NRAMP=3, except that a fourth, separate hyperbolic tangent ramp of specified duration (DRAMPElev, in days relative to cold start time plus FluxSettlingTime) specifically for elevation specified boundary forcing (e.g., tidal boundaries) will also be read on the DRAMP line. See the description of the DRAMP line for further information.

   = 5 Same as NRAMP=4, except that a fifth, separate hyperbolic tangent ramp of specified duration (DRAMPTip, in days relative to cold start time plus FluxSettlingTime) specifically for tidal potential forcing will also be read on the DRAMP line. See the description of the DRAMP line for further information.

   = 6 Same as NRAMP=5, except that a sixth, separate hyperbolic tangent ramp of specified duration (DRAMPMete, in days relative to cold start time plus FluxSettlingTime) specifically for meteorological forcing (i.e., wind and atmospheric pressure) will also be read on the DRAMP line. See the description of the DRAMP line for further information.

   = 7 Same as NRAMP=6, except that a seventh, separate hyperbolic tangent ramp of specified duration (DRAMPWRad, in days relative to cold start time plus FluxSettlingTime) specifically for wave radiation stress forcing will also be read on the DRAMP line. See the description of the DRAMP line for further information.

   = 8 Same as NRAMP=7, except that a delay parameter (DUnRampMete, in days relative to cold start time plus FluxSettlingTime) will also be read from the DRAMP line. The meteorological ramp delay parameter DUnRampMete is useful in cases where a meteorologically-forced run will be hotstarted from a long term meteorologically-free tidal spinup from cold start. The meteorological ramp delay parameter delays the start of the application of the meteorological ramp for the specified length of time, relative to the ADCIRC cold start time. See the description of the DRAMP line for further information.

.. _G:

G
   Gravitational constant. The units of this constant determine the distance units that ADCIRC operates with. (ADCIRC always operates in seconds and therefore the time units for G must be seconds.) When ICS = 2, it is required that G = 9.81 m/sec2. Regardless of ICS, when either NTIP = 1 or NCOR = 1, it is required that G = 9.81 m/sec2.

.. _TAU0:

TAU0
   Generalized Wave-Continuity Equation (GWCE) weighting factor that weights the relative contribution of the primitive and wave portions of the GWCE. If "primitive_weighting_in_continuity_equation" is specified as a nodal attribute in the fort.15 file above, this line will be read in but ignored. If a nodal attribute file is not used or "primitive_weighting_in_continuity_equation" is in the nodal attribute (fort.13) file, but not specified in the fort.15 file this TAU0 parameter will be used.

   = 0 the GWCE is a pure wave equation.

   < 1 the GWCE behaves like a pure primitive continuity equation. A good rule of thumb for setting TAU0 is to set it equal to the largest value of an equivalent linear friction factor (e.g, for linear friction TAU0 = TAU; for quadratic friction TAU0 = maximum (speed*CF/depth). Typical values for TAU0 are in the range of 0.005 – 0.1.

   = -1 the TAU0 is spatially varying but constant in time; it is calculated according to depth as follows: If the depth is >=10 TAU0 is set to 0.005, if the depth is < 10, TAU0 is set to 0.020.

   = -2 the TAU0 is spatially varying but constant in time; it is calculated according to depth as follows: if the depth is >=200 TAU0 is set to 0.005, if the depth is < 200 but > 1, then TAU0 is set to 1/depth, and if depth < 1, TAU0 is set to 1.0.

   = -3 the TAU0 varies spatially and in time. TAU0 is computed from TAU0Base read in from nodal attribute file.

      if TAU0Base < 0.025; TAU0 = TAU0Base (constant in time)

      if TAU0Base >= 0.025; TAU0 = TAU0Base + 1.5 TK(i) where TK(i)=Cd\|U\|/H

      TAU0Base values can be generated with the ADCIRC utitlity program tau0_gen.f. The program bases generation on the following logic:

         If the avg. distance between a node and its neighbors < 1750 m TAU0Base = 0.03

         If the avg. distance between a node and its neighbors > 1750 m AND depth < 10m; TAU0Base = 0.02

         If the avg. distance between a node and its neighbors > 1750 m AND depth > 10m; TAU0Base = 0.005

   = -5, FullDomainTimeVaryingTau0 = .True. the TAU0 varies spatially and in time, and is dependent on the local friction; it is limited to a range specified by Tau0FullDomainMin and Tau0FullDomainMax.

      Tau0=Tau0Min+1.5*TK(i)

   For tau0 formulations that vary spatially and temporally, ADCIRC is capable of writing out the tau0 values that it calculates internally. These values are written to a fort.90 file, which has the same format and output frequency as the water surface elevation output file (fort.63). The production of a fort.90 file is specified by placing a 1 in the tenths place of the tau0 input valuie in the fort.15 file. For example, if tau0=-3.1, the calculation of tau0 is still carried out according to the description of tau0=-3 above, and the fort.90 output file will also be produced.

.. _Tau0FullDomainMin:
.. _Tau0FullDomainMax:

Tau0FullDomainMin, Tau0FullDomainMax
   Include this line only if TAU0= -5. Specified values that the spatially and time varying TAU0 scheme must stay between. Suggested values are Tau0FullDomainMin = 0.005 and Tau0FullDomainMax = 0.2.

.. _DTDP:

DTDP
   ADCIRC time step (in seconds). Note: time in the model is computed as: TIME = STATIM*86400.+DTDP*IT.

      > 0 = The predictor-corrector algorithm is not used.

      < 0 = The predictor-corrector algorithm is used.

.. _STATIM:

STATIM
   Starting simulation time (in days). The first time step computes results at: TIME = STATIM*86400+DTDP. A nonzero value may be useful, for example, to align model output times with a specific time reference.

.. _REFTIM:

REFTIM
   Reference time (in days). This is used only to compute time for the harmonic forcing and analysis terms. A nonzero value allows equilibrium arguments to be used that have been calculated for a time other than TIME0 = STATIM*86400. The time used for harmonic terms is compute as: TIMEH = (STATIM–REFTIM)*86400.+DTDP*IT.

.. _POAN:

POAN(k)
   Parameter that weights bridge pier drag between adjacent and centerline nodes.

   = 2 if node represents a centerline node

   = 1 if node represents an adjacent node

.. _WTIMINC:

WTIMINC
   Time increment between meteorological forcing data sets (in seconds). This parameter and the line on which it appears depends on the value of the NWS parameter. See the Supplemental Meteorological/Wave/Ice Parameters table for details.

.. _YYYY:

YYYY,MM,DD,HH24,StormNumber,BLAdj
   For the Dynamic Holland model (NWS=8) and the Dynamic Asymetric Model (NWS=19), this is the coldstart datetime, the number of the storm in forecast ensemble, and the boundary layer adjustment factor. The datetime tells ADCIRC what time corresponds to t=0. For example, if the datetime is specified as 2005 08 29 06 on a cold start, then ADCIRC will find that time in the Single File Meteorological Forcing Input File, linearly interpolating if necessary, to get its initial wind state. For a hotstart, the time in the hotstart file will be added to this datetime before seeking the proper place in the Single File Meteorological Forcing Input File. For example, if the coldstart datetime is still specified as 2005 08 29 06, and the time in the hotstart file is 86400 seconds, ADCIRC will start up and interpolate in the Single File Meteorological Forcing Input File for the conditions at 6:00 am on August 30, 2005 to get its hotstart wind state. One limitation is that an ADCIRC run cannot cross the boundary of the calendar year, i.e. start in December and end in January. StormNumber is an integer and should be set to 1. BLAdj is the adjustment factor between wind speed at 10m and the wind speed at the top of the atmospheric boundary layer (winds at top of atm. b.l.)=(winds at 10m)/BLAdj. Reasonable range is 0.7 to 0.9.

.. _Geofactor:

Geofactor
   Integer that controls the form of the equation in the Generalized Asymmetric Holland Model (GAHM): (geofactor =1) the full gradient wind equation is used; (geofactor=0) a simplifying cyclostrophic balance is assumed at radius of maximum wind as done in the original Holland (1980) model derivation.

   The full gradient wind equation (geofactor = 1) is preferred, particularly for large or weak storms. Using the GAHM with geofactor = 0 should give results that are similar to the older NWS=19 parametric vortex(dynamic asymmetric vortex model) which is based on the original Holland model derivation.

.. _RSTIMINC:

RSTIMINC
   Time interval (in seconds) between successive wave radiation stress values in the Wave Radiation Stress Forcing File. This value must be specified in the Model Parameter and Periodic Boundary Condition File if the absolute value of NWS >=100.

.. _IREFYR:

IREFYR, IREFMO, IREFDAY, IREFHR, IREFMIN, REFSEC
   Starting time parameters for a Single File Meteorological Forcing Input File in US Navy Fleet Numeric format (NWS = 3, 103). These values are used in ADCIRC to compute WREFTIM which is the start time of the simulation in seconds since the beginning of the calendar year. ADCIRC is configured to accept only 1 calendar year's data, i.e., it is not possible to combine Fleet Numeric met data from two different years into a single file and then run.

   IREFYR = Year of the start of the simulation

   IREFMO = Month of the start of the simulation

   IREFDAY = Day of the start of the simulation

   IREFHR = Hour of the start of the simulation

   IREFMIN = Minute of the start of the simulation

   REFSEC
      Second of the start of the simulation.

.. _NWLAT:
.. _NWLON:
.. _WLATMAX:
.. _WLONMIN:
.. _WLATINC:
.. _WLONINC:

NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC
   Parameters describing the spatial structure of a Single File Meteorological Forcing Input File where met data is set up on a simple rectangular grid (NWS = 3, 6, 103, 106).

   NWLAT = Number of latitude values in the met file.

   NWLON = Number of longitude values in met file.

   WLATMAX = Maximum latitude (decimal deg) of data in met file (< 0 south of the equator).

   WLONMIN = Minimum longitude (decimal deg) of data in the met file (< 0 west of Greenwich meridian).

   WLATINC = Latitude increment (decimal deg) of data in the met file (must be > 0).

   WLONINC = Longitude increment (decimal deg) of data in the met file (must be > 0).

.. _RNDAY:

RNDAY
   Length of the ADCIRC run (in decimal days)

.. _DRAMP:

DRAMP
   Value (in decimal days) used to compute the ramp function that ramps up ADCIRC forcings from zero (if NRAMP=1). The ramp function is computed as RAMP=tanh(2.0*IT*DTDP/(86400.*DRAMP)) where IT = the time step number since the beginning of the model run. DRAMP is equal to the number of days when RAMP=0.96.

.. _DRAMPExtFlux:

DRAMPExtFlux
   Value (in decimal days) used to compute the ramp function that ramps up the nonzero external flux boundary condition.

.. _FluxSettlingTime:

FluxSettlingTime
   Time in days that it takes for the river flux boundary condition and the river bottom friction to equilibrate so the water surface elevation can find its steady state.

   From the start of the simulation until FluxSettlingTime has passed, the only forcing that is active is the external boundary flux forcing. All other forcings are set to zero.

   Once the FluxSettlingTime has passed, the other forcing functions begin their ramp up.

   The new IBTYPE=52 only works with periodic flux boundary conditions. Non-periodic flux boundary conditions cannot be specified for IBTYPE=52 boundaries.

.. _DRAMPIntFlux:

DRAMPIntFlux
   Value (in decimal days) used to compute the ramp function that ramps up the nonzero internal flux boundary condition.

.. _DRAMPElev:

DRAMPElev
   Value (in decimal days) used to compute the ramp function that ramps up the elevation-specified boundary condition.

.. _DRAMPTip:

DRAMPTip
   Value (in decimal days) used to compute the ramp function that ramps up the tidal potential.

.. _DRAMPMete:

DRAMPMete
   Value (in decimal days) used to compute the ramp function that ramps up the wind and atmospheric pressure.

.. _DRAMPWRad:

DRAMPWRad
   Value (in decimal days) used to compute the ramp function that ramps up the wave radiation stress.

.. _DUnRampMete:

DUnRampMete
   The meteorological ramp delay parameter with units of decimal days. It simply delays the application of the meteorological ramp for the specified length of time, relative to the ADCIRC cold start time.

.. _A00:
.. _B00:
.. _C00:

A00, B00,C00
   Time weighting factors (at time levels k+1, k, k-1, respectively) in the GWCE

.. _H0:

H0
   Minimum water depth

   If NOLIFA = 0, 1, H0 = minimum bathymetric depth. All bathymetric depths in the Grid and Boundary Information File less than H0 are changed to be equal to H0.

   If NOLIFA = 2, H0 = nominal water depth for a node (and the accompanying elements) to be considered dry (typical value 0.01 – 0.1 m).

.. _INTEGER:

INTEGER
   In the past, the wetting and drying algorithm required two additional integers as input. These extra parameters are no longer needed by the code, but they are still present to maintain backward compatibility. Their values will be ignored.

.. _VELMIN:

VELMIN
   Minimum velocity for wetting. A dry node wets if a water surface slope exists that would drive water from a currently wet node to the dry node and the steady-state current velocity that resulted would have a velocity > VELMIN. This parameter helps to keep nodes/elements from repeatedly turning on and off during the wetting process. A typical value might be 0.05 m/s.

.. _SLAM0:
.. _SFEA0:

SLAM0,SFEA0
   Longitude and latitude on which the CPP coordinate projection is centered (in degrees) if ICS = 2.

.. _TAU:

TAU
   Bottom friction is a linear function of depth-averaged velocity and TAU is the corresponding linear friction coefficient (units of 1/sec). In this case it is strongly recommended that TAU0 = TAU (Used with NOLIBF = 0). If some type of spatially varying bottom friction is specified in the NWP section, this input is ignored, and the friction coefficients will be read in from the nodal attributes file.

.. _CF:

CF
   Quadratic bottom friction coefficient (dimensionless) with the following specific meanings (note, for clarity, the quadratic friction coefficient name in the ADCIRC source code is FFACTOR):

   If NOLIBF = 1, By default, FFACTOR = CF and is spatially constant at each node in the domain. Spatially varying quadratic friction coefficients, Manning's n or Chezy coefficients may be specified using nodal attributes in the fort.13 file, (see the NWP section). In these cases the friction coefficient values specified in the nodal attribute file are converted to equivalent FFACTOR values at each node in the domain and at every timestep (e.g., see the documentation on the Manning's n nodal attribute for the formula used to convert Manning's n to an equivalent quadratic friction coefficient). If a Manning's n formulation is specified using nodal attributes, then CF is read in and used as the minimum equivalent quadratic friction coefficient (i.e., FFACTOR minimum = CF). If another bottom friction factor is specified via nodal attributes, CF is read in but ignored.

   If NOLIBF = 2, the hybrid bottom friction formulation is used FFACTOR = CF*[1+(HBREAK/H)**FTHETA]**(FGAMMA/FTHETA and CF is as specified in this expression. Note, that FFACTOR approaches CF in deep water (H > HBREAK) and the hybrid friction formulation reverts to a standard quadratic formulation. This option is not available for friction coefficients specified using nodal attributes.

.. _HBREAK:

HBREAK
   Break depth (units of length) utilized for NOLIBF = 2. If the water depth (H) is greater than HBREAK, bottom friction approaches a quadratic function of depth-averaged velocity with FFACTOR = CF. If the water depth is less than HBREAK, the friction factor increases as the depth decreases (e.g. as in a Manning-type friction law). (HBREAK = 1 m is recommended).

.. _FTHETA:

FTHETA
   Parameter (dimensionless) utilized in the hybrid bottom friction relationship (NOLIBF = 2) that determines how rapidly the hybrid bottom friction relationship approaches its deep water and shallow water limits when the water depth is greater than or less than HBREAK. (FTHETA = 10 is recommended).

.. _FGAMMA:

FGAMMA
   Parameter (dimensionless) utilized in the hybrid bottom friction relationship (NOLIBF = 2) that determines how the friction factor increases as the water depth decreases. Setting this to 1/3 gives a manning friction law type of behavior (FGAMMA = 1/3 is recommended).

.. _ESLM:

ESLM
   Spatially constant horizontal eddy viscosity for the momentum equations (units of length2/time)

.. _ESLC:

ESLC
   Spatially constant horizontal eddy diffusivity for the transport equation (units of length2/time). This is only specified if IM = 10.

.. _CORI:

CORI
   Constant Coriolis coefficient. This value is always read in, however it is only used in the computations when NCOR = 0.

.. _NTIF:

NTIF
   Number of tidal potential constituents

.. _TIPOTAG:

TIPOTAG(I)
   See description of TPK(I),AMIGT(I),ETRF(I),FFT(I) and FACET(I)

.. _TPK:
.. _AMIGT:
.. _ETRF:
.. _FFT:
.. _FACET:

TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I), I=1,NTIF
   Tidal potential amplitude, frequency, earth tide potential reduction factor (generally taken to be 0.690 for all constituents (Hendershott) but for more precise calculations can take on slightly different values (e.g. see Wahr, 1981)), nodal factor and equilibrium argument in degrees. These values are preceded by TIPOTAG(I) which is an alphanumeric descriptor (i.e. the constituent name)

.. _NBFR:

NBFR
   Number of periodic forcing frequencies on elevation specified boundaries. if NBFR=0 and a nonzero number of elevation specified boundary segments are included in the Grid and Boundary Information File, the elevation boundary condition is assumed to be non-periodic and will be read in from the Non-periodic Elevation Boundary Condition File. For reasons of backward compatability, NBFR is included in the Model Parameter and Periodic Boundary Condition File regardless of whether any elevation specified boundaries (IBTYPE=0) are defined in the fort.14 input.

.. _BOUNTAG:

BOUNTAG(k)
   See description of AMIG(k),FF(k),FACE(k)

.. _AMIG:
.. _FF:
.. _FACE:

AMIG(k),FF(k),FACE(k) k=1,NBFR
   Forcing frequency, nodal factor, equilibrium argument in degrees for tidal forcing on elevation specified boundaries. These values are preceded by BOUNTAG(k), an alphanumeric descriptor (i.e. the constituent name)

.. _ALPHAE:

ALPHAE
   See description of EMO(k,j), EFA(k,j) (<= 10 characters)

.. _EMO:
.. _EFA:

EMO(k,j),EFA(k,j) k=1,NBFR , j=1,NETA
   Amplitude and phase (in degrees) of the harmonic forcing function at the elevation specified boundaries for frequency k and elevation specified boundary forcing node j. NOTE that the parameter NETA is defined and read in from Grid and Boundary Information File: the forcing values are preceded by an alphanumeric descriptor EALPHA to facilitate verifying that the correct data matches a given frequency.

.. _ANGINN:

ANGINN
   Flow boundary nodes which are set up to have a normal flow essential boundary condition and have an inner angle less than ANGINN (specified in degrees) will have the tangential velocity zeroed. In either case, the normal velocity will be determined from the essential boundary condition.

.. _NFFR:

NFFR
   Number of frequencies in the specified normal flow external boundary condition. If NFFR=0 or NFFR=-1, the normal flow boundary condition is assumed to be non-periodic and will be read in from the Non-periodic, Normal Flux Boundary Condition File. If NFFR=0, ADCIRC assumes that the flux data in that file start at the cold start time; but if NFFR=-1, ADCIRC assumes that the flux data in that file start at the hot start time. On the other hand, positive integer values of NFFR indicate the number of frequency components that make up the periodic flux boundaries.

.. _FBOUNTAG:

FBOUNTAG(k)
   See description of FAMIGT(k),FFF(k),FFACE(k)

.. _FAMIGT:
.. _FFF:
.. _FFACE:

FAMIGT(k),FFF(k),FFACE(k) k=1,NFFR
   Forcing frequency, nodal factor, equilibrium argument in degrees for periodic normal flow forcing on flow boundaries. These values are preceded by FBOUNTAG(k), an alphanumeric descriptor (i.e. the constituent name)

.. _ALPHAQ:

ALPHAQ
   See description of QNAM(k,j),QNPH(k,j) (<=10 characters)

.. _QNAM:
.. _QNPH:

QNAM(k,j) ,QNPH(k,j) k=1,NFFR, j=1,NFLBN
   Amplitude and phase (in degrees) of the periodic normal flow/unit width (e.g. m2/s) for frequency I and "specified normal flow" boundary node j. A positive flow/unit width is into the domain and a negative flow/unit width is out of the domain. Note: the forcing values are preceded by an alphanumeric descriptor ALPHA to facilitate verifying that the correct data matches a given frequency.

.. _ENAM:
.. _ENPH:

ENAM(k,j) ,ENPH(k,j) k=1,NFFR, j=1,NFLBN
   Amplitude and phase of outgoing wave in IBTYPE=32 boundary condition (in degrees).

.. _NOUTE:
.. _TOUTSE:
.. _TOUTFE:
.. _NSPOOLE:

NOUTE, TOUTSE, TOUTFE, NSPOOLE
   NOUTE = Number of elevation recording stations.

      NOUTE =-3 Output is provided at the selected elevation recording stations in netCDF format. Following a hot start, a new fort.61.nc file is created.

      NOUTE =-2 Output is provided at the selected elevation recording stations in binary format. Following a hot start, a new fort.61 file is created.

      NOUTE =-1 Output is provided at the selected elevation recording stations in standard ascii format. Following a hot start, a new fort.61 file is created.

      NOUTE = 0 No output is provided at the selected elevation recording stations.

      NOUTE = 1 Output is provided at the selected elevation recording stations in standard ascii format. Following a hot start, continued output is merged into the existing fort.61 file.

      NOUTE = 2 Output is provided at the selected elevation recording stations in binary format. Following a hot start, continued output is merged into the existing fort.61 file.

      NOUTE = 3 Output is provided at the selected elevation recording stations in netCDF format. Following a hot start, continued output is merged into the existing fort.61.nc file.

   TOUTSE = the number of days after which elevation station data is recorded to fort.61 (TOUTSE is relative to STATIM)

   TOUTFE = the number of days after which elevation station data ceases to be recorded to fort.61 (TOUTFE is relative to STATIM)

   NSPOOLE = the number of time steps at which information is written to fort.61; i.e. the output is written to fort.61 every NSPOOLE time steps after TOUTSE

.. _NSTAE:

NSTAE
   The number of elevation recording stations (this is always read in regardless of the value of NOUTE)

.. _XEL:
.. _YEL:

XEL(k), YEL(k) k=1,NSTAE
   The coordinates of the elevation recording station k, for all NSTAE stations.

   If ICS = 1, coordinates are input as standard cartesian

   If ICS = 2, coordinates are input as degrees longitude and latitude

   If an elevation recording station is input which does not lie within the computational domain, a non-fatal error message will appear. If NFOVER has been set equal to 1, the code will estimate the nearest element and use that as the basis of interpolation. A proximity index is also printed out, which indicates how close or far the station coordinates are from the nearest element. This index may be interpreted as the number of elements that the station lies from the nearest element

.. _NOUTV:
.. _TOUTSV:
.. _TOUTFV:
.. _NSPOOLV:

NOUTV, TOUTSV,TOUTFV, NSPOOLV
   NOUTV = Output parameters which control the time series output provided for velocity solutions at selected velocity recording stations (fort.62 output)

      NOUTV =-3 Output is provided at the selected velocity recording stations in netCDF format. Following a hot start, a new fort.62.nc file is created.

      NOUTV =-2 Output is provided at the selected velocity recording stations in binary format. Following a hot start, a new fort.62 file is created.

      NOUTV =-1 Output is provided at the selected velocity recording stations in standard ascii format. Following a hot start, a new fort.62 file is created.

      NOUTV = 0 No output is provided at the selected velocity recording stations

      NOUTV = 1 Output is provided at the selected velocity recording stations in standard ascii format. Following a hot start, continued output is merged into the existing fort.62 file.

      NOUTV = 2 Output is provided at the selected velocity recording stations in binary format. Following a hot start, continued output is merged into the existing fort.62 file.

      NOUTV = 3 Output is provided at the selected velocity recording stations in netCDF format. Following a hot start, continued output is merged into the existing fort.62.nc file.

   TOUTSV = The number of days after which velocity station data is recorded to fort.62 (TOUTSV is relative to STATIM)

   TOUTFV = The number of days after which velocity station data ceases to be recorded to fort.62 (TOUTFV is relative to STATIM)

   NSPOOLV = The number of time step at which information is written to fort.62; i.e. the output is written to fort.62 every NSPOOLV time steps after TOUTSV

.. _NSTAV:

NSTAV
   The number of velocity recording stations (this is always read in regardless of the value of NOUTV)

.. _XEV:
.. _YEV:

XEV(k), YEV(k) k=1,NSTAV
   The coordinates of the velocity recording station k, for all NSTAV stations

   If ICS = 1, coordinates are input as standard cartesian

   If ICS = 2, coordinates are input as degrees longitude and latitude

   If a velocity recording station is input which does not lie within the computational domain, a non-fatal error message will appear. If NFOVER has been set equal to 1, the code will estimate the nearest element and use that as the basis of interpolation. A proximity index is also printed out, which indicates how close or far the station coordinates are from the nearest element. This index may be interpreted as the number of elements that the station lies from the nearest element

.. _NOUTC:
.. _TOUTSC:
.. _TOUTFC:
.. _NSPOOLC:

NOUTC, TOUTSC, TOUTFC, NSPOOLC
   NOUTC = Output parameters which control the time series output provided for concentration solutions at selected concentration recording stations (fort.81 output)

      NOUTC =-2 Output is provided at the selected concentration recording stations in binary format. Following a hot start, a new fort.91 file is created.

      NOUTC =-1 Output is provided at the selected concentration recording stations in standard ascii format. Following a hot start, a new fort.91 file is created.

      NOUTC = 0 no output is provided at the selected concentration recording stations

      NOUTC = 1 output is provided at the selected concentration recording stations in standard ascii format. Following a hot start, continued output is merged into the existing fort.91 file.

      NOUTC = 2 output is provided at the selected concentration recording stations in binary format. Following a hot start, continued output is merged into the existing fort.91 file.

   TOUTSC = The number of days after which concentration station data is recorded to fort.91 (TOUTSC is relative to STATIM)

   TOUTFC = The number of days after which concentration station data ceases to be recorded to fort.91 (TOUTFC is relative to STATIM)

   NSPOOLC = The number of time steps at which information is written to fort.81; i.e. the output is written to fort.81 every NSPOOLC time steps after TOUTSC

   This line is only read in if transport is included in the model run (i.e. IM=10)

.. _NSTAC:

NSTAC
   The number of concentration recording stations Note: this line is only read in if transport is included in the model run (i.e. IM=10) Note: this is read in even if NOUTC=0

.. _XEC:
.. _YEC:

XEC(k),YEC(k) k=1,NSTAC
   The coordinates of the concentration recording station k, for all NSTAC stations.

   This line is only read in if transport is included in the model run (i.e. IM=10)

   The coordinates must be consistent (i.e. cartesian or spherical) with the Grid and Boundary Information File and the coordinate designation parameter, ICS, in the Model Parameter and Periodic Boundary Condition File.

   If a concentration recording station is input which does not lie within the computational domain, a non-fatal error message will appear. If NFOVER has been set equal to 1, the code will estimate the nearest element and use that as the basis of interpolation. A proximity index is printed out in the fort.16 file that indicates how close or far the station coordinates are from the nearest element. This index may be interpreted as the number of elements that the station lies from the nearest element.
 
.. _NOUTM:
.. _TOUTSM:
.. _TOUTFM:
.. _NSPOOLM:

NOUTM, TOUTSM, TOUTFM, NSPOOLM
   NOUTM = Output parameters which control the time series output provided for met data at selected met recording stations (units 71&72 output)

      NOUTM =-3 Output is provided at the selected met recording stations in netCDF format. Following a hot start, new fort.71.nc&72.nc files are created.

      NOUTM =-2 Output is provided at the selected met recording stations in binary format. Following a hot start, new fort.71&72 files are created.

      NOUTM =-1 Output is provided at the selected met recording stations in standard ascii format. Following a hot start, new fort.71&72 files are created.

      NOUTM = 0 No output is provided at the selected met recording stations.

      NOUTM = 1 Output is provided at the selected met recording stations in standard ascii format. Following a hot start, continued output is merged into the existing fort.71&72 files.

      NOUTM = 2 Output is provided at the selected met recording stations in binary format. Following a hot start, continued output is merged into the existing fort.71&72 files.

      NOUTM = 3 Output is provided at the selected met recording stations in netCDF format. Following a hot start, continued output is merged into the existing fort.71.nc&72.nc files.

   TOUTSM = The number of days after which met station data is recorded to units 71&72 (TOUTSM is relative to STATIM)

   TOUTFM = The number of days after which met station data ceases to be recorded to units 71&72 (is relative to STATIM)

   NSPOOLM = The number of time steps at which information is written to units 71&72; i.e., output is written to units 71&72 every NSPOOLM time steps after TOUTSM. Note: this line is only read in if meteorological forcing is included in the model run (i.e. NWS<>0 and NWS<>100)

.. _NSTAM:

NSTAM
   The number of meteorological recording stations. Note: this line is only read in if met forcing is included in the model run (i.e. NWS<>0 and NWS<>100). This is read in even if NOUTM=0.

.. _XEM:
.. _YEM:

XEM(k),YEM(k) k=1,NSTAM
   The coordinates of the meteorological recording station I, for all NSTAM stations.

   This line is only read in if met forcing is included in the model run (i.e. NWS<>0 and NWS<>100)

   The coordinates must be consistent (i.e. cartesian or spherical) with the Grid and Boundary Information File and the coordinate designation parameter, ICS, in the Model Parameter and Periodic Boundary Condition File.

   If a meteorological recording station is input which does not lie within the computational domain, a non-fatal error message will appear. If NFOVER has been set equal to 1, the code will estimate the nearest element and use that as the basis of interpolation. A proximity index is printed out in the fort.16 file that indicates how close or far the station coordinates are from the nearest element. This index may be interpreted as the number of elements that the station lies from the nearest element

.. _NOUTGE:
.. _TOUTSGE:
.. _TOUTFGE:
.. _NSPOOLGE:

NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE
   NOUTGE = Output parameters which control the time series output provided for global elevation solutions at all nodes within the domain (fort.63 output)

      NOUTGE =-3 Global elevation output is provided in netCDF format. Following a hot start, a new fort.63.nc file is created.

      NOUTGE =-2 Global elevation output is provided in binary format. Following a hot start, a new fort.63 file is created.

      NOUTGE =-1 Global elevation output is provided in standard ascii format. Following a hot start, a new fort.63 file is created.

      NOUTGE = 0 No global elevation output is provided

      NOUTGE = 1 Global elevation output is provided in standard ascii format. Following a hot start, continued output is merged into the existing fort.63 file.

      NOUTGE = 2 Global elevation output is provided in binary format. Following a hot start, continued output is merged into the existing fort.63 file.

      NOUTGE = 3 Global elevation output is provided in netCDF format. Following a hot start, continued output is merged into the existing fort.63.nc file.

      NOUTGE = 4 Global elevation output is provided in sparse ascii format. Following a hot start, continued output is merged into the existing fort.63 file.

   TOUTSGE = The number of days after which global elevation data is recorded to fort.63 (TOUTSGE is relative to STATIM)

   TOUTFGE = The number of days after which global elevation data ceases to be recorded to fort.63 (TOUTFGE is relative to STATIM)

   NSPOOLGE = The number of time steps at which information is written to fort.63; i.e. the output is written to fort.63 every NSPOOLGE time steps after TOUTSGE.

.. _NOUTGV:
.. _TOUTSGV:
.. _TOUTFGV:
.. _NSPOOLGV:

NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
   Output parameters which control the time series output provided for global velocity solutions at all nodes within the domain (fort.64 output)

      NOUTGV =-3 Global velocity output is provided in netCDF format. Following a hot start, a new fort.64.nc file is created.

      NOUTGV =-2 Global velocity output is provided in binary format. Following a hot start, a new fort.64 file is created.

      NOUTGV =-1 Global velocity output is provided in standard ascii format. Following a hot start, a new fort.64 file is created.

      NOUTGV = 0 No global velocity output is provided

      NOUTGV = 1 Global velocity output is provided in standard ascii format. Following a hot start, continued output is merged into the existing fort.64 file.

      NOUTGV = 2 Global velocity output is provided in binary format. Following a hot start, continued output is merged into the existing fort.64 file.

      NOUTGV = 3 Global velocity output is provided in netCDF format. Following a hot start, continued output is merged into the existing fort.64.nc file.

      NOUTGV = 4 Global velocity output is provided in sparse ascii format. Following a hot start, continued output is merged into the existing fort.64 file.

   TOUTSGV = The number of days after which global velocity data is recorded to fort.64 (TOUTSGV is relative to STATIM)

   TOUTFGV = The number of days after which global velocity data ceases to be recorded to fort.64 (TOUTFGV is relative to STATIM)

   NSPOOLGV = The number of time steps at which information is written to fort.64; i.e. the output is written to fort.64 every NSPOOLGV time steps after TOUTSGV

.. _NOUTGC:
.. _TOUTSGC:
.. _TOUTFGC:
.. _NSPOOLGC:

NOUTGC, TOUTSGC, TOUTFGC, NSPOOLGC
   NOUTGC = Output parameters which control the time series output provided for global concentration solutions at all nodes within the domain (fort.93 output)

      NOUTGC =-2 Global concentration output is provided in binary format. Following a hot start, a new fort.93 file is created.

      NOUTGC =-1 Global concentration output is provided in standard ascii format. Following a hot start, a new fort.93 file is created.

      NOUTGC = 0 No global concentration output is provided

      NOUTGC = 1 Global concentration output is provided in standard ascii format. Following a hot start, continued output is merged into the existing fort.93 file.

      NOUTGC = 2 Global concentration output is provided in binary format. Following a hot start, continued output is merged into the existing fort. 93 file.

      NOUTGC = 3 Global concentration output is provided in sparse ascii format. Following a hot start, continued output is merged into the existing fort.93 file.

   TOUTSGC = The number of days after which global concentration data is recorded to fort.93 (TOUTSGC is relative to STATIM)

   TOUTFGC = The number of days after which global concentration data ceases to be recorded to fort.93 (TOUTFGC is relative to STATIM)

   NSPOOLGC = The number of time steps at which information is written to fort.93; i.e. the output is written to fort.93 every NSPOOLGC time steps after TOUTSGC

   This line is only read in if transport is included in the model run (i.e. IM=10)

.. _NOUTGW:
.. _TOUTSGW:
.. _TOUTFGW:
.. _NSPOOLGW:

NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW
   NOUTGW = Output parameters which control the time series output provided for wind stress or velocity and atmospheric pressure at all nodes within the domain (fort.73 and 74 output)

      NOUTGW =-3 Global wind stress/velocity and atmospheric pressure outputs are provided in netCDF format. Following a hot start, new fort.73.nc and 74.nc files are created.

      NOUTGW =-2 Global wind stress/velocity and atmospheric pressure outputs are provided in binary format. Following a hot start, new fort.73 and 74 files are created.

      NOUTGW =-1 Global wind stress/velocity and atmospheric pressure outputs are provided in standard ascii format. Following a hot start, new fort.73 and 74 files. are created.

      NOUTGW = 0 no global wind stress/velocity or atmospheric pressure output is provided

      NOUTGW = 1 Global wind stress/velocity and atmospheric pressure output are provided in standard ascii format. Following a hot start, continued output is merged into the existing fort.73 and 74 files.

      NOUTGW = 2 Global wind stress/velocity and atmospheric pressure output are provided in binary format. Following a hot start, continued output is merged into the existing fort.73 and 74 files.

      NOUTGW = 3 Global wind stress/velocity and atmospheric pressure output are provided in netCDF format. Following a hot start, continued output is merged into the existing fort.73.nc and 74.nc files.

      NOUTGW = 4 Global wind stress/velocity and atmospheric pressure output are provided in sparse ascii format. Following a hot start, continued output is merged into the existing fort.73 and 74 files.

   TOUTSGW = The number of days after which global wind stress/velocity and atmospheric pressure data are recorded to units 73,74 (TOUTSGW is relative to STATIM)

   TOUTFGW = The number of days after which global wind stress/velocity and atmospheric pressure data cease to be recorded to units 73,74 (TOUTFGW is relative to STATIM)

   NSPOOLGW = The number of time steps at which information is written to units 73,74; i.e. the output is written to units 73,74 every NSPOOLGW time steps after TOUTSGW

   This line is only read in if meteorological forcing is included in the model run (i.e. NWS<>0 and NWS<>100)

.. _NFREQ:

NFREQ
   Number of frequencies included in harmonic analysis of model results. Note: harmonic output is only available for 2DDI elevation and velocity

.. _NAMEFR:

NAMEFR(k)
   An alphanumeric descriptor (i.e. the constituent name) whose length must be <= 10 characters

.. _HAFREQ:
.. _HAFF:
.. _HAFACE:

HAFREQ(k), HAFF(k), HAFACE(k) k=1,NFREQ
   Parameters describing the constituents to be included in the harmonic analysis of model results

   HAFREQ(k) = frequency (rad/s)

   HAFF(k) = nodal factor

   HAFACE(k) = equilibrium argument (degrees)

   If a steady component will be included in the harmonic analysis, this must be the first constituent listed (i.e., the constituent corresponding to k=1) 

.. _THAS:
.. _THAF:
.. _NHAINC:
.. _FMV:

THAS, THAF, NHAINC, FMV
   Parameters that control the calculation of harmonic constituents both at stations and globally

   THAS = The number of days after which data starts to be harmonically analysed (THAS is relative to STATIM)

   THAF = The number of days after which data ceases to be harmonically analysed (THAF is relative to STATIM)

   NHAINC = The number of time steps at which information is harmonically analysed (information every NHAINC time steps after THAS is used in harmonic analysis)

   FMV = Fraction of the harmonic analysis period (extending back from the end of the harmonic analysis period) to use for comparing the water elevation and velocity means and variances from the raw model time series with corresponding means and variances of a time series resynthesized from the harmonic constituents. This comparison is helpful for identifying numerical instabilities and for determining how complete the harmonic analysis was. Examples:

      FMV = 0. Do not compute any means and vars.

      FMV = 0.1 Compute means and vars. over final 10% of period used in harmonic analysis

      FMV = 1.0 Compute means and vars. over entire period used in harmonic analysis

      Note: the means and variance calculations are only done if global harmonic calculations are performed. Results are written out to fort.55. A summary of the poorest comparisons throughout the domain and the node numbers where these occurred is given at the end of the fort.16 output file.

      Note: the time series resysthesis from the harmonic constituents can use up a lot of CPU time since this is done for every time step during the specified part of the harmonic analysis period. If the harmonic analysis period extends for only a few days, it is practical to set FMV=1. Otherwise, it becomes unreasonably time consuming to compute means and variances for more than 10-20 days. Ultimately, the practical limit to these calculations depends on the number of nodes, the number of constituents in the harmonic analysis, and the size of the time step.

.. _NHASE:
.. _NHASV:
.. _NHAGE:
.. _NHAGV:

NHASE, NHASV, NHAGE, NHAGV
   Parameters that control the spatial locations where harmonic analysis is performed

   NHASE = 0 no harmonic analysis is performed at the selected elevation recording stations

   NHASE = 1 harmonic analysis is performed at the selected elevation recording stations (output on fort.51) Note: the stations are as specified in the section on time series station elevation output

   NHASV = 0 no harmonic analysis is performed at the selected velocity recording stations

   NHASV = 1 harmonic analysis is performed at the selected velocity recording stations (output on fort.52) Note: the stations are as specified in the section on time series station velocity output

   NHAGE = 0 no harmonic analysis is performed for global elevations

   NHAGE = 1 harmonic analysis is performed for global elevations (output on fort.53)

   NHAGV = 0 no harmonic analysis is performed for global velocities

   NHAGV = 1 harmonic analysis is performed for global velocities (output on fort.54)

.. _NHSTAR:
.. _NHSINC:

NHSTAR, NHSINC
   NHSTAR = Parameters that control the generation of hot start output.

      NHSTAR = 0 no hot start output files generated

      NHSTAR = 1 hot start output files generated in binary format

      NHSTAR = 3 hot start output files generated in netCDF format

      NHSTAR = 5 hot start output files generated in netCDF4 format

   NHSINC = The number of time steps at which hot start output file is generated (hot start file is generated every NHSINC time steps). The time step increments are always counted starting from the cold start time.

.. _ITITER:
.. _ISLDIA:
.. _CONVCR:
.. _ITMAX:

ITITER, ISLDIA, CONVCR, ITMAX
   ITITER = Parameters that provide information about the solver that will be used for the GWCE.

      ITITER = -1 only for lumped, explicit GWCE, matrix is diagonal and no external solver is needed

      ITITER = 1 use iterative JCG solver (from ITPACKV 2D)

   ISLDIA = Parameters that provide information about the solver that will be used for the GWCE.

      ISLDIA = 0 fatal error messgs only from ITPACKV 2D(fort.33)

      ISLDIA = 1 warning messgs and minimum output from ITPACKV 2D (fort.33)

      ISLDIA = 2 reasonable summary of algorithm progress from ITPACKV 2D (fort.33)

      ISLDIA = 3 parameter values and informative comments from ITPACKV 2D (fort.33)

      ISLDIA = 4 approximate solution after each iteration from ITPACKV 2D (fort.33)

      ISLDIA = 5 original system from ITPACKV 2D (fort.33)

   CONVCR = Absolute convergence criteria (should be no smaller than 500 times the machine precision)

   ITMAX = Maximum number of iterations each time step

      Note: all of the parameters must be input regardless of whether a diagonal or iterative solver is selected. However, ISLDIA, CONVCR and ITMAX are only used with the iterative solvers
      
      CONVCR and ITMAX are only used with the iterative solvers

      Note: we typically use CONVCR=1E-10. After the first few time steps, the solutions usually converge within 5-10 iterations.

.. _ISLIP:

ISLIP = 3D bottom friction code

   ISLIP = 0, no slip bottom b.c.

   ISLIP = 1, linear slip bottom b.c.

   ISLIP = 2, quadratic slip bottom b.c. where the quadratic slip coefficient is computed using Log Layer formula
   .. image:: ../_static/log_layer_formula_islip_2.avif

   ISLIP = 3, quadratic slip b.c.

.. _KP:

KP
   3D bottom friction coefficient used in ADCIRC

   If ISLIP = 0 the bottom friction coefficient is ignored

   If ISLIP = 1 the bottom friction is a linear function of bottom velocity and KP is the corresponding linear friction coefficient (units of velocity)

   If ISLIP = 2 bottom friction is is computed using the log layer formula and KP is the minimum quadratic bottom friction coefficient (dimensionless)

   If ISLIP = 3 bottom friction is a quadratic function of bottom velocity and KP is the corresponding quadratic friction coefficient (dimensionless)

.. _Z0S:
.. _Z0B:

Z0S, Z0B
   Free surface & bottom roughnesses

   For IEVC=50, Z0S is the spatially constant value. For IEVC=51, the surface roughness length is computed dynamically (seeIEVC=51), and Z0S is the minimum surface roughness length.

   For nodal attributes bottom_roughness_length and mannings_n_at_the_sea_floor, the bottom roughness length is read in from the nodal attribute file either directly or as a Manning's n roughness. Currently Z0B is not used for either of these two cases.

   If a Manning's n roughness is read in, the roughness length is expressed in terms of the water depth H and the Manning's n:

   .. image:: ../_static/z0s_z0b_eqn.avif

   where K = 0.4 is the von Karman constant, and g is the gravitational acceleration (Bretschneider et al., 1986). New roughness lengths are computed at each time step, based on the computed water depth and Manning's n value at each mesh vertex.

.. _ALP1:
.. _ALP2:
.. _ALP3:

ALP1, ALP2, ALP3
   Time weighting coefficients for the 3D velocity solution.

   0.= fully explicit, 0.5=time centered, 1.= fully implicit

   ALP1 weights the Coriolis term

   ALP2 weights the bottom friction terms

   ALP3 weights the vertical diffusion terms

.. _IGC:
.. _NFEN:

IGC, NFEN
   Vertical grid code, # nodes in the vertical grid

   IGC = 0, vertical grid read in

   IGC = 1, uniform vertical grid generated

   IGC = 2, log vertical grid generated

   IGC = 3, log linear vertical grid generated

   IGC = 4, double log vertical grid generated

   IGC = 5, P-grid generated

   IGC = 6, sine grid generated


.. _IEVC:
.. _EVMIN:
.. _EVCON:

IEVC, EVMIN, EVCON
   Vertical eddy viscosity code, vertical eddy viscosity minimum value and vertical eddy viscosity constant

   NOTE: EVCON is only used for some of the vertical eddy viscosity formulations as discussed below.

   NOTE: In cases where vertical eddy viscosity is specified to vary linearly over the lower 20% of the water column, it actually varies linearly with a constant slope up to the vertical FE grid node that is less than or equal to the 20% location. The value is constant as specified at all FE grid nodes above the 20% location. The vertical eddy viscosity above and below the 20% level is joined by one additional linearly varying segment.

   NOTE: The vertical eddy viscosity is constrained to always be greater than or equal to EVMIN.

   IEVC = 0-9, EV constant in time & horizontal space

      0 – vertical eddy viscosity read in – EVCON is not used

      1 – EV = EVCON

   IEVC =10-19, Vertical eddy viscosity proportional to omega*h*h (Lynch and Officer (1986) Lynch and Werner (1987, 1991))

      10 – EV = omega*h*h/10 over the entire water column

      11 – EV = omega*h*h/1000 at bottom varies linear over lower 20% of water column

              = omega*h*h/10 in upper 80% of water column

      NOTE:For this vertical eddy viscosity formulation, EVCON is not used and omega is hardwired for a 12.42 hour tide.

   IEVC =20-29, EV proportional to kappa U* z

      20 – EV = 0.41U*Zo at bottom

              = 0.41U*Z over entire water column

      21 – EV = 0.41U*Zo at bottom

              = 0.41U*Z in lower 20% of water column

              = 0.082U*h in upper 80% of water column

      WHERE: U* is the friction velocity

      NOTE: For this EV formulation, EVCON is not used.

   IEVC =30-39, EV proportional to Uh (Davies 1990)

      30 – EV = 0.025\|U\|h/9.001 over entire water column

   31 – EV = EVCON \|U\|h over entire water column

   32 – EV = 0.025\|U\|h/9.001 in upper 80% of water column

           = 0.000025h\|U\|/9.001 at bottom varies linear over lower 20% of water column

   33 – EV = EVCON \|U\|h in upper 80% of water column

           = EVCON \|U\|h/1000. at bottom varies linear over lower 20% of water column

      WHERE: U is depth averaged velocity

      NOTE: For this vertical eddy viscosity formulation, EVCON is used only for IEVC =31,33

   IEVC = 40-49, EV proportional to U*U (Davies 1990)

      40 – EV = 2\|UU\|/9.001 over entire water column

      41 – EV = EVCON \|UU\| over entire water column

      42 – EV = 2\|UU\|/9.001 in upper 80% of water column

             = 0.002\|UU\|/9.001 at bottom varies linear over lower 20% of water column

      43 – EV = EVCON \|UU\| in upper 80% of water column

              = EVCON \|UU\|/1000. at bottom varies linear over lower 20% of water column

      WHERE: U is depth averaged velocity

      NOTE: For this EV formulation, EVCON is used only for IEVC=41,43

   IEVC =50, EV computed from Mellor-Yamada L2.5 closure. NOTE: For this EV formulation EVCON is not used.

      NOTE: For this EV formulation, EVCON is not used.

   IEVC =51, EV computed from Mellor-Yamada L2.5 closure with parameterizations to include enhanced mixing in the surface layer.

      .. image:: ../_static/ievc_51_note.avif

      .. image:: ../_static/ievc_51_eqns.avif
         
.. _EVTOT:

EVTOT(K)
   Eddy viscosity associated with vertical grid node K

.. _THETA1:
.. _THETA2:

THETA1, THETA2
   Time weighting coefficients for the MY2.5 turbulence soln. (include this line only if IEVC = 50)

   0.= fully explicit, 0.5=time centered, 1.= fully implicit

   THETA1 weights the dissipation term

   THETA2 weights the vertical diffusion term

.. _I3DSD:

I3DSD
   I3DSD = 0 no station 3D temperature, salinity and/or density info is output to unit 41

   I3DSD = 1 station 3D temperature, salinity and/or density info is output in ascii format

   I3DSD = 2 station 3D temperature, salinity and/or density info is output in binary format

.. _TO3DSDS:

TO3DSDS
   The number of days after which station 3D temperature, salinity and/or density are written to unit 41.

.. _TO3DSDF:

TO3DSDF
   The number of days after which station 3D temperature, salinity and/or density cease to be written to unit 41.

.. _NSPO3DSD:

NSPO3DSD
   The number of time steps at which data is written to unit 41. (i.e., data is output to unit 41 every NSPO3DSD time steps after TO3DSSD.)

.. _NSTA3DD:

NSTA3DD
   Number of 3D density stations

.. _X3DS:
.. _Y3DS:

X3DS(k), Y3DS(k)
   The coordinates of the 3D temperature, salinity, and/or density recording station k, for all NSTA3DD, NSTA3DV or NSTA3DT stations (only include this line if I3DSD, I3DSV or I3DST is not = 0)

.. _I3DSV:

I3DSV
   I3DSV = 0 no station 3D velocities are output to unit 42

   I3DSV = 1 station 3D velocities are output in ascii forma

   I3DSV = 2 station 3D velocities are output in binary format

.. _TO3DSVS:

TO3DSVS
   The number of days after which station 3 D velocities are written to unit 42.

.. _TO3DFVF:

TO3DFVF
   The number of days after which station 3 D velocities cease to be written to unit 42.

.. _NSPO3DSV:

NSPO3DSV
   The number of time steps at which data is written to unit 42. (i.e., data is output to unit 42 every NSPO3DSV time steps after TO3DSSV.)

.. _NSTA3DV:

NSTA3DV
   Number of 3D velocity stations

.. _I3DST:

I3DST
   I3DST = 0 no station 3D turbulence variables output to unit 43

   I3DST = 1 station 3D turbulence variables output in ascii format

   I3DST = 2 station 3D turbulence variables output in binary format

.. _TO3DSTS:

TO3DSTS
   The number of days after which station 3D turbulence variables are written to unit 43.

.. _TO3DSTF:

TO3DSTF
   The number of days after which station 3D turbulence variables cease to be written to unit 43.

.. _NSPO3DST:

NSPO3DST
   The number of time steps at which data is written to unit 43. (i.e., data is output to unit 43 every NSPO3DSV time steps after TO3DSSV.)

.. _NSTA3DT:

NSTA3DT
   Number of 3D turbulence stations

.. _I3DGD:

I3DGD
   I3DGD = 0 no global 3D temperature, salinity, and/or density info is output to unit 44

   I3DGD = 1 global 3D temperature, salinity, and/or density info is output in ascii format

   I3DGD = 2 global 3D temperature, salinity, and/or density info is output in binary format

.. _TO3DGDS:

TO3DGDS
   The number of days after which global 3D temperature, salinity, and/or density are written to unit 44.

.. _TO3DGDF:

TO3DGDF
   The number of days after which global 3D temperature, salinity, and/or density cease to be written to unit 44.

.. _NSPO3DGD:

NSPO3DGD
   The number of time steps at which data is written to unit 44. (i.e., data is output to unit 44 every NSPO3DGD time steps after TO3DSGD.)

.. _I3DGV:

I3DGV
   I3DGV = 0 no global 3D velocities are output to unit 45

   I3DGV = 1 global 3D velocities are output in ascii format

   I3DGV = 2 global 3D velocities are output in binary format

.. _TO3DGVS:

TO3DGVS
   The number of days after which global 3D velocity data is written to unit 45.

.. _TO3DGVF:

TO3DGVF
   The number of days after which global 3D velocity data ceases to be written to unit 45.

.. _NSPO3DGV:

NSPO3DGV
   The number of time steps at which data is written to unit 45. (i.e., data is output to unit 45 every NSPO3DGV time steps after TO3DSGV.)

.. _I3DGT:

I3DGT
   I3DGT = 0 no global 3D turbulence variables output to unit 46

   I3DGT = 1 global 3D turbulence variables output in ascii format

   I3DGT = 2 global 3D turbulence variables output in binary format

.. _TO3DGTS:

TO3DGTS
   The number of days after which global 3D turbulence variables are written to unit 46.

.. _TO3DGTF:

TO3DGTF
   The number of days after which global 3D turbulence variables cease to be written to unit 46.

.. _NSPO3DGT:

NSPO3DGT
   The number of time steps at which data is written to unit 46. (i.e., data is output to unit 46 every NSPO3DGT time steps after TO3DSGT.)

.. _RES_BC_FLAG:

RES_BC_FLAG
   Controls the type of boundary conditions used in the 3D baroclinic simulations. Must be the same as IDEN.

   RES_BC_FLAG < 0 Diagnostic simulations, so the only boundary condition utilized is the levels of no motion (steric adjustments to the elevations – fort.35). In these cases, RBCTIMEINC (levels of no motion boundary time interval) and BCSTATIM (starting time for the level of no motion boundary condition information) are needed.

   RES_BC_FLAG = 1 Not implemented (not a valid value).

   RES_BC_FLAG = 2 Prognostic simulation using salinity field, so the boundary conditions utilize both the levels of no motion (steric adjustments to the elevations found in fort.35) and the salinity field of the outside ocean (found in fort.36). In this case, RBCTIMEINC, SBCTIMEINC (levels of no motion and salinity boundary time interval), BCSTATIM and SBCSTATIM (starting time for the level of no motion and salinity boundary condition information) are needed.

   RES_BC_FLAG = 3 Prognostic simulation using temperature field, so the boundary conditions utilize both the levels of no motion (steric adjustments to the elevations found in fort.35) and the temperature field of the outside ocean (found in fort.37). In this case, RBCTIMEINC, TBCTIMEINC (levels of no motion and temperature boundary time interval), BCSTATIM and TBCSTATIM (starting time for the level of no motion and temperature boundary condition information) are needed.

   RES_BC_FLAG = 4 Prognostic simulation using both salinity and temperature fields, so the boundary conditions utilize the levels of no motion (steric adjustments to the elevations found in fort.35), and the salinity field and the temperature field of the outside ocean (found in fort.36 and fort.37, respectively). In this case, RBCTIMEINC, SBCTIMEINC, TBCTIMEINC, BCSTATIM, SBCSTATIM and TBCSTATIMare needed.


.. _BCFLAG_LNM:

BCFLAG_LNM
   Flag to control the levels of no motion parameterization used in the 3D ADCIRC run.

   BCFLAG_LNM = 1 Levels of no motion are elevations provided from a file (fort.35) and that file is accessed based on RBCTIMEINC.

   BCFLAG_LNM = 2 Not a valid value.

   BCFLAG_LNM = 3 Levels of no motion are elevations provided from the initial values of the elevations. These values can be obtained from either the initial condition file (fort.17) or the hotstart file (fort.67 or fort.68). This option should be used more for a diagnostic simulation than for the prognostic simulation.

.. _RBCTIMEINC:

RBCTIMEINC
   Time interval between data sets for the level of no motion boundary condition data, in seconds.

.. _BCSTATIM:

BCSTATIM
   Starting time (in seconds since ADCIRC cold start) for boundary condition data for the level of no motion boundary condition.

.. _SBCTIMEINC:

SBCTIMEINC
   Time interval between data sets for the salinity boundary condition data, in seconds.

.. _SBCSTATIM:

SBCSTATIM
   Starting time (in seconds since ADCIRC cold start) for boundary condition data for the salinity boundary condition.

.. _TBCTIMEINC:

TBCTIMEINC
   Time interval between data sets for the temperature boundary condition data, in seconds.

.. _TBCSTATIM:

TBCSTATIM
   Starting time (in seconds since ADCIRC cold start) for boundary condition data for the temperature boundary condition.

.. _TTBCTIMEINC:

TTBCTIMEINC
   Time interval between data sets for the surface heat flux boundary condition data, in seconds.

.. _TTBCSTATIM:

TTBCSTATIM
   Starting time (in seconds since ADCIRC cold start) for boundary condition data for the surface heat flux boundary condition.

.. _BCFLAG_TEMP:

BCFLAG_TEMP
   Controls the surface heat flux parameterization used in the 3D ADCIRC simulation with temperature. ADCIRC ignores this parameter unless RES_BC_FLAG is 3 or 4. Currently, there are three options for the surface heat flux boundary (all surface heat flux values are read in from a fort.38 file).

   BCFLAG_TEMP = 1 Surface heat flux values are directly provided from a file (fort.38), and that file is accessed based on TTBCTIMEINC (surface heat flux boundary time interval) and starts reading the values based on TTBCSTATIM (starting time for the surface heat flux boundary condition information).

   BCFLAG_TEMP = 2 Surface heat flux values are calculated from information obtained from a file (fort.38), and that file is accessed based on TTBCTIMEINC and starts reading the values based on TTBCSTATIM. This option looks for 6 values (sensible heat flux, latent heat flux, downward shortwave radiation, downward longwave radiation, upward shortwave radiation, upward longwave radiation) to use in calculating the surface heat flux from the following equation (Mellor, 1996): qh = qs + ql + (Sraddown – Sradup) + (lraddown – lradup) where qs = sensible heat flux, ql = latent heat flux, sraddown=downward shortwave radiation, sradup=upward shortwave radiation, lraddown=downward longwave radiation, lradup=upward longwave radiation.

   BCFLAG_TEMP = 3 Surface heat flux values are calculated from information obtained from a file (fort.38), and that file is accessed based on TTBCTIMEINC and starts reading the values based on TTBCSTATIM. This option looks for 4 values (net shortwave radiation, net longwave radiation, latent heat flux, sensible heat flux,) to use in calculating the surface heat flux from the following equation (Mellor, 1996): qh = qs + ql + Sradnet + lradnet where qs = sensible heat flux, ql = latent heat flux, sradnet=net shortwave radiation, lradnet=net longwave radiation.

.. _SPONGEDIST:

SPONGEDIST
   Controls the sponge layer that allows for a spatial ramp to occur on the wind and advection terms, but it does not vary in time. It starts at zero on the boundary and linearly increases to one over a user-defined distance into the domain. The distance must be given in either m or ft, depending on the units on gravity. The sponge layer is utilized for both the diagnostic and prognostic simulations for the wind terms and all advective terms (momentum and transport equations, if applicable). If SPONGEDIST is not equal to zero, and NOUTGE is not equal to zero, ADCIRC will produce a fulldomain output file that shows the extent of the sponge layer (similar in format to the maxele.63 file) called fort.92.

.. _EQNSTATE:

EQNSTATE
   Indicates the equation of state used to convert the salinity and temperature values into density values.

   EQNSTATE = 1 Use the equation of state given in Cushman-Roisin, B., Introduction to Geophysical Fluid Dynamics, Prentice-Hall, 1994, 320 pp. (1994): rho = rho0 (1 – alpha(T-T0)+beta(S-S0)) where rho0=1028 kg/m3 is the reference density of seawater, alpha=0.00017/degC is the coefficient of thermal expansion, and beta=0.00076 is the coefficient of saline concentration. The reference values for temperature and salinity are: T0=10 degC and S0=35 psu. This option can be used for IDEN=2, -2, 3, -3, 4 or -4.

   EQNSTATE = 2 Use the equation of state given in McDougall, T. J., D.G. Wright, D. R. Jackett and R. Feistel, "Accurate and computationally efficient algorithms for potential temperature and density of seawater", Journal of Atmospheric and Oceanic Technology, 20 (5), 2003, pp. 730-741. This uses the temperature, salinity and pressure in determining the density field. This option can be used when IDEN = 4 or -4. The equation used with this option can be found in the theory manual.

   EQNSTATE = 3 Use the equation of state is the UNESCO equation and is given in Gill, A.E., Atmosphere-Ocean Dynamics, Academic Press, 1982, 662 pp. and Mellor, G.L., Introduction to Physical Oceanography, American Institute of Physics, 1996, 284 pp. It uses the temperature, salinity and pressure in determining the density field.
   
   This option can be used when IDEN = 4 or -4. The equation used with this option can be found in the theory manual.

.. _NLSD:

NLSD
   Lateral salinity diffusion coefficient.

.. _NVSD:

NVSD
   Vertical salinity diffusion coefficient. If IEVC=50, this coefficient is calculated by the Mellor-Yamada equations and this value is ignored.

.. _NLTD:

NLTD
   Lateral temperature diffusion coefficient.

.. _NVTD:

NVTD
   Vertical temperature diffusion coefficient. If IEVC=50, this coefficient is calculated by the Mellor-Yamada equations and this value is ignored.

.. _ALP4:

ALP4
   Time stepping coefficient associated with the transport equation terms.

fort.12
~~~~~~~

.. _AGRID2:

AGRID2
   Alphanumeric file identification (<=24 characters). To facilitate organization of files for individual model runs, it is suggested that AGRID2 match AGRID in the Grid and Boundary Information File.

.. _STARTDRY:

STARTDRY(JN)
   Start dry code; STARTDRY = -88888 at nodes which will be initialized as dry. It can have any value at other nodes.

.. _DUM1:
.. _DUM2:

DUM1, DUM2
   Dummy variables

.. _fort.10:

fort.10
~~~~~~~

.. _DACONC:

DACONC
   Generic passive scalar 2D depth-averaged concentration field

.. _CONC:

CONC
   Generic passive scalar 3D concentration field

fort.11
~~~~~~~

.. _NVN:

NVN
   Number of nodes in vertical, must match NFEN

.. _NVP:

NVP
   Number of nodes in the horizontal grid, must match NP

.. _jki:

jki
   Node number

.. _NHNN:

NHNN
   Horizontal node number

.. _NVNN:

NVNN
   Vertical node number

.. _DASIGT:

DASIGT
   Sigma T value (kg/m^3) (=density-1000) for a 2DDI run

.. _DATEMP:

DATEMP
   Temperature (DEG C) for a 2DDI run

.. _DASALT:

DASALT
   Salinity (PSU) for a 2DDI run

.. _SIGT:

SIGT(NHNN,NVNN)
   Sigma T value (kg/m^3) (=density-1000)

.. _TEMP:

TEMP(NHNN,NVNN)
   Temperature (DEG C)

.. _SAL:

SAL(NHNN,NVNN)
   Salinity (PSU)

.. _SIGMA:

SIGMA(K)
   Dimensionless level of the vertical grid node K from -1 (bottom) to +1 (surface)

fort.13
~~~~~~~

.. _NumOfNodes:

NumOfNodes
   Number of nodes, must match NP from grid file

.. _NAttr:

NAttr
   Number of attributes contained in the fort.13 file, must be equal to or greater than NWP from fort.15 file

.. _AttrName:

AttrName(i)
   Nodal attribute name, the fort.13 file must contain data for all AttrNames that appear in the fort.15 file. Valid names followed by the ADCIRC variable are:

   primitive_weighting_in_continuity_equation – Tau0

   surface_submergence_state – StartDry

   quadratic_friction_coefficient_at_sea_floor – Fric

   surface_directional_effective_roughness_length – z0Land (Note: this attribute has ValuesPerNode = 12)

   surface_canopy_coefficient – VCanopy

   bridge_pilings_friction_paramenters – BK, BAlpha, BDelX, POAN (Note: this attribute has ValuesPerNode = 4)

   mannings_n_at_sea_floor – ManningsN

   chezy_friction_coefficient_at_sea_floor – ChezyFric

   sea_surface_height_above_geoid – GeoidOffset

   bottom_roughness_length – Z0b_var

   wave_refraction_in_swan – SwanWaveRefrac

   average_horizontal_eddy_viscosity_in_sea_water_wrt_depth – ESLM

   elemental_slope_limiter – elemental_slope_limiter_grad_max

   advection_state – AdvectionState

   initial_river_elevation – Eta2

   Note: if the user selects quadratic_friction_coefficient_at_sea_floor, mannings_n_at_sea_floor, or chezy_friction_coefficient_at_sea_floor, then NOLIBF must be 1 (nonlinear friction formulation) since all those formulations are nonlinear. If the NOLIBF were anything other than 1, it is an error that will cause ADCIRC to stop.

.. _Units:

Units(i)
   Physical units (ft, m/s, 1 = unitless)

.. _ValuesPerNode:

ValuesPerNode(i)
   Number of values at each node for a particular attribute

.. _DefaultAttrVal:

DefaultAttrVal(i,k)
   Default value(s) for the nodal attribute

.. _NumNodesNotDefaultVal:

NumNodesNotDefaultVal(i)
   Number of nodes with non-default values

.. _n:

n
   Node number

.. _AttrVal:

AttrVal(n,k)
   Nodal attribute value(s)

fort.19
~~~~~~~

.. _ESBIN:

ESBIN(k)
   Elevation (referenced to the geoid) at specified elevation node k. The sequencing is assumed to match what is defined in the elevation specified boundary condition part of the Grid and Boundary Information File.

.. _ETIMINC:

ETIMINC
   Time increment (secs) between consecutive sets of elevation specified boundary condition values contained in this file. 

fort.20
~~~~~~~

.. _FTIMINC:

FTIMINC
   Time increment (secs) between consecutive sets of normal flow boundary condition values contained in this file.

.. _NFLBN:

NFLBN
   Total number of flow boundary nodes

.. _QNIN:

QNIN(k)
   Normal flow/unit width (e.g., m2/s) at specified normal flow node k. A positive flow/unit width is into the domain and a negative flow/unit width is out of the domain. The sequencing is assumed to match what is defined in the part of the Grid and Boundary Information File specifying non zero normal flow boundaries.

.. _FRIC:

FRIC(k)
   Nodal bottom friction coefficient.
 

fort.22
~~~~~~~

.. _IWTIME:

IWTIME
   Time of the wind field in the following integer format: YEAR*1000000 + MONTH*10000 + DAY*100 + HR

.. _PRN:

PRN
   Applied atmospheric pressure at the free surface. Units depend on the specific type of wind input file.

.. _WDIR:

WDIR
   (If NWS = 3, 103) direction wind blows from in deg cw from north

.. _WSPEED:

WSPEED
   (If NWS =3, 103) wind speed in m/s

.. _WSX:

WSX
   Applied horizontal wind stress in the x,y directions divided by the reference density of water (should be units (length/time)2). An oceanographic convention is used where velocity is positive when it is blowing toward positive coordinate directions.

.. _WSY:

WSY
   Applied horizontal wind stress in the x,y directions divided by the reference density of water (should be units (length/time)2). An oceanographic convention is used where velocity is positive when it is blowing toward positive coordinate directions.

.. _WVX:
.. _WVY:

WVX, WVY
   Applied horizontal wind velocity in the x,y directions. An oceanographic convention is used where velocity is positive when it is blowing toward positive coordinate directions. Units depend on the specific type of wind input file.

.. _WVNX:
.. _WVNY:

WVNX, WVNY
   (If NWS = 4, -4, 104, -104) applied horizontal wind velocity in the x,y directions. An oceanographic convention is used where velocity is positive when it is blowing toward positive coordinate directions. Units are knots.

.. _WVXFN:
.. _WVYFN:

WVXFN, WVYFN
   (If NWS = 6, 106, 7, -7) applied horizontal wind velocity in the x,y directions. An oceanographic convention is used where velocity is positive when it is blowing toward positive coordinate directions. Units are assumed to be m/s

.. _LONB:
.. _LATB:

LONB, LATB
   (If NWS = 10) number of longitude and latitudes in a global Gaussian Lon/Lat grid (NWS = 10) or ETA-29 grid (NWS = 11), these are specified in the ADCIRC program.

.. _PG:

PG
   (If NWS = 10) surface pressure in m H20.

.. _UG:

UG
   (If NWS = 10) 10 meter U velocity in m/s.

.. _VG:

VG
   (If NWS = 10) 10 meter V velocity in m/s.

.. _PE:

PE
   (If NWS = 11) surface pressure in mbars.

.. _UE:

UE
   (If NWS = 11) 10 meter U velocity in m/s.

.. _VE:

VE
   (If NWS = 11) 10 meter V velocity in m/s

.. _NWSET:

NWSET
   (If NWS = 12) the number of wind/pressure fields to be used; 1 indicates that only the basin scale meteorological field should be used; and 2 indicates that two fields be used (both the basin scale meteorological field and the regional scale meteorological field).

.. _NWBS:

NWBS
   (If NWS = 12) is useful in cases where the start time of the met data does not coincide with the start time of an ADCIRC run. NWBS represents the time delay between the start of an ADCIRC run and the start of the gridded wind field data, in units of wind time increments. Examples are provided in the notes section for this met file format.

.. _DWM:

DWM
   (If NWS = 12) specifies a multiplication factor for the wind velocities. If the velocities should be used as-is, set this factor to 1.0.

.. _comment_line:

comment line
   Provenance information, max 1024 characters

.. _pressureWindRelationship:

pressureWindRelationship
   This character string specifies the method for assigning the barometric pressure, since it is not provided in HWind data. Allowable values are dvorak, knaffzehr, specifiedPc, and background.

.. _hours:

hours(i)
   Specifies the model time corresponding to the data in file i, in hours since cold start if NWS=15, or in hours since hot start if NWS=-15

.. _centralPressure:

centralPressure(i)
   If the pressure-wind relationship was set to specifiedPc, this column specifies the central pressure corresponding to the wind velocity field in file i. Otherwise it is ignored and can be set to -1 to emphasize this.

.. _rampMult:

rampMult(i)
   Multiplier that will be applied to the wind velocities in file i; can be used to provide a meteorological ramp if the hwind data are to be applied long after coldstart when the usual ADCIRC ramp functions are unavailable. Normally set to 1.0.

.. _filename:

filename(i)
   The filename of the HWind file i, as downloaded from NOAA HRD. The filename should be enclosed in quotes; this prevents the Fortran library from interpreting embedded "/" characters in full path filenames as record separators.

.. _dvorak:

dvorak
   The central pressure will be assigned by searching the wind field for the maximum wind speed and then plugging it into the formula pc=1015-(Vmax/3.92)**(1/0.644).

.. _knaffzehr:

knaffzehr
   Same as dvorak except the formula is pc=1010-(Vmax/3.92)**(1/0.76).

.. _specifiedPc:

specifiedPc
   The central pressure will be linearly interpolated in time from the values provided in the centralPressure(i) column. This allows the use of the BEST track central pressures with the HWind velocities.

.. _background:

background
   This value will cause the entire barometric pressure field to the background pressure.

fort.23
~~~~~~~

.. _RSX:
.. _RSY:

RSX, RSY
   Applied wave radiation stress in the x,y directions divided by the reference density of water (should be units (length/time)2). An oceanographic convention is used where stress is positive when it is pointed in positive coordinate directions.

fort.24
~~~~~~~

.. _SALTAMP:

SALTAMP(k,JN)
   Amplitude of the self attraction/earth tide loading forcing for constituent k and node number JN. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m.

.. _SALTPHA:

SALTPHA(k,JN)
   Phase (degrees) of the self attraction/earth tide loading forcing for constituent k and node number JN.

Output Files
------------

fort.51,52,53,54
~~~~~~~~~~~~~~~~

.. _EMAG:

EMAG(j,k)
   Elevation amplitude for constituent k at station j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _PHASEDE:

PHASEDE(j,k)
   Elevation phase (in deg) for constituent k at node or station j. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _EMAGT:

EMAGT(j,k)
   Elevation amplitude for constituent k at node j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _UMAG:

UMAG(j,k)
   X direction velocity amplitude for constituent k at station j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _PHASEDU:

PHASEDU(j,k)
   X – direction velocity phase (in deg) for constituent k at node or station j. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _VMAG:

VMAG(j,k)
   Y direction velocity amplitude for constituent k at station j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _PHASEDV:

PHASEDV(j,k)
   Y – direction velocity phase (in deg) for constituent k at node or station j. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _UMAGT:

UMAGT(j,k)
   X direction velocity amplitude for constituent k at node or station j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _VMAGT:

VMAGT(j,k)
   Y direction velocity amplitude (in units of distance consistent with gravity) for constituent k at node or station j. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, amplitude is in m, if gravity is in ft/s2, amplitude is in ft.). If ICS=2, gravity must be specified in m/s2 and the amplitude is in m. This quantity is computed by the harmonic analysis routines in ADCIRC.

.. _fort.55:

fort.55
~~~~~~~

.. _EAV:

EAV(J)
   Mean elevation in the resynthesized time series

.. _EAVDIF:

EAVDIF(J)
   Elevation variance in the resynthesized time series

.. _ESQ:

ESQ(J)
   Elevation variance in the resynthesized time series

.. _EVADIF:

EVADIF(J)
   Elevation variance in the resynthesized time series

.. _UAV:

UAV(J)
   Mean x-velocity in the resynthesized time series

.. _UAVDIF:

UAVDIF(J)

.. _USQ:

USQ(J)
   X-velocity variance in the resynthesized time series

.. _UVADIF:

UVADIF(J)

.. _VAV:

VAV(J)
   Mean y-velocity in the resynthesized time series

.. _VAVDIF:

VAVDIF(J)

.. _VSQ:

VSQ(J)
   y-velocity variance in the resynthesized time series

.. _VVADIF:

VVADIF(J)

fort.61,62,63,64,71,72,73,74,75,81,83
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _TIME:

TIME
   Model time (in seconds) (TIME = STATIM*86400 + IT*DT)

.. _IRTYPE:

IRTYPE
   The record type (= 1 for elevation files, = 2 for velocity files, and = 3 for 3D velocity file)

.. _NDSETSE:

NDSETSE
   The number of data sets to be written to fort.63

.. _NDSETSV:

NDSETSV
   The number of data sets to be written to fort.64

.. _NDSETSW:

NDSETSW
   The number of data sets to be spooled to fort.73 & 74

.. _NDSETSC:

NDSETSC
   The number of data sets to be spooled to fort.83

.. _NTRSPE:

NTRSPE
   The number of data sets to be written to fort.61

.. _NTRSPV:

NTRSPV
   The number of data sets to be written to fort.62

.. _NTRSPM:

NTRSPM
   The number of data sets to be spooled to fort.71 or fort.72

.. _NTRSPC:

NTRSPC
   The number of data sets to be spooled to fort.81

.. _ET00:

ET00(k)
   Surface elevation at NSTAE elevation recording stations. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, elevation is in m, if gravity is in ft/s2, elevation is in ft.). If ICS=2, gravity must be specified in m/s2 and the elevation is in m.

.. _ETA2:

ETA2(k)
   Surface elevation at node k at the current time step. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, elevation is in m, if gravity is in ft/s2, elevation is in ft.). If ICS=2, gravity must be specified in m/s2 and the elevation is in m.

.. _ETAMAX:

ETAMAX(k)
   Maximum surface elevation over the entire run at node k. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, elevation is in m, if gravity is in ft/s2, elevation is in ft.). If ICS=2, gravity must be specified in m/s2 and the elevation is in m.

.. _UU00:

UU00(k)
   x, y velocity at NSTAV velocity recording stations. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _UU:
.. _VV:

UU(k), VV(k)
   Depth-averaged velocity in the x,y -coordinate direction at node k at the current time step. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _RMP00:

RMP00(k)
   Atmospheric surface pressure (m of water) output at NSTAM meteorological recording stations.

.. _RMU00:
.. _RMV00:

RMU00, RMV00(k)
   x,y wind stress (NWS =1, 2, -2) or velocity (NWS = 3, 4, -4, 5, -5, 6, 10, 11) at NSTAM meteorological recording stations. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _CC00:

CC00(k)
   Scalar concentration at NSTAC concentration recording stations.

.. _PR2:

PR2(k)
   Atmospheric surface pressure (m of water) output for all nodes in the domain.

.. _WVNXOUT:
.. _WVNYOUT:

WVNXOUT(k), WVNYOUT(k)
   x,y wind stress (NWS =1, 2, -2) or velocity (NWS = 3, 4, -4, 5, -5, 6, 10, 11) or converted from 1-min average wind to 10-min average wind velocity (one2ten=0.8928) multiplied by the wind ramp in m/s (NWS = 8, 19, 20) for all nodes in the domain. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _C1:

C1(k)
   Scalar concentration output for all nodes in domain

fort.67, 68
~~~~~~~~~~~

.. _IT:

IT
   Model time step number since the beginning of the model run.

.. _ETA1:

ETA1(k)
   Surface elevation at node k at the previous time step

.. _ICSTP:

ICSTP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Scalar Concentration Time Series at Specified Concentration Recording Stations output file.

.. _IESTP:

IESTP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Elevation Time Series at Specified Elevation Recording Stations output file.

.. _IPSTP:

IPSTP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Atmospheric Pressure Time Series at Specified Meteorological Recording Stations output file.

.. _IVSTP:

IVSTP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Depth-averaged Velocity Time Series at Specified Velocity Recording Stations output file.

.. _IWSTP:

IWSTP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Wind Velocity Time Series at Specified Meteorological Recording Stations output file.

.. _IGCP:

IGCP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Scalar Concentration Time Series at All Nodes in the Model Grid output file.

.. _IGEP:

IGEP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Elevation Time Series at All Nodes in the Model Grid output file.

.. _IGPP:

IGPP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Atmospheric Pressure Time Series at All Nodes in the Model Grid output file.

.. _IGVP:

IGVP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Depth-averaged Velocity Time Series at All Nodes in the Model Grid output file.

.. _IGWP:

IGWP
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the Wind Stress or Velocity Time Series at All Nodes in the Model Grid output file.

.. _noff:

noff(i)
   The wet/dry state of element i where 1 indicates an element is categorized as wet on the timestep that the dataset was written while a value of 0 indicates that an element is categorized as dry.

.. _EP:

EP
   Scaling parameter that is used to maximize the diagonal dominance of the GWCE system matrix. It is computed within ADCIRC at a cold start or when the GWCE matrix changes (e.g., wetting and drying has occurred).

.. _NSCOUC:

NSCOUC
   Time step counter to determine when the next entry will be written to the Scalar Concentration Time Series at Specified Concentration Recording Stations output file.

.. _NSCOUE:

NSCOUE
   Time step counter to determine when the next entry will be written to the Elevation Time Series at Specified Elevation Recording Stations output file.

.. _NSCOUM:

NSCOUM
   Time step counter to determine when the next entry will be written to the Atmospheric Pressure Time Series at Specified
   Meteorological Recording Stations and Wind Velocity Time Series at Specified Meteorological Recording Stations output files.

.. _NSCOUV:

NSCOUV
   Time step counter to determine when the next entry will be written to the Depth-averaged Velocity Time Series at Specified Velocity Recording Stations output file.

.. _NSCOUGC:

NSCOUGC
   Time step counter to determine when the next entry will be written to the Scalar Concentration Time Series at All Nodes in the Model Grid output file.

.. _NSCOUGE:

NSCOUGE
   Time step counter to determine when the next entry will be written to the Elevation Time Series at All Nodes in the Model Grid output file.

.. _NSCOUGW:

NSCOUGW
   Time step counter to determine when the next entry will be written to the Atmospheric Pressure Time Series at All Nodes in the Model Grid and Wind Stress or Velocity Time Series at All Nodes in the Model Grid output files.

.. _NSCOUGV:

NSCOUGV
   Time step counter to determine when the next entry will be written to the Depth-averaged Velocity Time Series at All Nodes in the Model Grid output file.

.. _N3DSD:

N3DSD
   Time step counter to determine when the next entry will be written to the 3D Density, Temperature and/or Salinity at Specified Recording Stations (fort.41) output file.

.. _I3DSDRec:

I3DSDRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Density, Temperature and/or Salinity at Specified Recording Stations (fort.41) output file.

.. _N3DSV:

N3DSV
   Time step counter to determine when the next entry will be written to the 3D Velocity at Specified Recording Stations (fort.42) output file.

.. _I3DSVRec:

I3DSVRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Velocity at Specified Recording Stations (fort.42) output file.

.. _N3DST:

N3DST
   Time step counter to determine when the next entry will be written to the 3D Turbulence at Specified Recording Stations (fort.43) output file.

.. _I3DSTRec:

I3DSTRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Turbulence at Specified Recording Stations (fort.43) output file.

.. _N3DGD:

N3DGD
   Time step counter to determine when the next entry will be written to the 3D Density, Temperature and/or Salinity at All Nodes in the Model Grid (fort.44) output file.

.. _I3DGDRec:

I3DGDRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Density, Temperature and/or Salinity at All Nodes in the Model Grid (fort.44) output file.

.. _N3DGV:

N3DGV
   Time step counter to determine when the next entry will be written to the 3D Velocity at All Nodes in the Model Grid (fort.45) output file.

.. _I3DGVRec:

I3DGVRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Velocity at All Nodes in the Model Grid (fort.45) output file.

.. _N3DGT:

N3DGT
   Time step counter to determine when the next entry will be written to the 3D Turbulence at All Nodes in the Model Grid (fort.46) output file.

.. _I3DGTRec:

I3DGTRec
   Line number (for ASCII output) or record number (for binary output) of the most recent entry in the 3D Turbulence at All Nodes in the Model Grid (fort.46) output file.

.. _DUU:
.. _DUV:
.. _DVV:

DUU(k), DUV(k), DVV(k)
   Dispersion terms

.. _UU, VV:

UU(k), VV(k)
   Depth-averaged horizontal velocity

.. _BSX:
.. _BSY:

BSX(k), BSY(k)
   x, y bottom stresses

.. _VIDBCPDX:
.. _VIDBCPDY:

VIDBCPDX(k), VIDBCPDY(k)
   Vertically integrated Baroclinic pressure

.. _REAL(Q(k,j)):

REAL(Q(k,j))
   u velocity component

.. _AIMAG(Q(k,j)):

AIMAG(Q(k,j))
   v velocity component

.. _WZ(k,j):

WZ(k,j)
   Vertical velocity

.. _Q20(k,j):

Q20(k,j)
   Previous time step value of Q2

.. _HA:

HA(k,j)
   Coefficients in the least squares matrix used for the harmonic analysis

.. _CH1:

CH1(k)
   Depth-averaged scalar concentration value at node k at the current time step

.. _ELAV:

ELAV(k)
   Sum of elevations computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

.. _ELVA:

ELVA(k)
   Sum of squares of elevations computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

.. _FNAM8(1):

FNAM8(1)
   The first 8 characters of FNAME(k)

.. _FNAM8(2):

FNAM8(2)
   The second 8 characters of FNAME(k)

.. _GLOELV:

GLOELV(k,j)
   Harmonic analysis load vectors for elevation at all nodes in the model grid

.. _GLOULV:

GLOULV(k,j)
   Harmonic analysis load vectors for depth-averaged u velocity at all nodes in the model grid

.. _GLOVLV:

GLOVLV(k,j)
   Harmonic analysis load vectors for depth-averaged v velocity at all nodes in the model grid

.. _ICALL:

ICALL
   Number of times the harmonic analysis has been updated

.. _ICHA:

ICHA
   Time step counter to determine when the next update will be made to the harmonic analysis.

.. _IHARIND:

IHARIND
   Indicator of whether any harmonic analysis will be performed during the model run.

.. _ITUD:

ITUD
   Model time step when the harmonic analysis was last updated

.. _MM:

MM
   2*NFREQ – NF

.. _NF:

NF
   Indicator of whether the steady frequency is included in the harmonic analysis (NF = 1, steady is included; NF = 0, steady is not included).

.. _NTSTEPS:

NTSTEPS
   Number of time steps since harmonic analysis means and variance checking has begun

.. _NZ:

NZ
   Indicator of whether the steady frequency is included in the harmonic analysis (NZ = 0, steady is included; NZ = 1, steady is not included).

.. _STAELV:

STAELV(j,k)
   Harmonic analysis load vectors for elevation at elevation recording stations

.. _STAULV:

STAULV(j,k)
   Harmonic analysis load vectors for depth-averaged u velocity at velocity recording stations

.. _STAVLV:

STAVLV(j,k)
   Harmonic analysis load vectors for depth-averaged v velocity at velocity recording stations

.. _TIMEUD:

TIMEUD
   Model time when the harmonic analysis was last updated

.. _XVELAV:

XVELAV(k)
   Sum of depth-averaged u velocities computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

.. _XVELVA:

XVELVA(k)
   Sum of squares of depth-averaged u velocities computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

.. _YVELAV:

YVELAV(k)
   Sum of depth-averaged v velocities computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

.. _YVELVA:

YVELVA(k)
   Sum of squares of depth-averaged v velocities computed by ADCIRC, at every node k in the model grid, over all time steps since harmonic analysis means and variance checking has begun

fort.41, fort.42, fort.43, fort.44, fort.45, fort.46
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _NDSET3DSD:

NDSET3DSD
   Number of data sets to be written to fort.41

.. _SIGTSTA:

SIGTSTA
   Sigma T value (kg/m^3) (=density-1000) at a specified recording station

.. _SALSTA:

SALSTA
   Salinity (PSU) at a specified recording station

.. _TEMPSTA:

TEMPSTA
   Temperature (DEG C) at a specified recording station

.. _NDSET3DSV:

NDSET3DSV
   Number of data sets to be written to fort.42

.. _REAL(QSTA(k,j)):

REAL(QSTA(k,j))
   u velocity component for station output

.. _AIMAG(QSTA(k,j)):

AIMAG(QSTA(k,j))
   v velocity component for station output

.. _WZSTA:

WZSTA(k,j)
   Vertical velocity for station output. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _NDSET3DST:

NDSET3DST
   Number of data sets to be written to fort.43

.. _q20STA:

q20STA(M)
   Turbulent kinetic energy for station output

.. _ISTA:

ISTA(M)
   Mixing length for station output

.. _EVSTA:

EVSTA(M)
   Spatially constant horizontal eddy viscosity for the momentum equations (units of length2/time) for station output

.. _REAL(Q(j,M)):

REAL(Q(j,M))
   u velocity component for global output

.. _AIMAG(Q(j,M)):

AIMAG(Q(j,M))
   v velocity component for global output

.. _WZ(j,M):

WZ(j,M)
   Vertical velocity for global output. If ICS=1, units are determined by the units of gravity (G) specified in the fort.15 file (eg. if G is in m/s2, velocity is m/s, if gravity is in ft/s2, velocity is in ft/s). If ICS=2, gravity must be specified in m/s2 and the velocity units are m/s.

.. _q20(j,M):

q20(j,M)
   Turbulent kinetic energy

.. _l(j,M):

l(j,M)
   Mixing length

.. _EV(j,M):

EV(j,M)
   Horizontal eddy viscosity for the momentum equations (units of length2/time) for global output

.. _NDSET3DGD:

NDSET3DGD
   The number of data sets to be written to fort.44

.. _NDSET3DGV:

NDSET3DGV
   The number of data sets to be written to fort.45

.. _NDSET3DGT:

NDSET3DGT
   The number of data sets to be written to fort.46

.. _NTF:

NTF
   Temperature boundary condition file type (this file is not supported yet).

.. _NCPROJ:

NCPROJ
   Project Title (what is in the file).

.. _NCINST:

NCINST
   Project Institution (where file was produced).

.. _NCSOUR:

NCSOUR
   Project Source (how it was produced, eg. instrument type).

.. _NCHIST:

NCHIST
   Project History (audit trail of processing operations).

.. _NCREF:

NCREF
   Project References (pointers to publications, web documentation).

.. _NCCOM:

NCCOM
   Project Comments (any other comments about the file).

.. _NCHOST:

NCHOST
   Project Host.

.. _NCCONV:

NCCONV
   Conventions.

.. _NCCONT:

NCCONT
   Contact Information.

.. _NCDATE:

NCDATE
   The format of NCDATE must be as follows so that ADCIRC can create netcdf files that comply with the CF standard: yyyy-MM-dd hh:mm:ss tz. For example, if the cold start date/time of the run is midnight UTC on 1 May 2010, the NCDATE parameter should be set to 2010-05-01 00:00:00 UTC.

.. _WindDragLimit:

WindDragLimit
   This parameter controls the ceiling on the wind drag coefficient. The default value is 0.0035. Another commonly used value is 0.002.

.. _DragLawString:

DragLawString
   This character string controls the formulation used to calculate the nodal wind drag coefficients on the water surface from the nodal wind speed. The default value "Garratt" specifies the formula WindDragCoefficient = 0.001d0 * (0.75d0 + 0.067d0 * WindSpeed). Alternatively, for tropical cyclones, setting the 'DragLawString' equal to "Powell" causes a more complex formulation to be used, where different formulas apply to different storm quadrants. If ice coverage data have been specified via appropriate specification of the NWS parameter, the wind drag coefficient formulation can also be set to IceCube (the default) or RaysIce. Complete details are available in the documentation of the Ice Coverage Input Files (fort.25, fort.225, fort.227).

.. _rhoAir:

rhoAir
   Used to set the density of air; the default value is 1.15 kg/m^3.

.. _waveCoupling:

waveCoupling
   Fortran namelist that can be used to control the magnitude of the winds that are passed to coupled wave models.

.. _swanoutputcontrol:

swanoutputcontrol
   Fortran namelist that can be used to control the SWAN model related output control.

.. _metControl:

metControl
   Fortran namelist that can be used to control the meteorological forcing related parameters.

.. _WindWaveMultiplier:

WindWaveMultiplier
   This multiplier is applied after the winds have been read in by ADCIRC and time interpolated (if necessary) but without ADCIRC's ramp function. Since ADCIRC normally requires 10 minute averaged winds, many of ADCIRC's meteorological reading subroutines also include a conversion from the time averaging period of the source data (if such time averaging is a standard for that data type) to a 10 minute average. In these cases, the multiplier is applied to the derived 10 minute averaged wind velocity values, rather than the raw values from the original meteorological forcing file.

.. _timeBathyControl:

timeBathyControl
   Fortran namelist used to specify parameters related to simulations that have time varying bathymetry.

.. _NDDT:

NDDT
   Controls whether time varying bathymetry will be used, and if- so, the spatial extent and time reference of the data in the- time varying bathymetry input file (fort.141). See description below for more details.

   NDDT=0 Time varying bathymetry should not be used.

   NDDT>0 Time varying bathymetry should be used, and the beginning of the data in the time varying bathymetry file (fort.141) corresponds to the cold start time.

   NDDT<0 Time varying bathymetry should be used, and the beginning of the data in the time varying bathymetry file (fort.141) corresponds to the hotstart time.

   NDDT=1 or NDDT=-1 The time varying bathymetry file (fort.141) is a fulldomain file, that is, all nodes in the mesh are specified at each time increment in the file.

   NDDT=2 or NDDT=-2 The time varying bathymetry file (fort.141) contains nodal values over a limited area or subsection of the fulldomain mesh.

.. _BTIMINC:

BTIMINC
   Time increment (in seconds) between time varying bathymetry datasets in the time varying bathymetry (fort.141) file.

.. _BCHGTIMINC:

BCHGTIMINC
   Time increment (seconds) over which bathymetry changes during a BTIMINC interval.

.. _tau0var:

tau0var
   tau0 value written out at every node in the fort.90 file

.. _subdomainModeling:

subdomainModeling
   Fortran namelist used to activate subdomain modeling, which gives analysts the capability to define a small subdomain nested within a much larger domain. The fulldomain can thenbe run once, generating boundary conditions for the nested subdomain. These boundary conditions can then be used over and over for different internal configurations of the subdomain without having to run the larger domain again. Please see the ADCIRC documentation on subdomain modeling for further details.

.. _subdomainOn:

subdomainOn
   This logical variable activates subdomain modeling when set to '.true.'. The default value of '.false.' will be used if the subdomainModeling namelist is not found in the fort.15 file.

.. _wetDryControl:

wetDryControl
   Fortran namelist that can be used to turn on output of (1) the nodal wet/dry state to the nodecode.63 file and/or (2) the elemental wet/dry state to the noff.100 file. It can also be used to deactivate the use of elemental wet/dry state in the continuity and momentum equations by forcing all elements of the NOFF array to be 1.

.. _inundationOutputControl:

inundationOutputControl
   Fortran namelist that can be used to activate the production of a set of inundation-related output files, including initiallydry.63, everdried.63, inundationtime.63, maxinundepth.63, and endrisinginun.63. It is also used to set the inunThresh variable, i.e., the threshold of water depth above local ground level at which normally dry land is considered inundated (to eliminate trivial inundation levels from consideration).

.. _TVWControl:

TVWControl
   Fortran namelist that can be used to activate the time varying weirs feature (with use_TVW), specify the name of the input file that specifies the behavior of the time varying weir(s) (TVW_file), and the schedule for output of weir height to the fort.77 file (nout_TVW, touts_TVW, toutf_TVW, and nspool_TVW).

.. _outputNodeCode:

outputNodeCode
   Logical variable activates the production of the nodecode.63 file for nodal wet/dry state when set to '.true.'. The default value is '.false.'. The nodecode.63 output file is produced in ASCII format with the same structure as the fort.63 file (with the exception that the node code value is an integer) and on the same schedule as the fort.63 file. A nodal value of 1 in the nodecode.63 file indicates that the node is wet on the timestep that the data set was written, while a value of 0 indicates that the node is dry on the timestep that the data set was written. In order to be included in the computations of water surface elevation and velocity, a node must be wet.

.. _outputNOFF:

outputNOFF
   Logical variable activates the production of the noff.100 file for elemental wet/dry state when set to '.true.'. The default value is '.false.'. The noff.100 output file is produced in ASCII format with a structure similar to the fort.63 file, with the following exceptions: (a) the noff value is an integer; and (b) the noff array is elementally based, so the number of values for each time snap is equal to the number of elements, rather than the number of nodes. When outputNOFF is set to .true., the noff.100 file will be produced on the same schedule as the fort.63 file. An elemental value of 1 in the noff.100 file indicates that the element is wet on the timestep that the data set was written, while a value of 0 indicates that the element is dry on the timestep that the data set was written. In order to be included in the computations of water surface elevation and velocity, an element and all three of its nodes must be wet.

.. _noffActive:

noffActive
   NOFF is the name of elemental array in ADCIRC that holds the elemental wet/dry state. The noffActive logical variable can be used to completely disable the elemental wet dry state by setting NOFF to 1 (wet), always and everywhere, when set to '.false.'. The default value is '.true.'.

.. _inundationOutput:

inundationOutput
   Logical variable activates the production of the the inundation output files (consisting of initiallydry.63, everdried.63, inundationtime.63, maxinundepth.63, and endrisinginun.63) when set to '.true.'. The default value is '.false.'.

.. _inunThresh:

inunThresh
   Numerical value of water surface elevation above local ground level that marks the threshold at which normally dry land is considered inundated. The default value is 0.6m.

.. _use_TVW:

use_TVW
   Logical variable that activates the time varying weirs capability when set to '.true.'. The default value is '.false.'.

.. _TVW_file:

TVW_file
   Character string used to specify the name of the time varying weirs input file. The default name is 'fort.tvw'. This file defines the location(s) and schedule(s) or trigger(s) for weir elevation changes. If the TVW_file is not found, ADCIRC writes an informational message to the screen and to the fort.16 (log) file that the file was not found. The run will then continue, but all weirs will be static.

.. _nout_TVW:

nout_TVW
   Integer that activates the production of the fort.77 output file as well as specifying its format; setting nout_TVW to zero turns off the production of the fort.77 output file while setting it to positive or negative 1, 3, 4, or 5 produce files in ASCII, netCDF3, sparse ASCII, or netCDF4, respectively. Due to the sparse nature of the weir height output data, sparse ASCII format (nout_TVW=4) is recommended. If use_TVW is set to .true. but this parameter nout_TVW is not set at all, then nout_TVW is set to the same value as the file format specifier for the fort.63 file (NOUTGE). The data stored in the fort.77 file represent the change in the weir height (in meters, positive upward) compared to the weir height specified in the ADCIRC mesh file, rather than the absolute height of the weir at a particular time.

.. _touts_TVW:

touts_TVW
   Time in days since cold start after which output to the fort.77 file will start.

.. _toutf_TVW:

toutf_TVW
   Time in days since cold start when output to the fort.77 file will end.

.. _nspool_TVW:

nspool_TVW
   Time step increment at which output will be written to the fort.77 file when the ADCIRC model time is between touts_TVW and toutf_TVW (in days since cold start).

.. _TVW:

TVW
   Nodal array that contains the difference between the current weir elevations and their canonical elevations as specified in the ADCIRC mesh file (fort.14).

.. _NSTAE2:

NSTAE2
   Number of stations specified in the Elevation Station location file (elev_sta.151)

.. _NSTAV2:

NSTAV2
   Number of stations specified in the Velocity Station location file (vel_sta.151)

.. _NSTAC2:

NSTAC2
   Number of stations specified in the Concentration Station location file (conc_sta.151)

.. _NSTAM2:

NSTAM2
   Number of stations specified in the Meteorological Station location file (met_sta.151)

.. _initiallydry:

initiallydry(k)
   A value given to node k at cold start based on whether the node was dry or wet. If it is dry the value is 1, if it is wet, the value is 0.

.. _everdried:

everdried(k)
   A value at node k of -99999.0 if node k is dried at any time during the simulations, and 1 if it remains wet.

.. _driedtime:

driedtime(k)
   Total time in seconds that node k was dry during the simulation (0 if it remains wet)

.. _endrisinginun:

endrisinginun(k)
   If node k had rising inundation levels at the end of a simulation it is flagged with an integer value of 1, otherwise it is given an integer value of 0

.. _maxinundepth:

maxinundepth(k)
   The peak inundation depth at node k (in meters) above ground that occurred during the simulation

.. _maxinundepth_time:

maxinundepth_time(k)
   The time of the the peak inundation depth at the node in seconds since the cold start

.. _inundationtime:

inundationtime(k)
   The total accumulated time in seconds that node k was inundated beyond the threshold. Periods of inundation are counted toward the total time, even if they are not contiguous.

.. _inundationtime_onset:

inundationtime_onset(k)
   The time of onset of inundation beyond the threshold at node k in seconds since cold start. The time of onset data are useful in the context of real time model guidance for decision making.

.. _nodecode:

nodecode(k)
   The wet/dry state of node k where 1 indicates a node is categorized as wet on the timestep that the dataset was written while a value of 0 indicates that a node is categorized as dry.
