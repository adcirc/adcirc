Fort.13: Nodal Attributes File
==============================

The fort.13 file contains nodal attributes that are constant in time but spatially variable. This file is only read when `NWP > 0` in the `Model Parameter and Periodic Boundary Condition File`.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input.

.. parsed-literal::

   :ref:`AGRID <AGRID>`
   :ref:`NumOfNodes <NumOfNodes>`
   :ref:`NAttr <NAttr>`
   for i=1 to :ref:`NAttr <NAttr>`
      :ref:`AttrName(i) <AttrName>`
      :ref:`Units(i) <Units>`
      :ref:`ValuesPerNode(i) <ValuesPerNode>`
      for k=1 to :ref:`ValuesPerNode(i) <ValuesPerNode>`
         :ref:`DefaultAttrVal(i,k) <DefaultAttrVal>`
      end k loop
   end i loop
   for i=1 to :ref:`NAttr <NAttr>`
      :ref:`AttrName(i) <AttrName>`
      :ref:`NumNodesNotDefaultVal(i) <NumNodesNotDefaultVal>`
      for j=1 to :ref:`NumNodesNotDefaultVal(i) <NumNodesNotDefaultVal>`
         :ref:`n <n>` :ref:`AttrVal(n,k) <AttrVal>` for k=1 to :ref:`ValuesPerNode(i) <ValuesPerNode>`
         end k loop
      end j loop
   end i loop

Nodal Attributes
----------------

The following nodal attributes are available in ADCIRC:

* :ref:`primitive_weighting_in_continuity_equation`
* :ref:`surface_submergence_state`
* :ref:`quadratic_friction_coefficient_at_sea_floor`
* :ref:`surface_directional_effective_roughness_length`
* :ref:`surface_canopy_coefficient`
* :ref:`bottom_roughness_length`
* :ref:`average_horizontal_eddy_viscosity_in_sea_water_wrt_depth`
* :ref:`elemental_slope_limiter`
* :ref:`advection_state`
* :ref:`initial_river_elevation`
* :ref:`condensed_nodes`


The following sections describe each available nodal attribute and its properties.


.. _primitive_weighting_in_continuity_equation:

primitive_weighting_in_continuity_equation – Tau0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Same as existing TAU0.
   
   **Requirement**: Required, but not necessarily in the fort.13 file as there is a choice of specification methods. Can be specified in the Model Parameter and Periodic Boundary Condition File as a positive constant, in which case it is spatially uniform; or a negative constant, in which case it is spatially varying according to a hardcoded scheme based on depth. Can also be specified as a nodal attribute in the Nodal Attributes File (fort.13), in which case any value specified in the Model Parameter and Periodic Boundary Condition File is ignored (nodal attributes take precedence).
   
   **Units**: Unitless. (Units for nodal attributes are specified by the user).
   
   **Number of values per node**: 1.

   **Values**: Suggested range specified in description of TAU0.

.. _surface_submergence_state:

surface_submergence_state – StartDry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: This attribute can force a node to be categorized as dry at the start of the simulation, even if it is below the geoid. This is important for simulating places like New Orleans that are below sea level but are not underwater.
   
   **Requirement**: Optional.
   
   **Units**: Unitless.
   
   **Number of values per node**: 1.
   
   **Values**: If set to 1, the node is categorized as dry at the cold start of the simulation. If set to zero, the node is categorized as wet or dry depending on whether its depth is below or above the geoid.

.. _quadratic_friction_coefficient_at_sea_floor:

quadratic_friction_coefficient_at_sea_floor – Fric
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: See the description of CF in the Model Parameter and Periodic Boundary Condition File. If the user elects to load this from the fort.13 file, NOLIBF must be set to 1 or the run will terminate.
   
   **Requirement**: Optional.
   
   **Units**: Unitless.
   
   **Number of values per node**: 1.
   
   **Values**: Same as CF.

.. _surface_directional_effective_roughness_length:

surface_directional_effective_roughness_length – z0Land
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: A measure of the "roughness" of the land that can impede wind flow and reduce the surface stress that the wind applies. The ocean would be considered very smooth, and skyscrapers would be considered very rough. This attribute is directional, and the twelve values represent the roughness lengths "seen" by winds blowing from twelve different compass directions at each node. The orientation of the twelve values follows the trigonometric convention, that is, zero degrees represents due east, and the values proceed counter clockwise. In other words, the first value at a node is applied to winds blowing from west to east, the second value applies to winds blowing East-Northeast, etc.
   
   **Requirement**: Optional.
   
   **Units**: Specified by the user, as is the case for all nodal attributes. The data we use is provided in meters.
   
   **Number of values per node**: 12.
   
   **Values**: Greater than or equal to zero.

.. _surface_canopy_coefficient:

surface_canopy_coefficient – VCanopy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: This attribute allows the user to turn off wind stress in heavily forested areas that have been flooded, like a swamp. The canopy shields the water from the effect of the wind.
   
   **Requirement**: Optional.
   
   **Units**: Unitless.
   
   **Number of values per node**: 1.
   
   **Values**: Zero if the wind stress should be zero because of a canopy. One otherwise.

.. _bridge_pilings_friction_paramenters:

bridge_pilings_friction_paramenters – BK, BAlpha, BDelX, POAN
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Extra friction from bridge pilings. See Note below.

   **Requirement**: Very optional.
   
   **Units**: Vary, see note below and BK, BALPHA, BDELX and POAN.

   **Number of values per node**: 4.
   
   **Values**: Vary, see note below and BK, BALPHA, BDELX and POAN.

   **Note**: Special considerations must be used when designing a grid for an ADCIRC application that includes the effects of bridge pilings. Specifically, it is necessary to build the grid to provide at least three rows of nodes that parallel the bridge span. One row of nodes (centerline nodes) should lie along the approximate centerline of the bridge while the second and third rows of nodes (adjacent nodes) should lie on either side of the centerline nodes in the along steam direction. An initial implementation of obstruction drag in ADCIRC placed this drag entirely at the row of centerline nodes. However, tests showed that this arrangement led to significant oscillations in the numerical solution. The oscillations abated when the obstruction drag was distributed in the along stream direction so that 25 percent was located at each row of adjacent nodes and 50 percent was located at the row of centerline nodes. Node numbers and coefficient values at all nodes on the centerline and two adjacent rows must be entered in this input file. It is not necessary for centerline nodes to correspond to actual piling positions, (i.e., in the cross stream direction), since the overall effect of the pilings on the large scale circulation is all that is being represented. It is important, however, to construct a grid that is as uniform as possible in the vicinity of the bridge.

.. _mannings_n_at_sea_floor:

mannings_n_at_sea_floor – ManningsN
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Manning's n. If the user elects to use this nodal attribute, NOLIBF must be set to 1 or the run will terminate. During execution, the Manning's n value specified here is converted to an equivalent quadratic friction coefficient before the bottom stress is calculated. The equivalent quadratic friction coefficient is calculated according to the following formula at each node at each time step: Cd(t)=(g*n^2)/cuberoot(depth[+eta(t)]) where depth is the bathymetric depth and [eta(t)] is the water surface elevation. The addition of the water surface elevation is conditional upon the setting of NOLIFA: eta(t) is treated as zero if NOLIFA is set to zero in the fort.15 file. Finally, the value of CF in the fort.15 is used to set a lower limit on the resulting equivalent quadratic friction coefficient, since the Cd calculated from this formula tends to become small in deep water.

   **Requirement**: Optional.
   
   **Units**: Specified by user.
   
   **Number of values per node**: 1.
   
   **Values**: Greater than zero.

.. _chezy_friction_coefficient_at_sea_floor:

chezy_friction_coefficient_at_sea_floor – ChezyFric
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Chezy friction coefficient. If the user elects to use this nodal attribute, NOLIBF must be set to 1 or the run will terminate. 

   **Requirement**: Optional.
   
   **Units**: Specified by user.
   
   **Number of values per node**: 1.
   
   **Values**: Greater than zero.

.. _sea_surface_height_above_geoid:

sea_surface_height_above_geoid – GeoidOffset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Creates an initial offset of the sea surface from the geoid. If the offset is also specified at the boundaries, it will remain throughout the simulation. This has been used to simulate a steric effect, where water levels are higher in warm seasons because of thermal expansion.

   **Requirement**: Optional.
   
   **Units**: Specified by the user (length).
   
   **Number of values per node**: 1.
   
   **Values**: Any.

.. _wave_refraction_in_swan:

wave_refraction_in_swan – SwanWaveRefrac
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: To turn wave refraction on or off in the SWAN calculations during a coupled ADCIRC+SWAN run. In an ADCIRC-only run, this nodal attribute is ignored.

   **Requirement**: Optional.
   
   **Units**: Unitless.

   **Number of values per node**: 1.
   
   **Values**: If set to 1 at a node, wave refraction will be active at that node in the SWAN calculations during a coupled ADCIRC+SWAN run. If set to 0, wave refraction will be deactivated at that node in a coupled ADCIRC+SWAN run.

.. _bottom_roughness_length:

bottom_roughness_length – Z0b_var
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: Roughness length for 3D bottom friction calculations. Has no effect on a 2DDI ADCIRC run.

   **Requirement**: Optional.
   
   **Units**: Length (m)

   **Number of values per node**: 1.

   **Values**: Greater than zero. A reasonable range for the bottom_roughness_length would be 0.001m -0.2m , (densely vegetated overland values would have an upper limit of about 0.2-0.5m and smooth muddy bottoms could have values as low as 0.0001 m).

.. _average_horizontal_eddy_viscosity_in_sea_water_wrt_depth:

average_horizontal_eddy_viscosity_in_sea_water_wrt_depth – EVC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: See description of ESLC in the Model Parameter and Periodic Boundary Condition File.

   **Requirement**: Optional.
   
   **Units**: Specified by the user. ((length**2)/time).

   **Number of values per node**: 1.
   
   **Values**: Greater than or equal to zero.

.. _elemental_slope_limiter:

elemental_slope_limiter
^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: This nodal attribute is used to selectively limit the maximum elevation gradient that can occur across an element, thus improving numerical stability. Alternatively, it is also capable of merely logging individual elements where a specified elemental slope is exceeded at some point during the simulation. When this nodal attribute is loaded, warning messages will be written to the screen and to the fort.16 log file whenever the elevation gradient meets or exceeds the threshold value for the first time at a particular node. Once the elevation gradient is met or exceeded at a node, the elemental slope limiter remains active at that node for the remainder of the ADCIRC run. When the elemental slope limiter is active at a node, the water surface elevation at that node is reset to the average of the water surface elevations of the surrounding nodes. At the end of the run, a file called ESLNodes.63 will be written to indicate the nodes where the threshold elevation gradient was met or exceeded. Furthermore, If there is an ESLNodes.63 file in the input directory when ADCIRC starts, ADCIRC will load the existing ESLNodes.63 file along with the hotstart file so that the simulation can pick up where it left off, in terms of the elements where the slope is actively limited. When such a simulation finishes, it will overwrite the existing ESLNodes.63 file with a new one that reflects the updated state of limited elemental slopes.

   **Requirement**: Optional.

   **Units**: length/length or unitless.

   **Number of values per node**: 1.

   **Values**: Zero indicates that slope limiting is always active at that node, because a zero elevation gradient will always be met or exceeded; a positive value indicates the maximum gradient to be allowed at that node, at or beyond which the slope limiter is activated; a negative number indicates that the elevation gradients should be compared to the absolute value of the nodal attribute at this node. In the case of a negative number, ADCIRC will log a warning to the screen and to the fort.16 file the first time the elemental slope is exceeded, but ESLNodes.63 file will not be affected, and the elemental slope limiting will not actually occur at the node. A suggested value for this nodal attribute is 0.001.

   **ADCIRC Variable**: elemental_slope_limiter_grad_max

.. _advection_state:

advection_state – AdvectionState
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: The advection_state nodal attribute is used to provide fine grained control over the NOLICA and NOLICAT parameters, so that they can be set on an element-by-element basis. The bathymetric depth at the three nodes of an element are compared to the corresponding nodal attribute values, and if the bathymetric depth at all three nodes is greater than or equal to the corresponding nodal attribute value, the values of NOLICA and NOLICAT will be set to the value indicated in the Model Parameter and Periodic Boundary Condition (fort.15) file. If the bathymetric depth at any of the three nodes of a particular element is less than the corresponding nodal attribute value, the values of NOLICA and NOLICAT will be set to zero on that element.

   **Requirement**: Optional.

   **Units**: Length.

   **Number of values per node**: 1.

   **Values**: Any.

   **ADCIRC Variable**: AdvectionState

.. _initial_river_elevation:

initial_river_elevation
^^^^^^^^^^^^^^^^^^^^^^^

   **Description**: The initial_river_elevation nodal attribute is used to set the initial water surface elevation in inland rivers that have some portion of the river bed above mean sea level as well as an upstream flux boundary condition. ADCIRC assumes by default that vertices that are above mean sea level (i.e., with negative depths) will be dry when the simulation starts. This assumption is violated when the elevation of a river bed is above mean sea level, and there is an upstream flux boundary. This nodal attribute is used in those cases to provide the initial water surface elevation of the river at cold start.

   **Requirement**: Required only if the domain contains a river with an upstream flux boundary condition and a bed elevation above mean sea level.

   **Units**: Length.

   **Values**: Any.

   **ADCIRC Variable**: Eta2.

.. _condensed_nodes:

condensed_nodes
^^^^^^^^^^^^^^^

   **Description**: The condensed_nodes nodal attribute is used to specify a group of nodes, the nodal equations of which are to be condensed, and thus the solutions of which become identical. This condensation technique is useful when two or more nodes are within a close proximity and thus the model violates the CFL condition unless a smaller time step is used. The condensation is conducted at the level of nodal equations, which effectively extends the stencil of the nodes to relax the CFL condition. One typical use of the condensation is to relax the CFL condition along very narrow channels. The node pairs on the sides of a channel are specified as the condensed_nodes nodal attributes so that the model becomes insensitive to the width of the channel.

   **Requirement**: Optional.

   **Units**: Unitless.

   **Number of values per node**: Maximum number of nodes in a condensed node group. See :ref:`Example Fort.13 with Condensed Nodes <example_fort13_condensed_nodes>` for more information.

   **Values**: Node numbers.

   **ADCIRC Variable**: CondensedNodes.

Example
-------

The following is a simple example of a fort.13 file with two attributes:

.. code-block:: none

   Example ADCIRC Grid
   800    ! Total number of nodes in the grid
   2      ! Number of attributes
   primitive_weighting_in_continuity_equation
   unitless
   1      ! Number of values per node
   0.05   ! Default value for primitive_weighting_in_continuity_equation
   surface_submergence_state
   unitless
   1      ! Number of values per node
   0      ! Default value for surface_submergence_state
   primitive_weighting_in_continuity_equation
   2      ! Number of nodes with non-default values
   1 0.01 ! Node 1 has a value of 0.01
   3 0.02 ! Node 3 has a value of 0.02
   surface_submergence_state
   1      ! Number of nodes with non-default values
   4 1    ! Node 4 has a value of 1

This example demonstrates how to specify two attributes with different default values and how to override the default values for specific nodes.
