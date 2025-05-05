Wetting and Drying
==================

The ability to simulate the dynamic wetting and drying of elements is crucial for accurately modeling coastal flooding, tidal flats, and storm surge inundation. ADCIRC implements sophisticated algorithms to handle these processes.

Conceptual Approach
-------------------

Wetting and drying in ADCIRC is handled through an elemental approach. Key concepts include:

1. **Wet/Dry Status**: Each node is classified as either wet or dry. An element is considered wet only when all of its nodes are wet.

2. **Minimum Depth**: A threshold depth (:math:`H_{min}`) below which a node is considered dry.

3. **Gradual Transition**: The transition between wet and dry states is managed with smoothing algorithms to prevent numerical instabilities.

Basic Algorithm
---------------

The core wetting and drying algorithm follows these steps:

1. **Initial Classification**: 
   * At the beginning of the simulation, nodes with water depths below :math:`H_{min}` are marked as dry.
   * Nodal velocities at dry nodes are set to zero.

2. **Elemental Status Update**:
   * Elements are classified as wet (all nodes wet), dry (all nodes dry), or partially wet (mixed).
   * Partially wet elements are treated specially to ensure mass conservation and numerical stability.

3. **Dynamic Reclassification**:
   * At each time step, nodes are reclassified based on the computed water elevation.
   * A node becomes wet when :math:`\zeta + h > H_{min}`.
   * A node becomes dry when :math:`\zeta + h < H_{min}`.

4. **Velocity Treatment**:
   * When a node transitions from dry to wet, velocities are interpolated from neighboring wet nodes.
   * When a node transitions from wet to dry, velocities are gradually reduced to zero.

Numerical Techniques
--------------------

Several numerical techniques are employed to maintain stability and accuracy during wetting and drying:

### Extrapolation of Elevation

For partially wet elements, the water elevation at dry nodes is extrapolated from wet nodes using linear or higher-order extrapolation to maintain continuity.

### Velocity Handling

The momentum equations require special handling for wetting and drying:

1. **Near-dry Treatment**: When the water depth approaches :math:`H_{min}`, the bottom friction term is gradually increased to slow the flow:

   .. math::
   
       C_f' = C_f \left( \frac{H_{min}}{H} \right)^2

   This naturally retards flow as an area approaches drying.

2. **Dry Element Velocities**: In completely dry elements, the velocities are set to zero, and the momentum equations are not solved.

3. **Shallow Element Treatment**: For very shallow elements (:math:`H < 1.1 H_{min}`), a simplified momentum equation may be used with enhanced friction to prevent unrealistic velocities.

### Mass Conservation Enhancement

Special techniques ensure mass conservation during wetting and drying:

1. **Flux Limiters**: Fluxes across element boundaries are limited based on the available water volume to prevent negative depths.

2. **Volume Correction**: Small corrections are applied to ensure that the total volume of water in the domain is conserved despite discrete wetting and drying.

Implementation Details
----------------------

The implementation in ADCIRC involves several practical considerations:

### Threshold Selection

The minimum depth parameter :math:`H_{min}` is critical:
* Too small: May lead to numerical instabilities due to very shallow elements
* Too large: May artificially block water from reaching shallow areas

Typical values range from 0.01 to 0.1 meters, depending on the application and grid resolution.

### Nodal Attributes

ADCIRC allows spatial variation of wetting and drying parameters through nodal attributes:
* :math:`H_{min}` can vary spatially
* Minimum barrier heights for flow over topographic features
* Directional flow limitations for channelized flow

### Stability Considerations

To maintain stability during wetting and drying:

1. The time step may be automatically reduced when a large number of nodes are transitioning between wet and dry states.

2. Spatial filtering is applied more frequently near wetting/drying fronts to smooth potential oscillations.

3. Local modification of the tau parameter in the GWCE formulation helps maintain stability in shallow regions. 