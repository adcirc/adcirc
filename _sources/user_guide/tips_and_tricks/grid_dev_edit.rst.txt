.. meta::
   :description: Grid Development and Editing in ADCIRC
   :keywords: adcirc, grid development and editing

.. _grid_dev_edit:

Grid Development and Editing
============================

A mesh in the form of a :ref:`fort.14 file <fort14>` is required to run
ADCIRC. This page presents basic info on options available to mesh construction
end editing. For information on viewing meshes or other ADCIRC data, see
:doc:`visualization`. In many cases, a mesh may already exist to
suite one's needs. Given the complexities and challenges associated with mesh
construction, users should carefully weigh the merits of using an existing mesh,
revising an existing mesh, or building a new one.

.. _tools_for_adcirc_meshing:

Tools for ADCIRC Meshing
------------------------

SMS
~~~

SMS is a commercial Windows-based graphical program for creating, editing,
viewing, and running various models, including ADCIRC. It is popular with ADCIRC
users due to its features and support for ADCIRC file formats. More information
is available from Aquaveo
`here <https://www.aquaveo.com/software/sms-adcirc>`__ and
`here <https://xmswiki.com/wiki/SMS:ADCIRC>`__.

OceanMesh2D
~~~~~~~~~~~

Completely open source (GNU GPLv3.0) MATLAB-based software for automated
triangular mesh generation. The latest version is available from
https://github.com/CHLNDDEV/OceanMesh2D  [1]_  [2]_  [3]_.

It can be used as an end-to-end pre-processor for ADCIRC, including:

-  Generating the mesh.
-  Checking for and editing the mesh to satisfy Courant constraints.
-  Generating :ref:`fort.13 file <fort13>` attributes.
-  Generating :ref:`fort.15 file <fort15>`, including automatic generation of
   tidal potential information and tidal elevation boundary conditions.
-  Generating other input files such as :ref:`fort.24 file <fort24>`.
-  Storing the mesh and its attributes into a ``msh`` class container that can
   be saved as an efficient .mat binary.
-  Plotting the mesh and its attributes using the generalized ``plot`` command.
-  Writing the from the ``msh`` class container into the ASCII fort.xx input
   files.

References
==========

.. raw:: html

   <references />

.. [1]
   Pringle, W.J., Roberts, K.J., 2020. CHLNDDEV/OceanMesh2D: OceanMesh2D V3.0.0.
   doi:10.5281/zenodo.3721137

.. [2]
   Roberts, K.J., Pringle, W.J., Westerink, J.J., 2019. OceanMesh2D 1.0:
   MATLAB-based software for two-dimensional unstructured mesh generation in
   coastal ocean modeling. Geosci. Model Dev. 12, 1847–1868.
   doi:10.5194/gmd-12-1847-2019

.. [3]
   Roberts, K.J., Pringle, W.J., 2018. OceanMesh2D: User guide - Precise
   distance-based two-dimensional automated mesh generation.
   doi:10.13140/RG.2.2.21840.61446/1



.. raw:: html

   <style>
   .wrap-table th, .wrap-table td {
     white-space: normal !important;
     word-wrap: break-word !important;
     max-width: 100% !important;
     overflow-wrap: break-word !important;
     hyphens: auto !important;
   }
   </style>
