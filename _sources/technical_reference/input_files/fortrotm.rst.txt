.. meta::
   :description: Fort.rotm in ADCIRC
   :keywords: adcirc, fort.rotm

.. _fortrotm:

Fort.rotm: Rotation of the geographical coordinates
===================================================

The fort.rotm file is used to rotate the geographical
coordinates on the Spherical Earth (see :ref:`ICS <ics_parameter>` for details on how to
specify that rotation is desired). Rotation is typically desired in order to
move the North Pole from the ocean onto land (and keep the South Pole on land)
so that the singularity in the Spherical coordinate form of the governing
equations is removed (see the File Format Examples section below for examples of
valid rotations that achieve this goal).

File Format
-----------

File format is short but there are three variants. General format is: first line
is a string indicating format type, following lines are the values.

#. znorth_in_spherical_coors

   lon0 lat0 ! Comment: lon0 lat0 is the center of the new north pole

#. z-x-z

   alpha beta gamma ! Comment: intrinsic rotation of alpha, beta, gamma

#. rotation_matrix

   rot11 rot12 rot13
   rot21 rot22 rot23
   rot31 rot32 rot33
   ! Comment: a rotation matrix R: x-> x' is given

File Format Examples
--------------------

The following examples ensure that both poles are rotated onto land. The results
should stay effectively the same no matter the choice.

-  znorth_in_spherical_coors

   -42.8906 72.3200 ! Greenland-Antarctica

-  znorth_in_spherical_coors

   112.8516 40.3289 ! China-Argentina

-  znorth_in_spherical_coors

   114.16991 0.77432 ! Borneo-Brazil
