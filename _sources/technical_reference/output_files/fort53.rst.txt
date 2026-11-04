.. _fort53:

Fort.53: Elevation Harmonic Constituents at All Nodes in the Model Grid
=======================================================================

Amplitude and phase information obtained from harmonic analysis of the elevation solution at all nodes in the grid for the frequencies specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`NFREQ <NFREQ>`
   for k = 1, :ref:`NFREQ <NFREQ>`
      :ref:`HAFREQ(k) <HAFREQ>`, :ref:`HAFF(k) <HAFF>`, :ref:`HAFACE(k) <HAFACE>`, :ref:`NAMEFR(k) <NAMEFR>`
   end k loop

   :ref:`NP <NP>`
   for k=1, :ref:`NP <NP>`
      k
      for j=1, :ref:`NFREQ <NFREQ>`
         :ref:`EMAGT(k,j) <EMAGT>`, :ref:`PHASEDE(k,j) <PHASEDE>`
      end j loop
   end k loop

Note
----

* Output format is ascii 