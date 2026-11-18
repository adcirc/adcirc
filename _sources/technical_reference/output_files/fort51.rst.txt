Fort.51: Elevation Harmonic Constituents at Specified Elevation Recording Stations
==================================================================================

Amplitude and phase information obtained from harmonic analysis of the elevation solution at the elevation recording stations and for the frequencies specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`NFREQ <NFREQ>`
   for k = 1, :ref:`NFREQ <NFREQ>`
      :ref:`HAFREQ(k) <HAFREQ>`, :ref:`HAFF(k) <HAFF>`, :ref:`HAFACE(k) <HAFACE>`, :ref:`NAMEFR(k) <NAMEFR>`
   end k loop

   :ref:`NSTAE <NSTAE>`
   for k=1, :ref:`NSTAE <NSTAE>`
      for j=1, :ref:`NFREQ <NFREQ>`
         :ref:`EMAG(k,j) <EMAG>`, :ref:`PHASEDE(k,j) <PHASEDE>`
      end j loop
   end k loop

Note
----

* Output format is ascii 