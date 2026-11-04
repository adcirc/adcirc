.. _fort10:

Fort.10: Passive Scalar Transport Input File
============================================

File Structure
--------------

.. parsed-literal::

   Header Line 1
   Header Line 2
   :ref:`NVP <NVP>`
   for k=1 to :ref:`NVP <NVP>`
      :ref:`jki <jki>`, :ref:`DACONC(jki) <DACONC>`
   end k loop

File structure for a 3D run:

.. parsed-literal::

   Header Line 1
   Header Line 2
   :ref:`NVN <NVN>`, :ref:`NVP <NVP>`
   for k=1 to :ref:`NVP <NVP>`
      for j=1 to :ref:`NVN <NVN>`
         :ref:`NHNN <NHNN>`, :ref:`NVNN <NVNN>`, :ref:`CONC(NHNN,NVNN) <CONC>`
      end j loop
   end k loop 