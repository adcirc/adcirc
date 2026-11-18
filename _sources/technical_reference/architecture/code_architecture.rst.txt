.. meta::
   :description: Code Architecture in ADCIRC
   :keywords: adcirc, code architecture

Architecture
============

Code architecture comprises the "big-picture" structure and flow of the source
code. Since the ADCIRC model's source code has been under a state of refactoring
in recent years, there is not a 100% accurate source of information on the
code's current architecture. However, some `detailed diagrams of code
architecture <http://adcirc.org/home/documentation/users-manual-v50/adcirc-architecture/>`__
still capture the most important elements of the code's flow. Notable items that
are out-of-date include: modularization of code in timestep.F, fort.22 file
information is decomposed at run time instead of in adcprep, and adcpost is
generally no longer needed since model outputs can be provided globally by the
main code.



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
