This is modified from v45.02 

The code is compiled for a windows platform.

adcirc_sp.exe is the single precision executable
adcirc_dp.exe is the double precision executable

iterative solver convergence criteria no larger than 1x10e-5 for single precision,
  1x10e-8 for double precision. Note, the smaller the criterion the larger the number
  of iterations and run time.  Also, it may be necessary to increase the maximum number
  of iterations per time step to >25, depending on the tolerance that has been set.  You
  can judge this based on the number of iterations that is output to the screen every 
  time step.


If wetting and drying is activated, I recommend using:
  hmin = 0.1
  velmin = 0.1
  nodedrymin & nodewetmin = 0 NO LONGER ACTIVE PARAMETERS
