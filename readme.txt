
To run parallel ADCIRC on the UNC cluster:
written by RL 6/05/03


1.)  get access to the Intel Fortran & C compilers and MPICH-GM libraries by typing imp add intel_fortran intel_CC mpich-gm
     this need only to be done once, when a user account is first set up.  Alternative compilers and MPI libraries can be used
     as outlined at: http://www.unc.edu/atn/asg/scitool/parallel_setup.html

2.)  in directory work, file cmplrflags.mk, uncomment the LINUX-Intel section marked for the UNC cluster

3.)  in directory work, type:  make all
     This should have compiled and created 5 executable files: adcprep, adcprep2, adcpost, adcirc, padcirc

4.)  to run a serial ADCIRC job, create a directory with the required input files (e.g., fort.14, fort.15, etc.) and copy
     in adcirc.  For a very short job, simply type adcirc to run interactively.  For a longer job, submit a batch run:
     bsub -o screenout.out -q batch adcirc

5.)  to run a paralle ADCIRC job, create a directory with the required input files (e.g., fort.14, fort.15, etc) and copy
     in adcprep, adcpost, padcirc.  Run adcprep (interactively) and respond to various interactive questions.  Submit a batch run:
     bsub -o screenout.out -a mpich_gm -q parallel -n X padcirc   (Note, the X should be replaced by the number of processors 
     to be used).  For a very short run you can use the express queue rather than the parallel queue.  Once the job is completed,
     run adcpost (interactively) to gather output from each processor into a single, global set of files. 



Note: the cmplrflags.mk in this directory works on the UNC cluster