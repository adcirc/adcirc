# 09/13/2019 - COASTAL ACT Project
# This is a script to set the environment for a typical ESMF appliction
# How to use: from linux command line type: 
# source esmf-impi-env.sh
#

# Path to libraries and includes and binaries on Hera
# /apps/intel/compilers_and_libraries_2018.5.274/linux/bin/intel64/ifort
# ifor version: ifort (IFORT) 18.0.5 20180823
module load intel/18.0.5.274

# /apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/lib64
# ifor version: ifort (IFORT) 18.0.5 20180823
# mpi: libmpifort.so.12.0
module load impi/2018.0.4

# /apps/szip/2.1/lib
module load szip/2.1

# /apps/hdf5/1.10.4/intel/16.1.150/impi/5.1.2.150/
module load hdf5/1.10.4

# /apps/netcdf/4.6.1/intel/16.1.150/lib
module load netcdf/4.6.1



# Environment for ESMF v8.0.0 beta snapshot 48g
# /scratch1/NCEPDEV/nems/emc.nemspara/soft/esmf/8.0.0bs48g-intel18.0.5.274-impi2018.0.4-netcdf4.6.1
module use /home/emc.nemspara/SOFT-hera/modulefiles
module load esmf/8.0.0bs48g

# Alternative Environment for ESMF instead of use of module
# ESMF must have been compiled with these options:
# export ESMF_BOPT='g'
# export ESMF_COMM=intelmpi      # acceptable MPIs: mpich, mpich2, lam, openmpi or intelmpi
# export ESMF_COMPILER='intel'
# export ESMF_ABI='64'
# export ESMF_OS='Linux'         # uname -s
#
# Un comment below line and comment out above two lines for ESMF v8.0.0 modules 
# export ESMFMKFILE=/scratch1/NCEPDEV/nems/emc.nemspara/soft/esmf/8.0.0bs48g-intel18.0.5.274-impi2018.0.4-netcdf4.6.1/lib/esmf.mk
