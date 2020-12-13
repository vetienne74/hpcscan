
# Load needed modules
source /opt/nec/ve/nlc/2.1.0/bin/nlcvars.sh mpi
source /opt/nec/ve/mpi/2.11.0/bin/necmpivars.sh

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
export HPCSCAN_NTHREADS=8
export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpinc++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-O3 -g -fopenmp -w -fdiag-vector=3 -report-all'
export HPCSCAN_CPP_LIB=-lsca_openmp
#export HPCSCAN_CPP_LIB=/home/lgatineau/sca_201111/libsca_openmp.a

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for Aurora with $HPCSCAN_CPP.
