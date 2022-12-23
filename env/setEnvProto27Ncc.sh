
# Load needed modules
source /opt/nec/ve/nlc/2.3.0/bin/nlcvars.sh mpi
source /opt/nec/ve/mpi/2.21.0/bin/necmpivars.sh

# MPI config
#export HPCSCAN_MPI_INVOKER='mpirun -ve 0-7'
export HPCSCAN_MPI_INVOKER='mpirun'
export HPCSCAN_NTHREADS=8
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpinc++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-O3 -g -fopenmp -w -fdiag-vector=3 -report-all'
export HPCSCAN_CPP_LIB=-lsca_openmp

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB=

# HIP compiler
export HPCSCAN_HIP=
export HPCSCAN_HIP_FLAGCOMP=
export HPCSCAN_HIP_LIB=

# display Hpcscan settings
sh ./displayHpcscanEnv.sh
