
# Load needed modules
. /opt/intel/oneapi/setvars.sh --force

# MPI config
#export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
export HPCSCAN_NTHREADS=2
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP='mpiicpc -cxx=dpcpp'
export HPCSCAN_CPP_OPENACC_FLAG=
#export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -std=c++11 -xHost -ffast-math'
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -xHost -ffast-math'
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for Neptune with:
$HPCSCAN_CPP --version

