
# Load needed modules
module swap PrgEnv-cray PrgEnv-intel

# MPI config
export HPCSCAN_MPI_INVOKER='srun'

# OpenMP config
export HPCSCAN_NTHREADS=32
export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=CC
export HPCSCAN_CPP_FLAGCOMP='-O3 -g -openmp'
export HPCSCAN_CPP_FLAGLINK=-openmp
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_FLAGLINK=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for Shaheen with $HPCSCAN_CPP.
