
# Load needed modules
# no module needed on Mars

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'

# OpenMP config
export HPCSCAN_NTHREADS=1
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG='-fopenacc'
export HPCSCAN_CPP_FLAGCOMP='-g -O3 -mavx2 -fopenmp'
export HPCSCAN_CPP_LIB=

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
