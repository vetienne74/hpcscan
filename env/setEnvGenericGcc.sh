
# This is an generic script that needs to be sourced on your system prior to use hpcscan:
# - compilation with gcc
# - execution only on CPU

# Add below the modules that needed to be loaded for compilation on your system
# module load xxx/xxx

# MPI config
#export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
# define the number of threads to be used
export HPCSCAN_NTHREADS=16
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG=
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

