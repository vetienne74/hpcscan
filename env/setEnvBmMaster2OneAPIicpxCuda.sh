
# Load needed modules
. /opt/intel/oneapi/setvars.sh --force
export PATH=/usr/local/cuda/bin:$PATH

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
export HPCSCAN_NTHREADS=36
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP='mpiicpc -cxx=icpx'
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -xHost -ffast-math -I /usr/local/cuda/include/'
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=nvcc
export HPCSCAN_CUDA_FLAGCOMP='--std c++11'
export HPCSCAN_CUDA_LIB='-L/usr/local/cuda/lib64 -lcuda -lcudart -lnvidia-ml'

# HIP compiler
export HPCSCAN_HIP=
export HPCSCAN_HIP_FLAGCOMP=
export HPCSCAN_HIP_LIB=

# display Hpcscan settings
sh ./displayHpcscanEnv.sh
