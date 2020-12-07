
# Load needed modules
module load cuda/11.0.1
module load intel/2020 intelmpi/2020

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'

# OpenMP config
export HPCSCAN_NTHREADS=32
export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpiicpc
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -std=c++11 -xHost -ffast-math -fpermissive'
export HPCSCAN_CPP_FLAGLINK=-fopenmp
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=nvcc
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_FLAGLINK=
export HPCSCAN_CUDA_LIB='-L/usr/local/cuda/lib64 -lcuda -lcudart'

echo HPCSCAN set for Ibex with $HPCSCAN_CPP and $HPCSCAN_CUDA.
