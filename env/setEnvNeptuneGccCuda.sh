
# Load needed modules
# no module needed on Neptune

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'

# OpenMP config
export HPCSCAN_NTHREADS=2
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -mavx2 -fopenmp'
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=nvcc
export HPCSCAN_CUDA_FLAGCOMP='-gencode arch=compute_50,code=sm_50 -I /usr/lib/x86_64-linux-gnu/openmpi/include/'
export HPCSCAN_CUDA_LIB='-L/usr/local/cuda/lib64 -lcuda -lcudart'

echo HPCSCAN set for Neptune with $HPCSCAN_CPP and $HPCSCAN_CUDA.