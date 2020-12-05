
# MPI config
export HPCSCAN_MPI_INVOKER="mpirun --oversubscribe"

# OpenMP config
export HPCSCAN_NTHREADS=1
export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -mavx2'
export HPCSCAN_CPP_FLAGLINK=-fopenmp
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_FLAGLINK=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for Mars with $HPCSCAN_CPP.
