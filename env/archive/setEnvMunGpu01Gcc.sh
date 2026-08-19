
# Load needed modules
# no module needed on MunGpu01

#export PATH=/opt/mpi/ompi/bin:$PATH
export PATH=/opt/mvapich2/gdr/2.3.5/mcast/no-openacc/rocm4.1.1/mofed5.0/mpirun/gnu8.3.0/bin/:$PATH

# MPI config
#export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
export HPCSCAN_NTHREADS=128
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-g -O3 -mavx2 -fopenmp'
export HPCSCAN_CPP_LIB='-L/opt/rocm/lib -lamdhip64'

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
