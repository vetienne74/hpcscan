
# Load needed modules
module load pgi/20.1
module load cuda/11.0.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/csgv/pgi/20.1/el7.7_binary/linux86-64/20.1/mpi/openmpi-3.1.3/lib
export PATH=/sw/csgv/pgi/20.1/el7.7_binary/linux86-64/20.1/mpi/openmpi-3.1.3/bin:$PATH
#export PGI_ACC_NOTIFY=2

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun --oversubscribe'

# OpenMP config
export HPCSCAN_NTHREADS=32
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG='-acc -ta=tesla:cc70'
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -std=c++11 -mp -Minfo=accel'
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=nvcc
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB='-L/usr/local/cuda/lib64 -lcuda -lcudart -lnvidia-ml'

# HIP compiler
export HPCSCAN_HIP=
export HPCSCAN_HIP_FLAGCOMP=
export HPCSCAN_HIP_LIB=

# display Hpcscan settings
sh ./displayHpcscanEnv.sh
