# Environment for a GitHub Actions runner:
# - gcc + OpenMPI (single rank)
# - CUDA with PTX-only code for the simulated H100 (compute_90)
# Source this file from the env directory.

export HPCSCAN_MPI_INVOKER='mpirun'

export HPCSCAN_NTHREADS=2
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

export HPCSCAN_CPP=mpic++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP="-g -O2 -fopenmp -I${CUDA_HOME}/include"
export HPCSCAN_CPP_LIB=

# CUDA_HOME must point to the CUDA toolkit
export HPCSCAN_CUDA=nvcc
# nvcc also needs the MPI headers (Open MPI: mpic++ --showme:compile)
export HPCSCAN_CUDA_FLAGCOMP="-arch=compute_90 -cudart shared $(mpic++ --showme:compile)"
export HPCSCAN_CUDA_LIB="-L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs -lcuda -lcudart -lnvidia-ml"

export HPCSCAN_HIP=
export HPCSCAN_HIP_FLAGCOMP=
export HPCSCAN_HIP_LIB=

export HPCSCAN_CUSTOM_PATH=
export HPCSCAN_CUSTOM_OPT=

sh ./displayHpcscanEnv.sh
