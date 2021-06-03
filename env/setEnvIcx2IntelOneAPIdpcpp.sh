
# Load needed modules
#source /raid/opt/intel/parallel_studio_xe_2020.4.912/psxevars.sh
source /raid/opt/intel/oneapi/setvars.sh
#module load cuda/11.0.1
#module load intel/2020 intelmpi/2020

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun'

# OpenMP config
export HPCSCAN_NTHREADS=64
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP='mpiicpc -cxx=dpcpp'
export HPCSCAN_CPP_OPENACC_FLAG=
#export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -std=c++11 -xHost -ffast-math -fpermissive'
export HPCSCAN_CPP_FLAGCOMP='-w -g -O3 -fopenmp -std=c++11 -xHost -ffast-math'
#export HPCSCAN_CPP_FLAGCOMP='-xICELAKE-SERVER -w -g -O3 -fopenmp -std=c++11 -qopt-zmm-usage=high -no-prec-div -no-prec-sqrt -fno-fnalias -fno-nalias -ftz -fma -fp-model fast=2  -qopt-dynamic-alig -fimf-use-svml=true -fimf-precision -qopt-report=3 -qopt-report-annotate=text -qopt-report-phase=vecq -fpermissive'
export HPCSCAN_CPP_LIB=

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for Icx2 with $HPCSCAN_CPP.
