
# Load needed modules
#source /opt/nec/ve/nlc/2.1.0/bin/nlcvars.sh mpi
#source /opt/nec/ve/mpi/2.12.0/bin/necmpivars.sh
#export NMPI_CXX=/opt/nec/ve/ncc/3.1.0/bin/nc++

source /opt/nec/ve/nlc/2.3.0/bin/nlcvars.sh mpi
source /opt/nec/ve/mpi/2.15.0/bin/necmpivars.sh

# MPI config
export HPCSCAN_MPI_INVOKER='mpirun -ve 0-7'
export HPCSCAN_NTHREADS=8
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

# C++ compiler
export HPCSCAN_CPP=mpinc++
export HPCSCAN_CPP_OPENACC_FLAG=
export HPCSCAN_CPP_FLAGCOMP='-O3 -g -fopenmp -w -fdiag-vector=3 -report-all'
export HPCSCAN_CPP_LIB=-lsca_openmp

# CUDA compiler
export HPCSCAN_CUDA=
export HPCSCAN_CUDA_FLAGCOMP=
export HPCSCAN_CUDA_LIB=

echo HPCSCAN set for NEC SX-Aurora with $HPCSCAN_CPP.
