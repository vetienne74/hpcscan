
# environment for Ibex

module load cuda/11.0.1
module load intel/2020 intelmpi/2020

export HPCSCAN_MPI_INVOKER="srun"
export HPCSCAN_NTHREADS=32

alias sq='squeue -u x_etiennv'

export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

echo My environment is set for Ibex!
