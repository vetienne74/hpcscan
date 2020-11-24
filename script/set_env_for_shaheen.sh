
# environment for Shaheen

module swap PrgEnv-cray PrgEnv-intel

export HPCSCAN_MPI_INVOKER="srun"
export HPCSCAN_NTHREADS=32

alias sq='squeue -u x_etiennv'

export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

echo My environment is set for Shaheen!
