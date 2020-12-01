
# environment for Ibex

module load pgi/20.1
module load cuda/11.0.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/csgv/pgi/20.1/el7.7_binary/linux86-64/20.1/mpi/openmpi-3.1.3/lib
export PATH=/sw/csgv/pgi/20.1/el7.7_binary/linux86-64/20.1/mpi/openmpi-3.1.3/bin:$PATH
export PGI_ACC_NOTIFY=2
export HPCSCAN_MPI_INVOKER="srun"
export HPCSCAN_NTHREADS=1

alias sq='squeue -u x_etiennv'

export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

echo My environment is set for Ibex acc!
