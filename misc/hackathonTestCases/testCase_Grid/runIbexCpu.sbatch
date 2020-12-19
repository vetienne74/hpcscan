#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J hpcscan
#SBATCH -o hpcscan.%J.out
#SBATCH -e hpcscan.%J.err
#SBATCH --time=00:30:00
##SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
##SBATCH --cpus-per-task=32
#SBATCH --constraint=[v100]
#SBATCH --res=hackathon2020

hostname
source ../../script/set_env_for_ibex.sh
export KMP_AFFINITY=scatter,1,0,granularity=fine
export OMP_NUM_THREADS=32
mpirun -n 1 ../../bin/hpcscan -testCase Grid -dim 3 -n1 500 -n2 500 -n3 500 -testMode Baseline

