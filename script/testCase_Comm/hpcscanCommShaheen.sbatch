#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=8
#SBATCH --time=00:30:00

source ../../env/setEnvShaheenIntel.sh

script='hpcscanComm'
sh ${script}.sh > ${script}Shaheen.out

perfLog='hpcscan.perf.Comm.log'
mv ${perfLog} ${script}Shaheen.log
