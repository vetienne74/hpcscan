#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00

source ../../env/setEnvShaheenIntel.sh

script='hpcscanGrid'
sh ${script}.sh > ${script}Shaheen.out

perfLog='hpcscan.perf.Grid.log'
mv ${perfLog} ${script}Shaheen.log
