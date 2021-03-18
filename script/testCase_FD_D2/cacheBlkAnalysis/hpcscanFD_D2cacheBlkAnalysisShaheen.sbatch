#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00

source ../../../env/setEnvShaheenIntel.sh

script='hpcscanFD_D2cacheBlkAnalysis'
sh ${script}.sh > ${script}Shaheen.out

perfLog='hpcscan.perf.FD_D2.log'
mv ${perfLog} ${script}Shaheen.log
