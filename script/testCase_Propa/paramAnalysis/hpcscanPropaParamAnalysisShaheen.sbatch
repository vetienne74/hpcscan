#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00

source ../../../env/setEnvShaheenIntel.sh

script='hpcscanPropaParamAnalysis'
sh ${script}.sh > ${script}Shaheen.out

perfLog='hpcscan.perf.Propa.log'
mv ${perfLog} ${script}Shaheen.log
