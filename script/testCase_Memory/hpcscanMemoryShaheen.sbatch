#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=8
#SBATCH --time=00:30:00

source ../../env/setEnvShaheenIntel.sh

sh hpcscanMemory.sh > hpcscanCommMemory.out

perfLog='hpcscan.perf.Memory.log'
mv ${perfLog} ./hpcscanMemory.${perfLog}
