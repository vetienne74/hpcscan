#!/bin/bash
#SBATCH --account=k1256
#SBATCH --job-name=hpcscan
#SBATCH --output=hpcscan.out
#SBATCH --error=hpcscan.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00

source ../../env/setEnvShaheenIntel.sh

sh runValidationTests.sh
