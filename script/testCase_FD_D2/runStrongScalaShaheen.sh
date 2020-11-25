
source ../set_env_for_shaheen.sh

export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=32 

start_time=$(date)

n1=1000
n2=1000
n3=1000
fdOrder=8

sh ../clean_dir.sh

srun -n 1  ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder ${fdOrder} 

srun -n 2 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder ${fdOrder} -nsub2 2

srun -n 4 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder ${fdOrder} -nsub2 2 -nsub3 2

srun -n 6 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder ${fdOrder} -nsub2 3 -nsub3 2

srun -n 8 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder ${fdOrder} -nsub2 4 -nsub3 2 

perfLog='hpcscan.perf.FD_D2.log'
mv ${perfLog} ./results/runStrongScalaShaheen.${perfLog}


end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`

