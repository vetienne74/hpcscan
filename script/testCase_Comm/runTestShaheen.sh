
source ../set_env_for_shaheen.sh

export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=32 

start_time=$(date)

export n1=1000
export n2=1000
export n3=1000

sh ../clean_dir.sh

srun -n 8 ../../bin/hpcscan -testCase Comm -nsub2 4 -nsub3 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} 

srun -n 8 ../../bin/hpcscan -testCase Comm -nsub2 2 -nsub3 4 -n1 ${n1} -n2 ${n2} -n3 ${n3} 

srun -n 8  ../../bin/hpcscan -testCase Comm -nsub1 2 -nsub2 2 -nsub3 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} 

perfLog='hpcscan.perf.Comm.log'
mv ${perfLog} ./results/runTestShaheen.${perfLog}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`

