
# grid size
n1=1000
n2=1000
n3=1000

sh ../clean_dir.sh

# set 32 threads
export OMP_NUM_THREADS=32
export KMP_AFFINITY=scatter,1,0,granularity=fine
start_time=$(date)

# FD_D2 (Baseline)
srun -n 1 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2
srun -n 1 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4
srun -n 1 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8
srun -n 1 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12
srun -n 1 ../../bin/hpcscan -testCase FD_D2 -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16
perfLog='hpcscan.perf.FD_D2.log'
mv ${perfLog} ./results/runMediumGridShaheen.${perfLog}

# FD_D2_CacheBlk
srun -n 1 ../../bin/hpcscan -testCase FD_D2_CacheBlk -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2
srun -n 1 ../../bin/hpcscan -testCase FD_D2_CacheBlk -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4
srun -n 1 ../../bin/hpcscan -testCase FD_D2_CacheBlk -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8
srun -n 1 ../../bin/hpcscan -testCase FD_D2_CacheBlk -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12
srun -n 1 ../../bin/hpcscan -testCase FD_D2_CacheBlk -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16
perfLog='hpcscan.perf.FD_D2_CacheBlk.log'
mv ${perfLog} ./results/runMediumGridShaheen.${perfLog}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
