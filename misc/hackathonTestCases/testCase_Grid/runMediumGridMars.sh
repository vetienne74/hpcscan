
# grid size
n1=1000
n2=1000
n3=1000

sh ../../script/clean_dir.sh

start_time=$(date)

# 32 threads
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=2
mpirun -n 1 ../../bin/hpcscan -testCase Grid -dim 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} -testMode Baseline
mpirun -n 1 ../../bin/hpcscan -testCase Grid -dim 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} -testMode GPU1
mpirun -n 1 ../../bin/hpcscan -testCase Grid -dim 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} -testMode GPU2
mpirun -n 1 ../../bin/hpcscan -testCase Grid -dim 2 -n1 ${n1} -n2 ${n2} -n3 ${n3} -testMode GPU3

cp hpcscan.perf.Grid.log runMediumGridMars.out

cat runMediumGridMars.out

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
