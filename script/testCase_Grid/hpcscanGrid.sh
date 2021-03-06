
# grid size
n1=1000
n2=1000
n3=1000

export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

start_time=$(date)

sh ../clean_dir.sh

# grid size
n1=500
n2=500
n3=500
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Grid -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

# grid size
n1=1000
n2=1000
n3=1000
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Grid -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`

