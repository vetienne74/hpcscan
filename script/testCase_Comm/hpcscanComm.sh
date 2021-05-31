
# 8 MPI procs are considered
# all decomposition geometries are evaluated (total 10 configs)

# grid size
n1=1000
n2=1000
n3=1000

export KMP_AFFINITY=granularity=fine,compact 
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

start_time=$(date)

sh ../clean_dir.sh

# nsub1=1
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 1 -nsub2 1 -nsub3 8 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 1 -nsub2 2 -nsub3 4 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 1 -nsub2 4 -nsub3 2 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 1 -nsub2 8 -nsub3 1 -n1 ${n1} -n2 ${n2} -n3 ${n3}

# nsub1=2
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 2 -nsub2 1 -nsub3 4 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 2 -nsub2 2 -nsub3 2 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 2 -nsub2 4 -nsub3 1 -n1 ${n1} -n2 ${n2} -n3 ${n3}

# nsub1=4
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 4 -nsub2 1 -nsub3 2 -n1 ${n1} -n2 ${n2} -n3 ${n3}
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 4 -nsub2 2 -nsub3 1 -n1 ${n1} -n2 ${n2} -n3 ${n3}

# nsub1=8
$HPCSCAN_MPI_INVOKER -n 8  ../../bin/hpcscan -testCase Comm -nsub1 8 -nsub2 1 -nsub3 1 -n1 ${n1} -n2 ${n2} -n3 ${n3}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`

