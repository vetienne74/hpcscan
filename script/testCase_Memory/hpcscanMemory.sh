
# grid size
n1=1000
n2=1000
n3=1000

start_time=$(date)

sh ../clean_dir.sh

# from 1 to 8 threads
export KMP_AFFINITY=scatter,1,0,granularity=fine

export OMP_NUM_THREADS=1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=2
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=4
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=8
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=12
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=16
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=20
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=24
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=28
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=32
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
