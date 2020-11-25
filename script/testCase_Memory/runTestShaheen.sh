
# grid size
n1=500
n2=500
n3=500

start_time=$(date)

# from 1 to 8 threads
export KMP_AFFINITY=scatter,1,0,granularity=fine

export OMP_NUM_THREADS=1
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=2
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=4
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=8
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=12
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=16
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=20
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=24
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=28
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

export OMP_NUM_THREADS=32
srun -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
