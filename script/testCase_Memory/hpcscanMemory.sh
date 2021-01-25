
# grid size
n1=1000
n2=1000
n3=1000

start_time=$(date)

sh ../clean_dir.sh

# set thread affinity
export KMP_AFFINITY=scatter,1,0,granularity=fine

nthread=1
while [ $nthread -le $HPCSCAN_NTHREADS ]
do
    echo "nthread $nthread"
    export OMP_NUM_THREADS=$nthread
    # run hpcscan
    $HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase Memory -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3}
    
    if [ $nthread -lt 8 ]
    then
	# from 1 to 8, increment is 1
	nthread=$(( $nthread + 1 ))
    else
	# from 8 to $HPCSCAN_NTHREADS, increment is 4
	nthread=$(( $nthread + 4 ))
    fi
done

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
