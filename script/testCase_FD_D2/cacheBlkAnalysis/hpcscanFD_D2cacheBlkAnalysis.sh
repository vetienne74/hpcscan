
export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

start_time=$(date)

if [ -z "$1" ] 
then
    list_file=config*.sh
else
    list_file=$1
fi

# if not specified, use mpirun to launch MPI
if [ -z $HPCSCAN_MPI_INVOKER ]
then
    export HPCSCAN_MPI_INVOKER=mpirun
fi
echo Launch MPI with $HPCSCAN_MPI_INVOKER

# clean dir
sh ./clean_dir.sh
rm -f hpcscanFD_D2cacheBlkAnalysis.out

# build scripts
python3 buildConfigScript.py

# run hpcscan
# loop on all script files present in the directory
for scriptFile in $( ls ${list_file}); do

    #sh clean.sh
    
    echo '* Run hpcscan with' ${scriptFile}
    sh ${scriptFile}
    
done

end_time=$(date)
tester=$(whoami)
machine=$(hostname)

echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " $machine
echo "Done by    : " $tester
