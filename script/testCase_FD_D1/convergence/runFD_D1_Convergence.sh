
export KMP_AFFINITY=granularity=fine,compact 
export OMP_NUM_THREADS=$HPCSCAN_NTHREADS

start_time=$(date)
machine=`hostname`
today=`date +%F`

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

cp hpcscan.hwCounter.Propa.log propa.hwCounter.${machine}.${today}.log
cp hpcscan.perf.Propa.log propa.perf.${machine}.${today}.log

echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " $machine
echo "Done by    : " $tester
