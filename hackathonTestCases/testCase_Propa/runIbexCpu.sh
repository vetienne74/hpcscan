
export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=32

start_time=$(date)

paramAnalysisFile=runIbexCpu.out
rm -f ${paramAnalysisFile}

if [ -z "$1" ] 
then
    list_file=configIbexCpu*.sh
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
sh clean_dir.sh

# build scripts
python3 buildConfigScriptIbexCpu.py

# run hpcscan
# loop on all script files present in the directory
for scriptFile in $( ls ${list_file}); do

    #sh clean.sh
    
    echo '* Run hpcscan with' ${scriptFile}
    rm -f hpcscan.perf.Propa.log
    sh ${scriptFile} > tmp
    cat hpcscan.perf.Propa.log >> ${paramAnalysisFile}
    
done

end_time=$(date)
tester=$(whoami)
machine=$(hostname)

echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " $machine
echo "Done by    : " $tester
