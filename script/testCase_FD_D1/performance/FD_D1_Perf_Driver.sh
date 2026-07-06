
# grid size
n1=1000
n2=1000
n3=1000

export OMP_NUM_THREADS=$HPCSCAN_NTHREADS
export KMP_AFFINITY=granularity=fine,compact

# if not specified, testMode is set to Baseline
if [ -z $1 ]
then
    tM="Baseline"
else
    tM=$1
fi

echo '>>> Running with Test Mode' ${tM} '... <<<'

# if not specified, use mpirun to launch MPI
if [ "$HPCSCAN_MPI_INVOKER" = "" ]
then
    HPCSCAN_MPI_INVOKER=mpirun
fi
#echo Launch MPI with $HPCSCAN_MPI_INVOKER

#==========================================================================================================
# test case FD D1 Perf
#==========================================================================================================
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 6 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 10 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 14 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ${HPCSCAN_DIR}/bin/hpcscan -testCase FD_D1 -testMode $tM -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16 >> ${report_file}
mv hpcscan.perf.FD_D1.log hpcscanPerfFD_D1${tM}.log
