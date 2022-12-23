
# grid size
n1=1000
n2=1000
n3=1000

export OMP_NUM_THREADS=$HPCSCAN_NTHREADS
export KMP_AFFINITY=granularity=fine,compact

start_time=$(date)

sh ../clean_dir.sh

# FD_D2 Baseline
testMode='Baseline'
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 6
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 10
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 14
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16
mv hpcscan.perf.FD_D2.log hpcscanPerfFD_D2${testMode}.log

# FD_D2 NEC
testMode='NEC'
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4 -autoPad 
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 6 -autoPad -n1Offset 1 -n2Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8 -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 10 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12 -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 14 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16 -autoPad
mv hpcscan.perf.FD_D2.log hpcscanPerfFD_D2${testMode}.log

# FD_D2 NEC_SCA
testMode='NEC_SCA'
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 2 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 4 -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 6 -autoPad -n1Offset 1 -n2Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 8 -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 10 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 12 -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 14 -autoPad -n1Offset 1
$HPCSCAN_MPI_INVOKER -n 1 ../../bin/hpcscan -testCase FD_D2 -testMode $testMode -dim 3 -n1 ${n1} -n2 ${n2} -n3 ${n3} -fdOrder 16 -autoPad
mv hpcscan.perf.FD_D2.log hpcscanPerfFD_D2${testMode}.log

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`
