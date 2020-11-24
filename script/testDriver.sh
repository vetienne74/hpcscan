
#------------------------------------------------------------------------------------
# this script is used to assess good behavior of Xbenchmark in several configurations
# it is not intended to be used for performance measurements
# the test cases are very light and all can be run on a laptop within few minutes
# OMP_NUM_THREADS is set to 2 to allow efficient execution on every platforms
#------------------------------------------------------------------------------------

export OMP_NUM_THREADS=2

echo '------------------------------------------------------------------------'
echo Start testDriver.sh

# if not specified, testMode is set to Baseline
if [ -z $1 ]
then
    tM="Baseline"
else
    tM=$1
fi

echo '>>> Test Mode' ${tM} '<<<'

# if not specified, use mpirun to launch MPI
if [ -z $XBENCHMARK_MPI_INVOKER ]
then
    XBENCHMARK_MPI_INVOKER=mpirun
fi
echo Launch MPI with $XBENCHMARK_MPI_INVOKER

# run all test cases with default parameters
echo Running -testCase All ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase All >> ${report_file}

# test case Util
echo Running -testCase Util ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub2 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub2 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub1 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub1 2 -n1AddPad 5 -n2AddPad 5 -n3AddPad 5 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 -dim 1 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 -dim 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub2 2 -dim 2 >> ${report_file}

# test case Template
echo Running -testCase Template ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Template >> ${report_file}

# test case Memory
echo Running -testCase Memory ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Memory >> ${report_file}

# test case FD_D2
echo Running -testCase FD_D2 ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 1 -n1 1000 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 1 -n1 1000 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 1 -n1 1000 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 1 -n1 1000 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 1 -n1 1000 >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 2 -n1 20 -n2 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 2 -n1 20 -n2 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 2 -n1 20 -n2 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 2 -n1 20 -n2 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 2 -n1 20 -n2 20 >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 3 -n1 20 -n2 20 -n3 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 3 -n1 20 -n2 20 -n3 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 3 -n1 20 -n2 20 -n3 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 3 -n1 20 -n2 20 -n3 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 3 -n1 20 -n2 20 -n3 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 3 -n1 27 -n2 36 -n3 41 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1AddPad 5 -n2AddPad 5 -n3AddPad 5 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1MulPad 5 -n2MulPad 5 -n3MulPad 5 >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1MulPad 5 -n2MulPad 5 -n3MulPad 5 -autoPad >> ${report_file}

# test case Grid
echo Running -testCase Grid ...
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 1 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 1 -nsub1 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 -nsub2 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 -nsub2 2 -nsub1 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 -nsub3 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 -nsub3 2 -nsub2 2 >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 1 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 >> ${report_file}

# test case Comm
echo Running -testCase Comm ...
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub1 2 -n1 200 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 3 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub1 3 -n1 300 >> ${report_file}

# test case Propa
echo Running -testCase Propa ...

# -propaName Ac2Standard (default)
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 -snapInc 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 -param2 2 -param1 2 -snapInc 20 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 -param3 2 -snapInc 20 >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 -snapInc 20 -nsub1 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 -param2 2 -param1 2 -snapInc 20 -nsub2 2 >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 -param3 2 -snapInc 20 -nsub3 2 >> ${report_file}

# -propaName Ac2SplitComp
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 -snapInc 20 -propaName Ac2SplitComp >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 -param2 2 -param1 2 -snapInc 20 -propaName Ac2SplitComp >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 -param3 2 -snapInc 20 -propaName Ac2SplitComp >> ${report_file}

$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 -snapInc 20 -nsub1 2 -propaName Ac2SplitComp >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 -param2 2 -param1 2 -snapInc 20 -nsub2 2 -propaName Ac2SplitComp >> ${report_file}
$XBENCHMARK_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 -param3 2 -snapInc 20 -nsub3 2 -propaName Ac2SplitComp >> ${report_file}

echo End testDriver.sh
echo '------------------------------------------------------------------------'
