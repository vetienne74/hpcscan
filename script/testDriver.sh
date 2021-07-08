
#------------------------------------------------------------------------------------
# This script is used to assess good behavior of hpcscan in several configurations
# It is not intended to be used for performance measurements
# The test cases are very light and all can run on a laptop within few minutes
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
if [ "$HPCSCAN_MPI_INVOKER" = "" ]
then
    HPCSCAN_MPI_INVOKER=mpirun
fi
echo Launch MPI with $HPCSCAN_MPI_INVOKER

#==========================================================================================================
# test case Util
#==========================================================================================================

echo Running -testCase Util ...
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 4 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 6 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 8 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 10 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 12 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 14 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -fdOrder 16 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub2 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub2 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub1 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub3 2 -nsub1 2 -n1AddPad 5 -n2AddPad 5 -n3AddPad 5 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub1 2 -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Util -nsub2 2 -dim 2 >> ${report_file}

#==========================================================================================================
# test case Template
#==========================================================================================================

echo Running -testCase Template ...
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Template >> ${report_file}

#==========================================================================================================
# test case Memory
#==========================================================================================================

echo Running -testCase Memory ...
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Memory >> ${report_file}

#==========================================================================================================
# test case FD_D2
#==========================================================================================================

echo Running -testCase FD_D2 ...

# 1D all orders
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 6  -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 10 -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 14 -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 1 >> ${report_file}

# 2D all orders
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 6  -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 10 -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 14 -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 2 >> ${report_file}

# 3D all orders
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 6  -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 10 -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 14 -dim 3 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 3 >> ${report_file}

# 3D all orders with -autoPad
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 2  -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 4  -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 6  -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 8  -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 10 -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 12 -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 14 -dim 3 -autoPad >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 3 -autoPad >> ${report_file}

# custom grid size
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -fdOrder 16 -dim 3 -n1 27 -n2 36 -n3 41 >> ${report_file}

# Padding
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1AddPad 5 -n2AddPad 5 -n3AddPad 5 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1MulPad 5 -n2MulPad 5 -n3MulPad 5 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1MulPad 5 -n2MulPad 5 -n3MulPad 5 >> ${report_file}

# Offset
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1Offset 1 -n2Offset 2 -n3Offset 3 >> ${report_file}

# Offset + Padding
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase FD_D2 -n1Offset 1 -n2Offset 2 -n3Offset 3 \
		     -n1AddPad 5 -n2AddPad 5 -n3AddPad 5 >> ${report_file}

#==========================================================================================================
# test case Grid
#==========================================================================================================

echo Running -testCase Grid ...
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 1 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 1 -nsub1 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 -nsub2 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 2 -nsub2 2 -nsub1 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 -nsub3 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 4 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Grid -dim 3 -nsub3 2 -nsub2 2 >> ${report_file}

#==========================================================================================================
# test case Comm
#==========================================================================================================

echo Running -testCase Comm ...
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub1 2 -n1 200 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 3 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub1 3 -n1 300 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub2 2 -n2 200 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 3 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub2 3 -n2 300 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub3 2 -n3 200 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 3 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Comm -dim 3 -nsub3 3 -n3 300 >> ${report_file}

#==========================================================================================================
# test case Propa
#==========================================================================================================

echo Running -testCase Propa ...

# -propagator Ac2Standard (default)
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 2  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 4  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 6  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 8  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 10 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 12 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 14 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 16 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 2  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 4  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 6  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 8  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 10 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 12 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 14 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 16 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 2  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 4  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 6  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 8  >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 10 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 12 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 14 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 16 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 \
		     -snapInc 20 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 \
		     -param2 2 -param1 2 -snapInc 20 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 -param3 2 \
		     -snapInc 20 >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 \
		     -snapInc 20 -nsub1 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 \
		     -param2 2 -param1 2 -snapInc 20 -nsub2 2 >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 \
		     -param3 2 -snapInc 20 -nsub3 2 >> ${report_file}

# -propagator Ac2SplitComp
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 2  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 4  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 6  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 8  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 10 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 12 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 14 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -fdOrder 16 -propagator Ac2SplitComp >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 2  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 4  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 6  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 8  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 10 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 12 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 14 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -fdOrder 16 -propagator Ac2SplitComp >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 2  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 4  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 6  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 8  -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 10 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 12 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 14 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -fdOrder 16 -propagator Ac2SplitComp >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -param1 2 \
		     -snapInc 20 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 -n2 100 \
		     -param2 2 -param1 2 -snapInc 20 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 \
		     -param3 2 -snapInc 20 -propagator Ac2SplitComp >> ${report_file}

$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 1 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 \
		     -param1 2 -snapInc 20 -nsub1 2 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 2 -nt 50 -dt 0.001 -fdOrder 8 -n1 100 \
		     -n2 100 -param2 2 -param1 2 -snapInc 20 -nsub2 2 -propagator Ac2SplitComp >> ${report_file}
$HPCSCAN_MPI_INVOKER -n 2 ../bin/hpcscan -ntry 1 -testMode ${tM} -testCase Propa -dim 3 -nt 50 -dt 0.001 -param1 2 -param1 2 \
		     -param3 2 -snapInc 20 -nsub3 2 -propagator Ac2SplitComp >> ${report_file}

echo End testDriver.sh
echo '------------------------------------------------------------------------'
