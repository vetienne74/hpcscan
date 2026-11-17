# Short subset of runValidationTests.sh for one test mode, single MPI rank.
# Usage: sh runSmokeTests.sh <testMode>   (for example CUDA, CUDA_Opt or CUDA_Ref)
# Run it from the script directory with the hpcscan env set.
# Every case is checked against its reference solution and prints PASSED or FAILED.

if [ -z "$1" ]
then
    echo "usage: sh runSmokeTests.sh <testMode>"
    exit 1
fi
tM=$1

export OMP_NUM_THREADS=2
if [ "$HPCSCAN_MPI_INVOKER" = "" ]
then
    HPCSCAN_MPI_INVOKER=mpirun
fi

export report_file=runSmokeTests.${tM}.out
rm -f ${report_file}

run() {
    echo "hpcscan -testMode ${tM} $*"
    $HPCSCAN_MPI_INVOKER -n 1 ../bin/hpcscan -ntry 1 -testMode ${tM} "$@" >> ${report_file} 2>&1
    if [ $? -ne 0 ]
    then
        echo "  -> command failed"
        rc=1
    fi
}

rc=0

run -testCase Memory

for dim in 1 2 3
do
    for order in 2 8 16
    do
        run -testCase FD_D2 -dim ${dim} -fdOrder ${order}
        run -testCase FD_D1 -dim ${dim} -fdOrder ${order}
    done
    run -testCase Grid -dim ${dim}
done

run -testCase FD_D2 -dim 3 -fdOrder 8 -autoPad
run -testCase FD_D1 -dim 3 -fdOrder 8 -autoPad
run -testCase FD_D2 -dim 3 -fdOrder 16 -n1 27 -n2 36 -n3 41
run -testCase FD_D2 -n1Offset 1 -n2Offset 2 -n1AddPad 5 -n2AddPad 5 -n3AddPad 5

run -testCase Propa -dim 1 -fdOrder 8 -nt 50 -dt 0.001 -n1 100 -param1 2 -snapInc 20
run -testCase Propa -dim 2 -fdOrder 8 -nt 50 -dt 0.001 -n1 100 -n2 100 -param2 2 -param1 2 -snapInc 20
run -testCase Propa -dim 3 -fdOrder 4 -propagator Ac2SplitComp

sh clean_dir.sh
echo "report: ${report_file}"
exit ${rc}
