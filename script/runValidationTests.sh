
if [ -f $HPCSCAN_CPP ]
then
    echo It seems hpcscan environment has not been set.
    echo Source one of the env files in ../env or create a new one for your machine.
    exit
fi

# get testMode from command line
if [ -z "$1" ]
then
    testMode=All
else
    testMode=$1
fi

export tester=`whoami`
export machine=`hostname`
export report_file=runValidationTests.${machine}.${tester}.out

sh clean_dir.sh
rm -f ${report_file}

echo "\n========================================================================"
echo "                         START VALIDATION TESTS"
echo "========================================================================\n"
start_time=$(date)

if [ $testMode = "All" ]
then
    echo "--- All enabled testModes will be tested ---\n"
else
    echo "--- Only $testMode will be tested ---\n"
fi

if [ $testMode = "All" ] || [ $testMode = "Baseline" ]
then
    # run all tests with -testMode Baseline
    sh testDriver.sh Baseline
fi

if [ $testMode = "All" ] || [ $testMode = "CacheBlk" ]
then
    # run all tests with -testMode CacheBlk
    sh testDriver.sh CacheBlk
fi

# run all tests with -testMode OpenAcc
#sh testDriver.sh OpenAcc

if [ "$HPCSCAN_CUDA" = "nvcc" ]
then
    if [ $testMode = "All" ] || [ $testMode = "CUDA" ]
    then
	# run all tests with -testMode CUDA
	sh testDriver.sh CUDA
    else
	echo 'SKIP testMode CUDA'
    fi

    if [ $testMode = "All" ] || [ $testMode = "CUDA_Opt" ]
    then
	# run all tests with -testMode CUDA_Opt
	sh testDriver.sh CUDA_Opt
    else
	echo 'SKIP testMode CUDA_Opt'
    fi

    if [ $testMode = "All" ] || [ $testMode = "CUDA_Ref" ]
    then
	# run all tests with -testMode CUDA_Ref
	sh testDriver.sh CUDA_Ref
    else
	echo 'SKIP testMode CUDA_Ref'
    fi
else
    echo 'SKIP testMode CUDA'
    echo 'SKIP testMode CUDA_Opt'
    echo 'SKIP testMode CUDA_Ref'
fi

if [ $testMode = "All" ] || [ $testMode = "DPC++" ]
then
    # run all tests with -testMode DPC++
    if [ "$HPCSCAN_CPP" = "mpiicpc -cxx=dpcpp" ]
    then
	sh testDriver.sh DPC++
    else
	echo 'SKIP testMode DPC++'
    fi
else
    echo 'SKIP testMode DPC++'
fi

if [ "$HPCSCAN_HIP" = "hipcc" ]
then
    if [ $testMode = "All" ] || [ $testMode = "HIP" ]
    then
	# run all tests with -testMode HIP
	sh testDriver.sh HIP
    else
	echo 'SKIP testMode HIP'
    fi

    if [ $testMode = "All" ] || [ $testMode = "HIP_Opt" ]
    then
	# run all tests with -testMode HIP_Opt
	sh testDriver.sh HIP_Opt
    else
	echo 'SKIP testMode HIP_Opt'
    fi
else
    echo 'SKIP testMode HIP'
    echo 'SKIP testMode HIP_Opt'
fi

if [ "$HPCSCAN_CPP" = "mpinc++" ]
then
    if [ $testMode = "All" ] || [ $testMode = "NEC" ]
    then
	# run all tests with -testMode NEC
	sh testDriver.sh NEC
    else
	echo 'SKIP testMode NEC'
    fi

    if [ $testMode = "All" ] || [ $testMode = "NEC_SCA" ]
    then
	# run all tests with -testMode NEC_SCA
	sh testDriver.sh NEC_SCA
    else
	echo 'SKIP testMode NEC_SCA'
    fi
else
    echo 'SKIP testMode NEC'
    echo 'SKIP testMode NEC_SCA'
fi

end_time=$(date)
echo "\n# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " $machine
echo "Done by    : " $tester

echo '* SUMMARY'
echo '# PASSED  :' $(grep PASSED ${report_file}  | wc -l)
echo '# FAILED  :' $(grep FAILED ${report_file}  | wc -l)
echo '# ERROR   :' $(grep 'E R R O R' ${report_file}  | wc -l)
echo '# WARNING :' $(grep 'W A R N I N G' ${report_file} | wc -l)
echo "==> Results of tests are saved in" ${report_file}

echo "\n========================================================================"
echo "                         END VALIDATION TESTS"
echo "========================================================================\n"
