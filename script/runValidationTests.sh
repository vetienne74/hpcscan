
export tester=`whoami`
export machine=`hostname`
export report_file=runValidationTests.${machine}.${tester}.out

sh clean_dir.sh
rm -f ${report_file}

echo "\n========================================================================"
echo "                         START VALIDATION TESTS"
echo "========================================================================\n"
start_time=$(date)

# run all tests with -testMode Baseline
sh testDriver.sh Baseline

# run all tests with -testMode CacheBlk
sh testDriver.sh CacheBlk

# run all tests with -testMode NEC_SCA
if [ -z ${__NEC__} ]
then
    echo 'skip testDriver.sh NEC_SCA'
else
    sh testDriver.sh NEC_SCA
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
