
#------------------------------------------------------------------------------------
# this script is used to run the test cases for performance measurements
# the test cases are quite heavy and require HPC computing platforms
#------------------------------------------------------------------------------------

sh clean_dir.sh

rep='.'

echo "========================================================================"
echo "                         START TEST CASES"
echo "========================================================================"
start_time=$(date)

#####################################################################################
testCase1='testCase_Comm'
cd ${rep}/${testCase1}
script1='runTestShaheen'
mkdir -p results
sh ${script1}.sh | tee ./results/${script1}.out
cd -

#####################################################################################
testCase2='testCase_FD_D2'
cd ${rep}/${testCase2}
script2='runSmallGridShaheen'
mkdir -p results
sh ${script2}.sh | tee ./results/${script2}.out
script2b='runMediumGridShaheen'
sh ${script2b}.sh | tee ./results/${script2b}.out
script2c='runStrongScalaShaheen'
sh ${script2c}.sh | tee ./results/${script2c}.out
script2d='runWeakScalaShaheen'
sh ${script2d}.sh | tee ./results/${script2d}.out
cd -

#####################################################################################
testCase3='testCase_Grid'
cd ${rep}/${testCase3}
script3='runSmallGridShaheen'
mkdir -p results
sh ${script3}.sh | tee ./results/${script3}.out
script3b='runMediumGridShaheen'
sh ${script3b}.sh | tee ./results/${script3b}.out
cd -

#####################################################################################
testCase4='testCase_Memory'
script4='runTestShaheen'
cd ${rep}/${testCase4}
mkdir -p results
sh ${script4}.sh | tee ./results/${script4}.out
cd -

end_time=$(date)
echo "========================================================================"
echo "                         END TEST CASES"
echo "========================================================================"

echo '------------------------------------------------------------------------'
echo 'SUMMARY:'
echo ${testCase1}/${script1}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase1}/results/${script1}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase1}/results/${script1}.out | wc -l)
echo ${testCase2}/${script2}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase2}/results/${script2}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase2}/results/${script2}.out | wc -l)
echo ${testCase2}/${script2b}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase2}/results/${script2b}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase2}/results/${script2b}.out | wc -l)
echo ${testCase2}/${script2c}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase2}/results/${script2c}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase2}/results/${script2c}.out | wc -l)
echo ${testCase2}/${script2d}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase2}/results/${script2d}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase2}/results/${script2d}.out | wc -l)
echo ${testCase3}/${script3}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase3}/results/${script3}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase3}/results/${script3}.out | wc -l)
echo ${testCase3}/${script3b}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase3}/results/${script3b}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase3}/results/${script3b}.out | wc -l)
echo ${testCase4}/${script4}
echo '# PASSED  :' $(grep PASSED ${rep}/${testCase4}/results/${script4}.out | wc -l)
echo '# FAILED  :' $(grep FAILED ${rep}/${testCase4}/results/${script4}.out | wc -l)

echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname` 
echo "Done by    : " `whoami`
echo '------------------------------------------------------------------------'

