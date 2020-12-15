
# build
cd ../../build
make clean
make | tee ../misc/fileForReadme/make.txt
cd -

# command line parameters
../../bin/hpcscan -h | tee commandLineParam.txt

# hpcscan version
../../bin/hpcscan -v |tee version.txt

# run Propa test case
cd ../../script/
mpirun -n 1 ../bin/hpcscan -testCase Propa -testMode CacheBlk | tee ../misc/fileForReadme/runPropaTestCase.txt

# run validation tests
sh runValidationTests.sh | tee ../misc/fileForReadme/runValidationTests.txt
cd -
