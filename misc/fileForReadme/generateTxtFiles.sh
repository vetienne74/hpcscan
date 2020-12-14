
../../bin/hpcscan -h | tee commandLineParam.txt

../../bin/hpcscan -v |tee version.txt

cd ../../script/
sh runValidationTests.sh | tee ../misc/fileForReadme/runValidationTests.txt
cd -
