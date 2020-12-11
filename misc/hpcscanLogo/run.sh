
sh ../../script/clean_dir.sh

mpirun -n 1 ../../bin/hpcscan -testCase Propa -writeGrid \
       -tmax 0.2 -snapDt 0.1 \
       -dim 2 -n1 200 -n2 600 \
       -param1 4 -param2 8
