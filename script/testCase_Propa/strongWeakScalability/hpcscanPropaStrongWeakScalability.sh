
export KMP_AFFINITY=scatter,1,0,granularity=fine 
export OMP_NUM_THREADS=32 

start_time=$(date)

sh ../../clean_dir.sh

# Number of tries
ntry=10

# Strong Scalability
#-------------------
fdOrder=4
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2
fdOrder=8
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2

fdOrder=12
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2

# Weak Scalability
#-----------------
fdOrder=4
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 3000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 4000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2

fdOrder=8
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 3000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 4000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2

fdOrder=12
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 1000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 1 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 2 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 1000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10  \
			-nsub1 1 -nsub2 2 -nsub3 1

$HPCSCAN_MPI_INVOKER -n 4 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 2000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 2 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 6 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 3000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 3 -nsub3 2

$HPCSCAN_MPI_INVOKER -n 8 ../../../bin/hpcscan \
			-ntry ${ntry} -testMode CacheBlk -testCase Propa -propagator Ac2Standard \
			-dim 3 -nt 100 -dt 0.0 -ratioCFL 1.0 \
			-fdOrder ${fdOrder} -n1 1000 -n2 4000 -n3 2000 \
			-param1 130 -param2 130 -param3 130 -snapInc 10 \
			-nsub1 1 -nsub2 4 -nsub3 2

end_time=$(date)

echo '*** TEST ENDED ***'
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time
echo "On machine : " `hostname`
echo "Done by    : " `whoami`

