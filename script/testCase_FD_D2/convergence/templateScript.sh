
$HPCSCAN_MPI_INVOKER -n 1 ../../../bin/hpcscan \
			-ntry _ntry_ -testMode _testMode_ -testCase FD_D2 \
			-dim _dim_ \
			-fdOrder _fdOrder_ -n1 _n1_ -n2 _n2_ -n3 _n3_ \
			-param1 _nmode_ -param2 _nmode_ -param3 _nmode_ \
                        -autoPad -hwCounterDt 0.5
