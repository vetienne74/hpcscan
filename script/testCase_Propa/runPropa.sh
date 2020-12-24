
$HPCSCAN_CPP --version >& perf.out
mpirun -n 1 ../../bin/hpcscan -v >> perf.out

# 500^3 Baseline Ac2Standard
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode Baseline -autoPad -nt 100 -propagator Ac2Standard'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode Baseline -autoPad -nt 100 -propagator Ac2Standard -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 500^3 Baseline Ac2SplitComp
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode Baseline -autoPad -nt 100 -propagator Ac2SplitComp'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode Baseline -autoPad -nt 100 -propagator Ac2SplitComp -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 500^3 NEC Ac2Standard
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC -autoPad -nt 100 -propagator Ac2Standard'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC -autoPad -nt 100 -propagator Ac2Standard -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 500^3 NEC Ac2SplitComp
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC -autoPad -nt 100 -propagator Ac2SplitComp'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC -autoPad -nt 100 -propagator Ac2SplitComp -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 500^3 NEC_SCA Ac2SplitComp
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC_SCA -autoPad -nt 100 -propagator Ac2SplitComp'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 500 -n2 500 -n3 500 -testMode NEC_SCA -autoPad -nt 100 -propagator Ac2SplitComp -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 1000^3 NEC Ac2Standard
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 1000 -n2 1000 -n3 1000 -testMode NEC -autoPad -nt 100 -propagator Ac2Standard'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 1000 -n2 1000 -n3 1000 -testMode NEC -autoPad -nt 100 -propagator Ac2Standard -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

# 1000^3 NEC_SCA Ac2SplitComp
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 1000 -n2 1000 -n3 1000 -testMode NEC_SCA -autoPad -nt 100 -propagator Ac2SplitComp'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out
hpcscanCmd='mpirun -n 1 ../../bin/hpcscan -testCase Propa -n1 1000 -n2 1000 -n3 1000 -testMode NEC_SCA -autoPad -nt 100 -propagator Ac2SplitComp -boundary None'
echo $hpcscanCmd | tee >> perf.out
$hpcscanCmd > tmp 
grep 'errTestCase2' tmp | tee >> perf.out
grep 'Gpoint/s eff' tmp | tee >> perf.out

cat perf.out
