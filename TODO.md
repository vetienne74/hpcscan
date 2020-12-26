
# Current version

## Test mode OpenAcc

Finalize implementation

## Test mode Cuda

Finalize implementation

## propagator_Ac2.cpp

Optimal time step is not optimal when spacing grid sampling is not the same in the different axis

## testCase_Propa.cpp

Does not measure time spent in MPI communication (propaGpointFD = propaGpointEff, because testCase_time_com = 0)

# Additions for future versions

## Computation of spatial derivatives with FFT

## Propagator with pseudo spectral method (based on FFT)

## More complex propagator (such as elastic wave equation)
