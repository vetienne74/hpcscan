
# Current version

## Test mode OpenAcc

Finalize implementation

## Test mode Cuda

Finalize implementation

## propagator_Ac2.cpp

Optimal time step is not optimal when spacing grid sampling is not the same in the different axis

## testCase_Propa.cpp

Elapse time includes MPI comm. but it does not allow to know the time spent in MPI communication (propaGpointFD = propaGpointEff, because testCase_time_com = 0)

## testCase_Grid.cpp

ApplyBoundaryConditionGB is not exact but a rough estimate

Add an option to select the type of grid points to consider: inner or all points

For accelerators: add function to measure transfer between device to CPU and vice versa

# Additions for future versions

## FD computations

Computation of spatial derivatives with FFT

## Propagator

Propagator with pseudo spectral method (based on FFT)

More complex propagator (such as elastic wave equation)

Implement higher order in time than 2nd order
