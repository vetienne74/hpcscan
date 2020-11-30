
LIBS = -fopenmp

CUDAFLAGS = -L/usr/local/cuda/lib64 -lcuda -lcudart

# compiler options for INTEL
COPTION = -w -g -O3 -fopenmp -std=c++11 -xHost -ffast-math -fpermissive

# INTEL
$(info compile with Intel)
CPP     = mpiicpc
CC      = mpiicc
FC      = mpiifort
NVCC    = nvcc
