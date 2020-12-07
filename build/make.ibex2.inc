
#LIBS = -L/sw/csg/openmpi/3.0.0/el7.5_pgi17.10/executables/lib/libmpi.a

CUDAFLAGS = -L/usr/local/cuda/lib64 -lcuda -lcudart

# compiler options for INTEL
COPTION = -w -g -O3 -std=c++11 -acc -ta=tesla:cc70 -mp -Minfo=accel 

# INTEL
$(info compile with Intel)
CPP     = mpic++
CC      = mpicc
FC      = mpif90
NVCC    = nvcc
