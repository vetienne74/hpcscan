#
# libraries
#
INC = -openmp

# MKL
LIB_MKL= -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm

#
# compiler options for INTEL
# 
#COPTION = -C -g -openmp
#COPTION = -O3 -g -openmp -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -xavx -DMKL_ILP64 -qopenmp -I${MKLROOT}/include
#COPTION = -O3 -g -openmp -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -xavx
#COPTION = -O3 -g -openmp -no-inline-max-size -no-inline-max-per-routine -no-inline-max-per-compile
#-D NDEBUG is used to avoid conflict with assert and nr3.h
#COPTION = -O3 -g -openmp -DNDEBUG
COPTION = -O3 -g -openmp

#
# compilers
#
CPP     = CC
CC      = cc 
FC      = ftn

