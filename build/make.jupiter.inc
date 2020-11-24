
LIBS = -fopenmp 

#
# compiler options for INTEL
# 
#COPTION = -C -g -openmp
#COPTION = -O3 -g -openmp -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -xavx -DMKL_ILP64 -qopenmp -I${MKLROOT}/include
#COPTION = -O3 -g -openmp -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -xavx

#
# compiler options for g++
# 
#COPTION = -w -g -O3 -fopenmp -ftree-vectorize -fpermissive -mavx2
#COPTION = -w -g -O3 -fopenmp -ffast-math -fpermissive -mavx2
# -ffast-math prevents detection of NaN values
COPTION = -w -g -O3 -fopenmp -mavx2
#COPTION = -w -g -O0 -finstrument-functions -fopenmp -fpermissive

ifdef ONEAPI_ROOT
# INTEL
$(info compile with Intel ONEAPI)
CPP     = mpiicpc
CC      = mpiicc
FC      = mpiifort
#CPP     = mpicxx
#CC      = mpicxx
#FC      = mpiifort
else
# GNU
$(info compile with GNU compilers)
CPP     = mpic++
CC      = mpicc
FC      = mpif90
endif
