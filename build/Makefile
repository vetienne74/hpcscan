
#========================#
# check env has been set #
#========================#
ifeq ($(strip $(HPCSCAN_CPP)),)
$(info It seems hpcscan environment has not been set.)
$(info Source one of the env files in ../env or create a new one for your machine.)
$(error HPCSCAN_CPP is $(HPCSCAN_CPP))
else
$(info hpcscan will be compiled with $(HPCSCAN_CPP))
endif

#=====================#
#     directories     #
#=====================#

# source files
SRC_DIR        = ../src

# header files
INC_DIR        = ../src

# executable files
BIN_DIR        = ../bin

#=====================#
#       vpaths        #
#=====================#

VPATH = $(SRC_DIR)

#=====================#
#       options       #
#=====================#

# SINGLE / DOUBLE PRECISION
ifeq ($(findstring double,$(precision)), double)
HPCSCAN_CPP_FLAGCOMP+=-D_DOUBLE_PRECISION_
BUILD_MESSAGE='*** BUILD OK (DOUBLE PRECISION) *** '
else
BUILD_MESSAGE='*** BUILD OK (SINGLE PRECISION) *** '
endif

# WITH CUDA
ifneq ($(strip $(HPCSCAN_CUDA)),)
HPCSCAN_CPP_FLAGCOMP+=-D__CUDA__
endif

# WITH OPENACC
ifneq ($(strip $(HPCSCAN_CPP_OPENACC_FLAG)),)
HPCSCAN_CPP_FLAGCOMP+=-D__OPENACC__
HPCSCAN_CPP_FLAGCOMP+=$(HPCSCAN_CPP_OPENACC_FLAG)
endif

# WITH NEC
ifeq ($(HPCSCAN_CPP), mpinc++)
HPCSCAN_CPP_FLAGCOMP+=-D__NEC__
endif

$(info HPCSCAN_CPP_FLAGCOMP="$(HPCSCAN_CPP_FLAGCOMP)")

#=====================#
#     objects         #
#=====================#

# cpp standard objects
CPP_OBJ = config.o grid.o grid_CacheBlk.o grid_Factory.o \
	main.o output_report.o propagator_Factory.o propagator_Ac2.o propagator_Ac2SplitComp.o \
	testCase.o testCase_Comm.o testCase_FD_D2.o testCase_Grid.o \
	testCase_Propa.o testCase_Memory.o testCase_Template.o testCase_Util.o

# cpp OpenACC objects
CPP_OPENACC_OBJ = grid_OpenAcc.o

# Cuda objects
CUDA_OBJ = grid_Cuda.o

# NEC objects
CPP_NEC_OBJ = grid_NEC.o grid_NEC_SCA.o

# build list of needed objects
#-----------------------------
# common C++ objects
ALL_OBJ = ${CPP_OBJ}
ifneq ($(strip $(HPCSCAN_CUDA)),)
# CUDA objects
ALL_OBJ += ${CUDA_OBJ}
endif
# OpenACC objects
ifneq ($(strip $(HPCSCAN_CPP_OPENACC_FLAG)),)
ALL_OBJ += ${CPP_OPENACC_OBJ}
endif
# NEC objects
ifeq ($(HPCSCAN_CPP), mpinc++)
ALL_OBJ += ${CPP_NEC_OBJ}
endif

.SUFFIXES: .o .cpp .cu

.cu.o:
	$(HPCSCAN_CUDA) -c $(HPCSCAN_CUDA_FLAGCOMP) $<
.cpp.o:
	$(HPCSCAN_CPP) -c $(HPCSCAN_CPP_FLAGCOMP) $<

all: version hpcscan

hpcscan: $(ALL_OBJ)
	@mkdir -p $(BIN_DIR)
	$(HPCSCAN_CPP) $(HPCSCAN_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ $(ALL_OBJ) $(HPCSCAN_CPP_LIB) $(HPCSCAN_CUDA_LIB)
	@echo $(BUILD_MESSAGE)
	@$(HPCSCAN_CPP) --version

clean:
	rm -f *.L *.o *.mod *~ ./version $(BIN_DIR)/hpcscan

version:
	@if [ ${HPCSCAN_GIT_VER} ] ; then \
		./git_info.sh ../ > $@ ;\
		if test ! -e version ; then \
			echo "No version found, initialization" ;\
			cp version $(SRC_DIR)/version_hpcscan.h ; \
		fi;\
		if ! diff -q version $(SRC_DIR)/version_hpcscan.h ; then \
			echo "version has changed, regenerating version info"; \
			cp version $(SRC_DIR)/version_hpcscan.h ; \
		else \
			echo "version info is up to date"; \
		fi;\
	fi









