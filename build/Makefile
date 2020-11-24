#----------------------------------------------------------------#
#                                                                #  
# DO NOT FORGET TO CREATE THE INCLUDE FILE "make.MY_HOST.inc"    #  
# ACCORDING TO THE CONFIGURATION OF YOUR COMPUTING MACHINE       # 
# AND THEN INCLUDE THIS FILE IN THE MAKEFILE AS IT IS DONE BELOW #
#                                                                #
#----------------------------------------------------------------# 


#======================#
#     include file     #
#======================#

HOST=$(shell hostname)

# matches aurora (NEC)
ifeq ($(findstring aurora,$(HOST)), aurora)
include make.nec.inc
endif
ifeq ($(findstring robin,$(HOST)), robin)
include make.nec.inc
endif

# matches cn (Ibex)
ifeq ($(findstring cn,$(HOST)), cn)
include make.ibex.inc
endif

# matches cdl (Shaheen)                                                                                       
ifeq ($(findstring cdl,$(HOST)), cdl)
include make.shaheen.inc
endif

# matches gateway (Shaheen)
ifeq ($(findstring gateway,$(HOST)), gateway)
include make.shaheen.inc
endif

# matches jupiter
ifeq ($(findstring jupiter,$(HOST)), jupiter)
include make.jupiter.inc
GIT_VER=1
endif

# matches mars
ifeq ($(findstring mars,$(HOST)), mars)
include make.jupiter.inc
GIT_VER=1
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
COPTION+=-D_DOUBLE_PRECISION_
BUILD_MESSAGE='*** BUILD OK (DOUBLE PRECISION) *** '
else
BUILD_MESSAGE='*** BUILD OK (SINGLE PRECISION) *** '
endif

CFLAGS = -I. -I${INC_DIR} $(COPTION)

#=====================#
#     objects         #
#=====================#

# all objects
ALL_OBJ = config.o grid.o grid_CacheBlk.o grid_Factory.o grid_NEC_SCA.o \
	main.o output_report.o propagator_Factory.o propagator_Ac2.o propagator_Ac2SplitComp.o \
	testCase.o testCase_Comm.o testCase_FD_D2.o testCase_Grid.o \
	testCase_Propa.o testCase_Memory.o testCase_Template.o testCase_Util.o

.SUFFIXES: .o .cpp
.cpp.o:
	$(CPP) -c $(CFLAGS) $<

all: version hpcscan

hpcscan: $(ALL_OBJ)
	@mkdir -p $(BIN_DIR)
	$(CPP) $(CFLAGS) -o $(BIN_DIR)/$@ $(ALL_OBJ) $(LIBS)
	@echo $(BUILD_MESSAGE)
	@$(CPP) --version

clean:
	rm -f *.L *.o *.mod *~ ./version $(BIN_DIR)/hpcscan

version:
	@if [ ${GIT_VER} ] ; then \
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









