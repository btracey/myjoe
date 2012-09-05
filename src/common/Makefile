# Build makefile for common files
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep

COMMON_OBJ = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o \
             MshFilter.o UgpWithCv.o tc_vec3d.o
#COMMON_SRC = $(patsubst %.o, %.cpp,$(COMMON_OBJ))

.PHONY: default common help clean

default: common

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'all :    compiles common files used by joe, ray and charles'
	@echo 'clean :  clean common files'
	@echo ' '

common: $(COMMON_OBJ) common_comb

common_comb:
	$(MAKE) all -C mixing_combustion 

# Object files <file.o> dependencies
# Dependencies of the include files <file.h> are explicitly given in the above included files Makefile.dep
mpi_tommie.o: mpi_tommie.cpp $(mpi_tommie_H)
	$(CC) $(CFLAGS) -c $<

MpiStuff.o: MpiStuff.cpp $(MpiStuff_H)
	$(CC) $(CFLAGS) -c $<

Gp.o: Gp.cpp $(Gp_H)
	$(CC) $(CFLAGS) -c $<

MiscUtils.o: MiscUtils.cpp $(MiscUtils_H)
	$(CC) $(CFLAGS) -c $<

Ugp.o: Ugp.cpp $(Ugp_H)
	$(CC) $(CFLAGS) -c $<

Param.o: Param.cpp $(Param_H)
	$(CC) $(CFLAGS) -c $<

UgpWithTools.o: UgpWithTools.cpp $(UgpWithTools_H)
	$(CC) $(CFLAGS) -c $<

CdpFilter.o: CdpFilter.cpp $(CdpFilter_H)
	$(CC) $(CFLAGS) -c $<

MshFilter.o: MshFilter.cpp $(MshFilter_H)
	$(CC) $(CFLAGS) -c $<

UgpWithCv.o: UgpWithCv.cpp $(UgpWithCv_H)
	$(CC) $(CFLAGS) -c $<

UgpWithCvFake.o: UgpWithCvFake.cpp $(UgpWithCvFake_H)
	$(CC) $(CFLAGS) -c $<

tc_vec3d.o: tc_vec3d.cpp $(tc_vec3d_H)

Logging.o: Logging.cpp $(CTR_TOMMIETOO_LOGGING_H)
	$(CC) $(CFLAGS) -c $<

# Cleaning
clean:
	rm *.o
	$(MAKE) clean -C mixing_combustion
