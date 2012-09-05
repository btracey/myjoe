# Build makefile for joe
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep
include $(MUM_HOME)/src/joe/Makefile.dep

COMMON = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv.o Logging.o
COMMON_OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON))

# joe.o not included on purpose since user can have own joe.cpp
JOE_OBJS = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o
OBJS_AUX = $(COMMON_OBJS) $(JOE_OBJS) joe.o
JOE_LIBS = $(CLIBS) -lm

#Ronan
OBJS = $(OBJS_AUX) \
       $(MUM_HOME)/src/common/mixing_combustion/CombustionGeneralDefinitions.o\
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/ludcmp.o \
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/qrdcmp.o \
       $(MUM_HOME)/src/common/mixing_combustion/ChemtableCartesianLinear_single.o \
       $(MUM_HOME)/src/common/mixing_combustion/Table4D.o \
       $(MUM_HOME)/src/common/mixing_combustion/PressureTable.o

.PHONY: default help clean

default: joe

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'joe :      compiles joe (RANS solver)'
	@echo 'clean :    clean joe (included common files!)'
	@echo ' '

joe: obj joe.o
	$(CC) -o $@ $(CLD) $(OBJS) $(JOE_LIBS); mv $@ $(MUM_HOME)/bin

# joe.o treated separately because can be changed by user in case directory
joe.o: joe.cpp $(joe_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

# Other objects are both joe files and common files
obj: common_obj $(JOE_OBJS)

# Call "common makefile" recursively for all common files (update COMMON variable above if list has changed!)
common_obj:
	$(MAKE) $(COMMON) common_comb -C $(MUM_HOME)/src/common 

UgpWithCvCompFlow.o: UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $<

Scalars.o: Scalars.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $<

JoeWithModels.o: JoeWithModels.cpp $(JoeWithModels_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $<


clean:
	rm $(OBJS) $(MUM_HOME)/bin/joe 
