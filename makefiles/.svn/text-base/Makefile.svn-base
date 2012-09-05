# Case makefile for joe, ray and charles
# Author: Vincent Terrapon
# Version: 1.1
# Date: 12/2009

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep
include $(MUM_HOME)/src/joe/Makefile.dep
include $(MUM_HOME)/src/ray/Makefile.dep
#include $(MUM_HOME)/src/charles/Makefile.dep


.PHONY: default help clean

default: help

help:
	@echo ' '
	@echo 'You must specify a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: '
	@echo 'joe :        compiles our RANS code'
	@echo 'ray :        compiles our UQ manager'
	@echo 'charles :    compiles our LES code'
	@echo 'clean :      clean only object and executable files in case directory'
	@echo ' '
	@echo 'Do not forget to set up the environment variable MUM_HOME !'
	@echo ' '

##########
# JOE
##########

COMMON_JOE = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv.o Logging.o
JOE_OBJS = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o
OBJS_AUX = $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) $(addprefix $(MUM_HOME)/src/joe/, $(JOE_OBJS)) joe.o
JOE_LIBS = $(CLIBS) -lm

#Ronan
OBJS = $(OBJS_AUX) \
       $(MUM_HOME)/src/common/mixing_combustion/CombustionGeneralDefinitions.o\
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/ludcmp.o \
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/qrdcmp.o \
       $(MUM_HOME)/src/common/mixing_combustion/ChemtableCartesianLinear_single.o \
       $(MUM_HOME)/src/common/mixing_combustion/Table4D.o \
       $(MUM_HOME)/src/common/mixing_combustion/PressureTable.o \
       $(MUM_HOME)/src/common/mixing_combustion/Table3D.o \
       $(MUM_HOME)/src/common/mixing_combustion/ChemtableNormalized.o

joe: obj_joe joe.o
	$(CC) -o $@ $(CLD) $(OBJS) $(JOE_LIBS)

obj_joe:
	$(MAKE) obj -C $(MUM_HOME)/src/joe

joe.o: joe.cpp $(joe_H)
	$(CC) $(CFLAGS) -I $(MUM_HOME)/src/common -I $(MUM_HOME)/src/common/mixing_combustion -I $(MUM_HOME)/src/joe -c $<

joe.cpp:
	cp $(MUM_HOME)/src/joe/joe.cpp .


##########
# RAY
##########
COMMON_RAY = $(COMMON_JOE)
RAY_OBJS = 
OBJS_RAY = $(addprefix $(MUM_HOME)/src/common/, $(COMMON_RAY)) $(addprefix $(MUM_HOME)/src/ray/, $(RAY_OBJS)) ray.o
RAY_LIBS = $(CLIBS) -lm

ray: obj_ray ray.o
	$(CC) -o $@ $(CLD) $(OBJS_RAY) $(JOE_LIBS) $(RAY_LIBS)

obj_ray:
	$(MAKE) obj -C $(MUM_HOME)/src/ray

ray.o: ray.cpp $(ray_H)
	$(CC) $(CFLAGS) -I $(MUM_HOME)/src/common -I $(MUM_HOME)/src/ray -c $<

ray.cpp: 
	cp $(MUM_HOME)/src/ray/ray.cpp .


##########
# CHARLES
##########

charles:
	@echo 'Not yet implemented!'


##########
# ALL
##########

clean:
	rm joe ray charles *.o
