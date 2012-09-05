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

COMMON_JOE = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv2.o Logging.o
JOE_OBJS = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o
OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) $(addprefix $(MUM_HOME)/src/joe/, $(JOE_OBJS)) joe.o
JOE_LIBS = $(CLIBS) -lm

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
RAY_OBJS = 
OBJS_RAY = $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) $(addprefix $(MUM_HOME)/src/joe/, $(JOE_OBJS)) \
	   $(addprefix $(MUM_HOME)/src/ray/, $(RAY_OBJS)) ray.o
RAY_LIBS = $(CLIBS) -lm

ray: obj_joe obj_ray ray.o
	$(CC) -o $@ $(CLD) $(OBJS_RAY) $(JOE_LIBS) $(RAY_LIBS)

obj_ray:
	$(MAKE) obj -C $(MUM_HOME)/src/ray

ray.o: ray.cpp $(ray_H)
	$(CC) $(CFLAGS) -I $(MUM_HOME)/src/common -I $(MUM_HOME)/src/joe -I $(MUM_HOME)/src/ray -c $<

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
