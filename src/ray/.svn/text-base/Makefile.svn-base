# Build makefile for ray
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

COMMON = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv.o Logging.o
COMMON_OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON))

# ray.o not included on purpose since user can have own ray.cpp
RAY_OBJS = 

OBJS = $(COMMON_OBJS) $(RAY_OBJS) ray.o
RAY_LIBS = $(CLIBS) -lm

.PHONY: default help clean

default: ray

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'ray :      compiles ray (UQ manage)'
	@echo 'clean :    clean ray (included common/joe files!)'
	@echo ' '

ray: obj ray.o
	$(CC) -o $@ $(CLD) $(OBJS) $(RAY_LIBS); mv $@ $(MUM_HOME)/bin

# ray.o treated separately because can be changed by user in case directory
ray.o: ray.cpp $(ray_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ray -c $<

# Other objects are both joe files and common files
obj: common_obj $(RAY_OBJS)

# Call "common makefile" recursively for all common files (update COMMON variable above if list has changed!)
common_obj:
	$(MAKE) $(COMMON) -C $(MUM_HOME)/src/common

clean:
	rm $(OBJS) $(MUM_HOME)/bin/ray
