# Overall makefile for mum
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL = /bin/sh

# Home directory of mum, please define environment variable MUM_HOME!
#MUM_HOME = ./../..

# Includes machine specific compilation options
include $(MUM_HOME)/Makefile.in

.PHONY: default help all clean_common clean_joe clean_ray clean_tools clean

default: help

help:
	@echo ' '
	@echo 'You must specify a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: '
	@echo 'all :             builds common files, joe, ray, charles and tools'
	@echo 'joe :             compiles our RANS code'
	@echo 'ray :             compiles our UQ manager'
	@echo 'common :          compiles all common files'
	@echo 'tools :           compile only the tools'
	@echo 'clean_joe :       clean all joe objects and executable files (included in common)'
	@echo 'clean_ray :       clean all ray objects and executable files (included in common)'
	@echo 'clean_tools :     clean all tools'
	@echo 'clean_common :    clean common files'
	@echo 'clean :           clean all objects and executable files (common and tools included)'
	@echo ' '
	@echo 'Do not forget to set up the environment variable MUM_HOME !'
	@echo ' '


# Each target is built by calling the respective makefiles within the src directories
all: common joe ray charles tools

common:
	$(MAKE) $@ -C $(MUM_HOME)/src/common

joe:
	$(MAKE) $@ -C $(MUM_HOME)/src/joe

ray:
	$(MAKE) $@ -C $(MUM_HOME)/src/ray

tools:
	$(MAKE) $@ -C $(MUM_HOME)/src/tools

clean_common:
	$(MAKE) clean -i -C $(MUM_HOME)/src/common

clean_joe:
	$(MAKE) clean -i -C $(MUM_HOME)/src/joe

clean_ray:
	$(MAKE) clean -i -C $(MUM_HOME)/src/ray	

clean_charles:
	@echo 'Not yet implemented'

clean_tools:
	$(MAKE) clean -i -C $(MUM_HOME)/src/tools

clean: clean_common clean_joe clean_ray clean_charles clean_tools
	rm $(MUM_HOME)/bin/*

