# Build makefile for tools
# Author: Vincent Terrapon
# Version: 1.0
# Date: 08/2009

# NOTICE: if you add tools, don't forget to 
#         * update the variables TOOLS with the name of your executable
#         * update the variable COMMON with the common object files which are required and not yet in the list
#         * update the Makefile.dep with the dependencies for your tool
#         * add the compilation rules for your tool

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h> and object files <file.o>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep
include $(MUM_HOME)/src/joe/Makefile.dep
include $(MUM_HOME)/src/tools/Makefile.dep

TOOLS = TestCombustion TestFlamelet CreateChemTable CreateChemTable_new mum_hold ThermoProperties CreatePressureChemTable TestPressureTable ChemLinToTecplot prepro
TOOLS_EXECS = $(addprefix $(MUM_HOME)/bin/, $(TOOLS))
TOOLS_OBJS = $(addsuffix .o, $(TOOLS))
TOOL_LIBS = $(CLIBS) -lm

COMMON = Param.o MpiStuff.o MiscUtils.o mpi_tommie.o Gp.o Ugp.o UgpWithTools.o \
	CdpFilter.o MshFilter.o
COMMON_OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON))

.PHONY: default help clean

default: tools

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'tools :      compile all tools'
	@echo '<toolname>:  compile specific tool <toolname>'
	@echo 'clean :      clean tools (executable and object files)'
	@echo ' '
	@echo 'The available tools are: $(TOOLS)'
	@echo ' '

tools: $(TOOLS)

##################################################################################################################
### CreateChemtable4d ###

CreateChemtable4d: CreateChemtable4d.o $(CreateChemtable4d_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

CreateChemtable4d.o: CreateChemtable4d.cpp $(Combustion_H) $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### TestCombustion ###

TestCombustion: TestCombustion.o $(TestCombustion_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

TestCombustion.o: TestCombustion.cpp $(Combustion_H) $(Flamelet_H) $(ChemistryTable_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### TestFlamelet ###

TestFlamelet: TestFlamelet.o $(TestFlamelet_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

TestFlamelet.o: TestFlamelet.cpp $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### CreateChemTable ###

CreateChemTable: CreateChemTable.o $(CreateChemTable_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

CreateChemTable.o: CreateChemTable.cpp CreateChemTable.h $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### CreateChemTable_new ###

CreateChemTable_new: CreateChemTable_new.o $(CreateChemTable_new_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

CreateChemTable_new.o: CreateChemTable_new.cpp CreateChemTable_new.h $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### CreatePressureChemTable ###

CreatePressureChemTable: CreatePressureChemTable.o $(CreatePressureChemTable_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

CreatePressureChemTable.o: CreatePressureChemTable.cpp CreatePressureChemTable.h $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### ChemLinToTecplot ###

ChemLinToTecplot: ChemLinToTecplot.o $(ChemLinToTecplot_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

ChemLinToTecplot.o: ChemLinToTecplot.cpp ChemLinToTecplot.h
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<
##################################################################################################################
### TestPressureTable ###

TestPressureTable: TestPressureTable.o $(TestPressureTable_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

TestPressureTable.o: TestPressureTable.cpp TestPressureTable.h $(Flamelet_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### mum_hold ###

mum_hold: mum_hold.o $(mum_hold_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

mum_hold.o: mum_hold.cpp $(MpiStuff_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $<

##################################################################################################################
### ThermoProperties ###

ThermoProperties: ThermoProperties.o $(ThermoProperties_DEP)
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

ThermoProperties.o: ThermoProperties.cpp $(Combustion_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

##################################################################################################################
### 2dto3d ###

2dto3d: 2dto3d.cpp
	$(CC) -o $@ $(CLD) $^ $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

##################################################################################################################
### preprocessing ###

prepro: prepro.o $(Prepro_DEP)
	$(CC) -o $@ $(CLD) $(Prepro_DEP) prepro.o $(TOOL_LIBS); mv $@ $(MUM_HOME)/bin

prepro.o: prepro.cpp $(Prepro_H)
	$(CC) $(CFLAGS) -c prepro.cpp -I $(MUM_HOME)/src/common

##################################################################################################################


# Call "common makefile" recursively for all common files (update COMMON variable above if list has changed!)
$(COMMON_OBJS):
	$(MAKE) $(COMMON) -C $(MUM_HOME)/src/common 


clean:
	rm $(TOOLS_OBJS) $(TOOLS_EXECS) 
