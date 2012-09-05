include $(MUM_HOME)/Makefile.in

############################################################################
# REGLES DE COMPILATION
############################################################################

.SUFFIXES: .o .cpp 
.PHONY: all clean default

############################################################################
# Cibles
############################################################################
default: all

all : common_num CombustionGeneralDefinitions.o ChemtableCartesianLinear_single.o Table4D.o PressureTable.o Table3D.o ChemtableNormalized.o pdf_functions.o 

common_num:
	$(MAKE) all -C Numerical_Recipes 

# Dependencies
CombustionGeneralDefinitions.o : CombustionGeneralDefinitions.h
Table4D.o : CombustionGeneralDefinitions.h Table4D.h
Table3D.o : CombustionGeneralDefinitions.h Table3D.h
PressureTable.o : CombustionGeneralDefinitions.h Table4D.h PressureTable.h
ChemtableNormalized.o : CombustionGeneralDefinitions.h Table3D.h ChemtableNormalized.h 
ChemtableCartesianLinear_single.o : ChemtableCartesianLinear_single.cpp
pdf_functions.o : pdf_functions.cpp pdf_functions.h

clean :
	rm -f *.o *.mod *.e ;
	$(MAKE) clean -C Numerical_Recipes

.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp -I../ ;
