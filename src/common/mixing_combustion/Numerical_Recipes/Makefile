#########################################################
#                                                       #
# Makefile for Numerical Recipes                        #
#                                                       #
#########################################################

include $(MUM_HOME)/Makefile.in

.SUFFIXES: .o .cpp

.PHONY: default clean

default: all

all : qrdcmp.o ludcmp.o 

clean :
	rm -f *.o ;
                                                                                                                                                
.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp ;
