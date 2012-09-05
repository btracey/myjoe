#!/usr/bin/env python
#####################################################
# read Flamemaster files and                        #
# write in tecplot format                           #
# 12/02/09                                          #
# Note: take .out as input and create .dat          #
#       input is a file created by LT tool          #
#       file .out is a function of chi_st at        #
#       specific Z                                  #
# Use: FM2Tec360_bis *.out                          #
#####################################################
import sys

def readData(filename):
    file = open(filename,'r')
    lines = file.readlines()
    file.close()
    variables = lines[1].split()
    cleanVariableNames(variables)
    data = lines[2:]
    return variables,data
    

def cleanVariableNames(list):
    for entry in list:
        if entry.startswith('[')and entry.endswith(']'):
            list[(list.index(entry)-1)] = list[(list.index(entry)-1)]+' '+entry
            list.pop(list.index(entry))

def writeTecplotFile(filename,variables,data):
    file = open(filename,'w')
    file.write('TITLE=\"'+filename+'\"\n')
    file.write('VARIABLES=')
    for entry in variables:
        file.write('\"'+entry+'\"\n')
    file.write('ZONE T=\"Zone 1\"\n')
    file.write('I='+str(len(data))+', J=1, K=1, \n')
    file.write('ZONETYPE=Ordered \n')
    file.write('DATAPACKING=POINT \n')
    for line in data:
        file.write(line+'\n')
    file.close()


def tecplotFilename(origName):
    return origName.replace('.out','.dat')

#def getChi(filename):
#    chi = filename.split('chi')[1].split('tf')[0]
#    return chi


if len(sys.argv)==1:
    print 'Give at least 1 file!'
else:
    commandLineArguments = sys.argv[1:]
    for argument in commandLineArguments:
        variables,data = readData(argument)
        writeTecplotFile(tecplotFilename(argument),variables,data)
        print 'Done!'
