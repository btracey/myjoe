#ifndef CREATEPRESSURETABLE_H
#define CREATEPRESSURETABLE_H

//#include "ChemtableCartesianLinear_single.h"
#include "ChemtableNormalized.h"
#include "PressureTable.h"

#include <fstream>
using std::ifstream;

int       nZm,nZv,nC,nPress,nVar;
string    ChemTableFilename = "multiP_test";
double    *Pressure;
PressureTable FinalTable;

void GetDatabase(string filename,int ifile);

#endif
