#include "ChemLinToTecplot.h"

int main(int argc, char * argv[]) {

  char *tablename = argv[1];
  Table.Load(tablename);
  Table.WriteTableTecplot();

  int nZ   = Table.GetDimension1();
  int nZv  = Table.GetDimension2();
  int nC   = Table.GetDimension3();
  int nVar = Table.GetNumberVariables();

  double *Z  = new double[nZ];
  double *Zv = new double[nZv];
  double *C  = new double[nC];

  Table.CopyCoordinate1(Z);
  Table.CopyCoordinate2(Zv);
  Table.CopyCoordinate3(C);

  int n_E0 = Table.GetVarNameIndex("E0"); 
  int n_RoM = Table.GetVarNameIndex("RoM");
  int n_rho = Table.GetVarNameIndex("rho0");

  ofstream fout;
  fout.open("E0.dat");
  fout << "VARIABLES=\"ZMean\" \"E0\" \"RoM\" \"rho\" "<< endl;
  for (int j=0; j<nZv; j++) {
    fout << "ZONE T=Sz_"<< Zv[j] << ", I=" << nZ << endl;
    for (int i=0; i<nZ; i++) {
      fout << Z[i] << " " 
           << Table.GetTableValue(i,j,0,n_E0) << " " 
           << Table.GetTableValue(i,j,0,n_RoM) << " "
           << Table.GetTableValue(i,j,0,n_rho) << " "
           << endl;
    }
  }
  fout.close();

  delete [] Z;
  delete [] Zv;
  delete [] C;

  return 0;
}

