#include "CreatePressureChemTable.h"

void GetDatabase(string filename,int ifile) {

  ChemtableNormalized myChemTable;
  myChemTable.Load(filename);
  int n1,n2,n3,nvar;
  n1 = myChemTable.GetDimension1();
  n2 = myChemTable.GetDimension2();
  n3 = myChemTable.GetDimension3();
  nvar = myChemTable.GetNumberVariables();

  if (FinalTable.IsAllocated() == false)
  {
    // Allocation
    cout << "Allocate 5D Table" << endl;
    FinalTable.Allocate(nPress,n1,n2,n3,nvar);
    FinalTable.SetTableFileName(ChemTableFilename);
    nZm = FinalTable.GetDimension2();
    nZv = FinalTable.GetDimension3();
    nC  = FinalTable.GetDimension4();
    nVar= FinalTable.GetNumberVariables();
    for (int i=0; i<nVar;   ++i) FinalTable.AddVarName( myChemTable.GetVarName(i) );
    for (int i=0; i<nPress; ++i) FinalTable.SetCoordinate1(i, myChemTable.GetReferencePressure() );
    for (int i=0; i<nZm;    ++i) FinalTable.SetCoordinate2(i, myChemTable.GetCoordinate1(i) );
    for (int i=0; i<nZv;    ++i) FinalTable.SetCoordinate3(i, myChemTable.GetCoordinate2(i) );
    for (int i=0; i<nC;     ++i) FinalTable.SetCoordinate4(i, myChemTable.GetCoordinate3(i) );

  } else {
    // Check if the dimensions are consistent
    if ((n1 != nZm)||(n2 != nZv)||(n3 != nC)||(nvar != nVar))
    {
      cout << "Inconsistency in databases dimensions" <<endl;
      throw(-1);
    }
  }

  for (int i=0; i<nZm; i++)
    for (int j=0; j<nZv; j++)
      for (int k=0; k<nC; k++)
        for (int l=0; l<nVar; l++)
        {  
          FinalTable.SetTableValue(ifile,i,j,k,l,myChemTable.GetTableValue(i, j, k, l));
        }

  myChemTable.Unload();

}

int main(int argc, char * argv[]) {

  char * opt1 = argv[1];
  char * listfile;

  if (argc != 3) {
    cout << "Wrong number of arguments"<<endl;
    cout << "CreatePressureTable.e -dirlist <file>" << endl;
    throw(-1);
  }
  if (strcmp(opt1,"-dirlist") == 0)  {listfile = argv[2];}
  else {
    cout<<"Wrong option : "<<opt1<<endl;
    cout << "CreatePressureTable.e -dirlist <file>" << endl;
    throw(-1);
  }

  cout << "Opening dirlist : " << listfile << endl;
  ifstream fh(listfile);
  int rank;
  int nfiles;
  string dir;
  double myPress;
  vector<string> databaseNames;
  vector<double> PressVec;
  while (fh.eof() == 0)
  {
    if (fh  >> rank >> dir >> myPress)
    {
      cout << " " << rank << " " << dir <<  endl;
      databaseNames.push_back(dir+"/database");
      PressVec.push_back(myPress);
    }

  }
  fh.close();

  nPress = databaseNames.size();
  // Sort PressVec in Pressure
  Pressure = new double[nPress];
  vector<string> databaseNames_cpy = databaseNames;

  for (int i=0; i<nPress; ++i) {
		int iloc = MinIndex(&PressVec[0],nPress,0);
		Pressure[i] = PressVec[iloc];
		databaseNames[i] = databaseNames_cpy[iloc];
		PressVec[iloc] = 1.0e40;
  }

  for (int ifile = 0 ; ifile < nPress; ++ifile)
  {
    cout << databaseNames[ifile] << " " << Pressure[ifile] << endl;
    GetDatabase(databaseNames[ifile],ifile);
  }

  FinalTable.WriteTable("debug");
  FinalTable.Unload();

  return 0;
}
