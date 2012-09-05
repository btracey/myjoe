#include "ChemtableNormalized.h"

//------------------------------//
//         Constructors         //
//------------------------------//
ChemtableNormalized::ChemtableNormalized()
{
  CombustionModel = "FPVA";
  ChemTableType = "COEFF";
  version = 0.1;
  TableFilename = "";
}

ChemtableNormalized::ChemtableNormalized(string filename, int nZm, int nZv, int nC, int nVar) : Table3D(nZm, nZv, nC, nVar)
{
  CombustionModel = "FPVA";
  ChemTableType = "COEFF";
  version = 1.0;
  TableFilename = filename;
}

//------------------------------//
//         Destructor           //
//------------------------------//
ChemtableNormalized::~ChemtableNormalized()
{
  Unload();
}
//------------------------------//
//            Methods           //
//------------------------------//
void ChemtableNormalized::AddVarName(string var)
{
  VarNames.push_back(var);
  int size = VarNames.size();

  pair<string,int> newElem(var,size-1);
  pair<map<string, int>::iterator, bool> test = VarNamesMap.insert( newElem );
  if (test.second == false) { cerr << "Variable "<< var <<" already added"<< endl; throw(-1); }
}

void ChemtableNormalized::WriteTable(string output)
{
  FILE   *inFp = 0;
  char   filename[kStrLong], buffer_m[kStrMed];
  size_t dum;

  if (mpi_rank == 0) cout << endl << "***** Write table to file *****" << endl << endl;

  // ***************
  // Open table file
  // ***************
  if (output == "debug") cout << "Opening file "<< TableFilename << endl;
  TableFilename.copy(&filename[0], kStrLong);
  dum = TableFilename.size();
  if (dum < kStrLong) { strcpy(&filename[dum], "\0"); }
  else { cerr << "### Chemistry table file name is too long! ###" << endl; throw(-1); }

  if (!(inFp = fopen(filename, "wb"))) { cerr << "### Could not open chemistry table file " << filename << " ###" << endl; throw(-1); }

  // ************************
  // Write table general info
  // ************************

  // Version
  if (output == "debug") cout << "Write version " << version <<endl;
  dum = fwrite(&version, sizeof(double), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: version ###" << endl; throw(-1); }

  // CombustionModel
  if (output == "debug") cout << "Write model " << CombustionModel << endl;
  CombustionModel.copy(&buffer_m[0], kStrMed);
  dum = CombustionModel.size();
  if (dum < kStrMed) {strcpy(&buffer_m[dum], "\0");}
  else { cerr << "### Combustion model name is too long! ###" << endl; throw(-1); }

  dum = fwrite(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed) { cerr << "### Error writing file: Combustion model ###" << endl; throw(-1);}

  // ChemTableType
  if (output == "debug") cout << "Write ChemTableType " << ChemTableType << endl;
  ChemTableType.copy(&buffer_m[0], kStrMed);
  dum = ChemTableType.size();
  if (dum < kStrMed) { strcpy(&buffer_m[dum], "\0"); }
  else { cerr << "### Chem Table type is too long! ###" << endl; throw(-1); }

  dum = fwrite(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed) {cerr << "### Error writing file: Chem Table type ###" << endl; throw(-1);}

  // Preference
  dum = fwrite(&Preference, sizeof(double), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: Preference ###" << endl; throw(-1); }

  // **********************
  // Write table dimensions
  // **********************
  int nZm = GetDimension1();
  if (output == "debug") cout << "Write nZm "<< nZm << endl;
  dum = fwrite(&nZm, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nZm ###" << endl; throw(-1); }

  int nZv = GetDimension2();
  if (output == "debug") cout << "Write nZv "<< nZv << endl;
  dum = fwrite(&nZv, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nZv ###" << endl; throw(-1); }

  int nC = GetDimension3();
  if (output == "debug") cout << "Write nC "<< nC << endl;
  dum = fwrite(&nC, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nC ###" << endl; throw(-1); }

  int nVar = GetNumberVariables();
  if (output == "debug") cout << "Write nVar "<< nVar << endl;
  dum = fwrite(&nVar, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nVar ###" << endl; throw(-1);}

  // ***********************
  // Write table coordinates 
  // ***********************
  if (output == "debug") cout << "Write Zm vector "<< endl;
  double *Zm = new double[nZm];
  CopyCoordinate1(Zm);
  dum = fwrite(Zm, sizeof(double), nZm, inFp);
  if (dum != nZm) { cerr << "### Error writing file: Zm ###" << endl; throw(-1);}

  if (output == "debug") cout << "Write Zv vector "<< endl;
  double *Zv = new double[nZv];
  CopyCoordinate2(Zv);
  dum = fwrite(Zv, sizeof(double), nZv, inFp);
  if (dum != nZv) { cerr << "### Error writing file: Zv ###" << endl; throw(-1);}

  if (output == "debug") cout << "Write C vector "<< endl;
  double *CC = new double[nC]; 
  CopyCoordinate3(CC);
  dum = fwrite(CC, sizeof(double), nC, inFp);
  if (dum != nC) { cerr << "### Error writing file: CC ###" << endl; throw(-1);}

  // ***********************
  // Write name of variables
  // ***********************
  for (int l = 0; l < nVar; l++)
  {
    if (output == "debug") cout << "Write VarName "<< VarNames[l] << endl;
    VarNames[l].copy(&buffer_m[0], kStrMed);
    dum = VarNames[l].size();
    if (dum < kStrMed) {strcpy(&buffer_m[dum], "\0");}
    else { cerr << "### Variable name " << VarNames[l] << " is too long! ###" << endl; throw(-1);}

    dum = fwrite(&buffer_m, sizeof(char), kStrMed, inFp);
    if (dum != kStrMed) { cerr << "### Error writing file: Variable name " << VarNames[l] << " ###" << endl; throw(-1);}
  }

  // ****************
  // Write table data 
  // ****************
  if (output == "debug") cout << "Write Table" << endl;
  for (int i = 0; i < nZm; i++)
    for (int j = 0; j < nZv; j++)
      for (int k = 0; k < nC; k++)
        for (int n = 0; n < nVar; n++)
        {
          float value = GetTableValue(i,j,k,n);
          dum = fwrite(&value, sizeof(float), 1, inFp);
          if (dum != 1) { cerr << "### Error reading file: Table -> " << i << " " << j << " " << k << " " << n << " ### " << endl; throw(-1);}
        }
  // ****************
  // Close table file
  // ****************
  if (output == "debug") cout << "Close file" << endl;
  fclose(inFp);
}

void ChemtableNormalized::Load(string inputfilename)
{
  TableFilename = inputfilename;
  if (mpi_rank == 0) cout << endl << "Loading Database: " << inputfilename << endl;

  FILE *inFp = NULL;
  char filename[kStrLong], buffer_m[kStrMed];
  size_t dum;
    
  /* Read input file for chemistry table information */
  TableFilename.copy(&filename[0], kStrLong);
  dum = TableFilename.size();
  if (dum < kStrLong) strcpy(&filename[dum], "\0");

  // Open chemistry table file
  if (!(inFp = fopen(filename, "rb"))) { cerr << "### Could not open input file " << filename << " ###" << endl; throw(-1); }
    
  // *****************
  // Read general info
  // *****************

  // Version
  double tableVersion;
  dum = fread(&tableVersion, sizeof(double), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: tableVersion ###" << endl; throw(-1); }
  if (tableVersion != version)
  {
    cerr << "### Error with version of chemistry table! Table was created with version " << tableVersion
         << " of CreateChemTable and cannot be read with version " << version << " of interpolation function! ###" << endl;
    throw(-1);
  }

  // CombustionModel
  for (int i=0; i<kStrMed-1; i++) strcpy(&buffer_m[i], " ");
  dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed) { cerr << "### Error reading file: Combustion model ###" << endl; throw(-1); }
  strcpy(&buffer_m[kStrMed - 1], "\0");
  CombustionModel.assign(buffer_m);

  // ChemTableType
  for (int i=0; i<kStrMed-1; i++) strcpy(&buffer_m[i], " ");
  dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed) { cerr << "### Error reading file: Chemistry table type ###" << endl; throw(-1); }
  strcpy(&buffer_m[kStrMed - 1], "\0");
  ChemTableType.assign(buffer_m);

  // Preference
  dum = fread(&Preference, sizeof(double), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: Preference ###" << endl; throw(-1); }

  // *********************
  // Read table dimensions
  // *********************
  int nZm;
  dum = fread(&nZm, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n1 ###" << endl; throw(-1); }

  int nZv;
  dum = fread(&nZv, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n2 ###" << endl; throw(-1); }

  int nC;
  dum = fread(&nC, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n3 ###" << endl; throw(-1); }

  int nVar;
  dum = fread(&nVar, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: nvar ###" << endl; throw(-1); }

  // **********************
  // Allocate table
  // **********************
  Allocate(nZm,nZv,nC,nVar);

  // **********************
  // Read table coordinates 
  // **********************
  double aux;

  // Mixture fraction dimension
  for (int i=0; i<nZm; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) {cerr << "### Error reading file: x1 ###" << endl; throw(-1); }
    SetCoordinate1(i, aux);
  }

  // Mixture fraction variance
  for (int i=0; i<nZv; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) { cerr << "### Error reading file: x2 ###" << endl; throw(-1); }
    SetCoordinate2(i, aux);
  }

  // Progress variable
  for (int i=0; i<nC; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) { cerr << "### Error reading file: x3 ###" << endl; throw(-1); }
    SetCoordinate3(i, aux);
  }

  // **********************
  // Read name of variables
  // **********************
  for (int l = 0; l < nVar; l++)
  {
    for (int i=0; i<kStrMed-1; i++) strcpy(&buffer_m[i], " ");
    dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
    if (dum != kStrMed) { cerr << "### Error reading file: VarNames no. " << l << " ###" << endl; throw(-1); }
    strcpy(&buffer_m[kStrMed - 1], "\0");
    AddVarName(buffer_m);
  }

  // ***************
  // Read table data 
  // ***************
  float aux2;
  for (int j=0; j<nZm; ++j)
    for (int k=0; k<nZv; ++k)
      for (int l=0; l<nC; ++l)
        for (int p=0; p<nVar; ++p)
        {
          dum = fread(&aux2, sizeof(float), 1, inFp);
          if (dum != 1) { cerr << "### Error reading file: data ->" << j << " " << k << " " << l << " " << p << " ### " << endl; throw(-1);}
          SetTableValue(j, k, l, p, aux2); 
        }

  fclose(inFp);
}

void ChemtableNormalized::Unload(void) {
  DestroyTable();
  TableFilename.clear();
  CombustionModel.clear();
  ChemTableType.clear();
  VarNames.clear();
  VarNamesMap.clear();
}

double ChemtableNormalized::Lookup(double Zm, double Zvar, double C, string variable) {
  int ivar = GetVarNameIndex(variable);
  return Interpolate(Zm, Zvar, C, ivar);
}

double ChemtableNormalized::Lookup(double Zm, double Zvar, double C, int ivar) { return Interpolate(Zm, Zvar, C, ivar); }

double ChemtableNormalized::Lookup(string variable) {
  int ivar = GetVarNameIndex(variable);
  return Interpolate(ivar);
}

double ChemtableNormalized::Lookup(int ivar) { return Interpolate(ivar); }

void ChemtableNormalized::LookupCoeff(double &RoM,      double &T0,      double &E0,
                                double &Gam0,     double &a_Gam,   double &mu0,
                                double &a_mu,     double &lamOcp0, double &a_lamOcp,
                                double &src_prog, double Zm,
                                double Zvar,      double C) {
  InterpolatePoint(Zm, Zvar, C);
  RoM      = Interpolate(0);
  T0       = Interpolate(1);
  E0       = Interpolate(2);
  Gam0     = Interpolate(3);
  a_Gam    = Interpolate(4);
  mu0      = Interpolate(5);
  a_mu     = Interpolate(6);
  lamOcp0  = Interpolate(7);
  a_lamOcp = Interpolate(8);
  src_prog = Interpolate(9);
}

void ChemtableNormalized::LookupSelectedCoeff(double &RoM,  double &T0,    double &E0,
                                        double &Gam0, double &a_Gam,
                                        double Zm,    double Zvar,   double C) {
  InterpolatePoint(Zm, Zvar, C);
  RoM      = Interpolate(0);
  T0       = Interpolate(1);
  E0       = Interpolate(2);
  Gam0     = Interpolate(3);
  a_Gam    = Interpolate(4);
}

void ChemtableNormalized::print(void) {
	int nZm    = GetDimension1();
	int nZv    = GetDimension2();
	int nC     = GetDimension3();
	int nVar   = GetNumberVariables();

	cout << endl;
	cout << "------------------------------------------------------" << endl;
	cout << "  Multi-Pressure Chemtable: " << TableFilename << endl;
	cout << "------------------------------------------------------" << endl;
  cout << "  - Chemistry table type : " << ChemTableType << endl;
  cout << "  - Combustion model     : " << CombustionModel << endl;
  cout << "  - Reference pressure is: " << Preference << endl;
  cout << "  - Chemistry table size : " << nVar   << " variables" << endl;
  cout << "                           " << nZm    << " (Mean Mixture Fraction) "<< endl;
  cout << "                           " << nZv    << " (Mixture Fraction Variance) " << endl;
  cout << "                           " << nC     << " (Progress Variable) " << endl;
  cout << endl;

  cout << "  - Input coordinates    : " << endl;
  //cout << setiosflags(ios::fixed);
  cout << setprecision(4);
  //cout << dec;
	// Mixture Fraction Vec
	double *vec = new double[nZm];
	CopyCoordinate1(vec);
	cout  << "    x1 : Zmean      =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nZm-1] << "]" << endl;
	delete [] vec;
	// Mixture Fraction Variance vec
	vec = new double[nZv];
	CopyCoordinate2(vec);
	cout  << "    x2 : Zvar       =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nZv-1] << "]" << endl;
	delete [] vec;
	// Progress Variable vec
	vec = new double[nC];
	CopyCoordinate3(vec);
  cout  << "    x3 : Prog       =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nC-1] << "]" << endl;
  delete [] vec;
  cout << endl;

  // Display variables names, min and max
  cout << "  - Stored variables     : " << endl;
	float minval, maxval;
	int mi1, mi2, mi3, ma1, ma2, ma3;
	for (int l = 0; l < nVar; l++)
	{
		GetArrayInfo(l, minval, maxval, mi1, mi2, mi3, ma1, ma2, ma3);
		cout.width(6);
		cout << l << ". ";
		cout.width(12);
		cout << VarNames[l] << "  =  ["; cout.width(12); cout << minval << " ... "; cout.width(12); cout << maxval << "]" << endl;
	}
	cout << "------------------------------------------------------" << endl;
	cout << endl;
}

double ChemtableNormalized::ComputeNormalizedVariance(const double Zm, const double Zvar) {
  double Sz;
  double eps = 1.0e-5;
  if ((Zm <= 1.0-eps)&&(Zm >= eps)) {
    Sz = Zvar/(Zm*(1.0-Zm));
  } else {
    Sz = 0.0;
  }
  Sz = max(0.0, min(1.0, Sz));
  return Sz;
}

double ChemtableNormalized::ComputeNormalizedProg(const double Zm, const double Sz, const double Yc) {
  double eps = 1.0e-5;
  double Yc_eq = Lookup(Zm, Sz, 1.0, "PROG");
  double C_out;
  if (Yc_eq >= eps) {
    C_out = Yc/Yc_eq;
  } else {
    C_out = 0.0;
  }
  return C_out;
}

void ChemtableNormalized::BCAST_table(int mpi_root, MPI_Comm mpi_comm) {
  // Broadcast Table3D attributes
  BCAST(mpi_root, mpi_comm);

  // TableFileName, version, CombustionModel, ChemTableType are not used during the computation and are therefore not broadcast
  // TODO: Broadcast them anyway
  char buffer_m[kStrMed];
  int nVar = GetNumberVariables();
  size_t dum;
  for (int i=0; i<nVar; ++i) {
    if (mpi_rank == mpi_root) {
      for (int j=0; j<kStrMed-1; ++j) strcpy(&buffer_m[j]," ");
      VarNames[i].copy(&buffer_m[0], kStrMed);
      dum = VarNames[i].size();
      strcpy(&buffer_m[dum],"\0");
    }
    MPI_Bcast(buffer_m,kStrMed,MPI_CHAR,mpi_root,mpi_comm);    
    if (mpi_rank != mpi_root) AddVarName(buffer_m);
  }
  // Broadcast Preference
  double Pref;
  if (mpi_rank == mpi_root) Pref = GetReferencePressure(); 
  MPI_Bcast(&Pref,1,MPI_DOUBLE,mpi_root,mpi_comm);
  if (mpi_rank != mpi_root) SetReferencePressure(Pref);
}

void ChemtableNormalized::MPI_Load(string inputfilename) {
  int mpi_root = 0;
  // Only master process reads database file
  if (mpi_rank == mpi_root) Load(inputfilename);
  // Everyone waits for master to finish reading
  MPI_Barrier(mpi_comm);
  // Broadcast database 
  BCAST_table(mpi_root, mpi_comm);
  MPI_Barrier(mpi_comm);
}

void ChemtableNormalized::WriteTableTecplot(void)
/* Write table to file in Tecplot ASCII format for visualization */
{
  ofstream fout;
  char     filename[kStrLong] = "database.dat";
  cout << endl << " Write table to tecplot format file ... ";

  int nZm    = GetDimension1();
  int nZv    = GetDimension2();
  int nC     = GetDimension3();
  int nVar   = GetNumberVariables();
  // Mixture Fraction Vec
  double *Zm = new double[nZm];
  CopyCoordinate1(Zm);
  // Mixture Fraction Variance vec
  double *Zv = new double[nZv];
  CopyCoordinate2(Zv);
  // Progress Variable vec
  double *C = new double[nC];
  CopyCoordinate3(C);

  // *********
  // Open file
  // *********
  fout.open(filename);

  // ********************
  // Write tecplot header
  // ********************
  fout << "TITLE=\"Chemistry table\"" << endl;
  fout << "VARIABLES=\"I1\" \"I2\" \"I3\" \"ZMean\" \"ZVar\" \"Progress\" ";
  for (int j=0; j<nVar; j++)
    fout << "\"" << VarNames[j] << "\" ";
  fout << endl;

  fout << "ZONE I=" << nZm << ", J=" << nZv << ", K=" << nC << ", DATAPACKING=POINT" << endl;

  // *******************************
  // Write data for each table point
  // *******************************
  for (int k=0; k<nC; k++)
    for (int j=0; j<nZv; j++)
      for (int i=0; i<nZm; i++)
      {
        fout << i << "  " << j << "  " << k << "  ";
        fout << Zm[i] << "   " << Zv[j] << "   " << C[k];
        for (int n=0; n<nVar; n++)
          fout << "   " << GetTableValue(i,j,k,n);
        fout << endl;
      }

  fout.close();

  delete [] Zm;
  delete [] Zv;
  delete [] C;
  cout << "done !" << endl;
}

