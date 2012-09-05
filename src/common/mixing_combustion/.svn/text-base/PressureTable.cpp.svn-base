#include "PressureTable.h"

// static/extern/global variables
double PressureTable::P_scale = 1.0e5;
PressureTable* PressureTable::ptr_for_newton;

//------------------------------//
//         Constructors         //
//------------------------------//
PressureTable::PressureTable()
{
  CombustionModel = "FPVA";
  ChemTableType = "COEFF";
  version = 0.1;
  TableFilename = "";
  Z_newton = -1.0;
  Zv_newton = -1.0;
  rho_newton = -1.0;
  Yc_newton = -1.0;
  C_newton = -1.0;
  P_newton = -1.0;
  // default behavior: the table used in the newton wrapper is the last table object that has been declared
  ptr_for_newton = this;
}

PressureTable::PressureTable(string filename, int nPress, int nZm, int nZv, int nC, int nVar) : Table4D(nPress, nZm, nZv, nC, nVar)
{
  CombustionModel = "FPVA";
  ChemTableType = "COEFF";
  version = 0.1;
  TableFilename = filename;
  Z_newton = -1.0;
  Zv_newton = -1.0;
  rho_newton = -1.0;
  Yc_newton = -1.0;
  C_newton = -1.0;
  P_newton = -1.0;
  // default behavior: the table used in the newton wrapper is the last table object that has been declared
  ptr_for_newton = this;
}

//------------------------------//
//         Destructor           //
//------------------------------//
PressureTable::~PressureTable()
{
  Unload();
}
//------------------------------//
//            Methods           //
//------------------------------//
void PressureTable::AddVarName(string var)
{
  VarNames.push_back(var);
  int size = VarNames.size();

  pair<string,int> newElem(var,size-1);
  pair<map<string, int>::iterator, bool> test = VarNamesMap.insert( newElem );
  if (test.second == false) { cerr << "Variable "<< var <<" already added"<< endl; throw(-1); }
}

void PressureTable::WriteTable(string output)
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

  // **********************
  // Write table dimensions
  // **********************
  int nPress = GetDimension1();
  if (output == "debug") cout << "Write nPress "<< nPress << endl;
  dum = fwrite(&nPress, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nPress ###" << endl; throw(-1); }

  int nZm = GetDimension2();
  if (output == "debug") cout << "Write nZm "<< nZm << endl;
  dum = fwrite(&nZm, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nZm ###" << endl; throw(-1); }

  int nZv = GetDimension3();
  if (output == "debug") cout << "Write nZv "<< nZv << endl;
  dum = fwrite(&nZv, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error writing file: nZv ###" << endl; throw(-1); }

  int nC = GetDimension4();
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
  if (output == "debug") cout << "Write Pressure vector "<< endl;
  double *Pressure = new double[nPress];
  CopyCoordinate1(Pressure);
  dum = fwrite(Pressure, sizeof(double), nPress, inFp);
  if (dum != nPress) { cerr << "### Error writing file: Pressure ###" << endl; throw(-1);}

  if (output == "debug") cout << "Write Zm vector "<< endl;
  double *Zm = new double[nZm];
  CopyCoordinate2(Zm);
  dum = fwrite(Zm, sizeof(double), nZm, inFp);
  if (dum != nZm) { cerr << "### Error writing file: Zm ###" << endl; throw(-1);}

  if (output == "debug") cout << "Write Zv vector "<< endl;
  double *Zv = new double[nZv];
  CopyCoordinate3(Zv);
  dum = fwrite(Zv, sizeof(double), nZv, inFp);
  if (dum != nZv) { cerr << "### Error writing file: Zv ###" << endl; throw(-1);}

  if (output == "debug") cout << "Write C vector "<< endl;
  double *CC = new double[nC]; 
  CopyCoordinate4(CC);
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
  for (int ipress=0; ipress <nPress; ipress++)
    for (int i = 0; i < nZm; i++)
      for (int j = 0; j < nZv; j++)
        for (int k = 0; k < nC; k++)
          for (int n = 0; n < nVar; n++)
          {
            float value = GetTableValue(ipress,i,j,k,n);
            dum = fwrite(&value, sizeof(float), 1, inFp);
            if (dum != 1) { cerr << "### Error reading file: Table -> " << ipress << " " << i << " " << j << " " << k << " " << n << " ### " << endl; throw(-1);}
          }
  // ****************
  // Close table file
  // ****************
  if (output == "debug") cout << "Close file" << endl;
  fclose(inFp);
}

void PressureTable::Load(string inputfilename)
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

  // *********************
  // Read table dimensions
  // *********************
  int nPress;
  dum = fread(&nPress, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n1 ###" << endl; throw(-1); }

  int nZm;
  dum = fread(&nZm, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n2 ###" << endl; throw(-1); }

  int nZv;
  dum = fread(&nZv, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n3 ###" << endl; throw(-1); }

  int nC;
  dum = fread(&nC, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: n4 ###" << endl; throw(-1); }

  int nVar;
  dum = fread(&nVar, sizeof(int), 1, inFp);
  if (dum != 1) { cerr << "### Error reading file: nvar ###" << endl; throw(-1); }

  // **********************
  // Allocate table
  // **********************
  Allocate(nPress,nZm,nZv,nC,nVar);

  // **********************
  // Read table coordinates 
  // **********************
  double aux;

  // Pressure dimension
  for (int i=0; i<nPress; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) {cerr << "### Error reading file: x1 ###" << endl; throw(-1); }
    SetCoordinate1(i, aux);
  }

  // Mixture fraction dimension
  for (int i=0; i<nZm; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) {cerr << "### Error reading file: x2 ###" << endl; throw(-1); }
    SetCoordinate2(i, aux);
  }

  // Mixture fraction variance
  for (int i=0; i<nZv; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) { cerr << "### Error reading file: x3 ###" << endl; throw(-1); }
    SetCoordinate3(i, aux);
  }

  // Progress variable
  for (int i=0; i<nC; ++i)
  {
    dum = fread(&aux, sizeof(double), 1, inFp);
    if (dum != 1) { cerr << "### Error reading file: x4 ###" << endl; throw(-1); }
    SetCoordinate4(i, aux);
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
  for (int i=0; i<nPress; ++i)
    for (int j=0; j<nZm; ++j)
      for (int k=0; k<nZv; ++k)
        for (int l=0; l<nC; ++l)
          for (int p=0; p<nVar; ++p)
          {
            dum = fread(&aux2, sizeof(float), 1, inFp);
            if (dum != 1) { cerr << "### Error reading file: data ->" << i << " " << j << " " << k << " " << l << " " << p << " ### " << endl; throw(-1);}
            SetTableValue(i, j, k, l, p, aux2); 
          }

  fclose(inFp);
}

void PressureTable::Unload(void) {
  DestroyTable();
  TableFilename.clear();
  CombustionModel.clear();
  ChemTableType.clear();
  VarNames.clear();
  VarNamesMap.clear();
}

double PressureTable::Lookup(double P, double Zm, double Zvar, double C, string variable) {
  int ivar = GetVarNameIndex(variable); 
  return Interpolate(P, Zm, Zvar, C, ivar);
}

double PressureTable::Lookup(double P, double Zm, double Zvar, double C,int ivar) {return Interpolate(P, Zm, Zvar, C, ivar);}

double PressureTable::Lookup(string variable) {
  int ivar = GetVarNameIndex(variable); 
  return Interpolate(ivar);
}

double PressureTable::Lookup(int ivar) {return Interpolate(ivar);}

void PressureTable::LookupCoeff(double &RoM,      double &T0,      double &E0,
                                double &Gam0,     double &a_Gam,   double &mu0,
                                double &a_mu,     double &lamOcp0, double &a_lamOcp,
                                double &src_prog, double P,        double Zm,
                                double Zvar,      double C) {
  InterpolatePoint(P, Zm, Zvar, C);
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

void PressureTable::LookupSelectedCoeff(double &RoM,  double &T0,    double &E0,
                                        double &Gam0, double &a_Gam, double P,
                                        double Zm,    double Zvar,   double C) {
  InterpolatePoint(P, Zm, Zvar, C);
  RoM      = Interpolate(0);
  T0       = Interpolate(1);
  E0       = Interpolate(2);
  Gam0     = Interpolate(3);
  a_Gam    = Interpolate(4);
}

void PressureTable::print(void) {
	int nPress = GetDimension1();
	int nZm    = GetDimension2();
	int nZv    = GetDimension3();
	int nC     = GetDimension4();
	int nVar   = GetNumberVariables();

	cout << endl;
	cout << "------------------------------------------------------" << endl;
	cout << "  Multi-Pressure Chemtable: " << TableFilename << endl;
	cout << "------------------------------------------------------" << endl;
  cout << "  - Chemistry table type : " << ChemTableType << endl;
  cout << "  - Combustion model     : " << CombustionModel << endl;
  cout << "  - Chemistry table size : " << nVar   << " variables" << endl;
  cout << "                           " << nPress << " (Pressure) " << endl;
  cout << "                           " << nZm    << " (Mean Mixture Fraction) "<< endl;
  cout << "                           " << nZv    << " (Mixture Fraction Variance) " << endl;
  cout << "                           " << nC     << " (Progress Variable) " << endl;
  cout << endl;

  cout << "  - Input coordinates    : " << endl;
  //cout << setiosflags(ios::fixed);
  cout << setprecision(4);
  //cout << dec;
  // Pressure vec
	double *vec = new double[nPress];
	CopyCoordinate1(vec);
	cout  <<"    x1 : Pressure   =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nPress-1] << "]" << endl;
	delete [] vec;
	// Mixture Fraction Vec
	vec = new double[nZm];
	CopyCoordinate2(vec);
	cout  << "    x2 : Zmean      =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nZm-1] << "]" << endl;
	delete [] vec;
	// Mixture Fraction Variance vec
	vec = new double[nZv];
	CopyCoordinate3(vec);
	cout  << "    x3 : Zvar       =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nZv-1] << "]" << endl;
	delete [] vec;
	// Progress Variable vec
	vec = new double[nC];
	CopyCoordinate4(vec);
  cout  << "    x4 : Prog       =" << "\t ["; cout.width(6); cout << vec[0] << " ... "; cout.width(6); cout << vec[nC-1] << "]" << endl;
  delete [] vec;
  cout << endl;

  // Display variables names, min and max
  cout << "  - Stored variables     : " << endl;
	float minval, maxval;
	int mi1, mi2, mi3, ma1, ma2, ma3, mi4, ma4;
	for (int l = 0; l < nVar; l++)
	{
		GetArrayInfo(l, minval, maxval, mi1, mi2, mi3, mi4, ma1, ma2, ma3, ma4);
		cout.width(6);
		cout << l << ". ";
		cout.width(12);
		cout << VarNames[l] << "  =  ["; cout.width(12); cout << minval << " ... "; cout.width(12); cout << maxval << "]" << endl;
	}
	cout << "------------------------------------------------------" << endl;
	cout << endl;
}

void PressureTable::ComputeDatabaseInputs( const double rho, const double Z, const double Zv, const double Yc,
                                           double &P,  double &Sz , double &C,
                                           double P0,  double C0) {
  // Compute normalized variance
  double eps = 1.0e-3;
  if ((Z <= 1-eps)||(Z >= eps)) {
    Sz = Zv/(Z*(1.0-Z));
  } else {
    Sz = 0.0;
  }
  Sz = max(0.0, min(1.0, Sz));

  // Compute normalized progress variable
  Z_newton = Z;
  Zv_newton = Sz;
  Yc_newton = Yc;
  rho_newton = rho;
  // I need to solve ( Yc(C,P) - Yc_target , rho(C,P) - rho_target) = (0 , 0)
  C_newton = C0;
  P_newton = P0/P_scale;// Pressure is rescaled
  verbose_newton = false;
  Bool check;
  VecDoub_IO sol2(1, 0.0);
  VecDoub_IO sol(2, 0.0);
  int Newton_mode;

  try {
    //Newton  : newt(sol, check, usrfunc_wrapper);
    //Broyden : broydn(sol, check, usrfunc_wrapper);

    // Get Approximation of Yceq from highest pressure
    int nPress = GetDimension1();
    double P_max = GetCoordinate1(nPress-1);
    int nC = GetDimension4();
    double C_max = GetCoordinate4(nC-2); // Maximum C value before one
    double Yc_eq  = Lookup(P_max,Z_newton,Zv_newton,1.0,"PROG");
    double Yc_max = Lookup(P_max,Z_newton,Zv_newton,C_max,"PROG");
    // If Yc_eq is small, C is set to 1 and we look for pressure only
    if ((Yc_eq < eps)||(Yc > Yc_max)) {
      Newton_mode = 1;
      sol2[0] = P_newton;
      C_newton = 1.0;
      newt(sol2, check, usrfunc_wrapper2);
    } else {
      // usual
      Newton_mode = 2;
      sol[0] = C_newton;
      sol[1] = P_newton;
      try {newt(sol, check, usrfunc_wrapper);}
      catch (int e) {
       cout << "relaunched with Newton_mode = 1" <<  endl;
       Newton_mode = 1;
       C_newton = sol[0];
       sol2[0] = P0/P_scale;
       newt(sol2, check, usrfunc_wrapper2);
      }
    }
  }
  catch (int e) {
    int tol_failure = 0;
    if (Newton_mode == 1) {
      VecDoub out(1, 0.0);
      out = DensityFromP(sol2);
      if (abs(out[0]) > 1.0e-5) tol_failure = 1;
    } else if (Newton_mode == 2) {
      VecDoub out(2, 0.0);
      out = YcAndDensityFromCP(sol);
     if ((abs(out[0])>1.0e-4)||(abs(out[1]))>1.0e-5) tol_failure = 1;
    }
    cout << "Newton error. Test tol = " << tol_failure << endl;
    if (tol_failure == 1) { 
    cout << "Newton inputs : " << Z_newton   << " "
                               << Zv_newton  << " "
                               << Yc_newton  << " "
                               << rho_newton << " "
                               << C_newton   << " "
                               << P_newton*P_scale   << endl;
    cout << "C0 : " << C0 << " P0 : " << P0 << endl;
    cout << " Check : " << check << endl;
    if (Newton_mode == 2) {
      cout << " sol : C = " << sol[0] << " , P = " << sol[1]*P_scale << endl;
    } else if (Newton_mode == 1) {
      cout << " sol :  P = " << sol2[0]*P_scale << endl;
    }

    // Relaunch with verbose
    verbose_newton = true;
    try {
      Int MaxIter = 200;
      if (Newton_mode == 1) {
        sol2[0] = P_newton;
        C_newton = 1.0;
        cout << "Reinit with " <<  sol2[0] << " " << MaxIter <<  endl;
        newt(sol2, check, usrfunc_wrapper2, MaxIter);
      } else if (Newton_mode == 2) {
        // usual
        sol[0] = C_newton;
        sol[1] = P_newton;
        cout << "Reinit with " <<  sol[0] << " " << sol[1] << " " << MaxIter << endl;
        newt(sol, check, usrfunc_wrapper, MaxIter);
      } else {cerr << "Impossible value of Newton_mode: " << Newton_mode << endl;}

    }
    catch (int e2) {
      cout << "Bis Newton inputs : " << Z_newton   << " "
                                 << Zv_newton  << " "
                                 << Yc_newton  << " "
                                 << rho_newton << " "
                                 << C_newton   << " "
                                 << P_newton*P_scale   << endl;
      cout << "Bis Check : " << check << endl;
      if (Newton_mode == 2) {
        cout << "Bis sol : C = " << sol[0] << " , P = " << sol[1]*P_scale << endl;
      } else if (Newton_mode == 1) {
        cout << "Bis sol :  P = " << sol2[0]*P_scale << endl;
      }  
      throw(-1);
    }
    throw(-1);
    }
  }

  if (Newton_mode == 2) {
    C_newton = sol[0];
    P_newton = sol[1];
  } else if (Newton_mode == 1) {
    P_newton = sol2[0];
  }
  C = C_newton;
  P = P_newton*P_scale;
}

void PressureTable::ComputeDatabaseInputs_new( const double rho, const double Z, const double Zv, const double Yc, const double rhoE,
                                           double &P,  double &Sz , double &C,
                                           double P0,  double C0) {
  // Compute normalized variance
  double eps = 1.0e-3;
  Sz = ComputeNormalizedVariance(Z, Zv);

  // Compute normalized progress variable
  Z_newton = Z;
  Zv_newton = Sz;
  Yc_newton = Yc;
  rho_newton = rho;
  E_newton = rhoE/rho;
  // I need to solve ( Yc(C,P) - Yc_target , rho(C,P) - rho_target) = (0 , 0)
  C_newton = C0;
  P_newton = P0/P_scale;// Pressure is rescaled
  verbose_newton = false;
  Bool check;
  VecDoub_IO sol2(1, 0.0);
  VecDoub_IO sol(2, 0.0);
  int Newton_mode;


  try {
    // Get Approximation of Yceq from highest pressure
    int nPress = GetDimension1();
    double P_max = GetCoordinate1(nPress-1);
    int nC = GetDimension4();
    double C_max = GetCoordinate4(nC-2); // Maximum C value before one
    double Yc_eq  = Lookup(P_max,Z_newton,Zv_newton,1.0,"PROG");
    double Yc_max = Lookup(P_max,Z_newton,Zv_newton,C_max,"PROG");
    // If Yc_eq is small, C is set to 1 and we look for pressure only
    if ((Yc_eq < eps)||(Yc > Yc_max)) {
      Newton_mode = 1;
      sol2[0] = P_newton;
      C_newton = 1.0;
      newt(sol2, check, usrfunc_wrapper2);
    } else {
      // usual
      Newton_mode = 2;
      sol[0] = C_newton;
      sol[1] = P_newton;
      try {newt(sol, check, usrfunc_wrapper);}
      catch (int e) {
       cout << "relaunched with Newton_mode = 1" <<  endl;
       Newton_mode = 1;
       C_newton = sol[0];
       sol2[0] = P0/P_scale;
       newt(sol2, check, usrfunc_wrapper2);
      }
    }
  }
  catch (int e) {
    int tol_failure = 0;
    if (Newton_mode == 1) {
      VecDoub out(1, 0.0);
      out = DensityFromP(sol2);
      if (abs(out[0]) > 1.0e-5) tol_failure = 1;
    } else if (Newton_mode == 2) {
      VecDoub out(2, 0.0);
      out = YcAndDensityFromCP(sol);
     if ((abs(out[0])>1.0e-4)||(abs(out[1]))>1.0e-5) tol_failure = 1;
    }
    cout << "Newton error. Test tol = " << tol_failure << endl;
    if (tol_failure == 1) { 
    cout << "Newton inputs : " << Z_newton   << " "
                               << Zv_newton  << " "
                               << Yc_newton  << " "
                               << rho_newton << " "
                               << C_newton   << " "
                               << P_newton*P_scale   << endl;
    cout << "C0 : " << C0 << " P0 : " << P0 << endl;
    cout << " Check : " << check << endl;
    if (Newton_mode == 2) {
      cout << " sol : C = " << sol[0] << " , P = " << sol[1]*P_scale << endl;
    } else if (Newton_mode == 1) {
      cout << " sol :  P = " << sol2[0]*P_scale << endl;
    }
    }
  }
  if (Newton_mode == 2) {
    C_newton = sol[0];
    P_newton = sol[1];
  } else if (Newton_mode == 1) {
    P_newton = sol2[0];
  }
  C = C_newton;
  P = P_newton*P_scale;

}

void PressureTable::ComputeDatabaseInputs_FromT( const double rho, const double Z, const double Zv, const double Yc, const double T,
                                           double &P,  double &Sz , double &C,
                                           double P0,  double C0) {
  // Compute normalized variance
  double eps = 1.0e-3;
  Sz = ComputeNormalizedVariance(Z, Zv);

  // Compute normalized progress variable
  Z_newton = Z;
  Zv_newton = Sz;
  Yc_newton = Yc;
  rho_newton = rho;
  T_newton = T;
  // I need to solve ( Yc(C,P) - Yc_target , rho(C,P) - rho_target) = (0 , 0)
  C_newton = C0;
  P_newton = P0/P_scale;// Pressure is rescaled
  verbose_newton = false;
  Bool check;
  VecDoub_IO sol2(1, 0.0);
  VecDoub_IO sol(2, 0.0);
  int Newton_mode;


  try {
    // Get Approximation of Yceq from highest pressure
    int nPress = GetDimension1();
    double P_max = GetCoordinate1(nPress-1);
    int nC = GetDimension4();
    double C_max = GetCoordinate4(nC-2); // Maximum C value before one
    double Yc_eq  = Lookup(P_max,Z_newton,Zv_newton,1.0,"PROG");
    double Yc_max = Lookup(P_max,Z_newton,Zv_newton,C_max,"PROG");
    // If Yc_eq is small, C is set to 1 and we look for pressure only
    if ((Yc_eq < eps)||(Yc > Yc_max)) {
      Newton_mode = 1;
      sol2[0] = P_newton;
      C_newton = 1.0;
      newt(sol2, check, usrfunc_wrapper4);
    } else {
      // usual
      Newton_mode = 2;
      sol[0] = C_newton;
      sol[1] = P_newton;
      try {newt(sol, check, usrfunc_wrapper3);}
      catch (int e) {
       cout << "relaunched with Newton_mode = 1" <<  endl;
       Newton_mode = 1;
       C_newton = sol[0];
       sol2[0] = P0/P_scale;
       newt(sol2, check, usrfunc_wrapper4);
      }
    }
  }
  catch (int e) {
    int tol_failure = 0;
    if (Newton_mode == 1) {
      VecDoub out(1, 0.0);
      out = DensityFromP(sol2);
      if (abs(out[0]) > 1.0e-5) tol_failure = 1;
    } else if (Newton_mode == 2) {
      VecDoub out(2, 0.0);
      out = YcAndDensityFromCP(sol);
     if ((abs(out[0])>1.0e-4)||(abs(out[1]))>1.0e-5) tol_failure = 1;
    }
    cout << "Newton error. Test tol = " << tol_failure << endl;
    if (tol_failure == 1) {
    cout << "Newton inputs : " << Z_newton   << " "
                               << Zv_newton  << " "
                               << Yc_newton  << " "
                               << rho_newton << " "
                               << C_newton   << " "
                               << P_newton*P_scale   << endl;
    cout << "C0 : " << C0 << " P0 : " << P0 << endl;
    cout << " Check : " << check << endl;
    if (Newton_mode == 2) {
      cout << " sol : C = " << sol[0] << " , P = " << sol[1]*P_scale << endl;
    } else if (Newton_mode == 1) {
      cout << " sol :  P = " << sol2[0]*P_scale << endl;
    }
    }
  }
  if (Newton_mode == 2) {
    C_newton = sol[0];
    P_newton = sol[1];
  } else if (Newton_mode == 1) {
    P_newton = sol2[0];
  }
  C = C_newton;
  P = P_newton*P_scale;

}

VecDoub PressureTable::YcAndDensityFromCP(VecDoub_I vec) {
  double eps = 1.0e-4;
  double C = vec[0];
  double Press = vec[1]*P_scale;
  double Yc;
  double rho;
  VecDoub vec_out(2, 0.0);

  double Yc_eq  = Lookup(Press,Z_newton,Zv_newton,1.0,"PROG");
  Yc = C*Yc_eq;
  // Define Yc to find the solution C=0 when Yc_eq = 0 (necessary for Z=0 or 1, Sz = 1)
  if (Yc_eq == 0.0) {
    Yc = Yc_newton + C;
  } else if (Yc_eq < eps) {
  //                                 or when Yc_eq < eps
    Yc = Yc_newton*C;
  }

  int nPress = GetDimension1();
  double P_min = GetCoordinate1(0);
  double P_max = GetCoordinate1(nPress-1);

  if (Press <= P_max) {
    if (Press >= P_min) {
      rho = Lookup(Press,Z_newton,Zv_newton,C,"rho0");
    } else {
      double rho_min = Lookup(P_min,Z_newton,Zv_newton,C,"rho0");
      rho = rho_min*Press/P_min;
    }
  } else {
    double rho_max = Lookup(P_max,Z_newton,Zv_newton,C,"rho0");
    rho = rho_max*Press/P_max;
  }


  vec_out[0] = Yc-Yc_newton;
  vec_out[1] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << C << " , " << Press << " ) => ( "
         << vec_out[0] << " , " << vec_out[1] << " )"
         << endl;
    cout << "Yc : " << Yc << " vs " << Yc_newton << endl;
    cout << "rho : " << rho << " vs " << rho_newton << endl;
  }
  return vec_out;
}

VecDoub PressureTable::YcAndDensityFromCP_new(VecDoub_I vec) {
  double eps = 1.0e-4;
  double C = vec[0];
  double Press = vec[1]*P_scale;
  double Yc;
  double rho;
  VecDoub vec_out(2, 0.0);

  double Yc_eq  = Lookup(Press,Z_newton,Zv_newton,1.0,"PROG");
  Yc = C*Yc_eq;
  // Define Yc to find the solution C=0 when Yc_eq = 0 (necessary for Z=0 or 1, Sz = 1)
  if (Yc_eq == 0.0) {
    Yc = Yc_newton + C;
  } else if (Yc_eq < eps) {
  //                                 or when Yc_eq < eps
    Yc = Yc_newton*C;
  }

  double RoM, T0, E0, GAMMA0, AGAMMA;
  LookupSelectedCoeff(RoM, T0, E0, GAMMA0, AGAMMA, Press, Z_newton,Zv_newton,C);
  double Temp;
  if (AGAMMA == 0.0)
    Temp = T0 + (GAMMA0 - 1.0) / RoM * (E_newton - E0);
  else
    Temp = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E_newton - E0) / RoM) - 1.0);

  rho = Press/(RoM*Temp);

  vec_out[0] = Yc-Yc_newton;
  vec_out[1] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << C << " , " << Press << " ) => ( "
         << vec_out[0] << " , " << vec_out[1] << " )"
         << endl;
    cout << "Yc : " << Yc << " vs " << Yc_newton << endl;
    cout << "rho : " << rho << " vs " << rho_newton << endl;
  }
  return vec_out;
}

VecDoub PressureTable::YcAndDensityFromCPT(VecDoub_I vec) {
  double eps = 1.0e-4;
  double C = vec[0];
  double Press = vec[1]*P_scale;
  double Yc;
  double rho;
  VecDoub vec_out(2, 0.0);

  double Yc_eq  = Lookup(Press,Z_newton,Zv_newton,1.0,"PROG");
  Yc = C*Yc_eq;
  // Define Yc to find the solution C=0 when Yc_eq = 0 (necessary for Z=0 or 1, Sz = 1)
  if (Yc_eq == 0.0) {
    Yc = Yc_newton + C;
  } else if (Yc_eq < eps) {
  //                                 or when Yc_eq < eps
    Yc = Yc_newton*C;
  }

  double RoM = Lookup(Press, Z_newton,Zv_newton,C, "RoM");
  rho = Press/(RoM*T_newton);

  vec_out[0] = Yc-Yc_newton;
  vec_out[1] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << C << " , " << Press << " ) => ( "
         << vec_out[0] << " , " << vec_out[1] << " )"
         << endl;
    cout << "Yc : " << Yc << " vs " << Yc_newton << endl;
    cout << "rho : " << rho << " vs " << rho_newton << endl;
  }
  return vec_out;
}

VecDoub PressureTable::DensityFromP(VecDoub_I vec) {
  double Press = vec[0]*P_scale;
  double rho;
  VecDoub vec_out(1, 0.0);

  int nPress = GetDimension1();
  double P_min = GetCoordinate1(0);
  double P_max = GetCoordinate1(nPress-1);

  if (Press <= P_max) {
    if (Press >= P_min) {
      rho = Lookup(Press,Z_newton,Zv_newton,C_newton,"rho0");
    } else {
      double rho_min = Lookup(P_min,Z_newton,Zv_newton,C_newton,"rho0");
      rho = rho_min*Press/P_min;
    }
  } else {
    double rho_max = Lookup(P_max,Z_newton,Zv_newton,C_newton,"rho0");
    rho = rho_max*Press/P_max;
  }

  vec_out[0] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << Press << " ) => ( " 
         << rho << " vs " << rho_newton << " )"
         << endl;
  }
  return vec_out;
}

VecDoub PressureTable::DensityFromP_new(VecDoub_I vec) {
  double Press = vec[0]*P_scale;
  double rho;
  VecDoub vec_out(1, 0.0);

  double RoM, T0, E0, GAMMA0, AGAMMA;
  LookupSelectedCoeff(RoM, T0, E0, GAMMA0, AGAMMA, Press, Z_newton,Zv_newton,C_newton);
  double Temp;
  if (AGAMMA == 0.0)
    Temp = T0 + (GAMMA0 - 1.0) / RoM * (E_newton - E0);
  else
    Temp = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E_newton - E0) / RoM) - 1.0);

  rho = Press/(RoM*Temp);

  vec_out[0] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << Press << " ) => ( "
         << rho << " vs " << rho_newton << " )"
         << endl;
  }
  return vec_out;
}

VecDoub PressureTable::DensityFromPT(VecDoub_I vec) {
  double Press = vec[0]*P_scale;
  double rho;
  VecDoub vec_out(1, 0.0);

  double RoM = Lookup(Press, Z_newton,Zv_newton,C_newton, "RoM");
  rho = Press/(RoM*T_newton);

  vec_out[0] = rho-rho_newton;
  if (verbose_newton) {
    cout << "( " << Press << " ) => ( "
         << rho << " vs " << rho_newton << " )"
         << endl;
  }
  return vec_out;
}

VecDoub usrfunc_wrapper(VecDoub_I x) {
  return (PressureTable::ptr_for_newton)->YcAndDensityFromCP_new(x);
}

VecDoub usrfunc_wrapper2(VecDoub_I x) {
  return (PressureTable::ptr_for_newton)->DensityFromP_new(x);
}

VecDoub usrfunc_wrapper3(VecDoub_I x) {
  return (PressureTable::ptr_for_newton)->YcAndDensityFromCPT(x);
}

VecDoub usrfunc_wrapper4(VecDoub_I x) {
  return (PressureTable::ptr_for_newton)->DensityFromPT(x);
}

double PressureTable::ComputeNormalizedVariance(const double Zm, const double Zvar) {
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

double PressureTable::ComputeNormalizedProg(const double P, const double Zm, const double Sz, const double Yc) {
  double eps = 1.0e-5;
  double Yc_eq = Lookup(P, Zm, Sz, 1.0, "PROG");
  double C_out;
  if (Yc_eq >= eps) {
    C_out = Yc/Yc_eq;
  } else {
    C_out = 0.0;
  }
  return C_out;
}

void PressureTable::BCAST_pressure(int mpi_root, MPI_Comm mpi_comm) {
  // Broadcast Table4D attributes
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

}

void PressureTable::MPI_Load(string inputfilename) {
  int mpi_root = 0;
  // Only master process reads database file
  if (mpi_rank == mpi_root) Load(inputfilename);
  // Everyone waits for master to finish reading
  MPI_Barrier(mpi_comm);
  // Broadcast database 
  BCAST_pressure(mpi_root, mpi_comm);
  MPI_Barrier(mpi_comm);
}

