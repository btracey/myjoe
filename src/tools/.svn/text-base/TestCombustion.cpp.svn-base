/**
 * \brief Tool that check the combustion subroutines contained in Combustion.h.
 *
 * Perform interpolation of chemistry table and compare with table values;
 * Compute thermo properties of a mixture given the composition as function of temperature;
 * Compare thermo properties with flamelet solution for given Zmean, Zvar and Cmean or Chi.
 * \author Vincent Terrapon 
 * \date August 2009
 * \version 1.0
 */

#include "Combustion.h"
#include "Flamelet.h"
#include "ChemistryTableCartesianLinear.h"
#include "ChemistryTableCartesianCubic.h"
#include "ChemistryTableAdaptiveLinear.h"


template <class Chemtable>
class TestCombustion
{
private:
  Mixture   myMixture;            ///< Mixture containing the different species
  Chemtable myChemTable;          ///< Chemistry table used (rely on templates, should be incorporated in input file)
  string    CombustionModel;      ///< Type of combustion model (FPVA, STEADY)
  int       CombustionRegime;     ///< Type of combustion model (1: Steady Flamelet;  2:  FPVA)
  
public:
  /* Read thermo input file and load chemistry table */
  TestCombustion(ParamMap *myInputPtr)
  {
    string tablename;
    
    cout << endl << "***** TestCombustion initialization *****" << endl;
    
    // Open chemtable file and load chemistry table in memory
    tablename = myInputPtr->getStringParam("CHEMTABLE_FILE");    
    myChemTable.Load(tablename);
    CombustionModel = myChemTable.GetCombustionModel();
    if (CombustionModel == "FPVA")
      CombustionRegime = 2;
    else if (CombustionModel == "STEADY")
      CombustionRegime = 1;
    else if (CombustionModel == "ANALYTIC")
      CombustionRegime = 2;
    else
    {
      cerr << "### Wrong combustion model! should be FPVA, STEADY or ANALYTIC ###" << endl;
      throw(-1);
    }
    
    // Read thermodynamic input file (potentially containing more data than needed)
    myMixture.Load(CombustionRegime, myInputPtr, "YES");

    cout << endl << "***** TestCombustion initialized *****" << endl << endl;
  }
  
  /**********************************************************************************************/
  /* Clean memory */
  virtual ~TestCombustion()
  {
    myMixture.Unload();
    myChemTable.Unload();
    
    cout << endl << "***** TestCombustion finalized *****" << endl << endl;
  
    return;
  }
  
  /**********************************************************************************************/

  void TestChemtableInterpolation(ParamMap *myInputPtr)
  /* Check table interpolation: interpolated values for each direction for fixed other 
   * direction values and node values of table are output to file for comparison */
  {
    cout << "#####################################" << endl;
    cout << "Test interpolation of chemistry table" << endl;
    cout << "#####################################" << endl << endl;
    
    string chemtableType = myInputPtr->getStringParam("CHEMTABLE_TYPE");

    TestIndexSearch();
    TestNodalPointsInterpolation();
    
    if (CombustionModel == "ANALYTIC")
      TestAnalyticFunctionInterpolation();
    else
      cout << "Not performing 'TestAnalyticFunction()' since loaded table is not analytic table!" << endl;
    
    if ((chemtableType != "ADAPTIVE_LINEAR") && (chemtableType != "ADAPTIVE_CUBIC"))
    {
      double x1 = myInputPtr->getDoubleParam("1DPLOT_COORD_1");
      double x2 = myInputPtr->getDoubleParam("1DPLOT_COORD_2");
      double x3 = myInputPtr->getDoubleParam("1DPLOT_COORD_3");
      
      Test1DPlotInterpolation(x1, x2, x3);
    }
    else
      cout << "Not performing 'Test1DPlotInterpolation()' since loaded table is an adaptive table!" << endl;

    cout << endl;
    cout << "Chemistry interpolation test done!" << endl;
    cout << endl;
  }

/**********************************************************************************************/

  void TestIndexSearch()
  /* Check binary search */
  {
    cout << "Test indexing (binary search)" << endl;
    cout << "*****************************" << endl;
    cout << "Should not give any error message, except when value is less/greater than minimum/maximum value of table dimension;" << endl;
    cout << "this is due to clipping to these extreme values. Typical error message is:" << endl;
    cout << "=> ERROR!!! wrong index i from binary search at index ii : xi(ii) <= x < xi(ii+1)" << endl << endl;    
    
    double A1, A2, A3;    
    int nmax = 1000;
    
    int n1 = myChemTable.GetChemtableDimension1();
    int n2 = myChemTable.GetChemtableDimension2();
    int n3 = myChemTable.GetChemtableDimension3();
    
    double x1min = myChemTable.GetChemtableCoordinate1(1, 1, 1);
    double x1max = myChemTable.GetChemtableCoordinate1(n1, 1, 1);
    double x2min = myChemTable.GetChemtableCoordinate2(1, 1, 1);
    double x2max = myChemTable.GetChemtableCoordinate2(1, n2, 1);
    double x3min = myChemTable.GetChemtableCoordinate3(1, 1, 1);
    double x3max = myChemTable.GetChemtableCoordinate3(1, 1, n3);
   
    for (int i=0; i<nmax; i++)
    {
      A1 = x1min + double(i) * (x1max - x1min) / double(nmax-1);
      A2 = x2min + double(i) * (x2max - x2min) / double(nmax-1);
      A3 = x3min + double(i) * (x3max - x3min) / double(nmax-1);
      
      myChemTable.SetCoordinates(A1, A2, A3);
      
      int ii = myChemTable.GetInterpolationIndex(0);
      int jj = myChemTable.GetInterpolationIndex(1);
      int kk = myChemTable.GetInterpolationIndex(2);
      
      if ((A1 < myChemTable.GetChemtableCoordinate1(ii, jj, kk)) || (A1 >= myChemTable.GetChemtableCoordinate1(ii+1, jj, kk)))
        cout << "ERROR!!! wrong index 1 from binary search at index " << ii << " : " << myChemTable.GetChemtableCoordinate1(ii, jj, kk) << " <= " 
             << A1 << " < " << myChemTable.GetChemtableCoordinate1(ii+1, jj, kk) << endl; 
      if ((A2 < myChemTable.GetChemtableCoordinate2(ii, jj, kk)) || (A2 >= myChemTable.GetChemtableCoordinate2(ii, jj+1, kk)))
        cout << "ERROR!!! wrong index 2 from binary search at index " << jj << " : " << myChemTable.GetChemtableCoordinate2(ii, jj, kk) << " <= " 
             << A2 << " < " << myChemTable.GetChemtableCoordinate2(ii, jj+1, kk) << endl; 
      if ((A3 < myChemTable.GetChemtableCoordinate3(ii, jj, kk)) || (A3 >= myChemTable.GetChemtableCoordinate3(ii, jj, kk+1)))
        cout << "ERROR!!! wrong index 3 from binary search at index " << kk << " : " << myChemTable.GetChemtableCoordinate3(ii, jj, kk) << " <= " 
             << A3 << " < " << myChemTable.GetChemtableCoordinate3(ii, jj, kk+1) << endl; 
    }
    
    cout << endl;
  }

  /**********************************************************************************************/

  void TestNodalPointsInterpolation()
  /* Check binary search */
  {
    cout << "Test nodal points" << endl;
    cout << "*****************" << endl;
    cout << "Should not give any error message. Typical error message is:" << endl;
    cout << "=> ERROR!!! wrong nodal value at (x1, x2, x3): y_interp != ytable(ii,jj,kk)" << endl << endl;    
    
    double A1, A2, A3;    
    
    int n1 = myChemTable.GetChemtableDimension1();
    int n2 = myChemTable.GetChemtableDimension2();
    int n3 = myChemTable.GetChemtableDimension3();
       
    for (int i=1; i<n1+1; i++)
      for (int j=1; j<n2+1; j++)
        for (int k=1; k<n3+1; k++)
        {
          A1 = myChemTable.GetChemtableCoordinate1(i, j, k);
          A2 = myChemTable.GetChemtableCoordinate2(i, j, k);
          A3 = myChemTable.GetChemtableCoordinate3(i, j, k);
          
          for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
          {
            string var = myChemTable.GetChemtableVariables(l);
            double yinterp = myChemTable.Lookup(A1, A2, A3, var);
            double ytable  = myChemTable.GetChemtableValue(i, j, k, var);
            
            if (fabs(yinterp - ytable) > 1.0e-12)
              cout << "ERROR!!! wrong nodal value for variable " << var << " at (" << A1 << ", " << A2 << ", " << A3 << "): " 
                   << yinterp << " != " << ytable << endl; 
          }
        }
    
    cout << endl;
  }

  /**********************************************************************************************/

  void TestAnalyticFunctionInterpolation()
  /* Check binary search */
  {
    cout << "Test analytical functions" << endl;
    cout << "*************************" << endl;
    cout << "Compute error in L0-, L1- and L2-norms for a cubic function" << endl << endl;  
    
    double A1, A2, A3; 
    int nmax = 213;
    
    int n1 = myChemTable.GetChemtableDimension1();
    int n2 = myChemTable.GetChemtableDimension2();
    int n3 = myChemTable.GetChemtableDimension3();
    
    // bounds for x1, x2, x3 without the first cell to avoid errors at boundaries
    double x1min = myChemTable.GetChemtableCoordinate1(2, 1, 1);
    double x1max = myChemTable.GetChemtableCoordinate1(n1-1, 1, 1);
    double x2min = myChemTable.GetChemtableCoordinate2(1, 2, 1);
    double x2max = myChemTable.GetChemtableCoordinate2(1, n2-1, 1);
    double x3min = myChemTable.GetChemtableCoordinate3(1, 1, 2);
    double x3max = myChemTable.GetChemtableCoordinate3(1, 1, n3-1);
       
    
    double error0[4] = {0.0, 0.0, 0.0, 0.0}, error1[4] = {0.0, 0.0, 0.0, 0.0}, error2[4] = {0.0, 0.0, 0.0, 0.0};
    double yinterp, ytable, err;
    
    for (int i=0; i<=nmax; i++)
      for (int j=0; j<=nmax; j++)
        for (int k=0; k<=nmax; k++)
        {
          A1 = x1min + double(i) * (x1max - x1min) / double(nmax);
          A2 = x2min + double(j) * (x2max - x2min) / double(nmax);
          A3 = x3min + double(k) * (x3max - x3min) / double(nmax);
          
          yinterp = myChemTable.Lookup(A1, A2, A3, "Linear_Fcn");
          ytable  = linearFunction(A1, A2, A3);            
          err = yinterp - ytable;
          error0[0] = max(error0[0], fabs(err));
          error1[0] += fabs(err);
          error2[0] += err*err;
          
          yinterp = myChemTable.Lookup(A1, A2, A3, "Parabolic_Fcn");
          ytable  = parabolicFunction(A1, A2, A3);            
          err = yinterp - ytable;
          error0[1] = max(error0[1], fabs(err));
          error1[1] += fabs(err);
          error2[1] += err*err;
          
          yinterp = myChemTable.Lookup(A1, A2, A3, "Cubic_Fcn");
          ytable  = cubicFunction(A1, A2, A3);            
          err = yinterp - ytable;
          error0[2] = max(error0[2], fabs(err));
          error1[2] += fabs(err);
          error2[2] += err*err;
          
          yinterp = myChemTable.Lookup(A1, A2, A3, "General_Fcn");
          ytable  = otherAnalyticalFunction(A1, A2, A3);            
          err = yinterp - ytable;
          error0[3] = max(error0[3], fabs(err));
          error1[3] += fabs(err);
          error2[3] += err*err;
        }
    
    for (int i=0; i<4; i++)
    {
      error1[i] /= (nmax * nmax * nmax);
      error2[i] = sqrt(error2[i]) / (nmax * nmax * nmax);
    }
    
    cout << "Errors for linear function:               ";
    cout << "L0-norm = " << error0[0] << ",     L1-norm = " << error1[0] << ",      L2-norm = " << error2[0] << endl;
    cout << "Errors for parabolic function:            ";
    cout << "L0-norm = " << error0[1] << ",     L1-norm = " << error1[1] << ",      L2-norm = " << error2[1] << endl;
    cout << "Errors for cubic function:                ";
    cout << "L0-norm = " << error0[2] << ",     L1-norm = " << error1[2] << ",      L2-norm = " << error2[2] << endl;
    cout << "Errors for general function:              ";
    cout << "L0-norm = " << error0[3] << ",     L1-norm = " << error1[3] << ",      L2-norm = " << error2[3] << endl;
     
    cout << endl;
  }

  /**********************************************************************************************/

  void Test1DPlotInterpolation(double x1, double x2, double x3)
  {
    cout << "Test 1D plot interpolation" << endl;
    cout << "**************************" << endl;

    ofstream fout;
    
    int m = 2;
    
    int n1 = myChemTable.GetChemtableDimension1();
    int n2 = myChemTable.GetChemtableDimension2();
    int n3 = myChemTable.GetChemtableDimension3();

    double x1min = myChemTable.GetChemtableCoordinate1(1,  1, 1);
    double x1max = myChemTable.GetChemtableCoordinate1(n1, 1, 1);
    double x2min = myChemTable.GetChemtableCoordinate2(1, 1,  1);
    double x2max = myChemTable.GetChemtableCoordinate2(1, n2, 1);
    double x3min = myChemTable.GetChemtableCoordinate3(1, 1, 1 );
    double x3max = myChemTable.GetChemtableCoordinate3(1, 1, n3);
    
    int ii = myChemTable.GetChemtableIndex1(x1, x2, x3);
    int jj = myChemTable.GetChemtableIndex2(x1, x2, x3);
    int kk = myChemTable.GetChemtableIndex3(x1, x2, x3);
    
    double A1, A2, A3, A11, A22, A33;
    
    A11 = myChemTable.GetChemtableCoordinate1(ii, jj, kk);
    A22 = myChemTable.GetChemtableCoordinate2(ii, jj, kk);
    A33 = myChemTable.GetChemtableCoordinate3(ii, jj, kk);
    
    int nmax = 2500;
    
    
    
    // Open output file
    fout.open("Test1DPlotInterpolation.dat");
    if (fout.fail())
    {
      cerr << "### Cannot open output Test1DPlotInterpolation.dat file ###" << endl;
      throw(-1);
    }
    
    fout << "TITLE=\"1D comparison between interpolated data and nodal values\"" << endl;
    fout << "VARIABLES=\"Zmm\" \"Zvv\" \"Cmm\" ";
    for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      fout << "\"" << myChemTable.GetChemtableVariables(l) << "\" ";
    fout << endl;
    fout.setf(ios::fixed, ios::floatfield);
    fout.precision(15);
    
    // First direction, interpolated values
    fout << "ZONE" << endl;
    for (int i=0; i<=nmax; i++)
    {
      A1 = x1min + double(i) * (x1max - x1min) / double(nmax);

      fout << A1 << " " << A22 << " " << A33 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.Lookup(A1, A22, A33, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
    
    // Second direction, interpolated values
    fout << "ZONE" << endl;
    for (int j=0; j<=nmax; j++)
    {
      A2 = x2min + double(j) * (x2max - x2min) / double(nmax);

      fout << A11 << " " << A2 << " " << A33 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.Lookup(A11, A2, A33, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
    
    // Third direction, interpolated values
    fout << "ZONE" << endl;
    for (int k=0; k<=nmax; k++)
    {
      A3 = x3min + double(k) * (x3max - x3min) / double(nmax);

      fout << A11 << " " << A22 << " " << A3 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.Lookup(A11, A22, A3, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
        
    
    // First direction, nodal values
    fout << "ZONE" << endl;
    for (int i=1; i<=n1; i++)
    {
      A1 = myChemTable.GetChemtableCoordinate1(i, jj, kk);
      
      fout << A1 << " " << A22 << " " << A33 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.GetChemtableValue(i, jj, kk, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
    
    // Second direction, nodal values
    fout << "ZONE" << endl;
    for (int j=1; j<=n2; j++)
    {
      A2 = myChemTable.GetChemtableCoordinate2(ii, j, kk);
      
      fout << A11 << " " << A2 << " " << A33 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.GetChemtableValue(ii, j, kk, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
    
    // Third direction, nodal values
    fout << "ZONE" << endl;
    for (int k=1; k<=n3; k++)
    {
      A3 = myChemTable.GetChemtableCoordinate3(ii, jj, k);
      
      fout << A11 << " " << A22 << " " << A3 << " ";
      for (int l=0; l<myChemTable.GetChemtableDimension4(); l++)
      {
        string var = myChemTable.GetChemtableVariables(l);
        double yinterp = myChemTable.GetChemtableValue(ii, jj, k, var);
        fout << " " << yinterp << " ";
      }
      fout << endl;
    }
  
    fout.close();
  }

  /**********************************************************************************************/ 
  
  void TestChemtableMixtureThermoInterpolation(int type, ParamMap *myInputPtr)
  /* Compute thermoproperties as function of Z for zero variance and a progress variable from the flamelet files */
  {
    Flamelet        *myFlamelet;
    vector <string> VarProg;
    vector <double> WeightProg;
    vector <string> FlameletList; 
    int             nProgVar, nfiles, nkspecies, nZPoints;
    string          speciesname, ProgName, variablename, file, buffer_str;
    ofstream        fout;
    ifstream        fin;
    Param           *var=NULL;
    double          weight;
    double          chi_st, Z, C, chi, ent, Yk, W, temperature, ener;
    double          temp, E, rho, RoM, enthalpy, press, cp, muLam, LambdaOverCp, gamma, rhoDkZ, sos, CSource;
    
    
    if (CombustionModel == "ANALYTIC")
      return;
    
    cout << endl;
    cout << "######################" << endl;
    cout << "Test thermodyanmics " << type << endl;
    cout << "######################" << endl << endl;
    
    
    // Get some parameters and input values
    if (CombustionRegime == 2)
    {
      var = myInputPtr->getParam("VAR_PROG");
      nProgVar = var->getSize() - 1;
      if (nProgVar <= 0)
      {
        cerr << "### Variables for progress variable are missing in input file! ###" << endl;
        throw(-1);
      }
      for (int nP=0; nP<nProgVar; nP++)
      {
        speciesname = var->getString(nP + 1);
        VarProg.push_back(speciesname);
      }
      var = myInputPtr->getParam("WEIGHT_PROG");
      if (nProgVar != (var->getSize() - 1))
      {
        cout << "Weights for progress variable are missing in input file! Weights are automatically set to 1.0" << endl;
        for (int nP=0; nP<nProgVar; nP++)
          WeightProg.push_back(1.0);
      }
      for (int nP=0; nP<nProgVar; nP++)
      {
        weight = var->getDouble(nP + 1);
        WeightProg.push_back(weight);
      }
    }
    else if (CombustionRegime == 1)
    {
      VarProg.push_back("chi");
      WeightProg.push_back(1.0);
    }
    else
    {
      cerr << "### Combustion regime " << CombustionRegime << " not implemented! ###" << endl;
      throw(-1);
    }
    
    nkspecies = myMixture.GetNkspecies();
    
    // Open flamelet files
    string FlameletListFilename = myInputPtr->getStringParam("FLAMELET_LIST_FILENAME", "FlameletList.txt");
    char filename[kStrLong];
    size_t dum = FlameletListFilename.copy(&filename[0], kStrLong);
    dum = FlameletListFilename.size();
    if (dum < kStrLong)
      strcpy(&filename[dum], "\0");
    fin.open(filename);
    if (fin.fail())
    {
      cerr << "### Cannot open flamelet list input file " << filename << " ! ###" << endl;
      throw(-1);
    }
    while(getline(fin, buffer_str))
    {
      istringstream buf(buffer_str);
      buf >> file;
      if (buffer_str.empty())
        break;
      FlameletList.push_back(file);
    }
    fin.close();
    nfiles = FlameletList.size();
      
    // Open output file
    if (type == 1)
      fout.open("TestThermoPropertiesInterpolation.dat");
    else
      fout.open("TestThermoProperties.dat");
    if (fout.fail())
    {
      cerr << "### Cannot open output flamelet file ###" << endl;
      throw(-1);
    }
    fout << "TITLE=\"Thermo properties given Z at cst chi_st\"" << endl;
    fout << "VARIABLES=\"Z\" \"Chi\" \"Rgas\" \"W\" \"T\" \"E\" \"Enthalpy\" \"Cp\" \"rho\" \"P\" ";
    fout << "\"mu\" \"LamOCp\" \"gam\" \"rhoDk\" \"sos\" \"Csource\" "; 
    for (int j=0; j<nkspecies; j++)
      fout << "\"" << myMixture.GetMixSpecies(j) << "\" ";
    fout << endl;
  
    
    // Loop over all flamelets
    for (int nf=0; nf<nfiles; nf++)
    {    
      myFlamelet = new Flamelet;
      myFlamelet->LoadFlamelet(CombustionRegime, myInputPtr, FlameletList[nf], VarProg, WeightProg, "NO");
      chi_st = myFlamelet->GetFlameletChiSt();
      nZPoints = myFlamelet->GetFlameletSize();
      
      fout << "ZONE, T=\"Thermo properties case at chi_st=" << chi_st << "\"" << endl;
      
      // Loop over Z within flamelet
      for (int nZ=0; nZ<nZPoints; nZ++)
      {
        // Get flamelet values
        Z = myFlamelet->GetFlameletZ(nZ);
        ent = myFlamelet->GetFlameletVariable(nZ, "TotalEnthalpy");
        chi = myFlamelet->GetFlameletVariable(nZ, "chi");
        rho = myFlamelet->GetFlameletVariable(nZ, "density");
        temperature = myFlamelet->GetFlameletVariable(nZ, "temperature");
        W = myFlamelet->GetFlameletVariable(nZ, "W");
        ener = ent - kR / W * temperature;
        if (CombustionRegime == 2)
          C = myFlamelet->GetFlameletVariable(nZ, "PROG");
        
        if (myChemTable.GetChemtableType() != "COEFF")
        {        
          // Compute Yk from table through interpolation ...
          if (type == 1)
          {
            if (CombustionRegime == 2)
              myMixture.GetSpeciesMassFraction(CSource, Z, 0.0, C, myChemTable);
            else
            {
              myMixture.GetSpeciesMassFraction(Z, 0.0, chi, myChemTable);
              CSource = 0.0;
            }
          }
          // ... or get Yk directly from flamelet solution (by pass interpolation issues)
          else
          {
            for (int j=0; j<nkspecies; j++)
            {
              speciesname = myMixture.GetMixSpecies(j);
              Yk = myFlamelet->GetFlameletVariable(nZ, speciesname);
              myMixture.SetSpeciesMassFraction(Yk, speciesname);
            }
            myMixture.ComputeMixM();
          }  
          
          // Compute thermo properties given Yk
          CSource *= rho;
          RoM = myMixture.GetMixR_o_M();
          temp = myMixture.ComputeMixTemperature(ener);
          E = myMixture.ComputeMixEnergy(temp);
          //temp = T;
          enthalpy = myMixture.ComputeMixEnthalpy(temp);
          press = rho * RoM * temp;
          cp = myMixture.ComputeMixCp(temp);
          myMixture.ComputeMixLambdaMu(temp);
          muLam = myMixture.GetMixMul();
          LambdaOverCp = myMixture.GetMixLambda() / cp;
          gamma = cp / (cp - RoM);
          rhoDkZ = LambdaOverCp / myMixture.GetMixLewis();
          sos = sqrt(gamma*press/rho);
        }
        else
        {
          double T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC;
          myChemTable.LookupCoeff(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CSource, Z, 0.0, C);
          CSource *= rho;
          if (AGAMMA == 0.0)
            temp = T0 + (GAMMA0 - 1.0) / RoM * (ener - E0);
          else
            temp = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (ener - E0) / RoM) - 1.0);
          if (AGAMMA == 0.0)
            E = E0 + RoM / (GAMMA0 - 1.0) * (temp - T0);
          else
            E = E0 + RoM / AGAMMA * log(1.0 + AGAMMA * (temp - T0) / (GAMMA0 - 1.0));
          enthalpy = E + RoM * temp;
          press = rho * RoM * temp;
          muLam = MU0 + AMU * (temp - T0);
          LambdaOverCp = LOC0 + ALOC * (temp - T0);
          gamma = GAMMA0 + AGAMMA * (temp - T0);
          cp = gamma * RoM / (gamma - 1.0);
          sos = sqrt(gamma * RoM * temp);
          rhoDkZ = LambdaOverCp / 1.0;
        }
      
        fout << Z << " " << chi << " " << RoM << " " << kR/RoM << " " << temp << " " << E << " " << enthalpy << " " << cp << " " << rho << " " << press << " " 
             << muLam << " " << LambdaOverCp << " " << gamma << " " << rhoDkZ << " " << sos << " " << CSource;
        for (int j=0; j<nkspecies; j++)
        {
          speciesname = myMixture.GetMixSpecies(j);
          fout << " " << myMixture.GetMixYk(speciesname);
        }
        fout << endl;
      }
      if (myFlamelet != NULL)  delete myFlamelet;
      myFlamelet = NULL;
    }
    fout << endl;
  
  }
};

/**********************************************************************************************/

int main(int argc, char *argv[])
{
  char InputFileName[kStrLong];
  sprintf(InputFileName, "TestCombustion.in");
  for (int i=1; i<argc; i++)
    strcpy(InputFileName, argv[i]);

  ParamMap myInput(InputFileName);
  string chemtableType = myInput.getStringParam("CHEMTABLE_TYPE");
  
  if (chemtableType == "CARTESIAN_LINEAR")
  {
    TestCombustion<ChemtableCartesianLinear> testcomb(&myInput);
    
    testcomb.TestChemtableInterpolation(&myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(1, &myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(2, &myInput);
  }
  else if (chemtableType == "CARTESIAN_CUBIC")
  {
    TestCombustion<ChemtableCartesianCubic> testcomb(&myInput);
    
    testcomb.TestChemtableInterpolation(&myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(1, &myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(2, &myInput);
  }
  else if (chemtableType == "ADAPTIVE_LINEAR")
  {
    TestCombustion<ChemtableAdaptiveLinear> testcomb(&myInput);
    
    testcomb.TestChemtableInterpolation(&myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(1, &myInput);
    testcomb.TestChemtableMixtureThermoInterpolation(2, &myInput);
  }
  else 
  {
    cerr << "### Wrong chemistry table type! Should be CARTESIAN_LINEAR, CARTESIAN_CUBIC. ###" << endl;
    throw(-1);
  }
  
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
