/**
 * \brief Flamelet class containing flamelet objects from FlameMaster
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */

#ifndef FLAMELET_H
#define FLAMELET_H

//#define OLD_BETA

#include "Combustion.h"
#include <fstream>
using std::ifstream;

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Class Flamelet objects containing solutions of flamelets from FlameMaster.
 * 
 *  Load flamelet solution from FlameMaster file into a flamelet object and calculate progress variable and its source term.
 */
class Flamelet
{
private:
  
  double          Chi_st;                 ///< Stoichiometric scalar dissipation rate of flamelet.
  double          Tmax;                   ///< Maximum temperature of flamelet.
  double          Pressure;               ///< Pressure at which flamelet was computed.
  double          Toxi, Tfuel;            ///< Temperature boundary condition with which flamelet was computed.
  int             nSpecies;               ///< Number of species in flamelet.
  int             nVariables;             ///< Number of variables loaded.
  vector <string> VariableNames;          ///< Vector containing the name of all variables loaded.
  int             nPoints;                ///< Number of points in Z (also number of data for each variable).
  double          **myData;               ///< Data stored as 2D array of size Npoints x NVariables (dim 1: Z, dim 2: variables).
  double          *pdf;                   ///< Beta PDF used for convolution parametrized by a mean and variance.
  double          BetaMean;               ///< Mean of Beta PDF (used to avoid recomputing the pdf many times.
  double          BetaVar;                ///< Variance of Beta PDF (used to avoid recomputing the pdf many times.
  Mixture         *myMixture;             ///< Mixture for thermodynamics (array, 1 mixture for each Z).
  int             CombustionRegime;       ///< Combustion regime (1->steady flamelet, 2->fpva)
  
public:
  
  // Constructor
  Flamelet() 
  {
    myData = NULL;
    pdf = NULL;
    Chi_st = 0.0;
    Tmax = 0.0;
    Pressure = 0.0;
    nSpecies = 0;
    nVariables = 0;
    nPoints = 0;
    BetaMean = -1.0;
    BetaVar = -1.0;
  }
  
  // Destructor
  ~Flamelet()
  {
    if (myData != NULL) freeMem2D(myData, 0, nPoints-1, 0, nVariables-1); myData = NULL;
    VariableNames.clear();
    if (pdf != NULL) delete [] pdf; pdf = NULL;
    for (int k=0; k<nPoints; k++)
      myMixture[k].Unload();
    if (myMixture != NULL) delete [] myMixture; myMixture = NULL;
  }
  
  //Accessors
  int GetFlameletSize() {return nPoints;}
  int GetFlameletNSpecies() {return nSpecies;}
  int GetFlameletNVar() {return nVariables;}
  double GetFlameletChiSt() {return Chi_st;}
  double GetFlameletTmax() {return Tmax;}
  double GetFlameletPressure() {return Pressure;}
  double GetFlameletToxi() {return Toxi;}
  double GetFlameletTfuel() {return Tfuel;}
  double GetFlameletZ(int i) {return myData[i][0];}
  int GetFlameletVariableIndex(string var_str)
  {
    int j = -1;
    for (int var=0; var<VariableNames.size(); var++)
      if (!(var_str.compare(VariableNames.at(var))))
        j = var;
    if ((j < 0) || (j > nVariables))
    {
      cerr << "Variable " << var_str << " not found in loaded flamelet data!" << endl;
      throw(-1);
    }
    return j;
  }
  double GetFlameletVariable(int i, int j) 
  {
    return myData[i][j];
  }
  double GetFlameletVariable(int i, string var_str)
  {
    int j = GetFlameletVariableIndex(var_str);
    return myData[i][j];
  }
  double GetFlameletVariable(double Z, int j)
  {
    int i;
    int J = GetFlameletVariableIndex("Z");
    double w;
    if (Z < 0.0)
    {
      i = 0;
      w = 1.0;
    }
    else if (Z > 1.0)
    {
      i = nPoints - 2;
      w = 0.0;
    } 
    else
    {
      for (j=0; j<nPoints - 1; j++)
      {
        if (Z < myData[j+1][J])
        {
          i = j;
          break;
        }
      }
      w = (myData[i+1][J] - Z) / (myData[i+1][J] - myData[i][J]);
    }
    return w * myData[i][j] + (1.0 - w) * myData[i+1][j];
  }
  double GetFlameletVariable(double Z, string var_str)
  {
    int j = GetFlameletVariableIndex(var_str);
    return GetFlameletVariable(Z, j);
  }
  double GetFlameletVariableMax(int j)
  {
    if ((j < 0) || (j > nVariables))
    {
      cerr << "Index " << j << " out of bound in GetFlameletVariableMax" << endl;
      throw(-1);
    }
    double varMax = myData[0][j];
    for (int i=1; i<nPoints; i++)
      if (myData[i][j] > varMax)
        varMax = myData[i][j];
    return varMax;
  }
  double GetFlameletVariableMax(string var_str)
  {
    int j = GetFlameletVariableIndex(var_str);
    return GetFlameletVariableMax(j);
  }
  
  /** \brief Load data into flamelet object from FlameMaster solution file. 
   *  The input parameters are:
   *  filename: the name of the flamelet solution file
   *  VarProg: a vector of variables to build the progress variable
   *  WeightProg: a vector containing the corresponding weights to build the progress variable
   *  If steady flamelet is used, chi is considered the progress variable and the source term is set to 0.
   */
  void LoadFlamelet(int CombustionRegime, ParamMap *myInputPtr, string filename, const vector <string>& VarProg, const vector <double>& WeightProg, const string diag_output = "NO",
                    const string DefinitionProg = "USUAL", const double Zst = 1.0)
  {
    char      file[kStrLong];
    ifstream  fin;
    string    buffer_str, dummy, head;
    size_t    dum;
    int       nLines, nRest;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Check if VarProg and WeightProg have the same size
    if (VarProg.size() != WeightProg.size())
    {
      cerr << "### The number of progress variables does not match the number of weights! ###" << endl;
      throw(-1);
    }

    // Open flamelet file
    dum = filename.copy(&file[0], kStrLong);
    dum = filename.size();
    if (dum < kStrLong)
      strcpy(&file[dum], "\0");
    cout << "Flamelet " << filename;
    fin.open(file);
    if (fin.fail())
    {
      cerr << "### Cannot open input flamelet file " << filename << " ###" << endl;
      throw(-1);
    }
    
    // Read header part of the file to get general informations
    int Tfuelread = 0;
    while(getline(fin,buffer_str))
    {
      istringstream buf(buffer_str);
      buf >> head;
      if (!(head.compare("chi_st")))
      {
        buf >> dummy >> Chi_st;
      }
      if (!(head.compare("Tmax")))
      {
        buf >> dummy >> Tmax;
      }
      if (!(head.compare("pressure")))
      {
        buf >> dummy >> Pressure;
        Pressure *= 100000.0;                      // convert from bar to Pa
      }
      if (!(head.compare("Temperature")) && (Tfuelread != 1))          // first comes fuel temperature
      {
        buf >> dummy >> Tfuel;
        Tfuelread = 1;
      }
      if (!(head.compare("Temperature")) && (Tfuelread != 0))          // then comes oxidizer temperature
      {
        buf >> dummy >> Toxi;
      }
      if (!(head.compare("numOfSpecies")))
      {
        buf >> dummy >> nSpecies;
      }
      if (!(head.compare("gridPoints")))
      {
        buf >> dummy >> nPoints;
      }
      if (!(head.compare("body")))
      {
        break;
      }
    }
    
    // Calculate number of data lines -1 for each variable (5 data per line) and the number of data on last line
    nLines = nPoints / 5;
    nRest = nPoints % 5;
    
    // Read the file to get the variable names
    while(getline(fin,buffer_str))
    {
      istringstream buf(buffer_str);
      buf >> head;
      if (!(head.compare("trailer")))
        break;
      VariableNames.push_back(head);
      // Skip data
      for (int nL=0; nL<nLines; nL++)
      {
        getline(fin,buffer_str);
      }
      if (nRest != 0)
      {
        getline(fin,buffer_str);
      }
    }
    fin.close();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    nVariables = VariableNames.size() + 14;  // adding progress variable, its source term, R0, T0, R0*T0, rho0, energy0, mu0, lamOcp0, T1, R0*T1, energy1, mu1, lamOcp1
    nVariables = nVariables + 10;
    // Allocate memory for myData
    getMem2D(&myData, 0, nPoints-1, 0, nVariables-1, "LoadFlamelet: myData", true);
    pdf = new double[nPoints];
    
    myMixture = new Mixture[nPoints];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Reopen file to read variables
    fin.open(file);
    while(getline(fin,buffer_str))
    {
      istringstream buf(buffer_str);
      buf >> head;
      if (!(head.compare("body")))
      {
        break;
      }
    }
    // Read data
    int j=0;  // variable coordinate
    while(getline(fin,buffer_str))
    {
      istringstream buf(buffer_str);
      buf >> head;
      if (!(head.compare("trailer")))
        break;
      int i=0;
      for (int nL=0; nL<nLines; nL++)
      {
        getline(fin,buffer_str);
        istringstream buf(buffer_str);
        for (int nD=0; nD<5; nD++)
        {
          buf >> myData[i][j];
//            cout << i << "   " << myData[i][j] << endl;
          i++;
        }
      }
      if (nRest != 0)
      {
        getline(fin,buffer_str);
        istringstream buf(buffer_str);
        for (int nD=0; nD<nRest; nD++)
        {
          buf >> myData[i][j];
//            cout << i << "   " << myData[i][j] << endl;
          i++;
        }
      }
      j++;
    }
    fin.close();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Divide Species reaction rates and Heat Release by density
    
    int jrho = GetFlameletVariableIndex("density");
    // Species reaction rates
    string auxstring;
    for (int j=0; j<VariableNames.size(); j++) {
      auxstring = VariableNames[j].substr(0,8);
      if (!auxstring.compare("ProdRate")) {
        for (int k=0; k<nPoints; k++)
          myData[k][j] /= myData[k][jrho];
      }
    }
    // Heat Release
    int jHR  = GetFlameletVariableIndex("HeatRelease");
    for (int k=0; k<nPoints; k++)
      myData[k][jHR] /= myData[k][jrho];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Rename variables for species (i.e., remove 'massfraction-') for compatibility with main code
    for (int i=0; i<VariableNames.size(); i++)
    {
      size_t pos;
      pos = VariableNames[i].find("massfraction-");
      if ((int)pos != -1)
        VariableNames[i].erase(pos,13);
    }
//    for (int i=0; i<VariableNames.size(); i++)
//      cout << VariableNames[i] << endl;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Load mixture at each Z points
    
    for (int k=0; k<nPoints; k++)
    {
      if ((k == 0) && (diag_output == "YES"))
        myMixture[k].Load(CombustionRegime, myInputPtr, "YES");
      else
        myMixture[k].Load(CombustionRegime, myInputPtr, "NO");
      int nkspecies = myMixture[k].GetNkspecies();
//      if (nSpecies != nkspecies)
//      {
//        cerr << "### The number of species loaded in the mixture from input file, " << nkspecies 
//             << ", does not match the number of species in the flamelet solution, " << nSpecies << " ! ###" << endl;
//        throw(-1);
//      }
      for (int j=0; j<nkspecies; j++)
      {
        string speciesname = myMixture[k].GetMixSpecies(j);
        double Yk = GetFlameletVariable(k, speciesname);
        myMixture[k].SetSpeciesMassFraction(Yk, speciesname);
      }
      myMixture[k].ComputeMixM();
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Compute and add useful quantities: R0, T0, R0*T0, rho0, e0, mu0, lamOcp0
    int jH    = GetFlameletVariableIndex("TotalEnthalpy");
    int jT    = GetFlameletVariableIndex("temperature");
    int jR0   = VariableNames.size();
    int jT0   = jR0+1;
    int jRT0  = jT0+1;
    int jrho0 = jRT0+1;
    int jE0   = jrho0+1;
    int jmu0  = jE0+1;
    int jlam0 = jmu0+1;
  
    for (int k=0; k<nPoints; k++)
    {
      // gas constant R0
      myData[k][jR0] = myMixture[k].GetMixR_o_M();
      
      // temperature T0
      myData[k][jT0] = myMixture[k].ComputeMixTemperature_H(myData[k][jH], myData[k][jT]);  
      
      // gas constant * temperature R0*T0
      myData[k][jRT0] = myData[k][jR0] * myData[k][jT0];
      
      // density = P / RT
      myData[k][jrho0] = Pressure / myData[k][jRT0];
      
      // energy
      myData[k][jE0] = myData[k][jH] - myData[k][jRT0];
      
      myMixture[k].ComputeMixLambdaMu(myData[k][jT0]);      
      // viscosity
      myData[k][jmu0] = myMixture[k].GetMixMul();
      
      // thermal diffusion lambda0 / cp0
      myData[k][jlam0] = myMixture[k].GetMixLambda() / myMixture[k].ComputeMixCp(myData[k][jT0]);
    }
    
    VariableNames.push_back("R0");
    VariableNames.push_back("T0");
    VariableNames.push_back("RT0");
    VariableNames.push_back("rho0");
    VariableNames.push_back("E0");
    VariableNames.push_back("mu0");
    VariableNames.push_back("lamOcp0");      
    
    // Save storage for energy perturbution plus and minus (T, R*T, E, mu, lamOcp)
    VariableNames.push_back("Tpl");
    VariableNames.push_back("RTpl");
    VariableNames.push_back("Epl");
    VariableNames.push_back("mupl");
    VariableNames.push_back("lamOcppl"); 
    VariableNames.push_back("Tmi");
    VariableNames.push_back("RTmi");
    VariableNames.push_back("Emi");
    VariableNames.push_back("mumi");
    VariableNames.push_back("lamOcpmi"); 
    
    // Save storage for working variables T1, R0*T1, E1, mu1, lamOcp1
    VariableNames.push_back("T1");
    VariableNames.push_back("RT1");
    VariableNames.push_back("E1");
    VariableNames.push_back("mu1");
    VariableNames.push_back("lamOcp1"); 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    // Compute progress variable and source term of progress variable
    for (int i=0; i<VarProg.size(); i++)
    {
      // Progress variable
      int jV = -1;
      for (int j=0; j<VariableNames.size(); j++)
      {
        if (!VariableNames[j].compare(VarProg[i]))
          jV = j;
      }      
      if (jV == -1)
      {
        cerr << "### Error: the progress variable " << VarProg[i] << " was not found! ###" << endl;
        throw(-1);
      }
      for (int k=0; k<nPoints; k++)
        myData[k][nVariables-2] += WeightProg[i] * myData[k][jV];

      if (DefinitionProg == "GENERALIZED")
      {
        //Progress variable is set as C|Zst to be independent of mixture fraction
        double C_at_Zst = GetFlameletVariable(Zst, "PROG");
        for (int k=0; k<nPoints; k++)
          myData[k][nVariables-2] = C_at_Zst; 
      }
      
      // Source term
      jV = -1;
      string var = VarProg[i];
      if (!var.compare("chi"))
        for (int k=0; k<nPoints; k++)
          myData[k][nVariables-1] += 0.0;
      else
      {      
        if (!var.compare("temperature"))
          var = "HeatRelease";
        else
          var.insert(0,"ProdRate");
        
        for (int j=0; j<VariableNames.size(); j++)
        {
          if (!VariableNames[j].compare(var))
            jV = j;
        }      
        if (jV == -1)
        {
          cerr << "### Error: the source term for the progress variable " << var << " was not found! ###" << endl;
          throw(-1);
        }
        for (int k=0; k<nPoints; k++)
          myData[k][nVariables-1] += WeightProg[i] * myData[k][jV];
     }
    }
    
    VariableNames.push_back("PROG");
    VariableNames.push_back("SRC_PROG");
    
     cout << "    nPoints: " << nPoints << "    nSpecies: " << nSpecies << "    nVariables: " << nVariables 
         << "     Chi_st: " << Chi_st << "     Tmax: " << Tmax << "    Loaded!" << endl;
    
    return;
  }
  
  /***********************************************************************************************************/  

  /** \brief Given a perturbation in energy Delta_e, compute E1, T1, RT1, mu1 and lamOcp1 */
  void ComputeEnergyPerturbation(double rho_0_bar_Delta_e_tild)
  {
    int jE0    = GetFlameletVariableIndex("E0");
    int jR0    = GetFlameletVariableIndex("R0");
    int jT0    = GetFlameletVariableIndex("T0");
    int jrho0  = GetFlameletVariableIndex("rho0");
    int jEpl   = GetFlameletVariableIndex("Epl");
    int jTpl   = GetFlameletVariableIndex("Tpl");
    int jRTpl  = GetFlameletVariableIndex("RTpl");
    int jmupl  = GetFlameletVariableIndex("mupl");
    int jlampl = GetFlameletVariableIndex("lamOcppl");
    int jEmi   = GetFlameletVariableIndex("Emi");
    int jTmi   = GetFlameletVariableIndex("Tmi");
    int jRTmi  = GetFlameletVariableIndex("RTmi");
    int jmumi  = GetFlameletVariableIndex("mumi");
    int jlammi = GetFlameletVariableIndex("lamOcpmi");
    
    // Positive perturbation
    for (int k=0; k<nPoints; k++)
    {
      myData[k][jEpl]   = myData[k][jE0] + rho_0_bar_Delta_e_tild / myData[k][jrho0];
      myData[k][jTpl]   = myMixture[k].ComputeMixTemperature(myData[k][jEpl], myData[k][jT0]);
      myData[k][jRTpl]  = myData[k][jTpl] * myData[k][jR0];
      
      myMixture[k].ComputeMixLambdaMu(myData[k][jTpl]);
      myData[k][jmupl]  = myMixture[k].GetMixMul();
      myData[k][jlampl] = myMixture[k].GetMixLambda() / myMixture[k].ComputeMixCp(myData[k][jTpl]);
    }
    
    // Negative perturbation
    for (int k=0; k<nPoints; k++)
    {
      myData[k][jEmi]   = myData[k][jE0] - rho_0_bar_Delta_e_tild / myData[k][jrho0];
      myData[k][jTmi]   = myMixture[k].ComputeMixTemperature(myData[k][jEmi], myData[k][jT0]);
      myData[k][jRTmi]  = myData[k][jTmi] * myData[k][jR0];
      
      myMixture[k].ComputeMixLambdaMu(myData[k][jTmi]);
      myData[k][jmumi]  = myMixture[k].GetMixMul();
      myData[k][jlammi] = myMixture[k].GetMixLambda() / myMixture[k].ComputeMixCp(myData[k][jTmi]);
    }
  }
  
  /***********************************************************************************************************/  

  /** \brief Given a perturbation in energy Delta_e, compute E1, T1, RT1, mu1 and lamOcp1 */
  void ComputeEnergyPerturbation(double rho_0_bar, double Delta_e_tild)
  {
    int jE0   = GetFlameletVariableIndex("E0");
    int jR0   = GetFlameletVariableIndex("R0");
    int jT0   = GetFlameletVariableIndex("T0");
    int jrho0 = GetFlameletVariableIndex("rho0");
    int jE1   = GetFlameletVariableIndex("E1");
    int jT1   = GetFlameletVariableIndex("T1");
    int jRT1  = GetFlameletVariableIndex("RT1");
    int jmu1  = GetFlameletVariableIndex("mu1");
    int jlam1 = GetFlameletVariableIndex("lamOcp1");
    
    for (int k=0; k<nPoints; k++)
    {
      myData[k][jE1]   = myData[k][jE0] + rho_0_bar / myData[k][jrho0] * Delta_e_tild;
      myData[k][jT1]   = myMixture[k].ComputeMixTemperature(myData[k][jE1], myData[k][jT0]);
      myData[k][jRT1]  = myData[k][jT1] * myData[k][jR0];
      
      myMixture[k].ComputeMixLambdaMu(myData[k][jT1]);
      myData[k][jmu1]  = myMixture[k].GetMixMul();
      myData[k][jlam1] = myMixture[k].GetMixLambda() / myMixture[k].ComputeMixCp(myData[k][jT1]);
    }
  }
  
  /***********************************************************************************************************/  
 
  /** \brief Convolute the flamelet with a Beta-PDF given PDF parameters mean and variance for given variable */
  void ComputeBetaPDF(double mean, double var)
  {
    // Beta distribution is parametrized by a and b so that
    // pdf = Gamma(a+b)/(Gamma(a)*Gamma(b)) * x^(a-1) * (1-x)^(b-1)     where
    // mean = a / (a+b)
    // variance = (a*b) / ((a+b)^2 * (a+b+1))
    // which gives
    // a = mean * ( (mean*(1-mean))/variance - 1)
    // b = a / mean - a
    // The pdf is computed as e^[log(pdf)] to increase accuracy
    
    // Recompute the Beta PDF only if mean and var have changed
    if ((mean != BetaMean) || (var != BetaVar))
    {
      // reset BetaMean and BetaVar
      BetaMean = mean;
      BetaVar = var;
      
      // Get index of Z variable
      int J = GetFlameletVariableIndex("Z");
      
      // Initialize
      for (int k=0; k<nPoints; k++) {
        pdf[k] = 0.0;
      }
      
      // Zero mean: delta pdf at Z=0
      if (mean <= 1.0e-10)
      {
        pdf[0] = 1.0;
        return;
      }
      
      // Max mean: delta pdf at Z=1
      if (mean >= (1.0 - 1.0e-10))
      {
        pdf[nPoints-1] = 1.0;
        return;
      }
      
      // Zero variance: delta pdf at Z=mean
      if (var <= 1.0e-6) // RV change tolerance
      {
        double w;
        int k;
        for (int i=0; i<nPoints-1; i++)
        {
          if (mean < myData[i+1][J])
          {
            k = i;
            break;
          }
        }
        w = (myData[k+1][J] - mean) / (myData[k+1][J] - myData[k][J]);
        pdf[k] = w;
        pdf[k+1] = 1.0 - w;       
        return;
      }
      
      // Impossible cases (i.e., var > mean*(1-mean)): two delta pdf at Z=0 and Z=1
      if (var >= mean * (1.0 - mean))
      {
        pdf[0] = 1.0 - mean;
        pdf[nPoints-1] = mean;
        return;
      }
      
      // Otherwise
      double a = mean * (mean * (1.0 - mean) / var - 1.0);
      double b = a / mean - a;
      double factor = lgamma(a+b) - lgamma(a) - lgamma(b);
      double dz, tmp;
      
      // Left BC: explicit integration
#ifdef OLD_BETA
      dz = 0.5 * (myData[1][J] - myData[0][J]);
      tmp = a * log(dz) + factor;
      pdf[0] = exp(tmp) / a;
#else
      // RV change beta computation
      double z1 = myData[0][J];
      double z2 = 0.5 * (myData[1][J] + myData[0][J]); 
      pdf[0] = betai(a,b,z2)-betai(a,b,z1);
#endif

      // Right BC: explicit integration
#ifdef OLD_BETA
      dz = 0.5 * (myData[nPoints-1][J] - myData[nPoints-2][J]);
      tmp = b * log(dz) + factor;
      pdf[nPoints-1] = exp(tmp) / b;
#else
      // RV change beta computation
      z1 = 0.5 * (myData[nPoints-1][J] + myData[nPoints-2][J]);
      z2 = myData[nPoints-1][J];
      pdf[nPoints-1] = betai(a,b,z2)-betai(a,b,z1); 
#endif
      // Other points
      for (int k=1; k<nPoints-1; k++)
      {
#ifdef OLD_BETA
        dz = 0.5 * (myData[k+1][J] - myData[k-1][J]);
        tmp = (a - 1.0) * log(myData[k][J]) + (b - 1.0) * log(1.0 - myData[k][J]) + factor;
        pdf[k] = exp(tmp) * dz;
#else
        // RV change beta computation
        z1 = 0.5 * (myData[k][J] + myData[k-1][J]);
        z2 = 0.5 * (myData[k+1][J] + myData[k][J]);
        pdf[k] = betai(a,b,z2)-betai(a,b,z1);  
#endif
      }

      // Normalize pdf
      double integral = 0.0;
      for  (int k=0; k<nPoints; k++)
        integral += pdf[k];
//      // For debugging
      //cout << "Integral " << integral << endl;
      if (fabs(1.0-integral)>1.e-2) {
        cerr << "Beta PDF ("<<mean<<" , "<<var<<") not well discretized, Integral: " << integral <<endl;
#ifndef OLD_BETA
        throw(-1);
#endif
      }  
      integral = 1.0 / integral;
      for  (int k=0; k<nPoints; k++)
        pdf[k] *= integral;

    }
    
    return;
  } 

  /** \brief Convolute the flamelet with a Beta-PDF given PDF parameters mean and variance for given variable */
  double FlameletConvoluteWithPDF(double mean, double var, string variable)
  {
    // Compute first Beta PDF given mean and var (this updates array pdf)
    ComputeBetaPDF(mean, var);
    
    int J = GetFlameletVariableIndex(variable);
    int J2 = GetFlameletVariableIndex("Z");
    
    double value = 0.0;
    if ((variable == "density") || (variable == "rho0"))
    {
       for (int k=0; k<nPoints; k++)
          value += pdf[k] / myData[k][J];
       return 1.0 / value;
    }
    else
    {
      for (int k=0; k<nPoints; k++)
         value += pdf[k] * myData[k][J];
      return value;
    }
  }

  /** \brief Convolute the vector x with a Beta-PDF given PDF parameters mean and variance */
  double FlameletConvoluteWithPDF(double mean, double var, double *x)
  {
    // Compute first Beta PDF given mean and var (this updates array pdf)
    ComputeBetaPDF(mean, var);

    double value = 0.0;
    for (int k=0; k<nPoints; k++)
      value += pdf[k] * x[k];
    return value;
  }


// Functions to compute beta PDF
  double betai(double a, double b, double x)
  {
    double bt,betai;
    if ((x<0.0)||(x>1.0)) {cout << "bad argument X in BETAI. X = " << x << endl;
                           throw(-1);}

    if ((x == 0.0)||(x == 1.0)) { bt=0.0; }
    else { bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x)); }
 
    if (x < (a+1.0)/(a+b+2.0)) { betai = bt*betacf(a,b,x)/a; }
    else { betai = 1.0 - bt*betacf(b,a,1.0-x)/b; }

    return betai;
  }

  double betacf(double a, double b, double x)
  {
    double eps=3.0e-7,am,bm,az,qab,qap,qam,bz,aold,em,tem,d,ap,bp,app,bpp;
    int m,itmax=10000;

    am = 1.0;
    bm = 1.0;
    az = 1.0;
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    bz = 1.0 - qab*x/qap;

    m=0;
    aold=0.0;
    while (m <= itmax)
    {
      m = m + 1;
      em = (double) m;
      tem = em + em;
      d = em*(b-(double) m)*x/((qam+tem)*(a+tem));
      ap = az + d*am;
      bp = bz + d*bm;
  
      d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem));
      app = ap + d*az;
      bpp = bp + d*bz;
  
      aold = az;
      am = ap/bpp;
      bm = bp/bpp;
      az = app/bpp;
      bz = 1.0;
  
      if (fabs(az-aold)<eps*fabs(az)) break;
    }

    if (m>itmax) { cout << "A or B too big, or ITMAX too small" << endl;
                   cout << "A =" << a <<endl;
                   cout << "B =" << b <<endl;
                   cout << "Number of iterations: "<< m <<endl;
                   cout << "Error: "<< fabs(az-aold)/fabs(az) << " vs tolerance: "<< eps<<endl;
                   throw(-1);}
    return az;
  }

};
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
