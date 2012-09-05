#include "CreateChemTable.h"
/**********************************************************************************************/

/* Test if flamelet progress variable is larger or smaller by checking temperature in the filename */
/* Used to order the flamelet in increasing value of progress variable */
bool IsProgSmaller(string fl1, string fl2)
{
  size_t pos;
  string Tst;
  double T1, T2;
  
  pos = fl1.find("Tst");
  Tst = fl1.substr(pos+3);
  istringstream buf1(Tst);
  buf1 >> T1;
  
  pos = fl2.find("Tst");
  Tst = fl2.substr(pos+3);
  istringstream buf2(Tst);
  buf2 >> T2;
  
  return (T2 > T1);
}

/* Test if flamelet scalar dissipation rate is larger or smaller by checking the filename */
/* Used to order the flamelet in increasing value of scalar dissipation rate */
bool IsChiLarger(string fl1, string fl2)
{
  size_t pos1, pos2;
  string chist;
  double chi1, chi2;
  
  pos1 = fl1.find("chi") + 3;
  pos2 = fl1.find("tf");
  chist = fl1.substr(pos1, pos2 - pos1);
  istringstream buf1(chist);
  buf1 >> chi1;
  
  pos1 = fl2.find("chi") + 3;
  pos2 = fl2.find("tf");
  chist = fl2.substr(pos1, pos2 - pos1);
  istringstream buf2(chist);
  buf2 >> chi2;
  
  return (chi2 > chi1);
}

void CreateChemTableInitialize(ParamMap *myInputPtr)
/* Read parameters from input file */
{
  Param    *p = NULL;  
  string   FlameletListFilename, file, buffer_str;
  ifstream fin;
  
  version = Ver;
  //////////////////////////////////////////////////////////////////////////////////////////
  // Read input parameters
  //////////////////////////////////////////////////////////////////////////////////////////
  
  ChemTableFilename = myInputPtr->getStringParam("CHEMTABLE_FILENAME", "table.out");
  ChemTableType     = myInputPtr->getStringParam("CHEMTABLE_TYPE", "SPECIES");
  CombustionModel   = myInputPtr->getStringParam("COMBUSTION_MODEL", "FPVA");
  TecplotOutput     = myInputPtr->getStringParam("TECPLOT_OUTPUT", "YES");
  ScaleThirdDim     = myInputPtr->getStringParam("SCALE_THIRDDIM", "ADAPTIVE");
  NormalizedCoordinates = myInputPtr->getStringParam("CHEMTABLE_NORMALIZATION","NO");
  DefinitionProg        = myInputPtr->getStringParam("PROG_DEFINITION","USUAL");
  WriteSinglePrecision  = myInputPtr->getStringParam("WRITE_SINGLE_PRECISION","NO");
  
  if (ChemTableType == "COEFF")
  {
    VarToOutput.push_back("ROM");
    VarToOutput.push_back("T0");
    VarToOutput.push_back("E0");
    VarToOutput.push_back("GAMMA0");
    VarToOutput.push_back("AGAMMA");
    VarToOutput.push_back("MU0");
    VarToOutput.push_back("AMU");
    VarToOutput.push_back("LOC0");
    VarToOutput.push_back("ALOC");
    spvar = 9;
  }
  
  nZm = myInputPtr->getIntParam("N_ZMEAN", "100");
  nZv = myInputPtr->getIntParam("N_ZVAR", "50");
  nCC = myInputPtr->getIntParam("N_THIRDDIM", "100");
  
  Zst = myInputPtr->getDoubleParam("Z_ST", "0.03");

  // Read species and weights in progress variable
  if (CombustionModel == "FPVA")
  {
    CombustionRegime = 2;
    p = myInputPtr->getParam("VAR_PROG");
    if (p == NULL)
    {
      cerr << "### Missing input VAR_PROG in input file! ###" << endl;
      throw(-1);
    }
    int nVarProg = p->getSize() - 1;
    if (nVarProg <= 0)
    {
      cerr << "### Missing variables for progress variable in input file! ###" << endl;
      throw(-1);
    }
    for (int i = 0; i < nVarProg; i++)
      VarProg.push_back(p->getString(i + 1));
    
    p = myInputPtr->getParam("WEIGHT_PROG");
    if (p == NULL)
    {
      cout << "No weights found in input file for progress variable, set to default value 1.0" << endl;
      for (int i = 0; i < nVarProg; i++)
        WeightProg.push_back(1.0);
    }
    else
    {
      if (nVarProg != p->getSize()-1)
      {
        cerr << "### Number of variables for progress variable does not match the number of weights in input file! ###" << endl;
        throw(-1);
      }
      for (int i = 0; i < nVarProg; i++)
        WeightProg.push_back(p->getDouble(i + 1));
    }
  }
  else if (CombustionModel == "STEADY")
  {
    CombustionRegime = 1;
    VarProg.push_back("chi");
    WeightProg.push_back(1.0);
  }
  else if (CombustionModel == "ANALYTIC")
  {
    CombustionRegime = -1;
    VarToOutput.push_back("ZM");
    VarToOutput.push_back("ZV");
    VarToOutput.push_back("CM");
    VarToOutput.push_back("Linear_Fcn");
    VarToOutput.push_back("Parabolic_Fcn");
    VarToOutput.push_back("Cubic_Fcn");
    VarToOutput.push_back("General_Fcn");
  }
  else
  {
    cerr << "### Combustion model not recognized! ###" << endl;
    throw(-1);
  }
  
  // Add source term for progress variable and progress variable
  if (CombustionModel == "FPVA")
  {
    VarToOutput.push_back("SRC_PROG");
    VarToOutput.push_back("PROG");
    VarToOutput.push_back("HeatRelease");
    VarToOutput.push_back("rho0");
  }
  else if (CombustionModel == "STEADY")
    VarToOutput.push_back("PROG");
  
    
  // Read and sort list of flamelet files
  if (CombustionModel != "ANALYTIC")
  {    
    // Read species to store in table
    p = myInputPtr->getParam("SPECIES_TABLE");
    if (p != NULL)
    {
      int nVarTable = p->getSize() - 1;
      for (int i = 0; i < nVarTable; i++)
        VarToOutput.push_back(p->getString(i + 1));
    }

    FlameletListFilename = myInputPtr->getStringParam("FLAMELET_LIST_FILENAME", "FlameletList.txt");
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
    if (CombustionModel == "FPVA")        sort(FlameletList.begin(), FlameletList.end(), IsProgSmaller);
    else if (CombustionModel == "STEADY") sort(FlameletList.begin(), FlameletList.end(), IsChiLarger);
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Summarize inputs to screen
  //////////////////////////////////////////////////////////////////////////////////////////
  
  cout << endl;
  cout << "****************************************************************" << endl;
  cout << "**************** CreateChemTable Initialization ****************" << endl;
  cout << "****************************************************************" << endl;
  cout << endl;
  
  cout << "CHEMTABLE_FILENAME:      " << ChemTableFilename << endl;
  cout << "CHEMTABLE_TYPE:          " << ChemTableType << endl;
  cout << "CHEMTABLE_NORMALIZATION: " << NormalizedCoordinates << endl; 
  cout << "DEFINITION_PROG:         " << DefinitionProg <<endl;
  cout << "WRITE_SINGLE_PRECISION   " << WriteSinglePrecision << endl;
  cout << "COMBUSTION_MODEL:        " << CombustionModel << endl;
  cout << "TECPLOT_OUTPUT:          " << TecplotOutput << endl;
  cout << "SCALE_THIRDDIM:          " << ScaleThirdDim << endl;
  cout << "N_ZMEAN:                 " << nZm << endl;
  cout << "N_ZVAR:                  " << nZv << endl;
  cout << "N_THIRDDIM:              " << nCC << endl;
  cout << "Z_ST:                    " << Zst << endl;
  
  if (CombustionModel != "ANALYTIC")
  {
    cout << "PROGRESS VARIABLE:       " << WeightProg[0] << "    " << VarProg[0] << endl;
    for (int i=1; i<VarProg.size(); i++)
      cout << "                         " << WeightProg[i] << "    " << VarProg[i] << endl;
    cout << "STORED VARIABLES:        " << VarToOutput[0] << endl;
    for (int i=1; i<VarToOutput.size(); i++)
      cout << "                         " << VarToOutput[i] << endl;
    cout << "FLAMELET FILES:          " << FlameletList[0] << endl;
    for (int i=1; i<FlameletList.size(); i++)
      cout << "                         " << FlameletList[i] << endl;
  }
  
  cout << endl;
  cout << "****************************************************************" << endl << endl;
  
  return;
}

/**********************************************************************************************/
void PrepareTable()
/* Allocate arrays and create axis of table */
{
  //////////////////////////////////////////////////////////////////////////////////////////
  // ZMean
  //////////////////////////////////////////////////////////////////////////////////////////
  Zm = new double[nZm];
  
  // One third of the points between 0 and Zst : uniform mesh
  // Two third between Zst and 1 : linear growth
  // Create such mesh only if Zst sufficiently small
  if (Zst > 0.3)
  {
     double dz = 1.0 / (double) (nZm-1);
     for (int i=0; i<nZm; i++)
        Zm[i] = (double) (i) * dz;
  }
  else
  {
    // Linear mesh for [0,Zst]
    int nZcut = nZm / 3;
    double dz = Zst / (double) (nZcut);
    for (int i=0; i<=nZcut; i++)
      Zm[i] = (double) (i) * dz;
  
    // Mesh with linear growth to reach Z=1
    // Zm(i) = (i-nZcut)*dz + 0.5*(i-nZcut)(i-nZcut-1)*alpha + Zst
    // alpha from: Zm(nZm-1) = 1.0
    double m = (double) (nZm - 1 - nZcut);
    double r = 1.0 - Zst;
    double alpha = 2.0 * (r - m * dz) / (m * (m - 1.0));
    for (int i=nZcut+1; i< nZm; i++)
    {
      double a = (double) (i - nZcut);
      Zm[i] = a * dz + 0.5 * a * (a - 1.0) * alpha + Zst;
    }
  }
//  // Debugging check
//    cout << endl;
//    for (int i=0; i<nZm; i++)
//      cout << i << " " << Zm[i] << endl;
  

  //////////////////////////////////////////////////////////////////////////////////////////
  // ZVar
  //////////////////////////////////////////////////////////////////////////////////////////
  Zv = new double[nZv];
  if (NormalizedCoordinates != "YES")
  { 
    // Quadratic mesh between 0 and 0.25 (variance cannot be larger than Z*(1-Z))
    for (int j=0; j<nZv; j++)
      Zv[j] = 0.25 * pow((double) (j) / (double) (nZv-1), 2.0);
      // zV[j] = 0.25 * ((double) (j) / (double) (nZv-1));
  } else {
    // Zv is Zv/(z*(1-z)) in fact
    for (int j=0; j<nZv; j++)
      Zv[j] = pow((double) (j) / (double) (nZv-1),2.7);
  }
  
// //Debugging check
//  cout << endl;
//  for (int i=0; i<nZv; i++)
//    cout << i << " " << Zv[i] << endl;

  //////////////////////////////////////////////////////////////////////////////////////////
  // Number of variables
  //////////////////////////////////////////////////////////////////////////////////////////
  nVar = VarToOutput.size();
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Third dimension
  //////////////////////////////////////////////////////////////////////////////////////////  
  if (CombustionModel == "ANALYTIC")
    nCA = nCC;
  else
    nCA = FlameletList.size();
  getMem3D(&CA, 0, nZm-1, 0, nZv-1, 0, nCA-1, "CreateChemTableInitialize: CA", true);
  nC = nCA;
  if (ScaleThirdDim != "ADAPTIVE")
  {
    CC = new double[nCC];
    nC = nCC;
  }
 
  //////////////////////////////////////////////////////////////////////////////////////////
  // Tables
  //////////////////////////////////////////////////////////////////////////////////////////    
  getMem4D(&TabCA, 0, nZm-1, 0, nZv-1, 0, nCA-1, 0, nVar-1, "CreateChemTableInitialize: adaptive table", true);
  if (ScaleThirdDim != "ADAPTIVE") 
    getMem4D(&TabCC, 0, nZm-1, 0, nZv-1, 0, nCC-1, 0, nVar-1, "CreateChemTableInitialize: Cartesian table", true);
    
  return;
}

/**********************************************************************************************/

void CreateCartesianTable()
/* Use adaptive table TabCA to create Cartesian table TabCC */
{
  double min3 = +1.0e10;
  double max3 = -1.0e10;
  if (NormalizedCoordinates != "YES")
  {
    // Determine minimum and maximum value for third dimension
    for (int i=0; i<nZm; i++)
      for (int j=0; j<nZv; j++)
        for (int k=0; k<nCA; k++)
        {
          if (CA[i][j][k] > max3)
          max3 = CA[i][j][k];
          if (CA[i][j][k] < min3)
            min3 = CA[i][j][k];
        }
 } else {
   max3 = 1.0;
   min3 = 0.0;
 }

  // Construct third dimension vector CC
  if (ScaleThirdDim == "LIN")
    for (int k=0; k<nCC; k++)
      CC[k] = min3 + (max3 - min3) * ((double) k  / (double) (nCC-1));
  else if (ScaleThirdDim == "QUAD")
    for (int k=0; k<nCC; k++)
      CC[k] = min3 + (max3 - min3) * pow((double) k  / (double) (nCC-1), 2.0);
  else if (ScaleThirdDim == "LOG")
    for (int k=0; k<nCC; k++)
      CC[k] = min3 * pow((max3 / min3), ((double) k  / (double) (nCC-1)));    
  else
  {
    cerr << "### " << ScaleThirdDim << " scale for Cartesian table is unknown! ###" << endl;
    throw(-1);
  }


  // Loop over the three mapping directions in Cartesian table
  for (int i=0; i<nZm; i++)
    for (int j=0; j<nZv; j++)
      for (int k=0; k<nCC; k++)
      {
        double alpha_up, alpha_down;
        
        // Find the two flamelets right above and right below
        double err_up   = +1.0e15;
        double err_down = -1.0e15;

        int kk_up;
        int kk_down;

#ifndef INTERP_OLD
        BinarySearch(kk_down, alpha_down, CA[i][j], 0, nCA-1, CC[k]);
        kk_up = kk_down + 1;
        alpha_up = 1.0 - alpha_down;
        int kk_aux = kk_down;
        double aux =  alpha_down; 
#else
        kk_up   = -1;
        kk_down = -1;
        // Loop over the third dimension of the adaptive table
        for (int kk=0; kk<nCA; kk++)
        {
          
          double err = CA[i][j][kk] - CC[k];

          if ((err >= 0.0) && (err <= err_up))
          {
            kk_up  = kk;
            err_up = err;
          }
          if ((err <= 0.0) && (err >= err_down))
          {
            kk_down  = kk;
            err_down = err;
          }
        }
        // Interpolate
        if ((kk_up == -1) || (kk_down == -1))
        {
          if (kk_up == -1)
          {
            alpha_up   = 0.0;
            alpha_down = 1.0;
            kk_up = 0;
          }
          else
          {
            alpha_up   = 1.0;
            alpha_down = 0.0;
            kk_down = 0;
          }
        }
        else
        {
          if (kk_up == kk_down)
          {
            alpha_up   = 1.0;
            alpha_down = 0.0;
          }
          else
          {
            alpha_up   = (CC[k] - CA[i][j][kk_down]) / (CA[i][j][kk_up] - CA[i][j][kk_down]);
            alpha_down = (CA[i][j][kk_up] - CC[k])   / (CA[i][j][kk_up] - CA[i][j][kk_down]);
          }
        }
#endif
        //if (kk_aux != kk_down) {cout<<"Not same kk_down "<<kk_aux<<" "<<kk_down<<endl;throw(-1);}
        //if (aux != alpha_down) {cout<<"Not same alpha_down "<<aux<<" "<<alpha_down<<endl;throw(-1);}
        
        for (int n=0; n<nVar; n++)
          TabCC[i][j][k][n] = alpha_up   * TabCA[i][j][kk_up][n]
                            + alpha_down * TabCA[i][j][kk_down][n];

      }

  return;
}

/**********************************************************************************************/

void ComputeEnergyPerturbation(double &T_p_hat, double &mu_p_tild, double &lamOcp_p_tild, 
                               double Zm, double Zv, Flamelet *myFlamelet, 
                               double R_0_tild, double rho_0_bar, double Delta_e_tild)
{
  myFlamelet->ComputeEnergyPerturbation(rho_0_bar, Delta_e_tild);
  
  T_p_hat       = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "RT1") / R_0_tild;
  mu_p_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "mu1");
  lamOcp_p_tild = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "lamOcp1");  
}

/**********************************************************************************************/

void CalculateCoefficients(double &R_0_tild, double &T_0_hat, double &e_0_tild, double &gamma_0_tild, double &a_gamma, 
                           double &mu_0_tild, double &a_mu, double &lamOcp_0_tild, double &a_lamOcp, double &rho_0_bar,
                           double Zm, double Zv, Flamelet *myFlamelet)
/* Compute coefficients for energy, viscosity and thermal diffusion */
{
  double RT_0_tild;                                                                               // Convoluted variables at flamelet solution  
  double Delta_e_tild, T_p_hat, T_m_hat, mu_p_tild, mu_m_tild, lamOcp_p_tild, lamOcp_m_tild;      // Perturbed variables (for Delta_e_tild perturbation: positive -> ...p..., negative -> ...m...)
  
  // Compute convoluted values at flamelet solution
  R_0_tild      = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"R0");
  rho_0_bar     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"rho0");
  RT_0_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"RT0");
  e_0_tild      = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"E0");
  mu_0_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"mu0");
  lamOcp_0_tild = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"lamOcp0");
  T_0_hat       = RT_0_tild / R_0_tild;

  Delta_e_tild = 5000.0;
  ComputeEnergyPerturbation(T_p_hat, mu_p_tild, lamOcp_p_tild, Zm, Zv, myFlamelet, R_0_tild, rho_0_bar,  Delta_e_tild);
  ComputeEnergyPerturbation(T_m_hat, mu_m_tild, lamOcp_m_tild, Zm, Zv, myFlamelet, R_0_tild, rho_0_bar, -Delta_e_tild);
  
  // Derivatives: de/dT and d^2e/dT^2 at T=T0
  double dTm = T_0_hat - T_m_hat;
  double dTp = T_p_hat - T_0_hat;
  double dT  = T_p_hat - T_m_hat;
  double dedT = ( dTm * dTm * (e_0_tild + Delta_e_tild) + (dTp - dTm) * dT * e_0_tild - dTp * dTp * (e_0_tild - Delta_e_tild) ) / (dTp * dTm * dT);
  double d2edT2 = 2.0 * ( dTm * (e_0_tild + Delta_e_tild) - dT * e_0_tild + dTp * (e_0_tild - Delta_e_tild) ) / (dTp * dTm * dT);
  
  // Compute coefficients
  gamma_0_tild = R_0_tild / dedT + 1.0;
  a_gamma      = - d2edT2 * (gamma_0_tild - 1.0) * (gamma_0_tild - 1.0) / R_0_tild;
  a_mu         = ( dTm * dTm * mu_p_tild + (dTp - dTm) * dT * mu_0_tild - dTp * dTp * mu_m_tild ) / (dTp * dTm * dT);
  a_lamOcp     = ( dTm * dTm * lamOcp_p_tild + (dTp - dTm) * dT * lamOcp_0_tild - dTp * dTp * lamOcp_m_tild ) / (dTp * dTm * dT);
}

/**********************************************************************************************/

void CalculateCoefficients(double &R_0_tild, double &T_0_hat, double &e_0_tild, double &gamma_0_tild, double &a_gamma, 
                           double &mu_0_tild, double &a_mu, double &lamOcp_0_tild, double &a_lamOcp, double &rho_0_bar,
                           double rho0barDeltaEtild, double Zm, double Zv, Flamelet *myFlamelet)
/* Compute coefficients for energy, viscosity and thermal diffusion */
{
  double RT_0_tild;                                                                               // Convoluted variables at flamelet solution  
  double Delta_e_tild, T_p_hat, T_m_hat, mu_p_tild, mu_m_tild, lamOcp_p_tild, lamOcp_m_tild;      // Perturbed variables (for Delta_e_tild perturbation: positive -> ...p..., negative -> ...m...)
  
  // Compute convoluted values at flamelet solution
  R_0_tild      = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"R0");
  rho_0_bar     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"rho0");
  RT_0_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"RT0");
  e_0_tild      = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"E0");
  mu_0_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"mu0");
  lamOcp_0_tild = myFlamelet->FlameletConvoluteWithPDF(Zm,Zv,"lamOcp0");
  T_0_hat       = RT_0_tild / R_0_tild;

  // Compute quantities if energy perturbed by +/- rho_0_bar * Delta_E_tild
  T_p_hat       = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "RTpl") / R_0_tild;
  mu_p_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "mupl");
  lamOcp_p_tild = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "lamOcppl");  
  T_m_hat       = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "RTmi") / R_0_tild;
  mu_m_tild     = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "mumi");
  lamOcp_m_tild = myFlamelet->FlameletConvoluteWithPDF(Zm, Zv, "lamOcpmi");  

  // Compute quantities if energy perturbed by +/- Delta_e
  Delta_e_tild = rho0barDeltaEtild / rho_0_bar;
  
  // Derivatives: de/dT and d^2e/dT^2 at T=T0
  double dTm = T_0_hat - T_m_hat;
  double dTp = T_p_hat - T_0_hat;
  double dT  = T_p_hat - T_m_hat;
  double dedT = ( dTm * dTm * (e_0_tild + Delta_e_tild) + (dTp - dTm) * dT * e_0_tild - dTp * dTp * (e_0_tild - Delta_e_tild) ) / (dTp * dTm * dT);
  double d2edT2 = 2.0 * ( dTm * (e_0_tild + Delta_e_tild) - dT * e_0_tild + dTp * (e_0_tild - Delta_e_tild) ) / (dTp * dTm * dT);
  
  // Compute coefficients
  gamma_0_tild = R_0_tild / dedT + 1.0;
  a_gamma      = - d2edT2 * (gamma_0_tild - 1.0) * (gamma_0_tild - 1.0) / R_0_tild;
  a_mu         = ( dTm * dTm * mu_p_tild + (dTp - dTm) * dT * mu_0_tild - dTp * dTp * mu_m_tild ) / (dTp * dTm * dT);
  a_lamOcp     = ( dTm * dTm * lamOcp_p_tild + (dTp - dTm) * dT * lamOcp_0_tild - dTp * dTp * lamOcp_m_tild ) / (dTp * dTm * dT);
}

/**********************************************************************************************/

void TestCoeffApprox(double nZ, double rho0barDeltaEtild, double Zm, double Zv, Flamelet *myFlamelet)
/* Test coefficients by comparing with exact solution */
{
  double R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar;
  
  for (int k=0; k<nZ; k+=50)
  {
    cout << endl;
    cout << "Z           : " << myFlamelet->GetFlameletVariable(k, "Z")            << endl;
    cout << "Density     : " << myFlamelet->GetFlameletVariable(k, "density")      << "   " << myFlamelet->GetFlameletVariable(k, "rho0")    << "   " << myFlamelet->GetFlameletVariable(k, "density")      - myFlamelet->GetFlameletVariable(k, "rho0")    << endl;
    cout << "Temperature : " << myFlamelet->GetFlameletVariable(k, "temperature")  << "   " << myFlamelet->GetFlameletVariable(k, "T0")      << "   " << myFlamelet->GetFlameletVariable(k, "temperature")  - myFlamelet->GetFlameletVariable(k, "T0")      << endl;
    cout << "Viscosity   : " << myFlamelet->GetFlameletVariable(k, "mu")           << "   " << myFlamelet->GetFlameletVariable(k, "mu0")     << "   " << myFlamelet->GetFlameletVariable(k, "mu")           - myFlamelet->GetFlameletVariable(k, "mu0")     << endl;
    cout << "LambdaOcp   : " << myFlamelet->GetFlameletVariable(k, "lambdaOverCp") << "   " << myFlamelet->GetFlameletVariable(k, "lamOcp0") << "   " << myFlamelet->GetFlameletVariable(k, "lambdaOverCp") - myFlamelet->GetFlameletVariable(k, "lamOcp0") << endl;
    cout << endl;
  }
 
  // Open output file
  ofstream fout;

  fout.open("TestCoeff.dat");
  fout << "TITLE=\"Check coefficients for chemistry table\"" << endl;
  fout << "VARIABLES=\"de\"  ";
  fout << "\"Tex\" \"Tap\" \"Terror\" ";
  fout << "\"MUex\" \"MUap\" \"MUerror\" ";
  fout << "\"LAMOCPex\" \"LAMOCPap\" \"LAMOCPerror\" ";
  fout << endl;

  // Compute coefficients
  CalculateCoefficients(R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar, rho0barDeltaEtild, Zm, Zv, myFlamelet);
  
  cout << "Coefficients calculated" << endl;
  
  for (int i=-20; i<21; i++)
  {
    double Delta_e_tild = double(i) * 5000.0;
    
    // "Exact" solution
    double T_ex, mu_ex, lamOcp_ex;
    ComputeEnergyPerturbation(T_ex, mu_ex, lamOcp_ex, Zm, Zv, myFlamelet, R_0_tild, rho_0_bar, Delta_e_tild);
    
    // Approximation
    double T_ap, mu_ap, lamOcp_ap;
    if (a_gamma == 0.0)
      T_ap = T_0_hat + (gamma_0_tild - 1.0) / R_0_tild * Delta_e_tild;
    else
      T_ap = T_0_hat + (gamma_0_tild - 1.0) / a_gamma * (exp(a_gamma * Delta_e_tild / R_0_tild) - 1.0);
    mu_ap = mu_0_tild + a_mu * (T_ap - T_0_hat);
    lamOcp_ap = lamOcp_0_tild + a_lamOcp * (T_ap - T_0_hat);
    
    // Output to file
    fout << Delta_e_tild << " ";
    fout << T_ex      << " " << T_ap      << " " << T_ap-T_ex           << " ";
    fout << mu_ex     << " " << mu_ap     << " " << mu_ap-mu_ex         << " ";
    fout << lamOcp_ex << " " << lamOcp_ap << " " << lamOcp_ap-lamOcp_ex << " ";
    fout << endl;
  }

  fout.close();
  
}

/**********************************************************************************************/

void CreateAnalyticalTable()
/* Create table constructed analytically to test interpolation */
{
  double Zm_min = 0.0, Zm_max = 1.0;
  double ZV_min = 0.0, Zv_max = 0.25;
  double Cm_min = 0.0, Cm_max = 0.5;
  
  
  for (int i=0; i<nZm; i++)
    for (int j=0; j<nZv; j++)
      for (int k=0; k<nCA; k++)
      {
        CA[i][j][k] = Cm_min + (Cm_max - Cm_min) / double(nCA - 1) * double(k);
        TabCA[i][j][k][0] = Zm[i];
        TabCA[i][j][k][1] = Zv[j];
        TabCA[i][j][k][2] = CA[i][j][k];
        TabCA[i][j][k][3] = linearFunction(Zm[i], Zv[j], CA[i][j][k]);
        TabCA[i][j][k][4] = parabolicFunction(Zm[i], Zv[j], CA[i][j][k]);
        TabCA[i][j][k][5] = cubicFunction(Zm[i], Zv[j], CA[i][j][k]);
        TabCA[i][j][k][6] = otherAnalyticalFunction(Zm[i], Zv[j], CA[i][j][k]);
      }
  
}

/**********************************************************************************************/

void OutputTableStatistics()
/* Write statistics to standard output */
{
  cout << endl << endl << "***** Statistics *****" << endl << endl;
  
  cout << "Version of CreateChemTable:  " << version << endl;
  cout << "Table type:                  " << ChemTableType << " " << ScaleThirdDim << endl;
  cout << "Combustion model:            " << CombustionModel << endl;
  cout << "Reference pressure:          " << Preference << endl;
  cout << "Oxidizer temperature:        " << Toxi << endl;
  cout << "Fuel temperature:            " << Tfuel << endl;
  cout << "Progress variable:         ";
  for (int l=0; l<VarProg.size(); l++)
    cout << "+(" << WeightProg[l] << ")*" << VarProg[l];
  cout << endl;
  cout << "Dimensions - ZMean:          " << nZm << endl;
  cout << "             ZVar:           " << nZv << endl;
  cout << "             Prog:           " << nC << endl;
  
  cout << endl << endl;
  
  // Get minimum and maximum value of progress variable
  double ProgMin = +1.0e10;
  double ProgMax = -1.0e10;
  if (ScaleThirdDim == "ADAPTIVE")   // 3D vector
  {
    for (int i=0; i<nZm; i++)
      for (int j=0; j<nZv; j++)
        for (int k=0; k<nCA; k++)
        {
          if (CA[i][j][k] > ProgMax)
            ProgMax = CA[i][j][k];
          if (CA[i][j][k] < ProgMin)
            ProgMin = CA[i][j][k];
        }
  }
  else   // only 1D vector
  {
    ProgMax = CC[nC-1];
    ProgMin = CC[0];
  }    
  
  cout << "ZMean:   ";
  cout.width(8);
  cout << left << Zm[0];
  cout << " -> ";
  cout.width(8);
  cout << Zm[nZm-1] << endl;
  cout << "ZVar:    ";
  cout.width(8);
  cout << Zv[0];
  cout << " -> ";
  cout.width(8);
  cout << Zv[nZv-1] << endl;
  cout << "Prog:    ";
  cout.width(8);
  cout << ProgMin;
  cout << " -> ";
  cout.width(8); 
  cout << ProgMax << endl;
  
  cout << endl << endl;
  
  for (int n=0; n<nVar; n++)
  {
    double VarMin = +1.0e10;
    double VarMax = -1.0e10;
    for (int i=0; i<nZm; i++)
      for (int j=0; j<nZv; j++)
        for (int k=0; k<nC; k++)
        {
          if (Table[i][j][k][n] > VarMax)
            VarMax = Table[i][j][k][n];
          if (Table[i][j][k][n] < VarMin)
            VarMin = Table[i][j][k][n];
        } 
    cout.width(15);
    cout << VarToOutput[n];
    cout << ":       ";
    cout.width(8);
    cout << left << VarMin;
    cout << " -> ";
    cout.width(8);
    cout << VarMax << endl;
  }
  
  
  return;
}

/**********************************************************************************************/

void WriteTable()
/* Write table to file */
{
  FILE   *inFp = 0;
  char   filename[kStrLong], buffer_m[kStrMed];
  char   *p = 0, dummy;
  size_t dum;
  int    mask = 0;

  cout << endl << "***** Write table to file *****" << endl << endl;

  // ***************
  // Open table file
  // ***************
  dum = ChemTableFilename.copy(&filename[0], kStrLong);
  dum = ChemTableFilename.size();
  if (dum < kStrLong)
    strcpy(&filename[dum], "\0");
  else
  {
    cerr << "### Chemistry table file name is too long! ###" << endl;
    throw(-1);
  }
  if (!(inFp = fopen(filename, "wb")))
  {
    cerr << "### Could not open chemistry table file " << filename << " ###" << endl;
    throw(-1);
  }
  
  // ************************
  // Write table general info
  // ************************
  dum = fwrite(&version, sizeof(double), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: version ###" << endl;
    throw(-1);
  }
  
  dum = CombustionModel.copy(&buffer_m[0], kStrMed);
  dum = CombustionModel.size();
  if (dum < kStrMed)
    strcpy(&buffer_m[dum], "\0");
  else
  {
    cerr << "### Combustion model name is too long! ###" << endl;
    throw(-1);
  }
  dum = fwrite(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed)
  {
    cerr << "### Error writing file: Combustion model ###" << endl;
    throw(-1);
  }
  
  dum = ChemTableType.copy(&buffer_m[0], kStrMed);
  dum = ChemTableType.size();
  if (dum < kStrMed)
    strcpy(&buffer_m[dum], "\0");
  else
  {
    cerr << "### Chem Table type is too long! ###" << endl;
    throw(-1);
  }
  dum = fwrite(&buffer_m[0], sizeof(char), kStrMed, inFp);
  if (dum != kStrMed)
  {
    cerr << "### Error writing file: Chem Table type ###" << endl;
    throw(-1);
  }
  
  dum = fwrite(&Preference, sizeof(double), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: Preference ###" << endl;
    throw(-1);
  }
  dum = fwrite(&Toxi, sizeof(double), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: Toxi ###" << endl;
    throw(-1);
  }
  dum = fwrite(&Tfuel, sizeof(double), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: Tfuel ###" << endl;
    throw(-1);
  }
  
  // **********************
  // Write table dimensions
  // **********************
  dum = fwrite(&nZm, sizeof(int), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: nZm ###" << endl;
    throw(-1);
  }
  dum = fwrite(&nZv, sizeof(int), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: nZv ###" << endl;
    throw(-1);
  }
  
  dum = fwrite(&nC, sizeof(int), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: nC ###" << endl;
    throw(-1);
  }

  dum = fwrite(&nVar, sizeof(int), 1, inFp);
  if (dum != 1)
  {
    cerr << "### Error writing file: nVar ###" << endl;
    throw(-1);
  }

  // ***********************
  // Write table coordinates 
  // ***********************
  dum = fwrite(Zm, sizeof(double), nZm, inFp);
  if (dum != nZm)
  {
    cerr << "### Error writing file: Zm ###" << endl;
    throw(-1);
  }
  dum = fwrite(Zv, sizeof(double), nZv, inFp);
  if (dum != nZv)
  {
    cerr << "### Error writing file: Zv ###" << endl;
    throw(-1);
  }
    
  if (ScaleThirdDim == "ADAPTIVE")   // 3D vector
  {
    for (int i=0; i<nZm; i++)
      for (int j=0; j<nZv; j++)
        for (int k=0; k<nC; k++)
        {
          dum = fwrite(&CA[i][j][k], sizeof(double), 1, inFp);
          if (dum != 1)
          {
            cerr << "### Error writing file: CA -> " << i << " " << j << " " << k << " ###" << endl;
            throw(-1);
          }
        }
  }
  else   // only 1D vector
  {
    dum = fwrite(CC, sizeof(double), nC, inFp);
    if (dum != nC)
    {
      cerr << "### Error writing file: CC ###" << endl;
      throw(-1);
    }
  }    

  // ***********************
  // Write name of variables
  // ***********************
  for (int l = 0; l < nVar; l++)
  {
    dum = VarToOutput[l].copy(&buffer_m[0], kStrMed);
    dum = VarToOutput[l].size();
    if (dum < kStrMed)
      strcpy(&buffer_m[dum], "\0");
    else
    {
      cerr << "### Variable name " << VarToOutput[l] << " is too long! ###" << endl;
      throw(-1);
    }
    dum = fwrite(&buffer_m, sizeof(char), kStrMed, inFp);
    if (dum != kStrMed)
    {
      cerr << "### Error writing file: Variable name " << VarToOutput[l] << " ###" << endl;
      throw(-1);
    }
  }

  // ****************
  // Write table data 
  // ****************
  if (WriteSinglePrecision != "YES")
  {

  for (int i = 0; i < nZm; i++)
    for (int j = 0; j < nZv; j++)
      for (int k = 0; k < nC; k++)
        for (int n = 0; n < nVar; n++)
        {
          dum = fwrite(&Table[i][j][k][n], sizeof(double), 1, inFp);
          if (dum != 1)
          {
            cerr << "### Error reading file: Table -> " << i << " " << j << " " << k << " " << n << " ### " << endl;
            throw(-1);
          }
        }

  } else {

  // Allocate
  float  ****Table_single;
  getMem4D(&Table_single, 0, nZm-1, 0, nZv-1, 0, nCC-1, 0, nVar-1, "CreateChemTableInitialize: Single precision Cartesian table", true);
  // Copy
   for (int i = 0; i < nZm; i++)
    for (int j = 0; j < nZv; j++)
      for (int k = 0; k < nC; k++)
        for (int n = 0; n < nVar; n++)
        {Table_single[i][j][k][n] = (float) Table[i][j][k][n];}
   // Write
   for (int i = 0; i < nZm; i++)
     for (int j = 0; j < nZv; j++)
       for (int k = 0; k < nC; k++)
         for (int n = 0; n < nVar; n++)
         {
           dum = fwrite(&Table_single[i][j][k][n], sizeof(float), 1, inFp);
           if (dum != 1)
           {
             cerr << "### Error reading file: Table_single -> " << i << " " << j << " " << k << " " << n << " ### " << endl;
             throw(-1);
           }
         }

   // Destroy
   if (Table_single != NULL) freeMem4D(Table_single, 0, nZm-1, 0, nZv-1, 0, nCC-1, 0, nVar-1); Table_single = NULL;
  }
  // ****************
  // Close table file
  // ****************
  fclose(inFp);
}

/**********************************************************************************************/

void WriteTableTecplot()
/* Write table to file in Tecplot ASCII format for visualization */
{
  ofstream fout;
  char     filename[kStrLong];
  size_t   dum;

  cout << endl << "***** Write table to tecplot format file *****" << endl << endl;

  // *********
  // Open file
  // *********
  ChemTableFilename.append(".dat");
  dum = ChemTableFilename.copy(&filename[0], kStrLong);
  dum = ChemTableFilename.size();
  if (dum < kStrLong)
    strcpy(&filename[dum], "\0");
  else
  {
    cerr << "### Tecplot file name is too long! ###" << endl;
    throw(-1);
  }
  fout.open(filename);
  if (fout.fail())
  {
    cerr << "### Cannot open output tecplot file ###" << endl;
    throw(-1);
  }

  // ********************
  // Write tecplot header
  // ********************
  fout << "TITLE=\"Chemistry table\"" << endl;
  fout << "VARIABLES=\"I1\" \"I2\" \"I3\" \"ZMean\" \"ZVar\" \"Progress\" ";
  for (int j=0; j<nVar; j++)
    fout << "\"" << VarToOutput[j] << "\" ";
  fout << endl;

  fout << "ZONE I=" << nZm << ", J=" << nZv << ", K=" << nC << ", DATAPACKING=POINT" << endl;
  
  // *******************************
  // Write data for each table point
  // *******************************
  for (int k=0; k<nC; k++)
    for (int j=0; j<nZv; j++)
      for (int i=0; i<nZm; i++)
      {
        double prog;
        if (ScaleThirdDim == "ADAPTIVE")   // 3D vector
          prog = CA[i][j][k];
        else
          prog = CC[k];

        fout << i << "  " << j << "  " << k << "  ";
        fout << Zm[i] << "   " << Zv[j] << "   " << prog;       
        for (int n=0; n<nVar; n++)
          fout << "   " << Table[i][j][k][n];
        fout << endl;
      }
  
  fout.close();
  
  cout << endl << "***** Table written to tecplot format file *****" << endl << endl;
  
  return;
}

/**********************************************************************************************/

void CreateChemTableFinalize()
/* Clean memory */
{
  if (Zm != NULL)    delete [] Zm; Zm = NULL;
  if (Zv != NULL)    delete [] Zv; Zv = NULL;
  if (CC != NULL)    delete [] CC; CC = NULL;
  if (CA != NULL)    freeMem3D(CA, 0, nZm-1, 0, nZv-1, 0, nCA-1); CA = NULL;
  Table = NULL;
  if (TabCA != NULL) freeMem4D(TabCA, 0, nZm-1, 0, nZv-1, 0, nCA-1, 0, nVar-1); TabCA = NULL;
  if (TabCC != NULL) freeMem4D(TabCC, 0, nZm-1, 0, nZv-1, 0, nCC-1, 0, nVar-1); TabCC = NULL;

  cout << endl << "***** CreateChemTable finalized *****" << endl << endl;
  
  return;
}

/**********************************************************************************************/

int main(int argc, char *argv[])
{
  char InputFileName[kStrLong];
  sprintf(InputFileName, "CreateChemTable.in");
  for (int i=1; i<argc; i++)
    strcpy(InputFileName, argv[i]);
  
  
  ParamMap myInput(InputFileName);
  ParamMap *myInputPtr = &myInput.getThisParamMap();

  double prog_max = 0.0;                 // Maximum value of progress variable to check in really sorted
  int    nZPoints;                       // Number of points in flamelet
  string previousFl = "NONE";
  double rho0barDeltaEtild = 5000.0;     // Overall Delta E * rho0_bar for all flamelets
  
  CreateChemTableInitialize(myInputPtr);
  
#ifdef TESTING
  // Test computation of coefficients for 1 specific flamelet (first in list) and for specific ZMean and ZVar (define at top)
  myFlamelet = new Flamelet;
  myFlamelet->LoadFlamelet(CombustionRegime, myInputPtr, FlameletList[0], VarProg, WeightProg, "YES");
  nZPoints = myFlamelet->GetFlameletSize();
  
  // Compute energy perturbation (both positive and negative, saved in flamelet)
  myFlamelet->ComputeEnergyPerturbation(rho0barDeltaEtild);

  TestCoeffApprox(nZPoints, rho0barDeltaEtild, ZMEAN, ZVAR, myFlamelet);
  
  if (myFlamelet != NULL) delete myFlamelet;
 
#else  
  
  PrepareTable();
  
  if (CombustionModel != "ANALYTIC")
  {
      
    // Loop over flamelets
    for (int fl=0; fl<FlameletList.size(); fl++)
    {
      // Load flamelet
      myFlamelet = new Flamelet;
      myFlamelet->LoadFlamelet(CombustionRegime, myInputPtr, FlameletList[fl], VarProg, WeightProg, "NO", DefinitionProg, Zst);
      nZPoints = myFlamelet->GetFlameletSize();
      Preference = myFlamelet->GetFlameletPressure();
      Toxi = myFlamelet->GetFlameletToxi();
      Tfuel = myFlamelet->GetFlameletTfuel();
      double myMax = myFlamelet->GetFlameletVariableMax("PROG");
//    if (myMax <= prog_max)
//    {
//      cerr << "### The list of flamelets is not in ascending order of progress variable! ###" << endl;
//      cerr << "Flamelet " << FlameletList[fl] << " compared to flamelet " << previousFl << endl;
//      throw(-1);
//    }
          
      // Compute energy perturbation (both positive and negative, saved in flamelet)
      myFlamelet->ComputeEnergyPerturbation(rho0barDeltaEtild);
      
      // Loop over Zm and Zv
      for (int im=0; im<nZm; im++)
        for (int iv=0; iv<nZv; iv++)
        {
          double Zvar = Zv[iv];
          if (NormalizedCoordinates == "YES")
          {
            if ((Zm[im] != 0.0)&&(Zm[im] != 1.0))
            {
              Zvar = Zvar*(Zm[im]*(1.0-Zm[im])); // Conversion from Sz to Zvar
            } else {
              Zvar = 0.0;
            }
          }
  
          // For debuggin:
          // cout << "Zm = "<<Zm[im] << ", Zvar = " << Zvar << ", Sv = " << Zv[iv] <<"   ==> "<<myFlamelet->FlameletConvoluteWithPDF(Zm[im],Zvar,"Z") << endl;

          // Compute coefficients for energy, viscosity and thermal diffusion and store them in table
          if (ChemTableType == "COEFF")
          {
            double R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar;
//          CalculateCoefficients(R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar, Zm[im], Zv[iv], myFlamelet);
            CalculateCoefficients(R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar, rho0barDeltaEtild, Zm[im], Zvar, myFlamelet);
            TabCA[im][iv][fl][0] = R_0_tild;
            TabCA[im][iv][fl][1] = T_0_hat;
            TabCA[im][iv][fl][2] = e_0_tild;
            TabCA[im][iv][fl][3] = gamma_0_tild;
            TabCA[im][iv][fl][4] = a_gamma;
            TabCA[im][iv][fl][5] = mu_0_tild;
            TabCA[im][iv][fl][6] = a_mu;
            TabCA[im][iv][fl][7] = lamOcp_0_tild;
            TabCA[im][iv][fl][8] = a_lamOcp;
          }
          // Store in table all species/variables 
          for (int nv=spvar; nv<nVar; nv++)
          {
            TabCA[im][iv][fl][nv] = myFlamelet->FlameletConvoluteWithPDF(Zm[im],Zvar,VarToOutput[nv]);
            if (VarToOutput[nv] == "PROG")
            {
              CA[im][iv][fl] = TabCA[im][iv][fl][nv];  // Coordinate in progress variable
              // Check if mapping is unique (for adaptive only)
              if ((fl > 0) && (ScaleThirdDim == "ADAPTIVE"))
              {
                if ((CA[im][iv][fl] <= CA[im][iv][fl-1]) && (Zv[iv] <= Zm[im] * (1 - Zm[im])) && (im > 0) && (im < nZm-1))
                {
                  cout << "Mapping in the progress variable is not unique !!!!!!!!!!!" << endl;
                  cout << "Flamelet " << FlameletList[fl] << " and " << FlameletList[fl-1] << " at Zm=" << Zm[im] << " and Zv=" << Zv[iv] << endl;
                  cout << CA[im][iv][fl] << " " << CA[im][iv][fl-1] << endl;
                }
              }
            }  
          }
        }
      
      
      if (myFlamelet != NULL) delete    myFlamelet;
      prog_max = myMax;
      previousFl = FlameletList[fl];
    }
  }
  else
    CreateAnalyticalTable();
  
  Table = TabCA;

  if (NormalizedCoordinates == "YES")
  {
    // Determine minimum and maximum value for third dimension for each (Z,Sz) point and then normalize CA
    for (int i=0; i<nZm; i++)
    {
      for (int j=0; j<nZv; j++)
      {
        double min3 = +1.0e10;
        double max3 = -1.0e10;
        for (int k=0; k<nCA; k++)
        {
          if (CA[i][j][k] > max3)
          max3 = CA[i][j][k];
          if (CA[i][j][k] < min3)
            min3 = CA[i][j][k];
        }
        // Normalize CA
        for (int k=0; k<nCA; k++)
        {
          CA[i][j][k] = (CA[i][j][k]-min3)/(max3-min3);
        }
      }
    }
  }
  
  // If table is structured Cartesian, TabCA must be converted to TabCC in the progress variable dimension
  if (ScaleThirdDim != "ADAPTIVE")
  {
    cout << "Creating Cartesian table" << endl;
    CreateCartesianTable();
    Table = TabCC;
  }
  
  OutputTableStatistics();
  WriteTable();
  if (TecplotOutput == "YES")
    WriteTableTecplot();
  CreateChemTableFinalize();
  
#endif
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
