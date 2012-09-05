#include "CreateChemTable_new.h"
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
  if (pos == string::npos) {T1 = -1.0;} // Tst not found ==> quenched solution
  
  pos = fl2.find("Tst");
  Tst = fl2.substr(pos+3);
  istringstream buf2(Tst);
  buf2 >> T2;
  if (pos == string::npos) {T2 = -1.0;} // Tst not found ==> quenched solution
  
  return (T2 > T1);
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
  DefinitionProg        = myInputPtr->getStringParam("PROG_DEFINITION","USUAL");
  

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

  
  nZm = myInputPtr->getIntParam("N_ZMEAN", "100");
  nZv = myInputPtr->getIntParam("N_ZVAR", "20");
  nCC = myInputPtr->getIntParam("N_THIRDDIM", "100");
  
  Zst = myInputPtr->getDoubleParam("Z_ST", "0.03");

  // Read species and weights in progress variable
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

  
  // Add source term for progress variable and progress variable
  i_SRC = spvar;
  VarToOutput.push_back("SRC_PROG");
  VarToOutput.push_back("PROG");
  i_HR = i_SRC +2;
  VarToOutput.push_back("HeatRelease");
  VarToOutput.push_back("rho0");

  
    
  // Read and sort list of flamelet files  
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
  sort(FlameletList.begin(), FlameletList.end(), IsProgSmaller);
 
  //////////////////////////////////////////////////////////////////////////////////////////
  // Summarize inputs to screen
  //////////////////////////////////////////////////////////////////////////////////////////
  
  cout << endl;
  cout << "****************************************************************" << endl;
  cout << "**************** CreateChemTable Initialization ****************" << endl;
  cout << "****************************************************************" << endl;
  cout << endl;
  
  cout << "CHEMTABLE_FILENAME:      " << ChemTableFilename << endl;
  cout << "N_ZMEAN:                 " << nZm << endl;
  cout << "N_ZVAR:                  " << nZv << endl;
  cout << "N_THIRDDIM:              " << nCC << endl;
  cout << "Z_ST:                    " << Zst << endl;
  cout << "PROGRESS VARIABLE:       " << WeightProg[0] << "    " << VarProg[0] << endl;
  for (int i=1; i<VarProg.size(); i++)
    cout << "                         " << WeightProg[i] << "    " << VarProg[i] << endl;
  cout << "STORED VARIABLES:        " << VarToOutput[0] << endl;
  for (int i=1; i<VarToOutput.size(); i++)
    cout << "                         " << VarToOutput[i] << endl;
  cout << "FLAMELET FILES:          " << FlameletList[0] << endl;
  for (int i=1; i<FlameletList.size(); i++)
    cout << "                         " << FlameletList[i] << endl;

  
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
  // Zv is Zv/(z*(1-z)) in fact
  for (int j=0; j<nZv; j++)
    Zv[j] = pow((double) (j) / (double) (nZv-1),2.7);

  
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
  nCA = FlameletList.size();
  getMem3D(&CA, 0, nZm-1, 0, nZv-1, 0, nCA-1, "CreateChemTableInitialize: CA", true);
  // Construct third dimension vector CC
  CC = new double[nCC];
  for (int k=0; k<nCC; k++)
    CC[k] = ((double) k  / (double) (nCC-1));

  //////////////////////////////////////////////////////////////////////////////////////////
  // Tables
  //////////////////////////////////////////////////////////////////////////////////////////    
  getMem4D(&TabCA, 0, nZm-1, 0, nZv-1, 0, nCA-1, 0, nVar-1, "CreateChemTableInitialize: adaptive table", true);
    
  return;
}

/**********************************************************************************************/

void CreateCartesianTable()
/* Use adaptive table TabCA to create Cartesian table */
{
  Table.Allocate(nZm,nZv,nCC,nVar);
  Table.SetTableFileName(ChemTableFilename);
  Table.SetReferencePressure(Preference);
  for (int i=0; i<nVar; ++i) Table.AddVarName( VarToOutput[i] );
  for (int i=0; i<nZm;  ++i) Table.SetCoordinate1(i, Zm[i] );
  for (int i=0; i<nZv;  ++i) Table.SetCoordinate2(i, Zv[i] );
  for (int i=0; i<nCC;  ++i) Table.SetCoordinate3(i, CC[i] );
  
  // Loop over the three mapping directions in Cartesian table
  for (int i=0; i<nZm; i++)
    for (int j=0; j<nZv; j++)
      for (int k=0; k<nCC; k++)
      {
        double alpha_up;
        double alpha_down;
        int kk_up;
        int kk_down;

        BinarySearch(kk_down, alpha_down, CA[i][j], 0, nCA-1, CC[k]);
        kk_up = kk_down + 1;
        alpha_up = 1.0 - alpha_down;
        if ((CA[i][j][kk_down]>CC[k])||(CC[k]>CA[i][j][kk_up])) {
        	cout << "Failed "<<CC[k] << " in [" << CA[i][j][kk_down] << " " << CA[i][j][kk_up] << "]" << endl;
        }

        for (int n=0; n<nVar; n++) {
          double val = alpha_up   * TabCA[i][j][kk_up][n]
                     + alpha_down * TabCA[i][j][kk_down][n];
          Table.SetTableValue(i,j,k,n,val);
        }
      }

  return;
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
  R_0_tild      = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"R0");
  rho_0_bar     = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"rho0");
  RT_0_tild     = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"RT0");
  e_0_tild      = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"E0");
  mu_0_tild     = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"mu0");
  lamOcp_0_tild = myFlamelet->FlameletConvoluteWithPDF_new(Zm,Zv,"lamOcp0");
  T_0_hat       = RT_0_tild / R_0_tild;

  // Compute quantities if energy perturbed by +/- rho_0_bar * Delta_E_tild
  T_p_hat       = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "RTpl") / R_0_tild;
  mu_p_tild     = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "mupl");
  lamOcp_p_tild = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "lamOcppl");  
  T_m_hat       = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "RTmi") / R_0_tild;
  mu_m_tild     = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "mumi");
  lamOcp_m_tild = myFlamelet->FlameletConvoluteWithPDF_new(Zm, Zv, "lamOcpmi");  

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

int main(int argc, char *argv[])
{
  char InputFileName[kStrLong];
  sprintf(InputFileName, "CreateChemTable_new.in");
  for (int i=1; i<argc; i++)
    strcpy(InputFileName, argv[i]);
  
  
  ParamMap myInput(InputFileName);
  ParamMap *myInputPtr = &myInput.getThisParamMap();

  int    nZPoints;                       // Number of points in flamelet
  double rho0barDeltaEtild = 5000.0;     // Overall Delta E * rho0_bar for all flamelets
  
  CreateChemTableInitialize(myInputPtr);
  
  PrepareTable();
	  
  ofstream fout1;
  ofstream fout2;
  fout1.open("out1.dat");
  fout2.open("out2.dat");

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
	  // Compute energy perturbation (both positive and negative, saved in flamelet)
	  myFlamelet->ComputeEnergyPerturbation(rho0barDeltaEtild);
	  
	  // Loop over Zm and Zv
	  for (int im=0; im<nZm; im++)
	    for (int iv=0; iv<nZv; iv++)
	    {
	      double Zvar = Zv[iv]*(Zm[im]*(1.0-Zm[im])); // Conversion from Sz to Zvar
	      // Compute coefficients for energy, viscosity and thermal diffusion and store them in table
	      double R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar;
	      CalculateCoefficients(R_0_tild, T_0_hat, e_0_tild, gamma_0_tild, a_gamma, mu_0_tild, a_mu, lamOcp_0_tild, a_lamOcp, rho_0_bar, rho0barDeltaEtild, Zm[im], Zvar, myFlamelet);
              if (iv == 0) {
                fout1 << Zm[im] << " " << e_0_tild << endl;
              } else if (iv == 1) {
                fout2 << Zm[im] << " " << e_0_tild << endl;
              }
	      TabCA[im][iv][fl][0] = R_0_tild;
	      TabCA[im][iv][fl][1] = T_0_hat;
	      TabCA[im][iv][fl][2] = e_0_tild;
	      TabCA[im][iv][fl][3] = gamma_0_tild;
	      TabCA[im][iv][fl][4] = a_gamma;
	      TabCA[im][iv][fl][5] = mu_0_tild;
	      TabCA[im][iv][fl][6] = a_mu;
	      TabCA[im][iv][fl][7] = lamOcp_0_tild;
	      TabCA[im][iv][fl][8] = a_lamOcp;
	      // Store in table all species/variables 
	      for (int nv=spvar; nv<nVar; nv++)
	      {
	        TabCA[im][iv][fl][nv] = myFlamelet->FlameletConvoluteWithPDF_new(Zm[im],Zvar,VarToOutput[nv]);
	      }
	    }
	  if (myFlamelet != NULL) delete    myFlamelet;
          //throw(-1);
	}

  
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
        if (max3 != min3) {
        	CA[i][j][k] = (CA[i][j][k]-min3)/(max3-min3);
        } else {
        	CA[i][j][k] = ((double) k  / (double) (nCA-1));
        }
      }
      // Make CA well behaved
      int threshold =0;
      for (int k=0; k<nCA; k++)
      {
      	if (CA[i][j][k] == 1.0) threshold = 1;
      	if (threshold == 1) CA[i][j][k] == 1.0;
      }
      // Check
      for (int k=1; k<nCA; k++)
      {
      	if (CA[i][j][k] < CA[i][j][k-1]) {
      	  cout << " wait ! " << CA[i][j][k-1] <<" > "<< CA[i][j][k]<< endl;
      	  cout << "Zm="<<Zm[i] << " Zv="<<Zv[j]<<endl;
      	  throw(-1);
      	}
      }
      // Set src_prog to zero in C=0 and 1
      for (int k=0; k<nCA; k++) {
        if ( (CA[i][j][k] == 0.0)||(CA[i][j][k] == 1.0)) {
          TabCA[i][j][k][i_SRC] = 0.0;
          TabCA[i][j][k][i_HR]  = 0.0;
        }
      }
    }
  } 

  CreateCartesianTable();
  Table.WriteTable("no_debug");
  Table.print();
  Table.CheckNan();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
