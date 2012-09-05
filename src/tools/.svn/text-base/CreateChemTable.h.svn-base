/**
 * \brief Tool that creates the chemistry table
 *
 * Reads different input flamelet solution files and create chemistry table for steady flamelet model or for FPVA.
 * 
 * Implemented:
 * Structured Cartesian table (all 3 dimensions given by a vector);
 * Structured adaptive table (3rd dimension given by a matrix, i.e., dependent on first 2 other dimensions);
 * 
 * Either species only or species and coefficients for energy, viscosity and thermal diffusion are stored.
 * 
 * Version number is used to ensure that the table was created with the latest version.
 *
 * To be implemented:
 * Artificial Neural Network table (ANN)
 * 
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */

#include "Flamelet.h"

#include <algorithm>

//#define TESTING
//#define INTERP_OLD
 
#define ZMEAN 0.03
#define ZVAR 0.02

#define Ver 2.0

/**********************************************************************************************/
/* Global variables */
double          version;                                                   ///< Version of table creation tool      
int             nZm=0, nZv=0, nCC=0, nCA=0, nC= 0, nVar=0;                 ///< Dimensions of the chemistry table
double          Zst=1.0;                                                   ///< Z value between uniform and stretched mesh
int             spvar=0;                                                   ///< Index of first species variable in VarToOutput
int             CombustionRegime;                                          ///< Combustion regime (1->steady flamelet, 2->fpva)
double          Preference;                                                ///< Reference pressure at which flamelets were computed
double          Toxi, Tfuel;                                               ///< Temperature boundary conditions with which flamelets were computed
Flamelet        *myFlamelet;                                               ///< Flamelet object containing flamelet solution from file
vector <string> VarToOutput, VarProg;
vector <double> WeightProg;                                                ///< Vector of weight values for progress variable
vector <string> FlameletList;                                              ///< List of flamelet files
string          ScaleThirdDim, CombustionModel;
string          ChemTableFilename, ChemTableType;
string          TecplotOutput;                                             ///< "YES" if table also output in tecplot format
string          NormalizedCoordinates;                                     ///< "YES" if tabulated as function of normalized progress variable
                                                                           ///        and mixture fraction segregation factor
string          DefinitionProg;                                            ///< "GENERALIZED" gives c|Zst 
string          WriteSinglePrecision;                                      ///< "YES" 
double          ****Table;                                                 ///< Table data (Var, Zm, Zv, chi/C)
double          ****TabCC;                                                 ///< Table for Cartesian mesh
double          ****TabCA;                                                 ///< Table for adaptive mesh
double          *Zm;                                                       ///< First dimension ZMean
double          *Zv;                                                       ///< Second dimension ZVar
double          *CC;                                                       ///< Third dimension if Cartesian mesh
double          ***CA;                                                     ///< Third dimension if adaptive mesh
double          *BetaPDF;                                                  ///< Beta distribution
double          **WA;                                                      ///< Working array

/**********************************************************************************************/

/* Test if flamelet progress variable is larger or smaller by checking temperature in the filename */
/* Used to order the flamelet in increasing value of progress variable */
bool IsProgSmaller(string fl1, string fl2);

/* Test if flamelet scalar dissipation rate is larger or smaller by checking the filename */
/* Used to order the flamelet in increasing value of scalar dissipation rate */
bool IsChiLarger(string fl1, string fl2);

void CreateChemTableInitialize(ParamMap *myInputPtr);

/**********************************************************************************************/
void PrepareTable();

/**********************************************************************************************/

void CreateCartesianTable();

/**********************************************************************************************/

void ComputeEnergyPerturbation(double &T_p_hat, double &mu_p_tild, double &lamOcp_p_tild, 
                               double Zm, double Zv, Flamelet *myFlamelet, 
                               double R_0_tild, double rho_0_bar, double Delta_e_tild);

/**********************************************************************************************/

void CalculateCoefficients(double &R_0_tild, double &T_0_hat, double &e_0_tild, double &gamma_0_tild, double &a_gamma, 
                           double &mu_0_tild, double &a_mu, double &lamOcp_0_tild, double &a_lamOcp, double &rho_0_bar,
                           double Zm, double Zv, Flamelet *myFlamelet);
/* Compute coefficients for energy, viscosity and thermal diffusion */

/**********************************************************************************************/

void CalculateCoefficients(double &R_0_tild, double &T_0_hat, double &e_0_tild, double &gamma_0_tild, double &a_gamma, 
                           double &mu_0_tild, double &a_mu, double &lamOcp_0_tild, double &a_lamOcp, double &rho_0_bar,
                           double rho0barDeltaEtild, double Zm, double Zv, Flamelet *myFlamelet);
/* Compute coefficients for energy, viscosity and thermal diffusion */

/**********************************************************************************************/

void TestCoeffApprox(double nZ, double rho0barDeltaEtild, double Zm, double Zv, Flamelet *myFlamelet);
/* Test coefficients by comparing with exact solution */

/**********************************************************************************************/

void CreateAnalyticalTable();
/* Create table constructed analytically to test interpolation */

/**********************************************************************************************/

void OutputTableStatistics();
/* Write statistics to standard output */

/**********************************************************************************************/

void WriteTable();
/* Write table to file */

/**********************************************************************************************/

void WriteTableTecplot();
/* Write table to file in Tecplot ASCII format for visualization */

/**********************************************************************************************/

void CreateChemTableFinalize();
/* Clean memory */

/**********************************************************************************************/
