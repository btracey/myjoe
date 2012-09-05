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

#include "Flamelet_new.h"
#include "ChemtableNormalized.h"

#include <algorithm>

#define Ver 3.0

/**********************************************************************************************/
/* Global variables */
double          version;                                                   ///< Version of table creation tool      
int             nZm=0, nZv=0, nCC=0, nCA=0, nVar=0;                 ///< Dimensions of the chemistry table
double          Zst=1.0;                                                   ///< Z value between uniform and stretched mesh
int             spvar=0, i_SRC=0, i_HR = 0;                                                  ///< Index of first species variable in VarToOutput
int             CombustionRegime;                                          ///< Combustion regime (1->steady flamelet, 2->fpva)
double          Preference;                                                ///< Reference pressure at which flamelets were computed
double          Toxi, Tfuel;                                               ///< Temperature boundary conditions with which flamelets were computed
Flamelet        *myFlamelet;                                               ///< Flamelet object containing flamelet solution from file
vector <string> VarToOutput, VarProg;
vector <double> WeightProg;                                                ///< Vector of weight values for progress variable
vector <string> FlameletList;                                              ///< List of flamelet files
string          ChemTableFilename;
string          DefinitionProg;                                            ///< "GENERALIZED" gives c|Zst 
ChemtableNormalized Table;
double          ****TabCA;                                                 ///< Table for adaptive mesh
double          *Zm;                                                       ///< First dimension ZMean
double          *Zv;                                                       ///< Second dimension ZVar
double          *CC;                                                       ///< Third dimension if Cartesian mesh
double          ***CA;                                                     ///< Third dimension if adaptive mesh

/**********************************************************************************************/

/* Test if flamelet progress variable is larger or smaller by checking temperature in the filename */
/* Used to order the flamelet in increasing value of progress variable */
bool IsProgSmaller(string fl1, string fl2);

void CreateChemTableInitialize(ParamMap *myInputPtr);

/**********************************************************************************************/
void PrepareTable();

/**********************************************************************************************/

void CreateCartesianTable();

/**********************************************************************************************/

void CalculateCoefficients(double &R_0_tild, double &T_0_hat, double &e_0_tild, double &gamma_0_tild, double &a_gamma, 
                           double &mu_0_tild, double &a_mu, double &lamOcp_0_tild, double &a_lamOcp, double &rho_0_bar,
                           double rho0barDeltaEtild, double Zm, double Zv, Flamelet *myFlamelet);
/* Compute coefficients for energy, viscosity and thermal diffusion */


/**********************************************************************************************/
