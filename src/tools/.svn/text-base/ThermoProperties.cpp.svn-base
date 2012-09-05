/**
 * \brief Tool to compute different thermo properties.
 *
 * Use thermo binary input file and the functions in combustion.h.
 * \author Vincent Terrapon 
 * \date August 2010
 * \version 1.0
 */

#include "Combustion.h"

Mixture   myMixture;            ///< Mixture containing the different species.

// Overall properties
double    rho;                  ///< Density.
double    press;                ///< Pressure.
double    temp;                 ///< Temperature.
double    M;                    ///< Molar mass.
double    RoM;                  ///< Specific gas constant R/M.
double    vel;                  ///< Velocity.
double    cp;                   ///< Specific heat capacity at constant pressure.
double    cv;                   ///< Specific heat capacity at constant volume.
double    enthalpy;             ///< Specific enthalpy.
double    tot_enthalpy;         ///< Specific total enthalpy.
double    energy;               ///< Specific energy.
double    tot_energy;           ///< Specific total energy.
double    tot_press;            ///< Total pressure.
double    tot_temp;             ///< Total temperature.
double    gama;                 ///< Heat capacity ratio cp/cv.
double    a_sound;              ///< Sound speed.
double    Mach;                 ///< Mach number.
double    mu;                   ///< Laminar viscosity.
double    lambda;               ///< Thermal conductivity.
double    lamOcp;               ///< Thermal conductivity divided by heat capacity.

bool      p = false;
bool      T = false;
bool      H = false;
bool      r = false;


void computeBasicProperties_PT(double p, double T)
{
  rho = p / RoM / T;
  enthalpy = myMixture.ComputeMixEnthalpy(T);
}

void computeBasicProperties_PRHO(double p, double r)
{
  temp = p / RoM / r;
  enthalpy = myMixture.ComputeMixEnthalpy(temp);
}

void computeBasicProperties_PH(double p, double H)
{
  temp = myMixture.ComputeMixTemperature_H(H);
  rho = p / RoM / temp;
}

void computeBasicProperties_RHOT(double r, double T)
{
  press = r * RoM * T;
  enthalpy = myMixture.ComputeMixEnthalpy(T);
}

void computeBasicProperties_RHOH(double r, double H)
{
  temp = myMixture.ComputeMixTemperature_H(H);
  press = r * RoM * temp;
}

void computeOtherProperties()
{
  energy = myMixture.ComputeMixEnergy(T);

  cp = myMixture.ComputeMixCp(temp);
  cv = cp - RoM;
  gama = cp / cv;
  
  tot_enthalpy = enthalpy + 0.5 * vel * vel;
  tot_energy = energy + 0.5 * vel * vel;
  tot_temp = myMixture.ComputeMixTemperature_H(tot_enthalpy, temp);
  tot_press = press / pow(temp / tot_temp, gama / (gama - 1.0));
  
  a_sound = sqrt(gama * RoM * temp);
  Mach = vel / a_sound;
  
  myMixture.ComputeMixLambdaMu(temp);
  mu = myMixture.GetMixMul();
  lambda = myMixture.GetMixLambda();
  lamOcp = lambda / cp;
}

void outputResults()
{
  printf("\n");
  printf("THERMODYNAMIC PROPERTIES\n");
  printf("************************\n\n");
  
  printf("Pressure: \t\t\t\t\t%12.4e [Pa]\n", press);
  printf("Temperature: \t\t\t\t\t%12.4e [K]\n", temp);
  printf("Density]: \t\t\t\t\t%12.4e [kg/m^3\n", rho);
  printf("Velocity: \t\t\t\t\t%12.4e [m/s]\n", vel);
  
  printf("Molar mass: \t\t\t\t\t%12.4e [kg/kmol]\n", M);
  printf("Specific gas constant: \t\t\t\t%12.4e [J/kg/K]\n", RoM);
  
  printf("Specific heat capacity at constant pressure: \t%12.4e [J/kg/K]\n", cp);
  printf("Specific heat capacity at constant volume: \t%12.4e [J/kg/K]\n", cv);
  printf("Heat capacity ratio: \t\t\t\t%12.4e [-]\n", gama);
  printf("Specific enthalpy: \t\t\t\t%12.4e [J/kg]\n", enthalpy);
  printf("Specific energy: \t\t\t\t%12.4e [J/kg]\n", energy);
  
  printf("Specific total enthalpy: \t\t\t%12.4e [J/kg]\n", tot_enthalpy);
  printf("Specific total energy: \t\t\t\t%12.4e [J/kg]\n", tot_energy);
  printf("Total temperature: \t\t\t\t%12.4e [K]\n", tot_temp);
  printf("Total pressure: \t\t\t\t%12.4e [Pa]\n", tot_press);
  
  printf("Sound speed: \t\t\t\t\t%12.4e [m/s]\n", a_sound);
  printf("Mach number: \t\t\t\t\t%12.4e [-]\n", Mach);
  
  printf("Laminar dynamic viscosity: \t\t\t%12.4e [kg/m/s]\n", mu);
  printf("Laminar heat conductivity: \t\t\t%12.4e [J/K/m/s]\n", lambda);
  printf("Ratio conductivity to heat capacity: \t\t%12.4e [kg/m/s]\n", lamOcp);
  
  printf("\n");
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  initMpiStuff();
  
  char InputFileName[kStrLong];
  sprintf(InputFileName, "ThermoProperties.in");
  for (int i=1; i<argc; i++)
    strcpy(InputFileName, argv[i]);
  
  
  ParamMap myInput(InputFileName);
  ParamMap *myInputPtr = &myInput.getThisParamMap();
  
  // Load mixture
  int CombustionRegime = -1;
  myMixture.Load(CombustionRegime, myInputPtr, "YES");
  myMixture.GetSpeciesMassFraction(0.0);
  
  // Read input data
  vel   = myInputPtr->getDoubleParam("VELOCITY", 0.0);

  if (myInputPtr->checkParam("PRESSURE"))
  {
    p = true;
    press = myInputPtr->getDoubleParam("PRESSURE");
  }
  
  if (myInputPtr->checkParam("TEMPERATURE"))
  {
    T = true;
    temp = myInputPtr->getDoubleParam("TEMPERATURE");
  }
  
  if (myInputPtr->checkParam("DENSITY"))
  {
    r = true;
    rho = myInputPtr->getDoubleParam("DENSITY");
  }
  
  if (myInputPtr->checkParam("ENTHALPY"))
  {
    H = true;
    enthalpy = myInputPtr->getDoubleParam("ENTHALPY");
  }
    
  // Compute properties
  M = myMixture.GetMixM();
  RoM = myMixture.GetMixR_o_M();
  
  if (p && T)
    computeBasicProperties_PT(press, temp);
  else if (p && H)
    computeBasicProperties_PH(press, enthalpy);
  else if (p && r)
    computeBasicProperties_PRHO(press, rho);
  else if (r && T)
    computeBasicProperties_RHOT(rho, temp);
  else if (r && H)
    computeBasicProperties_RHOH(rho, enthalpy);
  else
  {
    cerr << "### Missing input data! ###" << endl;
    throw(-1);
  }
  
  computeOtherProperties();
  
  // Output results
  outputResults();
  
  return 0;
}
