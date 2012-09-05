#ifndef ATMOSPERICCOND_H
#define ATMOSPERICCOND_H


#include <math.h>

/**
 * Standard Atmosphere Computations: 
 * computes properties related to the 1976 standard atmosphere up to 84852 meters
 * 
 * input: 
 * \param height height in [m]
 * \param g Acceleration of gravity  in (m/s/s)
 * \param gamma Ratio of specific heats
 * \param R Gas constant for air in (J/kg-K)
 * 
 * output: 
 * \param temp temperature 
 * \param press pressure 
 * \param rho density 
 * \param sound sound speed 
 */
int getAthmospericConditions(double &temp, double &press, double &rho, 
    double height, double g, double gamma, double R);

void IsoThermal(double &TNew, double &PNew, double &RhoNew, 
    double Z0, double Z1, double T0, double P0, double Rho0, double g, double gamma, double R);

void Gradient(double &TNew, double &PNew, double &RhoNew, 
    double Z0, double Z1, double T0, double P0, double Rho0, double Lapse, double g, double gamma, double R);









inline int getAthmospericConditions(double &temp, double &press, double &rho, 
    double height, double g, double gamma, double R)
{
  // double g = 9.81;       // Acceleration of gravity (m/s/s)
  // double gamma = 1.4;    // Ratio of specific heats
  // double R = 287.0;      // Gas constant for air (J/kg-K)
  
  // Altitudes (m)
  double Start = 0.0;
  double H1 = 11000.0;
  double H2 = 20000.0;
  double H3 = 32000.0;
  double H4 = 47000.0;
  double H5 = 51000.0;
  double H6 = 71000.0;
  double H7 = 84852.0;
  
  //Lapse Rates (K/m)
  double L1 = -0.0065;
  double L2 =  0.0;
  double L3 =  0.001;
  double L4 =  0.0028;
  double L5 =  0.0;
  double L6 = -0.0028;
  double L7 = -0.002;

  // Initial Values
  double T0 = 288.16;     // (k)
  double P0 = 1.01325e5;  // (pa)
  double Rho0 = 1.225;    // (kg/m^3)

  double TNew, PNew, RhoNew;

  if (height > Start && height <= H1)
  {
    Gradient(TNew, PNew, RhoNew, Start, height, T0, P0, Rho0, L1, g, gamma, R);
  }
  else if (height > H1 && height <= H2)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, height, TNew, PNew, RhoNew, g, gamma, R);
  }
  else if (height > H2 && height <= H3)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, H2, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H2, height, TNew, PNew, RhoNew, L3, g, gamma, R);
  }
  else if (height > H3 && height <= H4)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, H2, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H2, H3, TNew, PNew, RhoNew, L3, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H3, height, TNew, PNew, RhoNew, L4, g, gamma, R);
  }
  else if (height > H4 && height <= H5)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, H2, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H2, H3, TNew, PNew, RhoNew, L3, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H3, H4, TNew, PNew, RhoNew, L4, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H4, height, TNew, PNew, RhoNew, g, gamma, R);
  }
  else if (height > H5 && height <= H6)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, H2, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H2, H3, TNew, PNew, RhoNew, L3, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H3, H4, TNew, PNew, RhoNew, L4, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H4, H5, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H5, height, TNew, PNew, RhoNew, L6, g, gamma, R);
  }
  else if (height > H6 && height <= H7)
  {
    Gradient(TNew, PNew, RhoNew, Start, H1, T0, P0, Rho0, L1, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H1, H2, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H2, H3, TNew, PNew, RhoNew, L3, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H3, H4, TNew, PNew, RhoNew, L4, g, gamma, R);
    IsoThermal(TNew, PNew, RhoNew, H4, H5, TNew, PNew, RhoNew, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H5, H6, TNew, PNew, RhoNew, L6, g, gamma, R);
    Gradient(TNew, PNew, RhoNew, H6, height, TNew, PNew, RhoNew, L7, g, gamma, R);
  }
  else
    return 0; // retrun -1 for error

  temp = TNew;
  press = PNew;
  rho = RhoNew;

  return 1; // return 0 for success 
}


inline void Gradient(double &TNew, double &PNew, double &RhoNew, 
    double Z0, double Z1, double T0, double P0, double Rho0, double Lapse, double g, double gamma, double R)
{
  TNew = T0 + Lapse*(Z1 - Z0);
  PNew = P0*pow(TNew/T0, -g/(Lapse*R));
  RhoNew = PNew/(R*TNew);
}

inline void IsoThermal(double &TNew, double &PNew, double &RhoNew, 
    double Z0, double Z1, double T0, double P0, double Rho0, double g, double gamma, double R)
{
  TNew = T0;
  PNew = P0*exp(-(g/(R*TNew))*(Z1-Z0));
  RhoNew = PNew/(R*TNew);
}







#endif /* ATMOSPERICCOND_H_ */
