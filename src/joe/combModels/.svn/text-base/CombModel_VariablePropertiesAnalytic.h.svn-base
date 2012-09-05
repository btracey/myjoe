/*
 * CombModel_VariablePropertiesAnalytic.h
 *
 *  Created on: September 08, 2010
 *      Author: Vincent Terrapon
 *      Update: 
 */

#ifndef COMBMODEL_VARIABLEPROPERTIESANALYTIC_H_
#define COMBMODEL_VARIABLEPROPERTIESANALYTIC_H_

#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                               RansCombVariablePropertiesAnalytic                    ----
// ------                                                                                 ------
// #############################################################################################

class RansCombVarPropertiesAnalytic: public RansCombBase {
public:
  // member variables
  double cp_1;
  double cp_2;
  double cp_0;
  double dcpdT;
  double h_0;
  double cv_0;
  
public:
  // constructors

  RansCombVarPropertiesAnalytic()
  {
    if (mpi_rank == 0)
      cout << "RansCombVarPropertiesAnalytic()" << endl;
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombVarPropertiesAnalytic()
  {
    if (mpi_rank == 0)
      cout << endl << "***** Variable Properties finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    double T1 = 300.0;
    double T2 = 1000.0;
    
    double R = getDoubleParam("GAS_CONSTANT", "288.1");
    for (int icv = 0; icv < ncv; icv++)
      RoM[icv] = R;
    updateCvData(RoM, REPLACE_DATA);
    
    
    double gamma_1 = getDoubleParam("GAMMA_AT_300",  "1.4");
    double gamma_2 = getDoubleParam("GAMMA_AT_1000", "1.25");

    // cp(T) = cp1 + (cp2 - cp1)/(T2 - T1) * (T - T1) := cp0 + dcpdT * T
    // cv(T) = cp(T) - R = cp0 - R + dcpdT * T := cv0 + dcpdT * T
    // h(T) = 0.5 * dcpdT * T^2 + cp0 * T + h0
    // e(T) = 0.5 * dcpdT * T^2 + cv0 * T + h0
    cp_1 = gamma_1 * R / (gamma_1 - 1.0);
    cp_2 = gamma_2 * R / (gamma_2 - 1.0);
    dcpdT = (cp_2 - cp_1) / (T2 - T1);
    cp_0 = cp_1 - dcpdT * T1;
    cv_0 = cp_0 - R;
    h_0 = - 0.5 * dcpdT * T1 * T1 - cp_0 * T1;
    
    if (mpi_rank == 0)
      cout << endl << "***** Variable Properties initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties.
   *
   *  Compute for the given composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double E, cp;
    double kinecv = 0.0;

    for (int icv = 0; icv < ncv; icv++)
    {
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      temp[icv] = (-cv_0 + sqrt(cv_0 * cv_0 - 2.0 * dcpdT * (h_0 - E))) / dcpdT;

      // Compute mixture gas constant, enthalpy and pressure
      enthalpy[icv] = E - RoM[icv] * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];

      // Compute mixture viscosity, heat conductivity, gamma
      cp = cp_0 + dcpdT * temp[icv];
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);
    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
  }
    
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    R = RoM[0];
    h = 0.5 * dcpdT * T * T + cp_0 * T + h_0;
    double cp = cp_0 + dcpdT * T;
    gam = cp / (cp - R);
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    R = RoM[0];
    T = p / (rho * R);
    h = 0.5 * dcpdT * T * T + cp_0 * T + h_0;
    double cp = cp_0 + dcpdT * T;
    gam = cp / (cp - R);
  }
  
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_e(double &T, double &p, double &R, double &gam, double &rho, double &e, double *Scal, int nScal)
  {
    R = RoM[0];
    T = (-cv_0 + sqrt(cv_0 * cv_0 - 2.0 * dcpdT * (h_0 - e))) / dcpdT;
    p = rho * R * T;
    double cp = cp_0 + dcpdT * T;
    gam = cp / (cp - R);
  }
 
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_T(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa]);
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_H(T_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],h_bfa[ifa]);
  }
  
  /*! \brief Compute for a given temperature the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  T Temperature in [K].
   */
  void ComputeProperties_T(double &h, double &R, double &gam, double &mu, double &lamOcp, double T)
  {
    R = RoM[0];
    h = 0.5 * dcpdT * T * T + cp_0 * T + h_0;
    double cp = cp_0 + dcpdT * T;
    gam = cp / (cp - R);
    if (mu_ref > 0.0)
    {
      mu = calcMuLam(T);
      lamOcp = mu / Pr;
    }
  }

  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] T Temperature in [K].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  H Enthalpy (chemical + sensible) in [J/kg].
   */
  void ComputeProperties_H(double &T, double &R, double &gam, double &mu, double &lamOcp, double h)
  {
    R = RoM[0];
    T = (-cp_0 + sqrt(cp_0 * cp_0 - 2.0 * dcpdT * (h_0 - h))) / dcpdT;
    double cp = cp_0 + dcpdT * T;
    gam = cp / (cp - R);
    if (mu_ref > 0.0)
    {
      mu = calcMuLam(T);
      lamOcp = mu / Pr;
    }
  }
};

#endif /* COMBMODEL_VARIABLEPROPERTIES_H_ */
