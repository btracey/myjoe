/*
 * CombModel_VariableProperties.h
 *
 *  Created on: July 16, 2009
 *      Author: Vincent Terrapon
 *      Update: Jan 11, 2010
 */

#ifndef COMBMODEL_VARIABLEPROPERTIES_H_
#define COMBMODEL_VARIABLEPROPERTIES_H_

#define MIXING_VISCOSITY

#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                               RansCombVariableProperties                            ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Mixing model based on mixture fraction. 
 *
 * The species mass fractions Species::Yk are fixed. There is neither combustion, nor mixing. But
 * the different thermodynamics properties (e.g., cp, gamma) can vary with temperature.
 * 
 */
class RansCombVarProperties: public RansCombBase {
public:
  // member variables
  Mixture myMixture;                    ///< Mixture containing all the species information.
  double  *muLam;                       ///< Laminar viscosity at cell center.
  double  *LambdaOverCp;                ///< Thermal diffusion coefficient at cell center.

public:
  // constructors

  RansCombVarProperties()
  {
    if (mpi_rank == 0)
      cout << "RansCombVarProperties()" << endl;
    
    // Other cell-centered quantities
    muLam        = NULL;        registerScalar(muLam, "muLam" , CV_DATA);
    LambdaOverCp = NULL;        registerScalar(LambdaOverCp, "LambdaOverCp", CV_DATA);
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombVarProperties()
  {
    myMixture.Unload();
    if (mpi_rank == 0)
      cout << endl << "***** Variable Properties finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    int CombustionRegime = -1;
    myMixture.Load(CombustionRegime, &getThisParamMap());
    myMixture.GetSpeciesMassFraction(0.0);
    
    for (int icv = 0; icv < ncv; icv++)
      RoM[icv] = myMixture.GetMixR_o_M();
    updateCvData(RoM, REPLACE_DATA);
   
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

      // Compute mixture temperature from energy and species mass fractions
      temp[icv] = myMixture.ComputeMixTemperature(E, temp[icv]);

      // Compute mixture gas constant, enthalpy and pressure
      enthalpy[icv] = myMixture.ComputeMixEnthalpy(temp[icv]);
      press[icv] = rho[icv] * RoM[icv] * temp[icv];

      if (isnan(press[icv]) || (press[icv]<= 0.0) || (rho[icv]<= 0.0))
      {
        cout << "WARNING! :" << endl;
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1]  << " / " << x_cv[icv][2] << endl;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / (rho[icv] * rho[icv]) << endl;
      }

      // Compute mixture viscosity, heat conductivity, gamma
      cp = myMixture.ComputeMixCp(temp[icv]);
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);

      myMixture.ComputeMixLambdaMu(temp[icv]);
      muLam[icv] = myMixture.GetMixMul();
      LambdaOverCp[icv] = myMixture.GetMixLambda() / cp;
    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
    
    updateCvData(muLam, REPLACE_DATA);
    updateCvData(LambdaOverCp, REPLACE_DATA);
  }
  
#ifdef MIXING_VISCOSITY  
  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *  
   */
  virtual void calcMaterialProperties()
  {
    if (mu_ref > 0.0)
    {
      // internal faces
      for (int ifa = nfa_b; ifa < nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];
  
        double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));
 
        mul_fa[ifa] = (w1*muLam[icv0] + w0*muLam[icv1])/(w0+w1);
        lamOcp_fa[ifa] = (w1*LambdaOverCp[icv0] + w0*LambdaOverCp[icv1])/(w0+w1);
      }    
      
      // boundary faces computed in setBC
    }
  }
#endif
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction();
    R = myMixture.GetMixR_o_M();
    h = myMixture.ComputeMixEnthalpy(T);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - R);
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction();
    R = myMixture.GetMixR_o_M();
    T = p / (rho * R);
    h = myMixture.ComputeMixEnthalpy(T);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - R);
  }
  
  // \brief Compute for a given energy, density and scalars: pressure, temperature, gas constant and ratio of specific heat
  virtual void calcThermoProp_e(double &T, double &p, double &R, double &gam, double &rho, double &e, double *Scal, int nScal)
  {
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction();
    R = myMixture.GetMixR_o_M();
    T = myMixture.ComputeMixTemperature(e, 1000.0);
    p = rho * R * T;
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - R);
  }
 
  virtual double calcMuLam(int icv)
  {
    return muLam[icv];
  }
  
  virtual double calcMuLam(double temp)
  {
    cout << "### RansCombVariableProperties::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
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
  void ComputeProperties_T(double &h, double &RoM, double &gam, double &mu, double &lamOcp, double T)
  {
    RoM = myMixture.GetMixR_o_M();
    h = myMixture.ComputeMixEnthalpy(T);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - RoM);
    myMixture.ComputeMixLambdaMu(T);
    mu = myMixture.GetMixMul();
    lamOcp = myMixture.GetMixLambda() / cp;
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
  void ComputeProperties_H(double &T, double &RoM, double &gam, double &mu, double &lamOcp, double h)
  {
    RoM = myMixture.GetMixR_o_M();
    T = myMixture.ComputeMixTemperature_H(h);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - RoM);
    myMixture.ComputeMixLambdaMu(T);
    mu = myMixture.GetMixMul();
    lamOcp = myMixture.GetMixLambda() / cp;
  }
};

#endif /* COMBMODEL_VARIABLEPROPERTIES_H_ */
