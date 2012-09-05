/*
 * CombModel_BinaryMixing.h
 *
 *  Created on: May 6, 2009
 *      Author: Vincent Terrapon
 *      Update: Jan 11, 2010
 */

#ifndef COMBMODEL_BINARYMIXING_H_
#define COMBMODEL_BINARYMIXING_H_

// Analytical deviation from linear mixing rules for single analytical flamelet
//#define ROM_DEF    R   = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
//#define CP_GAM_DEF cp  = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) ); gam = cp / (cp - R);
//#define H_DEF      h   = cp * T;
//#define T_DEF      T   = h / cp;      
//#define T_E_DEF    T   = E / (cp - R);

#define ROM_DEF    R   = zm * RoM_1 + (1.0 - zm) * RoM_2;
#define CP_GAM_DEF cp  = zm * Cp_1  + (1.0 - zm) * Cp_2;   gam = cp / (cp - R);
#define H_DEF      h   = gam * R * (T - Treference) / (gam - 1.0);
#define T_DEF      T   = (gam - 1.0) * h / (gam * R) + Treference;
#define T_E_DEF    T   = (gam - 1.0) * E / R + gam * Treference;

#include "UgpWithCvCompFlow.h"
#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                         RansCombBinaryMixing                                        ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Binary Mixing model based on mixture fraction. 
 *
 *  The species mass fractions Species::Yk are functions of the mean scalar mixture fraction \p Zmean
 *  only through specific mixing rules. It can represent pure mixing (linear mixing rule) or combustion
 *  (non-linear rule). Viscosity and thermal diffusion use simple temperature dependence rule. 
 *  Specific heat coefficients are not a direct function of the temperature but of the mixture fraction.
 *  Only 1 flamelet can be represented.
 */
//class RansCombBinaryMixing : virtual public UgpWithCvCompFlow
class RansCombBinaryMixing: public RansCombBase
{
public:
  // member variables

  double  *ZMean;                       ///< Mixture fraction Z.
//  double  *rhoDkZ;                      ///< Diffusion coefficient of mixture fraction Z.
  double  Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  double  Gamma_1;                      ///< Gamma of gas 1 (Z=1).
  double  Gamma_2;                      ///< Gamma of gas 2 (Z=0).
  double  MolMass_1;                    ///< Molecular mass of gas 1 (Z=1).
  double  MolMass_2;                    ///< Molecular mass of gas 2 (Z=0).
  double  RoM_1;                        ///< Gas constant of gas 1 (Z=1).
  double  RoM_2;                        ///< Gas constant of gas 2 (Z=0).
  double  Cp_1;                         ///< Constant specific heat capacity of gas 1 (Z=1).
  double  Cp_2;                         ///< Constant specific heat capacity of gas 2 (Z=0).
  double  rho_1;                        ///< Density of gas 1 (Z=1).
  double  rho_2;                        ///< Density of gas 2 (Z=0).
  int     ZMean_Index;                  ///< Index in scalar vector of ZMean.
  double  Treference;                   ///< Reference temperature for which enthalpy is zero ( H(Treference) = 0 ).
  
  double  coeff_1;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_2;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_3;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_4;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  dcp_0;                        ///< Coefficient to describe amplitude of the cp deviation in non-linear mixing rule
  double  dRoM_0;                       ///< Coefficient to describe amplitude of the RoM deviation in non-linear mixing rule

public:
  // constructors

  RansCombBinaryMixing()
  {
    if (mpi_rank == 0)
      cout << "RansCombBinaryMixing()" << endl;

    // Overall constants
    Schmidt_turb_Z = getDoubleParam("SCTURB_ZMean", "1.0");
    if (mpi_rank == 0) 
      cout << "Turbulent Schmidt number for Z Sc_t=" << Schmidt_turb_Z << endl;
  
    // Scalar transport equation for ZMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-4;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombBinaryMixing()
  {
    if (scalarTranspEqVector[ZMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZMean_Index].dpress_dphi; 
      scalarTranspEqVector[ZMean_Index].dpress_dphi = NULL;
    }

    if (mpi_rank == 0)
      cout << endl << "***** Binary Mixing finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    MolMass_1 = getDoubleParam("MOLMASS_1", "4.0026");
    MolMass_2 = getDoubleParam("MOLMASS_2", "2.016");
    
    Gamma_1 = getDoubleParam("GAMMA_1", "1.66");
    Gamma_2 = getDoubleParam("GAMMA_2", "1.41");
    
    RoM_1 = 8314.4 / MolMass_1;
    RoM_2 = 8314.4 / MolMass_2;
    
    Cp_1 = Gamma_1 * RoM_1 / (Gamma_1 - 1.0);
    Cp_2 = Gamma_2 * RoM_2 / (Gamma_2 - 1.0);   

//    // Non-linear mixnig rules
//    coeff_1 = getDoubleParam("NON_LIN_MIX_COEFF_1", "50.0");
//    coeff_2 = getDoubleParam("NON_LIN_MIX_COEFF_2", "10.0");
//    coeff_3 = getDoubleParam("NON_LIN_MIX_COEFF_3", "10.0");
//    coeff_4 = getDoubleParam("NON_LIN_MIX_COEFF_4", "10.0");
//    
//    // Default values correspond to linear mixing rule
//    dcp_0  = getDoubleParam("NON_LIN_MIX_D_CP", "0.0");
//    dRoM_0 = getDoubleParam("NON_LIN_MIX_D_R", "0.0");
    
    Treference = getDoubleParam("T_REFERENCE", "0.0");
  
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");               ZMean = eq->phi;               ZMean_Index = getScalarTransportIndex("ZMean");
    
    // For coupled solution
    if ((getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_COUPLED") || (getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_SEMICOUPLED")
        || (getStringParam("TIME_INTEGRATION") == "BDF2_SEMICOUPLED"))
      scalarTranspEqVector[ZMean_Index].dpress_dphi = new double[ncv_g];
      
    if (mpi_rank == 0)
      cout << endl << "***** Binary Mixing initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double zm, R, E, cp, gam, h, T;
    double kinecv = 0.0;

    for (int icv = 0; icv < ncv; icv++)
    {
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      // Compute mixture gas constant and specific heat coefficient
      zm = ZMean[icv];
      ROM_DEF
      RoM[icv] = R;
      CP_GAM_DEF
      gamma[icv] = gam;

      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      T_E_DEF
      temp[icv] = T;

      // Compute mixture enthalpy, pressure, gamma and speed of sound
      H_DEF
      enthalpy[icv] = h;
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
   }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);

    
    // viscosity and thermal diffusivity computed in UgpWithCvCompFlow::calcMaterialProperties() at cell faces
  }
  
  // \brief Compute derivative of pressure with respect to scalar for coupled solution (Jacobi matrix)
  virtual void pressureDerivativeHookScalarRansComb() 
  {
    double dcp_dscal = Cp_1 - Cp_2;
    double dR_dscal = RoM_1 - RoM_2;
    
    for (int icv = 0; icv < ncv; icv++)
    {
      double E = enthalpy[icv] - RoM[icv] * temp[icv];
      double dgam_dscal = (gamma[icv] - 1.0) / RoM[icv] * ((1.0 - gamma[icv]) * dcp_dscal + gamma[icv] * dR_dscal);
      double dRT_dscal = (E + RoM[icv] * Treference) * dgam_dscal + gamma[icv] * Treference * dR_dscal;
      scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = rho[icv] * dRT_dscal;
    }    
  } 
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    double cp;
    double zm = Scal[ZMean_Index];
    ROM_DEF
    CP_GAM_DEF
    H_DEF
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    double cp;
    double zm = Scal[ZMean_Index];
    ROM_DEF
    CP_GAM_DEF
    T = p / (rho * R);
    H_DEF
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi_bfa;
    
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_T(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa],zM_fa[ifa],0.0,0.0);
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi_bfa;
    
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_H(T_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],h_bfa[ifa],zM_fa[ifa],0.0,0.0);
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
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_T(double &h, double &R, double &gam, double &mu, double &lamOcp, double T, double zm, double zv, double chi)
  {
    double cp;
    ROM_DEF
    CP_GAM_DEF
    H_DEF
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
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_H(double &T, double &R, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double chi)
  {
    double cp;
    ROM_DEF
    CP_GAM_DEF
    T_DEF
    if (mu_ref > 0.0)
    {
      mu = calcMuLam(T);
      lamOcp = mu / Pr;
    }
  }

  virtual void diffusivityHookScalarRansComb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "ZMean")
    {
      eq = getScalarTransportData(name);
 
      // internal faces
      for (int ifa = nfa_b; ifa < nfa; ifa++)
        eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            eq->diff[ifa] = lamOcp_fa[ifa];
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
           eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
      }
    }
  }

};


#endif  /* COMBMODEL_BINARYMIXING_H_ */
