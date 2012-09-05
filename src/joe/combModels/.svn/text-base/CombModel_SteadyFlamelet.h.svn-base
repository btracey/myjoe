/*
 * CombModel_SteadyFlamelet.h
 *
 *  Created on: May 6, 2009
 *      Author: Vincent Terrapon
 */

#ifndef COMBMODEL_STEADYFLAMELET_H_
#define COMBMODEL_STEADYFLAMELET_H_

#define MIXING_VISCOSITY

#include "CombModel_Base.h"
#include "ChemistryTableCartesianLinear.h"
#include "ChemistryTableCartesianCubic.h"
#include "ChemistryTableAdaptiveLinear.h"


// #############################################################################################
// ------                                                                                 ------
// ----                            RansCombSteadyFlamelet                                   ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Combustion model based on steady flamelet. 
 *
 *  The species mass fractions Species::Yk are functions of the mean scalar mixture fraction \p Zmean,
 *  the variance of the scalar mixture fraction \p Zvar and the scalar dissipation rate \p chi.
 *  The BACKWARD_EULER_COUPLED time integration scheme is not implemented!!!!!!!!!!!!!!!!
 */
template <class Chemtable>
class RansCombSteadyFlamelet: public RansCombBase {
public:
  // member variables

  Chemtable        myChemTable;                  ///< Chemistry table based on flamelets.
  Mixture          myMixture;                    ///< Mixture containing all the species information.
  double           *ZMean;                       ///< Mixture fraction Z.
  double           *ZVar;                        ///< Variance of mixture fraction.
//  double           *rhoDkZ;                      ///< Diffusion coefficient of mixture fraction Z.
  double           *muLam;                       ///< Laminar viscosity at cell center.
  double           *LambdaOverCp;                ///< Thermal diffusion coefficient at cell center.
  double           *chi;                         ///< Scalar dissipation rate (without ZVar)
  double           *chi_bfa;                     ///< Scalar dissipation rate at boundary faces (without ZVar)
  double           Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  double           Schmidt_turb_Z2;              ///< Turbulent Schmidt number for variance of mixture fraction.
  double           Cchi;                         ///< Constant parameter for scalar dissipation rate.
  double           Cmu;                          ///< Constant parameter for scalar dissipation rate with k-omega model
  int              ZMean_Index;                  ///< Index in scalar vector of ZMean.
  int              ZVar_Index;                   ///< Index in scalar vector of ZVar.
  int              omega_Index;                  ///< Index in scalar vector of omega.
  int              kine_Index;                   ///< Index in scalar vector of kine.
  int              eps_Index;                    ///< Index in scalar vector of eps.
  int              nuSA_Index;                   ///< Index in scalar vector of nuSA.
  vector<double*>  OutputSpecies;                ///< Vector of pointers to scalars/species to be saved (e.g., for output).
  vector<string>   OutputSpeciesName;            ///< Vector with name of scalars/species to be saved (e.g., for output).

public:
  // constructors

  RansCombSteadyFlamelet()
  {
    if (mpi_rank == 0)
      cout << "RansCombSteadyFlamelet()" << endl;

    // Overall constants
    Cmu = 0.09;
    Cchi            = getDoubleParam("Cchi", "2.0");
    Schmidt_turb_Z  = getDoubleParam("SCTURB_ZMean", "1.0");
    Schmidt_turb_Z2 = getDoubleParam("SCTURB_ZVar", "1.0");
    if (mpi_rank == 0) 
    {
      cout << "Turbulent Schmidt number for ZMean Sc_t=" << Schmidt_turb_Z << endl;
      cout << "Turbulent Schmidt number for ZVar Sc_t2=" << Schmidt_turb_Z2 << endl;
    }
    
    // Read name of additional scalars/species to save (e.g., for output) and register corresponding scalar
    Param * p;
    if (getParam(p, "OUTPUT_SPECIES"))
    {
      OutputSpeciesName.resize(p->getSize()-1);
      OutputSpecies.resize(p->getSize()-1);
      for (int i=0; i<p->getSize()-1; i++)
      {
        OutputSpecies[i] = NULL;
        OutputSpeciesName[i] = p->getString(i + 1);
        registerScalar(OutputSpecies[i], OutputSpeciesName[i], CV_DATA);
      }
    }

    // Scalar transport equations for ZMean and ZVar
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-4;
    eq->phiMaxiter = 200;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");

    eq = registerScalarTransport("ZVar", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZVar", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-4;
    eq->phiMaxiter = 200;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    
    // Other cell-centered quantities
//    rhoDkZ       = NULL;        registerScalar(rhoDkZ, "rhoDkZ", CV_DATA);
    muLam        = NULL;        registerScalar(muLam, "muLam" , CV_DATA);
    LambdaOverCp = NULL;        registerScalar(LambdaOverCp, "LambdaOverCp", CV_DATA);
    chi          = NULL;        registerScalar(chi, "chi", CV_DATA);
  }

  /*! \brief Destructor: unload chemistry table and clean memory. */
  virtual ~RansCombSteadyFlamelet()
  {
    myChemTable.Unload();
    myMixture.Unload();
    if (chi_bfa != NULL) delete []chi_bfa;
    
    if (mpi_rank == 0)
      cout << endl << "***** Steady Flamelet finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load chemistry table. */
  virtual void initialHookScalarRansCombModel()
  {
    int CombustionRegime = 1;
    myMixture.Load(CombustionRegime, &getThisParamMap());
    myChemTable.Load(getStringParam("CHEMTABLE_FILE"));

    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");
    eq = getScalarTransportData("ZVar");          ZVar  = eq->phi;           ZVar_Index  = getScalarTransportIndex("ZVar");

    // save indices for the turbulent variables (needed to compute chi)
    if ((turbModel == KOM) || (turbModel == KOMSST))
    {
      omega_Index = getScalarTransportIndex("omega");
      kine_Index  = getScalarTransportIndex("kine");
    }
    else if (turbModel == KEPS)
    {
      eps_Index   = getScalarTransportIndex("eps");
      kine_Index  = getScalarTransportIndex("kine");
    }
    else if (turbModel == SA)
      nuSA_Index  = getScalarTransportIndex("nuSA");
    
    chi_bfa = new double[nfa_b];

    if (mpi_rank == 0)
      cout << endl << "***** Steady Flamelet initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean, \p Zvar, \p chi.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, diffusion coefficient, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double E, cp;
    double kinecv = 0.0;

    ComputeScalarDissipation(chi);

    for (int icv = 0; icv < ncv; icv++)
    {
      // Read and interpolate species mass fractions Yk from chemistry table 
      myMixture.GetSpeciesMassFraction(ZMean[icv], ZVar[icv], chi[icv]*ZVar[icv], myChemTable);

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      temp[icv] = myMixture.ComputeMixTemperature(E, temp[icv]);

      // Compute mixture gas constant, enthalpy and pressure
      RoM[icv] = myMixture.GetMixR_o_M();
      enthalpy[icv] = myMixture.ComputeMixEnthalpy(temp[icv]);
      press[icv] = rho[icv] * RoM[icv] * temp[icv];

      if (isnan(press[icv]) || (press[icv]<= 0.0))
      {
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1]  << " / " << x_cv[icv][2] << endl;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(vel[icv], vel[icv]) << "  Z = " << ZMean[icv] << endl;
      }

      // Compute mixture viscosity, heat conductivity, gamma, ...
      cp = myMixture.ComputeMixCp(temp[icv]);
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
      
      myMixture.ComputeMixLambdaMu(temp[icv]);
      muLam[icv] = myMixture.GetMixMul();
      LambdaOverCp[icv] = myMixture.GetMixLambda() / cp;
//      rhoDkZ[icv] = LambdaOverCp[icv] / myMixture.GetMixLewis();
      
      // Save additional species/scalars for output
      for (int i=0; i<OutputSpeciesName.size(); i++)
      {
        string SpeciesName = OutputSpeciesName[i];
        double *Species = getScalar(SpeciesName.c_str());  
        Species[icv] = myMixture.GetMixYk(SpeciesName);
      }

    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
    updateCvData(muLam, REPLACE_DATA);
    updateCvData(LambdaOverCp, REPLACE_DATA);
//    updateCvData(rhoDkZ, REPLACE_DATA);
    
    for (int i=0; i<OutputSpeciesName.size(); i++)
    {
      string SpeciesName = OutputSpeciesName[i];
      double *Species = getScalar(SpeciesName.c_str());  
      updateCvData(Species, REPLACE_DATA);
    }
   
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
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double chi;
    
    // Compute scalar dissipation rate
    if ((turbModel == KOM) || (turbModel == KOMSST))
      chi = Cchi * Cmu * Scal[omega_Index];
    else if (turbModel == KEPS)
      chi = Cchi * Scal[eps_Index] / Scal[kine_Index];
    else if (turbModel == SA)
    {
      cout << "### RansCombSteadyFlamelet::caclThermoProp_T NOT YET IMPLEMENTED ###" << endl;
      //chi = Cchi * sqrt(Cmu) * vortMag[icv0];       // pass icv0
    }    
    else
      chi = 0.0;
    
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm, zv, chi, myChemTable);
    R = myMixture.GetMixR_o_M();
    h = myMixture.ComputeMixEnthalpy(T);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - R);
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double chi;
    
    // Compute scalar dissipation rate
    if ((turbModel == KOM) || (turbModel == KOMSST))
      chi = Cchi * Cmu * Scal[omega_Index];
    else if (turbModel == KEPS)
      chi = Cchi * Scal[eps_Index] / Scal[kine_Index];
    else if (turbModel == SA)
    {
      cout << "### RansCombSteadyFlamelet::caclThermoProp_T NOT YET IMPLEMENTED ###" << endl;
      //chi = Cchi * sqrt(Cmu) * vortMag[icv0];        // pass icv0
    }    
    else
      chi = 0.0;
    
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm, zv, chi, myChemTable);
    R = myMixture.GetMixR_o_M();
    T = p / (rho * R);
    h = myMixture.ComputeMixEnthalpy(T);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - R);
  }
  
  virtual double calcMuLam(int icv)
  {
    return muLam[icv];
  }
  
  virtual double calcMuLam(double temp)
  {
    cout << "### RansCombSteadyFlamelet::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
   
    ComputeScalarDissipation(chi_bfa, zone);

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_T(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa],zM_fa[ifa],zV_fa[ifa],chi_bfa[ifa]);
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {    
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
   
    ComputeScalarDissipation(chi_bfa, zone);

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_H(T_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],h_bfa[ifa],zM_fa[ifa],zV_fa[ifa],chi_bfa[ifa]);
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
  void ComputeProperties_T(double &h, double &RoM, double &gam, double &mu, double &lamOcp, double T, double zm, double zv, double chi)
  {
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm, zv, chi, myChemTable);
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
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_H(double &T, double &RoM, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double chi)
  {
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm, zv, chi, myChemTable);
    RoM = myMixture.GetMixR_o_M();
    T = myMixture.ComputeMixTemperature_H(h);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - RoM);
    myMixture.ComputeMixLambdaMu(T);
    mu = myMixture.GetMixMul();
    lamOcp = myMixture.GetMixLambda() / cp;
  }
  
  void ComputeScalarDissipation(double *chi)
  {
    if ((turbModel == KOM) || (turbModel == KOMSST))
    {
      double *om = getR1("omega");
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * Cmu * om[icv];
    }
    
    else if (turbModel == KEPS)
    {
      double *epsi = getR1("eps");
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * epsi[icv] / kine[icv];
    }
    
    else if (turbModel == SA)
    {
      calcVorticity();
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * sqrt(Cmu) * vortMag[icv];
    }
    
    else
    {
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = 0.0;
    }
    updateCvData(chi, REPLACE_DATA);
  }
  
  void ComputeScalarDissipation(double *chi_bfa, FaZone *zone)
  {
    Param *param;
    if (getParam(param, zone->getName()))
      if (param->getString() == "WALL")
      {
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          chi_bfa[ifa] = chi[icv0];
        }
      }

      else
      {            
        if ((turbModel == KOM) || (turbModel == KOMSST))
        {
          ScalarTranspEq *eq;
          eq = getScalarTransportData("omega");
          double *omega_fa = eq->phi_bfa;
         
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            chi_bfa[ifa] = Cchi * Cmu * omega_fa[ifa];
        }
        
        else if (turbModel == KEPS)
        {
          ScalarTranspEq *eq;
          eq = getScalarTransportData("kine");
          double *kine_fa = eq->phi_bfa;
          eq = getScalarTransportData("eps");
          double *eps_fa = eq->phi_bfa;
        
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            chi_bfa[ifa] = Cchi * eps_fa[ifa] / kine_fa[ifa];
        }
        
        else if (turbModel == SA)
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            chi_bfa[ifa] = Cchi * sqrt(Cmu) * vortMag[icv0];
          }
        }
        
        else
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            chi_bfa[ifa] = 0.0;
        }
     }
  }

  virtual void diffusivityHookScalarRansComb(const string &name)
  {
    ScalarTranspEq *eq;

    if ((name == "ZMean") || (name == "ZVar") || (name == "CMean"))
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
  
  virtual void sourceHookScalarRansComb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "ZVar")
    {
      ScalarTranspEq *eq = getScalarTransportData("ZMean");       
      double *ZMean_bfa = eq->phi_bfa;  
      double *ZMean_grad[3] = eq->grad_phi;

      for (int icv = 0; icv < ncv; icv++)
      {
        double src = 2.0*InterpolateAtCellCenterFromFaceValues(mut_fa, icv)/Schmidt_turb_Z2 * vecDotVec3d(ZMean_grad[icv], ZMean_grad[icv]) - chi[icv] * rho[icv] * ZVar[icv];
        rhs[icv] += src * cv_volume[icv];
      }

      if (flagImplicit)
      {
        for (int icv = 0; icv < ncv; icv++)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - chi[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

};


#endif /* COMBMODEL_STEADYFLAMELET_H_ */
