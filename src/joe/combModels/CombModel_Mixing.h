/*
 * CombModel_Mixing.h
 *
 *  Created on: May 6, 2009
 *      Author: Vincent Terrapon
 */

#ifndef COMBMODEL_MIXING_H_
#define COMBMODEL_MIXING_H_

#define MIXING_VISCOSITY

#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                               RansCombMixing                                        ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Mixing model based on mixture fraction. 
 *
 *  The species mass fractions Species::Yk are functions of the mean scalar mixture fraction \p Zmean
 *  only through linear mixing rules. There is no combustion, species only mix. The amount of fuel in
 *  the fuel stream can be increased over time to help convergence.
 */
class RansCombMixing: public RansCombBase {
public:
  // member variables

  Mixture myMixture;                    ///< Mixture containing all the species information.
  double  *ZMean;                       ///< Mixture fraction Z.
//  double  *rhoDkZ;                      ///< Diffusion coefficient of mixture fraction Z.
  double  *muLam;                       ///< Laminar viscosity at cell center.
  double  *LambdaOverCp;                ///< Thermal diffusion coefficient at cell center.
  double  Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  int     ZMean_Index;                  ///< Index in scalar vector of ZMean.

public:
  // constructors

  RansCombMixing()
  {
    if (mpi_rank == 0)
      cout << "RansCombMixing()" << endl;

    // Overall constants
    Schmidt_turb_Z = getDoubleParam("SCTURB_ZMean", "0.5");
    if (mpi_rank == 0) 
      cout << "Turbulent Schmidt number for Z Sc_t=" << Schmidt_turb_Z << endl;
  
    // Scalar transport equation for ZMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -1.0;
    eq->upperBound = 2.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";
    
    // Other cell-centered quantities
//    rhoDkZ       = NULL;        registerScalar(rhoDkZ, "rhoDkZ", CV_DATA);
    muLam        = NULL;        registerScalar(muLam, "muLam" , CV_DATA);
    LambdaOverCp = NULL;        registerScalar(LambdaOverCp, "LambdaOverCp", CV_DATA);
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombMixing()
  {
    if (scalarTranspEqVector[ZMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZMean_Index].dpress_dphi; 
      scalarTranspEqVector[ZMean_Index].dpress_dphi = NULL;
    }

    myMixture.Unload();
    if (mpi_rank == 0)
      cout << endl << "***** Mixing finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    int CombustionRegime = 0;
    myMixture.Load(CombustionRegime, &getThisParamMap());

    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");          ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");
    
    // For coupled solution
    if ((getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_COUPLED") || (getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_SEMICOUPLED")
				|| (getStringParam("TIME_INTEGRATION") == "BDF2_SEMICOUPLED"))
      scalarTranspEqVector[ZMean_Index].dpress_dphi = new double[ncv_g];

    if (mpi_rank == 0)
      cout << endl << "***** Mixing initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double E, cp;
    double kinecv = 0.0;

    for (int icv = 0; icv < ncv; icv++)
    {
      // Read and interpolate species mass fractions Yk from chemistry table 
      if ( isnan(ZMean[icv]) ) {
        cout << "WARNING -- ZMean : " ;
        cout << "x = " << x_cv[icv][0] << " , " << x_cv[icv][1] << " , " << x_cv[icv][2] << "  " ;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / (rho[icv] * rho[icv]) << "  Z = " << ZMean[icv] << endl;
        fstream fin; fin.open("killjoe",fstream::in); if (!fin.is_open()){ fstream fout; fout.open("killjoe",fstream::out); fout.close();} else fin.close();
        cout.flush() ;
      }
      myMixture.GetSpeciesMassFraction(ZMean[icv]);

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      if (isnan(temp[icv]) || (temp[icv] < 0.001))
        temp[icv] = myMixture.ComputeMixTemperature(E);
      else
        temp[icv] = myMixture.ComputeMixTemperature(E, temp[icv]);
      if (isnan(temp[icv]) || (temp[icv]<= 0.0)) {
        cout << "WARNING -- temp : " ;
        cout << "x = " << x_cv[icv][0] << " , " << x_cv[icv][1] << " , " << x_cv[icv][2] << "  " ;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / (rho[icv] * rho[icv]) << "  Z = " << ZMean[icv] << endl;
        fstream fin; fin.open("killjoe",fstream::in); if (!fin.is_open()){ fstream fout; fout.open("killjoe",fstream::out); fout.close();} else fin.close();
        cout.flush() ;
      }

      // Compute mixture gas constant, enthalpy and pressure
      RoM[icv] = myMixture.GetMixR_o_M();
      enthalpy[icv] = myMixture.ComputeMixEnthalpy(temp[icv]);
      press[icv] = rho[icv] * RoM[icv] * temp[icv];

      if (isnan(press[icv]) || (press[icv]<= 0.0))
      {
        cout << "WARNING -- press : " ;
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1] << " / " << x_cv[icv][2] << "  " ;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / (rho[icv] * rho[icv]) << "  Z = " << ZMean[icv] << "  kine = " << kine[icv] << endl ;
        fstream fin; fin.open("killjoe",fstream::in); if (!fin.is_open()){ fstream fout; fout.open("killjoe",fstream::out); fout.close();} else fin.close();
      }

      // Compute mixture viscosity, heat conductivity, gamma
      cp = myMixture.ComputeMixCp(temp[icv]);
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);

      myMixture.ComputeMixLambdaMu(temp[icv]);
      muLam[icv] = myMixture.GetMixMul();
      LambdaOverCp[icv] = myMixture.GetMixLambda() / cp;
//      rhoDkZ[icv] = LambdaOverCp[icv] / myMixture.GetMixLewis();
   }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);

    updateCvData(LambdaOverCp, REPLACE_DATA);
    updateCvData(muLam, REPLACE_DATA);
//    updateCvData(rhoDkZ, REPLACE_DATA);
    
//    // Just for debugging: create table of energy and enthalpy as function of T and Z
//    if (mpi_rank == 0)
//      {
//        double T, zm;
//        Mixing_Factor = 1.0;
//        cout << endl << "Temperature table: Enthalpy" << endl;
//        for (int i = 0; i < 100; i++)
//          {
//            T = (double)i * 10.0 + 10.0;
//            cout << T << "  ";
//            for (int j = 0; j < 11; j++)
//              {
//                zm = (double) j * 0.1;
//                myMixture.GetSpeciesMassFraction(zm, Mixing_Factor);
//                cout << myMixture.ComputeMixEnthalpy(T) << "  ";
//              }
//            cout << endl;
//          }
//        cout << endl << "Temperature table: Energy" << endl;
//        for (int i = 0; i < 300; i++)
//          {
//            T = (double)i * 10.0 + 10.0;
//            cout << T << "  ";
//            for (int j = 0; j < 11; j++)
//              {
//                zm = (double) j * 0.1;
//                myMixture.GetSpeciesMassFraction(zm, Mixing_Factor);
//                cout << myMixture.ComputeMixEnergy(T) << "  ";
//              }
//            cout << endl;
//          }
//      }
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
 
        mul_fa[ifa] = (w1 * muLam[icv0] + w0 * muLam[icv1]) / (w0 + w1);
        lamOcp_fa[ifa] = (w1 * LambdaOverCp[icv0] + w0 * LambdaOverCp[icv1]) / (w0 + w1);
        
//        // New interpolation for the viscosity: interpolate Z and T, not mu
//        double Zmm = (w1 * ZMean[icv0] + w0 * ZMean[icv1]) / (w0 + w1);
//        double TT  = (w1 * temp[icv0] + w0 * temp[icv1]) / (w0 + w1);
//        myMixture.GetSpeciesMassFraction(Zmm);
//        double cp = myMixture.ComputeMixCp(TT);
//        myMixture.ComputeMixLambdaMu(TT);
//        mul_fa[ifa] = myMixture.GetMixMul() * 10.0;
//        lamOcp_fa[ifa] = myMixture.GetMixLambda() / cp;        
      }    
      
      // boundary faces computed in setBC
    }
  }
#endif  
  
  // \brief Compute derivative of pressure with respect to scalar for coupled solution (Jacobi matrix)
  virtual void pressureDerivativeHookScalarRansComb() 
  {
    double delta_ZM = 1.0e-5;
    double pressp, pressm;
    
    for (int icv = 0; icv < ncv; icv++)
    {
      double E = enthalpy[icv] - RoM[icv] * temp[icv];
      
      // Derivative with respect to ZMean
      double ZMp = ZMean[icv] + delta_ZM;
      double ZMm = ZMean[icv] - delta_ZM;
      pressure_scalSource(pressp, rho[icv], E, temp[icv], ZMp);
      pressure_scalSource(pressm, rho[icv], E, temp[icv], ZMm);
      scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = (pressp - pressm) / (2.0 * delta_ZM);
    }    
    updateCvData(scalarTranspEqVector[ZMean_Index].dpress_dphi, REPLACE_DATA);
  } 
  
  void pressure_scalSource(double &pp, double &rrho, double &E, double &ttemp, double &ZM)
  {
    myMixture.GetSpeciesMassFraction(ZM);
    pp = rrho * myMixture.GetMixR_o_M() * myMixture.ComputeMixTemperature(E, ttemp);
  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    double zm = Scal[ZMean_Index];
    
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm);
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
    
    // Read and interpolate species mass fractions Yk from chemistry table 
    myMixture.GetSpeciesMassFraction(zm);
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
    cout << "### RansCombMixing::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");        double *zM_fa = eq->phi_bfa;
    
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_T(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa],zM_fa[ifa],0.0,0.0);
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");        double *zM_fa = eq->phi_bfa;
    
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
  void ComputeProperties_T(double &h, double &RoM, double &gam, double &mu, double &lamOcp, double T, double zm, double zv, double chi)
  {
    // Read and interpolate species mass fractions Yk from linear mixing rules
    myMixture.GetSpeciesMassFraction(zm);
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
    // Read and interpolate species mass fractions Yk from linear mixing rules
    myMixture.GetSpeciesMassFraction(zm);
    RoM = myMixture.GetMixR_o_M();
    T = myMixture.ComputeMixTemperature_H(h);
    double cp = myMixture.ComputeMixCp(T);
    gam = cp / (cp - RoM);
    myMixture.ComputeMixLambdaMu(T);
    mu = myMixture.GetMixMul();
    lamOcp = myMixture.GetMixLambda() / cp;
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


#endif  /* COMBMODEL_MIXING_H_ */
