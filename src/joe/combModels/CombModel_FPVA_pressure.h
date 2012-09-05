/*
 * CombModel_FPVA_pressure.h
 *
 *  Created on: Mar 11, 2011
 *      Author: Ronan Vicquelin
 */

#ifndef COMBMODEL_FPVA_PRESSURE_H_
#define COMBMODEL_FPVA_PRESSURE_H_

#define PRESSURE_SCALING
#define MIXING_VISCOSITY

#include "CombModel_Base.h"
#include "PressureTable.h"

#define Tminimum 20.0


// #############################################################################################
// ------                                                                                 ------
// ----                            RansCombFPVA_new                                         ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Combustion model based on Flamelet Progress Variable Approach (FPVA).
 *
 *  The coefficients R0, T0, E0, ... are functions of the mean scalar mixture fraction \p Zmean,
 *  the variance of the scalar mixture fraction \p Zvar and the progress variable \p Cmean.
 */
class RansCombFPVA_pressure: public RansCombBase {
public:
  // member variables

  PressureTable  myChemTable;                  ///< Chemistry table based on flamelets.
  double             *ZMean;                       ///< Mixture fraction Z (mean).
  double             *ZVar;                        ///< Variance of mixture fraction.
  double             *CMean;                       ///< Progress variable C (mean).
  double             *muLam;                       ///< Laminar viscosity at cell center.
  double             *LambdaOverCp;                ///< Thermal diffusion coefficient at cell center.
  double             *chi;                         ///< Scalar dissipation rate (without ZVar).
  double             *CmeanSource;                 ///< Source term for progress variable C.
  double             *Zv_norm;                     ///< Normalized mixture fraction variance
  double             *C_norm;                      ///< Normalized progress variable 
  double             *Pressure_out;                       ///< Pressure find by the database
  double             Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  double             Schmidt_turb_Z2;              ///< Turbulent Schmidt number for variance of mixture fraction.
  double             Schmidt_turb_C;               ///< Turbulent Schmidt number for progress variable C.
  double             Cchi;                         ///< Constant parameter for scalar dissipation rate.
  double             Cmu;                          ///< Constant parameter for scalar dissipation rate with k-omega model.
  int                ZMean_Index;                  ///< Index in scalar vector of ZMean.
  int                ZVar_Index;                   ///< Index in scalar vector of ZVar.
  int                CMean_Index;                  ///< Index in scalar vector of CMean.
  vector<double*>    OutputSpecies;                ///< Vector of pointers to scalars/species to be saved (e.g., for output).
  vector<string>     OutputSpeciesName;            ///< Vector with name of scalars/species to be saved (e.g., for output).
  string             Comb_Reg;                     ///< Combustion regime (cold/hot/ignition) to set the progress variable.
  double             *dCmeanSource_dZM;            ///< Derivative of progress variable source term with respect to mean mixture fraction.
  double             *dCmeanSource_dZV;            ///< Derivative of progress variable source term with respect to mixture fraction variance.
  double             *dCmeanSource_dCM;            ///< Derivative of progress variable source term with respect to mean progress variable.

private:
  double T0;
  double E0;
  double GAMMA0;
  double AGAMMA;
  double MU0;
  double AMU;
  double LOC0;
  double ALOC;

public:
  // constructors

  RansCombFPVA_pressure()
  {
    if (mpi_rank == 0)
      cout << "RansCombFPVA_Coeff()" << endl;

    Comb_Reg = getStringParam("COMBUSTION_REGIME", "HOT");

    // Overall constants
    Cmu = 0.09;
    Cchi            = getDoubleParam("Cchi", "2.0");
    Schmidt_turb_Z  = getDoubleParam("SCTURB_ZMean", "0.5");
    Schmidt_turb_Z2 = getDoubleParam("SCTURB_ZVar", "0.5");
    Schmidt_turb_C  = getDoubleParam("SCTURB_CMean", "0.5");
    if (mpi_rank == 0)
    {
      cout << "Turbulent Schmidt number for ZMean Sc_t = " << Schmidt_turb_Z  << endl;
      cout << "Turbulent Schmidt number for ZVar Sc_t2 = " << Schmidt_turb_Z2 << endl;
      cout << "Turbulent Schmidt number for CMean Sc_t = " << Schmidt_turb_C  << endl;
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

    // Scalar transport equations for ZMean, ZVar and CMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";

    eq = registerScalarTransport("ZVar", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZVar", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";

    eq = registerScalarTransport("CMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_CMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";

    // Other cell-centered quantities
    muLam        = NULL;        registerScalar(muLam, "muLam" , CV_DATA);
    LambdaOverCp = NULL;        registerScalar(LambdaOverCp, "LambdaOverCp", CV_DATA);
    chi          = NULL;        registerScalar(chi, "chi", CV_DATA);
    CmeanSource  = NULL;        registerScalar(CmeanSource, "CmeanSource", CV_DATA);
    Zv_norm      = NULL;        registerScalar(Zv_norm, "Zv_norm", CV_DATA);
    C_norm       = NULL;        registerScalar(C_norm, "C_norm", CV_DATA);
    Pressure_out  = NULL;        registerScalar(Pressure_out, "P_out", CV_DATA);

    // Test if any BCs are not implemented yet
    bool testBC = check_CBC_INLET_SUBSONIC();
    if (testBC) throw(-1);
  }

  /*! \brief Destructor: unload chemistry table and clean memory. */
  virtual ~RansCombFPVA_pressure()
  {
    myChemTable.Unload();

    if (scalarTranspEqVector[ZMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZMean_Index].dpress_dphi;
      scalarTranspEqVector[ZMean_Index].dpress_dphi = NULL;
    }
    if (scalarTranspEqVector[ZVar_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZVar_Index].dpress_dphi;
      scalarTranspEqVector[ZVar_Index].dpress_dphi = NULL;
    }
    if (scalarTranspEqVector[CMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[CMean_Index].dpress_dphi;
      scalarTranspEqVector[CMean_Index].dpress_dphi = NULL;
    }

    if (dCmeanSource_dZM != NULL)
    {
      delete [] dCmeanSource_dZM;
      dCmeanSource_dZM = NULL;
    }
    if (dCmeanSource_dZV != NULL)
    {
      delete [] dCmeanSource_dZV;
      dCmeanSource_dZV = NULL;
    }
    if (dCmeanSource_dCM != NULL)
    {
      delete [] dCmeanSource_dCM;
      dCmeanSource_dCM = NULL;
    }

    if (mpi_rank == 0)
      cout << endl << "***** FPVA_Coeff Combustion finalized *****" << endl << endl;
  }

public:
  // member functions

  /*! \brief Read thermo input file and load chemistry table. */
  virtual void initialHookScalarRansCombModel()
  {
    //myChemTable.Load(getStringParam("CHEMTABLE_FILE"));
    myChemTable.MPI_Load(getStringParam("CHEMTABLE_FILE"));
    if (mpi_rank == 0) myChemTable.print();

    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");          ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");
    eq = getScalarTransportData("ZVar");           ZVar  = eq->phi;           ZVar_Index  = getScalarTransportIndex("ZVar");
    eq = getScalarTransportData("CMean");          CMean = eq->phi;           CMean_Index = getScalarTransportIndex("CMean");

    // For coupled solution
    if ((getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_COUPLED") || (getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_SEMICOUPLED"))
    {
      scalarTranspEqVector[ZMean_Index].dpress_dphi = new double[ncv_g];
      scalarTranspEqVector[ZVar_Index ].dpress_dphi = new double[ncv_g];
      scalarTranspEqVector[CMean_Index].dpress_dphi = new double[ncv_g];
      dCmeanSource_dZM = new double[ncv_g];
      dCmeanSource_dZV = new double[ncv_g];
      dCmeanSource_dCM = new double[ncv_g];
    }

    if (mpi_rank == 0)
      cout << endl << "***** FPVA_Coeff Combustion initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean, \p Zvar, \p chi.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, diffusion coefficient, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    //cout << " < calcStateVariables()" << endl;
    double E, cp;
    double kinecv = 0.0;

    ComputeScalarDissipation(chi);
    double P_out, Sz_out, C_out;
    for (int icv = 0; icv < ncv; icv++)
    {
      // Read and interpolate coefficients from chemistry table and save source term for progress variable

      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMean[icv], ZVar[icv], CMean[icv],rhoE[icv],
                                        P_out, Sz_out, C_out);
      Zv_norm[icv] = Sz_out;
      C_norm[icv]  = C_out;
      Pressure_out[icv]   = P_out;

      myChemTable.LookupCoeff(RoM[icv], T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CmeanSource[icv],
                              P_out, ZMean[icv], Sz_out, C_out);

      if (Comb_Reg == "COLD")
        CmeanSource[icv] = 0.0;

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];

      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      if (AGAMMA == 0.0)
        temp[icv] = T0 + (GAMMA0 - 1.0) / RoM[icv] * (E - E0);
      else
        temp[icv] = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E - E0) / RoM[icv]) - 1.0);

      // Temperature clipping
      temp[icv] = max(temp[icv], Tminimum);

      // Compute mixture enthalpy and pressure
      enthalpy[icv] = E + RoM[icv] * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
      if (press[icv] != P_out) {
        cout << "Error Pressure " << P_out << " " << press[icv] << endl;
        throw(-1);
      }

      // Compute mixture viscosity, heat conductivity, gamma, speed of sound
      muLam[icv] = MU0 * pow(temp[icv] / T0, 0.7);
      LambdaOverCp[icv] = LOC0 * pow(temp[icv] / T0, 0.62);
      gamma[icv] = GAMMA0 + AGAMMA * (temp[icv] - T0);
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);

      if (isnan(press[icv]) || (press[icv]<= 0.0))
      {
        cout << "WARNING! :" << endl;
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1]  << " / " << x_cv[icv][2] << endl;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(vel[icv], vel[icv]) << "  Z = " << ZMean[icv] << endl;
        throw(-1);
      }

      // Save additional species/scalars for output
      for (int i=0; i<OutputSpeciesName.size(); i++)
      {
        string SpeciesName = OutputSpeciesName[i];
        double *Species = getScalar(SpeciesName.c_str());
        Species[icv] = myChemTable.Lookup(SpeciesName);
      }
      //cout << icv << " : done!" << endl;

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
    updateCvData(CmeanSource, REPLACE_DATA);
//    updateCvData(rhoDkZ, REPLACE_DATA);

    for (int i=0; i<OutputSpeciesName.size(); i++)
    {
      string SpeciesName = OutputSpeciesName[i];
      double *Species = getScalar(SpeciesName.c_str());
      updateCvData(Species, REPLACE_DATA);
    }
    //cout << " < calcStateVariables : End" << endl;
  }

#ifdef MIXING_VISCOSITY
  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *
   */
  virtual void calcMaterialProperties()
  {
    //cout << " < calcMaterialProperties()" << endl;
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
    //cout << " < calcMaterialProperties : End" << endl;
  }
#endif

  // \brief Compute derivative of pressure with respect to scalar for coupled solution (Jacobi matrix)
  virtual void pressureDerivativeHookScalarRansComb()
  {
    //cout << " < pressureDerivativeHookScalarRansComb()" << endl;
    double delta_ZM = 1.0e-5;
    double delta_ZV = 1.0e-8;
    double delta_CM = 1.0e-5;

    double pressp, pressm, CmeanSourcep, CmeanSourcem;
    double P_out, Sz_out, C_out;
    for (int icv = 0; icv < ncv; icv++)
    {
      double E = enthalpy[icv] - RoM[icv] * temp[icv];

      // Derivative with respect to ZMean
      double ZMp = ZMean[icv] + delta_ZM;
      double ZMm = ZMean[icv] - delta_ZM;
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMp, ZVar[icv], CMean[icv], rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcep = myChemTable.Lookup(P_out, ZMp, Sz_out, C_out, "SRC_PROG");
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMm, ZVar[icv], CMean[icv], rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcem = myChemTable.Lookup(P_out, ZMm, Sz_out, C_out, "SRC_PROG");

      scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = 0.0;
      dCmeanSource_dZM[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZM);

      // Derivative with respect to ZVar
      double ZVp = ZVar[icv] + delta_ZV;
      double ZVm = ZVar[icv] - delta_ZV;
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMean[icv], ZVp, CMean[icv], rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcep = myChemTable.Lookup(P_out, ZMean[icv], Sz_out, C_out, "SRC_PROG");
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMean[icv], ZVm, CMean[icv], rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcem = myChemTable.Lookup(P_out, ZMean[icv], Sz_out, C_out, "SRC_PROG");

      scalarTranspEqVector[ZVar_Index].dpress_dphi[icv] = 0.0;
      dCmeanSource_dZV[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZV);

      // Derivative with respect to CMean
      double CMp = CMean[icv] + delta_CM;
      double CMm = CMean[icv] - delta_CM;
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMean[icv], ZVar[icv], CMp, rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcep = myChemTable.Lookup(P_out, ZMean[icv], Sz_out, C_out, "SRC_PROG");
      myChemTable.ComputeDatabaseInputs_new(rho[icv], ZMean[icv], ZVar[icv], CMm, rhoE[icv], P_out, Sz_out, C_out);
      CmeanSourcem = myChemTable.Lookup(P_out, ZMean[icv], Sz_out, C_out, "SRC_PROG");

      scalarTranspEqVector[CMean_Index].dpress_dphi[icv] = 0.0;
      dCmeanSource_dCM[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_CM);
    }

    updateCvData(scalarTranspEqVector[ZMean_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(scalarTranspEqVector[ZVar_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(scalarTranspEqVector[CMean_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(dCmeanSource_dZM, REPLACE_DATA);
    updateCvData(dCmeanSource_dZV, REPLACE_DATA);
    updateCvData(dCmeanSource_dCM, REPLACE_DATA);
    //cout << " < pressureDerivativeHookScalarRansComb : End" << endl;
  }


  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    //cout << " < calcThermoProp_T()" << endl;
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double cm = Scal[CMean_Index];
    double P_out, Sz_out, C_out;
    myChemTable.ComputeDatabaseInputs_FromT(rho, zm, zv, cm, T, P_out, Sz_out, C_out);
    myChemTable.LookupSelectedCoeff(R, T0, E0, GAMMA0, AGAMMA, P_out, zm, Sz_out, C_out);
    if (AGAMMA == 0.0)
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    p = rho * R * T;
    //cout << " < calcThermoProp_T : End" << endl;
  }

  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    //cout << " < calcThermoProp_p()" << endl;
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double cm = Scal[CMean_Index];
    // Sz
    double Sz = myChemTable.ComputeNormalizedVariance(zm, zv);
    // c
    double C_out = myChemTable.ComputeNormalizedProg(p, zm, Sz, cm);
    //
    myChemTable.LookupSelectedCoeff(R, T0, E0, GAMMA0, AGAMMA, p, zm, Sz, C_out);
    T = p / (rho * R);
    if (AGAMMA == 0.0)
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    //cout << " < calcThermoProp_p : End" << endl;
  }

  virtual double calcMuLam(int icv)
  {
    return muLam[icv];
  }

  virtual double calcMuLam(double temp)
  {
    //cout << "### RansCombFPVA::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
  }

  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    //cout << " < ComputeBCProperties_T()" << endl;
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
    eq = getScalarTransportData("CMean");       double *cM_fa = eq->phi_bfa;

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_TP(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa],p_bfa[ifa],zM_fa[ifa],zV_fa[ifa],cM_fa[ifa]);
    //cout << " < ComputeBCProperties_T : End" << endl;
  }

  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
    eq = getScalarTransportData("CMean");       double *cM_fa = eq->phi_bfa;

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_H_old(T_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],h_bfa[ifa],zM_fa[ifa],zV_fa[ifa],cM_fa[ifa]);
  }

  /*! \brief Compute for a given temperature & pressure the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition.
   *  \param[out] H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  T Temperature in [K].
   *  \param[in]  P Pressure in [Pa].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv normalized  Variance of mixture fraction.
   *  \param[in]  cm normalized progress variable.
   */
  void ComputeProperties_TP(double &h, double &RoM, double &gam, double &mu, double &lamOcp, 
                            const double T, const double P, const double zm, const double zv, const double cm)
  {
    //cout << " < ComputeProperties_TP()" << endl;
    // Sz
    double Sz = myChemTable.ComputeNormalizedVariance(zm, zv);
    // c
    double C_out = myChemTable.ComputeNormalizedProg(P, zm, Sz, cm);
    //
    double dummy;
    myChemTable.LookupCoeff(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, P, zm, Sz, C_out);
    if (AGAMMA == 0.0)
      h = E0 + RoM / (GAMMA0 - 1.0) * (T - T0) + RoM * T;
    else
      h = E0 + RoM / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + RoM * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    mu = MU0 * pow(T / T0, 0.7);
    lamOcp = LOC0 * pow(T / T0, 0.62);
    //cout << " < ComputeProperties_TP : End" << endl;
  }

  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   */
  void ComputeProperties_H(double &T, double &RoM, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double cm)
  {
    //cout << "ComputeProperties_H not implemented yet" << endl;
    //throw(-1);
    // zv and cm are assumed null
    // Without Newton ? Simple relation between h-h0 and T-T0 ?
  }

  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   */
  void ComputeProperties_H_old(double &T, double &RoM, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double cm)
  {
    //cout << "ComputeProperties_H not implemented yet" << endl;
    //throw(-1);
    // zv and cm are assumed null
    double Pref = 1.5e5;
    double dummy;
    myChemTable.LookupCoeff(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, Pref, zm, 0.0, 0.0);
    T = (h - E0) / RoM;
    T = SolveNewtonTemperature_H(RoM, T0, E0, GAMMA0, AGAMMA, h, T);
    T = max(T, Tminimum);
    gam = GAMMA0 + AGAMMA * (T - T0);
    mu = MU0 * pow(T / T0, 0.7);
    lamOcp = LOC0 * pow(T / T0, 0.62);
  }

  double SolveNewtonTemperature_H(double R, double t0, double e0, double gamma0, double agamma, double h, double Tguess)
  {
    int    itermax = 50, iter;
    double errormax = 1.0e-7, error;
    double Tresidual, hh, cpp;

    iter = 0;
    do
    {
      ++iter;
      if (agamma == 0.0)
        hh = e0 + R / (gamma0 - 1.0) * (Tguess - t0) + R * Tguess;
      else
        hh = e0 + R / agamma * log(1.0 + agamma * (Tguess - t0) / (gamma0 - 1.0)) + R * Tguess;
      cpp = R * (1.0 + 1.0 / (gamma0 + agamma * (Tguess - t0) - 1.0));

      Tresidual = Tguess - (hh - h) / cpp;

      error = Tguess - Tresidual;
      if (fabs(error) <= errormax)
        return Tresidual;

      Tguess = Tresidual;

    }
    while (iter < itermax);

    /* if the iteration has not converged, exit with error */
    cerr << "### Computation of temperature in Newton iteration has not converged: Tguess=" << Tguess + error
         << ", Tresidual=" << Tresidual << ", enthalpy=" << h << ", #iterations=" << iter << " ###" << endl;
    throw(-1);
  }

  void ComputeScalarDissipation(double *chi)
  {
    //if (mpi_rank == 0) cout << " < ComputeScalarDissipation()" << endl;
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
    //if (mpi_rank == 0) cout << " < ComputeScalarDissipation : End" << endl;
  }

  virtual void diffusivityHookScalarRansComb(const string &name)
  {
    //if (mpi_rank == 0) cout << " < diffusivityHookScalarRansComb()" << endl;
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
    //if (mpi_rank == 0) cout << " < diffusivityHookScalarRansComb : End" << endl;
  }

  virtual void sourceHookScalarRansComb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    //if (mpi_rank == 0) cout << " < sourceHookScalarRansComb_new()" << endl;
    if (name == "ZVar")
    {
      ScalarTranspEq *eq = getScalarTransportData("ZMean");
      double *ZMean_bfa = eq->phi_bfa;
      double (*ZMean_grad)[3] = eq->grad_phi;

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

    if (name == "CMean")
    {
      for (int icv = 0; icv < ncv; icv++)
      {
        rhs[icv] += rho[icv] * CmeanSource[icv] * cv_volume[icv];
      }

      if (flagImplicit)
      {
        for (int icv = 0; icv < ncv; icv++)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = 0.0;
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
    //if (mpi_rank == 0) cout << " < sourceHookScalarRansComb_new : End" << endl;
  }

  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    //if (mpi_rank == 0) cout << " < sourceHookRansCombCoupled()" << endl;
    ScalarTranspEq *eq = getScalarTransportData("ZMean");
    double (*ZMean_grad)[3] = eq->grad_phi;


    int ZMean_Coupled_Index = getScalarTransportCoupledIndex("ZMean");
    if (ZMean_Coupled_Index == -1)
    {
      cerr << "### Error ZMean_Coupled_Index is " << ZMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int ZVar_Coupled_Index = getScalarTransportCoupledIndex("ZVar");
    if (ZVar_Coupled_Index == -1)
    {
      cerr << "### Error ZVar_Coupled_Index is " << ZVar_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int CMean_Coupled_Index = getScalarTransportCoupledIndex("CMean");
    if (CMean_Coupled_Index == -1)
    {
      cerr << "### Error CMean_Coupled_Index is " << CMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      rhs[icv][5+ZVar_Coupled_Index ] += (2.0 * InterpolateAtCellCenterFromFaceValues(mut_fa, icv) / Schmidt_turb_Z2 * vecDotVec3d(ZMean_grad[icv], ZMean_grad[icv]) - chi[icv] * rho[icv] * ZVar[icv]) * cv_volume[icv];
      rhs[icv][5+CMean_Coupled_Index] += rho[icv] * CmeanSource[icv] * cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
        A[noc00][5+ZVar_Coupled_Index ][5+ZVar_Coupled_Index ] -= - chi[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+ZMean_Coupled_Index] -= dCmeanSource_dZM[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+ZVar_Coupled_Index ] -= dCmeanSource_dZV[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+CMean_Coupled_Index] -= dCmeanSource_dCM[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][0] -= (CmeanSource[icv] - ZMean[icv] * dCmeanSource_dZM[icv] - ZVar[icv] * dCmeanSource_dZV[icv] - CMean[icv] * dCmeanSource_dCM[icv]) * cv_volume[icv];
     }
    }
    //if (mpi_rank == 0) cout << " < sourceHookRansCombCoupled : End" << endl;
  }

  bool check_CBC_INLET_SUBSONIC(void) {
    bool test = false;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY) {
        Param *param;
        if (getParam(param, zone->getName())) {
          if (param->getString() == "CBC_SUBSONIC_INLET") {
             string FaName = zone->getName();
             Param *p;
             double phiBC;
             getParam(p,FaName + "." + "ZVar");
             phiBC = p->getDouble(1);
             if (phiBC != 0.0) {
               cerr << "CBC_SUBSONIC_INLET not implmented for ZVar <> 0.0" << endl;
               test = true;
             }
             getParam(p,FaName + "." + "CMean");
             phiBC = p->getDouble(1);
             if (phiBC != 0.0) {
               cerr << "CBC_SUBSONIC_INLET not implmented for CMean <> 0.0" << endl;
               throw(-1);
               test = true;
             }
          }
        }
      }
    return test;
  }
};

#endif /* COMBMODEL_FPVA_PRESSURE_H_ */
