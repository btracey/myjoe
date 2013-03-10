#ifndef RANSTURBMODEL_V2F_half_H
#define RANSTURBMODEL_V2F_half_H

#include "UgpWithCvCompFlow.h"


class RansTurbV2F_half : virtual public UgpWithCvCompFlow
{
public:
  RansTurbV2F_half()
  {
    if (mpi_rank == 0)
      cout << "RansTurbV2F_half()" << endl;
    
    turbModel = V2F;

    C_MU  = getDoubleParam("C_MU",  "0.22");
    SIG_K = getDoubleParam("SIG_K", "1.0");
    SIG_D = getDoubleParam("SIG_D", "1.3");
    CEPS1 = getDoubleParam("CEPS1", "1.4");
    CEPS2 = getDoubleParam("CEPS2", "1.9");
    C1    = getDoubleParam("C1",    "1.4");
    C2    = getDoubleParam("C2",    "0.3");
    CETA  = getDoubleParam("CETA",  "70.0");
    CL    = getDoubleParam("CL",    "0.23");
    ALPHA = getDoubleParam("ALPHA", "0.6");
    ENN   = getDoubleParam("ENN",   "6.0");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine",  CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.7");
    eq->phiZero = 1.0e-9;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-06;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-06;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_D;

    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    turbTS = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
    turbLS = NULL;       registerScalar(turbLS, "turbLS", CV_DATA);
    muT    = NULL;       registerScalar(muT,    "muT",    CV_DATA);
    v2     = NULL;       registerScalar(v2,     "v2",     CV_DATA);

    tturb  = NULL; registerScalar(tturb, "tturb", CV_DATA);
    tkol   = NULL; registerScalar(tkol, "tkol", CV_DATA);
    trel   = NULL; registerScalar(trel, "trel", CV_DATA);
    lturb  = NULL; registerScalar(lturb, "lturb", CV_DATA);
    lkol   = NULL; registerScalar(lkol, "lkol", CV_DATA);
    lrel   = NULL; registerScalar(lrel, "lrel", CV_DATA);

    omega = NULL; registerScalar(omega, "omega", CV_DATA);
    eps_rhs = NULL; registerScalar(eps_rhs, "eps_rhs", CV_DATA);
  }

  virtual ~RansTurbV2F_half() {}

public:

  double *eps, *v2;                           ///< turbulent scalars, introduced to have access to variables, results in to more readable code
  double *kine_bfa, *eps_bfa;                 ///< turbulent scalars at the boundary
  double *muT;                                ///< turbulent viscosity at cell center for output
  
  double C_MU, SIG_K, SIG_D, CEPS1, CEPS2, C1, C2, CETA, CL, ALPHA, ENN;

  double *tturb, *tkol, *trel, *lkol, *lrel, *lturb;
  double *omega, *eps_rhs;

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook V2F model" << endl;

    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");     kine = eq->phi;      kine_bfa = eq->phi_bfa;
    eq = getScalarTransportData("eps");      eps = eq->phi;       eps_bfa = eq->phi_bfa;
    
    // initialize from *.in file
    double turb[2];
    Param *pmy;
    if (getParam(pmy, "INITIAL_CONDITION_TURB"))
    {
      turb[0] = pmy->getDouble(1);
      turb[1] = pmy->getDouble(2);
    }
    else
    {
      cerr << " Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
      throw(-1);
    }

    if (!checkScalarFlag("kine"))
      for (int icv = 0; icv < ncv; icv++)
        kine[icv] = turb[0];

    if (!checkScalarFlag("eps"))
      for (int icv = 0; icv < ncv; icv++)
        eps[icv] = 0.09*kine[icv]*omega[icv];

    if (!checkScalarFlag("v2"))
      for (int icv = 0; icv < ncv; icv++)
        v2[icv] = 2.0/3.0*kine[icv];

    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("eps", REPLACE_DATA);
    updateCvDataByName("v2", REPLACE_DATA);
  }

  virtual void calcRansTurbViscMuet()
  {
    calcGradVel();
    // update strain rate tensor 
    calcStrainRateAndDivergence();

    // update v2 variable
    for (int icv = 0; icv < ncv; icv++)
      v2[icv] = -1.0*rij_diag[icv][1]/rho[icv];

    //--------------------------------
    // calculate turbulent viscosity
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      double invwtot = 1.0/(w0+w1);

      double rho_fa = (w1*rho[icv0] + w0*rho[icv1])*invwtot;
      double kine_fa = (w1*kine[icv0] + w0*kine[icv1])*invwtot;
      double eps_fa  = (w1*eps[icv0]  + w0*eps[icv1])*invwtot;
      double v2_fa   = (w1*v2[icv0]   + w0*v2[icv1])*invwtot;
      double str_fa  = (w1*strMag[icv0] + w0*strMag[icv1])*invwtot;
      double nuLam_fa = (w1*calcMuLam(icv0)/rho[icv0] + w0*calcMuLam(icv1)/rho[icv1])*invwtot;

      double TS  = calcTurbTimeScale(kine_fa, eps_fa, v2_fa, str_fa, nuLam_fa, 1);
      mut_fa[ifa] = min(C_MU*rho_fa*v2_fa*TS, 10000.0);
      //mut_fa[ifa] = min(0.09*rho_fa*kine_fa*kine_fa/eps_fa, 10000.0);
    }

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))                             // if wall ensure nu_t = 0.0
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            mut_fa[ifa] = 0.0;
        }
        else
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)     // otherwise make first order extrapolation to face
          {
            int icv0 = cvofa[ifa][0];
            double TS = calcTurbTimeScale(kine[icv0], eps[icv0], v2[icv0], strMag[icv0], calcMuLam(icv0)/rho[icv0], 1);
            mut_fa[ifa] = min(C_MU*rho[icv0]*v2[icv0]*TS, 10000.0);
            //mut_fa[ifa] = min(0.09*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 10000.0);
          }
      }
    }

    for (int icv=0; icv<ncv; icv++)
    {
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
      turbTS[icv] = calcTurbTimeScale(kine[icv], eps[icv], v2[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);
      turbLS[icv] = calcTurbLengthScale(kine[icv], eps[icv], v2[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);

      tturb[icv] = kine[icv]/eps[icv];
      tkol[icv] = 6.0*sqrt(calcMuLam(icv)/rho[icv]/eps[icv]);
      trel[icv] = kine[icv]/max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-08);  //log(exp(sqrt(3.0)*v2[icv]*C_MU*strMag[icv]) + exp(1.0e-8));
      lturb[icv] = CL*pow(kine[icv],1.5)/eps[icv];
      lkol[icv] = CL*CETA*pow(calcMuLam(icv)/rho[icv],0.75)/pow(eps[icv],0.25);
      lrel[icv] = pow(kine[icv],1.5)/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-08));
    }
  }

  /**
   * calculate diffusivity scalars
   */
  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *scal;

    if (name == "kine")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }

    if (name == "eps")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }
  }

  /*inline double calcTurbTimeScale(const double &kine, const double &eps, const double &v2,
                                  const double &str, const double &nu, int realizable)
  {
    double TimeScale = max(kine/eps, 6.0*sqrt(nu/eps));
    if (realizable)
      TimeScale = min(TimeScale, kine/(sqrt(3.0)*v2*C_MU*str));
    return TimeScale;
  }

  inline double calcTurbLengthScale(const double &kine, const double &eps, const double &v2,
                                    const double &str, const double &nu, int realizable)
  {
    double LengthScale = CL*max(pow(kine,1.5)/eps, CETA*pow(nu,0.75)/pow(eps,0.25));
    if (realizable)
      LengthScale = min(LengthScale, pow(kine,1.5)/(sqrt(3.0)*v2*C_MU*str));
    return LengthScale;
  }*/

  inline double calcTurbTimeScale(const double &kine, const double &eps, const double &v2,
                                  const double &str, const double &nu, int realizable)
  {
    double TimeScale = kine/eps;
    double KolScale = 6.0*sqrt(nu/eps);
    TimeScale = max(TimeScale, KolScale);
    //TimeScale = sqrt(TimeScale*TimeScale + 36.0*nu/eps);

    if (realizable){
      double RealScale = kine/(max(sqrt(3.0)*v2*C_MU*str,1.0e-08));
      //if (RealScale > KolScale)
        TimeScale = min(TimeScale, RealScale);
      //else
      //  TimeScale = KolScale;
    }

    //TimeScale = max(TimeScale, 6.0*sqrt(nu/eps));

    return TimeScale;
  }

  inline double calcTurbLengthScale(const double &kine, const double &eps, const double &v2,
                                    const double &str, const double &nu, int realizable)
  {
    double LengthScale = CL*pow(kine,1.5)/eps;
    double KolScale = CL*CETA*pow(nu,0.75)/pow(eps,0.25);
    LengthScale = max(LengthScale, KolScale);
    //LengthScale = sqrt(LengthScale*LengthScale + CL*CL*CETA*CETA*pow(nu,1.5)/pow(eps,0.5));

    if (realizable){
      double RealScale = pow(kine,1.5)/(max(sqrt(3.0)*v2*C_MU*str,1.0e-08));
      //if (RealScale > KolScale)
        LengthScale = min(LengthScale, RealScale);
      //else
      //  LengthScale = KolScale;
    }

    //LengthScale = max(LengthScale, CL*CETA*pow(nu,0.75)/pow(eps,0.25));

    return LengthScale;
  }

  double getTurbProd(int icv, int realizable)
  {
    double mu_t = C_MU*rho[icv]*v2[icv]*calcTurbTimeScale(kine[icv], eps[icv], v2[icv], strMag[icv], calcMuLam(icv)/rho[icv], realizable);
    //double mu_t = 0.09*rho[icv]*kine[icv]*kine[icv]/eps[icv];
    return mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")     // nu_t*str*str - k/k*eps
    for (int icv = 0; icv < ncv; icv++)
    {
      double src  = getTurbProd(icv, 1)-rho[icv]*eps[icv];
      rhs[icv]   += src*cv_volume[icv];

      // rho*eps/(rho*kine)
      if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = -eps[icv]/kine[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
    }

    if (name == "eps")      // (ce1*getTurbProd-ce2*rho*eps)/TS
    for (int icv = 0; icv < ncv; icv++)
    {
      double TS  = calcTurbTimeScale(kine[icv], eps[icv], v2[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);
      double ce1 = CEPS1*(1.0 + 0.045*pow(kine[icv]/v2[icv], 0.5));

      double src = (ce1*getTurbProd(icv, 1) - rho[icv]*CEPS2*eps[icv])/TS;
      //double Pe = CEPS1*C_MU*rho[icv]*v2[icv]*strMag[icv]*strMag[icv] - 1./TS*CEPS1*2./3.*rho[icv]*kine[icv]*diverg[icv];
      //double src = Pe - rho[icv]*CEPS2*eps[icv]/TS;
      rhs[icv]  += src*cv_volume[icv];
      eps_rhs[icv] = rhs[icv];

      // d(ce2*rho*eps/TS)/d(rho*eps)
      if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = -CEPS2/TS;
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (name == "eps")
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int icv0 = cvofa[ifa][0];

        double nVec[3], s_half[3];
        normVec3d(nVec, fa_normal[ifa]);
        vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
        double wallDist = fabs(vecDotVec3d(s_half, nVec));
        double nuLamCV = calcMuLam(icv0)/rho[icv0];
        double epsWall = 2.0*nuLamCV*kine[icv0]/(wallDist*wallDist);

        phi_fa[ifa] = epsWall;
      }
    }
  }

  virtual void finalHookScalarRansTurbModel()
  {}
};


#endif
