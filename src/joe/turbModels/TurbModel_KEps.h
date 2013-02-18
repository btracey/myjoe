#ifndef RANSTURBMODEL_KEps_H
#define RANSTURBMODEL_KEps_H

#include "UgpWithCvCompFlow.h"



class RansTurbKEps : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbKEps()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKEps()" << endl;
    
    turbModel = KEPS;

    C_MU       = getDoubleParam("C_MU", "0.09");
    sigma_kine = getDoubleParam("sigma_kine", "1.0");
    sigma_eps  = getDoubleParam("sigma_eps", "1.3");
    ceps1      = getDoubleParam("ceps1",  "1.44");
    ceps2      = getDoubleParam("ceps2",  "1.92");
    CETA       = getDoubleParam("CETA",  "70.0");
    CL         = getDoubleParam("CL",    "0.23");
    
    realizable = getIntParam("REALIZABILITY", "0");
    
    if (mpi_rank == 0)
      cout << "REALIZABILITY: " << realizable << endl;

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "1.0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-3;
    eq->phiMaxiter = 100;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "1.0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-3;
    eq->phiMaxiter = 100;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;

    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    turbTS = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
    turbLS = NULL;       registerScalar(turbLS, "turbLS", CV_DATA);
    muT    = NULL;       registerScalar(muT, "muT", CV_DATA);
  }

  virtual ~RansTurbKEps() {}

public:   // member vars

  double *eps;                          ///< turbulent scalars, introduced to have access to variables, results in to more readable code
  double *kine_bfa, *eps_bfa;           ///< turbulent scalars at the boundary
  double *muT;                          ///< turbulent viscosity at cell center for output
  
  int realizable;   ///< switch for realizability: 0... no 
                    ///                            1... simple c_mu limiter
                    ///                            2... time scale limiter

  double C_MU, sigma_kine, sigma_eps, ceps1, ceps2, CETA, CL;

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook KEPS model" << endl;

    // connect pointers 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");    kine = eq->phi;     kine_bfa = eq->phi_bfa;    
    eq = getScalarTransportData("eps");     eps = eq->phi;      eps_bfa = eq->phi_bfa;

    double kineInit, epsInit;
    Param *pmy;
    if (getParam(pmy, "INITIAL_CONDITION_TURB"))
    {
      kineInit = pmy->getDouble(1);
      epsInit = pmy->getDouble(2);
    }
    else
    {
      cerr << " Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
      throw(-1);
    }

    if (!checkScalarFlag("kine"))
      for (int icv=0; icv<ncv; icv++)
        kine[icv] = kineInit;

    if (!checkScalarFlag("eps"))
      for (int icv=0; icv<ncv; icv++)
        eps[icv] = epsInit;

    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("eps", REPLACE_DATA);
  }

  virtual void calcRansTurbViscMuet()
  {
    calcGradVel();
    // update strain rate tensor 
    calcStrainRateAndDivergence();


    // internal faces
    for (int ifa=nfa_b; ifa<nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      double rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
      double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
      double eps_fa = (w1*eps[icv0] + w0*eps[icv1])/(w0+w1);

      if (realizable == 1)
      {
        double strain = (w1*strMag[icv0] + w0*strMag[icv1])/(w0+w1);
        double s = strain*kine_fa/eps_fa;
        double cmu = min(C_MU, sqrt(C_MU)/max(s,1.0e-12));
        mut_fa[ifa] = min(cmu*rho_fa*kine_fa*kine_fa/eps_fa, 1000.0);
      }
      else
        mut_fa[ifa] = min(C_MU*rho_fa*kine_fa*kine_fa/eps_fa, 1000.0);
    }

    // boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (zoneIsWall(zone->getName()))
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          mut_fa[ifa] = 0.0;                                     // set mut zero at walls
      else
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          
          if (realizable == 1)
          {
            double s = strMag[icv0]*kine[icv0]/eps[icv0];
            double cmu = min(C_MU, sqrt(C_MU)/max(s,1.0e-12));
            mut_fa[ifa] = min(cmu*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 1000.0);    // zero order extrapolation for others
          }
          else
            mut_fa[ifa] = min(C_MU*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 1000.0);    // zero order extrapolation for others
        }
    }

    // just for output 
    for (int icv=0; icv<ncv; icv++)
    {
      muT[icv]    = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
      turbTS[icv] = calcTurbTimeScale(kine[icv], eps[icv], strMag[icv], calcMuLam(icv)/rho[icv], realizable);
      turbLS[icv] = calcTurbLengthScale(kine[icv], eps[icv], strMag[icv], calcMuLam(icv)/rho[icv], realizable);
    }
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "kine")
    {
      eq = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_kine;
    }
    
    if (name == "eps")
    {
      eq = getScalarTransportData(name);
      for (int ifa=nfa_b; ifa<nfa; ifa++)
        eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_eps;
    }
  }
  
  inline double calcTurbTimeScale(const double &kine, const double &eps,
                                  const double &str, const double &nu, int realizable)
  {
    double TimeScale = max(kine/eps, 6.0*sqrt(nu/eps));
    if (realizable)
      TimeScale = min(TimeScale, sqrt(3.0)/(2.0*C_MU*str));
    return TimeScale;
  }

  inline double calcTurbLengthScale(const double &kine, const double &eps,
                                    const double &str, const double &nu, int realizable)
  {
    double LengthScale = CL*max(pow(kine,1.5)/eps, CETA*pow(nu,0.75)/pow(eps,0.25));
    if (realizable)
      LengthScale = min(LengthScale, sqrt(kine)*sqrt(3.0)/(2*C_MU*str));
    return LengthScale;
  }

  virtual double calcTurbProd(int icv)
  {    
    if (realizable == 1)
    {
      double s = strMag[icv]*kine[icv]/eps[icv];
      double cmu = min(C_MU, sqrt(C_MU)/s);
      double mu_t = cmu*rho[icv]*kine[icv]*kine[icv]/eps[icv];
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    }
    else
    {
      double mu_t = C_MU*rho[icv]*kine[icv]*kine[icv]/eps[icv];
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")     // Pk - rho*eps
      for (int icv=0; icv<ncv; icv++)
      {
        double src = calcTurbProd(icv)-rho[icv]*eps[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          double dsrcdphi = eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] += dsrcdphi*cv_volume[icv];
        }
      }

    if (name == "eps")    //
      for (int icv=0; icv<ncv; icv++)
      {
        double TS = calcTurbTimeScale(kine[icv], eps[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);
        double src = eps[icv]/kine[icv]*(ceps1*calcTurbProd(icv) - ceps2*rho[icv]*eps[icv]);
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          double dsrcdphi = ceps2*eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] += dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if ((param->getString() == "WALL") && (name == "eps"))
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv = cvofa[ifa][0];

            double nVec[3], s_half[3];
            normVec3d(nVec,fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);
            double wallDist = fabs(vecDotVec3d(s_half, nVec));
            double muLamCV = calcMuLam(icv);
            phi_fa[ifa] = 2.0*muLamCV*kine[icv]/(rho[icv]*wallDist*wallDist);
          }
        }
    }
  }


};



#endif








