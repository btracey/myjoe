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
    
    realizable = getIntParam("REALIZABILITY", "0");
    
    if (mpi_rank == 0)
      cout << "REALIZABILITY: " << realizable << endl;

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "1.0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "1.0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;

    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
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



  double C_MU, sigma_kine, sigma_eps, ceps1, ceps2;

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
        double cmu = min(C_MU, sqrt(C_MU)/s);
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
            double cmu = min(C_MU, sqrt(C_MU)/s);
            mut_fa[ifa] = min(cmu*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 1000.0);    // zero order extrapolation for others
          }
          else
            mut_fa[ifa] = min(C_MU*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 1000.0);    // zero order extrapolation for others
        }
    }

    // just for output 
    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
//    {
//      muT[icv] = 0.0;
//
//      int foc_f = faocv_i[icv];
//      int foc_l = faocv_i[icv + 1] - 1;
//      for (int foc=foc_f; foc<=foc_l; foc++)
//      {
//        int ifa = faocv_v[foc];
//        muT[icv] += mut_fa[ifa];
//      }
//      muT[icv] /= (double)(foc_l-foc_f+1);
//    }
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "kine")
    {
      eq = getScalarTransportData(name);
      
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

        eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_kine;
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_kine;
          }
      }
    }
    
    if (name == "eps")
    {
      eq = getScalarTransportData(name);
      
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

        eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_eps;
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/sigma_eps;
          }
      }
    }
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

  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)
  {
    if (name == "kine")     // Pk - rho*eps
    for (int icv=0; icv<ncv; icv++)
    {
      double src = calcTurbProd(icv);
      rhs[icv] += src*cv_volume[icv];

      if (A != NULL)
      {
        double dsrcdphi = rho[icv]*eps[icv]/kine[icv];
        int noc00 = nbocv_i[icv];
        A[noc00] += dsrcdphi*cv_volume[icv];
      }
    }

    if (name == "eps")    // 
    for (int icv=0; icv<ncv; icv++)
    {
      double src = eps[icv]/kine[icv]*ceps1*calcTurbProd(icv);
      rhs[icv] += src*cv_volume[icv];


      if (A != NULL)
      {
        double dsrcdphi = ceps2*rho[icv]*eps[icv]/kine[icv];
        int noc00 = nbocv_i[icv];
        A[noc00] += dsrcdphi*cv_volume[icv];
      }
    }
  }
  
/*  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)
  {
    if (name == "kine")     // Pk - rho*eps
    for (int icv=0; icv<ncv; icv++)
    {
      double src = 0.0;
      if (SinhaFix == 1)
      {
        double b1Inf = 0.4;
        double b1Prime = b1Inf*(1.0-exp(1.0-Mach1));
        src = -2.0/3.0*rho[icv]*kine[icv]*(grad_u[icv][0][0])*(1.0-b1Prime);
      }
      else if (SinhaFix == 2)
      {
        double b1Inf = 0.4;
        double b1Prime = b1Inf*(1.0-exp(1.0-Mach1));
        src = calcTurbProd(icv)*(1.0-b1Prime);
      }
      else
        src = calcTurbProd(icv);
      
      rhs[icv] += src*cv_volume[icv];

      double dsrcdphi = rho[icv]*eps[icv]/kine[icv];
      int noc00 = nbocv_i[icv];
      A[noc00] += dsrcdphi*cv_volume[icv];
    }

    if (name == "eps")    // 
    for (int icv=0; icv<ncv; icv++)
    {
      double prod = calcTurbProd(icv);
      double src = 0.0;
      
      if ((SinhaFix == 1) || (SinhaFix == 2))
      {
        double ceps1Mod = 1.25 + 0.2*(Mach1-1.0);
        src = eps[icv]/kine[icv]*ceps1Mod*prod;
      }
      else
        src = eps[icv]/kine[icv]*ceps1*prod;
            
      rhs[icv] += src*cv_volume[icv];

      double dsrcdphi = ceps2*rho[icv]*eps[icv]/kine[icv];
      int noc00 = nbocv_i[icv];
      A[noc00] += dsrcdphi*cv_volume[icv];
    }
  }*/

  virtual void sourceHookScalarRansTurbExpl(double *rhs, const string &name)
  {
    if (name == "kine")     // Pk - rho*eps
    for (int icv=0; icv<ncv; icv++)
    {
      double src = calcTurbProd(icv)-rho[icv]*eps[icv];
      rhs[icv] += src*cv_volume[icv];
    }

    if (name == "eps")    // 
    for (int icv=0; icv<ncv; icv++)
    {
      double src = eps[icv]/kine[icv]*(ceps1*calcTurbProd(icv) - ceps2*rho[icv]*eps[icv]);
      rhs[icv] += src*cv_volume[icv];
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")     // Pk - rho*eps
    for (int icv=0; icv<ncv; icv++)
    {
      double src = calcTurbProd(icv)-rho[icv]*eps[icv];
      rhs[icv] += src*cv_volume[icv];
    }

    if (name == "eps")    //
    for (int icv=0; icv<ncv; icv++)
    {
      double src = eps[icv]/kine[icv]*(ceps1*calcTurbProd(icv) - ceps2*rho[icv]*eps[icv]);
      rhs[icv] += src*cv_volume[icv];
    }

    if (A!=NULL)
    {
      if (mpi_rank==0)
        cerr << "implizit part for sourceHookScalarRansTurb_new in keps not implemented yet!" << endl;
      throw(-1);
    }
  }
};



#endif








