#ifndef RANSTURBMODEL_KOMSST_H
#define RANSTURBMODEL_KOMSST_H

#include "UgpWithCvCompFlow.h"



class RansTurbKOmSST : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbKOmSST()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOmSST()" << endl;

    turbModel = KOMSST;

    sigma_k1  = getDoubleParam("sigma_k1",   "0.85");
    sigma_k2  = getDoubleParam("sigma_k2",   "1.0");
    sigma_om1 = getDoubleParam("sigma_om1",  "0.5");
    sigma_om2 = getDoubleParam("sigma_om2",  "0.856");
    beta_1    = getDoubleParam("beta_1",     "0.075");
    beta_2    = getDoubleParam("beta_2",     "0.0828");
    betaStar  = getDoubleParam("betaStar",   "0.09");
    a1        = getDoubleParam("a1", "0.31");
    alfa_1 = beta_1/betaStar - sigma_om1*pow(0.41, 2.0)/sqrt(betaStar);
    alfa_2 = beta_2/betaStar - sigma_om2*pow(0.41, 2.0)/sqrt(betaStar);


    SST_limitPK = getIntParam("SST_LIMIT_PK", "1");  // limiter for turbulent production: 0 ... no limiter
                                                     //                                   1 ... with Pk=min(Pk, 10*betastar*rho*kine*omega);

    SST_limiterShearStress = getIntParam("SST_LIMITER_SHEARSTRESS", "0");  // limit shear stress with: 0 ... vort magnitude,
                                                                           //                          1 ... strain rate magnitude

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.6");
    eq->phiZero = getDoubleParam("ZERO_kine", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_kine", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = getDoubleParam("ZERO_omega", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_omega", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    strMag      = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    vortMag     = NULL;       registerScalar(vortMag, "vortMag", CV_DATA);
    diverg      = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT         = NULL;       registerScalar(muT, "muT", CV_DATA);
    crossDiff   = NULL;       registerScalar(crossDiff, "crossDiff", CV_DATA);
    blendFuncF1 = NULL;       registerScalar(blendFuncF1, "blendFuncF1", CV_DATA);
    blendFuncF2 = NULL;       registerScalar(blendFuncF2, "blendFuncF2", CV_DATA);
    wallDist    = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);
    wallConn    = NULL;       // array of integers
    turbTS      = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
    turbLS      = NULL;       registerScalar(turbLS, "turbLS", CV_DATA);
  }

  virtual ~RansTurbKOmSST() {}

public:   // member vars

  double *omega;                                    ///< specific dissipation, introduced to have access to variables, results in to more readable code
  double (*grad_kine)[3], (*grad_omega)[3];
  double *kine_bfa, *omega_bfa;                     ///< turbulent scalars at the boundary
  double *muT;                                      ///< turbulent viscosity at cell center for output
  double *wallDist;                                 ///< distance to closest wall face
  int *wallConn;                                    ///< index of closest wall face
  double *crossDiff, *blendFuncF1, *blendFuncF2;
  
  // limiter for shear stress
  int SST_limiterShearStress;
  double *limiterFunc;

  double sigma_k1, sigma_k2, sigma_om1, sigma_om2, alfa_1, alfa_2, beta_1, beta_2, betaStar, a1;

  // limiter for Pk
  int SST_limitPK;

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook KOM SST model" << endl;

    wallConn = new int[ncv];
    calcWallDistance(wallConn, wallDist);

    if (SST_limiterShearStress == 1)      limiterFunc = strMag;   // shear stress limiter = strain rate magnitude
    else                                  limiterFunc = vortMag;  // shear stress limiter = vorticity magnitude
      
    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");      kine = eq->phi;    kine_bfa = eq->phi_bfa;   grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");     omega = eq->phi;   omega_bfa = eq->phi_bfa;  grad_omega = eq->grad_phi;
    

    // set initial condition if parameter given in Joe.in

    double kineInit, omegaInit;
    if (checkParam("INITIAL_CONDITION_TURB"))
    {
      kineInit = getParam("INITIAL_CONDITION_TURB")->getDouble(1);
      omegaInit = getParam("INITIAL_CONDITION_TURB")->getDouble(2);
    }
    else
    {
      cout << "Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
      throw(-1);
    }

    if (!checkScalarFlag("kine"))
      for (int icv=0; icv<ncv; icv++)
        kine[icv] = kineInit;

    if (!checkScalarFlag("omega"))
      for (int icv=0; icv<ncv; icv++)
        omega[icv] = omegaInit;


    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("omega", REPLACE_DATA);
  }
  
  virtual void calcRansTurbViscMuet()
  {
    calcGradVel();
    
    // update strain rate tensor 
    calcStrainRateAndDivergence();
    calcVorticity();    
    
    calcMenterBlendingFunctions();

    
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
      double om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);
      double limitFunc_fa = (w1*limiterFunc[icv0] + w0*limiterFunc[icv1])/(w0+w1);
      double blendFuncF2_fa = (w1*blendFuncF2[icv0] + w0*blendFuncF2[icv1])/(w0+w1);

      double zeta = min(1.0/om_fa, a1/(limitFunc_fa*blendFuncF2_fa));

      mut_fa[ifa] = min(max(rho_fa*kine_fa*zeta, 0.0), 1.0);
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

            double zeta = min(1.0/omega[icv0], a1/(limiterFunc[icv0]*blendFuncF2[icv0]));
            mut_fa[ifa] = min(max(rho[icv0]*kine[icv0]*zeta, 0.0), 1.0);
          }
      }
    
    // just for output and visulization
    for (int icv=0; icv<ncv; icv++)
    {
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
      turbTS[icv] = calcTurbTimeScale(kine[icv], omega[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);
      turbLS[icv] = calcTurbLengthScale(kine[icv], omega[icv], strMag[icv], calcMuLam(icv)/rho[icv], 1);
    }
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

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_k1*f1_fa + (1. - f1_fa)*sigma_k2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
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
            double coeff = sigma_k1*blendFuncF1[icv0] + (1. - blendFuncF1[icv0])*sigma_k2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }

    if (name == "omega")
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

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_om1*f1_fa + (1. - f1_fa)*sigma_om2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
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
            double coeff = sigma_om1*blendFuncF1[icv0] + (1. - blendFuncF1[icv0])*sigma_om2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0e5);
        double Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
        if (SST_limitPK == 1)
          Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);

        Pk = max(Pk, 0.0);  // for stabilization
        double src = Pk - betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - betaStar*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }

    if (name == "omega")
    {
      for (int icv=0; icv<ncv; icv++)
      {
        double F1 = blendFuncF1[icv];
        double alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
        double beta = F1*beta_1 + (1.0 - F1)*beta_2;

        double zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
        double Pk = strMag[icv]*strMag[icv] - 2./3.*zeta*diverg[icv];

        Pk = max(Pk, 0.0);  // for stabilization
        double src = alfa*rho[icv]*Pk + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");
    
    for (int icv = 0; icv < ncv; icv++)
    {
      double F1 = blendFuncF1[icv];
      double alfa = F1 * alfa_1 + (1.0 - F1) * alfa_2;
      double beta = F1 * beta_1 + (1.0 - F1) * beta_2;

      double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
      double mut  = min(max(rho[icv] * kine[icv] * zeta, 1.0e-8), 1.0e5);

      double Pk = max(mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
      if (SST_limitPK == 1)
        Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);

      rhs[icv][5+kine_Index]  += (Pk - betaStar*rho[icv]*omega[icv]*kine[icv])*cv_volume[icv];

      zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
      Pk = max(strMag[icv]*strMag[icv] - 2./3.*zeta*diverg[icv], 0.0);
      rhs[icv][5+omega_Index] += (alfa*rho[icv]*Pk + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv])*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
//        A[noc00][5+kine_Index][0]              -=   betaStar * kine[icv] * omega[icv] * cv_volume[icv];
//        A[noc00][5+kine_Index][5+omega_Index]  -= - betaStar * kine[icv] * cv_volume[icv];
        A[noc00][5+kine_Index][5+kine_Index]   -= - betaStar*omega[icv]*cv_volume[icv];
//        A[noc00][5+omega_Index][0]             -=   beta * omega[icv] * omega[icv] * cv_volume[icv];
        A[noc00][5+omega_Index][5+omega_Index] -= - 2.0*beta*omega[icv]*cv_volume[icv];
     }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if ((param->getString() == "WALL") && (name == "omega"))
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv = cvofa[ifa][0];
            double muLamCV = calcMuLam(icv);
            phi_fa[ifa] = 60.0*muLamCV/(rho[icv]*beta_1*wallDist[icv]*wallDist[icv]);
          }
        }
    }
  }

  virtual void calcMenterBlendingFunctions()
  {
    calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
    calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

    for (int icv=0; icv<ncv; icv++)
    {
      double d = wallDist[icv];
      double mue = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);

      crossDiff[icv] = max(2.0*rho[icv]*sigma_om2/omega[icv]*vecDotVec3d(grad_kine[icv], grad_omega[icv]),  1.0e-20);
//      crossDiff[icv] = min(crossDiff[icv], 10000.0);

      double gamma1 = 500.0*mue/(pow(d, 2.0)*rho[icv]*omega[icv]);
      double gamma2 = 4.0*sigma_om2*rho[icv]*kine[icv]/(d*d*crossDiff[icv]);
      double gamma3 = sqrt(kine[icv])/(betaStar*omega[icv]*d);

      double gamma = min(max(gamma1, gamma3), gamma2);
      blendFuncF1[icv] = tanh(pow(gamma,4.0));
      gamma = max(2.0*gamma3, gamma1);
      blendFuncF2[icv] = tanh(pow(gamma,2.0));
    }

    updateCvData(crossDiff, REPLACE_DATA);
    updateCvData(blendFuncF1, REPLACE_DATA);
    updateCvData(blendFuncF2, REPLACE_DATA);
  }

  inline double calcTurbTimeScale(const double &kine, const double &omega,
                                  const double &str, const double &nu, int realizable)
  {
    double eps = betaStar*kine*omega;
    double TimeScale = max(kine/eps, 2.0*sqrt(nu/eps));  //v2f has 6.0 instead of 2.0
    if (realizable)
      TimeScale = min(TimeScale, 1.0/(sqrt(6.0)*betaStar*str)); //this is not consistent with v2f
    return TimeScale;
  }

  inline double calcTurbLengthScale(const double &kine, const double &omega,
                                    const double &str, const double &nu, int realizable)
  {
    double CL = 0.23, CETA = 70.0;
    double eps = betaStar*kine*omega;
    double LengthScale = CL*max(pow(kine,1.5)/eps, CETA*pow(nu,0.75)/pow(eps,0.25));
    if (realizable)
      LengthScale = min(LengthScale, pow(kine,0.5)/(2.0/3.0*sqrt(3.0)*betaStar*str));
    return LengthScale;
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (wallConn != NULL) {delete [] wallConn; wallConn = NULL;}

    if (checkParam("CALC_RS_BOUSS"))
    {
      if (mpi_rank == 0)
        cout << "calculating Rs Bouss" << endl;
      for (int icv = 0; icv < ncv; icv++)
      {
        double term1 = (2.0/3.0) * rho[icv] * kine[icv];
        rij_diag[icv][0] = -term1 + muT[icv] * 2.0 * (grad_u[icv][0][0] - 1.0/3.0*diverg[icv]);
        rij_diag[icv][1] = -term1 + muT[icv] * 2.0 * (grad_u[icv][1][1] - 1.0/3.0*diverg[icv]);
        rij_diag[icv][2] = -term1 + muT[icv] * 2.0 * (grad_u[icv][2][2] - 1.0/3.0*diverg[icv]);
        rij_offdiag[icv][0] =     + muT[icv] * (grad_u[icv][0][1] + grad_u[icv][1][0]);
        rij_offdiag[icv][1] =     + muT[icv] * (grad_u[icv][0][2] + grad_u[icv][2][0]);
        rij_offdiag[icv][2] =     + muT[icv] * (grad_u[icv][1][2] + grad_u[icv][2][1]);
      }
    }
  }

};



#endif


