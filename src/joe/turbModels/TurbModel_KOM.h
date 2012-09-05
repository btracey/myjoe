#ifndef RANSTURBMODEL_KOM_H
#define RANSTURBMODEL_KOM_H

#include "UgpWithCvCompFlow.h"



/**
 * wilcox k-omega model as in Turbulence Modeling for CFD, Third Edition, Wilcox 2006
 */

class RansTurbKOm : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbKOm()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOm()" << endl;
    
    turbModel = KOM;

    betaStar   = getDoubleParam("betaStar",   "0.09");
    sigmaOmega = getDoubleParam("sigmaOmega", "0.5");
    sigmaStar  = getDoubleParam("sigmaStar",  "0.6");
    alfa       = getDoubleParam("alfa",       "0.52");
    sigmad0    = getDoubleParam("sigmad0",    "0.125");
    beta0      = getDoubleParam("beta0",      "0.0708");
    cLim       = getDoubleParam("cLim",       "0.875");

    KOM_RealizableConstraint = getIntParam("KOM_RealizableConstraint",  "1");   // this is the one from Wilcox 2006

    ScalarTranspEq *eq;

    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = 1.0e-6;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    strMag   = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    vortMag  = NULL;       registerScalar(vortMag, "vortMag", CV_DATA);
    diverg   = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT      = NULL;       registerScalar(muT, "muT", CV_DATA);
    wallDist = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);
  }

  virtual ~RansTurbKOm() {}

public:   // member vars

  double *kine_bfa, (*grad_kine)[3];
  double *omega, *omega_bfa, (*grad_omega)[3];
  double *muT;                      ///< turbulent viscosity at cell center for output
  double *wallDist;                 ///< wall distance
  
  double betaStar, sigmaOmega, sigmaStar, alfa, sigmad0, beta0, cLim;
  int KOM_RealizableConstraint;

public:   // member functions
  
  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook WILCOX KOM model" << endl;

    if (!checkParam("DO_NOT_CALC_WALLDIST"))
    if (!checkDataFlag(wallDist))
    {
      for (int icv=0; icv<ncv; icv++)
        wallDist[icv] = 0.0;

      calcWallDistance(NULL, wallDist);
    }

    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");    kine = eq->phi;     kine_bfa = eq->phi_bfa;   grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");   omega = eq->phi;    omega_bfa = eq->phi_bfa;  grad_omega = eq->grad_phi;

    double kineInit, omegaInit;
    Param *pmy;
    if (getParam(pmy, "INITIAL_CONDITION_TURB"))
    {
      kineInit = pmy->getDouble(1);
      omegaInit = pmy->getDouble(2);
    }
    else
    {
      if (mpi_rank == 0)
        cerr << " Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
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
      double strMag_fa = (w1*strMag[icv0] + w0*strMag[icv1])/(w0+w1);

      
      if (KOM_RealizableConstraint == 1)
      {
        double omega_tilde = max(om_fa, cLim*strMag_fa/sqrt(betaStar));
        mut_fa[ifa] = min(rho_fa*kine_fa/omega_tilde, 100.0);
      }
      else if (KOM_RealizableConstraint == 2)
      {
        double TS = min(1.0/om_fa, 0.6/(sqrt(6.0)*strMag_fa));
        mut_fa[ifa] = rho_fa*kine_fa*TS;
      }
      else
      {
        mut_fa[ifa] = min(rho_fa*kine_fa/om_fa, 100.0);
      }
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
          
          if (KOM_RealizableConstraint == 1)
          {
            double omega_tilde = max(omega[icv0], cLim*strMag[icv0]/sqrt(betaStar));
            mut_fa[ifa] = min(rho[icv0]*kine[icv0]/omega_tilde, 100.0);    // zero order extrapolation for others
          }
          else if (KOM_RealizableConstraint == 2)
          {
            double TS = min(1.0/omega[icv0], 0.6/(sqrt(6.0)*strMag[icv0]));
            mut_fa[ifa] = rho[icv0]*kine[icv0]*TS;
          }
          else
          {
            mut_fa[ifa] = min(rho[icv0]*kine[icv0]/omega[icv0], 100.0);    // zero order extrapolation for others
          }
        }
    }

    // just calculate mut in cell center by averaging over all face values
    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }



  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if ((name == "kine") || (name == "omega"))
    {
      double sigma;
      if (name == "kine")       sigma = sigmaStar;
      if (name == "omega")      sigma = sigmaOmega;

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

        double rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
        double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
        double om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);

        eq->diff[ifa] = mul_fa[ifa] + sigma*rho_fa*kine_fa/om_fa;//mut_fa[ifa];//
        //eq->diff[ifa] = mul_fa[ifa] + sigma*mut_fa[ifa];
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];            // set mut zero at walls
          }
        }
        else
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa] + sigma*rho_bfa[ifa]*kine_bfa[ifa]/omega_bfa[ifa]; //mut_fa[ifa];//
            //eq->diff[ifa] = mul_fa[ifa] + sigma*rho[icv0]*kine[icv0]/omega[icv0];//mut_fa[ifa];//
            //eq->diff[ifa] = mul_fa[ifa] + sigma*mut_fa[ifa];//
          }
        }
      }
    }
  }

  virtual double calcTurbProd(int icv)
  {
    if (KOM_RealizableConstraint == 1)
    {
      double omega_tilde = max(omega[icv], cLim*strMag[icv]/sqrt(betaStar));
      double mu_t = rho[icv]*kine[icv]/omega_tilde;
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    }
    else if (KOM_RealizableConstraint == 2)
    {
      double TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
      double mu_t = rho[icv]*kine[icv]*TS;
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    }
    else
    {
      double mu_t = rho[icv]*kine[icv]/omega[icv];
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    }
  }

  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)
  {
    if (name == "kine")
    for (int icv=0; icv<ncv; icv++)
    {
      rhs[icv] += calcTurbProd(icv)*cv_volume[icv];

      int noc00 = nbocv_i[icv];
      double dsrcdphi = - betaStar*rho[icv]*omega[icv];
      A[noc00] -= dsrcdphi*cv_volume[icv];
    }

    if (name == "omega")
    {
      double OM[3][3], STR_hat[3][3];

      for (int icv=0; icv<ncv; icv++)
      {
        double sigmad; 
        double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);
        
        double TS = 1.0;
        if (KOM_RealizableConstraint == 2)
        {
          TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
        }

        if (crossDiff <= 0.0)   sigmad = 0.0;
        else                    sigmad = sigmad0;

        double src = TS*alfa*omega[icv]/kine[icv]*calcTurbProd(icv) + sigmad*rho[icv]/omega[icv]*crossDiff;
        rhs[icv] += src*cv_volume[icv];

        double divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);
          
          // strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

        double chiOm = 0.0;
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));
        
        double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        double beta = beta0*fbeta;

        int noc00 = nbocv_i[icv];
        double dsrcdphi = -TS*beta*rho[icv]*omega[icv];
        A[noc00] -= dsrcdphi*cv_volume[icv];
      }
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
    for (int icv=0; icv<ncv; icv++)
    {
      double src = calcTurbProd(icv) - betaStar*rho[icv]*omega[icv]*kine[icv];
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
      double OM[3][3], STR_hat[3][3];

      for (int icv=0; icv<ncv; icv++)
      {

        double TS = 1.0;
        if (KOM_RealizableConstraint == 2)
          TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));

        double sigmad;
        double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);

        if (crossDiff <= 0.0)   sigmad = 0.0;
        else                    sigmad = sigmad0;

        double divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

          // strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

        double chiOm = 0.0;
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

        double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        double beta = beta0*fbeta;

        double src =  TS*alfa*omega[icv]/kine[icv]*calcTurbProd(icv)
                    - TS*beta*rho[icv]*omega[icv]*omega[icv]
                    + sigmad*rho[icv]/omega[icv]*crossDiff;

        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = -2.0*TS*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");
    
    double OM[3][3], STR_hat[3][3];

    for (int icv = 0; icv < ncv; icv++)
    {
      double TS = 1.0;
      if (KOM_RealizableConstraint == 2)
        TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));

      double sigmad;
      double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);

      if (crossDiff <= 0.0)   sigmad = 0.0;
      else                    sigmad = sigmad0;

      double divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

        // strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
        if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
        else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
      }

      double chiOm = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
        chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
      chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

      double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
      double beta = beta0*fbeta;

      double src =  TS*alfa*omega[icv]/kine[icv]*calcTurbProd(icv)
                  - TS*beta*rho[icv]*omega[icv]*omega[icv]
                  + sigmad*rho[icv]/omega[icv]*crossDiff;

      rhs[icv][5+kine_Index]  += (calcTurbProd(icv) - betaStar * rho[icv] * omega[icv] * kine[icv]) * cv_volume[icv];
      rhs[icv][5+omega_Index] += src * cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
        A[noc00][5+kine_Index][5+kine_Index]   -= - betaStar * omega[icv] * cv_volume[icv];
        A[noc00][5+omega_Index][5+omega_Index] -= - 2.0 * TS * beta * omega[icv] * cv_volume[icv];
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
            phi_fa[ifa] = 6.0*muLamCV/(rho[icv]*beta0*wallDist[icv]*wallDist[icv]);
          }
        }
    }
  }

};



#endif
