#ifndef RANSTURBMODEL_ASBMKEps_H
#define RANSTURBMODEL_ASBMKEps_H

#include "UgpWithCvCompFlow.h"

class RansTurbASBMKEps : virtual public UgpWithCvCompFlow
{
public:
  RansTurbASBMKEps()
  {
    if (mpi_rank == 0)
      cout << "RansTurbASBMKEps()" << endl;

    turbModel = ASBM;

    C_MU  = getDoubleParam("C_MU",  "0.22");
    SIG_K = getDoubleParam("SIG_K", "1.0");
    SIG_D = getDoubleParam("SIG_D", "1.3");
    CEPS1 = getDoubleParam("CEPS1", "1.4");
    CEPS2 = getDoubleParam("CEPS2", "1.9");
    CETA  = getDoubleParam("CETA",  "70.0");
    CL    = getDoubleParam("CL",    "0.23");

    LIMIT_PK     = getIntParam("LIMIT_PK",     "0");
    VEL_SCALE    = getIntParam("VEL_SCALE",    "0");
    RIJ_BASED_PK = getIntParam("RIJ_BASED_PK", "0");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine",  CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.7");
    eq->phiZero = 1.0e-10;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_D;

    v2     = NULL;       registerScalar(v2,     "v2"    , CV_DATA);
    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    turbTS = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
    turbLS = NULL;       registerScalar(turbLS, "turbLS", CV_DATA);
    muT    = NULL;       registerScalar(muT,    "muT",    CV_DATA);
    tturb  = NULL; registerScalar(tturb, "tturb", CV_DATA);
    tkol   = NULL; registerScalar(tkol, "tkol", CV_DATA);
    trel   = NULL; registerScalar(trel, "trel", CV_DATA);
    lturb  = NULL; registerScalar(lturb, "lturb", CV_DATA);
    lkol   = NULL; registerScalar(lkol, "lkol", CV_DATA);
    lrel   = NULL; registerScalar(lrel, "lrel", CV_DATA);
  }

  virtual ~RansTurbASBMKEps() {}

public:

  double *eps, *v2;                                     ///< introduced to have access to variables, results into more readable code
  double *kine_bfa, *eps_bfa;                           ///< turbulent scalars at the boundary
  double *muT;                                          ///< turbulent viscosity at cell center for output

  double C_MU, SIG_K, SIG_D, CEPS1, CEPS2, CETA, CL;

  double *tturb, *tkol, *trel, *lkol, *lrel, *lturb;

  int LIMIT_PK, VEL_SCALE, RIJ_BASED_PK;

public:
  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHook() for ASBMKEps" << endl;

    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");     kine = eq->phi;      kine_bfa = eq->phi_bfa;
    eq = getScalarTransportData("eps");      eps = eq->phi;       eps_bfa = eq->phi_bfa;
  }

  virtual void calcRansTurbViscMuet()
  {
    // compute wall-normal stress
    if(VEL_SCALE==0)
    	for (int icv = 0; icv < ncv; icv++)  {
      double fluc_norm = fabs(rij_diag[icv][1]);
      v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
    }

    if(VEL_SCALE==1)
    	for (int icv = 0; icv < ncv; icv++)  {
      double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
      double ny = vel[icv][1]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
      double fluc_norm = fabs(rij_diag[icv][0]*ny-rij_diag[icv][1]*nx);
      v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
    }

    if(VEL_SCALE==11)
    	for (int icv = 0; icv < ncv; icv++)  {
      double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
      double ny = vel[icv][1]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
      double fluc_norm = fabs(rij_diag[icv][0]*ny*ny+rij_diag[icv][1]*nx*nx-2*nx*ny*rij_offdiag[icv][0]);
      v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
    }

    if(VEL_SCALE==2)
    	for (int icv = 0; icv < ncv; icv++)  {
      double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][2]*vel[icv][2]+1.e-12);
      double nz = vel[icv][2]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][2]*vel[icv][2]+1.e-12);
      double fluc_norm = fabs(rij_diag[icv][0]*nz-rij_diag[icv][2]*nx);
      //double fluc_norm = fabs(rij_diag[icv][0]*nz*nz+rij_diag[icv][2]*nx*nx-2*nx*nz*rij_offdiag[icv][1]);
      v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
    }

    if(VEL_SCALE==3)
    	for (int icv = 0; icv < ncv; icv++)  {
      v2[icv] = 2./3.*kine[icv];
    }

    updateCvData(v2, REPLACE_DATA);

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

      double rho_fa    = (w1*rho[icv0] + w0*rho[icv1])*invwtot;
      double turbTS_fa = (w1*turbTS[icv0] + w0*turbTS[icv1])*invwtot;
      double v2_fa     = (w1*v2[icv0] + w0*v2[icv1])*invwtot;

      mut_fa[ifa] = min(max(C_MU*rho_fa*v2_fa*turbTS_fa, 0.0), 1.0);
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
            mut_fa[ifa] = min(max(C_MU*rho[icv0]*v2[icv0]*turbTS[icv0], 0.0), 1.0);
          }
      }
    }

    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
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

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv = 0; icv < ncv; icv++)
      {
        double src  = getTurbProd(icv) - rho[icv]*eps[icv];
        rhs[icv] += src*cv_volume[icv];

        //d(rho*kine*eps/kine)/d(rho*kine)
        if (flagImplicit)
        {
          double dsrcdphi = -eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    if (name == "eps")
      for (int icv = 0; icv < ncv; icv++)
      {
        // another good value is 0.0175
        double ce1 = CEPS1*(1.0 + 0.045*pow(kine[icv]/v2[icv], 0.5));

        double src = (ce1*getTurbProd(icv) - CEPS2*rho[icv]*eps[icv])/turbTS[icv];
        rhs[icv]  += src*cv_volume[icv];

        // d(ce2*rho*eps/TS)/d(rho*eps)
        if (flagImplicit)
        {
          double dsrcdphi = -CEPS2/turbTS[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = 5+getScalarTransportIndex("kine");
    int eps_Index =  5+getScalarTransportIndex("eps");

    double dsrcdphi;

    // mu_t*str*str - rho*k/k*eps
    for (int icv = 0; icv < ncv; icv++)
    {
      double kine_src  = getTurbProd(icv)-rho[icv]*eps[icv];
      double ce1 = CEPS1*(1.0 + 0.045*pow(kine[icv]/v2[icv], 0.5));
      double eps_src = (ce1*getTurbProd(icv) - CEPS2*rho[icv]*eps[icv])/turbTS[icv];

      rhs[icv][kine_Index] += kine_src*cv_volume[icv];
      rhs[icv][eps_Index]  += eps_src*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];

        dsrcdphi = -eps[icv]/kine[icv];
        A[noc00][kine_Index][kine_Index] -= dsrcdphi*cv_volume[icv];
        dsrcdphi = -1.;
        A[noc00][kine_Index][eps_Index]  -= dsrcdphi*cv_volume[icv];

        dsrcdphi = -CEPS2/turbTS[icv];
        A[noc00][eps_Index][eps_Index] -= dsrcdphi*cv_volume[icv];
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

  void calcTurbTimeScale()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];
      double TimeScale = max(kine[icv]/eps[icv], 6.0*sqrt(nu/eps[icv]));
      tturb[icv] = kine[icv]/eps[icv];
      tkol[icv] = 6.0*sqrt(nu/eps[icv]);
      trel[icv] = kine[icv]/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-14));
      bool realizable = true;
      if (realizable)
      {
        double RealScale = kine[icv]/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-14));
        TimeScale = min(TimeScale, RealScale);
      }

      turbTS[icv] = TimeScale;
    }
    updateCvData(turbTS,REPLACE_DATA);
  }

  void calcTurbLengthScale()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];
      double LengthScale = CL*max(pow(kine[icv],1.5)/eps[icv], CETA*pow(nu,0.75)/pow(eps[icv],0.25));

      bool realizable = true;
      if (realizable)
      {
        double RealScale = pow(kine[icv],1.5)/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-14));
        LengthScale = min(LengthScale, RealScale);
      }

      turbLS[icv] = LengthScale;
    }
    updateCvData(turbLS,REPLACE_DATA);
  }

  double getTurbProd(int icv)
  {
    double Pk, mu_t;

    switch(RIJ_BASED_PK)
    {
    case 0:
      mu_t = min(max(C_MU*rho[icv]*v2[icv]*turbTS[icv], 0.0), 1.0);
      Pk = mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
      break;
    case 1:
      Pk = rij_diag[icv][0]*grad_u[icv][0][0]    + rij_offdiag[icv][0]*grad_u[icv][0][1] + rij_offdiag[icv][1]*grad_u[icv][0][2] +
           rij_offdiag[icv][0]*grad_u[icv][1][0] + rij_diag[icv][1]*grad_u[icv][1][1]    + rij_offdiag[icv][2]*grad_u[icv][1][2] +
           rij_offdiag[icv][1]*grad_u[icv][2][0] + rij_offdiag[icv][2]*grad_u[icv][2][1] + rij_diag[icv][2]*grad_u[icv][2][2];
      break;
    }

    Pk = max(Pk,0.0);
    if (LIMIT_PK == 1)
      Pk = min(Pk, 20.0*rho[icv]*eps[icv]);
    return Pk;
  }

};

#endif

