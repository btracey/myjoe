#ifndef RANSTRSANSMODEL_GARET
#define RANSTRSANSMODEL_GARET

#include "TurbModel_KOMSST.h"


/*
 * transition model based on an intermittency transport equation
 * and a transport equation for the momentum-thickness Reynolds number
 *
 * Ref: "Correlation-Based Transition Modeling for Unstructured Parallelized Computational Fluid Dynamics Codes"
 * AIAA JOURNAL, Vol. 47, No. 12, December 2009,
 * Robin B. Langtry, Florian R. Menter
 *
 * this class inherits from the SST model as the transition model has been developed with it,
 * even though other combinations are possible (define TURB_MOD_FOR_TRANS, see below)
 *
 */



typedef RansTurbKOmSST TURB_MOD_FOR_TRANS; // define the turbulence model used with the transition model



class RansTransGaReT : virtual public TURB_MOD_FOR_TRANS
{
public:   // constructors

  RansTransGaReT()
  {
    if (mpi_rank == 0)
      cout << "RansTransGaReT()" << endl;

    turbModel = KOMSST;

    ce1 = 1.0;
    ce2 = 50.0;
    ctht = 0.03;
    ca1 = 2.0;
    ca2 = 0.06;

    ScalarTranspEq *eq;
    eq = registerScalarTransport("gam", CV_DATA);
    eq->phiZero = getDoubleParam("ZERO_gam", "1.0e-10");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_gam", "1.0e-4");
    eq->phiMaxiter = 20;
    eq->lowerBound = 1.0e-8;
    eq->upperBound = 1.0;

    eq = registerScalarTransport("ReT", CV_DATA);
    eq->phiZero = getDoubleParam("ZERO_ReT", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_ReT", "1.0e-2");
    eq->phiMaxiter = 20;
    eq->lowerBound = 1.0e-8;
    eq->upperBound = 1.0e8;

    Ftht  = NULL;       registerScalar(Ftht, "Ftht", CV_DATA);
    ReThC = NULL;       registerScalar(ReThC, "ReThC", CV_DATA);
    Mais = NULL;       registerScalar(Mais, "MA_IS", CV_DATA);

  }

  virtual ~RansTransGaReT() {}

public:   // member vars

  double *gam, *ReT;                                  ///< turbulent scalars, introduced to have access to variables, results in to more readable code
  double *gam_bfa, *ReT_bfa;                         ///< turbulent scalars at the boundary
  double *Ftht, *ReThC;
  double *Mais;

  double  ce1, ce2, ctht, ca1, ca2;

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    TURB_MOD_FOR_TRANS::initialHookScalarRansTurbModel();

    if (mpi_rank == 0) 
      cout << "initialHook gamma-ReTheta model" << endl;

    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("gam");   gam = eq->phi;    gam_bfa = eq->phi_bfa;
    eq = getScalarTransportData("ReT");   ReT = eq->phi;    ReT_bfa = eq->phi_bfa;
    
    if (!checkScalarFlag("gam"))
      for (int icv=0; icv<ncv; icv++)
        gam[icv] = 1.0;

    if (!checkScalarFlag("ReT"))
      for (int icv=0; icv<ncv; icv++)
        ReT[icv] = 200.0;

    updateCvDataByName("gam", REPLACE_DATA);
    updateCvDataByName("ReT", REPLACE_DATA);
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    TURB_MOD_FOR_TRANS::diffusivityHookScalarRansTurb(name);

    ScalarTranspEq *eq;

    if (name == "gam")
    {
      eq = getScalarTransportData(name);
      for (int ifa=0; ifa<nfa; ifa++)
        eq->diff[ifa] = mul_fa[ifa] + mut_fa[ifa];
    }

    if (name == "ReT")
    {
      eq = getScalarTransportData(name);
      for (int ifa=0; ifa<nfa; ifa++)
        eq->diff[ifa] = 2.0*(mul_fa[ifa] + mut_fa[ifa]);
    }
  }

  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)
  {
    TURB_MOD_FOR_TRANS::sourceHookScalarRansTurb(rhs, A, name);

    if (name == "gam")
    {
      if (mpi_rank == 0 )
        cerr << " explicit gam source not coded yet!" << endl;
      throw(-1);
    }

    if (name == "ReT")
    {
      if (mpi_rank == 0 )
        cerr << " explicit gam source not coded yet!" << endl;
      throw(-1);
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

      crossDiff[icv] = max(2.0*rho[icv]*sigma_om2/omega[icv]*vecDotVec3d(grad_kine[icv], grad_omega[icv]),  1.0e-10);
//      crossDiff[icv] = min(crossDiff[icv], 10000.0);

      double gamma1 = 500.0*mue/(pow(d, 2.0)*rho[icv]*omega[icv]);
      double gamma2 = 4.0*sigma_om2*rho[icv]*kine[icv]/(d*d*crossDiff[icv]);
      double gamma3 = sqrt(kine[icv])/(betaStar*omega[icv]*d);

      double gamma0 = min(max(gamma1, gamma3), gamma2);
      blendFuncF1[icv] = tanh(pow(gamma0,4.0));

      double Rey = rho[icv]*d*sqrt(kine[icv])/mue;
      double F3 = exp(-pow(Rey/120.0, 8.0));
      blendFuncF1[icv] = max(blendFuncF1[icv], F3);

      gamma0 = max(2.0*gamma3, gamma1);
      blendFuncF2[icv] = tanh(pow(gamma0,2.0));

      double U      = max(sqrt(vecDotVec3d(vel[icv], vel[icv])), 1.0e-6);
      double ThBL   = ReT[icv]*mue/(rho[icv]*U);
      double DeBL   = 7.5*ThBL;
      double Delta  = 50.0*vortMag[icv]*d*DeBL/U;
      double ReOm   = rho[icv]*omega[icv]*d*d/mue;
      double FWake  = exp(-pow(ReOm/1.0e+5, 2.0));
      double zeta1  = FWake*exp(-pow(d/Delta, 4.0));
      double zeta2  = 1.0-pow((gam[icv]-1./ce2)/(1.-1./ce2), 2.0);
      Ftht[icv]     = min(max(zeta1, zeta2), 1.0);

      if (ReT[icv]<=1870.0)
        ReThC[icv] = ReT[icv]-( 3.96035-1.20656e-2*ReT[icv]+8.68230e-4*pow(ReT[icv],2.0)
                               -6.96506e-7*pow(ReT[icv], 3.0)+1.74105e-10*pow(ReT[icv], 4.0));
      else
        ReThC[icv] = ReT[icv]-(593.11+ (ReT[icv]-1870.0)*0.482);

      Mais[icv] = sqrt(2.0/(gamma[icv]-1.0)*(-1.0+pow(max(p_ref/press[icv], 1.0), (gamma[icv]-1.0)/gamma[icv])));
    }

    updateCvData(crossDiff, REPLACE_DATA);
    updateCvData(blendFuncF1, REPLACE_DATA);
    updateCvData(blendFuncF2, REPLACE_DATA);
    updateCvData(Ftht, REPLACE_DATA);
    updateCvData(ReThC, REPLACE_DATA);
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double muLam   = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double wallD2  = wallDist[icv]*wallDist[icv];
        double ReNu    = rho[icv]*wallD2*strMag[icv]/muLam;
        double ReTurb  = rho[icv]*kine[icv]/(omega[icv]*muLam);
        double Freatt  = exp(-pow(ReTurb/20.0, 4.0));

        double gam_sep = min(2.0*max(ReNu/(3.235*ReThC[icv])-1., 0.0)*Freatt, 2.0)*Ftht[icv];
        double gam_eff = max(gam[icv], gam_sep);
        double faktDestr = min(max(gam_eff, 0.1), 1.0);

        double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0e5);
        double Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
        if (SST_limitPK == 1)
          Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);

        Pk = max(Pk, 0.0);  // for stabilization
        double src = gam_eff*Pk - faktDestr*betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - faktDestr*betaStar*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }

    if (name == "omega")
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

    if (name == "gam")
    {
      double fakt_nSigma = getDoubleParam("FAKT_NSIGMA", "0.0");

      for (int icv=0; icv<ncv; icv++)
      {
        double muLam  = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double wallD2 = wallDist[icv]*wallDist[icv];
        double ReNu   = rho[icv]*wallD2*strMag[icv]/muLam;
        double ReTurb = rho[icv]*kine[icv]/(omega[icv]*muLam);
        double Fturb  = exp(-pow(ReTurb/4.0, 4.0));

        double Fons1  = ReNu/(2.193*ReThC[icv]);
        double Fons2  = min(max(Fons1, pow(Fons1, 4.0)), 2.0);
        double Fons3  = max(1.0-pow(ReTurb/2.5, 3.0), 0.0);
        double Fons   = max(Fons2-Fons3, 0.0);

        double Flen = 0.3188;
        if (ReT[icv]<400.0)       Flen = 39.8189 - 1.19270e-2*ReT[icv] - 1.32567e-4*pow(ReT[icv], 2.0);
        else if (ReT[icv]<596.0)  Flen = 263.404 - 1.23939*ReT[icv] + 1.94548e-3*pow(ReT[icv], 2.0) - 1.01685e-6*pow(ReT[icv], 3.0);
        else if (ReT[icv]<1200.0) Flen = 0.5 - 3.0e-4*(ReT[icv] - 596.0);

        double ReOm  = rho[icv]*wallD2*omega[icv]/(500.0*muLam);
        double FsubL = exp(-pow(ReOm/0.4, 2.0));
        Flen = Flen*(1.0-FsubL) + 40.0*FsubL;

        double nsigma_mach = max(1.0 - (0.14 + 0.7*fakt_nSigma)*Mais[icv], 0.001);
        double PgConst = /*sqrt*/(nsigma_mach)*Flen*ca1*strMag[icv]*sqrt(Fons);
        double Pg = PgConst*rho[icv]*sqrt(gam[icv])*(1.0-ce1*gam[icv]);

        double EgConst = ca2*vortMag[icv]*Fturb;
        double Eg = EgConst*rho[icv]*gam[icv]*(ce2*gam[icv]-1.0);

        double src = Pg - Eg;
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = PgConst*(/*0.5/sqrt(gam[icv])*/-1.5*ce1*sqrt(gam[icv]))-EgConst*(2.*ce2*gam[icv]-1.);
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }


    if (name == "ReT")
    {
      int myCountNewtonFailure = 0;

//      FILE *fp = fopen("dummy","wt");

      double fakt_Retht = getDoubleParam("FAKT_RETHT", "0.0");


      for (int icv=0; icv<ncv; icv++)
      {
        double U      = max(sqrt(vecDotVec3d(vel[icv], vel[icv])), 1.0e-6);
        double muLam  = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double t      = 500.0*muLam/(rho[icv]*U*U);

        double dUdx = vel[icv][0]*grad_u[icv][0][0] + vel[icv][1]*grad_u[icv][1][0] + vel[icv][2]*grad_u[icv][2][0];
        double dUdy = vel[icv][0]*grad_u[icv][0][1] + vel[icv][1]*grad_u[icv][1][1] + vel[icv][2]*grad_u[icv][2][1];
        double dUdz = vel[icv][0]*grad_u[icv][0][2] + vel[icv][1]*grad_u[icv][1][2] + vel[icv][2]*grad_u[icv][2][2];
        double dUds = 1.0/(U*U)*(vel[icv][0]*dUdx + vel[icv][1]*dUdy + vel[icv][2]*dUdz);

        double TuLoc = max(100.0*sqrt(2.0/3.0*kine[icv])/U, 0.027);
        double Retht0;
        if (TuLoc <= 1.3)   Retht0 = 1173.51 - 589.428*TuLoc + 0.2196/(TuLoc*TuLoc);
        else                Retht0 = 331.50*pow(TuLoc - 0.5658, -0.671);

        double lam0 = rho[icv]*dUds/muLam; // later multiplied with theta**2
        double rhoUmu = rho[icv]*U/muLam;  // Re/Theta

        double theta = Retht0/rhoUmu;   // get initial value of theta for ZPG

//        fprintf(fp, "%.6le\t%.6le\t", x_cv[icv][0], theta);

        double NewtonIter = 0, maxNewtonIter = 10;  // start Newton interation, count the iterations
        do{
          double numeps = 1.0e-12;                  // define small number for numerical gradient estimation
          double fx     = fcorr(theta,        Retht0, lam0, TuLoc, rhoUmu);
          double fxeps  = fcorr(theta-numeps, Retht0, lam0, TuLoc, rhoUmu);

          theta = theta - fx/((fx-fxeps)/numeps);   // Newton step

          NewtonIter++;
        }while((fabs(fcorr(theta, Retht0, lam0, TuLoc, rhoUmu)) > 1.0e-10) && (NewtonIter < maxNewtonIter)); // check

        if (NewtonIter >= maxNewtonIter) myCountNewtonFailure++;

        double Retht = theta*rhoUmu;
//        Retht *= fakt_Retht*sqrt(1.0+0.38*pow(Mais[icv], 0.6)); // it was 1.6 max
        Retht *= (1.0 + 0.6*fakt_Retht*Mais[icv]);

//        fprintf(fp, "%.6le\t%.6le\t%.6le\t%.6le\t%.6le\n", theta, lam0*theta*theta, TuLoc, (double)NewtonIter, Retht);

        double Ptht = ctht*rho[icv]/t*(Retht-ReT[icv])*(1.0-Ftht[icv]);
        rhs[icv] += Ptht*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = -ctht/t*(1.0-Ftht[icv]);
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }

      int countNewtonFailure;
      MPI_Allreduce(&myCountNewtonFailure, &countNewtonFailure, 1, MPI_INT, MPI_SUM, mpi_comm);
      if (mpi_rank == 0)
        if (countNewtonFailure > 0)
          cout << "Newton iter failed in Ret eq: " << countNewtonFailure << endl;


      /*FILE *fp = fopen("lam.dat","wt");
      fprintf(fp, "\t");
      for (int iTu = 1; iTu <= 10; ++iTu)
      {
        double Tu =  (double)iTu/10.0*3.0;
        fprintf(fp, "%.6le\t", Tu);
      }
      fprintf(fp, "\n");

      for (int ilam = 0; ilam <= 50; ++ilam)
      {
        double lam = -0.1 + (double)ilam/50.0*0.2;

        fprintf(fp, "%.6le\t", lam);

        for (int iTu = 1; iTu <= 10; ++iTu)
        {
          double Tu =  (double)iTu/10.0*3.0;
          double flam = fcorr(lam, Tu);
          fprintf(fp, "%.6le\t", flam);
        }
        fprintf(fp, "\n");
      }

      fclose(fp);
      throw(-1);*/
    }
  }
/*
  inline double fcorr(const double lam, const double Tu)
  {
    double flam;
    if (lam <= 0.0) flam = 1. - (-12.986*lam-123.66*lam*lam-405.689*pow(lam,3.0))*exp(-pow(Tu/1.5, 1.5));
    else             flam = 1. + 0.275*(1.0-exp(-35.0*lam))*exp(-2.0*Tu);

    return flam;
  }*/

  inline double fcorr(double theta, const double Retht0, const double lam0, const double Tu, const double rhoUmu)
  {
    double lam = max(min(lam0*theta*theta, 0.1), -0.1);  // another bound
    double flam;
    if (lam0 <= 0.0) flam = 1. - (-12.986*lam-123.66*lam*lam-405.689*pow(lam,3.0))*exp(-pow(Tu/1.5, 1.5));
    else             flam = 1. + 0.275*(1.0-exp(-35.0*lam))*exp(-2.0*Tu);

    return (rhoUmu*theta - Retht0*flam);
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    TURB_MOD_FOR_TRANS::boundaryHookScalarRansTurb(phi_fa, zone, name);

    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if (name == "ReT")
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            double U     = max(sqrt(vecDotVec3d(vel_bfa[ifa], vel_bfa[ifa])), 1.0e-6);
            double TuLoc = max(100.0*sqrt(2.0*kine_bfa[ifa]/3.0)/U, 0.027);

            double Retht;
            if (TuLoc <= 1.3)   Retht = 1173.51 - 589.428*TuLoc + 0.2196/(TuLoc*TuLoc);
            else                Retht = 331.50*pow(TuLoc - 0.5658, -0.671);

            phi_fa[ifa] = Retht;
          }
        }
    }
  }
};



#endif


