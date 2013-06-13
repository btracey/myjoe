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

    C_MU  = getDoubleParam("C_MU",  "0.09");
    SIG_K = getDoubleParam("SIG_K", "1.0");
    SIG_D = getDoubleParam("SIG_D", "1.3");
    CEPS1 = getDoubleParam("CEPS1", "1.35");
    CEPS2 = getDoubleParam("CEPS2", "1.80");
    CETA  = getDoubleParam("CETA",  "70.0");
    CL    = getDoubleParam("CL",    "0.23");
    
    realizable = getIntParam("REALIZABILITY", "1");
    KEps_limitPK = getIntParam("KEps_LIMIT_PK", "1");

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
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "1.0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-3;
    eq->phiMaxiter = 100;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_D;

    strMag   = NULL;     registerScalar(strMag,   "strMag",   CV_DATA);
    diverg   = NULL;     registerScalar(diverg,   "diverg",   CV_DATA);
    turbTS   = NULL;     registerScalar(turbTS,   "turbTS",   CV_DATA);
    turbLS   = NULL;     registerScalar(turbLS,   "turbLS",   CV_DATA);
    muT      = NULL;     registerScalar(muT,      "muT",      CV_DATA);
    fmu      = NULL;     registerScalar(fmu,      "fmu",      CV_DATA);
    f2       = NULL;     registerScalar(f2,       "f2",       CV_DATA);

    yplus    = NULL;     registerScalar(yplus,    "yplus",    CV_DATA);
    wallDist = NULL;     registerScalar(wallDist, "wallDist", CV_DATA);
    wallConn = NULL;     // array of integers

    invDeltaNu    = NULL;
    my_invDeltaNu = NULL;

    wall_send_count = NULL;     wall_send_displ = NULL;
    wall_recv_count = NULL;     wall_recv_displ = NULL;

    debug1 = NULL; registerScalar(debug1, "debug1", CV_DATA);
    debug2 = NULL; registerScalar(debug2, "debug2", CV_DATA);
    debug3 = NULL; registerScalar(debug3, "debug3", CV_DATA);
    debug4 = NULL; registerScalar(debug4, "debug4", CV_DATA);
  }

  virtual ~RansTurbKEps() {}

public:   // member vars

  double *eps;                      ///< dissipation, introduced to have access to variables, results in to more readable code
  double *kine_bfa, *eps_bfa;       ///< turbulent scalars at the boundary
  double *muT;                      ///< turbulent viscosity at cell center for output
  double *fmu, *f2;                 ///< Chien damping functions
  double *yplus;                    ///< distance to closest wall in wall units
  double *wallDist;                 ///< distance to closest wall face
  int *wallConn;                    ///< index of closest wall face
  
  double *invDeltaNu;               ///< global inverse viscous length scale
  double *my_invDeltaNu;            ///< local inverse viscous length scale

  int *wall_send_count, *wall_send_displ;
  int *wall_recv_count, *wall_recv_displ;

  int realizable;     ///< upper limiter for time scale: 0... no
                      ///                                1... yes
  int KEps_limitPK;   ///< limiter for tke production: 0 ... no
                      ///                              1 ... yes

  double C_MU, SIG_K, SIG_D, CEPS1, CEPS2, CETA, CL;

  double *debug1, *debug2, *debug3, *debug4;
public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook KEPS model" << endl;

    // connect pointers 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");  kine = eq->phi;   kine_bfa = eq->phi_bfa;
    eq = getScalarTransportData("eps");   eps = eq->phi;    eps_bfa = eq->phi_bfa;

    // wall distance and connectivity
    wallConn = new int[ncv];
    calcWallDistance(wallConn, wallDist);

    // wall arrays
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
          {
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
              my_wall_faces++;
          }
      }

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    wall_send_count = new int[mpi_size];
    wall_send_displ = new int[mpi_size];
    for (int r=0; r<mpi_size; r++)
    {
      wall_send_count[r] = my_wall_faces;
      wall_send_displ[r] = 0;
    }

    wall_recv_count = new int[mpi_size];
    wall_recv_displ = new int[mpi_size];
    MPI_Allgather(&my_wall_faces, 1, MPI_INT, wall_recv_count, 1, MPI_INT, mpi_comm);

    wall_recv_displ[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      wall_recv_displ[i] = wall_recv_count[i-1] + wall_recv_displ[i-1];

    invDeltaNu    = new double[tot_wall_faces];
    my_invDeltaNu = new double[my_wall_faces];
  }

  virtual void calcRansTurbViscMuet()
  {
    calcGradVel();
    calcStrainRateAndDivergence();
    calcYplus();
    calcChienDampingFunctions();

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

      double rho_fa  = (w1*rho[icv0] + w0*rho[icv1])/(w0 + w1);
      double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0 + w1);
      double eps_fa  = (w1*eps[icv0] + w0*eps[icv1])/(w0 + w1);
      double fmu_fa  = (w1*fmu[icv0] + w0*fmu[icv1])/(w0 + w1);

      mut_fa[ifa] = min(max(C_MU*fmu_fa*rho_fa*kine_fa*kine_fa/eps_fa,0.0), 1.0);
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
            mut_fa[ifa] = min(max(C_MU*fmu[icv0]*rho[icv0]*kine[icv0]*kine[icv0]/eps[icv0], 0.0), 1.0);
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
      for (int ifa=nfa_b; ifa<nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }
  }
  
  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double muLamCV = calcMuLam(icv);
        double src = calcTurbProd(icv) - rho[icv]*eps[icv] - 2.0*muLamCV*kine[icv]/(wallDist[icv]*wallDist[icv]);
        rhs[icv] += src*cv_volume[icv];

        //d(-rho*kine*eps/kine)/d(rho*kine)
        if (flagImplicit)
        {
          double dsrcdphi = -eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }

    if (name == "eps")
      for (int icv=0; icv<ncv; icv++)
      {
        double muLamCv = calcMuLam(icv);
        /*double TS = calcTurbTimeScale(kine[icv], eps[icv], strMag[icv], muLamCv/rho[icv], realizable);
        double src = (CEPS1*calcTurbProd(icv) - CEPS2*f2[icv]*rho[icv]*eps[icv])/TS
                   - 2.0*muLamCv*eps[icv]*exp(-0.5*yplus[icv])/(wallDist[icv]*wallDist[icv]);*/

        double src = CEPS1*calcTurbProd(icv)*eps[icv]/kine[icv] - CEPS2*f2[icv]*rho[icv]*eps[icv]*eps[icv]/kine[icv]
                   - 2.0*muLamCv*eps[icv]*exp(-0.5*yplus[icv])/(wallDist[icv]*wallDist[icv]);

        rhs[icv] += src*cv_volume[icv];

        debug1[icv] = src;
        debug2[icv] = CEPS1*calcTurbProd(icv)*eps[icv]/kine[icv];
        debug3[icv] = -CEPS2*f2[icv]*rho[icv]*eps[icv]*eps[icv]/kine[icv];
        debug4[icv] = 2.0*muLamCv*eps[icv]*exp(-0.5*yplus[icv])/(wallDist[icv]*wallDist[icv]);
        // d(-ce2*rho*eps/TS)/d(rho*eps)
        if (flagImplicit)
        {
          double dsrcdphi = -CEPS2*eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {/*
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
    }*/
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
    double mu_t = min(max(C_MU*fmu[icv]*rho[icv]*kine[icv]*kine[icv]/eps[icv], 0.0), 1.0);

    double Pk = max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
    if (KEps_limitPK == 1)
      Pk = min(Pk, 20.0*rho[icv]*eps[icv]);

    return Pk;
  }

  virtual void calcYplus()
  {
    // compute my_deltaNu
    int count = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            {
              int icv = cvofa[ifa][0];

              double n[3], s_half[3], vel[3], velTang[3];

              normVec3d(n, fa_normal[ifa]);
              vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);

              vel[0] = rhou[icv][0]/rho[icv];
              vel[1] = rhou[icv][1]/rho[icv];
              vel[2] = rhou[icv][2]/rho[icv];
              double un = vecDotVec3d(n, vel);

              velTang[0] = vel[0] - un*n[0];
              velTang[1] = vel[1] - un*n[1];
              velTang[2] = vel[2] - un*n[2];

              double velMag = sqrt(vecDotVec3d(velTang, velTang));
              double walld = fabs(vecDotVec3d(s_half, n));

              double mulam = mul_fa[ifa];

              double tau = mulam*velMag/walld*sign(vel[0]);
              double utau = sqrt(fabs(tau)/rho[icv]);

              my_invDeltaNu[count] = utau/(mulam/rho[icv]);
              count++;
            }
      }

    // assemble deltaNu
    MPI_Alltoallv(my_invDeltaNu, wall_send_count, wall_send_displ, MPI_DOUBLE,
                  invDeltaNu,    wall_recv_count, wall_recv_displ, MPI_DOUBLE, mpi_comm);

    // compute yplus
    for (int icv = 0; icv < ncv; icv++)
      yplus[icv] = wallDist[icv]*invDeltaNu[wallConn[icv]];

    updateCvData(yplus, REPLACE_DATA);
  }

  virtual void calcChienDampingFunctions()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      fmu[icv] = 1.0 - exp(-0.0115*yplus[icv]);

      double muLamCV = calcMuLam(icv);
      double Ret = rho[icv]*kine[icv]*kine[icv]/(muLamCV*eps[icv]);
      f2[icv] = 1.0 - 0.4/1.8*exp(-Ret*Ret/36.0);
    }

    updateCvData(fmu, REPLACE_DATA);
    updateCvData(f2, REPLACE_DATA);
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (wallConn != NULL) {delete [] wallConn; wallConn = NULL;}

    if (invDeltaNu    != NULL) {delete [] invDeltaNu;    invDeltaNu    = NULL;}
    if (my_invDeltaNu != NULL) {delete [] my_invDeltaNu; my_invDeltaNu = NULL;}

    if (wall_send_count != NULL) {delete [] wall_send_count; wall_send_count = NULL;}
    if (wall_send_displ != NULL) {delete [] wall_send_displ; wall_send_displ = NULL;}
    if (wall_recv_count != NULL) {delete [] wall_recv_count; wall_recv_count = NULL;}
    if (wall_recv_displ != NULL) {delete [] wall_recv_displ; wall_recv_displ = NULL;}

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








