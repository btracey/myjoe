#ifndef RANSTURBMODEL_SA_H
#define RANSTURBMODEL_SA_H

#include "UgpWithCvCompFlow.h"




class RansTurbSA : virtual public UgpWithCvCompFlow
{
public:   // constructors
  RansTurbSA()
  {
    if (mpi_rank == 0)
      cout << "RansTurbSA()" << endl;

    turbModel = SA;


    // model constants
    cv1_3         = getDoubleParam("cv1_3",   "357.911");
    cw1           = getDoubleParam("cw1",     "3.239067817");
    cw2           = getDoubleParam("cw2",     "0.3");
    cw3_6         = getDoubleParam("cw3_6",   "64.0");
    sigma         = getDoubleParam("sigma",   "0.6666666");
    cb1           = getDoubleParam("cb1",     "0.1355");
    cb2           = getDoubleParam("cb2",     "0.622");
    KarmanConst_2 = getDoubleParam("KarmanConst_2", "0.1681");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("nuSA", CV_DATA);
    eq->relax = getDoubleParam("RELAX_SA", "0.5");
    eq->phiZero = 1.0e-9;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0;
    eq->turbSchmidtNumber = 1.0/sigma;
//    eq->reconstruction = "CONSERVATIVE";

    strMag   = NULL;      registerScalar(strMag, "strMag", CV_DATA);
    vortMag  = NULL;      registerScalar(vortMag, "vortMag", CV_DATA);
    muT      = NULL;      registerScalar(muT, "muT", CV_DATA);
    wallDist = NULL;      registerScalar(wallDist, "wallDist",  CV_DATA);
  }

  virtual ~RansTurbSA() {}

public:   // member vars

  double *nuSA;                             ///< turbulent scalars, introduced to have access to variables, results in to more readable code
  double *nuSA_bfa;                         ///< turbulent scalars at the boundary
  double *muT;                              ///< turbulent viscosity at cell center for output
  double *wallDist;                         ///< wall distance
  
  double cv1_3, cw1, cw2, cw3_6, sigma, cb1, cb2, KarmanConst_2;

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook SA model" << endl;
    
    if (!checkParam("DO_NOT_CALC_WALLDIST"))
    {
      if (!checkDataFlag(wallDist))
      {
        for (int icv=0; icv<ncv; icv++)
          wallDist[icv] = 0.0;
      
        calcWallDistance(NULL, wallDist);
      }    
    }
    
    // connect pointers 
    ScalarTranspEq *eq = getScalarTransportData("nuSA");      nuSA = eq->phi;       nuSA_bfa = eq->phi_bfa;

    double init_nuSA = getDoubleParam("INIT_NU_SA", "1.0e-6");
    
    if (!checkScalarFlag("nuSA"))
      for (int icv=0; icv<ncv; icv++)
        nuSA[icv] = init_nuSA; // arbitrary value  

    updateCvDataByName("nuSA", REPLACE_DATA);
  }

  virtual void calcRansTurbViscMuet()
  {
    // update grad_u with current state vector
    calcGradVel();
    // update strain rate tensor 
    calcStrainRate();
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

      double rho_fa  = (w1*rho[icv0]   + w0*rho[icv1])/(w0+w1);
      double nuSA_fa = (w1*nuSA[icv0]  + w0*nuSA[icv1])/(w0+w1);

      double muSA_fa = rho_fa*nuSA_fa;
      double chi_3   = pow(muSA_fa/mul_fa[ifa], 3.0);
      double fv1     = chi_3/(chi_3 + cv1_3);
      
      mut_fa[ifa]  = max(muSA_fa*fv1, 0.0);              
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

          double muSA  = rho[icv0]*nuSA[icv0];
          double mul_icv0 = InterpolateAtCellCenterFromFaceValues(mul_fa, icv0);
          double chi_3 = pow(muSA/mul_icv0, 3.0);
          double fv1   = chi_3/(chi_3 + cv1_3);
          mut_fa[ifa] = max(muSA*fv1, 0.0);    // zero order extrapolation for others
        }
    }

    // compute mut at cell center for output only 
    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);

  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    if (name == "nuSA")
    {
      ScalarTranspEq *eq = getScalarTransportData(name);
      double invSigma = 1.0/sigma;
      
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

        double rho_fa  = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
        double nuSA_fa = (w1*nuSA[icv0] + w0*nuSA[icv1])/(w0+w1);

        eq->diff[ifa] = invSigma*(mul_fa[ifa] + rho_fa*nuSA_fa);
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = invSigma*mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = invSigma*(mul_fa[ifa] + rho_bfa[ifa]*nuSA_bfa[ifa]);
          }
      }
    }
  }

  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)
  {   
    if (name == "nuSA")
    {
      // calculate gradient of nuSA
      double (*grad_nuSA)[3] = new double[ncv_g][3];
      calcCvScalarGrad(grad_nuSA, nuSA, nuSA_bfa, gradreconstruction, limiterNavierS, nuSA, epsilonSDWLS);
      
      for (int icv=0; icv<ncv; icv++)
      {
        // model variables
        double mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double Om         = max(vortMag[icv], 1.e-10);   // vorticity
        double Strainrate = max(strMag[icv], 1.e-10);    // strainrate
        double d = max(wallDist[icv], 1.0e-10);

        // model functions
        double chi   = nuSA[icv]*rho[icv]/mue;
        double chi_3 = pow(chi, 3.0);

        double inv_KarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        double fv1 = chi_3/(chi_3 + cv1_3);
        double fv2 = 1.0 - chi/(1.0 + chi*fv1);

        double BlendTurbulentProduction = 0.0;
        double S = max(Om + min(0.0, Strainrate - Om)*BlendTurbulentProduction, 0.0);
        double Shat = S + nuSA[icv]*fv2*inv_KarmanConst2_d2;

        double inv_Shat;
        if(fabs(Shat) < 1.0e-10) 
          inv_Shat = 1.0/(1.0e-10*sign(Shat));
        else 
          inv_Shat = 1.0/Shat;

        double r = min(nuSA[icv]*inv_Shat*inv_KarmanConst2_d2, 10.0);
        //double r = nuSA[icv]*inv_Shat*inv_KarmanConst2_d2;
        double g = r + cw2*(pow(r, 6.0)-r);
        double fw = g*pow((1.0+cw3_6)/(pow(g, 6.0)+cw3_6), 1.0/6.0);

        // model sources
        double H1 = rho[icv]*cb1*Shat*nuSA[icv];
        double H3 = rho[icv]*cb2/sigma*vecDotVec3d(grad_nuSA[icv], grad_nuSA[icv]);
        
        rhs[icv] += cv_volume[icv]*(H1 + H3);

        // destruction term divided by rho*nuSA
        double H2 = rho[icv]*cw1*fw*nuSA[icv]/pow(d, 2.0);
        int noc00 = nbocv_i[icv];
        A[noc00] += H2*cv_volume[icv];
      }

      delete [] grad_nuSA;
    }
  }

  virtual void sourceHookScalarRansTurbExplicit(double *rhs, const string &name)
  {   
    if (name == "nuSA")
    {
      // calculate gradient of nuSA
      double (*grad_nuSA)[3] = new double[ncv_g][3];
      calcCvScalarGrad(grad_nuSA, nuSA, NULL, gradreconstruction, limiterNavierS, nuSA, epsilonSDWLS);
      
      for (int icv=0; icv<ncv; icv++)
      {
        // model variables
        double mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double Om         = max(vortMag[icv], 1.e-8);   // vorticity
        double Strainrate = max(strMag[icv], 1.e-8);    // strainrate
        double d = wallDist[icv];

        // model functions
        double chi   = nuSA[icv]*rho[icv]/mue;
        double chi_3 = pow(chi, 3.0);

        double inv_rhoKarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        double fv1 = chi_3/(chi_3+cv1_3);
        double fv2 = 1.0 - chi/(1.0+chi*fv1);

        double BlendTurbulentProduction = 0.0;
        double S = max(Om + min(0.0, Strainrate- Om)*BlendTurbulentProduction, 0.0);
        double Shat = S + nuSA[icv]*fv2*inv_rhoKarmanConst2_d2;

        double inv_Shat;
        if(fabs(Shat) < 1.0e-10) 
          inv_Shat = 1.0/(1.0e-10*sign(Shat));
        else 
          inv_Shat = 1.0/Shat;

        double r = min(nuSA[icv]*inv_Shat*inv_rhoKarmanConst2_d2, 10.0);
        double g = r + cw2*(pow(r, 6.0)-r);
        double fw = g*pow((1.0+cw3_6)/(pow(g, 6.0)+cw3_6), 1.0/6.0);

        // model sources
        double H1 = rho[icv]*cb1*Shat*nuSA[icv];
        double H2 = rho[icv]*cw1*fw*nuSA[icv]*nuSA[icv]/pow(d, 2.0);
        double H3 = rho[icv]*cb2/sigma*vecDotVec3d(grad_nuSA[icv], grad_nuSA[icv]);

        rhs[icv] += cv_volume[icv]*(H1 - H2 + H3);
      }
      
      delete [] grad_nuSA;
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "nuSA")
    {
      // get gradient of nuSA
      double (*grad_nuSA)[3] = getScalarTransportData("nuSA")->grad_phi;

      for (int icv=0; icv<ncv; icv++)
      {
        // model variables
        double mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double Om         = max(vortMag[icv], 1.e-8);   // vorticity
        double Strainrate = max(strMag[icv], 1.e-8);    // strainrate
        double d = wallDist[icv];

        // model functions
        double chi   = nuSA[icv]*rho[icv]/mue;
        double chi_3 = pow(chi, 3.0);

        double inv_KarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        double fv1 = chi_3/(chi_3+cv1_3);
        double fv2 = 1.0 - chi/(1.0+chi*fv1);

        double BlendTurbulentProduction = 0.0;
        double S = max(Om + min(0.0, Strainrate- Om)*BlendTurbulentProduction, 0.0);
        double Shat = S + nuSA[icv]*fv2*inv_KarmanConst2_d2;

        double inv_Shat = 1.0/max(Shat, 1.0e-8);

/*        if(fabs(Shat) < 1.0e-10)
          inv_Shat = 1.0/(1.0e-10*sign(Shat));
        else
          inv_Shat = 1.0/Shat;*/

        double r   = min(nuSA[icv]*inv_Shat*inv_KarmanConst2_d2, 10.0);
        double g   = r + cw2*(pow(r, 6.0)-r);
        double g_6 = pow(g, 6.0);
        double fw  = g*pow((1.0+cw3_6)/(g_6+cw3_6), 1.0/6.0);

        // model sources
        double H1 =  cb1*Shat*nuSA[icv];
        double H2 = -cw1*fw*nuSA[icv]*nuSA[icv]/pow(d, 2.0);
        double H3 =  cb2/sigma*vecDotVec3d(grad_nuSA[icv], grad_nuSA[icv]);

        rhs[icv] += rho[icv]*cv_volume[icv]*(H1 + H2 + H3);

        if (flagImplicit)
        {
/*          double d_chi = 1./mue;
          double chi_3_p_cv1_3_min1 = 1./(chi_3 + cv1_3);
          double d_fv1 = 3.*chi*chi*chi_3_p_cv1_3_min1*(1. - chi_3*chi_3_p_cv1_3_min1)*d_chi;

          double chi_fv1_p_1_min1 = 1./(1.+chi*fv1);
          double d_fv2 = - d_chi*chi_fv1_p_1_min1 + chi*chi_fv1_p_1_min1*chi_fv1_p_1_min1*(d_chi*fv1 + chi*d_fv1);
          double d_S_s = (fv2/rho[icv] + nuSA[icv]*d_fv2)*inv_KarmanConst2_d2;

          double d_r  = inv_KarmanConst2_d2*(1./rho[icv] - nuSA[icv]*inv_Shat*d_S_s)*inv_Shat;
          double d_g  = d_r*( 1. + cw2*(6.*pow(r,5.) - 1) );

          double g_6_cw3_6_exp = pow(g_6+ cw3_6, -1.0/6.0);
          double cw3_6_p1_exp  = pow(1.0+cw3_6, 1.0/6.0);
          double d_fw = cw3_6_p1_exp*g_6_cw3_6_exp*d_g*( 1. - g_6/(g_6+cw3_6));

          double d_H1 =  cb1*(Shat + rho[icv]*nuSA[icv]*d_S_s);
          double d_H2 = -cw1*nuSA[icv]/d/d*(d_fw*rho[icv]*nuSA[icv] + 2.0*fw);
*/
          double d_H1 =  0.0;//cb1*Shat;
          double d_H2 = -2.0*cw1*fw*nuSA[icv]/pow(d, 2.0);

          int noc00 = nbocv_i[icv];
          A[noc00] -= (d_H1+d_H2)*cv_volume[icv];
        }
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int nuSA_Index = getScalarTransportIndex("nuSA");
    
    // get gradient of nuSA
    double (*grad_nuSA)[3] = getScalarTransportData("nuSA")->grad_phi;

    for (int icv = 0; icv < ncv; icv++)
    {
      // model variables
      double mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
      double Om         = max(vortMag[icv], 1.e-8);   // vorticity
      double Strainrate = max(strMag[icv], 1.e-8);    // strainrate
      double d = wallDist[icv];

      // model functions
      double chi   = nuSA[icv] * rho[icv] / mue;
      double chi_3 = pow(chi, 3.0);

      double inv_KarmanConst2_d2 = 1.0 / (KarmanConst_2 * d * d);

      double fv1 = chi_3 / (chi_3 + cv1_3);
      double fv2 = 1.0 - chi / (1.0 + chi * fv1);

      double BlendTurbulentProduction = 0.0;
      double S = max(Om + min(0.0, Strainrate- Om) * BlendTurbulentProduction, 0.0);
      double Shat = S + nuSA[icv] * fv2 * inv_KarmanConst2_d2;

      double inv_Shat = 1.0 / max(Shat, 1.0e-8);

/*        if(fabs(Shat) < 1.0e-10)
        inv_Shat = 1.0/(1.0e-10*sign(Shat));
      else
        inv_Shat = 1.0/Shat;*/

      double r   = min(nuSA[icv] * inv_Shat * inv_KarmanConst2_d2, 10.0);
      double g   = r + cw2 * (pow(r, 6.0) - r);
      double g_6 = pow(g, 6.0);
      double fw  = g * pow((1.0 + cw3_6) / (g_6 + cw3_6), 1.0 / 6.0);

      // model sources
      double H1 =   cb1 * Shat * nuSA[icv];
      double H2 = - cw1 * fw * nuSA[icv] * nuSA[icv] / pow(d, 2.0);
      double H3 =   cb2 / sigma * vecDotVec3d(grad_nuSA[icv], grad_nuSA[icv]);

      rhs[icv][5+nuSA_Index] += rho[icv] * cv_volume[icv] * (H1 + H2 + H3);

      if (flagImplicit)
      {
/*          double d_chi = 1./mue;
        double chi_3_p_cv1_3_min1 = 1./(chi_3 + cv1_3);
        double d_fv1 = 3.*chi*chi*chi_3_p_cv1_3_min1*(1. - chi_3*chi_3_p_cv1_3_min1)*d_chi;

        double chi_fv1_p_1_min1 = 1./(1.+chi*fv1);
        double d_fv2 = - d_chi*chi_fv1_p_1_min1 + chi*chi_fv1_p_1_min1*chi_fv1_p_1_min1*(d_chi*fv1 + chi*d_fv1);
        double d_S_s = (fv2/rho[icv] + nuSA[icv]*d_fv2)*inv_KarmanConst2_d2;

        double d_r  = inv_KarmanConst2_d2*(1./rho[icv] - nuSA[icv]*inv_Shat*d_S_s)*inv_Shat;
        double d_g  = d_r*( 1. + cw2*(6.*pow(r,5.) - 1) );

        double g_6_cw3_6_exp = pow(g_6+ cw3_6, -1.0/6.0);
        double cw3_6_p1_exp  = pow(1.0+cw3_6, 1.0/6.0);
        double d_fw = cw3_6_p1_exp*g_6_cw3_6_exp*d_g*( 1. - g_6/(g_6+cw3_6));

        double d_H1 =  cb1*(Shat + rho[icv]*nuSA[icv]*d_S_s);
        double d_H2 = -cw1*nuSA[icv]/d/d*(d_fw*rho[icv]*nuSA[icv] + 2.0*fw);
*/
        double d_H1 =   0.0; //cb1*Shat;
        double d_H2 = - 2.0 * cw1 * fw * nuSA[icv] / pow(d, 2.0);

        int noc00 = nbocv_i[icv];
//        A[noc00][5+nuSA_Index][0] -= - H2 * cv_volume[icv];
        A[noc00][5+nuSA_Index][5+nuSA_Index] -= (d_H1 + d_H2) * cv_volume[icv];
      }
    }
  }

};


#endif
