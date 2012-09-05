#include "UgpWithCvCompFlow.h"



void UgpWithCvCompFlow::solveScalars(double *massFlux)
{
  double *rhs = new double[ncv];
  double *A = new double[nbocv_s];

  // solve the scalars

  for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    solveScalarTransport(&(*data), rhs, A, massFlux);

  delete []rhs;
  delete []A;
}

void UgpWithCvCompFlow::solveScalarsImplUnsteady(double *massFlux)
{
  double *rhs = new double[ncv];
  double *A = new double[nbocv_s];

  // solve the scalars

  for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    solveScalarTransportImplUnsteady(&(*data), rhs, A, massFlux);

  delete []rhs;
  delete []A;
}



void UgpWithCvCompFlow::solveScalarsExplicit(double *massFlux)
{
  double *rhs = new double[ncv];

  // solve the scalars

  for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    solveScalarTransportExplicit(&(*data), rhs, massFlux);

  delete []rhs;
}



void UgpWithCvCompFlow::setScalarBC(FaZone *zone)
{

  for (ScalarTranspEqIterator scal = scalarTranspEqVector.begin(); scal < scalarTranspEqVector.end(); scal++)
  {
    string scalName(scal->getName());
    double phiBCval = 0.0, phiBCflux = 0.0;

    double *phi_fa = scal->phi_bfa;

    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      // .............................................................................................
      // HOOK BOUNDARY CONDITION
      // .............................................................................................
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        // user defined boundary hook
        boundaryHookScalarRansTurb(phi_fa, &(*zone), scalName);
        boundaryHookScalarRansComb(phi_fa, &(*zone), scalName);
      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          if (vecDotVec3d(vel_bfa[ifa], fa_normal[ifa]) > 1.0e-8)
            phi_fa[ifa] = scal->phi[icv0];
          else
            phi_fa[ifa] = phiBCval;
        }
      }
      // .............................................................................................
      // FLUX BOUNDARY CONDITION: NOT IMPLEMENTED
      // .............................................................................................
      else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
      {
        cerr << "scalarZoneIsFlux not implemented yet!" << endl;
        throw(-1);
      }
      // .............................................................................................
      // OTHER BOUNDARY CONDITIONS: NO FLUX
      // .............................................................................................
      else // apply zero flux if no boundary condition is specified
      {
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          phi_fa[ifa] = scal->phi[icv0];
        }
      }
    }
  }
}


/**
 *  solveScalarTransport: steady (no time term) implicit
 */
void UgpWithCvCompFlow::solveScalarTransport(ScalarTranspEq *scal, double *rhs, double *A, double *massFlux_fa)
{
  // initialize LHS and RHS and everything else
  for (int noc=0; noc<nbocv_s; noc++)     A[noc] = 0.0;
  for (int icv=0; icv<ncv; icv++)         rhs[icv] = 0.0;

  // -----------------------------------------------------------
  // define values at the boundary first (rho*phi)_fa
  // -----------------------------------------------------------
  string scalName(scal->getName());
  double phiBCval = 0.0, phiBCflux = 0.0;
  
  double *phi_fa = scal->phi_bfa;

  // -----------------------------------------------------------
  // then calculate gradient with boundary conditions
  // -----------------------------------------------------------

  // calc gradient of phi for viscous flux and second order
  static double (*gradPhi)[3] = new double[ncv_g][3];
  calcCvScalarGrad(gradPhi, scal->phi, phi_fa, gradreconstruction, limiterNavierS, scal->phi, epsilonSDWLS);

  // -----------------------------------------------------------
  // calc user defined diffusivity
  // -----------------------------------------------------------
  if (scal->diffTerm)
  {
    diffusivityHookScalarRansTurb(string(scal->getName()));
    diffusivityHookScalarRansComb(string(scal->getName()));
  }

  // ====================================================================
  // Part 1: interior face-based fluxes [nfa_b:nfa] ...
  // ====================================================================
  for (int ifa=nfa_b; ifa<nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];

    // ####################################### convective flux #######################################
    double lamPl = 0.0, lamMi = 0.0, explConv = 0.0;
    if (scal->convTerm)
    {
      lamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
      lamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));

      if ((sndOrder == true) && (stOrderScalar == false))
      {
        double dx0[3],dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

        explConv = sndOrder*(lamPl*vecDotVec3d(gradPhi[icv0], dx0) + lamMi*vecDotVec3d(gradPhi[icv1], dx1));
      }
    }

    // ####################################### viscous flux #######################################
    double explVisc = 0.0, implVisc = 0.0;
    if (scal->diffTerm)
    {
      double n[3], s[3];
      double nmag = normVec3d(n, fa_normal[ifa]);

      vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(s);
      double alpha = vecDotVec3d(n, s);

      double diff = scal->diff[ifa];
      implVisc = -diff*nmag*alpha/smag;

      for (int i=0; i<3; i++)
        explVisc += (gradPhi[icv0][i] + gradPhi[icv1][i])*(n[i] - alpha*s[i]);
      explVisc *= -0.5*diff*nmag;
    }

    //####################################### build equation system #######################################
    int noc00, noc01, noc11, noc10;
    getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    A[noc00]  += lamPl - implVisc;
    A[noc01]  += lamMi + implVisc;
    rhs[icv0] -= explConv + explVisc;

    if (icv1 < ncv)
    {
      A[noc11]  += -lamMi - implVisc;
      A[noc10]  += -lamPl + implVisc;
      rhs[icv1] += explConv + explVisc;
    }
  }



  // ====================================================================
  // Part 2: boundary face-based fluxes
  // ====================================================================

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int noc00 = nbocv_i[icv0];

          if (scal->convTerm)
          {
            if (massFlux_fa[ifa] > 0.0)     A[noc00]  += massFlux_fa[ifa];
            else                            rhs[icv0] -= massFlux_fa[ifa]*phi_fa[ifa];
          }

          if (scal->diffTerm)
          {
            double n[3], s_half[3];
            double nmag = normVec3d(n, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
            double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);
            
            double diff = scal->diff[ifa];

            double impl = -diff*nmag/smag_half;
            double expl = 0.0;  // no need for a tangential correction

            A[noc00]  -=  impl;
            rhs[icv0] +=  expl - impl*phi_fa[ifa];
          }
        }
      }
      else {
        if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int noc00 = nbocv_i[icv0];

            double phiBCvalTMP = phiBCval;

            if (scal->convTerm)
            {
              if (massFlux_fa[ifa] > 1.0e-8)
              {
                phiBCvalTMP = scal->phi[icv0];
                A[noc00] += massFlux_fa[ifa];
              }
              else
              {
                rhs[icv0] -= massFlux_fa[ifa]*phiBCvalTMP;
              }
            }

            if (scal->diffTerm)
            {
              double n[3], s_half[3];
              double nmag = normVec3d(n, fa_normal[ifa]);
              vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
              double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);

              double diff = scal->diff[ifa];

              double coeff = diff*nmag/smag_half;
              double impl = -coeff;
              double expl = 0.0;  // no need for a tangential correction

              A[noc00]  -=  impl;
              rhs[icv0] +=  expl + coeff*phiBCvalTMP;
            }
          }
        }
        else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
        {
          cerr << "scalarZoneIsFlux not implemented yet!" << endl;
          throw(-1);
        }
        else // apply zero flux if no boundary condition is specified
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is symmetry / zero flux" << endl;
        }
      }
    }
  }


  // ====================================================================
  // Part 3: add source term .......
  // ====================================================================

  sourceHookScalarRansTurb(rhs, A, scal->getName());
  sourceHookScalarRansComb(rhs, A, scal->getName());

  calcResidual(scal->resAve, scal->resMax, scal->phi, A, rhs);


  // correction for steady state
  // d(rho phi)/dt + d(rho phi u)/dx = ....
  // rho d(phi)/dt + phi d(rho)/dt = +d(rho phi u)/dt = ....
  //
  //      ^                  ^
  //      |                  |
  //                         |
  //  zero (steady state)    |
  //                         |
  //
  //                  d(rho)/dt = -d(rho u_i)/dx_i
  //
  for (int ifa=0; ifa<nfa; ifa++)
  {
    if (ifa >= nfa_b)
    {
      int icv0 = cvofa[ifa][0];
      int noc00 = nbocv_i[icv0];

      int icv1 = cvofa[ifa][1];
      int noc11 = nbocv_i[icv1];

      A[noc00] -= massFlux_fa[ifa];

      if (icv1 < ncv)
        A[noc11] += massFlux_fa[ifa];
    }
    else
    {
      int icv0 = cvofa[ifa][0];
      int noc00 = nbocv_i[icv0];

      A[noc00] -= massFlux_fa[ifa];
    }
  }

  // under-relaxation...
  for (int icv = 0; icv < ncv; icv++)
  {
    int noc = nbocv_i[icv];
    A[noc] /= scal->relax;
    rhs[icv] += (1.0 - scal->relax)*A[noc]*scal->phi[icv];
  }


  // solve the system
  solveCvScalarBcgstab(scal->phi, A, rhs, scal->phiZero, scal->phiZeroRel, scal->phiMaxiter, scal->getName());



  // bound the scalar for stability reasons
  for (int icv=0; icv<ncv; icv++)
    scal->phi[icv] = min(max(scal->phi[icv], scal->lowerBound), scal->upperBound);


  UserDefinedScalarClipping(scal->getName());

  // update into halo cells
  updateCvData(scal->phi, REPLACE_DATA);

}


/**
 *  solveScalarTransport: steady (no time term) implicit
 */
void UgpWithCvCompFlow::solveScalarTransportImplUnsteady(ScalarTranspEq *scal, double *rhs, double *A, double *massFlux_fa)
{
  // initialize LHS and RHS and everything else
  for (int noc=0; noc<nbocv_s; noc++)     A[noc] = 0.0;
  for (int icv=0; icv<ncv; icv++)         rhs[icv] = 0.0;

  // -----------------------------------------------------------
  // define values at the boundary first (rho*phi)_fa
  // -----------------------------------------------------------
  string scalName(scal->getName());
  double phiBCval = 0.0, phiBCflux = 0.0;
  
  double *phi_fa = scal->phi_bfa;

  // -----------------------------------------------------------
  // then calculate gradient with boundary conditions
  // -----------------------------------------------------------

  // calc gradient of phi for viscous flux and second order
  static double (*gradPhi)[3] = new double[ncv_g][3];
  calcCvScalarGrad(gradPhi, scal->phi, phi_fa, gradreconstruction, limiterNavierS, scal->phi, epsilonSDWLS);

  // -----------------------------------------------------------
  // calc user defined diffusivity
  // -----------------------------------------------------------
  if (scal->diffTerm)
  {
    diffusivityHookScalarRansTurb(string(scal->getName()));
    diffusivityHookScalarRansComb(string(scal->getName()));
  }

  // ====================================================================
  // Part 1: interior face-based fluxes [nfa_b:nfa] ...
  // ====================================================================
  for (int ifa=nfa_b; ifa<nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];

    // ####################################### convective flux #######################################
    double lamPl = 0.0, lamMi = 0.0, explConv = 0.0;
    if (scal->convTerm)
    {
      lamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
      lamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));

      if ((sndOrder == true) && (stOrderScalar == false))
      {
        double dx0[3],dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

        explConv = sndOrder*(lamPl*vecDotVec3d(gradPhi[icv0], dx0) + lamMi*vecDotVec3d(gradPhi[icv1], dx1));
      }
    }

    // ####################################### viscous flux #######################################
    double explVisc = 0.0, implVisc = 0.0;
    if (scal->diffTerm)
    {
      double n[3], s[3];
      double nmag = normVec3d(n, fa_normal[ifa]);

      vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(s);
      double alpha = vecDotVec3d(n, s);

      double diff = scal->diff[ifa];
      implVisc = -diff*nmag*alpha/smag;

      for (int i=0; i<3; i++)
        explVisc += (gradPhi[icv0][i] + gradPhi[icv1][i])*(n[i] - alpha*s[i]);
      explVisc *= -0.5*diff*nmag;
    }

    //####################################### build equation system #######################################
    int noc00, noc01, noc11, noc10;
    getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    A[noc00]  += lamPl - implVisc;
    A[noc01]  += lamMi + implVisc;
    rhs[icv0] -= explConv + explVisc;

    if (icv1 < ncv)
    {
      A[noc11]  += -lamMi - implVisc;
      A[noc10]  += -lamPl + implVisc;
      rhs[icv1] += explConv + explVisc;
    }
  }



  // ====================================================================
  // Part 2: boundary face-based fluxes
  // ====================================================================

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int noc00 = nbocv_i[icv0];

          if (scal->convTerm)
          {
            if (massFlux_fa[ifa] > 0.0)     A[noc00]  += massFlux_fa[ifa];
            else                            rhs[icv0] -= massFlux_fa[ifa]*phi_fa[ifa];
          }

          if (scal->diffTerm)
          {
            double n[3], s_half[3];
            double nmag = normVec3d(n, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
            double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);
            
            double diff = scal->diff[ifa];

            double impl = -diff*nmag/smag_half;
            double expl = 0.0;  // no need for a tangential correction

            A[noc00]  -=  impl;
            rhs[icv0] +=  expl - impl*phi_fa[ifa];
          }
        }
      }
      else {
        if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int noc00 = nbocv_i[icv0];

            double phiBCvalTMP = phiBCval;

            if (scal->convTerm)
            {
              if (massFlux_fa[ifa] > 1.0e-8)
              {
                phiBCvalTMP = scal->phi[icv0];
                A[noc00] += massFlux_fa[ifa];
              }
              else
              {
                rhs[icv0] -= massFlux_fa[ifa]*phiBCvalTMP;
              }
            }

            if (scal->diffTerm)
            {
              double n[3], s_half[3];
              double nmag = normVec3d(n, fa_normal[ifa]);
              vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
              double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);

              double diff = scal->diff[ifa];

              double coeff = diff*nmag/smag_half;
              double impl = -coeff;
              double expl = 0.0;  // no need for a tangential correction

              A[noc00]  -=  impl;
              rhs[icv0] +=  expl + coeff*phiBCvalTMP;
            }
          }
        }
        else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
        {
          cerr << "scalarZoneIsFlux not implemented yet!" << endl;
          throw(-1);
        }
        else // apply zero flux if no boundary condition is specified
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is symmetry / zero flux" << endl;
        }
      }
    }
  }


  // ====================================================================
  // Part 3: add source term .......
  // ====================================================================

  sourceHookScalarRansTurb(rhs, A, scal->getName());
  sourceHookScalarRansComb(rhs, A, scal->getName());

  calcResidual(scal->resAve, scal->resMax, scal->phi, A, rhs);


  // correction for steady state
  // d(rho phi)/dt + d(rho phi u)/dx = ....
  // rho d(phi)/dt + phi d(rho)/dt = +d(rho phi u)/dt = ....
  //
  //      ^                  ^
  //      |                  |
  //                         |
  //  zero (steady state)    |
  //                         |
  //
  //                  d(rho)/dt = -d(rho u_i)/dx_i
  //
  for (int ifa=0; ifa<nfa; ifa++)
  {
    if (ifa >= nfa_b)
    {
      int icv0 = cvofa[ifa][0];
      int noc00 = nbocv_i[icv0];

      int icv1 = cvofa[ifa][1];
      int noc11 = nbocv_i[icv1];

      A[noc00] -= massFlux_fa[ifa];

      if (icv1 < ncv)
        A[noc11] += massFlux_fa[ifa];
    }
    else
    {
      int icv0 = cvofa[ifa][0];
      int noc00 = nbocv_i[icv0];

      A[noc00] -= massFlux_fa[ifa];
    }
  }
  
  // time terms 
  for (int icv = 0; icv < ncv; icv++)
  {
    int noc = nbocv_i[icv];
    A[noc] += rho[icv]*cv_volume[icv]/local_dt[icv];
    rhs[icv] += rho[icv]*cv_volume[icv]/local_dt[icv]*scal->phi[icv];
  }

/*  // under-relaxation...
  for (int icv = 0; icv < ncv; icv++)
  {
    int noc = nbocv_i[icv];
    A[noc] /= scal->relax;
    rhs[icv] += (1.0 - scal->relax)*A[noc]*scal->phi[icv];
  }*/


  // solve the system
  solveCvScalarBcgstab(scal->phi, A, rhs, scal->phiZero, scal->phiZeroRel, scal->phiMaxiter, scal->getName());



  // bound the scalar for stability reasons
  for (int icv=0; icv<ncv; icv++)
    scal->phi[icv] = min(max(scal->phi[icv], scal->lowerBound), scal->upperBound);


  UserDefinedScalarClipping(scal->getName());

  // update into halo cells
  updateCvData(scal->phi, REPLACE_DATA);

}




/**
 *  solveScalarTransport:
 */
void UgpWithCvCompFlow::solveScalarTransportExplicit(ScalarTranspEq *scal, double *rhs, double *massFlux_fa)
{
  // initialize LHS and RHS and everything else
  for (int icv=0; icv<ncv; icv++)         rhs[icv] = 0.0;
  
  // -----------------------------------------------------------
  // define values at the boundary first phi_fa
  // -----------------------------------------------------------
  string scalName(scal->getName());
  double phiBCval = 0.0, phiBCflux = 0.0;
  
  double *phi_fa = scal->phi_bfa;

  // -----------------------------------------------------------
  // then calculate gradient with boundary conditions
  // -----------------------------------------------------------
  
  // calc gradient of phi for viscous flux and second order 
  static double (*gradPhi)[3] = new double[ncv_g][3];
  calcCvScalarGrad(gradPhi, scal->phi, phi_fa, gradreconstruction, limiterNavierS, scal->phi, epsilonSDWLS);


  // user defined diffusivity
  if (scal->diffTerm)
  {
    diffusivityHookScalarRansTurb(string(scal->getName()));
    diffusivityHookScalarRansComb(string(scal->getName()));
  }


  // ====================================================================
  // Part 1: interior face-based fluxes [nfa_b:nfa]...
  // ====================================================================

  double factSndOrder = 0.0;
  if ((sndOrder == true) && (stOrderScalar == false))
    factSndOrder = 1.0;

  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert((icv0 >= 0) && (icv0 < ncv));
    assert((icv1 >= 0) && (icv1 < ncv_g));

    double n[3], s[3];
    double nmag = normVec3d(n, fa_normal[ifa]);
    vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
    double smag = normVec3d(s);
    double alpha = vecDotVec3d(n, s);
    assert((alpha > 0.0)&&(alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...

    double dx0[3],dx1[3];
    vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
    vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

    // ####################################### convective flux #######################################
    double rhoLamPl = 0.0, rhoLamMi = 0.0, explConv = 0.0;
    if (scal->convTerm)
    {
      rhoLamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
      rhoLamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));
      explConv = rhoLamPl*(scal->phi[icv0]+factSndOrder*vecDotVec3d(gradPhi[icv0], dx0))
               + rhoLamMi*(scal->phi[icv1]+factSndOrder*vecDotVec3d(gradPhi[icv1], dx1));
    }

    // ####################################### viscous flux #######################################
    double explVisc = 0.0, implVisc = 0.0;
    if (scal->diffTerm)
    {
      double diff = scal->diff[ifa];
      explVisc = diff*nmag*(alpha*(scal->phi[icv1] - scal->phi[icv0])/smag
                     + 0.5*(gradPhi[icv0][0] + gradPhi[icv1][0])*(n[0] - alpha*s[0])
                     + 0.5*(gradPhi[icv0][1] + gradPhi[icv1][1])*(n[1] - alpha*s[1])
                     + 0.5*(gradPhi[icv0][2] + gradPhi[icv1][2])*(n[2] - alpha*s[2]));
    }

    //####################################### build equation system #######################################
    rhs[icv0] += -explConv + explVisc;

    if (icv1 < ncv)
      rhs[icv1] += explConv - explVisc;

  }


  // ====================================================================
  // Part 2: boundary face-based fluxes
  // ====================================================================

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        // user defined boundary hook
        boundaryHookScalarRansTurb(phi_fa, &(*zone), scalName);
        boundaryHookScalarRansComb(phi_fa, &(*zone), scalName);

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int noc00 = nbocv_i[icv0];

          if (scal->convTerm)
          {
            double rhoLamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
            double rhoLamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));
            double explConv = rhoLamPl*(scal->phi[icv0]) + rhoLamMi*(phi_fa[ifa]);
            rhs[icv0] -= explConv;
          }

          if (scal->diffTerm)
          {
            double n[3], s_half[3];
            double nmag = normVec3d(n, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
            double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);
            //double alpha = vecDotVec3d(n, s_half);    // HACK alfa should be always 1 at the boundary (provide such a grid!)

            double diff = scal->diff[ifa];
            double viscFlux = diff*nmag*(phi_fa[ifa]-scal->phi[icv0])/smag_half;
            rhs[icv0] +=  viscFlux;
          }
        }
      }
      else {
        if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int noc00 = nbocv_i[icv0];

            if (scal->convTerm)
            {
              double rhoLamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
              double rhoLamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));
              double explConv = rhoLamPl*(scal->phi[icv0]) + rhoLamMi*(phiBCval);
              rhs[icv0] -= explConv;
            }

            if (scal->diffTerm)
            {
              double n[3], s_half[3];
              double nmag = normVec3d(n, fa_normal[ifa]);
              vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
              double smag_half = fabs(vecDotVec3d(s_half, n));//normVec3d(s_half);
              //double alpha = vecDotVec3d(n, s_half);    // HACK alfa should be always 1 at the boundary (provide such a grid!)

              double diff = scal->diff[ifa];
              double viscFlux = diff*nmag*(phiBCval-scal->phi[icv0])/smag_half;
              rhs[icv0] +=  viscFlux;
            }
          }
        }
        else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
        {
          cerr << "scalarZoneIsFlux not implemented yet!" << endl;
          throw(-1);
        }
        else //(zoneIsSymmetry(zone->getName()))
        {
          if ((step == 1) && (mpi_rank == 0))
            cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int noc00 = nbocv_i[icv0];

            double dx0[3];
            vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
            double bv = scal->phi[icv0] + factSndOrder*vecDotVec3d(gradPhi[icv0], dx0);

            if (scal->convTerm)
            {
              double rhoLamPl = 0.5*(massFlux_fa[ifa] + fabs(massFlux_fa[ifa]));
              double rhoLamMi = 0.5*(massFlux_fa[ifa] - fabs(massFlux_fa[ifa]));
              double explConv = rhoLamPl*(scal->phi[icv0]) + rhoLamMi*(bv);
              rhs[icv0] -= explConv;
            }
          }
        }
      }
    }
  }

  // ====================================================================
  // Part 3: add source term .......
  // ====================================================================

  sourceHookScalarRansTurbExplicit(rhs, scal->getName());
  sourceHookScalarRansCombExplicit(rhs, scal->getName());

  // add    -phi*drho/dt  to rhs
  for (int ifa=0; ifa<nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    rhs[icv0] += scal->phi[icv0]*massFlux_fa[ifa];

//    cout << rhs[icv0] << "\t";
 
    if (ifa >= nfa_b)
    {
      int icv1 = cvofa[ifa][1];
      if (icv1 < ncv)  
        rhs[icv1] -= scal->phi[icv1]*massFlux_fa[ifa];
    }
  }


  // advance scalar and compute RESIDUAL
  scal->resAve = scal->resMax = 0.0;
  for (int icv = 0; icv < ncv; icv++)
  {
    scal->resAve = scal->resMax = fabs(rhs[icv]);

    double tmp = local_dt[icv]/cv_volume[icv];
    scal->phi[icv] += tmp/rho[icv]*rhs[icv];
  }

  // bound the scalar for stability reasons
  for (int icv=0; icv<ncv; icv++)
    scal->phi[icv] = min(max(scal->phi[icv], scal->lowerBound), scal->upperBound);

  // update into halo cells
  updateCvData(scal->phi, REPLACE_DATA);
}

void UgpWithCvCompFlow::calcViscousFluxScalar_new(double *rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, int flagImplicit)
{
  if (!transpScal.diffTerm)
    return;

  string scalName(transpScal.getName());
  double phiBCval = 0.0, phiBCflux = 0.0;

  double *phi = transpScal.phi;
  double *phi_bfa = transpScal.phi_bfa;
  double (*grad_phi)[3] = transpScal.grad_phi;
  double *diff = transpScal.diff;

  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  calcCvScalarGrad(grad_phi, phi, phi_bfa, gradreconstruction, limiterNavierS, phi, epsilonSDWLS);

  // =============================================================================================
  // compute user-defined diffusivity
  // =============================================================================================
  diffusivityHookScalarRansTurb(scalName);
  diffusivityHookScalarRansComb(scalName);

  // =============================================================================================
  // cycle trough internal faces and compute viscous flux and implicit matrix
  // =============================================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert((icv0 >= 0) && (icv0 < ncv));
    assert((icv1 >= 0) && (icv1 < ncv_g));
    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    double n[3] = {0.0, 0.0, 0.0}, s[3] = {0.0, 0.0, 0.0};
    double nmag = normVec3d(n, fa_normal[ifa]);
    vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
    double smag = normVec3d(s);
    double alpha = vecDotVec3d(n, s);
    assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...

    double viscFlux = diff[ifa] * nmag * (alpha * (phi[icv1] - phi[icv0]) / smag
                    + 0.5 * (grad_phi[icv0][0] + grad_phi[icv1][0]) * (n[0] - alpha * s[0])
                    + 0.5 * (grad_phi[icv0][1] + grad_phi[icv1][1]) * (n[1] - alpha * s[1])
                    + 0.5 * (grad_phi[icv0][2] + grad_phi[icv1][2]) * (n[2] - alpha * s[2]));

    rhs_rhoScal[icv0] += viscFlux;

    if (icv1 < ncv)
      rhs_rhoScal[icv1] -= viscFlux;

    if (flagImplicit)
    {
      AScal[noc00] +=   diff[ifa] * nmag * alpha / smag / rho[icv0];
      AScal[noc01] += - diff[ifa] * nmag * alpha / smag / rho[icv1];

      if (icv1 < ncv)
      {
        AScal[noc11] +=   diff[ifa] * nmag * alpha / smag / rho[icv1];
        AScal[noc10] += - diff[ifa] * nmag * alpha / smag / rho[icv0];
      }
    }


  }

  // =============================================================================================
  // cycle trough boundary faces and compute viscous flux and implicit matrix
  // =============================================================================================
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      // .............................................................................................
      // HOOK BOUNDARY CONDITION
      // .............................................................................................
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int noc00 = nbocv_i[icv0];
          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          double viscFlux = diff[ifa] * nmag * (phi_bfa[ifa] - phi[icv0]) / smag_half;
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit)
            AScal[noc00] += diff[ifa] * nmag / smag_half / rho[icv0];
        }
      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int noc00 = nbocv_i[icv0];
          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          double viscFlux = diff[ifa] * nmag * (phi_bfa[ifa] - phi[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit)
            AScal[noc00] += diff[ifa] * nmag / smag_half / rho[icv0];
        }
      }
      // .............................................................................................
      // FLUX BOUNDARY CONDITION: NOT IMPLEMENTED YET
      // .............................................................................................
      else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
      {
        cerr << "scalarZoneIsFlux not implemented yet!" << endl;
        throw(-1);
      }
      // .............................................................................................
      // OTHER BOUNDARY CONDITIONS: NO VISCOUS FLUX
      // .............................................................................................
      else
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
      }
    }
  }

}



