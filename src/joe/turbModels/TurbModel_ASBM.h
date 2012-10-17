#ifndef RANSTURBMODEL_ASBM_H
#define RANSTURBMODEL_ASBM_H

#include "TurbModel_KOMSST.h"
#include "TurbModel_V2F.h"


typedef RansTurbKOmSST TURB_MOD_FOR_ASBM; // define the turbulence model used with the ASBM
//typedef RansTurbV2F TURB_MOD_FOR_ASBM;

extern "C"{
  void asbm_(
      double*, double*, double*, double*, double*, double*,
      double*, double*, double*, int*, int*);
}

class RansTurbASBM : virtual public TURB_MOD_FOR_ASBM
{
public:   // constructors

  RansTurbASBM()
  {
    if (mpi_rank == 0)
      cout << "RansTurbASBM()" << endl;

    // General variables
    st_diag       = NULL;     registerVector(st_diag,       "st_diag",       CV_DATA);
    st_offdiag    = NULL;     registerVector(st_offdiag,    "st_offdiag",    CV_DATA);
    wt_offdiag    = NULL;     registerVector(wt_offdiag,    "wt_offdiag",    CV_DATA);

    as_diag       = NULL;     registerVector(as_diag,       "as_diag",       CV_DATA);
    as_offdiag    = NULL;     registerVector(as_offdiag,    "as_offdiag",    CV_DATA);
    ar_diag       = NULL;     registerVector(ar_diag,       "ar_diag",       CV_DATA);
    ar_offdiag    = NULL;     registerVector(ar_offdiag,    "ar_offdiag",    CV_DATA);

    etar          = NULL;     registerScalar(etar,          "etar",          CV_DATA);
    etaf          = NULL;     registerScalar(etaf,          "etaf",          CV_DATA);
    a2            = NULL;     registerScalar(a2,            "a2",            CV_DATA);

    scal_phi      = NULL;     registerScalar(scal_phi,      "scal_phi",      CV_DATA);
    scal_bet      = NULL;     registerScalar(scal_bet,      "scal_bet",      CV_DATA);
    scal_chi      = NULL;     registerScalar(scal_chi,      "scal_chi",      CV_DATA);

    // For debugging only
    hatwt         = NULL;     registerScalar(hatwt,         "hatwt",         CV_DATA);
    hatst         = NULL;     registerScalar(hatst,         "hatst",         CV_DATA);
    ststa         = NULL;     registerScalar(ststa,         "ststa",         CV_DATA);

    // Blocking variables
    block_diag    = NULL;     registerVector(block_diag,    "block_diag",    CV_DATA);
    block_offdiag = NULL;     registerVector(block_offdiag, "block_offdiag", CV_DATA);
    bphi          = NULL;     registerScalar(bphi,          "bphi",          CV_DATA);
    bphi_bfa      = NULL;     // this is a face array
    grad_bphi     = NULL;     // this array includes ghost cells
  }

  virtual ~RansTurbASBM() {}

protected:   // member vars

  // General variables
  int      start_asbm;            // iteration number at which ASBM starts

  double   (*st_diag)[3];         // diagonal rate of strain tensor times tau
  double   (*st_offdiag)[3];      // off-diagonal rate of strain tensor times tau
  double   (*wt_offdiag)[3];      // off-diagonal rate of rotation tensor times tau

  double   (*as_diag)[3];         // diagonal eddy axis tensor from st
  double   (*as_offdiag)[3];      // off-diagonal eddy axis tensor from st
  double   (*ar_diag)[3];         // diagonal eddy axis tensor from st and wt
  double   (*ar_offdiag)[3];      // off-diagonal eddy axis tensor from st and wt

  double   *etar;
  double   *etaf;
  double   *a2;

  double   *scal_phi;
  double   *scal_bet;
  double   *scal_chi;

  // For debugging only
  double   *hatwt;
  double   *hatst;
  double   *ststa;

  // Blocking variables
  double   (*block_diag)[3];      // diagonal blockage tensor
  double   (*block_offdiag)[3];   // off-diagonal blockage tensor
  double   *bphi;                 // phi for wall blocking
  double   *bphi_bfa;             // phi for wall blocking at boundaries
  double   (*grad_bphi)[3];       // gradients for wall blocking phi

  int      block_frq;             // frequency of blocking computation
  int      block_maxIter;         // linear solver parameters
  double   block_zeroAbs;
  double   block_zeroRel;

  // I/O variables
  FILE *finfo;              // output file for ASBM quantities

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    TURB_MOD_FOR_ASBM::initialHookScalarRansTurbModel();

    if (mpi_rank == 0)
      cout << "initialHook ASBM" << endl;

    start_asbm = getIntParam("START_ASBM", "0");
    block_frq = getIntParam("BLOCK_FRQ","10");

    if (!checkParam("LINEAR_SOLVER_BLOCK_TRESHOLDS"))
    {
      // add default values
      ParamMap::add("LINEAR_SOLVER_BLOCK_TRESHOLDS  MAX_ITER=500  ABS_RESID=1.0e-8  REL_RESID=1.0e-8");
      if (mpi_rank == 0)
        cout << "WARNING: added keyword"
             << "\"LINEAR_SOLVER_BLOCK_TRESHOLDS  MAX_ITER=500  ABS_RESID=1.0e-8  REL_RESID=1.0e-8\""
             << " to parameter map!" << endl;
    }
    block_maxIter = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getInt("MAX_ITER");
    block_zeroAbs = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getDouble("ABS_RESID");
    block_zeroRel = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getDouble("REL_RESID");

    if (mpi_rank == 0)
    {
      if ( (finfo=fopen("asbmInfo.dat","wt")) == NULL )
      {
        cerr << "Error: cannot open file asbmInfo.dat" << endl;
        throw(-1);
      }
    }

    // Initialize bphi, bphi_bfa, grad_bphi, and B tensor
    for (int icv = 0; icv < ncv; ++icv)
    {
        bphi[icv] = 0.0;

        block_diag[icv][0] = 0.0;
        block_diag[icv][1] = 0.0;
        block_diag[icv][2] = 0.0;

        block_offdiag[icv][0] = 0.0;
        block_offdiag[icv][1] = 0.0;
        block_offdiag[icv][2] = 0.0;
    }

    bphi_bfa  = new double[nfa_b];
    grad_bphi = new double[ncv_g][3];
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              bphi_bfa[ifa] = 1.0;
            }
          }
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              bphi_bfa[ifa] = bphi[icv0];
            }
          }
        }
      }
    }
    calcCvScalarGrad(grad_bphi, bphi, bphi_bfa, gradreconstruction, limiterNavierS, bphi, epsilonSDWLS);
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (bphi_bfa != NULL)  {delete [] bphi_bfa;     bphi_bfa = NULL;}
    if (grad_bphi != NULL) {delete [] grad_bphi;    grad_bphi = NULL;}
    if (mpi_rank == 0) fclose(finfo);
  }

  virtual void calcRansTurbViscMuet()
  {
    TURB_MOD_FOR_ASBM::calcRansTurbViscMuet();

    // Non-linear domain partitioning
    if (step == start_asbm) setNonLinearDomain();

    // Computation of ASBM variables
    if (step >= start_asbm){
        if ( step%block_frq == 0 ) calcBlockTensor();
        calcRsCenterASBM();
    }
  }

  virtual void setNonLinearDomain(){
    // default
    if (mpi_rank == 0)
      cout << "NON-LINEAR DOMAIN SET AS DEFAULT." << endl;
    for (int ifa = 0; ifa < nfa; ifa++)
      nonLinear[ifa] = 1.0;
  }

  void calcRsCenterASBM(){
    // ====================================================================
    //        variable declaration and initialization
    // ====================================================================

    // Input variables
    double WFT[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
    double ST[3][3], WT[3][3], BL[3][3];
    double EMPTY[3][3] = {{0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0}};

    // In and out variables
    double AS[3][3], AR[3][3];

    // Output variables
    double REY[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
    double DIM[3][3], CIR[3][3];

    // Error flags
    int ierr = 0, myNerr = 0, Nerr = 0, maxitrs = 99;

    // ====================================================================
    // compute inputs
    // ====================================================================
    double divg, tau;
    for (int icv = 0; icv < ncv; icv++){
        divg = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];
        tau  = turbTS[icv];

        st_diag[icv][0] = (grad_u[icv][0][0] - 1.0/3.0*divg)*tau;
        st_diag[icv][1] = (grad_u[icv][1][1] - 1.0/3.0*divg)*tau;
        st_diag[icv][2] = (grad_u[icv][2][2] - 1.0/3.0*divg)*tau;

        st_offdiag[icv][0] = 0.5*(grad_u[icv][0][1] + grad_u[icv][1][0])*tau;
        st_offdiag[icv][1] = 0.5*(grad_u[icv][0][2] + grad_u[icv][2][0])*tau;
        st_offdiag[icv][2] = 0.5*(grad_u[icv][1][2] + grad_u[icv][2][1])*tau;

        wt_offdiag[icv][0] = 0.5*(grad_u[icv][0][1] - grad_u[icv][1][0])*tau;
        wt_offdiag[icv][1] = 0.5*(grad_u[icv][0][2] - grad_u[icv][2][0])*tau;
        wt_offdiag[icv][2] = 0.5*(grad_u[icv][1][2] - grad_u[icv][2][1])*tau;
    }

    updateCvData(st_diag, REPLACE_ROTATE_DATA);
    updateCvData(st_offdiag, REPLACE_ROTATE_DATA);
    updateCvData(wt_offdiag, REPLACE_ROTATE_DATA);

    // Smooth inputs
    //smoothingVec(st_diag);
    //smoothingVec(st_offdiag);
    //smoothingVec(wt_offdiag);

    // ====================================================================
    // Compute cell-centered Reynolds stresses
    // ====================================================================
    for (int icv = 0; icv < ncv; icv++){

        // Anisotropic rate of strain tensor, sij = 0.5*(dui_dxj + duj_dxi) - 1/3*divg
        ST[0][0] = st_diag[icv][0];     ST[0][1] = st_offdiag[icv][0];     ST[0][2] = st_offdiag[icv][1];
        ST[1][0] = ST[0][1];            ST[1][1] = st_diag[icv][1];        ST[1][2] = st_offdiag[icv][2];
        ST[2][0] = ST[0][2];            ST[2][1] = ST[1][2];               ST[2][2] = st_diag[icv][2];

        // Rate of mean rotation tensor, C TO FORTRAN, NEED TRANSPOSE: wij = -0.5*(dui_dxj - duj_dxi)
        WT[0][0] = 0.0;                 WT[0][1] = -wt_offdiag[icv][0];    WT[0][2] = -wt_offdiag[icv][1];
        WT[1][0] = -WT[0][1];           WT[1][1] = 0.0;                    WT[1][2] = -wt_offdiag[icv][2];
        WT[2][0] = -WT[0][2];           WT[2][1] = -WT[1][2];              WT[2][2] = 0.0;

        // Blockage tensor
        BL[0][0] = block_diag[icv][0];  BL[0][1] = block_offdiag[icv][0];  BL[0][2] = block_offdiag[icv][1];
        BL[1][0] = BL[0][1];            BL[1][1] = block_diag[icv][1];     BL[1][2] = block_offdiag[icv][2];
        BL[2][0] = BL[0][2];            BL[2][1] = BL[1][2];               BL[2][2] = block_diag[icv][2];

        // Strained and rotated A's
        AS[0][0] = as_diag[icv][0];     AS[0][1] = as_offdiag[icv][0];     AS[0][2] = as_offdiag[icv][1];
        AS[1][0] = AS[0][1];            AS[1][1] = as_diag[icv][1];        AS[1][2] = as_offdiag[icv][2];
        AS[2][0] = AS[0][2];            AS[2][1] = AS[1][2];               AS[2][2] = as_diag[icv][2];

        AR[0][0] = ar_diag[icv][0];     AR[0][1] = ar_offdiag[icv][0];     AR[0][2] = ar_offdiag[icv][1];
        AR[1][0] = AR[0][1];            AR[1][1] = ar_diag[icv][1];        AR[1][2] = ar_offdiag[icv][2];
        AR[2][0] = AR[0][2];            AR[2][1] = AR[1][2];               AR[2][2] = ar_diag[icv][2];

        // Call the ASBM
        ierr = 0;

        asbm_(&(ST[0][0]), &(WT[0][0]), &(WFT[0][0]), &(BL[0][0]), &(AS[0][0]), &(AR[0][0]),
              &(REY[0][0]), &(DIM[0][0]), &(CIR[0][0]), &maxitrs, &ierr);

        if (ierr != 0)
        {
          myNerr ++;
          cout << "ASBM error number " << ierr << " in cell " << icv << endl;
          cout << "Cell-x: " << x_cv[icv][0]
               << " Cell-y: " << x_cv[icv][1]
               << " Cell-z: " << x_cv[icv][2] << endl;
        }

        rij_diag[icv][0] = -REY[0][0]*2*kine[icv]*rho[icv];
        rij_diag[icv][1] = -REY[1][1]*2*kine[icv]*rho[icv];
        rij_diag[icv][2] = -REY[2][2]*2*kine[icv]*rho[icv];

        rij_offdiag[icv][0] = -REY[0][1]*2*kine[icv]*rho[icv];
        rij_offdiag[icv][1] = -REY[0][2]*2*kine[icv]*rho[icv];
        rij_offdiag[icv][2] = -REY[1][2]*2*kine[icv]*rho[icv];

        as_diag[icv][0] = AS[0][0];       as_offdiag[icv][0] = AS[0][1];
        as_diag[icv][1] = AS[1][1];       as_offdiag[icv][1] = AS[0][2];
        as_diag[icv][2] = AS[2][2];       as_offdiag[icv][2] = AS[1][2];

        ar_diag[icv][0] = AR[0][0];       ar_offdiag[icv][0] = AR[0][1];
        ar_diag[icv][1] = AR[1][1];       ar_offdiag[icv][1] = AR[0][2];
        ar_diag[icv][2] = AR[2][2];       ar_offdiag[icv][2] = AR[1][2];

        scal_phi[icv] = CIR[0][0];
        scal_bet[icv] = CIR[1][0];
        scal_chi[icv] = CIR[2][0];

        etar[icv] = CIR[0][1];
        etaf[icv] = CIR[1][1];
        a2[icv]   = CIR[2][1];

        hatwt[icv] = CIR[0][2];
        hatst[icv] = CIR[1][2];
        ststa[icv] = CIR[2][2];
    }
    updateCvData(etar, REPLACE_DATA);
    updateCvData(etaf, REPLACE_DATA);
    updateCvData(a2, REPLACE_DATA);

    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);

    // Smooth Outputs
    smoothingVec(rij_diag);
    smoothingVec(rij_offdiag);

    // ====================================================================
    // Compute internal face Reynolds stresses
    // ====================================================================
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert( icv0 >= 0 );
      assert( icv1 >= 0 );

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      // Update Reynolds stresses
      rij_diag_fa[ifa][0] = (w1*rij_diag[icv0][0] + w0*rij_diag[icv1][0])/(w0 + w1);
      rij_diag_fa[ifa][1] = (w1*rij_diag[icv0][1] + w0*rij_diag[icv1][1])/(w0 + w1);
      rij_diag_fa[ifa][2] = (w1*rij_diag[icv0][2] + w0*rij_diag[icv1][2])/(w0 + w1);

      rij_offdiag_fa[ifa][0] = (w1*rij_offdiag[icv0][0] + w0*rij_offdiag[icv1][0])/(w0 + w1);
      rij_offdiag_fa[ifa][1] = (w1*rij_offdiag[icv0][1] + w0*rij_offdiag[icv1][1])/(w0 + w1);
      rij_offdiag_fa[ifa][2] = (w1*rij_offdiag[icv0][2] + w0*rij_offdiag[icv1][2])/(w0 + w1);
    }

    // ====================================================================
    // Compute boundary face Reynolds stresses
    // ====================================================================
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          // .............................................................................................
          // SYMMETRY BOUNDARY CONDITION
          // .............................................................................................
          if (param->getString() == "SYMMETRY")
          {
            // No viscous flux in this case (or yes!)

            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );

              // define Reynolds stresses
              rij_diag_fa[ifa][0] = rij_diag[icv0][0];
              rij_diag_fa[ifa][1] = rij_diag[icv0][1];
              rij_diag_fa[ifa][2] = rij_diag[icv0][2];

              rij_offdiag_fa[ifa][0] = rij_offdiag[icv0][0];
              rij_offdiag_fa[ifa][1] = rij_offdiag[icv0][1];
              rij_offdiag_fa[ifa][2] = rij_offdiag[icv0][2];

            }
          }
          // .............................................................................................
          // WALL BOUNDARY CONDITION
          // .............................................................................................
          else if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              rij_diag_fa[ifa][0] = 0.0;
              rij_diag_fa[ifa][1] = 0.0;
              rij_diag_fa[ifa][2] = 0.0;

              rij_offdiag_fa[ifa][0] = 0.0;
              rij_offdiag_fa[ifa][1] = 0.0;
              rij_offdiag_fa[ifa][2] = 0.0;
            }
          }
          // .............................................................................................
          // OTHER BOUNDARY CONDITIONS
          // .............................................................................................
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );

              // define Reynolds stresses
              rij_diag_fa[ifa][0] = rij_diag[icv0][0];
              rij_diag_fa[ifa][1] = rij_diag[icv0][1];
              rij_diag_fa[ifa][2] = rij_diag[icv0][2];

              rij_offdiag_fa[ifa][0] = rij_offdiag[icv0][0];
              rij_offdiag_fa[ifa][1] = rij_offdiag[icv0][1];
              rij_offdiag_fa[ifa][2] = rij_offdiag[icv0][2];
            }
          }
        }
      }
    }

    // ====================================================================
    // Output error file
    // ====================================================================
    MPI_Reduce(&myNerr,&Nerr,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
    {
      if (finfo == NULL)
      {
        cerr << "Error: cannot write to file asbmInfo.dat" << endl;
        throw(-1);
      }
      else
      {
        fprintf(finfo,"%8d\t%8d\n",step,Nerr);
      }
    }
  }

  void calcBlockTensor()
  {
    // =========================================================================================
    // allocate static memory
    // =========================================================================================

    static double *A_block   = new double[nbocv_s];
    static double *rhs_block = new double[ncv];

    for (int icv = 0; icv < ncv; ++icv)   rhs_block[icv] = 0.0;
    for (int noc = 0; noc < nbocv_s; noc++) A_block[noc] = 0.0;

    // =========================================================================================
    // compute A_block and rhs_block
    // =========================================================================================

    // ..........................................................................................
    // cycle trough internal faces first
    // ..........................................................................................
    for (int ifa = nfa_b; ifa < nfa; ++ifa)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert( icv0 >= 0 );
      assert( icv1 >= 0 );

      // viscous flux
      double sVec[3], nVec[3];
      double area = normVec3d(nVec, fa_normal[ifa]);
      vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
      double ds = normVec3d(sVec);

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      double grad_f[3], lim_grad[3];
      for (int i = 0; i < 3; i++)
        grad_f[i] = (w1*grad_bphi[icv0][i] + w0*grad_bphi[icv1][i])/(w0 + w1);
      double proj_grad = vecDotVec3d(grad_f,sVec);
      for (int i = 0; i < 3; i++)
        lim_grad[i] = grad_f[i] - proj_grad*sVec[i];

      double explVisc = (lim_grad[0]*nVec[0] + lim_grad[1]*nVec[1] + lim_grad[2]*nVec[2])*area;
      double implVisc = vecDotVec3d(sVec, nVec)*area/ds;

      // build equation system
      int noc00, noc01, noc11, noc10;
      getImplDependencyIndex(noc00,noc01,noc11,noc10,icv0,icv1);

      A_block[noc00] += implVisc;
      A_block[noc01] -= implVisc;
      rhs_block[icv0] += explVisc;

      if (icv1 < ncv)
      {
        A_block[noc11] += implVisc;
        A_block[noc10] -= implVisc;
        rhs_block[icv1] -= explVisc;
      }
    }

    // ..........................................................................................
    // cycle trough boundary faces
    // ..........................................................................................
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );
              int noc00 = nbocv_i[icv0];

              double sVec[3], nVec[3];
              double area = normVec3d(nVec, fa_normal[ifa]);
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double ds = normVec3d(sVec);

              double proj_grad = vecDotVec3d(grad_bphi[icv0],sVec);
              double lim_grad[3];
              for (int i = 0; i < 3; i++)
                lim_grad[i] = grad_bphi[icv0][i] - proj_grad*sVec[i];

              double implVisc = vecDotVec3d(sVec, nVec)*area/ds;
              double explVisc = (lim_grad[0]*nVec[0] +
                                 lim_grad[1]*nVec[1] +
                                 lim_grad[2]*nVec[2])*area +
                                 bphi_bfa[ifa]*implVisc;
              //rhs_block[icv0] += phi_asbm_bfa[ifa]*impl_coeff;

              // build equation system
              A_block[noc00] += implVisc;
              rhs_block[icv0] += explVisc;
            }
          }
          else
          {
            // do nothing as flux is zero
          }
        }
      }
    }

    // ..........................................................................................
    // cycle through cell centroids to compute source terms
    // ..........................................................................................

    /*cout << "before centroids" << endl;
    // norm of A matrix
    double myAnorm = 0.0, Anorm = 0.0;
    for (int noc = 0; noc < nbocv_s; noc++)
      myAnorm += A_block[noc];
    MPI_Reduce(&myAnorm,&Anorm,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0){
        if (asbmInfo.is_open())
          asbmInfo << "Anorm before: " << Anorm << endl;
        else
          cout << "Unable to write to asbmInfo.dat file." << endl;
    }*/

    for (int icv = 0; icv < ncv; icv++)
    {
      int noc00 = nbocv_i[icv];
      A_block[noc00] += cv_volume[icv]/(turbLS[icv]*turbLS[icv]);
    }

    // =========================================================================================
    // under-relaxation
    // =========================================================================================
    //for (int icv = 0; icv < ncv; icv++){
    //    int noc = nbocv_i[icv];
    //    A_block[noc] /= 1.0;
    //    rhs_block[icv] += (1.0 - 1.0)*A_block[noc]*phi_asbm[icv];
    //}

    // =========================================================================================
    // solve the linear system
    // =========================================================================================

    char scalName[] = "bphi";

    static double *bphi_temp = new double[ncv_g];
    for (int icv = 0; icv < ncv; ++icv)
      bphi_temp[icv] = bphi[icv];

    /*// norm of A matrix
    myAnorm = 0.0; Anorm = 0.0;
    for (int noc = 0; noc < nbocv_s; noc++)
      myAnorm += A_block[noc];
    MPI_Reduce(&myAnorm,&Anorm,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0){
        if (asbmInfo.is_open())
          asbmInfo << "Anorm after: " << Anorm << endl;
        else
          cout << "Unable to write to asbmInfo.dat file." << endl;
    }*/

    solveLinSysScalar(bphi_temp, A_block, rhs_block, block_zeroAbs, block_zeroRel, block_maxIter, scalName);

    // clipping
    for (int icv = 0; icv < ncv; ++icv)
    {
      if (bphi_temp[icv] > 0)
        bphi[icv] = bphi_temp[icv];
      else
        bphi[icv] = 0.0;
    }
    updateCvData(bphi, REPLACE_DATA);

    // compute residual
    double thisRes, myRes = 0.0, Res = 0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      thisRes = rhs_block[icv];

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1] - 1;

      thisRes -= A_block[noc_f]*bphi[icv];
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        thisRes -= A_block[noc]*bphi[icv_nbr];
      }

      myRes += fabs(thisRes);
    }
    MPI_Reduce(&myRes,&Res,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

    // output for info
    if (mpi_rank == 0)
    {
      if (finfo == NULL)
      {
        cerr << "Error: cannot write to file asbmInfo.dat" << endl;
        throw(-1);
      }
      else
      {
        fprintf(finfo,"bphi resid: %12.6e\n",Res);
      }
    }

    // =========================================================================================
    // set the boundary conditions
    // =========================================================================================

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              bphi_bfa[ifa] = 1.0;
            }
          }
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              bphi_bfa[ifa] = bphi[icv0];
            }
          }
        }
      }
    }

    // =========================================================================================
    // compute gradients, with boundary values
    // =========================================================================================

    calcCvScalarGrad(grad_bphi, bphi, bphi_bfa, gradreconstruction, limiterNavierS, bphi, epsilonSDWLS);

    // =========================================================================================
    // compute the blockage tensor
    // =========================================================================================
    double div_bphi;
    for (int icv = 0; icv < ncv; icv++){
      div_bphi = 0.0;
      div_bphi += grad_bphi[icv][0]*grad_bphi[icv][0] +
                  grad_bphi[icv][1]*grad_bphi[icv][1] +
                  grad_bphi[icv][2]*grad_bphi[icv][2];

      if (div_bphi > 1.0e-5){
        block_diag[icv][0] = grad_bphi[icv][0]*grad_bphi[icv][0]/div_bphi*bphi[icv];
        block_diag[icv][1] = grad_bphi[icv][1]*grad_bphi[icv][1]/div_bphi*bphi[icv];
        block_diag[icv][2] = grad_bphi[icv][2]*grad_bphi[icv][2]/div_bphi*bphi[icv];

        block_offdiag[icv][0] = grad_bphi[icv][0]*grad_bphi[icv][1]/div_bphi*bphi[icv];
        block_offdiag[icv][1] = grad_bphi[icv][0]*grad_bphi[icv][2]/div_bphi*bphi[icv];
        block_offdiag[icv][2] = grad_bphi[icv][1]*grad_bphi[icv][2]/div_bphi*bphi[icv];
      }
      else {
        block_diag[icv][0] = 0.0;
        block_diag[icv][1] = 0.0;
        block_diag[icv][2] = 0.0;

        block_offdiag[icv][0] = 0.0;
        block_offdiag[icv][1] = 0.0;
        block_offdiag[icv][2] = 0.0;
      }
    }

    updateCvData(block_diag, REPLACE_DATA);
    updateCvData(block_offdiag, REPLACE_DATA);

  }

  void smoothingVec(double (*input)[3])
  {
    double vol_tot, vol_icv;
    for (int icv = 0; icv < ncv; icv++)
    {
      vol_tot = cv_volume[icv];

      input[icv][0] *= vol_tot;
      input[icv][1] *= vol_tot;
      input[icv][2] *= vol_tot;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];

        vol_icv = cv_volume[icv_nbr];
        vol_tot += vol_icv;

        input[icv][0] += input[icv_nbr][0]*vol_icv;
        input[icv][1] += input[icv_nbr][1]*vol_icv;
        input[icv][2] += input[icv_nbr][2]*vol_icv;
      }

      input[icv][0] /= vol_tot;
      input[icv][1] /= vol_tot;
      input[icv][2] /= vol_tot;
    }
  }

  void smoothingScal(double *input)
  {
    double vol_tot, vol_icv;
    for (int icv = 0; icv < ncv; icv++)
    {
      vol_tot = cv_volume[icv];

      input[icv] *= vol_tot;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];

        vol_icv = cv_volume[icv_nbr];
        vol_tot += vol_icv;

        input[icv] += input[icv_nbr]*vol_icv;
      }
      input[icv] /= vol_tot;
    }
  }

};



#endif
