#ifndef RANSTURBMODEL_ASBM_H
#define RANSTURBMODEL_ASBM_H

#include "TurbModel_KOMSST.h"
#include "TurbModel_V2F.h"


typedef RansTurbKOmSST TURB_MOD_FOR_ASBM; // define the turbulence model used with the ASBM

extern "C"{
  void asbm_(
      double*, double*, double*, double*, double*, double*, double*,
      double*, double*, double*, double*, int*, int*);
}

class RansTurbASBM : virtual public TURB_MOD_FOR_ASBM
{
public:   // constructors

  RansTurbASBM()
  {
    if (mpi_rank == 0)
      cout << "RansTurbASBM()" << endl;

    st_diag       = NULL;     registerVector(st_diag,       "st_diag",       CV_DATA);
    st_offdiag    = NULL;     registerVector(st_offdiag,    "st_offdiag",    CV_DATA);
    wt_offdiag    = NULL;     registerVector(wt_offdiag,    "wt_offdiag",    CV_DATA);

    as_diag       = NULL;     registerVector(as_diag,       "as_diag",       CV_DATA);
    as_offdiag    = NULL;     registerVector(as_offdiag,    "as_offdiag",    CV_DATA);
    ar_diag       = NULL;     registerVector(ar_diag,       "ar_diag",       CV_DATA);
    ar_offdiag    = NULL;     registerVector(ar_offdiag,    "ar_offdiag",    CV_DATA);

    etar          = NULL;     registerScalar(etar,          "etar",       CV_DATA);
    etaf          = NULL;     registerScalar(etaf,          "etaf",       CV_DATA);
    a2            = NULL;     registerScalar(a2,            "a2",         CV_DATA);
  }

  virtual ~RansTurbASBM() {}

protected:   // member vars

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

  FILE *finfo;              // output file for ASBM quantities

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    TURB_MOD_FOR_ASBM::initialHookScalarRansTurbModel();

    if (mpi_rank == 0)
      cout << "initialHook ASBM" << endl;

    start_asbm = getIntParam("START_ASBM", "0");

    if (mpi_rank == 0)
    {
      if ( (finfo=fopen("asbmInfo.dat","wt")) == NULL )
      {
        cerr << "Error: cannot open file asbmInfo.dat" << endl;
        throw(-1);
      }
    }
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) fclose(finfo);
  }

  virtual void calcRansTurbViscMuet()
  {
    TURB_MOD_FOR_ASBM::calcRansTurbViscMuet();

    // Non-linear domain partitioning
    if (step == start_asbm) setNonLinearDomain();

    // Computation of ASBM variables
    if (step >= start_asbm){
        //if ( step%block_frq == 0 ) calcBlockTensor();
        calcRsCenterASBM();
    }
  }

  virtual void setNonLinearDomain(){
    // default
    cout << "non-linear domain set as default" << endl;
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
    double ST[3][3], WT[3][3], BX[3][3];
    double EMPTY[3][3] = {{0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0}};

    // In and out variables
    double AS[3][3] = {{0.3, 0.0, 0.0},
                       {0.0, 0.3, 0.0},
                       {0.0, 0.0, 0.4}};
    double AR[3][3] = {{0.3, 0.0, 0.0},
                       {0.0, 0.3, 0.0},
                       {0.0, 0.0, 0.4}};

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
        //BB[0][0] = block_diag[icv][0];  BB[0][1] = block_offdiag[icv][0];  BB[0][2] = block_offdiag[icv][1];
        //BB[1][0] = BB[0][1];            BB[1][1] = block_diag[icv][1];     BB[1][2] = block_offdiag[icv][2];
        //BB[2][0] = BB[0][2];            BB[2][1] = BB[1][2];               BB[2][2] = block_diag[icv][2];

        // Call the ASBM
        ierr = 0;

        asbm_(
            &(ST[0][0]), &(WT[0][0]), &(WFT[0][0]), &(EMPTY[0][0]), &(EMPTY[0][0]), &(EMPTY[0][0]),
            &(AS[0][0]), &(AR[0][0]), &(REY[0][0]), &(DIM[0][0]), &(CIR[0][0]), &maxitrs, &ierr);
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

        etar[icv] = CIR[1][0];
        etaf[icv] = CIR[1][1];
        a2[icv]   = CIR[1][2];

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
          cerr << "Error: cannot open file asbmInfo.dat" << endl;
          throw(-1);
        }
        else
        {
          fprintf(finfo,"%8d\t%8d\n",step,Nerr);
        }
    }
  }

  void smoothingVec(double (*input)[3]){
    double vol_tot, vol_icv;
    for (int icv = 0; icv < ncv; icv++){
        vol_tot = cv_volume[icv];

        input[icv][0] *= vol_tot;
        input[icv][1] *= vol_tot;
        input[icv][2] *= vol_tot;

        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;

        for (int noc = noc_f + 1; noc <= noc_l; noc++){
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
     for (int icv = 0; icv < ncv; icv++){
         vol_tot = cv_volume[icv];

         input[icv] *= vol_tot;

         int noc_f = nbocv_i[icv];
         int noc_l = nbocv_i[icv + 1] - 1;

         for (int noc = noc_f + 1; noc <= noc_l; noc++){
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
