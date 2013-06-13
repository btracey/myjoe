#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"
#include "turbModels/TurbModel_ASBM.h"
#include "turbModels/TransModel_GaReT.h"


#include "combModels/CombModel_BinaryMixing.h"
#include "combModels/CombModel_VariableProperties.h"
#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_FPVA.h"
#include "combModels/CombModel_FPVA_Coeff.h"




// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels
{
public:

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;
  }

  virtual ~MyJoe()  {}

  void initialHook() { JoeWithModels::initialHook();}

  
  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == name)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "%s.dat", name.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"rho\" \"press\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(fname,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone != faZoneList.end(); faZone++)  // loop over all boundary zones
      if (faZone->getKind() == FA_ZONE_BOUNDARY)
        if (faZone->getNameString() == name)                            // identify zone = name
        {
          for (int ifa = faZone->ifa_f; ifa <= faZone->ifa_l; ifa++)    // loop over boundary faces
          {
            int icv = cvofa[ifa][0];                                    // get cv of boundary face

            double n[3], s_half[3], vel[3], velTang[3];

            normVec3d(n, fa_normal[ifa]);                               // get area weighted face normal (outward pointing)
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

            double wallTemp = 300.0;
            double mulam = mul_fa[ifa];

            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv], temp[icv], tau, qDot, yplus, mul_fa[ifa]);
          }
        }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }

  /*
   * write all values, WARNING: should be used only for small grids.
   */
  virtual void writeAllValues()
  {
    FILE *fp;
    if ( mpi_rank ==0)
    {
      if ( (fp = fopen("flowfield.txt", "wt")) == NULL)
      {
        cout << "Error: cannot open file flowfield.txt" << endl;
        throw(-1);
      }
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,123,mpi_comm,&status);
      if ( (fp = fopen("flowfield.txt", "a")) == NULL)
      {
        cout << "Error: cannot open file flowfield.txt" << endl;
        throw(-1);
      }
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      fprintf(fp, "%.12le\t%.12le\t%.12le\t%.12le\t%.12le\t%.12le\t",
          x_cv[icv][0], x_cv[icv][1], rho[icv], rhou[icv][0], rhou[icv][1],press[icv]);

      for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
        fprintf(fp, "%.12le\t", data->phi[icv]);

      fprintf(fp, "%.12le\t", InterpolateAtCellCenterFromFaceValues(mul_fa, icv));

      fprintf(fp, "%.12le\t", grad_u[icv][0][1]);

      double *muT = getR1("muT");
      fprintf(fp, "%.12le\t", muT[icv]);

      // Reynolds stresses
      fprintf(fp, "%.12le\t", rij_diag[icv][0]);
      fprintf(fp, "%.12le\t", rij_diag[icv][1]);
      fprintf(fp, "%.12le\t", rij_diag[icv][2]);
      fprintf(fp, "%.12le\t", rij_offdiag[icv][0]);

      fprintf(fp, "\n");
    }
    fclose(fp);

    if ( mpi_rank < mpi_size-1)
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,123,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }

  void transformMeshHook()
  {
    if (!checkParam("TRANSFORM_MESH")) return;

    double scl_mesh = getDoubleParam("SCL_MESH", "1.0");
    for (int ino = 0; ino < getNno(); ++ino)
      x_no[ino][0] *= scl_mesh;

    // clearFlag for wall distance to recompute the wallDist when mesh deformed
    DoubleScalar *wallD = getScalarData("wallDist");
    wallD->clearFlag();
  }

  virtual void temporalHook()
  {
    if (step%10000 == 0)
    {
      if (checkParam("WRITE_ALLVALUES"))
        writeAllValues();

      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY)
        {
          Param *param;
          if (getParam(param, zone->getName()))
            if (param->getString() == "WALL")
              writeWallValues(zone->getName());
        }
    }
  }

  virtual void finalHook()
  {
    if (checkParam("WRITE_ALLVALUES"))
      writeAllValues();

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            writeWallValues(zone->getName());
      }
  }
  
/*
 *  Following virtual functions can be overloaded by user.
 *  Comment: some functions require that a scalar is defined.
 */

//  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {/*empty*/}  
//  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {/*empty*/}
//
//  virtual void UserDefinedScalarClipping(const string &name)  {/*empty*/}
//
//  virtual void initialHookScalarRansTurbModel() {/*empty*/}
//  virtual void calcRansTurbViscMuet() {/*empty*/}
//  virtual double calcTurbProd(int icv) {/*empty*/}
//  virtual double calcDissipProd(int icv) {/*empty*/}
//  virtual void diffusivityHookScalarRansTurb(const string &name) {/*empty*/}
//  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}
//  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name) {/*empty*/}
//  virtual void sourceHookScalarRansTurbExpl(double *rhs, const string &name) {/*empty*/}
//
//  virtual void initialHookScalarRansCombModel() {/*empty*/} 
//  virtual void sourceHookRansComb(double *rhsRho, double(*rhsRhou)[3], double *rhsRhoE, double(*A)[5][5]) {/*empty*/}
//  virtual void sourceHookScalarRansComb(double *rhs, double *A, const string &name) {/*empty*/}
//  virtual void boundaryHookScalarRansComb(double *phi_ph, FaZone *zone, const string &name) {/*empty*/}  
//  virtual void UserDefinedScalarClipping(const string &name) {/*empty*/}
  

  
  
// ##################################################################
// possible other constructors 
// ##################################################################
  
//  MyJoe(ParamMap &p, int i) : JoeWithModels(p), UgpWithCvCompFlow(p), ownID(i) {
//    if (mpi_rank == 0)      cout << "MyJoe()" << endl;
//    init();
//  }

//  MyJoe(char *name, int i) : JoeWithModels(name), UgpWithCvCompFlow(name), ownID(i) {
//    if (mpi_rank == 0)
//      cout << "MyJoe()" << endl;
//    init();
//  }
  
};


/*
 * MyJoe with Spalart & Allmaras Model
 */
class MyJoeSA : public MyJoe, public RansTurbSA
{
public:
  MyJoeSA(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSA()" << endl;
  }

  virtual ~MyJoeSA() {}
};

/*
 * MyJoe with Menter SST Model
 */
class MyJoeSST : public MyJoe, public RansTurbKOmSST
{
public:
  MyJoeSST(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSST()" << endl;
  }

  virtual ~MyJoeSST() {}
};

/*
 * MyJoe with Wilcox k-omega Model
 */
class MyJoeWX : public MyJoe, public RansTurbKOm
{
public:
  MyJoeWX(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWX()" << endl;
  }

  virtual ~MyJoeWX() {}
};

/*
 * MyJoe with k-epsilon Model
 */
class MyJoeKEps : public MyJoe, public RansTurbKEps
{
public:
  MyJoeKEps(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeKEps()" << endl;
  }

  virtual ~MyJoeKEps() {}
};

/*
 * MyJoe with Durbin v2-f Model
 */
class MyJoeV2F : public MyJoe, public RansTurbV2F
{
public:
  MyJoeV2F(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeV2F()" << endl;
  }

  virtual ~MyJoeV2F() {}
};

/*
 * MyJoe with ASBM model
 */
class MyJoeASBM : public MyJoe, public RansTurbASBM
{
public:
  MyJoeASBM(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeASBM()" << endl;
  }

  virtual ~MyJoeASBM() {}
};

/*
 * Compute the Barycentric maps
 */
class BaryMaps: public MyJoeASBM{
public:
  double (*baryCo)[3];          /// Barycentric map coordinates
  double (*colorC)[3];          /// Barycentric map colors

public:
  BaryMaps(char *name) : MyJoeASBM(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "BaryMaps()" << endl;

    baryCo = NULL;  registerVector(baryCo,  "baryCo",       CV_DATA);
    colorC = NULL;  registerVector(colorC,  "colorC",       CV_DATA);
  }

  virtual ~BaryMaps()  {}

  void perturbStress()
  {
    for (int icv=0; icv<ncv; icv++)
    {
      //####################################################################################
      //for every cell compute anisotropy tensor to get barycentric coords
      double aij[3][3],eigv[3][3],eigvt[3][3],eigs[3],eigsnew[3][3];
      int order[3];

      // open variables for eigenvalue and anisotropy tensor
      for (int i=0; i < 3; i++)
      {
        for (int j=0; j < 3; j++)
        {
          aij[i][j] = 0.0;
          eigv[i][j] = 0.0;
          eigvt[i][j] = 0.0;
          eigsnew[i][j] = 0.0;
        }
        eigs[i] = 0.0;
        order[i] = 0;
      }

      // anisotropy tensor (unordered)
//      double temp = 2.0*rho[icv]*kine[icv];
//      aij[0][0] = -1.0*rij_diag[icv][0]/temp - 1.0/3.0;
//      aij[1][1] = -1.0*rij_diag[icv][1]/temp - 1.0/3.0;
//      aij[2][2] = -1.0*rij_diag[icv][2]/temp - 1.0/3.0;
//      aij[0][1] = -1.0*rij_offdiag[icv][0]/temp;
//      aij[0][2] = -1.0*rij_offdiag[icv][1]/temp;
//      aij[1][2] = -1.0*rij_offdiag[icv][2]/temp;
//      aij[1][0] = aij[0][1];
//      aij[2][0] = aij[0][2];
//      aij[2][1] = aij[1][2];

      aij[0][0] = rij_diag_nd[icv][0] - 1.0/3.0;
      aij[1][1] = rij_diag_nd[icv][1] - 1.0/3.0;
      aij[2][2] = rij_diag_nd[icv][2] - 1.0/3.0;
      aij[0][1] = rij_offdiag_nd[icv][0];
      aij[0][2] = rij_offdiag_nd[icv][1];
      aij[1][2] = rij_offdiag_nd[icv][2];
      aij[1][0] = aij[0][1];
      aij[2][0] = aij[0][2];
      aij[2][1] = aij[1][2];

      //aij tensor now built - build barycentric map locations
      //compute eigenvalues of aij
      eigen_decomposition(aij,eigv,eigs);
      //generate transpose
      for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++)
          eigvt[j][i] = eigv[i][j];

      //sort eigenvalues by magnitude, save ordering in 'order'
      if ((eigs[0] >= eigs[1]) && (eigs[0] >= eigs[2]))
      {
        order[0] = 0;
        order[1] = 1;
        order[2] = 2;
        if (eigs[2] >= eigs[1])
        {
          order[1] = 2;
          order[2] = 1;
        }
      }
      if ((eigs[1] >= eigs[0]) && (eigs[1] >= eigs[2]))
      {
        order[0] = 1;
        order[1] = 0;
        order[2] = 2;
        if (eigs[2] >= eigs[0])
        {
          order[1] = 2;
          order[2] = 0;
        }
      }
      if ((eigs[2] >= eigs[0]) && (eigs[2] >= eigs[1]))
      {
        order[0] = 2;
        order[1] = 0;
        order[2] = 1;
        if (eigs[1] >= eigs[0])
        {
          order[1] = 1;
          order[2] = 0;
        }
      }

      //####################################################################################
      //Generate barycentric map coordinates
      double c1c = eigs[order[0]] - eigs[order[1]];
      double c2c = 2.0 * (eigs[order[1]] - eigs[order[2]]);
      double c3c = 3.0 * eigs[order[2]] + 1.0;

      //save barycentric coordinates
      baryCo[icv][0] = 1.0 * c1c + 0.0 * c2c + 0.5 * c3c;
      baryCo[icv][1] = 0.0 * c1c + 0.0 * c2c + 0.866025 * c3c;
      baryCo[icv][2] = 0.0;

      colorC[icv][0] = c1c;
      colorC[icv][1] = c2c;
      colorC[icv][2] = c3c;
    }
    updateCvData(baryCo, REPLACE_ROTATE_DATA);
    updateCvData(colorC, REPLACE_ROTATE_DATA);
  }

  virtual void temporalHook()
  {
    if (step%10000 == 0)
    {
      perturbStress();

      if (checkParam("WRITE_ALLVALUES"))
        writeAllValues();

      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY)
        {
          Param *param;
          if (getParam(param, zone->getName()))
            if (param->getString() == "WALL")
              writeWallValues(zone->getName());
        }
    }
  }

  virtual void finalHook()
  {
    perturbStress();

    if (checkParam("WRITE_ALLVALUES"))
      writeAllValues();

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            writeWallValues(zone->getName());
      }
  }

  // ###################################################################
  // The files below are used for computing eigenvalues and eigenvectors
  // ###################################################################

  void matMult(double (*res)[3], double (*m1)[3], double (*m2)[3])  // m1*m2 matrices
   {
     for (int i=0; i<3; i++)
       for (int j=0; j<3; j++)
         {
           res[i][j] = 0.0;
           for (int k=0; k<3; k++)
             res[i][j] += m1[i][k]*m2[k][j];
         }
   }

   //these routiens are for computing the eigenvalue decomposition
   static double hypot2(double x, double y)
   {
     return sqrt(x*x+y*y);
   }

   static void tred2(double V[3][3], double d[3], double e[3]) {

     //  This is derived from the Algol procedures tred2 by
     //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
     //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
     //  Fortran subroutine in EISPACK.

     for (int j = 0; j < 3; j++) {
       d[j] = V[3-1][j];
     }

     // Householder reduction to tridiagonal form.

     for (int i = 3-1; i > 0; i--) {

       // Scale to avoid under/overflow.

       double scale = 0.0;
       double h = 0.0;
       for (int k = 0; k < i; k++) {
         scale = scale + fabs(d[k]);
       }
       if (scale == 0.0) {
         e[i] = d[i-1];
         for (int j = 0; j < i; j++) {
           d[j] = V[i-1][j];
           V[i][j] = 0.0;
           V[j][i] = 0.0;
         }
       } else {

         // Generate Householder vector.

         for (int k = 0; k < i; k++) {
           d[k] /= scale;
           h += d[k] * d[k];
         }
         double f = d[i-1];
         double g = sqrt(h);
         if (f > 0) {
           g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (int j = 0; j < i; j++) {
           e[j] = 0.0;
         }

         // Apply similarity transformation to remaining columns.

         for (int j = 0; j < i; j++) {
           f = d[j];
           V[j][i] = f;
           g = e[j] + V[j][j] * f;
           for (int k = j+1; k <= i-1; k++) {
             g += V[k][j] * d[k];
             e[k] += V[k][j] * f;
           }
           e[j] = g;
         }
         f = 0.0;
         for (int j = 0; j < i; j++) {
           e[j] /= h;
           f += e[j] * d[j];
         }
         double hh = f / (h + h);
         for (int j = 0; j < i; j++) {
           e[j] -= hh * d[j];
         }
         for (int j = 0; j < i; j++) {
           f = d[j];
           g = e[j];
           for (int k = j; k <= i-1; k++) {
             V[k][j] -= (f * e[k] + g * d[k]);
           }
           d[j] = V[i-1][j];
           V[i][j] = 0.0;
         }
       }
       d[i] = h;
     }

     // Accumulate transformations.

     for (int i = 0; i < 3-1; i++) {
       V[3-1][i] = V[i][i];
       V[i][i] = 1.0;
       double h = d[i+1];
       if (h != 0.0) {
         for (int k = 0; k <= i; k++) {
           d[k] = V[k][i+1] / h;
         }
         for (int j = 0; j <= i; j++) {
           double g = 0.0;
           for (int k = 0; k <= i; k++) {
             g += V[k][i+1] * V[k][j];
           }
           for (int k = 0; k <= i; k++) {
             V[k][j] -= g * d[k];
           }
         }
       }
       for (int k = 0; k <= i; k++) {
         V[k][i+1] = 0.0;
       }
     }
     for (int j = 0; j < 3; j++) {
       d[j] = V[3-1][j];
       V[3-1][j] = 0.0;
     }
     V[3-1][3-1] = 1.0;
     e[0] = 0.0;
   }
         // Symmetric tridiagonal QL algorithm.

   static void tql2(double V[3][3], double d[3], double e[3]) {

                 //  This is derived from the Algol procedures tql2, by
                 //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
                 //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
                 //  Fortran subroutine in EISPACK.

                 for (int i = 1; i < 3; i++) {
                         e[i-1] = e[i];
                 }
                 e[3-1] = 0.0;

                 double f = 0.0;
                 double tst1 = 0.0;
                 double eps = pow(2.0,-52.0);
                 for (int l = 0; l < 3; l++) {

                         // Find small subdiagonal element

                         tst1 = max(tst1,(fabs(d[l]) + fabs(e[l])));
                         int m = l;
                         while (m < 3) {
                                 if (fabs(e[m]) <= eps*tst1) {
                                         break;
                                 }
                                 m++;
                         }

                         // If m == l, d[l] is an eigenvalue,
                         // otherwise, iterate.

                         if (m > l) {
                                 int iter = 0;
                                 do {
                                         iter = iter + 1;  // (Could check iteration count here.)

                                         // Compute implicit shift

                                         double g = d[l];
                                         double p = (d[l+1] - g) / (2.0 * e[l]);
                                         double r = hypot2(p,1.0);
                                         if (p < 0) {
                                                 r = -r;
                                         }
                                         d[l] = e[l] / (p + r);
                                         d[l+1] = e[l] * (p + r);
                                         double dl1 = d[l+1];
                                         double h = g - d[l];
                                         for (int i = l+2; i < 3; i++) {
                                                 d[i] -= h;
                                         }
                                         f = f + h;

                                         // Implicit QL transformation.

                                         p = d[m];
                                         double c = 1.0;
                                         double c2 = c;
                                         double c3 = c;
                                         double el1 = e[l+1];
                                         double s = 0.0;
                                         double s2 = 0.0;
                                         for (int i = m-1; i >= l; i--) {
                                                 c3 = c2;
                                                 c2 = c;
                                                 s2 = s;
                                                 g = c * e[i];
                                                 h = c * p;
                                                 r = hypot2(p,e[i]);
                                                 e[i+1] = s * r;
                                                 s = e[i] / r;
                                                 c = p / r;
                                                 p = c * d[i] - s * g;
                                                 d[i+1] = h + s * (c * g + s * d[i]);

                                                 // Accumulate transformation.

                                                 for (int k = 0; k < 3; k++) {
                                                         h = V[k][i+1];
                                                         V[k][i+1] = s * V[k][i] + c * h;
                                                         V[k][i] = c * V[k][i] - s * h;
                                                 }
                                         }
                                         p = -s * s2 * c3 * el1 * e[l] / dl1;
                                         e[l] = s * p;
                                         d[l] = c * p;

                                         // Check for convergence.

                                 } while (fabs(e[l]) > eps*tst1);
                         }
                         d[l] = d[l] + f;
                         e[l] = 0.0;
                 }

                 // Sort eigenvalues and corresponding vectors.

                 for (int i = 0; i < 3-1; i++) {
                         int k = i;
                         double p = d[i];
                         for (int j = i+1; j < 3; j++) {
                                 if (d[j] < p) {
                                         k = j;
                                         p = d[j];
                                 }
                         }
                         if (k != i) {
                                 d[k] = d[i];
                                 d[i] = p;
                                 for (int j = 0; j < 3; j++) {
                                         p = V[j][i];
                                         V[j][i] = V[j][k];
                                         V[j][k] = p;
                                 }
                         }
                 }
         }

   void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
                 double e[3];
                 for (int i = 0; i < 3; i++) {
                         for (int j = 0; j < 3; j++) {
                                 V[i][j] = A[i][j];
                         }
                 }
                 tred2(V, d, e);
                 tql2(V, d, e);
         }
};

/*
 * Flat channel with inlet and outlet bcs, SST
 */
class NonPerChanSST: public MyJoeSST{
protected:

public:
  NonPerChanSST(char *name) : MyJoeSST(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "NonPerChanSST()" << endl;
  }

  virtual ~NonPerChanSST() {}

  void readReynoldsStress()
  {
    int    nx;             // number of nodes in streamwise direction
    int    ny;             // number of nodes in cross-stream direction
    double ***reStress;     // holder for input profile data

    /*** Read the stresses ***/
    FILE *ifile;
    if ((ifile=fopen("./reStress.txt", "rt")) == NULL)
      {
        cout << "could not open reStress.txt" << endl;
        throw(-1);
      }

    // file values: uu vv ww uv
    fscanf(ifile, "nx=%d\tny=%d", &nx, &ny);
    reStress = new double** [nx];
    for (int i = 0; i < nx; i++){
      reStress[i] = new double* [ny];
      for (int j = 0; j < ny; j++)
        reStress[i][j] = new double [6];
    }

    for (int i = 0; i<nx; i++)
      for (int j = 0; j < ny; j++)
        for (int v=0; v<6; v++)
          fscanf(ifile, "%lf", &reStress[i][j][v]);

    fclose(ifile);

    // Interpolate
    double re11, re12, re21, re22;
    double rey1, rey2;

    for (int icv=0; icv<ncv; icv++)
    {
      int posx = 1, posy = 1;
      double fx, fy;

      // while posx and posx-1 don't sandwich x_cv, keep increasing posx
      while(reStress[posx][0][0] < x_cv[icv][0] && (posx < nx-1))
        posx++;
      // if profilesIC doesn't have a node high enough to sandwich x_cv[icv]
      if (x_cv[icv][0] > reStress[posx][0][0])
        fx = 1.0;
      // if profilesIC doesn't have a node low enough to sandwich x_cv[icv]
      else if (x_cv[icv][0] < reStress[posx-1][0][0])
        fx = 0.0;
      else
        fx = (x_cv[icv][0] - reStress[posx-1][0][0])/(reStress[posx][0][0] - reStress[posx-1][0][0]);

      // while posy and posy-1 don't sandwich x_cv, keep increasing posy
      while(reStress[posx][posy][1] < x_cv[icv][1] && (posy < ny-1))
        posy++;
      // if profilesIC doesn't have a node high enough to sandwich x_cv[icv]
      if (x_cv[icv][1] > reStress[posx][posy][1])
        fy = 1.0;
      // if profilesIC doesn't have a node low enough to sandwich x_cv[icv]
      else if (x_cv[icv][1] < reStress[posx][posy-1][1])
        fy = 0.0;
      else
        fy = (x_cv[icv][1] - reStress[posx][posy-1][1])/(reStress[posx][posy][1] - reStress[posx][posy-1][1]);

      re11 = reStress[posx-1][posy-1][0];
      re12 = reStress[posx][posy-1][0];
      re21 = reStress[posx-1][posy][0];
      re22 = reStress[posx][posy][0];
      rey1 = re11 + fx*(re12 - re11);
      rey2 = re21 + fx*(re22 - re21);
      rij_diag[icv][0] = rey1 + fy*(rey2 - rey1);

      re11 = reStress[posx-1][posy-1][1];
      re12 = reStress[posx][posy-1][1];
      re21 = reStress[posx-1][posy][1];
      re22 = reStress[posx][posy][1];
      rey1 = re11 + fx*(re12 - re11);
      rey2 = re21 + fx*(re22 - re21);
      rij_diag[icv][1] = rey1 + fy*(rey2 - rey1);

      re11 = reStress[posx-1][posy-1][2];
      re12 = reStress[posx][posy-1][2];
      re21 = reStress[posx-1][posy][2];
      re22 = reStress[posx][posy][2];
      rey1 = re11 + fx*(re12 - re11);
      rey2 = re21 + fx*(re22 - re21);
      rij_diag[icv][2] = rey1 + fy*(rey2 - rey1);

      re11 = reStress[posx-1][posy-1][3];
      re12 = reStress[posx][posy-1][3];
      re21 = reStress[posx-1][posy][3];
      re22 = reStress[posx][posy][3];
      rey1 = re11 + fx*(re12 - re11);
      rey2 = re21 + fx*(re22 - re21);
      rij_offdiag[icv][0] = rey1 + fy*(rey2 - rey1);
      rij_offdiag[icv][1] = 0.0;
      rij_offdiag[icv][2] = 0.0;
    }

    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);
  }

};

/*
 * Main Main Main Main Main Main
 */
int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of JoeWithModels
  int run = 0;

  // set input name to default "Joe.in"
  char inputFileName[50];
  sprintf(inputFileName, "Joe.in");

  for (int i=1; i<argc; i++)
  {
    string str(argv[i]);
    if (from_string<int>(run, str, std::dec))
    {
      if (mpi_rank == 0)
        cout << "You have specified run number = " << run << endl;
    }
    else
      strcpy(inputFileName, argv[i]);
  }

  if (mpi_rank == 0)
  {
    cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
    cout << "SPECIFIED RUN = " << run << endl;
  }


  try {

    // declare pointer to JoeWithModels
    JoeWithModels *joe;
    
    switch (run)
    {
    case 0:   joe = new MyJoe(inputFileName);           break;
    case 1:   joe = new MyJoeSA(inputFileName);         break;
    case 2:   joe = new MyJoeSST(inputFileName);        break;
    case 3:   joe = new MyJoeWX(inputFileName);         break;
    case 4:   joe = new MyJoeKEps(inputFileName);       break;
    case 5:   joe = new MyJoeV2F(inputFileName);        break;
    case 6:   joe = new MyJoeASBM(inputFileName);       break;
    case 7:   joe = new BaryMaps(inputFileName);        break;
    case 8:   joe = new NonPerChanSST(inputFileName);   break;
    default: 
      if (mpi_rank == 0)
        cerr << "ERROR: run number not available!" << endl;
      throw(-1);
    }
    
    // provide total runtime 
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    // run joe
    joe->run();
    

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

    // delete joe (make sure memory is deallocated in destructors
    delete joe;
  }
  catch (int e) {
    cerr << "Exception: " << e << endl;
    MPI_Finalize();
    return(-1);
  }
  catch(...) {
    cerr << "unhandled Exception.\n" << endl;
    MPI_Finalize();
    return(-1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return (0);

}



