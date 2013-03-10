#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"
#include "turbModels/TurbModel_V2F_half.h"
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
  void writeAllValues()
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
      fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t",
          x_cv[icv][0], x_cv[icv][1], rho[icv], rhou[icv][0], rhou[icv][1],press[icv]);

      for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
        fprintf(fp, "%.8le\t", data->phi[icv]);

      fprintf(fp, "%.8le\t", InterpolateAtCellCenterFromFaceValues(mul_fa, icv));

      fprintf(fp, "%.8le\t", grad_u[icv][0][1]);

      double *muT = getR1("muT");
      fprintf(fp, "%.8le\t", muT[icv]);

      // Reynolds stresses
      fprintf(fp, "%.8le\t", rij_diag[icv][0]);
      fprintf(fp, "%.8le\t", rij_diag[icv][1]);
      fprintf(fp, "%.8le\t", rij_diag[icv][2]);
      fprintf(fp, "%.8le\t", rij_offdiag[icv][0]);

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
    if (step%1000 == 0)
    {
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
 * MyJoe with half of the v2-f Model
 */
class MyJoeV2F_half : public MyJoe, public RansTurbV2F_half
{
public:
  MyJoeV2F_half(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeV2F_half()" << endl;
  }

  virtual ~MyJoeV2F_half() {}
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
 * Flat channel with periodic bc's, SST
 */
class PerChanSST: public MyJoeSST{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  PerChanSST(char *name) : MyJoeSST(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "PerChanSST()" << endl;
    boundVal = NULL;
  }

  virtual ~PerChanSST()
  {
    if (boundVal != NULL) delete []boundVal;
  }

  void initialHook()
  {
    JoeWithModels::initialHook();

    if (checkParam("SET_INIT_PROFILE"))
    {
      // Read inlet variable profile
      // file has variables y, rho, u, v, press, ...
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
      {
        cout << "could not open profiles.dat, apply boundary from input file" << endl;
        throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      for (int i=0; i<nn; i++)
        for (int v = 0; v < nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);

      // Specify initial condition over whole flow

      if(!checkDataFlag(rho))
      {
        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double f;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    f = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  f = 0.0;
          else    f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+f*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+f*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+f*(press2-press1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv])
              + rho[icv]*kine[icv];

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+f*(kine2-kine1);

          double omega1 = boundVal[pos-1][6];
          double omega2 = boundVal[pos][6];
          omega[icv] = omega1+f*(omega2-omega1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(omega, REPLACE_DATA);
      }
    }
  }

  void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5])
  {
    double grav  = getDoubleParam("grav",  "1.00");
    for (int icv = 0; icv < ncv; icv++)
    {
      rhs_rhou[icv][0] += rho[icv]*cv_volume[icv]*grav;
      rhs_rhoE[icv] += rhou[icv][0]*cv_volume[icv]*grav;

      if (A != NULL){
        A[nbocv_i[icv]][1][0] -= cv_volume[icv]*grav;
        A[nbocv_i[icv]][4][1] -= cv_volume[icv]*grav;
      }
    }
  }

  virtual void finalHook()
  {
    // Extract all data
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
};

/*
 * Flat channel with periodic bc's, V2F
 */
class PerChanV2F: public MyJoeV2F{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  PerChanV2F(char *name) : MyJoeV2F(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "PerChanV2F()" << endl;
    boundVal = NULL;
  }

  virtual ~PerChanV2F()
  {
    if (boundVal != NULL) delete []boundVal;
  }

  void initialHook()
  {
    JoeWithModels::initialHook();

    if (checkParam("SET_INIT_PROFILE"))
    {
      // Read inlet variable profile
      // file has variables y, rho, u, v, press, ...
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
      {
        cout << "could not open profiles.dat, apply boundary from input file" << endl;
          throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      for (int i=0; i<nn; i++)
        for (int v=0; v<nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);

      // Specify initial condition over whole flow

      if(!checkDataFlag(rho))
      {
        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double fi;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    fi = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  fi = 0.0;
          else    fi = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+fi*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+fi*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+fi*(press2-press1);

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+fi*(kine2-kine1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv]) +
              rho[icv]*kine[icv];

          double eps1 = boundVal[pos-1][6];
          double eps2 = boundVal[pos][6];
          eps[icv] = eps1+fi*(eps2-eps1);

          double v21 = boundVal[pos-1][7];
          double v22 = boundVal[pos][7];
          v2[icv] = v21 + fi*(v22 - v21);

          double f1 = boundVal[pos-1][8];
          double f2 = boundVal[pos][8];
          f[icv] = f1 + fi*(f2 - f1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(eps, REPLACE_DATA);
        updateCvData(v2, REPLACE_DATA);
        updateCvData(f, REPLACE_DATA);
      }
    }
  }

  void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]){
    double grav  = getDoubleParam("grav",  "1.00");

    for (int icv = 0; icv < ncv; icv++){
        rhs_rhou[icv][0] += rho[icv]*cv_volume[icv]*grav;
        rhs_rhoE[icv] += rhou[icv][0]*cv_volume[icv]*grav;

        if (A != NULL){
            A[nbocv_i[icv]][1][0] -= cv_volume[icv]*grav;
            A[nbocv_i[icv]][4][1] -= cv_volume[icv]*grav;
        }
    }
  }

  virtual void finalHook()
  {
    // Extract all data
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
};

/*
 * Flat channel with periodic bc's, KEps
 */
class PerChanKEps: public MyJoeKEps{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  PerChanKEps(char *name) : MyJoeKEps(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "PerChanKEps()" << endl;
    boundVal = NULL;
  }

  virtual ~PerChanKEps()
  {
    if (boundVal != NULL) delete []boundVal;
  }

  void initialHook()
  {
    JoeWithModels::initialHook();

    if (checkParam("SET_INIT_PROFILE"))
    {
      // Read inlet variable profile
      // file has variables y, rho, u, v, press, ...
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
      {
        cout << "could not open profiles.dat, apply boundary from input file" << endl;
          throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      for (int i=0; i<nn; i++)
        for (int v=0; v<nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);

      // Specify initial condition over whole flow

      if(!checkDataFlag(rho))
      {
        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double fi;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    fi = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  fi = 0.0;
          else    fi = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+fi*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+fi*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+fi*(press2-press1);

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+fi*(kine2-kine1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv]) +
              rho[icv]*kine[icv];

          double eps1 = boundVal[pos-1][6];
          double eps2 = boundVal[pos][6];
          eps[icv] = eps1+fi*(eps2-eps1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(eps, REPLACE_DATA);
      }
    }
  }

  void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]){
    double grav  = getDoubleParam("grav",  "1.00");
    for (int icv = 0; icv < ncv; icv++){
        rhs_rhou[icv][0] += rho[icv]*cv_volume[icv]*grav;
        rhs_rhoE[icv] += rhou[icv][0]*cv_volume[icv]*grav;

        if (A != NULL){
            A[nbocv_i[icv]][1][0] -= cv_volume[icv]*grav;
            A[nbocv_i[icv]][4][1] -= cv_volume[icv]*grav;
        }
    }
  }

  virtual void finalHook()
  {
    // Extract all data
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
};

/*
 * Boundary Layer on Flat Plate with SST
 */
class BlayerSST: public MyJoeSST {
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  BlayerSST(char *name) : MyJoeSST(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "BlayerSST()" << endl;
  }

  virtual ~BlayerSST()  {}

  void initialHook()
  {
    JoeWithModels::initialHook();

    boundVal = NULL;

    // Read inlet variable profile
    if (checkParam("READ_PROFILE"))
    {
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
      {
        cout << "could not open the file profiles.dat" << endl;
        throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      // file has variables y, rho, u, v, press, ...
      for (int i=0; i<nn; i++)
        for (int v = 0; v < nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);
      if (mpi_rank == 0)
        cout << "FINISHED READING INLET PROFILE" << endl;
    }

    // Specify initial condition over whole flow
    if (checkParam("SET_INIT_PROFILE"))
    {
      if(!checkDataFlag(rho))
      {
        if (boundVal == NULL)
        {
          if (mpi_rank == 0)
            cerr << "Have not read profiles.dat" << endl;
          throw(-1);
        }

        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double f;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    f = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  f = 0.0;
          else    f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+f*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+f*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+f*(press2-press1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv])
                      + rho[icv]*kine[icv];

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+f*(kine2-kine1);

          double omega1 = boundVal[pos-1][6];
          double omega2 = boundVal[pos][6];
          omega[icv] = omega1+f*(omega2-omega1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(omega, REPLACE_DATA);
      }
    }
  }

  virtual void boundaryHook(double *T_input, double (*vel_input)[3], double *p_input, FaZone *zone)
  {
    if (zone->getNameString() == getStringParam("INLET_NAME"))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }

      double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;
      gamma_bc = 1.4;
      double gm1 = gamma_bc - 1.0;
      double ovgm1 = 1.0/gm1;

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double uvel1 = boundVal[pos-1][2];
        double uvel2 = boundVal[pos][2];
        u_bc[0] = uvel1 + f*(uvel2 - uvel1);

        double vvel1 = boundVal[pos-1][3];
        double vvel2 = boundVal[pos][3];
        u_bc[1] = 0.0; //vvel1 + f*(vvel2 - vvel1);
        u_bc[2] = 0.0;

        double press1 = boundVal[pos-1][4];
        double press2 = boundVal[pos][4];
        p_bc = press1 + f*(press2-press1);

        double temp1 = press1/(boundVal[pos-1][1]*R_gas);
        double temp2 = press2/(boundVal[pos][1]*R_gas);
        T_bc = temp1 + f*(temp2-temp1);

        rho_bc = p_bc/(R_gas*T_bc);
        c_bc = sqrt(gamma_bc*p_bc/rho_bc);
        s_bc = pow(rho_bc,gamma_bc)/p_bc;

        int icv0 = cvofa[ifa][0];
        assert( icv0 >= 0 );

        double nVec[3];
        double area = normVec3d(nVec, fa_normal[ifa]);

        //Compute normal velocity of freestream
        double qn_bc = vecDotVec3d(nVec,u_bc);

        //Subtract from normal velocity of mesh
        double vn_bc = qn_bc;// - normalvel_bfa[ifa];

        //Compute normal velocity of internal cell
        double qne = vecDotVec3d(nVec,vel[icv0]);

        //Compute speed of sound in internal cell
        double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

        double ac1, ac2;

        //Compute the Riemann invariants in the halo cell
        if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
          ac1 = qne   + 2.0*ovgm1*ce;
        else                     // Supersonic inflow.
          ac1 = qn_bc + 2.0*ovgm1*c_bc;

        if(vn_bc > c_bc)             // Supersonic outflow.
          ac2 = qne   - 2.0*ovgm1*ce;
        else                     // Inflow or subsonic outflow.
          ac2 = qn_bc - 2.0*ovgm1*c_bc;


        double qnf = 0.5*(ac1 + ac2);
        double cf  = 0.25*(ac1 - ac2)*gm1;

        double velf[3], sf;

        if (vn_bc > 0)                                                       // Outflow
        {
          for (int i=0; i<3; i++)
            velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
          sf = pow( rho[icv0], gamma[icv0])/press[icv0];
        }
        else                                                                 // Inflow
        {
          for (int i=0; i<3; i++)
            velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
          sf = s_bc;
        }
        //Compute density, pressure and velocity at boundary face
        double rho_int = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
        for (int i=0; i<3; i++)
          vel_input[ifa][i] = velf[i];
        p_input[ifa] = rho_int*cf*cf/gamma[icv0];
        T_input[ifa] = p_input[ifa]/(R_gas*rho_int);
      }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbKOmSST::boundaryHookScalarRansTurb(phi_ph, zone, name);

    Param *param;

    if ((name == "kine") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + f*(kine2-kine1);
      }
    }

    if ((name == "omega") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double omega1 = boundVal[pos-1][6];
        double omega2 = boundVal[pos][6];
        phi_ph[ifa] =  (omega1 + f*(omega2-omega1));
      }
    }
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

  virtual void finalHook()
  {
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
};

/*
 * Boundary Layer on Flat Plate with V2F
 */
class BlayerV2F: public MyJoeV2F{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  BlayerV2F(char *name) : MyJoeV2F(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "BlayerV2F()" << endl;
  }

  virtual ~BlayerV2F()  {}

  void initialHook()
  {
    JoeWithModels::initialHook();

    boundVal = NULL;

    // Read inlet variable profile
    if (checkParam("READ_PROFILE"))
    {
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL) {
        cout << "could not open the file profiles.dat" << endl;
        throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      // file has variables y, rho, u, v, press, ...
      for (int i=0; i<nn; i++)
        for (int v=0; v<nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);
      if (mpi_rank == 0)
        cout << "FINISHED READING INLET PROFILE" << endl;
    }

    // Specify initial condition over whole flow
    if (checkParam("SET_INIT_PROFILE"))
    {
      if(!checkDataFlag(rho))
      {
        if (boundVal == NULL)
        {
          if (mpi_rank == 0)
            cerr << "Have not read profiles.dat" << endl;
          throw(-1);
        }

        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double fi;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    fi = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  fi = 0.0;
          else    fi = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+fi*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+fi*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+fi*(press2-press1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv])
                      + rho[icv]*kine[icv];

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+fi*(kine2-kine1);

          double eps1 = boundVal[pos-1][6];
          double eps2 = boundVal[pos][6];
          eps[icv] = eps1+fi*(eps2-eps1);

          double v21 = boundVal[pos-1][7];
          double v22 = boundVal[pos][7];
          v2[icv] = v21 + fi*(v22 - v21);

          double f1 = boundVal[pos-1][8];
          double f2 = boundVal[pos][8];
          f[icv] = f1 + fi*(f2 - f1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(eps, REPLACE_DATA);
        updateCvData(v2, REPLACE_DATA);
        updateCvData(f, REPLACE_DATA);
      }
    }
  }

  virtual void boundaryHook(double *T_input, double (*vel_input)[3], double *p_input, FaZone *zone)
  {
    if (zone->getNameString() == getStringParam("INLET_NAME"))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }

        double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;
        gamma_bc = 1.4;
        double gm1 = gamma_bc - 1.0;
        double ovgm1 = 1.0/gm1;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int pos=1;
            while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

            double fi;
            if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
            else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
            else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

            double uvel1 = boundVal[pos-1][2];
            double uvel2 = boundVal[pos][2];
            u_bc[0] = uvel1 + fi*(uvel2 - uvel1);

            double vvel1 = boundVal[pos-1][3];
            double vvel2 = boundVal[pos][3];
            u_bc[1] = 0.0; //vvel1 + f*(vvel2 - vvel1);
            u_bc[2] = 0.0;

            double press1 = boundVal[pos-1][4];
            double press2 = boundVal[pos][4];
            p_bc = press1 + fi*(press2-press1);

            double temp1 = press1/(boundVal[pos-1][1]*R_gas);
            double temp2 = press2/(boundVal[pos][1]*R_gas);
            T_bc = temp1 + fi*(temp2-temp1);

            rho_bc = p_bc/(R_gas*T_bc);
            c_bc = sqrt(gamma_bc*p_bc/rho_bc);
            s_bc = pow(rho_bc,gamma_bc)/p_bc;

            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            //Compute normal velocity of freestream
            double qn_bc = vecDotVec3d(nVec,u_bc);

            //Subtract from normal velocity of mesh
            double vn_bc = qn_bc;// - normalvel_bfa[ifa];

            //Compute normal velocity of internal cell
            double qne = vecDotVec3d(nVec,vel[icv0]);

            //Compute speed of sound in internal cell
            double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

            double ac1, ac2;

            //Compute the Riemann invariants in the halo cell
            if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
              ac1 = qne   + 2.0*ovgm1*ce;
            else                     // Supersonic inflow.
              ac1 = qn_bc + 2.0*ovgm1*c_bc;


            if(vn_bc > c_bc)             // Supersonic outflow.
              ac2 = qne   - 2.0*ovgm1*ce;
            else                     // Inflow or subsonic outflow.
              ac2 = qn_bc - 2.0*ovgm1*c_bc;


            double qnf = 0.5*(ac1 + ac2);
            double cf  = 0.25*(ac1 - ac2)*gm1;

            double velf[3], sf;

            if (vn_bc > 0)                                                       // Outflow
              {
                for (int i=0; i<3; i++)
                  velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
                sf = pow( rho[icv0], gamma[icv0])/press[icv0];
              }
            else                                                                 // Inflow
              {
                for (int i=0; i<3; i++)
                  velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
                sf = s_bc;
              }
            //Compute density, pressure and velocity at boundary face
            double rho_int = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
            for (int i=0; i<3; i++)
              vel_input[ifa][i] = velf[i];
            p_input[ifa] = rho_int*cf*cf/gamma[icv0];
            T_input[ifa] = p_input[ifa]/(R_gas*rho_int);
          }
      }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbV2F::boundaryHookScalarRansTurb(phi_ph, zone, name);

    Param *param;

    if ((name == "kine") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + fi*(kine2-kine1);
      }
    }

    if ((name == "eps") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double eps1 = boundVal[pos-1][6];
        double eps2 = boundVal[pos][6];
        phi_ph[ifa] =  (eps1 + fi*(eps2-eps1));
      }
    }
    if ((name == "v2") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double v21 = boundVal[pos-1][7];
        double v22 = boundVal[pos][7];
        phi_ph[ifa] =  (v21 + fi*(v22-v21));
      }
    }
    if ((name == "f") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double f1 = boundVal[pos-1][8];
        double f2 = boundVal[pos][8];
        phi_ph[ifa] =  (f1 + fi*(f2-f1));
      }
    }
  }

  void transformMeshHook()
  {
    if (!checkParam("TRANSFORM_MESH")) return;

    double scl_mesh = getDoubleParam("SCL_MESH", "1.0");
    for (int ino = 0; ino < getNno(); ++ino)
      x_no[ino][0] *= scl_mesh;

    // clearFlag for wall distance to recompute the wallDist when mesh deformed
    //DoubleScalar *wallD = getScalarData("wallDist");
    //wallD->clearFlag();
  }

  virtual void temporalHook()
  {
    if (step%10000 == 0)
    {
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
};

/*
 * Boundary Layer on Flat Plate with KEps
 */
class BlayerKEps: public MyJoeASBM{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  BlayerKEps(char *name) : MyJoeASBM(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "BlayerKEps()" << endl;
  }

  virtual ~BlayerKEps()  {}

  void initialHook()
    {
      JoeWithModels::initialHook();

      boundVal = NULL;

      // Read inlet variable profile
      if (checkParam("READ_PROFILE"))
      {
        FILE *ifile;
        if ((ifile=fopen("./profiles.dat", "rt")) == NULL) {
          cout << "could not open the file profiles.dat" << endl;
          throw(-1);
        }

        fscanf(ifile, "n=%d\td=%d", &nn, &nval);
        boundVal = new double *[nn];
        for (int i = 0; i < nn; i++)
          boundVal[i] = new double [nval];

        // file has variables y, rho, u, v, press, ...
        for (int i=0; i<nn; i++)
          for (int v=0; v<nval; v++)
            fscanf(ifile, "%lf", &boundVal[i][v]);

        fclose(ifile);
        if (mpi_rank == 0)
          cout << "FINISHED READING INLET PROFILE" << endl;
      }

      // Specify initial condition over whole flow
      if (checkParam("SET_INIT_PROFILE"))
      {
        if(!checkDataFlag(rho))
        {
          if (boundVal == NULL)
          {
            if (mpi_rank == 0)
              cerr << "Have not read profiles.dat" << endl;
            throw(-1);
          }
          for (int icv=0; icv<ncv; icv++)
          {
            int pos=1;
            // while pos and pos-1 dont sandwich x_cv, keep increasing pos
            while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

            double fi;
            // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
            if      (x_cv[icv][1] > boundVal[pos][0])    fi = 1.0;
            // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
            else if (x_cv[icv][1] < boundVal[pos-1][0])  fi = 0.0;
            else    fi = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

            double rho1 = boundVal[pos-1][1];
            double rho2 = boundVal[pos][1];
            rho[icv] = rho1+fi*(rho2-rho1);

            double uvel1 = boundVal[pos-1][2];
            double uvel2 = boundVal[pos][2];
            rhou[icv][0] = uvel1*rho1+fi*(uvel2*rho2-uvel1*rho1);
            rhou[icv][1] = rhou[icv][2]=0.0;

            double press1 = boundVal[pos-1][4];
            double press2 = boundVal[pos][4];
            press[icv] = press1+fi*(press2-press1);

            rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv])
                      + rho[icv]*kine[icv];

            double kine1 = boundVal[pos-1][5];
            double kine2 = boundVal[pos][5];
            kine[icv] = kine1+fi*(kine2-kine1);

            double eps1 = 0.09*boundVal[pos-1][5]*boundVal[pos-1][6];
            double eps2 = 0.09*boundVal[pos][5]*boundVal[pos][6];
            eps[icv] = eps1+fi*(eps2-eps1);
          }

          updateCvData(rhou, REPLACE_ROTATE_DATA);
          updateCvData(rho, REPLACE_DATA);
          updateCvData(rhoE, REPLACE_DATA);
          updateCvData(kine, REPLACE_DATA);
          updateCvData(eps, REPLACE_DATA);
        }
      }
    }

  virtual void boundaryHook(double *T_input, double (*vel_input)[3], double *p_input, FaZone *zone)
  {
    if (zone->getNameString() == getStringParam("INLET_NAME"))
      {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }

        double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;
        gamma_bc = 1.4;
        double gm1 = gamma_bc - 1.0;
        double ovgm1 = 1.0/gm1;

        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int pos=1;
            while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

            double fi;
            if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
            else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
            else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

            double uvel1 = boundVal[pos-1][2];
            double uvel2 = boundVal[pos][2];
            u_bc[0] = uvel1 + fi*(uvel2 - uvel1);

            double vvel1 = boundVal[pos-1][3];
            double vvel2 = boundVal[pos][3];
            u_bc[1] = 0.0; //vvel1 + f*(vvel2 - vvel1);
            u_bc[2] = 0.0;

            double press1 = boundVal[pos-1][4];
            double press2 = boundVal[pos][4];
            p_bc = press1 + fi*(press2-press1);

            double temp1 = press1/(boundVal[pos-1][1]*R_gas);
            double temp2 = press2/(boundVal[pos][1]*R_gas);
            T_bc = temp1 + fi*(temp2-temp1);

            rho_bc = p_bc/(R_gas*T_bc);
            c_bc = sqrt(gamma_bc*p_bc/rho_bc);
            s_bc = pow(rho_bc,gamma_bc)/p_bc;

            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            //Compute normal velocity of freestream
            double qn_bc = vecDotVec3d(nVec,u_bc);

            //Subtract from normal velocity of mesh
            double vn_bc = qn_bc;// - normalvel_bfa[ifa];

            //Compute normal velocity of internal cell
            double qne = vecDotVec3d(nVec,vel[icv0]);

            //Compute speed of sound in internal cell
            double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

            double ac1, ac2;

            //Compute the Riemann invariants in the halo cell
            if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
              ac1 = qne   + 2.0*ovgm1*ce;
            else                     // Supersonic inflow.
              ac1 = qn_bc + 2.0*ovgm1*c_bc;


            if(vn_bc > c_bc)             // Supersonic outflow.
              ac2 = qne   - 2.0*ovgm1*ce;
            else                     // Inflow or subsonic outflow.
              ac2 = qn_bc - 2.0*ovgm1*c_bc;


            double qnf = 0.5*(ac1 + ac2);
            double cf  = 0.25*(ac1 - ac2)*gm1;

            double velf[3], sf;

            if (vn_bc > 0)                                                       // Outflow
              {
                for (int i=0; i<3; i++)
                  velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
                sf = pow( rho[icv0], gamma[icv0])/press[icv0];
              }
            else                                                                 // Inflow
              {
                for (int i=0; i<3; i++)
                  velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
                sf = s_bc;
              }
            //Compute density, pressure and velocity at boundary face
            double rho_int = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
            for (int i=0; i<3; i++)
              vel_input[ifa][i] = velf[i];
            p_input[ifa] = rho_int*cf*cf/gamma[icv0];
            T_input[ifa] = p_input[ifa]/(R_gas*rho_int);
          }
      }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbV2F_half::boundaryHookScalarRansTurb(phi_ph, zone, name);

    Param *param;

    if ((name == "kine") && (zone->getNameString() == getStringParam("INLET_NAME")))
      {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int pos=1;
            while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

            double fi;
            if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
            else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
            else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

            double kine1 = boundVal[pos-1][5];
            double kine2 = boundVal[pos][5];
            phi_ph[ifa] = kine1 + fi*(kine2-kine1);
          }
      }

    if ((name == "eps") && (zone->getNameString() == getStringParam("INLET_NAME")))
      {
      if (boundVal == NULL)
      {
        if (mpi_rank == 0)
          cerr << "Have not read profiles.dat" << endl;
        throw(-1);
      }
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int pos=1;
            while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

            double fi;
            if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
            else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
            else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

            double eps1 = 0.09*boundVal[pos-1][5]*boundVal[pos-1][6];
            double eps2 = 0.09*boundVal[pos][5]*boundVal[pos][6];
            phi_ph[ifa] =  (eps1 + fi*(eps2-eps1));

          }
      }
  }

  void transformMeshHook()
  {
    if (!checkParam("TRANSFORM_MESH")) return;

    double scl_mesh = getDoubleParam("SCL_MESH", "1.0");
    for (int ino = 0; ino < getNno(); ++ino)
      x_no[ino][0] *= scl_mesh;

    // clearFlag for wall distance to recompute the wallDist when mesh deformed
    //DoubleScalar *wallD = getScalarData("wallDist");
    //wallD->clearFlag();
  }

  virtual void temporalHook()
  {
    if (step%10000 == 0)
    {
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
};
/*
 * Flat channel with inlet and outlet bcs, SST
 */
class NonPerChanSST: public MyJoeSST{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  NonPerChanSST(char *name) : MyJoeSST(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "NonPerChanSST()" << endl;
    boundVal = NULL;
  }

  virtual ~NonPerChanSST()
  {
    if (boundVal != NULL) delete []boundVal;
  }

  void initialHook()
  {
    JoeWithModels::initialHook();

    // Read inlet variable profile
    // file has variables y, rho, u, v, press, ...
    FILE *ifile;
    if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
    {
      cout << "could not open profiles.dat, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(ifile, "n=%d\td=%d", &nn, &nval);
    boundVal = new double *[nn];
    for (int i = 0; i < nn; i++)
      boundVal[i] = new double [nval];

    for (int i=0; i<nn; i++)
      for (int v = 0; v < nval; v++)
        fscanf(ifile, "%lf", &boundVal[i][v]);

    fclose(ifile);

    // Specify initial condition over whole flow
    if (checkParam("SET_INIT_PROFILE"))
    {
      if(!checkDataFlag(rho))
      {
        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double f;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    f = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  f = 0.0;
          else    f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+f*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+f*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+f*(press2-press1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv])
              + rho[icv]*kine[icv];

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+f*(kine2-kine1);

          double omega1 = boundVal[pos-1][6];
          double omega2 = boundVal[pos][6];
          omega[icv] = omega1+f*(omega2-omega1);
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(omega, REPLACE_DATA);
      }
    }
  }

  virtual void boundaryHook(double *T_input, double (*vel_input)[3], double *p_input, FaZone *zone)
  {
    if (zone->getNameString() == getStringParam("INLET_NAME"))
    {
      double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;
      gamma_bc = 1.4;
      double gm1 = gamma_bc - 1.0;
      double ovgm1 = 1.0/gm1;

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double uvel1 = boundVal[pos-1][2];
        double uvel2 = boundVal[pos][2];
        u_bc[0] = uvel1 + f*(uvel2 - uvel1);

        double vvel1 = boundVal[pos-1][3];
        double vvel2 = boundVal[pos][3];
        u_bc[1] = 0.0; //vvel1 + f*(vvel2 - vvel1);
        u_bc[2] = 0.0;

        double press1 = boundVal[pos-1][4];
        double press2 = boundVal[pos][4];
        p_bc = press1 + f*(press2-press1);

        double temp1 = press1/(boundVal[pos-1][1]*R_gas);
        double temp2 = press2/(boundVal[pos][1]*R_gas);
        T_bc = temp1 + f*(temp2-temp1);

        rho_bc = p_bc/(R_gas*T_bc);
        c_bc = sqrt(gamma_bc*p_bc/rho_bc);
        s_bc = pow(rho_bc,gamma_bc)/p_bc;

        int icv0 = cvofa[ifa][0];
        assert( icv0 >= 0 );

        double nVec[3];
        double area = normVec3d(nVec, fa_normal[ifa]);

        //Compute normal velocity of freestream
        double qn_bc = vecDotVec3d(nVec,u_bc);

        //Subtract from normal velocity of mesh
        double vn_bc = qn_bc;// - normalvel_bfa[ifa];

        //Compute normal velocity of internal cell
        double qne = vecDotVec3d(nVec,vel[icv0]);

        //Compute speed of sound in internal cell
        double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

        double ac1, ac2;

        //Compute the Riemann invariants in the halo cell
        if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
          ac1 = qne   + 2.0*ovgm1*ce;
        else                     // Supersonic inflow.
          ac1 = qn_bc + 2.0*ovgm1*c_bc;


        if(vn_bc > c_bc)             // Supersonic outflow.
          ac2 = qne   - 2.0*ovgm1*ce;
        else                     // Inflow or subsonic outflow.
          ac2 = qn_bc - 2.0*ovgm1*c_bc;


        double qnf = 0.5*(ac1 + ac2);
        double cf  = 0.25*(ac1 - ac2)*gm1;

        double velf[3], sf;

        if (vn_bc > 0)                                                       // Outflow
        {
          for (int i=0; i<3; i++)
            velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
          sf = pow( rho[icv0], gamma[icv0])/press[icv0];
        }
        else                                                                 // Inflow
        {
          for (int i=0; i<3; i++)
            velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
          sf = s_bc;
        }
        //Compute density, pressure and velocity at boundary face
        double rho_int = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
        for (int i=0; i<3; i++)
          vel_input[ifa][i] = velf[i];
        p_input[ifa] = rho_int*cf*cf/gamma[icv0];
        T_input[ifa] = p_input[ifa]/(R_gas*rho_int);
      }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbKOmSST::boundaryHookScalarRansTurb(phi_ph, zone, name);

    Param *param;

    if ((name == "kine") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + f*(kine2-kine1);
      }
    }

    if ((name == "omega") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double omega1 = boundVal[pos-1][6];
        double omega2 = boundVal[pos][6];
        phi_ph[ifa] =  (omega1 + f*(omega2-omega1));
      }
    }
  }

  virtual void finalHook()
  {
    // Extract all data
    //writeAllValues();

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            writeWallValues(zone->getName());
      }
  }

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
      // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
      if (x_cv[icv][0] > reStress[posx][0][0])
        fx = 1.0;
      // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
      else if (x_cv[icv][0] < reStress[posx-1][0][0])
        fx = 0.0;
      else
        fx = (x_cv[icv][0] - reStress[posx-1][0][0])/(reStress[posx][0][0] - reStress[posx-1][0][0]);

      // while posy and posy-1 don't sandwich x_cv, keep increasing posy
      while(reStress[posx][posy][1] < x_cv[icv][1] && (posy < ny-1))
        posy++;
      // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
      if (x_cv[icv][1] > reStress[posx][posy][1])
        fy = 1.0;
      // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
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
 * Flat channel with inlet and outlet bc's, V2F
 */
class NonPerChanV2F: public MyJoeV2F{
protected:
  int    nn;             // number of nodes in input profile
  int    nval;           // number of variables in input profile
  double **boundVal;     // holder for input profile data

public:
  NonPerChanV2F(char *name) : MyJoeV2F(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "NonPerChanV2F()" << endl;
    boundVal = NULL;
  }

  virtual ~NonPerChanV2F()
  {
    if (boundVal != NULL) delete []boundVal;
  }

  void initialHook()
  {
    JoeWithModels::initialHook();

    if (checkParam("SET_INIT_PROFILE"))
    {
      // Read inlet variable profile
      // file has variables y, rho, u, v, press, ...
      FILE *ifile;
      if ((ifile=fopen("./profiles.dat", "rt")) == NULL)
      {
        cout << "could not open profiles.dat, apply boundary from input file" << endl;
          throw(-1);
      }

      fscanf(ifile, "n=%d\td=%d", &nn, &nval);
      boundVal = new double *[nn];
      for (int i = 0; i < nn; i++)
        boundVal[i] = new double [nval];

      for (int i=0; i<nn; i++)
        for (int v=0; v<nval; v++)
          fscanf(ifile, "%lf", &boundVal[i][v]);

      fclose(ifile);

      // Specify initial condition over whole flow

      if(!checkDataFlag(rho))
      {
        for (int icv=0; icv<ncv; icv++)
        {
          int pos=1;
          // while pos and pos-1 dont sandwich x_cv, keep increasing pos
          while(boundVal[pos][0] < x_cv[icv][1] && (pos<nn-1))       pos++;

          double fi;
          // if boundVal doesn't have a node high enough to sandwich x_cv[icv]
          if      (x_cv[icv][1] > boundVal[pos][0])    fi = 1.0;
          // if boundVal doesn't have a node low enough to sandwich x_cv[icv]
          else if (x_cv[icv][1] < boundVal[pos-1][0])  fi = 0.0;
          else    fi = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

          double rho1 = boundVal[pos-1][1];
          double rho2 = boundVal[pos][1];
          rho[icv] = rho1+fi*(rho2-rho1);

          double uvel1 = boundVal[pos-1][2];
          double uvel2 = boundVal[pos][2];
          rhou[icv][0] = uvel1*rho1+fi*(uvel2*rho2-uvel1*rho1);
          rhou[icv][1] = rhou[icv][2]=0.0;

          double press1 = boundVal[pos-1][4];
          double press2 = boundVal[pos][4];
          press[icv] = press1+fi*(press2-press1);

          double kine1 = boundVal[pos-1][5];
          double kine2 = boundVal[pos][5];
          kine[icv] = kine1+fi*(kine2-kine1);

          rhoE[icv] = press[icv]/(gamma[icv]-1.0) + 0.5/rho[icv]*vecDotVec3d(rhou[icv],rhou[icv]) +
              rho[icv]*kine[icv];

          /*double eps1 = boundVal[pos-1][6];
          double eps2 = boundVal[pos][6];
          eps[icv] = eps1+fi*(eps2-eps1);*/

          /*double v21 = boundVal[pos-1][7];
          double v22 = boundVal[pos][7];
          v2[icv] = v21 + fi*(v22 - v21);*/

          /*double f1 = boundVal[pos-1][8];
          double f2 = boundVal[pos][8];
          f[icv] = f1 + fi*(f2 - f1);*/
        }

        updateCvData(rhou, REPLACE_ROTATE_DATA);
        updateCvData(rho, REPLACE_DATA);
        updateCvData(rhoE, REPLACE_DATA);
        updateCvData(kine, REPLACE_DATA);
        updateCvData(eps, REPLACE_DATA);
        updateCvData(v2, REPLACE_DATA);
        updateCvData(f, REPLACE_DATA);
      }
    }
  }

  virtual void boundaryHook(double *T_input, double (*vel_input)[3], double *p_input, FaZone *zone)
  {
    if (zone->getNameString() == getStringParam("INLET_NAME"))
    {
      double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;
      gamma_bc = 1.4;
      double gm1 = gamma_bc - 1.0;
      double ovgm1 = 1.0/gm1;

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double f;
        if      (x_fa[ifa][1] > boundVal[pos][0])   f = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) f = 0.0;
        else    f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double uvel1 = boundVal[pos-1][2];
        double uvel2 = boundVal[pos][2];
        u_bc[0] = uvel1 + f*(uvel2 - uvel1);

        double vvel1 = boundVal[pos-1][3];
        double vvel2 = boundVal[pos][3];
        u_bc[1] = 0.0; //vvel1 + f*(vvel2 - vvel1);
        u_bc[2] = 0.0;

        double press1 = boundVal[pos-1][4];
        double press2 = boundVal[pos][4];
        p_bc = press1 + f*(press2-press1);

        double temp1 = press1/(boundVal[pos-1][1]*R_gas);
        double temp2 = press2/(boundVal[pos][1]*R_gas);
        T_bc = temp1 + f*(temp2-temp1);

        rho_bc = p_bc/(R_gas*T_bc);
        c_bc = sqrt(gamma_bc*p_bc/rho_bc);
        s_bc = pow(rho_bc,gamma_bc)/p_bc;

        int icv0 = cvofa[ifa][0];
        assert( icv0 >= 0 );

        double nVec[3];
        double area = normVec3d(nVec, fa_normal[ifa]);

        //Compute normal velocity of freestream
        double qn_bc = vecDotVec3d(nVec,u_bc);

        //Subtract from normal velocity of mesh
        double vn_bc = qn_bc;// - normalvel_bfa[ifa];

        //Compute normal velocity of internal cell
        double qne = vecDotVec3d(nVec,vel[icv0]);

        //Compute speed of sound in internal cell
        double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

        double ac1, ac2;

        //Compute the Riemann invariants in the halo cell
        if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
          ac1 = qne   + 2.0*ovgm1*ce;
        else                     // Supersonic inflow.
          ac1 = qn_bc + 2.0*ovgm1*c_bc;


        if(vn_bc > c_bc)             // Supersonic outflow.
          ac2 = qne   - 2.0*ovgm1*ce;
        else                     // Inflow or subsonic outflow.
          ac2 = qn_bc - 2.0*ovgm1*c_bc;


        double qnf = 0.5*(ac1 + ac2);
        double cf  = 0.25*(ac1 - ac2)*gm1;

        double velf[3], sf;

        if (vn_bc > 0)                                                       // Outflow
        {
          for (int i=0; i<3; i++)
            velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
          sf = pow( rho[icv0], gamma[icv0])/press[icv0];
        }
        else                                                                 // Inflow
        {
          for (int i=0; i<3; i++)
            velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
          sf = s_bc;
        }
        //Compute density, pressure and velocity at boundary face
        double rho_int = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
        for (int i=0; i<3; i++)
          vel_input[ifa][i] = velf[i];
        p_input[ifa] = rho_int*cf*cf/gamma[icv0];
        T_input[ifa] = p_input[ifa]/(R_gas*rho_int);
      }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbV2F::boundaryHookScalarRansTurb(phi_ph, zone, name);

    Param *param;

    if ((name == "kine") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + fi*(kine2-kine1);
      }
    }

    if ((name == "eps") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double eps1 = boundVal[pos-1][6];
        double eps2 = boundVal[pos][6];
        phi_ph[ifa] =  (eps1 + fi*(eps2-eps1));
      }
    }
    if ((name == "v2") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double v21 = boundVal[pos-1][7];
        double v22 = boundVal[pos][7];
        phi_ph[ifa] =  (v21 + fi*(v22-v21));
      }
    }
    if ((name == "f") && (zone->getNameString() == getStringParam("INLET_NAME")))
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while ((boundVal[pos][0] < x_fa[ifa][1]) && (pos < nn-1))      pos++;

        double fi;
        if      (x_fa[ifa][1] > boundVal[pos][0])   fi = 1.0;
        else if (x_fa[ifa][1] < boundVal[pos-1][0]) fi = 0.0;
        else    fi = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double f1 = boundVal[pos-1][8];
        double f2 = boundVal[pos][8];
        phi_ph[ifa] =  (f1 + fi*(f2-f1));
      }
    }
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
    case 6:   joe = new PerChanSST(inputFileName);      break;
    case 7:   joe = new PerChanV2F(inputFileName);      break;
    case 8:   joe = new PerChanKEps(inputFileName);     break;
    case 9:   joe = new BlayerSST(inputFileName);       break;
    case 10:  joe = new BlayerV2F(inputFileName);       break;
    case 11:  joe = new BlayerKEps(inputFileName);      break;
    case 12:  joe = new NonPerChanSST(inputFileName);   break;
    case 13:  joe = new NonPerChanV2F(inputFileName);   break;
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



