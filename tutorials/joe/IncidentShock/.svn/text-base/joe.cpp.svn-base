#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

#include "myMem.h"




// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels
{
public:

  int npos, nval;
  double **boundVal;

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;

    // read inlet profile
    FILE *fp;

    if ((fp=fopen("inlet.txt", "rt")) == NULL)
    {
      if (mpi_rank == 0)
        cerr << "could not open inlet.txt, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(fp, "n=%d\td=%d", &npos, &nval);

    getMem2D(&boundVal, 0, npos-1, 0, nval-1, "boundVal");

    // file has values
    // y, rho, u, v, p, kine, omega
    for (int i=0; i<npos; i++)
      for (int v=0; v<nval; v++)
        fscanf(fp, "%lf", &boundVal[i][v]);

    if (mpi_rank == 0)
    {
      printf("n=%d\td=%d\n", npos, nval);

      for (int i=0; i<npos; i++)
      {
        for (int v=0; v<nval; v++)
          printf("%.6le\t", boundVal[i][v]);

        printf("\n");
      }
    }
  }

  virtual ~MyJoe()  {}

  void initialHook()
  {
    if (!checkDataFlag(rho))
    for (int icv=0; icv<ncv; icv++)
    {
      int pos=0;
      while(boundVal[pos][0] < x_cv[icv][1])      pos++;
      double f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

      double rho1 = boundVal[pos-1][1];
      double rho2 = boundVal[pos][1];
      rho[icv] = rho1 + f*(rho2-rho1);

      double rhou1 = boundVal[pos-1][2]*rho1;
      double rhou2 = boundVal[pos][2]*rho2;
      rhou[icv][0] = rhou1 + f*(rhou2-rhou1);

      double rhov1 = boundVal[pos-1][3]*rho1;
      double rhov2 = boundVal[pos][3]*rho2;
      rhou[icv][1] = rhov1 + f*(rhov2-rhov1);

      rhou[icv][2] = 0.0;

      double rhoe1 = boundVal[pos-1][4]/(gamma[icv] - 1.0) + 0.5*(rhou1*rhou1+rhov1*rhov1)/rho1;
      double rhoe2 = boundVal[pos][4]/(gamma[icv] - 1.0) + 0.5*(rhou2*rhou2+rhov2*rhov2)/rho2;
      rhoE[icv] = rhoe1 + f*(rhoe2-rhoe1);
    }

    // update all data
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }

  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone)
  {
    if (zone->getNameString() == "inlet")   // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;
        while(boundVal[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double rho1 = boundVal[pos-1][1];
        double rho2 = boundVal[pos][1];
        double rho = rho1 + f*(rho2-rho1);

        double uvel1 = boundVal[pos-1][2];
        double uvel2 = boundVal[pos][2];
        vel[ifa][0] = uvel1 + f*(uvel2-uvel1);

        double vvel1 = boundVal[pos-1][3];
        double vvel2 = boundVal[pos][3];
        vel[ifa][1] = vvel1 + f*(vvel2-vvel1);
        vel[ifa][2] = 0.0;

        double press1 = boundVal[pos-1][4];
        double press2 = boundVal[pos][4];
        press[ifa] = press1 + f*(press2-press1);

        temp[ifa] = press[ifa]/(rho*R_gas);
      }
    }
  }
  
  
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

    for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone != faZoneList.end(); faZone++)
      if (faZone->getKind() == FA_ZONE_BOUNDARY)
        if (faZone->getNameString() == name)
        {
          for (int ifa = faZone->ifa_f; ifa <= faZone->ifa_l; ifa++)
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

  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall");
  }

  virtual void finalHook()
  {
    writeWallValues("wall");
  }
};

// ###########################################################################################
/*
 * MyJoe with Spalart & Allmaras Model
 */
// ###########################################################################################
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

// ###########################################################################################
/*
 * MyJoe with Menter SST Model
 */
// ###########################################################################################
class MyJoeSST : public MyJoe, public RansTurbKOmSST
{
public:
  MyJoeSST(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSST()" << endl;
  }

  virtual ~MyJoeSST() {}

  virtual void initialHookScalarRansTurbModel()
  {
    RansTurbKOmSST::initialHookScalarRansTurbModel();

    if (!checkScalarFlag("kine"))
    for (int icv=0; icv<ncv; icv++)
    {
      int pos=0;
      while(boundVal[pos][0] < x_cv[icv][1])      pos++;
      double f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

      double kine1 = boundVal[pos-1][5];
      double kine2 = boundVal[pos][5];
      kine[icv] = kine1 + f*(kine2-kine1);
    }

    if (!checkScalarFlag("omega"))
    for (int icv=0; icv<ncv; icv++)
    {
      int pos=0;
      while(boundVal[pos][0] < x_cv[icv][1])      pos++;
      double f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

      double omega1 = boundVal[pos-1][6];
      double omega2 = boundVal[pos][6];
      omega[icv] = omega1 + f*(omega2-omega1);
    }

    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("omega", REPLACE_DATA);
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbKOmSST::boundaryHookScalarRansTurb(phi_ph, zone, name);

    if ((name == "kine") && (zone->getNameString() == "inlet"))   // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;
        while(boundVal[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + f*(kine2-kine1);
      }
    }

    if ((name == "omega") && (zone->getNameString() == "inlet"))  // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;
        while(boundVal[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double omega1 = boundVal[pos-1][6];
        double omega2 = boundVal[pos][6];
        phi_ph[ifa] = omega1 + f*(omega2-omega1);
      }
    }
  }
};

// ###########################################################################################
/*
 * MyJoe with Wolcox k-omega Model
 */
// ###########################################################################################
class MyJoeWX : public MyJoe, public RansTurbKOm
{
public:
  MyJoeWX(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWX()" << endl;
  }

  virtual ~MyJoeWX() {}

  virtual void initialHookScalarRansTurbModel()
  {
    RansTurbKOm::initialHookScalarRansTurbModel();

    if (!checkScalarFlag("kine"))
    for (int icv=0; icv<ncv; icv++)
    {
      int pos=0;
      while(boundVal[pos][0] < x_cv[icv][1])      pos++;
      double f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

      double kine1 = boundVal[pos-1][5];
      double kine2 = boundVal[pos][5];
      kine[icv] = kine1 + f*(kine2-kine1);
    }

    if (!checkScalarFlag("omega"))
    for (int icv=0; icv<ncv; icv++)
    {
      int pos=0;
      while(boundVal[pos][0] < x_cv[icv][1])      pos++;
      double f = (x_cv[icv][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

      double omega1 = boundVal[pos-1][6];
      double omega2 = boundVal[pos][6];
      omega[icv] = omega1 + f*(omega2-omega1);
    }

    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("omega", REPLACE_DATA);
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)
  {
    RansTurbKOm::boundaryHookScalarRansTurb(phi_ph, zone, name);

    if ((name == "kine") && (zone->getNameString() == "inlet"))   // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;
        while(boundVal[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double kine1 = boundVal[pos-1][5];
        double kine2 = boundVal[pos][5];
        phi_ph[ifa] = kine1 + f*(kine2-kine1);
      }
    }

    if ((name == "omega") && (zone->getNameString() == "inlet"))  // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;
        while(boundVal[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal[pos-1][0])/(boundVal[pos][0]-boundVal[pos-1][0]);

        double omega1 = boundVal[pos-1][6];
        double omega2 = boundVal[pos][6];
        phi_ph[ifa] = omega1 + f*(omega2-omega1);
      }
    }
  }
};




// ###########################################################################################

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
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
    case 0:   joe = new MyJoe(inputFileName);        break;
    case 1:   joe = new MyJoeSA(inputFileName);      break;
    case 2:   joe = new MyJoeSST(inputFileName);     break;
    case 3:   joe = new MyJoeWX(inputFileName);      break;
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



