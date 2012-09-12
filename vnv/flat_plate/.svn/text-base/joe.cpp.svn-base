#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

//#include "combModels/CombModel_BaseScalar.h"
//#include "combModels/CombModel_Mixing.h"
//#include "combModels/CombModel_SteadyFlamelet.h"


// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels 
{
  double *grad;
public:

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;   
    
    // register grad
    grad = NULL;        registerScalar(grad, "grad", CV_DATA);
  }

  virtual ~MyJoe()  {}

  void initialHook() {JoeWithModels::initialHook();}  
  
  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name, string fname)
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
    char charfname[200];
    sprintf(charfname, "%s", fname.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(charfname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << charfname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"rho\" \"press\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"Rex\" \"cflam\" \"cfturb\" \"cf\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(charfname,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    
    double muref = getDoubleParam("MU_REF");
    double rhotot = getDoubleParam("RHO_REF");
    double ptot = getDoubleParam("P_REF");

    double pstat;
    Param *p;
    if (getParam(p,"outlet"))
      pstat = p->getDouble(2);

    double Mach = sqrt(2.0/(GAMMA-1.0)*(pow(pstat/ptot, -(GAMMA-1)/GAMMA)-1.0));
    double rhostat = rhotot*pow(1.0+(GAMMA-1.0)/2.0*Mach*Mach, -1.0/(GAMMA-1.0));

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
            double mulam = calcMuLam(wallTemp);
            
            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;
            
            
            double Rex = x_fa[ifa][0]*rhostat*Mach/muref; 
            double cflam = 0.664/sqrt(Rex);
            double cfturb = 0.445/pow(log(0.06*Rex), 2.0);

            double cf = tau/(0.5*rhostat*Mach*Mach);
            
            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv], temp[icv], tau, qDot, yplus, Rex, cflam, cfturb, cf);
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
      writeWallValues("wall", "wallLam.dat");
  }

  virtual void finalHook()
  {
    writeWallValues("wall", "wallLam.dat");
  }
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
  
  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall", "wallSA.dat");
  }

  virtual void finalHook()
  {
    writeWallValues("wall", "wallSA.dat");
  }
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
  
  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall", "wallSST.dat");
  }

  virtual void finalHook()
  {
    writeWallValues("wall", "wallSST.dat");
  }
};

/*
 * MyJoe with Wolcox k-omega Model
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
  
  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall", "wallWX.dat");
  }

  virtual void finalHook()
  {
    writeWallValues("wall", "wallWX.dat");
  }
};






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



