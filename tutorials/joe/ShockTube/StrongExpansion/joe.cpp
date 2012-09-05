#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
//#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
//#include "turbModels/TurbModel_V2F.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels //, public RansTurbKOmSST, public RansCombSteadyFlamelet
{
public:
  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoe()" << endl;
  }

  virtual ~MyJoe()  {}

  virtual void initialHook()
  {
    double pL   = getDoubleParam("P_L", "1.0");
    double rhoL = getDoubleParam("RHO_L", "1.0");
    double uL   = getDoubleParam("U_L", "0.0");

    double pR   = getDoubleParam("P_R", "0.1");
    double rhoR = getDoubleParam("RHO_R", "0.125");
    double uR   = getDoubleParam("U_R", "0.0");

    for (int icv = 0; icv < ncv; icv++)           // loop over all cvs
    {
      if (x_cv[icv][0] < 0.0)                     // set left shock tube side based on cell center coordinates (xcv[icv][0] = x)
      {
        rho[icv]     = rhoL;
        rhou[icv][0] = rhoL * uL;
        rhou[icv][1] = rhou[icv][2] = 0.0;
        rhoE[icv]    = pL/(gamma[icv]-1.0) + 0.5 * rhoL * uL * uL;
      }
      else                                        // set right side
      {
        rho[icv]     = rhoR;
        rhou[icv][0] = rhoR * uR;
        rhou[icv][1] = rhou[icv][2] = 0.0;
        rhoE[icv]    = pR/(gamma[icv]-1.0) + 0.5 * rhoR * uR * uR;
      }
    }

    // update variables of inter-processor ghost cells
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }

  /*
   * write wall values to file in tecplot format
   */
  void finalHook()
  {
    int cv_global = 0;
    MPI_Allreduce(&ncv, &cv_global, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "shocktube.dat");

    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"x\" \"rho\" \"Ma\" \"press\" \"temp\" \n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", cv_global);
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

    for (int icv = 0; icv < ncv; icv++)
    {
      double Mach = vel[icv][0]/sos[icv];
      fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n", x_cv[icv][0], rho[icv], Mach, press[icv], temp[icv]);
    }
    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }
};


// ###########################################################################################
// ------                                                                               ------
// ------                    joe with scalar                                            ------
// ------                                                                               ------
// ###########################################################################################
class MyJoeWithScal : public MyJoe
{
public:
  double *phi;

public:
  MyJoeWithScal(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWithScal()" << endl;

    // register scalar transport equation name="phi"
    ScalarTranspEq *eq;
    eq = registerScalarTransport("phi", CV_DATA);
    eq->phiZero = 1.0e-5;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 20;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
//    eq->turbSchmidtNumber = 1.0; // not needed for Euler flow
  }

  virtual ~MyJoeWithScal()  {}

  void initialHook()
  {
    // call initial hook from MyJoe to initialize rho, rhou and rhoE
    MyJoe::initialHook();

    ScalarTranspEq *eq = getScalarTransportData("phi");
    phi    = eq->phi;

    double PhiL = getDoubleParam("PHI_L", "1.0");
    double PhiR = getDoubleParam("PHI_R", "0.5");

    for (int icv = 0; icv < ncv; icv++)
    {
      if (x_cv[icv][0] < 0.0)       phi[icv] = PhiL;
      else                          phi[icv] = PhiR;
    }
    //updateCv2DataByName("phi", REPLACE_DATA);
  }

  void finalHook()
  {
    int cv_global = 0;
    MPI_Allreduce(&ncv, &cv_global, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "shocktube.dat");

    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"rho\" \"M\" \"P\" \"T\" \"phi\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", cv_global);
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

    for (int icv = 0; icv < ncv; icv++)
    {
      double Mach = vel[icv][0]/sos[icv];
      fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n", x_cv[icv][0], rho[icv], Mach, press[icv], temp[icv], phi[icv]);
    }
    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }

};


int main(int argc, char * argv[])
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

    JoeWithModels *joe;
    
    switch (run)
    {
    case 0:   joe = new MyJoe(inputFileName);    	   break;
    case 1:   joe = new MyJoeWithScal(inputFileName);      break;
    default: 
      if (mpi_rank == 0)
        cerr << "ERROR: run number not available!" << endl;
    }
    
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();
    
    joe->run();
    
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

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



