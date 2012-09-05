#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"
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

  void initialHook()
  {

    double delta = getDoubleParam("SHOCK_THICKNESS", "1.0");    

    double Mach1 = 1.29;
    double rho1 = 1.0;
    double T1 = 1.0;

    double a1 = 1.0;
    double u1 = a1*Mach1;

    double M2 = sqrt((1.0+((GAMMA-1.0)/2.0)*Mach1*Mach1)/(GAMMA*Mach1*Mach1-(GAMMA-1.0)/2.0));
    double T2 = T1*(1.0+((2.0*GAMMA)/(GAMMA+1.0))*(Mach1*Mach1-1.0))*((2.0+(GAMMA-1.0)*Mach1*Mach1)/((GAMMA+1.0)*Mach1*Mach1));
    double rho2 = rho1*(Mach1*Mach1/(1.0+(GAMMA-1.0)/(GAMMA+1.0)*(Mach1*Mach1-1.0)));
    double u2 = M2*a1*sqrt(T2/T1);

    for (int icv = 0; icv < ncv; icv++)
    {
      double x = x_cv[icv][0]-2.0;
      double u = ((u1-u2)/2.0)*tanh(-x/(delta/2.0))+((u1+u2)/2.0);
      double rh = rho1*u1/u;
      double T = ((T2-T1)/2.0)*tanh(x/(delta/2.0))+((T1+T2)/2.0);

      rho[icv] = rh;
      rhou[icv][0] = rh*u;
      rhou[icv][1] = rhou[icv][2] = 0.0;
      rhoE[icv] = rh*R_gas*T/(GAMMA-1.0) + 0.5*rh*u*u;
    }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
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



