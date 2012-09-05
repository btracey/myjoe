#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"
//#include "turbModels/TransModel_GaReT.h"


#include "combModels/CombModel_BinaryMixing.h"
#include "combModels/CombModel_VariableProperties.h"
#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_FPVA.h"
#include "combModels/CombModel_FPVA_Coeff.h"





class MyJoe : public JoeWithModels
{
public:

  double *Mach;

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    Mach = NULL;   registerScalar(Mach, "Mach", CV_DATA);
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;
  }

  virtual ~MyJoe()  {}

  void transformMeshHook()
  {
    if (!checkParam("TRANSFORM_MESH")) return;

    double Amax = getDoubleParam("AMAX", "0.2");
    double Amin = getDoubleParam("AMIN", "0.02");

    for (int ino = 0; ino < getNno(); ++ino)
    {	
      if (x_no[ino][1] > 0.0)        x_no[ino][1] =  (Amax-Amin)*x_no[ino][0]*x_no[ino][0]+Amin;
      if (x_no[ino][1] < 0.0)        x_no[ino][1] = -(Amax-Amin)*x_no[ino][0]*x_no[ino][0]-Amin;
    }	
  }

  virtual void temporalHook() 
  {
    for (int icv = 0; icv < ncv; icv++)
      Mach[icv] = sqrt(vel[icv][0]*vel[icv][0] + vel[icv][1]*vel[icv][1] + vel[icv][1]*vel[icv][1]) / sos[icv];
    updateCvData(Mach, REPLACE_DATA);
  }

  virtual void finalHook() {}
};


/*
 * MyJoe with RansCombVarProperties
 */
class MyJoeVP : public MyJoe, public RansCombVarProperties
{
public:
  MyJoeVP(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeVP()" << endl;
  }

  virtual void initialHook()
  {
    Param *p;
    double T_in;
    double p_in;
    double u_in, v_in, w_in;
    double h_in, RoM_in, Gam_in, dum1, dum2;

    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;

    getParam(p, "inlet");
    if (p->getString() == "CBC_SUBSONIC_INLET")
    {
      u_in = 0.0;
      v_in = 0.0;
      w_in = 0.0;
      T_in = p->getDouble(5);
      p_in = p->getDouble(6);
    }
    else
    {
      cerr << "### Missing/wrong inflow boundary conditions ###" << endl;
      throw(-1);
    }

    ComputeProperties_T(h_in, RoM_in, Gam_in, dum1, dum2, T_in);

    if (!checkDataFlag(rho))
    {
      for (int icv = 0; icv < ncv; icv++)
        rho[icv] = p_in / (RoM_in * T_in);
    }

    if (!checkDataFlag(rhou))
    {
      for (int icv = 0; icv < ncv; icv++)
      {
        rhou[icv][0] = rho[icv] * u_in;
        rhou[icv][1] = rho[icv] * v_in;
        rhou[icv][2] = rho[icv] * w_in;
      }
    }

    if (!checkDataFlag(rhoE))
    {
      for (int icv = 0; icv < ncv; icv++)
        rhoE[icv] = h_in * rho[icv] - p_in +
                    0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
    }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }

  virtual ~MyJoeVP() {}
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
    case 1:   joe = new MyJoeVP(inputFileName);      break;
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



