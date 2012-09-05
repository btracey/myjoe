#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"
#include "turbModels/TransModel_GaReT.h"


#include "combModels/CombModel_BinaryMixing.h"
#include "combModels/CombModel_VariableProperties.h"
#include "combModels/CombModel_VariablePropertiesAnalytic.h"
#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_FPVA.h"
#include "combModels/CombModel_FPVA_Coeff.h"




// ###########################################################################################
// ------                                                                               ------
// ------                     MyJoe for Cylinder case - Vanilla version                 ------
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

  void initialHook() {JoeWithModels::initialHook();}  
  
};


// ###########################################################################################
// ------                                                                               ------
// ------            MyJoe for Cylinder case - Variable properties analytic             ------
// ------                                                                               ------
// ###########################################################################################
class MyJoeVarPropAnal : public JoeWithModels, public RansCombVarPropertiesAnalytic
{
public:

  MyJoeVarPropAnal(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)      cout << "MyJoeVarPropAnal()" << endl;
  }

  virtual ~MyJoeVarPropAnal()  {}

  virtual void initialHook()
  {
    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;

    double rho_init, p_init, u_init[3] = {0.0, 0.0, 0.0}, T_init, h_init, R, gamma_init;

    Param *p;
    if (getParam(p, "RHO_INITIAL"))       rho_init = p->getDouble(1);
    else                                  rho_init = rho_ref;

    if (getParam(p, "P_INITIAL"))         p_init = p->getDouble(1);
    else                                  p_init = p_ref;

    if (getParam(p, "U_INITIAL"))
    {
      u_init[0] = p->getDouble(1);
      u_init[1] = p->getDouble(2);
      u_init[2] = p->getDouble(3);
    }
    calcThermoProp_p(T_init, h_init, R, gamma_init, rho_init, p_init, NULL, 0);


    if (!checkDataFlag(rho))
      for (int icv=0; icv<ncv; icv++)
        rho[icv] = rho_init;

    if (!checkDataFlag(rhou))
      for (int icv=0; icv<ncv; icv++)
        for (int i=0; i<3; i++)
          rhou[icv][i] = rho_init*u_init[i];

    if (!checkDataFlag(gamma)) {
      for (int icv=0; icv<ncv; icv++)
        gamma[icv] = gamma_init ;
    }

    if (!checkDataFlag(rhoE))
      for (int icv=0; icv<ncv; icv++)
        rhoE[icv] = h_init - p_init/rho_init + 0.5*rho_init*vecDotVec3d(u_init, u_init);


    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

  }
};

// ###########################################################################################
// ------                                                                               ------
// ------            MyJoe for Cylinder case - Variable properties                      ------
// ------                                                                               ------
// ###########################################################################################
class MyJoeVarProp : public JoeWithModels, public RansCombVarProperties
{
public:

  MyJoeVarProp(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)      cout << "MyJoeVarProp()" << endl;
  }

  virtual ~MyJoeVarProp()  {}

  virtual void initialHook()
  {
    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;

    double rho_init, p_init, u_init[3] = {0.0, 0.0, 0.0}, T_init, h_init, R, gamma_init;

    Param *p;
    if (getParam(p, "RHO_INITIAL"))       rho_init = p->getDouble(1);
    else                                  rho_init = rho_ref;

    if (getParam(p, "P_INITIAL"))         p_init = p->getDouble(1);
    else                                  p_init = p_ref;

    if (getParam(p, "U_INITIAL"))
    {
      u_init[0] = p->getDouble(1);
      u_init[1] = p->getDouble(2);
      u_init[2] = p->getDouble(3);
    }
    calcThermoProp_p(T_init, h_init, R, gamma_init, rho_init, p_init, NULL, 0);


    if (!checkDataFlag(rho))
      for (int icv=0; icv<ncv; icv++)
        rho[icv] = rho_init;

    if (!checkDataFlag(rhou))
      for (int icv=0; icv<ncv; icv++)
        for (int i=0; i<3; i++)
          rhou[icv][i] = rho_init*u_init[i];

    if (!checkDataFlag(gamma)) {
      for (int icv=0; icv<ncv; icv++)
        gamma[icv] = gamma_init ;
    }

    if (!checkDataFlag(rhoE))
      for (int icv=0; icv<ncv; icv++)
        rhoE[icv] = h_init - p_init/rho_init + 0.5*rho_init*vecDotVec3d(u_init, u_init);


    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

  }
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
    case 0:   joe = new MyJoe(inputFileName);                 break;
    case 1:   joe = new MyJoeVarPropAnal(inputFileName);      break;
    case 2:   joe = new MyJoeVarProp(inputFileName);          break;
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



