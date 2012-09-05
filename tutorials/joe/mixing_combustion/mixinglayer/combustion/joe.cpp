#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
//#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
//#include "turbModels/TurbModel_V2F.h"

#include "combModels/CombModel_VariableProperties.h"
#include "combModels/CombModel_BinaryMixing.h"
#include "combModels/CombModel_Mixing.h"
//#include "combModels/CombModel_SteadyFlamelet.h"
#include "combModels/CombModel_FPVA.h"
#include "combModels/CombModel_FPVA_Coeff.h"

#define COEFF


// ###########################################################################################
// ------                                                                               ------
// ------                    Joe for mixing layer combustion                            ------
// ------                                                                               ------
// ###########################################################################################
#ifdef COEFF
//class mixlayercomb : public JoeWithModels, public RansCombFPVA_Coeff<ChemtableCartesianCubic>
class mixlayercomb : public JoeWithModels, public RansCombFPVA_Coeff<ChemtableAdaptiveLinear>
#else
//class mixlayercomb : public JoeWithModels, public RansCombFPVA<ChemtableCartesianCubic>
class mixlayercomb : public JoeWithModels, public RansCombFPVA<ChemtableAdaptiveLinear>
#endif

{
public:

  mixlayercomb(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "mixlayercomb()" << endl;
  }

  virtual ~mixlayercomb()  {}


  void initialHook()
  {
    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;
    
    // First initialize scalar ZMean
    double ZM_upper = getDoubleParam("upperinlet.ZMean", "0.0");
    double ZM_lower = getDoubleParam("lowerinlet.ZMean", "1.0");

    if (!checkScalarFlag("ZMean"))
    {   
      for (int icv = 0; icv < ncv; icv++)
        if (x_cv[icv][1] > 0.0)       
          ZMean[icv] = ZM_upper;
        else                      
	  ZMean[icv] = ZM_lower;
    }
    updateCvDataByName("ZMean", REPLACE_DATA);

    if (!checkScalarFlag("ZVar"))
    {   
      for (int icv = 0; icv < ncv; icv++)
        ZVar[icv] = 0.0;
    }
    updateCvDataByName("ZVar", REPLACE_DATA);

    if (!checkScalarFlag("CMean"))
    {   
      for (int icv = 0; icv < ncv; icv++)
#ifdef COEFF
	CMean[icv] = myChemTable.Lookup(ZMean[icv], ZVar[icv], 0.05, "PROG");
#else
        myMixture.GetProgressVariable(CMean[icv], ZMean[icv], ZVar[icv], 0.05, myChemTable);
#endif
    }
    updateCvDataByName("CMean", REPLACE_DATA);
    
//     for (int icv = 0; icv < ncv; icv++)
//     {
//       if (!checkScalarFlag("ZMeanRelaxed"))
// 	ZMeanRelaxed[icv] = ZMean[icv];
//       if (!checkScalarFlag("CMeanRelaxed"))
// 	CMeanRelaxed[icv] = CMean[icv];
//     }
//     updateCvDataByName("ZMeanRelaxed", REPLACE_DATA);
//     updateCvDataByName("CMeanRelaxed", REPLACE_DATA);


    // Then initialize flow field
    Param  *p;
    double dum;
    double U_upper, V_upper, W_upper, T_upper, P_upper;
    getParam(p, "upperinlet");
    if (p->getString() == "CBC")
    {
      U_upper = p->getDouble(2);
      V_upper = p->getDouble(3);
      W_upper = p->getDouble(4);
      T_upper = p->getDouble(5);
      P_upper = p->getDouble(6);
    }
    else
    {
      cerr << "### Missing CBC upperinlet conditions! ###" << endl;
      throw(-1);
    }

    double U_lower, V_lower, W_lower, T_lower, P_lower;
    getParam(p, "lowerinlet");
    if (p->getString() == "CBC")
    {
      U_lower = p->getDouble(2);
      V_lower = p->getDouble(3);
      W_lower = p->getDouble(4);
      T_lower = p->getDouble(5);
      P_lower = p->getDouble(6);
    }
    else
    {
      cerr << "### Missing CBC lowerinlet conditions! ###" << endl;
      throw(-1);
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      double U, V, W, T, P, h, R;

      if (x_cv[icv][1] > 0.0)       
      {
	U = U_upper;
        V = V_upper;
        W = W_upper;
        T = T_upper;
        P = P_upper;
      }
      else                      
      {
	U = U_lower;
        V = V_lower;
        W = W_lower;
        T = T_lower;
        P = P_lower;
      }

      ComputeProperties_T(h, R, dum, dum, dum, T, ZMean[icv], ZVar[icv], CMean[icv]);
      
      if (!checkDataFlag(rho)) 
        rho[icv] = P / (R * T);  
          
      if (!checkDataFlag(rhou))
      {
        rhou[icv][0] = rho[icv] * U;
        rhou[icv][1] = rho[icv] * V;
        rhou[icv][2] = rho[icv] * W;
      }

      double kin = 0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
      if (!checkDataFlag(rhoE))
        rhoE[icv] = h * rho[icv] - P + kin;
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

    if (mpi_rank == 0)
      cout << "Initial Conditions Set!"<< endl;
  }

//  virtual void temporalHook(){}

//  virtual void finalHook(){}
  
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
//  virtual void calcRansTurbViscMuet(double *rho, double (*rhou)[3]) {/*empty*/}
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
  
//  myJoe(ParamMap &p, int i) : JoeWithModels(p), UgpWithCvCompFlow(p), ownID(i) {
//    if (mpi_rank == 0)      cout << "myJoe()" << endl;
//    init();
//  }

//  myJoe(char *name, int i) : JoeWithModels(name), UgpWithCvCompFlow(name), ownID(i) {
//    if (mpi_rank == 0)
//      cout << "myJoe()" << endl;
//    init();
//  }
  
};


/*
 * mixlayercomb with Spalart & Allmaras Model
 */
class mixlayercomb_SA : public mixlayercomb, public RansTurbSA
{
public:
  mixlayercomb_SA(char *name) : mixlayercomb(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "mixlayercomb_SA()" << endl;
  }

  virtual ~mixlayercomb_SA() {}
};

/*
 * mixlayercomb with Menter SST Model
 */
class mixlayercomb_SST : public mixlayercomb, public RansTurbKOmSST
{
public:
  mixlayercomb_SST(char *name) : mixlayercomb(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "mixlayercomb_SST()" << endl;
  }

  virtual ~mixlayercomb_SST() {}
};

/*
 * mixlayercomb with Wilcox k-omega Model
 */
class mixlayercomb_WX : public mixlayercomb, public RansTurbKOm
{
public:
  mixlayercomb_WX(char *name) : mixlayercomb(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "mixlayercomb_WX()" << endl;
  }

  virtual ~mixlayercomb_WX() {}
};




int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
  int run = 2;

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
      case 1:   joe = new mixlayercomb_SA(inputFileName);      break;
      case 2:   joe = new mixlayercomb_SST(inputFileName);     break;
      case 3:   joe = new mixlayercomb_WX(inputFileName);      break;
      case 4:   joe = new mixlayercomb(inputFileName);         break;
      default:  joe = new mixlayercomb_SST(inputFileName);
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



