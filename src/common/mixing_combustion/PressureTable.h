#ifndef PRESSURETABLE_H
#define PRESSURETABLE_H
#include "Table4D.h"
#include <string>
#include <vector>
#include <map>
using std::string;
using std::vector;
using std::map;
using std::pair;
#include "./Numerical_Recipes/roots_multidim.h"

/*! /brief Class PressureTable describes a multiple pressure flamelet database
 *
 */
class PressureTable : public Table4D
{
public:

  //------------------------------//
  //         Constructors         //
  //------------------------------//
  PressureTable();
  PressureTable(string filename, int nPress, int nZm, int nZv, int nC, int nVar);
  //------------------------------//
  //         Destructor           //
  //------------------------------//
  ~PressureTable();
  //------------------------------//
  //     Get/Set Attributes       //
  //------------------------------//
  inline double GetTableVersion()  {return version;}
  inline string GetTableFileName() {return TableFilename;}
  inline string GetTableCombustionModel() {return CombustionModel;}
  inline string GetTableType() {return ChemTableType;}
  inline string GetVarName(int i) {return VarNames[i];}
  inline int    GetVarNameIndex(string Var) {return VarNamesMap.find(Var)->second;}
  inline static double GetPscale() {return P_scale;}

  inline void SetTableVersion(double ver) {version = ver;}
  inline void SetTableFileName(string name) {TableFilename = name;}
  inline void SetTableCombustionModel(string model) {CombustionModel = model;}
  inline void SetTableType(string type) {ChemTableType = type;}
  void AddVarName(string var);

  //------------------------------//
  //            Methods           //
  //------------------------------//

  //-----------------------------------------
  /*! \brief Write the database to the binary format
   *  All class attributes must be filled
   *  set output to= "debug" to display verbose and debugging
   */
  void WriteTable(string output);
  //-----------------------------------------
  /*! \brief Load a database binary file
  *   All class attributes are initialized
  *   filename is the name of the file to read
  *   Example: Load("my_pretty_database");
  */
  void Load(string filename);
  //-----------------------------------------
  // Deallocate and reset the class object
  //-----------------------------------------
  void Unload(void);
  //-----------------------------------------
  /*! \brief Interpolate the quantity "variable" in the 4D database
   * Example: Lookup(pressure_value, mixture_value, mixture_variance_value, prog_value, "my_variable");
   */
  double Lookup(double P, double Zm, double Zvar, double C, string variable);
  double Lookup(double P, double Zm, double Zvar, double C, int ivar);
  //-----------------------------------------
  /*! \brief Interpolate the quantity "variable" in the 4D database
   *  Interpolation index and weight are assumed to be already set
   *  Example: Lookup("my_variable");
   */
  double Lookup(string variable);
  double Lookup(int ivar);
  //-----------------------------------------
  /*! \brief Interpolate the first 10 variables in the database
   */
  void LookupCoeff(double &RoM,      double &T0,      double &E0,
                   double &Gam0,     double &a_Gam,   double &mu0,
                   double &a_mu,     double &lamOcp0, double &a_lamOcp,
                   double &src_prog, double P,        double Zm,
                   double Zvar,      double C);
  //-----------------------------------------
  /*! \brief Interpolate the first 5 variables in the database
   */
  void LookupSelectedCoeff(double &RoM,  double &T0,    double &E0,
                           double &Gam0, double &a_Gam, double P,
                           double Zm,    double Zvar,   double C);
  //-----------------------------------------
  /*! \brief Print Table information
   */
  void print(void);
  //-----------------------------------------
  /*! \brief compute normalized variance 
   */
  double ComputeNormalizedVariance(const double Zm, const double Zvar);
  //-----------------------------------------
  /*! \brief compute normalized progress variable 
   */
  double ComputeNormalizedProg(const double P, const double Zm, const double Sz, const double Yc);
  //-----------------------------------------
  /*! \brief compute database input from conservative variables
   */
  void ComputeDatabaseInputs(const double rho, const double Z, const double Zv, const double Yc,
                             double &P, double &Sz ,  double &C,
                             double P0=1.0e5, double C0 = 0.5);

  void ComputeDatabaseInputs_new(const double rho, const double Z, const double Zv, const double Yc, const double rhoE,
                             double &P, double &Sz ,  double &C,
                             double P0=1.0e5, double C0 = 0.5);

  void ComputeDatabaseInputs_FromT(const double rho, const double Z, const double Zv, const double Yc, const double T,
                             double &P, double &Sz ,  double &C,
                             double P0=1.0e5, double C0 = 0.5);
  //-----------------------------------------
  /*! \brief compute Yc and density from C and pressure
   * used for Newton/Broyden iteration to find C and P
   */
  VecDoub YcAndDensityFromCP(VecDoub_I vec);
  VecDoub DensityFromP(VecDoub_I vec);
  VecDoub YcAndDensityFromCP_new(VecDoub_I vec);
  VecDoub DensityFromP_new(VecDoub_I vec);
  VecDoub YcAndDensityFromCPT(VecDoub_I vec);
  VecDoub DensityFromPT(VecDoub_I vec);

  //-----------------------------------------
  /*! \brief MPI database reading 
   */
  void BCAST_pressure(int mpi_root, MPI_Comm mpi_comm);
  void MPI_Load(string inputfilename);

  // Pointer to database to be used in the Newton user-function wrapper
  static PressureTable* ptr_for_newton;
private:
  string    TableFilename;
  double    version;
  string    CombustionModel;
  string    ChemTableType;
  vector <string> VarNames;
  map<string,int> VarNamesMap;

  // Parameters to compute C and P through Newton iteration
  //  They are declared here because the user-defined  function
  //  does not take other inputs than C and P
  double Z_newton;
  double Zv_newton;
  double rho_newton;
  double Yc_newton;
  double C_newton;
  double P_newton;
  double E_newton;
  double T_newton;
  bool verbose_newton;
  static double P_scale; // scale pressure for better convergence in Newton iteration

};

VecDoub usrfunc_wrapper(VecDoub_I x);
VecDoub usrfunc_wrapper2(VecDoub_I x);
VecDoub usrfunc_wrapper3(VecDoub_I x);
VecDoub usrfunc_wrapper4(VecDoub_I x);

#endif

