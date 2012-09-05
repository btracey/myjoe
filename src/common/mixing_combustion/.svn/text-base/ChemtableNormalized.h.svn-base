#ifndef CHEMTABLENORMALIZED_H
#define CHEMTABLENORMALIZED_H
#include "Table3D.h"
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;

/*! /brief Class PressureTable describes a multiple pressure flamelet database
 *
 */
class ChemtableNormalized : public Table3D
{
public:

  //------------------------------//
  //         Constructors         //
  //------------------------------//
  ChemtableNormalized();
  ChemtableNormalized(string filename, int nZm, int nZv, int nC, int nVar);
  //------------------------------//
  //         Destructor           //
  //------------------------------//
  ~ChemtableNormalized();
  //------------------------------//
  //     Get/Set Attributes       //
  //------------------------------//
  inline double GetTableVersion()  {return version;}
  inline string GetTableFileName() {return TableFilename;}
  inline string GetTableCombustionModel() {return CombustionModel;}
  inline string GetTableType() {return ChemTableType;}
  inline double GetReferencePressure() {return Preference;}
  inline string GetVarName(int i) {return VarNames[i];}
  inline int    GetVarNameIndex(string Var) {return VarNamesMap.find(Var)->second;}

  inline void SetTableVersion(double ver) {version = ver;}
  inline void SetTableFileName(string name) {TableFilename = name;}
  inline void SetTableCombustionModel(string model) {CombustionModel = model;}
  inline void SetTableType(string type) {ChemTableType = type;}
  inline void SetReferencePressure(double P) {Preference = P;}
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
  void WriteTableTecplot();
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
  /*! \brief Interpolate the quantity "variable" in the 3D database
   * Example: Lookup(mixture_value, mixture_variance_value, prog_value, "my_variable");
   */
  double Lookup(double Zm, double Zvar, double C, string variable);
  double Lookup(double Zm, double Zvar, double C, int ivar);
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
                   double &src_prog, double Zm,
                   double Zvar,      double C);
  //-----------------------------------------
  /*! \brief Interpolate the first 5 variables in the database
   */
  void LookupSelectedCoeff(double &RoM,  double &T0,    double &E0,
                           double &Gam0, double &a_Gam,
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
  double ComputeNormalizedProg(const double Zm, const double Sz, const double Yc);

  //-----------------------------------------
  /*! \brief MPI database reading 
   */
  void BCAST_table(int mpi_root, MPI_Comm mpi_comm);
  void MPI_Load(string inputfilename);

private:
  string    TableFilename;
  double    version;
  string    CombustionModel;
  string    ChemTableType;
  double    Preference; 
  vector <string> VarNames;
  map<string,int> VarNamesMap;
};

#endif

