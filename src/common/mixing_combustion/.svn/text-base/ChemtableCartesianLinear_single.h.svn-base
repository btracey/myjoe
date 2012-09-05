/**
 * \brief Chemtable class for chemistry table: load table and lookup value based on 3 
 * coordinates using linear interpolation.
 * 
 *  This version of the chemistry table is based on a 3D structured Cartesian tabulation,
 *  i.e., the 3 dimensions are simple vectors. Adaptive is not implemented!
 *  The interpolation is linear.
 * 
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */

#ifndef CHEMTABLECARTESIANLINEAR_SINGLE_H
#define CHEMTABLECARTESIANLINEAR_SINGLE_H

#include "CombustionGeneralDefinitions.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;


#define Ver 2.0

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Class containing the chemistry table and associated functions.
 * 
 *  Read and load 3D table from chemistry table input file.
 *  Lookup values based on 3 dimensions by interpolating structured Cartesian table
 *  using linear interpolation 
 *  
 *  Class assumes a structured Cartesian table: x1=x1(i), x2=x2(j), x3=x3(k).
 */
class ChemtableCartesianLinear_single
{
private:
  double             version;                          ///< Version of the table interpolation routines (to check with table version)
  int                n1, n2, n3;                       ///< Dimensions of the table (corresponding to Zmean, Zvar and chi/progress variable)
  int                nvar;                             ///< Number of variables stored in table
  string             CombustionModel;                  ///< Combustion model used: Steady Flamelet or FPVA
  string             ChemTableType;                    ///< Type of chemistry table ("SPECIES" or "COEFF")
  double             Preference;                       ///< Reference pressure at which table was computed
  double             Toxi;                             ///< Temperature boundary condition for oxidizer (Z=0) with which table was computed
  double             Tfuel;                            ///< Temperature boundary condition for fuel (Z=1) with which table was computed
  vector <string>    VarNames;                         ///< Names of the stored variables
  double             *x1, *x2, *x3;                    ///< Axis coordinates of the table (Zmean, Zvar, chi/progress variable)
  float             ****Data;                         ///< Actual table data: Data(i,j,k,l), where i->x1, j->x2, k->x3, l->variable
  double             Zm, Zv, Cm;                       ///< Coordinates of point in table to interpolate (saved to avoid duplicate work)
  InterpolationIndex *interp;                          ///< Interpolation indices and weights (set by SetCoordinates once to avoid recomputing it for each species).
  
  map<string, int>                       myTableVar;   ///< Map containing the variables linking the name (string) to table last index l (int)
  map<string, int>::iterator             itTp;         ///< Iterators to loop over the variables stored in table
  pair<map<string, int>::iterator, bool> retT;         ///< Return value for map insertion to check if variable is already in map
  
public:  
  ChemtableCartesianLinear_single();
  
  /* Accessors */
  /*************/
  
  int GetChemtableDimension1() {return n1;}
  int GetChemtableDimension2() {return n2;}
  int GetChemtableDimension3() {return n3;}
  int GetChemtableDimension4() {return nvar;}
  
  double GetReferencePressure() {return Preference;}
  double GetToxi()  {return Toxi;}
  double GetTfuel() {return Tfuel;}
  
  string GetCombustionModel() {return CombustionModel;}
  string GetChemtableType()   {return ChemTableType;}
  string GetChemtableVariables(int i) {return VarNames[i];}

  double GetChemtableCoordinate1(int i1, int i2, int i3) {return x1[i1];}
  double GetChemtableCoordinate2(int i1, int i2, int i3) {return x2[i2];}
  double GetChemtableCoordinate3(int i1, int i2, int i3) {return x3[i3];}
  double GetChemtableCoordinate1(int i1) {return x1[i1];}
  double GetChemtableCoordinate2(int i2) {return x2[i2];}
  double GetChemtableCoordinate3(int i3) {return x3[i3];}
  
  int    GetChemtableIndex1(double A1, double A2, double A3) {int ii; BinarySearch(ii, x1, 1, n1, A1); return ii;}
  int    GetChemtableIndex2(double A1, double A2, double A3) {int jj; BinarySearch(jj, x2, 1, n2, A2); return jj;}
  int    GetChemtableIndex3(double A1, double A2, double A3) {int kk; BinarySearch(kk, x3, 1, n3, A3); return kk;}
  int    GetChemtableIndex1(double A1) {int ii; BinarySearch(ii, x1, 1, n1, A1); return ii;}
  int    GetChemtableIndex2(double A2) {int jj; BinarySearch(jj, x2, 1, n2, A2); return jj;}
  int    GetChemtableIndex3(double A3) {int kk; BinarySearch(kk, x3, 1, n3, A3); return kk;}

  float GetChemtableValue(int i1, int i2, int i3, int ivar) {return Data[i1][i2][i3][ivar];}
  float GetChemtableValue(int i1, int i2, int i3, string var) {return Data[i1][i2][i3][myTableVar.find(var)->second];}
  
  int    GetInterpolationIndex(int i) {return interp->index[i];}
  double GetInterpolationWeight(int i, int j) {return interp->weight[i][j];}

  /***********************************************************************************************************/

  /*! \brief Open, read and load chemistry table.
   *
   *  Read input file to determine table file name and some parameters, then read table from file and allocate 
   *  large array to contain table. Table is a 4D array with each species contiguous in memory.
   *  Also output information on table for verification.
   *  
   *  Table might use blanking if created with \p CreateChemtable tool from NGA. Blanking are dummy values
   *  if table is created with tool \p CreateChemTable in the tool directory.
   *  
   *  \param myInputPtr A pointer of type ParamMap pointing on the general input file.
   */
  void Load(string tablename);

  /***********************************************************************************************************/

  /*! \brief Clean table and deallocate memory */
  void Unload();

  /***********************************************************************************************************/

  /*! \brief Pad table for cubic interpolation. */
  void PadTable();
 
  /***********************************************************************************************************/  
  
  /*! \brief Compute index and interpolation weights based on 3 input parameters.
   *
   *  Computes lower indices of cell containing the 3 input parameters (usually: Zmean, Zvar and chi/progress variable)
   *  and interpolation weights to compute linear/cubic interpolation. If coordinates (input) are out of table bounds,
   *  the coordinates are projected to boundary of table: if x < x(pad), then x=x(pad).
   *  Data (indices + weights) is saved in variable #interp.
   *  \param[in]  A1  Value along first dimension (Zmean).
   *  \param[in]  A2  Value along second dimension (Zvar).
   *  \param[in]  A3  Value along third dimension (chi/progress variable).
   */
  void SetCoordinates(double A1, double A2, double A3);
 
  /***********************************************************************************************************/  
  
  /*! \brief Compute weights for cubic interpolation. */
  void Compute_Weights(double x, double x0, double xp1, double *wi);

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table based on indices, weights and name of variable.
   *
   *  Interpolation is linear. Any variable contained in table can be looked up. Indices and weights must be
   *  first computed using Chemtable::SetCoordinates . Returns the interpolated value of the variable.
   *  \param[in]  A1, A2, A3                             Coordinates.
   *  \param[in] tag                                     Name of variable to lookup.
   *  \return    Interpolated value of variable \p tag.
   */
  double Lookup(double A1, double A2, double A3, string tag);

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table to retrieve all coefficients.
   *
   *  Interpolation is linear. All coefficients contained in table are looked up. Indices and weights are
   *  first computed using Chemtable::SetCoordinates. Returns the interpolated value of the coefficients.
   *  \param[in]  A1, A2, A3                                                           Coordinates.
   *  \param[out] RoM, T0, E0, Gam0, a_Gam, mu0, a_mu, lamOcp0, a_lamOcp, src_prog     Interpolated value of all coefficients.
   */
  void LookupCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double &mu0, double &a_mu,
                   double &lamOcp0, double &a_lamOcp, double &src_prog, double A1, double A2, double A3);

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table to retrieve selected coefficients.
   *
   *  Interpolation is linear. Selected coefficients contained in table are looked up. Indices and weights are
   *  first computed using Chemtable::SetCoordinates. Returns the interpolated value of the coefficients.
   *  \param[in]  A1, A2, A3                   Coordinates.
   *  \param[out] RoM, T0, E0, Gam0, a_Gam     Interpolated value of all coefficients.
   */
  void LookupSelectedCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double A1, double A2, double A3);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

