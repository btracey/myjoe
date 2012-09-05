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

#ifndef CHEMISTRYTABLECARTESIANLINEAR_H
#define CHEMISTRYTABLECARTESIANLINEAR_H

#include "CombustionGeneralDefinitions.h"

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
class ChemtableCartesianLinear
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
  vector<string>     VarNames;                         ///< Names of the stored variables
  double             *x1, *x2, *x3;                    ///< Axis coordinates of the table (Zmean, Zvar, chi/progress variable)
  double             ****Data;                         ///< Actual table data: Data(i,j,k,l), where i->x1, j->x2, k->x3, l->variable
  double             Zm, Zv, Cm;                       ///< Coordinates of point in table to interpolate (saved to avoid duplicate work)
  InterpolationIndex *interp;                          ///< Interpolation indices and weights (set by SetCoordinates once to avoid recomputing it for each species).
  
  map<string, int>                       myTableVar;   ///< Map containing the variables linking the name (string) to table last index l (int)
  map<string, int>::iterator             itTp;         ///< Iterators to loop over the variables stored in table
  pair<map<string, int>::iterator, bool> retT;         ///< Return value for map insertion to check if variable is already in map
  
public:  
  ChemtableCartesianLinear()
  {
    if (mpi_rank == 0)
      cout << "ChemtableCartesianLinear()" << endl;
    
    x1   = NULL;
    x2   = NULL;
    x3   = NULL;
    Data = NULL;
    
    Zm = -1.0;
    Zv = -1.0;
    Cm = -1.0;
    
    version = Ver;

    interp = new InterpolationIndex(3,2);
  }
  
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

  double GetChemtableValue(int i1, int i2, int i3, int ivar) {return Data[i1][i2][i3][ivar];}
  double GetChemtableValue(int i1, int i2, int i3, string var) {return Data[i1][i2][i3][myTableVar.find(var)->second];}
  
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
  void Load(string tablename)
  {
    FILE *inFp = NULL;
    int i, j, k, l, mi1, mi2, mi3, ma1, ma2, ma3;
    double minval, maxval;
    string key;
    char filename[kStrLong], buffer_m[kStrMed];
    char *p = 0, dummy;
    size_t dum;

    /* Read input file for chemistry table information */
    dum = tablename.copy(&filename[0], kStrLong);
    dum = tablename.size();
    if (dum < kStrLong)
      strcpy(&filename[dum], "\0");

    /* Open chemistry table file */
    if (mpi_rank == 0)
      cout << endl << "Chemtable is " << filename << endl;
    
    if (!(inFp = fopen(filename, "rb")))
    {
      cerr << "### Could not open input file " << filename << " ###" << endl;
      throw(-1);
    }
    
    // *****************
    // Read general info
    // *****************
    double tableVersion;
    dum = fread(&tableVersion, sizeof(double), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: tableVersion ###" << endl;
      throw(-1);
    }
    if (tableVersion != version)
    {
      cerr << "### Error with version of chemistry table! Table was created with version " << tableVersion 
           << " of CreateChemTable and cannot be read with version " << version << " of interpolation function! ###" << endl;
      throw(-1);
    }

    for (i=0; i<kStrMed-1; i++)
      strcpy(&buffer_m[i], " ");
    dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
    if (dum != kStrMed)
    {
      cerr << "### Error reading file: Combustion model ###" << endl;
      throw(-1);
    }
    strcpy(&buffer_m[kStrMed - 1], "\0");
    CombustionModel.assign(buffer_m);

    for (i=0; i<kStrMed-1; i++)
      strcpy(&buffer_m[i], " ");
    dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
    if (dum != kStrMed)
    {
      cerr << "### Error reading file: Chemistry table type ###" << endl;
      throw(-1);
    }
    strcpy(&buffer_m[kStrMed - 1], "\0");
    ChemTableType.assign(buffer_m);
    
    dum = fread(&Preference, sizeof(double), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: Preference ###" << endl;
      throw(-1);
    }
    dum = fread(&Toxi, sizeof(double), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: Toxi ###" << endl;
      throw(-1);
    }
    dum = fread(&Tfuel, sizeof(double), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: Tfuel ###" << endl;
      throw(-1);
    }
      
    if (mpi_rank == 0)
    {
      cout << "Chemistry table type is: " << ChemTableType << endl;      
      cout << "The combustion model is: " << CombustionModel << endl;
      cout << "Reference pressure is:   " << Preference << endl;
      cout << "Toxidizer is:            " << Toxi << endl;
      cout << "Tfuel is:                " << Tfuel << endl;
    }
    
    
    // *********************
    // Read table dimensions
    // *********************
    dum = fread(&n1, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n1 ###" << endl;
      throw(-1);
    }
    dum = fread(&n2, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n2 ###" << endl;
      throw(-1);
    }
    dum = fread(&n3, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n3 ###" << endl;
      throw(-1);
    }
    dum = fread(&nvar, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: nvar ###" << endl;
      throw(-1);
    }
    if (mpi_rank == 0)
    {
      cout << "Chemistry table has size " << n1 << " x " << n2 << " x " << n3 << endl;
    }

    // **********************
    // Read table coordinates 
    // **********************
    getMem1D(&x1, 0, n1+1, "Chemtable::Load x1", true);          // padding: add 1 point at beginning and end of table for interpolation
    getMem1D(&x2, 0, n2+1, "Chemtable::Load x2", true);          // padding: add 1 point at beginning and end of table for interpolation
    getMem1D(&x3, 0, n3+1, "Chemtable::Load x3", true);          // padding: add 1 point at beginning and end of table for interpolation
    dum = fread(&x1[1], sizeof(double), n1, inFp);
    if (dum != n1)
    {
      cerr << "### Error reading file: x1 ###" << endl;
      throw(-1);
    }
    dum = fread(&x2[1], sizeof(double), n2, inFp);
    if (dum != n2)
    {
      cerr << "### Error reading file: x2 ###" << endl;
      throw(-1);
    }
    dum = fread(&x3[1], sizeof(double), n3, inFp);
    if (dum != n3)
    {
      cerr << "### Error reading file: x3 ###" << endl;
      throw(-1);
    }
    
    // temporary pad grid boundaries with 0's
    x1[0]    = 0.0;
    x2[0]    = 0.0;
    x3[0]    = 0.0;
    x1[n1+1] = 0.0;
    x2[n2+1] = 0.0;
    x3[n3+1] = 0.0;
    
    // output some information
    if (mpi_rank == 0)
    {
      cout.width(10);
      cout << "x1 = Zmean:      " << "\t min=" << MinVal(x1, n1+2, 1) << "\t\t max=" << MaxVal(x1, n1+2, 1) << endl;      
      cout.width(10);
      cout << "x2 = Zvar:       " << "\t min=" << MinVal(x2, n2+2, 1) << "\t\t max=" << MaxVal(x2, n2+2, 1) << endl;      
      cout.width(10);
      cout << "x3 = Prog/chi:   " << "\t min=" << MinVal(x3, n3+2, 1) << "\t\t max=" << MaxVal(x3, n3+2, 1) << endl;
      
      cout << "and contains " << nvar << " variables: " << endl;
    }


    // **********************
    // Read name of variables
    // **********************
    for (l = 0; l < nvar; l++)
    {
      for (i=0; i<kStrMed-1; i++)
        strcpy(&buffer_m[i], " ");
      dum = fread(&buffer_m[0], sizeof(char), kStrMed, inFp);
      if (dum != kStrMed)
      {
        cerr << "### Error reading file: VarNames no. " << l << " ###" << endl;
        throw(-1);
      }
      strcpy(&buffer_m[kStrMed - 1], "\0");
      VarNames.push_back(buffer_m);
      
      /* Create map to link variable names and index */
      key.assign(VarNames[l]);
      retT = myTableVar.insert(pair<string, int> (key, l));
      if (retT.second == false)
      {
        cerr << "### Variable " << key << " already in myTableVar in loop ###" << endl;
        throw(-1);
      }
    }

    // ***************
    // Read table data 
    // ***************
    getMem4D(&Data, 0, n1+1, 0, n2+1, 0, n3+1, 0, nvar-1, "Chemtable::Load Data", true);
    for (i=1; i<n1+1; i++)
      for (j=1; j<n2+1; j++)
        for (k=1; k<n3+1; k++)
          for (l=0; l<nvar; l++)
          {
            dum = fread(&Data[i][j][k][l], sizeof(double), 1, inFp);
            if (dum != 1)
            {
              cerr << "### Error reading file: data ->" << i << " " << j << " " << k << " " << l << " ### " << endl;
              throw(-1);
            }
          }

    fclose(inFp);

    
    // *********************************
    // Pad table for cubic interpolation
    // *********************************
    PadTable();
    
    
    // **********************************
    // Output useful information on table
    // **********************************
    if (mpi_rank == 0)
    {
      for (l = 0; l < nvar; l++)
      {
        GetArrayInfo(Data, l, n1+2, n2+2, n3+2, minval, maxval, mi1, mi2, mi3, ma1, ma2, ma3, 1);
        cout.width(3);
        cout << l << ". ";
        cout.width(15);
        cout << VarNames[l];
        cout << "\t min=" << minval << "\t at i=" << mi1 << "  j=" << mi2 << "  k=" << mi3;
        cout << "\t max=" << maxval << "\t at i=" << ma1 << "  j=" << ma2 << "  k=" << ma3 << endl << endl;
      }
    }
  }

  /***********************************************************************************************************/

  /*! \brief Clean table and deallocate memory */
  void Unload()
  {
    if (x1 != NULL)       {freeMem1D(x1, 0, n1+1);                                  x1 = NULL;}
    if (x2 != NULL)       {freeMem1D(x2, 0, n2+1);                                  x2 = NULL;}
    if (x3 != NULL)       {freeMem1D(x3, 0, n3+1);                                  x3 = NULL;}
    if (Data != NULL)     {freeMem4D(Data, 0, n1+1, 0, n2+1, 0, n3+1, 0, nvar-1);   Data = NULL;}
  }

  /***********************************************************************************************************/

  /*! \brief Pad table for cubic interpolation. */
  void PadTable()
  // Copied and adapted from CDP chemtable_m.f90
  {
    double wextrap_1[3][3], wextrap_n[3][3];
    double xi, x0, xp1, xp2;

    // pad grid boundaries
    x1[0]    = 2.0 * x1[1]  - x1[2];
    x2[0]    = 2.0 * x2[1]  - x2[2];
    x3[0]    = 2.0 * x3[1]  - x3[2];
    x1[n1+1] = 2.0 * x1[n1] - x1[n1-1];
    x2[n2+1] = 2.0 * x2[n2] - x2[n2-1];
    x3[n3+1] = 2.0 * x3[n3] - x3[n3-1];

    // pad all data
    // x1 dimension
    for (int j=1; j<n2+1; j++)
      for (int k=1; k<n3+1; k++)
        for (int l=0; l<nvar; l++)
        {
          Data[0][j][k][l]    = 2.0 * Data[1][j][k][l]  - Data[2][j][k][l];
          Data[n1+1][j][k][l] = 2.0 * Data[n1][j][k][l] - Data[n1-1][j][k][l];
        }                                                                    
    // x2 dimension
    for (int i=1; i<n1+1; i++)
      for (int k=1; k<n3+1; k++)
        for (int l=0; l<nvar; l++)
        {
          Data[i][0][k][l]    = 2.0 * Data[i][1][k][l]  - Data[i][2][k][l];
          Data[i][n2+1][k][l] = 2.0 * Data[i][n2][k][l] - Data[i][n2-1][k][l];
        }     
    // x3 dimension
    for (int i=1; i<n1+1; i++)
      for (int j=1; j<n2+1; j++)
        for (int l=0; l<nvar; l++)
        {
          Data[i][j][0][l]    = 2.0 * Data[i][j][1][l]  - Data[i][j][2][l];
          Data[i][j][n3+1][l] = 2.0 * Data[i][j][n3][l] - Data[i][j][n3-1][l];
        }

  }
  
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
  void SetCoordinates(double A1, double A2, double A3)
  // Find index ii such that x1[ii] <= A1 < x1[ii+1]
  // Compute weights
  {
    double x0, xp1;
    
    // First direction
    BinarySearch(interp->index[0], x1, 1, n1, A1);
    
    x0  = x1[interp->index[0]];
    xp1 = x1[interp->index[0] + 1];
    
    Compute_Weights(A1, x0, xp1, interp->weight[0]);
    
    // Second direction
    BinarySearch(interp->index[1], x2, 1, n2, A2);
    
    x0  = x2[interp->index[1]];
    xp1 = x2[interp->index[1] + 1];
    
    Compute_Weights(A2, x0, xp1, interp->weight[1]);
    
    // Third direction
    BinarySearch(interp->index[2], x3, 1, n3, A3);
    
    x0  = x3[interp->index[2]];
    xp1 = x3[interp->index[2] + 1];
    
    Compute_Weights(A3, x0, xp1, interp->weight[2]);

  }
  
  /***********************************************************************************************************/  
  
  /*! \brief Compute weights for cubic interpolation. */
  void Compute_Weights(double x, double x0, double xp1, double *wi)
  {
    double xp1x0;
    
    // short-hand notation
    xp1x0  = xp1 - x0;
    
    // weight for f(i)
    wi[0] = (xp1 - x) / xp1x0;
    
    // weight for f(i+1)
    wi[1] = (x - x0) / xp1x0;
   
  }

  int GetVarNameIndex( string tag ) {
    int ll;
    // Determine the index of the variable to lookup
    itTp = myTableVar.find(tag);
    if (itTp == myTableVar.end()) {
      cerr << "### Unknown variable in chemistry table (Lookup): " << tag << " ###" << endl;
      throw(-1);
    }
    return ll = itTp->second;
  }

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table based on indices, weights and name of variable.
   *
   *  Interpolation is linear. Any variable contained in table can be looked up. Indices and weights must be
   *  first computed using Chemtable::SetCoordinates . Returns the interpolated value of the variable.
   *  \param[in]  A1, A2, A3                             Coordinates.
   *  \param[in] tag                                     Name of variable to lookup.
   *  \return    Interpolated value of variable \p tag.
   */
  double Lookup(double A1, double A2, double A3, string tag)
  {
    int ll;
    double val, g[4], f[4];

    // Clip coordinates
    A1 = min(max(A1, x1[1]), x1[n1]);
    A2 = min(max(A2, x2[1]), x2[n2]);
    A3 = min(max(A3, x3[1]), x3[n3]);
    
    // Check if coordinates and weights need to be updated
    if ((Zm != A1) || (Zv != A2) || (Cm != A3))
    {
      Zm = A1;
      Zv = A2;
      Cm = A3;
      SetCoordinates(Zm, Zv, Cm);
    }
    
    // Determine the index of the variable to lookup
    ll = GetVarNameIndex(tag);

    // Interpolate value based on index and weights for given variable
    val = 0.0;    
    for (int k = 0; k < 2; k++)
    {
      g[k] = 0.0;
      for (int j = 0; j < 2; j++)
      {
        f[j] = 0.0;
        for (int i = 0; i < 2; i++)
          f[j] += interp->weight[0][i] * Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
        g[k] += interp->weight[1][j] * f[j];
      }
      val += interp->weight[2][k] * g[k];
    }
    
    return val;
  }

 /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table based on indices, weights and index of variable.
   *
   *  Interpolation is linear. Any variable contained in table can be looked up. Indices and weights must be
   *  first computed using Chemtable::SetCoordinates . Returns the interpolated value of the variable.
   *  \param[in]  A1, A2, A3                             Coordinates.
   *  \param[in] tag                                     Name of variable to lookup.
   *  \return    Interpolated value of variable \p tag.
   */
  double Lookup(double A1, double A2, double A3, int ll)
  {
    double val, g[4], f[4];

    // Clip coordinates
    A1 = min(max(A1, x1[1]), x1[n1]);
    A2 = min(max(A2, x2[1]), x2[n2]);
    A3 = min(max(A3, x3[1]), x3[n3]);

    // Check if coordinates and weights need to be updated
    if ((Zm != A1) || (Zv != A2) || (Cm != A3))
    {
      Zm = A1;
      Zv = A2;
      Cm = A3;
      SetCoordinates(Zm, Zv, Cm);
    }

    // Interpolate value based on index and weights for given variable
    val = 0.0;
    for (int k = 0; k < 2; k++)
    {
      g[k] = 0.0;
      for (int j = 0; j < 2; j++)
      {
        f[j] = 0.0;
        for (int i = 0; i < 2; i++)
          f[j] += interp->weight[0][i] * Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
        g[k] += interp->weight[1][j] * f[j];
      }
      val += interp->weight[2][k] * g[k];
    }

    return val;
  }


  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table to retrieve all coefficients.
   *
   *  Interpolation is linear. All coefficients contained in table are looked up. Indices and weights are
   *  first computed using Chemtable::SetCoordinates. Returns the interpolated value of the coefficients.
   *  \param[in]  A1, A2, A3                                                           Coordinates.
   *  \param[out] RoM, T0, E0, Gam0, a_Gam, mu0, a_mu, lamOcp0, a_lamOcp, src_prog     Interpolated value of all coefficients.
   */
  void LookupCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double &mu0, double &a_mu,
                   double &lamOcp0, double &a_lamOcp, double &src_prog, double A1, double A2, double A3)
  {
    int ll;
    double coeff[10], g[4], f[4];

    // Clip coordinates
    A1 = min(max(A1, x1[1]), x1[n1]);
    A2 = min(max(A2, x2[1]), x2[n2]);
    A3 = min(max(A3, x3[1]), x3[n3]);
    
    // Check if coordinates and weights need to be updated
    if ((Zm != A1) || (Zv != A2) || (Cm != A3))
    {
      Zm = A1;
      Zv = A2;
      Cm = A3;
      SetCoordinates(Zm, Zv, Cm);
    }
    
    for (int ll = 0; ll < 10; ll++)
    {
      // Interpolate value based on index and weights for given variable
      coeff[ll] = 0.0;    
      for (int k = 0; k < 2; k++)
      {
        g[k] = 0.0;
        for (int j = 0; j < 2; j++)
        {
          f[j] = 0.0;
          for (int i = 0; i < 2; i++)
            f[j] += interp->weight[0][i] * Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
          g[k] += interp->weight[1][j] * f[j];
        }
        coeff[ll] += interp->weight[2][k] * g[k];
      }
    }
    
    RoM      = coeff[0];
    T0       = coeff[1];
    E0       = coeff[2];
    Gam0     = coeff[3];
    a_Gam    = coeff[4];
    mu0      = coeff[5];
    a_mu     = coeff[6];
    lamOcp0  = coeff[7];
    a_lamOcp = coeff[8];
    src_prog = coeff[9];
  }

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table to retrieve selected coefficients.
   *
   *  Interpolation is linear. Selected coefficients contained in table are looked up. Indices and weights are
   *  first computed using Chemtable::SetCoordinates. Returns the interpolated value of the coefficients.
   *  \param[in]  A1, A2, A3                   Coordinates.
   *  \param[out] RoM, T0, E0, Gam0, a_Gam     Interpolated value of all coefficients.
   */
  void LookupSelectedCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double A1, double A2, double A3)
  {
    int ll;
    double coeff[5], g[4], f[4];

    // Clip coordinates
    A1 = min(max(A1, x1[1]), x1[n1]);
    A2 = min(max(A2, x2[1]), x2[n2]);
    A3 = min(max(A3, x3[1]), x3[n3]);
    
    // Check if coordinates and weights need to be updated
    if ((Zm != A1) || (Zv != A2) || (Cm != A3))
    {
      Zm = A1;
      Zv = A2;
      Cm = A3;
      SetCoordinates(Zm, Zv, Cm);
    }
    
    for (int ll = 0; ll < 5; ll++)
    {
      // Interpolate value based on index and weights for given variable
      coeff[ll] = 0.0;    
      for (int k = 0; k < 2; k++)
      {
        g[k] = 0.0;
        for (int j = 0; j < 2; j++)
        {
          f[j] = 0.0;
          for (int i = 0; i < 2; i++)
            f[j] += interp->weight[0][i] * Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
          g[k] += interp->weight[1][j] * f[j];
        }
        coeff[ll] += interp->weight[2][k] * g[k];
      }
    }
    
    RoM      = coeff[0];
    T0       = coeff[1];
    E0       = coeff[2];
    Gam0     = coeff[3];
    a_Gam    = coeff[4];
  }

/**********************************************************************************************/

  void WriteTableTecplot()
  /* Write table to file in Tecplot ASCII format for visualization */
  {
    ofstream fout;
    char     filename[kStrLong] = "database.dat";
    size_t   dum;
  
    cout << endl << " Write table to tecplot format file ... ";
  
    // *********
    // Open file
    // *********
    fout.open(filename);
  
    // ********************
    // Write tecplot header
    // ********************
    fout << "TITLE=\"Chemistry table\"" << endl;
    fout << "VARIABLES=\"I1\" \"I2\" \"I3\" \"ZMean\" \"ZVar\" \"Progress\" ";
    for (int j=0; j<nvar; j++)
      fout << "\"" << VarNames[j] << "\" ";
    fout << endl;
  
    fout << "ZONE I=" << n1 << ", J=" << n2 << ", K=" << n3 << ", DATAPACKING=POINT" << endl;
  
    // *******************************
    // Write data for each table point
    // *******************************
    for (int k=1; k<n3+1; k++)
      for (int j=1; j<n2+1; j++)
        for (int i=1; i<n1+1; i++)
        {
          double prog = x3[k];
          fout << i << "  " << j << "  " << k << "  ";
          fout << x1[i] << "   " << x2[j] << "   " << prog;
          for (int n=0; n<nvar; n++)
            fout << "   " << Data[i][j][k][n];
          fout << endl;
        }
  
    fout.close();
  
    cout << "done !" << endl;
  
    return;
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
};

#endif
