#include "ChemtableCartesianLinear_single.h"

ChemtableCartesianLinear_single::ChemtableCartesianLinear_single()
  {
    //if (mpi_rank == 0)
      cout << "ChemtableCartesianLinear_single()" << endl;
    
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

void ChemtableCartesianLinear_single::Load(string tablename)
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
    //if (mpi_rank == 0)
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
      
    //if (mpi_rank == 0)
    //{
      cout << "Chemistry table type is: " << ChemTableType << endl;      
      cout << "The combustion model is: " << CombustionModel << endl;
      cout << "Reference pressure is:   " << Preference << endl;
      cout << "Toxidizer is:            " << Toxi << endl;
      cout << "Tfuel is:                " << Tfuel << endl;
    //}
    
    
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
    //if (mpi_rank == 0)
    //{
      cout << "Chemistry table has size " << n1 << " x " << n2 << " x " << n3 << endl;
    //}

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
    //if (mpi_rank == 0)
    //{
      cout.width(10);
      cout << "x1 = Zmean:      " << "\t min=" << MinVal(x1, n1+2, 1) << "\t\t max=" << MaxVal(x1, n1+2, 1) << endl;      
      cout.width(10);
      cout << "x2 = Zvar:       " << "\t min=" << MinVal(x2, n2+2, 1) << "\t\t max=" << MaxVal(x2, n2+2, 1) << endl;      
      cout.width(10);
      cout << "x3 = Prog/chi:   " << "\t min=" << MinVal(x3, n3+2, 1) << "\t\t max=" << MaxVal(x3, n3+2, 1) << endl;
      
      cout << "and contains " << nvar << " variables: " << endl;
    //}


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
            dum = fread(&Data[i][j][k][l], sizeof(float), 1, inFp);
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
    //if (mpi_rank == 0)
    //{
      for (l = 0; l < nvar; l++)
      {
        GetArrayInfo_float(Data, l, n1+2, n2+2, n3+2, minval, maxval, mi1, mi2, mi3, ma1, ma2, ma3, 1);
        cout.width(3);
        cout << l << ". ";
        cout.width(15);
        cout << VarNames[l];
        cout << "\t min=" << minval << "\t at i=" << mi1 << "  j=" << mi2 << "  k=" << mi3;
        cout << "\t max=" << maxval << "\t at i=" << ma1 << "  j=" << ma2 << "  k=" << ma3 << endl << endl;
      }
    //}
  }

  /***********************************************************************************************************/

  /*! \brief Clean table and deallocate memory */
  void ChemtableCartesianLinear_single::Unload()
  {
    if (x1 != NULL)       {freeMem1D(x1, 0, n1+1);                                  x1 = NULL;}
    if (x2 != NULL)       {freeMem1D(x2, 0, n2+1);                                  x2 = NULL;}
    if (x3 != NULL)       {freeMem1D(x3, 0, n3+1);                                  x3 = NULL;}
    if (Data != NULL)     {freeMem4D(Data, 0, n1+1, 0, n2+1, 0, n3+1, 0, nvar-1);   Data = NULL;}
  }

  /***********************************************************************************************************/

  /*! \brief Pad table for cubic interpolation. */
  void ChemtableCartesianLinear_single::PadTable()
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
  void ChemtableCartesianLinear_single::SetCoordinates(double A1, double A2, double A3)
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
  void ChemtableCartesianLinear_single::Compute_Weights(double x, double x0, double xp1, double *wi)
  {
    double xp1x0;
    
    // short-hand notation
    xp1x0  = xp1 - x0;
    
    // weight for f(i)
    wi[0] = (xp1 - x) / xp1x0;
    
    // weight for f(i+1)
    wi[1] = (x - x0) / xp1x0;
   
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
  double ChemtableCartesianLinear_single::Lookup(double A1, double A2, double A3, string tag)
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
    itTp = myTableVar.find(tag);
    if (itTp == myTableVar.end())
    {
      cerr << "### Unknown variable in chemistry table (Lookup): " << tag << " ###" << endl;
      throw(-1);
    }
    ll = itTp->second;

    // Interpolate value based on index and weights for given variable
    val = 0.0;    
    for (int k = 0; k < 2; k++)
    {
      g[k] = 0.0;
      for (int j = 0; j < 2; j++)
      {
        f[j] = 0.0;
        for (int i = 0; i < 2; i++)
          f[j] += interp->weight[0][i] * (double) Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
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
  void ChemtableCartesianLinear_single::LookupCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double &mu0, double &a_mu,
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
            f[j] += interp->weight[0][i] * (double) Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
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
  void ChemtableCartesianLinear_single::LookupSelectedCoeff(double &RoM, double &T0, double &E0, double &Gam0, double &a_Gam, double A1, double A2, double A3)
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
            f[j] += interp->weight[0][i] * (double) Data[interp->index[0]+i][interp->index[1]+j][interp->index[2]+k][ll];
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

