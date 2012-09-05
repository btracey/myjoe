/**
 * \brief Artificial Neural Network (ANN) class for chemistry table: load table, compute value
 * 
 * \author Vincent Terrapon 
 * \date August 2009
 * \version 1.0
 * 
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!! NOT IMPLEMENTED YET !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

#ifndef CHEMISTRYTABLE_H
#define CHEMISTRYTABLE_H

#include "CombustionGeneralDefinitions.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Classes definition */
/**********************/

/*! \brief Neuron represent a neuron in the ANN with a certain number of inputs and 1 output.
 *  A neuron computes the weighted sum of its inputs with some specific weight and propagate this sum
 *  through a transfer function to give 1 output. The input signals come usually from the output signal of 
 *  other neurons upstream in the network or from the overall input signal.
 *  y = F(sum_i x_i * w_i + w_0)
 *  Neuron is the base class and NonlinearNeuron, HeadNeuron are derived classes.
 */
class Neuron
{
protected:
  int    myId;        ///< Id of the neuron.
  int    nInputs;     ///< The number of input signals (= number of connections with upstream neurons).
  double Output;      ///< The output signal y.
  double a_cst;       ///< Constant used by sigmoidal function (transfer function).
  double b_cst;       ///< Constant used by sigmoidal function (transfer function).
  
public:
  Neuron()  {}
  ~Neuron() {}
  
  void   SetmyId(int i) {myId = i;}
  int    GetmyId() {return myId;}
  void   SetNumberInputs(int i) {nInputs = i;}
  int    GetNumberInputs() {return nInputs;}
  void   SetTranserFunctionCst(double a, double b) {a_cst = a; b_cst = b;}
};

/*! \brief NonLinearNeuron is derived from Neuron for the neurons between head and tail neurons.
 *  The Inputs point to upstream neurons.
 */
class NonlinearNeuron: public Neuron
{
protected:
  Neuron **Inputs;     ///< The array of neurons sending the input signals x_i.
  double *Weights;     ///< The weight w_i of the different input signals.
  double Threshold;    ///< The threshold value w_0.
  
public:
  NonlinearNeuron()  {Inputs = NULL; Weights = NULL;}
  ~NonlinearNeuron() {CleanNonlinearNeuron;}
  
  void SetThreshold(double thresh) {Threshold = thresh;}
  void SetWeights(double *wi, int nwi)
  {
    if (nwi != nInputs)
    {
      cerr << "Number of weights does not match number of inputs in NonlinearNeuron::SetWeights()" << endl;
      throw(-1);
    }
    if (Weights == NULL) Weights = new double(nInputs);
    for (int i=0; i<nInputs; i++)
      Weights[i] = wi[i];
  }
  
  void SetInputs(Neuron **link, int nlink)
  {
    if (nlink != nInputs)
    {
      cerr << "Number of upstream neurons does not match number of inputs in NonlinearNeuron::SetInputs()" << endl;
      throw(-1);
    }
    if (Inputs == NULL) Inputs = new Neuron *[nInputs];
    for (int i=0; i<nInputs; i++)
      Inputs[i] = link[i];
  }
  
  void ComputeOutput()
  {
    // Weighted sum
    double sum = -1.0 * Threshold;
    for (int i=0; i<nInputs; i++)
    {
      sum += Weights[i] * Inputs[i]->GetOutput();
    }
    
    // Sigmoidal function
    return 2.0 * a_cst / (1.0 + exp(-b_cst * sum)) - a_cst;
  }
  
  double GetOutput() {return Output;} 
  
  void CleanNonlinearNeuron()
  {
    if (Inputs != NULL) delete [] Inputs;
    if (Weights != NULL) delete [] Weights;
  }
};

/*! \brief Tailneuron is the same as NonlinearNeuron but the transfer function is the identity function.
 *  y = sum_i x_i * w_i + w_0
 *  The output signal is also rescaled
 */
class TailNeuron: public NonlinearNeuron
{
protected:
  double ymin;    ///< Minimum value of signal
  double ymax;    ///< Maximum value of signal
  
public:
  TailNeuron()  {}
  ~TailNeuron() {}
  
  void SetSignalMinMax(double ymi, double yma) {ymin = ymi; ymax = yma;}

  void ComputeOutput()
  {
    double sum = Threshold;
    for (int i=0; i<nInputs; i++)
    {
      sum += Weights[i] * Inputs[i].GetOutput();
    }
    
    // No transfer function but unscaling of signal
    Output = 0.5 * (sum + 1.0) * (ymax - ymin) + ymin;
  } 
  double GetOutput() {return Output;}    // For TailNeuron output is unscaled signal
};

class HeadNeuron: public Neuron
{
protected:
  Signal *Inputs;     ///< Input signal for the entire ANN (here a signal and not a neuron).
  
public:
  HeadNeuron()  {Inputs = NULL;}
  ~HeadNeuron() {CleanHeadNeuron;}
  
  void SetInputs(Signal *link, int nlink)
  {
    if (nlink != nInputs)
    {
      cerr << "Number of upstream neurons does not match number of inputs in NonlinearNeuron::SetInputs()" << endl;
      throw(-1);
    }
    if (Inputs == NULL) Inputs = new Neuron *[nInputs];
    for (int i=0; i<nInputs; i++)
      Inputs[i] = link[i];
  }
  
  void ComputeOutput()
  {
    // Weighted sum
    double sum = Threshold;
    for (int i=0; i<nInputs; i++)
    {
      sum += Weights[i] * Inputs[i]->GetOutput();
    }
    
    // Sigmoidal function
    return 2.0 * a_cst / (1.0 + exp(-b_cst * sum)) - a_cst;
  }
  
  double GetOutput() {return Output.GetScaledSignal();}    // For NonlinearNeuron output is scaled signal
  
  void CleanHeadNeuron()
  {
    if (Inputs != NULL) delete [] Inputs;
  }
};

/***********************************************************************************************************/

/*! \brief Layer base class; a layer contains a certain number of neurons.
 *  The neurons within a layer act in parallel while different layers act in a serial way.
 */
class Layer
{
protected:
  int    myId;
  int    nNeurons;     ///< The number of input signals (= number of connections with upstream neurons).
  
public:
  Layer()  {}
  ~Layer() {}

  void SetmyId(int i) {myId = i;}
  int  GetmyId() {return myId;}
  void SetNumberNeurons(int i) {nNeurons = i;}
  int  GetNumberNeurons() {return nNeurons;}
  
};

/*! \brief Layers containing the nonlinear neurons, i.e., all layers between head layer and tail layer.
 */
class NonlinearLayer
{
  
};

/***********************************************************************************************************/

/*! \brief Class containing the chemistry table and associated functions.
 * 
 *  Load table and interpolate chemistry table. Read table from chemistry table input file.
 *  
 *  Base class assumes a structured Cartesian table: x1=x1(i), x2=x2(j), x3=x3(k).
 */
class Chemtable
{
private:
  int n1, n2, n3;                      ///< Dimensions of the table (corresponding to Zmean, Zvar and chi/progress variable)
  int nvar;                            ///< Number of variables stored in table
  char CombustionModel[kStrMed];       ///< Combustion model used: Steady Flamelet, FPVA or Enthalpy Flamelet
  char (*VarNames)[kStrMed];           ///< Names of the stored variables
  double *x1, *x2, *x3;                ///< Axis coordinates of the table (Zmean, Zvar, chi/progress variable)
  int ***iMask;                        ///< Indicating cells of table with no physical meaning: iMask(k,j,i) where k->x3, j->x2, i->x1
  bool is_iMask;                       ///< True if iMask (blanked cells) should be loaded
  double ****Data;                     ///< Actual table data: Data(l,k,j,i) where l->variable, k->x3, j->x2, i->x1

public:
  map<string, int> myTableVar;                   ///< Map containing the species linking the name (string) to table first index (int)
  map<string, int>::iterator itTp;               ///< Iterators to loop over the species
  pair<map<string, int>::iterator, bool> retT;   ///< Return value for map insertion to check if species is already in map
  
  /* Accessors */
  /*************/
  
  int GetChemtableDimension1() {return n1;}
  int GetChemtableDimension2() {return n2;}
  int GetChemtableDimension3() {return n3;}
  int GetChemtableDimension4() {return nvar;}

  double GetChemtableCoordinate1(int i1) {return x1[i1];}
  double GetChemtableCoordinate2(int i2) {return x2[i2];}
  double GetChemtableCoordinate3(int i3) {return x3[i3];}

  double GetChemtableValue(int i1, int i2, int i3, int ivar) {return Data[ivar][i3][i2][i1];}


  /***********************************************************************************************************/

  /*! \brief Open, read and load chemistry table.
   *
   *  Chemistry table is created with \p CreateChemtable tool from NGA. Read input file to determine table 
   *  file name and some parameters, then read table from file and allocate large array to contain table. 
   *  Table is a 4D array with each species contiguous in memory.
   *  Also output information on table for verification.
   *  \param myInputPtr A pointer of type ParamMap pointing on the general input file.
   */
  void Load(ParamMap *myInputPtr)
  {
    FILE *inFp = 0;
    int i, j, k, l, mi1, mi2, mi3, ma1, ma2, ma3;
    string key, buffer_str;
    char filename[kStrLong], buffer_m[kStrMed];
    char *p = 0, dummy;
    bool is_iMask;
    size_t dum;
    double minval, maxval;

    /* Read input file for chemistry table information */
    buffer_str = myInputPtr->getStringParam("CHEMTABLE_FILE");
    dum = buffer_str.copy(&filename[0], kStrLong);
    dum = buffer_str.size();
    if (dum < kStrLong)
      strcpy(&filename[dum], "\0");
    is_iMask = (bool) myInputPtr->getIntParam("LOAD_MASK");

    /* Read chemistry table file */
    if (mpi_rank == 0)
    {
      cout << endl << "Chemtable is " << filename << endl;
    }
    if (!(inFp = fopen(filename, "rb")))
    {
      cerr << "### Could not open input file " << filename << " ###\n";
      throw(-1);
    }
    /* Dimensions */
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

    /* Coordinates */
    getMem1D(&x1, 0, n1 - 1, "Chemtable::Load x1", true);
    getMem1D(&x2, 0, n2 - 1, "Chemtable::Load x2", true);
    getMem1D(&x3, 0, n3 - 1, "Chemtable::Load x3", true);
    dum = fread(x1, sizeof(double), n1, inFp);
    if (dum != n1)
    {
      cerr << "### Error reading file: x1 ###" << endl;
      throw(-1);
    }
    dum = fread(x2, sizeof(double), n2, inFp);
    if (dum != n2)
    {
      cerr << "### Error reading file: x2 ###" << endl;
      throw(-1);
    }
    dum = fread(x3, sizeof(double), n3, inFp);
    if (dum != n3)
    {
      cerr << "### Error reading file: x3 ###" << endl;
      throw(-1);
    }

    /* Masks */
    getMem3D(&iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1, "Chemtable::Load iMask", false);
    for (k = 0; k < n3; k++)
    {
      for (j = 0; j < n2; j++)
      {
        for (i = 0; i < n1; i++)
        {
          dum = fread(&iMask[k][j][i], sizeof(int), 1, inFp);
          if (dum != 1)
          {
            cerr << "### Error reading file: mask (" << i << j << k
            << ") ###" << endl;
            throw(-1);
          }
        }
      }
    }
    if (!(is_iMask))
    {
      freeMem3D(iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);
    } /* iMask not used, save memory */

    /* Other input data */
    for (i = 0; i < kStrMed - 1; i++)
    {
      strcpy(&CombustionModel[i], " ");
    }
    dum = fread(&CombustionModel[0], sizeof(char), kStrMed, inFp);
    if (dum != kStrMed)
    {
      cerr << "### Error reading file: Combustion model ###" << endl;
      throw(-1);
    }
    strcpy(&CombustionModel[kStrMed - 1], "\0");

    VarNames = new char[nvar][kStrMed];
    for (l = 0; l < nvar; l++)
    {
      strcpy(&VarNames[l][0], "\0");
      dum = fread(&buffer_m, sizeof(char), kStrMed, inFp);
      if (dum != kStrMed)
      {
        cerr << "### Error reading file: VarNames no. " << l << " ###"
        << endl;
        throw(-1);
      }
      p = strchr(buffer_m, ' ');
      if (p != 0)
        strcpy(p, "\0");
      strcpy(&VarNames[l][0], buffer_m);
      /* Create map to link variable names and index */
      key.assign(buffer_m);
      retT = myTableVar.insert(pair<string, int> (string(key), l));
      if (retT.second == false)
      {
        cout << "### Variable " << key
        << " already in myTableVar in loop ###" << endl;
      }
    }

    /* Table data */
    getMem4D(&Data, 0, nvar - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
        "Chemtable::Load Data", false);
    for (l = 0; l < nvar; l++)
    {
      for (k = 0; k < n3; k++)
      {
        for (j = 0; j < n2; j++)
        {
          for (i = 0; i < n1; i++)
          {
            dum = fread(&Data[l][k][j][i], sizeof(double), 1, inFp);
            if (dum != 1)
            {
              cerr << "### Error reading file: data (" << i << j << k
              << l << ") ### " << endl;
              throw(-1);
            }
          }
        }
      }
    }
    fclose(inFp);

    /* Output useful information on table */
    if (mpi_rank == 0)
    {
      if (is_iMask)
      {
        cout << "Chemistry table loaded with masks" << endl;
      }
      else
      {
        cout << "Chemistry table loaded without masks" << endl;
      }
      cout << "The combustion model is: " << CombustionModel << endl;
      cout << "Chemistry table has size " << n1 << " x " << n2 << " x " << n3
      << endl;
      cout << "x1 = Zmean:  " << "\t min=" << MinVal(x1, n1) << "  at i="
      << MinIndex(x1, n1) << "\t\t max=" << MaxVal(x1, n1) << "  at i="
      << MaxIndex(x1, n1) << endl;
      cout << "x2 = Zvar:  " << "\t min=" << MinVal(x2, n2) << "  at j="
      << MinIndex(x2, n2) << "\t\t max=" << MaxVal(x2, n2) << "  at j="
      << MaxIndex(x2, n2) << endl;
      strcpy(buffer_m, &VarNames[nvar - 1][0]);
      cout << "x3 = " << buffer_m << ":  " << "\t min=" << MinVal(x3, n3)
      << "  at k=" << MinIndex(x3, n3) << "\t\t max=" << MaxVal(x3, n3)
      << "  at k=" << MaxIndex(x3, n3) << endl;
      cout << "and contains " << nvar << " variables: " << endl;
      for (l = 0; l < nvar; l++)
      {
        strcpy(buffer_m, &VarNames[l][0]);
        GetArrayInfo(Data[l], n3, n2, n1, minval, maxval, mi3, mi2, mi1,
            ma3, ma2, ma1);
        cout.width(3);
        cout << l << ". ";
        cout.width(10);
        cout << buffer_m;
        cout << "\t min=" << minval << "\t at i=" << mi1 << "  j=" << mi2
        << "  k=" << mi3;
        cout << "\t max=" << maxval << "\t at i=" << ma1 << "  j=" << ma2
        << "  k=" << ma3 << endl;
      }
    }
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Clean table and deallocate memory */
  void Unload()
  {
    freeMem1D(x1, 0, n1 - 1);
    x1 = 0;
    freeMem1D(x2, 0, n2 - 1);
    x2 = 0;
    freeMem1D(x3, 0, n3 - 1);
    x3 = 0;
    if (is_iMask)
    { /* if iMask is used */
      freeMem3D(iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);
      is_iMask = 0;
    }
    delete[] VarNames;
    VarNames = 0;
    freeMem4D(Data, 0, nvar - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);
    Data = 0;

    return;
  }

  /***********************************************************************************************************/

  /***********************************************************************************************************/

  /*! \brief Compute value based on N inputs (e.g., ZMean, Zvar, C) and name of variable.
   *
   *  Propagate inputs through ANN and returns output
   *  \param[in] tag Name of variable to lookup.
   *  \param[in] xIn Pointer to array of inputs
   *  \return Interpolated value of variable \p tag.
   */
  double Lookup(string tag, double *xIn, int nIn)
  {
    return val;
  }

  /***********************************************************************************************************/
};

/***********************************************************************************************************/
/***********************************************************************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
