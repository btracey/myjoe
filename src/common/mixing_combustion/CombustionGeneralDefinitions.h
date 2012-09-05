/**
 * \brief General definitions of constants, classes, structures and functions used in combustion classes.
 *
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */

#ifndef COMBUSTIONGENERALDEFINITIONS_H
#define COMBUSTIONGENERALDEFINITIONS_H

#include <math.h>
#include "myMem.h"
#include <iostream>
#include <fstream>
using std::cout;
using std::cerr;
using std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////////////////


/** 
 * Overall constants
 */
#define kNameLen      32                 ///< String length identifiers for name of chemical species.
#define kStrMed       64                 ///< Size of strings in Chemtable class (for compatibility with input file).
#define kStrLong      256                ///< Size of strings for file names.

///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*! \brief Contains indices and weights to interpolate the 3 dimensional chemistry table.
 *
 *  The 3 indices correspond to coordinates just lower than input value: x(ii) <= x < x(ii+1).
 *  The 3 weights are used for interpolation (size depends if linear or cubic).
 *  The 3 directions corresponds to i,j,k indices of table (Zmean, Zvar, chi/progress variable).
 */
class InterpolationIndex
{
  
  public:
    
    int     TableDim;          ///< Dimension of table.
    int     StencilSize;       ///< Size of stencil (i.e., number of weights).
    int     *index;            ///< Index corresponding to lower bound: x(index[0]) <= x < x(index[0]+1).
    double  **weight;          ///< Weights for interpolation.
    
    // Constructor
    InterpolationIndex(void);
    InterpolationIndex(int n, int m);
  
    // Destructor
    ~InterpolationIndex();

    // Copy Method
    int copy(InterpolationIndex &newcopy);
    void print(void);
    bool compare(InterpolationIndex &other);
  
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * Useful functions
 */

void	GetArrayInfo(double ***A, int n1, int n2, int n3, double &minval, double &maxval, 
                     int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad = 0); 

void    GetArrayInfo_float(float ***A, int n1, int n2, int n3, double &minval, double &maxval,
                           int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad = 0);

void  GetArrayInfo(double ****A, int l, int n1, int n2, int n3, double &minval, double &maxval, 
                   int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad = 0); 

void  GetArrayInfo_float(float ****A, int l, int n1, int n2, int n3, double &minval, double &maxval,
                   int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad = 0); 

// Template declaration and implementation MUST be in header file (sniff ..)
template <class T>
void  GetArrayInfo4D(T *****A, int ivar, int n1, int n2, int n3, int n4,
                     T &minval, T &maxval,
                     int &mi1, int &mi2, int &mi3, int &mi4,
                     int &ma1, int &ma2, int &ma3, int &ma4)
{
  minval = A[0][0][0][0][ivar];
  mi1 = 0;
  mi2 = 0;
  mi3 = 0;
  mi4 = 0;

  maxval = A[0][0][0][0][ivar];
  ma1 = 0;
  ma2 = 0;
  ma3 = 0;
  ma4 = 0;

  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j)
      for (int k = 0; k < n3; ++k)
        for (int l =0; l < n4; ++l)
        {
          if ( A[i][j][k][l][ivar] < minval)
          {
            minval = A[i][j][k][l][ivar];
            mi1 = i;
            mi2 = j;
            mi3 = k;
            mi4 = l;
          }
          if (A[i][j][k][l][ivar] > maxval)
          {
            maxval = A[i][j][k][l][ivar];
            ma1 = i;
            ma2 = j;
            ma3 = k;
            ma4 = l;
          }
        }
  return;
}

template <class T>
void  GetArrayInfo3D(T ****A, int ivar, int n1, int n2, int n3,
                     T &minval, T &maxval,
                     int &mi1, int &mi2, int &mi3,
                     int &ma1, int &ma2, int &ma3)
{
  minval = A[0][0][0][ivar];
  mi1 = 0;
  mi2 = 0;
  mi3 = 0;

  maxval = A[0][0][0][ivar];
  ma1 = 0;
  ma2 = 0;
  ma3 = 0;

  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j)
      for (int k = 0; k < n3; ++k)
        {
          if ( A[i][j][k][ivar] < minval)
          {
            minval = A[i][j][k][ivar];
            mi1 = i;
            mi2 = j;
            mi3 = k;
          }
          if (A[i][j][k][ivar] > maxval)
          {
            maxval = A[i][j][k][ivar];
            ma1 = i;
            ma2 = j;
            ma3 = k;
          }
        }
  return;
}


void BinarySearch(int &i, double &w, double *Xarr, int i1, int i2, double x);
void BinarySearch(int &i, double *Xarr, int i1, int i2, double &x);


inline double MaxVal(double *A, int nmax, int pad)
{
        double maxval;
        maxval = A[pad];
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] > maxval) maxval = A[i];
        return maxval;
}

inline int MaxVal(int *A, int nmax, int pad)
{
        int maxval;
        maxval = A[pad];
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] > maxval) maxval = A[i];
        return maxval;
}

inline double  MinVal(double *A, int nmax, int pad)
{
        double minval;
        minval = A[pad];
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] < minval) minval = A[i];
        return minval;
}

inline int     MinVal(int *A, int nmax, int pad)
{
        int minval;
        minval = A[pad];
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] < minval) minval = A[i];
        return minval;
}

inline int     MaxIndex(double *A, int nmax, int pad)
{       
        int     maxindex = pad;
        for (int i = pad; i < nmax-pad; i++) 
                if (A[i] > A[maxindex]) maxindex = i;
        return maxindex;
}

inline int     MaxIndex(int *A, int nmax, int pad)
{       
        int     maxindex = pad;
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] > A[maxindex]) maxindex = i;
        return maxindex;
}

inline int     MinIndex(double *A, int nmax, int pad)
{
        int     minindex = pad;
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] < A[minindex]) minindex = i;
        return minindex;
}

inline int     MinIndex(int *A, int nmax, int pad)
{
        int     minindex = pad;
        for (int i = pad; i < nmax-pad; i++)
                if (A[i] < A[minindex]) minindex = i;
        return minindex;
}


inline double linearFunction(double x, double y, double z)
{
   return    2.0  * x*y*z
           + 0.25 * x*y
           - 0.7  * x*z
           + 0.2  * y*z
           - 0.3  * x
           + 1.0  * y
           - 0.5  * z
           + 1.25;
}

inline double parabolicFunction(double x, double y, double z)
{
  return (x-0.5)*(x-0.5)*(y-0.25)*y*(z-0.25);

//   return  - 2.0  * x*x*y*y*z*z 
//           + 0.25 * x*x*y*y*z
//           - 0.7  * x*y*y*z*z                                               
//           + 0.2  * x*y*z*z
//           - 0.3  * y*y*z
//           + 1.0  * y*z*z
//           - 0.5  * x
//           + 0.3  * y*z
//           - 1.0  * z
//           + 0.75 * x*z*z
//           + 2.5  * y*y
//           + 1.25;
}

inline double cubicFunction(double x, double y, double z)
{
  return 20.0*x*x*(x-0.95)*(y-0.25)*(y-0.25)*(y-0.25)*z*(z-0.3);

//   return    2.0  * x*x*x*y*y*y*z*z*z 
//           + 0.25 * x*x*x*y*y*z
//           - 0.7  * x*x*y*y*z*z                                               
//           + 0.2  * x*y*y*z*z
//           - 0.3  * y*y*z
//           + 1.0  * y*z*z*z
//           - 0.5  * x
//           + 0.3  * y*z
//           - 1.0  * z
//           + 0.75 * x*z*z
//           + 2.5  * (y-0.1)*(y-0.1)
//           + 1.25;
}

inline double otherAnalyticalFunction(double x, double y, double z)
{
   return exp(x*y-z) - z*z + sqrt(x) - y*y + 10.0;
}


#endif

