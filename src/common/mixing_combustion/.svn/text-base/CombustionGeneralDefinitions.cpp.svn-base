#include "CombustionGeneralDefinitions.h"
/**
 * \brief General definitions of constants, classes, structures and functions used in combustion classes.
 *
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */
InterpolationIndex::InterpolationIndex(void) {
  TableDim = 0;
  StencilSize = 0;
  index = NULL;
  weight = NULL;
}

InterpolationIndex::InterpolationIndex(int n, int m) 
    {
      TableDim = n;
      StencilSize = m;
      index = new int[n];
      for (int i=0; i<n; ++i) index[i] = 0;
      getMem2D(&weight, 0, TableDim-1, 0, StencilSize-1, "InterpolationIndex: weight", true);
    }
  
InterpolationIndex::~InterpolationIndex()
    {
      if (index  != NULL) {delete [] index; index = NULL;}
      if (weight != NULL) {freeMem2D(weight, 0, TableDim-1, 0, StencilSize-1); weight = NULL;}
    }
  
int InterpolationIndex::copy(InterpolationIndex &newcopy){
  if (newcopy.TableDim != TableDim) return -1;
  if (newcopy.StencilSize != StencilSize) return -1;
  for (int i=0; i<TableDim; ++i) {
    newcopy.index[i] = index[i];
    for (int j=0; j<StencilSize; ++j) {
      newcopy.weight[i][j] = weight[i][j];
    }
  }
  return 0;
}

void InterpolationIndex::print(void) {
  cout << "#--------------------------#"<<endl;
  cout << "  - TableDim:    " << TableDim  <<endl;
  cout << "  - StencilSize: " << StencilSize<< endl;
  cout << "  - Index:       ";
  for (int i=0; i<TableDim; ++i) cout << index[i] << " ";
  cout << endl;
  cout << "  - Weights:     ";
  for (int i=0; i<TableDim; ++i) cout << weight[i][0] << " ";
  cout << endl;
  cout << "#--------------------------#"<<endl;
}

bool InterpolationIndex::compare(InterpolationIndex &other) {
  bool test = true;
  test = test&&(other.TableDim == TableDim);
  test = test&&(other.StencilSize == StencilSize);
  if (!test) return test;
  for (int i=0; i<TableDim; ++i) {
    test = test&&(other.index[i] == index[i]);
    for (int j=0; j<StencilSize; ++j) {
          test = test&&(other.weight[i][j] == weight[i][j]);
        }
  }
  return test;
}
/**
 * Useful functions
 */


void  GetArrayInfo(double ***A, int n1, int n2, int n3, double &minval, double &maxval,
                          int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad)
{
        minval = A[pad][pad][pad];
        mi1 = pad; mi2 = pad; mi3 = pad;
        maxval = A[pad][pad][pad];
        ma1 = pad; ma2 = pad; ma3 = pad;
        for (int i = pad; i < n1-pad; i++)
                for (int j = pad; j < n2-pad; j++)
                        for (int k = pad; k < n3-pad; k++)
                        {
                                if (A[i][j][k] < minval)
                                {
                                        minval = A[i][j][k];
                                        mi1 = i; mi2 = j; mi3 = k;
                                }
                                if (A[i][j][k] > maxval)
                                {
                                        maxval = A[i][j][k];
                                        ma1 = i; ma2 = j; ma3 = k;
                                }
                        }
        return;
}

void  GetArrayInfo_float(float ***A, int n1, int n2, int n3, double &minval, double &maxval, 
                          int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad)
{
	minval = (double) A[pad][pad][pad];
	mi1 = pad; mi2 = pad; mi3 = pad;
	maxval = (double) A[pad][pad][pad];
	ma1 = pad; ma2 = pad; ma3 = pad;
	for (int i = pad; i < n1-pad; i++) 
		for (int j = pad; j < n2-pad; j++) 
			for (int k = pad; k < n3-pad; k++) 
			{
				if ((double) A[i][j][k] < minval) 
				{
					minval = (double) A[i][j][k];
					mi1 = i; mi2 = j; mi3 = k;
				}
				if ((double) A[i][j][k] > maxval) 
				{
					maxval = (double) A[i][j][k];
					ma1 = i; ma2 = j; ma3 = k;
				}
			}
	return;
}

void  GetArrayInfo(double ****A, int l, int n1, int n2, int n3, double &minval, double &maxval, 
                          int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad)
{
  minval = A[pad][pad][pad][l];
  mi1 = pad; mi2 = pad; mi3 = pad;
  maxval = A[pad][pad][pad][l];
  ma1 = pad; ma2 = pad; ma3 = pad;
  for (int i = pad; i < n1-pad; i++) 
    for (int j = pad; j < n2-pad; j++) 
      for (int k = pad; k < n3-pad; k++) 
      {
        if (A[i][j][k][l] < minval) 
        {
          minval = A[i][j][k][l];
          mi1 = i; mi2 = j; mi3 = k;
        }
        if (A[i][j][k][l] > maxval) 
        {
          maxval = A[i][j][k][l];
          ma1 = i; ma2 = j; ma3 = k;
        }
      }
  return;
}

void  GetArrayInfo_float(float ****A, int l, int n1, int n2, int n3, double &minval, double &maxval,
                          int &mi1, int &mi2, int &mi3, int &ma1, int &ma2, int &ma3, int pad)
{
  minval = (double) A[pad][pad][pad][l];
  mi1 = pad; mi2 = pad; mi3 = pad;
  maxval = (double) A[pad][pad][pad][l];
  ma1 = pad; ma2 = pad; ma3 = pad;
  for (int i = pad; i < n1-pad; i++)
    for (int j = pad; j < n2-pad; j++)
      for (int k = pad; k < n3-pad; k++)
      {
        if ((double) A[i][j][k][l] < minval)
        {
          minval = (double) A[i][j][k][l];
          mi1 = i; mi2 = j; mi3 = k;
        }
        if ((double) A[i][j][k][l] > maxval)
        {
          maxval = (double) A[i][j][k][l];
          ma1 = i; ma2 = j; ma3 = k;
        }
      }
  return;
}

void BinarySearch(int &i, double &w, double *Xarr, int i1, int i2, double x)
// return the index i for which Xarr[i] <= x < Xarr[i+1] and corresponding weight
{
  // check if out of array bounds
  if (x <= Xarr[i1])
  {
    i = i1;
    w = 1.0;
    return;
  }
  
  if (x >= Xarr[i2])
  {
    i = i2 - 1;
    w = 0.0;
    return;
  }
  
  // find index and weight
  while (i1 < i2-1)
  {
    int mid = (i1 + i2) / 2;
    if (x > Xarr[mid])
      i1 = mid;
    else if (x < Xarr[mid])
      i2 = mid;
    else
    {
      i = mid;
      w = 1.0;
      return;
    }
  }
  i = i1;
  w = (Xarr[i+1] - x) / (Xarr[i+1] - Xarr[i]);
}

void BinarySearch(int &i, double *Xarr, int i1, int i2, double &x)
// return the index i for which Xarr[i] <= x < Xarr[i+1]
// NOTE: Upon return, x is clipped to min/max value!!!
{
  // check if out of array bounds
  if (x <= Xarr[i1])
  {
    i = i1;
    x = Xarr[i1];
    return;
  }
  
  if (x >= Xarr[i2])
  {
    i = i2 - 1;
    x = Xarr[i2];
    return;
  }
  
  // find index and weight
  while (i1 < i2-1)
  {
    int mid = (i1 + i2) / 2;
    if (x > Xarr[mid])
      i1 = mid;
    else if (x < Xarr[mid])
      i2 = mid;
    else
    {
      i = mid;
      return;
    }
  }
  i = i1;
}
