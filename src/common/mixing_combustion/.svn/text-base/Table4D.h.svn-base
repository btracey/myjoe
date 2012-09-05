#ifndef TABLE4D_H
#define TABLE4D_H

#include "CombustionGeneralDefinitions.h"
#include <algorithm>
using std::min;
using std::max;

#include "MpiStuff.h"
using namespace MpiStuff;

/*! \brief Class to describe 4D databases
 * 
 *  Class assumes a structured Cartesian table: x1=x1(i), x2=x2(j), x3=x3(k), x4=x4(k)
 *
 */
class Table4D 
{
public:  

  //------------------------------//
  //         Constructors         //
  //------------------------------//
  Table4D(void);
  Table4D(int i1, int i2, int i3, int i4, int i5);
 
  //------------------------------//
  //         Destructor           //
  //------------------------------// 
  ~Table4D();
  
  //------------------------------//
  //     Get/Set Attributes       //
  //------------------------------//

  // Get
  inline int GetDimension1(void) {return n1;}
  inline int GetDimension2(void) {return n2;}
  inline int GetDimension3(void) {return n3;}
  inline int GetDimension4(void) {return n4;}
  inline int GetNumberVariables(void) {return nvar;}

  inline double GetCoordinate1(int i1) {return x1[i1];}
  inline double GetCoordinate2(int i2) {return x2[i2];}
  inline double GetCoordinate3(int i3) {return x3[i3];}
  inline double GetCoordinate4(int i4) {return x4[i4];}

  void GetInterpPoint(double &coor1, double &coor2, double &coor3, double &coor4) {
    coor1 = point[0];
    coor2 = point[1];
    coor3 = point[2];
    coor4 = point[3];
  }
  
  void GetInterpolation(InterpolationIndex &copy) { interp->copy(copy); }

  float GetTableValue(int i1, int i2, int i3, int i4, int ivar) {return Data[i1][i2][i3][i4][ivar];}
  double GetTableValueDouble(int i1, int i2, int i3, int i4, int ivar) {return (double) Data[i1][i2][i3][i4][ivar];}

  // Set
  void SetDimensions(int i1, int i2, int i3, int i4, int i5) {
    n1 = i1;
    n2 = i2;
    n3 = i3;
    n4 = i4;
    nvar = i5;
  }
  void SetCoordinate1(int i, double x) {x1[i] = x;}
  void SetCoordinate2(int i, double x) {x2[i] = x;}
  void SetCoordinate3(int i, double x) {x3[i] = x;}
  void SetCoordinate4(int i, double x) {x4[i] = x;}

  void SetInterpPoint(double coor1, double coor2, double coor3, double coor4) {
    // Clip coordinates
    coor1 = min(max(coor1, x1[0]), x1[n1-1]);
    coor2 = min(max(coor2, x2[0]), x2[n2-1]);
    coor3 = min(max(coor3, x3[0]), x3[n3-1]);
    coor4 = min(max(coor4, x4[0]), x4[n4-1]);
    // Set Point
    point[0] = coor1;
    point[1] = coor2;
    point[2] = coor3;
    point[3] = coor4;
  }

  void SetInterpolation(InterpolationIndex &my_interp) {
  	int test = my_interp.copy(*interp);
  	if (test != 0) {
  		cout << "could not set Interpolation Index"<<endl;
  		throw(-1);
  	}
  }

  void SetTableValue(int i1, int i2, int i3, int i4, int ivar, float value) {
    Data[i1][i2][i3][i4][ivar] = value;
  }
  void SetTableValue(int i1, int i2, int i3, int i4, int ivar, double value) {
    Data[i1][i2][i3][i4][ivar] = (float) value;
  }

  void SetAllocated(); 

  //------------------------------// 
  //    De-/Allocation Methods    // 
  //------------------------------// 
  void Allocate(int i1, int i2, int i3, int i4, int i5);
  void Allocate();
  void DestroyTable();

  //------------------------------//
  //       Other Methods          //
  //------------------------------//

  //-----------------------------
  // Check if the database and its
  // vector coordinates have already
  // been allocated
  //-----------------------------
  bool IsAllocated();
  //-----------------------------
  // Return min and max value of the ivar-th
  // variable and the location of the min and
  // max value in the database
  //-----------------------------
  void GetArrayInfo(int ivar, float &minval, float &maxval,
                    int &mi1, int &mi2, int &mi3, int &mi4,
                    int &ma1, int &ma2, int &ma3, int &ma4);

  //-----------------------------
  // Return a copy of x1
  //   the copy array is assumed already allocated
  //   at the right size
  //-----------------------------
  void CopyCoordinate1(double *x);
  //-----------------------------
  // Return a copy of x2
  //   the copy array is assumed already allocated
  //   at the right size
  //-----------------------------
  void CopyCoordinate2(double *x);
  //-----------------------------
  // Return a copy of x3
  //   the copy array is assumed already allocated
  //   at the right size
  //-----------------------------
  void CopyCoordinate3(double *x);
  //-----------------------------
  // Return a copy of x4
  //   the copy array is assumed already allocated
  //   at the right size
  //-----------------------------
  void CopyCoordinate4(double *x);
  //--------------------------------------
  // Set interpolation point in the database
  // and compute the corresponding interpolation
  // index and weights in the class
  //--------------------------------------
  void InterpolatePoint(double coor1, double coor2, double coor3, double coor4);
  void InterpolatePoint_old(double coor1, double coor2, double coor3, double coor4);
  //------------------------------------
  // Interpolate ivar-th variable in the database
  //   Example: var = Interpolate(val1, val2, val3, val4, i);
  //   val1, val2, val3, val4 are the databse coordinates
  //------------------------------------
  double Interpolate(double coor1, double coor2, double coor3, double coor4, int ivar);
  //-----------------------------------------
  // Return the ivar-th interpolated value of the database
  //   interpolation point and its interpolation index and weights
  //   are assumed already computed
  //-----------------------------------------
  double Interpolate(int ivar);
  //-----------------------------------------
  // Master process broadcasrs database to other ones 
  //-----------------------------------------
  void BCAST(int mpi_root, MPI_Comm mpi_comm);

private:
  int                n1, n2, n3, n4;
  int                nvar; 
  double             *x1, *x2, *x3, *x4;
  float              *****Data;
  double             point[4];
  int                interp_mode[4];
  InterpolationIndex *interp;
  bool               Allocated;

  void Init(void);
  //-----------------------------------------
  // Copy a coordinate vector to an other one
  //   the pointer *copy is assumed to have
  //   been already allocated with the right size
  //-----------------------------------------
  void CopyCoordinate(int n, double* x, double* copy);
  //-----------------------------------------
  // Compute interpolation index and weights
  //   interpolation point is assumed to be already set
  //-----------------------------------------
  void InterpolatePoint_old(void);
  void InterpolatePoint(void);
  //-----------------------------------------
  // Compute interpolation and weight in a 1D vector
  //  with different mode (0: BinarySearch
  //                       1: Linear
  //                       2: power law
  //                       3: two linear zone
  //                       4: Linear then linear stretch)
  //-----------------------------------------
  void ComputeIndexAndWeight(int &index, double &weight, double *vec, int size, double coor, int mode);
  //-----------------------------------------
  // Compute weight and index from linear mesh
  //   0 < coor < 1
  //-----------------------------------------
  void ComputeIndexAndWeightLin(int &index, double &weight, double *vec, int size, double coor);
  //-----------------------------------------
  // Compute weight and index from power-law mesh
  //   0 < coor < 1
  //-----------------------------------------
  void ComputeIndexAndWeightPow(int &index, double &weight, double *vec, int size, double coor, double power);
  //-----------------------------------------
  // Compute weight and index from two-linear zones mesh
  //   NOT VALIDATED
  //   z_min < coor < z_max
  //-----------------------------------------
  void ComputeIndexAndWeightTwoLin(int &index, double &weight, double *vec, int size, double coor,
                                   double z_min, double z_max, double z_cut, int i_cut);
  //-----------------------------------------
  // Compute weight and index from linear mesh followed by linear strecthing
  //   0 < coor < 1
  //-----------------------------------------
  void ComputeIndexAndWeightLinAndStretched(int &index, double &weight, double *vec, int size, double coor, double z_cut, int i_cut);
};
#endif

