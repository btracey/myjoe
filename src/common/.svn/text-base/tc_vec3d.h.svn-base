#ifndef TC_VEC3D_H
#define TC_VEC3D_H

#include <math.h> // cmath on g++ ?

void set_vec3d(int * v1,const int * v2);
void set_vec3d(double * v1,const double * v2);
void add_vec3d(double * v1,const double * v2);
void subtract_vec3d(double * v1,const double * v2);
void normalize_vec3d(double * v1);
double dotproduct_vec3d(const double * v1,const double * v2);
void divide_vec3d(double * v1,const double scalar);

inline double vecDotVec3d(const double * v1, const double * v2) {

  double dotproduct = 0.0;

  for (int i = 0; i < 3; i++)
    dotproduct += v1[i]*v2[i];

  return(dotproduct);
}

inline double normVec3d(double * v1) {

  double mag = vecDotVec3d(v1, v1);

  mag = sqrt(mag);

  for (int i = 0; i < 3; i++)
    v1[i] /= mag;

  return mag;
}

inline double normVec3d(double *vNorm, const double *v1) {

  double mag = vecDotVec3d(v1, v1);

  mag = sqrt(mag);

  for (int i = 0; i < 3; i++)
    vNorm[i] = v1[i]/mag;

  return mag;
}

inline void vecMinVec3d(double *vRes, const double *v1, const double *v2) {

  for (int i = 0; i < 3; i++)
    vRes[i] = v1[i] - v2[i];
}

inline void matMultMatR5(double (*res)[5], double (*m1)[5], double (*m2)[5])  // m*m matrices
{
  for (int i=0; i<5; i++)
  for (int j=0; j<5; j++)
  {
    res[i][j] = 0.0;
    for (int k=0; k<5; k++)
      res[i][j] += m1[i][k]*m2[k][j];
  }
}

inline void matMultDiagMatR5(double (*res)[5], double (*mat)[5], double *diag)  // m*m matrices
{
  for (int i=0; i<5; i++)
  for (int j=0; j<5; j++)
    res[i][j] = mat[i][j]*diag[j];
}

inline void matMultVecR5(double *res, double (*mat)[5], double *vec)  // m*m matrices !!!!!!
{
  for (int i=0; i<5; i++)
  {
    res[i] = 0.0;

    for (int j=0; j<5; j++)
      res[i] += mat[i][j]*vec[j];
  }
}







#endif




