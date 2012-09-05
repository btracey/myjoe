#include "tc_vec3d.h"

void set_vec3d(int * v1,const int * v2) {
  int i;
  for (i = 0; i < 3; i++) 
    v1[i] = v2[i];
}

void set_vec3d(double * v1,const double * v2) {
  int i;
  for (i = 0; i < 3; i++) 
    v1[i] = v2[i];
}

void add_vec3d(double * v1,const double * v2) {
  int i;
  for (i = 0; i < 3; i++) 
    v1[i] += v2[i];
}

void subtract_vec3d(double * v1,const double * v2) {
  int i;
  for (i = 0; i < 3; i++) 
    v1[i] -= v2[i];
}

void normalize_vec3d(double * v1) {
  double mag = 0.0;
  int i;
  for (i = 0; i < 3; i++) 
    mag += v1[i]*v1[i];
  mag = sqrt(mag);
  for (i = 0; i < 3; i++) 
    v1[i] /= mag;
}

void divide_vec3d(double * v1,const double scalar) {
  int i;
  for (i = 0; i < 3; i++) 
    v1[i] /= scalar;
}

        
