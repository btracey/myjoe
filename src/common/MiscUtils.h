#ifndef MISCUTILS_H
#define MISCUTILS_H

#include <iostream>
#include <cstring>    /* Added header to compile with GNU compiler >= 4.3 */
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>

#ifdef NO_ASSERT
#include <assert.h>
#endif

#include "MpiStuff.h"
using namespace MpiStuff;
using namespace std;

/**
 * MiscUtils 
 * \author Frank Ham and Gianluca Iaccarino
 */
namespace MiscUtils {

  /*
   * calcTaylorVelocityAndPressure
   */

  void calcTaylorVelocityAndPressure(double * u, double * p_ptr, 
				     const double * x, const double time, const double nu);

  /*
   * inviscid (Euler) vortex
   */
  
  void calcInviscidVortex(double * rho,double (*rhou)[3],double * rhoE,
			  const double (*x)[3],const int n,
			  const double t,const double gamma,
			  const double rho_inf,const double p_inf,const double u_inf,
			  const int periodic_flag);

  int atox(char * token);

  double atod(char * token);

  void nullTerminate(char *name);

  void nullTerminate(char *name, int len);

  void calcUniformDist(int * xod, const int nx, const int ndist);

  /**
   *  dumping a scalar range...
   */
  void dumpScalarRange(const double * s, const int n, const string& message);

  void dumpScalarRange(const double * s, const int n, const char * message);

  void dumpScalarRange(const int * s, const int n, char * message);

  void dumpScalarBins(const int * s, const int n, char * message);

  /**
   *  dumping a vector range...
   */

  void dumpVectorRange(const double (*v)[3], const int n, char * message);

  void dumpTensorRange(const double (*t)[3][3], const int n, char * message);

  /**
   *  byte swapping...
   */
  void byteSwap(int * ptr, const int n);

  void byteSwap(int (*ptr)[2], const int n1, const int n2);

  int byteSwap(const int old_value);

  MPI_Offset byteSwap(const MPI_Offset old_value);

  void byteSwap(double * ptr, const int n);

  void byteSwap(double (*ptr)[3], const int n1, const int n2);

  void byteSwap(double (*ptr)[3][3], const int n1, const int n2, const int n3);

  double byteSwap(const double old_value);

  void reorder_csr(int * x_i, int * x_v, int * order, const int n);

  void reorder_csr(int * x_i, int * x_v1, int * x_v2, int * order, const int n);

  void reorder(double * var,const int * order,const int n);

  void reorder(double (*var)[3],const int * order,const int n);

  void reorder(int * v,const int * order,const int n);
  
  void reorder(int (*v2)[2],const int * order,const int n);
  
  /**
   * checks if a token is in string
   * \param string 
   * \return 1 ... if token on lhs
   * \return 0 ... if no token
   * \return -1 ... if no token is within non-whitespace string
   */
  int checkHash(const string &str);

  /**
   * takes a string and seperates it into tokens
   * \return tokens stored in vector !!!
   * \param str ... processed string
   * \param str ... delimiters, standard delimiter is whitespace
   */
  void tokenizeString(vector<string> &tokens, const string &str,
                      const string &delimiters = " ");

  string getFileNameExtension(const string &str, const string &delimiters);

  // does a file exist?
  
  bool fileExists(const string& filename);

  bool fileExists(const char * filename);

  // resizing utilities.
  void resize(int * &scalar,const int n_new);  
  void resize(int * &scalar,const int n_old,const int n_new);  

  void resize(int (* &scalar)[2],const int n_new);
  void resize(int (* &scalar)[2],const int n_old,const int n_new);

  void resize(int (* &scalar)[4],const int n_new);
  void resize(int (* &scalar)[8],const int n_new);

  void resize(int (* &scalar)[4],const int n_old,const int n_new);
  void resize(int (* &scalar)[4][2],const int n_old,const int n_new);

  void resize(int (* &scalar)[6],const int n_old,const int n_new);
  void resize(int (* &scalar)[6][4],const int n_old,const int n_new);

  void resize(double * &scalar,const int n_new);
  void resize(double * &scalar,const int n_old,const int n_new);

  void resize(double (* &d3)[3],const int n_new);
  void resize(double (* &d3)[3],const int n_old,const int n_new);

}

#endif

