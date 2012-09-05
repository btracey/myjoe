#include "MiscUtils.h"
#include <sys/stat.h> // for fileExists
#include <stdio.h>
#include <string.h>

namespace MiscUtils {
  
  void calcTaylorVelocityAndPressure(double * u, double * p_ptr,
				     const double * x,const double time, const double nu) {
    
    double taylor_dx = 0.5;
    double taylor_dy = 0.5;
    double taylor_u = 0.25;
    double taylor_v = 0.15;
    //double taylor_u = 0.0;
    //double taylor_v = 0.0;
    //double taylor_sign = 1.0;

    double taylor_x = x[0] - taylor_dx - taylor_u*time;
    double taylor_y = x[1] - taylor_dy - taylor_v*time;
    double tmp = exp(-2.0*M_PI*M_PI*time*nu);
    u[0] = -cos(M_PI*taylor_x)*sin(M_PI*taylor_y)*tmp + taylor_u;
    u[1] = sin(M_PI*taylor_x)*cos(M_PI*taylor_y)*tmp + taylor_v;
    u[2] = 0.0;
    tmp = exp(-4.0*M_PI*M_PI*time*nu);
    *p_ptr = -0.25*(cos(2.0*M_PI*taylor_x) + cos(2.0*M_PI*taylor_y))*tmp;
  }
  
  void calcInviscidVortex(double * rho,double (*rhou)[3],double * rhoE,
			  const double (*x)[3],const int n,
			  const double t,const double gamma,
			  const double rho_inf,const double p_inf,const double u_inf,
			  const int periodic_flag) {
  
    // x,y position of the vortex
    // centered
    const double x0 = 0.0;
    const double y0 = 0.0;
    // near top right - for checking bc quickly...
    //double x0 = 3.0;
    //double y0 = 2.0;

    // direction of propagation
    const double theta = M_PI/3.0;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const double Ma_inf = u_inf/sqrt(gamma*p_inf/rho_inf);
    const double rc = 1.0;
    
    // circulation parameter...
    //const double e_twopi = 0.001; // very weak
    //double e_twopi = 0.005;
    const double e_twopi = 0.08; // normal
    //double e_twopi = 0.1;
    //double e_twopi = 0.4; //very strong 

    // setup...
    const double coeff = 0.5 * e_twopi*e_twopi * (gamma-1.0) * Ma_inf*Ma_inf;
    for (int i = 0; i < n; i++) {
      
      double dx = x[i][0] - x0 - u_inf*cos_theta*t;
      double dy = x[i][1] - y0 - u_inf*sin_theta*t;
    
      // for periodic -5 < x < 5, -5 < y < 5, bring the vortex back
      // into the domain...
      if (periodic_flag == 1) {
        while(dx > 5.0) dx -= 10.0;
        while(dx < -5.0) dx += 10.0;
        while(dy > 5.0) dy -= 10.0;
        while(dy < -5.0) dy += 10.0;
      }

      const double f0 = 1.0 - (( dx*dx ) + ( dy*dy ))/( rc*rc );
      rho[i] = rho_inf*pow( 1.0 - coeff * exp( f0 ) , 1.0/(gamma-1.0) );
      rhou[i][0] = rho[i]*u_inf*( cos_theta - e_twopi * ( dy )/rc * exp( f0 / 2.0 ) );
      rhou[i][1] = rho[i]*u_inf*( sin_theta + e_twopi * ( dx )/rc * exp( f0 / 2.0 ) );
      rhou[i][2] = 0.0;
    
      const double p = p_inf*pow( 1.0 - coeff * exp( f0 ) , gamma/(gamma-1.0) );
      rhoE[i] = p/(gamma-1.0) + 0.5*(rhou[i][0]*rhou[i][0] + rhou[i][1]*rhou[i][1] + rhou[i][2]*rhou[i][2])/rho[i];

    }

  }

  int atox(char * token) {
    int i;
    sscanf(token, "%x", &i);
    return (i);
  }

  double atod(char * token) {
    double d;
    sscanf(token, "%lf", &d);
    return (d);
  }

  void nullTerminate(char *name) {
    // f90 strings stored in binary files are space-terminated. This
    // adds a NULL termination after the last non-space character.
    int i = strlen(name)-1;

    while ((i > 0) && (name[i-1] == ' ')) {
      i--;
    }
    name[i] = '\0';
  }

  void nullTerminate(char *name, int len) {
    // f90 strings stored in binary files are space-terminated. This
    // adds a NULL termination after the last non-space character.
    int i = len;

    while ((i > 0) && (name[i-1] == ' ')) {
      i--;
    }
    name[i] = '\0';
  }

  void calcUniformDist(int * xod, const int nx, const int ndist) {
    xod[0] = 0;
    for (int id = 1; id <= ndist; id++) {
      xod[id] = (int)((double)id/(double)ndist*(double)nx + 0.5);
    }
  }
  
  void dumpScalarRange(const double * s, const int n, const string& message) {
    dumpScalarRange(s,n,message.c_str());
  }

  void dumpScalarRange(const double * s, const int n,const char * message) {

    double my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = 1.0E+20; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
      // nan check...
      assert( s[i] == s[i] ); 
    }

    double buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << " > dumpScalarRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

  }

  void dumpScalarRange(const int * s, const int n, char * message) {

    int my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = 1000000000; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
    }

    int buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_INT, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << " > dumpScalarRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

  }

  void dumpScalarBins(const int * s, const int n, char * message) {

    int my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = 1000000000; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
    }

    int buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_INT, MPI_MIN, mpi_comm);

    int nbins = -buf[1] - buf[0]+ 1;
    assert( (nbins > 0)&&(nbins < 100));

    int * my_bin_count = new int[nbins];
    for (int i = 0; i < nbins; i++)
      my_bin_count[i] = 0;

    for (int i = 0; i < n; i++)
      my_bin_count[s[i]-buf[0]] += 1;

    int * bin_count = new int[nbins];
    MPI_Reduce(my_bin_count, bin_count, nbins, MPI_INT, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > dumpScalarBins: "<< message << endl;
      for (int i = 0; i < nbins; i++)
        cout << " > value: "<< i+buf[0]<< ", count: "<< bin_count[i]<< endl;
    }

    delete[] my_bin_count;
    delete[] bin_count;

  }

  void dumpVectorRange(const double (*v)[3], const int n, char * message) {

    double my_buf[6];
    for (int j = 0; j < 6; j++)
      my_buf[j] = 1.0E+20; // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        my_buf[2*j] = min(my_buf[2*j], v[i][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1], -v[i][j]);
	// nan check...
	assert( v[i][j] == v[i][j] ); 
      }
    }

    double buf[6];
    MPI_Reduce(my_buf, buf, 6, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > dumpVectorRange: "<< message;
      for (int j = 0; j < 3; j++) {
        cout << ", "<< j << ": "<< buf[2*j]<< ":"<< -buf[2*j+1];
      }
      cout << endl;
    }

  }

  void dumpTensorRange(const double (*t)[3][3], const int n, char * message) {

    double my_buf[18];
    for (int j = 0; j < 18; j++)
      my_buf[j] = 1.0E+20; // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
	  my_buf[2*(3*j+k)] = min(my_buf[2*(3*j+k)], t[i][j][k]);
	  my_buf[2*(3*j+k)+1] = min(my_buf[2*(3*j+k)+1], -t[i][j][k]);
	  // nan check...
	  assert( t[i][j][k] == t[i][j][k] );
	}
      }
    }
    
    double buf[18];
    MPI_Reduce(my_buf, buf, 18, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0) {
      for (int j = 0; j < 3; j++) {
	cout << " > dumpTensorRange: "<< message;
	for (int k = 0; k < 3; k++) {
	  cout << ", " << j << "," << k << ": "<< buf[2*(3*j+k)]<< ":"<< -buf[2*(3*j+k)+1];
	}
	cout << endl;
      }
    }
    
  }
  
  void byteSwap(int * ptr, const int n) {
    for (int i = 0; i < n; i++) {
      ptr[i] = byteSwap(ptr[i]);
    }
  }

  void byteSwap(int (*ptr)[2], const int n1, const int n2) {
    assert(n2 == 2);
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < 2; j++) {
        ptr[i][j] = byteSwap(ptr[i][j]);
      }
    }
  }

  int byteSwap(const int old_value) {
    int new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[3];
    cnew[1] = cold[2];
    cnew[2] = cold[1];
    cnew[3] = cold[0];
    return (new_value);
  }

  MPI_Offset byteSwap(const MPI_Offset old_value) {
    MPI_Offset new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[7];
    cnew[1] = cold[6];
    cnew[2] = cold[5];
    cnew[3] = cold[4];
    cnew[4] = cold[3];
    cnew[5] = cold[2];
    cnew[6] = cold[1];
    cnew[7] = cold[0];
    return (new_value);
  }

  void byteSwap(double * ptr, const int n) {
    for (int i = 0; i < n; i++) {
      ptr[i] = byteSwap(ptr[i]);
    }
  }

  void byteSwap(double (*ptr)[3], const int n1, const int n2) {
    assert(n2 == 3);
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < 3; j++) {
        ptr[i][j] = byteSwap(ptr[i][j]);
      }
    }
  }

  void byteSwap(double (*ptr)[3][3], const int n1, const int n2, const int n3) {
    assert(n2 == 3);
    assert(n3 == 3);
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
	  ptr[i][j][k] = byteSwap(ptr[i][j][k]);
	}
      }
    }
  }

  double byteSwap(const double old_value) {
    double new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[7];
    cnew[1] = cold[6];
    cnew[2] = cold[5];
    cnew[3] = cold[4];
    cnew[4] = cold[3];
    cnew[5] = cold[2];
    cnew[6] = cold[1];
    cnew[7] = cold[0];
    return (new_value);
  }

  void reorder_csr(int * x_i, int * x_v, int * order, const int n) {

    // backup the old data...
    int * x_i_old = new int[n+1];
    for (int i = 0; i <= n; i++)
      x_i_old[i] = x_i[i];

    int * x_v_old = new int[x_i[n]];
    for (int j = 0; j < x_i[n]; j++)
      x_v_old[j] = x_v[j];

    // put the counts into x_i...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      x_i[i_new+1] = x_i_old[i+1] - x_i_old[i]; // put the count
    }

    // build new CSR...
    x_i[0] = 0;
    for (int i = 0; i < n; i++)
      x_i[i+1] += x_i[i];

    // copy in the values...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      int j_new_f = x_i[i_new];
      int j_f = x_i_old[i];
      int j_l = x_i_old[i+1]-1;
      for (int j = j_f; j <= j_l; j++) {
        x_v[j-j_f+j_new_f] = x_v_old[j];
      }
    }

    delete[] x_i_old;
    delete[] x_v_old;
  }

  void reorder_csr(int * x_i, int * x_v1, int * x_v2, int * order, const int n) {
    
    // backup the old data...
    int * x_i_old = new int[n+1];
    for (int i = 0; i <= n; i++)
      x_i_old[i] = x_i[i];
    
    int * x_v1_old = new int[x_i[n]];
    int * x_v2_old = new int[x_i[n]];
    for (int j = 0; j < x_i[n]; j++) {
      x_v1_old[j] = x_v1[j];
      x_v2_old[j] = x_v2[j];
    }

    // put the counts into x_i...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      x_i[i_new+1] = x_i_old[i+1] - x_i_old[i]; // put the count
    }

    // build new CSR...
    x_i[0] = 0;
    for (int i = 0; i < n; i++)
      x_i[i+1] += x_i[i];

    // copy in the values...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      int j_new_f = x_i[i_new];
      int j_f = x_i_old[i];
      int j_l = x_i_old[i+1]-1;
      for (int j = j_f; j <= j_l; j++) {
        x_v1[j-j_f+j_new_f] = x_v1_old[j];
        x_v2[j-j_f+j_new_f] = x_v2_old[j];
      }
    }

    delete[] x_i_old;
    delete[] x_v1_old;
    delete[] x_v2_old;

  }

  int checkHash(const string &str) 
  {
    size_t i = 0;
    if (str.length() == 0)     return 1;
  
    // ignore whitespace and tab
    while (str[i] == ' ' || str[i] == '\t')  
      i++;
  
    if (str.at(i) == '#')      return 1;
    {
      int hash_position = (int)str.find('#');

      if (hash_position == -1)    return str.length();    //return length of string to be used for the parameter 
      else                        return hash_position;   //only string until hash in line appears will be used 
    }
  }

  void tokenizeString(vector<string> &tokens, const string &str, const string &delimiters) {
    string::size_type lastPos = str.find_first_not_of(delimiters, 0); // skip delimiters at beginning
    string::size_type pos = str.find_first_of(delimiters, lastPos); // find first "non-delimiter"

    while (string::npos != pos || string::npos != lastPos) {
      tokens.push_back(str.substr(lastPos, pos - lastPos)); // add token to the vector<string>
      lastPos = str.find_first_not_of(delimiters, pos); // skip delimiters
      pos = str.find_first_of(delimiters, lastPos); // find next "non-delimiter"
    }
  }

  string getFileNameExtension(const string &str, const string &delimiters) {
    vector<string> tokens;
    tokenizeString(tokens, str, delimiters);
    return tokens[tokens.size()-1];
  }

  bool fileExists(const string& filename) {
    return( fileExists(filename.c_str()) );
  }

  bool fileExists(const char * filename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;
    // Attempt to get the file attributes
    intStat = stat(filename,&stFileInfo);
    if(intStat == 0) {
      // We were able to get the file attributes
      // so the file obviously exists.
      blnReturn = true;
    } else {
      // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.
      blnReturn = false;
    }
    return(blnReturn);
  }

  void resize(int * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(int * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      assert(scalar == NULL);
      scalar = new int[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int * tmp = new int[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }
  
  void resize(int (* &scalar)[2],const int n_new) {
    resize(scalar,0,n_new);
  }

  void resize(int (* &scalar)[2],const int n_old,const int n_new) {
    if (n_old == 0) {
      assert(scalar == NULL);
      scalar = new int[n_new][2];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[2] = new int[n_new][2];
      for (int i = 0; i < n_old; ++i) {
	tmp[i][0] = scalar[i][0];
	tmp[i][1] = scalar[i][1];
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(double * &scalar,const int n_new) {
    resize(scalar,0,n_new);
  }
  
  void resize(double * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      assert(scalar == NULL);
      scalar = new double[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      double * tmp = new double[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(double (*&d3)[3],const int n_new) {
    resize(d3,0,n_new);
  }
    
  void resize(double (*&d3)[3],const int n_old,const int n_new) {
    if (n_old == 0) {
      assert(d3 == NULL);
      d3 = new double[n_new][3];
    }
    else if (n_new > n_old) {
      double (*dtmp)[3] = new double[n_new][3];
      for (int i = 0; i < n_old; ++i) {
	dtmp[i][0] = d3[i][0];
	dtmp[i][1] = d3[i][1];
	dtmp[i][2] = d3[i][2];
      }
      delete[] d3;
      d3 = dtmp;
    }
  }
  
}

