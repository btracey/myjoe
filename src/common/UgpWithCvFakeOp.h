#ifndef UGPWITHCVFAKEOP_H
#define UGPWITHCVFAKEOP_H

#include "UgpWithCvFake.h"
#include "MiscUtils.h"
#include <math.h>

class UgpWithCvFakeOp : public UgpWithCvFake {
  
 private:
  
 protected:
  
  // face-based operators...
  int * cvofa0_i;
  int cvofa0_s;
  int * cvofa0_v;
  
  int * cvofa1_i;
  int cvofa1_s;
  int * cvofa1_v;

  double * cvofa0_value;
  double * cvofa1_value;
  double (*cvofa0_grad)[3];
  double (*cvofa1_grad)[3];
  
  double * fa_alpha;
  double * fa_alpha_no;
  
 public:
  
  // contructor for this class. Needs to be called
  // with a parameter filename to initialize the 
  // parameters (stored in a ParamMap)...

  UgpWithCvFakeOp() {
    
    if (mpi_rank == 0)
      cout << "UgpWithCvFakeOp()"<< endl;
    
    cvofa0_i = cvofa0_v = NULL;
    cvofa0_s = 0;

    cvofa1_i = cvofa1_v = NULL;
    cvofa1_s = 0;

    cvofa0_value = NULL;
    cvofa1_value = NULL;

    cvofa0_grad = NULL;
    cvofa1_grad = NULL;
    
    fa_alpha = NULL;
    fa_alpha_no = NULL;
    
  }
  
  ~UgpWithCvFakeOp() {

    delete[] cvofa0_i; 
    delete[] cvofa0_v;
    
    delete[] cvofa1_i; 
    delete[] cvofa1_v;
    
    delete[] cvofa0_value;
    delete[] cvofa1_value;

    delete[] cvofa0_grad;
    delete[] cvofa1_grad;
    
    delete[] fa_alpha;

    if (fa_alpha_no != NULL)
      unregisterScalar(fa_alpha_no);
    
  }
  
 protected:
  
  void init(const double fa_alpha_coeff) {

    // this init routin etake the fa_alpha_coeff, which multiplies the
    // dissipation term in the operators. a value of 1.0 results in a very low
    // dissipation operator that can handle grid transitions and skewness, although
    // bad skewness combined with anisotropy mayt require a higher value - I had
    // to use 5.0 on one grid, which is terrible, although the advantage of the
    // method is that in good regions, the grid is still dissipation-free.

    // initialize the unstructured grid partition...

    UgpWithCvFake::init();

    // and build our local operators for face-based reconstructions...

    buildOperators();
    checkOperators();
    calcSkewNormStuff(fa_alpha_coeff);

    // HACK
    //simpleOperators();

    /*
    cout << "Warning: going first order..." << endl;
    firstOrderOperators();
    FOR_IFA fa_alpha[ifa] = 1.0;
    */

  }
  
  void initFaAlphaVis() {
    
    registerScalar(fa_alpha_no,"FA_ALPHA",NO_DATA);

    FOR_INO fa_alpha_no[ino] = 0.0;
    FOR_ICV {
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; foc++) {
	const int ifa = faocv_v[foc];
	const int nof_f = noofa_i[ifa];
	const int nof_l = noofa_i[ifa+1]-1;
	for (int nof = nof_f; nof <= nof_l; ++nof) {
	  const int ino = noofa_v[nof];
	  fa_alpha_no[ino] = max(fa_alpha_no[ino],fa_alpha[ifa]);
	}
      }
    }
    updateNoR1(fa_alpha_no,MAX_DATA);
    
  }

private:
  
  void buildOperators();

  void combineOperatorsGR(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
			  double * cvofa0_value_q,double * cvofa1_value_q,
			  double * cvofa0_value_l1,double * cvofa1_value_l1,
			  double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]);
  
  void combineOperatorsLinear(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
			      double * cvofa0_value_q,double * cvofa1_value_q,
			      double * cvofa0_value_l1,double * cvofa1_value_l1,
			      double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]);
  
  void combineOperatorsHam(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
			   double * cvofa0_value_q,double * cvofa1_value_q,
			   double * cvofa0_value_l1,double * cvofa1_value_l1,
			   double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]);
  
  void checkOperators();

  void calcSkewNormStuff(const double fa_alpha_coeff);

  // inline this...
  void addMatrixEntry(double (*DpDt)[3],const double * coeff,const int i_global,const int j_global,
		      const int * nb2ocv_i,int * nb2ocv_v) {
    
    const int icv = i_global - cvora[mpi_rank];
    assert( (icv >= 0)&&(icv < ncv) );
    
    const int n2oc_f = nb2ocv_i[icv];
    const int n2oc_l = nb2ocv_i[icv+1]-1;
    for (int n2oc = n2oc_f; n2oc <= n2oc_l; ++n2oc) {
      if (nb2ocv_v[n2oc] == j_global) {
	FOR_I3 DpDt[n2oc][i] += coeff[i];
	return;
      }
      else if (nb2ocv_v[n2oc] == -1) {
	nb2ocv_v[n2oc] = j_global;
	FOR_I3 DpDt[n2oc][i] += coeff[i];
	return;
      }
    }
    
    // should never get here...
    cout << "Error: nb2ocv_i/v problem." << endl;
    throw(-1);
    
  }

 protected:

  void simpleOperators();

  void firstOrderOperators();

  void firstOrderOperatorsFlaggedCvs();

};

#endif // #ifndef UGPWITHCVFAKEOP_H
