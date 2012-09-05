#ifndef UGPWITHCVFAKE_H
#define UGPWITHCVFAKE_H

#include "UgpWithTools.h"

// for initialization routines...

//#include "CdpFilter.h"
//#include "MshFilter.h"

#include "tc_vec3d.h"


#define CV_G_DATA           10          // add a new registered data type


//#define epsilonSDWLS        1.0e-1      // constant to prevent division by zero

//#define limiterWeight(a) ((a*a))      // van Albada limiter
#define limiterWeight(a) (fabs(a))      // van Leer limiter

#define limSharpen   0.0

#ifndef sign
#define sign(a)       ( ((a) > 0) ? (1.0) : ((a) < 0.0 ? -1.0 : 0.0))
#endif

#define FOR_ICV_G for (int icv = 0; icv < ncv_g; ++icv)
#define FOR_ICV_GF for (int icv = 0; icv < ncv_gf; ++icv)

void ludeco(double(*A)[5], int order);
void lusolv(double(*A)[5], double *c, int order);

class ScalarTranspEq: public Data
{
public:

  double *phi;
  double *diff;

  double turbSchmidtNumber;
  bool convTerm, diffTerm;

  int phiMaxiter;
  double relax;
  double phiZero;
  double resMax, resAve;
  double lowerBound, upperBound;

  ScalarTranspEq(const char * name, int datatype) : Data(name, datatype)
  {
    phi = NULL;
    diff = NULL;

    turbSchmidtNumber = 1.0;

    convTerm = true;
    diffTerm = true;

    phiMaxiter = 500;
    relax = 0.8;
    phiZero = 1.0e-6;

    lowerBound = 0.0;
    upperBound = 1.0e10;
  }

  ~ScalarTranspEq() {}
};



class DoubleTranspVector: public Data
{
public:
  double (*phi)[3];
  DoubleTranspVector(const char * name, int datatype) : Data(name, datatype), phi(0) {}
};

class UgpWithCvFake: public UgpWithTools {

public:   // constructors/destructors

  UgpWithCvFake() {

    if (mpi_rank == 0)
      cout << "UgpWithCvFake()" << endl;

    ncv_g = 0;
    ncv_gf = 0;

    fa_normal = NULL;
    cv_volume = NULL;
    x_fa = NULL;
    x_cv = NULL;

    nbocv_s = 0;
    nbocv_i = nbocv_v = NULL;

    nbocv_lsg_coeff = NULL;
    boundary_fa_lsg_coeff = NULL;
  }

  virtual ~UgpWithCvFake()
  {
    delete []fa_normal;
    delete []cv_volume;
    delete []x_fa;
    delete []x_cv;

    delete []nbocv_i;
    delete []nbocv_v;

    delete []nbocv_lsg_coeff;
    delete []boundary_fa_lsg_coeff;
  }



protected:  // member variables

  enum {NOLIMITER, SDWLS, MINMOD, DOUBLEMINMOD, HARMONIC, SUPERBEE};

  enum LinearSolvers {
    HYPRE_AMG,
    HYPRE_PCG_AMG,
    HYPRE_PCG,
    HYPRE_BCG,
    HYPRE_GMRES_PILUT,
    HYPRE_GMRES_EUCLID,
    PCG_LOCAL,
    PCG,
    BCGSTAB,
    JACOBI,
    PETSC_GMRES,
    NONE // use for skipping pressure solve for scalability testing 
  };

  /**
   * ghost cvs - associated with interprocessor and periodic data,
   * internal cvs are [1:ncv-1], with ghost [ncv:ncv_g-1], and fake [ncv_g:ncv_gf-1]
   */
  int ncv_g,ncv_gf;

  double (*fa_normal)[3];     ///< face normal area weighted
  double *cv_volume;          ///< cell volume
  double (*x_fa)[3];          ///< face center coordinate
  double (*x_cv)[3];          ///< cell center coordinate

  int * nbocv_i;
  int nbocv_s;
  int * nbocv_v;

  /// least squares gradient operators...
  double (*nbocv_lsg_coeff)[3];
  double (*boundary_fa_lsg_coeff)[3];

  // data registration for transport variables ...
  list<ScalarTranspEq> scalarTranspEqList;            ///< list of scalar transport equations
  list<DoubleTranspVector> doubleTranspVectorList;    ///< list of vectors

public:   // member functions


  void init(const int lsg_flag = 0) {

    // standard initialization from a minimum Ugp...

    // call UgpWithTools's init routine...
    UgpWithTools::init();

    // UgpWithCvFake-based specific initializations...
    //addGhostCvs();
    addGhostAndFakeCvs();

    buildNbocv();

    calcGeometry();

    // note: no LSG grads specified here
    if (lsg_flag == 1) {
      calcLsgCoeff();
      checkLsgCoeff();
    }
  }

/*  void initializeFromRestartFile(const string &name) {

    string ext = getFileNameExtension(name, ".");

    if (ext == "cdp")
    {
      CdpFilter cf(name);
      cf.minInitUgp(this);
    }
    else if ((ext == "msh") || (ext == "cas"))
    {
      MshFilter mf(name);
      mf.initUgpMin(this);
    }
    else
    {
      // assume a generic restart...
      readRestart(name);
    }
    
    // calc volume, face vectors and attach ghost and fake cells 
    init(1);
  }*/


  ScalarTranspEq * registerScalarTransport(double *pScal, const char * name, int datatype)
  {
    // check that the name does not conflict with any other registered data...
    list<ScalarTranspEq>::iterator data;
    for (data = scalarTranspEqList.begin(); data != scalarTranspEqList.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
      {
        if (mpi_rank == 0)
          cerr << "Error: transport scalar already registered with name: " << name << endl;
        throw(-1);
      }
    }

    // push scalar into list...
    scalarTranspEqList.push_back(ScalarTranspEq(name, datatype));

    // connect pointer
    ScalarTranspEq *eq = getScalarTransportData(name);
    registerScalar(eq->phi, name, datatype);
    pScal = eq->phi;

    eq = &(scalarTranspEqList.back());
    return eq;
  }

  
  ScalarTranspEq * registerScalarTransport(const char * name, int datatype)
  {
    // check that the name does not conflict with any other registered data...
    list<ScalarTranspEq>::iterator data;
    for (data = scalarTranspEqList.begin(); data != scalarTranspEqList.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
      {
        if (mpi_rank == 0)
          cerr << "Error: transport scalar already registered with name: " << name << endl;
        throw(-1);
      }
    }

    // push it into list...
    scalarTranspEqList.push_back(ScalarTranspEq(name, datatype));
    ScalarTranspEq *transScal = &scalarTranspEqList.back();
    registerScalar(transScal->phi, name, datatype);
    
    return transScal;
  }

  ScalarTranspEq * getScalarTransportData(const string &name)
  {
    return getScalarTransportData(name.c_str());
  }

  ScalarTranspEq * getScalarTransportData(const char * name)
  {
    list<ScalarTranspEq>::iterator data;
    for (data = scalarTranspEqList.begin(); data != scalarTranspEqList.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
        return (&(*data));
    }
    return (NULL);
  }

  void registerTransportVector(const char * name, int datatype)
  {
    // check that the name does not conflict with any other registered data...
    list<DoubleTranspVector>::iterator data;
    for (data = doubleTranspVectorList.begin(); data != doubleTranspVectorList.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
      {
        if (mpi_rank == 0)
          cerr << "Error: transport vector already registered with name: " << name << endl;
        throw(-1);
      }
    }

    // push it into list...
    doubleTranspVectorList.push_back(DoubleTranspVector(name, datatype));

    DoubleTranspVector *transVec = getTranpVectorData(name);
    registerVector(transVec->phi, name, datatype);
  }

  DoubleTranspVector * getTranpVectorData(const char * name)
  {
    list<DoubleTranspVector>::iterator data;
    for (data = doubleTranspVectorList.begin(); data != doubleTranspVectorList.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
        return (&(*data));
    }
    return (NULL);
  }

  /**
   * scalar transport equation -> calc residual
   */
  void calcResidual(double &res_ave, double &res_max, double *phi, double *A, double *rhs) 
  {
    double this_res;
    res_max = -1.0;
    res_ave =  0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      this_res = rhs[icv];

      int noc_f = nbocv_i[icv]; // diagonal
      int noc_l = nbocv_i[icv+1]-1;

      this_res -= A[noc_f]*phi[icv];

      for (int noc = noc_f+1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        this_res -= A[noc]*phi[icv_nbr];
      }
      
      res_ave += fabs(this_res);//*this_res;
      res_max = max(res_max, this_res);
    }
//    res_ave /= (double)(ncv);
  }

  /**
   * vector transport equation -> calc residual
   */
  void calcResidual(double *res_ave, double *res_max, double (*phi)[3], double *A, double (*rhs)[3]) {

    double this_res[3];
    res_max[0] = res_max[1] = res_max[2] = -1.0;
    res_ave[0] = res_ave[1] = res_ave[2] =  0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      for (int i=0; i<3; i++) this_res[i] = rhs[icv][i];

      int noc_f = nbocv_i[icv]; // diagonal
      int noc_l = nbocv_i[icv+1]-1;

      for (int i=0; i<3; i++) this_res[i] -= A[noc_f]*phi[icv][i];
      for (int noc = noc_f+1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        for (int i=0; i<3; i++) this_res[i] -= A[noc]*phi[icv_nbr][i];
      }
      for (int i=0; i<3; i++)  {
        res_ave[i] += this_res[i]*this_res[i];
        res_max[i] = max(res_max[i], this_res[i]);
      }
    }
    for (int i=0; i<3; i++) res_ave[i] /= (double)(ncv);
  }



private:
  
  void addGhostCvs(); // do not use this one any more
  
  void addGhostAndFakeCvs();
  void buildNbocv();
  void calcGeometry();
  void calcLsgCoeff();
  void checkLsgCoeff();

protected:

  /**
   * search algorithm for icvm1,
   * sketch of face and cv alignment
   *
   \verbatim
                    ifa01
          |     |     |     |
          |  X  |  X  |  X  |
          |     |     |     |
           icvm1  icv0  cv1
   \endverbatim
   *  if icvm1 and icv1 are periodic cvs then icvm1=icv0 ?!?!?!? (should be)
   */
  inline int find_ICVm1(int icv0, int ifa01, double maxAngle)
  {
    if (icv0 >= ncv) return -1;  // return (-1) if icv0 is phantom cell

    double nvec1[3], nvec2[3];
    int notfound;

    int foc_f = faocv_i[icv0];
    int foc_l = faocv_i[icv0+1]-1;

    int icvm1 = -1;
    notfound = 1;

    for (int foc=foc_f; foc<=foc_l && notfound; foc++)
    {
      int ifa = faocv_v[foc];
      normVec3d(nvec1, fa_normal[ifa]);
      normVec3d(nvec2, fa_normal[ifa01]);
      double dot = vecDotVec3d(nvec1, nvec2);
      if ((fabs(dot) > maxAngle) && (ifa != ifa01) && (ifa >= nfa_b))
      {
        icvm1 = cvofa[ifa][0];
        notfound = 0;
      }
    }
    return icvm1;
  }

  /**
   * search algorithm for icvp2,
   * sketch of face and cv alignment
   *
   \verbatim
             ifa01
         |     |     |     |
         |  X  |  X  |  X  |
         |     |     |     |
          icv0  icv1  icvp2
   \endverbatim
   *
   *  if icvm1 and icv1 are periodic cvs then icvm1=icv0 ?!?!?!? (should be)
   *
   */
  inline int find_ICVp2(int icv1, int ifa01, double maxAngle)
  {
    if (icv1 >= ncv) return -1;  // return (-1) if icv0 is phantom cell

    double nvec1[3], nvec2[3];
    int notfound;

    int foc_f = faocv_i[icv1];
    int foc_l = faocv_i[icv1+1]-1;

    int icvp2 = -1;
    notfound = 1;

    for (int foc=foc_f; foc<=foc_l && notfound; foc++)
    {
      int ifa = faocv_v[foc];
      normVec3d(nvec1, fa_normal[ifa]);
      normVec3d(nvec2, fa_normal[ifa01]);
      double dot = vecDotVec3d(nvec1, nvec2);
      if ((fabs(dot) > maxAngle) && (ifa != ifa01) && (ifa >= nfa_b))
      {
        icvp2 = cvofa[ifa][1];
        notfound = 0;
      }
    }
    return icvp2;
  }

  /**
   *    noc00 -> diagonal element of the icv0 row
   *    noc11 -> diagonal element of the icv1 row
   *    noc01 -> flux contribution of icv1 to the icv0 balance
   *    noc10 -> flux contribution of icv0 to the icv1 balance
   */
  inline void getImplDependencyIndex(int &noc00, int &noc01, int &noc11, int &noc10, int icv0, int icv1) {
    noc00 = nbocv_i[icv0];
    noc01 = noc00;
    while (nbocv_v[++noc01] != icv1);

    if (icv1 < ncv) {
      noc11 = nbocv_i[icv1];
      noc10 = noc11;
      while (nbocv_v[++noc10] != icv0);
    }
  }

  
  
  /**
   * 
   * 
   *  
   * GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES
   *  
   * 
   *  
   */
  void calcCvScalarGrad(double (*gradPhi)[3], const double *phi, const double *phi_fa, int limiter, double *refValue, double epsilonSDWLS) 
  {
    switch (limiter)
    {
    case NOLIMITER:
      calcCvScalarGrad(gradPhi, phi);//, phi_fa);
      break;

    case SDWLS:
      calcCvScalarGradSDWLS(gradPhi, phi, phi_fa, refValue, epsilonSDWLS);
      break;

    case MINMOD:
    case DOUBLEMINMOD:
    case HARMONIC:
    case SUPERBEE:
      calcCvScalarGrad(gradPhi, phi);//, phi_fa);
      updateCvData(gradPhi, REPLACE_ROTATE_DATA);
      limitCvScalarGrad(gradPhi, limiter);
      break;
    }
    updateCvData(gradPhi, REPLACE_ROTATE_DATA);
  }


  void calcCvVectorGrad(double(*gradPhi)[3][3], const double (*phi)[3], const double (*phi_fa)[3], int limiter, double *refValue, double epsilonSDWLS) 
  {
    switch (limiter)
    {
    case NOLIMITER:
      calcCvVectorGrad(gradPhi, phi);//, phi_fa);
      break;

    case SDWLS:
      calcCvVectorGradSDWLS(gradPhi, phi, phi_fa, refValue, epsilonSDWLS);
      break;

    case MINMOD:
    case DOUBLEMINMOD:
    case HARMONIC:
    case SUPERBEE:
      calcCvVectorGrad(gradPhi, phi);//, phi_fa);
      updateCvData(gradPhi, REPLACE_ROTATE_DATA);
      limitCvVectorGrad(gradPhi, limiter);
      break;
    }
    updateCvData(gradPhi, REPLACE_ROTATE_DATA);
  }
  
  /** 
   * conventional limiters
   */
  double minmod(double a, double b)
  {
    if ((a * b) < 0.)               return 0.;
    else if (fabs(a) < fabs(b))     return a;
    else                            return b;
  }

  double maxmod(double a, double b)
  {
    if ((a * b) < 0.)               return 0.;
    else if (fabs(a) < fabs(b))     return b;
    else                            return a;
  }

  double doubleMinmod(double a, double b)
  {
    if ((a * b) < 0.)               return 0.;
    else                            return sign(a + b) * min(fabs(2. * a), min(fabs(2. * b), fabs(0.5 * (a + b))));
  }

  double harmonic(double a, double b)
  {
    if ((a * b) < 0.)               return 0.;
    else
    {
      if (fabs(a + b) < 1.e-12)     return 0.5 * (a + b);

      return sign(a + b) * min(fabs(0.5 * (a + b)), fabs(2. * a * b / (a + b)));
    }
  }

  double superbee(double a, double b)
  {
    return minmod(maxmod(a, b), minmod(2. * a, 2. * b));
  }
  
  
  
  /** 
   * 
   * calc gradients 
   *
   */
  void calcCvScalarGrad(double(*grad_p)[3], const double *p)
  {
    assert(nbocv_lsg_coeff != NULL);

    for (int icv = 0; icv < ncv; icv++)
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int i = 0; i < 3; i++)
        grad_p[icv][i] = nbocv_lsg_coeff[noc_f][i] * p[icv];

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
          grad_p[icv][i] += nbocv_lsg_coeff[noc][i] * p[icv_nbr];
      }
    }
  }

  /**
   * with boundary conditions p_bfa
   */
  void calcCvScalarGrad(double(*grad_p)[3], const double *p, const double *p_fa) 
  {
    calcCvScalarGrad(grad_p, p);

    assert(boundary_fa_lsg_coeff != NULL);

//    if (p_bfa != NULL)
    for (int ifa = 0; ifa < nfa_b; ifa++) 
    {
      int icv     = cvofa[ifa][0];
      int icvFake = cvofa[ifa][1];
      for (int i=0; i<3; i++)
        grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*p_fa[icvFake];// - p[icv]);
    }
  }

  /**
   * apply conventional limiters
   */
  void limitCvScalarGrad(double (*grad_p)[3], int limiter)
  {
    double maxGrad[3] = {0., 0., 0.};
    double minGrad[3] = {0., 0., 0.};
    
    static double (*tmp_grad)[3] = new double[ncv_g][3];
    
    // grad_p already updated in ghost cells 
    for (int icv=0; icv<ncv_g; icv++)
      for (int i=0; i<3; i++)
        tmp_grad[icv][i] = grad_p[icv][i];
//    updateCvData(tmp_grad, REPLACE_ROTATE_DATA);

    for (int icv=0; icv<ncv; icv++) 
    {
      // init to own cell values
      for (int i=0;i<3;i++) 
        maxGrad[i] = minGrad[i] = tmp_grad[icv][i];

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
        {
          maxGrad[i] = max(maxGrad[i], (1.-limSharpen)*tmp_grad[icv_nbr][i] + limSharpen*tmp_grad[icv][i]);
          minGrad[i] = min(minGrad[i], (1.-limSharpen)*tmp_grad[icv_nbr][i] + limSharpen*tmp_grad[icv][i]);
        }
      }

      switch (limiter)
      {
      case MINMOD:
        for (int i=0; i<3; i++)
          grad_p[icv][i] = minmod(maxGrad[i], minGrad[i]);
        break;
      case DOUBLEMINMOD:
        for (int i=0; i<3; i++)
          grad_p[icv][i] = doubleMinmod(maxGrad[i], minGrad[i]);
        break;
      case HARMONIC:
        for (int i=0; i<3; i++)
          grad_p[icv][i] = harmonic(maxGrad[i], minGrad[i]);
        break;
      case SUPERBEE:
        for (int i=0; i<3; i++)
          grad_p[icv][i] = superbee(maxGrad[i], minGrad[i]);
        break;
      }
    }
  }


  /** 
   * calculate gradients with weights acting as limiter (SDWLS ... solution dependent weighted least square)
   */
  void calcCvScalarGradSDWLS(double(*grad_p)[3], const double *p, const double *p_fa, double *refValue, double epsilonSDWLS) 
  {
    for (int icv = 0; icv < ncv; icv++) 
    {
      double swdx2 = 0.0; // sum weight dx2, etc...
      double swdy2 = 0.0;
      double swdz2 = 0.0;
      double swdxdy = 0.0;
      double swdxdz = 0.0;
      double swdydz = 0.0;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      // skip the diagonal in this loop...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        double dp = p[icv_nbr] - p[icv];
        double weight = 1.0/(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], 1.0e-8));

        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        swdx2 += weight * dx[0] * dx[0];
        swdy2 += weight * dx[1] * dx[1];
        swdz2 += weight * dx[2] * dx[2];
        swdxdy += weight * dx[0] * dx[1];
        swdxdz += weight * dx[0] * dx[2];
        swdydz += weight * dx[1] * dx[2];
      }

      // this cell may also touch some boundary faces that
      // do not have neighbours (no fake cells). Their
      // influence on the gradient has to be considered as well...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv + 1] - 1;
      for (int foc = foc_f; foc <= foc_l; foc++) 
      {
        int ifa = faocv_v[foc];
        if (ifa < nfa_b) 
        {
          // this is a boundary face...
          assert(cvofa[ifa][0] == icv);
          int icvFake = cvofa[ifa][1];

          double dp = 0.0; 
          dp = p[icvFake] - p[icv];
          double weight = 1.0/(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], 1.0e-8));

          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_cv[icvFake][i] - x_cv[icv][i];
         
          swdx2 += weight * dx[0] * dx[0];
          swdy2 += weight * dx[1] * dx[1];
          swdz2 += weight * dx[2] * dx[2];
          swdxdy += weight * dx[0] * dx[1];
          swdxdz += weight * dx[0] * dx[2];
          swdydz += weight * dx[1] * dx[2];
        }
      }

      double denom = 2.0 * swdxdy * swdxdz * swdydz
                          + swdx2 * swdy2 * swdz2
                          - swdx2 * swdydz * swdydz
                          - swdy2 * swdxdz * swdxdz
                          - swdz2 * swdxdy * swdxdy;

      // check non-singular...
      //assert(fabs(denom) > 1.0E-12 * cv_volume[icv] * cv_volume[icv]);
      // now the gradient...
      for (int i = 0; i < 3; i++)
        grad_p[icv][i] = 0.0;
      
      // build gradient...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        double dp = p[icv_nbr] - p[icv];
        double weight = 1.0/(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], 1.0e-8));

        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];
        
        // following coeff's multiply the delta across nbr - cv...
        grad_p[icv][0] += dp * weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                                       + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                                       + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;
        
        grad_p[icv][1] += dp * weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                                       + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                                       + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;
        
        grad_p[icv][2] += dp * weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                                       + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                                       + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;
      }

      // add boundary
      if (p_fa != NULL)
      for (int foc = foc_f; foc <= foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        if (ifa < nfa_b)
        {
          // this is a boundary face...
          assert(cvofa[ifa][0] == icv);
          int icvFake = cvofa[ifa][1];

          double dp = 0.0;
          dp = p[icvFake] - p[icv];
          double weight = 1.0/(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], 1.0e-8));

          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_cv[icvFake][i] - x_cv[icv][i];

          // following coeff's multiply the delta across nbr - cv...
          grad_p[icv][0] += dp * weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                                         + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                                         + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;
          
          grad_p[icv][1] += dp * weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                                         + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                                         + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;
          
          grad_p[icv][2] += dp * weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                                         + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                                         + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;
        }
      }
    }    
  }
  


  void calcCvVectorGradSDWLS(double(*grad_u)[3][3], const double(*u)[3], const double(*u_fa)[3], double *refValue, double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv; icv++) 
    {
      double swdx2 = 0.0; // sum weight dx2, etc...
      double swdy2 = 0.0;
      double swdz2 = 0.0;
      double swdxdy = 0.0;
      double swdxdz = 0.0;
      double swdydz = 0.0;
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      
      // skip the diagonal in this loop...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        const int icv_nbr = nbocv_v[noc];

        double du2 = 0.0;
        for (int i = 0; i < 3; i++) 
        {
          const double du = u[icv_nbr][i] - u[icv][i];
          du2 += du*du;
        }
        double weight = 1.0/(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]);

        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];
                
        swdx2 += weight * dx[0] * dx[0];
        swdy2 += weight * dx[1] * dx[1];
        swdz2 += weight * dx[2] * dx[2];
        swdxdy += weight * dx[0] * dx[1];
        swdxdz += weight * dx[0] * dx[2];
        swdydz += weight * dx[1] * dx[2];
      }
      
      // this cell may also touch some boundary faces that
      // do not have neighbours (no fake cells). Their
      // influence on the gradient has to be considered as well...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv + 1] - 1;
      for (int foc = foc_f; foc <= foc_l; foc++) 
      {
        int ifa = faocv_v[foc];
        if (ifa < nfa_b) 
        {
          // this is a boundary face...
          assert(cvofa[ifa][0] == icv);
          int icvFake = cvofa[ifa][1];

          double du2 = 0.0;
          for (int i = 0; i < 3; i++)
          {
            const double du = u[icvFake][i] - u[icv][i];
            du2 += du*du;
          }
          double weight = 1.0/(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]);

          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_cv[icvFake][i] - x_cv[icv][i];
          
          swdx2  += weight * dx[0] * dx[0];
          swdy2  += weight * dx[1] * dx[1];
          swdz2  += weight * dx[2] * dx[2];
          swdxdy += weight * dx[0] * dx[1];
          swdxdz += weight * dx[0] * dx[2];
          swdydz += weight * dx[1] * dx[2];
        }
      }
      double denom = 2.0 * swdxdy * swdxdz * swdydz
                         + swdx2 * swdy2  * swdz2
                         - swdx2 * swdydz * swdydz
                         - swdy2 * swdxdz * swdxdz
                         - swdz2 * swdxdy * swdxdy;

      // check non-singular...
      //assert(fabs(denom) > 1.0E-12 * cv_volume[icv] * cv_volume[icv]);
      // now the gradient...
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          grad_u[icv][i][j] = 0.0;

      // build gradient...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        double du[3];
        double du2 = 0.0;
        for (int i = 0; i < 3; i++) 
        {
          du[i] = u[icv_nbr][i] - u[icv][i];
          du2 += du[i]*du[i];
        }
        double weight = 1.0/(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]);

        double dx[3], tmp;
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        tmp = weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                      + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                      + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;

        for (int i = 0; i < 3; i++)
          grad_u[icv][i][0] += du[i] * tmp;

        tmp = weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                      + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                      + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;

        for (int i = 0; i < 3; i++)
          grad_u[icv][i][1] += du[i] * tmp;

        tmp = weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                      + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                      + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;

        for (int i = 0; i < 3; i++)
          grad_u[icv][i][2] += du[i] * tmp;
      }

      // add boundary faces 
      if (u_fa != NULL)
      for (int foc = foc_f; foc <= foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        if (ifa < nfa_b)
        {
          // this is a bounray face...
          assert(cvofa[ifa][0] == icv);
          int icvFake = cvofa[ifa][1];
          
          double du[3];
          double du2 = 0.0;
          for (int i = 0; i < 3; i++) 
          {
            du[i] = u[icvFake][i] - u[icv][i];
            du2 += du[i]*du[i];
          }
          double weight = 1.0/(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]);

          double dx[3], tmp;
          for (int i = 0; i < 3; i++)
            dx[i] = x_cv[icvFake][i] - x_cv[icv][i];

          tmp = weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                        + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                        + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;

          for (int i = 0; i < 3; i++)
            grad_u[icv][i][0] += du[i] * tmp;

          tmp = weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                        + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                        + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;

          for (int i = 0; i < 3; i++)
            grad_u[icv][i][1] += du[i] * tmp;

          tmp = weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                        + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                        + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;

          for (int i = 0; i < 3; i++)
            grad_u[icv][i][2] += du[i] * tmp;
        }
      }
    }
  }


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  void calcCvVectorGrad(double(*grad_u)[3][3], const double(*u)[3]) {

    assert(nbocv_lsg_coeff != NULL);

    for (int icv = 0; icv < ncv; icv++)
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          grad_u[icv][i][j] = nbocv_lsg_coeff[noc_f][j] * u[icv][i];

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            grad_u[icv][i][j] += nbocv_lsg_coeff[noc][j] * u[icv_nbr][i];
      }
    }
  }
  
  
  /**
   * with boundary conditions p_bfa
   */
  void calcCvVectorGrad(double(*grad_u)[3][3], const double(*u)[3], const double (*u_fa)[3]) 
  {
    calcCvVectorGrad(grad_u, u);
    
    assert(boundary_fa_lsg_coeff != NULL);
    for (int ifa = 0; ifa < nfa_b; ifa++) 
    {
      int icv = cvofa[ifa][0];
      int icvFake = cvofa[ifa][1];

      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          grad_u[icv][i][j] += boundary_fa_lsg_coeff[ifa][j]*u[icvFake][i];
    }
  }
  
  void limitCvVectorGrad(double (*grad_u)[3][3], int limiter)
  {
    double maxGrad[3][3] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double minGrad[3][3] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    
    static double (*tmp_grad)[3][3] = new double[ncv_g][3][3];
    // grad_u already updated in ghost cells
    for (int icv=0; icv<ncv_g; icv++)
      for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
          tmp_grad[icv][c][i] = grad_u[icv][c][i];
      
    for (int icv=0; icv<ncv; icv++) 
    {
      // init to own cell values
      for (int c=0; c<3; c++)
      for (int i=0; i<3; i++) 
        maxGrad[c][i] = minGrad[c][i] = tmp_grad[icv][c][i];

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;

      for (int noc = noc_f+1; noc<=noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
        {
          maxGrad[c][i] = max(maxGrad[c][i], (1.-limSharpen)*tmp_grad[icv_nbr][c][i] + limSharpen*tmp_grad[icv][c][i]);
          minGrad[c][i] = min(minGrad[c][i], (1.-limSharpen)*tmp_grad[icv_nbr][c][i] + limSharpen*tmp_grad[icv][c][i]);
        }
      }

      switch (limiter)
      {
      case MINMOD:
        for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
          grad_u[icv][c][i] = minmod(maxGrad[c][i], minGrad[c][i]);
        break;
      case DOUBLEMINMOD:
        for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
          grad_u[icv][c][i] = doubleMinmod(maxGrad[c][i], minGrad[c][i]);
        break;
      case HARMONIC:
        for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
          grad_u[icv][c][i] = harmonic(maxGrad[c][i], minGrad[c][i]);
        break;
      case SUPERBEE:
        for (int c=0; c<3; c++)
        for (int i=0; i<3; i++)
          grad_u[icv][c][i] = superbee(maxGrad[c][i], minGrad[c][i]);
        break;
      }
    }
  }
  


  int solveCvScalarJacobi(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double relax,
      const double zero, const int maxiter);

  int solveCvScalarCg(double * phi, double * Ap, double * rhs, const int mode, const double zero, const int maxiter);

  int solveCvScalarBcgstab(double * phi, double * Ap, double * rhs, const int mode, const double zero,
      const int maxiter, char *scalName);

  int solveCvScalarBcgstab(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double zero,
      const int maxiter);

  void matTimesVecOverCVs(double(*res)[5], double(*Ap)[5][5], double(*phi)[5]);
  void UpdateCvDataStateVec(double(*phi)[5]);

  int solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], 
      const double zeroAbs, const double zeroRel, const int maxiter, char *scalarName, int check_interval);


protected:

  void updateCvDataByName(char *name, const int action) {
    double *scal = getScalar(name);
    updateCvData(scal, action);
  }
};

#endif

