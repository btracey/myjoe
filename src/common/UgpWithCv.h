#ifndef UGPWITHCV_H
#define UGPWITHCV_H

#include "UgpWithTools.h"

// for initialization routines...
#include "CdpFilter.h"
#include "MshFilter.h"
#include "tc_vec3d.h"

#define CV_G_DATA           10          // add a new registered data type



// define one of these limiters: 
//#define limiterWeight(a) ((a*a))        // van Albada limiter
#define limiterWeight(a) (fabs(a))        // van Leer limiter

// define clipping value for sdwls gradient calculation of scalars:
#define SDWLSClippingScalar   1.0e-8


#define geomWeighting


// sharpen limiter for conventional limiters (limit gradients afterwards)
// limiterSharpen = 0.0 ... no sharpening
// limiterSharpen = 1.0 ... full sharpening / no limiting 
#define limSharpen   0.0



#ifndef sign
#define sign(a)       ( ((a) > 0) ? (1.0) : ((a) < 0.0 ? -1.0 : 0.0))
#endif

#define GG_WITH_DIST

//#define WITH_FAKE

#define FOR_ICV_G for (int icv = 0; icv < ncv_g; ++icv)
#define FOR_ICV_GF for (int icv = 0; icv < ncv_gf; ++icv)


void ludeco(double(*A)[5], int order);
void lusolv(double(*A)[5], double *c, int order);
void luinv(double(*A)[5], double Ainv[5][5], int order);


class ScalarTranspEq: public Data
{
public:

  double *phi;                            ///< scalar field
  double *phi_bfa;                        ///< scalar boundary values 
  double (*grad_phi)[3];                  ///< scalar gradients
  double *resid;                          ///< residual field

  double *rhophi;                         ///< conserved scalar field (rho * phi)
  double *rhophi_bfa;                     ///< conserved scalar boundary values  (rho * phi)
  double (*grad_rhophi)[3];               ///< conserved scalar gradients (rho * phi)

  double *diff;                           ///< scalar diffusion coefficient D ... d/dxi(D dphi/dxi)

  double *dpress_dphi;                    ///< derivative of pressure with respect to scalar, used for coupled solver

  double turbSchmidtNumber;
  bool convTerm, diffTerm;

  int phiMaxiter;
  double relax;
  double phiZero, phiZeroRel;
  double resMax, resAve;
  double lowerBound, upperBound;
  string reconstruction;                  ///< determine if conservative or standard reconstruction is used ("CONSERVATIVE" or "STANDARD")
  string coupling;                        ///< determine if scalar is solved coupled or uncoupled from NSE ("COUPLED" or "UNCOUPLED")

  ScalarTranspEq(const char * name, int datatype) : Data(name, datatype)
  {
    phi      = NULL;
    phi_bfa  = NULL;
    grad_phi = NULL;
    resid    = NULL;

    rhophi      = NULL;
    rhophi_bfa  = NULL;
    grad_rhophi = NULL;

    diff = NULL;

    dpress_dphi = NULL;
    
    turbSchmidtNumber = 1.0;

    convTerm = true;
    diffTerm = true;

    phiMaxiter = 100;
    relax = 1.0;
    phiZero = 1.0e-6;
    phiZeroRel = 1.0e-3;

    lowerBound = 0.0;
    upperBound = 1.0e10;

    reconstruction = "STANDARD";
    coupling = "UNCOUPLED";
  }

  ~ScalarTranspEq() {}
};



class DoubleTranspVector: public Data
{
public:
  double (*phi)[3];
  DoubleTranspVector(const char * name, int datatype) : Data(name, datatype), phi(0) {}
};



class UgpWithCv: public UgpWithTools {

private:

  list<Prcomm> cvPrcommList;

  Prcomm * getCvPrcomm(int rank) {
    list<Prcomm>::iterator prcomm = cvPrcommList.begin();
    while (prcomm != cvPrcommList.end()) {
      if (prcomm->getNbrRank() == rank)
        return (&(*prcomm));
      prcomm++;
    }
    cvPrcommList.push_back(Prcomm(rank));
    return (&(cvPrcommList.back()));
  }

public:   // constructors/destructors

  UgpWithCv() {

    if (mpi_rank == 0)
      cout << "UgpWithCv()" << endl;

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

  virtual ~UgpWithCv()
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

  enum {GRAD_GG, GRAD_LS, GRAD_SDWLS}; // gradient reconstruction (Least square, Green-Gauss, or solution depended least square)
  enum {NOLIMITER, BARTH_JESPERSEN_MOD};

  enum LinearSolvers {BCGSTAB, PETSC_GMRES, LUSGS, BCGSTAB_LINELET};


  /**
   * ghost cvs - associated with interprocessor and periodic data,
   * internal cvs are [1:ncv-1], with ghost [ncv:ncv_g-1]
   */
  int ncv_g;
  int ncv_gf;

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
  vector<ScalarTranspEq> scalarTranspEqVector;                               ///< vector of scalar transport equations to be treated segregated
  typedef vector<ScalarTranspEq>::iterator ScalarTranspEqIterator;           ///< iterator for vector of scalar transport equations
  list<DoubleTranspVector> doubleTranspVectorList;                           ///< list of vectors

	// (begin) Linelet preconditioner
	bool *Linelet;
	vector<int> *linelet_cv;
	int nSources;
	// (end) Linelet preconditioner

	// (begin) Multigrid variables
  int nbocv_s_mgLevel1;
	int ncv_g_mgLevel1;
	int ncv_gg_mgLevel1;
	double (*x_cv_mgLevel1)[3];
	double (*x_fa_mgLevel1)[3];
	double (*fa_normal_mgLevel1)[3];
	int ncv_ggf_mgLevel1;
	int nfa_b2_mgLevel1;
	int ncv_ggff_mgLevel1;
	double *cv_volume_mgLevel1;
	int *nbocv_i_mgLevel1;
  int *nbocv_v_mgLevel1;
	int **Children_CV_mgLevel1;
	int *nChildren_CV_mgLevel1;
	int *Parent_CV;
	bool *Agglomerate;
	int * d_interpI1;
	double * d_interpR1;
	double (*d_interpR2)[3];
	double (*d_interpR3)[3][3];
	
	int nbocv_s_mgLevel2;
	int ncv_g_mgLevel2;
	int ncv_gg_mgLevel2;
	double (*x_cv_mgLevel2)[3];
	double (*x_fa_mgLevel2)[3];
	double (*fa_normal_mgLevel2)[3];
	int ncv_ggf_mgLevel2;
	int nfa_b2_mgLevel2;
	int ncv_ggff_mgLevel2;
	double *cv_volume_mgLevel2;
	int *nbocv_i_mgLevel2;
  int *nbocv_v_mgLevel2;
	int **Children_CV_mgLevel2;
	int *nChildren_CV_mgLevel2;
	int *Parent_CV_mgLevel1;
	bool *Agglomerate_mgLevel1;
	
	int nbocv_s_mgLevel3;
	int ncv_g_mgLevel3;
	int ncv_gg_mgLevel3;
	double (*x_cv_mgLevel3)[3];
	double (*x_fa_mgLevel3)[3];
	double (*fa_normal_mgLevel3)[3];
	int ncv_ggf_mgLevel3;
	int nfa_b2_mgLevel3;
	int ncv_ggff_mgLevel3;
	double *cv_volume_mgLevel3;
	int *nbocv_i_mgLevel3;
  int *nbocv_v_mgLevel3;
	int **Children_CV_mgLevel3;
	int *nChildren_CV_mgLevel3;
	int *Parent_CV_mgLevel2;
	bool *Agglomerate_mgLevel2;
	
	// (end) Multigrid variables
	
public:   // member functions



  void init(const int lsg_flag = 0) {

    // standard initialization from a minimum Ugp...

    // call UgpWithTools's init routine...
    UgpWithTools::init();

    // UgpWithCv-based specific initializations...

    // complete csr structure...
    MPI_Allgather(&ncv,1,MPI_INT,&(cvora[1]),1,MPI_INT,mpi_comm);
    cvora[0] = 0;
    for (int i = 0; i < mpi_size; i++)
      cvora[i+1] += cvora[i];


#ifndef WITH_FAKE
    addGhostCvs();
    buildNbocv();
    calcGeometry();
#else
    addGhostAndFakeCvs();
    buildNbocvFake();
    calcGeometryFake();
#endif

    // note: no LSG grads specified here
    if (lsg_flag) {
      calcLsgCoeff();
      checkLsgCoeff();
    }

  }

  void initializeFromRestartFile(const string &name) {

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
        
    init();
    

#ifndef WITH_FAKE
    calcLsgCoeff();
//    checkLsgCoeff();
#else
    calcGradCoeffHam();
    checkGradCoeffHam();
#endif

  }

  void checkGradCoeffHam()
  {
    if (mpi_rank == 0)
      cout << " > checkGradCoeffHam(): ";

    const double phi0 = 1.01;
    const double phix = 2.012;
    const double phiy = 3.0123;
    const double phiz = 4.01234;

    double * phi = new double[ncv_gf];
    FOR_ICV_GF phi[icv] = phi0 + phix*x_cv[icv][0] + phiy*x_cv[icv][1] + phiz*x_cv[icv][2];

    double (*grad_phi)[3] = new double[ncv][3];

    calcCvScalarGradHam(grad_phi,phi);

    // what to see where it is failing?...
    //FOR_ICV cv_flag[icv] = 0;

    double my_max_error = 0.0;
    FOR_ICV {
      my_max_error = max( my_max_error, fabs(grad_phi[icv][0]-phix) );
      my_max_error = max( my_max_error, fabs(grad_phi[icv][1]-phiy) );
      my_max_error = max( my_max_error, fabs(grad_phi[icv][2]-phiz) );

      //if ( (fabs(grad_phi[icv][0]-phix) > 1.0E-6) ||
      //   (fabs(grad_phi[icv][1]-phiy) > 1.0E-6) ||
      //   (fabs(grad_phi[icv][2]-phiz) > 1.0E-6) )
      //cv_flag[icv] = 1;

    }

    //writeFlaggedCvsTecplot("test.dat");
    //throw(-1);

    double max_error;
    MPI_Reduce(&my_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0)
      cout << "max_error (should be small): " << max_error <<endl;

    delete[] phi;
    delete[] grad_phi;

  }


  void calcCvScalarGradHam(double (*dphidxi)[3], const double * phi, const double * phi_fa) {

    // populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
    // to ghost cvs, and there is nothing set in fake data (if present)...

    FOR_ICV {
      int noc_f = nbocv_i[icv];
      FOR_I3 dphidxi[icv][i] = nbocv_lsg_coeff[noc_f][i]*phi[icv];
      int noc_l = nbocv_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
        int icv_nbr = nbocv_v[noc];
        FOR_I3 dphidxi[icv][i] += nbocv_lsg_coeff[noc][i]*phi[icv_nbr];
      }
    }

    if (phi_fa != NULL)
    FOR_IFA_B {
      int icv0 = cvofa[ifa][0];
      FOR_I3 dphidxi[icv0][i] += boundary_fa_lsg_coeff[ifa][i]*phi_fa[ifa];
    }
  }

  void calcCvScalarGradHam(double (*dphidxi)[3], const double * phi) {

    // populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
    // to ghost cvs, and there is nothing set in fake data (if present)...

    FOR_ICV {
      int noc_f = nbocv_i[icv];
      FOR_I3 dphidxi[icv][i] = nbocv_lsg_coeff[noc_f][i]*phi[icv];
      int noc_l = nbocv_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
        int icv_nbr = nbocv_v[noc];
        FOR_I3 dphidxi[icv][i] += nbocv_lsg_coeff[noc][i]*phi[icv_nbr];
      }
    }

    FOR_IFA_B {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      FOR_I3 dphidxi[icv0][i] += boundary_fa_lsg_coeff[ifa][i]*phi[icv1];
    }
  }

  void calcCvVectorGradHam(double (*duidxj)[3][3], const double (*u)[3], const double (*u_fa)[3]) {

    // populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
    // to ghost cvs, and there is nothing set in fake data (if present)...

    FOR_ICV {
      int noc_f = nbocv_i[icv];
      FOR_I3 FOR_J3 duidxj[icv][i][j] = nbocv_lsg_coeff[noc_f][j]*u[icv][i];
      int noc_l = nbocv_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
        int icv_nbr = nbocv_v[noc];
        FOR_I3 FOR_J3 duidxj[icv][i][j] += nbocv_lsg_coeff[noc][j]*u[icv_nbr][i];
      }
    }

    if (u_fa != NULL)
    FOR_IFA_B {
      int icv0 = cvofa[ifa][0];
      FOR_I3 FOR_J3 duidxj[icv0][i][j] += boundary_fa_lsg_coeff[ifa][j]*u_fa[ifa][i];
    }

  }

  void calcCvVectorGradHam(double (*duidxj)[3][3], const double (*u)[3]) {

    // populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
    // to ghost cvs, and there is nothing set in fake data (if present)...

    FOR_ICV {
      int noc_f = nbocv_i[icv];
      FOR_I3 FOR_J3 duidxj[icv][i][j] = nbocv_lsg_coeff[noc_f][j]*u[icv][i];
      int noc_l = nbocv_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
        int icv_nbr = nbocv_v[noc];
        FOR_I3 FOR_J3 duidxj[icv][i][j] += nbocv_lsg_coeff[noc][j]*u[icv_nbr][i];
      }
    }

    FOR_IFA_B {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      FOR_I3 FOR_J3 duidxj[icv0][i][j] += boundary_fa_lsg_coeff[ifa][j]*u[icv1][i];
    }
  }

  ScalarTranspEq * registerScalarTransport(const char * name, int datatype)
  {
    // check that the name does not conflict with any other registered data...
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
      {
        if (mpi_rank == 0)
          cerr << "Error: transport scalar already registered with name: " << name << endl;
        throw(-1);
      }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // push it into vector,
    // ATTENTION: Reallocations of a stl vector invalidate all previously obtained iterators, references and pointers.
    // see: http://www.cplusplus.com/reference/stl/vector/push_back/
    scalarTranspEqVector.push_back(ScalarTranspEq(name, datatype));

    //-----------------------------------------------------------------------------------------------------------------
    // therefore, we need to re-connect the doubleScalarListPointer (phi) to the pointer sitting in ScalarList (**ptr)
    for (int i = 0; i < scalarTranspEqVector.size()-1; i++)
    {
      char residscal[30] = "resid_";
      strcat(residscal,scalarTranspEqVector[i].getName());
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++)
      {
        if (strcmp(scalarTranspEqVector[i].getName(), data->getName()) == 0)
          data->ptr = &scalarTranspEqVector[i].phi; //data->connectPointers(scalarTranspEqVector[i].phi);
        if (strcmp(residscal, data->getName()) == 0)
          data->ptr = &scalarTranspEqVector[i].resid;
      }
    }

    char residscal[30] = "resid_";
    ScalarTranspEq *transScal = &(scalarTranspEqVector.back());

    registerScalar(transScal->phi, name, datatype);
    registerScalar(transScal->resid, strcat(residscal,name), datatype);

    return transScal;
  }

  ScalarTranspEq * getScalarTransportData(const string &name)
  {
    return getScalarTransportData(name.c_str());
  }

  ScalarTranspEq * getScalarTransportData(const char * name)
  {
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
        return (&(*data));
    }
    return (NULL);
  }

  int getScalarTransportIndex(const string &name)
  {
    return getScalarTransportIndex(name.c_str());
  }

  int getScalarTransportIndex(const char * name)
  {
    for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    {
      if (strcmp(name, scalarTranspEqVector[iScal].getName()) == 0)
        return iScal;
    }
    return (-1);
  }

  int getScalarTransportCoupledIndex(const string &name)
  {
    return getScalarTransportIndex(name.c_str());
  }

  int getScalarTransportCoupledIndex(const char * name)
  {
    int iScalCoupled = 0;
    for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    {
      if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      {
        if (strcmp(name, scalarTranspEqVector[iScal].getName()) == 0)
          return iScalCoupled;
        iScalCoupled++;
      }
    }
    return (-1);
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
      
      res_ave += fabs(this_res);
      res_max = max(res_max, this_res);
    }
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

  void addGhostCvs();
  void addGhostAndFakeCvs();
  void buildNbocv();
  void buildNbocvFake();
  void calcGeometry();
  void calcGeometryFake();
  void calcLsgCoeff();
  void calcGradCoeffHam();
  void checkLsgCoeff();


  int getNextNodeOfFaceCCW(const int ifa,const int ino,const int icv)
  {
    int nof_f = noofa_i[ifa];
    int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof <= nof_l; ++nof)
    {
      int this_ino = noofa_v[nof];
      if (this_ino == ino) {
        if (cvofa[ifa][0] == icv) {
          if (nof == nof_l)
            return( noofa_v[nof_f] );
          else
            return( noofa_v[nof+1] );
        }
        else {
          assert( cvofa[ifa][1] == icv );
          if (nof == nof_f)
            return( noofa_v[nof_l] );
          else
            return( noofa_v[nof-1] );
        }
      }
    }
    cerr << "Error: could not find node " << ino << " in face " << ifa << endl;
    throw(-1);
  }

  int getPrevNodeOfFaceCCW(const int ifa,const int ino,const int icv)
  {
    int nof_f = noofa_i[ifa];
    int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof <= nof_l; ++nof) {
      int this_ino = noofa_v[nof];
      if (this_ino == ino) {
        if (cvofa[ifa][0] == icv) {
          if (nof == nof_f)
            return( noofa_v[nof_l] );
          else
            return( noofa_v[nof-1] );
        }
        else {
          assert( cvofa[ifa][1] == icv );
          if (nof == nof_l)
            return( noofa_v[nof_f] );
          else
            return( noofa_v[nof+1] );
        }
      }
    }

    cerr << "Error: could not find node " << ino << " in face " << ifa << endl;
    throw(-1);
  }

protected:


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

	inline void getImplDependencyIndex_mg(int &noc00, int &noc01, int &noc11, int &noc10, int icv0, int icv1, int iMesh) {
		if (iMesh == 1) {
			noc00 = nbocv_i_mgLevel1[icv0];
			noc01 = noc00;
			while (nbocv_v_mgLevel1[++noc01] != icv1);
			
			if (icv1 < ncv_mgLevel1) {
				noc11 = nbocv_i_mgLevel1[icv1];
				noc10 = noc11;
				while (nbocv_v_mgLevel1[++noc10] != icv0);
			}
		}
		if (iMesh == 2) {
			noc00 = nbocv_i_mgLevel2[icv0];
			noc01 = noc00;
			while (nbocv_v_mgLevel2[++noc01] != icv1);
			
			if (icv1 < ncv_mgLevel2) {
				noc11 = nbocv_i_mgLevel2[icv1];
				noc10 = noc11;
				while (nbocv_v_mgLevel2[++noc10] != icv0);
			}
		}
		if (iMesh == 3) {
			noc00 = nbocv_i_mgLevel3[icv0];
			noc01 = noc00;
			while (nbocv_v_mgLevel3[++noc01] != icv1);
			
			if (icv1 < ncv_mgLevel3) {
				noc11 = nbocv_i_mgLevel3[icv1];
				noc10 = noc11;
				while (nbocv_v_mgLevel3[++noc10] != icv0);
			}
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
  void calcCvScalarGrad(double (*gradPhi)[3], const double *phi, const double *phi_fa, 
      int type, int limiter, double *refValue, double epsilonSDWLS, double *alpha = NULL) 
  {
    switch (type)
    {
    case GRAD_GG:
#ifndef WITH_FAKE
      calcCvScalarGradGreenGauss(gradPhi, phi, phi_fa);
#else
      calcCvScalarGradHam(gradPhi, phi, phi_fa);
#endif
      if (limiter > NOLIMITER)
      {
        updateCvData(gradPhi, REPLACE_ROTATE_DATA);
        limitCvScalarGradBJ(gradPhi, phi, refValue, epsilonSDWLS, alpha);
      }
      break;

    case GRAD_LS:
      calcCvScalarGrad(gradPhi, phi, phi_fa);
      if (limiter > NOLIMITER)
      {
        updateCvData(gradPhi, REPLACE_ROTATE_DATA);
        limitCvScalarGradBJ(gradPhi, phi, refValue, epsilonSDWLS, alpha);
      }
      break;

    case GRAD_SDWLS:
      calcCvScalarGradSDWLS(gradPhi, phi, phi_fa, refValue, epsilonSDWLS);
      break;
    }

    updateCvData(gradPhi, REPLACE_ROTATE_DATA);
  }
  
  void calcCvVectorGrad(double(*gradPhi)[3][3], const double (*phi)[3], const double (*phi_fa)[3], 
      int type, int limiter, double *refValue, double epsilonSDWLS) 
  {
    switch (type)
    {
    case GRAD_GG:
#ifndef WITH_FAKE
      calcCvVectorGradGreenGauss(gradPhi, phi, phi_fa);
#else
      calcCvVectorGradHam(gradPhi, phi, phi_fa);
#endif

      if (limiter > NOLIMITER)
      {
        updateCvData(gradPhi, REPLACE_ROTATE_DATA);
        limitCvVectorGradBJ(gradPhi, phi, refValue, epsilonSDWLS);
      }
      break;

    case GRAD_LS:
      calcCvVectorGrad(gradPhi, phi, phi_fa);
      if (limiter > NOLIMITER)
      {
        updateCvData(gradPhi, REPLACE_ROTATE_DATA);
        limitCvVectorGradBJ(gradPhi, phi, refValue, epsilonSDWLS);
      }
      break;

    case GRAD_SDWLS:
      calcCvVectorGradSDWLS(gradPhi, phi, phi_fa, refValue, epsilonSDWLS);
      break;
    }

    updateCvData(gradPhi, REPLACE_ROTATE_DATA);
  }
  
  
  

//#define BJ
  /**
   * apply smooth barth jespersen limiters
   */
  void limitCvScalarGradBJ(double (*grad_p)[3], const double *phi, double *refValue, double epsilonSDWLS, double *alpha)
  {
    for (int icv = 0; icv < ncv; icv++) 
    {
      double phiMax = phi[icv];
      double phiMin = phiMax;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      // skip the diagonal in this loop...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) 
      {
        int icv_nbr = nbocv_v[noc];
        phiMax = max(phiMax, phi[icv_nbr]);
        phiMin = min(phiMin, phi[icv_nbr]);
      }

      double alfa = 1.0;

      if ((phiMax-phiMin) > 1.0e-12)
      {
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv + 1] - 1;
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];

          if (ifa >= nfa_b)
          {
            double r0[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
            double phifa = phi[icv] + vecDotVec3d(r0, grad_p[icv]);
#ifdef BJ
            if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
            else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
#else 
            double dp;
            double eps2 = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
            eps2 = pow(eps2, 2.0);
            double dm = phifa - phi[icv];

            if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
            else                                dp = phiMin - phi[icv];

            alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
          }
        }
      }
      else alfa = 0.0;

      for (int i=0; i<3; i++)
        grad_p[icv][i] *= alfa;

        // save alfa in alpha, if defined
      if (alpha != NULL)    alpha[icv] = alfa;
    }
  }

  /**
   * apply smooth barth jespersen limiters
   */
  void limitCvVectorGradBJ(double (*grad_p)[3][3], const double (*phi)[3], double *refValue, double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv; icv++) 
    {
      for (int i=0; i<3; i++)
      {
        double phiMax = phi[icv][i];
        double phiMin = phiMax;
  
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;
        // skip the diagonal in this loop...
        for (int noc = noc_f + 1; noc <= noc_l; noc++) 
        {
          int icv_nbr = nbocv_v[noc];
          phiMax = max(phiMax, phi[icv_nbr][i]);
          phiMin = min(phiMin, phi[icv_nbr][i]);
        }

        double alfa = 1.0;

        if ((phiMax-phiMin) > 1.0e-8)
        {
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv + 1] - 1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];

            if (ifa >= nfa_b)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
              double phifa = phi[icv][i] + vecDotVec3d(r0, grad_p[icv][i]);

#ifdef BJ
              if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
              else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
#else
              double dp;
              double eps2 = pow(epsilonSDWLS*refValue[icv], 2.0);
              double dm = phifa - phi[icv][i];

              if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
              else                                dp = phiMin - phi[icv][i];

              alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
            }
          }
        }
        else alfa = 0.0;

        for (int j=0; j<3; j++)
          grad_p[icv][i][j] *= alfa;
      }
    }
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

  void calcCvVectorGrad(double(*grad_u)[3][3], const double(*u)[3])
  {
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
  void calcCvScalarGrad(double(*grad_p)[3], const double *p, const double *p_fa) 
  {
    calcCvScalarGrad(grad_p, p);

    assert(boundary_fa_lsg_coeff != NULL);

    if (p_fa != NULL)
    for (int ifa = 0; ifa < nfa_b; ifa++) {
      int icv = cvofa[ifa][0];
      for (int i=0; i<3; i++)
        grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*p_fa[ifa];
    }
  }

  void calcCvVectorGrad(double(*grad_u)[3][3], const double(*u)[3], const double (*u_fa)[3])
  {
    calcCvVectorGrad(grad_u, u);

    assert(boundary_fa_lsg_coeff != NULL);

    if (u_fa != NULL)
    for (int ifa = 0; ifa < nfa_b; ifa++) {
      int icv = cvofa[ifa][0];
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          grad_u[icv][i][j] += boundary_fa_lsg_coeff[ifa][j]*u_fa[ifa][i];
    }
  }


  /** 
   * calculate gradients with Green-Gauss
   */


  void calcCvScalarGradGreenGauss(double(*grad_p)[3], const double *p, const double *p_fa) 
  {      
    for (int icv=0; icv<ncv; ++icv)
      for (int i = 0; i < 3; ++i)
        grad_p[icv][i] = 0.0;  

    for (int ifa = nfa_b; ifa < nfa; ++ifa)
    {
      // these faces have both icv0 and icv1 local...
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
#ifdef GG_WITH_DIST
      double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
#else
      double w0 = 0.5, w1 = 0.5;
#endif

      double faceValue = (w1*p[icv0]+w0*p[icv1])/(w0+w1);

      for (int i = 0; i < 3; ++i) 
      {
        grad_p[icv0][i] += faceValue*fa_normal[ifa][i];
        grad_p[icv1][i] -= faceValue*fa_normal[ifa][i];
      }
    }

    if (p_fa != NULL)
    for (int ifa = 0; ifa < nfa_b; ++ifa)
    {
      const int icv = cvofa[ifa][0];
      for (int i = 0; i < 3; ++i)
        grad_p[icv][i] += p_fa[ifa]*fa_normal[ifa][i];
    }

    for (int icv = 0; icv < ncv; ++icv)
      for (int i = 0; i < 3; ++i)
        grad_p[icv][i] /= cv_volume[icv];
  }
  
  /** 
   * calculate gradients with Green-Gauss
   */
  void calcCvVectorGradGreenGauss(double(*grad_phi)[3][3], const double (*phi)[3], const double (*phi_fa)[3]) 
  {      
    for (int icv=0; icv<ncv; ++icv)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          grad_phi[icv][i][j] = 0.0;  

    for (int ifa = nfa_b; ifa < nfa; ++ifa)
    {
      // these faces have both icv0 and icv1 local...
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
#ifdef GG_WITH_DIST
      double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
#else
      double w0 = 0.5, w1 = 0.5;
#endif

      for (int i = 0; i < 3; ++i)
      {
        double faceValue = (w1*phi[icv0][i] + w0*phi[icv1][i])/(w0+w1);

        for (int j = 0; j < 3; ++j)
        {
          grad_phi[icv0][i][j] += faceValue*fa_normal[ifa][j];
          grad_phi[icv1][i][j] -= faceValue*fa_normal[ifa][j];
        }
      }
    }

    if (phi_fa != NULL)
    for (int ifa = 0; ifa < nfa_b; ++ifa)
    {
      const int icv = cvofa[ifa][0];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          grad_phi[icv][i][j] += phi_fa[ifa][i]*fa_normal[ifa][j];
    }

    for (int icv = 0; icv < ncv; ++icv)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          grad_phi[icv][i][j] /= cv_volume[icv];
  }
  
	
	/** 
   * calculate gradients with Green-Gauss for the multigrid levels
   */
	void calcCvScalarGradGreenGauss_mg(double(*grad_p)[3], const double *p, const double *p_fa, int iMesh) {
		if (iMesh == 1) {
			for (int icv=0; icv<ncv_mgLevel1; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ++ifa)
			{
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel1[ifa][0];
				const int icv1 = cvofa_mgLevel1[ifa][1];				
				double faceValue = 0.5*(p[icv0]+p[icv1]);
				
				for (int i = 0; i < 3; ++i) 
				{
					grad_p[icv0][i] += faceValue*fa_normal_mgLevel1[ifa][i];
					grad_p[icv1][i] -= faceValue*fa_normal_mgLevel1[ifa][i];
				}
			}
			
			if (p_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel1; ++ifa)
				{
					const int icv = cvofa_mgLevel1[ifa][0];
					for (int i = 0; i < 3; ++i)
						grad_p[icv][i] += p_fa[ifa]*fa_normal_mgLevel1[ifa][i];
				}
			
			for (int icv = 0; icv < ncv_mgLevel1; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] /= cv_volume_mgLevel1[icv];
		}
		
		if (iMesh == 2) {
			for (int icv=0; icv<ncv_mgLevel2; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ++ifa)
			{
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel2[ifa][0];
				const int icv1 = cvofa_mgLevel2[ifa][1];				
				double faceValue = 0.5*(p[icv0]+p[icv1]);
				
				for (int i = 0; i < 3; ++i) 
				{
					grad_p[icv0][i] += faceValue*fa_normal_mgLevel2[ifa][i];
					grad_p[icv1][i] -= faceValue*fa_normal_mgLevel2[ifa][i];
				}
			}
			
			if (p_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel2; ++ifa)
				{
					const int icv = cvofa_mgLevel2[ifa][0];
					for (int i = 0; i < 3; ++i)
						grad_p[icv][i] += p_fa[ifa]*fa_normal_mgLevel2[ifa][i];
				}
			
			for (int icv = 0; icv < ncv_mgLevel2; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] /= cv_volume_mgLevel2[icv];		
		}
		
		if (iMesh == 3) {
			for (int icv=0; icv<ncv_mgLevel3; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ++ifa)
			{
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel3[ifa][0];
				const int icv1 = cvofa_mgLevel3[ifa][1];				
				double faceValue = 0.5*(p[icv0]+p[icv1]);
				
				for (int i = 0; i < 3; ++i) 
				{
					grad_p[icv0][i] += faceValue*fa_normal_mgLevel3[ifa][i];
					grad_p[icv1][i] -= faceValue*fa_normal_mgLevel3[ifa][i];
				}
			}
			
			if (p_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel3; ++ifa)
				{
					const int icv = cvofa_mgLevel3[ifa][0];
					for (int i = 0; i < 3; ++i)
						grad_p[icv][i] += p_fa[ifa]*fa_normal_mgLevel3[ifa][i];
				}
			
			for (int icv = 0; icv < ncv_mgLevel3; ++icv)
				for (int i = 0; i < 3; ++i)
					grad_p[icv][i] /= cv_volume_mgLevel3[icv];		
		}
		
  }

	
	/** 
   * calculate gradients with Green-Gauss for the multigrid levels
   */
  void calcCvVectorGradGreenGauss_mg(double(*grad_phi)[3][3], const double (*phi)[3], const double (*phi_fa)[3], int iMesh) 
  {
		
		if (iMesh == 1) {
			for (int icv=0; icv<ncv_mgLevel1; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ++ifa) {
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel1[ifa][0];
				const int icv1 = cvofa_mgLevel1[ifa][1];
				
				for (int i = 0; i < 3; ++i)
				{
					double faceValue = 0.5*(phi[icv0][i] + phi[icv1][i]);
					
					for (int j = 0; j < 3; ++j)
					{
						grad_phi[icv0][i][j] += faceValue*fa_normal_mgLevel1[ifa][j];
						grad_phi[icv1][i][j] -= faceValue*fa_normal_mgLevel1[ifa][j];
					}
				}
			}
			
			if (phi_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel1; ++ifa)
				{
					const int icv = cvofa_mgLevel1[ifa][0];
					for (int i = 0; i < 3; ++i)
						for (int j = 0; j < 3; ++j)
							grad_phi[icv][i][j] += phi_fa[ifa][i]*fa_normal_mgLevel1[ifa][j];
				}
			
			for (int icv = 0; icv < ncv_mgLevel1; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] /= cv_volume_mgLevel1[icv];
			
		}
		
		if (iMesh == 2) {
			for (int icv=0; icv<ncv_mgLevel2; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ++ifa) {
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel2[ifa][0];
				const int icv1 = cvofa_mgLevel2[ifa][1];
				
				for (int i = 0; i < 3; ++i)
				{
					double faceValue = 0.5*(phi[icv0][i] + phi[icv1][i]);
					
					for (int j = 0; j < 3; ++j)
					{
						grad_phi[icv0][i][j] += faceValue*fa_normal_mgLevel2[ifa][j];
						grad_phi[icv1][i][j] -= faceValue*fa_normal_mgLevel2[ifa][j];
					}
				}
			}
			
			if (phi_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel2; ++ifa)
				{
					const int icv = cvofa_mgLevel2[ifa][0];
					for (int i = 0; i < 3; ++i)
						for (int j = 0; j < 3; ++j)
							grad_phi[icv][i][j] += phi_fa[ifa][i]*fa_normal_mgLevel2[ifa][j];
				}
			
			for (int icv = 0; icv < ncv_mgLevel2; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] /= cv_volume_mgLevel2[icv];
			
		}
		
		if (iMesh == 3) {
			for (int icv=0; icv<ncv_mgLevel3; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] = 0.0;  
			
			for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ++ifa) {
				// these faces have both icv0 and icv1 local...
				const int icv0 = cvofa_mgLevel3[ifa][0];
				const int icv1 = cvofa_mgLevel3[ifa][1];
				
				for (int i = 0; i < 3; ++i)
				{
					double faceValue = 0.5*(phi[icv0][i] + phi[icv1][i]);
					
					for (int j = 0; j < 3; ++j)
					{
						grad_phi[icv0][i][j] += faceValue*fa_normal_mgLevel3[ifa][j];
						grad_phi[icv1][i][j] -= faceValue*fa_normal_mgLevel3[ifa][j];
					}
				}
			}
			
			if (phi_fa != NULL)
				for (int ifa = 0; ifa < nfa_b_mgLevel3; ++ifa)
				{
					const int icv = cvofa_mgLevel3[ifa][0];
					for (int i = 0; i < 3; ++i)
						for (int j = 0; j < 3; ++j)
							grad_phi[icv][i][j] += phi_fa[ifa][i]*fa_normal_mgLevel3[ifa][j];
				}
			
			for (int icv = 0; icv < ncv_mgLevel3; ++icv)
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						grad_phi[icv][i][j] /= cv_volume_mgLevel3[icv];
			
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
        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        double geomW = 1.0;

#ifdef geomWeighting
        double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        //geomW = pow(cv_volume[icv_nbr], 2.0/3.0)/dist2;
        geomW = dist2;///pow(cv_volume[icv_nbr], 2.0/3.0);
#endif
        //double weight = 1.0/(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], 1.0e-8));
        double weight = 1.0/(geomW*(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], SDWLSClippingScalar)));

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
          assert(cvofa[ifa][1] == -1);

          // if p_bfa != NULL, then consider boundary values in calc gradient
          double dp = 0.0;
          if (p_fa != NULL)   dp = p_fa[ifa] - p[icv];

          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_fa[ifa][i] - x_cv[icv][i];

          double geomW = 1.0;
#ifdef geomWeighting
          double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
          //geomW = pow(cv_volume[icv], 2.0/3.0)/dist2;
          geomW = dist2;///pow(cv_volume[icv], 2.0/3.0);
#endif
          //double weight = 1.0/(limiterWeight(dp) + geomW + max(epsilonSDWLS*refValue[icv], 1.0e-8));
          double weight = 1.0/(geomW*(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], SDWLSClippingScalar)));

          swdx2  += weight * dx[0] * dx[0];
          swdy2  += weight * dx[1] * dx[1];
          swdz2  += weight * dx[2] * dx[2];
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
        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        double geomW = 1.0;
#ifdef geomWeighting
        double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        //geomW = pow(cv_volume[icv_nbr], 2.0/3.0)/dist2;
        geomW = dist2;///pow(cv_volume[icv_nbr], 2.0/3.0);
#endif
        //double weight = 1.0/(limiterWeight(dp) + geomW + max(epsilonSDWLS*refValue[icv], 1.0e-8));
        double weight = 1.0/(geomW*(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], SDWLSClippingScalar)));

        
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
          assert(cvofa[ifa][1] == -1);

          double dp = 0.0;
          dp = p_fa[ifa] - p[icv];
          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_fa[ifa][i] - x_cv[icv][i];

          double geomW = 1.0;
#ifdef geomWeighting
          double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
          //geomW = pow(cv_volume[icv], 2.0/3.0)/dist2;
          geomW = dist2;///pow(cv_volume[icv], 2.0/3.0);
#endif
          //double weight = 1.0/(limiterWeight(dp) + geomW + max(epsilonSDWLS*refValue[icv], 1.0e-8));
          double weight = 1.0/(geomW*(limiterWeight(dp) + max(epsilonSDWLS*refValue[icv], SDWLSClippingScalar)));

        
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

        double dx[3];
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        double geomW = 1.0;
#ifdef geomWeighting
        double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        //geomW = pow(cv_volume[icv_nbr], 2.0/3.0)/dist2;
        geomW = dist2;///pow(cv_volume[icv_nbr], 2.0/3.0);
#endif
        //double weight = 1.0/(limiterWeight(sqrt(du2)) + geomW + epsilonSDWLS*refValue[icv]);
        double weight = 1.0/(geomW*(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]));



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
          assert(cvofa[ifa][1] == -1);
          // use tiny weight here, du2 == 0

          double du2 = 0.0;
          if (u_fa != NULL)
          for (int i = 0; i < 3; i++)
          {
            const double du = u_fa[ifa][i] - u[icv][i];
            du2 += du*du;
          }
          double dx[3];
          for (int i = 0; i < 3; i++)
            dx[i] = x_fa[ifa][i] - x_cv[icv][i];
      
          double geomW = 1.0;
#ifdef geomWeighting
          double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
          //geomW = pow(cv_volume[icv], 2.0/3.0)/dist2;
          geomW = dist2;///pow(cv_volume[icv], 2.0/3.0);
#endif
          //double weight = 1.0/(limiterWeight(sqrt(du2)) + geomW + epsilonSDWLS*refValue[icv]);
          double weight = 1.0/(geomW*(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]));
        
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
        double dx[3], tmp;
        for (int i = 0; i < 3; i++)
          dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

        double geomW = 1.0;
#ifdef geomWeighting
        double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        //geomW = pow(cv_volume[icv_nbr], 2.0/3.0)/dist2;
        geomW = dist2;///pow(cv_volume[icv_nbr], 2.0/3.0);
#endif
        //double weight = 1.0/(limiterWeight(sqrt(du2)) + geomW + epsilonSDWLS*refValue[icv]);
        double weight = 1.0/(geomW*(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]));
        
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
          assert(cvofa[ifa][1] == -1);
          
          double du[3];
          double du2 = 0.0;
          for (int i = 0; i < 3; i++) 
          {
            du[i] = u_fa[ifa][i] - u[icv][i];
            du2 += du[i]*du[i];
          }

          double dx[3], tmp;
          for (int i = 0; i < 3; i++)
            dx[i] = x_fa[ifa][i] - x_cv[icv][i];

          double geomW = 1.0;
#ifdef geomWeighting
          double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
          //geomW = pow(cv_volume[icv], 2.0/3.0)/dist2;
          geomW = dist2;///pow(cv_volume[icv], 2.0/3.0);
#endif
          //double weight = 1.0/(limiterWeight(sqrt(du2)) + geomW + epsilonSDWLS*refValue[icv]);
          double weight = 1.0/(geomW*(limiterWeight(sqrt(du2)) + epsilonSDWLS*refValue[icv]));
          
        
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


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  



  int solveCvScalarJacobi(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double relax,
      const double zero, const int maxiter);

  int solveCvScalarCg(double * phi, double * Ap, double * rhs, const int mode, const double zero, const int maxiter);

  int solveCvScalarBcgstab(double * phi, double * Ap, double * rhs,
      const double zero, const double zeroRel, const int maxiter, char *scalName);
	
	int solveCvScalarBcgstabLine(double * phi, double * Ap, double * rhs,
													 const double zero, const double zeroRel, const int maxiter, char *scalName);

  int solveCvScalarBcgstab(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double zero,
      const int maxiter);
	
	int solveCvScalarLusgs(double * phi, double * Ap, double * rhs,
												 const double zero, const double zeroRel, const int maxiter, char *scalName);

  void matTimesVecOverCVs(double(*res)[5], double(*Ap)[5][5], double(*phi)[5]);
  
  void UpdateCvDataStateVec(double(*phi)[5]);
	
  void UpdateCvDataStateVecScal(double **phi, int nScal);

  int solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], 
      const double zeroAbs, const double zeroRel, const int maxiter, char *scalarName);
	
	int solveCvVectorR5BcgstabLine(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], 
														 const double zeroAbs, const double zeroRel, const int maxiter, char *scalarName);

	int solveCvVectorR5Lusgs(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
													 const double zeroRel, const int maxiter, char *scalarName);
	
	void LowerProductVectorR5(double(*LowProd)[5], double(*Ap)[5][5], double(*Phi)[5]);
	
	void UpperProductVectorR5(double(*UpperProd)[5], double(*Ap)[5][5], double(*Phi)[5]);
	
	void DiagonalProductVectorR5(double(*DiagProd)[5], double(*Ap)[5][5], double(*Phi)[5]);
	
	void LowerProductScalar(double *LowProd, double *Ap, double *Phi);
	
	void UpperProductScalar(double *UpperProd, double *Ap, double *Phi);
	
	void DiagonalProductScalar(double *DiagProd, double *Ap, double *Phi);
	
	
	void getMultMatrix(double c[5][5], double a[5][5], double b[5][5]);
	
	void getMultMatrix(double c[5][5], double a[5][5], double **b);
	
	void getMultMatrix(double c[5], double a[5][5], double b[5]);
	
	void getMultMatrix(double c[5], double b[5], double a[5][5]);
	
	void getSubsMatrix(double c[5][5], double a[5][5], double b[5][5]);
	
	void getSubsMatrix(double c[5][5], double a[5][5], double **b);
	
	void getSubsMatrix(double c[5], double a[5], double b[5]);
	
	void matTimesVecOverCVs_mg(double(*res)[5], double(*Ap)[5][5], double(*phi)[5], int iMesh);

	void UpdateCvDataStateVec_mg(double(*phi)[5], int iMesh);

	int solveCvVectorR5Bcgstab_mg(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], 
																const double zeroAbs, const double zeroRel, const int maxiter, char *scalarName, int iMesh);
	
	void LowerProduct_mg(double(*LowProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh);
	
	void UpperProduct_mg(double(*UpperProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh);
	
	void DiagonalProduct_mg(double(*DiagProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh);
	
	int solveCvVectorR5Lusgs_mg(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
															const double zeroRel, const int maxiter, char *scalarName, int iMesh);

	void updateCvData_mg(int * d, const int action, int iMesh) {
		
		if (iMesh == 1) {
			for (int icv_mgLevel1 = 0; icv_mgLevel1 < ncv_mgLevel1; icv_mgLevel1++)
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
					int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
					d_interpI1[icv_fine] = d[icv_mgLevel1];
				}	
			
			updateI1(d_interpI1, action, cvPrcommList);
			
			for (int icv_mgLevel1 = ncv_mgLevel1; icv_mgLevel1 < ncv_g_mgLevel1; icv_mgLevel1++) {
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel1] = d_interpI1[icv_fine];
			}
		}
		
		if (iMesh == 2) {
			for (int icv_mgLevel2 = 0; icv_mgLevel2 < ncv_mgLevel2; icv_mgLevel2++)
				for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
					int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
					for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1< nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
						int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
						d_interpI1[icv_fine] = d[icv_mgLevel2];
					}
				}
			
			updateI1(d_interpI1, action, cvPrcommList);
			
			for (int icv_mgLevel2 = ncv_mgLevel2; icv_mgLevel2 < ncv_g_mgLevel2; icv_mgLevel2++) {
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel2] = d_interpI1[icv_fine];
			}
		}
		
		if (iMesh == 3) {
			for (int icv_mgLevel3 = 0; icv_mgLevel3 < ncv_mgLevel3; icv_mgLevel3++)
				for (int iChildren_mgLevel3 = 0; iChildren_mgLevel3 < nChildren_CV_mgLevel3[icv_mgLevel3]; iChildren_mgLevel3++) {
					int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][iChildren_mgLevel3];
					for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
						int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
						for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1 < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
							int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
							d_interpI1[icv_fine] = d[icv_mgLevel3];
						}
					}
				}
			
			updateI1(d_interpI1, action, cvPrcommList);
			
			for (int icv_mgLevel3 = ncv_mgLevel3; icv_mgLevel3 < ncv_g_mgLevel3; icv_mgLevel3++) {
				int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][0];
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel3] = d_interpI1[icv_fine];
			}
		}
	}
	
	void updateCvData_mg(double *d, const int action, int iMesh) {

		if (iMesh == 1) {
			for (int icv_mgLevel1 = 0; icv_mgLevel1 < ncv_mgLevel1; icv_mgLevel1++)
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
					int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
					d_interpR1[icv_fine] = d[icv_mgLevel1];
				}	
			
			updateR1(d_interpR1, action, cvPrcommList);
			
			for (int icv_mgLevel1 = ncv_mgLevel1; icv_mgLevel1 < ncv_g_mgLevel1; icv_mgLevel1++) {
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel1] = d_interpR1[icv_fine];
			}
		}
		
		if (iMesh == 2) {
			for (int icv_mgLevel2 = 0; icv_mgLevel2 < ncv_mgLevel2; icv_mgLevel2++)
				for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
					int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
					for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1< nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
						int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
						d_interpR1[icv_fine] = d[icv_mgLevel2];
					}
				}
			
			updateR1(d_interpR1, action, cvPrcommList);
			
			for (int icv_mgLevel2 = ncv_mgLevel2; icv_mgLevel2 < ncv_g_mgLevel2; icv_mgLevel2++) {
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel2] = d_interpR1[icv_fine];
			}
		}
		
		if (iMesh == 3) {
			for (int icv_mgLevel3 = 0; icv_mgLevel3 < ncv_mgLevel3; icv_mgLevel3++)
				for (int iChildren_mgLevel3 = 0; iChildren_mgLevel3 < nChildren_CV_mgLevel3[icv_mgLevel3]; iChildren_mgLevel3++) {
					int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][iChildren_mgLevel3];
					for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
						int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
						for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1 < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
							int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
							d_interpR1[icv_fine] = d[icv_mgLevel3];
						}
					}
				}
			
			updateR1(d_interpR1, action, cvPrcommList);
			
			for (int icv_mgLevel3 = ncv_mgLevel3; icv_mgLevel3 < ncv_g_mgLevel3; icv_mgLevel3++) {
				int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][0];
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				d[icv_mgLevel3] = d_interpR1[icv_fine];
			}
		}
		
	}
	
	void updateCvData_mg(double(*d)[3], const int action, int iMesh) {
		
		if (iMesh == 1) {
			for (int icv_mgLevel1 = 0; icv_mgLevel1 < ncv_mgLevel1; icv_mgLevel1++)
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
					int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
					for (int iVar = 0; iVar < 3; iVar++)
						d_interpR2[icv_fine][iVar] = d[icv_mgLevel1][iVar];
				}	
			
			updateR2(d_interpR2, action, cvPrcommList);
			
			for (int icv_mgLevel1 = ncv_mgLevel1; icv_mgLevel1 < ncv_g_mgLevel1; icv_mgLevel1++) {
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					d[icv_mgLevel1][iVar] = d_interpR2[icv_fine][iVar];
			}
		}
		
		if (iMesh == 2) {
			for (int icv_mgLevel2 = 0; icv_mgLevel2 < ncv_mgLevel2; icv_mgLevel2++)
				for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
					int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
					for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1< nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
						int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
						for (int iVar = 0; iVar < 3; iVar++)
							d_interpR2[icv_fine][iVar] = d[icv_mgLevel2][iVar];
					}
				}
			
			updateR2(d_interpR2, action, cvPrcommList);
			
			for (int icv_mgLevel2 = ncv_mgLevel2; icv_mgLevel2 < ncv_g_mgLevel2; icv_mgLevel2++) {
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					d[icv_mgLevel2][iVar] = d_interpR2[icv_fine][iVar];
			}
		}
		
		if (iMesh == 3) {
			for (int icv_mgLevel3 = 0; icv_mgLevel3 < ncv_mgLevel3; icv_mgLevel3++)
				for (int iChildren_mgLevel3 = 0; iChildren_mgLevel3 < nChildren_CV_mgLevel3[icv_mgLevel3]; iChildren_mgLevel3++) {
					int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][iChildren_mgLevel3];
					for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
						int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
						for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1 < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
							int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
							for (int iVar = 0; iVar < 3; iVar++)
								d_interpR2[icv_fine][iVar] = d[icv_mgLevel3][iVar];
						}
					}
				}
			
			updateR2(d_interpR2, action, cvPrcommList);
			
			for (int icv_mgLevel3 = ncv_mgLevel3; icv_mgLevel3 < ncv_g_mgLevel3; icv_mgLevel3++) {
				int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][0];
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					d[icv_mgLevel3][iVar] = d_interpR2[icv_fine][iVar];
			}
		}
	}
	
	void updateCvData_mg(double(*d)[3][3], const int action, int iMesh) {
		if (iMesh == 1) {
			for (int icv_mgLevel1 = 0; icv_mgLevel1 < ncv_mgLevel1; icv_mgLevel1++)
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
					int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
					for (int iVar = 0; iVar < 3; iVar++)
						for (int jVar = 0; jVar < 3; jVar++)
							d_interpR3[icv_fine][iVar][jVar] = d[icv_mgLevel1][iVar][jVar] ;
				}	
			
			updateR3(d_interpR3, action, cvPrcommList);
			
			for (int icv_mgLevel1 = ncv_mgLevel1; icv_mgLevel1 < ncv_g_mgLevel1; icv_mgLevel1++) {
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					for (int jVar = 0; jVar < 3; jVar++)
						d[icv_mgLevel1][iVar][jVar]  = d_interpR3[icv_fine][iVar][jVar] ;
			}
		}
		
		if (iMesh == 2) {
			for (int icv_mgLevel2 = 0; icv_mgLevel2 < ncv_mgLevel2; icv_mgLevel2++)
				for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
					int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
					for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1< nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
						int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
						for (int iVar = 0; iVar < 3; iVar++)
							for (int jVar = 0; jVar < 3; jVar++)
								d_interpR3[icv_fine][iVar][jVar]  = d[icv_mgLevel2][iVar][jVar] ;
					}
				}
			
			updateR3(d_interpR3, action, cvPrcommList);
			
			for (int icv_mgLevel2 = ncv_mgLevel2; icv_mgLevel2 < ncv_g_mgLevel2; icv_mgLevel2++) {
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					for (int jVar = 0; jVar < 3; jVar++)
						d[icv_mgLevel2][iVar][jVar]  = d_interpR3[icv_fine][iVar][jVar] ;
			}
		}
		
		if (iMesh == 3) {
			for (int icv_mgLevel3 = 0; icv_mgLevel3 < ncv_mgLevel3; icv_mgLevel3++)
				for (int iChildren_mgLevel3 = 0; iChildren_mgLevel3 < nChildren_CV_mgLevel3[icv_mgLevel3]; iChildren_mgLevel3++) {
					int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][iChildren_mgLevel3];
					for (int iChildren_mgLevel2 = 0; iChildren_mgLevel2 < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren_mgLevel2++) {
						int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren_mgLevel2];
						for (int iChildren_mgLevel1 = 0; iChildren_mgLevel1 < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren_mgLevel1++) {
							int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren_mgLevel1];
							for (int iVar = 0; iVar < 3; iVar++)
								for (int jVar = 0; jVar < 3; jVar++)
									d_interpR3[icv_fine][iVar][jVar]  = d[icv_mgLevel3][iVar][jVar] ;
						}
					}
				}
			
			updateR3(d_interpR3, action, cvPrcommList);
			
			for (int icv_mgLevel3 = ncv_mgLevel3; icv_mgLevel3 < ncv_g_mgLevel3; icv_mgLevel3++) {
				int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][0];
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
				int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
				for (int iVar = 0; iVar < 3; iVar++)
					for (int jVar = 0; jVar < 3; jVar++)
						d[icv_mgLevel3][iVar][jVar]  = d_interpR3[icv_fine][iVar][jVar] ;
			}
		}
	}
		
protected:

  // ---------------------------
  // ghost/fake data updates
  // ---------------------------

  void updateCvData(int * d, const int action) {
    updateI1(d, action, cvPrcommList);
  }

  void updateCvData(double * d, const int action) {
    updateR1(d, action, cvPrcommList);
  }

  void updateCvData(double(*d)[3], const int action) {
    updateR2(d, action, cvPrcommList);
  }

  void updateCvData(double(*d)[3][3], const int action) {
    updateR3(d, action, cvPrcommList);
  }

  void updateCvSymmetricR3(double(*diag)[3], double(*offd)[3], const int action) {
    updateSymmetricR3(diag, offd, action, cvPrcommList);
  }

  void updateCvSymmetricR3(double(*d)[6], const int action) {
    // here the 6 tensor elements are stored as
    // xx,xy,xz,yy,yz,zz
    updateSymmetricR3(d, action, cvPrcommList);
  }


  void updateCvDataByName(char *name, const int action) {
    double *scal = getScalar(name);
    updateCvData(scal, action);
  }
};

#endif

