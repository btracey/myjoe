#include "UgpWithCvFakeOp.h"
#include "Param.h"
#include "CdpFilter.h"
#include "MshFilter.h"
#include "MiscUtils.h"
#include <math.h>

//#include "shapeFunctions.h"

// cf4 - my 4th attempt at a cv-based compressible flow solver for LES
// F. Ham - April 2009

enum SgsModels {
  SMAGORINSKY_SGS,
  CSM_SGS,
  VREMAN_SGS,
  NO_SGS
};      

class Charles : public UgpWithCvFakeOp, public ParamMap {

private:

  double wtime;
  int ifa_debug;

  int runtime_flag;
  double runtime;

protected:
  
  double time; /// simulation time
  double dt;   /// time step
  int step;    /// step index for simulation
  int nsteps;  /// number of steps to run
  double dt_target;

  // shock-capturing parameters... 
  double dilCoeff;
  double swCoeff;

  // gas props...
  double gamma,p_ref,rho_ref,T_ref,R_gas,mu_ref,mu_power_law,Pr;
  
  // cv field data...
  double *rho;
  double (*rhou)[3];
  double *rhoE;
  
  // for visualizion...
  double *dil;

  int sgs_model;
  double *mu_sgs;
  double *k_sgs;

  // primitive data...
  double (*u)[3];
  double *p;
  double *T;
  double *h;
  
  // how often to look closely at output...
  int check_interval;

  // operator mods for shock-capturing...
  double * fa_ss;
  double * fa_ss_no;
  double * fa_alpha_ss;

  // boundary conditions...
  /*
  double * rho_bc;
  double (*rhou_bc)[3];
  double * rhoE_bc;
  */

public:
  
  // contructor for this class. Needs to be called
  // with a parameter filename to initialize the 
  // parameters (stored in a ParamMap)...

  Charles(char *name) : ParamMap(name) {

    if (mpi_rank == 0)
      cout << "Charles()"<< endl;

    // for timing...
    
    wtime = MPI_Wtime();
    
    // -------------------------
    // gas props...
    // -------------------------

    gamma   = getDoubleParam("GAMMA",1.4);
    p_ref   = getDoubleParam("P_REF");
    rho_ref = getDoubleParam("RHO_REF");
    T_ref   = getDoubleParam("T_REF");

    // ideal gas law is P / rho = R * T ...
    R_gas = p_ref/rho_ref/T_ref;
    
    // viscosity at Tref...
    mu_ref = getDoubleParam("MU_REF");
    
    // note: specify 0 here for constant visc = mu_ref
    mu_power_law = getDoubleParam("MU_POWER_LAW",0.76);
    Pr = getDoubleParam("PR",0.7);

    // shock-capturing scheme will be applied for local
    // absolute dilitations (i.e. both compressions AND expansions) 
    // above this value. For reason of rocustness, its default value is
    // set to slightly negative, so it is applied everywhere.
    dilCoeff = getDoubleParam("DIL_COEFF",-1.0);

    // solution-weighted coeff: a larger value here reduces the amount of solution-dependent
    // weighting...
    swCoeff = getDoubleParam("SW_COEFF",0.1); 

    // --------------------------
    // data registration. Registers our
    // class's (Charles class) data with the
    // Ugp, so we can read/write/visualize
    // our field data by name. Also handles
    // memory allocation.
    // --------------------------

    registerValue(dt,  "DT");
    registerValue(step,"STEP");
    registerValue(time,"TIME");

    // note these may be overwritten when restart file is read in...

    dt = 0.0; // gets set in calcDt()
    step = 0;
    time = 0.0;

    dt_target = getDoubleParam("DT");
    //upwindFactor = getDoubleParam("UPWIND_FACTOR",0.25);
    nsteps = getIntParam("NSTEPS",-1);

    runtime_flag = 0;
    if (checkParam("RUNTIME")) {
      runtime_flag = 1;
      if (mpi_rank == 0) {
	runtime = getDoubleParam("RUNTIME");
	cout << " > RUNTIME [hours]: " << runtime << endl;
	// convert to the time in seconds when we have to stop...
	runtime = runtime*3600.0 + wtime;
      }
    }
    
    // register data...
    // does the memory allocation as well.
    
    rho  = NULL; registerScalar(rho, "RHO", CV_DATA);
    rhou = NULL; registerVector(rhou,"RHOU",CV_DATA);
    rhoE = NULL; registerScalar(rhoE,"RHOE",CV_DATA);

    u = NULL; registerVector(u,"U",CV_DATA);
    p = NULL; registerScalar(p,"P",CV_DATA);
    T = NULL; registerScalar(T,"T",CV_DATA);

    /*
    // Interp hack...
    rho  = NULL; registerScalar(rho, "RHO_CV", CV_DATA);
    rhou = NULL; registerVector(rhou,"RHOU_CV",CV_DATA);
    rhoE = NULL; registerScalar(rhoE,"RHOE_CV",CV_DATA);
    */

    dil    = NULL; registerScalar(dil,    "DIL",    CV_DATA);
    mu_sgs = NULL; registerScalar(mu_sgs, "MU_SGS", CV_DATA);
    k_sgs  = NULL; registerScalar(k_sgs,  "K_SGS",  CV_DATA);

    fa_ss_no = NULL; registerScalar(fa_ss_no,"FA_SS",NO_DATA);

    fa_ss = NULL;
    fa_alpha_ss = NULL;

    // bcs...
    
    /*
    rho_bc = NULL;
    rhou_bc = NULL;
    rhoE_bc = NULL;
    */

    // sgs model...

    Param *p;
    if (getParam(p,"SGS_MODEL")) {
      string name = p->getString();
      if (name == "SMAGORINSKY") {
	if (mpi_rank == 0)
	  cout << "SGS_MODEL: SMAGORINSKY" << endl;
	sgs_model = SMAGORINSKY_SGS;
      }
      else if (name == "CSM") {
	if (mpi_rank == 0)
	  cout << "SGS_MODEL: CSM (Coherent Structure Model)" << endl;
	sgs_model = CSM_SGS;
      }
      else if (name == "VREMAN") {
	if (mpi_rank == 0)
	  cout << "SGS_MODEL: VREMAN" << endl;
	sgs_model = VREMAN_SGS;
      }
      else if (name == "NONE") {
	if (mpi_rank == 0)
	  cout << "SGS_MODEL: NONE" << endl;
	sgs_model = NO_SGS;
      }
      else {      
	if (mpi_rank == 0)
	  cerr << "Error: unrecognized SGS_MODEL: " << name << endl;
	throw(-1);
      }
    }
    else {
      if (mpi_rank == 0)
	cout << "SGS_MODEL: NONE (default)" << endl;
      sgs_model = NO_SGS;
    }

    // stats - after all registration...
    
    if (getParam(p,"STATS"))
      initStats(p);

    // --------------------------
    // report...
    // --------------------------
    
    if (mpi_rank == 0) {
      cout << "Gas properties: " << endl;
      cout << "  gamma        : " << gamma << endl;
      cout << "  p_ref        : " << p_ref << endl;
      cout << "  rho_ref      : " << rho_ref << endl;
      cout << "  T_ref        : " << T_ref << endl;
      cout << "  R_gas        : " << R_gas << endl;
      cout << "  mu_ref       : " << mu_ref << endl;
      cout << "  mu_power_law : " << mu_power_law << endl;
      cout << "  Pr           : " << Pr << endl;
      cout << "Solver settings: " << endl;
      cout << "  dt           : " << dt << endl;
      cout << "  nsteps       : " << nsteps << endl;
      cout << "  dilCoeff     : " << dilCoeff << endl;
      cout << "  swCoeff      : " << swCoeff << endl;
    }

    // other stuff...
    
    check_interval = getIntParam("CHECK_INTERVAL",1);

  }

  virtual ~Charles() {
   
    delete[] u;
    delete[] p;
    delete[] T;
    delete[] h;
    
  }

  /*
  void transformMeshHook() {
    cout << "WARNING: **************** transforming mesh ********************" << endl;
    srand(1); // seed the random number generator so the grids are repeatable...    
    FOR_INO FOR_I3 x_no[ino][i] +=  0.01*(2.0*rand()/(double)RAND_MAX - 1.0);
  }
  */

  void init() {
    
    // first step is to read a restart file...
    // may or may not contain data. Multiple formats
    // are also supported...
    
    string restart = getStringParam("RESTART");
    string ext = getFileNameExtension(restart, ".");
    if (ext == "cdp") {
      CdpFilter cf(restart);
      cf.minInitUgp(this);
    } 
    else if ((ext == "msh")||(ext == "cas")) {
      MshFilter mf(restart);
      mf.initUgpMin(this);
    } 
    else {
      // assume this is a standard cdp file...
      readRestart(restart);
    }
    
    // initialize the unstructured grid partition...
    
    double fa_alpha_coeff = getDoubleParam("FA_ALPHA_COEFF",2.0);
    UgpWithCvFakeOp::init(fa_alpha_coeff);
    
    if (checkParam("FIRST_ORDER"))
      firstOrderOperators();
    
    if (checkParam("RESET_STATS"))
      resetStats();

    // allocate primitive data...
    
    // u p T can now be used for stats
    //u = new double[ncv_gf][3];
    //p = new double[ncv_gf];
    //T = new double[ncv_gf];
    h = new double[ncv_gf];

    // the shock-sensor lives at the faces...
    fa_ss = new double[nfa];

    // this is the fa_alpha value including the shock sensor info...
    fa_alpha_ss = new double[nfa];

    // let the user set an initial condition, and
    // modify operators if desired...
    
    initialHook();
    
    // if you wanr to register FA_ALPHA and visualize...
    // this is in UgpWithCvFakeOp

    initFaAlphaVis();

    // update ghost data, make sure everything 
    // is consistent between processors...

    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou,REPLACE_ROTATE_DATA);
    updateCvData(rhoE,REPLACE_DATA);
    
    /*
    rho_bc = new double[nfa_b];
    rhou_bc = new double[nfa_b][3];
    rhoE_bc = new double[nfa_b];
    */

    updateFakeBCData(rho,rhou,rhoE,time,1); // verbose

    // look for WRITE_DATA lines...

    vector<Param> *paramVec;
    if (getParam(paramVec,"WRITE_DATA"))
      initWriteData(paramVec);

    // and write the data if this is step 0...

    if (step == 0)
      writeData(step);

  }
  
protected:
  
  // the hook routines are normally made "virtual" to
  // allow the user to automatically override these
  // with their own hooks...
  
  virtual void initialHook() {
    
    // look for default initialization...
    // for now, force the user to specifiy either ALL of RHO0, U0 and P0,
    // or the initial field of step 0 will NOT be initialized...
    
    int init_flag = 0;
    double rho0,u0[3],p0;
    
    if (checkParam("RHO0")) {
      init_flag = 1;
      rho0 = getDoubleParam("RHO0");
    }
    
    Param *p;
    if (getParam(p,"U0")) {
      assert( init_flag == 1 );
      u0[0] = p->getDouble(1);
      u0[1] = p->getDouble(2);
      u0[2] = p->getDouble(3);
    }
    
    if (checkParam("P0")) {
      assert( init_flag == 1 );
      p0 = getDoubleParam("P0");
    }
    
    if (init_flag == 1) {
      if (mpi_rank == 0)
	cout << "setting initial condition: RHO0, U0, P0: " << rho0 << " " <<
	  u0[0] << " " << u0[1] << " " << u0[2] << " " << p0 << endl;
      double rhoE0 = p0/(gamma - 1.0) + 
	0.5*rho0*( u0[0]*u0[0] + 
		   u0[1]*u0[1] + 
		   u0[2]*u0[2] );
      FOR_ICV {
	rho[icv] = rho0;
	FOR_I3 rhou[icv][i] = rho0*u0[i];
	rhoE[icv] = rhoE0;
      }
    }
    
    
  }
  
  virtual void temporalHook() {
    
  }

  virtual void finalHook() {
    
  }

public:
  
  void run() {
    
    // standard explicit 3rd-order RK solver...

    if (mpi_rank == 0)
      cout << "run()" << endl;

    // data req'd by the RK solver - no ghosts req'd
    // here (just locally owned cvs), so no need to use 
    // registration. Just allocate them here...
    
    double *drho1       = new double[ncv];
    double (*drhou1)[3] = new double[ncv][3];
    double *drhoE1      = new double[ncv];
    
    double *drho2       = new double[ncv];
    double (*drhou2)[3] = new double[ncv][3];
    double *drhoE2      = new double[ncv];
    
    double *drho3       = new double[ncv];
    double (*drhou3)[3] = new double[ncv][3];
    double *drhoE3      = new double[ncv];
  
    double *rho0        = new double[ncv];
    double (*rhou0)[3]  = new double[ncv][3];
    double *rhoE0       = new double[ncv];

    int done = doneSolver();
    while (done != 1) {
      
      step++;
      calcDt();
      time += dt;
      
      if ((check_interval > 0)&&(step%check_interval == 0)&&(mpi_rank == 0)) {
	cout << 
	  "\n----------------------------------------------------------\n" << 
	  " starting step: " << step << " time: " << time << " dt: " << dt << 
	  "\n----------------------------------------------------------" << endl;
      }
      
      calcSgsStuff(step%check_interval == 0); // verbosity
      
      // -------------------------------------------------
      // copy the current solution into rho0,rhou0, etc...
      // -------------------------------------------------
      
      FOR_ICV {
	rho0[icv]            = rho[icv];
	FOR_I3 rhou0[icv][i] = rhou[icv][i];
	rhoE0[icv]           = rhoE[icv];
      }
    
      // ------------
      // rk step 1...
      // ------------
      // return a rhs (integrated over volume) at 
      // the local cv's...
      calcRhs(drho1,drhou1,drhoE1,rho,rhou,rhoE,time-dt,1); // also include the step number
      // invert...
      FOR_ICV {
	double tmp = dt/cv_volume[icv];
	drho1[icv]            *= tmp;
	FOR_I3 drhou1[icv][i] *= tmp;
	drhoE1[icv]           *= tmp;
	// ----------------------------
	// and do the update...
	// ----------------------------
	rho[icv]              += 0.5*drho1[icv];
	FOR_I3 rhou[icv][i]   += 0.5*drhou1[icv][i];
	rhoE[icv]             += 0.5*drhoE1[icv];
      }

      // update Ghosts...
      updateCvData(rho, REPLACE_DATA);
      updateCvData(rhou,REPLACE_ROTATE_DATA);
      updateCvData(rhoE,REPLACE_DATA);
      
      updateFakeBCData(rho,rhou,rhoE,time-0.5*dt);
      
      // ------------
      // rk step 2...
      // ------------
      calcRhs(drho2,drhou2,drhoE2,rho,rhou,rhoE,time-0.5*dt,2);
      FOR_ICV {
	double tmp = dt/cv_volume[icv];
	drho2[icv]            *= tmp;
	FOR_I3 drhou2[icv][i] *= tmp;
	drhoE2[icv]           *= tmp;
	// ----------------------------
	// and do the update...
	// ----------------------------
	rho[icv]            = rho0[icv]     - drho1[icv]     + 2.0*drho2[icv];
	FOR_I3 rhou[icv][i] = rhou0[icv][i] - drhou1[icv][i] + 2.0*drhou2[icv][i];
	rhoE[icv]           = rhoE0[icv]    - drhoE1[icv]    + 2.0*drhoE2[icv];
      }
      updateCvData(rho, REPLACE_DATA);
      updateCvData(rhou,REPLACE_ROTATE_DATA);
      updateCvData(rhoE,REPLACE_DATA);

      updateFakeBCData(rho,rhou,rhoE,time);
      
      // ------------
      // rk step 3...
      // ------------
      calcRhs(drho3,drhou3,drhoE3,rho,rhou,rhoE,time,3);
      FOR_ICV {
	double tmp = dt/cv_volume[icv];
	drho3[icv]            *= tmp;
	drhoE3[icv]           *= tmp;
	FOR_I3 drhou3[icv][i] *= tmp;
	// ----------------------------
	// and do the update...
	// ----------------------------
	rho[icv]            = rho0[icv]     + (drho1[icv]     + 4.0*drho2[icv]     + drho3[icv])/6.0; 
	FOR_I3 rhou[icv][i] = rhou0[icv][i] + (drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i])/6.0; 
	rhoE[icv]           = rhoE0[icv]    + (drhoE1[icv]    + 4.0*drhoE2[icv]    + drhoE3[icv])/6.0;
      }
      updateCvData(rho, REPLACE_DATA);
      updateCvData(rhou,REPLACE_ROTATE_DATA);
      updateCvData(rhoE,REPLACE_DATA);

      updateFakeBCData(rho,rhou,rhoE,time);
      
      // ------------
      // take a look...
      // ------------ 
      if (step%check_interval == 0) {
	
	dumpScalarRange(rho,ncv,"RHO");
	dumpVectorRange(rhou,ncv,"RHOU");
	dumpScalarRange(rhoE,ncv,"RHOE");

	calcDilitation(dil);
	dumpScalarRange(dil,ncv,"DIL");
	
	// compute current rho,p,T for stats...
	FOR_ICV {
	  FOR_I3 u[icv][i] = rhou[icv][i]/rho[icv];
	  p[icv] = (gamma-1.0)*( rhoE[icv] - 0.5*( rhou[icv][0]*u[icv][0] +
						   rhou[icv][1]*u[icv][1] + 
						   rhou[icv][2]*u[icv][2] ) );
	  T[icv] = p[icv]/rho[icv]/R_gas;
	}
	
	// statistics...
	updateStats((double)check_interval*dt);

	if (mpi_rank == 0) {
	  double wtime0 = wtime;
	  wtime = MPI_Wtime();
	  cout << " > time since last check[s]: " << wtime - wtime0 << endl;
	}

	done = doneSolver();
	
      }

      temporalHook();
      
      writeData(step);
      
    }
    
    finalHook();

    writeRestart();
    
  }
  
  void calcDilitation(double * dil) {
    
    double * undA = new double[nfa];
    FOR_IFA undA[ifa] = 0.0;
    
    // left states...
    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      const int cof_f = cvofa0_i[ifa];
      const int cof_l = cvofa0_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	double coeff = 0.5*cvofa0_value[cof];
	FOR_I3 undA[ifa] += coeff*rhou[icv][i]/rho[icv]*fa_normal[ifa][i];
      }
    }
    
    // here we subtract across processors because 
    //normals are in opposite directions...
    updateFaR1(undA,SUBTRACT_DATA);

    // right states...
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
      const int cof_f = cvofa1_i[ifa];
      const int cof_l = cvofa1_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];
	double coeff = 0.5*cvofa1_value[cof];
	FOR_I3 undA[ifa] += coeff*rhou[icv][i]/rho[icv]*fa_normal[ifa][i];
      }
    }

    // for boundaries, just use left reconstruction...
    for (int ifa = 0; ifa < nfa_b; ++ifa) {
      const int cof_f = cvofa0_i[ifa];
      const int cof_l = cvofa0_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	double coeff = cvofa0_value[cof];
	FOR_I3 undA[ifa] += coeff*rhou[icv][i]/rho[icv]*fa_normal[ifa][i];
      }
    }

    FOR_ICV dil[icv] = 0.0;
    for (int ifa = 0; ifa < nfa_bpi; ++ifa) {
      // these faces have only icv0 local...
      const int icv0 = cvofa[ifa][0];
      dil[icv0] += undA[ifa];
    }
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
      // these faces have both icv0 and icv1 local...
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      dil[icv0] += undA[ifa];
      dil[icv1] -= undA[ifa];
    }
    FOR_ICV dil[icv] /= cv_volume[icv];

    delete[] undA;

  }
  
  void updateFakeBCData(double * rho,double (*rhou)[3],double * rhoE,
			const double rk_time,const int verbose = 0) {
    
    int bc_err = 0;
    
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
      if (zone->getKind() == FA_ZONE_BOUNDARY) {
	
	Param * param;
	if ( getParam(param,zone->getName()) ) {
	  
	  if ( param->getString() == "WALL_ISOTHERMAL" ) {

	    double Twall = param->getDouble(2);
	    
	    if ((mpi_rank == 0)&&(verbose))
	      cout << "Applying WALL_ISOTHERMAL bc to zone: " << zone->getName() << ", T: " << Twall << endl;
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      const int icv0 = cvofa[ifa][0];
	      assert( (icv0 >= 0)&&(icv0 < ncv) ); // icv0 is local
	      
	      const int icv1 = cvofa[ifa][1];
	      assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) ); // icv1 is "fake"
	      
	      // here we want the velocity to pass through 0, and the temperature to
	      // pass through Twall...
	      
	      double p0 = (gamma-1.0)*( rhoE[icv0] - 0.5*( rhou[icv0][0]*rhou[icv0][0] +
							   rhou[icv0][1]*rhou[icv0][1] +
							   rhou[icv0][2]*rhou[icv0][2] )/rho[icv0] );

	      // T0 = p0/rho0/R_gas...
	      double T0 = p0/rho[icv0]/R_gas;
	      
	      // use p1 == p0, 0.5*(T0 + T1) == Twall...
	      
	      double T1 = 2.0*Twall - T0;
	      rho[icv1] = p0/T1/R_gas;

	      // uwall == 0.0...
	      
	      rhou[icv1][0] = -rho[icv1]*rhou[icv0][0]/rho[icv0];
	      rhou[icv1][1] = -rho[icv1]*rhou[icv0][1]/rho[icv0];
	      rhou[icv1][2] = -rho[icv1]*rhou[icv0][2]/rho[icv0];
	      
	      // p0 == p1 
	      
	      rhoE[icv1] = p0/(gamma-1.0) + 0.5*( rhou[icv1][0]*rhou[icv1][0] +
						  rhou[icv1][1]*rhou[icv1][1] +
						  rhou[icv1][2]*rhou[icv1][2] )/rho[icv1];
	      
	    }

	  }
	  else if ( param->getString() == "WALL_ADIABATIC" ) {
	    
	    if ((mpi_rank == 0)&&(verbose))
	      cout << "Applying WALL_ADIABATIC bc to zone: " << zone->getName() << endl;
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      const int icv0 = cvofa[ifa][0];
	      assert( (icv0 >= 0)&&(icv0 < ncv) ); // icv0 is local
	      
	      const int icv1 = cvofa[ifa][1];
	      assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) ); // icv1 is "fake"
	      
	      rho[icv1]            = rho[icv0];
	      FOR_I3 rhou[icv1][i] = -rhou[icv0][i]; // reflect rhou
	      rhoE[icv1]           = rhoE[icv0];
	      
	    }
	    
	  }	      
	  else if (  param->getString() == "CBC" ) {
	    
	    double rho_bc = param->getDouble(2);
	    double u_bc[3];
	    u_bc[0] = param->getDouble(3);
	    u_bc[1] = param->getDouble(4);
	    u_bc[2] = param->getDouble(5);
	    double p_bc = param->getDouble(6);
	    
	    if ((mpi_rank == 0)&&(verbose))
	      cout << "Applying CBC (characteristic bc) to zone: " << zone->getName() << ", rho, u, p: " << 
		rho_bc << " " << u_bc[0] << " " << u_bc[1] << " " << u_bc[2] << " " << p_bc << endl;
	    
	    double rhoE_bc = p_bc/(gamma - 1.0) + 
	      0.5*rho_bc*( u_bc[0]*u_bc[0] + 
			   u_bc[1]*u_bc[1] + 
			   u_bc[2]*u_bc[2] );

	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      const int icv0 = cvofa[ifa][0];
	      assert( (icv0 >= 0)&&(icv0 < ncv) ); // icv0 is local
	      
	      const int icv1 = cvofa[ifa][1];
	      assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) ); // icv1 is "fake"
	      
	      rho[icv1]            = rho_bc;
	      FOR_I3 rhou[icv1][i] = rho_bc*u_bc[i];
	      rhoE[icv1]           = rhoE_bc;
	      
	    }
	    
	  }
	  else if (  param->getString() == "SYMMETRY" ) {
	    
	    if ((mpi_rank == 0)&&(verbose))
	      cout << "Applying SYMMETRY bc to zone: " << zone->getName() << endl;
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      const int icv0 = cvofa[ifa][0];
	      assert( (icv0 >= 0)&&(icv0 < ncv) ); // icv0 is local
	      
	      const int icv1 = cvofa[ifa][1];
	      assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) ); // icv1 is "fake"
	      
	      // get the outward-pointing unit normal...
	      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] + 
				  fa_normal[ifa][1]*fa_normal[ifa][1] + 
				  fa_normal[ifa][2]*fa_normal[ifa][2] );
	      double unitNormal[3] = { fa_normal[ifa][0]/area,
			      fa_normal[ifa][1]/area,
			      fa_normal[ifa][2]/area };
	      
	      // the right state is the same as the left state with normal components
	      // of the velocity reflected...
	      rho[icv1]            = rho[icv0];
	      double rhoun = rhou[icv0][0]*unitNormal[0] + rhou[icv0][1]*unitNormal[1] + rhou[icv0][2]*unitNormal[2];
	      FOR_I3 rhou[icv1][i] = rhou[icv0][i] - 2.0*rhoun*unitNormal[i];
	      rhoE[icv1]           = rhoE[icv0];
	      
	    }
	    
	  }
	  else if (  param->getString() == "HOOK" ) {

	    if ((mpi_rank == 0)&&(verbose))
	      cout << "Applying HOOK bc to zone: " << zone->getName() << endl;

	    // pass control to the user to populate the ghost cv's associated with
	    // this boundary zone. note that the ghost cv's are ordered the same as the
	    // boundary faces, except they are at the end of the cv list (boundary faces are
	    // at the start)...
	    if ( bcHook(&(rho[ncv_g]),&(rhou[ncv_g]),&(rhoE[ncv_g]),*zone,rk_time) == -1 ) {
	      if (mpi_rank == 0)
		cerr << "Error: bcHook not setup for zone: " <<  zone->getName() << endl;
	      bc_err = -1;
	    }
	    
	  }
	  else {
	    
	    if (mpi_rank == 0)
	      cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
	    bc_err = -1;
	    
	  }
	}
	else {
	  
	  if (mpi_rank == 0)
	    cerr << "Error: cannot determine bc for zone: " << zone->getName() << endl;
	  bc_err = -1;
	  
	}
	
      }
    }
    
    if (bc_err != 0) {
      if (mpi_rank == 0)
	cout << "Error: there was a bc problem" << endl;
      throw(-1);
    }
    
  }
  
  virtual int bcHook(double * rho_bc,double (*rhou_bc)[3],double * rhoE_bc,
		     FaZone& zone,const double rk_time) {
    
    // return 0 if you handle the bc, otherwise returning -1 will throw a bc error... 

    /*
    // here is a sample snippet...
    if (zone.getNameString() == "y1-FLUID") {
      for (int ifa = zone.ifa_f; ifa <= zone.ifa_l; ++ifa) {
	rho_bc[ifa] = rho_ref;
	rhou_bc[ifa][0] = rho_ref*3.0;
	rhou_bc[ifa][1] = 0.0;
	rhou_bc[ifa][2] = 0.0;
	rhoE_bc[ifa] = p_ref/(gamma - 1.0) + 
	  0.5*( rhou_bc[ifa][0]*rhou_bc[ifa][0] + 
		rhou_bc[ifa][1]*rhou_bc[ifa][1] + 
		rhou_bc[ifa][2]*rhou_bc[ifa][2] )/rho_bc[ifa];
      }
      return(0);
    }
    */

    // note that you can also access the face coordinates using x_fa[ifa][0..2] 
    
    return(-1);
    
  }

  void calcDt() {

    dt = dt_target;

  }

  int doneSolver() {

    // ---------------------------------------------
    // returns 1 if we are done, 0 otherwise...
    // ---------------------------------------------
    
    int done = 0;
  
    // if the current completed step is greater or equal to the user
    // specified nsteps, then we are done...
    
    if ((nsteps >= 0)&&(step >= nsteps))
      done = 1;
  
    // some stuff needs to be decided by rank 0 only...
    
    if (mpi_rank == 0) {
      
      // RUNTIME...
      
      if ((runtime_flag)&&(MPI_Wtime() > runtime)) {
	cout << " > reached RUNTIME" << endl;
	done = 1;
      }
      
      // killcf4 file...
      
      int ierr = MPI_File_delete("killcf4",MPI_INFO_NULL);
      if (ierr == 0) {
	cout << " > found file killcf4" << endl;
	done = 1;
      }
      
    }
    MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
    
    return(done);
    
  }

  void calcSgsStuff(const int verbose) {

    // define cv-based mu_sgs and k_sgs 
    // WARNING: note that the sgs quantities are NOT set in the fake cells!...
    
    switch (sgs_model) {
    case NO_SGS:
      
      if ((mpi_rank == 0)&&(verbose))
	cout << " > calcSgsStuff: zeroing mu and k" << endl;
      
      FOR_ICV_G mu_sgs[icv] = k_sgs[icv] = 0.0;
      return;
      
    case VREMAN_SGS:
      
      if ((mpi_rank == 0)&&(verbose))
	cout << " > calcSgsStuff: VREMAN" << endl;
      
      {
	
	// -------------------------------------------
	// VR_PI = sqrt(B/alal)
	// 
	// where,
	// 
	// B      = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23
	// alij   = duj_dxi
	// betaij = delta_1**2 * al1i*al1j + delta_2**2 * al2i*al2j + delta_3**2 * al3i*al3j
	// alal   = sum(sij_d(1:3,icv)**2) + 2.0_WP*sum(sij_od(1:3,icv)**2)
	// -------------------------------------------
	
	const double vreman_coeff= 0.07; 
	const double Pr_t = 0.9;
	
	/*
	// hack the velocity field...
	FOR_ICV_GF {
	  rho[icv] = 1.0;
	  rhou[icv][0] = 1.1 + 1.2*x_cv[icv][0] + 1.32*x_cv[icv][1] + 1.432*x_cv[icv][2];
	  rhou[icv][1] = 2.1 + 2.2*x_cv[icv][0] + 2.32*x_cv[icv][1] + 2.432*x_cv[icv][2];
	  rhou[icv][2] = 3.1 + 3.2*x_cv[icv][0] + 3.32*x_cv[icv][1] + 3.432*x_cv[icv][2];
	}
	*/
	
	double (*duidxj)[3][3] = new double[ncv][3][3];
	calcCvStrainRate(duidxj,rho,rhou);
	
	// if you are hacking the velocity field... 
	//dumpTensorRange(duidxj,ncv,"DUIDXJ");
	//throw(-1);
	
	FOR_ICV {
	  
	  double dx2 = pow( cv_volume[icv], 2.0/3.0 );
	  double dy2 = dx2;
	  double dz2 = dx2;
	  
	  double alpha11 = duidxj[icv][0][0];
	  double alpha22 = duidxj[icv][1][1];
	  double alpha33 = duidxj[icv][2][2];

	  double alpha12 = duidxj[icv][0][1];
	  double alpha13 = duidxj[icv][0][2];
	  double alpha23 = duidxj[icv][1][2];
	  
	  double alpha21 = duidxj[icv][1][0];
	  double alpha31 = duidxj[icv][2][0];
	  double alpha32 = duidxj[icv][2][1];
	  
	  double beta11  = dx2*alpha11*alpha11+dy2*alpha21*alpha21+dz2*alpha31*alpha31;
	  double beta12  = dx2*alpha11*alpha12+dy2*alpha21*alpha22+dz2*alpha31*alpha32;
	  double beta13  = dx2*alpha11*alpha13+dy2*alpha21*alpha23+dz2*alpha31*alpha33;
	  double beta22  = dx2*alpha12*alpha12+dy2*alpha22*alpha22+dz2*alpha32*alpha32;
	  double beta23  = dx2*alpha12*alpha13+dy2*alpha22*alpha23+dz2*alpha32*alpha33;
	  double beta33  = dx2*alpha13*alpha13+dy2*alpha23*alpha23+dz2*alpha33*alpha33;

	  double B       = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23;
	  B              = (B + fabs(B))*0.5;
	  double alal    = 
	    alpha11*alpha11+alpha22*alpha22+alpha33*alpha33 +
	    alpha12*alpha12+alpha13*alpha13+alpha23*alpha23 + 
	    alpha21*alpha21+alpha31*alpha31+alpha32*alpha32;
	  
	  double s_mag = sqrt(B/alal); // includes lengthscale squared too
	  
	  // mu_sgs...
	  mu_sgs[icv] = rho[icv]*vreman_coeff*s_mag;

	  // clip at the constant smagorinsky level with c2 = 0.5^2...
	  //mu_sgs[icv] = min( mu_sgs[icv], rho_no[icv]*0.05*dx2*sqrt(2.0*alal) );
	  
	}

	delete[] duidxj;
	
	// ------------------------------------------------
	// filter a few times using a cv->no->cv filter...
	// ------------------------------------------------

	// 1. pre-compute the normalization. Equal weighting, so use no_flag...

	FOR_INO no_flag[ino] = 0;
	
	FOR_ICV {
	  const int noc_f = noocv_i[icv];
	  const int noc_l = noocv_i[icv+1]-1;
	  for (int noc = noc_f; noc <= noc_l; ++noc) {
	    const int ino = noocv_v[noc];
	    no_flag[ino] += 1;
	  }
	}
	
	updateNoI1(no_flag,ADD_DATA);
	
	// 2. now do 3 cycles of cv->no->cv...
	
	double * mu_sgs_no = new double[nno];
	
	for (int iter = 0; iter < 3; ++iter) {
	  
	  FOR_INO mu_sgs_no[ino] = 0.0;
	  
	  FOR_ICV {
	    const int noc_f = noocv_i[icv];
	    const int noc_l = noocv_i[icv+1]-1;
	    for (int noc = noc_f; noc <= noc_l; ++noc) {
	      const int ino = noocv_v[noc];
	      mu_sgs_no[ino] += mu_sgs[icv];
	    }
	  }

	  updateNoR1(mu_sgs_no,ADD_DATA);
	  
	  FOR_INO mu_sgs_no[ino] /= (double)no_flag[ino];
	  
	  FOR_ICV {
	    mu_sgs[icv] = 0.0;
	    const int noc_f = noocv_i[icv];
	    const int noc_l = noocv_i[icv+1]-1;
	    for (int noc = noc_f; noc <= noc_l; ++noc) {
	      const int ino = noocv_v[noc];
	      mu_sgs[icv] += mu_sgs_no[ino];
	    }
	    mu_sgs[icv] /= (double)(noc_l-noc_f+1);
	  }

	}

	delete[] mu_sgs_no;

	// ------------------------------------------------
	// update filtered mu_sgs into ghost cells...
	// ------------------------------------------------
	
	updateCvData(mu_sgs,REPLACE_DATA);
	
	// ------------------------------------------------
	// compute k_sgs based on filtered mu_sgs and turbulent Pr...
	// ------------------------------------------------
	
	FOR_ICV_G k_sgs[icv] = R_gas * gamma / ( gamma - 1.0) * mu_sgs[icv] / Pr_t;
	
      }
      break;
      
    default:
      
      if (mpi_rank == 0)
	cerr << "Error: unrecognized SGS_MODEL" << endl;
      throw(-1);
      
    }
    
    if (verbose) {
      dumpScalarRange(mu_sgs,ncv,"MU_SGS");
      dumpScalarRange(k_sgs,ncv,"K_SGS");
    }
    
  }
  
  void calcRhs(double * rho_rhs,double (*rhou_rhs)[3],double * rhoE_rhs,
	       const double *rho,const double (*rhou)[3],const double * rhoE,
	       const double rk_time,const int rk_step) {

    // update the cv-based primitive data...
    
    for (int icv = 0; icv < ncv_gf; icv++) {
      FOR_I3 u[icv][i] = rhou[icv][i]/rho[icv];
      p[icv] = (gamma-1.0)*( rhoE[icv] - 0.5*( rhou[icv][0]*u[icv][0] +
					       rhou[icv][1]*u[icv][1] + 
					       rhou[icv][2]*u[icv][2] ) );
      T[icv] = p[icv]/rho[icv]/R_gas;
      h[icv] = ( rhoE[icv] + p[icv] )/rho[icv];
    }

    // for shock detection, we are going to use the Ducros? switch:
    // theta = du_i/dx_i
    // vort = eps_ijk du_k/dx_j
    // s = -theta/(fabs(theta) + sqrt(vort_i*vort_i))
    // s > s_lim, use shock-capturing scheme.

    double (*grad_u0)[3][3] = new double[nfa][3][3];
    double (*grad_T0)[3] = new double[nfa][3];

    double (*grad_u1)[3][3] = new double[nfa][3][3];
    double (*grad_T1)[3] = new double[nfa][3];
    
    // left side first (note we skip boundary faces)...

    for (int ifa = nfa_b; ifa < nfa; ++ifa) {

      FOR_J3 {
	FOR_I3 grad_u0[ifa][i][j] = 0.0;
	grad_T0[ifa][j]           = 0.0;
      }
      
      const int cof_f = cvofa0_i[ifa];
      const int cof_l = cvofa0_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	FOR_J3 {
	  double coeff = cvofa0_grad[cof][j];
	  FOR_I3 grad_u0[ifa][i][j] += coeff*u[icv][i];
	  grad_T0[ifa][j]           += coeff*T[icv];
	}
      }
      
    }

    // copy the computed left state to the right state in
    // boundary/processor boundary faces...
    for (int ifa = nfa_b; ifa < nfa_bpi; ++ifa) {
      FOR_J3 {
	FOR_I3 grad_u1[ifa][i][j] = grad_u0[ifa][i][j];
	grad_T1[ifa][j]           = grad_T0[ifa][j];
      }
    }
    
    // exchange...
    updateFaR3(grad_u1,REPLACE_ROTATE_DATA);
    updateFaR2(grad_T1,REPLACE_ROTATE_DATA);
    
    // complete right states...
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
      
      FOR_J3 {
	FOR_I3 grad_u1[ifa][i][j] = 0.0;
	grad_T1[ifa][j]           = 0.0;
      }
      
      // cycle through the cell nbrs...
      const int cof_f = cvofa1_i[ifa];
      const int cof_l = cvofa1_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];
	FOR_J3 {
	  double coeff = cvofa1_grad[cof][j];
	  FOR_I3 grad_u1[ifa][i][j] += coeff*u[icv][i];
	  grad_T1[ifa][j]           += coeff*T[icv];
	}
      }
      
    }

    // ================================================================================
    // if this is the first of the rk steps, recompute the shock sensor - we want to
    // hold this on for all sub-steps of the rk cycle...
    // ================================================================================
    
    if (rk_step == 1) {

      // turn off on boundary faces...
      for (int ifa = 0; ifa < nfa_b; ++ifa) fa_ss[ifa] = 0.0; 
      
      // internal faces...
      for (int ifa = nfa_b; ifa < nfa; ++ifa) {

	// divergence...
	double theta = 0.5*( grad_u0[ifa][0][0] + grad_u0[ifa][1][1] + grad_u0[ifa][2][2] +
			     grad_u1[ifa][0][0] + grad_u1[ifa][1][1] + grad_u1[ifa][2][2] );

	// vorticity...
	// XXXXXX: should probably average this a bit...  
	double w_x = 0.25*( grad_u0[ifa][1][2] - grad_u0[ifa][2][1] +
			    grad_u1[ifa][1][2] - grad_u1[ifa][2][1] );
	double w_y = 0.25*( grad_u0[ifa][2][0] - grad_u0[ifa][0][2] +
			    grad_u1[ifa][2][0] - grad_u1[ifa][0][2] );
	double w_z = 0.25*( grad_u0[ifa][0][1] - grad_u0[ifa][1][0] +
			    grad_u1[ifa][0][1] - grad_u1[ifa][1][0] );

	// compute sensor...
	/*
	fa_ss[ifa] = -theta/( fabs(theta) + max(sqrt( w_x*w_x + w_y*w_y + w_z*w_z ),ducroisCoeff) );
	*/
	
	// use dillitation... 
	fa_ss[ifa] = theta;
	
      }

    }

    // set the fa_flag to 1 in faces that require shock-capturing reconstuction. This is
    // a bit safer than checking the fa_ss threshold locally because it synchronizes
    // interprocessor faces...
    
    FOR_IFA {
      // shock capturing based on abs mag of local dilitation - handles strong
      // expansion as well...
      if (fabs(fa_ss[ifa]) > dilCoeff) 
	fa_flag[ifa] = 1;
      else
	fa_flag[ifa] = 0;      
    }
    updateFaI1(fa_flag,ADD_DATA);
    
    //if (step < 2) FOR_IFA fa_flag[ifa] = 1;

    /*
    // grow faces by one layer?...
    // XXXXX - is this neccessary?...
    // extend faces by nbrs cells...
    FOR_ICV cv_flag[icv] = 0;
    for (int ifa = 0; ifa < nfa_bpi; ++ifa) if (fa_flag[ifa] > 0) {
      // these faces have only icv0 local...
      const int icv0 = cvofa[ifa][0];
      cv_flag[icv0] = 1;
    }
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa) if (fa_flag[ifa] > 0) {
      // these faces have both icv0 and icv1 local...
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      cv_flag[icv0] = 1;
      cv_flag[icv1] = 1;
    }
    updateCvData(cv_flag,REPLACE_DATA);
    for (int ifa = 0; ifa < nfa_b; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      if (cv_flag[icv0] > 0)
	fa_flag[ifa] = 1;
    }
    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      if ((cv_flag[icv0] > 0)||(cv_flag[icv1] > 0))
	fa_flag[ifa] = 1;
    }
    */

    // now any face with fa_flag > 0 should be treated with shock-capturing...
    
    for (int ifa = 0; ifa < nfa_b; ++ifa)
      fa_alpha_ss[ifa] = fa_alpha[ifa];
    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      if (fa_flag[ifa] > 0) {
	fa_alpha_ss[ifa] = 1.0;
      }
      else {
	fa_alpha_ss[ifa] = fa_alpha[ifa];
      }
    }

    // ========================================================
    // visualization
    // ========================================================
    if ((rk_step == 1)&&(step%check_interval == 0)) {
      
      // visualization...
      // update fa_ss_no to allow visualization by user...
      FOR_INO fa_ss_no[ino] = -1.0;
      FOR_ICV {
	const int foc_f = faocv_i[icv];
	const int foc_l = faocv_i[icv+1]-1;
	for (int foc = foc_f; foc <= foc_l; foc++) {
	  const int ifa = faocv_v[foc];
	  const int nof_f = noofa_i[ifa];
	  const int nof_l = noofa_i[ifa+1]-1;
	  for (int nof = nof_f; nof <= nof_l; ++nof) {
	    const int ino = noofa_v[nof];
	    fa_ss_no[ino] = max(fa_ss_no[ino],fa_ss[ifa]);
	  }
	}
      }
      updateNoR1(fa_ss_no,MAX_DATA);
      
    }    
    // ========================================================
    
    // -----------------------------------------------------------------
    // now face-data reconstruction including effect of shock sensor...
    // -----------------------------------------------------------------

    double *rho0 = new double[nfa];
    double (*u0)[3] = new double[nfa][3];
    double *p0 = new double[nfa];
    double *h0 = new double[nfa];

    double *rho1 = new double[nfa];
    double (*u1)[3] = new double[nfa][3];
    double *p1 = new double[nfa];
    double *h1 = new double[nfa];

    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      
      if (fa_flag[ifa] != 0) {
	
	// solution-weighted reconstruction for shock capturing...
	
	const int cof_f = cvofa0_i[ifa];
	const int cof_l = cvofa0_i[ifa+1]-1;
      
	// set the upwind (left) state - this is 
	// always the first cell in the list...

	const int icv0 = cvofa0_v[cof_f];
	rho0[ifa]         = rho[icv0];
	FOR_I3 u0[ifa][i] = u[icv0][i];
	p0[ifa]           = p[icv0];
	h0[ifa]           = h[icv0];

	// use some fraction of the square of the local value...
	// as swCoeff approaches zero, the scheme approaches first
	// order. By using this approach, traditionally problematic 
	// parts of the flux - i.e. zero density and pressure - 
	// become first order as well...
      
	double eps_drho2 = swCoeff*rho[icv0]*rho[icv0];
	double eps_du2   = swCoeff*gamma*p[icv0]/rho[icv0]; // use c here
	double eps_dp2   = swCoeff*p[icv0]*p[icv0];
	double eps_dh2   = swCoeff*h[icv0]*h[icv0]; 
      
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip first
	  const int icv = cvofa0_v[cof];
	  double coeff = cvofa0_value[cof];
	  // rho...
	  double delta = (rho[icv] - rho[icv0]); 
	  rho0[ifa] += eps_drho2/( delta*delta + eps_drho2 )*coeff*delta;
	  // u...
	  FOR_I3 {
	    delta = (u[icv][i] - u[icv0][i]); 
	    u0[ifa][i] += eps_du2/( delta*delta + eps_du2 )*coeff*delta;
	  }
	  // p...
	  delta = (p[icv] - p[icv0]); 
	  p0[ifa] += eps_dp2/( delta*delta + eps_dp2 )*coeff*delta;
	  // h...
	  delta = (h[icv] - h[icv0]); 
	  h0[ifa] += eps_dh2/( delta*delta + eps_dh2 )*coeff*delta;
	}

      }
      else {

	// normal reconstruction...

	rho0[ifa]         = 0.0;
	FOR_I3 u0[ifa][i] = 0.0;
	p0[ifa]           = 0.0;
	h0[ifa]           = 0.0;
      
	const int cof_f = cvofa0_i[ifa];
	const int cof_l = cvofa0_i[ifa+1]-1;
	for (int cof = cof_f; cof <= cof_l; ++cof) {
	  const int icv = cvofa0_v[cof];
	  double coeff = cvofa0_value[cof];
	  rho0[ifa]         += coeff*rho[icv];
	  FOR_I3 u0[ifa][i] += coeff*u[icv][i];
	  p0[ifa]           += coeff*p[icv];
	  h0[ifa]           += coeff*h[icv];
	}
	
      }
      
    }

    // copy the computed left state to the right state in
    // boundary/processor boundary faces...
    for (int ifa = nfa_b; ifa < nfa_bpi; ++ifa) {
      rho1[ifa]         = rho0[ifa];
      FOR_I3 u1[ifa][i] = u0[ifa][i];
      p1[ifa]           = p0[ifa];
      h1[ifa]           = h0[ifa];
    }
    
    // exchange...
    updateFaR1(rho1,REPLACE_DATA);
    updateFaR2(u1,REPLACE_ROTATE_DATA);
    updateFaR1(p1,REPLACE_DATA);
    updateFaR1(h1,REPLACE_DATA);
    
    // complete right states...
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
      
      if (fa_flag[ifa] != 0) {

	// solution-weighted reconstruction for shock capturing...
	
	const int cof_f = cvofa1_i[ifa];
	const int cof_l = cvofa1_i[ifa+1]-1;
      
	// set the upwind (right) state - this is 
	// always the first cell in the list...

	const int icv1 = cvofa1_v[cof_f];
	rho1[ifa]         = rho[icv1];
	FOR_I3 u1[ifa][i] = u[icv1][i];
	p1[ifa]           = p[icv1];
	h1[ifa]           = h[icv1];

	// use some fraction of the square of the local value...
	// as swCoeff approaches zero, the scheme approaches first
	// order. By using this approach, traditionally problematic 
	// parts of the flux - i.e. zero density and pressure - 
	// become first order as well...
      
	double eps_drho2 = swCoeff*rho[icv1]*rho[icv1];
	double eps_du2   = swCoeff*gamma*p[icv1]/rho[icv1]; // use c here
	double eps_dp2   = swCoeff*p[icv1]*p[icv1];
	double eps_dh2   = swCoeff*h[icv1]*h[icv1]; 
      
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip first
	  const int icv = cvofa1_v[cof];
	  double coeff = cvofa1_value[cof];
	  // rho...
	  double delta = (rho[icv] - rho[icv1]); 
	  rho1[ifa] += eps_drho2/( delta*delta + eps_drho2 )*coeff*delta;
	  // u...
	  FOR_I3 {
	    delta = (u[icv][i] - u[icv1][i]); 
	    u1[ifa][i] += eps_du2/( delta*delta + eps_du2 )*coeff*delta;
	  }
	  // p...
	  delta = (p[icv] - p[icv1]); 
	  p1[ifa] += eps_dp2/( delta*delta + eps_dp2 )*coeff*delta;
	  // h...
	  delta = (h[icv] - h[icv1]); 
	  h1[ifa] += eps_dh2/( delta*delta + eps_dh2 )*coeff*delta;
	}
	
      }
      else {
      
	// zero the left state...
	rho1[ifa]         = 0.0;
	FOR_I3 u1[ifa][i] = 0.0;
	p1[ifa]           = 0.0;
	h1[ifa]           = 0.0;
	
	// cycle through the cell nbrs...
	const int cof_f = cvofa1_i[ifa];
	const int cof_l = cvofa1_i[ifa+1]-1;
	for (int cof = cof_f; cof <= cof_l; ++cof) {
	  const int icv = cvofa1_v[cof];
	  double coeff = cvofa1_value[cof];
	  rho1[ifa]         += coeff*rho[icv];
	  FOR_I3 u1[ifa][i] += coeff*u[icv][i];
	  p1[ifa]           += coeff*p[icv];
	  h1[ifa]           += coeff*h[icv];
	}

      }
      
    }

    // now all left AND right face states are set on all
    // faces EXCEPT boundary faces (where a BC is required)...

    // ======================================
    // let the user specify a rhs source - if
    // there is no source, then this routine
    // zeros everything...
    // ======================================
    
    calcSourceHook(rho_rhs,rhou_rhs,rhoE_rhs,
		   rho,rhou,rhoE,rk_time,rk_step);

    // ======================================
    // the face based mu_sgs and k_sgs...
    // ======================================
    
    double * mu_sgs_fa = new double[nfa];
    double * k_sgs_fa = new double[nfa];

    // on boundary faces, take the internal cell value...
    for (int ifa = 0; ifa < nfa_b; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      mu_sgs_fa[ifa] = mu_sgs[icv0];
      k_sgs_fa[ifa]  = k_sgs[icv0];
    }
    // on inter-processor and internal faces, use the simple average...
    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      mu_sgs_fa[ifa] = 0.5*( mu_sgs[icv0] + mu_sgs[icv1] );
      k_sgs_fa[ifa]  = 0.5*( k_sgs[icv0]  + k_sgs[icv1] );
    }
    
    // ======================================
    // compute and assemble fluxes...
    // internal faces...
    // ======================================
    
    for (int ifa = nfa_b; ifa < nfa; ++ifa) {
      
      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
			  fa_normal[ifa][1]*fa_normal[ifa][1] +
			  fa_normal[ifa][2]*fa_normal[ifa][2] );
      
      double nx = fa_normal[ifa][0]/area;
      double ny = fa_normal[ifa][1]/area;
      double nz = fa_normal[ifa][2]/area;

      // ===============================
      // Jameson...
      // ===============================

      double rhoun = 0.25*(rho0[ifa]+rho1[ifa])*( nx*(u0[ifa][0]+u1[ifa][0]) +
						  ny*(u0[ifa][1]+u1[ifa][1]) +
						  nz*(u0[ifa][2]+u1[ifa][2]) );
      double Frho = (1.0-fa_alpha_ss[ifa])*rhoun;
      double Frhou[3];
      Frhou[0] = (1.0-fa_alpha_ss[ifa])*0.5*( rhoun*(u0[ifa][0]+u1[ifa][0]) + (p0[ifa]+p1[ifa])*nx );
      Frhou[1] = (1.0-fa_alpha_ss[ifa])*0.5*( rhoun*(u0[ifa][1]+u1[ifa][1]) + (p0[ifa]+p1[ifa])*ny );
      Frhou[2] = (1.0-fa_alpha_ss[ifa])*0.5*( rhoun*(u0[ifa][2]+u1[ifa][2]) + (p0[ifa]+p1[ifa])*nz );
      double FrhoE = (1.0-fa_alpha_ss[ifa])*0.5*( rhoun*( h0[ifa] + h1[ifa] ) );
      
      // ===============================
      // Roe flux...
      // ===============================

      // first part part of the convective (Euler) flux...
      double rhoun0 = rho0[ifa]*( u0[ifa][0]*nx + u0[ifa][1]*ny + u0[ifa][2]*nz );
      double rhoun1 = rho1[ifa]*( u1[ifa][0]*nx + u1[ifa][1]*ny + u1[ifa][2]*nz );
      Frho     += fa_alpha_ss[ifa]*0.5*(rhoun0 + rhoun1);
      Frhou[0] += fa_alpha_ss[ifa]*0.5*(rhoun0*u0[ifa][0] + rhoun1*u1[ifa][0] + (p0[ifa] + p1[ifa])*nx);
      Frhou[1] += fa_alpha_ss[ifa]*0.5*(rhoun0*u0[ifa][1] + rhoun1*u1[ifa][1] + (p0[ifa] + p1[ifa])*ny);
      Frhou[2] += fa_alpha_ss[ifa]*0.5*(rhoun0*u0[ifa][2] + rhoun1*u1[ifa][2] + (p0[ifa] + p1[ifa])*nz);
      FrhoE    += fa_alpha_ss[ifa]*0.5*(rhoun0*h0[ifa] + rhoun1*h1[ifa]);
      
      // original Roe (no entropy fix)...
      double sqrt_rho0 = sqrt(rho0[ifa]);
      double sqrt_rho1 = sqrt(rho1[ifa]);
      double tmp = 1.0/(sqrt_rho0 + sqrt_rho1 );
      double uAvg = tmp*(u0[ifa][0]*sqrt_rho0 + u1[ifa][0]*sqrt_rho1 );
      double vAvg = tmp*(u0[ifa][1]*sqrt_rho0 + u1[ifa][1]*sqrt_rho1 );
      double wAvg = tmp*(u0[ifa][2]*sqrt_rho0 + u1[ifa][2]*sqrt_rho1 );
      double hAvg = tmp*(h0[ifa]*sqrt_rho0 + h1[ifa]*sqrt_rho1 );
      
      double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
      double a2Avg = (gamma-1.0)*(hAvg - alphaAvg);
      if ( a2Avg < 1.0E-12 ) {
	cout << "Error: internal Roe flux failed at: " << x_fa[ifa][0] << " " << x_fa[ifa][1] << " " << x_fa[ifa][2] << endl;
	throw(-1);
      }
      double aAvg  = sqrt(a2Avg);
      double unAvg = uAvg*nx+ vAvg*ny+ wAvg*nz;

      double dr =  rho1[ifa] - rho0[ifa];
      double dru = rho1[ifa]*u1[ifa][0] - rho0[ifa]*u0[ifa][0];
      double drv = rho1[ifa]*u1[ifa][1] - rho0[ifa]*u0[ifa][1];
      double drw = rho1[ifa]*u1[ifa][2] - rho0[ifa]*u0[ifa][2];
      double drE = (rho1[ifa]*h1[ifa] - p1[ifa]) - (rho0[ifa]*h0[ifa] - p0[ifa]);
      
      double lam1 = fabs(unAvg);
      double lam2 = fabs(unAvg + aAvg);
      double lam3 = fabs(unAvg - aAvg);
      
      // entropy "fix"...
      {
	double eps = 0.5*( fabs(rhoun0/rho0[ifa] - rhoun1/rho1[ifa]) + fabs(sqrt(gamma*p0[ifa]/rho0[ifa]) - sqrt(gamma*p1[ifa]/rho1[ifa])));
	if (lam1 < 2.*eps)    lam1 = 0.25*lam1*lam1/eps + eps;
	if (lam2 < 2.*eps)    lam2 = 0.25*lam2*lam2/eps + eps;
	if (lam3 < 2.*eps)    lam3 = 0.25*lam3*lam3/eps + eps;
      }
      
      double abv1 = 0.5*(lam2 + lam3);
      double abv2 = 0.5*(lam2 - lam3);
      double abv3 = abv1 - lam1;
      double abv4 = (gamma-1.0)*(alphaAvg*dr - uAvg*dru - vAvg*drv - wAvg*drw + drE );
      double abv5 = unAvg*dr - nx*dru - ny*drv - nz*drw;
      double abv6 = abv3*abv4/a2Avg - abv2*abv5/aAvg;
      double abv7 = abv3*abv5 - abv2*abv4/aAvg;
      
      Frho     -= fa_alpha_ss[ifa]*0.5*(lam1*dr  + abv6);
      Frhou[0] -= fa_alpha_ss[ifa]*0.5*(lam1*dru + uAvg*abv6 - nx*abv7);
      Frhou[1] -= fa_alpha_ss[ifa]*0.5*(lam1*drv + vAvg*abv6 - ny*abv7);
      Frhou[2] -= fa_alpha_ss[ifa]*0.5*(lam1*drw + wAvg*abv6 - nz*abv7);
      FrhoE    -= fa_alpha_ss[ifa]*0.5*(lam1*drE + hAvg*abv6 - unAvg*abv7);
      
      // ===============================
      // skew-symmetric form...
      // ===============================

      /*
      double rhoun0 = rho0[ifa]*( u0[ifa][0]*nx + u0[ifa][1]*ny + u0[ifa][2]*nz );
      double rhoun1 = rho1[ifa]*( u1[ifa][0]*nx + u1[ifa][1]*ny + u1[ifa][2]*nz );
      
      double Frho     = 0.5*(rhoun0 + rhoun1);
      // add 1/2 the convetive term for both momentum and energy...
      double Frhou[3];
      Frhou[0] = 0.25*Frho*(u0[ifa][0]+u1[ifa][0]) + 0.5*(p0[ifa]+p1[ifa])*nx;
      Frhou[1] = 0.25*Frho*(u0[ifa][1]+u1[ifa][1]) + 0.5*(p0[ifa]+p1[ifa])*ny;
      Frhou[2] = 0.25*Frho*(u0[ifa][2]+u1[ifa][2]) + 0.5*(p0[ifa]+p1[ifa])*nz;
      double FrhoE = 0.25*Frho*(h0[ifa]+h1[ifa]);
      */

      // ===============================
      // viscous fluxes...
      // ===============================

      const double T_f = 0.5*(p0[ifa]/rho0[ifa] + p1[ifa]/rho1[ifa])/R_gas;
      const double mu_f = mu_ref*pow( T_f/T_ref, mu_power_law );
      const double mu_f_total = mu_f + mu_sgs_fa[ifa];
      const double k_f_total = R_gas * gamma / ( gamma - 1.0) * mu_f / Pr + k_sgs_fa[ifa]; // cp = R * gamma / (gamma - 1)
      
      // compute tau without the viscosity multiplier for now  because this term has the 
      // total viscosity in the momentum equation, but only the molecular viscosity in the
      // energy equation...

      double tau[3][3];
      FOR_I3 {
	FOR_J3 tau[i][j] = 0.5*( grad_u0[ifa][i][j] + grad_u0[ifa][j][i] +
				 grad_u1[ifa][i][j] + grad_u1[ifa][j][i] );
	// note 2/3 * 0.5 == 1/3...
	tau[i][i] -= 1.0/3.0*( grad_u0[ifa][0][0] + grad_u0[ifa][1][1] + grad_u0[ifa][2][2] +
			       grad_u1[ifa][0][0] + grad_u1[ifa][1][1] + grad_u1[ifa][2][2] );
      }
      
      double tauijnj[3] = {
	nx*tau[0][0] + ny*tau[0][1] + nz*tau[0][2],
	nx*tau[1][0] + ny*tau[1][1] + nz*tau[1][2],
	nx*tau[2][0] + ny*tau[2][1] + nz*tau[2][2] };
      
      Frhou[0] -= mu_f_total*tauijnj[0];
      Frhou[1] -= mu_f_total*tauijnj[1];
      Frhou[2] -= mu_f_total*tauijnj[2];
      
      // -k*dTdj*nj - ui*mu*tauij*nj
      FrhoE -= 
	0.5*k_f_total*( nx*(grad_T0[ifa][0]+grad_T1[ifa][0]) + 
			ny*(grad_T0[ifa][1]+grad_T1[ifa][1]) +
			nz*(grad_T0[ifa][2]+grad_T1[ifa][2]) ) +
	0.5*mu_f*( (u0[ifa][0]+u1[ifa][0])*tauijnj[0] + 
		   (u0[ifa][1]+u1[ifa][1])*tauijnj[1] + 
		   (u0[ifa][2]+u1[ifa][2])*tauijnj[2] );


      // ===============================
      // regular flux update...
      // ===============================

      int icv0 = cvofa[ifa][0];
      rho_rhs[icv0]            -= area*Frho;
      FOR_I3 rhou_rhs[icv0][i] -= area*Frhou[i];
      rhoE_rhs[icv0]           -= area*FrhoE;
      
      // icv1 may be ghost...
      int icv1 = cvofa[ifa][1];
      if (icv1 < ncv) {
	rho_rhs[icv1]            += area*Frho;
	FOR_I3 rhou_rhs[icv1][i] += area*Frhou[i];
	rhoE_rhs[icv1]           += area*FrhoE;
      }

      // ===============================
      // update fluxes including remaining skew-symmetric terms...
      // ===============================

      /*
      int icv0 = cvofa[ifa][0];
      rho_rhs[icv0]            -= area*Frho;
      FOR_I3 rhou_rhs[icv0][i] -= area*( Frhou[i] + 
					 0.5*u[icv0][i]*Frho + 
					 0.25*(rhou[icv0][0]*nx + rhou[icv0][1]*ny + rhou[icv0][2]*nz)*
					 (u0[ifa][i]+u1[ifa][i]) );
      rhoE_rhs[icv0]           -= area*( FrhoE +
					 0.5*h[icv0]*Frho + 
					 0.25*(rhou[icv0][0]*nx + rhou[icv0][1]*ny + rhou[icv0][2]*nz)*
					 (h0[ifa]+h1[ifa]) );
      
      // icv1 may be ghost...
      int icv1 = cvofa[ifa][1];
      if (icv1 < ncv) {
	rho_rhs[icv1]            += area*Frho;
	FOR_I3 rhou_rhs[icv1][i] += area*( Frhou[i] + 
					   0.5*u[icv1][i]*Frho + 
					   0.25*(rhou[icv1][0]*nx + rhou[icv1][1]*ny + rhou[icv1][2]*nz)*
					   (u0[ifa][i]+u1[ifa][i]) );
	rhoE_rhs[icv1]           += area*( FrhoE +
					   0.5*h[icv1]*Frho + 
					   0.25*(rhou[icv1][0]*nx + rhou[icv1][1]*ny + rhou[icv1][2]*nz)*
					   (h0[ifa]+h1[ifa]) );
      }
      */
      
    }

    // =========================================================
    // boundary conditions...
    // =========================================================

    int bc_err = 0;
    
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
      if (zone->getKind() == FA_ZONE_BOUNDARY) {
	
	Param * param;
	if ( getParam(param,zone->getName()) ) {
	  
	  if ( param->getString() == "WALL_ISOTHERMAL" ) {
	    
	    // the fake cell properties should already be set such that they are
	    // consistent with the isothermal bc...
	    
	    double Twall = param->getDouble(2);
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      // geometry...
	      
	      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
				  fa_normal[ifa][1]*fa_normal[ifa][1] +
				  fa_normal[ifa][2]*fa_normal[ifa][2] );
	      
	      double nx = fa_normal[ifa][0]/area;
	      double ny = fa_normal[ifa][1]/area;
	      double nz = fa_normal[ifa][2]/area;

	      // zero the left state...
	      p0[ifa]           = 0.0;
	      
	      // and the left gradients...
	      FOR_J3 {
		FOR_I3 grad_u0[ifa][i][j] = 0.0;
		grad_T0[ifa][j]           = 0.0;
	      }
	      
	      // cycle through the cell nbrs...
	      const int cof_f = cvofa0_i[ifa];
	      const int cof_l = cvofa0_i[ifa+1]-1;
	      for (int cof = cof_f; cof <= cof_l; ++cof) {
		const int icv = cvofa0_v[cof];
		double coeff = cvofa0_value[cof];
		p0[ifa]           += coeff*p[icv];
		FOR_J3 {
		  coeff = cvofa0_grad[cof][j];
		  FOR_I3 grad_u0[ifa][i][j] += coeff*u[icv][i];
		  grad_T0[ifa][j]           += coeff*T[icv];
		}
	      }

	      // euler part is just the pressure...

	      double Frho,Frhou[3],FrhoE;
	      Frho = 0.0;
	      Frhou[0] = p0[ifa]*nx;
	      Frhou[1] = p0[ifa]*ny;
	      Frhou[2] = p0[ifa]*nz;
	      FrhoE = 0.0;

	      // viscous...
	      
	      const double mu_f = mu_ref*pow( Twall/T_ref, mu_power_law );
	      const double mu_f_total = mu_f + mu_sgs_fa[ifa];
	      const double k_f_total = R_gas * gamma / ( gamma - 1.0) * mu_f / Pr + k_sgs_fa[ifa]; // cp = R * gamma / (gamma - 1)
	      
	      // compute tau without the viscosity multiplier for now - because this term has the 
	      // total viscosity in the momentum equation, but only the molecular viscosity in the
	      // energy equation...
      
	      double tau[3][3];
	      FOR_I3 {
		FOR_J3 tau[i][j] = ( grad_u0[ifa][i][j] + grad_u0[ifa][j][i] );
		tau[i][i] -= 2.0/3.0*( grad_u0[ifa][0][0] + grad_u0[ifa][1][1] + grad_u0[ifa][2][2] );
	      }
	      
	      double tauijnj[3] = {
		nx*tau[0][0] + ny*tau[0][1] + nz*tau[0][2],
		nx*tau[1][0] + ny*tau[1][1] + nz*tau[1][2],
		nx*tau[2][0] + ny*tau[2][1] + nz*tau[2][2] };
	      
	      Frhou[0] -= mu_f_total*tauijnj[0];
	      Frhou[1] -= mu_f_total*tauijnj[1];
	      Frhou[2] -= mu_f_total*tauijnj[2];
	      
	      // -k*dTdj*nj - ui*mu*tauij*nj, ui == 0...
	      FrhoE -= 
		k_f_total*( nx*(grad_T0[ifa][0]) + 
			    ny*(grad_T0[ifa][1]) +
			    nz*(grad_T0[ifa][2]) );
	      // update...	    
	      int icv0 = cvofa[ifa][0];
	      rho_rhs[icv0]            -= area*Frho;
	      FOR_I3 rhou_rhs[icv0][i] -= area*Frhou[i];
	      rhoE_rhs[icv0]           -= area*FrhoE;
	      
	    }
	    
	  }
	  else if ( param->getString() == "WALL_ADIABATIC" ) {
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      // geometry...
	      
	      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
				  fa_normal[ifa][1]*fa_normal[ifa][1] +
				  fa_normal[ifa][2]*fa_normal[ifa][2] );
	      
	      double nx = fa_normal[ifa][0]/area;
	      double ny = fa_normal[ifa][1]/area;
	      double nz = fa_normal[ifa][2]/area;

	      // zero the left state...
	      rho0[ifa]         = 0.0;
	      p0[ifa]           = 0.0;
	      
	      // and the left gradients...
	      FOR_J3 {
		FOR_I3 grad_u0[ifa][i][j] = 0.0;
	      }
	      
	      // cycle through the cell nbrs...
	      const int cof_f = cvofa0_i[ifa];
	      const int cof_l = cvofa0_i[ifa+1]-1;
	      for (int cof = cof_f; cof <= cof_l; ++cof) {
		const int icv = cvofa0_v[cof];
		double coeff = cvofa0_value[cof];
		rho0[ifa]         += coeff*rho[icv];
		p0[ifa]           += coeff*p[icv];
		FOR_J3 {
		  coeff = cvofa0_grad[cof][j];
		  FOR_I3 grad_u0[ifa][i][j] += coeff*u[icv][i];
		}
	      }

	      // euler part is just the pressure...

	      double Frho,Frhou[3],FrhoE;
	      Frho = 0.0;
	      Frhou[0] = p0[ifa]*nx;
	      Frhou[1] = p0[ifa]*ny;
	      Frhou[2] = p0[ifa]*nz;
	      FrhoE = 0.0;
	      
	      // viscous...

	      const double T_f = p0[ifa]/rho0[ifa]/R_gas;
	      const double mu_f = mu_ref*pow( T_f/T_ref, mu_power_law );
	      const double mu_f_total = mu_f + mu_sgs_fa[ifa];
	      
	      // compute tau without the viscosity multiplier for now - because this term has the 
	      // total viscosity in the momentum equation, but only the molecular viscosity in the
	      // energy equation...
      
	      double tau[3][3];
	      FOR_I3 {
		FOR_J3 tau[i][j] = ( grad_u0[ifa][i][j] + grad_u0[ifa][j][i] );
		tau[i][i] -= 2.0/3.0*( grad_u0[ifa][0][0] + grad_u0[ifa][1][1] + grad_u0[ifa][2][2] );
	      }
	      
	      double tauijnj[3] = {
		nx*tau[0][0] + ny*tau[0][1] + nz*tau[0][2],
		nx*tau[1][0] + ny*tau[1][1] + nz*tau[1][2],
		nx*tau[2][0] + ny*tau[2][1] + nz*tau[2][2] };
	      
	      Frhou[0] -= mu_f_total*tauijnj[0];
	      Frhou[1] -= mu_f_total*tauijnj[1];
	      Frhou[2] -= mu_f_total*tauijnj[2];
	      
	      // update...	    
	      int icv0 = cvofa[ifa][0];
	      rho_rhs[icv0]            -= area*Frho;
	      FOR_I3 rhou_rhs[icv0][i] -= area*Frhou[i];
	      rhoE_rhs[icv0]           -= area*FrhoE;
	      
	    }
	    
	  }	
	  else if (( param->getString() == "CBC" )||( param->getString() == "HOOK" )) {
	    
	    // for a CBC (characteristic BC) or a hook bc, the fake data is already populated with 
	    // the desired values...
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
				  fa_normal[ifa][1]*fa_normal[ifa][1] +
				  fa_normal[ifa][2]*fa_normal[ifa][2] );
	      double unitNormal[3] = { fa_normal[ifa][0]/area,
				       fa_normal[ifa][1]/area,	
				       fa_normal[ifa][2]/area };
      
	      const int icv0 = cvofa[ifa][0];
	      assert( (icv0 >= 0)&&(icv0 < ncv) ); // icv0 is local
	      
	      const int icv1 = cvofa[ifa][1];
	      assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) ); // icv1 is "fake"
	      
	      double Frho,Frhou[3],FrhoE;
	      if ( calcEulerFluxRoe(Frho,Frhou,FrhoE,
				    rho[icv0],rhou[icv0],rhoE[icv0],
				    rho[icv1],rhou[icv1],rhoE[icv1],
				    unitNormal) != 0 ) {
		cerr << "Error: Flux failed for boundary face at: " << x_fa[ifa][0] << " " << x_fa[ifa][1] << " " << x_fa[ifa][2] << endl;
		throw(-1);
	      }
	      
	      Frho            *= area;
	      FOR_I3 Frhou[i] *= area;
	      FrhoE           *= area;
	      
	      rho_rhs[icv0]            -= Frho;
	      FOR_I3 rhou_rhs[icv0][i] -= Frhou[i];
	      rhoE_rhs[icv0]           -= FrhoE;
	      
	    }
	    
	  }
	  else if ( param->getString() == "SYMMETRY" ) {
	    
	    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	      
	      // geometry...
	      
	      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
				  fa_normal[ifa][1]*fa_normal[ifa][1] +
				  fa_normal[ifa][2]*fa_normal[ifa][2] );
	      
	      double nx = fa_normal[ifa][0]/area;
	      double ny = fa_normal[ifa][1]/area;
	      double nz = fa_normal[ifa][2]/area;
	      
	      int icv0 = cvofa[ifa][0]; // needed for cv pressure...
	      
	      double Frho = 0.0;
	      double Frhou[3];
	      Frhou[0] = p[icv0]*nx;
	      Frhou[1] = p[icv0]*ny;
	      Frhou[2] = p[icv0]*nz;
	      double FrhoE = 0.0;

	      Frho            *= area;
	      FOR_I3 Frhou[i] *= area;
	      FrhoE           *= area;
	      
	      rho_rhs[icv0]            -= Frho;
	      FOR_I3 rhou_rhs[icv0][i] -= Frhou[i];
	      rhoE_rhs[icv0]           -= FrhoE;
	      
	    }

	  }
	  else {
	    
	    if (mpi_rank == 0)
	      cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
	    bc_err = -1;
	    
	  }
	}
	else {
	  
	  if (mpi_rank == 0)
	    cerr << "Error: cannot determine bc for zone: " << zone->getName() << endl;
	  bc_err = -1;
	  
	}
	
      }
    }
    
    if (bc_err != 0) {
      cout << "there was a bc problem" << endl;
      throw(-1);
    }    

    delete[] mu_sgs_fa;
    delete[] k_sgs_fa;
    
    delete[] rho0;
    delete[] u0;
    delete[] p0;
    delete[] h0;

    delete[] rho1;
    delete[] u1;
    delete[] p1;
    delete[] h1;
    
    delete[] grad_u0;
    delete[] grad_T0;
    
    delete[] grad_u1;
    delete[] grad_T1;
    
  }

  virtual void calcSourceHook(double * rho_rhs,double (*rhou_rhs)[3],double * rhoE_rhs,
			      const double *rho,const double (*rhou)[3],const double * rhoE,
			      const double rk_time,const int rk_step) {

    
    // if you overload this routine, you are responsible for zeroing everthing
    // that does not have a source. i.e.
    
    FOR_ICV rho_rhs[icv] = 0.0;
    FOR_ICV FOR_I3 rhou_rhs[icv][i] = 0.0;
    FOR_ICV rhoE_rhs[icv] = 0.0;
    
  }
  
  int calcEulerFluxRoe(double &Frho, double *Frhou, double &FrhoE,
		       const double& rho0, const double *rhou0, const double& rhoE0,
		       const double& rho1, const double *rhou1, const double& rhoE1,
		       const double *nVec) {
    
    double gm1 = gamma - 1.0;
    
    double u0[3] = { rhou0[0]/rho0,
		     rhou0[1]/rho0,
		     rhou0[2]/rho0 };

    double u1[3] = { rhou1[0]/rho1,
		     rhou1[1]/rho1,
		     rhou1[2]/rho1 };

    double rhoun0 = vecDotVec3d(rhou0, nVec);
    double rhoun1 = vecDotVec3d(rhou1, nVec);

    double p0 = gm1*( rhoE0 - 0.5*vecDotVec3d(rhou0, u0) );
    double p1 = gm1*( rhoE1 - 0.5*vecDotVec3d(rhou1, u1) );
    
    /*
    p0 = max(1.0e-12,p0);
    p1 = max(1.0e-12,p1);
    */
    
    double H0 = (rhoE0 + p0)/rho0;
    double H1 = (rhoE1 + p1)/rho1;

    // first part part of the convective (Euler) flux...
    Frho     = 0.5*(rhoun0 + rhoun1 );
    Frhou[0] = 0.5*(rhoun0*u0[0]+ rhoun1*u1[0]+ (p0 + p1 )*nVec[0]);
    Frhou[1] = 0.5*(rhoun0*u0[1]+ rhoun1*u1[1]+ (p0 + p1 )*nVec[1]);
    Frhou[2] = 0.5*(rhoun0*u0[2]+ rhoun1*u1[2]+ (p0 + p1 )*nVec[2]);
    FrhoE    = 0.5*(rhoun0*H0 + rhoun1*H1 );

    // original Roe (no entropy fix)...
    double sqrt_rho0 = sqrt(rho0);
    double sqrt_rho1 = sqrt(rho1);
    double tmp = 1.0/(sqrt_rho0 + sqrt_rho1 );
    double uAvg = tmp*(u0[0]*sqrt_rho0 + u1[0]*sqrt_rho1 );
    double vAvg = tmp*(u0[1]*sqrt_rho0 + u1[1]*sqrt_rho1 );
    double wAvg = tmp*(u0[2]*sqrt_rho0 + u1[2]*sqrt_rho1 );
    double hAvg = tmp*(H0*sqrt_rho0 + H1*sqrt_rho1 );
    // modified steger warming
    /*double uAvg = 0.5*(u0[0] + u1[0] );
    double vAvg = 0.5*(u0[1] + u1[1] );
    double wAvg = 0.5*(u0[2] + u1[2] );
    double hAvg = 0.5*(H0 + H1 );*/
    
    double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
    double a2Avg = gm1*(hAvg - alphaAvg);
    if ( a2Avg < 1.0E-12 )
      return(-1);
    //a2Avg = max(1.0E-12,a2Avg);
    double aAvg  = sqrt(a2Avg);
    double unAvg = uAvg*nVec[0]+ vAvg*nVec[1]+ wAvg*nVec[2];

    double dr = rho1 - rho0;
    double dru = rhou1[0] - rhou0[0];
    double drv = rhou1[1] - rhou0[1];
    double drw = rhou1[2] - rhou0[2];
    double drE = rhoE1 - rhoE0;

    double lam1 = fabs(unAvg);
    double lam2 = fabs(unAvg + aAvg);
    double lam3 = fabs(unAvg - aAvg);

    // entropy "fix"...
    {
      double eps = 0.5*( fabs(rhoun0/rho0 - rhoun1/rho1) + fabs(sqrt(gamma*p0/rho0) - sqrt(gamma*p1/rho1)));
      if (lam1 < 2.*eps)    lam1 = 0.25*lam1*lam1/eps + eps;
      if (lam2 < 2.*eps)    lam2 = 0.25*lam2*lam2/eps + eps;
      if (lam3 < 2.*eps)    lam3 = 0.25*lam3*lam3/eps + eps;
    }

    double abv1 = 0.5*(lam2 + lam3);
    double abv2 = 0.5*(lam2 - lam3);
    double abv3 = abv1 - lam1;
    double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv - wAvg*drw + drE );
    double abv5 = unAvg*dr - nVec[0]*dru - nVec[1]*drv - nVec[2]*drw;
    double abv6 = abv3*abv4/a2Avg - abv2*abv5/aAvg;
    double abv7 = abv3*abv5 - abv2*abv4/aAvg;
    
    Frho     -= 0.5*(lam1*dr  + abv6);
    Frhou[0] -= 0.5*(lam1*dru + uAvg*abv6 - nVec[0]*abv7);
    Frhou[1] -= 0.5*(lam1*drv + vAvg*abv6 - nVec[1]*abv7);
    Frhou[2] -= 0.5*(lam1*drw + wAvg*abv6 - nVec[2]*abv7);
    FrhoE    -= 0.5*(lam1*drE + hAvg*abv6 - unAvg*abv7);

    return(0);

  }

      
};

class EulerCf4 : public Charles {
  
private:
  
  int EULER_VORTEX_PERIODIC_FLAG;
  double EULER_VORTEX_U_INF;

  double * rho_error;
  
public:
  
  // Constructor that takes input file...
  EulerCf4(char *name) : Charles(name) {
    
    if (mpi_rank == 0)
      cout << "EulerCf4()" << endl;
    
    EULER_VORTEX_PERIODIC_FLAG = 0; // turn off for yaser prob...
    EULER_VORTEX_U_INF = 0.5;
    
    rho_error = NULL; registerScalar(rho_error, "RHO_ERROR", CV_DATA);
    
  }
  
  void initialHook() {
    
    if (step == 0) {
      
      if (mpi_rank == 0)
	cout << " > EulerCf4: initialHook(): " << endl;

      calcInviscidVortexLocal(rho,rhou,rhoE,
			 x_cv,ncv,0.0,gamma,rho_ref,p_ref,
			 EULER_VORTEX_U_INF,EULER_VORTEX_PERIODIC_FLAG);
      
    }

  }

  int bcHook(double * rho_bc,double (*rhou_bc)[3],double * rhoE_bc,
	     FaZone& zone,const double rk_time) {

    calcInviscidVortexLocal(&(rho_bc[zone.ifa_f]),&(rhou_bc[zone.ifa_f]),&(rhoE_bc[zone.ifa_f]),
			    &(x_fa[zone.ifa_f]),zone.ifa_l-zone.ifa_f+1,rk_time,gamma,rho_ref,p_ref,
			    EULER_VORTEX_U_INF,EULER_VORTEX_PERIODIC_FLAG);
    
    return(0);
    
  }

  void transformMeshHook() {
    
    cout << "WARNING: **************** transforming mesh ********************" << endl;
    
    // yaser's test mesh is x: 0-4, y: -1..1, z: -1..1
    
    FOR_INO {
      // shift and scale x...
      x_no[ino][0] = (x_no[ino][0] - 2.0)*20.0;
      x_no[ino][1] = x_no[ino][1]*20.0;
      x_no[ino][2] = x_no[ino][2]*20.0;
    }
    
  }

  /*
  void transformMeshHook() {
    
    if (mpi_rank == 0)
      cout << "XXXXX added shear to grid" << endl;
    
    //FOR_INO x_no[ino][0] += x_no[ino][1];

    // twisted sinusoidal transform...
    // ARB 2004...

    FOR_INO {
      double x = x_no[ino][0]/5.0;
      double y = x_no[ino][1]/5.0;
      x_no[ino][0] = 5.0*( x + 0.15*sin(M_PI*y) );
      x_no[ino][1] = 5.0*( y + 0.15*sin(M_PI*x) );
    }
    
  }
  */
  
  // this is in MiscUtils also, but does not inline for some reason...
  
  inline void calcInviscidVortexLocal(double * rho,double (*rhou)[3],double * rhoE,
			       const double (*x)[3],const int n,
			       const double t,const double gamma,
			       const double rho_inf,const double p_inf,const double u_inf,
			       const int periodic_flag) {

    // most grids have -5..5 in x...
    const double xmin = -5.0;
    const double xmax =  5.0;

    // except the prism grids, which have (-5..5)*cos(30)...
    //const double xmin = -4.330127019;
    //const double xmax =  4.330127019;

    // x,y position of the vortex
    // centered
    const double x0 = -32.0;
    const double y0 = 10.0;
    // near top right - for checking bc quickly...
    //double x0 = 3.0;
    //double y0 = 2.0;

    // direction of propagation
    //const double theta = M_PI/3.0;
    const double theta = 0.0;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const double Ma_inf = u_inf/sqrt(gamma*p_inf/rho_inf);
    const double rc = 1.0;

    // circulation parameter...
    //const double e_twopi = 0.001; // very weak
    //const double e_twopi = 0.005;
    const double e_twopi = 0.08; // normal
    //double e_twopi = 0.1;
    //double e_twopi = 0.4; //very strong 
    //double e_twopi = 1.0; //very very strong 

    // setup...
    const double coeff = 0.5 * e_twopi*e_twopi * (gamma-1.0) * Ma_inf*Ma_inf;
    for (int i = 0; i < n; i++) {
      
      double dx = x[i][0] - x0 - u_inf*cos_theta*t;
      double dy = x[i][1] - y0 - u_inf*sin_theta*t;
    
      // for periodic -5 < x < 5, -5 < y < 5, bring the vortex back
      // into the domain...
      if (periodic_flag == 1) {
        while(dx > xmax) dx -= (xmax-xmin);
        while(dx < xmin) dx += (xmax-xmin);
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

      // HACK - uniform flow...
      /*
      rho[i] = rho_inf;
      rhou[i][0] = rho_inf*u_inf*cos_theta;
      rhou[i][1] = rho_inf*u_inf*sin_theta;
      rhou[i][2] = rho_inf*u_inf*0.01;
      rhoE[i] = p_inf/(gamma-1.0) + 0.5*(rhou[i][0]*rhou[i][0] + rhou[i][1]*rhou[i][1] + rhou[i][2]*rhou[i][2])/rho[i];
      */

    }

  }

  void temporalHook() {
		     
    if (step%check_interval == 0) {
      
      if (mpi_rank == 0)
	cout << "Errors:" << endl;

      double * rho_exact = new double[ncv];
      double (*rhou_exact)[3] = new double[ncv][3];
      double * rhoE_exact = new double[ncv];
      
      double my_l2[6],my_linf[5];
      for (int i = 0; i < 6; i++)
	my_l2[i] = 0.0;
      for (int i = 0; i < 5; i++)
	my_linf[i] = 0.0;

      calcInviscidVortexLocal(rho_exact,rhou_exact,rhoE_exact,
			 x_cv,ncv,time,gamma,rho_ref,p_ref,
			 EULER_VORTEX_U_INF,EULER_VORTEX_PERIODIC_FLAG);
      
      for (int icv = 0; icv < ncv; icv++) {
	
	double delta = rho[icv] - rho_exact[icv];
	
	rho_error[icv] = sqrt(delta*delta);

	my_l2[0] += cv_volume[icv]*delta*delta;
	my_linf[0] = max( fabs(delta), my_linf[0] );
	
	for (int i = 0; i < 3; i++) {
	  delta = rhou[icv][i] - rhou_exact[icv][i];
	  my_l2[i+1] += cv_volume[icv]*delta*delta;
	  my_linf[i+1] = max( fabs(delta), my_linf[i+1] );
	}
	
	delta = rhoE[icv] - rhoE_exact[icv];
	my_l2[4] += cv_volume[icv]*delta*delta;
	my_linf[4] = max( fabs(delta), my_linf[4] );

	my_l2[5] += cv_volume[icv];

      }
      double l2[6],linf[5];
      MPI_Reduce(my_l2,l2,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      MPI_Reduce(my_linf,linf,5,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
	cout << " > time, L2: " << time << " " << 
	  sqrt(l2[0]/l2[5]) << " " <<
	  sqrt(l2[1]/l2[5]) << " " <<
	  sqrt(l2[2]/l2[5]) << " " <<
	  sqrt(l2[3]/l2[5]) << " " <<
	  sqrt(l2[4]/l2[5]) << endl;
	cout << " > time, Linf: " << time << " " <<
	  linf[0] << " " <<
	  linf[1] << " " <<
	  linf[2] << " " <<
	  linf[3] << " " <<
	  linf[4] << endl;
      }
      
      delete[] rho_exact;
      delete[] rhou_exact;
      delete[] rhoE_exact;
      
    }

  }
  
};


class ChannelCf4 : public Charles {
  
private:

  double (*u_avg)[3]; 
  double (*u_rms)[3]; 
  double (*u_rey)[3]; 
  double (*p_avg);
  double (*p_rms);
  double (*T_avg);
  double (*T_rms);
  double stats_time;
 
  double * part;
 
public:
  
  // Constructor that takes input file...
  ChannelCf4(char *name) : Charles(name) {
    
    if (mpi_rank == 0)
      cout << "ChannelCf4()" << endl;
    
    registerValue(stats_time,"STATS_TIME");
    
    u_avg = NULL; registerVector(u_avg, "U_AVG", CV_DATA);
    u_rms = NULL; registerVector(u_rms, "U_RMS", CV_DATA);
    u_rey = NULL; registerVector(u_rey, "U_REY", CV_DATA);
    
    p_avg = NULL; registerScalar(p_avg, "P_AVG", CV_DATA);
    p_rms = NULL; registerScalar(p_rms, "P_RMS", CV_DATA);

    T_avg = NULL; registerScalar(T_avg, "T_AVG", CV_DATA);
    T_rms = NULL; registerScalar(T_rms, "T_RMS", CV_DATA);

    part = NULL; registerScalar(part,"PART",NO_DATA);
    
  }

  
  void initialHook() {
    
    if (step == 0) {

      if (mpi_rank == 0)
	cout << " > ChannelCf4: initialHook()" << endl;
      
      /*
      // laminar...
      FOR_ICV laminarChannelIC(rho[icv],rhou[icv],rhoE[icv],
      x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
      */
      
      // turbulent...
      FOR_ICV turbChannelIC(rho[icv],rhou[icv],rhoE[icv],
			    x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
      
    }

    if (checkParam("RESET_STATS")) {
      
      if (mpi_rank == 0)
	cout << " > resetting stats" << endl;

      stats_time = 0.0;

      FOR_ICV p_avg[icv] = 0.0;
      FOR_ICV p_rms[icv] = 0.0;

      FOR_ICV T_avg[icv] = 0.0;
      FOR_ICV T_rms[icv] = 0.0;

      FOR_ICV FOR_I3 u_avg[icv][i] = 0.0;
      FOR_ICV FOR_I3 u_rms[icv][i] = 0.0;
      FOR_ICV FOR_I3 u_rey[icv][i] = 0.0;

    }

    FOR_INO part[ino] = (double)mpi_rank;
    
    // ==========================================
    // initialize the x/z averaging operator...
    // ==========================================
    
    if (average == NULL) {

      if (mpi_rank == 0)
	cout << " > setting up average for y stats." << endl;

      double * y = new double[ncv];
      FOR_ICV y[icv] = x_cv[icv][1];
      average = new Average(y,cv_volume,ncv,1.0E-7);
      average->setKind(AVERAGE_Y);
      delete[] y;
    }
    else {
      
      // is some other part of the code has set up this averaging operator,
      // make sure it is as expected...
      assert( average->getKind() == AVERAGE_Y );
    }
    
  }
  
  void laminarChannelIC(double &rho,double *rhou,double &rhoE,
			const double x,const double y,const double z) {
    
    double p = R_gas*rho_ref*T_ref;
    double T = T_ref + (gamma-1.0) * Pr / ( 12.0 * pow( mu_ref, 2.0 ) * R_gas * gamma ) * ( 1.0 - y*y*y*y );
    rho = p / ( R_gas * T );
    rhou[0] = rho*( 0.5 / mu_ref * ( 1.0 - y*y ) );
    rhou[1] = 0.0;
    rhou[2] = 0.0;
    rhoE = p/(gamma - 1.0) + 
      0.5*( rhou[0]*rhou[0] + 
	    rhou[1]*rhou[1] + 
	    rhou[2]*rhou[2] )/rho;
    
  }

  void turbChannelIC(double &rho,double *rhou,double &rhoE,
		     const double x,const double y,const double z) {
    
    // approximate turbulent mean profile...
    rho = rho_ref;
    rhou[0] = rho_ref*25.0*(1.0 + y*y)*(1.0 - y*y);
    rhou[1] = 0.0;
    rhou[2] = 0.0;
    
    // big perturbation...
    rhou[0] += 10*sin(10.2*x + 13.1*y + 11.7*z);
    rhou[1] += 10*sin(12.2*x + 11.1*y + 12.7*z);
    rhou[2] += 10*sin(14.2*x + 10.1*y + 9.7*z);

    rhoE = p_ref/(gamma - 1.0) + 
      0.5*( rhou[0]*rhou[0] + 
	    rhou[1]*rhou[1] + 
	    rhou[2]*rhou[2] )/rho;
  }
  
  void calcSourceHook(double * rho_rhs,double (*rhou_rhs)[3],double * rhoE_rhs,
		      const double *rho,const double (*rhou)[3],const double * rhoE,
		      const double rk_time,const int rk_step) {
    
    // unit forcing in the x-direction...
    
    FOR_ICV {
      rho_rhs[icv]     = 0.0;
      rhou_rhs[icv][0] = rho[icv]*cv_volume[icv];
      rhou_rhs[icv][1] = 0.0;
      rhou_rhs[icv][2] = 0.0;
      rhoE_rhs[icv]    = rhou[icv][0]*cv_volume[icv];
    }
    
  }
      
  void temporalHook() {
    
    // only add stats every 10 time steps...
    
    if (step%10 == 0) {
      
      // do stats...

      // ========================
      // Pressure
      // ========================
      
      double * r1 = new double[ncv];
      
      // put pressure in r1...
      FOR_ICV r1[icv] = (gamma-1.0)*( rhoE[icv] - 0.5*( rhou[icv][0]*rhou[icv][0] +
							rhou[icv][1]*rhou[icv][1] +
							rhou[icv][2]*rhou[icv][2] )/rho[icv] );      
      addStats(p_avg,p_rms,r1,ncv,stats_time,dt);
      
      // ========================
      // T = p / (rho*R_gas)
      // ========================
      
      FOR_ICV r1[icv] /= rho[icv]*R_gas;
      
      addStats(T_avg,T_rms,r1,ncv,stats_time,dt);
      
      delete[] r1;
      
      // ========================
      // velocity...
      // ========================
      
      double (*r2)[3] = new double[ncv][3];
      
      FOR_ICV FOR_I3 r2[icv][i] = rhou[icv][i]/rho[icv];
      
      addStats(u_avg,u_rms,u_rey,r2,ncv,stats_time,dt);
      
      delete[] r2;

      // ========================
      // update stats time
      // ========================
      
      stats_time += dt;

    }


  }
  
  void addStats(double * r1_avg,double * r1_rms,
		const double * r1,const int n,const double old_stats_time,const double dt) {
    
    double stats_time = old_stats_time + dt;
    
    for (int i = 0; i < n; ++i) {
      // put sum( phi**2 * dt )_new / sum( dt ) in r1_rms[i]...
      r1_rms[i] = ( old_stats_time*( r1_rms[i]*r1_rms[i] + r1_avg[i]*r1_avg[i] ) + dt*(r1[i]*r1[i]) ) / stats_time;
      // put new means in r1_avg...
      r1_avg[i] = ( old_stats_time*r1_avg[i] + dt*(r1[i]) ) / stats_time;
    }

    // apply additional averaging...
    if (average != NULL) {
      average->apply( r1_avg );
      average->apply( r1_rms );
    }
    
    for (int i = 0; i < n; ++i) {
      // complete rms...
      r1_rms[i] = sqrt( fabs( r1_rms[i] - r1_avg[i]*r1_avg[i] ) );
    }
    
  }
  
  void addStats(double (*r2_avg)[3],double (*r2_rms)[3],double (*r2_rey)[3], 
		const double (*r2)[3],const int n,const double old_stats_time,const double dt) {
    
    double stats_time = old_stats_time + dt;
    
    for (int i = 0; i < n; ++i) {
      // put sum( ui**2 * dt )_new / sum( dt ) in r2_rms(i)...
      r2_rms[i][0] = ( old_stats_time*( r2_rms[i][0]*r2_rms[i][0] + r2_avg[i][0]*r2_avg[i][0] ) + dt*(r2[i][0]*r2[i][0]) ) / stats_time;
      r2_rms[i][1] = ( old_stats_time*( r2_rms[i][1]*r2_rms[i][1] + r2_avg[i][1]*r2_avg[i][1] ) + dt*(r2[i][1]*r2[i][1]) ) / stats_time;
      r2_rms[i][2] = ( old_stats_time*( r2_rms[i][2]*r2_rms[i][2] + r2_avg[i][2]*r2_avg[i][2] ) + dt*(r2[i][2]*r2[i][2]) ) / stats_time;
      // put sum( ui*uj * dt )_new / sum( dt ) in r2_rey(i)...
      r2_rey[i][0] = ( old_stats_time*( r2_rey[i][0] + r2_avg[i][1]*r2_avg[i][2] ) + dt*(r2[i][1]*r2[i][2]) ) / stats_time;
      r2_rey[i][1] = ( old_stats_time*( r2_rey[i][1] + r2_avg[i][2]*r2_avg[i][0] ) + dt*(r2[i][2]*r2[i][0]) ) / stats_time;
      r2_rey[i][2] = ( old_stats_time*( r2_rey[i][2] + r2_avg[i][0]*r2_avg[i][1] ) + dt*(r2[i][0]*r2[i][1]) ) / stats_time;
      // put new means in r2_avg...
      r2_avg[i][0] = ( old_stats_time*r2_avg[i][0] + dt*(r2[i][0]) ) / stats_time;
      r2_avg[i][1] = ( old_stats_time*r2_avg[i][1] + dt*(r2[i][1]) ) / stats_time;
      r2_avg[i][2] = ( old_stats_time*r2_avg[i][2] + dt*(r2[i][2]) ) / stats_time;
    }

    // apply additional averaging...
    if (average != NULL) {
      average->apply( r2_avg );
      average->apply( r2_rms );
      average->apply( r2_rey );
    }

    for (int i = 0; i < n; ++i) {
      // complete rms...
      r2_rms[i][0] = sqrt( fabs( r2_rms[i][0] - r2_avg[i][0]*r2_avg[i][0] ) );
      r2_rms[i][1] = sqrt( fabs( r2_rms[i][1] - r2_avg[i][1]*r2_avg[i][1] ) );
      r2_rms[i][2] = sqrt( fabs( r2_rms[i][2] - r2_avg[i][2]*r2_avg[i][2] ) );
      // complete rey...
      r2_rey[i][0] = r2_rey[i][0] - r2_avg[i][1]*r2_avg[i][2];
      r2_rey[i][1] = r2_rey[i][1] - r2_avg[i][2]*r2_avg[i][0];
      r2_rey[i][2] = r2_rey[i][2] - r2_avg[i][0]*r2_avg[i][1];
    }

  }
      


};  
  
class HITCf4 : public Charles {
  
private:

public:
  
  // Constructor that takes input file...
  HITCf4(char *name) : Charles(name) {
    
    if (mpi_rank == 0)
      cout << "HITCf4()" << endl;
    
  }
  
  void initialHook() {
    
    if (step == 0) {

      if (mpi_rank == 0)
	cout << " > HITCf4: initialHook()" << endl;
      
      FOR_ICV taylorGreenIC(rho[icv],rhou[icv],rhoE[icv],
			    x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
      
    }
    
  }
    
  void taylorGreenIC(double &rho,double *rhou,double &rhoE,
		     const double x,const double y,const double z) {
    rho = rho_ref;
    rhou[0] = rho_ref*(sin(x)*cos(y)*cos(z));
    rhou[1] = rho_ref*(-cos(x)*sin(y)*cos(z));
    rhou[2] = 0.0;
    const double p = p_ref + ( (cos(2.0*z) + 2.0)*(cos(2.0*x) + cos(2.0*y)) - 2.0 )/16.0;
    rhoE = p/(gamma - 1.0) + 
      0.5*( rhou[0]*rhou[0] + 
	    rhou[1]*rhou[1] + 
	    rhou[2]*rhou[2] )/rho;
  }

  
  /*
  void calcSourceHook(double * rho_rhs,double (*rhou_rhs)[3],double * rhoE_rhs,
		      const double *rho,const double (*rhou)[3],const double * rhoE,
		      const double rk_time,const int rk_step) {
    
    // unit forcing in the x-direction...

    
    FOR_ICV {
      rho_rhs[icv]     = 0.0;
      rhou_rhs[icv][0] = rho[icv]*cv_volume[icv];
      rhou_rhs[icv][1] = 0.0;
      rhou_rhs[icv][2] = 0.0;
      rhoE_rhs[icv]    = rhou[icv][0]*cv_volume[icv];
    }
    
  }
  */


};  
  


class ShockTubeCf4 : public Charles {

private:
  
public:
  
  // Constructor that takes input file...
  ShockTubeCf4(char *name) : Charles(name) {
    
    if (mpi_rank == 0)
      cout << "ShockTubeCf4()" << endl;
    
  }
  
  void initialHook() {
    
    if (mpi_rank == 0)
      cout << " > ShockTubeCf4: initialHook()" << endl;
    
    double pL = 1.0;
    double rhoL = 1.0;
 
    double pR = 0.1;
    double rhoR = 0.125;
 
    FOR_ICV {
      if (x_cv[icv][0] < 0.0) {
	rho[icv] = rhoL;
	rhou[icv][0] = rhou[icv][1] = rhou[icv][2] = 0.0;
	rhoE[icv] = pL/(gamma-1.0);
      }
      else {
	rho[icv] = rhoR;
	rhou[icv][0] = rhou[icv][1] = rhou[icv][2] = 0.0;
	rhoE[icv] = pR/(gamma-1.0);
      }
    }
    
  }

};



class JetLinearCf4 : public Charles {
  
private:
  
public:
  
  // Constructor that takes input file...
  JetLinearCf4 (char *name) : Charles(name) {
    
    if (mpi_rank == 0)
      cout << "JetLinearCf4()" << endl;
  }

  // hooks are the virtual functions in Charles...

  void initialHook() {

    //if (step == 0) {

      if (mpi_rank == 0)
	cout << "Applying initial condition" << endl;
      
      FOR_ICV {
	rho[icv] = rho_ref;
	rhou[icv][0] = 0.0;
	rhou[icv][1] = 0.0;
	rhou[icv][2] = 0.0;
	rhoE[icv] = p_ref/(gamma-1.0);
      }
      
      //}

    // =========================================
    // and set some operators to first order...
    // =========================================
      /*
    FOR_ICV cv_flag[icv] = 0;
    double xp = -0.4;
    FOR_ICV {
      // sphere near nozzle...
      double dx = x_cv[icv][0] - xp;
      double dy = x_cv[icv][1]; 
      double dz = x_cv[icv][2];
      // sphere near nozzle...
      if ( sqrt(dx*dx + dy*dy + dz*dz) < 0.25 )
	cv_flag[icv] = 1;
      // outside x-cylinder, r = 0.5... 
      if ( sqrt(dy*dy + dz*dz) > 0.5 )
	cv_flag[icv] = 1;
      // outlet region...
      if ( x_cv[icv][0] > 1.75 )
	cv_flag[icv] = 1;
      
    }
    firstOrderOperatorsFlaggedCvs();
      */
    
    // ============================
    // as a test, set fa_alpha to 1 everywhere - fully upwinded
    // ============================
    //FOR_IFA fa_alpha[ifa] = 1.0;

    // ============================
    // also reset the time...
    // ============================
    time = 0.0;
    
  }

  // apply source term...
  
  void calcSourceHook(double * rho_rhs,double (*rhou_rhs)[3],double * rhoE_rhs,
		      const double *rho,const double (*rhou)[3],const double * rhoE,
		      const double rk_time,const int rk_step) {
    
    double xp[3] = { 0.635, 0.0, 0.0 };
    double freq = 13385.8;
    double angular_freq = 2 * M_PI * freq;
    // speed of sound (squared) is...
    double c2 = gamma*p_ref/rho_ref;
    double one_over_sigma2 = 10.0*freq*freq/c2;
    double strength = 10.0*sin(angular_freq*rk_time);
    
    if ((rk_step == 1)&&(mpi_rank == 0))
      cout << "time, strength: " << rk_time << " " << strength << endl;
    
    FOR_ICV {
      
      // calc the radius of this cv centroid...
      double dr2 = 0.0;
      FOR_I3 {
        double dx = x_cv[icv][i] - xp[i];
        dr2 += dx*dx;
      }
      
      double force = strength*exp(-dr2*one_over_sigma2)*cv_volume[icv];
      
      rho_rhs[icv]            = rho[icv]*force;
      FOR_I3 rhou_rhs[icv][i] = rhou[icv][i]*force;
      rhoE_rhs[icv]           = rhoE[icv]*force;

      /*
	double H0 =  rhoE[icv] + (gamma-1.0)*(rhoE[icv]-0.5*( 
	rhou[icv][0]*rhou[icv][0] +
	rhou[icv][1]*rhou[icv][1] +
	rhou[icv][2]*rhou[icv][2] )/rho[icv]);
	rhoE_rhs[icv] = rho_rhs[icv] * H0;
      */
      
    }
    
  }

};


    
int main(int argc, char * argv[]) {

  MPI_Init(&argc, &argv);

  try {
      
    // mpi is in the namespace "MpiStuff" that contains 
    // commonly used things like 
    // mpi_rank, mpi_size, mpi_comm, etc...

    initMpiStuff();
    
    // Base...
    //Charles cf4("cf4.in");
    
    // Joe's SFD class...
    //SfdCf4 cf4("cf4_sfd.in");
    
    // Shocktube class
    //ShockTubeCf4 cf4("cf4_shocktube.in");
    
    // ShuOsher...
    //ShuOsherCf4 cf4("cf4_ShuOsher.in");

    // Euler vortex...
    //EulerCf4 cf4("cf4_euler.in");

    // Taylor-Green
    //TaylorGreenCf4 cf4("cf4_tg.in");
    
    // 1D testing...
    //SmallCf4 cf4("cf4_small.in");
    
    // Woodward-Collela bump... 
    WCCf4 cf4("cf4_wc.in");
    
    // Laminar or Turbulent Channel...
    //ChannelCf4 cf4("cf4_channel.in");

    // Laminar or Turbulent Channel...
    //HITCf4 cf4("cf4_hit.in");

    // interp version...
    //InterpCf4 cf4("cf4.in");

    // include regions of first order for certain meshes...
    //ModCf4 cf4("cf4.in");

    // do something - also Honeywell...
    //MiscCf4 cf4("cf4.in");
    
    // G's blayer stuff...
    //BlayerCf4 cf4("cf4_bl.in");
   
    // Sun model...
    //SunCf4 cf4("cf4_sun.in");

    // for Simon...
    //JetLinearCf4 cf4("cf4_jet_linear.in");
    JetCf4 cf4("cf4.in");
    //UtrcCf4 cf4("cf4.in");

    //RinglebCf4 cf4("cf4.in");

    cf4.init();
    cf4.run();
    
  }
  catch (int e) {
    cerr << "Exception: " << e << endl;
    MPI_Finalize();
    return(-1);
  }
  catch(...) {
    cerr << "unhandled Exception.\n" << endl;
    MPI_Finalize();
    return(-1);
  }
    
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return(0);
  
}



