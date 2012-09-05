#ifndef UGPWITHCVCOMPFLOW_H
#define UGPWITHCVCOMPFLOW_H

#include "MiscUtils.h"
#include "tc_vec3d.h"

#include "UgpWithCv.h"
#include "CdpFilter.h"
#include "MshFilter.h"
#include "Param.h"

#include "AtmosphericCond.h"

#include "HypreSolver.h"
#include "PetscSolver.h"

#include "Logging.h"
using namespace logging;


//#define temp_reconstruction   // determine whether temperature or pressure reconstruction is used at faces
//#define alpha_limiter         // determine if rho*Phi is limited based on rho or rho*Phi


/**
 * Class contains the main routines for unstructured Euler for now
 *
 * \author Frank Ham, Rene Pecnik and Vincent Terrapon
 * \date March, 2008
 */


class UgpWithCvCompFlow : public UgpWithCv, public ParamMap
{
public:   // constructors/destructors

  /**
   * standard constructor needed for Models which virtual inherit from UgpWithCvCompFlow
   */
  UgpWithCvCompFlow() {    init();  }

  /**
   * constructor, pass ParamMap
   */
  UgpWithCvCompFlow(ParamMap &p) : ParamMap(p) {    init();  }

  /**
   * constructor, pass name of joe's input file
   */
  UgpWithCvCompFlow(char *name) : ParamMap(name) {    init();  }


  /**
   * registers scalars and vectors for main equations (Navier-Stokes)
   * and set main parameters (solver: nsteps; NS-eq: gravity...)
   */
  void init()
  {
    Param *param;
    if (getParam(param, "LOGGING_LEVEL"))
    {
      string name = param->getString();
      if      (name == "EVERYTHING") setLoggingLevel(EVERYTHING);
      else if (name == "DEBUG_HI"  ) setLoggingLevel(DEBUG_HI  );
      else if (name == "DEBUG_LO"  ) setLoggingLevel(DEBUG_LO  );
      else if (name == "INFO_HI"   ) setLoggingLevel(INFO_HI   );
      else if (name == "INFO_LO"   ) setLoggingLevel(INFO_LO   );
      else if (name == "WARNING"   ) setLoggingLevel(WARNING   );
      else if (name == "ERROR"     ) setLoggingLevel(ERROR     );
      else if (name == "CRITICAL"  ) setLoggingLevel(CRITICAL  );
      else if (name == "FATAL"     ) setLoggingLevel(FATAL     );
      else if (name == "SILENT"    ) setLoggingLevel(SILENT    );
      else {
        lout(WARNING) << "Warning: Unknown logging level, using default.\n";
      }
    }


    if (mpi_rank == 0)    cout << "UgpWithCvCompFlow()"<< endl;

    // ----------------------------------------------------------------------------------------
    // register data...
    // ----------------------------------------------------------------------------------------

	// multigrid sensors
    MGSensor  = NULL;         registerScalar(MGSensor, "MG_SENSOR", CV_DATA);
	  
    // conservative variables
    rho   = NULL;         registerScalar(rho,  "RHO",  CV_DATA);
    rhou  = NULL;         registerVector(rhou, "RHOU", CV_DATA);
    rhoE  = NULL;         registerScalar(rhoE, "RHOE", CV_DATA);

    // interpolated conservative variables
    rho_interp   = NULL;         registerScalar(rho_interp,  "RHO_INTERP",  CV_DATA);
    rhou_interp  = NULL;         registerVector(rhou_interp, "RHOU_INTERP", CV_DATA);
    rhoE_interp  = NULL;         registerScalar(rhoE_interp, "RHOE_INTERP", CV_DATA);
	  
    // primitive variables
    vel      = NULL;      registerVector(vel,      "vel",      CV_DATA);
    press    = NULL;      registerScalar(press,    "press",    CV_DATA);
    temp     = NULL;      registerScalar(temp,     "temp",     CV_DATA);
    enthalpy = NULL;      registerScalar(enthalpy, "enthalpy", CV_DATA);

    // material properties
    gamma    = NULL;      registerScalar(gamma,    "gamma",    CV_DATA);
    RoM      = NULL;      registerScalar(RoM,      "RoM",      CV_DATA);
    sos      = NULL;      registerScalar(sos,      "sos",      CV_DATA);
    
    // diffusion coefficient for momentum and energy at cell faces
    mul_fa = NULL;
    lamOcp_fa = NULL;
    
    massFlux_fa = NULL;

    // turbulence 
    mut_fa = NULL;          // defined at faces, allocated in initializeFromRestartFile
    kine   = NULL;          // registerScalar(kine, "kine", CV_DATA);

    // just to take a look at the residual field of energy
    residField  = NULL;     registerScalar(residField,  "residField",  CV_DATA);
    residField2 = NULL;     registerScalar(residField2, "residField2", CV_DATA); 

    local_dt = NULL;        registerScalar(local_dt, "local_dt", CV_DATA);


    grad_rho = NULL;        // register only for second order euler flux or if viscous
    grad_u = NULL;          // register only for second order euler flux or if viscous
    grad_p = NULL;          // register only for second order euler flux, if pressure reconstruction is used
    grad_temp = NULL;       // register only for second order euler flux, if temperature reconstruction is used
    grad_enthalpy = NULL;   // register only if viscous
    
    alpha_rho = NULL;       // register only for second order euler flux to limit scalars

    strMag = NULL;
    vortMag = NULL;
    diverg = NULL;


    // ----------------------------------------------------------------------------------------
    // register flow parameters
    // ----------------------------------------------------------------------------------------

    gravity[0] = getDoubleParam("GX", "0.0");
    gravity[1] = getDoubleParam("GY", "0.0");
    gravity[2] = getDoubleParam("GZ", "0.0");


    turbModel = NONE;

    viscMode      = getStringParam("MU_MODE", "POWERLAW");
    mu_ref        = getDoubleParam("MU_REF", "0.0");           // NOTE: if MU_MODE == SUTHERLAND then set MU_REF = 1.716e-5
    mu_power_law  = getDoubleParam("MU_POWER_LAW", "0.76");    // note: specify 0 here for constant visc = mu_ref
    SL_Tref       = getDoubleParam("SL_TREF", "273.15");       // Sutherland law, Tref should be 273.15
    SL_Sref       = getDoubleParam("SL_SREF", "110.4");        // Sutherland law, S value should be 110.4

    // if HEIGHT is specified in the input file T_ref, p_ref, rho_ref and mu_ref are recalculated using getAthmospericConditions
    if (checkParam("ALTITUDE"))
    {
      if (!getAthmospericConditions(T_ref, p_ref, rho_ref, getDoubleParam("ALTITUDE"), 9.81, 1.4, 287.0))
      {
        cerr << "ERROR in getAthmospericConditions occurred" << endl;
        throw(-1);
      }

      viscMode = "SUTHERLAND";    // change to Sutherland law for viscosity if height is specified 
      mu_ref = 1.716e-5;          // hard coded for reference viscosity at conditions given below
      SL_Tref = 273.15;           // hard coded for reference temperature in Sutherland law
      SL_Sref = 110.4;            // hard coded for constant in Sutherland law
    }
    else
    {
      T_ref    = getDoubleParam("T_REF");
      p_ref    = getDoubleParam("P_REF");
      rho_ref  = getDoubleParam("RHO_REF");
    }

    // ideal gas law is P / rho = R * T ...
    GAMMA  = getDoubleParam("GAMMA", "1.4");
    R_gas  = p_ref/rho_ref/T_ref;
    Pr     = getDoubleParam("PR", "0.72");
    PrTurb = getDoubleParam("PRTURB", "0.9");

    // ----------------------------------------------------------------------------------------
    // register solver parameters
    // ----------------------------------------------------------------------------------------

    // check if second order for Navier-Stokes
    if (checkParam("SPATIAL_SECOND_ORDER"))       sndOrder = true;
    else                                          sndOrder = false;

    if (checkParam("SCALAR_FIRST_ORDER"))         stOrderScalar = true;
    else                                          stOrderScalar = false;

    // approximate Riemann solver
    // for coupled or semi-coupled simulations only HLLC, HAENEL, AUSMD, AUSMV, or AUSMDV
    SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");  // "HLLC", "ROE", "AUSM", "AUSMD", "AUSMV", "AUSMDV", "HAENEL", "JST"
    if (SpaceIntName == "AUSMD")
    {
      SpaceIntName = "AUSMDV";
      AUSM_Type = "D";
    }
    else if (SpaceIntName == "AUSMV")
    {
      SpaceIntName = "AUSMDV";
      AUSM_Type = "V";
    }
    else
      AUSM_Type = "DV";

    ShockFixFlagCell = NULL;                  // array containing the flagged cells

    if (SpaceIntName == "AUSMDV")
    {
      Shock_Fix_Type = getStringParam("SHOCK_FIX_TYPE", "HAENEL");         // "NONE", "HAENEL", "ALPHA_LR"
      Shock_Fix_Flag = getStringParam("SHOCK_FIX_FLAG", "CELL_FLAG_HEX");  // "CELL_FLAG_HEX", "CELL_FLAG_TET"
      if (Shock_Fix_Flag == "CELL_FLAG_HEX")
        Shock_Fix_Flag_Sum = 1.1;             // both cells on each side of the interface must be flagged to use the shock fix
      else
        Shock_Fix_Flag_Sum = 0.9;             // only one cell on either side of the interface must be flagged to use the shock fix
      Entropy_Fix = getStringParam("ENTROPY_FIX", "EXPANSION_FAN");  // "NONE", "EXPANSION_FAN"
    }
    else
    {
      Entropy_Fix    = "NONE";
      Shock_Fix_Type = "NONE";
      Shock_Fix_Flag = "NONE";
      Shock_Fix_Flag_Sum = 0.0;
    }
    if (Shock_Fix_Type != "NONE")
      registerScalar(ShockFixFlagCell, "ShockFixFlagCell", CV_DATA);

    // register type of gradient reconstruction and limiters
    if (!checkParam("GRAD_RECONSTRUCTION"))
    {
      ParamMap::add("GRAD_RECONSTRUCTION  SCHEME=GREENGAUSS  LIMITER=YES  EPS=0.001");    // add default values
      if (mpi_rank == 0)
        cout << "WARNING: added keyword \"GRAD_RECONSTRUCTION  SCHEME=GREENGAUSS  LIMITER=YES  EPS=0.001\"" <<
                " to parameter map!" << endl;
    }
    string gradReconstr = getParam("GRAD_RECONSTRUCTION")->getString("SCHEME");

    if      (gradReconstr == "GREENGAUSS")       gradreconstruction = GRAD_GG;
    else if (gradReconstr == "LEASTSQUARES")     gradreconstruction = GRAD_LS;
    else if (gradReconstr == "SDWLS")            gradreconstruction = GRAD_SDWLS;
    else
    {
      if (mpi_rank == 0)
        cerr << "ERROR: gradient reconstruction type "
             << "\"" << gradReconstr << "\"" << " not recognized, choose GRAD_RECONSTRUCTION = \"GREENGAUSS\", \"LEASTSQUARE\" or \"SDWLS\"" << endl;
      throw(-1);
    }

    string limiterstr;
    if ((gradreconstruction == GRAD_GG) || (gradreconstruction == GRAD_LS))
    {
      limiterstr   = getParam("GRAD_RECONSTRUCTION")->getString("LIMITER");

      if      ((limiterstr == "NOLIMITER") || (limiterstr == "NO"))              limiterNavierS = NOLIMITER;
      else if ((limiterstr == "BARTH_JESPERSEN_MOD") || (limiterstr == "YES"))   limiterNavierS = BARTH_JESPERSEN_MOD;
      else 
      {
        if (mpi_rank == 0)
          cerr << "ERROR: Grandient limiter (" << limiterstr <<
              ") not recognized, choose GRAD_RECONSTRUCTION = SCHEME=<YES, NO, NOLIMITER or BARTH_JESPERSEN_MOD>" << endl;
        throw(-1);
      }

      epsilonSDWLS = getParam("GRAD_RECONSTRUCTION")->getDouble("EPS");
    }
    else if (gradreconstruction == GRAD_SDWLS)
    {
      limiterstr = "VANLEER";
      limiterNavierS = 0;
      epsilonSDWLS = getParam("GRAD_RECONSTRUCTION")->getDouble("EPS");
    }

    // register CFL condition
    
    if (checkParam("CFL"))
    {
      timeStepMode = "CFL";
      cfl = getDoubleParam("CFL");
    }
    else if (checkParam("CFL_MIN"))
    {
      timeStepMode = "CFL_MIN";
      cfl = getDoubleParam("CFL_MIN");
    }
    else if (checkParam("DT"))
    {
      timeStepMode = "DT";
      cfl = 0.0;
    }
    else
    {
      cout << "ERROR: DT, CFL or CFL_MIN not specified" << endl;
      throw(-1);
    }


    nsteps = getIntParam("NSTEPS", "100");
    registerValue(step,"STEP");
    registerValue(time,"TIME");

    runtime_flag = 0;
    wtime = MPI_Wtime();
    if (checkParam("RUNTIME"))
    {
      runtime_flag = 1;
      if (mpi_rank==0)
      {
        runtime = getDoubleParam("RUNTIME");
        cout<<" > RUNTIME [hours]: "<<runtime<<endl;
        // convert to the time in seconds when we have to stop...
        runtime = runtime*3600.0+wtime;
      }
    }

    write_restart  = getIntParam("WRITE_RESTART", "-1");     // default is no writing
    check_interval = getIntParam("CHECK_INTERVAL", "10");

    resid_energ_th = getDoubleParam("RESID_ENERG_TH", 0.0);



    // ----------------------------------------------------------------------------------------
    // set Navier-Stokes solver
    // ----------------------------------------------------------------------------------------

    Param * p;
    string tIntName = getStringParam("TIME_INTEGRATION", "BACKWARD_EULER");
    string sIntName = getStringParam("SPACE_INTEGRATION", "HLLC");


    // ----------------------------------------------------------------------------------------
    // pick linear solver for Navier-Stokes
    // ----------------------------------------------------------------------------------------

    petscSolver = NULL;
    nbocv_v_global = NULL;

    if (getParam(p,"LINEAR_SOLVER_NS"))
    {
      string name = p->getString();

      if (name == "PETSC_GMRES")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: PETSC_GMRES" << endl;
        linearSolverNS = PETSC_GMRES;
      }
      else if (name == "BCGSTAB")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: BCGSTAB" << endl;
        linearSolverNS = BCGSTAB;
      }
			else if (name == "LUSGS")
			{
				if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: LUSGS" << endl;
				linearSolverNS = LUSGS;
			}
			else if (name == "BCGSTAB_LINELET")
			{
				if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: BCGSTAB_LINELET" << endl;
				linearSolverNS = BCGSTAB_LINELET;
			}
      else
      {
        if (mpi_rank == 0)        cerr << "Error: unrecognized LINEAR_SOLVER_NS: " << name << endl;
        throw(-1);
      }
    }
    else
    {
      if (mpi_rank == 0)          cout << "LINEAR_SOLVER_NS: BCGSTAB (default)" << endl;
      linearSolverNS = BCGSTAB;
    }


    // ----------------------------------------------------------------------------------------
    // pick linear solver for scalar
    // ----------------------------------------------------------------------------------------
    petscSolverScalars = NULL;

    if (getParam(p,"LINEAR_SOLVER_SCALARS"))
    {
      string name = p->getString();

      if (name == "PETSC_GMRES")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: PETSC_GMRES" << endl;
        linearSolverScalars = PETSC_GMRES;
      }
      else if (name == "BCGSTAB")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: BCGSTAB" << endl;
        linearSolverScalars = BCGSTAB;
      }
			else if (name == "LUSGS")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: LUSGS" << endl;
        linearSolverScalars = LUSGS;
      }
			else if (name == "BCGSTAB_LINELET")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: BCGSTAB_LINELET" << endl;
        linearSolverScalars = BCGSTAB_LINELET;
      }
      else
      {
        if (mpi_rank == 0)        cerr << "Error: unrecognized LINEAR_SOLVER_NS: " << name << endl;
        throw(-1);
      }
    }
    else
    {
      if (linearSolverNS == PETSC_GMRES)
      {
        linearSolverScalars = PETSC_GMRES;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: PETSC_GMRES (default)" << endl;
      }
			else if (linearSolverNS == LUSGS)
      {
        linearSolverScalars = LUSGS;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: LUSGS (default)" << endl;
      }
			else if (linearSolverNS == BCGSTAB_LINELET)
      {
        linearSolverScalars = BCGSTAB_LINELET;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: BCGSTAB_LINELET (default)" << endl;
      }
      else
      {
        linearSolverScalars = BCGSTAB;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: BCGSTAB (default)" << endl;
      }
    }

    initial_flowfield_output = getStringParam("INITIAL_FLOWFIELD_OUTPUT", "YES");

    // ----------------------------------------------------------------------------------------
    // write some parameters on the screen
    // ----------------------------------------------------------------------------------------

    if (mpi_rank == 0)
    {
      cout << "Gas properties         " << endl;
      cout << "    GAMMA            : " << GAMMA << endl;
      cout << "    R_GAS            : " << R_gas << endl;
      cout << "    P_REF            : " << p_ref << endl;
      cout << "    RHO_REF          : " << rho_ref << endl;
      cout << "    T_REF            : " << T_ref << endl;
      cout << "    SOS_REF          : " << sqrt(GAMMA*R_gas*T_ref) << endl;
      cout << "Material properties    " << endl;
      cout << "    MU_MODE          : " << viscMode << endl;
      cout << "      mu_ref         : " << mu_ref << endl;
      cout << "      mu_power_law   : " << mu_power_law << endl;
      cout << "      SL_Tref        : " << SL_Tref << endl;
      cout << "      SL_Sref        : " << SL_Sref << endl;
      cout << "    Pr               : " << Pr << endl;
      cout << "    PrTurb           : " << PrTurb << endl;
      cout << "Body Force           : " << endl;
      cout << "    gravity          : " << gravity[0] << " " << gravity[1] << " " << gravity[2] << endl;
      cout << "Solver settings        " << endl;
      cout << "    second order     : " << sndOrder << endl;
	    cout << "    1st order scalars: " << stOrderScalar << endl;
      cout << "    gradient recon   : " << gradReconstr << ", value = " << gradreconstruction << endl;
      cout << "    gradient limiter : " << limiterstr << ", value = " << limiterNavierS << endl;
	    cout << "    viscous limiter  : " << getStringParam("VISC_LIMITER","YES") << endl;
      if (gradreconstruction == GRAD_SDWLS) 
        cout << "                     : " << "epsilonSDWLS = " << epsilonSDWLS << endl;
      cout << "TIME_INTEGRATION     : " << tIntName << endl;
      cout << "    nsteps           : " << nsteps << endl;
      cout << "    timeStepMode     : " << timeStepMode << endl;
      cout << "    timeStepMethod   : " << getStringParam("TIME_STEP","MAX_BASED") << endl;
      cout << "SPACE_INTEGRATION    : " << sIntName << endl;
      cout << "    SymBCMethod      : " << getStringParam("SYMM_BC","ORIGINAL_NUMERIC") << endl;
      cout << "    Shock Fix        : " << Shock_Fix_Type << endl;
      cout << "    Shock Fig Flag   : " << Shock_Fix_Flag << endl;
      cout << "    Entropy Fix      : " << Entropy_Fix << endl;
      cout << "COMBUSTION_REGIME    : " << getStringParam("COMBUSTION_REGIME", "HOT") << endl;
      cout << "    SmoothSource     : " << getIntParam("NSMOOTH_MEANSOURCE","1") << endl;
      cout << "CHEMICAL_TABLE AUQ   : " << getStringParam("UQ_TABLE", "NO") << endl;
      cout << "    Source Term      : " << getDoubleParam("SOURCE_PERT_FACTOR","0.0") << endl;
      cout << "    Temperature Term : " << getDoubleParam("TEMP_PERT_FACTOR","0.0") << endl;

      cout << "--------------------------------------------" << endl << endl;
    }


    // ----------------------------------------------------------------------------------------
    // parse any WRITE_DATA parameters...
    // ----------------------------------------------------------------------------------------

    // init writedata
    initWriteData(this);

    // init probes
    initProbes(this);
  }

  virtual void initializeFromRestartFile(const string &name)
  {
    // ----------------------------------------------------------------------------------------
    // read file
    // ----------------------------------------------------------------------------------------
    
    UgpWithCv::initializeFromRestartFile(name);
    
    
    // ----------------------------------------------------------------------------------------
    // boundary faces, allocate memory, TODO: should be incorporated in to the associated arrays ???
    // ----------------------------------------------------------------------------------------

    rho_bfa = new double[nfa_b];
    T_bfa   = new double[nfa_b];
    vel_bfa = new double[nfa_b][3];
    p_bfa   = new double[nfa_b];
    h_bfa   = new double[nfa_b];
    gam_bfa = new double[nfa_b];
    RoM_bfa = new double[nfa_b];

    // ----------------------------------------------------------------------------------------
    // init memory for face-based data
    // ----------------------------------------------------------------------------------------
    massFlux_fa = new double[nfa];
    mut_fa    = new double[nfa];
    mul_fa    = new double[nfa];
    lamOcp_fa = new double[nfa];
    
    // ----------------------------------------------------------------------------------------
    // init memory for scalar diffusion coeff
    // ----------------------------------------------------------------------------------------
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      data->phi_bfa  = new double[nfa_b];
      data->grad_phi = new double[ncv_g][3];
      if (data->reconstruction == "CONSERVATIVE")
      {
        data->rhophi      = new double[ncv_g];
        data->rhophi_bfa  = new double[nfa_b];
        data->grad_rhophi = new double[ncv_g][3];        
      }
      data->diff     = new double[nfa];
    }
    
    // ----------------------------------------------------------------------------------------
    // gradients
    // ----------------------------------------------------------------------------------------
    if (sndOrder == true)
    {
      grad_rho = new double[ncv_g][3];
      registerScalar(alpha_rho,  "alpha_rho",  CV_DATA);
#ifdef temp_reconstruction
      grad_temp = new double[ncv_g][3];
#else
      grad_p = new double[ncv_g][3];
#endif
//      grad_gam = new double[ncv_g][3];
    }

    grad_u = new double[ncv_g][3][3];             // allocate always

    if ((mu_ref > 0.0) || (sndOrder == true))     // for now, will change
      grad_enthalpy = new double[ncv_g][3];
  }

  virtual ~UgpWithCvCompFlow()
  {
    if (rho_bfa != NULL)         {delete [] rho_bfa;        rho_bfa = NULL;}
    if (T_bfa   != NULL)         {delete [] T_bfa;          T_bfa   = NULL;}
    if (vel_bfa != NULL)         {delete [] vel_bfa;        vel_bfa = NULL;}
    if (p_bfa   != NULL)         {delete [] p_bfa;          p_bfa   = NULL;}
    if (h_bfa   != NULL)         {delete [] h_bfa;          h_bfa   = NULL;}
    if (gam_bfa != NULL)         {delete [] gam_bfa;        gam_bfa = NULL;}
    if (RoM_bfa != NULL)         {delete [] RoM_bfa;        RoM_bfa = NULL;}
    
    
    if (massFlux_fa != NULL)       {delete [] massFlux_fa;       massFlux_fa = NULL;}
    if (mut_fa      != NULL)       {delete [] mut_fa;            mut_fa      = NULL;}
    if (mul_fa      != NULL)       {delete [] mul_fa;            mul_fa      = NULL;}
    if (lamOcp_fa   != NULL)       {delete [] lamOcp_fa;         lamOcp_fa   = NULL;}

    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      if (data->phi_bfa     != NULL)     {delete [] data->phi_bfa;            data->phi_bfa     = NULL;}
      if (data->grad_phi    != NULL)     {delete [] data->grad_phi;           data->grad_phi    = NULL;}
      if (data->rhophi      != NULL)     {delete [] data->rhophi;             data->rhophi      = NULL;}
      if (data->rhophi_bfa  != NULL)     {delete [] data->rhophi_bfa;         data->rhophi_bfa  = NULL;}
      if (data->grad_rhophi != NULL)     {delete [] data->grad_rhophi;        data->grad_rhophi = NULL;}
      if (data->diff        != NULL)     {delete [] data->diff;               data->diff        = NULL;}
    }

    if (grad_rho      != NULL)          {delete [] grad_rho;             grad_rho      = NULL;}
#ifdef temp_reconstruction
    if (grad_temp     != NULL)          {delete [] grad_temp;            grad_temp     = NULL;}
#else
    if (grad_p        != NULL)          {delete [] grad_p;               grad_p        = NULL;}
#endif
    if (grad_u        != NULL)          {delete [] grad_u;               grad_u        = NULL;}
    if (grad_enthalpy != NULL)          {delete [] grad_enthalpy;        grad_enthalpy = NULL;}
//    if (grad_gam      != NULL)          {delete [] grad_gam;             grad_gam      = NULL;}

    if (ShockFixFlagCell != NULL)       {delete [] ShockFixFlagCell;     ShockFixFlagCell = NULL;}
  }

public:   // member variables

  double  *rho;             ///< density \f$\rho\f$
  double (*rhou)[3];        ///< momentum \f$\rho u\f$
  double  *rhoE;            ///< energy \f$\rho E\f$
	
  double  *rho_interp;             ///< density \f$\rho\f$
  double (*rhou_interp)[3];        ///< momentum \f$\rho u\f$
  double  *rhoE_interp;            ///< energy \f$\rho E\f$

  double (*vel)[3];         ///< velocity
  double *press;            ///< pressure
  double *temp;             ///< temperature
  double *enthalpy;         ///< chemical + sensible enthalpy (no kinetic energy!)

  double *RoM;              ///< gas constant at cell center
  double *gamma;            ///< ratio of specific heat at cell center
  double *sos;              ///< speed of sound at cell center
  
  double *mul_fa;           ///< laminar viscosity at cell faces
  double *lamOcp_fa;        ///< heat conductivity at cell faces
  
  double *mut_fa;           ///< turbulent viscosity at cell faces
  double *kine;     
  //  double *kine, *eps;         ///< turbulent kinetic energy and turbulent dissipation 

  double *massFlux_fa;      ///< mass flux for scalar solver

  double *residField;
  double *residField2;
  double resid_energ_th;    ///< energy residual treshold

  double *strMag;           ///< magnitude strain rate
  double *vortMag;          ///< magnitude vorticity
  double *diverg;           ///< divergence

  // ----------------------------------------------------------------------------------------
  // boundary faces, allocate memory, TODO: should be incorporated in to the associated arrays ???
  // ----------------------------------------------------------------------------------------
  double *rho_bfa, *T_bfa, (*vel_bfa)[3], *p_bfa, *h_bfa, *gam_bfa, *RoM_bfa;

  double (*grad_enthalpy)[3]; ///< enthalpy gradient, registered only if 2nd order euler flux calc turned on
  double (*grad_rho)[3];      ///< density gradient, registered only if 2nd order euler flux calc turned on
  double (*grad_u)[3][3];     ///< velocity gradient tensor, registered only 2nd order euler flux calc turned on or viscous calc
  double (*grad_p)[3];        ///< pressure gradient, registered only if 2nd order euler flux calc turned on and temp_reconstruction NOT defined
  double (*grad_temp)[3];     ///< temperature gradient, registered only if 2nd order euler flux calc turned on and temp_reconstruction defined
//  double (*grad_gam)[3];
  
  double *alpha_rho;          ///< limiter factor for gradients of rho, used to limit scalars in computation of Euler flux

  double p_ref, rho_ref, T_ref, R_gas;
  double GAMMA, Pr, PrTurb;

  string viscMode;
  double mu_ref;
  double mu_power_law;
  double SL_Tref, SL_Sref;    ///< Sutherland law, Tref and S value, should be 273.15 and 110.4

  double gravity[3];

  double cfl, *local_dt;

  string timeStepMode;
  int nsteps, step;
  double time;

  double runtime;    ///< user specified runtime in hours
  int runtime_flag;  ///< if flag==1 then use runtime to determine the job cancellation
  double wtime;

  int write_restart, check_interval;

  int gradreconstruction;
  bool sndOrder;
  bool stOrderScalar;
  int limiterNavierS;
  double epsilonSDWLS;
  
  int turbModel;
  enum TurbModel{NONE, SA, KEPS, KOM, KOMSST, V2F};


  // Approximate Riemann solver
  string SpaceIntName;         ///< Approximate Riemann solver: "HLLC", "ROE", "AUSM", "AUSMD", "AUSMV", "AUSMDV", "HAENEL", or "JST"
  string AUSM_Type;            ///< Type of AUSMDV: "DV" -> AUSMDV, "D" -> AUSMD, "V" -> AUSMV
  string Shock_Fix_Type;       ///< Type of shock fix: "NONE", "HAENEL", "ALPHA_LR"  (only for AUSMDV)
  string Shock_Fix_Flag;       ///< Type of flagging for shock fix: "CELL_FLAG_HEX", "CELL_FLAG_TET" (only for AUSMDV)
  double Shock_Fix_Flag_Sum;   ///< Sum of flagged cells: > 1.0 for CELL_FLAG_HEX, >= 1.0 for CELL_FLAG_TET
  string Entropy_Fix;          ///< Type of entropy fix: "NONE", "EXPANSION_FAN"  (only for AUSMDV)
  double *ShockFixFlagCell;    ///< Cells flagged for special Euler flux treatment (shock fix) at compressible sonic points

  string SymmPressBC;          ///< Type of treatment for symmetry boundary condition: "ANALYTICAL", "ORIGINAL_NUMERIC", or "NEW_NUMERIC"


	// (begin) JST variables
	double  *PressSensor;
	double  (*Conserv_Und_Lapl)[5];
	double  *Lambda;
	double *p1_Und_Lapl; double *p2_Und_Lapl; 
	double *t1_Und_Lapl; double *t2_Und_Lapl; 
	double *mut1_Und_Lapl; double *mut2_Und_Lapl;
	bool *BoundaryCV; int *NeighborCV;
	double *Lambda_mgLevel1; int *NeighborCV_mgLevel1;
	double *Lambda_mgLevel2; int *NeighborCV_mgLevel2;
	double *Lambda_mgLevel3; int *NeighborCV_mgLevel3;
	
	// (end) JST variables

	// (begin) Multigrid variables
  double  *MGSensor;
	
	int *cvora_mgLevel1;
  double  *rho_mgLevel1;
  double (*rhou_mgLevel1)[3]; 
  double  *rhoE_mgLevel1;
  double (*vel_mgLevel1)[3];
  double *press_mgLevel1; 
  double *temp_mgLevel1;
  double *enthalpy_mgLevel1;
  double *RoM_mgLevel1;
	double *muT_mgLevel1;
  double *gamma_mgLevel1;
  double *sos_mgLevel1;
  double *kine_mgLevel1;
	double *mul_fa_mgLevel1;
	double *mut_fa_mgLevel1;
	double *lamOcp_fa_mgLevel1;
	double *local_dt_mgLevel1;
	double (*grad_enthalpy_mgLevel1)[3];
  double (*grad_u_mgLevel1)[3][3];
	double *muLam_mgLevel1;
	double *LambdaOverCp_mgLevel1;			
	double  *rho_old_mgLevel1;
  double (*rhou_old_mgLevel1)[3]; 
  double  *rhoE_old_mgLevel1;
  double  *rho_TruncError_mgLevel1;
  double (*rhou_TruncError_mgLevel1)[3]; 
  double  *rhoE_TruncError_mgLevel1;
  double  *rho_TruncError;
  double (*rhou_TruncError)[3]; 
  double  *rhoE_TruncError;
  double  *rho_res_old;
  double (*rhou_res_old)[3]; 
  double  *rhoE_res_old;
  double  *rho_res_sum;
  double (*rhou_res_sum)[3]; 
  double  *rhoE_res_sum;
  double *rho_bfa_mgLevel1; 
	double *T_bfa_mgLevel1; 
	double (*vel_bfa_mgLevel1)[3]; 
	double *p_bfa_mgLevel1; 
	double *h_bfa_mgLevel1; 
	double *gam_bfa_mgLevel1; 
	double *RoM_bfa_mgLevel1;
	double  *MGSensor_mgLevel1;
	
	int *cvora_mgLevel2;
	double  *rho_mgLevel2;
  double (*rhou_mgLevel2)[3]; 
  double  *rhoE_mgLevel2;
  double (*vel_mgLevel2)[3];
  double *press_mgLevel2; 
  double *temp_mgLevel2;
  double *enthalpy_mgLevel2;
  double *RoM_mgLevel2;
	double *muT_mgLevel2;
  double *gamma_mgLevel2;
  double *sos_mgLevel2;
  double *kine_mgLevel2;
	double *mul_fa_mgLevel2;
	double *mut_fa_mgLevel2;
	double *lamOcp_fa_mgLevel2;
	double *local_dt_mgLevel2;
	double (*grad_enthalpy_mgLevel2)[3];
  double (*grad_u_mgLevel2)[3][3];
	double *muLam_mgLevel2;
	double *LambdaOverCp_mgLevel2;			
	double  *rho_old_mgLevel2;
  double (*rhou_old_mgLevel2)[3]; 
  double  *rhoE_old_mgLevel2;
  double  *rho_TruncError_mgLevel2;
  double (*rhou_TruncError_mgLevel2)[3]; 
  double  *rhoE_TruncError_mgLevel2;
  double *rho_bfa_mgLevel2; 
	double *T_bfa_mgLevel2; 
	double (*vel_bfa_mgLevel2)[3]; 
	double *p_bfa_mgLevel2; 
	double *h_bfa_mgLevel2; 
	double *gam_bfa_mgLevel2; 
	double *RoM_bfa_mgLevel2;
	double  *MGSensor_mgLevel2;
	
	int *cvora_mgLevel3;
	double  *rho_mgLevel3;
  double (*rhou_mgLevel3)[3]; 
  double  *rhoE_mgLevel3;
  double (*vel_mgLevel3)[3];
  double *press_mgLevel3; 
  double *temp_mgLevel3;
  double *enthalpy_mgLevel3;
  double *RoM_mgLevel3;
	double *muT_mgLevel3;
  double *gamma_mgLevel3;
  double *sos_mgLevel3;
  double *kine_mgLevel3;
	double *mul_fa_mgLevel3;
	double *mut_fa_mgLevel3;
	double *lamOcp_fa_mgLevel3;
	double *local_dt_mgLevel3;
	double (*grad_enthalpy_mgLevel3)[3];
  double (*grad_u_mgLevel3)[3][3];
	double *muLam_mgLevel3;
	double *LambdaOverCp_mgLevel3;			
	double  *rho_old_mgLevel3;
  double (*rhou_old_mgLevel3)[3]; 
  double  *rhoE_old_mgLevel3;
  double  *rho_TruncError_mgLevel3;
  double (*rhou_TruncError_mgLevel3)[3]; 
  double  *rhoE_TruncError_mgLevel3;
  double *rho_bfa_mgLevel3; 
	double *T_bfa_mgLevel3; 
	double (*vel_bfa_mgLevel3)[3]; 
	double *p_bfa_mgLevel3; 
	double *h_bfa_mgLevel3; 
	double *gam_bfa_mgLevel3; 
	double *RoM_bfa_mgLevel3;
	double  *MGSensor_mgLevel3;
	// (end) Multigrid variables
  
protected: 

  int linearSolverNS;         ///< linear solver for Navier-Stokes equations
  int linearSolverScalars;    ///< linear solver for scalar equations

  PetscSolver *petscSolver, *petscSolverScalars;   ///< petsc's linear solvers
  int *cvora5;
  int *nbocv_v_global;
  
	
	// (begin) Multigrid variables
	PetscSolver *petscSolver_mgLevel1, *petscSolverScalars_mgLevel1;   ///< petsc's linear solvers
	int *cvora5_mgLevel1;
	PetscSolver *petscSolver_mgLevel2, *petscSolverScalars_mgLevel2;   ///< petsc's linear solvers
	int *cvora5_mgLevel2;
	PetscSolver *petscSolver_mgLevel3, *petscSolverScalars_mgLevel3;   ///< petsc's linear solvers
	int *cvora5_mgLevel3;
	int *nbocv_v_global_mgLevel1;
	int *nbocv_v_global_mgLevel2;
	int *nbocv_v_global_mgLevel3;
	// (end) Multigrid variables
	
  string initial_flowfield_output;       ///< tag (YES or NO) for output of initial field


public:   // member functions

  /**
   *  solve coupled linear system for Navier-Stokes equations
   */
  void solveCoupledLinSysNS(double (*phi)[5],
                            double (*Ap)[5][5],
                            double (*rhs)[5],
                            const double zeroAbs,
                            const double zeroRel,
                            const int maxIter)
  {
    switch (linearSolverNS)
    {
      case PETSC_GMRES:

        // on the first time, instantiate the petsc solver...
        if (petscSolver == NULL)
        {
          if (nbocv_v_global == NULL)
          {
            nbocv_v_global = new int[ncv_g];
            for (int icv = 0; icv < ncv; icv++)
              nbocv_v_global[icv] = cvora[mpi_rank] + icv;
            updateCvData(nbocv_v_global, REPLACE_DATA);
          }

          petscSolver = new PetscSolver(cvora, nbocv_i, nbocv_v, 5);
          petscSolver->setTresholds(zeroAbs, zeroRel, maxIter);
        }

        petscSolver->solveGMRES(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global, 5);

        break;

      case BCGSTAB:

        solveCvVectorR5Bcgstab(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq");   // solve the system
        break;
				
			case LUSGS:
				
				solveCvVectorR5Lusgs(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq");
				break;
				
			case BCGSTAB_LINELET:
				
        solveCvVectorR5BcgstabLine(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq");   // solve the system
				break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }

  void solveLinSysScalar(double *phi, double *Ap, double *rhs,
      const double zeroAbs, const double zeroRel, const int maxIter, char *scalarName)
  {
    switch (linearSolverScalars)
    {
      case PETSC_GMRES:

        // on the first time, instantiate the petsc solver...
        if (petscSolverScalars == NULL)
        {
          if (nbocv_v_global == NULL)
          {
            nbocv_v_global = new int[ncv_g];
            for (int icv = 0; icv < ncv; icv++)
              nbocv_v_global[icv] = cvora[mpi_rank] + icv;
            updateCvData(nbocv_v_global, REPLACE_DATA);
          }
          petscSolverScalars = new PetscSolver(cvora, nbocv_i, nbocv_v, 1);
        }

        petscSolverScalars->setTresholds(zeroAbs, zeroRel, maxIter);
        petscSolverScalars->solveGMRES(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global);

        break;

      case BCGSTAB:
        solveCvScalarBcgstab(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, scalarName);
        break;
				
			case LUSGS:
        solveCvScalarLusgs(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, scalarName);
        break;
				
			case BCGSTAB_LINELET:
        solveCvScalarBcgstabLine(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, scalarName);
        break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver for scalars: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }

  /**
   *  solve coupled linear system for Navier-Stokes equations and scalars together
   */
  void solveCoupledLinSysNSCoupled(double **phi,
                                   double ***Ap,
                                   double **rhs,
                                   const double zeroAbs,
                                   const double zeroRel,
                                   const int maxIter,
                                   int nScal)
  {
    switch (linearSolverNS)
    {
      case PETSC_GMRES:

        // on the first time, instantiate the petsc solver...
        if (petscSolver == NULL)
        {
          if (nbocv_v_global == NULL)
          {
            nbocv_v_global = new int[ncv_g];
            for (int icv = 0; icv < ncv; icv++)
              nbocv_v_global[icv] = cvora[mpi_rank] + icv;
            updateCvData(nbocv_v_global, REPLACE_DATA);
          }

          petscSolver = new PetscSolver(cvora, nbocv_i, nbocv_v, 5+nScal);
          petscSolver->setTresholds(zeroAbs, zeroRel, maxIter);
        }

        petscSolver->solveGMRESCoupled(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global, 5+nScal);

        break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }

  /**
   *    calc time step
   */
  virtual double calcDt(double cfl = 0.0);


  /**
   *
   *  explicit euler flux HLLC
   *
   */
  virtual int calcEulerFlux_HLLC(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc);


  /**
   *
   *  implicit euler flux HLLC supersonic
   *
   */
  virtual int calcEulerFluxMatrices_HLLC(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc);
  

  /**
   *
   *  implicit euler flux for HLLC subsonic
   *
   */
  void calcSubSonicJacobeanHLLC(double (*AL)[5], double (*AR)[5],
      double rhoL, const double *uL, double pL, double eL, double qL, double psiL, double SL,
      double rhoSL, double *rhouSL, double eSL,
      double *dSMdUL, double *dSMdUR, double *dpsdUL, double *dpsdUR, double SM, double pS,
      double gamma, const double *nV);

	/**
   *
   *  explicit euler flux Roe
   *
   */
  virtual int calcEulerFlux_Roe(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
																const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
																const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
																const double area, const double *nVec, const int nScal, const double surfVeloc);
	
	/**
	 *
	 *  implicit euler flux for Roe
	 *
	 */
	virtual int calcEulerFluxMatrices_Roe(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
										  const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
										  const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
										  const double area, const double *nVec, const int nScal, const double surfVeloc);
	
	/**
	 *
	 *  explicit euler flux JST
	 *
	 */
	virtual int calcEulerFlux_JST(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
								  const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
								  const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
								  const double area, const double *nVec, const int nScal, const double surfVeloc,
								  const double SensorL, const double LambdaL, const double *Und_LaplL, const int NeighborL, const double SensorR, 
								  const double LambdaR, const double *Und_LaplR, const int NeighborR);
	
	virtual int calcEulerFlux_Lax(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
								  const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
								  const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
								  const double area, const double *nVec, const int nScal, const double surfVeloc,
								  const double LambdaL, const int NeighborL, const double LambdaR, const int NeighborR);
	
	virtual void calcPMatrix(double P[5][5], const double *vel, double c, double rrho, 
													 const double *nV, double gamma, double surfVeloc);
	
	virtual void calcInvPMatrix(double invP[5][5], const double *vel, double c, double rrho, 
															const double *nV, double gamma, double surfVeloc);

	
  // nV is not normalized
  void calcJacobianA(double (*A)[5], const double *vel, double pp, double rrho,
      const double *nV, double gamma, double surfVeloc);
  
	/**
   *
   *  explicit euler flux AUSM
   *
   */
  virtual int calcEulerFlux_AUSM(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
																 const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
																 const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
																 const double area, const double *nVec, const int nScal, const double surfVeloc);

  /**
   *
   *  explicit euler flux AUSMDV
   *
   */
  virtual int calcEulerFlux_AUSMDV(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
                                   const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
                                   const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
                                   const double area, const double *nVec, const int nScal, const double surfVeloc,
                                   const string AUSMType, const string ShockFix, const string EntropyFix);

  /**
   *
   *  implicit euler flux AUSMDV
   *
   */
  virtual int calcEulerFluxMatrices_AUSMDV(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc,
      const string AUSMType, const string ShockFix, const string EntropyFix);

  /**
   *
   *  explicit euler flux based on Haenel & Schwane
   *
   */
  virtual int calcEulerFlux_HAENEL(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
                                   const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
                                   const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
                                   const double area, const double *nVec, const int nScal, const double surfVeloc);

  /**
   *
   *  implicit euler flux based on Haenel & Schwane
   *
   */
  virtual int calcEulerFluxMatrices_HAENEL(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc);
	
  /**
   *
   * viscous flux routine for implicit
   * needs lsg_coeff0 and lsg_coeff1 ---> LS gradient coefficients
   *
   */
  virtual void addViscFlux(double *Frhou, double &FrhoE, double (*A0)[5], double (*A1)[5],
      const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double kine0,
      const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double kine1,
      const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, 
      const double area, const double *nVec, const double smag, const double *sVec);
  
  
 
  
  

  /**
   *  explicit euler flux HLLC, coupled NSE and scalars
   */
  virtual void calcEulerFluxCoupled_HLLC(double *EulerFlux,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);

  /**
   *  matrices for implicit euler flux HLLC, coupled NSE and scalars
   */
  virtual void calcEulerFluxMatricesCoupled_HLLC(double **A_L, double **A_R,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double *dpress_dscalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double *dpress_dscalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);
  
  /**
   *  supersonic Jacobian matrix for Euler flux, coupled NSE and scalars
   */  
  void calcJacobianCoupled_HLLC(double **A, const double *vel, const double pp, const double rrho, const double hh, const double *scal, const double *dpress_dscal, 
      const double *nV, const double gamma, const double surfVeloc, const int nScal);

  /**
   *  derivative of contact surface speed
   */  
  void calcdSMdU_HLLC(double *dSMdU, const double *vel, const double pp, const double rrho, const double hh, const double *dpress_dscal,
                 const double *nV, const double gamma, const double SM, const double S, const double invertrho, const double Factor, const int nScal);

  /**
   *  subsonic Jacobian for Euler flux, coupled NSE and scalars
   */  
  void calcSubSonicJacobianCoupled_HLLC(double **A, const double *dSMdU, const double rhoS, const double *rhouS, const double rhoeS, const double *scal,
                                        const double pS, const double Omega, const double sM, const double rhoSmq, const double *nV, const double area, const int nScal);
  

  /**
   *  explicit euler flux AUSMDV, coupled NSE and scalars
   */
  virtual void calcEulerFluxCoupled_AUSMDV(double *EulerFlux,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment,
      const string AUSMType, const string ShockFix, const string EntropyFix);

  /**
   *  matrices for implicit euler flux AUSMDV, coupled NSE and scalars
   */
  virtual void calcEulerFluxMatricesCoupled_AUSMDV(double **A_L, double **A_R,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double *dpress_dscalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double *dpress_dscalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment,
      const string AUSMType, const string ShockFix, const string EntropyFix);

  /**
   *  explicit euler flux based on Haenel & Schwane, coupled NSE and scalars
   */
  virtual void calcEulerFluxCoupled_HAENEL(double *EulerFlux,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);

  /**
   *  matrices for implicit euler flux based on Haenel & Schwane, coupled NSE and scalars
   */
  virtual void calcEulerFluxMatricesCoupled_HAENEL(double **A_L, double **A_R,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double *dpress_dscalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double *dpress_dscalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);

  /**
   *  explicit and implicit viscous flux, coupled NSE and scalars
   */
  virtual void calcViscousFluxCoupled(double *ViscousFlux, double **A0, double **A1,
      const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double *Scal0, const double (*gradScal0)[3], const double *dpress_dscal0, const double kine0,
      const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double *Scal1, const double (*gradScal1)[3], const double *dpress_dscal1, const double kine1,
      const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, const double *diff, const double *DiffTerm,
      const double area, const double *nVec, const double smag, const double *sVec, const double alpha, const int nScal);
  
  
  

  
  /**
   * loop over scalars and solve them
   */
  void setScalarBC(FaZone *zone);
  
  void solveScalars(double *massFlux);
  
  void solveScalarsImplUnsteady(double *massFlux);
  
  void solveScalarsExplicit(double *massFlux);

  void solveScalarTransport(ScalarTranspEq *scal, double *rhs, double *A, double *massFlux);
  
  void solveScalarTransportImplUnsteady(ScalarTranspEq *scal, double *rhs, double *A, double *massFlux);
  
  void solveScalarTransportExplicit(ScalarTranspEq *scal, double *rhs, double *massFlux_fa);
  
  virtual void UserDefinedScalarClipping(const string &name)  {/*empty*/}
  
  void calcViscousFluxScalar_new(double *rhs_rhoScal, double *Ascal, ScalarTranspEq &transpScal, int flagImplicit);

  
  

  
  virtual void calcGradVel()
  {
    calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, limiterNavierS, sos, epsilonSDWLS);
  }

  virtual void calcStrainRateAndDivergence()
  {
    if ((strMag == NULL) || (diverg == NULL))
    {
      cout << "you have tried to calculate the magnitude of the strain rate and the divergence of U, but one of them has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv; icv++)
    {
      diverg[icv] = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      strMag[icv] = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i])-1./3.*diverg[icv], 2.0);
        else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);

      strMag[icv] = sqrt(2.0*strMag[icv]);
    }
    updateCvData(strMag, REPLACE_DATA);
    updateCvData(diverg, REPLACE_DATA);
  }

  virtual void calcDivergence()
  {
    if (diverg == NULL)
    {
      cout << "you have tried to calculate the divergence of U, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv; icv++)
      diverg[icv] = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

    updateCvData(diverg, REPLACE_DATA);
  }
  
  virtual void calcStrainRate()
  {
    if (strMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the strain rate, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv; icv++)
    {
      strMag[icv] = 0.0;
      double diverg = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]) - 1./3.*diverg, 2.0);
        else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);

      strMag[icv] = sqrt(2.0*strMag[icv]);
    }
    updateCvData(strMag, REPLACE_DATA);
  }

  virtual void calcVorticity()
  {
    if (vortMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the vorticity, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv; icv++)
    {
      vortMag[icv] = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        vortMag[icv] += pow(0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]), 2.0);

      vortMag[icv] = sqrt(2.0*vortMag[icv]);
    }
    updateCvData(vortMag, REPLACE_DATA);
  }


  /** 
   * calculate wall distance: minimum distance of cv to wall face by looping over all cv's and faces 
   */
  void calcWallDistance(double *phi, double *wd)
  {
    if (mpi_rank == 0)
      cout << "calc wall distance" << endl;

    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    for (int icv=0; icv<ncv; icv++)
      wd[icv] = 1.0e20;

    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
          {
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
              my_wall_faces++;
          }
      }

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    if (tot_wall_faces == 0)
    {
      if (mpi_rank == 0)
        cerr << "ERROR: you have been trying to calculate the wall distance, but there are no walls in you domain ?!?" << endl;
      throw(-1);
    }

    // set up send side
    int send_count[mpi_size], send_displ[mpi_size];
    for (int r=0; r<mpi_size; r++)
    {
      send_count[r] = my_wall_faces*3;
      send_displ[r] = 0;
    }

    // set up receive side
    int recv_count[mpi_size], recv_displ[mpi_size];
    int wallcoord = my_wall_faces*3;
    MPI_Allgather(&wallcoord, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    
    recv_displ[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_displ[i] = recv_count[i-1] + recv_displ[i-1];

    // fill up send buffer
    int count = 0;
    double *my_wall_fa        = new double[my_wall_faces*3];
    double *my_wall_fa_normal = new double[my_wall_faces*3];
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            {
              for (int i=0; i<3; i++)
              {
                my_wall_fa[3*count + i] = x_fa[ifa][i];
                my_wall_fa_normal[3*count + i] = fa_normal[ifa][i];
              }
  
              count++;
            }
      }
    
    // receive all wall face coordinates
    double *wall_fa   = new double[tot_wall_faces*3];
    double *wall_fa_n = new double[tot_wall_faces*3];
    
    MPI_Alltoallv(my_wall_fa, send_count, send_displ, MPI_DOUBLE, 
                  wall_fa,    recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

    MPI_Alltoallv(my_wall_fa_normal, send_count, send_displ, MPI_DOUBLE, 
                  wall_fa_n,         recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

    delete [] my_wall_fa;
    delete [] my_wall_fa_normal;


    // find minimum wall distance 
    for (int icv=0; icv<ncv; icv++)
    {
      double minDist = 1.0e20;
      int ifaSave;

      for (int ifa=0; ifa<tot_wall_faces; ifa++)
      {
        double dist = sqrt(   pow(x_cv[icv][0]-wall_fa[3*ifa  ], 2.0)
                            + pow(x_cv[icv][1]-wall_fa[3*ifa+1], 2.0)
                            + pow(x_cv[icv][2]-wall_fa[3*ifa+2], 2.0));
        
        if (dist < minDist)
        {
          ifaSave = ifa;
          minDist = dist;
        }
      }
      
/*      double fa_x[3] = {wall_fa[3*ifaSave], wall_fa[3*ifaSave+1], wall_fa[3*ifaSave+2]};
      double fa_n[3] = {wall_fa_n[3*ifaSave], wall_fa_n[3*ifaSave+1], wall_fa_n[3*ifaSave+2]};
      double n[3], s_half[3];
      
      normVec3d(n, fa_n);
      vecMinVec3d(s_half, fa_x, x_cv[icv]);
      wd[icv] = fabs(vecDotVec3d(s_half, n));*/
      
      double fa_x[3] = {wall_fa[3*ifaSave], wall_fa[3*ifaSave+1], wall_fa[3*ifaSave+2]};
      double fa_n[3] = {wall_fa_n[3*ifaSave], wall_fa_n[3*ifaSave+1], wall_fa_n[3*ifaSave+2]};
      double n[3], s_half[3];
      
      double nmag = normVec3d(n, fa_n);
      vecMinVec3d(s_half, fa_x, x_cv[icv]);

      double sMag2 = vecDotVec3d(s_half, s_half);
      double dNorm2 = vecDotVec3d(s_half, n)*vecDotVec3d(s_half, n);
      double dTang2 = sMag2-dNorm2;

      double specRadius2 = 0.5*nmag;
      

      if (dTang2 > specRadius2)   wd[icv] = minDist;
      else                        wd[icv] = fabs(vecDotVec3d(s_half, n));


//      wd[icv] = minDist;
    }
    /*
    // correct closest wall
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if (param->getString() == "WALL")
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          
          double n[3], s_half[3];
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          wd[icv0] = fabs(vecDotVec3d(s_half, n));  //normVec3d(s_half);          
        }
    }*/
    
    MPI_Barrier(mpi_comm);
    
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > time to compute wall distance [s]: " << wtime - wtime0 << endl;
    }
    
    delete [] wall_fa;
    delete [] wall_fa_n;    
  }
  
  // Interpolate variable at the cell center from values at the faces (e.g., for viscosity)
  double InterpolateAtCellCenterFromFaceValues(double *phi, int icv)
  {
    double phiC = 0.0;

    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv + 1] - 1;
    for (int foc=foc_f; foc<=foc_l; foc++)
    {
      int ifa = faocv_v[foc];
      phiC += phi[ifa];
    }
    phiC /= (double)(foc_l-foc_f+1);    
    
    return phiC;
  }
  
  
  
  
  
  
  
  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // turbmodels turbmodels turbmodels turbmodels turbmodels turbmodels 
  //
  // ---------------------------------------------------------------------------
  // ###########################################################################
  
  virtual void calcRansTurbViscMuet()
  {
    static int flag = 1;
    
    // provide zero mut for laminar calculations
    if (flag == 1)
    {
      flag = 0;
      for (int ifa=0; ifa<nfa; ifa++)
        mut_fa[ifa] = 0.0;
    }
  }

  virtual void initialHookScalarRansTurbModel() {/*empty*/}
  virtual void diffusivityHookScalarRansTurb(const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansTurbExplicit(double *rhs, const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit) {/*empty*/}
  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}  
  virtual void sourceHookRansTurb(double *rhsRho, double (*rhsRhou)[3], double *rhsRhoE, double (*A)[5][5]) {/*empty*/}
  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit) {/*empty*/}
  virtual void pressureDerivativeHookScalarRansTurb() {/*empty*/}
  
  

  
  
  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // thermodynamics thermodynamics thermodynamics thermodynamics 
  //
  // ---------------------------------------------------------------------------  
  // ###########################################################################
  
  virtual double calcMuLam(double temp)
  {
    double muLam = 0.0;
    
    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp + SL_Sref); 
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }
  
  virtual double calcMuLam(int icv)
  {
    double muLam = 0.0;
    
    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp[icv]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv] + SL_Sref); 
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp[icv]/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }
  
  
  /*! \brief Compute constants gas constant R_gas and heat capacity ratio gamma.
   *
   *  Compute R_gas and gamma and initialize the arrays RoM and gam.
   */
  virtual void initialHookScalarRansCombModel()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      RoM[icv] = R_gas;
      gamma[icv] = GAMMA;
    }

    updateCvData(RoM, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
  }



  /*! \brief Compute for each cell the states of variables and material properties based on constant gamma and R/M.
   *
   *  Compute temperature temp, viscosity mul, heat coefficient lambda,
   *  enthalpy enthalp and pressure press.
   */
  virtual void calcRansStateVarAndMaterialProperties(double *rho, double(*rhou)[3], double *rhoE)
  {
    // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
    calcStateVariables();
    
    // Compute viscosity and thermal diffusivity at cell faces
    calcMaterialProperties();
  }

  /*! \brief Compute for each cell the states of variables based on constant gamma and R/M.
   *
   *  Compute velocity vel, temperature temp, enthalpy enthalpy, pressure press, specific heat ratio gamma, speed of sounds sos and gas constant R.
   *  R and gamma not actually computed since constant
   *  
   */
  virtual void calcStateVariables()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      if (rho[icv] <= 0.0)
        cout << "negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      double kinecv = 0.0;
      if (kine != NULL)
        kinecv = kine[icv];

      double pr = (gamma[icv] - 1.0) * (rhoE[icv] - 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / rho[icv] - rho[icv] * kinecv);
      if (pr <= 0.0)
      {
        cout << "negative pressure at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
      }
      else
        press[icv] = pr;

      temp[icv] = press[icv] / (rho[icv] * RoM[icv]);
      enthalpy[icv] = gamma[icv] * RoM[icv] / (gamma[icv] - 1.0) * temp[icv];
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);
    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
  }

  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *  
   */
  virtual void calcMaterialProperties()
  {
    if (mu_ref > 0.0)
    {
      if (viscMode == "SUTHERLAND")
      {
        // internal faces
        for (int ifa = nfa_b; ifa < nfa; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
          vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
          vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
          double w0 = sqrt(vecDotVec3d(dx0, dx0));
          double w1 = sqrt(vecDotVec3d(dx1, dx1));

          double temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
          mul_fa[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
          lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
        }

          // boundary faces computed next in setBC
      } 
      else if (viscMode == "POWERLAW")
      {
        // internal faces
        for (int ifa = nfa_b; ifa < nfa; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
          vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
          vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
          double w0 = sqrt(vecDotVec3d(dx0, dx0));
          double w1 = sqrt(vecDotVec3d(dx1, dx1));
   
          double temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
          mul_fa[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
          lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
        }
  
        // boundary faces computed next in setBC
      }
      else
      {
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1);
      }
    }
  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    p = rho * R * T;
    h = gam * R / (gam - 1.0) * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    T = p / (rho * R);
    h = gam * R / (gam - 1.0) * T;
  }
  
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_e(double &T, double &p, double &R, double &gam, double &rho, double &e, double *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    T = e * (gam - 1.0) / R;
    p = rho * R * T;
  }

  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    {
      gam_bfa[ifa] = GAMMA;
      RoM_bfa[ifa] = R_gas;
      h_bfa[ifa] = GAMMA * R_gas / (GAMMA - 1.0) * T_bfa[ifa];
    }

    if (mu_ref > 0.0)
    {
      if (viscMode == "SUTHERLAND")        
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          mul_fa[ifa] = mu_ref*pow(T_bfa[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa[ifa] + SL_Sref);
          //mul_fa[ifa] = 1.458e-6 * pow(T_bfa[ifa], 1.5) / (T_bfa[ifa] + 110.4);
      else if (viscMode == "POWERLAW")     
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          mul_fa[ifa] = mu_ref * pow(T_bfa[ifa]/T_ref, mu_power_law);
      else
      {
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1);
      }
      
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
    }
  }  
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    {
      gam_bfa[ifa] = GAMMA;
      RoM_bfa[ifa] = R_gas;
      T_bfa[ifa] = h_bfa[ifa] * (GAMMA - 1.0) / (GAMMA * R_gas);
    }

    if (mu_ref > 0.0)
    {
      if (viscMode == "SUTHERLAND")        
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          mul_fa[ifa] = mu_ref*pow(T_bfa[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa[ifa] + SL_Sref);
          //mul_fa[ifa] = 1.458e-6 * pow(T_bfa[ifa], 1.5) / (T_bfa[ifa] + 110.4);
      else if (viscMode == "POWERLAW")     
        for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          mul_fa[ifa] = mu_ref * pow(T_bfa[ifa]/T_ref, mu_power_law);
      else
      {
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
    }
  }



  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // combustion combustion combustion combustion combustion combustion
  //
  // ---------------------------------------------------------------------------  
  // ###########################################################################

  virtual void diffusivityHookScalarRansComb(const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansComb(double *rhs, double *A, const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansCombExplicit(double *rhs, const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansComb_new(double *rhs, double *A, const string &name, int flagImplicit)  {/*empty*/}
  virtual void boundaryHookScalarRansComb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}  
  virtual void sourceHookRansComb(double *rhsRho, double(*rhsRhou)[3], double *rhsRhoE, double(*A)[5][5])  {/*empty*/}  
  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit)  {/*empty*/}  
  virtual void pressureDerivativeHookScalarRansComb() {/*empty*/}
  
  

  int zoneIsWall(const string & faName) 
  {
    Param *p;
    string paraName;
    int retVal = 0;

    if (getParam(p, faName)) {
      paraName = p->getString();
      retVal = (paraName == "WALL") || (paraName == "Wall") || (paraName == "wall") ||
               (paraName == "WALL.T") || (paraName == "Wall.T") || (paraName == "wall.T");
    }
    return retVal;
  }

  // #########################################################################################
  //
  //    scalar boundaries / scalar boundaries / scalar boundaries
  //
  // #########################################################################################
  int scalarZoneIsHook(const string &faName, const string &scalName) {

    Param *p;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName)) {
      string name = p->getString(1);
      retVal = (name == "HOOK") || (name == "Hook")|| (name == "hook");
    } else
      retVal = 0;

    return retVal;
  }

  int scalarZoneIsDirichlet(double &phiBC, const string &faName, const string &scalName) {

    Param *p;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName)) {
      phiBC = p->getDouble(1);
      retVal = 1;
    } else
      retVal = 0;

    return retVal;
  }

  int scalarZoneIsFlux(double phiBCflux, const string &faName, const string &scalName) {

    Param *p;
    string paraName;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName + ".flux")) {
      phiBCflux = p->getDouble(1);
      retVal = 1;
    } else
      retVal = 0;

    return retVal;
  }

// Begin Multigrid
	
	void solveCoupledLinSysNS_mg(double (*phi)[5], double (*Ap)[5][5],
															 double (*rhs)[5], const double zeroAbs,
															 const double zeroRel, const int maxIter, int iMesh) {
		switch (linearSolverNS)
		{
			case PETSC_GMRES:
				
				if (iMesh == 1) {
					
					if (petscSolver_mgLevel1 == NULL) {
						if (nbocv_v_global_mgLevel1 == NULL) {
							nbocv_v_global_mgLevel1 = new int[ncv_g_mgLevel1];
							for (int icv = 0; icv < ncv_mgLevel1; icv++) {
								nbocv_v_global_mgLevel1[icv] = cvora_mgLevel1[mpi_rank] + icv;
							}
							
							updateCvData_mg(nbocv_v_global_mgLevel1, REPLACE_DATA, 1);
						}
						
						
						petscSolver_mgLevel1 = new PetscSolver(cvora_mgLevel1, nbocv_i_mgLevel1, nbocv_v_mgLevel1, 5);
						petscSolver_mgLevel1->setTresholds(zeroAbs, zeroRel, maxIter);
					}
					
					petscSolver_mgLevel1->solveGMRES(Ap, phi, rhs, cvora_mgLevel1, nbocv_i_mgLevel1, nbocv_v_mgLevel1, nbocv_v_global_mgLevel1, 5);
				}
				
				if (iMesh == 2) {
					
					if (petscSolver_mgLevel2 == NULL) {
						if (nbocv_v_global_mgLevel2 == NULL) {
							nbocv_v_global_mgLevel2 = new int[ncv_g_mgLevel2];
							for (int icv = 0; icv < ncv_mgLevel2; icv++)
								nbocv_v_global_mgLevel2[icv] = cvora_mgLevel2[mpi_rank] + icv;
							updateCvData_mg(nbocv_v_global_mgLevel2, REPLACE_DATA, 2);
						}
						
						petscSolver_mgLevel2 = new PetscSolver(cvora_mgLevel2, nbocv_i_mgLevel2, nbocv_v_mgLevel2, 5);
						petscSolver_mgLevel2->setTresholds(zeroAbs, zeroRel, maxIter);
					}
					
					petscSolver_mgLevel2->solveGMRES(Ap, phi, rhs, cvora_mgLevel2, nbocv_i_mgLevel2, nbocv_v_mgLevel2, nbocv_v_global_mgLevel2, 5);
				}
				
				if (iMesh == 3) {
					if (petscSolver_mgLevel3 == NULL) {
						if (nbocv_v_global_mgLevel3 == NULL) {
							nbocv_v_global_mgLevel3 = new int[ncv_g_mgLevel3];
							for (int icv = 0; icv < ncv_mgLevel3; icv++)
								nbocv_v_global_mgLevel3[icv] = cvora_mgLevel3[mpi_rank] + icv;
							updateCvData_mg(nbocv_v_global_mgLevel3, REPLACE_DATA, 3);
						}
						
						petscSolver_mgLevel3 = new PetscSolver(cvora_mgLevel3, nbocv_i_mgLevel3, nbocv_v_mgLevel3, 5);
						petscSolver_mgLevel3->setTresholds(zeroAbs, zeroRel, maxIter);
					}
					
					petscSolver_mgLevel3->solveGMRES(Ap, phi, rhs, cvora_mgLevel3, nbocv_i_mgLevel3, nbocv_v_mgLevel3, nbocv_v_global_mgLevel3, 5);
				}
				
				break;
				
			case BCGSTAB:
				
				solveCvVectorR5Bcgstab_mg(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq", iMesh);   // solve the system
				break;
				
			case LUSGS:
				
				solveCvVectorR5Lusgs_mg(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq", iMesh);   // solve the system
				break;
				
			default:
				if (mpi_rank == 0)
					cerr << "Error: unrecognized solver: " << linearSolverNS << endl;
				throw(-1);
				break;
		}
	}
	
	virtual double calcDt_mg(double cfl = 0.0, int iMesh = 1);
	
	virtual void addViscFlux_mg(double *Frhou, double &FrhoE, double (*A0)[5], double (*A1)[5],
															const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double kine0,
															const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double kine1,
															const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, 
															const double area, const double *nVec);
	
	double InterpolateAtCellCenterFromFaceValues_mg(double *phi, int icv, int iMesh)
	{
		double phiC = 0.0;
		
		if (iMesh == 1) {
			int foc_f = faocv_i_mgLevel1[icv];
			int foc_l = faocv_i_mgLevel1[icv + 1] - 1;
			for (int foc=foc_f; foc<=foc_l; foc++)
			{
				int ifa = faocv_v_mgLevel1[foc];
				phiC += phi[ifa];
			}
			phiC /= (double)(foc_l-foc_f+1);  
		}
		
		if (iMesh == 2) {
			int foc_f = faocv_i_mgLevel2[icv];
			int foc_l = faocv_i_mgLevel2[icv + 1] - 1;
			for (int foc=foc_f; foc<=foc_l; foc++)
			{
				int ifa = faocv_v_mgLevel2[foc];
				phiC += phi[ifa];
			}
			phiC /= (double)(foc_l-foc_f+1);  
		}
		
		return phiC;
	}
	
	virtual void calcRansTurbViscMuet_mg(int iMesh) {
		
		if (iMesh == 1) {
			for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
				int icv0 = cvofa_mgLevel1[ifa][0];
				int icv1 = cvofa_mgLevel1[ifa][1];
				mut_fa_mgLevel1[ifa] = 0.5*(muT_mgLevel1[icv0]+muT_mgLevel1[icv1]);
			}
			
			for (int ifa = 0; ifa < nfa_b_mgLevel1; ifa++) {
				int icv0 = cvofa_mgLevel1[ifa][0];
				mut_fa_mgLevel1[ifa] = muT_mgLevel1[icv0];
			}
			
			for (list<FaZone>::iterator zone = faZoneList_mgLevel1.begin(); zone != faZoneList_mgLevel1.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY)  {
					Param *param;
					if (getParam(param, zone->getName()))
						if (param->getString() == "WALL")
							for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
								mut_fa_mgLevel1[ifa] = 0.0;
				}
		}
		if (iMesh == 2) {
			for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
				int icv0 = cvofa_mgLevel2[ifa][0];
				int icv1 = cvofa_mgLevel2[ifa][1];
				mut_fa_mgLevel2[ifa] = 0.5*(muT_mgLevel2[icv0]+muT_mgLevel2[icv1]);
			}
			
			for (int ifa = 0; ifa < nfa_b_mgLevel2; ifa++) {
				int icv0 = cvofa_mgLevel2[ifa][0];
				mut_fa_mgLevel2[ifa] = muT_mgLevel2[icv0];
			}
			for (list<FaZone>::iterator zone = faZoneList_mgLevel2.begin(); zone != faZoneList_mgLevel2.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY)  {
					Param *param;
					if (getParam(param, zone->getName()))
						if (param->getString() == "WALL")
							for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
								mut_fa_mgLevel2[ifa] = 0.0;
				}
		}
		if (iMesh == 3) {
			for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
				int icv0 = cvofa_mgLevel3[ifa][0];
				int icv1 = cvofa_mgLevel3[ifa][1];
				mut_fa_mgLevel3[ifa] = 0.5*(muT_mgLevel3[icv0]+muT_mgLevel3[icv1]);
			}
			
			for (int ifa = 0; ifa < nfa_b_mgLevel3; ifa++) {
				int icv0 = cvofa_mgLevel3[ifa][0];
				mut_fa_mgLevel3[ifa] = muT_mgLevel3[icv0];
			}
			for (list<FaZone>::iterator zone = faZoneList_mgLevel3.begin(); zone != faZoneList_mgLevel3.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY)  {
					Param *param;
					if (getParam(param, zone->getName()))
						if (param->getString() == "WALL")
							for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
								mut_fa_mgLevel3[ifa] = 0.0;
				}
		}
	}
	
	virtual void calcStateVariables_mg(int iMesh)  {
		if (iMesh == 1) {
			for (int icv = 0; icv < ncv_mgLevel1; icv++) {
				
				if (rho_mgLevel1[icv] <= 0.0) {
					cout << "negative density at xcv (mg): " << rho_mgLevel1[icv] <<"-"<< x_cv_mgLevel1[icv][0] << ", " << x_cv_mgLevel1[icv][1] << ", " << x_cv_mgLevel1[icv][2] << endl;
				}
				
				for (int i=0; i<3; i++)
					vel_mgLevel1[icv][i] = rhou_mgLevel1[icv][i]/rho_mgLevel1[icv];
				
				double kinecv = 0.0;
				
				double pr = (gamma_mgLevel1[icv] - 1.0) * (rhoE_mgLevel1[icv] - 0.5 * vecDotVec3d(rhou_mgLevel1[icv], rhou_mgLevel1[icv]) / rho_mgLevel1[icv] - rho_mgLevel1[icv] * kinecv);
				//				if (pr <= 0.0) {
				//					cout << "negative pressure at xcv (mg): " << gamma_mgLevel1[icv]<< ", " << rhoE_mgLevel1[icv] << ", " << rho_mgLevel1[icv] << ", " << rhou_mgLevel1[icv][0] << ", " << rhou_mgLevel1[icv][1] << ", " << rhou_mgLevel1[icv][2] << endl;
				//				}
				//				else 
				press_mgLevel1[icv] = pr;
				
				temp_mgLevel1[icv] = press_mgLevel1[icv] / (rho_mgLevel1[icv] * RoM_mgLevel1[icv]);
				enthalpy_mgLevel1[icv] = gamma_mgLevel1[icv] * RoM_mgLevel1[icv] / (gamma_mgLevel1[icv] - 1.0) * temp_mgLevel1[icv];
				sos_mgLevel1[icv] = sqrt(gamma_mgLevel1[icv] * press_mgLevel1[icv] / rho_mgLevel1[icv]);
				
			}
		}
		if (iMesh == 2) {
			for (int icv = 0; icv < ncv_mgLevel2; icv++) {
				
				if (rho_mgLevel2[icv] <= 0.0)
					cout << "negative density at xcv (mg): " << x_cv_mgLevel2[icv][0] << ", " << x_cv_mgLevel2[icv][1] << ", " << x_cv_mgLevel2[icv][2] << endl;
				
				for (int i=0; i<3; i++)
					vel_mgLevel2[icv][i] = rhou_mgLevel2[icv][i]/rho_mgLevel2[icv];
				
				double kinecv = 0.0;
				
				double pr = (gamma_mgLevel2[icv] - 1.0) * (rhoE_mgLevel2[icv] - 0.5 * vecDotVec3d(rhou_mgLevel2[icv], rhou_mgLevel2[icv]) / rho_mgLevel2[icv] - rho_mgLevel2[icv] * kinecv);
				//				if (pr <= 0.0) {
				//					cout << "negative pressure at xcv (mg): " << gamma_mgLevel2[icv]<< ", " << rhoE_mgLevel2[icv] << ", " << rho_mgLevel2[icv] << ", " << rhou_mgLevel2[icv][0] << ", " << rhou_mgLevel2[icv][1] << ", " << rhou_mgLevel2[icv][2] << endl;
				//				}
				//				else 
				press_mgLevel2[icv] = pr;
				
				temp_mgLevel2[icv] = press_mgLevel2[icv] / (rho_mgLevel2[icv] * RoM_mgLevel2[icv]);
				enthalpy_mgLevel2[icv] = gamma_mgLevel2[icv] * RoM_mgLevel2[icv] / (gamma_mgLevel2[icv] - 1.0) * temp_mgLevel2[icv];
				sos_mgLevel2[icv] = sqrt(gamma_mgLevel2[icv] * press_mgLevel2[icv] / rho_mgLevel2[icv]);
				
			}
		}
		if (iMesh == 3) {
			for (int icv = 0; icv < ncv_mgLevel3; icv++) {
				
				if (rho_mgLevel3[icv] <= 0.0)
					cout << "negative density at xcv (mg): " << x_cv_mgLevel3[icv][0] << ", " << x_cv_mgLevel3[icv][1] << ", " << x_cv_mgLevel3[icv][2] << endl;
				
				for (int i=0; i<3; i++)
					vel_mgLevel3[icv][i] = rhou_mgLevel3[icv][i]/rho_mgLevel3[icv];
				
				double kinecv = 0.0;
				
				double pr = (gamma_mgLevel3[icv] - 1.0) * (rhoE_mgLevel3[icv] - 0.5 * vecDotVec3d(rhou_mgLevel3[icv], rhou_mgLevel3[icv]) / rho_mgLevel3[icv] - rho_mgLevel3[icv] * kinecv);
				//				if (pr <= 0.0) {
				//					cout << "negative pressure at xcv (mg): " << gamma_mgLevel3[icv]<< ", " << rhoE_mgLevel3[icv] << ", " << rho_mgLevel3[icv] << ", " << rhou_mgLevel3[icv][0] << ", " << rhou_mgLevel3[icv][1] << ", " << rhou_mgLevel3[icv][2] << endl;
				//				}
				//				else 
				press_mgLevel3[icv] = pr;
				
				temp_mgLevel3[icv] = press_mgLevel3[icv] / (rho_mgLevel3[icv] * RoM_mgLevel3[icv]);
				enthalpy_mgLevel3[icv] = gamma_mgLevel3[icv] * RoM_mgLevel3[icv] / (gamma_mgLevel3[icv] - 1.0) * temp_mgLevel3[icv];
				sos_mgLevel3[icv] = sqrt(gamma_mgLevel3[icv] * press_mgLevel3[icv] / rho_mgLevel3[icv]);
				
			}
		}
	}
	
	virtual void calcMaterialProperties_mg(int iMesh) {
		if (iMesh == 1) {
			if (mu_ref > 0.0) {
				if (viscMode == "SUTHERLAND") {
					// internal faces
					for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
						
						int icv0 = cvofa_mgLevel1[ifa][0];
						int icv1 = cvofa_mgLevel1[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel1[icv0] + temp_mgLevel1[icv1]);
						mul_fa_mgLevel1[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
						lamOcp_fa_mgLevel1[ifa] = mul_fa_mgLevel1[ifa] / Pr;
					}
					// boundary faces computed next in setBC
				} 
				else if (viscMode == "POWERLAW") {
					// internal faces
					for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
						
						int icv0 = cvofa_mgLevel1[ifa][0];
						int icv1 = cvofa_mgLevel1[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel1[icv0] + temp_mgLevel1[icv1]);
						mul_fa_mgLevel1[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
						lamOcp_fa_mgLevel1[ifa] = mul_fa_mgLevel1[ifa] / Pr;
					}
					
					// boundary faces computed next in setBC
				}
				else {
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
			}
		}
		if (iMesh == 2) {
			if (mu_ref > 0.0) {
				if (viscMode == "SUTHERLAND") {
					// internal faces
					for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
						
						int icv0 = cvofa_mgLevel2[ifa][0];
						int icv1 = cvofa_mgLevel2[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel2[icv0] + temp_mgLevel2[icv1]);
						mul_fa_mgLevel2[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
						lamOcp_fa_mgLevel2[ifa] = mul_fa_mgLevel2[ifa] / Pr;
					}
					// boundary faces computed next in setBC
				} 
				else if (viscMode == "POWERLAW") {
					// internal faces
					for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
						
						int icv0 = cvofa_mgLevel2[ifa][0];
						int icv1 = cvofa_mgLevel2[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel2[icv0] + temp_mgLevel2[icv1]);
						mul_fa_mgLevel2[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
						lamOcp_fa_mgLevel2[ifa] = mul_fa_mgLevel2[ifa] / Pr;
						
					}
					
					// boundary faces computed next in setBC
				}
				else {
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
			}
		}
		
		if (iMesh == 3) {
			if (mu_ref > 0.0) {
				if (viscMode == "SUTHERLAND") {
					// internal faces
					for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
						
						int icv0 = cvofa_mgLevel3[ifa][0];
						int icv1 = cvofa_mgLevel3[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel3[icv0] + temp_mgLevel3[icv1]);
						mul_fa_mgLevel3[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
						lamOcp_fa_mgLevel3[ifa] = mul_fa_mgLevel3[ifa] / Pr;
					}
					// boundary faces computed next in setBC
				} 
				else if (viscMode == "POWERLAW") {
					// internal faces
					for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
						
						int icv0 = cvofa_mgLevel3[ifa][0];
						int icv1 = cvofa_mgLevel3[ifa][1];
						
						double temperature  = 0.5*(temp_mgLevel3[icv0] + temp_mgLevel3[icv1]);
						mul_fa_mgLevel3[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
						lamOcp_fa_mgLevel3[ifa] = mul_fa_mgLevel3[ifa] / Pr;
					}
					
					// boundary faces computed next in setBC
				}
				else {
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
			}
		}
	}
	
	virtual void ComputeBCProperties_T_mg(FaZone *zone, int iMesh)
	{
		if (iMesh == 1) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel1[ifa] = GAMMA;
				RoM_bfa_mgLevel1[ifa] = R_gas;
				h_bfa_mgLevel1[ifa] = GAMMA * R_gas / (GAMMA - 1.0) * T_bfa_mgLevel1[ifa];
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND") {    
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel1[ifa] = mu_ref*pow(T_bfa_mgLevel1[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel1[ifa] + SL_Sref);
				}
				else if (viscMode == "POWERLAW") {
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel1[ifa] = mu_ref * pow(T_bfa_mgLevel1[ifa]/T_ref, mu_power_law);
				}
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel1[ifa] = mul_fa_mgLevel1[ifa] / Pr;
			}
		}
		if (iMesh == 2) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel2[ifa] = GAMMA;
				RoM_bfa_mgLevel2[ifa] = R_gas;
				h_bfa_mgLevel2[ifa] = GAMMA * R_gas / (GAMMA - 1.0) * T_bfa_mgLevel2[ifa];
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND") {    
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel2[ifa] = mu_ref*pow(T_bfa_mgLevel2[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel2[ifa] + SL_Sref);
				}
				else if (viscMode == "POWERLAW") {
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel2[ifa] = mu_ref * pow(T_bfa_mgLevel2[ifa]/T_ref, mu_power_law);
				}
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel2[ifa] = mul_fa_mgLevel2[ifa] / Pr;
			}
		}
		if (iMesh == 3) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel3[ifa] = GAMMA;
				RoM_bfa_mgLevel3[ifa] = R_gas;
				h_bfa_mgLevel3[ifa] = GAMMA * R_gas / (GAMMA - 1.0) * T_bfa_mgLevel3[ifa];
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND") {    
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel3[ifa] = mu_ref*pow(T_bfa_mgLevel3[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel3[ifa] + SL_Sref);
				}
				else if (viscMode == "POWERLAW") {
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel3[ifa] = mu_ref * pow(T_bfa_mgLevel3[ifa]/T_ref, mu_power_law);
				}
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel3[ifa] = mul_fa_mgLevel3[ifa] / Pr;
			}
		}		
	}  
	
	virtual void ComputeBCProperties_H_mg(FaZone *zone, int iMesh) {
		
		if (iMesh == 1) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel1[ifa] = GAMMA;
				RoM_bfa_mgLevel1[ifa] = R_gas;
				T_bfa_mgLevel1[ifa] = h_bfa_mgLevel1[ifa] * (GAMMA - 1.0) / (GAMMA * R_gas);
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND")        
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel1[ifa] = mu_ref*pow(T_bfa_mgLevel1[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel1[ifa] + SL_Sref);
				//mul_fa[ifa] = 1.458e-6 * pow(T_bfa[ifa], 1.5) / (T_bfa[ifa] + 110.4);
				else if (viscMode == "POWERLAW")     
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel1[ifa] = mu_ref * pow(T_bfa_mgLevel1[ifa]/T_ref, mu_power_law);
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel1[ifa] = mul_fa_mgLevel1[ifa] / Pr;
			}
		}
		if (iMesh == 2) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel2[ifa] = GAMMA;
				RoM_bfa_mgLevel2[ifa] = R_gas;
				T_bfa_mgLevel2[ifa] = h_bfa_mgLevel2[ifa] * (GAMMA - 1.0) / (GAMMA * R_gas);
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND")        
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel2[ifa] = mu_ref*pow(T_bfa_mgLevel2[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel2[ifa] + SL_Sref);
				//mul_fa[ifa] = 1.458e-6 * pow(T_bfa[ifa], 1.5) / (T_bfa[ifa] + 110.4);
				else if (viscMode == "POWERLAW")     
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel2[ifa] = mu_ref * pow(T_bfa_mgLevel2[ifa]/T_ref, mu_power_law);
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel2[ifa] = mul_fa_mgLevel2[ifa] / Pr;
			}
		}
		if (iMesh == 3) {
			for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			{
				gam_bfa_mgLevel3[ifa] = GAMMA;
				RoM_bfa_mgLevel3[ifa] = R_gas;
				T_bfa_mgLevel3[ifa] = h_bfa_mgLevel3[ifa] * (GAMMA - 1.0) / (GAMMA * R_gas);
			}
			
			if (mu_ref > 0.0)
			{
				if (viscMode == "SUTHERLAND")        
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel3[ifa] = mu_ref*pow(T_bfa_mgLevel3[ifa]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(T_bfa_mgLevel3[ifa] + SL_Sref);
				//mul_fa[ifa] = 1.458e-6 * pow(T_bfa[ifa], 1.5) / (T_bfa[ifa] + 110.4);
				else if (viscMode == "POWERLAW")     
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						mul_fa_mgLevel3[ifa] = mu_ref * pow(T_bfa_mgLevel3[ifa]/T_ref, mu_power_law);
				else
				{
					cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
					throw(-1);
				}
				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
					lamOcp_fa_mgLevel3[ifa] = mul_fa_mgLevel3[ifa] / Pr;
			}
		}		
	}
};

#endif
