#ifndef JOEWITHMODELS_H
#define JOEWITHMODELS_H

#include "MiscUtils.h"
#include "tc_vec3d.h"
#include <iomanip>

#include "UgpWithCvCompFlow.h"


class JoeWithModels: virtual public UgpWithCvCompFlow
{
public:

  /**
   * constructor, pass ParamMap
   */
  JoeWithModels(ParamMap &p) : UgpWithCvCompFlow(p) {    init();  }

  JoeWithModels(ParamMap &p, int s) : UgpWithCvCompFlow(p) {    init();  }

  /**
   * constructor, pass name of joe's input file
   */
  JoeWithModels(char *name) : UgpWithCvCompFlow(name) {    init();  }

  JoeWithModels(char *name, int s) : UgpWithCvCompFlow(name) {    init();  }

  virtual void init()
  {
    if (mpi_rank == 0)
      cout << "JoeWithModels virtual init()"<< endl;
  }

  virtual ~JoeWithModels() {}

public:

  /**
   * check runtime and file
   */
  int doneSolver(const double residEnergy)
  {
    // ---------------------------------------------
    // returns 1 if we are done, 0 otherwise...
    // ---------------------------------------------
    int done = 0;

    if (mpi_rank==0)
    {

      // if the current completed step is greater or equal to the user
      // specified nsteps, then we are done...

      if ((nsteps >= 0) && (step >= nsteps))
      {
        cout << "step: " << step << endl;
        cout << "nsteps: " << nsteps << endl;
        cout << " > step from restart file is greater than NSTEP in Joe.in, you can reset the step by using RESET_STEP in Joe.in. stopping."<<endl;
        done = 1;
      }

      // check energy residual
      if (residEnergy < resid_energ_th)
      {
        cout<<" > reached threshold for energy residual. stopping."<<endl;
        done = 1;
      }

      // check if energy residual is NAN (not a number)
      if (residEnergy != residEnergy)
      {
        cout<<" > nan in energy residual detected. stopping."<<endl;
        done = 1;
      }

      // RUNTIME...
      if ((runtime_flag) && (MPI_Wtime() > runtime))
      {
        cout<<" > reached specified RUNTIME. stopping."<<endl;
        done = 1;
      }

      // killjoe file...
      if (fileExists("killjoe"))
      {
        MPI_File_delete("killjoe", MPI_INFO_NULL);
        cout<<" > Found file killjoe. stopping."<<endl;
        done = 1;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    return (done);
  }



  /**
   * hook for initializing Navier-Stokes equations
   */
  virtual void initialHook()
  {
    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;

    double rho_init, p_init, u_init[3] = {0.0, 0.0, 0.0};

    Param *p;
    if (getParam(p, "RHO_INITIAL"))       rho_init = p->getDouble(1);
    else                                  rho_init = rho_ref;

    if (getParam(p, "P_INITIAL"))         p_init = p->getDouble(1);
    else                                  p_init = p_ref;

    if (getParam(p, "U_INITIAL")) 
    {
      u_init[0] = p->getDouble(1);
      u_init[1] = p->getDouble(2);
      u_init[2] = p->getDouble(3);
    }

    if (!checkDataFlag(rho))
      for (int icv=0; icv<ncv; icv++)
        rho[icv] = rho_init;

    if (!checkDataFlag(rhou))
      for (int icv=0; icv<ncv; icv++)
        for (int i=0; i<3; i++)
          rhou[icv][i] = rho_init*u_init[i];

    if (!checkDataFlag(gamma)) {
      double gamma_init = GAMMA ;
      for (int icv=0; icv<ncv; icv++)
        gamma[icv] = gamma_init ;
    }

    if (!checkDataFlag(rhoE))
      for (int icv=0; icv<ncv; icv++)
        rhoE[icv] = p_init/(gamma[icv]-1.0)+ 0.5*rho_init*vecDotVec3d(u_init, u_init);


    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

//    writeData(0);
  }


  /**
   * boundary HOOK
   */
  virtual void boundaryHook(double *rho, double (*vel)[3], double *press, FaZone *zone)  { /*empty*/ }


  /**
   * source HOOK, rhs already contains fluxes, don't overwrite !!!!!
   * e.g.: gravity
   *
   *  for (int icv = 0; icv < ncv; icv++)
   *  {
   *    for (int i = 0; i < 3; i++)
   *    rhs_rhou[icv][i] += cv_volume[icv]*rho[icv]*gravity[i];
   *    rhs_rhoE[icv] += cv_volume[icv]*(gravity[0]*rhou[icv][0]+gravity[1]*rhou[icv][1]+gravity[2]*rhou[icv][2]);
   *  }
   */
  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5])  { /*empty*/ }

  /**
   * source HOOK, rhs already contains fluxes, don't overwrite !!!!!
   */
  virtual void sourceHookCoupled(double **rhs, double ***A, int nScal, int flagImplicit)  { /*empty*/ }


  /**
   * minimum stuff to run a compressible Navier-Stokes simulation
   */
  void run();

  /**
   * explicit backward Euler
   */
  void runForwardEuler();

  /**
   * explicit 3-step Runge-Kutta
   */
  void runRK();

  /**
   * global implicit fully-coupled relaxation 
   */
  void runBackwardEuler();

  /**
   * global implicit fully-coupled relaxation including scalars
   */
  void runBackwardEulerCoupled();

  /**
   * global implicit semi-coupled relaxation (only some scalars are solved fully coupled)
   */
  void runBackwardEulerSemiCoupled();

  /**
   * global implicit fully-coupled relaxation for unsteady calcs 
   */
  void runBDF2();
	
	/**
   * global implicit semi-coupled relaxation for unsteady calcs 
   */
  void runBDF2SemiCoupled();

  /**
   * set boundary conditions for Navier-Stokes and scalars 
   */
  void setBC();
  
  /**
   * calculate RHS, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcRhs(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit);
  
  /**
   * calculate Euler flux for both NSE and scalars, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcEulerFlux(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit);

  /**
   * calculate viscous flux for NSE only, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcViscousFluxNS(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit);
  
  /**
   * calculate RHS, routine used for fully coupled implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcRhsCoupled(double **rhs, double ***A, int nScal, int flagImplicit);
  
  /**
   * calculate both Euler and viscous fluxes for both NSE and scalars, routine used for fully coupled implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcFluxCoupled(double **rhs, double ***A, int nScal, int flagImplicit);
  
  /**
   * calculate RHS, routine used for semi-coupled implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcRhsSemiCoupled(double **rhs, double **rhs_rhoScal, double ***A, double ***AScal, int nScalCoupled, int nScalUncoupled, int flagImplicit);
  
  /**
   * calculate both Euler and viscous fluxes for both NSE and scalars, routine used for semi-coupled implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcFluxSemiCoupled(double **rhs, double **rhs_rhoScal, double ***A, double ***AScal, int nScalCoupled, int nScalUncoupled, int flagImplicit);
  
  /**
   * reshape RHS and separate between coupled and uncoupled scalars
   */
  void reshapeRHSSemiCoupled(double **rhs, double **rhsScal, double *Flux, int icv, double factor, int nScal);
 
  /**
   * reshape Jacobi matrix and separate between coupled and uncoupled scalars
   */
  void reshapeJacobianSemiCoupled(double ***A, double ***AScal, double **A0, double **A1, int noc0, int noc1, double factor, int nScal);
  
  /**
   * show residual
   */
  void showResidue(double *rhsResid);
  
  /*
   * add your own temporal output here
   */
  virtual void temporalHook() { /*empty*/ }

  /*
   * add your own finalization here
   */
  virtual void finalHook() { /*empty*/ }
	
	/*
   * Find the face between two control volumes
   */
	int findFace(int cv_first, int cv_second);

	/*
   * Find the linelet for the line preconditioning
   */
	void initializeLinelet();
	
	/**
   * global implicit fully-coupled relaxation with multigrid 
   */
	void runBackwardEulerMultigrid();
	
	/**
   * compute viscous fluxes for coarse multigrid levels
   */
	void calcViscousFluxNS_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh);
	
	/**
   * compute inviscid fluxes for coarse multigrid levels
   */
  void calcEulerFlux_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh);
	
	/**
   * time integration for the finest grid level
   */
	void calcTimeInt(double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double (*A)[5][5], double (*dq)[5], double (*rhs)[5], double underRelax, int maxIterLS, double zeroAbsLS, double zeroRelLS);
	
	/**
   * time integration for the coarse grid level
   */
	void calcTimeInt_mg(double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double (*A)[5][5], double (*dq)[5], double (*rhs)[5], 
											double underRelax, int maxIterLS, double zeroAbsLS, double zeroRelLS, int iMesh);
	
	/**
   * Compute Rhs of the multigrid method
   */
	void calcRhs_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh);
	
	/**
   * set boundary conditions for Navier-Stokes multigrid method
   */
	void setBC_mg(int iMesh);
	
	/**
   * Multigrid operators
   */
	void Set_SaveSolution(int index, int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetProlongated_Correction(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetResidual_Smoothing(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int nSmooth, double val_smooth_coeff, int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetSolution(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh);
	
	/**
   * Compute damping sensor for multigrid
   */
	void calcDampingSensors(void);
	
	/**
	 * Compute parameters for the JST method
	 */
	void calcJSTCoeff_Const(int iMesh);
	void calcJSTCoeff_Var(int iMesh);
	
	/**
	 * Adapt initial grid
	 */
	void AdaptInitialGrid(void);
	
	/**
	 * Adapt initial grid
	 */
	int FindEdgeCV(int ip_0, int ip_1, int **Edge, int *nbono_i, int *nbono_v);
	int FindFaceCV(int icv, int ip_0, int ip_1, int ip_3, int ip_4);
	int CheckTetCode(bool *AdaptCode);
	int CheckHexaCode(bool *AdaptCode);
	int CheckPyramCode(bool *AdaptCode);
	void TetDivision(int code, int *nodes, int **Division, int *nPart);
	void HexaDivision(int code, int *nodes, int **Division, int *nPart);
	void PyramDivision(int code, int *nodes, int **Division, int *nPart);


	
	/**
   * Multigrid operators
   */
	void SetResidual_Term(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetForcing_Term(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double *rhs_rho_mgLevel1, double (*rhs_rhou_mgLevel1)[3], 
											 double *rhs_rhoE_mgLevel1, int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetRestricted_Solution(int iMesh);
	
	/**
	 * Multigrid operators
	 */
	void SetUpdateGhost_Solution(int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetProjected_Solution(int iMesh);
	
	/**
   * Multigrid operators
   */
	void SetProjected_Residual(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double *rhs_rho_mgLevel1, double (*rhs_rhou_mgLevel1)[3], double *rhs_rhoE_mgLevel1);
	
	/*
   * Find the face between two control volumes
   */
	int findFace_mgLevel1(int cv_first, int cv_second);
	
	/*
   * Find the face between two control volumes
   */
	int findFace_mgLevel2(int cv_first, int cv_second);
	
	/*
   * Find the face between two control volumes
   */
	int findFace_mgLevel3(int cv_first, int cv_second);
	
	/*
   * Check the orientation of the control volume
   */
	double checkFaceOrientation(int ifa, int cv_first, int cv_second);
	
	/*
   * Check the orientation of the control volume
   */
	double checkFaceOrientation_mgLevel1(int ifa, int cv_first, int cv_second);
	
	/*
   * Check the orientation of the control volume
   */
	double checkFaceOrientation_mgLevel2(int ifa, int cv_first, int cv_second);
	
	/*
   * Create first coarse level from the initial grid
   */
	void initializeFromFineGrid();
	
	/*
   * Create first coarsest grid levels
   */
	void initializeFromMultiGrid(int iMesh);
	
	/*
   * Residual visualization in serial
   */
	void setHistoryFile();
	
	/*
   * Residual visualization in serial
   */
	void writeHistoryFile(double *rhsResid);
  
	virtual void boundaryHook_mg(double *rho, double (*vel)[3], double *press, FaZone *zone, int iMesh)  { /*empty*/ }

};


#endif



