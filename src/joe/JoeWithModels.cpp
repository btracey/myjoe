#include "JoeWithModels.h"
#include "myMem.h"

#define SHOCKFIX_ACROSS_SHOCK

void JoeWithModels::run()
{
	// read mesh or restart file
	initializeFromRestartFile(getStringParam("RESTART"));
	
	// adapt initial grid
	//	AdaptInitialGrid();
	
	if (checkParam("RESET_STEP"))
	{
	  step = 0;
	  if (mpi_rank == 0)
	    cout << "RESET_STEP: " << step << endl;
	}
	
	// initialize models
	initialHookScalarRansTurbModel();
	initialHookScalarRansCombModel();
	
	// initialize Navier-Stokes equations
	initialHook(); 
	
	// initialize Linele for preconditioning
	string SolverName = getStringParam("LINEAR_SOLVER_NS");
	string SemiName   = getStringParam("SEMICOARSENING","NO");
	if ((SolverName == "BCGSTAB_LINELET") || (SemiName == "YES")) initializeLinelet();
	
	string tIntName = getStringParam("TIME_INTEGRATION");
	
	if      (tIntName == "FORWARD_EULER")               runForwardEuler();
	else if (tIntName == "RK")                          runRK();
	else if (tIntName == "BACKWARD_EULER")              runBackwardEuler();
	else if (tIntName == "BACKWARD_EULER_MULTIGRID")    runBackwardEulerMultigrid();
	else if (tIntName == "BACKWARD_EULER_COUPLED")      runBackwardEulerCoupled();
	else if (tIntName == "BACKWARD_EULER_SEMICOUPLED")  runBackwardEulerSemiCoupled();
	else if (tIntName == "BDF2")                        runBDF2();
	else if (tIntName == "BDF2_SEMICOUPLED")		    runBDF2SemiCoupled();
	else
	  if (mpi_rank == 0)
	  {
	    cerr << "ERROR: wrong time integration scheme specified !" << endl;
	    cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2, BACKWARD_EULER_COUPLED, BACKWARD_EULER_SEMICOUPLED" << endl;
	  }

}

void JoeWithModels::runForwardEuler()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  double *Residual = new double[5+nScal];

  double *drho = new double[ncv];
  double (*drhou)[3] = new double[ncv][3];
  double *drhoE = new double[ncv];
  double **drhoScal = NULL;
  if (nScal > 0) getMem2D(&drhoScal, 0, nScal-1, 0, ncv-1, "drhoScal");

  double (*dummyA)[5][5] = NULL;       // provide dummy pointer for explicit mode!
  double ***dummyAScal   = NULL;       // provide dummy pointer for explicit mode!

  if (checkParam("WRITE_HISTORY")) setHistoryFile();

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // Loop over time steps
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  step = 0;   // set to zero, even if read from restart file, as for steady state time step doesn't matter

  int done = 0;
  if (nsteps == 0)    done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  while (done != 1)
  {
    step++;
    double dt = calcDt(cfl);
    time += dt;


    // =========================================================================================
    // EXPLICIT FORWARD STEP
    // =========================================================================================

    calcRhs(drho, drhou, drhoE, drhoScal, dummyA, dummyAScal, 0);

    for (int icv = 0; icv < ncv; icv++)
    {
      double tmp = local_dt[icv]/cv_volume[icv];
      rho[icv] += tmp * drho[icv];
      rhoE[icv] += tmp * drhoE[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] += tmp * drhou[icv][i];
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;

      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        phi[icv] = ((rho[icv] - tmp * drho[icv]) * phi[icv] + tmp * drhoScal[iScal][icv]) / rho[icv];
      }
      updateCvData(phi, REPLACE_DATA);
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(drho[icv]);
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(drhou[icv][i]);
      myResidual[4] += fabs(drhoE[icv]);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(drhoScal[iScal][icv]);

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt << " time: " << time << endl;

      showResidue(Residual);

      if (checkParam("WRITE_HISTORY"))
    	  writeHistoryFile(Residual);
    }

    temporalHook();

    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation

  }

  writeRestart();

  finalHookScalarRansTurbModel();
  finalHook();

  delete [] drho;
  delete [] drhou;
  delete [] drhoE;
  if (nScal > 0) freeMem2D(drhoScal, 0, nScal-1, 0, ncv-1);

  delete [] myResidual;
  delete [] Residual;
}

void JoeWithModels::runRK()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  double *Residual = new double[5+nScal];

  double *drho1 = new double[ncv];
  double (*drhou1)[3] = new double[ncv][3];
  double *drhoE1 = new double[ncv];
  double **drhoScal1 = NULL;
  if (nScal > 0) getMem2D(&drhoScal1, 0, nScal-1, 0, ncv-1, "drhoScal1");

  double *drho2 = new double[ncv];
  double (*drhou2)[3] = new double[ncv][3];
  double *drhoE2 = new double[ncv];
  double **drhoScal2 = NULL;
  if (nScal > 0) getMem2D(&drhoScal2, 0, nScal-1, 0, ncv-1, "drhoScal2");

  double *drho3 = new double[ncv];
  double (*drhou3)[3] = new double[ncv][3];
  double *drhoE3 = new double[ncv];
  double **drhoScal3 = NULL;
  if (nScal > 0) getMem2D(&drhoScal3, 0, nScal-1, 0, ncv-1, "drhoScal3");

  double *rho0 = new double[ncv];
  double (*rhou0)[3] = new double[ncv][3];
  double *rhoE0 = new double[ncv];
  double **rhoScal0 = NULL;
  if (nScal > 0) getMem2D(&rhoScal0, 0, nScal-1, 0, ncv-1, "rhoScal0");

  double (*dummyA)[5][5] = NULL;       // provide dummy pointer for explicit mode!
  double ***dummyAScal   = NULL;       // provide dummy pointer for explicit mode!

  if (checkParam("WRITE_HISTORY")) setHistoryFile();

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // Loop over time steps
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  while (done != 1)
  {
    step++;
    double dt = calcDt(cfl);
    time += dt;


    // =========================================================================================
    // copy the current solution into rho0, rhou0, etc...
    // =========================================================================================
    for (int icv = 0; icv < ncv; icv++)
    {
      rho0[icv] = rho[icv];
      rhoE0[icv] = rhoE[icv];
      for (int i = 0; i < 3; i++)
        rhou0[icv][i] = rhou[icv][i];
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        rhoScal0[iScal][icv] = rho[icv] * phi[icv];
    }


    // =========================================================================================
    // RK STEP 1
    // =========================================================================================
    calcRhs(drho1, drhou1, drhoE1, drhoScal1, dummyA, dummyAScal, 0);

    for (int icv = 0; icv < ncv; icv++)
    {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho1[icv] *= tmp;
      drhoE1[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou1[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal1[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++)
    {
      rho[icv] = rho0[icv] + 0.5*drho1[icv];
      rhoE[icv] = rhoE0[icv] + 0.5*drhoE1[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] + 0.5*drhou1[icv][i];
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = (rhoScal0[iScal][icv] + 0.5*drhoScal1[iScal][icv]) / rho[icv];
      updateCvData(phi, REPLACE_DATA);
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // RK STEP 2
    // =========================================================================================
    calcRhs(drho2, drhou2, drhoE2, drhoScal2, dummyA, dummyAScal, 0);

    for (int icv = 0; icv < ncv; icv++) {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho2[icv] *= tmp;
      drhoE2[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou2[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal2[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++) {
      rho[icv] = rho0[icv] - drho1[icv] + 2.0*drho2[icv];
      rhoE[icv] = rhoE0[icv] - drhoE1[icv] + 2.0*drhoE2[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] - drhou1[icv][i] + 2.0*drhou2[icv][i];
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = (rhoScal0[iScal][icv] - drhoScal1[iScal][icv] + 2.0*drhoScal2[iScal][icv]) / rho[icv];
      updateCvData(phi, REPLACE_DATA);
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // RK STEP 3
    // =========================================================================================
    calcRhs(drho3, drhou3, drhoE3, drhoScal3, dummyA, dummyAScal, 0);
    for (int icv = 0; icv < ncv; icv++) {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho3[icv] *= tmp;
      drhoE3[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou3[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal3[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++) {
      rho[icv] = rho0[icv] + (drho1[icv] + 4.0*drho2[icv] + drho3[icv])/6.0;
      rhoE[icv] = rhoE0[icv] + (drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv])/6.0;
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] + (drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i])/6.0;
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = (rhoScal0[iScal][icv] + (drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0) / rho[icv];
      updateCvData(phi, REPLACE_DATA);
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(drho1[icv] + 4.0*drho2[icv] + drho3[icv]) / 6.0;
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i]) / 6.0;
      myResidual[4] += fabs(drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv]) / 6.0;
      residField[icv] = fabs(drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv]) / 6.0;
    }

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv]
                                    + drhoScal3[iScal][icv])/6.0;

      switch (iScal)
      {
      case 0:
        for (int icv = 0; icv < ncv; icv++)
          residField0[icv] = fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv]
    			          + drhoScal3[iScal][icv])/6.0;
        break;
      case 1:
        for (int icv = 0; icv < ncv; icv++)
          residField1[icv] = fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv]
    	                         + drhoScal3[iScal][icv])/6.0;
        break;
      case 2:
        for (int icv = 0; icv < ncv; icv++)
          residField2[icv] = fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv]
    	                         + drhoScal3[iScal][icv])/6.0;
        break;
      case 3:
        for (int icv = 0; icv < ncv; icv++)
          residField3[icv] = fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv]
    	                         + drhoScal3[iScal][icv])/6.0;
        break;
      }
    }

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt << " time: " << time << endl;

      showResidue(Residual);

      if (checkParam("WRITE_HISTORY"))
    	  writeHistoryFile(Residual);
    }

    temporalHook();

    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation

  }

  writeRestart();

  finalHookScalarRansTurbModel();
  finalHook();

  delete [] drho1;
  delete [] drhou1;
  delete [] drhoE1;
  if (nScal > 0) freeMem2D(drhoScal1, 0, nScal-1, 0, ncv-1);

  delete [] drho2;
  delete [] drhou2;
  delete [] drhoE2;
  if (nScal > 0) freeMem2D(drhoScal2, 0, nScal-1, 0, ncv-1);

  delete [] drho3;
  delete [] drhou3;
  delete [] drhoE3;
  if (nScal > 0) freeMem2D(drhoScal3, 0, nScal-1, 0, ncv-1);

  delete [] rho0;
  delete [] rhou0;
  delete [] rhoE0;
  if (nScal > 0) freeMem2D(rhoScal0, 0, nScal-1, 0, ncv-1);

  delete [] myResidual;
  delete [] Residual;
}


void JoeWithModels::runBackwardEuler()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  double *Residual = new double[5+nScal];

  double *RHSrho = new double[ncv];
  double (*RHSrhou)[3] = new double[ncv][3];
  double *RHSrhoE = new double[ncv];

  double (*A)[5][5] = new double[nbocv_s][5][5];
  double (*dq)[5] = new double[ncv_g][5];
  double (*rhs)[5] = new double[ncv][5];

  double ***AScal  = NULL;  if (nScal > 0) getMem3D(&AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
  double **dScal   = NULL;  if (nScal > 0) getMem2D(&dScal,   0, nScal-1, 0, ncv_g-1, "dScal");
  double **rhsScal = NULL;  if (nScal > 0) getMem2D(&rhsScal, 0, nScal-1, 0, ncv-1, "rhsScal");

  //------------------------------------
  // some parameters
  //------------------------------------
  if (checkParam("WRITE_HISTORY")) setHistoryFile();

  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
  double underRelaxScalars = getDoubleParam("UNDER_RELAXATION_SCALARS", "0.2");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

	string SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");
	if (SpaceIntName == "JST") { 
		PressSensor = new double[ncv_g]; Conserv_Und_Lapl = new double[ncv_g][5]; Lambda = new double[ncv_g];
		BoundaryCV = new bool[ncv_g]; NeighborCV = new int[ncv_g];
		p1_Und_Lapl = new double[ncv_g]; p2_Und_Lapl   = new double[ncv_g]; 
		calcJSTCoeff_Const(0);
	}

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
  int done = doneSolver(1.0e20);

  if (initial_flowfield_output == "YES")     writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

    // ---------------------------------------------------------------------------------
    // Compute RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          A[noc][i][j] = 0.0;
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 5; i++)
        dq[icv][i] = 0.0;

    for (int iScal = 0; iScal < nScal; iScal++)                         // set AScal, dScal to zero! rhs is set zero in calcRHS
    {
      for (int i = 0; i <= 5; i++)
        for (int noc = 0; noc < nbocv_s; noc++)
          AScal[iScal][i][noc] = 0.0;
      for (int icv = 0; icv < ncv_g; icv++)
        dScal[iScal][icv] = 0.0;
    }


    // ---------------------------------------------------------------------------------
    // calculate RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------

    calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, true);

    // ---------------------------------------------------------------------------------
    // solve linear system for the NSE
    // ---------------------------------------------------------------------------------

    for (int icv=0; icv<ncv; ++icv)                                 // prepare rhs and A
    {
      rhs[icv][0] = underRelax*RHSrho[icv];
      rhs[icv][1] = underRelax*RHSrhou[icv][0];
      rhs[icv][2] = underRelax*RHSrhou[icv][1];
      rhs[icv][3] = underRelax*RHSrhou[icv][2];
      rhs[icv][4] = underRelax*RHSrhoE[icv];

      residField[icv] = RHSrhou[icv][0];

      double tmp = cv_volume[icv]/(local_dt[icv]);
      for (int i = 0; i < 5; i++)
        A[nbocv_i[icv]][i][i] += tmp;                               // diagonal part ( vol/dt + A )
    }

    solveCoupledLinSysNS(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS);  // solve linear system

    for (int icv=0; icv<ncv; icv++)                                 // update solution: Q_new = Q_old + Delta_Q
    {
      rho[icv]     += dq[icv][0];
      rhou[icv][0] += dq[icv][1];
      rhou[icv][1] += dq[icv][2];
      rhou[icv][2] += dq[icv][3];
      rhoE[icv]    += dq[icv][4];
    }

    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

    UpdateCvDataStateVec(dq);                                       // update dq since neighbors needed to compute RHS of scalars


    // ---------------------------------------------------------------------------------
    // solve linear system for the scalars
    // ---------------------------------------------------------------------------------

    // the scalars are solved separately from the NSE but in order to ensure consistency with
    // the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou)
    // are used on the RHS of the equations. This means that AScal[iScal][4] is the only implicit
    // term on the LHS, while AScal[iScal][0-3] are put back to the right side.

    for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
    {
      string scalname = scalarTranspEqVector[iScal].getName();
      for (int icv = 0; icv < ncv; ++icv)
      {
        rhsScal[iScal][icv] *= underRelaxScalars ;
              
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;

        double tmp = cv_volume[icv]/(local_dt[icv]);
        AScal[iScal][5][noc_f] += tmp;                                 // diagonal part ( vol/dt + A )        
        

        // move the other implicit terms to the RHS
        for (int noc = noc_f; noc <= noc_l; noc++)
          rhsScal[iScal][icv] = rhsScal[iScal][icv]
                              - AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
                              - AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
                              - AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
                              - AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
                              - AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];
      }
      
      switch (iScal)
      {
      case 0: for (int icv = 0; icv < ncv; icv++) residField0[icv] = rhsScal[iScal][icv]; break;
      case 1: for (int icv = 0; icv < ncv; icv++) residField1[icv] = rhsScal[iScal][icv]; break;
      case 2: for (int icv = 0; icv < ncv; icv++) residField2[icv] = rhsScal[iScal][icv]; break;
      case 3: for (int icv = 0; icv < ncv; icv++) residField3[icv] = rhsScal[iScal][icv]; break;
      }

      solveLinSysScalar(dScal[iScal], AScal[iScal][5], rhsScal[iScal],
                        scalarTranspEqVector[iScal].phiZero,
                        scalarTranspEqVector[iScal].phiZeroRel,
                        scalarTranspEqVector[iScal].phiMaxiter,
                        scalarTranspEqVector[iScal].getName());

      // update scalars and clip
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScal][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
      updateCvData(phi, REPLACE_DATA);
    }



    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }
	  
    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(RHSrho[icv]);
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(RHSrhou[icv][i]);
      myResidual[4] += fabs(RHSrhoE[icv]);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(rhsScal[iScal][icv]/underRelax);

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);

      if (checkParam("WRITE_HISTORY"))
    	  writeHistoryFile(Residual);
    }

    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  //temporalHook();
  finalHookScalarRansTurbModel();
  finalHook();

  writeData(step);
  writeRestart();


  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------

  delete [] A;
  delete [] rhs;

  delete [] RHSrho;
  delete [] RHSrhou;
  delete [] RHSrhoE;

  delete [] dq;

  if (nScal > 0) freeMem3D(AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1);
  if (nScal > 0) freeMem2D(dScal,   0, nScal-1, 0, ncv_g-1);
  if (nScal > 0) freeMem2D(rhsScal, 0, nScal-1, 0, ncv-1);
}


void JoeWithModels::runBDF2()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  double *Residual = new double[5+nScal];
  double *fstResidual = new double[5+nScal];
  double *relResidual = new double[5+nScal];

  double *RHSrho = new double[ncv];
  double (*RHSrhou)[3] = new double[ncv][3];
  double *RHSrhoE = new double[ncv];

  double (*A)[5][5] = new double[nbocv_s][5][5];
  double (*dq)[5] = new double[ncv_g][5];
  double (*rhs)[5] = new double[ncv][5];
  double (*qn)[5] = new double[ncv][5];
  double (*qnm1)[5] = new double[ncv][5];

  double ***AScal   = NULL;  if (nScal > 0) getMem3D(&AScal,    0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
  double **dScal    = NULL;  if (nScal > 0) getMem2D(&dScal,    0, nScal-1, 0, ncv_g-1, "dScal");
  double **rhsScal  = NULL;  if (nScal > 0) getMem2D(&rhsScal,  0, nScal-1, 0, ncv-1, "rhsScal");
  double **qnScal   = NULL;  if (nScal > 0) getMem2D(&qnScal,   0, nScal-1, 0, ncv-1, "qnScal");
  double **qnm1Scal = NULL;  if (nScal > 0) getMem2D(&qnm1Scal, 0, nScal-1, 0, ncv-1, "qnm1Scal");

  //------------------------------------
  // some parameters
  //------------------------------------
  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.9");
  int mnNewton = getIntParam("N_NEWTON_ITER", "4");
  int bdf2_TVD_limiter = getIntParam("BDF2_TVD_LIMITER", "1");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");


  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------

  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();



  //double bdf2Alfa = 1.5, bdf2Beta = 2.0, bdf2Gamma = 0.5;
  double bdf2Alfa = 0.5, bdf2Beta = 1.0, bdf2Gamma = 0.5;

  // store n and n-1 and solve n+1 using second order
  for (int icv=0; icv<ncv; icv++)
  {
    qn[icv][0] = rho[icv];
    qn[icv][1] = rhou[icv][0];
    qn[icv][2] = rhou[icv][1];
    qn[icv][3] = rhou[icv][2];
    qn[icv][4] = rhoE[icv];
  }

  for (int iScal = 0; iScal < nScal; iScal++)
  {
    double *phi = scalarTranspEqVector[iScal].phi;
    for (int icv=0; icv<ncv; icv++)
      qnScal[iScal][icv] = rho[icv]*phi[icv];
  }

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);
    time += dt_min;


    // store n and n-1 and solve n+1 using second order
    for (int icv=0; icv<ncv; icv++)
    {
      qnm1[icv][0] = qn[icv][0];
      qnm1[icv][1] = qn[icv][1];
      qnm1[icv][2] = qn[icv][2];
      qnm1[icv][3] = qn[icv][3];
      qnm1[icv][4] = qn[icv][4];

      qn[icv][0] = rho[icv];
      qn[icv][1] = rhou[icv][0];
      qn[icv][2] = rhou[icv][1];
      qn[icv][3] = rhou[icv][2];
      qn[icv][4] = rhoE[icv];
    }

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;

      for (int icv=0; icv<ncv; icv++)
      {
        qnm1Scal[iScal][icv] = qnScal[iScal][icv];
        qnScal[iScal][icv] = rho[icv]*phi[icv];
      }
    }


    // ---------------------------------------------------------------------------------
    // start Newton iterations
    // ---------------------------------------------------------------------------------


    if ((mpi_rank == 0) && (step%check_interval == 0))
      printf("dt: %.6le\n", dt_min);

    relResidual[4] = 1.0e20;

    int pN = 0;
    while ((pN < mnNewton) && (relResidual[4] > 1.0e-6))              // newton steps!!!
    {
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] = 0.0;

      for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 5; j++)
            A[noc][i][j] = 0.0;
			
      for (int icv = 0; icv < ncv_g; icv++)
        for (int i = 0; i < 5; i++)
          dq[icv][i] = 0.0;

      for (int iScal = 0; iScal < nScal; iScal++)                         // set AScal, dScal to zero! rhs is set zero in calcRHS
      {
        for (int i = 0; i <= 5; i++)
          for (int noc = 0; noc < nbocv_s; noc++)
            AScal[iScal][i][noc] = 0.0;
        for (int icv = 0; icv < ncv_g; icv++)
          dScal[iScal][icv] = 0.0;
      }


      // ---------------------------------------------------------------------------------
      // calculate RHS for both NSE and scalars
      // ---------------------------------------------------------------------------------

      calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, true);

      // ---------------------------------------------------------------------------------
      // solve linear system for the NSE
      // ---------------------------------------------------------------------------------

      for (int icv=0; icv<ncv; ++icv)                                 // prepare rhs and A
      {
        double tmp = cv_volume[icv]/(local_dt[icv]);

        double psi[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
        //double psi[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        if (bdf2_TVD_limiter == 1)
        {
          double tvdR[5], minmod[5];

          tvdR[0] = (rho[icv]-    qn[icv][0])/(qn[icv][0]-qnm1[icv][0]);
          tvdR[1] = (rhou[icv][0]-qn[icv][1])/(qn[icv][1]-qnm1[icv][1]);
          tvdR[2] = (rhou[icv][1]-qn[icv][2])/(qn[icv][2]-qnm1[icv][2]);
          tvdR[3] = (rhou[icv][2]-qn[icv][3])/(qn[icv][3]-qnm1[icv][3]);
          tvdR[4] = (rhoE[icv]-   qn[icv][4])/(qn[icv][4]-qnm1[icv][4]);

          for (int i=0; i<5; i++)
          {
            if (tvdR[i] > 0.0)
            {
              if (1.0 < fabs(tvdR[i]))    minmod[i] = 1.0;
              else                        minmod[i] = tvdR[i];
            }
            else                          minmod[i] = 0.0;

            psi[i] = pow(minmod[i]/max(1.0, fabs(tvdR[i])), 0.5);
          }
        }

        rhs[icv][0] = underRelax*(RHSrho[icv]     - ((1.0+psi[0]*bdf2Alfa)*rho[icv]     - (1.0+psi[0]*bdf2Beta)*qn[icv][0] + psi[0]*bdf2Gamma*qnm1[icv][0])*tmp);
        rhs[icv][1] = underRelax*(RHSrhou[icv][0] - ((1.0+psi[1]*bdf2Alfa)*rhou[icv][0] - (1.0+psi[1]*bdf2Beta)*qn[icv][1] + psi[1]*bdf2Gamma*qnm1[icv][1])*tmp);
        rhs[icv][2] = underRelax*(RHSrhou[icv][1] - ((1.0+psi[2]*bdf2Alfa)*rhou[icv][1] - (1.0+psi[2]*bdf2Beta)*qn[icv][2] + psi[2]*bdf2Gamma*qnm1[icv][2])*tmp);
        rhs[icv][3] = underRelax*(RHSrhou[icv][2] - ((1.0+psi[3]*bdf2Alfa)*rhou[icv][2] - (1.0+psi[3]*bdf2Beta)*qn[icv][3] + psi[3]*bdf2Gamma*qnm1[icv][3])*tmp);
        rhs[icv][4] = underRelax*(RHSrhoE[icv]    - ((1.0+psi[4]*bdf2Alfa)*rhoE[icv]    - (1.0+psi[4]*bdf2Beta)*qn[icv][4] + psi[4]*bdf2Gamma*qnm1[icv][4])*tmp);

        residField[icv] = RHSrhoE[icv];

        for (int i=0; i<5; i++)
          myResidual[i] += fabs(rhs[icv][i]);

        for (int i = 0; i < 5; i++)
          A[nbocv_i[icv]][i][i] += (1.0+psi[0]*bdf2Alfa)*tmp;             // diagonal part ( vol/dt + A )
      }

      solveCoupledLinSysNS(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS);  // solve linear system

      for (int icv=0; icv<ncv; icv++)                                     // update solution: Q_new = Q_old + Delta_Q
      {
        rho[icv]     += dq[icv][0];
        rhou[icv][0] += dq[icv][1];
        rhou[icv][1] += dq[icv][2];
        rhou[icv][2] += dq[icv][3];
        rhoE[icv]    += dq[icv][4];
      }

      updateCvData(rho,  REPLACE_DATA);
      updateCvData(rhou, REPLACE_ROTATE_DATA);
      updateCvData(rhoE, REPLACE_DATA);

      UpdateCvDataStateVec(dq);                                           // update dq since neighbors needed to compute RHS of scalars


      // ---------------------------------------------------------------------------------
      // solve linear system for the scalars
      // ---------------------------------------------------------------------------------

      // the scalars are solved separately from the NSE but in order to ensure consistency with
      // the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou, rhoE)
      // are used on the RHS of the equations. This means that AScal[iScal][5] is the only implicit
      // term on the LHS, while AScal[iScal][0-4] are put back to the right side.

      for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
      {
        double *phi = scalarTranspEqVector[iScal].phi;

        for (int icv = 0; icv < ncv; ++icv)
        {
          double tmp = cv_volume[icv]/(local_dt[icv]);
          double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];

          double psi = 1.0;
          if (bdf2_TVD_limiter == 1)
          {
            double minmod, tvdR = (rhoPhi-qnScal[iScal][icv])/(qnScal[iScal][icv]-qnm1Scal[iScal][icv]);
            if (tvdR > 0.0)
            {
              if (1.0 < fabs(tvdR))    minmod= 1.0;
              else                     minmod = tvdR;
            }
            else                       minmod = 0.0;

            psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5);
          }

          rhsScal[iScal][icv] = underRelax*( rhsScal[iScal][icv]
                                   -(  (1.0+psi*bdf2Alfa)*rhoPhi
                                     - (1.0+psi*bdf2Beta)*qnScal[iScal][icv]
                                         + psi*bdf2Gamma*qnm1Scal[iScal][icv])*tmp);

          int noc_f = nbocv_i[icv];
          int noc_l = nbocv_i[icv + 1] - 1;

          AScal[iScal][5][noc_f] += (1.0+psi*bdf2Alfa)*tmp;                                 // diagonal part ( vol/dt + A )

          // move the other implicit terms to the RHS
          for (int noc = noc_f; noc <= noc_l; noc++)
            rhsScal[iScal][icv] = rhsScal[iScal][icv]
                                - AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
                                - AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
                                - AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
                                - AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
                                - AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];

          myResidual[5+iScal] += fabs(rhsScal[iScal][icv]);
        }

        solveLinSysScalar(dScal[iScal], AScal[iScal][5], rhsScal[iScal],
                          scalarTranspEqVector[iScal].phiZero,
                          scalarTranspEqVector[iScal].phiZeroRel,
                          scalarTranspEqVector[iScal].phiMaxiter,
                          scalarTranspEqVector[iScal].getName());

        // update scalars and clip
        for (int icv = 0; icv < ncv; icv++)
          phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScal][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);

        updateCvData(phi, REPLACE_DATA);
      }

      calcStateVariables();

      calcMaterialProperties();

      setBC();

      calcRansTurbViscMuet();

      MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

      if (mpi_rank == 0)
      {
        if (pN == 0)
          for (int i=0; i<5+nScal; i++)
            fstResidual[i] = Residual[i];

        for (int i=0; i<5+nScal; i++)
          relResidual[i] = Residual[i]/(fstResidual[i]+1.0e-10);

        if (step%check_interval == 0)
        {
          printf("Newton step residual: %d:\t", pN+1);
          for (int i=0; i<5+nScal; i++)
            printf("%.4le\t", relResidual[i]);
          printf("\n");
        }
      }

      pN++;
    }

    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);
    }

    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  temporalHook();
  finalHookScalarRansTurbModel();
  finalHook();

  writeRestart();


  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------
  delete [] RHSrho;
  delete [] RHSrhou;
  delete [] RHSrhoE;

  delete [] A;
  delete [] dq;
  delete [] rhs;
  delete [] qn;
  delete [] qnm1;

  if (nScal > 0) freeMem3D(AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1);
  if (nScal > 0) freeMem2D(dScal,   0, nScal-1, 0, ncv_g-1);
  if (nScal > 0) freeMem2D(rhsScal, 0, nScal-1, 0, ncv-1);
  if (nScal > 0) freeMem2D(qnScal, 0, nScal-1, 0, ncv-1);
  if (nScal > 0) freeMem2D(qnm1Scal, 0, nScal-1, 0, ncv-1);

}

void JoeWithModels::runBDF2SemiCoupled()
{
  int nScal = scalarTranspEqVector.size();
  int nScalCoupled = 0;
  int nScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      nScalCoupled++;
    else
      nScalUncoupled++;
  }
	
  // All: NSE + scalars
  double *myResidual = new double[5+nScal];
  double *Residual   = new double[5+nScal];
	double *fstResidual = new double[5+nScal];
  double *relResidual = new double[5+nScal];
	
  // NSE + coupled scalars
  double ***A;           getMem3D(&A,   0, nbocv_s-1, 0, 5+nScalCoupled-1, 0, 5+nScalCoupled-1, "JoeWithModels::runBackwardEulerCoupled -> A",   true);
  double **rhs;          getMem2D(&rhs, 0, ncv-1,     0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> rhs", true);
  double **dq;           getMem2D(&dq,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> dq",  true);
	double **qn;           getMem2D(&qn,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> qn",  true);
  double **qnm1;         getMem2D(&qnm1,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> qnm1",  true);
	
  // Uncoupled scalars
  double ***AScal  = NULL;  if (nScalUncoupled > 0) getMem3D(&AScal,   0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1, "AScal", true);
  double **rhsScal = NULL;  if (nScalUncoupled > 0) getMem2D(&rhsScal, 0, nScalUncoupled-1, 0, ncv-1, "rhsScal", true);
  double **dScal   = NULL;  if (nScalUncoupled > 0) getMem2D(&dScal,   0, nScalUncoupled-1, 0, ncv_g-1, "dScal", true);
	double **qnScal   = NULL;  if (nScalUncoupled > 0) getMem2D(&qnScal,   0, nScalUncoupled-1, 0, ncv-1, "qnScal", true);
  double **qnm1Scal = NULL;  if (nScalUncoupled > 0) getMem2D(&qnm1Scal, 0, nScalUncoupled-1, 0, ncv-1, "qnm1Scal", true);
	
  //------------------------------------
  // some parameters
  //------------------------------------
  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
	int mnNewton = getIntParam("N_NEWTON_ITER", "4");
	int bdf2_TVD_limiter = getIntParam("BDF2_TVD_LIMITER", "1");
	
  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
			" to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");
	
  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }
	
  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");
	
  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();
	
  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();
	
  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();
	
  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();
	
	
  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;
	
  if (initial_flowfield_output == "YES") writeData(0);
	
  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();
	
	double bdf2Alfa = 0.5, bdf2Beta = 1.0, bdf2Gamma = 0.5;
	
  // store n and n-1 and solve n+1 using second order
  for (int icv=0; icv<ncv; icv++)
  {
    qn[icv][0] = rho[icv];
    qn[icv][1] = rhou[icv][0];
    qn[icv][2] = rhou[icv][1];
    qn[icv][3] = rhou[icv][2];
    qn[icv][4] = rhoE[icv];
  }
	
	int iScalCoupled = 0; int iScalUncoupled = 0;
	for (int iScal = 0; iScal < nScal; iScal++)	{
		double *phi = scalarTranspEqVector[iScal].phi;
		if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
			for (int icv = 0; icv < ncv; icv++)
				qn[icv][5 + iScalCoupled] = rho[icv]*phi[icv];
			iScalCoupled++;
		}
		else {
			for (int icv = 0; icv < ncv; icv++)
				qnScal[iScalUncoupled][icv] = rho[icv]*phi[icv];
			iScalUncoupled++;	
		}
	}
	
  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);
		time += dt_min;
		
		
		// store n and n-1 and solve n+1 using second order
    for (int icv=0; icv<ncv; icv++) {
      qnm1[icv][0] = qn[icv][0]; qn[icv][0] = rho[icv];
      qnm1[icv][1] = qn[icv][1]; qn[icv][1] = rhou[icv][0];
      qnm1[icv][2] = qn[icv][2]; qn[icv][2] = rhou[icv][1];
      qnm1[icv][3] = qn[icv][3]; qn[icv][3] = rhou[icv][2];
      qnm1[icv][4] = qn[icv][4]; qn[icv][4] = rhoE[icv];
	  }
		
		int iScalCoupled = 0; int iScalUncoupled = 0; int icv;
		for (int iScal = 0; iScal < nScal; iScal++) {
			double *phi = scalarTranspEqVector[iScal].phi;
			if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
				for (icv = 0; icv < ncv; ++icv)
					qnm1[icv][5+iScalCoupled] = qn[icv][5+iScalCoupled]; qn[icv][5+iScalCoupled] = rho[icv]*phi[icv];
				iScalCoupled++;
			}
			else {
				for (icv=0; icv<ncv; icv++)
					qnm1Scal[iScalUncoupled][icv] = qnScal[iScalUncoupled][icv]; qnScal[iScalUncoupled][icv] = rho[icv]*phi[icv];
				iScalUncoupled++;
			}
		}
		
		// ---------------------------------------------------------------------------------
    // start Newton iterations
    // ---------------------------------------------------------------------------------
		
		
    if ((mpi_rank == 0) && (step%check_interval == 0))
      printf("dt: %.6le\n", dt_min);
		
    relResidual[4] = 1.0e20;
		
    int pN = 0;
    while ((pN < mnNewton) && (relResidual[4] > 1.0e-6))              // newton steps!!!
    {
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] = 0.0;
			
			
			for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
				for (int i = 0; i < 5+nScalCoupled; i++)
					for (int j = 0; j < 5+nScalCoupled; j++)
						A[noc][i][j] = 0.0;
			
			for (int icv = 0; icv < ncv_g; icv++)
				for (int i = 0; i < 5+nScalCoupled; i++)
					dq[icv][i] = 0.0;
			
			for (int iScal = 0; iScal < nScalUncoupled; iScal++)                // set AScal, dScal to zero! rhs is set zero in calcRHS
			{
				for (int i = 0; i <= 5; i++)
					for (int noc = 0; noc < nbocv_s; noc++)
						AScal[iScal][i][noc] = 0.0;
				for (int icv = 0; icv < ncv_g; icv++)
					dScal[iScal][icv] = 0.0;
			}
			
			// ---------------------------------------------------------------------------------
			// calculate RHS for both NSE and scalars
			// ---------------------------------------------------------------------------------
			
			calcRhsSemiCoupled(rhs, rhsScal, A, AScal, nScalCoupled, nScalUncoupled, true);
			
			// ---------------------------------------------------------------------------------
			// solve linear system for the NSE
			// ---------------------------------------------------------------------------------
			
			
			for (int icv = 0; icv < ncv; icv++)                                     // prepare rhs and A
			{
				
				double tmp = cv_volume[icv]/(local_dt[icv]);
				
				double psi[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
				
				if (bdf2_TVD_limiter == 1) {
					double tvdR[5], minmod[5];
					tvdR[0] = (rho[icv]-    qn[icv][0])/(qn[icv][0]-qnm1[icv][0]);
					tvdR[1] = (rhou[icv][0]-qn[icv][1])/(qn[icv][1]-qnm1[icv][1]);
					tvdR[2] = (rhou[icv][1]-qn[icv][2])/(qn[icv][2]-qnm1[icv][2]);
					tvdR[3] = (rhou[icv][2]-qn[icv][3])/(qn[icv][3]-qnm1[icv][3]);
					tvdR[4] = (rhoE[icv]-   qn[icv][4])/(qn[icv][4]-qnm1[icv][4]);
					
					for (int i=0; i<5; i++) {
						if (tvdR[i] > 0.0) {
							if (1.0 < fabs(tvdR[i]))    minmod[i] = 1.0;
							else                        minmod[i] = tvdR[i];
						}
						else                          minmod[i] = 0.0;
						psi[i] = pow(minmod[i]/max(1.0, fabs(tvdR[i])), 0.5);
					}
				}
				
				rhs[icv][0] = underRelax*(rhs[icv][0] - ((1.0+psi[0]*bdf2Alfa)*rho[icv]     - (1.0+psi[0]*bdf2Beta)*qn[icv][0] + psi[0]*bdf2Gamma*qnm1[icv][0])*tmp);
				rhs[icv][1] = underRelax*(rhs[icv][1] - ((1.0+psi[1]*bdf2Alfa)*rhou[icv][0] - (1.0+psi[1]*bdf2Beta)*qn[icv][1] + psi[1]*bdf2Gamma*qnm1[icv][1])*tmp);
				rhs[icv][2] = underRelax*(rhs[icv][2] - ((1.0+psi[2]*bdf2Alfa)*rhou[icv][1] - (1.0+psi[2]*bdf2Beta)*qn[icv][2] + psi[2]*bdf2Gamma*qnm1[icv][2])*tmp);
				rhs[icv][3] = underRelax*(rhs[icv][3] - ((1.0+psi[3]*bdf2Alfa)*rhou[icv][2] - (1.0+psi[3]*bdf2Beta)*qn[icv][3] + psi[3]*bdf2Gamma*qnm1[icv][3])*tmp);
				rhs[icv][4] = underRelax*(rhs[icv][4] - ((1.0+psi[4]*bdf2Alfa)*rhoE[icv]    - (1.0+psi[4]*bdf2Beta)*qn[icv][4] + psi[4]*bdf2Gamma*qnm1[icv][4])*tmp);
				
				residField[icv] = rhs[icv][4];
				
				for (int i=0; i<5; i++)
					myResidual[i] += fabs(rhs[icv][i]);
				
				for (int i = 0; i < 5; i++)
					A[nbocv_i[icv]][i][i] += (1.0+psi[0]*bdf2Alfa)*tmp;                                       // diagonal part ( vol/dt + A )
			}  
			
			// Add the 2nd order derivative contgribution to the residual for coupled scalars
			int iScalCoupled = 0;
			for (int iScal = 0; iScal < nScal; iScal++)  
				if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
					double *phi = scalarTranspEqVector[iScal].phi;
					
					for (int icv = 0; icv < ncv; ++icv) {
						double tmp = cv_volume[icv]/(local_dt[icv]);
						double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];
						
						double psi = 1.0;
						if (bdf2_TVD_limiter == 1) {
							double minmod, tvdR = (rhoPhi-qn[icv][5+iScalCoupled])/(qn[icv][5+iScalCoupled]-qnm1[icv][5+iScalCoupled]);
							if (tvdR > 0.0) {
								if (1.0 < fabs(tvdR))    minmod= 1.0;
								else                     minmod = tvdR;
							}
							else                       minmod = 0.0;
							psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5);
						}
						
						rhs[5+iScalCoupled][icv] = underRelax*( rhs[icv][5 + iScalCoupled] - ((1.0+psi*bdf2Alfa)*phi[icv] - (1.0+psi*bdf2Beta)*qn[icv][5+iScalCoupled] + psi*bdf2Gamma*qnm1[icv][5+iScalCoupled])*tmp);
						
						myResidual[5+iScalCoupled] += fabs(rhs[icv][5+iScalCoupled]);
						
						A[nbocv_i[icv]][5+iScalCoupled][5+iScalCoupled] += (1.0+psi*bdf2Alfa)*tmp;             // diagonal part ( vol/dt + A )
						
					}
					iScalCoupled++;
				}
			
			// solve linear system    
			solveCoupledLinSysNSCoupled(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, nScalCoupled); 
			
			for (int icv = 0; icv < ncv; icv++)                                     // update solution for NSE: Q_new = Q_old + Delta_Q
			{
				rho[icv]     += dq[icv][0];
				rhou[icv][0] += dq[icv][1];
				rhou[icv][1] += dq[icv][2];
				rhou[icv][2] += dq[icv][3];
				rhoE[icv]    += dq[icv][4];
			}
			
			updateCvData(rho,  REPLACE_DATA);
			updateCvData(rhou, REPLACE_ROTATE_DATA);
			updateCvData(rhoE, REPLACE_DATA);
			
			UpdateCvDataStateVecScal(dq, nScalCoupled);                             // update dq since neighbors needed to compute RHS of scalars
			
			
			iScalCoupled = 0;
			for (int iScal = 0; iScal < nScal; iScal++)                             // update + clip solution for coupled scalars: Q_new = (rho_old * Q_old + Delta_Q ) / rho_new
			{
				if (scalarTranspEqVector[iScal].coupling == "COUPLED")
				{
					double *phi = scalarTranspEqVector[iScal].phi;
					for (int icv = 0; icv < ncv; icv++)
						phi[icv] = min(max((phi[icv] * (rho[icv] - dq[icv][0]) + dq[icv][5+iScalCoupled]) / rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
					updateCvData(phi, REPLACE_DATA);
					iScalCoupled++;
				}
			}
			
			
			// ---------------------------------------------------------------------------------
			// solve linear system for the uncoupled scalars
			// ---------------------------------------------------------------------------------
			
			// the uncoupled scalars are solved separately from the NSE but in order to ensure consistency with
			// the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou, rhoE)
			// are used on the RHS of the equations. This means that AScal[iScal][5] is the only implicit
			// term on the LHS, while AScal[iScal][0-4] are put back to the right side.
			
			int iScalUncoupled = 0;
			for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
				if (scalarTranspEqVector[iScal].coupling != "COUPLED")
				{
					double *phi = scalarTranspEqVector[iScal].phi;
					for (int icv = 0; icv < ncv; ++icv)
					{
						double tmp = cv_volume[icv]/(local_dt[icv]);
						double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];
						
						double psi = 1.0;
						if (bdf2_TVD_limiter == 1)
						{
							double minmod, tvdR = (rhoPhi-qnScal[iScalUncoupled][icv])/(qnScal[iScalUncoupled][icv]-qnm1Scal[iScalUncoupled][icv]);
							if (tvdR > 0.0)
							{
								if (1.0 < fabs(tvdR))    minmod= 1.0;
								else                     minmod = tvdR;
							}
							else                       minmod = 0.0;
							
							psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5);
						}
						
						rhsScal[iScalUncoupled][icv] = underRelax*( rhsScal[iScalUncoupled][icv]
																											 -(  (1.0+psi*bdf2Alfa)*rhoPhi
																												 - (1.0+psi*bdf2Beta)*qnScal[iScalUncoupled][icv]
																												 + psi*bdf2Gamma*qnm1Scal[iScalUncoupled][icv])*tmp);
						
						int noc_f = nbocv_i[icv];
						int noc_l = nbocv_i[icv + 1] - 1;
						
						AScal[iScalUncoupled][5][noc_f] += (1.0+psi*bdf2Alfa)*tmp;                                 // diagonal part ( vol/dt + A )
						
						
						// move the other implicit terms to the RHS
						for (int noc = noc_f; noc <= noc_l; noc++)
							rhsScal[iScalUncoupled][icv] = rhsScal[iScalUncoupled][icv]
							- AScal[iScalUncoupled][0][noc] * dq[nbocv_v[noc]][0]
							- AScal[iScalUncoupled][1][noc] * dq[nbocv_v[noc]][1]
							- AScal[iScalUncoupled][2][noc] * dq[nbocv_v[noc]][2]
							- AScal[iScalUncoupled][3][noc] * dq[nbocv_v[noc]][3]
							- AScal[iScalUncoupled][4][noc] * dq[nbocv_v[noc]][4];
						
						myResidual[5+nScalCoupled+iScalUncoupled] += fabs(rhsScal[iScalUncoupled][icv]);
						
					}
					
					solveLinSysScalar(dScal[iScalUncoupled], AScal[iScalUncoupled][5], rhsScal[iScalUncoupled],
														scalarTranspEqVector[iScal].phiZero,
														scalarTranspEqVector[iScal].phiZeroRel,
														scalarTranspEqVector[iScal].phiMaxiter,
														scalarTranspEqVector[iScal].getName());
					
					
					// update scalars and clip
					for (int icv = 0; icv < ncv; icv++)
						phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScalUncoupled][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
					
					updateCvData(phi, REPLACE_DATA);
					
					iScalUncoupled++;
					
				}
			
			calcStateVariables();
			
			calcMaterialProperties();
			
			setBC();
			
			calcRansTurbViscMuet();
			
			MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			if (mpi_rank == 0)
			{
				if (pN == 0)
					for (int i=0; i<5+nScal; i++)
						fstResidual[i] = Residual[i];
				
				for (int i=0; i<5+nScal; i++)
					relResidual[i] = Residual[i]/(fstResidual[i]+1.0e-10);
				
				if (step%check_interval == 0)
				{
					printf("Newton step residual: %d:\t", pN+1);
					for (int i=0; i<5+nScal; i++)
						printf("%.4le\t", relResidual[i]);
					printf("\n");
				}
			}
			
			pN++;
		}
		
		// =========================================================================================
		// calculate and show residual
		// =========================================================================================
		if (step%check_interval == 0)
		{
			if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
				cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;
			
			showResidue(Residual);
		}
		
		temporalHook();
		dumpProbes(step, 0.0);
		writeData(step);
		
		if ((write_restart > 0) && (step % write_restart == 0))
			writeRestart(step);
		
		done = doneSolver(Residual[4]);   // pass the energy residual
	}
	
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
	{
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
	}
	
	
	// ---------------------------------------------------------------------------------
	// output
	// ---------------------------------------------------------------------------------
	
	temporalHook();
	finalHookScalarRansTurbModel();
	finalHook();
	
	writeRestart();
	
	
	// ---------------------------------------------------------------------------------
	// delete memory
	// ---------------------------------------------------------------------------------
	
	delete [] A;
	delete [] dq;
	delete [] rhs;
	delete [] qn;
	delete [] qnm1;
	
	if (nScalUncoupled > 0) freeMem3D(AScal,   0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1);
	if (nScalUncoupled > 0) freeMem2D(dScal,   0, nScalUncoupled-1, 0, ncv_g-1);
	if (nScalUncoupled > 0) freeMem2D(rhsScal, 0, nScalUncoupled-1, 0, ncv-1);
	if (nScalUncoupled > 0) freeMem2D(qnScal, 0, nScalUncoupled-1, 0, ncv-1);
	if (nScalUncoupled > 0) freeMem2D(qnm1Scal, 0, nScalUncoupled-1, 0, ncv-1);
	
}
	
#ifdef CHECK
void JoeWithModels::runBDF2SemiCoupled()
{
  // Count the number of coupled equations.
	int nScal = scalarTranspEqVector.size();
	int nScalCoupled = 0; int nScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++) {
    if (scalarTranspEqVector[iScal].coupling == "COUPLED") nScalCoupled++;
    else nScalUncoupled++;
  }
	
	// All: NSE + scalars
  double *myResidual = new double[5+nScal];
  double *Residual = new double[5+nScal];
  double *fstResidual = new double[5+nScal];
  double *relResidual = new double[5+nScal];
		
	// NSE + coupled scalars
  double ***A;           getMem3D(&A,   0, nbocv_s-1, 0, 5+nScalCoupled-1, 0, 5+nScalCoupled-1, "JoeWithModels::runBackwardEulerCoupled -> A",   true);
  double **rhs;          getMem2D(&rhs, 0, ncv-1,     0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> rhs", true);
  double **dq;           getMem2D(&dq,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> dq",  true);
  double **qn;           getMem2D(&qn,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> qn",  true);
  double **qnm1;         getMem2D(&qnm1,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> qnm1",  true);
		
	// Uncoupled scalars
  double ***AScal   = NULL;  if (nScalUncoupled > 0) getMem3D(&AScal,    0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1, "AScal", true);
  double **dScal    = NULL;  if (nScalUncoupled > 0) getMem2D(&dScal,    0, nScalUncoupled-1, 0, ncv_g-1, "dScal", true);
  double **rhsScal  = NULL;  if (nScalUncoupled > 0) getMem2D(&rhsScal,  0, nScalUncoupled-1, 0, ncv-1, "rhsScal", true);
  double **qnScal   = NULL;  if (nScalUncoupled > 0) getMem2D(&qnScal,   0, nScalUncoupled-1, 0, ncv-1, "qnScal", true);
  double **qnm1Scal = NULL;  if (nScalUncoupled > 0) getMem2D(&qnm1Scal, 0, nScalUncoupled-1, 0, ncv-1, "qnm1Scal", true);
	
	
  //------------------------------------
  // some parameters
  //------------------------------------
  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
  int mnNewton = getIntParam("N_NEWTON_ITER", "4");
	int bdf2_TVD_limiter = getIntParam("BDF2_TVD_LIMITER", "1");
	
  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
			" to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");
	
  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }
	
  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");
	
	
  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();
	
  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();
	
  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();
	
  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();
	
	
  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
	
  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;
	
  if (initial_flowfield_output == "YES") writeData(0);
	
  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();
	
  double bdf2Alfa = 0.5, bdf2Beta = 1.0, bdf2Gamma = 0.5;
	
  // store n and n-1 and solve n+1 using second order
  for (int icv=0; icv<ncv; icv++)
  {
    qn[icv][0] = rho[icv];
    qn[icv][1] = rhou[icv][0];
    qn[icv][2] = rhou[icv][1];
    qn[icv][3] = rhou[icv][2];
    qn[icv][4] = rhoE[icv];
  }
	
	int iScalCoupled = 0; int iScalUncoupled = 0;
	for (int iScal = 0; iScal < nScal; iScal++)	{
		double *phi = scalarTranspEqVector[iScal].phi;
		if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
			for (int icv = 0; icv < ncv; icv++)
				qn[icv][5 + iScalCoupled] = rho[icv]*phi[icv];
			iScalCoupled++;
		}
		else {
			for (int icv = 0; icv < ncv; icv++)
				qnScal[iScalUncoupled][icv] = rho[icv]*phi[icv];
			iScalUncoupled++;	
		}
	}
	
	
  while (done != 1) {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);
    time += dt_min;
		
		
    // store n and n-1 and solve n+1 using second order
    for (int icv=0; icv<ncv; icv++) {
      qnm1[icv][0] = qn[icv][0]; qn[icv][0] = rho[icv];
      qnm1[icv][1] = qn[icv][1]; qn[icv][1] = rhou[icv][0];
      qnm1[icv][2] = qn[icv][2]; qn[icv][2] = rhou[icv][1];
      qnm1[icv][3] = qn[icv][3]; qn[icv][3] = rhou[icv][2];
      qnm1[icv][4] = qn[icv][4]; qn[icv][4] = rhoE[icv];
	  }
		
		int iScalCoupled = 0; int iScalUncoupled = 0;
		for (int iScal = 0; iScal < nScal; iScal++) {
			double *phi = scalarTranspEqVector[iScal].phi;
			if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
				for (int icv = 0; icv < ncv; ++icv) {
					qnm1[icv][5+iScalCoupled] = qn[icv][5+iScalCoupled]; qn[icv][5+iScalCoupled] = rho[icv]*phi[icv];
				}
				iScalCoupled++;
			}
			else {
				for (int icv=0; icv<ncv; icv++) {
					qnm1Scal[iScalUncoupled][icv] = qnScal[iScalUncoupled][icv]; qnScal[iScalUncoupled][icv] = rho[icv]*phi[icv];
				}
				iScalUncoupled++;
			}
		}

    // ---------------------------------------------------------------------------------
    // start Newton iterations
    // ---------------------------------------------------------------------------------
		
		
    if ((mpi_rank == 0) && (step%check_interval == 0))
      printf("dt: %.6le\n", dt_min);
		
    relResidual[4] = 1.0e20;
		
    int pN = 0;
    while ((pN < mnNewton) && (relResidual[4] > 1.0e-6))              // newton steps!!!
    {
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] = 0.0;
			
      for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
        for (int i = 0; i < 5+nScalCoupled; i++)
          for (int j = 0; j < 5+nScalCoupled; j++)
            A[noc][i][j] = 0.0;
			
      for (int icv = 0; icv < ncv_g; icv++)
        for (int i = 0; i < 5+nScalCoupled; i++)
          dq[icv][i] = 0.0;
			
      for (int iScal = 0; iScal < nScalUncoupled; iScal++) {                         // set AScal, dScal to zero! rhs is set zero in calcRHS
        for (int i = 0; i <= 5; i++)
          for (int noc = 0; noc < nbocv_s; noc++)
            AScal[iScal][i][noc] = 0.0;
        for (int icv = 0; icv < ncv_g; icv++)
          dScal[iScal][icv] = 0.0;
      }
			
			
      // ---------------------------------------------------------------------------------
      // calculate RHS for both NSE and scalars
      // ---------------------------------------------------------------------------------
			
			calcRhsSemiCoupled(rhs, rhsScal, A, AScal, nScalCoupled, nScalUncoupled, true);
			
      // ---------------------------------------------------------------------------------
      // solve linear system for the NSE
      // ---------------------------------------------------------------------------------

			// Add the 2nd order derivative contgribution to the residual for the Navier-Stokes Equations 
      for (int icv=0; icv<ncv; ++icv)                                 // prepare rhs and A
      {
        double tmp = cv_volume[icv]/(local_dt[icv]);
				
				double psi[5] = {1.0, 1.0, 1.0, 1.0, 1.0};

        if (bdf2_TVD_limiter == 1) {
          double tvdR[5], minmod[5];
          tvdR[0] = (rho[icv]-    qn[icv][0])/(qn[icv][0]-qnm1[icv][0]);
          tvdR[1] = (rhou[icv][0]-qn[icv][1])/(qn[icv][1]-qnm1[icv][1]);
          tvdR[2] = (rhou[icv][1]-qn[icv][2])/(qn[icv][2]-qnm1[icv][2]);
          tvdR[3] = (rhou[icv][2]-qn[icv][3])/(qn[icv][3]-qnm1[icv][3]);
          tvdR[4] = (rhoE[icv]-   qn[icv][4])/(qn[icv][4]-qnm1[icv][4]);
					
          for (int i=0; i<5; i++) {
            if (tvdR[i] > 0.0) {
              if (1.0 < fabs(tvdR[i]))    minmod[i] = 1.0;
              else                        minmod[i] = tvdR[i];
            }
            else                          minmod[i] = 0.0;
            psi[i] = pow(minmod[i]/max(1.0, fabs(tvdR[i])), 0.5);
          }
        }
				
//				rhs[icv][0] = underRelax*(rhs[icv][0] - ((1.0+psi[0]*bdf2Alfa)*rho[icv]     - (1.0+psi[0]*bdf2Beta)*qn[icv][0] + psi[0]*bdf2Gamma*qnm1[icv][0])*tmp);
 //       rhs[icv][1] = underRelax*(rhs[icv][1] - ((1.0+psi[1]*bdf2Alfa)*rhou[icv][0] - (1.0+psi[1]*bdf2Beta)*qn[icv][1] + psi[1]*bdf2Gamma*qnm1[icv][1])*tmp);
  //      rhs[icv][2] = underRelax*(rhs[icv][2] - ((1.0+psi[2]*bdf2Alfa)*rhou[icv][1] - (1.0+psi[2]*bdf2Beta)*qn[icv][2] + psi[2]*bdf2Gamma*qnm1[icv][2])*tmp);
 //       rhs[icv][3] = underRelax*(rhs[icv][3] - ((1.0+psi[3]*bdf2Alfa)*rhou[icv][2] - (1.0+psi[3]*bdf2Beta)*qn[icv][3] + psi[3]*bdf2Gamma*qnm1[icv][3])*tmp);
 //       rhs[icv][4] = underRelax*(rhs[icv][4] - ((1.0+psi[4]*bdf2Alfa)*rhoE[icv]    - (1.0+psi[4]*bdf2Beta)*qn[icv][4] + psi[4]*bdf2Gamma*qnm1[icv][4])*tmp);
				
				residField[icv] = rhs[icv][4];

        for (int i=0; i<5; i++)
          myResidual[i] += fabs(rhs[icv][i]);
				
 //       for (int i = 0; i < 5; i++)
 //         A[nbocv_i[icv]][i][i] += (1.0+psi[0]*bdf2Alfa)*tmp;             // diagonal part ( vol/dt + A )
      }
			
			// Add the 2nd order derivative contgribution to the residual for coupled scalars
			int iScalCoupled = 0;
			for (int iScal = 0; iScal < nScal; iScal++)  
				if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
					double *phi = scalarTranspEqVector[iScal].phi;
					
					for (int icv = 0; icv < ncv; ++icv) {
						double tmp = cv_volume[icv]/(local_dt[icv]);
						double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];

						double psi = 1.0;
						if (bdf2_TVD_limiter == 1) {
							double minmod, tvdR = (rhoPhi-qn[icv][5+iScalCoupled])/(qn[icv][5+iScalCoupled]-qnm1[icv][5+iScalCoupled]);
							if (tvdR > 0.0) {
								if (1.0 < fabs(tvdR))    minmod= 1.0;
								else                     minmod = tvdR;
							}
							else                       minmod = 0.0;
							psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5);
						}
						
	//					rhs[5+iScalCoupled][icv] = underRelax*( rhs[icv][5 + iScalCoupled] - ((1.0+psi*bdf2Alfa)*phi[icv] - (1.0+psi*bdf2Beta)*qn[icv][5+iScalCoupled] + psi*bdf2Gamma*qnm1[icv][5+iScalCoupled])*tmp);
						
						myResidual[5+iScalCoupled] += fabs(rhs[icv][5+iScalCoupled]);
						
	//					A[nbocv_i[icv]][5+iScalCoupled][5+iScalCoupled] += (1.0+psi*bdf2Alfa)*tmp;             // diagonal part ( vol/dt + A )
						
					}
					iScalCoupled++;
				}

			
			solveCoupledLinSysNSCoupled(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, nScalCoupled); 
			
			
      for (int icv=0; icv<ncv; icv++)                                     // update solution: Q_new = Q_old + Delta_Q
      {
        rho[icv]     += dq[icv][0];
        rhou[icv][0] += dq[icv][1];
        rhou[icv][1] += dq[icv][2];
        rhou[icv][2] += dq[icv][3];
        rhoE[icv]    += dq[icv][4];
      }
			
      updateCvData(rho,  REPLACE_DATA);
      updateCvData(rhou, REPLACE_ROTATE_DATA);
      updateCvData(rhoE, REPLACE_DATA);
			
			UpdateCvDataStateVecScal(dq, nScalCoupled);                             // update dq since neighbors needed to compute RHS of scalars
			
			iScalCoupled = 0;
			for (int iScal = 0; iScal < nScal; iScal++) {                             // update + clip solution for coupled scalars: Q_new = (rho_old * Q_old + Delta_Q ) / rho_new
				if (scalarTranspEqVector[iScal].coupling == "COUPLED") {
					double *phi = scalarTranspEqVector[iScal].phi;
					for (int icv = 0; icv < ncv; icv++)
						phi[icv] = min(max((phi[icv] * (rho[icv] - dq[icv][0]) + dq[icv][5+iScalCoupled]) / rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
					updateCvData(phi, REPLACE_DATA);
					iScalCoupled++;
				}
			}
			
      // ---------------------------------------------------------------------------------
      // solve linear system for the uncoupled scalars
      // ---------------------------------------------------------------------------------
			
      // the uncoupled scalars are solved separately from the NSE but in order to ensure consistency with
      // the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou, rhoE)
      // are used on the RHS of the equations. This means that AScal[iScal][5] is the only implicit
      // term on the LHS, while AScal[iScal][0-4] are put back to the right side.
			
			int iScalUncoupled = 0;
      for (int iScal = 0; iScal < nScal; iScal++)  
				if (scalarTranspEqVector[iScal].coupling != "COUPLED")
      {
        double *phi = scalarTranspEqVector[iScal].phi;
				
        for (int icv = 0; icv < ncv; ++icv)
        {
          double tmp = cv_volume[icv]/(local_dt[icv]);
          double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];
					
					
					double psi = 1.0;
          if (bdf2_TVD_limiter == 1)
          {
            double minmod, tvdR = (rhoPhi-qnScal[iScalUncoupled][icv])/(qnScal[iScalUncoupled][icv]-qnm1Scal[iScalUncoupled][icv]);
            if (tvdR > 0.0)
            {
              if (1.0 < fabs(tvdR))    minmod= 1.0;
              else                     minmod = tvdR;
            }
            else                       minmod = 0.0;
						
            psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5);
          }
					
          rhsScal[iScalUncoupled][icv] = underRelax*( rhsScal[iScalUncoupled][icv]
																						-(  (1.0+psi*bdf2Alfa)*rhoPhi
																							- (1.0+psi*bdf2Beta)*qnScal[iScalUncoupled][icv]
																							+ psi*bdf2Gamma*qnm1Scal[iScalUncoupled][icv])*tmp);
					
          int noc_f = nbocv_i[icv];
          int noc_l = nbocv_i[icv + 1] - 1;
					
          AScal[iScalUncoupled][5][noc_f] += (1.0+psi*bdf2Alfa)*tmp;                                 // diagonal part ( vol/dt + A )
					
          // move the other implicit terms to the RHS
          for (int noc = noc_f; noc <= noc_l; noc++)
            rhsScal[iScalUncoupled][icv] = rhsScal[iScalUncoupled][icv]
						- AScal[iScalUncoupled][0][noc] * dq[nbocv_v[noc]][0]
						- AScal[iScalUncoupled][1][noc] * dq[nbocv_v[noc]][1]
						- AScal[iScalUncoupled][2][noc] * dq[nbocv_v[noc]][2]
						- AScal[iScalUncoupled][3][noc] * dq[nbocv_v[noc]][3]
						- AScal[iScalUncoupled][4][noc] * dq[nbocv_v[noc]][4];
					
          myResidual[5+nScalCoupled+iScalUncoupled] += fabs(rhsScal[iScalUncoupled][icv]);
        }
				
        solveLinSysScalar(dScal[iScalUncoupled], AScal[iScalUncoupled][5], rhsScal[iScalUncoupled],
                          scalarTranspEqVector[iScal].phiZero,
                          scalarTranspEqVector[iScal].phiZeroRel,
                          scalarTranspEqVector[iScal].phiMaxiter,
                          scalarTranspEqVector[iScal].getName());
				
        // update scalars and clip
        for (int icv = 0; icv < ncv; icv++)
          phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScalUncoupled][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
				
        updateCvData(phi, REPLACE_DATA);
				
				iScalUncoupled++;

      }
			
      calcStateVariables();
			
      calcMaterialProperties();
			
      setBC();
			
      calcRansTurbViscMuet();
			
      MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
      if (mpi_rank == 0)
      {
        if (pN == 0)
          for (int i=0; i<5+nScal; i++)
            fstResidual[i] = Residual[i];
				
        for (int i=0; i<5+nScal; i++)
          relResidual[i] = Residual[i]/(fstResidual[i]+1.0e-10);
				
        if (step%check_interval == 0)
        {
          printf("Newton step residual: %d:\t", pN+1);
          for (int i=0; i<5+nScal; i++)
            printf("%.4le\t", relResidual[i]);
          printf("\n");
        }
      }
			
      pN++;
    }
		
    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;
			
      showResidue(Residual);
    }
		
    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);
		
    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);
		
    done = doneSolver(Residual[4]);   // pass the energy residual
  }
	
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }
	
	
  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------
	
  temporalHook();
  finalHook();
	
  writeRestart();
	
	
  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------
	
  delete [] A;
  delete [] dq;
  delete [] rhs;
  delete [] qn;
  delete [] qnm1;
	
  if (nScalUncoupled > 0) freeMem3D(AScal,   0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1);
  if (nScalUncoupled > 0) freeMem2D(dScal,   0, nScalUncoupled-1, 0, ncv_g-1);
  if (nScalUncoupled > 0) freeMem2D(rhsScal, 0, nScalUncoupled-1, 0, ncv-1);
  if (nScalUncoupled > 0) freeMem2D(qnScal, 0, nScalUncoupled-1, 0, ncv-1);
  if (nScalUncoupled > 0) freeMem2D(qnm1Scal, 0, nScalUncoupled-1, 0, ncv-1);
	
}
#endif

void JoeWithModels::calcRhs(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit)
{
  // set RHS to zero
  for (int icv = 0; icv < ncv; icv++)
  {
    rhs_rho[icv] = 0.0;
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv][i] = 0.0;
    rhs_rhoE[icv] = 0.0;
  }

  // set scalars RHS to zero
  for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    for (int icv = 0; icv < ncv; icv++)
      rhs_rhoScal[iScal][icv] = 0.0;

  // =======================================================================================
  // NAVIER-STOKES
  // =======================================================================================

  // compute Euler Flux for NS and scalars
  calcEulerFlux(rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, A, AScal, flagImplicit);

  // compute viscous Flux for NS
  if (mu_ref > 0.0)
    calcViscousFluxNS(rhs_rho, rhs_rhou, rhs_rhoE, A, flagImplicit);

  // add source terms to RHS of Navier-Stokes equations
  sourceHook(rhs_rho, rhs_rhou, rhs_rhoE, A);
  sourceHookRansTurb(rhs_rho, rhs_rhou, rhs_rhoE, A);
  sourceHookRansComb(rhs_rho, rhs_rhou, rhs_rhoE, A);


  // =======================================================================================
  // SCALARS
  // =======================================================================================
  for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
  {
    // compute viscous Flux for scalars and add source terms to RHS of scalar equations
    if (AScal == NULL)
    {
      if (mu_ref > 0.0)
        calcViscousFluxScalar_new(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal], flagImplicit);

      sourceHookScalarRansTurb_new(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
      sourceHookScalarRansComb_new(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
    }
    else
    {
      if (mu_ref > 0.0)
        calcViscousFluxScalar_new(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal], flagImplicit);

      sourceHookScalarRansTurb_new(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
      sourceHookScalarRansComb_new(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
    }
  }
}

void JoeWithModels::calcEulerFlux(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit)
{
  SymmPressBC = getStringParam("SYMM_BC","ORIGINAL_NUMERIC");

  int nScal = scalarTranspEqVector.size();

  double (*Apl)[5] = NULL;
  double (*Ami)[5] = NULL;

  if (flagImplicit)
  {
    Apl = new double[5][5];
    Ami = new double[5][5];

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        Apl[i][j] = Ami[i][j] = 0.0;
  }

  // Implicit matrix for scalars: the definition is not easy, so ask an expert (e.g., R. Pecnik)
  double (*AplScal)[6] = NULL;
  double (*AmiScal)[6] = NULL;

  if (flagImplicit)
  {
    if (nScal > 0) AplScal = new double[nScal][6];
    if (nScal > 0) AmiScal = new double[nScal][6];

    for (int iScal = 0; iScal < nScal; iScal++)
      for (int i = 0; i <= 5; i++)
        AplScal[iScal][i] = AmiScal[iScal][i] = 0.0;
  }

  double Frho, Frhou[3], FrhoE;
  double *FrhoScal       = NULL;         if (nScal > 0) FrhoScal       = new double[nScal];
  double *Scalar0        = NULL;         if (nScal > 0) Scalar0        = new double[nScal];            // cell face if second order
  double *Scalar1        = NULL;         if (nScal > 0) Scalar1        = new double[nScal];            // cell face if second order
  double *ScalCV0        = NULL;         if (nScal > 0) ScalCV0        = new double[nScal];            // cell center
  double *ScalCV1        = NULL;         if (nScal > 0) ScalCV1        = new double[nScal];            // cell center
  double *ScalConvTerm   = NULL;         if (nScal > 0) ScalConvTerm   = new double[nScal];            // 0 if convective term not considered, otherwise 1


  // If JST scheme, it is necessary to do some preprocessing
  if (SpaceIntName == "JST") calcJSTCoeff_Var(0);
	
  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
    
  for (int iScal = 0; iScal < nScal; iScal++)
    ScalConvTerm[iScal] = (double)scalarTranspEqVector[iScal].convTerm;

  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
	if (sndOrder == true) {
		
		// Second order reconstruction for mean-flow
		calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, limiterNavierS, sos, epsilonSDWLS);

		// Second order reconstruction for scalars
		if (stOrderScalar == false) {
			
			for (int icv = 0; icv < ncv; icv++)
				alpha_rho[icv] = 1.0;
			
			calcCvScalarGrad(grad_rho, rho, rho_bfa, gradreconstruction, limiterNavierS, rho, epsilonSDWLS, alpha_rho);
#ifdef temp_reconstruction
			calcCvScalarGrad(grad_temp, temp, T_bfa, gradreconstruction, limiterNavierS, temp, epsilonSDWLS);
#else
			calcCvScalarGrad(grad_p, press, p_bfa, gradreconstruction, limiterNavierS, press, epsilonSDWLS);
#endif
			
			for (int iScal = 0; iScal < nScal; iScal++) {
				
				if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") {
					// For the scalar reconstruction, Grad(rho*Phi) is required
					// Gradients of rho*Phi are saved in grad_phi temporarily
					// Gradients of rho*Phi are also limited like rho with alpha_rho
					// Boundary face values rho_bfa*Phi_fa are saved in Phi_fa temporarily
					double *rhoPhi = new double[ncv_g];
					double *phi = scalarTranspEqVector[iScal].phi;
					double *phi_bfa = scalarTranspEqVector[iScal].phi_bfa;
					double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
					
					// Compute rho*Phi
					for (int icv = 0; icv < ncv; icv++)
						rhoPhi[icv] = rho[icv] * phi[icv];
					updateCvData(rhoPhi, REPLACE_DATA);
					
					// Compute rho*Phi at boundary faces
					for (int ifa = 0; ifa < nfa_b; ifa++)
						phi_bfa[ifa] *= rho_bfa[ifa];
					
#ifdef alpha_limiter
					// Compute gradients of rho*Phi and limit based on rho with alpha_rho
					calcCvScalarGrad(grad_phi, rhoPhi, phi_bfa, gradreconstruction, NOLIMITER, rhoPhi, epsilonSDWLS);
					for (int icv = 0; icv < ncv; icv++)
						for (int i = 0; i < 3; i++)
							grad_phi[icv][i] *= alpha_rho[icv];
					updateCvData(grad_phi, REPLACE_ROTATE_DATA);
#else
					// Compute gradients of rho*Phi and limit based on rho*Phi
					calcCvScalarGrad(grad_phi, rhoPhi, phi_bfa, gradreconstruction, limiterNavierS, rhoPhi, epsilonSDWLS);
#endif
					
					// Compute back Phi at boundary faces
					for (int ifa = 0; ifa < nfa_b; ifa++)
						phi_bfa[ifa] /= rho_bfa[ifa];
					
					delete [] rhoPhi;
				}
				else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD") {
					double *phi = scalarTranspEqVector[iScal].phi;
					double *phi_bfa = scalarTranspEqVector[iScal].phi_bfa;
					double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
					calcCvScalarGrad(grad_phi, phi, phi_bfa, gradreconstruction, limiterNavierS, phi, epsilonSDWLS);
				}
				
				else {
					cerr << "### JoeWithModels::calcEulerFlux_new => Wrong reconstruction type for scalars! ###" << endl;
					throw(-1);
				}
			}
		}
	}
  // ===============================================================================================
  // Shock fix with cell flagging for AUSMDV Euler fluxes
  // ===============================================================================================
	if (Shock_Fix_Type != "NONE")
	{
	  // set flag variable ShockFixFlagCell to zero first for all cells
	  for (int icv=0; icv<ncv; ++icv)
	    ShockFixFlagCell[icv] = 0.0;                    // ideally should be integer and not float
	  updateCvData(ShockFixFlagCell, REPLACE_DATA);

	  // loop over all internal faces to locate shock location (use cell center values)
	  for (int ifa = nfa_b; ifa < nfa; ifa++)
	  {
	    int icv0 = cvofa[ifa][0];
	    int icv1 = cvofa[ifa][1];
	    double nVec[3] = {0.0, 0.0, 0.0};
	    double area = normVec3d(nVec, fa_normal[ifa]);
	    double surfv = 0.0;
	    double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
      double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
      double cL   = sqrt(gamma[icv0] * press[icv0] / rho[icv0]);
      double cR   = sqrt(gamma[icv1] * press[icv1] / rho[icv1]);

      // check if sonic point, if yes, flag the two cells on both sides of the interface
      if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
      {
        ShockFixFlagCell[icv0] = 1.0;
        ShockFixFlagCell[icv1] = 1.0;
      }
	  }
    updateCvData(ShockFixFlagCell, REPLACE_DATA);
	}

  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    double rho0 = rho[icv0];
    double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    double p0 = press[icv0];
    double T0 = temp[icv0];
    double h0 = enthalpy[icv0];
    double gam0 = gamma[icv0];
    double R0 = RoM[icv0];
    double kineCV0 = 0.0;          // cell center
    double kineFA0 = 0.0;          // cell face if second order

    double rho1 = rho[icv1];
    double u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    double p1 = press[icv1];
    double T1 = temp[icv1];
    double h1 = enthalpy[icv1];
    double gam1 = gamma[icv1];
    double R1 = RoM[icv1];
    double kineCV1 = 0.0;
    double kineFA1 = 0.0;

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

	  
	  if (sndOrder == true) {
		  
		  double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
		  vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
		  vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);
		  
		  // ----------------------------------------
		  // left side
		  // ----------------------------------------
		  rho0 += vecDotVec3d(r0, grad_rho[icv0]);
		  
#ifdef temp_reconstruction
		  T0 += vecDotVec3d(r0, grad_temp[icv0]);
		  if ((T0 <= 0.0) || (rho0 <= 0.0)) {
			  T0 = temp[icv0];
			  rho0 = rho[icv0];
			  myCountReducedOrder++;
		  }
#else
		  p0 += vecDotVec3d(r0, grad_p[icv0]);
		  if ((p0 <= 0.0) || (rho0 <= 0.0)) {
			  p0 = press[icv0];
			  rho0 = rho[icv0];
			  myCountReducedOrder++;
		  }
#endif
		  else {
			  for (int i = 0; i < 3; i++)
				  u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
			  
			  if (stOrderScalar == false) {
				  for (int iScal = 0; iScal < nScal; iScal++) {
					  double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
					  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
						  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
					  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
						  Scalar0[iScal] = Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0]);
				  }
			  }
			  
		  }
		  
		  // ----------------------------------------
		  // right side
		  // ----------------------------------------
		  rho1 += vecDotVec3d(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
		  T1 += vecDotVec3d(r1, grad_temp[icv1]);
		  if ((T1 <= 0.0) || (rho1 <= 0.0)) {
			  T1 = temp[icv1];
			  rho1 = rho[icv1];
			  myCountReducedOrder++;
		  }
#else
		  p1 += vecDotVec3d(r1, grad_p[icv1]);
		  if ((p1 <= 0.0) || (rho1 <= 0.0)) {
			  p1 = press[icv1];
			  rho1 = rho[icv1];
			  myCountReducedOrder++;
		  }
#endif
		  else {
			  for (int i = 0; i < 3; i++)
				  u1[i] += vecDotVec3d(r1, grad_u[icv1][i]);
			  
			  if (stOrderScalar == false) {
				  for (int iScal = 0; iScal < nScal; iScal++) {
					  double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
					  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
						  Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d(r1, grad_phi[icv1])) / rho1;
					  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
						  Scalar1[iScal] = Scalar1[iScal] + vecDotVec3d(r1, grad_phi[icv1]);
				  }
			  }
			  
		  }
		  
		  
		  // .............................................................................................
		  // calculation of other variables at faces: p/T, h, R, gam
		  // .............................................................................................
#ifdef temp_reconstruction
		  calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
		  calcThermoProp_T(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
		  calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
		  calcThermoProp_p(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif
		  
	  }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];         // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];         // cell face right,
    }


    // Shock Fix for AUSMDV, identify sonic points to use a more dissipative numerical scheme (Haenel & Schwane)
    string ShockFix = "NONE";
    if (Shock_Fix_Type != "NONE")
    {
      if ((ShockFixFlagCell[icv0] + ShockFixFlagCell[icv1]) > Shock_Fix_Flag_Sum) // either one cell is flagged (CELL_FLAG_TET), or both (CELL_FLAG_HEX)
      {
#ifdef SHOCKFIX_ACROSS_SHOCK
        // use dissipative scheme also across shock
#else
        // use dissipative scheme only in directions parallel to shock, not across shock
        double nVec[3] = {0.0, 0.0, 0.0};
        double area = normVec3d(nVec, fa_normal[ifa]);
        double surfv = 0.0;
        double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
        double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
        double cL   = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);
        double cR   = sqrt(gamma[icv1]*press[icv1]/rho[icv1]);

        if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
          // do nothing: ShockFix is already "NONE"
        else
#endif
          ShockFix = Shock_Fix_Type;
      }
    }

    // .............................................................................................
    // calculation of Euler Flux explicit using HLLC
    // .............................................................................................

	  if (SpaceIntName == "HLLC")
		  calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
							 rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
							 rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
							 area, nVec, nScal, 0.0);
	  else if (SpaceIntName == "AUSMDV")
    {
      if (ShockFix != "HAENEL")
        calcEulerFlux_AUSMDV(Frho, Frhou, FrhoE, FrhoScal,
               rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
               rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
               area, nVec, nScal, 0.0, AUSM_Type, ShockFix, Entropy_Fix);
      else
        calcEulerFlux_HAENEL(Frho, Frhou, FrhoE, FrhoScal,
               rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
               rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
               area, nVec, nScal, 0.0);
    }
	  else if (SpaceIntName == "ROE")
		  calcEulerFlux_Roe(Frho, Frhou, FrhoE, FrhoScal,
							rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
							rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
							area, nVec, nScal, 0.0);
	  else if (SpaceIntName == "AUSM")
		  calcEulerFlux_AUSM(Frho, Frhou, FrhoE, FrhoScal,
							 rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
							 rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
							 area, nVec, nScal, 0.0);
	  else if (SpaceIntName == "HAENEL")
      calcEulerFlux_HAENEL(Frho, Frhou, FrhoE, FrhoScal,
               rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
               rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
               area, nVec, nScal, 0.0);
	  else if (SpaceIntName == "JST") {
		  double PressSens0 = PressSensor[icv0]; double PressSens1 = PressSensor[icv1];
		  double Lambda0 = Lambda[icv0]; double Lambda1 = Lambda[icv1];
		  int Neighbor0 = NeighborCV[icv0]; double Neighbor1 = NeighborCV[icv1];
		  double *Und_Lapl0 = Conserv_Und_Lapl[icv0]; double *Und_Lapl1 = Conserv_Und_Lapl[icv1];
		  calcEulerFlux_JST(Frho, Frhou, FrhoE, FrhoScal,
							rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
							rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
							area, nVec, nScal, 0.0,
							PressSens0, Lambda0, Und_Lapl0, Neighbor0, 
							PressSens1, Lambda1, Und_Lapl1, Neighbor1);
	  }

    // icv0 is always valid...
    rhs_rho[icv0] -= Frho;
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv0][i] -= Frhou[i];
    rhs_rhoE[icv0] -= FrhoE;
    for (int iScal = 0; iScal < nScal; iScal++)
      rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

    // icv1 can be ghost... ? we are doing a loop (int ifa = nfa_b; ifa < nfa; ifa++) 
	// so the faces woth ghost cells are not included!!
    if (icv1 < ncv)
    {
      rhs_rho[icv1] += Frho;
      for (int i = 0; i < 3; i++)
        rhs_rhou[icv1][i] += Frhou[i];
      rhs_rhoE[icv1] += FrhoE;
      for (int iScal = 0; iScal < nScal; iScal++)
        rhs_rhoScal[iScal][icv1] += ScalConvTerm[iScal] * FrhoScal[iScal];
    }

    // .............................................................................................
    // calculate implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit)
    {
			if (SpaceIntName == "HLLC")
				calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
																	 rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
																	 rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
																	 area, nVec, nScal, 0.0);
			else if (SpaceIntName == "AUSMDV")
			{
        if (ShockFix != "HAENEL")
          calcEulerFluxMatrices_AUSMDV(Apl, Ami, AplScal, AmiScal,
                                  rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
                                  rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
                                  area, nVec, nScal, 0.0, AUSM_Type, ShockFix, Entropy_Fix);
        else
          calcEulerFluxMatrices_HAENEL(Apl, Ami, AplScal, AmiScal,
                                  rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
                                  rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
                                  area, nVec, nScal, 0.0);
			}
			else if (SpaceIntName == "ROE")
				calcEulerFluxMatrices_Roe(Apl, Ami, AplScal, AmiScal,
																	rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
																	rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
																	area, nVec, nScal, 0.0);
			else if (SpaceIntName == "AUSM")
				calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
																	rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
																	rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
																	area, nVec, nScal, 0.0);
			else if (SpaceIntName == "HAENEL")
          calcEulerFluxMatrices_HAENEL(Apl, Ami, AplScal, AmiScal,
                                  rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
                                  rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
                                  area, nVec, nScal, 0.0);
			else if (SpaceIntName == "JST")
				calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
									  rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
									  rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
									  area, nVec, nScal, 0.0);

      for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
      {
        A[noc00][i][j] += Apl[i][j];
        A[noc01][i][j] += Ami[i][j];
      }

      if (icv1 < ncv)  // if icv1 is internal...
      {
        for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
          A[noc11][i][j] -= Ami[i][j];
          A[noc10][i][j] -= Apl[i][j];
        }
      }

      for (int iScal = 0; iScal < nScal; iScal++)
        for (int i = 0; i <= 5; i++)
        {
          AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
          AScal[iScal][i][noc01] += ScalConvTerm[iScal] * AmiScal[iScal][i];

          if (icv1 < ncv)
          {
            AScal[iScal][i][noc11] -= ScalConvTerm[iScal] * AmiScal[iScal][i];
            AScal[iScal][i][noc10] -= ScalConvTerm[iScal] * AplScal[iScal][i];
          }
        }
    }


  }


  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY")
        {
			
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              double *phi_bfa = scalarTranspEqVector[iScal].phi_bfa;
              Scalar0[iScal] = phi_bfa[ifa];
            }

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi_bfa[ifa];

            if (SpaceIntName == "HLLC")
              calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSMDV")
              calcEulerFlux_AUSMDV(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "ROE")
              calcEulerFlux_Roe(Frho, Frhou, FrhoE, FrhoScal,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                      area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSM")
              calcEulerFlux_AUSM(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "HAENEL")
              calcEulerFlux_HAENEL(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "JST") {
              double PressSens0 = PressSensor[icv0]; double PressSens1 = PressSensor[icv0];
              double Lambda0 = Lambda[icv0]; double Lambda1 = Lambda[icv0];
              int Neighbor0 = NeighborCV[icv0]; double Neighbor1 = NeighborCV[icv0];
              double *Und_Lapl0 = Conserv_Und_Lapl[icv0]; double *Und_Lapl1 = Conserv_Und_Lapl[icv0];
              calcEulerFlux_JST(Frho, Frhou, FrhoE, FrhoScal,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                      area, nVec, nScal, 0.0,
                      PressSens0, Lambda0, Und_Lapl0, Neighbor0,
                      PressSens1, Lambda1, Und_Lapl1, Neighbor1);
            }



            if ((SymmPressBC == "ORIGINAL_NUMERIC") || (SymmPressBC == "NEW_NUMERIC")) {
              rhs_rho[icv0] -= Frho;
              for (int i = 0; i < 3; i++)
                rhs_rhou[icv0][i] -= Frhou[i];
              rhs_rhoE[icv0] -= FrhoE;
            }
            if (SymmPressBC == "ANALYTICAL") {
              rhs_rho[icv0] -= 0.0;
              for (int i = 0; i < 3; i++)
                rhs_rhou[icv0][i] -= p_bfa[ifa]*nVec[i]*area;
              rhs_rhoE[icv0] -= 0.0;
            }

            for (int iScal = 0; iScal < nScal; iScal++)
              rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

            if (flagImplicit)
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
                             rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                             rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                             area, nVec, nScal, 0.0);
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatrices_AUSMDV(Apl, Ami, AplScal, AmiScal,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "ROE")
                calcEulerFluxMatrices_Roe(Apl, Ami, AplScal, AmiScal,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            area, nVec, nScal, 0.0);
              else if (SpaceIntName == "AUSM")
                calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            area, nVec, nScal, 0.0);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatrices_HAENEL(Apl, Ami, AplScal, AmiScal,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            area, nVec, nScal, 0.0);
              else if (SpaceIntName == "JST")
                calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                            area, nVec, nScal, 0.0);


              if (SymmPressBC == "ORIGINAL_NUMERIC") {
                int noc00 = nbocv_i[icv0];
                for (int i = 0; i < 5; i++)
                  for (int j = 0; j < 5; j++)
                    A[noc00][i][j] += Apl[i][j];
              }

              if (SymmPressBC == "NEW_NUMERIC") {
                int noc00 = nbocv_i[icv0];
                for (int i = 0; i < 5; i++)
                  for (int j = 0; j < 5; j++)
                    A[noc00][i][j] += Apl[i][j]+Ami[i][j];
              }

              if (SymmPressBC == "ANALYTICAL") {
                double SimmJac[5][5];
                double a2 = gam_bfa[ifa]-1.0;
                double phi = 0.5*a2*(vel_bfa[ifa][0]*vel_bfa[ifa][0]+vel_bfa[ifa][1]*vel_bfa[ifa][1]+vel_bfa[ifa][2]*vel_bfa[ifa][2]);
                for (int i = 0; i<5; i++) {
                  SimmJac[0][i] = 0.0;
                  SimmJac[4][i] = 0.0;
                }
                for (int i = 0; i<3; i++) {
                  SimmJac[i+1][0] = phi*nVec[i]*area;
                  for (int j = 0; j<3; j++)
                    SimmJac[i+1][j+1] = -a2*vel_bfa[ifa][j]*nVec[i]*area;
                  SimmJac[i+1][4] = a2*nVec[i]*area;
                }

                int noc00 = nbocv_i[icv0]; // icv0's diagonal
                for (int i = 0; i < 5; i++)
                  for (int j = 0; j < 5; j++)
                    A[noc00][i][j] += SimmJac[i][j];


              }

              int noc00 = nbocv_i[icv0]; // icv0's diagonal
              for (int iScal = 0; iScal < nScal; iScal++)
                for (int i = 0; i <= 5; i++)
                  AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
            }
          }
        }
		  // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              double *phi_bfa = scalarTranspEqVector[iScal].phi_bfa;
              Scalar0[iScal] = phi_bfa[ifa];
            }

            if (SpaceIntName == "HLLC")
              calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSMDV")
              calcEulerFlux_AUSMDV(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "ROE")
              calcEulerFlux_Roe(Frho, Frhou, FrhoE, FrhoScal,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                      area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSM")
              calcEulerFlux_AUSM(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "HAENEL")
              calcEulerFlux_HAENEL(Frho, Frhou, FrhoE, FrhoScal,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "JST")
            {
              double PressSens0 = PressSensor[icv0]; double PressSens1 = PressSensor[icv0];
              double Lambda0 = Lambda[icv0]; double Lambda1 = Lambda[icv0];
              int Neighbor0 = NeighborCV[icv0]; double Neighbor1 = NeighborCV[icv0];
              double *Und_Lapl0 = Conserv_Und_Lapl[icv0]; double *Und_Lapl1 = Conserv_Und_Lapl[icv0];
              calcEulerFlux_JST(Frho, Frhou, FrhoE, FrhoScal,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                      area, nVec, nScal, 0.0,
                      PressSens0, Lambda0, Und_Lapl0, Neighbor0,
                      PressSens1, Lambda1, Und_Lapl1, Neighbor1);
            }


            rhs_rho[icv0] -= Frho;
            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
            for (int iScal = 0; iScal < nScal; iScal++)
              rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

            if (flagImplicit)
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
                               rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                               rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                               area, nVec, nScal, 0.0);
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatrices_AUSMDV(Apl, NULL, AplScal, AmiScal,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "ROE")
                calcEulerFluxMatrices_Roe(Apl, NULL, AplScal, AmiScal,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              area, nVec, nScal, 0.0);
              else if (SpaceIntName == "AUSM")
                calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              area, nVec, nScal, 0.0);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatrices_HAENEL(Apl, NULL, AplScal, AmiScal,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              area, nVec, nScal, 0.0);
              else if (SpaceIntName == "JST")
                calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, 0.0,
                              area, nVec, nScal, 0.0);
							
              int noc00 = nbocv_i[icv0]; // icv0's diagonal
              for (int i = 0; i < 5; i++)
              for (int j = 0; j < 5; j++)
                A[noc00][i][j] += Apl[i][j];

              for (int iScal = 0; iScal < nScal; iScal++)
                for (int i = 0; i <= 5; i++)
                  AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
            }
          }
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
        // .............................................................................................
        else
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            double rho0 = rho[icv0];
            double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double p0 = press[icv0];
            double T0 = temp[icv0];
            double h0 = enthalpy[icv0];
            double gam0 = gamma[icv0];
            double R0 = RoM[icv0];
            double kineCV0 = 0.0;           // cell center
            double kineFA0 = 0.0;           // cell face

            double kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              double *phi = scalarTranspEqVector[iScal].phi;
              double *phi_bfa = scalarTranspEqVector[iScal].phi_bfa;
              ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
              Scalar1[iScal] = phi_bfa[ifa];
            }

			  if (sndOrder == true) {
				  double r0[3] = {0.0, 0.0, 0.0};
				  vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
				  
				  // left side
				  rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
				  T0 += vecDotVec3d(r0, grad_temp[icv0]);
				  if ((T0 <= 0.0) || (rho0 <= 0.0)) {
					  T0 = temp[icv0];
					  rho0 = rho[icv0];
					  myCountReducedOrder++;
				  }
#else
				  p0 += vecDotVec3d(r0, grad_p[icv0]);
				  if ((p0 <= 0.0) || (rho0 <= 0.0)) {
					  p0 = press[icv0];
					  rho0 = rho[icv0];
					  myCountReducedOrder++;
				  }
#endif
				  else {
					  for (int i = 0; i < 3; i++)
						  u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
					  
					  if (stOrderScalar == false) {
						  for (int iScal = 0; iScal < nScal; iScal++)
						  {
							  double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
							  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
								  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
							  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
								  Scalar0[iScal] = Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0]);
						  }
					  }
				  }
				  
				  // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
				  calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
				  calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
			  }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];          // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            if (SpaceIntName == "HLLC")
              calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
                       rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSMDV")
              calcEulerFlux_AUSMDV(Frho, Frhou, FrhoE, FrhoScal,
                       rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                       area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "ROE")
              calcEulerFlux_Roe(Frho, Frhou, FrhoE, FrhoScal,
                      rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                      area, nVec, nScal, 0.0);
            else if (SpaceIntName == "AUSM")
              calcEulerFlux_AUSM(Frho, Frhou, FrhoE, FrhoScal,
                       rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "HAENEL")
              calcEulerFlux_HAENEL(Frho, Frhou, FrhoE, FrhoScal,
                       rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                       area, nVec, nScal, 0.0);
            else if (SpaceIntName == "JST") {
              double PressSens0 = PressSensor[icv0]; double PressSens1 = PressSensor[icv0];
              double Lambda0 = Lambda[icv0]; double Lambda1 = Lambda[icv0];
              int Neighbor0 = NeighborCV[icv0]; double Neighbor1 = NeighborCV[icv0];
              double *Und_Lapl0 = Conserv_Und_Lapl[icv0]; double *Und_Lapl1 = Conserv_Und_Lapl[icv0];
              calcEulerFlux_JST(Frho, Frhou, FrhoE, FrhoScal,
                      rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                      area, nVec, nScal, 0.0,
                      PressSens0, Lambda0, Und_Lapl0, Neighbor0,
                      PressSens1, Lambda1, Und_Lapl1, Neighbor1);
            }
			  
						
            rhs_rho[icv0] -= Frho;
            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
            for (int iScal = 0; iScal < nScal; iScal++)
              rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

            if (flagImplicit)
            {
							if (SpaceIntName == "HLLC")
								calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
																					 rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
																					 rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
																					 area, nVec, nScal, 0.0);
							//              calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
							//                        rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
							//                        rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
							//                        area, nVec, nScal, 0.0);
							else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatrices_AUSMDV(Apl, NULL, AplScal, AmiScal,
                                          rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
                                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                                          area, nVec, nScal, 0.0, AUSM_Type, "NONE", Entropy_Fix);
							else if (SpaceIntName == "ROE")
								calcEulerFluxMatrices_Roe(Apl, NULL, AplScal, AmiScal,
																					 rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
																					 rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
																					 area, nVec, nScal, 0.0);
							else if (SpaceIntName == "AUSM")
								calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
																					rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
																					rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
																					area, nVec, nScal, 0.0);
							else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatrices_HAENEL(Apl, NULL, AplScal, AmiScal,
                                          rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
                                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                                          area, nVec, nScal, 0.0);
							else if (SpaceIntName == "JST")
                calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
                              rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, kineCV0,
                              rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                              area, nVec, nScal, 0.0);

              int noc00 = nbocv_i[icv0]; // icv0's diagonal
              for (int i = 0; i < 5; i++)
              for (int j = 0; j < 5; j++)
                A[noc00][i][j] += Apl[i][j];

              for (int iScal = 0; iScal < nScal; iScal++)
                for (int i = 0; i <= 5; i++)
                  AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
            }
          }
        }
      }
    }
  }


  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INTEGER, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL)  delete [] Apl;
  if (Ami != NULL)  delete [] Ami;

  if (AplScal != NULL) delete [] AplScal;
  if (AmiScal != NULL) delete [] AmiScal;

  if (nScal > 0) delete [] FrhoScal;
  if (nScal > 0) delete [] Scalar0;
  if (nScal > 0) delete [] Scalar1;
  if (nScal > 0) delete [] ScalCV0;
  if (nScal > 0) delete [] ScalCV1;
  if (nScal > 0) delete [] ScalConvTerm;
}

void JoeWithModels::calcViscousFluxNS(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit)
{
  double (*A0)[5];
  double (*A1)[5];

  if (flagImplicit)
  {
    A0 = new double[5][5];
    A1 = new double[5][5];

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        A0[i][j] = A1[i][j] = 0.0;
  }
  else
    A0 = A1 = NULL;

  double Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;

  // save the index of kine if defined
  int kine_index = getScalarTransportIndex("kine");

  // ====================================================================
  //        compute gradients, with boundary values
  // ====================================================================
	
	string ViscLim = getStringParam("VISC_LIMITER","YES");
	if (ViscLim == "YES") {
		if (sndOrder != true)
			calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, limiterNavierS, sos, epsilonSDWLS);
		calcCvScalarGrad(grad_enthalpy, enthalpy, h_bfa, gradreconstruction, limiterNavierS, enthalpy, epsilonSDWLS);
	}
	else {
		calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, NOLIMITER, sos, epsilonSDWLS);
		calcCvScalarGrad(grad_enthalpy, enthalpy, h_bfa, gradreconstruction, NOLIMITER, enthalpy, epsilonSDWLS);
	}
	
  // ====================================================================
  // cycle through internal faces, assembling flux and matrix
  // ====================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);
    double sVec[3] = {0.0, 0.0, 0.0};
    vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
    double smag = normVec3d(sVec);
    double uAux_fa[3] = { 0.5*(vel[icv0][0]+ vel[icv1][0]),
                          0.5*(vel[icv0][1]+ vel[icv1][1]),
                          0.5*(vel[icv0][2]+ vel[icv1][2])};

    // kinetic energy, if defined
    double kine0 = 0.0;
    double kine1 = 0.0;
    double kine_fa = 0.0;
    if (kine_index > -1)
    {
      double *phi = scalarTranspEqVector[kine_index].phi;
      kine0 = phi[icv0];
      kine1 = phi[icv1];
      double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      kine_fa  = (w1*kine0 + w0*kine1) / (w0+w1);
    }

    // calculate viscous flux
    addViscFlux(Frhou, FrhoE, A0, A1,
              rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
              rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
              nonLinear[ifa], rij_diag_fa[ifa], rij_offdiag_fa[ifa], mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
              area, nVec, smag, sVec);

    if (flagImplicit)
    {
      for (int i=1; i<5; i++)
      for (int j=0; j<5; j++)
      {
        A[noc00][i][j] -= A0[i][j];
        A[noc01][i][j] -= A1[i][j];
      }

      if (icv1 < ncv)  // if icv1 is internal...
      {
        for (int i=1; i<5; i++)
        for (int j=0; j<5; j++)
        {
          A[noc11][i][j] += A1[i][j];
          A[noc10][i][j] += A0[i][j];
        }
      }
    }

    // icv0 is always valid...
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv0][i] -= Frhou[i];
    rhs_rhoE[icv0] -= FrhoE;

    // icv1 can be ghost...
    if (icv1 < ncv)
    {
      for (int i = 0; i < 3; i++)
        rhs_rhou[icv1][i] += Frhou[i];
      rhs_rhoE[icv1] += FrhoE;
    }
  }

  // ====================================================================
  // cycle through boundary faces, assembling flux and matrix
  // ====================================================================
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION
        // .............................................................................................
        if ((param->getString() == "SYMMETRY") || (param->getString() == "NOTHING"))
        {
          if (kine_index > -1)
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );
              double nVec[3] = {0.0, 0.0, 0.0};
              double area = normVec3d(nVec, fa_normal[ifa]);

              double *phi_bfa = scalarTranspEqVector[kine_index].phi_bfa;
              double kine_fa  = phi_bfa[ifa];

              //double tmp = (1.0 - nonLinear[ifa])*1.0/3.0*(rho[icv0] + rho_bfa[ifa])*kine_fa;
              double tmp = 1.0/3.0*(rho[icv0] + rho_bfa[ifa])*kine_fa;

              for (int i = 0; i < 3; i++)
                rhs_rhou[icv0][i] -= tmp*fa_normal[ifa][i];

              /*rhs_rhou[icv0][0] += nonLinear[ifa]*area*
                  (rij_diag_fa[ifa][0]    * nVec[0] + rij_offdiag_fa[ifa][0] * nVec[1] + rij_offdiag_fa[ifa][1] * nVec[2]);
              rhs_rhou[icv0][1] += nonLinear[ifa]*area*
                  (rij_offdiag_fa[ifa][0] * nVec[0] + rij_diag_fa[ifa][1]    * nVec[1] + rij_offdiag_fa[ifa][2] * nVec[2]);
              rhs_rhou[icv0][2] += nonLinear[ifa]*area*
                  (rij_offdiag_fa[ifa][1] * nVec[0] + rij_offdiag_fa[ifa][2] * nVec[1] + rij_diag_fa[ifa][2]    * nVec[2]);*/
            }
          }
        }
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);
            double sVec[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
            double smag = fabs(vecDotVec3d(sVec, nVec));  // project sVec to wall face normal

            double kine0 = 0.0;
            double kine1 = 0.0;
            double kine_fa = 0.0;
            if (kine_index > -1)
            {
              double *phi = scalarTranspEqVector[kine_index].phi;
              kine0 = phi[icv0];
            }

            // calculate viscous flux
            addViscFlux(Frhou, FrhoE, A0, NULL,
                      rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
                      rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], kine1,
                      nonLinear[ifa], rij_diag_fa[ifa], rij_offdiag_fa[ifa], mul_fa[ifa], 0.0, lamOcp_fa[ifa], kine_fa, vel_bfa[ifa],
                      area, nVec, smag, nVec);  /* <- use nVec here instead of sVec, to avoid inaccurate correction*/

            if (flagImplicit)
            {
              int noc00 = nbocv_i[icv0]; // icv0's diagonal

              for (int i=1; i<5; i++)
                for (int j=0; j<5; j++)
                  A[noc00][i][j] -= A0[i][j];
            }

            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
          }
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS
        // .............................................................................................
        else
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);
            double sVec[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
            double smag = normVec3d(sVec);

            double kine0 = 0.0;
            double kine1 = 0.0;
            double kine_fa = 0.0;
            if (kine_index > -1)
            {
              double *phi = scalarTranspEqVector[kine_index].phi;
              double *phi_bfa = scalarTranspEqVector[kine_index].phi_bfa;
              kine0 = phi[icv0];
              kine1 = phi_bfa[ifa];
              kine_fa = kine1;
            }

            // calculate viscous flux
            addViscFlux(Frhou, FrhoE, A0, NULL,
                      rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
                      rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], kine1,
                      nonLinear[ifa], rij_diag_fa[ifa], rij_offdiag_fa[ifa], mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, vel_bfa[ifa],
                      area, nVec, smag, sVec);

            if (flagImplicit)
            {
              int noc00 = nbocv_i[icv0]; // icv0's diagonal

              for (int i=1; i<5; i++)
                for (int j=0; j<5; j++)
                  A[noc00][i][j] -= A0[i][j];
            }

            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
          }
        }
      }
    }
  }

  if (A0  != NULL)  delete [] A0;
  if (A1  != NULL)  delete [] A1;
}

void JoeWithModels::setBC()
{
  static int first = 1;
  int bc_err = 0;

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // HOOK BOUNDARY CONDITION
        // .............................................................................................
        if (param->getString() == "HOOK")
        {
          if ((first) && (mpi_rank == 0))
            cout << "Applying HOOK                to zone: "<< zone->getName() << endl;

          boundaryHook(T_bfa, vel_bfa, p_bfa, &(*zone));

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // CBC BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC")
        {
          double u_bc[3], T_bc, p_bc;

          for (int i=0; i<3; i++)
            u_bc[i] = param->getDouble(i+2);
          T_bc = param->getDouble(5);
          p_bc = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC                 to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0]<< " "<< u_bc[1]<< " "
                 << u_bc[2] <<" T_bc: " << T_bc << " p_bc: " << p_bc << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
            {
              double velMagN = vecDotVec3d(vel[icv0], nVec);
              double mach = fabs(velMagN)/sos[icv0];

              if (mach >= 1.0)
              {
                T_bfa[ifa] = temp[icv0];
                for (int i=0; i<3; i++)
                  vel_bfa[ifa][i] = vel[icv0][i];
                p_bfa[ifa] = press[icv0];
              }
              else
              {
                T_bfa[ifa] = temp[icv0];
                for (int i=0; i<3; i++)
                  vel_bfa[ifa][i] = vel[icv0][i];
                p_bfa[ifa] = p_bc;
              }
            }
            else      // inlet
            {
              T_bfa[ifa] = T_bc;
              for (int i=0; i<3; i++)
                vel_bfa[ifa][i] = u_bc[i];
              p_bfa[ifa] = p_bc;
            }
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
       }
        // .............................................................................................
        // CBC SUBSONIC INLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_INLET")
        {
          double angleU[3], Ttot, htot, ptot;

          for (int i=0; i<3; i++)
            angleU[i] = param->getDouble(i+2);
          Ttot = param->getDouble(5);
          ptot = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC        to zone: "<< zone->getName() <<"\t angleU: "<<angleU[0]<<" "<<angleU[1]<<" "<<angleU[2]
                 << " Ttot: " << Ttot << " Ptot: "<< ptot << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double u[3] = {rhou[icv0][0]/rho[icv0], rhou[icv0][1]/rho[icv0], rhou[icv0][2]/rho[icv0]};
            double wPow2 = vecDotVec3d(u, u);               // velocity squared
            double vel = sqrt(wPow2);                       // apply angle to extrapolated velocity
            for (int i=0; i<3; i++)
              vel_bfa[ifa][i] = angleU[i]*vel;
            T_bfa[ifa] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
            p_bfa[ifa] = ptot;                              // save temporary total pressure
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));                  // total enthalpy from total temperature

          
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            double wPow2 = vecDotVec3d(vel_bfa[ifa], vel_bfa[ifa]);         // velocity squared          
            h_bfa[ifa] -= 0.5* wPow2;                                       // static enthalpy
          }
          
          ComputeBCProperties_H(&(*zone));                  // static temperature and thermo properties from static enthalpy
          
          // Assumes isentropic relations to determine static pressure (= constant cp)
          // At first approximation ok, but could be improved; should for now be considered in defining p_bc
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            p_bfa[ifa] = ptot * pow(T_bfa[ifa]/Ttot, gam_bfa[ifa]/(gam_bfa[ifa]-1.0));
        }
        // .............................................................................................
        // CBC SUBSONIC OUTLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_OUTLET")
        {
          double p_bc = param->getDouble(2);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC_OUTLET to zone: "<< zone->getName() << "\t pOut: "<< p_bc << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            // Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
            p_bfa[ifa] = p_bc;
            for (int i=0; i<3; i++)
              vel_bfa[ifa][i] = rhou[icv0][i]/rho[icv0];

            T_bfa[ifa] = temp[icv0];
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "SYMMETRY")
        {
          if ((first)&&(mpi_rank == 0))
            cout << "Applying SYMMETRY            to zone: "<< zone->getName() << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            // flip u, APPROXIMATION ---> take velocity at the cell center
            double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double un = vecDotVec3d(nVec, u0);
            for (int i = 0; i < 3; i++)
              vel_bfa[ifa][i] = u0[i] - 1.0*un*nVec[i];
            //assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

            T_bfa[ifa] = temp[icv0];
            p_bfa[ifa] = press[icv0];
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          if ((first)&&(mpi_rank == 0))
            cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            for (int i = 0; i < 3; i++)
              vel_bfa[ifa][i] = rhou[icv0][i]/rho[icv0];

            T_bfa[ifa] = temp[icv0];
            p_bfa[ifa] = press[icv0];
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          int i=0;
          double T_bc = 0.0;
          if ((i = param->findString("TEMP")) != 0)
            T_bc = param->getDouble(i+1);

          if ((first)&&(mpi_rank == 0))
          {
            if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
            else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            for (int i = 0; i < 3; i++)
              vel_bfa[ifa][i] = 0.0;

            if (T_bc > 0.0)   T_bfa[ifa] = T_bc;           // wall temperature
            else              T_bfa[ifa] = temp[icv0];     // adiabatic wall

            p_bfa[ifa] = press[icv0];                      // Assumes zero pressure gradient at the wall
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // FARFIELD BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "FARFIELD")
        {
          double u_bc[3], T_bc, p_bc, rho_bc, gamma_bc, c_bc, s_bc;

          for (int i=0; i<3; i++)
            u_bc[i] = param->getDouble(i+2);
          p_bc = param->getDouble(5);
          rho_bc = param->getDouble(6);
          gamma_bc = param->getDouble(7);

          double gm1 = gamma_bc - 1.0;
          double ovgm1 = 1.0/gm1;

          c_bc = sqrt(gamma_bc*p_bc/rho_bc);
          s_bc = pow(rho_bc,gamma_bc)/p_bc ;

          T_bc = p_bc/(R_gas*rho_bc);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying FARFIELD      to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0] << " "<< u_bc[1]
                 << " "<< u_bc[2] <<" rho_bc: " << rho_bc << " p_bc: " << p_bc << " gamma_bc: " << gamma_bc << " T_bc: "
                 << T_bc << endl;

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            //Compute normal velocity of freestream
            double qn_bc = vecDotVec3d(nVec,u_bc);

            //Subtract from normal velocity of mesh
            double vn_bc = qn_bc;// - normalvel_bfa[ifa];

            //Compute normal velocity of internal cell
            double qne = vecDotVec3d(nVec,vel[icv0]);

            //Compute speed of sound in internal cell
            double ce = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);

            double ac1, ac2;

            //Compute the Riemann invariants in the halo cell
            if (vn_bc > -c_bc)           // Outflow or subsonic inflow.
              ac1 = qne   + 2.0*ovgm1*ce;
            else                     // Supersonic inflow.
              ac1 = qn_bc + 2.0*ovgm1*c_bc;


            if(vn_bc > c_bc)             // Supersonic outflow.
              ac2 = qne   - 2.0*ovgm1*ce;
            else                     // Inflow or subsonic outflow.
              ac2 = qn_bc - 2.0*ovgm1*c_bc;


            double qnf = 0.5*(ac1 + ac2);
            double cf  = 0.25*(ac1 - ac2)*gm1;

            double velf[3], sf;

            if (vn_bc > 0)                                                       // Outflow
            {
              for (int i=0; i<3; i++)
                velf[i] = vel[icv0][i] + (qnf - qne)*nVec[i];
              sf = pow( rho[icv0], gamma[icv0])/press[icv0];
            }
            else                                                                 // Inflow
            {
              for (int i=0; i<3; i++)
                velf[i] = u_bc[i] + (qnf - qn_bc)*nVec[i];
              sf = s_bc;
            }

            //Compute density, pressure and velocity at boundary face
            rho_bfa[ifa] = pow( (sf*cf*cf/gamma[icv0]), ovgm1);
            for (int i=0; i<3; i++)
              vel_bfa[ifa][i] = velf[i];
            p_bfa[ifa] = rho_bfa[ifa]*cf*cf/gamma[icv0];
            T_bfa[ifa] = p_bfa[ifa]/(R_gas*rho_bfa[ifa]);

          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS
        // .............................................................................................
        else
        {
          if (mpi_rank == 0)
            cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
          bc_err = 1;
        }
      }
      else
      {
        if (mpi_rank == 0)
          cerr << "Error: no bc set for: "<< zone->getName() << endl;
        bc_err = 1;
      }
    }

  // update density at boundary using EOS
  for (int ifa = 0; ifa < nfa_b; ifa++)
    rho_bfa[ifa] = p_bfa[ifa] / (RoM_bfa[ifa] * T_bfa[ifa]);


  if (bc_err != 0)
    throw(-1);

  first = 0;
}

void JoeWithModels::showResidue(double *rhsResid)
{
  int nScal = scalarTranspEqVector.size();

  // residual label at every 10 output steps
  if ((mpi_rank == 0) && ((step%(check_interval*10) == 0) || (step == 1)))
  {
    printf("                rho          rhou-X       rhou-Y       rhou-Z       rhoE      ");
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12s", scalarTranspEqVector[iScal].getName());
    cout << endl;
  }

  // residual value at each output step
  if (mpi_rank == 0)
  {
    printf("RESID: %6d %12.4e %12.4e %12.4e %12.4e %12.4e", step, rhsResid[0], rhsResid[1], rhsResid[2], rhsResid[3], rhsResid[4]);
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12.4e", rhsResid[5+iScal]);
    cout << endl;
  }
}


void JoeWithModels::runBackwardEulerCoupled()
{
  int nScal = scalarTranspEqVector.size();
  for (int iScal = 0; iScal < nScal; iScal++)
    scalarTranspEqVector[iScal].coupling = "COUPLED";

  double *myResidual = new double[5+nScal];
  double *Residual   = new double[5+nScal];
  
  double ***A;           getMem3D(&A,   0, nbocv_s-1, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::runBackwardEulerCoupled -> A",   true);
  double **rhs;          getMem2D(&rhs, 0, ncv-1,     0, 5+nScal-1,               "JoeWithModels::runBackwardEulerCoupled -> rhs", true);
  double **dq;           getMem2D(&dq,  0, ncv_g-1,   0, 5+nScal-1,               "JoeWithModels::runBackwardEulerCoupled -> dq",  true);

  //------------------------------------
  // some parameters
  //------------------------------------
  if (checkParam("WRITE_HISTORY")) setHistoryFile();

  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
  double underRelaxScalars = getDoubleParam("UNDER_RELAXATION_SCALARS", "0.2");
	
  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

	string SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");
	if (SpaceIntName == "JST") { 
		PressSensor = new double[ncv_g]; Conserv_Und_Lapl = new double[ncv_g][5]; Lambda = new double[ncv_g];
		BoundaryCV = new bool[ncv_g]; NeighborCV = new int[ncv_g];
		p1_Und_Lapl = new double[ncv_g]; p2_Und_Lapl   = new double[ncv_g]; 
		calcJSTCoeff_Const(0);
	}

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
  int done = doneSolver(1.0e20);

  if (initial_flowfield_output == "YES")     writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

    // ---------------------------------------------------------------------------------
    // Compute RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
          A[noc][i][j] = 0.0;
    
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 5+nScal; i++)
        dq[icv][i] = 0.0;

    // ---------------------------------------------------------------------------------
    // calculate RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    
    calcRhsCoupled(rhs, A, nScal, true);

    
    // ---------------------------------------------------------------------------------
    // calculate residual
    // ---------------------------------------------------------------------------------

    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }
    
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] += fabs(rhs[icv][i]);
      
    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    for (int icv = 0; icv < ncv; icv++) residField[icv] = rhs[icv][4];

    for (int iScal = 0; iScal < nScal; iScal ++)
    {
      switch (iScal)
      {
      case 0:
    	for (int icv = 0; icv < ncv; icv++) residField0[icv] = rhs[icv][5+iScal];
    	break;
      case 1:
  	    for (int icv = 0; icv < ncv; icv++) residField1[icv] = rhs[icv][5+iScal];
  	    break;
      case 2:
  	    for (int icv = 0; icv < ncv; icv++) residField2[icv] = rhs[icv][5+iScal];
  	    break;
      case 3:
  	    for (int icv = 0; icv < ncv; icv++) residField3[icv] = rhs[icv][5+iScal];
  	    break;
      }
    }
    // ---------------------------------------------------------------------------------
    // solve linear system for the NSE
    // ---------------------------------------------------------------------------------
    
    for (int icv = 0; icv < ncv; icv++)              // prepare rhs and A
    {
      for (int i = 0; i < 5+nScal; i++)              // relaxation
        rhs[icv][i] *= underRelax;

      double tmp = cv_volume[icv]/(local_dt[icv]);
      for (int i = 0; i < 5+nScal; i++)
        A[nbocv_i[icv]][i][i] += tmp;                // diagonal part ( vol/dt + A )
    }      
    
//    For debugging Jacobi matrix A
//    cout << endl;
//    int np;
//    for (int icv = 0; icv < ncv; icv++)                                     // prepare rhs and A
//    {
//      if ((x_cv[icv][0] > 0.93) && (x_cv[icv][0] < 0.94))
//        if ((x_cv[icv][1] > 0.0) && (x_cv[icv][1] < 3.88e-5))
//          np = icv;
//    }
//    for (int i = 0; i < 5+nScal; i++)
//    {
//      printf("%5d\t%15.10lf\t", i, rhs[np][i]);
//      for (int j = 0; j < 5+nScal; j++)
//        printf("%15.10lf\t", A[nbocv_i[np]][i][j]);
//      cout << endl;
//    }
//    cout << endl;

    solveCoupledLinSysNSCoupled(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, nScal);  // solve linear system
             
    for (int icv = 0; icv < ncv; icv++)                                     // update solution for NSE: Q_new = Q_old + Delta_Q
    {
      rho[icv]     += dq[icv][0];
      rhou[icv][0] += dq[icv][1];
      rhou[icv][1] += dq[icv][2];
      rhou[icv][2] += dq[icv][3];
      rhoE[icv]    += dq[icv][4];
    }

    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)                         // update + clip solution for scalars: Q_new = (rho_old * Q_old + Delta_Q ) / rho_new
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = min(max((phi[icv] * (rho[icv] - dq[icv][0]) + dq[icv][5+iScal]) / rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
      updateCvData(phi, REPLACE_DATA);
    }


    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // show residual
    // =========================================================================================
    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);

      if (checkParam("WRITE_HISTORY"))
    	  writeHistoryFile(Residual);
    }

    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  temporalHook();
  finalHookScalarRansTurbModel();
  finalHook();

  writeRestart();



  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------
  delete [] myResidual;
  delete [] Residual;
  
  freeMem3D(A,   0, nbocv_s-1, 0, 5+nScal-1, 0, 5+nScal-1);         A   = NULL;
  freeMem2D(rhs, 0, ncv-1,     0, 5+nScal-1);                       rhs = NULL;
  freeMem2D(dq,  0, ncv_g-1,   0, 5+nScal-1);                       dq  = NULL;
}


void JoeWithModels::calcRhsCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
{
  // set RHS to zero
  for (int icv = 0; icv < ncv; icv++)
    for (int i = 0; i < 5+nScal; i++)
      rhs[icv][i] = 0.0;
  
  // compute Euler and viscous fluxes for NS and scalars
  calcFluxCoupled(rhs, A, nScal, flagImplicit);
  
  // add source terms to RHS of Navier-Stokes equations
  sourceHookCoupled(rhs, A, nScal, flagImplicit);
  sourceHookRansTurbCoupled(rhs, A, nScal, flagImplicit);
  sourceHookRansCombCoupled(rhs, A, nScal, flagImplicit);
}

void JoeWithModels::calcFluxCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
{
  SymmPressBC = getStringParam("SYMM_BC","ANALYTICAL");

  double **Apl = NULL;
  double **Ami = NULL;
  double **A0  = NULL;
  double **A1  = NULL;
  
  if (flagImplicit)
  {
    getMem2D(&Apl, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Apl", true);
    getMem2D(&Ami, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Ami", true);
    for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        Apl[i][j] = 0.0;
        Ami[i][j] = 0.0;
      }
    if (mu_ref > 0.0)
    {
      getMem2D(&A0,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A0 ", true);
      getMem2D(&A1,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A1 ", true);
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A0[i][j] = 0.0;
          A1[i][j] = 0.0;
        }
    }
  }

  double *EulerFlux   = new double[5+nScal];
  double *ViscousFlux = new double[5+nScal];
  
  for (int i = 0; i < 5+nScal; i++)
  {
    EulerFlux[i]   = 0.0;
    ViscousFlux[i] = 0.0;
  }
  
  double *Scalar0        = NULL;         if (nScal > 0) Scalar0       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *Scalar1        = NULL;         if (nScal > 0) Scalar1       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *ScalCV0        = NULL;         if (nScal > 0) ScalCV0       = new double[nScal];            // cell center
  double *ScalCV1        = NULL;         if (nScal > 0) ScalCV1       = new double[nScal];            // cell center
  double (*gradScal0)[3] = NULL;         if (nScal > 0) gradScal0     = new double[nScal][3];         // gradient of scalars
  double (*gradScal1)[3] = NULL;         if (nScal > 0) gradScal1     = new double[nScal][3];         // gradient of scalars
  double *dpress_dscal0  = NULL;         if (nScal > 0) dpress_dscal0 = new double[nScal];            // derivative of pressure with respect to scalars
  double *dpress_dscal1  = NULL;         if (nScal > 0) dpress_dscal1 = new double[nScal];            // derivative of pressure with respect to scalars
  double *diffScal       = NULL;         if (nScal > 0) diffScal      = new double[nScal];            // diffusivity of scalars
  double *ConvTerm       = NULL;         if (nScal > 0) ConvTerm      = new double[nScal];            // 0 if convective term not considered, otherwise 1
  double *DiffTerm       = NULL;         if (nScal > 0) DiffTerm      = new double[nScal];            // 0 if diffusive  term not considered, otherwise 1
  
  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
    
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    ConvTerm[iScal] = (double)scalarTranspEqVector[iScal].convTerm;
    DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
    dpress_dscal0[iScal] = 0.0;
    dpress_dscal1[iScal] = 0.0;
  }
  
  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  if ((sndOrder == true) || (mu_ref > 0.0))
    calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, limiterNavierS, sos, epsilonSDWLS);                      // velocity gradients
  
	if (sndOrder == true) {
		if (stOrderScalar == false) {
#ifdef alpha_limiter
			for (int icv = 0; icv < ncv; icv++) alpha_rho[icv] = 1.0;
			calcCvScalarGrad(grad_rho, rho, rho_bfa, gradreconstruction, limiterNavierS, rho, epsilonSDWLS, alpha_rho);
#else
			calcCvScalarGrad(grad_rho, rho, rho_bfa, gradreconstruction, limiterNavierS, rho, epsilonSDWLS);                    // density gradients
#endif
#ifdef temp_reconstruction
			calcCvScalarGrad(grad_temp, temp, T_bfa, gradreconstruction, limiterNavierS, temp, epsilonSDWLS);                   // temperature gradients
#else
			calcCvScalarGrad(grad_p, press, p_bfa, gradreconstruction, limiterNavierS, press, epsilonSDWLS);                    // pressure gradients
#endif
		}
	}
  
  if (mu_ref > 0.0)
    calcCvScalarGrad(grad_enthalpy, enthalpy, h_bfa, gradreconstruction, limiterNavierS, enthalpy, epsilonSDWLS);       // enthalpy gradients
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    string scalName(scalarTranspEqVector[iScal].getName());

    if ((scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") && (sndOrder == true) && (stOrderScalar == false))
    {
      // For the scalar reconstruction, Grad(rho * phi) is required
      double *phi              = scalarTranspEqVector[iScal].phi;
      double *rhoPhi           = scalarTranspEqVector[iScal].rhophi;
      double *phi_bfa          = scalarTranspEqVector[iScal].phi_bfa;
      double *rhophi_bfa       = scalarTranspEqVector[iScal].rhophi_bfa;
      double (*grad_rhophi)[3] = scalarTranspEqVector[iScal].grad_rhophi;

      // Compute rho * phi
      for (int icv = 0; icv < ncv; icv++)
        rhoPhi[icv] = rho[icv] * phi[icv];
      updateCvData(rhoPhi, REPLACE_DATA);

      // Compute rho * phi at boundary faces
      for (int ifa = 0; ifa < nfa_b; ifa++)
        rhophi_bfa[ifa] = rho_bfa[ifa] * phi_bfa[ifa];

#ifdef alpha_limiter
      // Compute gradients of rho*Phi and limit based on rho with alpha_rho
      calcCvScalarGrad(grad_rhophi, rhoPhi, rhophi_bfa, gradreconstruction, NOLIMITER, rhoPhi, epsilonSDWLS);          // scalars gradients
      for (int icv = 0; icv < ncv; icv++)
        for (int i = 0; i < 3; i++)
          grad_rhophi[icv][i] *= alpha_rho[icv];
      updateCvData(grad_rhophi, REPLACE_ROTATE_DATA);
#else
      // Compute gradients of rho*Phi and limit based on rho*Phi
      calcCvScalarGrad(grad_rhophi, rhoPhi, rhophi_bfa, gradreconstruction, limiterNavierS, rhoPhi, epsilonSDWLS);     // scalars gradients
#endif
    }

    if ( ((scalarTranspEqVector[iScal].reconstruction == "STANDARD") && (sndOrder == true) && (stOrderScalar == false)) || ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0) */) )
    {
      double *phi           = scalarTranspEqVector[iScal].phi;
      double *phi_bfa       = scalarTranspEqVector[iScal].phi_bfa;
      double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
      calcCvScalarGrad(grad_phi, phi, phi_bfa, gradreconstruction, limiterNavierS, phi, epsilonSDWLS);                 // scalars gradients
    }
    
    if ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0)*/)
    {
      diffusivityHookScalarRansTurb(scalName);
      diffusivityHookScalarRansComb(scalName);
    }
  }
  
  if (flagImplicit)
  {
    pressureDerivativeHookScalarRansTurb();      // compute pressure derivative dP/dScal
    pressureDerivativeHookScalarRansComb();
  }
  // ===============================================================================================
  // Shock fix with cell flagging for AUSMDV euler fluxes
  // ===============================================================================================
  if (Shock_Fix_Type != "NONE")
  {
    // set flag variable ShockFixFlagCell to zero first for all cells
    for (int icv=0; icv<ncv; ++icv)
      ShockFixFlagCell[icv] = 0.0;
    updateCvData(ShockFixFlagCell, REPLACE_DATA);

    // loop over all internal faces to locate shock location (use cell center values)
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      double nVec[3] = {0.0, 0.0, 0.0};
      double area = normVec3d(nVec, fa_normal[ifa]);
      double surfv = 0.0;
      double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
      double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
      double cL   = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);
      double cR   = sqrt(gamma[icv1]*press[icv1]/rho[icv1]);

      // check if sonic point, if yes, flag the two cells on both sides of the interface
      if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
      {
        ShockFixFlagCell[icv0] = 1.0;
        ShockFixFlagCell[icv1] = 1.0;
      }
    }
    updateCvData(ShockFixFlagCell, REPLACE_DATA);
  }


  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    double rho0 = rho[icv0];
    double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    double p0 = press[icv0];
    double T0 = temp[icv0];
    double h0 = enthalpy[icv0];
    double gam0 = gamma[icv0];
    double R0 = RoM[icv0];
    double kineCV0 = 0.0;          // cell center
    double kineFA0 = 0.0;          // cell face if second order

    double rho1 = rho[icv1];
    double u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    double p1 = press[icv1];
    double T1 = temp[icv1];
    double h1 = enthalpy[icv1];
    double gam1 = gamma[icv1];
    double R1 = RoM[icv1];
    double kineCV1 = 0.0;
    double kineFA1 = 0.0;
    
    double kine_fa = 0.0;          // interpolated kinetic energy at face

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

    if (sndOrder == true)
    {
      double r0[3] = {0.0, 0.0, 0.0};
      double r1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

      // ----------------------------------------
      // left side
      // ----------------------------------------
      rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
      T0 += vecDotVec3d(r0, grad_temp[icv0]);
      if ((T0 <= 0.0) || (rho0 <= 0.0))
      {
        T0 = temp[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#else
      p0 += vecDotVec3d(r0, grad_p[icv0]);
      if ((p0 <= 0.0) || (rho0 <= 0.0))
      {
        p0 = press[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
		  
		  if (stOrderScalar == false) {
			  for (int iScal = 0; iScal < nScal; iScal++)
			  {
				  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
					  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
				  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
					  Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
			  }
		  }
      }

      // ----------------------------------------
      // right side
      // ----------------------------------------
      rho1 += vecDotVec3d(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
      T1 += vecDotVec3d(r1, grad_temp[icv1]);
      if ((T1 <= 0.0) || (rho1 <= 0.0))
      {
        T1 = temp[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#else
      p1 += vecDotVec3d(r1, grad_p[icv1]);
      if ((p1 <= 0.0) || (rho1 <= 0.0))
      {
        p1 = press[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u1[i] += vecDotVec3d(r1, grad_u[icv1][i]);
		  
		  if (stOrderScalar == false) {
			  for (int iScal = 0; iScal < nScal; iScal++) {
				  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
					  Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_rhophi[icv1])) / rho1;
				  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
					  Scalar1[iScal] += vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_phi[icv1]);
			  }
		  }
      }
     
      
      // .............................................................................................
      // calculation of other variables at faces: p/T, h, R, gam
      // .............................................................................................
#ifdef temp_reconstruction
      calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
      calcThermoProp_T(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
      calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
      calcThermoProp_p(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif

    }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];                     // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];                     // cell face right,
      
      double dx0[3] = {0.0, 0.0, 0.0};
      double dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      kine_fa  = (w1 * kineCV0 + w0 * kineCV1) / (w0 + w1);  // cell face interpolated
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Shock Fix for AUSMDV, identify sonic points to use a more dissipative numerical scheme (e.g., Haenel & Schwane)
    string ShockFix = "NONE";
    if (Shock_Fix_Type != "NONE")
    {
      if ((ShockFixFlagCell[icv0] + ShockFixFlagCell[icv1]) > Shock_Fix_Flag_Sum) // either one cell is flagged (CELL_FLAG_TET), or both (CELL_FLAG_HEX)
      {
#ifdef SHOCKFIX_ACROSS_SHOCK
        // use dissipative scheme also across shock
#else
        double nVec[3] = {0.0, 0.0, 0.0};
        double area = normVec3d(nVec, fa_normal[ifa]);
        double surfv = 0.0;
        double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
        double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
        double cL   = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);
        double cR   = sqrt(gamma[icv1]*press[icv1]/rho[icv1]);

        // use dissipative scheme only in directions parallel to shock, not across shock
        if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
          // do nothing: ShockFix is already "NONE"
        else
#endif
          ShockFix = Shock_Fix_Type;
      }
    }

    // .............................................................................................
    // calculate Euler Flux explicit using HLLC
    // .............................................................................................
    if (SpaceIntName == "HLLC")
      calcEulerFluxCoupled_HLLC(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
    else if (SpaceIntName == "AUSMDV")
    {
      if (ShockFix != "HAENEL")
        calcEulerFluxCoupled_AUSMDV(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, ShockFix, Entropy_Fix);
      else
        calcEulerFluxCoupled_HAENEL(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
    }
    else if (SpaceIntName == "HAENEL")
      calcEulerFluxCoupled_HAENEL(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

    // icv0 is always valid...
    for (int i = 0; i < 5+nScal; i++)
      rhs[icv0][i] -= EulerFlux[i];
    
    // icv1 can be ghost...
    if (icv1 < ncv)
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv1][i] += EulerFlux[i];

    // .............................................................................................
    // calculate Euler implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit)
    {
      for (int iScal = 0; iScal < nScal; iScal++)
        if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
        {
          dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];
          dpress_dscal1[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv1];
        }

      if (SpaceIntName == "HLLC")
        calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
      else if (SpaceIntName == "AUSMDV")
      {
        if (ShockFix != "HAENEL")
          calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, ShockFix, Entropy_Fix);
        else
          calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
      }
      else if (SpaceIntName == "HAENEL")
        calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

      for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        A[noc00][i][j] += Apl[i][j];
        A[noc01][i][j] += Ami[i][j];
      }

      if (icv1 < ncv)  // if icv1 is internal...
        for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc11][i][j] -= Ami[i][j];
          A[noc10][i][j] -= Apl[i][j];
        }
    }

    // .............................................................................................
    // calculate viscous Flux explicit and implicit matrix
    // .............................................................................................
    if (mu_ref > 0.0)
    {
      double sVec[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(sVec);
      double alpha = vecDotVec3d(nVec, sVec);
      assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));             // alpha should now contain s dot n...
      double uAux_fa[3] = { 0.5 * (vel[icv0][0] + vel[icv1][0]),
                            0.5 * (vel[icv0][1] + vel[icv1][1]),
                            0.5 * (vel[icv0][2] + vel[icv1][2])};
      
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        for (int i = 0; i < 3; i++)
        {
          gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
          gradScal1[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv1][i];
        }
        diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
      }      
      
      calcViscousFluxCoupled(ViscousFlux, A0, A1,
                 rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                 rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], ScalCV1, gradScal1, dpress_dscal1, kineCV1,
                 mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa, diffScal, DiffTerm, 
                 area, nVec, smag, sVec, alpha, nScal);
  
      // icv0 is always valid...
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv0][i] -= ViscousFlux[i];
      
      // icv1 can be ghost...
      if (icv1 < ncv)
        for (int i = 0; i < 5+nScal; i++)
          rhs[icv1][i] += ViscousFlux[i];
      
      if (flagImplicit)
      {
        for (int i = 0; i < 5+nScal; i++) 
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc00][i][j] += A0[i][j];
          A[noc01][i][j] += A1[i][j];
        }
  
        if (icv1 < ncv)  // if icv1 is internal...
          for (int i = 0; i < 5+nScal; i++) 
          for (int j = 0; j < 5+nScal; j++)
          {
            A[noc11][i][j] -= A1[i][j];
            A[noc10][i][j] -= A0[i][j];
          }
      }
    }

  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi_bfa[ifa];

            // Euler flux
            if ((SymmPressBC == "ORIGINAL_NUMERIC") || (SymmPressBC == "NEW_NUMERIC"))
            {
              // Euler flux
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            }
            if (SymmPressBC == "ANALYTICAL")
            {
              // Euler flux
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);   // more "correct" but apparently less stable...
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
            }

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SymmPressBC == "ORIGINAL_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              }

              if (SymmPressBC == "NEW_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "Haenel")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              }

              if (SymmPressBC == "ANALYTICAL")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);          // more "correct" but apparently less stable...
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
              }

              if (SymmPressBC == "NEW_NUMERIC")
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j]+Ami[i][j];
              }
              else
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j];
              }
            }

            // Viscous flux: only 2/3 viscosity times trace of strain rate tensor and 2/3 * rho * kine
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double tmp = 2.0 / 3.0 * ((mul_fa[ifa] + mut_fa[ifa]) * (grad_u[icv0][0][0] + grad_u[icv0][1][1] + grad_u[icv0][2][2]) + rho_bfa[ifa] * kineFA);
              rhs[icv0][1] -= area * tmp * nVec[0];
              rhs[icv0][2] -= area * tmp * nVec[1];
              rhs[icv0][3] -= area * tmp * nVec[2];
              
              if (flagImplicit)
              {
                // No implicit term considered here!
              }
            }
            
          }
        }
        
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // ASSUMES THAT IF NSE HAS NEUMANN BC, ALL SCALARS HAVE NEUMANN BC!!!
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi_bfa[ifa];

            // Euler flux
            if (SpaceIntName == "HLLC")
              calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            else if (SpaceIntName == "AUSMDV")
              calcEulerFluxCoupled_AUSMDV(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "HAENEL")
              calcEulerFluxCoupled_HAENEL(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SpaceIntName == "HLLC")
                calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              if (SymmPressBC == "NEW_NUMERIC")
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j]+Ami[i][j];
              }
              else
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j];
              }
            }

            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = 0.0;
              }
              
              calcViscousFluxCoupled(ViscousFlux, A0, A1,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa], grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa], grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA, vel_bfa[ifa], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j]; // + A1[i][j];   // !!!! Apl + Ami for Neumann boundary condition: Qr = Ql -> dF/dQ = dF/dQr + dF/dQl
              }
            }            
          }
        }
        
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalName(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalName) || scalarZoneIsDirichlet(dummy, zone->getName(), scalName)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];            // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal]  = scalarTranspEqVector[iScal].phi[icv0];
              Scalar0[iScal]  = scalarTranspEqVector[iScal].phi_bfa[ifa];
            }
            
            double kineCV0 = 0.0;
            double kineFA = 0.0;
            if (kine_Index > -1)
              kineCV0 = scalarTranspEqVector[kine_Index].phi[icv0];
                        
            // Euler flux
            if ((SymmPressBC == "ORIGINAL_NUMERIC") || (SymmPressBC == "NEW_NUMERIC"))
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            }
            if (SymmPressBC == "ANALYTICAL")
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
           }

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SymmPressBC == "ORIGINAL_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              }

              if (SymmPressBC == "NEW_NUMERIC") {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              }

              if (SymmPressBC == "ANALYTICAL") {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
              }

              if (SymmPressBC == "NEW_NUMERIC")
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j]+Ami[i][j];
              }
              else
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += Apl[i][j];
              }
            }
            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], 0.0, lamOcp_fa[ifa], kineFA, vel_bfa[ifa], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
        
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), ...)
        // .............................................................................................
        else
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalname(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalname) || scalarZoneIsDirichlet(dummy, zone->getName(), scalname)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0]; // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            double rho0 = rho[icv0];
            double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double p0 = press[icv0];
            double T0 = temp[icv0];
            double h0 = enthalpy[icv0];
            double gam0 = gamma[icv0];
            double R0 = RoM[icv0];
            double kineCV0 = 0.0;           // cell center
            double kineFA0 = 0.0;           // cell face

            double kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv0];
              Scalar1[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];
            }

            if (sndOrder == true)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

              // left side
              rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
              T0 += vecDotVec3d(r0, grad_temp[icv0]);
              if ((T0 <= 0.0) || (rho0 <= 0.0))
              {
                T0 = temp[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#else
              p0 += vecDotVec3d(r0, grad_p[icv0]);
              if ((p0 <= 0.0) || (rho0 <= 0.0))
              {
                p0 = press[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#endif
              else
              {
                for (int i = 0; i < 3; i++)
                  u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
				  
				  if (stOrderScalar == false) {
					  for (int iScal = 0; iScal < nScal; iScal++)
					  {
						  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
						  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							  Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
					  }
				  }
              }

              // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
              calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
              calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
            }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];          // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            if (SpaceIntName == "HLLC")
              calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            else if (SpaceIntName == "AUSMDV")
              calcEulerFluxCoupled_AUSMDV(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "HAENEL")
              calcEulerFluxCoupled_HAENEL(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];
            
            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              // HACK:: Q_bfa might depend on Q0 so that the Jacobi matrix Apl should be adapted, not yet implemented
              if (SpaceIntName == "HLLC")
                calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                        rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                        rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                      rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                      area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                      rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                      area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[i][j];
            }

            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
//              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
//              double alpha = 1.0;
              double smag = normVec3d(sVec);
              double alpha = vecDotVec3d(nVec, sVec);
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA1,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA1, vel_bfa[ifa], diffScal, DiffTerm, 
//                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
                         area, nVec, smag, sVec, alpha, nScal);  

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
      }
    }
  }


  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INTEGER, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL) {freeMem2D(Apl, 0, 5+nScal-1, 0, 5+nScal-1);   Apl = NULL;}
  if (Ami != NULL) {freeMem2D(Ami, 0, 5+nScal-1, 0, 5+nScal-1);   Ami = NULL;}
  if (A0  != NULL) {freeMem2D(A0,  0, 5+nScal-1, 0, 5+nScal-1);   A0  = NULL;}
  if (A1  != NULL) {freeMem2D(A1,  0, 5+nScal-1, 0, 5+nScal-1);   A1  = NULL;}
  
  if (EulerFlux     != NULL) {delete [] EulerFlux;        EulerFlux     = NULL;}
  if (ViscousFlux   != NULL) {delete [] ViscousFlux;      ViscousFlux   = NULL;}
  
  if (Scalar0       != NULL) {delete [] Scalar0;          Scalar0       = NULL;}
  if (Scalar1       != NULL) {delete [] Scalar1;          Scalar1       = NULL;}
  if (ScalCV0       != NULL) {delete [] ScalCV0;          ScalCV0       = NULL;}
  if (ScalCV1       != NULL) {delete [] ScalCV1;          ScalCV1       = NULL;}
  if (gradScal0     != NULL) {delete [] gradScal0;        gradScal0     = NULL;}
  if (gradScal1     != NULL) {delete [] gradScal1;        gradScal1     = NULL;}
  if (dpress_dscal0 != NULL) {delete [] dpress_dscal0;    dpress_dscal0 = NULL;}
  if (dpress_dscal1 != NULL) {delete [] dpress_dscal1;    dpress_dscal1 = NULL;}
  if (diffScal      != NULL) {delete [] diffScal;         diffScal      = NULL;}
  if (ConvTerm      != NULL) {delete [] ConvTerm;         ConvTerm      = NULL;}
  if (DiffTerm      != NULL) {delete [] DiffTerm;         DiffTerm      = NULL;}
}

void JoeWithModels::runBackwardEulerSemiCoupled()
{
  int nScal = scalarTranspEqVector.size();
  int nScalCoupled = 0;
  int nScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      nScalCoupled++;
    else
      nScalUncoupled++;
  }

  // All: NSE + scalars
  double *myResidual = new double[5+nScal];
  double *Residual   = new double[5+nScal];

  // NSE + coupled scalars
  double ***A;           getMem3D(&A,   0, nbocv_s-1, 0, 5+nScalCoupled-1, 0, 5+nScalCoupled-1, "JoeWithModels::runBackwardEulerCoupled -> A",   true);
  double **rhs;          getMem2D(&rhs, 0, ncv-1,     0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> rhs", true);
  double **dq;           getMem2D(&dq,  0, ncv_g-1,   0, 5+nScalCoupled-1,                      "JoeWithModels::runBackwardEulerCoupled -> dq",  true);

  // Uncoupled scalars
  double ***AScal  = NULL;  if (nScalUncoupled > 0) getMem3D(&AScal,   0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1, "AScal", true);
  double **rhsScal = NULL;  if (nScalUncoupled > 0) getMem2D(&rhsScal, 0, nScalUncoupled-1, 0, ncv-1, "rhsScal", true);
  double **dScal   = NULL;  if (nScalUncoupled > 0) getMem2D(&dScal,   0, nScalUncoupled-1, 0, ncv_g-1, "dScal", true);
  
  //------------------------------------
  // some parameters
  //------------------------------------
  if (checkParam("WRITE_HISTORY")) setHistoryFile();

  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
  double underRelaxScalars = getDoubleParam("UNDER_RELAXATION_SCALARS", "0.2");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

	string SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");
	if (SpaceIntName == "JST") { 
		PressSensor = new double[ncv_g]; Conserv_Und_Lapl = new double[ncv_g][5]; Lambda = new double[ncv_g];
		BoundaryCV = new bool[ncv_g]; NeighborCV = new int[ncv_g];
		p1_Und_Lapl = new double[ncv_g]; p2_Und_Lapl   = new double[ncv_g]; 
		calcJSTCoeff_Const(0);
	}
	
  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();

	
  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
  int done = doneSolver(1.0e20);

  if (initial_flowfield_output == "YES") writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
				
    double dt_min = calcDt(cfl);

    // ---------------------------------------------------------------------------------
    // Compute RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
      for (int i = 0; i < 5+nScalCoupled; i++)
        for (int j = 0; j < 5+nScalCoupled; j++)
          A[noc][i][j] = 0.0;
    
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 5+nScalCoupled; i++)
        dq[icv][i] = 0.0;

    for (int iScal = 0; iScal < nScalUncoupled; iScal++)                // set AScal, dScal to zero! rhs is set zero in calcRHS
    {
      for (int i = 0; i <= 5; i++)
        for (int noc = 0; noc < nbocv_s; noc++)
          AScal[iScal][i][noc] = 0.0;
      for (int icv = 0; icv < ncv_g; icv++)
        dScal[iScal][icv] = 0.0;
    }

    // ---------------------------------------------------------------------------------
    // calculate RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    
    calcRhsSemiCoupled(rhs, rhsScal, A, AScal, nScalCoupled, nScalUncoupled, true);
    
    for (int icv = 0; icv < ncv; icv++) residField[icv] = rhs[icv][4];

    for (int iScal = 0; iScal < nScalCoupled; iScal ++)
    {
      switch (iScal)
      {
      case 0:
    	for (int icv = 0; icv < ncv; icv++) residField0[icv] = rhs[icv][5+iScal];
    	break;
      case 1:
  	    for (int icv = 0; icv < ncv; icv++) residField1[icv] = rhs[icv][5+iScal];
  	    break;
      case 2:
  	    for (int icv = 0; icv < ncv; icv++) residField2[icv] = rhs[icv][5+iScal];
  	    break;
      case 3:
  	    for (int icv = 0; icv < ncv; icv++) residField3[icv] = rhs[icv][5+iScal];
  	    break;
      }
    }
    // ---------------------------------------------------------------------------------
    // solve linear system for the NSE
    // ---------------------------------------------------------------------------------
    
    for (int icv = 0; icv < ncv; icv++)               // prepare rhs and A
    {
      for (int i = 0; i < 5+nScalCoupled; i++)        // relaxation
        rhs[icv][i] *= underRelax;

      double tmp = cv_volume[icv]/(local_dt[icv]);
      for (int i = 0; i < 5+nScalCoupled; i++)
        A[nbocv_i[icv]][i][i] += tmp;                 // diagonal part ( vol/dt + A )
    }      
    
    // solve linear system    
    solveCoupledLinSysNSCoupled(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, nScalCoupled); 
             
    for (int icv = 0; icv < ncv; icv++)                                     // update solution for NSE: Q_new = Q_old + Delta_Q
    {
      rho[icv]     += dq[icv][0];
      rhou[icv][0] += dq[icv][1];
      rhou[icv][1] += dq[icv][2];
      rhou[icv][2] += dq[icv][3];
      rhoE[icv]    += dq[icv][4];
    }

    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);

    UpdateCvDataStateVecScal(dq, nScalCoupled);                             // update dq since neighbors needed to compute RHS of scalars

    int iScalCoupled = 0;
    for (int iScal = 0; iScal < nScal; iScal++)                             // update + clip solution for coupled scalars: Q_new = (rho_old * Q_old + Delta_Q ) / rho_new
    {
      if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      {
        double *phi = scalarTranspEqVector[iScal].phi;
        for (int icv = 0; icv < ncv; icv++)
          phi[icv] = min(max((phi[icv] * (rho[icv] - dq[icv][0]) + dq[icv][5+iScalCoupled]) / rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
        updateCvData(phi, REPLACE_DATA);
        iScalCoupled++;
      }
    }


    // ---------------------------------------------------------------------------------
    // solve linear system for the uncoupled scalars
    // ---------------------------------------------------------------------------------

    // the uncoupled scalars are solved separately from the NSE but in order to ensure consistency with
    // the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou, rhoE)
    // are used on the RHS of the equations. This means that AScal[iScal][5] is the only implicit
    // term on the LHS, while AScal[iScal][0-4] are put back to the right side.

    int iScalUncoupled = 0;
    for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
    {
      if (scalarTranspEqVector[iScal].coupling != "COUPLED")
      {
        string scalname = scalarTranspEqVector[iScal].getName();
        for (int icv = 0; icv < ncv; ++icv)
        {
          rhsScal[iScalUncoupled][icv] *= underRelax;
                
          int noc_f = nbocv_i[icv];
          int noc_l = nbocv_i[icv + 1] - 1;
  
          double tmp = cv_volume[icv]/(local_dt[icv]);
          AScal[iScalUncoupled][5][noc_f] += tmp;                                 // diagonal part ( vol/dt + A )        
  
          // move the other implicit terms to the RHS
          for (int noc = noc_f; noc <= noc_l; noc++)
            rhsScal[iScalUncoupled][icv] = rhsScal[iScalUncoupled][icv]
                                           - AScal[iScalUncoupled][0][noc] * dq[nbocv_v[noc]][0]
                                           - AScal[iScalUncoupled][1][noc] * dq[nbocv_v[noc]][1]
                                           - AScal[iScalUncoupled][2][noc] * dq[nbocv_v[noc]][2]
                                           - AScal[iScalUncoupled][3][noc] * dq[nbocv_v[noc]][3]
                                           - AScal[iScalUncoupled][4][noc] * dq[nbocv_v[noc]][4];
        }

        switch (iScal)
        {
        case 0:
        	for (int icv = 0; icv < ncv; icv++) residField0[icv] = rhsScal[iScalUncoupled][icv];
        	break;
        case 1:
        	for (int icv = 0; icv < ncv; icv++) residField1[icv] = rhsScal[iScalUncoupled][icv];
        	break;
        case 2:
        	for (int icv = 0; icv < ncv; icv++) residField2[icv] = rhsScal[iScalUncoupled][icv];
        	break;
        case 3:
        	for (int icv = 0; icv < ncv; icv++) residField3[icv] = rhsScal[iScalUncoupled][icv];
        	break;
        }

        solveLinSysScalar(dScal[iScalUncoupled], AScal[iScalUncoupled][5], rhsScal[iScalUncoupled],
                          scalarTranspEqVector[iScal].phiZero,
                          scalarTranspEqVector[iScal].phiZeroRel,
                          scalarTranspEqVector[iScal].phiMaxiter,
                          scalarTranspEqVector[iScal].getName());
  
        // update scalars and clip
        double *phi = scalarTranspEqVector[iScal].phi;
        for (int icv = 0; icv < ncv; icv++)
          phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScalUncoupled][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
        updateCvData(phi, REPLACE_DATA);
        
        iScalUncoupled++;
      }
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();

    
    // ---------------------------------------------------------------------------------
    // calculate and show residual
    // ---------------------------------------------------------------------------------

    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }
    
    // NSE
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        myResidual[i] += fabs(rhs[icv][i] / underRelax);
    
    // Scalars
    iScalCoupled = iScalUncoupled = 0;
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      {
        for (int icv = 0; icv < ncv; icv++)
          myResidual[5+iScal] += fabs(rhs[icv][5+iScalCoupled] / underRelax);
        iScalCoupled++;
      }
      else
      {
        for (int icv = 0; icv < ncv; icv++)
          myResidual[5+iScal] += fabs(rhsScal[iScalUncoupled][icv] / underRelax);
        iScalUncoupled++;
      }
    }
    
    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);
		
      if (checkParam("WRITE_HISTORY"))
    	  writeHistoryFile(Residual);
    }


    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  temporalHook();
  finalHookScalarRansTurbModel();
  finalHook();

  writeRestart();



  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------
  delete [] myResidual;
  delete [] Residual;
  
  freeMem3D(A,   0, nbocv_s-1, 0, 5+nScalCoupled-1, 0, 5+nScalCoupled-1);         A   = NULL;
  freeMem2D(dq,  0, ncv_g-1,   0, 5+nScalCoupled-1);                              dq  = NULL;
  freeMem2D(rhs, 0, ncv-1,     0, 5+nScalCoupled-1);                              rhs = NULL;

  if (nScalUncoupled > 0) freeMem3D(AScal,   0, nScalUncoupled-1, 0, 5, 0, nbocv_s-1);         AScal   = NULL;
  if (nScalUncoupled > 0) freeMem2D(dScal,   0, nScalUncoupled-1, 0, ncv_g-1);                 dScal   = NULL;
  if (nScalUncoupled > 0) freeMem2D(rhsScal, 0, nScalUncoupled-1, 0, ncv-1);                   rhsScal = NULL;
}


void JoeWithModels::calcRhsSemiCoupled(double **rhs, double **rhs_rhoScal, double ***A, double ***AScal, int nScalCoupled, int nScalUncoupled, int flagImplicit)
{
  int nScal = nScalCoupled + nScalUncoupled;
  
  // set RHS to zero (NSE + coupled scalars)
  for (int icv = 0; icv < ncv; icv++)
    for (int i = 0; i < 5+nScalCoupled; i++)
      rhs[icv][i] = 0.0;

  // set uncoupled scalars RHS to zero
  for (int iScalUncoupled = 0; iScalUncoupled < nScalUncoupled; iScalUncoupled++)
    for (int icv = 0; icv < ncv; icv++)
      rhs_rhoScal[iScalUncoupled][icv] = 0.0;
  
  // compute Euler and viscous fluxes for NS and scalars
  calcFluxSemiCoupled(rhs, rhs_rhoScal, A, AScal, nScalCoupled, nScalUncoupled, flagImplicit);
  
  // add source terms to RHS of Navier-Stokes equations
  sourceHookCoupled(rhs, A, nScal, flagImplicit);
//  sourceHookRansTurbCoupled(rhs, A, nScal, flagImplicit);              // HACK:: needs to make sure that if turbulence scalars coupled, modified!    
  sourceHookRansCombCoupled(rhs, A, nScal, flagImplicit);                // HACK:: needs to make sure that if uncoupled, not called!          
  
  // =======================================================================================
  // SCALARS
  // =======================================================================================
  int iScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    if (scalarTranspEqVector[iScal].coupling != "COUPLED")
    {
      // compute source terms of uncoupled scalar equations
      if (AScal == NULL)
      {
        sourceHookScalarRansTurb_new(rhs_rhoScal[iScalUncoupled], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
        sourceHookScalarRansComb_new(rhs_rhoScal[iScalUncoupled], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
      }
      else
      {
        sourceHookScalarRansTurb_new(rhs_rhoScal[iScalUncoupled], AScal[iScalUncoupled][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
        sourceHookScalarRansComb_new(rhs_rhoScal[iScalUncoupled], AScal[iScalUncoupled][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
      }
      iScalUncoupled++;
    }
  }
}

void JoeWithModels::reshapeRHSSemiCoupled(double **rhs, double **rhsScal, double *Flux, int icv, double factor, int nScal)
{
  // First NSE
  for (int i = 0; i < 5; i++)
    rhs[icv][i] += factor * Flux[i];
  
  // Then scalars
  int iScalCoupled = 0;
  int iScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    if (scalarTranspEqVector[iScal].coupling == "COUPLED")
    {
      rhs[icv][5+iScalCoupled] += factor * Flux[5+iScal];
      iScalCoupled++;
    }
    else
    {
      rhsScal[iScalUncoupled][icv] += factor * Flux[5+iScal];
      iScalUncoupled++;        
    }
  }
}

void JoeWithModels::reshapeJacobianSemiCoupled(double ***A, double ***AScal, double **A0, double **A1, int noc0, int noc1, double factor, int nScal)
{
  // Upper left quadrant (only NSE)
  for (int i = 0; i < 5; i++)
  for (int j = 0; j < 5; j++)
  {
    A[noc0][i][j] += factor * A0[i][j];
    if (A1 != NULL)
      A[noc1][i][j] += factor * A1[i][j];
  }      

  // Upper right quadrant (NSE + scalars)
  for (int i = 0; i < 5; i++)
  {
    int iScalCoupled = 0;
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      {
        A[noc0][i][5+iScalCoupled] += factor * A0[i][5+iScal];
        if (A1 != NULL)
          A[noc1][i][5+iScalCoupled] += factor * A1[i][5+iScal];            
        iScalCoupled++;
      }
    }
  }
  
  // Lower quadrants (scalars)
  int iScalCoupled = 0;
  int iScalUncoupled = 0;
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    if (scalarTranspEqVector[iScal].coupling == "COUPLED")
    {
      // left quadrant (with NSE)
      for (int i = 0; i < 5; i++)
      {
        A[noc0][5+iScalCoupled][i] += factor * A0[5+iScal][i];
        if (A1 != NULL)
          A[noc1][5+iScalCoupled][i] += factor * A1[5+iScal][i]; 
      }
      
      // right quadrant (with scalars)
      int jScalCoupled = 0;
      for (int jScal = 0; jScal < nScal; jScal++)
      {
        if (scalarTranspEqVector[jScal].coupling == "COUPLED")
        {
          A[noc0][5+iScalCoupled][5+jScalCoupled] += factor * A0[5+iScal][5+jScal];
          if (A1 != NULL)
            A[noc1][5+iScalCoupled][5+jScalCoupled] += factor * A1[5+iScal][5+jScal];               
          jScalCoupled++;
        }
      }          
      iScalCoupled++;
    }
    else
    {
      for (int i = 0; i < 5; i++)
      {
        AScal[iScalUncoupled][i][noc0] += factor * A0[5+iScal][i];
        if (A1 != NULL)
          AScal[iScalUncoupled][i][noc1] += factor * A1[5+iScal][i];
      }
      AScal[iScalUncoupled][5][noc0] += factor * A0[5+iScal][5+iScal];
      if (A1 != NULL)
        AScal[iScalUncoupled][5][noc1] += factor * A1[5+iScal][5+iScal];
      iScalUncoupled++;
    } 
  }

}

void JoeWithModels::calcFluxSemiCoupled(double **rhs, double **rhs_rhoScal, double ***A, double ***AScal, int nScalCoupled, int nScalUncoupled, int flagImplicit)
{
  SymmPressBC = getStringParam("SYMM_BC","ANALYTICAL");

  int nScal = nScalCoupled + nScalUncoupled;
  
  double **Apl = NULL;
  double **Ami = NULL;
  double **A0  = NULL;
  double **A1  = NULL;
  
  if (flagImplicit)
  {
    getMem2D(&Apl, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Apl", true);
    getMem2D(&Ami, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Ami", true);
    for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        Apl[i][j] = 0.0;
        Ami[i][j] = 0.0;
      }
    if (mu_ref > 0.0)
    {
      getMem2D(&A0,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A0 ", true);
      getMem2D(&A1,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A1 ", true);
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A0[i][j] = 0.0;
          A1[i][j] = 0.0;
        }
    }
  }

  double *EulerFlux   = new double[5+nScal];
  double *ViscousFlux = new double[5+nScal];
  
  for (int i = 0; i < 5+nScal; i++)
  {
    EulerFlux[i]   = 0.0;
    ViscousFlux[i] = 0.0;
  }
  
  double *Scalar0        = NULL;         if (nScal > 0) Scalar0       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *Scalar1        = NULL;         if (nScal > 0) Scalar1       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *ScalCV0        = NULL;         if (nScal > 0) ScalCV0       = new double[nScal];            // cell center
  double *ScalCV1        = NULL;         if (nScal > 0) ScalCV1       = new double[nScal];            // cell center
  double (*gradScal0)[3] = NULL;         if (nScal > 0) gradScal0     = new double[nScal][3];         // gradient of scalars
  double (*gradScal1)[3] = NULL;         if (nScal > 0) gradScal1     = new double[nScal][3];         // gradient of scalars
  double *dpress_dscal0  = NULL;         if (nScal > 0) dpress_dscal0 = new double[nScal];            // derivative of pressure with respect to scalars
  double *dpress_dscal1  = NULL;         if (nScal > 0) dpress_dscal1 = new double[nScal];            // derivative of pressure with respect to scalars
  double *diffScal       = NULL;         if (nScal > 0) diffScal      = new double[nScal];            // diffusivity of scalars
  double *ConvTerm       = NULL;         if (nScal > 0) ConvTerm      = new double[nScal];            // 0 if convective term not considered, otherwise 1
  double *DiffTerm       = NULL;         if (nScal > 0) DiffTerm      = new double[nScal];            // 0 if diffusive  term not considered, otherwise 1
  
  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
    
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    ConvTerm[iScal] = (double)scalarTranspEqVector[iScal].convTerm;
    DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
    dpress_dscal0[iScal] = 0.0;
    dpress_dscal1[iScal] = 0.0;
  }
  
  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  if ((sndOrder == true) || (mu_ref > 0.0))
    calcCvVectorGrad(grad_u, vel, vel_bfa, gradreconstruction, limiterNavierS, sos, epsilonSDWLS);                      // velocity gradients
  
	if (sndOrder == true) {
		if (stOrderScalar == false) {
#ifdef alpha_limiter
			for (int icv = 0; icv < ncv; icv++)
				alpha_rho[icv] = 1.0;
			calcCvScalarGrad(grad_rho, rho, rho_bfa, gradreconstruction, limiterNavierS, rho, epsilonSDWLS, alpha_rho);
#else
			calcCvScalarGrad(grad_rho, rho, rho_bfa, gradreconstruction, limiterNavierS, rho, epsilonSDWLS);                    // density gradients
#endif
#ifdef temp_reconstruction
			calcCvScalarGrad(grad_temp, temp, T_bfa, gradreconstruction, limiterNavierS, temp, epsilonSDWLS);                   // temperature gradients
#else
			calcCvScalarGrad(grad_p, press, p_bfa, gradreconstruction, limiterNavierS, press, epsilonSDWLS);                    // pressure gradients
#endif
		}
	}
  
  if (mu_ref > 0.0)
    calcCvScalarGrad(grad_enthalpy, enthalpy, h_bfa, gradreconstruction, limiterNavierS, enthalpy, epsilonSDWLS);       // enthalpy gradients
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    string scalName(scalarTranspEqVector[iScal].getName());

    if ((scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") && (sndOrder == true) && (stOrderScalar == false))
    {
      // For the scalar reconstruction, Grad(rho * phi) is required
      double *phi              = scalarTranspEqVector[iScal].phi;
      double *rhoPhi           = scalarTranspEqVector[iScal].rhophi;
      double *phi_bfa          = scalarTranspEqVector[iScal].phi_bfa;
      double *rhophi_bfa       = scalarTranspEqVector[iScal].rhophi_bfa;
      double (*grad_rhophi)[3] = scalarTranspEqVector[iScal].grad_rhophi;

      // Compute rho * phi
      for (int icv = 0; icv < ncv; icv++)
        rhoPhi[icv] = rho[icv] * phi[icv];
      updateCvData(rhoPhi, REPLACE_DATA);

      // Compute rho * phi at boundary faces
      for (int ifa = 0; ifa < nfa_b; ifa++)
        rhophi_bfa[ifa] = rho_bfa[ifa] * phi_bfa[ifa];

#ifdef alpha_limiter
      // Compute gradients of rho*Phi and limit based on rho with alpha_rho
      calcCvScalarGrad(grad_rhophi, rhoPhi, rhophi_bfa, gradreconstruction, NOLIMITER, rhoPhi, epsilonSDWLS);          // scalars gradients
      for (int icv = 0; icv < ncv; icv++)
        for (int i = 0; i < 3; i++)
          grad_rhophi[icv][i] *= alpha_rho[icv];
      updateCvData(grad_rhophi, REPLACE_ROTATE_DATA);
#else
      // Compute gradients of rho*Phi and limit based on rho*Phi
      calcCvScalarGrad(grad_rhophi, rhoPhi, rhophi_bfa, gradreconstruction, limiterNavierS, rhoPhi, epsilonSDWLS);     // scalars gradients
#endif
    }

    if ( ((scalarTranspEqVector[iScal].reconstruction == "STANDARD") && (sndOrder == true) && (stOrderScalar == false)) || ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0) */) )
    {
      double *phi           = scalarTranspEqVector[iScal].phi;
      double *phi_bfa       = scalarTranspEqVector[iScal].phi_bfa;
      double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
      calcCvScalarGrad(grad_phi, phi, phi_bfa, gradreconstruction, limiterNavierS, phi, epsilonSDWLS);                 // scalars gradients
    }
    
    if ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0)*/)
    {
      diffusivityHookScalarRansTurb(scalName);
      diffusivityHookScalarRansComb(scalName);
    }
  }
  
  if (flagImplicit)
  {
    pressureDerivativeHookScalarRansTurb();      // compute pressure derivative dP/dScal
    pressureDerivativeHookScalarRansComb();
  }
  // ===============================================================================================
  // Shock fix with cell flagging for AUSMDV euler fluxes
  // ===============================================================================================
  if (Shock_Fix_Type != "NONE")
  {
    // set flag variable ShockFixFlagCell to zero first for all cells
    for (int icv=0; icv<ncv; ++icv)
      ShockFixFlagCell[icv] = 0.0;
    updateCvData(ShockFixFlagCell, REPLACE_DATA);

    // loop over all internal faces to locate shock location (use cell center values)
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      double nVec[3] = {0.0, 0.0, 0.0};
      double area = normVec3d(nVec, fa_normal[ifa]);
      double surfv = 0.0;
      double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
      double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
      double cL   = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);
      double cR   = sqrt(gamma[icv1]*press[icv1]/rho[icv1]);

      // check if sonic point, if yes, flag the two cells on both sides of the interface
      if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
      {
        ShockFixFlagCell[icv0] = 1.0;
        ShockFixFlagCell[icv1] = 1.0;
      }
    }
    updateCvData(ShockFixFlagCell, REPLACE_DATA);
  }

  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    double rho0 = rho[icv0];
    double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    double p0 = press[icv0];
    double T0 = temp[icv0];
    double h0 = enthalpy[icv0];
    double gam0 = gamma[icv0];
    double R0 = RoM[icv0];
    double kineCV0 = 0.0;          // cell center
    double kineFA0 = 0.0;          // cell face if second order

    double rho1 = rho[icv1];
    double u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    double p1 = press[icv1];
    double T1 = temp[icv1];
    double h1 = enthalpy[icv1];
    double gam1 = gamma[icv1];
    double R1 = RoM[icv1];
    double kineCV1 = 0.0;
    double kineFA1 = 0.0;
    
    double kine_fa = 0.0;          // interpolated kinetic energy at face

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

    if (sndOrder == true)
    {
      double r0[3] = {0.0, 0.0, 0.0};
      double r1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

      // ----------------------------------------
      // left side
      // ----------------------------------------
      rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
      T0 += vecDotVec3d(r0, grad_temp[icv0]);
      if ((T0 <= 0.0) || (rho0 <= 0.0))
      {
        T0 = temp[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#else
      p0 += vecDotVec3d(r0, grad_p[icv0]);
      if ((p0 <= 0.0) || (rho0 <= 0.0))
      {
        p0 = press[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
		  
		  if (stOrderScalar == false) {
			  for (int iScal = 0; iScal < nScal; iScal++)
			  {
				  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
					  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
				  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
					  Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
			  }
		  }
      }

      // ----------------------------------------
      // right side
      // ----------------------------------------
      rho1 += vecDotVec3d(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
      T1 += vecDotVec3d(r1, grad_temp[icv1]);
      if ((T1 <= 0.0) || (rho1 <= 0.0))
      {
        T1 = temp[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#else
      p1 += vecDotVec3d(r1, grad_p[icv1]);
      if ((p1 <= 0.0) || (rho1 <= 0.0))
      {
        p1 = press[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u1[i] += vecDotVec3d(r1, grad_u[icv1][i]);
		  
		  if (stOrderScalar == false) {
			  for (int iScal = 0; iScal < nScal; iScal++)
			  {
				  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
					  Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_rhophi[icv1])) / rho1;
				  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
					  Scalar1[iScal] += vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_phi[icv1]);
			  }
		  }
      }
     
      
      // .............................................................................................
      // calculation of other variables at faces: p/T, h, R, gam
      // .............................................................................................
#ifdef temp_reconstruction
      calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
      calcThermoProp_T(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
      calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
      calcThermoProp_p(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif

    }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];                     // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];                     // cell face right,
      
      double dx0[3] = {0.0, 0.0, 0.0};
      double dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      kine_fa  = (w1 * kineCV0 + w0 * kineCV1) / (w0 + w1);  // cell face interpolated
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Shock Fix for AUSMDV, identify sonic points to use a more dissipative numerical scheme (e.g., Haenel & Schwane)
    string ShockFix = "NONE";
    if (Shock_Fix_Type != "NONE")
    {
      if ((ShockFixFlagCell[icv0] + ShockFixFlagCell[icv1]) > Shock_Fix_Flag_Sum) // either one cell is flagged (CELL_FLAG_TET), or both (CELL_FLAG_HEX)
      {
#ifdef SHOCKFIX_ACROSS_SHOCK
        // use dissipative scheme also across shock
#else
        // use dissipative scheme only in directions parallel to shock, not across shock
        double nVec[3] = {0.0, 0.0, 0.0};
        double area = normVec3d(nVec, fa_normal[ifa]);
        double surfv = 0.0;
        double unL  = vecDotVec3d(vel[icv0], nVec) - surfv;
        double unR  = vecDotVec3d(vel[icv1], nVec) - surfv;
        double cL   = sqrt(gamma[icv0]*press[icv0]/rho[icv0]);
        double cR   = sqrt(gamma[icv1]*press[icv1]/rho[icv1]);

        if ((((unL-cL) > 0.0) && ((unR-cR) < 0.0)) || (((unL+cL) > 0.0) && ((unR+cR) < 0.0)))
          // do nothing: ShockFix is already "NONE"
        else
#endif
          ShockFix = Shock_Fix_Type;
      }
    }

    // .............................................................................................
    // calculate Euler Flux explicit using HLLC or AUSMDV
    // .............................................................................................
    if (SpaceIntName == "HLLC")
      calcEulerFluxCoupled_HLLC(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
    else if (SpaceIntName == "AUSMDV")
    {
      if (ShockFix != "HAENEL")
        calcEulerFluxCoupled_AUSMDV(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, ShockFix, Entropy_Fix);
      else
        calcEulerFluxCoupled_HAENEL(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
    }
    else if (SpaceIntName == "HAENEL")
      calcEulerFluxCoupled_HAENEL(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

    // icv0 is always valid...
    reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);
    // icv1 can be ghost...
    if (icv1 < ncv)
      reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv1, 1.0, nScal);

    
    // .............................................................................................
    // calculate Euler implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit)
    {
      for (int iScal = 0; iScal < nScal; iScal++)
        if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
        {
          dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];
          dpress_dscal1[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv1];
        }

      if (SpaceIntName == "HLLC")
        calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
      else if (SpaceIntName == "AUSMDV")
      {
        if (ShockFix != "HAENEL")
          calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, ShockFix, Entropy_Fix);
        else
          calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
      }
      else if (SpaceIntName == "HAENEL")
        calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
      
      // reshape different Jacobi matrix for coupled and uncoupled equations      
      reshapeJacobianSemiCoupled(A, AScal, Apl, Ami, noc00, noc01, 1.0, nScal);
      if (icv1 < ncv)  // if icv1 is internal...
        reshapeJacobianSemiCoupled(A, AScal, Apl, Ami, noc10, noc11, -1.0, nScal);
      
    }
    
    // .............................................................................................
    // calculate viscous Flux explicit and implicit matrix
    // .............................................................................................
    if (mu_ref > 0.0)
    {
      double sVec[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(sVec);
      double alpha = vecDotVec3d(nVec, sVec);
      assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));             // alpha should now contain s dot n...
      double uAux_fa[3] = { 0.5 * (vel[icv0][0] + vel[icv1][0]),
                            0.5 * (vel[icv0][1] + vel[icv1][1]),
                            0.5 * (vel[icv0][2] + vel[icv1][2])};
      
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        for (int i = 0; i < 3; i++)
        {
          gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
          gradScal1[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv1][i];
        }
        diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
      }      
      
      calcViscousFluxCoupled(ViscousFlux, A0, A1,
                 rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                 rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], ScalCV1, gradScal1, dpress_dscal1, kineCV1,
                 mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa, diffScal, DiffTerm, 
                 area, nVec, smag, sVec, alpha, nScal);
  
      // icv0 is always valid...
      reshapeRHSSemiCoupled(rhs, rhs_rhoScal, ViscousFlux, icv0, -1.0, nScal);
      // icv1 can be ghost...
      if (icv1 < ncv)
        reshapeRHSSemiCoupled(rhs, rhs_rhoScal, ViscousFlux, icv1, 1.0, nScal);

      
      if (flagImplicit)
      {   
        // reshape different Jacobi matrix for coupled and uncoupled equations      
        reshapeJacobianSemiCoupled(A, AScal, A0, A1, noc00, noc01, 1.0, nScal);
        if (icv1 < ncv)  // if icv1 is internal...
          reshapeJacobianSemiCoupled(A, AScal, A0, A1, noc10, noc11, -1.0, nScal);
      }
    }

  }

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi_bfa[ifa];
			  
            if ((SymmPressBC == "ORIGINAL_NUMERIC") || (SymmPressBC == "NEW_NUMERIC"))
            {
              // Euler flux
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            }
            if (SymmPressBC == "ANALYTICAL")
            {
              // Euler flux
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);          // more "correct" but apparently less stable...
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
           }
            reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);
            
            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SymmPressBC == "ORIGINAL_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
              }
				
              if (SymmPressBC == "NEW_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                reshapeJacobianSemiCoupled(A, AScal, Apl, Ami, noc00, 0, 1.0, nScal);
              }
				
              if (SymmPressBC == "ANALYTICAL")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
													  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);          // more "correct" but apparently less stable...
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                            area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...
                reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
              }
				
            }

            // Viscous flux: only 2/3 viscosity times trace of strain rate tensor and 2/3 * rho * kine
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double tmp = 2.0 / 3.0 * ((mul_fa[ifa] + mut_fa[ifa]) * (grad_u[icv0][0][0] + grad_u[icv0][1][1] + grad_u[icv0][2][2]) + rho_bfa[ifa] * kineFA);
              rhs[icv0][1] -= area * tmp * nVec[0];
              rhs[icv0][2] -= area * tmp * nVec[1];
              rhs[icv0][3] -= area * tmp * nVec[2];
              
              if (flagImplicit)
              {
                // No implicit term considered here!
              }
            }
            
          }
        }
        
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // ASSUMES THAT IF NSE HAS NEUMANN BC, ALL SCALARS HAVE NEUMANN BC!!!
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi_bfa[ifa];

            // Euler flux
            if (SpaceIntName == "HLLC")
              calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            else if (SpaceIntName == "AUSMDV")
              calcEulerFluxCoupled_AUSMDV(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "HAENEL")
              calcEulerFluxCoupled_HAENEL(EulerFlux,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
           reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SpaceIntName == "HLLC")
                calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
             reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
            }

            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = 0.0;
              }
              
              calcViscousFluxCoupled(ViscousFlux, A0, A1,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa], grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa], grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA, vel_bfa[ifa], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
              
              reshapeRHSSemiCoupled(rhs, rhs_rhoScal, ViscousFlux, icv0, -1.0, nScal);
              
              if (flagImplicit)
              {
                reshapeJacobianSemiCoupled(A, AScal, A0, NULL, noc00, 0, 1.0, nScal);
              }
            }            
          }
        }
        
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalName(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalName) || scalarZoneIsDirichlet(dummy, zone->getName(), scalName)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];            // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal]  = scalarTranspEqVector[iScal].phi[icv0];
              Scalar0[iScal]  = scalarTranspEqVector[iScal].phi_bfa[ifa];
            }
            
            double kineCV0 = 0.0;
            double kineFA = 0.0;
            if (kine_Index > -1)
              kineCV0 = scalarTranspEqVector[kine_Index].phi[icv0];
                    
            // Euler flux
            if ((SymmPressBC == "ORIGINAL_NUMERIC") || (SymmPressBC == "NEW_NUMERIC"))
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);
            }
            if (SymmPressBC == "ANALYTICAL")
            {
              if (SpaceIntName == "HLLC")
                calcEulerFluxCoupled_HLLC(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxCoupled_AUSMDV(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxCoupled_HAENEL(EulerFlux,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, kineFA,
                          area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
              reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);
            }

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              if (SymmPressBC == "ORIGINAL_NUMERIC")
              {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
               reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
              }

              if (SymmPressBC == "NEW_NUMERIC") {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, Ami,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                reshapeJacobianSemiCoupled(A, AScal, Apl, Ami, noc00, 0, 1.0, nScal);
              }

              if (SymmPressBC == "ANALYTICAL") {
                if (SpaceIntName == "HLLC")
                  calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
                else if (SpaceIntName == "AUSMDV")
                  calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE", AUSM_Type, "NONE", Entropy_Fix);
                else if (SpaceIntName == "HAENEL")
                  calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, dpress_dscal0, kineFA,
                                  area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
                reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
              }
				
            }
            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], 0.0, lamOcp_fa[ifa], kineFA, vel_bfa[ifa], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
              
              reshapeRHSSemiCoupled(rhs, rhs_rhoScal, ViscousFlux, icv0, -1.0, nScal);
              
              if (flagImplicit)
              {
                reshapeJacobianSemiCoupled(A, AScal, A0, NULL, noc00, 0, 1.0, nScal);
              }
              
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
        
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), ...)
        // .............................................................................................
        else
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalname(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalname) || scalarZoneIsDirichlet(dummy, zone->getName(), scalname)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0]; // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            double rho0 = rho[icv0];
            double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double p0 = press[icv0];
            double T0 = temp[icv0];
            double h0 = enthalpy[icv0];
            double gam0 = gamma[icv0];
            double R0 = RoM[icv0];
            double kineCV0 = 0.0;           // cell center
            double kineFA0 = 0.0;           // cell face

            double kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv0];
              Scalar1[iScal] = scalarTranspEqVector[iScal].phi_bfa[ifa];
            }

            if (sndOrder == true)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

              // left side
              rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
              T0 += vecDotVec3d(r0, grad_temp[icv0]);
              if ((T0 <= 0.0) || (rho0 <= 0.0))
              {
                T0 = temp[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#else
              p0 += vecDotVec3d(r0, grad_p[icv0]);
              if ((p0 <= 0.0) || (rho0 <= 0.0))
              {
                p0 = press[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#endif
              else
              {
                for (int i = 0; i < 3; i++)
                  u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
				  
				  if (stOrderScalar == false) {
					  for (int iScal = 0; iScal < nScal; iScal++)
					  {
						  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							  Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
						  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							  Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
					  }
				  }
              }

              // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
              calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
              calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
            }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];          // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            if (SpaceIntName == "HLLC")
              calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            else if (SpaceIntName == "AUSMDV")
              calcEulerFluxCoupled_AUSMDV(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
            else if (SpaceIntName == "HAENEL")
              calcEulerFluxCoupled_HAENEL(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa], T_bfa[ifa], h_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
            reshapeRHSSemiCoupled(rhs, rhs_rhoScal, EulerFlux, icv0, -1.0, nScal);
            
            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              // HACK:: Q_bfa might depend on Q0 so that the Jacobi matrix Apl should be adapted, not yet implemented
              if (SpaceIntName == "HLLC")
                calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                        rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                        rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              else if (SpaceIntName == "AUSMDV")
                calcEulerFluxMatricesCoupled_AUSMDV(Apl, NULL,
                      rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                      area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS", AUSM_Type, "NONE", Entropy_Fix);
              else if (SpaceIntName == "HAENEL")
                calcEulerFluxMatricesCoupled_HAENEL(Apl, NULL,
                      rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, dpress_dscal0, kineCV0,
                      rho_bfa[ifa], vel_bfa[ifa], p_bfa[ifa],  T_bfa[ifa], h_bfa[ifa],     RoM_bfa[ifa], gam_bfa[ifa], Scalar1, dpress_dscal0, kineFA1,
                      area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
              reshapeJacobianSemiCoupled(A, AScal, Apl, NULL, noc00, 0, 1.0, nScal);
            }

            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
//              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
//              double alpha = 1.0;
              double smag = normVec3d(sVec);
              double alpha = vecDotVec3d(nVec, sVec);
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho_bfa[ifa], vel_bfa[ifa], grad_u[icv0], h_bfa[ifa],     grad_enthalpy[icv0], T_bfa[ifa], RoM_bfa[ifa], gam_bfa[ifa], Scalar0, gradScal0, dpress_dscal0, kineFA1,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA1, vel_bfa[ifa], diffScal, DiffTerm, 
//                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
                         area, nVec, smag, sVec, alpha, nScal);  
              
              reshapeRHSSemiCoupled(rhs, rhs_rhoScal, ViscousFlux, icv0, -1.0, nScal);
              
              if (flagImplicit)
              {
                reshapeJacobianSemiCoupled(A, AScal, A0, NULL, noc00, 0, 1.0, nScal);
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
      }
    }
  }


  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INTEGER, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL) {freeMem2D(Apl, 0, 5+nScal-1, 0, 5+nScal-1);   Apl = NULL;}
  if (Ami != NULL) {freeMem2D(Ami, 0, 5+nScal-1, 0, 5+nScal-1);   Ami = NULL;}
  if (A0  != NULL) {freeMem2D(A0,  0, 5+nScal-1, 0, 5+nScal-1);   A0  = NULL;}
  if (A1  != NULL) {freeMem2D(A1,  0, 5+nScal-1, 0, 5+nScal-1);   A1  = NULL;}
  
  if (EulerFlux     != NULL) {delete [] EulerFlux;        EulerFlux     = NULL;}
  if (ViscousFlux   != NULL) {delete [] ViscousFlux;      ViscousFlux   = NULL;}
  
  if (Scalar0       != NULL) {delete [] Scalar0;          Scalar0       = NULL;}
  if (Scalar1       != NULL) {delete [] Scalar1;          Scalar1       = NULL;}
  if (ScalCV0       != NULL) {delete [] ScalCV0;          ScalCV0       = NULL;}
  if (ScalCV1       != NULL) {delete [] ScalCV1;          ScalCV1       = NULL;}
  if (gradScal0     != NULL) {delete [] gradScal0;        gradScal0     = NULL;}
  if (gradScal1     != NULL) {delete [] gradScal1;        gradScal1     = NULL;}
  if (dpress_dscal0 != NULL) {delete [] dpress_dscal0;    dpress_dscal0 = NULL;}
  if (dpress_dscal1 != NULL) {delete [] dpress_dscal1;    dpress_dscal1 = NULL;}
  if (diffScal      != NULL) {delete [] diffScal;         diffScal      = NULL;}
  if (ConvTerm      != NULL) {delete [] ConvTerm;         ConvTerm      = NULL;}
  if (DiffTerm      != NULL) {delete [] DiffTerm;         DiffTerm      = NULL;}
}

void JoeWithModels::initializeLinelet() {
	
	double alpha = 0.5, weight, nVec[3];
	int noc00, noc01, noc11, noc10, noc_f, noc_l;
	Param *p;
	bool *check_cv;
	
	Linelet = new bool[ncv];
	check_cv = new bool [ncv]; 
	for (int icv = 0; icv<ncv; icv++) 
		check_cv[icv] = true;
	
	// The source control volume for the linelets are the wall control volumes
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
		if (zone->getKind() == FA_ZONE_BOUNDARY)
			if (getParam(p, zone->getName()))
				if (p->getString() == "WALL") {
					nSources = zone->ifa_l-zone->ifa_f+1;
					linelet_cv = new vector<int>[nSources];
					int iSources = 0;
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
						int icv0 = cvofa[ifa][0];
						linelet_cv[iSources].push_back(icv0);
						check_cv[icv0] = false;
						iSources++;
					}
				}
	
	// Add new control volumes to the source linelets
	for (int iSources = 0; iSources < nSources; iSources++) {
		bool add = true;
		int index_cv = 0;
		do {
			
			// Set the source term
			int icv = linelet_cv[iSources][index_cv];
			
			// Compute max_weight
			double max_weight = 0.0;
			noc_f = nbocv_i[icv];
			noc_l = nbocv_i[icv+1]-1;
			for (int noc = noc_f+1; noc <= noc_l; noc++) {
				int ino = nbocv_v[noc];
				int ifa_new = findFace(icv, ino);
				double area = normVec3d(nVec, fa_normal[ifa_new]);
				double volume_cv1 = cv_volume[icv];
				double volume_cno = cv_volume[ino];					
				weight = 0.5*area*((1.0/volume_cv1)+(1.0/volume_cno));
				max_weight = max(max_weight, weight);
			}
			
			// Verify if any face of the control volume must be added
			noc_f = nbocv_i[icv];
			noc_l = nbocv_i[icv+1]-1;
			add = false;
			int counter = 0;
			int ino, next_cv;
			for (int noc = noc_f+1; noc <= noc_l; noc++) {
				ino = nbocv_v[noc];
				int ifa_new = findFace(icv, ino);
				double area = normVec3d(nVec, fa_normal[ifa_new]);
				double volume_cv1 = cv_volume[icv];
				double volume_cno = cv_volume[ino];					
				weight = 0.5*area*((1.0/volume_cv1)+(1.0/volume_cno));
				
				if (check_cv[ino])
					if (weight/max_weight > alpha)
						if ((index_cv == 0) || ((index_cv > 0) && (ino != linelet_cv[iSources][index_cv-1]))) {
							add = true;
							next_cv = ino;
							counter++;
						}
			}
			
			// We have arrived to a isotropic zone (there is no line any more).
			if (counter > 1) add = false;
			
			if (add) {
				linelet_cv[iSources].push_back(next_cv);
				check_cv[next_cv] = false;
				index_cv++;
			}
		} while (add);
	}
	delete [] check_cv; 
	
	for (int ifa = 0; ifa < nfa; ifa++) {
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];
		Linelet[icv0] = false;
		Linelet[icv1] = false;
	}
	
	for (int iSources = 0; iSources < nSources; iSources++)
		for (int index_cv = 0; index_cv < linelet_cv[iSources].size()-1; index_cv++){
			int icv0 = linelet_cv[iSources][index_cv];
			int icv1 = linelet_cv[iSources][index_cv+1];
			Linelet[icv0] = true;
			Linelet[icv1] = true;
		}
}

int JoeWithModels::findFace(int cv_first, int cv_second) {
	
	if (cv_first == cv_second) return -1;
	if ((cv_first >= ncv) && (cv_second >= ncv)) return -1;
	int faocv_first_f = faocv_i[cv_first];
	int faocv_first_l = faocv_i[cv_first + 1] - 1;
	for (int faocv_first = faocv_first_f; faocv_first <= faocv_first_l; faocv_first++) {
		int faocv_second_f = faocv_i[cv_second];
		int faocv_second_l = faocv_i[cv_second + 1] - 1;
		for (int faocv_second = faocv_second_f; faocv_second <= faocv_second_l; faocv_second++)
			if (faocv_v[faocv_first] == faocv_v[faocv_second]) {
				return faocv_v[faocv_first];
			}
	}
	return -1;
}


// Multigrid implementation

void JoeWithModels::runBackwardEulerMultigrid() {
	
#define MESH_0 0
#define MESH_1 1
#define MESH_2 2
#define MESH_3 3
	
	// Fine grid level definitions
	int nScal = scalarTranspEqVector.size();
	int PreSmoothing = getIntParam("PRE_SMOOTHINGS", "2");
	int AgglomSmoothing = getIntParam("AGGLOMERATION_SMOOTHINGS", "1");
	string Simplify_MG = getStringParam("SIMPLIFY_MG","YES");
	int nMesh = 0;
	nMesh = getIntParam("MULTIGRID_LEVELS", "3");
	
	double *myResidual = new double[5+nScal];
	double *Residual   = new double[5+nScal];
	
	double *RHSrho       = new double[ncv];
	double (*RHSrhou)[3] = new double[ncv][3];
	double *RHSrhoE      = new double[ncv];
	
	p1_Und_Lapl   = new double[ncv_g]; p2_Und_Lapl   = new double[ncv_g]; 
	t1_Und_Lapl   = new double[ncv_g]; t2_Und_Lapl   = new double[ncv_g]; 
	mut1_Und_Lapl   = new double[ncv_g]; mut2_Und_Lapl   = new double[ncv_g]; 	
	
	double (*A)[5][5] = new double[nbocv_s][5][5];
	double (*dq)[5]   = new double[ncv_g][5];
	double (*rhs)[5]  = new double[ncv][5];
	
	if (nMesh > 0) {
		initializeFromFineGrid();
		cvora_mgLevel1 = new int[mpi_size+1];
    	MPI_Allgather(&ncv_mgLevel1, 1, MPI_INT, cvora_mgLevel1+1, 1, MPI_INT, mpi_comm);
    	cvora_mgLevel1[0] = 0;
    	for (int rank = 0; rank<mpi_size; ++rank)
			cvora_mgLevel1[rank+1] += cvora_mgLevel1[rank];
    	assert(cvora_mgLevel1[mpi_rank+1]-cvora_mgLevel1[mpi_rank]==ncv_mgLevel1);
	}
	
	if (nMesh > 1) {
		initializeFromMultiGrid(MESH_2);
		cvora_mgLevel2 = new int[mpi_size+1];
    	MPI_Allgather(&ncv_mgLevel2, 1, MPI_INT, cvora_mgLevel2+1, 1, MPI_INT, mpi_comm);
    	cvora_mgLevel2[0] = 0;
    	for (int rank = 0; rank<mpi_size; ++rank)
			cvora_mgLevel2[rank+1] += cvora_mgLevel2[rank];
    	assert(cvora_mgLevel2[mpi_rank+1]-cvora_mgLevel2[mpi_rank]==ncv_mgLevel2);
	}
	
	if (nMesh > 2) {
		initializeFromMultiGrid(MESH_3);
		cvora_mgLevel3 = new int[mpi_size+1];
    	MPI_Allgather(&ncv_mgLevel3, 1, MPI_INT, cvora_mgLevel3+1, 1, MPI_INT, mpi_comm);
    	cvora_mgLevel3[0] = 0;
    	for (int rank = 0; rank<mpi_size; ++rank)
			cvora_mgLevel3[rank+1] += cvora_mgLevel3[rank];
    	assert(cvora_mgLevel3[mpi_rank+1]-cvora_mgLevel3[mpi_rank]==ncv_mgLevel3);
	}
	
	MPI_Barrier(mpi_comm);
	
	double *RHSrho_mgLevel1, (*RHSrhou_mgLevel1)[3], *RHSrhoE_mgLevel1;
	double *RHSrho_mgLevel2, (*RHSrhou_mgLevel2)[3], *RHSrhoE_mgLevel2;
	double *RHSrho_mgLevel3, (*RHSrhou_mgLevel3)[3], *RHSrhoE_mgLevel3;
	double (*A_mgLevel1)[5][5], (*dq_mgLevel1)[5], (*rhs_mgLevel1)[5];
	double (*A_mgLevel2)[5][5], (*dq_mgLevel2)[5], (*rhs_mgLevel2)[5];	
	double (*A_mgLevel3)[5][5], (*dq_mgLevel3)[5], (*rhs_mgLevel3)[5];
	
	if (nMesh > 0) {
		SetRestricted_Solution(MESH_0);
		RHSrho_mgLevel1 = new double[ncv_mgLevel1];
		RHSrhou_mgLevel1 = new double[ncv_mgLevel1][3];
		RHSrhoE_mgLevel1 = new double[ncv_mgLevel1];
		
		A_mgLevel1 = new double[nbocv_s_mgLevel1][5][5];
		dq_mgLevel1 = new double[ncv_mgLevel1][5];	
		rhs_mgLevel1 = new double[ncv_mgLevel1][5];
		
		if (nMesh > 1) {
			SetRestricted_Solution(MESH_1);
			RHSrho_mgLevel2       = new double[ncv_mgLevel2];
			RHSrhou_mgLevel2 = new double[ncv_mgLevel2][3];
			RHSrhoE_mgLevel2      = new double[ncv_mgLevel2];
			
			A_mgLevel2 = new double[nbocv_s_mgLevel2][5][5];
			dq_mgLevel2   = new double[ncv_mgLevel2][5];
			rhs_mgLevel2  = new double[ncv_mgLevel2][5];
			
			if (nMesh > 2) {
				SetRestricted_Solution(MESH_2);
				RHSrho_mgLevel3       = new double[ncv_mgLevel3];
				RHSrhou_mgLevel3 = new double[ncv_mgLevel3][3];
				RHSrhoE_mgLevel3      = new double[ncv_mgLevel3];
				
				A_mgLevel3 = new double[nbocv_s_mgLevel3][5][5];
				dq_mgLevel3   = new double[ncv_mgLevel3][5];	
				rhs_mgLevel3  = new double[ncv_mgLevel3][5];
			}}}
	
	MPI_Barrier(mpi_comm);
	
	// Scalars definitions
	double ***AScal  = NULL;  if (nScal > 0) getMem3D(&AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
	double **dScal   = NULL;  if (nScal > 0) getMem2D(&dScal,   0, nScal-1, 0, ncv_g-1, "dScal");
	double **rhsScal = NULL;  if (nScal > 0) getMem2D(&rhsScal, 0, nScal-1, 0, ncv-1, "rhsScal");
	
	// Some parameters that are only relevant for backward euler
	double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");
	
	if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS")) {
		ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
			" to parameter map!" << endl;
	}
	int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
	double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
	double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");
	
	if (!checkParam("CFL_RAMP")) {
		ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
		if (mpi_rank == 0)
			cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
	}
	
  	int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  	int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  	double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  	double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");
	double dt_min_mg;
	
	string SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");
	if (SpaceIntName == "JST") { 
		PressSensor = new double[ncv_g]; Conserv_Und_Lapl = new double[ncv_g][5]; Lambda = new double[ncv_g];
		BoundaryCV = new bool[ncv_g]; NeighborCV = new int[ncv_g];
		p1_Und_Lapl = new double[ncv_g]; p2_Und_Lapl   = new double[ncv_g]; 
		calcJSTCoeff_Const(0);

		if (nMesh > 0) { 
			Lambda_mgLevel1 = new double[ncv_g_mgLevel1]; NeighborCV_mgLevel1 = new int[ncv_g_mgLevel1];
			calcJSTCoeff_Const(1);
			if (nMesh > 1) { 
				Lambda_mgLevel2 = new double[ncv_g_mgLevel2]; NeighborCV_mgLevel2 = new int[ncv_g_mgLevel2];
				calcJSTCoeff_Const(2);
				if (nMesh > 2) { 
					Lambda_mgLevel3 = new double[ncv_g_mgLevel3]; NeighborCV_mgLevel3 = new int[ncv_g_mgLevel3];
					calcJSTCoeff_Const(3);
				}
			}
		}
		
	
	}
	
	// Update state properties: velocity, pressure, temperature, enthalpy, gamma and R
	// Update material properties: laminar viscosity and heat conductivity
	// Set BC's for NS and scalars.
	// Update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars.
  	calcStateVariables(); calcMaterialProperties(); setBC(); calcRansTurbViscMuet();
	
	if (nMesh > 0) { 
		calcStateVariables_mg(MESH_1);	calcMaterialProperties_mg(MESH_1);	
		setBC_mg(MESH_1); calcRansTurbViscMuet_mg(MESH_1);
		if (nMesh > 1) { 
			calcStateVariables_mg(MESH_2);	calcMaterialProperties_mg(MESH_2);	
			setBC_mg(MESH_2); calcRansTurbViscMuet_mg(MESH_2);
			if (nMesh > 2) { 
				calcStateVariables_mg(MESH_3);	calcMaterialProperties_mg(MESH_3);	
				setBC_mg(MESH_3); calcRansTurbViscMuet_mg(MESH_3);
			}}}
	
	//   Loop over time steps
	int done = doneSolver(1.0e20);
	
	if (initial_flowfield_output == "YES") writeData(0);
	
	// provide total runtime
	double wtime, wtime0;
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();
	
	while (done != 1) {
		
		if ( step >= getIntParam("SWITCH_ON_MG", "50") ) {
			if ( step < getIntParam("SWITCH_OFF_MG", "1500") )
				nMesh = getIntParam("MULTIGRID_LEVELS", "3");
			else
				nMesh = 0;
		}
		
		// Set to zero the trucation error
		if (nMesh > 0){
			for (int icv = 0; icv < ncv_mgLevel1; icv++) {
				rho_TruncError_mgLevel1[icv] = 0.0; rhou_TruncError_mgLevel1[icv][0] = 0.0;
				rhou_TruncError_mgLevel1[icv][1] = 0.0; rhou_TruncError_mgLevel1[icv][2] = 0.0;
				rhoE_TruncError_mgLevel1[icv] = 0.0; }
			if (nMesh > 1){
				for (int icv = 0; icv < ncv_mgLevel2; icv++) {
					rho_TruncError_mgLevel2[icv] = 0.0; rhou_TruncError_mgLevel2[icv][0] = 0.0;
					rhou_TruncError_mgLevel2[icv][1] = 0.0; rhou_TruncError_mgLevel2[icv][2] = 0.0;
					rhoE_TruncError_mgLevel2[icv] = 0.0; }
				if (nMesh > 2){
					for (int icv = 0; icv < ncv_mgLevel3; icv++) {
						rho_TruncError_mgLevel3[icv] = 0.0; rhou_TruncError_mgLevel3[icv][0] = 0.0;
						rhou_TruncError_mgLevel3[icv][1] = 0.0; rhou_TruncError_mgLevel3[icv][2] = 0.0;
						rhoE_TruncError_mgLevel3[icv] = 0.0; }
				}}}
		
    	step++;
		
		// Compute time step
		if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL)) cfl *= incCFL;
		
		// Space integration and time integration on the finest level
		calcStateVariables(); calcMaterialProperties(); setBC();  calcRansTurbViscMuet();
		double dt_min = calcDt(cfl);
		
		
		// Compute RHS for both NSE and scalars
		for (int noc = 0; noc < nbocv_s; noc++)	// set A, dq to zero! rhs is set zero in calcRHS
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					A[noc][i][j] = 0.0;
		
		for (int icv = 0; icv < ncv_g; icv++)
			for (int i = 0; i < 5; i++)
				dq[icv][i] = 0.0;
		
		for (int iScal = 0; iScal < nScal; iScal++)	// set AScal, dScal to zero! rhs is set zero in calcRHS
		{
			for (int i = 0; i <= 5; i++)
				for (int noc = 0; noc < nbocv_s; noc++)
					AScal[iScal][i][noc] = 0.0;
			for (int icv = 0; icv < ncv_g; icv++)
				dScal[iScal][icv] = 0.0;
		}
		
		
		
		calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, true);
		
		calcTimeInt(RHSrho, RHSrhou, RHSrhoE, A, dq, rhs, underRelax, maxIterLS, zeroAbsLS, zeroRelLS);
		
		// Calculate and show residual on the finest grid
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = 0.0;
			Residual[i] = 0.0;
		}
		
		for (int icv = 0; icv < ncv; icv++) {
			myResidual[0] += fabs(RHSrho[icv]);
			for (int i=0; i<3; i++)
				myResidual[i+1] += fabs(RHSrhou[icv][i]);
			myResidual[4] += fabs(RHSrhoE[icv]);
		}
		
		// solve linear system for the scalars
		// the scalars are solved separately from the NSE but in order to ensure consistency with
		// the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou)
		// are used on the RHS of the equations. This means that AScal[iScal][4] is the only implicit
		// term on the LHS, while AScal[iScal][0-3] are put back to the right side.
		
		for (int iScal = 0; iScal < nScal; iScal++) {
			string scalname = scalarTranspEqVector[iScal].getName();
			for (int icv = 0; icv < ncv; ++icv) {
				rhsScal[iScal][icv] *= underRelax;
				
				int noc_f = nbocv_i[icv];
				int noc_l = nbocv_i[icv + 1] - 1;
				
				double tmp = cv_volume[icv]/(local_dt[icv]);
				AScal[iScal][5][noc_f] += tmp;	// diagonal part ( vol/dt + A )        
				
				// move the other implicit terms to the RHS
				for (int noc = noc_f; noc <= noc_l; noc++)
					rhsScal[iScal][icv] = rhsScal[iScal][icv]
					- AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
					- AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
					- AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
					- AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
					- AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];
			}
			
			solveLinSysScalar(dScal[iScal], AScal[iScal][5], rhsScal[iScal],
							  scalarTranspEqVector[iScal].phiZero,
							  scalarTranspEqVector[iScal].phiZeroRel,
							  scalarTranspEqVector[iScal].phiMaxiter,
							  scalarTranspEqVector[iScal].getName());
			
			// update scalars and clip
			double *phi = scalarTranspEqVector[iScal].phi;
			for (int icv = 0; icv < ncv; icv++)
				phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScal][icv])/rho[icv], 
								   scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
			updateCvData(phi, REPLACE_DATA);
			
		}
		
		for (int iScal = 0; iScal < nScal; iScal++)
			for (int icv = 0; icv < ncv; icv++)
				myResidual[5+iScal] += fabs(rhsScal[iScal][icv]/underRelax);
		
		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
		calcStateVariables(); calcMaterialProperties(); setBC(); calcRansTurbViscMuet();
		
		if (nMesh >= 1) {
			// LEVEL 1 - Compute $r^*_(k+1) = Damp(I^(k+1)_k(P_k+F_k(u_k)))$
			calcDampingSensors();
			calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, false);
			SetResidual_Term(RHSrho, RHSrhou, RHSrhoE, MESH_0);
			
			
			// LEVEL 1 - Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$
			if (Simplify_MG == "NO") Set_SaveSolution(1, MESH_1);
			SetRestricted_Solution(MESH_0);
			calcStateVariables_mg(MESH_1); calcMaterialProperties_mg(MESH_1); setBC_mg(MESH_1); calcRansTurbViscMuet_mg(MESH_1);
			calcRhs_mg(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, A_mgLevel1, false, MESH_1);
	    	if (Simplify_MG == "NO") Set_SaveSolution(-1, MESH_1);
			
			// LEVEL 1 - Compute $P_(k+1) = r^*_(k+1) - r_(k+1)
			SetForcing_Term(RHSrho, RHSrhou, RHSrhoE, RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, MESH_1);
			
			// LEVEL 1 - Presmoothing
			for (int iPreSmooth = 0; iPreSmooth < PreSmoothing; iPreSmooth++) {
				calcStateVariables_mg(MESH_1); calcMaterialProperties_mg(MESH_1); calcRansTurbViscMuet_mg(MESH_1); setBC_mg(MESH_1); 
				dt_min_mg = calcDt_mg(cfl, MESH_1);
				calcRhs_mg(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, A_mgLevel1, true, MESH_1);
				calcTimeInt_mg(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, A_mgLevel1, dq_mgLevel1, 
							   rhs_mgLevel1, underRelax, maxIterLS, zeroAbsLS, zeroRelLS, MESH_1);
			}
			
			if (nMesh >= 2) {
				// LEVEL 2 - Compute $r^*_(k+1) = Damp(I^(k+1)_k(P_k+F_k(u_k)))$
				calcStateVariables_mg(MESH_1); calcMaterialProperties_mg(MESH_1); setBC_mg(MESH_1); calcRansTurbViscMuet_mg(MESH_1);
				calcRhs_mg(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, A_mgLevel1, false, MESH_1);
				SetResidual_Term(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, MESH_1);

				// LEVEL 2 - Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$
				if (Simplify_MG == "NO") Set_SaveSolution(1, MESH_2);
				SetRestricted_Solution(MESH_1);
				calcStateVariables_mg(MESH_2); calcMaterialProperties_mg(MESH_2); setBC_mg(MESH_2); calcRansTurbViscMuet_mg(MESH_2);
				calcRhs_mg(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, A_mgLevel2, false, MESH_2);
				if (Simplify_MG == "NO") Set_SaveSolution(-1, MESH_2);

				// LEVEL 2 - Compute $P_(k+1) = r^*_(k+1) - r_(k+1)$
				SetForcing_Term(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, MESH_2);

				// LEVEL 2 - Presmoothing
				for (int iPreSmooth = 0; iPreSmooth < PreSmoothing; iPreSmooth++) {
					calcStateVariables_mg(MESH_2);
					dt_min_mg = calcDt_mg(cfl, MESH_2);
					calcMaterialProperties_mg(MESH_2); setBC_mg(MESH_2); calcRansTurbViscMuet_mg(MESH_2);
					calcRhs_mg(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, A_mgLevel2, true, MESH_2);
					calcTimeInt_mg(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, A_mgLevel2, dq_mgLevel2, 
								   rhs_mgLevel2, underRelax, maxIterLS, zeroAbsLS, zeroRelLS, MESH_2); 
				}

				if (nMesh == 3) {
					// LEVEL 3 - Compute $r^*_(k+1) = Damp(I^(k+1)_k(P_k+F_k(u_k)))$
					calcStateVariables_mg(MESH_2); calcMaterialProperties_mg(MESH_2); setBC_mg(MESH_2); calcRansTurbViscMuet_mg(MESH_2);
					calcRhs_mg(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, A_mgLevel2, false, MESH_2);
					SetResidual_Term(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, MESH_2);
					
					// LEVEL 3 - Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$
					if (Simplify_MG == "NO") Set_SaveSolution(1, MESH_3);
					SetRestricted_Solution(MESH_2);
					calcStateVariables_mg(MESH_3); calcMaterialProperties_mg(MESH_3); setBC_mg(MESH_3); calcRansTurbViscMuet_mg(MESH_3);
					calcRhs_mg(RHSrho_mgLevel3, RHSrhou_mgLevel3, RHSrhoE_mgLevel3, A_mgLevel3, false, MESH_3);
					if (Simplify_MG == "NO") Set_SaveSolution(-1, MESH_3);
					
					// LEVEL 3 - Compute $P_(k+1) = r^*_(k+1) - r_(k+1)
					SetForcing_Term(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, RHSrho_mgLevel3, RHSrhou_mgLevel3, RHSrhoE_mgLevel3, MESH_3);
					
					// LEVEL 3 - Presmoothing
					for (int iPreSmooth = 0; iPreSmooth < PreSmoothing; iPreSmooth++) {
						calcStateVariables_mg(MESH_3);
						dt_min_mg = calcDt_mg(cfl, MESH_3);
						calcMaterialProperties_mg(MESH_3); setBC_mg(MESH_3); calcRansTurbViscMuet_mg(MESH_3);
						calcRhs_mg(RHSrho_mgLevel3, RHSrhou_mgLevel3, RHSrhoE_mgLevel3, A_mgLevel3, true, MESH_3);
						calcTimeInt_mg(RHSrho_mgLevel3, RHSrhou_mgLevel3, RHSrhoE_mgLevel3, A_mgLevel3, dq_mgLevel3, 
									   rhs_mgLevel3, underRelax, maxIterLS, zeroAbsLS, zeroRelLS, MESH_3);
					}
					
					// LEVEL 3 - Compute prolongated solution with implicit residual smoothing $u^(new)_k = u_k +  I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k)$
					SetProlongated_Correction(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, MESH_3);
					SetResidual_Smoothing(RHSrho, RHSrhou, RHSrhoE, AgglomSmoothing, 1.25, MESH_2);
					SetSolution(RHSrho_mgLevel2, RHSrhou_mgLevel2, RHSrhoE_mgLevel2, MESH_2);
				}
				// LEVEL 2 - Compute prolongated solution with implicit residual smoothing $u^(new)_k = u_k +  I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k)$
				SetProlongated_Correction(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, MESH_2);
				SetResidual_Smoothing(RHSrho, RHSrhou, RHSrhoE, AgglomSmoothing, 1.25, MESH_1);
				SetSolution(RHSrho_mgLevel1, RHSrhou_mgLevel1, RHSrhoE_mgLevel1, MESH_1);
			}
			
			// LEVEL 1 - Compute prolongated solution with implicit residual smoothing $u^(new)_k = u_k +  I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k)$
			SetProlongated_Correction(RHSrho, RHSrhou, RHSrhoE, MESH_1);
			SetResidual_Smoothing(RHSrho, RHSrhou, RHSrhoE, AgglomSmoothing, 1.25, MESH_0);
			SetSolution(RHSrho, RHSrhou, RHSrhoE, MESH_0);
		}
		
		
		if (step%check_interval == 0) {
			if ((mpi_rank == 0) && (step%(check_interval*10) == 0)) {
				cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;
			}
			
			showResidue(Residual);
			
		}
		
		
		temporalHook();
		dumpProbes(step, 0.0);
		
		writeData(step);
		
		if ((write_restart > 0) && (step % write_restart == 0))
			writeRestart(step);
		
		
		done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation
		
	}
	
	
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
	{
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
	}
	
	
	// output
	temporalHook();
	finalHook();
	writeRestart();
	
	// delete memory
	MPI_Barrier(mpi_comm);
	delete [] A; delete [] rhs; delete [] RHSrho; delete [] RHSrhou; delete [] RHSrhoE; delete [] dq;
	/*  if (nMesh > 0) { delete [] A_mgLevel1; delete [] rhs_mgLevel1; delete [] RHSrho_mgLevel1; delete [] RHSrhou_mgLevel1; delete [] RHSrhoE_mgLevel1; delete [] dq_mgLevel1; }
	 if (nMesh > 1) { delete [] A_mgLevel2; delete [] rhs_mgLevel2; delete [] RHSrho_mgLevel2; delete [] RHSrhou_mgLevel2; delete [] RHSrhoE_mgLevel2; delete [] dq_mgLevel2; }
	 if (nMesh > 2) { delete [] A_mgLevel3; delete [] rhs_mgLevel3; delete [] RHSrho_mgLevel3; delete [] RHSrhou_mgLevel3; delete [] RHSrhoE_mgLevel3; delete [] dq_mgLevel3; } */
	
	if (nScal > 0) freeMem3D(AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1);
	if (nScal > 0) freeMem2D(dScal,   0, nScal-1, 0, ncv_g-1);
	if (nScal > 0) freeMem2D(rhsScal, 0, nScal-1, 0, ncv-1);
	
	MPI_Finalize();
	exit(1);
	
}
	
	
void JoeWithModels::calcTimeInt_mg(double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double (*A)[5][5], double (*dq)[5], 
								   double (*rhs)[5], double underRelax, int maxIterLS, double zeroAbsLS, double zeroRelLS, int iMesh)
{
	double CFLRelax = getDoubleParam("CFL_RELAXATION", "0.75");
	
	if (iMesh == 1) {
		for (int icv = 0; icv < ncv_g_mgLevel1; icv++)
			for (int i = 0; i < 5; i++)
				dq[icv][i] = 0.0;
		
		// solve linear system for the NSE
		for (int icv=0; icv<ncv_mgLevel1; ++icv) {
			rhs[icv][0] = underRelax*(RHSrho[icv]+rho_TruncError_mgLevel1[icv]);
			rhs[icv][1] = underRelax*(RHSrhou[icv][0]+rhou_TruncError_mgLevel1[icv][0]);
			rhs[icv][2] = underRelax*(RHSrhou[icv][1]+rhou_TruncError_mgLevel1[icv][1]);
			rhs[icv][3] = underRelax*(RHSrhou[icv][2]+rhou_TruncError_mgLevel1[icv][2]);
			rhs[icv][4] = underRelax*(RHSrhoE[icv]+rhoE_TruncError_mgLevel1[icv]);
			
			double tmp = pow(double (iMesh), CFLRelax)*(cv_volume_mgLevel1[icv]/(local_dt_mgLevel1[icv]));
			for (int i = 0; i < 5; i++)
				A[nbocv_i_mgLevel1[icv]][i][i] += tmp;
		}
		
		// solve linear system
 		solveCoupledLinSysNS_mg(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, iMesh);  
		
		for (int icv=0; icv<ncv_mgLevel1; icv++) {
			rho_mgLevel1[icv]     += dq[icv][0];
			rhou_mgLevel1[icv][0] += dq[icv][1];
			rhou_mgLevel1[icv][1] += dq[icv][2];
			rhou_mgLevel1[icv][2] += dq[icv][3];
			rhoE_mgLevel1[icv]    += dq[icv][4];
			
		}
		
		// Update the solution at the ghost cells
//		updateCvData_mg(rho_mgLevel1,  REPLACE_DATA, MESH_1);
//		updateCvData_mg(rhou_mgLevel1, REPLACE_ROTATE_DATA, MESH_1);
//		updateCvData_mg(rhoE_mgLevel1, REPLACE_DATA, MESH_1);
		
	}
	
	if (iMesh == 2) {
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++)
			for (int i = 0; i < 5; i++)
				dq[icv][i] = 0.0;
		
		// solve linear system for the NSE
		for (int icv=0; icv<ncv_mgLevel2; ++icv) {
			rhs[icv][0] = underRelax*(RHSrho[icv]+rho_TruncError_mgLevel2[icv]);
			rhs[icv][1] = underRelax*(RHSrhou[icv][0]+rhou_TruncError_mgLevel2[icv][0]);
			rhs[icv][2] = underRelax*(RHSrhou[icv][1]+rhou_TruncError_mgLevel2[icv][1]);
			rhs[icv][3] = underRelax*(RHSrhou[icv][2]+rhou_TruncError_mgLevel2[icv][2]);
			rhs[icv][4] = underRelax*(RHSrhoE[icv]+rhoE_TruncError_mgLevel2[icv]);
			
			double tmp = pow(double (iMesh), CFLRelax)*(cv_volume_mgLevel2[icv]/(local_dt_mgLevel2[icv]));
			for (int i = 0; i < 5; i++)
				A[nbocv_i_mgLevel2[icv]][i][i] += tmp; 
			
			
		}
		
		// solve linear system
		solveCoupledLinSysNS_mg(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, iMesh);  
		for (int icv=0; icv<ncv_mgLevel2; icv++) {
			rho_mgLevel2[icv]     += dq[icv][0];
			rhou_mgLevel2[icv][0] += dq[icv][1];
			rhou_mgLevel2[icv][1] += dq[icv][2];
			rhou_mgLevel2[icv][2] += dq[icv][3];
			rhoE_mgLevel2[icv]    += dq[icv][4];
		}
		
		// Update the solution at the ghost cells
//		updateCvData_mg(rho_mgLevel2,  REPLACE_DATA, MESH_2);
//		updateCvData_mg(rhou_mgLevel2, REPLACE_ROTATE_DATA, MESH_2);
//		updateCvData_mg(rhoE_mgLevel2, REPLACE_DATA, MESH_2);
	}
	
	if (iMesh == 3) {
		for (int icv = 0; icv < ncv_g_mgLevel3; icv++)
			for (int i = 0; i < 5; i++)
				dq[icv][i] = 0.0;
		
		// solve linear system for the NSE
		for (int icv=0; icv<ncv_mgLevel3; ++icv) {
			rhs[icv][0] = underRelax*(RHSrho[icv]+rho_TruncError_mgLevel3[icv]);
			rhs[icv][1] = underRelax*(RHSrhou[icv][0]+rhou_TruncError_mgLevel3[icv][0]);
			rhs[icv][2] = underRelax*(RHSrhou[icv][1]+rhou_TruncError_mgLevel3[icv][1]);
			rhs[icv][3] = underRelax*(RHSrhou[icv][2]+rhou_TruncError_mgLevel3[icv][2]);
			rhs[icv][4] = underRelax*(RHSrhoE[icv]+rhoE_TruncError_mgLevel3[icv]);
			
			double tmp = pow(double (iMesh), CFLRelax)*(cv_volume_mgLevel3[icv]/(local_dt_mgLevel3[icv]));
			for (int i = 0; i < 5; i++)
				A[nbocv_i_mgLevel3[icv]][i][i] += tmp;
		}
		
		// solve linear system
		solveCoupledLinSysNS_mg(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, iMesh);  
		for (int icv=0; icv<ncv_mgLevel3; icv++) {
			rho_mgLevel3[icv]     += dq[icv][0];
			rhou_mgLevel3[icv][0] += dq[icv][1];
			rhou_mgLevel3[icv][1] += dq[icv][2];
			rhou_mgLevel3[icv][2] += dq[icv][3];
			rhoE_mgLevel3[icv]    += dq[icv][4];
		}
		
		// Update the solution at the ghost cells
//		updateCvData_mg(rho_mgLevel3,  REPLACE_DATA, MESH_3);
//		updateCvData_mg(rhou_mgLevel3, REPLACE_ROTATE_DATA, MESH_3);
//		updateCvData_mg(rhoE_mgLevel3, REPLACE_DATA, MESH_3);
	}
	
}

void JoeWithModels::calcRhs_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh)
{
	
	if (iMesh == 1) {
		// Compute RHS for both NSE and scalars
		for (int noc = 0; noc < nbocv_s_mgLevel1; noc++)
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					A[noc][i][j] = 0.0;
		
		// set RHS to zero
		for (int icv = 0; icv < ncv_mgLevel1; icv++) {
			rhs_rho[icv] = 0.0;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv][i] = 0.0;
			rhs_rhoE[icv] = 0.0;
		}
	}
	
	if (iMesh == 2) {
		// Compute RHS for both NSE and scalars
		for (int noc = 0; noc < nbocv_s_mgLevel2; noc++)
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					A[noc][i][j] = 0.0;
		
		// set RHS to zero
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			rhs_rho[icv] = 0.0;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv][i] = 0.0;
			rhs_rhoE[icv] = 0.0;
		}	
	}
	
	if (iMesh == 3) {
		// Compute RHS for both NSE and scalars
		for (int noc = 0; noc < nbocv_s_mgLevel3; noc++)
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					A[noc][i][j] = 0.0;
		
		// set RHS to zero
		for (int icv = 0; icv < ncv_mgLevel3; icv++) {
			rhs_rho[icv] = 0.0;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv][i] = 0.0;
			rhs_rhoE[icv] = 0.0;
		}	
	}
	
  // Compute Euler Flux for NS and scalars
  calcEulerFlux_mg(rhs_rho, rhs_rhou, rhs_rhoE, A, flagImplicit, iMesh);
	
  // Compute viscous Flux for NS
	if (mu_ref > 0.0)
    calcViscousFluxNS_mg(rhs_rho, rhs_rhou, rhs_rhoE, A, flagImplicit, iMesh);
	
}

void JoeWithModels::calcEulerFlux_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh) {
	
	string SpaceIntName = getStringParam("SPACE_INTEGRATION","HLLC");
	
	double (*Apl)[5] = NULL;
	double (*Ami)[5] = NULL;
	
	if (flagImplicit) {
		Apl = new double[5][5];
		Ami = new double[5][5];
		
		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Apl[i][j] = Ami[i][j] = 0.0;
	}
	
	double Frho, Frhou[3], FrhoE;
	
	// If JST scheme, it is necessary to do some preprocessing
	if (SpaceIntName == "JST") calcJSTCoeff_Var(iMesh);
	
	// count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
	int CountReducedOrder = 0;
	int myCountReducedOrder = 0;
	
	
	if (iMesh == 1) {
		
		// cycle through internal faces, assembling flux to both sides
		for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
			
			int icv0 = cvofa_mgLevel1[ifa][0];
			int icv1 = cvofa_mgLevel1[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
			
			// reconstruction of variables at faces: rho, u, T or P, scalars
			double rho0 = rho_mgLevel1[icv0];
			double u0[3] = {vel_mgLevel1[icv0][0], vel_mgLevel1[icv0][1], vel_mgLevel1[icv0][2]};
			double p0 = press_mgLevel1[icv0];
			double T0 = temp_mgLevel1[icv0];
			double h0 = enthalpy_mgLevel1[icv0];
			double gam0 = gamma_mgLevel1[icv0];
			double R0 = RoM_mgLevel1[icv0];
			double kineCV0 = 0.0;
			double kineFA0 = 0.0;
			
			double rho1 = rho_mgLevel1[icv1];
			double u1[3] = {vel_mgLevel1[icv1][0], vel_mgLevel1[icv1][1], vel_mgLevel1[icv1][2]};
			double p1 = press_mgLevel1[icv1];
			double T1 = temp_mgLevel1[icv1];
			double h1 = enthalpy_mgLevel1[icv1];
			double gam1 = gamma_mgLevel1[icv1];
			double R1 = RoM_mgLevel1[icv1];
			double kineCV1 = 0.0;
			double kineFA1 = 0.0;
			
			// calculation of Euler Flux explicit using HLLC
			if (SpaceIntName == "HLLC")
				calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "AUSM")
				calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "ROE")
				calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "JST") {
				double Lambda0 = Lambda_mgLevel1[icv0]; double Lambda1 = Lambda_mgLevel1[icv1];
				int Neighbor0 = NeighborCV_mgLevel1[icv0]; double Neighbor1 = NeighborCV_mgLevel1[icv1];
				calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
								  rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0,
								  rho1, u1, p1, T1, h1, R1, gam1, NULL, kineFA1,
								  area, nVec, 0, 0.0,
								  Lambda0, Neighbor0, Lambda1, Neighbor1);
			}
			
			// icv0 is always valid...
			rhs_rho[icv0] -= Frho;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel1) {
				rhs_rho[icv1] += Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
			
			// calculate implicit matrix using HLLC
			if (flagImplicit) {
				if (SpaceIntName == "HLLC")
					calcEulerFluxMatrices_HLLC(Apl, Ami, NULL, NULL,
											   rho_mgLevel1[icv0], vel_mgLevel1[icv0], press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0], NULL, kineCV0,
											   rho_mgLevel1[icv1], vel_mgLevel1[icv1], press_mgLevel1[icv1], temp_mgLevel1[icv1], enthalpy_mgLevel1[icv1], RoM_mgLevel1[icv1], gamma_mgLevel1[icv1], NULL, kineCV1,
											   area, nVec, 0, 0.0);
				if (SpaceIntName == "ROE")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel1[icv0], vel_mgLevel1[icv0], press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0], NULL, kineCV0,
											  rho_mgLevel1[icv1], vel_mgLevel1[icv1], press_mgLevel1[icv1], temp_mgLevel1[icv1], enthalpy_mgLevel1[icv1], RoM_mgLevel1[icv1], gamma_mgLevel1[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "AUSM")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel1[icv0], vel_mgLevel1[icv0], press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0], NULL, kineCV0,
											  rho_mgLevel1[icv1], vel_mgLevel1[icv1], press_mgLevel1[icv1], temp_mgLevel1[icv1], enthalpy_mgLevel1[icv1], RoM_mgLevel1[icv1], gamma_mgLevel1[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "JST")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel1[icv0], vel_mgLevel1[icv0], press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0], NULL, kineCV0,
											  rho_mgLevel1[icv1], vel_mgLevel1[icv1], press_mgLevel1[icv1], temp_mgLevel1[icv1], enthalpy_mgLevel1[icv1], RoM_mgLevel1[icv1], gamma_mgLevel1[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++) {
						A[noc00][i][j] += Apl[i][j];
						A[noc01][i][j] += Ami[i][j];
					}
				
				if (icv1 < ncv_mgLevel1) { // if icv1 is internal...
					
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc11][i][j] -= Ami[i][j];
							A[noc10][i][j] -= Apl[i][j];
						}
				}
			}
		}
		
		
		// cycle through boundary faces, assembling flux
		for (list<FaZone>::iterator zone = faZoneList_mgLevel1.begin(); zone != faZoneList_mgLevel1.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName()))  {
					
					// SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
					if (param->getString() == "SYMMETRY")  {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)  {
							
							
							int icv0 = cvofa_mgLevel1[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							double kineFA = 0.0;
							
							rhs_rho[icv0] -= 0.0;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= p_bfa_mgLevel1[ifa]*nVec[i]*area;
							rhs_rhoE[icv0] -= 0.0;							
							
							if (flagImplicit)
							{
								double SimmJac[5][5];
								double a2 = gam_bfa_mgLevel1[ifa]-1.0;
								double phi = 0.5*a2*(vel_bfa_mgLevel1[ifa][0]*vel_bfa_mgLevel1[ifa][0]+vel_bfa_mgLevel1[ifa][1]*vel_bfa_mgLevel1[ifa][1]+vel_bfa_mgLevel1[ifa][2]*vel_bfa_mgLevel1[ifa][2]);
								for (int i = 0; i<5; i++) {
									SimmJac[0][i] = 0.0;
									SimmJac[4][i] = 0.0;
								}
								for (int i = 0; i<3; i++) {
									SimmJac[i+1][0] = phi*nVec[i]*area;
									for (int j = 0; j<3; j++)
										SimmJac[i+1][j+1] = -a2*vel_bfa_mgLevel1[ifa][j]*nVec[i]*area;
									SimmJac[i+1][4] = a2*nVec[i]*area;
								}
								
								int noc00 = nbocv_i_mgLevel1[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += SimmJac[i][j];
							}
						}
					}
					
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel1[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST")
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0,
												  Lambda_mgLevel1[icv0], NeighborCV_mgLevel1[icv0], Lambda_mgLevel1[icv0], NeighborCV_mgLevel1[icv0]);
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit)
							{
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
																						 rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						 rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						 area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
																						rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
																						rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
																						area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
															  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel1[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
					
					// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel1[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							double rho0 = rho_mgLevel1[icv0];
							double u0[3] = {vel_mgLevel1[icv0][0], vel_mgLevel1[icv0][1], vel_mgLevel1[icv0][2]};
							double p0 = press_mgLevel1[icv0];
							double T0 = temp_mgLevel1[icv0];
							double h0 = enthalpy_mgLevel1[icv0];
							double gam0 = gamma_mgLevel1[icv0];
							double R0 = RoM_mgLevel1[icv0];
							double kineCV0 = 0.0;           // cell center
							double kineFA0 = 0.0;           // cell face
							double kineFA1 = 0.0;
							
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST") {
								double Lambda0 = Lambda_mgLevel1[icv0]; double Lambda1 = Lambda_mgLevel1[icv0];
								int Neighbor0 = NeighborCV_mgLevel1[icv0]; double Neighbor1 = NeighborCV_mgLevel1[icv0];
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa], T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0,
												  Lambda0, Neighbor0, Lambda1, Neighbor1);
							}
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit) {
								
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
															   rho_mgLevel1[icv0],    vel_mgLevel1[icv0],    press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0],    gamma_mgLevel1[icv0],  NULL, kineCV0,
															   rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa],  T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa],     RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
															   area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel1[icv0],    vel_mgLevel1[icv0],    press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0],    gamma_mgLevel1[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa],  T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa],     RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel1[icv0],    vel_mgLevel1[icv0],    press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0],    gamma_mgLevel1[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa],  T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa],     RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel1[icv0],    vel_mgLevel1[icv0],    press_mgLevel1[icv0], temp_mgLevel1[icv0], enthalpy_mgLevel1[icv0], RoM_mgLevel1[icv0],    gamma_mgLevel1[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], p_bfa_mgLevel1[ifa],  T_bfa_mgLevel1[ifa], h_bfa_mgLevel1[ifa],     RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel1[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
				}
			}
		}
	}
	
	if (iMesh == 2) {
		
		// cycle through internal faces, assembling flux to both sides
		for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
			
			int icv0 = cvofa_mgLevel2[ifa][0];
			int icv1 = cvofa_mgLevel2[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
			
			// reconstruction of variables at faces: rho, u, T or P, scalars
			double rho0 = rho_mgLevel2[icv0];
			double u0[3] = {vel_mgLevel2[icv0][0], vel_mgLevel2[icv0][1], vel_mgLevel2[icv0][2]};
			double p0 = press_mgLevel2[icv0];
			double T0 = temp_mgLevel2[icv0];
			double h0 = enthalpy_mgLevel2[icv0];
			double gam0 = gamma_mgLevel2[icv0];
			double R0 = RoM_mgLevel2[icv0];
			double kineCV0 = 0.0;
			double kineFA0 = 0.0;
			
			double rho1 = rho_mgLevel2[icv1];
			double u1[3] = {vel_mgLevel2[icv1][0], vel_mgLevel2[icv1][1], vel_mgLevel2[icv1][2]};
			double p1 = press_mgLevel2[icv1];
			double T1 = temp_mgLevel2[icv1];
			double h1 = enthalpy_mgLevel2[icv1];
			double gam1 = gamma_mgLevel2[icv1];
			double R1 = RoM_mgLevel2[icv1];
			double kineCV1 = 0.0;
			double kineFA1 = 0.0;
			
			// calculation of Euler Flux explicit using HLLC
			if (SpaceIntName == "HLLC")
				calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "AUSM")
				calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "ROE")
				calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "JST") {
				double Lambda0 = Lambda_mgLevel2[icv0]; double Lambda1 = Lambda_mgLevel2[icv1];
				int Neighbor0 = NeighborCV_mgLevel2[icv0]; double Neighbor1 = NeighborCV_mgLevel2[icv1];
				calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
								  rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0,
								  rho1, u1, p1, T1, h1, R1, gam1, NULL, kineFA1,
								  area, nVec, 0, 0.0,
								  Lambda0, Neighbor0, Lambda1, Neighbor1);
			}
			
			// icv0 is always valid...
			rhs_rho[icv0] -= Frho;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel2) {
				rhs_rho[icv1] += Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
			
			// calculate implicit matrix using HLLC
			if (flagImplicit) {
				if (SpaceIntName == "HLLC")
					calcEulerFluxMatrices_HLLC(Apl, Ami, NULL, NULL,
											   rho_mgLevel2[icv0], vel_mgLevel2[icv0], press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0], gamma_mgLevel2[icv0], NULL, kineCV0,
											   rho_mgLevel2[icv1], vel_mgLevel2[icv1], press_mgLevel2[icv1], temp_mgLevel2[icv1], enthalpy_mgLevel2[icv1], RoM_mgLevel2[icv1], gamma_mgLevel2[icv1], NULL, kineCV1,
											   area, nVec, 0, 0.0);
				if (SpaceIntName == "ROE")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel2[icv0], vel_mgLevel2[icv0], press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0], gamma_mgLevel2[icv0], NULL, kineCV0,
											  rho_mgLevel2[icv1], vel_mgLevel2[icv1], press_mgLevel2[icv1], temp_mgLevel2[icv1], enthalpy_mgLevel2[icv1], RoM_mgLevel2[icv1], gamma_mgLevel2[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "AUSM")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel2[icv0], vel_mgLevel2[icv0], press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0], gamma_mgLevel2[icv0], NULL, kineCV0,
											  rho_mgLevel2[icv1], vel_mgLevel2[icv1], press_mgLevel2[icv1], temp_mgLevel2[icv1], enthalpy_mgLevel2[icv1], RoM_mgLevel2[icv1], gamma_mgLevel2[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "JST")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel2[icv0], vel_mgLevel2[icv0], press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0], gamma_mgLevel2[icv0], NULL, kineCV0,
											  rho_mgLevel2[icv1], vel_mgLevel2[icv1], press_mgLevel2[icv1], temp_mgLevel2[icv1], enthalpy_mgLevel2[icv1], RoM_mgLevel2[icv1], gamma_mgLevel2[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++) {
						A[noc00][i][j] += Apl[i][j];
						A[noc01][i][j] += Ami[i][j];
					}
				
				if (icv1 < ncv_mgLevel2) { // if icv1 is internal...
					
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc11][i][j] -= Ami[i][j];
							A[noc10][i][j] -= Apl[i][j];
						}
				}
			}
		}
		
		
		// cycle through boundary faces, assembling flux
		for (list<FaZone>::iterator zone = faZoneList_mgLevel2.begin(); zone != faZoneList_mgLevel2.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName()))  {
					
					// SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
					if (param->getString() == "SYMMETRY")  {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)  {
							
							
							int icv0 = cvofa_mgLevel2[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							double kineFA = 0.0;
							
							rhs_rho[icv0] -= 0.0;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= p_bfa_mgLevel2[ifa]*nVec[i]*area;
							rhs_rhoE[icv0] -= 0.0;							
							
							if (flagImplicit)
							{
								double SimmJac[5][5];
								double a2 = gam_bfa_mgLevel2[ifa]-1.0;
								double phi = 0.5*a2*(vel_bfa_mgLevel2[ifa][0]*vel_bfa_mgLevel2[ifa][0]+vel_bfa_mgLevel2[ifa][1]*vel_bfa_mgLevel2[ifa][1]+vel_bfa_mgLevel2[ifa][2]*vel_bfa_mgLevel2[ifa][2]);
								for (int i = 0; i<5; i++) {
									SimmJac[0][i] = 0.0;
									SimmJac[4][i] = 0.0;
								}
								for (int i = 0; i<3; i++) {
									SimmJac[i+1][0] = phi*nVec[i]*area;
									for (int j = 0; j<3; j++)
										SimmJac[i+1][j+1] = -a2*vel_bfa_mgLevel2[ifa][j]*nVec[i]*area;
									SimmJac[i+1][4] = a2*nVec[i]*area;
								}
								
								int noc00 = nbocv_i_mgLevel2[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += SimmJac[i][j];
							}
						}
					}
					
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel2[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST")
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0,
												  Lambda_mgLevel2[icv0], NeighborCV_mgLevel2[icv0], Lambda_mgLevel2[icv0], NeighborCV_mgLevel2[icv0]);
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit)
							{
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
															   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															   area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel2[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
					
					// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel2[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							double rho0 = rho_mgLevel2[icv0];
							double u0[3] = {vel_mgLevel2[icv0][0], vel_mgLevel2[icv0][1], vel_mgLevel2[icv0][2]};
							double p0 = press_mgLevel2[icv0];
							double T0 = temp_mgLevel2[icv0];
							double h0 = enthalpy_mgLevel2[icv0];
							double gam0 = gamma_mgLevel2[icv0];
							double R0 = RoM_mgLevel2[icv0];
							double kineCV0 = 0.0;           // cell center
							double kineFA0 = 0.0;           // cell face
							double kineFA1 = 0.0;
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST") {
								double Lambda0 = Lambda_mgLevel2[icv0]; double Lambda1 = Lambda_mgLevel2[icv0];
								int Neighbor0 = NeighborCV_mgLevel2[icv0]; double Neighbor1 = NeighborCV_mgLevel2[icv0];
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa], T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0,
												  Lambda0, Neighbor0, Lambda1, Neighbor1);
							}
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit) {
								
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
															   rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  NULL, kineCV0,
															   rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa],  T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa],     RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
															   area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa],  T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa],     RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa],  T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa],     RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    press_mgLevel2[icv0], temp_mgLevel2[icv0], enthalpy_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], p_bfa_mgLevel2[ifa],  T_bfa_mgLevel2[ifa], h_bfa_mgLevel2[ifa],     RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel2[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
				}
			}
		}
		 
	}
	
	if (iMesh == 3) {
		
		// cycle through internal faces, assembling flux to both sides
		for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
			
			int icv0 = cvofa_mgLevel3[ifa][0];
			int icv1 = cvofa_mgLevel3[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
			
			// reconstruction of variables at faces: rho, u, T or P, scalars
			double rho0 = rho_mgLevel3[icv0];
			double u0[3] = {vel_mgLevel3[icv0][0], vel_mgLevel3[icv0][1], vel_mgLevel3[icv0][2]};
			double p0 = press_mgLevel3[icv0];
			double T0 = temp_mgLevel3[icv0];
			double h0 = enthalpy_mgLevel3[icv0];
			double gam0 = gamma_mgLevel3[icv0];
			double R0 = RoM_mgLevel3[icv0];
			double kineCV0 = 0.0;
			double kineFA0 = 0.0;
			
			double rho1 = rho_mgLevel3[icv1];
			double u1[3] = {vel_mgLevel3[icv1][0], vel_mgLevel3[icv1][1], vel_mgLevel3[icv1][2]};
			double p1 = press_mgLevel3[icv1];
			double T1 = temp_mgLevel3[icv1];
			double h1 = enthalpy_mgLevel3[icv1];
			double gam1 = gamma_mgLevel3[icv1];
			double R1 = RoM_mgLevel3[icv1];
			double kineCV1 = 0.0;
			double kineFA1 = 0.0;
			
			// calculation of Euler Flux explicit using HLLC
			if (SpaceIntName == "HLLC")
				calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "AUSM")
				calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													 p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "ROE")
				calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL, rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0, rho1, u1, 
													p1, T1, h1, R1, gam1, NULL, kineFA1, area, nVec, 0, 0.0);
			if (SpaceIntName == "JST") {
				double Lambda0 = Lambda_mgLevel3[icv0]; double Lambda1 = Lambda_mgLevel3[icv1];
				int Neighbor0 = NeighborCV_mgLevel3[icv0]; double Neighbor1 = NeighborCV_mgLevel3[icv1];
				calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
								  rho0, u0, p0, T0, h0, R0, gam0, NULL, kineFA0,
								  rho1, u1, p1, T1, h1, R1, gam1, NULL, kineFA1,
								  area, nVec, 0, 0.0,
								  Lambda0, Neighbor0, Lambda1, Neighbor1);
			}
			
			// icv0 is always valid...
			rhs_rho[icv0] -= Frho;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel3) {
				rhs_rho[icv1] += Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
			
			// calculate implicit matrix using HLLC
			if (flagImplicit) {
				if (SpaceIntName == "HLLC")
					calcEulerFluxMatrices_HLLC(Apl, Ami, NULL, NULL,
											   rho_mgLevel3[icv0], vel_mgLevel3[icv0], press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0], gamma_mgLevel3[icv0], NULL, kineCV0,
											   rho_mgLevel3[icv1], vel_mgLevel3[icv1], press_mgLevel3[icv1], temp_mgLevel3[icv1], enthalpy_mgLevel3[icv1], RoM_mgLevel3[icv1], gamma_mgLevel3[icv1], NULL, kineCV1,
											   area, nVec, 0, 0.0);
				if (SpaceIntName == "ROE")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel3[icv0], vel_mgLevel3[icv0], press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0], gamma_mgLevel3[icv0], NULL, kineCV0,
											  rho_mgLevel3[icv1], vel_mgLevel3[icv1], press_mgLevel3[icv1], temp_mgLevel3[icv1], enthalpy_mgLevel3[icv1], RoM_mgLevel3[icv1], gamma_mgLevel3[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "AUSM")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel3[icv0], vel_mgLevel3[icv0], press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0], gamma_mgLevel3[icv0], NULL, kineCV0,
											  rho_mgLevel3[icv1], vel_mgLevel3[icv1], press_mgLevel3[icv1], temp_mgLevel3[icv1], enthalpy_mgLevel3[icv1], RoM_mgLevel3[icv1], gamma_mgLevel3[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				if (SpaceIntName == "JST")
					calcEulerFluxMatrices_Roe(Apl, Ami, NULL, NULL,
											  rho_mgLevel3[icv0], vel_mgLevel3[icv0], press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0], gamma_mgLevel3[icv0], NULL, kineCV0,
											  rho_mgLevel3[icv1], vel_mgLevel3[icv1], press_mgLevel3[icv1], temp_mgLevel3[icv1], enthalpy_mgLevel3[icv1], RoM_mgLevel3[icv1], gamma_mgLevel3[icv1], NULL, kineCV1,
											  area, nVec, 0, 0.0);
				
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++) {
						A[noc00][i][j] += Apl[i][j];
						A[noc01][i][j] += Ami[i][j];
					}
				
				if (icv1 < ncv_mgLevel3) { // if icv1 is internal...
					
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc11][i][j] -= Ami[i][j];
							A[noc10][i][j] -= Apl[i][j];
						}
				}
			}
		}
		
		// cycle through boundary faces, assembling flux
		for (list<FaZone>::iterator zone = faZoneList_mgLevel3.begin(); zone != faZoneList_mgLevel3.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName()))  {
					
					// SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
					if (param->getString() == "SYMMETRY")  {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)  {
							
							
							int icv0 = cvofa_mgLevel3[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							double kineFA = 0.0;
							
							rhs_rho[icv0] -= 0.0;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= p_bfa_mgLevel3[ifa]*nVec[i]*area;
							rhs_rhoE[icv0] -= 0.0;							
							
							if (flagImplicit)
							{
								double SimmJac[5][5];
								double a2 = gam_bfa_mgLevel3[ifa]-1.0;
								double phi = 0.5*a2*(vel_bfa_mgLevel3[ifa][0]*vel_bfa_mgLevel3[ifa][0]+vel_bfa_mgLevel3[ifa][1]*vel_bfa_mgLevel3[ifa][1]+vel_bfa_mgLevel3[ifa][2]*vel_bfa_mgLevel3[ifa][2]);
								for (int i = 0; i<5; i++) {
									SimmJac[0][i] = 0.0;
									SimmJac[4][i] = 0.0;
								}
								for (int i = 0; i<3; i++) {
									SimmJac[i+1][0] = phi*nVec[i]*area;
									for (int j = 0; j<3; j++)
										SimmJac[i+1][j+1] = -a2*vel_bfa_mgLevel3[ifa][j]*nVec[i]*area;
									SimmJac[i+1][4] = a2*nVec[i]*area;
								}
								
								int noc00 = nbocv_i_mgLevel3[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += SimmJac[i][j];
							}
						}
					}
					
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel3[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST")
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
												  area, nVec, 0, 0.0,
												  Lambda_mgLevel3[icv0], NeighborCV_mgLevel3[icv0], Lambda_mgLevel3[icv0], NeighborCV_mgLevel3[icv0]);
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit)
							{
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
															   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															   area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, 0.0,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel3[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
					
					// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							
							int icv0 = cvofa_mgLevel3[ifa][0];
							
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							double rho0 = rho_mgLevel3[icv0];
							double u0[3] = {vel_mgLevel3[icv0][0], vel_mgLevel3[icv0][1], vel_mgLevel3[icv0][2]};
							double p0 = press_mgLevel3[icv0];
							double T0 = temp_mgLevel3[icv0];
							double h0 = enthalpy_mgLevel3[icv0];
							double gam0 = gamma_mgLevel3[icv0];
							double R0 = RoM_mgLevel3[icv0];
							double kineCV0 = 0.0;           // cell center
							double kineFA0 = 0.0;           // cell face
							double kineFA1 = 0.0;
							
							
							if (SpaceIntName == "HLLC")
								calcEulerFlux_HLLC(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "AUSM")
								calcEulerFlux_AUSM(Frho, Frhou, FrhoE, NULL,
												   rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
												   area, nVec, 0, 0.0);
							if (SpaceIntName == "ROE")
								calcEulerFlux_Roe(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0);
							if (SpaceIntName == "JST") {
								double Lambda0 = Lambda_mgLevel3[icv0]; double Lambda1 = Lambda_mgLevel3[icv0];
								int Neighbor0 = NeighborCV_mgLevel3[icv0]; double Neighbor1 = NeighborCV_mgLevel3[icv0];
								calcEulerFlux_Lax(Frho, Frhou, FrhoE, NULL,
												  rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         NULL, kineFA0,
												  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa], T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
												  area, nVec, 0, 0.0,
												  Lambda0, Neighbor0, Lambda1, Neighbor1);
							}
							
							rhs_rho[icv0] -= Frho;
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
							if (flagImplicit) {
								
								if (SpaceIntName == "HLLC")
									calcEulerFluxMatrices_HLLC(Apl, NULL, NULL, NULL,
															   rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  NULL, kineCV0,
															   rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa],  T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa],     RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
															   area, nVec, 0, 0.0);
								if (SpaceIntName == "ROE")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa],  T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa],     RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "AUSM")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa],  T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa],     RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								if (SpaceIntName == "JST")
									calcEulerFluxMatrices_Roe(Apl, NULL, NULL, NULL,
															  rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    press_mgLevel3[icv0], temp_mgLevel3[icv0], enthalpy_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  NULL, kineCV0,
															  rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], p_bfa_mgLevel3[ifa],  T_bfa_mgLevel3[ifa], h_bfa_mgLevel3[ifa],     RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], NULL, kineFA1,
															  area, nVec, 0, 0.0);
								
								int noc00 = nbocv_i_mgLevel3[icv0]; // icv0's diagonal
								for (int i = 0; i < 5; i++)
									for (int j = 0; j < 5; j++)
										A[noc00][i][j] += Apl[i][j];
								
							}
						}
					}
				}
			}
		}
	}
	
	
  if (Apl != NULL)  delete [] Apl;
  if (Ami != NULL)  delete [] Ami;
	
}

void JoeWithModels::calcViscousFluxNS_mg(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit, int iMesh)
{
  double (*A0)[5];
  double (*A1)[5];
	
  if (flagImplicit) {
    A0 = new double[5][5];
    A1 = new double[5][5];
		
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        A0[i][j] = A1[i][j] = 0.0;
  }
  else A0 = A1 = NULL;
	
  double Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;
	
	if (iMesh == 1) {
		
		// compute gradients, with boundary values
		calcCvVectorGradGreenGauss_mg(grad_u_mgLevel1, vel_mgLevel1, vel_bfa_mgLevel1, iMesh);
		calcCvScalarGradGreenGauss_mg(grad_enthalpy_mgLevel1, enthalpy_mgLevel1, h_bfa_mgLevel1, iMesh);
		
		// cycle through internal faces, assembling flux and matrix
		for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
			
			int icv0 = cvofa_mgLevel1[ifa][0];
			int icv1 = cvofa_mgLevel1[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
			
			double uAux_fa[3] = { 0.5*(vel_mgLevel1[icv0][0]+ vel_mgLevel1[icv1][0]),
				0.5*(vel_mgLevel1[icv0][1]+ vel_mgLevel1[icv1][1]),
				0.5*(vel_mgLevel1[icv0][2]+ vel_mgLevel1[icv1][2])};
			
			// kinetic energy, if defined
			double kine0 = 0.0;
			double kine1 = 0.0;
			double kine_fa = 0.0;
			
			// calculate viscous flux
			addViscFlux_mg(Frhou, FrhoE, A0, A1,
										 rho_mgLevel1[icv0], vel_mgLevel1[icv0], grad_u_mgLevel1[icv0], enthalpy_mgLevel1[icv0], 
										 grad_enthalpy_mgLevel1[icv0], temp_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0], kine0,
										 rho_mgLevel1[icv1], vel_mgLevel1[icv1], grad_u_mgLevel1[icv1], enthalpy_mgLevel1[icv1], 
										 grad_enthalpy_mgLevel1[icv1], temp_mgLevel1[icv1], RoM_mgLevel1[icv1], gamma_mgLevel1[icv1], kine1,
										 mul_fa_mgLevel1[ifa], mut_fa_mgLevel1[ifa], lamOcp_fa_mgLevel1[ifa], kine_fa, uAux_fa, area, nVec);
			
			if (flagImplicit) {
				for (int i=1; i<5; i++)
					for (int j=0; j<5; j++) {
						A[noc00][i][j] -= A0[i][j];
						A[noc01][i][j] -= A1[i][j];
					}
				
				if (icv1 < ncv_mgLevel1) {
					for (int i=1; i<5; i++)
						for (int j=0; j<5; j++) {
							A[noc11][i][j] += A1[i][j];
							A[noc10][i][j] += A0[i][j];
						}
				}
			}
			
			// icv0 is always valid...
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel1) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
		}
		
		
		// cycle through boundary faces, assembling flux and matrix
		for (list<FaZone>::iterator zone = faZoneList_mgLevel1.begin(); zone != faZoneList_mgLevel1.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				if (getParam(param, zone->getName())) {
					// SYMMETRY BOUNDARY CONDITION
					if (param->getString() == "SYMMETRY") {
						// No viscous flux in this case
					}
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel1[icv0], vel_mgLevel1[icv0],    grad_u_mgLevel1[icv0], enthalpy_mgLevel1[icv0], grad_enthalpy_mgLevel1[icv0], temp_mgLevel1[icv0], RoM_mgLevel1[icv0], gamma_mgLevel1[icv0],  kine0,
														 rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], grad_u_mgLevel1[icv0], h_bfa_mgLevel1[ifa], grad_enthalpy_mgLevel1[icv0], T_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], kine1,
														 mul_fa_mgLevel1[ifa], 0.0, lamOcp_fa_mgLevel1[ifa], kine_fa, vel_bfa_mgLevel1[ifa], area, nVec); 
							
							if (flagImplicit) {
								int noc00 = nbocv_i_mgLevel1[icv0]; // icv0's diagonal
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
						}
					}
					
					// OTHER BOUNDARY CONDITIONS
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel1[icv0],    vel_mgLevel1[icv0],    grad_u_mgLevel1[icv0], enthalpy_mgLevel1[icv0], grad_enthalpy_mgLevel1[icv0], temp_mgLevel1[icv0], RoM_mgLevel1[icv0],    gamma_mgLevel1[icv0],  kine0,
														 rho_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa], grad_u_mgLevel1[icv0], h_bfa_mgLevel1[ifa],     grad_enthalpy_mgLevel1[icv0], T_bfa_mgLevel1[ifa], RoM_bfa_mgLevel1[ifa], gam_bfa_mgLevel1[ifa], kine1,
														 mul_fa_mgLevel1[ifa], mut_fa_mgLevel1[ifa], lamOcp_fa_mgLevel1[ifa], kine_fa, vel_bfa_mgLevel1[ifa],
														 area, nVec);
							
							if (flagImplicit) {
								int noc00 = nbocv_i_mgLevel1[icv0]; // icv0's diagonal
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
						}
					}
				}
			}
		}
	}
	if (iMesh == 2) {
		// compute gradients, with boundary values
		calcCvVectorGradGreenGauss_mg(grad_u_mgLevel2, vel_mgLevel2, vel_bfa_mgLevel2, iMesh);
		calcCvScalarGradGreenGauss_mg(grad_enthalpy_mgLevel2, enthalpy_mgLevel2, h_bfa_mgLevel2, iMesh);
		
		// cycle through internal faces, assembling flux and matrix
		for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++)
		{
			int icv0 = cvofa_mgLevel2[ifa][0];
			int icv1 = cvofa_mgLevel2[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
			
			double uAux_fa[3] = { 0.5*(vel_mgLevel2[icv0][0]+ vel_mgLevel2[icv1][0]),
				0.5*(vel_mgLevel2[icv0][1]+ vel_mgLevel2[icv1][1]),
				0.5*(vel_mgLevel2[icv0][2]+ vel_mgLevel2[icv1][2])};
			
			// kinetic energy, if defined
			double kine0 = 0.0;
			double kine1 = 0.0;
			double kine_fa = 0.0;
			
			// calculate viscous flux
			addViscFlux_mg(Frhou, FrhoE, A0, A1,
										 rho_mgLevel2[icv0], vel_mgLevel2[icv0], grad_u_mgLevel2[icv0], enthalpy_mgLevel2[icv0], grad_enthalpy_mgLevel2[icv0], temp_mgLevel2[icv0], RoM_mgLevel2[icv0], gamma_mgLevel2[icv0], kine0,
										 rho_mgLevel2[icv1], vel_mgLevel2[icv1], grad_u_mgLevel2[icv1], enthalpy_mgLevel2[icv1], grad_enthalpy_mgLevel2[icv1], temp_mgLevel2[icv1], RoM_mgLevel2[icv1], gamma_mgLevel2[icv1], kine1,
										 mul_fa_mgLevel2[ifa], mut_fa_mgLevel2[ifa], lamOcp_fa_mgLevel2[ifa], kine_fa, uAux_fa, area, nVec);
			
			if (flagImplicit) {
				for (int i=1; i<5; i++)
					for (int j=0; j<5; j++) {
						A[noc00][i][j] -= A0[i][j];
						A[noc01][i][j] -= A1[i][j];
					}
				
				if (icv1 < ncv_mgLevel2) {
					for (int i=1; i<5; i++)
						for (int j=0; j<5; j++) {
							A[noc11][i][j] += A1[i][j];
							A[noc10][i][j] += A0[i][j];
						}
				}
			}
			
			// icv0 is always valid...
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel2) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
		}
		
		
		// cycle through boundary faces, assembling flux and matrix
		for (list<FaZone>::iterator zone = faZoneList_mgLevel2.begin(); zone != faZoneList_mgLevel2.end(); zone++)
		{
			if (zone->getKind() == FA_ZONE_BOUNDARY)
			{
				Param *param;
				if (getParam(param, zone->getName()))
				{
					// SYMMETRY BOUNDARY CONDITION
					if (param->getString() == "SYMMETRY") {
						// No viscous flux in this case
					}
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    grad_u_mgLevel2[icv0], enthalpy_mgLevel2[icv0], grad_enthalpy_mgLevel2[icv0], temp_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  kine0,
														 rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], grad_u_mgLevel2[icv0], h_bfa_mgLevel2[ifa],     grad_enthalpy_mgLevel2[icv0], T_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], kine1,
														 mul_fa_mgLevel2[ifa], 0.0, lamOcp_fa_mgLevel2[ifa], kine_fa, vel_bfa_mgLevel2[ifa],
														 area, nVec); 
							
							if (flagImplicit) {
								int noc00 = nbocv_i_mgLevel2[icv0]; // icv0's diagonal
								
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
						}
					}
					
					// OTHER BOUNDARY CONDITIONS
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel2[icv0],    vel_mgLevel2[icv0],    grad_u_mgLevel2[icv0], enthalpy_mgLevel2[icv0], grad_enthalpy_mgLevel2[icv0], temp_mgLevel2[icv0], RoM_mgLevel2[icv0],    gamma_mgLevel2[icv0],  kine0,
														 rho_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa], grad_u_mgLevel2[icv0], h_bfa_mgLevel2[ifa],     grad_enthalpy_mgLevel2[icv0], T_bfa_mgLevel2[ifa], RoM_bfa_mgLevel2[ifa], gam_bfa_mgLevel2[ifa], kine1,
														 mul_fa_mgLevel2[ifa], mut_fa_mgLevel2[ifa], lamOcp_fa_mgLevel2[ifa], kine_fa, vel_bfa_mgLevel2[ifa],
														 area, nVec);
							
							if (flagImplicit)
							{
								int noc00 = nbocv_i_mgLevel2[icv0]; // icv0's diagonal
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
						}
					}
				}
			}
		}
	}
	
	if (iMesh == 3) {
		// compute gradients, with boundary values
		calcCvVectorGradGreenGauss_mg(grad_u_mgLevel3, vel_mgLevel3, vel_bfa_mgLevel3, iMesh);
		calcCvScalarGradGreenGauss_mg(grad_enthalpy_mgLevel3, enthalpy_mgLevel3, h_bfa_mgLevel3, iMesh);
		
		// cycle through internal faces, assembling flux and matrix
		for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++)
		{
			int icv0 = cvofa_mgLevel3[ifa][0];
			int icv1 = cvofa_mgLevel3[ifa][1];
			assert( icv0 >= 0 );
			assert( icv1 >= 0 );
			
			int noc00, noc01, noc11, noc10;
			if (flagImplicit)
				getImplDependencyIndex_mg(noc00, noc01, noc11, noc10, icv0, icv1, iMesh);
			
			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
			
			double uAux_fa[3] = { 0.5*(vel_mgLevel3[icv0][0]+ vel_mgLevel3[icv1][0]),
				0.5*(vel_mgLevel3[icv0][1]+ vel_mgLevel3[icv1][1]),
				0.5*(vel_mgLevel3[icv0][2]+ vel_mgLevel3[icv1][2])};
			
			// kinetic energy, if defined
			double kine0 = 0.0;
			double kine1 = 0.0;
			double kine_fa = 0.0;
			
			// calculate viscous flux
			addViscFlux_mg(Frhou, FrhoE, A0, A1,
										 rho_mgLevel3[icv0], vel_mgLevel3[icv0], grad_u_mgLevel3[icv0], enthalpy_mgLevel3[icv0], grad_enthalpy_mgLevel3[icv0], temp_mgLevel3[icv0], RoM_mgLevel3[icv0], gamma_mgLevel3[icv0], kine0,
										 rho_mgLevel3[icv1], vel_mgLevel3[icv1], grad_u_mgLevel3[icv1], enthalpy_mgLevel3[icv1], grad_enthalpy_mgLevel3[icv1], temp_mgLevel3[icv1], RoM_mgLevel3[icv1], gamma_mgLevel3[icv1], kine1,
										 mul_fa_mgLevel3[ifa], mut_fa_mgLevel3[ifa], lamOcp_fa_mgLevel3[ifa], kine_fa, uAux_fa, area, nVec);
			
			if (flagImplicit) {
				for (int i=1; i<5; i++)
					for (int j=0; j<5; j++) {
						A[noc00][i][j] -= A0[i][j];
						A[noc01][i][j] -= A1[i][j];
					}
				
				if (icv1 < ncv_mgLevel3) {
					for (int i=1; i<5; i++)
						for (int j=0; j<5; j++) {
							A[noc11][i][j] += A1[i][j];
							A[noc10][i][j] += A0[i][j];
						}
				}
			}
			
			// icv0 is always valid...
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv0][i] -= Frhou[i];
			rhs_rhoE[icv0] -= FrhoE;
			
			// icv1 can be ghost...
			if (icv1 < ncv_mgLevel3) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[icv1][i] += Frhou[i];
				rhs_rhoE[icv1] += FrhoE;
			}
		}
		
		
		// cycle through boundary faces, assembling flux and matrix
		for (list<FaZone>::iterator zone = faZoneList_mgLevel3.begin(); zone != faZoneList_mgLevel3.end(); zone++)
		{
			if (zone->getKind() == FA_ZONE_BOUNDARY)
			{
				Param *param;
				if (getParam(param, zone->getName()))
				{
					// SYMMETRY BOUNDARY CONDITION
					if (param->getString() == "SYMMETRY") {
						// No viscous flux in this case
					}
					// WALL BOUNDARY CONDITION
					else if (param->getString() == "WALL") {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    grad_u_mgLevel3[icv0], enthalpy_mgLevel3[icv0], grad_enthalpy_mgLevel3[icv0], temp_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  kine0,
														 rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], grad_u_mgLevel3[icv0], h_bfa_mgLevel3[ifa],     grad_enthalpy_mgLevel3[icv0], T_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], kine1,
														 mul_fa_mgLevel3[ifa], 0.0, lamOcp_fa_mgLevel3[ifa], kine_fa, vel_bfa_mgLevel3[ifa],
														 area, nVec); 
							
							if (flagImplicit) {
								int noc00 = nbocv_i_mgLevel3[icv0]; // icv0's diagonal
								
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
							
						}
					}
					
					// OTHER BOUNDARY CONDITIONS
					else {
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							double kine0 = 0.0;
							double kine1 = 0.0;
							double kine_fa = 0.0;
							
							// calculate viscous flux
							addViscFlux_mg(Frhou, FrhoE, A0, NULL,
														 rho_mgLevel3[icv0],    vel_mgLevel3[icv0],    grad_u_mgLevel3[icv0], enthalpy_mgLevel3[icv0], grad_enthalpy_mgLevel3[icv0], temp_mgLevel3[icv0], RoM_mgLevel3[icv0],    gamma_mgLevel3[icv0],  kine0,
														 rho_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa], grad_u_mgLevel3[icv0], h_bfa_mgLevel3[ifa],     grad_enthalpy_mgLevel3[icv0], T_bfa_mgLevel3[ifa], RoM_bfa_mgLevel3[ifa], gam_bfa_mgLevel3[ifa], kine1,
														 mul_fa_mgLevel3[ifa], mut_fa_mgLevel3[ifa], lamOcp_fa_mgLevel3[ifa], kine_fa, vel_bfa_mgLevel3[ifa],
														 area, nVec);
							
							if (flagImplicit)
							{
								int noc00 = nbocv_i_mgLevel3[icv0]; // icv0's diagonal
								for (int i=1; i<5; i++)
									for (int j=0; j<5; j++)
										A[noc00][i][j] -= A0[i][j];
							}
							
							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= Frhou[i];
							rhs_rhoE[icv0] -= FrhoE;
						}
					}
				}
			}
		}
	}
	
  if (A0  != NULL)  delete [] A0;
  if (A1  != NULL)  delete [] A1;
}

void JoeWithModels::setBC_mg(int iMesh) {
	
	if (iMesh == 1) {
		static int first = 1;
		int bc_err = 0;
		
		for (list<FaZone>::iterator zone = faZoneList_mgLevel1.begin(); zone != faZoneList_mgLevel1.end(); zone++)
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName())) {
					
					// HOOK BOUNDARY CONDITION
					if (param->getString() == "HOOK") {
						if ((first) && (mpi_rank == 0))
							cout << "Applying HOOK                to zone: "<< zone->getName() << endl;
						
						boundaryHook_mg(T_bfa_mgLevel1, vel_bfa_mgLevel1, p_bfa_mgLevel1, &(*zone), iMesh);
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					
					// CBC BOUNDARY CONDITION
					else if (param->getString() == "CBC") {
						double u_bc[3], T_bc, p_bc;
						
						for (int i=0; i<3; i++)
							u_bc[i] = param->getDouble(i+2);
						T_bc = param->getDouble(5);
						p_bc = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
							{
								double velMagN = vecDotVec3d(vel_mgLevel1[icv0], nVec);
								double mach = fabs(velMagN)/sos_mgLevel1[icv0];
								
								if (mach >= 1.0) {
									T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel1[ifa][i] = vel_mgLevel1[icv0][i];
									p_bfa_mgLevel1[ifa] = press_mgLevel1[icv0];
								}
								else
								{
									T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel1[ifa][i] = vel_mgLevel1[icv0][i];
									p_bfa_mgLevel1[ifa] = p_bc;
								}
							}
							else      // inlet
							{
								T_bfa_mgLevel1[ifa] = T_bc;
								for (int i=0; i<3; i++)
									vel_bfa_mgLevel1[ifa][i] = u_bc[i];
								p_bfa_mgLevel1[ifa] = p_bc;
							}
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// CBC SUBSONIC INLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_INLET")
					{
						double angleU[3], Ttot, htot, ptot;
						
						for (int i=0; i<3; i++)
							angleU[i] = param->getDouble(i+2);
						Ttot = param->getDouble(5);
						ptot = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							double u[3] = {rhou_mgLevel1[icv0][0]/rho_mgLevel1[icv0], rhou_mgLevel1[icv0][1]/rho_mgLevel1[icv0], rhou_mgLevel1[icv0][2]/rho_mgLevel1[icv0]};
							double wPow2 = vecDotVec3d(u, u);               // velocity squared
							double vel = sqrt(wPow2);                       // apply angle to extrapolated velocity
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel1[ifa][i] = angleU[i]*vel;
							T_bfa_mgLevel1[ifa] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);                  // total enthalpy from total temperature
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							double wPow2 = vecDotVec3d(vel_bfa_mgLevel1[ifa], vel_bfa_mgLevel1[ifa]);         // velocity squared          
							h_bfa_mgLevel1[ifa] -= 0.5* wPow2;                                       // static enthalpy
						}
						
						ComputeBCProperties_H_mg(&(*zone), iMesh);                  // static temperature and thermo properties from static enthalpy
						
						// Assumes isentropic relations to determine static pressure (= constant cp)
						// At first approximation ok, but could be improved; should for now be considered in defining p_bc
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
							p_bfa_mgLevel1[ifa] = ptot * pow(T_bfa_mgLevel1[ifa]/Ttot, gam_bfa_mgLevel1[ifa]/(gam_bfa_mgLevel1[ifa]-1.0));
					}
					// .............................................................................................
					// CBC SUBSONIC OUTLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_OUTLET")  {
						double p_bc = param->getDouble(2);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
							p_bfa_mgLevel1[ifa] = p_bc;
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel1[ifa][i] = rhou_mgLevel1[icv0][i]/rho_mgLevel1[icv0];
							
							T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// SYMMETRY BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "SYMMETRY")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
							
							// flip u, APPROXIMATION ---> take velocity at the cell center
							double u0[3] = {vel_mgLevel1[icv0][0], vel_mgLevel1[icv0][1], vel_mgLevel1[icv0][2]};
							double un = vecDotVec3d(nVec, u0);
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel1[ifa][i] = u0[i] - 1.0*un*nVec[i];
							//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);
							
							T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];
							p_bfa_mgLevel1[ifa] = press_mgLevel1[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// NEUMANN BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "NEUMANN")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel1[ifa][i] = rhou_mgLevel1[icv0][i]/rho_mgLevel1[icv0];
							
							T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];
							p_bfa_mgLevel1[ifa] = press_mgLevel1[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// WALL BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "WALL")
					{
						int i=0;
						double T_bc = 0.0;
						if ((i = param->findString("TEMP")) != 0)
							T_bc = param->getDouble(i+1);
						
						if ((first)&&(mpi_rank == 0))
						{
							if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
							else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
						}
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel1[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel1[ifa][i] = 0.0;
							
							if (T_bc > 0.0)   T_bfa_mgLevel1[ifa] = T_bc;           // wall temperature
							else              T_bfa_mgLevel1[ifa] = temp_mgLevel1[icv0];     // adiabatic wall
							
							p_bfa_mgLevel1[ifa] = press_mgLevel1[icv0];                      // Assumes zero pressure gradient at the wall
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// OTHER BOUNDARY CONDITIONS
					// .............................................................................................
					else {
						if (mpi_rank == 0)
							cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
						bc_err = 1;
					}
				}
				else {
					if (mpi_rank == 0)
						cerr << "Error: no bc set for: "<< zone->getName() << endl;
					bc_err = 1;
				}
			}
		
		// update density at boundary using EOS
		for (int ifa = 0; ifa < nfa_b_mgLevel1; ifa++)
			rho_bfa_mgLevel1[ifa] = p_bfa_mgLevel1[ifa] / (RoM_bfa_mgLevel1[ifa] * T_bfa_mgLevel1[ifa]);
		
		
		if (bc_err != 0)
			throw(-1);
		
		first = 0;
	}
	
	if (iMesh == 2) {
		
	  static int first = 1;
		int bc_err = 0;
		
		for (list<FaZone>::iterator zone = faZoneList_mgLevel2.begin(); zone != faZoneList_mgLevel2.end(); zone++)
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName())) {
					
					// HOOK BOUNDARY CONDITION
					if (param->getString() == "HOOK") {
						if ((first) && (mpi_rank == 0))
							cout << "Applying HOOK                to zone: "<< zone->getName() << endl;
						
						boundaryHook_mg(T_bfa_mgLevel2, vel_bfa_mgLevel2, p_bfa_mgLevel2, &(*zone), iMesh);
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					
					// CBC BOUNDARY CONDITION
					else if (param->getString() == "CBC") {
						double u_bc[3], T_bc, p_bc;
						
						for (int i=0; i<3; i++)
							u_bc[i] = param->getDouble(i+2);
						T_bc = param->getDouble(5);
						p_bc = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
							{
								double velMagN = vecDotVec3d(vel_mgLevel2[icv0], nVec);
								double mach = fabs(velMagN)/sos_mgLevel2[icv0];
								
								if (mach >= 1.0) {
									T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel2[ifa][i] = vel_mgLevel2[icv0][i];
									p_bfa_mgLevel2[ifa] = press_mgLevel2[icv0];
								}
								else
								{
									T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel2[ifa][i] = vel_mgLevel2[icv0][i];
									p_bfa_mgLevel2[ifa] = p_bc;
								}
							}
							else      // inlet
							{
								T_bfa_mgLevel2[ifa] = T_bc;
								for (int i=0; i<3; i++)
									vel_bfa_mgLevel2[ifa][i] = u_bc[i];
								p_bfa_mgLevel2[ifa] = p_bc;
							}
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// CBC SUBSONIC INLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_INLET")
					{
						double angleU[3], Ttot, htot, ptot;
						
						for (int i=0; i<3; i++)
							angleU[i] = param->getDouble(i+2);
						Ttot = param->getDouble(5);
						ptot = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							double u[3] = {rhou_mgLevel2[icv0][0]/rho_mgLevel2[icv0], rhou_mgLevel2[icv0][1]/rho_mgLevel2[icv0], rhou_mgLevel2[icv0][2]/rho_mgLevel2[icv0]};
							double wPow2 = vecDotVec3d(u, u);               // velocity squared
							double vel = sqrt(wPow2);                       // apply angle to extrapolated velocity
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel2[ifa][i] = angleU[i]*vel;
							T_bfa_mgLevel2[ifa] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);                  // total enthalpy from total temperature
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							double wPow2 = vecDotVec3d(vel_bfa_mgLevel2[ifa], vel_bfa_mgLevel2[ifa]);         // velocity squared          
							h_bfa_mgLevel2[ifa] -= 0.5* wPow2;                                       // static enthalpy
						}
						
						ComputeBCProperties_H_mg(&(*zone), iMesh);                  // static temperature and thermo properties from static enthalpy
						
						// Assumes isentropic relations to determine static pressure (= constant cp)
						// At first approximation ok, but could be improved; should for now be considered in defining p_bc
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
							p_bfa_mgLevel2[ifa] = ptot * pow(T_bfa_mgLevel2[ifa]/Ttot, gam_bfa_mgLevel2[ifa]/(gam_bfa_mgLevel2[ifa]-1.0));
					}
					// .............................................................................................
					// CBC SUBSONIC OUTLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_OUTLET")  {
						double p_bc = param->getDouble(2);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
							p_bfa_mgLevel2[ifa] = p_bc;
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel2[ifa][i] = rhou_mgLevel2[icv0][i]/rho_mgLevel2[icv0];
							
							T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// SYMMETRY BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "SYMMETRY")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
							
							// flip u, APPROXIMATION ---> take velocity at the cell center
							double u0[3] = {vel_mgLevel2[icv0][0], vel_mgLevel2[icv0][1], vel_mgLevel2[icv0][2]};
							double un = vecDotVec3d(nVec, u0);
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel2[ifa][i] = u0[i] - 1.0*un*nVec[i];
							//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);
							
							T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];
							p_bfa_mgLevel2[ifa] = press_mgLevel2[icv0];
							
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// NEUMANN BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "NEUMANN")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel2[ifa][i] = rhou_mgLevel2[icv0][i]/rho_mgLevel2[icv0];
							
							T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];
							p_bfa_mgLevel2[ifa] = press_mgLevel2[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// WALL BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "WALL")
					{
						int i=0;
						double T_bc = 0.0;
						if ((i = param->findString("TEMP")) != 0)
							T_bc = param->getDouble(i+1);
						
						if ((first)&&(mpi_rank == 0))
						{
							if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
							else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
						}
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel2[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel2[ifa][i] = 0.0;
							
							if (T_bc > 0.0)   T_bfa_mgLevel2[ifa] = T_bc;           // wall temperature
							else              T_bfa_mgLevel2[ifa] = temp_mgLevel2[icv0];     // adiabatic wall
							
							p_bfa_mgLevel2[ifa] = press_mgLevel2[icv0];                      // Assumes zero pressure gradient at the wall
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// OTHER BOUNDARY CONDITIONS
					// .............................................................................................
					else {
						if (mpi_rank == 0)
							cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
						bc_err = 1;
					}
				}
				else {
					if (mpi_rank == 0)
						cerr << "Error: no bc set for: "<< zone->getName() << endl;
					bc_err = 1;
				}
			}
		
		// update density at boundary using EOS
		for (int ifa = 0; ifa < nfa_b_mgLevel2; ifa++)
			rho_bfa_mgLevel2[ifa] = p_bfa_mgLevel2[ifa] / (RoM_bfa_mgLevel2[ifa] * T_bfa_mgLevel2[ifa]);
		
		if (bc_err != 0)
			throw(-1);
		
		first = 0;
		
		
	}
	
	if (iMesh == 3) {
	  static int first = 1;
		int bc_err = 0;
		
		for (list<FaZone>::iterator zone = faZoneList_mgLevel3.begin(); zone != faZoneList_mgLevel3.end(); zone++)
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;
				
				if (getParam(param, zone->getName())) {
					
					// HOOK BOUNDARY CONDITION
					if (param->getString() == "HOOK") {
						if ((first) && (mpi_rank == 0))
							cout << "Applying HOOK                to zone: "<< zone->getName() << endl;
						
						boundaryHook_mg(T_bfa_mgLevel3, vel_bfa_mgLevel3, p_bfa_mgLevel3, &(*zone), iMesh);
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					
					// CBC BOUNDARY CONDITION
					else if (param->getString() == "CBC") {
						double u_bc[3], T_bc, p_bc;
						
						for (int i=0; i<3; i++)
							u_bc[i] = param->getDouble(i+2);
						T_bc = param->getDouble(5);
						p_bc = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
							{
								double velMagN = vecDotVec3d(vel_mgLevel3[icv0], nVec);
								double mach = fabs(velMagN)/sos_mgLevel3[icv0];
								
								if (mach >= 1.0) {
									T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel3[ifa][i] = vel_mgLevel3[icv0][i];
									p_bfa_mgLevel3[ifa] = press_mgLevel3[icv0];
								}
								else
								{
									T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];
									for (int i=0; i<3; i++)
										vel_bfa_mgLevel3[ifa][i] = vel_mgLevel3[icv0][i];
									p_bfa_mgLevel3[ifa] = p_bc;
								}
							}
							else      // inlet
							{
								T_bfa_mgLevel3[ifa] = T_bc;
								for (int i=0; i<3; i++)
									vel_bfa_mgLevel3[ifa][i] = u_bc[i];
								p_bfa_mgLevel3[ifa] = p_bc;
							}
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// CBC SUBSONIC INLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_INLET")
					{
						double angleU[3], Ttot, htot, ptot;
						
						for (int i=0; i<3; i++)
							angleU[i] = param->getDouble(i+2);
						Ttot = param->getDouble(5);
						ptot = param->getDouble(6);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							double u[3] = {rhou_mgLevel3[icv0][0]/rho_mgLevel3[icv0], rhou_mgLevel3[icv0][1]/rho_mgLevel3[icv0], rhou_mgLevel3[icv0][2]/rho_mgLevel3[icv0]};
							double wPow2 = vecDotVec3d(u, u);               // velocity squared
							double vel = sqrt(wPow2);                       // apply angle to extrapolated velocity
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel3[ifa][i] = angleU[i]*vel;
							T_bfa_mgLevel3[ifa] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);                  // total enthalpy from total temperature
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							double wPow2 = vecDotVec3d(vel_bfa_mgLevel3[ifa], vel_bfa_mgLevel3[ifa]);         // velocity squared          
							h_bfa_mgLevel3[ifa] -= 0.5* wPow2;                                       // static enthalpy
						}
						
						ComputeBCProperties_H_mg(&(*zone), iMesh);                  // static temperature and thermo properties from static enthalpy
						
						// Assumes isentropic relations to determine static pressure (= constant cp)
						// At first approximation ok, but could be improved; should for now be considered in defining p_bc
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
							p_bfa_mgLevel3[ifa] = ptot * pow(T_bfa_mgLevel3[ifa]/Ttot, gam_bfa_mgLevel3[ifa]/(gam_bfa_mgLevel3[ifa]-1.0));
					}
					// .............................................................................................
					// CBC SUBSONIC OUTLET BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "CBC_SUBSONIC_OUTLET")  {
						double p_bc = param->getDouble(2);
						
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
							p_bfa_mgLevel3[ifa] = p_bc;
							for (int i=0; i<3; i++)
								vel_bfa_mgLevel3[ifa][i] = rhou_mgLevel3[icv0][i]/rho_mgLevel3[icv0];
							
							T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// SYMMETRY BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "SYMMETRY")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							double nVec[3];
							double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
							
							// flip u, APPROXIMATION ---> take velocity at the cell center
							double u0[3] = {vel_mgLevel3[icv0][0], vel_mgLevel3[icv0][1], vel_mgLevel3[icv0][2]};
							double un = vecDotVec3d(nVec, u0);
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel3[ifa][i] = u0[i] - 1.0*un*nVec[i];
							//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);
							
							T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];
							p_bfa_mgLevel3[ifa] = press_mgLevel3[icv0];
							
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// NEUMANN BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "NEUMANN")
					{
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel3[ifa][i] = rhou_mgLevel3[icv0][i]/rho_mgLevel3[icv0];
							
							T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];
							p_bfa_mgLevel3[ifa] = press_mgLevel3[icv0];
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// WALL BOUNDARY CONDITION
					// .............................................................................................
					else if (param->getString() == "WALL")
					{
						int i=0;
						double T_bc = 0.0;
						if ((i = param->findString("TEMP")) != 0)
							T_bc = param->getDouble(i+1);
						
						if ((first)&&(mpi_rank == 0))
						{
							if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
							else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
						}
						
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa_mgLevel3[ifa][0];
							assert( icv0 >= 0 );
							
							for (int i = 0; i < 3; i++)
								vel_bfa_mgLevel3[ifa][i] = 0.0;
							
							if (T_bc > 0.0)   T_bfa_mgLevel3[ifa] = T_bc;           // wall temperature
							else              T_bfa_mgLevel3[ifa] = temp_mgLevel3[icv0];     // adiabatic wall
							
							p_bfa_mgLevel3[ifa] = press_mgLevel3[icv0];                      // Assumes zero pressure gradient at the wall
						}
						
						ComputeBCProperties_T_mg(&(*zone), iMesh);
					}
					// .............................................................................................
					// OTHER BOUNDARY CONDITIONS
					// .............................................................................................
					else {
						if (mpi_rank == 0)
							cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
						bc_err = 1;
					}
				}
				else {
					if (mpi_rank == 0)
						cerr << "Error: no bc set for: "<< zone->getName() << endl;
					bc_err = 1;
				}
			}
		
		// update density at boundary using EOS
		for (int ifa = 0; ifa < nfa_b_mgLevel3; ifa++)
			rho_bfa_mgLevel3[ifa] = p_bfa_mgLevel3[ifa] / (RoM_bfa_mgLevel3[ifa] * T_bfa_mgLevel3[ifa]);
		
		if (bc_err != 0)
			throw(-1);
		
		first = 0;
		
		
	}
	
}


void JoeWithModels::Set_SaveSolution(int index, int iMesh) {
	unsigned long icv;
	if (iMesh == 1) {
		if (index > 0)
			for(icv = 0; icv < ncv_mgLevel1; icv++) {
				rho_old_mgLevel1[icv] = rho_mgLevel1[icv];
				rhou_old_mgLevel1[icv][0] = rhou_mgLevel1[icv][0];
				rhou_old_mgLevel1[icv][1] = rhou_mgLevel1[icv][1];
				rhou_old_mgLevel1[icv][2] = rhou_mgLevel1[icv][2];
				rhoE_old_mgLevel1[icv] = rhoE_mgLevel1[icv];
			}
		if (index < 0)
			for(icv = 0; icv < ncv_mgLevel1; icv++) {
				rho_mgLevel1[icv] = rho_old_mgLevel1[icv];
				rhou_mgLevel1[icv][0] = rhou_old_mgLevel1[icv][0];
				rhou_mgLevel1[icv][1] = rhou_old_mgLevel1[icv][1];
				rhou_mgLevel1[icv][2] = rhou_old_mgLevel1[icv][2];
				rhoE_mgLevel1[icv] = rhoE_old_mgLevel1[icv];
			}
	}
	if (iMesh == 2) {
		if (index > 0)
			for(icv = 0; icv < ncv_mgLevel2; icv++) {
				rho_old_mgLevel2[icv] = rho_mgLevel2[icv];
				rhou_old_mgLevel2[icv][0] = rhou_mgLevel2[icv][0];
				rhou_old_mgLevel2[icv][1] = rhou_mgLevel2[icv][1];
				rhou_old_mgLevel2[icv][2] = rhou_mgLevel2[icv][2];
				rhoE_old_mgLevel2[icv] = rhoE_mgLevel2[icv];
			}
		if (index < 0)
			for(icv = 0; icv < ncv_mgLevel2; icv++) {
				rho_mgLevel2[icv] = rho_old_mgLevel2[icv];
				rhou_mgLevel2[icv][0] = rhou_old_mgLevel2[icv][0];
				rhou_mgLevel2[icv][1] = rhou_old_mgLevel2[icv][1];
				rhou_mgLevel2[icv][2] = rhou_old_mgLevel2[icv][2];
				rhoE_mgLevel2[icv] = rhoE_old_mgLevel2[icv];
			}
	}
	if (iMesh == 3) {
		if (index > 0)
			for(icv = 0; icv < ncv_mgLevel3; icv++) {
				rho_old_mgLevel3[icv] = rho_mgLevel3[icv];
				rhou_old_mgLevel3[icv][0] = rhou_mgLevel3[icv][0];
				rhou_old_mgLevel3[icv][1] = rhou_mgLevel3[icv][1];
				rhou_old_mgLevel3[icv][2] = rhou_mgLevel3[icv][2];
				rhoE_old_mgLevel3[icv] = rhoE_mgLevel3[icv];
			}
		if (index < 0)
			for(icv = 0; icv < ncv_mgLevel3; icv++) {
				rho_mgLevel3[icv] = rho_old_mgLevel3[icv];
				rhou_mgLevel3[icv][0] = rhou_old_mgLevel3[icv][0];
				rhou_mgLevel3[icv][1] = rhou_old_mgLevel3[icv][1];
				rhou_mgLevel3[icv][2] = rhou_old_mgLevel3[icv][2];
				rhoE_mgLevel3[icv] = rhoE_old_mgLevel3[icv];
			}
	}
}

void JoeWithModels::SetProlongated_Correction(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	double icv_volume_parent, icv_volume_children;
	
	if (iMesh == 1) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel1[icv_coarse];
			
			rho_old_mgLevel1[icv_coarse] = rho_mgLevel1[icv_coarse];
			rhou_old_mgLevel1[icv_coarse][0] = rhou_mgLevel1[icv_coarse][0];
			rhou_old_mgLevel1[icv_coarse][1] = rhou_mgLevel1[icv_coarse][1];
			rhou_old_mgLevel1[icv_coarse][2] = rhou_mgLevel1[icv_coarse][2];
			rhoE_old_mgLevel1[icv_coarse] = rhoE_mgLevel1[icv_coarse];
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				icv_volume_children = cv_volume[icv_fine];
				
				rho_old_mgLevel1[icv_coarse] -= rho[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel1[icv_coarse][0] -= rhou[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel1[icv_coarse][1] -= rhou[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel1[icv_coarse][2] -= rhou[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_old_mgLevel1[icv_coarse] -= rhoE[icv_fine]*icv_volume_children/icv_volume_parent;
			}
			
		}
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++)
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				
				rhs_rho[icv_fine] = rho_old_mgLevel1[icv_coarse];
				rhs_rhou[icv_fine][0] = rhou_old_mgLevel1[icv_coarse][0];
				rhs_rhou[icv_fine][1] = rhou_old_mgLevel1[icv_coarse][1];
				rhs_rhou[icv_fine][2] = rhou_old_mgLevel1[icv_coarse][2];
				rhs_rhoE[icv_fine] = rhoE_old_mgLevel1[icv_coarse];
			}
	}
	
	if (iMesh == 2) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel2[icv_coarse];
			
			rho_old_mgLevel2[icv_coarse] = rho_mgLevel2[icv_coarse];
			rhou_old_mgLevel2[icv_coarse][0] = rhou_mgLevel2[icv_coarse][0];
			rhou_old_mgLevel2[icv_coarse][1] = rhou_mgLevel2[icv_coarse][1];
			rhou_old_mgLevel2[icv_coarse][2] = rhou_mgLevel2[icv_coarse][2];
			rhoE_old_mgLevel2[icv_coarse] = rhoE_mgLevel2[icv_coarse];
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				icv_volume_children = cv_volume_mgLevel1[icv_fine];
				
				rho_old_mgLevel2[icv_coarse] -= rho_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel2[icv_coarse][0] -= rhou_mgLevel1[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel2[icv_coarse][1] -= rhou_mgLevel1[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel2[icv_coarse][2] -= rhou_mgLevel1[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_old_mgLevel2[icv_coarse] -= rhoE_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
			}
			
		}
		
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++)
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				
				rhs_rho[icv_fine] = rho_old_mgLevel2[icv_coarse];
				rhs_rhou[icv_fine][0] = rhou_old_mgLevel2[icv_coarse][0];
				rhs_rhou[icv_fine][1] = rhou_old_mgLevel2[icv_coarse][1];
				rhs_rhou[icv_fine][2] = rhou_old_mgLevel2[icv_coarse][2];
				rhs_rhoE[icv_fine] = rhoE_old_mgLevel2[icv_coarse];
			}
	}
	
	if (iMesh == 3) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel3[icv_coarse];
			
			rho_old_mgLevel3[icv_coarse] = rho_mgLevel3[icv_coarse];
			rhou_old_mgLevel3[icv_coarse][0] = rhou_mgLevel3[icv_coarse][0];
			rhou_old_mgLevel3[icv_coarse][1] = rhou_mgLevel3[icv_coarse][1];
			rhou_old_mgLevel3[icv_coarse][2] = rhou_mgLevel3[icv_coarse][2];
			rhoE_old_mgLevel3[icv_coarse] = rhoE_mgLevel3[icv_coarse];
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				icv_volume_children = cv_volume_mgLevel2[icv_fine];
				
				rho_old_mgLevel3[icv_coarse] -= rho_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel3[icv_coarse][0] -= rhou_mgLevel2[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel3[icv_coarse][1] -= rhou_mgLevel2[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_old_mgLevel3[icv_coarse][2] -= rhou_mgLevel2[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_old_mgLevel3[icv_coarse] -= rhoE_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
			}
			
		}
		
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++)
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				
				rhs_rho[icv_fine] = rho_old_mgLevel3[icv_coarse];
				rhs_rhou[icv_fine][0] = rhou_old_mgLevel3[icv_coarse][0];
				rhs_rhou[icv_fine][1] = rhou_old_mgLevel3[icv_coarse][1];
				rhs_rhou[icv_fine][2] = rhou_old_mgLevel3[icv_coarse][2];
				rhs_rhoE[icv_fine] = rhoE_old_mgLevel3[icv_coarse];
			}
	}
}


void JoeWithModels::SetResidual_Smoothing(double *rhs_rho, double (*rhs_rhou)[3], 
										  double *rhs_rhoE, int nSmooth, double val_smooth_coeff, 
										  int iMesh) {
	
	if (iMesh == 0) {
		
		for (int icv = 0; icv < ncv; icv++) {
			rho_res_old[icv] = rhs_rho[icv];
			rhou_res_old[icv][0] = rhs_rhou[icv][0];
			rhou_res_old[icv][1] = rhs_rhou[icv][1];
			rhou_res_old[icv][2] = rhs_rhou[icv][2];
			rhoE_res_old[icv] = rhs_rhoE[icv];
		}
		
		for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
			for (int icv = 0; icv < ncv; icv++) {
				rho_res_sum[icv] = 0.0;
				rhou_res_sum[icv][0] = 0.0;
				rhou_res_sum[icv][1] = 0.0;
				rhou_res_sum[icv][2] = 0.0;
				rhoE_res_sum[icv] = 0.0;
			}
			
			for (int ifa = nfa_b; ifa < nfa; ifa++) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				if (icv1 < ncv) {
					rho_res_sum[icv0] += rhs_rho[icv1]; rho_res_sum[icv1] += rhs_rho[icv0];
					rhou_res_sum[icv0][0] += rhs_rhou[icv1][0]; rhou_res_sum[icv1][0] += rhs_rhou[icv0][0];
					rhou_res_sum[icv0][1] += rhs_rhou[icv1][1]; rhou_res_sum[icv1][1] += rhs_rhou[icv0][1];
					rhou_res_sum[icv0][2] += rhs_rhou[icv1][2]; rhou_res_sum[icv1][2] += rhs_rhou[icv0][2];
					rhoE_res_sum[icv0] += rhs_rhoE[icv1]; rhoE_res_sum[icv1] += rhs_rhoE[icv0];
				}
			}
			
			for (int icv = 0; icv < ncv; icv++) {
				int non_f = nbocv_i[icv];
				int non_l = nbocv_i[icv+1]-1;
				double nneigh = double(non_l-non_f);
				rhs_rho[icv] = (rho_res_old[icv] + val_smooth_coeff*rho_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][0] = (rhou_res_old[icv][0] + val_smooth_coeff*rhou_res_sum[icv][0]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][1] = (rhou_res_old[icv][1] + val_smooth_coeff*rhou_res_sum[icv][1]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][2] = (rhou_res_old[icv][2] + val_smooth_coeff*rhou_res_sum[icv][2]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhoE[icv] = (rhoE_res_old[icv] + val_smooth_coeff*rhoE_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				
			}
			
			for (int ifa = 0; ifa < nfa_b; ifa++) {
				int icv0 = cvofa[ifa][0];
				rhs_rho[icv0] = rho_res_old[icv0];
				rhs_rhou[icv0][0] = rhou_res_old[icv0][0];
				rhs_rhou[icv0][1] = rhou_res_old[icv0][1];
				rhs_rhou[icv0][2] = rhou_res_old[icv0][2];
				rhs_rhoE[icv0] = rhoE_res_old[icv0];
			}
		}
	}
	
	if (iMesh == 1) {
		
		for (int icv = 0; icv < ncv_mgLevel1; icv++) {
			rho_res_old[icv] = rhs_rho[icv];
			rhou_res_old[icv][0] = rhs_rhou[icv][0];
			rhou_res_old[icv][1] = rhs_rhou[icv][1];
			rhou_res_old[icv][2] = rhs_rhou[icv][2];
			rhoE_res_old[icv] = rhs_rhoE[icv];
		}
		
		for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
			for (int icv = 0; icv < ncv_mgLevel1; icv++) {
				rho_res_sum[icv] = 0.0;
				rhou_res_sum[icv][0] = 0.0;
				rhou_res_sum[icv][1] = 0.0;
				rhou_res_sum[icv][2] = 0.0;
				rhoE_res_sum[icv] = 0.0;
			}
			
			for (int ifa = nfa_b; ifa < nfa_mgLevel1; ifa++) {
				int icv0 = cvofa_mgLevel1[ifa][0];
				int icv1 = cvofa_mgLevel1[ifa][1];
				if (icv1 < ncv_mgLevel1) {
					rho_res_sum[icv0] += rhs_rho[icv1]; rho_res_sum[icv1] += rhs_rho[icv0];
					rhou_res_sum[icv0][0] += rhs_rhou[icv1][0]; rhou_res_sum[icv1][0] += rhs_rhou[icv0][0];
					rhou_res_sum[icv0][1] += rhs_rhou[icv1][1]; rhou_res_sum[icv1][1] += rhs_rhou[icv0][1];
					rhou_res_sum[icv0][2] += rhs_rhou[icv1][2]; rhou_res_sum[icv1][2] += rhs_rhou[icv0][2];
					rhoE_res_sum[icv0] += rhs_rhoE[icv1]; rhoE_res_sum[icv1] += rhs_rhoE[icv0];
				}
			}
			
			for (int icv = 0; icv < ncv_mgLevel1; icv++) {
				int non_f = nbocv_i[icv];
				int non_l = nbocv_i[icv+1]-1;
				double nneigh = double(non_l-non_f);
				rhs_rho[icv] = (rho_res_old[icv] + val_smooth_coeff*rho_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][0] = (rhou_res_old[icv][0] + val_smooth_coeff*rhou_res_sum[icv][0]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][1] = (rhou_res_old[icv][1] + val_smooth_coeff*rhou_res_sum[icv][1]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][2] = (rhou_res_old[icv][2] + val_smooth_coeff*rhou_res_sum[icv][2]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhoE[icv] = (rhoE_res_old[icv] + val_smooth_coeff*rhoE_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				
			}
			
			for (int ifa = 0; ifa < nfa_b_mgLevel1; ifa++) {
				int icv0 = cvofa[ifa][0];
				rhs_rho[icv0] = rho_res_old[icv0];
				rhs_rhou[icv0][0] = rhou_res_old[icv0][0];
				rhs_rhou[icv0][1] = rhou_res_old[icv0][1];
				rhs_rhou[icv0][2] = rhou_res_old[icv0][2];
				rhs_rhoE[icv0] = rhoE_res_old[icv0];
			}
		}
	}
	
	if (iMesh == 2) {
		
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			rho_res_old[icv] = rhs_rho[icv];
			rhou_res_old[icv][0] = rhs_rhou[icv][0];
			rhou_res_old[icv][1] = rhs_rhou[icv][1];
			rhou_res_old[icv][2] = rhs_rhou[icv][2];
			rhoE_res_old[icv] = rhs_rhoE[icv];
		}
		
		for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
			for (int icv = 0; icv < ncv_mgLevel2; icv++) {
				rho_res_sum[icv] = 0.0;
				rhou_res_sum[icv][0] = 0.0;
				rhou_res_sum[icv][1] = 0.0;
				rhou_res_sum[icv][2] = 0.0;
				rhoE_res_sum[icv] = 0.0;
			}
			
			for (int ifa = nfa_b; ifa < nfa_mgLevel2; ifa++) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				if (icv1 < ncv_mgLevel2) {
					rho_res_sum[icv0] += rhs_rho[icv1]; rho_res_sum[icv1] += rhs_rho[icv0];
					rhou_res_sum[icv0][0] += rhs_rhou[icv1][0]; rhou_res_sum[icv1][0] += rhs_rhou[icv0][0];
					rhou_res_sum[icv0][1] += rhs_rhou[icv1][1]; rhou_res_sum[icv1][1] += rhs_rhou[icv0][1];
					rhou_res_sum[icv0][2] += rhs_rhou[icv1][2]; rhou_res_sum[icv1][2] += rhs_rhou[icv0][2];
					rhoE_res_sum[icv0] += rhs_rhoE[icv1]; rhoE_res_sum[icv1] += rhs_rhoE[icv0];
				}	
			}
			
			for (int icv = 0; icv < ncv_mgLevel2; icv++) {
				int non_f = nbocv_i[icv];
				int non_l = nbocv_i[icv+1]-1;
				double nneigh = double(non_l-non_f);
				rhs_rho[icv] = (rho_res_old[icv] + val_smooth_coeff*rho_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][0] = (rhou_res_old[icv][0] + val_smooth_coeff*rhou_res_sum[icv][0]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][1] = (rhou_res_old[icv][1] + val_smooth_coeff*rhou_res_sum[icv][1]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhou[icv][2] = (rhou_res_old[icv][2] + val_smooth_coeff*rhou_res_sum[icv][2]) /(1.0 + val_smooth_coeff*nneigh);
				rhs_rhoE[icv] = (rhoE_res_old[icv] + val_smooth_coeff*rhoE_res_sum[icv]) /(1.0 + val_smooth_coeff*nneigh);
				
			}
			
			for (int ifa = 0; ifa < nfa_b_mgLevel2; ifa++) {
				int icv0 = cvofa[ifa][0];
				rhs_rho[icv0] = rho_res_old[icv0];
				rhs_rhou[icv0][0] = rhou_res_old[icv0][0];
				rhs_rhou[icv0][1] = rhou_res_old[icv0][1];
				rhs_rhou[icv0][2] = rhou_res_old[icv0][2];
				rhs_rhoE[icv0] = rhoE_res_old[icv0];
			}
		}
	}
}

void JoeWithModels::SetSolution(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh) {
	
	double relaxation = getDoubleParam("CORRECTION_RELAXATION", "0.75");
	
	if (iMesh == 0) {
		for (int icv = 0; icv < ncv; icv++) {
			rho[icv] += relaxation*rhs_rho[icv];
			rhou[icv][0] += relaxation*rhs_rhou[icv][0];
			rhou[icv][1] += relaxation*rhs_rhou[icv][1];
			rhou[icv][2] += relaxation*rhs_rhou[icv][2];
			rhoE[icv] += relaxation*rhs_rhoE[icv];
		}
	}
	if (iMesh == 1) {
		for (int icv = 0; icv < ncv_mgLevel1; icv++) {
			rho_mgLevel1[icv] += relaxation*rhs_rho[icv];
			rhou_mgLevel1[icv][0] += relaxation*rhs_rhou[icv][0];
			rhou_mgLevel1[icv][1] += relaxation*rhs_rhou[icv][1];
			rhou_mgLevel1[icv][2] += relaxation*rhs_rhou[icv][2];
			rhoE_mgLevel1[icv] += relaxation*rhs_rhoE[icv];
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			rho_mgLevel2[icv] += relaxation*rhs_rho[icv];
			rhou_mgLevel2[icv][0] += relaxation*rhs_rhou[icv][0];
			rhou_mgLevel2[icv][1] += relaxation*rhs_rhou[icv][1];
			rhou_mgLevel2[icv][2] += relaxation*rhs_rhou[icv][2];
			rhoE_mgLevel2[icv] += relaxation*rhs_rhoE[icv];
		}
	}
	
}


void JoeWithModels::calcDampingSensors() {
	
	double factor = getDoubleParam("MG_DAMPING_FACTOR", "0.0");
	
	
	double MinMuT = 1.0, MaxMuT = 0.0;
	for (int icv = 0; icv < ncv; icv++) {
		MGSensor[icv] = 1.0;
		p1_Und_Lapl[icv] = 0.0;
		p2_Und_Lapl[icv] = 0.0;
		t1_Und_Lapl[icv] = 0.0;
		t2_Und_Lapl[icv] = 0.0;
		mut1_Und_Lapl[icv] = 0.0;
		mut2_Und_Lapl[icv] = 0.0;
	}
	
	for (int ifa = nfa_b; ifa < nfa; ifa++) {
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];
		if ((icv0 < ncv) && (icv1 < ncv)) {
			
			double Pressure_0 = press[icv0];
			double Pressure_1 = press[icv1];
			
			double Temp_0 = temp[icv0];
			double Temp_1 = temp[icv1];
			
			double muT_0 = 0.0;
			double muT_1 = 0.0;
			if (mu_ref > 0.0) {
				muT_0 = InterpolateAtCellCenterFromFaceValues(mut_fa, icv0); // muT[icv0];
				muT_1 = InterpolateAtCellCenterFromFaceValues(mut_fa, icv1); // muT[icv1];
			}
			
			MaxMuT = max(MaxMuT, muT_0);
			MinMuT = min(MinMuT, muT_0);
			MaxMuT = max(MaxMuT, muT_1);
			MinMuT = min(MinMuT, muT_1);
			
			p1_Und_Lapl[icv0] += (Pressure_1 - Pressure_0);
			p1_Und_Lapl[icv1] += (Pressure_0 - Pressure_1);
			
			p2_Und_Lapl[icv0] += (Pressure_0 + Pressure_1);
			p2_Und_Lapl[icv1] += (Pressure_0 + Pressure_1);
			
			t1_Und_Lapl[icv0] += (Temp_1 - Temp_0);
			t1_Und_Lapl[icv1] += (Temp_0 - Temp_1);
			
			t2_Und_Lapl[icv0] += (Temp_0 + Temp_1);
			t2_Und_Lapl[icv1] += (Temp_0 + Temp_1);
			
			if (mu_ref > 0.0) {
				mut1_Und_Lapl[icv0] += (muT_1 - muT_0);
				mut1_Und_Lapl[icv1] += (muT_0 - muT_1);
				
				mut2_Und_Lapl[icv0] += (muT_0 + muT_1);
				mut2_Und_Lapl[icv1] += (muT_0 + muT_1);
			}
		}
	}
	
	for (int icv = 0; icv < ncv; icv++) {
		double PressSensor = 1.0-min(1.0,factor*fabs(p1_Und_Lapl[icv])/p2_Und_Lapl[icv]);
		double TempSensor = 1.0-min(1.0,factor*fabs(t1_Und_Lapl[icv])/t2_Und_Lapl[icv]);
		
		double MuTSensor;
		if (mu_ref > 0.0) MuTSensor = 1.0-min(1.0,factor*fabs(mut1_Und_Lapl[icv])/mut2_Und_Lapl[icv]);
		else MuTSensor = 1.0;
		MGSensor[icv] = PressSensor*TempSensor*MuTSensor;
	}
	
	double max_sensor;
	for (int icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
		max_sensor = 0.0;
		for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			max_sensor = max(max_sensor, MGSensor[icv_fine]);
		}
		MGSensor_mgLevel1[icv_coarse] = max_sensor;
	}
	
	for (int icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++) {
		max_sensor = 0.0;
		for (int iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
			int icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
			max_sensor = max(max_sensor, MGSensor_mgLevel1[icv_fine]);
		}
		MGSensor_mgLevel2[icv_coarse] = max_sensor;
	}
	
	for (int icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++) {
		max_sensor = 0.0;
		for (int iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
			int icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
			max_sensor = max(max_sensor, MGSensor_mgLevel2[icv_fine]);
		}
		MGSensor_mgLevel3[icv_coarse] = max_sensor;
	}
}

void JoeWithModels::SetResidual_Term(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, int iMesh) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	double relaxation = getDoubleParam("FORCING_TERM_RELAXATION", "0.75");
	
	// Add t previous forcing term and do the restriction to the coarse level damping several zones.
	if (iMesh == 0) {
		for (int icv = 0; icv < ncv; icv++) {
			rhs_rho[icv] = rhs_rho[icv];
			rhs_rhou[icv][0] = rhs_rhou[icv][0];
			rhs_rhou[icv][1] = rhs_rhou[icv][1];
			rhs_rhou[icv][2] = rhs_rhou[icv][2];
			rhs_rhoE[icv] = rhs_rhoE[icv];
			
		}
		
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			 rho_TruncError_mgLevel1[icv_coarse] = 0.0;
			 rhou_TruncError_mgLevel1[icv_coarse][0] = 0.0;
			 rhou_TruncError_mgLevel1[icv_coarse][1] = 0.0;
			 rhou_TruncError_mgLevel1[icv_coarse][2] = 0.0;
			 rhoE_TruncError_mgLevel1[icv_coarse] = 0.0;
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				rho_TruncError_mgLevel1[icv_coarse] += relaxation*MGSensor_mgLevel1[icv_coarse]*rhs_rho[icv_fine];
				rhou_TruncError_mgLevel1[icv_coarse][0] += relaxation*MGSensor_mgLevel1[icv_coarse]*rhs_rhou[icv_fine][0];
				rhou_TruncError_mgLevel1[icv_coarse][1] += relaxation*MGSensor_mgLevel1[icv_coarse]*rhs_rhou[icv_fine][1];
				rhou_TruncError_mgLevel1[icv_coarse][2] += relaxation*MGSensor_mgLevel1[icv_coarse]*rhs_rhou[icv_fine][2];
				rhoE_TruncError_mgLevel1[icv_coarse] += relaxation*MGSensor_mgLevel1[icv_coarse]*rhs_rhoE[icv_fine];
			}
		}
	}
	if (iMesh == 1) {
		for (int icv = 0; icv < ncv_mgLevel1; icv++) {
			rhs_rho[icv] = rhs_rho[icv] + rho_TruncError_mgLevel1[icv];
			rhs_rhou[icv][0] = rhs_rhou[icv][0] + rhou_TruncError_mgLevel1[icv][0];
			rhs_rhou[icv][1] = rhs_rhou[icv][1] + rhou_TruncError_mgLevel1[icv][1];
			rhs_rhou[icv][2] = rhs_rhou[icv][2] + rhou_TruncError_mgLevel1[icv][2];
			rhs_rhoE[icv] = rhs_rhoE[icv] + rhoE_TruncError_mgLevel1[icv];
		}
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++) {
			 rho_TruncError_mgLevel2[icv_coarse] = 0.0;
			 rhou_TruncError_mgLevel2[icv_coarse][0] = 0.0;
			 rhou_TruncError_mgLevel2[icv_coarse][1] = 0.0;
			 rhou_TruncError_mgLevel2[icv_coarse][2] = 0.0;
			 rhoE_TruncError_mgLevel2[icv_coarse] = 0.0;
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				rho_TruncError_mgLevel2[icv_coarse] += relaxation*MGSensor_mgLevel2[icv_coarse]*rhs_rho[icv_fine];
				rhou_TruncError_mgLevel2[icv_coarse][0] += relaxation*MGSensor_mgLevel2[icv_coarse]*rhs_rhou[icv_fine][0];
				rhou_TruncError_mgLevel2[icv_coarse][1] += relaxation*MGSensor_mgLevel2[icv_coarse]*rhs_rhou[icv_fine][1];
				rhou_TruncError_mgLevel2[icv_coarse][2] += relaxation*MGSensor_mgLevel2[icv_coarse]*rhs_rhou[icv_fine][2];
				rhoE_TruncError_mgLevel2[icv_coarse] += relaxation*MGSensor_mgLevel2[icv_coarse]*rhs_rhoE[icv_fine];
			}
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			rhs_rho[icv] = rhs_rho[icv] + rho_TruncError_mgLevel2[icv];
			rhs_rhou[icv][0] = rhs_rhou[icv][0] + rhou_TruncError_mgLevel2[icv][0];
			rhs_rhou[icv][1] = rhs_rhou[icv][1] + rhou_TruncError_mgLevel2[icv][1];
			rhs_rhou[icv][2] = rhs_rhou[icv][2] + rhou_TruncError_mgLevel2[icv][2];
			rhs_rhoE[icv] = rhs_rhoE[icv] + rhoE_TruncError_mgLevel2[icv];
		}
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++) {
			 rho_TruncError_mgLevel3[icv_coarse] = 0.0;
			 rhou_TruncError_mgLevel3[icv_coarse][0] = 0.0;
			 rhou_TruncError_mgLevel3[icv_coarse][1] = 0.0;
			 rhou_TruncError_mgLevel3[icv_coarse][2] = 0.0;
			 rhoE_TruncError_mgLevel3[icv_coarse] = 0.0;
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				rho_TruncError_mgLevel3[icv_coarse] += relaxation*MGSensor_mgLevel3[icv_coarse]*rhs_rho[icv_fine];
				rhou_TruncError_mgLevel3[icv_coarse][0] += relaxation*MGSensor_mgLevel3[icv_coarse]*rhs_rhou[icv_fine][0];
				rhou_TruncError_mgLevel3[icv_coarse][1] += relaxation*MGSensor_mgLevel3[icv_coarse]*rhs_rhou[icv_fine][1];
				rhou_TruncError_mgLevel3[icv_coarse][2] += relaxation*MGSensor_mgLevel3[icv_coarse]*rhs_rhou[icv_fine][2];
				rhoE_TruncError_mgLevel3[icv_coarse] += relaxation*MGSensor_mgLevel3[icv_coarse]*rhs_rhoE[icv_fine];
			}
		}
	}
}

void JoeWithModels::SetForcing_Term(double *rhs_rho_fine, double (*rhs_rhou_fine)[3], double *rhs_rhoE_fine, double *rhs_rho_coarse, 
																		double (*rhs_rhou_coarse)[3], double *rhs_rhoE_coarse, int iMesh) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	
	if (iMesh == 1) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			rho_TruncError_mgLevel1[icv_coarse] = rho_TruncError_mgLevel1[icv_coarse] - rhs_rho_coarse[icv_coarse];
			rhou_TruncError_mgLevel1[icv_coarse][0] = rhou_TruncError_mgLevel1[icv_coarse][0] - rhs_rhou_coarse[icv_coarse][0];
			rhou_TruncError_mgLevel1[icv_coarse][1] = rhou_TruncError_mgLevel1[icv_coarse][1] - rhs_rhou_coarse[icv_coarse][1];
			rhou_TruncError_mgLevel1[icv_coarse][2] = rhou_TruncError_mgLevel1[icv_coarse][2] - rhs_rhou_coarse[icv_coarse][2];
			rhoE_TruncError_mgLevel1[icv_coarse] = rhoE_TruncError_mgLevel1[icv_coarse] - rhs_rhoE_coarse[icv_coarse];		
		}
	}
	if (iMesh == 2) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++) {
			rho_TruncError_mgLevel2[icv_coarse] = rho_TruncError_mgLevel2[icv_coarse] - rhs_rho_coarse[icv_coarse];
			rhou_TruncError_mgLevel2[icv_coarse][0] = rhou_TruncError_mgLevel2[icv_coarse][0] - rhs_rhou_coarse[icv_coarse][0];
			rhou_TruncError_mgLevel2[icv_coarse][1] = rhou_TruncError_mgLevel2[icv_coarse][1] - rhs_rhou_coarse[icv_coarse][1];
			rhou_TruncError_mgLevel2[icv_coarse][2] = rhou_TruncError_mgLevel2[icv_coarse][2] - rhs_rhou_coarse[icv_coarse][2];
			rhoE_TruncError_mgLevel2[icv_coarse] = rhoE_TruncError_mgLevel2[icv_coarse] - rhs_rhoE_coarse[icv_coarse];		
		}	
	}
	if (iMesh == 3) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++) {
			rho_TruncError_mgLevel3[icv_coarse] = rho_TruncError_mgLevel3[icv_coarse] - rhs_rho_coarse[icv_coarse];
			rhou_TruncError_mgLevel3[icv_coarse][0] = rhou_TruncError_mgLevel3[icv_coarse][0] - rhs_rhou_coarse[icv_coarse][0];
			rhou_TruncError_mgLevel3[icv_coarse][1] = rhou_TruncError_mgLevel3[icv_coarse][1] - rhs_rhou_coarse[icv_coarse][1];
			rhou_TruncError_mgLevel3[icv_coarse][2] = rhou_TruncError_mgLevel3[icv_coarse][2] - rhs_rhou_coarse[icv_coarse][2];
			rhoE_TruncError_mgLevel3[icv_coarse] = rhoE_TruncError_mgLevel3[icv_coarse] - rhs_rhoE_coarse[icv_coarse];		
			
		}
	}
}

void JoeWithModels::SetRestricted_Solution(int iMesh) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	double icv_volume_parent, icv_volume_children;
	
	if (iMesh == 0) {
		
		// Compute the interpolation on the coarse grid
		for (icv_coarse = 0; icv_coarse < ncv_g_mgLevel1; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel1[icv_coarse];
			rho_mgLevel1[icv_coarse] = 0.0;
			rhou_mgLevel1[icv_coarse][0] = 0.0;
			rhou_mgLevel1[icv_coarse][1] = 0.0;
			rhou_mgLevel1[icv_coarse][2] = 0.0;
			rhoE_mgLevel1[icv_coarse] = 0.0;
			gamma_mgLevel1[icv_coarse] = 0.0;
			RoM_mgLevel1[icv_coarse] = 0.0;
			muT_mgLevel1[icv_coarse] = 0.0;
			sos_mgLevel1[icv_coarse] = 0.0;
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				icv_volume_children = cv_volume[icv_fine];
				
				rho_mgLevel1[icv_coarse] += rho[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel1[icv_coarse][0] += rhou[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel1[icv_coarse][1] += rhou[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel1[icv_coarse][2] += rhou[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_mgLevel1[icv_coarse] += rhoE[icv_fine]*icv_volume_children/icv_volume_parent;
				gamma_mgLevel1[icv_coarse] += gamma[icv_fine]*icv_volume_children/icv_volume_parent;
				RoM_mgLevel1[icv_coarse] += RoM[icv_fine]*icv_volume_children/icv_volume_parent;
				sos_mgLevel1[icv_coarse] += sos[icv_fine]*icv_volume_children/icv_volume_parent;
				muT_mgLevel1[icv_coarse] += InterpolateAtCellCenterFromFaceValues(mut_fa, icv_fine)*icv_volume_children/icv_volume_parent;
			}
		}
		
		
		string ProjParallel = getStringParam("PROJECT_PARALLEL","NO");
		if (ProjParallel == "YES") {
			
			SetUpdateGhost_Solution(MESH_1);
			/*
			 // 	Project the interpolation to the fine auxiliar interpolated grid for doing the Send/Receive
			 for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			 for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
			 icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			 rho_interp[icv_fine] = rho_mgLevel1[icv_coarse];
			 rhou_interp[icv_fine][0] = rhou_mgLevel1[icv_coarse][0];
			 rhou_interp[icv_fine][1] = rhou_mgLevel1[icv_coarse][1];
			 rhou_interp[icv_fine][2] = rhou_mgLevel1[icv_coarse][2];
			 rhoE_interp[icv_fine] = rhoE_mgLevel1[icv_coarse];
			 }
			 }
			 
			 // Transfer the interpolated solution
			 updateCvData(rhou_interp, REPLACE_ROTATE_DATA);
			 updateCvData(rho_interp,  REPLACE_DATA);
			 updateCvData(rhoE_interp, REPLACE_DATA);
			 
			 // Update the solution the boundary elements
			 for (icv_coarse = ncv_mgLevel1; icv_coarse < ncv_g_mgLevel1; icv_coarse++) {
			 icv_fine = Children_CV_mgLevel1[icv_coarse][0];
			 rho_mgLevel1[icv_coarse] = rho_interp[icv_fine];
			 rhou_mgLevel1[icv_coarse][0] = rhou_interp[icv_fine][0];
			 rhou_mgLevel1[icv_coarse][1] = rhou_interp[icv_fine][1];
			 rhou_mgLevel1[icv_coarse][2] = rhou_interp[icv_fine][2];
			 rhoE_mgLevel1[icv_coarse] = rhoE_interp[icv_fine];
			 }*/
			
		}
		
	}
	
	if (iMesh == 1) {
		for (icv_coarse = 0; icv_coarse < ncv_g_mgLevel2; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel2[icv_coarse];
			
			rho_mgLevel2[icv_coarse] = 0.0;
			rhou_mgLevel2[icv_coarse][0] = 0.0;
			rhou_mgLevel2[icv_coarse][1] = 0.0;
			rhou_mgLevel2[icv_coarse][2] = 0.0;
			rhoE_mgLevel2[icv_coarse] = 0.0;
			gamma_mgLevel2[icv_coarse] = 0.0;
			RoM_mgLevel2[icv_coarse] = 0.0;
			muT_mgLevel2[icv_coarse] = 0.0;
			sos_mgLevel2[icv_coarse] = 0.0;
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				icv_volume_children = cv_volume_mgLevel1[icv_fine];
				
				rho_mgLevel2[icv_coarse] += rho_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel2[icv_coarse][0] += rhou_mgLevel1[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel2[icv_coarse][1] += rhou_mgLevel1[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel2[icv_coarse][2] += rhou_mgLevel1[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_mgLevel2[icv_coarse] += rhoE_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				gamma_mgLevel2[icv_coarse] += gamma_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				RoM_mgLevel2[icv_coarse] += RoM_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				muT_mgLevel2[icv_coarse] += muT_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
				sos_mgLevel2[icv_coarse] += sos_mgLevel1[icv_fine]*icv_volume_children/icv_volume_parent;
			}
		}
	}
	if (iMesh == 2) {
		for (icv_coarse = 0; icv_coarse < ncv_g_mgLevel3; icv_coarse++) {
			icv_volume_parent = cv_volume_mgLevel3[icv_coarse];
			
			rho_mgLevel3[icv_coarse] = 0.0;
			rhou_mgLevel3[icv_coarse][0] = 0.0;
			rhou_mgLevel3[icv_coarse][1] = 0.0;
			rhou_mgLevel3[icv_coarse][2] = 0.0;
			rhoE_mgLevel3[icv_coarse] = 0.0;
			gamma_mgLevel3[icv_coarse] = 0.0;
			RoM_mgLevel3[icv_coarse] = 0.0;
			muT_mgLevel3[icv_coarse] = 0.0;
			sos_mgLevel3[icv_coarse] = 0.0;
			
			
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				icv_volume_children = cv_volume_mgLevel2[icv_fine];
				
				rho_mgLevel3[icv_coarse] += rho_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel3[icv_coarse][0] += rhou_mgLevel2[icv_fine][0]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel3[icv_coarse][1] += rhou_mgLevel2[icv_fine][1]*icv_volume_children/icv_volume_parent;
				rhou_mgLevel3[icv_coarse][2] += rhou_mgLevel2[icv_fine][2]*icv_volume_children/icv_volume_parent;
				rhoE_mgLevel3[icv_coarse] += rhoE_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				gamma_mgLevel3[icv_coarse] += gamma_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				RoM_mgLevel3[icv_coarse] += RoM_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				muT_mgLevel3[icv_coarse] += muT_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
				sos_mgLevel3[icv_coarse] += sos_mgLevel2[icv_fine]*icv_volume_children/icv_volume_parent;
			}
		}
	}
}


void JoeWithModels::SetUpdateGhost_Solution(int iMesh) {
	
	if (iMesh == 1) {
		// 	Project the interpolation to the fine auxiliar interpolated grid for doing the Send/Receive
		for (int icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				rho_interp[icv_fine] = rho_mgLevel1[icv_coarse];
				rhou_interp[icv_fine][0] = rhou_mgLevel1[icv_coarse][0];
				rhou_interp[icv_fine][1] = rhou_mgLevel1[icv_coarse][1];
				rhou_interp[icv_fine][2] = rhou_mgLevel1[icv_coarse][2];
				rhoE_interp[icv_fine] = rhoE_mgLevel1[icv_coarse];
			}
		}
		
		// Transfer the interpolated solution
    	updateCvData(rhou_interp, REPLACE_ROTATE_DATA);
    	updateCvData(rho_interp,  REPLACE_DATA);
    	updateCvData(rhoE_interp, REPLACE_DATA);
		
		// Update the solution the boundary elements
		for (int icv_coarse = ncv_mgLevel1; icv_coarse < ncv_g_mgLevel1; icv_coarse++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][0];
			rho_mgLevel1[icv_coarse] = rho_interp[icv_fine];
			rhou_mgLevel1[icv_coarse][0] = rhou_interp[icv_fine][0];
			rhou_mgLevel1[icv_coarse][1] = rhou_interp[icv_fine][1];
			rhou_mgLevel1[icv_coarse][2] = rhou_interp[icv_fine][2];
			rhoE_mgLevel1[icv_coarse] = rhoE_interp[icv_fine];
		}
	}
	
	if (iMesh == 2) {
		// 	Project the interpolation to the fine auxiliar interpolated grid for doing the Send/Receive
		for (int icv_mgLevel2 = 0; icv_mgLevel2 < ncv_mgLevel2; icv_mgLevel2++) {
			for (int iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren++) {
				int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren];
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
					int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
					rho_interp[icv_fine] = rho_mgLevel1[icv_mgLevel2];
					rhou_interp[icv_fine][0] = rhou_mgLevel1[icv_mgLevel2][0];
					rhou_interp[icv_fine][1] = rhou_mgLevel1[icv_mgLevel2][1];
					rhou_interp[icv_fine][2] = rhou_mgLevel1[icv_mgLevel2][2];
					rhoE_interp[icv_fine] = rhoE_mgLevel1[icv_mgLevel2];
				}
			}
		}
		
		// Transfer the interpolated solution
    	updateCvData(rhou_interp, REPLACE_ROTATE_DATA);
    	updateCvData(rho_interp,  REPLACE_DATA);
    	updateCvData(rhoE_interp, REPLACE_DATA);
		
		// Update the solution the boundary elements
		for (int icv_mgLevel2 = ncv_mgLevel2; icv_mgLevel2 < ncv_g_mgLevel2; icv_mgLevel2++) {
			int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
			int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
			rho_mgLevel1[icv_mgLevel2] = rho_interp[icv_fine];
			rhou_mgLevel1[icv_mgLevel2][0] = rhou_interp[icv_fine][0];
			rhou_mgLevel1[icv_mgLevel2][1] = rhou_interp[icv_fine][1];
			rhou_mgLevel1[icv_mgLevel2][2] = rhou_interp[icv_fine][2];
			rhoE_mgLevel1[icv_mgLevel2] = rhoE_interp[icv_fine];
		}
	}
	
	if (iMesh == 3) {
		// 	Project the interpolation to the fine auxiliar interpolated grid for doing the Send/Receive
		for (int icv_mgLevel3 = 0; icv_mgLevel3 < ncv_mgLevel3; icv_mgLevel3++) {
			for (int iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_mgLevel3]; iChildren++) {
				int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][iChildren];
				for (int iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_mgLevel2]; iChildren++) {
					int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][iChildren];
					for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_mgLevel1]; iChildren++) {
						int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][iChildren];
						rho_interp[icv_fine] = rho_mgLevel1[icv_mgLevel3];
						rhou_interp[icv_fine][0] = rhou_mgLevel1[icv_mgLevel3][0];
						rhou_interp[icv_fine][1] = rhou_mgLevel1[icv_mgLevel3][1];
						rhou_interp[icv_fine][2] = rhou_mgLevel1[icv_mgLevel3][2];
						rhoE_interp[icv_fine] = rhoE_mgLevel1[icv_mgLevel3];
					}
				}
			}
		}
		
		// Transfer the interpolated solution
    	updateCvData(rhou_interp, REPLACE_ROTATE_DATA);
    	updateCvData(rho_interp,  REPLACE_DATA);
    	updateCvData(rhoE_interp, REPLACE_DATA);
		
		// Update the solution the boundary elements
		for (int icv_mgLevel3 = ncv_mgLevel3; icv_mgLevel3 < ncv_g_mgLevel3; icv_mgLevel3++) {
			int icv_mgLevel2 = Children_CV_mgLevel3[icv_mgLevel3][0];
			int icv_mgLevel1 = Children_CV_mgLevel2[icv_mgLevel2][0];
			int icv_fine = Children_CV_mgLevel1[icv_mgLevel1][0];
			rho_mgLevel1[icv_mgLevel3] = rho_interp[icv_fine];
			rhou_mgLevel1[icv_mgLevel3][0] = rhou_interp[icv_fine][0];
			rhou_mgLevel1[icv_mgLevel3][1] = rhou_interp[icv_fine][1];
			rhou_mgLevel1[icv_mgLevel3][2] = rhou_interp[icv_fine][2];
			rhoE_mgLevel1[icv_mgLevel3] = rhoE_interp[icv_fine];
		}
	}
}



void JoeWithModels::SetProjected_Solution(int iMesh) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	
	if (iMesh == 1) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
				
				rho[icv_fine] = rho_mgLevel1[icv_coarse];
				rhou[icv_fine][0] = rhou_mgLevel1[icv_coarse][0];
				rhou[icv_fine][1] = rhou_mgLevel1[icv_coarse][1];
				rhou[icv_fine][2] = rhou_mgLevel1[icv_coarse][2];
				rhoE[icv_fine] = rhoE_mgLevel1[icv_coarse];
			}
		}
	}
	if (iMesh == 2) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse++) {
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel2[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				
				rho_mgLevel1[icv_fine] = rho_mgLevel2[icv_coarse];
				rhou_mgLevel1[icv_fine][0] = rhou_mgLevel2[icv_coarse][0];
				rhou_mgLevel1[icv_fine][1] = rhou_mgLevel2[icv_coarse][1];
				rhou_mgLevel1[icv_fine][2] = rhou_mgLevel2[icv_coarse][2];
				rhoE_mgLevel1[icv_fine] = rhoE_mgLevel2[icv_coarse];
			}
		}
	}
	if (iMesh == 3) {
		for (icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse++) {
			for (iChildren = 0; iChildren < nChildren_CV_mgLevel3[icv_coarse]; iChildren++) {
				icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				
				rho_mgLevel2[icv_fine] = rho_mgLevel3[icv_coarse];
				rhou_mgLevel2[icv_fine][0] = rhou_mgLevel3[icv_coarse][0];
				rhou_mgLevel2[icv_fine][1] = rhou_mgLevel3[icv_coarse][1];
				rhou_mgLevel2[icv_fine][2] = rhou_mgLevel3[icv_coarse][2];
				rhoE_mgLevel2[icv_fine] = rhoE_mgLevel3[icv_coarse];
			}
		}
	}	
}

void JoeWithModels::SetProjected_Residual(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double *rhs_rho_mgLevel1, double (*rhs_rhou_mgLevel1)[3], double *rhs_rhoE_mgLevel1) {
	unsigned long icv_fine, icv_coarse;
	unsigned short iChildren;
	
	for (icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse++) {
		for (iChildren = 0; iChildren < nChildren_CV_mgLevel1[icv_coarse]; iChildren++) {
			icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			
			rhs_rho[icv_fine] = rhs_rho_mgLevel1[icv_coarse];
			rhs_rhou[icv_fine][0] = rhs_rhou_mgLevel1[icv_coarse][0];
			rhs_rhou[icv_fine][1] = rhs_rhou_mgLevel1[icv_coarse][1];
			rhs_rhou[icv_fine][2] = rhs_rhou_mgLevel1[icv_coarse][2];
			rhs_rhoE[icv_fine] = rhs_rhoE_mgLevel1[icv_coarse];
			
		}
	}	
}

void JoeWithModels::initializeFromFineGrid() {
	
	double my_volume_sum = 0.0;
	for (int icv = 0; icv<ncv; icv++)
		my_volume_sum += cv_volume[icv];
	double TotalVolume;
	MPI_Allreduce(&my_volume_sum, &TotalVolume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	
	int max_neighbor = 150;
	
	// Relation between the domain length and the aglomerated element.
	double Maxlength = getDoubleParam("MAX_LENGTH_AGGLOMERATION", "1.0E8");
	
	// Allocate some structures for doing the agglomeration of the fine grid
	Parent_CV = new int[ncv_g];
	Agglomerate = new bool[ncv_g];
	for (int icv = 0; icv < ncv_g; icv++)
		Agglomerate[icv] = false;
	
	// The ghost control volume can not be aglomerated using the information of this domain
	for (int icv = ncv; icv < ncv_g; icv++)
		Agglomerate[icv] = true;
	
	nChildren_CV_mgLevel1 = new int[ncv_g];
	Children_CV_mgLevel1 = new int *[ncv_g];
	for (int icv = 0; icv < ncv_g; icv++)
		Children_CV_mgLevel1[icv] = new int[max_neighbor];
	
	list<FaZone>::iterator zone_fine;
	list<FaZone>::iterator zone_coarse;
	
	string twoDimwithSymmetry = getStringParam("2D_WITH_SYMMETRY","NO");
	
	// Compute the agglomeration struture without compact storage method
	ncv_mgLevel1 = 0;
	
	// Agglomerate the linelets
	string SemiName = getStringParam("SEMICOARSENING","NO");
	if (SemiName == "YES") {
		for (int iSources = 0; iSources < nSources-3; iSources++) {
			for (int index_cv = 0; index_cv < linelet_cv[iSources].size()-5; index_cv++) {
				int icv0 = linelet_cv[iSources][index_cv];
				int icv1 = linelet_cv[iSources][index_cv+1];
				int icv2 = linelet_cv[iSources][index_cv+2];
				int icv3 = linelet_cv[iSources][index_cv+3];
				if ((Agglomerate[icv0] == false) && (Agglomerate[icv1] == false) &&
					(Agglomerate[icv2] == false) && (Agglomerate[icv3] == false)) {
					Parent_CV[icv0] = ncv_mgLevel1; Agglomerate[icv0] = true;
					Parent_CV[icv1] = ncv_mgLevel1; Agglomerate[icv1] = true;
					Parent_CV[icv2] = ncv_mgLevel1; Agglomerate[icv2] = true;
					Parent_CV[icv3] = ncv_mgLevel1; Agglomerate[icv3] = true; 
					nChildren_CV_mgLevel1[ncv_mgLevel1] = 4;
					Children_CV_mgLevel1[ncv_mgLevel1][0] = icv0;
					Children_CV_mgLevel1[ncv_mgLevel1][1] = icv1;
					Children_CV_mgLevel1[ncv_mgLevel1][2] = icv2;
					Children_CV_mgLevel1[ncv_mgLevel1][3] = icv3;
					ncv_mgLevel1++;
				}
			}
		}
	}
	
	// Second agglomerate the inner elements
	string HexaAgglomeration = getStringParam("ONLY_HEXA_AGGLOMERATION","NO");
	if (HexaAgglomeration == "NO") {
		
		for (int icv = 0; icv < ncv; icv++)
			
			if (Agglomerate[icv] == false) {
				
				// Add the seed control volume
				Parent_CV[icv] = ncv_mgLevel1; 
				Agglomerate[icv] = true;
				nChildren_CV_mgLevel1[ncv_mgLevel1] = 1;
				Children_CV_mgLevel1[ncv_mgLevel1][0] = icv;
				
				// Add the indirect neighbor control volume (neighbor that shares a point with the original element)
				int noc_f = nbocv_i[icv];
				int noc_l = nbocv_i[icv+1]-1;
				for (int noc = noc_f; noc <= noc_l; noc++) {
					int nbocv = nbocv_v[noc];
					if (Agglomerate[nbocv] == false) {
						Parent_CV[nbocv] = ncv_mgLevel1;
						Agglomerate[nbocv] = true;
						Children_CV_mgLevel1[ncv_mgLevel1][nChildren_CV_mgLevel1[ncv_mgLevel1]] = nbocv;
						nChildren_CV_mgLevel1[ncv_mgLevel1]++;
					}
				}
				
				vector<int> First_Neighbor_Points;
				vector<int> Second_Neighbor_Points;
				vector<int> Third_Neighbor_Points;
				vector<int> Fourth_Neighbor_Points;
				vector<int> Fifth_Neighbor_Points;
				vector<int> Sixth_Neighbor_Points;
				
				First_Neighbor_Points.push_back(icv);
				for (int jNode = nbocv_i[icv]+1; jNode <= nbocv_i[icv+1]-1; jNode++) {
					int kPoint = nbocv_v[jNode];
					if (kPoint < ncv) First_Neighbor_Points.push_back(kPoint);
				}
				
				
				for (int iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++) {
					int jPoint = First_Neighbor_Points[iNeighbor];
					for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
						int kPoint = nbocv_v[jNode];
						if (kPoint < ncv) Second_Neighbor_Points.push_back(kPoint);
					}
				}
				
				
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++) {
					int jPoint = Second_Neighbor_Points[iNeighbor];
					for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
						int kPoint = nbocv_v[jNode];
						if (kPoint < ncv) Third_Neighbor_Points.push_back(kPoint);
					}
				}
				
				
				for (int iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++) {
					int jPoint = Third_Neighbor_Points[iNeighbor];
					for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
						int kPoint = nbocv_v[jNode];
						if (kPoint < ncv) Fourth_Neighbor_Points.push_back(kPoint);
					}
				}
				
				
				/*				for (int iNeighbor = 0; iNeighbor <	Fourth_Neighbor_Points.size(); iNeighbor ++) {
				 int jPoint = Fourth_Neighbor_Points[iNeighbor];
				 for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
				 int kPoint = nbocv_v[jNode];
				 if (kPoint < ncv) Fifth_Neighbor_Points.push_back(kPoint);
				 }
				 }
				 
				 
				 for (int iNeighbor = 0; iNeighbor <	Fifth_Neighbor_Points.size(); iNeighbor ++) {
				 int jPoint = Fifth_Neighbor_Points[iNeighbor];
				 for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
				 int kPoint = nbocv_v[jNode];
				 if (kPoint < ncv) Sixth_Neighbor_Points.push_back(kPoint);
				 }
				 }*/
				
				// Storage the points of the seed node, add the nodes of the original element
				vector<int> Seed_Points;
				for (int foc=faocv_i[icv]; foc<=faocv_i[icv + 1] - 1; foc++)
					for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
						Seed_Points.push_back(noofa_v[nof]);
				
				// Add the elements that shares a vertex with the new control volume
				vector<int> add_cv;
				for (int iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++) {
					int firstcv = First_Neighbor_Points[iNeighbor];
					for (int foc=faocv_i[firstcv]; foc<=faocv_i[firstcv + 1] - 1; foc++)
						for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
							for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
								if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(firstcv);
				}
				
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++) {
					int secondcv = Second_Neighbor_Points[iNeighbor];
					for (int foc=faocv_i[secondcv]; foc<=faocv_i[secondcv + 1] - 1; foc++)
						for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
							for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
								if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(secondcv);
				}
				
				for (int iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++) {
					int thirdcv = Third_Neighbor_Points[iNeighbor];
					for (int foc=faocv_i[thirdcv]; foc<=faocv_i[thirdcv + 1] - 1; foc++)
						for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
							for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
								if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(thirdcv);
				}
				
				for (int iNeighbor = 0; iNeighbor <	Fourth_Neighbor_Points.size(); iNeighbor ++) {
					int fourthcv = Fourth_Neighbor_Points[iNeighbor];
					for (int foc=faocv_i[fourthcv]; foc<=faocv_i[fourthcv + 1] - 1; foc++)
						for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
							for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
								if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(fourthcv);
				}
				
				/*  			for (int iNeighbor = 0; iNeighbor <	Fifth_Neighbor_Points.size(); iNeighbor ++) {
				 int fifthcv = Fifth_Neighbor_Points[iNeighbor];
				 for (int foc=faocv_i[fifthcv]; foc<=faocv_i[fifthcv + 1] - 1; foc++)
				 for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
				 for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
				 if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(fifthcv);
				 }
				 
				 for (int iNeighbor = 0; iNeighbor <	Sixth_Neighbor_Points.size(); iNeighbor ++) {
				 int sixthcv = Sixth_Neighbor_Points[iNeighbor];
				 for (int foc=faocv_i[sixthcv]; foc<=faocv_i[sixthcv + 1] - 1; foc++)
				 for (int nof = noofa_i[faocv_v[foc]]; nof<noofa_i[faocv_v[foc]+1]; nof++)
				 for (int iSeedPoint = 0; iSeedPoint < Seed_Points.size(); iSeedPoint ++)
				 if (Seed_Points[iSeedPoint] == noofa_v[nof]) add_cv.push_back(sixthcv);
				 }*/
				
				std::sort( add_cv.begin(), add_cv.end() );
				std::vector<int>::iterator new_end_pos;
				new_end_pos = std::unique( add_cv.begin(), add_cv.end() );
				add_cv.erase( new_end_pos, add_cv.end() );
				
				
				double TempVol = cv_volume[icv];	
				
				for (int kNode = 0; kNode <	add_cv.size(); kNode ++) {
					for (int iNode = 0; iNode <	add_cv.size(); iNode ++) {
						int nbocv = add_cv[iNode];
						
						// Check the order for adding new elements (we add elements that have a 
						// common face, with the lelements that have been previously added )
						bool check = false;
						for (int iChildren = 0; iChildren < nChildren_CV_mgLevel1[ncv_mgLevel1]; iChildren++) {
							int TestCV = Children_CV_mgLevel1[ncv_mgLevel1][iChildren];
							for (int jNode = nbocv_i[TestCV]+1; jNode <= nbocv_i[TestCV+1]-1; jNode++)
								if (nbocv_v[jNode] == nbocv) { check = true; break; }
						}
						
						if ((Agglomerate[nbocv] == false) && (TotalVolume * Maxlength > TempVol) && (check)) {
							TempVol += cv_volume[nbocv];
							Parent_CV[nbocv] = ncv_mgLevel1;
							Agglomerate[nbocv] = true;
							Children_CV_mgLevel1[ncv_mgLevel1][nChildren_CV_mgLevel1[ncv_mgLevel1]] = nbocv;
							nChildren_CV_mgLevel1[ncv_mgLevel1]++;
						}
					}
				}	
				
				// Update number coarse grid control volumes
				ncv_mgLevel1++;
			}
	}
	
	// Fast algorithm for hexa agglomeration
	else { 
		
		for (int icv = 0; icv < ncv; icv++) {
			if (Agglomerate[icv] == false) {
				
				// Add the seed control volume
				Parent_CV[icv] = ncv_mgLevel1; 
				Agglomerate[icv] = true;
				nChildren_CV_mgLevel1[ncv_mgLevel1] = 1;
				Children_CV_mgLevel1[ncv_mgLevel1][0] = icv;
				
				// Add the first neighbor control volume
				int noc_f = nbocv_i[icv];
				int noc_l = nbocv_i[icv+1]-1;
				for (int noc = noc_f; noc <= noc_l; noc++) {
					int nbocv = nbocv_v[noc];
					if (Agglomerate[nbocv] == false) {
						Parent_CV[nbocv] = ncv_mgLevel1;
						Agglomerate[nbocv] = true;
						Children_CV_mgLevel1[ncv_mgLevel1][nChildren_CV_mgLevel1[ncv_mgLevel1]] = nbocv;
						nChildren_CV_mgLevel1[ncv_mgLevel1]++;
					}
				}
				
				// Add the second neighbor control volume
				vector<unsigned long> First_Neighbor_Points;
				First_Neighbor_Points.push_back(icv);
				
				// Localize first neighbors and create a list
				for (int noc = nbocv_i[icv]+1; noc <= nbocv_i[icv+1]-1; noc++)
					First_Neighbor_Points.push_back(nbocv_v[noc]);
				
				// Identify a list with the second neighbors and its origin
				vector<int> Second_Neighbor_Points;
				vector<int> Second_Origin_Points;
				
				for (int iNode = nbocv_i[icv]+1; iNode <= nbocv_i[icv+1]-1; iNode++) {
					int jPoint = nbocv_v[iNode];
					if (jPoint < ncv) {
						for (int jNode = nbocv_i[jPoint]+1; jNode <= nbocv_i[jPoint+1]-1; jNode++) {
							int kPoint = nbocv_v[jNode];
							
							// Verify if the element belongs to the list of the first neighbors
							bool SecondNeighborSeed = true;
							for (int iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
								if (kPoint == First_Neighbor_Points[iNeighbor]) {
									SecondNeighborSeed = false;
									break;
								}
							if (SecondNeighborSeed) {
								Second_Neighbor_Points.push_back(kPoint);
								Second_Origin_Points.push_back(jPoint);
							}
						}
					}
				}
				
				vector<int> Repeated_Neighbor_Points;					
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
					for (int jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
						if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) && (iNeighbor != jNeighbor))
							Repeated_Neighbor_Points.push_back(Second_Neighbor_Points[iNeighbor]);
				
				for (int iNode = 0; iNode <	Repeated_Neighbor_Points.size(); iNode ++) {
					int nbocv = Repeated_Neighbor_Points[iNode];
					if (Agglomerate[nbocv] == false) {
						Parent_CV[nbocv] = ncv_mgLevel1;
						Agglomerate[nbocv] = true;
						Children_CV_mgLevel1[ncv_mgLevel1][nChildren_CV_mgLevel1[ncv_mgLevel1]] = nbocv;
						nChildren_CV_mgLevel1[ncv_mgLevel1]++;
						
					}
				}
				
				// Update number coarse grid control volumes
				ncv_mgLevel1++;
			}
		}
	}
	
	
	// In case there is a control volume which has not being agglomerated
	for (int icv = 0; icv < ncv; icv++)
		if (Agglomerate[icv] == false) {
			Parent_CV[icv] = ncv_mgLevel1;
			Agglomerate[icv] = true;
			nChildren_CV_mgLevel1[ncv_mgLevel1] = 1;
			Children_CV_mgLevel1[ncv_mgLevel1][0] = icv;
			ncv_mgLevel1++;
		}
	
	
	// Store the number of ghost control volumes of the coarse grid
	ncv_g_mgLevel1 = ncv_mgLevel1;
	
	// Add the ghost control volume at the end of the agglomerated volumes
	for (int icv = ncv; icv < ncv_g; icv++) {
		Parent_CV[icv] = ncv_g_mgLevel1;
		Agglomerate[icv] = true;
		nChildren_CV_mgLevel1[ncv_g_mgLevel1] = 1;
		Children_CV_mgLevel1[ncv_g_mgLevel1][0] = icv;
		ncv_g_mgLevel1++;
	}
	
	double MeanChildren = 0.0;
	int MaxChildren = 0, MinChildren = max_neighbor;
	for (int icv = 0; icv < ncv_mgLevel1; icv++) {
		MeanChildren += nChildren_CV_mgLevel1[icv];
		MaxChildren = max(nChildren_CV_mgLevel1[icv], MaxChildren);
		MinChildren = min(nChildren_CV_mgLevel1[icv], MinChildren);
	}
	
	int my_rank; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	cout << "Domain: " << my_rank << ". CV of the MG level: " << ncv_mgLevel1 << ". Ghost CV of the MG level: " 
	<< ncv_g_mgLevel1-ncv_mgLevel1 << ". Reduction factor: 1/" <<MeanChildren/double(ncv_mgLevel1)<< "." <<endl;
	MPI_Barrier( MPI_COMM_WORLD );
//	cout << "Domain: " << my_rank << ". MaxChildren: " << MaxChildren << ". MinChildren: " 
//	<< MinChildren << ". MeanChildren: " <<int(MeanChildren/double(ncv_mgLevel1))<< "." <<endl;
	
	MPI_Barrier( MPI_COMM_WORLD );
	
	// Compute the neighbors of a control volume without compact storage method
	int **neighbor, *nNeighbor;
	neighbor = new int* [ncv_g_mgLevel1];
	for (int icv = 0; icv < ncv_g_mgLevel1; icv++)
		neighbor[icv] = new int[max_neighbor];
	nNeighbor = new int[ncv_g_mgLevel1];
	bool add_cv;
	
	
	for (int icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse ++) {
		nNeighbor[icv_coarse] = 0;
		for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel1[icv_coarse]; iChildren ++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			int nboc_f = nbocv_i[icv_fine];
			int nboc_l = nbocv_i[icv_fine + 1] - 1;
			
			for (int noc = nboc_f + 1; noc <= nboc_l; noc++) {
				int icv_fine_neighbor = nbocv_v[noc];
				int icv_parent = Parent_CV[icv_fine_neighbor];
				if (icv_parent != icv_coarse) {
					add_cv = true;
					for (int iNeighbor = 0; iNeighbor < nNeighbor[icv_coarse]; iNeighbor++)
						if (neighbor[icv_coarse][iNeighbor] == icv_parent) {
							add_cv = false; break; 
						}
					if (add_cv) {
						neighbor[icv_coarse][nNeighbor[icv_coarse]] = icv_parent;
						nNeighbor[icv_coarse]++;
					}
				}
			}
		}
	}
	
	
	// Write the neighbors struture using a compact storage method
	assert(nbocv_i_mgLevel1==NULL);
	assert(nbocv_v_mgLevel1==NULL);
	
	nbocv_i_mgLevel1 = new int[ncv_mgLevel1+1];
	int nnbocv_v_mgLevel1 = ncv_mgLevel1;
	for (int icv = 0; icv < ncv_mgLevel1; icv++) {
		nbocv_i_mgLevel1[icv+1] = 1;
		nnbocv_v_mgLevel1 += nNeighbor[icv];
	}
	nbocv_v_mgLevel1 = new int[nnbocv_v_mgLevel1];
	
	int inbocv_v_mgLevel1 = 0;
	for (int icv = 0; icv < ncv_mgLevel1; icv++) {
		nbocv_i_mgLevel1[icv] = inbocv_v_mgLevel1;
		nbocv_i_mgLevel1[icv+1] = nbocv_i_mgLevel1[icv] + nNeighbor[icv] + 1;
		nbocv_v_mgLevel1[inbocv_v_mgLevel1] = icv; inbocv_v_mgLevel1++;
		for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++) {
			nbocv_v_mgLevel1[inbocv_v_mgLevel1] = neighbor[icv][iNeighbor]; inbocv_v_mgLevel1++;
		}
	}	
	nbocv_s_mgLevel1 = nbocv_i_mgLevel1[ncv_mgLevel1];
	
	
	
	// We are doing the  dimensionalization using the fine grid information
	assert(cvofa_mgLevel1==NULL);
	cvofa_mgLevel1 = new int[nfa][2];
	
	// First boundary faces using the cvofa vector
	bool check_coarse[ncv_mgLevel1];
	
	nfa_mgLevel1 = 0;
	for (zone_fine = faZoneList.begin(); zone_fine != faZoneList.end(); zone_fine++)
		if (zone_fine->getKind() == FA_ZONE_BOUNDARY) {
			faZoneList_mgLevel1.push_back(FaZone());
			FaZone * zone_coarse = &(faZoneList_mgLevel1.back());
			zone_coarse->setIndex(zone_fine->getIndex());
			zone_coarse->setName(zone_fine->getName());
			zone_coarse->setKind(zone_fine->getKind());
			zone_coarse->ifa_f = nfa_mgLevel1;
			for (int icv = 0; icv < ncv_mgLevel1; icv++) 
				check_coarse[icv] = true;
			for (int ifa = zone_fine->ifa_f; ifa<=zone_fine->ifa_l; ifa++) {
				int icv_fine =  cvofa[ifa][0];
				int icv_coarse = Parent_CV[icv_fine];
				// The control volume has not been previously added, and it is a real one (no ghost)
				if ((check_coarse[icv_coarse]) && (icv_coarse < ncv_mgLevel1)) {
					cvofa_mgLevel1[nfa_mgLevel1][0] = icv_coarse;
					cvofa_mgLevel1[nfa_mgLevel1][1] = -1;
					nfa_mgLevel1++;
					check_coarse[icv_coarse] = false;
					if ((zone_coarse->getNameString() == "symm") && (twoDimwithSymmetry == "YES")){
						cvofa_mgLevel1[nfa_mgLevel1][0] = icv_coarse;
						cvofa_mgLevel1[nfa_mgLevel1][1] = -1;
						nfa_mgLevel1++;
					}
				}
			}
			zone_coarse->ifa_l = nfa_mgLevel1-1;
		}
	
	nfa_b_mgLevel1 = nfa_bp_mgLevel1 = nfa_bpi_mgLevel1 = nfa_mgLevel1;
	
	// Then inner faces using the neighbor information (the sortering is from lower to higher values),
	// We don't add ghost control volumen at this stage. 
	for (int icv = 0; icv < ncv_mgLevel1; icv++)
		for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
			if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] < ncv_mgLevel1)) {
				cvofa_mgLevel1[nfa_mgLevel1][0] = icv;
				cvofa_mgLevel1[nfa_mgLevel1][1] = neighbor[icv][iNeighbor];
				nfa_mgLevel1++;
			}
	
	
	
	// Create the face of control volume structure using the buildFaocv() code
	assert(faocv_i_mgLevel1 == NULL);
	assert(faocv_v_mgLevel1 == NULL);
	
	// we don't assume anything about ordering here, simply that if 
	// cvofa[ifa][0,1] is >= 0, then it represents a valid connection...
	
	faocv_i_mgLevel1 = new int[ncv_mgLevel1+1];
	for (int icv = 0; icv < ncv_mgLevel1; icv++)
		faocv_i_mgLevel1[icv+1] = 0;
	
	for (int iter = 0; iter<2; iter++) {
		for (int ifa = 0; ifa<nfa_mgLevel1; ifa++) {
			// icv0 is always valid... 
			int icv = cvofa_mgLevel1[ifa][0];
			assert((icv>=0)&&(icv<ncv_mgLevel1));
			if (iter==0) { faocv_i_mgLevel1[icv+1] += 1; }
			else {
				faocv_v_mgLevel1[faocv_i_mgLevel1[icv]] = ifa;
				faocv_i_mgLevel1[icv] += 1;
			}
			// icv1 may or may not be valid - here we assume that
			// the storage of a periodic face index has been removed
			// and -1 indicates a boundary icv...
			icv = cvofa_mgLevel1[ifa][1];
			if (icv>=0) {
				assert(icv<ncv_mgLevel1);
				if (iter==0) { faocv_i_mgLevel1[icv+1] += 1; }
				else {
					faocv_v_mgLevel1[faocv_i_mgLevel1[icv]] = ifa;
					faocv_i_mgLevel1[icv] += 1;
				}
			}
		}
		if (iter==0) {
			faocv_i_mgLevel1[0] = 0;
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				faocv_i_mgLevel1[icv+1] += faocv_i_mgLevel1[icv];
			faocv_s_mgLevel1 = faocv_i_mgLevel1[ncv_mgLevel1];
			faocv_v_mgLevel1 = new int[faocv_s_mgLevel1];
		}
		else {
			for (int icv = ncv_mgLevel1; icv>0; icv--)
				faocv_i_mgLevel1[icv] = faocv_i_mgLevel1[icv-1];
			faocv_i_mgLevel1[0] = 0;
		}
	}
	
	int nfa_r_mgLevel1 = nfa_mgLevel1;
	
	// Add the ghost control volume to the cvofa structure
	for (int icv = 0; icv < ncv_mgLevel1; icv++)
		for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
			if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] >= ncv_mgLevel1)) {
				cvofa_mgLevel1[nfa_mgLevel1][0] = icv;
				cvofa_mgLevel1[nfa_mgLevel1][1] = neighbor[icv][iNeighbor];
				nfa_mgLevel1++;
			}
	
	// Set inner normal faces, and cv_volume		
	assert(fa_normal_mgLevel1==NULL);
	assert(x_fa_mgLevel1==NULL);
	fa_normal_mgLevel1 = new double[nfa_mgLevel1][3];
	x_fa_mgLevel1 = new double[nfa_mgLevel1][3];
	
	assert(cv_volume_mgLevel1==NULL);
	assert(x_cv_mgLevel1==NULL);
	cv_volume_mgLevel1 = new double[ncv_g_mgLevel1];
	x_cv_mgLevel1 = new double[ncv_g_mgLevel1][3];
	
	int fa_counter[nfa_mgLevel1];
	
	// Compute the area of the coarse volume
	for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel1; icv_coarse ++) {
		cv_volume_mgLevel1[icv_coarse] = 0.0;
		for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel1[icv_coarse]; iChildren ++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			cv_volume_mgLevel1[icv_coarse] += cv_volume[icv_fine];
		}
	}
	
	// Compute the cog of the coarse volume
	for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel1; icv_coarse ++) {
		x_cv_mgLevel1[icv_coarse][0] = 0.0;
		x_cv_mgLevel1[icv_coarse][1] = 0.0;
		x_cv_mgLevel1[icv_coarse][2] = 0.0;
		for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel1[icv_coarse]; iChildren ++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			x_cv_mgLevel1[icv_coarse][0] += x_cv[icv_fine][0];
			x_cv_mgLevel1[icv_coarse][1] += x_cv[icv_fine][1];
			x_cv_mgLevel1[icv_coarse][2] += x_cv[icv_fine][2];
		}
		x_cv_mgLevel1[icv_coarse][0] = x_cv_mgLevel1[icv_coarse][0]/nChildren_CV_mgLevel1[icv_coarse];
		x_cv_mgLevel1[icv_coarse][1] = x_cv_mgLevel1[icv_coarse][1]/nChildren_CV_mgLevel1[icv_coarse];
		x_cv_mgLevel1[icv_coarse][2] = x_cv_mgLevel1[icv_coarse][2]/nChildren_CV_mgLevel1[icv_coarse];
	}
	
	// Compute the normal face of the coarse grid for the internal faces
	for (int ifa = 0; ifa < nfa_mgLevel1; ifa ++) {
		fa_normal_mgLevel1[ifa][0] = 0.0;
		fa_normal_mgLevel1[ifa][1] = 0.0;
		fa_normal_mgLevel1[ifa][2] = 0.0;
		x_fa_mgLevel1[ifa][0] = 0.0;
		x_fa_mgLevel1[ifa][1] = 0.0;
		x_fa_mgLevel1[ifa][2] = 0.0;
		fa_counter[ifa] = 0;
	}
	
	bool *AddFace;
	AddFace = new bool[nfa];
	for (int ifa = 0; ifa < nfa; ifa++) AddFace[ifa] = true;
	
	for (int icv_coarse = 0; icv_coarse < ncv_mgLevel1; icv_coarse ++)
		for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel1[icv_coarse]; iChildren ++) {
			int icv_fine = Children_CV_mgLevel1[icv_coarse][iChildren];
			int nofinecv_f = nbocv_i[icv_fine];
			int nofinecv_l = nbocv_i[icv_fine + 1] - 1;
			for (int nocv = nofinecv_f+1; nocv <= nofinecv_l; nocv++) {
				int icv_fine_neighbor = nbocv_v[nocv];
				int icv_coarse_neighbor = Parent_CV[icv_fine_neighbor];
				
				if ((icv_fine < ncv) && (icv_fine_neighbor < ncv)) {
					
					// Find the face that conect the fine volumes (icv_fine and icv_fine_neighbor)
					int ifa_fine = findFace(icv_fine, icv_fine_neighbor);
					if (ifa_fine == -1) { cout << "There is an error detecting a fine face" << endl; cin.get(); }
					
					// We only add once each face of the finest grid
					if ((AddFace[ifa_fine]) && (icv_coarse != icv_coarse_neighbor)) {
						
						// Find the face that conect the coarse volumes (icv_coarse and icv_coarse_neighbor)...
						int ifa_coarse = findFace_mgLevel1(icv_coarse, icv_coarse_neighbor);
						if (ifa_coarse == -1) { cout << "There is an error detecting a coarse face" << icv_coarse <<", "<< icv_coarse_neighbor <<", "<< ifa_fine <<", "<< icv_fine_neighbor <<endl; cin.get(); }
						
						// Compute the face normal and check that the orientation is the same as in the finest grid
						double FaceOrientation = checkFaceOrientation(ifa_fine, icv_fine, icv_fine_neighbor);
						fa_normal_mgLevel1[ifa_coarse][0] += fa_normal[ifa_fine][0]*FaceOrientation;
						fa_normal_mgLevel1[ifa_coarse][1] += fa_normal[ifa_fine][1]*FaceOrientation;
						fa_normal_mgLevel1[ifa_coarse][2] += fa_normal[ifa_fine][2]*FaceOrientation;
						AddFace[ifa_fine] = false;
					}
				}
			}
		}
	
	
	// Compute the normal face of the coarse grid for the ghost faces
	for (int ifa_coarse = nfa_r_mgLevel1; ifa_coarse < nfa_mgLevel1; ifa_coarse++) {
		int real_icv_coarse = cvofa_mgLevel1[ifa_coarse][0];
		int ghost_icv_coarse = cvofa_mgLevel1[ifa_coarse][1];
		int ghost_icv_fine = Children_CV_mgLevel1[ghost_icv_coarse][0];
		for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel1[real_icv_coarse]; iChildren ++) {
			int real_icv_fine = Children_CV_mgLevel1[real_icv_coarse][iChildren];
			for (int ifa_fine = 0; ifa_fine < nfa; ifa_fine ++) {
				if (cvofa[ifa_fine][1] >= ncv) {
					if ((cvofa[ifa_fine][0] == real_icv_fine) && (cvofa[ifa_fine][1] == ghost_icv_fine)) {
						fa_normal_mgLevel1[ifa_coarse][0] += fa_normal[ifa_fine][0];
						fa_normal_mgLevel1[ifa_coarse][1] += fa_normal[ifa_fine][1];
						fa_normal_mgLevel1[ifa_coarse][2] += fa_normal[ifa_fine][2];
						
					}	
				}
			}
		}
	}
	
	// Compute the normal face of the coarse grid for the boundary faces
	int (*faceIndex)[10]; faceIndex = new int[ncv_mgLevel1][10];
	string (*faceName)[10]; faceName = new string[ncv_mgLevel1][10];
	int *nface; nface = new int[ncv_mgLevel1];
	int icv_coarse, ifa_coarse;
	
	for (int icv = 0; icv < ncv_mgLevel1; icv++) nface[icv] = 0;
	
	for (zone_coarse = faZoneList_mgLevel1.begin(); zone_coarse != faZoneList_mgLevel1.end(); zone_coarse++)
		if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
			for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
				icv_coarse = cvofa_mgLevel1[ifa_coarse][0];
				faceIndex[icv_coarse][nface[icv_coarse]] = ifa_coarse;
				faceName[icv_coarse][nface[icv_coarse]].assign(zone_coarse->getNameString());
				nface[icv_coarse]++;
			}	
	
	for (zone_fine = faZoneList.begin(); zone_fine != faZoneList.end(); zone_fine++)
		if (zone_fine->getKind() == FA_ZONE_BOUNDARY)
			for (int ifa_fine = zone_fine->ifa_f; ifa_fine<=zone_fine->ifa_l; ifa_fine++) {
				int icv_fine =  cvofa[ifa_fine][0];
				int icv_coarse = Parent_CV[icv_fine];
				for (int iface = 0; iface < nface[icv_coarse]; iface++) {
					if (faceName[icv_coarse][iface] == zone_fine->getNameString()) {
						ifa_coarse = faceIndex[icv_coarse][iface];
						goto Compute_face;
					}
				}
				
			Compute_face:
				if (AddFace[ifa_fine]) {					
					if ((zone_fine->getNameString() == "symm" ) && (twoDimwithSymmetry == "YES")) {
						if (fa_normal[ifa_fine][2] > 0.0) {
							fa_normal_mgLevel1[ifa_coarse][0] += fa_normal[ifa_fine][0];
							fa_normal_mgLevel1[ifa_coarse][1] += fa_normal[ifa_fine][1];
							fa_normal_mgLevel1[ifa_coarse][2] += fa_normal[ifa_fine][2];
							AddFace[ifa_fine] = false;	
						}
						else {
							fa_normal_mgLevel1[ifa_coarse+1][0] += fa_normal[ifa_fine][0];
							fa_normal_mgLevel1[ifa_coarse+1][1] += fa_normal[ifa_fine][1];
							fa_normal_mgLevel1[ifa_coarse+1][2] += fa_normal[ifa_fine][2];
							AddFace[ifa_fine] = false;	
						}
					}
					else {
						fa_normal_mgLevel1[ifa_coarse][0] += fa_normal[ifa_fine][0];
						fa_normal_mgLevel1[ifa_coarse][1] += fa_normal[ifa_fine][1];
						fa_normal_mgLevel1[ifa_coarse][2] += fa_normal[ifa_fine][2];
						AddFace[ifa_fine] = false;
					}
					if (zone_fine->getNameString() == "inlet" ) {
						x_fa_mgLevel1[ifa_coarse][0] += x_fa[ifa_fine][0];
						x_fa_mgLevel1[ifa_coarse][1] += x_fa[ifa_fine][1];
						x_fa_mgLevel1[ifa_coarse][2] += x_fa[ifa_fine][2];
						fa_counter[ifa_coarse]++;
					}
				}
			}
	
	for (zone_coarse = faZoneList_mgLevel1.begin(); zone_coarse != faZoneList_mgLevel1.end(); zone_coarse++)
		if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
			if (zone_coarse->getNameString() == "inlet" )
				for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
					x_fa_mgLevel1[ifa_coarse][0] = x_fa_mgLevel1[ifa_coarse][0]/fa_counter[ifa_coarse];
					x_fa_mgLevel1[ifa_coarse][1] = x_fa_mgLevel1[ifa_coarse][1]/fa_counter[ifa_coarse];
					x_fa_mgLevel1[ifa_coarse][2] = x_fa_mgLevel1[ifa_coarse][2]/fa_counter[ifa_coarse];
				}
	
	delete [] AddFace;
	
	MGSensor_mgLevel1 = new double[ncv_g_mgLevel1];
	
	rho_mgLevel1 = new double[ncv_g_mgLevel1];
  	rhou_mgLevel1 = new double[ncv_g_mgLevel1][3];
  	rhoE_mgLevel1 = new double[ncv_g_mgLevel1];
	
	rho_old_mgLevel1 = new double[ncv_g_mgLevel1];
  	rhou_old_mgLevel1 = new double[ncv_g_mgLevel1][3];
  	rhoE_old_mgLevel1 = new double[ncv_g_mgLevel1];
	
  	vel_mgLevel1 = new double[ncv_g_mgLevel1][3];
  	press_mgLevel1 = new double[ncv_g_mgLevel1];
  	temp_mgLevel1 = new double[ncv_g_mgLevel1];
  	enthalpy_mgLevel1 = new double[ncv_g_mgLevel1];
  	RoM_mgLevel1 = new double[ncv_g_mgLevel1];
	muT_mgLevel1 = new double[ncv_g_mgLevel1];
  	gamma_mgLevel1 = new double[ncv_g_mgLevel1];
  	sos_mgLevel1 = new double[ncv_g_mgLevel1];
	kine_mgLevel1 = new double[ncv_g_mgLevel1];
	local_dt_mgLevel1 = new double[ncv_g_mgLevel1];
	grad_enthalpy_mgLevel1 = new double[ncv_g_mgLevel1][3];
  	grad_u_mgLevel1 = new double[ncv_g_mgLevel1][3][3];
	muLam_mgLevel1 = new double[ncv_g_mgLevel1];
	LambdaOverCp_mgLevel1 = new double[ncv_g_mgLevel1];			
	
	rho_TruncError_mgLevel1 = new double[ncv_g_mgLevel1];
  	rhou_TruncError_mgLevel1 = new double[ncv_g_mgLevel1][3]; 
  	rhoE_TruncError_mgLevel1 = new double[ncv_g_mgLevel1];
	
	rho_res_old = new double[ncv_g];
  	rhou_res_old = new double[ncv_g][3]; 
  	rhoE_res_old = new double[ncv_g];
	
	rho_res_sum = new double[ncv_g];
  	rhou_res_sum = new double[ncv_g][3]; 
  	rhoE_res_sum = new double[ncv_g];
	
	rho_bfa_mgLevel1 = new double[nfa_mgLevel1];
	T_bfa_mgLevel1 = new double[nfa_mgLevel1];
	vel_bfa_mgLevel1 = new double[nfa_mgLevel1][3];
	p_bfa_mgLevel1 = new double[nfa_mgLevel1];
	h_bfa_mgLevel1 = new double[nfa_mgLevel1];
	gam_bfa_mgLevel1 = new double[nfa_mgLevel1];
	RoM_bfa_mgLevel1 = new double[nfa_mgLevel1];
	mul_fa_mgLevel1 = new double[nfa_mgLevel1];
	mut_fa_mgLevel1 = new double[nfa_mgLevel1];
	lamOcp_fa_mgLevel1 = new double[nfa_mgLevel1];
	
	d_interpI1 = new int[ncv_g];
	d_interpR1 = new double[ncv_g];
	d_interpR2 = new double[ncv_g][3];
	d_interpR3 = new double[ncv_g][3][3];
	
}	

void JoeWithModels::initializeFromMultiGrid(int iMesh) {
	string twoDimwithSymmetry = getStringParam("2D_WITH_SYMMETRY","NO");
	
	double my_volume_sum = 0.0;
	for (int icv = 0; icv<ncv; icv++)
		my_volume_sum += cv_volume[icv];
	double TotalVolume;
	MPI_Allreduce(&my_volume_sum, &TotalVolume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	
	
	int max_neighbor = 150;
	
	// Relation between the domain length and the aglomerated element.
	double Maxlength = getDoubleParam("MAX_LENGTH_AGGLOMERATION", "1.0E8");
	
	

	
	if (iMesh == 2) {
		
		list<FaZone>::iterator zone_fine;
		list<FaZone>::iterator zone_coarse;
		
		// Allocate some structures for doing the agglomeration
		Parent_CV_mgLevel1 = new int[ncv_g_mgLevel1];
		Agglomerate_mgLevel1 = new bool[ncv_g_mgLevel1];
		for (int icv = 0; icv < ncv_g_mgLevel1; icv++)
			Agglomerate_mgLevel1[icv] = false;
		
		// The ghost control volume can not be aglomerated using the information of this domain
		for (int icv = ncv_mgLevel1; icv < ncv_g_mgLevel1; icv++)
			Agglomerate_mgLevel1[icv] = true;
		
		nChildren_CV_mgLevel2 = new int[ncv_g_mgLevel1];
		Children_CV_mgLevel2 = new int *[ncv_g_mgLevel1];
		for (int icv = 0; icv < ncv_g_mgLevel1; icv++)
			Children_CV_mgLevel2[icv] = new int[max_neighbor];
		
		// Compute the agglomeration struture without compact storage method
		ncv_mgLevel2 = 0;
		
		
		// Second Agglomerate_mgLevel1 the inner elements
		for (int icv = 0; icv < ncv_mgLevel1; icv++)
			
			if (Agglomerate_mgLevel1[icv] == false) { 
				
				// Add the seed control volume
				Parent_CV_mgLevel1[icv] = ncv_mgLevel2; 
				Agglomerate_mgLevel1[icv] = true;
				nChildren_CV_mgLevel2[ncv_mgLevel2] = 1;
				Children_CV_mgLevel2[ncv_mgLevel2][0] = icv;
				
//				vector<int> Repeated_Neighbor_Points;
//				Repeated_Neighbor_Points.clear();
//				Repeated_Neighbor_Points.push_back(icv);
				
				// Add the first neighbor control volume
				int noc_f = nbocv_i_mgLevel1[icv];
				int noc_l = nbocv_i_mgLevel1[icv+1]-1;
				for (int noc = noc_f; noc <= noc_l; noc++) {
					int nbocv = nbocv_v_mgLevel1[noc];
					if (Agglomerate_mgLevel1[nbocv] == false) {
//						Repeated_Neighbor_Points.push_back(nbocv);
						Parent_CV_mgLevel1[nbocv] = ncv_mgLevel2;
						Agglomerate_mgLevel1[nbocv] = true;
						Children_CV_mgLevel2[ncv_mgLevel2][nChildren_CV_mgLevel2[ncv_mgLevel2]] = nbocv;
						nChildren_CV_mgLevel2[ncv_mgLevel2]++;
					}
				}
				
				// Add the second neighbor control volume
				vector<unsigned long> First_Neighbor_Points;
				First_Neighbor_Points.push_back(icv);
				
				// Localize first neighbors and create a list
				for (int noc = nbocv_i_mgLevel1[icv]+1; noc <= nbocv_i_mgLevel1[icv+1]-1; noc++)
					First_Neighbor_Points.push_back(nbocv_v_mgLevel1[noc]);
				
				// Identify a list with the second neighbors and its origin
				vector<int> Second_Neighbor_Points;
				vector<int> Second_Origin_Points;
				
				for (int iNode = nbocv_i_mgLevel1[icv]+1; iNode <= nbocv_i_mgLevel1[icv+1]-1; iNode++) {
					int jPoint = nbocv_v_mgLevel1[iNode];
					if (jPoint < ncv_mgLevel1) {
						for (int jNode = nbocv_i_mgLevel1[jPoint]+1; jNode <= nbocv_i_mgLevel1[jPoint+1]-1; jNode++) {
							int kPoint = nbocv_v_mgLevel1[jNode];
							
							// Verify if the element belongs to the list of the first neighbors
							bool SecondNeighborSeed = true;
							for (int iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
								if (kPoint == First_Neighbor_Points[iNeighbor]) {
									SecondNeighborSeed = false;
									break;
								}
							if (SecondNeighborSeed) {
								Second_Neighbor_Points.push_back(kPoint);
								Second_Origin_Points.push_back(jPoint);
							}
						}
					}
				}
				
				// Control volume that is neighbor of two first neighbors is added
				vector<int> Repeated_Neighbor_Points;					
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
					for (int jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
						if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) && (iNeighbor != jNeighbor))
							Repeated_Neighbor_Points.push_back(Second_Neighbor_Points[iNeighbor]);
				
/*				// Second neighbor that is neighbor of another second neighbor is added too
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
					for (int iNode = nbocv_i_mgLevel1[Second_Neighbor_Points[iNeighbor]]+1; 
						 iNode <= nbocv_i_mgLevel1[Second_Neighbor_Points[iNeighbor]+1]-1; iNode++) {
						int Third_Neighbor = nbocv_v_mgLevel1[iNode];
						
						for (int jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++) {
							if ((Second_Neighbor_Points[jNeighbor] == Third_Neighbor) && 
								(Second_Neighbor_Points[jNeighbor] != Second_Neighbor_Points[iNeighbor]))
								Repeated_Neighbor_Points.push_back(Third_Neighbor);
						}
					}
				
				
				double TempVol = cv_volume_mgLevel1[icv];	
				
				for (int kNode = 0; kNode <	Repeated_Neighbor_Points.size(); kNode ++) {
					for (int iNode = 0; iNode <	Repeated_Neighbor_Points.size(); iNode ++) {
						int nbocv = Repeated_Neighbor_Points[iNode];
						
						// Check the order for adding new elements (we add elements that have a 
						// common face, with the lelements that have been previously added )
						bool check = false;
						int init;
						if (nChildren_CV_mgLevel2[ncv_mgLevel2] == 1) init = 0;
						else init = 1;
						for (int iChildren = init; iChildren < nChildren_CV_mgLevel2[ncv_mgLevel2]; iChildren++) {
							int TestCV = Children_CV_mgLevel2[ncv_mgLevel2][iChildren];
							for (int jNode = nbocv_i_mgLevel1[TestCV]+1; jNode <= nbocv_i_mgLevel1[TestCV+1]-1; jNode++)
								if (nbocv_v_mgLevel1[jNode] == nbocv) { check = true; break; }
						}
						
						
						if ((Agglomerate_mgLevel1[nbocv] == false) && (TotalVolume * Maxlength > TempVol) && (check)) {
							TempVol += cv_volume_mgLevel1[nbocv];
							Parent_CV_mgLevel1[nbocv] = ncv_mgLevel2;
							Agglomerate_mgLevel1[nbocv] = true;
							Children_CV_mgLevel2[ncv_mgLevel2][nChildren_CV_mgLevel2[ncv_mgLevel2]] = nbocv;
							nChildren_CV_mgLevel2[ncv_mgLevel2]++;
						}
					}
				}*/
				
				for (int iNode = 0; iNode <	Repeated_Neighbor_Points.size(); iNode ++) {
					int nbocv = Repeated_Neighbor_Points[iNode];
					if (Agglomerate_mgLevel1[nbocv] == false) {
						Parent_CV_mgLevel1[nbocv] = ncv_mgLevel2;
						Agglomerate_mgLevel1[nbocv] = true;
						Children_CV_mgLevel2[ncv_mgLevel2][nChildren_CV_mgLevel2[ncv_mgLevel2]] = nbocv;
						nChildren_CV_mgLevel2[ncv_mgLevel2]++;
						
					}
				}
				
				// Update number coarse grid control volumes
				ncv_mgLevel2++; 
			}
		
		
		// In case there is a control volume which has not being Agglomerate_mgLevel1d
		for (int icv = 0; icv < ncv_mgLevel1; icv++)
			if (Agglomerate_mgLevel1[icv] == false) {
				Parent_CV_mgLevel1[icv] = ncv_mgLevel2;
				Agglomerate_mgLevel1[icv] = true;
				nChildren_CV_mgLevel2[ncv_mgLevel2] = 1;
				Children_CV_mgLevel2[ncv_mgLevel2][0] = icv;
				ncv_mgLevel2++;
			}
		
		// Store the number of control volumes of the coarse grid
		ncv_g_mgLevel2 = ncv_mgLevel2;
		
		// Add the ghost control volume at the end of the agglomerated volumes
		for (int icv = ncv_mgLevel1; icv < ncv_g_mgLevel1; icv++) {
			Parent_CV_mgLevel1[icv] = ncv_g_mgLevel2;
			Agglomerate_mgLevel1[icv] = true;
			nChildren_CV_mgLevel2[ncv_g_mgLevel2] = 1;
			Children_CV_mgLevel2[ncv_g_mgLevel2][0] = icv;
			ncv_g_mgLevel2++;
		}
		
		double MeanChildren = 0.0;
		int MaxChildren = 0; int MinChildren = max_neighbor;
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			MeanChildren += nChildren_CV_mgLevel2[icv];
			MaxChildren = max(nChildren_CV_mgLevel2[icv], MaxChildren);
			MinChildren = min(nChildren_CV_mgLevel2[icv], MinChildren);
		}
		
		int my_rank; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
		cout << "Domain: " << my_rank << ". CV of the MG level: " << ncv_mgLevel2 << ". Ghost CV of the MG level: " 
		<< ncv_g_mgLevel2-ncv_mgLevel2 << ". Reduction factor: 1/" <<MeanChildren/double(ncv_mgLevel2)<< "." <<endl;
		MPI_Barrier( MPI_COMM_WORLD );
//		cout << "Domain: " << my_rank << ". MaxChildren: " << MaxChildren << ". MinChildren: " 
//		<< MinChildren << ". MeanChildren: " <<int(MeanChildren/double(ncv_mgLevel2))<< "." <<endl;
		
		// Compute the neighbors of a control volume without compact storage method
		int **neighbor, *nNeighbor;
		neighbor = new int* [ncv_g_mgLevel2];
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++)
			neighbor[icv] = new int[max_neighbor];
		nNeighbor = new int[ncv_g_mgLevel2];
		bool add_cv;
		
		for (int icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse ++) {
			nNeighbor[icv_coarse] = 0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel2[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				int nboc_f = nbocv_i_mgLevel1[icv_fine];
				int nboc_l = nbocv_i_mgLevel1[icv_fine + 1] - 1;
				for (int noc = nboc_f + 1; noc <= nboc_l; noc++) {
					int icv_fine_neighbor = nbocv_v_mgLevel1[noc];
					int icv_parent = Parent_CV_mgLevel1[icv_fine_neighbor];
					if (icv_parent != icv_coarse) {
						add_cv = true;
						for (int iNeighbor = 0; iNeighbor < nNeighbor[icv_coarse]; iNeighbor++)
							if (neighbor[icv_coarse][iNeighbor] == icv_parent) {
								add_cv = false; break; 
							}
						if (add_cv) {
							neighbor[icv_coarse][nNeighbor[icv_coarse]] = icv_parent;
							nNeighbor[icv_coarse]++;
						}
					}
				}
			}
		}
		
		// Write the neighbors struture using a compact storage method
		assert(nbocv_i_mgLevel2==NULL);
		assert(nbocv_v_mgLevel2==NULL);
		
		nbocv_i_mgLevel2 = new int[ncv_mgLevel2+1];
		int nnbocv_v_mgLevel2 = ncv_mgLevel2;
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			nbocv_i_mgLevel2[icv+1] = 1;
			nnbocv_v_mgLevel2 += nNeighbor[icv];
		}
		nbocv_v_mgLevel2 = new int[nnbocv_v_mgLevel2];
		
		int inbocv_v_mgLevel2 = 0;
		for (int icv = 0; icv < ncv_mgLevel2; icv++) {
			nbocv_i_mgLevel2[icv] = inbocv_v_mgLevel2;
			nbocv_i_mgLevel2[icv+1] = nbocv_i_mgLevel2[icv] + nNeighbor[icv] + 1;
			nbocv_v_mgLevel2[inbocv_v_mgLevel2] = icv; inbocv_v_mgLevel2++;
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++) {
				nbocv_v_mgLevel2[inbocv_v_mgLevel2] = neighbor[icv][iNeighbor]; inbocv_v_mgLevel2++;
			}
		}	
		nbocv_s_mgLevel2 = nbocv_i_mgLevel2[ncv_mgLevel2];
		
		
		// We are doing the  dimensionalization using the fine grid information
		assert(cvofa_mgLevel2==NULL);
		cvofa_mgLevel2 = new int[nfa][2];
		
		bool check_coarse[ncv_mgLevel2];
		
		
		nfa_mgLevel2 = 0;
		for (zone_fine = faZoneList_mgLevel1.begin(); zone_fine != faZoneList_mgLevel1.end(); zone_fine++)
			if (zone_fine->getKind() == FA_ZONE_BOUNDARY) {
				faZoneList_mgLevel2.push_back(FaZone());
				FaZone * zone_coarse = &(faZoneList_mgLevel2.back());
				zone_coarse->setIndex(zone_fine->getIndex());
				zone_coarse->setName(zone_fine->getName());
				zone_coarse->setKind(zone_fine->getKind());
				zone_coarse->ifa_f = nfa_mgLevel2;
				for (int icv = 0; icv < ncv_mgLevel2; icv++) 
					check_coarse[icv] = true;
				for (int ifa = zone_fine->ifa_f; ifa<=zone_fine->ifa_l; ifa++) {
					int icv_fine =  cvofa_mgLevel1[ifa][0];
					int icv_coarse = Parent_CV_mgLevel1[icv_fine];
					if ((check_coarse[icv_coarse]) && (icv_coarse < ncv_mgLevel2))  { // We are checking if the face has been previously added to a coarse element
						cvofa_mgLevel2[nfa_mgLevel2][0] = icv_coarse;
						cvofa_mgLevel2[nfa_mgLevel2][1] = -1;
						nfa_mgLevel2++;
						check_coarse[icv_coarse] = false;
						if ((zone_coarse->getNameString() == "symm" ) && (twoDimwithSymmetry == "YES")) {
							cvofa_mgLevel2[nfa_mgLevel2][0] = icv_coarse;
							cvofa_mgLevel2[nfa_mgLevel2][1] = -1;
							nfa_mgLevel2++;
						}
					}
				}
				zone_coarse->ifa_l = nfa_mgLevel2-1;
			}
		
		nfa_b_mgLevel2 = nfa_bp_mgLevel2 = nfa_bpi_mgLevel2 = nfa_mgLevel2;
		
		// Then inner faces using the neighbor information (the sortering is from lower to higher values)
		for (int icv = 0; icv < ncv_mgLevel2; icv++)
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
				if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] < ncv_mgLevel2)) {
					cvofa_mgLevel2[nfa_mgLevel2][0] = icv;
					cvofa_mgLevel2[nfa_mgLevel2][1] = neighbor[icv][iNeighbor];
					nfa_mgLevel2++;
				}
		
		// Create the face of control volume structure using the buildFaocv() code
		assert(faocv_i_mgLevel2 == NULL);
		assert(faocv_v_mgLevel2 == NULL);
		
		// we don't assume anything about ordering here, simply that if 
		// cvofa_mgLevel1[ifa][0,1] is >= 0, then it represents a valid connection...
		
		faocv_i_mgLevel2 = new int[ncv_mgLevel2+1];
		for (int icv = 0; icv < ncv_mgLevel2; icv++)
			faocv_i_mgLevel2[icv+1] = 0;
		
		for (int iter = 0; iter<2; iter++) {
			for (int ifa = 0; ifa<nfa_mgLevel2; ifa++) {
				// icv0 is always valid... 
				int icv = cvofa_mgLevel2[ifa][0];
				assert((icv>=0)&&(icv<ncv_mgLevel2));
				if (iter==0) { faocv_i_mgLevel2[icv+1] += 1; }
				else {
					faocv_v_mgLevel2[faocv_i_mgLevel2[icv]] = ifa;
					faocv_i_mgLevel2[icv] += 1;
				}
				// icv1 may or may not be valid - here we assume that
				// the storage of a periodic face index has been removed
				// and -1 indicates a boundary icv...
				icv = cvofa_mgLevel2[ifa][1];
				if (icv>=0) {
					assert(icv<ncv_mgLevel2);
					if (iter==0) { faocv_i_mgLevel2[icv+1] += 1; }
					else {
						faocv_v_mgLevel2[faocv_i_mgLevel2[icv]] = ifa;
						faocv_i_mgLevel2[icv] += 1;
					}
				}
			}
			if (iter==0) {
				faocv_i_mgLevel2[0] = 0;
				for (int icv = 0; icv<ncv_mgLevel2; icv++)
					faocv_i_mgLevel2[icv+1] += faocv_i_mgLevel2[icv];
				faocv_s_mgLevel2 = faocv_i_mgLevel2[ncv_mgLevel2];
				faocv_v_mgLevel2 = new int[faocv_s_mgLevel2];
			}
			else {
				for (int icv = ncv_mgLevel2; icv>0; icv--)
					faocv_i_mgLevel2[icv] = faocv_i_mgLevel2[icv-1];
				faocv_i_mgLevel2[0] = 0;
			}
		}
		
		int nfa_r_mgLevel2 = nfa_mgLevel2;
		
		// Add the ghost control volume to the cvofa structure
		for (int icv = 0; icv < ncv_mgLevel2; icv++)
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
				if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] >= ncv_mgLevel2)) {
					cvofa_mgLevel2[nfa_mgLevel2][0] = icv;
					cvofa_mgLevel2[nfa_mgLevel2][1] = neighbor[icv][iNeighbor];
					nfa_mgLevel2++;
				}
		
		
		// Set inner normal faces, and cv_volume		
		assert(fa_normal_mgLevel2==NULL);
		assert(x_fa_mgLevel2==NULL);
		fa_normal_mgLevel2 = new double[nfa_mgLevel2][3];
		x_fa_mgLevel2 = new double[nfa_mgLevel2][3];
		
		assert(cv_volume_mgLevel2==NULL);
		assert(x_cv_mgLevel2==NULL);
		cv_volume_mgLevel2 = new double[ncv_g_mgLevel2];
		x_cv_mgLevel2 = new double[ncv_g_mgLevel2][3];
		
		int fa_counter[nfa_mgLevel2];
		
		// Compute the area of the coarse volume
		for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel2; icv_coarse ++) {
			cv_volume_mgLevel2[icv_coarse] = 0.0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel2[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				cv_volume_mgLevel2[icv_coarse] += cv_volume_mgLevel1[icv_fine];
			}
		}
		
		// Compute the cog of the coarse volume
		for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel2; icv_coarse ++) {
			x_cv_mgLevel2[icv_coarse][0] = 0.0;
			x_cv_mgLevel2[icv_coarse][1] = 0.0;
			x_cv_mgLevel2[icv_coarse][2] = 0.0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel2[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				x_cv_mgLevel2[icv_coarse][0] += x_cv_mgLevel1[icv_fine][0];
				x_cv_mgLevel2[icv_coarse][1] += x_cv_mgLevel1[icv_fine][1];
				x_cv_mgLevel2[icv_coarse][2] += x_cv_mgLevel1[icv_fine][2];
			}
			x_cv_mgLevel2[icv_coarse][0] = x_cv_mgLevel2[icv_coarse][0]/nChildren_CV_mgLevel2[icv_coarse];
			x_cv_mgLevel2[icv_coarse][1] = x_cv_mgLevel2[icv_coarse][1]/nChildren_CV_mgLevel2[icv_coarse];
			x_cv_mgLevel2[icv_coarse][2] = x_cv_mgLevel2[icv_coarse][2]/nChildren_CV_mgLevel2[icv_coarse];
		}
		
		
		// Compute the normal face of the coarse grid for the internal faces
		for (int ifa = 0; ifa < nfa_mgLevel2; ifa ++) {
			fa_normal_mgLevel2[ifa][0] = 0.0;
			fa_normal_mgLevel2[ifa][1] = 0.0;
			fa_normal_mgLevel2[ifa][2] = 0.0;
			x_fa_mgLevel2[ifa][0] = 0.0;
			x_fa_mgLevel2[ifa][1] = 0.0;
			x_fa_mgLevel2[ifa][2] = 0.0;
			fa_counter[ifa] = 0;
		}
		
		bool *AddFace;
		AddFace = new bool[nfa_mgLevel1];
		for (int ifa = 0; ifa < nfa_mgLevel1; ifa++) AddFace[ifa] = true;
		
		for (int icv_coarse = 0; icv_coarse < ncv_mgLevel2; icv_coarse ++)
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel2[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel2[icv_coarse][iChildren];
				int nofinecv_f = nbocv_i_mgLevel1[icv_fine];
				int nofinecv_l = nbocv_i_mgLevel1[icv_fine + 1] - 1;
				for (int nocv = nofinecv_f+1; nocv <= nofinecv_l; nocv++) {
					int icv_fine_neighbor = nbocv_v_mgLevel1[nocv];
					int icv_coarse_neighbor = Parent_CV_mgLevel1[icv_fine_neighbor];
					
					if ((icv_fine < ncv_mgLevel1) && (icv_fine_neighbor < ncv_mgLevel1)) {
						
						// Find the face that conect the fine volumes (icv_fine and icv_fine_neighbor)
						int ifa_fine = findFace_mgLevel1(icv_fine, icv_fine_neighbor);
						if (ifa_fine == -1) { cout << "There is an error detecting a fine face" << endl; cin.get(); }
						
						// We only add once each face of the finest grid
						if ((AddFace[ifa_fine]) && (icv_coarse != icv_coarse_neighbor)) {
							
							// Find the face that conect the coarse volumes (icv_coarse and icv_coarse_neighbor)...
							int ifa_coarse = findFace_mgLevel2(icv_coarse, icv_coarse_neighbor);
							if (ifa_coarse == -1) { cout << "There is an error detecting a coarse face" << icv_coarse<<" "<< icv_coarse_neighbor <<" "<< icv_fine <<" "<< icv_fine_neighbor << endl; cin.get(); }
							
							// Compute the face normal and check that the orientation is the same as in the finest grid
							double FaceOrientation = checkFaceOrientation_mgLevel1(ifa_fine, icv_fine, icv_fine_neighbor);
							fa_normal_mgLevel2[ifa_coarse][0] += fa_normal_mgLevel1[ifa_fine][0]*FaceOrientation;
							fa_normal_mgLevel2[ifa_coarse][1] += fa_normal_mgLevel1[ifa_fine][1]*FaceOrientation;
							fa_normal_mgLevel2[ifa_coarse][2] += fa_normal_mgLevel1[ifa_fine][2]*FaceOrientation;
							AddFace[ifa_fine] = false;
						}
					}
				}
			}
		
		
		// Compute the normal face of the coarse grid for the ghost faces
		for (int ifa_coarse = nfa_r_mgLevel2; ifa_coarse < nfa_mgLevel2; ifa_coarse++) {
			int real_icv_coarse = cvofa_mgLevel2[ifa_coarse][0];
			int ghost_icv_coarse = cvofa_mgLevel2[ifa_coarse][1];
			int ghost_icv_fine = Children_CV_mgLevel2[ghost_icv_coarse][0];
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel2[real_icv_coarse]; iChildren ++) {
				int real_icv_fine = Children_CV_mgLevel2[real_icv_coarse][iChildren];
				for (int ifa_fine = 0; ifa_fine < nfa_mgLevel1; ifa_fine ++) {
					if (cvofa_mgLevel1[ifa_fine][1] >= ncv_mgLevel1) {
						if ((cvofa_mgLevel1[ifa_fine][0] == real_icv_fine) && (cvofa_mgLevel1[ifa_fine][1] == ghost_icv_fine)) {
							fa_normal_mgLevel2[ifa_coarse][0] += fa_normal_mgLevel1[ifa_fine][0];
							fa_normal_mgLevel2[ifa_coarse][1] += fa_normal_mgLevel1[ifa_fine][1];
							fa_normal_mgLevel2[ifa_coarse][2] += fa_normal_mgLevel1[ifa_fine][2];
						}	
					}
				}
			}
		}
		
		
		int (*faceIndex)[10]; faceIndex = new int[ncv_mgLevel2][10];
		string (*faceName)[10]; faceName = new string[ncv_mgLevel2][10];
		int *nface; nface = new int[ncv_mgLevel2];
		int icv_coarse, ifa_coarse;
		
		for (int icv = 0; icv < ncv_mgLevel2; icv++) nface[icv] = 0;
		
		for (zone_coarse = faZoneList_mgLevel2.begin(); zone_coarse != faZoneList_mgLevel2.end(); zone_coarse++)
			if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
				for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
					icv_coarse = cvofa_mgLevel2[ifa_coarse][0];
					faceIndex[icv_coarse][nface[icv_coarse]] = ifa_coarse;
					faceName[icv_coarse][nface[icv_coarse]].assign(zone_coarse->getNameString());
					nface[icv_coarse]++;
				}	
		
		for (zone_fine = faZoneList_mgLevel1.begin(); zone_fine != faZoneList_mgLevel1.end(); zone_fine++)
			if (zone_fine->getKind() == FA_ZONE_BOUNDARY)
				for (int ifa_fine = zone_fine->ifa_f; ifa_fine<=zone_fine->ifa_l; ifa_fine++) {
					int icv_fine =  cvofa_mgLevel1[ifa_fine][0];
					int icv_coarse = Parent_CV_mgLevel1[icv_fine];
					for (int iface = 0; iface < nface[icv_coarse]; iface++) {
						if (faceName[icv_coarse][iface] == zone_fine->getNameString()) {
							ifa_coarse = faceIndex[icv_coarse][iface];
							goto Compute_face;
						}
					}
					
				Compute_face:
					if (AddFace[ifa_fine]) {
						if ((zone_fine->getNameString() == "symm" ) && (twoDimwithSymmetry == "YES")) {
							if (fa_normal_mgLevel1[ifa_fine][2] > 0.0) {
								fa_normal_mgLevel2[ifa_coarse][0] += fa_normal_mgLevel1[ifa_fine][0];
								fa_normal_mgLevel2[ifa_coarse][1] += fa_normal_mgLevel1[ifa_fine][1];
								fa_normal_mgLevel2[ifa_coarse][2] += fa_normal_mgLevel1[ifa_fine][2];
								AddFace[ifa_fine] = false;	
							}
							else {
								fa_normal_mgLevel2[ifa_coarse+1][0] += fa_normal_mgLevel1[ifa_fine][0];
								fa_normal_mgLevel2[ifa_coarse+1][1] += fa_normal_mgLevel1[ifa_fine][1];
								fa_normal_mgLevel2[ifa_coarse+1][2] += fa_normal_mgLevel1[ifa_fine][2];
								AddFace[ifa_fine] = false;	
							}
						}
						else {
							fa_normal_mgLevel2[ifa_coarse][0] += fa_normal_mgLevel1[ifa_fine][0];
							fa_normal_mgLevel2[ifa_coarse][1] += fa_normal_mgLevel1[ifa_fine][1];
							fa_normal_mgLevel2[ifa_coarse][2] += fa_normal_mgLevel1[ifa_fine][2];
							AddFace[ifa_fine] = false;
						}
						if (zone_fine->getNameString() == "inlet" ) {
							x_fa_mgLevel2[ifa_coarse][0] += x_fa_mgLevel1[ifa_fine][0];
							x_fa_mgLevel2[ifa_coarse][1] += x_fa_mgLevel1[ifa_fine][1];
							x_fa_mgLevel2[ifa_coarse][2] += x_fa_mgLevel1[ifa_fine][2];
							fa_counter[ifa_coarse]++;
						}
					}
				}
		
		for (zone_coarse = faZoneList_mgLevel2.begin(); zone_coarse != faZoneList_mgLevel2.end(); zone_coarse++)
			if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
				if (zone_coarse->getNameString() == "inlet" )
					for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
						x_fa_mgLevel2[ifa_coarse][0] = x_fa_mgLevel2[ifa_coarse][0]/fa_counter[ifa_coarse];
						x_fa_mgLevel2[ifa_coarse][1] = x_fa_mgLevel2[ifa_coarse][1]/fa_counter[ifa_coarse];
						x_fa_mgLevel2[ifa_coarse][2] = x_fa_mgLevel2[ifa_coarse][2]/fa_counter[ifa_coarse];
					}
		
		
		delete [] AddFace;
		
		MGSensor_mgLevel2 = new double[ncv_g_mgLevel2];
		
		rho_mgLevel2 = new double[ncv_g_mgLevel2];
		rhou_mgLevel2 = new double[ncv_g_mgLevel2][3];
		rhoE_mgLevel2 = new double[ncv_g_mgLevel2];
		
		rho_old_mgLevel2 = new double[ncv_g_mgLevel2];
		rhou_old_mgLevel2 = new double[ncv_g_mgLevel2][3];
		rhoE_old_mgLevel2 = new double[ncv_g_mgLevel2];
		
		vel_mgLevel2 = new double[ncv_g_mgLevel2][3];
		press_mgLevel2 = new double[ncv_g_mgLevel2];
		temp_mgLevel2 = new double[ncv_g_mgLevel2];
		enthalpy_mgLevel2 = new double[ncv_g_mgLevel2];
		RoM_mgLevel2 = new double[ncv_g_mgLevel2];
		muT_mgLevel2 = new double[ncv_g_mgLevel2];
		gamma_mgLevel2 = new double[ncv_g_mgLevel2];
		sos_mgLevel2 = new double[ncv_g_mgLevel2];
		kine_mgLevel2 = new double[ncv_g_mgLevel2];
		local_dt_mgLevel2 = new double[ncv_g_mgLevel2];
		grad_enthalpy_mgLevel2 = new double[ncv_g_mgLevel2][3];
		grad_u_mgLevel2 = new double[ncv_g_mgLevel2][3][3];
		muLam_mgLevel2 = new double[ncv_g_mgLevel2];
		LambdaOverCp_mgLevel2 = new double[ncv_g_mgLevel2];			
		
		
		rho_TruncError_mgLevel2 = new double[ncv_g_mgLevel2];
		rhou_TruncError_mgLevel2 = new double[ncv_g_mgLevel2][3]; 
		rhoE_TruncError_mgLevel2 = new double[ncv_g_mgLevel2];
		
		rho_bfa_mgLevel2 = new double[nfa_mgLevel2];
		T_bfa_mgLevel2 = new double[nfa_mgLevel2];
		vel_bfa_mgLevel2 = new double[nfa_mgLevel2][3];
		p_bfa_mgLevel2 = new double[nfa_mgLevel2];
		h_bfa_mgLevel2 = new double[nfa_mgLevel2];
		gam_bfa_mgLevel2 = new double[nfa_mgLevel2];
		RoM_bfa_mgLevel2 = new double[nfa_mgLevel2];
		mul_fa_mgLevel2 = new double[nfa_mgLevel2];
		mut_fa_mgLevel2 = new double[nfa_mgLevel2];
		lamOcp_fa_mgLevel2 = new double[nfa_mgLevel2];
		
	}
	
	if (iMesh == 3) {
		
		list<FaZone>::iterator zone_fine;
		list<FaZone>::iterator zone_coarse;
		
		// Allocate some structures for doing the agglomeration
		Parent_CV_mgLevel2 = new int[ncv_g_mgLevel2];
		Agglomerate_mgLevel2 = new bool[ncv_g_mgLevel2];
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++)
			Agglomerate_mgLevel2[icv] = false;
		
		// The ghost control volume can not be aglomerated using the information of this domain
		for (int icv = ncv_mgLevel2; icv < ncv_g_mgLevel2; icv++)
			Agglomerate_mgLevel2[icv] = true;
		
		nChildren_CV_mgLevel3 = new int[ncv_g_mgLevel2];
		Children_CV_mgLevel3 = new int *[ncv_g_mgLevel2];
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++)
			Children_CV_mgLevel3[icv] = new int[max_neighbor];
		
		// Compute the agglomeration struture without compact storage method
		ncv_mgLevel3 = 0;
		
		
		// Second Agglomerate_mgLevel2 the inner elements
		for (int icv = 0; icv < ncv_mgLevel2; icv++)
			
			if (Agglomerate_mgLevel2[icv] == false) { 
				
				// Add the seed control volume
				Parent_CV_mgLevel2[icv] = ncv_mgLevel3; 
				Agglomerate_mgLevel2[icv] = true;
				nChildren_CV_mgLevel3[ncv_mgLevel3] = 1;
				Children_CV_mgLevel3[ncv_mgLevel3][0] = icv;
				
				vector<int> Repeated_Neighbor_Points;
				Repeated_Neighbor_Points.clear();
				Repeated_Neighbor_Points.push_back(icv);
				
				// Add the first neighbor control volume
				int noc_f = nbocv_i_mgLevel2[icv];
				int noc_l = nbocv_i_mgLevel2[icv+1]-1;
				for (int noc = noc_f; noc <= noc_l; noc++) {
					int nbocv = nbocv_v_mgLevel2[noc];
					if (Agglomerate_mgLevel2[nbocv] == false) {
						
						Repeated_Neighbor_Points.push_back(nbocv);
						Parent_CV_mgLevel2[nbocv] = ncv_mgLevel3;
						Agglomerate_mgLevel2[nbocv] = true;
						Children_CV_mgLevel3[ncv_mgLevel3][nChildren_CV_mgLevel3[ncv_mgLevel3]] = nbocv;
						nChildren_CV_mgLevel3[ncv_mgLevel3]++;
					}
				}
				
				// Add the second neighbor control volume
				vector<unsigned long> First_Neighbor_Points;
				First_Neighbor_Points.push_back(icv);
				
				// Localize first neighbors and create a list
				for (int noc = nbocv_i_mgLevel2[icv]+1; noc <= nbocv_i_mgLevel2[icv+1]-1; noc++)
					First_Neighbor_Points.push_back(nbocv_v_mgLevel2[noc]);
				
				// Identify a list with the second neighbors and its origin
				vector<int> Second_Neighbor_Points;
				vector<int> Second_Origin_Points;
				
				for (int iNode = nbocv_i_mgLevel2[icv]+1; iNode <= nbocv_i_mgLevel2[icv+1]-1; iNode++) {
					int jPoint = nbocv_v_mgLevel2[iNode];
					if (jPoint < ncv_mgLevel2) {
						for (int jNode = nbocv_i_mgLevel2[jPoint]+1; jNode <= nbocv_i_mgLevel2[jPoint+1]-1; jNode++) {
							int kPoint = nbocv_v_mgLevel2[jNode];
							
							// Verify if the element belongs to the list of the first neighbors
							bool SecondNeighborSeed = true;
							for (int iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
								if (kPoint == First_Neighbor_Points[iNeighbor]) {
									SecondNeighborSeed = false;
									break;
								}
							if (SecondNeighborSeed) {
								Second_Neighbor_Points.push_back(kPoint);
								Second_Origin_Points.push_back(jPoint);
							}
						}
					}
				}
				
				// Control volume that is neighbor of two first neighbors is added
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
					for (int jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++) {
						if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) && (iNeighbor != jNeighbor)) {
							Repeated_Neighbor_Points.push_back(Second_Neighbor_Points[iNeighbor]);
						}
					}
				
				// Second neighbor that is neighbor of another second neighbor is added too
				for (int iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
					for (int iNode = nbocv_i_mgLevel2[Second_Neighbor_Points[iNeighbor]]+1; 
						 iNode <= nbocv_i_mgLevel2[Second_Neighbor_Points[iNeighbor]+1]-1; iNode++) {
						int Third_Neighbor = nbocv_v_mgLevel2[iNode];
						
						for (int jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++) {
							if ((Second_Neighbor_Points[jNeighbor] == Third_Neighbor) && 
								(Second_Neighbor_Points[jNeighbor] != Second_Neighbor_Points[iNeighbor]))
								Repeated_Neighbor_Points.push_back(Third_Neighbor);
						}
					}
				
				
				double TempVol = cv_volume_mgLevel2[icv];	
				
				for (int kNode = 0; kNode <	Repeated_Neighbor_Points.size(); kNode ++) {
					for (int iNode = 0; iNode <	Repeated_Neighbor_Points.size(); iNode ++) {
						int nbocv = Repeated_Neighbor_Points[iNode];
						
						// Check the order for adding new elements (we add elements that have a 
						// common face, with the lelements that have been previously added )
						bool check = false;
						int init;
						if (nChildren_CV_mgLevel3[ncv_mgLevel3] == 1) init = 0;
						else init = 1;
						for (int iChildren = init; iChildren < nChildren_CV_mgLevel3[ncv_mgLevel3]; iChildren++) {
							int TestCV = Children_CV_mgLevel3[ncv_mgLevel3][iChildren];
							for (int jNode = nbocv_i_mgLevel2[TestCV]+1; jNode <= nbocv_i_mgLevel2[TestCV+1]-1; jNode++)
								if (nbocv_v_mgLevel2[jNode] == nbocv) { check = true; break; }
						}
						
						
						if ((Agglomerate_mgLevel2[nbocv] == false) && (TotalVolume * Maxlength > TempVol) && (check)) {
							TempVol += cv_volume_mgLevel2[nbocv];
							Parent_CV_mgLevel2[nbocv] = ncv_mgLevel3;
							Agglomerate_mgLevel2[nbocv] = true;
							Children_CV_mgLevel3[ncv_mgLevel3][nChildren_CV_mgLevel3[ncv_mgLevel3]] = nbocv;
							nChildren_CV_mgLevel3[ncv_mgLevel3]++;
						}
					}
				}
				
				// Update number coarse grid control volumes
				ncv_mgLevel3++; 
			}
		
		
		// In case there is a control volume which has not being Agglomerate_mgLevel2d
		for (int icv = 0; icv < ncv_mgLevel2; icv++)
			if (Agglomerate_mgLevel2[icv] == false) {
				Parent_CV_mgLevel2[icv] = ncv_mgLevel3;
				Agglomerate_mgLevel2[icv] = true;
				nChildren_CV_mgLevel3[ncv_mgLevel3] = 1;
				Children_CV_mgLevel3[ncv_mgLevel3][0] = icv;
				ncv_mgLevel3++;
			}
		
		// Store the number of control volumes of the coarse grid
		ncv_g_mgLevel3 = ncv_mgLevel3;
		
		// Add the ghost control volume at the end of the agglomerated volumes
		for (int icv = ncv_mgLevel2; icv < ncv_g_mgLevel2; icv++) {
			Parent_CV_mgLevel2[icv] = ncv_g_mgLevel3;
			Agglomerate_mgLevel2[icv] = true;
			nChildren_CV_mgLevel3[ncv_g_mgLevel3] = 1;
			Children_CV_mgLevel3[ncv_g_mgLevel3][0] = icv;
			ncv_g_mgLevel3++;
		}
		
		double MeanChildren = 0.0;
		int MaxChildren = 0; int MinChildren = max_neighbor;
		for (int icv = 0; icv < ncv_mgLevel3; icv++) {
			MeanChildren += nChildren_CV_mgLevel3[icv];
			MaxChildren = max(nChildren_CV_mgLevel3[icv], MaxChildren);
			MinChildren = min(nChildren_CV_mgLevel3[icv], MinChildren);
		}
		
		int my_rank; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
		cout << "Domain: " << my_rank << ". CV of the MG level: " << ncv_mgLevel3 << ". Ghost CV of the MG level: " 
		<< ncv_g_mgLevel3-ncv_mgLevel3 << ". Reduction factor: 1/" <<MeanChildren/double(ncv_mgLevel3)<< "." <<endl;
		MPI_Barrier( MPI_COMM_WORLD );
		//		cout << "Domain: " << my_rank << ". MaxChildren: " << MaxChildren << ". MinChildren: " 
		//		<< MinChildren << ". MeanChildren: " <<int(MeanChildren/double(ncv_mgLevel3))<< "." <<endl;
		
		// Compute the neighbors of a control volume without compact storage method
		int **neighbor, *nNeighbor;
		neighbor = new int* [ncv_g_mgLevel3];
		for (int icv = 0; icv < ncv_g_mgLevel3; icv++)
			neighbor[icv] = new int[max_neighbor];
		nNeighbor = new int[ncv_g_mgLevel3];
		bool add_cv;
		
		for (int icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse ++) {
			nNeighbor[icv_coarse] = 0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel3[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				int nboc_f = nbocv_i_mgLevel2[icv_fine];
				int nboc_l = nbocv_i_mgLevel2[icv_fine + 1] - 1;
				for (int noc = nboc_f + 1; noc <= nboc_l; noc++) {
					int icv_fine_neighbor = nbocv_v_mgLevel2[noc];
					int icv_parent = Parent_CV_mgLevel2[icv_fine_neighbor];
					if (icv_parent != icv_coarse) {
						add_cv = true;
						for (int iNeighbor = 0; iNeighbor < nNeighbor[icv_coarse]; iNeighbor++)
							if (neighbor[icv_coarse][iNeighbor] == icv_parent) {
								add_cv = false; break; 
							}
						if (add_cv) {
							neighbor[icv_coarse][nNeighbor[icv_coarse]] = icv_parent;
							nNeighbor[icv_coarse]++;
						}
					}
				}
			}
		}
		
		// Write the neighbors struture using a compact storage method
		assert(nbocv_i_mgLevel3==NULL);
		assert(nbocv_v_mgLevel3==NULL);
		
		nbocv_i_mgLevel3 = new int[ncv_mgLevel3+1];
		int nnbocv_v_mgLevel3 = ncv_mgLevel3;
		for (int icv = 0; icv < ncv_mgLevel3; icv++) {
			nbocv_i_mgLevel3[icv+1] = 1;
			nnbocv_v_mgLevel3 += nNeighbor[icv];
		}
		nbocv_v_mgLevel3 = new int[nnbocv_v_mgLevel3];
		
		int inbocv_v_mgLevel3 = 0;
		for (int icv = 0; icv < ncv_mgLevel3; icv++) {
			nbocv_i_mgLevel3[icv] = inbocv_v_mgLevel3;
			nbocv_i_mgLevel3[icv+1] = nbocv_i_mgLevel3[icv] + nNeighbor[icv] + 1;
			nbocv_v_mgLevel3[inbocv_v_mgLevel3] = icv; inbocv_v_mgLevel3++;
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++) {
				nbocv_v_mgLevel3[inbocv_v_mgLevel3] = neighbor[icv][iNeighbor]; inbocv_v_mgLevel3++;
			}
		}	
		nbocv_s_mgLevel3 = nbocv_i_mgLevel3[ncv_mgLevel3];
		
		
		// We are doing the  dimensionalization using the fine grid information
		assert(cvofa_mgLevel3==NULL);
		cvofa_mgLevel3 = new int[nfa][2];
		
		bool check_coarse[ncv_mgLevel3];
		
		
		nfa_mgLevel3 = 0;
		for (zone_fine = faZoneList_mgLevel2.begin(); zone_fine != faZoneList_mgLevel2.end(); zone_fine++)
			if (zone_fine->getKind() == FA_ZONE_BOUNDARY) {
				faZoneList_mgLevel3.push_back(FaZone());
				FaZone * zone_coarse = &(faZoneList_mgLevel3.back());
				zone_coarse->setIndex(zone_fine->getIndex());
				zone_coarse->setName(zone_fine->getName());
				zone_coarse->setKind(zone_fine->getKind());
				zone_coarse->ifa_f = nfa_mgLevel3;
				for (int icv = 0; icv < ncv_mgLevel3; icv++) 
					check_coarse[icv] = true;
				for (int ifa = zone_fine->ifa_f; ifa<=zone_fine->ifa_l; ifa++) {
					int icv_fine =  cvofa_mgLevel2[ifa][0];
					int icv_coarse = Parent_CV_mgLevel2[icv_fine];
					if ((check_coarse[icv_coarse]) && (icv_coarse < ncv_mgLevel3))  { // We are checking if the face has been previously added to a coarse element
						cvofa_mgLevel3[nfa_mgLevel3][0] = icv_coarse;
						cvofa_mgLevel3[nfa_mgLevel3][1] = -1;
						nfa_mgLevel3++;
						check_coarse[icv_coarse] = false;
						if ((zone_coarse->getNameString() == "symm" ) && (twoDimwithSymmetry == "YES")) {
							cvofa_mgLevel3[nfa_mgLevel3][0] = icv_coarse;
							cvofa_mgLevel3[nfa_mgLevel3][1] = -1;
							nfa_mgLevel3++;
						}
					}
				}
				zone_coarse->ifa_l = nfa_mgLevel3-1;
			}
		
		nfa_b_mgLevel3 = nfa_bp_mgLevel3 = nfa_bpi_mgLevel3 = nfa_mgLevel3;
		
		// Then inner faces using the neighbor information (the sortering is from lower to higher values)
		for (int icv = 0; icv < ncv_mgLevel3; icv++)
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
				if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] < ncv_mgLevel3)) {
					cvofa_mgLevel3[nfa_mgLevel3][0] = icv;
					cvofa_mgLevel3[nfa_mgLevel3][1] = neighbor[icv][iNeighbor];
					nfa_mgLevel3++;
				}
		
		// Create the face of control volume structure using the buildFaocv() code
		assert(faocv_i_mgLevel3 == NULL);
		assert(faocv_v_mgLevel3 == NULL);
		
		// we don't assume anything about ordering here, simply that if 
		// cvofa_mgLevel2[ifa][0,1] is >= 0, then it represents a valid connection...
		
		faocv_i_mgLevel3 = new int[ncv_mgLevel3+1];
		for (int icv = 0; icv < ncv_mgLevel3; icv++)
			faocv_i_mgLevel3[icv+1] = 0;
		
		for (int iter = 0; iter<2; iter++) {
			for (int ifa = 0; ifa<nfa_mgLevel3; ifa++) {
				// icv0 is always valid... 
				int icv = cvofa_mgLevel3[ifa][0];
				assert((icv>=0)&&(icv<ncv_mgLevel3));
				if (iter==0) { faocv_i_mgLevel3[icv+1] += 1; }
				else {
					faocv_v_mgLevel3[faocv_i_mgLevel3[icv]] = ifa;
					faocv_i_mgLevel3[icv] += 1;
				}
				// icv1 may or may not be valid - here we assume that
				// the storage of a periodic face index has been removed
				// and -1 indicates a boundary icv...
				icv = cvofa_mgLevel3[ifa][1];
				if (icv>=0) {
					assert(icv<ncv_mgLevel3);
					if (iter==0) { faocv_i_mgLevel3[icv+1] += 1; }
					else {
						faocv_v_mgLevel3[faocv_i_mgLevel3[icv]] = ifa;
						faocv_i_mgLevel3[icv] += 1;
					}
				}
			}
			if (iter==0) {
				faocv_i_mgLevel3[0] = 0;
				for (int icv = 0; icv<ncv_mgLevel3; icv++)
					faocv_i_mgLevel3[icv+1] += faocv_i_mgLevel3[icv];
				faocv_s_mgLevel3 = faocv_i_mgLevel3[ncv_mgLevel3];
				faocv_v_mgLevel3 = new int[faocv_s_mgLevel3];
			}
			else {
				for (int icv = ncv_mgLevel3; icv>0; icv--)
					faocv_i_mgLevel3[icv] = faocv_i_mgLevel3[icv-1];
				faocv_i_mgLevel3[0] = 0;
			}
		}
		
		int nfa_r_mgLevel3 = nfa_mgLevel3;
		
		// Add the ghost control volume to the cvofa structure
		for (int icv = 0; icv < ncv_mgLevel3; icv++)
			for (int iNeighbor = 0; iNeighbor < nNeighbor[icv]; iNeighbor++)
				if ((neighbor[icv][iNeighbor] > icv) && (neighbor[icv][iNeighbor] >= ncv_mgLevel3)) {
					cvofa_mgLevel3[nfa_mgLevel3][0] = icv;
					cvofa_mgLevel3[nfa_mgLevel3][1] = neighbor[icv][iNeighbor];
					nfa_mgLevel3++;
				}
		
		
		// Set inner normal faces, and cv_volume		
		assert(fa_normal_mgLevel3==NULL);
		assert(x_fa_mgLevel3==NULL);
		fa_normal_mgLevel3 = new double[nfa_mgLevel3][3];
		x_fa_mgLevel3 = new double[nfa_mgLevel3][3];
		
		assert(cv_volume_mgLevel3==NULL);
		assert(x_cv_mgLevel3==NULL);
		cv_volume_mgLevel3 = new double[ncv_g_mgLevel3];
		x_cv_mgLevel3 = new double[ncv_g_mgLevel3][3];
		
		int fa_counter[nfa_mgLevel3];
		
		// Compute the area of the coarse volume
		for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel3; icv_coarse ++) {
			cv_volume_mgLevel3[icv_coarse] = 0.0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel3[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				cv_volume_mgLevel3[icv_coarse] += cv_volume_mgLevel2[icv_fine];
			}
		}
		
		// Compute the cog of the coarse volume
		for (int icv_coarse = 0; icv_coarse < ncv_g_mgLevel3; icv_coarse ++) {
			x_cv_mgLevel3[icv_coarse][0] = 0.0;
			x_cv_mgLevel3[icv_coarse][1] = 0.0;
			x_cv_mgLevel3[icv_coarse][2] = 0.0;
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel3[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				x_cv_mgLevel3[icv_coarse][0] += x_cv_mgLevel2[icv_fine][0];
				x_cv_mgLevel3[icv_coarse][1] += x_cv_mgLevel2[icv_fine][1];
				x_cv_mgLevel3[icv_coarse][2] += x_cv_mgLevel2[icv_fine][2];
			}
			x_cv_mgLevel3[icv_coarse][0] = x_cv_mgLevel3[icv_coarse][0]/nChildren_CV_mgLevel3[icv_coarse];
			x_cv_mgLevel3[icv_coarse][1] = x_cv_mgLevel3[icv_coarse][1]/nChildren_CV_mgLevel3[icv_coarse];
			x_cv_mgLevel3[icv_coarse][2] = x_cv_mgLevel3[icv_coarse][2]/nChildren_CV_mgLevel3[icv_coarse];
		}
		
		
		// Compute the normal face of the coarse grid for the internal faces
		for (int ifa = 0; ifa < nfa_mgLevel3; ifa ++) {
			fa_normal_mgLevel3[ifa][0] = 0.0;
			fa_normal_mgLevel3[ifa][1] = 0.0;
			fa_normal_mgLevel3[ifa][2] = 0.0;
			x_fa_mgLevel3[ifa][0] = 0.0;
			x_fa_mgLevel3[ifa][1] = 0.0;
			x_fa_mgLevel3[ifa][2] = 0.0;
			fa_counter[ifa] = 0;
		}
		
		bool *AddFace;
		AddFace = new bool[nfa_mgLevel2];
		for (int ifa = 0; ifa < nfa_mgLevel2; ifa++) AddFace[ifa] = true;
		
		for (int icv_coarse = 0; icv_coarse < ncv_mgLevel3; icv_coarse ++)
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel3[icv_coarse]; iChildren ++) {
				int icv_fine = Children_CV_mgLevel3[icv_coarse][iChildren];
				int nofinecv_f = nbocv_i_mgLevel2[icv_fine];
				int nofinecv_l = nbocv_i_mgLevel2[icv_fine + 1] - 1;
				for (int nocv = nofinecv_f+1; nocv <= nofinecv_l; nocv++) {
					int icv_fine_neighbor = nbocv_v_mgLevel2[nocv];
					int icv_coarse_neighbor = Parent_CV_mgLevel2[icv_fine_neighbor];
					
					if ((icv_fine < ncv_mgLevel2) && (icv_fine_neighbor < ncv_mgLevel2)) {
						
						// Find the face that conect the fine volumes (icv_fine and icv_fine_neighbor)
						int ifa_fine = findFace_mgLevel2(icv_fine, icv_fine_neighbor);
						if (ifa_fine == -1) { cout << "There is an error detecting a fine face" << endl; cin.get(); }
						
						// We only add once each face of the finest grid
						if ((AddFace[ifa_fine]) && (icv_coarse != icv_coarse_neighbor)) {
							
							// Find the face that conect the coarse volumes (icv_coarse and icv_coarse_neighbor)...
							int ifa_coarse = findFace_mgLevel3(icv_coarse, icv_coarse_neighbor);
							if (ifa_coarse == -1) { cout << "There is an error detecting a coarse face" << icv_coarse<<" "<< icv_coarse_neighbor <<" "<< icv_fine <<" "<< icv_fine_neighbor << endl; cin.get(); }
							
							// Compute the face normal and check that the orientation is the same as in the finest grid
							double FaceOrientation = checkFaceOrientation_mgLevel2(ifa_fine, icv_fine, icv_fine_neighbor);
							fa_normal_mgLevel3[ifa_coarse][0] += fa_normal_mgLevel2[ifa_fine][0]*FaceOrientation;
							fa_normal_mgLevel3[ifa_coarse][1] += fa_normal_mgLevel2[ifa_fine][1]*FaceOrientation;
							fa_normal_mgLevel3[ifa_coarse][2] += fa_normal_mgLevel2[ifa_fine][2]*FaceOrientation;
							AddFace[ifa_fine] = false;
						}
					}
				}
			}
		
		
		// Compute the normal face of the coarse grid for the ghost faces
		for (int ifa_coarse = nfa_r_mgLevel3; ifa_coarse < nfa_mgLevel3; ifa_coarse++) {
			int real_icv_coarse = cvofa_mgLevel3[ifa_coarse][0];
			int ghost_icv_coarse = cvofa_mgLevel3[ifa_coarse][1];
			int ghost_icv_fine = Children_CV_mgLevel3[ghost_icv_coarse][0];
			for (int iChildren = 0; iChildren <  nChildren_CV_mgLevel3[real_icv_coarse]; iChildren ++) {
				int real_icv_fine = Children_CV_mgLevel3[real_icv_coarse][iChildren];
				for (int ifa_fine = 0; ifa_fine < nfa_mgLevel2; ifa_fine ++) {
					if (cvofa_mgLevel2[ifa_fine][1] >= ncv_mgLevel2) {
						if ((cvofa_mgLevel2[ifa_fine][0] == real_icv_fine) && (cvofa_mgLevel2[ifa_fine][1] == ghost_icv_fine)) {
							fa_normal_mgLevel3[ifa_coarse][0] += fa_normal_mgLevel2[ifa_fine][0];
							fa_normal_mgLevel3[ifa_coarse][1] += fa_normal_mgLevel2[ifa_fine][1];
							fa_normal_mgLevel3[ifa_coarse][2] += fa_normal_mgLevel2[ifa_fine][2];
						}	
					}
				}
			}
		}
		
		
		int (*faceIndex)[10]; faceIndex = new int[ncv_mgLevel3][10];
		string (*faceName)[10]; faceName = new string[ncv_mgLevel3][10];
		int *nface; nface = new int[ncv_mgLevel3];
		int icv_coarse, ifa_coarse;
		
		for (int icv = 0; icv < ncv_mgLevel3; icv++) nface[icv] = 0;
		
		for (zone_coarse = faZoneList_mgLevel3.begin(); zone_coarse != faZoneList_mgLevel3.end(); zone_coarse++)
			if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
				for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
					icv_coarse = cvofa_mgLevel3[ifa_coarse][0];
					faceIndex[icv_coarse][nface[icv_coarse]] = ifa_coarse;
					faceName[icv_coarse][nface[icv_coarse]].assign(zone_coarse->getNameString());
					nface[icv_coarse]++;
				}	
		
		for (zone_fine = faZoneList_mgLevel2.begin(); zone_fine != faZoneList_mgLevel2.end(); zone_fine++)
			if (zone_fine->getKind() == FA_ZONE_BOUNDARY)
				for (int ifa_fine = zone_fine->ifa_f; ifa_fine<=zone_fine->ifa_l; ifa_fine++) {
					int icv_fine =  cvofa_mgLevel2[ifa_fine][0];
					int icv_coarse = Parent_CV_mgLevel2[icv_fine];
					for (int iface = 0; iface < nface[icv_coarse]; iface++) {
						if (faceName[icv_coarse][iface] == zone_fine->getNameString()) {
							ifa_coarse = faceIndex[icv_coarse][iface];
							goto Compute_face_;
						}
					}
					
				Compute_face_:
					if (AddFace[ifa_fine]) {
						if ((zone_fine->getNameString() == "symm" ) && (twoDimwithSymmetry == "YES")) {
							if (fa_normal_mgLevel2[ifa_fine][2] > 0.0) {
								fa_normal_mgLevel3[ifa_coarse][0] += fa_normal_mgLevel2[ifa_fine][0];
								fa_normal_mgLevel3[ifa_coarse][1] += fa_normal_mgLevel2[ifa_fine][1];
								fa_normal_mgLevel3[ifa_coarse][2] += fa_normal_mgLevel2[ifa_fine][2];
								AddFace[ifa_fine] = false;	
							}
							else {
								fa_normal_mgLevel3[ifa_coarse+1][0] += fa_normal_mgLevel2[ifa_fine][0];
								fa_normal_mgLevel3[ifa_coarse+1][1] += fa_normal_mgLevel2[ifa_fine][1];
								fa_normal_mgLevel3[ifa_coarse+1][2] += fa_normal_mgLevel2[ifa_fine][2];
								AddFace[ifa_fine] = false;	
							}
						}
						else {
							fa_normal_mgLevel3[ifa_coarse][0] += fa_normal_mgLevel2[ifa_fine][0];
							fa_normal_mgLevel3[ifa_coarse][1] += fa_normal_mgLevel2[ifa_fine][1];
							fa_normal_mgLevel3[ifa_coarse][2] += fa_normal_mgLevel2[ifa_fine][2];
							AddFace[ifa_fine] = false;
						}
						if (zone_fine->getNameString() == "inlet" ) {
							x_fa_mgLevel3[ifa_coarse][0] += x_fa_mgLevel2[ifa_fine][0];
							x_fa_mgLevel3[ifa_coarse][1] += x_fa_mgLevel2[ifa_fine][1];
							x_fa_mgLevel3[ifa_coarse][2] += x_fa_mgLevel2[ifa_fine][2];
							fa_counter[ifa_coarse]++;
						}
					}
				}
		
		for (zone_coarse = faZoneList_mgLevel3.begin(); zone_coarse != faZoneList_mgLevel3.end(); zone_coarse++)
			if (zone_coarse->getKind() == FA_ZONE_BOUNDARY)
				if (zone_coarse->getNameString() == "inlet" )
					for (int ifa_coarse = zone_coarse->ifa_f; ifa_coarse<=zone_coarse->ifa_l; ifa_coarse++) {
						x_fa_mgLevel3[ifa_coarse][0] = x_fa_mgLevel3[ifa_coarse][0]/fa_counter[ifa_coarse];
						x_fa_mgLevel3[ifa_coarse][1] = x_fa_mgLevel3[ifa_coarse][1]/fa_counter[ifa_coarse];
						x_fa_mgLevel3[ifa_coarse][2] = x_fa_mgLevel3[ifa_coarse][2]/fa_counter[ifa_coarse];
					}
		
		
		delete [] AddFace;
		
		MGSensor_mgLevel3 = new double[ncv_g_mgLevel3];
		
		rho_mgLevel3 = new double[ncv_g_mgLevel3];
		rhou_mgLevel3 = new double[ncv_g_mgLevel3][3];
		rhoE_mgLevel3 = new double[ncv_g_mgLevel3];
		
		rho_old_mgLevel3 = new double[ncv_g_mgLevel3];
		rhou_old_mgLevel3 = new double[ncv_g_mgLevel3][3];
		rhoE_old_mgLevel3 = new double[ncv_g_mgLevel3];
		
		vel_mgLevel3 = new double[ncv_g_mgLevel3][3];
		press_mgLevel3 = new double[ncv_g_mgLevel3];
		temp_mgLevel3 = new double[ncv_g_mgLevel3];
		enthalpy_mgLevel3 = new double[ncv_g_mgLevel3];
		RoM_mgLevel3 = new double[ncv_g_mgLevel3];
		muT_mgLevel3 = new double[ncv_g_mgLevel3];
		gamma_mgLevel3 = new double[ncv_g_mgLevel3];
		sos_mgLevel3 = new double[ncv_g_mgLevel3];
		kine_mgLevel3 = new double[ncv_g_mgLevel3];
		local_dt_mgLevel3 = new double[ncv_g_mgLevel3];
		grad_enthalpy_mgLevel3 = new double[ncv_g_mgLevel3][3];
		grad_u_mgLevel3 = new double[ncv_g_mgLevel3][3][3];
		muLam_mgLevel3 = new double[ncv_g_mgLevel3];
		LambdaOverCp_mgLevel3 = new double[ncv_g_mgLevel3];			
		
		
		rho_TruncError_mgLevel3 = new double[ncv_g_mgLevel3];
		rhou_TruncError_mgLevel3 = new double[ncv_g_mgLevel3][3]; 
		rhoE_TruncError_mgLevel3 = new double[ncv_g_mgLevel3];
		
		rho_bfa_mgLevel3 = new double[nfa_mgLevel3];
		T_bfa_mgLevel3 = new double[nfa_mgLevel3];
		vel_bfa_mgLevel3 = new double[nfa_mgLevel3][3];
		p_bfa_mgLevel3 = new double[nfa_mgLevel3];
		h_bfa_mgLevel3 = new double[nfa_mgLevel3];
		gam_bfa_mgLevel3 = new double[nfa_mgLevel3];
		RoM_bfa_mgLevel3 = new double[nfa_mgLevel3];
		mul_fa_mgLevel3 = new double[nfa_mgLevel3];
		mut_fa_mgLevel3 = new double[nfa_mgLevel3];
		lamOcp_fa_mgLevel3 = new double[nfa_mgLevel3];
		
	}
}	

int JoeWithModels::findFace_mgLevel1(int cv_first, int cv_second) {
	
	if (cv_first == cv_second) return -1;
	int faocv_first_f = faocv_i_mgLevel1[cv_first];
	int faocv_first_l = faocv_i_mgLevel1[cv_first + 1] - 1;
	for (int faocv_first = faocv_first_f; faocv_first <= faocv_first_l; faocv_first++) {
		int faocv_second_f = faocv_i_mgLevel1[cv_second];
		int faocv_second_l = faocv_i_mgLevel1[cv_second + 1] - 1;
		for (int faocv_second = faocv_second_f; faocv_second <= faocv_second_l; faocv_second++)
			if (faocv_v_mgLevel1[faocv_first] == faocv_v_mgLevel1[faocv_second]) 
				return faocv_v_mgLevel1[faocv_first];
	}
	return -1;
}

int JoeWithModels::findFace_mgLevel2(int cv_first, int cv_second) {
	
	if (cv_first == cv_second) return -1;
	int faocv_first_f = faocv_i_mgLevel2[cv_first];
	int faocv_first_l = faocv_i_mgLevel2[cv_first + 1] - 1;
	for (int faocv_first = faocv_first_f; faocv_first <= faocv_first_l; faocv_first++) {
		int faocv_second_f = faocv_i_mgLevel2[cv_second];
		int faocv_second_l = faocv_i_mgLevel2[cv_second + 1] - 1;
		for (int faocv_second = faocv_second_f; faocv_second <= faocv_second_l; faocv_second++)
			if (faocv_v_mgLevel2[faocv_first] == faocv_v_mgLevel2[faocv_second]) 
				return faocv_v_mgLevel2[faocv_first];
	}
	return -1;
}

int JoeWithModels::findFace_mgLevel3(int cv_first, int cv_second) {
	
	if (cv_first == cv_second) return -1;
	int faocv_first_f = faocv_i_mgLevel3[cv_first];
	int faocv_first_l = faocv_i_mgLevel3[cv_first + 1] - 1;
	for (int faocv_first = faocv_first_f; faocv_first <= faocv_first_l; faocv_first++) {
		int faocv_second_f = faocv_i_mgLevel3[cv_second];
		int faocv_second_l = faocv_i_mgLevel3[cv_second + 1] - 1;
		for (int faocv_second = faocv_second_f; faocv_second <= faocv_second_l; faocv_second++)
			if (faocv_v_mgLevel3[faocv_first] == faocv_v_mgLevel3[faocv_second]) 
				return faocv_v_mgLevel3[faocv_first];
	}
	return -1;
}

double JoeWithModels::checkFaceOrientation(int ifa, int cv_first, int cv_second) {
	if ((cvofa[ifa][0] == cv_first) && (cvofa[ifa][1] == cv_second)) return 1.0;
	else return -1.0;
}

double JoeWithModels::checkFaceOrientation_mgLevel1(int ifa, int cv_first, int cv_second) {
	if ((cvofa_mgLevel1[ifa][0] == cv_first) && (cvofa_mgLevel1[ifa][1] == cv_second)) return 1.0;
	else return -1.0;
}

double JoeWithModels::checkFaceOrientation_mgLevel2(int ifa, int cv_first, int cv_second) {
	if ((cvofa_mgLevel2[ifa][0] == cv_first) && (cvofa_mgLevel2[ifa][1] == cv_second)) return 1.0;
	else return -1.0;
}

void JoeWithModels::calcTimeInt(double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double (*A)[5][5], double (*dq)[5], 
								double (*rhs)[5], double underRelax, int maxIterLS, double zeroAbsLS, double zeroRelLS)
{
	for (int icv = 0; icv < ncv_g; icv++)
		for (int i = 0; i < 5; i++)
			dq[icv][i] = 0.0;
	
	// solve linear system for the NSE
	for (int icv=0; icv < ncv; ++icv) {
		rhs[icv][0] = underRelax*(RHSrho[icv]);
		rhs[icv][1] = underRelax*(RHSrhou[icv][0]);
		rhs[icv][2] = underRelax*(RHSrhou[icv][1]);
		rhs[icv][3] = underRelax*(RHSrhou[icv][2]);
		rhs[icv][4] = underRelax*(RHSrhoE[icv]);
		
		residField[icv] = RHSrhoE[icv];
		
		double tmp = cv_volume[icv]/(local_dt[icv]);
		for (int i = 0; i < 5; i++)
			A[nbocv_i[icv]][i][i] += tmp;                               // diagonal part ( vol/dt + A )
	}
	
	// solve linear system
	solveCoupledLinSysNS(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS);
	
	// update solution: Q_new = Q_old + Delta_Q
	for (int icv=0; icv<ncv; icv++)  {
		rho[icv]     += dq[icv][0];
		rhou[icv][0] += dq[icv][1];
		rhou[icv][1] += dq[icv][2];
		rhou[icv][2] += dq[icv][3];
		rhoE[icv]    += dq[icv][4];
	}
	
	updateCvData(rho,  REPLACE_DATA);
  	updateCvData(rhou, REPLACE_ROTATE_DATA);
	updateCvData(rhoE, REPLACE_DATA);
	
	UpdateCvDataStateVec(dq);                                       // update dq since neighbors needed to compute RHS of scalars
	
}

void JoeWithModels::setHistoryFile()
{
  FILE *fp;
  char fname[200] = "history.dat";
  int nScal = scalarTranspEqVector.size();

  if (mpi_rank == 0)
  {
    if ( (fp = fopen(fname,"wt"))== NULL )
    {
      cerr << "Error: cannot open file" << fname << endl;
      throw(-1);
    }
    fprintf(fp, "VARIABLES = \"Iter\" \"Res_rho\" \"Res_rhou-X\" \"Res_rhou-Y\" \"Res_rhou-Z\" \"Res_rhoE\"");
    for (int iScal = 0; iScal < nScal; iScal++)
      fprintf(fp," \"Res_%s\"",scalarTranspEqVector[iScal].getName());
    fprintf(fp,"\nZone T = \"residuals\" F = point\n");
    
    fclose(fp);
  }
}

void JoeWithModels::writeHistoryFile(double *rhsResid)
{
  FILE *fp;
  char fname[200] = "history.dat";
  int nScal = scalarTranspEqVector.size();

  if (mpi_rank == 0)
  {
    if ( (fp = fopen(fname,"a"))== NULL )
    {
      cerr << "Error: cannot open file" << fname << endl;
      throw(-1);
    }
    fprintf(fp,"%8d",step);
    fprintf(fp,"%14.4e%14.4e%14.4e%14.4e%14.4e",
        log10(rhsResid[0]), log10(rhsResid[1]), log10(rhsResid[2]),
        log10(rhsResid[3]), log10(rhsResid[4]));

    for (int iScal = 0; iScal < nScal; iScal++)
      fprintf(fp,"%14.4e",log10(rhsResid[5 + iScal]));

    fprintf(fp,"\n");
    fclose(fp);
  }
}

void JoeWithModels::calcJSTCoeff_Const(int iMesh) {
	
	if (iMesh == 0) {
		for (int icv = 0; icv < ncv_g; icv++) {
			BoundaryCV[icv] = false;
			NeighborCV[icv] = 0;
		}
		
		for (int ifa = nfa_b; ifa < nfa; ifa++) {
			BoundaryCV[cvofa[ifa][0]] = true;
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			NeighborCV[icv0]++;
			NeighborCV[icv1]++;
		}
	}
	
	if (iMesh == 1) {
		for (int icv = 0; icv < ncv_g_mgLevel1; icv++)
			NeighborCV_mgLevel1[icv] = 0;
		
		for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
			int icv0 = cvofa_mgLevel1[ifa][0];
			int icv1 = cvofa_mgLevel1[ifa][1];
			NeighborCV_mgLevel1[icv0]++;
			NeighborCV_mgLevel1[icv1]++;
		}
	}
	
	if (iMesh == 2) {
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++)
			NeighborCV_mgLevel2[icv] = 0;
		
		for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
			int icv0 = cvofa_mgLevel2[ifa][0];
			int icv1 = cvofa_mgLevel2[ifa][1];
			NeighborCV_mgLevel2[icv0]++;
			NeighborCV_mgLevel2[icv1]++;
		}
	}
	
	if (iMesh == 3) {
		for (int icv = 0; icv < ncv_g_mgLevel3; icv++)
			NeighborCV_mgLevel3[icv] = 0;
		
		for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
			int icv0 = cvofa_mgLevel3[ifa][0];
			int icv1 = cvofa_mgLevel3[ifa][1];
			NeighborCV_mgLevel3[icv0]++;
			NeighborCV_mgLevel3[icv1]++;
		}
	}
	
}

void JoeWithModels::calcJSTCoeff_Var(int iMesh) {
	
	if (iMesh == 0) {
		
		double Diff[5];
		
		for (int icv = 0; icv < ncv_g; icv++) {
			p1_Und_Lapl[icv] = 0.0;
			p2_Und_Lapl[icv] = 0.0;
			for (int iVar = 0; iVar < 5; iVar++)
				Conserv_Und_Lapl[icv][iVar] = 0.0;
			Lambda[icv] = 0.0;
		}
		
		for (int ifa = nfa_b; ifa < nfa; ifa++) {
			
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal[ifa]);
			double Pressure_0 = press[icv0]; double Pressure_1 = press[icv1];
			double rho0 = rho[icv0]; double rho1 = rho[icv1];
			double gamma0 = gamma[icv0]; double gamma1 = gamma[icv1];
			double SoundSpeed0 = sqrt(gamma0*Pressure_0/rho0); 
			double SoundSpeed1 = sqrt(gamma1*Pressure_1/rho1);
			
			double u0[3] = {rhou[icv0][0]/rho[icv0], rhou[icv0][1]/rho[icv0], rhou[icv0][2]/rho[icv0]};
			double u1[3] = {rhou[icv1][0]/rho[icv1], rhou[icv1][1]/rho[icv1], rhou[icv1][2]/rho[icv1]};
			double un0  = vecDotVec3d(u0, nVec);
			double un1  = vecDotVec3d(u1, nVec);
			double Lambda0 = 0.5*(fabs(un0) + SoundSpeed0)*area;
			double Lambda1 = 0.5*(fabs(un1) + SoundSpeed1)*area;
			bool boundary_0 = BoundaryCV[icv0]; 
			bool boundary_1 = BoundaryCV[icv1]; 
			
			Lambda[icv0] += Lambda0; Lambda[icv0] += Lambda1;
			Lambda[icv1] += Lambda0; Lambda[icv1] += Lambda1;
			
			p1_Und_Lapl[icv0] += (Pressure_1 - Pressure_0); p1_Und_Lapl[icv1] += (Pressure_0 - Pressure_1);
			p2_Und_Lapl[icv0] += (Pressure_0 + Pressure_1); p2_Und_Lapl[icv1] += (Pressure_0 + Pressure_1);
			
			Diff[0] = rho[icv0] - rho[icv1];		
			Diff[1] = rhou[icv0][0] - rhou[icv1][0];	
			Diff[2] = rhou[icv0][1] - rhou[icv1][1];	
			Diff[3] = rhou[icv0][2] - rhou[icv1][2];		
			Diff[4] = (rhoE[icv0]+Pressure_0) - (rhoE[icv1]+Pressure_1);
			
			// Both points inside Omega
			if ((!boundary_0 && !boundary_1) || (boundary_0 && boundary_1)) {
				for (int iVar = 0; iVar < 5; iVar++) Conserv_Und_Lapl[icv0][iVar] -= Diff[iVar];
				for (int iVar = 0; iVar < 5; iVar++) Conserv_Und_Lapl[icv1][iVar] += Diff[iVar];
			}
			// icv0 inside Omega, icv1 on the boundary 
			if (!boundary_0 && boundary_1)
				for (int iVar = 0; iVar < 5; iVar++) Conserv_Und_Lapl[icv0][iVar] -= Diff[iVar];
			// icv1 inside Omega, icv0 on the boundary
			if (boundary_0 && !boundary_1)
				for (int iVar = 0; iVar < 5; iVar++) Conserv_Und_Lapl[icv1][iVar] += Diff[iVar];
		}
		
		for (int icv = 0; icv < ncv_g; icv++)
			PressSensor[icv] = fabs(p1_Und_Lapl[icv])/p2_Und_Lapl[icv];
	}
	
	if (iMesh == 1) {
		for (int icv = 0; icv < ncv_g_mgLevel1; icv++) Lambda_mgLevel1[icv] = 0.0;
		
		for (int ifa = nfa_b_mgLevel1; ifa < nfa_mgLevel1; ifa++) {
			int icv0 = cvofa_mgLevel1[ifa][0];
			int icv1 = cvofa_mgLevel1[ifa][1];
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel1[ifa]);
			double Pressure_0 = press_mgLevel1[icv0]; double Pressure_1 = press_mgLevel1[icv1];
			double rho0 = rho_mgLevel1[icv0]; double rho1 = rho_mgLevel1[icv1];
			double gamma0 = gamma_mgLevel1[icv0]; double gamma1 = gamma_mgLevel1[icv1];
			double SoundSpeed0 = sqrt(gamma0*Pressure_0/rho0); 
			double SoundSpeed1 = sqrt(gamma1*Pressure_1/rho1);
			double u0[3] = {rhou_mgLevel1[icv0][0]/rho_mgLevel1[icv0], rhou_mgLevel1[icv0][1]/rho_mgLevel1[icv0], rhou_mgLevel1[icv0][2]/rho_mgLevel1[icv0]};
			double u1[3] = {rhou_mgLevel1[icv1][0]/rho_mgLevel1[icv1], rhou_mgLevel1[icv1][1]/rho_mgLevel1[icv1], rhou_mgLevel1[icv1][2]/rho_mgLevel1[icv1]};
			double un0  = vecDotVec3d(u0, nVec); double un1  = vecDotVec3d(u1, nVec);
			double Lambda0 = 0.5*(fabs(un0) + SoundSpeed0)*area; double Lambda1 = 0.5*(fabs(un1) + SoundSpeed1)*area;

			Lambda_mgLevel1[icv0] += Lambda0; Lambda_mgLevel1[icv0] += Lambda1;
			Lambda_mgLevel1[icv1] += Lambda0; Lambda_mgLevel1[icv1] += Lambda1;
		}
	}
	
	if (iMesh == 2) {
		for (int icv = 0; icv < ncv_g_mgLevel2; icv++) Lambda_mgLevel2[icv] = 0.0;
		
		for (int ifa = nfa_b_mgLevel2; ifa < nfa_mgLevel2; ifa++) {
			int icv0 = cvofa_mgLevel2[ifa][0];
			int icv1 = cvofa_mgLevel2[ifa][1];
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel2[ifa]);
			double Pressure_0 = press_mgLevel2[icv0]; double Pressure_1 = press_mgLevel2[icv1];
			double rho0 = rho_mgLevel2[icv0]; double rho1 = rho_mgLevel2[icv1];
			double gamma0 = gamma_mgLevel2[icv0]; double gamma1 = gamma_mgLevel2[icv1];
			double SoundSpeed0 = sqrt(gamma0*Pressure_0/rho0); 
			double SoundSpeed1 = sqrt(gamma1*Pressure_1/rho1);
			double u0[3] = {rhou_mgLevel2[icv0][0]/rho_mgLevel2[icv0], rhou_mgLevel2[icv0][1]/rho_mgLevel2[icv0], rhou_mgLevel2[icv0][2]/rho_mgLevel2[icv0]};
			double u1[3] = {rhou_mgLevel2[icv1][0]/rho_mgLevel2[icv1], rhou_mgLevel2[icv1][1]/rho_mgLevel2[icv1], rhou_mgLevel2[icv1][2]/rho_mgLevel2[icv1]};
			double un0  = vecDotVec3d(u0, nVec); double un1  = vecDotVec3d(u1, nVec);
			double Lambda0 = 0.5*(fabs(un0) + SoundSpeed0)*area; double Lambda1 = 0.5*(fabs(un1) + SoundSpeed1)*area;
			
			Lambda_mgLevel2[icv0] += Lambda0; Lambda_mgLevel2[icv0] += Lambda1;
			Lambda_mgLevel2[icv1] += Lambda0; Lambda_mgLevel2[icv1] += Lambda1;
		}
	}
	
	if (iMesh == 3) {
		for (int icv = 0; icv < ncv_g_mgLevel3; icv++) Lambda_mgLevel3[icv] = 0.0;
		
		for (int ifa = nfa_b_mgLevel3; ifa < nfa_mgLevel3; ifa++) {
			int icv0 = cvofa_mgLevel3[ifa][0];
			int icv1 = cvofa_mgLevel3[ifa][1];
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal_mgLevel3[ifa]);
			double Pressure_0 = press_mgLevel3[icv0]; double Pressure_1 = press_mgLevel3[icv1];
			double rho0 = rho_mgLevel3[icv0]; double rho1 = rho_mgLevel3[icv1];
			double gamma0 = gamma_mgLevel3[icv0]; double gamma1 = gamma_mgLevel3[icv1];
			double SoundSpeed0 = sqrt(gamma0*Pressure_0/rho0); 
			double SoundSpeed1 = sqrt(gamma1*Pressure_1/rho1);
			double u0[3] = {rhou_mgLevel3[icv0][0]/rho_mgLevel3[icv0], rhou_mgLevel3[icv0][1]/rho_mgLevel3[icv0], rhou_mgLevel3[icv0][2]/rho_mgLevel3[icv0]};
			double u1[3] = {rhou_mgLevel3[icv1][0]/rho_mgLevel3[icv1], rhou_mgLevel3[icv1][1]/rho_mgLevel3[icv1], rhou_mgLevel3[icv1][2]/rho_mgLevel3[icv1]};
			double un0  = vecDotVec3d(u0, nVec); double un1  = vecDotVec3d(u1, nVec);
			double Lambda0 = 0.5*(fabs(un0) + SoundSpeed0)*area; double Lambda1 = 0.5*(fabs(un1) + SoundSpeed1)*area;
			
			Lambda_mgLevel3[icv0] += Lambda0; Lambda_mgLevel3[icv0] += Lambda1;
			Lambda_mgLevel3[icv1] += Lambda0; Lambda_mgLevel3[icv1] += Lambda1;
		}
	}
	
}

int JoeWithModels::FindEdgeCV(int ip_0, int ip_1, int **Edge, int *nbono_i, int *nbono_v) {
	
	int ino = 0, iNode;
	for (iNode = nbono_i[ip_0]; iNode <= nbono_i[ip_0+1]-1; iNode++) {
		ino = nbono_v[iNode];
		if (ino == ip_1) break;
	}
	
	if (ino == ip_1) return Edge[ip_0][iNode - nbono_i[ip_0]];	
	else {
		cout <<"Ups! I don't find the edge that connect "<< ip_0 <<" and "<< ip_1 <<"."<< endl;
		exit(1);
		return -1;
	}
}

int JoeWithModels::FindFaceCV(int icv, int ip_0, int ip_1, int ip_2, int ip_3) {
	
	int counter, ifa;
	
	for (int foc = faocv_i[icv]; foc <= faocv_i[icv + 1]-1; foc++) {
		ifa = faocv_v[foc];
		counter = 0;
    for (int nof = noofa_i[ifa]; nof <= noofa_i[ifa+1]-1; nof++)  {
      int ino = noofa_v[nof];
			if ((ino == ip_0) || (ino == ip_1) || 
					(ino == ip_2) || (ino == ip_3) ) counter++;
		}
		if (counter == 4) break;
	}
	
	if (counter == 4) return ifa;	
	else {
		cout <<"Ups! I don't find the face that connect "<< ip_0 <<", "<< ip_1 <<", "<< ip_2 <<", "<< ip_3 <<"."<< endl;
		exit(1);
		return -1;
	}
}

int JoeWithModels::CheckTetCode(bool *AdaptCode) {
	
	int Code = -1;
	
	// Default
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true) || (AdaptCode[4] == true) ||
			(AdaptCode[5] == true) ) { Code = 0; }
	
	// Combination 1:8
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) ) {Code = 1; return Code;}
	
	// Combination 1:4
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) ) {Code = 2; return Code;}

	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) ) {Code = 3; return Code;}

	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) ) {Code = 4; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) ) {Code = 5; return Code;}

	// Combination 1:2
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) ) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) ) {Code = 7; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) ) {Code = 8; return Code; }
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) ) {Code = 9; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) ) {Code = 10; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) ) {Code = 11; return Code; }
	
	return Code;
	
}

int JoeWithModels::CheckHexaCode(bool *AdaptCode) {
	
	int Code = -1;
	
	// Default

	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true) || (AdaptCode[4] == true) ||
			(AdaptCode[5] == true) || (AdaptCode[6] == true) || (AdaptCode[7] == true) || (AdaptCode[8] == true) || (AdaptCode[9] == true) ||
			(AdaptCode[10] == true) || (AdaptCode[11] == true)) { Code = 0; }
	
	// Combination 1:8

	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == true) && (AdaptCode[7] == true) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 1; return Code;}
	
	// Combination 1:4

	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == true) && (AdaptCode[7] == false) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == false) && (AdaptCode[7] == true) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == true) && (AdaptCode[7] == true) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 4; return Code;}
	
	// Combination 1:2

	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == false) && (AdaptCode[7] == false) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 5; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == false) && (AdaptCode[7] == true) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == true) && (AdaptCode[7] == false) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 7; return Code; }
			
		return Code;
	
}

int JoeWithModels::CheckPyramCode(bool *AdaptCode) {
	
	int Code = -1;
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true)) {Code = 0;}
	
	// Combination 1P -> 1P
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 1; return Code;}
	
	// Combination 1P -> 3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 4; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 5; return Code;}
	
	// Combination 1P -> 4T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 7; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 8; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 9; return Code;}
	
	// Combination 1P -> 2P
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 10; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 11; return Code;}
	
	// Combination 1P -> 1P+3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 12; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 13; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 14; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 15; return Code;}
	
	
	// Combination 1P -> 4P
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 16; return Code;}
	
	return Code;
	
}

void JoeWithModels::TetDivision(int code , int *nodes, int **Division, int *nPart) {
	
	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		Division[4][0] = 5; Division[5][0] = 5; Division[6][0] = 5; Division[7][0] = 5;
		*nPart = 8;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[6];	Division[0][3] = nodes[4]; Division[0][4] = nodes[0]; 
		Division[1][1] = nodes[7];	Division[1][2] = nodes[4];	Division[1][3] = nodes[5]; Division[1][4] = nodes[1]; 
		Division[2][1] = nodes[8];	Division[2][2] = nodes[7];	Division[2][3] = nodes[9]; Division[2][4] = nodes[3]; 
		Division[3][1] = nodes[9];	Division[3][2] = nodes[5];	Division[3][3] = nodes[6]; Division[3][4] = nodes[2]; 
		Division[4][1] = nodes[8];	Division[4][2] = nodes[7];	Division[4][3] = nodes[5]; Division[4][4] = nodes[9]; 
		Division[5][1] = nodes[8];	Division[5][2] = nodes[5];	Division[5][3] = nodes[6]; Division[5][4] = nodes[9]; 
		Division[6][1] = nodes[5];	Division[6][2] = nodes[7];	Division[6][3] = nodes[8]; Division[6][4] = nodes[4]; 
		Division[7][1] = nodes[6];	Division[7][2] = nodes[5];	Division[7][3] = nodes[8]; Division[7][4] = nodes[4]; 
		
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[8];	Division[0][4] = nodes[2]; 
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[7]; Division[1][4] = nodes[2]; 
		Division[2][1] = nodes[8];	Division[2][2] = nodes[4];	Division[2][3] = nodes[7]; Division[2][4] = nodes[2]; 
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[2]; 
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[1];	Division[0][2] = nodes[7];	Division[0][3] = nodes[5];	Division[0][4] = nodes[0]; 
		Division[1][1] = nodes[7];	Division[1][2] = nodes[3];	Division[1][3] = nodes[9]; Division[1][4] = nodes[0]; 
		Division[2][1] = nodes[7];	Division[2][2] = nodes[9];	Division[2][3] = nodes[5]; Division[2][4] = nodes[0]; 
		Division[3][1] = nodes[9];	Division[3][2] = nodes[2];	Division[3][3] = nodes[5]; Division[3][4] = nodes[0]; 
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[6];	Division[0][3] = nodes[0];	Division[0][4] = nodes[1]; 
		Division[1][1] = nodes[8];	Division[1][2] = nodes[9];	Division[1][3] = nodes[6]; Division[1][4] = nodes[1]; 
		Division[2][1] = nodes[9];	Division[2][2] = nodes[2];	Division[2][3] = nodes[6]; Division[2][4] = nodes[1]; 
		Division[3][1] = nodes[3];	Division[3][2] = nodes[9];	Division[3][3] = nodes[8]; Division[3][4] = nodes[1]; 
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[6];	Division[0][4] = nodes[3]; 
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[5]; Division[1][4] = nodes[3]; 
		Division[2][1] = nodes[6];	Division[2][2] = nodes[4];	Division[2][3] = nodes[5]; Division[2][4] = nodes[3]; 
		Division[3][1] = nodes[6];	Division[3][2] = nodes[5];	Division[3][3] = nodes[2]; Division[3][4] = nodes[3]; 
		
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[2];	Division[0][4] = nodes[3]; 
		Division[1][1] = nodes[2];	Division[1][2] = nodes[4];	Division[1][3] = nodes[1]; Division[1][4] = nodes[3]; 
		
	}	

	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[2];	Division[0][2] = nodes[0];	Division[0][3] = nodes[5];	Division[0][4] = nodes[3]; 
		Division[1][1] = nodes[5];	Division[1][2] = nodes[0];	Division[1][3] = nodes[1]; Division[1][4] = nodes[3]; 
		
	}	

	if (code == 8) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[6];	Division[0][2] = nodes[0];	Division[0][3] = nodes[1];	Division[0][4] = nodes[3]; 
		Division[1][1] = nodes[2];	Division[1][2] = nodes[6];	Division[1][3] = nodes[1]; Division[1][4] = nodes[3]; 
		
	}	

	if (code == 9) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[7];	Division[0][2] = nodes[1];	Division[0][3] = nodes[0];	Division[0][4] = nodes[2]; 
		Division[1][1] = nodes[3];	Division[1][2] = nodes[7];	Division[1][3] = nodes[0]; Division[1][4] = nodes[2]; 
		
	}	

	if (code == 10) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[8]; 
		Division[1][1] = nodes[1];	Division[1][2] = nodes[3];	Division[1][3] = nodes[2]; Division[1][4] = nodes[8]; 
		
	}	

	if (code == 11) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[9]; 
		Division[1][1] = nodes[0];	Division[1][2] = nodes[1];	Division[1][3] = nodes[9]; Division[1][4] = nodes[3]; 
		
	}		
	
}

void JoeWithModels::HexaDivision(int code , int *nodes, int **Division, int *nPart) {
	
	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		Division[4][0] = 9; Division[5][0] = 9; Division[6][0] = 9; Division[7][0] = 9;
		*nPart = 8;
		
		// nodes that compose each element
		Division[0][1] = nodes[20];	Division[0][2] = nodes[8];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[26]; Division[0][6] = nodes[25]; Division[0][7] = nodes[17]; Division[0][8] = nodes[22];
		
		Division[1][1] = nodes[20];	Division[1][2] = nodes[9];	Division[1][3] = nodes[2]; Division[1][4] = nodes[10]; 
		Division[1][5] = nodes[26]; Division[1][6] = nodes[22]; Division[1][7] = nodes[18]; Division[1][8] = nodes[23];
		
		Division[2][1] = nodes[20];	Division[2][2] = nodes[10];	Division[2][3] = nodes[3]; Division[2][4] = nodes[11]; 
		Division[2][5] = nodes[26]; Division[2][6] = nodes[23]; Division[2][7] = nodes[19]; Division[2][8] = nodes[24];
		
		Division[3][1] = nodes[20];	Division[3][2] = nodes[11];	Division[3][3] = nodes[0]; Division[3][4] = nodes[8]; 
		Division[3][5] = nodes[26]; Division[3][6] = nodes[24]; Division[3][7] = nodes[16]; Division[3][8] = nodes[25];
		
		Division[4][1] = nodes[26];	Division[4][2] = nodes[25];	Division[4][3] = nodes[17];	Division[4][4] = nodes[22]; 
		Division[4][5] = nodes[21]; Division[4][6] = nodes[12]; Division[4][7] = nodes[5]; Division[4][8] = nodes[13];
		
		Division[5][1] = nodes[26];	Division[5][2] = nodes[22];	Division[5][3] = nodes[18]; Division[5][4] = nodes[23]; 
		Division[5][5] = nodes[21]; Division[5][6] = nodes[13]; Division[5][7] = nodes[6]; Division[5][8] = nodes[14];
		
		Division[6][1] = nodes[26];	Division[6][2] = nodes[23];	Division[6][3] = nodes[19]; Division[6][4] = nodes[24]; 
		Division[6][5] = nodes[21]; Division[6][6] = nodes[14]; Division[6][7] = nodes[7]; Division[6][8] = nodes[15];
		
		Division[7][1] = nodes[26];	Division[7][2] = nodes[24];	Division[7][3] = nodes[16]; Division[7][4] = nodes[25]; 
		Division[7][5] = nodes[21]; Division[7][6] = nodes[15]; Division[7][7] = nodes[4]; Division[7][8] = nodes[12];
		
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[10]; 
		Division[0][5] = nodes[25]; Division[0][17] = nodes[25]; Division[0][7] = nodes[18]; Division[0][8] = nodes[23];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[8];	Division[1][3] = nodes[10]; Division[1][4] = nodes[3]; 
		Division[1][5] = nodes[16]; Division[1][6] = nodes[25]; Division[1][7] = nodes[23]; Division[1][8] = nodes[19];
		
		Division[2][1] = nodes[25];	Division[2][2] = nodes[17];	Division[2][3] = nodes[18]; Division[2][4] = nodes[23]; 
		Division[2][5] = nodes[12]; Division[2][6] = nodes[5]; Division[2][7] = nodes[6]; Division[2][8] = nodes[14];
		
		Division[3][1] = nodes[16];	Division[3][2] = nodes[25];	Division[3][3] = nodes[23]; Division[3][4] = nodes[19]; 
		Division[3][5] = nodes[4]; Division[3][6] = nodes[12]; Division[3][7] = nodes[14]; Division[3][8] = nodes[7];
				
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		*nPart = 4;

		// nodes that compose each element
		Division[0][1] = nodes[11];	Division[0][2] = nodes[0];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[24]; Division[0][6] = nodes[16]; Division[0][7] = nodes[17]; Division[0][8] = nodes[22];
		
		Division[1][1] = nodes[3];	Division[1][2] = nodes[11];	Division[1][9] = nodes[2]; Division[1][4] = nodes[2]; 
		Division[1][5] = nodes[19]; Division[1][6] = nodes[24]; Division[1][7] = nodes[22]; Division[1][8] = nodes[18];
		
		Division[2][1] = nodes[24];	Division[2][2] = nodes[16];	Division[2][3] = nodes[17]; Division[2][4] = nodes[22]; 
		Division[2][5] = nodes[15]; Division[2][6] = nodes[4]; Division[2][7] = nodes[5]; Division[2][8] = nodes[13];
		
		Division[3][1] = nodes[19];	Division[3][2] = nodes[24];	Division[3][3] = nodes[22]; Division[3][4] = nodes[18]; 
		Division[3][5] = nodes[7]; Division[3][6] = nodes[15]; Division[3][7] = nodes[13]; Division[3][8] = nodes[6];
		
	}	

	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9;
		*nPart = 4;

		// nodes that compose each element
		Division[0][1] = nodes[20];	Division[0][2] = nodes[8];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[21]; Division[0][6] = nodes[12]; Division[0][7] = nodes[5]; Division[0][8] = nodes[13];
		
		Division[1][1] = nodes[20];	Division[1][2] = nodes[9];	Division[1][3] = nodes[2]; Division[1][4] = nodes[10]; 
		Division[1][5] = nodes[21]; Division[1][6] = nodes[13]; Division[1][7] = nodes[6]; Division[1][8] = nodes[14];
		
		Division[2][1] = nodes[20];	Division[2][2] = nodes[10];	Division[2][3] = nodes[0]; Division[2][4] = nodes[8]; 
		Division[2][5] = nodes[21]; Division[2][6] = nodes[15]; Division[2][7] = nodes[4]; Division[2][8] = nodes[12];
		
		Division[3][1] = nodes[20];	Division[3][2] = nodes[10];	Division[3][3] = nodes[3]; Division[3][4] = nodes[11]; 
		Division[3][5] = nodes[21]; Division[3][6] = nodes[14]; Division[3][7] = nodes[7]; Division[3][8] = nodes[15];
		
	}	

	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;

		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[3]; 
		Division[0][5] = nodes[16]; Division[0][6] = nodes[17]; Division[0][7] = nodes[18]; Division[0][8] = nodes[19];
		
		Division[1][1] = nodes[16];	Division[1][2] = nodes[17];	Division[1][3] = nodes[18]; Division[1][4] = nodes[19]; 
		Division[1][5] = nodes[4]; Division[1][6] = nodes[5]; Division[1][7] = nodes[6]; Division[1][8] = nodes[7];
				
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;

		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[9];	Division[0][4] = nodes[11]; 
		Division[0][5] = nodes[4]; Division[0][6] = nodes[5]; Division[0][7] = nodes[13]; Division[0][8] = nodes[15];
		
		Division[1][1] = nodes[9];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3]; Division[1][4] = nodes[11]; 
		Division[1][5] = nodes[13]; Division[1][6] = nodes[6]; Division[1][7] = nodes[7]; Division[1][8] = nodes[15];
		
	}	
	
	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;

		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[10]; 
		Division[0][5] = nodes[12]; Division[0][5] = nodes[5]; Division[0][7] = nodes[6]; Division[0][8] = nodes[14];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[8];	Division[1][3] = nodes[10]; Division[1][4] = nodes[3]; 
		Division[1][5] = nodes[4]; Division[1][6] = nodes[12]; Division[1][7] = nodes[14]; Division[1][8] = nodes[7];
		
	}	
	
}

void JoeWithModels::PyramDivision(int code , int *nodes, int **Division, int *nPart) {
	
	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 6;
		*nPart = 1;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[3]; Division[0][5] = nodes[4]; 
	
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[5];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3]; Division[1][4] = nodes[4];
				
		Division[2][1] = nodes[5];	Division[2][2] = nodes[3];	Division[2][3] = nodes[0]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[6];	Division[1][3] = nodes[3]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[3];	Division[2][2] = nodes[6];	Division[2][3] = nodes[2]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[7]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[0];	Division[2][2] = nodes[7];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[8];	Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[2];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[3]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[6]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[6];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[6];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[6];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[6];	Division[2][2] = nodes[2];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[0];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 8) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[8]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[8];	Division[1][2] = nodes[1];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[1];	Division[2][2] = nodes[2];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 9) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[8]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[5];	Division[2][3] = nodes[2]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 10) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[7];	Division[0][4] = nodes[3]; Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2];  Division[1][4] = nodes[7]; Division[1][5] = nodes[4];
				
	}	
	
	if (code == 11) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6];	Division[0][4] = nodes[8]; Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[8];	Division[1][2] = nodes[6];	Division[1][3] = nodes[2];  Division[1][4] = nodes[3]; Division[1][5] = nodes[4];
		
	}	
	
	if (code == 12) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[7];	Division[0][4] = nodes[3];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[6]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[6];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[6];	Division[3][2] = nodes[2];	Division[3][3] = nodes[7]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 13) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[8];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[6];	Division[1][2] = nodes[2];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[6];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 14) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[5];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[7];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[5];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 15) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[6];	Division[0][3] = nodes[2]; Division[0][4] = nodes[3];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[1];	Division[2][3] = nodes[6]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[5];	Division[3][2] = nodes[6];	Division[3][3] = nodes[8]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 16) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6; Division[2][0] = 6; Division[3][0] = 6;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[9];	Division[0][2] = nodes[8];	Division[0][3] = nodes[0]; Division[0][4] = nodes[5];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[9];	Division[1][2] = nodes[5];	Division[1][3] = nodes[1]; Division[1][4] = nodes[6]; Division[1][5] = nodes[4];
		
		Division[2][1] = nodes[9];	Division[2][2] = nodes[6];	Division[2][3] = nodes[2]; Division[2][4] = nodes[7]; Division[2][5] = nodes[4]; 
		
		Division[3][1] = nodes[9];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[8]; Division[3][5] = nodes[4];
		
	}	

}




void JoeWithModels::AdaptInitialGrid()
{
		
	// Plot non-adapted solution using tecplot
	ofstream tecplot_file;
	tecplot_file.open("original.plt", ios::out);
	
	tecplot_file << " VARIABLES = \"x\",\"y\",\"z\" " << endl;	
	tecplot_file << " ZONE T = \"Time = 0.0\", N= "<< nno <<" , E = "<< ncv <<" , F = FEPOINT, ET = BRICK"<< endl;
	
	
	for (int ino = 0; ino < nno; ino++)
		tecplot_file << x_no[ino][0]<< " " << x_no[ino][1] << " " << x_no[ino][2] << endl;	
	
	for (int icv = 0; icv < ncv; icv++) {
		int noc_f = noocv_i[icv];
		int noc_l = noocv_i[icv+1]-1;
		
		for (int noc = noc_f; noc <= noc_l; noc++) {
			int ino = noocv_v[noc];
			tecplot_file << ino+1;
			if (noc != noc_l) tecplot_file << "\t";
			else tecplot_file<< "\n";
		}
	}
	
	
	// build a nodal connectivity 
	int * cvono_i_tmp = new int[nno+1];
	for (int ino = 0; ino < nno; ino++)
    cvono_i_tmp[ino+1] = 0;
	int * cvono_v_tmp = NULL;
	for (int iter = 0; iter < 2; iter++)
	{
		// put -1 in all nodes...
		for (int ino = 0; ino < nno; ino++)
      no_flag[ino] = -1;
		for (int icv = 0; icv < getNcv(); icv++)
		{
			const int foc_f = faocv_i[icv];
			const int foc_l = faocv_i[icv+1]-1;
			for (int foc = foc_f; foc <= foc_l; foc++)
			{
				const int ifa = faocv_v[foc];
				const int nof_f = noofa_i[ifa];
				const int nof_l = noofa_i[ifa+1]-1;
				for (int nof = nof_f; nof <= nof_l; nof++)
				{
					const int ino = noofa_v[nof];
					if (no_flag[ino] != icv)
					{
						no_flag[ino] = icv;
						if (iter == 0)
						{
							cvono_i_tmp[ino+1] += 1;
						}
						else
						{
							cvono_v_tmp[cvono_i_tmp[ino]] = icv;
							cvono_i_tmp[ino] += 1;
						}
					}
				}
			}
		}
		if (iter == 0)
		{
			cvono_i_tmp[0] = 0;
			for (int ino = 0; ino < nno; ino++)
        cvono_i_tmp[ino+1] += cvono_i_tmp[ino];
			cvono_v_tmp = new int[cvono_i_tmp[nno]];
		}
		else
		{
			for (int ino = nno; ino > 0; ino--)
        cvono_i_tmp[ino] = cvono_i_tmp[ino-1];
			cvono_i_tmp[0] = 0;
		}
	}
	
	int * nbono_i = new int[nno+1];
	for (int ino = 0; ino < nno; ino++)
    nbono_i[ino+1] = 0;
	int * nbono_v = NULL;
	for (int iter = 0; iter < 2; iter++)
	{
		// put -1 in all nodes...
		for (int ino = 0; ino < nno; ino++)
      no_flag[ino] = -1;
		for (int ino = 0; ino < nno; ino++)
		{
			// no need to include diagonal...
			assert( no_flag[ino] != ino );
			no_flag[ino] = ino;
			// get nbrs from cvono structure...
			int con_f = cvono_i_tmp[ino];
			int con_l = cvono_i_tmp[ino+1]-1;
			for (int con = con_f; con <= con_l; con++)
			{
				const int icv = cvono_v_tmp[con];
				const int nooc_f = noocv_i[icv];
				const int nooc_l = noocv_i[icv+1]-1;
				
				int nooc;
				for (nooc = nooc_f; nooc <= nooc_l; nooc++)
					if (ino == noocv_v[nooc]) break;
				
				int HexaCode[3] = {0,0,0};
				int nodeIndex = nooc - nooc_f;
				
				if(nodeIndex == 0) { HexaCode[0] = 1; HexaCode[1] = 3; HexaCode[2] = 4; }
				if(nodeIndex == 1) { HexaCode[0] = 2; HexaCode[1] = 0; HexaCode[2] = 5; }
				if(nodeIndex == 2) { HexaCode[0] = 3; HexaCode[1] = 1; HexaCode[2] = 6; }
				if(nodeIndex == 3) { HexaCode[0] = 2; HexaCode[1] = 0; HexaCode[2] = 7; }
				if(nodeIndex == 4) { HexaCode[0] = 0; HexaCode[1] = 5; HexaCode[2] = 7; }
				if(nodeIndex == 5) { HexaCode[0] = 4; HexaCode[1] = 6; HexaCode[2] = 1; }
				if(nodeIndex == 6) { HexaCode[0] = 7; HexaCode[1] = 5; HexaCode[2] = 2; }
				if(nodeIndex == 7) { HexaCode[0] = 4; HexaCode[1] = 6; HexaCode[2] = 3; }

				for (int iVar = 0; iVar < 3; iVar++)
				{
					nooc = nooc_f + HexaCode[iVar];
					const int ino_nbr = noocv_v[nooc];
					if (no_flag[ino_nbr] != ino)
					{
						no_flag[ino_nbr] = ino;
						if (iter == 0)
						{
							nbono_i[ino+1] += 1;
						}
						else
						{
							nbono_v[nbono_i[ino]] = ino_nbr;
							nbono_i[ino] += 1;
						}
					}
				}
				
				
			}
		}
		if (iter == 0)
		{
			nbono_i[0] = 0;
			for (int ino = 0; ino < nno; ino++)
        nbono_i[ino+1] += nbono_i[ino];
			nbono_v = new int[nbono_i[nno]];
		}
		else
		{
			for (int ino = nno; ino > 0; ino--)
        nbono_i[ino] = nbono_i[ino-1];
			nbono_i[0] = 0;
		}
	}
	
	delete[] cvono_i_tmp;
	delete[] cvono_v_tmp;
	
	
// Create the edge structure	
	int ino, jno, iEdge, jNode, iNode, iVar, TestEdge = 0;
	
	int **Edge;
	Edge = new  int*[nno];
	for (ino = 0; ino < nno; ino++)
		Edge[ino] = new int [50];
	
	for(ino = 0; ino < nno; ino++)
		for(iVar = 0; iVar < 50; iVar++)
			Edge[ino][iVar] = -1;

	int nEdge = 0;
	for(ino = 0; ino < nno; ino++)
		for (iNode = nbono_i[ino]; iNode <= nbono_i[ino+1]-1; iNode++) {
			jno = nbono_v[iNode];
			
			for (jNode = nbono_i[jno]; jNode <= nbono_i[jno+1]-1; jNode++)
				if ((nbono_v[jNode] == ino) && (jno != ino)) {
					TestEdge = Edge[jno][jNode - nbono_i[jno]];
					break;
				}
			
			if (TestEdge == -1) {
				Edge[ino][iNode - nbono_i[ino]] = nEdge;
				Edge[jno][jNode - nbono_i[jno]] = nEdge;
				nEdge++;
			}
		}

	int **EdgeCV;
	EdgeCV = new  int*[nEdge];
	for (int iEdge = 0; iEdge < nEdge; iEdge++)
		EdgeCV[iEdge] = new int [2];
	
	for(ino = 0; ino < nno; ino++)
		for (iNode = nbono_i[ino]; iNode <= nbono_i[ino+1]-1; iNode++) {
			jno = nbono_v[iNode];
			iEdge = FindEdgeCV(ino, jno, Edge, nbono_i, nbono_v);
			if (ino < jno) {
				EdgeCV[iEdge][0] = ino;
				EdgeCV[iEdge][1] = jno;
			}
		}

	// Create a edge structure method for the hexa grid.
	int *HexaAdaptCode;
	int **HexaEdgeIndex;
	bool **HexaEdgeCode;
	int **HexaEdgeNode;
	int **HexaFaceIndex;
	bool **HexaFaceCode;
	int **HexaFaceNode;
	int **HexaElemIndex;
	bool **HexaElemCode;
	int **HexaElemNode;
	
	HexaAdaptCode = new int[ncv];
	HexaEdgeIndex = new int*[ncv];
	HexaEdgeCode = new bool*[ncv];
	HexaEdgeNode = new int*[ncv];
	HexaFaceIndex = new int*[ncv];
	HexaFaceCode = new bool*[ncv];
	HexaFaceNode = new int*[ncv];
	HexaElemIndex = new int*[ncv];
	HexaElemCode = new bool*[ncv];
	HexaElemNode = new int*[ncv];

	for (int icv = 0; icv < ncv; icv++) {
		HexaEdgeIndex[icv] = new int [12];
		HexaEdgeCode[icv] = new bool [12];
		HexaEdgeNode[icv] = new int [12];
		HexaFaceIndex[icv] = new int [6];
		HexaFaceCode[icv] = new bool [6];
		HexaFaceNode[icv] = new int [6];
		HexaElemIndex[icv] = new int [1];
		HexaElemCode[icv] = new bool [1];
		HexaElemNode[icv] = new int [1];
	}
	
	
	
	for (int icv = 0; icv < ncv; icv ++) {
		int noc_f = noocv_i[icv];
		
		int ip_0 = noocv_v[noc_f];		int ip_1 = noocv_v[noc_f+1];
		int ip_2 = noocv_v[noc_f+2];	int ip_3 = noocv_v[noc_f+3];
		int ip_4 = noocv_v[noc_f+4];	int ip_5 = noocv_v[noc_f+5];
		int ip_6 = noocv_v[noc_f+6];	int ip_7 = noocv_v[noc_f+7];			
			
		HexaEdgeIndex[icv][0] = FindEdgeCV(ip_0,ip_1, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][0] = false;		HexaEdgeNode[icv][0] = -1;
		HexaEdgeIndex[icv][1] = FindEdgeCV(ip_1,ip_2, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][1] = false;		HexaEdgeNode[icv][1] = -1;
		HexaEdgeIndex[icv][2] = FindEdgeCV(ip_2,ip_3, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][2] = false;		HexaEdgeNode[icv][2] = -1;
		HexaEdgeIndex[icv][3] = FindEdgeCV(ip_3,ip_0, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][3] = false;		HexaEdgeNode[icv][3] = -1;
		HexaEdgeIndex[icv][4] = FindEdgeCV(ip_4,ip_5, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][4] = false;		HexaEdgeNode[icv][4] = -1;
		HexaEdgeIndex[icv][5] = FindEdgeCV(ip_5,ip_6, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][5] = false;		HexaEdgeNode[icv][5] = -1;
		HexaEdgeIndex[icv][6] = FindEdgeCV(ip_6,ip_7, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][6] = false;		HexaEdgeNode[icv][6] = -1;
		HexaEdgeIndex[icv][7] = FindEdgeCV(ip_7,ip_4, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][7] = false;		HexaEdgeNode[icv][7] = -1;
		HexaEdgeIndex[icv][8] = FindEdgeCV(ip_0,ip_4, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][8] = false;		HexaEdgeNode[icv][8] = -1;
		HexaEdgeIndex[icv][9] = FindEdgeCV(ip_1,ip_5, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][9] = false;		HexaEdgeNode[icv][9] = -1;
		HexaEdgeIndex[icv][10] = FindEdgeCV(ip_2,ip_6, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][10] = false;	HexaEdgeNode[icv][10] = -1;	
		HexaEdgeIndex[icv][11] = FindEdgeCV(ip_3,ip_7, Edge, nbono_i, nbono_v);		HexaEdgeCode[icv][11] = false;	HexaEdgeNode[icv][11] = -1;
		
		HexaFaceIndex[icv][0] = FindFaceCV(icv, ip_0, ip_1, ip_2, ip_3);		HexaFaceCode[icv][0] = false;		HexaFaceNode[icv][0] = -1;
		HexaFaceIndex[icv][1] = FindFaceCV(icv, ip_4, ip_5, ip_6, ip_7);		HexaFaceCode[icv][1] = false;		HexaFaceNode[icv][1] = -1;
		HexaFaceIndex[icv][2] = FindFaceCV(icv, ip_1, ip_2, ip_6, ip_5);		HexaFaceCode[icv][2] = false;		HexaFaceNode[icv][2] = -1;
		HexaFaceIndex[icv][3] = FindFaceCV(icv, ip_3, ip_2, ip_6, ip_7);		HexaFaceCode[icv][3] = false;		HexaFaceNode[icv][3] = -1;
		HexaFaceIndex[icv][4] = FindFaceCV(icv, ip_0, ip_3, ip_7, ip_4);		HexaFaceCode[icv][4] = false;		HexaFaceNode[icv][4] = -1;
		HexaFaceIndex[icv][5] = FindFaceCV(icv, ip_0, ip_1, ip_5, ip_4);		HexaFaceCode[icv][5] = false;		HexaFaceNode[icv][5] = -1;

		HexaElemIndex[icv][0] = icv; HexaElemCode[icv][0] = false;	HexaElemNode[icv][0] = -1;

	}
	
	// Initial edges that are going to be divided
	bool DivEdge[nEdge]; bool DivElem[ncv];
	for (int iEdge = 0; iEdge < nEdge; iEdge ++) DivEdge[iEdge] = false;
	for (int icv = 0; icv < ncv; icv ++) DivElem[icv] = false;
	for (int icv = 1; icv < 2; icv ++) DivElem[icv] = true;
	
	
	for (int icv = 0; icv < ncv; icv ++) {
		if (DivElem[icv] == true) {
			for (int iIndex = 0; iIndex < 12; iIndex++) {
				
				int iEdge = HexaEdgeIndex[icv][iIndex];
				int ino = EdgeCV[iEdge][0];
				int jno = EdgeCV[iEdge][1];
				
				double i_no_z = x_no[ino][2];
				double j_no_z = x_no[jno][2];

				if ((i_no_z > 0.005) && (j_no_z > 0.005)) {
					DivEdge[HexaEdgeIndex[icv][iIndex]] = true;
					HexaEdgeCode[icv][iIndex] = true;
					cout << x_no[ino][0] << " " << x_no[ino][1] << endl;
				}
				
				if ((i_no_z < 0.005) && (j_no_z < 0.005)) {
					DivEdge[HexaEdgeIndex[icv][iIndex]] = true;
					HexaEdgeCode[icv][iIndex] = true;
					cout << x_no[ino][0] << " " << x_no[ino][1] << endl;
				}
				
			}
			
		}
	}

	// We must verify that all the cv have the right edges marked (semiadaptation)
	for (int icv = 0; icv < ncv; icv ++) {
		for (int iIndex = 0; iIndex < 12; iIndex++) {
			if (DivEdge[HexaEdgeIndex[icv][iIndex]] == true) {
				HexaEdgeCode[icv][iIndex] = true;
			}
		}
	}
	
	// Only those elements that verify certain rules will be marked for hexa adaptation... 
	// the others will need a new point in the middle and pyrams
	for (int icv = 0; icv < ncv; icv ++) {
		HexaAdaptCode[icv] = CheckHexaCode(HexaEdgeCode[icv]);
		if (HexaAdaptCode[icv] != -1) cout << HexaAdaptCode[icv] <<" "<< icv << endl;
		
		// Set the HexaFaceCode, and HexaElemCode
		if (HexaAdaptCode[icv] == 1) {
			HexaFaceCode[icv][0] = true; HexaFaceCode[icv][1] = true; HexaFaceCode[icv][2] = true; 
			HexaFaceCode[icv][3] = true; HexaFaceCode[icv][4] = true; HexaFaceCode[icv][5] = true;
			HexaElemCode[icv][0] = true;
		}
		if (HexaAdaptCode[icv] == 2) {
			HexaFaceCode[icv][3] = true; HexaFaceCode[icv][5] = true;
		}
		if (HexaAdaptCode[icv] == 3) {
			HexaFaceCode[icv][2] = true; HexaFaceCode[icv][4] = true;
		}
		if (HexaAdaptCode[icv] == 4) {
			HexaFaceCode[icv][0] = true; HexaFaceCode[icv][1] = true;
		}
	}

	// Create the new nodes on the edges, on the faces, and in the element.
	int NodeAtEdges[nEdge]; int NodeAtFaces[nfa]; int NodeAtElem[ncv];
	for (int iEdge = 0; iEdge < nEdge; iEdge++) NodeAtEdges[iEdge] = -1;
	for (int ifa = 0; ifa < nfa; ifa++) NodeAtFaces[ifa] = -1;
	for (int icv = 0; icv < ncv; icv++) NodeAtElem[icv] = -1;

	int nNewNodes = nno;
	
	
	double **NewNodeCoord;
	NewNodeCoord = new double *[10*nno];
	for (int ino=0; ino <10*nno; ino++)
		NewNodeCoord[ino] = new double[3];
	
	
	int no_0, no_1, no_2, no_3, no_4, no_5, no_6, no_7;
	
	for (int icv = 0; icv < ncv; icv ++) {
		int noc_f = noocv_i[icv];
		int ip_0 = noocv_v[noc_f];		int ip_1 = noocv_v[noc_f+1];
		int ip_2 = noocv_v[noc_f+2];	int ip_3 = noocv_v[noc_f+3];
		int ip_4 = noocv_v[noc_f+4];	int ip_5 = noocv_v[noc_f+5];
		int ip_6 = noocv_v[noc_f+6];	int ip_7 = noocv_v[noc_f+7];	
				
		for (int iIndex = 0; iIndex < 12; iIndex++) {

			if (HexaEdgeCode[icv][iIndex] == true) {
				
				if (NodeAtEdges[HexaEdgeIndex[icv][iIndex]] != -1)
					HexaEdgeNode[icv][iIndex] = NodeAtEdges[HexaEdgeIndex[icv][iIndex]];

				if (NodeAtEdges[HexaEdgeIndex[icv][iIndex]] == -1) {
					
					NodeAtEdges[HexaEdgeIndex[icv][iIndex]] = nNewNodes;
					HexaEdgeNode[icv][iIndex] = nNewNodes;
					// Compute the coordinates of the new node					
					if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1;}
					if (iIndex == 1) {no_0 = ip_1; no_1 = ip_2;}
					if (iIndex == 2) {no_0 = ip_2; no_1 = ip_3;}
					if (iIndex == 3) {no_0 = ip_3; no_1 = ip_0;}
					if (iIndex == 4) {no_0 = ip_4; no_1 = ip_5;}
					if (iIndex == 5) {no_0 = ip_5; no_1 = ip_6;}
					if (iIndex == 6) {no_0 = ip_6; no_1 = ip_7;}
					if (iIndex == 7) {no_0 = ip_7; no_1 = ip_4;}
					if (iIndex == 8) {no_0 = ip_0; no_1 = ip_4;}
					if (iIndex == 9) {no_0 = ip_1; no_1 = ip_5;}
					if (iIndex == 10) {no_0 = ip_2; no_1 = ip_6;}
					if (iIndex == 11) {no_0 = ip_3; no_1 = ip_7;}
					NewNodeCoord[nNewNodes][0] = 0.5*(x_no[no_0][0]+x_no[no_1][0]);
					NewNodeCoord[nNewNodes][1] = 0.5*(x_no[no_0][1]+x_no[no_1][1]);
					NewNodeCoord[nNewNodes][2] = 0.5*(x_no[no_0][2]+x_no[no_1][2]);
					nNewNodes++;
				}
			}
		}
		
		for (int iIndex = 0; iIndex < 6; iIndex++) {
			if (HexaFaceCode[icv][iIndex] == true) {
				
				if (NodeAtFaces[HexaFaceIndex[icv][iIndex]] != -1)
					HexaFaceNode[icv][iIndex] = NodeAtFaces[HexaFaceIndex[icv][iIndex]];

				if (NodeAtFaces[HexaFaceIndex[icv][iIndex]] == -1) {
					NodeAtFaces[HexaFaceIndex[icv][iIndex]] = nNewNodes;
					HexaFaceNode[icv][iIndex] = nNewNodes;
					// Compute the coordinates of the new node
					if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_2; no_3 = ip_3;}
					if (iIndex == 1) {no_0 = ip_4; no_1 = ip_5; no_2 = ip_6; no_3 = ip_7;}
					if (iIndex == 2) {no_0 = ip_1; no_1 = ip_2; no_2 = ip_6; no_3 = ip_5;}
					if (iIndex == 3) {no_0 = ip_3; no_1 = ip_2; no_2 = ip_6; no_3 = ip_7;}
					if (iIndex == 4) {no_0 = ip_0; no_1 = ip_3; no_2 = ip_7; no_3 = ip_4;}
					if (iIndex == 5) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_5; no_3 = ip_4;}
					NewNodeCoord[nNewNodes][0] = 0.25*(x_no[no_0][0]+x_no[no_1][0]+x_no[no_2][0]+x_no[no_3][0]);
					NewNodeCoord[nNewNodes][1] = 0.25*(x_no[no_0][1]+x_no[no_1][1]+x_no[no_2][1]+x_no[no_3][1]);
					NewNodeCoord[nNewNodes][2] = 0.25*(x_no[no_0][2]+x_no[no_1][2]+x_no[no_2][2]+x_no[no_3][2]);
					nNewNodes++;
					
				}
			}
		}
		
		for (int iIndex = 0; iIndex < 1; iIndex++) {
			if (HexaElemCode[icv][iIndex] == true) {
				if (NodeAtElem[HexaElemIndex[icv][iIndex]] != -1) 
					HexaElemNode[icv][iIndex] = NodeAtElem[HexaElemIndex[icv][iIndex]];
					
				if (NodeAtElem[HexaElemIndex[icv][iIndex]] == -1) {
					NodeAtElem[HexaElemIndex[icv][iIndex]] = nNewNodes;
					HexaElemNode[icv][iIndex] = nNewNodes;
					// Compute the coordinates of the new node
					if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_2; no_3 = ip_3; no_4 = ip_4; no_5 = ip_5; no_6 = ip_6; no_7 = ip_7;}
					NewNodeCoord[nNewNodes][0] = 0.125*(x_no[no_0][0]+x_no[no_1][0]+x_no[no_2][0]+x_no[no_3][0]+x_no[no_4][0]+x_no[no_5][0]+x_no[no_6][0]+x_no[no_7][0]);
					NewNodeCoord[nNewNodes][1] = 0.125*(x_no[no_0][1]+x_no[no_1][1]+x_no[no_2][1]+x_no[no_3][1]+x_no[no_4][1]+x_no[no_5][1]+x_no[no_6][1]+x_no[no_7][1]);
					NewNodeCoord[nNewNodes][2] = 0.125*(x_no[no_0][2]+x_no[no_1][2]+x_no[no_2][2]+x_no[no_3][2]+x_no[no_4][2]+x_no[no_5][2]+x_no[no_6][2]+x_no[no_7][2]);
					nNewNodes++;
				}
			}
		}
	}
	
	// if Hexa adapt code equals 0, then a semidivision is applied
	int nSemiDivided = 0;
	for (int icv = 0; icv < ncv; icv ++) {
		if (HexaAdaptCode[icv] == 0)
			nSemiDivided++;
	}
	
	// If semidivision, then divide add a new point, divide the hexa into pyramids, 
	// and find the right combination, it also create the new node.
	int nPyram = nSemiDivided*6;
	int nPyramNode = nSemiDivided;
		
	
	int *PyramAdaptCode;
	int **PyramNode;	
	int **HexaPyramIndex;
	int *PyramHexaIndex;
	
	int **PyramEdgeIndex; 
	bool **PyramEdgeCode;
	int **PyramEdgeNode; 
	
	int **PyramFaceNode;
	
	int **PyramElemNode; 
	
	PyramAdaptCode = new int [nPyram];
	PyramNode = new int *[nPyram];	
	HexaPyramIndex = new int *[ncv];	
	PyramHexaIndex = new int [nPyram];	
	PyramEdgeIndex = new int *[nPyram]; 
	PyramEdgeCode = new bool *[nPyram];
	PyramEdgeNode = new int *[nPyram]; 
	PyramFaceNode = new int *[nPyram];
	PyramElemNode = new int *[nPyram]; 
	
	for (int iPyram = 0; iPyram < nPyram; iPyram++) {
		
		PyramNode[iPyram] = new int [4];
		PyramEdgeIndex[iPyram] = new int [4];
		PyramEdgeCode[iPyram] = new bool [4];
		PyramEdgeNode[iPyram] = new int [4]; 
		PyramFaceNode[iPyram] = new int [1];
		PyramElemNode[iPyram] = new int [1];
		
	}
	
	for (int icv = 0; icv < ncv; icv++) {
		HexaPyramIndex[icv] = new int [6];
	}
	
	nPyram = 0;
	for (int icv = 0; icv < ncv; icv ++) {
		if (HexaAdaptCode[icv] == 0) {
			
			// Write the edge combination on the base.
			int noc_f = noocv_i[icv];
		
			int ip_0 = noocv_v[noc_f];		int ip_1 = noocv_v[noc_f+1];
			int ip_2 = noocv_v[noc_f+2];	int ip_3 = noocv_v[noc_f+3];
			int ip_4 = noocv_v[noc_f+4];	int ip_5 = noocv_v[noc_f+5];
			int ip_6 = noocv_v[noc_f+6];	int ip_7 = noocv_v[noc_f+7];
			
			NewNodeCoord[nNewNodes][0] = (1.0/8.0)*(x_no[ip_0][0] + x_no[ip_1][0] + x_no[ip_2][0] + x_no[ip_3][0] + x_no[ip_4][0] + x_no[ip_5][0] + x_no[ip_6][0] + x_no[ip_7][0]);
			NewNodeCoord[nNewNodes][1] = (1.0/8.0)*(x_no[ip_0][1] + x_no[ip_1][1] + x_no[ip_2][1] + x_no[ip_3][1] + x_no[ip_4][1] + x_no[ip_5][1] + x_no[ip_6][1] + x_no[ip_7][1]);
			NewNodeCoord[nNewNodes][2] = (1.0/8.0)*(x_no[ip_0][2] + x_no[ip_1][2] + x_no[ip_2][2] + x_no[ip_3][2] + x_no[ip_4][2] + x_no[ip_5][2] + x_no[ip_6][2] + x_no[ip_7][2]);
			
			// Create the 1st pyramid.			
			HexaPyramIndex[icv][0] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_1;
			PyramNode[nPyram][2] = ip_2; PyramNode[nPyram][3] = ip_3;

			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][0]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][1];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][2]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][3];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][0]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][1];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][2]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][3];
		
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;
			
			// Create the 2nd pyramid.			
			HexaPyramIndex[icv][1] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_4; PyramNode[nPyram][1] = ip_7;
			PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_5;
			
			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][7]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][6];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][5]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][4];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][7]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][6];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][5]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][4];
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;
			
			// Create the 3th pyramid.			
			HexaPyramIndex[icv][2] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_4;
			PyramNode[nPyram][2] = ip_5; PyramNode[nPyram][3] = ip_1;
			
			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][8]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][4];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][9]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][0];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][8]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][4];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][9]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][0];
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;
			
			// Create the 4th pyramid.		
			HexaPyramIndex[icv][3] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_3; PyramNode[nPyram][1] = ip_2;
			PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_7;
			
			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][2]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][10];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][6]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][11];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][2]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][10];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][6]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][11];
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;
			
			// Create the 5th pyramid.			
			HexaPyramIndex[icv][4] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_1; PyramNode[nPyram][1] = ip_5;
			PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_2;
			
			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][9]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][5];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][10]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][1];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][9]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][5];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][10]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][1];
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;

			// Create the 6th pyramid.	
			HexaPyramIndex[icv][5] = nPyram; PyramHexaIndex[nPyram] = icv; PyramElemNode[nPyram][0] = nNewNodes; 
			
			PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_3;
			PyramNode[nPyram][2] = ip_7; PyramNode[nPyram][3] = ip_4;
			
			PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[icv][3]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[icv][11];
			PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[icv][7]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[icv][8];
			
			PyramEdgeNode[nPyram][0] = HexaEdgeNode[icv][3]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[icv][11];
			PyramEdgeNode[nPyram][2] = HexaEdgeNode[icv][7]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[icv][8];
			
			for (int iIndex = 0; iIndex < 4; iIndex++)
				if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
			nPyram++;
			
			nNewNodes++; 
	
		}
	}
	
	
//	Check the kind of Pyram partitioning that should be applied
	for (int iPyram = 0; iPyram < nPyram; iPyram ++) {
		PyramAdaptCode[iPyram] = CheckPyramCode(PyramEdgeCode[iPyram]);
		if (PyramAdaptCode[iPyram] == 0) cout << "There is a problem with one Pyram" << endl;
	}
		
	// Find the node that is in the base of the Pyram
	for (int iPyram = 0; iPyram < nPyram; iPyram ++) {
		int icv = PyramHexaIndex[iPyram];
		int ifa = FindFaceCV(icv, PyramNode[iPyram][0], PyramNode[iPyram][1], PyramNode[iPyram][2], PyramNode[iPyram][3]);
		PyramFaceNode[iPyram][0] = NodeAtFaces[ifa];
	}
		
	// Plot adapted solution using tecplot
	ofstream tecplot_file_;
	tecplot_file_.open("adapted.plt", ios::out);
	tecplot_file_ << " VARIABLES = \"x\",\"y\",\"z\" " << endl;
	
	int nNewcv = ncv;
	for (int icv = 0; icv < ncv; icv ++) {
		if (HexaAdaptCode[icv] == 1) nNewcv = nNewcv + 7;
		if (HexaAdaptCode[icv] == 2) nNewcv = nNewcv + 3;
		if (HexaAdaptCode[icv] == 3) nNewcv = nNewcv + 3;
		if (HexaAdaptCode[icv] == 4) nNewcv = nNewcv + 3;
		if (HexaAdaptCode[icv] == 5) nNewcv = nNewcv + 1;
		if (HexaAdaptCode[icv] == 6) nNewcv = nNewcv + 1;
		if (HexaAdaptCode[icv] == 7) nNewcv = nNewcv + 1;
		if (HexaAdaptCode[icv] == 0) {
			int iPyram;
			for (int iIndex = 0; iIndex < 6; iIndex++) {
				iPyram = HexaPyramIndex[icv][iIndex];
				if (PyramAdaptCode[iPyram] == 1) nNewcv = nNewcv + 1;
				if (PyramAdaptCode[iPyram] == 2) nNewcv = nNewcv + 3;
				if (PyramAdaptCode[iPyram] == 3) nNewcv = nNewcv + 3;
				if (PyramAdaptCode[iPyram] == 4) nNewcv = nNewcv + 3;
				if (PyramAdaptCode[iPyram] == 5) nNewcv = nNewcv + 3;
				if (PyramAdaptCode[iPyram] == 6) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 7) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 8) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 9) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 10) nNewcv = nNewcv + 2;
				if (PyramAdaptCode[iPyram] == 11) nNewcv = nNewcv + 2;
				if (PyramAdaptCode[iPyram] == 12) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 13) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 14) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 15) nNewcv = nNewcv + 4;
				if (PyramAdaptCode[iPyram] == 16) nNewcv = nNewcv + 4;
			}
			nNewcv = nNewcv-1;
		}
	}
		
	tecplot_file_ << " ZONE T = \"Time = 0.0\", N= "<< nNewNodes <<" , E = "<< nNewcv<<" , F = FEPOINT, ET = BRICK"<< endl;
	
	
	for (int ino = 0; ino < nno; ino++)
		tecplot_file_ << x_no[ino][0]<< " " << x_no[ino][1] << " " << x_no[ino][2] << endl;	
	for (int ino = nno; ino < nNewNodes; ino++)
		tecplot_file_ << NewNodeCoord[ino][0]<< " " << NewNodeCoord[ino][1] << " " << NewNodeCoord[ino][2] << endl;	

	for (int icv = 0; icv < ncv; icv++) {
		if (HexaAdaptCode[icv] == -1) {
			int noc_f = noocv_i[icv];
			int noc_l = noocv_i[icv+1]-1;
			
			for (int noc = noc_f; noc <= noc_l; noc++) {
				int ino = noocv_v[noc];
				tecplot_file_ << ino+1;
				if (noc != noc_l) tecplot_file_ << "\t";
				else tecplot_file_<< "\n";
			}
		}
	}
	
	int nodes[27]; int **Division; int nPart;
	
	Division = new int*[100];
	for (int iVar = 0; iVar < 100; iVar++)
		Division[iVar] = new int[100];
	
	for (int icv = 0; icv < ncv; icv ++) {
		
		// Hexa elements...
		if (HexaAdaptCode[icv] > 0) {

			// First the corners
			int noc_f = noocv_i[icv];
			nodes[0] = noocv_v[noc_f];		nodes[1] = noocv_v[noc_f+1];
			nodes[2] = noocv_v[noc_f+2];	nodes[3] = noocv_v[noc_f+3];
			nodes[4] = noocv_v[noc_f+4];	nodes[5] = noocv_v[noc_f+5];
			nodes[6] = noocv_v[noc_f+6];		nodes[7] = noocv_v[noc_f+7];
			
			// Next the points that correspond to the broken edges.
			nodes[8] = HexaEdgeNode[icv][0]; nodes[9] = HexaEdgeNode[icv][1];
			nodes[10] = HexaEdgeNode[icv][2]; nodes[11] = HexaEdgeNode[icv][3];
			nodes[12] = HexaEdgeNode[icv][4]; nodes[13] = HexaEdgeNode[icv][5];
			nodes[14] = HexaEdgeNode[icv][6]; nodes[15] = HexaEdgeNode[icv][7];
			nodes[16] = HexaEdgeNode[icv][8]; nodes[17] = HexaEdgeNode[icv][9];
			nodes[18] = HexaEdgeNode[icv][10]; nodes[19] = HexaEdgeNode[icv][11];
			
			// Next the points that correspond to the faces.
			nodes[20] = HexaFaceNode[icv][0];
			nodes[21] = HexaFaceNode[icv][1];
			nodes[22] = HexaFaceNode[icv][2];
			nodes[23] = HexaFaceNode[icv][3];
			nodes[24] = HexaFaceNode[icv][4];
			nodes[25] = HexaFaceNode[icv][5];
			
			// Next the points that correspond to the element.
			nodes[26] = HexaElemNode[icv][0];
			
			
			HexaDivision(HexaAdaptCode[icv], nodes, Division, &nPart);
			for (int iPart = 0; iPart < nPart; iPart++) {
				for (int iDivision = 1; iDivision < Division[iPart][0]; iDivision++) {
					tecplot_file_ << Division[iPart][iDivision]+1 << "\t";;
				}
				tecplot_file_ << endl;
			}
		}
		
		// Piram elements...
		if (HexaAdaptCode[icv] == 0) {
			for (int iIndex = 0; iIndex < 6; iIndex ++) {
				int iPyram = HexaPyramIndex[icv][iIndex];
				// First the corners
				nodes[0] = PyramNode[iPyram][0];	
				nodes[1] = PyramNode[iPyram][1];	
				nodes[2] = PyramNode[iPyram][2];	
				nodes[3] = PyramNode[iPyram][3];	
				nodes[4] = PyramElemNode[iPyram][0];
				
				
				// Next the points that correspond to the broken edges.
				nodes[5] = PyramEdgeNode[iPyram][0];
				nodes[6] = PyramEdgeNode[iPyram][1];
				nodes[7] = PyramEdgeNode[iPyram][2];
				nodes[8] = PyramEdgeNode[iPyram][3];
				
				// Next the points that correspond to the base face.
				nodes[9] = PyramFaceNode[iPyram][0];

				PyramDivision(PyramAdaptCode[iPyram], nodes, Division, &nPart);
				for (int iPart = 0; iPart < nPart; iPart++) {
					// Tets case
					if (Division[iPart][0] == 5) {
						tecplot_file_ << 
						Division[iPart][1]+1 << "\t"<< Division[iPart][2]+1 << "\t"<< 
						Division[iPart][3]+1 << "\t"<< Division[iPart][3]+1 << "\t"<< 
						Division[iPart][4]+1 << "\t"<< Division[iPart][4]+1 << "\t"<< 
						Division[iPart][4]+1 << "\t"<< Division[iPart][4]+1; 
					
					}
					
					if (Division[iPart][0] == 6) {
						tecplot_file_ << 
						Division[iPart][1]+1 << "\t"<< Division[iPart][2]+1 << "\t"<< 
						Division[iPart][3]+1 << "\t"<< Division[iPart][4]+1 << "\t"<< 
						Division[iPart][5]+1 << "\t"<< Division[iPart][5]+1 << "\t"<< 
						Division[iPart][5]+1 << "\t"<< Division[iPart][5]+1;
					
					}
					
					tecplot_file_ << endl;
				}
			}
		}
	}
}
	
