#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

#include "combModels/CombModel_Base.h"
#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_SteadyFlamelet.h"




// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class PlateLaminar : public JoeWithModels 
{

  int npos, nval;
  double **boundVal1;
public:
  
  double (*vgrad_diag)[3];        ///< velocity gradient
  double (*vgrad_offdiag)[3];        ///< velocity gradient

public:  
  PlateLaminar(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {    
    if (mpi_rank == 0)
      cout << "PlateLaminar()" << endl;
	  
	vgrad_diag = NULL;     registerVector(vgrad_diag,    "vgrad_diag",    CV_DATA);
	vgrad_offdiag = NULL;     registerVector(vgrad_offdiag,    "vgrad_offdiag",    CV_DATA);

  }

  virtual ~PlateLaminar()  {}

  void transformMeshHook()
  {
	if (!checkParam("TRANSFORM_MESH")) return;
	double scale_mesh_X = getDoubleParam("SCALE_MESH_X", "1.0");
	double scale_mesh_Y = getDoubleParam("SCALE_MESH_Y", "1.0");
	for (int ino = 0; ino < getNno(); ++ino)
	{	
		x_no[ino][0] *= scale_mesh_X;
		x_no[ino][1] *= scale_mesh_Y;
	}
		
		
	// clearFlag for wall distance to recompute the wallDist when mesh deformed
	DoubleScalar *wallD = getScalarData("wallDist");
	wallD->clearFlag();
  }
	
  void initialHook() 
  {
    JoeWithModels::initialHook();
    
    double fakt_bl = getDoubleParam("SCALE_BL", "1.0");
    
    FILE *fp;

    // ======================================================
    // read laminar profile
    if ((fp=fopen("laminar_inlet.dat", "rt")) == NULL)
    {
      if (mpi_rank == 0)
        cerr << "could not open laminar_inlet.dat, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(fp, "n=%d\td=%d", &npos, &nval);
    getMem2D(&boundVal1, 0, npos-1, 0, nval-1, "boundVal1");

    // file has values
    // x, y, z,press, rho, u, v, w
    for (int i=0; i<npos; i++)
    {
      for (int v=0; v<nval; v++)
      {
        fscanf(fp, "%lf", &boundVal1[i][v]);
        if (v==1)
	  boundVal1[i][v] *= fakt_bl;
        
	printf("%.6le\t", boundVal1[i][v]);
      }
      printf("\n");
    }

    fclose(fp);

    if (!checkDataFlag(rho))
      for (int icv=0; icv<ncv; icv++)
      {
        int pos=1;
		
        while((boundVal1[pos][1] < x_cv[icv][1]) && (pos < npos-1))      pos++;
        double f = (x_cv[icv][1]-boundVal1[pos-1][1])/(boundVal1[pos][1]-boundVal1[pos-1][1]);

        double rho1 = boundVal1[pos-1][4];
        double rho2 = boundVal1[pos][4];
        rho[icv] = rho1 + f*(rho2-rho1);

        double uvel1 = boundVal1[pos-1][5];
        double uvel2 = boundVal1[pos][5];
        rhou[icv][0] = uvel1*rho1 + f*(uvel2*rho2-uvel1*rho1);
        rhou[icv][1] = rhou[icv][2] = 0.0;

        double press1 = boundVal1[pos-1][3];
        double press2 = boundVal1[pos][3];
        press[icv] = press1 + f*(press2-press1);
		  
        rhoE[icv] = press[icv]/(gamma[icv]-1.0)+ 0.5/rho[icv]*vecDotVec3d(rhou[icv], rhou[icv]);    
      }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho,  REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    
  }  
  
  
  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone)
  {
    if (zone->getNameString() == "Inlet-5")   // if more HOOK boundaries are defined
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=1;
        while((boundVal1[pos][1] < x_fa[ifa][1]) && (pos < npos-1))      pos++;
        double f = (x_fa[ifa][1]-boundVal1[pos-1][1])/(boundVal1[pos][1]-boundVal1[pos-1][1]);
        
        double uvel1 = boundVal1[pos-1][5];
        double uvel2 = boundVal1[pos][5];
        vel[ifa][0] = uvel1 + f*(uvel2-uvel1);
        
        double vvel1 = boundVal1[pos-1][6];
        double vvel2 = boundVal1[pos][6];
        vel[ifa][1] = 0.0;//vvel1 + f*(vvel2-vvel1);
        vel[ifa][2] = 0.0;
        
        double rho1 = boundVal1[pos-1][4];
        double rho2 = boundVal1[pos][4];
        double rho = rho1 + f*(rho2-rho1);
        
        double press1 = boundVal1[pos-1][3];
        double press2 = boundVal1[pos][3];
        press[ifa] = press1 + f*(press2-press1);

		temp[ifa] = press[ifa]/(rho*R_gas);
      }
    }
  }
	
  void temporalHook()
  {    
	
	calcGradVel();
	for (int icv=0; icv<ncv; icv++)
	{
		  vgrad_diag[icv][0] = grad_u[icv][0][0];
		  vgrad_diag[icv][1] = grad_u[icv][1][1];
		  vgrad_diag[icv][2] = grad_u[icv][2][2];
		  vgrad_offdiag[icv][0] = grad_u[icv][0][1];
		  vgrad_offdiag[icv][1] = grad_u[icv][0][2];
		  vgrad_offdiag[icv][2] = grad_u[icv][1][2];
	}
	updateCvData(vgrad_diag, REPLACE_ROTATE_DATA);	
	updateCvData(vgrad_offdiag, REPLACE_ROTATE_DATA);
	
	
/*    if (step%500 == 0)
	{
		writeRampWall(step);
		//writeProfile(step);
	}
*/
  }

  void finalHook() 
  {
	
	  FILE *fp;
    if (turbModel == KOMSST)
		fp = fopen("KOMSST_profile.txt", "wt");
	else if (turbModel == KOM)
		fp = fopen("KOM_profile.txt", "wt");
	else if (turbModel == KEPS)
		fp = fopen("KEps_profile.txt", "wt");
	else if (turbModel == SA)
		fp = fopen("SA_profile.txt", "wt");
	else
	{
		cout << "No Turbulence model was specified, writing to newprofile.txt" << endl;
		fp = fopen("newprofile.txt", "wt");
    }
    

    for (int icv=0; icv<ncv; icv++)
    //if ((x_cv[icv][0] > 0.265) && (x_cv[icv][0] < 0.27)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.385) && (x_cv[icv][0] < 0.39)) // ReTheta2800
    //if ((x_cv[icv][0] > 0.352) && (x_cv[icv][0] < 0.358)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.288) && (x_cv[icv][0] < 0.294)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.32) && (x_cv[icv][0] < 0.324)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.3) && (x_cv[icv][0] < 0.306)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.312) && (x_cv[icv][0] < 0.32)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.306) && (x_cv[icv][0] < 0.312)) // ReTheta2000
    //if ((x_cv[icv][0] > 0.3) && (x_cv[icv][0] < 0.302)) // ReTheta2000
    {
      fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t", x_cv[icv][0], x_cv[icv][1], rho[icv], rhou[icv][0], rhou[icv][1], press[icv]);

      for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
        fprintf(fp, "%.8le\t", data->phi[icv]);
      
      fprintf(fp, "%.8le\t", InterpolateAtCellCenterFromFaceValues(mul_fa, icv));
      
      if (turbModel == KOMSST)
      {
        fprintf(fp, "%.8le\t", grad_u[icv][0][1]);
        
        double *muT = getR1("muT");
        fprintf(fp, "%.8le\t", muT[icv]);
      }

      fprintf(fp, "\n");
    }

    fclose(fp);
  }

  void writeProfile(int step)
  {
	  // write to file
	  FILE *fp;
	  char fname[200];
	  
	  if (turbModel == KOMSST)
		  sprintf(fname, "SST.%06d.txt",step);
	  else if (turbModel == KOM)
		  sprintf(fname, "KOM.%06d.txt",step);
	  else if (turbModel == SA)
		  sprintf(fname, "SA.%06d.txt",step);
	  else
	  {
		  cout << "No Turbulence model was specified, writing to newprofile.txt" << endl;
		  sprintf(fname, "profile.%06d.dat",step);
	  }
	  
	  if ( mpi_rank == 0 )
	  {
		  if ( (fp=fopen(fname,"w"))==NULL )
		  {
			  cerr << "Error: cannot open file " << fname << endl;
			  throw(-1);
		  }
	  }
	  else
	  {
		  int dummy;
		  MPI_Status status;
		  MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
		  if ( (fp=fopen(fname,"a"))==NULL )
		  {
			  cerr << "Error: cannot open file " << fname << endl;
			  throw(-1);
		  }
	  }
		  
	  for (int icv=0; icv<ncv; icv++)
	  {
		  fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t", x_cv[icv][0], x_cv[icv][1], rho[icv], rhou[icv][0], rhou[icv][1], press[icv]);
		  
                  for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
			  fprintf(fp, "%.8le\t", data->phi[icv]);
		  
		  fprintf(fp, "%.8le\t", InterpolateAtCellCenterFromFaceValues(mul_fa, icv));
		  
		  if (turbModel == KOMSST)
		  {
			  fprintf(fp, "%.8le\t", grad_u[icv][0][1]);
			  
			  double *muT = getR1("muT");
			  fprintf(fp, "%.8le\t", muT[icv]);
		  }
		  
		  fprintf(fp, "\n");
	  }
	  
	  fclose(fp);
  }
	
  void writeRampWall(int step)
  {
    ////////////////////////////////////////////////////////////
    // define these values to get dim.-less quantities
    ////////////////////////////////////////////////////////////
    double delta = 6.925e-3;
    double pinf = p_ref;

    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == "wall")
          my_wall_faces = zone->ifa_l - zone->ifa_f + 1;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "rampWall.%06d.dat",step);
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"w"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X/delta\" \"Y\" \"Z\" \"rho\" \"press/pinf\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(fname,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone != faZoneList.end(); faZone++)
      if (faZone->getKind() == FA_ZONE_BOUNDARY)
        if (faZone->getNameString() == "wall")
        {
          for (int ifa = faZone->ifa_f; ifa <= faZone->ifa_l; ifa++)
          {
            int icv = cvofa[ifa][0];

            double n[3], s_half[3], vel[3], velTang[3];

            normVec3d(n, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);

            vel[0] = rhou[icv][0]/rho[icv];
            vel[1] = rhou[icv][1]/rho[icv];
            vel[2] = rhou[icv][2]/rho[icv];
            double un = vecDotVec3d(vel, vel);

            velTang[0] = vel[0] - un*n[0];
            velTang[1] = vel[1] - un*n[1];
            velTang[2] = vel[2] - un*n[2];

            double velMag = sqrt(vecDotVec3d(vel, vel));
            double wallDist = fabs(vecDotVec3d(n, s_half));

            double tau = mul_fa[ifa]*velMag/wallDist*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = wallDist*utau/(mul_fa[ifa]/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mul_fa[ifa]/Pr;
            double wallTemp = 307.0;
            double qDot = kTotal*(temp[icv]-wallTemp)/wallDist;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t",
                (x_fa[ifa][0])/delta, x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv]/pinf, temp[icv], tau, qDot, yplus, mul_fa[ifa]);

            fprintf(fp, "\n");
          }
        }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }
/*
 *  Following virtual functions can be overloaded by user.
 *  Comment: some functions require that a scalar is defined.
 */

//  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {/*empty*/}  
//  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {/*empty*/}
//  void temporalHook() {/*empty*/}  
//  void finalHook() {/*empty*/}
//
//  virtual void UserDefinedScalarClipping(const string &name)  {/*empty*/}
//
//  virtual void initialHookScalarRansTurbModel() {/*empty*/}
//  virtual void calcRansTurbViscMuet(double *rho, double (*rhou)[3]) {/*empty*/}
//  virtual double calcTurbProd(int icv) {/*empty*/}
//  virtual double calcDissipProd(int icv) {/*empty*/}
//  virtual void diffusivityHookScalarRansTurb(const string &name) {/*empty*/}
//  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}
//  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name) {/*empty*/}
//  virtual void sourceHookScalarRansTurbExpl(double *rhs, const string &name) {/*empty*/}
//
//  virtual void initialHookScalarRansCombModel() {/*empty*/} 
//  virtual void sourceHookRansComb(double *rhsRho, double(*rhsRhou)[3], double *rhsRhoE, double(*A)[5][5]) {/*empty*/}
//  virtual void sourceHookScalarRansComb(double *rhs, double *A, const string &name) {/*empty*/}
//  virtual void boundaryHookScalarRansComb(double *phi_ph, FaZone *zone, const string &name) {/*empty*/}  
//  virtual void UserDefinedScalarClipping(const string &name) {/*empty*/}

};


// ###########################################################################################
// ------                                                                               ------
// ###########################################################################################

class PlateTurbSA : public PlateLaminar, public RansTurbSA 
{
public:
  PlateTurbSA(char *name) : PlateLaminar(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)
      cout << "PlateTurbSA()" << endl;
  }

  virtual ~PlateTurbSA()  {}
  
//  void initialHook() {JoeWithModels::initialHook();}    
};


// ###########################################################################################
// ------                                                                               ------
// ###########################################################################################

class PlateTurbSST : public PlateLaminar, public RansTurbKOmSST
{
public:
  PlateTurbSST(char *name) : PlateLaminar(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)
      cout << "PlateTurbSST()" << endl;
  }

  virtual ~PlateTurbSST()  {}
  
//  void initialHook() {JoeWithModels::initialHook();}
};

// ###########################################################################################
// ------                                                                               ------
// ###########################################################################################

class PlateTurbWXKOM : public PlateLaminar, public RansTurbKOm
{
public:
  PlateTurbWXKOM(char *name) : PlateLaminar(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)
      cout << "PlateTurbWXKOM()" << endl;
  }

  virtual ~PlateTurbWXKOM()  {}
  
//  void initialHook() {JoeWithModels::initialHook();}    
  
  virtual void initialHookScalarRansTurbModel() 
  {
    
    RansTurbKOm::initialHookScalarRansTurbModel();
    
/*    double kineInit, omegaInit;
    Param *pmy;
    if (getParam(pmy, "INITIAL_CONDITION_TURB"))
    {
      kineInit = pmy->getDouble(1);
      omegaInit = pmy->getDouble(2);
    }
    else
    {
      cerr << " Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
      throw(-1);
    }

  //  if (!checkScalarFlag("kine"))
      for (int icv=0; icv<ncv; icv++)
        kine[icv] = kineInit;

   // if (!checkScalarFlag("omega"))
      for (int icv=0; icv<ncv; icv++)
        omega[icv] = omegaInit;

    updateCvDataByName("kine", REPLACE_DATA);
    updateCvDataByName("omega", REPLACE_DATA);*/
  }

};



int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
  int run = 0;

  // set input name to default "Joe.in"
  char inputFileName[50];
  sprintf(inputFileName, "Joe.in");

  for (int i=1; i<argc; i++)
  {
    string str(argv[i]);
    if (from_string<int>(run, str, std::dec))
    {
      if (mpi_rank == 0)
        cout << "You have specified run number = " << run << endl;
    }
    else
      strcpy(inputFileName, argv[i]);
  }

  if (mpi_rank == 0)
  {
    cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
    cout << "SPECIFIED RUN = " << run << endl;
  }


  try {

    JoeWithModels *joe;

    switch (run)
    {
    case 0:   joe = new PlateLaminar(inputFileName);    break;
    case 1:   joe = new PlateTurbSA(inputFileName);     break;
    case 2:   joe = new PlateTurbSST(inputFileName);    break;
    case 3:   joe = new PlateTurbWXKOM(inputFileName);    break;
    default: 
      if (mpi_rank == 0)
        cerr << "ERROR: run number not available!" << endl;
    }
    
    
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    joe->run();

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

    delete joe;

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
  return (0);
}



