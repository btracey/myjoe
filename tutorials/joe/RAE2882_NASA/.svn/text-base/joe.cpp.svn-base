#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"


#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_SteadyFlamelet.h"
#include "combModels/CombModel_FPVA.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels
{
public:
  
  class POINT
  {
  public :
    double mx, my, mz;
    
    bool operator==(const POINT &other) 
    {
      if ((mx == other.mx) && (my == other.my) && (mz == other.mz))   return true;
      else return false;
    }
  };
  

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;
    
    init();
  }

  virtual ~MyJoe()  {}
  
  
  void init() 
  {
    // do your own initialization here

    cout << "init RAE()" << endl;

    Remove("farfield");

    double Ma  = getDoubleParam("MACH");

    double ReNum = getDoubleParam("RENUM", "6.5e6");

    if (ReNum != 0.0)    mu_ref = Ma*1.4/ReNum;
    else                 mu_ref = 0.0;

    cout << "MU_REF (mu_ref) changed to: " << mu_ref << endl;

    double aoa = M_PI/180.0*getDoubleParam("AOA");
    double gamma = getDoubleParam("GAMMA");

    double sos = sqrt(gamma*p_ref/rho_ref);
    double speed = Ma*sos;

    double velx = speed*cos(aoa);
    double vely = speed*sin(aoa);

    char name[200];
    sprintf(name, "farfield = CBC %.6lf %.6lf 0.0 %.6lf %.6lf", velx, vely, T_ref, p_ref);
    ParamMap::add(string(name));

    //call init from JoeWithModels afterwards
    JoeWithModels::init();
  }


  void transformMeshHook()
  {

    if (!checkParam("TRANSFORM_MESH")) return;

    double thickness_mean = 0.1211;
    double thickness = getDoubleParam("PROFILETHICKNESS", "0.1211");
    
    // sanity check 
//    if (thickness_mean == thickness) return;
    
    
    // ##################################################
    // get the wall nodes
    for (int ino = 0; ino < getNno(); ++ino) 
      no_flag[ino] = 0;


    double (*coordWall)[3] = NULL;

    // count the wall nodes first 
    // a bit tricky as the profile is a closed loop, 
    // which means that on one cpu you get wallNodes=nfa whereas on multiple cores you get wallNodes=nfa+1 per core 
    int myWallNodes = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
              for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) 
              {
                int ino = noofa_v[nof];
                
                if (no_flag[ino] == 0)
                {
                  myWallNodes++;
                  no_flag[ino] = 1;
                }
              }
      }
    
    // then get the actual coordinates of this mesh-nodes/cpu
    double *myWallNodeCo = new double[myWallNodes*3];
    
    for (int ino = 0; ino < getNno(); ++ino) 
      no_flag[ino] = 0;
    
    int iter = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
              for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) 
              {
                int ino = noofa_v[nof];
  
                if (no_flag[ino] == 0)
                {
                  for (int i=0; i<3; ++i)
                    myWallNodeCo[3*iter+i] = x_no[ino][i];
                  
                  no_flag[ino] = 1;
                  iter++;
                }
              }
      }

    // then do the parallel stuff
    int totWallNodes = 0;
    MPI_Allreduce(&myWallNodes, &totWallNodes, 1, MPI_INT, MPI_SUM, mpi_comm);

    // set up send side
    int send_count[mpi_size], send_displ[mpi_size];
    for (int r=0; r<mpi_size; r++)
    {
      send_count[r] = myWallNodes*3;
      send_displ[r] = 0;
    }

    // set up receive side
    int recv_count[mpi_size], recv_displ[mpi_size];
    int wallcoord = myWallNodes*3;
    MPI_Allgather(&wallcoord, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

    recv_displ[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_displ[i] = recv_count[i-1] + recv_displ[i-1];

    double *totWallNodeCo   = new double[totWallNodes*3];

    MPI_Alltoallv(myWallNodeCo,  send_count, send_displ, MPI_DOUBLE, 
                  totWallNodeCo, recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

    delete [] myWallNodeCo;

    // #####################################################################################################################
    // get rid of the double entries and put the points into a vector<POINT> -> with the vector it is much slower though
    POINT p;
    vector<POINT> profile;
    for (int i = 0; i < totWallNodes; ++i)
    {
      p.mx = totWallNodeCo[3*i];
      p.my = totWallNodeCo[3*i+1];
      p.mz = totWallNodeCo[3*i+2];
      profile.push_back(p);
    }
    
    totWallNodes = 0;
    delete [] totWallNodeCo; totWallNodeCo = NULL;    // delete array

    if (mpi_rank == 0)
      cout << "vector size: " << profile.size() << endl;

    vector<POINT>::iterator end_pos(profile.end());
    vector<POINT>::iterator start_pos(profile.begin());
    POINT value = *start_pos; 
    while (++start_pos != end_pos)
    {
       end_pos = std::remove(start_pos, end_pos, value);
       value = *start_pos;
    }
    profile.erase( end_pos, profile.end() );

    if (mpi_rank == 0)
      cout << "vector size: after " << profile.size() << endl;



    // ######################################################################
    // compute wall displacements 
    // load Camber line first
    FILE *fp;
    if ((fp = fopen("rae2822_coord_camberline_10001.dat", "rt")) == NULL)
    {
      if (mpi_rank == 0)
        cout << "you should provide the file rae2822_coord_camberline_10001.dat" << endl;
      throw(-1);
    }
    int nPoints = 10001;
    vector<POINT> camber;
    for (int i = 0; i < nPoints; ++i) 
    {
      fscanf(fp, "%lf %lf", &p.mx, &p.my);
      camber.push_back(p);
    }
    fclose(fp);
    
    // the compute the displacement
    double *displ = new double[profile.size()];

    for (int i=0; i<profile.size(); ++i)
    {
      int pos=1;
      while((camber[pos].mx < profile[i].mx) && (pos<nPoints-1))      pos++;

      double f = (profile[i].mx-camber[pos-1].mx)/(camber[pos].mx-camber[pos-1].mx);
      double yCamberLine = camber[pos-1].my + f*(camber[pos].my-camber[pos-1].my);
      
      double localThickness = profile[i].my - yCamberLine;
      double newLocalWallCoord = yCamberLine + localThickness*thickness/thickness_mean;
      
      displ[i] = newLocalWallCoord - profile[i].my;

      if (mpi_rank == 0)
        printf("%d\t%d\t%d\t%lf\t%lf\t%lf\n", i, pos, nPoints, displ[i], camber[pos].my, profile[i].mx);
    }

    
    if (mpi_rank == 0)
      cout << " finished i in loop : " << endl;

    
    // ######################################################################
    // deform mesh
    
    double c = getDoubleParam("EXPONENT_FOR_DISTANCEWEIGHT", "-4.0");
    
    int nWallPoints = profile.size();
    double *phi = new double[nWallPoints];

    for (int ino = 0; ino < getNno(); ++ino)
    {
      int iwallNode = -1;
      double phiSum = 0.0;

      for (int i=0; i<nWallPoints; ++i)
      {
        if ((profile[i].mx==x_no[ino][0]) && (profile[i].my==x_no[ino][1]) && (profile[i].mz==x_no[ino][2]))
        {
          iwallNode = i;
        }
        else
        {
          double dist = sqrt(  pow(profile[i].mx-x_no[ino][0], 2.0) 
                             + pow(profile[i].my-x_no[ino][1], 2.0)
                             + pow(profile[i].mz-x_no[ino][2], 2.0));
          phi[i] = pow(dist, c);
          phiSum += phi[i];
        }
      }
      
      if (iwallNode != -1)
        x_no[ino][1] += displ[iwallNode];
      else
      {
        double displMesh = 0.0;

        for (int i=0; i<nWallPoints; ++i)
          displMesh += phi[i]*displ[i]/phiSum;
        
        x_no[ino][1] += displMesh;
      }
   
      if ((mpi_rank == 0) && (ino%1000 == 0))
        cout << " ino in loop : " << ino << endl;
    }
    
    
    delete [] phi;
    delete [] displ;
    
    
    // clearFlag for wall distance to recompute the wallDist when mesh deformed
    DoubleScalar *wallD = getScalarData("wallDist");
    wallD->clearFlag();

  }
  

  void initialHook() {JoeWithModels::initialHook();}  

  
  
  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == name)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "%s.dat", name.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"nx\" \"ny\" \"nz\" \"rho\" \"press\"");
      fprintf(fp, " \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\" \"vel-x\" \"walldist\"\n");
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
        if (faZone->getNameString() == name)
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
            double un = vecDotVec3d(n, vel);

            velTang[0] = vel[0] - un*n[0];
            velTang[1] = vel[1] - un*n[1];
            velTang[2] = vel[2] - un*n[2];

            double velMag = sqrt(vecDotVec3d(velTang, velTang));
            double walld = fabs(vecDotVec3d(s_half, n));
            //walld = sqrt(vecDotVec3d(s_half,s_half));

            double wallTemp = 300.0;
            double mulam = mul_fa[ifa];

            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], fa_normal[ifa][0], fa_normal[ifa][1], fa_normal[ifa][2], 
                rho[icv], press[icv], temp[icv], tau, qDot, yplus, mul_fa[ifa], velMag, walld);
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

  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall");
  }

  virtual void finalHook()
  {
    writeWallValues("wall");
  }
  
/*
 *  Following virtual functions can be overloaded by user.
 *  Comment: some functions require that a scalar is defined.
 */

//  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {/*empty*/}  
//  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {/*empty*/}
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
  

  
  
// ##################################################################
// possible other constructors 
// ##################################################################
  
//  MyJoe(ParamMap &p, int i) : JoeWithModels(p), UgpWithCvCompFlow(p), ownID(i) {
//    if (mpi_rank == 0)      cout << "MyJoe()" << endl;
//    init();
//  }

//  MyJoe(char *name, int i) : JoeWithModels(name), UgpWithCvCompFlow(name), ownID(i) {
//    if (mpi_rank == 0)
//      cout << "MyJoe()" << endl;
//    init();
//  }
  
};


/*
 * MyJoe with Spalart & Allmaras Model
 */
class MyJoeSA : public MyJoe, public RansTurbSA
{
public:
  MyJoeSA(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSA()" << endl;
  }

  virtual ~MyJoeSA() {}
  
  
  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == name)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "%s.dat", name.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"nx\" \"ny\" \"nz\" \"rho\" \"press\" \"temp\"");
      fprintf(fp, " \"tau_wall\" \"qdot\" \"yplus\" \"muLam\" \"vel-x\" \"vel-y\" \"walldist\" \"muT\" \"mutfa\"");
      fprintf(fp, " \"u11\" \"u12\" \"u21\" \"u2ex2\"\n");
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
        if (faZone->getNameString() == name)
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
            double un = vecDotVec3d(n, vel);

            velTang[0] = vel[0] - un*n[0];
            velTang[1] = vel[1] - un*n[1];
            velTang[2] = vel[2] - un*n[2];

            double velMag = sqrt(vecDotVec3d(velTang, velTang));
            double walld = fabs(vecDotVec3d(s_half, n));
            //walld = sqrt(vecDotVec3d(s_half,s_half));

            double wallTemp = 300.0;
            double mulam = calcMuLam(wallTemp);

            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], fa_normal[ifa][0], fa_normal[ifa][1], fa_normal[ifa][2], 
                rho[icv], press[icv], temp[icv], tau, T_bfa[ifa], yplus, mul_fa[ifa], vel[0], vel[1], walld, muT[icv], mut_fa[ifa]);
            
            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\n", grad_u[icv][0][0], grad_u[icv][0][1], grad_u[icv][1][0], grad_u[icv][1][1]);
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

  virtual void temporalHook()
  {
    if (step%10 == 0)
      writeWallValues("wall");
  }

  virtual void finalHook()
  {
    writeWallValues("wall");
  }


};

/*
 * MyJoe with Menter SST Model
 */
class MyJoeSST : public MyJoe, public RansTurbKOmSST
{
public:
  MyJoeSST(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSST()" << endl;
  }

  virtual ~MyJoeSST() {}
};

/*
 * MyJoe with Wolcox k-omega Model
 */
class MyJoeWX : public MyJoe, public RansTurbKOm
{
public:
  MyJoeWX(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWX()" << endl;
  }

  virtual ~MyJoeWX() {}
};






int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
  int run = 1;

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

    // declare pointer to JoeWithModels
    JoeWithModels *joe;
    
    switch (run)
    {
    case 0:   joe = new MyJoe(inputFileName);        break;
    case 1:   joe = new MyJoeSA(inputFileName);      break;
    case 2:   joe = new MyJoeSST(inputFileName);     break;
    case 3:   joe = new MyJoeWX(inputFileName);      break;
    default: 
      if (mpi_rank == 0)
        cerr << "ERROR: run number not available!" << endl;
      throw(-1);
    }
    
    // provide total runtime 
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    // run joe
    joe->run();
    

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

    // delete joe (make sure memory is deallocated in destructors
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



