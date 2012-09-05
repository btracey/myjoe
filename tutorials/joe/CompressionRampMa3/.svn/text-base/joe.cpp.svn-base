#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

#include "combModels/CombModel_Mixing.h"
#include "combModels/CombModel_SteadyFlamelet.h"


// ###########################################################################################
// ------                                                                               ------
// ------                                                                               ------
// ------                                                                               ------
// ###########################################################################################

// define turbmodel here; e.g. RansTurbSA

class CompRamp : public JoeWithModels, public RansTurbSA
{

  // memory of inlet profiles (2 inlet profiles at different ReTheta)
#define NVALUES 6
  double (*boundVal1)[NVALUES];
  double (*boundVal2)[NVALUES];

public:
  CompRamp(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)
      cout << "CompRamp()" << endl;
  }

  virtual ~CompRamp()  {}
  
  
  /**
   * initial hook reads two inlet profiles 
   */
  void initialHook()
  {
    JoeWithModels::initialHook();
        
    FILE *fp;

    // ======================================================    
    // read first profile
    if ((fp=fopen("inletReTheta2000.txt", "rt")) == NULL)
    {
      cout << "could not open inletReTheta2000.txt, apply boundary from input file" << endl;
      throw(-1);
    }

    int nPos;
    fscanf(fp, "n=%d", &nPos);

    boundVal1 = new double[nPos][NVALUES];

    // file has values 
    // y, rho, rhou, rhov, rhoe, nuSA
    for (int i=0; i<nPos; i++)
      for (int v=0; v<NVALUES; v++)
        fscanf(fp, "%lf", &boundVal1[i][v]);
    
    fclose(fp);


    // ======================================================
    // read second profile
    if ((fp=fopen("inletReTheta2800.txt", "rt")) == NULL)
    {
      cout << "could not open inletReTheta2800.txt, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(fp, "n=%d", &nPos);

    boundVal2 = new double[nPos][NVALUES];

    // file has values 
    // y, rho, rhou, rhov, rhoe, nuSA
    for (int i=0; i<nPos; i++)
      for (int v=0; v<NVALUES; v++)
        fscanf(fp, "%lf", &boundVal2[i][v]);
    
    fclose(fp);
  }
  
  double getValue(int pos, int val, double fakt)
  {
    return (boundVal1[pos][val] + fakt*(boundVal2[pos][val]-boundVal1[pos][val]));
  }

  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone)
  {
    if (zone->getNameString() == "inlet")   // if more HOOK boundaries are defined
    {
      double interpol = getDoubleParam("FAKT_INTERPOL_INLET_PROFILES", "0.0");
      
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      { 
        int pos=0;
        while(boundVal1[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal1[pos-1][0])/(boundVal1[pos][0]-boundVal1[pos-1][0]);
  
        double rho1 = getValue(pos-1, 1, interpol);
        double rho2 = getValue(pos,   1, interpol);
        double rho = rho1 + f*(rho2-rho1);
  
        double uvel1 = getValue(pos-1, 2, interpol)/rho1;//boundVal[pos-1][2]/boundVal[pos-1][1];
        double uvel2 = getValue(pos,   2, interpol)/rho2;
        vel[ifa][0] = uvel1 + f*(uvel2-uvel1);
  
        double vvel1 = getValue(pos-1, 3, interpol)/rho1;//boundVal[pos-1][3]/boundVal[pos-1][1];
        double vvel2 = getValue(pos,   3, interpol)/rho2;//boundVal[pos][3]/boundVal[pos][1];
        vel[ifa][1] = vvel1 + f*(vvel2-vvel1);
        vel[ifa][2] = 0.0;

        double rhoe1 = getValue(pos-1, 4, interpol);
        double rhoe2 = getValue(pos,   4, interpol);
        double press1 = (GAMMA-1.0)*(rhoe1-0.5*(uvel1*uvel1+vvel1*vvel1)*rho1);
        double press2 = (GAMMA-1.0)*(rhoe2-0.5*(uvel2*uvel2+vvel2*vvel2)*rho2);
        press[ifa] = press1 + f*(press2-press1);

        temp[ifa] = press[ifa]/(rho*R_gas);
      }
    }

  /*  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      printf("%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n", x_fa[ifa][1], rho[ifa], vel[ifa][0], vel[ifa][1], press[ifa]);*/
  }

  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name) 
  {
    if (zone->getNameString() == "inlet")   // if more HOOK boundaries are defined
    {
      double fakt = getDoubleParam("FAKT_NU_SA", "1.0");
      double interpol = getDoubleParam("FAKT_INTERPOL_INLET_PROFILES", "0.0");

      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int pos=0;  
        while(boundVal1[pos][0] < x_fa[ifa][1])      pos++;
        double f = (x_fa[ifa][1]-boundVal1[pos-1][0])/(boundVal1[pos][0]-boundVal1[pos-1][0]);
  
        double nuSA1 = getValue(pos-1, 5, interpol);
        double nuSA2 = getValue(pos,   5, interpol);
        phi_ph[ifa] = fakt*(nuSA1 + f*(nuSA2-nuSA1));
      }
    }
  }

  void temporalHook()
  {    
    if (step%10 == 0)
    {
      ////////////////////////////////////////////////////////////
      // define these values to get dim.-less quantities  
      ////////////////////////////////////////////////////////////
      double delta = 6.7e-3;
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
      sprintf(fname, "rampWall.dat");
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

              double tau = mul_fa[ifa]*velMag/wallDist;
              double utau = sqrt(fabs(tau)/rho[icv]);
              double yplus = wallDist*utau/(mul_fa[ifa]/rho[icv]);

              double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mul_fa[ifa]/Pr;
              double wallTemp = 307.0;
              double qDot = kTotal*(temp[icv]-wallTemp)/wallDist;

              fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t",
                  x_fa[ifa][0]/delta, x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv]/pinf, temp[icv], tau, qDot, yplus, mul_fa[ifa]);

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
  }
};











////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/// 
//   MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN 
///
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

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

    
    // instantiate JOE
    joe = new CompRamp(inputFileName);


    
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    // run JOE
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



