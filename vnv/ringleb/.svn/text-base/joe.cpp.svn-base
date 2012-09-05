#include "JoeWithModels.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels
{
public:

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe()" << endl;

    init();
    
    strMag = NULL;        registerScalar(strMag, "strMag", CV_DATA);

  }

  virtual ~MyJoe()  {}

  void init() {/*empty*/}

  void initialHook() {
   JoeWithModels::initialHook();
   writeFluentMsh("fluent.msh");
  }  

  void Ringleb(double xx, double yy, double &press, double *velo, double &tempe)
  {
    //This routine spits out pressure, velocity and temperature for a
    //Ringleb Flow for a given spatial location (xx, yy)
    //You can use this to set boundary conditions at the inflow and outflow
    //and to compare with the exact solution

    //Adapted from Masatsuka

    int imax = 500;
    double tol = 1.e-15;
    double Vmin = 0.;
    double Vmax = 2.;
    double b,rho,L,V,psi,Vpsi,theta;        
    double Vp = 0.5*( Vmin + Vmax );

    for(int i=1 ; i<=imax ; i++)
    {
      b  = sqrt(1. - 0.2*Vp*Vp);
      rho = pow(b,5);
 
      L  = 1/b + 1/(3*pow(b,3)) + 1/(5*pow(b,5))
          - 0.5*log((1+b)/(1-b));
      V   = sqrt(1/sqrt(fabs((xx-0.5*L)*(xx-0.5*L)+yy*yy) )/2/rho);

      if (fabs(V-Vp) < tol) break;
  
      if (i == imax) 
        cout<<"did not converge... Sorry."<<endl;

      Vp = V;

    }

    b = sqrt(1. - 0.2*V*V);
    rho = pow(b,5);
    press = rho*b*b/GAMMA;
    L = 1/b + 1/(3*pow(b,3)) + 1/(5*pow(b,5))-0.5*log((1+b)/(1-b));
    psi = sqrt(0.5/V/V-(xx-0.5*L)*rho);
    Vpsi= sqrt(0.5-V*V*(xx-0.5*L)*rho);
 
    if ( fabs( Vpsi - 1. ) < tol  )
         theta = 0.5*3.141592653589793238;
    else
         theta = asin( psi*V );

    velo[0] = -V*cos(theta);
    velo[1] = -V*sin(theta);
    velo[2] = 0.;

    tempe = press/rho;
  }


  virtual void boundaryHook(double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone)
  {
    // Ringleb Flow
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    {
      int icv1 = cvofa[ifa][1];
      Ringleb(x_fa[ifa][0], x_fa[ifa][1], p_fa[icv1], vel_fa[icv1], T_fa[icv1]);
    }
  }
  
  void finalHook()
  {
    double T_exact,vel_exact[3],p_exact,rho_exact;
    double Error=0.,Volume=0.,total_error=0.,total_volume=0.;
  
  
    // Ringleb Flow
    for (int icv = 0; icv <ncv; icv++)
    {
      Ringleb(x_cv[icv][0], x_cv[icv][1], p_exact, vel_exact, T_exact);
      rho_exact = p_exact/T_exact;
      Error+= pow(rho_exact-rho[icv],2)*cv_volume[icv];
      Volume+= cv_volume[icv];
    }
  
    MPI_Allreduce(&Error, &total_error, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&Volume, &total_volume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    
    if (mpi_rank == 0)
    {
      cout << "============================================" << endl;
      cout << " Final L2 Error in density  = XXX YYY " << sqrt(total_error)*96.0 << endl;
      cout << " Total Volume  = " << total_volume << endl;
      cout << "============================================" << endl;
    }
  
  }
};




int main(int argc, char *argv[])
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

    // declare pointer to JoeWithModels
    JoeWithModels *joe;
    
    switch (run)
    {
    case 0:   joe = new MyJoe(inputFileName);        break;



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



