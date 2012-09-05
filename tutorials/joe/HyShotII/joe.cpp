#include "JoeWithModels.h"

//#include "turbModels/TurbModel_KOM.h"
//#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
//#include "turbModels/TurbModel_V2F.h"

//#include "combModels/CombModel_VariableProperties.h"
#include "combModels/CombModel_Mixing.h"
//#include "combModels/CombModel_SteadyFlamelet.h"
#include "combModels/CombModel_FPVA.h"
#include "combModels/CombModel_FPVA_Coeff.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Basic Joe with k-omega SST and:                            ------
// ------                    - transition for (x,y) < (X_TRANSITION,Y_TRANSITION)       ------
// ------                    - writeWallValues to Tecplot                               ------
// ------                    - massflux of Zmean, Zvar and Cmean into and out of domain ------
// ------                                                                               ------
// ###########################################################################################

class JoeKOMbasic : public JoeWithModels, public RansTurbKOmSST
{
public:

  double x_transition, y_transition;

public:
  JoeKOMbasic(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "JoeKOMbasic()" << endl;
  }

  virtual ~JoeKOMbasic() {}

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")  for (int icv=0; icv<ncv; icv++) {
      bool setToLaminar = ( (x_cv[icv][0] < x_transition) && (x_cv[icv][1] < y_transition) ) ;
      double Pk = 0.0;
      if ( ! setToLaminar ) {
        double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0e5);
        Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
        if (SST_limitPK == 1)  Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);
        Pk = max(Pk, 0.0);
      }
      double src = Pk - betaStar*rho[icv]*omega[icv]*kine[icv];
      rhs[icv] += src*cv_volume[icv];
      if (flagImplicit) {
        int noc00 = nbocv_i[icv];
        double dsrcdphi = - betaStar*omega[icv];
        A[noc00] -= dsrcdphi*cv_volume[icv];
      }
    }
    if (name == "omega")  for (int icv=0; icv<ncv; icv++) {
      bool setToLaminar = ( (x_cv[icv][0] < x_transition) && (x_cv[icv][1] < y_transition) ) ;
      double F1 = blendFuncF1[icv];
      double alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
      double beta = F1*beta_1 + (1.0 - F1)*beta_2;
      double Pk = 0.0;
      if ( ! setToLaminar ) {
        double zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
        Pk = strMag[icv]*strMag[icv] - 2./3.*zeta*diverg[icv];
        Pk = max(Pk, 0.0);
      }
      double src = alfa*rho[icv]*Pk + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
      rhs[icv] += src*cv_volume[icv];
      if (flagImplicit) {
        int noc00 = nbocv_i[icv];
        double dsrcdphi =  - 2.0*beta*omega[icv];
        A[noc00] -= dsrcdphi*cv_volume[icv];
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");

    for (int icv = 0; icv < ncv; icv++) {
      bool setToLaminar = ( (x_cv[icv][0] < x_transition) && (x_cv[icv][1] < y_transition) ) ;

      double F1 = blendFuncF1[icv];
      double alfa = F1 * alfa_1 + (1.0 - F1) * alfa_2;
      double beta = F1 * beta_1 + (1.0 - F1) * beta_2;

      double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
      double mut  = min(max(rho[icv] * kine[icv] * zeta, 1.0e-8), 1.0e5);

      double Pk = 0.0 ;
      if ( ! setToLaminar ) {
        Pk = max(mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv], 0.0);
        if (SST_limitPK == 1)  Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);
      }
      rhs[icv][5+kine_Index]  += (Pk - betaStar*rho[icv]*omega[icv]*kine[icv])*cv_volume[icv];

      Pk = 0.0 ;
      if ( ! setToLaminar ) {
	zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
	Pk = max(strMag[icv]*strMag[icv] - 2./3.*zeta*diverg[icv], 0.0);
      }
      rhs[icv][5+omega_Index] += (alfa*rho[icv]*Pk + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv])*cv_volume[icv];

      if (flagImplicit) {
        int noc00 = nbocv_i[icv];
        A[noc00][5+kine_Index][5+kine_Index]   -= - betaStar*omega[icv]*cv_volume[icv];
        A[noc00][5+omega_Index][5+omega_Index] -= - 2.0*beta*omega[icv]*cv_volume[icv];
     }
    }
  }

  virtual void initialHook()
  {
    JoeWithModels::initialHook() ;
    //  Read transition info
    x_transition = getDoubleParam("X_TRANSITION", "-1.0e8");
    y_transition = getDoubleParam("Y_TRANSITION", "-1.0e8");
    if ( ( x_transition > -1.0e7 ) && ( y_transition < -1.0e7 ) )  y_transition=1.0e8 ;
    if ( ( y_transition > -1.0e7 ) && ( x_transition < -1.0e7 ) )  x_transition=1.0e8 ;
    if (mpi_rank == 0)
      cout << "JoeKOMbasic::initialHook: setToLaminar for (x,y)<(" << x_transition << "," << y_transition << ")" << endl;
  }
  
  virtual void writeWallValues(string zonename, char *fname)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == zonename)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    if ( mpi_rank == 0 ) {
      if ( (fp=fopen(fname,"wt"))==NULL ) {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"rho\" \"press\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\" \"ZMean\" \"ZVar\" \"CMean\" \n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else {
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
    if (faZone->getNameString() == zonename)
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

        double velMag   = sqrt(vecDotVec3d(velTang, velTang));
        double walld    = fabs(vecDotVec3d(s_half, n));

        double wallTemp = T_bfa[ifa];
        double mulam    = mul_fa[ifa];

        double tau      = mulam * velMag / walld * sign(vel[0]);
        double utau     = sqrt(fabs(tau) / rho[icv]);
        double yplus    = walld * utau / (mulam / rho[icv]);

        double kTotal   = RoM[icv] * gamma[icv] / (gamma[icv] - 1.0) * lamOcp_fa[ifa];
        double qDot     = kTotal * (temp[icv] - wallTemp) / walld;

        double ZM       = getZMean(ifa);
        double ZV       = getZVar(ifa);
        double CM       = getCMean(ifa);

        fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
            x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv], wallTemp, tau, qDot, yplus, mulam, ZM, ZV, CM);
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

  virtual double getZMean(int ifa) { return 0.0; }
  virtual double getZVar(int ifa)  { return 0.0; }
  virtual double getCMean(int ifa) { return 0.0; }

  virtual void temporalHook()
  {
    //  Write Tecplot files for all walls
    if (true) {
      int interval = 10000000 ;
      Param *param;
      if ( getParam(param,"WRITE_DATA") ) {
        int i = 1;
        while ( i < param->getSize() )
          if ( param->getString(i++) == "INTERVAL" )
            interval = param->getInt(i++);
      }
      if ( step%interval == 0 ) {
        //  Find all walls, write tecplot file for each...
        for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zoneIsWall(zone->getName())) {
          char fname[200] ;  sprintf(fname,"%s.%06d.dat",zone->getName(),step) ;
          writeWallValues(zone->getName(),fname) ;
        }
      }
    }
  }


  virtual void finalHook()
  {
    if (true) {
      //  Find all walls, write tecplot file for each...
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      if (zoneIsWall(zone->getName())) {
        char fname[200] ;  sprintf(fname,"%s.%06d.dat",zone->getName(),step) ;
        writeWallValues(zone->getName(),fname) ;
      }
    }
  }

};


// ###########################################################################################
// ------                                                                               ------
// ###########################################################################################

class JoeKOM_mixing : public JoeKOMbasic, public RansCombMixing
{
public:

  JoeKOM_mixing(char *name) : JoeKOMbasic(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0) cout << "JoeKOM_mixing()" << endl;
  }

  virtual ~JoeKOM_mixing()  {}

  void initialHook()
  {
    JoeKOMbasic::initialHook();
    if (mpi_rank == 0)  cout << "JoeKOM_mixing::initialHook()"<< endl;

    if (!checkScalarFlag("ZMean")) {   
      if (mpi_rank==0)  cout << "Set ZMean=0 everywhere" << endl ;
      for (int icv = 0; icv < ncv; icv++)
        ZMean[icv] = 0.0;
      updateCvDataByName("ZMean", REPLACE_DATA);
    }
  }

  virtual double getZMean(int ifa) 
  { 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");          
    return  eq->phi_bfa[ifa];
  }

};


// ###########################################################################################
// ------                                                                               ------
// ###########################################################################################

//class JoeKOM_fpva_coeff : public JoeKOMbasic, public RansCombFPVA_Coeff<ChemtableAdaptiveLinear>
class JoeKOM_fpva_coeff : public JoeKOMbasic, public RansCombFPVA_Coeff<ChemtableCartesianLinear>
{
public:
  JoeKOM_fpva_coeff(char *name) : JoeKOMbasic(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0) {
      cout << "JoeKOM_fpva_coeff()" << endl;
#ifdef PRESSURE_SCALING_FPVA_COEFF
      cout << "PRESSURE_SCALING_FPVA_COEFF == " << PRESSURE_SCALING_FPVA_COEFF << endl ;
#endif
    }
  }

  virtual ~JoeKOM_fpva_coeff()  {}

  void initialHook()
  {
    JoeKOMbasic::initialHook();
    if (mpi_rank == 0)  cout << "JoeKOM_fpva_coeff::initialHook()"<< endl;

    if (!checkScalarFlag("ZMean")) {   
      if (mpi_rank==0)  cout << "Set ZMean=0 everywhere" << endl ;
      for (int icv = 0; icv < ncv; icv++)
        ZMean[icv] = 0.0;
      updateCvDataByName("ZMean", REPLACE_DATA);
    }

    if (!checkScalarFlag("ZVar")) {   
      if (mpi_rank==0)  cout << "Set ZVar =0 everywhere" << endl ;
      for (int icv = 0; icv < ncv; icv++)
        ZVar[icv] = 0.0;
      updateCvDataByName("ZVar", REPLACE_DATA);
    }

    bool specifiedCold = (false) ;
    Param *param; if ( getParam(param,"COMBUSTION_REGIME") ) if ( param->getString() == "COLD" ) specifiedCold = (true) ;
    if (specifiedCold) {
      if (mpi_rank==0)  cout << "Set CMean=0 everywhere" << endl ;
      for (int icv = 0; icv < ncv; icv++)  CMean[icv] = 0.0 ;
      updateCvDataByName("CMean", REPLACE_DATA);
    }
    else {
      double CMeanMax = 0.0 ;
      if ( checkScalarFlag("CMean") ) {
        for (int icv = 0; icv < ncv; icv++)  CMeanMax = max( CMean[icv] , CMeanMax ) ;
        double dummy = CMeanMax ;  MPI_Allreduce(&dummy,&CMeanMax,1,MPI_DOUBLE,MPI_MAX, mpi_comm);   
      }
      if ( CMeanMax < 1.0e-5 ) {
        if (mpi_rank==0)  cout << "Set CMean<=0.05 everywhere" << endl ;
        for (int icv = 0; icv < ncv; icv++)  CMean[icv] = myChemTable.Lookup(ZMean[icv],ZVar[icv],0.05,"PROG");
        updateCvDataByName("CMean", REPLACE_DATA);
      }
    }
  }

  virtual double getZMean(int ifa) 
  { 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");          
    return  eq->phi_bfa[ifa];
  }
  virtual double getZVar(int ifa)
  { 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZVar");          
    return  eq->phi_bfa[ifa];
  }
  virtual double getCMean(int ifa)
  { 
    ScalarTranspEq *eq;
    eq = getScalarTransportData("CMean");          
    return  eq->phi_bfa[ifa];
  } 

  virtual void temporalHook()
  {
    JoeKOMbasic::temporalHook() ;

    //  Integrate mass-flux of fuel and progress-variable into and out of domain
    //  Also integrate source term of progress-variable over full domain
    //  Find fraction of reverse and subsonic flow at outlet
    if (true) {
      int interval = getIntParam("CHECK_INTERVAL","1") ;
      if ( step%interval == 0 ) {
        double sum[11] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0} ;
        //  Boundaries
        for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY) {
          //  Decide if boundary is inlet, outlet, or anything else...
          //  Note that outflows are sometimes specified as NEUMANN...
          int bndkind=0 ;
          Param *param;
          if (getParam(param, zone->getName())) {
            if (param->getString() == "HOOK")  bndkind=1 ;
            if (param->getString() == "CBC")  bndkind=1 ;
            if (param->getString() == "CBC_SUBSONIC_INLET")  bndkind=1 ;
            if (param->getString() == "CBC_SUBSONIC_OUTLET")  bndkind=2 ;
            if (param->getString() == "NEUMANN")  bndkind=2 ;
          }
          //  Inlet of any kind -- get ZMean flux
          if (bndkind==1)  for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++) {
            double n[3] ;
            double farea = normVec3d(n, fa_normal[ifa]);
            double un = vecDotVec3d(n, vel_bfa[ifa]);
            sum[0] += - rho_bfa[ifa] * un * farea * getZMean(ifa) ;
            sum[9] += - rho_bfa[ifa] * un * farea ;
          }
          //  Outlet of any kind -- get ZMean, ZVar, CMean fluxes; also get fraction of reverse/subsonic faces
          if (bndkind==2)  for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++) {
            int icv = cvofa[ifa][0];
            double n[3] ;
            double farea = normVec3d(n, fa_normal[ifa]);
            double un = vecDotVec3d(n, vel_bfa[ifa]);
            sum[1] += rho_bfa[ifa] * un * farea * getZMean(ifa) ;
            sum[2] += rho_bfa[ifa] * un * farea * getZVar(ifa) ;
            sum[3] += rho_bfa[ifa] * un * farea * getCMean(ifa) ;
            if ( fabs(un) > 1.0e-3 ) {   //  i.e., actual outlet and not symmetry plane...
              sum[6] += farea ;
              if ( un < 0.0 )   sum[7] += farea ;   //  reverse flow
              if ( (0.0<un) & (un<sos[icv]) )   sum[8] += farea ;   //  subsonic flow
            }
            sum[10] += rho_bfa[ifa] * un * farea ;
          }
        }
        //  Full domain
        double CmeanSourceFactor = 1.0 ;
        for (int icv=0; icv<ncv; icv++) {
#ifdef PRESSURE_SCALING_FPVA_COEFF
          CmeanSourceFactor = pow( press[icv]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
          sum[4] += cv_volume[icv] * rho[icv] * myChemTable.Lookup(ZMean[icv],ZVar[icv],CMean[icv],"SRC_PROG") * CmeanSourceFactor ;
          sum[5] += cv_volume[icv] * rho[icv] * myChemTable.Lookup(ZMean[icv],ZVar[icv],CMean[icv],"HeatRelease") * CmeanSourceFactor ;
        }
        double tmp[11] = {sum[0],sum[1],sum[2],sum[3],sum[4],sum[5],sum[6],sum[7],sum[8],sum[9],sum[10]} ;
        MPI_Allreduce(tmp,sum, 11, MPI_DOUBLE, MPI_SUM, mpi_comm);
        if (mpi_rank==0) printf("Monit: %6d %12.4e   %12.4e %12.4e %12.4e   %12.4e %12.4e   %12.4e %12.4e   %12.4e %12.4e\n", step,sum[0] , sum[1],sum[2],sum[3] , sum[4],sum[5] , sum[7]/sum[6],sum[8]/sum[6] , sum[9],sum[10] ) ;
      }
    }
  }

  //  Integrate source terms over full domain in transverse directions by binning, dump to file
  void integrateSources(double dxbin, double dybin, double dzbin)
  {
    if (mpi_rank==0)  cout << "JoeKOM_fpva_coeff::integrateSources(" << dxbin << "," << dybin << "," << dzbin << ")" << endl ;
    //  Find min/max x,y,z
    double xmin[3] = {  1.0e10 ,  1.0e10 ,  1.0e10 } ;
    double xmax[3] = { -1.0e10 , -1.0e10 , -1.0e10 } ;
    for (int icv=0; icv<ncv; icv++) {
      for (int d=0 ; d<3 ; d++)  xmin[d] = min( x_cv[icv][d] , xmin[d] ) ;
      for (int d=0 ; d<3 ; d++)  xmax[d] = max( x_cv[icv][d] , xmax[d] ) ;
    }
    double tmp[3] = { xmin[0] , xmin[1] , xmin[2] } ;
    MPI_Allreduce(tmp,xmin,3,MPI_DOUBLE,MPI_MIN,mpi_comm) ;
    tmp[0] = xmax[0] ; tmp[1] = xmax[1] ; tmp[2] = xmax[2] ;
    MPI_Allreduce(tmp,xmax,3,MPI_DOUBLE,MPI_MAX,mpi_comm) ;
    //  Set nbins[], adjust xmin[]
    int nbins[3] ;
    double dbin[3] = { dxbin , dybin , dzbin } ;
    for (int d=0 ; d<3 ; d++) {
      nbins[d] = (int) ( ( xmax[d]-xmin[d] + dbin[d] - 1.0e-8 ) / dbin[d] ) ;
      nbins[d] = max(1,nbins[d]) ;
      xmin[d] -= 1.0e-8 ;
      xmax[d] += 1.0e-8 ;
      dbin[d] = ( xmax[d]-xmin[d] ) / (double) nbins[d] ;
    }
    MPI_Bcast(nbins,3,MPI_INT,0,mpi_comm) ;
    MPI_Bcast(xmin,3,MPI_DOUBLE,0,mpi_comm) ;
    MPI_Bcast(dbin,3,MPI_DOUBLE,0,mpi_comm) ;
    //  Allocate integrals
#define  NSOURCES   2
    double *xsum[3][1+NSOURCES] ;
    for (int d=0 ; d<3 ; d++) {
      for (int j=0 ; j < 1+NSOURCES ; j++) {
        xsum[d][j] = new double[ nbins[d] ] ;
        for (int i=0 ; i < nbins[d] ; i++)
          xsum[d][j][i] = 0.0 ;
      }
    }
    //  Integrate by binning
    //if (mpi_rank==0)  cout << "Preference (for integrateSources) == " << Preference << endl ;
    double CmeanSourceFactor = 1.0 ;
    for (int icv=0; icv<ncv; icv++) {
#ifdef PRESSURE_SCALING_FPVA_COEFF
      CmeanSourceFactor = pow( press[icv]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
      for (int d=0 ; d<3 ; d++) {
        int bin=0 ;
        while ( xmin[d]+dbin[d]*(bin+1) < x_cv[icv][d] )  bin++ ;
        bin = max(0,min(nbins[d]-1,bin)) ;
        xsum[d][0][bin] += cv_volume[icv] ;
        xsum[d][1][bin] += cv_volume[icv] * rho[icv] * myChemTable.Lookup(ZMean[icv],ZVar[icv],CMean[icv],"SRC_PROG") * CmeanSourceFactor ;
        xsum[d][2][bin] += cv_volume[icv] * rho[icv] * myChemTable.Lookup(ZMean[icv],ZVar[icv],CMean[icv],"HeatRelease") * CmeanSourceFactor ;
        //  add more stuff here, but increase NSOURCES...
      }
    }
    //  Collect across all ranks
    double *tmpsum = new double[ max(nbins[0],max(nbins[1],nbins[2])) ] ;
    for (int d=0 ; d<3 ; d++) {
      for (int j=0 ; j < 1+NSOURCES ; j++) {
        MPI_Allreduce(xsum[d][j],tmpsum, nbins[d], MPI_DOUBLE, MPI_SUM, mpi_comm);
        for (int i=0 ; i < nbins[d] ; i++) {
          xsum[d][j][i] = 0.0 ;
          if ( fabs(tmpsum[i]) > 1.0e-15 )  xsum[d][j][i] = tmpsum[i] ;
        }
      }
    }
    if (mpi_rank==0) {
      //  Write to screen
      if (false) {
        for (int d=0 ; d<3 ; d++) {
          cout << "Integrated sources in direction " << d << " : xbin, sum(dV) [m^3], sum(rho*omega_H2O*dV) [kg/s], sum(rho*omega_Heat*dV) [W]" << endl ;
          for (int i=0 ; i < nbins[d] ; i++) {
            cout << "   " ;  cout.width(8) ;  cout << xmin[d]+dbin[d]*(i+0.5) ;
            for (int j=0 ; j < 1+NSOURCES ; j++) {
              cout << " \t" ;  cout.width(9) ;  cout << xsum[d][j][i] ;
            }
            cout << endl ;
          }
        }
      }
      //  Dump to file
      char fname[200] ;
      FILE *fp ;
      for (int d=0 ; d<3 ; d++) {
        sprintf(fname, "integratedSources%d.%06d.txt", d, step) ;
        fp = fopen(fname,"w") ;
        if (fp==NULL)  cout << "ERROR when opening file " << fname << endl ;
        fprintf(fp,"%%  Integrated sources in direction %d : xbin, sum(dV) [m^3], sum(rho*omega_H2O*dV) [kg/s], sum(rho*omega_Heat*dV) [W]\n", d) ;
        for (int i=0 ; i < nbins[d] ; i++) {
          fprintf(fp," %15.8g ", xmin[d]+dbin[d]*(i+0.5) ) ;
          for (int j=0 ; j < 1+NSOURCES ; j++)
            fprintf(fp," %15.8g", xsum[d][j][i] ) ;
          fprintf(fp,"\n") ;
        }
        fclose(fp);
      }
    }
    //  Clean up
    delete [] tmpsum ;
    for (int d=0 ; d<3 ; d++) {
      for (int j=0 ; j < 1+NSOURCES ; j++) {
        if ( xsum[d][j] != NULL )
          delete [] xsum[d][j] ;
      }
    }
#undef  NSOURCES
  }

  //  Compute massfluxes across all boundaries
  void computeBoundaryMassfluxes()
  {
#define _FIX_OUTLET_IN_VINCEGRID_    //  Vince's grid has unified outlet (bleed and actual outlet) -- fix by adding extra bnd...

    if (mpi_rank==0)  cout << "JoeKOM_fpva_coeff::computeBoundaryMassfluxes()" << endl ;
    //  Number of "species" : 0=total, 1-3={ZMean,ZVar,CMean}, 4->n+3={all species}
    int nspec = 4 + OutputSpeciesName.size() ;
    //  Number of boundaries
    int nbnds = 0 ;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)  if (zone->getKind() == FA_ZONE_BOUNDARY)
      nbnds++ ;
#ifdef  _FIX_OUTLET_IN_VINCEGRID_
      nbnds++ ;   //  Add extra bnd for bleed outlet...
#endif
    //  Allocate storage
    double **massflux = new double*[nbnds] ;
    for (int bnd=0 ; bnd < nbnds ; bnd++) {
      massflux[bnd] = new double[nspec] ;
      for (int i=0 ; i < nspec ; i++)
        massflux[bnd][i] = 0.0 ;
    }
    //  Loop over all boundaries
    int bnd=0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)  if (zone->getKind() == FA_ZONE_BOUNDARY) {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
        int icv = cvofa[ifa][0];
        double n[3];
        double area = normVec3d(n,fa_normal[ifa]) ;
        double mftot = area * rho_bfa[ifa] * vecDotVec3d(n, vel_bfa[ifa]);
        massflux[bnd][0] += mftot ;
        double ZM = getZMean(ifa);
        double ZV = getZVar(ifa);
        double CM = getCMean(ifa);
        massflux[bnd][1] += mftot * ZM ;
        massflux[bnd][2] += mftot * ZV ;
        massflux[bnd][3] += mftot * CM ;
#ifdef  _FIX_OUTLET_IN_VINCEGRID_
        if ( (zone->getNameString() == "outlet" ) && ( x_fa[ifa][0] < 0.4 ) ) {
          massflux[nbnds-1][0] += mftot ;
          massflux[nbnds-1][1] += mftot * ZM ;
          massflux[nbnds-1][2] += mftot * ZV ;
          massflux[nbnds-1][3] += mftot * CM ;\
        }
#endif
        for (int i=0; i<OutputSpeciesName.size(); i++) {
          string SpeciesName = OutputSpeciesName[i];
          double speciesmassfraction = myChemTable.Lookup(ZM, ZV, CM, SpeciesName);
          if (speciesmassfraction < 1.0e-27)  speciesmassfraction = 0.0;
          massflux[bnd][4+i] += mftot * speciesmassfraction ;
#ifdef  _FIX_OUTLET_IN_VINCEGRID_
          if ( (zone->getNameString() == "outlet" ) && ( x_fa[ifa][0] < 0.4 ) ) {
            massflux[nbnds-1][4+i] += mftot * speciesmassfraction ;
          }
#endif
        }
      }
      bnd++;
    }
    //  Collect across all ranks
    double *tmp = new double[nspec] ;
    for (bnd=0 ; bnd < nbnds ; bnd++) {
      MPI_Allreduce(massflux[bnd],tmp,nspec, MPI_DOUBLE, MPI_SUM, mpi_comm);
      for (int i=0 ; i < nspec ; i++) {
        massflux[bnd][i] = 0.0 ;
        if ( fabs(tmp[i]) > 1.0e-10 )  massflux[bnd][i] = tmp[i] ;
      }
    }
    //  Dump to file
    if (mpi_rank==0) {
      char fname[200] ;
      sprintf(fname, "boundaryMassfluxes.%06d.txt", step) ;
      FILE *fp ;  fp = fopen(fname,"w") ;
      if (fp==NULL)  cout << "ERROR when opening file " << fname << endl ;
      fprintf(fp,"%%  Massfluxes across all boundaries, directed out of domain, [kg/s]\n") ;
      //  Line with boundary names
        bnd=0;
        fprintf(fp,"%%") ;
        for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)  if (zone->getKind() == FA_ZONE_BOUNDARY) {
          fprintf(fp,"\t  %s    ", zone->getNameString().data() ) ;
          bnd++ ;
        }
#ifdef  _FIX_OUTLET_IN_VINCEGRID_
        fprintf(fp,"\t  bleed_outlet") ;
#endif
        fprintf(fp,"\n") ;
      //  Each species
      for (int i=0 ; i < nspec ; i++) {
        bnd=0;
        fprintf(fp,"  ") ;
        for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)  if (zone->getKind() == FA_ZONE_BOUNDARY) {
          fprintf(fp,"\t %11.4e", massflux[bnd][i] ) ;
          bnd++ ;
        }
#ifdef  _FIX_OUTLET_IN_VINCEGRID_
        fprintf(fp,"\t %11.4e", massflux[nbnds-1][i] ) ;
#endif
        if      (i==0)  fprintf(fp,"\t  %%   Total") ;
        else if (i==1)  fprintf(fp,"\t  %%   ZMean") ;
        else if (i==2)  fprintf(fp,"\t  %%   ZVar") ;
        else if (i==3)  fprintf(fp,"\t  %%   CMean") ;
        else            fprintf(fp,"\t  %%   %s", OutputSpeciesName[i-4].data() ) ;
        fprintf(fp,"\n") ;
      }
      fclose(fp);
    }
    //  Clean up
    if ( tmp != NULL ) {
      delete [] tmp ;
      tmp = NULL ;
    }
    if ( massflux != NULL ) {
      for (int bnd=0 ; bnd < nbnds ; bnd++)
        if ( massflux[bnd] != NULL )
          delete [] massflux[bnd] ;
      delete [] massflux ;
      massflux = NULL ;
    }
  }

  //  Extract profiles -- use Vince's ugly method for parallel writing...
  void extractProfiles()
  {
    if (mpi_rank==0)  cout << "JoeKOM_fpva_coeff::extractProfiles()" << endl ;
    double xmin[3] = { -1.0e10 , -1.0e10 , -1.0e10 } ;
    double xmax[3] = {  1.0e10 ,  1.0e10 ,  1.0e10 } ;

    //  y-profiles for inflow to hyshot LES
    xmin[0] = 0.3742 ;  xmax[0] = 0.3772 ;  xmax[2] = 0.0006 ;

    //  Open file at rank 0, stop and wait at all other ranks
    FILE *fp ;
    char fname[200] ;  sprintf(fname,"extractProfiles.txt") ;
    if (mpi_rank==0) {
      fp = fopen(fname,"w") ;
      if (fp==NULL)  cout << "ERROR when opening file " << fname << endl ;
      fprintf(fp,"%%  Dumping all cells with center within bin\n") ;
      fprintf(fp,"%%  Columns  1- 9 :  x,y,z , rho,u,v,w,p,T\n") ;
      fprintf(fp,"%%          10-16 :  enthalphy, gamma, sos, RoM, kine, omega, diverg\n") ;
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(fname,"a"))==NULL )
      {
        cout << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    //  Process cells in bin
    for (int icv=0; icv<ncv; icv++) {
      if ( ( xmin[0] < x_cv[icv][0] ) && ( x_cv[icv][0] < xmax[0] )
      &&   ( xmin[1] < x_cv[icv][1] ) && ( x_cv[icv][1] < xmax[1] )
      &&   ( xmin[2] < x_cv[icv][2] ) && ( x_cv[icv][2] < xmax[2] ) ) {
	fprintf(fp," %15.8g %15.8g %15.8g", x_cv[icv][0], x_cv[icv][1], x_cv[icv][2] ) ;
	fprintf(fp," %15.8g %15.8g %15.8g %15.8g", rho[icv], rhou[icv][0]/rho[icv], rhou[icv][1]/rho[icv], rhou[icv][2]/rho[icv] ) ;
	fprintf(fp," %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g", press[icv], temp[icv], enthalpy[icv], gamma[icv], sos[icv], RoM[icv] ) ;
	fprintf(fp," %15.8g %15.8g %15.8g\n", kine[icv], omega[icv], diverg[icv] ) ;
      }
    }

    //  Close file and pass baton
    fclose(fp) ;
    double dummy[1000] ;  for (int i=1 ; i < 5000 ; i++)  dummy[i] = 1.0 / sqrt( (double) i ) ;
    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }

#undef _FIX_OUTLET_IN_VINCEGRID_
};





// ###########################################################################################
// ------              HyShot (reacting)                                                ------
// ###########################################################################################

class hyshot_fpva_coeff : public JoeKOM_fpva_coeff
{
public:

#define NVAL 7
  int    npos;
  double (*bval)[NVAL];

  hyshot_fpva_coeff(char *name) : JoeKOM_fpva_coeff(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0) {
      cout << "hyshot_fpva_coeff()" << endl;
    }

    // read inlet profile
    // Y press temp kine omega vel-X vel-Y
    // 0 1     2    3    4     5     6
    FILE *fp;
    if ((fp=fopen("inflow.profile", "rt")) == NULL) {
      if (mpi_rank == 0)
        cerr << "could not open inflow.profile, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(fp, "n=%d\n", &npos);
    bval = new double[npos][NVAL];
    for (int i=0; i<npos; i++)
      for (int v=0; v<NVAL; v++)
        fscanf(fp, "%lf", &bval[i][v]);
    fclose(fp);

    if (mpi_rank == 0) {
      printf("n=%d\td=%d\n", npos, NVAL);
      for (int i=0; i<npos; i++) {
        for (int v=0; v<NVAL; v++)
          printf("%.6le\t", bval[i][v]);
        printf("\n");
      }
    }
  }

  ~hyshot_fpva_coeff()
  {
    if ( bval != NULL ) {
      delete [] bval ;
      bval = NULL ;
    }
  }

  double getValues(double ycoord, int ival)
  {
    int pos=1;
    while ((bval[pos][0] < ycoord) && (pos<npos-1))     pos++;
    double f = (ycoord-bval[pos-1][0])/(bval[pos][0]-bval[pos-1][0]);
    f = max(0.0,min(1.0,f)) ;
    return ( bval[pos-1][ival] + f*(bval[pos][ival]-bval[pos-1][ival]) ) ;
  }
    
  // Y press temp kine omega vel-X vel-Y
  // 0 1     2    3    4     5     6
  virtual void boundaryHook(double *T, double (*vel)[3], double *press, FaZone *zone)
  {
    //  Note that we here assume that only a single boundary has been defined as HOOK!!!
    //  Otherwise should do "if" around the boundary name...
    //if (zone->getNameString() == "air_inflow" )   // if more HOOK boundaries are defined
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    {
      T[ifa]      = getValues(x_fa[ifa][1], 2);
      vel[ifa][0] = getValues(x_fa[ifa][1], 5);
      vel[ifa][1] = getValues(x_fa[ifa][1], 6);
      vel[ifa][2] = 0.0;
      press[ifa]  = getValues(x_fa[ifa][1], 1);
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    RansTurbKOmSST::boundaryHookScalarRansTurb(phi_fa, zone, name);

    //  Note that this function is called both for boundaries defined as HOOK and for all walls
    //  Hence we must do the stuff below only if it is NOT a wall
    Param *param;
    if (getParam(param, zone->getName()))
      if (param->getString() == "HOOK") {
        if (name == "kine")
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            phi_fa[ifa] = getValues(x_fa[ifa][1], 3);
        if (name == "omega")
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            phi_fa[ifa] = getValues(x_fa[ifa][1], 4);
      }
  }
 
  virtual void finalHook()
  {
    JoeKOMbasic::finalHook() ;
    computeBoundaryMassfluxes() ;
    computeWallNorms(0.075/8.0/10.0) ;
    integrateSources(0.250/32.0,0.0098/16.0,0.075/8.0/16.0) ;
    //extractProfiles() ;
  }

  virtual void initialHook()
  {
    JoeKOM_fpva_coeff::initialHook() ;
    //  Verify that only a single boundary has been defined as HOOK
    //  Otherwise will get trouble in boundaryHook below...
    if (mpi_rank==0)  cout << "hyshot_fpva_coeff::initialHook: counting number of HOOK boundaries..." << endl ;
    int nhooks=0;
    Param *param;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (getParam(param, zone->getName()))
          if (param->getString() == "HOOK")  nhooks++ ;
    if ( nhooks != 1 ) {
      if (mpi_rank==0)  cout << "ERROR : " << nhooks << " boundaries defined as HOOK..." << endl ;
      exit(-1) ;
    }

    //  Set initial field from inflow profile
    if (!checkScalarFlag("RHO")) {   
      if (mpi_rank==0)  cout << "Setting initial N-S fields from inflow profile..." << endl ;
      for (int icv = 0; icv < ncv; icv++) {
        double Pprof = getValues(x_cv[icv][1], 1);
        double Tprof = getValues(x_cv[icv][1], 2);
        double Uprof = getValues(x_cv[icv][1], 5);
        double Vprof = getValues(x_cv[icv][1], 6);
        double Wprof = 0.0;
        double Hprof, dum1, dum2 ;
        ComputeProperties_T(Hprof, RoM[icv], gamma[icv], dum1, dum2, Tprof, ZMean[icv], ZVar[icv], CMean[icv]);
        rho[icv] = Pprof / (RoM[icv] * Tprof) ;
        rhou[icv][0] = rho[icv] * Uprof ;
        rhou[icv][1] = rho[icv] * Vprof ;
        rhou[icv][2] = rho[icv] * Wprof ;
        rhoE[icv] = Hprof*rho[icv] - Pprof + 0.5*vecDotVec3d(rhou[icv],rhou[icv])/rho[icv] ;
      }
      updateCvData(rhou, REPLACE_ROTATE_DATA);
      updateCvData(rho,  REPLACE_DATA);
      updateCvData(rhoE, REPLACE_DATA);
    }
  }


  //  Compute different norms of wall pressure on lower/upper walls as unstart proxies
  //  To avoid having to give name of wall, do for all walls in x=[0.45,0.66]
  //  Compute norms with powers 1,2,4,6,8
  void computeWallNorms(double dzbin)
  {
#define  NTERMS   7    // area; powers 1,2,4,6,8; max
    const double xbeg = 0.45 ;
    const double xend = 0.66 ;
    if (mpi_rank==0)  cout << "hyshot_fpva_coeff::computeWallNorms(" << dzbin << "), [xbeg,xend]==[" << xbeg << "," << xend << "]" << endl ;
    //  Number of boundaries
    int nbnds = 0 ;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    if (zoneIsWall(zone->getName()))
      nbnds++ ;
    //  Find min/max z
    double zmin = 1.0e10 ;
    double zmax = -1.0e10 ;
    for (int icv=0; icv<ncv; icv++) {
      zmin = min( x_cv[icv][2] , zmin ) ;
      zmax = max( x_cv[icv][2] , zmax ) ;
    }
    double ztmp=zmin ;
    MPI_Allreduce(&ztmp,&zmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm) ;
    ztmp=zmax ;
    MPI_Allreduce(&ztmp,&zmax,1,MPI_DOUBLE,MPI_MAX,mpi_comm) ;
    //  Set nbins, adjust zmin
    int nbins ;
    nbins = (int) ( ( zmax-zmin + dzbin - 2.0e-8 ) / dzbin ) ;
    nbins = max(1,nbins) ;
    zmin -= 1.0e-8 ;
    zmax += 1.0e-8 ;
    MPI_Bcast(&nbins,1,MPI_INT,0,mpi_comm) ;
    MPI_Bcast(&zmin,1,MPI_DOUBLE,0,mpi_comm) ;
    double dbin = ( zmax-zmin ) / (double) nbins ;
    MPI_Bcast(&dbin,1,MPI_DOUBLE,0,mpi_comm) ;
    //  Allocate storage
    double **norm[NTERMS] ;
    for (int n=0 ; n < NTERMS ; n++) {
      norm[n] = new double*[nbnds] ;
      for (int bnd=0 ; bnd < nbnds ; bnd++) {
        norm[n][bnd] = new double[nbins] ;
        for (int i=0 ; i < nbins ; i++)
          norm[n][bnd][i] = 0.0 ;
      }
    }
    //  Loop over all boundaries
    int bnd=0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    if (zoneIsWall(zone->getName())) {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
        int icv = cvofa[ifa][0];
        if ( ( xbeg < x_cv[icv][0] ) && ( x_cv[icv][0] < xend ) ) {
          int bin=0 ;
          while ( zmin+dbin*(bin+1) < x_cv[icv][2] )  bin++ ;
          bin = max(0,min(nbins-1,bin)) ;
          double dummy[3];
          double area = normVec3d(dummy,fa_normal[ifa]) ;
          norm[0][bnd][bin] += area ;
          norm[1][bnd][bin] += area * press[icv] ;
          norm[2][bnd][bin] += area * pow( press[icv] , 2 ) ;
          norm[3][bnd][bin] += area * pow( press[icv] , 4 ) ;
          norm[4][bnd][bin] += area * pow( press[icv] , 6 ) ;
          norm[5][bnd][bin] += area * pow( press[icv] , 8 ) ;
          norm[6][bnd][bin] = max( press[icv] , norm[6][bnd][bin] ) ;
        }
      }
      bnd++;
    }
    //  Collect across all ranks
    double *tmp = new double[nbins] ;
    for (int bnd=0 ; bnd < nbnds ; bnd++) {
      for (int n=0 ; n < 6 ; n++) {
        MPI_Allreduce(norm[n][bnd],tmp, nbins, MPI_DOUBLE, MPI_SUM, mpi_comm);
        for (int i=0 ; i < nbins ; i++)
          norm[n][bnd][i] = tmp[i] ;
      }
      for (int n=6 ; n < NTERMS ; n++) {
        MPI_Allreduce(norm[n][bnd],tmp, nbins, MPI_DOUBLE, MPI_MAX, mpi_comm);
        for (int i=0 ; i < nbins ; i++)
          norm[n][bnd][i] = tmp[i] ;
      }
    }
    //  Dump to file
    if (mpi_rank==0) {
      char fname[200] ;
      FILE *fp ;
      sprintf(fname, "wallNorms.%06d.txt", step) ;
      fp = fopen(fname,"w") ;
      if (fp==NULL)  cout << "ERROR when opening file " << fname << endl ;
      fprintf(fp,"%%  Norms at walls, binned in z-direction : zbin, sum(dA) [m^2], sum(p^{1,2,4,6,8}*dA) [Pa^?*m^2], max(p) [Pa]\n") ;
      fprintf(fp,"%%  Number of walls == %d\n", nbnds) ;
      //  Loop over walls
      int bnd=0;
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      if (zoneIsWall(zone->getName())) {
        fprintf(fp,"%%  Wall == %s\n", zone->getNameString().data() ) ;
        for (int i=0 ; i < nbins ; i++) {
          fprintf(fp," %15.8g ", zmin+dbin*(i+0.5) ) ;
          for (int n=0 ; n < NTERMS ; n++)
            fprintf(fp," %15.8g", norm[n][bnd][i] ) ;
          fprintf(fp,"\n") ;
        }
        bnd++;
      }
      fclose(fp);
    }
    //  Clean up
    delete [] tmp ;
    for (int n=0 ; n < NTERMS ; n++) {
      for (int bnd=0 ; bnd < nbnds ; bnd++)
        if ( norm[n][bnd] != NULL )
          delete [] norm[n][bnd] ;
      delete [] norm[n] ;
      norm[n] = NULL ;
    }
#undef  NTERMS
  }



};




// ###########################################################################################
// ------              HyShot (cold, using mixing class) -- for Francisco's testing...  ------
// ###########################################################################################

class hyshot_mixing : public JoeKOM_mixing
{
public:

#define NVAL 7
  int    npos;
  double (*bval)[NVAL];

  hyshot_mixing(char *name) : JoeKOM_mixing(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0) {
      cout << "hyshot_mixing()" << endl;
    }

    // read inlet profile
    // Y press temp kine omega vel-X vel-Y
    // 0 1     2    3    4     5     6
    FILE *fp;
    if ((fp=fopen("inflow.profile", "rt")) == NULL) {
      if (mpi_rank == 0)
        cerr << "could not open inflow.profile, apply boundary from input file" << endl;
      throw(-1);
    }

    fscanf(fp, "n=%d\n", &npos);
    bval = new double[npos][NVAL];
    for (int i=0; i<npos; i++)
      for (int v=0; v<NVAL; v++)
        fscanf(fp, "%lf", &bval[i][v]);
    fclose(fp);

    if (mpi_rank == 0) {
      printf("n=%d\td=%d\n", npos, NVAL);
      for (int i=0; i<npos; i++) {
        for (int v=0; v<NVAL; v++)
          printf("%.6le\t", bval[i][v]);
        printf("\n");
      }
    }
  }

  ~hyshot_mixing()
  {
    if ( bval != NULL ) {
      delete [] bval ;
      bval = NULL ;
    }
  }

  double getValues(double ycoord, int ival)
  {
    int pos=1;
    while ((bval[pos][0] < ycoord) && (pos<npos-1))     pos++;
    double f = (ycoord-bval[pos-1][0])/(bval[pos][0]-bval[pos-1][0]);
    f = max(0.0,min(1.0,f)) ;
    return ( bval[pos-1][ival] + f*(bval[pos][ival]-bval[pos-1][ival]) ) ;
  }
    
  // Y press temp kine omega vel-X vel-Y
  // 0 1     2    3    4     5     6
  virtual void boundaryHook(double *T, double (*vel)[3], double *press, FaZone *zone)
  {
    //  Note that we here assume that only a single boundary has been defined as HOOK!!!
    //  Otherwise should do "if" around the boundary name...
    //if (zone->getNameString() == "air_inflow" )   // if more HOOK boundaries are defined
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    {
      T[ifa]      = getValues(x_fa[ifa][1], 2);
      vel[ifa][0] = getValues(x_fa[ifa][1], 5);
      vel[ifa][1] = getValues(x_fa[ifa][1], 6);
      vel[ifa][2] = 0.0;
      press[ifa]  = getValues(x_fa[ifa][1], 1);
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    RansTurbKOmSST::boundaryHookScalarRansTurb(phi_fa, zone, name);

    //  Note that this function is called both for boundaries defined as HOOK and for all walls
    //  Hence we must do the stuff below only if it is NOT a wall
    Param *param;
    if (getParam(param, zone->getName()))
      if (param->getString() == "HOOK") {
        if (name == "kine")
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            phi_fa[ifa] = getValues(x_fa[ifa][1], 3);
        if (name == "omega")
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            phi_fa[ifa] = getValues(x_fa[ifa][1], 4);
      }
  }
 
  virtual void initialHook()
  {
    JoeKOM_mixing::initialHook() ;
    //  Verify that only a single boundary has been defined as HOOK
    //  Otherwise will get trouble in boundaryHook below...
    if (mpi_rank==0)  cout << "hyshot_mixing::initialHook: counting number of HOOK boundaries..." << endl ;
    int nhooks=0;
    Param *param;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (getParam(param, zone->getName()))
          if (param->getString() == "HOOK")  nhooks++ ;
    if ( nhooks != 1 ) {
      if (mpi_rank==0)  cout << "ERROR : " << nhooks << " boundaries defined as HOOK..." << endl ;
      exit(-1) ;
    }

    //  Set initial field from inflow profile
    if (!checkScalarFlag("RHO")) {   
      if (mpi_rank==0)  cout << "Setting initial N-S fields from inflow profile..." << endl ;
      for (int icv = 0; icv < ncv; icv++) {
        double Pprof = getValues(x_cv[icv][1], 1);
        double Tprof = getValues(x_cv[icv][1], 2);
        double Uprof = getValues(x_cv[icv][1], 5);
        double Vprof = getValues(x_cv[icv][1], 6);
        double Wprof = 0.0;
        double Hprof, dum1, dum2 ;
        ComputeProperties_T(Hprof, RoM[icv], gamma[icv], dum1, dum2, Tprof, ZMean[icv], 0.0, 0.0);
        rho[icv] = Pprof / (RoM[icv] * Tprof) ;
        rhou[icv][0] = rho[icv] * Uprof ;
        rhou[icv][1] = rho[icv] * Vprof ;
        rhou[icv][2] = rho[icv] * Wprof ;
        rhoE[icv] = Hprof*rho[icv] - Pprof + 0.5*vecDotVec3d(rhou[icv],rhou[icv])/rho[icv] ;
      }
      updateCvData(rhou, REPLACE_ROTATE_DATA);
      updateCvData(rho,  REPLACE_DATA);
      updateCvData(rhoE, REPLACE_DATA);
    }
  }

};






// ###########################################################################################
// ###########################################################################################

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
  int run = 0 ;

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
      case 1:  joe = new JoeKOM_mixing(inputFileName);  break;
      case 2:  joe = new JoeKOM_fpva_coeff(inputFileName);  break;
      case 3:  joe = new hyshot_fpva_coeff(inputFileName);  break;
      case 4:  joe = new hyshot_mixing(inputFileName);  break;
      default:  if (mpi_rank==0)  cout << "ERROR : MUST SPECIFY RUN NUMBER..." << endl ;
    }
    
    // provide total runtime 
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    // run joe
    if ( joe != NULL )  joe->run();
    

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

    // delete joe (make sure memory is deallocated in destructors
    if ( joe != NULL )  delete joe;
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



