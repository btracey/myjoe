#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

#include "combModels/CombModel_BinaryMixing.h"
#include "combModels/CombModel_Mixing.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Advection of Binary Mixture                                ------
// ------                                                                               ------
// ###########################################################################################
class MyJoeBinaryMixing : public JoeWithModels, public RansCombBinaryMixing
{
public:

  MyJoeBinaryMixing(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoeBinaryMixing()" << endl;
  }

  virtual ~MyJoeBinaryMixing()  {}

  void initialHook() 
  {
    // First compute and update mixture profile ZMean 
    string init_profile = getStringParam("INIT_PROFILE", "SMOOTH_HAT");
    
    // Sinusoidal function
    if (init_profile == "SINE")
    {
      for (int icv = 0; icv < ncv; icv++)
        ZMean[icv] = 0.5 * (1.0 - cos(2.0 * M_PI * x_cv[icv][0]));
    }
    
    // Smoothed top-hat function
    else if (init_profile == "SMOOTH_HAT")
    {
      const double delta = 1.5*10.0/256.0; // smooth initial condition 10 points in 256
      for (int icv = 0; icv < ncv; icv++)
      {
        double rx = fabs(x_cv[icv][0]-0.5);
        if(rx <= 0.25-delta)
          ZMean[icv] = 1.0;
        else if (rx >= 0.25+delta)
          ZMean[icv] = 0.0;
        else
        {
          double alpha_He = (1.0 + cos(M_PI * (rx - (0.25 - delta)) / (2.0*delta))) / 2.0;
          ZMean[icv] = alpha_He * 1.0 + (1.0 - alpha_He) * 0.0;
        }
      }
    }
    else
    {
      cerr << "### Init profile not recognized, should be SINE or SMOOTH_HAT! ###" << endl;
      throw(-1);
    }
    
    updateCvDataByName("ZMean", REPLACE_DATA);

    // Read initial conditions and case
    double Uin   = getDoubleParam("U_INITIAL", "1000");
    double Pin   = getDoubleParam("P_REF", "100000.0");
    double Tin   = getDoubleParam("T_REF", "1300.0");
    double RHOin = getDoubleParam("RHO_REF", "0.26738");
    
    double Tfuel = getDoubleParam("T_FUEL", "300.0");
    double Toxy  = getDoubleParam("T_OXY", "1300.0");
    double hfuel, hoxy, RoMfuel, RoMoxy, gamfuel, gamoxy, dum1, dum2;

    rho_1 = Pin / (RoM_1 * Tfuel);
    rho_2 = Pin / (RoM_2 * Toxy);
    double rhofuel, rhooxy;
    calcThermoProp_Tp(hfuel, rhofuel, RoMfuel, gamfuel, Tfuel, Pin, 1.0);
    calcThermoProp_Tp(hoxy,  rhooxy,  RoMoxy,  gamoxy,  Toxy,  Pin, 0.0);
    
    string caseoption = getStringParam("CASE", "T_CST");

    // Compute and update all conserved variables
    for (int icv = 0; icv < ncv; icv++)
    {
      double pressure = Pin;

      if (caseoption == "T_CST")
      {
        calcThermoProp_Tp(enthalpy[icv], rho[icv], RoM[icv], gamma[icv], Tin, Pin, ZMean[icv]);
      }
      else if (caseoption == "RHO_CST")
      {
        rho[icv] = RHOin;
        calcThermoProp_prho(temp[icv], enthalpy[icv], RoM[icv], gamma[icv], pressure, rho[icv], ZMean[icv]);
      } 
      else if (caseoption == "T_FIXED")
      {
        enthalpy[icv] = ZMean[icv] * hfuel + (1.0 - ZMean[icv]) * hoxy; 
        calcThermoProp_Hp(temp[icv], rho[icv], RoM[icv], gamma[icv], enthalpy[icv], pressure, ZMean[icv]);
      }
      else
      {
        cerr << "### Case option not recognized! ###" << endl;
        throw(-1);
      }

      rhou[icv][0] = rho[icv] * Uin;
      rhou[icv][1] = 0.0;
      rhou[icv][2] = 0.0;
      rhoE[icv] = rho[icv] * enthalpy[icv] - pressure + 0.5 * rho[icv] * Uin * Uin;
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);

    if (mpi_rank == 0) cout << "Field initialized" << endl;
  }  
  
  void finalHook()
  {
    // Compute L1, L2 and Linf error for T, P, rho, U, and Z

    double myErr_T_L1   = 0.0, myErr_T_L2   = 0.0, myErr_T_Linf   = 0.0;
    double myErr_P_L1   = 0.0, myErr_P_L2   = 0.0, myErr_P_Linf   = 0.0;
    double myErr_rho_L1 = 0.0, myErr_rho_L2 = 0.0, myErr_rho_Linf = 0.0;
    double myErr_U_L1   = 0.0, myErr_U_L2   = 0.0, myErr_U_Linf   = 0.0;
    double myErr_Z_L1   = 0.0, myErr_Z_L2   = 0.0, myErr_Z_Linf   = 0.0;
    
    double Err_T_L1   = 0.0, Err_T_L2   = 0.0, Err_T_Linf   = 0.0;
    double Err_P_L1   = 0.0, Err_P_L2   = 0.0, Err_P_Linf   = 0.0;
    double Err_rho_L1 = 0.0, Err_rho_L2 = 0.0, Err_rho_Linf = 0.0;
    double Err_U_L1   = 0.0, Err_U_L2   = 0.0, Err_U_Linf   = 0.0;
    double Err_Z_L1   = 0.0, Err_Z_L2   = 0.0, Err_Z_Linf   = 0.0;
   

    // Read initial conditions and case
    double Uin   = getDoubleParam("U_INITIAL", "1000");
    double Pin   = getDoubleParam("P_REF", "100000.0");
    double Tin   = getDoubleParam("T_REF", "1300.0");
    double RHOin = getDoubleParam("RHO_REF", "0.26738");
    
    double Tfuel = getDoubleParam("T_FUEL", "300.0");
    double Toxy  = getDoubleParam("T_OXY", "1300.0");
    double hfuel, hoxy, RoMfuel, RoMoxy, gamfuel, gamoxy, dum1, dum2;
    ComputeProperties_T(hfuel, RoMfuel, gamfuel, dum1, dum2, Tfuel, 1.0, 0.0, 0.0);
    ComputeProperties_T(hoxy,  RoMoxy,  gamoxy,  dum1, dum2, Toxy,  0.0, 0.0, 0.0);
  
    string caseoption = getStringParam("CASE", "T_CST");
    
    double Pstar = Pin;
    double Ustar = Uin;
    double Zstar = 1.0;
    
    double Tstar = 0.0, myTstar = 0.0;
    double rhostar = 0.0, myrhostar = 0.0;
    
    string init_profile = getStringParam("INIT_PROFILE", "SMOOTH_HAT");
    const double delta = 1.5*10.0/256.0; // smooth initial condition 10 points in 256
    for (int icv = 0; icv < ncv; icv++)
    {
      double T, P=Pin, RHO, U=Uin, Z;
      double Hin, RoMin, gamin, dum1, dum2;
      
      // Sinusoidal function
      if (init_profile == "SINE")
      {
        Z = 0.5 * (1.0 - cos(2.0 * M_PI * x_cv[icv][0]));
      }
      else if (init_profile == "SMOOTH_HAT")
      {      
        double rx = fabs(x_cv[icv][0]-0.5);
        if(rx <= 0.25-delta)
          Z = 1.0;
        else if (rx >= 0.25+delta)
          Z = 0.0;
        else
        {
          double alpha_He = (1.0 + cos(M_PI * (rx - (0.25 - delta)) / (2.0*delta))) / 2.0;
          Z = alpha_He * 1.0 + (1.0 - alpha_He) * 0.0;
        }
      }

      if (caseoption == "T_CST")
      {
        T = Tin;
        ComputeProperties_T(Hin, RoMin, gamin, dum1, dum2, Tin, Z, 0.0, 0.0);
        RHO = P / (RoMin * Tin);  
        myTstar = Tin;
        myrhostar = max(RHO, myrhostar);
      }
      else if (caseoption == "RHO_CST")
      {
        RHO = RHOin;
        // just to get RoM first
        ComputeProperties_T(Hin, RoMin, gamin, dum1, dum2, Tin, Z, 0.0, 0.0);
        T = P / (RHO * RoMin);
        myTstar = max(T, myTstar);
        myrhostar = RHO;
      }
      else if (caseoption == "T_FIXED")
      {
        Hin = Z * hfuel + (1.0 - Z) * hoxy; 
        ComputeProperties_H(T, RoMin, gamin, dum1, dum2, Hin, Z, 0.0, 0.0);
        RHO = P / (RoMin * T); 
        myTstar = max(T, myTstar);
        myrhostar = max(RHO, myrhostar);
      }
      
      // compute errors
      double Err_T   = temp[icv] - T;
      double Err_P   = press[icv] - P;
      double Err_rho = rho[icv] - RHO;
      double Err_U   = vel[icv][0] - U;
      double Err_Z   = ZMean[icv] - Z;
      
      myErr_T_L1   += fabs(Err_T);        myErr_T_L2   += (Err_T   * Err_T);        myErr_T_Linf   = max(myErr_T_Linf,   fabs(Err_T));
      myErr_P_L1   += fabs(Err_P);        myErr_P_L2   += (Err_P   * Err_P);        myErr_P_Linf   = max(myErr_P_Linf,   fabs(Err_P));
      myErr_rho_L1 += fabs(Err_rho);      myErr_rho_L2 += (Err_rho * Err_rho);      myErr_rho_Linf = max(myErr_rho_Linf, fabs(Err_rho));
      myErr_U_L1   += fabs(Err_U);        myErr_U_L2   += (Err_U   * Err_U);        myErr_U_Linf   = max(myErr_U_Linf,   fabs(Err_U));
      myErr_Z_L1   += fabs(Err_Z);        myErr_Z_L2   += (Err_Z   * Err_Z);        myErr_Z_Linf   = max(myErr_Z_Linf,   fabs(Err_Z));
    }
    
    // add over all processors and reduce
    MPI_Allreduce(&myErr_T_L1,   &Err_T_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_T_L2,   &Err_T_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_T_Linf, &Err_T_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_P_L1,   &Err_P_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_P_L2,   &Err_P_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm); 
    MPI_Allreduce(&myErr_P_Linf, &Err_P_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_rho_L1,   &Err_rho_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_rho_L2,   &Err_rho_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_rho_Linf, &Err_rho_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_U_L1,   &Err_U_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_U_L2,   &Err_U_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_U_Linf, &Err_U_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_Z_L1,   &Err_Z_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_Z_L2,   &Err_Z_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_Z_Linf, &Err_Z_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    int ncv_tot;
    MPI_Allreduce(&ncv, &ncv_tot, 1, MPI_INT, MPI_SUM, mpi_comm);
    
    MPI_Allreduce(&myTstar,   &Tstar,   1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    MPI_Allreduce(&myrhostar, &rhostar, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    // rescale
    Err_T_L1     = Err_T_L1     / (double(ncv_tot) * Tstar);
    Err_T_L2     = Err_T_L2     / (double(ncv_tot));                  Err_T_L2   = sqrt(Err_T_L2)   / Tstar;
    Err_T_Linf   = Err_T_Linf   / (Tstar);
    Err_P_L1     = Err_P_L1     / (double(ncv_tot) * Pstar);
    Err_P_L2     = Err_P_L2     / (double(ncv_tot));                  Err_P_L2   = sqrt(Err_P_L2)   / Pstar;
    Err_P_Linf   = Err_P_Linf   / (Pstar);
    Err_rho_L1   = Err_rho_L1   / (double(ncv_tot) * rhostar);
    Err_rho_L2   = Err_rho_L2   / (double(ncv_tot));                  Err_rho_L2 = sqrt(Err_rho_L2) / rhostar;
    Err_rho_Linf = Err_rho_Linf / (rhostar);
    Err_U_L1     = Err_U_L1     / (double(ncv_tot) * Ustar);
    Err_U_L2     = Err_U_L2     / (double(ncv_tot));                  Err_U_L2   = sqrt(Err_U_L2)   / Ustar;
    Err_U_Linf   = Err_U_Linf   / (Ustar);
    Err_Z_L1     = Err_Z_L1     / (double(ncv_tot) * Zstar);
    Err_Z_L2     = Err_Z_L2     / (double(ncv_tot));                  Err_Z_L2   = sqrt(Err_Z_L2)   / Zstar;
    Err_Z_Linf   = Err_Z_Linf   / (Zstar);
    
    // output to file
    if (mpi_rank == 0)
    {
      FILE *fp;
      fp = fopen("error.dat", "at");
      fprintf(fp, "%.4i\t%.8e\t", ncv_tot, 1.0/int(ncv_tot));
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t", Err_T_L1,   Err_P_L1,   Err_rho_L1,   Err_U_L1,   Err_Z_L1);
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t", Err_T_L2,   Err_P_L2,   Err_rho_L2,   Err_U_L2,   Err_Z_L2);
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Err_T_Linf, Err_P_Linf, Err_rho_Linf, Err_U_Linf, Err_Z_Linf);
      fclose(fp);
    }
  }
};


// ###########################################################################################
// ------                                                                               ------
// ------                    Advection of Real 2 Components Mixture                     ------
// ------                                                                               ------
// ###########################################################################################
class MyJoeMixing : public JoeWithModels, public RansCombMixing
{
public:

  MyJoeMixing(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "RansCombMixing()" << endl;
  }

  virtual ~MyJoeMixing()  {}
  
  void initialHook() 
  {
    // First compute and update mixture profile ZMean 
    string init_profile = getStringParam("INIT_PROFILE", "SMOOTH_HAT");
    
    // Sinusoidal function
    if (init_profile == "SINE")
    {
      for (int icv = 0; icv < ncv; icv++)
        ZMean[icv] = 0.5 * (1.0 - cos(2.0 * M_PI * x_cv[icv][0]));
    }
    
    // Smoothed top-hat function
    else if (init_profile == "SMOOTH_HAT")
    {
      const double delta = 1.5*10.0/256.0; // smooth initial condition 10 points in 256
      for (int icv = 0; icv < ncv; icv++)
      {
        double rx = fabs(x_cv[icv][0]-0.5);
        if(rx <= 0.25-delta)
          ZMean[icv] = 1.0;
        else if (rx >= 0.25+delta)
          ZMean[icv] = 0.0;
        else
        {
          double alpha_He = (1.0 + cos(M_PI * (rx - (0.25 - delta)) / (2.0*delta))) / 2.0;
          ZMean[icv] = alpha_He * 1.0 + (1.0 - alpha_He) * 0.0;
        }
      }
    }
    else
    {
      cerr << "### Init profile not recognized, should be SINE or SMOOTH_HAT! ###" << endl;
      throw(-1);
    }
    
    updateCvDataByName("ZMean", REPLACE_DATA);

    // Read initial conditions and case
    double Uin   = getDoubleParam("U_INITIAL", "1000");
    double Pin   = getDoubleParam("P_REF", "100000.0");
    double Tin   = getDoubleParam("T_REF", "1300.0");
    double RHOin = getDoubleParam("RHO_REF", "0.26738");
    
    double Tfuel = getDoubleParam("T_FUEL", "300.0");
    double Toxy  = getDoubleParam("T_OXY", "1300.0");
    double hfuel, hoxy, RoMfuel, RoMoxy, gamfuel, gamoxy, dum1, dum2;
    ComputeProperties_T(hfuel, RoMfuel, gamfuel, dum1, dum2, Tfuel, 1.0, 0.0, 0.0);
    ComputeProperties_T(hoxy,  RoMoxy,  gamoxy,  dum1, dum2, Toxy,  0.0, 0.0, 0.0);
    
    string caseoption = getStringParam("CASE", "T_CST");

    // Compute and update all conserved variables
    for (int icv = 0; icv < ncv; icv++)
    {
      double pressure = Pin;

      if (caseoption == "T_CST")
      {
        ComputeProperties_T(enthalpy[icv], RoM[icv], gamma[icv], dum1, dum2, Tin, ZMean[icv], 0.0, 0.0);
        rho[icv] = pressure / (RoM[icv] * Tin);  
      }
      else if (caseoption == "RHO_CST")
      {
        rho[icv] = RHOin;
        // just to get RoM first
        ComputeProperties_T(enthalpy[icv], RoM[icv], gamma[icv], dum1, dum2, Tin, ZMean[icv], 0.0, 0.0);
        double temper = pressure / (rho[icv] * RoM[icv]);
        ComputeProperties_T(enthalpy[icv], RoM[icv], gamma[icv], dum1, dum2, temper, ZMean[icv], 0.0, 0.0);
      }
      else if (caseoption == "T_FIXED")
      {
        enthalpy[icv] = ZMean[icv] * hfuel + (1.0 - ZMean[icv]) * hoxy; 
        ComputeProperties_H(temp[icv], RoM[icv], gamma[icv], dum1, dum2, enthalpy[icv], ZMean[icv], 0.0, 0.0);
        rho[icv] = pressure / (RoM[icv] * temp[icv]);
      }
      else
      {
        cerr << "### Case option not recognized! ###" << endl;
        throw(-1);
      }

      rhou[icv][0] = rho[icv] * Uin;
      rhou[icv][1] = 0.0;
      rhou[icv][2] = 0.0;
      rhoE[icv] = rho[icv] * enthalpy[icv] - pressure + 0.5 * rho[icv] * Uin * Uin;
    }
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);

    if (mpi_rank == 0) cout << "Field initialized" << endl;
  }  
  
  void finalHook()
  {
    // Compute L1, L2 and Linf error for T, P, rho, U, and Z

    double myErr_T_L1   = 0.0, myErr_T_L2   = 0.0, myErr_T_Linf   = 0.0;
    double myErr_P_L1   = 0.0, myErr_P_L2   = 0.0, myErr_P_Linf   = 0.0;
    double myErr_rho_L1 = 0.0, myErr_rho_L2 = 0.0, myErr_rho_Linf = 0.0;
    double myErr_U_L1   = 0.0, myErr_U_L2   = 0.0, myErr_U_Linf   = 0.0;
    double myErr_Z_L1   = 0.0, myErr_Z_L2   = 0.0, myErr_Z_Linf   = 0.0;
    
    double Err_T_L1   = 0.0, Err_T_L2   = 0.0, Err_T_Linf   = 0.0;
    double Err_P_L1   = 0.0, Err_P_L2   = 0.0, Err_P_Linf   = 0.0;
    double Err_rho_L1 = 0.0, Err_rho_L2 = 0.0, Err_rho_Linf = 0.0;
    double Err_U_L1   = 0.0, Err_U_L2   = 0.0, Err_U_Linf   = 0.0;
    double Err_Z_L1   = 0.0, Err_Z_L2   = 0.0, Err_Z_Linf   = 0.0;
   

    // Read initial conditions and case
    double Uin   = getDoubleParam("U_INITIAL", "1000");
    double Pin   = getDoubleParam("P_REF", "100000.0");
    double Tin   = getDoubleParam("T_REF", "1300.0");
    double RHOin = getDoubleParam("RHO_REF", "0.26738");
    
    double Tfuel = getDoubleParam("T_FUEL", "300.0");
    double Toxy  = getDoubleParam("T_OXY", "1300.0");
    double hfuel, hoxy, RoMfuel, RoMoxy, gamfuel, gamoxy, dum1, dum2;
    ComputeProperties_T(hfuel, RoMfuel, gamfuel, dum1, dum2, Tfuel, 1.0, 0.0, 0.0);
    ComputeProperties_T(hoxy,  RoMoxy,  gamoxy,  dum1, dum2, Toxy,  0.0, 0.0, 0.0);
  
    string caseoption = getStringParam("CASE", "T_CST");
    
    double Pstar = Pin;
    double Ustar = Uin;
    double Zstar = 1.0;
    
    double Tstar = 0.0, myTstar = 0.0;
    double rhostar = 0.0, myrhostar = 0.0;
    
    string init_profile = getStringParam("INIT_PROFILE", "SMOOTH_HAT");
    const double delta = 1.5*10.0/256.0; // smooth initial condition 10 points in 256
    for (int icv = 0; icv < ncv; icv++)
    {
      double T, P=Pin, RHO, U=Uin, Z;
      double Hin, RoMin, gamin, dum1, dum2;
      
      // Sinusoidal function
      if (init_profile == "SINE")
      {
        Z = 0.5 * (1.0 - cos(2.0 * M_PI * x_cv[icv][0]));
      }
      else if (init_profile == "SMOOTH_HAT")
      {      
        double rx = fabs(x_cv[icv][0]-0.5);
        if(rx <= 0.25-delta)
          Z = 1.0;
        else if (rx >= 0.25+delta)
          Z = 0.0;
        else
        {
          double alpha_He = (1.0 + cos(M_PI * (rx - (0.25 - delta)) / (2.0*delta))) / 2.0;
          Z = alpha_He * 1.0 + (1.0 - alpha_He) * 0.0;
        }
      }

      if (caseoption == "T_CST")
      {
        T = Tin;
        ComputeProperties_T(Hin, RoMin, gamin, dum1, dum2, Tin, Z, 0.0, 0.0);
        RHO = P / (RoMin * Tin);  
        myTstar = Tin;
        myrhostar = max(RHO, myrhostar);
      }
      else if (caseoption == "RHO_CST")
      {
        RHO = RHOin;
        // just to get RoM first
        ComputeProperties_T(Hin, RoMin, gamin, dum1, dum2, Tin, Z, 0.0, 0.0);
        T = P / (RHO * RoMin);
        myTstar = max(T, myTstar);
        myrhostar = RHO;
      }
      else if (caseoption == "T_FIXED")
      {
        Hin = Z * hfuel + (1.0 - Z) * hoxy; 
        ComputeProperties_H(T, RoMin, gamin, dum1, dum2, Hin, Z, 0.0, 0.0);
        RHO = P / (RoMin * T); 
        myTstar = max(T, myTstar);
        myrhostar = max(RHO, myrhostar);
      }
      
      // compute errors
      double Err_T   = temp[icv] - T;
      double Err_P   = press[icv] - P;
      double Err_rho = rho[icv] - RHO;
      double Err_U   = vel[icv][0] - U;
      double Err_Z   = ZMean[icv] - Z;
      
      myErr_T_L1   += fabs(Err_T);        myErr_T_L2   += Err_T   * Err_T;        myErr_T_Linf   = max(myErr_T_Linf,   fabs(Err_T));
      myErr_P_L1   += fabs(Err_P);        myErr_P_L2   += Err_P   * Err_P;        myErr_P_Linf   = max(myErr_P_Linf,   fabs(Err_P));
      myErr_rho_L1 += fabs(Err_rho);      myErr_rho_L2 += Err_rho * Err_rho;      myErr_rho_Linf = max(myErr_rho_Linf, fabs(Err_rho));
      myErr_U_L1   += fabs(Err_U);        myErr_U_L2   += Err_U   * Err_U;        myErr_U_Linf   = max(myErr_U_Linf,   fabs(Err_U));
      myErr_Z_L1   += fabs(Err_Z);        myErr_Z_L2   += Err_Z   * Err_Z;        myErr_Z_Linf   = max(myErr_Z_Linf,   fabs(Err_Z));
    }
    
    // add over all processors and reduce
    MPI_Allreduce(&myErr_T_L1,   &Err_T_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_T_L2,   &Err_T_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_T_Linf, &Err_T_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_P_L1,   &Err_P_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_P_L2,   &Err_P_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm); 
    MPI_Allreduce(&myErr_P_Linf, &Err_P_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_rho_L1,   &Err_rho_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_rho_L2,   &Err_rho_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_rho_Linf, &Err_rho_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_U_L1,   &Err_U_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_U_L2,   &Err_U_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_U_Linf, &Err_U_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    MPI_Allreduce(&myErr_Z_L1,   &Err_Z_L1,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_Z_L2,   &Err_Z_L2,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myErr_Z_Linf, &Err_Z_Linf, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    int ncv_tot;
    MPI_Allreduce(&ncv, &ncv_tot, 1, MPI_INT, MPI_SUM, mpi_comm);
    
    MPI_Allreduce(&myTstar,   &Tstar,   1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    MPI_Allreduce(&myrhostar, &rhostar, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    
    // rescale
    Err_T_L1     = Err_T_L1     / (double(ncv_tot) * Tstar);
    Err_T_L2     = Err_T_L2     / (double(ncv_tot));                  Err_T_L2   = sqrt(Err_T_L2)   / Tstar;
    Err_T_Linf   = Err_T_Linf   / (Tstar);
    Err_P_L1     = Err_P_L1     / (double(ncv_tot) * Pstar);
    Err_P_L2     = Err_P_L2     / (double(ncv_tot));                  Err_P_L2   = sqrt(Err_P_L2)   / Pstar;
    Err_P_Linf   = Err_P_Linf   / (Pstar);
    Err_rho_L1   = Err_rho_L1   / (double(ncv_tot) * rhostar);
    Err_rho_L2   = Err_rho_L2   / (double(ncv_tot));                  Err_rho_L2 = sqrt(Err_rho_L2) / rhostar;
    Err_rho_Linf = Err_rho_Linf / (rhostar);
    Err_U_L1     = Err_U_L1     / (double(ncv_tot) * Ustar);
    Err_U_L2     = Err_U_L2     / (double(ncv_tot));                  Err_U_L2   = sqrt(Err_U_L2)   / Ustar;
    Err_U_Linf   = Err_U_Linf   / (Ustar);
    Err_Z_L1     = Err_Z_L1     / (double(ncv_tot) * Zstar);
    Err_Z_L2     = Err_Z_L2     / (double(ncv_tot));                  Err_Z_L2   = sqrt(Err_Z_L2)   / Zstar;
    Err_Z_Linf   = Err_Z_Linf   / (Zstar);
    
    // output to file
    if (mpi_rank == 0)
    {
      FILE *fp;
      fp = fopen("error.dat", "at");
      fprintf(fp, "%.4i\t%.8e\t", ncv_tot, 1.0/int(ncv_tot));
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t", Err_T_L1,   Err_P_L1,   Err_rho_L1,   Err_U_L1,   Err_Z_L1);
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t", Err_T_L2,   Err_P_L2,   Err_rho_L2,   Err_U_L2,   Err_Z_L2);
      fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Err_T_Linf, Err_P_Linf, Err_rho_Linf, Err_U_Linf, Err_Z_Linf);
      fclose(fp);
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
    case 0:   joe = new MyJoeBinaryMixing(inputFileName);        break;
    case 1:   joe = new MyJoeMixing(inputFileName);              break;
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



