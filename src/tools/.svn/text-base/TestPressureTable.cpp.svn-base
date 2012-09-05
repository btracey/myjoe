/*
 * TestPressureTable.cpp
 *
 *  Created on: Feb 18, 2011
 *      Author: ronanv
 */

#include "TestPressureTable.h"

int main(int argc, char * argv[]) {

  // Read Table filename
  if (argc > 0) {
    ChemTableFilename = argv[1];
  }


  time_t start0;
  time(&start0);
  PressureTable NewTable;
  NewTable.Load(ChemTableFilename);
  NewTable.print();

  int nPress = NewTable.GetDimension1();
  int nZm = NewTable.GetDimension2();
  int nZv = NewTable.GetDimension3();
  int nC  = NewTable.GetDimension4();
  // Check Pressure array
  //double *myP = new double[nPress];
  //NewTable.CopyCoordinate1(myP);
  //for(int i=0; i<nPress; ++i) cout << myP[i] << " ";
  //cout << endl;

  //ChemtableCartesianLinear_single ChemTableFPVA;
  //ChemTableFPVA.Load("RUN_P2.9000e05/database");

  double Zm = 0.2;
  double Zvar = 0.1;
  double C = 0.233;
  cout << "Interpolating T0" << endl;
  //cout << NewTable.Lookup(2.9e5, Zm, Zvar, C, "T0")<< endl;
  // Get Interpolation Index used previously
  InterpolationIndex toto(4,2);
  //NewTable.GetInterpolation(toto);
  //toto.print();
  // Copy toto to tata and add one to Prog dimension
  InterpolationIndex tata(4,2);
  //toto.copy(tata);
  //tata.index[3] += 1;
  //tata.print();

  //NewTable.SetInterpolation(tata);
  //double temp = NewTable.Lookup("T0");
  //cout << "ic + 1 : " << temp <<endl;

  /*
  // --------------------------
  // Test change of interpolation mode
  // --------------------------
  time_t start, end;
  time(&start);
  for (int i=0; i< 1e6 ; ++i) {
    C = ( double (rand()) )/( double (RAND_MAX) );
    //cout << "C = " << C << endl;
    Zm = ( double (rand()) )/( double (RAND_MAX) );
    //cout << "Z = " << Zm << endl;
    Zvar = ( double (rand()) )/( double (RAND_MAX) );
    //cout << "Zvar = " << Zvar << endl;

    NewTable.InterpolatePoint_old(2.9e5,Zm, Zvar,C);
    NewTable.GetInterpolation(toto);
    //toto.print();
    NewTable.InterpolatePoint(2.9e5,Zm, Zvar,C);
    NewTable.GetInterpolation(tata);
    //tata.print();
    if (!toto.compare(tata)) {
      cout << "OLD and NEW disagree" << endl;
      throw(-1);
    }
  }
  time(&end);
  cout << difftime(end, start) << " seconds have passed in test 0" << endl;
  // Time for old interp only
  time(&start);
  for (int i=0; i< 4e7 ; ++i) {
    C = ( double (rand()) )/( double (RAND_MAX) );
    Zm = ( double (rand()) )/( double (RAND_MAX) );
    Zvar = ( double (rand()) )/( double (RAND_MAX) );
    NewTable.InterpolatePoint_old(2.9e5,Zm, Zvar,C);
  }
  time(&end);
  cout << difftime(end, start) << " seconds have passed in test 1" << endl;
  // Time for new interp only
  time(&start);
  for (int i=0; i< 4e7 ; ++i) {
    C = ( double (rand()) )/( double (RAND_MAX) );
    Zm = ( double (rand()) )/( double (RAND_MAX) );
    Zvar = ( double (rand()) )/( double (RAND_MAX) );
    NewTable.InterpolatePoint(2.9e5,Zm, Zvar,C);
  }
  time(&end);
  cout << difftime(end, start) << " seconds have passed in test 2" << endl;
  cout << difftime(end, start0) << " seconds have passed in total" << endl;
*/

/*
  // --------------------------
  // Test Newton/Broyden ierations to compute C,P
  // --------------------------

  // Here is the solution I want to retrieve
  C = 0.0;
  double P = 0.879425145e5;
  Zm = 0.000621758;
  double Sz = 0.068987; // normalized
  // What are the corresponding rho, rhoZ, rhoZv, RhoYc
  //double rho = NewTable.Lookup(P, Zm, Sz, C,"rho0");
  //double rho = NewTable.Lookup(P, Zm, Sz, 1.0,"rho0");
  double rho = 0.612164;
  //double Zv = Sz*Zm*(1.0-Zm); // non-normalized
  double Zv =  0.00220893;
  //double Yc = NewTable.Lookup(P, Zm, Sz, C,"PROG");
  //double Yc_eq  = NewTable.Lookup(P,Zm,Sz,1.0,"PROG");
  //double  Yc = C*Yc_eq;
  double Yc = 0.00501301;
  cout << "Inputs :" << endl;
  cout << "  - rho   : " << rho << endl;
  cout << "  - Z  : " << Zm << endl;
  cout << "  - Sz : " << Sz << endl;
  cout << "  - Yc : " << Yc << endl;
  cout << endl;
  // Initial guess:
  double C_ini = 0.5;
  double P_ini = 1.0e5;
  // Call Newton/Broyden method
  double Sz_out, C_out, P_out;
  NewTable.ComputeDatabaseInputs(rho, Zm, Zv, Yc,
                                 P_out, Sz_out, C_out,
                                 P_ini, C_ini);
  // Ouput:
  cout << "Solution :" << endl;
  cout << "  - P  : " << P_out  << " vs " << P    << endl;
  cout << "  - Z  : " << Zm     << endl;
  cout << "  - Sz : " << Sz_out << " vs " << Sz << endl;
  cout << "  - C  : " << C_out  << " vs " << C    << endl;
  cout << endl;

  double rho_test = NewTable.Lookup(P_out, Zm, Sz_out, C_out,"rho0");
  double Yc_test  = C_out*NewTable.Lookup(P,Zm,Sz,1.0,"PROG");
  cout << "Outputs :" << endl;
  cout << "  - rho : " << rho_test << endl;
  cout << "  - Yc  : " << Yc_test << endl;
  cout << endl;

  double P_scale = PressureTable::GetPscale();
  cout << "L2 error : " << sqrt( (P-P_out)/P_scale*(P-P_out)/P_scale
                                 + (Sz-Sz_out)*(Sz-Sz_out)
                                 + (C-C_out)*(C-C_out) ) << endl;

  */
  
  // --------------------------
  // Compare Pressure scaling
  // --------------------------
  ofstream tecplotfile;
  tecplotfile.open("Compare_pressure_fixed_Z_SZ.dat");
  tecplotfile << "VARIABLES = Pressure HR HR_approx P_scaled HR_scaled HR_approx_scaled"
              << " SRC SRC_approx SRC_scaled SRC_appprox_scaled" << endl;
  Zm = 0.03;
  Zvar = 0.0;
  //C = 0.5;
  double press0 = 1.5e5;
  int nP = 20;
  double HR0, HR, HR_approx;
  double SRC0, SRC, SRC_approx;

  int nzone = 10;
  for (int izone = 0; izone < nzone ; ++izone) {
    C = ((double) izone)/((double) nzone - 1.0);
    tecplotfile << "ZONE T=\"C=" << C << "\", I=" << nP << endl;
    HR0 = NewTable.Lookup(press0, Zm, Zvar, C,"HeatRelease");
    SRC0 = NewTable.Lookup(press0, Zm, Zvar, C,"SRC_PROG");
    double rho0 = NewTable.Lookup(press0, Zm, Zvar, C,"rho0");
    for (int i=0; i<nP; ++i) {
      double press = 5.0e4 + (5.0e5-5.0e4)*((double) i)/((double) nP - 1.0);
      HR = NewTable.Lookup(press, Zm, Zvar, C,"HeatRelease");
      double rho = NewTable.Lookup(press, Zm, Zvar, C,"rho0");
      HR_approx = pow(press/press0,2.0)*HR0;
      SRC = NewTable.Lookup(press, Zm, Zvar, C,"SRC_PROG");
      SRC_approx = pow(press/press0,2.0)*SRC0;
      tecplotfile << press << " " << rho*HR << " " << rho*HR_approx << " "
                  << press/press0 << " " << rho*HR/(rho0*HR0) << " " << pow(press/press0,2.0) << " "
                  << SRC << " " << rho*SRC_approx << " "
                  << rho*SRC/(rho0*SRC0) << " " << pow(press/press0,2.0)
                  << endl;
    }
  }
  tecplotfile.close();

  //
  tecplotfile.open("Compare_pressure_fixed_C_SZ.dat");
  tecplotfile << "VARIABLES = Pressure HR HR_approx P_scaled HR_scaled HR_approx_scaled"
              << " SRC SRC_approx SRC_scaled SRC_appprox_scaled" << endl;
  Zvar = 0.0;
  C = 0.5;
  nzone = 30;
  for (int izone = 1; izone < nzone-1 ; ++izone) {
    Zm = ((double) izone)/((double) nzone - 1.0);
    tecplotfile << "ZONE T=\"Zm=" << Zm << "\", I=" << nP << endl;
    HR0 = NewTable.Lookup(press0,  Zm, Zvar, C,"HeatRelease");
    SRC0 = NewTable.Lookup(press0, Zm, Zvar, C,"SRC_PROG");
    double rho0 = NewTable.Lookup(press0, Zm, Zvar, C,"rho0");
    for (int i=0; i<nP; ++i) {
      double press = 5.0e4 + (5.0e5-5.0e4)*((double) i)/((double) nP - 1.0);
      HR = NewTable.Lookup(press, Zm, Zvar, C,"HeatRelease");
      double rho = NewTable.Lookup(press, Zm, Zvar, C,"rho0");
      HR_approx = pow(press/press0,2.0)*HR0;
      SRC = NewTable.Lookup(press, Zm, Zvar, C,"SRC_PROG");
      SRC_approx = pow(press/press0,2.0)*SRC0;
      tecplotfile << press << " " << rho*HR << " " << rho*HR_approx << " "
                  << press/press0 << " " << rho*HR/(rho0*HR0) << " " << pow(press/press0,2.0) << " "
                  << SRC << " " << rho*SRC_approx << " "
                  << rho*SRC/(rho0*SRC0) << " " << pow(press/press0,2.0)
                  << endl;
    }
  }
  tecplotfile.close();


  // Compare Yceq at different pressure
  ofstream Yceqfile;
  Yceqfile.open("Yceq.dat");
  Yceqfile << "VARIABLES = Z Yceq T" << endl;
  for (int i=0; i<nPress; ++i) {
    double myP = NewTable.GetCoordinate1(i);
    Yceqfile << "Zone T = \"P=" << myP << "\", I =" << nZm << endl;
    for (int j=0; j<nZm; ++j) {
      double myZ = NewTable.GetCoordinate2(j);
      Yceqfile << myZ << " "
               << NewTable.Lookup(myP, myZ, 0.0, 1.0, "PROG") << " "
               << NewTable.Lookup(myP, myZ, 0.0, 1.0, "T0") << endl;
    }
  }
  Yceqfile.close();

  // Compare Yceq at different Zvar
  Yceqfile.open("Yceq_Sz.dat");
  Yceqfile << "VARIABLES = Z Yceq T" << endl;
  double myP = NewTable.GetCoordinate1(nPress);
  for (int i=0; i<nZv; ++i) {
    double mySz = NewTable.GetCoordinate3(i);
    Yceqfile << "Zone T = \"Sz=" << mySz << "\", I =" << nZm << endl;
    for (int j=0; j<nZm; ++j) {
      double myZ = NewTable.GetCoordinate2(j);
      Yceqfile << myZ << " "
               << NewTable.Lookup(myP, myZ, mySz, 1.0, "PROG") << " "
               << NewTable.Lookup(myP, myZ, mySz, 1.0, "T0") << endl;
    }
  }
  Yceqfile.close();

  // Plot surfaces rho_min, rho_max (Zvar = 0.0)
  ofstream rhofile;
  rhofile.open("Rho_min_max.dat");
  rhofile << "VARIABLES = z c rho_min rho_max" << endl;
  rhofile << "ZONE I=" << nZm << " ,J=" << nC << " F=POINT" << endl;
  double myC, myZ, rho;
  double rho_min, rho_max;
  for (int i=0; i<nC; ++i) {
    myC = NewTable.GetCoordinate4(i);
    for (int j=0; j<nZm; ++j) {
      myZ = NewTable.GetCoordinate2(j);
      rho_min = 1.e6;
      rho_max = -1.0;
      for (int k=0; k<nPress; ++k) {
        double P = NewTable.GetCoordinate1(k);
        double r = NewTable.Lookup(P, myZ, 0.0, myC, "ROM");
        double temp = NewTable.Lookup(P, myZ, 0.0, myC, "T0");
        rho = P/r/temp;
        rho_min = min(rho_min,rho);
        rho_max = max(rho_max,rho);
      }
      rhofile << myZ << " " << myC << " " << rho_min << " " << rho_max << endl;
    }
  }
  rhofile.close();

  // Pressure scaling oh max OH
  ofstream ohfile;
  ohfile.open("max_Yoh.dat");
  ohfile << "VARIABLES = pressure YOH_max omega_max" << endl;
  double myZv;
  double oh_max, omega_max;
  for (int i=0; i<nPress; ++i) {
    myP = NewTable.GetCoordinate1(i);
    oh_max = -1.0;
    omega_max = -1.0;

    for (int j=0; j<nZm; ++j) {
      myZ = NewTable.GetCoordinate2(j);
      for (int k=0; k<nZv; ++k) {
        myZv = NewTable.GetCoordinate3(k);
        //for (int l=0; l<nC; ++l) {
          //myC = NewTable.GetCoordinate4(l);
          myC = 0.7;
          oh_max = max(oh_max,NewTable.Lookup(myP, myZ, myZv, myC,"OH"));
          omega_max = max(omega_max, NewTable.Lookup(myP, myZ, myZv, myC,"SRC_PROG"));
        //}
      }
    }
    ohfile << myP << " " << oh_max << " " << omega_max << endl;
  }
  ohfile.close();

  // Compare YOH at different pressure
  ofstream YOHfile;
  YOHfile.open("Maps_ZC_P.dat");
  YOHfile << "VARIABLES = Z C YOH HR rhoHR YH2O T YH2" << endl;
  for (int i=0; i<nPress; ++i) {
    double myP = NewTable.GetCoordinate1(i);
    YOHfile << "Zone T = \"P=" << myP << "\", I =" << nZm << " ,J=" << nC << " F=POINT" << endl;
    for (int k=0; k<nC; ++k) {
      double myC = NewTable.GetCoordinate4(k);

      for (int j=0; j<nZm; ++j) {
        double myZ = NewTable.GetCoordinate2(j);
        HR0 = NewTable.Lookup(myP, myZ, 0.0, myC,"HeatRelease");
        double rho = NewTable.Lookup(myP, myZ, 0.0, myC,"rho0");
        YOHfile << myZ << " " << myC << " "
               << NewTable.Lookup(myP, myZ, 0.0, myC, "OH")
               << HR0 << rho*HR0
               << NewTable.Lookup(myP, myZ, 0.0, myC, "H2O")
               << NewTable.Lookup(myP, myZ, 0.0, myC, "T0")
               << NewTable.Lookup(myP, myZ, 0.0, myC, "H2") << endl;
      }
    }
  }
  YOHfile.close();
  
  return 0;
}
