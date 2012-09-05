#ifndef CHEMTABLE3D_H
#define CHEMTABLE3D_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

#ifdef NO_ASSERT
#include <assert.h>
#endif

#include <myMem.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::fstream;
using std::ios;

class ChemtableVar {
public:
  ChemtableVar() {
    flag = 0;
    data = NULL;
  }
  int flag;
  string name;
  double  ***data;
};

class Chemtable3d {

private:

  double p_ref;
  double * Z;  
  double * aZv;  
  double * aC;  
//  double *** CA;
  int n; 
  int nzv; 
  int nC;
  vector<ChemtableVar> varVec;
  
public:
  
  Chemtable3d() {
    n = 0;
    nzv = 0;
    nC =0;
    Z = NULL;
    aZv = NULL;
    aC = NULL;
  }
  
  ~Chemtable3d() {
    if (Z != NULL) delete[] Z;
    if (aZv != NULL) delete[] aZv;
    if (aC != NULL) delete[] aC;
    for (int i = 0; i < varVec.size(); ++i) {
      if (varVec[i].data != NULL) 
         freeMem3D(varVec[i].data, 0, nC-1, 0, nzv-1, 0, n-1); 
         varVec[i].data = NULL;
    }
    //  freeMem3D(CA, 0, nC-1, 0, nzv-1, 0, n-1); 

  }

  void setSize(const int nvar,const int nC, const int nzv, const int n) {
    assert( this->nC == 0 );
    assert( this->n == 0 );
    assert( this->nzv == 0 );
    assert(Z == NULL);
    assert(aZv == NULL);
    assert(aC == NULL);
    assert(varVec.size() == 0);

    this->n = n;
    this->nzv = nzv;
    this->nC = nC;
    Z = new double[n];
    aZv = new double[nzv];
    aC = new double[nC];
    
    varVec.resize(nvar);
    for (int i = 0; i < nvar; ++i)
    getMem3D(&varVec[i].data,0, nC-1, 0, nzv-1, 0, n-1, "LoadTabC: myData", true);

//    getMem3D(&CA, 0, nC-1, 0, nzv-1, 0, n-1, "LoadTabC: myData", true);

  }

  int getN() const { return( n ); }
  int getNZV() const { return( nzv ); }
  int getNC() const { return( nC ); }

  void setPressure(const double p_ref) {
    this->p_ref = p_ref;
  }

  double getPressure() const {
    return( p_ref );
  }
  
  double getZ(const int i) const {
    assert( (i >= 0)&&(i < n) );
    assert( Z != NULL );
    return( Z[i] );
  }

  double getVar(const int tag,const int k, const int j, const int i) const {
    assert( (tag >= 0)&&(tag < varVec.size()) );
    assert( (i >= 0)&&(i < n) );
    assert( (j >= 0)&&(j < nzv) );
    assert( (k >= 0)&&(k < nC) );
    return( varVec[tag].data[k][j][i] );
  }

  void scaleVar(const int tag,const double scale) {
    assert( (tag >= 0)&&(tag < varVec.size()) );
     
    for (int k = 0; k < nC; ++k) {
    for (int j = 0; j < nzv; ++j) {
      for (int i = 0; i < n; ++i)
        varVec[tag].data[k][j][i] *= scale;
    }
    }
  }

  int getTag(const char * name) const {
    for (int i = 0; i < varVec.size(); ++i)
      if (varVec[i].name == name)
	return(i);
    cerr << "Error: could not match name in chemtable: \"" << name << "\"" << endl; 
    throw(-1);
  }
    
  void setVarName(const int tag,const char * name) {
    assert( (tag >= 0)&&(tag < varVec.size()) );
    assert( varVec[tag].flag == 0 );
    assert( strlen(name) < 32 );
    varVec[tag].flag = 1;
    varVec[tag].name = name;
  }

  void setZ(const int i,const double value) {
    assert( (i >= 0)&&(i < n) );
    assert( Z != NULL );
    Z[i] = value;
  }

  void setZv(const int i,const double value) {
    assert( (i >= 0)&&(i < nzv) );
    assert( aZv != NULL );
    aZv[i] = value;
  }

  void setaC(const int i,const double value) {
    assert( (i >= 0)&&(i < nC) );
    assert( aC != NULL );
    aC[i] = value;
  }
/*
  void setCA(const int k, const int j, const int i, const double value) {
    assert( (k >= 0)&&(k < nC) );
    assert( (j >= 0)&&(j < nzv) );
    assert( (i >= 0)&&(i < n) );
    assert( CA != NULL );

    CA[k][j][i] = value;
  }
*/
  
  void setVar(const int tag,const int k, const int j, const int i, const double value) {
    assert( (tag >= 0)&&(tag < varVec.size()) );
    assert( (i >= 0)&&(i < n) );
    assert( (j >= 0)&&(j < nzv) );
    assert( (k >= 0)&&(k < nC) );
    varVec[tag].data[k][j][i] = value;
  }

   #define lookupLinear lookup

  double lookupLinear(const int tag,const double thisZ, const double thisZv, const double thisaC) {

    assert( (tag >= 0)&&(tag < varVec.size()) );
    assert( n >= 2 );
    assert( nzv >= 2 );
    assert( nC >= 2 );
   // double *** data = varVec[tag].data;
    
    vector<double>num0(aC,aC+nC);
    vector<double>num1(aZv,aZv+nzv);
    vector<double>num2(Z,Z+n);

    for (int j=0; j<nzv;j++) {
       cout << "aZv is" <<aZv[j]<<endl;
     }

    vector<double>::const_iterator largest = 
     max_element(num2.begin(), num2.end());

    for (int j=0; j<nC;j++) {
       cout << "aC is" << aC[j]<<endl;
     }

    cout << "The largest number is " << *largest << endl;
    vector<double>::const_iterator tinest = 
     min_element( num2.begin(), num2.end() );
    cout << "The smallest number is " << *tinest << endl;

 
   // eight corner points
   // one interior 
   // six faces 
   // 12 edges
    int iZ, iZv, iC;
    if ( thisZ <= Z[0] ) {
        iZ = 0;
        if (thisZv<= aZv[0]) {
            iZv = 0;
            if(thisaC <= aC[0]) {
              iC =0;
              return( varVec[tag].data[0][0][0] );
            } else if (thisaC >= aC[nC-1]) {
              return( varVec[tag].data[nC-1][0][0] );
            }else {
               vector<double>::const_iterator one =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int j = one - num0.begin();
               double w1 = (aC[j]-thisaC)/(aC[j]-aC[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[j-1][0][0] + w2*varVec[tag].data[j][0][0] );
            }
        }else if (thisZv >=aZv[nzv-1]) {
            iZv = nzv-1;        

            if(thisaC <= aC[0]) {
              iC =0;
              return( varVec[tag].data[0][iZv][0] );
            } else if (thisaC >= aC[nC-1]) {
              return( varVec[tag].data[nC-1][iZv][0] );
            }else {
               vector<double>::const_iterator one =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int j = one - num0.begin();
               double w1 = (aC[j]-thisaC)/(aC[j]-aC[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[j-1][iZv][0] + w2*varVec[tag].data[j][iZv][0] );
            }

        }else { 

            if(thisaC <= aC[0]) {
               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();
               double w1 = (aZv[j]-thisZv)/(aZv[j]-aZv[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[0][j-1][0] + w2*varVec[tag].data[0][j][0] );
            } else if (thisaC >= aC[nC-1]) {
               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();
               double w1 = (aZv[j]-thisZv)/(aZv[j]-aZv[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[nC-1][j-1][0] + w2*varVec[tag].data[nC-1][j][0] );
            }else {   //bilinear interpolation over C & Zv

               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();

               vector<double>::const_iterator zero =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int k = zero - num0.begin();

               double deno = (aZv[j]-aZv[j-1])*(aC[k]-aC[k-1]); 
               double w1= (aZv[j]-thisZv) * (aC[k]-thisaC)/deno;
               double w2= (thisZv-aZv[j-1]) * (aC[k]-thisaC)/deno;
               double w3= (aZv[j]-thisZv) * (thisaC-aC[k-1])/deno;
               double w4= (thisZv-aZv[j-1]) * (thisaC-aC[k-1])/deno;
               return( w1*varVec[tag].data[k-1][j-1][0] + w2*varVec[tag].data[k-1][j][0] 
                   +w3*varVec[tag].data[k][j-1][0] + w4*varVec[tag].data[k][j][0]);
            }
               
        }
    }
    else if ( thisZ >= Z[n-1]) {
        iZ = n-1;
        if (thisZv<= aZv[0]) {
            if(thisaC <= aC[0]) {
              return( varVec[tag].data[0][0][n-1] );
            } else if (thisaC >= aC[nC-1]) {
              return( varVec[tag].data[nC-1][0][n-1] );
            }else {
               vector<double>::const_iterator one =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int j = one - num0.begin();
               double w1 = (aC[j]-thisaC)/(aC[j]-aC[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[j-1][0][n-1] + w2*varVec[tag].data[j][0][n-1] );
            }
        }else if (thisZv >=aZv[nzv-1]) {
            iZv = nzv-1;        

            if(thisaC <= aC[0]) {
              iC =0;
              return( varVec[tag].data[0][iZv][n-1] );
            } else if (thisaC >= aC[nC-1]) {
              return( varVec[tag].data[nC-1][iZv][n-1] );
            }else {
               vector<double>::const_iterator one =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int j = one - num0.begin();
               double w1 = (aC[j]-thisaC)/(aC[j]-aC[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[j-1][iZv][n-1] + w2*varVec[tag].data[j][iZv][n-1] );
            }

        }else { 

            if(thisaC <= aC[0]) {
               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();
               double w1 = (aZv[j]-thisZv)/(aZv[j]-aZv[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[0][j-1][n-1] + w2*varVec[tag].data[0][j][n-1] );
            } else if (thisaC >= aC[nC-1]) {
               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();
               double w1 = (aZv[j]-thisZv)/(aZv[j]-aZv[j-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[nC-1][j-1][n-1] + w2*varVec[tag].data[nC-1][j][n-1] );
            }else {   //bilinear interpolation over C & Zv

               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();

               vector<double>::const_iterator zero =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int k = zero - num0.begin();

               double deno = (aZv[j]-aZv[j-1])*(aC[k]-aC[k-1]); 
               double w1= (aZv[j]-thisZv) * (aC[k]-thisaC)/deno;
               double w2= (thisZv-aZv[j-1]) * (aC[k]-thisaC)/deno;
               double w3= (aZv[j]-thisZv) * (thisaC-aC[k-1])/deno;
               double w4= (thisZv-aZv[j-1]) * (thisaC-aC[k-1])/deno;
               return( w1*varVec[tag].data[k-1][j-1][n-1] + w2*varVec[tag].data[k-1][j][n-1] 
                   +w3*varVec[tag].data[k][j-1][n-1] + w4*varVec[tag].data[k][j][n-1]);
            }
               
        }



    } else {
      //  trilinear 

        if (thisZv<= aZv[0]) {
            if(thisaC <= aC[0]) {

               vector<double>::const_iterator one =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int i = one - num2.begin();
               double w1 = (Z[i]-thisZ)/(Z[i]-Z[i-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[0][0][i-1] + w2*varVec[tag].data[0][0][i] );

            } else if (thisaC >= aC[nC-1]) {
               vector<double>::const_iterator one =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int i = one - num2.begin();
               double w1 = (Z[i]-thisZ)/(Z[i]-Z[i-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[nC-1][0][i-1] + w2*varVec[tag].data[nC-1][0][i] );

            }else {
               vector<double>::const_iterator two =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int j = two - num2.begin();

               vector<double>::const_iterator zero =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int k = zero - num0.begin();

               double deno = (Z[j]-Z[j-1])*(aC[k]-aC[k-1]); 
               double w1= (Z[j]-thisZ) * (aC[k]-thisaC)/deno;
               double w2= (thisZ-Z[j-1]) * (aC[k]-thisaC)/deno;
               double w3= (Z[j]-thisZ) * (thisaC-aC[k-1])/deno;
               double w4= (thisZ-Z[j-1]) * (thisaC-aC[k-1])/deno;
               return( w1*varVec[tag].data[k-1][0][j-1] + w2*varVec[tag].data[k-1][0][j] 
                   +w3*varVec[tag].data[k][0][j-1] + w4*varVec[tag].data[k][0][j]);
            }
        }else if (thisZv >=aZv[nzv-1]) {
            if(thisaC <= aC[0]) {
               vector<double>::const_iterator one =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int i = one - num2.begin();
               double w1 = (Z[i]-thisZ)/(Z[i]-Z[i-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[0][nzv-1][i-1] + w2*varVec[tag].data[0][nzv-1][i] );

            } else if (thisaC >= aC[nC-1]) {
               vector<double>::const_iterator one =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int i = one - num2.begin();
               double w1 = (Z[i]-thisZ)/(Z[i]-Z[i-1]);
               double w2 = 1.0 - w1;
               return( w1*varVec[tag].data[nC-1][nzv-1][i-1] + w2*varVec[tag].data[nC-1][nzv-1][i] );

            }else {
               vector<double>::const_iterator two =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int j = two - num2.begin();

               vector<double>::const_iterator zero =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int k = zero - num0.begin();

               double deno = (Z[j]-Z[j-1])*(aC[k]-aC[k-1]); 
               double w1= (Z[j]-thisZ) * (aC[k]-thisaC)/deno;
               double w2= (thisZ-Z[j-1]) * (aC[k]-thisaC)/deno;
               double w3= (Z[j]-thisZ) * (thisaC-aC[k-1])/deno;
               double w4= (thisZ-Z[j-1]) * (thisaC-aC[k-1])/deno;
               return( w1*varVec[tag].data[k-1][nzv-1][j-1] + w2*varVec[tag].data[k-1][nzv-1][j] 
                   +w3*varVec[tag].data[k][nzv-1][j-1] + w4*varVec[tag].data[k][nzv-1][j]);
            }
        

        } else {
            if(thisaC <= aC[0]) {

               vector<double>::const_iterator two =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int j = two - num2.begin();

               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int k = one - num1.begin();

               double deno = (Z[j]-Z[j-1])*(aZv[k]-aZv[k-1]); 
               double w1= (Z[j]-thisZ) * (aZv[k]-thisZv)/deno;
               double w2= (thisZ-Z[j-1]) * (aZv[k]-thisZv)/deno;
               double w3= (Z[j]-thisZ) * (thisZv-aZv[k-1])/deno;
               double w4= (thisZ-Z[j-1]) * (thisZv-aZv[k-1])/deno;
               return( w1*varVec[tag].data[0][k-1][j-1] + w2*varVec[tag].data[0][k-1][j] 
                   +w3*varVec[tag].data[0][k][j-1] + w4*varVec[tag].data[0][k][j]);
            } else if (thisaC >= aC[nC-1]) {

               vector<double>::const_iterator two =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int j = two - num2.begin();

               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int k = one - num1.begin();

               double deno = (Z[j]-Z[j-1])*(aZv[k]-aZv[k-1]); 
               double w1= (Z[j]-thisZ) * (aZv[k]-thisZv)/deno;
               double w2= (thisZ-Z[j-1]) * (aZv[k]-thisZv)/deno;
               double w3= (Z[j]-thisZ) * (thisZv-aZv[k-1])/deno;
               double w4= (thisZ-Z[j-1]) * (thisZv-aZv[k-1])/deno;
               return( w1*varVec[tag].data[nC-1][k-1][j-1] + w2*varVec[tag].data[nC-1][k-1][j] 
                   +w3*varVec[tag].data[nC-1][k][j-1] + w4*varVec[tag].data[nC-1][k][j]);

            }else {
            // trilinear interpolation 
               cout << "doing trilinear interp" <<endl;
               vector<double>::const_iterator two =
                  upper_bound(num2.begin(),num2.end(),thisZ);              
               int i = two - num2.begin();

               vector<double>::const_iterator one =
                  upper_bound(num1.begin(),num1.end(),thisZv);              
               int j = one - num1.begin();

               vector<double>::const_iterator zero =
                  upper_bound(num0.begin(),num0.end(),thisaC);              
               int k = zero - num0.begin();

               double deno = (Z[i]-Z[i-1])*(aZv[j]-aZv[j-1])*(aC[k]-aC[k-1]); 

               double w111 = (Z[i]-thisZ)*(aZv[j]-thisZv)*(aC[k]-thisaC)/deno;   
               double w112 = (Z[i]-thisZ)*(aZv[j]-thisZv)*(thisaC-aC[k-1])/deno;   
               double w211 = (thisZ-Z[i-1])*(aZv[j]-thisZv)*(aC[k]-thisaC)/deno;   
               double w212 = (thisZ-Z[i-1])*(aZv[j]-thisZv)*(thisaC-aC[k-1])/deno;   
               double w121 = (Z[i]-thisZ)*(thisZv-aZv[j-1])*(aC[k]-thisaC)/deno;   
               double w122 = (Z[i]-thisZ)*(thisZv-aZv[j-1])*(thisaC-aC[k-1])/deno;   
               double w221 = (thisZ-Z[i-1])*(thisZv-aZv[j-1])*(aC[k]-thisaC)/deno;   
               double w222 = (thisZ-Z[i-1])*(thisZv-aZv[j-1])*(thisaC-aC[k-1])/deno;   
               return( w111*varVec[tag].data[k-1][j-1][i-1] + w112*varVec[tag].data[k-1][j-1][i] 
                   +w211*varVec[tag].data[k][j-1][i-1] + w212*varVec[tag].data[k][j-1][i]
                   +w121*varVec[tag].data[k-1][j][i-1] + w122*varVec[tag].data[k-1][j][i]
                   +w221*varVec[tag].data[k][j][i-1] + w222*varVec[tag].data[k][j][i]);
            }
        }
    }
  }


  void write(const string& filename,const bool verbose = false) {
    write(filename.c_str(),verbose);
  }

  void write(const char * filename,const bool verbose = false) {

    if (verbose)
      cout << "writing chemtable1d to file: " << filename << endl;
    
    fstream fout(filename, ios::out|ios::binary|ios::trunc);
    
    if (fout.fail()) {
      cerr << "Error: could not open file: " << filename << endl;
      throw(-1);
    }
    
    if (verbose)
      cout << " > nvar,n: " << varVec.size() << " " << n << endl;
    
    int ivar[9];
    ivar[0] = 123581321; // magic number
    ivar[1] = 3; // chemtable dimension
    ivar[2] = varVec.size(); // nvars + ... + SRC_PROG
    ivar[3] = n;
    ivar[4] = nzv;
    ivar[5] = nC;
    ivar[6] = 0;
    ivar[7] = 0;
    ivar[8] = 0;
    fout.write((char*)ivar,sizeof(int)*9);

    double dvar[9];
    dvar[0] = p_ref;
    dvar[1] = 0.0;
    dvar[2] = 0.0;
    dvar[3] = 0.0;
    dvar[4] = 0.0;
    dvar[5] = 0.0;
    dvar[6] = 0.0;
    dvar[7] = 0.0;
    dvar[8] = 0.0;

    double tabvar[n];
    fout.write((char*)dvar,sizeof(double)*9);
    // Z...
    fout.write((char*)Z,sizeof(double)*n);
    // Zv...
    fout.write((char*)aZv,sizeof(double)*nzv);
    // C...
    for (int i=0; i< nC; i++) { 
        cout<<"aC write is " << aC[i] << endl;
    }
    fout.write((char*)aC,sizeof(double)*nC);

    // vars...
    for (int i = 0; i < varVec.size(); ++i) {
      char name[32];
      strcpy(name,varVec[i].name.c_str());
      if (verbose)
	cout << " > var: \"" << name << "\"" << endl;
      fout.write(name,32);
      for (int k=0; k <nC; ++k) {
         for (int j=0; j <nzv; ++j){ 
             for (int l=0; l<n;++l) tabvar[l]=varVec[i].data[k][j][l]; 
             fout.write((char*)tabvar,sizeof(double)*n);
         }
      }
    }
    
    fout.close();

  }

  void writetec(const string& filename,const bool verbose = false) {
    writetec(filename.c_str(),verbose);
  }

  void writetec(const char * filename,const bool verbose = false) {

    if (verbose)
      cout << "writing chemtable1d to file: " << filename << endl;
    
    ofstream outf;

    outf.open(filename);

    if (outf.fail()) {
      cerr << "Error: could not open file: " << filename << endl;
      throw(-1);
    }
    
    if (verbose)
      cout << " > nvar,n: " << varVec.size() << " " << n << endl;
    
    int ivar[9];
    ivar[0] = 123581321; // magic number
    ivar[1] = 3; // chemtable dimension
    ivar[2] = varVec.size()+1; // nvars
    ivar[3] = n;
    ivar[4] = nzv;
    ivar[5] = nC;
    ivar[6] = 0;
    ivar[7] = 0;
    ivar[8] = 0;
    //outf << ivar[0] << ivar[1] << ivar[2] <<endl;

    double dvar[9];
    dvar[0] = p_ref;
    dvar[1] = 0.0;
    dvar[2] = 0.0;
    dvar[3] = 0.0;
    dvar[4] = 0.0;
    dvar[5] = 0.0;
    dvar[6] = 0.0;
    dvar[7] = 0.0;
    dvar[8] = 0.0;
   // outf<<dvar[0] << dvar[1] <<dvar[2] <<endl;

    // Z...
  //  outf<<Z[0] << Z[1] << Z[2]<<endl;

    // vars...
// TECPLOT variables
//    outf<< "Z "<<"  " <<varVec[0].name.c_str() <<"  "<< varVec[1].name.c_str()<<endl;
  outf << "TITLE=\"Chemistry table\"" << endl;
  outf << "VARIABLES=\"Zm\"  \"ZVar\"  \"Progress\" ";


  for (int j=0; j<varVec.size(); ++j)
    outf << "\"" << varVec[j].name.c_str()<< "\" ";
    outf << endl;
   outf<< " ZONE T=\"Floor\", I=" <<nC 
              <<" J=" << nzv <<" K=" << n <<" F=POINT " <<endl;

    for (int k=0; k<n; ++k){
    for (int kk=0; kk<nzv; ++kk){
    for (int kc=0; kc<nC; ++kc){
       for (int i = 0; i < varVec.size(); ++i) {
         char name[32];
         strcpy(name,varVec[i].name.c_str());
         if (verbose) {
	    cout << " > var: \"" << name << "\"" << endl;
         }
        }
           outf<< Z[k] <<"  ";
           outf<< aZv[kk] <<"  ";
          // outf<< CA[kc][kk][k] <<"  ";
           outf<< aC[kc] <<"  ";
           for (int j=0; j<varVec.size(); ++j) 
           outf <<varVec[j].data[kc][kk][k] <<"  ";
           outf<<endl;
    }
    }
    }
    outf.close();

  }

  void read(const string& filename,const bool verbose = false) {
    read(filename.c_str(),verbose);
  }
    
  void read(const char * filename,const bool verbose = false) {
    
    if (verbose)
      cout << "reading chemtable1d from file: " << filename << endl;
    
    fstream fin(filename,ios::in|ios::binary);

    if (fin.fail()) {
      cerr << "Error: could not open file: " << filename << endl;
      throw(-1);
    }
    
    int ivar[9];
    fin.read((char*)ivar,sizeof(int)*9);
    if (ivar[0] != 123581321) {
      if (verbose)
	cout << " > file requires byte swap: not implemented yet" << endl;
      throw(-1);
    }
    assert( ivar[1] == 3 ); // chemtable dimension
    setSize(ivar[2],ivar[5],ivar[4],ivar[3]);

    if (verbose)
      cout << " > nvar,n: " << varVec.size() << " " << n << endl;

    double dvar[9];
    double tabvar[n];
    fin.read((char*)dvar,sizeof(double)*9);
    setPressure(dvar[0]);
    
    if (verbose)
      cout << " > p_ref: " << p_ref << endl;
    
    // Z...
    fin.read((char*)Z,sizeof(double)*n);
    fin.read((char*)aZv,sizeof(double)*nzv);
    fin.read((char*)aC,sizeof(double)*nC);

    for (int j=0; j<nC;j++) {
         cout << "aC read is " <<aC[j]<<endl;
    }

    // vars...
    for (int i = 0; i < varVec.size(); ++i) {
      char name[32];
      fin.read(name,32);
      if (verbose)
	cout << " > var: \"" << name << "\" range: ";
      setVarName(i,name);
      // fin.read((char*)varVec[i].data,sizeof(double)*nC*nzv*n);
      for (int k=0; k <nC; ++k){ 
         for (int j=0; j <nzv; ++j){ 
             fin.read((char*)tabvar,sizeof(double)*n);
             for (int l=0; l<n;++l) varVec[i].data[k][j][l] = tabvar[l]; 
         }
      }

      if (verbose) {
	double minval = varVec[i].data[0][0][0];
	double maxval = varVec[i].data[0][0][0];
	for (int ii = 1; ii < n; ++ii) {
	  minval = min( varVec[i].data[0][0][ii], minval);
	  maxval = max( varVec[i].data[0][0][ii], maxval);
	}
	cout << minval << " " << maxval << endl;
      }
    }
    
    fin.close();

  }

};
    
#endif

