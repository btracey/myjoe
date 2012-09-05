#ifndef UGP_H
#define UGP_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Adt.h"
#include "Param.h"

#include "MpiStuff.h"
using namespace MpiStuff;

#include "MiscUtils.h"
using namespace MiscUtils;

#include "Gp.h"

// #############################################
// #############################################
// #############################################

#define FOR_INO for (int ino = 0; ino < nno; ++ino)
#define FOR_IFA for (int ifa = 0; ifa < nfa; ++ifa)
#define FOR_IFA_B for (int ifa = 0; ifa < nfa_b; ++ifa)
#define FOR_IFA_NONB for (int ifa = nfa_b; ifa < nfa; ++ifa)
#define FOR_ICV for (int icv = 0; icv < ncv; ++icv)
#define FOR_I2 for (int i = 0; i < 2; ++i)
#define FOR_I3 for (int i = 0; i < 3; ++i)
#define FOR_J3 for (int j = 0; j < 3; ++j)
#define FOR_K3 for (int k = 0; k < 3; ++k)
#define FOR_ICV_B for (int icv = 0; icv < ncv_b; ++icv)
#define FOR_ICV_NONB for (int icv = ncv_b; icv < ncv; ++icv)

#define NNOF_MAX 8
#define NNOC_MAX 81
#define NFOC_MAX 24 // 4x6
#define ZONE_NAME_LEN         64

// face zone kind's...
// reserve -1 for unknown
#define FA_ZONE_UNKNOWN          -1
#define FA_ZONE_PERIODIC_UNKNOWN -2

#define FA_ZONE_BOUNDARY          1

#define FA_ZONE_PERIODIC_CART     2
#define FA_ZONE_PERIODIC_CYL_X    3
#define FA_ZONE_PERIODIC_CYL_Y    4
#define FA_ZONE_PERIODIC_CYL_Z    5

#define FA_ZONE_INTERNAL          6

// range for all and periodic zones (tommie needs this?)...
#define FA_ZONE_FIRST             1
#define FA_ZONE_LAST              6
#define FA_ZONE_PERIODIC_FIRST    2
#define FA_ZONE_PERIODIC_LAST     5

// cv zone kinds...
#define CV_ZONE_UNKNOWN          -1
#define CV_ZONE_FLUID             1

// types of elements...
// limit this to 8 types max for use in adaptation
#define HEX_TYPE     0
#define PRISM_TYPE   1
#define PYRAMID_TYPE 2
#define TET_TYPE     3
#define UNKNOWN_TYPE 4

// types of faces...
// limit this to 4 types max for use in adaptation
#define QUAD_TYPE    0
#define TRI_TYPE     1

class ScalarBcData {

private:

  int kind;

public:

  string name;
  double *data;
  int * flag;

  ScalarBcData(const string& name) {
    this->name = name;
    kind = -1;
    flag = NULL;
    data = NULL;
  }

  int getKind() const {
    return (kind);
  }
  void setKind(const int kind) {
    this->kind = kind;
  }

};

class VectorBcData {

private:

  int kind;

public:

  string name;
  double (*data)[3];
  int * flag;

  VectorBcData(const string& name) {
    this->name = name;
    kind = -1;
    flag = NULL;
    data = NULL;
  }

  int getKind() const {
    return (kind);
  }
  void setKind(const int kind) {
    this->kind = kind;
  }

};

class FaZone {

private:

  char name[ZONE_NAME_LEN];
  int index, periodicIndex, kind;

public:

  bool ib_flag;

  // for boundary conditions...
  list<ScalarBcData> scalarBcDataList;
  list<VectorBcData> vectorBcDataList;

  int flag;

  // lists of faces...
  int ifa_f, ifa_l;
  int ifa_i_f, ifa_i_l;

  // periodicity...
  double periodic_data[3];
  int bits;

  // overload operator for sorting face zones...
  // sort based on kind, where kind is ordered as
  // above...
  bool operator <(FaZone &fzcompare) {
    return ((kind<fzcompare.kind)||((kind==fzcompare.kind)&&(strcmp(name, fzcompare.name)<0)));
  }

  FaZone() {

    ib_flag = false;

    strcpy(name, "UNKNOWN");
    index = periodicIndex = 0;
    kind = FA_ZONE_UNKNOWN;
    flag = 0;

    ifa_f = 0;
    ifa_l = -1;
    ifa_i_f = 0;
    ifa_i_l = -1;

    bits = 0;
    FOR_I3
      periodic_data[i] = 0.0;

  }

  int getKind() const {
    return (kind);
  }
  void setKind(const int kind) {
    this->kind = kind;
  }

  bool isPeriodic() const {
    return ((kind>=FA_ZONE_PERIODIC_FIRST)&&(kind<=FA_ZONE_PERIODIC_LAST));
  }

  int getIndex() const {
    return (index);
  }
  void setIndex(const int index) {
    this->index = index;
  }

  int getPeriodicIndex() const {
    return (periodicIndex);
  }
  void setPeriodicIndex(const int index) {
    periodicIndex = index;
  }

  int getBits() const {
    return (bits);
  }
  void setBits(const int new_bits) {
    bits = new_bits;
  }

  char * getName() {
    return (name);
  }
  string getNameString() const {
    return (string(name));
  }

  void setName(const string newName) {
    assert(newName.length()<=ZONE_NAME_LEN);
    strcpy(name, newName.c_str());
  }

  void setName(char * new_name) {
    assert(strlen(new_name)<=ZONE_NAME_LEN);
    strcpy(name, new_name);
  }

  void setPeriodicData(double * data) {
    for (int id = 0; id<3; id++)
      periodic_data[id] = data[id];
  }

  void getPeriodicData(double * data) {
    for (int id = 0; id<3; id++)
      data[id] = periodic_data[id];
  }

  void scalePeriodicData(const double factor) {
    assert(kind==FA_ZONE_PERIODIC_CART);
    for (int id = 0; id<3; id++)
      periodic_data[id] *= factor;
  }

  void shearXYPeriodicData(const double factor) {
    assert(kind==FA_ZONE_PERIODIC_CART);
    periodic_data[0] += factor*periodic_data[1];
  }

  void scalePeriodicData(const double * factor) {
    assert(kind==FA_ZONE_PERIODIC_CART);
    for (int id = 0; id<3; id++)
      periodic_data[id] *= factor[id];
  }

  void rotatePeriodicDataX(const double degrees) {
    assert(kind==FA_ZONE_PERIODIC_CART);
    double cos_t = cos(M_PI*degrees/180.0);
    double sin_t = sin(M_PI*degrees/180.0);
    double dy = cos_t*periodic_data[1]-sin_t*periodic_data[2];
    double dz = cos_t*periodic_data[2]+sin_t*periodic_data[1];
    periodic_data[1] = dy;
    periodic_data[2] = dz;
  }

  void rotatePeriodicDataZ(const double degrees) {
    assert(kind==FA_ZONE_PERIODIC_CART);
    double cos_t = cos(M_PI*degrees/180.0);
    double sin_t = sin(M_PI*degrees/180.0);
    double dx = cos_t*periodic_data[0]-sin_t*periodic_data[1];
    double dy = cos_t*periodic_data[1]+sin_t*periodic_data[0];
    periodic_data[0] = dx;
    periodic_data[1] = dy;
  }

  void dump() {
    cout<<"FaZone::dump(), name: \""<<name<<"\", index, kind, bits: "<<" "<<index<<" "<<kind<<" ";
    for (int ii = 0; ii<12; ii++) {
      if (bits&(1<<ii)) {
        cout<<"1";
      }
      else {
        cout<<"0";
      }
    }
    cout<<endl;
  }

  VectorBcData * getVectorBcData(const string& name) {
    for (list<VectorBcData>::iterator i = vectorBcDataList.begin(); i!=vectorBcDataList.end(); i++) {
      if (i->name==name) return (&(*i));
    }
    cerr<<"Error: getVectorBcData: could not find data: "<<name<<endl;
    throw(-1);
  }

  ScalarBcData * getScalarBcData(const string& name) {
    for (list<ScalarBcData>::iterator i = scalarBcDataList.begin(); i!=scalarBcDataList.end(); i++) {
      if (i->name==name) return (&(*i));
    }
    cerr<<"Error: getScalarBcData: could not find data: "<<name<<endl;
    throw(-1);
  }

};

class CvZone {

private:

  char name[ZONE_NAME_LEN];
  int index, kind;
  int icv_f, icv_l;

public:

  int flag;

  CvZone() {
    strcpy(name, "UNKNOWN");
    index = 0;
    kind = CV_ZONE_UNKNOWN;
    flag = 0;
  }

  int getKind() const {
    return (kind);
  }
  void setKind(const int kind) {
    this->kind = kind;
  }

  int getIndex() const {
    return (index);
  }
  void setIndex(const int index) {
    this->index = index;
  }

  char * getName() {
    return (name);
  }

  void setName(const string newName) {
    assert(newName.length()<=ZONE_NAME_LEN);
    strcpy(name, newName.c_str());
  }

  void setName(char * new_name) {
    assert(strlen(new_name)<=ZONE_NAME_LEN);
    strcpy(name, new_name);
  }

  void dump() {
    cout<<"CvZone::dump(), name: \""<<name<<"\", index, kind: "<<" "<<index<<" "<<kind<<endl;
  }

};

// various datatypes...
#define VALUE_DATA 0
#define CV_DATA    1
#define FA_DATA    2
#define NO_DATA    3

#define DATA_NAME_LEN 32

class Data {

private:

  char name[DATA_NAME_LEN];
  int datatype, flag;

public:

  Data(const char * name, int datatype) {
    strcpy(this->name, name);
    this->datatype = datatype;
    flag = 0;
  }

  Data(const char * name) {
    strcpy(this->name, name);
    this->datatype = VALUE_DATA;
    flag = 0;
  }

  char * getName() {
    return (name);
  }

  void setName(const char * new_name) {
    assert(strlen(new_name)<=DATA_NAME_LEN);
    strcpy(name, new_name);
  }

  int getDatatype() const {
    return (datatype);
  }

  int getFlag() const {
    return (flag);
  }

  int checkFlag() const {
    return (flag!=0);
  }

  void setFlag() {
    setFlag(1);
  }

  void clearFlag() {
    setFlag(0);
  }

  void setFlag(const int value) {
    flag = value;
  }

};

class DoubleValue: public Data {

public:

  double * ptr;

  DoubleValue(double &value, const char * name) :
    Data(name) {
    ptr = &value;
  }

};

class IntValue: public Data {

public:

  int * ptr;

  IntValue(int &value, const char * name) :
    Data(name) {
    ptr = &value;
  }

};

class IntScalar: public Data {

public:

  int ** ptr;

  IntScalar(int *&scalar, const char * name, int datatype) :
    Data(name, datatype) {
    ptr = &scalar;
  }

  ~IntScalar() {
    delete[] (*ptr);
  }

};

class DoubleScalar: public Data {

public:

  double ** ptr;

  DoubleScalar(double *&scalar, const char * name, int datatype) :
    Data(name, datatype) {
    ptr = &scalar;
  }

  ~DoubleScalar() {
    delete[] (*ptr);
  }

};

class DoubleVector: public Data {

public:

  double (**ptr)[3];

  DoubleVector(double(*&vector)[3], const char * name, int datatype) :
    Data(name, datatype) {
    ptr = &vector;
  }

  ~DoubleVector() {
    delete[] (*ptr);
  }

};

class DoubleTensor: public Data {

public:

  double (**ptr)[3][3];

  DoubleTensor(double(*&tensor)[3][3], const char * name, int datatype) :
    Data(name, datatype) {
    ptr = &tensor;
  }

  ~DoubleTensor() {
    delete[] (*ptr);
  }

};

// at times it is neccessary to sort pairs or triples of ints according to one
// of their indices. For example, faces need to be matched based on global index
// during the building of the face communication structure...

typedef struct {
  int i1, i2;
} IntPair;

class IntPairCompare1: public std::binary_function<IntPair&, IntPair&, bool> {
public:
  bool operator()(const IntPair& a, const IntPair& b) {
    return (a.i1<b.i1);
  }
};

// same as above, but produces a stable sort behavior...

class IntPairCompare12: public std::binary_function<IntPair&, IntPair&, bool> {
public:
  bool operator()(const IntPair& a, const IntPair& b) {
    return ((a.i1<b.i1)||((a.i1==b.i1)&&(a.i2<b.i2)));
  }
};

typedef struct {
  int i1, i2, i3;
} IntTriple;

class IntTripleCompare1: public std::binary_function<IntTriple&, IntTriple&, bool> {
public:
  bool operator()(const IntTriple& a, const IntTriple& b) {
    return (a.i1<b.i1);
  }
};

class IntTripleCompare3: public std::binary_function<IntTriple&, IntTriple&, bool> {
public:
  bool operator()(const IntTriple& a, const IntTriple& b) {
    return (a.i3<b.i3);
  }
};

class IntTripleCompare12: public std::binary_function<IntTriple&, IntTriple&, bool> {
public:
  bool operator()(const IntTriple& a, const IntTriple& b) {
    return ((a.i1<b.i1)||((a.i1==b.i1)&&(a.i2<b.i2)));
  }
};

class IntTripleCompare23: public std::binary_function<IntTriple&, IntTriple&, bool> {
public:
  bool operator()(const IntTriple& a, const IntTriple& b) {
    return ((a.i2<b.i2)||((a.i2==b.i2)&&(a.i3<b.i3)));
  }
};

class IntTripleCompare123: public std::binary_function<IntTriple&, IntTriple&, bool> {
public:
  bool operator()(const IntTriple& a, const IntTriple& b) {
    return ((a.i1<b.i1)||((a.i1==b.i1)&&(a.i2<b.i2))||((a.i1==b.i1)&&(a.i2==b.i2)&&(a.i3<b.i3)));
  }
};

typedef struct {
  int i1, i2, i3, i4;
} IntQuad;

class IntQuadCompare1: public std::binary_function<IntQuad&, IntQuad&, bool> {
public:
  bool operator()(const IntQuad& a, const IntQuad& b) {
    return (a.i1<b.i1);
  }
};

class IntQuadCompare4: public std::binary_function<IntQuad&, IntQuad&, bool> {
public:
  bool operator()(const IntQuad& a, const IntQuad& b) {
    return (a.i4<b.i4);
  }
};

class IntQuadCompare12: public std::binary_function<IntQuad&, IntQuad&, bool> {
public:
  bool operator()(const IntQuad& a, const IntQuad& b) {
    return ((a.i1<b.i1)||((a.i1==b.i1)&&(a.i2<b.i2)));
  }
};

class IntQuadCompare123: public std::binary_function<IntQuad&, IntQuad&, bool> {
public:
  bool operator()(const IntQuad& a, const IntQuad& b) {
    return ((a.i1<b.i1)||((a.i1==b.i1)&&(a.i2<b.i2))||((a.i1==b.i1)&&(a.i2==b.i2)&&(a.i3<b.i3)));
  }
};

// Ugp implements some generic solvers. These
// are the modes suported in some cases...

/// \def solver modes:
#define RELATIVE_RESIDUAL 1
#define ABSOLUTE_RESIDUAL 2

#define UGP_NNOC_MAX 80

class Ugp: public Gp {

  // the various supported filters get access to the private data...
  friend class CdpFilter;

private:

  // the Ugp registration allows for data-free registration
  // using these lists to store the actual data that would normally
  // be stored by the calling process...
  list<int> intList;
  list<double> doubleList;
  list<double*> doublePtrList;
  list<double(*)[3]> doublePtr3List;

public:

  int ncv_b;
  int ncv_bp;
  int ncv; ///< number of local control volumes

  int nno; ///< number of local nodes
  int nno_i; // internal nodes - i.e. nodes that have all their dist-1 nbrs local

  // face subdivisions...
  int nfa_b; // 1. boundary zones         [0:nfa_b-1]
  int nfa_bp; // 2. periodic boundary      [nfa_b:nfa_bp-1]
  int nfa_bpi; // 3. processor boundary     [nfa_bp:nfa_bpi-1]
  int nfa; // 4. internal               [nfa_bpi:nfa]

  int * noofa_i; ///< [nfa+1] index range for node-of-face CSR structure
  int noofa_s;
  int * noofa_v; ///< [noofa_i[nfa]] values for node-of-face CSR structure
  int (*cvofa)[2]; ///< [nfa][2] cv-of-face: cvofa[*][1] is -1 on boundary faces

  int * faocv_i;
  int faocv_s;
  int * faocv_v;

  int * noocv_i;
  int noocv_s;
  int * noocv_v;

  // partition info...
  int * cv_part; // the cv partition - zero-indexed, possibly read in from file
  int npart; // current number of partitions in cv_part.

  // face and cv zone lists...
  list<FaZone> faZoneList;
  list<CvZone> cvZoneList;

  int * noora;
  int * faora;
  int * cvora;

  int * cv_flag; /* [ncv] flag for cv's */
  int * fa_flag; /* [nfa] flag for faces */
  int * no_flag;

  double (*x_no)[3];

	// (begin) Multigrid variables
  int ncv_mgLevel1;
	int nfa_b_mgLevel1;
	int nfa_mgLevel1;
	int (*cvofa_mgLevel1)[2];
	list<FaZone> faZoneList_mgLevel1;
  int * faocv_i_mgLevel1;
  int faocv_s_mgLevel1;
  int * faocv_v_mgLevel1;
	int nfa_bp_mgLevel1;
  int nfa_bpi_mgLevel1;
	
	
	int ncv_mgLevel2;
	int nfa_b_mgLevel2;
	int nfa_mgLevel2;
	int (*cvofa_mgLevel2)[2];
	list<FaZone> faZoneList_mgLevel2;
  int * faocv_i_mgLevel2;
  int faocv_s_mgLevel2;
  int * faocv_v_mgLevel2;
	int nfa_bp_mgLevel2;
  int nfa_bpi_mgLevel2;
	
	int ncv_mgLevel3;
	int nfa_b_mgLevel3;
	int nfa_mgLevel3;
	int (*cvofa_mgLevel3)[2];
	list<FaZone> faZoneList_mgLevel3;
  int * faocv_i_mgLevel3;
  int faocv_s_mgLevel3;
  int * faocv_v_mgLevel3;
	int nfa_bp_mgLevel3;
  int nfa_bpi_mgLevel3;
	// (end) Multigrid variables
	
  // data registration...
  list<IntValue> intValueList;
  list<IntScalar> intScalarList;
  list<DoubleValue> doubleValueList;
  list<DoubleScalar> doubleScalarList;
  list<DoubleVector> doubleVectorList;
  list<DoubleTensor> doubleTensorList;

  // registered int for cv groupings for "linelets"...
  int * cvGroup;

  void buildXora(int * &Xora, int nX) {
    if (Xora==NULL) Xora = new int[mpi_size+1];
    MPI_Allgather(&nX, 1, MPI_INT, Xora+1, 1, MPI_INT, mpi_comm);
    Xora[0] = 0;
    for (int i = 0; i<mpi_size; i++)
      Xora[i+1] += Xora[i];
  }

  // allow addition of another min initialized ugp...
  Ugp& operator+=(Ugp& ugp) {

    // not set up to handle registered data yet...
    assert(intValueList.size()==0);
    assert(intScalarList.size()==0);
    assert(doubleValueList.size()==0);
    assert(doubleScalarList.size()==0);
    assert(doubleVectorList.size()==0);
    assert(doubleTensorList.size()==0);

    // update the counts...
    int nno_old = nno;
    nno += ugp.nno;
    int nfa_old = nfa;
    nfa += ugp.nfa;
    int ncv_old = ncv;
    ncv += ugp.ncv;

    // because the global numbering is striped across all processors, it
    // is neccessary to update all indices (existing and new)...

    // nodes...
    int * noora_old = NULL;
    buildXora(noora_old, nno_old);
    int * noora_ugp = NULL;
    buildXora(noora_ugp, ugp.nno);

    // faces...
    int * faora_old = NULL;
    buildXora(faora_old, nfa_old);
    int * faora_ugp = NULL;
    buildXora(faora_ugp, ugp.nfa);

    // cvs...
    int * cvora_old = NULL;
    buildXora(cvora_old, ncv_old);
    int * cvora_ugp = NULL;
    buildXora(cvora_ugp, ugp.ncv);

    // ======================
    // nodes...
    // ======================
    if (ugp.nno>0) {
      // ----------------------------
      // node coords x_no...
      // ----------------------------
      resize(x_no, nno_old, nno);
      for (int ino = 0; ino<ugp.nno; ++ino)
        FOR_I3
          x_no[ino+nno_old][i] = ugp.x_no[ino][i];
      // ----------------------------
      // no_flag: just copy it over...
      // ----------------------------
      resize(no_flag, nno_old, nno);
      for (int ino = 0; ino<ugp.nno; ++ino)
        no_flag[ino+nno_old] = ugp.no_flag[ino];
    }

    // ======================
    // faces...
    // ======================
    // adjust the existing global node indexing...
    if (nfa_old>0) {
      assert(noofa_s==noofa_i[nfa_old]);
      for (int nof = 0; nof<noofa_s; ++nof) {
        int ino_old = noofa_v[nof];
        assert((ino_old>=0)&&(ino_old<noora_old[mpi_size]));
        // figure out where this node lives based on its global index...
        int rank = 0;
        while (ino_old>=noora_old[rank+1])
          ++rank;
        // update its index...
        noofa_v[nof] = ino_old+noora_ugp[rank];
      }
    }
    // adjust the global cv indexing in cvofa...
    int * fa_check = new int[nfa_old];
    for (int ifa = 0; ifa<nfa_old; ++ifa)
      fa_check[ifa] = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); ++zone) {
      int zone_index = zone->getIndex();
      for (int ifa = 0; ifa<nfa_old; ++ifa) {
        if (fa_flag[ifa]==zone_index) {
          // check that every face gts touched once...
          assert(fa_check[ifa]==0);
          fa_check[ifa] = 1;
          // icv0 is always a global cv number...
          int icv0_old = cvofa[ifa][0];
          assert((icv0_old>=0)&&(icv0_old<cvora_old[mpi_size]));
          int rank = 0;
          while (icv0_old>=cvora_old[rank+1])
            ++rank;
          // update its index...
          cvofa[ifa][0] = icv0_old+cvora_ugp[rank];
          // icv1 may be -1, in which case it remains unchanged. Another possibility is
          // that it is a global face indexing associated with a periodic pairing...
          int icv1_old = cvofa[ifa][1];
          if (icv1_old>=0) {
            if (zone->isPeriodic()) {
              // icv1_old is a global face index...
              assert((icv1_old>=0)&&(icv1_old<faora_old[mpi_size]));
              rank = 0;
              while (icv1_old>=faora_old[rank+1])
                ++rank;
              cvofa[ifa][1] = icv1_old+faora_ugp[rank];
            }
            else {
              // icv1_old is a global cv index...
              assert((icv1_old>=0)&&(icv1_old<cvora_old[mpi_size]));
              rank = 0;
              while (icv1_old>=cvora_old[rank+1])
                ++rank;
              cvofa[ifa][1] = icv1_old+cvora_ugp[rank];
            }
          }
          else {
            // should be -1...
            assert(icv1_old==-1);
          }
        }
      }
    }
    // make sure all faces were touched...
    for (int ifa = 0; ifa<nfa_old; ++ifa)
      assert(fa_check[ifa]==1);
    delete[] fa_check;

    // add new global node numbers...
    if (ugp.nfa>0) {
      // ----------------------------
      // the noofa_i/v CSR struct requires
      // special initialization the first time...
      // ----------------------------
      if (nfa_old==0) {
        assert(noofa_s==0);
        assert(noofa_i==NULL);
        assert(noofa_v==NULL);
        assert(ugp.noofa_s==ugp.noofa_i[ugp.nfa]);
        noofa_s = ugp.noofa_s;
        resize(noofa_i, nfa+1);
        resize(noofa_v, noofa_s);
        noofa_i[0] = 0;
      }
      else {
        assert(noofa_s==noofa_i[nfa_old]);
        assert(ugp.noofa_s==ugp.noofa_i[ugp.nfa]);
        int noofa_s_old = noofa_s;
        noofa_s += ugp.noofa_s;
        resize(noofa_i, nfa_old+1, nfa+1);
        resize(noofa_v, noofa_s_old, noofa_s);
      }
      for (int ifa = 0; ifa<ugp.nfa; ++ifa) {
        noofa_i[ifa+nfa_old+1] = noofa_i[ifa+nfa_old]+ugp.noofa_i[ifa+1]-ugp.noofa_i[ifa];
        for (int nof = ugp.noofa_i[ifa]; nof<ugp.noofa_i[ifa+1]; ++nof) {
          int ino_ugp = ugp.noofa_v[nof];
          assert((ino_ugp>=0)&&(ino_ugp<noora_ugp[mpi_size]));
          int rank = 0;
          while (ino_ugp>=noora_ugp[rank+1])
            ++rank;
          int nof_new = nof-ugp.noofa_i[ifa]+noofa_i[ifa+nfa_old];
          noofa_v[nof_new] = ino_ugp+noora_old[rank+1];
        }
      }
      // ----------------------------
      // fa_flag: contains the faZone index, and
      // gets set below. Set it to -1 for now...
      // ----------------------------
      resize(fa_flag, nfa_old, nfa);
      for (int ifa = 0; ifa<ugp.nfa; ++ifa)
        fa_flag[ifa+nfa_old] = -1;
      // ----------------------------
      // cvofa: icv0 is always valid and contains the
      // global 0-based index of icv0. icv1 must be interpreted
      // based on the zone info -- see below.
      // ----------------------------
      resize(cvofa, nfa_old, nfa);
    }

    // ======================
    // cvs...
    // ======================
    if (ugp.ncv>0) {
      // ----------------------------
      // cv_flag...
      // ----------------------------
      resize(cv_flag, ncv_old, ncv);
      for (int icv = 0; icv<ugp.ncv; ++icv)
        cv_flag[icv+ncv_old] = ugp.cv_flag[icv];
    }

    // ======================
    // face zones...
    // ======================
    for (list<FaZone>::iterator zone = ugp.faZoneList.begin(); zone!=ugp.faZoneList.end(); ++zone) {
      // look for a match in the existing faZoneList...
      int zone_index = -1;
      FaZone * existing_zone = getFaZone(zone->getName());
      if (existing_zone!=NULL) {
        // a zone exists with the same name - make sure it is the same kind (e.g. INTERNAL, BOUNDARY, etc)...
        assert(existing_zone->getKind()==zone->getKind());
        // and has the same periodicity (if present)...
        FOR_I3
          assert(existing_zone->periodic_data[i]==zone->periodic_data[i]);
        zone_index = existing_zone->getIndex();
      }
      else {
        // adding a new zone automagically returns a new unique index...
        zone_index = addFaZone(&(*zone));
      }
      // now modify the fa_flag in the new range and set cvofa...
      assert(zone_index>=0);
      for (int ifa = 0; ifa<ugp.nfa; ++ifa) {
        if (ugp.fa_flag[ifa]==zone->getIndex()) {
          // check...
          assert(fa_flag[ifa+nfa_old]==-1);
          fa_flag[ifa+nfa_old] = zone_index;
          // icv0 is always a global cv number...
          int icv0_ugp = ugp.cvofa[ifa][0];
          assert((icv0_ugp>=0)&&(icv0_ugp<cvora_ugp[mpi_size]));
          int rank = 0;
          while (icv0_ugp>=cvora_ugp[rank+1])
            ++rank;
          // update its index...
          cvofa[ifa+nfa_old][0] = icv0_ugp+cvora_old[rank+1];
          // icv1 may be -1, in which case it remains unchanged. Another possibility is
          // that it is a global face indexing associated with a periodic pairing...
          int icv1_ugp = ugp.cvofa[ifa][1];
          if (icv1_ugp>=0) {
            if (zone->isPeriodic()) {
              // icv1_ugp is a global face index...
              assert((icv1_ugp>=0)&&(icv1_ugp<faora_ugp[mpi_size]));
              rank = 0;
              while (icv1_ugp>=faora_ugp[rank+1])
                ++rank;
              cvofa[ifa+nfa_old][1] = icv1_ugp+faora_old[rank+1];
            }
            else {
              // icv1_ugp is a global cv index...
              assert((icv1_ugp>=0)&&(icv1_ugp<cvora_ugp[mpi_size]));
              rank = 0;
              while (icv1_ugp>=cvora_ugp[rank+1])
                ++rank;
              cvofa[ifa+nfa_old][1] = icv1_ugp+cvora_old[rank+1];
            }
          }
          else {
            // should be -1...
            assert(icv1_ugp==-1);
            cvofa[ifa+nfa_old][1] = -1;
          }
        }
      }
    }
    // check...
    for (int ifa = 0; ifa<ugp.nfa; ++ifa)
      assert(fa_flag[ifa+nfa_old]>=0);

    /*
     for (list<CvZone>::iterator zone = ugp.cvZoneList.begin(); zone != ugp.cvZoneList.end(); ++zone) {
     cout << "got cv zone: " << zone->getName() << endl;
     */

    // should have a corresponding check for cvs...

    // cleanup...
    delete[] noora_old;
    delete[] noora_ugp;
    delete[] faora_old;
    delete[] faora_ugp;
    delete[] cvora_old;
    delete[] cvora_ugp;

    return *this;
  }

protected:

  // this is a bit ugly, and could be removed if faZoneList were
  // a vector instead...
  int * fa_kind_table;

  int n_periodic_bit_table;
  int * periodic_bit_table;

private:

  // these are global counts?
  int no_count;
  int fa_count;
  int cv_count;

  // io stuff...
  int io_common_flag;
  int io_read_flag;
  int io_write_flag;
  int byte_swap;
  MPI_File fh;
  MPI_Offset offset, int_size, double_size, header_size;
  int my_no_count, // local and global node counts
      my_fa_count, // local and global face counts
      my_noofa_count, noofa_count, // ...
      my_cv_count;
  int * my_no_list; // size my_no_count
  int * my_fa_list; // size my_fa_count
  int * my_cv_list; // size my_cv_count
  MPI_Datatype MPI_Header, no_i1_type, no_d1_type, no_d3_type, fa_i1_type, fa_i2_type, fa_d1_type, noofa_i1_type,
      cv_i1_type, cv_d1_type, cv_d3_type, cv_d33_type;

public:

  Ugp() {

    if (mpi_rank==0) cout<<"Ugp()"<<endl;

    ncv_b = ncv_bp = ncv = 0;

    nno = 0;
    nno_i = 0;

    nfa_b = nfa_bp = nfa_bpi = nfa = 0;

    noofa_i = NULL;
    noofa_v = NULL;
    noofa_s = 0;
    cvofa = NULL;

    faocv_i = NULL;
    faocv_s = 0;
    faocv_v = NULL;

    noocv_i = NULL;
    noocv_s = 0;
    noocv_v = NULL;

    cv_flag = NULL;
    fa_flag = NULL;
    no_flag = NULL;

    noora = NULL;
    faora = NULL;
    cvora = NULL;

    // partitioning info...
    npart = 0;
    cv_part = NULL;

    // node coordinates...
    x_no = NULL;

    // get rid of this eventually...
    fa_kind_table = NULL;

    // a table of periodic bit pairs: each unique
    // periodic transform is associated with a bit, stored in
    // the faZone->bits, and has an inverse transform, associated with
    // another face zone. This matching is not demanded of the user
    // when they specify the grid, but is instead worked out by the Ugp
    // redistReconnectReorder routine. This allows a developer to more
    // easily add and/or remove periodicity from a given mesh.

    n_periodic_bit_table = 0;
    periodic_bit_table = NULL;

    // io_status...
    io_common_flag = 0;
    io_read_flag = 0;
    io_write_flag = 0;
    my_no_list = NULL;
    my_fa_list = NULL;
    my_cv_list = NULL;

    // cvGroup - this will be used to modify the partition if specified...
    cvGroup = NULL;

  }

  ~Ugp() {

    if (mpi_rank==0) cout<<"~Ugp()"<<endl;

    if (noofa_i!=NULL) delete[] noofa_i;
    if (noofa_v!=NULL) delete[] noofa_v;
    noofa_s = 0;

    if (cvofa!=NULL) delete[] cvofa;

    if (faocv_i!=NULL) delete[] faocv_i;
    if (faocv_v!=NULL) delete[] faocv_v;
    faocv_s = 0;

    if (noocv_i!=NULL) delete[] noocv_i;
    if (noocv_v!=NULL) delete[] noocv_v;
    noocv_s = 0;

    if (cv_part!=NULL) delete[] cv_part;

    if (noora!=NULL) delete[] noora;
    if (faora!=NULL) delete[] faora;
    if (cvora!=NULL) delete[] cvora;

    if (cv_flag!=NULL) delete[] cv_flag;
    if (fa_flag!=NULL) delete[] fa_flag;
    if (no_flag!=NULL) delete[] no_flag;

    if (x_no!=NULL) delete[] x_no;

    if (fa_kind_table!=NULL) delete[] fa_kind_table;
    if (periodic_bit_table!=NULL) delete[] periodic_bit_table;

    if (my_no_list!=NULL) delete[] my_no_list;
    if (my_fa_list!=NULL) delete[] my_fa_list;
    if (my_cv_list!=NULL) delete[] my_cv_list;

  }

  // ===================================================
  // ===================================================
  // ===================================================
  // ===================================================

  void calcCvVolumeAndCentroid(double& vol_cv, double * x_cv, const int icv) {

    /// HACK XXXXX fix this eventually

    FOR_I3
      x_cv[i] = 0.0;
    const int noc_f = noocv_i[icv];
    const int noc_l = noocv_i[icv+1]-1;
    for (int noc = noc_f; noc<=noc_l; ++noc) {
      const int ino = noocv_v[noc];
      FOR_I3
        x_cv[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (noc_l-noc_f+1);
    FOR_I3
      x_cv[i] *= tmp;

    vol_cv = 1.0;

  }

  void calcCvCenter(double * x_cv, const int icv) {

    // returns an "approximate center" of the cv based on the
    // simple average of cv nodes...

    assert((icv>=0)&(icv<ncv));

    FOR_I3
      x_cv[i] = 0.0;
    const int noc_f = noocv_i[icv];
    const int noc_l = noocv_i[icv+1]-1;
    for (int noc = noc_f; noc<=noc_l; ++noc) {
      const int ino = noocv_v[noc];
      FOR_I3
        x_cv[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (noc_l-noc_f+1);
    FOR_I3
      x_cv[i] *= tmp;

  }

  void calcFaceCenter(double * x_fa, const int ifa) {

    // returns an "approximate center" of the face based on the
    // simple average of face nodes...

    assert((ifa>=0)&(ifa<nfa));

    FOR_I3
      x_fa[i] = 0.0;
    const int nof_f = noofa_i[ifa];
    const int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof<=nof_l; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3
        x_fa[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (nof_l-nof_f+1);
    FOR_I3
      x_fa[i] *= tmp;

  }

  void calcFaceCenterAndNormal(double * x_fa, double * n_fa, const int ifa) {

    // returns an "approximate center" of the face based on the
    // simple average of face nodes, and then the correct normal (the
    // normal with area magnitude does not depend on the definition of
    // the face center).

    assert((ifa>=0)&(ifa<nfa));

    FOR_I3
      x_fa[i] = 0.0;
    const int nof_f = noofa_i[ifa];
    const int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof<=nof_l; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3
        x_fa[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (nof_l-nof_f+1);
    FOR_I3
      x_fa[i] *= tmp;

    FOR_I3
      n_fa[i] = 0.0;
    int ino1 = noofa_v[nof_l];
    for (int nof = nof_f; nof<=nof_l; ++nof) {
      int ino0 = ino1;
      ino1 = noofa_v[nof];
      double dx0[3];
      FOR_I3
        dx0[i] = x_no[ino0][i]-x_fa[i];
      double dx1[3];
      FOR_I3
        dx1[i] = x_no[ino1][i]-x_fa[i];
      n_fa[0] += dx0[1]*dx1[2]-dx0[2]*dx1[1];
      n_fa[1] += dx0[2]*dx1[0]-dx0[0]*dx1[2];
      n_fa[2] += dx0[0]*dx1[1]-dx0[1]*dx1[0];
    }
    FOR_I3
      n_fa[i] *= 0.5;

  }

  void calcFaceCenterAndUnitNormal(double * x_fa, double * n_fa, const int ifa) {

    // returns an "approximate center" of the face based on the
    // simple average of face nodes, and then the correct normal (the
    // normal with area magnitude does not depend on the definition of
    // the face center).

    assert((ifa>=0)&(ifa<nfa));

    FOR_I3
      x_fa[i] = 0.0;
    const int nof_f = noofa_i[ifa];
    const int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof<=nof_l; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3
        x_fa[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (nof_l-nof_f+1);
    FOR_I3
      x_fa[i] *= tmp;

    FOR_I3
      n_fa[i] = 0.0;
    int ino1 = noofa_v[nof_l];
    for (int nof = nof_f; nof<=nof_l; ++nof) {
      int ino0 = ino1;
      ino1 = noofa_v[nof];
      double dx0[3];
      FOR_I3
        dx0[i] = x_no[ino0][i]-x_fa[i];
      double dx1[3];
      FOR_I3
        dx1[i] = x_no[ino1][i]-x_fa[i];
      n_fa[0] += dx0[1]*dx1[2]-dx0[2]*dx1[1];
      n_fa[1] += dx0[2]*dx1[0]-dx0[0]*dx1[2];
      n_fa[2] += dx0[0]*dx1[1]-dx0[1]*dx1[0];
    }
    double inv_mag = 1.0/sqrt(n_fa[0]*n_fa[0]+n_fa[1]*n_fa[1]+n_fa[2]*n_fa[2]);
    FOR_I3
      n_fa[i] *= inv_mag;
  }

  // these should apply to minimum initializations...
  // transformations...

  void rotateX(const double degrees) {

    if (mpi_rank==0) cout<<"Ugp::rotateX(): "<<degrees<<endl;

    double rad = degrees*M_PI/180.0;
    double cos_t = cos(rad);
    double sin_t = sin(rad);
    for (int ino = 0; ino<nno; ino++) {
      double y = x_no[ino][1]*cos_t-x_no[ino][2]*sin_t;
      double z = x_no[ino][2]*cos_t+x_no[ino][1]*sin_t;
      x_no[ino][1] = y;
      x_no[ino][2] = z;
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      int kind = zone->getKind();
      switch (kind) {
      case FA_ZONE_PERIODIC_CART:
        zone->rotatePeriodicDataX(degrees);
        break;
      case FA_ZONE_PERIODIC_CYL_Y:
      case FA_ZONE_PERIODIC_CYL_Z:
        if (mpi_rank==0) cerr<<"Error: cannot x-rotate a mesh with CYL_Y or CYL_Z periodicity"<<endl;
        throw(-1);
      }
    }

  }

  void rotateZ(const double degrees) {

    if (mpi_rank==0) cout<<"Ugp::rotateZ(): "<<degrees<<endl;

    double rad = degrees*M_PI/180.0;
    double cos_t = cos(rad);
    double sin_t = sin(rad);
    for (int ino = 0; ino<nno; ino++) {
      double x = x_no[ino][0]*cos_t-x_no[ino][1]*sin_t;
      double y = x_no[ino][1]*cos_t+x_no[ino][0]*sin_t;
      x_no[ino][0] = x;
      x_no[ino][1] = y;
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      int kind = zone->getKind();
      switch (kind) {
      case FA_ZONE_PERIODIC_CART:
        zone->rotatePeriodicDataZ(degrees);
        break;
      case FA_ZONE_PERIODIC_CYL_X:
      case FA_ZONE_PERIODIC_CYL_Y:
        if (mpi_rank==0) cerr<<"Error: cannot z-rotate a mesh with CYL_X or CYL_Y periodicity"<<endl;
        throw(-1);
      }
    }

  }

  void scale(const double factor) {

    if (mpi_rank==0) cout<<"Ugp::scale: "<<factor<<endl;

    for (int ino = 0; ino<nno; ino++) {
      x_no[ino][0] *= factor;
      x_no[ino][1] *= factor;
      x_no[ino][2] *= factor;
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      if (zone->getKind()==FA_ZONE_PERIODIC_CART) zone->scalePeriodicData(factor);
    }

  }

  void scale(const double * factor) {

    if (mpi_rank==0) cout<<"Ugp::scale: "<<factor[0]<<" "<<factor[1]<<" "<<factor[2]<<endl;

    for (int ino = 0; ino<nno; ino++) {
      x_no[ino][0] *= factor[0];
      x_no[ino][1] *= factor[1];
      x_no[ino][2] *= factor[2];
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      switch (zone->getKind()) {
      case FA_ZONE_PERIODIC_CART:
        zone->scalePeriodicData(factor);
        break;
      case FA_ZONE_PERIODIC_CYL_X:
        // aniotropic scaling is fine for CYL_X periodicity as long as y-z is isotropic...
        assert(factor[1]==factor[2]);
        break;
      case FA_ZONE_PERIODIC_CYL_Y:
        assert(factor[0]==factor[2]);
        break;
      case FA_ZONE_PERIODIC_CYL_Z:
        assert(factor[0]==factor[1]);
        break;
      }
    }

  }

  void shearXY(const double factor) {

    if (mpi_rank==0) cout<<"Ugp::shearXY: "<<factor<<endl;

    for (int ino = 0; ino<nno; ino++) {
      x_no[ino][0] += factor*x_no[ino][1];
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      switch (zone->getKind()) {
      case FA_ZONE_PERIODIC_CART:
        zone->shearXYPeriodicData(factor);
        break;
      case FA_ZONE_PERIODIC_CYL_X:
      case FA_ZONE_PERIODIC_CYL_Y:
      case FA_ZONE_PERIODIC_CYL_Z:
        cerr<<"Error :shear cannot be applied to grid with cylindrical periodicity"<<endl;
        throw(-1);
      }
    }

  }

  void translate(const double * dx) {

    if (mpi_rank==0) cout<<"Ugp::translate: "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;

    for (int ino = 0; ino<nno; ino++) {
      x_no[ino][0] += dx[0];
      x_no[ino][1] += dx[1];
      x_no[ino][2] += dx[2];
    }

    // and modify any periodic transforms...
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      switch (zone->getKind()) {
      case FA_ZONE_PERIODIC_CYL_X:
        // we cannot tranlate in the plane of the axis...
        assert((dx[1]==0.0)&&(dx[2]==0.0));
        break;
      case FA_ZONE_PERIODIC_CYL_Y:
        assert((dx[0]==0.0)&&(dx[2]==0.0));
        break;
      case FA_ZONE_PERIODIC_CYL_Z:
        assert((dx[0]==0.0)&&(dx[1]==0.0));
        break;
      }
    }

  }

  void setCvGroup() {

    // ================================================
    // building of unique cv groups is used to introduce an
    // additional constraint on the cv partition (to handle linelets
    // or other preconditioning for example). The task of producing
    // the group is simplified for the user to the task of specifying the
    // vistual function cvofaGroupHook to return an int of 1 if the
    // nbring cvs should be in the same group, and zero if not. The value
    // for boundary faces can be either 1 or zero, and will not matter.
    // ================================================

    // it is possible that the group is NULL, i.e. it has not been registered. If this
    // is the case, then register it here...
    if (cvGroup==NULL) registerData(cvGroup, "CV_GROUP", CV_DATA);

    // clear the cvGroup int...
    FOR_ICV
      cvGroup[icv] = -1;

    // set the fa_flag to indicate the barriers between groups: 1 = implicit, 0 = explicit
    setFaFlagForCvGroupHook();

    // check...
    FOR_IFA
      assert((fa_flag[ifa]==0)||(fa_flag[ifa]==1));

    // to ensure consistency, reduce using MAX_DATA just in case the user's
    // cvofaGroupHook was not symmetric - this errs on the side of grouping...
    updateFaData(fa_flag, MAX_DATA);

    // now group...


    MPI_Pause("about to group");

  }

  virtual void setFaFlagForCvGroupHook() {
    // set the fa_flag to 1 if the nbring cvs should
    // be put in the same group, otherwise zero...
    FOR_IFA
      fa_flag[ifa] = 0;
  }

  // ===================================================
  // ===================================================
  // ===================================================
  // ===================================================

  FaZone * getFaZone(const string& name) {
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++) {
      if (zone->getName()==name) {
        return (&(*zone));
      }
    }
    return (NULL);
  }

  int getCvOfFa(const int ifa, const int i) const {
#ifdef DEBUG
    assert( (ifa >= 0)&&(ifa < nfa) );
    assert( (i == 0)||(i == 1) );
#endif
    return (cvofa[ifa][i]);
  }

  int getFaFlag(const int ifa) const {
#ifdef DEBUG
    assert( (ifa >= 0)&&(ifa < nfa) );
#endif
    return (fa_flag[ifa]);
  }

  void add(Ugp& ugp) {

    cout<<"wow"<<endl;

  }

  void preinit();

  void splitNonPlanarFaces(const double tol);

  // ===================================================
  // ===================================================
  // ===================================================
  // ===================================================

  void init() {

    if (mpi_rank==0) cout<<" > Ugp::init()"<<endl;

    // ----------------------------------------------------
    // we need valid versions of noora,faora,cvora...
    // these should have been constructed by preinit...
    // ----------------------------------------------------
    assert(noora!=NULL);
    assert(faora!=NULL);
    assert(cvora!=NULL);

    // in some cases we need to repart...
    if (cvGroup!=NULL) setCvGroupPartMetis(); // uses Metis/ParMetis to partition cvs based on cvGroup[icv]
    else setCvPartMetis(); // uses Metis/ParMetis to partition cvs

    // setCvPartCheap();           // puts the first ncv/np cells on rank0, next, etc...
    // setCvPartBad();             // salt and pepper - cvs alternating

    // and check...
    checkCvPart();

    // we always need to redistribute the grid based on the part...
    redistReconnectReorder();

    // add a transformation at this point if desired...
    transformMeshHook();

    // Reconstruct cvora - because the partition is cv-based, this
    // still has meaning. Note in the MPI_Allgather call below,
    // cvora+1 is a simple way of saying &(cvora[1])... i.e. we put the result
    // into the buffer starting at the second position. This simplifies the
    // subsequent build of cvora...

    MPI_Allgather(&ncv, 1, MPI_INT, cvora+1, 1, MPI_INT, mpi_comm);
    cvora[0] = 0;
    for (int rank = 0; rank<mpi_size; ++rank)
      cvora[rank+1] += cvora[rank];
    assert(cvora[mpi_rank+1]-cvora[mpi_rank]==ncv);

    // this check fails for certain filters, so skip. If you need to know the
    // global cv count, use cvora[mpi_size]...
    /*
     cout << "cv_count: " << cv_count << endl;
     assert( cvora[mpi_size] == cv_count );
     */

  }

  void postinit() {

    if (mpi_rank==0) cout<<" > Ugp::postinit()"<<endl;

    // call this to sync nodes when new periodic connections have been created...
    syncNodes();
    checkFaceCorrespondence();
    checkNodeCorrespondence();

  }

  int getInversePeriodicBits(const int bits) {

    // for the case of no bits set, simply return 0...
    if (bits==0) return (0);

    // to use this function, the periodic_bit_table must exist...
    assert(n_periodic_bit_table>0);
    assert(periodic_bit_table!=NULL);

    int inverse_bits = 0;
    for (int bit = 0; bit<n_periodic_bit_table; bit++) {
      if (bits&(1<<bit)) {
        inverse_bits |= (1<<(periodic_bit_table[bit]));
      }
    }

    return (inverse_bits);

  }

  int addPeriodicBits(const int b1, const int b2) {

    // adding periodic bits is a bit complex because if
    // one has a bits set and the other has an inverse bit set, then
    // the result should have neither set. Also, we should never be
    // adding the same bits (implication: 2X the transform)...

    if (b1==0) return (b2);
    else if (b2==0) return (b1);

    // to use this function, the periodic_bit_table must exist...
    assert(n_periodic_bit_table>0);
    assert(periodic_bit_table!=NULL);

    int add_bits = 0;
    for (int bit = 0; bit<n_periodic_bit_table; bit++) {
      assert(bit!=periodic_bit_table[bit]);
      int b1_set = (b1&(1<<bit));
      int b2_set = (b2&(1<<bit));
      int b1i_set = (b1&(1<<(periodic_bit_table[bit])));
      int b2i_set = (b2&(1<<(periodic_bit_table[bit])));
      assert(!((b1_set)&&(b2_set)));
      assert(!((b1i_set)&&(b2i_set)));
      assert(!((b1_set)&&(b1i_set)));
      assert(!((b2_set)&&(b2i_set)));
      if (((b1_set)&&(!b2i_set))||((b2_set)&&(!b1i_set))) {
        add_bits |= (1<<bit);
      }
    }

    return (add_bits);

  }

  void dumpPeriodicBits(const int bits) {

    if (bits==0) {
      cout<<"0"<<endl;
      return;
    }

    assert(n_periodic_bit_table>0);
    assert(periodic_bit_table!=NULL);

    for (int bit = 0; bit<n_periodic_bit_table; bit++) {
      assert(bit!=periodic_bit_table[bit]);
      if (bits&(1<<bit)) {
        cout<<"1";
      }
      else {
        cout<<"0";
      }
    }
    cout<<endl;

  }

  void dumpPeriodicBitsTable() {

    assert(n_periodic_bit_table>0);
    assert(periodic_bit_table!=NULL);

    for (int bit = 0; bit<n_periodic_bit_table; bit++)
      cout<<" bit: "<<bit<<" periodic_bit_table[bit]: "<<periodic_bit_table[bit]<<endl;
  }

  // ==========================================================================
  // data-free registration routines: used to register data based on names
  // when we don't explicitly have a pointer available
  // ==========================================================================

  void registerData(ParamMap& params) {

    Param * p;
    if (params.getParam(p, "REGISTER_INT_VALUE")||params.getParam(p, "REGISTER_I0")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerIntValue(p->getString(i));
      }
    }
    if (params.getParam(p, "REGISTER_DOUBLE_VALUE")||params.getParam(p, "REGISTER_R0")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerDoubleValue(p->getString(i));
      }
    }
    if (params.getParam(p, "REGISTER_CV_DOUBLE_SCALAR")||params.getParam(p, "REGISTER_CV_R1")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerDoubleScalar(p->getString(i), CV_DATA);
      }
    }
    if (params.getParam(p, "REGISTER_NO_DOUBLE_SCALAR")||params.getParam(p, "REGISTER_NO_R1")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerDoubleScalar(p->getString(i), NO_DATA);
      }
    }
    if (params.getParam(p, "REGISTER_CV_DOUBLE_VECTOR")||params.getParam(p, "REGISTER_CV_R2")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerDoubleVector(p->getString(i), CV_DATA);
      }
    }
    if (params.getParam(p, "REGISTER_NO_DOUBLE_VECTOR")||params.getParam(p, "REGISTER_NO_R2")) {
      for (int i = 1; i<p->getSize(); i++) {
        registerDoubleVector(p->getString(i), NO_DATA);
      }
    }

  }

  void registerIntValue(const string& name) {
    registerIntValue(name.c_str());
  }

  void registerIntValue(const char * name) {
    int value = 0;
    intList.push_back(value);
    registerData(intList.back(), name);
  }

  void registerDoubleValue(const string& name) {
    registerDoubleValue(name.c_str());
  }

  void registerDoubleValue(const char * name) {
    double value = 0.0;
    doubleList.push_back(value);
    registerData(doubleList.back(), name);
  }

  void registerDoubleScalar(const string& name, const int datatype) {
    registerDoubleScalar(name.c_str(), datatype);
  }

  void registerDoubleScalar(const char * name, const int datatype) {
    double * ptr = NULL;
    doublePtrList.push_back(ptr);
    registerData(doublePtrList.back(), name, datatype);
  }

  void registerDoubleVector(const string& name, const int datatype) {
    registerDoubleVector(name.c_str(), datatype);
  }

  void registerDoubleVector(const char * name, const int datatype) {
    double (*ptr)[3] = NULL;
    doublePtr3List.push_back(ptr);
    registerData(doublePtr3List.back(), name, datatype);
  }

  // ==========================================================================
  // these are provided for backward compatibility. Use the standard registerData
  // with the type of registration inferred from the arguments. These will be
  // phased out evenytually.
  // ==========================================================================
  void registerValue(int &value, const char * name) {
    registerData(value, name);
  }

  void registerValue(int &value, const string &name) {
    registerData(value, name.c_str());
  }

  void registerValue(double &value, const char * name) {
    registerData(value, name);
  }

  void registerValue(double &value, const string &name) {
    registerData(value, name.c_str());
  }

  void registerScalar(double *&scalar, const char * name, int datatype) {
    registerData(scalar, name, datatype);
  }

  void registerScalar(double *&scalar, const string &name, int datatype) {
    registerData(scalar, name.c_str(), datatype);
  }

  void registerVector(double(*&vector)[3], const char * name, int datatype) {
    registerData(vector, name, datatype);
  }

  void registerVector(double(*&vector)[3], const string &name, int datatype) {
    registerData(vector, name.c_str(), datatype);
  }

  void registerTensor(double(*&tensor)[3][3], const char * name, int datatype) {
    registerData(tensor, name, datatype);
  }

  void registerTensor(double(*&tensor)[3][3], const string &name, int datatype) {
    registerData(tensor, name.c_str(), datatype);
  }

  // ================================================
  // registration...
  // TODO: XXXXX still need to provide common access functions
  // for registered data. Unfortunately, because the return value
  // of an argument cannot be used to determine which interface to
  // use, this limits our flexibility a bit.
  // ================================================

  void registerData(int &value, const char * name) {

    if (mpi_rank==0) cout<<"registerValue (int)   : "<<name<<endl;

    // check that the name does not conflict with any other registered data...
    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: value already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // push it into list...

    intValueList.push_back(IntValue(value, name));

    // and zero?...

    value = 0;
  }

  IntValue * getIntValueData(const string& name) {
    return (getIntValueData(name.c_str()));
  }

  IntValue * getIntValueData(const char * name) {

    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  // ########################################################

  void registerData(double &value, const char * name) {

    if (mpi_rank==0) cout<<"registerValue (double): "<<name<<endl;

    // check that the name does not conflict with any other registered data...
    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: value already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // push it into list...

    doubleValueList.push_back(DoubleValue(value, name));

    // and zero...

    value = 0.0;

  }

  DoubleValue * getDoubleValueData(const string& name) {
    return (getDoubleValueData(name.c_str()));
  }

  DoubleValue * getDoubleValueData(const char * name) {

    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  DoubleScalar * getDoubleScalarData(const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  DoubleVector * getDoubleVectorData(const char * name) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  void unregisterScalar(double *&scalar) {
    // no-op for now
  }

  void registerData(int *&scalar, const char * name, int datatype) {

    // check that the name does not conflict with any other registered data...
    for (list<IntScalar>::iterator data = intScalarList.begin(); data!=intScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: scalar already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // user is required to NULLIFY so they know that they are not to manage memory...

    assert(scalar==NULL);

    // push it into list...

    intScalarList.push_back(IntScalar(scalar, name, datatype));

    // you can register when the data size is still unknown...

    int n = 0;
    switch (datatype) {
    case CV_DATA:
      if (mpi_rank==0) cout<<"registerData int CV_DATA: "<<name<<endl;
      n = ncv;
      break;
    case FA_DATA:
      if (mpi_rank==0) cout<<"registerData int FA_DATA: "<<name<<endl;
      n = nfa;
      break;
    case NO_DATA:
      if (mpi_rank==0) cout<<"registerData int NO_DATA: "<<name<<endl;
      n = nno;
      break;
    default:
      if (mpi_rank==0) cout<<"Error: unsupported datatype: "<<datatype<<endl;
      throw(-1);
    }

    if (n>0) {
      scalar = new int[n];
      for (int i = 0; i<n; i++) {
        scalar[i] = 0;
      }
    }

  }

  void registerData(double *&scalar, const char * name, int datatype) {

    // check that the name does not conflict with any other registered data...
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: scalar already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // user is required to NULLIFY so they know that they are not to manage memory...

    assert(scalar==NULL);

    // push it into list...

    doubleScalarList.push_back(DoubleScalar(scalar, name, datatype));

    // you can register when the data size is still unknown...

    int n = 0;
    switch (datatype) {
    case CV_DATA:
      if (mpi_rank==0) cout<<"registerScalar CV_DATA: "<<name<<endl;
      n = ncv;
      break;
    case FA_DATA:
      if (mpi_rank==0) cout<<"registerScalar FA_DATA: "<<name<<endl;
      n = nfa;
      break;
    case NO_DATA:
      if (mpi_rank==0) cout<<"registerScalar NO_DATA: "<<name<<endl;
      n = nno;
      break;
    default:
      if (mpi_rank==0) cout<<"Error: unsupported datatype: "<<datatype<<endl;
      throw(-1);
    }

    if (n>0) {
      scalar = new double[n];
      for (int i = 0; i<n; i++) {
        scalar[i] = 0.0;
      }
    }

  }

  void registerData(double(*&vector)[3], const char * name, int datatype) {

    // check that the name does not conflict with any other registered data...

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: vector already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // user is required to NULLIFY so they know that they are not to manage memory...

    assert(vector==NULL);

    // push it into list...

    doubleVectorList.push_back(DoubleVector(vector, name, datatype));

    // you can register when the data size is still unknown...

    int n = 0;
    switch (datatype) {
    case CV_DATA:
      if (mpi_rank==0) cout<<"registerVector CV_DATA: "<<name<<endl;
      n = ncv;
      break;
    case FA_DATA:
      if (mpi_rank==0) cout<<"registerVector FA_DATA: "<<name<<endl;
      n = nfa;
      break;
    case NO_DATA:
      if (mpi_rank==0) cout<<"registerVector NO_DATA: "<<name<<endl;
      n = nno;
      break;
    default:
      if (mpi_rank==0) cout<<"Error: unsupported datatype: "<<datatype<<endl;
      throw(-1);
    }

    if (n>0) {
      vector = new double[n][3];
      for (int i = 0; i<n; i++) {
        for (int j = 0; j<3; j++) {
          vector[i][j] = 0.0;
        }
      }
    }

  }

  void registerData(double(*&tensor)[3][3], const char * name, int datatype) {

    // check that the name does not conflict with any other registered data...

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        if (mpi_rank==0) cerr<<"Error: tensor already registered with name: "<<name<<endl;
        throw(-1);
      }
    }

    // user is required to NULLIFY so they know that they are not to manage memory...

    assert(tensor==NULL);

    // push it into list...

    doubleTensorList.push_back(DoubleTensor(tensor, name, datatype));

    // you can register when the data size is still unknown...

    int n = 0;
    switch (datatype) {
    case CV_DATA:
      if (mpi_rank==0) cout<<"registerTensor CV_DATA: "<<name<<endl;
      n = ncv;
      break;
    case FA_DATA:
      if (mpi_rank==0) cout<<"registerTensor FA_DATA: "<<name<<endl;
      n = nfa;
      break;
    case NO_DATA:
      if (mpi_rank==0) cout<<"registerTensor NO_DATA: "<<name<<endl;
      n = nno;
      break;
    default:
      if (mpi_rank==0) cout<<"Error: unsupported datatype: "<<datatype<<endl;
      throw(-1);
    }

    if (n>0) {
      tensor = new double[n][3][3];
      for (int i = 0; i<n; i++) {
        for (int j = 0; j<3; j++) {
          for (int k = 0; k<3; k++) {
            tensor[i][j][k] = 0.0;
          }
        }
      }
    }

  }

  int getI0(const char * name) {

    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find registered I0: "<<name<<endl;
    throw(-1);

  }

  DoubleScalar * getScalarData(const string& name) {
    return (getScalarData(name.c_str()));
  }

  DoubleScalar * getScalarData(const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  double * getScalar(const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find registered scalar: "<<name<<endl;
    throw(-1);
  }

  double * getR1(const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find registered R1: "<<name<<endl;
    throw(-1);
  }

  double * getRegisteredDoubleScalar(const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // throw an error because a NULL could be ambiguous here...
    cerr<<"Error: could not find registered double scalar with name: "<<name<<endl;
    throw(-1);

  }

  DoubleVector * getVectorData(const string& name) {
    return (getVectorData(name.c_str()));
  }

  DoubleVector * getVectorData(const char * name) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  double (* getR2(const char * name))[3] {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find registered R2: "<<name<<endl;
    throw(-1);

  }

  double (* getRegisteredDoubleVector(const char * name))[3] {
    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // throw an error because a NULL could be ambiguous here...
    cerr<<"Error: could not find registered double vector with name: "<<name<<endl;
    throw(-1);

  }

  DoubleTensor * getTensorData(const string& name) {
    return (getTensorData(name.c_str()));
  }

  DoubleTensor * getTensorData(const char * name) {

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return (&(*data));
    }
    return (NULL);
  }

  double (* getR3(const char * name))[3][3]
  {

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (strcmp(name, data->getName())==0) return ((*data->ptr));
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find registered R3: "<<name<<endl;
    throw(-1);

  }

  void setR1Name(const double * ptr, const char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if ((*data->ptr)==ptr) {
        if (mpi_rank==0) cout<<" > changing R1 name "<<data->getName()<<" to "<<name<<endl;
        data->setName(name);
        return;
      }
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find R1 in setR1Name."<<endl;
    throw(-1);

  }

  void setR2Name(const double(*ptr)[3], const char * name) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if ((*data->ptr)==ptr) {
        if (mpi_rank==0) cout<<" > changing R2 name "<<data->getName()<<" to "<<name<<endl;
        data->setName(name);
        return;
      }
    }

    // Error - a NULL could be ambiguous here...
    cerr<<"Error: could not find R2 in setR2Name."<<endl;
    throw(-1);

  }

              /*
               double * getVector(char * name) {

               for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++) {
               if ( strcmp( name, data->getName() ) == 0 )
               return( (*data->ptr) );
               }
               // this routine has to error, because a NULL is ambiguous - i.e. could be
               // no registered data, or a zero-sized pointer...
               cerr << "Error: could not find registered vector: " << name << endl;
               throw(-1);
               }
               */

  void clearDataFlags() {

    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++)
      data->clearFlag();

    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++)
      data->clearFlag();

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
      data->clearFlag();

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
      data->clearFlag();

  }

  void setDataFlags() {

    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++)
      data->setFlag();

    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++)
      data->setFlag();

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
      data->setFlag();

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
      data->setFlag();

  }

  void setDataFlag(const char * name) {

    for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }

    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }

    if (mpi_rank==0) cout<<"Warning: could not set data flag for name: "<<name<<endl;

  }

  void setScalarFlag(char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cout<<"Warning: setScalarFlag: could not find: "<<name<<", skipping."<<endl;

  }

  void setDataFlag(double &value) {

    for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++) {
      if (data->ptr==&value) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cerr<<"Error: setDataFlag: could not find matching pointer."<<endl;
    throw(-1);

  }

  void setDataFlag(double *&scalar) {

    setScalarFlag(scalar);

  }

  void setScalarFlag(double *&scalar) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (data->ptr==&scalar) {
        data->setFlag();
        return;
      }
    }
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if ((*data->ptr!=NULL)&&(*data->ptr==scalar)) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cerr<<"Error: setScalarFlag: could not find matching pointer."<<endl;
    throw(-1);

  }

  int checkDataFlag(double *&scalar) {

    return checkScalarFlag(scalar);

  }

  int checkScalarFlag(char * name) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        return data->checkFlag();
      }
    }
    if (mpi_rank==0) cout<<"Warning: setScalarFlag: could not find: "<<name<<", skipping."<<endl;
    throw(-1);
  }

  int checkScalarFlag(double *&scalar) {

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
      if (data->ptr==&scalar) {
        return data->checkFlag();
      }
    }
    if (mpi_rank==0) cerr<<"Error: checkScalarFlag: could not find matching pointer."<<endl;
    throw(-1);

  }

  void setVectorFlag(char * name) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cout<<"Warning: setVectorFlag: could not find: "<<name<<", skipping."<<endl;

  }

  void setDataFlag(double(*&vector)[3]) {

    setVectorFlag(vector);

  }

  void setVectorFlag(double(*&vector)[3]) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (data->ptr==&vector) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cout<<"Warning: setVectorFlag: could not find matching pointer, skipping."<<endl;

  }

  int checkDataFlag(double(*&vector)[3]) {

    return checkVectorFlag(vector);

  }

  int checkVectorFlag(char * name) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        return data->checkFlag();
      }
    }
    if (mpi_rank==0) cout<<"Warning: setVectorFlag: could not find: "<<name<<", skipping."<<endl;
    throw(-1);

  }

  int checkVectorFlag(double(*&vector)[3]) {

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
      if (data->ptr==&vector) {
        return (data->checkFlag());
      }
    }
    if (mpi_rank==0) cerr<<"Error: checkVectorFlag: could not find matching pointer."<<endl;
    throw(-1);

  }

  void setTensorFlag(char * name) {

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (strcmp(name, data->getName())==0) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cout<<"Warning: doubleTensorList: could not find: "<<name<<", skipping."<<endl;

  }

  void setDataFlag(double(*&tensor)[3][3]) {

    setTensorFlag(tensor);

  }

  void setTensorFlag(double(*&tensor)[3][3]) {

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (data->ptr==&tensor) {
        data->setFlag();
        return;
      }
    }
    if (mpi_rank==0) cout<<"Warning: setTensorFlag: could not find matching pointer, skipping."<<endl;

  }

  int checkDataFlag(double(*&tensor)[3][3]) {

    return checkTensorFlag(tensor);

  }

  int checkTensorFlag(double(*&tensor)[3][3]) {

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
      if (data->ptr==&tensor) {
        return data->checkFlag();
      }
    }
    if (mpi_rank==0) cerr<<"Error: checkTensorFlag: could not find matching pointer."<<endl;
    throw(-1);

  }

  void allocateRegisteredNoData() {

    assert(nno>=0);

    if (nno>0) {

      // double scalars: [nno]...
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
        if (data->getDatatype()==NO_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nno];
          for (int ino = 0; ino<nno; ino++) {
            (*(data->ptr))[ino] = 0.0;
          }
        }
      }

      // double vectors: [nno][3]...
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
        if (data->getDatatype()==NO_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nno][3];
          for (int ino = 0; ino<nno; ino++) {
            for (int j = 0; j<3; j++) {
              (*(data->ptr))[ino][j] = 0.0;
            }
          }
        }
      }

      // double tensors: [nno][3][3]...
      for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
        if (data->getDatatype()==NO_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nno][3][3];
          for (int ino = 0; ino<nno; ino++) {
            for (int j = 0; j<3; j++) {
              for (int k = 0; k<3; k++) {
                (*(data->ptr))[ino][j][k] = 0.0;
              }
            }
          }
        }
      }

    }

  }

  void allocateRegisteredCvData() {

    assert(ncv>=0);

    if (ncv>0) {

      // double scalars: [ncv]...
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
        if (data->getDatatype()==CV_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[ncv];
          for (int icv = 0; icv<ncv; icv++) {
            (*(data->ptr))[icv] = 0.0;
          }
        }
      }

      // double vectors: [ncv][3]...
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
        if (data->getDatatype()==CV_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[ncv][3];
          for (int icv = 0; icv<ncv; icv++) {
            for (int j = 0; j<3; j++) {
              (*(data->ptr))[icv][j] = 0.0;
            }
          }
        }
      }

      // double tensors: [ncv][3][3]...
      for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
        if (data->getDatatype()==CV_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[ncv][3][3];
          for (int icv = 0; icv<ncv; icv++) {
            for (int j = 0; j<3; j++) {
              for (int k = 0; k<3; k++) {
                (*(data->ptr))[icv][j][k] = 0.0;
              }
            }
          }
        }
      }

    }

  }

  void allocateRegisteredFaData() {

    assert(nfa>=0);

    if (nfa>0) {

      // double scalars: [nfa]...
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
        if (data->getDatatype()==FA_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nfa];
          for (int ifa = 0; ifa<nfa; ifa++) {
            (*(data->ptr))[ifa] = 0.0;
          }
        }
      }

      // double vectors: [nfa][3]...
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
        if (data->getDatatype()==FA_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nfa][3];
          for (int ifa = 0; ifa<nfa; ifa++) {
            for (int j = 0; j<3; j++) {
              (*(data->ptr))[ifa][j] = 0.0;
            }
          }
        }
      }

      // double tensors: [nfa][3][3]...
      for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
        if (data->getDatatype()==FA_DATA) {
          assert(*(data->ptr)==NULL);
          *(data->ptr) = new double[nfa][3][3];
          for (int ifa = 0; ifa<nfa; ifa++) {
            for (int j = 0; j<3; j++) {
              for (int k = 0; k<3; k++) {
                (*(data->ptr))[ifa][j][k] = 0.0;
              }
            }
          }
        }
      }

    }

  }

private:

  int addFaZone(char * name, const int kind) {

    // cannot be periodic...
    assert((kind<FA_ZONE_PERIODIC_FIRST)||(kind>FA_ZONE_PERIODIC_LAST));

    double periodic_data[3];
    for (int i = 0; i<3; i++)
      periodic_data[i] = 0.0;

    // returns the face zone index...
    return addFaZone(name, kind, periodic_data);

  }

  int addFaZone(FaZone * zone) {

    return addFaZone(zone->getName(), zone->getKind(), zone->periodic_data);

  }

  int addFaZone(char * name, const int kind, double * periodic_data) {

    // must be valid...
    assert((kind>=FA_ZONE_FIRST)&&(kind<=FA_ZONE_LAST));

    // first make sure that a zone of this name does not already exist...
    int new_index = -1;
    list<FaZone>::iterator iz;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++) {
      if (strcmp(iz->getName(), name)==0) {
        // there is an existing zone with this name - make sure its characteristics
        // are the same...
        cout<<"XXXX found zone - need to check if same"<<endl;
        throw(-1);
      }
      new_index = max(iz->getIndex(), new_index);
    }

    faZoneList.push_back(FaZone());
    iz = faZoneList.end();
    iz--;

    iz->setName(name);
    iz->setKind(kind);

    new_index++;
    iz->setIndex(new_index);

    iz->setPeriodicData(periodic_data);

    if (mpi_rank==0) cout<<" > adding FaZone: "<<iz->getName()<<endl;

    // returns the face zone index...
    return (new_index);

  }

protected:

  // we need versions of these virtual functions to support translations/rotations
  // across periodic boundaries...

  void rangeTranslateVector(double * xi, const int n, const int flag, double * transform) {
    // for Ugp, the range flag stores the KIND of the associated zone. In the
    // future it may be required to combine radial and cartesian periodicity, in which
    // case we will need to modify this slightly, likely with bits or something, i.e.
    // allow cart + cyl_x, y or z.
    switch (flag) {
    case FA_ZONE_INTERNAL:
      return;
    case FA_ZONE_PERIODIC_CART:
      // Cartesian translation...
      for (int i = 0; i<n; i++) {
        for (int j = 0; j<3; j++) {
          xi[3*i+j] += transform[j];
        }
      }
      return;
    case FA_ZONE_PERIODIC_CYL_X:
      // rotate about x...
      for (int i = 0; i<n; i++) {
        double this_y = xi[3*i+1];
        double this_z = xi[3*i+2];
        xi[3*i+1] = transform[0]*this_y-transform[1]*this_z;
        xi[3*i+2] = transform[0]*this_z+transform[1]*this_y;
      }
      return;
    case FA_ZONE_PERIODIC_CYL_Z:
      // rotate about z...
      for (int i = 0; i<n; i++) {
        double this_x = xi[3*i+0];
        double this_y = xi[3*i+1];
        xi[3*i+0] = transform[0]*this_x-transform[1]*this_y;
        xi[3*i+1] = transform[0]*this_y+transform[1]*this_x;
        // and include any additional translation provided in transform[2]...
        // to support the case of multiple periodicity in annular sectors
        xi[3*i+2] += transform[2];
      }
      return;
    case FA_ZONE_PERIODIC_CYL_Y:
    default:
      cerr<<"Error: unsupported range flag: "<<flag<<endl;
      throw(-1);
    }
  }

  void rangeRotateVector(double * vi, const int n, const int flag, double * transform) {
    switch (flag) {
    case FA_ZONE_INTERNAL:
    case FA_ZONE_PERIODIC_CART:
      return;
    case FA_ZONE_PERIODIC_CYL_X:
      // rotate about x...
      for (int i = 0; i<n; i++) {
        double this_vy = vi[3*i+1];
        double this_vz = vi[3*i+2];
        vi[3*i+1] = transform[0]*this_vy-transform[1]*this_vz;
        vi[3*i+2] = transform[0]*this_vz+transform[1]*this_vy;
      }
      return;
    case FA_ZONE_PERIODIC_CYL_Z:
      // rotate about z...
      for (int i = 0; i<n; i++) {
        double this_vx = vi[3*i+0];
        double this_vy = vi[3*i+1];
        vi[3*i+0] = transform[0]*this_vx-transform[1]*this_vy;
        vi[3*i+1] = transform[0]*this_vy+transform[1]*this_vx;
      }
      return;
    case FA_ZONE_PERIODIC_CYL_Y:
    default:
      cerr<<"Error: unsupported range flag: "<<flag<<endl;
      throw(-1);
    }
  }

  void rangeRotateTensor(double * tij, const int n, const int flag, double * transform) {

    // in this case, the tensor is packed as t11, t12, t13, t21, t22, t23, t31, t32, t33, next...

    switch (flag) {
    case FA_ZONE_INTERNAL:
    case FA_ZONE_PERIODIC_CART:
      // no transformation neccessary...
      return;
    case FA_ZONE_PERIODIC_CYL_X: {

      // sector periodicity about x...
      double ct = transform[0]; // cos(theta)
      double st = transform[1]; // sin(theta)
      double ct2 = ct*ct;
      double st2 = st*st;
      double ctst = ct*st;

      // the following is the expansion of [R][T][Rt]
      // where [R] = [[1,  0,  0]
      //              [0, ct, -st]
      //              [0, st,  ct]]

      for (int i = 0; i<n; i++) {

        //double t11 = tij[9*i  ]; // unchanged
        double t12 = tij[9*i+1];
        double t13 = tij[9*i+2];

        double t21 = tij[9*i+3];
        double t22 = tij[9*i+4];
        double t23 = tij[9*i+5];

        double t31 = tij[9*i+6];
        double t32 = tij[9*i+7];
        double t33 = tij[9*i+8];

        //tij[9*i  ] = t11;
        tij[9*i+1] = ct*t12-st*t13;
        tij[9*i+2] = ct*t13+st*t12;

        tij[9*i+3] = ct*t21-st*t31;
        tij[9*i+4] = ct2*t22+st2*t33-ctst*(t23+t32);
        tij[9*i+5] = ct2*t23-st2*t32+ctst*(t22-t33);

        tij[9*i+6] = ct*t31+st*t21;
        tij[9*i+7] = ct2*t32-st2*t23+ctst*(t22-t33);
        tij[9*i+8] = ct2*t33+st2*t22+ctst*(t23+t32);

      }
    }
      return;
    case FA_ZONE_PERIODIC_CYL_Z: {

      // sector periodicity about z...
      double ct = transform[0]; // cos(theta)
      double st = transform[1]; // sin(theta)
      double ct2 = ct*ct;
      double st2 = st*st;
      double ctst = ct*st;

      // the following is the expansion of [R][T][Rt]
      // where [R] = [[ct, -st,  0]
      //              [st,  ct,  0]
      //              [ 0,   0,  1]]

      for (int i = 0; i<n; i++) {

        double t11 = tij[9*i];
        double t12 = tij[9*i+1];
        double t13 = tij[9*i+2];

        double t21 = tij[9*i+3];
        double t22 = tij[9*i+4];
        double t23 = tij[9*i+5];

        double t31 = tij[9*i+6];
        double t32 = tij[9*i+7];
        //double t33 = tij[9*i+8]; // unchanged

        tij[9*i] = ct2*t11-ctst*(t12+t21)+st2*t22;
        tij[9*i+1] = ctst*(t11-t22)+ct2*t12-st2*t21;
        tij[9*i+2] = ct*t13-st*t23;

        tij[9*i+3] = ctst*(t11-t22)-st2*t12+ct2*t21;
        tij[9*i+4] = st2*t11+ctst*(t12+t21)+ct2*t22;
        tij[9*i+5] = ct*t23+st*t13;

        tij[9*i+6] = ct*t31-st*t32;
        tij[9*i+7] = ct*t32+st*t31;
        //tij[9*i+8] = t33;

      }
    }
      return;
    default:
      cerr<<"Error: unsupported range flag: "<<flag<<endl;
      throw(-1);
    }
  }

  void rangeRotateSymmetricTensor(double * tij, const int n, const int flag, double * transform) {

    // To save space, symmetric tensors can be stored as six components.
    // in this case, the tensor is packed as t11, t22, t33, t23, t13, t12, next...
    // i.e. diagonal part first, then off-diagonal. Note that the ordering for the
    // off-diagonal part is based on cycling the missing index.

    switch (flag) {
    case FA_ZONE_INTERNAL:
    case FA_ZONE_PERIODIC_CART:
      // no transformation neccessary...
      return;

    case FA_ZONE_PERIODIC_CYL_X: {

      // sector periodicity about x...
      double ct = transform[0]; // cos(theta)
      double st = transform[1]; // sin(theta)
      double ct2 = ct*ct;
      double st2 = st*st;
      double ctst = ct*st;

      // the following is the expansion of [R][T][Rt]
      // where [R] = [[1,  0,  0]
      //              [0, ct, -st]
      //              [0, st,  ct]]

      for (int i = 0; i<n; i++) {

        //double t11 = tij[6*i  ];
        double t22 = tij[6*i+1];
        double t33 = tij[6*i+2];

        double t23 = tij[6*i+3];
        double t13 = tij[6*i+4];
        double t12 = tij[6*i+5];

        //tij[6*i  ] = t11;
        tij[6*i+1] = ct2*t22+st2*t33-2.0*ctst*t23;
        tij[6*i+2] = ct2*t33+st2*t22+2.0*ctst*t23;

        tij[6*i+3] = (ct2-st2)*t23+ctst*(t22-t33);
        tij[6*i+4] = ct*t13+st*t12;
        tij[6*i+5] = ct*t12-st*t13;

      }
    }
      return;

    case FA_ZONE_PERIODIC_CYL_Z: {

      // sector periodicity about z...
      double ct = transform[0]; // cos(theta)
      double st = transform[1]; // sin(theta)
      double ct2 = ct*ct;
      double st2 = st*st;
      double ctst = ct*st;

      // the following is the expansion of [R][T][Rt]
      // where [R] = [[ct, -st,  0]
      //              [st,  ct,  0]
      //              [ 0,   0,  1]]

      for (int i = 0; i<n; i++) {

        double t11 = tij[6*i];
        double t22 = tij[6*i+1];
        //double t33 = tij[6*i+2];

        double t23 = tij[6*i+3];
        double t13 = tij[6*i+4];
        double t12 = tij[6*i+5];

        tij[6*i] = ct2*t11-2.0*ctst*t12+st2*t22;
        tij[6*i+1] = st2*t11+2.0*ctst*t12+ct2*t22;
        //tij[6*i+2] = t33;

        tij[6*i+3] = st*t13+ct*t23;
        tij[6*i+4] = ct*t13-st*t23;
        tij[6*i+5] = ctst*(t11-t22)+(ct2-st2)*t12;

      }
    }
      return;

    default:
      cerr<<"Error: unsupported range flag: "<<flag<<endl;
      throw(-1);
    }
  }

public:

  int getNcv() const {
    return (ncv);
  }
  int getNfa() const {
    return (nfa);
  }
  int getNno() const {
    return (nno);
  }

  //void reconnectAndRemove();

protected:

  virtual void transformMeshHook() {
    // overload this to modify the mesh in some way
    // if you do overload it and you do transform the mesh, consider calling
    // rebuildPeriodicData()
  }

private:

  void setCvPartMetis();
  void setCvGroupPartMetis();
  void setCvPartCheap();
  void setCvPartBad();
  void checkCvPart();
  void redistReconnectReorder();

private:

  void syncNoofa();

  // ============================================
  // node reconnection stuff, as well as the common routine "buildPrcommsAndRanges"
  // that builds parallel connectivity with potentially multiple
  // periodicity from a Doano vector
  // ============================================

private:

  void reconnectNodes();

protected:

  // stuff used for node reconnection, but also for Cv2 reconnection, or any other
  // time when multiple periodicity reconnection is performed...

  typedef struct {
    int ino, ino_nbr, rank_nbr, bits;
  } Daono;

  // eventually it is neccessary to sort the Daono vector by rank, then
  // flag (i.e. transformation), then local ino to making packing as contiguous
  // as possible.

  class DaonoCompare: public std::binary_function<Daono&, Daono&, bool> {
  public:
    bool operator()(const Daono& a, const Daono& b) {
      return ((a.rank_nbr<b.rank_nbr)||((a.rank_nbr==b.rank_nbr)&&(a.bits<b.bits))||((a.rank_nbr
          ==b.rank_nbr)&&(a.bits==b.bits)&&(a.ino<b.ino)));
    }
  };

  void buildPrcommsAndRanges(vector<Daono>& daono_v, list<Prcomm>& prcommList);

private:

  void checkCvonoCounts();

public:

  void checkFaceCorrespondence();
  void checkNodeCorrespondence();
  void syncNodes();

protected:

  void buildFaocv();
  void buildNoocv();
  void reorderFaocvNoocv();

private:

  int reorderElementNodesAndFaces(int * no_flag_cv, int * fa_flag_cv, const int nno_cv, const int nfa_cv,
      const double(* x_no_cv)[3], const int * noofa_i_cv, const int * noofa_v_cv);
  int getQuadFaceMatching(int &ino0, int &ino1, int &ino2, int &ino3, const int nfa_ss,
      const int * noofa_i_ss, const int * noofa_v_ss);
  int getTriFaceMatching(int &ino0, int &ino1, int &ino2, const int nfa_ss, const int * noofa_i_ss,
      const int * noofa_v_ss);

public:

  void writeBoundaryFacesTecplot(const char * prefix);
  void writeBoundaryFacesStl();
  void writeFlaggedFacesTecplot(char * filename);
  void writeFlaggedCvsTecplotASCII(char * filename);
  void writeFlaggedCvsTecplot(char * filename);

private:

  void buildFaKindTable() {

    assert(fa_kind_table==NULL);

    // recall that fa_flag now contains the face zone of each face - we
    // cannot assume anything about face ordering at present.

    int n_fa_kind_table = -1;
    list<FaZone>::iterator iz;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
      n_fa_kind_table = max(n_fa_kind_table, iz->getIndex());
    assert(n_fa_kind_table>=0);
    n_fa_kind_table++; // zero indexing

    // make sure all fa_flag's fall in this range...
    for (int ifa = 0; ifa<nfa; ifa++) {
      assert((fa_flag[ifa]>=0)&&(fa_flag[ifa]<n_fa_kind_table));
    }

    fa_kind_table = new int[n_fa_kind_table];
    for (int i = 0; i<n_fa_kind_table; i++)
      fa_kind_table[i] = -1;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++) {
      int index = iz->getIndex();
      int kind = iz->getKind();
      assert(fa_kind_table[index]==-1);
      fa_kind_table[index] = kind;
    }

  }

  void clearFaKindTable() {

    if (fa_kind_table!=NULL) {
      delete[] fa_kind_table;
      fa_kind_table = NULL;
    }

  }

public:

  void writeFluentMsh(char * filename);

  void flagCvsIntersectingPlane(const double x, const double y, const double z, const double nx,
      const double ny, const double nz, const double tol) {

    for (int ino = 0; ino<getNno(); ino++) {
      double delta = (x_no[ino][0]-x)*nx+(x_no[ino][1]-y)*ny+(x_no[ino][2]-z)*nz;
      if (delta>tol) {
        no_flag[ino] = 1;
      }
      else if (delta<-tol) {
        no_flag[ino] = -1;
      }
      else {
        no_flag[ino] = 0;
      }
    }

    for (int icv = 0; icv<ncv; icv++) {
      // use the flag to store the status of nodes so far...
      // flag =  0 : no status
      // flag = -1 : one or more -1 nodes found
      // flag = +1 : one or more +1 nodes found
      // flag = +2 : flag this cv
      int flag = 0;
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; ++foc) {
        int ifa = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; ++nof) {
          int ino = noofa_v[nof];
          switch (no_flag[ino]) {
          case 0:
            flag = 2;
            break;
          case -1:
            switch (flag) {
            case 0:
              flag = -1;
              break;
            case 1:
              flag = 2;
              break;
            }
            break;
          case 1:
            switch (flag) {
            case 0:
              flag = 1;
              break;
            case -1:
              flag = 2;
              break;
            }
            break;
          }
          if (flag==2) break;
        }
        if (flag==2) break;
      }
      if (flag==2) cv_flag[icv] = 1;
    }

  }

protected:

  int deleteFileCollective(char * filename) {

    // This routine returns 0 if the file was found and deleted
    // successfully, non-zero otherwise...

    int ierr;
    if (mpi_rank==0) {
      // only root process tries to delete...
      ierr = MPI_File_delete(filename, MPI_INFO_NULL);
      if (ierr==0) cout<<"deleteFileCollective: "<<filename<<endl;
    }
    MPI_Bcast(&ierr, 1, MPI_INT, 0, mpi_comm);
    return (ierr);

  }

  int solveScalarCgLocal(double * phi, double * A, double * rhs, const int nno, int * nbono_i,
      int * nbono_v, const double zero, const int maxiter) {

    // purely local cg solve - I hope you know what you are doing...
    int mode = ABSOLUTE_RESIDUAL;

    // we need the following work arrays...
    double * res = new double[nno];
    double * v = new double[nno];
    double * p = new double[nno];

    // initialize...
    for (int ino = 0; ino<nno; ino++)
      p[ino] = 0.0;

    double rho = 1.0;

    // calculate the residual in rhs format...
    // also include a few assertions to check if everything
    // is set up properly for local solutions. If these
    // fail, check your nbono_i and nbono_v carefully - are you sure you shouldn't
    // be using a parallel solver?

    for (int ino = 0; ino<nno; ino++) {
      int non_f = nbono_i[ino];
      // diagonal check...
      assert(nbono_v[non_f]==ino);
      int non_l = nbono_i[ino+1]-1;
      res[ino] = rhs[ino]-A[non_f]*phi[ino]; // diagonal
      // neighbors...
      for (int non = non_f+1; non<=non_l; non++) {
        int ino_nbr = nbono_v[non];
        assert((ino_nbr>=0)&&(ino_nbr<nno));
        res[ino] -= A[non]*phi[ino_nbr];
      }
    }

    double res0_max = 0.0;
    if (mode==RELATIVE_RESIDUAL) for (int ino = 0; ino<nno; ino++)
      res0_max = max(res0_max, fabs(res[ino]/A[nbono_i[ino]]));

    int iter = 0;
    int done = 0;
    while (done==0) {

      iter++;

      // diagonal precon...
      for (int ino = 0; ino<nno; ino++)
        v[ino] = res[ino]/A[nbono_i[ino]];

      double rho_prev = rho;
      rho = 0.0;
      for (int ino = 0; ino<nno; ino++)
        rho += res[ino]*v[ino];

      if (fabs(rho_prev)<1.0E-25) rho_prev = 1.0E-25;
      double beta = rho/rho_prev;

      for (int ino = 0; ino<nno; ino++)
        p[ino] = v[ino]+beta*p[ino];

      // v = [A]{p}...
      for (int ino = 0; ino<nno; ino++) {
        int non_f = nbono_i[ino];
        int non_l = nbono_i[ino+1]-1;
        v[ino] = A[non_f]*p[ino]; // diagonal
        for (int non = non_f+1; non<=non_l; non++)
          v[ino] += A[non]*p[nbono_v[non]];
      }

      double gamma = 0.0;
      for (int ino = 0; ino<nno; ino++)
        gamma += p[ino]*v[ino];

      if (fabs(gamma)<1.0E-25) gamma = 1.0E-25;

      double alpha = rho/gamma;
      for (int ino = 0; ino<nno; ino++)
        phi[ino] += alpha*p[ino];

      // check if we are done...
      if (iter%5==0) {

        // on the other iterations, use this approximation...
        /*
         for (int ino = 0; ino < nno; ino++)
         res[ino] -= alpha*v[ino];
         */

        // recompute the residual...
        double res_max = 0.0;
        for (int ino = 0; ino<nno; ino++) {
          int non_f = nbono_i[ino];
          int non_l = nbono_i[ino+1]-1;
          res[ino] = rhs[ino]-A[non_f]*phi[ino]; // diagonal
          // neighbors...
          for (int non = non_f+1; non<=non_l; non++)
            res[ino] -= A[non]*phi[nbono_v[non]];
          res_max = max(res_max, fabs(res[ino]/A[nbono_i[ino]]));
        }

        //cout << "cg iter, res_max: " << iter << " " << res_max << endl;
        if ((mode==ABSOLUTE_RESIDUAL)&&(res_max<=zero)) {
          done = 1;
        }
        else if ((mode==RELATIVE_RESIDUAL)&&(res_max/(res0_max+1.0E-20)<=zero)) {
          done = 1;
        }
        else if (iter>maxiter) {
          cout<<"Warning: solveScalarCgLocal did not converge after "<<maxiter<<" iters, res_max: "
              <<res_max<<endl;
          done = 2;
        }
        else { //if (iter > maxiter/2) {
          cout<<" > solveScalarCgLocal: iter, res_max: "<<iter<<" "<<res_max<<endl;
        }

      }
      else {

        // on the other iterations, use this approximation...
        for (int ino = 0; ino<nno; ino++)
          res[ino] -= alpha*v[ino];
      }

    }

    //delete[]
    delete[] res;
    delete[] v;
    delete[] p;

    // let the calling routine know if we were successful...
    return (done==1);

  }

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // i/o stuff...
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

public:

  void writeRestart() {
    writeRestart("restart.out");
  }
  void writeRestart(const int index) {
    char filename[32];
    sprintf(filename, "restart.%06d.out", index);
    writeRestart(filename);
  }
  void writeRestart(const string& filename) {
    writeRestart(filename.c_str());
  }
  void writeRestart(const char * filename);

  int getIoWriteFlag() const {
    return (io_write_flag);
  }

  void writeSnapshot(const int index);
  void readSnapshot(const int index);

private:

#define UGP_IO_HEADER_NAME_LEN    52 // select this to make header size 256 bytes
  /*
   typedef struct {
   char name[UGP_IO_HEADER_NAME_LEN];
   int id;
   MPI_Offset skip;
   int idata[16];
   double rdata[16];
   } Header;
   */

  // use the following header class to avoid unitialized errors in valgrind...

  class Header {
  public:
    char name[UGP_IO_HEADER_NAME_LEN];
    int id;
    MPI_Offset skip;
    int idata[16];
    double rdata[16];
    Header() {
      for (int i = 0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        name[i] = 0;
      id = 0;
      skip = 0;
      for (int i = 0; i<16; ++i) {
        idata[i] = 0;
        rdata[i] = 0.0;
      }
    }
  };

  void initIoCommon();

  void byteSwapHeader(Header * header, const int n) {
    for (int i = 0; i<n; i++) {
      header[i].id = byteSwap(header[i].id);
      header[i].skip = byteSwap(header[i].skip);
      byteSwap(header[i].idata, 16);
      byteSwap(header[i].rdata, 16);
    }
  }

  void writeI0(const int data, const int id, char * name);
  void writeD0(const double data, const int id, char * name);

  void writeCVsI1(int *data, const int id, char * name);
  void writeCVsD1(double *data, const int id, char * name);
  void writeCVsD3(double(*data)[3], const int id, char * name);
  void writeCVsD33(double(*data)[3][3], const int id, char * name);

  void writeNodeI1(int *data, const int id, char * name);
  void writeNodeD1(double *data, const int id, char * name);
  void writeNodeD3(double(*data)[3], const int id, char * name);

  void writeFaceI1(int *data, const int id, char * name);
  void writeFaceD1(double *data, const int id, char * name);
  void writeNoofa();
  void writeCvofa();

  void writeCvI1(int *data, const int id, char * name);
  void writeCvPart();

public:

  void readRestart() {
    readRestart("restart.in");
  }
  void readRestart(const string& filename) {
    readRestart(filename.c_str());
  }
  void readRestart(const char * filename);

  // these routines can read data node data from a restart file. The user
  // requests a certain range of data, and the routine returns the actual
  // range.

  int readRestartNoX(double(*x_no)[3], const int i0, const int i1, const char * filename);

  int readRestartNoXCollective(double(*x_no)[3], const int i0, const int i1, const char * filename);

  int readRestartNoR1(double * data, const int i0, const int i1, const char * name, const char * filename);

  int readRestartNoR2(double(*data)[3], const int i0, const int i1, const char * name,
      const char * filename);

  int readRestartNoR2Collective(double(*data)[3], const int i0, const int i1, const char * name,
      const char * filename);

private:

  void initReadRestart();
  void readNodeCheck(Header& header);
  void readFaceCheck(Header& header);
  void readCvCheck(Header& header);
  void readNoofa(Header& header);
  void readCvofa(Header& header);
  void readFaZone(Header& header);
  void readFaD1(double *data, Header& header);
  void readCvPart(Header& header);
  void readCvD1(double *data, Header& header);
  void readCvD3(double(*data)[3], Header& header);
  void readCvD33(double(*data)[3][3], Header& header);
  void readNoD1(double *data, Header& header);
  void readNoD3(double(*data)[3], Header& header);
  void finalizeReadRestart();

  void cleanupReadRestart();

};

#endif
