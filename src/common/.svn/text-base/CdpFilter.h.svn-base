#ifndef CDPFILTER_H
#define CDPFILTER_H

#include "MiscUtils.h"
using namespace MiscUtils;

#include "Ugp.h"

// =====================================================================
// note that the CdpFilter class is a friend of the Ugp class, so it can
// set Ugp's private data, as it must during the building process...
// =====================================================================

#define CDP_IO_MAGIC_NUMBER  1235813
#define CDP_IO_VERSION             4

#define CDP_IO_R0                 10
#define CDP_IO_I0                 11
#define CDP_IO_NO_COUNTS         101
#define CDP_IO_NODE_CC           102
#define CDP_IO_NO_FLAG_TEST      103
#define CDP_IO_NO_R1             104
#define CDP_IO_NO_R2             105
#define CDP_IO_NO_BA_R1          106
#define CDP_IO_NO_BA_R2          107
#define CDP_IO_NO_I1             108
#define CDP_IO_FA_COUNTS         201
#define CDP_IO_FA_FLAG_TEST      202
#define CDP_IO_FA_ZONE           203
#define CDP_IO_NOOFA_AND_EDOFA   204
#define CDP_IO_CVOFA             205
#define CDP_IO_FA_R1             206
#define CDP_IO_FA_R2             207
#define CDP_IO_FA_BA_R1          208
#define CDP_IO_FA_BA_R2          209
#define CDP_IO_FA_BA_I1          210
#define CDP_IO_ED_COUNTS         301
#define CDP_IO_ED_FLAG_TEST      302
#define CDP_IO_NOOED             303
#define CDP_IO_ED_R1             304
#define CDP_IO_ED_R2             305
#define CDP_IO_CV_COUNTS         401
#define CDP_IO_CV_FLAG_TEST      402
#define CDP_IO_CV_ZONE           403
#define CDP_IO_CV_PART           404
#define CDP_IO_CV_R1             405
#define CDP_IO_CV_R2             406
#define CDP_IO_LP_X              501
#define CDP_IO_LP_R0             502
#define CDP_IO_LP_R1             503
#define CDP_IO_LP_R2             504       
#define CDP_IO_LAST             1001

#define CDP_IO_HEADER_NAME_LEN    52 // select this to make header size 256 bytes
#define CDP_FA_ZONE_I             1
#define CDP_FA_ZONE_BI            2
#define CDP_FA_ZONE_BP            3
#define CDP_FA_ZONE_BA            4

#define CDP_PERIODIC_CART         1 
#define CDP_PERIODIC_CYL_X        2
#define CDP_PERIODIC_CYL_Y        3
#define CDP_PERIODIC_CYL_Z        4

class CdpFilter {

private:

  string filename;

  typedef struct {
    char name[CDP_IO_HEADER_NAME_LEN];
    int id;
    MPI_Offset skip;
    int idata[16];
    double rdata[16];
  } CdpHeader;

  MPI_Datatype MPI_CdpHeader;
  MPI_Aint header_extent;

  MPI_Datatype no_scalar_int_t;
  MPI_Datatype no_scalar_double_t;
  MPI_Datatype no_vector_double_t;

  MPI_Datatype fa_scalar_int_t;
  int ifa_global;

  MPI_Datatype cv_scalar_int_t;
  int icv_global;

  int byte_swap;

  Ugp * ugp;
  MPI_File fh;
  MPI_Offset offset;
  CdpHeader header;

public:

  CdpFilter(const string filename) {

    if (mpi_rank == 0)
      cout << "CdpFilter(): filename: "<< filename << endl;

    // nullify the ugp pointer...
    ugp = NULL;

    // store the filename...
    assert(filename.length() < 128);
    this->filename = filename;

    // check some sizes...
    if ((sizeof(int) != 4) || (sizeof(MPI_Offset) != 8)|| (sizeof(float) != 4)
        || (sizeof(double) != 8)) {
      if (mpi_rank == 0) {
        cerr << "Error: unexpected sizes:"<< endl;
        cerr << "int       : "<< sizeof(int)<< endl;
        cerr << "MPI_Offset: "<< sizeof(MPI_Offset)<< endl;
        cerr << "float     : "<< sizeof(float)<< endl;
        cerr << "double    : "<< sizeof(double)<< endl;
      }
      throw(-1);
    }

  }

  ~CdpFilter() {
    if (mpi_rank == 0)
      cout << "~CdpFilter()"<< endl;
  }

  void minInitUgp(Ugp * ugp) {
    // store a pointer to the ugp...
    this->ugp = ugp;

    // read the specified restart file and 
    // produce a minimum initialization of the Ugp...
    readCdpRestart();
  }

private:

  void readCdpRestart();

  // to produce an approximately uniform distribution of nX. The
  // user should have preallocated/newed xod with size [ndist+1]...
  // the zero-indexed global range for each zero-indexed id is then 
  // xod[id] .. xod[id+1]-1 inclusive 

  void byteSwapCdpHeader(CdpHeader * header, const int n) {
    for (int i = 0; i < n; i++) {
      header[i].id = byteSwap(header[i].id);
      header[i].skip = byteSwap(header[i].skip);
      byteSwap(header[i].idata, 16);
      byteSwap(header[i].rdata, 16);
    }
  }

private:

  void build_no_structs();

  void no_flag_test();

  void build_fa_structs();

  void fa_flag_test();

  void build_cv_structs();

  void cv_flag_test();

  void read_cv_part();

  void add_face_zone();

  void read_noofa_iv();

  void read_cvofa();

  void read_no_scalar(double *s);

  void read_no_vector(double (*v)[3]);

  void finalize_and_check();

};

#endif

