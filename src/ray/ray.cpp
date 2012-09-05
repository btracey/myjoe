
#include "postuq.h"

#if 0
#include "Param.h"
#include "UgpWithTools.h"
#include "tc_vec3d.h"


class PostUQ: public UgpWithTools
{
private:

  int n_cv_pack;
  int * cv_pack_count;
  int * cv_pack_disp;
  int * cv_pack_index;
  int * cv_pack_buffer_int;
  double * cv_pack_buffer_double;

  int n_cv_unpack;
  int * cv_unpack_count;
  int * cv_unpack_disp;
  int * cv_unpack_index;
  int * cv_unpack_buffer_int;
  double * cv_unpack_buffer_double;

protected:

  double (*x_cv)[3];
  PostUQ *ugp2;

public:

  PostUQ()
  {
    if (mpi_rank==0) cout<<"PostUQ()"<<endl;

    n_cv_pack = 0;
    cv_pack_count = NULL;
    cv_pack_disp = NULL;
    cv_pack_index = NULL;
    cv_pack_buffer_int = NULL;
    cv_pack_buffer_double = NULL;

    n_cv_unpack = 0;
    cv_unpack_count = NULL;
    cv_unpack_disp = NULL;
    cv_unpack_index = NULL;
    cv_unpack_buffer_int = NULL;
    cv_unpack_buffer_double = NULL;

    // used to match cv-data...
    x_cv = NULL;

    // ugp2 carries the restart file that might have different partitioning
    ugp2 = NULL;
  }

  ~PostUQ()
  {
    // cleanup all the cv_pack and cv_unpack stuff...
    deleteNoInitCvDataMap();

    // x_cv...
    if (x_cv!=NULL) {
      delete[] x_cv;
      x_cv = NULL;
    }
  }

  void buildNoInitCvDataMap(double(*x_cv2)[3], int * cvora2)
  {
    if (mpi_rank==0) cout<<"buildNoInitCvDataMap()"<<endl;

    // cleanup from any previous map...
    deleteNoInitCvDataMap();

    assert(cvora[mpi_size]==cvora2[mpi_size]);
    int ncv2 = cvora2[mpi_rank+1]-cvora2[mpi_rank];

    // prepare to send the x_cv2's to ranks where there may be
    // a match. The send side is considered the data layout associated
    // with x_cv2, ncv2, cvora2...

    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    double * send_buffer_x_cv2 = NULL;
    int send_count_sum;
    for (int i = 0; i<mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter<2; iter++) {
      for (int icv2 = 0; icv2<ncv2; ++icv2) {
        for (int rank = 0; rank<mpi_size; ++rank) {
          if ((x_cv2[icv2][0]>=cvAdtBBMinOfRank[rank][0])&&(x_cv2[icv2][0]<=cvAdtBBMaxOfRank[rank][0])
            &&(x_cv2[icv2][1]>=cvAdtBBMinOfRank[rank][1])&&(x_cv2[icv2][1]<=cvAdtBBMaxOfRank[rank][1])
            &&(x_cv2[icv2][2]>=cvAdtBBMinOfRank[rank][2])&&(x_cv2[icv2][2]<=cvAdtBBMaxOfRank[rank][2])) {
            if (iter==0) {
              send_count[rank] += 3;
            }
            else {
              send_buffer_x_cv2[send_disp[rank]] = x_cv2[icv2][0];
              send_buffer_x_cv2[send_disp[rank]+1] = x_cv2[icv2][1];
              send_buffer_x_cv2[send_disp[rank]+2] = x_cv2[icv2][2];
              send_disp[rank] += 3;
            }
          }
        }
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i<mpi_size; i++)
        send_disp[i] = send_count[i-1]+send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter==0) {
        send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
        send_buffer_x_cv2 = new double[send_count_sum];
      }
    }

    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
    double * recv_buffer_x_cv2 = new double[recv_count_sum];
    MPI_Alltoallv(send_buffer_x_cv2, send_count, send_disp, MPI_DOUBLE, recv_buffer_x_cv2, recv_count, recv_disp, MPI_DOUBLE, mpi_comm);

    // the x_cv2 data is now on the recv side. Use the cvAdt to identify point matches
    // and return the closest match we can find. since we are now going to send one
    // value per request rather than 3, remove the factor of 3...
    send_count_sum /= 3;
    recv_count_sum /= 3;
    for (int i = 0; i<mpi_size; i++) {
      send_count[i] /= 3;
      send_disp[i] /= 3;
      recv_count[i] /= 3;
      recv_disp[i] /= 3;
    }

    delete[] send_buffer_x_cv2;

    double * recv_buffer_d2 = new double[recv_count_sum];
    int * recv_buffer_icv = new int[recv_count_sum];

    int bbListSize, bbList[ADT_LIST_MAX];
    for (int irecv = 0; irecv<recv_count_sum; ++irecv) {
      double this_x_cv2[3];
      FOR_I3
        this_x_cv2[i] = recv_buffer_x_cv2[3*irecv+i];
      recv_buffer_icv[irecv] = -1;
      // use the Adt to get a list of possible candidates...
      cvAdt->buildListForPoint(bbListSize, bbList, this_x_cv2);
      // cycle through candidates and store the closest...
      for (int ibb = 0; ibb<bbListSize; ++ibb) {
        const int icv = bbList[ibb];
        double d2 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv][i]-this_x_cv2[i];
          d2 += dx*dx;
        }
        if ((recv_buffer_icv[irecv]==-1)||(d2<recv_buffer_d2[irecv])) {
          recv_buffer_icv[irecv] = icv; // the local icv is fine
          recv_buffer_d2[irecv] = d2;
        }
      }
    }

    delete[] recv_buffer_x_cv2;

    double * send_buffer_d2 = new double[send_count_sum];
    int * send_buffer_icv = new int[send_count_sum];

    // send the recv stuff nack to the send data layout...
    MPI_Alltoallv(recv_buffer_icv, recv_count, recv_disp, MPI_INT, send_buffer_icv, send_count, send_disp, MPI_INT,
        mpi_comm);
    MPI_Alltoallv(recv_buffer_d2, recv_count, recv_disp, MPI_DOUBLE, send_buffer_d2, send_count, send_disp, MPI_DOUBLE,
        mpi_comm);

    // HERE - do a double loop, although there is likely some shorter way
    // to do this because we know that the data layouts are 1:1 matched. Other topologies
    // are not single-valued (nodes, edges, faces), so this is more general...
    assert(cv_pack_count==NULL);
    cv_pack_count = new int[mpi_size];
    assert(cv_pack_disp==NULL);
    cv_pack_disp = new int[mpi_size];
    // these get set below...
    assert(cv_pack_index==NULL);
    assert(cv_pack_buffer_int==NULL);
    for (int i = 0; i<mpi_size; i++)
      cv_pack_count[i] = 0;
    double my_d2_max = 0.0;
    for (int iter = 0; iter<2; ++iter) {
      for (int icv2 = 0; icv2<ncv2; ++icv2) {
        int rank_min = -1;
        int icv_min;
        double d2_min;
        for (int rank = 0; rank<mpi_size; ++rank) {
          if ((x_cv2[icv2][0]>=cvAdtBBMinOfRank[rank][0])&&(x_cv2[icv2][0]<=cvAdtBBMaxOfRank[rank][0])&&(x_cv2[icv2][1]
              >=cvAdtBBMinOfRank[rank][1])&&(x_cv2[icv2][1]<=cvAdtBBMaxOfRank[rank][1])&&(x_cv2[icv2][2]
              >=cvAdtBBMinOfRank[rank][2])&&(x_cv2[icv2][2]<=cvAdtBBMaxOfRank[rank][2])) {
            double this_d2 = send_buffer_d2[send_disp[rank]];
            int this_icv = send_buffer_icv[send_disp[rank]];
            send_disp[rank] += 1;
            if (((rank_min==-1)||(this_d2<d2_min))&&(this_icv>=0)) {
              rank_min = rank;
              icv_min = this_icv;
              d2_min = this_d2;
            }
          }
        }
        // if we don't find even one, then this is a problem...
        assert(rank_min>=0);
        // we need to build a pack/unpack buffer that knows how to
        // pack the data quickly...
        if (iter==0) cv_pack_count[rank_min] += 1;
        else {
          // store the max d2 (of all our min's) for reporting...
          my_d2_max = max(my_d2_max, d2_min);
          // fill index and buffer, inc disp...
          cv_pack_index[cv_pack_disp[rank_min]] = icv2; // local value
          cv_pack_buffer_int[cv_pack_disp[rank_min]] = icv_min; // send to rank_min below
          cv_pack_disp[rank_min] += 1;
        }
      }
      // calculate cv_pack_disp on both iterations...
      cv_pack_disp[0] = 0;
      for (int i = 1; i<mpi_size; i++)
        cv_pack_disp[i] = cv_pack_count[i-1]+cv_pack_disp[i-1];
      // on the first time through...
      if (iter==0) {
        // rebuild send_disp...
        send_disp[0] = 0;
        for (int i = 1; i<mpi_size; i++)
          send_disp[i] = send_count[i-1]+send_disp[i-1];
        // allocate pack buffers...
        n_cv_pack = cv_pack_disp[mpi_size-1]+cv_pack_count[mpi_size-1];
        cv_pack_index = new int[n_cv_pack];
        cv_pack_buffer_int = new int[n_cv_pack];
      }
    }

    double d2_max;
    MPI_Reduce(&my_d2_max, &d2_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank==0) cout<<" > all cvs matched within tol: "<<sqrt(d2_max)<<endl;

    // cleanup...
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    delete[] send_buffer_d2;
    delete[] send_buffer_icv;
    delete[] recv_buffer_d2;
    delete[] recv_buffer_icv;

    // now build the unpack stuff...
    assert(cv_unpack_count==NULL);
    cv_unpack_count = new int[mpi_size];
    MPI_Alltoall(cv_pack_count, 1, MPI_INT, cv_unpack_count, 1, MPI_INT, mpi_comm);
    assert(cv_unpack_disp==NULL);
    cv_unpack_disp = new int[mpi_size];
    cv_unpack_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      cv_unpack_disp[i] = cv_unpack_count[i-1]+cv_unpack_disp[i-1];
    n_cv_unpack = cv_unpack_disp[mpi_size-1]+cv_unpack_count[mpi_size-1];
    assert(cv_unpack_index==NULL);
    cv_unpack_index = new int[n_cv_unpack];
    // exchange directly into the unpack index...
    MPI_Alltoallv(cv_pack_buffer_int, cv_pack_count, cv_pack_disp, MPI_INT, cv_unpack_index, cv_unpack_count,
        cv_unpack_disp, MPI_INT, mpi_comm);

    // allocate remaining buffers...
    assert(cv_unpack_buffer_int==NULL);
    cv_unpack_buffer_int = new int[n_cv_unpack];
    assert(cv_pack_buffer_double==NULL);
    cv_pack_buffer_double = new double[n_cv_pack];
    assert(cv_unpack_buffer_double==NULL);
    cv_unpack_buffer_double = new double[n_cv_unpack];

    // check: try comparing components of x_cv2...
    FOR_I3 {
      for (int ii = 0; ii<n_cv_pack; ++ii) {
        int icv2 = cv_pack_index[ii];
        cv_pack_buffer_double[ii] = x_cv2[icv2][i];
      }
      MPI_Alltoallv(cv_pack_buffer_double, cv_pack_count, cv_pack_disp, MPI_DOUBLE, cv_unpack_buffer_double,
          cv_unpack_count, cv_unpack_disp, MPI_DOUBLE, mpi_comm);
      double my_d2_max = 0.0;
      for (int ii = 0; ii<n_cv_unpack; ++ii) {
        int icv = cv_unpack_index[ii];
        double dx = cv_unpack_buffer_double[ii]-x_cv[icv][i];
        my_d2_max = max(my_d2_max, dx*dx);
      }
      double d2_max;
      MPI_Reduce(&my_d2_max, &d2_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank==0) cout<<" > check for component: "<<i<<" tol (should be zero): "<<sqrt(d2_max)<<endl;
    }

    // also check that all cv's are in fact touched and touched
    // once during uppacking...
    FOR_ICV
      cv_flag[icv] = 0;
    for (int ii = 0; ii<n_cv_unpack; ++ii) {
      int icv = cv_unpack_index[ii];
      assert(cv_flag[icv]==0);
      cv_flag[icv] = 1;
    }
    FOR_ICV
      assert(cv_flag[icv]==1);

  }

  void deleteNoInitCvDataMap()
  {
    n_cv_pack = 0;
    if (cv_pack_count!=NULL) {
      delete[] cv_pack_count;
      cv_pack_count = NULL;
    }
    if (cv_pack_disp!=NULL) {
      delete[] cv_pack_disp;
      cv_pack_disp = NULL;
    }
    if (cv_pack_index!=NULL) {
      delete[] cv_pack_index;
      cv_pack_index = NULL;
    }
    if (cv_pack_buffer_int!=NULL) {
      delete[] cv_pack_buffer_int;
      cv_pack_buffer_int = NULL;
    }
    if (cv_pack_buffer_double!=NULL) {
      delete[] cv_pack_buffer_double;
      cv_pack_buffer_double = NULL;
    }

    n_cv_unpack = 0;
    if (cv_unpack_count!=NULL) {
      delete[] cv_unpack_count;
      cv_unpack_count = NULL;
    }
    if (cv_unpack_disp!=NULL) {
      delete[] cv_unpack_disp;
      cv_unpack_disp = NULL;
    }
    if (cv_unpack_index!=NULL) {
      delete[] cv_unpack_index;
      cv_unpack_index = NULL;
    }
    if (cv_unpack_buffer_int!=NULL) {
      delete[] cv_unpack_buffer_int;
      cv_unpack_buffer_int = NULL;
    }
    if (cv_unpack_buffer_double!=NULL) {
      delete[] cv_unpack_buffer_double;
      cv_unpack_buffer_double = NULL;
    }
  }

  void updateNoInitCvData(double * scalar, const double * noInitScalar)
  {
    for (int ii = 0; ii<n_cv_pack; ++ii) {
      int icv2 = cv_pack_index[ii];
      cv_pack_buffer_double[ii] = noInitScalar[icv2];
    }

    MPI_Alltoallv(cv_pack_buffer_double, cv_pack_count, cv_pack_disp, MPI_DOUBLE, cv_unpack_buffer_double,
        cv_unpack_count, cv_unpack_disp, MPI_DOUBLE, mpi_comm);

    for (int ii = 0; ii<n_cv_unpack; ++ii) {
      int icv = cv_unpack_index[ii];
      scalar[icv] = cv_unpack_buffer_double[ii];
    }

  }

  void updateNoInitCvData(double(*vector)[3], const double(*noInitVector)[3])
  {
    FOR_I3 {
      for (int ii = 0; ii<n_cv_pack; ++ii) {
        int icv2 = cv_pack_index[ii];
        cv_pack_buffer_double[ii] = noInitVector[icv2][i];
      }

      MPI_Alltoallv(cv_pack_buffer_double, cv_pack_count, cv_pack_disp, MPI_DOUBLE, cv_unpack_buffer_double,
          cv_unpack_count, cv_unpack_disp, MPI_DOUBLE, mpi_comm);

      for (int ii = 0; ii<n_cv_unpack; ++ii) {
        int icv = cv_unpack_index[ii];
        vector[icv][i] = cv_unpack_buffer_double[ii];
      }
    }

  }

  void calcCvCenterInit(double(* x_cv_approx)[3])
  {
    FOR_ICV {
      x_cv_approx[icv][0] = 0.0;
      x_cv_approx[icv][1] = 0.0;
      x_cv_approx[icv][2] = 0.0;
    }

    FOR_IFA {
      double x_fa_approx[3];
      x_fa_approx[0] = 0.0;
      x_fa_approx[1] = 0.0;
      x_fa_approx[2] = 0.0;
      //Loop over all nodes of each face and sum coordinates
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1];
      for (int nof = nof_f; nof<nof_l; ++nof) {
        int ino = noofa_v[nof];
        x_fa_approx[0] += x_no[ino][0];
        x_fa_approx[1] += x_no[ino][1];
        x_fa_approx[2] += x_no[ino][2];
      }
      //Divide by number of nodes per face
      double tmp = 1.0/(double) (nof_l-nof_f);
      x_fa_approx[0] *= tmp;
      x_fa_approx[1] *= tmp;
      x_fa_approx[2] *= tmp;
      // send to both cv's...
      int icv0 = cvofa[ifa][0];
      assert((icv0>=0)&&(icv0<ncv));
      x_cv_approx[icv0][0] += x_fa_approx[0];
      x_cv_approx[icv0][1] += x_fa_approx[1];
      x_cv_approx[icv0][2] += x_fa_approx[2];
      int icv1 = cvofa[ifa][1];
      if ((icv1>=0)&&(icv1<ncv)) {
        x_cv_approx[icv1][0] += x_fa_approx[0];
        x_cv_approx[icv1][1] += x_fa_approx[1];
        x_cv_approx[icv1][2] += x_fa_approx[2];
      }
    }

    // divide by the number of faces per cell...
    FOR_ICV {
      double tmp = 1.0/(double) (faocv_i[icv+1]-faocv_i[icv]);
      x_cv_approx[icv][0] *= tmp;
      x_cv_approx[icv][1] *= tmp;
      x_cv_approx[icv][2] *= tmp;
    }

  }

  void calcCvCenterNoInit(double(* x_cv_approx)[3])
  {
    if (mpi_rank==0) cout<<"calcCvCenterNoInit()"<<endl;

    // we use the noora and cvora here...
    assert(noora[mpi_rank+1]-noora[mpi_rank]==nno);
    assert(cvora[mpi_rank+1]-cvora[mpi_rank]==ncv);

    // Since the grid hasn't been redistributed, we need the calculate
    // approximate CV centers based on approximate face centers

    // Buffers..
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * send_buffer_int = NULL;
    int send_count_sum;
    for (int i = 0; i<mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter<2; iter++) {
      FOR_IFA {
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; ++nof) {
          // at this point, noofa_v stores the GLOBAL node index...
          const int ino = noofa_v[nof];
          assert((ino>=0)&&(ino<noora[mpi_size]));
          int rank = 0;
          while (noora[rank+1]<=ino)
            ++rank;
          if (iter==0) send_count[rank] += 1;
          else {
            send_buffer_int[send_disp[rank]] = ino-noora[rank]; // local node index
            send_disp[rank] += 1;
          }
        }
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i<mpi_size; i++)
        send_disp[i] = send_count[i-1]+send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter==0) {
        send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
        send_buffer_int = new int[send_count_sum];
      }
    }

    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
    int * recv_buffer_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buffer_int, send_count, send_disp, MPI_INT, recv_buffer_int, recv_count, recv_disp, MPI_INT,
        mpi_comm);

    // some cleanup...
    delete[] send_buffer_int;

    // now load up the nodal coordinates...
    double (*recv_buffer_x)[3] = new double[recv_count_sum][3];
    for (int i = 0; i<recv_count_sum; ++i) {
      const int ino = recv_buffer_int[i];
      assert((ino>=0)&&(ino<nno));
      recv_buffer_x[i][0] = x_no[ino][0];
      recv_buffer_x[i][1] = x_no[ino][1];
      recv_buffer_x[i][2] = x_no[ino][2];
    }

    delete[] recv_buffer_int;

    // and send back: here we have three doubles per previous int...
    for (int i = 0; i<mpi_size; i++) {
      send_count[i] *= 3;
      send_disp[i] *= 3;
      recv_count[i] *= 3;
      recv_disp[i] *= 3;
    }
    double (*send_buffer_x)[3] = new double[send_count_sum][3];
    MPI_Alltoallv(recv_buffer_x, recv_count, recv_disp, MPI_DOUBLE, send_buffer_x, send_count, send_disp, MPI_DOUBLE,
        mpi_comm);

    // some more cleanup...
    delete[] recv_buffer_x;

    // return send_disp...
    for (int i = 0; i<mpi_size; i++)
      send_disp[i] /= 3;

    // x_fa...
    double (*x_fa_approx)[3] = new double[nfa][3];
    FOR_IFA {
      // zero this face's x_fa_approx...
      x_fa_approx[ifa][0] = 0.0;
      x_fa_approx[ifa][1] = 0.0;
      x_fa_approx[ifa][2] = 0.0;
      // loop on nodes...
      const int nof_f = noofa_i[ifa];
      const int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof<=nof_l; ++nof) {
        // at this point, noofa_v stores the GLOBAL node index...
        const int ino = noofa_v[nof];
        assert((ino>=0)&&(ino<noora[mpi_size]));
        int rank = 0;
        while (noora[rank+1]<=ino)
          ++rank;
        x_fa_approx[ifa][0] += send_buffer_x[send_disp[rank]][0];
        x_fa_approx[ifa][1] += send_buffer_x[send_disp[rank]][1];
        x_fa_approx[ifa][2] += send_buffer_x[send_disp[rank]][2];
        send_disp[rank] += 1;
      }
      double tmp = 1.0/(double) (nof_l-nof_f+1);
      x_fa_approx[ifa][0] *= tmp;
      x_fa_approx[ifa][1] *= tmp;
      x_fa_approx[ifa][2] *= tmp;
    }

    // more cleanup...
    delete[] send_buffer_x;

    // ===========================
    // now reduce to cv's...
    // ===========================

    for (int i = 0; i<mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter<2; iter++) {
      FOR_IFA {
        // icv0...
        const int icv0 = cvofa[ifa][0];
        assert((icv0>=0)&&(icv0<cvora[mpi_size]));
        int rank = 0;
        while (cvora[rank+1]<=icv0)
          ++rank;
        if (iter==0) send_count[rank] += 1;
        else {
          send_buffer_int[send_disp[rank]] = icv0-cvora[rank]; // local cv index
          send_buffer_x[send_disp[rank]][0] = x_fa_approx[ifa][0];
          send_buffer_x[send_disp[rank]][1] = x_fa_approx[ifa][1];
          send_buffer_x[send_disp[rank]][2] = x_fa_approx[ifa][2];
          send_disp[rank] += 1;
        }
        // icv1...
        // the fa_kind_table should be built at teh end of readRestart...
        if (fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL) {
          const int icv1 = cvofa[ifa][1];
          assert((icv1>=0)&&(icv1<cvora[mpi_size]));
          rank = 0;
          while (cvora[rank+1]<=icv1)
            ++rank;
          if (iter==0) send_count[rank] += 1;
          else {
            send_buffer_int[send_disp[rank]] = icv1-cvora[rank]; // local cv index
            send_buffer_x[send_disp[rank]][0] = x_fa_approx[ifa][0];
            send_buffer_x[send_disp[rank]][1] = x_fa_approx[ifa][1];
            send_buffer_x[send_disp[rank]][2] = x_fa_approx[ifa][2];
            send_disp[rank] += 1;
          }
        }
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i<mpi_size; i++)
        send_disp[i] = send_count[i-1]+send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter==0) {
        send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
        send_buffer_int = new int[send_count_sum];
        send_buffer_x = new double[send_count_sum][3];
      }
    }

    delete[] x_fa_approx;

    // now build the recv side stuff...
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    recv_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
    recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
    recv_buffer_int = new int[recv_count_sum];

    MPI_Alltoallv(send_buffer_int, send_count, send_disp, MPI_INT, recv_buffer_int, recv_count, recv_disp, MPI_INT,
        mpi_comm);

    delete[] send_buffer_int;

    for (int i = 0; i<mpi_size; i++) {
      send_count[i] *= 3;
      send_disp[i] *= 3;
      recv_count[i] *= 3;
      recv_disp[i] *= 3;
    }

    recv_buffer_x = new double[recv_count_sum][3];

    MPI_Alltoallv(send_buffer_x, send_count, send_disp, MPI_DOUBLE, recv_buffer_x, recv_count, recv_disp, MPI_DOUBLE,
        mpi_comm);

    delete[] send_buffer_x;

    FOR_ICV {
      x_cv_approx[icv][0] = 0.0;
      x_cv_approx[icv][1] = 0.0;
      x_cv_approx[icv][2] = 0.0;
      cv_flag[icv] = 0; // used for normalization
    }
    for (int i = 0; i<recv_count_sum; ++i) {
      const int icv = recv_buffer_int[i];
      assert((icv>=0)&&(icv<ncv));
      x_cv_approx[icv][0] += recv_buffer_x[i][0];
      x_cv_approx[icv][1] += recv_buffer_x[i][1];
      x_cv_approx[icv][2] += recv_buffer_x[i][2];
      cv_flag[icv] += 1;
    }
    dumpScalarRange(cv_flag, ncv, "CV_FLAG");
    FOR_ICV {
      double tmp = 1.0/(double) cv_flag[icv];
      x_cv_approx[icv][0] *= tmp;
      x_cv_approx[icv][1] *= tmp;
      x_cv_approx[icv][2] *= tmp;
    }

    delete[] recv_buffer_int;
    delete[] recv_buffer_x;
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

    dumpVectorRange(x_cv_approx, ncv, "X_CV_APPROX");

  }

  void calcCvCenterNoInitPC(double(* x_cv_approx)[3])
  {
    if (mpi_rank==0) cout<<"calcCvCenterNoInitPC()"<<endl;

    //Since the grid hasn't been redistributed, we need the calculate approximate CV centers based on approximate face centers

    //Buffers..
    int * send_count = new int[mpi_size];
    int send_count_sum = 0;
    int * send_disp = new int[mpi_size];
    int * send_buffer = NULL;
    int * face_buffer = NULL;

    int mrecv;
    int nof_f;
    int nof_l;

    //Go and get nodal coordinates for each node associated with each face
    for (int iter = 0; iter<2; iter++) {
      if (iter==1) {
        send_count_sum = send_count[mpi_size-1]+send_disp[mpi_size-1];
        send_buffer = new int[send_count_sum];
        face_buffer = new int[send_count_sum];
      }

      send_disp[0] = 0;

      for (int m = 0; m<mpi_size; ++m) {
        send_count[m] = 0;
        FOR_IFA {
          nof_f = noofa_i[ifa];
          nof_l = noofa_i[ifa+1];
          // Loop over all nodes of faces and pack the node number
          for (int nof = nof_f; nof<nof_l; nof++) {
            int ino = noofa_v[nof];
            mrecv = 0;
            while (ino>=noora[mrecv]) {
              mrecv++;
            }
            mrecv--;
            if (mrecv==m) {
              if (iter==1) {
                send_buffer[send_disp[m]+send_count[m]] = ino;
                face_buffer[send_disp[m]+send_count[m]] = ifa; //We keep this for later so that we can fill face centers
              }
              send_count[m]++;
            }
          }
        }

        if (m>0) send_disp[m] = send_count[m-1]+send_disp[m-1];
      }
    }

    //Standard MPI send/receive
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;

    for (int m = 1; m<mpi_size; m++) {
      recv_disp[m] = recv_count[m-1]+recv_disp[m-1];
    }

    int recv_count_sum = recv_count[mpi_size-1]+recv_disp[mpi_size-1];
    int * recv_buffer = new int[recv_count_sum];

    MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

    double * send_buffer_double = new double[3*recv_count_sum];

    //Pack up nodal coordinates in the order we received them
    for (int irecv = 0; irecv<recv_count_sum; irecv++) {
      int ino = recv_buffer[irecv]-noora[mpi_rank];
      send_buffer_double[3*irecv] = x_no[ino][0];
      send_buffer_double[3*irecv+1] = x_no[ino][1];
      send_buffer_double[3*irecv+2] = x_no[ino][2];
    }

    for (int m = 0; m<mpi_size; ++m) {
      send_count[m] = recv_count[m];
      send_disp[m] = recv_disp[m];
    }

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

    recv_disp[0] = 0;

    for (int m = 1; m<mpi_size; m++) {
      recv_disp[m] = recv_count[m-1]+recv_disp[m-1];
    }

    recv_count_sum = recv_count[mpi_size-1]+recv_disp[mpi_size-1];
    double * recv_buffer_double = new double[3*recv_count_sum];

    //Resize the counts and offsets for sending coordinate data
    send_disp[0] = 0;
    recv_disp[0] = 0;

    for (int m = 0; m<mpi_size; ++m) {
      send_count[m] = send_count[m]*3;
      recv_count[m] = recv_count[m]*3;
      if (m>0) {
        send_disp[m] = send_count[m-1]+send_disp[m-1];
        recv_disp[m] = recv_count[m-1]+recv_disp[m-1];
      }
    }

    MPI_Alltoallv(send_buffer_double, send_count, send_disp, MPI_DOUBLE, recv_buffer_double, recv_count, recv_disp,
        MPI_DOUBLE, mpi_comm);

    //Initialize face centers
    double (* x_fa_approx)[3] = new double[nfa][3]; //Face center as geometric average of nodes
    FOR_IFA {
      x_fa_approx[ifa][0] = 0;
      x_fa_approx[ifa][1] = 0;
      x_fa_approx[ifa][2] = 0;
    }

    //Populate face centers using face information generated before
    for (int irecv = 0; irecv<recv_count_sum; irecv++) {
      x_fa_approx[face_buffer[irecv]][0] += recv_buffer_double[3*irecv];
      x_fa_approx[face_buffer[irecv]][1] += recv_buffer_double[3*irecv+1];
      x_fa_approx[face_buffer[irecv]][2] += recv_buffer_double[3*irecv+2];
    }

    //Divide by number of nodes
    FOR_IFA {
      x_fa_approx[ifa][0] /= (noofa_i[ifa+1]-noofa_i[ifa]);
      x_fa_approx[ifa][1] /= (noofa_i[ifa+1]-noofa_i[ifa]);
      x_fa_approx[ifa][2] /= (noofa_i[ifa+1]-noofa_i[ifa]);
    }

    //Pack up x_fa's to send to their corresponding CVs
    for (int iter = 0; iter<2; iter++) {
      if (iter==1) {
        send_count_sum = send_count[mpi_size-1]+send_disp[mpi_size-1];
        send_buffer_double = new double[send_count_sum*3];
        send_buffer = new int[send_count_sum];
      }

      send_disp[0] = 0;

      for (int m = 0; m<mpi_size; ++m) {
        send_count[m] = 0;
        FOR_IFA {
          for (int i = 0; i<2; i++) {
            int icv = cvofa[ifa][i];
            //If the face is periodic, we don't want to send it's center to a CV far away...
            if ((i==1)&&(fa_kind_table[fa_flag[ifa]]!=FA_ZONE_INTERNAL)) {
              icv = -1;
            }
            if (icv!=-1) {
              mrecv = 0;
              while (icv>=cvora[mrecv]) {
                mrecv++;
              }
              mrecv--;
              if (mrecv==m) {
                if (iter==1) {
                  send_buffer[send_disp[m]+send_count[m]] = icv;
                  send_buffer_double[3*(send_disp[m]+send_count[m])] = x_fa_approx[ifa][0];
                  send_buffer_double[3*(send_disp[m]+send_count[m])+1] = x_fa_approx[ifa][1];
                  send_buffer_double[3*(send_disp[m]+send_count[m])+2] = x_fa_approx[ifa][2];
                }
                send_count[m]++;
              }
            }
          }
        }
        if (m>0) send_disp[m] = send_count[m-1]+send_disp[m-1];
      }
    }

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

    recv_disp[0] = 0;

    for (int m = 1; m<mpi_size; m++) {
      recv_disp[m] = recv_count[m-1]+recv_disp[m-1];
    }

    recv_count_sum = recv_count[mpi_size-1]+recv_disp[mpi_size-1];
    recv_buffer = new int[recv_count_sum];

    MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

    recv_buffer_double = new double[3*recv_count_sum];

    //Resize the counts and offsets for sending coordinate data
    send_disp[0] = 0;
    recv_disp[0] = 0;

    for (int m = 0; m<mpi_size; ++m) {
      send_count[m] = send_count[m]*3;
      recv_count[m] = recv_count[m]*3;
      if (m>0) {
        send_disp[m] = send_count[m-1]+send_disp[m-1];
        recv_disp[m] = recv_count[m-1]+recv_disp[m-1];
      }
    }

    MPI_Alltoallv(send_buffer_double, send_count, send_disp, MPI_DOUBLE, recv_buffer_double, recv_count, recv_disp,
        MPI_DOUBLE, mpi_comm);

    //Initialize CV centers
    FOR_ICV {
      x_cv_approx[icv][0] = 0;
      x_cv_approx[icv][1] = 0;
      x_cv_approx[icv][2] = 0;
    }

    int * nfaocv = new int[ncv]; //Number of faces per CV

    FOR_ICV {
      nfaocv[icv] = 0;
    }

    for (int irecv = 0; irecv<recv_count_sum; irecv++) {
      int icv = recv_buffer[irecv]-cvora[mpi_rank];
      x_cv_approx[icv][0] += recv_buffer_double[3*irecv];
      x_cv_approx[icv][1] += recv_buffer_double[3*irecv+1];
      x_cv_approx[icv][2] += recv_buffer_double[3*irecv+2];
      nfaocv[icv]++;
    }

    FOR_ICV {
      x_cv_approx[icv][0] /= nfaocv[icv];
      x_cv_approx[icv][1] /= nfaocv[icv];
      x_cv_approx[icv][2] /= nfaocv[icv];
    }

    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    delete[] send_buffer;
    delete[] recv_buffer;
    delete[] recv_buffer_double;
    delete[] send_buffer_double;
    delete[] x_fa_approx;

  }

  virtual void init(ParamMap& params)
  {
    // init the stats BEFORE reading the restart - it may have stat info.
    initStatsUQ(&params);

    // read the RESTART...
    string restart = params.getStringParam("RESTART");
    readRestart(restart);

    // initialize...
    UgpWithTools::init();

    // init any WRITE_DATA requests...
    initWriteData(&params);

    // calculate the restart file cell centers
    x_cv = new double[ncv][3];
    calcCvCenterInit(x_cv);

    dumpVectorRange(x_no, nno, "X_NO");
    dumpVectorRange(x_cv, ncv, "X_CV");

    // we also need the cvAdt. This builds it if it hasn't been built already...
    useCvAdt();

//    writeData(0);
  }

  virtual double readAndMapRestart(char *fname, ParamMap& params)
  {
    // look for the file "restart.XXXXXX.out". This is the default restart
    // file for the snapshots. If it exists, then read it an initialize/
    // reinitialize the communication pattern...

    // we may already have one in memory. If so, get rid of it...
    if (ugp2!=NULL) delete ugp2;

    // fresh realloc...
    ugp2 = new PostUQ();
    ugp2->registerData(params);
    ugp2->readRestart(fname);
    // NO init: leave the data striped across processors.

    // building the map between the 2 data layouts requires the ugp2 data
    // to have the same face-based centroids calculated: thanks PC...
    ugp2->x_cv = new double[ugp2->ncv][3];
    ugp2->calcCvCenterNoInit(ugp2->x_cv);
    dumpVectorRange(ugp2->x_no, ugp2->nno, "X_NO2");
    dumpVectorRange(ugp2->x_cv, ugp2->ncv, "X_CV2");
    assert(ugp2->cvora[mpi_rank+1]-ugp2->cvora[mpi_rank]==ugp2->ncv);

    // build the map: its structures live here in ugp...
    buildNoInitCvDataMap(ugp2->x_cv, ugp2->cvora);

    // copy over registered data2 that was read in from the snapshot. Following the same convention
    // used in reading full restarts, we know that data was read in if its flag was set...
    for (list<IntValue>::iterator data2 = ugp2->intValueList.begin(); data2!=ugp2->intValueList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"IntValue data2 was found in snapshot: "<<data2->getName()<<endl;
        // this PostUQ should have the same data registered...
        IntValue * data = getIntValueData(data2->getName());
        assert(data!=NULL);
        *(data->ptr) = *(data2->ptr);
      }
    }

    for (list<DoubleValue>::iterator data2 = ugp2->doubleValueList.begin(); data2!=ugp2->doubleValueList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleValue data2 was found in snapshot: "<<data2->getName()<<endl;
        DoubleValue * data = getDoubleValueData(data2->getName());
        assert(data!=NULL);
        *(data->ptr) = *(data2->ptr);
      }
    }

    for (list<DoubleScalar>::iterator data2 = ugp2->doubleScalarList.begin(); data2!=ugp2->doubleScalarList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleScalar data2 was found in snapshot: "<<data2->getName()<<endl;
        DoubleScalar * data = getDoubleScalarData(data2->getName());
        assert(data!=NULL);
        updateNoInitCvData(*(data->ptr), *(data2->ptr)); // note order of args!
      }
    }

    for (list<DoubleVector>::iterator data2 = ugp2->doubleVectorList.begin(); data2!=ugp2->doubleVectorList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleVector data2 was found in snapshot: "<<data2->getName()<<endl;
        DoubleVector * data = getDoubleVectorData(data2->getName());
        assert(data!=NULL);
        updateNoInitCvData(*(data->ptr), *(data2->ptr)); // note order of args!
      }
    }

    return *ugp2->getDoubleValueData("CC_WEIGHT")->ptr;
  }

  virtual void run(ParamMap& params)
  {
    if (mpi_rank==0) cout<<"PostUQ::run()"<<endl;

    // specify the ranges and increments of the restart files
    int ensemble_first = params.getParam("ENSEMBLES")->getInt(1);
    int ensemble_inc   = params.getParam("ENSEMBLES")->getInt(2);
    int ensemble_last  = params.getParam("ENSEMBLES")->getInt(3);

    // get from user the location of the run_ samples
    char home[200];
    sprintf(home, "%s", params.getStringParam("ENSEMBLES_HOME").c_str());

    if (mpi_rank==0)
    {
      cout<<"ENSEMBLES_HOME folder mode: first, inc, last: "<<ensemble_first<<" "<<ensemble_inc<<" "<<ensemble_last<<endl;
      cout<<"from location: "<<home<<endl;
    }

    for (int ensemble = ensemble_first; ensemble<=ensemble_last; ensemble += ensemble_inc)
    {
      // read the restart that might have a different partitioning
      char filename[64];
      sprintf(filename, "%s/case%d/restart.out", home, ensemble);
      if (mpi_rank==0) cout<<"ENSEMBLES_HOME mode: reading restart.out in case"<<ensemble<<"..."<<endl;
      double weight = readAndMapRestart(filename, params);

      if (mpi_rank==0) cout<<"the weight is: " << weight << endl;

      // let the user do something with this new data...
      temporalHook();

      // stats...
      updateStatsUQ(weight);
    }
  }

  virtual void temporalHook()  {  }

  virtual void finalHook(ParamMap& params)
  {
    writeData(0);

    if (params.checkParam("WRITE_RESULT"))      writeRestart("result.out");
  }

};

#endif


// ===============================================================
// ===============================================================
// ===============================================================

class MyPostUQ: public PostUQ {

public:

  MyPostUQ() {
    if (mpi_rank==0) cout<<"PostUQ()"<<endl;

    // you can register additional data here
  }

  virtual double readAndMapRestart(char *fname, ParamMap& params)
  {
    // look for the file "restart.XXXXXX.out". This is the default restart
    // file for the snapshots. If it exists, then read it an initialize/
    // reinitialize the communication pattern...

    // we may already have one in memory. If so, get rid of it...
    if (ugp2!=NULL) delete ugp2;

    // fresh realloc...
    ugp2 = new PostUQ();
    ugp2->registerData(params);
    ugp2->readRestart(fname);
    // NO init: leave the data striped across processors.

    // building the map between the 2 data layouts requires the ugp2 data
    // to have the same face-based centroids calculated: thanks PC...
    ugp2->x_cv = new double[ugp2->ncv][3];
    ugp2->calcCvCenterNoInit(ugp2->x_cv);
    dumpVectorRange(ugp2->x_no, ugp2->nno, "X_NO2");
    dumpVectorRange(ugp2->x_cv, ugp2->ncv, "X_CV2");
    assert(ugp2->cvora[mpi_rank+1]-ugp2->cvora[mpi_rank]==ugp2->ncv);

    // build the map: its structures live here in ugp...
    buildNoInitCvDataMap(ugp2->x_cv, ugp2->cvora);

    // copy over registered data2 that was read in from the snapshot. Following the same convention
    // used in reading full restarts, we know that data was read in if its flag was set...
    for (list<IntValue>::iterator data2 = ugp2->intValueList.begin(); data2!=ugp2->intValueList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"IntValue data2 was found in sample: "<<data2->getName()<<endl;
        // this PostUQ should have the same data registered...
        IntValue * data = getIntValueData(data2->getName());
        assert(data!=NULL);
        *(data->ptr) = *(data2->ptr);
      }
    }

    for (list<DoubleValue>::iterator data2 = ugp2->doubleValueList.begin(); data2!=ugp2->doubleValueList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleValue data2 was found in sample: "<<data2->getName()<<endl;
        DoubleValue * data = getDoubleValueData(data2->getName());
        assert(data!=NULL);
        *(data->ptr) = *(data2->ptr);
      }
    }

    for (list<DoubleScalar>::iterator data2 = ugp2->doubleScalarList.begin(); data2!=ugp2->doubleScalarList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleScalar data2 was found in sample: "<<data2->getName()<<endl;
        DoubleScalar * data = getDoubleScalarData(data2->getName());
        assert(data!=NULL);
        updateNoInitCvData(*(data->ptr), *(data2->ptr)); // note order of args!
      }
    }

    for (list<DoubleVector>::iterator data2 = ugp2->doubleVectorList.begin(); data2!=ugp2->doubleVectorList.end(); data2++) {
      if (data2->getFlag()!=0) {
        if (mpi_rank==0) cout<<"DoubleVector data2 was found in sample: "<<data2->getName()<<endl;
        DoubleVector * data = getDoubleVectorData(data2->getName());
        assert(data!=NULL);
        updateNoInitCvData(*(data->ptr), *(data2->ptr)); // note order of args!
      }
    }

    return *ugp2->getDoubleValueData("CC_WEIGHT")->ptr;
  }

  virtual void run(ParamMap& params)
  {
    if (mpi_rank==0) cout<<"PostUQ::run()"<<endl;

    // specify the ranges and increments of the restart files
    int ensemble_first = params.getParam("ENSEMBLES")->getInt(1);
    int ensemble_inc   = params.getParam("ENSEMBLES")->getInt(2);
    int ensemble_last  = params.getParam("ENSEMBLES")->getInt(3);

    // get from user the location of the run_ samples
    char home[200];
    sprintf(home, "%s", params.getStringParam("ENSEMBLES_HOME").c_str());

    if (mpi_rank==0)
    {
      cout<<"ENSEMBLES_HOME folder mode: first, inc, last: "<<ensemble_first<<" "<<ensemble_inc<<" "<<ensemble_last<<endl;
      cout<<"from location: "<<home<<endl;
    }

    for (int ensemble = ensemble_first; ensemble<=ensemble_last; ensemble += ensemble_inc)
    {
      // read the restart that might have a different partitioning
      char filename[64];
      sprintf(filename, "%s/case%d/restart.out", home, ensemble);
      if (mpi_rank==0) cout<<"ENSEMBLES_HOME mode: reading restart.out in case"<<ensemble<<"..."<<endl;
      double weight = readAndMapRestart(filename, params);

      if (mpi_rank==0) cout<<"the weight is: " << weight << endl;

      // let the user do something with this new data...
      temporalHook();

      // stats...
      updateStatsUQ(weight);
    }


  }
};

// ===========================================================
// ===========================================================
// ===========================================================

int main(int argc, char * argv[])
{
  // intialize the MPI environment - on a single processor,
  // this does nothing.

  MPI_Init(&argc, &argv);

  // try/catch is for catching errors that the code "throws"

  try {

    // this also does nothing in a serial case, in parallel
    // is gets some MPI stuff ready for us to use...

    initMpiStuff();

    ParamMap params;
    params.addParamsFromArgs(argc, argv);
    params.addParamsFromFile("ray.in");

    PostUQ ugp;

    ugp.registerData(params);
    ugp.init(params);
    ugp.run(params);
    ugp.finalHook(params);

  }
  catch (...) {
    cerr<<"unhandled Exception.\n"<<endl;
    MPI_Finalize();
    return (-1);
  }

  // shut down MPI stuff...

  MPI_Barrier( MPI_COMM_WORLD);
  MPI_Finalize();
  return (0);
}



