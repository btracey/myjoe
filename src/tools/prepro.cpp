#include "Param.h"
#include "UgpWithTools.h"
#include "CdpFilter.h"
#include "MshFilter.h"

//#include <algorithm> // required for std::sort on uP

class PreUgp : public UgpWithTools {

private:

  // XXXXX this stuff should be unified with PostUgp some day! 
  // F. Ham, Sept 2009

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

  // no data matching...

  int n_no_pack;
  int * no_pack_count;
  int * no_pack_disp;
  int * no_pack_index;
  int * no_pack_buffer_int;
  double * no_pack_buffer_double;
  
  int n_no_unpack;
  int * no_unpack_count;
  int * no_unpack_disp;
  int * no_unpack_index;
  int * no_unpack_buffer_int;
  double * no_unpack_buffer_double;

protected:
  
  double fa_delta_max;
  
public:

  double (*x_cv)[3];
  
  PreUgp() {
    
    if (mpi_rank == 0)
      cout << "PreUgp()" << endl;
    
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

    n_no_pack = 0;
    no_pack_count = NULL;
    no_pack_disp = NULL;
    no_pack_index = NULL;
    no_pack_buffer_int = NULL;
    no_pack_buffer_double = NULL;
    
    n_no_unpack = 0;
    no_unpack_count = NULL;
    no_unpack_disp = NULL;
    no_unpack_index = NULL;
    no_unpack_buffer_int = NULL;
    no_unpack_buffer_double = NULL;

    // used to match cv-data...
    x_cv = NULL;
    
  }
  
  ~PreUgp() {
    
    // cleanup all the cv_pack and cv_unpack stuff...
    deleteNoInitCvDataMap();
    
    // cleanup all the cv_pack and cv_unpack stuff...
    deleteNoInitNoDataMap();
    
    // x_cv...
    if (x_cv != NULL) {
      delete[] x_cv; 
      x_cv = NULL;
    }
    
  }

  void init() {
    
    UgpWithTools::init();

  }

  void calcCvCenterInit(double (* x_cv_approx)[3]) {
    
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
      for (int nof = nof_f; nof < nof_l; ++nof) {
	int ino = noofa_v[nof];
	x_fa_approx[0] += x_no[ino][0];
	x_fa_approx[1] += x_no[ino][1];
	x_fa_approx[2] += x_no[ino][2];
      }
      //Divide by number of nodes per face
      double tmp = 1.0/(double)(nof_l-nof_f);
      x_fa_approx[0] *= tmp;
      x_fa_approx[1] *= tmp;
      x_fa_approx[2] *= tmp;
      // send to both cv's...
      int icv0 = cvofa[ifa][0];
      assert( (icv0 >= 0)&&(icv0 < ncv) );
      x_cv_approx[icv0][0] += x_fa_approx[0]; 
      x_cv_approx[icv0][1] += x_fa_approx[1]; 
      x_cv_approx[icv0][2] += x_fa_approx[2];
      int icv1 = cvofa[ifa][1];
      if ( (icv1 >= 0)&&(icv1 < ncv) ) {
	x_cv_approx[icv1][0] += x_fa_approx[0]; 
	x_cv_approx[icv1][1] += x_fa_approx[1]; 
	x_cv_approx[icv1][2] += x_fa_approx[2];
      }
    }

    // divide by the number of faces per cell... 
    FOR_ICV {
      double tmp = 1.0/(double)(faocv_i[icv+1]-faocv_i[icv]);
      x_cv_approx[icv][0] *= tmp;
      x_cv_approx[icv][1] *= tmp;
      x_cv_approx[icv][2] *= tmp;
    }

  }
  
  void calcCvCenterNoInit(double (* x_cv_approx)[3]) {
    
    if (mpi_rank == 0) 
      cout << "calcCvCenterNoInit()" << endl;
    
    // we use the noora and cvora here... 
    assert( noora[mpi_rank+1]-noora[mpi_rank] == nno );
    assert( cvora[mpi_rank+1]-cvora[mpi_rank] == ncv );

    // Since the grid hasn't been redistributed, we need the calculate 
    // approximate CV centers based on approximate face centers
    
    // Buffers..
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * send_buffer_int = NULL;
    int send_count_sum;
    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      FOR_IFA {
	const int nof_f = noofa_i[ifa];
	const int nof_l = noofa_i[ifa+1]-1;
	for (int nof = nof_f; nof <= nof_l; ++nof) {
	  // at this point, noofa_v stores the GLOBAL node index...
	  const int ino = noofa_v[nof];
	  assert( (ino >= 0)&&(ino < noora[mpi_size]) );
	  int rank = 0;
	  while (noora[rank+1] <= ino) ++rank;
	  if (iter == 0)
	    send_count[rank] += 1;
	  else {
	    send_buffer_int[send_disp[rank]] = ino - noora[rank]; // local node index
	    send_disp[rank] += 1;
	  }
	}
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	send_disp[i] = send_count[i-1] + send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter == 0) {
	send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
	send_buffer_int = new int[send_count_sum];
      }
    }
    
    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    int * recv_buffer_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buffer_int,send_count,send_disp,MPI_INT,
		  recv_buffer_int,recv_count,recv_disp,MPI_INT,
		  mpi_comm);

    // some cleanup...
    delete[] send_buffer_int;
    
    // now load up the nodal coordinates...
    double (*recv_buffer_x)[3] = new double[recv_count_sum][3];
    for (int i = 0; i < recv_count_sum; ++i) {
      const int ino = recv_buffer_int[i];
      assert( (ino >= 0)&&(ino < nno) );
      recv_buffer_x[i][0] = x_no[ino][0]; 
      recv_buffer_x[i][1] = x_no[ino][1]; 
      recv_buffer_x[i][2] = x_no[ino][2];
    }

    delete[] recv_buffer_int;
    
    // and send back: here we have three doubles per previous int...
    for (int i = 0; i < mpi_size; i++) {
      send_count[i] *= 3;
      send_disp[i]  *= 3;
      recv_count[i] *= 3;
      recv_disp[i]  *= 3;
    }
    double (*send_buffer_x)[3] = new double[send_count_sum][3];
    MPI_Alltoallv(recv_buffer_x,recv_count,recv_disp,MPI_DOUBLE,
		  send_buffer_x,send_count,send_disp,MPI_DOUBLE,
		  mpi_comm);
    
    // some more cleanup...
    delete[] recv_buffer_x;
    
    // return send_disp... 
    for (int i = 0; i < mpi_size; i++)
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
      for (int nof = nof_f; nof <= nof_l; ++nof) {
	// at this point, noofa_v stores the GLOBAL node index...
	const int ino = noofa_v[nof];
	assert( (ino >= 0)&&(ino < noora[mpi_size]) );
	int rank = 0;
	while (noora[rank+1] <= ino) ++rank;
	x_fa_approx[ifa][0] += send_buffer_x[send_disp[rank]][0]; 
	x_fa_approx[ifa][1] += send_buffer_x[send_disp[rank]][1]; 
	x_fa_approx[ifa][2] += send_buffer_x[send_disp[rank]][2];
	send_disp[rank] += 1;
      }
      double tmp = 1.0/(double)(nof_l-nof_f+1);
      x_fa_approx[ifa][0] *= tmp;
      x_fa_approx[ifa][1] *= tmp;
      x_fa_approx[ifa][2] *= tmp;
    }
    
    // more cleanup...
    delete[] send_buffer_x;
    
    // ===========================
    // now reduce to cv's...
    // ===========================

    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      FOR_IFA {
	// icv0...
	const int icv0 = cvofa[ifa][0];
	assert( (icv0 >= 0)&&(icv0 < cvora[mpi_size]) );
	int rank = 0;
	while (cvora[rank+1] <= icv0) ++rank;
	if (iter == 0)
	  send_count[rank] += 1;
	else {
	  send_buffer_int[send_disp[rank]] = icv0 - cvora[rank]; // local cv index
	  send_buffer_x[send_disp[rank]][0] = x_fa_approx[ifa][0];
	  send_buffer_x[send_disp[rank]][1] = x_fa_approx[ifa][1];
	  send_buffer_x[send_disp[rank]][2] = x_fa_approx[ifa][2];
	  send_disp[rank] += 1;
	}
	// icv1...
	// the fa_kind_table should be built at teh end of readRestart...
        if (fa_kind_table[fa_flag[ifa]] == FA_ZONE_INTERNAL) {
	  const int icv1 = cvofa[ifa][1];
	  assert( (icv1 >= 0)&&(icv1 < cvora[mpi_size]) );
	  rank = 0;
	  while (cvora[rank+1] <= icv1) ++rank;
	  if (iter == 0)
	    send_count[rank] += 1;
	  else {
	    send_buffer_int[send_disp[rank]] = icv1 - cvora[rank]; // local cv index
	    send_buffer_x[send_disp[rank]][0] = x_fa_approx[ifa][0];
	    send_buffer_x[send_disp[rank]][1] = x_fa_approx[ifa][1];
	    send_buffer_x[send_disp[rank]][2] = x_fa_approx[ifa][2];
	    send_disp[rank] += 1;
	  }
	}
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	send_disp[i] = send_count[i-1] + send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter == 0) {
	send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
	send_buffer_int = new int[send_count_sum];
	send_buffer_x = new double[send_count_sum][3];
      }
    }
    
    delete[] x_fa_approx;
    
    // now build the recv side stuff...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    recv_buffer_int = new int[recv_count_sum];

    MPI_Alltoallv(send_buffer_int,send_count,send_disp,MPI_INT,
		  recv_buffer_int,recv_count,recv_disp,MPI_INT,
		  mpi_comm);

    delete[] send_buffer_int;

    for (int i = 0; i < mpi_size; i++) {
      send_count[i] *= 3;
      send_disp[i]  *= 3;
      recv_count[i] *= 3;
      recv_disp[i]  *= 3;
    }
    
    recv_buffer_x = new double[recv_count_sum][3];
    
    MPI_Alltoallv(send_buffer_x,send_count,send_disp,MPI_DOUBLE,
		  recv_buffer_x,recv_count,recv_disp,MPI_DOUBLE,
		  mpi_comm);

    delete[] send_buffer_x;

    FOR_ICV {
      x_cv_approx[icv][0] = 0.0;
      x_cv_approx[icv][1] = 0.0;
      x_cv_approx[icv][2] = 0.0;
      cv_flag[icv] = 0; // used for normalization
    }
    for (int i = 0; i < recv_count_sum; ++i) {
      const int icv = recv_buffer_int[i];
      assert( (icv >= 0)&&(icv < ncv) );
      x_cv_approx[icv][0] += recv_buffer_x[i][0];
      x_cv_approx[icv][1] += recv_buffer_x[i][1];
      x_cv_approx[icv][2] += recv_buffer_x[i][2];
      cv_flag[icv] += 1;
    }
    dumpScalarRange(cv_flag,ncv,"CV_FLAG");
    FOR_ICV {
      double tmp = 1.0/(double)cv_flag[icv];
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
    
    dumpVectorRange(x_cv_approx,ncv,"X_CV_APPROX");

  }

  void buildNoInitCvDataMap2(double (*x_cv2)[3],int * cvora2) {

    // somewhat different than the cv data map built in Postpro for
    // 1:1 snapshot interpolation, this one allows for non-coincident data
    // and returns the cv's that are unset in cv_flag...

    if (mpi_rank == 0)
      cout << "buildNoInitCvDataMap2()" << endl;

    // cleanup from any previous map...
    deleteNoInitCvDataMap();

    //assert( cvora[mpi_size] == cvora2[mpi_size] );
    int ncv2 = cvora2[mpi_rank+1] - cvora2[mpi_rank];

    // prepare to send the x_cv2's to ranks where there may be 
    // a match. The send side is considered the data layout associated
    // with x_cv2, ncv2, cvora2...
    
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    double * send_buffer_x_cv2 = NULL;
    int send_count_sum;
    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int icv2 = 0; icv2 < ncv2; ++icv2) {
	for (int rank = 0; rank < mpi_size; ++rank) {
	  if ( (x_cv2[icv2][0] >= cvAdtBBMinOfRank[rank][0])&&(x_cv2[icv2][0] <= cvAdtBBMaxOfRank[rank][0])&&
	       (x_cv2[icv2][1] >= cvAdtBBMinOfRank[rank][1])&&(x_cv2[icv2][1] <= cvAdtBBMaxOfRank[rank][1])&&
	       (x_cv2[icv2][2] >= cvAdtBBMinOfRank[rank][2])&&(x_cv2[icv2][2] <= cvAdtBBMaxOfRank[rank][2]) ) {
	    if (iter == 0) {
	      send_count[rank] += 3;
	    }
	    else {
	      send_buffer_x_cv2[send_disp[rank]  ] = x_cv2[icv2][0];
	      send_buffer_x_cv2[send_disp[rank]+1] = x_cv2[icv2][1];
	      send_buffer_x_cv2[send_disp[rank]+2] = x_cv2[icv2][2];
	      send_disp[rank] += 3;
	    }
	  }
	}
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	send_disp[i] = send_count[i-1] + send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter == 0) {
	send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
	send_buffer_x_cv2 = new double[send_count_sum];
      }
    }

    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    double * recv_buffer_x_cv2 = new double[recv_count_sum];
    MPI_Alltoallv(send_buffer_x_cv2,send_count,send_disp,MPI_DOUBLE,
		  recv_buffer_x_cv2,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    
    // the x_cv2 data is now on the recv side. Use the cvAdt to identify point matches
    // and return the closest match we can find. since we are now going to send one
    // value per request rather than 3, remove the factor of 3...
    send_count_sum /= 3;
    recv_count_sum /= 3;
    for (int i = 0; i < mpi_size; i++) {
      send_count[i] /= 3;
      send_disp[i] /= 3;
      recv_count[i] /= 3;
      recv_disp[i] /= 3;
    }

    delete[] send_buffer_x_cv2;
    
    double * recv_buffer_d2 = new double[recv_count_sum];
    int * recv_buffer_icv = new int[recv_count_sum];
    
    int bbListSize,bbList[ADT_LIST_MAX];
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      double this_x_cv2[3];
      FOR_I3 this_x_cv2[i] = recv_buffer_x_cv2[3*irecv+i];
      recv_buffer_icv[irecv] = -1;
      // use the Adt to get a list of possible candidates...
      cvAdt->buildListForPoint(bbListSize,bbList,this_x_cv2);
      // cycle through candidates and store the closest...
      for (int ibb = 0; ibb < bbListSize; ++ibb) {
	const int icv = bbList[ibb];
	double d2 = 0.0;
	FOR_I3 {
	  double dx = x_cv[icv][i] - this_x_cv2[i];
	  d2 += dx*dx;
	}
	if ( (recv_buffer_icv[irecv] == -1)||(d2 < recv_buffer_d2[irecv]) ) {
	  recv_buffer_icv[irecv] = icv; // the local icv is fine
	  recv_buffer_d2[irecv] = d2;
	}
      }
    }
    
    delete[] recv_buffer_x_cv2;
    
    double * send_buffer_d2 = new double[send_count_sum];
    int * send_buffer_icv = new int[send_count_sum];

    // send the recv stuff nack to the send data layout...
    MPI_Alltoallv(recv_buffer_icv,recv_count,recv_disp,MPI_INT,
		  send_buffer_icv,send_count,send_disp,MPI_INT,
		  mpi_comm);
    MPI_Alltoallv(recv_buffer_d2,recv_count,recv_disp,MPI_DOUBLE,
		  send_buffer_d2,send_count,send_disp,MPI_DOUBLE,
		  mpi_comm);
    
    // HERE - do a double loop, although there is likely some shorter way
    // to do this because we know that the data layouts are 1:1 matched. Other topologies
    // are not single-valued (nodes, edges, faces), so this is more general...
    assert( cv_pack_count == NULL );
    cv_pack_count = new int[mpi_size];
    assert( cv_pack_disp == NULL );
    cv_pack_disp = new int[mpi_size];
    // these get set below...
    assert( cv_pack_index == NULL );
    assert( cv_pack_buffer_int == NULL );
    for (int i = 0; i < mpi_size; i++)
      cv_pack_count[i] = 0;
    double my_d2_max = 0.0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icv2 = 0; icv2 < ncv2; ++icv2) {
	int rank_min = -1;
	int icv_min;
	double d2_min;
	for (int rank = 0; rank < mpi_size; ++rank) {
	  if ( (x_cv2[icv2][0] >= cvAdtBBMinOfRank[rank][0])&&(x_cv2[icv2][0] <= cvAdtBBMaxOfRank[rank][0])&&
	       (x_cv2[icv2][1] >= cvAdtBBMinOfRank[rank][1])&&(x_cv2[icv2][1] <= cvAdtBBMaxOfRank[rank][1])&&
	       (x_cv2[icv2][2] >= cvAdtBBMinOfRank[rank][2])&&(x_cv2[icv2][2] <= cvAdtBBMaxOfRank[rank][2]) ) {
	    double this_d2 = send_buffer_d2[send_disp[rank]];
	    int this_icv = send_buffer_icv[send_disp[rank]];
	    send_disp[rank] += 1;
	    if ( ((rank_min == -1)||(this_d2 < d2_min)) && (this_icv >= 0) ) {
	      rank_min = rank;
	      icv_min  = this_icv;
	      d2_min   = this_d2;
	    }
	  }
	}
	// for non-coincident data, we may not find even 1...
	if ( rank_min >= 0 ) {
	  // we need to build a pack/unpack buffer that knows how to
	  // pack the data quickly...
	  if (iter == 0)
	    cv_pack_count[rank_min] += 1;
	  else {
	    // store the max d2 (of all our min's) for reporting...
	    my_d2_max = max(my_d2_max,d2_min);
	    // fill index and buffer, inc disp... 
	    cv_pack_index[cv_pack_disp[rank_min]] = icv2; // local value
	    cv_pack_buffer_int[cv_pack_disp[rank_min]] = icv_min; // send to rank_min below
	    cv_pack_disp[rank_min] += 1;
	  }
	}
      }
      // calculate cv_pack_disp on both iterations...
      cv_pack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	cv_pack_disp[i] = cv_pack_count[i-1] + cv_pack_disp[i-1];
      // on the first time through...
      if (iter == 0) {
	// rebuild send_disp...
	send_disp[0] = 0;
	for (int i = 1; i < mpi_size; i++)
	  send_disp[i] = send_count[i-1] + send_disp[i-1];
	// allocate pack buffers...
	n_cv_pack = cv_pack_disp[mpi_size-1] + cv_pack_count[mpi_size-1];
	cv_pack_index = new int[n_cv_pack];
	cv_pack_buffer_int = new int[n_cv_pack];
      }
    }
    
    double d2_max;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);    
    if (mpi_rank == 0)
      cout << " > cvs matched within tol: " << sqrt(d2_max) << endl;
    
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
    assert( cv_unpack_count == NULL );
    cv_unpack_count = new int[mpi_size];
    MPI_Alltoall(cv_pack_count,1,MPI_INT,cv_unpack_count,1,MPI_INT,mpi_comm);
    assert( cv_unpack_disp == NULL );
    cv_unpack_disp = new int[mpi_size];
    cv_unpack_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      cv_unpack_disp[i] = cv_unpack_count[i-1] + cv_unpack_disp[i-1];
    n_cv_unpack = cv_unpack_disp[mpi_size-1] + cv_unpack_count[mpi_size-1];
    assert(cv_unpack_index == NULL );
    cv_unpack_index = new int[n_cv_unpack];
    // exchange directly into the unpack index...
    MPI_Alltoallv(cv_pack_buffer_int,cv_pack_count,cv_pack_disp,MPI_INT,
		  cv_unpack_index,cv_unpack_count,cv_unpack_disp,MPI_INT,mpi_comm);

    // allocate remaining buffers...
    assert( cv_unpack_buffer_int == NULL );
    cv_unpack_buffer_int = new int[n_cv_unpack];
    assert( cv_pack_buffer_double == NULL );
    cv_pack_buffer_double = new double[n_cv_pack];
    assert( cv_unpack_buffer_double == NULL );
    cv_unpack_buffer_double = new double[n_cv_unpack];
    
    // check: try comparing components of x_cv2...
    /*
      FOR_I3 {
      for (int ii = 0; ii < n_cv_pack; ++ii) {
      int icv2 = cv_pack_index[ii];
      cv_pack_buffer_double[ii] = x_cv2[icv2][i];
      }
      MPI_Alltoallv(cv_pack_buffer_double,cv_pack_count,cv_pack_disp,MPI_DOUBLE,
      cv_unpack_buffer_double,cv_unpack_count,cv_unpack_disp,MPI_DOUBLE,mpi_comm);
      double my_d2_max = 0.0;
      for (int ii = 0; ii < n_cv_unpack; ++ii) {
      int icv = cv_unpack_index[ii];
      double dx = cv_unpack_buffer_double[ii] - x_cv[icv][i];
      my_d2_max = max(my_d2_max,dx*dx);
      }
      double d2_max;
      MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);    
      if (mpi_rank == 0)
      cout << " > check for component: " << i << " tol (should be zero): " << sqrt(d2_max) << endl;
      }
    */

    // also check how many cv's are in fact touched and touched
    // once during uppacking...
    FOR_ICV cv_flag[icv] = 0;
    for (int ii = 0; ii < n_cv_unpack; ++ii) {
      int icv = cv_unpack_index[ii];
      cv_flag[icv] += 1;
    }

    int my_buf[2] = { ncv , 0 };
    FOR_ICV 
      if (cv_flag[icv] > 0)
	my_buf[1] += 1;
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);    
    if (mpi_rank == 0)
      cout << " > buildNoInitCvDataMap2: set " << buf[1] << " out of " << buf[0] << " cvs" << endl;
  
  }
  
  void deleteNoInitCvDataMap() {
    n_cv_pack = 0;
    if (cv_pack_count != NULL) {
      delete[] cv_pack_count; cv_pack_count = NULL;
    }
    if (cv_pack_disp != NULL) {
      delete[] cv_pack_disp; cv_pack_disp = NULL;
    }
    if (cv_pack_index != NULL) {
      delete[] cv_pack_index; cv_pack_index = NULL;
    }
    if (cv_pack_buffer_int != NULL) {
      delete[] cv_pack_buffer_int; cv_pack_buffer_int = NULL;
    }
    if (cv_pack_buffer_double != NULL) {
      delete[] cv_pack_buffer_double; cv_pack_buffer_double = NULL;
    }
    
    n_cv_unpack = 0;
    if (cv_unpack_count != NULL) {
      delete[] cv_unpack_count; cv_unpack_count = NULL;
    }
    if (cv_unpack_disp != NULL) {
      delete[] cv_unpack_disp; cv_unpack_disp = NULL;
    }
    if (cv_unpack_index != NULL) {
      delete[] cv_unpack_index; cv_unpack_index = NULL;
    }
    if (cv_unpack_buffer_int != NULL) {
      delete[] cv_unpack_buffer_int; cv_unpack_buffer_int = NULL;
    }
    if (cv_unpack_buffer_double != NULL) {
      delete[] cv_unpack_buffer_double; cv_unpack_buffer_double = NULL;
    }
  }

  void updateNoInitCvData(double * scalar,const double * noInitScalar) {

    for (int ii = 0; ii < n_cv_pack; ++ii) {
      int icv2 = cv_pack_index[ii];
      cv_pack_buffer_double[ii] = noInitScalar[icv2];
    }
    
    MPI_Alltoallv(cv_pack_buffer_double,  cv_pack_count,  cv_pack_disp,  MPI_DOUBLE,
		  cv_unpack_buffer_double,cv_unpack_count,cv_unpack_disp,MPI_DOUBLE,mpi_comm);
    
    for (int ii = 0; ii < n_cv_unpack; ++ii) {
      int icv = cv_unpack_index[ii];
      scalar[icv] = cv_unpack_buffer_double[ii];
    }

  }
  
  void updateNoInitCvData(double (*vector)[3],const double (*noInitVector)[3]) {

    FOR_I3 {
      for (int ii = 0; ii < n_cv_pack; ++ii) {
	int icv2 = cv_pack_index[ii];
	cv_pack_buffer_double[ii] = noInitVector[icv2][i];
      }
      
      MPI_Alltoallv(cv_pack_buffer_double,  cv_pack_count,  cv_pack_disp,  MPI_DOUBLE,
		    cv_unpack_buffer_double,cv_unpack_count,cv_unpack_disp,MPI_DOUBLE,mpi_comm);
    
      for (int ii = 0; ii < n_cv_unpack; ++ii) {
	int icv = cv_unpack_index[ii];
	vector[icv][i] = cv_unpack_buffer_double[ii];
      }
    }

  }


  void buildNoInitNoDataMap(double (*x_no2)[3],int * noora2) {
    
    if (mpi_rank == 0)
      cout << "buildNoInitNoDataMap()" << endl;
    
    // cleanup from any previous map...
    deleteNoInitNoDataMap();
    
    // for node matching, we send the node value from the snapshot to
    // the first matched location, then do a reduction to get the correct
    // value everywhere...
    int nno2 = noora2[mpi_rank+1] - noora2[mpi_rank];

    // prepare to send the x_no2's to ranks where there may be 
    // a match. The send side is considered the data layout associated
    // with x_no2, nno2, noora2...
    
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    double * send_buffer_x_no2 = NULL;
    int send_count_sum;
    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int ino2 = 0; ino2 < nno2; ++ino2) {
	for (int rank = 0; rank < mpi_size; ++rank) {
	  if ( (x_no2[ino2][0] >= cvAdtBBMinOfRank[rank][0])&&(x_no2[ino2][0] <= cvAdtBBMaxOfRank[rank][0])&&
	       (x_no2[ino2][1] >= cvAdtBBMinOfRank[rank][1])&&(x_no2[ino2][1] <= cvAdtBBMaxOfRank[rank][1])&&
	       (x_no2[ino2][2] >= cvAdtBBMinOfRank[rank][2])&&(x_no2[ino2][2] <= cvAdtBBMaxOfRank[rank][2]) ) {
	    if (iter == 0) {
	      send_count[rank] += 3;
	    }
	    else {
	      send_buffer_x_no2[send_disp[rank]  ] = x_no2[ino2][0];
	      send_buffer_x_no2[send_disp[rank]+1] = x_no2[ino2][1];
	      send_buffer_x_no2[send_disp[rank]+2] = x_no2[ino2][2];
	      send_disp[rank] += 3;
	    }
	  }
	}
      }
      // calculate send_disp on both iterations...
      send_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	send_disp[i] = send_count[i-1] + send_disp[i-1];
      // on the first time through, allocate the send_buffer...
      if (iter == 0) {
	send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
	send_buffer_x_no2 = new double[send_count_sum];
      }
    }
    
    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    double * recv_buffer_x_no2 = new double[recv_count_sum];
    MPI_Alltoallv(send_buffer_x_no2,send_count,send_disp,MPI_DOUBLE,
		  recv_buffer_x_no2,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    
    // the x_no2 data is now on the recv side. Use the cvAdt to identify point matches
    // and return the closest match we can find. since we are now going to send one
    // value per request rather than 3, remove the factor of 3...
    send_count_sum /= 3;
    recv_count_sum /= 3;
    for (int i = 0; i < mpi_size; i++) {
      send_count[i] /= 3;
      send_disp[i] /= 3;
      recv_count[i] /= 3;
      recv_disp[i] /= 3;
    }

    delete[] send_buffer_x_no2;
    
    double * recv_buffer_d2 = new double[recv_count_sum];
    int * recv_buffer_ino = new int[recv_count_sum];
    
    int bbListSize,bbList[ADT_LIST_MAX];
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      double this_x_no2[3];
      FOR_I3 this_x_no2[i] = recv_buffer_x_no2[3*irecv+i];
      recv_buffer_ino[irecv] = -1;
      // use the Adt to get a list of possible cv candidates...
      cvAdt->buildListForPoint(bbListSize,bbList,this_x_no2);
      // cycle through candidate cvs and store the closest node...
      for (int ibb = 0; ibb < bbListSize; ++ibb) {
	const int icv = bbList[ibb];
	for (int noc = noocv_i[icv]; noc < noocv_i[icv+1]; ++noc) {
	  int ino = noocv_v[noc];
	  double d2 = 0.0;
	  FOR_I3 {
	    double dx = x_no[ino][i] - this_x_no2[i];
	    d2 += dx*dx;
	  }
	  if ( (recv_buffer_ino[irecv] == -1)||(d2 < recv_buffer_d2[irecv]) ) {
	    recv_buffer_ino[irecv] = ino; // the local ino is fine
	    recv_buffer_d2[irecv] = d2;
	  }
	}
      }
    }
    
    delete[] recv_buffer_x_no2;
    
    double * send_buffer_d2 = new double[send_count_sum];
    int * send_buffer_ino = new int[send_count_sum];
    
    // send the recv stuff nack to the send data layout...
    MPI_Alltoallv(recv_buffer_ino,recv_count,recv_disp,MPI_INT,
		  send_buffer_ino,send_count,send_disp,MPI_INT,
		  mpi_comm);
    MPI_Alltoallv(recv_buffer_d2,recv_count,recv_disp,MPI_DOUBLE,
		  send_buffer_d2,send_count,send_disp,MPI_DOUBLE,
		  mpi_comm);
    
    assert( no_pack_count == NULL );
    no_pack_count = new int[mpi_size];
    assert( no_pack_disp == NULL );
    no_pack_disp = new int[mpi_size];
    // these get set below...
    assert( no_pack_index == NULL );
    assert( no_pack_buffer_int == NULL );
    for (int i = 0; i < mpi_size; i++)
      no_pack_count[i] = 0;
    double my_d2_max = 0.0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ino2 = 0; ino2 < nno2; ++ino2) {
	int rank_min = -1;
	int ino_min;
	double d2_min;
	for (int rank = 0; rank < mpi_size; ++rank) {
	  if ( (x_no2[ino2][0] >= cvAdtBBMinOfRank[rank][0])&&(x_no2[ino2][0] <= cvAdtBBMaxOfRank[rank][0])&&
	       (x_no2[ino2][1] >= cvAdtBBMinOfRank[rank][1])&&(x_no2[ino2][1] <= cvAdtBBMaxOfRank[rank][1])&&
	       (x_no2[ino2][2] >= cvAdtBBMinOfRank[rank][2])&&(x_no2[ino2][2] <= cvAdtBBMaxOfRank[rank][2]) ) {
	    double this_d2 = send_buffer_d2[send_disp[rank]];
	    int this_ino = send_buffer_ino[send_disp[rank]];
	    send_disp[rank] += 1;
	    if ( ((rank_min == -1)||(this_d2 < d2_min)) && (this_ino >= 0) ) {
	      rank_min = rank;
	      ino_min  = this_ino;
	      d2_min   = this_d2;
	    }
	  }
	}
	// if we don't find even one, then this is a problem...
	assert( rank_min >= 0 );
	// we need to build a pack/unpack buffer that knows how to
	// pack the data quickly...
	if (iter == 0)
	  no_pack_count[rank_min] += 1;
	else {
	  // store the max d2 (of all our min's) for reporting...
	  my_d2_max = max(my_d2_max,d2_min);
	  // fill index and buffer, inc disp... 
	  no_pack_index[no_pack_disp[rank_min]] = ino2; // local value
	  no_pack_buffer_int[no_pack_disp[rank_min]] = ino_min; // send to rank_min below
	  no_pack_disp[rank_min] += 1;
	}
      }
      // calculate no_pack_disp on both iterations...
      no_pack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
	no_pack_disp[i] = no_pack_count[i-1] + no_pack_disp[i-1];
      // on the first time through...
      if (iter == 0) {
	// rebuild send_disp...
	send_disp[0] = 0;
	for (int i = 1; i < mpi_size; i++)
	  send_disp[i] = send_count[i-1] + send_disp[i-1];
	// allocate pack buffers...
	n_no_pack = no_pack_disp[mpi_size-1] + no_pack_count[mpi_size-1];
	no_pack_index = new int[n_no_pack];
	no_pack_buffer_int = new int[n_no_pack];
      }
    }

    double d2_max;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);    
    if (mpi_rank == 0)
      cout << " > all nos matched within tol: " << sqrt(d2_max) << endl;
    
    // cleanup...
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    delete[] send_buffer_d2;
    delete[] send_buffer_ino;
    delete[] recv_buffer_d2;
    delete[] recv_buffer_ino;

    // now build the unpack stuff...
    assert( no_unpack_count == NULL );
    no_unpack_count = new int[mpi_size];
    MPI_Alltoall(no_pack_count,1,MPI_INT,no_unpack_count,1,MPI_INT,mpi_comm);
    assert( no_unpack_disp == NULL );
    no_unpack_disp = new int[mpi_size];
    no_unpack_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      no_unpack_disp[i] = no_unpack_count[i-1] + no_unpack_disp[i-1];
    n_no_unpack = no_unpack_disp[mpi_size-1] + no_unpack_count[mpi_size-1];
    assert(no_unpack_index == NULL );
    no_unpack_index = new int[n_no_unpack];
    // exchange directly into the unpack index...
    MPI_Alltoallv(no_pack_buffer_int,no_pack_count,no_pack_disp,MPI_INT,
		  no_unpack_index,no_unpack_count,no_unpack_disp,MPI_INT,mpi_comm);

    // allocate remaining buffers...
    assert( no_unpack_buffer_int == NULL );
    no_unpack_buffer_int = new int[n_no_unpack];
    assert( no_pack_buffer_double == NULL );
    no_pack_buffer_double = new double[n_no_pack];
    assert( no_unpack_buffer_double == NULL );
    no_unpack_buffer_double = new double[n_no_unpack];
    
    // check: try comparing components of x_no2...
    FOR_I3 {
      for (int ii = 0; ii < n_no_pack; ++ii) {
	int ino2 = no_pack_index[ii];
	no_pack_buffer_double[ii] = x_no2[ino2][i];
      }
      MPI_Alltoallv(no_pack_buffer_double,no_pack_count,no_pack_disp,MPI_DOUBLE,
		    no_unpack_buffer_double,no_unpack_count,no_unpack_disp,MPI_DOUBLE,mpi_comm);
      double my_d2_max = 0.0;
      for (int ii = 0; ii < n_no_unpack; ++ii) {
	int ino = no_unpack_index[ii];
	double dx = no_unpack_buffer_double[ii] - x_no[ino][i];
	my_d2_max = max(my_d2_max,dx*dx);
      }
      double d2_max;
      MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);    
      if (mpi_rank == 0)
	cout << " > check for component: " << i << " tol (should be zero): " << sqrt(d2_max) << endl;
    }

    // also check that all no's are in fact touched and touched
    // once during unpacking (for point-matched grids, otherwise
    // comment out the assert below for general node interpolation between
    // different grids)...
    FOR_INO no_flag[ino] = 0;
    for (int ii = 0; ii < n_no_unpack; ++ii) {
      int ino = no_unpack_index[ii];
      //assert( no_flag[ino] == 0 );
      no_flag[ino] = 1;
    }

    /*
      int my_buf[2] = { nno , 0 };
    FOR_INO
      if (no_flag[ino] > 0)
	my_buf[1] += 1;
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);    
    if (mpi_rank == 0)
      cout << " > buildNoInitCvDataMap2: set " << buf[1] << " out of " << buf[0] << " cvs" << endl;
    updateNoData(no_flag,ADD_NO_PERIODIC_DATA);
    FOR_INO assert( no_flag[ino] == 1 );
    */

  }
  
  void deleteNoInitNoDataMap() {
    n_no_pack = 0;
    if (no_pack_count != NULL) {
      delete[] no_pack_count; no_pack_count = NULL;
    }
    if (no_pack_disp != NULL) {
      delete[] no_pack_disp; no_pack_disp = NULL;
    }
    if (no_pack_index != NULL) {
      delete[] no_pack_index; no_pack_index = NULL;
    }
    if (no_pack_buffer_int != NULL) {
      delete[] no_pack_buffer_int; no_pack_buffer_int = NULL;
    }
    if (no_pack_buffer_double != NULL) {
      delete[] no_pack_buffer_double; no_pack_buffer_double = NULL;
    }
    
    n_no_unpack = 0;
    if (no_unpack_count != NULL) {
      delete[] no_unpack_count; no_unpack_count = NULL;
    }
    if (no_unpack_disp != NULL) {
      delete[] no_unpack_disp; no_unpack_disp = NULL;
    }
    if (no_unpack_index != NULL) {
      delete[] no_unpack_index; no_unpack_index = NULL;
    }
    if (no_unpack_buffer_int != NULL) {
      delete[] no_unpack_buffer_int; no_unpack_buffer_int = NULL;
    }
    if (no_unpack_buffer_double != NULL) {
      delete[] no_unpack_buffer_double; no_unpack_buffer_double = NULL;
    }
  }
  
  void updateNoInitNoData(double * scalar,const double * noInitScalar) {

    FOR_INO scalar[ino] = 0.0;

    for (int ii = 0; ii < n_no_pack; ++ii) {
      int ino2 = no_pack_index[ii];
      no_pack_buffer_double[ii] = noInitScalar[ino2];
    }
    
    MPI_Alltoallv(no_pack_buffer_double,  no_pack_count,  no_pack_disp,  MPI_DOUBLE,
		  no_unpack_buffer_double,no_unpack_count,no_unpack_disp,MPI_DOUBLE, mpi_comm);
    
    for (int ii = 0; ii < n_no_unpack; ++ii) {
      int ino = no_unpack_index[ii];
      scalar[ino] = no_unpack_buffer_double[ii];
    }

    updateNoData(scalar,ADD_NO_PERIODIC_DATA);

  }
  
  void updateNoInitNoData(double (*vector)[3],const double (*noInitVector)[3]) {

    FOR_INO FOR_I3 vector[ino][i] = 0.0;

    FOR_I3 {
      for (int ii = 0; ii < n_no_pack; ++ii) {
	int ino2 = no_pack_index[ii];
	no_pack_buffer_double[ii] = noInitVector[ino2][i];
      }
      
      MPI_Alltoallv(no_pack_buffer_double,  no_pack_count,  no_pack_disp,  MPI_DOUBLE,
		    no_unpack_buffer_double,no_unpack_count,no_unpack_disp,MPI_DOUBLE,mpi_comm);
    
      for (int ii = 0; ii < n_no_unpack; ++ii) {
	int ino = no_unpack_index[ii];
	vector[ino][i] = no_unpack_buffer_double[ii];
      }
    }

    updateNoData(vector,ADD_NO_PERIODIC_DATA);

  }

  void setFaFlagForCvGroupHook() {

    if (mpi_rank == 0)
      cout << "PreUgp::setFaFlagForCvGroupHook(), fa_delta_max: " << fa_delta_max << endl;
    
    // this hook is called by the Ugp infrastructure to 
    // setup the cvGroup[icv] integer. Basically set the fa_flag
    // to 1 between cells belonging in the same group, or 0
    // if they are potentially different groups...
    
    double * fa_delta = new double[nfa];
    FOR_IFA fa_delta[ifa] = 0.0;
    
    FOR_ICV {
      
      // we need the "approximate" cv center...
      double x_cv[3];
      calcCvCenter(x_cv,icv);
      
      // now loop on the faces and compute the distance between cv's 
      // projected on the face normal...
      const int foc_f = faocv_i[icv]; 
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	double x_fa[3],n_fa[3];
	calcFaceCenterAndUnitNormal(x_fa,n_fa,ifa);
	double delta = 0.0;
	if (cvofa[ifa][0] == icv) {
	  FOR_I3 delta += (x_fa[i] - x_cv[i])*n_fa[i];
	}
	else {
	  assert( cvofa[ifa][1] == icv );
	  FOR_I3 delta -= (x_fa[i] - x_cv[i])*n_fa[i];
	}
	assert( delta > 0.0 );
	fa_delta[ifa] += delta;
      }

    }
      
    // add across inter-processor/periodic faces...
    updateFaData(fa_delta,ADD_DATA);
    dumpScalarRange(fa_delta,nfa,"FA_DELTA");
    
    FOR_IFA
      if (fa_delta[ifa] <= fa_delta_max)
	fa_flag[ifa] = 1;
      else 
	fa_flag[ifa] = 0;
    
    delete[] fa_delta;

    // as a hack, flip the faces and have a look...
    FOR_IFA fa_flag[ifa] = 1-fa_flag[ifa];
    writeFlaggedFacesTecplot("explicit_faces.dat");
    FOR_IFA fa_flag[ifa] = 1-fa_flag[ifa];

  }
  
  void extendRegisteredNoData() {
    
    // for no data, reduce along edges - 
    
    int * no_count = new int[nno];
    
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++) {
      if (data->getDatatype() == NO_DATA) {
	
	FOR_INO {
	  no_count[ino] = no_flag[ino];
	  // expect zero or 1...
	  assert( (no_count[ino] == 0)||(no_count[ino] == 1) );
	}
	
	int done = 0;
	int iter = 0;
	while (done != 1) {
	  
	  int my_done = 1;
	  
	  iter++;
	  if (mpi_rank == 0)
	    cout << " > extend scalar: iter: " << iter << endl;
      
	  FOR_IFA {
	    int nof_f = noofa_i[ifa];
	    int nof_l = noofa_i[ifa+1]-1;
	    int ino1 = noofa_v[nof_l];
	    for (int nof = nof_f; nof <= nof_l; ++nof) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      if ((no_count[ino0] <= 0)&&(no_count[ino1] == 1)) {
		(*(data->ptr))[ino0] += (*(data->ptr))[ino1];
		no_count[ino0] -= 1;
		my_done = 0;
	      }
	      else if ((no_count[ino1] <= 0)&&(no_count[ino0] == 1)) {
		(*(data->ptr))[ino1] += (*(data->ptr))[ino0];
		no_count[ino1] -= 1;
		my_done = 0;
	      }
	    }
	  }
	  
	  // switch no_count to positives...
	  
	  FOR_INO {
	    if (no_count[ino] < 0)
	      no_count[ino] = -no_count[ino];
	  }

	  // reduce...
	  updateNoData((*(data->ptr)),ADD_DATA);
	  updateNoData(no_count,ADD_DATA);
	  
	  // now set cv's using the simple average of set face data...
	  FOR_INO {
	    if (no_count[ino] > 1) {
	      (*(data->ptr))[ino] /= (double)no_count[ino];
	      no_count[ino] = 1;
	    }
	  }
	  
	  MPI_Allreduce(&my_done, &done, 1, MPI_INT, MPI_MIN, mpi_comm);

	}
      }
    }
    
    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++) {
      if (data->getDatatype() == NO_DATA) {


	FOR_INO {
	  no_count[ino] = no_flag[ino];
	  // expect zero or 1...
	  assert( (no_count[ino] == 0)||(no_count[ino] == 1) );
	}
	
	int done = 0;
	int iter = 0;
	while (done != 1) {
	  
	  int my_done = 1;
	  
	  iter++;
	  if (mpi_rank == 0)
	    cout << " > extend vector: iter: " << iter << endl;
	  
	  FOR_IFA {
	    int nof_f = noofa_i[ifa];
	    int nof_l = noofa_i[ifa+1]-1;
	    int ino1 = noofa_v[nof_l];
	    for (int nof = nof_f; nof <= nof_l; ++nof) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      if ((no_count[ino0] <= 0)&&(no_count[ino1] == 1)) {
		FOR_I3 (*(data->ptr))[ino0][i] += (*(data->ptr))[ino1][i];
		no_count[ino0] -= 1;
		my_done = 0;
	      }
	      else if ((no_count[ino1] <= 0)&&(no_count[ino0] == 1)) {
		FOR_I3 (*(data->ptr))[ino1][i] += (*(data->ptr))[ino0][i];
		no_count[ino1] -= 1;
		my_done = 0;
	      }
	    }
	  }
	  
	  // switch no_count to positives...
	  
	  FOR_INO {
	    if (no_count[ino] < 0)
	      no_count[ino] = -no_count[ino];
	  }

	  // reduce...
	  updateNoData((*(data->ptr)),ADD_ROTATE_DATA);
	  updateNoData(no_count,ADD_DATA);
	  
	  // now set cv's using the simple average of set face data...
	  FOR_INO {
	    if (no_count[ino] > 1) {
	      FOR_I3 (*(data->ptr))[ino][i] /= (double)no_count[ino];
	      no_count[ino] = 1;
	    }
	  }
	  
	  MPI_Allreduce(&my_done, &done, 1, MPI_INT, MPI_MIN, mpi_comm);
	  
	}
      }
    }
    
    //
    // now check if all data is set. For 2 cases, it is possible that the above routine
    // successfully exits but the cv dat ahas not been set:
    //
    // 1. there was no data set to start with
    // 2. the domain has disconnected regions and one or more of these regions
    //    had no data set.
    //
    // here we simply write a worning, and let the user fix it if they want to.
    //

    int my_unset_count = 0;
    FOR_INO
      if (no_count[ino] <= 0)
	my_unset_count += 1;
    int unset_count;
    MPI_Reduce(&my_unset_count,&unset_count,1,MPI_INT,MPI_SUM,0,mpi_comm);    
    if (mpi_rank == 0) {
      if (unset_count == 0)
	cout << " > all nos successfully set." << endl;
      else
	cout << " > Warning: ******* not all nos were set by extension *******: " << unset_count << endl;
    }

    // cleanup...
    delete[] no_count;
    
  }

  void extendRegisteredCvData() {
    
    double * fa_scalar = NULL;
    double (*fa_vector)[3] = NULL;
    
    int done = 0;
    int iter = 0;
    while (done != 1) {

      iter++;
      if (mpi_rank == 0)
	cout << " > extend iter: " << iter << endl;
      
      // find faces that have one cv set and the other unset. This represents the
      // front of known data...
      int my_done = 1;
      
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++) {
	if (data->getDatatype() == CV_DATA) {

	  // only allocate if we need it...
	  if (fa_scalar == NULL)
	    fa_scalar = new double[nfa];
	  
	  // push set cv data into faces... 
	  FOR_IFA {
	    // zero...
	    fa_scalar[ifa] = 0.0;
	    fa_flag[ifa] = 0;
	    // icv0...
	    const int icv0 = cvofa[ifa][0];
	    if (cv_flag[icv0] > 0) {
	      fa_scalar[ifa] += (*(data->ptr))[icv0];
	      fa_flag[ifa] += 1;
	    }
	    // icv1...
	    const int icv1 = cvofa[ifa][1];
	    if ((icv1 >= 0)&&(icv1 < ncv)&&(cv_flag[icv1] > 0)) {
	      fa_scalar[ifa] += (*(data->ptr))[icv1];
	      fa_flag[ifa] += 1;
	    }
	  }
	  
	  // update across processors...
	  updateFaData(fa_scalar,ADD_DATA);
	  updateFaData(fa_flag,ADD_DATA);
	  
	  // now set cv's using the simple average of set face data...
	  FOR_ICV {
	    // because we are cycling through multiple data here,
	    // this cv_flag could be 0 or -1 here...
	    if (cv_flag[icv] <= 0) {
	      double value = 0.0;
	      int count = 0;
	      // cycle through faces, looking for faces with a value...
	      const int foc_f = faocv_i[icv];
	      const int foc_l = faocv_i[icv+1]-1;
	      for (int foc = foc_f; foc <= foc_l; ++foc) {
		const int ifa = faocv_v[foc];
		if (fa_flag[ifa] > 0) {
		  // must be 1...
		  assert( fa_flag[ifa] == 1 );
		  value += fa_scalar[ifa];
		  count += 1;
		}
	      }
	      // if we got atleast one, then set...
	      if (count > 0) {
		(*(data->ptr))[icv] = value/(double)count;
		cv_flag[icv] = -1;
		// note - put my_done here rather then simply based on whether we have any
		// unset cv's because for certain pathological cases (complex domains and/or 
		// no matched interpolation data) we will not iterate forever...
		my_done = 0;
	      }
	    }
	  }
	  
	}
      }
      
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++) {
	if (data->getDatatype() == CV_DATA) {
	  
	  // only allocate if we need it...
	  if (fa_vector == NULL)
	    fa_vector = new double[nfa][3];
	  
	  // push set cv data into faces... 
	  FOR_IFA {
	    // zero...
	    FOR_I3 fa_vector[ifa][i] = 0.0;
	    fa_flag[ifa] = 0;
	    // icv0...
	    const int icv0 = cvofa[ifa][0];
	    if (cv_flag[icv0] > 0) {
	      FOR_I3 fa_vector[ifa][i] += (*(data->ptr))[icv0][i];
	      fa_flag[ifa] += 1;
	    }
	    // icv1...
	    const int icv1 = cvofa[ifa][1];
	    if ((icv1 >= 0)&&(icv1 < ncv)&&(cv_flag[icv1] > 0)) {
	      FOR_I3 fa_vector[ifa][i] += (*(data->ptr))[icv1][i];
	      fa_flag[ifa] += 1;
	    }
	  }
	  
	  // update across processors...
	  updateFaData(fa_vector,ADD_ROTATE_DATA);
	  updateFaData(fa_flag,ADD_DATA);
	  
	  // now set cv's using the simple average of set face data...
	  FOR_ICV {
	    // because we are cycling through multiple data here,
	    // this cv_flag could be 0 or -1 here...
	    if (cv_flag[icv] <= 0) {
	      double value[3] = { 0.0, 0.0, 0.0 };
	      int count = 0;
	      // cycle through faces, looking for faces with a value...
	      const int foc_f = faocv_i[icv];
	      const int foc_l = faocv_i[icv+1]-1;
	      for (int foc = foc_f; foc <= foc_l; ++foc) {
		const int ifa = faocv_v[foc];
		if (fa_flag[ifa] > 0) {
		  // must be 1...
		  assert( fa_flag[ifa] == 1 );
		  FOR_I3 value[i] += fa_vector[ifa][i];
		  count += 1;
		}
	      }
	      // if we got atleast one, then set...
	      if (count > 0) {
		FOR_I3 (*(data->ptr))[icv][i] = value[i]/(double)count;
		cv_flag[icv] = -1;
		// note - put my_done here rather then simply based on whether we have any
		// unset cv's because for certain pathological cases (complex domains and/or 
		// no matched interpolation data) we will not iterate forever...
		my_done = 0;
	      }
	    }
	  }

	}
      }

      // if we reset my_done, then there were some cv's set...
      if (my_done == 0) 
	FOR_ICV
	  if (cv_flag[icv] == -1)
	    cv_flag[icv] = 1;
   
      MPI_Allreduce(&my_done, &done, 1, MPI_INT, MPI_MIN, mpi_comm);
      
    }
    
    //
    // now check if all data is set. For 2 cases, it is possible that the above routine
    // successfully exits but the cv dat ahas not been set:
    //
    // 1. there was no data set to start with
    // 2. the domain has disconnected regions and one or more of these regions
    //    had no data set.
    //
    // here we simply write a worning, and let the user fix it if they want to.
    //

    int my_unset_count = 0;
    FOR_ICV
      if (cv_flag[icv] <= 0)
	my_unset_count += 1;
    int unset_count;
    MPI_Reduce(&my_unset_count,&unset_count,1,MPI_INT,MPI_SUM,0,mpi_comm);    
    if (mpi_rank == 0) {
      if (unset_count == 0)
	cout << " > all cvs successfully set." << endl;
      else
	cout << " > Warning: ******* not all cvs were set by extension *******: " << unset_count << endl;
    }

    // cleanup...

    if (fa_scalar != NULL) delete[] fa_scalar;
    if (fa_vector != NULL) delete[] fa_vector;
    
  }
  
  void dumpBoundaryZoneSummary() {

    if (mpi_rank == 0)
      cout << "dumpBoundaryZoneSummary()" << endl;

    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
   
      if (zone->getKind() == FA_ZONE_BOUNDARY) {

	int my_nfa_bc = zone->ifa_l - zone->ifa_f + 1;
	int nfa_bc;
	MPI_Reduce(&my_nfa_bc,&nfa_bc,1,MPI_INT,MPI_SUM,0,mpi_comm);    

	if (mpi_rank == 0)
	  cout << " > " << zone->getName() << ", face count: " << nfa_bc << ", area: ";
	
	double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
	for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	  double x_fa[3],n_fa[3];
	  calcFaceCenterAndNormal(x_fa,n_fa,ifa);
	  FOR_I3 my_buf[i] += n_fa[i];
	  my_buf[3] += sqrt( n_fa[0]*n_fa[0] + n_fa[1]*n_fa[1] + n_fa[2]*n_fa[2] );
	}

	double buf[4];
	MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);    
	
	if (mpi_rank == 0)
	  cout << buf[3] << ", outward normal: " << buf[0] << " " << buf[1] << " " << buf[2] << endl;
	
      }

    }

  }

  void preproHook() { }

  /*
    void preproHook() {

    if (mpi_rank == 0)
    cout << "preproHook()" << endl;
    
    double (*x_cv)[3] = new double[ncv][3];
    calcCvCenterInit(x_cv);
    
    double (*rhou)[3] = getRegisteredDoubleVector("RHOU");
    FOR_ICV {
    // remove the tangential component. A unit vector in the tangential direction is...
    double v[3];
    v[0] = 0.0;
    v[1] = -x_cv[icv][2];
    v[2] = x_cv[icv][1];
    double inv_mag = 1.0/sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    FOR_I3 v[i] *= inv_mag;
    double rhouv = rhou[icv][0]*v[0] + rhou[icv][1]*v[1] + rhou[icv][2]*v[2];
    FOR_I3 rhou[icv][i] -= rhouv*v[i];
    }

    delete[] x_cv;
    
    }
  */
  
};

void dumpPreproHelp() {
  
  cout << "\nprepro: Parallel mesh and solution preprocessing\n================================================\n\nUsage:\n\n  mpirun -np 4 ./prepro" << endl;
  cout << "\nor run in single-processor mode: \n\n  ./prepro" << endl;
  cout << "\nPrepro's behavior is controlled from parameters parsed from the default input file \"prepro.in\"." << endl;
  cout << "Input parameters can also be added at the commandline using double hyphen notation as follows:" << endl;
  cout << "\n  ./prepro --RESTART filename.cas --TECPLOT" << endl;
  cout << "\nNote that commandline parameters take precidence over input file parameters." << endl; 
  cout << "\nParameters\n------------" << endl;
  cout << "\nRESTART <filename>" << endl;
  cout << "\nRESTART defines the source file for mesh (and potentially solution data). Supported filters include" << endl;
  cout << "Fluent *.cas or Gambit *.msh files (ASCII only), and cdp2.4 *.cdp files. Other extensions will be assumed" << endl;
  cout << "Mum restart/result files." << endl;
  cout << "\nMultiple RESTART files can be read in and combined by simply specifying multiple restart parameters." << endl;
  cout << "\nRENAME <old_zone_name> <new_zone_name> [<periodic_kind> <periodic_data...>]" << endl;
  cout << "\nexample:\n" << endl;
  cout << "\nRENAME can be used to change face zone names, introduce or remove mesh periodicity," << endl;
  cout << "and reconnect zones from single or multiple input meshes. example:" << endl;
  cout << "RENAME y0-inside periodic_y0 CYL_X 22.5" << endl;
  cout << "\nSupported periodicities include:" << endl;
  cout << "\n  CART <dx> <dy> <dz>" << endl; 
  cout << "\n  CYL_X <degrees>" << endl; 
  cout << "\n  CYL_Y <degrees>" << endl; 
  cout << "\n  CYL_Z <degrees>" << endl; 
  cout << "\nUsing the RENAME param:" << endl;
  cout << endl;
  cout << "When renaming a zone to a periodic zone, you must" << endl;
  cout << "provide the type of periodicity and its associated transform on" << endl;
  cout << "the same line as the RENAME parameter." << endl;
  cout << "For example, to RENAME a zone named \"x0\" to \"periodic_x0\" and" << endl;
  cout << "introduce Cartesian periodicity with another zone (like periodic_x1)" << endl;
  cout << "located +1 units in x, use the folowing param:" << endl; 
  cout << endl;
  cout << "RENAME x0 periodic_x0 CART 1.0 0.0 0.0" << endl;
  cout << "RENAME x1 periodic_x1 CART -1.0 0.0 0.0" << endl;
  cout << endl;
  cout << "for cylindrical periodicity, use the keyword CYL_Z (or CYL_X or CYL_Y)" << endl;
  cout << "followed by an angle in degrees. e.g.:" << endl;
  cout << endl;
  cout << "RENAME y0 periodic_theta0 CYL_Z 22.5" << endl;
  cout << "RENAME y1 periodic_theta1 CYL_Z -22.5" << endl;
  cout << "\nThe transformation you specify should be that required to get from the" <<endl;
  cout << "zone you are renaming to the periodic pair." << endl;
  cout << "\nUsing the TRANSFORM param:" << endl;
  cout << endl;
  cout << "TRANSFORM TRANSLATE <dx> <dy> <dz>" << endl;
  cout << "TRANSFORM ROTATE_X <degrees>" << endl;
  cout << "TRANSFORM ROTATE_Y <degrees>" << endl;
  cout << "TRANSFORM ROTATE_Z <degrees>" << endl;
  cout << "TRANSFORM SCALE <factor>" << endl;
  cout << "TRANSFORM SCALE <factor_x> <factor_y> <factor_z>" << endl;
  cout << "\nINTERP_FROM_RESTART <old_restart_file>" << endl;
  cout << "\nREGISTER_INT_VALUE STEP REGISTER_DOUBLE_VALUE TIME" << endl;
  
}

int main(int argc,char * argv[]) {

  // intialize the MPI environment - on a single processor,
  // this does nothing.

  MPI_Init(&argc,&argv);

  // try/catch is for catching errors that the code "throws"
  
  try {
    
    // this also does nothing in a serial case, in parallel 
    // is gets some MPI stuff ready for us to use...
    
    initMpiStuff();

    //
    // a ParamMap is a way to manage parameters that the user
    // specifies. For example, RESTART = filename, these can 
    // come from a file (typically call "prepro.in") or can be 
    // parsed form the command line arguments using -- notation, e.g.
    // 
    // ./prepro --RESTART filename
    // 
    
    ParamMap params;
    params.addParamsFromArgs(argc,argv);
    params.addParamsFromFile("prepro.in");

    if (params.checkParam("HELP")) {
      if (mpi_rank == 0)
	dumpPreproHelp();
    }
    else {

      if (mpi_rank == 0)
	params.dump();

      // PC HACK...
      //cout << "here is a as a string: " << params.getStringParam("a") << endl;
      //cout << "here is a as a double: " << params.getDoubleParam("a") << endl;
    
      PreUgp ugp;
      ugp.registerData(params);
      bool first = true;

      // =====================================================
      // step 1 - minimum initialization...
      // =====================================================
    
      vector<Param> *paramVec;
      assert(params.getParam(paramVec,"RESTART"));
    
      // allow for multiple RESTART's
    
      for (vector<Param>::iterator p = paramVec->begin(); p != paramVec->end(); p++) {
	
	if (first) {
	  
	  first = false;
	  
	  string restart = p->getString();
	  string ext = getFileNameExtension(restart,".");
	  if (ext == "cdp") {
	    // old CDP format...
	    CdpFilter cf(restart);
	    cf.minInitUgp(&ugp);
	  } 
	  else if ((ext == "msh")||(ext == "cas")) {
	    // Fluent format...
	    MshFilter mf(restart);
	    mf.initUgpMin(&ugp);
	  } 
	  else {
	    // assume this is a standard mum file...
	    ugp.readRestart(restart);
	  }

	}
	else {

	  Ugp this_ugp;
	  this_ugp.registerData(params);

	  // ===========================================================
	  // read the restart file...
	  // ===========================================================
	  string restart = p->getString();
	  string ext = getFileNameExtension(restart,".");
	  if (ext == "cdp") {
	    // old CDP format...
	    CdpFilter cf(restart);
	    cf.minInitUgp(&this_ugp);
	  } 
	  else if ((ext == "msh")||(ext == "cas")) {
	    // Fluent format...
	    MshFilter mf(restart);
	    mf.initUgpMin(&this_ugp);
	  } 
	  else {
	    // assume this is a standard mum file...
	    this_ugp.readRestart(restart);
	  }
	  
	  // ===========================================================
	  // look for additional arguments associated with this RESTART...
	  // ===========================================================
	  int i = 2;
	  while (i < p->getSize()) {
	    string token = p->getString(i++);
	    if (token == "TRANSFORM") {
	      token = p->getString(i++);
	      if (token == "SCALE") {
		cerr << "Error: TRANSFORM SCALE not implemented for individual RESTARTs" << endl;
		throw(-1);
	      }
	      else {
		if (mpi_rank == 0)
		  cerr << "Error: unsupported TRANSFORM token: " << token << endl;
		throw(-1);
	      }
	    } 
	    else {
	      if (mpi_rank == 0)
		cerr << "Error: unsupported RESTART token: " << token << endl;
	      throw(-1);
	    }
	  }
	  
	  // add to the PreUgp...
	  ugp += this_ugp;
	  
	}

      }
	
      // ====================================
      // we have applied all RESTART args...
      // now look for global transforms...
      // ====================================
    
      if (params.getParam(paramVec,"TRANSFORM")) {
      
	for (vector<Param>::iterator p = paramVec->begin(); p != paramVec->end(); p++) {
	  
	  string transformName = p->getString();
	    
	  if (transformName.compare("TRANSLATE") == 0) {

	    double dx[3];
	    dx[0] = p->getDouble(2);
	    dx[1] = p->getDouble(3);
	    dx[2] = p->getDouble(4);
	      
	    if (mpi_rank == 0)
	      cout << " > applying TRANSFORM TRANSLATE " << dx[0] << " " << dx[1] << " " << dx[2] << endl;
	    ugp.translate(dx);
	      
	  }
	  else if (transformName.compare("SHEAR_XY") == 0) {
	      
	    double alpha = p->getDouble(2);
	      
	    if (mpi_rank == 0)
	      cout << " > applying TRANSFORM SHEAR_XY " << alpha << endl;
	    ugp.shearXY(alpha);
	      
	  }
	  else if (transformName.compare("ROTATE_X") == 0) {

	    double alpha = p->getDouble(2);

	    if (mpi_rank == 0)
	      cout << " > applying TRANSFORM ROTATE_X " << alpha << " degrees" << endl;
	    ugp.rotateX(alpha);
	      
	  }
	  else if (transformName.compare("ROTATE_Z") == 0) {

	    double alpha = p->getDouble(2);

	    if (mpi_rank == 0)
	      cout << " > applying TRANSFORM ROTATE_Z " << alpha << " degrees" << endl;
	    ugp.rotateZ(alpha);
	      
	  }
	  else if (transformName.compare("SCALE") == 0) {
	      
	    // the scale comand can be followed by either 1 or 3 doubles...
	    switch (p->getSize()) {
	    case 3:
	      {
		double factor = p->getDouble(2);
		if (mpi_rank == 0)
		  cout << " > applying TRANSFORM SCALE " << factor << endl;
		ugp.scale(factor);
	      }
	      break;
	    case 5:
	      {
		double factor[3];
		factor[0] = p->getDouble(2);
		factor[1] = p->getDouble(3);
		factor[2] = p->getDouble(4);
		if (mpi_rank == 0)
		  cout << " > applying TRANSFORM SCALE " << factor[0] << " " << factor[1] << " " << factor[2] << endl;
		ugp.scale(factor);
	      }
	      break;
	    default:
	      if (mpi_rank == 0)
		cerr << "Error: TRANSFORM SCALE expects either 1 (isotropic) or 3 (anisotropic) doubles" << endl;
	      throw(-1);
	    }
	      

	  }
	  else {

	    cerr << "Error: unrecognized TRANSFORM: " << transformName << endl;
	    throw(-1);

	  }

	}

      }

      // ===================================
    
      // zone renaming, periodicity, etc...
    
      if (params.getParam(paramVec,"RENAME")) {
      
	for (vector<Param>::iterator p = paramVec->begin(); p != paramVec->end(); p++) {
	  
	  string oldName = p->getString();
	  
	  FaZone * zone = ugp.getFaZone(oldName);
	  if (zone == NULL) {
	    if (mpi_rank == 0)
	      cerr << "Error: RENAME: could not find zone named: " << oldName << endl;
	    throw(-1);
	  }

	  // the next string should be the new name...
	    
	  string newName = p->getString(2);
	    
	  if (mpi_rank == 0) 
	    cout << " > RENAMING zone " << oldName << " to " << newName << endl;
	    
	  // set the new name...

	  zone->setName(newName);

	  // and pull the zone index. Recall that the min initialization has the face zone
	  // index in fa_flag...

	  int index = zone->getIndex();
	    
	  // the user can add periodicity to a grid by using the RENAME command
	  // and renaming the zone so that is starts with "periodic" or "PERIODIC"...
	    
	  if ((newName.compare(0,8,"periodic") == 0)||(newName.compare(0,8,"PERIODIC") == 0)) {
	    
	    // look for the periodic transform - it will tell us what kind of periodicity
	    // is being introduced...
	    
	    string periodicType = p->getString(3);
	    if (periodicType == "CART") {
	      if (mpi_rank == 0) 
		cout << " > applying CART periodicity to zone: " << newName << endl;
	      zone->setKind(FA_ZONE_PERIODIC_CART);
	      // should be 3 doubles that describe the transformation...
	      double periodic_data[3];
	      periodic_data[0] = p->getDouble(4);
	      periodic_data[1] = p->getDouble(5);
	      periodic_data[2] = p->getDouble(6);
	      zone->setPeriodicData(periodic_data);
	    }

	    else if (periodicType == "CYL_X") {
	      if (mpi_rank == 0) 
		cout << " > applying CYL_X periodicity to zone: " << newName << endl;
	      zone->setKind(FA_ZONE_PERIODIC_CYL_X);
	      // should be one double describing the angle in degrees...
	      double alpha = p->getDouble(4);
	      double periodic_data[3];
	      periodic_data[0] = cos(M_PI*alpha/180.0);
	      periodic_data[1] = sin(M_PI*alpha/180.0);
	      periodic_data[2] = 0.0;
	      zone->setPeriodicData(periodic_data);
	    }


	    else if (periodicType == "CYL_Y") {
	      if (mpi_rank == 0) 
		cout << " > applying CYL_Y periodicity to zone: " << newName << endl;
	      zone->setKind(FA_ZONE_PERIODIC_CYL_Y);
	      // should be one double describing the angle in degrees...
	      double alpha = p->getDouble(4);
	      double periodic_data[3];
	      periodic_data[0] = cos(M_PI*alpha/180.0);
	      periodic_data[1] = sin(M_PI*alpha/180.0);
	      periodic_data[2] = 0.0;
	      zone->setPeriodicData(periodic_data);
	    }


	    else if (periodicType == "CYL_Z") {
	      if (mpi_rank == 0) 
		cout << " > applying CYL_Z periodicity to zone: " << newName << endl;
	      zone->setKind(FA_ZONE_PERIODIC_CYL_Z);
	      // should be one double describing the angle in degrees...
	      double alpha = p->getDouble(4);
	      double periodic_data[3];
	      periodic_data[0] = cos(M_PI*alpha/180.0);
	      periodic_data[1] = sin(M_PI*alpha/180.0);
	      periodic_data[2] = 0.0;
	      zone->setPeriodicData(periodic_data);
	    }
	    else {
	      if (mpi_rank == 0)
		cerr << "Error: RENAME: unrecognized type of periodicity: " << periodicType << endl;
	      throw(-1);
	    }

	    // check any face index in cvofa[ifa][1]...
	    for (int ifa = 0; ifa < ugp.getNfa(); ++ifa)
	      if (ugp.getFaFlag(ifa) == index) 
		assert( ugp.getCvOfFa(ifa,1) == -1 );
	      
	  }
	  else if ( zone->isPeriodic() ) {
	    
	    // if the zone was periodic, we can delete periodicity as follows...
	    zone->setKind(FA_ZONE_BOUNDARY);
	    
	    // check any face index in cvofa[ifa][1]...
	    for (int ifa = 0; ifa < ugp.getNfa(); ++ifa)
	      if (ugp.getFaFlag(ifa) == index) 
		ugp.cvofa[ifa][1] = -1;
	      
	  }
	  else {
	    
	    // ensure the zone is a boundary zone...
	    assert( zone->getKind() == FA_ZONE_BOUNDARY );
	    
	  }
	  
	}
	
      }
      
      // ===================================
    
      // preinit ensures all faces are set as expected - e.g. if any
      // periodicity has been added, or if zones need to be reconnected,
      // or if formerly periodic zones have been turned into boundaries...
      ugp.preinit();
    
      // normal init...
      ugp.init();
    
      // postinit synchronizes the nodes...
      ugp.postinit();
    
      // =====================================================
      // now look for data to interpolate...
      // =====================================================
    
      Param * p;
      if (params.getParam(p,"INTERP_FROM_RESTART")) {
      
	// interpolate any registered data from another restart file...
	if ( mpi_rank == 0)
	  cout << "INTERP_FROM_RESTART: " << p->getString() << endl;
      
	PreUgp ugp2;
	ugp2.registerData(params);
	ugp2.readRestart(p->getString());
	// NO init here for ugp2!
	
	// to build the cv data map, we need x_cv for both ugp and ugp2... 
	if (ugp.x_cv == NULL) {
	  ugp.x_cv = new double[ugp.ncv][3];
	  ugp.calcCvCenterInit(ugp.x_cv);
	}
	
	assert(ugp2.x_cv == NULL);
	ugp2.x_cv = new double[ugp2.ncv][3];
	ugp2.calcCvCenterNoInit(ugp2.x_cv);

	/*
	// HACK - transform the ugp2 cv's (and nodes?) to interpolate data into
	// another part of a full geometry and be sure to comment out the 
	// extend data routoines below...
	cout << "WARNING: transforming cv centroids" << endl;
	for (int icv = 0; icv < ugp2.ncv; ++icv) {
	  // flip z...
	  ugp2.x_cv[icv][2] *= -1.0;
	}
	*/

	// and ugp's adt...
	ugp.useCvAdt();
	
	// build the map: its structures lives in ugp...
	// Note that this routine returns cv_flag with 0 in cv's 
	// that will not recieve data, and >= 1 in cells getting data...
	ugp.buildNoInitCvDataMap2(ugp2.x_cv,ugp2.cvora);
	
	// build a nodal data map too...
	ugp.buildNoInitNoDataMap(ugp2.x_no,ugp2.noora);	  
	
	// ============================================
	// now interpolate
	// ============================================
	if (mpi_rank == 0)
	  cout << "First-order interpolation of registered data..." << endl;
	
	// now loop thru registered data and update in Ugp...
	for (list<IntValue>::iterator data2 = ugp2.intValueList.begin(); data2 != ugp2.intValueList.end(); data2++) {
	  if (mpi_rank == 0)
	    cout << " > setting int value: " << data2->getName() << endl;
	  IntValue * data = ugp.getIntValueData(data2->getName());
	  assert( data != NULL );
	  *(data->ptr) = *(data2->ptr);
	}

	for (list<DoubleValue>::iterator data2 = ugp2.doubleValueList.begin(); data2 != ugp2.doubleValueList.end(); data2++) {
	  if (mpi_rank == 0)
	    cout << " > setting double value: " << data2->getName() << endl;
	  DoubleValue * data = ugp.getDoubleValueData(data2->getName());
	  assert( data != NULL );
	  *(data->ptr) = *(data2->ptr);
	}
	 
	for (list<DoubleScalar>::iterator data2 = ugp2.doubleScalarList.begin(); data2 != ugp2.doubleScalarList.end(); data2++) {
	  if (data2->getDatatype() == CV_DATA) {
	    if (mpi_rank == 0)
	      cout << " > interpolating cv double scalar: " << data2->getName() << endl;
	    DoubleScalar * data = ugp.getDoubleScalarData(data2->getName());
	    assert( data != NULL );
	    ugp.updateNoInitCvData(*(data->ptr),*(data2->ptr)); // note order of args!
	  }
	  else if (data2->getDatatype() == NO_DATA) {
	    if (mpi_rank == 0)
	      cout << " > interpolating no double scalar: " << data2->getName() << endl;
	    DoubleScalar * data = ugp.getDoubleScalarData(data2->getName());
	    assert( data != NULL );
	    ugp.updateNoInitNoData(*(data->ptr),*(data2->ptr)); // note order of args!
	  }
	  else {
	    if (mpi_rank == 0)
	      cerr << "Error: unsupported data type: " << data2->getDatatype() << endl;
	    throw(-1);
	  }
	}
	
	for (list<DoubleVector>::iterator data2 = ugp2.doubleVectorList.begin(); data2 != ugp2.doubleVectorList.end(); data2++) {
	  if (data2->getDatatype() == CV_DATA) {
	    if (mpi_rank == 0)
	      cout << " > interpolating cv double vector: " << data2->getName() << endl;
	    DoubleVector * data = ugp.getDoubleVectorData(data2->getName());
	    assert( data != NULL );
	    ugp.updateNoInitCvData(*(data->ptr),*(data2->ptr)); // note order of args!
	  }
	  else if (data2->getDatatype() == NO_DATA) {
	    if (mpi_rank == 0)
	      cout << " > interpolating no double vector: " << data2->getName() << endl;
	    DoubleVector * data = ugp.getDoubleVectorData(data2->getName());
	    assert( data != NULL );
	    ugp.updateNoInitNoData(*(data->ptr),*(data2->ptr)); // note order of args!
	  }
	  else {
	    if (mpi_rank == 0)
	      cerr << "Error: unsupported data type: " << data2->getDatatype() << endl;
	    throw(-1);
	  }
	}
	
	// now do data extension into neighboring cv's...
	
	if (mpi_rank == 0)
	  cout << "Extending registered data into unset neighbors..." << endl;
	
	// HACK - skip the extend if you have multiple data
	//cout << "WARNING: skipping extend" << endl;
	ugp.extendRegisteredCvData();
	ugp.extendRegisteredNoData();

      }

      // =====================================================
      // step 3 - dump result
      // =====================================================
      
      //if (params.checkParam("ZONE_SUMMARY"))
      ugp.dumpBoundaryZoneSummary();
      
      if (params.checkParam("TECPLOT"))
	ugp.writeBoundaryFacesTecplot("faces");
      
      if (params.checkParam("MSH"))
	ugp.writeFluentMsh("fluent.msh");
      
      if (params.checkParam("STL"))
	ugp.writeBoundaryFacesStl();
      
      ugp.preproHook();
      
      // setCvGroup calls the virtual function setFaFlagForCvGroupHook() - see above...
      /*
	ugp.fa_delta_max = params.getDoubleParam("FA_DELTA_MAX");
	ugp.setCvGroup();
      */

      if (!params.checkParam("NO_RESTART"))
	ugp.writeRestart();
      
      // let the user know we are done...
      if (mpi_rank == 0) 
	cout << "prepro finished successfully" << endl;


    }

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
  
  // shut down MPI stuff...
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return(0);
  
}

