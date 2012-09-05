#ifndef AVERAGE_H
#define AVERAGE_H

#include <iostream>
#include <vector>

#ifdef NO_ASSERT
#include <assert.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::vector;

#include <algorithm> // required for std::sort on uP

#include "MpiStuff.h"
using namespace MpiStuff;

// kinds of averaging
enum AverageKinds {
  AVERAGE_NONE,
  AVERAGE_Z,
  AVERAGE_X_Z,
  AVERAGE_THETAX
};

class Average1dCompare : public std::binary_function<int&, int&, bool> {

private:
  
  const double * v1;
  double tol;
  
public: 
  
  Average1dCompare(const double * v1,const double tol) {
    this->v1 = v1;
    this->tol = tol;
  }
  
  bool operator()(const int& a,const int& b) {
    return( v1[a] < v1[b]-tol );
  }
  
};

class Average2dCompare : public std::binary_function<int&, int&, bool> {
  
private:
  
  const double * v1;
  const double * v2;
  double tol;
  
public: 
  
  Average2dCompare(const double * v1,const double * v2,const double tol) {
    this->v1 = v1;
    this->v2 = v2;
    this->tol = tol;
  }

  bool operator()(const int& a,const int& b) {
    return( ( v1[a] < v1[b]-tol )||( ( v1[a] < v1[b]+tol )&&( v2[a] < v2[b]-tol ) ) );
  }

};

class Average {
  
private:
  
  int n,nbins;
  int * bin_of_index; // [n] dimension of variables, typically nno, or ncv
  double * weight_of_index; // [n]
  double * bin_weight; // [nbins] number of bins
  double * bin_value;
  
  int * pack_count;
  int * pack_disp;
  int npack;
  int * pack_index;
  double * pack_buffer;
  
  int * unpack_count;
  int * unpack_disp;
  int nunpack;
  int * unpack_index;
  double * unpack_buffer;

  int kind;

public:
  
  Average(const double * v1,const double * weight,const int n,const double tol) {
    
    init0();
    init(v1,weight,n,tol);
    
  }
  
  Average(const double * v1,const double * v2,const double * weight,
	  const int n,const double tol) {

    init0();
    init(v1,v2,weight,n,tol);
    
  }

 public:

  void setKind(const int newkind) {
    kind = newkind;
  }

  int getKind() const {
    return(kind);
  }
  
  int size() const {
    return(n);
  }

 private:

  void init0() {

    n = 0;
    bin_of_index = NULL;
    weight_of_index = NULL;
    
    nbins = 0;
    bin_weight = NULL;
    bin_value = NULL;

    pack_count = NULL;
    pack_disp = NULL;
    npack = 0;
    pack_index = NULL;
    pack_buffer = NULL;
    
    unpack_count = NULL;
    unpack_disp = NULL;
    nunpack = 0;
    unpack_index = NULL;
    unpack_buffer = NULL;

    kind = AVERAGE_NONE;

  }
  
  void init(const double * v1,const double * weight,const int n,const double tol) {
	   
    // initialization of averaging class for 1d average - basically bins common points
    // that have the same 1d location within "tol"...
    
    if (mpi_rank == 0)
      cout << "Average::init() for 1D bins, tol: " << tol << endl;
    
    this->n = n;
    double * v1_bin = NULL;
    
    {
    
      vector<int> indexVec(n);
      for (int i = 0; i < n; i++)
        indexVec[i] = i;
      
      std::sort(indexVec.begin(),indexVec.end(),Average1dCompare(v1,tol));
      
      assert( bin_of_index == NULL );
      bin_of_index = new int[n];

      for (int iter = 0; iter < 2; iter++) {
        nbins = 0;
        double this_v1;
        for (int i = 0; i < n; i++) {
          int index = indexVec[i];
          if ( (nbins == 0)||( fabs(this_v1-v1[index]) > tol ) ) {
            this_v1 = v1[index];
            if (iter == 1) {
              v1_bin[nbins] = this_v1;
            }
            nbins++;
          }
          if (iter == 1)
            bin_of_index[index] = nbins-1;
        }
        if (iter == 0) {
          v1_bin = new double[nbins];
        }
      }
      
    }
    
    // assign bin_weight and weight_of_index, do some checking...
    assert( bin_weight == NULL );
    bin_weight = new double[nbins];
    for (int ibin = 0; ibin < nbins; ibin++)
      bin_weight[ibin] = 0.0;
    assert( weight_of_index == NULL );
    weight_of_index = new double[n];
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      assert( fabs(v1[i] - v1_bin[ibin]) < tol );
      weight_of_index[i] = weight[i];
      bin_weight[ibin] += weight[i];
    }
    
    //cout << "n, nbins: " << n << " " << nbins << endl;
    
    // we also need to exchange with processor neighbors. Everyone compute and 
    // exchange their v1 bounding boxes first...
    
    double local_bbox[2];
    local_bbox[0] = 1.0E+20; // v1 min
    local_bbox[1] = -1.0E+20; // v1 max
    for (int i = 0; i < n; i++) {
      local_bbox[0] = min(local_bbox[0],v1[i]);
      local_bbox[1] = max(local_bbox[1],v1[i]);
    }
    double *global_bbox = new double[2*mpi_size];
    MPI_Allgather(local_bbox,2,MPI_DOUBLE,global_bbox,2,MPI_DOUBLE,mpi_comm);

    // now decide how many v1's we are going to send to each rank...
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * send_buffer_bin = NULL;
    double * send_buffer_double = NULL;
    int send_count_sum;
    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) if (rank != mpi_rank) {
        for (int i = 0; i < nbins; i++) {
          if ( (v1_bin[i] > global_bbox[2*rank+0]-tol)&&(v1_bin[i] < global_bbox[2*rank+1]+tol) ) {
            if (iter == 0) {
              send_count[rank] += 1;
            }
            else {
              send_buffer_bin[send_disp[rank]] = i; // store the sent bin
              send_buffer_double[send_disp[rank]  ] = v1_bin[i];
              send_disp[rank] += 1;
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
        send_buffer_bin = new int[send_count_sum];
        send_buffer_double = new double[send_count_sum];
      }
    }

    // we stored the bin in send_buffer_bin, so no need to keep global bbox's...
    delete[] global_bbox;
    
    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    double * recv_buffer_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buffer_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buffer_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    
    delete[] send_buffer_double;
    
    // will hold a number >= 0 if we got a match in one of our bins, -1 otherwise...
    int * recv_buffer_int = new int[recv_count_sum];
    for (int i = 0; i < recv_count_sum; i++) {
      // following is n^2 in nbins -- sorry. Could do a 
      // second sort and save some time here, but it
      // is just pre-processing, so n^2 with small n is OK?!...
      recv_buffer_int[i] = -1; 
      double this_v1 = recv_buffer_double[i];
      for (int ibin = 0; ibin < nbins; ibin++) {
        if ( fabs(this_v1-v1_bin[ibin]) < tol ) {
          //if (mpi_rank == 0) 
          //cout << "got match" << endl;
          recv_buffer_int[i] = ibin;
          break;
        }
      }
    }
    
    delete[] recv_buffer_double;
    
    // send this matching info back to the sending processes...
    
    int * send_buffer_int = new int[send_count_sum];
    MPI_Alltoallv(recv_buffer_int,recv_count,recv_disp,MPI_INT,
                  send_buffer_int,send_count,send_disp,MPI_INT,mpi_comm);
    
    // now recv and send exchange roles. During averaging, recv will actually 
    // collect its various bins and send them to the send processes for addition
    // to their bin sums, then normalization to produce averaging...
    
    // ---------------------------------------------
    // the recv side becomes the "pack" side...
    // ---------------------------------------------
    pack_count = new int[mpi_size];
    pack_disp = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++)
      pack_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) {
        for (int i = recv_disp[rank]; i < recv_disp[rank]+recv_count[rank]; i++) {
          if (recv_buffer_int[i] >= 0) {
            if (iter == 0) {              
              pack_count[rank] += 1;
            }
            else {
              pack_index[pack_disp[rank]] = recv_buffer_int[i];
              pack_disp[rank] += 1;
            }
          }
        }
      }
      // calculate pack_disp on both iterations...
      pack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
        pack_disp[i] = pack_count[i-1] + pack_disp[i-1];
      // on the first time through, allocate the pack stuff...
      if (iter == 0) {
        npack = pack_disp[mpi_size-1] + pack_count[mpi_size-1];
        pack_index = new int[npack];
        pack_buffer = new double[npack];
      }
    }

    delete[] recv_count;
    delete[] recv_disp;
    delete[] recv_buffer_int;
    
    // ---------------------------------------------
    // and the send side becomes the unpack side...
    // ---------------------------------------------
    unpack_count = new int[mpi_size];
    unpack_disp = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++)
      unpack_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) {
        for (int i = send_disp[rank]; i < send_disp[rank]+send_count[rank]; i++) {
          if (send_buffer_int[i] >= 0) {
            if (iter == 0) {              
              unpack_count[rank] += 1;
            }
            else {
              unpack_index[unpack_disp[rank]] = send_buffer_bin[i]; // recall we stored above
              unpack_disp[rank] += 1;
            }
          }
        }
      }
      // calculate pack_disp on both iterations...
      unpack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
        unpack_disp[i] = unpack_count[i-1] + unpack_disp[i-1];
      // on the first time through, allocate the unpack stuff...
      if (iter == 0) {
        nunpack = unpack_disp[mpi_size-1] + unpack_count[mpi_size-1];
        unpack_index = new int[nunpack];
        unpack_buffer = new double[nunpack];
      }
    }
    
    delete[] send_count;
    delete[] send_disp;
    delete[] send_buffer_bin;
    delete[] send_buffer_int;
    
    // ---------------------------------------------------------
    // now check bin correspondence by exchanging v1...
    // ---------------------------------------------------------
    
    // v1...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = v1_bin[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      assert( fabs(unpack_buffer[i] - v1_bin[unpack_index[i]]) < tol );
    
    // finally sum the bin_weight's...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = bin_weight[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      bin_weight[unpack_index[i]] += unpack_buffer[i];
    
    // final stuff...
    assert( bin_value == NULL );
    bin_value = new double[nbins];
    
  }
  
  void init(const double * v1,const double * v2,const double * weight,const int n,const double tol) {

    // initialization of averaging class for 2d average - basically bins common points
    // that have the same 2d location within "tol"...

    if (mpi_rank == 0)
      cout << "Average::init() for 2D bins, tol: " << tol << endl;
    
    this->n = n;
    double * v1_bin = NULL;
    double * v2_bin = NULL;
    
    {
      
      vector<int> indexVec(n);
      for (int i = 0; i < n; i++)
        indexVec[i] = i;
      
      std::sort(indexVec.begin(),indexVec.end(),Average2dCompare(v1,v2,tol));
      
      assert( bin_of_index == NULL );
      bin_of_index = new int[n];

      for (int iter = 0; iter < 2; iter++) {
        nbins = 0;
        double this_v1,this_v2;
        for (int i = 0; i < n; i++) {
          int index = indexVec[i];
          if ( (nbins == 0)||( fabs(this_v1-v1[index]) > tol )||( fabs(this_v2-v2[index]) > tol ) ) {
            this_v1 = v1[index];
            this_v2 = v2[index];
            if (iter == 1) {
              v1_bin[nbins] = this_v1;
              v2_bin[nbins] = this_v2;
            }
            nbins++;
          }
          if (iter == 1)
            bin_of_index[index] = nbins-1;
        }
        if (iter == 0) {
          v1_bin = new double[nbins];
          v2_bin = new double[nbins];
        }
      }
      
    }
    
    // assign bin_weight and weight_of_index, do some checking...
    assert( bin_weight == NULL );
    bin_weight = new double[nbins];
    for (int ibin = 0; ibin < nbins; ibin++)
      bin_weight[ibin] = 0.0;
    assert( weight_of_index == NULL );
    weight_of_index = new double[n];
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      assert( fabs(v1[i] - v1_bin[ibin]) < tol );
      assert( fabs(v2[i] - v2_bin[ibin]) < tol );
      weight_of_index[i] = weight[i];
      bin_weight[ibin] += weight[i];
    }
    
    //cout << "n, nbins: " << n << " " << nbins << endl;
    
    // we also need to exchange with processor neighbors. Everyone compute and 
    // exchange their v1,v2 bounding boxes first...
    
    double local_bbox[4];
    local_bbox[0] = 1.0E+20; // v1 min
    local_bbox[1] = -1.0E+20; // v1 max
    local_bbox[2] = 1.0E+20; // v2 min
    local_bbox[3] = -1.0E+20; // v2 max
    for (int i = 0; i < n; i++) {
      local_bbox[0] = min(local_bbox[0],v1[i]);
      local_bbox[1] = max(local_bbox[1],v1[i]);
      local_bbox[2] = min(local_bbox[2],v2[i]);
      local_bbox[3] = max(local_bbox[3],v2[i]);
    }
    double *global_bbox = new double[4*mpi_size];
    MPI_Allgather(local_bbox,4,MPI_DOUBLE,global_bbox,4,MPI_DOUBLE,mpi_comm);

    // now decide how many v1,v2 pairs we are going to send to each rank...
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * send_buffer_bin = NULL;
    double * send_buffer_double = NULL;
    int send_count_sum;
    for (int i = 0; i < mpi_size; i++)
      send_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) if (rank != mpi_rank) {
        for (int i = 0; i < nbins; i++) {
          if ( (v1_bin[i] > global_bbox[4*rank+0]-tol)&&(v1_bin[i] < global_bbox[4*rank+1]+tol)&&
               (v2_bin[i] > global_bbox[4*rank+2]-tol)&&(v2_bin[i] < global_bbox[4*rank+3]+tol) ) {
            if (iter == 0) {
              send_count[rank] += 2;
            }
            else {
              send_buffer_bin[send_disp[rank]/2] = i; // store the sent bin
              send_buffer_double[send_disp[rank]  ] = v1_bin[i];
              send_buffer_double[send_disp[rank]+1] = v2_bin[i];
              send_disp[rank] += 2;
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
        assert( send_count_sum%2 == 0 );
        send_buffer_bin = new int[send_count_sum/2];
        send_buffer_double = new double[send_count_sum];
      }
    }

    // we stored the bin in send_buffer_bin, so no need to keep global bbox's...
    delete[] global_bbox;
    
    // now build the recv side stuff...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    double * recv_buffer_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buffer_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buffer_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

    delete[] send_buffer_double;
      
    // now split everything by 2...
    recv_count_sum /= 2;
    send_count_sum /= 2;
    for (int i = 0; i < mpi_size; i++) {
      recv_count[i] /= 2;
      recv_disp[i] /= 2;
      send_count[i] /= 2;
      send_disp[i] /= 2;
    }

    // will hold a number >= 0 if we got a match in one of our bins, -1 otherwise...
    int * recv_buffer_int = new int[recv_count_sum];
    for (int i = 0; i < recv_count_sum; i++) {
      // following is n^2 in nbins -- sorry. Could do a 
      // second sort and save some time here, but it
      // is just pre-processing, so n^2 with small n is OK?!...
      recv_buffer_int[i] = -1; 
      double this_v1 = recv_buffer_double[2*i];
      double this_v2 = recv_buffer_double[2*i+1];
      for (int ibin = 0; ibin < nbins; ibin++) {
        //if (mpi_rank == 0) 
        //cout << "comparing " << this_v1 << "," << this_v2 << " to " << v1_bin[ibin] << "," << v2_bin[ibin] << endl;
        if ( (fabs(this_v1-v1_bin[ibin]) < tol)&&(fabs(this_v2-v2_bin[ibin]) < tol) ) {
          //if (mpi_rank == 0) 
          //cout << "got match" << endl;
          recv_buffer_int[i] = ibin;
          break;
        }
      }
    }

    delete[] recv_buffer_double;
    
    // send this matching info back to the sending processes...
    
    int * send_buffer_int = new int[send_count_sum];
    MPI_Alltoallv(recv_buffer_int,recv_count,recv_disp,MPI_INT,
                  send_buffer_int,send_count,send_disp,MPI_INT,mpi_comm);

    // now recv and send exchange roles. During averaging, recv will actually 
    // collect its various bins and send them to the send processes for addition
    // to their bin sums, then normalization to produce averaging...

    // ---------------------------------------------
    // the recv side becomes the "pack" side...
    // ---------------------------------------------
    pack_count = new int[mpi_size];
    pack_disp = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++)
      pack_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) {
        for (int i = recv_disp[rank]; i < recv_disp[rank]+recv_count[rank]; i++) {
          if (recv_buffer_int[i] >= 0) {
            if (iter == 0) {              
              pack_count[rank] += 1;
            }
            else {
              pack_index[pack_disp[rank]] = recv_buffer_int[i];
              pack_disp[rank] += 1;
            }
          }
        }
      }
      // calculate pack_disp on both iterations...
      pack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
        pack_disp[i] = pack_count[i-1] + pack_disp[i-1];
      // on the first time through, allocate the pack stuff...
      if (iter == 0) {
        npack = pack_disp[mpi_size-1] + pack_count[mpi_size-1];
        pack_index = new int[npack];
        pack_buffer = new double[npack];
      }
    }

    delete[] recv_count;
    delete[] recv_disp;
    delete[] recv_buffer_int;
    
    // ---------------------------------------------
    // and the send side becomes the unpack side...
    // ---------------------------------------------
    unpack_count = new int[mpi_size];
    unpack_disp = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++)
      unpack_count[i] = 0;
    for (int iter = 0; iter < 2; iter++) {
      for (int rank = 0; rank < mpi_size; rank++) {
        for (int i = send_disp[rank]; i < send_disp[rank]+send_count[rank]; i++) {
          if (send_buffer_int[i] >= 0) {
            if (iter == 0) {              
              unpack_count[rank] += 1;
            }
            else {
              unpack_index[unpack_disp[rank]] = send_buffer_bin[i]; // recall we stored above
              unpack_disp[rank] += 1;
            }
          }
        }
      }
      // calculate pack_disp on both iterations...
      unpack_disp[0] = 0;
      for (int i = 1; i < mpi_size; i++)
        unpack_disp[i] = unpack_count[i-1] + unpack_disp[i-1];
      // on the first time through, allocate the unpack stuff...
      if (iter == 0) {
        nunpack = unpack_disp[mpi_size-1] + unpack_count[mpi_size-1];
        unpack_index = new int[nunpack];
        unpack_buffer = new double[nunpack];
      }
    }

    delete[] send_count;
    delete[] send_disp;
    delete[] send_buffer_bin;
    delete[] send_buffer_int;
    
    // ---------------------------------------------------------
    // now check bin correspondence by exchanging v1 and v2...
    // ---------------------------------------------------------

    // v1...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = v1_bin[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      assert( fabs(unpack_buffer[i] - v1_bin[unpack_index[i]]) < tol );
    
    // v2...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = v2_bin[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      assert( fabs(unpack_buffer[i] - v2_bin[unpack_index[i]]) < tol );
    
    // finally sum the bin_weight's...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = bin_weight[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      bin_weight[unpack_index[i]] += unpack_buffer[i];

    // final stuff...
    assert( bin_value == NULL );
    bin_value = new double[nbins];
    
  }

 public:

  void apply(double * var) {
    
    for (int ibin = 0; ibin < nbins; ibin++)
      bin_value[ibin] = 0.0;
    
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      double weight = weight_of_index[i];
      bin_value[ibin] += weight*var[i];
    }
    
    // reduce values...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = bin_value[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      bin_value[unpack_index[i]] += unpack_buffer[i];
    
    for (int ibin = 0; ibin < nbins; ibin++)
      bin_value[ibin] /= bin_weight[ibin];
    
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      var[i] = bin_value[ibin];
    }
    
  }
  
  void applySum(double * var) {
    
    for (int ibin = 0; ibin < nbins; ibin++)
      bin_value[ibin] = 0.0;
    
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      double weight = weight_of_index[i];
      bin_value[ibin] += weight*var[i];
    }
    
    // reduce values...
    for (int i = 0; i < npack; i++)
      pack_buffer[i] = bin_value[pack_index[i]];
    MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                  unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
    for (int i = 0; i < nunpack; i++)
      bin_value[unpack_index[i]] += unpack_buffer[i];
    
    // remove this normalization...
    //for (int ibin = 0; ibin < nbins; ibin++)
    //  bin_value[ibin] /= bin_weight[ibin];
    
    for (int i = 0; i < n; i++) {
      int ibin = bin_of_index[i];
      var[i] = bin_value[ibin];
    }
    
  }
  
  void apply(double (*var)[3]) {
    
    for (int j = 0; j < 3; j++) {
      
      for (int ibin = 0; ibin < nbins; ibin++)
        bin_value[ibin] = 0.0;
      
      for (int i = 0; i < n; i++) {
        int ibin = bin_of_index[i];
        double weight = weight_of_index[i];
        bin_value[ibin] += weight*var[i][j];
      }
    
      // reduce values...
      for (int i = 0; i < npack; i++)
        pack_buffer[i] = bin_value[pack_index[i]];
      MPI_Alltoallv(pack_buffer,pack_count,pack_disp,MPI_DOUBLE,
                    unpack_buffer,unpack_count,unpack_disp,MPI_DOUBLE,mpi_comm);
      for (int i = 0; i < nunpack; i++)
        bin_value[unpack_index[i]] += unpack_buffer[i];
      
      for (int ibin = 0; ibin < nbins; ibin++)
        bin_value[ibin] /= bin_weight[ibin];
    
      for (int i = 0; i < n; i++) {
        int ibin = bin_of_index[i];
        var[i][j] = bin_value[ibin];
      }
      
    }
    
  }

};

#endif

