#ifndef MPI_TOMMIE_H
#define MPI_TOMMIE_H

#include <stdio.h>
#include <iostream>

/// dummy MPI defines
/// defines are size of types !?!?!?!?!
#define MPI_CHAR             1
#define MPI_INT              4
#define MPI_BYTE             1
#define MPI_INTEGER          4
#define MPI_FLOAT            4
#define MPI_DOUBLE           8
#define MPI_DOUBLE_PRECISION 8
#define MPI_LONG             8
#define MPI_LONG_LONG        8

#define MPI_MAX 1
#define MPI_SUM 2
#define MPI_MIN 3

#define MPI_COMM_WORLD 1
#define MPI_COMM_SELF  1
#define MPI_INFO_NULL 1
#define MPI_COMM_NULL -1

#define MPI_MODE_RDONLY              2  // ADIO_RDONLY
#define MPI_MODE_RDWR                8  // ADIO_RDWR
#define MPI_MODE_WRONLY              4  // ADIO_WRONLY
#define MPI_MODE_CREATE              1  // ADIO_CREATE
#define MPI_MODE_EXCL               64  // ADIO_EXCL
#define MPI_MODE_DELETE_ON_CLOSE    16  // ADIO_DELETE_ON_CLOSE
#define MPI_MODE_UNIQUE_OPEN        32  // ADIO_UNIQUE_OPEN
#define MPI_MODE_APPEND            128  // ADIO_APPEND
#define MPI_MODE_SEQUENTIAL        256  // ADIO_SEQUENTIAL
#define MPI_SEEK_SET  1
#define MPI_SEEK_CUR  1

/// typedefs
typedef int MPI_Comm;
typedef int MPI_Status;
typedef int MPI_Request;
typedef int MPI_Datatype;
typedef int MPI_Info;
typedef long int MPI_Aint;
typedef int MPI_Op;

#ifdef MPI_OFFSET_IS_LONG_LONG_INT
typedef long long int MPI_Offset;
#else
typedef long int MPI_Offset;
#endif

typedef FILE * MPI_File;

int MPI_Init(int * argc, char *** argv);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm * mycomm);
int MPI_Comm_split(MPI_Comm comm,int color,int key,MPI_Comm * comm_out);
int MPI_Comm_rank(MPI_Comm comm, int * mpi_rank);
int MPI_Comm_size(MPI_Comm mycomm, int * mpi_size);
int MPI_Barrier(MPI_Comm comm);
int MPI_Finalize();
int MPI_Waitall(int count, MPI_Request *requests, MPI_Status *status);
int MPI_Waitany(int count, MPI_Request *requests, int *index, MPI_Status *status);

int MPI_File_open(MPI_Comm comm, char *fname, int mode, int info, MPI_File *fp);
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
    char *datarep, int info);
int MPI_File_close(MPI_File *fh);
int MPI_File_delete(char *fname,MPI_Info info);

int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Type_indexed(int count, int *blocklens, int *indices, MPI_Datatype old_type,
    MPI_Datatype *newtype);
int MPI_Type_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type,
    MPI_Datatype *newtype);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types,
    MPI_Datatype *newtype);
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent);
int MPI_Type_size(MPI_Datatype datatype, int *size);
int MPI_File_seek(MPI_File fh, MPI_Offset offset, int status);
int MPI_Address(void *location, MPI_Aint *address);

double MPI_Wtime();

int MPI_File_set_size(MPI_File fh,MPI_Offset size);

//
// implementation of DUMMY-MPI functions when compiling without mpi,
// function templates have to be implemented in the header file
//

///
/// MPI imitator fuction is using C style for io stream
/// instead of using sizeof(T) for size specification in fread -
/// had to use DEFINE of datatype, see functions 'read_no_vector' and 'read_cvofa' ...
template<class T>
int MPI_File_read_all(MPI_File fp, T *buf, int count, MPI_Datatype datatype,
    MPI_Status *status) {
  fread(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_read(MPI_File fp, T *buf, int count, MPI_Datatype datatype,
    MPI_Status *status) {
  fread(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write_all(MPI_File fp, T *buf, int count, MPI_Datatype datatype,
    MPI_Status *status) {
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write(MPI_File fp, T *buf, int count, MPI_Datatype datatype,
		   MPI_Status *status) {
  fwrite(buf, datatype, count, fp);
  return 0;
}

//MPI_File_write_at(fh,offset,ibuf,2,MPI_INT,&status);

template<class T>
int MPI_File_write_at(MPI_File fp, MPI_Offset offset, T *buf, int count, MPI_Datatype datatype,
		      MPI_Status *status) {
  MPI_File_seek(fp,offset,*status);
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_Reduce(T* local, T * global, int n, int type, int action, int rank, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    global[i] = local[i];
  return (0);
}

template<class T>
int MPI_Reduce(T (*local)[3], T (*global)[3], int n, int type, int action, int rank, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      global[i][j] = local[i][j];
  return (0);
}

template<class T>
int MPI_Allreduce(T *local, T *global, int n, int type, int action, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    global[i] = local[i];
  return (0);
}

template<class T>
int MPI_Send(T *data, int count, int type, int rank, int tag, MPI_Comm comm) {
  //    std::cerr << "Error: should not be sending"<< std::endl;
  //    throw(-1);
  return 0;
}

template<class T>
int MPI_Recv(T *data, int count, MPI_Datatype datatype, int rank, int tag,
    MPI_Comm comm, MPI_Status * status) {
  //    std::cerr << "Error: should not be receiving"<< std::endl;
  //    throw(-1);
  return 0;
}

template<class T>
int MPI_Sendrecv(T * sbuf, int ns, int type_s, int rank_s, int tag_s, T * rbuf,
    int nr, int type_r, int rank_r, int tag_r, MPI_Comm
    comm, MPI_Status * status) {
  //    std::cerr << "Error: should not be send/receiving"<< std::endl;
  //    throw(-1);
  return 0;
}

template<class T>
int MPI_Alltoall(T *sendbuf, int sendcount, MPI_Datatype sendtype, T *recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  for (int i = 0; i < sendcount; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Alltoallv(T *sendbuf, int *sendcnts, int *sdispls, MPI_Datatype sendtype,
		  T *recvbuf, int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
  for (int i = 0; i < sendcnts[0]; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Alltoallv(T (*sendbuf)[3], int *sendcnts, int *sdispls, MPI_Datatype sendtype,
		  T (*recvbuf)[3], int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
  //assert( sendcnts[0]%3 == 0 );
  for (int i = 0; i < sendcnts[0]/3; i++) {
    recvbuf[i][0] = sendbuf[i][0];
    recvbuf[i][1] = sendbuf[i][1];
    recvbuf[i][2] = sendbuf[i][2];
  }
  return (0);
}

template<class T>
int MPI_Gather(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
               T *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm) {
  //assert(sendcnt == recvcount);
  //assert(root == 0);
  for (int i = 0; i < sendcnt; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Gatherv(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
                T *recvbuf, int *recvcnts, int *displs,
                MPI_Datatype recvtype,
                int root, MPI_Comm comm) {
  //assert(sendcnt == recvcnts[0]);
  //assert(root == 0);
  for (int i = 0; i < sendcnt; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Allgatherv(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
		   T *recvbuf, int *recvcnts, int *displs,
		   MPI_Datatype recvtype,
		   MPI_Comm comm) {
  //assert(sendcnt == recvcnts[0]);
  //assert(root == 0);
  for (int i = 0; i < sendcnt; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Allgather(T * buf1, int n1, int type1, T * buf2, int n2, int type2,
    MPI_Comm comm) {
  for (int i = 0; i < n1; i++)
    buf2[i] = buf1[i];
  return (0);
}

template<class T>
int MPI_Allgather(T (buf1)[3], int n1, int type1, T (* buf2)[3], int n2, int type2,
    MPI_Comm comm) {
  //assert( n1%3 == 0 );
  for (int i = 0; i < n1/3; i++) {
    buf2[i][0] = buf1[0];
    buf2[i][1] = buf1[1];
    buf2[i][2] = buf1[2];
  }
  return (0);
}

template<class T>
int MPI_Scatterv(T *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype,
                 T *recvbuf, int recvcnt, MPI_Datatype recvtype,
                 int root, MPI_Comm comm) {
  //assert(sendcnts[0] == recvcnt);
  //assert(root == 0);
  for (int i = 0; i < recvcnt; i++)
    recvbuf[i] = sendbuf[i];
  return (0);
}

template<class T>
int MPI_Scan(T * buf1, T * buf2, int n, int type, int action, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    buf2[i] = buf1[i];
  return (0);
}

template<class T>
int MPI_Bcast(T * data, int n, int type, int rank, MPI_Comm comm) {
  return (0);
}

template<class T>
int MPI_Issend(T *data, int count, MPI_Datatype datatype, int rank, int tag,
    MPI_Comm comm, MPI_Request * request) {
  return (0);
}

template<class T>
int MPI_Irecv(T *data, int count, MPI_Datatype datatype, int rank, int tag,
    MPI_Comm comm, MPI_Request * request) {
  return (0);
}

#endif


