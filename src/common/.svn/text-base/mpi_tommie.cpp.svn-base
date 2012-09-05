#ifdef NO_MPI

#include "stdio.h"    // io file stream in "C" style
#include "mpi_tommie.h"
#include <iostream>

int MPI_Init(int * argc, char *** argv)
{
    return (0);
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm * mycomm)
{
    *mycomm = comm;
    return (0);
}

int MPI_Comm_split(MPI_Comm comm,int color,int key,MPI_Comm * comm_out) 
{
  return MPI_Comm_dup(comm,comm_out);
}

int MPI_Comm_rank(MPI_Comm comm, int * mpi_rank)
{
    *mpi_rank = 0;
    return (0);
}

int MPI_Comm_size(MPI_Comm mycomm, int * mpi_size)
{
    *mpi_size = 1;
    return (0);
}

int MPI_Barrier(MPI_Comm comm)
{
    return (0);
}

int MPI_Finalize()
{
    return (0);
}

int MPI_Waitall(int count, MPI_Request * requests, MPI_Status * statuses)
{
    return (0);
}

int MPI_Waitany(int count, MPI_Request * requests, int *index, MPI_Status * statuses)
{
    return (0);
}

/// MPI imitator fuction is using C style for io stream,
/// opening files in binary mode!!!
int MPI_File_open(MPI_Comm comm, char *fname, int mode, int info, MPI_File *fp)
{
    char cmode[10];
    switch (mode)
    {
    case MPI_MODE_RDONLY:
      sprintf(cmode, "rb");
      break;
    case MPI_MODE_WRONLY:
      sprintf(cmode, "wb");   // only binary
      break;
    case MPI_MODE_WRONLY | MPI_MODE_CREATE:
      sprintf(cmode, "wb");   // only binary
      break;
    default:
      std::cerr << "Error: cannot figure out mode" << std::endl;
      break;
    }

    if ((*fp=fopen(fname, cmode)) == NULL)
    {
        printf("could not open %s\n", fname);
        return -1;
    }

    return 0;
}

/// MPI imitator fuction is using C style for io stream
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype, char *datarep, int info)
{
    fseek (fh ,disp ,SEEK_SET);
    return 0;
}

int MPI_File_close(MPI_File *fp)
{
    fclose(*fp);
    return 0;
}

int MPI_File_delete(char *fname,MPI_Info info) {
  // remove returns 0 if successfull...
  return remove(fname);
}

int MPI_Type_free(MPI_Datatype *datatype)
{
    return 0;
}

int MPI_Type_indexed(int count, int *blocklens, int *indices, MPI_Datatype old_type, MPI_Datatype *newtype)
{
    return 0;
}

int MPI_Type_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type, MPI_Datatype *newtype)
{
    return 0;
}

int MPI_Type_commit(MPI_Datatype *datatype)
{
    return 0;
}

int MPI_Type_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types, MPI_Datatype *newtype)
{
    *newtype = 0;
    for (int i=0; i<count; i++)
        *newtype += blocklens[i]*old_types[i];
    return 0;
}

int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent)
{
    *extent = (long int)datatype;
    return 0;
}

int MPI_Type_size(MPI_Datatype datatype, int *size)
{
    *size = (long int)datatype;
    return 0;
}

double MPI_Wtime() {

  return(0.0);

}

int MPI_File_set_size(MPI_File fh,MPI_Offset size) {

  return 0;
}

int MPI_File_seek(MPI_File fh, MPI_Offset offset, int status)
{
  return 0;
}


#endif


