#ifndef MPISTUFF_H
#define MPISTUFF_H

#ifdef NO_MPI
# include "mpi_tommie.h"
#else
# include <mpi.h>
#endif

#ifdef MPI_OFFSET_IS_LONG_LONG_INT
# define MPI_OFFSET_DATATYPE MPI_LONG_LONG
# define MPI_INT8 MPI_LONG_LONG
typedef long long int int8;
#else
# define MPI_OFFSET_DATATYPE MPI_LONG
# define MPI_INT8 MPI_LONG
typedef long int int8;
#endif

namespace MpiStuff {

  // these namespace members are declared "extern" so we can include
  // the namespace definition in a header file included by
  // multiple routines
  
  extern int mpi_rank;
  extern int mpi_size;
  extern MPI_Comm mpi_comm;

  // call this method if you really are running mpi...
  extern void initMpiStuff();
  extern void initMpiStuff(MPI_Comm& comm);
  extern void MPI_Pause(char * message);
  extern void MPI_Type_indexed_clean(int n_zones, int *my_zone_count, int *my_zone_disp, MPI_Datatype old_type, MPI_Datatype *new_type);
  
};

#endif
