#include "MpiStuff.h"
#include <iostream>
#include <stdio.h>
#include <string.h>

using std::cout;
using std::endl;

namespace MpiStuff {

// give single processor initial values to these
// in practice, the user should call initMpiStuff
// near the start of there cade, after MPI_Init(*,*)...

int mpi_rank = 0;
int mpi_size = 1;
MPI_Comm mpi_comm = MPI_COMM_WORLD;

void initMpiStuff() {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
}

void initMpiStuff(MPI_Comm& comm) {
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_dup(comm, &mpi_comm);
}

void MPI_Pause(char * message) {
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) {
        cout << message << endl;
        getchar();
    }
    MPI_Barrier(mpi_comm);
}

}

