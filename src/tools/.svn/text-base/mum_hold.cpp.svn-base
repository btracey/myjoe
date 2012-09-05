#include <iostream>
#include "MpiStuff.h"
using namespace MpiStuff;

using std::cout;
using std::cerr;
using std::endl;

int main(int argc,char * argv[]) {
  
  MPI_Init(&argc,&argv);
  
  try {

    initMpiStuff();

    if (mpi_rank == 0)
      cout << "mum_hold" << endl;
    
    int done = 0;
    while (done != 1) {
      // do something costly...
      const double wait_time = MPI_Wtime();
      while ( wait_time + 5.0 > MPI_Wtime() );
      if (mpi_rank == 0) {
	cout << " > checking for killhold..." << endl;
	const int ierr = MPI_File_delete("killhold",MPI_INFO_NULL);
	if (ierr == 0)
	  done = 1;
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
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
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return(0);
  
}

