#include <fstream>
#include <iostream>

#include "MpiStuff.h"
#include "Logging.h"

#ifndef DEFAULT_LOGGING_LEVEL
#define DEFAULT_LOGGING_LEVEL INFO_HI
#endif

namespace logging
{
  int current_level = DEFAULT_LOGGING_LEVEL;

  std::fstream null_stream;
  
  std::ostream & lout(int level)
  {
    if (MpiStuff::mpi_rank == 0 && current_level <= level)
      return std::cout;
    else return null_stream;
  }

  std::ostream & lerr(int level)
  {
    if (MpiStuff::mpi_rank == 0 && current_level <= level)
      return std::cerr;
    else return null_stream;
  }

  void setLoggingLevel(int level)
  {
    current_level = level;
  }
}

