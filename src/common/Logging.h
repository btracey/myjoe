#ifndef CTR_TOMMIETOO_LOGGING_H
#define CTR_TOMMIETOO_LOGGING_H

#include <string>

namespace logging
{
  enum LoggingLevel{
    EVERYTHING =   0,
    DEBUG_LO   =  10,
    DEBUG_HI   =  20,
    INFO_LO    =  30,
    INFO       =  40,
    INFO_HI    =  50,
    WARNING    =  60,
    ERROR      =  70,
    CRITICAL   =  80,
    FATAL      =  90,
    SILENT     = 100
  };

  std::ostream& lout(int level);
  std::ostream& lerr(int level);

  void setLoggingLevel(int level);
}

#endif

