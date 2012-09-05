#ifndef MSH_FILTER_H
#define MSH_FILTER_H

#include "MiscUtils.h"
using namespace MiscUtils;

#include "Ugp.h"

#define MSH_BUFFER_SIZE      65536
#define MSH_TOKEN_SIZE          63
#define MAX_NO_PER_FA            8

class MshFilter {

private:

  FILE * fp;
  int pos, max_pos, level;
  char buf[MSH_BUFFER_SIZE];
  char token[MSH_TOKEN_SIZE+1];
  int (*noofa_v_tmp)[MAX_NO_PER_FA];
  Ugp * ugp;

public:

  MshFilter(string filename) {

    if (mpi_rank == 0)
      cout << "MshFilter(): filename: "<< filename << endl;

    ugp = NULL;
    noofa_v_tmp = NULL;

    if ( (fp = fopen(filename.c_str(), "rb")) == NULL) {
      cerr << "Error: cannot open msh file: "<< filename << endl;
      throw(-1);
    }

    pos = max_pos = 0;
    level = 0;
  }

  void initUgpMin(Ugp * ugp) {

    this->ugp = ugp;

    int token_level;
    while ((token_level = getNextToken(token)) != -1) {
      if (token_level == 1) {
        //if (mpi_rank == 0)
        //  cout << "level 1 token: " << token << endl;
        if (strcmp(token, "0") == 0) {
          if (getNextToken(token) != 1) {
            cerr << "Error: expect another level 1 token."<< endl;
            throw(-1);
          }
          //cout << "0: " << token << endl;
        } else if (strcmp(token, "2") == 0) {
          processDimension();
        } else if (strcmp(token, "10") == 0) {
          processNodes();
        } else if (strcmp(token, "101") == 0) {
          processProjectedNodes();
        } else if (strcmp(token, "12") == 0) {
          processCvs();
        } else if (strcmp(token, "13") == 0) {
          processFaces();
        } else if (strcmp(token, "18") == 0) {
          processPeriodic();
        } else if ( (strcmp(token, "39") == 0) || (strcmp(token, "45") == 0)) {
          processZones();
        } else {
          if (mpi_rank == 0)
            cout << "Warning: unrecognized level 1 token: \""<< token << "\""<< endl;
        }
      }
    }

    fclose(fp);

    // check that the bracket levels closed properly...
    if (level != 0) {
      cerr << "Error: msh bracket levels did not close: "<< level << endl;
      throw(-1);
    }

    finalizeAndCheck();

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      cout << "MshFilter::initUgpMin(): OK"<< endl;

  }

private:

  int getNextToken(char * token);

  void processDimension();

  void processNodes();

  void processProjectedNodes();

  void processCvs();

  void processFaces();

  void processPeriodic();

  void processZones();

  void finalizeAndCheck();

};

#endif // #ifndef MSHFILTER_H
