#include "UgpWithTools.h"

// needed on uP for std::sort...
//#include <algorithm>

void UgpWithTools1::initStats(ParamMap * params) {
  if (mpi_rank==0) cout<<"UgpWithTools1::initStats()"<<endl;

  Param *p;
  if (params->getParam(p, "STATS")) {
    for (int i = 1; i<p->getSize(); i++) {
      DoubleScalar * ds = getScalarData(p->getString(i));
      if (ds!=NULL) {
        scalarStatsList.push_back(ScalarStats(ds));
        ScalarStats * stats = &(scalarStatsList.back());
        // it would be much nicer to use strings here...
        char name[DATA_NAME_LEN];
        switch (ds->getDatatype()) {
        case NO_DATA:
          sprintf(name, "%s_WGT", ds->getName());
          registerValue(stats->weight, name);
          sprintf(name, "%s_AVG", ds->getName());
          registerScalar(stats->avg, name, NO_DATA);
          sprintf(name, "%s_RMS", ds->getName());
          registerScalar(stats->rms, name, NO_DATA);
          break;
        case CV_DATA:
          sprintf(name, "%s_WGT", ds->getName());
          registerValue(stats->weight, name);
          sprintf(name, "%s_AVG", ds->getName());
          registerScalar(stats->avg, name, CV_DATA);
          sprintf(name, "%s_RMS", ds->getName());
          registerScalar(stats->rms, name, CV_DATA);
          break;
        default:
          if (mpi_rank==0) cerr<<"Error: initStats: only CV_DATA and NO_DATA supported in stats for now."<<endl;
          throw(-1);
        }
      }
      else {
        DoubleVector * dv = getVectorData(p->getString(i));
        if (dv!=NULL) {
          vectorStatsList.push_back(VectorStats(dv));
          VectorStats * stats = &(vectorStatsList.back());
          // it would be much nicer to use strings here...
          char name[DATA_NAME_LEN];
          switch (dv->getDatatype()) {
          case NO_DATA:
            sprintf(name, "%s_WGT", dv->getName());
            registerValue(stats->weight, name);
            sprintf(name, "%s_AVG", dv->getName());
            registerVector(stats->avg, name, NO_DATA);
            sprintf(name, "%s_RMS", dv->getName());
            registerVector(stats->rms, name, NO_DATA);
            sprintf(name, "%s_REY", dv->getName());
            registerVector(stats->rey, name, NO_DATA);
            break;
          case CV_DATA:
            sprintf(name, "%s_WGT", dv->getName());
            registerValue(stats->weight, name);
            sprintf(name, "%s_AVG", dv->getName());
            registerVector(stats->avg, name, CV_DATA);
            sprintf(name, "%s_RMS", dv->getName());
            registerVector(stats->rms, name, CV_DATA);
            sprintf(name, "%s_REY", dv->getName());
            registerVector(stats->rey, name, CV_DATA);
            break;
          default:
            if (mpi_rank==0) cerr<<"Error: initStats: only CV_DATA and NO_DATA supported in stats for now."<<endl;
            throw(-1);
          }
        }
        else {
          if (mpi_rank==0) cerr<<"Error: initStats: cannot find registered data matching: "<<p->getString(i)<<endl;
          throw(-1);
        }
      }
    }

    // finally, look for some averaging...
    assert(stats_average_kind==AVERAGE_NONE);

    if (params->getParam(p, "STATS_AVERAGE")) {
      string name = p->getString();
      if (name=="X_Z") {
        if (mpi_rank==0) cout<<" > STATS_AVERAGE: X_Z"<<endl;
        stats_average_kind = AVERAGE_X_Z;
      }
      else if (name=="THETAX") {
        if (mpi_rank==0) cout<<" > STATS_AVERAGE: THETAX"<<endl;
        stats_average_kind = AVERAGE_THETAX;
      }
      else {
        if (mpi_rank==0) cout<<" > Warning: unrecognized STATS_AVERAGE: "<<name<<", skipping."<<endl;
      }
    }
  }
}

void UgpWithTools1::initStatsUQ(ParamMap * params) {
  if (mpi_rank==0) cout<<"UgpWithTools1::initStats()"<<endl;

  Param *p;
  if (params->getParam(p, "STATS")) {
    for (int i = 1; i<p->getSize(); i++) {
      DoubleScalar * ds = getScalarData(p->getString(i));
      if (ds!=NULL) {
        scalarStatsList.push_back(ScalarStats(ds));
        ScalarStats * stats = &(scalarStatsList.back());
        // it would be much nicer to use strings here...
        char name[DATA_NAME_LEN];
        switch (ds->getDatatype()) {
        case NO_DATA:
          sprintf(name, "%s_WGT", ds->getName());
          registerValue(stats->weight, name);
          sprintf(name, "%s_AVG", ds->getName());
          registerScalar(stats->avg, name, NO_DATA);
          sprintf(name, "%s_STD", ds->getName());
          registerScalar(stats->rms, name, NO_DATA);
          break;
        case CV_DATA:
          sprintf(name, "%s_WGT", ds->getName());
          registerValue(stats->weight, name);
          sprintf(name, "%s_AVG", ds->getName());
          registerScalar(stats->avg, name, CV_DATA);
          sprintf(name, "%s_STD", ds->getName());
          registerScalar(stats->rms, name, CV_DATA);
          break;
        default:
          if (mpi_rank==0) cerr<<"Error: initStats: only CV_DATA and NO_DATA supported in stats for now."<<endl;
          throw(-1);
        }
      }
      else {
        DoubleVector * dv = getVectorData(p->getString(i));
        if (dv!=NULL) {
          vectorStatsList.push_back(VectorStats(dv));
          VectorStats * stats = &(vectorStatsList.back());
          // it would be much nicer to use strings here...
          char name[DATA_NAME_LEN];
          switch (dv->getDatatype()) {
          case NO_DATA:
            sprintf(name, "%s_WGT", dv->getName());
            registerValue(stats->weight, name);
            sprintf(name, "%s_AVG", dv->getName());
            registerVector(stats->avg, name, NO_DATA);
            sprintf(name, "%s_STD", dv->getName());
            registerVector(stats->rms, name, NO_DATA);
            break;
          case CV_DATA:
            sprintf(name, "%s_WGT", dv->getName());
            registerValue(stats->weight, name);
            sprintf(name, "%s_AVG", dv->getName());
            registerVector(stats->avg, name, CV_DATA);
            sprintf(name, "%s_STD", dv->getName());
            registerVector(stats->rms, name, CV_DATA);

            break;
          default:
            if (mpi_rank==0) cerr<<"Error: initStats: only CV_DATA and NO_DATA supported in stats for now."<<endl;
            throw(-1);
          }
        }
        else {
          if (mpi_rank==0) cerr<<"Error: initStats: cannot find registered data matching: "<<p->getString(i)<<endl;
          throw(-1);
        }
      }
    }

    // finally, look for some averaging...
    assert(stats_average_kind==AVERAGE_NONE);

    if (params->getParam(p, "STATS_AVERAGE")) {
      string name = p->getString();
      if (name=="X_Z") {
        if (mpi_rank==0) cout<<" > STATS_AVERAGE: X_Z"<<endl;
        stats_average_kind = AVERAGE_X_Z;
      }
      else if (name=="THETAX") {
        if (mpi_rank==0) cout<<" > STATS_AVERAGE: THETAX"<<endl;
        stats_average_kind = AVERAGE_THETAX;
      }
      else {
        if (mpi_rank==0) cout<<" > Warning: unrecognized STATS_AVERAGE: "<<name<<", skipping."<<endl;
      }
    }
  }
}

void UgpWithTools::initWriteData(ParamMap * params) {
  if (mpi_rank==0) cout<<"initWriteData()"<<endl;

  vector<Param> *paramVec;
  if (params->getParam(paramVec, "WRITE_DATA")) {
    for (vector<Param>::iterator ip = paramVec->begin(); ip!=paramVec->end(); ip++)
      writeDataList.push_back(WriteData(&(*ip)));
  }
}

void UgpWithTools::initProbes(ParamMap * params) {
  if (mpi_rank==0) cout<<"initProbes()"<<endl;

  vector<Param> *paramVec;
  if (params->getParam(paramVec, "PROBE")) {
    // probes will need the cv adt...
    useCvAdt();

    // cycle through...
    for (vector<Param>::iterator ip = paramVec->begin(); ip!=paramVec->end(); ip++)
      probeList.push_back(Probe(&(*ip)));

    // make sure names are unique...
    int ierr = 0;
    for (list<Probe>::iterator pr1 = probeList.begin(); pr1!=probeList.end(); pr1++) {
      for (list<Probe>::iterator pr2 = probeList.begin(); pr2!=pr1; pr2++) {
        if (pr1->getName()==pr2->getName()) {
          if (mpi_rank==0) cerr<<"Error: PROBE NAME: "<<pr1->getName()<<" duplicated."<<endl;
          ierr = 1;
        }
      }
    }
    if (ierr==1) throw(-1);
  }
}

