#include "MshFilter.h"

int MshFilter::getNextToken(char * token) {

  int token_pos = 0;
  int token_level = level;
  int quote_mode = 0;

  while (1) {

    if (pos >= max_pos) {
      max_pos = fread(buf, sizeof(char), MSH_BUFFER_SIZE, fp);
      if (max_pos == 0) {
        if (token_pos > 0) {
          break;
        } else {
          return (-1);
        }
      }
      pos = 0;
    }

    char c = buf[pos++];
    if (c == '"') {
      if (quote_mode == 1) {
        break;
      }
      quote_mode = 1;
      continue;
    } else if (quote_mode == 1) {
      if (token_pos >= MSH_TOKEN_SIZE) {
        cerr << "Error: increase MSH_TOKEN_SIZE."<< endl;
        throw(-1);
      }
      token[token_pos++] = c;
    } else if (c == ' '|| c == '\t'|| c == '\n') {
      if (token_pos > 0) {
        break;
      } else {
        continue;
      }
    } else if (c == '(') {
      level++;
      if (token_pos > 0) {
        break;
      } else {
        token_level++;
        continue;
      }
    } else if (c == ')') {
      level--;
      if (token_pos > 0) {
        break;
      } else {
        token_level--;
        continue;
      }
    } else {
      if (token_pos >= MSH_TOKEN_SIZE) {
        cerr << "Error: increase MSH_TOKEN_SIZE."<< endl;
      }
      token[token_pos++] = c;
    }

  }

  token[token_pos] = '\0';
  return (token_level);

}

void MshFilter::processDimension() {

  // next token should be level 1 also...
  if (getNextToken(token) != 1) {
    cerr << "Error: expect level == 1."<< endl;
    throw(-1);
  }

  // require 3D mesh...
  if (strcmp(token, "3") != 0) {
    cerr << "Error: dimension != 3: "<< token << endl;
    throw(-1);
  }

}

void MshFilter::processNodes() {

  if (getNextToken(token) != 2) {
    cerr << "Error: expect level == 2"<< endl;
    throw(-1);
  }

  if (strcmp(token, "0") == 0) {

    // node count...
    getNextToken(token);
    int istart = atox(token);
    assert(istart == 1);

    getNextToken(token);
    int iend = atox(token);

    // the global node count is iend...
    // we will distribute the reading of the nodes across the processors...

    // just lump all nodes together for now...
    int nno_global = iend;
    if (mpi_rank == 0)
      cout << "nno_global: "<< nno_global << endl;

    // nodes are striped as evenly as possible across all processors...
    assert(ugp->noora == NULL);
    ugp->noora = new int[mpi_size + 1];
    calcUniformDist(ugp->noora, nno_global, mpi_size);

    assert(ugp->nno == 0);
    ugp->nno = ugp->noora[mpi_rank+1] - ugp->noora[mpi_rank];

    // the no_flag is not used for anything yet, but must be
    // allocated as part of the minimum initialization... 
    assert(ugp->no_flag == NULL );
    ugp->no_flag = new int[ugp->nno];
    for (int ino = 0; ino < ugp->nno; ino++)
      ugp->no_flag[ino] = 0;

    // the nodal coordinates are not registered data - 
    // allocate them mannually...   
    assert(ugp->x_no == NULL );
    ugp->x_no = new double[ugp->nno][3];

    // allocate any other registered nodal data...
    ugp->allocateRegisteredNoData();

  } else {

    int zone_index = atox(token);
    assert(zone_index != 0);

    // next are the global indices of this node range...

    getNextToken(token);
    int istart = atox(token) - 1;

    getNextToken(token);
    int iend = atox(token) - 1;

    // only parse the remainder if we own nodes in this range...
    if ((istart < ugp->noora[mpi_rank+1])&&(iend >= ugp->noora[mpi_rank])) {

      getNextToken(token);
      if ((strcmp(token, "1") != 0)&&(strcmp(token, "2") != 0)) {
        cerr << "Error: expecting 1 or 2"<< endl;
        throw(-1);
      }

      getNextToken(token);
      if (strcmp(token, "3") != 0) {
        cerr << "Error: expecting 3"<< endl;
        throw(-1);
      }

      assert(ugp->no_flag != NULL );
      assert(ugp->x_no != NULL );

      for (int ino_global = istart; ino_global <= iend; ino_global++) {
        int ino = ino_global - ugp->noora[mpi_rank];
        if ((ino >= 0)&&(ino < ugp->nno)) {
          assert(ugp->no_flag[ino] == 0);
          ugp->no_flag[ino] = zone_index;
          for (int i = 0; i < 3; i++) {
            // parse 3 doubles...
            getNextToken(token);
            ugp->x_no[ino][i] = atod(token);
          }
        } else {
          // advance 3 spaces...
          for (int i = 0; i < 3; i++)
            getNextToken(token);
        }
      }

    }

  }

}

void MshFilter::processProjectedNodes() {

  DoubleVector * x_no_project = ugp->getVectorData("X_NO_PROJECT");
  if (x_no_project == NULL) {
    if (mpi_rank == 0)
      cout
          << "x_no_project data will be ignored: vector data X_NO_PROJECT not registered."
          << endl;
    return;
  }
  if (mpi_rank == 0)
    cout << "reading x_no_project..."<< endl;

  assert(x_no_project->getDatatype() == NO_DATA);

  if (getNextToken(token) != 2) {
    cerr << "Error: expect level == 2"<< endl;
    throw(-1);
  }

  if (strcmp(token, "0") == 0) {

    cerr << "Error: unexpected token in projected nodes."<< endl;
    throw(-1);

  } else {

    int zone_index = atox(token);
    assert(zone_index != 0);

    // next are the global indices of this node range...

    getNextToken(token);
    int istart = atox(token) - 1;

    getNextToken(token);
    int iend = atox(token) - 1;

    // only parse the remainder if we own nodes in this range...
    if ((istart < ugp->noora[mpi_rank+1])&&(iend >= ugp->noora[mpi_rank])) {

      getNextToken(token);
      if ((strcmp(token, "1") != 0)&&(strcmp(token, "2") != 0)) {
        cerr << "Error: expecting 1 or 2"<< endl;
        throw(-1);
      }

      getNextToken(token);
      if (strcmp(token, "3") != 0) {
        cerr << "Error: expecting 3"<< endl;
        throw(-1);
      }

      for (int ino_global = istart; ino_global <= iend; ino_global++) {
        int ino = ino_global - ugp->noora[mpi_rank];
        if ((ino >= 0)&&(ino < ugp->nno)) {
          for (int i = 0; i < 3; i++) {
            // parse 3 doubles...
            getNextToken(token);
            (*x_no_project->ptr)[ino][i] = atod(token);
          }
        } else {
          // advance 3 spaces...
          for (int i = 0; i < 3; i++)
            getNextToken(token);
        }
      }

    }

  }

}

void MshFilter::processCvs() {

  if (getNextToken(token) != 2) {
    cerr << "Error: expect level == 2"<< endl;
    throw(-1);
  }

  if (strcmp(token, "0") == 0) {

    // cell count...
    getNextToken(token);
    int istart = atox(token);
    assert(istart == 1);

    getNextToken(token);
    int iend = atox(token);

    int ncv_global = iend;
    if (mpi_rank == 0)
      cout << "ncv_global: "<< ncv_global << endl;

    assert(ugp->cvora == NULL);
    ugp->cvora = new int[mpi_size + 1];
    calcUniformDist(ugp->cvora, ncv_global, mpi_size);

    assert(ugp->ncv == 0);
    ugp->ncv = ugp->cvora[mpi_rank+1] - ugp->cvora[mpi_rank];

    // ugp->cv_flag is used to store the cv zone index...
    assert(ugp->cv_flag == NULL );
    ugp->cv_flag = new int[ugp->ncv];
    for (int icv = 0; icv < ugp->ncv; icv++)
      ugp->cv_flag[icv] = 0;

    // allocate any other registered cv data...
    ugp->allocateRegisteredCvData();

  } else {

    int zone_index = atox(token);
    assert(zone_index != 0);

    getNextToken(token);
    int istart = atox(token) - 1; // zero indexing

    getNextToken(token);
    int iend = atox(token) - 1; // zero indexing

    ugp->cvZoneList.push_back(CvZone());
    CvZone * zone = &(ugp->cvZoneList.back());
    zone->setIndex(zone_index);
    zone->setKind(CV_ZONE_UNKNOWN);

    // only parse the rest if it contains cv's we own...
    if ((istart < ugp->cvora[mpi_rank+1])&&(iend >= ugp->cvora[mpi_rank])) {

      getNextToken(token);
      int celltype = atox(token);

      // put the associated zone number into the cv_flag...
      assert(ugp->cv_flag != NULL );

      getNextToken(token);
      int dummy = atox(token);
      if (dummy == 0) {
        // these cells listed individually for some reason?...
        for (int icv_global = istart; icv_global <= iend; icv_global++) {
          getNextToken(token);
          int icv = icv_global - ugp->cvora[mpi_rank];
          if ((icv >= 0)&&(icv < ugp->ncv)) {
            assert(ugp->cv_flag[icv] == 0);
            ugp->cv_flag[icv] = zone_index;
          }
        }
      } else {
        for (int icv_global = istart; icv_global <= iend; icv_global++) {
          int icv = icv_global - ugp->cvora[mpi_rank];
          if ((icv >= 0)&&(icv < ugp->ncv)) {
            assert(ugp->cv_flag[icv] == 0);
            ugp->cv_flag[icv] = zone_index;
          }
        }
      }

    }

  }

}

void MshFilter::processFaces() {

  if (getNextToken(token) != 2) {
    cerr << "Error: expect level == 2"<< endl;
    throw(-1);
  }

  if (strcmp(token, "0") == 0) {

    // face count...
    getNextToken(token);
    int istart = atox(token);
    assert(istart == 1);

    getNextToken(token);
    int iend = atox(token);

    // faces are striped as evenly as possible across all processors...
    int nfa_global = iend;
    if (mpi_rank == 0)
      cout << "nfa_global: "<< nfa_global << endl;

    assert(ugp->faora == NULL);
    ugp->faora = new int[mpi_size + 1];
    calcUniformDist(ugp->faora, nfa_global, mpi_size);

    assert( (ugp->nfa) == 0);
    ugp->nfa = ugp->faora[mpi_rank+1] - ugp->faora[mpi_rank];

    // fa_flag is used to store the face zone index...
    assert(ugp->fa_flag == NULL );
    ugp->fa_flag = new int[(ugp->nfa)];
    for (int ifa = 0; ifa < (ugp->nfa); ifa++)
      ugp->fa_flag[ifa] = 0;

    // faced-based connectivity structures...
    ugp->noofa_i = new int[ugp->nfa+1];
    ugp->cvofa = new int[ugp->nfa][2];
    noofa_v_tmp = new int[ugp->nfa][MAX_NO_PER_FA];

  } else {

    int zone_index = atox(token);
    assert(zone_index != 0);

    getNextToken(token);
    int istart = atox(token) - 1;

    getNextToken(token);
    int iend = atox(token) - 1;

    //cout << "face zone " << zone_index << " nfa: " << iend-istart+1 << endl;

    // everybody adds the face zones...
    ugp->faZoneList.push_back(FaZone());
    FaZone * zone = &(ugp->faZoneList.back());
    zone->setIndex(zone_index);
    zone->setKind(FA_ZONE_UNKNOWN);

    // only parse the rest if you need to...
    if ((istart < ugp->faora[mpi_rank+1])&&(iend >= ugp->faora[mpi_rank])) {

      getNextToken(token);
      int facetype = atox(token);

      getNextToken(token);
      int nno_per_face = atox(token);

      // support 6.3 polygons...

      if ((nno_per_face == 0)||(nno_per_face == 5)) {

        // this version of the faces has the number
        // of nodes on each line...
        for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

          getNextToken(token);
          nno_per_face = atox(token);
          assert(nno_per_face <= MAX_NO_PER_FA);

          int no_list[MAX_NO_PER_FA];
          for (int i = 0; i < nno_per_face; i++) {
            getNextToken(token);
            no_list[i] = atox(token) - 1; // use zero-indexing
          }

          getNextToken(token);
          int icv0 = atox(token) - 1;

          getNextToken(token);
          int icv1 = atox(token) - 1;

          // get the local face indexing...
          int ifa = ifa_global - ugp->faora[mpi_rank];
          if ((ifa >= 0)&&(ifa < ugp->nfa)) {

            assert(ugp->fa_flag[ifa] == 0);
            ugp->fa_flag[ifa] = zone_index;

            // flip the face right away so that any face with a -1 cv
            // has it in the second position, i.e. cvofa[ifa][1]...
            if (icv0 >= 0) {
              // put no count in ifa+1 position for CSR struct 
              ugp->noofa_i[ifa+1] = nno_per_face;
              // reverse the order of the nodes - msh has left-handed rule...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_tmp[ifa][i] = no_list[nno_per_face-i-1];
              ugp->cvofa[ifa][0] = icv0;
              ugp->cvofa[ifa][1] = icv1;
            } else {
              assert(icv1 >= 0);
              // put no count in ifa+1 position for CSR struct 
              ugp->noofa_i[ifa+1] = nno_per_face;
              // use the current face ordering...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_tmp[ifa][i] = no_list[i];
              ugp->cvofa[ifa][0] = icv1;
              ugp->cvofa[ifa][1] = icv0;
            }

          }

        }

      } else {

        // these faces have the same number of nodes...
        assert(nno_per_face <= MAX_NO_PER_FA);

        // this version of the faces has the number
        // of nodes on each line...
        for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

          int no_list[MAX_NO_PER_FA];
          for (int i = 0; i < nno_per_face; i++) {
            getNextToken(token);
            no_list[i] = atox(token) - 1; // use zero-indexing
          }

          getNextToken(token);
          int icv0 = atox(token) - 1;

          getNextToken(token);
          int icv1 = atox(token) - 1;

          // get the local face indexing...
          int ifa = ifa_global - ugp->faora[mpi_rank];
          if ((ifa >= 0)&&(ifa < ugp->nfa)) {

            assert(ugp->fa_flag[ifa] == 0);
            ugp->fa_flag[ifa] = zone_index;

            // flip the face right away so that any face with a -1 cv
            // has it in the second position, i.e. cvofa[ifa][1]...
            if (icv0 >= 0) {
              // put no count in ifa+1 position for CSR struct 
              ugp->noofa_i[ifa+1] = nno_per_face;
              // reverse the order of the nodes - msh has left-handed rule...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_tmp[ifa][i] = no_list[nno_per_face-i-1];
              ugp->cvofa[ifa][0] = icv0;
              ugp->cvofa[ifa][1] = icv1;
            } else {
              assert(icv1 >= 0);
              // put no count in ifa+1 position for CSR struct 
              ugp->noofa_i[ifa+1] = nno_per_face;
              // use the current face ordering...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_tmp[ifa][i] = no_list[i];
              ugp->cvofa[ifa][0] = icv1;
              ugp->cvofa[ifa][1] = icv0;
            }

          }

        }

      }

    }

  }

}

void MshFilter::processPeriodic() {

  getNextToken(token);
  int istart = atox(token) - 1;

  getNextToken(token);
  int iend = atox(token) - 1;

  getNextToken(token);
  int zone1_index = atox(token);

  getNextToken(token);
  int zone2_index = atox(token);

  // get these zones in the zone list...
  int found1 = 0;
  int found2 = 0;
  for (list<FaZone>::iterator fzone = ugp->faZoneList.begin(); fzone
      != ugp->faZoneList.end(); fzone++) {
    if (fzone->getIndex() == zone1_index) {
      assert(found1 == 0);
      found1 = 1;
      assert(fzone->getKind() == FA_ZONE_UNKNOWN);
      fzone->setKind(FA_ZONE_PERIODIC_UNKNOWN);
      assert(fzone->getPeriodicIndex() == 0);
      fzone->setPeriodicIndex(zone2_index );
    }
    if (fzone->getIndex() == zone2_index) {
      assert(found2 == 0);
      found2 = 1;
      assert(fzone->getKind() == FA_ZONE_UNKNOWN);
      fzone->setKind(FA_ZONE_PERIODIC_UNKNOWN); // assume Cartesian periodicity for now, but this gets cleaned up below
      assert(fzone->getPeriodicIndex() == 0);
      fzone->setPeriodicIndex(zone1_index );
    }
  }
  assert( (found1 == 1)&&(found2 == 1));

  /*
   cout << "Setting zones " << zone1_index << " and " << zone2_index << " to periodic." << endl;
   cout << "istart, iend: " << istart << " " << iend << endl;
   */

  // we all need to do this, just in case the face pair is in our local faces...
  for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

    getNextToken(token);
    int ifa1_global = atox(token) - 1;

    getNextToken(token);
    int ifa2_global = atox(token) - 1;

    if ((ifa1_global < ugp->faora[mpi_rank+1])&&(ifa1_global>= ugp->faora[mpi_rank])) {
      int ifa = ifa1_global - ugp->faora[mpi_rank];
      // use cvofa[ifa][1] top store the global face...
      assert(ugp->cvofa[ifa][1] == -1);
      ugp->cvofa[ifa][1] = ifa2_global;
    }

    if ((ifa2_global < ugp->faora[mpi_rank+1])&&(ifa2_global>= ugp->faora[mpi_rank])) {
      int ifa = ifa2_global - ugp->faora[mpi_rank];
      // use cvofa[ifa][1] top store the global face...
      assert(ugp->cvofa[ifa][1] == -1);
      ugp->cvofa[ifa][1] = ifa1_global;
    }

  }

}

void MshFilter::processZones() {

  getNextToken(token);
  int zone_index = atoi(token);
  assert(zone_index != 0);

  // fluent's name...
  getNextToken(token);
  char name1[MSH_TOKEN_SIZE];
  strcpy(name1, token);

  // user name...
  getNextToken(token);
  char name2[MSH_TOKEN_SIZE];
  strcpy(name2, token);

  //if (mpi_rank == 0)
  //  cout << "zone " << zone_index << ", name1: " << name1 << ", name2: " << name2 << endl;

  // try face zones first...
  for (list<FaZone>::iterator fzone = ugp->faZoneList.begin(); fzone
      != ugp->faZoneList.end(); fzone++) {
    if (fzone->getIndex() == zone_index) {
      // we found a match - make sure this zone has not been visited...
      assert(fzone->flag == 0);
      fzone->flag = 1;
      // name the zone based on the name2...
      fzone->setName(name2);
      return;
    }
  }

  // if we are still here, we did not find a metch...
  for (list<CvZone>::iterator czone = ugp->cvZoneList.begin(); czone
      != ugp->cvZoneList.end(); czone++) {
    if (czone->getIndex() == zone_index) {
      // we found a match - make sure this zone has not been visited...
      assert(czone->flag == 0);
      czone->flag = 1;
      // name the zone based on the name2...
      czone->setName(name2);
      return;
    }
  }

  if (mpi_rank == 0)
    cout << "Warning: could not match zone "<< name1 << ":"<< name2
        << " in cv or fa zone list. Skipping."<< endl;

}

void MshFilter::finalizeAndCheck() {

  // all elements should have their flags set...

  for (int ino = 0; ino < ugp->nno; ino++)
    assert(ugp->no_flag[ino] > 0);

  for (int icv = 0; icv < ugp->ncv; icv++)
    assert(ugp->cv_flag[icv] > 0);
  
  // the noofa_i/v CSR structure can now be built...

  int my_ierr = 0;
  ugp->noofa_i[0] = 0;
  for (int ifa = 0; ifa < ugp->nfa; ifa++) {
    //assert(ugp->fa_flag[ifa] > 0);
    if (ugp->fa_flag[ifa] <= 0) {
      my_ierr += 1;
    }
    ugp->noofa_i[ifa+1] += ugp->noofa_i[ifa];
  }
  ugp->noofa_s = ugp->noofa_i[ugp->nfa];
  assert(ugp->noofa_v == NULL );
  ugp->noofa_v = new int[ugp->noofa_s];
  for (int ifa = 0; ifa < ugp->nfa; ifa++) {
    int nof_f = ugp->noofa_i[ifa];
    int nof_l = ugp->noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof <= nof_l; nof++) {
      ugp->noofa_v[nof] = noofa_v_tmp[ifa][nof-nof_f];
    }
  }
  delete[] noofa_v_tmp;

  int ierr;
  MPI_Allreduce(&my_ierr, &ierr, 1, MPI_INT, MPI_SUM, mpi_comm);
  if (ierr > 0) {
    if (mpi_rank == 0)
      cout << "Warning: some faces are not associated with a zone: " << ierr << endl;
  }
  
  //MPI_Pause("XXXX");
  
  // now make sure that all zones have been flagged, and set the zone
  // types based on whether any periodicity has been set previously...

  for (list<FaZone>::iterator fzone = ugp->faZoneList.begin(); 
       fzone != ugp->faZoneList.end(); fzone++) {
    // make sure the zone's name has been set...
    assert(fzone->flag == 1);
    // to set the zone kind, we need to check the faces...
    if (fzone->getKind() == FA_ZONE_UNKNOWN) {
      // if the zone is unknown, then it has not been associated with any periodicity,
      // and if will be either FA_ZONE_BOUNDARY or FA_ZONE_INTERNAL, depending on 
      // whether it has icv1 defined or not. Sometimes a single zone can have both
      // defined and undefined faces, and this should be handled by splitting in 
      // the future. For now, just error:
      int my_zone_kind = FA_ZONE_UNKNOWN;
      for (int ifa = 0; ifa < ugp->nfa; ifa++) {
        if (ugp->fa_flag[ifa] == fzone->getIndex()) {
          // check icv1...
          int icv1 = ugp->cvofa[ifa][1];
          if (icv1 == -1) {
            if (my_zone_kind == FA_ZONE_UNKNOWN) {
              my_zone_kind = FA_ZONE_BOUNDARY;
            } 
	    else if (my_zone_kind != FA_ZONE_BOUNDARY) {
              cout << "Error: zone "<< fzone->getName()
		   << " has both defined and undefined icv1. Implement zone splitting."
		   << endl;
              //throw(-1);
            }
          } 
	  else {
            // make sure icv1 is in the valid global range of cv's...
            assert( (icv1 >= 0)&&(icv1 < ugp->cvora[mpi_size]));
            if (my_zone_kind == FA_ZONE_UNKNOWN) {
              my_zone_kind = FA_ZONE_INTERNAL;
            } 
	    else if (my_zone_kind != FA_ZONE_INTERNAL) {
              cout << "Error: zone "<< fzone->getName()
		   << " has both defined and undefined icv1. Implement zone splitting."
		   << endl;
              //throw(-1);
            }
          }
        }
      }
      // because we do not have all faces, we need to reduce...
      int zone_kind;
      MPI_Allreduce(&my_zone_kind, &zone_kind, 1, MPI_INT, MPI_MAX, mpi_comm);
      if (my_zone_kind != zone_kind) {
        if (my_zone_kind != FA_ZONE_UNKNOWN) {
          cout << "Error: zone "<< fzone->getName()
	       << " has both defined and undefined icv1. Implement zone splitting."<< endl;
          //throw(-1);
        }
      } else if (zone_kind == FA_ZONE_UNKNOWN) {
        cout << "Error: zone "<< fzone->getName() << " has no faces."<< endl;
        //throw(-1);
      }
      fzone->setKind(zone_kind);
      if (mpi_rank == 0)
        fzone->dump();
    } 
    else if (fzone->getKind() == FA_ZONE_PERIODIC_UNKNOWN) {
      // make sure all our faces were set with a global face index in their cvofa position...
      for (int ifa = 0; ifa < ugp->nfa; ifa++) {
        if (ugp->fa_flag[ifa] == fzone->getIndex()) {
          int icv1 = ugp->cvofa[ifa][1];
          assert( (icv1 >= 0)&&(icv1 < ugp->faora[mpi_size]));
        }
      }
      // find our periodic pair. 
      list<FaZone>::iterator fzone2;
      for (fzone2 = ugp->faZoneList.begin(); fzone2 != ugp->faZoneList.end(); fzone2++) {
        if (fzone2->getIndex() == fzone->getPeriodicIndex())
          break;
      }
      // make sure we found one...
      assert(fzone2 != ugp->faZoneList.end() );
      // make sure it is tied to us...
      assert(fzone2->getPeriodicIndex() == fzone->getIndex() );
      // check fzone2...
      /*
       for (int ifa = 0; ifa < ugp->nfa; ifa++) {
       if (ugp->fa_flag[ifa] == fzone2->getIndex()) {
       int icv1 = ugp->cvofa[ifa][1];
       assert( (icv1 >= 0)&&(icv1 < ugp->faora[mpi_size]) );
       }
       }
       */
      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      // for now, clear any periodicity...
      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      fzone->setKind(FA_ZONE_BOUNDARY);
      for (int ifa = 0; ifa < ugp->nfa; ifa++) {
        if (ugp->fa_flag[ifa] == fzone->getIndex()) {
          ugp->cvofa[ifa][1] = -1;
        }
      }
      // and dump...
      if (mpi_rank == 0)
        fzone->dump();
    } 
    else {
      cerr << "Error: unexpected fzone->getKind(): "<< fzone->getKind() << endl;
      throw(-1);
    }

  }

  for (list<CvZone>::iterator czone = ugp->cvZoneList.begin(); czone
      != ugp->cvZoneList.end(); czone++) {
    // make sure the zone has been set...
    assert(czone->flag == 1);
    // assume fluid for now...
    czone->setKind(CV_ZONE_FLUID);
    if (mpi_rank == 0)
      czone->dump();

  }
}

