#include "CdpFilter.h"

void CdpFilter::readCdpRestart() {

  // build the MPI_CdpHeader...

  MPI_Datatype oldtypes[5];
  int blockcounts[5];
  MPI_Aint offsets[5];
  int header_size;

  // the name...
  offsets[0] = 0;
  oldtypes[0] = MPI_CHAR;
  blockcounts[0] = CDP_IO_HEADER_NAME_LEN;

  // the id...
  offsets[1] = offsets[0] + 1*CDP_IO_HEADER_NAME_LEN;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 1;

  // the offset skip...
  offsets[2] = offsets[1] + 4;
#ifdef MPI_OFFSET_IS_LONG_LONG_INT
  oldtypes[2] = MPI_LONG_LONG;
#else
  oldtypes[2] = MPI_LONG;
#endif
  blockcounts[2] = 1;

  // the 16 ints...
  offsets[3] = offsets[2] + 8;
  oldtypes[3] = MPI_INT;
  blockcounts[3] = 16;

  // the 16 reals...
  offsets[4] = offsets[3] + 4*16;
  oldtypes[4] = MPI_DOUBLE;
  blockcounts[4] = 16;

  MPI_Type_struct(5, blockcounts, offsets, oldtypes, &MPI_CdpHeader);
  MPI_Type_commit(&MPI_CdpHeader);
  MPI_Type_extent(MPI_CdpHeader, &header_extent);
  MPI_Type_size(MPI_CdpHeader, &header_size);

  if ((header_extent != header_size) || (header_size != 256)) {
    cerr << "Error: header size/extent problem: "<< endl;
    cout << "header_extent : "<< header_extent << endl;
    cout << "header_size   : "<< header_size << endl;
    throw(-1);
  }

  // reset the ifa_global to the offset of new faces (zero for now)...
  ifa_global = 0;

  // reset the icv_global to the offset of new cvs (zero for now)...
  icv_global = 0;

  // open file...
  char name[200];
  sprintf(name, filename.c_str());
  if (MPI_File_open(mpi_comm, name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != 0) {
    cerr << "Error: could not open restart file: "<< filename << endl;
    throw(-1);
  }

  // first 4 ints are:
  // 0. magic nummber
  // 1. i/o version
  // 2. WP - should be 8 bytes
  // 3. ?

  // assume we do not have to byte_swap...
  byte_swap = 0;

  int itmp[4];
  MPI_Status status;
  MPI_File_read_all(fh, itmp, 4, MPI_INT, &status);

  if (itmp[0] != CDP_IO_MAGIC_NUMBER) {
    byteSwap(itmp, 4);
    if (itmp[0] != CDP_IO_MAGIC_NUMBER) {
      cerr << "Error: restart file does not start as expected. aborting."<< endl;
      throw(-1);
    }
    if (mpi_rank == 0)
      cout << "File requires byte swapping."<< endl;
    byte_swap = 1;
  }

  if (itmp[1] != CDP_IO_VERSION) {
    if (mpi_rank == 0)
      cerr << "Error: io version differs: "<< itmp[1]<< endl;
    throw(-1);
  }

  if (itmp[2] != 8) {
    if (mpi_rank == 0)
      cerr << "Error: WP not 8: "<< itmp[2]<< endl;
    throw(-1);
  }

  // remaining file consists of headers, potentially followed by data...
  // CDP writes the file in 2 "sections" each terminated by CDP_IO_LAST, 
  // so that data registration can be
  // done in between knowing the sizes and reading the data. Because C/C++
  // can store pointers to pointers, this multiple should not be required...

  offset = 4*4;
  int section = 0;
  while (section < 2) {

    MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, &header, 1, MPI_CdpHeader, &status);

    // also, header's string comes from a fortran file, so null terminate....
    nullTerminate(header.name, CDP_IO_HEADER_NAME_LEN);

    if (byte_swap)
      byteSwapCdpHeader(&header, 1);

    /*
     if (mpi_rank == 0)
     {
     cout <<"got header name: |"<<header.name<<"|"<<"\tid: "<<header.id<<"   \tswap:"<<byte_swap<<endl;
     for (int i = 0; i < 4; i++)
     {
     cout << "header.idata["<< i << "]: "<< header.idata[i]<< endl;
     }
     for (int i = 0; i < 16; i++) {
     cout << "header.rdata["<< i << "]: "<< header.rdata[i]<< endl;
     }
     }
     */

    switch (header.id) {
    case CDP_IO_LAST:
      section++;
      break;
    case CDP_IO_NO_COUNTS:
      build_no_structs();
      break;
    case CDP_IO_NO_FLAG_TEST:
      no_flag_test();
      break;
    case CDP_IO_FA_COUNTS:
      build_fa_structs();
      break;
    case CDP_IO_FA_FLAG_TEST:
      fa_flag_test();
      break;
    case CDP_IO_CV_COUNTS:
      build_cv_structs();
      break;
    case CDP_IO_CV_FLAG_TEST:
      cv_flag_test();
      break;
    case CDP_IO_CV_PART:
      read_cv_part();
      break;
    case CDP_IO_FA_ZONE:
      add_face_zone();
      break;
    case CDP_IO_NOOFA_AND_EDOFA:
      read_noofa_iv();
      break;
    case CDP_IO_CVOFA:
      read_cvofa();
      break;
    case CDP_IO_NODE_CC:
      read_no_vector(ugp->x_no);
      dumpVectorRange(ugp->x_no, ugp->nno, "X_NO");
      break;
    case CDP_IO_NO_R1: 
      {
        DoubleScalar * data = ugp->getScalarData(header.name);
        if ((data != NULL)&&(data->getDatatype() == NO_DATA)) {
          read_no_scalar(*(data->ptr));
          dumpScalarRange(*(data->ptr), ugp->nno, header.name);
        } 
        else {
          if (mpi_rank == 0)
            cout << "skipping NO_DATA scalar: "<< header.name<< endl;
        }
      }
      break;

    case CDP_IO_NO_R2: 
      {
        DoubleVector * data = ugp->getVectorData(header.name);
        if ((data != NULL)&&(data->getDatatype() == NO_DATA)) {
          read_no_vector(*(data->ptr));
          dumpVectorRange(*(data->ptr), ugp->nno, header.name);
        } 
        else {
          if (mpi_rank == 0)
            cout << "skipping NO_DATA vector: "<< header.name<< endl;
        }
      }
      break;

    default:
      //                        if (mpi_rank == 0)
      //                                cout << "skipping CDP record: "<< header.name<< endl;
      break;

    }
    offset += header.skip;

  }

  // cleanup...
  MPI_File_close(&fh);
  MPI_Type_free(&MPI_CdpHeader);
  MPI_Type_free(&no_scalar_int_t);
  MPI_Type_free(&no_scalar_double_t);
  MPI_Type_free(&no_vector_double_t);
  MPI_Type_free(&fa_scalar_int_t);
  MPI_Type_free(&cv_scalar_int_t);

  // produce the minimum partition...
  finalize_and_check();

}

void CdpFilter::build_no_structs() {

  // just lump all nodes together for now...
  int nno_global = header.idata[0] + header.idata[1]+ header.idata[2]+ header.idata[3];
  if (mpi_rank == 0)
    cout << "nno_global: "<< nno_global << endl;

  // nodes are striped as evenly as possible across all processors...

  assert(ugp->noora == NULL);
  ugp->noora = new int[mpi_size + 1];
  calcUniformDist(ugp->noora, nno_global, mpi_size);

  assert(ugp->nno == 0);
  ugp->nno = ugp->noora[mpi_rank+1] - ugp->noora[mpi_rank];
  int my_disp = ugp->noora[mpi_rank];
  MPI_Type_indexed(1, &(ugp->nno), &my_disp, MPI_INT, &no_scalar_int_t);
  MPI_Type_commit(&no_scalar_int_t);

  MPI_Type_indexed(1, &(ugp->nno), &my_disp, MPI_DOUBLE, &no_scalar_double_t);
  MPI_Type_commit(&no_scalar_double_t);

  int nno3 = 3*(ugp->nno);
  int my_disp3 = 3*my_disp;
  MPI_Type_indexed(1, &nno3, &my_disp3, MPI_DOUBLE, &no_vector_double_t);
  MPI_Type_commit(&no_vector_double_t);

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

}

void CdpFilter::build_fa_structs() {

  // faces are striped as evenly as possible across all processors...
  int nfa_global = header.idata[0];
  if (mpi_rank == 0)
    cout << "nfa_global: "<< nfa_global << endl;

  assert(ugp->faora == NULL);
  ugp->faora = new int[mpi_size + 1];
  calcUniformDist(ugp->faora, nfa_global, mpi_size);

  assert( (ugp->nfa) == 0);
  (ugp->nfa) = ugp->faora[mpi_rank+1] - ugp->faora[mpi_rank];
  int my_disp = ugp->faora[mpi_rank];
  MPI_Type_indexed(1, &(ugp->nfa), &my_disp, MPI_INT, &fa_scalar_int_t);
  MPI_Type_commit(&fa_scalar_int_t);

  // fa_flag is used to store the face zone index...
  assert(ugp->fa_flag == NULL );
  ugp->fa_flag = new int[(ugp->nfa)];
  for (int ifa = 0; ifa < (ugp->nfa); ifa++)
    ugp->fa_flag[ifa] = -1;

}

void CdpFilter::build_cv_structs() {

  // cvs are striped as evenly as possible across all processors...
  int ncv_global = header.idata[0] + header.idata[1]+ header.idata[2]+ header.idata[3];
  if (mpi_rank == 0)
    cout << "ncv_global: "<< ncv_global << endl;

  assert(ugp->cvora == NULL);
  ugp->cvora = new int[mpi_size + 1];
  calcUniformDist(ugp->cvora, ncv_global, mpi_size);

  assert( (ugp->ncv) == 0);
  (ugp->ncv) = ugp->cvora[mpi_rank+1] - ugp->cvora[mpi_rank];
  int my_disp = ugp->cvora[mpi_rank];
  MPI_Type_indexed(1, &(ugp->ncv), &my_disp, MPI_INT, &cv_scalar_int_t);
  MPI_Type_commit(&cv_scalar_int_t);

  // ugp->cv_flag is used to store the cv zone index...
  assert(ugp->cv_flag == NULL );
  ugp->cv_flag = new int[(ugp->ncv)];
  for (int icv = 0; icv < (ugp->ncv); icv++)
    ugp->cv_flag[icv] = -1;

  // allocate any other registered cv data...
  ugp->allocateRegisteredCvData();

}

void CdpFilter::no_flag_test() {

  assert(header.idata[0] == ugp->noora[mpi_size]);
  int * my_no_flag = new int[(ugp->nno)];
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, no_scalar_int_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, my_no_flag, (ugp->nno), MPI_INT, &status);
  if (byte_swap)
    byteSwap(my_no_flag, (ugp->nno));
  for (int i = 0; i < (ugp->nno); i++) {
    if (i + ugp->noora[mpi_rank]+ 1!= my_no_flag[i]) {
      cerr << "Error: failed no_flag_test"<< endl;
      throw(-1);
    }
  }
  delete[] my_no_flag;

  //MPI_Barrier(mpi_comm);
  //if (mpi_rank == 0) 
  //  cout << "no_flag_test: OK" << endl;

}

void CdpFilter::fa_flag_test() {

  assert(header.idata[0] == ugp->faora[mpi_size]);
  int * my_fa_flag = new int[(ugp->nfa)];
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, fa_scalar_int_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, my_fa_flag, (ugp->nfa), MPI_INT, &status);
  if (byte_swap)
    byteSwap(my_fa_flag, (ugp->nfa));
  for (int i = 0; i < (ugp->nfa); i++) {
    if (i + ugp->faora[mpi_rank]+ 1!= my_fa_flag[i]) {
      cerr << "Error: failed fa_flag_test"<< endl;
      throw(-1);
    }
  }
  delete[] my_fa_flag;

  //MPI_Barrier(mpi_comm);
  //if (mpi_rank == 0)
  //  cout << "fa_flag_test: OK" << endl;
}

void CdpFilter::cv_flag_test() {

  assert(header.idata[0] == ugp->cvora[mpi_size]);
  int * my_cv_flag = new int[(ugp->ncv)];
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, cv_scalar_int_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, my_cv_flag, (ugp->ncv), MPI_INT, &status);
  if (byte_swap)
    byteSwap(my_cv_flag, (ugp->ncv));
  for (int i = 0; i < (ugp->ncv); i++) {
    if (i + ugp->cvora[mpi_rank]+ 1!= my_cv_flag[i]) {
      cerr << "Error: failed cv_flag_test"<< endl;
      throw(-1);
    }
  }
  delete[] my_cv_flag;

  //MPI_Barrier(mpi_comm);
  //if (mpi_rank == 0)
  //  cout << "cv_flag_test: OK" << endl;

}

void CdpFilter::read_cv_part() {

  assert(header.idata[0] == ugp->cvora[mpi_size]);
  assert(ugp->cv_part == NULL);
  ugp->cv_part = new int[(ugp->ncv)];
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, cv_scalar_int_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, ugp->cv_part, (ugp->ncv), MPI_INT, &status);
  if (byte_swap)
    byteSwap(ugp->cv_part, (ugp->ncv));
  int my_max_part = -1;
  for (int icv = 0; icv < (ugp->ncv); icv++) {
    assert(ugp->cv_part[icv] >= 0);
    my_max_part = max(my_max_part, ugp->cv_part[icv]);
  }
  my_max_part++;
  MPI_Allreduce(&my_max_part, &(ugp->npart), 1, MPI_INT, MPI_MAX, mpi_comm);
  if (mpi_rank == 0)
    cout << " > file contains partition info for npart: " << ugp->npart << endl;

}

void CdpFilter::add_face_zone() {

  int index = -1;
  switch (header.idata[0]) {
  case CDP_FA_ZONE_I:
  case CDP_FA_ZONE_BI:
    index = ugp->addFaZone(header.name, FA_ZONE_INTERNAL);
    break;
  case CDP_FA_ZONE_BP:
    switch (header.idata[1]) {
    case CDP_PERIODIC_CART:
      // first 3 reals of header.rdata are CDP's periodic transform...
      index = ugp->addFaZone(header.name, FA_ZONE_PERIODIC_CART, header.rdata);
      break;
    case CDP_PERIODIC_CYL_X:
      index = ugp->addFaZone(header.name, FA_ZONE_PERIODIC_CYL_X, header.rdata);
      break;
    case CDP_PERIODIC_CYL_Y:
      index = ugp->addFaZone(header.name, FA_ZONE_PERIODIC_CYL_Y, header.rdata);
      break;
    case CDP_PERIODIC_CYL_Z:
      index = ugp->addFaZone(header.name, FA_ZONE_PERIODIC_CYL_Z, header.rdata);
      break;
    default:
      cerr << "Error: unrecognized CDP_PERIODIC_TYPE: "<< header.idata[1]<< endl;
      throw(-1);
    }
    break;
  case CDP_FA_ZONE_BA:
    index = ugp->addFaZone(header.name, FA_ZONE_BOUNDARY);
    break;
  default:
    cerr << "Error: unrecognized CDP_FA_ZONE: "<< header.idata[0]<< endl;
    throw(-1);
  }

  // we now should have a unique index...
  assert(index >= 0);

  // we should always have some faces...
  assert(header.idata[4] > 0);

  // CDP restart stores the face count only with the zone contiguous, so set the
  // fa_flag of the current faces to the face zone index...
  for (int i = 0; i < header.idata[4]; i++) {
    // if this ifa_global is part of our local faces, then set the fa_zone, which
    // should have its current zone set to -1... 
    int ifa_local = ifa_global - ugp->faora[mpi_rank]; // faora[mpi_rank] contains my disp
    if ( (ifa_local >= 0)&&(ifa_local < (ugp->nfa))) {
      assert(ugp->fa_flag[ifa_local] == -1);
      ugp->fa_flag[ifa_local] = index;
    }
    ifa_global++;
  }

}

void CdpFilter::read_noofa_iv() {

  assert(header.idata[0] == ugp->faora[mpi_size]);
  assert(ugp->noofa_i == NULL);
  ugp->noofa_i = new int[(ugp->nfa)+1];
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, fa_scalar_int_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, ugp->noofa_i+1, (ugp->nfa), MPI_INT, &status);
  if (byte_swap)
    byteSwap(ugp->noofa_i+1, (ugp->nfa));
  ugp->noofa_i[0] = 0;
  for (int ifa = 0; ifa < (ugp->nfa); ifa++)
    ugp->noofa_i[ifa+1] += ugp->noofa_i[ifa];
  ugp->noofa_s = ugp->noofa_i[(ugp->nfa)];
  assert(ugp->noofa_v == NULL);
  ugp->noofa_v = new int[ugp->noofa_s];
  int my_disp;
  MPI_Scan(&(ugp->noofa_s), &my_disp, 1, MPI_INT, MPI_SUM, mpi_comm);
  // the last processor should have the global node-of-face count...
  assert( (mpi_rank < mpi_size - 1) || (my_disp == header.idata[1]));
  // change to a 0-indexed displacement
  my_disp -= ugp->noofa_s;
  // build an indexed type and read ugp->noofa_v...
  MPI_Datatype noofa_v_t;
  MPI_Type_indexed(1, &(ugp->noofa_s), &my_disp, MPI_INT, &noofa_v_t);
  MPI_Type_commit(&noofa_v_t);
  MPI_File_set_view(fh, offset+header_extent+4*ugp->faora[mpi_size], MPI_INT, noofa_v_t,
      "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, ugp->noofa_v, ugp->noofa_s, MPI_INT, &status);
  if (byte_swap)
    byteSwap(ugp->noofa_v, ugp->noofa_s);
  // convert to zero-based global indexing.
  for (int nof = 0; nof < ugp->noofa_s; nof++)
    ugp->noofa_v[nof] -= 1;
  MPI_Type_free(&noofa_v_t);

}

void CdpFilter::read_cvofa() {

  assert(ugp->cvofa == NULL);
  ugp->cvofa = new int[(ugp->nfa)][2];

  int my_size = 2*(ugp->nfa);
  int my_disp = 2*ugp->faora[mpi_rank];
  MPI_Datatype cvofa_t;
  MPI_Type_indexed(1, &my_size, &my_disp, MPI_INT, &cvofa_t);
  MPI_Type_commit(&cvofa_t);
  MPI_File_set_view(fh, offset+header_extent, MPI_INT, cvofa_t, "native", MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, ugp->cvofa, 2*(ugp->nfa), MPI_INT, &status);
  MPI_Type_free(&cvofa_t);

  if (byte_swap)
    byteSwap(ugp->cvofa, (ugp->nfa), 2);

  // convert to zero-based global indexing, check...
  for (int ifa = 0; ifa < (ugp->nfa); ifa++) {
    ugp->cvofa[ifa][0] -= 1;
    ugp->cvofa[ifa][1] -= 1;
  }

}

void CdpFilter::read_no_scalar(double *s) {

  MPI_File_set_view(fh, offset+header_extent, MPI_DOUBLE, no_scalar_double_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, s, ugp->nno, MPI_DOUBLE, &status);
  if (byte_swap)
    byteSwap(s, ugp->nno);

}

void CdpFilter::read_no_vector(double (*v)[3]) {

  MPI_File_set_view(fh, offset+header_extent, MPI_DOUBLE, no_vector_double_t, "native",
      MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, v, 3*(ugp->nno), MPI_DOUBLE, &status);
  if (byte_swap)
    byteSwap(v, (ugp->nno), 3);

}

void CdpFilter::finalize_and_check() {

  // check that all faces have their face zone set in fa_flag...
  for (int ifa = 0; ifa < ugp->nfa; ifa++) {
    assert(ugp->fa_flag[ifa] >= 0);
  }

  // check that the nodes of all faces are valid and reasonable...
  for (int ifa = 0; ifa < ugp->nfa; ifa++) {
    int nnof = 0;
    for (int nof = ugp->noofa_i[ifa]; nof < ugp->noofa_i[ifa+1]; nof++) {
      int ino = ugp->noofa_v[nof];
      // increment the number of nodes...
      nnof++;
      // node numbering should be global zero-indexed...
      assert( (ino >= 0)&&(ino < ugp->noora[mpi_size]));
      // ensure there are no nodal duplications...
      for (int nof2 = ugp->noofa_i[ifa]; nof2 < nof; nof2++) {
        assert(ugp->noofa_v[nof2] != ino );
      }
    }
    // each face should have between 3 and 8 nodes...
    assert( (nnof >= 3)&&(nnof <= 8));
  }

  {

    // it would be much better to use a vector of faZone's and avoid this
    // table, but for now...
    ugp->buildFaKindTable();

    // check that the cvs are reasonable and valid...
    for (int ifa = 0; ifa < (ugp->nfa); ifa++) {
      // first cv should always be valid...
      assert( (ugp->cvofa[ifa][0] >= 0)&&(ugp->cvofa[ifa][0]< ugp->cvora[mpi_size]));
      // second cv can be -1 for boundary faces or very negative for periodic faces
      // recall that this was a trick used for periodic reconnection in CDP restarts...
      assert(ugp->cvofa[ifa][1] < ugp->cvora[mpi_size]);
      // recall that the fa_flag contains the faZone index, so we can get the face zone
      // kind from the table we built above...
      if (ugp->cvofa[ifa][1] >= 0) {
        // if the second cv is valid, then this must be an internal face...
        assert(ugp->fa_kind_table[ugp->fa_flag[ifa]]== FA_ZONE_INTERNAL);
      } 
      else if (ugp->cvofa[ifa][1] == -1) {
        assert(ugp->fa_kind_table[ugp->fa_flag[ifa]]== FA_ZONE_BOUNDARY);
      } 
      else {
        // must be one of the periodic zones...
        assert( (ugp->fa_kind_table[ugp->fa_flag[ifa]]>= FA_ZONE_PERIODIC_FIRST)
            &&(ugp->fa_kind_table[ugp->fa_flag[ifa]]<= FA_ZONE_PERIODIC_LAST));
        // modify the cvofa to contain the matching global face index...
        ugp->cvofa[ifa][1] = -ugp->cvofa[ifa][1] - 2;
        assert( (ugp->cvofa[ifa][1] >= 0)&&(ugp->cvofa[ifa][1]< ugp->faora[mpi_size]));
        // for one rank, we can see if this makes sense...
        //assert( (ugp->nfa) == ugp->faora[mpi_size] );
        //cout << "ifa, ifa2: " << ugp->fa_flag[ifa] << " " << ugp->fa_flag[ugp->cvofa[ifa][1]] << endl; 
      }
    }

  }

  // finally, CDP does not support cell zones, so add one cell zone
  // and put all cvs in it...

  {

    assert(ugp->cvZoneList.size() == 0);
    ugp->cvZoneList.push_back(CvZone());
    list<CvZone>::iterator iz = ugp->cvZoneList.begin();
    iz->setName("fluid");
    iz->setIndex(0);
    iz->setKind(CV_ZONE_FLUID);

    // for a minimum initialization of Ugp, the cv_flag contains the
    // index of the associated cvZone - always 0 in this case...
    for (int icv = 0; icv < (ugp->ncv); icv++) {
      assert(ugp->cv_flag[icv] == -1);
      ugp->cv_flag[icv] = 0;
    }

  }

  // at this point, the minimum initialization of the Ugp is complete.

}
