#include "Ugp.h"
#include <math.h>
#include <algorithm> // required for std::sort on uP
#ifdef WITH_PARMETIS
#include "parmetis.h"
#endif

void Ugp::preinit()
{

  // call this to build relevant connections between
  // either newly periodic faces or newly internal faces...

  if (mpi_rank==0) cout<<" > Ugp::preinit()"<<endl;

  // do not assume noora, faora, cvora exist...
  buildXora(noora, nno);
  buildXora(faora, nfa);
  buildXora(cvora, ncv);

  // split any non-planar faces first...
  //splitNonPlanarFaces(1.0E-8);

  // clear any fa_kind_table...
  clearFaKindTable();

  // merge zones with common names...
  list<FaZone>::iterator zone2 = faZoneList.begin();
  while (zone2!=faZoneList.end())
  {
    // cycle zone1 up to but not including zone2
    list<FaZone>::iterator zone1;
    for (zone1 = faZoneList.begin(); zone1!=zone2; zone1++)
    {
      if (zone2->getNameString()==zone1->getNameString())
      {
        if (mpi_rank==0) cout<<" > merging zones with common name: "<<zone1->getNameString()<<endl;
        // should be the same kind...
        assert(zone2->getKind()==zone1->getKind());
        // indices should be different...
        int index2 = zone2->getIndex();
        int index1 = zone1->getIndex();
        assert(index1!=index2);
        // convert any faces associated with index2 to index1...
        for (int ifa = 0; ifa<nfa; ifa++)
          if (fa_flag[ifa]==index2) fa_flag[ifa] = index1;
        // and remove zone2...
        list<FaZone>::iterator zone2_copy = zone2;
        zone2++;
        faZoneList.erase(zone2_copy);
        // and break out of zone1 loop...
        break;
      }
    }
    // if we completed the above loop, then inc zone2...
    if (zone1==zone2) zone2++;
  }

  // ============================================================
  // do any reconnection now...
  // ============================================================

  int * fa_tmp = new int[nfa];
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_tmp[ifa] = -1;

  int nReconnect = 0;

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
  {
    int kind = zone->getKind();
    int index = zone->getIndex();
    if (kind==FA_ZONE_INTERNAL)
    {
      // all faces of this kind should have a valid icv1...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]==index)
        {
          if (cvofa[ifa][1]<0)
          {
            cerr<<"Error: cvofa[ifa][1] < 0"<<endl;
            throw(-1);
          }
        }
      }
    }
    else if (kind==FA_ZONE_BOUNDARY)
    {
      // all faces of this kind should have a -1 in there icv1...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]==index)
        {
          if (cvofa[ifa][1]!=-1)
          {
            cvofa[ifa][1] = -1;
          }
        }
      }
    }
    else if ((kind>=FA_ZONE_PERIODIC_FIRST)&&(kind<=FA_ZONE_PERIODIC_LAST))
    {
      // all faces of this kind should have a global face 
      // index indicating their associated global face...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]==index)
        {
          if (cvofa[ifa][1]==-1)
          {
            // this is a face in a periodic zone, but it does not
            // have the face index in its cvofa[ifa][1]...
            fa_tmp[ifa] = nReconnect++;
          }
          else
          {
            // the cvofa[ifa][1] should already contain the periodic
            // face...
            assert(cvofa[ifa][1]>=0);
          }
        }
      }
    }
  }

  int nReconnectGlobal;
  MPI_Allreduce(&nReconnect, &nReconnectGlobal, 1, MPI_INT, MPI_SUM, mpi_comm);
  if (nReconnectGlobal==0) return;

  if (mpi_rank==0) cout<<" > reconnecting "<<nReconnectGlobal<<" faces..."<<endl;

  // build a reverse-lookup to get the global face number from the reconnect index...
  int * fa_local_rev = new int[nReconnect];
  for (int ifa = 0; ifa<nfa; ifa++)
    if (fa_tmp[ifa]>=0) fa_local_rev[fa_tmp[ifa]] = ifa;

  // we also need these...
  double (*xReconnect)[3] = new double[nReconnect][3];
  double (*bbmin)[3] = new double[nReconnect][3];
  double (*bbmax)[3] = new double[nReconnect][3];

  assert(noora!=NULL);

  // =====================================================================
  // we need to send out a request for the nodal coordinates associated with
  // the reconnect faces so we can build the face centers...
  // =====================================================================
  int * send_count = new int[mpi_size];
  int * send_disp = new int[mpi_size];
  int * send_buffer = NULL;
  int send_count_sum;
  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_tmp[ifa]>=0)
      {
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++)
        {
          const int ino = noofa_v[nof];
          assert((ino>=0)&&(ino<noora[mpi_size]));
          // where is this node located?...
          int rank;
          for (rank = 0; rank<mpi_size; rank++)
          {
            if (ino<noora[rank+1]) break;
          }
          if (iter==0)
          {
            send_count[rank] += 1;
          }
          else
          {
            send_buffer[send_disp[rank]] = ino;
            send_disp[rank] += 1;
          }
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff...
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  int * recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  delete[] send_buffer;

  // on the recv side, cycle through the requested nodes - check
  // that we own them, then pack their coordinates...

  double * recv_buffer_double = new double[3*recv_count_sum];
  for (int i = 0; i<recv_count_sum; i++)
  {
    int ino = recv_buffer[i];
    int ino_local = ino-noora[mpi_rank];
    assert((ino_local>=0)&&(ino_local<nno));
    recv_buffer_double[3*i] = x_no[ino_local][0];
    recv_buffer_double[3*i+1] = x_no[ino_local][1];
    recv_buffer_double[3*i+2] = x_no[ino_local][2];
  }

  delete[] recv_buffer;

  double * send_buffer_double = new double[3*send_count_sum];
  for (int i = 0; i<mpi_size; i++)
  {
    recv_count[i] *= 3;
    recv_disp[i] *= 3;
    send_count[i] *= 3;
    send_disp[i] *= 3;
  }

  MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
      MPI_DOUBLE, mpi_comm);

  delete[] recv_buffer_double;

  double x_no_local[NNOF_MAX][3];

  for (int ifa = 0; ifa<nfa; ifa++)
  {

    if (fa_tmp[ifa]>=0)
    {

      int iReconnect = fa_tmp[ifa];

      // use the noofa_i/v to unpack the relevant nodes, which are
      // coming from all over...
      FOR_I3
        xReconnect[iReconnect][i] = 0.0;
      const int nof_f = noofa_i[ifa];
      const int nof_l = noofa_i[ifa+1]-1;
      assert(nof_l-nof_f+1<=NNOF_MAX);
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        const int ino = noofa_v[nof];
        assert((ino>=0)&&(ino<noora[mpi_size]));
        // where is this node located?...
        int rank;
        for (rank = 0; rank<mpi_size; rank++)
        {
          if (ino<noora[rank+1]) break;
        }
        FOR_I3
        {
          xReconnect[iReconnect][i] += send_buffer_double[send_disp[rank]];
          x_no_local[nof-nof_f][i] = send_buffer_double[send_disp[rank]];
          send_disp[rank] += 1;
        }
      }
      // normalize...
      FOR_I3
        xReconnect[iReconnect][i] /= (double) (nof_l-nof_f+1);

      // now compute the minimum edge length for this face - this
      // will get used as an appropriate local length scale for
      // populating the adt...
      double d2_min = 0.0; // gets set below
      int nof1 = nof_l-nof_f;
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        int nof0 = nof1;
        nof1 = nof-nof_f;
        double d2 = 0.0;
        for (int i = 0; i<3; i++)
        {
          double dx = x_no_local[nof1][i]-x_no_local[nof0][i];
          d2 += dx*dx;
        }
        if ((nof==nof_f)||(d2<d2_min)) d2_min = d2;
      }
      // build the bbmin,bbmax as the sphere around the face center with
      // radius a fraction of the minimum edge length...
      for (int i = 0; i<3; i++)
      {
        bbmin[iReconnect][i] = xReconnect[iReconnect][i]-0.25*sqrt(d2_min);
        bbmax[iReconnect][i] = xReconnect[iReconnect][i]+0.25*sqrt(d2_min);
      }
    }
  }

  delete[] send_buffer_double;

  /*
   // take a look...
   cout << "[" << mpi_rank << "]" << nReconnect << endl;
   for (int i = 0; i < nReconnect; ++i)
   cout << "[" << mpi_rank << "]" << i << " X: " << xReconnect[i][0] << " " << xReconnect[i][1] << " " << xReconnect[i][2] << endl;
   */

  // ========================================================
  // nodes have been exchanged and we have local xReconnect, 
  // bbmin and bbmax for each face.
  // ========================================================

  // build the local Adt...
  Adt<double> * adt = new Adt<double> (nReconnect, bbmin, bbmax);

  // everyone gets ALL adt bounding boxes...
  // these are used in making decisions on who to send our points to...
  double local_bbox[6] = { 1.0E+20, -1.0E+20, // xmin, xmax, 
      1.0E+20, -1.0E+20, // ymin, ymax...
      1.0E+20, -1.0E+20 };
  for (int ir = 0; ir<nReconnect; ++ir)
  {
    FOR_I3
    {
      local_bbox[2*i] = min(local_bbox[2*i], bbmin[ir][i]);
      local_bbox[2*i+1] = max(local_bbox[2*i+1], bbmax[ir][i]);
    }
  }
  double (*global_bbox)[6] = new double[mpi_size][6];
  MPI_Allgather(local_bbox, 6, MPI_DOUBLE, &(global_bbox[0][0]), 6, MPI_DOUBLE, mpi_comm);

  // don't need these any more...
  delete[] bbmin;
  delete[] bbmax;

  // we need the transformed coordinates to do the reconnection. Transformations
  // are virtual functions accessed using the "kind" of transformation and 
  // the zone's periodic data...
  double (*xReconnectT)[3] = new double[nReconnect][3];
  for (int ir = 0; ir<nReconnect; ir++)
    for (int i = 0; i<3; i++)
      xReconnectT[ir][i] = xReconnect[ir][i];

  int iReconnect = 0;
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
  {
    int kind = zone->getKind();
    int index = zone->getIndex();
    if (kind==FA_ZONE_INTERNAL)
    {
      // all faces of this kind should have a valid icv1...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]==index)
        {
          if (cvofa[ifa][1]<0)
          {
            cerr<<"Error: cvofa[ifa][1] < 0"<<endl;
            throw(-1);
          }
        }
      }
    }
    else if ((kind>=FA_ZONE_PERIODIC_FIRST)&&(kind<=FA_ZONE_PERIODIC_LAST))
    {
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]==index)
        {
          if (cvofa[ifa][1]==-1)
          {
            assert(fa_tmp[ifa]==iReconnect);
            rangeTranslateVector(xReconnectT[iReconnect], 1, kind, zone->periodic_data);
            iReconnect++;
          }
        }
      }
    }
  }
  assert(iReconnect==nReconnect);

  delete[] fa_tmp;

  // send out the transformed points to all processors with possible matches...

  // a way to get the reconnect index from the send index...
  int * ir_of_isend = NULL;

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;

  for (int iter = 0; iter<2; iter++)
  {
    for (int rank = 0; rank<mpi_size; ++rank)
    {
      for (int ir = 0; ir<nReconnect; ir++)
      {
        if ((xReconnectT[ir][0]>=global_bbox[rank][0])&&(xReconnectT[ir][0]<=global_bbox[rank][1])&&(xReconnectT[ir][1]
            >=global_bbox[rank][2])&&(xReconnectT[ir][1]<=global_bbox[rank][3])&&(xReconnectT[ir][2]
            >=global_bbox[rank][4])&&(xReconnectT[ir][2]<=global_bbox[rank][5]))
        {
          // this is point that should be sent to this rank...
          if (iter==0)
          {
            send_count[rank] += 3;
          }
          else
          {
            FOR_I3
              send_buffer_double[send_disp[rank]+i] = xReconnectT[ir][i];
            ir_of_isend[send_disp[rank]/3] = ir;
            send_disp[rank] += 3;
          }
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer_double = new double[send_count_sum];
      ir_of_isend = new int[send_count_sum/3];
    }
  }

  delete[] xReconnectT;
  delete[] global_bbox;

  // now build the recv side stuff...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];

  recv_buffer_double = new double[recv_count_sum];
  MPI_Alltoallv(send_buffer_double, send_count, send_disp, MPI_DOUBLE, recv_buffer_double, recv_count, recv_disp,
      MPI_DOUBLE, mpi_comm);

  delete[] send_buffer_double;

  // now use the adt on the recv side to find the closest local match, if
  // one exists...

  recv_count_sum /= 3;

  int * recv_ifa_pair = new int[recv_count_sum];
  double * recv_ifa_d2 = new double[recv_count_sum];
  int bbListSize, bbList[ADT_LIST_MAX];

  for (int irecv = 0; irecv<recv_count_sum; ++irecv)
  {
    recv_ifa_pair[irecv] = -1;
    recv_ifa_d2[irecv] = 0.0;
    double point[3];
    point[0] = recv_buffer_double[3*irecv];
    point[1] = recv_buffer_double[3*irecv+1];
    point[2] = recv_buffer_double[3*irecv+2];
    // query the Adt for a match...
    adt->buildListForPoint(bbListSize, bbList, point);
    // XXX when reconnecting internal faces, we will need to 
    // skip a perfect match - LATER...
    for (int ibb = 0; ibb<bbListSize; ++ibb)
    {
      double d2 = 0.0;
      int ir = bbList[ibb]; // reconnect index 
      FOR_I3
      {
        double dx = xReconnect[ir][i]-point[i];
        d2 += dx*dx;
      }
      if ((recv_ifa_pair[irecv]==-1)||(d2<recv_ifa_d2[irecv]))
      {
        recv_ifa_pair[irecv] = fa_local_rev[ir]+faora[mpi_rank]; // add the global offset
        recv_ifa_d2[irecv] = d2;
      }
    }
  }

  // cleanup...

  delete adt;
  delete[] xReconnect;
  delete[] recv_buffer_double;

  // now send this result back to the send process...

  for (int i = 0; i<mpi_size; i++)
  {
    recv_count[i] /= 3;
    recv_disp[i] /= 3;
    send_count[i] /= 3;
    send_disp[i] /= 3;
  }

  send_count_sum /= 3;

  int * send_ifa_pair = new int[send_count_sum];
  double * send_ifa_d2 = new double[send_count_sum];
  MPI_Alltoallv(recv_ifa_pair, recv_count, recv_disp, MPI_INT, send_ifa_pair, send_count, send_disp, MPI_INT, mpi_comm);
  MPI_Alltoallv(recv_ifa_d2, recv_count, recv_disp, MPI_DOUBLE, send_ifa_d2, send_count, send_disp, MPI_DOUBLE,
      mpi_comm);

  delete[] recv_ifa_pair;
  delete[] recv_ifa_d2;

  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  // use this fa_pair to store the face pair with the closest d2...

  int * fa_pair = new int[nReconnect];
  double * fa_d2 = new double[nReconnect];
  for (int ir = 0; ir<nReconnect; ++ir)
    fa_pair[ir] = -1;

  for (int isend = 0; isend<send_count_sum; ++isend)
  {
    // we can get the reconnect index from the sedn index...
    int ir = ir_of_isend[isend];
    assert((ir>=0)&&(ir<nReconnect));
    // the ifa_pair and ifa_d2 contain the global face match and 
    // its d2 (distance squared) respectively...
    if ((fa_pair[ir]==-1)||(fa_d2[ir]<send_ifa_d2[isend]))
    {
      fa_pair[ir] = send_ifa_pair[isend];
      fa_d2[ir] = send_ifa_d2[isend];
    }
  }

  delete[] send_ifa_pair;
  delete[] send_ifa_d2;
  delete[] ir_of_isend;

  // now we have all the matched faces back. 

  double my_d2_max = 0.0;
  for (int ir = 0; ir<nReconnect; ir++)
  {
    int ifa_local = fa_local_rev[ir];
    assert(fa_pair[ir]>=0); // make sure a match was found
    assert(cvofa[ifa_local][1]==-1);
    cvofa[ifa_local][1] = fa_pair[ir]; // global face index
    my_d2_max = max(fa_d2[ir], my_d2_max);
  }

  delete[] fa_pair;
  delete[] fa_d2;
  delete[] fa_local_rev;

  double d2_max;
  MPI_Reduce(&my_d2_max, &d2_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank==0) cout<<" > reconnected faces matched within tol: "<<sqrt(d2_max)<<endl;

}

void Ugp::splitNonPlanarFaces(const double tol)
{

  if (mpi_rank==0) cout<<"splitNonPlanarFaces()"<<endl;

  // called BEFORE grid is redistributed.

  assert(faora[mpi_rank+1]-faora[mpi_rank]==nfa);
  assert(noora[mpi_rank+1]-noora[mpi_rank]==nno);

  // We need the face nodes to determine whether face is non-coplanar enough
  // to split...

  // Buffers..
  int * send_count = new int[mpi_size];
  int * send_disp = new int[mpi_size];
  int * send_buffer_int = NULL;
  int send_count_sum;
  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    FOR_IFA
    {
      const int nof_f = noofa_i[ifa];
      const int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof<=nof_l; ++nof)
      {
        // at this point, noofa_v stores the GLOBAL node index...
        const int ino = noofa_v[nof];
        assert((ino>=0)&&(ino<noora[mpi_size]));
        int rank = 0;
        while (noora[rank+1]<=ino)
          ++rank;
        if (iter==0) send_count[rank] += 1;
        else
        {
          send_buffer_int[send_disp[rank]] = ino-noora[rank]; // local node index
          send_disp[rank] += 1;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer_int = new int[send_count_sum];
    }
  }

  // now build the recv side stuff...
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  int * recv_buffer_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer_int, send_count, send_disp, MPI_INT, recv_buffer_int, recv_count, recv_disp, MPI_INT,
      mpi_comm);

  // some cleanup...
  delete[] send_buffer_int;

  // now load up the nodal coordinates...
  double (*recv_buffer_x)[3] = new double[recv_count_sum][3];
  for (int i = 0; i<recv_count_sum; ++i)
  {
    const int ino = recv_buffer_int[i];
    assert((ino>=0)&&(ino<nno));
    recv_buffer_x[i][0] = x_no[ino][0];
    recv_buffer_x[i][1] = x_no[ino][1];
    recv_buffer_x[i][2] = x_no[ino][2];
  }

  delete[] recv_buffer_int;

  // and send back: here we have three doubles per previous int...
  for (int i = 0; i<mpi_size; i++)
  {
    send_count[i] *= 3;
    send_disp[i] *= 3;
    recv_count[i] *= 3;
    recv_disp[i] *= 3;
  }
  double (*send_buffer_x)[3] = new double[send_count_sum][3];
  MPI_Alltoallv(recv_buffer_x, recv_count, recv_disp, MPI_DOUBLE, send_buffer_x, send_count, send_disp, MPI_DOUBLE,
      mpi_comm);

  // some more cleanup...
  delete[] recv_buffer_x;

  // return send_disp... 
  for (int i = 0; i<mpi_size; i++)
    send_disp[i] /= 3;

  // x_fa...
  int * fa_split = new int[nfa];
  int nfa_new = 0;
  int noofa_s_new = 0;
  double myMaxComponent = 0.0;
  double x_noofa[NNOF_MAX][3];
  FOR_IFA
  {
    // put the action to take in fa_split...
    fa_split[ifa] = 0;
    // loop on nodes...
    const int nof_f = noofa_i[ifa];
    const int nof_l = noofa_i[ifa+1]-1;
    assert(nof_l-nof_f+1<=NNOF_MAX);
    for (int nof = nof_f; nof<=nof_l; ++nof)
    {
      // at this point, noofa_v stores the GLOBAL node index...
      const int ino = noofa_v[nof];
      assert((ino>=0)&&(ino<noora[mpi_size]));
      int rank = 0;
      while (noora[rank+1]<=ino)
        ++rank;
      FOR_I3
        x_noofa[nof-nof_f][i] = send_buffer_x[send_disp[rank]][i];
      send_disp[rank] += 1;
    }
    // =================================
    // check how planar the face is...
    // =================================
    // 1. approximate center...
    const int nnof = nof_l-nof_f+1;
    double x_fa_approx[3] = { 0.0, 0.0, 0.0 };
    for (int nof = 0; nof<nnof; ++nof)
      FOR_I3
        x_fa_approx[i] += x_noofa[nof][i];
    double tmp = 1.0/(double) nnof;
    FOR_I3
      x_fa_approx[i] *= tmp;
    // 2. face normal...
    double fa_normal[3] = { 0.0, 0.0, 0.0 };
    int nof1 = nnof-1;
    for (int nof = 0; nof<nnof; ++nof)
    {
      int nof0 = nof1;
      nof1 = nof;
      double v0[3], v1[3];
      FOR_I3
      {
        v0[i] = x_noofa[nof0][i]-x_fa_approx[i];
        v1[i] = x_noofa[nof1][i]-x_fa_approx[i];
      }
      fa_normal[0] += 0.5*(v0[1]*v1[2]-v0[2]*v1[1]);
      fa_normal[1] += 0.5*(v0[2]*v1[0]-v0[0]*v1[2]);
      fa_normal[2] += 0.5*(v0[0]*v1[1]-v0[1]*v1[0]);
    }
    // 3. compare component of vector from approx face center to
    // the node (should be zero) against an appropriate length scale...
    double maxComponent = 0.0;
    for (int nof = 0; nof<nnof; ++nof)
    {
      maxComponent = max(maxComponent, fabs((x_noofa[nof][0]-x_fa_approx[0])*fa_normal[0]+(x_noofa[nof][1]
          -x_fa_approx[1])*fa_normal[1]+(x_noofa[nof][2]-x_fa_approx[2])*fa_normal[2]));
    }
    // normlize
    double area2 = fa_normal[0]*fa_normal[0]+fa_normal[1]*fa_normal[1]+fa_normal[2]*fa_normal[2];
    maxComponent /= pow(area2, 0.75);
    // store the global max for diagnostics...
    myMaxComponent = max(myMaxComponent, maxComponent);
    // take action at non-coplanar faces...
    if (maxComponent>=tol)
    {
      // decide what action to take based on the type of face...
      if (nnof==4)
      {
        // we will need one new face...
        nfa_new += 1;
        // we need a position for 3 new nodes, minus 1 already present to
        // describe the quad.
        noofa_s_new += 2;
        // quad face - insert diaagonal and form 2 faces...
        // choose minimum distance diagonal...
        double dist02 = 0.0;
        double dist13 = 0.0;
        FOR_I3
        {
          double dx = x_noofa[2][i]-x_noofa[0][i];
          dist02 += dx*dx;
          dx = x_noofa[3][i]-x_noofa[1][i];
          dist13 += dx*dx;
        }
        if (dist02<=dist13) fa_split[ifa] = 1;
        else fa_split[ifa] = 2;
      }
      else
      {
        cout<<"unsupported non-planar face: "<<nnof<<endl;
        throw(-1);
      }
    }
  }

  // cleanup...
  delete[] send_buffer_x;
  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  int my_ibuf[2] = { nfa, nfa_new };
  int ibuf[2];
  MPI_Allreduce(my_ibuf, ibuf, 2, MPI_INT, MPI_SUM, mpi_comm);

  double globalMaxComponent;
  MPI_Reduce(&myMaxComponent, &globalMaxComponent, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank==0)
  {
    cout<<" > adding new faces due to non-planar splitting: "<<ibuf[1]<<" out of "<<ibuf[0]<<" faces total"<<endl;
    cout<<" > globalMaxComponent (max normalized measure of face non-planarity): "<<globalMaxComponent<<endl;
  }

  // we used Allreduce so we can all make this decision together...

  if (ibuf[1]>0)
  {

    // however only ranks with new faces need do the following...

    if (nfa_new>0)
    {

      assert(noofa_s==noofa_i[nfa]);

      // resize all data, copying existing values...
      resize(noofa_i, nfa+1, nfa+nfa_new+1);
      resize(noofa_v, noofa_s, noofa_s+noofa_s_new);
      resize(cvofa, nfa, nfa+nfa_new);
      resize(fa_flag, nfa, nfa+nfa_new);

      // loop backward thru original faces...
      int ifa_offset = nfa_new;
      int nof_offset = noofa_s_new;
      for (int ifa = nfa-1; ifa>=0; --ifa)
      {
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        noofa_i[ifa+ifa_offset+1] = nof_l+nof_offset+1;
        switch (fa_split[ifa])
        {
        case 0:
          // just shift face-based data using offsets...
          // use reverse order here...
          for (int nof = nof_l; nof>=nof_f; --nof)
            noofa_v[nof+nof_offset] = noofa_v[nof];
          cvofa[ifa+ifa_offset][0] = cvofa[ifa][0];
          cvofa[ifa+ifa_offset][1] = cvofa[ifa][1];
          fa_flag[ifa+ifa_offset] = fa_flag[ifa];
          break;
        case 1:
        {
          // should be quad...
          assert(nof_l-nof_f+1==4);
          // take the quad nodes (global indexing, but doesn't matter)...
          int ino0 = noofa_v[nof_f];
          int ino1 = noofa_v[nof_f+1];
          int ino2 = noofa_v[nof_f+2];
          int ino3 = noofa_v[nof_f+3];
          // face1 is tri 0-1-2...
          noofa_i[ifa+ifa_offset] = nof_f+nof_offset+1;
          noofa_v[nof_f+nof_offset+1] = ino0;
          noofa_v[nof_f+nof_offset+2] = ino1;
          noofa_v[nof_f+nof_offset+3] = ino2;
          // face0 is tri 0-2-3...
          noofa_v[nof_f+nof_offset-2] = ino0;
          noofa_v[nof_f+nof_offset-1] = ino2;
          noofa_v[nof_f+nof_offset] = ino3;
          // cvofa's are same. This violates the one-face-per-cv pair
          // assumption, but we need to be able to handle this for
          // GG. accuracy...
          cvofa[ifa+ifa_offset][0] = cvofa[ifa][0];
          cvofa[ifa+ifa_offset][1] = cvofa[ifa][1];
          cvofa[ifa+ifa_offset-1][0] = cvofa[ifa][0];
          cvofa[ifa+ifa_offset-1][1] = cvofa[ifa][1];
          // same zone...
          fa_flag[ifa+ifa_offset] = fa_flag[ifa];
          fa_flag[ifa+ifa_offset-1] = fa_flag[ifa];
          // adjust nof_offset by 2 only...
          nof_offset -= 2;
          ifa_offset -= 1;
        }
          break;
        case 2:
        {
          // should be quad...
          assert(nof_l-nof_f+1==4);
          // take the quad nodes (global indexing, but doesn't matter)...
          int ino0 = noofa_v[nof_f];
          int ino1 = noofa_v[nof_f+1];
          int ino2 = noofa_v[nof_f+2];
          int ino3 = noofa_v[nof_f+3];
          // face1 is tri 0-1-3...
          noofa_i[ifa+ifa_offset] = nof_f+nof_offset+1;
          noofa_v[nof_f+nof_offset+1] = ino0;
          noofa_v[nof_f+nof_offset+2] = ino1;
          noofa_v[nof_f+nof_offset+3] = ino3;
          // face0 is tri 1-2-3...
          noofa_v[nof_f+nof_offset-2] = ino1;
          noofa_v[nof_f+nof_offset-1] = ino2;
          noofa_v[nof_f+nof_offset] = ino3;
          // cvofa's are same. This violates the one-face-per-cv pair
          // assumption, but we need to be able to handle this for
          // Green-Gauss accuracy...
          cvofa[ifa+ifa_offset][0] = cvofa[ifa][0];
          cvofa[ifa+ifa_offset][1] = cvofa[ifa][1];
          cvofa[ifa+ifa_offset-1][0] = cvofa[ifa][0];
          cvofa[ifa+ifa_offset-1][1] = cvofa[ifa][1];
          // same zone...
          fa_flag[ifa+ifa_offset] = fa_flag[ifa];
          fa_flag[ifa+ifa_offset-1] = fa_flag[ifa];
          // adjust nof_offset by 2 only...
          nof_offset -= 2;
          ifa_offset -= 1;
        }
          break;
        default:
          cout<<"Error: unsupported face splitting: "<<fa_split[ifa]<<endl;
          throw(-1);
        }
      }

      // make sure we got zero offsets...
      assert(ifa_offset==0);
      assert(nof_offset==0);

      // now adjust nfa...
      nfa += nfa_new;
      noofa_s += noofa_s_new;
      assert(noofa_i[nfa]==noofa_s);

    }

    // we ALL adjust faora...
    MPI_Allgather(&nfa, 1, MPI_INT, faora+1, 1, MPI_INT, mpi_comm);
    faora[0] = 0;
    for (int rank = 0; rank<mpi_size; ++rank)
      faora[rank+1] += faora[rank];

    // we did this at the top...
    assert(faora[mpi_rank+1]-faora[mpi_rank]==nfa);

  }

  delete[] fa_split;

} // void Ugp::splitNonPlanarFaces() {

void Ugp::setCvGroupPartMetis()
{

  MPI_Pause("setCvGroupPartMetis");

}

void Ugp::setCvPartMetis()
{

  // recomputes this ugp's cv_part and npart to be load-balanced and 
  // consistent with the current mpi_size...

  if (mpi_rank==0) cout<<"setCvPartMetis()"<<endl;

  dumpScalarRange(fa_flag, nfa, "FA_FLAG");
  //dumpScalarBins(fa_flag,nfa,"FA_FLAG");

  if (cv_part==NULL) cv_part = new int[ncv];

  // could make this more than mpi_size, but for now...
  npart = mpi_size;

  // for the case of just 1 partition, there is no 
  // need to call ParMetis/Metis...
  if (mpi_size==1)
  {
    for (int icv = 0; icv<ncv; icv++)
      cv_part[icv] = 0;
    return;
  }

  if (fa_kind_table==NULL) buildFaKindTable();

  // now we need to build a temporary nbocv_i/v structure to 
  // compute the Metis/ParMetis partition. include only
  // neighbours connected through internal faces...

  int * send_count = new int[mpi_size];
  int * send_disp = new int[mpi_size];
  int * send_buffer = NULL;
  int send_count_sum;
  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      assert(fa_kind_table[fa_flag[ifa]]>=0);
      if (fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL)
      {
        // we will be sending a message to the processor(s) of each of
        // these cv's...
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];
        // make sure both cv's are valid global cv indices...
        assert((icv0>=0)&&(icv0<cvora[mpi_size]));
        assert((icv1>=0)&&(icv1<cvora[mpi_size]));
        // icv0...
        int rank0;
        for (rank0 = 0; rank0<mpi_size; rank0++)
        {
          if (icv0<cvora[rank0+1]) break;
        }
        // icv1...
        int rank1;
        for (rank1 = 0; rank1<mpi_size; rank1++)
        {
          if (icv1<cvora[rank1+1]) break;
        }
        // to save some buffer size, treat same-rank cv's slightly differently...
        if (rank0==rank1)
        {
          // when the processor is the same, we send the pair once
          // and use logic on the recv processor to figure this out...
          if (iter==0)
          {
            send_count[rank0] += 2;
          }
          else
          {
            send_buffer[send_disp[rank0]] = icv0;
            send_buffer[send_disp[rank0]+1] = icv1;
            send_disp[rank0] += 2;
          }
        }
        else
        {
          // send the cv pair to both processors...
          if (iter==0)
          {
            send_count[rank0] += 2;
            send_count[rank1] += 2;
          }
          else
          {
            send_buffer[send_disp[rank0]] = icv0;
            send_buffer[send_disp[rank0]+1] = icv1;
            send_disp[rank0] += 2;
            send_buffer[send_disp[rank1]] = icv1;
            send_buffer[send_disp[rank1]+1] = icv0;
            send_disp[rank1] += 2;
          }
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff...
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  int * recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // and build the nbocv_tmp_i/v struct...
  int * nbocv_tmp_i = new int[ncv+1];
  int * nbocv_tmp_v = NULL;
  for (int icv = 0; icv<ncv; icv++)
    nbocv_tmp_i[icv+1] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int i = 0; i<recv_count_sum; i += 2)
    {
      // the two global cvs are...
      int icv0 = recv_buffer[i];
      int icv1 = recv_buffer[i+1];
      // icv0 should always be local...
      int icv0_local = icv0-cvora[mpi_rank];
      assert((icv0_local>=0)&&(icv0_local<ncv));
      if (iter==0)
      {
        nbocv_tmp_i[icv0_local+1] += 1;
      }
      else
      {
        nbocv_tmp_v[nbocv_tmp_i[icv0_local]] = icv1; // immediately to global indexing as required by Metis...        
        nbocv_tmp_i[icv0_local] += 1;
      }
      // icv1 may or may not be local...
      int icv1_local = icv1-cvora[mpi_rank];
      if ((icv1_local>=0)&&(icv1_local<ncv))
      {
        if (iter==0)
        {
          nbocv_tmp_i[icv1_local+1] += 1;
        }
        else
        {
          nbocv_tmp_v[nbocv_tmp_i[icv1_local]] = icv0; // immediately to global indexing as required by Metis...        
          nbocv_tmp_i[icv1_local] += 1;
        }
      }
    }
    if (iter==0)
    {
      // convert counts to CSR structure...
      nbocv_tmp_i[0] = 0;
      for (int icv = 0; icv<ncv; icv++)
        nbocv_tmp_i[icv+1] += nbocv_tmp_i[icv];
      int nbocv_tmp_s = nbocv_tmp_i[ncv];
      nbocv_tmp_v = new int[nbocv_tmp_s];
    }
    else
    {
      // restore CSR structure...
      for (int icv = ncv; icv>0; icv--)
        nbocv_tmp_i[icv] = nbocv_tmp_i[icv-1];
      nbocv_tmp_i[0] = 0;
    }
  }

  // secondly, compute a cv-weight that is the sum of the 
  // face node counts. In this case, include all faces...

  delete[] send_buffer;
  delete[] recv_buffer;

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // all faces have a valid icv0...
      int icv0 = cvofa[ifa][0];
      assert((icv0>=0)&&(icv0<cvora[mpi_size]));
      int rank0;
      for (rank0 = 0; rank0<mpi_size; rank0++)
      {
        if (icv0<cvora[rank0+1]) break;
      }
      if (iter==0)
      {
        send_count[rank0] += 2;
      }
      else
      {
        send_buffer[send_disp[rank0]] = icv0;
        send_buffer[send_disp[rank0]+1] = noofa_i[ifa+1]-noofa_i[ifa];
        send_disp[rank0] += 2;
      }
      // only internal faces have a valid icv1...
      if (fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL)
      {
        int icv1 = cvofa[ifa][1];
        assert((icv1>=0)&&(icv1<cvora[mpi_size]));
        int rank1;
        for (rank1 = 0; rank1<mpi_size; rank1++)
        {
          if (icv1<cvora[rank1+1]) break;
        }
        if (iter==0)
        {
          send_count[rank1] += 2;
        }
        else
        {
          send_buffer[send_disp[rank1]] = icv1;
          send_buffer[send_disp[rank1]+1] = noofa_i[ifa+1]-noofa_i[ifa];
          send_disp[rank1] += 2;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // and build the cv_weights used in the partitioning...
  int * cv_weights = new int[ncv]; // for future dual constraint, double this.
  for (int icv = 0; icv<ncv; icv++)
    cv_weights[icv] = 0;
  for (int i = 0; i<recv_count_sum; i += 2)
  {
    // the global cv is first...
    int icv = recv_buffer[i];
    // the nnof is second...
    int nnof = recv_buffer[i+1];
    // icv should always be local...
    int icv_local = icv-cvora[mpi_rank];
    assert((icv_local>=0)&&(icv_local<ncv));
    cv_weights[icv_local] += nnof;
  }

  // or, go to uniform weighting...
  for (int icv = 0; icv<ncv; icv++)
    cv_weights[icv] = 1;

  dumpScalarRange(cv_weights, ncv, "CV_WEIGHTS");

  delete[] send_count;
  delete[] send_disp;
  delete[] send_buffer;

  delete[] recv_count;
  delete[] recv_disp;
  delete[] recv_buffer;

  // and partition...

  {

    int wgtflag = 2; // 0: no weights, 1: edge weights only, 2: vertex weights only, 3: both 
    int numflag = 0; // c-style indexing
    int ncon = 1; // number of constraints

    float * tpwgts = new float[ncon*npart];
    for (int i = 0; i<ncon*npart; i++)
      tpwgts[i] = 1.0/(float) npart;

    float ubvec[2]; // MAXNCON
    for (int i = 0; i<ncon; i++)
      ubvec[i] = 1.05;

    int options[3];
    options[0] = 1; // 0: following 2 options are skipped, 1: use following options 
    options[1] = 3; // print more timing info: default: 0  
    options[2] = 15; // random seed, default: 15

    int edgecut;

#ifndef NO_MPI       // preliminary like this
    //#ifdef WITH_PARMETIS
    ParMETIS_V3_PartKway(cvora, nbocv_tmp_i, nbocv_tmp_v, cv_weights, NULL, &wgtflag, &numflag, &ncon, &npart, tpwgts,
        ubvec, options, &edgecut, cv_part, &mpi_comm);
    ///#endif
#endif

    delete[] tpwgts;

    /*
     // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     // bad partition HACK...
     // use this to test partition-independence of your Ugp-based solver
     if (mpi_rank == 0)
     cout << "Warning: making very bad partition..." << endl;
     for (int icv = 0; icv < ncv; icv++)
     cv_part[icv] = icv%npart;
     // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     */

    dumpScalarRange(cv_part, ncv, "CV_PART");
    //dumpScalarBins(cv_part,ncv,"CV_PART");

  }

  delete[] nbocv_tmp_i;
  delete[] nbocv_tmp_v;
  delete[] cv_weights;

}

void Ugp::setCvPartCheap()
{

  if (mpi_rank==0) cout<<"setCvPartCheap()"<<endl;

  if (cv_part==NULL) cv_part = new int[ncv];

  for (int icv = 0; icv<ncv; icv++)
    cv_part[icv] = mpi_rank;

  npart = mpi_size;

}

void Ugp::setCvPartBad()
{

  if (mpi_rank==0) cout<<"setCvPartBad()"<<endl;

  if (cv_part==NULL) cv_part = new int[ncv];

  int this_part = (mpi_rank+1)%mpi_size;
  for (int icv = 0; icv<ncv; icv++)
  {
    cv_part[icv] = this_part;
    this_part = (this_part+1)%mpi_size;
  }

  npart = mpi_size;

}

void Ugp::checkCvPart()
{

  assert(cv_part!=NULL);

  int * my_cv_part_count = new int[npart];
  for (int i = 0; i<npart; i++)
    my_cv_part_count[i] = 0;

  for (int icv = 0; icv<ncv; icv++)
  {
    assert((cv_part[icv]>=0)&&(cv_part[icv]<npart));
    my_cv_part_count[cv_part[icv]] += 1;
  }

  int * cv_part_count = new int[npart];
  MPI_Reduce(my_cv_part_count, cv_part_count, npart, MPI_INT, MPI_SUM, 0, mpi_comm);
  if (mpi_rank==0)
  {
    int ncv_sum = cv_part_count[0];
    int ncv_max = cv_part_count[0];
    int ncv_min = cv_part_count[0];
    for (int i = 1; i<npart; i++)
    {
      ncv_sum += cv_part_count[i];
      ncv_max = max(ncv_max, cv_part_count[i]);
      ncv_min = min(ncv_min, cv_part_count[i]);
    }
    cout<<"checkCvPart(), npart: "<<npart<<", ncv global: "<<ncv_sum<<" min/avg/max per proc: "<<ncv_min<<" "
        <<(double) ncv_sum/(double) npart<<" "<<ncv_max<<endl;
    /*
     for (int i = 0; i < npart; i++)
     cout << "rank: " << i << ", cv count: " << cv_part_count[i] << endl;
     */
  }
  delete[] cv_part_count;

  delete[] my_cv_part_count;

}

void Ugp::redistReconnectReorder()
{

  // performs the tedious task of redistributing the grid given a new 
  // cv partition. This routine assumes the following minimum conditions
  // exist:
  //
  // 1. faces are uniquely distributed across the processors according 
  //    to faora[0..mpi_size], and also have local cvofa and noofa_i/v
  //    structures set in terms of global (and not neccessarily on the
  //    same processor cells and nodes respectively). fa_flag also contains
  //    a zero or positive index of the associated face zone, from which the 
  //    face kind (INTERNAL, BOUNDARY. PERIODIC, ) can be determined.
  //
  // 2. cells are uniquely distributed across the processors according to 
  //    cvora. the cv_part contains the 0-indexed partition info.
  // 
  // 3. nodes are uniquely distributed across the processors according to
  //    noora. Mainly data and specifically x_no is associated with them.
  // 
  // for ALL faces in ALL zones, icv0 (cvofa[ifa][0]) is always a valid
  // global cv index. For boundary faces, icv1 should be -1. For periodic faces
  // icv1 should contain the matching FACE index as a -2 indexed number ( i.e. 
  // 0 == -2, 1 == -3, etc).
  //
  // also, all cell and face zones should be synchronized across all processors. 

  if (mpi_rank==0) cout<<"redistReconnectReorder()"<<endl;

  // check that the partition is as expected...
  assert(npart==mpi_size);
  assert(cv_part!=NULL);
  for (int icv = 0; icv<ncv; icv++)
  {
    assert((cv_part[icv]>=0)&&(cv_part[icv]<npart));
  }

  // we also need the fa_kind_table...
  if (fa_kind_table==NULL) buildFaKindTable();

  // looks good, so now do redistribution of cv's and all associated faces/nodes...
  // note that we completely skip edges for now. They will have to be
  // added back in prior to developing a consistent unstructured refinement
  // capability, but they are not required in general...

  // start by having the face-partition request the cv-partition for 
  // each of its cv's...

  int * send_count = new int[mpi_size];
  int * send_disp = new int[mpi_size];
  int * send_buffer;
  int send_count_sum;
  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // all faces have a valid icv0...
      int icv0 = cvofa[ifa][0];
      assert((icv0>=0)&&(icv0<cvora[mpi_size]));
      int rank0;
      for (rank0 = 0; rank0<mpi_size; rank0++)
      {
        if (icv0<cvora[rank0+1]) break;
      }
      if (iter==0)
      {
        send_count[rank0] += 1;
      }
      else
      {
        send_buffer[send_disp[rank0]] = icv0;
        send_disp[rank0] += 1;
      }
      // only internal faces have a valid icv1...
      if (fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL)
      {
        int icv1 = cvofa[ifa][1];
        assert((icv1>=0)&&(icv1<cvora[mpi_size]));
        int rank1;
        for (rank1 = 0; rank1<mpi_size; rank1++)
        {
          if (icv1<cvora[rank1+1]) break;
        }
        if (iter==0)
        {
          send_count[rank1] += 1;
        }
        else
        {
          send_buffer[send_disp[rank1]] = icv1;
          send_disp[rank1] += 1;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  int recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  int * recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // replace the contents of the recv buffer with the rank from cv_part...
  for (int i = 0; i<recv_count_sum; i++)
  {
    int icv = recv_buffer[i];
    int icv_local = icv-cvora[mpi_rank];
    assert((icv_local>=0)&&(icv_local<ncv));
    recv_buffer[i] = cv_part[icv_local];
  }

  MPI_Alltoallv(recv_buffer, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);

  delete[] recv_buffer;

  // now build the rank-of-cv-of-face... 
  int (*raocvofa)[2] = new int[nfa][2];
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    // all faces have a valid icv0...
    int icv0 = cvofa[ifa][0];
    assert((icv0>=0)&&(icv0<cvora[mpi_size]));
    int rank0;
    for (rank0 = 0; rank0<mpi_size; rank0++)
    {
      if (icv0<cvora[rank0+1]) break;
    }
    raocvofa[ifa][0] = send_buffer[send_disp[rank0]];
    send_disp[rank0] += 1;
    // only internal faces have a valid icv1...
    if (fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL)
    {
      int icv1 = cvofa[ifa][1];
      assert((icv1>=0)&&(icv1<cvora[mpi_size]));
      int rank1;
      for (rank1 = 0; rank1<mpi_size; rank1++)
      {
        if (icv1<cvora[rank1+1]) break;
      }
      raocvofa[ifa][1] = send_buffer[send_disp[rank1]];
      send_disp[rank1] += 1;
    }
    else
    {
      // it is critical that if the second cv is not valid, set this
      // quantity to -1...
      raocvofa[ifa][1] = -1;
    }
  }

  delete[] send_buffer;

  // -------------------------------------------------
  // now send the faces to the appropriate ranks...
  // to do this, we need to count and pack the face flags (zones), 
  // their global cvs, and global nodes (ordered properly)...
  // -------------------------------------------------

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // how many nodes...
      int nnof = noofa_i[ifa+1]-noofa_i[ifa];
      assert(nnof>=3);
      // all faces have a valid icv0...
      int icv0 = cvofa[ifa][0];
      int rank0 = raocvofa[ifa][0]; // the rank where this cv will finally be
      assert(rank0>=0);
      // not all faces have a valid icv1...
      int icv1 = cvofa[ifa][1];
      int rank1 = raocvofa[ifa][1]; // the rank where this cv will finally be
      if (rank1==-1)
      {
        // this is a boundary or periodic face, so icv1 contains either -1
        // or the global face index associated with the periodicity...
        // check this...
        if (icv1==-1)
        {
          assert(fa_kind_table[fa_flag[ifa]]==FA_ZONE_BOUNDARY);
        }
        else
        {
          assert((fa_kind_table[fa_flag[ifa]]>=FA_ZONE_PERIODIC_FIRST)&&(fa_kind_table[fa_flag[ifa]]
              <=FA_ZONE_PERIODIC_LAST));
          assert((icv1>=0)&&(icv1<faora[mpi_size]));
        }
        // in either case, send the icv1 as the second data...
        if (iter==0)
        {
          send_count[rank0] += 5+nnof; // icv0, icv1 [-1 or periodic face], fa_zone, old ifa, nnof, nodes.
        }
        else
        {
          send_buffer[send_disp[rank0]] = icv0;
          send_buffer[send_disp[rank0]+1] = icv1; // actually a global face index or -1
          send_buffer[send_disp[rank0]+2] = fa_flag[ifa];
          send_buffer[send_disp[rank0]+3] = ifa+faora[mpi_rank]; // original global ifa - for future data.
          send_buffer[send_disp[rank0]+4] = nnof;
          send_disp[rank0] += 5;
          // the nodes get packed in their current r-h-rule order... 
          for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
          {
            send_buffer[send_disp[rank0]] = noofa_v[nof];
            send_disp[rank0] += 1;
          }
        }
      }
      else if (rank0==rank1)
      {
        // this is an internal face that will remain internal in the new
        // partition - both its cv's are in the same rank...
        assert(fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL);
        if (iter==0)
        {
          send_count[rank0] += 5+nnof;
        }
        else
        {
          send_buffer[send_disp[rank0]] = icv0;
          send_buffer[send_disp[rank0]+1] = icv1;
          send_buffer[send_disp[rank0]+2] = fa_flag[ifa];
          send_buffer[send_disp[rank0]+3] = ifa+faora[mpi_rank]; // original global ifa
          send_buffer[send_disp[rank0]+4] = nnof;
          send_disp[rank0] += 5;
          // the nodes get packed in their current r-h-rule order... 
          for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
          {
            send_buffer[send_disp[rank0]] = noofa_v[nof];
            send_disp[rank0] += 1;
          }
        }
      }
      else
      {
        // this is an internal face that will be on the inter-processor
        // boundary of the new partition - its cv's are different ranks.
        // here we need to send it to both.
        assert(rank1>=0);
        assert(fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL);
        if (iter==0)
        {
          send_count[rank0] += 5+nnof;
          send_count[rank1] += 5+nnof;
        }
        else
        {
          // rank0...
          send_buffer[send_disp[rank0]] = icv0;
          send_buffer[send_disp[rank0]+1] = ifa+faora[mpi_rank]; // put a global face index here
          send_buffer[send_disp[rank0]+2] = -fa_flag[ifa]-1; // negative face zone indicates IB
          send_buffer[send_disp[rank0]+3] = ifa+faora[mpi_rank]; // original global ifa
          send_buffer[send_disp[rank0]+4] = nnof;
          send_disp[rank0] += 5;
          // the nodes get packed in their current r-h-rule order... 
          for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
          {
            send_buffer[send_disp[rank0]] = noofa_v[nof];
            send_disp[rank0] += 1;
          }
          // rank1...
          send_buffer[send_disp[rank1]] = icv1;
          send_buffer[send_disp[rank1]+1] = ifa+faora[mpi_rank]; // global face index
          send_buffer[send_disp[rank1]+2] = -fa_flag[ifa]-1; // make this negative to indicate INTERNAL-BOUNDARY
          send_buffer[send_disp[rank1]+3] = -(ifa+faora[mpi_rank])-1; // neg global face index minus 1
          send_buffer[send_disp[rank1]+4] = nnof;
          send_disp[rank1] += 5;
          // the nodes get packed in the opposite direction...
          for (int nof = noofa_i[ifa+1]-1; nof>=noofa_i[ifa]; nof--)
          {
            send_buffer[send_disp[rank1]] = noofa_v[nof];
            send_disp[rank1] += 1;
          }
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  delete[] send_buffer;

  // rebuild the new cvofa, noofa_i, noofa_v, for the new partition...
  // also, store the zone in fa_zone (was fa_flag), and use fa_flag
  // to store the original global face index...

  int nfa_old = nfa;
  int * fa_zone;
  delete[] cvofa;
  delete[] fa_flag;
  delete[] noofa_i;
  delete[] noofa_v;
  for (int iter = 0; iter<3; iter++)
  {
    int ifa = 0;
    int i = 0;
    while (i<recv_count_sum)
    {
      int icv0 = recv_buffer[i++];
      int icv1 = recv_buffer[i++];
      int this_fa_zone = recv_buffer[i++]; // can be negative, indicating inter-processor
      int ifa_global = recv_buffer[i++]; // can be negative, indicating face flipping
      int nnof = recv_buffer[i++];
      if (iter==0)
      {
        // on the first iter, we are just counting
        // the new faces and there is no need to read the nodes...
        i += nnof;
      }
      else if (iter==1)
      {
        // on the second iter, we are counting nodes of face,
        // but there is still no reason to read the nodes...
        cvofa[ifa][0] = icv0;
        cvofa[ifa][1] = icv1;
        fa_zone[ifa] = this_fa_zone;
        fa_flag[ifa] = ifa_global;
        noofa_i[ifa+1] += nnof;
        i += nnof;
      }
      else
      {
        // on the third iter, read the nodes and build noofa_v...
        for (int nof = 0; nof<nnof; nof++)
        {
          noofa_v[noofa_i[ifa]] = recv_buffer[i++];
          noofa_i[ifa] += 1;
        }
      }
      // increment ifa...
      ifa++;
    }
    assert(i==recv_count_sum);
    if (iter==0)
    {
      nfa = ifa;
      fa_zone = new int[nfa];
      fa_flag = new int[nfa];
      cvofa = new int[nfa][2];
      noofa_i = new int[nfa+1];
      for (ifa = 0; ifa<nfa; ifa++)
        noofa_i[ifa+1] = 0;
    }
    else if (iter==1)
    {
      // make noofa_i CSR...
      noofa_i[0] = 0;
      for (ifa = 0; ifa<nfa; ifa++)
        noofa_i[ifa+1] += noofa_i[ifa];
      noofa_s = noofa_i[nfa];
      noofa_v = new int[noofa_s];
    }
    else
    {
      // return noofa_i to CSR...
      for (ifa = nfa; ifa>0; ifa--)
        noofa_i[ifa] = noofa_i[ifa-1];
      noofa_i[0] = 0;
    }
  }

  //cout << "mpi_rank, nfa_old, nfa: " << mpi_rank << " " << nfa_old << " " << nfa << endl;

  delete[] recv_buffer;

  // ----------------------------------------------
  // now it is time to collect the nodes. The current
  // node numbering on the faces is global node numbering,
  // so we can use noora[] to form requests for these nodes.
  // A critical problem, however, is how to deal with the 
  // node duplication currently present. Avoiding sorting
  // we can solve this problem in linear time with the
  // use of a bit more memory...
  // ----------------------------------------------

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // how many nodes...
      for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
      {
        int ino = noofa_v[nof]; // should be a global node index...
        assert((ino>=0)&&(ino<noora[mpi_size]));
        int rank;
        for (rank = 0; rank<mpi_size; rank++)
        {
          if (ino<noora[rank+1]) break;
        }
        if (iter==0)
        {
          send_count[rank] += 1;
        }
        else
        {
          send_buffer[send_disp[rank]] = ino;
          send_disp[rank] += 1;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // now cycle through our nodes, and return a modified recv_buffer 
  // indicating the nodal duplication within the rank...
  assert(no_flag!=NULL);
  for (int rank = 0; rank<mpi_size; rank++)
  {
    // we need to do this separately for each rank...
    for (int i = recv_disp[rank]; i<recv_disp[rank]+recv_count[rank]; i++)
    {
      int ino = recv_buffer[i]-noora[mpi_rank]; // local node numbering...
      assert((ino>=0)&&(ino<nno));
      no_flag[ino] = -1;
    }
    // now loop through the relevant range of the recv buffer...
    for (int i = recv_disp[rank]; i<recv_disp[rank]+recv_count[rank]; i++)
    {
      int ino = recv_buffer[i]-noora[mpi_rank]; // local node numbering...
      if (no_flag[ino]==-1)
      {
        // this node has not been visited before - store the
        // location of the node in the sub-list relative to the
        // start of the rank's requests...
        no_flag[ino] = i-recv_disp[rank];
        // set this position of the recv_buffer to
        // -1, indicating that a new node has been requested for the first time...
        recv_buffer[i] = -1;
      }
      else
      {
        // this node has already been requested by this processor...
        // return the request in the recv_buffer...
        recv_buffer[i] = no_flag[ino];
      }
    }
  }

  // now send the modified recv_buffer back to the calling process... 
  MPI_Alltoallv(recv_buffer, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);

  delete[] recv_buffer;

  // modify the send buffer to contain a new local node indexing...
  int nno_old = nno;
  delete[] no_flag;

  nno = 0;
  for (int rank = 0; rank<mpi_size; rank++)
  {
    // we need to do this separately for each rank...
    for (int i = send_disp[rank]; i<send_disp[rank]+send_count[rank]; i++)
    {
      if (send_buffer[i]==-1)
      {
        send_buffer[i] = nno++;
      }
      else
      {
        int i_prev = send_buffer[i]+send_disp[rank];
        assert(i_prev<i);
        send_buffer[i] = send_buffer[i_prev];
      }
    }
  }

  //cout << "mpi_rank, nno_old, nno: " << mpi_rank << " " << nno_old << " " << nno << endl;

  // reallocate the no_flag and set to -1... 
  no_flag = new int[nno];
  for (int ino = 0; ino<nno; ino++)
    no_flag[ino] = -1;

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    // how many nodes...
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof]; // should be a global node index...
      assert((ino>=0)&&(ino<noora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (ino<noora[rank+1]) break;
      }
      // pull the local node number from the modified send_buffer
      int ino_local = send_buffer[send_disp[rank]];
      send_disp[rank] += 1;

      // modify noofa_v to contain this new local indexing...
      assert((ino_local>=0)&&(ino_local<nno));
      noofa_v[nof] = ino_local;

      // store the old global indexing in the newly allocated no_flag...
      if (no_flag[ino_local]==-1)
      {
        no_flag[ino_local] = ino; // old global
      }
      else
      {
        assert(no_flag[ino_local]==ino);
      }
    }
  }

  // check no_flag was completely set...
  for (int ino = 0; ino<nno; ino++)
  {
    assert((no_flag[ino]>=0)&&(no_flag[ino]<noora[mpi_size]));
  }

  delete[] send_buffer;

  // =======================================
  // cv's...
  // =======================================

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // all faces hav a valid icv0, currently global...
      int icv0 = cvofa[ifa][0];
      assert((icv0>=0)&&(icv0<cvora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (icv0<cvora[rank+1]) break;
      }
      if (iter==0)
      {
        send_count[rank] += 1;
      }
      else
      {
        send_buffer[send_disp[rank]] = icv0;
        send_disp[rank] += 1;
      }
      // icv1 is valid if the cv is internal AND not inter-processor...
      if ((fa_zone[ifa]>=0)&&(fa_kind_table[fa_zone[ifa]]==FA_ZONE_INTERNAL))
      {
        int icv1 = cvofa[ifa][1];
        assert((icv1>=0)&&(icv1<cvora[mpi_size]));
        int rank;
        for (rank = 0; rank<mpi_size; rank++)
        {
          if (icv1<cvora[rank+1]) break;
        }
        if (iter==0)
        {
          send_count[rank] += 1;
        }
        else
        {
          send_buffer[send_disp[rank]] = icv1;
          send_disp[rank] += 1;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // on the receive side, we should own all the cv's being requested...
  // here we cannot use cv_flag because it contains the cv zone... 
  int * cv_tmp = new int[ncv];
  // now cycle through our nodes, and return a modified recv_buffer 
  // indicating the cv duplication within the rank...
  for (int rank = 0; rank<mpi_size; rank++)
  {
    // we need to do this separately for each rank...
    for (int i = recv_disp[rank]; i<recv_disp[rank]+recv_count[rank]; i++)
    {
      int icv = recv_buffer[i]-cvora[mpi_rank]; // local cv numbering...
      assert((icv>=0)&&(icv<ncv));
      cv_tmp[icv] = -1;
    }
    // now loop through the relevant range of the recv buffer...
    for (int i = recv_disp[rank]; i<recv_disp[rank]+recv_count[rank]; i++)
    {
      int icv = recv_buffer[i]-cvora[mpi_rank]; // local cv numbering...
      assert((icv>=0)&&(icv<ncv));
      if (cv_tmp[icv]==-1)
      {
        // this cv has not been visited before - store the
        // location of the cv in the sub-list relative to the
        // start of the rank's requests...
        cv_tmp[icv] = i-recv_disp[rank];
        // set this position of the recv_buffer to
        // -1, indicating that a new cv has been requested for the first time...
        recv_buffer[i] = -1;
      }
      else
      {
        // this cv has already been requested by this processor...
        // return the request in the recv_buffer...
        recv_buffer[i] = cv_tmp[icv];
      }
    }
  }
  delete[] cv_tmp;

  // now send the modified recv_buffer back to the calling process... 
  MPI_Alltoallv(recv_buffer, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);

  delete[] recv_buffer;

  // modify the send buffer to contain a new local cv indexing...
  int ncv_old = ncv;
  ncv = 0;
  for (int rank = 0; rank<mpi_size; rank++)
  {
    // we need to do this separately for each rank...
    for (int i = send_disp[rank]; i<send_disp[rank]+send_count[rank]; i++)
    {
      if (send_buffer[i]==-1)
      {
        send_buffer[i] = ncv++;
      }
      else
      {
        int i_prev = send_buffer[i]+send_disp[rank];
        assert(i_prev<i);
        send_buffer[i] = send_buffer[i_prev];
      }
    }
  }

  //cout << "mpi_rank, ncv_old, ncv: " << mpi_rank << " " << ncv_old << " " << ncv << endl;

  // reallocate the cv_tmp with the new ncv to store the old 
  // cv index - we cannot deallocate and reallocate cv_flag because it
  // still contains the cv zone...

  cv_tmp = new int[ncv];
  for (int icv = 0; icv<ncv; icv++)
    cv_tmp[icv] = -1;

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    // all faces have a valid icv0, currently global...
    int icv0 = cvofa[ifa][0];
    assert((icv0>=0)&&(icv0<cvora[mpi_size]));
    int rank;
    for (rank = 0; rank<mpi_size; rank++)
    {
      if (icv0<cvora[rank+1]) break;
    }
    // the send buffer now has the local numbering in it...
    int icv0_local = send_buffer[send_disp[rank]];
    assert((icv0_local>=0)&&(icv0_local<ncv));
    send_disp[rank] += 1;
    // go for local indexing... 
    cvofa[ifa][0] = icv0_local;
    // and store the old global indexing...
    if (cv_tmp[icv0_local]==-1)
    {
      cv_tmp[icv0_local] = icv0; // old global
    }
    else
    {
      assert(cv_tmp[icv0_local]==icv0);
    }
    // icv1 is valid if the cv is internal AND not inter-processor...
    if ((fa_zone[ifa]>=0)&&(fa_kind_table[fa_zone[ifa]]==FA_ZONE_INTERNAL))
    {
      int icv1 = cvofa[ifa][1];
      assert((icv1>=0)&&(icv1<cvora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (icv1<cvora[rank+1]) break;
      }
      // the send buffer now has the local numbering in it...
      int icv1_local = send_buffer[send_disp[rank]];
      assert((icv1_local>=0)&&(icv1_local<ncv));
      send_disp[rank] += 1;
      // go for local indexing... 
      cvofa[ifa][1] = icv1_local;
      // and store the old global indexing...
      if (cv_tmp[icv1_local]==-1)
      {
        cv_tmp[icv1_local] = icv1; // old global
      }
      else
      {
        assert(cv_tmp[icv1_local]==icv1);
      }
    }
  }

  // check cv_tmp was completely set...
  for (int icv = 0; icv<ncv; icv++)
  {
    assert((cv_tmp[icv]>=0)&&(cv_tmp[icv]<cvora[mpi_size]));
  }

  delete[] send_buffer;

  // -------------------------------------------------------------------
  // At this point, the fa_flag contains the original global face index
  // for each new face (the zone index is now in the temp array fa_zone).
  // Global index fa_flag can be -ve, indicating a flipped 
  // face orientation with respect to the original face. This is
  // ofcourse important for getting the sign of directional data associated 
  // with faces correct.
  //
  // The no_flag contains the original global indexing for nodes.
  // Nodes had no original "zone", so there is no analogous no_zone.
  // 
  // And cv_tmp contains the original cv global indexing for cv's, and the
  // original cv zones has not yet been transfered.
  // -------------------------------------------------------------------

  // ===================================================================
  // complete cv zone transfer and do data exchanges...
  // ===================================================================

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  int * send_buffer2;
  for (int iter = 0; iter<2; iter++)
  {
    for (int icv = 0; icv<ncv; icv++)
    {
      int icv_global = cv_tmp[icv];
      assert((icv_global>=0)&&(icv_global<cvora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (icv_global<cvora[rank+1]) break;
      }
      if (iter==0)
      {
        send_count[rank] += 1;
      }
      else
      {
        send_buffer[send_disp[rank]] = icv_global;
        send_buffer2[send_disp[rank]] = icv;
        send_disp[rank] += 1;
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
      send_buffer2 = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // -------------------------
  // 1. cv_part as a check...
  // -------------------------
  int * recv_buffer2 = new int[recv_count_sum];
  for (int i = 0; i<recv_count_sum; i++)
  {
    int icv = recv_buffer[i]-cvora[mpi_rank];
    assert((icv>=0)&&(icv<ncv_old));
    recv_buffer2[i] = cv_part[icv];
  }
  MPI_Alltoallv(recv_buffer2, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);
  for (int i = 0; i<send_count_sum; i++)
  {
    int icv = send_buffer2[i]; // recall we put the local cv number in send_buffer2...
    assert((icv>=0)&&(icv<ncv));
    assert(send_buffer[i]==mpi_rank);
  }

  // -------------------------
  // 2. old cv_flag contains the cv_zone...
  // -------------------------
  for (int i = 0; i<recv_count_sum; i++)
  {
    int icv = recv_buffer[i]-cvora[mpi_rank];
    assert((icv>=0)&&(icv<ncv_old));
    recv_buffer2[i] = cv_flag[icv];
  }
  delete[] cv_flag;
  MPI_Alltoallv(recv_buffer2, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);
  int * cv_zone = new int[ncv];
  for (int i = 0; i<send_count_sum; i++)
  {
    int icv = send_buffer2[i];
    assert((icv>=0)&&(icv<ncv));
    cv_zone[icv] = send_buffer[i];
    // expect one cv zone for now...
    //assert( cv_zone[icv] == 0 );
  }
  // and point cv_flag to cv_tmp
  cv_flag = cv_tmp;

  // so now, like the faces, cv_flag has the global cv indexing, and cv_zone
  // has the cv zone numbering -- just one zone for now, but can be 
  // any number in the future.

  // registered cv scalars...
  {
    double * recv_buffer_double = NULL;
    double * send_buffer_double = NULL;
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if (data->getDatatype()==CV_DATA)
      {
        if (recv_buffer_double==NULL) recv_buffer_double = new double[recv_count_sum];
        for (int i = 0; i<recv_count_sum; i++)
        {
          int icv = recv_buffer[i]-cvora[mpi_rank];
          assert((icv>=0)&&(icv<ncv_old));
          recv_buffer_double[i] = (*data->ptr)[icv];
        }
        if (send_buffer_double==NULL) send_buffer_double = new double[send_count_sum];
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[ncv];
        for (int i = 0; i<send_count_sum; i++)
        {
          int icv = send_buffer2[i];
          assert((icv>=0)&&(icv<ncv));
          (*data->ptr)[icv] = send_buffer_double[i];
        }
      }
    }
    if (recv_buffer_double!=NULL) delete[] recv_buffer_double;
    if (send_buffer_double!=NULL) delete[] send_buffer_double;
  }

  // registered cv vectors...
  for (int i = 0; i<mpi_size; i++)
  {
    recv_count[i] *= 3;
    recv_disp[i] *= 3;
    send_count[i] *= 3;
    send_disp[i] *= 3;
  }

  {
    double * recv_buffer_double = NULL;
    double * send_buffer_double = NULL;
    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
    {
      if (data->getDatatype()==CV_DATA)
      {
        if (recv_buffer_double==NULL) recv_buffer_double = new double[3*recv_count_sum];
        for (int i = 0; i<recv_count_sum; i++)
        {
          int icv = recv_buffer[i]-cvora[mpi_rank];
          assert((icv>=0)&&(icv<ncv_old));
          for (int j = 0; j<3; j++)
            recv_buffer_double[3*i+j] = (*data->ptr)[icv][j];
        }
        if (send_buffer_double==NULL) send_buffer_double = new double[3*send_count_sum];
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[ncv][3];
        for (int i = 0; i<send_count_sum; i++)
        {
          int icv = send_buffer2[i];
          assert((icv>=0)&&(icv<ncv));
          for (int j = 0; j<3; j++)
            (*data->ptr)[icv][j] = send_buffer_double[3*i+j];
        }
      }
    }
    if (recv_buffer_double!=NULL) delete[] recv_buffer_double;
    if (send_buffer_double!=NULL) delete[] send_buffer_double;
  }

  // registered cv tensors...
  for (int i = 0; i<mpi_size; i++)
  {
    recv_count[i] *= 3;
    recv_disp[i] *= 3;
    send_count[i] *= 3;
    send_disp[i] *= 3;
  }

  {
    double * recv_buffer_double = NULL;
    double * send_buffer_double = NULL;
    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
    {
      if (data->getDatatype()==CV_DATA)
      {
        if (recv_buffer_double==NULL) recv_buffer_double = new double[9*recv_count_sum];
        for (int i = 0; i<recv_count_sum; i++)
        {
          int icv = recv_buffer[i]-cvora[mpi_rank];
          assert((icv>=0)&&(icv<ncv_old));
          for (int j = 0; j<3; j++)
            for (int k = 0; k<3; k++)
              recv_buffer_double[9*i+3*j+k] = (*data->ptr)[icv][j][k];
        }
        if (send_buffer_double==NULL) send_buffer_double = new double[9*send_count_sum];
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[ncv][3][3];
        for (int i = 0; i<send_count_sum; i++)
        {
          int icv = send_buffer2[i];
          assert((icv>=0)&&(icv<ncv));
          for (int j = 0; j<3; j++)
            for (int k = 0; k<3; k++)
              (*data->ptr)[icv][j][k] = send_buffer_double[9*i+3*j+k];
        }
      }
    }
    if (recv_buffer_double!=NULL) delete[] recv_buffer_double;
    if (send_buffer_double!=NULL) delete[] send_buffer_double;
  }

  delete[] recv_buffer;
  delete[] recv_buffer2;
  delete[] send_buffer;
  delete[] send_buffer2;

  // ===================================================================
  // nodal data exchanges...
  // ===================================================================

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ino = 0; ino<nno; ino++)
    {
      int ino_global = no_flag[ino];
      assert((ino_global>=0)&&(ino_global<noora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (ino_global<noora[rank+1]) break;
      }
      if (iter==0)
      {
        send_count[rank] += 1;
      }
      else
      {
        send_buffer[send_disp[rank]] = ino_global;
        send_buffer2[send_disp[rank]] = ino;
        send_disp[rank] += 1;
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
      send_buffer2 = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  delete[] send_buffer;

  // allocate the buffers required for the vector sends - we always have to
  // send x_no, so might as well do it all at once...

  // do NOT get confused - recv (the original nodal striping) is sending to 
  // send (the final partitioning)...

  {

    double * recv_buffer_double = new double[3*recv_count_sum];
    double * send_buffer_double = new double[3*send_count_sum];

    // ----------------------
    // scalar nodal data...
    // ----------------------

    // registered nodal scalars...
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if (data->getDatatype()==NO_DATA)
      {
        for (int i = 0; i<recv_count_sum; i++)
        {
          int ino = recv_buffer[i]-noora[mpi_rank];
          recv_buffer_double[i] = (*data->ptr)[ino];
        }
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[nno];
        for (int i = 0; i<send_count_sum; i++)
        {
          int ino = send_buffer2[i];
          (*data->ptr)[ino] = send_buffer_double[i];
        }
      }
    }

    // ----------------------
    // vector nodal data...
    // ----------------------

    // multiply counts and disps by 3...
    for (int i = 0; i<mpi_size; i++)
    {
      recv_count[i] *= 3;
      recv_disp[i] *= 3;
      send_count[i] *= 3;
      send_disp[i] *= 3;
    }

    // x_no...
    for (int i = 0; i<recv_count_sum; i++)
    {
      int ino = recv_buffer[i]-noora[mpi_rank];
      assert((ino>=0)&&(ino<nno_old));
      for (int j = 0; j<3; j++)
        recv_buffer_double[3*i+j] = x_no[ino][j];
    }
    MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
        MPI_DOUBLE, mpi_comm);
    delete[] x_no;
    x_no = new double[nno][3];
    for (int i = 0; i<send_count_sum; i++)
    {
      // recall we put the local node number in send_buffer2...
      int ino = send_buffer2[i];
      for (int j = 0; j<3; j++)
        x_no[ino][j] = send_buffer_double[3*i+j];
    }

    // registered nodal vectors...
    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
    {
      if (data->getDatatype()==NO_DATA)
      {
        for (int i = 0; i<recv_count_sum; i++)
        {
          int ino = recv_buffer[i]-noora[mpi_rank];
          for (int j = 0; j<3; j++)
            recv_buffer_double[3*i+j] = (*data->ptr)[ino][j];
        }
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[nno][3];
        for (int i = 0; i<send_count_sum; i++)
        {
          int ino = send_buffer2[i];
          for (int j = 0; j<3; j++)
            (*data->ptr)[ino][j] = send_buffer_double[3*i+j];
        }
      }
    }

    // ----------------------
    // nodal tensor data
    // ----------------------

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
    {
      if (data->getDatatype()==NO_DATA)
      {
        cout<<"not setup to handle nodal tensor reorder here"<<endl;
        throw(-1);
      }
    }

    delete[] send_buffer_double;
    delete[] recv_buffer_double;

  }

  // done nodal data...

  delete[] send_buffer2;
  delete[] recv_buffer;

  // ===================================================================
  // face data exchanges...
  // ===================================================================

  // recall that the face global indexing in fa_flag is 
  // signed, to indicate when flipping of direction data is required.

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      int ifa_global = max(fa_flag[ifa], -fa_flag[ifa]-1);
      assert((ifa_global>=0)&&(ifa_global<faora[mpi_size]));
      int rank;
      for (rank = 0; rank<mpi_size; rank++)
      {
        if (ifa_global<faora[rank+1]) break;
      }
      if (iter==0)
      {
        send_count[rank] += 1;
      }
      else
      {
        send_buffer[send_disp[rank]] = ifa_global;
        send_buffer2[send_disp[rank]] = ifa;
        send_disp[rank] += 1;
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
      send_buffer2 = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  delete[] send_buffer;

  // do NOT get confused - recv (the original nodal striping) is sending to 
  // send (the final partitioning)...

  {

    double * recv_buffer_double = NULL;
    double * send_buffer_double = NULL;

    // ----------------------
    // scalar face data...
    // ----------------------

    // registered scalars...
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if (data->getDatatype()==FA_DATA)
      {
        // allocate...
        if (recv_buffer_double==NULL) recv_buffer_double = new double[recv_count_sum];
        if (send_buffer_double==NULL) send_buffer_double = new double[send_count_sum];
        // pack...
        for (int i = 0; i<recv_count_sum; i++)
        {
          int ifa = recv_buffer[i]-faora[mpi_rank];
          recv_buffer_double[i] = (*data->ptr)[ifa];
        }
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[nfa];
        for (int i = 0; i<send_count_sum; i++)
        {
          int ifa = send_buffer2[i];
          // recall face data is signed by default...
          if (fa_flag[ifa]>=0) (*data->ptr)[ifa] = send_buffer_double[i];
          else (*data->ptr)[ifa] = -send_buffer_double[i];
        }
      }
    }

    if (recv_buffer_double!=NULL)
    {
      delete[] recv_buffer_double;
      recv_buffer_double = NULL;
    }
    if (send_buffer_double!=NULL)
    {
      delete[] send_buffer_double;
      send_buffer_double = NULL;
    }

    // ----------------------
    // vector face data...
    // ----------------------

    // multiply counts and disps by 3...
    for (int i = 0; i<mpi_size; i++)
    {
      recv_count[i] *= 3;
      recv_disp[i] *= 3;
      send_count[i] *= 3;
      send_disp[i] *= 3;
    }

    // registered face vectors...
    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
    {
      if (data->getDatatype()==FA_DATA)
      {
        // allocate...
        if (recv_buffer_double==NULL) recv_buffer_double = new double[3*recv_count_sum];
        if (send_buffer_double==NULL) send_buffer_double = new double[3*send_count_sum];
        // pack...
        for (int i = 0; i<recv_count_sum; i++)
        {
          int ifa = recv_buffer[i]-faora[mpi_rank];
          for (int j = 0; j<3; j++)
            recv_buffer_double[3*i+j] = (*data->ptr)[ifa][j];
        }
        MPI_Alltoallv(recv_buffer_double, recv_count, recv_disp, MPI_DOUBLE, send_buffer_double, send_count, send_disp,
            MPI_DOUBLE, mpi_comm);
        delete[] (*data->ptr);
        *data->ptr = new double[nfa][3];
        for (int i = 0; i<send_count_sum; i++)
        {
          int ifa = send_buffer2[i];
          // recall face data is signed by default...
          if (fa_flag[ifa]>=0) for (int j = 0; j<3; j++)
            (*data->ptr)[ifa][j] = send_buffer_double[3*i+j];
          else for (int j = 0; j<3; j++)
            (*data->ptr)[ifa][j] = -send_buffer_double[3*i+j];
        }
      }
    }

    if (recv_buffer_double!=NULL)
    {
      delete[] recv_buffer_double;
      recv_buffer_double = NULL;
    }
    if (send_buffer_double!=NULL)
    {
      delete[] send_buffer_double;
      send_buffer_double = NULL;
    }

    // ----------------------
    // face tensor data - not yet
    // ----------------------

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
    {
      if (data->getDatatype()==FA_DATA)
      {
        cout<<"not setup to handle face tensor reorder here"<<endl;
        throw(-1);
      }
    }

    // done face data...

  }

  delete[] send_buffer2;
  delete[] recv_buffer;

  // ===================================================================
  // face connectivity...
  // ===================================================================

  for (int i = 0; i<mpi_size; i++)
    send_count[i] = 0;
  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // pack the global face number if this is a ib/periodic face...
      if (fa_zone[ifa]<0)
      {
        assert(fa_kind_table[-fa_zone[ifa]-1]==FA_ZONE_INTERNAL);
        // note that a face from an internal boundary can be positive or negative,
        // depending on its default outward-pointing orientation on the 
        // new processor relative to its original orientation on the 
        // original global face distribution.
        int ifa_global = max(fa_flag[ifa], -fa_flag[ifa]-1);
        int rank = 0;
        while (ifa_global>=faora[rank+1])
          rank++;
        if (iter==0)
        {
          send_count[rank] += 1;
        }
        else
        {
          send_buffer[send_disp[rank]] = fa_flag[ifa]; // send the signed value
          send_disp[rank] += 1;
        }
      }
      else if ((fa_kind_table[fa_zone[ifa]]>=FA_ZONE_PERIODIC_FIRST)&&(fa_kind_table[fa_zone[ifa]]
          <=FA_ZONE_PERIODIC_LAST))
      {
        // recall that periodic faces have their cvofa[][1] == ifa nbr
        int ifa_global = cvofa[ifa][1];
        assert((ifa_global>=0)&&(ifa_global<faora[mpi_size]));
        int rank = 0;
        while (ifa_global>=faora[rank+1])
          rank++;
        if (iter==0)
        {
          send_count[rank] += 1;
        }
        else
        {
          send_buffer[send_disp[rank]] = ifa_global;
          send_disp[rank] += 1;
        }
      }
    }
    // calculate send_disp on both iterations...
    send_disp[0] = 0;
    for (int i = 1; i<mpi_size; i++)
      send_disp[i] = send_count[i-1]+send_disp[i-1];
    // on the first time through, allocate the send_buffer...
    if (iter==0)
    {
      send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
      send_buffer = new int[send_count_sum];
    }
  }

  // now build the recv side stuff and exchange...
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
  recv_disp[0] = 0;
  for (int i = 1; i<mpi_size; i++)
    recv_disp[i] = recv_count[i-1]+recv_disp[i-1];
  recv_count_sum = recv_disp[mpi_size-1]+recv_count[mpi_size-1];
  recv_buffer = new int[recv_count_sum];
  MPI_Alltoallv(send_buffer, send_count, send_disp, MPI_INT, recv_buffer, recv_count, recv_disp, MPI_INT, mpi_comm);

  // on the recv side, determine where the requested faces were sent in the 
  // final partition and send back the rank...
  for (int i = 0; i<recv_count_sum; i++)
  {
    int ifa_global = recv_buffer[i];
    if (ifa_global<0)
    {
      // a negative global face index indicates an internal boundary face 
      // flipped on the new processor, so the new face should be...
      int ifa_local = -ifa_global-1-faora[mpi_rank];
      assert((ifa_local>=0)&&(ifa_local<nfa_old));
      // both ranks should be valid AND different...
      int rank0 = raocvofa[ifa_local][0];
      int rank1 = raocvofa[ifa_local][1];
      assert((rank0>=0)&&(rank0<mpi_size));
      assert((rank1>=0)&&(rank1<mpi_size));
      assert(rank0!=rank1);
      // because this is the flipped face, the rank of the face
      // pair we are searching for is the first rank: rank0
      recv_buffer[i] = rank0;
    }
    else
    {
      // this face could either be internal unflipped or periodic...
      // if it is periodic, then the raocvofa[ifa][1] should be -1, 
      // and we need to return the rank0. If it is an internal face with
      // both ranks valid in its cv's, then it is an internal face that
      // has not been flipped, and we need to return rank1... 
      int ifa_local = ifa_global-faora[mpi_rank];
      assert((ifa_local>=0)&&(ifa_local<nfa_old));
      int rank1 = raocvofa[ifa_local][1];
      if (rank1>=0)
      {
        recv_buffer[i] = rank1;
      }
      else
      {
        int rank0 = raocvofa[ifa_local][0];
        assert((rank0>=0)&&(rank0<mpi_size));
        recv_buffer[i] = rank0;
      }
    }
  }

  // done with this...
  delete[] raocvofa;

  // note that it would also be possible to reduce the local face indexing
  // to the old face partition with two ints per internal face, and we
  // could return from this routine with a rank and local face index on
  // that rank to determine the face matching. for now, I prefer to 
  // send a request to the processor, perhaps even after ordering, with
  // the global face indexing, and the processor can then build the 
  // pack and unpack structures with a bit of sorting to get the
  // face correspondence correct...

  MPI_Alltoallv(recv_buffer, recv_count, recv_disp, MPI_INT, send_buffer, send_count, send_disp, MPI_INT, mpi_comm);

  delete[] recv_buffer;

  // now the send buffer has the ranks of the internal boundary and 
  // periodic faces where inter (or intra, for the case of periodic self-
  // matching) processor faces are located, along with global indexing.

  // take these out of the send buffer and put them into a rank of face...
  int * fa_nbr_rank = new int[nfa];
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_nbr_rank[ifa] = -1;
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    // pack the global face number if this is a ib/periodic face...
    if (fa_zone[ifa]<0)
    {
      assert(fa_kind_table[-fa_zone[ifa]-1]==FA_ZONE_INTERNAL);
      // note that a face from an internal boundary can be positive or negative,
      // depending on its default outward-pointing orientation on the 
      // new processor relative to its original orientation on the 
      // original global face distribution.
      int ifa_global = max(fa_flag[ifa], -fa_flag[ifa]-1);
      int rank = 0;
      while (ifa_global>=faora[rank+1])
        rank++;
      fa_nbr_rank[ifa] = send_buffer[send_disp[rank]];
      send_disp[rank] += 1;
      assert((fa_nbr_rank[ifa]>=0)&&(fa_nbr_rank[ifa]<mpi_size));
      assert(fa_nbr_rank[ifa]!=mpi_rank); // IB cannot be self-matching
    }
    else if ((fa_kind_table[fa_zone[ifa]]>=FA_ZONE_PERIODIC_FIRST)
        &&(fa_kind_table[fa_zone[ifa]]<=FA_ZONE_PERIODIC_LAST))
    {
      // recall that periodic faces have their cvofa[][1] == ifa nbr
      int ifa_global = cvofa[ifa][1];
      assert((ifa_global>=0)&&(ifa_global<faora[mpi_size]));
      int rank = 0;
      while (ifa_global>=faora[rank+1])
        rank++;
      fa_nbr_rank[ifa] = send_buffer[send_disp[rank]];
      send_disp[rank] += 1;
      assert((fa_nbr_rank[ifa]>=0)&&(fa_nbr_rank[ifa]<mpi_size));
      // self-matching is possible for periodic
    }
  }
  delete[] send_buffer;
  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  // now we should be able to build the paired communicators for faces
  // using the Gp constructors...
  {
    // cycle through all faces and count the number of possible ranges,
    // as well as the number of ranks we need to consider. For a good 
    // partitioning, this will be small compared to the total number of ranks,
    // making this operation more scalable, even for 1000's of processors...

    int * rank_count = new int[mpi_size];
    for (int rank = 0; rank<mpi_size; rank++)
      rank_count[rank] = 0;
    int np = 0; // number of face pairs
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_nbr_rank[ifa]>=0)
      {
        rank_count[fa_nbr_rank[ifa]] += 1;
        np++;
      }
    }

    // as a check, face exchanges should be symmetric... i.e. every sent
    // request should have a matching face that has just sent to us...
    int * check_count = new int[mpi_size];
    MPI_Alltoall(rank_count, 1, MPI_INT, check_count, 1, MPI_INT, mpi_comm);
    for (int rank = 0; rank<mpi_size; rank++)
    {
      assert(rank_count[rank]==check_count[rank]);
    }
    delete[] check_count;

    // passed symmetric check - neccessary but not sufficient. 

    // build and sort a list of local/global face pairings that includes 
    // just the faces associated with exchanges - improves the sort time...
    vector<IntPair> pairVec;
    if (np>0)
    {
      pairVec.resize(np);
      np = 0;
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_nbr_rank[ifa]>=0)
        {
          pairVec[np].i2 = ifa;
          // the sorted index is the first: i1, which should be
          // based on the global face index in this case...
          pairVec[np].i1 = max(fa_flag[ifa], -fa_flag[ifa]-1);
          np++;
        }
      }
      // sort using i1...
      std::sort(pairVec.begin(), pairVec.end(), IntPairCompare1());
    }

    // Now build the face prcomm for each range...    

    // start by setting the "bits" associated with each face range. 
    // note that face ranges are synchronized across processors, so this 
    // will be the same on all ranks...

    // here we sort the faZoneList for aesthetic reasons only :)
    faZoneList.sort();

    assert(periodic_bit_table==NULL);
    n_periodic_bit_table = 0;
    list<FaZone>::iterator iz;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
    {
      if ((iz->getKind()>=FA_ZONE_PERIODIC_FIRST)&&(iz->getKind()<=FA_ZONE_PERIODIC_LAST))
      {
        assert(n_periodic_bit_table<31); // if this fails, there is way too much periodicity
        iz->bits = 1<<n_periodic_bit_table;
        n_periodic_bit_table++;
      }
      else
      {
        iz->bits = 0;
      }
    }
    int * my_periodic_bit_table = NULL;
    if (n_periodic_bit_table>0)
    {
      periodic_bit_table = new int[n_periodic_bit_table];
      my_periodic_bit_table = new int[n_periodic_bit_table];
      for (int i = 0; i<n_periodic_bit_table; i++)
        my_periodic_bit_table[i] = -1;
    }

    // look at the sorted face zones...
    if (mpi_rank==0)
    {
      cout<<"faZoneList dump after sort: "<<endl;
      for (list<FaZone>::iterator iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
      {
        iz->dump();
      }
    }

    // now add face paired communicators for each non-zero rank_count...

    for (int rank = 0; rank<mpi_size; rank++)
    {
      if (rank_count[rank]>0)
      {
        Prcomm * prcomm = getFacePrcomm(rank);
        // size the pack and unpack vecs right away...
        prcomm->packIndexVec.resize(rank_count[rank]);
        prcomm->unpackIndexVec.resize(rank_count[rank]);
        // the packBufferInt gets 2 times this size because it will pass the global
        // face index and the bits for each face...
        int npack = 2*rank_count[rank];
        int * ptr;
        if ((ptr = (int*) realloc(prcomm->packBufferInt, sizeof(int)*npack))==NULL)
        {
          cerr<<"Error: realloc failed in initial build"<<endl;
          throw(-1);
        }
        prcomm->packBufferInt = ptr;
        prcomm->packBufferIntSize = npack;
        // now build ranges - internal first...
        int i_f = 0;
        int i_l = i_f;
        for (int ifa = 0; ifa<nfa; ifa++)
        {
          if (fa_nbr_rank[ifa]==rank)
          {
            // this is a face that requires exchange with rank. Check if
            // it is an internal face that has been split. This can be determined
            // if its zone is an internal zone. 
            if (fa_zone[ifa]<0)
            {
              assert(mpi_rank!=rank); // no self-passing of non-periodic data.
              // this must be an internal zone...
              assert(fa_kind_table[-fa_zone[ifa]-1]==FA_ZONE_INTERNAL);
              prcomm->packIndexVec[i_l] = ifa;
              prcomm->packBufferInt[2*i_l] = max(fa_flag[ifa], -fa_flag[ifa]-1); // get the positive global index
              prcomm->packBufferInt[2*i_l+1] = 0; // bits = 0 for internal processor boundary faces.
              i_l++;
            }
            else
            {
              // any other should be periodic...
              assert((fa_kind_table[fa_zone[ifa]]>=FA_ZONE_PERIODIC_FIRST)&&(fa_kind_table[fa_zone[ifa]]
                  <=FA_ZONE_PERIODIC_LAST));
            }
          }
        }
        // if there were any internals added, then add a new internal range...
        if (i_l>i_f)
        {
          prcomm->addPackRange(i_f, i_l-1, 0, FA_ZONE_INTERNAL, NULL);
          i_f = i_l;
        }
        // now do the periodic ranges. Each periodic range has bits. The bits
        // are set to be unique per periodic face zone, but with no specific ordering
        // i.e. the precise pairing of periodic faces is not known, but can be worked
        // out because we have the global face index pairs... 
        for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
        {
          if ((iz->getKind()>=FA_ZONE_PERIODIC_FIRST)&&(iz->getKind()<=FA_ZONE_PERIODIC_LAST))
          {
            // this is a periodic range. Count if there are any faces in it...
            for (int ifa = 0; ifa<nfa; ifa++)
            {
              if ((fa_nbr_rank[ifa]==rank)&&(fa_zone[ifa]==iz->getIndex()))
              {
                prcomm->packIndexVec[i_l] = ifa;
                assert((cvofa[ifa][1]>=0)&&(cvofa[ifa][1]<faora[mpi_size]));
                prcomm->packBufferInt[2*i_l] = cvofa[ifa][1]; // the positive global index is in cvofa[ifa][1] 
                assert(iz->bits!=0);
                prcomm->packBufferInt[2*i_l+1] = iz->bits; // unique bits set in each periodic zone
                i_l++;
              }
            }
          }
          // if there were any added, then add a new periodic range...
          if (i_l>i_f)
          {
            prcomm->addPackRange(i_f, i_l-1, iz->bits, iz->getKind(), iz->periodic_data);
            i_f = i_l;
          }
        }
        // we should have accounted for all faces...
        assert(i_l==rank_count[rank]);
      }
    }
    delete[] rank_count;

    // now exchange the face buffers - here we want to put the current global face indices
    // and bits from packBufferInt to unpackBufferInt...

    int requestCount = 0;
    list<Prcomm>::iterator prcomm;
    for (prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
      if (prcomm->getNbrRank()!=mpi_rank) requestCount++;
    MPI_Request * sendRequestArray;
    MPI_Request * recvRequestArray;
    MPI_Status * statusArray;
    if (requestCount>0)
    {
      sendRequestArray = new MPI_Request[requestCount];
      recvRequestArray = new MPI_Request[requestCount];
      statusArray = new MPI_Status[requestCount];
    }
    requestCount = 0;
    for (prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
    {
      // ensure size of unpackBufferInt...
      int nunpack = 2*prcomm->unpackIndexVec.size();
      if (nunpack>prcomm->unpackBufferIntSize)
      {
        int * ptr;
        if ((ptr = (int*) realloc(prcomm->unpackBufferInt, sizeof(int)*nunpack))==NULL)
        {
          cerr<<"Error: realloc failed"<<endl;
          throw(-1);
        }
        prcomm->unpackBufferInt = ptr;
        prcomm->unpackBufferIntSize = nunpack;
      }
      // post the Irecv...
      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Irecv(prcomm->unpackBufferInt, nunpack, MPI_INT, prcomm->getNbrRank(), 5000, mpi_comm,
            &(recvRequestArray[requestCount]));
      }
      // assert size of packBufferInt...
      int npack = 2*prcomm->packIndexVec.size();
      assert(npack==prcomm->packBufferIntSize);
      // already packed, so send...
      if (prcomm->getNbrRank()==mpi_rank)
      {
        assert(npack==nunpack);
        for (int i = 0; i<npack; i++)
          prcomm->unpackBufferInt[i] = prcomm->packBufferInt[i];
      }
      else
      {
        MPI_Issend(prcomm->packBufferInt, npack, MPI_INT, prcomm->getNbrRank(), 5000, mpi_comm,
            &(sendRequestArray[requestCount]));
        // inc request count...
        requestCount++;
      }
    }
    // now we wait for all messages to be received...
    if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);
    // now unpack and figure out face matches...
    vector<IntTriple> tripleVec;
    for (prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
    {
      int nunpack_half = prcomm->unpackIndexVec.size();
      assert(nunpack_half>0);
      tripleVec.resize(nunpack_half);
      for (int i = 0; i<nunpack_half; i++)
      {
        tripleVec[i].i2 = i; // use i2 to store the original position...
        tripleVec[i].i1 = prcomm->unpackBufferInt[2*i]; // the global face index -- i1 is sorted on
        tripleVec[i].i3 = prcomm->unpackBufferInt[2*i+1]; // bits of forward transform
      }
      // sort using i1...
      std::sort(tripleVec.begin(), tripleVec.end(), IntTripleCompare1());
      int ii = 0;
      for (int i = 0; i<tripleVec.size(); i++)
      {
        while (tripleVec[i].i1!=pairVec[ii].i1)
        {
          ii++;
          assert(ii<pairVec.size());
        }
        prcomm->unpackIndexVec[tripleVec[i].i2] = pairVec[ii].i2; // the local index of this global face
        prcomm->unpackBufferInt[tripleVec[i].i2] = tripleVec[i].i3; // bits, so we can build unpack ranges
      }
    }
    // cleanup after sends...
    if (requestCount>0)
    {
      MPI_Waitall(requestCount, sendRequestArray, statusArray);
      delete[] sendRequestArray;
      delete[] recvRequestArray;
      delete[] statusArray;
    }

    // reset cvofa[ifa][1] to -1 for all boundary cv's...
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_nbr_rank[ifa]>=0)
      {
        cvofa[ifa][1] = -1;
      }
      // clear it...
      fa_nbr_rank[ifa] = -1;
    }

    // and build the unpack ranges based on the bits...
    for (prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
    {
      // put the pack range bits into fa_nbr_rank...
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int ifa = prcomm->packIndexVec[i];
          // every face should appear once and only once in all exchanges...
          assert(fa_nbr_rank[ifa]==-1);
          // use fa_nbr_rank to store the pack range bits...
          fa_nbr_rank[ifa] = ri->getBits();
        }
      }
      // now check unpackIndexVec and build unpack ranges...
      int i_l = 0;
      while (i_l<prcomm->unpackIndexVec.size())
      {
        int i_f = i_l;
        int unpack_bits = prcomm->unpackBufferInt[i_l];
        int ifa = prcomm->unpackIndexVec[i_l];
        // this face should be one of the faces also in the packBuffer...
        assert(fa_nbr_rank[ifa]>=0);
        int pack_bits = fa_nbr_rank[ifa];
        while ((i_l<prcomm->unpackIndexVec.size())&&(prcomm->unpackBufferInt[i_l]==unpack_bits))
        {
          // check...
          ifa = prcomm->unpackIndexVec[i_l];
          // this face should be one of the faces also in the packBuffer...
          assert(fa_nbr_rank[ifa]>=0);
          assert(fa_nbr_rank[ifa]==pack_bits);
          i_l++;
        }
        // we have an unpack range that matches a pack range...
        //cout << "face range: " << i_f << ":" << i_l-1 << " bits: " << pack_bits << ":" << unpack_bits << endl;
        assert(i_l>i_f);
        if (unpack_bits==0)
        {
          assert(pack_bits==0);
          assert(i_f==0); // should be the first
          prcomm->addUnpackRange(i_f, i_l-1, 0, FA_ZONE_INTERNAL, NULL);

        }
        else
        {
          assert(pack_bits!=0);
          assert(pack_bits!=unpack_bits);
          // we can now insert or check 2 entries in the periodic_bit_table...
          // note that we are actually setting /checking my_periodic_bit_table, because
          // not all processors will be able to fill any/all entries, so we need to reduce
          // below - a final synchronization chck...
          int pack_bit;
          for (pack_bit = 0; pack_bit<n_periodic_bit_table; pack_bit++)
          {
            if (pack_bits==(1<<pack_bit)) break;
          }
          assert(pack_bit<n_periodic_bit_table);
          int unpack_bit;
          for (unpack_bit = 0; unpack_bit<n_periodic_bit_table; unpack_bit++)
          {
            if (unpack_bits==(1<<unpack_bit)) break;
          }
          assert(unpack_bit<n_periodic_bit_table);
          // set/check table entries...
          if (my_periodic_bit_table[pack_bit]==-1)
          {
            my_periodic_bit_table[pack_bit] = unpack_bit;
          }
          else
          {
            assert(my_periodic_bit_table[pack_bit]==unpack_bit);
          }
          if (my_periodic_bit_table[unpack_bit]==-1)
          {
            my_periodic_bit_table[unpack_bit] = pack_bit;
          }
          else
          {
            assert(my_periodic_bit_table[unpack_bit]==pack_bit);
          }
          // find the zone that matches the pack_bits (we need it for the periodic transform)...
          for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
          {
            if (iz->getBits()==pack_bits) break;
          }
          // make sure we found one...
          assert(iz!=faZoneList.end());
          // add an unpack range... 
          prcomm->addUnpackRange(i_f, i_l-1, pack_bits, iz->getKind(), iz->periodic_data);
        }
      }
    }
    if (n_periodic_bit_table>0)
    {
      MPI_Allreduce(my_periodic_bit_table, periodic_bit_table, n_periodic_bit_table, MPI_INT, MPI_MAX, mpi_comm);
      // make sure atleast one person found every entry, and that out found entries match
      // the global found entries...
      for (int i = 0; i<n_periodic_bit_table; i++)
      {
        assert(periodic_bit_table[i]>=0);
        if (my_periodic_bit_table[i]>=0)
        {
          assert(my_periodic_bit_table[i]==periodic_bit_table[i]);
        }
      }
      delete[] my_periodic_bit_table;
    }

  }

  // this has been reused 
  delete[] fa_nbr_rank;

  // now that cvofa[ifa][1]'s have been reset to -1, 
  // build the face of cv structure...
  buildFaocv();

  // do a discrete GCL check...
  // this makes sure that every edge appearing in the cv is
  // matched once and only once in the opposite direction. We
  // can use this check in the future to reduce the node numbering
  // and eliminate common nodes introduced by reconnection. This 
  // will require some sending and recv with the original nodal
  // distribution noora to ensure eliminated nodes are comprehended 
  // on all processors, independent of partition...

  {
    int nodelist1[64];
    int nodelist2[64];
    for (int icv = 0; icv<ncv; icv++)
    {
      int nlist = 0;
      for (int foc = faocv_i[icv]; foc<faocv_i[icv+1]; foc++)
      {
        int ifa = faocv_v[foc];
        if (cvofa[ifa][0]==icv)
        {
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          int ino2 = noofa_v[nof_l];
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino1 = ino2;
            ino2 = noofa_v[nof];
            // this is an outward pointing face, so look for 
            // ino2-ino1 in the edge list...
            int i;
            for (i = 0; i<nlist; i++)
            {
              if ((nodelist1[i]==ino2)&&(nodelist2[i]==ino1)) break;
            }
            if (i<nlist)
            {
              // got a match at i...
              nodelist1[i] = nodelist2[i] = -1;
            }
            else
            {
              // no match was found, so add ino1-ino2...
              assert(nlist<64);
              nodelist1[nlist] = ino1;
              nodelist2[nlist] = ino2;
              nlist++;
            }
          }
        }
        else
        {
          assert(cvofa[ifa][1]==icv);
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          int ino2 = noofa_v[nof_l];
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino1 = ino2;
            ino2 = noofa_v[nof];
            // this is an inward pointing face, so look for 
            // ino1-ino2 in the edge list...
            int i;
            for (i = 0; i<nlist; i++)
            {
              if ((nodelist1[i]==ino1)&&(nodelist2[i]==ino2)) break;
            }
            if (i<nlist)
            {
              // got a match at i...
              nodelist1[i] = nodelist2[i] = -1;
            }
            else
            {
              // no match was found, so add ino2-ino1...
              assert(nlist<64);
              nodelist1[nlist] = ino2;
              nodelist2[nlist] = ino1;
              nlist++;
            }
          }
        }
      }
      // at this point, the node lists should be all negative, indicating that
      // each edge was found twice...
      assert(nlist>=6); // 6 is the smallest number of edges associated with a valid cv (tet)...
      for (int i = 0; i<nlist; i++)
      {
        assert(nodelist1[i]==-1);
        assert(nodelist2[i]==-1);
      }
    }

    // let the user know
    MPI_Barrier(mpi_comm);
    if (mpi_rank==0) cout<<"passed discrete GCL."<<endl;

  }

  // ====================================
  // reorder...
  // ====================================

  {

    // to reorder cells, faces, and (eventually) nodes, we require the following:
    // 1. fa_flag contains the coarse-grain ordering of the faces - i.e. a zero-based
    // index indicating a global face ordering. The ultimate face ordering is then 
    // the best we can do using RCM still obeying these bounds
    //
    // 2. cv_flag contains the coarse-grain ordering of the cells...
    //
    // 3. no_flag contains the coarse-grain ordering of the nodes...
    // 
    // the new ordering is returned in cv_flag, fa_flag, no_flag

    // for faces, right now "fa_zone" contains the index of the
    // FaZone that owns the face. 

    // -------------------
    // faces...
    // faces are ordered according to zone first:
    // 1. boundary zones         [0:nfa_b-1]
    // 2. periodic boundary      [nfa_b:nfa_bp-1]
    // 3. processor boundary     [nfa_bp:nfa_bpi-1]
    // 4. internal               [nfa_bpi:nfa-1]
    // -------------------

    // build a local lookup table to get new zone list order for a given face index...
    int n_fa_order_table = -1;
    list<FaZone>::iterator iz;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
      n_fa_order_table = max(n_fa_order_table, iz->getIndex());
    assert(n_fa_order_table>=0);
    n_fa_order_table++; // zero indexing
    // make sure all fa_flag's fall in this range...
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // recall that face zones can be negative, indicating interprocessor boundary faces...
      int this_fa_zone = max(fa_zone[ifa], -fa_zone[ifa]-1);
      assert((this_fa_zone>=0)&&(this_fa_zone<n_fa_order_table));
    }
    int * fa_order_table = new int[n_fa_order_table];
    for (int i = 0; i<n_fa_order_table; i++)
      fa_order_table[i] = -1;
    int n_zone = 0;
    int n_internal = 0;
    for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
    {
      int index = iz->getIndex();
      assert(fa_order_table[index]==-1);
      fa_order_table[index] = n_zone++;
      // internal zones are treated specially because their faces get split into two 
      // contiguous parts: the processor boundary part (if any), and the fully
      // interior part (if any)...
      if (iz->getKind()==FA_ZONE_INTERNAL) n_internal++;
    }
    // now cycle through all faces and count...
    //cout << "n_internal, n_zone: " << n_internal << " " << n_zone << endl;
    int * fa_count = new int[n_zone+n_internal+1];
    for (int i = 0; i<n_zone+n_internal; i++)
      fa_count[i+1] = 0;
    for (int iter = 0; iter<2; iter++)
    {
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        assert(cvofa[ifa][0]>=0);
        int izone;
        if (fa_zone[ifa]<0)
        {
          int this_fa_zone = -fa_zone[ifa]-1;
          assert(fa_kind_table[this_fa_zone]==FA_ZONE_INTERNAL);
          assert(cvofa[ifa][1]==-1);
          izone = fa_order_table[this_fa_zone];
        }
        else
        {
          int this_fa_zone = fa_zone[ifa];
          if (fa_kind_table[this_fa_zone]==FA_ZONE_INTERNAL)
          {
            assert(cvofa[ifa][1]>=0);
            izone = fa_order_table[this_fa_zone]+n_internal;
          }
          else
          {
            assert(cvofa[ifa][1]==-1);
            izone = fa_order_table[this_fa_zone];
          }
        }
        if (iter==0)
        {
          fa_count[izone+1] += 1;
        }
        else
        {
          fa_flag[ifa] = fa_count[izone];
          fa_count[izone] += 1;
        }
      }
      if (iter==0)
      {
        // build csr...
        fa_count[0] = 0;
        for (int i = 0; i<n_zone+n_internal; i++)
          fa_count[i+1] += fa_count[i];
        // face subdivisions...
        nfa_b = -1;
        nfa_bp = -1;
        nfa_bpi = -1;
        // set faZoneList numbering...
        int izone = 0;
        for (iz = faZoneList.begin(); iz!=faZoneList.end(); iz++)
        {
          switch (iz->getKind())
          {
          case FA_ZONE_BOUNDARY:
            iz->ifa_f = fa_count[izone];
            iz->ifa_l = fa_count[izone+1]-1;
            break;
          case FA_ZONE_PERIODIC_CART:
          case FA_ZONE_PERIODIC_CYL_X:
          case FA_ZONE_PERIODIC_CYL_Y:
          case FA_ZONE_PERIODIC_CYL_Z:
            if (nfa_b==-1) nfa_b = fa_count[izone];
            iz->ifa_f = fa_count[izone];
            iz->ifa_l = fa_count[izone+1]-1;
            break;
          case FA_ZONE_INTERNAL:
            if (nfa_b==-1) nfa_b = fa_count[izone];
            if (nfa_bp==-1) nfa_bp = fa_count[izone];
            if (nfa_bpi==-1) nfa_bpi = fa_count[izone+n_internal];
            iz->ifa_f = fa_count[izone];
            iz->ifa_l = fa_count[izone+1]-1;
            // internal zones have two ranges of faces - the 
            // totally internal range and the interprocessor range...
            iz->ifa_i_f = fa_count[izone+n_internal];
            iz->ifa_i_l = fa_count[izone+n_internal+1]-1;
            break;
          default:
            cerr<<"Error: do not know how to handle this kind: "<<iz->getKind()<<endl;
            throw(-1);
          }
          izone++;
        }
      }
    }
    delete[] fa_order_table;
    delete[] fa_count;

    // fa_flag now holds the new ordering for faces...
    // so modify all face data indexing...

    // faocv_v...
    {
      vector<int> fvec;
      for (int icv = 0; icv<ncv; icv++)
      {
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv+1]-1;
        fvec.resize(foc_l-foc_f+1);
        for (int foc = foc_f; foc<=foc_l; foc++)
        {
          fvec[foc-foc_f] = fa_flag[faocv_v[foc]];
        }
        // sort so the faces are ordered based on new numbering...
        std::sort(fvec.begin(), fvec.end());
        for (int foc = foc_f; foc<=foc_l; foc++)
        {
          faocv_v[foc] = fvec[foc-foc_f];
        }
      }
    }

    // noofa_i/v...
    reorder_csr(noofa_i, noofa_v, fa_flag, nfa);

    // facePrcomm's - sort ranges?...
    for (list<Prcomm>::iterator prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
    {
      for (int i = 0; i<prcomm->packIndexVec.size(); i++)
      {
        int ifa_old = prcomm->packIndexVec[i];
        assert((ifa_old>=0)&&(ifa_old<nfa));
        prcomm->packIndexVec[i] = fa_flag[ifa_old];
      }
      for (int i = 0; i<prcomm->unpackIndexVec.size(); i++)
      {
        int ifa_old = prcomm->unpackIndexVec[i];
        assert((ifa_old>=0)&&(ifa_old<nfa));
        prcomm->unpackIndexVec[i] = fa_flag[ifa_old];
      }
    }

    // cvofa...
    {
      int (*i2)[2] = new int[nfa][2];
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        i2[ifa][0] = cvofa[ifa][0];
        i2[ifa][1] = cvofa[ifa][1];
      }
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        int ifa_new = fa_flag[ifa];
        cvofa[ifa_new][0] = i2[ifa][0];
        cvofa[ifa_new][1] = i2[ifa][1];
      }
      delete[] i2;
    }

    // any registered face data?...
    {
      double * fa_scalar = NULL;
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
      {
        if (data->getDatatype()==FA_DATA)
        {
          if (fa_scalar==NULL) fa_scalar = new double[nfa];
          for (int ifa = 0; ifa<nfa; ifa++)
            fa_scalar[ifa] = (*data->ptr)[ifa];
          for (int ifa = 0; ifa<nfa; ifa++)
            (*data->ptr)[fa_flag[ifa]] = fa_scalar[ifa];
        }
      }
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
      {
        if (data->getDatatype()==FA_DATA)
        {
          if (fa_scalar==NULL) fa_scalar = new double[nfa];
          for (int i = 0; i<3; i++)
          {
            for (int ifa = 0; ifa<nfa; ifa++)
              fa_scalar[ifa] = (*data->ptr)[ifa][i];
            for (int ifa = 0; ifa<nfa; ifa++)
              (*data->ptr)[fa_flag[ifa]][i] = fa_scalar[ifa];
          }
        }
      }
      for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
      {
        if (data->getDatatype()==FA_DATA)
        {
          cerr<<"Error: reordering of registered face tensors not implemented"<<endl;
        }
      }
      if (fa_scalar!=NULL) delete[] fa_scalar;
    }

    // -------------------
    // cvs are ordered as follows:
    // 1. cvs with atleast 1 boundary face and any number of inter-processor faces [0:ncv_b-1]
    //    - this set of cells are listed first to allow certain boundary conditions to
    //      be applied (e.g. gradient correction) because they correspond to the boundary
    //      faces being listed first as well...
    // 2. cvs with no boundary faces, atleast 1 inter-processor face               [ncv_b:ncv_bp-1]
    //    - this set of cvs is broken out to allow hiding of latency. The only cvs with
    //      distance-1 neighbors that could be on other processors (ghost) are [0:ncv_bp-1]
    // 3. cvs with only internal faces...                                          [ncv_bp:ncv-1] 
    //    - this range of cvs has no boundary/inter-processor dependence.
    // -------------------

    // XXX for now, disregard cv zoning...
    for (int icv = 0; icv<ncv; icv++)
      cv_flag[icv] = -1;
    for (int ifa = 0; ifa<nfa_b; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert((icv0>=0)&&(icv0<ncv));
      assert(icv1==-1);
      cv_flag[icv0] = 0;
    }
    for (int ifa = nfa_b; ifa<nfa_bpi; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert((icv0>=0)&&(icv0<ncv));
      assert(icv1==-1);
      if (cv_flag[icv0]==-1) cv_flag[icv0] = 1;
    }
    for (int ifa = nfa_bpi; ifa<nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert((icv0>=0)&&(icv0<ncv));
      assert((icv1>=0)&&(icv1<ncv));
      if (cv_flag[icv0]==-1) cv_flag[icv0] = 2;
      if (cv_flag[icv1]==-1) cv_flag[icv1] = 2;
    }

    {
      int cv_counts[4];
      for (int i = 0; i<3; i++)
        cv_counts[i+1] = 0;
      for (int icv = 0; icv<ncv; icv++)
      {
        // make sure all cv's were counted basedon their face connections...
        assert((cv_flag[icv]>=0)&&(cv_flag[icv]<3));
        cv_counts[cv_flag[icv]+1] += 1;
      }
      // make cv_counts hold the offset...
      cv_counts[0] = 0;
      for (int i = 0; i<3; i++)
        cv_counts[i+1] += cv_counts[i];

      // set the sorted cv ranges...
      ncv_b = cv_counts[1];
      ncv_bp = cv_counts[2];
      assert(ncv==cv_counts[3]);

      for (int icv = 0; icv<ncv; icv++)
      {
        assert((cv_flag[icv]>=0)&&(cv_flag[icv]<3));
        int i = cv_flag[icv];
        cv_flag[icv] = cv_counts[i];
        assert((cv_flag[icv]>=0)&&(cv_flag[icv]<ncv));
        cv_counts[i] += 1;
      }
    }

    // cv_flag now contains a new ordering, so reorder all cv stuff...

    for (int ifa = 0; ifa<nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      cvofa[ifa][0] = cv_flag[icv0];
      if (icv1!=-1)
      {
        assert((icv1>=0)&&(icv1<ncv));
        cvofa[ifa][1] = cv_flag[icv1];
      }
    }

    // cv_zone...
    {
      int *i1 = new int[ncv];
      for (int icv = 0; icv<ncv; icv++)
        i1[icv] = cv_zone[icv];
      for (int icv = 0; icv<ncv; icv++)
        cv_zone[cv_flag[icv]] = i1[icv];
      delete[] i1;
    }

    // faocv csr struct...
    reorder_csr(faocv_i, faocv_v, cv_flag, ncv);

    // registered cv data...
    {
      double * cv_scalar = NULL;
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
      {
        if (data->getDatatype()==CV_DATA)
        {
          if (cv_scalar==NULL) cv_scalar = new double[ncv];
          for (int icv = 0; icv<ncv; icv++)
            cv_scalar[icv] = (*data->ptr)[icv];
          for (int icv = 0; icv<ncv; icv++)
            (*data->ptr)[cv_flag[icv]] = cv_scalar[icv];
        }
      }
      for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
      {
        if (data->getDatatype()==CV_DATA)
        {
          if (cv_scalar==NULL) cv_scalar = new double[ncv];
          for (int i = 0; i<3; i++)
          {
            for (int icv = 0; icv<ncv; icv++)
              cv_scalar[icv] = (*data->ptr)[icv][i];
            for (int icv = 0; icv<ncv; icv++)
              (*data->ptr)[cv_flag[icv]][i] = cv_scalar[icv];
          }
        }
      }
      for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
      {
        if (data->getDatatype()==CV_DATA)
        {
          if (cv_scalar==NULL) cv_scalar = new double[ncv];
          for (int i = 0; i<3; i++)
          {
            for (int j = 0; j<3; j++)
            {
              for (int icv = 0; icv<ncv; icv++)
                cv_scalar[icv] = (*data->ptr)[icv][i][j];
              for (int icv = 0; icv<ncv; icv++)
                (*data->ptr)[cv_flag[icv]][i][j] = cv_scalar[icv];
            }
          }
        }
      }
      if (cv_scalar!=NULL) delete[] cv_scalar;
    }

  }

  delete[] fa_zone;
  delete[] cv_zone;

  // --------------------------------------
  // nodes...
  // --------------------------------------

  //#define NODE_REORDER
#ifdef NODE_REORDER
  {

    // build a temporary nodal connectivity to facilitate this reordering...

    // build cvono_i/v...

    int * cvono_i_tmp = new int[getNno()+1];
    for (int ino = 0; ino < getNno(); ino++)
    cvono_i_tmp[ino+1] = 0;
    int * cvono_v_tmp = NULL;
    for (int iter = 0; iter < 2; iter++)
    {
      // put -1 in all nodes...
      for (int ino = 0; ino < getNno(); ino++)
      no_flag[ino] = -1;
      for (int icv = 0; icv < getNcv(); icv++)
      {
        const int foc_f = faocv_i[icv];
        const int foc_l = faocv_i[icv+1]-1;
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          const int ifa = faocv_v[foc];
          const int nof_f = noofa_i[ifa];
          const int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof <= nof_l; nof++)
          {
            const int ino = noofa_v[nof];
            if (no_flag[ino] != icv)
            {
              no_flag[ino] = icv;
              if (iter == 0)
              {
                cvono_i_tmp[ino+1] += 1;
              }
              else
              {
                cvono_v_tmp[cvono_i_tmp[ino]] = icv;
                cvono_i_tmp[ino] += 1;
              }
            }
          }
        }
      }
      if (iter == 0)
      {
        cvono_i_tmp[0] = 0;
        for (int ino = 0; ino < getNno(); ino++)
        cvono_i_tmp[ino+1] += cvono_i_tmp[ino];
        cvono_v_tmp = new int[cvono_i_tmp[getNno()]];
      }
      else
      {
        for (int ino = getNno(); ino > 0; ino--)
        cvono_i_tmp[ino] = cvono_i_tmp[ino-1];
        cvono_i_tmp[0] = 0;
      }
    }

    int * nbono_i_tmp = new int[getNno()+1];
    for (int ino = 0; ino < getNno(); ino++)
    nbono_i_tmp[ino+1] = 0;
    int * nbono_v_tmp = NULL;
    for (int iter = 0; iter < 2; iter++)
    {
      // put -1 in all nodes...
      for (int ino = 0; ino < getNno(); ino++)
      no_flag[ino] = -1;
      for (int ino = 0; ino < getNno(); ino++)
      {
        // no need to include diagonal...
        assert( no_flag[ino] != ino );
        no_flag[ino] = ino;
        // get nbrs from cvono structure...
        int con_f = cvono_i_tmp[ino];
        int con_l = cvono_i_tmp[ino+1]-1;
        for (int con = con_f; con <= con_l; con++)
        {
          const int icv = cvono_v_tmp[con];
          const int foc_f = faocv_i[icv];
          const int foc_l = faocv_i[icv+1]-1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            const int ifa = faocv_v[foc];
            const int nof_f = noofa_i[ifa];
            const int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              const int ino_nbr = noofa_v[nof];
              if (no_flag[ino_nbr] != ino)
              {
                no_flag[ino_nbr] = ino;
                if (iter == 0)
                {
                  nbono_i_tmp[ino+1] += 1;
                }
                else
                {
                  nbono_v_tmp[nbono_i_tmp[ino]] = ino_nbr;
                  nbono_i_tmp[ino] += 1;
                }
              }
            }
          }
        }
      }
      if (iter == 0)
      {
        nbono_i_tmp[0] = 0;
        for (int ino = 0; ino < getNno(); ino++)
        nbono_i_tmp[ino+1] += nbono_i_tmp[ino];
        nbono_v_tmp = new int[nbono_i_tmp[getNno()]];
      }
      else
      {
        for (int ino = getNno(); ino > 0; ino--)
        nbono_i_tmp[ino] = nbono_i_tmp[ino-1];
        nbono_i_tmp[0] = 0;
      }
    }

    delete[] cvono_i_tmp;
    delete[] cvono_v_tmp;

    // now use no_flag to store a "levelset" from the nodal data that must
    // be packed and exchanged...

    for (int ino = 0; ino < nno; ino++)
    no_flag[ino] = -1;

    // now put a 0 in every node that requires communication...
    for (list<Prcomm>::iterator prcomm = facePrcommList.begin(); prcomm != facePrcommList.end(); prcomm++)
    {
      for (int i = 0; i < prcomm->packIndexVec.size(); i++)
      {
        const int ifa = prcomm->packIndexVec[i];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof <= nof_l; nof++)
        {
          const int ino = noofa_v[nof];
          no_flag[ino] = 0;
        }
      }
    }

    // now advance the "front" through the remaining nodes to produce a positive
    // int in every no_flag where the processor boundary is all 0, first row
    // of nodes in is 1, etc...

    int done = 0;
    while (done != 1)
    {
      // assume we are done...
      done = 1;
      for (int ino = 0; ino < nno; ino++)
      {
        // for each node, find the minimum no_flag amongst the
        // nbrs...
        int no_flag_min = -1;
        const int non_f = nbono_i_tmp[ino];
        const int non_l = nbono_i_tmp[ino+1]-1;
        // recal tmp has no diagonal...
        for (int non = non_f; non <= non_l; non++)
        {
          const int ino_nbr = nbono_v_tmp[non];
          if (no_flag[ino_nbr] >= 0)
          {
            if (no_flag_min == -1)
            no_flag_min = no_flag[ino_nbr];
            else
            no_flag_min = min( no_flag_min, no_flag[ino_nbr] );
          }
        }
        if ( no_flag_min >= 0 )
        {
          if ( no_flag[ino] == -1 )
          {
            done = 0;
            no_flag[ino] = no_flag_min+1;
          }
          else if ( no_flag[ino] > no_flag_min+1 )
          {
            done = 0;
            no_flag[ino] = no_flag_min+1;
          }
        }
      }
    }

    // now count how many members are in each level...
    // note that it is still possible that some/all nodes
    // are -1, so set these to 1...

    int nlevels = 0;
    for (int ino = 0; ino < nno; ino++)
    {
      if (no_flag[ino] == -1)
      no_flag[ino] = 1;
      nlevels = max( nlevels, no_flag[ino] );
    }
    nlevels += 1;

    int * level = new int[nlevels+1];
    for (int i = 0; i < nlevels; i++)
    level[i+1] = 0;

    // we actually want to pack the max level first, so flip the
    // numbering, so the exchanged boundary is last. It could go either
    // way, but this is more consistent with active node numbering, where
    // additional, neccessary nodes will be available locally, but
    // do not need to be iterated on (see UgpWithToolsNo.h/cpp)...

    for (int ino = 0; ino < nno; ino++)
    {
      const int ilevelp1 = nlevels - no_flag[ino]; // 1..nlevels
      assert( ilevelp1 >= 1 );
      no_flag[ino] = -ilevelp1;
      level[ilevelp1] += 1;
    }

    level[0] = 0;
    for (int i = 0; i < nlevels; i++)
    level[i+1] += level[i];
    assert( level[nlevels] == nno );

    // we can now set nno_i, the internal node range...
    nno_i = level[nlevels-1];

    // we need a stack...
    int * stack = new int[nno];
    for (int i = 0; i < nno; i++)
    stack[i] = -1;

    // now loop through the stack...
    for (int i = 0; i < nno; i++)
    {
      if (stack[i] == -1)
      {
        // we need to find the next available node for this
        // position in the stack...
        int ino_seed = -1;
        int nun_seed = 0; // nun = number of ungrouped nbrs, gets reset below
        for (int ino = 0; ino < nno; ino++)
        {
          if (no_flag[ino] < 0)
          {
            // compute the complexity...
            int nun = 0;
            const int non_f = nbono_i_tmp[ino];
            const int non_l = nbono_i_tmp[ino+1]-1;
            for (int non = non_f; non <= non_l; non++)
            {
              const int ino_nbr = nbono_v_tmp[non];
              if( no_flag[ino_nbr] < 0 )
              nun++;
            }
            // this is a valid candidate...
            if ( (ino_seed == -1) ||
                (no_flag[ino] > no_flag[ino_seed]) ||
                ((no_flag[ino] == no_flag[ino_seed])&&(nun < nun_seed)) )
            {
              ino_seed = ino;
              nun_seed = nun;
            }
          }
        }
        // we should have found one...
        assert( ino_seed >= 0 );
        // the level from the seed should match the stack index...
        const int ilevel = -no_flag[ino_seed]-1;
        assert( level[ilevel] == i );
        stack[i] = ino_seed;
        no_flag[ino_seed] = i;
        level[ilevel] += 1;
      }

      // pop this node off the stack...
      const int ino = stack[i];
      assert( no_flag[ino] == i );

      // now cycle through its nbrs and push into the stack...
      const int non_f = nbono_i_tmp[ino];
      const int non_l = nbono_i_tmp[ino+1]-1;
      // note that the diagonal is not in these tmp connectivities...
      for (int non = non_f; non <= non_l; non++)
      {
        int ino_nbr = nbono_v_tmp[non];
        if (no_flag[ino_nbr] < 0)
        {
          const int ilevel = -no_flag[ino_nbr]-1;
          const int j = level[ilevel];
          level[ilevel] += 1;
          assert( stack[j] == -1 );
          stack[j] = ino_nbr;
          no_flag[ino_nbr] = j;
        }
      }

    }

    // cleanup...
    delete[] level;
    delete[] stack;

    /*

     // =====================================
     // report matrix bandwidth...
     // =====================================

     // old/current...

     int max_bw = 0;
     double avg_bw = 0.0;
     for (int ino = 0; ino < nno; ino++) {
     int ino_min = ino;
     int ino_max = ino;
     const int non_f = nbono_i_tmp[ino];
     const int non_l = nbono_i_tmp[ino+1]-1;
     for (int non = non_f; non <= non_l; non++) {
     const int ino_nbr = nbono_v_tmp[non];
     ino_min = min( ino_min, ino_nbr );
     ino_max = max( ino_max, ino_nbr );
     }
     max_bw = max( max_bw, ino_max - ino_min + 1 );
     avg_bw += (double)(ino_max - ino_min + 1);
     }
     cout << "old matrix bw: max: " << max_bw << ", avg: " << avg_bw/(double)nno << endl;

     // new...

     max_bw = 0;
     avg_bw = 0.0;
     for (int ino = 0; ino < nno; ino++) {
     int ino_min = no_flag[ino];
     int ino_max = no_flag[ino];
     const int non_f = nbono_i_tmp[ino];
     const int non_l = nbono_i_tmp[ino+1]-1;
     for (int non = non_f; non <= non_l; non++) {
     const int ino_nbr = no_flag[nbono_v_tmp[non]];
     ino_min = min( ino_min, ino_nbr );
     ino_max = max( ino_max, ino_nbr );
     }
     max_bw = max( max_bw, ino_max - ino_min + 1 );
     avg_bw += (double)(ino_max - ino_min + 1);
     }
     cout << "new matrix bw: max: " << max_bw << ", avg: " << avg_bw/(double)nno << endl;

     */

    delete[] nbono_i_tmp;
    delete[] nbono_v_tmp;

    // ---------------------------------------------------------------------------
    // now no_flag contains the new node ordering, so reorder all things nodal...
    // ---------------------------------------------------------------------------

    // noofa_v...
    for (int nof = 0; nof < noofa_s; nof++)
    noofa_v[nof] = no_flag[noofa_v[nof]];

    // noocv_v - does not exist here yet...
    /*
     for (int noc = 0; noc < noocv_s; noc++)
     noocv_v[noc] = no_flag[noocv_v[noc]];
     */

    // node prcomm - does not exist yet...
    /*
     for (list<Prcomm>::iterator prcomm = nodePrcommList.begin(); prcomm != nodePrcommList.end(); prcomm++) {
     for (int i = 0; i < prcomm->packIndexVec.size(); i++)
     prcomm->packIndexVec[i] = no_flag[prcomm->packIndexVec[i]];
     for (int i = 0; i < prcomm->unpackIndexVec.size(); i++)
     prcomm->unpackIndexVec[i] = no_flag[prcomm->unpackIndexVec[i]];
     }
     */

    // we need a tmp vector...
    double * tmp = new double[nno];

    // x_no coords...

    for (int i = 0; i < 3; i++)
    {
      for (int ino = 0; ino < nno; ino++)
      tmp[ino] = x_no[ino][i];
      for (int ino = 0; ino < nno; ino++)
      x_no[no_flag[ino]][i] = tmp[ino];
    }

    // registered node data...

    // double scalars...

    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++)
    {
      if (data->getDatatype() == NO_DATA)
      {
        for (int ino = 0; ino < nno; ino++)
        tmp[ino] = (*(data->ptr))[ino];
        for (int ino = 0; ino < nno; ino++)
        (*(data->ptr))[no_flag[ino]] = tmp[ino];
      }
    }

    // double vectors: [nno][3]...

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++)
    {
      if (data->getDatatype() == NO_DATA)
      {
        for (int i = 0; i < 3; i++)
        {
          for (int ino = 0; ino < nno; ino++)
          tmp[ino] = (*(data->ptr))[ino][i];
          for (int ino = 0; ino < nno; ino++)
          (*(data->ptr))[no_flag[ino]][i] = tmp[ino];
        }
      }
    }

    // double tensors: [nno][3][3]...

    for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data != doubleTensorList.end(); data++)
    {
      if (data->getDatatype() == NO_DATA)
      {
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            for (int ino = 0; ino < nno; ino++)
            tmp[ino] = (*(data->ptr))[ino][i][j];
            for (int ino = 0; ino < nno; ino++)
            (*(data->ptr))[no_flag[ino]][i][j] = tmp[ino];
          }
        }
      }
    }

    delete[] tmp;

    //MPI_Pause("done RCM");

  }
#endif

  // finish off...

  checkFaceCorrespondence();

  syncNoofa();

  reconnectNodes();

  // synchronize nodes - only if the user demands it...

  // sync

#ifdef NODE_REORDER
  // check node ordering...
  for (int ino = 0; ino < nno; ino++)
  no_flag[ino] = 1;
  updateNoI1(no_flag,ADD_DATA);
  for (int ino = 0; ino < nno_i; ino++)
  assert( no_flag[ino] == 1 );
  for (int ino = nno_i; ino < nno; ino++)
  assert( no_flag[ino] > 1 );
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  cout << " > node ordering (nno_i) OK." << endl;
#endif

  checkNodeCorrespondence();

  //checkCvonoCounts();

  buildNoocv();
  reorderFaocvNoocv();

} // Ugp::redistReconnectReorder()

// ###########################################################

void Ugp::syncNoofa()
{

  if (mpi_rank==0) cout<<"syncNoofa()"<<endl;

  // step 1 - tie break. Only one of the 2 faces is going to
  // rotate its nodal ordering...

  // start by checking that the noofa count is the same across processors...
  int my_nnof_max = 0;
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    int nnof = noofa_i[ifa+1]-noofa_i[ifa];
    my_nnof_max = max(my_nnof_max, nnof);
    assert((nnof>=3)&&(nnof<=8));
    fa_flag[ifa] = nnof;
  }
  updateFaI1(fa_flag, REPLACE_DATA);
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    if (fa_flag[ifa]!=noofa_i[ifa+1]-noofa_i[ifa])
    {
      cerr<<"Error: node of face count does not match"<<endl;
      throw(-1);
    }
  }

  // now tie break to select one of the two faces to potentially
  // do node cycling...
  int my_offset;
  MPI_Scan(&nfa, &my_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  my_offset -= nfa;
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_flag[ifa] = ifa+my_offset;
  updateFaI1(fa_flag, MIN_DATA);
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    if (fa_flag[ifa]==ifa+my_offset)
    {
      fa_flag[ifa] = 0;
    }
    else
    {
      fa_flag[ifa] = 1;
    }
  }

  // exchange the coordinates of the first node...
  double (*x_fa)[3] = new double[nfa][3];

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    int ino = noofa_v[noofa_i[ifa]];
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] = x_no[ino][i];
  }
  updateFaR2(x_fa, REPLACE_TRANSLATE_DATA);
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    // only check matching node in higher face as determined by
    // the tie breaking algorithm above...
    if (fa_flag[ifa]==1)
    {
      // select the first node...
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1]-1;
      // select the first...
      int nof_min = nof_f;
      double d2_min = 0.0;
      for (int i = 0; i<3; i++)
      {
        double dx = x_fa[ifa][i]-x_no[noofa_v[nof_min]][i];
        d2_min += dx*dx;
      }
      // check against the rest...
      for (int nof = nof_f+1; nof<=nof_l; nof++)
      {
        double d2 = 0.0;
        for (int i = 0; i<3; i++)
        {
          double dx = x_fa[ifa][i]-x_no[noofa_v[nof]][i];
          d2 += dx*dx;
        }
        if (d2<d2_min)
        {
          nof_min = nof;
          d2_min = d2;
        }
      }
      // check that the selected node is atleast X (twice) as far from all
      // other nodes than from the selected node...
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        if (nof!=nof_min)
        {
          double d2 = 0.0;
          for (int i = 0; i<3; i++)
          {
            double dx = x_fa[ifa][i]-x_no[noofa_v[nof]][i];
            d2 += dx*dx;
          }
          if (d2<=d2_min)
          {
            cout<<"d2, d2_min: "<<d2<<" "<<d2_min<<endl;
            cout<<"nof_f, nof_l : "<<nof_f<<" "<<nof_l<<endl;
            for (int nof2 = nof_f; nof2<=nof_l; nof2++)
            {
              cout<<"nof2, ino, x: "<<nof2<<" "<<noofa_v[nof2]<<" "<<x_no[noofa_v[nof2]][0]<<" "
                  <<x_no[noofa_v[nof2]][1]<<" "<<x_no[noofa_v[nof2]][2]<<endl;
            }
          }
          assert(d2>d2_min);
          if (d2<2.0*d2_min)
          {
            cerr<<"Error: syncNoofa nodes are too close compared to tol: "<<sqrt(d2_min)<<" "<<sqrt(d2)<<endl;
            throw(-1);
          }
        }
      }
      // we want to rotate the node list until the selected
      // node at nof_min is the last node...
      int ino_l = noofa_v[nof_min];
      while (noofa_v[nof_l]!=ino_l)
      {
        //cout << "rotating nodes" << endl;
        int ino_f = noofa_v[nof_f];
        for (int nof = nof_f; nof<nof_l; nof++)
          noofa_v[nof] = noofa_v[nof+1];
        noofa_v[nof_l] = ino_f;
      }
    }
  }

  // ---------------------------------------------
  // check global correspondence of all nodes...
  // ---------------------------------------------

  int nnof_max;
  MPI_Allreduce(&my_nnof_max, &nnof_max, 1, MPI_INT, MPI_MAX, mpi_comm);

  // put a 2 in boundary faces...
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_flag[ifa] = 1;
  updateFaI1(fa_flag, ADD_DATA);

  double my_d2_max = 0.0;
  for (int nof_inc = 0; nof_inc<nnof_max; nof_inc++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_flag[ifa]==2)
      {
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nof = nof_f+nof_inc;
        if (nof<=nof_l)
        {
          int ino = noofa_v[nof];
          for (int i = 0; i<3; i++)
            x_fa[ifa][i] = x_no[ino][i];
        }
      }
    }
    updateFaR2(x_fa, REPLACE_TRANSLATE_DATA);
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_flag[ifa]==2)
      {
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nof = nof_l-nof_inc;
        if (nof>=nof_f)
        {
          int ino = noofa_v[nof];
          double d2 = 0.0;
          for (int i = 0; i<3; i++)
          {
            double dx = x_fa[ifa][i]-x_no[ino][i];
            d2 += dx*dx;
          }
          my_d2_max = max(my_d2_max, d2);
        }
      }
    }
  }

  delete[] x_fa;

  double d2_max;
  MPI_Reduce(&my_d2_max, &d2_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank==0) cout<<"syncNoofa: nodes of faces matched wihin tol (should be small): "<<sqrt(d2_max)<<endl;

}

void Ugp::reconnectNodes()
{

  // requires noofa to be synced across processors

  // nodal reconnection builds the nodal paired communicators
  // requires the face paired communicators...

  if (mpi_rank==0) cout<<"reconnectNodes()"<<endl;

  assert(nodePrcommList.size()==0);

  // every node needs to know what other nodes it is associated
  // with across inter-processor internal, periodic, and full 
  // boundaries. This is done through the already-existing 
  // face prcomm lists...

  // we want to store data associated with each node:
  // 1. transform index (1,2,3,etc and bitwise combinations?)
  // 2. associated rank for nbr - can be mpi_rank (self-periodic or self-full)
  // 3. local node number on nbr rank - cannot be ino

  // this is done with the following algorithm: first the node exchanges are seeded with
  // whatever our own face transformation is, and this data is stored on the
  // receiving proces in a data-of-node structure. On subsequent packs, the
  // data-of-node structure is cycled through and new entries are transformed and 
  // exchanged. Once an entry is transformed and exchanged, it is marked as done
  // by setting its local node number ino to -1


  // start by ensuring the face nodes are sync'd across processors - this means
  // that our first node it the recieving processors last node, etc...

  int * daono_i = new int[nno+1];
  int ino;
  for (ino = 0; ino<=nno; ino++)
    daono_i[ino] = 0;

  int * new_data_count = new int[nno];
  vector<Daono> daono_v;

  int first = 1;
  int done = 0;
  while (done!=1)
  {

    // ------------------------------------------
    // pack up the data...
    // ------------------------------------------

    for (int iter = 0; iter<2; iter++)
    {
      //cout << "pack iter: " << iter << endl;
      for (list<Prcomm>::iterator prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
      {
        prcomm->npack_v = 0;
        for (ino = 0; ino<nno; ino++)
          no_flag[ino] = -1;
        for (list<Range>::iterator range = prcomm->packRangeList.begin(); range!=prcomm->packRangeList.end(); range++)
        {
          int i_f = range->getIndexFirst();
          int i_l = range->getIndexLast();
          int periodic_bits = range->getBits(); // may be zero (INTERNAL) or a single specific bit indicating a FULL or PERIODIC range 
          int inverse_periodic_bits = getInversePeriodicBits(periodic_bits);
          for (int i = i_f; i<=i_l; i++)
          {
            int ifa = prcomm->packIndexVec[i];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof<=nof_l; nof++)
            {
              ino = noofa_v[nof];
              // use no_flag to prevent double accounting in a given range...
              if (no_flag[ino]!=periodic_bits)
              {
                no_flag[ino] = periodic_bits;
                if (first==1)
                {
                  // on the first time, simply pack our own node, bits, and rank
                  // to seed the exchanges...
                  if (iter==1)
                  {
                    prcomm->packBufferInt[prcomm->npack_v] = ino;
                    prcomm->packBufferInt[prcomm->npack_v+1] = periodic_bits;
                    prcomm->packBufferInt[prcomm->npack_v+2] = mpi_rank;
                  }
                  prcomm->npack_v += 3;
                }
                else
                {
                  // on subsequent passes, repackage and send back all data that is new
                  // as indicated by ino == -1...
                  int don_f = daono_i[ino];
                  int don_l = daono_i[ino+1]-1;
                  for (int don = don_f; don<=don_l; don++)
                  {
                    if (daono_v[don].ino==-1)
                    {
                      // this is a new member, so pack if our current transformation doesn't simply
                      // reverse the existing transformation...
                      if ((periodic_bits==0)||((daono_v[don].bits&inverse_periodic_bits)==0))
                      {
                        if (iter==1)
                        {
                          prcomm->packBufferInt[prcomm->npack_v] = daono_v[don].ino_nbr;
                          prcomm->packBufferInt[prcomm->npack_v+1] = daono_v[don].bits|periodic_bits;
                          prcomm->packBufferInt[prcomm->npack_v+2] = daono_v[don].rank_nbr;
                        }
                        prcomm->npack_v += 3;
                      }
                    }
                  }
                }
                // always include a -1 to indicate we should move on to the next node...
                if (iter==1) prcomm->packBufferInt[prcomm->npack_v] = -1;
                prcomm->npack_v += 1;
              }
            }
          }
        }
        if (iter==0)
        {
          // pack buffer needs space for npack...
          prcomm->ensurePackBufferIntSize(prcomm->npack_v);
        }
      }
    }

    // ------------------------------------------
    // as soon as we have finished packing, mark all current
    // data-of-node to indicate they have already been 
    // combined with our local transformations - this is done
    // by setting daono_v[].ino to the local ino (initially -1)...
    // ------------------------------------------

    for (ino = 0; ino<nno; ino++)
    {
      int don_f = daono_i[ino];
      int don_l = daono_i[ino+1]-1;
      int don;
      for (don = don_f; don<=don_l; don++)
      {
        if (daono_v[don].ino==-1)
        {
          daono_v[don].ino = ino;
        }
        else
        {
          // check that ino is properly set...
          assert(daono_v[don].ino==ino);
        }
      }
    }

    // ------------------------------------------
    // do a variable sized buffer exchange...
    // note that this sets prcomm%nunpack_v and ensures 
    // the unpackBufferInt size, then exchanges current
    // packBufferInt into unpackBufferInt...
    // ------------------------------------------

    exchangePrcommBufferIntV(facePrcommList);

    // ------------------------------------------
    // unpack...
    // ------------------------------------------

    for (ino = 0; ino<nno; ino++)
      new_data_count[ino] = 0;

    // assume we are done...
    int my_done = 1;

    for (int iter = 0; iter<2; iter++)
    {
      //cout << "unpack iter: " << iter << endl;
      for (list<Prcomm>::iterator prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
      {
        int nunpack_v = 0;
        int rank = prcomm->getNbrRank();
        for (ino = 0; ino<nno; ino++)
          no_flag[ino] = -1;
        for (list<Range>::iterator range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
        {
          int i_f = range->getIndexFirst();
          int i_l = range->getIndexLast();
          int periodic_bits = range->getBits();
          for (int i = i_f; i<=i_l; i++)
          {
            int ifa = prcomm->unpackIndexVec[i];
            // Note: here we assume the pack and unpack noofa ordering is flipped
            // across processors/boundaries - if you called syncNoofa, this should be the case...
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_l; nof>=nof_f; nof--)
            {
              ino = noofa_v[nof];
              // use no_flag to prevent double accounting...
              if (no_flag[ino]!=periodic_bits)
              {
                no_flag[ino] = periodic_bits;
                int ino_nbr = prcomm->unpackBufferInt[nunpack_v++];
                while (ino_nbr!=-1)
                {
                  int bits_nbr = prcomm->unpackBufferInt[nunpack_v++];
                  int rank_nbr = prcomm->unpackBufferInt[nunpack_v++];
                  // do something here with ino_nbr/bits_nbr/rank_nbr...
                  int found = 0;
                  int don = -1;
                  if ((ino_nbr==ino)&&(rank_nbr==mpi_rank))
                  {
                    // if we get ourself back, bits_nbr should be zero because
                    // we purposely avoid reversing transformations, however multiple
                    // non-transformations have to be allowed and cannot easily be
                    // detected, thus this little bit of extra messaging...
                    assert(bits_nbr==0);
                    found = 1;
                  }
                  else
                  {
                    // cycle through the data, 
                    // looking for matches...
                    int don_f = daono_i[ino];
                    int don_l = daono_i[ino+1]-1;
                    for (don = don_f; don<=don_l; don++)
                    {
                      // check...
                      if (daono_v[don].ino==-2)
                      {
                        // this is a free position in the vector, so we did not match. 
                        // anything in the current daono_v and we should insert here, 
                        // so break out with the don value...
                        assert(iter==1);
                        break;
                      }
                      else if ((daono_v[don].ino_nbr==ino_nbr)&&(daono_v[don].rank_nbr==rank_nbr))
                      {
                        // bits should match...
                        assert(daono_v[don].bits==bits_nbr);
                        found = 1;
                        break;
                      }
                    }
                  }
                  if (found==0)
                  {
                    if (iter==0)
                    {
                      // clear done...
                      my_done = 0;
                      // on the first iter, we should have searched the whole don range...
                      assert(don==daono_i[ino+1]);
                      new_data_count[ino] += 1;
                    }
                    else
                    {
                      // we did not find and need to add - this should have the
                      // don value in it where we need to add...
                      assert((don>=daono_i[ino])&&(don<=daono_i[ino+1]-1));
                      assert(daono_v[don].ino==-2);
                      // set the data...
                      daono_v[don].ino_nbr = ino_nbr;
                      daono_v[don].bits = bits_nbr;
                      daono_v[don].rank_nbr = rank_nbr;
                      daono_v[don].ino = -1; // switch to -1 to indicate set
                    }
                  }

                  // try for next ino_nbr...
                  ino_nbr = prcomm->unpackBufferInt[nunpack_v++];

                }
              }
            }
          }
        }
        assert(nunpack_v==prcomm->nunpack_v);
      }

      // if we are done, then there is no need to do second iter...
      if (my_done==1)
      {
        break;
      }
      else if (iter==0)
      {

        // make space in daono_i and daono_v[] vectors for 
        // new ino/bits/rank combos...

        int ishift = 0;
        for (ino = 0; ino<nno; ino++)
          ishift += new_data_count[ino];
        assert(ishift>0);

        // step 1 - resize...
        daono_v.resize(daono_v.size()+ishift);
        daono_i[nno] += ishift;
        for (ino = nno-1; ino>=0; ino--)
        {

          ishift -= new_data_count[ino];
          daono_i[ino] += ishift;
          int don_f = daono_i[ino];
          int don_l = daono_i[ino+1]-1;

          // clear the unset values...
          int don_l0 = don_l-new_data_count[ino];

          // clear the locations for the new data... 
          int don;
          for (don = don_l0+1; don<=don_l; don++)
          {
            daono_v[don].ino_nbr = -1;
            daono_v[don].bits = -1;
            daono_v[don].rank_nbr = -1;
            daono_v[don].ino = -2;
          }

          // do this backwards, to avoid overwriting...
          for (don = don_l0; don>=don_f; don--)
          {
            daono_v[don].ino_nbr = daono_v[don-ishift].ino_nbr;
            daono_v[don].bits = daono_v[don-ishift].bits;
            daono_v[don].rank_nbr = daono_v[don-ishift].rank_nbr;
            // use ino entry for some checking...
            assert(daono_v[don-ishift].ino==ino);
            daono_v[don-ishift].ino = -3; // should get overwritten 
            daono_v[don].ino = ino;
          }
        }

        // check...
        assert(ishift==0);

      }
      else
      {

        // on the second iter, there may have been some duplicate
        // entries in the new data. These will be indicated by their
        // daono_v[don].ino == -2:

        int ishift = 0;
        for (ino = 0; ino<nno; ino++)
        {
          int don_f = daono_i[ino];
          if (ishift>0) daono_i[ino] -= ishift;
          int don_l = daono_i[ino+1]-1;
          int don;
          for (don = don_f; don<=don_l; don++)
          {
            if (daono_v[don].ino==-2)
            {
              ishift++;
            }
            else
            {
              assert((daono_v[don].ino==ino)||(daono_v[don].ino==-1));
              if (ishift>0)
              {
                daono_v[don-ishift].ino_nbr = daono_v[don].ino_nbr;
                daono_v[don-ishift].bits = daono_v[don].bits;
                daono_v[don-ishift].rank_nbr = daono_v[don].rank_nbr;
                daono_v[don-ishift].ino = daono_v[don].ino;
              }
            }
          }
        }
        if (ishift>0)
        {
          assert(daono_v.size()==daono_i[nno]);
          daono_i[nno] -= ishift;
          daono_v.resize(daono_v.size()-ishift);
        }

      }

    }

    /*
     cout << "******* node data *********" << endl;
     for (ino = 0; ino < nno; ino++) {
     int don_f = daono_i[ino];
     int don_l = daono_i[ino+1]-1;
     cout << "don entries for ino: " << ino << endl;
     int don;
     for (don = don_f; don <= don_l; don++) {
     cout << " > ino_nbr, rank_nbr, bits, ino : " <<
     daono_v[don].ino_nbr << " " <<
     daono_v[don].rank_nbr << " " <<
     daono_v[don].bits << " " <<
     daono_v[don].ino << endl;
     }
     }
     cout << "******* node data *********" << endl;
     */

    // clear first...
    first = 0;

    // if anyone's not done, then non of us are done...
    MPI_Allreduce(&my_done, &done, 1, MPI_INT, MPI_MIN, mpi_comm);

  }

  /*
   cout << "******* final node data *********" << endl;
   for (ino = 0; ino < nno; ino++) {
   int don_f = daono_i[ino];
   int don_l = daono_i[ino+1]-1;
   cout << "don entries for ino: " << ino << endl;
   int don;
   for (don = don_f; don <= don_l; don++) {
   cout << " > ino_nbr, rank_nbr, bits, ino : " <<
   daono_v[don].ino_nbr << " " <<
   daono_v[don].rank_nbr << " " <<
   daono_v[don].bits << " " <<
   daono_v[don].ino << endl;
   }
   }
   cout << "******* final node data *********" << endl;
   */

  // preliminary cleanup...
  delete[] daono_i;
  delete[] new_data_count;

  // ------------------------------------------------
  // now sort the daono_v data according to 
  // 1. rank
  // 2. bits
  // 3. local node number
  // this will make packing contiguous.
  // ------------------------------------------------

  // sorting the daono_v vector as is will allow us to build the
  // unpack data sizes and ranges, but we cannot get the unpack
  // indices until after we exchange with the pack (see further below)...

  if (daono_v.size()>0) std::sort(daono_v.begin(), daono_v.end(), DaonoCompare());

  // take a look...
  /*
   {
   cout << "******* node data *********" << endl;
   for (int don = 0; don < daono_v.size(); don++) {
   cout << mpi_rank << " > ino, ino_nbr, rank_nbr, bits: " <<
   daono_v[don].ino << " " <<
   daono_v[don].ino_nbr << " " <<
   daono_v[don].rank_nbr << " ";
   for (int ii = 0; ii < 24; ii++) {
   if (daono_v[don].bits & (1<<ii)) {
   cout << "1";
   }
   else {
   cout << "0";
   }
   }
   cout << "\n";
   }
   cout << "******* node data *********" << endl;
   }
   */
  //MPI_Pause("unpack sort");

  // ------------------------------
  // unpack Prcomms and Ranges...
  // ------------------------------

  buildPrcommsAndRanges(daono_v, nodePrcommList);

  // daono_v cleaned automatically

}

void Ugp::buildPrcommsAndRanges(vector<Daono>& daono_v, list<Prcomm>& prcommList)
{

  assert(prcommList.size()==0);

  int don_l = 0;
  while (don_l<daono_v.size())
  {
    int don_f = don_l;
    int rank_nbr = daono_v[don_f].rank_nbr;
    // advance don_l to the end of the rank...
    while ((don_l<daono_v.size())&&(daono_v[don_l].rank_nbr==rank_nbr))
      don_l++;

    // now [don_f:don_l-1] contains the same rank, rank_nbr...
    // these are all the nodes going to rank_nbr.
    //cout << " > nodes to rank: " << rank_nbr << " range, count: " << 
    //don_f << ":" << don_l-1 << " " << don_l-don_f << endl;

    // add the new paired communicator and resize the unpackIndexVec...
    prcommList.push_back(Prcomm(rank_nbr));
    Prcomm * prcomm = &(prcommList.back());
    prcomm->unpackIndexVec.resize(don_l-don_f);
    prcomm->ensureUnpackBufferIntSize(2*(don_l-don_f));

    // the pack buffer will get the same size - i.e. the number of nodes 
    // we recv is the same as the number we pack and send...
    //prcomm->packIndexVec.resize(don_l-don_f); --- NOT ALWAYS.

    // put our localnode numbers in the unpack vec...
    // and we will put 2 ints into the unpackBufferInt:
    // first is the unpack node on the other side, second
    // is the periodic bits, so the other side knows how to pack.
    for (int don = don_f; don<don_l; don++)
    {
      prcomm->unpackIndexVec[don-don_f] = daono_v[don].ino;
      prcomm->unpackBufferInt[2*(don-don_f)] = daono_v[don].ino_nbr;
      prcomm->unpackBufferInt[2*(don-don_f)+1] = daono_v[don].bits;
    }

    // now add ranges...
    int don2_l = don_f;
    while (don2_l<don_l)
    {
      int don2_f = don2_l;
      int bits = daono_v[don2_f].bits;
      // advance don2_l to the end of the bits...
      while ((don2_l<don_l)&&(daono_v[don2_l].bits==bits))
        don2_l++;
      // add the range...
      if (bits==0)
      {
        assert(don2_l-don_f-1<don_l-don_f);
        prcomm->addUnpackRange(don2_f-don_f, don2_l-don_f-1, 0, FA_ZONE_INTERNAL, NULL);
      }
      else
      {

        // the bits contain how we got here, so the inverse is how to get back...
        int inverse_bits = getInversePeriodicBits(bits);

        // also add possibility to double a particular transform: this is required
        // for the case of neighbors of neighbors...
        int double_transform = 0;
        if (bits&(1<<30))
        {
          double_transform = 1;
          inverse_bits |= (1<<30);
        }

        // Multiple periodicity requires the combined transformation. For now, this requires
        // that the transformation be either a single transformation or a sum of
        // Cartesian transformations. In the future, we could support Cartesian+
        // cylindrical (e.g. periodic annular region), but this would require 
        // that the transforms commute (i.e. order of application does not matter), 
        // and the flag and virtual function machinery for transfromation and 
        // rotation can handle this - probably need to increase the amount of data a 
        // periodic Face zone and Range hold (currently 3 doubles).

        int zone_kind = -1;
        double periodic_data[3];
        for (int bit = 0; bit<n_periodic_bit_table; bit++)
        {
          if (inverse_bits&(1<<bit))
          {
            // get this zone from the faZoneList...
            list<FaZone>::iterator zone;
            for (zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
            {
              if (zone->getBits()==(1<<bit)) break;
            }
            // make sure we found one...
            assert(zone!=faZoneList.end());
            if (zone_kind==-1)
            {
              // this is the first we found - allow any supported periodicity...
              zone_kind = zone->getKind();
              for (int i = 0; i<3; i++)
                periodic_data[i] = zone->periodic_data[i];
              // check for doubling here...
              if (double_transform==1)
              {
                switch (zone_kind)
                {
                case FA_ZONE_PERIODIC_CART:
                  for (int i = 0; i<3; i++)
                    periodic_data[i] *= 2.0;
                  break;
                default:
                  cerr<<"Error: do not know how to double this kind of periodicity (yet): "<<zone_kind<<endl;
                  throw(-1);
                }
              }
            }
            else
            {
              // multiple periodicity should NOT involve doubling
              assert(double_transform==0);
              // this is multiple periodicity - treatment here depends on the
              // particular periodicity being combined...
              switch (zone_kind)
              {
              case FA_ZONE_PERIODIC_CART:
                switch (zone->getKind())
                {
                case FA_ZONE_PERIODIC_CART:
                  // both are Cartesian, so just add the transformations...
                  for (int i = 0; i<3; i++)
                    periodic_data[i] += zone->periodic_data[i];
                  break;
                case FA_ZONE_PERIODIC_CYL_Z:
                  // this is cylindrical periodicity, so check that
                  // the cartesian transforms are just in the z-direction...
                  assert(periodic_data[0]==0.0);
                  assert(periodic_data[1]==0.0);
                  // put the cos and sin into 0,1...
                  periodic_data[0] = zone->periodic_data[0]; // cos(thetaz)
                  periodic_data[1] = zone->periodic_data[1]; // sin(thetaz)
                  // the periodic_data[2] should already contain the dz, so
                  // just leave it...
                  assert(zone->periodic_data[2]==0.0);
                  // switch the zone to CYL_Z...
                  zone_kind = FA_ZONE_PERIODIC_CYL_Z;
                  break;
                default:
                  cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                  throw(-1);
                }
                break;
              case FA_ZONE_PERIODIC_CYL_Z:
                switch (zone->getKind())
                {
                case FA_ZONE_PERIODIC_CART:
                  // this is cartesian periodicity mixed with cyl_z, so there should
                  // be no x,y transformation...
                  assert(zone->periodic_data[0]==0.0);
                  assert(zone->periodic_data[1]==0.0);
                  // add our z-transformation...
                  periodic_data[2] += zone->periodic_data[2]; // cos(thetaz)
                  break;
                default:
                  cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                  throw(-1);
                }
                break;
              default:
                cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                throw(-1);
              }
            }
          }
        }
        assert(zone_kind!=-1);
        assert(don2_l-don_f-1<don_l-don_f);
        prcomm->addUnpackRange(don2_f-don_f, don2_l-don_f-1, inverse_bits, zone_kind, periodic_data);
      }
    }
  }

  // ------------------------------
  // pack Prcomms and Ranges...
  // ------------------------------

  // reverse exchange the unpackBufferInt into the packBufferInt...
  exchangePrcommBufferIntReverseFirst(prcommList, 2);

  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); ++prcomm)
  {

    // the packBufferInt (filled during the above reverse exchange) contains
    // the local node number for packing, and the local bits associated with any 
    // periodicity...
    for (int i = 0; i<prcomm->packIndexVec.size(); i++)
      prcomm->packIndexVec[i] = prcomm->packBufferInt[2*i];

    // for the ranges...
    int i_l = 0;
    while (i_l<prcomm->packIndexVec.size())
    {
      int i_f = i_l;
      int bits = prcomm->packBufferInt[2*i_f+1];
      // advance don2_l to the end of the bits...
      while ((i_l<prcomm->packIndexVec.size())&&(prcomm->packBufferInt[2*i_l+1]==bits))
        i_l++;
      // add the pack Range...
      if (bits==0)
      {
        prcomm->addPackRange(i_f, i_l-1, 0, FA_ZONE_INTERNAL, NULL);
      }
      else
      {
        int zone_kind = -1;
        double periodic_data[3];

        int double_transform = 0;
        if (bits&(1<<30))
        {
          double_transform = 1;
        }

        for (int bit = 0; bit<n_periodic_bit_table; bit++)
        {
          if (bits&(1<<bit))
          {
            // get this zone from the faZoneList...
            list<FaZone>::iterator zone;
            for (zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
            {
              if (zone->getBits()==(1<<bit)) break;
            }
            // make sure we found one...
            assert(zone!=faZoneList.end());
            if (zone_kind==-1)
            {
              // this is the first we found - allow any supported periodicity...
              zone_kind = zone->getKind();
              for (int i = 0; i<3; i++)
                periodic_data[i] = zone->periodic_data[i];

              // check for doubling here...
              if (double_transform==1)
              {
                switch (zone_kind)
                {
                case FA_ZONE_PERIODIC_CART:
                  for (int i = 0; i<3; i++)
                    periodic_data[i] *= 2.0;
                  break;
                default:
                  cerr<<"Error: do not know how to double this kind of periodicity (yet): "<<zone_kind<<endl;
                  throw(-1);
                }
              }

            }
            else
            {
              // multiple periodicity should NOT involve doubling
              assert(double_transform==0);
              // this is multiple periodicity - treatment here depends on the
              // particular periodicity being combined...
              switch (zone_kind)
              {
              case FA_ZONE_PERIODIC_CART:
                switch (zone->getKind())
                {
                case FA_ZONE_PERIODIC_CART:
                  // both are Cartesian, so just add the transformations...
                  for (int i = 0; i<3; i++)
                    periodic_data[i] += zone->periodic_data[i];
                  break;
                case FA_ZONE_PERIODIC_CYL_Z:
                  // this is cylindrical periodicity, so check that
                  // the cartesian transforms are just in the z-direction...
                  assert(periodic_data[0]==0.0);
                  assert(periodic_data[1]==0.0);
                  // put the cos and sin into 0,1...
                  periodic_data[0] = zone->periodic_data[0]; // cos(thetaz)
                  periodic_data[1] = zone->periodic_data[1]; // sin(thetaz)
                  // the periodic_data[2] should already contain the dz, so
                  // just leave it...
                  assert(zone->periodic_data[2]==0.0);
                  // switch the zone to CYL_Z...
                  zone_kind = FA_ZONE_PERIODIC_CYL_Z;
                  break;
                default:
                  cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                  throw(-1);
                }
                break;
              case FA_ZONE_PERIODIC_CYL_Z:
                switch (zone->getKind())
                {
                case FA_ZONE_PERIODIC_CART:
                  // this is cartesian periodicity mixed with cyl_z, so there should
                  // be no x,y transformation...
                  assert(zone->periodic_data[0]==0.0);
                  assert(zone->periodic_data[1]==0.0);
                  // add our z-transformation...
                  periodic_data[2] += zone->periodic_data[2]; // cos(thetaz)
                  break;
                default:
                  cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                  throw(-1);
                }
                break;
              default:
                cerr<<"Error: this multiple periodicity not yet supported."<<endl;
                throw(-1);
              }
            }
          }
        }
        assert(zone_kind!=-1);
        prcomm->addPackRange(i_f, i_l-1, bits, zone_kind, periodic_data);
      }
    }

  }

}

void Ugp::checkFaceCorrespondence()
{

  // do face data exchange to test transformations...

  double my_dx_max[6];
  double (*x_fa)[3] = new double[nfa][3];

  // try the forward (normal) exchanges...

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        x_fa[ifa][i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);
  }
  updateFaR2(x_fa, REPLACE_TRANSLATE_DATA);
  for (int i = 0; i<3; i++)
    my_dx_max[i] = 0.0;
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    double this_x_fa[3];
    for (int i = 0; i<3; i++)
      this_x_fa[i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        this_x_fa[i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
    {
      this_x_fa[i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);
      my_dx_max[i] = max(my_dx_max[i], fabs(this_x_fa[i]-x_fa[ifa][i]));
    }
  }

  // try the reverse exchanges...

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        x_fa[ifa][i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);
  }
  updateFaR2Reverse(x_fa, REPLACE_TRANSLATE_DATA);
  for (int i = 0; i<3; i++)
    my_dx_max[3+i] = 0.0;
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    double this_x_fa[3];
    for (int i = 0; i<3; i++)
      this_x_fa[i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        this_x_fa[i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
    {
      this_x_fa[i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);
      my_dx_max[3+i] = max(my_dx_max[3+i], fabs(this_x_fa[i]-x_fa[ifa][i]));
    }
  }

  delete[] x_fa;

  double dx_max[6];
  MPI_Reduce(my_dx_max, dx_max, 6, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank==0)
  {
    cout<<"checkFaceCorrespondence (should be zero), forward: "<<dx_max[0]<<" "<<dx_max[1]<<" "<<dx_max[2]
        <<", reverse: "<<dx_max[3]<<" "<<dx_max[4]<<" "<<dx_max[5]<<endl;
    /*
     // just in case this goes unnoticed by the user...
     for (int i = 0; i < 6; i++) {
     if ( dx_max[i] > 1.0E-10) {
     cerr << "Error: face correspondence failed" << endl;
     throw(-1);
     }
     }
     */
  }

}

void Ugp::checkNodeCorrespondence()
{

  // do node data exchange to test transformations...

  double my_dx_max[6];
  double (*x_no_copy)[3] = new double[nno][3];

  // try the forward (normal) exchanges...

  for (int ino = 0; ino<nno; ino++)
    for (int i = 0; i<3; i++)
      x_no_copy[ino][i] = x_no[ino][i];
  updateNoR2(x_no_copy, REPLACE_TRANSLATE_DATA);
  for (int i = 0; i<3; i++)
    my_dx_max[i] = 0.0;
  for (int ino = 0; ino<nno; ino++)
  {
    for (int i = 0; i<3; i++)
    {
      my_dx_max[i] = max(my_dx_max[i], fabs(x_no_copy[ino][i]-x_no[ino][i]));
    }
  }

  // try the reverse exchanges...

  for (int ino = 0; ino<nno; ino++)
    for (int i = 0; i<3; i++)
      x_no_copy[ino][i] = x_no[ino][i];
  updateNoR2Reverse(x_no_copy, REPLACE_TRANSLATE_DATA);
  for (int i = 0; i<3; i++)
    my_dx_max[3+i] = 0.0;
  for (int ino = 0; ino<nno; ino++)
  {
    for (int i = 0; i<3; i++)
    {
      my_dx_max[3+i] = max(my_dx_max[3+i], fabs(x_no_copy[ino][i]-x_no[ino][i]));
    }
  }

  delete[] x_no_copy;

  double dx_max[6];
  MPI_Reduce(my_dx_max, dx_max, 6, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank==0)
  {
    cout<<"checkNodeCorrespondence (should be zero), forward: "<<dx_max[0]<<" "<<dx_max[1]<<" "<<dx_max[2]
        <<", reverse: "<<dx_max[3]<<" "<<dx_max[4]<<" "<<dx_max[5]<<endl;
    /*
     // just in case this goes unnoticed by the user...
     for (int i = 0; i < 6; i++) {
     if ( dx_max[i] > 1.0E-10) {
     cerr << "Error: node correspondence failed" << endl;
     throw(-1);
     }
     }
     */
  }

}

void Ugp::syncNodes()
{

  if (mpi_rank==0) cout<<"syncNodes()"<<endl;

  for (int ino = 0; ino<nno; ino++)
    no_flag[ino] = 1;

  updateNoR2(x_no, ADD_TRANSLATE_DATA);
  updateNoI1(no_flag, ADD_DATA);

  for (int ino = 0; ino<nno; ino++)
    for (int i = 0; i<3; i++)
      x_no[ino][i] /= (double) no_flag[ino];

}

void Ugp::buildFaocv()
{

  if (mpi_rank==0) cout<<"buildFaocv()"<<endl;

  assert(faocv_i==NULL);
  assert(faocv_v==NULL);

  // we don't assume anything about ordering here, simply that if 
  // cvofa[ifa][0,1] is >= 0, then it represents a valid connection...

  faocv_i = new int[ncv+1];
  for (int icv = 0; icv<ncv; icv++)
    faocv_i[icv+1] = 0;

  for (int iter = 0; iter<2; iter++)
  {
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // icv0 is always valid... 
      int icv = cvofa[ifa][0];
      assert((icv>=0)&&(icv<ncv));
      if (iter==0)
      {
        faocv_i[icv+1] += 1;
      }
      else
      {
        faocv_v[faocv_i[icv]] = ifa;
        faocv_i[icv] += 1;
      }
      // icv1 may or may not be valid - here we assume that
      // the storage of a periodic face index has been removed
      // and -1 indicates a boundary icv...
      icv = cvofa[ifa][1];
      if (icv>=0)
      {
        assert(icv<ncv);
        if (iter==0)
        {
          faocv_i[icv+1] += 1;
        }
        else
        {
          faocv_v[faocv_i[icv]] = ifa;
          faocv_i[icv] += 1;
        }
      }
    }
    if (iter==0)
    {
      faocv_i[0] = 0;
      for (int icv = 0; icv<ncv; icv++)
        faocv_i[icv+1] += faocv_i[icv];
      faocv_s = faocv_i[ncv];
      faocv_v = new int[faocv_s];
    }
    else
    {
      // return csr...
      for (int icv = ncv; icv>0; icv--)
        faocv_i[icv] = faocv_i[icv-1];
      faocv_i[0] = 0;
    }
  }

}

void Ugp::buildNoocv()
{

  if (mpi_rank==0) cout<<"buildNoocv()"<<endl;

  assert(noocv_i==NULL);
  assert(noocv_v==NULL);

  noocv_i = new int[ncv+1];
  for (int icv = 0; icv<ncv; ++icv)
    noocv_i[icv+1] = 0;

  for (int iter = 0; iter<2; ++iter)
  {
    // put -1 in all nodes...
    for (int ino = 0; ino<nno; ino++)
      no_flag[ino] = -1;
    // cycle through cv's...
    for (int icv = 0; icv<ncv; icv++)
    {
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++)
        {
          int ino = noofa_v[nof];
          if (no_flag[ino]!=icv)
          {
            no_flag[ino] = icv;
            if (iter==0)
            {
              noocv_i[icv+1] += 1;
            }
            else
            {
              noocv_v[noocv_i[icv]] = ino;
              noocv_i[icv] += 1;
            }
          }
        }
      }
    }
    if (iter==0)
    {
      noocv_i[0] = 0;
      for (int icv = 0; icv<ncv; icv++)
        noocv_i[icv+1] += noocv_i[icv];
      noocv_s = noocv_i[ncv];
      noocv_v = new int[noocv_s];
    }
    else
    {
      for (int icv = ncv; icv>0; icv--)
        noocv_i[icv] = noocv_i[icv-1];
      noocv_i[0] = 0;
    }
  }

}

void Ugp::reorderFaocvNoocv()
{

  if (mpi_rank==0) cout<<"reorderFaocvNoocv()"<<endl;

  assert(noocv_i!=NULL);
  assert(noocv_v!=NULL);

  // this routine identifies basic primatives and sets their node ordering to the
  // default ordering. Other cells (cells with hanging nodes, for example) are left as is...

  int nno_cv, nfa_cv;
  double x_no_cv[NNOC_MAX][3];
  int ino_of_ino_cv[NNOC_MAX];
  int noofa_i_cv[NFOC_MAX+1];
  int noofa_v_cv[8*NFOC_MAX];

  int no_flag_cv[NNOC_MAX];
  int fa_flag_cv[NFOC_MAX];
  int ifa_of_ifa_cv[NFOC_MAX];

  for (int ino = 0; ino<nno; ino++)
    no_flag[ino] = -1;
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_flag[ifa] = -1;

  int my_type_count[5] = { 0, 0, 0, 0, 0 };

  for (int icv = 0; icv<ncv; icv++)
  {

    const int noc_f = noocv_i[icv];
    const int noc_l = noocv_i[icv+1]-1;
    nno_cv = noc_l-noc_f+1;
    assert(nno_cv<=NNOC_MAX);
    for (int noc = noc_f; noc<=noc_l; noc++)
    {
      const int ino = noocv_v[noc];
      assert(no_flag[ino]==-1);
      no_flag[ino] = noc-noc_f;
      x_no_cv[noc-noc_f][0] = x_no[ino][0];
      x_no_cv[noc-noc_f][1] = x_no[ino][1];
      x_no_cv[noc-noc_f][2] = x_no[ino][2];
      ino_of_ino_cv[noc-noc_f] = ino;
    }
    nfa_cv = 0;
    noofa_i_cv[0] = 0;

    // loop on faces...
    int nof_cv = 0;
    const int foc_f = faocv_i[icv];
    const int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc<=foc_l; foc++)
    {
      const int ifa = faocv_v[foc];
      assert(fa_flag[ifa]==-1);
      fa_flag[ifa] = foc-foc_f;
      nfa_cv++;
      assert(nfa_cv<=NFOC_MAX);
      noofa_i_cv[nfa_cv] = noofa_i_cv[nfa_cv-1]+noofa_i[ifa+1]-noofa_i[ifa];
      assert(nof_cv==noofa_i_cv[nfa_cv-1]);
      // add faces such that they are outward pointing...
      int nof_begin, nof_end, nof_inc;
      if (cvofa[ifa][0]==icv)
      {
        nof_begin = noofa_i[ifa];
        nof_end = noofa_i[ifa+1];
        nof_inc = 1;
      }
      else
      {
        assert(cvofa[ifa][1]==icv);
        nof_begin = noofa_i[ifa+1]-1;
        nof_end = noofa_i[ifa]-1;
        nof_inc = -1;
      }
      for (int nof = nof_begin; nof!=nof_end; nof += nof_inc)
      {
        const int ino = noofa_v[nof];
        assert(no_flag[ino]>=0);
        assert(nof_cv<8*NFOC_MAX);
        noofa_v_cv[nof_cv++] = no_flag[ino];
      }
      // also store for reordering faces...
      ifa_of_ifa_cv[foc-foc_f] = ifa;
    }

    int cv_type = reorderElementNodesAndFaces(no_flag_cv, fa_flag_cv, nno_cv, nfa_cv, x_no_cv, noofa_i_cv, noofa_v_cv);

    switch (cv_type)
    {
    case HEX_TYPE:
    case PRISM_TYPE:
    case PYRAMID_TYPE:
    case TET_TYPE:
    case UNKNOWN_TYPE:
      my_type_count[cv_type] += 1;
      break;
    default:
      cerr<<"Error: unknown type returned by reorderElementNodesAndFaces: "<<cv_type<<endl;
      throw(-1);
    }

    // reorder the nodes and reset the no_flag...
    for (int noc = noc_f; noc<=noc_l; ++noc)
    {
      const int ino_cv = no_flag_cv[noc-noc_f];
      const int ino = ino_of_ino_cv[ino_cv];
      assert(ino>=0);
      noocv_v[noc] = ino;
      ino_of_ino_cv[ino_cv] = -1;
      assert(no_flag[ino]>=0);
      no_flag[ino] = -1;
    }

    // reorder the faces and reset the fa_flag...
    // this is slightly different than above because the
    // fa_flag_cv contains the new location of that particular
    // face...
    for (int foc = foc_f; foc<=foc_l; ++foc)
    {
      const int foc_new = fa_flag_cv[foc-foc_f];
      const int ifa = ifa_of_ifa_cv[foc-foc_f];
      assert(ifa>=0);
      faocv_v[foc_new+foc_f] = ifa;
      assert(fa_flag[ifa]>=0);
      fa_flag[ifa] = -1;
    }

  }

  int type_count[5];
  MPI_Reduce(my_type_count, type_count, 5, MPI_INT, MPI_SUM, 0, mpi_comm);
  if (mpi_rank==0)
  {
    int type_count_sum = 0;
    for (int i = 0; i<5; i++)
      type_count_sum += type_count[i];
    cout<<"************* element summary **************"<<"\nhex:     "<<type_count[HEX_TYPE]<<"\nprism:   "
        <<type_count[PRISM_TYPE]<<"\npyramid: "<<type_count[PYRAMID_TYPE]<<"\ntet:     "<<type_count[TET_TYPE]
        <<"\nunknown: "<<type_count[UNKNOWN_TYPE]<<"\ntotal:   "<<type_count_sum
        <<"\n********** end element summary *************"<<endl;
  }

}

int Ugp::reorderElementNodesAndFaces(int * no_flag_cv, int * fa_flag_cv, const int nno_cv, const int nfa_cv,
    const double(* x_no_cv)[3], const int * noofa_i_cv, const int * noofa_v_cv)
{

  for (int ino = 0; ino<nno_cv; ++ino)
    no_flag_cv[ino] = -1;

  for (int ifa = 0; ifa<nfa_cv; ++ifa)
    fa_flag_cv[ifa] = -1;

  // count the number of quad and hex faces...
  double tri_count = 0;
  double quad_count = 0;
  for (int ifa = 0; ifa<nfa_cv; ifa++)
  {
    const int nnof = noofa_i_cv[ifa+1]-noofa_i_cv[ifa];
    if (nnof==3) ++tri_count;
    else if (nnof==4) ++quad_count;
  }

  // try to determine is this cv fits one of the standard types...
  int standard_type = -1;
  if ((nfa_cv==6)&&(nno_cv==8)&&(quad_count==6)&&(tri_count==0))
  {

    standard_type = HEX_TYPE;

    // ------------------------------------------
    // hex: standard node numbering is...
    //
    //             f3
    //        7------6   
    //       /  f5  /|
    //      4------5 | f1
    //   f0 | 3    | 2 
    //      |  f2  |/
    //      0------1
    //          f4
    // faces: use x0,x1,y0,y1,z0,z1 numbering 
    // ------------------------------------------
    // no_flag_cv is currently -1. use the getQuadFaceMatching routines to find the face
    // matching the nodes that are set and set the remaining nodes with its node values...
    // the node numbering uses an outward pointing right-hand rule...

    int ifa = getQuadFaceMatching(no_flag_cv[0], no_flag_cv[3], no_flag_cv[2], no_flag_cv[1], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 4;
    ifa = getQuadFaceMatching(no_flag_cv[0], no_flag_cv[1], no_flag_cv[5], no_flag_cv[4], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 2;
    ifa = getQuadFaceMatching(no_flag_cv[1], no_flag_cv[2], no_flag_cv[6], no_flag_cv[5], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 1;
    ifa = getQuadFaceMatching(no_flag_cv[2], no_flag_cv[3], no_flag_cv[7], no_flag_cv[6], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 3;
    ifa = getQuadFaceMatching(no_flag_cv[3], no_flag_cv[0], no_flag_cv[4], no_flag_cv[7], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 0;
    ifa = getQuadFaceMatching(no_flag_cv[4], no_flag_cv[5], no_flag_cv[6], no_flag_cv[7], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 5;

    // check...
    for (int ifa = 0; ifa<nfa_cv; ifa++)
      assert(fa_flag_cv[ifa]!=-1);

  }
  else if ((nfa_cv==5)&&(nno_cv==6)&&(quad_count==3)&&(tri_count==2))
  {

    standard_type = PRISM_TYPE;

    // ------------------------------------------
    // prism: standard node numbering is...
    //       
    //           __
    //       ___/  5
    //    __/     /|
    //   3-------4 |
    //   |       | 2
    //   |       |/
    //   0-------1
    //
    // ------------------------------------------
    int ifa = getTriFaceMatching(no_flag_cv[0], no_flag_cv[2], no_flag_cv[1], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 3;
    ifa = getQuadFaceMatching(no_flag_cv[0], no_flag_cv[1], no_flag_cv[4], no_flag_cv[3], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 0;
    ifa = getQuadFaceMatching(no_flag_cv[1], no_flag_cv[2], no_flag_cv[5], no_flag_cv[4], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 1;
    ifa = getQuadFaceMatching(no_flag_cv[2], no_flag_cv[0], no_flag_cv[3], no_flag_cv[5], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 2;
    ifa = getTriFaceMatching(no_flag_cv[3], no_flag_cv[4], no_flag_cv[5], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 4;

    // check...
    for (int ifa = 0; ifa<nfa_cv; ifa++)
      assert(fa_flag_cv[ifa]!=-1);

  }
  else if ((nfa_cv==5)&&(nno_cv==5)&&(quad_count==1)&&(tri_count==4))
  {

    standard_type = PYRAMID_TYPE;

    // ------------------------------------------
    // pyramid: standard node numbering is...
    // 
    //       4_
    //      / \\_
    //     / 3 \ \_2
    //    /     \ /
    //   0-------1
    //
    // -----------------------------------------
    int ifa = getQuadFaceMatching(no_flag_cv[0], no_flag_cv[3], no_flag_cv[2], no_flag_cv[1], nfa_cv, noofa_i_cv,
        noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 4;
    ifa = getTriFaceMatching(no_flag_cv[0], no_flag_cv[1], no_flag_cv[4], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 0;
    ifa = getTriFaceMatching(no_flag_cv[1], no_flag_cv[2], no_flag_cv[4], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 1;
    ifa = getTriFaceMatching(no_flag_cv[2], no_flag_cv[3], no_flag_cv[4], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 2;
    ifa = getTriFaceMatching(no_flag_cv[3], no_flag_cv[0], no_flag_cv[4], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 3;

    // check...
    for (int ifa = 0; ifa<nfa_cv; ifa++)
      assert(fa_flag_cv[ifa]!=-1);

  }
  else if ((nfa_cv==4)&&(nno_cv==4)&&(quad_count==0)&&(tri_count==4))
  {

    standard_type = TET_TYPE;

    // ------------------------------------------
    // tet: standard node numbering is...
    // 
    //       3_
    //      / \\_
    //     /   \ \_2
    //    /     \ /
    //   0-------1
    //
    // -----------------------------------------
    int ifa = getTriFaceMatching(no_flag_cv[0], no_flag_cv[2], no_flag_cv[1], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 3;
    ifa = getTriFaceMatching(no_flag_cv[0], no_flag_cv[1], no_flag_cv[3], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 2;
    ifa = getTriFaceMatching(no_flag_cv[1], no_flag_cv[2], no_flag_cv[3], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 0;
    ifa = getTriFaceMatching(no_flag_cv[2], no_flag_cv[0], no_flag_cv[3], nfa_cv, noofa_i_cv, noofa_v_cv);
    assert(ifa>=0);
    assert(fa_flag_cv[ifa]==-1);
    fa_flag_cv[ifa] = 1;

    // check...
    for (int ifa = 0; ifa<nfa_cv; ifa++)
      assert(fa_flag_cv[ifa]!=-1);

  }
  else
  {

    standard_type = UNKNOWN_TYPE;

    // don't do any reordering...

    for (int ino = 0; ino<nno_cv; ino++)
      no_flag_cv[ino] = ino;

    for (int ifa = 0; ifa<nfa_cv; ++ifa)
      fa_flag_cv[ifa] = ifa;

  }

  return (standard_type);

}

int Ugp::getQuadFaceMatching(int &ino0, int &ino1, int &ino2, int &ino3, const int nfa_ss, const int * noofa_i_ss,
    const int * noofa_v_ss)
{

  // find the quad face in the list matching the passed nodes

  for (int ifa = 0; ifa<nfa_ss; ifa++)
  {
    int nof_f = noofa_i_ss[ifa];
    int nof_l = noofa_i_ss[ifa+1]-1;
    // consider quads only...
    if (nof_l-nof_f+1==4)
    {
      for (int nof0 = nof_f; nof0<=nof_l; nof0++)
      {
        if ((ino0==-1)||(ino0==noofa_v_ss[nof0]))
        {
          int nof1 = nof0+1;
          if (nof1>nof_l) nof1 = nof_f;
          if ((ino1==-1)||(ino1==noofa_v_ss[nof1]))
          {
            int nof2 = nof1+1;
            if (nof2>nof_l) nof2 = nof_f;
            if ((ino2==-1)||(ino2==noofa_v_ss[nof2]))
            {
              int nof3 = nof2+1;
              if (nof3>nof_l) nof3 = nof_f;
              if ((ino3==-1)||(ino3==noofa_v_ss[nof3]))
              {
                // matched!
                // set all nodes (some/all may already be set -- doesn't matter)...
                ino0 = noofa_v_ss[nof0];
                ino1 = noofa_v_ss[nof1];
                ino2 = noofa_v_ss[nof2];
                ino3 = noofa_v_ss[nof3];
                // and return...
                return (ifa);
              }
            }
          }
        }
      }
    }
  }

  return (-1);

}

int Ugp::getTriFaceMatching(int &ino0, int &ino1, int &ino2, const int nfa_ss, const int * noofa_i_ss,
    const int * noofa_v_ss)
{

  // find the quad face in the list matching the passed nodes

  for (int ifa = 0; ifa<nfa_ss; ifa++)
  {
    int nof_f = noofa_i_ss[ifa];
    int nof_l = noofa_i_ss[ifa+1]-1;
    // consider quads only...
    if (nof_l-nof_f+1==3)
    {
      for (int nof0 = nof_f; nof0<=nof_l; nof0++)
      {
        if ((ino0==-1)||(ino0==noofa_v_ss[nof0]))
        {
          int nof1 = nof0+1;
          if (nof1>nof_l) nof1 = nof_f;
          if ((ino1==-1)||(ino1==noofa_v_ss[nof1]))
          {
            int nof2 = nof1+1;
            if (nof2>nof_l) nof2 = nof_f;
            if ((ino2==-1)||(ino2==noofa_v_ss[nof2]))
            {
              // matched!
              // set all nodes (some/all may already be set -- doesn't matter)...
              ino0 = noofa_v_ss[nof0];
              ino1 = noofa_v_ss[nof1];
              ino2 = noofa_v_ss[nof2];
              // and return...
              return (ifa);
            }
          }
        }
      }
    }
  }

  return (-1);

}

void Ugp::checkCvonoCounts()
{

  // the nodes of a triply-periodic hex mesh are surrounded 
  // by exactly 8 cv's. check this.

  // XXXX this result is not partition independent because it adds
  // counts from all nodes on all partitions, and of course there
  // are nodal duplicates...

  // this would be easier with noocv_i/v, but use faocv instead...
  assert(faocv_i!=NULL);
  assert(faocv_v!=NULL);

  for (int ino = 0; ino<nno; ino++)
    no_flag[ino] = 0;

  for (int icv = 0; icv<ncv; icv++)
  {
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc<=foc_l; foc++)
    {
      int ifa = faocv_v[foc];
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        int ino = noofa_v[nof];
        if (no_flag[ino]>=0) no_flag[ino] = -no_flag[ino]-1;
      }
    }
    for (int foc = foc_f; foc<=foc_l; foc++)
    {
      int ifa = faocv_v[foc];
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        int ino = noofa_v[nof];
        if (no_flag[ino]<0) no_flag[ino] = -no_flag[ino];
      }
    }
  }

  updateNoI1(no_flag, ADD_DATA);

  int my_cvono_count_max = 0;
  for (int ino = 0; ino<nno; ino++)
    my_cvono_count_max = max(my_cvono_count_max, no_flag[ino]);

  int cvono_count_max;
  MPI_Allreduce(&my_cvono_count_max, &cvono_count_max, 1, MPI_INT, MPI_MAX, mpi_comm);

  int * my_bin_count = new int[cvono_count_max+1];
  for (int i = 0; i<=cvono_count_max; i++)
    my_bin_count[i] = 0;

  for (int ino = 0; ino<nno; ino++)
    my_bin_count[no_flag[ino]] += 1;

  int * bin_count = new int[cvono_count_max+1];
  MPI_Reduce(my_bin_count, bin_count, cvono_count_max+1, MPI_INT, MPI_SUM, 0, mpi_comm);
  if (mpi_rank==0)
  {
    cout<<"checkCvonoCounts: "<<endl;
    for (int i = 0; i<=cvono_count_max; i++)
      cout<<" n, global count: "<<i<<", "<<bin_count[i]<<endl;
  }

  delete[] my_bin_count;
  delete[] bin_count;

}

void Ugp::writeBoundaryFacesTecplot(const char * prefix)
{

  // set the flag in all faces to zero...
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_flag[ifa] = 0;

  // loop on all face zones and set the flag for boundary zone faces to 1...
  for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone!=faZoneList.end(); faZone++)
  {
    if ((faZone->getKind()==FA_ZONE_BOUNDARY)||((faZone->getKind()>=FA_ZONE_PERIODIC_FIRST)&&(faZone->getKind()
        <=FA_ZONE_PERIODIC_LAST)))
    {

      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ifa++)
      {
        fa_flag[ifa] = 1;
      }

      char filename[128];
      sprintf(filename, "%s.%s.dat", prefix, faZone->getName());

      writeFlaggedFacesTecplot(filename);

      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ifa++)
      {
        fa_flag[ifa] = 0;
      }

    }
  }

}

void Ugp::writeBoundaryFacesStl()
{

  // recall format...
  /*
   solid
   facet normal 1 1 1
   outer loop
   vertex 0.570125 -0.820108 -0.015980
   vertex 0.720667 -0.689533 -0.050486
   vertex 0.602396 -0.770875 -0.200511
   endloop
   endfacet
   ...
   */

  for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone!=faZoneList.end(); faZone++)
  {
    if (faZone->getKind()==FA_ZONE_BOUNDARY)
    {

      char filename[128];
      sprintf(filename, "%s.stl", faZone->getName());

      FILE * fp;
      if (mpi_rank==0)
      {
        if ((fp = fopen(filename, "w"))==NULL)
        {
          cerr<<"Error: cannot open file "<<filename<<endl;
          throw(-1);
        }
        fprintf(fp, "solid\n");
      }
      else
      {
        int dummy;
        MPI_Status status;
        MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1234, mpi_comm, &status);
        if ((fp = fopen(filename, "a"))==NULL)
        {
          cerr<<"Error: cannot open file "<<filename<<endl;
          throw(-1);
        }
      }

      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ifa++)
      {
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        switch (nof_l-nof_f+1)
        {
        case 3:
          // a single tri...
          fprintf(fp, "facet normal 1 1 1\n");
          fprintf(fp, "outer loop\n");
          {
            int ino0 = noofa_v[nof_f];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino0][0], x_no[ino0][1], x_no[ino0][2]);
            int ino1 = noofa_v[nof_f+1];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino1][0], x_no[ino1][1], x_no[ino1][2]);
            int ino2 = noofa_v[nof_f+2];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino2][0], x_no[ino2][1], x_no[ino2][2]);
          }
          fprintf(fp, "endloop\n");
          fprintf(fp, "endfacet\n");
          break;
        case 4:
          // 2 tris
          fprintf(fp, "facet normal 1 1 1\n");
          fprintf(fp, "outer loop\n");
          {
            int ino0 = noofa_v[nof_f];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino0][0], x_no[ino0][1], x_no[ino0][2]);
            int ino1 = noofa_v[nof_f+1];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino1][0], x_no[ino1][1], x_no[ino1][2]);
            int ino2 = noofa_v[nof_f+2];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino2][0], x_no[ino2][1], x_no[ino2][2]);
          }
          fprintf(fp, "endloop\n");
          fprintf(fp, "endfacet\n");
          // number 2...
          fprintf(fp, "facet normal 1 1 1\n");
          fprintf(fp, "outer loop\n");
          {
            int ino0 = noofa_v[nof_f];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino0][0], x_no[ino0][1], x_no[ino0][2]);
            int ino1 = noofa_v[nof_f+2];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino1][0], x_no[ino1][1], x_no[ino1][2]);
            int ino2 = noofa_v[nof_f+3];
            fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino2][0], x_no[ino2][1], x_no[ino2][2]);
          }
          fprintf(fp, "endloop\n");
          fprintf(fp, "endfacet\n");
          break;
        default:
        {
          double xc[3] = { 0.0, 0.0, 0.0 };
          for (int nof = nof_f; nof<=nof_l; ++nof)
          {
            int ino = noofa_v[nof];
            FOR_I3
              xc[i] += x_no[ino][i];
          }
          FOR_I3
            xc[i] /= (double) (nof_l-nof_f+1);
          int ino1 = noofa_v[nof_l];
          for (int nof = nof_f; nof<=nof_l; ++nof)
          {
            int ino0 = ino1;
            ino1 = noofa_v[nof];
            fprintf(fp, "facet normal 1 1 1\n");
            fprintf(fp, "outer loop\n");
            {
              fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino0][0], x_no[ino0][1], x_no[ino0][2]);
              fprintf(fp, "vertex %lf %lf %lf\n", x_no[ino1][0], x_no[ino1][1], x_no[ino1][2]);
              fprintf(fp, "vertex %lf %lf %lf\n", xc[0], xc[1], xc[2]);
            }
            fprintf(fp, "endloop\n");
            fprintf(fp, "endfacet\n");
          }
        }
        }
      }

      if (mpi_rank==mpi_size-1)
      {
        fprintf(fp, "endsolid\n");
        fclose(fp);
      }
      else
      {
        fclose(fp);
        int dummy = 1;
        MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1234, mpi_comm);
      }
    }
  }

  MPI_Barrier(mpi_comm);

}

void Ugp::writeFlaggedFacesTecplot(char * filename)
{

  // a sequential tecplot file writing...
  // write also flagged registered data

  if (mpi_rank==0) cout<<"writeFlaggedFacesTecplot: "<<filename<<endl;

  FILE * fp;
  if (mpi_rank==0)
  {
    if ((fp = fopen(filename, "w"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }
    fprintf(fp, "TITLE = \"flagged faces\"\n");
    fprintf(fp, "VARIABLES = \"X\"\n");
    fprintf(fp, "\"Y\"\n");
    fprintf(fp, "\"Z\"\n");

    // add header for flagged registered data (scalars...add vectors later)
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if ((data->getDatatype()==NO_DATA)&&(data->checkFlag()))
      {
        // mixture fraction "Z" causes a problem with name correspondence...
        if ((strcmp(data->getName(), "X")==0)||(strcmp(data->getName(), "Y")==0)||(strcmp(data->getName(), "Z")==0))
        {
          fprintf(fp, "\"%s (scalar)\"\n", data->getName());
        }
        else
        {
          fprintf(fp, "\"%s\"\n", data->getName());
        }
      }
    }

  }
  else
  {
    int dummy;
    MPI_Status status;
    MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1234, mpi_comm, &status);
    if ((fp = fopen(filename, "a"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }
  }

  int ncount = 0; // node count
  int ecount = 0; // element count

  for (int ifa = 0; ifa<getNfa(); ifa++)
  {
    if (fa_flag[ifa]!=0)
    {
      int nnof = noofa_i[ifa+1]-noofa_i[ifa];
      assert(nnof>=3);
      ncount += nnof+1; // all nodes plus the center...
      ecount += nnof; // center-edge tri
    }
  }

  if (ncount>0)
  {

    fprintf(fp, "ZONE T=\"%s, rank = %d\"\n", filename, mpi_rank);
    fprintf(fp, "N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE\n", ncount, ecount);

    // nodes first...
    for (int i = 0; i<3; i++)
    { // the node coordinates
      int lines_written = 0;
      for (int ifa = 0; ifa<getNfa(); ifa++)
      {
        if (fa_flag[ifa]!=0)
        {
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          double x_fa = 0.0;
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino = noofa_v[nof];
            assert((ino>=0)&&(ino<nno));
            x_fa += x_no[ino][i];
          }
          x_fa /= (double) (nof_l-nof_f+1);
          fprintf(fp, "%lf ", x_fa);
          if (++lines_written>5)
          {
            fprintf(fp, "\n");
            lines_written = 0;
          }
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino = noofa_v[nof];
            fprintf(fp, "%lf ", x_no[ino][i]);
            if (++lines_written>5)
            {
              fprintf(fp, "\n");
              lines_written = 0;
            }
          }
        }
      }
      if (lines_written!=0) fprintf(fp, "\n");
    }

    // now the registered data (scalars)
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      int lines_written = 0;
      if ((data->getDatatype()==NO_DATA)&&(data->checkFlag()))
      {
        for (int ifa = 0; ifa<getNfa(); ifa++)
        {
          if (fa_flag[ifa]!=0)
          {
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            double phi_fa = 0.0;
            for (int nof = nof_f; nof<=nof_l; nof++)
            {
              int ino = noofa_v[nof];
              assert((ino>=0)&&(ino<nno));
              phi_fa += (*data->ptr)[ino];
            }
            phi_fa /= (double) (nof_l-nof_f+1);
            fprintf(fp, "%lf ", phi_fa);
            if (++lines_written>5)
            {
              fprintf(fp, "\n");
              lines_written = 0;
            }
            for (int nof = nof_f; nof<=nof_l; nof++)
            {
              int ino = noofa_v[nof];
              fprintf(fp, "%lf ", (*data->ptr)[ino]);
              if (++lines_written>5)
              {
                fprintf(fp, "\n");
                lines_written = 0;
              }
            }
          }
        }
      }
      if (lines_written!=0) fprintf(fp, "\n");
    }

    // connectivity...
    ncount = 0;
    for (int ifa = 0; ifa<getNfa(); ifa++)
    {
      if (fa_flag[ifa]!=0)
      {
        int ino_c = ncount;
        ncount += noofa_i[ifa+1]-noofa_i[ifa]+1; // all nodes plus the center...
        // cycle through the edges...
        int ino2 = ncount-1;
        for (int ino = ino_c+1; ino<ncount; ino++)
        {
          int ino1 = ino2;
          ino2 = ino;
          fprintf(fp, "%d %d %d\n", ino_c+1, ino1+1, ino2+1); // 1-indexing
        }
      }
    }

  }

  fclose(fp);

  if (mpi_rank<mpi_size-1)
  {
    int dummy = 1;
    MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1234, mpi_comm);
  }

  MPI_Barrier(mpi_comm);

}

void Ugp::writeFlaggedCvsTecplotASCII(char * filename)
{

  // a sequential tecplot file writing...
  // write also flagged registered data

  if (mpi_rank==0) cout<<"writeFlaggedCvsTecplot: "<<filename<<endl;

  // step 1: build a global node numbering...
  // for nodes associated with flagged cv's...

  // build a global "node" numbering for the cv's, faces, nodes of flagged cv's...

  int my_buf[2];
  my_buf[0] = 0; // local ncount
  my_buf[1] = 0; // local ecount

  // use the connectivity to determine what has to be dumped...

  for (int ino = 0; ino<nno; ++ino)
    no_flag[ino] = -1;

  for (int ifa = 0; ifa<nfa; ++ifa)
    fa_flag[ifa] = -1;

  for (int icv = 0; icv<ncv; ++icv)
  {
    if (cv_flag[icv]==0)
    {
      // a negative 2 means skip completely...
      cv_flag[icv] = -2;
    }
    else
    {
      // at first assume we do not have to dump a value
      // at this cv centroid...
      cv_flag[icv] = -1;
      // all nodes have to be dumped...
      const int noc_f = noocv_i[icv];
      const int noc_l = noocv_i[icv+1]-1;
      const int nnoc = noc_l-noc_f+1;
      for (int noc = noc_f; noc<=noc_l; ++noc)
      {
        const int ino = noocv_v[noc];
        if (no_flag[ino]==-1)
        {
          no_flag[ino] = 0; // count later
        }
      }
      // faces may or may not be dumped...
      int ntri = 0;
      int nquad = 0;
      int ne = 0;
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      const int nfoc = foc_l-foc_f+1;
      for (int foc = foc_f; foc<=foc_l; ++foc)
      {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l-nof_f+1;
        if (fa_flag[ifa]==-1)
        {
          // this face has not yet been visited, so we need to check
          // if there are 3 or 4 nodes, and if 4, if they are coplanar...
          switch (nnof)
          {
          case 3:
            // a tri face - never a problem...
            fa_flag[ifa] = -2;
            break;
          case 4:
            // a quad face - only a problem if non-coplanar...
          {
            const int ino0 = noofa_v[nof_f];
            const int ino1 = noofa_v[nof_f+1];
            const int ino2 = noofa_v[nof_f+2];
            const int ino3 = noofa_v[nof_f+3];
            double v1[3] = { x_no[ino1][0]-x_no[ino0][0], x_no[ino1][1]-x_no[ino0][1], x_no[ino1][2]-x_no[ino0][2] };
            double v2[3] = { x_no[ino2][0]-x_no[ino0][0], x_no[ino2][1]-x_no[ino0][1], x_no[ino2][2]-x_no[ino0][2] };
            double normal[3] = { v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] };
            // check the distance of the 4th node off the plane of the other 3...
            double delta = (x_no[ino3][0]-x_no[ino0][0])*normal[0]+(x_no[ino3][1]-x_no[ino0][1])*normal[1]
                +(x_no[ino3][2]-x_no[ino0][2])*normal[2];
            if (fabs(delta)<0.01*pow(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2], 0.75))
            {
              // coplanar within 1%...
              fa_flag[ifa] = -2;
            }
            else
            {
              // not coplanar enough...
              //cout << "Warning: got non-coplanar quad face" << endl;
              fa_flag[ifa] = 0;
            }
          }
            break;
          default:
            // this is a face with one or more hanging nodes? need to dump the
            // center and triangulate to the edges...
            fa_flag[ifa] = 0;
            break;
          }
        }

        // now whether we visited it this time or not, we need to
        // what kind of face it is...

        if (fa_flag[ifa]==-2)
        {
          // tri or valid quad...
          if (nnof==3) ntri++;
          else
          {
            assert(nnof==4);
            nquad++;
          }
          // only one element here...
          ne += 1;
        }
        else
        {
          // if the face is zero, it is getting dumped, and we
          // need to add one element for each edge...
          assert(fa_flag[ifa]==0);
          ne += nnof;
        }

      }

      // now decide whether we dump the cv... 
      if ((nnoc==8)&&(nfoc==6)&&(nquad==6)&&(ntri==0)&&(ne==6))
      {
        // hex cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==6)&&(nfoc==5)&&(nquad==3)&&(ntri==2)&&(ne==5))
      {
        // prism cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==5)&&(nfoc==5)&&(nquad==1)&&(ntri==4)&&(ne==5))
      {
        // pyramid cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==4)&&(nfoc==4)&&(nquad==0)&&(ntri==4)&&(ne==4))
      {
        // tet cv...
        my_buf[1] += 1;
      }
      else
      {
        cv_flag[icv] = 0;
        my_buf[1] += ne;
      }
    }
  }

  // now my_buf[1] constains the element count.
  // cycle through nodes, faces, cells to get the vertex count...

  for (int ino = 0; ino<nno; ++ino)
  {
    if (no_flag[ino]==0)
    {
      no_flag[ino] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  for (int ifa = 0; ifa<nfa; ++ifa)
  {
    if (fa_flag[ifa]==0)
    {
      fa_flag[ifa] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  for (int icv = 0; icv<ncv; ++icv)
  {
    if (cv_flag[icv]==0)
    {
      cv_flag[icv] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  // figure out our nodal displacement...

  int disp;
  MPI_Scan(&(my_buf[0]), &disp, 1, MPI_INT, MPI_SUM, mpi_comm);
  disp -= my_buf[0];

  for (int ino = 0; ino<nno; ++ino)
    if (no_flag[ino]>=0) no_flag[ino] += disp;

  for (int ifa = 0; ifa<nfa; ++ifa)
    if (fa_flag[ifa]>=0) fa_flag[ifa] += disp;

  for (int icv = 0; icv<ncv; ++icv)
    if (cv_flag[icv]>=0) cv_flag[icv] += disp;

  // send our counts to rank 0...
  int buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_INT, MPI_SUM, 0, mpi_comm);
  if (mpi_rank==0)
  {

    // rank 0 writes the header...

    FILE * fp;
    if ((fp = fopen(filename, "w"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }
    fprintf(fp, "TITLE = \"flagged cvs\"\n");
    fprintf(fp, "VARIABLES = \"X\"\n");
    fprintf(fp, "\"Y\"\n");
    fprintf(fp, "\"Z\"\n");

    // add header for flagged registered data...
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
        case NO_DATA:
        case CV_DATA:
          // mixtrure fraction "Z" causes a problem with name correspondence...
          if ((strcmp(data->getName(), "X")==0)||(strcmp(data->getName(), "Y")==0)||(strcmp(data->getName(), "Z")==0))
          {
            fprintf(fp, "\"%s (scalar)\"\n", data->getName());
          }
          else
          {
            fprintf(fp, "\"%s\"\n", data->getName());
          }
          break;
        default:
          if (mpi_rank==0) cout<<"Warning: writeFlaggedCvsTecplot does not support data type: "<<data->getDatatype()
              <<", skipping."<<endl;
        }
      }
    }

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
        case NO_DATA:
        case CV_DATA:
          fprintf(fp, "\"%s-X\"\n", data->getName());
          fprintf(fp, "\"%s-Y\"\n", data->getName());
          fprintf(fp, "\"%s-Z\"\n", data->getName());
          break;
        default:
          if (mpi_rank==0) cout<<"Warning: writeFlaggedCvsTecplot does not support data type: "<<data->getDatatype()
              <<", skipping."<<endl;
        }
      }
    }

    fprintf(fp, "ZONE T=\"writeFlaggedCvsTecplot, np=%d\"\n", mpi_size);
    fprintf(fp, "N=%d, E=%d, F=FEBLOCK, ET=BRICK\n", buf[0], buf[1]);
    fclose(fp);

  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // x,y,z...
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  for (int i = 0; i<3; i++)
  {

    int lines_written = 0;
    if (mpi_rank==0)
    {
      // on the first coordinate, start writing right away...
      if (i>0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
      }
    }
    else
    {
      // processes wait for previous to finish, and recv's lines_written...
      MPI_Status status;
      MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
    }

    FILE * fp;
    if ((fp = fopen(filename, "a"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }

    // nodes are first...
    for (int ino = 0; ino<nno; ino++)
    {
      if (no_flag[ino]>=0)
      {
        fprintf(fp, "%18.15le ", x_no[ino][i]);
        if (++lines_written>5)
        {
          fprintf(fp, "\n");
          lines_written = 0;
        }
      }
    }

    // then faces...
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_flag[ifa]>=0)
      {
        double x_fa = 0.0;
        double weight = 0.0;
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++)
        {
          int ino = noofa_v[nof];
          x_fa += x_no[ino][i];
          weight += 1.0;
        }
        x_fa /= weight;
        fprintf(fp, "%18.15le ", x_fa);
        if (++lines_written>5)
        {
          fprintf(fp, "\n");
          lines_written = 0;
        }
      }
    }

    // then cv's...
    for (int icv = 0; icv<ncv; icv++)
    {
      if (cv_flag[icv]>=0)
      {
        double x_cv = 0.0;
        double weight = 0.0;
        int noc_f = noocv_i[icv];
        int noc_l = noocv_i[icv+1]-1;
        for (int noc = noc_f; noc<=noc_l; noc++)
        {
          int ino = noocv_v[noc];
          x_cv += x_no[ino][i];
          weight += 1.0;
        }
        x_cv /= weight;
        fprintf(fp, "%18.15le ", x_cv);
        if (++lines_written>5)
        {
          fprintf(fp, "\n");
          lines_written = 0;
        }
      }

    }

    if (mpi_rank<mpi_size-1)
    {
      // for ranks that are not last, just close and send a message on...
      fclose(fp);
      MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
    }
    else
    {
      // the last rank writes a return if required, then sends a message to the first...
      if (lines_written!=0) fprintf(fp, "\n");
      fclose(fp);
      int dummy_int = 0;
      MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
    }

  }

  // used for cv data interpolation...
  double * no_scalar = NULL;
  double * no_weight = NULL;
  double * no_scalar_ptr = NULL;

  // now the registered data (scalars)...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
      case NO_DATA:
        // ======================== node data ==========================
        no_scalar_ptr = (*data->ptr);
        break;
      case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_scalar==NULL) no_scalar = new double[nno];
        if (no_weight==NULL) no_weight = new double[nno];
        for (int ino = 0; ino<nno; ino++)
        {
          no_scalar[ino] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv<ncv; icv++)
        {
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            no_scalar[ino] += (*data->ptr)[icv];
            no_weight[ino] += 1.0;
          }
        }
        updateNoR1(no_scalar, ADD_DATA);
        updateNoR1(no_weight, ADD_DATA);
        for (int ino = 0; ino<nno; ino++)
          no_scalar[ino] /= no_weight[ino];
        no_scalar_ptr = no_scalar;
        break;
      default:
        continue;
      }

      // if we made it here, then we have some node data pointed to in no_scalar_ptr,
      // and can dump using no, fa, then cv writing...

      int lines_written = 0;
      if (mpi_rank==0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
      }
      else
      {
        // processes wait for previous to finish, and recv's lines_written...
        MPI_Status status;
        MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
      }

      FILE * fp;
      if ((fp = fopen(filename, "a"))==NULL)
      {
        cerr<<"Error: cannot open file "<<filename<<endl;
        throw(-1);
      }

      // nodes are first...
      for (int ino = 0; ino<nno; ino++)
      {
        if (no_flag[ino]>=0)
        {
          fprintf(fp, "%18.15le ", no_scalar_ptr[ino]);
          if (++lines_written>5)
          {
            fprintf(fp, "\n");
            lines_written = 0;
          }
        }
      }

      // then faces...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]>=0)
        {
          double phi_fa = 0.0;
          double weight = 0.0;
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino = noofa_v[nof];
            phi_fa += no_scalar_ptr[ino];
            weight += 1.0;
          }
          phi_fa /= weight;
          fprintf(fp, "%18.15le ", phi_fa);
          if (++lines_written>5)
          {
            fprintf(fp, "\n");
            lines_written = 0;
          }
        }
      }

      // then cv's...
      for (int icv = 0; icv<ncv; icv++)
      {
        if (cv_flag[icv]>=0)
        {
          double phi_cv = 0.0;
          double weight = 0.0;
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            phi_cv += no_scalar_ptr[ino];
            weight += 1.0;
          }
          phi_cv /= weight;
          fprintf(fp, "%18.15le ", phi_cv);
          if (++lines_written>5)
          {
            fprintf(fp, "\n");
            lines_written = 0;
          }
        }
      }

      if (mpi_rank<mpi_size-1)
      {
        // for ranks that are not last, just close and send a message on...
        fclose(fp);
        MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
      }
      else
      {
        // the last rank writes a return if required, then sends a message to the first...
        if (lines_written!=0) fprintf(fp, "\n");
        fclose(fp);
        int dummy_int = 0;
        MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
      }

    }

  }

  // cleanup and prepart for vector interpolation...
  if (no_scalar!=NULL) delete[] no_scalar;

  double (*no_vector)[3] = NULL;
  double (*no_vector_ptr)[3] = NULL;

  // loop on vector data...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
      case NO_DATA:
        // ================= nodal data ======================
        no_vector_ptr = (*data->ptr);
        break;
      case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_vector==NULL) no_vector = new double[nno][3];
        if (no_weight==NULL) no_weight = new double[nno];
        for (int ino = 0; ino<nno; ino++)
        {
          for (int i = 0; i<3; i++)
            no_vector[ino][i] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv<ncv; icv++)
        {
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            for (int i = 0; i<3; i++)
              no_vector[ino][i] += (*data->ptr)[icv][i];
            no_weight[ino] += 1.0;
          }
        }
        updateNoR2(no_vector, ADD_ROTATE_DATA);
        updateNoR1(no_weight, ADD_DATA);
        for (int ino = 0; ino<nno; ino++)
          for (int i = 0; i<3; i++)
            no_vector[ino][i] /= no_weight[ino];
        no_vector_ptr = no_vector;
        break;
      default:
        continue;
      }

      // write the data...
      for (int i = 0; i<3; i++)
      {

        int lines_written = 0;
        if (mpi_rank==0)
        {
          MPI_Status status;
          int dummy_int;
          MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
        }
        else
        {
          // processes wait for previous to finish, and recv's lines_written...
          MPI_Status status;
          MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
        }

        FILE * fp;
        if ((fp = fopen(filename, "a"))==NULL)
        {
          cerr<<"Error: cannot open file "<<filename<<endl;
          throw(-1);
        }

        // nodes are first...
        for (int ino = 0; ino<nno; ino++)
        {
          if (no_flag[ino]>=0)
          {
            fprintf(fp, "%18.15le ", no_vector_ptr[ino][i]);
            if (++lines_written>5)
            {
              fprintf(fp, "\n");
              lines_written = 0;
            }
          }
        }

        // then faces...
        for (int ifa = 0; ifa<nfa; ifa++)
        {
          if (fa_flag[ifa]>=0)
          {
            double phi_fa = 0.0;
            double weight = 0.0;
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof<=nof_l; nof++)
            {
              int ino = noofa_v[nof];
              phi_fa += no_vector_ptr[ino][i];
              weight += 1.0;
            }
            phi_fa /= weight;
            fprintf(fp, "%18.15le ", phi_fa);
            if (++lines_written>5)
            {
              fprintf(fp, "\n");
              lines_written = 0;
            }
          }
        }

        // then cv's...
        for (int icv = 0; icv<ncv; icv++)
        {
          if (cv_flag[icv]>=0)
          {
            double phi_cv = 0.0;
            double weight = 0.0;
            int noc_f = noocv_i[icv];
            int noc_l = noocv_i[icv+1]-1;
            for (int noc = noc_f; noc<=noc_l; noc++)
            {
              int ino = noocv_v[noc];
              phi_cv += no_vector_ptr[ino][i];
              weight += 1.0;
            }
            phi_cv /= weight;
            fprintf(fp, "%18.15le ", phi_cv);
            if (++lines_written>5)
            {
              fprintf(fp, "\n");
              lines_written = 0;
            }
          }
        }

        if (mpi_rank<mpi_size-1)
        {
          // for ranks that are not last, just close and send a message on...
          fclose(fp);
          MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
        }
        else
        {
          // the last rank writes a return if required, then sends a message to the first...
          if (lines_written!=0) fprintf(fp, "\n");
          fclose(fp);
          int dummy_int = 0;
          MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
        }
      }
    }
  }

  if (no_vector!=NULL) delete[] no_vector;
  if (no_weight!=NULL) delete[] no_weight;

  // ======================================
  // connectivity...
  // ======================================

  if (mpi_rank==0)
  {
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
  }
  else
  {
    // processes wait for previous to finish, and recv's lines_written...
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
  }

  FILE * fp;
  if ((fp = fopen(filename, "a"))==NULL)
  {
    cerr<<"Error: cannot open file "<<filename<<endl;
    throw(-1);
  }

  for (int icv = 0; icv<ncv; icv++)
  {
    if (cv_flag[icv]==-1)
    {
      // this is a cv that is a standard primative. We can figure out
      // its type from its node count...
      int noc_f = noocv_i[icv];
      int noc_l = noocv_i[icv+1]-1;
      int nnoc = noc_l-noc_f+1;
      switch (nnoc)
      {
      case 8:
        // hex...
        fprintf(fp,
            "%d %d %d %d %d %d %d %d\n",
            no_flag[noocv_v[noc_f+0]]+1, // 1-indexed
            no_flag[noocv_v[noc_f+1]]+1, no_flag[noocv_v[noc_f+2]]+1, no_flag[noocv_v[noc_f+3]]+1,
            no_flag[noocv_v[noc_f+4]]+1, no_flag[noocv_v[noc_f+5]]+1, no_flag[noocv_v[noc_f+6]]+1,
            no_flag[noocv_v[noc_f+7]]+1);
        break;
      case 6:
        // prism...
        fprintf(fp,
            "%d %d %d %d %d %d %d %d\n",
            no_flag[noocv_v[noc_f+0]]+1, // 1-indexed
            no_flag[noocv_v[noc_f+1]]+1, no_flag[noocv_v[noc_f+2]]+1, no_flag[noocv_v[noc_f+2]]+1,
            no_flag[noocv_v[noc_f+3]]+1, no_flag[noocv_v[noc_f+4]]+1, no_flag[noocv_v[noc_f+5]]+1,
            no_flag[noocv_v[noc_f+5]]+1);
        break;
      case 5:
        // pyramid...
        fprintf(fp,
            "%d %d %d %d %d %d %d %d\n",
            no_flag[noocv_v[noc_f+0]]+1, // 1-indexed
            no_flag[noocv_v[noc_f+1]]+1, no_flag[noocv_v[noc_f+2]]+1, no_flag[noocv_v[noc_f+3]]+1,
            no_flag[noocv_v[noc_f+4]]+1, no_flag[noocv_v[noc_f+4]]+1, no_flag[noocv_v[noc_f+4]]+1,
            no_flag[noocv_v[noc_f+4]]+1);
        break;
      case 4:
        // tet...
        fprintf(fp,
            "%d %d %d %d %d %d %d %d\n",
            no_flag[noocv_v[noc_f+0]]+1, // 1-indexed
            no_flag[noocv_v[noc_f+1]]+1, no_flag[noocv_v[noc_f+2]]+1, no_flag[noocv_v[noc_f+2]]+1,
            no_flag[noocv_v[noc_f+3]]+1, no_flag[noocv_v[noc_f+3]]+1, no_flag[noocv_v[noc_f+3]]+1,
            no_flag[noocv_v[noc_f+3]]+1);
        break;
      default:
        cerr<<"Error: expecting a standard primative: "<<nnoc<<endl;
        throw(-1);
      }
    }
    else if (cv_flag[icv]>=0)
    {
      // this cv is getting dumped by tesselating to each face...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nnof = nof_l-nof_f+1;
        if (fa_flag[ifa]==-2)
        {
          // this face must be 3 of 4 nodes, thus forming a pyramid
          switch (nnof)
          {
          case 3:
            // joint to the cv forming a tet...
            if (cvofa[ifa][0]==icv)
            {
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[noofa_v[nof_f+2]]+1, // 1-indexed
                  no_flag[noofa_v[nof_f+1]]+1, no_flag[noofa_v[nof_f+0]]+1, no_flag[noofa_v[nof_f+0]]+1,
                  cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[noofa_v[nof_f+0]]+1, // 1-indexed
                  no_flag[noofa_v[nof_f+1]]+1, no_flag[noofa_v[nof_f+2]]+1, no_flag[noofa_v[nof_f+2]]+1,
                  cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1);
            }
            break;
          case 4:
            // joint to the cv forming a prism...
            if (cvofa[ifa][0]==icv)
            {
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[noofa_v[nof_f+3]]+1, // 1-indexed
                  no_flag[noofa_v[nof_f+2]]+1, no_flag[noofa_v[nof_f+1]]+1, no_flag[noofa_v[nof_f+0]]+1,
                  cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[noofa_v[nof_f+0]]+1, // 1-indexed
                  no_flag[noofa_v[nof_f+1]]+1, no_flag[noofa_v[nof_f+2]]+1, no_flag[noofa_v[nof_f+3]]+1,
                  cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1);
            }
            break;
          default:
            cerr<<"Error: expecting 3 or 4 nodes on this face: "<<nnof<<endl;
            throw(-1);
          }
        }
        else
        {
          assert(fa_flag[ifa]>=0);
          int ino1 = noofa_v[nof_l];
          for (int nof = nof_f; nof<=nof_l; ++nof)
          {
            int ino0 = ino1;
            ino1 = noofa_v[nof];
            // joint to the fa and cv forming a tet...
            if (cvofa[ifa][0]==icv)
            {
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[ino1]+1, // 1-indexed
                  no_flag[ino0]+1, fa_flag[ifa]+1, fa_flag[ifa]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1,
                  cv_flag[icv]+1);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                  no_flag[ino0]+1, // 1-indexed
                  no_flag[ino1]+1, fa_flag[ifa]+1, fa_flag[ifa]+1, cv_flag[icv]+1, cv_flag[icv]+1, cv_flag[icv]+1,
                  cv_flag[icv]+1);
            }
          }
        }
      }
    }
  }

  fclose(fp);

  if (mpi_rank<mpi_size-1)
  {
    int dummy = 1;
    MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
  }

  MPI_Barrier(mpi_comm);

}

void writeTecNames(FILE *fp, char *name)
{
  int i=0;
  do
  {
    int dummy = (int)name[i];
    fwrite( &dummy, sizeof(int), 1, fp);
    i++;
  }
  while(name[i] != '\0');

  int dummy = 0;
  fwrite( &dummy, sizeof(int), 1, fp);
}

void Ugp::writeFlaggedCvsTecplot(char * filename)
{

  // a sequential tecplot file writing...
  // write also flagged registered data

  if (mpi_rank==0) cout<<"writeFlaggedCvsTecplot: "<<filename<<endl;

  // step 1: build a global node numbering...
  // for nodes associated with flagged cv's...

  // build a global "node" numbering for the cv's, faces, nodes of flagged cv's...

  int my_buf[2];
  my_buf[0] = 0; // local ncount
  my_buf[1] = 0; // local ecount

  // use the connectivity to determine what has to be dumped...

  for (int ino = 0; ino<nno; ++ino)
    no_flag[ino] = -1;

  for (int ifa = 0; ifa<nfa; ++ifa)
    fa_flag[ifa] = -1;

  for (int icv = 0; icv<ncv; ++icv)
  {
    if (cv_flag[icv]==0)
    {
      // a negative 2 means skip completely...
      cv_flag[icv] = -2;
    }
    else
    {
      // at first assume we do not have to dump a value
      // at this cv centroid...
      cv_flag[icv] = -1;
      // all nodes have to be dumped...
      const int noc_f = noocv_i[icv];
      const int noc_l = noocv_i[icv+1]-1;
      const int nnoc = noc_l-noc_f+1;
      for (int noc = noc_f; noc<=noc_l; ++noc)
      {
        const int ino = noocv_v[noc];
        if (no_flag[ino]==-1)
        {
          no_flag[ino] = 0; // count later
        }
      }
      // faces may or may not be dumped...
      int ntri = 0;
      int nquad = 0;
      int ne = 0;
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      const int nfoc = foc_l-foc_f+1;
      for (int foc = foc_f; foc<=foc_l; ++foc)
      {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l-nof_f+1;
        if (fa_flag[ifa]==-1)
        {
          // this face has not yet been visited, so we need to check
          // if there are 3 or 4 nodes, and if 4, if they are coplanar...
          switch (nnof)
          {
          case 3:
            // a tri face - never a problem...
            fa_flag[ifa] = -2;
            break;
          case 4:
            // a quad face - only a problem if non-coplanar...
          {
            const int ino0 = noofa_v[nof_f];
            const int ino1 = noofa_v[nof_f+1];
            const int ino2 = noofa_v[nof_f+2];
            const int ino3 = noofa_v[nof_f+3];
            double v1[3] = { x_no[ino1][0]-x_no[ino0][0], x_no[ino1][1]-x_no[ino0][1], x_no[ino1][2]-x_no[ino0][2] };
            double v2[3] = { x_no[ino2][0]-x_no[ino0][0], x_no[ino2][1]-x_no[ino0][1], x_no[ino2][2]-x_no[ino0][2] };
            double normal[3] = { v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] };
            // check the distance of the 4th node off the plane of the other 3...
            double delta = (x_no[ino3][0]-x_no[ino0][0])*normal[0]+(x_no[ino3][1]-x_no[ino0][1])*normal[1]
                +(x_no[ino3][2]-x_no[ino0][2])*normal[2];
            if (fabs(delta)<0.01*pow(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2], 0.75))
            {
              // coplanar within 1%...
              fa_flag[ifa] = -2;
            }
            else
            {
              // not coplanar enough...
              //cout << "Warning: got non-coplanar quad face" << endl;
              fa_flag[ifa] = 0;
            }
          }
            break;
          default:
            // this is a face with one or more hanging nodes? need to dump the
            // center and triangulate to the edges...
            fa_flag[ifa] = 0;
            break;
          }
        }

        // now whether we visited it this time or not, we need to
        // what kind of face it is...

        if (fa_flag[ifa]==-2)
        {
          // tri or valid quad...
          if (nnof==3) ntri++;
          else
          {
            assert(nnof==4);
            nquad++;
          }
          // only one element here...
          ne += 1;
        }
        else
        {
          // if the face is zero, it is getting dumped, and we
          // need to add one element for each edge...
          assert(fa_flag[ifa]==0);
          ne += nnof;
        }

      }

      // now decide whether we dump the cv...
      if ((nnoc==8)&&(nfoc==6)&&(nquad==6)&&(ntri==0)&&(ne==6))
      {
        // hex cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==6)&&(nfoc==5)&&(nquad==3)&&(ntri==2)&&(ne==5))
      {
        // prism cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==5)&&(nfoc==5)&&(nquad==1)&&(ntri==4)&&(ne==5))
      {
        // pyramid cv...
        my_buf[1] += 1;
      }
      else if ((nnoc==4)&&(nfoc==4)&&(nquad==0)&&(ntri==4)&&(ne==4))
      {
        // tet cv...
        my_buf[1] += 1;
      }
      else
      {
        cv_flag[icv] = 0;
        my_buf[1] += ne;
      }
    }
  }

  // now my_buf[1] constains the element count.
  // cycle through nodes, faces, cells to get the vertex count...

  for (int ino = 0; ino<nno; ++ino)
  {
    if (no_flag[ino]==0)
    {
      no_flag[ino] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  for (int ifa = 0; ifa<nfa; ++ifa)
  {
    if (fa_flag[ifa]==0)
    {
      fa_flag[ifa] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  for (int icv = 0; icv<ncv; ++icv)
  {
    if (cv_flag[icv]==0)
    {
      cv_flag[icv] = my_buf[0];
      my_buf[0] += 1;
    }
  }

  // figure out our nodal displacement...

  int disp;
  MPI_Scan(&(my_buf[0]), &disp, 1, MPI_INT, MPI_SUM, mpi_comm);
  disp -= my_buf[0];

  for (int ino = 0; ino<nno; ++ino)
    if (no_flag[ino]>=0) no_flag[ino] += disp;

  for (int ifa = 0; ifa<nfa; ++ifa)
    if (fa_flag[ifa]>=0) fa_flag[ifa] += disp;

  for (int icv = 0; icv<ncv; ++icv)
    if (cv_flag[icv]>=0) cv_flag[icv] += disp;



  // send our counts to rank 0...
  int buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_INT, MPI_SUM, 0, mpi_comm);
  if (mpi_rank==0)
  {

    // rank 0 writes the header...

    // get number of variables to write first
    int numvars = 3; // x,y,z coordinates
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
      if (data->checkFlag())
        numvars++;

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
      if (data->checkFlag())
        numvars = numvars+3;

    FILE * fp;
    if ((fp = fopen(filename, "wb"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }


    char sdummy[300];
    sprintf(sdummy, "#!TDV75 ");
    fwrite(sdummy, sizeof(char), 8, fp);

    int idummy = 1;
    fwrite(&idummy, sizeof(int), 1, fp);          // byte order integer value of 1
    writeTecNames(fp, "flagged cvs");
    fwrite(&numvars, sizeof(int), 1, fp);
    writeTecNames(fp, "X");
    writeTecNames(fp, "Y");
    writeTecNames(fp, "Z");

    // write variable names
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
        case NO_DATA:
        case CV_DATA:
          // mixtrure fraction "Z" causes a problem with name correspondence...
          if ((strcmp(data->getName(), "X")==0)||(strcmp(data->getName(), "Y")==0)||(strcmp(data->getName(), "Z")==0))
          {
            char varname[200];
            sprintf(varname, "%s-scalar", data->getName());
            writeTecNames(fp, varname);
          }
          else
          {
            writeTecNames(fp, data->getName());
          }
          break;
        default:
          if (mpi_rank==0) cout<<"Warning: writeFlaggedCvsTecplot does not support data type: "<<data->getDatatype()
              <<", skipping."<<endl;
        }
      }
    }

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
        case NO_DATA:
        case CV_DATA:
          char varname[200];
          sprintf(varname, "%s-X", data->getName());
          writeTecNames(fp, varname);
          sprintf(varname, "%s-Y", data->getName());
          writeTecNames(fp, varname);
          sprintf(varname, "%s-Z", data->getName());
          writeTecNames(fp, varname);
          break;
        default:
          if (mpi_rank==0) cout<<"Warning: writeFlaggedCvsTecplot does not support data type: "<<data->getDatatype()
              <<", skipping."<<endl;
        }
      }
    }

    float zonemarker = 299.0;
    fwrite(&zonemarker, sizeof(float), 1, fp);
    sprintf(sdummy, "writeFlaggedCvsTecplot, np=%d", mpi_size);         // write zone name
    writeTecNames(fp, sdummy);
    idummy =  2;          fwrite(&idummy, sizeof(int), 1, fp);          // format = FEBLOCK
    idummy = -1;          fwrite(&idummy, sizeof(int), 1, fp);          // Zone color - tecplot choose
    int numPt = buf[0];   fwrite(&numPt, sizeof(int), 1, fp);           // number of nodes
    int numEl = buf[1];   fwrite(&numEl, sizeof(int), 1, fp);           // number of elements
    idummy = 3;           fwrite(&idummy, sizeof(int), 1, fp);          // Bricks

    float marker2  =  357.0;
    fwrite(&marker2,  sizeof(float), 1, fp);                            // end zone

    // starting data set here
    fwrite(&zonemarker, sizeof(float), 1, fp);                          // zone marker for data
    idummy = 0; fwrite(&idummy, sizeof(int), 1, fp);                    // repeat
    idummy = 1;                                                         // write format 1: float
    for (int ivar = 0; ivar < numvars; ++ivar)
      fwrite(&idummy, sizeof(int), 1, fp);


    fclose(fp);
  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // x,y,z...
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  for (int i = 0; i<3; i++)
  {

    int lines_written = 0;
    if (mpi_rank==0)
    {
      // on the first coordinate, start writing right away...
      if (i>0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
      }
    }
    else
    {
      // processes wait for previous to finish, and recv's lines_written...
      MPI_Status status;
      MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
    }

    FILE * fp;
    if ((fp = fopen(filename, "ab"))==NULL)
    {
      cerr<<"Error: cannot open file "<<filename<<endl;
      throw(-1);
    }

    // nodes are first...
    for (int ino = 0; ino<nno; ino++)
    {
      if (no_flag[ino]>=0)
      {
        float x = (float)x_no[ino][i];
        fwrite(&x, sizeof(float), 1, fp);
        if (++lines_written>5)
        {
          lines_written = 0;
        }
      }
    }

    // then faces...
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_flag[ifa]>=0)
      {
        double x_fa = 0.0;
        double weight = 0.0;
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++)
        {
          int ino = noofa_v[nof];
          x_fa += x_no[ino][i];
          weight += 1.0;
        }
        x_fa /= weight;
        float x = (float)x_fa;
        fwrite(&x, sizeof(float), 1, fp);
        if (++lines_written>5)
        {
          lines_written = 0;
        }
      }
    }

    // then cv's...
    for (int icv = 0; icv<ncv; icv++)
    {
      if (cv_flag[icv]>=0)
      {
        double x_cv = 0.0;
        double weight = 0.0;
        int noc_f = noocv_i[icv];
        int noc_l = noocv_i[icv+1]-1;
        for (int noc = noc_f; noc<=noc_l; noc++)
        {
          int ino = noocv_v[noc];
          x_cv += x_no[ino][i];
          weight += 1.0;
        }
        x_cv /= weight;
        float fdummy = (float)x_cv;
        fwrite(&fdummy, sizeof(float), 1, fp);
        if (++lines_written>5)
        {
          lines_written = 0;
        }
      }

    }

    if (mpi_rank<mpi_size-1)
    {
      // for ranks that are not last, just close and send a message on...
      fclose(fp);
      MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
    }
    else
    {
      // the last rank sends a message to the first...
      fclose(fp);
      int dummy_int = 0;
      MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
    }

  }

  // used for cv data interpolation...
  double * no_scalar = NULL;
  double * no_weight = NULL;
  double * no_scalar_ptr = NULL;

  // now the registered data (scalars)...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
      case NO_DATA:
        // ======================== node data ==========================
        no_scalar_ptr = (*data->ptr);
        break;
      case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_scalar==NULL) no_scalar = new double[nno];
        if (no_weight==NULL) no_weight = new double[nno];
        for (int ino = 0; ino<nno; ino++)
        {
          no_scalar[ino] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv<ncv; icv++)
        {
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            no_scalar[ino] += (*data->ptr)[icv];
            no_weight[ino] += 1.0;
          }
        }
        updateNoR1(no_scalar, ADD_DATA);
        updateNoR1(no_weight, ADD_DATA);
        for (int ino = 0; ino<nno; ino++)
          no_scalar[ino] /= no_weight[ino];
        no_scalar_ptr = no_scalar;
        break;
      default:
        continue;
      }

      // if we made it here, then we have some node data pointed to in no_scalar_ptr,
      // and can dump using no, fa, then cv writing...

      int lines_written = 0;
      if (mpi_rank==0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
      }
      else
      {
        // processes wait for previous to finish, and recv's lines_written...
        MPI_Status status;
        MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
      }

      FILE * fp;
      if ((fp = fopen(filename, "ab"))==NULL)
      {
        cerr<<"Error: cannot open file "<<filename<<endl;
        throw(-1);
      }

      // nodes are first...
      for (int ino = 0; ino<nno; ino++)
      {
        if (no_flag[ino]>=0)
        {
          float fdummy = (float)no_scalar_ptr[ino];
          fwrite(&fdummy, sizeof(float), 1, fp);
          if (++lines_written>5)
          {
            lines_written = 0;
          }
        }
      }

      // then faces...
      for (int ifa = 0; ifa<nfa; ifa++)
      {
        if (fa_flag[ifa]>=0)
        {
          double phi_fa = 0.0;
          double weight = 0.0;
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof<=nof_l; nof++)
          {
            int ino = noofa_v[nof];
            phi_fa += no_scalar_ptr[ino];
            weight += 1.0;
          }
          phi_fa /= weight;
          float fdummy = (float)phi_fa;
          fwrite(&fdummy, sizeof(float), 1, fp);
          if (++lines_written>5)
          {
            lines_written = 0;
          }
        }
      }

      // then cv's...
      for (int icv = 0; icv<ncv; icv++)
      {
        if (cv_flag[icv]>=0)
        {
          double phi_cv = 0.0;
          double weight = 0.0;
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            phi_cv += no_scalar_ptr[ino];
            weight += 1.0;
          }
          phi_cv /= weight;
          float fdummy = (float)phi_cv;
          fwrite(&fdummy, sizeof(float), 1, fp);
          if (++lines_written>5)
          {
            lines_written = 0;
          }
        }
      }

      if (mpi_rank<mpi_size-1)
      {
        // for ranks that are not last, just close and send a message on...
        fclose(fp);
        MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
      }
      else
      {
        // the last rank sends a message to the first...
        fclose(fp);
        int dummy_int = 0;
        MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
      }

    }

  }

  // cleanup and prepart for vector interpolation...
  if (no_scalar!=NULL) delete[] no_scalar;

  double (*no_vector)[3] = NULL;
  double (*no_vector_ptr)[3] = NULL;

  // loop on vector data...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
      case NO_DATA:
        // ================= nodal data ======================
        no_vector_ptr = (*data->ptr);
        break;
      case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_vector==NULL) no_vector = new double[nno][3];
        if (no_weight==NULL) no_weight = new double[nno];
        for (int ino = 0; ino<nno; ino++)
        {
          for (int i = 0; i<3; i++)
            no_vector[ino][i] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv<ncv; icv++)
        {
          int noc_f = noocv_i[icv];
          int noc_l = noocv_i[icv+1]-1;
          for (int noc = noc_f; noc<=noc_l; noc++)
          {
            int ino = noocv_v[noc];
            for (int i = 0; i<3; i++)
              no_vector[ino][i] += (*data->ptr)[icv][i];
            no_weight[ino] += 1.0;
          }
        }
        updateNoR2(no_vector, ADD_ROTATE_DATA);
        updateNoR1(no_weight, ADD_DATA);
        for (int ino = 0; ino<nno; ino++)
          for (int i = 0; i<3; i++)
            no_vector[ino][i] /= no_weight[ino];
        no_vector_ptr = no_vector;
        break;
      default:
        continue;
      }

      // write the data...
      for (int i = 0; i<3; i++)
      {

        int lines_written = 0;
        if (mpi_rank==0)
        {
          MPI_Status status;
          int dummy_int;
          MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
        }
        else
        {
          // processes wait for previous to finish, and recv's lines_written...
          MPI_Status status;
          MPI_Recv(&lines_written, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
        }

        FILE * fp;
        if ((fp = fopen(filename, "ab"))==NULL)
        {
          cerr<<"Error: cannot open file "<<filename<<endl;
          throw(-1);
        }

        // nodes are first...
        for (int ino = 0; ino<nno; ino++)
        {
          if (no_flag[ino]>=0)
          {
            float fdummy = (float)no_vector_ptr[ino][i];
            fwrite(&fdummy, sizeof(float), 1, fp);
            if (++lines_written>5)
            {
              lines_written = 0;
            }
          }
        }

        // then faces...
        for (int ifa = 0; ifa<nfa; ifa++)
        {
          if (fa_flag[ifa]>=0)
          {
            double phi_fa = 0.0;
            double weight = 0.0;
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof<=nof_l; nof++)
            {
              int ino = noofa_v[nof];
              phi_fa += no_vector_ptr[ino][i];
              weight += 1.0;
            }
            phi_fa /= weight;
            float fdummy = (float)phi_fa;
            fwrite(&fdummy, sizeof(float), 1, fp);
            if (++lines_written>5)
            {
              lines_written = 0;
            }
          }
        }

        // then cv's...
        for (int icv = 0; icv<ncv; icv++)
        {
          if (cv_flag[icv]>=0)
          {
            double phi_cv = 0.0;
            double weight = 0.0;
            int noc_f = noocv_i[icv];
            int noc_l = noocv_i[icv+1]-1;
            for (int noc = noc_f; noc<=noc_l; noc++)
            {
              int ino = noocv_v[noc];
              phi_cv += no_vector_ptr[ino][i];
              weight += 1.0;
            }
            phi_cv /= weight;
            float fdummy = (float)phi_cv;
            fwrite(&fdummy, sizeof(float), 1, fp);
            if (++lines_written>5)
            {
              lines_written = 0;
            }
          }
        }

        if (mpi_rank<mpi_size-1)
        {
          // for ranks that are not last, just close and send a message on...
          fclose(fp);
          MPI_Send(&lines_written, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
        }
        else
        {
          // the last rank sends a message to the first...
          fclose(fp);
          int dummy_int = 0;
          MPI_Send(&dummy_int, 1, MPI_INT, 0, 1111, mpi_comm);
        }
      }
    }
  }

  if (no_vector!=NULL) delete[] no_vector;
  if (no_weight!=NULL) delete[] no_weight;

  // ======================================
  // connectivity...
  // ======================================

  if (mpi_rank==0)
  {
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int, 1, MPI_INT, mpi_size-1, 1111, mpi_comm, &status);
  }
  else
  {
    // processes wait for previous to finish, and recv's lines_written...
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int, 1, MPI_INT, mpi_rank-1, 2222, mpi_comm, &status);
  }

  FILE * fp;
  if ((fp = fopen(filename, "ab"))==NULL)
  {
    cerr<<"Error: cannot open file "<<filename<<endl;
    throw(-1);
  }


  int idummy[8];
  if (mpi_rank==0)
  {
    idummy[0] = 0;
    fwrite(idummy, sizeof(int), 1, fp);   // repeat
  }


  for (int icv = 0; icv<ncv; icv++)
  {
    if (cv_flag[icv]==-1)
    {
      // this is a cv that is a standard primative. We can figure out
      // its type from its node count...
      int noc_f = noocv_i[icv];
      int noc_l = noocv_i[icv+1]-1;
      int nnoc = noc_l-noc_f+1;

      switch (nnoc)
      {
      case 8:
        // hex...
        idummy[0] = no_flag[noocv_v[noc_f+0]]+1;
        idummy[1] = no_flag[noocv_v[noc_f+1]]+1;
        idummy[2] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[3] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[4] = no_flag[noocv_v[noc_f+4]]+1;
        idummy[5] = no_flag[noocv_v[noc_f+5]]+1;
        idummy[6] = no_flag[noocv_v[noc_f+6]]+1;
        idummy[7] = no_flag[noocv_v[noc_f+7]]+1;
        fwrite(idummy, sizeof(int), 8, fp);
        break;
      case 6:
        // prism...
        idummy[0] = no_flag[noocv_v[noc_f+0]]+1;
        idummy[1] = no_flag[noocv_v[noc_f+1]]+1;
        idummy[2] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[3] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[4] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[5] = no_flag[noocv_v[noc_f+4]]+1;
        idummy[6] = no_flag[noocv_v[noc_f+5]]+1;
        idummy[7] = no_flag[noocv_v[noc_f+5]]+1;
        fwrite(idummy, sizeof(int), 8, fp);
        break;
      case 5:
        // pyramid...
        idummy[0] = no_flag[noocv_v[noc_f+0]]+1;
        idummy[1] = no_flag[noocv_v[noc_f+1]]+1;
        idummy[2] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[3] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[4] = no_flag[noocv_v[noc_f+4]]+1;
        idummy[5] = no_flag[noocv_v[noc_f+4]]+1;
        idummy[6] = no_flag[noocv_v[noc_f+4]]+1;
        idummy[7] = no_flag[noocv_v[noc_f+4]]+1;
        fwrite(idummy, sizeof(int), 8, fp);
        break;
      case 4:
        // tet...
        idummy[0] = no_flag[noocv_v[noc_f+0]]+1;
        idummy[1] = no_flag[noocv_v[noc_f+1]]+1;
        idummy[2] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[3] = no_flag[noocv_v[noc_f+2]]+1;
        idummy[4] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[5] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[6] = no_flag[noocv_v[noc_f+3]]+1;
        idummy[7] = no_flag[noocv_v[noc_f+3]]+1;
        fwrite(idummy, sizeof(int), 8, fp);
        break;
      default:
        cerr<<"Error: expecting a standard primative: "<<nnoc<<endl;
        throw(-1);
      }
    }
    else if (cv_flag[icv]>=0)
    {
      // this cv is getting dumped by tesselating to each face...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nnof = nof_l-nof_f+1;
        if (fa_flag[ifa]==-2)
        {
          // this face must be 3 of 4 nodes, thus forming a pyramid
          switch (nnof)
          {
          case 3:
            // joint to the cv forming a tet...
            if (cvofa[ifa][0]==icv)
            {
              idummy[0] = no_flag[noofa_v[nof_f+2]]+1;
              idummy[1] = no_flag[noofa_v[nof_f+1]]+1;
              idummy[2] = no_flag[noofa_v[nof_f+0]]+1;
              idummy[3] = no_flag[noofa_v[nof_f+0]]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              idummy[0] = no_flag[noofa_v[nof_f+0]]+1;
              idummy[1] = no_flag[noofa_v[nof_f+1]]+1;
              idummy[2] = no_flag[noofa_v[nof_f+2]]+1;
              idummy[3] = no_flag[noofa_v[nof_f+2]]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
            break;
          case 4:
            // joint to the cv forming a prism...
            if (cvofa[ifa][0]==icv)
            {
              idummy[0] = no_flag[noofa_v[nof_f+3]]+1;
              idummy[1] = no_flag[noofa_v[nof_f+2]]+1;
              idummy[2] = no_flag[noofa_v[nof_f+1]]+1;
              idummy[3] = no_flag[noofa_v[nof_f+0]]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              idummy[0] = no_flag[noofa_v[nof_f+0]]+1;
              idummy[1] = no_flag[noofa_v[nof_f+1]]+1;
              idummy[2] = no_flag[noofa_v[nof_f+2]]+1;
              idummy[3] = no_flag[noofa_v[nof_f+3]]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
            break;
          default:
            cerr<<"Error: expecting 3 or 4 nodes on this face: "<<nnof<<endl;
            throw(-1);
          }
        }
        else
        {
          assert(fa_flag[ifa]>=0);
          int ino1 = noofa_v[nof_l];
          for (int nof = nof_f; nof<=nof_l; ++nof)
          {
            int ino0 = ino1;
            ino1 = noofa_v[nof];
            // joint to the fa and cv forming a tet...
            if (cvofa[ifa][0]==icv)
            {
              idummy[0] = no_flag[ino1]+1;
              idummy[1] = no_flag[ino0]+1;
              idummy[2] = fa_flag[ifa]+1;
              idummy[3] = fa_flag[ifa]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
            else
            {
              assert(cvofa[ifa][1]==icv);
              idummy[0] = no_flag[ino0]+1;
              idummy[1] = no_flag[ino1]+1;
              idummy[2] = fa_flag[ifa]+1;
              idummy[3] = fa_flag[ifa]+1;
              idummy[4] = cv_flag[icv]+1;
              idummy[5] = cv_flag[icv]+1;
              idummy[6] = cv_flag[icv]+1;
              idummy[7] = cv_flag[icv]+1;
              fwrite(idummy, sizeof(int), 8, fp);
            }
          }
        }
      }
    }
  }

  fclose(fp);

  if (mpi_rank<mpi_size-1)
  {
    int dummy = 1;
    MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 2222, mpi_comm);
  }

  MPI_Barrier(mpi_comm);
}






#ifdef JUNKJUNK

void Ugp::writeFlaggedCvsTecplot(char * filename)
{

  // a sequential tecplot file writing...
  // write also flagged registered data

  if (mpi_rank == 0)
  cout << "writeFlaggedCvsTecplot: " << filename << endl;

  // step 1: build a global node numbering...
  // for nodes associated with flagged cv's...

  // build a global "node" numbering for the cv's, faces, nodes of flagged cv's...

  int my_buf[2];
  my_buf[0] = 0; // ncount
  my_buf[1] = 0; // ecount

  // use the connectivity to determine what has to be dumped...

  for (int ino = 0; ino < nno; ino++)
  no_flag[ino] = 0;

  for (int ifa = 0; ifa < nfa; ifa++)
  fa_flag[ifa] = 0;

  for (int icv = 0; icv < ncv; icv++)
  {
    if (cv_flag[icv] != 0)
    {
      cv_flag[icv] = 1;
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        if (fa_flag[ifa] == 0)
        fa_flag[ifa] = 1;
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof <= nof_l; nof++)
        {
          int ino = noofa_v[nof];
          if (no_flag[ino] == 0)
          no_flag[ino] = 1;
          ++(my_buf[1]);
        }
      }
    }
  }

  // number the nodes, faces and then cells - this makes
  // dumping a bit easier...

  for (int ino = 0; ino < nno; ino++)
  {
    if (no_flag[ino] == 1)
    {
      no_flag[ino] = ++(my_buf[0]); // 1-indexed to allow sign flipping
    }
  }
  for (int ifa = 0; ifa < nfa; ifa++)
  {
    if (fa_flag[ifa] == 1)
    {
      fa_flag[ifa] = ++(my_buf[0]);
    }
  }
  for (int icv = 0; icv < ncv; icv++)
  {
    if (cv_flag[icv] == 1)
    {
      cv_flag[icv] = ++(my_buf[0]);
    }
  }

  // figure out our nodal displacement...
  int disp;
  MPI_Scan(&(my_buf[0]),&disp,1,MPI_INT,MPI_SUM,mpi_comm);
  disp -= my_buf[0];

  for (int ino = 0; ino < nno; ino++)
  if (no_flag[ino] > 0)
  no_flag[ino] += disp;

  for (int ifa = 0; ifa < nfa; ifa++)
  if (fa_flag[ifa] > 0)
  fa_flag[ifa] += disp;

  for (int icv = 0; icv < ncv; icv++)
  if (cv_flag[icv] > 0)
  cv_flag[icv] += disp;

  // send our counts to rank 0...
  int buf[2];
  MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
  {

    // rank 0 writes the header...

    FILE * fp;
    if ( (fp=fopen(filename,"w"))==NULL )
    {
      cerr << "Error: cannot open file " << filename << endl;
      throw(-1);
    }
    fprintf(fp,"TITLE = \"flagged cvs\"\n");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    // add header for flagged registered data...
    for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
          case NO_DATA:
          case CV_DATA:
          // mixtrure fraction "Z" causes a problem with name correspondence...
          if ( (strcmp(data->getName(),"X") == 0) ||
              (strcmp(data->getName(),"Y") == 0) ||
              (strcmp(data->getName(),"Z") == 0) )
          {
            fprintf(fp,"\"%s (scalar)\"\n",data->getName());
          }
          else
          {
            fprintf(fp,"\"%s\"\n",data->getName());
          }
          break;
          default:
          if (mpi_rank == 0)
          cout << "Warning: writeFlaggedCvsTecplot does not support data type: " << data->getDatatype() <<
          ", skipping." << endl;
        }
      }
    }

    for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++)
    {
      if (data->checkFlag())
      {
        switch (data->getDatatype())
        {
          case NO_DATA:
          case CV_DATA:
          fprintf(fp,"\"%s-X\"\n",data->getName());
          fprintf(fp,"\"%s-Y\"\n",data->getName());
          fprintf(fp,"\"%s-Z\"\n",data->getName());
          break;
          default:
          if (mpi_rank == 0)
          cout << "Warning: writeFlaggedCvsTecplot does not support data type: " << data->getDatatype() <<
          ", skipping." << endl;
        }
      }
    }

    fprintf(fp,"ZONE T=\"writeFlaggedCvsTecplot, np=%d\"\n",mpi_size);
    fprintf(fp,"N=%d, E=%d, F=FEBLOCK, ET=TETRAHEDRON\n",buf[0],buf[1]);
    fclose(fp);

  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // x,y,z...
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  for (int i = 0; i < 3; i++)
  {

    int lines_written = 0;
    if (mpi_rank == 0)
    {
      // on the first coordinate, start writing right away...
      if (i > 0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int,1,MPI_INT,mpi_size-1,1111,mpi_comm,&status);
      }
    }
    else
    {
      // processes wait for previous to finish, and recv's lines_written...
      MPI_Status status;
      MPI_Recv(&lines_written,1,MPI_INT,mpi_rank-1,2222,mpi_comm,&status);
    }

    FILE * fp;
    if ( (fp=fopen(filename,"a"))==NULL )
    {
      cerr << "Error: cannot open file " << filename << endl;
      throw(-1);
    }

    // nodes are first...
    for (int ino = 0; ino < nno; ino++)
    {
      if (no_flag[ino] > 0)
      {
        fprintf(fp,"%le ",x_no[ino][i]);
        if (++lines_written > 5)
        {
          fprintf(fp,"\n");
          lines_written=0;
        }
      }
    }

    // then faces...
    for (int ifa = 0; ifa < nfa; ifa++)
    {
      if (fa_flag[ifa] > 0)
      {
        double x_fa = 0.0;
        double weight = 0.0;
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof <= nof_l; nof++)
        {
          int ino = noofa_v[nof];
          x_fa += x_no[ino][i];
          weight += 1.0;
        }
        x_fa /= weight;
        fprintf(fp,"%le ",x_fa);
        if (++lines_written > 5)
        {
          fprintf(fp,"\n");
          lines_written=0;
        }
      }
    }

    // then cv's...
    for (int icv = 0; icv < ncv; icv++)
    {
      if (cv_flag[icv] > 0)
      {
        double x_cv = 0.0;
        double weight = 0.0;
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv+1]-1;
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];
          assert( fa_flag[ifa] != 0 );
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof <= nof_l; nof++)
          {
            int ino = noofa_v[nof];
            assert( no_flag[ino] != 0 );
            // use a flipping ot the face node to include each node only once
            // in the cv average...
            if (no_flag[ino] > 0)
            {
              // here we can just flip the sign because it is never zero - 
              no_flag[ino] = -no_flag[ino];
              x_cv += x_no[ino][i];
              weight += 1.0;
            }
          }
        }
        // flip all the nodes back...
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof <= nof_l; nof++)
          {
            int ino = noofa_v[nof];
            if (no_flag[ino] < 0)
            {
              no_flag[ino] = -no_flag[ino];
            }
          }
        }
        // write the cv...
        x_cv /= weight;
        fprintf(fp,"%le ",x_cv);
        if (++lines_written > 5)
        {
          fprintf(fp,"\n");
          lines_written=0;
        }
      }
    }

    if (mpi_rank < mpi_size-1)
    {
      // for ranks that are not last, just close and send a message on...
      fclose(fp);
      MPI_Send(&lines_written,1,MPI_INT,mpi_rank+1,2222,mpi_comm);
    }
    else
    {
      // the last rank writes a return if required, then sends a message to the first...
      if (lines_written != 0)
      fprintf(fp,"\n");
      fclose(fp);
      int dummy_int = 0;
      MPI_Send(&dummy_int,1,MPI_INT,0,1111,mpi_comm);
    }

  }

  // used for cv data interpolation...
  double * no_scalar = NULL;
  double * no_weight = NULL;
  double * no_scalar_ptr = NULL;

  // now the registered data (scalars)...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
        case NO_DATA:
        // ======================== node data ==========================
        no_scalar_ptr = (*data->ptr);
        break;
        case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_scalar == NULL) no_scalar = new double[nno];
        if (no_weight == NULL) no_weight = new double[nno];
        for (int ino = 0; ino < nno; ino++)
        {
          no_scalar[ino] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv < ncv; icv++)
        {
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv+1]-1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              // include unselected nodes in this data set...
              if (no_flag[ino] >= 0)
              {
                no_flag[ino] = -no_flag[ino]-1;
                no_scalar[ino] += (*data->ptr)[icv];
                no_weight[ino] += 1.0;
              }
            }
          }
          // flip all the nodes back...
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              if (no_flag[ino] < 0)
              {
                no_flag[ino] = -no_flag[ino]-1;
              }
            }
          }
        }
        updateNoR1(no_scalar,ADD_DATA);
        updateNoR1(no_weight,ADD_DATA);
        for (int ino = 0; ino < nno; ino++)
        no_scalar[ino] /= no_weight[ino];
        no_scalar_ptr = no_scalar;
        break;
        default:
        continue;
      }

      // if we made it here, then we have some node data pointed to in no_scalar_ptr,
      // and can dump using no, fa, then cv writing...

      int lines_written = 0;
      if (mpi_rank == 0)
      {
        MPI_Status status;
        int dummy_int;
        MPI_Recv(&dummy_int,1,MPI_INT,mpi_size-1,1111,mpi_comm,&status);
      }
      else
      {
        // processes wait for previous to finish, and recv's lines_written...
        MPI_Status status;
        MPI_Recv(&lines_written,1,MPI_INT,mpi_rank-1,2222,mpi_comm,&status);
      }

      FILE * fp;
      if ( (fp=fopen(filename,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << filename << endl;
        throw(-1);
      }

      // nodes are first...
      for (int ino = 0; ino < nno; ino++)
      {
        if (no_flag[ino] > 0)
        {
          fprintf(fp,"%le ",no_scalar_ptr[ino]);
          if (++lines_written > 5)
          {
            fprintf(fp,"\n");
            lines_written=0;
          }
        }
      }

      // then faces...
      for (int ifa = 0; ifa < nfa; ifa++)
      {
        if (fa_flag[ifa] > 0)
        {
          double phi_fa = 0.0;
          double weight = 0.0;
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          for (int nof = nof_f; nof <= nof_l; nof++)
          {
            int ino = noofa_v[nof];
            phi_fa += no_scalar_ptr[ino];
            weight += 1.0;
          }
          phi_fa /= weight;
          fprintf(fp,"%le ",phi_fa);
          if (++lines_written > 5)
          {
            fprintf(fp,"\n");
            lines_written=0;
          }
        }
      }

      // then cv's...
      for (int icv = 0; icv < ncv; icv++)
      {
        if (cv_flag[icv] > 0)
        {
          double phi_cv = 0.0;
          double weight = 0.0;
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv+1]-1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            assert( fa_flag[ifa] != 0 );
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              assert( no_flag[ino] != 0 );
              // use a flipping ot the face node to include each node only once
              // in the cv average...
              if (no_flag[ino] > 0)
              {
                // here we can just flip the sign because it is never zero - 
                no_flag[ino] = -no_flag[ino];
                phi_cv += no_scalar_ptr[ino];
                weight += 1.0;
              }
            }
          }
          // flip all the nodes back...
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              if (no_flag[ino] < 0)
              {
                no_flag[ino] = -no_flag[ino];
              }
            }
          }
          // write the cv...
          phi_cv /= weight;
          fprintf(fp,"%le ",phi_cv);
          if (++lines_written > 5)
          {
            fprintf(fp,"\n");
            lines_written=0;
          }
        }
      }

      if (mpi_rank < mpi_size-1)
      {
        // for ranks that are not last, just close and send a message on...
        fclose(fp);
        MPI_Send(&lines_written,1,MPI_INT,mpi_rank+1,2222,mpi_comm);
      }
      else
      {
        // the last rank writes a return if required, then sends a message to the first...
        if (lines_written != 0)
        fprintf(fp,"\n");
        fclose(fp);
        int dummy_int = 0;
        MPI_Send(&dummy_int,1,MPI_INT,0,1111,mpi_comm);
      }

    }

  }

  // cleanup and prepart for vector interpolation...
  if (no_scalar != NULL) delete[] no_scalar;

  double (*no_vector)[3] = NULL;
  double (*no_vector_ptr)[3] = NULL;

  // loop on vector data...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++)
  {
    if (data->checkFlag())
    {

      switch (data->getDatatype())
      {
        case NO_DATA:
        // ================= nodal data ======================
        no_vector_ptr = (*data->ptr);
        break;
        case CV_DATA:
        // ======================== cv data ==========================
        // cv data is first averaged to the nodes, then dumped
        // using the nodal data algorithm...
        // include ALL cv data and ALL nodes at this point...
        if (no_vector == NULL) no_vector = new double[nno][3];
        if (no_weight == NULL) no_weight = new double[nno];
        for (int ino = 0; ino < nno; ino++)
        {
          for (int i = 0; i < 3; i++)
          no_vector[ino][i] = 0.0;
          no_weight[ino] = 0.0;
        }
        // distribute cv data to the nodes...
        for (int icv = 0; icv < ncv; icv++)
        {
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv+1]-1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              // include unselected nodes in this data set...
              if (no_flag[ino] >= 0)
              {
                no_flag[ino] = -no_flag[ino]-1;
                for (int i = 0; i < 3; i++)
                no_vector[ino][i] += (*data->ptr)[icv][i];
                no_weight[ino] += 1.0;
              }
            }
          }
          // flip all the nodes back...
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              if (no_flag[ino] < 0)
              {
                no_flag[ino] = -no_flag[ino]-1;
              }
            }
          }
        }
        updateNoR2(no_vector,ADD_ROTATE_DATA);
        updateNoR1(no_weight,ADD_DATA);
        for (int ino = 0; ino < nno; ino++)
        for (int i = 0; i < 3; i++)
        no_vector[ino][i] /= no_weight[ino];
        no_vector_ptr = no_vector;
        break;
        default:
        continue;
      }

      // write the data...
      for (int i = 0; i < 3; i++)
      {

        int lines_written = 0;
        if (mpi_rank == 0)
        {
          MPI_Status status;
          int dummy_int;
          MPI_Recv(&dummy_int,1,MPI_INT,mpi_size-1,1111,mpi_comm,&status);
        }
        else
        {
          // processes wait for previous to finish, and recv's lines_written...
          MPI_Status status;
          MPI_Recv(&lines_written,1,MPI_INT,mpi_rank-1,2222,mpi_comm,&status);
        }

        FILE * fp;
        if ( (fp=fopen(filename,"a"))==NULL )
        {
          cerr << "Error: cannot open file " << filename << endl;
          throw(-1);
        }

        // nodes are first...
        for (int ino = 0; ino < nno; ino++)
        {
          if (no_flag[ino] > 0)
          {
            fprintf(fp,"%le ",no_vector_ptr[ino][i]);
            if (++lines_written > 5)
            {
              fprintf(fp,"\n");
              lines_written=0;
            }
          }
        }

        // then faces...
        for (int ifa = 0; ifa < nfa; ifa++)
        {
          if (fa_flag[ifa] > 0)
          {
            double phi_fa = 0.0;
            double weight = 0.0;
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            for (int nof = nof_f; nof <= nof_l; nof++)
            {
              int ino = noofa_v[nof];
              phi_fa += no_vector_ptr[ino][i];
              weight += 1.0;
            }
            phi_fa /= weight;
            fprintf(fp,"%le ",phi_fa);
            if (++lines_written > 5)
            {
              fprintf(fp,"\n");
              lines_written=0;
            }
          }
        }

        // then cv's...
        for (int icv = 0; icv < ncv; icv++)
        {
          if (cv_flag[icv] > 0)
          {
            double phi_cv = 0.0;
            double weight = 0.0;
            int foc_f = faocv_i[icv];
            int foc_l = faocv_i[icv+1]-1;
            for (int foc = foc_f; foc <= foc_l; foc++)
            {
              int ifa = faocv_v[foc];
              assert( fa_flag[ifa] != 0 );
              int nof_f = noofa_i[ifa];
              int nof_l = noofa_i[ifa+1]-1;
              for (int nof = nof_f; nof <= nof_l; nof++)
              {
                int ino = noofa_v[nof];
                assert( no_flag[ino] != 0 );
                // use a flipping ot the face node to include each node only once
                // in the cv average...
                if (no_flag[ino] > 0)
                {
                  // here we can just flip the sign because it is never zero - 
                  no_flag[ino] = -no_flag[ino];
                  phi_cv += no_vector_ptr[ino][i];
                  weight += 1.0;
                }
              }
            }
            // flip all the nodes back...
            for (int foc = foc_f; foc <= foc_l; foc++)
            {
              int ifa = faocv_v[foc];
              int nof_f = noofa_i[ifa];
              int nof_l = noofa_i[ifa+1]-1;
              for (int nof = nof_f; nof <= nof_l; nof++)
              {
                int ino = noofa_v[nof];
                if (no_flag[ino] < 0)
                {
                  no_flag[ino] = -no_flag[ino];
                }
              }
            }
            // write the cv...
            phi_cv /= weight;
            fprintf(fp,"%le ",phi_cv);
            if (++lines_written > 5)
            {
              fprintf(fp,"\n");
              lines_written=0;
            }
          }
        }

        if (mpi_rank < mpi_size-1)
        {
          // for ranks that are not last, just close and send a message on...
          fclose(fp);
          MPI_Send(&lines_written,1,MPI_INT,mpi_rank+1,2222,mpi_comm);
        }
        else
        {
          // the last rank writes a return if required, then sends a message to the first...
          if (lines_written != 0)
          fprintf(fp,"\n");
          fclose(fp);
          int dummy_int = 0;
          MPI_Send(&dummy_int,1,MPI_INT,0,1111,mpi_comm);
        }
      }
    }
  }

  if (no_vector != NULL) delete[] no_vector;
  if (no_weight != NULL) delete[] no_weight;

  // ======================================
  // connectivity...
  // ======================================

  if (mpi_rank == 0)
  {
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int,1,MPI_INT,mpi_size-1,1111,mpi_comm,&status);
  }
  else
  {
    // processes wait for previous to finish, and recv's lines_written...
    MPI_Status status;
    int dummy_int;
    MPI_Recv(&dummy_int,1,MPI_INT,mpi_rank-1,2222,mpi_comm,&status);
  }

  FILE * fp;
  if ( (fp=fopen(filename,"a"))==NULL )
  {
    cerr << "Error: cannot open file " << filename << endl;
    throw(-1);
  }

  for (int icv = 0; icv < ncv; icv++)
  {
    if (cv_flag[icv] > 0)
    {
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; foc++)
      {
        int ifa = faocv_v[foc];
        assert( fa_flag[ifa] > 0 );
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        if (cvofa[ifa][0] == icv)
        {
          int ino2 = noofa_v[nof_l];
          for (int nof = nof_f; nof <= nof_l; nof++)
          {
            int ino1 = ino2;
            ino2 = noofa_v[nof];
            assert( no_flag[ino1] > 0 );
            assert( no_flag[ino2] > 0 );
            fprintf(fp,"%d %d %d %d\n",fa_flag[ifa],no_flag[ino2],no_flag[ino1],cv_flag[icv]);
          }
        }
        else
        {
          assert( cvofa[ifa][1] == icv );
          int ino2 = noofa_v[nof_f];
          for (int nof = nof_l; nof >= nof_f; nof--)
          {
            int ino1 = ino2;
            ino2 = noofa_v[nof];
            assert( no_flag[ino1] > 0 );
            assert( no_flag[ino2] > 0 );
            fprintf(fp,"%d %d %d %d\n",fa_flag[ifa],no_flag[ino2],no_flag[ino1],cv_flag[icv]);
          }
        }
      }
    }
  }

  fclose(fp);

  if ( mpi_rank < mpi_size-1 )
  {
    int dummy = 1;
    MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,2222,mpi_comm);
  }

  MPI_Barrier(mpi_comm);

}

#endif

// ================================================================
// Parallel Fluent writer
// ================================================================

#define FLUENT_TYPE_INTERIOR 2
#define FLUENT_TYPE_WALL 3

void Ugp::writeFluentMsh(char * filename)
{

  if (mpi_rank==0) cout<<"writeFluentMsh: "<<filename<<endl;

  // build global node, face, and cv indices and (on the first time only)
  // the mpi types required to write data...

  // what to do about this?...
  //assert( io_write_flag == 0 );

  // -------------------
  // nodes...
  // -------------------

  int my_no_offset;
  MPI_Scan(&nno, &my_no_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  my_no_offset -= nno;
  for (int ino = 0; ino<nno; ino++)
    no_flag[ino] = ino+my_no_offset;
  updateNoI1(no_flag, MIN_NO_PERIODIC_DATA);

  int my_no_count = 0;
  for (int ino = 0; ino<nno; ino++)
  {
    if (no_flag[ino]==ino+my_no_offset)
    {
      my_no_count++;
      no_flag[ino] = 1; // this one gets dumped
    }
    else
    {
      no_flag[ino] = 0;
    }
  }

  MPI_Scan(&my_no_count, &my_no_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  int no_count = my_no_offset;
  MPI_Bcast(&no_count, 1, MPI_INT, mpi_size-1, mpi_comm);
  my_no_offset -= my_no_count;

  int * ino_global = new int[nno];
  for (int ino = 0; ino<nno; ino++)
  {
    if (no_flag[ino]==1)
    {
      ino_global[ino] = my_no_offset++;
    }
    else
    {
      ino_global[ino] = -1;
    }
  }
  updateNoI1(ino_global, MAX_NO_PERIODIC_DATA);

  // check...
  for (int ino = 0; ino<nno; ino++)
    assert((ino_global[ino]>=0)&&(ino_global[ino]<no_count));

  // -------------------
  // faces...
  // -------------------

  // this stuff required below, so keep outside of braces...

  int my_fa_offset;
  MPI_Scan(&nfa, &my_fa_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  my_fa_offset -= nfa;
  for (int ifa = 0; ifa<nfa; ifa++)
    fa_flag[ifa] = ifa+my_fa_offset;
  updateFaI1(fa_flag, MIN_NO_PERIODIC_DATA);

  int my_fa_count = 0;
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    if (fa_flag[ifa]==ifa+my_fa_offset)
    {
      my_fa_count++;
      fa_flag[ifa] = 1; // this one gets dumped
    }
    else
    {
      fa_flag[ifa] = 0;
    }
  }

  MPI_Scan(&my_fa_count, &my_fa_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  int fa_count = my_fa_offset;
  MPI_Bcast(&fa_count, 1, MPI_INT, mpi_size-1, mpi_comm);

  int nzones = faZoneList.size();
  int * my_fa_zone_count = new int[nzones];

  int izone = -1;
  for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone!=faZoneList.end(); faZone++)
  {
    izone++;
    my_fa_zone_count[izone] = 0;
    switch (faZone->getKind())
    {
    case FA_ZONE_INTERNAL:
      // internal zones have 2 ranges...
      for (int ifa = faZone->ifa_i_f; ifa<=faZone->ifa_i_l; ++ifa)
      {
        if (fa_flag[ifa]==1)
        {
          my_fa_zone_count[izone] += 1;
        }
      }
      // no break;
    default:
      // all other zones have just one range... 
      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ++ifa)
      {
        if (fa_flag[ifa]==1)
        {
          my_fa_zone_count[izone] += 1;
        }
      }
    }
  }

  // all need fa_zone_count...
  int * fa_zone_count = new int[nzones];
  MPI_Allreduce(my_fa_zone_count, fa_zone_count, nzones, MPI_INT, MPI_SUM, mpi_comm);
  delete[] my_fa_zone_count;

  // all check...
  int fa_count_check = 0;
  for (int i = 0; i<nzones; ++i)
    fa_count_check += fa_zone_count[i];
  assert(fa_count_check==fa_count);

  // -------------------
  // cvs...
  // -------------------

  int my_cv_count = ncv;
  int my_cv_offset;
  MPI_Scan(&my_cv_count, &my_cv_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
  int cv_count = my_cv_offset;
  MPI_Bcast(&cv_count, 1, MPI_INT, mpi_size-1, mpi_comm);
  my_cv_offset -= my_cv_count;
  for (int icv = 0; icv<ncv; icv++)
    cv_flag[icv] = icv+my_cv_offset;

  int * cvofa_tmp = new int[nfa_bpi];
  for (int ifa = 0; ifa<nfa_b; ++ifa)
    cvofa_tmp[ifa] = -1;
  for (int ifa = nfa_b; ifa<nfa_bpi; ++ifa)
    cvofa_tmp[ifa] = cv_flag[cvofa[ifa][0]];
  updateFaI1(cvofa_tmp, REPLACE_DATA);

  // --------------------------------------------------
  // now ready to write the file...
  // --------------------------------------------------

  // timing...

  MPI_Barrier(mpi_comm);
  double wtime;
  if (mpi_rank==0) wtime = MPI_Wtime();

  // open...

  if (MPI_File_open(mpi_comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: writeFluentMsh: could not open : "<<filename<<endl;
    throw(-1);
  }

  // this normally works - if you need it bigger, I would be interested to know
  // F. Ham: fham@stanford.edu
  MPI_Offset cbuf_size = (MPI_Offset) 64*max(my_fa_count, my_no_count);
  char * cbuf = new char[cbuf_size];

  int n = sprintf(cbuf, "(0 \"Extrude to Fluent\")\n");
  n += sprintf(cbuf+n, "\n");
  n += sprintf(cbuf+n, "(0 \"Dimensions:\")\n");
  n += sprintf(cbuf+n, "(2 3)\n");
  n += sprintf(cbuf+n, "\n");

  //
  // header and global dimensions
  //
  n += sprintf(cbuf+n, "(0 \"Grid:\")\n");
  n += sprintf(cbuf+n, "\n");
  n += sprintf(cbuf+n, "(12 (0 %x %x 0))\n", 1, cv_count);
  n += sprintf(cbuf+n, "(13 (0 %x %x 0))\n", 1, fa_count);
  n += sprintf(cbuf+n, "(10 (0 %x %x 0 3))\n", 1, no_count);
  n += sprintf(cbuf+n, "\n");

  //
  // cells
  //
  int cv_element_type = 7; // polyhedra = 7
  int cv_type = 1; // active = 1 - inactive = 32
  n += sprintf(cbuf+n, "(0 \"Cells:\")\n");
  n += sprintf(cbuf+n, "(12 (%x %x %x %x %x))\n", 4001, 1, cv_count, cv_type, cv_element_type);
  n += sprintf(cbuf+n, "(45 (%d fluid %s)())\n", 4001, "fluid");
  n += sprintf(cbuf+n, "\n");

  //
  // start of faces
  //
  n += sprintf(cbuf+n, "(0 \"Faces:\")\n");
  n += sprintf(cbuf+n, "\n");

  if (mpi_rank==0)
  {
    cout<<" > header"<<endl;
    MPI_Status status;
    MPI_File_write(fh, cbuf, n, MPI_CHAR, &status);
  }

  // initialize the offset...
  MPI_Offset offset = n;

  // now the face zones...

  //
  // faces
  //
  int fa_element_type = 5; // 5 is polygonal; 4 is quadrilateral
  int this_zone_type;
  char this_zone_name[128];
  char this_fluent_name[128];
  int this_zone_id = 10;
  int ifa_global = 0;
  izone = -1;
  for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone!=faZoneList.end(); faZone++)
  {
    if (mpi_rank==0) cout<<" > faces for zone \""<<faZone->getName()<<"\""<<endl;
    izone++;
    if ((faZone->getKind()==FA_ZONE_BOUNDARY)||faZone->isPeriodic())
    {
      sprintf(this_zone_name, "%s", faZone->getName());
      sprintf(this_fluent_name, "wall");
      this_zone_type = FLUENT_TYPE_WALL;
      this_zone_id += 1;
      if (mpi_rank==0)
      {
        n = sprintf(cbuf, "(13 (%x %x %x %x %x)(\n", this_zone_id, ifa_global+1, ifa_global+fa_zone_count[izone],
            this_zone_type, fa_element_type);
      }
      else
      {
        n = 0;
      }
      ifa_global += fa_zone_count[izone];
      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ++ifa)
      {
        if (fa_flag[ifa]==1)
        {
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          // number of nodes...
          n += sprintf(cbuf+n, "%x ", nof_l-nof_f+1);
          // reverse the looping order (msh is left-handed)...
          for (int nof = nof_l; nof>=nof_f; nof--)
          {
            int ino = ino_global[noofa_v[nof]];
            n += sprintf(cbuf+n, "%x ", ino+1);
          }
          int icv0 = cvofa[ifa][0];
          assert((icv0>=0)&&(icv0<ncv));
          int icv1 = cvofa[ifa][1];
          assert(icv1==-1);
          n += sprintf(cbuf+n, "%x %x\n", cv_flag[icv0]+1, icv1+1);
          assert(n<(cbuf_size*4)/5);
        }
      }
      // last rank closes...
      if (mpi_rank==mpi_size-1)
      {
        n += sprintf(cbuf+n, "))\n");
      }
    }
    else if (faZone->getKind()==FA_ZONE_INTERNAL)
    {
      sprintf(this_zone_name, "%s", faZone->getName());
      sprintf(this_fluent_name, "interior");
      this_zone_type = FLUENT_TYPE_INTERIOR;
      this_zone_id += 1;
      if (mpi_rank==0)
      {
        n = sprintf(cbuf, "(13 (%x %x %x %x %x)(\n", this_zone_id, ifa_global+1, ifa_global+fa_zone_count[izone],
            this_zone_type, fa_element_type);
      }
      else
      {
        n = 0;
      }
      ifa_global += fa_zone_count[izone];
      for (int ifa = faZone->ifa_f; ifa<=faZone->ifa_l; ++ifa)
      {
        if (fa_flag[ifa]==1)
        {
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          // number of nodes...
          n += sprintf(cbuf+n, "%x ", nof_l-nof_f+1);
          // reverse the looping order (msh is left-handed)...
          for (int nof = nof_l; nof>=nof_f; nof--)
          {
            int ino = ino_global[noofa_v[nof]];
            n += sprintf(cbuf+n, "%x ", ino+1);
          }
          int icv0 = cvofa[ifa][0];
          assert((icv0>=0)&&(icv0<ncv));
          assert((ifa>=nfa_b)&&(ifa<nfa_bpi));
          int icv1 = cvofa_tmp[ifa];
          n += sprintf(cbuf+n, "%x %x\n", cv_flag[icv0]+1, icv1+1);
          assert(n<(cbuf_size*4)/5);
        }
      }
      for (int ifa = faZone->ifa_i_f; ifa<=faZone->ifa_i_l; ++ifa)
      {
        if (fa_flag[ifa]==1)
        {
          int nof_f = noofa_i[ifa];
          int nof_l = noofa_i[ifa+1]-1;
          // number of nodes...
          n += sprintf(cbuf+n, "%x ", nof_l-nof_f+1);
          // reverse the looping order (msh is left-handed)...
          for (int nof = nof_l; nof>=nof_f; nof--)
          {
            int ino = ino_global[noofa_v[nof]];
            n += sprintf(cbuf+n, "%x ", ino+1);
          }
          int icv0 = cvofa[ifa][0];
          assert((icv0>=0)&&(icv0<ncv));
          int icv1 = cvofa[ifa][1];
          assert((icv1>=0)&&(icv1<ncv));
          n += sprintf(cbuf+n, "%x %x\n", cv_flag[icv0]+1, cv_flag[icv1]+1);
          assert(n<(cbuf_size*4)/5);
        }
      }
      // last rank closes...
      if (mpi_rank==mpi_size-1)
      {
        n += sprintf(cbuf+n, "))\n");
      }
    }
    else
    {
      cerr<<"Error: unsupported zone for now."<<endl;
      throw(-1);
    }

    // apparently this is 4 on BGL. This can cause problems for
    // grids that are big because of the char buffering 
    // technique used.

    assert(sizeof(MPI_Aint)==8);

    MPI_Aint nA = (MPI_Aint) n;
    MPI_Aint dispA;
    MPI_Scan(&nA, &dispA, 1, MPI_OFFSET_DATATYPE, MPI_SUM, mpi_comm);
    MPI_Aint global_sizeA = dispA;
    MPI_Bcast(&global_sizeA, 1, MPI_OFFSET_DATATYPE, mpi_size-1, mpi_comm);
    dispA -= nA;

    MPI_Datatype char_type;
    MPI_Type_hindexed(1, &n, &dispA, MPI_CHAR, &char_type);
    MPI_Type_commit(&char_type);

    if (MPI_File_set_view(fh, offset, MPI_CHAR, char_type, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: mpi_char view at offset: "<<offset<<endl;
      throw(-1);
    }

    MPI_Status status;
    MPI_File_write_all(fh, cbuf, n, MPI_CHAR, &status);

    offset += global_sizeA;
    MPI_Type_free(&char_type);

    /*
     int disp;
     MPI_Scan(&n,&disp,1,MPI_INT,MPI_SUM,mpi_comm);
     int global_size = disp;
     MPI_Bcast(&global_size,1,MPI_INT,mpi_size-1,mpi_comm);
     if (mpi_rank == 0)
     cout << " > global_size: " << global_size << endl;
     disp -= n;

     MPI_Pause("OK?");

     MPI_Datatype char_type;
     MPI_Type_indexed(1,&n,&disp,MPI_CHAR,&char_type);
     MPI_Type_commit(&char_type);

     if (MPI_File_set_view(fh, offset, MPI_CHAR, char_type, "native", MPI_INFO_NULL) != 0) {
     cerr << "Error: mpi_char view at offset: " << offset << endl;
     throw(-1);
     }

     MPI_Status status;
     MPI_File_write_all(fh, cbuf, n, MPI_CHAR, &status);

     offset += global_size;
     MPI_Type_free(&char_type);
     */

  }

  // loop a second time to dump names...

  n = 0;
  this_zone_id = 10;
  izone = -1;
  for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone!=faZoneList.end(); faZone++)
  {
    izone++;
    if ((faZone->getKind()==FA_ZONE_BOUNDARY)||faZone->isPeriodic())
    {
      sprintf(this_zone_name, "%s", faZone->getName());
      sprintf(this_fluent_name, "wall");
      this_zone_type = FLUENT_TYPE_WALL;
      this_zone_id += 1;
      n += sprintf(cbuf+n, "(39 (%d %s %s 1)())\n", this_zone_id, this_fluent_name, this_zone_name);
    }
    else if (faZone->getKind()==FA_ZONE_INTERNAL)
    {
      sprintf(this_zone_name, "%s", faZone->getName());
      sprintf(this_fluent_name, "interior");
      this_zone_type = FLUENT_TYPE_INTERIOR;
      this_zone_id += 1;
      n += sprintf(cbuf+n, "(39 (%d %s %s 1)())\n", this_zone_id, this_fluent_name, this_zone_name);
    }
    else
    {
      cerr<<"Error: unsupported zone for now."<<endl;
      throw(-1);
    }
  }

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  if (mpi_rank==0)
  {
    cout<<" > face zone names"<<endl;
    MPI_Status status;
    MPI_File_write(fh, cbuf, n, MPI_CHAR, &status);
  }

  offset += n;

  if (mpi_rank==0) cout<<" > nodes"<<endl;

  //
  // nodes...
  //
  if (mpi_rank==0)
  {
    n = sprintf(cbuf, "(0 \"Nodes:\")\n");
    n += sprintf(cbuf+n, "(10 (1 %x %x 1 3 )(\n", 1, no_count);
  }
  else
  {
    n = 0;
  }

  for (int ino = 0; ino<nno; ino++)
  {
    if (no_flag[ino]==1)
    {
      n += sprintf(cbuf+n, "%15.8e %15.8e %15.8e\n", x_no[ino][0], x_no[ino][1], x_no[ino][2]);
      assert(n<(cbuf_size*4)/5);
    }
  }

  if (mpi_rank==mpi_size-1)
  {
    n += sprintf(cbuf+n, "))\n");
  }

  {

    assert(sizeof(MPI_Aint)==8);

    MPI_Aint nA = (MPI_Aint) n;
    MPI_Aint dispA;
    MPI_Scan(&nA, &dispA, 1, MPI_OFFSET_DATATYPE, MPI_SUM, mpi_comm);
    MPI_Aint global_sizeA = dispA;
    MPI_Bcast(&global_sizeA, 1, MPI_OFFSET_DATATYPE, mpi_size-1, mpi_comm);
    dispA -= nA;

    MPI_Datatype char_type;
    MPI_Type_hindexed(1, &n, &dispA, MPI_CHAR, &char_type);
    MPI_Type_commit(&char_type);

    if (MPI_File_set_view(fh, offset, MPI_CHAR, char_type, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: mpi_char view at offset: "<<offset<<endl;
      throw(-1);
    }

    MPI_Status status;
    MPI_File_write_all(fh, cbuf, n, MPI_CHAR, &status);

    offset += global_sizeA;
    MPI_Type_free(&char_type);

    /*
     int disp;
     MPI_Scan(&n,&disp,1,MPI_INT,MPI_SUM,mpi_comm);
     int global_size = disp;
     MPI_Bcast(&global_size,1,MPI_INT,mpi_size-1,mpi_comm);
     disp -= n;

     MPI_Datatype char_type;
     MPI_Type_indexed(1,&n,&disp,MPI_CHAR,&char_type);
     MPI_Type_commit(&char_type);

     if (MPI_File_set_view(fh, offset, MPI_CHAR, char_type, "native", MPI_INFO_NULL) != 0) {
     cerr << "Error: mpi_char view at offset: " << offset << endl;
     throw(-1);
     }

     MPI_Status status;
     MPI_File_write_all(fh, cbuf, n, MPI_CHAR, &status);

     offset += global_size;
     MPI_Type_free(&char_type);
     */

  }

  // trim the file to the current size and close
  MPI_File_set_size(fh, offset);
  MPI_File_close(&fh);

  if (mpi_rank==0) cout<<" > done"<<endl;

  // timing...
  MPI_Barrier(mpi_comm);
  if (mpi_rank==0)
  {
    wtime = MPI_Wtime()-wtime;
    cout<<" > writeFluentMsh timing: nno: "<<no_count<<", nfa: "<<fa_count<<", ncv: "<<cv_count<<", np: "<<mpi_size
        <<", size[B]: "<<offset<<", time[s]: "<<wtime<<", rate[MB/s]: "<<(double) offset/1.0E+6/wtime<<endl;
  }

  // cleanup...  
  delete[] cbuf;
  delete[] cvofa_tmp;
  delete[] fa_zone_count;
  delete[] ino_global;

}

// ######################################################################
// io routines
// ######################################################################

#define UGP_IO_MAGIC_NUMBER           123581321
#define UGP_IO_VERSION                1

// restart record indexing - do not change these...

#define UGP_IO_NO_FA_CV_NOOFA_COUNTS  10
#define UGP_IO_I0                     11
#define UGP_IO_D0                     12

#define UGP_IO_FA_CHECK               20
#define UGP_IO_NOOFA_I_AND_V          21
#define UGP_IO_CVOFA                  22
#define UGP_IO_FA_ZONE                23
#define UGP_IO_FA_D1                  24
#define UGP_IO_FA_D3                  25

#define UGP_IO_NO_CHECK               30
#define UGP_IO_X_NO                   31
#define UGP_IO_NO_D1                  32
#define UGP_IO_NO_D3                  33

#define UGP_IO_CV_CHECK               40
#define UGP_IO_CV_PART                41
#define UGP_IO_CV_D1                  42
#define UGP_IO_CV_D3                  43
#define UGP_IO_CV_D33                 44

#define UGP_IO_DATA                   50
#define UGP_IO_EOF                    51

void Ugp::initIoCommon()
{

  assert(io_common_flag==0);

  // check expected sizes...
  if ((sizeof(int)!=4)||(sizeof(MPI_Offset)!=8)||(sizeof(float)!=4)||(sizeof(double)!=8))
  {
    if (mpi_rank==0)
    {
      cerr<<"Error: unexpected sizes:"<<endl;
      cerr<<"int       : "<<sizeof(int)<<endl;
      cerr<<"MPI_Offset: "<<sizeof(MPI_Offset)<<endl;
      cerr<<"float     : "<<sizeof(float)<<endl;
      cerr<<"double    : "<<sizeof(double)<<endl;
    }
    throw(-1);
  }

  MPI_Datatype oldtypes[5];
  int blockcounts[5];
  MPI_Aint offsets[5];

  // the name...
  offsets[0] = 0;
  oldtypes[0] = MPI_CHAR;
  blockcounts[0] = UGP_IO_HEADER_NAME_LEN;

  // the id...
  offsets[1] = offsets[0]+1*UGP_IO_HEADER_NAME_LEN;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 1;

  // the offset skip...
  offsets[2] = offsets[1]+4;
  oldtypes[2] = MPI_OFFSET_DATATYPE;
  blockcounts[2] = 1;

  // the 16 ints...
  offsets[3] = offsets[2]+8;
  oldtypes[3] = MPI_INT;
  blockcounts[3] = 16;

  // the 16 reals...
  offsets[4] = offsets[3]+4*16;
  oldtypes[4] = MPI_DOUBLE;
  blockcounts[4] = 16;

  // build the header type...
  MPI_Type_struct(5, blockcounts, offsets, oldtypes, &MPI_Header);
  MPI_Type_commit(&MPI_Header);

  int header_size_int;
  MPI_Aint header_extent;
  MPI_Type_extent(MPI_Header, &header_extent);
  MPI_Type_size(MPI_Header, &header_size_int);

  // check that is is 256
  if ((header_extent!=header_size_int)||(header_size_int!=256))
  {
    cerr<<"Error: header size/extent problem: "<<endl;
    cout<<"header_extent : "<<header_extent<<endl;
    cout<<"header_size  : "<<header_size_int<<endl;

    int long_size;
    MPI_Type_size(MPI_LONG, &long_size);
    cout<<"long_size: "<<long_size<<endl;

    throw(-1);
  }

  // these are type MPI_Offset (8 byte int)...
  int_size = 4;
  double_size = 8;
  header_size = 256;

}

// #######################################################

void Ugp::writeRestart(const char * filename)
{

  if (mpi_rank==0) cout<<"writeRestart: "<<filename<<endl;

  // make sure the various io types are initialized...

  if (io_common_flag==0)
  {
    assert(io_write_flag==0);
    initIoCommon();
    io_common_flag = 1;
  }
  if (io_read_flag==1)
  {
    assert(io_write_flag==0);
    cleanupReadRestart();
    io_read_flag = 2; // NOT back to 0 - only allow one readRestart per unique Ugp for now
  }

  // build global node, face, and cv indices and (on the first time only)
  // the mpi types required to write data...

  // -------------------
  // nodes...
  // -------------------

  {
    // stripe the nodes in 2 zones - boundary, then the rest...
    int * no_zone = new int[nno];
    for (int ino = 0; ino<nno; ino++)
      no_zone[ino] = 1;
    // reset boundary nodes to zone 0...
    for (int ifa = 0; ifa<nfa_b; ifa++)
    {
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof<=nof_l; nof++)
      {
        int ino = noofa_v[nof];
        no_zone[ino] = 0;
      }
    }
    // boundary nodes (0's) win...
    updateNoI1(no_zone, MIN_DATA);

    // no_zone now contains 0 or 1.
    // now use no_flag to break inter-processor ties...

    int my_offset;
    MPI_Scan(&nno, &my_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
    my_offset -= nno;
    for (int ino = 0; ino<nno; ino++)
      no_flag[ino] = ino+my_offset;
    updateNoI1(no_flag, MIN_NO_PERIODIC_DATA);

    int n_no_zones = 2;
    int * my_no_zone_count = new int[n_no_zones];
    for (int i = 0; i<n_no_zones; i++)
      my_no_zone_count[i] = 0;

    for (int ino = 0; ino<nno; ino++)
    {
      if (no_flag[ino]==ino+my_offset)
      {
        no_flag[ino] = my_no_zone_count[no_zone[ino]];
        my_no_zone_count[no_zone[ino]] += 1;
      }
      else
      {
        no_flag[ino] = -1;
        no_zone[ino] = -1;
      }
    }

    // no_zone now contains 0,1 or -1 in nodes not to be dumped.
    // and no_flag contains a local node indexing. 

    // first, build a node list for local packing...

    my_no_count = 0; // this is part of Ugp
    int * my_no_zone_disp = new int[n_no_zones];
    for (int i = 0; i<n_no_zones; i++)
    {
      my_no_zone_disp[i] = my_no_count;
      my_no_count += my_no_zone_count[i];
    }

    if (io_write_flag==0)
    {
      assert(my_no_list==NULL);
      my_no_list = new int[my_no_count];
    }

    for (int ino = 0; ino<nno; ino++)
    {
      if (no_zone[ino]>=0)
      {
        my_no_list[my_no_zone_disp[no_zone[ino]]] = ino;
        my_no_zone_disp[no_zone[ino]] += 1;
      }
    }

    // Now build stuff neccessary for mpi types and turn 
    // no_flag into a global node indexing with the
    // correct offsets...

    MPI_Scan(my_no_zone_count, my_no_zone_disp, n_no_zones, MPI_INT, MPI_SUM, mpi_comm);

    int * no_zone_count = new int[n_no_zones];
    if (mpi_rank==mpi_size-1)
    {
      for (int i = 0; i<n_no_zones; i++)
      {
        no_zone_count[i] = my_no_zone_disp[i];
      }
    }
    MPI_Bcast(no_zone_count, n_no_zones, MPI_INT, mpi_size-1, mpi_comm);

    no_count = 0; // also part of Ugp
    for (int i = 0; i<n_no_zones; i++)
    {
      // adjust the disp so we can use it...
      my_no_zone_disp[i] = my_no_zone_disp[i]-my_no_zone_count[i]+no_count;
      // and count the global node sizes...
      no_count += no_zone_count[i];
    }

    // add disp to positive no_flag's to produce a global node numbering in no_flag...
    for (int ino = 0; ino<nno; ino++)
      if (no_zone[ino]>=0) no_flag[ino] += my_no_zone_disp[no_zone[ino]];
    updateNoI1(no_flag, MAX_NO_PERIODIC_DATA);

    // now all nodes contain their associated global index in no_flag. We can tell
    // the specific nodes that won the ties (and need to be dumped) by no_zone >= 0.

    //if (mpi_rank == 0)
    //  cout << "> nno: " << no_count << endl;

    // on the first time, build and commit the mpi types...

    if (io_write_flag==0)
    {

      MPI_Type_indexed(n_no_zones, my_no_zone_count, my_no_zone_disp, MPI_INT, &no_i1_type);
      MPI_Type_commit(&no_i1_type);

      MPI_Type_indexed(n_no_zones, my_no_zone_count, my_no_zone_disp, MPI_DOUBLE, &no_d1_type);
      MPI_Type_commit(&no_d1_type);

      for (int i = 0; i<n_no_zones; i++)
      {
        my_no_zone_count[i] *= 3;
        my_no_zone_disp[i] *= 3;
      }
      MPI_Type_indexed(n_no_zones, my_no_zone_count, my_no_zone_disp, MPI_DOUBLE, &no_d3_type);
      MPI_Type_commit(&no_d3_type);

    }

    delete[] no_zone;
    delete[] my_no_zone_count;
    delete[] my_no_zone_disp;
    delete[] no_zone_count;

  }

  // -------------------
  // faces...
  // -------------------

  // this stuff required below, so keep outside of braces...

  int * fa_zone_count = NULL;
  int n_fa_zones = 0;

  {
    int * fa_zone = new int[nfa];
    for (int ifa = 0; ifa<nfa; ifa++)
      fa_zone[ifa] = -1;

    // use the faZoneList to set the fa_zone...
    {
      int izone = -1;
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
      {
        izone++;
        switch (zone->getKind())
        {
        case FA_ZONE_INTERNAL:
          // internal faces have the additional face range...
          for (int ifa = zone->ifa_i_f; ifa<=zone->ifa_i_l; ifa++)
          {
            assert(fa_zone[ifa]==-1);
            fa_zone[ifa] = izone;
          }
          // no break...
        default:
          // all face zones have the normal range...
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            assert(fa_zone[ifa]==-1);
            fa_zone[ifa] = izone;
          }
        }
      }

      // check that all faces have been accounted for...
      for (int ifa = 0; ifa<nfa; ifa++)
        assert(fa_zone[ifa]>=0);
    }

    // the number of zones...
    n_fa_zones = faZoneList.size();

    int * my_fa_zone_count = new int[n_fa_zones];
    int * my_noofa_zone_count = new int[n_fa_zones];
    for (int i = 0; i<n_fa_zones; i++)
    {
      my_fa_zone_count[i] = 0;
      my_noofa_zone_count[i] = 0;
    }

    // now break ties...

    int my_offset;
    MPI_Scan(&nfa, &my_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
    my_offset -= nfa;
    for (int ifa = 0; ifa<nfa; ifa++)
      fa_flag[ifa] = ifa+my_offset;

    updateFaI1(fa_flag, MIN_NO_PERIODIC_DATA);

    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_flag[ifa]==ifa+my_offset)
      {
        fa_flag[ifa] = my_fa_zone_count[fa_zone[ifa]];
        my_fa_zone_count[fa_zone[ifa]] += 1;
        my_noofa_zone_count[fa_zone[ifa]] += noofa_i[ifa+1]-noofa_i[ifa];
      }
      else
      {
        fa_flag[ifa] = -1;
        fa_zone[ifa] = -1;
      }
    }

    // fist build the list for local packing...

    my_fa_count = 0; // this is part of Ugp
    int * my_fa_zone_disp = new int[n_fa_zones];
    for (int i = 0; i<n_fa_zones; i++)
    {
      my_fa_zone_disp[i] = my_fa_count;
      my_fa_count += my_fa_zone_count[i];
    }

    if (io_write_flag==0)
    {

      assert(my_fa_list==NULL);
      my_fa_list = new int[my_fa_count];

    }

    for (int ifa = 0; ifa<nfa; ifa++)
    {
      if (fa_zone[ifa]>=0)
      {
        my_fa_list[my_fa_zone_disp[fa_zone[ifa]]] = ifa;
        my_fa_zone_disp[fa_zone[ifa]] += 1;
      }
    }

    // now build the mpi struct...

    MPI_Scan(my_fa_zone_count, my_fa_zone_disp, n_fa_zones, MPI_INT, MPI_SUM, mpi_comm);

    int * my_noofa_zone_disp = new int[n_fa_zones];
    MPI_Scan(my_noofa_zone_count, my_noofa_zone_disp, n_fa_zones, MPI_INT, MPI_SUM, mpi_comm);

    fa_zone_count = new int[n_fa_zones]; // allocated outside face braces
    int * noofa_zone_count = new int[n_fa_zones];
    if (mpi_rank==mpi_size-1)
    {
      for (int i = 0; i<n_fa_zones; i++)
      {
        fa_zone_count[i] = my_fa_zone_disp[i];
        noofa_zone_count[i] = my_noofa_zone_disp[i];
      }
    }
    MPI_Bcast(fa_zone_count, n_fa_zones, MPI_INT, mpi_size-1, mpi_comm);
    MPI_Bcast(noofa_zone_count, n_fa_zones, MPI_INT, mpi_size-1, mpi_comm);

    fa_count = 0;
    my_noofa_count = 0;
    noofa_count = 0;
    for (int i = 0; i<n_fa_zones; i++)
    {
      // adjust the disp so we can use it...
      my_fa_zone_disp[i] = my_fa_zone_disp[i]-my_fa_zone_count[i]+fa_count;
      my_noofa_zone_disp[i] = my_noofa_zone_disp[i]-my_noofa_zone_count[i]+noofa_count;
      // and count the local and global face sizes...
      fa_count += fa_zone_count[i];
      my_noofa_count += my_noofa_zone_count[i];
      noofa_count += noofa_zone_count[i];
    }

    // now put a global index in fa_flag
    for (int ifa = 0; ifa<nfa; ifa++)
      if (fa_zone[ifa]>=0) fa_flag[ifa] += my_fa_zone_disp[fa_zone[ifa]];
    updateFaI1(fa_flag, MAX_NO_PERIODIC_DATA);

    //if (mpi_rank == 0)
    //  cout << "> nfa: " << fa_count << endl;

    if (io_write_flag==0)
    {

      MPI_Type_indexed(n_fa_zones, my_fa_zone_count, my_fa_zone_disp, MPI_INT, &fa_i1_type);
      MPI_Type_commit(&fa_i1_type);

      MPI_Type_indexed(n_fa_zones, my_fa_zone_count, my_fa_zone_disp, MPI_DOUBLE, &fa_d1_type);
      MPI_Type_commit(&fa_d1_type);

      for (int i = 0; i<n_fa_zones; i++)
      {
        my_fa_zone_count[i] *= 2;
        my_fa_zone_disp[i] *= 2;
      }
      MPI_Type_indexed(n_fa_zones, my_fa_zone_count, my_fa_zone_disp, MPI_INT, &fa_i2_type);
      MPI_Type_commit(&fa_i2_type);

      MPI_Type_indexed(n_fa_zones, my_noofa_zone_count, my_noofa_zone_disp, MPI_INT, &noofa_i1_type);
      MPI_Type_commit(&noofa_i1_type);

    }

    delete[] fa_zone;
    delete[] my_fa_zone_count;
    delete[] my_noofa_zone_count;
    delete[] my_fa_zone_disp;
    delete[] my_noofa_zone_disp;
    delete[] noofa_zone_count;

  }

  // -------------------
  // cvs - this is easier because they are not
  // duplicated across processors, but we still need
  // to allow for cv zones: for now, use 2...
  // -------------------

  {

    int n_cv_zones = 2;
    int * my_cv_zone_count = new int[n_cv_zones];
    for (int i = 0; i<n_cv_zones; i++)
      my_cv_zone_count[i] = 0;

    int * cv_zone = new int[ncv];
    for (int icv = 0; icv<ncv_b; icv++)
    {
      cv_zone[icv] = 0;
      cv_flag[icv] = my_cv_zone_count[0];
      my_cv_zone_count[0] += 1;
    }
    for (int icv = ncv_b; icv<ncv; icv++)
    {
      cv_zone[icv] = 1;
      cv_flag[icv] = my_cv_zone_count[1];
      my_cv_zone_count[1] += 1;
    }

    // build my_cv_list for cv packing...

    my_cv_count = 0; // this is part of Ugp
    int * my_cv_zone_disp = new int[n_cv_zones];
    for (int i = 0; i<n_cv_zones; i++)
    {
      my_cv_zone_disp[i] = my_cv_count;
      my_cv_count += my_cv_zone_count[i];
    }

    // we should nave counted all cvs...
    assert(my_cv_count==ncv);

    if (io_write_flag==0)
    {

      assert(my_cv_list==NULL);
      my_cv_list = new int[my_cv_count];

    }

    for (int icv = 0; icv<ncv; icv++)
    {
      if (cv_zone[icv]>=0)
      {
        my_cv_list[my_cv_zone_disp[cv_zone[icv]]] = icv;
        my_cv_zone_disp[cv_zone[icv]] += 1;
      }
      else
      {
        cerr<<"Error: for now all cvs should be packed"<<endl;
        throw(-1);
      }
    }

    // now the cv mpi types...

    MPI_Scan(my_cv_zone_count, my_cv_zone_disp, n_cv_zones, MPI_INT, MPI_SUM, mpi_comm);

    int * cv_zone_count = new int[n_cv_zones];
    if (mpi_rank==mpi_size-1)
    {
      for (int i = 0; i<n_cv_zones; i++)
      {
        cv_zone_count[i] = my_cv_zone_disp[i];
      }
    }
    MPI_Bcast(cv_zone_count, n_cv_zones, MPI_INT, mpi_size-1, mpi_comm);

    cv_count = 0;
    for (int i = 0; i<n_cv_zones; i++)
    {
      // adjust the disp so we can use it to build a type...
      my_cv_zone_disp[i] = my_cv_zone_disp[i]-my_cv_zone_count[i]+cv_count;
      // and count the local and global cv sizes...
      cv_count += cv_zone_count[i];
    }

    // now put a global index in cv_flag
    for (int icv = 0; icv<ncv; icv++)
    {
      assert(cv_zone[icv]>=0);
      cv_flag[icv] += my_cv_zone_disp[cv_zone[icv]];
    }

    //if (mpi_rank == 0)
    //  cout << "> ncv: " << cv_count << endl;

    if (io_write_flag==0)
    {

      MPI_Type_indexed(n_cv_zones, my_cv_zone_count, my_cv_zone_disp, MPI_INT, &cv_i1_type);
      MPI_Type_commit(&cv_i1_type);

      MPI_Type_indexed(n_cv_zones, my_cv_zone_count, my_cv_zone_disp, MPI_DOUBLE, &cv_d1_type);
      MPI_Type_commit(&cv_d1_type);

      for (int i = 0; i<n_cv_zones; i++)
      {
        my_cv_zone_count[i] *= 3;
        my_cv_zone_disp[i] *= 3;
      }
      MPI_Type_indexed(n_cv_zones, my_cv_zone_count, my_cv_zone_disp, MPI_DOUBLE, &cv_d3_type);
      MPI_Type_commit(&cv_d3_type);

      for (int i = 0; i<n_cv_zones; i++)
      {
        my_cv_zone_count[i] *= 3;
        my_cv_zone_disp[i] *= 3;
      }
      MPI_Type_indexed(n_cv_zones, my_cv_zone_count, my_cv_zone_disp, MPI_DOUBLE, &cv_d33_type);
      MPI_Type_commit(&cv_d33_type);

    }

    delete[] my_cv_zone_count;
    delete[] cv_zone;
    delete[] my_cv_zone_disp;
    delete[] cv_zone_count;

  }

  // --------------------------------------------------
  // now ready to write the file...
  // --------------------------------------------------

  io_write_flag = 1;

  // timing...

  MPI_Barrier(mpi_comm);
  double wtime;
  if (mpi_rank==0) wtime = MPI_Wtime();

  // open...

  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(mpi_comm, dummy, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: writeRestart: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // ------------------------------------------
  // the start of the file contains 2 ints:
  // 1. magic number
  // 2. io version...
  // ------------------------------------------

  if (mpi_rank==0)
  {

    cout<<" > header"<<endl;

    int ibuf[2];
    ibuf[0] = UGP_IO_MAGIC_NUMBER;
    ibuf[1] = UGP_IO_VERSION;

    MPI_Status status;
    MPI_File_write(fh, ibuf, 2, MPI_INT, &status);

  }

  offset = int_size*2;

  // ------------------------------------------
  // counts...
  // ------------------------------------------

  if (mpi_rank==0)
  {

    cout<<" > counts"<<endl;

    Header header;
    sprintf(header.name, "NO_FA_CV_NOOFA_COUNTS");
    header.id = UGP_IO_NO_FA_CV_NOOFA_COUNTS;
    header.skip = header_size;
    header.idata[0] = no_count;
    header.idata[1] = fa_count;
    header.idata[2] = cv_count;
    header.idata[3] = noofa_count;

    MPI_Status status;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }
  offset = offset+header_size;

  // ------------------------------------------
  // node check: write the global index currently in no_flag.
  // When we re-read this, we can check for 0,1,2,3,4...nno-1.
  // ------------------------------------------

  writeNodeI1(no_flag, UGP_IO_NO_CHECK, "NO_CHECK");

  // ------------------------------------------
  // face check...
  // ------------------------------------------

  writeFaceI1(fa_flag, UGP_IO_FA_CHECK, "FA_CHECK");

  // ------------------------------------------
  // cv check...
  // ------------------------------------------

  writeCvI1(cv_flag, UGP_IO_CV_CHECK, "CV_CHECK");

  // ------------------------------------------
  // connectivity
  // ------------------------------------------

  writeNoofa();
  writeCvofa();

  // ------------------------------------------
  // face zones...
  // ------------------------------------------

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the face zone headers...

  if (mpi_rank==0)
  {

    MPI_Status status;
    Header header;
    header.id = UGP_IO_FA_ZONE;
    header.skip = header_size;

    int izone = -1;
    int ifa_global = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone!=faZoneList.end(); zone++)
    {
      izone++;
      sprintf(header.name, zone->getName());
      // put the kind in idata[0]...
      header.idata[0] = zone->getKind();
      // put the 0-indexed face range in idata[1], idata[2]...
      header.idata[1] = ifa_global;
      ifa_global += fa_zone_count[izone];
      header.idata[2] = ifa_global-1;
      // put the periodic transform (if any) in rdata...
      zone->getPeriodicData(header.rdata);
      // and write...
      MPI_File_write(fh, &header, 1, MPI_Header, &status);
    }
    // make sure we came out alright...
    assert(izone+1==n_fa_zones);
    assert(ifa_global==fa_count);

  }

  // all update offset...

  offset += header_size*n_fa_zones;

  // ------------------------------------------
  // current cv partition...
  // ------------------------------------------

  writeCvPart();

  // ------------------------------------------
  // node coordinates: x_no...
  // ------------------------------------------

  writeNodeD3(x_no, UGP_IO_X_NO, "X_NO");

  // ------------------------------------------
  // start of data...
  // ------------------------------------------

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }
  if (mpi_rank==0)
  {
    MPI_Status status;
    Header header;
    header.id = UGP_IO_DATA;
    sprintf(header.name, "DATA");
    header.skip = header_size;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);
  }
  offset += header_size;

  // ------------------------------------------
  // data...
  // ------------------------------------------

  for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++)
  {
    writeI0(*data->ptr, UGP_IO_I0, data->getName());
  }

  for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++)
  {
    writeD0(*data->ptr, UGP_IO_D0, data->getName());
  }

  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
  {
    switch (data->getDatatype())
    {
    case NO_DATA:
      writeNodeD1(*data->ptr, UGP_IO_NO_D1, data->getName());
      break;

    case CV_DATA:
      writeCVsD1(*data->ptr, UGP_IO_CV_D1, data->getName());
      break;

    case FA_DATA:
      writeFaceD1(*data->ptr, UGP_IO_FA_D1, data->getName());
      break;

    default:
      if (mpi_rank==0) cout<<"Warning: Scalar io not setup for registered datatype: "<<data->getDatatype()<<endl;
    }
  }

  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
  {
    switch (data->getDatatype())
    {

    case NO_DATA:
      writeNodeD3(*data->ptr, UGP_IO_NO_D3, data->getName());
      break;

    case CV_DATA:
      writeCVsD3(*data->ptr, UGP_IO_CV_D3, data->getName());
      break;

    default:
      if (mpi_rank==0) cout<<"Warning: Vector io not setup for registered datatype: "<<data->getDatatype()<<endl;
    }
  }

  for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
  {
    switch (data->getDatatype())
    {

    case CV_DATA:
      writeCVsD33(*data->ptr, UGP_IO_CV_D33, data->getName());
      break;

    default:
      if (mpi_rank==0) cout<<"Warning: Tensor io not setup for registered datatype: "<<data->getDatatype()<<endl;
    }
  }

  // ------------------------------------------
  // eof...
  // ------------------------------------------

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }
  if (mpi_rank==0)
  {
    MPI_Status status;
    Header header;
    header.id = UGP_IO_EOF;
    sprintf(header.name, "EOF");
    header.skip = header_size;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);
  }
  offset += header_size;

  // trim the file to the current size and close
  MPI_File_set_size(fh, offset);
  MPI_File_close(&fh);

  // timing...
  MPI_Barrier(mpi_comm);
  if (mpi_rank==0)
  {
    wtime = MPI_Wtime()-wtime;
    cout<<" > writeRestart timing: nno: "<<no_count<<", nfa: "<<fa_count<<", ncv: "<<cv_count<<", np: "<<mpi_size
        <<", size[B]: "<<offset<<", time[s]: "<<wtime<<", rate[MB/s]: "<<(double) offset/1.0E+6/wtime<<endl;
  }

  // cleanup...  
  delete[] fa_zone_count;

  // note that we leave the io types, etc to
  // allow snapshots, etc, although we have now
  // set the io_write_flag to 1...

}

void Ugp::writeSnapshot(const int index)
{

  // if we don't have a consistent restart (connectivity, coordinates),
  // then the snapshot is useless, so...

  if (io_write_flag!=1) writeRestart(index);

  char filename[32];
  sprintf(filename, "snapshot.%06d.out", index);

  if (mpi_rank==0) cout<<"writeSnapshot: "<<filename<<endl;

  // timing...

  MPI_Barrier(mpi_comm);
  double wtime;
  if (mpi_rank==0) wtime = MPI_Wtime();

  // open...

  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(mpi_comm, dummy, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: writeSnapshot: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // ------------------------------------------
  // just as in a restart, the start of the file contains 2 ints:
  // 1. magic number
  // 2. io version
  // ------------------------------------------

  if (mpi_rank==0)
  {

    int ibuf[2];
    ibuf[0] = UGP_IO_MAGIC_NUMBER;
    ibuf[1] = UGP_IO_VERSION;

    MPI_Status status;
    MPI_File_write(fh, ibuf, 2, MPI_INT, &status);

  }

  offset = int_size*2;

  // ------------------------------------------
  // then straight to data...
  // ------------------------------------------

  for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++)
  {
    if (data->checkFlag())
    {
      if (mpi_rank==0) cout<<" > int value: "<<data->getName()<<": "<<*data->ptr<<endl;
      writeI0(*data->ptr, UGP_IO_I0, data->getName());
    }
  }

  for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++)
  {
    if (data->checkFlag())
    {
      if (mpi_rank==0) cout<<" > double value: "<<data->getName()<<": "<<*data->ptr<<endl;
      writeD0(*data->ptr, UGP_IO_D0, data->getName());
    }
  }

  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
  {
    if (data->checkFlag())
    {
      switch (data->getDatatype())
      {
      case NO_DATA:
        if (mpi_rank==0) cout<<" > node scalar: "<<data->getName()<<endl;
        writeNodeD1(*data->ptr, UGP_IO_NO_D1, data->getName());
        break;
      case CV_DATA:
        if (mpi_rank==0) cout<<" > cv scalar: "<<data->getName()<<endl;
        writeCVsD1(*data->ptr, UGP_IO_CV_D1, data->getName());
        break;
      default:
        if (mpi_rank==0) cout<<"Warning: Snapshot io not setup for registered datatype: "<<data->getDatatype()<<endl;
      }
    }
  }

  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
  {
    if (data->checkFlag())
    {
      switch (data->getDatatype())
      {
      case NO_DATA:
        if (mpi_rank==0) cout<<" > node vector: "<<data->getName()<<endl;
        writeNodeD3(*data->ptr, UGP_IO_NO_D3, data->getName());
        break;
      case CV_DATA:
        if (mpi_rank==0) cout<<" > cv vector: "<<data->getName()<<endl;
        writeCVsD3(*data->ptr, UGP_IO_CV_D3, data->getName());
        break;
      default:
        if (mpi_rank==0) cout<<"Warning: Snapshot io not setup for registered datatype: "<<data->getDatatype()<<endl;
      }
    }
  }

  // ------------------------------------------
  // eof...
  // ------------------------------------------

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }
  if (mpi_rank==0)
  {
    MPI_Status status;
    Header header;
    header.id = UGP_IO_EOF;
    sprintf(header.name, "EOF");
    header.skip = header_size;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);
  }
  offset += header_size;

  // trim the file to the current size and close
  MPI_File_set_size(fh, offset);
  MPI_File_close(&fh);

  // timing...
  MPI_Barrier(mpi_comm);
  if (mpi_rank==0)
  {
    wtime = MPI_Wtime()-wtime;
    cout<<" > writeSnapshot timing: np: "<<mpi_size<<", size[B]: "<<offset<<", time[s]: "<<wtime<<", rate[MB/s]: "
        <<(double) offset/1.0E+6/wtime<<endl;
  }

}

void Ugp::readSnapshot(const int index)
{

  char filename[32];
  sprintf(filename, "snapshot.%06d.out", index);

  if (mpi_rank==0) cout<<"readSnapshot: "<<filename<<endl;

  // when data has been read from the restart file...

  clearDataFlags();

  // timing...

  MPI_Barrier(mpi_comm);
  double wtime;
  if (mpi_rank==0) wtime = MPI_Wtime();

  // open...

  if (MPI_File_open(mpi_comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestart: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  byte_swap = 0;

  int itmp[2];
  MPI_Status status;
  MPI_File_read_all(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      if (mpi_rank==0) cerr<<"Error: readSnapshot: file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    if (mpi_rank==0) cout<<"File requires byte swapping."<<endl;
    byte_swap = 1;
  }

  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...

  offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  Header header;
  int done = 0;
  while (done!=1)
  {

    if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: readRestart: could not set view for header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    if (MPI_File_read_all(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestart: could not read all for header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // check...
    /*
     int count = 0;
     MPI_Get_count(&status,MPI_Header,&count);
     if (count != 1) {
     cerr << "Error: readRestart: could not get count for header expected at offset: " << offset << endl;
     throw(-1);
     }
     */

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    switch (header.id)
    {
    case UGP_IO_I0:
    {
      IntValue * data = getIntValueData(header.name);
      if (data!=NULL)
      {
        if (mpi_rank==0) cout<<" > reading int value: "<<header.name<<": "<<header.idata[0]<<endl;
        *data->ptr = header.idata[0];
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping int value: "<<header.name<<": "<<header.idata[0]<<endl;
      }
    }
      break;
    case UGP_IO_D0:
    {
      DoubleValue * data = getDoubleValueData(header.name);
      if (data!=NULL)
      {
        if (mpi_rank==0) cout<<" > reading double value: "<<header.name<<": "<<header.rdata[0]<<endl;
        *data->ptr = header.rdata[0];
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping double value: "<<header.name<<": "<<header.rdata[0]<<endl;
      }
    }
      break;
    case UGP_IO_CV_D1:
    {
      DoubleScalar * data = getScalarData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA scalar: "<<header.name<<endl;
        readCvD1(*data->ptr, header);
        dumpScalarRange(*data->ptr, ncv, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA scalar: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_CV_D3:
    {
      DoubleVector * data = getVectorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA vector: "<<header.name<<endl;
        readCvD3(*data->ptr, header);
        dumpVectorRange(*data->ptr, ncv, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA vector: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_CV_D33:
    {
      DoubleTensor * data = getTensorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA tensor: "<<header.name<<endl;
        readCvD33(*data->ptr, header);
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA tensor: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_NO_D1:
    {
      DoubleScalar * data = getScalarData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==NO_DATA))
      {
        if (mpi_rank==0) cout<<" > reading NO_DATA scalar: "<<header.name<<endl;
        readNoD1(*data->ptr, header);
        dumpScalarRange(*data->ptr, nno, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping NO_DATA scalar: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_NO_D3:
    {
      DoubleVector * data = getVectorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==NO_DATA))
      {
        if (mpi_rank==0) cout<<" > reading NO_DATA vector: "<<header.name<<endl;
        readNoD3(*data->ptr, header);
        dumpVectorRange(*data->ptr, nno, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping NO_DATA vector: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_EOF:
      done = 1;
      break;
    default:
      if (mpi_rank==0) cout<<" > skipping unrecognized header: "<<header.name<<endl;
    }

    offset += header.skip;

  }

  // and close...
  MPI_File_close(&fh);

  // timing...
  MPI_Barrier(mpi_comm);
  if (mpi_rank==0)
  {
    wtime = MPI_Wtime()-wtime;
    cout<<" > readSnapshot timing: nno: "<<no_count<<", nfa: "<<fa_count<<", ncv: "<<cv_count<<", np: "<<mpi_size
        <<", size[B]: "<<offset<<", time[s]: "<<wtime<<", rate[MB/s]: "<<(double) offset/1.0E+6/wtime<<endl;
  }

}

// ####################################################

void Ugp::writeI0(const int data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > int value: "<<name<<", value: "<<data<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeI0: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size;
    header.idata[0] = data;

    MPI_Status status;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // adjust offset...

  offset += header_size;

}

void Ugp::writeD0(const double data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > double value: "<<name<<", value: "<<data<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeD0: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size;
    header.rdata[0] = data;

    MPI_Status status;
    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // adjust offset...

  offset += header_size;

}

void Ugp::writeCVsD1(double *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > cv double scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*cv_count;
    header.idata[0] = cv_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);
  }

  // everybody pack...

  double * tmp = new double[my_cv_count];

  for (int i = 0; i<my_cv_count; i++)
  {
    int icv = my_cv_list[i];
    assert((icv>=0)&&(icv<ncv));
    tmp[i] = data[icv];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_cv_count, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*cv_count;

  // cleanup...

  delete[] tmp;

}

void Ugp::writeCVsD3(double(*data)[3], const int id, char * name)
{

  if (mpi_rank==0) cout<<" > cv double vector: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD3: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*cv_count*3;
    header.idata[0] = cv_count;
    header.idata[1] = 3;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  double (*tmp)[3] = new double[my_cv_count][3];

  for (int i = 0; i<my_cv_count; i++)
  {
    int icv = my_cv_list[i];
    for (int j = 0; j<3; j++)
    {
      tmp[i][j] = data[icv][j];
    }
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d3_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD3: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_cv_count*3, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*cv_count*3;

  // cleanup

  delete[] tmp;

}

void Ugp::writeCVsD33(double(*data)[3][3], const int id, char * name)
{

  if (mpi_rank==0) cout<<" > cv double tensor: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD33: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*cv_count*9;
    header.idata[0] = cv_count;
    header.idata[1] = 3;
    header.idata[2] = 3;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  double (*tmp)[3][3] = new double[my_cv_count][3][3];

  for (int i = 0; i<my_cv_count; i++)
  {
    int icv = my_cv_list[i];
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        tmp[i][j][k] = data[icv][j][k];
      }
    }
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d33_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCVsD33: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_cv_count*9, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*cv_count*9;

  // cleanup

  delete[] tmp;

}

void Ugp::writeNodeI1(int *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > no int scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+int_size*no_count;
    header.idata[0] = no_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  int * tmp = new int[my_no_count];

  for (int i = 0; i<my_no_count; i++)
  {
    int ino = my_no_list[i];
    assert((ino>=0)&&(ino<nno));
    tmp[i] = data[ino];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_INT, no_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_no_count, MPI_INT, &status);

  // update offset...

  offset += header_size+int_size*no_count;

  // cleanup...

  delete[] tmp;

}

void Ugp::writeNodeD1(double *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > no double scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*no_count;
    header.idata[0] = no_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  double * tmp = new double[my_no_count];

  for (int i = 0; i<my_no_count; i++)
  {
    int ino = my_no_list[i];
    assert((ino>=0)&&(ino<nno));
    tmp[i] = data[ino];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, no_d1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_no_count, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*no_count;

  // cleanup...

  delete[] tmp;

}

void Ugp::writeNodeD3(double(*data)[3], const int id, char * name)
{

  if (mpi_rank==0) cout<<" > no double vector: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeD3: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*no_count*3;
    header.idata[0] = no_count;
    header.idata[1] = 3;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  double (*tmp)[3] = new double[my_no_count][3];

  for (int i = 0; i<my_no_count; i++)
  {
    int ino = my_no_list[i];
    for (int j = 0; j<3; j++)
    {
      tmp[i][j] = data[ino][j];
    }
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, no_d3_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNodeD3: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_no_count*3, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*no_count*3;

  // cleanup

  delete[] tmp;

}

void Ugp::writeFaceI1(int *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > fa int scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeFaceI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+int_size*fa_count;
    header.idata[0] = fa_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  int * tmp = new int[my_fa_count];

  for (int i = 0; i<my_fa_count; i++)
  {
    int ifa = my_fa_list[i];
    assert((ifa>=0)&&(ifa<nfa));
    tmp[i] = data[ifa];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_INT, fa_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeFaceI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_fa_count, MPI_INT, &status);

  // update offset...

  offset += header_size+int_size*fa_count;

  // cleanup

  delete[] tmp;

}

void Ugp::writeFaceD1(double *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > fa double scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeFaceD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+double_size*fa_count;
    header.idata[0] = fa_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  double * tmp = new double[my_fa_count];

  for (int i = 0; i<my_fa_count; i++)
  {
    int ifa = my_fa_list[i];
    assert((ifa>=0)&&(ifa<nfa));
    tmp[i] = data[ifa];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, fa_d1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeFaceD1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_fa_count, MPI_DOUBLE, &status);

  // update offset...

  offset += header_size+double_size*fa_count;

  // cleanup...

  delete[] tmp;

}

// ##############################################################

void Ugp::writeNoofa()
{

  if (mpi_rank==0) cout<<" > noofa"<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNoofa: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, "NOOFA_I_AND_V");
    header.id = UGP_IO_NOOFA_I_AND_V;
    header.skip = header_size+int_size*fa_count+int_size*noofa_count;
    header.idata[0] = fa_count;
    header.idata[1] = noofa_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // increment offset...

  offset += header_size;

  // everybody pack...

  assert(my_noofa_count>my_fa_count);
  int * tmp = new int[my_noofa_count];

  for (int i = 0; i<my_fa_count; i++)
  {
    int ifa = my_fa_list[i];
    assert((ifa>=0)&&(ifa<nfa));
    tmp[i] = noofa_i[ifa+1]-noofa_i[ifa]; // just the node count...
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset, MPI_INT, fa_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNoofa: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_fa_count, MPI_INT, &status);

  // increment offset...

  offset += int_size*fa_count;

  // everybody pack...

  int ii = 0;
  for (int i = 0; i<my_fa_count; i++)
  {
    int ifa = my_fa_list[i];
    int nof_f = noofa_i[ifa];
    int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof<=nof_l; nof++)
    {
      // recall no_flag contains the 0-based global node index
      tmp[ii++] = no_flag[noofa_v[nof]];
    }
  }
  assert(ii==my_noofa_count);

  // everybody set the view...

  if (MPI_File_set_view(fh, offset, MPI_INT, noofa_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeNoofa: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_noofa_count, MPI_INT, &status);

  // increment offset...

  offset += int_size*noofa_count;

  // cleanup...

  delete[] tmp;

}

// ###########################################################

void Ugp::writeCvofa()
{

  if (mpi_rank==0) cout<<" > cvofa"<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvofa: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, "CVOFA");
    header.id = UGP_IO_CVOFA;
    header.skip = header_size+int_size*fa_count*2;
    header.idata[0] = fa_count;
    header.idata[1] = 2;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // increment offset...

  offset += header_size;

  // everybody pack...

  // the global cv index is currently in cv_flag, and the global face index is currently in 
  // fa_flag. For cvofa, we want to write the 2 valid cvs for all internal faces, the
  // valid cv and the valid face pair for periodic faces, and the valid cv and -1 for 
  // the boundary faces. To avoid ambiguity, periodic face pairings are written as -ifa-2, to 
  // allow -1 to have the special meaning of a boundary face.

  int * tmp = new int[nfa];

  // here we assume faces are sorted as expected in Ugp...
  for (int ifa = 0; ifa<nfa_b; ifa++)
    tmp[ifa] = -1; // assigned boundary
  for (int ifa = nfa_b; ifa<nfa_bp; ifa++)
    tmp[ifa] = -fa_flag[ifa]-2; // periodic boundary
  for (int ifa = nfa_bp; ifa<nfa_bpi; ifa++)
    tmp[ifa] = cv_flag[cvofa[ifa][0]];
  for (int ifa = nfa_bpi; ifa<nfa; ifa++)
    tmp[ifa] = -1;
  updateFaI1(tmp, REPLACE_DATA);
  // check...
  for (int ifa = 0; ifa<nfa_b; ifa++)
    assert(tmp[ifa]==-1); // assigned boundary
  for (int ifa = nfa_b; ifa<nfa_bp; ifa++)
    assert((tmp[ifa]<=-2)&&(tmp[ifa]!=-fa_flag[ifa]-2));
  for (int ifa = nfa_bp; ifa<nfa_bpi; ifa++)
    assert((tmp[ifa]>=0)&&(tmp[ifa]!=cv_flag[cvofa[ifa][0]]));
  for (int ifa = nfa_bpi; ifa<nfa; ifa++)
    assert(tmp[ifa]==-1);

  int (*tmp2)[2] = new int[my_fa_count][2];

  for (int i = 0; i<my_fa_count; i++)
  {
    int ifa = my_fa_list[i];
    assert((ifa>=0)&&(ifa<nfa));
    tmp2[i][0] = cv_flag[cvofa[ifa][0]];
    if (ifa<nfa_b) tmp2[i][1] = -1;
    else if (ifa<nfa_bpi) tmp2[i][1] = tmp[ifa];
    else tmp2[i][1] = cv_flag[cvofa[ifa][1]];
  }

  // don't need tmp anymore...

  delete[] tmp;

  // everybody set the view...

  if (MPI_File_set_view(fh, offset, MPI_INT, fa_i2_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvofa: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp2, my_fa_count*2, MPI_INT, &status);

  // increment offset...

  offset += int_size*fa_count*2;

  // cleanup...

  delete[] tmp2;

}

// ###########################################################

void Ugp::writeCvI1(int *data, const int id, char * name)
{

  if (mpi_rank==0) cout<<" > cv int scalar: "<<name<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, name);
    header.id = id;
    header.skip = header_size+int_size*cv_count;
    header.idata[0] = cv_count;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  int * tmp = new int[my_cv_count];

  for (int i = 0; i<my_cv_count; i++)
  {
    int icv = my_cv_list[i];
    assert((icv>=0)&&(icv<ncv));
    tmp[i] = data[icv];
  }

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_INT, cv_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_cv_count, MPI_INT, &status);

  // update offset...

  offset += header_size+int_size*cv_count;

  // cleanup...

  delete[] tmp;

}

void Ugp::writeCvPart()
{

  if (mpi_rank==0) cout<<" > cv part"<<endl;

  // do we all have to set the header view? ...

  if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvPart: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // just rank zero writes the header...

  MPI_Status status;
  if (mpi_rank==0)
  {

    Header header;
    sprintf(header.name, "CV_PART");
    header.id = UGP_IO_CV_PART;
    header.skip = header_size+int_size*cv_count;
    header.idata[0] = cv_count;
    header.idata[1] = mpi_size;

    MPI_File_write(fh, &header, 1, MPI_Header, &status);

  }

  // everybody pack...

  int * tmp = new int[my_cv_count];

  for (int i = 0; i<my_cv_count; i++)
    tmp[i] = mpi_rank;

  // everybody set the view...

  if (MPI_File_set_view(fh, offset+header_size, MPI_INT, cv_i1_type, "native", MPI_INFO_NULL)!=0)
  {
    cerr<<"Error: writeCvI1: could not set view at offset: "<<offset<<endl;
    throw(-1);
  }

  // everybody write...

  MPI_File_write_all(fh, tmp, my_cv_count, MPI_INT, &status);

  // update offset...

  offset += header_size+int_size*cv_count;

  // cleanup...

  delete[] tmp;

}

// #########################################################

int Ugp::readRestartNoX(double(*x_no)[3], const int i0, const int i1, const char * filename)
{

  assert((i0>=0)&&(i1>=i0));

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }

  // non-collective: open with MPI_COMM_SELF...

  MPI_File fh;
  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(MPI_COMM_SELF, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestartNoX: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  int byte_swap = 0;
  int itmp[2];
  MPI_Status status;
  MPI_File_read(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      cerr<<"Error: readRestartNoX: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    byte_swap = 1;
  }
  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: readRestartNoX: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...
  MPI_Offset offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  int i1_actual = -1;
  int done = 0;
  while (done!=1)
  {

    // read the next header...
    if (MPI_File_seek(fh, offset, MPI_SEEK_SET)!=0)
    {
      cerr<<"Error: readRestartNoX: could MPI_File_seek: "<<offset<<endl;
      throw(-1);
    }

    Header header;
    if (MPI_File_read(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestartNoX: could not read header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    if (header.id==UGP_IO_EOF)
    {
      cerr<<"Error: readRestartNoX: got UGP_IO_EOF"<<endl;
      throw(-1);
    }
    else if (header.id==UGP_IO_X_NO)
    {
      int this_nno = header.idata[0];
      assert(header.idata[1]==3);
      if (i0<this_nno)
      {
        i1_actual = min(i1, this_nno-1);
        if (i0>0) MPI_File_seek(fh, double_size*3*i0, MPI_SEEK_CUR);
        MPI_File_read(fh, x_no, 3*(i1_actual-i0+1), MPI_DOUBLE, &status);
        if (byte_swap) byteSwap(x_no, (i1_actual-i0+1), 3);
      }
      done = 1;
    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

  return (i1_actual);

}

int Ugp::readRestartNoXCollective(double(*x_no)[3], const int i0, const int i1, const char * filename)
{

  assert((i0>=0)&&(i1>=i0));

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }

  MPI_File fh;
  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(mpi_comm, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestartNoXCollective: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  int byte_swap = 0;
  int itmp[2];
  MPI_Status status;
  MPI_File_read_all(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      cerr<<"Error: readRestartNoXCollective: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    byte_swap = 1;
  }
  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: readRestartNoXCollective: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...
  MPI_Offset offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  int i1_actual = -1;
  int done = 0;
  while (done!=1)
  {

    // read the next header...
    if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: readRestartNoXCollective: could not set view: "<<offset<<endl;
      throw(-1);
    }

    Header header;
    if (MPI_File_read_all(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestartNoX: could not read header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    if (header.id==UGP_IO_EOF)
    {
      cerr<<"Error: readRestartNoX: got UGP_IO_EOF"<<endl;
      throw(-1);
    }
    else if (header.id==UGP_IO_X_NO)
    {
      int this_nno = header.idata[0];
      assert(header.idata[1]==3);
      if (i0<this_nno)
      {
        i1_actual = min(i1, this_nno-1);
        if (i0>0)
        {
          if (MPI_File_set_view(fh, offset+header_size+double_size*3*i0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)
              !=0)
          {
            cerr<<"Error: readRestartNoXCollective: could not set view: "<<offset<<endl;
            throw(-1);
          }
        }
        MPI_File_read_all(fh, x_no, 3*(i1_actual-i0+1), MPI_DOUBLE, &status);
        if (byte_swap) byteSwap(x_no, (i1_actual-i0+1), 3);
      }
      done = 1;
    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

  return (i1_actual);

}

int Ugp::readRestartNoR1(double *data, const int i0, const int i1, const char * name, const char * filename)
{

  assert((i0>=0)&&(i1>=i0));

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }

  // non-collective: open with MPI_COMM_SELF...

  MPI_File fh;
  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(MPI_COMM_SELF, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestartNoR1: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  int byte_swap = 0;
  int itmp[2];
  MPI_Status status;
  MPI_File_read(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      cerr<<"Error: readRestartNoR1: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    byte_swap = 1;
  }
  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: readRestartNoR1: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...
  MPI_Offset offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  int i1_actual = -1;
  int done = 0;
  while (done!=1)
  {

    // read the next header...
    if (MPI_File_seek(fh, offset, MPI_SEEK_SET)!=0)
    {
      cerr<<"Error: readRestartNoR1: could MPI_File_seek: "<<offset<<endl;
      throw(-1);
    }

    Header header;
    if (MPI_File_read(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestartNoR1: could not read header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    if (header.id==UGP_IO_EOF)
    {
      cerr<<"Error: readRestartNoR1: got UGP_IO_EOF, R1 data not in file: "<<name<<endl;
      throw(-1);
    }
    else if ((header.id==UGP_IO_NO_D1)&&(strcmp(header.name, name)==0))
    {
      int this_nno = header.idata[0];
      if (i0<this_nno)
      {
        i1_actual = min(i1, this_nno-1);
        if (i0>0) MPI_File_seek(fh, double_size*i0, MPI_SEEK_CUR);
        MPI_File_read(fh, data, (i1_actual-i0+1), MPI_DOUBLE, &status);
        if (byte_swap) byteSwap(data, (i1_actual-i0+1));
      }
      done = 1;
    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

  return (i1_actual);

}

int Ugp::readRestartNoR2(double(*data)[3], const int i0, const int i1, const char * name, const char * filename)
{

  assert((i0>=0)&&(i1>=i0));

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }

  // non-collective: open with MPI_COMM_SELF...

  MPI_File fh;
  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(MPI_COMM_SELF, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestartNoR2: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  int byte_swap = 0;
  int itmp[2];
  MPI_Status status;
  MPI_File_read(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      cerr<<"Error: readRestartNoR2: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    byte_swap = 1;
  }
  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: readRestartNoR2: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...
  MPI_Offset offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  int i1_actual = -1;
  int done = 0;
  while (done!=1)
  {

    // read the next header...
    if (MPI_File_seek(fh, offset, MPI_SEEK_SET)!=0)
    {
      cerr<<"Error: readRestartNoR2: could MPI_File_seek: "<<offset<<endl;
      throw(-1);
    }

    Header header;
    if (MPI_File_read(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestartNoR2: could not read header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    if (header.id==UGP_IO_EOF)
    {
      cerr<<"Error: readRestartNoR2: got UGP_IO_EOF, R2 data not in file: "<<name<<endl;
      throw(-1);
    }
    else if ((header.id==UGP_IO_NO_D3)&&(strcmp(header.name, name)==0))
    {
      int this_nno = header.idata[0];
      assert(header.idata[1]==3);
      if (i0<this_nno)
      {
        i1_actual = min(i1, this_nno-1);
        if (i0>0) MPI_File_seek(fh, double_size*3*i0, MPI_SEEK_CUR);
        MPI_File_read(fh, data, 3*(i1_actual-i0+1), MPI_DOUBLE, &status);
        if (byte_swap) byteSwap(data, (i1_actual-i0+1), 3);
      }
      done = 1;
    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

  return (i1_actual);

}

int Ugp::readRestartNoR2Collective(double(*data)[3], const int i0, const int i1, const char * name,
    const char * filename)
{

  assert((i0>=0)&&(i1>=i0));

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }

  MPI_File fh;
  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(mpi_comm, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestartNoR2: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  int byte_swap = 0;
  int itmp[2];
  MPI_Status status;
  MPI_File_read_all(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      cerr<<"Error: readRestartNoR2: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    byte_swap = 1;
  }
  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: readRestartNoR2: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...
  MPI_Offset offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  int i1_actual = -1;
  int done = 0;
  while (done!=1)
  {

    // read the next header...
    if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: readRestartNoXCollective: could not set view: "<<offset<<endl;
      throw(-1);
    }

    Header header;
    if (MPI_File_read_all(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestartNoR2: could not read header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    if (header.id==UGP_IO_EOF)
    {
      cerr<<"Error: readRestartNoR2: got UGP_IO_EOF, R2 data not in file: "<<name<<endl;
      throw(-1);
    }
    else if ((header.id==UGP_IO_NO_D3)&&(strcmp(header.name, name)==0))
    {
      int this_nno = header.idata[0];
      assert(header.idata[1]==3);
      if (i0<this_nno)
      {
        i1_actual = min(i1, this_nno-1);
        if (i0>0)
        {
          if (MPI_File_set_view(fh, offset+header_size+double_size*3*i0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)
              !=0)
          {
            cerr<<"Error: readRestartNoXCollective: could not set view: "<<offset<<endl;
            throw(-1);
          }
        }
        MPI_File_read_all(fh, data, 3*(i1_actual-i0+1), MPI_DOUBLE, &status);
        if (byte_swap) byteSwap(data, (i1_actual-i0+1), 3);
      }
      done = 1;
    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

  return (i1_actual);

}

// #########################################################

void Ugp::readRestart(const char * filename)
{

  if (mpi_rank==0) cout<<"readRestart: "<<filename<<endl;

  // io flags: we want to be able to read once, and write many... 

  if (io_common_flag==0)
  {
    initIoCommon();
    io_common_flag = 1;
  }
  if (io_read_flag!=0)
  {
    if (mpi_rank==0) cerr<<"Error: can only read restart once at present."<<endl;
    throw(-1);
  }
  if (io_write_flag!=0)
  {
    if (mpi_rank==0) cerr<<"Error: how can io_write_flag be set?"<<endl;
    throw(-1);
  }

  io_read_flag = 1;

  // clear the data flags - they are set to indicate
  // when data has been read from the restart file...

  clearDataFlags();

  // timing...

  MPI_Barrier(mpi_comm);
  double wtime;
  if (mpi_rank==0) wtime = MPI_Wtime();

  // open...

  char dummy[64];
  sprintf(dummy, "%s", filename);
  if (MPI_File_open(mpi_comm, dummy, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0)
  {
    cerr<<"Error: readRestart: could not open restart file: "<<filename<<endl;
    throw(-1);
  }

  // first 2 ints are:
  // 0. magic nummber
  // 1. i/o version

  // assume we do not have to byte_swap...
  byte_swap = 0;

  int itmp[2];
  MPI_Status status;
  MPI_File_read_all(fh, itmp, 2, MPI_INT, &status);
  if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
  {
    byteSwap(itmp, 2);
    if (itmp[0]!=UGP_IO_MAGIC_NUMBER)
    {
      if (mpi_rank==0) cerr<<"Error: readRestart: restart file does not start as expected. aborting."<<endl;
      throw(-1);
    }
    if (mpi_rank==0) cout<<"File requires byte swapping."<<endl;
    byte_swap = 1;
  }

  if (itmp[1]!=UGP_IO_VERSION)
  {
    cerr<<"Error: io version differs from current version: "<<itmp[1]<<endl;
    throw(-1);
  }

  // initialize the offset...

  offset = int_size*2;

  // remaining file consists of headers, potentially followed by data...

  Header header;
  int done = 0;
  while (done!=1)
  {

    if (MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL)!=0)
    {
      cerr<<"Error: readRestart: could not set view for header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    if (MPI_File_read_all(fh, &header, 1, MPI_Header, &status)!=0)
    {
      cerr<<"Error: readRestart: could not read all for header expected at offset: "<<offset<<endl;
      throw(-1);
    }

    // check...
    /*
     int count = 0;
     MPI_Get_count(&status,MPI_Header,&count);
     if (count != 1) {
     cerr << "Error: readRestart: could not get count for header expected at offset: " << offset << endl;
     throw(-1);
     }
     */

    // swap if neccessary...
    if (byte_swap) byteSwapHeader(&header, 1);

    switch (header.id)
    {
    case UGP_IO_NO_FA_CV_NOOFA_COUNTS:
      if (mpi_rank==0) cout<<" > global nno, nfa, ncv: "<<header.idata[0]<<" "<<header.idata[1]<<" "<<header.idata[2]
          <<endl;
      no_count = header.idata[0];
      fa_count = header.idata[1];
      cv_count = header.idata[2];
      noofa_count = header.idata[3];
      initReadRestart();
      break;
    case UGP_IO_NO_CHECK:
      readNodeCheck(header);
      break;
    case UGP_IO_FA_CHECK:
      readFaceCheck(header);
      break;
    case UGP_IO_CV_CHECK:
      readCvCheck(header);
      break;
    case UGP_IO_NOOFA_I_AND_V:
      readNoofa(header);
      break;
    case UGP_IO_CVOFA:
      readCvofa(header);
      break;
    case UGP_IO_FA_ZONE:
      readFaZone(header);
      break;
    case UGP_IO_CV_PART:
      readCvPart(header);
      break;
    case UGP_IO_X_NO:
      if (mpi_rank==0) cout<<" > reading NO_DATA vector: X_NO"<<endl;
      readNoD3(x_no, header);
      dumpVectorRange(x_no, nno, "X_NO");
      break;
    case UGP_IO_DATA:
      break;
    case UGP_IO_I0:
    {
      IntValue * data = getIntValueData(header.name);
      if (data!=NULL)
      {
        if (mpi_rank==0) cout<<" > reading int value: "<<header.name<<": "<<header.idata[0]<<endl;
        *data->ptr = header.idata[0];
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping int value: "<<header.name<<": "<<header.idata[0]<<endl;
      }
    }
      break;
    case UGP_IO_D0:
    {
      DoubleValue * data = getDoubleValueData(header.name);
      if (data!=NULL)
      {
        if (mpi_rank==0) cout<<" > reading double value: "<<header.name<<": "<<header.rdata[0]<<endl;
        *data->ptr = header.rdata[0];
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping double value: "<<header.name<<": "<<header.rdata[0]<<endl;
      }
    }
      break;
    case UGP_IO_FA_D1:
    {
      DoubleScalar * data = getScalarData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==FA_DATA))
      {
        if (mpi_rank==0) cout<<" > reading FA_DATA scalar: "<<header.name<<endl;
        readFaD1(*data->ptr, header);
        dumpScalarRange(*data->ptr, nfa, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping FA_DATA scalar: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_CV_D1:
    {
      DoubleScalar * data = getScalarData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA scalar: "<<header.name<<endl;
        readCvD1(*data->ptr, header);
        dumpScalarRange(*data->ptr, ncv, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA scalar: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_CV_D3:
    {
      DoubleVector * data = getVectorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA vector: "<<header.name<<endl;
        readCvD3(*data->ptr, header);
        dumpVectorRange(*data->ptr, ncv, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA vector: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_CV_D33:
    {
      DoubleTensor * data = getTensorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==CV_DATA))
      {
        if (mpi_rank==0) cout<<" > reading CV_DATA tensor: "<<header.name<<endl;
        readCvD33(*data->ptr, header);
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping CV_DATA tensor: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_NO_D1:
    {
      DoubleScalar * data = getScalarData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==NO_DATA))
      {
        if (mpi_rank==0) cout<<" > reading NO_DATA scalar: "<<header.name<<endl;
        readNoD1(*data->ptr, header);
        dumpScalarRange(*data->ptr, nno, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping NO_DATA scalar: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_NO_D3:
    {
      DoubleVector * data = getVectorData(header.name);
      if ((data!=NULL)&&(data->getDatatype()==NO_DATA))
      {
        if (mpi_rank==0) cout<<" > reading NO_DATA vector: "<<header.name<<endl;
        readNoD3(*data->ptr, header);
        dumpVectorRange(*data->ptr, nno, data->getName());
        data->setFlag();
      }
      else
      {
        if (mpi_rank==0) cout<<" > skipping NO_DATA vector: "<<header.name<<endl;
      }
    }
      break;
    case UGP_IO_EOF:
      done = 1;
      break;
    default:
      if (mpi_rank==0) cout<<" > skipping unrecognized header: "<<header.name<<endl;
    }

    offset += header.skip;

  }

  // and close...
  MPI_File_close(&fh);

  // now check that everything is OK...
  finalizeReadRestart();

  // timing...
  MPI_Barrier(mpi_comm);
  if (mpi_rank==0)
  {
    wtime = MPI_Wtime()-wtime;
    cout<<" > readRestart timing: nno: "<<no_count<<", nfa: "<<fa_count<<", ncv: "<<cv_count<<", np: "<<mpi_size
        <<", size[B]: "<<offset<<", time[s]: "<<wtime<<", rate[MB/s]: "<<(double) offset/1.0E+6/wtime<<endl;
  }

}

void Ugp::initReadRestart()
{

  // --------------------------------------
  // nodes...
  // --------------------------------------

  assert(noora==NULL);
  noora = new int[mpi_size+1];
  calcUniformDist(noora, no_count, mpi_size);

  assert(nno==0);
  nno = noora[mpi_rank+1]-noora[mpi_rank];

  int my_disp = noora[mpi_rank];
  MPI_Type_indexed(1, &nno, &my_disp, MPI_INT, &no_i1_type);
  MPI_Type_commit(&no_i1_type);

  MPI_Type_indexed(1, &nno, &my_disp, MPI_DOUBLE, &no_d1_type);
  MPI_Type_commit(&no_d1_type);

  int nno3 = 3*nno;
  int my_disp3 = 3*my_disp;
  MPI_Type_indexed(1, &nno3, &my_disp3, MPI_DOUBLE, &no_d3_type);
  MPI_Type_commit(&no_d3_type);

  // the no_flag is not used for anything yet, but must be
  // allocated as part of the minimum initialization... 

  assert(no_flag==NULL);
  no_flag = new int[nno];

  // the nodal coordinates are not registered data - 
  // allocate them mannually...   

  assert(x_no==NULL);
  x_no = new double[nno][3];

  // allocate any registered nodal data...

  allocateRegisteredNoData();

  // --------------------------------------
  // faces...
  // --------------------------------------

  assert(faora==NULL);
  faora = new int[mpi_size+1];
  calcUniformDist(faora, fa_count, mpi_size);

  assert(nfa==0);
  nfa = faora[mpi_rank+1]-faora[mpi_rank];
  my_disp = faora[mpi_rank];

  MPI_Type_indexed(1, &nfa, &my_disp, MPI_INT, &fa_i1_type);
  MPI_Type_commit(&fa_i1_type);

  MPI_Type_indexed(1, &nfa, &my_disp, MPI_DOUBLE, &fa_d1_type);
  MPI_Type_commit(&fa_d1_type);

  int nfa2 = nfa*2;
  int my_disp2 = my_disp*2;
  MPI_Type_indexed(1, &nfa2, &my_disp2, MPI_INT, &fa_i2_type);
  MPI_Type_commit(&fa_i2_type);

  assert(fa_flag==NULL);
  fa_flag = new int[nfa];

  // allocate any registered face data...

  allocateRegisteredFaData();

  // --------------------------------------
  // cvs...
  // --------------------------------------

  assert(cvora==NULL);
  cvora = new int[mpi_size+1];
  calcUniformDist(cvora, cv_count, mpi_size);

  assert(ncv==0);
  ncv = cvora[mpi_rank+1]-cvora[mpi_rank];
  my_disp = cvora[mpi_rank];

  MPI_Type_indexed(1, &ncv, &my_disp, MPI_INT, &cv_i1_type);
  MPI_Type_commit(&cv_i1_type);

  MPI_Type_indexed(1, &ncv, &my_disp, MPI_DOUBLE, &cv_d1_type);
  MPI_Type_commit(&cv_d1_type);

  int ncv3 = 3*ncv;
  my_disp3 = 3*my_disp;
  MPI_Type_indexed(1, &ncv3, &my_disp3, MPI_DOUBLE, &cv_d3_type);
  MPI_Type_commit(&cv_d3_type);

  int ncv33 = 9*ncv;
  int my_disp33 = 9*my_disp;
  MPI_Type_indexed(1, &ncv33, &my_disp33, MPI_DOUBLE, &cv_d33_type);
  MPI_Type_commit(&cv_d33_type);

  assert(cv_flag==NULL);
  cv_flag = new int[ncv];

  // allocate any other registered cv data...

  allocateRegisteredCvData();

}

void Ugp::readNodeCheck(Header &header)
{

  if (mpi_rank==0) cout<<" > readNodeCheck: ";

  assert(header.idata[0]==noora[mpi_size]);
  MPI_File_set_view(fh, offset+header_size, MPI_INT, no_i1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, no_flag, nno, MPI_INT, &status);
  if (byte_swap) byteSwap(no_flag, nno);

  for (int ino = 0; ino<nno; ino++)
  {
    if (ino+noora[mpi_rank]!=no_flag[ino])
    {
      cerr<<"Error: nodeCheck failed."<<endl;
      throw(-1);
    }
    // return no_flag to -1...
    no_flag[ino] = -1;
  }

  if (mpi_rank==0) cout<<"OK"<<endl;

}

void Ugp::readFaceCheck(Header &header)
{

  if (mpi_rank==0) cout<<" > readFaceCheck: ";

  assert(header.idata[0]==faora[mpi_size]);
  MPI_File_set_view(fh, offset+header_size, MPI_INT, fa_i1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, fa_flag, nfa, MPI_INT, &status);
  if (byte_swap) byteSwap(fa_flag, nfa);

  for (int ifa = 0; ifa<nfa; ifa++)
  {
    if (ifa+faora[mpi_rank]!=fa_flag[ifa])
    {
      cerr<<"Error: faceCheck failed."<<endl;
      throw(-1);
    }
    // return fa_flag to -1...
    fa_flag[ifa] = -1;
  }

  if (mpi_rank==0) cout<<"OK"<<endl;

}

void Ugp::readCvCheck(Header &header)
{

  if (mpi_rank==0) cout<<" > readCvCheck  : ";

  assert(header.idata[0]==cvora[mpi_size]);
  MPI_File_set_view(fh, offset+header_size, MPI_INT, cv_i1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, cv_flag, ncv, MPI_INT, &status);
  if (byte_swap) byteSwap(cv_flag, ncv);

  for (int icv = 0; icv<ncv; icv++)
  {
    if (icv+cvora[mpi_rank]!=cv_flag[icv])
    {
      cerr<<"Error: cvCheck failed."<<endl;
      throw(-1);
    }
    // return cv_flag to -1...
    cv_flag[icv] = -1;
  }

  if (mpi_rank==0) cout<<"OK"<<endl;

}

void Ugp::readNoofa(Header &header)
{

  assert(header.idata[0]==faora[mpi_size]);
  assert(header.idata[1]==noofa_count);

  assert(noofa_i==NULL);
  noofa_i = new int[nfa+1];

  MPI_File_set_view(fh, offset+header_size, MPI_INT, fa_i1_type, "native", MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, noofa_i+1, nfa, MPI_INT, &status);
  if (byte_swap) byteSwap(noofa_i+1, nfa);

  noofa_i[0] = 0;
  for (int ifa = 0; ifa<nfa; ifa++)
    noofa_i[ifa+1] += noofa_i[ifa];

  noofa_s = noofa_i[nfa];
  assert(noofa_v==NULL);
  noofa_v = new int[noofa_s];

  int my_disp;
  MPI_Scan(&noofa_s, &my_disp, 1, MPI_INT, MPI_SUM, mpi_comm);

  // the last processor should have the global node-of-face count...
  assert((mpi_rank<mpi_size-1)||(my_disp==header.idata[1]));

  // change to a 0-indexed displacement
  my_disp -= noofa_s;

  // build an indexed type and read noofa_i1_type...
  MPI_Type_indexed(1, &noofa_s, &my_disp, MPI_INT, &noofa_i1_type);
  MPI_Type_commit(&noofa_i1_type);

  MPI_File_set_view(fh, offset+header_size+int_size*faora[mpi_size], MPI_INT, noofa_i1_type, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, noofa_v, noofa_s, MPI_INT, &status);
  if (byte_swap) byteSwap(noofa_v, noofa_s);

}

void Ugp::readCvofa(Header &header)
{

  assert(header.idata[0]==faora[mpi_size]);
  assert(header.idata[1]==2);

  assert(cvofa==NULL);
  cvofa = new int[nfa][2];

  MPI_File_set_view(fh, offset+header_size, MPI_INT, fa_i2_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, cvofa, 2*nfa, MPI_INT, &status);
  if (byte_swap) byteSwap(cvofa, nfa, 2);

}

void Ugp::readFaZone(Header &header)
{

  // this requires the face zone kind (in idata[0]) and the periodic 
  // transform, if any (in rdata)...
  int index = addFaZone(header.name, header.idata[0], header.rdata);

  // we now should have a unique positive index...
  assert(index>=0);

  // and set the fa_flag to this index in associated faces...
  for (int ifa_global = header.idata[1]; ifa_global<=header.idata[2]; ifa_global++)
  {
    int ifa_local = ifa_global-faora[mpi_rank]; // faora[mpi_rank] contains my disp
    if ((ifa_local>=0)&&(ifa_local<nfa))
    {
      assert(fa_flag[ifa_local]==-1);
      fa_flag[ifa_local] = index;
    }
  }

}

void Ugp::readFaD1(double *data, Header& header)
{

  assert(header.idata[0]==faora[mpi_size]);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, fa_d1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, nfa, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, nfa);
}

void Ugp::readCvPart(Header &header)
{

  assert(header.idata[0]==cvora[mpi_size]);
  assert(cv_part==NULL);
  cv_part = new int[ncv];
  MPI_File_set_view(fh, offset+header_size, MPI_INT, cv_i1_type, "native", MPI_INFO_NULL);
  MPI_Status status;
  MPI_File_read_all(fh, cv_part, ncv, MPI_INT, &status);
  if (byte_swap) byteSwap(cv_part, ncv);

  int my_max_part = -1;
  for (int icv = 0; icv<ncv; icv++)
  {
    assert(cv_part[icv]>=0);
    my_max_part = max(my_max_part, cv_part[icv]);
  }
  my_max_part++;
  MPI_Allreduce(&my_max_part, &npart, 1, MPI_INT, MPI_MAX, mpi_comm);

  assert(npart==header.idata[1]);

  if (mpi_rank==0) cout<<" > file contains partition info for npart: "<<npart<<endl;
}

void Ugp::readCvD1(double *data, Header& header)
{

  assert(header.idata[0]==cvora[mpi_size]);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, ncv, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, ncv);
}

void Ugp::readCvD3(double(*data)[3], Header& header)
{

  assert(header.idata[0]==cvora[mpi_size]);
  assert(header.idata[1]==3);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d3_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, 3*ncv, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, ncv, 3);
}

void Ugp::readCvD33(double(*data)[3][3], Header& header)
{

  assert(header.idata[0]==cvora[mpi_size]);
  assert(header.idata[1]==3);
  assert(header.idata[2]==3);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, cv_d33_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, 9*ncv, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, ncv, 3, 3);
}

void Ugp::readNoD1(double *data, Header& header)
{

  assert(header.idata[0]==noora[mpi_size]);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, no_d1_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, nno, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, nno);
}

void Ugp::readNoD3(double(*data)[3], Header& header)
{

  assert(header.idata[0]==noora[mpi_size]);
  assert(header.idata[1]==3);

  MPI_File_set_view(fh, offset+header_size, MPI_DOUBLE, no_d3_type, "native", MPI_INFO_NULL);

  MPI_Status status;
  MPI_File_read_all(fh, data, 3*nno, MPI_DOUBLE, &status);
  if (byte_swap) byteSwap(data, nno, 3);
}

void Ugp::finalizeReadRestart()
{

  // check that all faces have their face zone set in fa_flag...
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    assert(fa_flag[ifa]>=0);
  }

  // check that the nodes of all faces are valid and reasonable...
  for (int ifa = 0; ifa<nfa; ifa++)
  {
    int nnof = 0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++)
    {
      int ino = noofa_v[nof];
      // increment the number of nodes...
      nnof++;
      // node numbering should be global zero-indexed...
      assert((ino>=0)&&(ino<noora[mpi_size]));
      // ensure there are no nodal duplications...
      for (int nof2 = noofa_i[ifa]; nof2<nof; nof2++)
      {
        assert(noofa_v[nof2]!=ino);
      }
    }
    // each face should have between 3 and 8 nodes...
    assert((nnof>=3)&&(nnof<=8));
  }

  {

    // it would be much better to use a vector of faZone's and avoid this
    // table, but for now...
    buildFaKindTable();

    // check that the cvs are reasonable and valid...
    for (int ifa = 0; ifa<nfa; ifa++)
    {
      // first cv should always be valid...
      assert((cvofa[ifa][0]>=0)&&(cvofa[ifa][0]<cvora[mpi_size]));
      // second cv can be -1 for boundary faces or very negative for periodic faces
      // recall that this was a trick used for periodic reconnection in CDP restarts...
      assert(cvofa[ifa][1]<cvora[mpi_size]);
      // recall that the fa_flag contains the faZone index, so we can get the face zone
      // kind from the table we built above...
      if (cvofa[ifa][1]>=0)
      {
        // if the second cv is valid, then this must be an internal face...
        assert(fa_kind_table[fa_flag[ifa]]==FA_ZONE_INTERNAL);
      }
      else if (cvofa[ifa][1]==-1)
      {
        assert(fa_kind_table[fa_flag[ifa]]==FA_ZONE_BOUNDARY);
      }
      else
      {
        // must be one of the periodic zones...
        assert((fa_kind_table[fa_flag[ifa]]>=FA_ZONE_PERIODIC_FIRST)&&(fa_kind_table[fa_flag[ifa]]
            <=FA_ZONE_PERIODIC_LAST));
        // modify the cvofa to contain the matching global face index...
        cvofa[ifa][1] = -cvofa[ifa][1]-2;
        assert((cvofa[ifa][1]>=0)&&(cvofa[ifa][1]<faora[mpi_size]));
      }
    }

  }

  {
    assert(cvZoneList.size()==0);
    cvZoneList.push_back(CvZone());
    list<CvZone>::iterator iz = cvZoneList.begin();
    iz->setName("fluid");
    iz->setIndex(0);
    iz->setKind(CV_ZONE_FLUID);

    // for a minimum initialization of Ugp, the cv_flag contains the
    // index of the associated cvZone - always 0 in this case...
    for (int icv = 0; icv<ncv; icv++)
    {
      assert(cv_flag[icv]==-1);
      cv_flag[icv] = 0;
    }
  }

  // at this point, the minimum initialization of the Ugp is complete.
  // and we can free up some stuff...

  // XXXXX in the future, we may want to leave some of this stuff
  // to allow SNAPSHOT reads, etc...

}

void Ugp::cleanupReadRestart()
{

  MPI_Type_free(&no_i1_type);
  MPI_Type_free(&no_d1_type);
  MPI_Type_free(&no_d3_type);
  MPI_Type_free(&fa_i1_type);
  MPI_Type_free(&fa_i2_type);
  MPI_Type_free(&cv_i1_type);
  MPI_Type_free(&cv_d1_type);
  MPI_Type_free(&cv_d3_type);
  MPI_Type_free(&cv_d33_type);
  MPI_Type_free(&noofa_i1_type);
}

