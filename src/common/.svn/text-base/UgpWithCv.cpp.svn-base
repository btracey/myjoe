#include "UgpWithCv.h"
#include <math.h>

#include "Logging.h"
using namespace logging;

void UgpWithCv::addGhostCvs() {

  // we have to add ghosts for both internal (inter-processor) and
  // periodic boundaries...

  ncv_g = ncv;
  for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); fa_prcomm!=facePrcommList.end(); fa_prcomm++) {
    Prcomm * cv_prcomm = getCvPrcomm(fa_prcomm->getNbrRank());
    // size the pack and unpack vecs the same as the faces - they could potentially be less
    // but for now...
    cv_prcomm->packIndexVec.resize(fa_prcomm->packIndexVec.size());
    for (int i = 0; i<fa_prcomm->packIndexVec.size(); i++) {
      int ifa = fa_prcomm->packIndexVec[i];
      // the internal face is always valid...
      assert((cvofa[ifa][0]>=0)&&(cvofa[ifa][0]<nfa));
      cv_prcomm->packIndexVec[i] = cvofa[ifa][0];
    }
    cv_prcomm->unpackIndexVec.resize(fa_prcomm->unpackIndexVec.size());
    for (int i = 0; i<fa_prcomm->unpackIndexVec.size(); i++) {
      int ifa = fa_prcomm->unpackIndexVec[i];
      // the ghost face should not be valid (yet)...
      assert(cvofa[ifa][1]==-1);
      cv_prcomm->unpackIndexVec[i] = ncv_g;
      cvofa[ifa][1] = ncv_g;
      ncv_g++;
    }
    // for proper transformations at periodic boundaries, it is neccessary (for the packing atleast)
    // to specify the Range's that divide up the Prcomm. Since this is an exact copy, this is trivial -
    // just copy the Face's packRanges...
    for (list<Range>::iterator ri = fa_prcomm->packRangeList.begin(); ri!=fa_prcomm->packRangeList.end(); ri++) {
      cv_prcomm->addPackRange(ri->getIndexFirst(), ri->getIndexLast(), ri->getBits(), ri->getFlag(), ri->dxyz);
    }
  }

  assert(ncv_g>=ncv);

  // --------------------------------------------------
  // adjust data sizes to include ghosts in CV_DATA
  // --------------------------------------------------

  // reallocate the cv_flag for this new cv size...
  if (ncv_g>ncv) {
    delete[] cv_flag;
    cv_flag = new int[ncv_g];
  }

  // double scalars: [ncv]...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g>ncv) {
        double * tmp = new double[ncv_g];
        for (int icv = 0; icv<ncv; icv++)
          tmp[icv] = (*(data->ptr))[icv];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_DATA);
    }
  }

  // double vectors: [ncv][3]...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g>ncv) {
        double (*tmp)[3] = new double[ncv_g][3];
        for (int icv = 0; icv<ncv; icv++)
          for (int j = 0; j<3; j++)
            tmp[icv][j] = (*(data->ptr))[icv][j];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

  // double tensors: [ncv][3][3]...
  for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g>ncv) {
        double (*tmp)[3][3] = new double[ncv_g][3][3];
        for (int icv = 0; icv<ncv; icv++)
          for (int j = 0; j<3; j++)
            for (int k = 0; k<3; k++)
              tmp[icv][j][k] = (*(data->ptr))[icv][j][k];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

  //cout << "ncv, ncv_g-ncv: " << ncv << " " << ncv_g-ncv << endl;

}

void UgpWithCv::addGhostAndFakeCvs() {

  if (mpi_rank==0) cout<<" > addGhostAndFakeCvs()"<<endl;

  // we have to add ghosts for both internal (inter-processor) and
  // periodic boundaries, AND fakes for local boundary faces...

  // put a global index into cv_flag...
  assert(cvora[mpi_rank+1]-cvora[mpi_rank]==ncv);
  FOR_ICV
    cv_flag[icv] = icv+cvora[mpi_rank];

  // now copy this global index into fa_flag along inter-processor boundaries...
  for (int ifa = 0; ifa<nfa_b; ++ifa)
    fa_flag[ifa] = -1;
  for (int ifa = nfa_b; ifa<nfa_bpi; ++ifa)
    fa_flag[ifa] = cv_flag[cvofa[ifa][0]];
  for (int ifa = nfa_bpi; ifa<nfa; ++ifa)
    fa_flag[ifa] = -2;
  updateFaData(fa_flag, REPLACE_DATA);

  // check faces are not overwritten...
  for (int ifa = 0; ifa<nfa_b; ++ifa)
    assert(fa_flag[ifa]==-1);
  for (int ifa = nfa_bpi; ifa<nfa; ++ifa)
    assert(fa_flag[ifa]==-2);

  // for indexing local stuff...
  for (int icv = 0; icv<ncv; ++icv)
    cv_flag[icv] = 0;

  // we also need a cv_tmp that is as large as the largest
  // cv count on a single processor
  int ncv_max = 0;
  for (int rank = 0; rank<mpi_size; ++rank)
    ncv_max = max(ncv_max, cvora[rank+1]-cvora[rank]);

  int * cv_flag_nbr = new int[ncv_max];
  for (int icv = 0; icv<ncv_max; ++icv)
    cv_flag_nbr[icv] = -1;

  // ncv_g starts at current ncv...
  ncv_g = ncv;
  for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); fa_prcomm!=facePrcommList.end(); fa_prcomm++) {

    // nbr rank stuff...
    int nbr_rank = fa_prcomm->getNbrRank();

    // new paired communicator associated with this nbr_rank...
    Prcomm * cv_prcomm = getCvPrcomm(nbr_rank);

    // cycle through face ranges. There will be a corresponding range in
    // each cv paired communicator, however the size may be slightly smaller
    // due to common cv's...

    // pack stuff...

    for (int iter = 0; iter<2; ++iter) {
      int cv_index_l = -1;
      for (list<Range>::iterator fa_range = fa_prcomm->packRangeList.begin(); fa_range!=fa_prcomm->packRangeList.end(); fa_range++) {
        // start cv_index_f at the last plus 1...
        int cv_index_f = cv_index_l+1;
        // the faces in this range are...
        int fa_index_f = fa_range->getIndexFirst();
        int fa_index_l = fa_range->getIndexLast();
        for (int fa_index = fa_index_f; fa_index<=fa_index_l; ++fa_index) {
          int ifa = fa_prcomm->packIndexVec[fa_index];
          // should be a periodic/interprocessor face...
          assert(fa_flag[ifa]>=0);
          // the internal cv should always be valid...
          int icv0 = cvofa[ifa][0];
          assert((icv0>=0)&&(icv0<ncv));
          // if it has not been packed yet for this range, then count it...
          if (cv_flag[icv0]==0) {
            cv_flag[icv0] = 1;
            cv_index_l += 1;
            if (iter==1) cv_prcomm->packIndexVec[cv_index_l] = icv0;
          }
        }
        // clear flag...
        for (int fa_index = fa_index_f; fa_index<=fa_index_l; ++fa_index) {
          int ifa = fa_prcomm->packIndexVec[fa_index];
          int icv0 = cvofa[ifa][0];
          cv_flag[icv0] = 0;
        }
        if (iter==0) {
          cv_prcomm->addPackRange(cv_index_f, cv_index_l, fa_range->getBits(), fa_range->getFlag(), fa_range->dxyz);
        }
      }
      if (iter==0) cv_prcomm->packIndexVec.resize(cv_index_l+1);
    }

    // unpack stuff...

    int ncv_g0 = ncv_g;
    for (int iter = 0; iter<2; ++iter) {
      ncv_g = ncv_g0;
      int cv_index_l = -1;
      for (list<Range>::iterator fa_range = fa_prcomm->unpackRangeList.begin(); fa_range
          !=fa_prcomm->unpackRangeList.end(); fa_range++) {
        // start cv_index_f at the last plus 1...
        int cv_index_f = cv_index_l+1;
        // the faces in this range are...
        int fa_index_f = fa_range->getIndexFirst();
        int fa_index_l = fa_range->getIndexLast();
        for (int fa_index = fa_index_f; fa_index<=fa_index_l; ++fa_index) {
          int ifa = fa_prcomm->unpackIndexVec[fa_index];
          // should be a periodic/interprocessor face...
          int icv_global = fa_flag[ifa];
          assert((icv_global>=cvora[nbr_rank])&&(icv_global<cvora[nbr_rank+1]));
          int icv1 = icv_global-cvora[nbr_rank];
          if (cv_flag_nbr[icv1]==-1) {
            cv_index_l += 1;
            if (iter==1) cv_prcomm->unpackIndexVec[cv_index_l] = ncv_g;
            cv_flag_nbr[icv1] = ncv_g;
            ncv_g++;
          }
          if (iter==1) {
            assert(cvofa[ifa][1]<0);
            cvofa[ifa][1] = cv_flag_nbr[icv1];
          }
        }
        // clear flag...
        for (int fa_index = fa_index_f; fa_index<=fa_index_l; ++fa_index) {
          int ifa = fa_prcomm->unpackIndexVec[fa_index];
          int icv_global = fa_flag[ifa];
          assert((icv_global>=cvora[nbr_rank])&&(icv_global<cvora[nbr_rank+1]));
          int icv1 = icv_global-cvora[nbr_rank];
          cv_flag_nbr[icv1] = -1;
        }
        if (iter==0) {
          cv_prcomm->addUnpackRange(cv_index_f, cv_index_l, fa_range->getBits(), fa_range->getFlag(), fa_range->dxyz);
        }
      }
      if (iter==0) cv_prcomm->unpackIndexVec.resize(cv_index_l+1);
    }

  }

  // check and cleanup...
  assert(ncv_g>=ncv);
  delete[] cv_flag_nbr;

  // --------------------------------------------------
  // fakes are really easy: there is always a 1:1 correspondence
  // with boundary faces...
  // --------------------------------------------------

  ncv_gf = ncv_g;
  for (int ifa = 0; ifa<nfa_b; ++ifa) {
    assert(cvofa[ifa][1]<0);
    cvofa[ifa][1] = ncv_gf++;
  }

  // --------------------------------------------------
  // adjust data sizes to include ghosts in CV_DATA
  // --------------------------------------------------

  // i think this must always be larger...
  assert(ncv_gf>ncv);

  // reallocate the cv_flag for this new cv size...
  if (ncv_gf>ncv) {
    delete[] cv_flag;
    cv_flag = new int[ncv_gf];
  }

  // double scalars: [ncv]...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf>ncv) {
        double * tmp = new double[ncv_gf];
        for (int icv = 0; icv<ncv; icv++)
          tmp[icv] = (*(data->ptr))[icv];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_DATA);
    }
  }

  // double vectors: [ncv][3]...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf>ncv) {
        double (*tmp)[3] = new double[ncv_gf][3];
        for (int icv = 0; icv<ncv; icv++)
          for (int j = 0; j<3; j++)
            tmp[icv][j] = (*(data->ptr))[icv][j];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

  // double tensors: [ncv][3][3]...
  for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
    if (data->getDatatype()==CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf>ncv) {
        double (*tmp)[3][3] = new double[ncv_gf][3][3];
        for (int icv = 0; icv<ncv; icv++)
          for (int j = 0; j<3; j++)
            for (int k = 0; k<3; k++)
              tmp[icv][j][k] = (*(data->ptr))[icv][j][k];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

}

void UgpWithCv::calcGeometryFake() {

  if (mpi_rank==0) cout<<" > calcGeometry()"<<endl;

  // ======================================
  // face normal and centroid...
  // normal has area magnitude
  // ======================================

  if (x_fa==NULL) x_fa = new double[nfa][3];

  if (fa_normal==NULL) fa_normal = new double[nfa][3];

  for (int ifa = 0; ifa<nfa; ifa++) {

    for (int i = 0; i<3; i++)
      x_fa[ifa][i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        x_fa[ifa][i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);

    // we can compute the face normal directly with this approx centroid...

    for (int i = 0; i<3; i++)
      fa_normal[ifa][i] = 0.0;
    int ino2 = noofa_v[noofa_i[ifa+1]-1]; // last node in the list
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
      int ino1 = ino2;
      ino2 = noofa_v[nof];
      double v1[3];
      double v2[3];
      for (int i = 0; i<3; i++) {
        v1[i] = x_no[ino1][i]-x_fa[ifa][i];
        v2[i] = x_no[ino2][i]-x_fa[ifa][i];
      }
      fa_normal[ifa][0] += 0.5*(v1[1]*v2[2]-v1[2]*v2[1]);
      fa_normal[ifa][1] += 0.5*(v1[2]*v2[0]-v1[0]*v2[2]);
      fa_normal[ifa][2] += 0.5*(v1[0]*v2[1]-v1[1]*v2[0]);
    }

    double fa_area = sqrt(fa_normal[ifa][0]*fa_normal[ifa][0]+fa_normal[ifa][1]*fa_normal[ifa][1]+fa_normal[ifa][2]
        *fa_normal[ifa][2]);
    double unit_normal[3];
    FOR_I3
      unit_normal[i] = fa_normal[ifa][i]/fa_area;

    // find true centroid...

    bool done = false;
    int iter = 0;
    while (!done) {
      double dx_fa[3] = { 0.0, 0.0, 0.0 };
      double x_fa_tmp[3] = { 0.0, 0.0, 0.0 };
      double area_sum = 0.0;
      int ino2 = noofa_v[noofa_i[ifa+1]-1]; // last node in the list
      for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
        int ino1 = ino2;
        ino2 = noofa_v[nof];
        double v1[3];
        double v2[3];
        for (int i = 0; i<3; i++) {
          v1[i] = x_no[ino1][i]-x_fa[ifa][i];
          v2[i] = x_no[ino2][i]-x_fa[ifa][i];
        }
        double this_normal[3];
        this_normal[0] = (v1[1]*v2[2]-v1[2]*v2[1]);
        this_normal[1] = (v1[2]*v2[0]-v1[0]*v2[2]);
        this_normal[2] = (v1[0]*v2[1]-v1[1]*v2[0]);
        /*
         double this_area = sqrt( this_normal[0]*this_normal[0] +
         this_normal[1]*this_normal[1] +
         this_normal[2]*this_normal[2] );
         */
        double this_area = this_normal[0]*unit_normal[0]+this_normal[1]*unit_normal[1]+this_normal[2]*unit_normal[2];
        FOR_I3
          x_fa_tmp[i] += this_area*(x_no[ino1][i]+x_no[ino2][i]+x_fa[ifa][i]);
        FOR_I3
          dx_fa[i] += this_area*(v1[i]+v2[i]);
        area_sum += this_area;
      }
      double d2 = 0.0;
      FOR_I3 {
        x_fa_tmp[i] /= 3.0*area_sum;
        dx_fa[i] /= 2.0*area_sum;
        x_fa[ifa][i] += dx_fa[i];
        d2 += dx_fa[i]*dx_fa[i];
      }
      // just set to the x_fa_tmp value...
      FOR_I3
        x_fa[ifa][i] = x_fa_tmp[i];
      if (++iter>200) {
        cerr<<"Error: face centroid did not converge: d2/area_sum: "<<d2/area_sum<<endl;
        throw(-1);
      }
      else if (iter>150) {
        cout<<" > face centroid not converging: "<<ifa<<" "<<iter<<" "<<d2/area_sum<<endl;
      }
      done = (d2/area_sum<1.0E-20);
    }

  }

  // ======================================
  // x_cv, cv_volume in ghost as well...
  // ======================================

  if (x_cv==NULL) x_cv = new double[ncv_gf][3];

  if (cv_volume==NULL) cv_volume = new double[ncv_gf];

  for (int icv = 0; icv<ncv; icv++) {
    // compute an approximate center as the mean of the face x_fa's...
    double x_cv_approx[3];
    for (int i = 0; i<3; i++)
      x_cv_approx[i] = 0.0;
    for (int foc = faocv_i[icv]; foc<faocv_i[icv+1]; foc++) {
      int ifa = faocv_v[foc];
      for (int i = 0; i<3; i++)
        x_cv_approx[i] += x_fa[ifa][i];
    }
    // divide by the number of faces...
    for (int i = 0; i<3; i++)
      x_cv_approx[i] /= (double) (faocv_i[icv+1]-faocv_i[icv]);
    // add tet volumes...
    cv_volume[icv] = 0.0;
    for (int i = 0; i<3; i++)
      x_cv[icv][i] = 0.0;
    // loop on the faces of this cv...
    for (int foc = faocv_i[icv]; foc<faocv_i[icv+1]; foc++) {
      int ifa = faocv_v[foc];
      double v1[3];
      for (int i = 0; i<3; i++)
        v1[i] = x_fa[ifa][i]-x_cv_approx[i];
      // check if the face is inward or outward wrt this cv...
      if (cvofa[ifa][0]==icv) {
        // face is outward - loop through edges in forward direction...
        int ino2 = noofa_v[noofa_i[ifa+1]-1];
        for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i<3; i++) {
            v2[i] = x_no[ino1][i]-x_cv_approx[i];
            v3[i] = x_no[ino2][i]-x_cv_approx[i];
          }
          // 2 nodes, the face, and the approx cv form a tet...
          double this_volume = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])+v1[1]*(v2[2]*v3[0]-v2[0]*v3[2])+v1[2]*(v2[0]*v3[1]-v2[1]
              *v3[0]);
          assert(this_volume>0.0); // check on the grid ordering/right-handedness
          cv_volume[icv] += this_volume;
          for (int i = 0; i<3; i++)
            x_cv[icv][i] += this_volume*(x_cv_approx[i]+x_fa[ifa][i]+x_no[ino1][i]+x_no[ino2][i]);
        }
      }
      else {
        assert(cvofa[ifa][1]==icv);
        // face is outward, loop through edges in backward direction...
        int ino2 = noofa_v[noofa_i[ifa]];
        for (int nof = noofa_i[ifa+1]-1; nof>=noofa_i[ifa]; nof--) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i<3; i++) {
            v2[i] = x_no[ino1][i]-x_cv_approx[i];
            v3[i] = x_no[ino2][i]-x_cv_approx[i];
          }
          double this_volume = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])+v1[1]*(v2[2]*v3[0]-v2[0]*v3[2])+v1[2]*(v2[0]*v3[1]-v2[1]
              *v3[0]);
          //assert(this_volume > 0.0);
          if (this_volume<=0.0) cout<<" > Warning: volume negative: "<<this_volume<<endl;

          cv_volume[icv] += this_volume;
          for (int i = 0; i<3; i++)
            x_cv[icv][i] += this_volume*(x_cv_approx[i]+x_fa[ifa][i]+x_no[ino1][i]+x_no[ino2][i]);
        }
      }
    }
    // normalize both...
    for (int i = 0; i<3; i++)
      x_cv[icv][i] /= 4.0*cv_volume[icv];
    cv_volume[icv] /= 6.0;
  }

  // update volumes and centroids across processors...
  updateCvData(cv_volume, REPLACE_DATA); // blocking exchanges
  updateCvData(x_cv, REPLACE_TRANSLATE_DATA); // "TRANSLATE" is for any periodicity

  // finally fake data...
  /*
   for (int ifa = 0; ifa < nfa_b; ifa++) {
   // copy the cv_volume...
   int icv0 = cvofa[ifa][0];
   assert( (icv0 >= 0)&&(icv0 < ncv) );
   int icv1 = cvofa[ifa][1];
   assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) );
   // copy the volume...
   cv_volume[icv1] = cv_volume[icv0];
   double nmag = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
   fa_normal[ifa][1]*fa_normal[ifa][1] +
   fa_normal[ifa][2]*fa_normal[ifa][2] );
   double shalf_dot_n = ( (x_fa[ifa][0] - x_cv[icv0][0])*fa_normal[ifa][0] +
   (x_fa[ifa][1] - x_cv[icv0][1])*fa_normal[ifa][1] +
   (x_fa[ifa][2] - x_cv[icv0][2])*fa_normal[ifa][2] )/nmag;
   assert( shalf_dot_n > 0.0 );

   // Method 1: relect the cv through the normal...
   //for (int i = 0; i < 3; i++)
   //   x_cv[icv1][i] = x_cv[icv0][i] + 2.0*shalf_dot_n*fa_normal[ifa][i]/nmag;

   // Method 2: put the cv off the face by the cv-fa normal distance. This
   // this is better for certain grids with hanging nodes at physical boundaries.
   for (int i = 0; i < 3; i++)
   x_cv[icv1][i] = x_fa[ifa][i] + shalf_dot_n*fa_normal[ifa][i]/nmag;
   }
   */

  // XXXXX: major mod - feb 16: put cv centroid right on face location, and give it zero volume...
  for (int ifa = 0; ifa<nfa_b; ifa++) {
    int icv1 = cvofa[ifa][1];
    assert((icv1>=ncv_g)&&(icv1<ncv_gf));
    // copy the volume...
    cv_volume[icv1] = 0.0;
    // and put the cv centroid right ontop of the face...
    FOR_I3
      x_cv[icv1][i] = x_fa[ifa][i];
  }

  // check...
  double my_volume_sum = 0.0;
  for (int icv = 0; icv<ncv; icv++)
    my_volume_sum += cv_volume[icv];
  double volume_sum;
  MPI_Reduce(&my_volume_sum, &volume_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  // also check s dot n...
  double my_buf[2];
  my_buf[0] = 1.0;
  my_buf[1] = 0.0;
  for (int ifa = 0; ifa<nfa_b; ifa++) {
    // the first faces are boundary faces, and should have icv1 >= ncv_g...
    assert((cvofa[ifa][1]>=ncv_g)&&(cvofa[ifa][1]<ncv_gf));
  }
  for (int ifa = nfa_b; ifa<nfa; ifa++) {
    // the remaining faces should have valid cvofa, with a value
    // in the ghost range for the periodic cvs, and in the
    // owned cv range for totally internal cvs...
    assert(cvofa[ifa][1]>=0);
    if (ifa<nfa_bpi) {
      assert((cvofa[ifa][1]>=ncv)&&(cvofa[ifa][1]<ncv_g));
    }
    else {
      assert(cvofa[ifa][1]<ncv);
    }
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    double s[3];
    double smag = 0.0;
    double nmag = 0.0;
    for (int i = 0; i<3; i++) {
      s[i] = x_cv[icv1][i]-x_cv[icv0][i];
      smag += s[i]*s[i];
      nmag += fa_normal[ifa][i]*fa_normal[ifa][i];
    }
    smag = sqrt(smag);
    nmag = sqrt(nmag);
    //cout << "ifa, smag, nmag: " << ifa << " " << smag << " " << nmag << endl;
    double n[3];
    for (int i = 0; i<3; i++) {
      s[i] /= smag;
      n[i] = fa_normal[ifa][i]/nmag;
    }
    double dp = s[0]*n[0]+s[1]*n[1]+s[2]*n[2];
    assert(dp>0.0);
    my_buf[0] = min(my_buf[0], dp);
    my_buf[1] = min(my_buf[1], -dp);
  }

  double buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);

  if (mpi_rank==0) cout<<" > calcGeometry(), cv volume: "<<volume_sum<<", s dot n min/max: "<<buf[0]<<" "<<-buf[1]
      <<endl;

  // also compute the min cv distance...
  double my_min_delta = 1.0E+20;
  FOR_IFA {
    const int icv0 = cvofa[ifa][0];
    assert((icv0>=0)&&(icv0<ncv));
    const int icv1 = cvofa[ifa][1];
    assert(icv1>=0);
    double d2 = 0.0;
    FOR_I3 {
      double dx = x_cv[icv1][i]-x_cv[icv0][i];
      d2 += dx*dx;
    }
    my_min_delta = min(my_min_delta, sqrt(d2));
  }
  double min_delta;
  MPI_Reduce(&my_min_delta, &min_delta, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
  if (mpi_rank==0) cout<<" > calcGeometry(), min distance between cvs: "<<min_delta<<endl;

  //MPI_Pause("done");

}

void UgpWithCv::buildNbocvFake() {

  //
  // nbocv_i/v (neighbor-of-cv) is the CSR (compact-storage-row) format
  // matrix that holds the immediate cv neighbours. The first element
  // of any row is the diagonal element, i.e.
  //
  // double * A = new double[nbocv_s];
  // for (int icv = 0; icv < ncv; icv++) {
  //   int noc_f = nbocv_i[icv];
  //   assert( nbocv_v[noc_f] == icv );
  //   A[noc_f] = ... // the diagonal
  //   ...
  //
  // requirements: face data must be ordered.
  //


  if (mpi_rank==0) cout<<" > buildNbocv()"<<endl;

  // cvora should be valid...
  assert(cvora[mpi_rank+1]-cvora[mpi_rank]==ncv);

  // these should not exist...
  assert(nbocv_i==NULL);
  assert(nbocv_v==NULL);

  nbocv_i = new int[ncv+1];
  FOR_ICV
    nbocv_i[icv+1] = 0; // used for counts in first iteration below...

  for (int icv = 0; icv<ncv; ++icv)
    cv_flag[icv] = -1; // locally owned
  for (int icv = ncv; icv<ncv_g; ++icv)
    cv_flag[icv] = -2; // ghost
  for (int icv = ncv_g; icv<ncv_gf; ++icv)
    cv_flag[icv] = -3; // fake

  for (int iter = 0; iter<2; iter++) {
    FOR_ICV {
      // diagonal first...
      assert(cv_flag[icv]==-1);
      cv_flag[icv] = 1;
      if (iter==0) {
        nbocv_i[icv+1] += 1;
      }
      else {
        nbocv_v[nbocv_i[icv]] = icv;
        nbocv_i[icv] += 1;
      }
      // get nbrs using faces...
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; ++foc) {
        const int ifa = faocv_v[foc];
        int icv_nbr;
        if (cvofa[ifa][0]==icv) {
          icv_nbr = cvofa[ifa][1];
        }
        else {
          assert(cvofa[ifa][1]==icv);
          icv_nbr = cvofa[ifa][0];
        }
        if (ifa<nfa_b) {
          // this is a boundary face, and should have a fake cv. It
          // is NOT included in the nbocv CSR structure...
          assert(cv_flag[icv_nbr]==-3);
          cv_flag[icv_nbr] = 3; // switch to positive for checking
        }
        else if (ifa<nfa_bpi) {
          // this is a ghost face. It may already be included in the
          // list because two faces can point to the same cv now:
          // (addresses non-coplanar face problem)
          if (cv_flag[icv_nbr]==-2) {
            cv_flag[icv_nbr] = 2; // flip
            if (iter==0) {
              nbocv_i[icv+1] += 1;
            }
            else {
              nbocv_v[nbocv_i[icv]] = icv_nbr;
              nbocv_i[icv] += 1;
            }
          }
          else {
            assert(cv_flag[icv_nbr]==2);
          }
        }
        else {
          // this is an internal face. It may already be included in the
          // list for the same reason as above...
          if (cv_flag[icv_nbr]==-1) {
            cv_flag[icv_nbr] = 1; // flip
            if (iter==0) {
              nbocv_i[icv+1] += 1;
            }
            else {
              nbocv_v[nbocv_i[icv]] = icv_nbr;
              nbocv_i[icv] += 1;
            }
          }
          else {
            assert(cv_flag[icv_nbr]==1);
          }
        }
      }
      // clear flags...
      for (int foc = foc_f; foc<=foc_l; ++foc) {
        const int ifa = faocv_v[foc];
        if (cv_flag[cvofa[ifa][0]]>0) cv_flag[cvofa[ifa][0]] = -cv_flag[cvofa[ifa][0]];
        if (cv_flag[cvofa[ifa][1]]>0) cv_flag[cvofa[ifa][1]] = -cv_flag[cvofa[ifa][1]];
      }
    }
    if (iter==0) {
      nbocv_i[0] = 0;
      FOR_ICV
        nbocv_i[icv+1] += nbocv_i[icv];
      nbocv_s = nbocv_i[ncv];
      nbocv_v = new int[nbocv_s];
    }
    else {
      // return nbocv_i...
      for (int icv = ncv-1; icv>0; --icv)
        nbocv_i[icv] = nbocv_i[icv-1];
      nbocv_i[0] = 0;
    }
  }

}

void UgpWithCv::buildNbocv() {

  //
  // nbocv_i/v (neighbor-of-cv) is the CSR (compact-storage-row) format
  // matrix that holds the immediate cv neighbours. The first element
  // of any row is the diagonal element, i.e.
  //
  // double * A = new double[nbocv_s];
  // for (int icv = 0; icv < ncv; icv++) {
  //   int noc_f = nbocv_i[icv];
  //   assert( nbocv_v[noc_f] == icv );
  //   A[noc_f] = ... // the diagonal
  //   ...
  //
  // requirements: face data must be ordered.
  //


  if (mpi_rank==0) cout<<"buildNbocv()"<<endl;

  assert(nbocv_i==NULL);
  assert(nbocv_v==NULL);

  nbocv_i = new int[ncv+1];
  for (int icv = 0; icv<ncv; icv++)
    nbocv_i[icv+1] = 1; // 1 for the diagonal

  for (int iter = 0; iter<2; iter++) {
    // skip boundary faces - these do not have icvs that
    // participate in the nbocv_i/v structure...
    for (int ifa = nfa_b; ifa<nfa; ifa++) {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert(icv1>=0);
      if (iter==0) {
        nbocv_i[icv0+1] += 1;
        if (icv1<ncv) nbocv_i[icv1+1] += 1;
      }
      else {
        nbocv_v[nbocv_i[icv0]] = icv1;
        nbocv_i[icv0] += 1;
        if (icv1<ncv) {
          nbocv_v[nbocv_i[icv1]] = icv0;
          nbocv_i[icv1] += 1;
        }
      }
    }
    if (iter==0) {
      nbocv_i[0] = 0;
      for (int icv = 0; icv<ncv; icv++)
        nbocv_i[icv+1] += nbocv_i[icv];
      nbocv_s = nbocv_i[ncv];
      nbocv_v = new int[nbocv_s];
      // add the diagonal...
      for (int icv = 0; icv<ncv; icv++) {
        nbocv_v[nbocv_i[icv]] = icv;
        nbocv_i[icv] += 1;
      }
    }
    else {
      for (int icv = ncv; icv>0; icv--)
        nbocv_i[icv] = nbocv_i[icv-1];
      nbocv_i[0] = 0;
    }
  }

}

/**
 *
 * face normal and centroid...
 * normal has area magnitude
 *
 */
void UgpWithCv::calcGeometry() {
  assert(x_fa==NULL);
  x_fa = new double[nfa][3]; // face centroid

  assert(fa_normal==NULL);
  fa_normal = new double[nfa][3]; // face normal

  for (int ifa = 0; ifa<nfa; ifa++) {
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] = 0.0;
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
      int ino = noofa_v[nof];
      for (int i = 0; i<3; i++)
        x_fa[ifa][i] += x_no[ino][i];
    }
    for (int i = 0; i<3; i++)
      x_fa[ifa][i] /= (double) (noofa_i[ifa+1]-noofa_i[ifa]);

    // we can compute the face normal directly with this approx centroid...

    for (int i = 0; i<3; i++)
      fa_normal[ifa][i] = 0.0;
    int ino2 = noofa_v[noofa_i[ifa+1]-1]; // last node in the list
    for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
      int ino1 = ino2;
      ino2 = noofa_v[nof];
      double v1[3];
      double v2[3];
      for (int i = 0; i<3; i++) {
        v1[i] = x_no[ino1][i]-x_fa[ifa][i];
        v2[i] = x_no[ino2][i]-x_fa[ifa][i];
      }
      fa_normal[ifa][0] += 0.5*(v1[1]*v2[2]-v1[2]*v2[1]);
      fa_normal[ifa][1] += 0.5*(v1[2]*v2[0]-v1[0]*v2[2]);
      fa_normal[ifa][2] += 0.5*(v1[0]*v2[1]-v1[1]*v2[0]);
    }

    double fa_area = sqrt(fa_normal[ifa][0]*fa_normal[ifa][0]+fa_normal[ifa][1]*fa_normal[ifa][1]+fa_normal[ifa][2]
        *fa_normal[ifa][2]);
    double unit_normal[3];
    FOR_I3
      unit_normal[i] = fa_normal[ifa][i]/fa_area;

    // find true centroid...

    bool done = false;
    int iter = 0;
    while (!done) {
      double dx_fa[3] = { 0.0, 0.0, 0.0 };
      double x_fa_tmp[3] = { 0.0, 0.0, 0.0 };
      double area_sum = 0.0;
      int ino2 = noofa_v[noofa_i[ifa+1]-1]; // last node in the list
      for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
        int ino1 = ino2;
        ino2 = noofa_v[nof];
        double v1[3];
        double v2[3];
        for (int i = 0; i<3; i++) {
          v1[i] = x_no[ino1][i]-x_fa[ifa][i];
          v2[i] = x_no[ino2][i]-x_fa[ifa][i];
        }
        double this_normal[3];
        this_normal[0] = (v1[1]*v2[2]-v1[2]*v2[1]);
        this_normal[1] = (v1[2]*v2[0]-v1[0]*v2[2]);
        this_normal[2] = (v1[0]*v2[1]-v1[1]*v2[0]);
        /*
         double this_area = sqrt( this_normal[0]*this_normal[0] +
         this_normal[1]*this_normal[1] +
         this_normal[2]*this_normal[2] );
         */
        double this_area = this_normal[0]*unit_normal[0]+this_normal[1]*unit_normal[1]+this_normal[2]*unit_normal[2];
        FOR_I3
          x_fa_tmp[i] += this_area*(x_no[ino1][i]+x_no[ino2][i]+x_fa[ifa][i]);
        FOR_I3
          dx_fa[i] += this_area*(v1[i]+v2[i]);
        area_sum += this_area;
      }
      double d2 = 0.0;
      FOR_I3 {
        x_fa_tmp[i] /= 3.0*area_sum;
        dx_fa[i] /= 2.0*area_sum;
        x_fa[ifa][i] += dx_fa[i];
        d2 += dx_fa[i]*dx_fa[i];
      }
      // just set to the x_fa_tmp value...
      FOR_I3
        x_fa[ifa][i] = x_fa_tmp[i];
      if (++iter>200) {
        cerr<<"Error: face centroid did not converge: d2/area_sum: "<<d2/area_sum<<endl;
        throw(-1);
      }
      else if (iter>150) {
        cout<<" > face centroid not converging: "<<ifa<<" "<<iter<<" "<<d2/area_sum<<endl;
      }
      done = (d2/area_sum<1.0E-15);
    }

  }

  // ======================================
  // x_cv, cv_volume in ghost as well...
  // ======================================

  x_cv = new double[ncv_g][3]; // cell center coordinate
  cv_volume = new double[ncv_g]; // its volume


  int count = 0;

  for (int icv = 0; icv<ncv; icv++) {
    // compute an approximate center as the mean of the face x_fa's...
    double x_cv_approx[3];
    for (int i = 0; i<3; i++)
      x_cv_approx[i] = 0.0;
    for (int foc = faocv_i[icv]; foc<faocv_i[icv+1]; foc++) {
      int ifa = faocv_v[foc];
      for (int i = 0; i<3; i++)
        x_cv_approx[i] += x_fa[ifa][i];
    }
    // divide by the number of faces...
    for (int i = 0; i<3; i++)
      x_cv_approx[i] /= (double) (faocv_i[icv+1]-faocv_i[icv]);
    // add tet volumes...
    cv_volume[icv] = 0.0;
    for (int i = 0; i<3; i++)
      x_cv[icv][i] = 0.0;
    // loop on the faces of this cv...
    for (int foc = faocv_i[icv]; foc<faocv_i[icv+1]; foc++) {
      int ifa = faocv_v[foc];
      double v1[3];
      for (int i = 0; i<3; i++)
        v1[i] = x_fa[ifa][i]-x_cv_approx[i];
      // check if the face is inward or outward wrt this cv...
      if (cvofa[ifa][0]==icv) {
        // face is outward - loop through edges in forward direction...
        int ino2 = noofa_v[noofa_i[ifa+1]-1];
        for (int nof = noofa_i[ifa]; nof<noofa_i[ifa+1]; nof++) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i<3; i++) {
            v2[i] = x_no[ino1][i]-x_cv_approx[i];
            v3[i] = x_no[ino2][i]-x_cv_approx[i];
          }
          // 2 nodes, the face, and the approx cv form a tet...
          double this_volume = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])+v1[1]*(v2[2]*v3[0]-v2[0]*v3[2])+v1[2]*(v2[0]*v3[1]-v2[1]
              *v3[0]);
          //          assert(this_volume > 0.0); // check on the grid ordering/right-handedness
          if (this_volume<=0.0) {
            cout<<" > Warning: subtet volume negative: "<<this_volume<<endl;
            count++;
          }
          cv_volume[icv] += this_volume;
          for (int i = 0; i<3; i++)
            x_cv[icv][i] += this_volume*(x_cv_approx[i]+x_fa[ifa][i]+x_no[ino1][i]+x_no[ino2][i]);
        }
      }
      else {
        assert(cvofa[ifa][1]==icv);
        // face is outward, loop through edges in backward direction...
        int ino2 = noofa_v[noofa_i[ifa]];
        for (int nof = noofa_i[ifa+1]-1; nof>=noofa_i[ifa]; nof--) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i<3; i++) {
            v2[i] = x_no[ino1][i]-x_cv_approx[i];
            v3[i] = x_no[ino2][i]-x_cv_approx[i];
          }
          double this_volume = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])+v1[1]*(v2[2]*v3[0]-v2[0]*v3[2])+v1[2]*(v2[0]*v3[1]-v2[1]
              *v3[0]);
          //assert(this_volume > 0.0);
          if (this_volume<=0.0) {
            cout<<" > Warning: subtet volume negative: "<<this_volume<<endl;
            count++;
          }

          cv_volume[icv] += this_volume;
          for (int i = 0; i<3; i++)
            x_cv[icv][i] += this_volume*(x_cv_approx[i]+x_fa[ifa][i]+x_no[ino1][i]+x_no[ino2][i]);
        }
      }
    }
    // normalize both...
    for (int i = 0; i<3; i++)
      x_cv[icv][i] /= 4.0*cv_volume[icv];
    cv_volume[icv] /= 6.0;
    if (cv_volume[icv]<=0.0) cout<<" > Warning: total cell volume negative: "<<cv_volume[icv]<<endl;
  }

  // update volumes and centroids across processors...
  updateCvData(cv_volume, REPLACE_DATA); // blocking exchanges
  updateCvData(x_cv, REPLACE_TRANSLATE_DATA); // "TRANSLATE" is for any periodicity

  // check...
  double my_volume_sum = 0.0;
  for (int icv = 0; icv<ncv; icv++)
    my_volume_sum += cv_volume[icv];
  double volume_sum;
  MPI_Reduce(&my_volume_sum, &volume_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  // also check s dot n...
  double my_buf[2];
  my_buf[0] = 1.0;
  my_buf[1] = 0.0;
  for (int ifa = 0; ifa<nfa_b; ifa++) {
    // the first faces are boundary faces, and should have no icv1...
    assert(cvofa[ifa][1]==-1);
  }
  for (int ifa = nfa_b; ifa<nfa; ifa++) {
    // the remaining faces should have valid cvofa, with a value
    // in the ghost range for the periodic cvs, and in the
    // owned cv range for totally internal cvs...
    assert(cvofa[ifa][1]>=0);
    if (ifa<nfa_bpi) {
      assert((cvofa[ifa][1]>=ncv)&&(cvofa[ifa][1]<ncv_g));
    }
    else {
      assert(cvofa[ifa][1]<ncv);
    }
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    double s[3];
    double smag = 0.0;
    double nmag = 0.0;
    for (int i = 0; i<3; i++) {
      s[i] = x_cv[icv1][i]-x_cv[icv0][i];
      smag += s[i]*s[i];
      nmag += fa_normal[ifa][i]*fa_normal[ifa][i];
    }
    smag = sqrt(smag);
    nmag = sqrt(nmag);
    //cout << "ifa, smag, nmag: " << ifa << " " << smag << " " << nmag << endl;
    double n[3];
    for (int i = 0; i<3; i++) {
      s[i] /= smag;
      n[i] = fa_normal[ifa][i]/nmag;
    }
    double dp = s[0]*n[0]+s[1]*n[1]+s[2]*n[2];
    assert(dp>0.0);
    my_buf[0] = min(my_buf[0], dp);
    my_buf[1] = min(my_buf[1], -dp);
  }

  double buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);

  if (mpi_rank==0) cout<<"calcGeometry(), cv volume: "<<volume_sum<<", s dot n min/max: "<<buf[0]<<" "<<-buf[1]<<endl;

  //MPI_Pause("done");

}

void UgpWithCv::calcGradCoeffHam() {
  if (mpi_rank==0) cout<<" > calcGradCoeffHam()"<<endl;

  // a Green-Gauss gradient reconstruction for cv's
  assert(nbocv_lsg_coeff==NULL);
  nbocv_lsg_coeff = new double[nbocv_s][3];

  // there should be exactly one fake cv for each boundary face...
  int ncv_gf = ncv_g;
  //  assert( ncv_gf - ncv_g == nfa_b );
  assert(boundary_fa_lsg_coeff==NULL);
  boundary_fa_lsg_coeff = new double[nfa_b][3];

  int faono_i_tmp[NNOC_MAX+1];
  const int FAONO_S_MAX = 4*NNOC_MAX; // assume 4 faces per node
  int faono_v_tmp[FAONO_S_MAX];
  double faono_weight_tmp[FAONO_S_MAX];

  FOR_INO
    no_flag[ino] = -1;
  FOR_IFA
    fa_flag[ifa] = -1;
  FOR_ICV {

    int noc_f = noocv_i[icv];
    int noc_l = noocv_i[icv+1]-1;
    for (int noc = noc_f; noc<=noc_l; ++noc) {
      int ino = noocv_v[noc];
      assert(no_flag[ino]==-1);
      no_flag[ino] = noc-noc_f;
    }
    int nnoc = noc_l-noc_f+1;
    assert(nnoc<=NNOC_MAX);

    // at the end of this routine, the grad coeff's for each nbr
    // cv will be stored in faocv_grad. This can accomodate the
    // fact that fake cvs are not cv neighbors...

    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc<=foc_l; ++foc) {
      int ifa = faocv_v[foc];
      assert(fa_flag[ifa]==-1);
      fa_flag[ifa] = foc-foc_f;
    }
    int nfoc = foc_l-foc_f+1;
    assert(nfoc<=NFOC_MAX);

    // build local faono for this cv...

    for (int ii = 0; ii<nnoc; ++ii)
      faono_i_tmp[ii+1] = 0;
    int faono_s_tmp;

    for (int iter = 0; iter<=1; ++iter) {
      for (int foc = foc_f; foc<=foc_l; ++foc) {
        int ifa = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; ++nof) {
          int ino = noofa_v[nof];
          int ii = no_flag[ino];
          assert(ii>=0);
          if (iter==0) {
            faono_i_tmp[ii+1] += 1;
          }
          else {
            faono_v_tmp[faono_i_tmp[ii]] = ifa;
            faono_i_tmp[ii] += 1;
          }
        }
      }
      if (iter==0) {
        faono_i_tmp[0] = 0;
        for (int ii = 0; ii<nnoc; ++ii)
          faono_i_tmp[ii+1] += faono_i_tmp[ii];
        faono_s_tmp = faono_i_tmp[nnoc];
        assert(faono_s_tmp<=FAONO_S_MAX);
      }
      else {
        for (int ii = nnoc; ii>0; --ii)
          faono_i_tmp[ii] = faono_i_tmp[ii-1];
        faono_i_tmp[0] = 0;
      }
    }

    // for nodes with just 2 nbrs, additional nbrs associated with edge-shared nodes must be
    // added...
    int vshift = 0;
    for (int noc = noc_f; noc<=noc_l; ++noc) {
      int ii = noc-noc_f;
      if (faono_i_tmp[ii+1]-faono_i_tmp[ii]==2) {
        int fon_f = faono_i_tmp[ii];
        int ifa0 = faono_v_tmp[fon_f];
        int ifa1 = faono_v_tmp[fon_f+1];
        int ino = noocv_v[noc];
        int ino0 = getPrevNodeOfFaceCCW(ifa0, ino, icv);
        assert(getNextNodeOfFaceCCW(ifa1, ino, icv)==ino0);
        int ii0 = no_flag[ino0];
        vshift += faono_i_tmp[ii0+1]-faono_i_tmp[ii0]-2;
        int ino1 = getNextNodeOfFaceCCW(ifa0, ino, icv);
        assert(getPrevNodeOfFaceCCW(ifa1, ino, icv)==ino1);
        int ii1 = no_flag[ino1];
        vshift += faono_i_tmp[ii1+1]-faono_i_tmp[ii1]-2;
      }
    }
    if (vshift>0) {
      faono_s_tmp += vshift;
      assert(faono_s_tmp<=FAONO_S_MAX);
      // cycle through in reverse...
      for (int noc = noc_l; noc>=noc_f; --noc) {
        int ii = noc-noc_f;
        int fon_f = faono_i_tmp[ii];
        int fon_l = faono_i_tmp[ii+1]-1;
        faono_i_tmp[ii+1] += vshift;
        if (fon_l-fon_f+1>2) {
          // for any node that already has 3 or more nbrs, just shift up...
          for (int fon = fon_l; fon>=fon_f; --fon)
            faono_v_tmp[fon+vshift] = faono_v_tmp[fon];
        }
        else {
          assert(fon_l-fon_f+1==2);
          int ifa0 = faono_v_tmp[fon_f];
          int ifa1 = faono_v_tmp[fon_f+1];
          int ino = noocv_v[noc];
          int ino0 = getPrevNodeOfFaceCCW(ifa0, ino, icv);
          assert(getNextNodeOfFaceCCW(ifa1, ino, icv)==ino0);
          int ii0 = no_flag[ino0];
          int fon0_f = faono_i_tmp[ii0];
          int fon0_l = faono_i_tmp[ii0+1]-1;
          for (int fon0 = fon0_l; fon0>=fon0_f; --fon0) {
            int ifa_nbr = faono_v_tmp[fon0];
            if ((ifa_nbr!=ifa0)&&(ifa_nbr!=ifa1)) {
              faono_v_tmp[fon_l+vshift] = ifa_nbr;
              --vshift;
            }
          }
          int ino1 = getNextNodeOfFaceCCW(ifa0, ino, icv);
          assert(getPrevNodeOfFaceCCW(ifa1, ino, icv)==ino1);
          int ii1 = no_flag[ino1];
          int fon1_f = faono_i_tmp[ii1];
          int fon1_l = faono_i_tmp[ii1+1]-1;
          for (int fon1 = fon1_l; fon1>=fon1_f; --fon1) {
            int ifa_nbr = faono_v_tmp[fon1];
            if ((ifa_nbr!=ifa0)&&(ifa_nbr!=ifa1)) {
              faono_v_tmp[fon_l+vshift] = ifa_nbr;
              --vshift;
            }
          }
          // and add original pair of faces...
          faono_v_tmp[fon_f+1+vshift] = ifa1;
          faono_v_tmp[fon_f+vshift] = ifa0;
        }
      }
      assert(vshift==0);
    }

    // now sort faces and build the geomtric interpolation factors
    // that determine each node's value as a function of surrounding cv's...

    for (int ii = 0; ii<faono_s_tmp; ++ii)
      faono_weight_tmp[ii] = -1000.0;

    for (int noc = noc_f; noc<=noc_l; ++noc) {
      int ino = noocv_v[noc];
      int ii = noc-noc_f;
      int fon_f = faono_i_tmp[ii];
      int fon_l = faono_i_tmp[ii+1]-1;

      switch (fon_l-fon_f+1) {

      case 3: {
        int ifa0 = faono_v_tmp[fon_f];
        int ifa1 = faono_v_tmp[fon_f+1];
        int ifa2 = faono_v_tmp[fon_f+2];
        int ino01 = getNextNodeOfFaceCCW(ifa0, ino, icv);
        if (getPrevNodeOfFaceCCW(ifa1, ino, icv)!=ino01) {
          // flip ifa1 and ifa2...
          faono_v_tmp[fon_f+1] = ifa2;
          faono_v_tmp[fon_f+2] = ifa1;
          ifa1 = faono_v_tmp[fon_f+1];
          ifa2 = faono_v_tmp[fon_f+2];
          assert(getPrevNodeOfFaceCCW(ifa1, ino, icv)==ino01);
        }
        // now should all be in order...
        int ino12 = getNextNodeOfFaceCCW(ifa1, ino, icv);
        assert(ino12==getPrevNodeOfFaceCCW(ifa2, ino, icv));
        int ino20 = getNextNodeOfFaceCCW(ifa2, ino, icv);
        assert(ino20==getPrevNodeOfFaceCCW(ifa0, ino, icv));
        // and now the interpolation values...
        int icv0;
        if (cvofa[ifa0][0]==icv) icv0 = cvofa[ifa0][1];
        else {
          assert(cvofa[ifa0][1]==icv);
          icv0 = cvofa[ifa0][0];
        }
        int icv1;
        if (cvofa[ifa1][0]==icv) icv1 = cvofa[ifa1][1];
        else {
          assert(cvofa[ifa1][1]==icv);
          icv1 = cvofa[ifa1][0];
        }
        int icv2;
        if (cvofa[ifa2][0]==icv) icv2 = cvofa[ifa2][1];
        else {
          assert(cvofa[ifa2][1]==icv);
          icv2 = cvofa[ifa2][0];
        }

        double dx0[3];
        FOR_I3
          dx0[i] = x_cv[icv0][i]-x_cv[icv][i];
        double dx1[3];
        FOR_I3
          dx1[i] = x_cv[icv1][i]-x_cv[icv][i];
        double dx2[3];
        FOR_I3
          dx2[i] = x_cv[icv2][i]-x_cv[icv][i];
        double dx[3];
        FOR_I3
          dx[i] = x_no[ino][i]-x_cv[icv][i];
        double w0 = dx[0]*(dx1[1]*dx2[2]-dx1[2]*dx2[1])+dx[1]*(dx1[2]*dx2[0]-dx1[0]*dx2[2])+dx[2]*(dx1[0]*dx2[1]-dx1[1]
            *dx2[0]);
        double w1 = dx[0]*(dx2[1]*dx0[2]-dx2[2]*dx0[1])+dx[1]*(dx2[2]*dx0[0]-dx2[0]*dx0[2])+dx[2]*(dx2[0]*dx0[1]-dx2[1]
            *dx0[0]);
        double w2 = dx[0]*(dx0[1]*dx1[2]-dx0[2]*dx1[1])+dx[1]*(dx0[2]*dx1[0]-dx0[0]*dx1[2])+dx[2]*(dx0[0]*dx1[1]-dx0[1]
            *dx1[0]);
        double w = dx0[0]*(dx1[1]*dx2[2]-dx1[2]*dx2[1])+dx0[1]*(dx1[2]*dx2[0]-dx1[0]*dx2[2])+dx0[2]*(dx1[0]*dx2[1]
            -dx1[1]*dx2[0]);
        assert(faono_weight_tmp[fon_f]==-1000.0);
        assert(faono_weight_tmp[fon_f+1]==-1000.0);
        assert(faono_weight_tmp[fon_f+2]==-1000.0);
        faono_weight_tmp[fon_f] = w0/w;
        faono_weight_tmp[fon_f+1] = w1/w;
        faono_weight_tmp[fon_f+2] = w2/w;
      }
        break;

      case 4: {
        // here we are going to use a inverse-distance weighted least squares,
        // so ordering doesn't matter...
        int ifa0 = faono_v_tmp[fon_f];
        int ifa1 = faono_v_tmp[fon_f+1];
        int ifa2 = faono_v_tmp[fon_f+2];
        int ifa3 = faono_v_tmp[fon_f+3];

        // and the interpolation values...
        int icv0;
        if (cvofa[ifa0][0]==icv) icv0 = cvofa[ifa0][1];
        else {
          assert(cvofa[ifa0][1]==icv);
          icv0 = cvofa[ifa0][0];
        }
        int icv1;
        if (cvofa[ifa1][0]==icv) icv1 = cvofa[ifa1][1];
        else {
          assert(cvofa[ifa1][1]==icv);
          icv1 = cvofa[ifa1][0];
        }
        int icv2;
        if (cvofa[ifa2][0]==icv) icv2 = cvofa[ifa2][1];
        else {
          assert(cvofa[ifa2][1]==icv);
          icv2 = cvofa[ifa2][0];
        }
        int icv3;
        if (cvofa[ifa3][0]==icv) icv3 = cvofa[ifa3][1];
        else {
          assert(cvofa[ifa3][1]==icv);
          icv3 = cvofa[ifa3][0];
        }
        // we need to get weights for the nodal value in terms of the
        // surrounding cv values. Weight these by inverse distance from the node...
        double w0 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv0][i]-x_no[ino][i];
          w0 += dx*dx;
        }
        w0 = 1.0/sqrt(w0);
        double w1 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv1][i]-x_no[ino][i];
          w1 += dx*dx;
        }
        w1 = 1.0/sqrt(w1);
        double w2 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv2][i]-x_no[ino][i];
          w2 += dx*dx;
        }
        w2 = 1.0/sqrt(w2);
        double w3 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv3][i]-x_no[ino][i];
          w3 += dx*dx;
        }
        w3 = 1.0/sqrt(w3);
        // dx's are relative to the icv, where the polynomial is constrained...
        double dx0[3];
        FOR_I3
          dx0[i] = x_cv[icv0][i]-x_cv[icv][i];
        double dx1[3];
        FOR_I3
          dx1[i] = x_cv[icv1][i]-x_cv[icv][i];
        double dx2[3];
        FOR_I3
          dx2[i] = x_cv[icv2][i]-x_cv[icv][i];
        double dx3[3];
        FOR_I3
          dx3[i] = x_cv[icv3][i]-x_cv[icv][i];
        // sums...
        double swxx = w0*dx0[0]*dx0[0]+w1*dx1[0]*dx1[0]+w2*dx2[0]*dx2[0]+w3*dx3[0]*dx3[0];
        double swxy = w0*dx0[0]*dx0[1]+w1*dx1[0]*dx1[1]+w2*dx2[0]*dx2[1]+w3*dx3[0]*dx3[1];
        double swzx = w0*dx0[0]*dx0[2]+w1*dx1[0]*dx1[2]+w2*dx2[0]*dx2[2]+w3*dx3[0]*dx3[2];
        double swyy = w0*dx0[1]*dx0[1]+w1*dx1[1]*dx1[1]+w2*dx2[1]*dx2[1]+w3*dx3[1]*dx3[1];
        double swyz = w0*dx0[1]*dx0[2]+w1*dx1[1]*dx1[2]+w2*dx2[1]*dx2[2]+w3*dx3[1]*dx3[2];
        double swzz = w0*dx0[2]*dx0[2]+w1*dx1[2]*dx1[2]+w2*dx2[2]*dx2[2]+w3*dx3[2]*dx3[2];
        // denom...
        double denom = swxy*swxy*swzz+swzx*swzx*swyy+swyz*swyz*swxx-swxx*swyy*swzz-2.0*swxy*swyz*swzx;
        double dx[3];
        FOR_I3
          dx[i] = x_no[ino][i]-x_cv[icv][i];
        double cx = (dx[0]*(swyz*swyz-swyy*swzz)+dx[1]*(swxy*swzz-swzx*swyz)+dx[2]*(swzx*swyy-swxy*swyz))/denom;
        double cy = (dx[1]*(swzx*swzx-swzz*swxx)+dx[2]*(swyz*swxx-swxy*swzx)+dx[0]*(swxy*swzz-swyz*swzx))/denom;
        double cz = (dx[2]*(swxy*swxy-swxx*swyy)+dx[0]*(swzx*swyy-swyz*swxy)+dx[1]*(swyz*swxx-swzx*swxy))/denom;
        w0 *= cx*dx0[0]+cy*dx0[1]+cz*dx0[2];
        w1 *= cx*dx1[0]+cy*dx1[1]+cz*dx1[2];
        w2 *= cx*dx2[0]+cy*dx2[1]+cz*dx2[2];
        w3 *= cx*dx3[0]+cy*dx3[1]+cz*dx3[2];
        assert(faono_weight_tmp[fon_f]==-1000.0);
        assert(faono_weight_tmp[fon_f+1]==-1000.0);
        assert(faono_weight_tmp[fon_f+2]==-1000.0);
        assert(faono_weight_tmp[fon_f+3]==-1000.0);
        faono_weight_tmp[fon_f] = w0;
        faono_weight_tmp[fon_f+1] = w1;
        faono_weight_tmp[fon_f+2] = w2;
        faono_weight_tmp[fon_f+3] = w3;
      }
        break;

      case 5: {
        // here we are going to use a inverse-distance weighted least squares,
        // so ordering doesn't matter. XXXX In the future, we should use a version
        // of this reconstruction that constrains the reconstructions associated
        // with shared faces, but for now, just use inverse distance and hope for the
        // best...
        int ifa0 = faono_v_tmp[fon_f];
        int ifa1 = faono_v_tmp[fon_f+1];
        int ifa2 = faono_v_tmp[fon_f+2];
        int ifa3 = faono_v_tmp[fon_f+3];
        int ifa4 = faono_v_tmp[fon_f+4];

        // and the interpolation values...
        int icv0;
        if (cvofa[ifa0][0]==icv) icv0 = cvofa[ifa0][1];
        else {
          assert(cvofa[ifa0][1]==icv);
          icv0 = cvofa[ifa0][0];
        }
        int icv1;
        if (cvofa[ifa1][0]==icv) icv1 = cvofa[ifa1][1];
        else {
          assert(cvofa[ifa1][1]==icv);
          icv1 = cvofa[ifa1][0];
        }
        int icv2;
        if (cvofa[ifa2][0]==icv) icv2 = cvofa[ifa2][1];
        else {
          assert(cvofa[ifa2][1]==icv);
          icv2 = cvofa[ifa2][0];
        }
        int icv3;
        if (cvofa[ifa3][0]==icv) icv3 = cvofa[ifa3][1];
        else {
          assert(cvofa[ifa3][1]==icv);
          icv3 = cvofa[ifa3][0];
        }
        int icv4;
        if (cvofa[ifa4][0]==icv) icv4 = cvofa[ifa4][1];
        else {
          assert(cvofa[ifa4][1]==icv);
          icv4 = cvofa[ifa4][0];
        }
        // we need to get weights for the nodal value in terms of the
        // surrounding cv values. Weight these by inverse distance from the node...
        double w0 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv0][i]-x_no[ino][i];
          w0 += dx*dx;
        }
        w0 = 1.0/sqrt(w0);
        double w1 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv1][i]-x_no[ino][i];
          w1 += dx*dx;
        }
        w1 = 1.0/sqrt(w1);
        double w2 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv2][i]-x_no[ino][i];
          w2 += dx*dx;
        }
        w2 = 1.0/sqrt(w2);
        double w3 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv3][i]-x_no[ino][i];
          w3 += dx*dx;
        }
        w3 = 1.0/sqrt(w3);
        double w4 = 0.0;
        FOR_I3 {
          double dx = x_cv[icv4][i]-x_no[ino][i];
          w4 += dx*dx;
        }
        w4 = 1.0/sqrt(w4);
        // dx's are relative to the icv, where the polynomial is constrained...
        double dx0[3];
        FOR_I3
          dx0[i] = x_cv[icv0][i]-x_cv[icv][i];
        double dx1[3];
        FOR_I3
          dx1[i] = x_cv[icv1][i]-x_cv[icv][i];
        double dx2[3];
        FOR_I3
          dx2[i] = x_cv[icv2][i]-x_cv[icv][i];
        double dx3[3];
        FOR_I3
          dx3[i] = x_cv[icv3][i]-x_cv[icv][i];
        double dx4[3];
        FOR_I3
          dx4[i] = x_cv[icv4][i]-x_cv[icv][i];
        // sums...
        double swxx = w0*dx0[0]*dx0[0]+w1*dx1[0]*dx1[0]+w2*dx2[0]*dx2[0]+w3*dx3[0]*dx3[0]+w4*dx4[0]*dx4[0];
        double swxy = w0*dx0[0]*dx0[1]+w1*dx1[0]*dx1[1]+w2*dx2[0]*dx2[1]+w3*dx3[0]*dx3[1]+w4*dx4[0]*dx4[1];
        double swzx = w0*dx0[0]*dx0[2]+w1*dx1[0]*dx1[2]+w2*dx2[0]*dx2[2]+w3*dx3[0]*dx3[2]+w4*dx4[0]*dx4[2];
        double swyy = w0*dx0[1]*dx0[1]+w1*dx1[1]*dx1[1]+w2*dx2[1]*dx2[1]+w3*dx3[1]*dx3[1]+w4*dx4[1]*dx4[1];
        double swyz = w0*dx0[1]*dx0[2]+w1*dx1[1]*dx1[2]+w2*dx2[1]*dx2[2]+w3*dx3[1]*dx3[2]+w4*dx4[1]*dx4[2];
        double swzz = w0*dx0[2]*dx0[2]+w1*dx1[2]*dx1[2]+w2*dx2[2]*dx2[2]+w3*dx3[2]*dx3[2]+w4*dx4[2]*dx4[2];
        // denom...
        double denom = swxy*swxy*swzz+swzx*swzx*swyy+swyz*swyz*swxx-swxx*swyy*swzz-2.0*swxy*swyz*swzx;
        double dx[3];
        FOR_I3
          dx[i] = x_no[ino][i]-x_cv[icv][i];
        double cx = (dx[0]*(swyz*swyz-swyy*swzz)+dx[1]*(swxy*swzz-swzx*swyz)+dx[2]*(swzx*swyy-swxy*swyz))/denom;
        double cy = (dx[1]*(swzx*swzx-swzz*swxx)+dx[2]*(swyz*swxx-swxy*swzx)+dx[0]*(swxy*swzz-swyz*swzx))/denom;
        double cz = (dx[2]*(swxy*swxy-swxx*swyy)+dx[0]*(swzx*swyy-swyz*swxy)+dx[1]*(swyz*swxx-swzx*swxy))/denom;
        w0 *= cx*dx0[0]+cy*dx0[1]+cz*dx0[2];
        w1 *= cx*dx1[0]+cy*dx1[1]+cz*dx1[2];
        w2 *= cx*dx2[0]+cy*dx2[1]+cz*dx2[2];
        w3 *= cx*dx3[0]+cy*dx3[1]+cz*dx3[2];
        w4 *= cx*dx4[0]+cy*dx4[1]+cz*dx4[2];
        assert(faono_weight_tmp[fon_f]==-1000.0);
        assert(faono_weight_tmp[fon_f+1]==-1000.0);
        assert(faono_weight_tmp[fon_f+2]==-1000.0);
        assert(faono_weight_tmp[fon_f+3]==-1000.0);
        assert(faono_weight_tmp[fon_f+4]==-1000.0);
        faono_weight_tmp[fon_f] = w0;
        faono_weight_tmp[fon_f+1] = w1;
        faono_weight_tmp[fon_f+2] = w2;
        faono_weight_tmp[fon_f+3] = w3;
        faono_weight_tmp[fon_f+4] = w4;
      }
        break;

      default:
        cerr<<"Error: strange number of faces per node: "<<fon_l-fon_f+1<<endl;
        throw(-1);
      }
    }

    // ====================================================================
    // we now have the interpolation weights associated with each node...
    // now loop through the faces and build the gradient...
    // ====================================================================

    double faocv_grad[NFOC_MAX][3];
    for (int ii = 0; ii<nfoc; ++ii)
      FOR_I3
        faocv_grad[ii][i] = 0.0;
    double cv_vol = 0.0;

    for (int foc = foc_f; foc<=foc_l; ++foc) {

      int ifa = faocv_v[foc];

      double fa_area = sqrt(fa_normal[ifa][0]*fa_normal[ifa][0]+fa_normal[ifa][1]*fa_normal[ifa][1]+fa_normal[ifa][2]
          *fa_normal[ifa][2]);
      double unit_normal[3];
      FOR_I3
        unit_normal[i] = fa_normal[ifa][i]/fa_area;

      // loop thru face nodes in CCW order...
      int nof_f = noofa_i[ifa];
      int nof_begin, nof_end, nof_inc, nnof;
      if (cvofa[ifa][0]==icv) {
        nof_begin = noofa_i[ifa];
        nof_end = noofa_i[ifa+1];
        nof_inc = 1;
        nnof = nof_end-nof_begin;
      }
      else {
        assert(cvofa[ifa][1]==icv);
        nof_begin = noofa_i[ifa+1]-1;
        nof_end = noofa_i[ifa]-1;
        nof_inc = -1;
        nnof = nof_begin-nof_end;
        FOR_I3
          unit_normal[i] = -unit_normal[i];
      }

      // compute the set of of weights required to interpolate to the
      // face centroid. Since the centroid is already calculated, just use the
      // following formula....
      //
      // (sum_{sub-tri}(mag))*x_fa = sum_{sub-tri}(mag*(x_fa+x_no0+x_no1)/3);
      //
      // and solve for x_fa.

      double noofa_weight[NNOF_MAX];
      for (int ii = 0; ii<nnof; ++ii)
        noofa_weight[ii] = 0.0;
      double noofa_weight_sum = 0.0;

      int nof1 = nof_end-nof_inc;
      int ino1 = noofa_v[nof1];
      for (int nof = nof_begin; nof!=nof_end; nof += nof_inc) {
        int nof0 = nof1;
        int ino0 = ino1;
        nof1 = nof;
        ino1 = noofa_v[nof1];
        double dx0[3];
        FOR_I3
          dx0[i] = x_no[ino0][i]-x_fa[ifa][i];
        double dx1[3];
        FOR_I3
          dx1[i] = x_no[ino1][i]-x_fa[ifa][i];
        double this_normal[3];
        this_normal[0] = dx0[1]*dx1[2]-dx0[2]*dx1[1];
        this_normal[1] = dx0[2]*dx1[0]-dx0[0]*dx1[2];
        this_normal[2] = dx0[0]*dx1[1]-dx0[1]*dx1[0];
        // note this is the fluent definition (I think) that always produces
        // positive weights, although the convergence is potentially very slow
        /*
         double mag = sqrt( this_normal[0]*this_normal[0] +
         this_normal[1]*this_normal[1] +
         this_normal[2]*this_normal[2] );
         */
        // dotting on the uniquely defined unit_normal (does not depend on the
        // face center/centroid) results in faster convergence, but some of these
        // weights could be negative if the face-handedness is right on the edge. The
        // only important thing is that this reconstruction is the same as that used in
        // the calcGeometry routine...
        double mag = this_normal[0]*unit_normal[0]+this_normal[1]*unit_normal[1]+this_normal[2]*unit_normal[2];
        noofa_weight[nof0-nof_f] += mag;
        noofa_weight[nof1-nof_f] += mag;
        noofa_weight_sum += mag*2.0;
      }
      // normalize and check...
      double x_fa_check[3] = { 0.0, 0.0, 0.0 };
      for (int nof = nof_begin; nof!=nof_end; nof += nof_inc) {
        noofa_weight[nof-nof_f] /= noofa_weight_sum;
        int ino = noofa_v[nof];
        FOR_I3
          x_fa_check[i] += noofa_weight[nof-nof_f]*x_no[ino][i];
      }
      double d2 = 0.0;
      FOR_I3 {
        double dx = x_fa_check[i]-x_fa[ifa][i];
        d2 += dx*dx;
      }
      //cout << "dist: " << sqrt(d2) << endl;
      assert(sqrt(d2)<1.0E-10*sqrt(noofa_weight_sum)); // noofa_weight_sum is proportional to the face area

      // now assemble gradient...
      nof1 = nof_end-nof_inc;
      ino1 = noofa_v[nof1];
      for (int nof = nof_begin; nof!=nof_end; nof += nof_inc) {
        int nof0 = nof1;
        int ino0 = ino1;
        nof1 = nof;
        ino1 = noofa_v[nof1];
        double dx0[3];
        FOR_I3
          dx0[i] = x_no[ino0][i]-x_fa[ifa][i];
        double dx1[3];
        FOR_I3
          dx1[i] = x_no[ino1][i]-x_fa[ifa][i];
        double this_normal[3];
        this_normal[0] = dx0[1]*dx1[2]-dx0[2]*dx1[1];
        this_normal[1] = dx0[2]*dx1[0]-dx0[0]*dx1[2];
        this_normal[2] = dx0[0]*dx1[1]-dx0[1]*dx1[0];
        // this is twice the outward pointing normal between x_fa, x_no0 and x_no1...
        // ino0...
        int fon_f = faono_i_tmp[no_flag[ino0]];
        int fon_l = faono_i_tmp[no_flag[ino0]+1]-1;
        for (int fon = fon_f; fon<=fon_l; ++fon) {
          int ifa_nbr = faono_v_tmp[fon];
          double this_weight = faono_weight_tmp[fon];
          FOR_I3
            faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
        }
        // ino1...
        fon_f = faono_i_tmp[no_flag[ino1]];
        fon_l = faono_i_tmp[no_flag[ino1]+1]-1;
        for (int fon = fon_f; fon<=fon_l; ++fon) {
          int ifa_nbr = faono_v_tmp[fon];
          double this_weight = faono_weight_tmp[fon];
          FOR_I3
            faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
        }
        // face...
        for (int nof2 = nof_begin; nof2!=nof_end; nof2 += nof_inc) {
          int ino2 = noofa_v[nof2];
          fon_f = faono_i_tmp[no_flag[ino2]];
          fon_l = faono_i_tmp[no_flag[ino2]+1]-1;
          for (int fon = fon_f; fon<=fon_l; ++fon) {
            int ifa_nbr = faono_v_tmp[fon];
            double this_weight = faono_weight_tmp[fon]*noofa_weight[nof2-nof_f];
            FOR_I3
              faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
          }
        }
        // to get the sub-tet volume, dot with the vector from the
        // cell centroid to the face centroid
        double this_vol = this_normal[0]*(x_fa[ifa][0]-x_cv[icv][0])+this_normal[1]*(x_fa[ifa][1]-x_cv[icv][1])
            +this_normal[2]*(x_fa[ifa][2]-x_cv[icv][2]);
        assert(this_vol>0.0);
        cv_vol += this_vol;

      }
    }

    assert(fabs(cv_vol-6.0*cv_volume[icv])<1.0E-10*cv_volume[icv]);

    double diag[3] = { 0.0, 0.0, 0.0 };
    int nboc_f = nbocv_i[icv]; // diagonal...
    for (int foc = foc_f; foc<=foc_l; ++foc) {
      int ifa = faocv_v[foc];
      FOR_I3
        diag[i] -= faocv_grad[foc-foc_f][i]/cv_vol;
      if (ifa<nfa_b) {
        // boundary face / fake cv...
        assert(cvofa[ifa][0]==icv);
        FOR_I3
          boundary_fa_lsg_coeff[ifa][i] = faocv_grad[foc-foc_f][i]/cv_vol;
      }
      else {
        // internal or inter-processor face / internal or ghost cv
        if (cvofa[ifa][0]==icv) {
          int nboc = nboc_f; // first is the diagonal...
          while (nbocv_v[++nboc]!=cvofa[ifa][1])
            assert(nboc<nbocv_i[icv+1]);
          FOR_I3
            nbocv_lsg_coeff[nboc][i] = faocv_grad[foc-foc_f][i]/cv_vol;
        }
        else {
          assert(cvofa[ifa][1]==icv);
          int nboc = nboc_f; // first is the diagonal...
          while (nbocv_v[++nboc]!=cvofa[ifa][0])
            assert(nboc<nbocv_i[icv+1]);
          FOR_I3
            nbocv_lsg_coeff[nboc][i] = faocv_grad[foc-foc_f][i]/cv_vol;
        }
      }
    }
    FOR_I3
      nbocv_lsg_coeff[nboc_f][i] = diag[i];

    // cleanup...

    for (int noc = noc_f; noc<=noc_l; ++noc) {
      int ino = noocv_v[noc];
      no_flag[ino] = -1;
    }
    for (int foc = foc_f; foc<=foc_l; ++foc) {
      int ifa = faocv_v[foc];
      fa_flag[ifa] = -1;
    }

  }

  //MPI_Pause("done Ham grad");
  // new
  // delete
}

void UgpWithCv::calcLsgCoeff() {

  if (mpi_rank==0) cout<<"calcLsgCoeff()"<<endl;

  // least squares gradient coeff's...
  assert(nbocv_lsg_coeff==NULL);
  nbocv_lsg_coeff = new double[nbocv_s][3];
  for (int noc = 0; noc<nbocv_s; noc++)
    for (int i = 0; i<3; i++)
      nbocv_lsg_coeff[noc][i] = 0.0;

  assert(boundary_fa_lsg_coeff==NULL);
  boundary_fa_lsg_coeff = new double[nfa_b][3];
  for (int ifa = 0; ifa<nfa_b; ifa++)
    for (int i = 0; i<3; i++)
      boundary_fa_lsg_coeff[ifa][i] = 0.0;

  for (int icv = 0; icv<ncv; icv++) {
    double swdx2 = 0.0; // sum weight dx2, etc...
    double swdy2 = 0.0;
    double swdz2 = 0.0;
    double swdxdy = 0.0;
    double swdxdz = 0.0;
    double swdydz = 0.0;
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    // skip the diagonal in this loop...
    for (int noc = noc_f+1; noc<=noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      double weight = 1.0;
      double dx[3];
      for (int i = 0; i<3; i++)
        dx[i] = x_cv[icv_nbr][i]-x_cv[icv][i];

      weight = 1.0/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

      swdx2 += weight*dx[0]*dx[0];
      swdy2 += weight*dx[1]*dx[1];
      swdz2 += weight*dx[2]*dx[2];
      swdxdy += weight*dx[0]*dx[1];
      swdxdz += weight*dx[0]*dx[2];
      swdydz += weight*dx[1]*dx[2];
    }

    // this cell may also touch some boundary faces that
    // do not have neighbours (no fake cells). Their
    // influence on the gradient has to be considered as well...
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc<=foc_l; foc++) {
      int ifa = faocv_v[foc];
      if (ifa<nfa_b) {
        // this is a bounray face...
        assert(cvofa[ifa][0]==icv);
        assert(cvofa[ifa][1]==-1);
        double weight = 1.0;
        double dx[3];
        for (int i = 0; i<3; i++)
          dx[i] = x_fa[ifa][i]-x_cv[icv][i];

        weight = 1.0/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

        swdx2 += weight*dx[0]*dx[0];
        swdy2 += weight*dx[1]*dx[1];
        swdz2 += weight*dx[2]*dx[2];
        swdxdy += weight*dx[0]*dx[1];
        swdxdz += weight*dx[0]*dx[2];
        swdydz += weight*dx[1]*dx[2];
      }
    }

    double denom = 2.0*swdxdy*swdxdz*swdydz+swdx2*swdy2*swdz2-swdx2*swdydz*swdydz-swdy2*swdxdz*swdxdz-swdz2*swdxdy
        *swdxdy;

    // check non-singular...
    assert(fabs(denom)>1.0E-12*cv_volume[icv]*cv_volume[icv]);

    // assemble coeff's...
    for (int noc = noc_f+1; noc<=noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      double weight = 1.0;
      double dx[3];
      for (int i = 0; i<3; i++)
        dx[i] = x_cv[icv_nbr][i]-x_cv[icv][i];

      weight = 1.0/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

      // following coeff's multiply the delta across nbr - cv...
      nbocv_lsg_coeff[noc][0] = weight*((swdy2*swdz2-swdydz*swdydz)*dx[0]+(swdxdz*swdydz-swdxdy*swdz2)*dx[1]+(swdxdy
          *swdydz-swdxdz*swdy2)*dx[2])/denom;
      nbocv_lsg_coeff[noc][1] = weight*((swdxdz*swdydz-swdxdy*swdz2)*dx[0]+(swdx2*swdz2-swdxdz*swdxdz)*dx[1]+(swdxdy
          *swdxdz-swdydz*swdx2)*dx[2])/denom;
      nbocv_lsg_coeff[noc][2] = weight*((swdxdy*swdydz-swdxdz*swdy2)*dx[0]+(swdxdy*swdxdz-swdydz*swdx2)*dx[1]+(swdx2
          *swdy2-swdxdy*swdxdy)*dx[2])/denom;
      // and subtract from the diagonal...
      for (int i = 0; i<3; i++)
        nbocv_lsg_coeff[noc_f][i] -= nbocv_lsg_coeff[noc][i];
    }

    for (int foc = foc_f; foc<=foc_l; foc++) {
      int ifa = faocv_v[foc];
      if (ifa<nfa_b) {
        // this is a boundary face...
        assert(cvofa[ifa][0]==icv);
        assert(cvofa[ifa][1]==-1);
        double weight = 1.0;
        double dx[3];
        for (int i = 0; i<3; i++)
          dx[i] = x_fa[ifa][i]-x_cv[icv][i];

        weight = 1.0/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

        boundary_fa_lsg_coeff[ifa][0] = weight*((swdy2*swdz2-swdydz*swdydz)*dx[0]+(swdxdz*swdydz-swdxdy*swdz2)*dx[1]
            +(swdxdy*swdydz-swdxdz*swdy2)*dx[2])/denom;
        boundary_fa_lsg_coeff[ifa][1] = weight*((swdxdz*swdydz-swdxdy*swdz2)*dx[0]+(swdx2*swdz2-swdxdz*swdxdz)*dx[1]
            +(swdxdy*swdxdz-swdydz*swdx2)*dx[2])/denom;
        boundary_fa_lsg_coeff[ifa][2] = weight*((swdxdy*swdydz-swdxdz*swdy2)*dx[0]+(swdxdy*swdxdz-swdydz*swdx2)*dx[1]
            +(swdx2*swdy2-swdxdy*swdxdy)*dx[2])/denom;

        // and subtract from the diagonal for Dirichlet boundary condition
        for (int i = 0; i<3; i++)
          nbocv_lsg_coeff[noc_f][i] -= boundary_fa_lsg_coeff[ifa][i];

        // do not subtract from the diagonal, so that the default
        // gradient boundary condition is d/ds = 0, figure
        // out how to make this d/dn = 0 later.
      }
    }
  }

}

void UgpWithCv::checkLsgCoeff() {

  if (mpi_rank==0) cout<<"checkLsgCoeff()...";

  double * p = new double[ncv_g];
  double * p_bc = new double[nfa_b];

  double constant_grad_p[3];
  constant_grad_p[0] = 1.12;
  constant_grad_p[1] = 2.13;
  constant_grad_p[2] = 3.14;

  double (*grad_p)[3] = new double[ncv_g][3];

  // put a linear grad in p - everywhere...
  for (int icv = 0; icv<ncv_g; icv++)
    p[icv] = 1.1+1.12*x_cv[icv][0]+2.13*x_cv[icv][1]+3.14*x_cv[icv][2];

  // including the boundary...
  for (int ifa = 0; ifa<nfa_b; ifa++)
    p_bc[ifa] = 1.1+1.12*x_fa[ifa][0]+2.13*x_fa[ifa][1]+3.14*x_fa[ifa][2];

  // compute gradients...
  for (int icv = 0; icv<ncv; icv++) {
    for (int i = 0; i<3; i++)
      grad_p[icv][i] = 0.0;
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    for (int noc = noc_f; noc<=noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      for (int i = 0; i<3; i++)
        grad_p[icv][i] += nbocv_lsg_coeff[noc][i]*p[icv_nbr];
    }
  }

  // for boundary closures, one of the following is possible...

  // if bc is dirichlet...

  for (int ifa = 0; ifa<nfa_b; ifa++) {
    int icv = cvofa[ifa][0];
    for (int i = 0; i<3; i++)
      grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*p_bc[ifa];
  }

  /**/

  // if you know the full gradient...
  /*
   for (int ifa = 0; ifa < nfa_b; ifa++) {
   int icv = cvofa[ifa][0];
   double grad_p_ds = 0.0;
   for (int i = 0; i < 3; i++) {
   double ds = x_fa[ifa][i] - x_cv[icv][i];
   grad_p_ds += ds*constant_grad_p[i];
   }
   for (int i = 0; i < 3; i++)
   grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*grad_p_ds;
   }
   */

  // if you know the normal component of the gradient...
  /*  double (*coeff_b)[3][3] = new double[ncv_b][3][3];

   // zero the coeff_b matrix...
   for (int icv = 0; icv < ncv_b; icv++) {
   for (int i = 0; i < 3; i++) {
   for (int j = 0; j < 3; j++) {
   coeff_b[icv][i][j] = 0.0;
   }
   }
   }

   for (int ifa = 0; ifa < nfa_b; ifa++) {
   int icv = cvofa[ifa][0];
   // this cv should be in the boundary cv's...
   assert((icv >= 0) && (icv < ncv_b));
   double n[3], s[3];
   double nmag = 0.0;
   double smag = 0.0;
   for (int i = 0; i < 3; i++) {
   n[i] = fa_normal[ifa][i];
   nmag += n[i] * n[i];
   s[i] = x_fa[ifa][i] - x_cv[icv][i];
   smag += s[i] * s[i];
   }
   nmag = sqrt(nmag);
   smag = sqrt(smag);
   double alpha = 0.0;
   for (int i = 0; i < 3; i++) {
   n[i] /= nmag;
   s[i] /= smag;
   alpha += n[i] * s[i];
   }
   assert((alpha > 0.0) && (alpha < 1.00000001));
   // the exact dpdn is then ...
   double dpdn = constant_grad_p[0] * n[0] + constant_grad_p[1] * n[1] + constant_grad_p[2] * n[2];
   // similar to above, we want to write the face-based value in terms of
   // the known dpdn and the implicit cell based full dpdx_i...
   // note that (see above):
   // grad_p_ds = grad_p[j].dot.s[j]*smag
   //           = alpha*grad_p[j].dot.n[j]*smag + grad_p[j].dot.(s[j]-alpha*n[j])*smag;
   //           = alpha*dpdn*smag + grad_p[j].dot.(n[j]-alpha*s[j])*smag;
   // now put the first part into grad_p, and the second into the 3x3 grad coeff matrix...
   for (int i = 0; i < 3; i++) {
   grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i] * alpha * dpdn * smag;
   for (int j = 0; j < 3; j++) {
   coeff_b[icv][i][j] -= boundary_fa_lsg_coeff[ifa][i] * (s[j] - alpha * n[j]) * smag;
   }
   }
   }

   // gradients in the boundary cells may have to be modified...
   for (int icv = 0; icv < ncv_b; icv++) {
   // include the unit diagonal here...
   for (int i = 0; i < 3; i++)
   coeff_b[icv][i][i] += 1.0;

   double denom = coeff_b[icv][0][0] * coeff_b[icv][1][1] * coeff_b[icv][2][2]
   + coeff_b[icv][0][2] * coeff_b[icv][2][1] * coeff_b[icv][1][0]
   + coeff_b[icv][2][0] * coeff_b[icv][0][1] * coeff_b[icv][1][2]
   - coeff_b[icv][0][0] * coeff_b[icv][1][2] * coeff_b[icv][2][1]
   - coeff_b[icv][1][1] * coeff_b[icv][0][2] * coeff_b[icv][2][0]
   - coeff_b[icv][2][2] * coeff_b[icv][0][1] * coeff_b[icv][1][0];

   double rhs0 = grad_p[icv][0];
   double rhs1 = grad_p[icv][1];
   double rhs2 = grad_p[icv][2];

   grad_p[icv][0] = ((coeff_b[icv][1][1] * coeff_b[icv][2][2] - coeff_b[icv][1][2] * coeff_b[icv][2][1]) * rhs0
   + (coeff_b[icv][2][1] * coeff_b[icv][0][2] - coeff_b[icv][2][2] * coeff_b[icv][0][1]) * rhs1
   + (coeff_b[icv][0][1] * coeff_b[icv][1][2] - coeff_b[icv][0][2] * coeff_b[icv][1][1]) * rhs2) / denom;
   grad_p[icv][1] = ((coeff_b[icv][1][2] * coeff_b[icv][2][0] - coeff_b[icv][1][0] * coeff_b[icv][2][2]) * rhs0
   + (coeff_b[icv][2][2] * coeff_b[icv][0][0] - coeff_b[icv][2][0] * coeff_b[icv][0][2]) * rhs1
   + (coeff_b[icv][0][2] * coeff_b[icv][1][0] - coeff_b[icv][0][0] * coeff_b[icv][1][2]) * rhs2) / denom;
   grad_p[icv][2] = ((coeff_b[icv][1][0] * coeff_b[icv][2][1] - coeff_b[icv][1][1] * coeff_b[icv][2][0]) * rhs0
   + (coeff_b[icv][2][0] * coeff_b[icv][0][1] - coeff_b[icv][2][1] * coeff_b[icv][0][0]) * rhs1
   + (coeff_b[icv][0][0] * coeff_b[icv][1][1] - coeff_b[icv][0][1] * coeff_b[icv][1][0]) * rhs2) / denom;
   }*/

  for (int icv = 0; icv<ncv; icv++) {
    if ((fabs(grad_p[icv][0]-1.12)>1.0E-8)||(fabs(grad_p[icv][1]-2.13)>1.0E-8)||(fabs(grad_p[icv][2]-3.14)>1.0E-8)) {
      cerr<<"Error: linear grad test failed: "<<grad_p[icv][0]-1.12<<" "<<grad_p[icv][1]-2.13<<" "<<grad_p[icv][2]-3.14
          <<endl;
      throw(-1);
    }
  }

  if (mpi_rank==0) cout<<"OK"<<endl;

  delete[] p;
  delete[] p_bc;
  delete[] grad_p;
  //  delete[] coeff_b;

}

int UgpWithCv::solveCvScalarJacobi(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double relax,
    const double zero, const int maxiter) {

  // need gradients in ghost as well...
  double (*grad_phi)[3] = new double[ncv_g][3];

  // res only in cvs we own...
  double * res = new double[ncv];

  int iter = 0;
  int done = 0;
  while (done==0) {

    iter++;

    // compute [Ap]{phi} + [Ap_grad]{grad_phi}...

    // start with {grad_phi}...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      for (int i = 0; i<3; i++)
        grad_phi[icv][i] = nbocv_lsg_coeff[noc_f][i]*phi[icv];
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i<3; i++)
          grad_phi[icv][i] += nbocv_lsg_coeff[noc][i]*phi[icv_nbr];
      }
    }

    // also need the gradient in the ghosts...
    updateCvData(grad_phi, REPLACE_ROTATE_DATA);

    // {res} = [Ap]{phi} + [Ap_grad]{grad_phi}
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      res[icv] = Ap[noc_f]*phi[icv]; // diagonal
      for (int i = 0; i<3; i++)
        res[icv] += Ap_grad[noc_f][i]*grad_phi[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        res[icv] += Ap[noc]*phi[icv_nbr];
        for (int i = 0; i<3; i++)
          res[icv] += Ap_grad[noc][i]*grad_phi[icv_nbr][i];
      }
    }

    // update phi and compute max (Linf) normalized residual...
    double my_res_max = 0.0;
    for (int icv = 0; icv<ncv; icv++) {
      double this_res = (rhs[icv]-res[icv])/Ap[nbocv_i[icv]];
      my_res_max = max(my_res_max, fabs(this_res));
      // update phi with relaxation...
      phi[icv] += relax*this_res;
    }

    // update phi in ghost cv's of all processors...
    updateCvData(phi, REPLACE_DATA);

    // are we done? decide on one processor...
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank==0) {
      //if (iter%10 == 0)
      //cout << "jacobi iter, res_max: " << iter << " " << res_max << endl;
      if (res_max<zero) {
        done = 1;
      }
      else if (iter>maxiter) {
        cout<<"Warning: solveCvScalarJacobi did not converge after "<<maxiter<<" iters, res_max: "<<res_max<<endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

  } // while (done == 0)

  // cleanup...
  delete[] grad_phi;
  delete[] res;

  // let the calling routine know if we were successful...
  return (done==1);

}

int UgpWithCv::solveCvScalarCg(double * phi, double * Ap, double * rhs, const int mode, const double zero,
    const int maxiter) {

  // we need the following work arrays...
  double * res = new double[ncv];
  double * v = new double[ncv];
  double * p = new double[ncv_g];

  // initialize...
  for (int icv = 0; icv<ncv; icv++)
    p[icv] = 0.0;

  double rho = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv<ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    res[icv] = Ap[noc_f]*phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f+1; noc<=noc_l; noc++)
      res[icv] += Ap[noc]*phi[nbocv_v[noc]];
    res[icv] = rhs[icv]-res[icv];
  }

  double res0_max;
  if (mode==RELATIVE_RESIDUAL) {
    double my_res0_max = 0.0;

    for (int icv = 0; icv<ncv; icv++)
      my_res0_max = max(my_res0_max, fabs(res[icv]/Ap[nbocv_i[icv]]));

    MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  }

  int iter = 0;
  int done = 0;
  while (done==0) {

    iter++;

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++)
      v[icv] = res[icv]/Ap[nbocv_i[icv]];

    double rho_prev = rho;

    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_rho += res[icv]*v[icv];

    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prev)<1.0E-20) rho_prev = 1.0E-20;
    double beta = rho/rho_prev;

    for (int icv = 0; icv<ncv; icv++)
      p[icv] = v[icv]+beta*p[icv];
    updateCvData(p, REPLACE_DATA);

    // v = [Ap]{p}...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      v[icv] = Ap[noc_f]*p[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        v[icv] += Ap[noc]*p[nbocv_v[noc]];
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_gamma += p[icv]*v[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;

    double alpha = rho/gamma;
    for (int icv = 0; icv<ncv_g; icv++)
      phi[icv] += alpha*p[icv];

    // check if we are done...
    if (iter%5==0) {

      double my_res_max = 0.0;

      // recompute the residual...
      updateCvData(phi, REPLACE_DATA);
      for (int icv = 0; icv<ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv+1]-1;
        res[icv] = Ap[noc_f]*phi[icv]; // diagonal
        for (int noc = noc_f+1; noc<=noc_l; noc++)
          res[icv] += Ap[noc]*phi[nbocv_v[noc]];
        // above is LHS. residual is then...
        res[icv] = rhs[icv]-res[icv];
        my_res_max = max(my_res_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank==0) {
        cout<<"cg iter, res_max: "<<iter<<" "<<res_max<<endl;

        if ((mode==ABSOLUTE_RESIDUAL)&&(res_max<=zero)) {
          done = 1;
        }
        else if ((mode==RELATIVE_RESIDUAL)&&(res_max/(res0_max+1.0E-20)<=zero)) {
          done = 1;
        }
        else if (iter>maxiter) {
          cout<<"Warning: solveCvScalarCg did not converge after "<<maxiter<<" iters, res_max: "<<res_max<<endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }
    else {

      // on the other iterations, use this approximation...
      for (int icv = 0; icv<ncv; icv++)
        res[icv] -= alpha*v[icv];

    }

  }

  //delete[]
  delete[] res;
  delete[] v;
  delete[] p;

  // let the calling routine know if we were successful...
  return (done==1);

}

int UgpWithCv::solveCvScalarBcgstab(double * phi, double * Ap, double * rhs, const double zero, const double zeroRel,
    const int maxiter, char *scalarName) {

  // we need the following work arrays...
  double *res = new double[ncv];
  double *res0 = new double[ncv];
  double *p = new double[ncv];
  double *v = new double[ncv];
  double *t = new double[ncv];
  double *s = new double[ncv];
  double *phat = new double[ncv_g];
  double *shat = new double[ncv_g];

  // initialize...
  for (int icv = 0; icv<ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  for (int icv = 0; icv<ncv_g; icv++) {
    phat[icv] = 0.0;
    shat[icv] = 0.0;
  }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv<ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    res[icv] = Ap[noc_f]*phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f+1; noc<=noc_l; noc++)
      res[icv] += Ap[noc]*phi[nbocv_v[noc]];

    res[icv] = rhs[icv]-res[icv];
    res0[icv] = res[icv];
  }

  double res0_max;
  double my_res0_max = 0.0;
  for (int icv = 0; icv<ncv; icv++)
    my_res0_max = max(my_res0_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

  int iter = 0;
  int done = 0;

  // check if we are already done
  if (mpi_rank==0) if (res0_max<=zero) {
    cout<<scalarName<<": "<<iter<<endl;
    done = 1;
  }
  MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

  while (done==0) {
    iter++;

    double rho_prime = -omega*rho;

    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_rho += res[icv]*res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
    double beta = alpha*rho/rho_prime;

    for (int icv = 0; icv<ncv; icv++)
      p[icv] = res[icv]-beta*(p[icv]-omega*v[icv]);

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++)
      phat[icv] = p[icv]/Ap[nbocv_i[icv]];
    updateCvData(phat, REPLACE_DATA);

    // v = [A]{phat}
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      v[icv] = Ap[noc_f]*phat[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        v[icv] += Ap[noc]*phat[nbocv_v[noc]];
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_gamma += v[icv]*res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;

    alpha = rho/gamma;

    for (int icv = 0; icv<ncv; icv++)
      s[icv] = res[icv]-alpha*v[icv];

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++)
      shat[icv] = s[icv]/Ap[nbocv_i[icv]];
    updateCvData(shat, REPLACE_DATA);

    // t = [A] shat...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      t[icv] = Ap[noc_f]*shat[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        t[icv] += Ap[noc]*shat[nbocv_v[noc]];
    }

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv<ncv; icv++) {
      my_buf[0] += s[icv]*t[icv];
      my_buf[1] += t[icv]*t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0]/(buf[1]+1.0E-20);

    // update phi...
    for (int icv = 0; icv<ncv; icv++)
      phi[icv] += alpha*phat[icv]+omega*shat[icv];
    updateCvData(phi, REPLACE_DATA);

    double my_res_max = 0.0;

    // recompute the residual...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      res[icv] = Ap[noc_f]*phi[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        res[icv] += Ap[noc]*phi[nbocv_v[noc]];
      // above is LHS. residual is then...
      res[icv] = rhs[icv]-res[icv];
      my_res_max = max(my_res_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
    }

    // check if we are done...
    if (iter%1==0) {
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

      if (mpi_rank==0) {

        lout(DEBUG_HI)<<scalarName<<":\t"<<"bcgstab iter, res_max: "<<iter<<" "<<res_max/(res0_max+1.0E-12)<<endl;

        if ((res_max<=zero)||(res_max/(res0_max+1.0E-12)<=zeroRel)||(iter>maxiter)) {
          lout(INFO_LO)<<scalarName<<"\t iter/resAbs/resRel/tresholds = "<<iter<<"\t"<<res_max<<"\t"<<res_max/(res0_max
              +1.0E-12)<<"\t"<<zero<<"\t"<<zeroRel<<endl;
          done = 1;
        }

        if (iter>maxiter) {
          cout<<"\nWarning: "<<scalarName<<" solveCvScalarBcgstab did not converge after "<<maxiter
              <<" iters, res_max: "<<res_max<<endl;
          done = 2;
        }
      }

      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }
  }

  delete[] res;
  delete[] res0;
  delete[] p;
  delete[] v;
  delete[] t;
  delete[] s;
  delete[] phat;
  delete[] shat;

  // let the calling routine know if we were successful...
  return (done==1);

}

//##############################################################################
//
//      math functions used for BCGSTAB,
//      ludeco ... LU decomposition
//      lusolv ... LU back substitution
//      luinv ... LU inverse
//
//##############################################################################
void ludeco(double(*A)[5], int order) {
  for (int jc = 1; jc<order; jc++)
    A[0][jc] /= A[0][0];

  int jrjc = 0;

  for (;;) {
    jrjc++;
    int jrjcm1 = jrjc-1;
    int jrjcp1 = jrjc+1;
    for (int jr = jrjc; jr<order; jr++) {
      double sum = A[jr][jrjc];
      for (int jm = 0; jm<=jrjcm1; jm++)
        sum -= A[jr][jm]*A[jm][jrjc];
      A[jr][jrjc] = sum;
    }
    if (jrjc==(order-1)) return;

    for (int jc = jrjcp1; jc<order; jc++) {
      double sum = A[jrjc][jc];

      for (int jm = 0; jm<=jrjcm1; jm++)
        sum -= A[jrjc][jm]*A[jm][jc];

      A[jrjc][jc] = sum/A[jrjc][jrjc];
    }
  }
}

void lusolv(double(*A)[5], double *c, int order) {
  //              ...First l(inv)*b...
  c[0] = c[0]/A[0][0];
  for (int jr = 1; jr<order; jr++) {
    int jrm1 = jr-1;
    double sum = c[jr];
    for (int jm = 0; jm<=jrm1; jm++)
      sum -= A[jr][jm]*c[jm];
    c[jr] = sum/A[jr][jr];
  }

  //             ...Next u(inv) of l(inv)*b...
  for (int jrjr = 1; jrjr<order; jrjr++) {
    int jr = (order-1)-jrjr;
    int jrp1 = jr+1;
    double sum = c[jr];
    for (int jmjm = jrp1; jmjm<order; jmjm++) {
      int jm = (order-1)-jmjm+jrp1;
      sum -= A[jr][jm]*c[jm];
    }
    c[jr] = sum;
  }
}

void luinv(double(*A)[5], double Ainv[5][5], int order) {
	double c[5];
	
	for (int i=0; i<order; i++) {
		if (i==0) {c[0] = 1.0; c[1] = 0.0; c[2] = 0.0; c[3] = 0.0; c[4] = 0.0;}
		if (i==1) {c[0] = 0.0; c[1] = 1.0; c[2] = 0.0; c[3] = 0.0; c[4] = 0.0;}
		if (i==2) {c[0] = 0.0; c[1] = 0.0; c[2] = 1.0; c[3] = 0.0; c[4] = 0.0;}
		if (i==3) {c[0] = 0.0; c[1] = 0.0; c[2] = 0.0; c[3] = 1.0; c[4] = 0.0;}
		if (i==4) {c[0] = 0.0; c[1] = 0.0; c[2] = 0.0; c[3] = 0.0; c[4] = 1.0;}
		
		//              ...First l(inv)*b...
		c[0] = c[0]/A[0][0];
		for (int jr = 1; jr<order; jr++) {
			int jrm1 = jr-1;
			double sum = c[jr];
			for (int jm = 0; jm<=jrm1; jm++)
				sum -= A[jr][jm]*c[jm];
			c[jr] = sum/A[jr][jr];
		}
		
		//             ...Next u(inv) of l(inv)*b...
		for (int jrjr = 1; jrjr<order; jrjr++) {
			int jr = (order-1)-jrjr;
			int jrp1 = jr+1;
			double sum = c[jr];
			for (int jmjm = jrp1; jmjm<order; jmjm++) {
				int jm = (order-1)-jmjm+jrp1;
				sum -= A[jr][jm]*c[jm];
			}
			c[jr] = sum;
		}
		
		Ainv[0][i] = c[0]; Ainv[1][i] = c[1]; Ainv[2][i] = c[2]; Ainv[3][i] = c[3]; Ainv[4][i] = c[4];
	}
	
}


void UgpWithCv::matTimesVecOverCVs(double(*res)[5], double(*Ap)[5][5], double(*phi)[5]) {
  for (int icv = 0; icv<ncv; icv++) {
    for (int i = 0; i<5; i++)
      res[icv][i] = 0.0;

    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    for (int noc = noc_f; noc<=noc_l; noc++) // cv diagonal + neighbors...
    {
      const int icv_nbr = nbocv_v[noc];

      for (int i = 0; i<5; i++)
        for (int j = 0; j<5; j++)
          res[icv][i] += Ap[noc][i][j]*phi[icv_nbr][j];
    }
  }
}

void UgpWithCv::UpdateCvDataStateVec(double(*phi)[5]) 
{
  static double *scal = new double[ncv_g];
  static double (*vec)[3] = new double[ncv_g][3];

  for (int icv = 0; icv<ncv_g; icv++)
    scal[icv] = phi[icv][0];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    phi[icv][0] = scal[icv];

  for (int icv = 0; icv<ncv_g; icv++)
    scal[icv] = phi[icv][4];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    phi[icv][4] = scal[icv];

  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<3; i++)
      vec[icv][i] = phi[icv][i+1];
  updateCvData(vec, REPLACE_ROTATE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<3; i++)
      phi[icv][i+1] = vec[icv][i];
}

void UgpWithCv::UpdateCvDataStateVecScal(double **phi, int nScal) 
{
  static double *scal = new double[ncv_g];
  static double (*vec)[3] = new double[ncv_g][3];

  for (int icv = 0; icv<ncv_g; icv++)
    scal[icv] = phi[icv][0];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    phi[icv][0] = scal[icv];

  for (int icv = 0; icv<ncv_g; icv++)
    scal[icv] = phi[icv][4];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    phi[icv][4] = scal[icv];

  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<3; i++)
      vec[icv][i] = phi[icv][i+1];
  updateCvData(vec, REPLACE_ROTATE_DATA);
  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<3; i++)
      phi[icv][i+1] = vec[icv][i];
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    for (int icv = 0; icv<ncv_g; icv++)
      scal[icv] = phi[icv][5+iScal];
    updateCvData(scal, REPLACE_DATA);
    for (int icv = 0; icv<ncv_g; icv++)
      phi[icv][5+iScal] = scal[icv];
  }
  
}

int UgpWithCv::solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
    const double zeroRel, const int maxiter, char *scalarName) {

  // we need the following work arrays...
  static double (*res)[5] = new double[ncv][5];
  static double (*res0)[5] = new double[ncv][5];
  static double (*p)[5] = new double[ncv][5];
  static double (*v)[5] = new double[ncv][5];
  static double (*t)[5] = new double[ncv][5];
  static double (*s)[5] = new double[ncv][5];
  static double (*phat)[5] = new double[ncv_g][5];
  static double (*shat)[5] = new double[ncv_g][5];

  static double (*LUDEC)[5][5] = new double[ncv][5][5]; // LU decomposed diagonal for speed up
  double b[5];

  // LU decompose diagonal and store
  for (int icv = 0; icv<ncv; icv++) {
    for (int i = 0; i<5; i++)
      for (int j = 0; j<5; j++)
        LUDEC[icv][i][j] = Ap[nbocv_i[icv]][i][j];

    ludeco(LUDEC[icv], 5);
  }

  // initialize...
  for (int icv = 0; icv<ncv; icv++)
    for (int i = 0; i<5; i++) {
      p[icv][i] = 0.0;
      v[icv][i] = 0.0;
    }

  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<5; i++) {
      phat[icv][i] = 0.0;
      shat[icv][i] = 0.0;
    }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format ...
  matTimesVecOverCVs(res, Ap, phi);

  for (int icv = 0; icv<ncv; icv++)
    for (int i = 0; i<5; i++) {
      res[icv][i] = (rhs[icv][i]-res[icv][i]);
      res0[icv][i] = res[icv][i];
    }

  // compute first residual  
  double res0_max, my_res0_max = 0.0;

  for (int icv = 0; icv<ncv; icv++) {
    double tmp[5];

    for (int i = 0; i<5; i++)
      tmp[i] = res[icv][i];

    lusolv(LUDEC[icv], tmp, 5);

    for (int i = 0; i<5; i++)
      my_res0_max = max(my_res0_max, fabs(tmp[i]));
  }

  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

  // start iteration

  int iter = 0;
  int done = 0;
  double res_max;

  while (done==0) {

    iter++;

    double rho_prime = -omega*rho;

    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        my_rho += res[icv][i]*res0[icv][i];

    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;

    double beta = alpha*rho/rho_prime;

    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        p[icv][i] = res[icv][i]-beta*(p[icv][i]-omega*v[icv][i]);

    // block-diagonal precon...
    // solve the system
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        phat[icv][i] = p[icv][i];
      lusolv(LUDEC[icv], phat[icv], 5);
    }

    UpdateCvDataStateVec(phat);

    // v = [A]{phat}
    matTimesVecOverCVs(v, Ap, phat);

    double my_gamma = 0.0;

    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        my_gamma += v[icv][i]*res0[icv][i];

    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    //    gamma = my_gamma;

    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;

    alpha = rho/gamma;

    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        s[icv][i] = res[icv][i]-alpha*v[icv][i];

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        shat[icv][i] = s[icv][i];
      lusolv(LUDEC[icv], shat[icv], 5);
    }

    UpdateCvDataStateVec(shat);

    // t = [A] shat...
    matTimesVecOverCVs(t, Ap, shat);

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++) {
        my_buf[0] += s[icv][i]*t[icv][i];
        my_buf[1] += t[icv][i]*t[icv][i];
      }

    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0]/(buf[1]+1.0E-20);

    // update phi...
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        phi[icv][i] += alpha*phat[icv][i]+omega*shat[icv][i];

    UpdateCvDataStateVec(phi);

    double my_res_max = 0.0;

    // recompute the residual...
    matTimesVecOverCVs(res, Ap, phi);

    for (int icv = 0; icv<ncv; icv++) // above is LHS. residual is then...
      for (int i = 0; i<5; i++)
        res[icv][i] = rhs[icv][i]-res[icv][i];

    // compute normalized res
    for (int icv = 0; icv<ncv; icv++) {
      double tmp[5];

      for (int i = 0; i<5; i++)
        tmp[i] = res[icv][i];

      lusolv(LUDEC[icv], tmp, 5);

      for (int i = 0; i<5; i++)
        my_res_max = max(my_res_max, fabs(tmp[i]));
    }

    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

    // check residual 
    if (mpi_rank==0) {
      //      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;

      if ((res_max<=zeroAbs)||(res_max/(res0_max+1.0E-12)<=zeroRel)) done = 1;
      else if (iter>maxiter) {
        /*        cout << "Warning: " << scalarName
         << " solveCvScalarBcgstab did not converge after " << maxiter
         << " iters, res_max: " << res_max << endl;*/

        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
  }

  //lout(INFO_LO)<<"\tbcgstab: iter: "<<iter<<" res_max: "<<res_max<<endl;
  lout(INFO_LO)<<"navier"<<"\t iter/resAbs/resRel/tresholds = "<<iter<<"\t"<<res_max<<"\t"<<res_max/(res0_max
              +1.0E-12)<<"\t"<<zeroAbs<<"\t"<<zeroRel<<endl;

  /*    delete[] res;
   delete[] res0;
   delete[] p;
   delete[] v;
   delete[] t;
   delete[] s;
   delete[] phat;
   delete[] shat;

   delete [] LU;*/

  // let the calling routine know if we were successful...
  return (done==1);
}

int UgpWithCv::solveCvScalarBcgstab(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double zero,
    const int maxiter) {

  // need gradients in ghost as well...
  double (*grad)[3] = new double[ncv_g][3];

  // we need the following work arrays...
  double * res = new double[ncv];
  double * res0 = new double[ncv];
  double * p = new double[ncv];
  double * v = new double[ncv];
  double * t = new double[ncv];
  double * s = new double[ncv];
  double * phat = new double[ncv_g];
  double * shat = new double[ncv_g];

  // initialize...
  for (int icv = 0; icv<ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv<ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    for (int i = 0; i<3; i++)
      grad[icv][i] = nbocv_lsg_coeff[noc_f][i]*phi[icv];
    for (int noc = noc_f+1; noc<=noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      for (int i = 0; i<3; i++)
        grad[icv][i] += nbocv_lsg_coeff[noc][i]*phi[icv_nbr];
    }
  }
  updateCvData(grad, REPLACE_ROTATE_DATA);
  for (int icv = 0; icv<ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    res[icv] = Ap[noc_f]*phi[icv]; // diagonal
    for (int i = 0; i<3; i++)
      res[icv] += Ap_grad[noc_f][i]*grad[icv][i]; // diagonal gradient
    // cv neighbors...
    for (int noc = noc_f+1; noc<=noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      res[icv] += Ap[noc]*phi[icv_nbr];
      for (int i = 0; i<3; i++)
        res[icv] += Ap_grad[noc][i]*grad[icv_nbr][i];
    }
    res[icv] = rhs[icv]-res[icv];
    res0[icv] = res[icv];
  }

  int iter = 0;
  int done = 0;
  while (done==0) {

    iter++;

    double rho_prime = -omega*rho;

    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_rho += res[icv]*res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
    double beta = alpha*rho/rho_prime;

    for (int icv = 0; icv<ncv; icv++)
      p[icv] = res[icv]-beta*(p[icv]-omega*v[icv]);

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++)
      phat[icv] = p[icv]/Ap[nbocv_i[icv]];
    updateCvData(phat, REPLACE_DATA);

    // v = [A]{phat}
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      for (int i = 0; i<3; i++)
        grad[icv][i] = nbocv_lsg_coeff[noc_f][i]*phat[icv];
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i<3; i++)
          grad[icv][i] += nbocv_lsg_coeff[noc][i]*phat[icv_nbr];
      }
    }
    updateCvData(grad, REPLACE_ROTATE_DATA);
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      v[icv] = Ap[noc_f]*phat[icv]; // diagonal
      for (int i = 0; i<3; i++)
        v[icv] += Ap_grad[noc_f][i]*grad[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        v[icv] += Ap[noc]*phat[icv_nbr];
        for (int i = 0; i<3; i++)
          v[icv] += Ap_grad[noc][i]*grad[icv_nbr][i];
      }
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_gamma += v[icv]*res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;

    alpha = rho/gamma;

    for (int icv = 0; icv<ncv; icv++)
      s[icv] = res[icv]-alpha*v[icv];

    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++)
      shat[icv] = s[icv]/Ap[nbocv_i[icv]];
    updateCvData(shat, REPLACE_DATA);

    // t = [A] shat...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      for (int i = 0; i<3; i++)
        grad[icv][i] = nbocv_lsg_coeff[noc_f][i]*shat[icv];
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i<3; i++)
          grad[icv][i] += nbocv_lsg_coeff[noc][i]*shat[icv_nbr];
      }
    }
    updateCvData(grad, REPLACE_ROTATE_DATA);
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      t[icv] = Ap[noc_f]*shat[icv]; // diagonal
      for (int i = 0; i<3; i++)
        t[icv] += Ap_grad[noc_f][i]*grad[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f+1; noc<=noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        t[icv] += Ap[noc]*shat[icv_nbr];
        for (int i = 0; i<3; i++)
          t[icv] += Ap_grad[noc][i]*grad[icv_nbr][i];
      }
    }

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv<ncv; icv++) {
      my_buf[0] += s[icv]*t[icv];
      my_buf[1] += t[icv]*t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0]/(buf[1]+1.0E-20);

    // update phi...
    for (int icv = 0; icv<ncv; icv++)
      phi[icv] += alpha*phat[icv]+omega*shat[icv];
    updateCvData(phi, REPLACE_DATA);

    // check if we are done...
    if (iter%5==0) {

      double my_res_max = 0.0;

      // recompute the residual...
      for (int icv = 0; icv<ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv+1]-1;
        for (int i = 0; i<3; i++)
          grad[icv][i] = nbocv_lsg_coeff[noc_f][i]*phi[icv];
        for (int noc = noc_f+1; noc<=noc_l; noc++) {
          int icv_nbr = nbocv_v[noc];
          for (int i = 0; i<3; i++)
            grad[icv][i] += nbocv_lsg_coeff[noc][i]*phi[icv_nbr];
        }
      }
      updateCvData(grad, REPLACE_ROTATE_DATA);
      for (int icv = 0; icv<ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv+1]-1;
        res[icv] = Ap[noc_f]*phi[icv]; // diagonal
        for (int i = 0; i<3; i++)
          res[icv] += Ap_grad[noc_f][i]*grad[icv][i]; // diagonal gradient
        // cv neighbors...
        for (int noc = noc_f+1; noc<=noc_l; noc++) {
          int icv_nbr = nbocv_v[noc];
          res[icv] += Ap[noc]*phi[icv_nbr];
          for (int i = 0; i<3; i++)
            res[icv] += Ap_grad[noc][i]*grad[icv_nbr][i];
        }
        res[icv] = rhs[icv]-res[icv];
        my_res_max = max(my_res_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank==0) {
        //cout << "bcgstab iter, res_max: " << iter << " " << res_max << endl;
        if (res_max<zero) {
          done = 1;
        }
        else if (iter>maxiter) {
          cout<<"Warning: solveCvScalarBcgstab did not converge after "<<maxiter<<" iters, res_max: "<<res_max<<endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }
    else {

      // on the other iterations, use this approximation...
      for (int icv = 0; icv<ncv; icv++)
        res[icv] = s[icv]-omega*t[icv];

    }

  }

  //delete[]
  delete[] grad;
  delete[] res;
  delete[] res0;
  delete[] p;
  delete[] v;
  delete[] t;
  delete[] s;
  delete[] phat;
  delete[] shat;

  // let the calling routine know if we were successful...
  return (done==1);

}

int UgpWithCv::solveCvScalarLusgs(double * phi, double * Ap, double * rhs, const double zero, const double zeroRel,
																	const int maxiter, char *scalarName) {
	
	
	static double (*LowProd) = new double[ncv];
	static double (*UpperProd) = new double[ncv];
	static double (*DiagProd) = new double[ncv];
	static double (*AuxVector) = new double[ncv];
	static double (*x_n) = new double[ncv];
	
	// First part of the symmetric iteration: (D+L).x^* = b
	// compute Lx^*
	LowerProductScalar(LowProd, Ap, phi);  
	
	// compute aux_vector = b-L.x^*, and solve Dx^* = aux_vector, using a LU method
	for (int icv = 0; icv<ncv; icv++) {
		x_n[icv] = (rhs[icv]-LowProd[icv])/Ap[nbocv_i[icv]];
	}
	
	// Second part of the symmetric iteration: (D+U).x^1 = D.x^*
	// compute D.x^*
	DiagonalProductScalar(DiagProd, Ap, x_n); 
	
	// compute aux_vector = D.x^*
	for (int icv = 0; icv<ncv; icv++)
		AuxVector[icv] = DiagProd[icv];
	
	// compute U.x^1
	UpperProductScalar(UpperProd, Ap, x_n);
	
	// compute aux_vector = D.x^*-U.x^1
	for (int icv = 0; icv<ncv; icv++)
		AuxVector[icv] -= UpperProd[icv];
	
	// solve D.x^1 = aux_vector, using a LU method
	double res_max;
	for (int icv = 0; icv<ncv; icv++) {
		res_max = max(res_max, fabs(AuxVector[icv]));
		x_n[icv] = AuxVector[icv]/Ap[nbocv_i[icv]];
	}
	
	// update phi...
	for (int icv = 0; icv<ncv; icv++)
		phi[icv] = x_n[icv];
	
	// let the calling routine know if we were successful...
  return 1;
}

int UgpWithCv::solveCvVectorR5Lusgs(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
																		const double zeroRel, const int maxiter, char *scalarName) {
	
	
	static double (*LowProd)[5] = new double[ncv][5];
	static double (*UpperProd)[5] = new double[ncv][5];
	static double (*DiagProd)[5] = new double[ncv][5];
	static double (*AuxVector)[5] = new double[ncv][5];
	static double (*x_n)[5] = new double[ncv][5];
	double tmp[5];
	static double (*LUDEC)[5][5] = new double[ncv][5][5];
	
	// LU decompose diagonal and store
	for (int icv = 0; icv<ncv; icv++) {
		for (int i = 0; i<5; i++)
			for (int j = 0; j<5; j++)
				LUDEC[icv][i][j] = Ap[nbocv_i[icv]][i][j];
		ludeco(LUDEC[icv], 5);
	}
	
	// First part of the symmetric iteration: (D+L).x^* = b
	// compute Lx^*
	LowerProductVectorR5(LowProd, Ap, phi);  
	
	// compute aux_vector = b-L.x^*, and solve Dx^* = aux_vector, using a LU method
	for (int icv = 0; icv<ncv; icv++) {
		for (int i = 0; i<5; i++) 
			tmp[i] = rhs[icv][i]-LowProd[icv][i];
		lusolv(LUDEC[icv], tmp, 5);
		for (int i = 0; i<5; i++) 
			x_n[icv][i] = tmp[i];
	}
	
	// Second part of the symmetric iteration: (D+U).x^1 = D.x^*
	// compute D.x^*
	DiagonalProductVectorR5(DiagProd, Ap, x_n); 
	
	// compute aux_vector = D.x^*
	for (int icv = 0; icv<ncv; icv++)
		for (int i = 0; i<5; i++)
			AuxVector[icv][i] = DiagProd[icv][i];
	
	// compute U.x^1
	UpperProductVectorR5(UpperProd, Ap, x_n);
	
	// compute aux_vector = D.x^*-U.x^1
	for (int icv = 0; icv<ncv; icv++)
		for (int i = 0; i<5; i++)
			AuxVector[icv][i] -= UpperProd[icv][i];
	
	// solve D.x^1 = aux_vector, using a LU method
	double res_max;
	for (int icv = 0; icv<ncv; icv++) {
		lusolv(LUDEC[icv], AuxVector[icv], 5);
		for (int i = 0; i<5; i++) {
			res_max = max(res_max, fabs(AuxVector[icv][i]));
			x_n[icv][i] = AuxVector[icv][i];
		}
	}
	
	// update phi...
	for (int icv = 0; icv<ncv; icv++)
		for (int i = 0; i<5; i++)
			phi[icv][i] = x_n[icv][i];
	
	// let the calling routine know if we were successful...
  return 1;
}

void UgpWithCv::LowerProductScalar(double *LowProd, double *Ap, double *Phi) {
	for (int icv = 0; icv<ncv; icv++) {
		LowProd[icv] = 0.0;
		int noc_f = nbocv_i[icv]+1;
		int noc_l = nbocv_i[icv+1]-1;
		for (int noc = noc_f; noc<=noc_l; noc++) {
			int icv_nbr = nbocv_v[noc];
			if (icv_nbr < icv) {
				LowProd[icv] += Ap[noc]*Phi[icv_nbr];
			}
		}
	}
}


void UgpWithCv::UpperProductScalar(double *UpperProd, double *Ap, double *Phi) {
	for (int icv = 0; icv<ncv; icv++) {
		UpperProd[icv] = 0.0;
		int noc_f = nbocv_i[icv]+1;
		int noc_l = nbocv_i[icv+1]-1;
		for (int noc = noc_f; noc<=noc_l; noc++) {
			int icv_nbr = nbocv_v[noc];
			if (icv_nbr > icv) {
				UpperProd[icv] += Ap[noc]*Phi[icv_nbr];
			}
		}
	}
}


void UgpWithCv::DiagonalProductScalar(double *DiagProd, double *Ap, double *Phi) {
	for (int icv = 0; icv<ncv; icv++) {
		DiagProd[icv] = 0.0;
		int icv_nbr = nbocv_v[nbocv_i[icv]];
		DiagProd[icv] += Ap[nbocv_i[icv]]*Phi[icv_nbr];
	}
}


void UgpWithCv::LowerProductVectorR5(double(*LowProd)[5], double(*Ap)[5][5], double(*Phi)[5]) {
	for (int icv = 0; icv<ncv; icv++) {
		for (int i = 0; i<5; i++)
			LowProd[icv][i] = 0.0;
		
		int noc_f = nbocv_i[icv]+1;
		int noc_l = nbocv_i[icv+1]-1;
		for (int noc = noc_f; noc<=noc_l; noc++) {
			int icv_nbr = nbocv_v[noc];
			if (icv_nbr < icv) {
				for (int i = 0; i<5; i++)
					for (int j = 0; j<5; j++)
						LowProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
			}
		}
	}
}



void UgpWithCv::UpperProductVectorR5(double(*UpperProd)[5], double(*Ap)[5][5], double(*Phi)[5]) {
	for (int icv = 0; icv<ncv; icv++) {
		for (int i = 0; i<5; i++)
			UpperProd[icv][i] = 0.0;
		
		int noc_f = nbocv_i[icv]+1;
		int noc_l = nbocv_i[icv+1]-1;
		for (int noc = noc_f; noc<=noc_l; noc++) {
			int icv_nbr = nbocv_v[noc];
			if (icv_nbr > icv) {
				for (int i = 0; i<5; i++)
					for (int j = 0; j<5; j++)
						UpperProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
			}
		}
	}
}


void UgpWithCv::DiagonalProductVectorR5(double(*DiagProd)[5], double(*Ap)[5][5], double(*Phi)[5]) {
	for (int icv = 0; icv<ncv; icv++) {
		for (int i = 0; i<5; i++)
			DiagProd[icv][i] = 0.0;
		
		int icv_nbr = nbocv_v[nbocv_i[icv]];
		for (int i = 0; i<5; i++)
			for (int j = 0; j<5; j++)
				DiagProd[icv][i] += Ap[nbocv_i[icv]][i][j]*Phi[icv_nbr][j];
	}
}

int UgpWithCv::solveCvScalarBcgstabLine(double * phi, double * Ap, double * rhs, const double zero, const double zeroRel,
																				const int maxiter, char *scalarName) {
	
  // we need the following work arrays...
  double *res = new double[ncv];
  double *res0 = new double[ncv];
  double *p = new double[ncv];
  double *v = new double[ncv];
  double *t = new double[ncv];
  double *s = new double[ncv];
  double *phat = new double[ncv_g];
  double *shat = new double[ncv_g];
	
	double b, m, aux_0, aux_1, aux_2, invA, aux_A;
	int noc00, noc01, noc11, noc10;
	
	int n=0;
	for (int iSources = 0; iSources < nSources; iSources++) 
		n = max (n,int(linelet_cv[iSources].size()));
	
	double (*mhat) = new double[n];
	double (*bhat) = new double[n];
	double (*LUbhat) = new double[n];
	double (*invbhat) = new double[n];
	double (*vhat) = new double[n];
	
  // initialize...
  for (int icv = 0; icv<ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  for (int icv = 0; icv<ncv_g; icv++) {
    phat[icv] = 0.0;
    shat[icv] = 0.0;
  }
	
  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
	
  // calculate the residual in rhs format...
  for (int icv = 0; icv<ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv+1]-1;
    res[icv] = Ap[noc_f]*phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f+1; noc<=noc_l; noc++)
      res[icv] += Ap[noc]*phi[nbocv_v[noc]];
		
    res[icv] = rhs[icv]-res[icv];
    res0[icv] = res[icv];
  }
	
  double res0_max;
  double my_res0_max = 0.0;
  for (int icv = 0; icv<ncv; icv++)
    my_res0_max = max(my_res0_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
	
  int iter = 0;
  int done = 0;
	
  // check if we are already done
  if (mpi_rank==0) if (res0_max<=zero) {
    cout<<scalarName<<": "<<iter<<endl;
    done = 1;
  }
  MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
	
  while (done==0) {
    iter++;
		
    double rho_prime = -omega*rho;
		
    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_rho += res[icv]*res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
    double beta = alpha*rho/rho_prime;
		
    for (int icv = 0; icv<ncv; icv++)
      p[icv] = res[icv]-beta*(p[icv]-omega*v[icv]);
		
		/*********** JACOBI ************/		
/*    for (int icv = 0; icv<ncv; icv++)
      phat[icv] = p[icv]/Ap[nbocv_i[icv]];
    updateCvData(phat, REPLACE_DATA); */
		/*********** JACOBI ************/		
    
		/*********** LINE IMPLICIT ************/
		for (int icv = 0; icv < ncv; icv++) 
			if (!Linelet[icv]) {
				for (int i = 0; i < 5; i++)
					phat[icv] = p[icv]/Ap[nbocv_i[icv]];
			}
		
		// Solve linelet
		for (int iSources = 0; iSources < nSources; iSources++) {
			int n = linelet_cv[iSources].size();
			
			// LU descompusition
			vhat[0]=p[linelet_cv[iSources][0]];
			bhat[0]=Ap[nbocv_i[linelet_cv[iSources][0]]];
			
			for (int i = 1; i < n; i++) {
				int cv_im1 = linelet_cv[iSources][i-1];
				int cv_i = linelet_cv[iSources][i];
				
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_im1, cv_i);
				invbhat[i-1]= 1.0/bhat[i-1];
				
				mhat[i] = Ap[noc10]*invbhat[i-1];
				aux_0 = mhat[i] * Ap[noc01];
				bhat[i] = Ap[noc11]-aux_0;
				
				// forward substituton
				aux_1 = mhat[i] * vhat[i-1];											
				vhat[i] = p[linelet_cv[iSources][i]] - aux_1;
			}
			
			// backward substituton
			int cv_nm1 = linelet_cv[iSources][n-1];
			invbhat[n-1] = 1.0/bhat[n-1];
			phat[cv_nm1] = invbhat[n-1]* vhat[n-1];
			
			for (int i = n-2; i >= 0; i--) {
				int cv_i = linelet_cv[iSources][i];
				int cv_ip1 = linelet_cv[iSources][i+1];
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_i, cv_ip1);
				aux_2 = Ap[noc01]*phat[cv_ip1];
				aux_1 = vhat[i]-aux_2;
				phat[cv_i] = invbhat[i]*aux_1;
			}
		}
		
		updateCvData(phat, REPLACE_DATA);
		/*********** LINE IMPLICIT ************/		
		
		// v = [A]{phat}
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      v[icv] = Ap[noc_f]*phat[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        v[icv] += Ap[noc]*phat[nbocv_v[noc]];
    }
		
    double my_gamma = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      my_gamma += v[icv]*res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;
		
    alpha = rho/gamma;
		
    for (int icv = 0; icv<ncv; icv++)
      s[icv] = res[icv]-alpha*v[icv];
		
		/*********** JACOBI ************/		
/*    for (int icv = 0; icv<ncv; icv++)
      shat[icv] = s[icv]/Ap[nbocv_i[icv]];
    updateCvData(shat, REPLACE_DATA);*/
		/*********** JACOBI ************/		
		
		/*********** LINE IMPLICIT ************/
		for (int icv = 0; icv < ncv; icv++) 
			if (!Linelet[icv]) {
				for (int i = 0; i < 5; i++)
					shat[icv] = s[icv]/Ap[nbocv_i[icv]];
			}
		
		// Solve linelet
		for (int iSources = 0; iSources < nSources; iSources++) {
			int n = linelet_cv[iSources].size();
			
			// LU descompusition
			vhat[0]=s[linelet_cv[iSources][0]];
			bhat[0]=Ap[nbocv_i[linelet_cv[iSources][0]]];
			
			for (int i = 1; i < n; i++) {
				int cv_im1 = linelet_cv[iSources][i-1];
				int cv_i = linelet_cv[iSources][i];
				
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_im1, cv_i);
				invbhat[i-1]= 1.0/bhat[i-1];
				
				mhat[i] = Ap[noc10]*invbhat[i-1];
				aux_0 = mhat[i] * Ap[noc01];
				bhat[i] = Ap[noc11]-aux_0;
				
				// forward substituton
				aux_1 = mhat[i] * vhat[i-1];											
				vhat[i] = s[linelet_cv[iSources][i]] - aux_1;
			}
			
			// backward substituton
			int cv_nm1 = linelet_cv[iSources][n-1];
			invbhat[n-1] = 1.0/bhat[n-1];
			shat[cv_nm1] = invbhat[n-1]* vhat[n-1];
			
			for (int i = n-2; i >= 0; i--) {
				int cv_i = linelet_cv[iSources][i];
				int cv_ip1 = linelet_cv[iSources][i+1];
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_i, cv_ip1);
				aux_2 = Ap[noc01]*shat[cv_ip1];
				aux_1 = vhat[i]-aux_2;
				shat[cv_i] = invbhat[i]*aux_1;
			}
		}
		
    updateCvData(shat, REPLACE_DATA);
		/*********** LINE IMPLICIT ************/		
				
    // t = [A] shat...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      t[icv] = Ap[noc_f]*shat[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        t[icv] += Ap[noc]*shat[nbocv_v[noc]];
    }
		
    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv<ncv; icv++) {
      my_buf[0] += s[icv]*t[icv];
      my_buf[1] += t[icv]*t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    omega = buf[0]/(buf[1]+1.0E-20);
		
    // update phi...
    for (int icv = 0; icv<ncv; icv++)
      phi[icv] += alpha*phat[icv]+omega*shat[icv];
    updateCvData(phi, REPLACE_DATA);
		
    double my_res_max = 0.0;
		
    // recompute the residual...
    for (int icv = 0; icv<ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1]-1;
      res[icv] = Ap[noc_f]*phi[icv]; // diagonal
      for (int noc = noc_f+1; noc<=noc_l; noc++)
        res[icv] += Ap[noc]*phi[nbocv_v[noc]];
      // above is LHS. residual is then...
      res[icv] = rhs[icv]-res[icv];
      my_res_max = max(my_res_max, fabs(res[icv]/Ap[nbocv_i[icv]]));
    }
		
    // check if we are done...
    if (iter%1==0) {
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
			
      if (mpi_rank==0) {
				
        lout(DEBUG_HI)<<scalarName<<":\t"<<"bcgstab iter, res_max: "<<iter<<" "<<res_max/(res0_max+1.0E-12)<<endl;
//       cout<<scalarName<<":\t"<<"bcgstab iter, res_max: "<<iter<<" "<<res_max/(res0_max+1.0E-12)<<endl;
				
        if ((res_max<=zero)||(res_max/(res0_max+1.0E-12)<=zeroRel)) {
          lout(INFO_LO)<<scalarName<<"\t iter/resAbs/resRel/tresholds = "<<iter<<"\t"<<res_max<<"\t"<<res_max/(res0_max
																																																							 +1.0E-12)<<"\t"<<zero<<"\t"<<zeroRel<<endl;
          done = 1;
        }
				
        if (iter>maxiter) {
          cout<<"\nWarning: "<<scalarName<<" solveCvScalarBcgstab did not converge after "<<maxiter
					<<" iters, res_max: "<<res_max<<endl;
          done = 2;
        }
      }
			
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }
  }
	
  delete[] res;
  delete[] res0;
  delete[] p;
  delete[] v;
  delete[] t;
  delete[] s;
  delete[] phat;
  delete[] shat;
	
  // let the calling routine know if we were successful...
  return iter;
}

int UgpWithCv::solveCvVectorR5BcgstabLine(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
																					const double zeroRel, const int maxiter, char *scalarName) {
	
  // we need the following work arrays...
  static double (*res)[5] = new double[ncv][5];
  static double (*res0)[5] = new double[ncv][5];
  static double (*p)[5] = new double[ncv][5];
  static double (*v)[5] = new double[ncv][5];
  static double (*t)[5] = new double[ncv][5];
  static double (*s)[5] = new double[ncv][5];
  static double (*phat)[5] = new double[ncv_g][5];
  static double (*shat)[5] = new double[ncv_g][5];
	static double (*LUDEC)[5][5] = new double[ncv][5][5]; // LU decomposed diagonal for speed up
	
	double b[5], m[5][5], aux_0[5][5], aux_1[5], aux_2[5], invA[5][5], aux_A[5][5];
	int noc00, noc01, noc11, noc10;
		
	int n=0;
	for (int iSources = 0; iSources < nSources; iSources++) 
		n = max (n,int(linelet_cv[iSources].size()));
	
	double (*mhat)[5][5] = new double[n][5][5];
	double (*bhat)[5][5] = new double[n][5][5];
	double (*LUbhat)[5][5] = new double[n][5][5];
	double (*invbhat)[5][5] = new double[n][5][5];
	double (*vhat)[5] = new double[n][5];
	
  // LU decompose diagonal and store
  for (int icv = 0; icv<ncv; icv++) {
    for (int i = 0; i<5; i++)
      for (int j = 0; j<5; j++)
        LUDEC[icv][i][j] = Ap[nbocv_i[icv]][i][j];
		
    ludeco(LUDEC[icv], 5);
  }
	
  // initialize...
  for (int icv = 0; icv<ncv; icv++)
    for (int i = 0; i<5; i++) {
      p[icv][i] = 0.0;
      v[icv][i] = 0.0;
    }
	
  for (int icv = 0; icv<ncv_g; icv++)
    for (int i = 0; i<5; i++) {
      phat[icv][i] = 0.0;
      shat[icv][i] = 0.0;
    }
	
  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
	
  // calculate the residual in rhs format ...
  matTimesVecOverCVs(res, Ap, phi);
	
  for (int icv = 0; icv<ncv; icv++)
    for (int i = 0; i<5; i++) {
      res[icv][i] = (rhs[icv][i]-res[icv][i]);
      res0[icv][i] = res[icv][i];
    }
	
  // compute first residual  
  double res0_max, my_res0_max = 0.0;
	
  for (int icv = 0; icv<ncv; icv++) {
    double tmp[5];
		
    for (int i = 0; i<5; i++)
      tmp[i] = res[icv][i];
		
    lusolv(LUDEC[icv], tmp, 5);
		
    for (int i = 0; i<5; i++)
      my_res0_max = max(my_res0_max, fabs(tmp[i]));
  }
	
  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
	
  // start iteration
	
  int iter = 0;
  int done = 0;
  double res_max;
	
  while (done==0) {
		
    iter++;
		
    double rho_prime = -omega*rho;
		
    double my_rho = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        my_rho += res[icv][i]*res0[icv][i];
		
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
		
    double beta = alpha*rho/rho_prime;
		
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        p[icv][i] = res[icv][i]-beta*(p[icv][i]-omega*v[icv][i]);
		
		
		
		/*********** NO PRECONDITIONING ************/		
/*    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        phat[icv][i] = p[icv][i];
    }
    UpdateCvDataStateVec(phat);*/
		/*********** NO PRECONDITIONING ************/		
		
		/*********** JACOBI ************/		
 /*   // block-diagonal precon...
    // solve the system
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        phat[icv][i] = p[icv][i];
      lusolv(LUDEC[icv], phat[icv], 5);
    }
    UpdateCvDataStateVec(phat);*/
		/*********** JACOBI ************/		
		
		/*********** LINE IMPLICIT ************/
		for (int icv = 0; icv < ncv; icv++) 
			if (!Linelet[icv]) {
				for (int i = 0; i < 5; i++)
					phat[icv][i] = p[icv][i];
				lusolv(LUDEC[icv], phat[icv], 5);
			}
		
		// Solve linelet
		for (int iSources = 0; iSources < nSources; iSources++) {
			int n = linelet_cv[iSources].size();
			
			// LU descompusition
			for (int i = 0; i < 5; i++) {
				vhat[0][i]=p[linelet_cv[iSources][0]][i];
				for (int j = 0; j < 5; j++)
					bhat[0][i][j]=Ap[nbocv_i[linelet_cv[iSources][0]]][i][j];
			}
			
			for (int i = 1; i < n; i++) {
				int cv_im1 = linelet_cv[iSources][i-1];
				int cv_i = linelet_cv[iSources][i];
				
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_im1, cv_i);
				for (int iVar = 0; iVar < 5; iVar++)
					for (int jVar = 0; jVar < 5; jVar++)
						LUbhat[i-1][iVar][jVar] = bhat[i-1][iVar][jVar];
				ludeco(LUbhat[i-1], 5); luinv(LUbhat[i-1], invbhat[i-1], 5);
				
				getMultMatrix(mhat[i], Ap[noc10], invbhat[i-1]);
				getMultMatrix(aux_0, mhat[i], Ap[noc01]);
				getSubsMatrix(bhat[i], Ap[noc11], aux_0);
				
				// forward substituton
				getMultMatrix(aux_1, mhat[i], vhat[i-1]);											
				getSubsMatrix(vhat[i], p[linelet_cv[iSources][i]], aux_1);
			}
			
			// backward substituton
			int cv_nm1 = linelet_cv[iSources][n-1];
			for (int iVar = 0; iVar < 5; iVar++)
				for (int jVar = 0; jVar < 5; jVar++)
					LUbhat[n-1][iVar][jVar] = bhat[n-1][iVar][jVar];
			ludeco(LUbhat[n-1], 5); luinv(LUbhat[n-1], invbhat[n-1], 5);
			getMultMatrix(phat[cv_nm1], invbhat[n-1], vhat[n-1]);
			
			for (int i = n-2; i >= 0; i--) {
				int cv_i = linelet_cv[iSources][i];
				int cv_ip1 = linelet_cv[iSources][i+1];
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_i, cv_ip1);
				getMultMatrix(aux_2, Ap[noc01], phat[cv_ip1]);
				getSubsMatrix(aux_1, vhat[i], aux_2);
				getMultMatrix(phat[cv_i], invbhat[i], aux_1);
			}
		}
		
		UpdateCvDataStateVec(phat);
		/*********** LINE IMPLICIT ************/		
		
    // v = [A]{phat}
    matTimesVecOverCVs(v, Ap, phat);
		
    double my_gamma = 0.0;
		
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        my_gamma += v[icv][i]*res0[icv][i];
		
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    //    gamma = my_gamma;
		
    if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;
		
    alpha = rho/gamma;
		
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        s[icv][i] = res[icv][i]-alpha*v[icv][i];
		
		/*********** NO PRECONDITIONING ************/				
/*    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        shat[icv][i] = s[icv][i];
    }
    UpdateCvDataStateVec(shat);*/
		/*********** NO PRECONDITIONING ************/	
		
		/*********** JACOBI ************/				
/*    // diagonal precon...
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<5; i++)
        shat[icv][i] = s[icv][i];
      lusolv(LUDEC[icv], shat[icv], 5);
    }
    UpdateCvDataStateVec(shat);*/
		/*********** JACOBI ************/		
		
		/*********** LINE IMPLICIT ************/
		for (int icv = 0; icv < ncv; icv++) 
			if (!Linelet[icv]) {
				for (int i = 0; i < 5; i++)
					shat[icv][i] = s[icv][i];
				lusolv(LUDEC[icv], shat[icv], 5);
			}
		
		// Solve linelet
		for (int iSources = 0; iSources < nSources; iSources++) {
			int n = linelet_cv[iSources].size();
			
			// LU descompusition
			for (int i = 0; i < 5; i++) {
				vhat[0][i]=s[linelet_cv[iSources][0]][i];
				for (int j = 0; j < 5; j++)
					bhat[0][i][j]=Ap[nbocv_i[linelet_cv[iSources][0]]][i][j];
			}
			
			for (int i = 1; i < n; i++) {
				int cv_im1 = linelet_cv[iSources][i-1];
				int cv_i = linelet_cv[iSources][i];
				
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_im1, cv_i);
				for (int iVar = 0; iVar < 5; iVar++)
					for (int jVar = 0; jVar < 5; jVar++)
						LUbhat[i-1][iVar][jVar] = bhat[i-1][iVar][jVar];
				ludeco(LUbhat[i-1], 5); luinv(LUbhat[i-1], invbhat[i-1], 5);
				
				getMultMatrix(mhat[i], Ap[noc10], invbhat[i-1]);
				getMultMatrix(aux_0, mhat[i], Ap[noc01]);
				getSubsMatrix(bhat[i], Ap[noc11], aux_0);
				
				// forward substituton
				getMultMatrix(aux_1, mhat[i], vhat[i-1]);											
				getSubsMatrix(vhat[i], s[linelet_cv[iSources][i]], aux_1);
			}
			
			// backward substituton
			int cv_nm1 = linelet_cv[iSources][n-1];
			for (int iVar = 0; iVar < 5; iVar++)
				for (int jVar = 0; jVar < 5; jVar++)
					LUbhat[n-1][iVar][jVar] = bhat[n-1][iVar][jVar];
			ludeco(LUbhat[n-1], 5); luinv(LUbhat[n-1], invbhat[n-1], 5);
			getMultMatrix(shat[cv_nm1], invbhat[n-1], vhat[n-1]);
			
			for (int i = n-2; i >= 0; i--) {
				int cv_i = linelet_cv[iSources][i];
				int cv_ip1 = linelet_cv[iSources][i+1];
				getImplDependencyIndex(noc00, noc01, noc11, noc10, cv_i, cv_ip1);
				getMultMatrix(aux_2, Ap[noc01], shat[cv_ip1]);
				getSubsMatrix(aux_1, vhat[i], aux_2);
				getMultMatrix(shat[cv_i], invbhat[i], aux_1);
				
			}
		}
		
		UpdateCvDataStateVec(shat);
		/*********** LINE IMPLICIT ************/	
		
    // t = [A] shat...
    matTimesVecOverCVs(t, Ap, shat);
		
    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++) {
        my_buf[0] += s[icv][i]*t[icv][i];
        my_buf[1] += t[icv][i]*t[icv][i];
      }
		
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		
    omega = buf[0]/(buf[1]+1.0E-20);
		
    // update phi...
    for (int icv = 0; icv<ncv; icv++)
      for (int i = 0; i<5; i++)
        phi[icv][i] += alpha*phat[icv][i]+omega*shat[icv][i];
		
    UpdateCvDataStateVec(phi);
		
    double my_res_max = 0.0;
		
    // recompute the residual...
    matTimesVecOverCVs(res, Ap, phi);
		
    for (int icv = 0; icv<ncv; icv++) // above is LHS. residual is then...
      for (int i = 0; i<5; i++)
        res[icv][i] = rhs[icv][i]-res[icv][i];
		
    // compute normalized res
    for (int icv = 0; icv<ncv; icv++) {
      double tmp[5];
			
      for (int i = 0; i<5; i++)
        tmp[i] = res[icv][i];
			
      lusolv(LUDEC[icv], tmp, 5);
			
      for (int i = 0; i<5; i++)
        my_res_max = max(my_res_max, fabs(tmp[i]));
    }
		
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
		
    // check residual 
    if (mpi_rank==0) {
			//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;
			
      if ((res_max<=zeroAbs)||(res_max/(res0_max+1.0E-12)<=zeroRel)) done = 1;
      else if (iter>maxiter) {
				//        cout << "Warning: " << scalarName
				//         << " solveCvScalarBcgstab did not converge after " << maxiter
				//         << " iters, res_max: " << res_max << endl;
				
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
  }
	
  lout(INFO_LO)<<"\tbcgstab: iter: "<<iter<<" res_max: "<<res_max<<endl;
	
  /*    delete[] res;
   delete[] res0;
   delete[] p;
   delete[] v;
   delete[] t;
   delete[] s;
   delete[] phat;
   delete[] shat;
	 
   delete [] LU;*/
	
  // let the calling routine know if we were successful...
  return iter;
}

void UgpWithCv::getMultMatrix(double c[5][5], double a[5][5], double b[5][5]) {
	for(int i = 0; i < 5; i++) 
		for(int j = 0; j < 5; j++) {
			c[i][j] = 0.0;
			for(int k = 0; k < 5; k++) 
				c[i][j] +=  a[i][k] * b[k][j];
		}
}

void UgpWithCv::getMultMatrix(double c[5][5], double a[5][5], double **b) {
	for(int i = 0; i < 5; i++) 
		for(int j = 0; j < 5; j++) {
			c[i][j] = 0.0;
			for(int k = 0; k < 5; k++) 
				c[i][j] +=  a[i][k] * b[k][j];
		}
}

void UgpWithCv::getMultMatrix(double c[5], double a[5][5], double b[5]) {
	for(int i = 0; i < 5; i++) {
		c[i] =  0.0;
		for(int j = 0; j < 5; j++)
			c[i] +=  a[i][j] * b[j];
	}
}

void UgpWithCv::getMultMatrix(double c[5], double b[5], double a[5][5]) {
	for(int i = 0; i < 5; i++) {
		c[i] =  0.0;
		for(int j = 0; j < 5; j++)
			c[i] +=  a[i][j] * b[j];
	}
}

void UgpWithCv::getSubsMatrix(double c[5][5], double a[5][5], double b[5][5]) {
	for(int i = 0; i < 5; i++) 
		for(int j = 0; j < 5; j++)
			c[i][j] =  a[i][j] - b[i][j];
}

void UgpWithCv::getSubsMatrix(double c[5][5], double a[5][5], double **b) {
	for(int i = 0; i < 5; i++) 
		for(int j = 0; j < 5; j++)
			c[i][j] =  a[i][j] - b[i][j];
}

void UgpWithCv::getSubsMatrix(double c[5], double a[5], double b[5]) {
	for(int i = 0; i < 5; i++) 
		c[i] =  a[i] - b[i];
}

// Begin Multigrid
int UgpWithCv::solveCvVectorR5Bcgstab_mg(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
																				 const double zeroRel, const int maxiter, char *scalarName, int iMesh) {
	
	if (iMesh == 1) {
		
		// we need the following work arrays...
		static double (*res)[5] = new double[ncv_mgLevel1][5];
		static double (*res0)[5] = new double[ncv_mgLevel1][5];
		static double (*p)[5] = new double[ncv_mgLevel1][5];
		static double (*v)[5] = new double[ncv_mgLevel1][5];
		static double (*t)[5] = new double[ncv_mgLevel1][5];
		static double (*s)[5] = new double[ncv_mgLevel1][5];
		static double (*phat)[5] = new double[ncv_g_mgLevel1][5];
		static double (*shat)[5] = new double[ncv_g_mgLevel1][5];
		
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel1][5][5]; // LU decomposed diagonal for speed up
		double b[5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel1[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		
		// initialize...
		for (int icv = 0; icv<ncv_mgLevel1; icv++)
			for (int i = 0; i<5; i++) {
				p[icv][i] = 0.0;
				v[icv][i] = 0.0;
			}
		
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			for (int i = 0; i<5; i++) {
				phat[icv][i] = 0.0;
				shat[icv][i] = 0.0;
			}
		
		double rho = 1.0;
		double alpha = 1.0;
		double omega = 1.0;
		
		// calculate the residual in rhs format ...
		matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
		
		for (int icv = 0; icv<ncv_mgLevel1; icv++)
			for (int i = 0; i<5; i++) {
				res[icv][i] = (rhs[icv][i]-res[icv][i]);
				res0[icv][i] = res[icv][i];
			}
		
		// compute first residual  
		double res0_max, my_res0_max = 0.0;
		
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			double tmp[5];
			
			for (int i = 0; i<5; i++)
				tmp[i] = res[icv][i];
			
			lusolv(LUDEC[icv], tmp, 5);
			
			for (int i = 0; i<5; i++)
				my_res0_max = max(my_res0_max, fabs(tmp[i]));
		}
		
		MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
		
		// start iteration
		
		int iter = 0;
		int done = 0;
		double res_max;
		
		while (done==0) {
			
			iter++;
			
			double rho_prime = -omega*rho;
			
			double my_rho = 0.0;
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++) {
					my_rho += res[icv][i]*res0[icv][i];
				}
			
			MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
			
			double beta = alpha*rho/rho_prime;
			
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++)
					p[icv][i] = res[icv][i]-beta*(p[icv][i]-omega*v[icv][i]);
			
			// block-diagonal precon...
			// solve the system
			for (int icv = 0; icv<ncv_mgLevel1; icv++) {
				for (int i = 0; i<5; i++)
					phat[icv][i] = p[icv][i];
				lusolv(LUDEC[icv], phat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(phat, iMesh);
			
			// v = [A]{phat}
			matTimesVecOverCVs_mg(v, Ap, phat, iMesh);
			
			double my_gamma = 0.0;
			
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++)
					my_gamma += v[icv][i]*res0[icv][i];
			
			double gamma;
			MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			//    gamma = my_gamma;
			
			if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;
			
			alpha = rho/gamma;
			
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++)
					s[icv][i] = res[icv][i]-alpha*v[icv][i];
			
			// diagonal precon...
			for (int icv = 0; icv<ncv_mgLevel1; icv++) {
				for (int i = 0; i<5; i++)
					shat[icv][i] = s[icv][i];
				lusolv(LUDEC[icv], shat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(shat, iMesh);
			
			// t = [A] shat...
			matTimesVecOverCVs_mg(t, Ap, shat, iMesh);
			
			double my_buf[2];
			my_buf[0] = my_buf[1] = 0.0;
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++) {
					my_buf[0] += s[icv][i]*t[icv][i];
					my_buf[1] += t[icv][i]*t[icv][i];
				}
			
			double buf[2];
			MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			omega = buf[0]/(buf[1]+1.0E-20);
			
			// update phi...
			for (int icv = 0; icv<ncv_mgLevel1; icv++)
				for (int i = 0; i<5; i++)
					phi[icv][i] += alpha*phat[icv][i]+omega*shat[icv][i];
			
//			UpdateCvDataStateVec_mg(phi, iMesh);
			
			double my_res_max = 0.0;
			
			// recompute the residual...
			matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
			
			for (int icv = 0; icv<ncv_mgLevel1; icv++) // above is LHS. residual is then...
				for (int i = 0; i<5; i++)
					res[icv][i] = rhs[icv][i]-res[icv][i];
			
			// compute normalized res
			for (int icv = 0; icv<ncv_mgLevel1; icv++) {
				double tmp[5];
				
				for (int i = 0; i<5; i++)
					tmp[i] = res[icv][i];
				
				lusolv(LUDEC[icv], tmp, 5);
				
				for (int i = 0; i<5; i++)
					my_res_max = max(my_res_max, fabs(tmp[i]));
			}
			
			MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
			
			// check residual 
			if (mpi_rank==0) {
				//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;
				
				if ((res_max<=zeroAbs)||(res_max/(res0_max+1.0E-12)<=zeroRel)) done = 1;
				else if (iter>maxiter) {
					/*        cout << "Warning: " << scalarName
					 << " solveCvScalarBcgstab did not converge after " << maxiter
					 << " iters, res_max: " << res_max << endl;*/
					
					done = 2;
				}
			}
			MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
		}
		
		lout(INFO_LO)<<"\tbcgstab: iter: "<<iter<<" res_max: "<<res_max<<endl;
		
		/*    delete[] res;
		 delete[] res0;
		 delete[] p;
		 delete[] v;
		 delete[] t;
		 delete[] s;
		 delete[] phat;
		 delete[] shat;
		 
		 delete [] LU;*/
		
		// let the calling routine know if we were successful...
		return (done==1);
		
	}
	
	if (iMesh == 2) {
		
		// we need the following work arrays...
		static double (*res)[5] = new double[ncv_mgLevel2][5];
		static double (*res0)[5] = new double[ncv_mgLevel2][5];
		static double (*p)[5] = new double[ncv_mgLevel2][5];
		static double (*v)[5] = new double[ncv_mgLevel2][5];
		static double (*t)[5] = new double[ncv_mgLevel2][5];
		static double (*s)[5] = new double[ncv_mgLevel2][5];
		static double (*phat)[5] = new double[ncv_g_mgLevel2][5];
		static double (*shat)[5] = new double[ncv_g_mgLevel2][5];
		
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel2][5][5]; // LU decomposed diagonal for speed up
		double b[5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel2[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		
		// initialize...
		for (int icv = 0; icv<ncv_mgLevel2; icv++)
			for (int i = 0; i<5; i++) {
				p[icv][i] = 0.0;
				v[icv][i] = 0.0;
			}
		
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			for (int i = 0; i<5; i++) {
				phat[icv][i] = 0.0;
				shat[icv][i] = 0.0;
			}
		
		double rho = 1.0;
		double alpha = 1.0;
		double omega = 1.0;
		
		// calculate the residual in rhs format ...
		matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
		
		for (int icv = 0; icv<ncv_mgLevel2; icv++)
			for (int i = 0; i<5; i++) {
				res[icv][i] = (rhs[icv][i]-res[icv][i]);
				res0[icv][i] = res[icv][i];
			}
		
		// compute first residual  
		double res0_max, my_res0_max = 0.0;
		
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			double tmp[5];
			
			for (int i = 0; i<5; i++)
				tmp[i] = res[icv][i];
			
			lusolv(LUDEC[icv], tmp, 5);
			
			for (int i = 0; i<5; i++)
				my_res0_max = max(my_res0_max, fabs(tmp[i]));
		}
		
		MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
		
		// start iteration
		
		int iter = 0;
		int done = 0;
		double res_max;
		
		while (done==0) {
			
			iter++;
			
			double rho_prime = -omega*rho;
			
			double my_rho = 0.0;
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++) {
					my_rho += res[icv][i]*res0[icv][i];
				}
			
			MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
			
			double beta = alpha*rho/rho_prime;
			
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++)
					p[icv][i] = res[icv][i]-beta*(p[icv][i]-omega*v[icv][i]);
			
			// block-diagonal precon...
			// solve the system
			for (int icv = 0; icv<ncv_mgLevel2; icv++) {
				for (int i = 0; i<5; i++)
					phat[icv][i] = p[icv][i];
				lusolv(LUDEC[icv], phat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(phat, iMesh);
			
			// v = [A]{phat}
			matTimesVecOverCVs_mg(v, Ap, phat, iMesh);
			
			double my_gamma = 0.0;
			
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++)
					my_gamma += v[icv][i]*res0[icv][i];
			
			double gamma;
			MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			//    gamma = my_gamma;
			
			if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;
			
			alpha = rho/gamma;
			
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++)
					s[icv][i] = res[icv][i]-alpha*v[icv][i];
			
			// diagonal precon...
			for (int icv = 0; icv<ncv_mgLevel2; icv++) {
				for (int i = 0; i<5; i++)
					shat[icv][i] = s[icv][i];
				lusolv(LUDEC[icv], shat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(shat, iMesh);
			
			// t = [A] shat...
			matTimesVecOverCVs_mg(t, Ap, shat, iMesh);
			
			double my_buf[2];
			my_buf[0] = my_buf[1] = 0.0;
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++) {
					my_buf[0] += s[icv][i]*t[icv][i];
					my_buf[1] += t[icv][i]*t[icv][i];
				}
			
			double buf[2];
			MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			omega = buf[0]/(buf[1]+1.0E-20);
			
			// update phi...
			for (int icv = 0; icv<ncv_mgLevel2; icv++)
				for (int i = 0; i<5; i++)
					phi[icv][i] += alpha*phat[icv][i]+omega*shat[icv][i];
			
//			UpdateCvDataStateVec_mg(phi, iMesh);
			
			double my_res_max = 0.0;
			
			// recompute the residual...
			matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
			
			for (int icv = 0; icv<ncv_mgLevel2; icv++) // above is LHS. residual is then...
				for (int i = 0; i<5; i++)
					res[icv][i] = rhs[icv][i]-res[icv][i];
			
			// compute normalized res
			for (int icv = 0; icv<ncv_mgLevel2; icv++) {
				double tmp[5];
				
				for (int i = 0; i<5; i++)
					tmp[i] = res[icv][i];
				
				lusolv(LUDEC[icv], tmp, 5);
				
				for (int i = 0; i<5; i++)
					my_res_max = max(my_res_max, fabs(tmp[i]));
			}
			
			MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
			
			// check residual 
			if (mpi_rank==0) {
				//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;
				
				if ((res_max<=zeroAbs)||(res_max/(res0_max+1.0E-12)<=zeroRel)) done = 1;
				else if (iter>maxiter) {
					/*        cout << "Warning: " << scalarName
					 << " solveCvScalarBcgstab did not converge after " << maxiter
					 << " iters, res_max: " << res_max << endl;*/
					
					done = 2;
				}
			}
			MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
		}
		
		lout(INFO_LO)<<"\tbcgstab: iter: "<<iter<<" res_max: "<<res_max<<endl;
		
		/*    delete[] res;
		 delete[] res0;
		 delete[] p;
		 delete[] v;
		 delete[] t;
		 delete[] s;
		 delete[] phat;
		 delete[] shat;
		 
		 delete [] LU;*/
		
		// let the calling routine know if we were successful...
		return (done==1);
	}
	
	if (iMesh == 3) {
		
		// we need the following work arrays...
		static double (*res)[5] = new double[ncv_mgLevel3][5];
		static double (*res0)[5] = new double[ncv_mgLevel3][5];
		static double (*p)[5] = new double[ncv_mgLevel3][5];
		static double (*v)[5] = new double[ncv_mgLevel3][5];
		static double (*t)[5] = new double[ncv_mgLevel3][5];
		static double (*s)[5] = new double[ncv_mgLevel3][5];
		static double (*phat)[5] = new double[ncv_g_mgLevel3][5];
		static double (*shat)[5] = new double[ncv_g_mgLevel3][5];
		
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel3][5][5]; // LU decomposed diagonal for speed up
		double b[5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel3[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		
		// initialize...
		for (int icv = 0; icv<ncv_mgLevel3; icv++)
			for (int i = 0; i<5; i++) {
				p[icv][i] = 0.0;
				v[icv][i] = 0.0;
			}
		
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			for (int i = 0; i<5; i++) {
				phat[icv][i] = 0.0;
				shat[icv][i] = 0.0;
			}
		
		double rho = 1.0;
		double alpha = 1.0;
		double omega = 1.0;
		
		// calculate the residual in rhs format ...
		matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
		
		for (int icv = 0; icv<ncv_mgLevel3; icv++)
			for (int i = 0; i<5; i++) {
				res[icv][i] = (rhs[icv][i]-res[icv][i]);
				res0[icv][i] = res[icv][i];
			}
		
		// compute first residual  
		double res0_max, my_res0_max = 0.0;
		
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			double tmp[5];
			
			for (int i = 0; i<5; i++)
				tmp[i] = res[icv][i];
			
			lusolv(LUDEC[icv], tmp, 5);
			
			for (int i = 0; i<5; i++)
				my_res0_max = max(my_res0_max, fabs(tmp[i]));
		}
		
		MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
		
		// start iteration
		
		int iter = 0;
		int done = 0;
		double res_max;
		
		while (done==0) {
			
			iter++;
			
			double rho_prime = -omega*rho;
			
			double my_rho = 0.0;
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++) {
					my_rho += res[icv][i]*res0[icv][i];
				}
			
			MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			if (fabs(rho_prime)<1.0E-20) rho_prime = 1.0E-20;
			
			double beta = alpha*rho/rho_prime;
			
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++)
					p[icv][i] = res[icv][i]-beta*(p[icv][i]-omega*v[icv][i]);
			
			// block-diagonal precon...
			// solve the system
			for (int icv = 0; icv<ncv_mgLevel3; icv++) {
				for (int i = 0; i<5; i++)
					phat[icv][i] = p[icv][i];
				lusolv(LUDEC[icv], phat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(phat, iMesh);
			
			// v = [A]{phat}
			matTimesVecOverCVs_mg(v, Ap, phat, iMesh);
			
			double my_gamma = 0.0;
			
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++)
					my_gamma += v[icv][i]*res0[icv][i];
			
			double gamma;
			MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			//    gamma = my_gamma;
			
			if (fabs(gamma)<1.0E-20) gamma = 1.0E-20;
			
			alpha = rho/gamma;
			
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++)
					s[icv][i] = res[icv][i]-alpha*v[icv][i];
			
			// diagonal precon...
			for (int icv = 0; icv<ncv_mgLevel3; icv++) {
				for (int i = 0; i<5; i++)
					shat[icv][i] = s[icv][i];
				lusolv(LUDEC[icv], shat[icv], 5);
			}
			
//			UpdateCvDataStateVec_mg(shat, iMesh);
			
			// t = [A] shat...
			matTimesVecOverCVs_mg(t, Ap, shat, iMesh);
			
			double my_buf[2];
			my_buf[0] = my_buf[1] = 0.0;
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++) {
					my_buf[0] += s[icv][i]*t[icv][i];
					my_buf[1] += t[icv][i]*t[icv][i];
				}
			
			double buf[2];
			MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
			
			omega = buf[0]/(buf[1]+1.0E-20);
			
			// update phi...
			for (int icv = 0; icv<ncv_mgLevel3; icv++)
				for (int i = 0; i<5; i++)
					phi[icv][i] += alpha*phat[icv][i]+omega*shat[icv][i];
			
//			UpdateCvDataStateVec_mg(phi, iMesh);
			
			double my_res_max = 0.0;
			
			// recompute the residual...
			matTimesVecOverCVs_mg(res, Ap, phi, iMesh);
			
			for (int icv = 0; icv<ncv_mgLevel3; icv++) // above is LHS. residual is then...
				for (int i = 0; i<5; i++)
					res[icv][i] = rhs[icv][i]-res[icv][i];
			
			// compute normalized res
			for (int icv = 0; icv<ncv_mgLevel3; icv++) {
				double tmp[5];
				
				for (int i = 0; i<5; i++)
					tmp[i] = res[icv][i];
				
				lusolv(LUDEC[icv], tmp, 5);
				
				for (int i = 0; i<5; i++)
					my_res_max = max(my_res_max, fabs(tmp[i]));
			}
			
			MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
			
			// check residual 
			if (mpi_rank==0) {
				//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;
				
				if ((res_max<=zeroAbs)||(res_max/(res0_max+1.0E-12)<=zeroRel)) done = 1;
				else if (iter>maxiter) {
					/*        cout << "Warning: " << scalarName
					 << " solveCvScalarBcgstab did not converge after " << maxiter
					 << " iters, res_max: " << res_max << endl;*/
					
					done = 2;
				}
			}
			MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
		}
		
		lout(INFO_LO)<<"\tbcgstab: iter: "<<iter<<" res_max: "<<res_max<<endl;
		
		/*    delete[] res;
		 delete[] res0;
		 delete[] p;
		 delete[] v;
		 delete[] t;
		 delete[] s;
		 delete[] phat;
		 delete[] shat;
		 
		 delete [] LU;*/
		
		// let the calling routine know if we were successful...
		return (done==1);
	}
	
}

void UgpWithCv::LowerProduct_mg(double(*LowProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh) {
	if (iMesh == 1) {
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				LowProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel1[icv]+1;
			int noc_l = nbocv_i_mgLevel1[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel1[noc];
				if (icv_nbr < icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							LowProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				LowProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel2[icv]+1;
			int noc_l = nbocv_i_mgLevel2[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel2[noc];
				if (icv_nbr < icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							LowProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
	if (iMesh == 3) {
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				LowProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel3[icv]+1;
			int noc_l = nbocv_i_mgLevel3[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel3[noc];
				if (icv_nbr < icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							LowProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
}

void UgpWithCv::UpperProduct_mg(double(*UpperProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh) {
	if (iMesh == 1) {
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				UpperProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel1[icv]+1;
			int noc_l = nbocv_i_mgLevel1[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel1[noc];
				if (icv_nbr > icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							UpperProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				UpperProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel2[icv]+1;
			int noc_l = nbocv_i_mgLevel2[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel2[noc];
				if (icv_nbr > icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							UpperProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
	if (iMesh == 3) {
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				UpperProd[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel3[icv]+1;
			int noc_l = nbocv_i_mgLevel3[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) {
				int icv_nbr = nbocv_v_mgLevel3[noc];
				if (icv_nbr > icv) {
					for (int i = 0; i<5; i++)
						for (int j = 0; j<5; j++)
							UpperProd[icv][i] += Ap[noc][i][j]*Phi[icv_nbr][j];
				}
			}
		}
	}
}

void UgpWithCv::DiagonalProduct_mg(double(*DiagProd)[5], double(*Ap)[5][5], double(*Phi)[5], int iMesh) {
	if (iMesh == 1) {
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				DiagProd[icv][i] = 0.0;
			
			int icv_nbr = nbocv_v_mgLevel1[nbocv_i_mgLevel1[icv]];
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					DiagProd[icv][i] += Ap[nbocv_i_mgLevel1[icv]][i][j]*Phi[icv_nbr][j];
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				DiagProd[icv][i] = 0.0;
			
			int icv_nbr = nbocv_v_mgLevel2[nbocv_i_mgLevel2[icv]];
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					DiagProd[icv][i] += Ap[nbocv_i_mgLevel2[icv]][i][j]*Phi[icv_nbr][j];
		}
	}
	if (iMesh == 3) {
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				DiagProd[icv][i] = 0.0;
			
			int icv_nbr = nbocv_v_mgLevel3[nbocv_i_mgLevel3[icv]];
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					DiagProd[icv][i] += Ap[nbocv_i_mgLevel3[icv]][i][j]*Phi[icv_nbr][j];
		}
	}
}

int UgpWithCv::solveCvVectorR5Lusgs_mg(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5], const double zeroAbs,
																			 const double zeroRel, const int maxiter, char *scalarName, int iMesh) {
	if (iMesh == 1) {
		static double (*LowProd)[5] = new double[ncv_mgLevel1][5];
		static double (*UpperProd)[5] = new double[ncv_mgLevel1][5];
		static double (*DiagProd)[5] = new double[ncv_mgLevel1][5];
		static double (*AuxVector)[5] = new double[ncv_mgLevel1][5];
		static double (*x_n)[5] = new double[ncv_mgLevel1][5];
		double tmp[5];
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel1][5][5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel1[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		// First part of the symmetric iteration: (D+L).x^* = b
		// compute Lx^*
		LowerProduct_mg(LowProd, Ap, phi, 1);  
		
		// compute aux_vector = b-L.x^*, and solve Dx^* = aux_vector, using a LU method
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++) 
				tmp[i] = rhs[icv][i]-LowProd[icv][i];
			lusolv(LUDEC[icv], tmp, 5);
			for (int i = 0; i<5; i++) 
				x_n[icv][i] = tmp[i];
		}
		
		// Second part of the symmetric iteration: (D+U).x^1 = D.x^*
		// compute D.x^*
		DiagonalProduct_mg(DiagProd, Ap, x_n, 1); 
		
		// compute aux_vector = D.x^*
		for (int icv = 0; icv<ncv_mgLevel1; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] = DiagProd[icv][i];
		
		// compute U.x^1
		UpperProduct_mg(UpperProd, Ap, x_n, 1);
		
		// compute aux_vector = D.x^*-U.x^1
		for (int icv = 0; icv<ncv_mgLevel1; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] -= UpperProd[icv][i];
		
		// solve D.x^1 = aux_vector, using a LU method
		double res_max;
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			lusolv(LUDEC[icv], AuxVector[icv], 5);
			for (int i = 0; i<5; i++) {
				res_max = max(res_max, fabs(AuxVector[icv][i]));
				x_n[icv][i] = AuxVector[icv][i];
			}
		}
		
		// update phi...
		for (int icv = 0; icv<ncv_mgLevel1; icv++)
			for (int i = 0; i<5; i++)
				phi[icv][i] = x_n[icv][i];
	}
	
	if (iMesh == 2) {
		static double (*LowProd)[5] = new double[ncv_mgLevel2][5];
		static double (*UpperProd)[5] = new double[ncv_mgLevel2][5];
		static double (*DiagProd)[5] = new double[ncv_mgLevel2][5];
		static double (*AuxVector)[5] = new double[ncv_mgLevel2][5];
		static double (*x_n)[5] = new double[ncv_mgLevel2][5];
		double tmp[5];
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel2][5][5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel2[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		// First part of the symmetric iteration: (D+L).x^* = b
		// compute Lx^*
		LowerProduct_mg(LowProd, Ap, phi, 2);  
		
		// compute aux_vector = b-L.x^*, and solve Dx^* = aux_vector, using a LU method
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++) 
				tmp[i] = rhs[icv][i]-LowProd[icv][i];
			lusolv(LUDEC[icv], tmp, 5);
			for (int i = 0; i<5; i++) 
				x_n[icv][i] = tmp[i];
		}
		
		// Second part of the symmetric iteration: (D+U).x^1 = D.x^*
		// compute D.x^*
		DiagonalProduct_mg(DiagProd, Ap, x_n, 2); 
		
		// compute aux_vector = D.x^*
		for (int icv = 0; icv<ncv_mgLevel2; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] = DiagProd[icv][i];
		
		// compute U.x^1
		UpperProduct_mg(UpperProd, Ap, x_n, 2);
		
		// compute aux_vector = D.x^*-U.x^1
		for (int icv = 0; icv<ncv_mgLevel2; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] -= UpperProd[icv][i];
		
		// solve D.x^1 = aux_vector, using a LU method
		double res_max;
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			lusolv(LUDEC[icv], AuxVector[icv], 5);
			for (int i = 0; i<5; i++) {
				res_max = max(res_max, fabs(AuxVector[icv][i]));
				x_n[icv][i] = AuxVector[icv][i];
			}
		}
		
		// update phi...
		for (int icv = 0; icv<ncv_mgLevel2; icv++)
			for (int i = 0; i<5; i++)
				phi[icv][i] = x_n[icv][i];
	}
	if (iMesh == 3) {
		static double (*LowProd)[5] = new double[ncv_mgLevel3][5];
		static double (*UpperProd)[5] = new double[ncv_mgLevel3][5];
		static double (*DiagProd)[5] = new double[ncv_mgLevel3][5];
		static double (*AuxVector)[5] = new double[ncv_mgLevel3][5];
		static double (*x_n)[5] = new double[ncv_mgLevel3][5];
		double tmp[5];
		static double (*LUDEC)[5][5] = new double[ncv_mgLevel3][5][5];
		
		// LU decompose diagonal and store
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				for (int j = 0; j<5; j++)
					LUDEC[icv][i][j] = Ap[nbocv_i_mgLevel3[icv]][i][j];
			ludeco(LUDEC[icv], 5);
		}
		
		// First part of the symmetric iteration: (D+L).x^* = b
		// compute Lx^*
		LowerProduct_mg(LowProd, Ap, phi, 3);  
		
		// compute aux_vector = b-L.x^*, and solve Dx^* = aux_vector, using a LU method
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++) 
				tmp[i] = rhs[icv][i]-LowProd[icv][i];
			lusolv(LUDEC[icv], tmp, 5);
			for (int i = 0; i<5; i++) 
				x_n[icv][i] = tmp[i];
		}
		
		// Second part of the symmetric iteration: (D+U).x^1 = D.x^*
		// compute D.x^*
		DiagonalProduct_mg(DiagProd, Ap, x_n, 3); 
		
		// compute aux_vector = D.x^*
		for (int icv = 0; icv<ncv_mgLevel3; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] = DiagProd[icv][i];
		
		// compute U.x^1
		UpperProduct_mg(UpperProd, Ap, x_n, 3);
		
		// compute aux_vector = D.x^*-U.x^1
		for (int icv = 0; icv<ncv_mgLevel3; icv++)
			for (int i = 0; i<5; i++)
				AuxVector[icv][i] -= UpperProd[icv][i];
		
		// solve D.x^1 = aux_vector, using a LU method
		double res_max;
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			lusolv(LUDEC[icv], AuxVector[icv], 5);
			for (int i = 0; i<5; i++) {
				res_max = max(res_max, fabs(AuxVector[icv][i]));
				x_n[icv][i] = AuxVector[icv][i];
			}
		}
		
		// update phi...
		for (int icv = 0; icv<ncv_mgLevel3; icv++)
			for (int i = 0; i<5; i++)
				phi[icv][i] = x_n[icv][i];
	}
	return (1);
}

void UgpWithCv::matTimesVecOverCVs_mg(double(*res)[5], double(*Ap)[5][5], double(*phi)[5], int iMesh) {
	if (iMesh == 1) {
		for (int icv = 0; icv<ncv_mgLevel1; icv++) {
			for (int i = 0; i<5; i++)
				res[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel1[icv];
			int noc_l = nbocv_i_mgLevel1[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) // cv diagonal + neighbors...
			{
				const int icv_nbr = nbocv_v_mgLevel1[noc];
				
				for (int i = 0; i<5; i++)
					for (int j = 0; j<5; j++)
						res[icv][i] += Ap[noc][i][j]*phi[icv_nbr][j];
			}
		}
	}
	if (iMesh == 2) {
		for (int icv = 0; icv<ncv_mgLevel2; icv++) {
			for (int i = 0; i<5; i++)
				res[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel2[icv];
			int noc_l = nbocv_i_mgLevel2[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) // cv diagonal + neighbors...
			{
				const int icv_nbr = nbocv_v_mgLevel2[noc];
				
				for (int i = 0; i<5; i++)
					for (int j = 0; j<5; j++)
						res[icv][i] += Ap[noc][i][j]*phi[icv_nbr][j];
			}
		}
	}	
	if (iMesh == 3) {
		for (int icv = 0; icv<ncv_mgLevel3; icv++) {
			for (int i = 0; i<5; i++)
				res[icv][i] = 0.0;
			
			int noc_f = nbocv_i_mgLevel3[icv];
			int noc_l = nbocv_i_mgLevel3[icv+1]-1;
			for (int noc = noc_f; noc<=noc_l; noc++) // cv diagonal + neighbors...
			{
				const int icv_nbr = nbocv_v_mgLevel3[noc];
				
				for (int i = 0; i<5; i++)
					for (int j = 0; j<5; j++)
						res[icv][i] += Ap[noc][i][j]*phi[icv_nbr][j];
			}
		}
	}	
	
}

void UgpWithCv::UpdateCvDataStateVec_mg(double(*phi)[5], int iMesh) 
{
	if (iMesh == 1) {
		static double *scal = new double[ncv_g_mgLevel1];
		static double (*vec)[3] = new double[ncv_g_mgLevel1][3];
		
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			scal[icv] = phi[icv][0];
		updateCvData_mg(scal, REPLACE_DATA, 1);
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			phi[icv][0] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			scal[icv] = phi[icv][4];
		updateCvData_mg(scal, REPLACE_DATA, 1);
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			phi[icv][4] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			for (int i = 0; i<3; i++)
				vec[icv][i] = phi[icv][i+1];
		updateCvData_mg(vec, REPLACE_ROTATE_DATA, 1);
		for (int icv = 0; icv<ncv_g_mgLevel1; icv++)
			for (int i = 0; i<3; i++)
				phi[icv][i+1] = vec[icv][i];
	}
	if (iMesh == 2) {
		static double *scal = new double[ncv_g_mgLevel2];
		static double (*vec)[3] = new double[ncv_g_mgLevel2][3];
		
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			scal[icv] = phi[icv][0];
		updateCvData_mg(scal, REPLACE_DATA, 2);
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			phi[icv][0] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			scal[icv] = phi[icv][4];
		updateCvData_mg(scal, REPLACE_DATA, 2);
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			phi[icv][4] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			for (int i = 0; i<3; i++)
				vec[icv][i] = phi[icv][i+1];
		updateCvData_mg(vec, REPLACE_ROTATE_DATA, 2);
		for (int icv = 0; icv<ncv_g_mgLevel2; icv++)
			for (int i = 0; i<3; i++)
				phi[icv][i+1] = vec[icv][i];
	}
	if (iMesh == 3) {
		static double *scal = new double[ncv_g_mgLevel3];
		static double (*vec)[3] = new double[ncv_g_mgLevel3][3];
		
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			scal[icv] = phi[icv][0];
		updateCvData_mg(scal, REPLACE_DATA, 3);
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			phi[icv][0] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			scal[icv] = phi[icv][4];
		updateCvData_mg(scal, REPLACE_DATA, 3);
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			phi[icv][4] = scal[icv];
		
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			for (int i = 0; i<3; i++)
				vec[icv][i] = phi[icv][i+1];
		updateCvData_mg(vec, REPLACE_ROTATE_DATA, 3);
		for (int icv = 0; icv<ncv_g_mgLevel3; icv++)
			for (int i = 0; i<3; i++)
				phi[icv][i+1] = vec[icv][i];
	}	
}
