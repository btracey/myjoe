#include "UgpWithCvFake.h"
#include <math.h>

void UgpWithCvFake::addGhostCvs() {

  // we have to add ghosts for both internal (inter-processor) and
  // periodic boundaries...
  
  // include fake cv's too
  cerr << "Error: include fake cvs too, even if you don't need them." << endl;
  cerr << "Use addGhostAndFakeCvs" << endl;
  throw(-1);
  
  ncv_g = ncv;
  for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); fa_prcomm != facePrcommList.end(); fa_prcomm++) {
    Prcomm * cv_prcomm = getCvPrcomm(fa_prcomm->getNbrRank());
    // size the pack and unpack vecs the same as the faces - they could potentially be less
    // but for now...
    cv_prcomm->packIndexVec.resize(fa_prcomm->packIndexVec.size());
    for (int i = 0; i < fa_prcomm->packIndexVec.size(); i++) {
      int ifa = fa_prcomm->packIndexVec[i];
      // the internal face is always valid...
      assert((cvofa[ifa][0] >= 0) && (cvofa[ifa][0] < nfa));
      cv_prcomm->packIndexVec[i] = cvofa[ifa][0];
    }
    cv_prcomm->unpackIndexVec.resize(fa_prcomm->unpackIndexVec.size());
    for (int i = 0; i < fa_prcomm->unpackIndexVec.size(); i++) {
      int ifa = fa_prcomm->unpackIndexVec[i];
      // the ghost face should not be valid (yet)...
      assert(cvofa[ifa][1] == -1);
      cv_prcomm->unpackIndexVec[i] = ncv_g;
      cvofa[ifa][1] = ncv_g;
      ncv_g++;
    }
    // for proper transformations at periodic boundaries, it is neccessary (for the packing atleast)
    // to specify the Range's that divide up the Prcomm. Since this is an exact copy, this is trivial -
    // just copy the Face's packRanges...
    for (list<Range>::iterator ri = fa_prcomm->packRangeList.begin(); ri != fa_prcomm->packRangeList.end(); ri++) {
      cv_prcomm->addPackRange(ri->getIndexFirst(), ri->getIndexLast(), ri->getBits(), ri->getFlag(), ri->dxyz);
    }
  }

  assert( ncv_g >= ncv );

  // --------------------------------------------------
  // adjust data sizes to include ghosts in CV_DATA
  // --------------------------------------------------

  // reallocate the cv_flag for this new cv size...
  if (ncv_g > ncv) {
    delete[] cv_flag;
    cv_flag = new int[ncv_g];
  }

  // double scalars: [ncv]...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g > ncv) {
        double * tmp = new double[ncv_g];
        for (int icv = 0; icv < ncv; icv++)
          tmp[icv] = (*(data->ptr))[icv];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_DATA);
    }
  }

  // double vectors: [ncv][3]...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g > ncv) {
        double (*tmp)[3] = new double[ncv_g][3];
        for (int icv = 0; icv < ncv; icv++)
          for (int j = 0; j < 3; j++)
            tmp[icv][j] = (*(data->ptr))[icv][j];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

  // double tensors: [ncv][3][3]...
  for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data != doubleTensorList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_g > ncv) {
        double (*tmp)[3][3] = new double[ncv_g][3][3];
        for (int icv = 0; icv < ncv; icv++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
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

void UgpWithCvFake::addGhostAndFakeCvs() {

  // we have to add ghosts for both internal (inter-processor) and
  // periodic boundaries, AND fakes for local boundary faces...
  
  ncv_g = ncv;
  for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); fa_prcomm != facePrcommList.end(); fa_prcomm++) {
    Prcomm * cv_prcomm = getCvPrcomm(fa_prcomm->getNbrRank());
    // size the pack and unpack vecs the same as the faces - they could potentially be less
    // but for now...
    cv_prcomm->packIndexVec.resize(fa_prcomm->packIndexVec.size());
    for (int i = 0; i < fa_prcomm->packIndexVec.size(); i++) {
      int ifa = fa_prcomm->packIndexVec[i];
      // the internal face is always valid...
      assert((cvofa[ifa][0] >= 0) && (cvofa[ifa][0] < nfa));
      cv_prcomm->packIndexVec[i] = cvofa[ifa][0];
    }
    cv_prcomm->unpackIndexVec.resize(fa_prcomm->unpackIndexVec.size());
    for (int i = 0; i < fa_prcomm->unpackIndexVec.size(); i++) {
      int ifa = fa_prcomm->unpackIndexVec[i];
      // the ghost face should not be valid (yet)...
      assert(cvofa[ifa][1] == -1);
      cv_prcomm->unpackIndexVec[i] = ncv_g;
      cvofa[ifa][1] = ncv_g;
      ncv_g++;
    }
    // for proper transformations at periodic boundaries, it is neccessary (for the packing atleast)
    // to specify the Range's that divide up the Prcomm. Since this is an exact copy, this is trivial -
    // just copy the Face's packRanges...
    for (list<Range>::iterator ri = fa_prcomm->packRangeList.begin(); ri != fa_prcomm->packRangeList.end(); ri++) {
      cv_prcomm->addPackRange(ri->getIndexFirst(), ri->getIndexLast(), ri->getBits(), ri->getFlag(), ri->dxyz);
    }
  }

  assert( ncv_g >= ncv );

  // fakes are really easy...

  ncv_gf = ncv_g;
  for (int ifa = 0; ifa < nfa_b; ++ifa) {
    assert( cvofa[ifa][1] < 0 );
    cvofa[ifa][1] = ncv_gf++;
  }

  // --------------------------------------------------
  // adjust data sizes to include ghosts in CV_DATA
  // --------------------------------------------------

  // i think is will always be larger...
  assert( ncv_gf > ncv ); 
  
  // reallocate the cv_flag for this new cv size...
  if (ncv_gf > ncv) {
    delete[] cv_flag;
    cv_flag = new int[ncv_gf];
  }

  // double scalars: [ncv]...
  for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf > ncv) {
        double * tmp = new double[ncv_gf];
        for (int icv = 0; icv < ncv; icv++)
          tmp[icv] = (*(data->ptr))[icv];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_DATA);
    }
  }

  // double vectors: [ncv][3]...
  for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data != doubleVectorList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf > ncv) {
        double (*tmp)[3] = new double[ncv_gf][3];
        for (int icv = 0; icv < ncv; icv++)
          for (int j = 0; j < 3; j++)
            tmp[icv][j] = (*(data->ptr))[icv][j];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }

  // double tensors: [ncv][3][3]...
  for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data != doubleTensorList.end(); data++) {
    if (data->getDatatype() == CV_DATA) {
      // only reallocate if this processor actually got some ghost cv's...
      if (ncv_gf > ncv) {
        double (*tmp)[3][3] = new double[ncv_gf][3][3];
        for (int icv = 0; icv < ncv; icv++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
              tmp[icv][j][k] = (*(data->ptr))[icv][j][k];
        delete[] (*(data->ptr));
        *(data->ptr) = tmp;
      }
      // fill ghost cvs...
      updateCvData(*(data->ptr), REPLACE_ROTATE_DATA);
    }
  }
}

void UgpWithCvFake::buildNbocv() {

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
  
  // rebuild cvora...
  MPI_Allgather(&ncv,1,MPI_INT,&(cvora[1]),1,MPI_INT,mpi_comm);
  cvora[0] = 0;
  for (int i = 0; i < mpi_size; i++)
    cvora[i+1] += cvora[i];
  assert( cvora[mpi_rank+1]-cvora[mpi_rank] == ncv );
  
  if (mpi_rank == 0)
    cout << "buildNbocv()" << endl;

  assert(nbocv_i == NULL);
  assert(nbocv_v == NULL);

  nbocv_i = new int[ncv + 1];
  for (int icv = 0; icv < ncv; icv++)
    nbocv_i[icv + 1] = 1; // 1 for the diagonal

  for (int iter = 0; iter < 2; iter++) {
    // skip boundary faces - these do not have icvs that
    // participate in the nbocv_i/v structure...
    for (int ifa = nfa_b; ifa < nfa; ifa++) {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert(icv1 >= 0);
      if (iter == 0) {
        nbocv_i[icv0 + 1] += 1;
        if (icv1 < ncv)
          nbocv_i[icv1 + 1] += 1;
      } else {
        nbocv_v[nbocv_i[icv0]] = icv1;
        nbocv_i[icv0] += 1;
        if (icv1 < ncv) {
          nbocv_v[nbocv_i[icv1]] = icv0;
          nbocv_i[icv1] += 1;
        }
      }
    }
    if (iter == 0) {
      nbocv_i[0] = 0;
      for (int icv = 0; icv < ncv; icv++)
        nbocv_i[icv + 1] += nbocv_i[icv];
      nbocv_s = nbocv_i[ncv];
      nbocv_v = new int[nbocv_s];
      // add the diagonal...
      for (int icv = 0; icv < ncv; icv++) {
        nbocv_v[nbocv_i[icv]] = icv;
        nbocv_i[icv] += 1;
      }
    } else {
      for (int icv = ncv; icv > 0; icv--)
        nbocv_i[icv] = nbocv_i[icv - 1];
      nbocv_i[0] = 0;
    }
  }

}

void UgpWithCvFake::calcGeometry() {

  // ======================================
  // face normal and centroid...
  // normal has area magnitude
  // ======================================

  assert(x_fa == NULL);
  x_fa = new double[nfa][3];

  assert(fa_normal == NULL);
  fa_normal = new double[nfa][3];

  for (int ifa = 0; ifa < nfa; ifa++) {
    double x_fa_approx[3];
    for (int i = 0; i < 3; i++)
      x_fa_approx[i] = 0.0;
    for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) {
      int ino = noofa_v[nof];
      for (int i = 0; i < 3; i++)
        x_fa_approx[i] += x_no[ino][i];
    }
    // normalize: divide by the number of nodes...
    for (int i = 0; i < 3; i++)
      x_fa_approx[i] /= (double) (noofa_i[ifa + 1] - noofa_i[ifa]);
    // we can compute the face normal directly with this approx centroid...
    for (int i = 0; i < 3; i++)
      fa_normal[ifa][i] = 0.0;
    int ino2 = noofa_v[noofa_i[ifa + 1] - 1]; // last node in the list
    for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) {
      int ino1 = ino2;
      ino2 = noofa_v[nof];
      double v1[3];
      double v2[3];
      for (int i = 0; i < 3; i++) {
        v1[i] = x_no[ino1][i] - x_fa_approx[i];
        v2[i] = x_no[ino2][i] - x_fa_approx[i];
      }
      fa_normal[ifa][0] += 0.5 * (v1[1] * v2[2] - v1[2] * v2[1]);
      fa_normal[ifa][1] += 0.5 * (v1[2] * v2[0] - v1[0] * v2[2]);
      fa_normal[ifa][2] += 0.5 * (v1[0] * v2[1] - v1[1] * v2[0]);
    }
    // now we can improve the centroid...
    double area_sum = 0.0;
    for (int i = 0; i < 3; i++)
      x_fa[ifa][i] = 0.0;
    // loop through the edges again to get the normal pieces...
    ino2 = noofa_v[noofa_i[ifa + 1] - 1];
    for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) {
      int ino1 = ino2;
      ino2 = noofa_v[nof];
      double v1[3];
      double v2[3];
      for (int i = 0; i < 3; i++) {
        v1[i] = x_no[ino1][i] - x_fa_approx[i];
        v2[i] = x_no[ino2][i] - x_fa_approx[i];
      }
      double this_area = fa_normal[ifa][0] * (v1[1] * v2[2] - v1[2] * v2[1]) + fa_normal[ifa][1] * (v1[2] * v2[0]
          - v1[0] * v2[2]) + fa_normal[ifa][2] * (v1[0] * v2[1] - v1[1] * v2[0]);
      //assert(this_area > 0.0);
      if (this_area <= 0.0)
        cout << " Warning: area negative -> " << this_area << endl;
      area_sum += this_area;
      for (int i = 0; i < 3; i++)
        x_fa[ifa][i] += this_area * (x_fa_approx[i] + x_no[ino1][i] + x_no[ino2][i]);
    }
    for (int i = 0; i < 3; i++)
      x_fa[ifa][i] /= 3.0 * area_sum;

  }

  // ======================================
  // x_cv, cv_volume in ghost as well...
  // ======================================

  x_cv = new double[ncv_gf][3];
  cv_volume = new double[ncv_gf];

  for (int icv = 0; icv < ncv; icv++) {
    // compute an approximate center as the mean of the face x_fa's...
    double x_cv_approx[3];
    for (int i = 0; i < 3; i++)
      x_cv_approx[i] = 0.0;
    for (int foc = faocv_i[icv]; foc < faocv_i[icv + 1]; foc++) {
      int ifa = faocv_v[foc];
      for (int i = 0; i < 3; i++)
        x_cv_approx[i] += x_fa[ifa][i];
    }
    // divide by the number of faces...
    for (int i = 0; i < 3; i++)
      x_cv_approx[i] /= (double) (faocv_i[icv + 1] - faocv_i[icv]);
    // add tet volumes...
    cv_volume[icv] = 0.0;
    for (int i = 0; i < 3; i++)
      x_cv[icv][i] = 0.0;
    // loop on the faces of this cv...
    for (int foc = faocv_i[icv]; foc < faocv_i[icv + 1]; foc++) {
      int ifa = faocv_v[foc];
      double v1[3];
      for (int i = 0; i < 3; i++)
        v1[i] = x_fa[ifa][i] - x_cv_approx[i];
      // check if the face is inward or outward wrt this cv...
      if (cvofa[ifa][0] == icv) {
        // face is outward - loop through edges in forward direction...
        int ino2 = noofa_v[noofa_i[ifa + 1] - 1];
        for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; nof++) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i < 3; i++) {
            v2[i] = x_no[ino1][i] - x_cv_approx[i];
            v3[i] = x_no[ino2][i] - x_cv_approx[i];
          }
          // 2 nodes, the face, and the approx cv form a tet...
          double this_volume = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
                             + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2])
                             + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
          assert(this_volume > 0.0); // check on the grid ordering/right-handedness
          cv_volume[icv] += this_volume;
          for (int i = 0; i < 3; i++)
            x_cv[icv][i] += this_volume * (x_cv_approx[i] + x_fa[ifa][i] + x_no[ino1][i] + x_no[ino2][i]);
        }
      } else {
        assert(cvofa[ifa][1] == icv);
        // face is outward, loop through edges in backward direction...
        int ino2 = noofa_v[noofa_i[ifa]];
        for (int nof = noofa_i[ifa + 1] - 1; nof >= noofa_i[ifa]; nof--) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          for (int i = 0; i < 3; i++) {
            v2[i] = x_no[ino1][i] - x_cv_approx[i];
            v3[i] = x_no[ino2][i] - x_cv_approx[i];
          }
          double this_volume = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
                             + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2])
                             + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
          //assert(this_volume > 0.0);
          if (this_volume <= 0.0)
            cout << " Warning: volume negative -> " << this_volume << endl;

          cv_volume[icv] += this_volume;
          for (int i = 0; i < 3; i++)
            x_cv[icv][i] += this_volume * (x_cv_approx[i] + x_fa[ifa][i] + x_no[ino1][i] + x_no[ino2][i]);
        }
      }
    }
    // normalize both...
    for (int i = 0; i < 3; i++)
      x_cv[icv][i] /= 4.0 * cv_volume[icv];
    cv_volume[icv] /= 6.0;
  }

  // update volumes and centroids across processors...
  updateCvData(cv_volume, REPLACE_DATA); // blocking exchanges
  updateCvData(x_cv, REPLACE_TRANSLATE_DATA); // "TRANSLATE" is for any periodicity
  
  // finally fake data... 
  for (int ifa = 0; ifa < nfa_b; ifa++) {
    // copy the cv_volume...
    int icv0 = cvofa[ifa][0];
    assert( (icv0 >= 0)&&(icv0 < ncv) );
    int icv1 = cvofa[ifa][1];
    assert( (icv1 >= ncv_g)&&(icv1 < ncv_gf) );
    // copy the volume...
    cv_volume[icv1] = cv_volume[icv0];
    // and relect the cv through the normal...
    double nmag = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
			fa_normal[ifa][1]*fa_normal[ifa][1] +
			fa_normal[ifa][2]*fa_normal[ifa][2] );
    double shalf_dot_n = ( (x_fa[ifa][0] - x_cv[icv0][0])*fa_normal[ifa][0] + 
			   (x_fa[ifa][1] - x_cv[icv0][1])*fa_normal[ifa][1] + 
			   (x_fa[ifa][2] - x_cv[icv0][2])*fa_normal[ifa][2] )/nmag;
    assert( shalf_dot_n > 0.0 );
    for (int i = 0; i < 3; i++) 
      x_cv[icv1][i] = x_cv[icv0][i] + 2.0*shalf_dot_n*fa_normal[ifa][i]/nmag;
  }
  
  // check...
  double my_volume_sum = 0.0;
  for (int icv = 0; icv < ncv; icv++)
    my_volume_sum += cv_volume[icv];
  double volume_sum;
  MPI_Reduce(&my_volume_sum, &volume_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
  
  // also check s dot n...
  double my_buf[2];
  my_buf[0] = 1.0;
  my_buf[1] = 0.0;
  for (int ifa = 0; ifa < nfa_b; ifa++) {
    // the first faces are boundary faces, and should have icv1 >= ncv_g...
    assert( (cvofa[ifa][1] >= ncv_g)&&(cvofa[ifa][1] < ncv_gf) );
  }
  for (int ifa = nfa_b; ifa < nfa; ifa++) {
    // the remaining faces should have valid cvofa, with a value
    // in the ghost range for the periodic cvs, and in the
    // owned cv range for totally internal cvs...
    assert(cvofa[ifa][1] >= 0);
    if (ifa < nfa_bpi) {
      assert((cvofa[ifa][1] >= ncv) && (cvofa[ifa][1] < ncv_g));
    } else {
      assert(cvofa[ifa][1] < ncv);
    }
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    double s[3];
    double smag = 0.0;
    double nmag = 0.0;
    for (int i = 0; i < 3; i++) {
      s[i] = x_cv[icv1][i] - x_cv[icv0][i];
      smag += s[i] * s[i];
      nmag += fa_normal[ifa][i] * fa_normal[ifa][i];
    }
    smag = sqrt(smag);
    nmag = sqrt(nmag);
    //cout << "ifa, smag, nmag: " << ifa << " " << smag << " " << nmag << endl;
    double n[3];
    for (int i = 0; i < 3; i++) {
      s[i] /= smag;
      n[i] = fa_normal[ifa][i] / nmag;
    }
    double dp = s[0] * n[0] + s[1] * n[1] + s[2] * n[2];
    assert(dp > 0.0);
    my_buf[0] = min(my_buf[0], dp);
    my_buf[1] = min(my_buf[1], -dp);
  }
  
  double buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
  
  if (mpi_rank == 0)
    cout << "calcGeometry(), cv volume: " << volume_sum << ", s dot n min/max: " << buf[0] << " " << -buf[1] << endl;
  
  //MPI_Pause("done");

}

void UgpWithCvFake::calcLsgCoeff() 
{
  if (mpi_rank == 0)
    cout << "calcLsgCoeff()" << endl;
  
  // least squares gradient coeff's...
  assert(nbocv_lsg_coeff == NULL);
  nbocv_lsg_coeff = new double[nbocv_s][3];
  for (int noc = 0; noc < nbocv_s; noc++)
    for (int i = 0; i < 3; i++)
      nbocv_lsg_coeff[noc][i] = 0.0;

  assert(boundary_fa_lsg_coeff == NULL);
  boundary_fa_lsg_coeff = new double[nfa_b][3];
  for (int ifa = 0; ifa < nfa_b; ifa++)
    for (int i = 0; i < 3; i++)
      boundary_fa_lsg_coeff[ifa][i] = 0.0;

  for (int icv = 0; icv < ncv; icv++) 
  {
    double swdx2 = 0.0; // sum weight dx2, etc...
    double swdy2 = 0.0;
    double swdz2 = 0.0;
    double swdxdy = 0.0;
    double swdxdz = 0.0;
    double swdydz = 0.0;
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    // skip the diagonal in this loop...
    for (int noc = noc_f + 1; noc <= noc_l; noc++) 
    {
      int icv_nbr = nbocv_v[noc];
      double weight = 1.0;
      double dx[3];
      for (int i = 0; i < 3; i++)
        dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];
      swdx2 += weight * dx[0] * dx[0];
      swdy2 += weight * dx[1] * dx[1];
      swdz2 += weight * dx[2] * dx[2];
      swdxdy += weight * dx[0] * dx[1];
      swdxdz += weight * dx[0] * dx[2];
      swdydz += weight * dx[1] * dx[2];
    }

    // this cell may also touch some boundary faces that
    // do not have neighbours (no fake cells). Their
    // influence on the gradient has to be considered as well...
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv + 1] - 1;
    for (int foc = foc_f; foc <= foc_l; foc++) 
    {
      int ifa = faocv_v[foc];
      if (ifa < nfa_b) 
      {
        // this is a bounray face...
        assert(cvofa[ifa][0] == icv);
        int icvFake = cvofa[ifa][1];

        double weight = 1.0;
        double dx[3];
        for (int i = 0; i < 3; i++)
          //dx[i] = x_fa[ifa][i] - x_cv[icv][i];
          dx[i] = x_cv[icvFake][i] - x_cv[icv][i];

        swdx2 += weight * dx[0] * dx[0];
        swdy2 += weight * dx[1] * dx[1];
        swdz2 += weight * dx[2] * dx[2];
        swdxdy += weight * dx[0] * dx[1];
        swdxdz += weight * dx[0] * dx[2];
        swdydz += weight * dx[1] * dx[2];
      }
    }

    double denom = 2.0 * swdxdy * swdxdz * swdydz
                       + swdx2 * swdy2 * swdz2
                       - swdx2 * swdydz * swdydz
                       - swdy2 * swdxdz * swdxdz
                       - swdz2 * swdxdy * swdxdy;

    // check non-singular...
    assert(fabs(denom) > 1.0E-12 * cv_volume[icv] * cv_volume[icv]);

    // assemble coeff's...
    for (int noc = noc_f + 1; noc <= noc_l; noc++) 
    {
      int icv_nbr = nbocv_v[noc];
      double weight = 1.0;
      double dx[3];
      
      for (int i = 0; i < 3; i++)
        dx[i] = x_cv[icv_nbr][i] - x_cv[icv][i];

      // following coeff's multiply the delta across nbr - cv...
      nbocv_lsg_coeff[noc][0] = weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                                        + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                                        + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;
      nbocv_lsg_coeff[noc][1] = weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                                        + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                                        + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;
      nbocv_lsg_coeff[noc][2] = weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                                        + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                                        + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;
      // and subtract from the diagonal...
      for (int i = 0; i < 3; i++)
        nbocv_lsg_coeff[noc_f][i] -= nbocv_lsg_coeff[noc][i];
    }

    for (int foc = foc_f; foc <= foc_l; foc++) 
    {
      int ifa = faocv_v[foc];
      if (ifa < nfa_b) 
      {
        // this is a bounray face...
        assert(cvofa[ifa][0] == icv);
        int icvFake = cvofa[ifa][1];
        
        double weight = 1.0;
        double dx[3];
        for (int i = 0; i < 3; i++)
          //dx[i] = x_fa[ifa][i] - x_cv[icv][i];
          dx[i] = x_cv[icvFake][i] - x_cv[icv][i];
        
        boundary_fa_lsg_coeff[ifa][0] = weight * ((swdy2 * swdz2 - swdydz * swdydz) * dx[0]
                                                + (swdxdz * swdydz - swdxdy * swdz2) * dx[1]
                                                + (swdxdy * swdydz - swdxdz * swdy2) * dx[2]) / denom;
        boundary_fa_lsg_coeff[ifa][1] = weight * ((swdxdz * swdydz - swdxdy * swdz2) * dx[0]
                                                + (swdx2 * swdz2 - swdxdz * swdxdz) * dx[1]
                                                + (swdxdy * swdxdz - swdydz * swdx2) * dx[2]) / denom;
        boundary_fa_lsg_coeff[ifa][2] = weight * ((swdxdy * swdydz - swdxdz * swdy2) * dx[0]
                                                + (swdxdy * swdxdz - swdydz * swdx2) * dx[1]
                                                + (swdx2 * swdy2 - swdxdy * swdxdy) * dx[2]) / denom;
        // and subtract from the diagonal...
        for (int i = 0; i < 3; i++)
          nbocv_lsg_coeff[noc_f][i] -= boundary_fa_lsg_coeff[ifa][i];
      }
    }
  }
}

void UgpWithCvFake::checkLsgCoeff() 
{
  if (mpi_rank == 0)
    cout << "checkLsgCoeff()...";

  double *p = new double[ncv_gf];

  double constant_grad_p[3];
  constant_grad_p[0] = 1.12;
  constant_grad_p[1] = 2.13;
  constant_grad_p[2] = 3.14;

  double (*grad_p)[3] = new double[ncv_g][3];

  // put a linear grad in p - everywhere...
  for (int icv = 0; icv < ncv_gf; icv++)
    p[icv] = 1.1 + 1.12 * x_cv[icv][0] + 2.13 * x_cv[icv][1] + 3.14 * x_cv[icv][2];

  // compute gradients...
  for (int icv = 0; icv < ncv; icv++) 
  {
    for (int i = 0; i < 3; i++)
      grad_p[icv][i] = 0.0;
    
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    for (int noc = noc_f; noc <= noc_l; noc++) 
    {
      int icv_nbr = nbocv_v[noc];
      for (int i = 0; i < 3; i++)
        grad_p[icv][i] += nbocv_lsg_coeff[noc][i]*p[icv_nbr];
    }
  }

  /*
  // add boundaries
  for (int ifa = 0; ifa < nfa_b; ifa++) 
  {
    int icv     = cvofa[ifa][0];
    int icvFake = cvofa[ifa][1];
    for (int i=0; i<3; i++)
      grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*(p[icvFake]-p[icv]);
  }*/

  // for boundary closures, one of the following is possible...

  // if bc is dirichlet...
  for (int ifa = 0; ifa < nfa_b; ifa++) 
  {
    int icv     = cvofa[ifa][0];
    int icvFake = cvofa[ifa][1];
    for (int i = 0; i < 3; i++)
      grad_p[icv][i] += boundary_fa_lsg_coeff[ifa][i]*(p[icvFake]);//-p[icv]);
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
/*
  // if you know the normal component of the gradient...
  double (*coeff_b)[3][3] = new double[ncv_b][3][3];

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
  for (int icv = 0; icv < ncv_b; icv++) 
  {
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
  }
  */
  
  double epsMax = 1.0E-9;
  for (int icv = 0; icv < ncv; icv++) 
  {
    if ((fabs(grad_p[icv][0] - 1.12) > epsMax) || (fabs(grad_p[icv][1] - 2.13) > epsMax) || (fabs(grad_p[icv][2] - 3.14) > epsMax)) 
    {
      cerr << "Error: linear grad test failed: " << grad_p[icv][0] - 1.12 << " " << grad_p[icv][1] - 2.13 << " " << grad_p[icv][2] - 3.14 << endl;
      throw(-1);
    }
  }

  if (mpi_rank == 0)
    cout << "OK" << endl;
  
  delete[] p;
  delete[] grad_p;
}

int UgpWithCvFake::solveCvScalarJacobi(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double relax,
    const double zero, const int maxiter) {

  // need gradients in ghost as well...
  double (*grad_phi)[3] = new double[ncv_g][3];

  // res only in cvs we own...
  double * res = new double[ncv];

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // compute [Ap]{phi} + [Ap_grad]{grad_phi}...

    // start with {grad_phi}...
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      for (int i = 0; i < 3; i++)
        grad_phi[icv][i] = nbocv_lsg_coeff[noc_f][i] * phi[icv];
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
          grad_phi[icv][i] += nbocv_lsg_coeff[noc][i] * phi[icv_nbr];
      }
    }

    // also need the gradient in the ghosts...
    updateCvData(grad_phi, REPLACE_ROTATE_DATA);

    // {res} = [Ap]{phi} + [Ap_grad]{grad_phi}
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      res[icv] = Ap[noc_f] * phi[icv]; // diagonal
      for (int i = 0; i < 3; i++)
        res[icv] += Ap_grad[noc_f][i] * grad_phi[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        res[icv] += Ap[noc] * phi[icv_nbr];
        for (int i = 0; i < 3; i++)
          res[icv] += Ap_grad[noc][i] * grad_phi[icv_nbr][i];
      }
    }

    // update phi and compute max (Linf) normalized residual...
    double my_res_max = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
      double this_res = (rhs[icv] - res[icv]) / Ap[nbocv_i[icv]];
      my_res_max = max(my_res_max, fabs(this_res));
      // update phi with relaxation...
      phi[icv] += relax * this_res;
    }

    // update phi in ghost cv's of all processors...
    updateCvData(phi, REPLACE_DATA);

    // are we done? decide on one processor...
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank == 0) {
      //if (iter%10 == 0)
      //cout << "jacobi iter, res_max: " << iter << " " << res_max << endl;
      if (res_max < zero) {
        done = 1;
      } else if (iter > maxiter) {
        cout << "Warning: solveCvScalarJacobi did not converge after " << maxiter << " iters, res_max: " << res_max
            << endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

  } // while (done == 0)

  // cleanup...
  delete[] grad_phi;
  delete[] res;

  // let the calling routine know if we were successful...
  return (done == 1);

}

int UgpWithCvFake::solveCvScalarCg(double * phi, double * Ap, double * rhs, const int mode, const double zero, const int maxiter) {

  // we need the following work arrays...
  double * res = new double[ncv];
  double * v = new double[ncv];
  double * p = new double[ncv_g];

  // initialize...
  for (int icv = 0; icv < ncv; icv++)
    p[icv] = 0.0;

  double rho = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    res[icv] = Ap[noc_f] * phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f + 1; noc <= noc_l; noc++)
      res[icv] += Ap[noc] * phi[nbocv_v[noc]];
    res[icv] = rhs[icv] - res[icv];
  }

  double res0_max;
  if (mode == RELATIVE_RESIDUAL)
  {
    double my_res0_max = 0.0;

    for (int icv = 0; icv < ncv; icv++)
      my_res0_max = max(my_res0_max, fabs(res[icv] / Ap[nbocv_i[icv]]));

    MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      v[icv] = res[icv] / Ap[nbocv_i[icv]];

    double rho_prev = rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_rho += res[icv] * v[icv];

    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prev) < 1.0E-20)
      rho_prev = 1.0E-20;
    double beta = rho / rho_prev;

    for (int icv = 0; icv < ncv; icv++)
      p[icv] = v[icv] + beta * p[icv];
    updateCvData(p, REPLACE_DATA);

    // v = [Ap]{p}...
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      v[icv] = Ap[noc_f] * p[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        v[icv] += Ap[noc] * p[nbocv_v[noc]];
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_gamma += p[icv] * v[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    double alpha = rho / gamma;
    for (int icv = 0; icv < ncv_g; icv++)
      phi[icv] += alpha * p[icv];

    // check if we are done...
    if (iter % 5 == 0) {

      double my_res_max = 0.0;

      // recompute the residual...
      updateCvData(phi, REPLACE_DATA);
      for (int icv = 0; icv < ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;
        res[icv] = Ap[noc_f] * phi[icv]; // diagonal
        for (int noc = noc_f + 1; noc <= noc_l; noc++)
          res[icv] += Ap[noc] * phi[nbocv_v[noc]];
        // above is LHS. residual is then...
        res[icv] = rhs[icv] - res[icv];
        my_res_max = max(my_res_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) 
      {
        cout << "cg iter, res_max: " << iter << " " << res_max << endl;

        if ((mode == ABSOLUTE_RESIDUAL) && (res_max <= zero)) {
          done = 1;
        } else if ((mode == RELATIVE_RESIDUAL) && (res_max / (res0_max + 1.0E-20) <= zero)) {
          done = 1;
        } else if (iter > maxiter) {
          cout << "Warning: solveCvScalarCg did not converge after " << maxiter << " iters, res_max: " << res_max
              << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    } else {

      // on the other iterations, use this approximation...
      for (int icv = 0; icv < ncv; icv++)
        res[icv] -= alpha * v[icv];

    }

  }

  //delete[]
  delete[] res;
  delete[] v;
  delete[] p;

  // let the calling routine know if we were successful...
  return (done == 1);

}

int UgpWithCvFake::solveCvScalarBcgstab(double * phi, double * Ap, double * rhs,
    const int mode, const double zero, const int maxiter, char *scalarName) {

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
  for (int icv = 0; icv < ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  for (int icv = 0; icv < ncv_g; icv++) {
    phat[icv] = 0.0;
    shat[icv] = 0.0;
  }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; icv++) 
  {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    res[icv] = Ap[noc_f] * phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f + 1; noc <= noc_l; noc++)
      res[icv] += Ap[noc] * phi[nbocv_v[noc]];

    res[icv] = rhs[icv] - res[icv];
    res0[icv] = res[icv];
  }

  double res0_max;
  if (mode == RELATIVE_RESIDUAL) 
  {
    double my_res0_max = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_res0_max = max(my_res0_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
    MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  }

  int iter = 0;
  int done = 0;
  while (done == 0) 
  {
    iter++;

    double rho_prime = -omega * rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_rho += res[icv] * res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime) < 1.0E-20)
      rho_prime = 1.0E-20;
    double beta = alpha * rho / rho_prime;

    for (int icv = 0; icv < ncv; icv++)
      p[icv] = res[icv] - beta * (p[icv] - omega * v[icv]);

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      phat[icv] = p[icv] / Ap[nbocv_i[icv]];
    updateCvData(phat, REPLACE_DATA);

    // v = [A]{phat}
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      v[icv] = Ap[noc_f] * phat[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        v[icv] += Ap[noc] * phat[nbocv_v[noc]];
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_gamma += v[icv] * res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    alpha = rho / gamma;

    for (int icv = 0; icv < ncv; icv++)
      s[icv] = res[icv] - alpha * v[icv];

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      shat[icv] = s[icv] / Ap[nbocv_i[icv]];
    updateCvData(shat, REPLACE_DATA);

    // t = [A] shat...
    for (int icv = 0; icv < ncv; icv++) 
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      t[icv] = Ap[noc_f] * shat[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        t[icv] += Ap[noc] * shat[nbocv_v[noc]];
    }

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv < ncv; icv++) 
    {
      my_buf[0] += s[icv] * t[icv];
      my_buf[1] += t[icv] * t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0] / (buf[1] + 1.0E-20);

    // update phi...
    for (int icv = 0; icv < ncv; icv++)
      phi[icv] += alpha * phat[icv] + omega * shat[icv];
    updateCvData(phi, REPLACE_DATA);

    double my_res_max = 0.0;

    // recompute the residual...
    for (int icv = 0; icv < ncv; icv++) 
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      res[icv] = Ap[noc_f] * phi[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        res[icv] += Ap[noc] * phi[nbocv_v[noc]];
      // above is LHS. residual is then...
      res[icv] = rhs[icv] - res[icv];
      my_res_max = max(my_res_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
    }

    // check if we are done...
    if (iter % 2 == 0) 
    {
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

      if (mpi_rank == 0) 
      {
/*        if (iter % 40 == 0)
          cout << "bcgstab iter, res_max: " << iter << " " << res_max << endl;*/
        
        if ((mode == ABSOLUTE_RESIDUAL) && (res_max <= zero)) 
        {
//          cout << scalarName << " -> bcgstab iter, res_max: " << iter << " " << res_max << endl;
          done = 1;
        } 
        else if ((mode == RELATIVE_RESIDUAL) && (res_max / (res0_max + 1.0E-12) <= zero)) 
        {
          done = 1;
        } 
        else if (iter > maxiter) 
        {
          cout << "\nWarning: " << scalarName << " solveCvScalarBcgstab did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }

      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }
  }
  /*
  if (mpi_rank == 0)
    cout << "\n" << scalarName << " solveCvScalarBcgstab: " << iter << endl;*/

  //delete[]
  delete[] res;
  delete[] res0;
  delete[] p;
  delete[] v;
  delete[] t;
  delete[] s;
  delete[] phat;
  delete[] shat;

  // let the calling routine know if we were successful...
  return (done == 1);

}

//##############################################################################
//
//      math functions used for BCGSTAB,
//      ludeco ... LU decomposition
//      lusolv ... LU back substitution
//
//##############################################################################
void ludeco(double(*A)[5], int order) {
  for (int jc = 1; jc < order; jc++)
    A[0][jc] /= A[0][0];

  int jrjc = 0;

  for (;;) {
    jrjc++;
    int jrjcm1 = jrjc - 1;
    int jrjcp1 = jrjc + 1;
    for (int jr = jrjc; jr < order; jr++) {
      double sum = A[jr][jrjc];
      for (int jm = 0; jm <= jrjcm1; jm++)
        sum -= A[jr][jm] * A[jm][jrjc];
      A[jr][jrjc] = sum;
    }
    if (jrjc == (order - 1))
      return;

    for (int jc = jrjcp1; jc < order; jc++) {
      double sum = A[jrjc][jc];

      for (int jm = 0; jm <= jrjcm1; jm++)
        sum -= A[jrjc][jm] * A[jm][jc];

      A[jrjc][jc] = sum / A[jrjc][jrjc];
    }
  }
}

void lusolv(double(*A)[5], double *c, int order) {
  //              ...First l(inv)*b...
  c[0] = c[0] / A[0][0];
  for (int jr = 1; jr < order; jr++) {
    int jrm1 = jr - 1;
    double sum = c[jr];
    for (int jm = 0; jm <= jrm1; jm++)
      sum -= A[jr][jm] * c[jm];
    c[jr] = sum / A[jr][jr];
  }

  //             ...Next u(inv) of l(inv)*b...
  for (int jrjr = 1; jrjr < order; jrjr++) {
    int jr = (order - 1) - jrjr;
    int jrp1 = jr + 1;
    double sum = c[jr];
    for (int jmjm = jrp1; jmjm < order; jmjm++) {
      int jm = (order - 1) - jmjm + jrp1;
      sum -= A[jr][jm] * c[jm];
    }
    c[jr] = sum;
  }
}

void UgpWithCvFake::matTimesVecOverCVs(double(*res)[5], double(*Ap)[5][5], double(*phi)[5])
{
  for (int icv = 0; icv < ncv; icv++) 
  {
    for (int i = 0; i < 5; i++)
      res[icv][i] = 0.0;
    
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    for (int noc = noc_f; noc <= noc_l; noc++)      // cv diagonal + neighbors... 
    {
      const int icv_nbr = nbocv_v[noc];

      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          res[icv][i] += Ap[noc][i][j] * phi[icv_nbr][j];
    }
  }
}

void UgpWithCvFake::UpdateCvDataStateVec(double(*phi)[5]) {
  static double *scal = new double[ncv_g];
  static double (*vec)[3] = new double[ncv_g][3];

  for (int icv = 0; icv < ncv_g; icv++)
    scal[icv] = phi[icv][0];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv < ncv_g; icv++)
    phi[icv][0] = scal[icv];

  for (int icv = 0; icv < ncv_g; icv++)
    scal[icv] = phi[icv][4];
  updateCvData(scal, REPLACE_DATA);
  for (int icv = 0; icv < ncv_g; icv++)
    phi[icv][4] = scal[icv];

  for (int icv = 0; icv < ncv_g; icv++)
    for (int i = 0; i < 3; i++)
      vec[icv][i] = phi[icv][i + 1];
  updateCvData(vec, REPLACE_ROTATE_DATA);
  for (int icv = 0; icv < ncv_g; icv++)
    for (int i = 0; i < 3; i++)
      phi[icv][i + 1] = vec[icv][i];
}

int UgpWithCvFake::solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5],
                                      const double zeroAbs, const double zeroRel, const int maxiter, 
                                      char *scalarName, int check_interval) {

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
  for (int icv = 0; icv < ncv; icv++) {
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        LUDEC[icv][i][j] = Ap[nbocv_i[icv]][i][j];

    ludeco(LUDEC[icv], 5);
  }

  // initialize...
  for (int icv = 0; icv < ncv; icv++)
    for (int i = 0; i < 5; i++) {
      p[icv][i] = 0.0;
      v[icv][i] = 0.0;
    }

  for (int icv = 0; icv < ncv_g; icv++)
    for (int i = 0; i < 5; i++) {
      phat[icv][i] = 0.0;
      shat[icv][i] = 0.0;
    }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format ...
  matTimesVecOverCVs(res, Ap, phi);

  for (int icv = 0; icv < ncv; icv++)
  for (int i = 0; i < 5; i++) 
  {
    res[icv][i] = (rhs[icv][i] - res[icv][i]);
    res0[icv][i] = res[icv][i];
  }

  // compute first residual  
  double res0_max, my_res0_max = 0.0;
  
  for (int icv = 0; icv < ncv; icv++) 
  {
    double tmp[5];

    for (int i = 0; i < 5; i++) 
      tmp[i] = res[icv][i];
    
    lusolv(LUDEC[icv], tmp, 5);
    
    for (int i = 0; i < 5; i++) 
      my_res0_max = max(my_res0_max, fabs(tmp[i]));
  }

  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

  // start iteration

  int iter = 0;
  int done = 0;
  double res_max;

  while (done == 0) {

    iter++;

    double rho_prime = -omega*rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        my_rho += res[icv][i] * res0[icv][i];
    
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime) < 1.0E-20)
      rho_prime = 1.0E-20;

    double beta = alpha * rho / rho_prime;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        p[icv][i] = res[icv][i] - beta * (p[icv][i] - omega * v[icv][i]);

    // block-diagonal precon...
    // solve the system
    for (int icv = 0; icv < ncv; icv++) {
      for (int i = 0; i < 5; i++)
        phat[icv][i] = p[icv][i];
      lusolv(LUDEC[icv], phat[icv], 5);
    }

    UpdateCvDataStateVec(phat);

    // v = [A]{phat}
    matTimesVecOverCVs(v, Ap, phat);

    double my_gamma = 0.0;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        my_gamma += v[icv][i] * res0[icv][i];

    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

//    gamma = my_gamma; 
    
    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    alpha = rho/gamma;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        s[icv][i] = res[icv][i] - alpha * v[icv][i];

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++) {
      for (int i = 0; i < 5; i++)
        shat[icv][i] = s[icv][i];
      lusolv(LUDEC[icv], shat[icv], 5);
    }

    UpdateCvDataStateVec(shat);

    // t = [A] shat...
    matTimesVecOverCVs(t, Ap, shat);

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++) {
        my_buf[0] += s[icv][i] * t[icv][i];
        my_buf[1] += t[icv][i] * t[icv][i];
      }

    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0] / (buf[1] + 1.0E-20);

    // update phi...
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        phi[icv][i] += alpha * phat[icv][i] + omega * shat[icv][i];

    UpdateCvDataStateVec(phi);


    double my_res_max = 0.0;

    // recompute the residual...
    matTimesVecOverCVs(res, Ap, phi);

    for (int icv = 0; icv < ncv; icv++)           // above is LHS. residual is then...
    for (int i = 0; i < 5; i++)
      res[icv][i] = rhs[icv][i] - res[icv][i];

    // compute normalized res
    for (int icv = 0; icv < ncv; icv++) 
    {
      double tmp[5];

      for (int i = 0; i < 5; i++) 
        tmp[i] = res[icv][i];

      lusolv(LUDEC[icv], tmp, 5);

      for (int i = 0; i < 5; i++) 
        my_res_max = max(my_res_max, fabs(tmp[i]));
    }

    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

    // check residual 
    if (mpi_rank == 0)
    {
//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;
      
      if ((res_max <= zeroAbs) || (res_max / (res0_max + 1.0E-12) <= zeroRel))
        done = 1;
      else if (iter > maxiter)
      {
/*        cout << "Warning: " << scalarName 
             << " solveCvScalarBcgstab did not converge after " << maxiter 
             << " iters, res_max: " << res_max << endl;*/
        
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
  }

  if (check_interval == 0)
    if (mpi_rank == 0)
      cout << "\tbcgstab: iter: " << iter << " res_max: " << res_max << endl;

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
  return (done == 1);
}

int UgpWithCvFake::solveCvScalarBcgstab(double * phi, double * Ap, double(*Ap_grad)[3], double * rhs, const double zero, const int maxiter) {

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
  for (int icv = 0; icv < ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    for (int i = 0; i < 3; i++)
      grad[icv][i] = nbocv_lsg_coeff[noc_f][i] * phi[icv];
    for (int noc = noc_f + 1; noc <= noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      for (int i = 0; i < 3; i++)
        grad[icv][i] += nbocv_lsg_coeff[noc][i] * phi[icv_nbr];
    }
  }
  updateCvData(grad, REPLACE_ROTATE_DATA);
  for (int icv = 0; icv < ncv; icv++) {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    res[icv] = Ap[noc_f] * phi[icv]; // diagonal
    for (int i = 0; i < 3; i++)
      res[icv] += Ap_grad[noc_f][i] * grad[icv][i]; // diagonal gradient
    // cv neighbors...
    for (int noc = noc_f + 1; noc <= noc_l; noc++) {
      int icv_nbr = nbocv_v[noc];
      res[icv] += Ap[noc] * phi[icv_nbr];
      for (int i = 0; i < 3; i++)
        res[icv] += Ap_grad[noc][i] * grad[icv_nbr][i];
    }
    res[icv] = rhs[icv] - res[icv];
    res0[icv] = res[icv];
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    double rho_prime = -omega * rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_rho += res[icv] * res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime) < 1.0E-20)
      rho_prime = 1.0E-20;
    double beta = alpha * rho / rho_prime;

    for (int icv = 0; icv < ncv; icv++)
      p[icv] = res[icv] - beta * (p[icv] - omega * v[icv]);

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      phat[icv] = p[icv] / Ap[nbocv_i[icv]];
    updateCvData(phat, REPLACE_DATA);

    // v = [A]{phat}
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      for (int i = 0; i < 3; i++)
        grad[icv][i] = nbocv_lsg_coeff[noc_f][i] * phat[icv];
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
          grad[icv][i] += nbocv_lsg_coeff[noc][i] * phat[icv_nbr];
      }
    }
    updateCvData(grad, REPLACE_ROTATE_DATA);
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      v[icv] = Ap[noc_f] * phat[icv]; // diagonal
      for (int i = 0; i < 3; i++)
        v[icv] += Ap_grad[noc_f][i] * grad[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        v[icv] += Ap[noc] * phat[icv_nbr];
        for (int i = 0; i < 3; i++)
          v[icv] += Ap_grad[noc][i] * grad[icv_nbr][i];
      }
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_gamma += v[icv] * res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    alpha = rho / gamma;

    for (int icv = 0; icv < ncv; icv++)
      s[icv] = res[icv] - alpha * v[icv];

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      shat[icv] = s[icv] / Ap[nbocv_i[icv]];
    updateCvData(shat, REPLACE_DATA);

    // t = [A] shat...
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      for (int i = 0; i < 3; i++)
        grad[icv][i] = nbocv_lsg_coeff[noc_f][i] * shat[icv];
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        for (int i = 0; i < 3; i++)
          grad[icv][i] += nbocv_lsg_coeff[noc][i] * shat[icv_nbr];
      }
    }
    updateCvData(grad, REPLACE_ROTATE_DATA);
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      t[icv] = Ap[noc_f] * shat[icv]; // diagonal
      for (int i = 0; i < 3; i++)
        t[icv] += Ap_grad[noc_f][i] * grad[icv][i]; // diagonal gradient
      // cv neighbors...
      for (int noc = noc_f + 1; noc <= noc_l; noc++) {
        int icv_nbr = nbocv_v[noc];
        t[icv] += Ap[noc] * shat[icv_nbr];
        for (int i = 0; i < 3; i++)
          t[icv] += Ap_grad[noc][i] * grad[icv_nbr][i];
      }
    }

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
      my_buf[0] += s[icv] * t[icv];
      my_buf[1] += t[icv] * t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0] / (buf[1] + 1.0E-20);

    // update phi...
    for (int icv = 0; icv < ncv; icv++)
      phi[icv] += alpha * phat[icv] + omega * shat[icv];
    updateCvData(phi, REPLACE_DATA);

    // check if we are done...
    if (iter % 5 == 0) {

      double my_res_max = 0.0;

      // recompute the residual...
      for (int icv = 0; icv < ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;
        for (int i = 0; i < 3; i++)
          grad[icv][i] = nbocv_lsg_coeff[noc_f][i] * phi[icv];
        for (int noc = noc_f + 1; noc <= noc_l; noc++) {
          int icv_nbr = nbocv_v[noc];
          for (int i = 0; i < 3; i++)
            grad[icv][i] += nbocv_lsg_coeff[noc][i] * phi[icv_nbr];
        }
      }
      updateCvData(grad, REPLACE_ROTATE_DATA);
      for (int icv = 0; icv < ncv; icv++) {
        int noc_f = nbocv_i[icv];
        int noc_l = nbocv_i[icv + 1] - 1;
        res[icv] = Ap[noc_f] * phi[icv]; // diagonal
        for (int i = 0; i < 3; i++)
          res[icv] += Ap_grad[noc_f][i] * grad[icv][i]; // diagonal gradient
        // cv neighbors...
        for (int noc = noc_f + 1; noc <= noc_l; noc++) {
          int icv_nbr = nbocv_v[noc];
          res[icv] += Ap[noc] * phi[icv_nbr];
          for (int i = 0; i < 3; i++)
            res[icv] += Ap_grad[noc][i] * grad[icv_nbr][i];
        }
        res[icv] = rhs[icv] - res[icv];
        my_res_max = max(my_res_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        //cout << "bcgstab iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        } else if (iter > maxiter) {
          cout << "Warning: solveCvScalarBcgstab did not converge after " << maxiter << " iters, res_max: " << res_max
              << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    } else {

      // on the other iterations, use this approximation...
      for (int icv = 0; icv < ncv; icv++)
        res[icv] = s[icv] - omega * t[icv];

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
  return (done == 1);

}

