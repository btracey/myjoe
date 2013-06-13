#ifndef UGPWITHTOOLS_H
#define UGPWITHTOOLS_H

#include <fstream>

#include "Ugp.h"
#include "Param.h"
#include "Adt.h"
#include "Average.h"

// Tools0 contains point lookup functions...

class UgpWithTools0: public Ugp {

protected:

  double (*cvAdtBBMinOfRank)[3];
  double (*cvAdtBBMaxOfRank)[3];

public:

  Adt<double> * cvAdt;
  Average * average;

  UgpWithTools0() {

    if (mpi_rank==0) cout<<"UgpWithTools0()"<<endl;

    cvAdt = NULL;
    cvAdtBBMinOfRank = NULL;
    cvAdtBBMaxOfRank = NULL;

    // used for spatial averaging for certain grids
    // by dynamic models, statistics, ...

    average = NULL;

  }

protected:

  void init() {

    if (mpi_rank==0) cout<<"UgpWithTools0::init()"<<endl;

    Ugp::init();

  }

public:

  bool pointIsInsideCv(const double * xp, const int icv) {

    int foc, ifa;
    return (pointIsInsideCv(foc, ifa, xp, icv));

  }

  bool pointIsInsideCv(int& foc, int& ifa, const double * xp, const int icv) {

    // this routine determines if a point is inside a cv...

    // get the approx cv center...

    const int noc_f = noocv_i[icv];
    const int noc_l = noocv_i[icv+1]-1;
    double x_cv[3] = { 0.0, 0.0, 0.0 };
    for (int noc = noc_f; noc<=noc_l; ++noc) {
      const int ino = noocv_v[noc];
      FOR_I3
        x_cv[i] += x_no[ino][i];
    }
    double tmp = 1.0/(double) (noc_l-noc_f+1);
    FOR_I3
      x_cv[i] *= tmp;

    // the point relative to the cv center...

    double dxp[3];
    FOR_I3
      dxp[i] = xp[i]-x_cv[i];

    // cycle through the faces...

    const int foc_f = faocv_i[icv];
    const int foc_l = faocv_i[icv+1]-1;
    for (foc = foc_f; foc<=foc_l; ++foc) {
      ifa = faocv_v[foc];
      // get the approx fa center...
      const int nof_f = noofa_i[ifa];
      const int nof_l = noofa_i[ifa+1]-1;
      double dx_fa[3] = { 0.0, 0.0, 0.0 };
      for (int nof = nof_f; nof<=nof_l; ++nof) {
        const int ino = noofa_v[nof];
        FOR_I3
          dx_fa[i] += x_no[ino][i];
      }
      tmp = 1.0/(double) (nof_l-nof_f+1);
      FOR_I3
        dx_fa[i] = dx_fa[i]*tmp-x_cv[i];

      // for this face to contain the point xp, the dot-product of dx_fa and dxp
      // should be positive...

      //if ( dxp[0]*dx_fa[0] + dxp[1]*dx_fa[1] + dxp[2]*dx_fa[2] >= 0.0 ) {
      if (1) {

        // the face normal orientation (nodes follow rh-rule)...

        double fa_direction;
        if (cvofa[ifa][0]==icv) fa_direction = 1.0;
        else fa_direction = -1.0;

        // now loop on face sub-tri's...

        int ino1 = noofa_v[nof_l];
        double dx_no1[3];
        FOR_I3
          dx_no1[i] = x_no[ino1][i]-x_cv[i];
        double cv_fa_no1_xp_vol = fa_direction*calcSignedTetVolume(dx_fa, dx_no1, dxp);

        for (int nof = nof_f; nof<=nof_l; ++nof) {

          // copy the ino1 data into ino0...
          const int ino0 = ino1;
          double dx_no0[3];
          FOR_I3
            dx_no0[i] = dx_no1[i];
          const double cv_fa_no0_xp_vol = cv_fa_no1_xp_vol;

          // next ino1 data...
          ino1 = noofa_v[nof];
          FOR_I3
            dx_no1[i] = x_no[ino1][i]-x_cv[i];
          cv_fa_no1_xp_vol = fa_direction*calcSignedTetVolume(dx_fa, dx_no1, dxp);

          // if cv_fa_no0_xp is positive and cv_fa_no1_xp is negative, then
          // this sub-tet may contain the point...
          if ((cv_fa_no0_xp_vol>=0.0)&&(cv_fa_no1_xp_vol<=0.0)) {

            const double cv_no0_no1_xp_vol = fa_direction*calcSignedTetVolume(dx_no0, dx_no1, dxp);

            if (cv_no0_no1_xp_vol>=0.0) {

              // this subface must contain the point, although it may be in the next cell...

              const double cv_no0_no1_fa_vol = fa_direction*calcSignedTetVolume(dx_no0, dx_no1, dx_fa);

              if (cv_fa_no0_xp_vol-cv_fa_no1_xp_vol+cv_no0_no1_xp_vol<=cv_no0_no1_fa_vol) return true;
              else return false;

            }
          }
        }
      }
    }

    cerr<<"Error: did not find a subtri that contains the point: "<<xp[0]<<" "<<xp[1]<<" "<<xp[2]<<endl;
    throw(-1);

  }

  inline double calcSignedTetVolume(const double * dx1, const double * dx2, const double * dx3) const {
    return (dx3[0]*(dx1[1]*dx2[2]-dx1[2]*dx2[1])+dx3[1]*(dx1[2]*dx2[0]-dx1[0]*dx2[2])+dx3[2]*(dx1[0]*dx2[1]-dx1[1]
        *dx2[0]));
  }

  void useCvAdt() {

    // in the future we should rebuild this here
    // if the mesh has been adapted. For now, just consider
    // whether it has been initialied once...

    if (cvAdt==NULL) initCvAdt();

  }

  int setCvInterpWeightsForPoint(int &cv_size, int * cv_list, double * cv_weight, const double * xp) {

    // for now, this returns one point with weight 1 if a containing cv 
    // is found...
    useCvAdt();

    int nList;
    int cvList[ADT_LIST_MAX];
    cvAdt->buildListForPoint(nList, cvList, xp);

    // now cycle through each cv and check...
    for (int ii = 0; ii<nList; ii++) {

      const int icv = cvList[ii];

      // look for a subtet in this cv that contains xp...

      const int noc_f = noocv_i[icv];
      const int noc_l = noocv_i[icv+1]-1;
      const int nnoc = noc_l-noc_f+1;
      double x_cv[3] = { 0.0, 0.0, 0.0 };
      for (int noc = noc_f; noc<=noc_l; ++noc) {
        const int ino = noocv_v[noc];
        for (int i = 0; i<3; ++i)
          x_cv[i] += x_no[ino][i];
      }
      double tmp = 1.0/(double) nnoc;
      for (int i = 0; i<3; ++i)
        x_cv[i] *= tmp;

      // now loop on faces...

      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; ++foc) {
        const int ifa = faocv_v[foc];
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
        }
        double dx_fa[3] = { 0.0, 0.0, 0.0 };
        for (int nof = nof_begin; nof!=nof_end; nof += nof_inc) {
          const int ino = noofa_v[nof];
          for (int i = 0; i<3; ++i)
            dx_fa[i] += x_no[ino][i];
        }
        tmp = 1.0/(double) nnof;
        for (int i = 0; i<3; ++i)
          dx_fa[i] = dx_fa[i]*tmp-x_cv[i];

        // now loop on subtets...

        int ino1 = noofa_v[nof_end-nof_inc];
        double dx_no1[3] = { x_no[ino1][0]-x_cv[0], x_no[ino1][1]-x_cv[1], x_no[ino1][2]-x_cv[2] };

        for (int nof = nof_begin; nof!=nof_end; nof += nof_inc) {

          // copy down previous ino1 stuff to ino0...

          const int ino0 = ino1;
          const double dx_no0[3] = { dx_no1[0], dx_no1[1], dx_no1[2] };

          // and take the next one...

          ino1 = noofa_v[nof];
          dx_no1[0] = x_no[ino1][0]-x_cv[0];
          dx_no1[1] = x_no[ino1][1]-x_cv[1];
          dx_no1[2] = x_no[ino1][2]-x_cv[2];

          // now this subtet is well defined. We will add a phi gradient to
          // ino0 and/or ino1 if the reversed velocity vector is contained in
          // the tet...

          // the full volume is...
          const double vol = dx_no1[0]*(dx_fa[1]*dx_no0[2]-dx_fa[2]*dx_no0[1])+dx_no1[1]*(dx_fa[2]*dx_no0[0]-dx_fa[0]
              *dx_no0[2])+dx_no1[2]*(dx_fa[0]*dx_no0[1]-dx_fa[1]*dx_no0[0]);
          assert(vol>0.0);

          // ==========================
          // get the location of the point relative to ino0...
          // ==========================

          const double dx_p0[3] = { xp[0]-x_no[ino0][0], xp[1]-x_no[ino0][1], xp[2]-x_no[ino0][2] };

          const double vol_ino1 = dx_p0[0]*(dx_fa[1]*dx_no0[2]-dx_fa[2]*dx_no0[1])+dx_p0[1]*(dx_fa[2]*dx_no0[0]
              -dx_fa[0]*dx_no0[2])+dx_p0[2]*(dx_fa[0]*dx_no0[1]-dx_fa[1]*dx_no0[0]);

          if (vol_ino1>=0.0) {

            const double vol_ifa = dx_p0[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1])+dx_p0[1]*(dx_no0[2]*dx_no1[0]
                -dx_no0[0]*dx_no1[2])+dx_p0[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);

            if (vol_ifa>=0.0) {

              const double dx_no0_[3] = { dx_no0[0]-dx_fa[0], dx_no0[1]-dx_fa[1], dx_no0[2]-dx_fa[2] };

              const double dx_no1_[3] = { dx_no1[0]-dx_fa[0], dx_no1[1]-dx_fa[1], dx_no1[2]-dx_fa[2] };

              const double vol_icv = dx_p0[0]*(dx_no1_[1]*dx_no0_[2]-dx_no1_[2]*dx_no0_[1])+dx_p0[1]*(dx_no1_[2]
                  *dx_no0_[0]-dx_no1_[0]*dx_no0_[2])+dx_p0[2]*(dx_no1_[0]*dx_no0_[1]-dx_no1_[1]*dx_no0_[0]);

              if (vol_icv>=0.0) {

                // we got 3 positive volumes...
                if (vol_ino1+vol_ifa+vol_icv<=vol) {

                  // this is it...
                  cv_size = 1;
                  cv_list[0] = icv;
                  cv_weight[0] = 1.0;
                  return (1);

                }
              }
            }
          }
        }
      }
    }

    return (0);

  }

  int setNodeInterpWeightsForPoint(int &no_size, int * no_list, double * no_weight, const double * xp) {

    // for this one we do not know the cv, so we need to use the cvAdt..
    useCvAdt();

    // XXXXX some day this should return a vector<int>...
    int nList;
    int cvList[ADT_LIST_MAX];
    cvAdt->buildListForPoint(nList, cvList, xp);

    // now cycle through each cv and compute the weights...
    for (int i = 0; i<nList; i++) {

      const int icv = cvList[i];

      // compute the cv center - start by zeroing
      // no_flag in all this cv's nodes...
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++) {
          const int ino = noofa_v[nof];
          no_flag[ino] = -1;
        }
      }
      // and now cycle through the nodes, visiting each
      // cv node once...
      int nnoc = 0;
      double x_cv[3] = { 0.0, 0.0, 0.0 };
      for (int foc = foc_f; foc<=foc_l; foc++) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++) {
          const int ino = noofa_v[nof];
          if (no_flag[ino]==-1) {
            no_flag[ino] = nnoc++;
            for (int i = 0; i<3; i++)
              x_cv[i] += x_no[ino][i];
          }
        }
      }
      for (int i = 0; i<3; i++)
        x_cv[i] /= (double) nnoc;

      // so the point location relative to this x_cv...

      double dx_p[3];
      for (int i = 0; i<3; i++)
        dx_p[i] = xp[i]-x_cv[i];

      // now find a face where all cross products are positive...

      for (int foc = foc_f; foc<=foc_l; foc++) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l-nof_f+1;
        double dx_fa[3] = { 0.0, 0.0, 0.0 };
        for (int nof = nof_f; nof<=nof_l; nof++) {
          const int ino = noofa_v[nof];
          for (int i = 0; i<3; i++)
            dx_fa[i] += x_no[ino][i];
        }
        for (int i = 0; i<3; i++)
          dx_fa[i] = dx_fa[i]/(double) nnof-x_cv[i];
        // and we need an fa sign...
        int fa_sign;
        if (cvofa[ifa][0]==icv) {
          fa_sign = 1;
        }
        else {
          assert(cvofa[ifa][1]==icv);
          fa_sign = -1;
        }
        // loop on edges...
        int ino1 = noofa_v[nof_l];
        double dx_no1[3];
        for (int i = 0; i<3; i++)
          dx_no1[i] = x_no[ino1][i]-x_cv[i];
        double vol_no0 = (double) fa_sign*(dx_p[0]*(dx_no1[1]*dx_fa[2]-dx_no1[2]*dx_fa[1])+dx_p[1]*(dx_no1[2]*dx_fa[0]
            -dx_no1[0]*dx_fa[2])+dx_p[2]*(dx_no1[0]*dx_fa[1]-dx_no1[1]*dx_fa[0]));
        for (int nof = nof_f; nof<=nof_l; nof++) {
          const int ino0 = ino1;
          const double vol_no1 = -vol_no0;
          ino1 = noofa_v[nof];
          double dx_no0[3];
          for (int i = 0; i<3; i++) {
            dx_no0[i] = dx_no1[i];
            dx_no1[i] = x_no[ino1][i]-x_cv[i];
          }
          vol_no0 = (double) fa_sign*(dx_p[0]*(dx_no1[1]*dx_fa[2]-dx_no1[2]*dx_fa[1])+dx_p[1]*(dx_no1[2]*dx_fa[0]
              -dx_no1[0]*dx_fa[2])+dx_p[2]*(dx_no1[0]*dx_fa[1]-dx_no1[1]*dx_fa[0]));
          if ((vol_no0>=0.0)&&(vol_no1>=0.0)) {
            double tmp[3];
            tmp[0] = dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1];
            tmp[1] = dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2];
            tmp[2] = dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0];
            double vol_fa = (double) fa_sign*(dx_p[0]*tmp[0]+dx_p[1]*tmp[1]+dx_p[2]*tmp[2]);
            if (vol_fa>=0.0) {
              // check the subtet volume and compare...
              double vol_st = (double) fa_sign*(dx_fa[0]*tmp[0]+dx_fa[1]*tmp[1]+dx_fa[2]*tmp[2]);
              assert(vol_st>0.0);
              double vol_cv = vol_st-vol_no0-vol_no1-vol_fa;
              if (vol_cv>-1.0E-10*vol_st) {
                no_size = nnoc;
                for (int i = 0; i<no_size; i++)
                  no_weight[i] = 0.0;
                no_weight[no_flag[ino0]] += vol_no0/vol_st;
                no_weight[no_flag[ino1]] += vol_no1/vol_st;
                for (int nof2 = nof_f; nof2<=nof_l; nof2++) {
                  const int ino2 = noofa_v[nof2];
                  no_weight[no_flag[ino2]] += vol_fa/(double) nnof/vol_st;
                }
                // XXXXX noocv_i/v would really help this routine...
                for (int foc2 = foc_f; foc2<=foc_l; foc2++) {
                  const int ifa2 = faocv_v[foc2];
                  const int nof2_f = noofa_i[ifa2];
                  const int nof2_l = noofa_i[ifa2+1]-1;
                  for (int nof2 = nof2_f; nof2<=nof2_l; nof2++) {
                    const int ino2 = noofa_v[nof2];
                    if (no_flag[ino2]>=0) {
                      no_list[no_flag[ino2]] = ino2;
                      no_weight[no_flag[ino2]] += vol_cv/(double) nnoc/vol_st;
                      no_flag[ino2] = -1;
                    }
                  }
                }
                return (1);
              }
            }
          }
        }
      }
    }
    return (0);
  }

protected:

  void initCvAdt() {

    assert(cvAdt==NULL);

    assert(faocv_i!=NULL);
    assert(faocv_v!=NULL);
    assert(noofa_i!=NULL);
    assert(noofa_v!=NULL);

    double (*bbmin)[3] = new double[ncv][3];
    double (*bbmax)[3] = new double[ncv][3];

    for (int icv = 0; icv<ncv; icv++) {
      double my_bbmin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
      double my_bbmax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++) {
          const int ino = noofa_v[nof];
          for (int i = 0; i<3; i++) {
            my_bbmin[i] = min(x_no[ino][i], my_bbmin[i]);
            my_bbmax[i] = max(x_no[ino][i], my_bbmax[i]);
          }
        }
      }
      // add a small fraction to the bounding box to ensure it
      // contains the cv...
      for (int i = 0; i<3; i++) {
        const double delta = 1.0E-6*(my_bbmax[i]-my_bbmin[i]);
        bbmin[icv][i] = my_bbmin[i]-delta;
        bbmax[icv][i] = my_bbmax[i]+delta;
      }
    }

    cvAdt = new Adt<double> (ncv, bbmin, bbmax);

    // also, everyone gets a copy of everyone's bounding box. This helps
    // when we need to decide where to send data...

    assert(cvAdtBBMinOfRank==NULL);
    assert(cvAdtBBMaxOfRank==NULL);
    cvAdtBBMinOfRank = new double[mpi_size][3];
    cvAdtBBMaxOfRank = new double[mpi_size][3];

    double cvAdtBBMin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
    double cvAdtBBMax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
    for (int icv = 0; icv<ncv; icv++) {
      for (int i = 0; i<3; i++) {
        cvAdtBBMin[i] = min(cvAdtBBMin[i], bbmin[icv][i]);
        cvAdtBBMax[i] = max(cvAdtBBMax[i], bbmax[icv][i]);
      }
    }
    MPI_Allgather(cvAdtBBMin, 3, MPI_DOUBLE, cvAdtBBMinOfRank, 3, MPI_DOUBLE, mpi_comm);
    MPI_Allgather(cvAdtBBMax, 3, MPI_DOUBLE, cvAdtBBMaxOfRank, 3, MPI_DOUBLE, mpi_comm);

    delete[] bbmin;
    delete[] bbmax;

  }

};

// =====================================================
// basic statistics classes for scalars and vectors...
// =====================================================

class ScalarStats {

public:

  DoubleScalar * ds;
  double *avg;
  double *rms;
  double weight;

  ScalarStats(DoubleScalar * ds) {

    //if (mpi_rank == 0)
    //cout << "adding stats for scalar: " << ds->getName() << endl;

    // take a copy of the data pointer...
    this->ds = ds;
    avg = NULL;
    rms = NULL;

    // XXXXX to restart with stats, this HAS to be a registered value...
    weight = 0.0;

  }

  virtual void update(const int n, const double new_weight, Average * average) {

    double * d = *(ds->ptr);
    double old_weight = weight;
    weight += new_weight;
    for (int i = 0; i<n; i++) {
      // put sum( phi**2 * new_weight )_new / sum( new_weight ) in rms[i]...
      rms[i] = (old_weight*(rms[i]*rms[i]+avg[i]*avg[i])+new_weight*(d[i]*d[i]))/weight;
      // put new means in r1_avg...
      avg[i] = (old_weight*avg[i]+new_weight*d[i])/weight;
    }

    if (average!=NULL) {
      assert(average->size()==n); // check consistency - only CV_DATA for now
      average->apply(avg);
      average->apply(rms);
    }

    for (int i = 0; i<n; i++) {
      // complete rms...
      rms[i] = sqrt(fabs(rms[i]-avg[i]*avg[i]));
    }

  }

  // quick fix
  void updateUQ(const int n, const double new_weight, Average * average) {

    double * d = *(ds->ptr);
    double old_weight = weight;
    weight += new_weight;
    for (int i = 0; i<n; i++)
    {
      rms[i] *= rms[i];
      rms[i] += avg[i]*avg[i];
      avg[i] += d[i]*new_weight;
      rms[i] += d[i]*d[i]*new_weight;
      rms[i] -= avg[i]*avg[i];
      rms[i] = sqrt(fabs(rms[i]));
    }

    if (average!=NULL) {
      assert(average->size()==n); // check consistency - only CV_DATA for now
      average->apply(avg);
      average->apply(rms);
    }
  }

  void reset(const int n) {

    weight = 0.0;
    for (int i = 0; i<n; i++)
      avg[i] = rms[i] = 0.0;

  }

};

class VectorStats {

public:

  DoubleVector * dv;
  double (*avg)[3];
  double (*rms)[3];
  double (*rey)[3];
  double weight;

  VectorStats(DoubleVector * dv) {

    //if (mpi_rank == 0)
    //cout << "adding stats for vector: " << dv->getName() << endl;

    // take a copy of the data pointer...
    this->dv = dv;
    avg = NULL;
    rms = NULL;
    rey = NULL;

    // XXXXX to restart with stats, this HAS to be a registered value...
    weight = 0.0;

  }

  void update(const int n, const double new_weight, Average * average) {

    double (*d)[3] = *(dv->ptr);
    double old_weight = weight;
    weight += new_weight;
    for (int i = 0; i<n; i++) {

      // put sum( ui**2 * new_weight )_new / sum( new_weight ) in rms(i)...
      for (int j = 0; j<3; j++)
        rms[i][j] = (old_weight*(rms[i][j]*rms[i][j]+avg[i][j]*avg[i][j])+new_weight*(d[i][j]*d[i][j]))/weight;
      // put sum( ui*uj * new_weight )_new / sum( new_weight ) in rey(i)...
      rey[i][0] = (old_weight*(rey[i][0]+avg[i][1]*avg[i][2])+new_weight*(d[i][1]*d[i][2]))/weight;
      rey[i][1] = (old_weight*(rey[i][1]+avg[i][2]*avg[i][0])+new_weight*(d[i][2]*d[i][0]))/weight;
      rey[i][2] = (old_weight*(rey[i][2]+avg[i][0]*avg[i][1])+new_weight*(d[i][0]*d[i][1]))/weight;
      // put new means in avg...
      for (int j = 0; j<3; j++)
        avg[i][j] = (old_weight*avg[i][j]+new_weight*d[i][j])/weight;

    }

    if (average!=NULL) {
      assert(average->size()==n); // check consistency - only CV_DATA support right now
      average->apply(avg);
      average->apply(rms);
      average->apply(rey);
    }

    for (int i = 0; i<n; i++) {

      // complete rms...
      for (int j = 0; j<3; j++)
        rms[i][j] = sqrt(fabs(rms[i][j]-avg[i][j]*avg[i][j]));
      // complete rey...
      rey[i][0] -= avg[i][1]*avg[i][2];
      rey[i][1] -= avg[i][2]*avg[i][0];
      rey[i][2] -= avg[i][0]*avg[i][1];

    }

  }

  // quick fix
  void updateUQ(const int n, const double new_weight, Average * average) {

    double (*d)[3] = *(dv->ptr);
    double old_weight = weight;
    weight += new_weight;
    for (int i = 0; i<n; i++)
      for (int j = 0; j<3; j++)
      {
        rms[i][j] *= rms[i][j];
        rms[i][j] += avg[i][j]*avg[i][j];
        avg[i][j] += d[i][j]*new_weight;
        rms[i][j] += d[i][j]*d[i][j]*new_weight;
        rms[i][j] -= avg[i][j]*avg[i][j];
        rms[i][j] = sqrt(fabs(rms[i][j]));
      }

    if (average!=NULL) {
      assert(average->size()==n); // check consistency - only CV_DATA support right now
      average->apply(avg);
      average->apply(rms);
      average->apply(rey);
    }
  }

  void reset(const int n) {

    weight = 0.0;
    for (int i = 0; i<n; i++)
      for (int j = 0; j<3; j++)
        avg[i][j] = rms[i][j] = rey[i][j] = 0.0;

  }

};

class UgpWithTools1: public UgpWithTools0 {

protected:

private:

  // statistics lists...
  list<ScalarStats> scalarStatsList;
  list<VectorStats> vectorStatsList;
  int stats_average_kind;

public:

  UgpWithTools1() {

    if (mpi_rank==0) cout<<"UgpWithTools1()"<<endl;

    // assume no special averaging of the stats: i.e. stats will remain
    // point-wise. You can always add averaging to postpro when you 
    // look at the restart file...
    stats_average_kind = AVERAGE_NONE;

  }

protected:

  void init() {

    if (mpi_rank==0) cout<<"UgpWithTools1::init()"<<endl;

    UgpWithTools0::init();

  }

  /**

   Call initStats() in a class that inherits UgpWithTools to make
   statistics capabilities available. It should be called BEFORE the
   restart file is read, so that the statistics themselves can be
   stored and restarted if desired, e.g.:

   \verbatim
   Param * p;
   if (getParam(p,"STATS"))
   initStats(p);
   \endverbatim

   where the parameter line containing STATS should simply be a list of
   registered data names, e.g.:

   \verbatim
   STATS U P PHI1 PHI2
   \endverbatim

   This initialization routine will register new data to support the
   statistics. For scalar PHI1 the registered data will be:

   \verbatim
   PHI1_WGT
   PHI1_AVG
   PHI1_RMS
   \endverbatim

   and for a vector U, the registered data will be:

   \verbatim
   U_WGT
   U_AVG
   U_RMS
   U_REY
   \endverbatim

   */

  //void initStats(Param * p);

  void initStats(ParamMap * params);

  void initStatsUQ(ParamMap * params);

  void updateStats(const double weight, const int verbose = 0) {

    Average * stats_average = NULL;

    switch (stats_average_kind) {
    case AVERAGE_NONE:
      break;
    case AVERAGE_X_Z:
      if (average==NULL) {
        if (mpi_rank==0) cout<<" > setting up average for STATS: X_Z"<<endl;
        double * y_cv = new double[ncv];
        double * vol_cv = new double[ncv];
        FOR_ICV {
          double x[3];
          calcCvVolumeAndCentroid(vol_cv[icv], x, icv);
          y_cv[icv] = x[1];
        }
        average = new Average(y_cv, vol_cv, ncv, 1.0E-7);
        average->setKind(AVERAGE_X_Z);
        delete[] y_cv;
        delete[] vol_cv;
      }
      else {
        // is some other part of the code has set up this averaging operator,
        // make sure it is as expected...
        assert(average->getKind()==AVERAGE_X_Z);
      }
      // set the stats_average to the average...
      stats_average = average;
      break;
    case AVERAGE_THETAX:
      if (average==NULL) {
        if (mpi_rank==0) cout<<" > setting up average for STATS: THETAX"<<endl;
        double * x_cv = new double[ncv];
        double * r_cv = new double[ncv];
        double * vol_cv = new double[ncv];
        FOR_ICV {
          double x[3];
          calcCvVolumeAndCentroid(vol_cv[icv], x, icv);
          x_cv[icv] = x[0];
          r_cv[icv] = sqrt(x[1]*x[1]+x[2]*x[2]);
        }
        average = new Average(x_cv, r_cv, vol_cv, ncv, 1.0E-7);
        average->setKind(AVERAGE_THETAX);
        delete[] x_cv;
        delete[] r_cv;
        delete[] vol_cv;
      }
      else {
        // is some other part of the code has set up this averaging operator,
        // make sure it is as expected...
        assert(average->getKind()==AVERAGE_THETAX);
      }
      // set the stats_average to the average...
      stats_average = average;
      break;
    default:
      if (mpi_rank==0) cout<<"Warning: unrecognized stats_average_kind. Skipping."<<endl;
    }

    // scalars first...

    for (list<ScalarStats>::iterator stats = scalarStatsList.begin(); stats!=scalarStatsList.end(); stats++) {

      if ((mpi_rank==0)&&(verbose==1)) cout<<" > updating stats for scalar: "<<stats->ds->getName()
          <<", new averaging time: ";

      switch (stats->ds->getDatatype()) {
      case NO_DATA:
        stats->update(nno, weight, stats_average);
        break;
      case CV_DATA:
        stats->update(ncv, weight, stats_average);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: updateStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }

      if ((mpi_rank==0)&&(verbose==1)) cout<<stats->weight<<endl;

    }

    // then vectors...

    for (list<VectorStats>::iterator stats = vectorStatsList.begin(); stats!=vectorStatsList.end(); stats++) {

      if ((mpi_rank==0)&&(verbose==1)) cout<<" > updating stats for vector: "<<stats->dv->getName()
          <<", new averaging time: ";

      switch (stats->dv->getDatatype()) {
      case NO_DATA:
        stats->update(nno, weight, stats_average);
        break;
      case CV_DATA:
        stats->update(ncv, weight, stats_average);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: updateStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }

      if ((mpi_rank==0)&&(verbose==1)) cout<<stats->weight<<endl;

    }

  }


  // update stats for uncertainty quantification
  // quick fix is copy code. otherwise big changes would be necessary: stats should have operators or templates.

  void updateStatsUQ(const double weight, const int verbose = 0) {

    Average * stats_average = NULL;

    switch (stats_average_kind) {
    case AVERAGE_NONE:
      break;
    case AVERAGE_X_Z:
      if (average==NULL) {
        if (mpi_rank==0) cout<<" > setting up average for STATS: X_Z"<<endl;
        double * y_cv = new double[ncv];
        double * vol_cv = new double[ncv];
        FOR_ICV {
          double x[3];
          calcCvVolumeAndCentroid(vol_cv[icv], x, icv);
          y_cv[icv] = x[1];
        }
        average = new Average(y_cv, vol_cv, ncv, 1.0E-7);
        average->setKind(AVERAGE_X_Z);
        delete[] y_cv;
        delete[] vol_cv;
      }
      else {
        // is some other part of the code has set up this averaging operator,
        // make sure it is as expected...
        assert(average->getKind()==AVERAGE_X_Z);
      }
      // set the stats_average to the average...
      stats_average = average;
      break;
    case AVERAGE_THETAX:
      if (average==NULL) {
        if (mpi_rank==0) cout<<" > setting up average for STATS: THETAX"<<endl;
        double * x_cv = new double[ncv];
        double * r_cv = new double[ncv];
        double * vol_cv = new double[ncv];
        FOR_ICV {
          double x[3];
          calcCvVolumeAndCentroid(vol_cv[icv], x, icv);
          x_cv[icv] = x[0];
          r_cv[icv] = sqrt(x[1]*x[1]+x[2]*x[2]);
        }
        average = new Average(x_cv, r_cv, vol_cv, ncv, 1.0E-7);
        average->setKind(AVERAGE_THETAX);
        delete[] x_cv;
        delete[] r_cv;
        delete[] vol_cv;
      }
      else {
        // is some other part of the code has set up this averaging operator,
        // make sure it is as expected...
        assert(average->getKind()==AVERAGE_THETAX);
      }
      // set the stats_average to the average...
      stats_average = average;
      break;
    default:
      if (mpi_rank==0) cout<<"Warning: unrecognized stats_average_kind. Skipping."<<endl;
    }

    // scalars first...

    for (list<ScalarStats>::iterator stats = scalarStatsList.begin(); stats!=scalarStatsList.end(); stats++) {

      if ((mpi_rank==0)&&(verbose==1)) cout<<" > updating stats for scalar: "<<stats->ds->getName()
          <<", new averaging time: ";

      switch (stats->ds->getDatatype()) {
      case NO_DATA:
        stats->updateUQ(nno, weight, stats_average);
        break;
      case CV_DATA:
        stats->updateUQ(ncv, weight, stats_average);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: updateStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }

      if ((mpi_rank==0)&&(verbose==1)) cout<<stats->weight<<endl;

    }

    // then vectors...

    for (list<VectorStats>::iterator stats = vectorStatsList.begin(); stats!=vectorStatsList.end(); stats++) {

      if ((mpi_rank==0)&&(verbose==1)) cout<<" > updating stats for vector: "<<stats->dv->getName()
          <<", new averaging time: ";

      switch (stats->dv->getDatatype()) {
      case NO_DATA:
        stats->updateUQ(nno, weight, stats_average);
        break;
      case CV_DATA:
        stats->updateUQ(ncv, weight, stats_average);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: updateStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }

      if ((mpi_rank==0)&&(verbose==1)) cout<<stats->weight<<endl;

    }

  }

  /**

   Resets the statistics - call this at the end of the initialization
   phase when all data is fully allocated. e.g.:

   \verbatim
   if (checkParam("RESET_STATS"))
   resetStats();
   \endverbatim

   */

  void resetStats() {

    if (mpi_rank==0) cout<<"resetStats()"<<endl;

    // scalars first...

    for (list<ScalarStats>::iterator stats = scalarStatsList.begin(); stats!=scalarStatsList.end(); stats++) {
      switch (stats->ds->getDatatype()) {
      case NO_DATA:
        stats->reset(nno);
        break;
      case CV_DATA:
        stats->reset(ncv);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: resetStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }
    }

    // then vectors...

    for (list<VectorStats>::iterator stats = vectorStatsList.begin(); stats!=vectorStatsList.end(); stats++) {
      switch (stats->dv->getDatatype()) {
      case NO_DATA:
        stats->reset(nno);
        break;
      case CV_DATA:
        stats->reset(ncv);
        break;
      default:
        if (mpi_rank==0) cerr<<"Error: resetStats: only nodal and cv data supported in stats for now."<<endl;
        throw(-1);
      }
    }
  }

public:

  void setStatsFlags() {

    // scalars first...

    for (list<ScalarStats>::iterator stats = scalarStatsList.begin(); stats!=scalarStatsList.end(); stats++) {
      setScalarFlag(stats->avg);
      setScalarFlag(stats->rms);
    }

    // then vectors...

    for (list<VectorStats>::iterator stats = vectorStatsList.begin(); stats!=vectorStatsList.end(); stats++) {
      setVectorFlag(stats->avg);
      setVectorFlag(stats->rms);
      setVectorFlag(stats->rey);
    }

  }

};

// ==========================================================
// these classes implement the virtual routine "flagCvs"
// to set the cv flag to 1 in a user-specified part of the
// geometry.
// ==========================================================

class WriteDataGeom {

private:

public:

  virtual ~WriteDataGeom() {
    ;
  }

  virtual void flagCvs(Ugp * ugp) = 0;

};

class AllWriteDataGeom: public WriteDataGeom {

public:

  virtual ~AllWriteDataGeom() {
  }

  void flagCvs(Ugp * ugp) {
    // All just flags everything...
    for (int icv = 0; icv<ugp->getNcv(); icv++)
      ugp->cv_flag[icv] = 1;
  }

};

class PlaneWriteDataGeom: public WriteDataGeom {

private:

  double x, y, z, nx, ny, nz;

public:

  virtual ~PlaneWriteDataGeom() {
  }

  PlaneWriteDataGeom(const double x, const double y, const double z, const double nx, const double ny, const double nz) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
  }

  void flagCvs(Ugp * ugp) {

    // flag the cv's in Ugp intersecting the plane defined
    // by the point and normal...

    double tol = 1.0E-6;

    for (int ino = 0; ino<ugp->getNno(); ino++) {
      double delta = (ugp->x_no[ino][0]-x)*nx+(ugp->x_no[ino][1]-y)*ny+(ugp->x_no[ino][2]-z)*nz;
      if (delta>tol) {
        ugp->no_flag[ino] = 1;
      }
      else if (delta<-tol) {
        ugp->no_flag[ino] = -1;
      }
      else {
        ugp->no_flag[ino] = 0;
      }
    }

    for (int icv = 0; icv<ugp->getNcv(); icv++) {
      // use the flag to store the status of nodes so far...
      // flag =  0 : no status
      // flag = -1 : one or more -1 nodes found
      // flag = +1 : one or more +1 nodes found
      // flag = +2 : flag this cv
      int flag = 0;
      int foc_f = ugp->faocv_i[icv];
      int foc_l = ugp->faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++) {
        int ifa = ugp->faocv_v[foc];
        int nof_f = ugp->noofa_i[ifa];
        int nof_l = ugp->noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++) {
          int ino = ugp->noofa_v[nof];
          switch (ugp->no_flag[ino]) {
          case 0:
            flag = 2;
            break;
          case -1:
            switch (flag) {
            case 0:
              flag = -1;
              break;
            case 1:
              flag = 2;
              break;
            }
            break;
          case 1:
            switch (flag) {
            case 0:
              flag = 1;
              break;
            case -1:
              flag = 2;
              break;
            }
            break;
          }
          if (flag==2) break;
        }
        if (flag==2) break;
      }
      if (flag==2) ugp->cv_flag[icv] = 1;
    }

  }

};

class CylinderWriteDataGeom: public WriteDataGeom {

private:

  double x, y, z, nx, ny, nz, r;

public:

  virtual ~CylinderWriteDataGeom() {
  }

  CylinderWriteDataGeom(const double x, const double y, const double z, const double nx, const double ny,
      const double nz, const double r) {

    // point on cylinder axis...
    this->x = x;
    this->y = y;
    this->z = z;

    double mag = sqrt(nx*nx+ny*ny+nz*nz);

    // unit normal...
    this->nx = nx/mag;
    this->ny = ny/mag;
    this->nz = nz/mag;

    // radius...
    this->r = r;
  }

  void flagCvs(Ugp * ugp) {

    // flag the cv's in Ugp intersecting the cylinder defined
    // by the point, normal, radius...

    double tol = 1.0E-6;

    for (int ino = 0; ino<ugp->getNno(); ino++) {
      double dx = ugp->x_no[ino][0]-x;
      double dy = ugp->x_no[ino][1]-y;
      double dz = ugp->x_no[ino][2]-z;
      // so the distance along the normal is...
      double delta = dx*nx+dy*ny+dz*nz;
      // so the final delta radius is...
      delta = sqrt(fabs(dx*dx+dy*dy+dz*dz-delta*delta))-r;

      if (delta>tol) {
        ugp->no_flag[ino] = 1;
      }
      else if (delta<-tol) {
        ugp->no_flag[ino] = -1;
      }
      else {
        ugp->no_flag[ino] = 0;
      }
    }

    for (int icv = 0; icv<ugp->getNcv(); icv++) {
      // use the flag to store the status of nodes so far...
      // flag =  0 : no status
      // flag = -1 : one or more -1 nodes found
      // flag = +1 : one or more +1 nodes found
      // flag = +2 : flag this cv
      int flag = 0;
      int foc_f = ugp->faocv_i[icv];
      int foc_l = ugp->faocv_i[icv+1]-1;
      for (int foc = foc_f; foc<=foc_l; foc++) {
        int ifa = ugp->faocv_v[foc];
        int nof_f = ugp->noofa_i[ifa];
        int nof_l = ugp->noofa_i[ifa+1]-1;
        for (int nof = nof_f; nof<=nof_l; nof++) {
          int ino = ugp->noofa_v[nof];
          switch (ugp->no_flag[ino]) {
          case 0:
            flag = 2;
            break;
          case -1:
            switch (flag) {
            case 0:
              flag = -1;
              break;
            case 1:
              flag = 2;
              break;
            }
            break;
          case 1:
            switch (flag) {
            case 0:
              flag = 1;
              break;
            case -1:
              flag = 2;
              break;
            }
            break;
          }
          if (flag==2) break;
        }
        if (flag==2) break;
      }
      if (flag==2) ugp->cv_flag[icv] = 1;
    }

  }

};

class BoxWriteDataGeom: public WriteDataGeom {

private:

  double xyzMin[3], xyzMax[3];

public:

  virtual ~BoxWriteDataGeom() {
  }

  BoxWriteDataGeom(const double x0, const double x1, const double y0, const double y1, const double z0, const double z1) {

    this->xyzMin[0] = x0;
    this->xyzMin[1] = y0;
    this->xyzMin[2] = z0;

    this->xyzMax[0] = x1;
    this->xyzMax[1] = y1;
    this->xyzMax[2] = z1;

  }

  void flagCvs(Ugp * ugp) {

    for (int icv = 0; icv<ugp->getNcv(); icv++) {

      // assume we are going to dump this cv...
      ugp->cv_flag[icv] = 1;

      const int noc_f = ugp->noocv_i[icv];
      const int noc_l = ugp->noocv_i[icv+1]-1;

      for (int i = 0; i<3; ++i) {
        int flag = 0;
        for (int noc = noc_f; noc<=noc_l; ++noc) {
          const int ino = ugp->noocv_v[noc];
          if (ugp->x_no[ino][i]<xyzMin[i]) {
            // this node is on the left of the min...
            switch (flag) {
            case 0:
            case -1:
              flag = -1;
              break;
            case +1:
              flag = 0;
              break;
            }
          }
          else if (ugp->x_no[ino][i]>xyzMax[i]) {
            // this node is on the right of the max...
            switch (flag) {
            case 0:
            case +1:
              flag = +1;
              break;
            case -1:
              flag = 0;
              break;
            }
          }
          else {
            // this node is in the interval...
            flag = 0;
          }
          // if the flag is zero, break out...
          if (flag==0) break;
        }
        // a non-zero flag at this point means that
        // the nodes are entirely on one size of the box, so
        if (flag!=0) {
          ugp->cv_flag[icv] = 0;
          break;
        }
      }

    }

  }

};

// #############################################

enum WriteDataFormats {
  NO_FORMAT, TECPLOT_FORMAT
};

class WriteData {

private:

  int format, interval;
  string name;
  WriteDataGeom *geom;
  list<string> varList;

public:

  virtual ~WriteData() {
  }

  WriteData(Param *p) {

    format = NO_FORMAT;
    interval = -1;
    int name_set = 0;
    int vars_set = 0;
    geom = NULL;

    int ierr = 0;
    int i = 1;
    while (i<p->getSize()) {
      string token = p->getString(i++);
      if (token=="FORMAT") {
        token = p->getString(i++);
        if (token=="TECPLOT") {
          format = TECPLOT_FORMAT;
        }
        else {
          if (mpi_rank==0) cerr<<"Error: unsuported WRITE_DATA FORMAT: \""<<token<<"\""<<endl;
          ierr = 1;
        }
      }
      else if (token=="NAME") {
        name = p->getString(i++);
        name_set = 1;
      }
      else if (token=="INTERVAL") {
        interval = p->getInt(i++);
      }
      else if (token=="GEOM") {
        token = p->getString(i++);
        if (token=="PLANE") {
          double x = p->getDouble(i++);
          double y = p->getDouble(i++);
          double z = p->getDouble(i++);
          double nx = p->getDouble(i++);
          double ny = p->getDouble(i++);
          double nz = p->getDouble(i++);
          geom = new PlaneWriteDataGeom(x, y, z, nx, ny, nz);
        }
        else if (token=="CYLINDER") {
          double x = p->getDouble(i++);
          double y = p->getDouble(i++);
          double z = p->getDouble(i++);
          double nx = p->getDouble(i++);
          double ny = p->getDouble(i++);
          double nz = p->getDouble(i++);
          double r = p->getDouble(i++);
          geom = new CylinderWriteDataGeom(x, y, z, nx, ny, nz, r);
        }
        else if (token=="BOX") {
          double x0 = p->getDouble(i++);
          double x1 = p->getDouble(i++);
          double y0 = p->getDouble(i++);
          double y1 = p->getDouble(i++);
          double z0 = p->getDouble(i++);
          double z1 = p->getDouble(i++);
          geom = new BoxWriteDataGeom(x0, x1, y0, y1, z0, z1);
        }
        else if (token=="ALL") {
          geom = new AllWriteDataGeom();
        }
        else {
          if (mpi_rank==0) cerr<<"Error: unsuported WRITE_DATA GEOM: \""<<token<<"\""<<endl;
          ierr = 1;
        }
      }
      else if (token=="VARS") {
        // the var list is a variable length, so must be
        // listed at the end...
        while (i<p->getSize()) {
          varList.push_back(p->getString(i++));
          vars_set++;
        }
      }
      else {
        if (mpi_rank==0) cerr<<"Error: unrecognized WRITE_DATA token: "<<token<<endl;
        ierr = 1;
      }

    }

    if (format==NO_FORMAT) {
      if (mpi_rank==0) cout<<" > WRITE_DATA: missing FORMAT. using default: TECPLOT"<<endl;
      format = TECPLOT_FORMAT;
    }

    if (interval<=0) {
      if (mpi_rank==0) cout<<" > WRITE_DATA: missing INTERVAL. using default: 1"<<endl;
      interval = 1;
    }

    if (name_set==0) {
      if (mpi_rank==0) cerr<<"Error: WRITE_DATA missing valid NAME"<<endl;
      ierr = 1;
    }

    if (geom==NULL) {
      if (mpi_rank==0) cout<<" > WRITE_DATA: missing GEOM. using default: ALL"<<endl;
      geom = new AllWriteDataGeom();
    }

    if (vars_set==0) {
      if (mpi_rank==0) cerr<<"Error: WRITE_DATA missing valid VARS"<<endl;
      ierr = 1;
    }

    if (ierr!=0) throw(-1);

  }

  void write(const int step, UgpWithTools1 * ugp) {
    if (step%interval==0) {
      // set the data flags...
      ugp->clearDataFlags();
      for (list<string>::iterator var = varList.begin(); var!=varList.end(); var++) {
        if (var->compare("STATS")==0) ugp->setStatsFlags();
        else ugp->setDataFlag(var->c_str());
      }

      // if a geom is specified, then we flag a bunch of cvs and write...
      if (geom!=NULL) {
        for (int icv = 0; icv<ugp->getNcv(); icv++)
          ugp->cv_flag[icv] = 0;
        geom->flagCvs(ugp);

        char filename[32];
        sprintf(filename, "%s.%06d.plt", name.c_str(), step);
        ugp->writeFlaggedCvsTecplot(filename);
      }
    }
  }
};

// ========================================================
// Probe class gets instantiated for each PROBE parameter
// line...
// ========================================================

enum WriteProbeFormats {
  ASCII, BINARY
};

class Probe {

private:

  string name;
  int np;
  double (*xp)[3];
  int interval;
  list<string> varList;
  int first;

  //added to write probe files as PLOT3D format
  int ni, nj, nk;
  int format;
  int plot3d_format;
  int write_separate;
  //

  int * cvopp_i;
  int cvopp_s;
  int * cvopp_v;
  double * cvopp_wgt;

  int * noopp_i;
  int noopp_s;
  int * noopp_v;
  double * noopp_wgt;

  /*
   int no_data_rank;
   int no_data_size;
   int no_data_index[UGP_NNOC_MAX];
   double no_data_weight[UGP_NNOC_MAX];
   int cv_data_rank;
   int cv_data_size;
   int cv_data_index[1];
   double cv_data_weight[1];
   */

public:

  Probe(Param *p) {

    // initialize structs...
    cvopp_i = NULL;
    cvopp_s = 0;
    cvopp_v = NULL;
    cvopp_wgt = NULL;

    noopp_i = NULL;
    noopp_s = 0;
    noopp_v = NULL;
    noopp_wgt = NULL;

    first = 1;
    interval = -1;
    int name_set = 0;
    int geom_set = 0;
    int vars_set = 0;

    format = ASCII;
    plot3d_format = 0;
    write_separate = 0;

    // use this to store the np probe points, no matter what the 
    // probe geom...

    int ierr = 0;
    int i = 1;
    while (i<p->getSize()) {
      string token = p->getString(i++);
      if (token=="NAME") {
        assert(name_set==0);
        name = p->getString(i++);
        name_set = 1;
      }
      else if (token=="INTERVAL") {
        interval = p->getInt(i++);
        assert(interval>0);
      }
      else if (token=="BINARY") {
        format = BINARY;
      }
      else if (token=="WRITE_SEPARATE") {
        write_separate = 1;
      }

      else if (token=="GEOM") {
        token = p->getString(i++);
        if (token=="POINT") {
          np = 1;
          xp = new double[np][3];
          xp[0][0] = p->getDouble(i++);
          xp[0][1] = p->getDouble(i++);
          xp[0][2] = p->getDouble(i++);
          geom_set = 1;
        }
        else if (token=="LINE") {
          // LINE <x0> <y0> <z0> <x1> <y1> <z1> <n>
          double x0 = p->getDouble(i++);
          double y0 = p->getDouble(i++);
          double z0 = p->getDouble(i++);
          double x1 = p->getDouble(i++);
          double y1 = p->getDouble(i++);
          double z1 = p->getDouble(i++);
          np = p->getInt(i++);
          assert(np>=2);
          xp = new double[np][3];
          for (int ip = 0; ip<np; ++ip) {
            xp[ip][0] = x0+(double) ip/(double) (np-1)*(x1-x0);
            xp[ip][1] = y0+(double) ip/(double) (np-1)*(y1-y0);
            xp[ip][2] = z0+(double) ip/(double) (np-1)*(z1-z0);
          }
          geom_set = 1;
        }
        else if (token=="CIRCLE_XRN") {
          double x = p->getDouble(i++);
          double r = p->getDouble(i++);
          np = p->getInt(i++);
          xp = new double[np][3];
          for (int ip = 0; ip<np; ++ip) {
            // use point spacing staggered off the axes...
            double theta = (double) (2*ip+1)/(double) np*M_PI;
            xp[ip][0] = x;
            xp[ip][1] = r*cos(theta);
            xp[ip][2] = r*sin(theta);
          }
          geom_set = 1;
        }
        else if (token=="CIRCLE_ZRN") {
          double z = p->getDouble(i++);
          double r = p->getDouble(i++);
          np = p->getInt(i++);
          xp = new double[np][3];
          for (int ip = 0; ip<np; ++ip) {
            // use point spacing staggered off teh axes...
            double theta = (double) (2*ip+1)/(double) np*M_PI;
            xp[ip][0] = r*cos(theta);
            xp[ip][1] = r*sin(theta);
            xp[ip][2] = z;
          }
          geom_set = 1;
        }
        else if (token=="BOX") {
          // the BOX probe samples the data in a structured array
          // BOX <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <nx> <ny> <nz>
          double xmin = p->getDouble(i++);
          double xmax = p->getDouble(i++);
          double ymin = p->getDouble(i++);
          double ymax = p->getDouble(i++);
          double zmin = p->getDouble(i++);
          double zmax = p->getDouble(i++);
          assert(xmax>xmin);
          assert(ymax>ymin);
          assert(zmax>zmin);
          int nx = p->getInt(i++);
          int ny = p->getInt(i++);
          int nz = p->getInt(i++);
          assert(nx>=1);
          assert(ny>=1);
          assert(nz>=1);
          np = nx*ny*nz;
          xp = new double[np][3];
          int ip = 0;
          for (int k = 0; k<nz; ++k) {
            double this_z;
            if (nz==1) this_z = 0.5*(zmin+zmax);
            else this_z = zmin+(double) k/(double) (nz-1)*(zmax-zmin);
            for (int j = 0; j<ny; ++j) {
              double this_y;
              if (ny==1) this_y = 0.5*(ymin+ymax);
              else this_y = ymin+(double) j/(double) (ny-1)*(ymax-ymin);
              for (int i = 0; i<nx; ++i) {
                double this_x;
                if (nx==1) this_x = 0.5*(xmin+xmax);
                else this_x = xmin+(double) i/(double) (nx-1)*(xmax-xmin);
                // set the point coords...
                xp[ip][0] = this_x;
                xp[ip][1] = this_y;
                xp[ip][2] = this_z;
                ++ip;
              }
            }
          }
          assert(ip==np);
          geom_set = 1;
        }
        else if (token=="FILE") {
          // an ASCII input file with plot3d format
          // example :
          // PROBE NAME=probes/p1 INTERVAL=20 WRITE_SEPARATE BINARY GEOM=FILE mesh.xyz VARS=RHO P T U
          // writes separate binary files in plot3d format at every 20 timestep
          plot3d_format = 1;

          FILE *fpout;

          string filename = p->getString(i++);

          if ((fpout = fopen(filename.c_str(), "r"))==NULL) {
            cerr<<"Error: cannot open probe file "<<filename<<". Check if directory exists."<<endl;
            ierr = 1;
          }

          fscanf(fpout, " %9i %9i %9i\n", &ni, &nj, &nk);

          np = ni*nj*nk;

          assert(np>=2);
          xp = new double[np][3];

          int ip = 0;
          for (int k = 0; k<nk; k++)
            for (int j = 0; j<nj; j++) {
              for (int i = 0; i<ni; i++) {
                fscanf(fpout, "%lf", &xp[ip][0]);
                ++ip;
              }
            }

          ip = 0;

          for (int k = 0; k<nk; k++)
            for (int j = 0; j<nj; j++) {
              for (int i = 0; i<ni; i++) {
                fscanf(fpout, "%lf", &xp[ip][1]);
                ++ip;
              }
            }

          ip = 0;

          for (int k = 0; k<nk; k++)
            for (int j = 0; j<nj; j++) {
              for (int i = 0; i<ni; i++) {
                fscanf(fpout, "%lf", &xp[ip][2]);
                ++ip;
              }
            }

          fclose(fpout);
          geom_set = 1;
        }
        else {
          if (mpi_rank==0) {
            cerr<<"Error: unrecognized PROBE GEOM: "<<token<<endl;
            cerr<<"valid GEOMS: POINT, CIRCLE_ZRN, CIRCLE_XRN"<<endl;
          }
          ierr = 1;
        }
      }
      else if (token=="VARS") {
        // the var list is a variable length, so must be
        // listed at the end...
        assert(name_set==1);
        assert(interval>0);
        assert(geom_set==1);
        while (i<p->getSize()) {
          varList.push_back(p->getString(i++));
          vars_set++;
        }
      }
      else {
        if (mpi_rank==0) cerr<<"Error: unrecognized PROBE token: "<<token<<endl;
        ierr = 1;
      }
    }

    if (name_set==0) {
      if (mpi_rank==0) cerr<<"Error: PROBE missing valid NAME"<<endl;
      ierr = 1;
    }

    if (interval<=0) {
      if (mpi_rank==0) cerr<<"Error: PROBE missing valid INTERVAL"<<endl;
      ierr = 1;
    }

    if (geom_set==0) {
      if (mpi_rank==0) cerr<<"Error: PROBE missing valid GEOM"<<endl;
      ierr = 1;
    }

    if (vars_set==0) {
      if (mpi_rank==0) cerr<<"Error: PROBE missing valid VARS"<<endl;
      ierr = 1;
    }

    // we should have some points...
    assert(np>0);

    // rank 0 should try to write the probe description...
    if (mpi_rank==0) {
      cout<<" > adding PROBE: "<<name<<", number of probe points: "<<np<<endl;
      // also write the description of this probe to a file...
      FILE * fp;
      char filename[64];
      sprintf(filename, "%s.README", name.c_str());
      if ((fp = fopen(filename, "a"))==NULL) {
        cerr<<"Error: cannot open probe file "<<filename<<". Check if directory exists."<<endl;
        ierr = 1;
      }
      else {
        time_t rawtime;
        time(&rawtime);
        fprintf(fp, "Local time: %s", ctime(&rawtime)); // adds its own \n
        fprintf(fp, "Job size: %d\n", mpi_size);
        fprintf(fp, "Parameter line: ");
        p->writeToFile(fp); // also adds its own parameter line
        fclose(fp);
      }
    }
    MPI_Bcast(&ierr, 1, MPI_INT, 0, mpi_comm);

    if (ierr!=0) throw(-1);

  }

public:

  string getName() const {
    return (name);
  }

public:

  void dump(const int step, const double time, UgpWithTools0 * ugp) {

    if (step%interval==0) {

      // for buffering the data...

      double * my_data = new double[np];
      double * data = NULL;
      if (mpi_rank==0) data = new double[np];

      // ==========================================
      // setup everything on the first time...
      // ==========================================

      if (first==1) {

        // reset first...
        first = 0;

        // on the very first time, write or append the probe_name.X,Y,Z files...
        if (mpi_rank==0) {

          FILE * fp;
          char filename[64];

          sprintf(filename, "%s.README", name.c_str());
          if ((fp = fopen(filename, "a"))==NULL) {
            cerr<<"Error: cannot open file "<<filename<<endl;
            throw(-1);
          }
          fprintf(fp, "Initial step: %d\n", step);
          fprintf(fp, "Initial time: %18.16e\n\n", time); // leave an extra space..
          fclose(fp);

          // X...
          sprintf(filename, "%s.X", name.c_str());
          if ((fp = fopen(filename, "a"))==NULL) {
            cerr<<"Error: cannot open file "<<filename<<endl;
            throw(-1);
          }
          fprintf(fp, "%d %18.16e %d", step, time, np);
          for (int ip = 0; ip<np; ++ip)
            fprintf(fp, " %18.16e", xp[ip][0]);
          fprintf(fp, "\n");
          fclose(fp);

          // Y...
          sprintf(filename, "%s.Y", name.c_str());
          if ((fp = fopen(filename, "a"))==NULL) {
            cerr<<"Error: cannot open file "<<filename<<endl;
            throw(-1);
          }
          fprintf(fp, "%d %18.16e %d", step, time, np);
          for (int ip = 0; ip<np; ++ip)
            fprintf(fp, " %18.16e", xp[ip][1]);
          fprintf(fp, "\n");
          fclose(fp);

          // Z...
          sprintf(filename, "%s.Z", name.c_str());
          if ((fp = fopen(filename, "a"))==NULL) {
            cerr<<"Error: cannot open file "<<filename<<endl;
            throw(-1);
          }
          fprintf(fp, "%d %18.16e %d", step, time, np);
          for (int ip = 0; ip<np; ++ip)
            fprintf(fp, " %18.16e", xp[ip][2]);
          fprintf(fp, "\n");
          fclose(fp);

        }

        // calc the cv offset...
        int my_cv_offset;
        MPI_Scan(&(ugp->ncv), &my_cv_offset, 1, MPI_INT, MPI_SUM, mpi_comm);
        my_cv_offset -= (ugp->ncv);

        // locate the points...
        int nList;
        int cvList[ADT_LIST_MAX];
        int * my_icv_flag = new int[np];
        for (int ip = 0; ip<np; ++ip) {
          my_icv_flag[ip] = -1;
          ugp->cvAdt->buildListForPoint(nList, cvList, xp[ip]);
          for (int i = 0; i<nList; i++) {
            const int icv = cvList[i];
            int foc, ifa;
            if (ugp->pointIsInsideCv(foc, ifa, xp[ip], icv)) {
              my_icv_flag[ip] = icv+my_cv_offset; // use a global cv index
              break;
            }
          }
        }

        // at this point, my_icv_flag should contain the global cv index of
        // the containing cell, or -1...

        int * icv_flag = new int[np];
        MPI_Allreduce(my_icv_flag, icv_flag, np, MPI_INT, MPI_MAX, mpi_comm);

        // any -1's in icv_flag indicate probe points that could not be
        // found...
        int ierr = 0;
        for (int ip = 0; ip<np; ++ip) {
          if (icv_flag[ip]==-1) {
            if (mpi_rank==0) cerr<<"Error: PROBE "<<name<<": following point could not be located: "<<xp[ip][0]<<" "
                <<xp[ip][1]<<" "<<xp[ip][2]<<endl;
            ierr = 1;
          }
        }

        // determine which interpolation data structures we need...

        int noData = 0;
        int cvData = 0;
        for (list<string>::iterator var = varList.begin(); var!=varList.end(); var++) {
          // DoubleScalars first...
          DoubleScalar * scalarData = ugp->getScalarData(var->c_str());
          if (scalarData!=NULL) {
            if (scalarData->getDatatype()==NO_DATA) noData = 1;
            else if (scalarData->getDatatype()==CV_DATA) cvData = 1;
            else {
              if (mpi_rank==0) cerr<<"Error: VAR: "<<*var<<" is not NO_DATA or CV_DATA."<<endl;
              ierr = 1;
            }
          }
          else {
            DoubleVector * vectorData = ugp->getVectorData(var->c_str());
            if (vectorData!=NULL) {
              if (vectorData->getDatatype()==NO_DATA) noData = 1;
              else if (vectorData->getDatatype()==CV_DATA) cvData = 1;
              else {
                if (mpi_rank==0) cerr<<"Error: VAR: "<<*var<<" is not NO_DATA or CV_DATA."<<endl;
                ierr = 1;
              }
            }
            else {
              if (mpi_rank==0) cerr<<"Error: VAR: "<<*var<<" cannot be found in registered data."<<endl;
              ierr = 1;
            }
          }
        }

        if (ierr!=0) throw(-1);

        // ===============================================
        // for cv's just report centroid data...
        // ===============================================

        if (cvData==1) {

          assert(cvopp_i==NULL);
          cvopp_i = new int[np+1]; // cv-of-probe-point...

          for (int iter = 0; iter<2; ++iter) {
            cvopp_i[0] = 0;
            for (int ip = 0; ip<np; ++ip) {
              int cop = cvopp_i[ip];
              if (icv_flag[ip]==my_icv_flag[ip]) {
                // we own this point...
                if (iter==0) {
                  cop++;
                }
                else {
                  cvopp_v[cop] = icv_flag[ip]-my_cv_offset;
                  cvopp_wgt[cop] = 1.0;
                  cop++;
                }
              }
              cvopp_i[ip+1] = cop;
            }
            if (iter==0) {
              cvopp_s = cvopp_i[np];
              cvopp_v = new int[cvopp_s];
              cvopp_wgt = new double[cvopp_s];
            }
          }

          // loop on coordinate index...

          double my_dx_max[3] = { 0.0, 0.0, 0.0 };

          FOR_I3 {

            for (int ip = 0; ip<np; ip++) {
              my_data[ip] = 0.0;
              const int cop_f = cvopp_i[ip];
              const int cop_l = cvopp_i[ip+1]-1;
              for (int cop = cop_f; cop<=cop_l; ++cop) {
                const int icv = cvopp_v[cop];
                assert((icv>=0)&&(icv<ugp->ncv));
                // use the simple average of nodes to report the "approx"
                // cv center - this should be ok for checking purposes, and
                // will certainly identify when overlapped cells are present...
                const int noc_f = ugp->noocv_i[icv];
                const int noc_l = ugp->noocv_i[icv+1]-1;
                const double this_wgt = cvopp_wgt[cop]/(double) (noc_l-noc_f+1);
                for (int noc = noc_f; noc<=noc_l; ++noc) {
                  const int ino = ugp->noocv_v[noc];
                  my_data[ip] += this_wgt*ugp->x_no[ino][i];
                }
                my_dx_max[i] = max(fabs(my_data[ip]-xp[ip][i]), my_dx_max[i]);
              }
            }

            MPI_Reduce(my_data, data, np, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

            if (mpi_rank==0) {
              // filename depends on "i"...
              char filename[64];
              switch (i) {
              case 0:
                sprintf(filename, "%s.X_CV", name.c_str());
                break;
              case 1:
                sprintf(filename, "%s.Y_CV", name.c_str());
                break;
              case 2:
                sprintf(filename, "%s.Z_CV", name.c_str());
                break;
              }
              // write ...
              FILE * fp;
              if ((fp = fopen(filename, "a"))==NULL) {
                cerr<<"Error: cannot open file "<<filename<<endl;
                throw(-1);
              }
              fprintf(fp, "%d %18.16e %d", step, time, np);
              for (int ip = 0; ip<np; ++ip)
                fprintf(fp, " %18.16e", data[ip]);
              fprintf(fp, "\n");
              fclose(fp);
            }

          }

          double dx_max[3];
          MPI_Reduce(my_dx_max, dx_max, 3, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
          if (mpi_rank==0) cout<<" > max delta in probe points relative to cv centers: "<<dx_max[0]<<" "<<dx_max[1]
              <<" "<<dx_max[2]<<endl;

        }

        // ===============================================
        // for node-based data, just report nearest-node data
        // ===============================================

        if (noData==1) {

          double my_d2_max = 0.0;

          assert(noopp_i==NULL);
          noopp_i = new int[np+1]; // node-of-probe-point...

          for (int iter = 0; iter<2; ++iter) {
            noopp_i[0] = 0;
            for (int ip = 0; ip<np; ++ip) {
              int nop = noopp_i[ip];
              if (icv_flag[ip]==my_icv_flag[ip]) {
                // we own this point...
                if (iter==0) {
                  // for now choose one point, but in the future, use linear interpolation
                  nop++;
                }
                else {
                  int icv = icv_flag[ip]-my_cv_offset;
                  int noc_f = ugp->noocv_i[icv];
                  int noc_l = ugp->noocv_i[icv+1]-1;
                  int ino_min = -1;
                  double d2_min;
                  for (int noc = noc_f; noc<=noc_l; ++noc) {
                    const int ino = ugp->noocv_v[noc];
                    double d2 = 0.0;
                    FOR_I3 {
                      double dx = ugp->x_no[ino][i]-xp[ip][i];
                      d2 += dx*dx;
                    }
                    if ((ino_min==-1)||(d2<d2_min)) {
                      ino_min = ino;
                      d2_min = d2;
                    }
                  }
                  assert(ino_min!=-1);
                  my_d2_max = max(my_d2_max, d2_min);
                  noopp_v[nop] = ino_min;
                  noopp_wgt[nop] = 1.0;
                  nop++;
                }
              }
              noopp_i[ip+1] = nop;
            }
            if (iter==0) {
              noopp_s = noopp_i[np];
              noopp_v = new int[noopp_s];
              noopp_wgt = new double[noopp_s];
            }
          }

          double d2_max;
          MPI_Reduce(&my_d2_max, &d2_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
          if (mpi_rank==0) cout<<" > max delta in probe points relative to nodes: "<<sqrt(d2_max)<<endl;

          // loop on coordinate index...

          double my_dx_max[3] = { 0.0, 0.0, 0.0 };

          FOR_I3 {

            for (int ip = 0; ip<np; ip++) {
              my_data[ip] = 0.0;
              const int nop_f = noopp_i[ip];
              const int nop_l = noopp_i[ip+1]-1;
              for (int nop = nop_f; nop<=nop_l; ++nop) {
                const int ino = noopp_v[nop];
                assert((ino>=0)&&(ino<ugp->nno));
                const double this_wgt = noopp_wgt[nop];
                my_data[ip] += this_wgt*ugp->x_no[ino][i];
                my_dx_max[i] = max(fabs(my_data[ip]-xp[ip][i]), my_dx_max[i]);
              }
            }

            MPI_Reduce(my_data, data, np, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

            if (mpi_rank==0) {
              // filename depends on "i"...
              char filename[64];
              switch (i) {
              case 0:
                sprintf(filename, "%s.X_NO", name.c_str());
                break;
              case 1:
                sprintf(filename, "%s.Y_NO", name.c_str());
                break;
              case 2:
                sprintf(filename, "%s.Z_NO", name.c_str());
                break;
              }
              // write ...
              FILE * fp;
              if ((fp = fopen(filename, "a"))==NULL) {
                cerr<<"Error: cannot open file "<<filename<<endl;
                throw(-1);
              }
              fprintf(fp, "%d %18.16e %d", step, time, np);
              for (int ip = 0; ip<np; ++ip)
                fprintf(fp, " %18.16e", data[ip]);
              fprintf(fp, "\n");
              fclose(fp);
            }

          }

          double dx_max[3];
          MPI_Reduce(my_dx_max, dx_max, 3, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
          if (mpi_rank==0) cout<<" > max delta in probe points relative to node interpolation: "<<dx_max[0]<<" "
              <<dx_max[1]<<" "<<dx_max[2]<<endl;

        }

        // cleanup...
        delete[] my_icv_flag;
        delete[] icv_flag;
        delete[] xp;

      }

      // ==============================================
      // now dump the data...
      // ==============================================

      for (list<string>::iterator var = varList.begin(); var!=varList.end(); var++) {
        // DoubleScalars first...
        DoubleScalar * scalarData = ugp->getScalarData(var->c_str());
        if (scalarData!=NULL) {
          if (scalarData->getDatatype()==NO_DATA) {
            for (int ip = 0; ip<np; ip++) {
              my_data[ip] = 0.0;
              const int nop_f = noopp_i[ip];
              const int nop_l = noopp_i[ip+1]-1;
              for (int nop = nop_f; nop<=nop_l; ++nop) {
                const int ino = noopp_v[nop];
                my_data[ip] += noopp_wgt[nop]*(*(scalarData->ptr))[ino];
              }
            }
          }
          else if (scalarData->getDatatype()==CV_DATA) {
            for (int ip = 0; ip<np; ip++) {
              my_data[ip] = 0.0;
              const int cop_f = cvopp_i[ip];
              const int cop_l = cvopp_i[ip+1]-1;
              for (int cop = cop_f; cop<=cop_l; ++cop) {
                const int icv = cvopp_v[cop];
                my_data[ip] += cvopp_wgt[cop]*(*(scalarData->ptr))[icv];
              }
            }
          }

          // my_data has been populated above...
          MPI_Reduce(my_data, data, np, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

          // only rank 0 writes...
          if (mpi_rank==0) {
            char filename[64];

            if (write_separate==1) {
              sprintf(filename, "%s.%s.%d", name.c_str(), var->c_str(), step);
            }
            else {
              sprintf(filename, "%s.%s", name.c_str(), var->c_str());
            }

            FILE * fp;
            if ((fp = fopen(filename, "a"))==NULL) {
              cerr<<"Error: cannot open file "<<filename<<endl;
              throw(-1);
            }
            if (format==ASCII) {

              if (plot3d_format==1) {
                fprintf(fp, "%d %d %d 1", ni, nj, nk);
                fprintf(fp, "\n");
              }

              // Yaser! what is up here?

              else if (write_separate==0) fprintf(fp, "%d %18.16e %d", step, time, np);

              for (int ip = 0; ip<np; ++ip)
                fprintf(fp, " %18.16e", data[ip]);
              fprintf(fp, "\n");

            }
            else {

              // binary output?...

              if (plot3d_format==1) {
                fwrite(&ni, sizeof(int), 1, fp);
                fwrite(&nj, sizeof(int), 1, fp);
                fwrite(&nk, sizeof(int), 1, fp);
                int dummy = 1;
                fwrite(&dummy, sizeof(int), 1, fp);
              }

              else if (write_separate==0) {
                fwrite(&step, sizeof(int), 1, fp);
                fwrite(&time, sizeof(double), 1, fp);
                fwrite(&np, sizeof(int), 1, fp);
              }
              fwrite(&data[0], sizeof(double), np, fp);
            }

            fclose(fp);
          }

        }
        else {

          DoubleVector * vectorData = ugp->getVectorData(var->c_str());
          if (vectorData!=NULL) {
            // loop on conponents...
            FOR_I3 {

              // pack local data buffer...
              if (vectorData->getDatatype()==NO_DATA) {
                for (int ip = 0; ip<np; ip++) {
                  my_data[ip] = 0.0;
                  const int nop_f = noopp_i[ip];
                  const int nop_l = noopp_i[ip+1]-1;
                  for (int nop = nop_f; nop<=nop_l; ++nop) {
                    const int ino = noopp_v[nop];
                    my_data[ip] += noopp_wgt[nop]*(*(vectorData->ptr))[ino][i];
                  }
                }
              }
              else if (vectorData->getDatatype()==CV_DATA) {
                for (int ip = 0; ip<np; ip++) {
                  my_data[ip] = 0.0;
                  const int cop_f = cvopp_i[ip];
                  const int cop_l = cvopp_i[ip+1]-1;
                  for (int cop = cop_f; cop<=cop_l; ++cop) {
                    const int icv = cvopp_v[cop];
                    my_data[ip] += cvopp_wgt[cop]*(*(vectorData->ptr))[icv][i];
                  }
                }
              }

              // my_data has been populated above...
              MPI_Reduce(my_data, data, np, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

              // only rank 0 writes...
              if (mpi_rank==0) {
                char filename[64];

                if (write_separate==1) {
                  sprintf(filename, "%s.%s.%d", name.c_str(), var->c_str(), step);
                }
                else {

                  switch (i) {
                  case 0:
                    sprintf(filename, "%s.%s-X", name.c_str(), var->c_str());
                    break;
                  case 1:
                    sprintf(filename, "%s.%s-Y", name.c_str(), var->c_str());
                    break;
                  case 2:
                    sprintf(filename, "%s.%s-Z", name.c_str(), var->c_str());
                    break;
                  }
                }
                FILE * fp;
                if ((fp = fopen(filename, "a"))==NULL) {
                  cerr<<"Error: cannot open file "<<filename<<endl;
                  throw(-1);
                }

                if (format==ASCII) {
                  if ((plot3d_format==1)&&(i==0)) {
                    fprintf(fp, "%d %d %d 3", ni, nj, nk);
                    fprintf(fp, "\n");
                  }

                  else if (write_separate==0) fprintf(fp, "%d %18.16e %d", step, time, np);

                  for (int ip = 0; ip<np; ++ip)
                    fprintf(fp, " %18.16e", data[ip]);
                  fprintf(fp, "\n");

                }
                else {
                  if ((plot3d_format==1)&&(i==0)) {
                    fwrite(&ni, sizeof(int), 1, fp);
                    fwrite(&nj, sizeof(int), 1, fp);
                    fwrite(&nk, sizeof(int), 1, fp);
                    int dummy = 3;
                    fwrite(&dummy, sizeof(int), 1, fp);
                  }

                  else if (write_separate==0) {
                    fwrite(&step, sizeof(int), 1, fp);
                    fwrite(&time, sizeof(double), 1, fp);
                    fwrite(&np, sizeof(int), 1, fp);
                  }
                  fwrite(&data[0], sizeof(double), np, fp);
                }

                fclose(fp);
              }

            }

          }

        }

      }

      delete[] my_data;
      if (mpi_rank==0) delete[] data;

    }

  }

#ifdef LATERLATER

private:

  void buildNoDataWeights(UgpWithTools0 * ugp)
  {

    int my_no_data_rank = mpi_size;
    if (ugp->setNodeInterpWeightsForPoint(no_data_size, no_data_index, no_data_weight, xyz))
    my_no_data_rank = mpi_rank;

    // min rank process that found it wins...
    MPI_Allreduce(&my_no_data_rank, &no_data_rank, 1, MPI_INT, MPI_MIN, mpi_comm);

    if (no_data_rank == mpi_size)
    {
      if (mpi_rank == 0)
      cerr << "Error: Probe " << name << ": could not locate point in grid: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
      throw(-1);
    }

    // do a check if it is us...
    if (no_data_rank == mpi_rank)
    {
      double x_test[3] =
      { 0.0, 0.0, 0.0};
      for (int i = 0; i < no_data_size; i++)
      for (int j = 0; j < 3; j++)
      x_test[j] += no_data_weight[i] * ugp->x_no[no_data_index[i]][j];
      assert(fabs(x_test[0] - xyz[0]) < 1.0E-8);
      assert(fabs(x_test[1] - xyz[1]) < 1.0E-8);
      assert(fabs(x_test[2] - xyz[2]) < 1.0E-8);
    }

  }

  void buildCvDataWeights(UgpWithTools0 * ugp)
  {

    int my_cv_data_rank = mpi_size;
    if (ugp->setCvInterpWeightsForPoint(cv_data_size, cv_data_index, cv_data_weight, xyz))
    my_cv_data_rank = mpi_rank;

    // min rank process that found it wins...
    MPI_Allreduce(&my_cv_data_rank, &cv_data_rank, 1, MPI_INT, MPI_MIN, mpi_comm);

    if (cv_data_rank == mpi_size)
    {
      if (mpi_rank == 0)
      cerr << "Error: Probe " << name << ": could not locate point in grid: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
      throw(-1);
    }

    if (cv_data_rank == mpi_rank)
    {
      assert( cv_data_size == 1 );
      const int icv = cv_data_index[0];
      double dx[3] =
      { 0.0,0.0,0.0};
      int noc_f = ugp->noocv_i[icv];
      int noc_l = ugp->noocv_i[icv+1]-1;
      for (int noc = noc_f; noc <= noc_l; ++noc)
      {
        const int ino = ugp->noocv_v[noc];
        for (int i = 0; i < 3; ++i)
        dx[i] += ugp->x_no[ino][i];
      }
      double d2 = 0.0;
      for (int i = 0; i < 3; ++i)
      {
        dx[i] = dx[i]/(double)(noc_l-noc_f+1) - xyz[i];
        d2 += dx[i]*dx[i];
      }
      cout << " > Warning: CV probes report closest CV values. Approx cv-probe distance: " << sqrt(d2) << endl;
    }

  }

  /*

   int interpRegisteredData(double * data,const string& name) {

   DoubleScalar * scalarData = getScalarData(name);
   if (scalarData != NULL) {

   // this is a scalar, check if it is node-based...
   switch (scalarData->datatype) {
   case NO_DATA:

   if (no_data_rank == -2)
   buildNoDataWeights();

   assert( no_data_rank >= 0 );

   if (no_data_rank == mpi_rank) {
   data[0] = 0.0;
   for (int i = 0; i < no_data_size; i++) {
   data[0] += no_data_weight[i]*(*scalarData->ptr)[no_data_index[i]];
   }
   }


   HEREHERE


   // returns the dimension of the interpolated data

   */

#endif

};

// ############################################################

class UgpWithTools: public UgpWithTools1 {

protected:

private:

  // for post-processing output...
  list<WriteData> writeDataList;

  // for probes...
  list<Probe> probeList;

public:

  UgpWithTools() {

    if (mpi_rank==0) cout<<"UgpWithTools()"<<endl;

  }

  void writeData(const int step) {
    for (list<WriteData>::iterator wd = writeDataList.begin(); wd!=writeDataList.end(); wd++)
      wd->write(step, this);

  }

  void initWriteData(ParamMap * params);

protected:

  void init() {

    if (mpi_rank==0) cout<<"UgpWithTools::init()"<<endl;

    UgpWithTools1::init();

  }

  /**

   Call initWriteData() in a class that inherits UgpWithTools to make
   WRITE_DATA capabilities available. For example:

   \verbatim
   vector<Param> *paramVec;
   if (getParam(paramVec,"WRITE_DATA"))
   initWriteData(paramVec);
   \endverbatim

   Where the paramVec contains the vector of WRITE_DATA params with
   general format:

   \verbatim
   WRITE_DATA FORMAT <format> NAME <name> INTERVAL <int> GEOM <geom> VARS <vars>
   \endverbatim

   A specific example is:

   \verbatim
   WRITE_DATA FORMAT TECPLOT NAME all INTERVAL 500 GEOM ALL VARS U P PHI1 PHI2 PHI3
   \endverbatim

   */

  void initProbes(ParamMap * params);

  void dumpProbes(const int step, const double time) {
    for (list<Probe>::iterator pr = probeList.begin(); pr!=probeList.end(); pr++)
      pr->dump(step, time, this);
  }

};

#endif

